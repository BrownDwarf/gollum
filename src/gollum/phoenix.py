r"""
PHOENIX Spectrum
-----------------

A container for a single Phoenix grid-point spectrum of wavelength and flux :math:`F(\lambda)`.

PHOENIXSpectrum
###############
"""

import os
import numpy as np

from copy import deepcopy
from itertools import product
from warnings import filterwarnings
from logging import getLogger
from tqdm import tqdm
from urllib.error import URLError
from gollum.precomputed_spectrum import PrecomputedSpectrum
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.io import fits
from astropy import units as u
from specutils import SpectrumCollection
from specutils.spectra.spectrum1d import Spectrum1D
from bokeh.io import show, output_notebook
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import Slider, Range1d
from bokeh.layouts import layout, Spacer
from bokeh.models.widgets import Button, Div

log = getLogger(__name__)

#  See Issue: https://github.com/astropy/specutils/issues/779
filterwarnings("ignore", category=AstropyDeprecationWarning)
# See Issue: https://github.com/astropy/specutils/issues/800
filterwarnings("ignore", category=RuntimeWarning)


class PHOENIXSpectrum(PrecomputedSpectrum):
    """
    A container for PHOENIX spectra

    Parameters
    ----------
    teff : int
        The teff label of the PHOENIX model to read in.  Must be on the PHOENIX grid.
    logg : float
        The logg label of the PHOENIX model to read in.  Must be on the PHOENIX grid.
    metallicity : float
        The metallicity label of the PHOENIX model to read in. Must be on the PHOENIX grid.
    path : str
        The path to your locally downloaded PHOENIX grid library. Default: "~/libraries/raw/PHOENIX/"
    download : bool
        [Experimental] Set to True if you want to download the spectra from the internet; requires an internet connection.
    wl_lo : float
        The shortest wavelength of the models to keep (\u212B)
    wl_hi : float
        The longest wavelength of the models to keep (\u212B)
    """

    def __init__(
        self,
        *args,
        teff=None,
        logg=None,
        metallicity=0.0,  # solar by default
        path="~/libraries/raw/PHOENIX/",
        download=False,
        wl_lo=8038,
        wl_hi=12849,
        **kwargs,
    ):

        if teff and logg:
            if not download:
                base_path = os.path.expanduser(path)
                assert os.path.exists(base_path), "Given path does not exist."

                wl_file = f"{base_path}/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
                assert os.path.exists(wl_file), f"PHOENIX models must be in {base_path}"
            else:
                site = "ftp://phoenix.astro.physik.uni-goettingen.de/v2.0/HiResFITS/"
                log.info("[WIP]Downloading PHOENIX models from the internet...")
                log.info(f"We are using this FTP site: {site}")
                wl_file = f"{site}WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
                base_path = f"{site}PHOENIX-ACES-AGSS-COND-2011/"

            Z_string = f"{metallicity:+0.1f}" if metallicity else "-0.0"

            fn = f"{base_path}/Z{Z_string}/lte{teff:05d}-{logg:0.2f}{Z_string}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
            if not os.path.exists(fn) and not download:
                # Have to add website requests, download check is a temporary fix
                raise FileExistsError(
                    "No PHOENIX Spectrum file exists for given parameters."
                )

            wl_orig = fits.open(wl_file)[0].data.astype(np.float64)
            mask = (wl_orig > wl_lo) & (wl_orig < wl_hi)
            wl_out = wl_orig[mask]

            flux_orig = fits.open(fn)[0].data.astype(np.float64)
            flux_native = flux_orig[mask]
            native_flux_unit = u.erg / u.s / u.cm ** 2 / u.cm

            super().__init__(
                spectral_axis=wl_out * u.AA,
                flux=flux_native * native_flux_unit,
                meta={
                    "teff": teff,
                    "logg": logg,
                    "metallicity": metallicity,
                    "native_flux_unit": native_flux_unit,
                },
                **kwargs,
            )

        else:
            super().__init__(*args, **kwargs)

    teff = property(lambda self: self.meta.get("teff"))
    logg = property(lambda self: self.meta.get("logg"))
    metallicity = property(lambda self: self.meta.get("metallicity"))


class PHOENIXGrid(SpectrumCollection):
    """
    A container for a grid of PHOENIX precomputed synthetic spectra of stars.

    Parameters
    ----------
    teff_range : tuple
        The teff limits of the grid model to read in.
    logg_range : tuple
        The logg limits of the grid model to read in.
    metallicity_range : tuple
        The metallicity limits of the grid model to read in.
    path : str
        The path to your locally downloaded PHOENIX grid library. Default: "~/libraries/raw/PHOENIX/"
    wl_lo : float
        The shortest wavelength of the models to keep (\u212B)
    wl_hi : float
        The longest wavelength of the models to keep (\u212B)
    """

    def __init__(
        self,
        teff_range=None,
        logg_range=None,
        metallicity_range=None,
        path="~/libraries/raw/PHOENIX/",
        wl_lo=8038,
        wl_hi=12849,
        **kwargs,
    ):

        if set(("flux", "spectral_axis", "meta")).issubset(kwargs):
            # Trigger a passthrough
            super().__init__(**kwargs)
        else:

            teff_points = np.hstack(
                (np.arange(2300, 7000, 100), np.arange(7000, 12001, 200))
            )
            # Todo: some T_eff ranges go to log(g) = 0.0, consider adding these
            logg_points = np.arange(2.0, 6.01, 0.5)

            metallicity_points = np.array([-4, -3, -2, -1.5, -1, -0.5, 0, 0.5, 1])

            if teff_range:
                subset = (teff_points >= teff_range[0]) & (teff_points <= teff_range[1])
                teff_points = teff_points[subset]

            if logg_range:
                subset = (logg_points >= logg_range[0]) & (logg_points <= logg_range[1])
                logg_points = logg_points[subset]

            if metallicity_range:
                subset = (metallicity_points >= metallicity_range[0]) & (
                    metallicity_points <= metallicity_range[1]
                )
                metallicity_points = metallicity_points[subset]

            wavelengths, fluxes, grid_points, missing = [], [], [], 0
            pbar = tqdm(
                product(teff_points, logg_points, metallicity_points),
                total=len(teff_points) * len(logg_points) * len(metallicity_points),
            )

            for teff, logg, Z in pbar:
                pbar.desc = f"Processing Teff={teff}K|log(g)={logg:0.2f}|Z={Z:+0.1f}"
                try:
                    spec = PHOENIXSpectrum(
                        teff=teff,
                        logg=logg,
                        metallicity=Z,
                        path=path,
                        wl_lo=wl_lo,
                        wl_hi=wl_hi,
                    )
                except (FileExistsError, URLError):
                    log.info(f"No file for Teff={teff}K|log(g)={logg:0.2f}|Z={Z:+0.1f}")
                    missing += 1
                wavelengths.append(spec.wavelength)
                fluxes.append(spec.flux)
                grid_points.append((teff, logg, Z))
            assert grid_points != [], "Empty grid; parameter limits out of range"
            print(
                f"{missing} files not found; grid may not cover given parameter ranges fully"
            ) if missing else None

            super().__init__(
                flux=np.array(fluxes) * fluxes[0].unit,
                spectral_axis=np.array(wavelengths) * wavelengths[0].unit,
                meta={
                    "teff_points": teff_points,
                    "logg_points": logg_points,
                    "metallicity_points": metallicity_points,
                    "grid_labels": ("T_eff", "log(g)", "Z"),
                    "n_spectra": len(grid_points),
                    "grid_points": grid_points,
                    "lookup_dict": {value: i for i, value in enumerate(grid_points)},
                },
            )

    def __getitem__(self, key):
        flux = self.flux[key]
        if flux.ndim != 1:
            raise ValueError(
                "Currently only 1D data structures may be returned from slice operations."
            )
        try:
            meta = self.meta[key]
        except (KeyError, TypeError):
            meta = self.meta

        return PHOENIXSpectrum(
            flux=flux,
            spectral_axis=self.spectral_axis[key],
            uncertainty=self.uncertainty[key] if self.uncertainty else None,
            wcs=self.wcs[key] if self.wcs else None,
            mask=self.mask[key] if self.mask else None,
            meta=meta,
        )

    grid_points = property(lambda self: self.meta["grid_points"])
    teff_points = property(lambda self: self.meta["teff_points"])
    metallicity_points = property(lambda self: self.meta["metallicity_points"])
    logg_points = property(lambda self: self.meta["logg_points"])
    grid_labels = property(lambda self: self.meta["grid_labels"])
    n_spectra = property(lambda self: self.meta["n_spectra"])
    lookup_dict = property(lambda self: self.meta["lookup_dict"])

    def truncate(self, wavelength_range=None, data=None):
        """Truncate the wavelength range of the grid

        Parameters
        ----------
        wavelength_range: list or tuple
            A pair of values that denote the shortest and longest wavelengths
            for truncating the grid.
        data: Spectrum1D-like
            A spectrum to which this method will match the wavelength limits
        
        Returns
        -------
        truncated_spectrum: Spectrum1D-like
            The spectrum after being truncated to the given wavelength range
        """
        fiducial_spectrum = deepcopy(self[0])
        wavelength_units = fiducial_spectrum.wavelength.unit
        flux_units = fiducial_spectrum.flux.unit

        if data and wavelength_range:
            wavelength_range = (
                fiducial_spectrum.wavelength.value.min() * wavelength_units,
                fiducial_spectrum.wavelength.value.max() * wavelength_units,
            )

        shortest_wavelength, longest_wavelength = wavelength_range

        wavelengths, fluxes = [], []
        for spectrum in self:
            mask = (spectrum.wavelength > shortest_wavelength) & (
                spectrum.wavelength < longest_wavelength
            )
            wavelengths.append(spectrum.wavelength.value[mask])
            fluxes.append(spectrum.flux.value[mask])

        assert fluxes and wavelengths

        return self.__class__(
            flux=np.array(fluxes) * flux_units,
            spectral_axis=np.array(wavelengths) * wavelength_units,
            meta=self.meta,
        )

    get_index = lambda self, grid_point: self.lookup_dict[grid_point]

    def find_nearest_teff(self, value):
        idx = np.abs(self.teff_points - value).argmin()
        return self.teff_points[idx]

    def find_nearest_metallicity(self, value):
        idx = np.abs(self.metallicity_points - value).argmin()
        return self.metallicity_points[idx]

    def show_dashboard(
        self, data=None, notebook_url="localhost:8888"
    ):  # pragma: no cover
        """Show an interactive dashboard for the PHOENIX grid;
        heavily inspired by the lightkurve .interact() method

        Parameters
        ----------
        data: Spectrum1D-like
            A normalized data spectrum over which to plot the models
        notebook_url: str
            Location of the Jupyter notebook page (default: "localhost:8888")
            When showing Bokeh applications, the Bokeh server must be
            explicitly configured to allow connections originating from
            different URLs. This parameter defaults to the standard notebook
            host and port. If you are running on a different location, you
            will need to supply this value for the application to display
            properly. If no protocol is supplied in the URL, e.g. if it is
            of the form "localhost:8888", then "http" will be used.
        """

        def create_interact_ui(doc):

            # Make the spectrum source
            scalar_norm = np.percentile(self[0].flux.value, 95)
            spec_source = ColumnDataSource(
                data={
                    "wavelength": self[0].wavelength.value,
                    "flux": self[0].flux.value / scalar_norm,
                    "native_flux": self[0].flux.value / scalar_norm,
                    "native_wavelength": self[0].wavelength.value,
                }
            )

            fig = figure(
                title="PHOENIX Interactive Dashboard",
                plot_height=500,
                plot_width=950,
                tools="pan,wheel_zoom,box_zoom,tap,reset",
                toolbar_location="below",
                border_fill_color="whitesmoke",
            )
            fig.title.offset = 350
            fig.yaxis.axis_label = "Flux"
            fig.xaxis.axis_label = "Wavelength (\u00B5m)"
            fig.y_range = Range1d(start=0, end=1.5)

            fig.step(
                "wavelength",
                "flux",
                line_width=1,
                color="red",
                source=spec_source,
                nonselection_line_color="red",
                nonselection_line_alpha=1.0,
            )
            wl_lo, wl_hi = (
                self[0].wavelength.value.min(),
                self[0].wavelength.value.max(),
            )

            if data:
                assert isinstance(
                    data, Spectrum1D
                ), "The data spectrum must be Spectrum1D-like"
                new_lo, new_hi = (
                    data.wavelength.value.min(),
                    data.wavelength.value.max(),
                )
                assert (new_lo < wl_hi) & (
                    new_hi > wl_lo
                ), "Data should overlap the models, double check your wavelength limits."
                wl_lo, wl_hi = new_lo, new_hi

                data_source = ColumnDataSource(
                    data={"wavelength": data.wavelength.value, "flux": data.flux.value}
                )
                fig.step(
                    "wavelength", "flux", line_width=1, color="blue", source=data_source
                )

            fig.x_range = Range1d(start=wl_lo, end=wl_hi)

            # Slider to decimate the data
            smoothing_slider = Slider(
                start=0.1,
                end=200,
                value=0.1,
                step=0.1,
                title="Rotational Broadening: v sin(i) [km/s]",
                width=700,
            )

            vz_slider = Slider(
                start=-200,
                end=200,
                value=0.00,
                step=0.05,
                title="Radial Velocity: RV [km/s]",
                width=700,
                format="0.000f",
            )

            teff_slider = Slider(
                start=min(self.teff_points),
                end=max(self.teff_points),
                value=self.teff_points[1],
                step=100,
                title="Effective Temperature: T_eff [K]",
                width=700,
            )
            teff_message = Div(
                text=f"Closest point: {self.teff_points[1]}K", width=100, height=10,
            )
            logg_slider = Slider(
                start=min(self.logg_points),
                end=max(self.logg_points),
                value=5.0,
                step=0.50,
                title="Surface Gravity: log(g) [cm/s^2]",
                width=700,
            )

            metallicity_slider = Slider(
                start=min(self.metallicity_points),
                end=max(self.metallicity_points),
                value=0.0,
                step=0.50,
                title="Metallicity: Z",
                width=700,
            )

            r_button = Button(label=">", button_type="default", width=32)
            l_button = Button(label="<", button_type="default", width=32)

            def update_upon_smooth(attr, old, new):
                """Callback to take action when smoothing slider changes"""
                new_spec = PHOENIXSpectrum(
                    spectral_axis=spec_source.data["native_wavelength"] * u.Angstrom,
                    flux=spec_source.data["native_flux"] * u.dimensionless_unscaled,
                ).rotationally_broaden(new)
                spec_source.data["flux"] = new_spec.flux.value

            def update_upon_vz(attr, old, new):
                """Callback to take action when vz slider changes"""
                new_spec = PHOENIXSpectrum(
                    spectral_axis=spec_source.data["native_wavelength"] * u.Angstrom,
                    flux=spec_source.data["native_flux"] * u.dimensionless_unscaled,
                ).rv_shift(new)
                spec_source.data["wavelength"] = new_spec.wavelength.value

            def update_upon_teff_selection(attr, old, new):
                """Callback to take action when teff slider changes"""
                teff = self.find_nearest_teff(new)
                if teff != old:
                    teff_message.text = f"Closest point: {teff}K"
                    point = (teff, logg_slider.value, metallicity_slider.value)
                    native_spec = self[self.get_index(point)].normalize(percentile=95)
                    new_spec = native_spec.rotationally_broaden(
                        smoothing_slider.value
                    ).rv_shift(vz_slider.value)

                    spec_source.data = {
                        "native_wavelength": native_spec.wavelength.value,
                        "native_flux": native_spec.flux.value,
                        "wavelength": new_spec.wavelength.value,
                        "flux": new_spec.flux.value,
                    }

            def update_upon_metallicity_selection(attr, old, new):
                """Callback to take action when teff slider changes"""
                metallicity = self.find_nearest_metallicity(new)
                if metallicity != old:
                    teff = self.find_nearest_teff(teff_slider.value)
                    point = (teff, logg_slider.value, metallicity)
                    native_spec = self[self.get_index(point)].normalize(percentile=95)
                    new_spec = native_spec.rotationally_broaden(
                        smoothing_slider.value
                    ).rv_shift(vz_slider.value)

                    spec_source.data = {
                        "native_wavelength": native_spec.wavelength.value,
                        "native_flux": native_spec.flux.value,
                        "wavelength": new_spec.wavelength.value,
                        "flux": new_spec.flux.value,
                    }

            def update_upon_logg_selection(attr, old, new):
                """Callback to take action when logg slider changes"""
                teff = self.find_nearest_teff(teff_slider.value)
                metallicity = self.find_nearest_metallicity(metallicity_slider.value)
                point = (teff, new, metallicity)
                native_spec = self[self.get_index(point)].normalize(percentile=95)
                new_spec = native_spec.rotationally_broaden(
                    smoothing_slider.value
                ).rv_shift(vz_slider.value)

                spec_source.data = {
                    "native_wavelength": native_spec.wavelength.value,
                    "native_flux": native_spec.flux.value,
                    "wavelength": new_spec.wavelength.value,
                    "flux": new_spec.flux.value,
                }

            def go_right_by_one():
                """Step forward by a single cadence"""
                new_index = np.abs(self.teff_points - teff_slider.value).argmin() + 1
                if new_index < len(self.teff_points):
                    teff_slider.value = self.teff_points[new_index]

            def go_left_by_one():
                """Step backward by a single cadence"""
                new_index = np.abs(self.teff_points - teff_slider.value).argmin() - 1
                if new_index >= 0:
                    teff_slider.value = self.teff_points[new_index]

            r_button.on_click(go_right_by_one)
            l_button.on_click(go_left_by_one)
            smoothing_slider.on_change("value", update_upon_smooth)
            vz_slider.on_change("value", update_upon_vz)
            teff_slider.on_change("value", update_upon_teff_selection)
            logg_slider.on_change("value", update_upon_logg_selection)
            metallicity_slider.on_change("value", update_upon_metallicity_selection)

            sp1, sp2, sp3, sp4 = (Spacer(width=w) for w in (5, 10, 20, 100))

            widgets_and_figures = layout(
                [fig],
                [l_button, sp1, r_button, sp2, teff_slider, sp3, teff_message],
                [sp4, logg_slider],
                [sp4, metallicity_slider],
                [sp4, smoothing_slider],
                [sp4, vz_slider],
            )
            doc.add_root(widgets_and_figures)

        output_notebook(verbose=False, hide_banner=True)
        return show(create_interact_ui, notebook_url=notebook_url)
