r"""
PHOENIX Spectrum
-----------------

A container for a single Phoenix grid-point spectrum of wavelength and flux :math:`F(\lambda)`.

PHOENIXSpectrum
###############
"""

import os

from itertools import product
from tqdm import tqdm
from urllib.error import URLError
from contextlib import suppress
from numpy.ma import compressed, masked_outside
from gollum.utilities import _truncate
from gollum.precomputed_spectrum import *
from astropy.utils.exceptions import AstropyWarning
from astropy.io import fits
from specutils import SpectrumCollection
from bokeh.io import show, output_notebook
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import Slider, Range1d, Toggle
from bokeh.layouts import layout, Spacer

log = getLogger(__name__)

#  See Issue: https://github.com/astropy/specutils/issues/779
filterwarnings("ignore", category=AstropyDeprecationWarning)
filterwarnings("ignore", category=AstropyWarning)
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
    Z : float
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
        Z=0.0,  # solar by default
        path="~/libraries/raw/PHOENIX/",
        download=False,
        wl_lo=8038,
        wl_hi=12849,
        **kwargs,
    ):
        if not (teff or logg):
            super().__init__(*args, **kwargs)
            return

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

        Z_string = f"{Z:+0.1f}" if Z else "-0.0"

        wl_orig = fits.open(wl_file)[0].data.astype(np.float64)
        mask = (wl_orig >= wl_lo) & (wl_orig <= wl_hi)
        wl_out = wl_orig[mask]

        fn = f"{base_path}/Z{Z_string}/lte{teff:05d}-{logg:0.2f}{Z_string}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
        flux_orig = fits.open(fn)[0].data.astype(np.float64)
        flux_native = flux_orig[mask]
        native_flux_unit = u.erg / u.s / u.cm ** 2 / u.cm

        super().__init__(
            spectral_axis=wl_out * u.AA,
            flux=flux_native * native_flux_unit,
            meta={
                "teff": teff,
                "logg": logg,
                "Z": Z,
                "native_flux_unit": native_flux_unit,
            },
            **kwargs,
        )

    teff = property(lambda self: self.meta.get("teff"))
    logg = property(lambda self: self.meta.get("logg"))
    Z = property(lambda self: self.meta.get("Z"))


class PHOENIXGrid(SpectrumCollection):
    """
    A container for a grid of PHOENIX precomputed synthetic spectra of stars.

    Parameters
    ----------
    teff_range : tuple
        The teff limits of the grid model to read in.
    logg_range : tuple
        The logg limits of the grid model to read in.
    Z_range : tuple
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
        Z_range=None,
        path="~/libraries/raw/PHOENIX/",
        wl_lo=8038,
        wl_hi=12849,
        download=False,
        experimental=False,
        **kwargs,
    ):
        if set(("flux", "spectral_axis", "meta")).issubset(kwargs):
            super().__init__(**kwargs)
        else:

            teff_points = np.hstack(
                (np.arange(2300, 7000, 100), np.arange(7000, 12001, 200))
            )
            # Todo: some T_eff ranges go to log(g) = 0.0, consider adding these
            logg_points = np.arange(2.0, 6.01, 0.5)

            Z_points = np.array([-4, -3, -2, -1.5, -1, -0.5, 0, 0.5, 1])

            if teff_range:
                teff_points = compressed(masked_outside(teff_points, *teff_range))

            if logg_range:
                logg_points = compressed(masked_outside(logg_points, *logg_range))

            if Z_range:
                Z_points = compressed(masked_outside(Z_points, *Z_range))

            wavelengths, fluxes, grid_points = [], [], []
            iterlen = len(teff_points) * len(logg_points) * len(Z_points)
            pbar = tqdm(product(teff_points, logg_points, Z_points), total=iterlen)

            for teff, logg, Z in pbar:
                pbar.desc = f"Processing Teff={teff}K|log(g)={logg:0.2f}|Z={Z:+0.1f}"
                with suppress(FileNotFoundError, URLError):
                    spec = PHOENIXSpectrum(
                        teff=teff,
                        logg=logg,
                        Z=Z,
                        path=path,
                        wl_lo=wl_lo,
                        wl_hi=wl_hi,
                        download=download,
                    ).normalize(95)
                    wavelengths.append(spec.wavelength)
                    fluxes.append(spec.flux)
                    grid_points.append((teff, logg, Z))

            assert grid_points != [], "Empty grid; parameter limits out of range"
            super().__init__(
                flux=np.array(fluxes) * fluxes[0].unit,
                spectral_axis=np.array(wavelengths) * wavelengths[0].unit,
                meta={
                    "teff_points": teff_points,
                    "logg_points": logg_points,
                    "Z_points": Z_points,
                    "grid_labels": ("T_eff", "log(g)", "Z"),
                    "n_spectra": len(grid_points),
                    "grid_points": grid_points,
                    "lookup_dict": {value: i for i, value in enumerate(grid_points)},
                },
            )
            if experimental:
                from gollum.experimental import ExpPHOENIXGrid

                self.__class__ = ExpPHOENIXGrid

    def __getitem__(self, key):
        flux = self.flux[key]
        if flux.ndim != 1:
            raise ValueError(
                "Currently only 1D data structures may be returned from slice operations."
            )
        meta = self.meta.get(key, self.meta)

        meta["teff"], meta["logg"], meta["Z"] = self.grid_points[key]
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
    Z_points = property(lambda self: self.meta["Z_points"])
    logg_points = property(lambda self: self.meta["logg_points"])
    grid_labels = property(lambda self: self.meta["grid_labels"])
    n_spectra = property(lambda self: self.meta["n_spectra"])
    lookup_dict = property(lambda self: self.meta["lookup_dict"])

    truncate = _truncate
    get_index = lambda self, grid_point: self.lookup_dict[grid_point]

    def find_nearest_teff(self, value):
        idx = np.abs(self.teff_points - value).argmin()
        return self.teff_points[idx]

    def find_nearest_Z(self, value):
        idx = np.abs(self.Z_points - value).argmin()
        return self.Z_points[idx]

    def show_dashboard(self, data=None, url="localhost:8888"):  # pragma: no cover
        """Show an interactive dashboard for the PHOENIX grid;
        heavily inspired by the lightkurve .interact() method.

        If data is used, we recommend that the grid first be truncated to it,
        with a margin of 50 Angstroms on either end of the spectral axis to allow for
        radial velocity shift and rotational broadening to operate properly.

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
            spec_source = ColumnDataSource(
                data={
                    "wavelength": self[0].wavelength.value,
                    "flux": self[0].flux.value,
                    "native_flux": self[0].flux.value,
                    "native_wavelength": self[0].wavelength.value,
                }
            )
            wl_lo, wl_hi = (
                self[0].wavelength.value.min(),
                self[0].wavelength.value.max(),
            )
            fig = figure(
                title="PHOENIX Interactive Dashboard",
                width=950,
                height=500,
                tools="pan,wheel_zoom,box_zoom,tap,reset",
                toolbar_location="below",
                border_fill_color="whitesmoke",
            )

            if data:
                assert isinstance(
                    data, Spectrum1D
                ), "The data spectrum must be Spectrum1D-like"
                new_lo, new_hi = (
                    data.wavelength.value.min(),
                    data.wavelength.value.max(),
                )
                assert (
                    wl_lo < new_lo < new_hi < wl_hi
                ), "Data should overlap the models, double check your wavelength limits."
                wl_lo, wl_hi = new_lo, new_hi

                fig.step(
                    "wavelength",
                    "flux",
                    line_width=1,
                    color="black",
                    legend_label=data.meta["header"]["OBJECT"],
                    source=ColumnDataSource(
                        data={
                            "wavelength": data.wavelength.value,
                            "flux": data.flux.value,
                        }
                    ),
                )

            fig.title.align = "center"
            fig.title.text_font_size = "16pt"
            fig.yaxis.axis_label = "Normalized Flux"
            fig.xaxis.axis_label = "Wavelength (\u212B)"
            fig.axis.axis_label_text_font_style = "bold"
            fig.x_range = Range1d(start=wl_lo, end=wl_hi)
            fig.y_range = Range1d(start=0, end=1.5)
            fig.legend.location = "top_right"
            fig.legend.click_policy = "hide"
            fig.step(
                "wavelength",
                "flux",
                line_width=1,
                color="crimson",
                source=spec_source,
                nonselection_line_color="red",
                nonselection_line_alpha=1.0,
                legend_label="PHOENIX Model",
            )

            smoothing_slider = Slider(
                start=0.1,
                end=200,
                value=0.1,
                step=0.1,
                title="Rotational Broadening: v sin(i) [km/s]",
                width=460,
                bar_color="blue",
            )
            rv_slider = Slider(
                start=-200,
                end=200,
                value=0.00,
                step=0.05,
                title="Radial Velocity: RV [km/s]",
                width=460,
                format="0.000f",
                bar_color="blue",
            )
            teff_slider = Slider(
                start=min(self.teff_points),
                end=max(self.teff_points),
                value=min(self.teff_points),
                step=100,
                title="Effective Temperature: T_eff [K]",
                width=460,
                bar_color="red",
            )
            logg_slider = Slider(
                start=min(self.logg_points),
                end=max(self.logg_points),
                value=min(self.logg_points),
                step=0.50,
                title="Surface Gravity: log(g) [cm/s^2]",
                width=460,
                bar_color="red",
            )
            Z_slider = Slider(
                start=min(self.Z_points),
                end=max(self.Z_points),
                value=min(self.Z_points),
                step=0.50,
                title="Metallicity: Z",
                width=460,
                bar_color="red",
            )
            scale_slider = Slider(
                start=0.1,
                end=2.0,
                value=1.0,
                step=0.005,
                title="Normalization Scalar",
                width=460,
                bar_color="black",
            )
            continuum_toggle = Toggle(
                label="Fit Continuum (disables normalization)", button_type="success"
            )

            def update_to_continuum(active):
                """Callback to take action when the continuum toggle is toggled"""
                if active:
                    new_spec = PHOENIXSpectrum(
                        spectral_axis=spec_source.data["wavelength"] * u.AA,
                        flux=spec_source.data["flux"] * DV,
                    ).tilt_to_data(data)
                    scale_slider.disabled = True
                    continuum_toggle.label = "Undo Continuum (enables normalization)"
                else:
                    new_spec = (
                        PHOENIXSpectrum(
                            spectral_axis=spec_source.data["native_wavelength"] * u.AA,
                            flux=spec_source.data["native_flux"] * DV,
                        )
                        .rotationally_broaden(smoothing_slider.value)
                        .multiply(scale_slider.value * DV)
                        .rv_shift(rv_slider.value)
                    )
                    scale_slider.disabled = False
                    continuum_toggle.label = "Fit Continuum (disables normalization)"

                spec_source.data["flux"] = new_spec.flux.value

            def update_upon_scale(attr, old, new):
                """Callback to take action when normalization slider changes"""
                new_spec = (
                    PHOENIXSpectrum(
                        spectral_axis=spec_source.data["native_wavelength"] * u.AA,
                        flux=spec_source.data["native_flux"] * DV,
                    )
                    .rotationally_broaden(smoothing_slider.value)
                    .multiply(new * DV)
                    .rv_shift(rv_slider.value)
                )
                spec_source.data["flux"] = new_spec.flux.value

            def update_upon_smooth(attr, old, new):
                """Callback to take action when smoothing slider changes"""
                new_spec = PHOENIXSpectrum(
                    spectral_axis=spec_source.data["native_wavelength"] * u.AA,
                    flux=spec_source.data["native_flux"] * DV,
                ).rotationally_broaden(new)
                new_spec = (
                    new_spec.tilt_to_data(data)
                    if continuum_toggle.active
                    else new_spec.multiply(scale_slider.value * DV)
                )

                spec_source.data["flux"] = new_spec.flux.value

            def update_upon_rv(attr, old, new):
                """Callback to take action when RV slider changes"""
                new_spec = PHOENIXSpectrum(
                    spectral_axis=spec_source.data["native_wavelength"] * u.AA,
                    flux=spec_source.data["native_flux"] * DV,
                ).rv_shift(new)
                spec_source.data["wavelength"] = new_spec.wavelength.value

            def update_upon_teff_selection(attr, old, new):
                """Callback to take action when teff slider changes"""
                teff = self.find_nearest_teff(new)
                if teff != old:
                    point = (teff, logg_slider.value, Z_slider.value)
                    native_spec = self[self.get_index(point)]
                    new_spec = native_spec.rotationally_broaden(
                        smoothing_slider.value
                    ).rv_shift(rv_slider.value)

                    new_spec = (
                        new_spec.tilt_to_data(data)
                        if continuum_toggle.active
                        else new_spec.multiply(scale_slider.value * DV)
                    )

                    spec_source.data = {
                        "native_wavelength": native_spec.wavelength.value,
                        "native_flux": native_spec.flux.value,
                        "wavelength": new_spec.wavelength.value,
                        "flux": new_spec.flux.value,
                    }
                teff_slider.value = teff

            def update_upon_Z_selection(attr, old, new):
                """Callback to take action when metallicity slider changes"""
                Z = self.find_nearest_Z(new)
                if Z != old:
                    point = (teff_slider.value, logg_slider.value, Z)
                    native_spec = self[self.get_index(point)]
                    new_spec = native_spec.rotationally_broaden(
                        smoothing_slider.value
                    ).rv_shift(rv_slider.value)

                    new_spec = (
                        new_spec.tilt_to_data(data)
                        if continuum_toggle.active
                        else new_spec.multiply(scale_slider.value * DV)
                    )

                    spec_source.data = {
                        "native_wavelength": native_spec.wavelength.value,
                        "native_flux": native_spec.flux.value,
                        "wavelength": new_spec.wavelength.value,
                        "flux": new_spec.flux.value,
                    }

            def update_upon_logg_selection(attr, old, new):
                """Callback to take action when logg slider changes"""
                Z = self.find_nearest_Z(Z_slider.value)
                point = (teff_slider.value, new, Z)
                native_spec = self[self.get_index(point)]
                new_spec = native_spec.rotationally_broaden(
                    smoothing_slider.value
                ).rv_shift(rv_slider.value)

                new_spec = (
                    new_spec.tilt_to_data(data)
                    if continuum_toggle.active
                    else new_spec.multiply(scale_slider.value * DV)
                )

                spec_source.data = {
                    "native_wavelength": native_spec.wavelength.value,
                    "native_flux": native_spec.flux.value,
                    "wavelength": new_spec.wavelength.value,
                    "flux": new_spec.flux.value,
                }

            continuum_toggle.on_click(update_to_continuum)
            smoothing_slider.on_change("value", update_upon_smooth)
            rv_slider.on_change("value", update_upon_rv)
            teff_slider.on_change("value", update_upon_teff_selection)
            logg_slider.on_change("value", update_upon_logg_selection)
            Z_slider.on_change("value", update_upon_Z_selection)
            scale_slider.on_change("value", update_upon_scale)

            sp = Spacer(width=20)
            doc.add_root(
                layout(
                    [fig],
                    [continuum_toggle],
                    [teff_slider, sp, smoothing_slider],
                    [logg_slider, sp, rv_slider],
                    [Z_slider, sp, scale_slider],
                    background="whitesmoke",
                )
            )

        output_notebook(verbose=False, hide_banner=True)
        return show(create_interact_ui, notebook_url=url)
