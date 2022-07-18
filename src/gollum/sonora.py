r"""
Sonora Spectrum
-----------------

A container for a single grid-point of the Sonora precomputed synthetic model spectrum of brown dwarfs and free-floating Gas Giant planets.  The spectrum is a vector with coordinates wavelength and flux :math:`F(\lambda)`.

SonoraSpectrum
###############
"""

import os
import numpy as np

from itertools import product
from warnings import filterwarnings
from logging import getLogger
from tqdm import tqdm
from gollum.utilities import _truncate
from gollum.precomputed_spectrum import PrecomputedSpectrum
from gollum.telluric import TelluricSpectrum
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy import units as u
from pandas import read_csv
from specutils import SpectrumCollection, Spectrum1D
from scipy.ndimage import gaussian_filter1d
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


class Sonora2017Spectrum(PrecomputedSpectrum):
    """
    A container for Sonora 2017 spectra

    Parameters
    ----------
    teff : int
        The teff label of the Sonora model to read in.  Must be on the Sonora grid.
    logg : float
        The logg label of the Sonora model to read in.  Must be on the Sonora grid.
    path : str
        The path to your local Sonora grid library.  You must have the Sonora grid downloaded locally.  Default: "~/libraries/raw/Sonora/"
    wl_lo : float
        The shortest wavelength of the models to keep (Angstroms)
    wl_hi : float
        The longest wavelength of the models to keep (Angstroms)
    """

    def __init__(
        self,
        *args,
        teff=None,
        logg=None,
        path="~/libraries/raw/Sonora/",
        wl_lo=8038,
        wl_hi=12849,
        **kwargs,
    ):

        teff_points = np.hstack(
            (
                np.arange(500, 600, 25),
                np.arange(600, 1000, 50),
                np.arange(1000, 2401, 100),
            )
        )
        logg_points = np.arange(4.0, 5.51, 0.25)

        # Map logg (cgs) to the gravity labels used in file names
        logg_par_dict = {
            4.0: "100",
            4.25: "178",
            4.5: "316",
            4.75: "562",
            5.0: "1000",
            5.25: "1780",
            5.5: "3160",
        }

        if teff and logg:
            base_path = os.path.expanduser(path)
            assert os.path.exists(base_path), "Given path does not exist."
            assert teff in teff_points, "teff must be a point on the grid"
            assert logg in logg_points, "logg must be a point on the grid"

            fn = f"{base_path}/sp_t{teff}g{logg_par_dict[logg]}nc_m0.0.gz"

            df_native = (
                read_csv(
                    fn,
                    skiprows=[0, 1],
                    delim_whitespace=True,
                    compression="gzip",
                    names=["wavelength_um", "flux"],
                )
                .sort_values("wavelength_um")
                .reset_index(drop=True)
            )

            # convert to Angstroms
            df_native["wavelength"] = df_native["wavelength_um"] * 10000.0
            mask = (df_native.wavelength > wl_lo) & (df_native.wavelength < wl_hi)
            df_trimmed = df_native[mask].reset_index(drop=True)

            super().__init__(
                spectral_axis=df_trimmed.wavelength.values * u.Angstrom,
                flux=df_trimmed.flux.values * u.erg / u.s / u.cm ** 2 / u.Hz,
                **kwargs,
            )

        else:
            super().__init__(*args, **kwargs)


class Sonora2021Spectrum(PrecomputedSpectrum):
    r"""
    A container for a single Sonora precomputed synthetic spectrum of a brown dwarfs or free-floating
    Gas Giant planet.

    Parameters
    ----------
    teff : int
        The teff label of the Sonora model to read in.  Must be on the Sonora grid.
    logg : float
        The logg label of the Sonora model to read in.  Must be on the Sonora grid.
    metallicity : float
        The metallicity label of the Sonora model to read in.  Must be on the Sonora grid.
    path : str
        The path to your local Sonora grid library.  You must have the Sonora grid downloaded locally.  Default: "~/libraries/raw/Sonora/"
    wl_lo : float
        The shortest wavelength of the models to keep (Angstroms)
    wl_hi : float
        The longest wavelength of the models to keep (Angstroms)
    """

    def __init__(
        self,
        *args,
        teff=None,
        logg=None,
        metallicity=0.0,
        path="~/libraries/raw/SonoraBobcat2021/",
        wl_lo=8038,
        wl_hi=12849,
        **kwargs,
    ):

        teff_points = np.hstack(
            (
                np.arange(500, 600, 25),
                np.arange(600, 1000, 50),
                np.arange(1000, 2401, 100),
            )
        )
        logg_points = np.arange(3.0, 5.51, 0.25)

        # Map logg (cgs) to the gravity labels used in file names
        logg_par_dict = {
            3.0: "10",
            3.25: "17",
            3.5: "31",
            3.75: "56",
            4.0: "100",
            4.25: "178",
            4.5: "316",
            4.75: "562",
            5.0: "1000",
            5.25: "1780",
            5.5: "3160",
        }

        metallicity_points = np.arange(-0.5, 0.51, 0.5)

        if teff and logg:
            base_path = os.path.expanduser(path)
            assert os.path.exists(base_path), "Given path does not exist."
            assert teff in teff_points, "teff must be a point on the grid"
            assert logg in logg_points, "logg must be a point on the grid"
            assert metallicity in metallicity_points, "Fe/H must be a point on the grid"

            Z_string = f"{metallicity:+0.1f}" if metallicity else "0.0.gz"
            fn = f"{base_path}sp_t{teff}g{logg_par_dict[logg]}nc_m{Z_string}"

            df_native = (
                read_csv(
                    fn,
                    skiprows=[0, 1],
                    delim_whitespace=True,
                    names=["wavelength_um", "flux"],
                )
                .sort_values("wavelength_um")
                .reset_index(drop=True)
            )

            # convert to Angstroms
            df_native["wavelength"] = df_native["wavelength_um"] * 10000.0
            mask = (df_native.wavelength > wl_lo) & (df_native.wavelength < wl_hi)
            df_trimmed = df_native[mask].reset_index(drop=True)

            super().__init__(
                spectral_axis=df_trimmed.wavelength.values * u.Angstrom,
                flux=df_trimmed.flux.values * u.erg / u.s / u.cm ** 2 / u.Hz,
                **kwargs,
            )

        else:
            super().__init__(*args, **kwargs)


SonoraSpectrum = Sonora2021Spectrum


class SonoraGrid(SpectrumCollection):
    """
    A container for a grid of Sonora precomputed synthetic spectra of brown dwarfs and free-floating
    Gas Giant planets.

    Parameters
    ----------
    teff_range : tuple
        The Teff limits of the grid model to read in.
    logg_range : tuple
        The logg limits of the Sonora model to read in.
    metallicity_range : tuple
        The metallicity limits of the Sonora model to read in
    path : str
        The path to your locally downloaded Sonora grid library. Default: "~/libraries/raw/SonoraBobcat2021/"
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
        path="~/libraries/raw/SonoraBobcat2021/",
        wl_lo=8038,
        wl_hi=12849,
        **kwargs,
    ):

        if set(("flux", "spectral_axis", "meta")).issubset(kwargs):
            super().__init__(**kwargs)
        else:
            teff_points = np.hstack(
                (
                    np.arange(500, 600, 25),
                    np.arange(600, 1000, 50),
                    np.arange(1000, 2401, 100),
                )
            )
            logg_points = np.arange(3.0, 5.51, 0.25)

            metallicity_points = np.arange(-0.5, 0.51, 0.5)

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
                pbar.desc = f"Processing Teff={teff}K|logg={logg:0.2f}|Z={Z:0.1f}"
                try:
                    spec = SonoraSpectrum(
                        teff=teff,
                        logg=logg,
                        metallicity=Z,
                        path=path,
                        wl_lo=wl_lo,
                        wl_hi=wl_hi,
                    )
                    wavelengths.append(spec.wavelength)
                    fluxes.append(spec.flux)
                    grid_points.append((teff, logg, Z))
                except FileNotFoundError:
                    log.info(f"No file for Teff={teff}K|logg={logg:0.2f}|Z={Z:0.1f}")
                    missing += 1

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

        return SonoraSpectrum(
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

    def instrumental_broaden(self, resolving_power):
        """Instrumental broaden the grid"""

        # Currently assumes they all have the same wavelength grid!
        angstroms_per_pixel = np.median(np.diff(self.wavelength[0, :].value))
        lam0 = np.median(self.wavelength[0, :].value)
        delta_lam = lam0 / resolving_power

        scale_factor = 2.355
        sigma = delta_lam / scale_factor / angstroms_per_pixel

        convolved_flux = (
            gaussian_filter1d(self.flux.value, sigma, axis=1) * self.flux.unit
        )
        return self.__class__(
            flux=convolved_flux, spectral_axis=self.wavelength, meta=self.meta
        )

    def decimate(self, decimation_fraction=0.1):
        """Decimate the grid"""

        ## Hmmm, this implementation may be brittle...
        fluxes = []
        wavelengths = []

        for spec in self:
            newspec = spec.resample_to_uniform_in_velocity(decimation_fraction)
            fluxes.append(newspec.flux.value)
            wavelengths.append(newspec.wavelength.value)

        output = SonoraGrid(
            flux=np.array(fluxes) * newspec.flux.unit,
            spectral_axis=np.array(wavelengths) * newspec.wavelength.unit,
            meta=self.meta,
        )
        return output

    truncate = _truncate
    get_index = lambda self, grid_point: self.lookup_dict[grid_point]

    #  Need to add a function to find the near grid point in the case it doesn't exist (find nearest point in a lattice)
    def find_nearest_grid_point(self, teff, logg, metallicity):
        current = np.array((teff, logg, metallicity))
        mindist = np.Inf
        for point in map(np.array, self.grid_points):
            if (current_dist := np.linalg.norm(current - point)) < mindist:
                mindist, minpoint = current_dist, point
        return tuple(minpoint)

    def find_nearest_teff(self, value):
        idx = (np.abs(self.teff_points - value)).argmin()
        return self.teff_points[idx]

    def show_dashboard(
        self, data=None, notebook_url="localhost:8888", show_telluric=True
    ):  # pragma: no cover
        """Show an interactive dashboard for the Sonora grid;
        heavily inspired by the lightkurve .interact() method.

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
            scalar_norm = np.percentile(self[0].flux.value, 95)
            spec_source = ColumnDataSource(
                data={
                    "wavelength": self[0].wavelength,
                    "flux": self[0].flux.value / scalar_norm,
                    "native_flux": self[0].flux.value / scalar_norm,
                    "native_wavelength": self[0].wavelength.value,
                }
            )

            fig = figure(
                title="Sonora Bobcat Interactive Dashboard",
                plot_height=340,
                plot_width=600,
                tools="pan,wheel_zoom,box_zoom,reset,save",
                toolbar_location="below",
                border_fill_color="whitesmoke",
            )
            fig.title.offset = -10
            fig.yaxis.axis_label = "Flux "
            fig.xaxis.axis_label = "Wavelength (Angstrom)"
            fig.y_range = Range1d(start=0, end=1.9)

            wl_lo, wl_hi = (
                self[0].wavelength.value.min(),
                self[0].wavelength.value.max(),
            )

            instrumental_resolution = 100000
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
                    data=dict(wavelength=data.wavelength.value, flux=data.flux.value,)
                )
                fig.step(
                    "wavelength",
                    "flux",
                    line_width=2,
                    legend_label="data",
                    source=data_source,
                )
                if hasattr(data, "instrumental_resolution"):
                    instrumental_resolution = data.instrumental_resolution

            fig.x_range = Range1d(start=wl_lo, end=wl_hi)

            if show_telluric:
                tell_spec = TelluricSpectrum(
                    path="default", wl_lo=wl_lo, wl_hi=wl_hi
                ).instrumental_broaden(instrumental_resolution * 2)
                tell_source = ColumnDataSource(
                    data=dict(
                        wavelength=tell_spec.wavelength.value,
                        flux=tell_spec.flux.value,
                    )
                )
                out_glyph = fig.step(
                    "wavelength",
                    "flux",
                    line_width=2,
                    color="#bdc3c7",
                    legend_label="Telluric",
                    source=tell_source,
                )
                out_glyph.level = "underlay"

            fig.step(
                "wavelength",
                "flux",
                line_width=2,
                color="Red",
                source=spec_source,
                legend_label="Sonora Model",
                nonselection_line_color="DarkOrange",
                nonselection_line_alpha=1.0,
            )

            fig.legend.location = "top_left"
            fig.legend.orientation = "horizontal"

            # Slider to decimate the data
            smoothing_slider = Slider(
                start=0.1,
                end=40,
                value=0.1,
                step=0.1,
                title="Rotational Broadening: v sin(i) [km/s]",
                width=490,
            )

            vz_slider = Slider(
                start=-200,
                end=200,
                value=0.00,
                step=0.05,
                title="Radial Velocity: RV [km/s]",
                width=490,
                format="0.000f",
            )

            teff_slider = Slider(
                start=min(self.teff_points),
                end=max(self.teff_points),
                value=1000,
                step=25,
                title="Effective Temperature: T_eff [Kelvin]",
                width=490,
            )
            teff_message = Div(
                text="Closest T_eff point: {}".format(1000), width=100, height=10
            )
            logg_slider = Slider(
                start=min(self.logg_points),
                end=max(self.logg_points),
                value=5.0,
                step=0.25,
                title="Surface Gravity: log(g) [cm/s^2]",
                width=490,
            )
            logg_message = Div(
                text="Closest log(g) point: {}".format(1000), width=100, height=10
            )
            metallicity_slider = Slider(
                start=min(self.metallicity_points),
                end=max(self.metallicity_points),
                value=0.0,
                step=0.5,
                title="Metallicity: Metallicity [Fe/H]",
                width=490,
            )
            metallicity_message = Div(
                text="Closest Metallicity point: {}".format(1000), width=100, height=10
            )
            scale_slider = Slider(
                start=0.1,
                end=2.0,
                value=1.0,
                step=0.005,
                title="Normalization Scalar",
                width=490,
            )
            r_button = Button(label=">", button_type="default", width=30)
            l_button = Button(label="<", button_type="default", width=30)

            def update_upon_scale(attr, old, new):
                """Callback to take action when smoothing slider changes"""
                new_spec = (
                    SonoraSpectrum(
                        spectral_axis=spec_source.data["native_wavelength"]
                        * u.Angstrom,
                        flux=spec_source.data["native_flux"] * u.dimensionless_unscaled,
                    )
                    .rotationally_broaden(smoothing_slider.value)
                    .multiply(new * u.dimensionless_unscaled)
                    .rv_shift(vz_slider.value)
                )
                spec_source.data["flux"] = new_spec.flux.value

            def update_upon_smooth(attr, old, new):
                """Callback to take action when smoothing slider changes"""
                new_spec = (
                    SonoraSpectrum(
                        spectral_axis=spec_source.data["native_wavelength"]
                        * u.Angstrom,
                        flux=spec_source.data["native_flux"] * u.dimensionless_unscaled,
                    )
                    .rotationally_broaden(new)
                    .multiply(scale_slider.value * u.dimensionless_unscaled)
                    .rv_shift(vz_slider.value)
                )
                spec_source.data["flux"] = new_spec.flux.value

            def update_upon_vz(attr, old, new):
                """Callback to take action when vz slider changes"""
                new_spec = SonoraSpectrum(
                    spectral_axis=spec_source.data["native_wavelength"] * u.Angstrom,
                    flux=spec_source.data["native_flux"] * u.dimensionless_unscaled,
                ).rv_shift(new)
                spec_source.data["wavelength"] = new_spec.wavelength.value

            def update_upon_teff_selection(attr, old, new):
                """Callback to take action when teff slider changes"""
                logg = logg_slider.value
                metallicity = metallicity_slider.value
                new_grid_point = self.find_nearest_grid_point(new, logg, metallicity)

                if new_grid_point[0] != old:
                    teff = new_grid_point[0]
                    logg = new_grid_point[1]
                    metallicity = new_grid_point[2]

                    teff_message.text = "Closest T_eff point: {}".format(teff)
                    logg_message.text = "Closest log(g) point: {}".format(logg)
                    metallicity_message.text = "Closest Metallicity point: {}".format(
                        metallicity
                    )
                    index = self.get_index(new_grid_point)

                    native_spec = self[index].normalize(percentile=95)
                    new_spec = (
                        native_spec.rotationally_broaden(smoothing_slider.value)
                        .multiply(scale_slider.value * u.dimensionless_unscaled)
                        .rv_shift(vz_slider.value)
                    )

                    spec_source.data = {
                        "native_wavelength": native_spec.wavelength.value,
                        "native_flux": native_spec.flux.value,
                        "wavelength": new_spec.wavelength.value,
                        "flux": new_spec.flux.value,
                    }
                else:
                    pass

            def update_upon_logg_selection(attr, old, new):
                """Callback to take action when logg slider changes"""
                teff = teff_slider.value
                metallicity = metallicity_slider.value
                new_grid_point = self.find_nearest_grid_point(teff, new, metallicity)

                if new_grid_point[1] != old:
                    teff = new_grid_point[0]
                    logg = new_grid_point[1]
                    metallicity = new_grid_point[2]

                    teff_message.text = "Closest T_eff point: {}".format(teff)
                    logg_message.text = "Closest log(g) point: {}".format(logg)
                    metallicity_message.text = "Closest Metallicity point: {}".format(
                        metallicity
                    )
                    index = self.get_index(new_grid_point)

                    native_spec = self[index].normalize(percentile=95)
                    new_spec = (
                        native_spec.rotationally_broaden(smoothing_slider.value)
                        .multiply(scale_slider.value * u.dimensionless_unscaled)
                        .rv_shift(vz_slider.value)
                    )

                    spec_source.data = {
                        "native_wavelength": native_spec.wavelength.value,
                        "native_flux": native_spec.flux.value,
                        "wavelength": new_spec.wavelength.value,
                        "flux": new_spec.flux.value,
                    }

                else:
                    pass

            def update_upon_metallicity_selection(attr, old, new):
                """Callback to take action when metallicity slider changes"""
                teff = teff_slider.value
                logg = logg_slider.value
                new_grid_point = self.find_nearest_grid_point(teff, logg, new)

                if new_grid_point[2] != old:
                    teff = new_grid_point[0]
                    logg = new_grid_point[1]
                    metallicity = new_grid_point[2]

                    teff_message.text = "Closest T_eff point: {}".format(teff)
                    logg_message.text = "Closest log(g) point: {}".format(logg)
                    metallicity_message.text = "Closest Metallicity point: {}".format(
                        metallicity
                    )
                    index = self.get_index(new_grid_point)

                    native_spec = self[index].normalize(percentile=95)
                    new_spec = (
                        native_spec.rotationally_broaden(smoothing_slider.value)
                        .multiply(scale_slider.value * u.dimensionless_unscaled)
                        .rv_shift(vz_slider.value)
                    )

                    spec_source.data = {
                        "native_wavelength": native_spec.wavelength.value,
                        "native_flux": native_spec.flux.value,
                        "wavelength": new_spec.wavelength.value,
                        "flux": new_spec.flux.value,
                    }
                else:
                    pass

            def go_right_by_one():
                """Step forward in time by a single cadence"""
                current_index = np.abs(self.teff_points - teff_slider.value).argmin()
                new_index = current_index + 1
                if new_index <= (len(self.teff_points) - 1):
                    teff_slider.value = self.teff_points[new_index]

            def go_left_by_one():
                """Step back in time by a single cadence"""
                current_index = np.abs(self.teff_points - teff_slider.value).argmin()
                new_index = current_index - 1
                if new_index >= 0:
                    teff_slider.value = self.teff_points[new_index]

            r_button.on_click(go_right_by_one)
            l_button.on_click(go_left_by_one)
            smoothing_slider.on_change("value", update_upon_smooth)
            vz_slider.on_change("value", update_upon_vz)
            teff_slider.on_change("value", update_upon_teff_selection)
            logg_slider.on_change("value", update_upon_logg_selection)
            metallicity_slider.on_change("value", update_upon_metallicity_selection)
            scale_slider.on_change("value", update_upon_scale)

            sp1, sp2, sp3, sp4, sp5 = (
                Spacer(width=5),
                Spacer(width=10),
                Spacer(width=20),
                Spacer(width=100),
                Spacer(width=25),
            )

            widgets_and_figures = layout(
                [fig],
                [l_button, sp1, r_button, sp2, teff_slider, sp5, teff_message],
                [sp4, logg_slider, sp3, logg_message],
                [sp4, metallicity_slider, sp3, metallicity_message],
                [sp4, smoothing_slider],
                [sp4, vz_slider],
                [sp4, scale_slider],
            )
            doc.add_root(widgets_and_figures)

        output_notebook(verbose=False, hide_banner=True)
        show(create_interact_ui, notebook_url=notebook_url)
