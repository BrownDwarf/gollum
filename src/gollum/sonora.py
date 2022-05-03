r"""
Sonora Spectrum
-----------------

A container for a single grid-point of the Sonora precomputed synthetic model spectrum of brown dwarfs and free-floating Gas Giant planets.  The spectrum is a vector with coordinates wavelength and flux :math:`F(\lambda)`.

SonoraSpectrum
###############
"""

import copy
import warnings
import logging
from gollum.precomputed_spectrum import PrecomputedSpectrum
from gollum.telluric import TelluricSpectrum
import numpy as np
import astropy
import pandas as pd
from astropy import units as u
from specutils import SpectrumCollection
import os
from specutils.spectra.spectrum1d import Spectrum1D
from tqdm import tqdm

from math import sqrt

from bokeh.io import show, output_notebook, push_notebook
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import Slider, Span, Range1d, Dropdown
from bokeh.layouts import layout, Spacer
from bokeh.models.widgets import Button, Div

log = logging.getLogger(__name__)

#  See Issue: https://github.com/astropy/specutils/issues/779
warnings.filterwarnings(
    "ignore", category=astropy.utils.exceptions.AstropyDeprecationWarning
)
# See Issue: https://github.com/astropy/specutils/issues/800
warnings.filterwarnings("ignore", category=RuntimeWarning)


class Sonora2017Spectrum(PrecomputedSpectrum):
    """
    A container for a single Sonora precomputed synthetic spectrum of a brown dwarfs or free-floating
    Gas Giant planet.

    Args:
        Teff (int): The Teff label of the Sonora model to read in.  Must be on the Sonora grid.
        logg (float): The logg label of the Sonora model to read in.  Must be on the Sonora grid.
        path (str): The path to your local Sonora grid library.  You must have the Sonora grid downloaded locally.  Default: "~/libraries/raw/Sonora/"
        wl_lo (float): the bluest wavelength of the models to keep (Angstroms)
        wl_hi (float): the reddest wavelength of the models to keep (Angstroms)
        """

    def __init__(
        self, *args, teff=None, logg=None, path=None, wl_lo=8038, wl_hi=12849, **kwargs
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

        if path is None:
            path = "~/libraries/raw/Sonora/"

        if (teff is not None) & (logg is not None):
            base_path = os.path.expanduser(path)
            assert os.path.exists(
                base_path
            ), "You must specify the path to local Sonora models: {}".format(base_path)

            print("successfully found the path to Sonora models: {}".format(base_path))

            assert teff in teff_points, "Teff must be on the grid points"
            assert logg in logg_points, "logg must be on the grid points"

            base_name = "sp_t{0:0>.0f}g{1:}nc_m0.0".format(
                float(teff), logg_par_dict[logg]
            )
            fn = base_path + "/" + base_name + ".gz"

            assert os.path.exists(fn), "Double check that the file {} exists".format(fn)

            # Units: micron, erg/cm^2/s/Hz
            df_native = (
                pd.read_csv(
                    fn,
                    skiprows=[0, 1],
                    delim_whitespace=True,
                    compression="gzip",
                    names=["wavelength_um", "flux"],
                )
                .sort_values("wavelength_um")
                .reset_index(drop=True)
            )

            # convert to Angstrom
            df_native["wavelength"] = df_native["wavelength_um"] * 10_000.0
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

    Args:
        Teff (int): The Teff label of the Sonora model to read in.  Must be on the Sonora grid.
        logg (float): The logg label of the Sonora model to read in.  Must be on the Sonora grid.
        path (str): The path to your local Sonora grid library.  You must have the Sonora grid downloaded locally.  Default: "~/libraries/raw/Sonora/"
        wl_lo (float): the bluest wavelength of the models to keep (Angstroms)
        wl_hi (float): the reddest wavelength of the models to keep (Angstroms)
        """

    def __init__(
        self,
        *args,
        teff=None,
        logg=None,
        metallicity=0.0,
        path=None,
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

        if path is None:
            path = "~/libraries/raw/SonoraBobcat2021/"

        if (teff is not None) & (logg is not None):
            base_path = os.path.expanduser(path)
            assert os.path.exists(
                base_path
            ), "You must specify the path to local Sonora models: {}".format(base_path)

            assert teff in teff_points, "Teff must be on the grid points"
            assert logg in logg_points, "logg must be on the grid points"
            assert metallicity in metallicity_points, "Fe/H must be a valid point"

            if metallicity <= 0:
                base_name = "sp_t{0:0>.0f}g{1:}nc_m{2:0.01f}".format(
                    float(teff), logg_par_dict[logg], float(metallicity)
                )
            else:
                base_name = "sp_t{0:0>.0f}g{1:}nc_m{2:+0.1f}".format(
                    float(teff), logg_par_dict[logg], float(metallicity)
                )
            fn = base_path + "/" + base_name

            assert os.path.exists(fn), "Double check that the file {} exists".format(fn)

            # Units: micron, erg/cm^2/s/Hz
            df_native = (
                pd.read_csv(
                    fn,
                    skiprows=[0, 1],
                    delim_whitespace=True,
                    names=["wavelength_um", "flux"],
                )
                .sort_values("wavelength_um")
                .reset_index(drop=True)
            )

            # convert to Angstrom
            df_native["wavelength"] = df_native["wavelength_um"] * 10_000.0
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
    r"""
    A container for a grid of Sonora precomputed synthetic spectra of brown dwarfs and free-floating
    Gas Giant planets.

    Args:
        Teff_range (tuple): The Teff limits of the grid model to read in.
        logg (tuple): The logg limits of the Sonora model to read in.
        metallicity_range (tuple): The metallicity limits of the Sonora model to read in
        path (str): The path to your local Sonora grid library.
            You must have the Sonora grid downloaded locally.
            Default: "~/libraries/raw/Sonora/"
        wl_lo (float): the bluest wavelength of the models to keep (Angstroms)
        wl_hi (float): the reddest wavelength of the models to keep (Angstroms)
        """

    def __init__(
        self,
        teff_range=None,
        logg_range=None,
        metallicity_range=None,
        path=None,
        wl_lo=8038,
        wl_hi=12849,
        **kwargs,
    ):

        if (
            ("flux" in kwargs.keys())
            and ("spectral_axis" in kwargs.keys())
            and ("meta" in kwargs.keys())
        ):
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

            if teff_range is not None:
                subset = (teff_points >= teff_range[0]) & (teff_points <= teff_range[1])
                teff_points = teff_points[subset]

            if logg_range is not None:
                subset = (logg_points >= logg_range[0]) & (logg_points <= logg_range[1])
                logg_points = logg_points[subset]

            if metallicity_range is not None:
                subset = (metallicity_points >= metallicity_range[0]) & (metallicity_points <= metallicity_range[1])
                metallicity_points = metallicity_points[subset]

            wavelengths, fluxes = [], []
            grid_points = []

            pbar = tqdm(teff_points)
            for teff in pbar:
                for logg in logg_points:
                    for metallicity in metallicity_points:
                        pbar.set_description(
                            "Processing Teff={} K, logg={:0.2f}, metallicity={:0.1f}".format(
                                teff, logg, metallicity
                            )
                        )
                        grid_point = (teff, logg, metallicity)
                        # "To do": See issue 31
                        # Temporary work around: pass if we can't read the file
                        try:
                            spec = SonoraSpectrum(
                                teff=teff,
                                logg=logg,
                                metallicity=metallicity,
                                path=path,
                                wl_lo=wl_lo,
                                wl_hi=wl_hi,
                            )
                            wavelengths.append(spec.wavelength)
                            fluxes.append(spec.flux)
                            grid_points.append(grid_point)
                        except:
                            log.info(
                                "Grid point Teff={} K, logg={:0.2f}, metallicity={:0.1f} does not exist".format(
                                    teff, logg, metallicity
                                )
                                + " in Sonora library"
                            )
            flux_out = np.array(fluxes) * fluxes[0].unit
            wave_out = np.array(wavelengths) * wavelengths[0].unit

            # Make a quick-access dictionary
            n_spectra = len(grid_points)
            lookup_dict = {grid_points[i]: i for i in range(n_spectra)}
            meta = {
                "teff_points": teff_points,
                "logg_points": logg_points,
                "metallicity_points": metallicity_points,
                "grid_labels": ("T_eff", "log(g)"),
                "n_spectra": n_spectra,
                "grid_points": grid_points,
                "lookup_dict": lookup_dict,
            }

            super().__init__(flux=flux_out, spectral_axis=wave_out, meta=meta)

    def __getitem__(self, key):
        flux = self.flux[key]
        if flux.ndim != 1:
            raise ValueError(
                "Currently only 1D data structures may be "
                "returned from slice operations."
            )
        spectral_axis = self.spectral_axis[key]
        uncertainty = None if self.uncertainty is None else self.uncertainty[key]
        mask = None if self.mask is None else self.mask[key]
        if self.meta is None:
            meta = None
        else:
            try:
                meta = self.meta[key]
            except KeyError:
                meta = self.meta

        return SonoraSpectrum(
            flux=flux,
            spectral_axis=spectral_axis,
            uncertainty=uncertainty,
            wcs=None,
            mask=mask,
            meta=meta,
        )

    @property
    def grid_points(self):
        """What are the coordinates of the grid?"""
        return self.meta["grid_points"]

    @property
    def teff_points(self):
        """What are the Teff points of the grid?"""
        return self.meta["teff_points"]

    @property
    def logg_points(self):
        """What are the logg points of the grid?"""
        return self.meta["logg_points"]

    @property
    def metallicity_points(self):
        """What are the metallicity points of the grid?"""
        return self.meta["metallicity_points"]

    @property
    def grid_labels(self):
        """What are the grid labels?"""
        return self.meta["grid_labels"]

    @property
    def n_spectra(self):
        """How many distinct spectra are in the grid?"""
        return self.meta["n_spectra"]

    @property
    def lookup_dict(self):
        """Lookup dictioary for spectra from their grid coordinates"""
        return self.meta["lookup_dict"]

    def truncate(self, wavelength_range=None, data=None):
        """Truncate the wavelength range of the grid

        Parameters
        ----------
        wavelength_range: List or Tuple of Quantities
            A pair of values that denote the shortest and longest wavelengths
            for truncating the grid.
        data: Spectrum1D-like
            A spectrum to which this method will match the wavelength limits

        """
        fiducial_spectrum = copy.deepcopy(self[0])
        wavelength_units = fiducial_spectrum.wavelength.unit
        flux_units = fiducial_spectrum.flux.unit

        if (data is not None) and (wavelength_range is None):
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

        fluxes = np.array(fluxes) * flux_units
        wavelengths = np.array(wavelengths) * wavelength_units
        assert fluxes is not None
        assert wavelengths is not None

        return self.__class__(flux=fluxes, spectral_axis=wavelengths, meta=self.meta)

    def get_index(self, grid_point):
        """Get the spectrum index associated with a given grid point
        """
        return self.lookup_dict[grid_point]

    def get_distance(self, gridpoint1, gridpoint2):
        return sqrt(((gridpoint1[0]-gridpoint2[0]) ** 2) + ((gridpoint1[1]-gridpoint2[1]) ** 2) + ((gridpoint1[2]-gridpoint2[2]) ** 2))

    #  Need to add a function to find the near grid point in the case it doesn't exist (find nearest point in a lattice)
    def find_nearest_grid_point(self, teff, logg, metallicity):
        current = (teff, logg, metallicity)

        distances = []
        for point in self.grid_points:
            distances.append(self.get_distance(current, point))
        shortest_distance = min(distances)
        idx = distances.index(shortest_distance)
        nearest_point = self.grid_points[idx]

        return nearest_point

    def find_nearest_teff(self, value):
        idx = (np.abs(self.teff_points - value)).argmin()
        return self.teff_points[idx]

    def show_dashboard(
        self, data=None, notebook_url="localhost:8888", show_telluric=True
    ):
        """Show an interactive dashboard for interacting with the Sonora grid
        Heavily inspired by the lightkurve .interact() method.

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
                data=dict(
                    wavelength=self[0].wavelength,
                    flux=self[0].flux.value / scalar_norm,
                    native_flux=self[0].flux.value / scalar_norm,
                    native_wavelength=self[0].wavelength.value,
                )
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

            instrumental_resolution = 100_000
            if data is not None:
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
                    metallicity_message.text = "Closest Metallicity point: {}".format(metallicity)
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
                    metallicity_message.text = "Closest Metallicity point: {}".format(metallicity)
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
                    metallicity_message.text = "Closest Metallicity point: {}".format(metallicity)
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

            sp1, sp2, sp3, sp4 = (
                Spacer(width=5),
                Spacer(width=10),
                Spacer(width=20),
                Spacer(width=100),
            )

            widgets_and_figures = layout(
                [fig],
                [l_button, sp1, r_button, sp2, teff_slider, sp3, teff_message],
                [sp4, logg_slider, sp3, logg_message],
                [sp4, metallicity_slider, sp3, metallicity_message], 
                [sp4, smoothing_slider],
                [sp4, vz_slider],
                [sp4, scale_slider],
            )
            doc.add_root(widgets_and_figures)

        output_notebook(verbose=False, hide_banner=True)
        show(create_interact_ui, notebook_url=notebook_url)
