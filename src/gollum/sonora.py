r"""
Sonora Spectrum
-----------------

A container for a single grid-point of the Sonora precomputed synthetic model spectrum of brown dwarfs and free-floating Gas Giant planets.  The spectrum is a vector with coordinates wavelength and flux :math:`F(\lambda)`.

SonoraSpectrum
###############
"""

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


class SonoraSpectrum(PrecomputedSpectrum):
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


class SonoraGrid(SpectrumCollection):
    r"""
    A container for a grid of Sonora precomputed synthetic spectra of brown dwarfs and free-floating 
    Gas Giant planets.  

    Args:
        Teff_range (tuple): The Teff limits of the grid model to read in.
        logg (tuple): The logg limits of the Sonora model to read in.
        path (str): The path to your local Sonora grid library.  
            You must have the Sonora grid downloaded locally.  
            Default: "~/libraries/raw/Sonora/"
        wl_lo (float): the bluest wavelength of the models to keep (Angstroms)
        wl_hi (float): the reddest wavelength of the models to keep (Angstroms)
        """

    def __init__(
        self, teff_range=None, logg_range=None, path=None, wl_lo=8038, wl_hi=12849,
    ):
        teff_points = np.hstack(
            (
                np.arange(500, 600, 25),
                np.arange(600, 1000, 50),
                np.arange(1000, 2401, 100),
            )
        )
        logg_points = np.arange(4.0, 5.51, 0.25)

        if teff_range is not None:
            subset = (teff_points >= teff_range[0]) & (teff_points <= teff_range[1])
            teff_points = teff_points[subset]

        if logg_range is not None:
            subset = (logg_points >= logg_range[0]) & (logg_points <= logg_range[1])
            logg_points = logg_points[subset]

        self.teff_points = teff_points
        self.logg_points = logg_points
        self.grid_labels = ("T_eff", "log(g)")

        wavelengths, fluxes = [], []
        grid_points = []

        pbar = tqdm(teff_points)
        for teff in pbar:
            for logg in logg_points:
                pbar.set_description(
                    "Processing Teff={} K, logg={:0.2f}".format(teff, logg)
                )
                grid_point = (teff, logg)
                spec = SonoraSpectrum(
                    teff=teff, logg=logg, path=path, wl_lo=wl_lo, wl_hi=wl_hi
                )
                wavelengths.append(spec.wavelength)
                fluxes.append(spec.flux)
                grid_points.append(grid_point)
        flux_out = np.array(fluxes) * fluxes[0].unit
        wave_out = np.array(wavelengths) * wavelengths[0].unit

        # Make a quick-access dictionary
        n_spectra = len(grid_points)
        self.n_spectra = n_spectra
        lookup_dict = {grid_points[i]: i for i in range(n_spectra)}
        self.lookup_dict = lookup_dict

        super().__init__(
            flux=flux_out, spectral_axis=wave_out, meta={"grid_points": grid_points}
        )

    def __getitem__(self, key):
        flux = self.flux[key]
        if flux.ndim != 1:
            raise ValueError(
                "Currently only 1D data structures may be "
                "returned from slice operations."
            )
        spectral_axis = self.spectral_axis[key]
        uncertainty = None if self.uncertainty is None else self.uncertainty[key]
        wcs = None if self.wcs is None else self.wcs[key]
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
            wcs=wcs,
            mask=mask,
            meta=meta,
        )

    @property
    def grid_points(self):
        """What are the grid points of the spectrum?"""
        return self.meta["grid_points"]

    def get_index(self, grid_point):
        """Get the spectrum index associated with a given grid point
        """
        return self.lookup_dict[grid_point]

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
                color="darkslateblue",
                source=spec_source,
                legend_label="Sonora Model",
                nonselection_line_color="darkslateblue",
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
                text="Closest grid point: {}".format(1000), width=100, height=10
            )
            logg_slider = Slider(
                start=min(self.logg_points),
                end=max(self.logg_points),
                value=5.0,
                step=0.25,
                title="Surface Gravity: log(g) [cm/s^2]",
                width=490,
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
                teff = self.find_nearest_teff(new)
                if teff != old:
                    teff_message.text = "Closest grid point: {}".format(teff)
                    logg = logg_slider.value
                    grid_point = (teff, logg)
                    index = self.get_index(grid_point)

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
                teff = self.find_nearest_teff(teff_slider.value)

                grid_point = (teff, new)
                index = self.get_index(grid_point)

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
                [sp4, logg_slider],
                [sp4, smoothing_slider],
                [sp4, vz_slider],
                [sp4, scale_slider],
            )
            doc.add_root(widgets_and_figures)

        output_notebook(verbose=False, hide_banner=True)
        show(create_interact_ui)
