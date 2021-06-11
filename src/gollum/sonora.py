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
import numpy as np
import astropy
import pandas as pd
from astropy import units as u
from specutils import SpectrumCollection
import os
from tqdm import tqdm

from bokeh.io import show, output_notebook, push_notebook
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import Slider, Span, Range1d, Dropdown
from bokeh.layouts import layout, Spacer
from bokeh.models.widgets import Button, Div

from scipy.ndimage import gaussian_filter1d
from collections import OrderedDict

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
        """What is the provenance of each spectrum?"""
        return self.meta["grid_points"]

    def get_index(self, grid_point):
        """Get the spectrum index associated with a given grid point
        """
        return self.lookup_dict[grid_point]

    def find_nearest_teff(self, value):
        idx = (np.abs(self.teff_points - value)).argmin()
        return self.teff_points[idx]

    def create_interact_ui(self, doc):

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
            title="Sonora Bobcat in Bokeh",
            plot_height=340,
            plot_width=600,
            tools="pan,wheel_zoom,box_zoom,tap,reset",
            toolbar_location="below",
            border_fill_color="whitesmoke",
        )
        fig.title.offset = -10
        fig.yaxis.axis_label = "Flux "
        fig.xaxis.axis_label = "Wavelength (micron)"
        fig.y_range = Range1d(start=0, end=1.5)
        xmin, xmax = (
            self.wavelength[0].value.min() * 0.995,
            self.wavelength[0].value.max() * 1.005,
        )
        fig.x_range = Range1d(start=xmin, end=xmax)

        fig.step(
            "wavelength",
            "flux",
            line_width=1,
            color="gray",
            source=spec_source,
            nonselection_line_color="gray",
            nonselection_line_alpha=1.0,
        )

        # Slider to decimate the data
        smoothing_slider = Slider(
            start=0.1,
            end=40,
            value=0.1,
            step=0.1,
            title="Spectral resolution kernel",
            width=490,
        )

        vz_slider = Slider(
            start=-0.009,
            end=0.009,
            value=0.00,
            step=0.0005,
            title="Radial Velocity",
            width=490,
            format="0.000f",
        )

        teff_slider = Slider(
            start=min(self.teff_points),
            end=max(self.teff_points),
            value=1000,
            step=25,
            title="Teff",
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
            title="logg",
            width=490,
        )
        r_button = Button(label=">", button_type="default", width=30)
        l_button = Button(label="<", button_type="default", width=30)

        def update_upon_smooth(attr, old, new):
            """Callback to take action when smoothing slider changes"""
            # spec_source.data["wavelength"] = df_nir.wavelength.values[::new]
            spec_source.data["flux"] = gaussian_filter1d(
                spec_source.data["native_flux"], new
            )

        def update_upon_vz(attr, old, new):
            """Callback to take action when vz slider changes"""
            spec_source.data["wavelength"] = (
                spec_source.data["native_wavelength"] - new * 10_000
            )
            # spec_source.data["flux"] = gaussian_filter1d(df_nir.flux.values, new)

        def update_upon_teff_selection(attr, old, new):
            """Callback to take action when teff slider changes"""
            teff = self.find_nearest_teff(new)
            if teff != old:
                teff_message.text = "Closest grid point: {}".format(teff)
                logg = logg_slider.value
                grid_point = (teff, logg)
                index = self.get_index(grid_point)
                spec = self[index]
                scalar_norm = np.percentile(spec.flux.value, 95)
                spec_source.data["native_wavelength"] = spec.wavelength.value
                spec_source.data["native_flux"] = spec.flux.value / scalar_norm
                spec_source.data["wavelength"] = spec.wavelength.value - vz_slider.value
                spec_source.data["flux"] = gaussian_filter1d(
                    spec.flux.value / scalar_norm, smoothing_slider.value
                )

            else:
                pass

        def update_upon_logg_selection(attr, old, new):
            """Callback to take action when logg slider changes"""
            teff = self.find_nearest_teff(teff_slider.value)

            logg = logg_slider.value
            grid_point = (teff, new)
            index = self.get_index(grid_point)
            spec = self[index]

            scalar_norm = np.percentile(spec.flux.value, 95)
            spec_source.data["native_wavelength"] = spec.wavelength.value
            spec_source.data["native_flux"] = spec.flux.value / scalar_norm
            spec_source.data["wavelength"] = spec.wavelength.value - vz_slider.value
            spec_source.data["flux"] = gaussian_filter1d(
                spec.flux.value / scalar_norm, smoothing_slider.value
            )

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
        )
        doc.add_root(widgets_and_figures)

    def show_dashboard(self):
        """Show an interactive dashboard for interacting with the models"""
        output_notebook(verbose=False, hide_banner=True)
        show(self.create_interact_ui)
        pass
