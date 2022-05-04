r"""
PHOENIX Spectrum
-----------------

A container for a single Phoenix grid-point spectrum of wavelength and flux :math:`F(\lambda)`.

PHOENIXSpectrum
###############
"""

import copy
import warnings
import logging
from gollum.precomputed_spectrum import PrecomputedSpectrum
import numpy as np
import astropy
from astropy.io import fits
from astropy import units as u
from specutils import SpectrumCollection
from specutils.spectra.spectrum1d import Spectrum1D
from tqdm import tqdm
import os

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


class PHOENIXSpectrum(PrecomputedSpectrum):
    r"""
    A container for PHOENIX spectra

    Args:
        Teff (int): The Teff label of the PHOENIX model to read in.  Must be on the PHOENIX grid.
        logg (float): The logg label of the PHOENIX model to read in.  Must be on the PHOENIX grid.
        path (str): The path to your local PHOENIX grid library.  You must have the PHOENIX
            grid downloaded locally.  Default: "~/libraries/raw/PHOENIX/"
        download (bool): **Experimental** Whether or not you want to download the spectra 
            from the internet.  Requires an internet connection to work.
        wl_lo (float): the bluest wavelength of the models to keep (Angstroms)
        wl_hi (float): the reddest wavelength of the models to keep (Angstroms)
    """

    def __init__(
        self,
        *args,
        teff=None,
        logg=None,
        metallicity=None,
        path=None,
        download=False,
        wl_lo=8038,
        wl_hi=12849,
        **kwargs,
    ):

        if (teff is not None) & (logg is not None):

            if metallicity is None:
                metallicity = 0.0  # solar by default

            if path is None:
                path = "~/libraries/raw/PHOENIX/"

            if download == False:
                base_path = os.path.expanduser(path)
                assert os.path.exists(
                    base_path
                ), "You must specify the path to local PHOENIX models"

                wl_filename = base_path + "/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
                assert os.path.exists(
                    wl_filename
                ), f"You need to place the PHOENIX models in {base_path}"
            else:
                log.info(
                    "Experimental feature! Attempting to download PHOENIX models from the internet..."
                )
                log.info(
                    "We are using this FTP site: ftp://phoenix.astro.physik.uni-goettingen.de/v2.0/HiResFITS/"
                )
                wl_filename = "ftp://phoenix.astro.physik.uni-goettingen.de/v2.0/HiResFITS/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
                base_path = "ftp://phoenix.astro.physik.uni-goettingen.de/v2.0/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/"

            wl_orig = fits.open(wl_filename)[0].data.astype(np.float64)

            mask = (wl_orig > wl_lo) & (wl_orig < wl_hi)
            wl_out = wl_orig[mask]

            # Deal with metallicity
            metallicity_string = f"{metallicity:+0.1f}"

            metallicity_string = "-0.0" if metallicity == 0.0 else metallicity_string

            fn = (
                base_path
                + "/Z{}/lte{:05d}-{:0.2f}{}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
            ).format(metallicity_string, teff, logg, metallicity_string)

            flux_orig = fits.open(fn)[0].data.astype(np.float64)
            # Units: erg/s/cm^2/cm
            flux_native = flux_orig[mask]

            native_flux_units = u.erg / u.s / u.cm ** 2 / u.cm

            meta_dict = {
                "teff": teff,
                "logg": logg,
                "metallicity": metallicity,
                "native_flux_unit": native_flux_units,
            }

            super().__init__(
                spectral_axis=wl_out * u.AA,
                flux=flux_native * native_flux_units,
                meta=meta_dict,
                **kwargs,
            )

        else:
            super().__init__(*args, **kwargs)

    @property
    def teff(self):
        """The input Effective Temperature associated with this model"""
        if "teff" in self.meta.keys():
            return self.meta["teff"]
        else:
            return None

    @property
    def logg(self):
        """The input surface gravity associated with this model"""
        if "logg" in self.meta.keys():
            return self.meta["logg"]
        else:
            return None

    @property
    def metallicity(self):
        """The input metallicity associated with this model"""
        if "metallicity" in self.meta.keys():
            return self.meta["metallicity"]
        else:
            return None


class PHOENIXGrid(SpectrumCollection):
    r"""
    A container for a grid of PHOENIX precomputed synthetic spectra of stars.

    Args:
        Teff_range (tuple): The Teff limits of the grid model to read in.
        logg (tuple): The logg limits of the grid model to read in.
        path (str): The path to your local PHOENIX grid library.
            You must have the PHOENIX grid downloaded locally.
            Default: "~/libraries/raw/PHOENIX/"
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
            # Trigger a passthrough
            super().__init__(**kwargs)
        else:

            teff_points = np.hstack(
                (np.arange(2300, 7000, 100), np.arange(7000, 12_001, 200))
            )
            # Todo: some T_eff ranges go to log(g) = 0.0, consider adding these
            logg_points = np.arange(2.0, 6.01, 0.5)

            metallicity_points = np.array([-4, -3, -2, -1.5, -1, -0.5, 0, 0.5, 1])

            if teff_range is not None:
                subset = (teff_points >= teff_range[0]) & (teff_points <= teff_range[1])
                teff_points = teff_points[subset]

            if logg_range is not None:
                subset = (logg_points >= logg_range[0]) & (logg_points <= logg_range[1])
                logg_points = logg_points[subset]

            if metallicity_range is not None:
                subset = (metallicity_points >= metallicity_range[0]) & (
                    metallicity_points <= metallicity_range[1]
                )
                metallicity_points = metallicity_points[subset]

            wavelengths, fluxes = [], []
            grid_points = []

            pbar = tqdm(teff_points)
            for teff in pbar:
                for logg in logg_points:
                    for metallicity in metallicity_points:
                        pbar.set_description(
                            "Processing Teff={} K, logg={:0.2f}, Z={:+0.1f}".format(
                                teff, logg, metallicity
                            )
                        )
                        grid_point = (teff, logg, metallicity)
                        spec = PHOENIXSpectrum(
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
            flux_out = np.array(fluxes) * fluxes[0].unit
            wave_out = np.array(wavelengths) * wavelengths[0].unit

            # Make a quick-access dictionary
            n_spectra = len(grid_points)
            lookup_dict = {grid_points[i]: i for i in range(n_spectra)}
            meta = {
                "teff_points": teff_points,
                "logg_points": logg_points,
                "metallicity_points": metallicity_points,
                "grid_labels": ("T_eff", "log(g)", "Z"),
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
        wcs = None if self.wcs is None else self.wcs[key]
        mask = None if self.mask is None else self.mask[key]
        if self.meta is None:
            meta = None
        else:
            try:
                meta = self.meta[key]
            except KeyError:
                meta = self.meta

        return PHOENIXSpectrum(
            flux=flux,
            spectral_axis=spectral_axis,
            uncertainty=uncertainty,
            wcs=wcs,
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
    def metallicity_points(self):
        """What are the metallicity points of the grid?"""
        return self.meta["metallicity_points"]

    @property
    def logg_points(self):
        """What are the logg points of the grid?"""
        return self.meta["logg_points"]

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
        """Get the spectrum index associated with a given grid point"""
        return self.lookup_dict[grid_point]

    def find_nearest_teff(self, value):
        idx = (np.abs(self.teff_points - value)).argmin()
        return self.teff_points[idx]

    def find_nearest_metallicity(self, value):
        idx = (np.abs(self.metallicity_points - value)).argmin()
        return self.metallicity_points[idx]

    def show_dashboard(self, data=None, notebook_url="localhost:8888"):
        """Show an interactive dashboard for interacting with the PHOENIX grid
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
                    wavelength=self[0].wavelength.value,
                    flux=self[0].flux.value / scalar_norm,
                    native_flux=self[0].flux.value / scalar_norm,
                    native_wavelength=self[0].wavelength.value,
                )
            )

            fig = figure(
                title="PHOENIX Interactive Dashboard",
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

            fig.step(
                "wavelength",
                "flux",
                line_width=1,
                color="gray",
                source=spec_source,
                nonselection_line_color="gray",
                nonselection_line_alpha=1.0,
            )
            wl_lo, wl_hi = (
                self[0].wavelength.value.min(),
                self[0].wavelength.value.max(),
            )

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
                value=self.teff_points[1],
                step=100,
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
                step=0.50,
                title="Surface Gravity: log(g) [cm/s^2]",
                width=490,
            )

            metallicity_slider = Slider(
                start=min(self.metallicity_points),
                end=max(self.metallicity_points),
                value=0.0,
                step=0.50,
                title="Metallicity: Z",
                width=490,
            )

            r_button = Button(label=">", button_type="default", width=30)
            l_button = Button(label="<", button_type="default", width=30)

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
                    teff_message.text = "Closest grid point: {}".format(teff)
                    logg = logg_slider.value
                    metallicity = metallicity_slider.value
                    grid_point = (teff, logg, metallicity)
                    index = self.get_index(grid_point)

                    native_spec = self[index].normalize(percentile=95)
                    new_spec = native_spec.rotationally_broaden(
                        smoothing_slider.value
                    ).rv_shift(vz_slider.value)

                    spec_source.data = {
                        "native_wavelength": native_spec.wavelength.value,
                        "native_flux": native_spec.flux.value,
                        "wavelength": new_spec.wavelength.value,
                        "flux": new_spec.flux.value,
                    }

                else:
                    pass

            def update_upon_metallicity_selection(attr, old, new):
                """Callback to take action when teff slider changes"""
                metallicity = self.find_nearest_metallicity(new)
                if metallicity != old:
                    teff = self.find_nearest_teff(teff_slider.value)
                    logg = logg_slider.value
                    grid_point = (teff, logg, metallicity)
                    index = self.get_index(grid_point)

                    native_spec = self[index].normalize(percentile=95)
                    new_spec = native_spec.rotationally_broaden(
                        smoothing_slider.value
                    ).rv_shift(vz_slider.value)

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
                metallicity = self.find_nearest_metallicity(metallicity_slider.value)
                grid_point = (teff, new, metallicity)
                index = self.get_index(grid_point)

                native_spec = self[index].normalize(percentile=95)
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
                [sp4, metallicity_slider],
                [sp4, smoothing_slider],
                [sp4, vz_slider],
            )
            doc.add_root(widgets_and_figures)

        output_notebook(verbose=False, hide_banner=True)
        return show(create_interact_ui, notebook_url=notebook_url)
