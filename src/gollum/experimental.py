r"""
PHOENIX Spectrum
-----------------

A container for a single Phoenix grid-point spectrum of wavelength and flux :math:`F(\lambda)`.

PHOENIXSpectrum
###############
"""

import numpy as np

from warnings import filterwarnings
from logging import getLogger
from gollum.phoenix import PHOENIXGrid, PHOENIXSpectrum
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy import units as u
from specutils import Spectrum1D
from bokeh.io import show, output_notebook
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import Slider, Range1d
from bokeh.layouts import layout, Spacer

log = getLogger(__name__)

#  See Issue: https://github.com/astropy/specutils/issues/779
filterwarnings("ignore", category=AstropyDeprecationWarning)
# See Issue: https://github.com/astropy/specutils/issues/800
filterwarnings("ignore", category=RuntimeWarning)


class ExperimentalPHOENIXGrid(PHOENIXGrid):
    """
    An experimental version of PHOENIXGrid.
    """

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
            scalar_norm = np.percentile(self[0].flux.value, 95)
            spec_source = ColumnDataSource(
                data={
                    "wavelength": self[0].wavelength.value,
                    "flux": self[0].flux.value / scalar_norm,
                    "native_flux": self[0].flux.value / scalar_norm,
                    "native_wavelength": self[0].wavelength.value,
                }
            )
            wl_lo, wl_hi = (
                self[0].wavelength.value.min(),
                self[0].wavelength.value.max(),
            )
            fig = figure(
                title="PHOENIX Interactive Dashboard",
                plot_height=500,
                plot_width=950,
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
                assert (new_lo < wl_hi) & (
                    new_hi > wl_lo
                ), "Data should overlap the models, double check your wavelength limits."
                wl_lo, wl_hi = new_lo, new_hi

                fig.step(
                    "wavelength",
                    "flux",
                    line_width=1,
                    color="blue",
                    source=ColumnDataSource(
                        data={
                            "wavelength": data.wavelength.value,
                            "flux": data.flux.value,
                        }
                    ),
                )

            fig.title.offset = 350
            fig.yaxis.axis_label = "Flux"
            fig.xaxis.axis_label = "Wavelength (\u00B5m)"
            fig.x_range = Range1d(start=wl_lo, end=wl_hi)
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

            smoothing_slider = Slider(
                start=0.1,
                end=200,
                value=0.1,
                step=0.1,
                title="Rotational Broadening: v sin(i) [km/s]",
                width=460,
            )
            vz_slider = Slider(
                start=-200,
                end=200,
                value=0.00,
                step=0.05,
                title="Radial Velocity: RV [km/s]",
                width=460,
                format="0.000f",
            )
            teff_slider = Slider(
                start=min(self.teff_points),
                end=max(self.teff_points),
                value=min(self.teff_points),
                step=100,
                title="Effective Temperature: T_eff [K]",
                width=460,
            )
            logg_slider = Slider(
                start=min(self.logg_points),
                end=max(self.logg_points),
                value=min(self.logg_points),
                step=0.50,
                title="Surface Gravity: log(g) [cm/s^2]",
                width=460,
            )
            metallicity_slider = Slider(
                start=min(self.metallicity_points),
                end=max(self.metallicity_points),
                value=min(self.metallicity_points),
                step=0.50,
                title="Metallicity: Z",
                width=460,
            )
            scale_slider = Slider(
                start=0.1,
                end=2.0,
                value=1.0,
                step=0.005,
                title="Normalization Scalar",
                width=460,
            )
            spot_temp_slider = Slider(
                start=2000,
                end=4000,
                value=2000,
                step=100,
                title="Starspot Temperature [K]",
                width=460,
            )
            spot_area_slider = Slider(
                start=0,
                end=1,
                value=0,
                step=0.05,
                title="Starspot Filling Factor",
                width=460,
            )

            def update_upon_scale(attr, old, new):
                """Callback to take action when normalization slider changes"""
                new_spec = (
                    PHOENIXSpectrum(
                        spectral_axis=spec_source.data["native_wavelength"]
                        * u.Angstrom,
                        flux=spec_source.data["native_flux"] *
                        u.dimensionless_unscaled,
                    )
                    .rotationally_broaden(smoothing_slider.value)
                    .multiply(new * u.dimensionless_unscaled)
                    .rv_shift(vz_slider.value)
                )
                spec_source.data["flux"] = new_spec.flux.value

            def update_upon_smooth(attr, old, new):
                """Callback to take action when smoothing slider changes"""
                new_spec = (
                    PHOENIXSpectrum(
                        spectral_axis=spec_source.data["native_wavelength"]
                        * u.Angstrom,
                        flux=spec_source.data["native_flux"] *
                        u.dimensionless_unscaled,
                    )
                    .rotationally_broaden(new)
                    .multiply(scale_slider.value * u.dimensionless_unscaled)
                )
                spec_source.data["flux"] = new_spec.flux.value

            def update_upon_vz(attr, old, new):
                """Callback to take action when vz slider changes"""
                new_spec = PHOENIXSpectrum(
                    spectral_axis=spec_source.data["native_wavelength"] * u.Angstrom,
                    flux=spec_source.data["native_flux"] *
                    u.dimensionless_unscaled,
                ).rv_shift(new)
                spec_source.data["wavelength"] = new_spec.wavelength.value

            def update_upon_teff_selection(attr, old, new):
                """Callback to take action when teff slider changes"""
                teff = self.find_nearest_teff(new)
                if teff != old:
                    point = (teff, logg_slider.value, metallicity_slider.value)
                    native_spec = self[self.get_index(point)]
                    # ambient absolute scalar value
                    scalar_value_amb = np.percentile(native_spec.flux, 95)
                    native_spec = native_spec.divide(scalar_value_amb)

                    # Todo: also read in the starspot spectrum.
                    # ...
                    # ...
                    # ...
                    # ...
                    # scalar_value_spot
                    # absolute_ratio = scalar_value_spot / scalar_value_amb
                    # ...
                    # ...

                    # Multiply by area ratio:

                    # spot_spec = f_spot * spot_spec_abs...

                    new_spec = (
                        native_spec.rotationally_broaden(
                            smoothing_slider.value)
                        .rv_shift(vz_slider.value)
                        .multiply(scale_slider.value * u.dimensionless_unscaled)
                    )

                    spec_source.data = {
                        "native_wavelength": native_spec.wavelength.value,
                        "native_flux": native_spec.flux.value,
                        "wavelength": new_spec.wavelength.value,
                        "flux": new_spec.flux.value,
                    }
                teff_slider.value = teff

            def update_upon_metallicity_selection(attr, old, new):
                """Callback to take action when metallicity slider changes"""
                metallicity = self.find_nearest_metallicity(new)
                if metallicity != old:
                    teff = self.find_nearest_teff(teff_slider.value)
                    point = (teff, logg_slider.value, metallicity)
                    native_spec = self[self.get_index(
                        point)].normalize(percentile=95)
                    new_spec = (
                        native_spec.rotationally_broaden(
                            smoothing_slider.value)
                        .rv_shift(vz_slider.value)
                        .multiply(scale_slider.value * u.dimensionless_unscaled)
                    )

                    spec_source.data = {
                        "native_wavelength": native_spec.wavelength.value,
                        "native_flux": native_spec.flux.value,
                        "wavelength": new_spec.wavelength.value,
                        "flux": new_spec.flux.value,
                    }

            def update_upon_logg_selection(attr, old, new):
                """Callback to take action when logg slider changes"""
                teff = self.find_nearest_teff(teff_slider.value)
                metallicity = self.find_nearest_metallicity(
                    metallicity_slider.value)
                point = (teff, new, metallicity)
                native_spec = self[self.get_index(
                    point)].normalize(percentile=95)
                new_spec = (
                    native_spec.rotationally_broaden(smoothing_slider.value)
                    .rv_shift(vz_slider.value)
                    .multiply(scale_slider.value * u.dimensionless_unscaled)
                )

                spec_source.data = {
                    "native_wavelength": native_spec.wavelength.value,
                    "native_flux": native_spec.flux.value,
                    "wavelength": new_spec.wavelength.value,
                    "flux": new_spec.flux.value,
                }

            def update_spot_temp(attr, old, new):
                """Callback to take action when spot temp slider changes"""
                pass

            def update_spot_area(attr, old, new):
                """"Callback to take action when spot area changes"""
                pass

            smoothing_slider.on_change("value", update_upon_smooth)
            vz_slider.on_change("value", update_upon_vz)
            teff_slider.on_change("value", update_upon_teff_selection)
            logg_slider.on_change("value", update_upon_logg_selection)
            metallicity_slider.on_change(
                "value", update_upon_metallicity_selection)
            scale_slider.on_change("value", update_upon_scale)

            sp = Spacer(width=20)
            doc.add_root(
                layout(
                    [fig],
                    [teff_slider, sp, smoothing_slider],
                    [logg_slider, sp, vz_slider],
                    [metallicity_slider, sp, scale_slider],
                    [spot_temp_slider, sp, spot_area_slider],
                )
            )

        output_notebook(verbose=False, hide_banner=True)
        return show(create_interact_ui, notebook_url=notebook_url)
