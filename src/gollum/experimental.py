r"""
Experimental PHOENIX Spectrum
-----------------

A container for an experimental Phoenix grid-point spectrum of wavelength and flux :math:`F(\lambda)`.

ExpPHOENIXSpectrum
###############
"""

import numpy as np

from warnings import filterwarnings
from logging import getLogger
from gollum.phoenix import PHOENIXSpectrum, PHOENIXGrid
from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyWarning
from astropy import units as u
from specutils import Spectrum1D
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


class ExpPHOENIXGrid(PHOENIXGrid):
    """A container for an experimental grid of PHOENIX precomputed synthetic spectra of stars."""
    def show_dashboard(
        self, data=None, notebook_url="localhost:8888"
    ):  # pragma: no cover
        """Show an interactive dashboard for the experimental PHOENIX grid;
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
            spec_source = ColumnDataSource(
                data={
                    "wavelength": self[0].wavelength.value,
                    "flux": self[0].flux.value / np.percentile(self[0].flux.value, 95),
                    "native_flux": self[0].flux.value,
                    "native_wavelength": self[0].wavelength.value,
                }
            )
            wl_lo, wl_hi = (
                self[0].wavelength.value[0],
                self[0].wavelength.value[-1],
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
                    data.wavelength.value[0],
                    data.wavelength.value[-1],
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

            fig.title.offset = 280
            fig.title.text_font_size = "16pt"
            fig.yaxis.axis_label = "Normalized Flux"
            fig.xaxis.axis_label = "Wavelength (\u212B)"
            fig.axis.axis_label_text_font_style = "bold"
            fig.x_range = Range1d(start=wl_lo, end=wl_hi)
            fig.y_range = Range1d(start=0, end=1.5)
            fig.legend.location = "top_right"
            fig.step(
                "wavelength",
                "flux",
                line_width=1,
                color="crimson",
                source=spec_source,
                nonselection_line_color="crimson",
                nonselection_line_alpha=1.0,
                legend_label="PHOENIX Model",
            )

            smoothing_slider = Slider(
                start=0.1,
                end=200,
                value=0.1,
                step=0.1,
                title="Rotational Broadening [km/s]",
                width=460,
                format='0.f',
                bar_color="blue",
            )
            rv_slider = Slider(
                start=-200,
                end=200,
                value=0.00,
                step=0.05,
                title="Radial Velocity [km/s]",
                width=460,
                bar_color="blue",
            )
            teff_slider = Slider(
                start=self.teff_points[0],
                end=self.teff_points[-1],
                value=self.teff_points[0],
                step=100,
                title="Effective Temperature [K]",
                width=460,
                bar_color="red",
                margin=(0, 20, 0, 0),
            )
            logg_slider = Slider(
                start=self.logg_points[0],
                end=self.logg_points[-1],
                value=self.logg_points[0],
                step=0.50,
                title="Surface Gravity [cm/s\u00b2]",
                width=460,
                bar_color="red",
                margin=(0, 20, 0, 0),
            )
            metallicity_slider = Slider(
                start=self.metallicity_points[0],
                end=self.metallicity_points[-1],
                value=self.metallicity_points[0],
                step=0.50,
                title="Metallicity [dex]",
                width=460,
                bar_color="red",
                margin=(0, 20, 0, 0),
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
            spot_temp_slider = Slider(
                start=2300,
                end=4000,
                value=2300,
                step=100,
                title="Starspot Temperature [K]",
                width=460,
                bar_color="maroon",
                margin=(0, 20, 0, 0),
            )
            fill_factor_slider = Slider(
                start=0,
                end=1,
                value=0,
                step=0.05,
                title="Starspot Filling Factor",
                width=460,
                bar_color="maroon",
            )
            continuum_toggle = Toggle(
                label="Fit Continuum (disables normalization)", button_type="success"
            )

            def update_continuum(active):
                """Callback that toggles continuum auto-fit"""
                if active:
                    spec_source.data["flux"] = (
                        PHOENIXSpectrum(
                            spectral_axis=spec_source.data["wavelength"] * u.AA,
                            flux=spec_source.data["flux"] * u.dimensionless_unscaled,
                        )
                        .tilt_to_data(data)
                        .flux.value
                    )
                    scale_slider.disabled = True
                    continuum_toggle.label = "Undo Continuum (enables normalization)"
                else:
                    spec_source.data["flux"] = (
                        PHOENIXSpectrum(
                            spectral_axis=spec_source.data["wavelength"] * u.AA,
                            flux=spec_source.data["native_flux"]
                            * u.dimensionless_unscaled,
                        )
                        .normalize(percentile=95)
                        .rotationally_broaden(smoothing_slider.value)
                        .flux.value
                        * scale_slider.value
                    )
                    scale_slider.disabled = False
                    continuum_toggle.label = "Fit Continuum (disables normalization)"

            def update_rv(attr, old, new):
                """Callback that RV shifts the spectrum"""
                spec_source.data["wavelength"] = (
                    PHOENIXSpectrum(
                        spectral_axis=spec_source.data["native_wavelength"] * u.AA,
                        flux=spec_source.data["flux"] * u.dimensionless_unscaled,
                    )
                    .rv_shift(new)
                    .wavelength.value
                )

            def update_smoothing(attr, old, new):
                """Callback that rotationally broadens the spectrum"""
                spec = (
                    PHOENIXSpectrum(
                        spectral_axis=spec_source.data["wavelength"] * u.AA,
                        flux=spec_source.data["native_flux"] * u.dimensionless_unscaled,
                    )
                    .normalize(percentile=95)
                    .rotationally_broaden(new)
                )
                spec_source.data["flux"] = (
                    spec.tilt_to_data(data).flux.value
                    if continuum_toggle.active
                    else spec.flux.value * scale_slider.value
                )

            def update_scale(attr, old, new):
                """Callback that scales the spectrum"""
                spec_source.data["flux"] *= new / old

            def update_native(attr, old, new):
                """Callback that updates the intrinsic parameters behind the spectrum"""
                teff_slider.value = self.find_nearest_teff(teff_slider.value)
                spot_temp_slider.value = self.find_nearest_teff(spot_temp_slider.value)
                metallicity_slider.value = self.find_nearest_metallicity(
                    metallicity_slider.value
                )
                spot = PHOENIXSpectrum(
                    teff=spot_temp_slider.value,
                    logg=logg_slider.value,
                    metallicity=metallicity_slider.value,
                    wl_lo=spec_source.data["native_wavelength"][0],
                    wl_hi=spec_source.data["native_wavelength"][-1],
                )
                base = PHOENIXSpectrum(
                    teff=teff_slider.value,
                    logg=logg_slider.value,
                    metallicity=metallicity_slider.value,
                    wl_lo=spec_source.data["native_wavelength"][0],
                    wl_hi=spec_source.data["native_wavelength"][-1],
                )
                spec_source.data["native_flux"] = (
                    base.flux.value * (1 - fill_factor_slider.value)
                    + spot.flux.value * fill_factor_slider.value
                )
                final = (
                    PHOENIXSpectrum(
                        spectral_axis=spec_source.data["wavelength"] * u.AA,
                        flux=spec_source.data["native_flux"] * u.dimensionless_unscaled,
                    )
                    .normalize(percentile=95)
                    .rotationally_broaden(smoothing_slider.value)
                )
                spec_source.data["flux"] = (
                    final.tilt_to_data(data).flux.value
                    if continuum_toggle.active
                    else final.flux.value * scale_slider.value
                )

            continuum_toggle.on_click(update_continuum)
            rv_slider.on_change("value", update_rv)
            smoothing_slider.on_change("value", update_smoothing)
            scale_slider.on_change("value", update_scale)
            teff_slider.on_change("value", update_native)
            logg_slider.on_change("value", update_native)
            metallicity_slider.on_change("value", update_native)
            spot_temp_slider.on_change("value", update_native)
            fill_factor_slider.on_change("value", update_native)

            sp = Spacer(width=20)
            doc.add_root(
                layout(
                    [sp, fig],
                    [sp, continuum_toggle],
                    [sp, teff_slider, smoothing_slider],
                    [sp, logg_slider, rv_slider],
                    [sp, metallicity_slider, scale_slider],
                    [sp, spot_temp_slider, fill_factor_slider],
                    background="whitesmoke",
                )
            )

        output_notebook(verbose=False, hide_banner=True)
        return show(create_interact_ui, notebook_url=notebook_url)
