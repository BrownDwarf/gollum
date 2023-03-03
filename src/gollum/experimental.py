r"""
Experimental PHOENIX Spectrum
-----------------

A container for an experimental Phoenix grid-point spectrum of wavelength and flux :math:`F(\lambda)`.

ExpPHOENIXSpectrum
###############
"""

from gollum.phoenix import *

log = getLogger(__name__)

filterwarnings("ignore", category=AstropyDeprecationWarning)
filterwarnings("ignore", category=AstropyWarning)
filterwarnings("ignore", category=RuntimeWarning)


class ExpPHOENIXGrid(PHOENIXGrid):
    """A container for an experimental grid of PHOENIX precomputed synthetic spectra of stars."""

    def show_dashboard(self, data=None, url="localhost:8888"):  # pragma: no cover
        """Show an interactive dashboard for the experimental PHOENIX grid;
        heavily inspired by the lightkurve .interact() method

        Parameters
        ----------
        data: Spectrum1D-like
            A normalized data spectrum over which to plot the models
        url: str
            Location of the Jupyter notebook page (default: "localhost:8888")
            If you are running on a different location, you
            will need to supply this value for the application to display
            properly. If no protocol is supplied in the URL, e.g. if it is
            of the form "localhost:8888", then "http" will be used.
        """

        def create_interact_ui(doc):
            self.cur_photo = self[0]
            self.cur_spot = self[0]
            self.cur_total = self[0]
            wl_i, flux_i = self[0].wavelength.value, self[0].normalize(95).flux.value
            cds = ColumnDataSource(
                data={
                    "wl": wl_i,
                    "nat_wl": wl_i,
                    "flux": flux_i,
                    "nat_flux": self[0].flux.value,
                    "photo_flux": flux_i,
                    "photo_nat": flux_i,
                    "spot_flux": flux_i * 0,
                    "spot_nat": flux_i,
                }
            )
            wl_lo, wl_hi = wl_i[0], wl_i[-1]
            fig = figure(
                title="PHOENIX Interactive Dashboard",
                width=950,
                height=500,
                tools="pan,wheel_zoom,box_zoom,tap,reset",
                toolbar_location="below",
                border_fill_color="whitesmoke",
            )

            if data:
                new_lo, new_hi = data.wavelength.value[0], data.wavelength.value[-1]
                assert (
                    wl_lo < new_lo < new_hi < wl_hi
                ), "Data must overlap models, expand your wavelength range."
                wl_lo, wl_hi = new_lo, new_hi

                fig.step(
                    x="wl",
                    y="flux",
                    color="black",
                    legend_label=data.meta.get("header").get("OBJECT"),
                    source=ColumnDataSource(
                        data={"wl": data.wavelength.value, "flux": data.flux.value,}
                    ),
                )
            fig.title.align = "center"
            fig.title.text_font_size = "18pt"
            fig.yaxis.axis_label = "Normalized Flux"
            fig.xaxis.axis_label = "Wavelength (\u212B)"
            fig.axis.axis_label_text_font_style = "bold"
            fig.x_range = Range1d(start=wl_lo, end=wl_hi)
            fig.y_range = Range1d(start=0, end=1.5)
            fig.legend.location = "top_right"
            fig.legend.click_policy = "hide"
            fig.step(
                x="wl",
                y="flux",
                color="crimson",
                source=cds,
                legend_label="PHOENIX Model: Total Flux",
            )
            fig.step(
                x="wl",
                y="photo_flux",
                color="violet",
                source=cds,
                legend_label="PHOENIX Model: Photosphere Flux",
                level="underlay",
            ).visible = False
            fig.step(
                x="wl",
                y="spot_flux",
                color="lavender",
                source=cds,
                legend_label="PHOENIX Model: Starspot Flux",
                level="underlay",
            ).visible = False

            smooths = Slider(
                start=0.1,
                end=200,
                value=0.1,
                step=0.1,
                title="Rotational Broadening [km/s]",
                width=460,
                format="0.f",
                bar_color="blue",
            )
            rvs = Slider(
                start=-200,
                end=200,
                value=0.00,
                title="Radial Velocity [km/s]",
                width=460,
                bar_color="blue",
            )
            teffs = Slider(
                start=self.teff_points[0],
                end=self.teff_points[-1],
                value=self.teff_points[0],
                step=100,
                title="Effective Temperature [K]",
                width=460,
                bar_color="red",
                margin=(0, 20, 0, 0),
            )
            loggs = Slider(
                start=self.logg_points[0],
                end=self.logg_points[-1],
                value=self.logg_points[0],
                step=0.50,
                title="Surface Gravity [cm/s\u00b2]",
                width=460,
                bar_color="red",
                margin=(0, 20, 0, 0),
            )
            Zs = Slider(
                start=self.Z_points[0],
                end=self.Z_points[-1],
                value=self.Z_points[0],
                step=0.50,
                title="Metallicity [dex]",
                width=460,
                bar_color="red",
                margin=(0, 20, 0, 0),
            )
            scales = Slider(
                start=0.1,
                end=2.0,
                value=1.0,
                step=0.005,
                title="Scale Factor",
                width=460,
                bar_color="black",
            )
            spot_temps = Slider(
                start=self.teff_points[0],
                end=self.teff_points[-1],
                value=self.teff_points[0],
                step=100,
                title="Starspot Temperature [K]",
                width=460,
                bar_color="maroon",
                margin=(0, 20, 0, 0),
            )
            fills = Slider(
                start=0,
                end=1,
                value=0,
                step=0.05,
                title="Starspot Filling Factor",
                width=460,
                bar_color="maroon",
            )
            continuum = Toggle(
                label="Fit Continuum (disables scaling)", button_type="success",
            )

            def toggle_continuum(active):
                """Callback that toggles continuum auto-fit"""
                if active:
                    cds.data["flux"] = self.cur_total.tilt_to_data(data).flux.value
                    scales.disabled = True
                    continuum.label = "Undo Continuum (enables scaling)"
                else:
                    cds.data["flux"] = self.cur_total.flux.value * scales.value
                    scales.disabled = False
                    continuum.label = "Fit Continuum (disables scaling)"

            def update_rv(attr, old, new):
                """Callback that RV shifts the spectrum"""
                cds.data["wl"] = self[0].rv_shift(new).wavelength.value

            def update_smoothing(attr, old, new):
                """Callback that rotationally broadens the spectrum"""
                photo = self.cur_photo.rotationally_broaden(new)
                spot = self.cur_spot.rotationally_broaden(new)
                self.cur_total = (
                    photo * (1 - fills.value) + spot * fills.value
                ).normalize(95)

                cds.data["photo_flux"] = photo.normalize(95).flux.value * (
                    1 - fills.value
                )
                cds.data["spot_flux"] = spot.normalize(95).flux.value * fills.value
                cds.data["flux"] = (
                    self.cur_total.flux.value * scales.value
                    if not continuum.active
                    else self.cur_total.tilt_to_data(data).flux.value
                )

            def update_scale(attr, old, new):
                """Callback that scales the spectra"""
                cds.data["flux"] *= new / old

            def update_spot(attr, old, new):
                """Callback that updates the starspot's temperature"""
                spot_temps.value = self.find_nearest_teff(new)
                self.cur_spot = self[self.get_index((new, loggs.value, Zs.value))]
                spot = self.cur_spot.rotationally_broaden(smooths.value)
                self.cur_total = (
                    self.cur_photo.rotationally_broaden(smooths.value)
                    * (1 - fills.value)
                    + spot * fills.value
                ).normalize(95)

                cds.data["spot_flux"] = spot.normalize(95).flux.value * fills.value
                cds.data["flux"] = (
                    self.cur_total.flux.value * scales.value
                    if not continuum.active
                    else self.cur_total.tilt_to_data(data).flux.value
                )

            def update_fill_factor(attr, old, new):
                if old not in (0, 1):
                    cds.data["spot_flux"] *= new / old
                    cds.data["photo_flux"] *= (1 - new) / (1 - old)
                    self.cur_total = self[0]._copy(
                        flux=(cds.data["photo_flux"] + cds.data["spot_flux"]) * DV
                    )
                    cds.data["flux"] = (
                        self.cur_total.flux.value * scales.value
                        if not continuum.active
                        else self.cur_total.tilt_to_data(data).flux.value
                    )
                else:
                    update_smoothing("value", 0, smooths.value)

            def update_native(attr, old, new):
                """Callback that updates the intrinsic parameters behind the spectrum"""
                teffs.value = self.find_nearest_teff(teffs.value)
                spot_temps.value = self.find_nearest_teff(spot_temps.value)
                Zs.value = self.find_nearest_Z(Zs.value)
                cds.data["photo_flux"] = PHOENIXSpectrum(
                    teff=teffs.value,
                    logg=loggs.value,
                    Z=Zs.value,
                    wl_lo=cds.data["nat_wl"][0],
                    wl_hi=cds.data["nat_wl"][-1],
                ).normalize(95).rotationally_broaden(smooths.value).flux.value * (
                    1 - fills.value
                )

                cds.data["spot_flux"] = (
                    PHOENIXSpectrum(
                        teff=spot_temps.value,
                        logg=loggs.value,
                        Z=Zs.value,
                        wl_lo=cds.data["nat_wl"][0],
                        wl_hi=cds.data["nat_wl"][-1],
                    )
                    .normalize(95)
                    .rv_shift(rvs.value)
                    .rotationally_broaden(smooths.value)
                    .flux.value
                    * fills.value
                )

                cds.data["nat_flux"] = cds.data["photo_flux"] + cds.data["spot_flux"]
                final = PrecomputedSpectrum(
                    spectral_axis=cds.data["wl"] * u.AA, flux=cds.data["nat_flux"] * DV,
                )
                cds.data["flux"] = (
                    final.tilt_to_data(data).flux.value
                    if continuum.active
                    else final.flux.value * scales.value
                )

            continuum.on_click(toggle_continuum)
            rvs.on_change("value", update_rv)
            smooths.on_change("value", update_smoothing)
            scales.on_change("value", update_scale)
            teffs.on_change("value", update_native)
            loggs.on_change("value", update_native)
            Zs.on_change("value", update_native)
            spot_temps.on_change("value", update_spot)
            fills.on_change("value", update_fill_factor)

            sp = Spacer(width=20)
            doc.add_root(
                layout(
                    [sp, fig],
                    [sp, continuum],
                    [sp, teffs, smooths],
                    [sp, loggs, rvs],
                    [sp, Zs, scales],
                    [sp, spot_temps, fills],
                    background="whitesmoke",
                )
            )

        output_notebook(verbose=False, hide_banner=True)
        return show(create_interact_ui, notebook_url=url)
