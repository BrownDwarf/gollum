r"""
Precomputed Spectrum
-----------------

A container for a single precomputed synthetic model spectrum at a single grid-point, with wavelength and flux :math:`F(\lambda)`.

PrecomputedSpectrum
###############
"""
import specutils
import numpy as np
import matplotlib.pyplot as plt

from copy import deepcopy
from logging import getLogger
from warnings import filterwarnings
from gollum.utilities import apply_numpy_mask
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from specutils import Spectrum1D
from specutils.manipulation import LinearInterpolatedResampler
from specutils.fitting import fit_generic_continuum
from astropy import units as u, constants as const
from astropy.modeling.physical_models import BlackBody
from astropy.utils.exceptions import AstropyDeprecationWarning

log = getLogger(__name__)

#  See Issue: https://github.com/astropy/specutils/issues/779
filterwarnings("ignore", category=AstropyDeprecationWarning)
# See Issue: https://github.com/astropy/specutils/issues/800
filterwarnings("ignore", category=RuntimeWarning)


class PrecomputedSpectrum(Spectrum1D):
    """
    An abstract container base class for a Precomputed spectrum
    """

    def __init__(self, *args, **kwargs):
        # Todo, we could put the wavelength limits in here.
        super().__init__(*args, **kwargs)

    @property
    def velocity_spacing(self):
        """The velocity sampling of the spectrum
        
        Returns
        -------
        velocity_spacing : np.array
            vector of per pixel velocity sampling
        """
        c_kmps = const.c.to(u.km / u.s).value
        per_pixel_velocity_sampling = (
            c_kmps * np.diff(self.wavelength.value) / self.wavelength.value[1:]
        )
        velocity_variation = np.std(per_pixel_velocity_sampling)

        # Problems may arise if pixel spacings jump around by more than 50 m/s:
        velocity_variation_threshold = 0.05  # km/s
        if velocity_variation > velocity_variation_threshold:
            log.info(
                f"""Your velocity sampling appears to be non-uniform at the {velocity_variation:0.4f} km/s level, 
                which could affect future convolution processes. Consider applying the `resample_to_uniform_in_velocity` method."""
            )
        return np.median(per_pixel_velocity_sampling) * u.km / u.s

    apply_boolean_mask = apply_numpy_mask

    def normalize(self, percentile=50):
        """Normalize spectrum by some given percentile

        Parameters
        ----------
        percentile : int
            The percentile to which the spectrum will be normalized (default: 50th percentile)

        Returns
        -------
        normalized_spec : PrecomputedSpectrum
            Normalized spectrum
        """
        spec = self._copy(
            spectral_axis=self.wavelength.value * self.wavelength.unit, wcs=None
        )
        scalar_flux = (
            np.nanpercentile(spec.flux.value, percentile) * spec.flux.unit
            if percentile
            else np.nanmedian(spec.flux.value) * spec.flux.unit
        )

        return spec.divide(scalar_flux, handle_meta="first_found")

    def rotationally_broaden(self, vsini, u1=0.0, u2=0.0):
        r"""Rotationally broaden the spectrum for a given :math:`v\sin{i}`
        Implementation inspired by https://github.com/HajimeKawahara/exojax

        Known limitation: If the wavelength sampling changes with wavelength,
          the convolution becomes inaccurate.  It may be better to FFT,
          following Starfish.

        Parameters
        ----------
        vsini : int
            :math:`v\sin{i}` in units of km/s
        u1 : float
            Limb-darkening coefficient 1
        u2 : float
            Limb-darkening coefficient 2

        Returns
        -------
        broadened_spec : PrecomputedSpectrum
            Rotationally Broadened Spectrum
        """
        lam0 = np.median(self.wavelength.value)
        x2 = (299792.458 * (self.wavelength.value - lam0) / (lam0 * vsini)) ** 2
        kernel = np.where(
            x2 < 1,
            np.pi / 2 * u1 * (1 - x2)
            + np.sqrt(1 - x2) * (2 - 2 * u1 - 4 / 3 * u2 * u2 * x2),
            0,
        )
        kernel, positive_elements = kernel / np.sum(kernel, axis=0), kernel > 0
        return (
            self._copy(
                flux=np.convolve(
                    self.flux.value, kernel[positive_elements], mode="same"
                )
                * self.flux.unit
            )
            if positive_elements.any()
            else self
        )

    def instrumental_broaden(self, resolving_power=55000):
        r"""Instrumentally broaden the spectrum for a given instrumental resolution

        Known limitation: If the wavelength sampling changes with wavelength,
          the convolution becomes inaccurate.  It may be better to FFT,
          following Starfish.

        Parameters
        ----------
        resolving_power : int
            Instrumental resolving power :math:`R = \frac{\lambda}{\delta \lambda}`

        Returns
        -------
        broadened_spec : PrecomputedSpectrum
            Instrumentally broadened spectrum
        """
        # In detail the spectral resolution is wavelength dependent...
        # For now we assume a constant resolving power
        angstroms_per_pixel = np.median(np.diff(self.wavelength.value))
        lam0 = np.median(self.wavelength.value)
        delta_lam = lam0 / resolving_power

        scale_factor = 2.355
        sigma = delta_lam / scale_factor / angstroms_per_pixel

        convolved_flux = gaussian_filter1d(self.flux.value, sigma) * self.flux.unit
        return self._copy(flux=convolved_flux)

    def rv_shift(self, rv):
        """Shift the spectrum by a radial velocity

        Parameters
        ----------
        rv : int
            Radial velocity in km/s

        Returns
        -------
        shifted_spec : PrecomputedSpectrum
            RV-Shifted Spectrum
        """
        output = deepcopy(self)
        output.radial_velocity = rv * u.km / u.s
        return self._copy(
            spectral_axis=output.wavelength.value * output.wavelength.unit, wcs=None
        )

    def resample(self, target_spectrum):
        """Resample spectrum at the wavelength points of another spectrum

        Parameters
        ----------
        target_spectrum : Spectrum1D
            Spectrum whose wavelength grid you seek to match

        Returns
        -------
        resampled_spec : PrecomputedSpectrum
            Resampled spectrum
        """
        output = LinearInterpolatedResampler()(self, target_spectrum.wavelength)

        return self._copy(
            spectral_axis=output.wavelength.value * output.wavelength.unit,
            flux=output.flux,
            wcs=None,
        )

    def resample_to_uniform_in_velocity(self, oversample=1.4):
        """Resample spectrum to a uniform-in-velocity pixel spacing

        Parameters
        ----------
        oversample : float
            The desired oversampling in velocity, compared to the typical velocity sampling at original pixel spacing.
            Typically, you will want to oversample to ensure the narrowest lines remain resolved at the new sampling,
            at the expense of more pixels than you started with.
            If your original spectrum has large gaps, you may end up with many more pixels than you started with.

        Returns
        -------
        resampled_spec : PrecomputedSpectrum
            Resampled spectrum
        """

        c_kmps = const.c.to(u.km / u.s).value
        lambda_0 = self.wavelength.value.min()
        lambda_max = self.wavelength.value.max()

        # Compute the per pixel resolving power
        # Note: assumes a reasonably contiguous sampling in wavelength
        # Major gaps in wavelength will mess up this approach.
        per_pixel_resolution = self.wavelength.value[1:] / np.diff(
            self.wavelength.value
        )
        median_pixel_resolution = np.median(per_pixel_resolution)
        max_pixel_resolution = np.max(per_pixel_resolution)
        new_pixel_resolution = median_pixel_resolution * oversample

        if new_pixel_resolution < max_pixel_resolution:
            log.info(
                f"""You are trying to oversample the spectrum by a factor of {oversample}.
                The highest existing per-pixel resolution of the spectrum was {max_pixel_resolution:0.1f}, whereas your new resolution is only {new_pixel_resolution:0.1f}.
                You may want to consider a higher oversample factor to avoid information loss."""
            )

        velocity_resolution_kmps = c_kmps / new_pixel_resolution

        velocity_max = c_kmps * np.log(lambda_max / lambda_0)
        velocity_vector = np.arange(0, velocity_max, velocity_resolution_kmps)
        new_wavelength_sampling = lambda_0 * np.exp(velocity_vector / c_kmps)

        output = LinearInterpolatedResampler()(self, new_wavelength_sampling * u.AA)

        return self._copy(
            spectral_axis=output.wavelength.value * output.wavelength.unit,
            flux=output.flux,
            wcs=None,
        )

    def decimate(self, decimation_factor=0.1, n_pixels=None, resolving_power=None):
        """Decimate the number of samples in the spectrum

        Parameters
        ----------
        decimation_factor : float
            The fraction of pixels to keep. Default: 0.1
        n_pixels  : int
            The number of pixels to keep. Default: 2,000
        resolving_power : int
            The resolving power of the new spectrum.  Default: 3,000

        Returns
        -------
        decimated_spec : PrecomputedSpectrum
            Decimated Spectrum
        """
        if n_pixels or resolving_power:
            raise NotImplementedError(
                "n_pixels and resolving_power are not implemented yet"
            )

        return self.resample_to_uniform_in_velocity(oversample=decimation_factor)

    def get_blackbody_spectrum(self, teff=None):
        """Get the blackbody spectrum associated with the input model"""
        if not teff:
            if hasattr(self, "teff") and self.teff:
                teff = self.teff
            else:
                raise NotImplementedError(
                    "Your subclass may not have implemented the teff attribute yet."
                )

        blackbody_model = BlackBody(temperature=teff * u.Kelvin)
        blackbody_flux_per_Hz = blackbody_model(self.wavelength) * np.pi * u.steradian

        flux_unit = self.flux.unit
        normalize = False
        if flux_unit == u.dimensionless_unscaled:
            if "native_flux_unit" in self.meta:
                flux_unit = self.meta["native_flux_unit"]
                normalize = True
            else:
                raise NotImplementedError(
                    "We do not yet support inferring units for this subclass"
                )

        blackbody_flux = blackbody_flux_per_Hz.to(
            flux_unit, equivalencies=u.spectral_density(self.wavelength)
        )

        output = PrecomputedSpectrum(flux=blackbody_flux, spectral_axis=self.wavelength)

        return output.normalize() if normalize else output

    def divide_by_blackbody(self, teff=None):
        """Divide the spectrum by a blackbody

        Parameters
        ----------
        teff : float
            The effective temperature of the blackbody to divide by.

        Returns
        -------
        divided_spec : PrecomputedSpectrum
            The spectrum after being divided by the blackbody.
        """
        return self.divide(self.get_blackbody_spectrum(teff=teff), handle_meta="ff")

    def tilt_to_data(self, target_spectrum, return_model=False):
        """Tilt the template towards a data spectrum by fitting and dividing by a low-order polynomial

        Parameters
        ----------
        target_spectrum : Spectrum1D
            A Spectrum1D object whose flux tilt to match against
        return_model : bool
            Whether or not to return the model

        Returns
        -------
        tilted_spec : PrecomputedSpectrum
            Tilted spectrum
        """
        resampled_self = self.resample(target_spectrum)
        model = fit_generic_continuum(
            target_spectrum.divide(resampled_self, handle_meta="ff")
        )
        tilted_spec = self.multiply(model(self.wavelength))

        return (tilted_spec, model) if return_model else tilted_spec

    def fit_continuum(self, pixel_distance=5001, polyorder=3, return_coeffs=False):
        """Finds the low frequency continuum trend using scipy's find_peaks filter and linear algebra.
        Currently broken

        Parameters
        ----------
        pixel_distance : int
            The minimum separation between peaks, in pixels. Default = 5001 px
        polyorder : int
            The polynomial degree to be used for peak fitting
        return_coeffs : bool
            If `True`, returns a 2-tuple (trend_spec, trend_coeffs) where trend_spec is the fitted trend.

        Returns
        -------
        spec_out : PrecomputedSpectrum
            New Spectrum object representing a fit of the continuum.
        If ``return_coeffs`` is set to ``True``, this method will also return:
            coeffs : np.array
                New vector of polynomial coefficients that reproduce the trend.
        """

        x_vector, y_vector = self.wavelength.value, self.flux.value

        if pixel_distance > len(x_vector):
            raise ValueError(
                "Please provide a pixel_distance smaller than the spectrum length."
            )

        peak_inds = find_peaks(y_vector, distance=pixel_distance)[0]

        x_peaks, y_peaks = x_vector[peak_inds], y_vector[peak_inds]

        A_matrix, A_full = np.vander(x_peaks, polyorder), np.vander(x_vector, polyorder)

        coeffs = np.linalg.lstsq(A_matrix, y_peaks, rcond=None)

        spec_out = self._copy(flux=np.dot(coeffs, A_full.T) * self.flux.unit)

        return (spec_out, coeffs) if return_coeffs else spec_out

    def to_pandas(self):
        """Export the spectrum to a pandas dataframe"""
        try:
            from pandas import DataFrame
        except ImportError:
            log.error("The to_pandas method requires the optional dependency pandas")

        return DataFrame(
            {
                "wavelength": (wl := self.wavelength.value),
                "flux": self.flux.value,
                "mask": self.mask if self.mask else np.zeros(len(wl), dtype=int),
            }
        )

    def plot(self, ax=None, ylo=0.6, yhi=1.2, figsize=(10, 4), **kwargs):
        """Plot a quick look of the spectrum"

        Parameters
        ----------
        ax : `~matplotlib.axes.Axes`
            A matplotlib axes object to plot into. If no axes is provided, a new one will be generated.
        ylo : float
            Y-axis lower bound
        yhi : float
            Y-axis upper bound
        figsize : tuple
            Dimensions of the figure
        label : str
            The label for plt.legend()

        Returns
        -------
        ax : `~matplotlib.axes.Axes`
            The axis to display and/or modify
        """
        if not ax:
            ax = plt.subplots(1, figsize=figsize)[1]
            ax.set_ylim(ylo, yhi)
            ax.set_xlabel(r"$\lambda \;(\AA)$")
            ax.set_ylabel("Flux")

        ax.step(self.wavelength, self.flux, where="mid", **kwargs)
        return ax
