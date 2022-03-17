r"""
Precomputed Spectrum
-----------------

A container for a single precomputed synthetic model spectrum at a single grid-point, with wavelength and flux :math:`F(\lambda)`.

PrecomputedSpectrum
###############
"""

import warnings
import logging
import numpy as np
import astropy
from astropy import units as u
from astropy.modeling.physical_models import BlackBody
import specutils

import matplotlib.pyplot as plt
import copy

from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import UnivariateSpline, interp1d
from scipy.signal import find_peaks, wavelets
from specutils.manipulation import LinearInterpolatedResampler
from specutils.fitting import fit_generic_continuum
from gollum.utilities import apply_numpy_mask


log = logging.getLogger(__name__)

#  See Issue: https://github.com/astropy/specutils/issues/779
warnings.filterwarnings(
    "ignore", category=astropy.utils.exceptions.AstropyDeprecationWarning
)
# See Issue: https://github.com/astropy/specutils/issues/800
warnings.filterwarnings("ignore", category=RuntimeWarning)

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from specutils import Spectrum1D
    from specutils import SpectrumList


class PrecomputedSpectrum(Spectrum1D):
    r"""
    An abstract container base class for a Precomputed spectrum

    """

    def __init__(self, *args, **kwargs):

        # Todo, we could put the wavelength limits in here.

        super().__init__(*args, **kwargs)

    def apply_boolean_mask(self, mask):
        """Apply a boolean mask to the spectrum

        Parameters
        ----------
        mask: boolean mask, typically a numpy array
            The boolean mask with numpy-style masking: True means "keep" that index and
            False means discard that index
        """

        spec = apply_numpy_mask(self, mask)

        return spec

    def normalize(self, percentile=None):
        """Normalize spectrum by its median value

        Args:
            percentile: The percentile to which the spectrum will be normalized
                (default: 50th percentile)

        Returns
        -------
        normalized_spec : (PrecomputedSpectrum)
            Normalized Spectrum
        """
        spec = self._copy(
            spectral_axis=self.wavelength.value * self.wavelength.unit, wcs=None
        )
        if percentile is None:
            scalar_flux = np.nanmedian(spec.flux.value) * spec.flux.unit
        else:
            scalar_flux = np.nanpercentile(spec.flux.value, percentile) * spec.flux.unit

        return spec.divide(scalar_flux, handle_meta="first_found")

    def rotationally_broaden(self, vsini, u1=0.0, u2=0.0):
        r"""Rotationally broaden the spectrum for a given :math:`v\sin{i}`
        Implementation inspired by https://github.com/HajimeKawahara/exojax

        Known limitation: If the wavelength sampling changes with wavelength,
          the convolution becomes inaccurate.  It may be better to FFT,
          following Starfish.

        Args:
            vsini: :math:`v\sin{i}` in units of km/s
            u1: Limb-darkening coefficient 1
            u2: Limb-darkening coefficient 2

        Returns
        -------
        broadened_spec : (PrecomputedSpectrum)
            Rotationally Broadened Spectrum
        """
        lam0 = np.median(self.wavelength.value)
        velocity_grid = 299792.458 * (self.wavelength.value - lam0) / lam0
        x = velocity_grid / vsini
        x2 = x * x
        kernel = np.where(
            x2 < 1.0,
            np.pi / 2.0 * u1 * (1.0 - x2)
            - 2.0 / 3.0 * np.sqrt(1.0 - x2) * (-3.0 + 3.0 * u1 + u2 * 2.0 * u2 * x2),
            0.0,
        )
        kernel = kernel / np.sum(kernel, axis=0)
        positive_elements = kernel > 0
        if positive_elements.any():
            kernel = kernel[positive_elements]
            convolved_flux = (
                np.convolve(self.flux.value, kernel, mode="same") * self.flux.unit
            )
            return self._copy(flux=convolved_flux)
        else:
            return self

    def instrumental_broaden(self, resolving_power=55_000):
        r"""Instrumentally broaden the spectrum for a given instrumental resolution, R

        Known limitation: If the wavelength sampling changes with wavelength,
          the convolution becomes inaccurate.  It may be better to FFT,
          following Starfish.

        Args:
            resolving_power: Instrumental resolving power :math:`R = \frac{\lambda}{\delta \lambda}`

        Returns
        -------
        broadened_spec : (PrecomputedSpectrum)
            Instrumentally Broadened Spectrum
        """
        # TODO: I think we want the Nadarya-Watson estimator here instead
        angstroms_per_pixel = np.median(np.diff(self.wavelength.value))
        lam0 = np.median(self.wavelength.value)
        delta_lam = lam0 / resolving_power

        scale_factor = 2.355
        sigma = delta_lam / scale_factor / angstroms_per_pixel

        convolved_flux = gaussian_filter1d(self.flux.value, sigma) * self.flux.unit
        return self._copy(flux=convolved_flux)

    def rv_shift(self, rv):
        """Shift the spectrum by a radial velocity, in units of km/s

        Args:
            rv: Radial velocity in units of km/s

        Returns
        -------
        shifted_spec : (PrecomputedSpectrum)
            RV Shifted Spectrum
        """
        try:
            output = copy.deepcopy(self)
            output.radial_velocity = rv * u.km / u.s
            return self._copy(
                spectral_axis=output.wavelength.value * output.wavelength.unit, wcs=None
            )
        except:
            log.error(
                "rv shift requires specutils version >= 1.2, you have: {}".format(
                    specutils.__version__
                )
            )
            raise

    def resample(self, target_spectrum):
        """Resample spectrum at the wavelength points of the other spectrum

        Args:
            target_spectrum: A Spectrum1D spectrum whose wavelength grid you seek to match

        Returns
        -------
        resampled_spec : (PrecomputedSpectrum)
            Resampled spectrum
        """
        fluxc_resample = LinearInterpolatedResampler()
        output = fluxc_resample(self, target_spectrum.wavelength)

        return self._copy(
            spectral_axis=output.wavelength.value * output.wavelength.unit,
            flux=output.flux,
            wcs=None,
        )

    def get_blackbody_spectrum(self, teff=None):
        """Get the blackbody spectrum associated with the input model"""
        if teff is None:
            # Try to look for a teff attribute
            if hasattr(self, "teff"):
                if self.teff is not None:
                    teff = self.teff
                else:
                    raise NotImplementedError(
                        "Your subclass may not have implemented the teff attribute yet."
                    )

        blackbody_model = BlackBody(temperature=teff * u.Kelvin)
        blackbody_flux_per_sr = blackbody_model(self.wavelength)

        blackbody_flux_per_Hz = blackbody_flux_per_sr * np.pi * u.steradian

        flux_unit = self.flux.unit
        normalize = False
        if flux_unit == u.dimensionless_unscaled:
            if "native_flux_unit" in self.meta.keys():
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
        if normalize:
            return output.normalize()
        else:
            return output

    def divide_by_blackbody(self, teff=None):
        """Divide the spectrum by a blackbody

        Args:
            target_spectrum: Spectrum1D
                A Spectrum1D spectrum whose flux tilt you seek to match
            return_model: (bool)
                Whether or not to return the model

        Returns
        -------
        tilted_spec : (PrecomputedSpectrum)
            Tilted spectrum
        """
        return self.divide(self.get_blackbody_spectrum(teff=teff), handle_meta="ff")

    def tilt_to_data(self, target_spectrum, return_model=False):
        """Tilt the template towards a data spectrum by fitting and dividing by a low-order polynomial

        Args:
            target_spectrum: Spectrum1D
                A Spectrum1D spectrum whose flux tilt you seek to match
            return_model: (bool)
                Whether or not to return the model

        Returns
        -------
        tilted_spec : (PrecomputedSpectrum)
            Tilted spectrum
        """
        ## Assume the template is resampled exactly to the data...
        ratio = target_spectrum.divide(self, handle_meta="ff")
        g1_fit = fit_generic_continuum(ratio)
        y_continuum_fitted = g1_fit(target_spectrum.wavelength)

        tilted_spectrum = self.multiply(y_continuum_fitted)

        if return_model:
            return (tilted_spectrum, g1_fit)
        else:
            return tilted_spectrum

    def fit_continuum(
        self, pixel_distance=5001, polyorder=3, return_coeffs=False,
    ):
        """Finds the low frequency continuum trend using scipy's find_peaks filter
        and linear algebra.

        Parameters
        ----------
        pixel_distance : int
            The minimum separation between peaks, in pixels.  Default = 5001 pixels
        polyorder : int
            The order of the polynomial used to fit the peaks
        return_coeffs : bool
            If `True`, the method will return a tuple of two elements
            (trend_spec, trend_coeffs) where trend_spec is the fitted trend.
        Returns
        -------
        continuum_spec : `PrecomputedSpectrum`
            New Spectrum object representing a fit of the continuum.
        If ``return_coeffs`` is set to ``True``, this method will also return:
        trend_coeffs : np.array
            New vector of polynomial coefficients containing to reproduce the trend.
        """

        x_vector = self.wavelength.value
        y_vector = self.flux.value

        if pixel_distance > len(x_vector):
            raise IndexError(
                "Your pixel_distance is larger than the spectrum."
                " Please provide a smaller pixel distance."
            )
        peak_inds, _ = find_peaks(y_vector, distance=pixel_distance)

        x_peaks = x_vector[peak_inds]
        y_peaks = y_vector[peak_inds]

        A_matrix = np.vander(x_peaks, polyorder)
        A_full = np.vander(x_vector, polyorder)

        ATA = np.dot(A_matrix.T, A_matrix)

        coeffs = np.linalg.solve(ATA, np.dot(A_matrix.T, y_peaks))
        y_full = np.dot(coeffs, A_full.T)

        spec_out = self._copy(flux=y_full * self.flux.unit)

        if return_coeffs:
            return (spec_out, coeffs)
        else:
            return spec_out

    def to_pandas(self):
        """Export the spectrum to a pandas dataframe"""
        try:
            import pandas as pd
        except ImportError:
            log.error("The to_pandas method requires the optional dependency pandas")
        if self.mask is not None:
            mask = self.mask
        else:
            mask = np.zeros(len(self.wavelength.value), dtype=int)

        return pd.DataFrame(
            {
                "wavelength": self.wavelength.value,
                "flux": self.flux.value,
                "mask": mask,
            }
        )

    def plot(self, ax=None, ylo=0.6, yhi=1.2, figsize=(10, 4), **kwargs):
        """Plot a quick look of the spectrum"

        Parameters
        ----------
        ax : `~matplotlib.axes.Axes`
            A matplotlib axes object to plot into. If no axes is provided,
            a new one will be generated.
        ylo : scalar
            Lower limit of the y axis
        yhi : scalar
            Upper limit of the y axis
        figsize : tuple
            The figure size for the plot
        label : str
            The legend label to for plt.legend()

        Returns
        -------
        ax : (`~matplotlib.axes.Axes`)
            The axis to display and/or modify
        """
        if ax is None:
            fig, ax = plt.subplots(1, figsize=figsize)
            ax.set_ylim(ylo, yhi)
            ax.set_xlabel(r"$\lambda \;(\AA)$")
            ax.set_ylabel("Flux")
            ax.step(self.wavelength, self.flux, where="mid", **kwargs)
        else:
            ax.step(self.wavelength, self.flux, where="mid", **kwargs)

        return ax
