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
from astropy.io import fits
from astropy import units as u
import specutils

import matplotlib.pyplot as plt
import os
import copy

from scipy.ndimage import gaussian_filter1d
from specutils.manipulation import LinearInterpolatedResampler


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
        if percentile is None:
            scalar_flux = np.median(self.flux)
        else:
            scalar_flux = np.percentile(self.flux, percentile)

        return self.divide(scalar_flux, handle_meta="first_found")

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
            output = copy.copy(self)
            output.radial_velocity = rv * u.km / u.s
            return self._copy(
                spectral_axis=output.wavelength.value * self.wavelength.unit
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

        return self._copy(spectral_axis=output.wavelength, flux=output.flux)

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
            ax.step(self.wavelength, self.flux, **kwargs)
        else:
            ax.step(self.wavelength, self.flux, **kwargs)

        return ax
