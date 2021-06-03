r"""
PHOENIX Spectrum
-----------------

A container for a single Phoenix grid-point spectrum of wavelength and flux :math:`F(\lambda)`.

PHOENIXSpectrum
###############
"""

import warnings
import logging
import numpy as np
import astropy
from astropy.io import fits
from astropy import units as u

import matplotlib.pyplot as plt
import os
import copy

from scipy.ndimage import gaussian_filter1d


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


class PHOENIXSpectrum(Spectrum1D):
    r"""
    A container for PHOENIX spectra

    Args:
        Teff (int): The Teff label of the PHOENIX model to read in.  Must be on the PHOENIX grid.
        logg (float): The logg label of the PHOENIX model to read in.  Must be on the PHOENIX grid.
        path (str): The path to your local PHOENIX grid library.  You must have the PHOENIX
            grid downloaded locally.  Default: "~/libraries/raw/PHOENIX/"
        wl_lo (float): the bluest wavelength of the models to keep (Angstroms)
        wl_hi (float): the reddest wavelength of the models to keep (Angstroms)
        """

    def __init__(
        self, *args, teff=None, logg=None, path=None, wl_lo=8038, wl_hi=12849, **kwargs
    ):

        if path is None:
            path = "~/libraries/raw/PHOENIX/"

        if (teff is not None) & (logg is not None):
            base_path = os.path.expanduser(path)
            assert os.path.exists(
                base_path
            ), "You must specify the path to local PHOENIX models"

            wl_filename = base_path + "/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
            assert os.path.exists(
                wl_filename
            ), "You need to place the PHOENIX models in {}".format(base_path)

            wl_orig = fits.open(wl_filename)[0].data.astype(np.float64)

            mask = (wl_orig > wl_lo) & (wl_orig < wl_hi)
            wl_out = wl_orig[mask]

            fn = (
                base_path
                + "/Z-0.0/lte{:05d}-{:0.2f}-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
            ).format(teff, logg)
            assert os.path.exists(fn), "Double check that the file {} exists".format(fn)

            flux_orig = fits.open(fn)[0].data.astype(np.float64)
            # Units: erg/s/cm^2/cm
            flux_native = flux_orig[mask]

            super().__init__(
                spectral_axis=wl_out * u.AA,
                flux=flux_native * u.erg / u.s / u.cm ** 2 / u.cm,
                **kwargs,
            )

        else:
            super().__init__(*args, **kwargs)

    def normalize(self):
        """Normalize spectrum by its median value

        Returns
        -------
        normalized_spec : (PHOENIXSpectrum)
            Normalized Spectrum
        """
        median_flux = np.median(self.flux)

        return self.divide(median_flux, handle_meta="first_found")

    def rotationally_broaden(self, vsini, u1=0.0, u2=0.0):
        """Rotationally broaden the spectrum for a given vsini
        Implementation inspired by https://github.com/HajimeKawahara/exojax 

        Known limitation: If the wavelength sampling changes with wavelength, 
          the convolution becomes inaccurate.  It may be better to FFT,
          following Starfish.

        Args:
            vsini: V sini for rotation 
            u1: Limb-darkening coefficient 1
            u2: Limb-darkening coefficient 2
            
        Returns
        -------
        broadened_spec : (PHOENIXSpectrum)
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
        kernel = kernel[kernel > 0]

        convolved_flux = np.convolve(self.flux, kernel, mode="same")
        return PHOENIXSpectrum(spectral_axis=self.wavelength, flux=convolved_flux)

    def instrumental_broaden(self, resolving_power=55_000):
        """Instrumentally broaden the spectrum for a given instrumental resolution, R

        Known limitation: If the wavelength sampling changes with wavelength, 
          the convolution becomes inaccurate.  It may be better to FFT,
          following Starfish.

        Args:
            resolving_power: Instrumental resolving power :math:`R\equiv \frac{\lambda}{\delta \lambda}` 
            
        Returns
        -------
        broadened_spec : (PHOENIXSpectrum)
            Instrumentally Broadened Spectrum
        """
        angstroms_per_pixel = np.median(np.diff(self.wavelength.value))
        lam0 = np.median(self.wavelength.value)
        delta_lam = lam0 / resolving_power

        scale_factor = 2.355  # Todo is this sigma or FWHM?
        sigma = delta_lam * scale_factor / angstroms_per_pixel

        convolved_flux = gaussian_filter1d(self.flux.value, sigma) * self.flux.unit
        return PHOENIXSpectrum(spectral_axis=self.wavelength, flux=convolved_flux)

    def rv_shift(self, rv=0):
        """Shift the spectrum by a radial velocity, in units of km/s

        Args:
            rv: Radial velocity in units of km/s
            
        Returns
        -------
        shifted_spec : (PHOENIXSpectrum)
            RV Shifted Spectrum
        """
        self.radial_velocity = rv * u.km / u.s
        return self

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
