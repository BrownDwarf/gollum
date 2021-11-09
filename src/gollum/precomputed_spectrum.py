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
from scipy.signal import savgol_filter
from specutils.manipulation import LinearInterpolatedResampler
from specutils.fitting import fit_generic_continuum


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

    def flatten(
        self,
        window_length=101,
        polyorder=2,
        return_trend=False,
        break_tolerance=5,
        niters=3,
        sigma=3,
        mask=None,
        **kwargs,
    ):
        """Removes the low frequency trend using scipy's Savitzky-Golay filter.
        This method wraps `scipy.signal.savgol_filter`.  Abridged from the
        `lightkurve` method with the same name for flux time series, and mirrored
        from the `muler` method with the same name.

        Parameters
        ----------
        window_length : int
            The length of the filter window (i.e. the number of coefficients).
            ``window_length`` must be a positive odd integer.
        polyorder : int
            The order of the polynomial used to fit the samples. ``polyorder``
            must be less than window_length.
        return_trend : bool
            If `True`, the method will return a tuple of two elements
            (flattened_spec, trend_spec) where trend_spec is the removed trend.
        break_tolerance : int
            If there are large gaps in wavelength, flatten will split the flux into
            several sub-spectra and apply `savgol_filter` to each
            individually. A gap is defined as a region in wavelength larger than
            `break_tolerance` times the median gap.  To disable this feature,
            set `break_tolerance` to None.
        niters : int
            Number of iterations to iteratively sigma clip and flatten. If more than one, will
            perform the flatten several times, removing outliers each time.
        sigma : int
            Number of sigma above which to remove outliers from the flatten
        mask : boolean array with length of self.wavelength
            Boolean array to mask data with before flattening. Flux values where
            mask is True will not be used to flatten the data. An interpolated
            result will be provided for these points. Use this mask to remove
            data you want to preserve, e.g. spectral regions of interest.
        **kwargs : dict
            Dictionary of arguments to be passed to `scipy.signal.savgol_filter`.
        Returns
        -------
        flatten_spec : `EchelleSpectrum`
            New light curve object with long-term trends removed.
        If ``return_trend`` is set to ``True``, this method will also return:
        trend_spec : `EchelleSpectrum`
            New light curve object containing the trend that was removed.
        """
        if mask is None:
            mask = np.ones(len(self.wavelength), dtype=bool)
        else:
            # Deep copy ensures we don't change the original.
            mask = copy.deepcopy(~mask)
        # No NaNs
        mask &= np.isfinite(self.flux)
        # No outliers
        mask &= np.nan_to_num(np.abs(self.flux - np.nanmedian(self.flux))) <= (
            np.nanstd(self.flux) * sigma
        )
        for iter in np.arange(0, niters):
            if break_tolerance is None:
                break_tolerance = np.nan
            if polyorder >= window_length:
                polyorder = window_length - 1
                log.warning(
                    "polyorder must be smaller than window_length, "
                    "using polyorder={}.".format(polyorder)
                )
            # Split the lightcurve into segments by finding large gaps in time
            dlam = self.wavelength.value[mask][1:] - self.wavelength.value[mask][0:-1]
            with warnings.catch_warnings():  # Ignore warnings due to NaNs
                warnings.simplefilter("ignore", RuntimeWarning)
                cut = np.where(dlam > break_tolerance * np.nanmedian(dlam))[0] + 1
            low = np.append([0], cut)
            high = np.append(cut, len(self.wavelength[mask]))
            # Then, apply the savgol_filter to each segment separately
            trend_signal = u.Quantity(
                np.zeros(len(self.wavelength[mask])), unit=self.flux.unit
            )
            for l, h in zip(low, high):
                # Reduce `window_length` and `polyorder` for short segments;
                # this prevents `savgol_filter` from raising an exception
                # If the segment is too short, just take the median
                if np.any([window_length > (h - l), (h - l) < break_tolerance]):
                    trend_signal[l:h] = np.nanmedian(self.flux[mask][l:h])
                else:
                    # Scipy outputs a warning here that is not useful, will be fixed in version 1.2
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", FutureWarning)
                        trsig = savgol_filter(
                            x=self.flux.value[mask][l:h],
                            window_length=window_length,
                            polyorder=polyorder,
                            **kwargs,
                        )
                        trend_signal[l:h] = u.Quantity(trsig, trend_signal.unit)
            # Ignore outliers;
            # Note that it's possible numerical noise can cause outliers...
            # If this happens you can add `1e-14` below to avoid detecting
            # outliers which are merely caused by numerical noise.
            mask1 = np.nan_to_num(np.abs(self.flux[mask] - trend_signal)) < (
                np.nanstd(self.flux[mask] - trend_signal)
                * sigma
                # + Quantity(1e-14, self.flux.unit)
            )
            f = interp1d(
                self.wavelength.value[mask][mask1],
                trend_signal[mask1],
                fill_value="extrapolate",
            )
            trend_signal = u.Quantity(f(self.wavelength.value), self.flux.unit)
            mask[mask] &= mask1

        flatten_spec = copy.deepcopy(self)
        trend_spec = self._copy(flux=trend_signal)
        with warnings.catch_warnings():
            # ignore invalid division warnings
            warnings.simplefilter("ignore", RuntimeWarning)
            flatten_spec = flatten_spec.divide(trend_spec, handle_meta="ff")
        if return_trend:
            return flatten_spec, trend_spec
        else:
            return flatten_spec

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
