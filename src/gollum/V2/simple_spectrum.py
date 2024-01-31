import matplotlib.pyplot as plt
import numpy as np

from dataclasses import dataclass, replace
from dotenv import set_key
from inspect import getfile
from numpy.typing import ArrayLike
from pathlib import Path
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import gaussian_filter1d
from typing import Self

@dataclass(slots=True)
class SimpleSpectrum:
    '''
    A minimal container for storing spectra.
    
    Parameters
    ----------
    `wavelength: ArrayLike`
        The spectral axis.
    `flux: ArrayLike`
        The spectral flux density.
    '''
    wavelength: ArrayLike
    flux: ArrayLike

    def __init__(self, wavelength: ArrayLike, flux: ArrayLike):
        if len(wavelength) != len(flux):
            raise ValueError("Wavelength and flux must have the same length.")
        self.wavelength = np.array(wavelength, dtype=np.float64)
        self.flux = np.array(flux, dtype=np.float64)
    
    def __add__(self, other: Self):
        return replace(self, flux=self.flux + other.flux)
    
    def __sub__(self, other: Self):
        return replace(self, flux=self.flux - other.flux)
    
    def __mul__(self, other: Self):
        return replace(self, flux=self.flux * other.flux)
    
    def __truediv__(self, other: Self):
        return replace(self, flux=self.flux / other.flux)
    
    def __iadd__(self, other: Self):
        self.flux += other.flux
        return self
    
    def __isub__(self, other: Self):
        self.flux -= other.flux
        return self
    
    def __imul__(self, other: Self):
        self.flux *= other.flux
        return self
    
    def __itruediv__(self, other: Self):
        self.flux /= other.flux
        return self
    
    def normalize(self, percentile: float):
        '''
        [TRANSFORM] Normalize flux to a given percentile.
        
        Parameters
        ----------
        `percentile: float`
            The percentile at which the normalized flux is to equal 1.
        '''
        return replace(self, flux=self.flux/np.percentile(self.flux, percentile))
    
    def rsmooth(self, vsini: float):
        '''
        [TRANSFORM] Apply simple rotational broadening to the spectrum.
        
        Parameters
        ----------
        `vsini: float`
            The projected rotational velocity in km/s.
        '''
        mid = np.median(self.wavelength)
        x = 299792 * (self.wavelength - mid) / (mid * vsini)
        conv = np.nan_to_num(np.sqrt(1 - x*x))
        return replace(self, flux=np.convolve(self.flux, conv[conv > 0]/conv.sum(), mode='same'))
    
    def ismooth(self, R: int):
        '''
        [TRANSFORM] Apply instrumental broadening to the spectrum.
        
        Parameters
        ----------
        `R: int`
            The resolving power of the instrument.
        '''
        return replace(self, flux=gaussian_filter1d(self.flux, np.median(self.wavelength) / (2.355 * R * np.median(np.diff(self.wavelength)))))
    
    def shift(self, rv: float):
        '''
        [TRANSFORM] Apply a non-relativistic Doppler shift to the spectrum.
        
        Parameters
        ----------
        `rv: float`
            The radial velocity in km/s.
        '''
        return replace(self, wavelength=self.wavelength * (1 + rv / 299792))

    def resample(self, spectral_axis: ArrayLike):
        '''
        [TRANSFORM] Resample the spectrum onto a new spectral axis.

        Parameters
        ----------
        `spectral_axis: ArrayLike`
            The new spectral axis in Angstroms.
        '''
        return replace(self, wavelength=spectral_axis, flux=PchipInterpolator(self.wavelength, self.flux)(spectral_axis))

    def view(self, **kwargs):
        '''
        [VISUALIZE] Plot a quick view of the spectrum.
        
        Parameters
        ----------
        `**kwargs`
            Passthrough for matplotlib.pyplot.plot.
        '''

        plt.figure(figsize=(12, 6))
        plt.plot(self.wavelength, self.flux, **kwargs)
        plt.xlabel('Wavelength')
        plt.ylabel('Spectral Flux Density')
        plt.tight_layout()
        plt.show()
    
    @staticmethod
    def blackbody(spectral_axis: ArrayLike, T: float):
        '''
        [ACCESS] Create a blackbody spectrum at a given temperature with units of W/m^2/nm.

        Parameters
        ----------
        `spectral_axis: ArrayLike`
            The spectral axis in Angstroms.
        `T: float`
            The temperature of the blackbody in K.
        '''
        return SimpleSpectrum(spectral_axis, 1.191043E25 / (spectral_axis**5 * (np.exp(1.43877688E8 / (T * spectral_axis)) - 1)))
    
    @staticmethod
    def _env():
        '''
        [ACCESS] The path to gollum's environment file.
        '''
        return Path(getfile(SimpleSpectrum)).parent / '.env'
    
    @classmethod
    def configure(cls, data_path: str):
        '''
        [SETUP] Configure gollum's environment for a model grid.

        Parameters
        ----------
        `data_path: str`
            The path to the data directory.
        '''
        set_key(cls._env(), cls.__name__.replace('Spectrum', '_PATH'), str(Path(data_path).expanduser().resolve()))
