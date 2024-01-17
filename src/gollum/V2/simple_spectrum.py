import numpy as np
import pandas as pd

from dataclasses import dataclass
from numpy.typing import ArrayLike
from scipy.ndimage import gaussian_filter1d
from typing import Self

@dataclass(slots=True)
class SimpleSpectrum:
    '''
    A minimal container for storing spectra.
    
    Parameters
    ----------
    wavelength: ArrayLike
        The spectral axis. Must be in nm.
    flux: ArrayLike
        The spectral flux density. Must be in W/m^2/nm.
    '''
    wavelength: ArrayLike
    flux: ArrayLike

    def __init__(self, wavelength: ArrayLike, flux: ArrayLike):
        if len(wavelength) != len(flux):
            raise ValueError("Wavelength and flux must have the same length.")
        self.wavelength = np.array(wavelength, dtype=np.float64)
        self.flux = np.array(flux, dtype=np.float64)
    
    def __add__(self, other: Self):
        return self.__class__(self.wavelength, self.flux + other.flux)
    
    def __sub__(self, other: Self):
        return self.__class__(self.wavelength, self.flux - other.flux)
    
    def __mul__(self, other: Self):
        return self.__class__(self.wavelength, self.flux * other.flux)
    
    def __truediv__(self, other: Self):
        return self.__class__(self.wavelength, self.flux / other.flux)
    
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
    
    def __eq__(self, other: Self):
        return np.allclose(self.flux, other.flux)
    
    def normalize(self, percentile: float):
        '''
        [TRANSFORM] Normalize flux to a given percentile.
        
        Parameters
        ----------
        percentile: float
            The percentile at which the normalized flux is to equal 1.
        '''
        return self.__class__(wavelength=self.wavelength, 
                                 flux=self.flux/np.percentile(self.flux, percentile))
    
    def rsmooth(self, vsini: float):
        '''
        [TRANSFORM] Apply simple rotational broadening to the spectrum.
        
        Parameters
        ----------
        vsini: float
            The projected rotational velocity in km/s.
        '''
        mid = np.median(self.wavelength)
        x = 299792 * (self.wavelength - mid) / (mid * vsini)
        conv = np.nan_to_num(np.sqrt(1 - x*x))
        return self.__class__(wavelength=self.wavelength, 
                                 flux=np.convolve(self.flux, conv[conv > 0]/conv.sum(), mode='same'))
    
    def ismooth(self, R: int):
        '''
        [TRANSFORM] Apply instrumental broadening to the spectrum.
        
        Parameters
        ----------
        R: int
            The resolving power of the instrument.
        '''
        return self.__class__(wavelength=self.wavelength, 
                                 flux=gaussian_filter1d(self.flux, np.median(self.wavelength) / 
                                                        (2.355 * R * np.median(np.diff(self.wavelength)))))
    
    