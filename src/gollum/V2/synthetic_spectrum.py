import numpy as np
import pandas as pd

from dataclasses import dataclass
from typing import Self

@dataclass(slots=True)
class SyntheticSpectrum:
    wavelength: np.ndarray
    flux: np.ndarray
    
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
    
    def __isub__(self, other: Self):
        self.flux -= other.flux
    
    def __imul__(self, other: Self):
        self.flux *= other.flux
    
    def __itruediv__(self, other: Self):
        self.flux /= other.flux
    
    def __eq__(self, other: Self):
        return np.allclose(self.flux, other.flux)
    
    