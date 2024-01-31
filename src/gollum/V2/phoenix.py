from gollum.V2.simple_spectrum import SimpleSpectrum

import pandas as pd

from dataclasses import dataclass
from dotenv import get_key

@dataclass(slots=True)
class PHOENIXSpectrum(SimpleSpectrum):
    '''
    A container for PHOENIX synthetic spectra (Husser et al. 2013).
    
    Parameters
    ----------
    `teff: int`
        The effective temperature.
    `logg: float`
        The surface gravity.
    `Z: float`
        The iron abundance.
    `alpha: float`
        The alpha abundance of the PHOENIX model.
    `extent: tuple[float, float]`
        The wavelength limits of the PHOENIX model.
    `**kwargs` (`wavelength` & `flux`)
        Bypasses the PHOENIX constructor and directly assigns the wavelength and flux arrays. \n
        Still requires fundamental stellar parameters; if you don't need the metadata use `SimpleSpectrum` directly.
    '''
    teff: int
    logg: float
    Z: float
    alpha: float

    def __init__(self, teff: int, logg: float, Z: float, alpha: float, extent: tuple[float, float] = None, **kwargs):
        self.teff, self.logg, self.Z, self.alpha = teff, logg, Z, alpha
        if kwargs:
            super(PHOENIXSpectrum, self).__init__(**kwargs)
            return
        if not self._env().exists() or not (path := get_key(self._env(), 'PHOENIX_PATH')):
            self._env().touch()
            raise RuntimeError('PHOENIX data path not set. Set it with: PHOENIXSpectrum.configure(<path>)')
        super(PHOENIXSpectrum, self).__init__([10, 20, 30], [40, 50, 60])


x = PHOENIXSpectrum(1, 2, 3, 7, (4, 5))
y = PHOENIXSpectrum(1, 2, 3, 7, wavelength=[1, 2, 3], flux=[4, 5, 6])
print(x + y)
