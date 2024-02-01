from gollum.V2.simple_spectrum import SimpleSpectrum

from dataclasses import dataclass
from dotenv import get_key
from pandas import read_parquet

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

    def __init__(self, teff: int = None, logg: float = None, Z: float = None, alpha: float = None, extent: tuple[float, float] = None, **kwargs):
        self.teff, self.logg, self.Z, self.alpha = point =  teff, logg, Z, alpha
        if kwargs:
            super(PHOENIXSpectrum, self).__init__(**kwargs)
            return

        path = get_key(self._env(), 'PHOENIX_PATH') or 'zenodo link'
        df = read_parquet(path, 'pyarrow', [str(point)], filters=[('index', '>=', extent[0]), ('index', '<=', extent[1])] if extent else None)
        super(PHOENIXSpectrum, self).__init__(df.index, df[str(point)])


import pandas as pd
pd.DataFrame(index=[1, 2, 3], data={'(6000, 4.5, 0.0, 0.0)': [4, 5, 6]}).to_parquet('test.parquet.gz', compression='gzip')

x = PHOENIXSpectrum(6000, 4.5, 0.0, 0.0)
print(x)

