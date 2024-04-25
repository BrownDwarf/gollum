from astropy.io import fits
from numpy import float64
from pandas import DataFrame
from pathlib import Path
from re import search
from tqdm import tqdm

def main():
    fluff = 'PHOENIX-ACES-AGSS-COND-2011'
    p = Path('~/data/PHOENIX').expanduser()
    wavelength = fits.open(p / f'WAVE_{fluff}.fits')[0].data.astype(float64)
    p /= fluff
    fluxes = dict()
    for f in tqdm(map(lambda x: x.name[:-39], sorted(p.glob('*/*')))):
        T = int(search(r'\d{5}', f)[0])
        GZ = search(r'(\d\.\d{2})([+-]\d\.\d)', f)
        G = float(GZ[1])
        Z = float(GZ[2])
        A = float(x[1]) if (x := search(r'Alpha=([+-]\d\.\d{2})', f)) else 0.0
        z_str = f'Z{Z:+.1f}' if Z else 'Z-0.0'
        a_str = f'.Alpha={A:+.2f}' if A else ''
        flux = fits.open(p / f'{z_str}{a_str}/{f}.{fluff}-HiRes.fits')[0].data.astype(float64)
        fluxes[str((T, G, Z, A))] = flux
    DataFrame(index=wavelength, data=fluxes, dtype='float64[pyarrow]').rename_axis('wavelength').to_parquet('PHOENIX.parquet.br', compression='brotli')

    #flux = fits.open(p / f'Z-0.0/lte02300+0.50-0.0.{fluff}-HiRes.fits')[0].data.astype(np.float64)
    #df = pd.DataFrame(index=wavelength, data={'(2300, 0.5, 0.0, 0.0)': flux}).rename_axis('wavelength')
    #df.to_parquet('test.parquet.gz', compression='gzip')

if __name__ == "__main__":
    main()

