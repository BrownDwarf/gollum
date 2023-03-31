r"""
coolTLUSTY Spectrum
-----------------

A container for a single grid-point of the coolTLUSTY precomputed synthetic model spectrum of brown dwarfs and free-floating Gas Giant planets.  The spectrum is a vector with coordinates wavelength and flux :math:`F(\lambda)`.

coolTLUSTYSpectrum
###############
"""

import os

from itertools import product
from tqdm import tqdm
from gollum.utilities import _truncate
from gollum.precomputed_spectrum import *
from gollum.telluric import TelluricSpectrum
from pandas import read_csv
from specutils import SpectrumCollection
from bokeh.io import show, output_notebook
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import Slider, Range1d
from bokeh.layouts import layout, Spacer
from bokeh.models.widgets import Button, Div

log = getLogger(__name__)

#  See Issue: https://github.com/astropy/specutils/issues/779
filterwarnings("ignore", category=AstropyDeprecationWarning)
# See Issue: https://github.com/astropy/specutils/issues/800
filterwarnings("ignore", category=RuntimeWarning)


class coolTLUSTYSpectrum(PrecomputedSpectrum):
    """
    A container for coolTLUSTY 2023 spectra

    Parameters
    ----------
    teff : int
        The teff label of the coolTLUSTY model to read in.  Must be on the coolTLUSTY grid.
    logg : float
        The logg label of the coolTLUSTY model to read in.  Must be on the coolTLUSTY grid.
    path : str
        The path to your local coolTLUSTY grid library.  You must have the coolTLUSTY grid downloaded locally.  Default: "~/libraries/raw/coolTLUSTY/"
    wl_lo : float
        The shortest wavelength of the models to keep (Angstroms)
    wl_hi : float
        The longest wavelength of the models to keep (Angstroms)
    """

    def __init__(
        self,
        *args,
        teff=None,
        logg=None,
        z=None,
        path="~/libraries/raw/coolTLUSTY/ClearEQ/",
        wl_lo=8038,
        wl_hi=12849,
        **kwargs,
    ):

        teff_points = np.arange(250, 601, 25)
        logg_points = np.arange(3.5, 5.01, 0.25)
        z_points = np.array([0.316,1.000, 3.160])


        if teff and logg:
            base_path = os.path.expanduser(path)
            assert os.path.exists(base_path), "Given path does not exist."
            assert teff in teff_points, "teff must be a point on the grid"
            assert logg in logg_points, "logg must be a point on the grid"
            assert z in z_points, "Fe/H must be a point on the grid"

            fn = "{}T{:3d}_g{:0.2f}_Z{:0.3f}.21".format(base_path, int(teff), logg, z)

            df_native = read_csv(fn, delim_whitespace=True, usecols=['LAMBDA(mic)', 'FLAM'])
            df_native['wavelength_um'] = df_native['LAMBDA(mic)'].str.replace('D', 'e').astype(float)
            df_native['flux'] = df_native['FLAM'].str.replace('D', 'e').astype(float)

            # convert to Angstroms
            df_native["wavelength"] = df_native["wavelength_um"] * 10000.0
            mask = (df_native.wavelength > wl_lo) & (df_native.wavelength < wl_hi)
            df_trimmed = df_native[mask].reset_index(drop=True)

            super().__init__(
                spectral_axis=df_trimmed.wavelength.values * u.AA,
                flux=df_trimmed.flux.values * u.erg / u.s / u.cm**2 / u.AA,
                **kwargs,
            )

        else:
            super().__init__(*args, **kwargs)
