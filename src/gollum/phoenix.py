r"""
PHOENIX Spectrum
-----------------

A container for a single Phoenix grid-point spectrum of wavelength and flux :math:`F(\lambda)`.

PHOENIXSpectrum
###############
"""

import warnings
import logging
from gollum.precomputed_spectrum import PrecomputedSpectrum
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


class PHOENIXSpectrum(PrecomputedSpectrum):
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
