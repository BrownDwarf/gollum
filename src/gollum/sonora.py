r"""
Sonora Spectrum
-----------------

A container for a single grid-point of the Sonora precomputed synthetic model spectrum of brown dwarfs and free-floating Gas Giant planets.  The spectrum is a vector with coordinates wavelength and flux :math:`F(\lambda)`.

SonoraSpectrum
###############
"""

import warnings
import logging
from gollum.precomputed_spectrum import PrecomputedSpectrum
import numpy as np
import astropy
import pandas as pd
from astropy import units as u
import os

log = logging.getLogger(__name__)

#  See Issue: https://github.com/astropy/specutils/issues/779
warnings.filterwarnings(
    "ignore", category=astropy.utils.exceptions.AstropyDeprecationWarning
)
# See Issue: https://github.com/astropy/specutils/issues/800
warnings.filterwarnings("ignore", category=RuntimeWarning)


class SonoraSpectrum(PrecomputedSpectrum):
    r"""
    A container for Sonora precomputed synthetic spectra of brown dwarfs and free-floating 
    Gas Giant planets.  This 

    Args:
        Teff (int): The Teff label of the Sonora model to read in.  Must be on the Sonora grid.
        logg (float): The logg label of the Sonora model to read in.  Must be on the Sonora grid.
        path (str): The path to your local Sonora grid library.  You must have the Sonora grid downloaded locally.  Default: "~/libraries/raw/Sonora/"
        wl_lo (float): the bluest wavelength of the models to keep (Angstroms)
        wl_hi (float): the reddest wavelength of the models to keep (Angstroms)
        """

    def __init__(
        self, *args, teff=None, logg=None, path=None, wl_lo=8038, wl_hi=12849, **kwargs
    ):

        teff_points = np.hstack(
            (
                np.arange(500, 600, 25),
                np.arange(600, 1000, 50),
                np.arange(1000, 2401, 100),
            )
        )
        logg_points = np.arange(4.0, 5.51, 0.25)

        # Map logg (cgs) to the gravity labels used in file names
        logg_par_dict = {
            4.0: "100",
            4.25: "178",
            4.5: "316",
            4.75: "562",
            5.0: "1000",
            5.25: "1780",
            5.5: "3160",
        }

        if path is None:
            path = "~/libraries/raw/Sonora/"

        if (teff is not None) & (logg is not None):
            base_path = os.path.expanduser(path)
            assert os.path.exists(
                base_path
            ), "You must specify the path to local Sonora models"
            assert teff in teff_points, "Teff must be on the grid points"
            assert logg in logg_points, "logg must be on the grid points"

            base_name = "sp_t{0:0>.0f}g{1:}nc_m0.0".format(
                float(teff), logg_par_dict[logg]
            )
            fn = base_path + "/" + base_name + ".gz"

            assert os.path.exists(fn), "Double check that the file {} exists".format(fn)

            # Units: micron, erg/cm^2/s/Hz
            df_native = (
                pd.read_csv(
                    fn,
                    skiprows=[0, 1],
                    delim_whitespace=True,
                    compression="gzip",
                    names=["wavelength_um", "flux"],
                )
                .sort_values("wavelength_um")
                .reset_index(drop=True)
            )

            # convert to Angstrom
            df_native["wavelength"] = df_native["wavelength_um"] * 10_000.0
            mask = (df_native.wavelength > wl_lo) & (df_native.wavelength < wl_hi)
            df_trimmed = df_native[mask].reset_index(drop=True)

            super().__init__(
                spectral_axis=df_trimmed.wavelength.values * u.Angstrom,
                flux=df_trimmed.flux.values * u.erg / u.s / u.cm ** 2 / u.Hz,
                **kwargs,
            )

        else:
            super().__init__(*args, **kwargs)
