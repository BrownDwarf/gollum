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
from specutils import SpectrumCollection
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
    A container for a single Sonora precomputed synthetic spectrum of a brown dwarfs or free-floating 
    Gas Giant planet. 

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
            ), "You must specify the path to local Sonora models: {}".format(base_path)
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


class SonoraGrid(SpectrumCollection):
    r"""
    A container for a grid of Sonora precomputed synthetic spectra of brown dwarfs and free-floating 
    Gas Giant planets.  

    Args:
        Teff_range (tuple): The Teff limits of the grid model to read in.
        logg (tuple): The logg limits of the Sonora model to read in.
        path (str): The path to your local Sonora grid library.  
            You must have the Sonora grid downloaded locally.  
            Default: "~/libraries/raw/Sonora/"
        wl_lo (float): the bluest wavelength of the models to keep (Angstroms)
        wl_hi (float): the reddest wavelength of the models to keep (Angstroms)
        """

    def __init__(
        self, teff_range=None, logg_range=None, path=None, wl_lo=8038, wl_hi=12849,
    ):
        teff_points = np.hstack(
            (
                np.arange(500, 600, 25),
                np.arange(600, 1000, 50),
                np.arange(1000, 2401, 100),
            )
        )
        logg_points = np.arange(4.0, 5.51, 0.25)

        if teff_range is not None:
            subset = (teff_points > teff_range[0]) & (teff_points < teff_range[1])
            teff_points = teff_points[subset]

        if logg_range is not None:
            subset = (logg_points > logg_range[0]) & (logg_points < logg_range[1])
            logg_points = logg_points[subset]

        wavelengths, fluxes = [], []
        for teff in teff_points:
            for logg in logg_points:
                spec = SonoraSpectrum(
                    teff=teff, logg=logg, path=path, wl_lo=wl_lo, wl_hi=wl_hi
                )
                wavelengths.append(spec.wavelength)
                fluxes.append(spec.flux)
        flux_out = np.array(fluxes) * fluxes[0].unit
        wave_out = np.array(wavelengths) * wavelengths[0].unit
        super().__init__(flux=flux_out, spectral_axis=wave_out)

    def __getitem__(self, key):
        flux = self.flux[key]
        if flux.ndim != 1:
            raise ValueError(
                "Currently only 1D data structures may be "
                "returned from slice operations."
            )
        spectral_axis = self.spectral_axis[key]
        uncertainty = None if self.uncertainty is None else self.uncertainty[key]
        wcs = None if self.wcs is None else self.wcs[key]
        mask = None if self.mask is None else self.mask[key]
        if self.meta is None:
            meta = None
        else:
            try:
                meta = self.meta[key]
            except KeyError:
                meta = self.meta

        return SonoraSpectrum(
            flux=flux,
            spectral_axis=spectral_axis,
            uncertainty=uncertainty,
            wcs=wcs,
            mask=mask,
            meta=meta,
        )

