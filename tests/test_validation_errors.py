import numpy as np
import pytest
from astropy import units as u

from gollum.coolTLUSTY import coolTLUSTYSpectrum
from gollum.phoenix import PHOENIXSpectrum
from gollum.sonora import SonoraSpectrum
from gollum.utilities import _truncate, apply_numpy_mask


def test_apply_numpy_mask_requires_spectrum1d():
    with pytest.raises(TypeError, match="Spectrum1D"):
        apply_numpy_mask(spec="not-a-spectrum", mask=np.array([True, False]))


def test_apply_numpy_mask_requires_non_empty_output():
    spec = PHOENIXSpectrum(
        spectral_axis=np.ones(3) * u.AA, flux=np.ones(3) * u.dimensionless_unscaled
    )
    with pytest.raises(ValueError, match="at least one pixel"):
        apply_numpy_mask(spec=spec, mask=np.array([False, False, False]))


def test_truncate_requires_exactly_one_selector():
    with pytest.raises(ValueError, match="only one of the following"):
        _truncate(grid=[], wavelength_range=None, data=None)


def test_phoenix_missing_path_raises_filenotfound():
    with pytest.raises(FileNotFoundError, match="Given path does not exist"):
        PHOENIXSpectrum(teff=5000, logg=4, path="/does/not/exist", download=False)


def test_sonora_missing_path_raises_filenotfound():
    with pytest.raises(FileNotFoundError, match="Given path does not exist"):
        SonoraSpectrum(teff=1000, logg=4.5, path="/does/not/exist")


def test_cooltlusty_missing_path_raises_filenotfound():
    with pytest.raises(FileNotFoundError, match="Given path does not exist"):
        coolTLUSTYSpectrum(teff=250, logg=3.5, z=0.316, path="/does/not/exist")
