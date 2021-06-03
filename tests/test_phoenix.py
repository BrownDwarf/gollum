# from astropy.nddata.ccddata import _uncertainty_unit_equivalent_to_parent
import pytest
import time
from gollum.phoenix import PHOENIXSpectrum
from specutils import Spectrum1D

# from astropy.nddata.nduncertainty import StdDevUncertainty
import numpy as np
import glob
import astropy


def test_basic():
    """Do the basic methods work?"""

    spec = PHOENIXSpectrum(teff=5000, logg=4)

    assert spec is not None
    assert isinstance(spec, Spectrum1D)
    assert isinstance(spec.flux, np.ndarray)
    assert len(spec.flux) == len(spec.wavelength)

    new_spec = spec.normalize()

    assert new_spec.shape[0] == spec.shape[0]
    assert np.median(new_spec.flux) == 1

    ax = new_spec.plot(label="demo", color="r")
    assert ax is not None
