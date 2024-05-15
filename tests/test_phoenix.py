from warnings import filterwarnings, catch_warnings
filterwarnings("ignore", category=DeprecationWarning)
from gollum.phoenix import PHOENIXSpectrum, PHOENIXGrid
from pytest import raises
from urllib.error import URLError
from specutils import Spectrum1D
import numpy as np
import astropy.units as u
from astropy.utils.exceptions import AstropyDeprecationWarning

def test_spectrum():
    """Testing the PHOENIXSpectrum class"""
    filterwarnings("ignore", category=AstropyDeprecationWarning)
    with raises(URLError):
        PHOENIXSpectrum(teff=5, logg=2, Z=2.0, download=True)

    spec = PHOENIXSpectrum(teff=5000, logg=4, download=True)

    assert spec
    assert isinstance(spec, Spectrum1D)
    assert isinstance(spec.flux, np.ndarray)
    assert len(spec.flux) == len(spec.wavelength)

    new_spec = spec.normalize()

    assert new_spec.shape[0] == spec.shape[0]
    assert np.median(new_spec.flux) == 1

    assert new_spec.plot(label="demo", color="r")

    new_spec = (
        spec.rotationally_broaden(28.8).rv_shift(10.1).instrumental_broaden(55000)
    )

    assert new_spec
    assert isinstance(new_spec, Spectrum1D)
    assert isinstance(new_spec.flux, np.ndarray)
    assert len(new_spec.flux) == len(new_spec.wavelength) == len(spec.wavelength)
    assert isinstance(spec, PHOENIXSpectrum)

    mask = (spec.wavelength.value < 14000) & (spec.wavelength.value > 12000)
    masked_spec = spec.apply_boolean_mask(mask)
    assert masked_spec.wavelength.value.min() > 12000
    assert masked_spec.wavelength.value.max() < 14000


def test_resample():
    """Testing resampling methods"""

    spec = PHOENIXSpectrum(teff=5000, logg=4, download=True)

    assert spec

    target_wavelength = [1.0, 1.01, 1.02, 1.04, 1.05] * u.micron
    target = Spectrum1D(
        spectral_axis=target_wavelength,
        flux=np.ones(len(target_wavelength)) * u.erg / u.s,
    )
    resampled_spec = spec.resample(target)

    assert resampled_spec
    assert isinstance(resampled_spec, Spectrum1D)
    assert isinstance(resampled_spec.flux, np.ndarray)
    assert len(resampled_spec.flux) == len(target.wavelength)
    assert isinstance(spec, PHOENIXSpectrum)


def test_grid():
    """Testing the PHOENIXGrid methods"""

    with catch_warnings(action='ignore', category=DeprecationWarning):
        grid = PHOENIXGrid(
            teff_range=(5000, 5100), logg_range=(2, 2.5), Z_range=(0, 0.5), experimental=True, download=True
        )
    assert grid
    assert len(grid)
    assert isinstance(grid[0], PHOENIXSpectrum)
    assert grid.find_nearest_teff(5080) == 5100
    assert grid.find_nearest_Z(0.4) == 0.5
    assert grid.find_nearest_teff(6000) == 5100
    assert isinstance(
        grid.truncate(wavelength_range=(10000 * u.AA, 10500 * u.AA)), PHOENIXGrid
    )
