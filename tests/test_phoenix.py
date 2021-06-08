from gollum.phoenix import PHOENIXSpectrum
from specutils import Spectrum1D
import numpy as np
import astropy.units as u


def test_basic():
    """Do the basic methods work?"""

    spec = PHOENIXSpectrum(teff=5000, logg=4)

    assert spec is not None
    assert isinstance(spec, Spectrum1D)
    assert isinstance(spec, PHOENIXSpectrum)
    assert isinstance(spec.flux, np.ndarray)
    assert len(spec.flux) == len(spec.wavelength)

    new_spec = spec.normalize()

    assert new_spec.shape[0] == spec.shape[0]
    assert np.median(new_spec.flux) == 1

    ax = new_spec.plot(label="demo", color="r")
    assert ax is not None

    new_spec = (
        spec.rotationally_broaden(28.8).rv_shift(10.1).instrumental_broaden(55_000)
    )

    assert new_spec is not None
    assert isinstance(new_spec, Spectrum1D)
    assert isinstance(new_spec.flux, np.ndarray)
    assert len(new_spec.flux) == len(new_spec.wavelength)
    assert len(new_spec.flux) == len(spec.wavelength)
    assert isinstance(spec, PHOENIXSpectrum)


def test_resample():
    """Do the basic methods work?"""

    spec = PHOENIXSpectrum(teff=5000, logg=4)

    assert spec is not None

    target_wavelength = [1.0, 1.01, 1.02, 1.04, 1.05] * u.micron
    target = Spectrum1D(
        spectral_axis=target_wavelength,
        flux=np.ones(len(target_wavelength)) * u.erg / u.s,
    )
    resampled_spec = spec.resample(target)

    assert resampled_spec is not None
    assert isinstance(resampled_spec, Spectrum1D)
    assert isinstance(resampled_spec.flux, np.ndarray)
    assert len(resampled_spec.flux) == len(target.wavelength)
    assert isinstance(spec, PHOENIXSpectrum)
