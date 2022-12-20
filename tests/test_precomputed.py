from gollum.phoenix import PrecomputedSpectrum, PHOENIXSpectrum
from specutils import Spectrum1D
from pandas import DataFrame
import numpy as np
import astropy.units as u


def test_basic():
    """Do the basic methods work?"""
    assert PrecomputedSpectrum(
        spectral_axis=np.ones(5) * u.AA, flux=np.ones(5) * u.dimensionless_unscaled
    )
    spec = PHOENIXSpectrum(teff=7000, logg=4, Z=0.5, download=True)

    assert isinstance(spec.velocity_spacing, np.ndarray)
    assert np.median(spec.normalize().flux.value) == 1
    assert isinstance(spec.resample_to_uniform_in_velocity(), PHOENIXSpectrum)
    assert np.all(
        spec.resample_to_uniform_in_velocity().wavelength
        == spec.decimate(1.4).wavelength
    )
    assert isinstance(spec.get_blackbody_spectrum(), PrecomputedSpectrum)
    assert isinstance(spec.divide_by_blackbody(spec.teff), PHOENIXSpectrum)

    tilted_spec, model = spec.tilt_to_data(
        target_spectrum=PHOENIXSpectrum(teff=3000, logg=4, Z=0.5, download=True), return_model=True,
    )
    assert np.all(tilted_spec.flux == spec.multiply(model(spec.wavelength)).flux)
    assert isinstance(spec.to_pandas(), DataFrame)

    assert isinstance(spec.fit_continuum(), PrecomputedSpectrum)


def test_resampling():
    """Does resampling work?"""

    n_pix = 1000
    x_vec = np.linspace(10_000, 10_100, n_pix) * u.Angstrom
    y_continuum = np.ones_like(x_vec.value)
    y_perturbation = 0.05 * np.sin(2 * np.pi * x_vec.value / 6.0)
    y_vec = (y_continuum + y_perturbation) * u.dimensionless_unscaled
    spec = PrecomputedSpectrum(spectral_axis=x_vec, flux=y_vec)

    new_spec = spec.resample_to_uniform_in_velocity()

    assert new_spec is not None
    assert isinstance(new_spec, Spectrum1D)
    assert isinstance(new_spec.flux, np.ndarray)
    assert len(new_spec.wavelength) != len(spec.wavelength)

    new_spec = spec.resample_to_uniform_in_velocity(oversample=3.5)

    assert new_spec is not None
    assert isinstance(new_spec, Spectrum1D)
    assert isinstance(new_spec.flux, np.ndarray)
    assert len(new_spec.wavelength) != len(spec.wavelength)

    # Test undersampling... should log a warning...
    new_spec = spec.resample_to_uniform_in_velocity(oversample=0.2)

    assert new_spec is not None
    assert isinstance(new_spec, Spectrum1D)
    assert isinstance(new_spec.flux, np.ndarray)
    assert len(new_spec.wavelength) != len(spec.wavelength)
