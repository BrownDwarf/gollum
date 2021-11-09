from gollum.phoenix import PrecomputedSpectrum
from specutils import Spectrum1D
import numpy as np
import astropy.units as u


def test_basic():
    """Do the basic methods work?"""

    n_pix = 1000
    x_vec = np.linspace(10_000, 10_100, n_pix) * u.Angstrom
    y_continuum = np.ones_like(x_vec.value)
    y_perturbation = 0.05 * np.sin(2 * np.pi * x_vec.value / 6.0)
    y_vec = (y_continuum + y_perturbation) * u.dimensionless_unscaled
    spec = PrecomputedSpectrum(spectral_axis=x_vec, flux=y_vec)

    assert spec is not None
    assert isinstance(spec, Spectrum1D)
    assert isinstance(spec, PrecomputedSpectrum)
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
    assert isinstance(spec, PrecomputedSpectrum)

    new_spec = new_spec.fit_continuum(pixel_distance=91)

    assert new_spec is not None
    assert isinstance(new_spec, Spectrum1D)
    assert isinstance(new_spec.flux, np.ndarray)
    assert len(new_spec.flux) == len(new_spec.wavelength)
    assert len(new_spec.flux) == len(spec.wavelength)
    assert isinstance(spec, PrecomputedSpectrum)
