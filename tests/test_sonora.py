from gollum.precomputed_spectrum import PrecomputedSpectrum
from gollum.sonora import SonoraGrid, SonoraSpectrum
from specutils import Spectrum1D
import numpy as np
import astropy.units as u
from specutils.spectra.spectrum_collection import SpectrumCollection


def test_basic():
    """Do the basic methods work?"""

    spec = SonoraSpectrum(teff=1200, logg=4.5)

    assert spec is not None
    assert isinstance(spec, Spectrum1D)
    assert isinstance(spec, PrecomputedSpectrum)
    assert isinstance(spec, SonoraSpectrum)
    assert isinstance(spec.flux, np.ndarray)
    assert len(spec.flux) == len(spec.wavelength)

    new_spec = spec.normalize()

    assert new_spec.shape[0] == spec.shape[0]
    assert np.median(new_spec.flux) == 1

    ax = new_spec.plot(label="demo", color="r")
    assert ax

    new_spec = (
        spec.rotationally_broaden(28.8).rv_shift(10.1).instrumental_broaden(55_000)
    )

    assert new_spec
    assert isinstance(new_spec, Spectrum1D)
    assert isinstance(new_spec.flux, np.ndarray)
    assert len(new_spec.flux) == len(new_spec.wavelength)
    assert len(new_spec.flux) == len(spec.wavelength)
    assert isinstance(spec, SonoraSpectrum)


def test_resample():
    """Do the basic methods work?"""

    spec = SonoraSpectrum(teff=850, logg=4.25)

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
    assert isinstance(spec, SonoraSpectrum)


def test_grid():
    """Do the basic methods work?"""

    grid = SonoraGrid(
        teff_range=(950, 1020),
        logg_range=(4.25, 4.75),
        metallicity_range=(0.0, 0.5),
    )

    assert grid
    assert isinstance(grid, SpectrumCollection)
    assert grid.nspectral > 0

    assert isinstance(grid[2], SonoraSpectrum)
    assert hasattr(grid[3], "rotationally_broaden")
    assert hasattr(grid, "teff_points")
    assert hasattr(grid, "grid_labels")
    assert hasattr(grid, "grid_points")

    teff_point = grid.find_nearest_teff(1001.5)
    assert teff_point == 1000
    teff_point = grid.find_nearest_teff(0)
    assert teff_point == grid.teff_points.min()
    teff_point = grid.find_nearest_teff(1e6)
    assert teff_point == grid.teff_points.max()

    newgrid = grid.truncate(wavelength_range=(10800.0 * u.Angstrom, 11800 * u.Angstrom))
    assert newgrid
    assert isinstance(newgrid, SpectrumCollection)
    assert newgrid.grid_points == grid.grid_points
    assert len(newgrid[0].wavelength.value) < len(grid[0].wavelength.value)
    assert len(newgrid) == len(grid)


def test_nearest_gridpoint():
    """Do the basic methods work?"""

    grid = SonoraGrid(
        teff_range=(950, 1200),
        logg_range=(4.25, 4.75),
        metallicity_range=(0.0, 0.5),
    )

    assert grid

    # Test if get_nearest_grid_point returns a valid point (itself)
    assert grid.find_nearest_grid_point(1100, 4.5, 0.0) == (1100, 4.5, 0.0)

    # Test if get_nearest_grid_point return a valid point for an invalid input
    assert grid.find_nearest_grid_point(1100.1, 4.501, 0.001) == (1100, 4.5, 0.0)
