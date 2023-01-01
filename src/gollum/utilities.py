from copy import deepcopy
from numpy import array, floor, ceil
from specutils import Spectrum1D


def apply_numpy_mask(spec, mask):
    """Applies a numpy-style boolean mask to an input spectrum (True=Keep, False=Discard)

    Parameters
    ----------
    spec : Spectrum1D object
        Object containing a spectrum
    mask : boolean mask, typically a numpy array
        The mask to apply to the spectrum

    Returns
    -------
    masked_spec: Spectrum1D object
        The spectrum with the mask applied
    """

    assert isinstance(spec, Spectrum1D), "Input must be a specutils Spectrum1D object"
    assert mask.any(), "The masked spectrum must have at least one pixel remaining"

    if (mask_length := len(mask)) != (npx := len(spec.spectral_axis.value)):
        raise IndexError(
            f"Your mask has {mask_length} entries and your spectrum has {npx} pixels. They should be the same shape"
        )

    return spec.__class__(
        spectral_axis=spec.wavelength.value[mask] * spec.wavelength.unit,
        flux=spec.flux[mask],
        mask=spec.mask[mask] if spec.mask else None,
        wcs=None,
        meta=deepcopy(spec.meta) if spec.meta else None,
    )


def _truncate(grid, wavelength_range=None, data=None):
    """Truncate the wavelength range of the grid

    Parameters
    ----------
    wavelength_range: list or tuple
        A pair of values (with units) that denote the shortest and longest wavelengths
        for truncating the grid.
    data: Spectrum1D-like
        A spectrum to which this method will match the wavelength limits

    Returns
    -------
    truncated_spectrum: Spectrum1D-like
        The spectrum after being truncated to the given wavelength range
    """
    assert (
        bool(data) + bool(wavelength_range) == 1
    ), "Please provide only one of the following: data OR wavelength_range"
    wl_lo, wl_hi = (
        (floor(data.wavelength.min()), ceil(data.wavelength.max()))
        if data
        else wavelength_range
    )

    wavelengths, fluxes = [], []
    for spec in grid:
        mask = (spec.wavelength >= wl_lo) & (spec.wavelength <= wl_hi)
        wavelengths.append(spec.wavelength.value[mask])
        fluxes.append(spec.flux.value[mask])

    assert fluxes and wavelengths
    return grid.__class__(
        flux=array(fluxes) * grid[0].flux.unit,
        spectral_axis=array(wavelengths) * grid[0].wavelength.unit,
        meta=grid.meta,
    )
