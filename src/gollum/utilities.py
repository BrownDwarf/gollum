from copy import deepcopy
from numpy import array
from specutils.spectra import Spectrum1D


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
        A pair of values that denote the shortest and longest wavelengths
        for truncating the grid.
    data: Spectrum1D-like
        A spectrum to which this method will match the wavelength limits
        
    Returns
    -------
    truncated_spectrum: Spectrum1D-like
        The spectrum after being truncated to the given wavelength range
    """
    fiducial_spec = deepcopy(grid[0])
    wavelength_units = fiducial_spec.wavelength.unit
    flux_units = fiducial_spec.flux.unit

    if data and not wavelength_range:
        wavelength_range = (
            fiducial_spec.wavelength.value.min() * wavelength_units,
            fiducial_spec.wavelength.value.max() * wavelength_units,
        )

    shortest_wavelength, longest_wavelength = wavelength_range

    wavelengths, fluxes = [], []
    for spectrum in grid:
        mask = (spectrum.wavelength > shortest_wavelength) & (
            spectrum.wavelength < longest_wavelength
        )
        wavelengths.append(spectrum.wavelength.value[mask])
        fluxes.append(spectrum.flux.value[mask])

    assert fluxes and wavelengths

    return grid.__class__(
        flux=array(fluxes) * flux_units,
        spectral_axis=array(wavelengths) * wavelength_units,
        meta=grid.meta,
    )
