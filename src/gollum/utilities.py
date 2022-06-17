from copy import deepcopy
from specutils.spectra import Spectrum1D


def apply_numpy_mask(spec, mask):
    """Applies a numpy-style boolean mask to an input spectrum (True=Keep, False=Discard)

    Parameters
    ----------
    spec : Spectrum1D object
           Object containing a spectrum
    mask : boolean mask, typically a numpy array
           The boolean mask to apply to the spectrum
    
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
