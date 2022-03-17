import copy
from specutils.spectra import Spectrum1D


def apply_numpy_mask(spec, mask):
    """Applies a boolean mask to an input spectrum, numpy-style (True=Keep, False=Discard)


    Parameters
    ----------
    spec: Spectrum1D-like object
        Object storing spectrum
    mask: boolean mask, typically a numpy array
        The boolean mask with numpy-style masking: True means "keep" that index and False means discard that index
    """

    assert isinstance(spec, Spectrum1D), "Input must be a specutils Spectrum1D object"

    assert mask.sum() > 0, "The masked spectrum must have at least one pixel remaining"

    if len(mask) != len(spec.spectral_axis.value):
        raise IndexError(
            "Your boolean mask has {} entries and your spectrum has {} pixels.  "
            " The boolean mask should have the same shape as the spectrum."
        )

    if spec.mask is not None:
        mask_out = spec.mask[mask]
    else:
        mask_out = None

    if spec.meta is not None:
        meta_out = copy.deepcopy(spec.meta)
    else:
        meta_out = None

    return spec.__class__(
        spectral_axis=spec.wavelength.value[mask] * spec.wavelength.unit,
        flux=spec.flux[mask],
        mask=mask_out,
        wcs=None,
        meta=meta_out,
    )
