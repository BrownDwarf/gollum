"""Centralized warning policy for gollum.

Keep this list short and specific. We only suppress known, non-actionable
third-party warnings that are outside this package's control.
"""

from warnings import filterwarnings
from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyWarning


def apply_gollum_warning_filters():
    # specutils currently emits this deprecation while gollum still subclasses Spectrum1D.
    filterwarnings(
        "ignore",
        category=AstropyDeprecationWarning,
        message=r"The Spectrum1D class is deprecated.*",
    )

    # Known astropy/specutils warnings tracked in project history.
    filterwarnings(
        "ignore",
        category=AstropyWarning,
        message=r".*Spectrum1D.*",
    )
