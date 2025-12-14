"""
wavelength.py

Author  : Kévin Walcarius
Date    : 2025-11-19
Version : 1.0
License : MIT
Summary : Base Wavelength class with pydentic, ensures ordering and
          positivity of the user input wavelength.
"""

from typing import List

from pydantic import BaseModel

from pymopsmap.utils import SortedPosFloat64List


# =================================================================================================
# Base class for wavelength
# =================================================================================================
class Wavelength(BaseModel):
    """
    A simple model for wavelengths in MOPSMAP. Type checking
    is performed with pydantic.
    """

    values: SortedPosFloat64List

    @property
    def command(self) -> str:
        if len(self.values) == 1:
            return f"wavelength {self.values[0]}"
        return "wavelength list " + " ".join(str(v) for v in self.values)


# =================================================================================================
# Wrapper
# =================================================================================================
def wavelength(wl: Wavelength | float | List[float]) -> Wavelength:
    """
    Convenience wrapper for creating a :class:`Wavelength` object.

    Parameters
    ----------
    wl : float | list[float] | Wavelength
        Input wavelength(s).

    Returns
    -------
    Wavelength
        Normalized wavelength object.
    """
    if isinstance(wl, Wavelength):
        return wl
    return Wavelength(values=wl)
