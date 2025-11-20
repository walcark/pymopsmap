"""
wavelength.py

Author  : Kévin Walcarius
Date    : 2025-11-19
Version : 1.0
License : MIT
Summary : Base Wavelength class with pydentic, ensures ordering and
          positivity of the user input wavelength.
"""

from pydantic import BaseModel, Field, field_validator, model_validator
from typing import List, Any


# =================================================================================================
# Base class for wavelength
# =================================================================================================
class Wavelength(BaseModel):
    """
    Unified wavelength representation for MOPSMAP.

    A `Wavelength` object always stores wavelengths as an ordered
    list of floats.

    Parameters
    ----------
    values : float | list[float]
        Wavelength(s) in microns.

    Examples
    --------
    Single wavelength
    >>> Wavelength(values=0.55).to_section()
    'wavelength 0.55'

    Multiple wavelengths
    >>> Wavelength([0.44, 0.55, 0.67]).to_section()
    'wavelength list 0.44 0.55 0.67'

    Usage inside the public API
    >>> from pymopsmap.api import wavelength
    >>> wl = wavelength([0.44, 0.55, 0.67])
    >>> wl.values
    [0.44, 0.55, 0.67]
    """

    values: List[float] = Field(..., min_length=1)

    @field_validator("values", mode="before")
    @classmethod
    def coerce_scalar_or_list(cls, v: Any) -> List[float]:
        """Transform float input to list."""
        if isinstance(v, (int, float)):
            return [float(v)]
        if isinstance(v, list):
            return v
        raise TypeError("wl must be a float or a list of floats.")

    @model_validator(mode="after")
    def validate_values(self):
        # positivity
        if any(w <= 0 for w in self.values):
            raise ValueError(f"Wavelengths must be > 0 (got {self.values})")
        # sorted ascending
        if self.values != sorted(self.values):
            raise ValueError(
                f"Wavelengths must be sorted ascending (got {self.values})"
            )
        return self

    def to_section(self) -> str:
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
