"""Input adapters — convert external data formats to pymopsmap models."""

from .cams import CamsAerosol, CamsVersion, read_aerosol_microphysical_parameters
from .opac import OpacMix, OpacMixName

__all__ = [
    "CamsAerosol",
    "CamsVersion",
    "read_aerosol_microphysical_parameters",
    "OpacMix",
    "OpacMixName",
]
