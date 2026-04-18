"""Input adapters — convert external data formats to pymopsmap models."""

from .cams import CamsAerosol, CamsVersion, read_aerosol_microphysical_parameters

__all__ = ["CamsAerosol", "CamsVersion", "read_aerosol_microphysical_parameters"]
