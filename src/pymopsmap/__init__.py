"""
Aerosol optical properties from microphysics, on the grid you ask for.

A wrapper around MOPSMAP (Gasteiger and Wiegner, 2018, GMD 11:2739). Load an
aerosol, compute its optical properties:

    import pymopsmap as pm

    op = pm.load(pm.CAMS.SULPHATE).compute(wl=[0.44, 0.55], rh=[0, 50, 90])

The result is a plain xarray Dataset, with a ``.mopsmap`` accessor for the
conversions this library adds.
"""

from __future__ import annotations

from pymopsmap import accessors as _accessors  # noqa: F401  registers .mopsmap
from pymopsmap import psd, shapes
from pymopsmap.engine.outputs import OutputType
from pymopsmap.exceptions import (
    CoverageError,
    DomainError,
    MopsmapError,
    SchemaError,
)
from pymopsmap.microparams import MicroParameters
from pymopsmap.species import CamsSpecie as CAMS
from pymopsmap.species import Mix, Specie, load, opac_mix
from pymopsmap.species import OpacSpecie as OPAC

__all__ = [
    # Loading a species, and mixing several
    "load",
    "Specie",
    "Mix",
    "CAMS",
    "OPAC",
    "opac_mix",
    # Building one by hand
    "shapes",
    "psd",
    "MicroParameters",
    # Asking for outputs
    "OutputType",
    # Errors a caller acts on
    "CoverageError",
    "DomainError",
    "MopsmapError",
    "SchemaError",
]
