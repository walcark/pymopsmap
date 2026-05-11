"""pymopsmap — Python wrapper for MOPSMAP aerosol optical properties."""

from __future__ import annotations

from dataclasses import dataclass

import xarray as xr

from pymopsmap.engine import compute_optical_properties
from pymopsmap.exceptions import (
    DatasetSourceNotConfiguredError,
    DownloadError,
    IndexFileError,
    MopsmapError,
)
from pymopsmap.models import (
    MicroParameters,
    OptiProps,
    ParametricSweep,
    ParticleMixture,
)

# Re-export shapes and PSDs for convenience
from pymopsmap.models.microparams import (
    DistrListPSD,
    DistrType,
    FileDefinedPSD,
    FixedPSD,
    Irregular,
    IrregularDistrFile,
    IrregularOverlay,
    LognormalPSD,
    ModifiedGammaPSD,
    Sphere,
    Spheroid,
    SpheroidDistrFile,
    SpheroidLognormal,
)
from pymopsmap.models.output_request import (
    DEFAULT_OUTPUT,
    OutputRequest,
    OutputType,
)

# ---------------------------------------------------------------------------
# Cache status
# ---------------------------------------------------------------------------


@dataclass
class CacheStatusReport:
    cached: list[str]
    missing: list[str]


def cache_status(
    mps: MicroParameters | ParticleMixture | ParametricSweep,
) -> CacheStatusReport:
    """Return which dataset files are cached and which are missing."""
    from pymopsmap.cache.optical import OpticalDatasetCache
    from pymopsmap.cache.resolver import NCFileResolver

    dataset_cache = OpticalDatasetCache()
    index_path = dataset_cache.full_path("index.nc")

    if not dataset_cache.is_cached("index.nc"):
        return CacheStatusReport(cached=[], missing=["index.nc"])

    resolver = NCFileResolver(index_path)
    modes = _flatten_modes(mps)
    required = resolver.resolve(modes)

    cached = [p for p in required if dataset_cache.is_cached(p)]
    missing = [p for p in required if not dataset_cache.is_cached(p)]
    return CacheStatusReport(cached=cached, missing=missing)


# ---------------------------------------------------------------------------
# Prefetch
# ---------------------------------------------------------------------------


def prefetch(
    mps: MicroParameters | ParticleMixture | ParametricSweep,
    quiet: bool = False,
) -> None:
    """Download all missing dataset files without running MOPSMAP."""
    from pymopsmap.cache.downloader import DatasetDownloader
    from pymopsmap.cache.optical import OpticalDatasetCache
    from pymopsmap.cache.resolver import NCFileResolver

    dataset_cache = OpticalDatasetCache()
    downloader = DatasetDownloader(cache=dataset_cache, quiet=quiet)

    if not dataset_cache.is_cached("index.nc"):
        downloader.download("index.nc")

    index_path = dataset_cache.full_path("index.nc")
    resolver = NCFileResolver(index_path)
    modes = _flatten_modes(mps)
    required = resolver.resolve(modes)
    downloader.download_missing(required)


# ---------------------------------------------------------------------------
# Primary compute entry point
# ---------------------------------------------------------------------------


def compute(
    mps: MicroParameters | ParticleMixture | ParametricSweep,
    output_types: OutputRequest = DEFAULT_OUTPUT,
    rh: float | None = None,
    quiet: bool = False,
) -> OptiProps:
    """
    Compute optical properties.

    Handles file resolution, caching, downloading, and MOPSMAP execution
    transparently.
    """
    if isinstance(mps, MicroParameters):
        mps = ParticleMixture([mps])
    return compute_optical_properties(
        target=mps,
        output_types=output_types,
        rh=rh,
        quiet=quiet,
    )


# ---------------------------------------------------------------------------
# Convenience wrappers
# ---------------------------------------------------------------------------


def kext(
    mps: MicroParameters | ParticleMixture | ParametricSweep,
    rh: float | None = None,
    quiet: bool = False,
) -> xr.DataArray:
    return compute(mps, rh=rh, quiet=quiet).sel("kext")


def ssa(
    mps: MicroParameters | ParticleMixture | ParametricSweep,
    rh: float | None = None,
    quiet: bool = False,
) -> xr.DataArray:
    return compute(mps, rh=rh, quiet=quiet).sel("ssa")


def phase(
    mps: MicroParameters | ParticleMixture | ParametricSweep,
    rh: float | None = None,
    quiet: bool = False,
) -> xr.DataArray:
    from pymopsmap.models.output_request import OutputType

    return compute(
        mps,
        output_types=frozenset(
            {OutputType.INTEGRATED, OutputType.PHASE_FUNCTION}
        ),
        rh=rh,
        quiet=quiet,
    ).sel("phase")


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _flatten_modes(
    mps: MicroParameters | ParticleMixture | ParametricSweep,
) -> list[MicroParameters]:
    if isinstance(mps, ParametricSweep):
        result = []
        for mixture in mps.mixtures:
            result.extend(mixture.modes)
        return result
    if isinstance(mps, ParticleMixture):
        return list(mps.modes)
    return [mps]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

__all__ = [
    # Main entry points
    "compute",
    "cache_status",
    "prefetch",
    "kext",
    "ssa",
    "phase",
    # Classes
    "MicroParameters",
    "ParticleMixture",
    "ParametricSweep",
    "OptiProps",
    "CacheStatusReport",
    # Output types
    "OutputType",
    "OutputRequest",
    "DEFAULT_OUTPUT",
    # Shape types
    "Sphere",
    "Spheroid",
    "SpheroidLognormal",
    "SpheroidDistrFile",
    "Irregular",
    "IrregularDistrFile",
    "IrregularOverlay",
    # PSD types
    "FixedPSD",
    "LognormalPSD",
    "ModifiedGammaPSD",
    "FileDefinedPSD",
    "DistrListPSD",
    "DistrType",
    # Exceptions
    "DatasetSourceNotConfiguredError",
    "DownloadError",
    "IndexFileError",
    "MopsmapError",
]
