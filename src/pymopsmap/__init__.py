"""pymopsmap — Python wrapper for MOPSMAP aerosol optical properties."""

from __future__ import annotations

from dataclasses import dataclass

from pymopsmap import species
from pymopsmap.exceptions import (
    DatasetSourceNotConfiguredError,
    DomainError,
    DownloadError,
    IndexFileError,
    MopsmapError,
)
from pymopsmap.models import MicroParameters

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
from pymopsmap.species import CamsSpecie as CAMS
from pymopsmap.species import Mix, Specie, load

# ---------------------------------------------------------------------------
# Dataset cache
# ---------------------------------------------------------------------------


@dataclass
class CacheStatusReport:
    cached: list[str]
    missing: list[str]


def cache_status(modes: list[MicroParameters]) -> CacheStatusReport:
    """Return which optical dataset files are cached and which are missing."""
    from pymopsmap.scatlib.cache import OpticalDatasetCache
    from pymopsmap.scatlib.resolver import NCFileResolver

    dataset_cache = OpticalDatasetCache()
    if not dataset_cache.is_cached("index.nc"):
        return CacheStatusReport(cached=[], missing=["index.nc"])

    resolver = NCFileResolver(dataset_cache.full_path("index.nc"))
    required = resolver.resolve(modes)
    return CacheStatusReport(
        cached=[p for p in required if dataset_cache.is_cached(p)],
        missing=[p for p in required if not dataset_cache.is_cached(p)],
    )


def prefetch(modes: list[MicroParameters], quiet: bool = False) -> None:
    """Download all missing optical dataset files without running MOPSMAP."""
    from pymopsmap.scatlib.cache import OpticalDatasetCache
    from pymopsmap.scatlib.downloader import DatasetDownloader
    from pymopsmap.scatlib.resolver import NCFileResolver

    dataset_cache = OpticalDatasetCache()
    downloader = DatasetDownloader(cache=dataset_cache, quiet=quiet)

    if not dataset_cache.is_cached("index.nc"):
        downloader.download("index.nc")

    resolver = NCFileResolver(dataset_cache.full_path("index.nc"))
    downloader.download_missing(resolver.resolve(modes))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

__all__ = [
    # Species
    "CAMS",
    "Mix",
    "Specie",
    "load",
    "species",
    # Dataset cache
    "cache_status",
    "prefetch",
    # Classes
    "MicroParameters",
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
    "DomainError",
    "DownloadError",
    "IndexFileError",
    "MopsmapError",
]
