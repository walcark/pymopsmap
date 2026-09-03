"""pymopsmap — Python wrapper for MOPSMAP aerosol optical properties."""

from __future__ import annotations

from dataclasses import dataclass

from pymopsmap import species
from pymopsmap.engine.outputs import DEFAULT_OUTPUT, OutputRequest, OutputType
from pymopsmap.exceptions import (
    CoverageError,
    DatasetSourceNotConfiguredError,
    DomainError,
    DownloadError,
    IndexFileError,
    MopsmapError,
)
from pymopsmap.microparams import MicroParameters

# Re-export shapes and PSDs for convenience
from pymopsmap.psd import (
    DistrListPSD,
    DistrType,
    FileDefinedPSD,
    FixedPSD,
    LognormalPSD,
    ModifiedGammaPSD,
)
from pymopsmap.shapes import (
    Irregular,
    IrregularDistrFile,
    IrregularOverlay,
    Sphere,
    Spheroid,
    SpheroidDistrFile,
    SpheroidLognormal,
)
from pymopsmap.species import OPAC_COMPOSITIONS, Mix, Specie, load, opac_mix
from pymopsmap.species import CamsSpecie as CAMS
from pymopsmap.species import OpacSpecie as OPAC

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
    "OPAC",
    "OPAC_COMPOSITIONS",
    "opac_mix",
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
    "CoverageError",
    "DomainError",
    "DownloadError",
    "IndexFileError",
    "MopsmapError",
]
