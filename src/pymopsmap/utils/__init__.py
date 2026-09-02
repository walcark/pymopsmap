"""Utility exports — types, environment paths, logging, temp files, caching."""

import os
from pathlib import Path

from .logging import get_logger
from .temp import get_tempdir, get_tempfile
from .types import Float64List, PosFloat64List, SortedPosFloat64List
from .validation import check_within_grid

ROOT_PATH = Path(__file__).resolve().parent.parent.parent.parent
DATA_PATH = Path(os.getenv("DATA_PATH", ROOT_PATH / "data"))
MOPSMAP_PATH = ROOT_PATH / "bin/mopsmap/mopsmap"

_default_cache = Path.home() / ".cache" / "pymopsmap"
CACHE_DIR = Path(os.getenv("PYMOPSMAP_CACHE_DIR", _default_cache))
DATASET_CACHE_DIR = CACHE_DIR / "dataset"
RESULT_CACHE_DIR = CACHE_DIR / "results"
DATASET_SOURCE: str | None = os.getenv("PYMOPSMAP_DATASET_SOURCE", None)

__all__ = [
    "get_logger",
    "get_tempdir",
    "get_tempfile",
    "Float64List",
    "PosFloat64List",
    "SortedPosFloat64List",
    "check_within_grid",
    "DATA_PATH",
    "MOPSMAP_PATH",
    "CACHE_DIR",
    "DATASET_CACHE_DIR",
    "RESULT_CACHE_DIR",
    "DATASET_SOURCE",
]
