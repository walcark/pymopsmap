from .logging import get_logger
from .temp import get_tempdir, get_tempfile
from .types import Float64List, PosFloat64List, SortedPosFloat64List

from pathlib import Path
import os

ROOT_PATH = Path(__file__).resolve().parent.parent.parent.parent
DATA_PATH = Path(os.getenv("DATA_PATH", ROOT_PATH / "data"))
MOPSMAP_PATH = ROOT_PATH / "bin/mopsmap/mopsmap"

__all__ = [
    "get_logger",
    "get_tempdir",
    "get_tempfile",
    "Float64List",
    "PosFloat64List",
    "SortedPosFloat64List",
    "DATA_PATH",
    "MOPSMAP_PATH",
]
