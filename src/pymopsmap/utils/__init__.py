from .logging import get_logger
from .temp import get_tempdir, get_tempfile
from .types import Float64List, PosFloat64List, SortedPosFloat64List

from pathlib import Path
import os

ROOT_PATH = Path(__file__).resolve().parent.parent.parent.parent
CAMS_DATA_PATH = Path(os.getenv("CAMS_DATA_PATH", ROOT_PATH / "data/cams"))

__all__ = [
    "get_logger",
    "get_tempdir",
    "get_tempfile",
    "Float64List",
    "PosFloat64List",
    "SortedPosFloat64List",
    "CAMS_DATA_PATH",
]
