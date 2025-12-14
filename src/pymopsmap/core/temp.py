"""
temp.py — Temporary workspace management for the whole Python session.

This module:
- Creates a unique temporary directory on import.
- Exposes a decorator `use_tempdir` so functions receive that directory.
- Registers an automatic cleanup on interpreter exit.
"""

from __future__ import annotations
from threading import get_ident
from pathlib import Path
import tempfile
import shutil
import atexit
import uuid
import os


# Global temporary directory
TEMP_DIR: Path = Path(tempfile.gettempdir())
TEMP_DIR /= f"pymopsmap-temp-{uuid.uuid4().hex}"
TEMP_DIR.mkdir(parents=True, exist_ok=True)


def get_tempdir() -> Path:
    tid = get_ident()
    d = TEMP_DIR / f"worker-{tid}"
    d.mkdir(exist_ok=True)
    return d


# Cleanup method called after program termination
DISABLE_TEMP_CLEANUP = bool(os.environ.get("PYMOPSMAP_KEEP_TEMP", 0))
print(DISABLE_TEMP_CLEANUP)
print(TEMP_DIR)


def _cleanup() -> None:
    """Delete the temporary directory safely at Python exit."""
    if DISABLE_TEMP_CLEANUP:
        return
    if TEMP_DIR.exists():
        shutil.rmtree(TEMP_DIR, ignore_errors=True)


atexit.register(_cleanup)
