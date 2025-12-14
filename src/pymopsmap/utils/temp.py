"""
temp.py

Author  : Kévin Walcarius
Date    : 2025-11-19
Version : 1.0
License : MIT
Summary : Tools to manage a temporary workspace.
"""

from __future__ import annotations

import atexit
import os
import shutil
import tempfile
import uuid
from pathlib import Path

# --------------------------------------------------------------------------
# Global temporary directory
# --------------------------------------------------------------------------

TEMP_DIR = Path(tempfile.gettempdir()) / f"pymopsmap-{uuid.uuid4().hex}"
TEMP_DIR.mkdir(parents=True, exist_ok=True)

TEMP_REGISTRY: dict[str, Path] = {}

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------


def get_tempdir() -> Path:
    """Return the global temporary directory (ensure it exists)."""
    TEMP_DIR.mkdir(exist_ok=True)
    return TEMP_DIR


def get_tempfile(filename: str) -> Path:
    """Return a stable path inside the tempdir for a given filename."""
    if filename not in TEMP_REGISTRY:
        TEMP_REGISTRY[filename] = get_tempdir() / filename
    return TEMP_REGISTRY[filename]


# --------------------------------------------------------------------------
# Cleanup
# --------------------------------------------------------------------------

DISABLE_TEMP_CLEANUP = os.environ.get("PYMOPSMAP_KEEP_TEMP", "").lower() in (
    "1",
    "true",
    "yes",
)


def _cleanup():
    if DISABLE_TEMP_CLEANUP:
        return
    if TEMP_DIR.exists():
        shutil.rmtree(TEMP_DIR, ignore_errors=True)


atexit.register(_cleanup)
