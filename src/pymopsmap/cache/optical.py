"""Optical dataset local cache — lookup and storage of NC files by key."""

import os
from pathlib import Path

from pymopsmap.utils import DATASET_CACHE_DIR


class OpticalDatasetCache:
    def __init__(self, root_dir: Path | None = None):
        self.root_dir = root_dir or DATASET_CACHE_DIR
        self.root_dir.mkdir(parents=True, exist_ok=True)

    def is_cached(self, relative_path: str) -> bool:
        p = self.full_path(relative_path)
        return p.exists() and p.stat().st_size > 0

    def full_path(self, relative_path: str) -> Path:
        return self.root_dir / relative_path

    def register(self, relative_path: str, tmp_path: Path) -> None:
        target = self.full_path(relative_path)
        target.parent.mkdir(parents=True, exist_ok=True)
        os.replace(tmp_path, target)

    def cached_files(self) -> list[str]:
        result = []
        for p in self.root_dir.rglob("*.nc"):
            result.append(str(p.relative_to(self.root_dir)))
        return sorted(result)
