"""Concurrent runs fetch the same dataset file without colliding."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from pymopsmap.scatlib.cache import OpticalDatasetCache
from pymopsmap.scatlib.downloader import DatasetDownloader


def _source(tmp_path: Path) -> Path:
    source = tmp_path / "source" / "spheres"
    source.mkdir(parents=True)
    (source / "a.nc").write_bytes(b"x" * 4096)
    return tmp_path / "source"


class TestParallelFetch:
    def test_the_same_file_fetched_at_once_lands_intact(self, tmp_path):
        """
        A shared scratch name meant two threads wrote the same bytes over
        each other, and one of them registered a half-written file.
        """
        cache = OpticalDatasetCache(root_dir=tmp_path / "cache")
        downloader = DatasetDownloader(cache=cache, source=_source(tmp_path))

        with ThreadPoolExecutor(max_workers=8) as pool:
            list(
                pool.map(
                    lambda _: downloader.download("spheres/a.nc"), range(8)
                )
            )

        assert cache.full_path("spheres/a.nc").stat().st_size == 4096

    def test_no_scratch_file_is_left_behind(self, tmp_path):
        cache = OpticalDatasetCache(root_dir=tmp_path / "cache")
        downloader = DatasetDownloader(cache=cache, source=_source(tmp_path))

        with ThreadPoolExecutor(max_workers=4) as pool:
            list(
                pool.map(
                    lambda _: downloader.download("spheres/a.nc"), range(4)
                )
            )

        assert not list((tmp_path / "cache").rglob("*.tmp*"))
