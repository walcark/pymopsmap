"""Inspecting and pre-fetching the dataset files a species needs."""

from __future__ import annotations

import pytest

import pymopsmap as pm

WL = [0.44, 0.55]


@pytest.fixture
def resolved(monkeypatch):
    """Record what the species asks the dataset for, without touching disk."""
    asked: dict = {}

    class FakeCache:
        def __init__(self, *a, **k) -> None:
            pass

        def is_cached(self, path: str) -> bool:
            return path == "index.nc"

        def full_path(self, path: str):
            return path

    class FakeDownloader:
        def __init__(self, *a, **k) -> None:
            pass

        def download(self, path: str) -> None:
            asked.setdefault("downloaded", []).append(path)

        def download_missing(self, paths: list[str]) -> None:
            asked.setdefault("downloaded", []).extend(paths)

    class FakeResolver:
        def __init__(self, *a, **k) -> None:
            pass

        def resolve(self, modes, rh=None):
            asked["modes"] = modes
            asked["rh"] = rh
            return ["index.nc", "spheres/a.nc", "spheres/b.nc"]

    monkeypatch.setattr(
        "pymopsmap.scatlib.cache.OpticalDatasetCache", FakeCache
    )
    monkeypatch.setattr(
        "pymopsmap.scatlib.downloader.DatasetDownloader", FakeDownloader
    )
    monkeypatch.setattr(
        "pymopsmap.scatlib.resolver.NCFileResolver", FakeResolver
    )
    return asked


class TestCacheStatus:
    def test_it_splits_cached_from_missing(self, resolved):
        report = pm.load(pm.CAMS.SULPHATE).cache_status(wl=WL, rh=50.0)

        assert report.cached == ["index.nc"]
        assert report.missing == ["spheres/a.nc", "spheres/b.nc"]

    def test_it_resolves_for_both_modes(self, resolved):
        pm.load(pm.CAMS.SULPHATE).cache_status(wl=WL, rh=50.0)

        assert len(resolved["modes"]) == 2

    def test_it_passes_the_humidity_the_engine_will_apply(self, resolved):
        """A kappa species is grown by MOPSMAP, which moves its files."""
        pm.load(pm.OPAC.WASO).cache_status(wl=WL, rh=80.0)

        assert resolved["rh"] == 80.0

    def test_a_tabulated_species_leaves_no_growth_to_the_engine(
        self, resolved
    ):
        pm.load(pm.CAMS.SULPHATE).cache_status(wl=WL, rh=50.0)

        assert resolved["rh"] is None


class TestPrefetch:
    def test_it_downloads_what_is_missing(self, resolved):
        pm.load(pm.CAMS.SULPHATE).prefetch(wl=WL, rh=50.0)

        assert "spheres/a.nc" in resolved["downloaded"]

    def test_it_returns_nothing(self, resolved):
        assert pm.load(pm.CAMS.SULPHATE).prefetch(wl=WL, rh=50.0) is None
