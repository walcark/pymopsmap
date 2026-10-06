"""Each MOPSMAP run owns its directory, so concurrent runs cannot collide."""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor

import pytest

from pymopsmap.engine.workspace import Workspace


class TestIsolation:
    def test_two_workspaces_have_distinct_directories(self):
        with Workspace() as first, Workspace() as second:
            assert first.path != second.path

    def test_the_same_file_name_lands_in_each_of_them(self):
        with Workspace() as first, Workspace() as second:
            assert first.file("mopsmap.txt") != second.file("mopsmap.txt")

    def test_concurrent_workspaces_do_not_share_a_path(self):
        def make(_):
            with Workspace() as workspace:
                path = workspace.file("mopsmap.txt")
                path.write_text(str(workspace.path))
                return path.read_text()

        with ThreadPoolExecutor(max_workers=8) as pool:
            seen = list(pool.map(make, range(16)))

        assert len(set(seen)) == 16


class TestLifetime:
    def test_the_directory_exists_inside_the_block(self):
        with Workspace() as workspace:
            assert workspace.path.is_dir()

    def test_it_is_removed_on_exit(self):
        with Workspace() as workspace:
            path = workspace.path

        assert not path.exists()

    def test_a_file_is_reused_within_one_run(self):
        with Workspace() as workspace:
            assert workspace.file("ri_1.txt") == workspace.file("ri_1.txt")

    def test_keeping_it_is_opt_in(self, monkeypatch):
        monkeypatch.setenv("PYMOPSMAP_KEEP_TEMP", "1")

        with Workspace() as workspace:
            path = workspace.path

        assert path.exists()
        os.rmdir(path) if not any(path.iterdir()) else None

    def test_it_is_removed_even_when_the_run_raises(self):
        with pytest.raises(RuntimeError):
            with Workspace() as workspace:
                path = workspace.path
                raise RuntimeError("boom")

        assert not path.exists()
