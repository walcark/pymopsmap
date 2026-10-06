"""MOPSMAP reports some failures on stdout while still exiting zero."""

from __future__ import annotations

from pathlib import Path

import pytest

from pymopsmap.exceptions import MopsmapError


class TestErrorDetection:
    def test_a_missing_dataset_file_is_raised(self, monkeypatch, tmp_path):
        """
        The binary prints the error, produces no data, and exits 0. The exit
        code alone therefore says nothing.
        """
        from pymopsmap.engine import launcher

        stdout = (
            " Error opening /somewhere/spheres/sphere_1.5600_0.008600.nc: "
            "No such file or directory\n"
        )
        _fake_run(monkeypatch, launcher, stdout=stdout, returncode=0)

        with pytest.raises(MopsmapError, match="sphere_1.5600"):
            launcher.launch_mopsmap(tmp_path / "in.txt")

    def test_a_partial_run_is_raised_too(self, monkeypatch, tmp_path):
        """Some wavelengths can succeed while others silently do not."""
        from pymopsmap.engine import launcher

        stdout = (
            "integrated 5.32000E-01 4.8E-06 9.5E-01 6.3E-01 1.1E-01 9.6E+08"
            " 6.0E-06 9.3E-13 1.5E-06 NaN NaN NaN\n"
            " Error opening /somewhere/spheres/sphere_1.4800_0.004300.nc:"
            " No such file or directory\n"
        )
        _fake_run(monkeypatch, launcher, stdout=stdout, returncode=0)

        with pytest.raises(MopsmapError, match="sphere_1.4800"):
            launcher.launch_mopsmap(tmp_path / "in.txt")

    def test_a_fortran_stop_is_raised(self, monkeypatch, tmp_path):
        from pymopsmap.engine import launcher

        stdout = "Error: Relative humidity > 99% not supported.\n"
        _fake_run(monkeypatch, launcher, stdout=stdout, returncode=0)

        with pytest.raises(MopsmapError, match="Relative humidity"):
            launcher.launch_mopsmap(tmp_path / "in.txt")

    def test_a_clean_run_passes(self, monkeypatch, tmp_path):
        from pymopsmap.engine import launcher

        stdout = (
            "integrated 5.32000E-01 4.8E-06 9.5E-01 6.3E-01 1.1E-01 9.6E+08"
            " 6.0E-06 9.3E-13 1.5E-06 NaN NaN NaN\n"
        )
        _fake_run(monkeypatch, launcher, stdout=stdout, returncode=0)

        assert launcher.launch_mopsmap(tmp_path / "in.txt")["stdout"] == stdout

    def test_a_warning_is_not_an_error(self, monkeypatch, tmp_path):
        from pymopsmap.engine import launcher

        stdout = (
            "Warning: size parameter close to the limit\n"
            "integrated 5.32000E-01 4.8E-06 9.5E-01 6.3E-01 1.1E-01 9.6E+08"
            " 6.0E-06 9.3E-13 1.5E-06 NaN NaN NaN\n"
        )
        _fake_run(monkeypatch, launcher, stdout=stdout, returncode=0)

        assert (
            "Warning" in launcher.launch_mopsmap(tmp_path / "in.txt")["stdout"]
        )

    def test_a_non_zero_exit_still_raises(self, monkeypatch, tmp_path):
        from pymopsmap.engine import launcher

        _fake_run(monkeypatch, launcher, stdout="", returncode=2)

        with pytest.raises(MopsmapError):
            launcher.launch_mopsmap(tmp_path / "in.txt")


def _fake_run(monkeypatch, launcher, stdout: str, returncode: int) -> None:
    class Result:
        def __init__(self) -> None:
            self.stdout = stdout
            self.stderr = ""
            self.returncode = returncode

    monkeypatch.setattr(launcher.subprocess, "run", lambda *a, **k: Result())
    monkeypatch.setattr(launcher, "MOPSMAP_PATH", Path("/fake/mopsmap"))
