"""Refractive index file generation for multi-mode particle mixtures."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from pymopsmap.engine.commands import microparams_command
from pymopsmap.models import LognormalPSD, MicroParameters, Sphere

WL = [0.44, 0.55, 0.67]


def _mode(n_real: float, n_imag: float, rm: float) -> MicroParameters:
    return MicroParameters(
        wavelength=WL,
        n_real=n_real,
        n_imag=n_imag,
        shape=Sphere(),
        psd=LognormalPSD(rm=rm, sigma=1.5, n=1e6, rmin=0.001, rmax=40.0),
    )


def _refr_paths(command: str) -> list[Path]:
    return [Path(m) for m in re.findall(r"refrac file '([^']+)'", command)]


def _n_real_column(path: Path) -> list[float]:
    return [float(line.split()[1]) for line in path.read_text().splitlines()]


class TestMultiModeRefractiveIndexFiles:
    def test_each_mode_gets_its_own_file(self):
        fine, coarse = _mode(1.40, 1e-3, 0.05), _mode(1.60, 5e-2, 1.0)

        paths = _refr_paths(microparams_command([fine, coarse]))

        assert len(paths) == 2
        assert paths[0] != paths[1]

    def test_each_file_holds_its_own_refractive_index(self):
        fine, coarse = _mode(1.40, 1e-3, 0.05), _mode(1.60, 5e-2, 1.0)

        fine_path, coarse_path = _refr_paths(
            microparams_command([fine, coarse])
        )

        assert _n_real_column(fine_path) == pytest.approx([1.40] * len(WL))
        assert _n_real_column(coarse_path) == pytest.approx([1.60] * len(WL))

    def test_single_mode_writes_one_file(self):
        paths = _refr_paths(microparams_command(_mode(1.45, 1e-4, 0.1)))

        assert len(paths) == 1
        assert paths[0].exists()
