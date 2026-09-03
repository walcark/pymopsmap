"""Refractive indices reach MOPSMAP without losing significant digits."""

from __future__ import annotations

import math

import pytest

from pymopsmap.engine.commands import refr_command, write_refr_file
from pymopsmap.microparams import MicroParameters
from pymopsmap.psd import LognormalPSD
from pymopsmap.shapes import Sphere

WL = [0.44, 0.55, 0.67]
# Weakly absorbing species reach these magnitudes: CAMS sulphate goes down to
# 8.8e-9, and fixed-point formatting with six decimals writes them as zero.
NI = [1.8461538e-07, 8.83360462e-09, 3.21e-05]
NR = [1.4529, 1.4312, 1.4204]


def _columns(path) -> list[list[float]]:
    return [
        [float(token) for token in line.split()]
        for line in path.read_text().splitlines()
    ]


class TestRefractiveIndexFile:
    def test_small_imaginary_parts_survive(self):
        path = write_refr_file(wl=WL, nr=NR, ni=NI, mode_index=1)

        written = _columns(path)
        for row, expected in zip(written, NI):
            assert row[2] == pytest.approx(expected, rel=1e-9)

    def test_no_value_is_flushed_to_zero(self):
        path = write_refr_file(wl=WL, nr=NR, ni=NI, mode_index=2)

        assert all(row[2] > 0 for row in _columns(path))

    def test_real_parts_and_wavelengths_survive(self):
        path = write_refr_file(wl=WL, nr=NR, ni=NI, mode_index=3)

        written = _columns(path)
        for row, wl, nr in zip(written, WL, NR):
            assert row[0] == pytest.approx(wl, rel=1e-9)
            assert row[1] == pytest.approx(nr, rel=1e-9)

    def test_wavelengths_stay_strictly_ascending(self):
        """MOPSMAP refuses a refractive index file that is not ascending."""
        close = [0.5500000, 0.5500001, 0.5500002]
        path = write_refr_file(
            wl=close, nr=[1.45] * 3, ni=[1e-4] * 3, mode_index=4
        )

        column = [row[0] for row in _columns(path)]
        assert column == sorted(column)
        assert len(set(column)) == 3


class TestConstantRefracCommand:
    def test_single_wavelength_keeps_its_imaginary_part(self):
        command = refr_command(wl=[0.55], nr=[1.4312], ni=[8.8e-09])

        assert float(command.split()[2]) == pytest.approx(8.8e-09, rel=1e-9)


class TestModelDoesNotRound:
    def test_small_values_are_not_rounded_to_ten_decimals(self):
        mp = MicroParameters(
            wavelength=[0.55],
            n_real=[1.4312],
            n_imag=[1.8461538461538464e-07],
            shape=Sphere(),
            psd=LognormalPSD(rm=0.1, sigma=1.5, n=1e6, rmin=0.001, rmax=40.0),
        )

        assert mp.n_imag[0] == pytest.approx(1.8461538461538464e-07, rel=1e-12)

    def test_non_finite_values_are_still_rejected(self):
        from pydantic import ValidationError

        with pytest.raises(ValidationError):
            MicroParameters(
                wavelength=[0.55],
                n_real=[math.nan],
                n_imag=[1e-4],
                shape=Sphere(),
                psd=LognormalPSD(
                    rm=0.1, sigma=1.5, n=1e6, rmin=0.001, rmax=40.0
                ),
            )
