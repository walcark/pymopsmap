"""Non-finite numeric inputs are rejected at validation time."""

from __future__ import annotations

import math

import pytest
from pydantic import ValidationError

from pymopsmap.microparams import MicroParameters
from pymopsmap.psd import LognormalPSD
from pymopsmap.shapes import Sphere


def _micro(n_real, n_imag) -> MicroParameters:
    return MicroParameters(
        wavelength=[0.55],
        n_real=n_real,
        n_imag=n_imag,
        shape=Sphere(),
        psd=LognormalPSD(rm=0.1, sigma=1.5, n=1e6, rmin=0.001, rmax=40.0),
    )


class TestNonFiniteRefractiveIndex:
    def test_nan_real_part_is_rejected(self):
        with pytest.raises(ValidationError):
            _micro([math.nan], [1e-4])

    def test_nan_imaginary_part_is_rejected(self):
        with pytest.raises(ValidationError):
            _micro([1.45], [math.nan])

    def test_infinite_value_is_rejected(self):
        with pytest.raises(ValidationError):
            _micro([math.inf], [1e-4])

    def test_finite_values_pass(self):
        assert _micro([1.45], [1e-4]).n_real == [1.45]


class TestNonFiniteWavelength:
    def test_nan_wavelength_is_rejected(self):
        with pytest.raises(ValidationError):
            MicroParameters(
                wavelength=[0.44, math.nan],
                n_real=1.45,
                n_imag=1e-4,
                shape=Sphere(),
                psd=LognormalPSD(
                    rm=0.1, sigma=1.5, n=1e6, rmin=0.001, rmax=40.0
                ),
            )
