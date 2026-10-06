"""Every output type a request can ask for, through the public path."""

from __future__ import annotations

import os

import pytest

import pymopsmap as pm

pytestmark = pytest.mark.skipif(
    os.getenv("PYMOPSMAP_DATASET_SOURCE") is None,
    reason="PYMOPSMAP_DATASET_SOURCE not set, skipping integration tests",
)

ANGLES = 181


@pytest.fixture(scope="module")
def specie():
    return pm.load(pm.CAMS.SULPHATE)


@pytest.mark.parametrize(
    ("output", "expected"),
    [
        (pm.OutputType.INTEGRATED, {"wl"}),
        (pm.OutputType.LIDAR, {"wl"}),
        (pm.OutputType.PHASE_FUNCTION, {"wl", "theta"}),
        (pm.OutputType.VOLUME_SCATTERING_FUNCTION, {"wl", "theta"}),
        (pm.OutputType.SCATTERING_MATRIX, {"wl", "theta", "element"}),
        (pm.OutputType.COEFF, {"wl", "l", "coeff_element"}),
    ],
)
def test_compute_returns_every_output_type(specie, output, expected):
    """
    A sweep has to declare the axes of what it collects.

    Every angular output used to fail on its first point with "cannot reshape
    array of size 2000 into shape (1,)", because the contract declared one
    number per wavelength for all of them.
    """
    result = specie.compute(
        wl=[0.55],
        rh=50,
        outputs=frozenset({output}),
        n_angles=ANGLES,
        quiet=True,
    )

    assert expected <= set(result.sizes)
    if "theta" in expected:
        assert result.sizes["theta"] == ANGLES


def test_the_angle_count_is_part_of_the_store_key(specie):
    """Two angle counts are two shapes, so they cannot share a store."""
    coarse = specie.compute(
        wl=[0.55],
        rh=50,
        outputs=frozenset({pm.OutputType.PHASE_FUNCTION}),
        n_angles=91,
        quiet=True,
    )
    fine = specie.compute(
        wl=[0.55],
        rh=50,
        outputs=frozenset({pm.OutputType.PHASE_FUNCTION}),
        n_angles=361,
        quiet=True,
    )

    assert coarse.sizes["theta"] == 91
    assert fine.sizes["theta"] == 361
