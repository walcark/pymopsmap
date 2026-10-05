"""
Check the pipeline against an independent Mie code.

Every other validation in this directory compares MOPSMAP against numbers the
MOPSMAP authors published, so it cannot catch an error they and this wrapper
would share. ``miepython`` is a separate implementation of Mie theory, which
makes this the one check that does not go through the MOPSMAP data set at all.

Only spheres can be checked this way: a spheroid or an aggregate needs a
T-matrix code, which is what the data set exists to avoid.

Run with:
  PYMOPSMAP_DATASET_SOURCE=... pixi run -e dev pytest \
      tests/integration/test_mie_reference.py
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from pymopsmap import MicroParameters
from pymopsmap.engine import run_point
from pymopsmap.engine.outputs import DEFAULT_OUTPUT
from pymopsmap.psd import LognormalPSD
from pymopsmap.shapes import Sphere

pytestmark = pytest.mark.skipif(
    os.getenv("PYMOPSMAP_DATASET_SOURCE") is None,
    reason="PYMOPSMAP_DATASET_SOURCE not set, skipping integration tests",
)

# The mode of Table 3 of Gasteiger and Wiegner (2018), whose "explicit" column
# is itself a Mie calculation, at the two refractive indices it uses.
MODE = dict(rm=0.5, sigma=2.0, rmin=0.001, rmax=4.0, n=1e9)
WAVELENGTH_UM = 0.62832

# The data set samples the size parameter in 1 % steps and the refractive
# index on a coarse grid, which is the error Table 3 is about: the article
# reports the deviations it causes as a few tenths of a percent.
SAMPLING_TOLERANCE = 3e-3


def _mie_ensemble(
    n_real: float, n_imag: float, points: int = 20000
) -> tuple[float, float, float]:
    """
    Integrate Mie theory over the lognormal mode.

    Returns
    -------
    tuple of float
        Extinction coefficient in m-1, single scattering albedo, and
        asymmetry parameter.
    """
    import miepython

    radii = np.logspace(np.log10(MODE["rmin"]), np.log10(MODE["rmax"]), points)
    # dN/dr of a lognormal, normalised over (0, inf), as log_distr.f90 writes
    # it. The quadrature runs on ln r, so the Jacobian cancels one 1/r.
    log_sigma = np.log(MODE["sigma"])
    density = (
        MODE["n"]
        / (np.sqrt(2.0 * np.pi) * log_sigma)
        * np.exp(-((np.log(radii / MODE["rm"])) ** 2) / (2.0 * log_sigma**2))
    )

    size = 2.0 * np.pi * radii / WAVELENGTH_UM
    # miepython takes m = n - ik, the opposite sign convention.
    q_ext, q_sca, _, asymmetry = miepython.efficiencies_mx(
        complex(n_real, -n_imag), size
    )

    # Radii are in um and the concentration in m-3, so the cross sections come
    # out in um2 m-3 and need the square of a micrometre in square metres.
    area = np.pi * radii**2 * 1e-12
    weight = density * area
    log_radii = np.log(radii)

    k_ext = np.trapezoid(weight * q_ext, log_radii)
    k_sca = np.trapezoid(weight * q_sca, log_radii)
    g = np.trapezoid(weight * q_sca * asymmetry, log_radii) / k_sca
    return float(k_ext), float(k_sca / k_ext), float(g)


@pytest.mark.parametrize(
    ("n_real", "n_imag"), [(1.52, 0.0), (1.54, 0.005), (1.45, 0.001)]
)
def test_pipeline_matches_an_independent_mie_code(
    n_real: float, n_imag: float
) -> None:
    """One lognormal mode of spheres, computed twice by different means."""
    mode = MicroParameters(
        wavelength=[WAVELENGTH_UM],
        n_real=n_real,
        n_imag=n_imag,
        shape=Sphere(),
        psd=LognormalPSD(**MODE),
    )
    result = run_point([mode], DEFAULT_OUTPUT, quiet=True)

    k_ext, ssa, g = _mie_ensemble(n_real, n_imag)
    assert float(result["kext"].values.ravel()[0]) == pytest.approx(
        k_ext, rel=SAMPLING_TOLERANCE
    )
    assert float(result["ssa"].values.ravel()[0]) == pytest.approx(
        ssa, rel=SAMPLING_TOLERANCE
    )
    assert float(result["g"].values.ravel()[0]) == pytest.approx(
        g, rel=SAMPLING_TOLERANCE
    )
