"""
Validation against the published numbers of Gasteiger and Wiegner (2018).

Reference: J. Gasteiger and M. Wiegner, "MOPSMAP v1.0: a versatile tool for
the modeling of aerosol optical properties", Geosci. Model Dev. 11:2739-2762,
doi:10.5194/gmd-11-2739-2018.

Every expected value here is read from the article or from the reference
output files the authors ship in ``bin/mopsmap/misc/paper_examples``, not from
a previous run of this package. A figure can be eyeballed; a table cannot, so
the tables are what the pipeline is held to.

Run with:
  PYMOPSMAP_DATASET_SOURCE=... pixi run -e dev pytest \
      tests/integration/test_gasteiger_2018.py
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest
import xarray as xr

from pymopsmap import MicroParameters
from pymopsmap.engine import run_point
from pymopsmap.engine.outputs import OutputType
from pymopsmap.psd import DistrListPSD, DistrType, LognormalPSD
from pymopsmap.shapes import (
    Irregular,
    Sphere,
    Spheroid,
    SpheroidDistrFile,
)

pytestmark = pytest.mark.skipif(
    os.getenv("PYMOPSMAP_DATASET_SOURCE") is None,
    reason="PYMOPSMAP_DATASET_SOURCE not set, skipping integration tests",
)

INTEGRATED_AND_LIDAR = frozenset({OutputType.INTEGRATED, OutputType.LIDAR})

# Table 3 and Table 4 are printed to four significant digits, so a relative
# agreement of 1e-3 is the most the published values can support. The lidar
# ratio and the depolarisation ratio are quotients of two interpolated
# quantities and lose a digit.
TIGHT = 1.5e-3
LOOSE = 1.0e-2


def _value(result, name: str) -> float:
    """The single wavelength point of one output variable."""
    return float(result[name].values.ravel()[0])


# --------------------------------------------------------------------------
# Table 3: one lognormal mode, size sampling and refractive index
# interpolation. Section 3.3 of the article.
# --------------------------------------------------------------------------

# The article gives no concentration; 1000 cm-3 is what reproduces its
# extinction coefficients, and the other four rows do not depend on it.
TABLE3_PSD = dict(rm=0.5, sigma=2.0, rmin=0.001, rmax=4.0, n=1e9)
TABLE3_WL = 0.62832

# Columns "data set" of Table 3: alpha_ext (km-1), omega_0, g, S (sr), delta_l.
TABLE3 = {
    ("sphere", 1.52, 0.0): (4.808, 1.0000, 0.7045, 10.52, 0.0000),
    ("spheroid", 1.52, 0.0): (4.863, 1.0000, 0.7018, 42.75, 0.3063),
    ("sphere", 1.54, 0.005): (4.793, 0.8845, 0.7331, 13.13, 0.0000),
    ("spheroid", 1.54, 0.005): (4.844, 0.8892, 0.7382, 58.25, 0.2502),
}


@pytest.mark.parametrize(("shape_name", "n_real", "n_imag"), list(TABLE3))
def test_table3(shape_name: str, n_real: float, n_imag: float) -> None:
    """Reproduce one column of Table 3."""
    shape = (
        Sphere()
        if shape_name == "sphere"
        else Spheroid(mode="prolate", aspect_ratio=2.0)
    )
    mode = MicroParameters(
        wavelength=[TABLE3_WL],
        n_real=n_real,
        n_imag=n_imag,
        shape=shape,
        psd=LognormalPSD(**TABLE3_PSD),
    )
    result = run_point([mode], INTEGRATED_AND_LIDAR, quiet=True)

    kext, ssa, g, lidar_ratio, depol = TABLE3[(shape_name, n_real, n_imag)]
    # MOPSMAP reports the extinction coefficient in m-1, the article in km-1.
    assert _value(result, "kext") * 1e3 == pytest.approx(kext, rel=TIGHT)
    assert _value(result, "ssa") == pytest.approx(ssa, rel=TIGHT)
    assert _value(result, "g") == pytest.approx(g, rel=TIGHT)
    assert _value(result, "lidar_ratio") == pytest.approx(
        lidar_ratio, rel=LOOSE
    )
    assert _value(result, "depol_ratio") == pytest.approx(depol, abs=2e-3)


# --------------------------------------------------------------------------
# Table 4: the five COSMO-MUSCAT dust size bins at 500 nm. Section 5.2.
#
# "The size bins are determined by the radius limits 0.1, 0.3, 0.9, 2.6, 8,
# and 24 um. We assumed constant dv/dlnr within each bin." A two-point
# dv/dlnr table with equal values is exactly that.
# --------------------------------------------------------------------------

BIN_EDGES = [0.1, 0.3, 0.9, 2.6, 8.0, 24.0]
DUST_WL = 0.5
DUST_INDEX = (1.53, 0.0078)

# The aspect ratio distribution of Kandler et al. (2009), as MOPSMAP ships it.
AR_KANDLER = str(Path("bin/mopsmap/data/ar_kandler").resolve())

# Table 4: omega_0 and g of each bin, spheres then prolate spheroids. "For the
# latter case we assumed volume-equivalent sizes to keep the particle mass
# constant", which is size_equ vol.
TABLE4 = {
    (1, "spheres"): (0.9632, 0.6567),
    (2, "spheres"): (0.9216, 0.6866),
    (3, "spheres"): (0.7903, 0.8088),
    (4, "spheres"): (0.6450, 0.8998),
    (5, "spheres"): (0.5561, 0.9442),
    (1, "spheroids"): (0.9628, 0.6585),
    (2, "spheroids"): (0.9264, 0.7111),
    (3, "spheroids"): (0.7934, 0.8109),
    (4, "spheroids"): (0.6485, 0.9017),
    (5, "spheroids"): (0.5601, 0.9419),
}


@pytest.mark.parametrize(("bin_index", "shape_name"), list(TABLE4))
def test_table4(bin_index: int, shape_name: str) -> None:
    """Reproduce one half of one Table 4 column."""
    low, high = BIN_EDGES[bin_index - 1], BIN_EDGES[bin_index]
    spherical = shape_name == "spheres"
    mode = MicroParameters(
        wavelength=[DUST_WL],
        n_real=DUST_INDEX[0],
        n_imag=DUST_INDEX[1],
        shape=(
            Sphere()
            if spherical
            else SpheroidDistrFile(distr_filename=AR_KANDLER)
        ),
        psd=DistrListPSD(
            radii=[low, high],
            concentrations=[1.0, 1.0],
            distr_type=DistrType.DVDLNR,
        ),
        size_equ="cs" if spherical else "vol",
    )
    result = run_point([mode], INTEGRATED_AND_LIDAR, quiet=True)

    ssa, g = TABLE4[(bin_index, shape_name)]
    assert _value(result, "ssa") == pytest.approx(ssa, rel=TIGHT)
    assert _value(result, "g") == pytest.approx(g, rel=TIGHT)


# --------------------------------------------------------------------------
# Table 6: the reference ensemble of the uncertainty example. Section 5.5,
# whose text gives the three reference values outright.
# --------------------------------------------------------------------------


def test_table6_reference_ensemble() -> None:
    """Reproduce the dust-like ensemble the Jacobian is taken around."""
    mode = MicroParameters(
        wavelength=[0.532],
        n_real=1.53,
        n_imag=0.0063,
        shape=Spheroid(mode="prolate", aspect_ratio=2.0),
        psd=LognormalPSD(rm=0.1, sigma=2.6, rmin=0.001, rmax=20.0, n=1e6),
    )
    result = run_point([mode], INTEGRATED_AND_LIDAR, quiet=True)

    assert _value(result, "ssa") == pytest.approx(0.9020, rel=TIGHT)
    assert _value(result, "g") == pytest.approx(0.7319, rel=TIGHT)
    assert _value(result, "lidar_ratio") == pytest.approx(69.95, rel=LOOSE)


# --------------------------------------------------------------------------
# Table 5: one lognormal mode seen through the three size equivalences.
# Section 5.4. The expected values come from the reference runs the authors
# ship in ``paper_examples/sect_54_size_equivalence``, which carry six digits
# where the article prints three.
#
# Those runs were made with the concentration in cm-3 where MOPSMAP reads m-3,
# so every extensive quantity in them is 1e6 below the article's. The factor
# is applied here rather than in the expectations, which stay as shipped.
# --------------------------------------------------------------------------

# "N0 = 10^3.66 cm-3, which results in a concentration of N = 100 cm-3 in the
# range from rmin to rmax": the two statements disagree, since the lognormal
# puts 96.5 % of its particles inside that range. The in-range concentration is
# the one that reproduces the table, and MOPSMAP echoes it back as n.
TABLE5_N = 100.0 / 0.9646465211511944 * 1e6
TABLE5_SCALE = 1e6

TABLE5_PSD = dict(rm=0.5, sigma=2.0, rmin=0.001, rmax=1.75, n=TABLE5_N)

# Columns of output_sphere, output_cs, output_vol and output_vol_cs: the
# integrated line (ext, ssa, g, reff, n, a, v, mass) and the lidar line
# (backscatter, S, delta_l). Shape D is the aggregate of Gasteiger et al.
# (2011b), of which the article uses xi_vc = 0.8708.
TABLE5: dict[str, tuple[str, dict[str, float]]] = {
    "sphere": (
        "cs",
        dict(
            kext=3.50095e-10,
            ssa=8.97033e-01,
            g=7.21522e-01,
            reff=9.84353e-01,
            cross_dens=1.41136e-10,
            vol_dens=1.85237e-16,
            mass_conc=4.81616e-10,
            backscatter=3.02728e-11,
            lidar_ratio=1.15647e01,
            depol_ratio=0.0,
        ),
    ),
    "aggregate-cs": (
        "cs",
        dict(
            kext=3.47313e-10,
            ssa=9.22396e-01,
            g=6.79080e-01,
            reff=9.83807e-01,
            cross_dens=1.41051e-10,
            vol_dens=1.22174e-16,
            mass_conc=3.17653e-10,
            backscatter=1.03345e-11,
            lidar_ratio=3.36072e01,
            depol_ratio=4.49701e-01,
        ),
    ),
    "aggregate-vol": (
        "vol",
        dict(
            kext=4.49223e-10,
            ssa=9.10603e-01,
            g=6.80427e-01,
            reff=9.83785e-01,
            cross_dens=1.86012e-10,
            vol_dens=1.85020e-16,
            mass_conc=4.81051e-10,
            backscatter=1.36772e-11,
            lidar_ratio=3.28447e01,
            depol_ratio=4.54431e-01,
        ),
    ),
    "aggregate-vol-cs": (
        "vol_cs_ratio",
        dict(
            kext=7.50807e-10,
            ssa=8.83109e-01,
            g=6.89094e-01,
            reff=9.83764e-01,
            cross_dens=3.23499e-10,
            vol_dens=4.24329e-16,
            mass_conc=1.10325e-09,
            backscatter=2.27816e-11,
            lidar_ratio=3.29566e01,
            depol_ratio=4.53986e-01,
        ),
    ),
}


@pytest.mark.parametrize("case", sorted(TABLE5))
def test_table5(case: str) -> None:
    """Reproduce one column of Table 5."""
    size_equ, expected = TABLE5[case]
    shape = Sphere() if case == "sphere" else Irregular(shape_id="D")
    mode = MicroParameters(
        wavelength=[0.532],
        n_real=1.54,
        n_imag=0.005,
        shape=shape,
        psd=LognormalPSD(**TABLE5_PSD),
        density=2.6,
        size_equ=size_equ,
    )
    result = run_point([mode], INTEGRATED_AND_LIDAR, quiet=True)

    for name, reference in expected.items():
        obtained = _value(result, name)
        if name in _INTENSIVE:
            assert obtained == pytest.approx(reference, rel=TIGHT, abs=2e-3), (
                name
            )
        else:
            assert obtained == pytest.approx(
                reference * TABLE5_SCALE, rel=TIGHT
            ), name


# The quantities that do not scale with the particle concentration.
_INTENSIVE = frozenset({"ssa", "g", "reff", "lidar_ratio", "depol_ratio"})


# --------------------------------------------------------------------------
# Table 6: the Jacobian around that ensemble, by the same central difference
# the authors use in paper_examples/sect_55_error_calc/jacobian.py, a relative
# perturbation of one percent.
#
# The article warns that "partial derivatives d zeta / d mr are constant
# between the mr grid points of the data set", so these are interpolation
# slopes rather than physical ones. They are also differences of nearly equal
# numbers printed to two or three significant digits: d omega_0 / d eps comes
# out of a change in omega_0 of 2e-4. Ten percent is what that supports.
# --------------------------------------------------------------------------

TABLE6 = {
    "n_real": {"ssa": -0.037, "g": -0.428, "lidar_ratio": -360.0},
    "n_imag": {"ssa": -11.0, "g": +3.69, "lidar_ratio": +2839.0},
    "aspect_ratio": {"ssa": +0.010, "g": +0.058, "lidar_ratio": +48.3},
}

TABLE6_REFERENCE = {"n_real": 1.53, "n_imag": 0.0063, "aspect_ratio": 2.0}


def _table6_run(**overrides: float) -> xr.Dataset:
    """The reference ensemble of section 5.5, with one parameter moved."""
    values = TABLE6_REFERENCE | overrides
    mode = MicroParameters(
        wavelength=[0.532],
        n_real=values["n_real"],
        n_imag=values["n_imag"],
        shape=Spheroid(mode="prolate", aspect_ratio=values["aspect_ratio"]),
        psd=LognormalPSD(rm=0.1, sigma=2.6, rmin=0.001, rmax=20.0, n=1e6),
    )
    return run_point([mode], INTEGRATED_AND_LIDAR, quiet=True)


@pytest.mark.parametrize("parameter", sorted(TABLE6))
def test_table6_jacobian(parameter: str) -> None:
    """Reproduce one row of the Jacobian matrix."""
    reference = TABLE6_REFERENCE[parameter]
    low = _table6_run(**{parameter: reference * 0.99})
    high = _table6_run(**{parameter: reference * 1.01})
    span = reference * 0.02

    for name, expected in TABLE6[parameter].items():
        derivative = (_value(high, name) - _value(low, name)) / span
        assert derivative == pytest.approx(expected, rel=0.1), name
