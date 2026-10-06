"""
Reproduce Figure 9 of Gasteiger and Wiegner (2018), GMD 11:2739.

    PYMOPSMAP_DATASET_SOURCE=... \
        pixi run -e dev python -m scripts.validation.gasteiger_fig9

"Volume scattering function of dust at lambda = 355 nm (arbitrary scale) using
either the mi distribution measured by Kandler et al. (2011), the average mi of
these measurements, or applying the non-absorbing fraction parameterization
with different X."

Four of the five curves. The red one needs the size-resolved imaginary index
distribution of Kandler et al. (2011), a supplement that MOPSMAP does not ship
with its own version of this example, so it is not here. What is here is the
black curve, every particle carrying the average index, and the three blue
ones, where a fraction X of the particles does not absorb at all and the rest
absorbs the more for it.

The ensemble is the OPAC desert type at RH = 0 %, with the mineral components
as prolate spheroids with the aspect ratio distribution of Kandler et al.
(2009), as the authors' own script
``paper_examples/sect_56_dust_refr_variability/dust_refr_variability.py``
sets it up.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

from pymopsmap import MicroParameters
from pymopsmap.engine import run_point
from pymopsmap.engine.outputs import OutputType
from pymopsmap.psd import LognormalPSD
from pymopsmap.shapes import Sphere, SpheroidDistrFile

FIGURES = Path(__file__).resolve().parents[2] / "docs" / "figures"
AR_KANDLER = "ar_kandler"

WAVELENGTH_UM = 0.355
N_ANGLES = 1801

# The average imaginary index of the Kandler et al. (2011) measurements, which
# section 5.6 gives as 0.0175, "close to 0.0166 given for the mineral
# components in OPAC at 355 nm".
MINERAL_N_IMAG = 0.0175
WASO_INDEX = (1.53, 0.005)
N_REAL = 1.53

# n [cm-3], r_mod [um], sigma, spheroidal.
MODES = [
    (6.667, 0.0212, 2.24, False),
    (0.898, 0.07, 1.95, True),
    (0.1016, 0.39, 2.00, True),
    (0.000472, 1.9, 2.15, True),
]
# The authors read the radius range off the Kandler bins. Without that file
# the OPAC range is used instead, which is close enough that the integrated
# properties still land on the published ones.
R_MIN_UM = 0.005
R_MAX_UM = 20.0

# The four cases the figure draws, and the colours the caption gives them.
CASES = {
    "Average $m_i$": (0.0, "black", "-", 2.0),
    "Average $m_i$ with nonabs. fraction $X$=0.25": (
        0.25,
        "blue",
        "--",
        0.9,
    ),
    "Average $m_i$ with nonabs. fraction $X$=0.50": (0.5, "blue", "-", 2.0),
    "Average $m_i$ with nonabs. fraction $X$=0.75": (0.75, "blue", "-", 0.9),
}

OUTPUTS = frozenset(
    {OutputType.INTEGRATED, OutputType.LIDAR, OutputType.PHASE_FUNCTION}
)

# What section 5.6 states for the two cases that do not need the measured
# index distribution: "omega_0 = 0.741 when using the average mi and
# omega_0 = 0.834 using the parameterization with X = 0.5", then "for the
# asymmetry parameter g, we obtain 0.744, 0.789, and 0.749 for the measured,
# averaged, and parameterized cases", "values of 41, 78, and 42 sr" for the
# lidar ratio and "0.241, 0.212, and 0.220" for the depolarisation.
PUBLISHED = {
    "Average $m_i$": "   0.741  0.789    78.00   0.2120",
    "Average $m_i$ with nonabs. fraction $X$=0.50": (
        "   0.834  0.749    42.00   0.2200"
    ),
}

# What the report prints in the left column, since the legend labels are the
# ones the published figure uses and do not fit a table.
SHORT = {
    "Average $m_i$": "average m_i",
    "Average $m_i$ with nonabs. fraction $X$=0.25": "X = 0.25",
    "Average $m_i$ with nonabs. fraction $X$=0.50": "X = 0.50",
    "Average $m_i$ with nonabs. fraction $X$=0.75": "X = 0.75",
}


def _modes(nonabs_fraction: float) -> list[MicroParameters]:
    """The four desert modes, with one non-absorbing fraction on the dust."""
    built = []
    for concentration, rm, sigma, spheroidal in MODES:
        built.append(
            MicroParameters(
                wavelength=[WAVELENGTH_UM],
                n_real=N_REAL,
                n_imag=MINERAL_N_IMAG if spheroidal else WASO_INDEX[1],
                shape=(
                    SpheroidDistrFile(distr_filename=AR_KANDLER)
                    if spheroidal
                    else Sphere()
                ),
                psd=LognormalPSD(
                    rm=rm,
                    sigma=sigma,
                    n=concentration * 1e6,
                    rmin=R_MIN_UM,
                    rmax=R_MAX_UM,
                ),
                # The water soluble mode keeps its own index: the variability
                # section 5.6 studies is the dust one.
                nonabs_fraction=nonabs_fraction if spheroidal else 0.0,
            )
        )
    return built


def compute() -> xr.Dataset:
    """
    The volume scattering function of each case.

    Returns
    -------
    xr.Dataset
        Dimensions ``(case, theta)`` for the scattering function, ``(case,)``
        for the four integrated quantities the script reports.
    """
    per_case = []
    for name, (fraction, *_) in CASES.items():
        point = run_point(
            _modes(fraction), OUTPUTS, quiet=True, n_angles=N_ANGLES
        )
        scalar = {
            key: float(point[key].values.ravel()[0])
            for key in ("ssa", "g", "lidar_ratio", "depol_ratio")
        }
        # The volume scattering function is the phase function scaled by the
        # scattering coefficient; the figure is on an arbitrary scale, so the
        # scattering coefficient is what makes the four comparable.
        scattering = float(point["ksca"].values.ravel()[0])
        per_case.append(
            xr.Dataset(
                {
                    "vsf": (
                        ("theta",),
                        point["phase"].values.ravel()
                        * scattering
                        / (4.0 * np.pi),
                    ),
                    **scalar,
                },
                coords={"theta": point["theta"].values},
            )
        )
        print(f"  {name}: done")
    return xr.concat(
        per_case,
        dim=xr.DataArray(list(CASES), dims="case", name="case"),
    )


def report(data: xr.Dataset) -> None:
    """Print what section 5.6 says, beside the computed values."""
    print("\n=== section 5.6 ===")
    print(f"{'':18s}{'omega_0':>7s}{'g':>7s}{'S [sr]':>9s}{'delta_l':>9s}")
    for case in data["case"].values:
        row = data.sel(case=case)
        print(
            f"{SHORT[str(case)]:18s}{float(row['ssa']):7.4f}"
            f"{float(row['g']):7.4f}{float(row['lidar_ratio']):9.2f}"
            f"{float(row['depol_ratio']):9.4f}"
        )
    print()
    for label, row in PUBLISHED.items():
        name = f"published, {SHORT[label].replace('average m_i', 'avg')}"
        print(f"{name:18s}{row[2:]}")

    backward = data["vsf"].sel(theta=slice(150.0, 180.0))
    spread = backward.max("case") / backward.min("case")
    print(
        f"\nspread over 150 to 180 deg: at most a factor"
        f" {float(spread.max()):.2f}"
    )


def plot(data: xr.Dataset, path: Path) -> None:
    """Draw the volume scattering functions."""
    import matplotlib.pyplot as plt

    # The published scale is arbitrary. Putting the forward value of the
    # black curve on its 30 lets the two be laid side by side.
    scale = 30.0 / float(data["vsf"].sel(case=list(CASES)[0]).isel(theta=0))

    fig, axis = plt.subplots(figsize=(8, 6))
    for name, (_, colour, style, width) in CASES.items():
        axis.semilogy(
            data["theta"],
            data["vsf"].sel(case=name) * scale,
            style,
            color=colour,
            linewidth=width,
            label=name,
        )
    axis.set_xlim(0, 180)
    axis.set_ylim(2e-3, 5e1)
    axis.set_xticks(np.arange(0, 181, 20))
    axis.set_xlabel("Scattering angle $\\theta$")
    axis.set_ylabel("Volume scattering function $\\tilde{a}_1$ (arb. scale)")
    axis.grid(which="both", alpha=0.3)
    axis.legend()
    fig.suptitle("Gasteiger and Wiegner (2018), Figure 9, recomputed")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    print(f"\nwrote {path}")


def main() -> None:
    print("Computing the dust ensemble under four index treatments:")
    data = compute()
    report(data)
    plot(data, FIGURES / "gasteiger_fig9.png")


if __name__ == "__main__":
    main()
