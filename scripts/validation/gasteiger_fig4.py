"""
Reproduce Figure 4 of Gasteiger and Wiegner (2018), GMD 11:2739.

    PYMOPSMAP_DATASET_SOURCE=... \
        pixi run -e dev python -m scripts.validation.gasteiger_fig4

"Examples illustrating the effect of the limited size resolution of the
MOPSMAP data set (a, c) and the effect of the interpolation between the
refractive index grid points of the data set (b, d)."

This is the figure about the data set's own error, so it is the one that needs
something outside MOPSMAP to compare against. ``miepython`` provides it.

Panels (b) and (d) are complete: the four grid points the data set holds
around m = 1.54 + 0.005i, the value MOPSMAP interpolates there, and Mie theory
computed explicitly on the same size grid.

Panels (a) and (c) show the spheres only. Their second pair of curves is a
T-matrix calculation at high size resolution for prolate spheroids, which no
Python package provides and which the data set exists to avoid.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

from pymopsmap import MicroParameters
from pymopsmap.engine import run_point
from pymopsmap.engine.outputs import DEFAULT_OUTPUT
from pymopsmap.psd import FixedPSD
from pymopsmap.shapes import Sphere
from pymopsmap.utils import DATASET_CACHE_DIR

FIGURES = Path(__file__).resolve().parents[2] / "docs" / "figures"

WAVELENGTH_UM = 0.5
CONCENTRATION_M3 = 1e6
X_MAX = 40.0

# The index of panels (a) and (c), a grid point of the data set, so nothing is
# interpolated and the only error left is the size sampling.
SAMPLING_INDEX = (1.52, 0.0)

# The index of panels (b) and (d), which sits between four grid points. The
# caption names all four.
TARGET_INDEX = (1.54, 0.005)
CORNERS = [
    (1.52, 0.0043),
    (1.52, 0.0060811),
    (1.56, 0.0043),
    (1.56, 0.0060811),
]

PANELS = {"qext": "Extinction efficiency $q_{ext}$", "g": "Asymmetry para. g"}


def dataset_sizes() -> np.ndarray:
    """
    The size parameter grid the data set itself uses, up to x = 40.

    Reading it from the file is what "the same x grid as used by the data set"
    means; it rises in one percent steps.
    """
    path = DATASET_CACHE_DIR / "spheres" / "sphere_1.5200_0.004300.nc"
    with xr.open_dataset(path) as ds:
        sizes = ds["sizepara"].values
    return sizes[(sizes >= 0.2) & (sizes <= X_MAX)]


def _curve(grid: np.ndarray, n_real: float, n_imag: float) -> xr.Dataset:
    """Run one refractive index over the size grid, one particle at a time."""
    values = {name: np.full(len(grid), np.nan) for name in PANELS}
    for index, x in enumerate(grid):
        point = run_point(
            [
                MicroParameters(
                    wavelength=[WAVELENGTH_UM],
                    n_real=n_real,
                    n_imag=n_imag,
                    shape=Sphere(),
                    psd=FixedPSD(
                        radius=x * WAVELENGTH_UM / (2.0 * np.pi),
                        n=CONCENTRATION_M3,
                    ),
                )
            ],
            DEFAULT_OUTPUT,
            quiet=True,
        )
        values["qext"][index] = float(point["kext"].values.ravel()[0]) / float(
            point["cross_dens"].values.ravel()[0]
        )
        values["g"][index] = float(point["g"].values.ravel()[0])
    return xr.Dataset(
        {name: (("x",), value) for name, value in values.items()},
        coords={"x": grid},
    )


def _mie(grid: np.ndarray, n_real: float, n_imag: float) -> xr.Dataset:
    """The same quantities from miepython, which shares no code with it."""
    import miepython

    # miepython takes m = n - ik, the opposite sign convention.
    q_ext, _, _, asymmetry = miepython.efficiencies_mx(
        complex(n_real, -n_imag), grid
    )
    return xr.Dataset(
        {
            "qext": (("x",), np.asarray(q_ext)),
            "g": (("x",), np.asarray(asymmetry)),
        },
        coords={"x": grid},
    )


def compute() -> dict[str, xr.Dataset]:
    """
    Every curve of the figure that can be drawn.

    Returns
    -------
    dict of str to xr.Dataset
        ``sampling_dataset`` and ``sampling_mie`` for the left column,
        ``interpolated``, ``explicit`` and one entry per grid corner for the
        right one.
    """
    grid = dataset_sizes()
    fine = np.arange(0.2, X_MAX + 0.002, 0.002)
    print(f"  {len(grid)} data set sizes, {len(fine)} for the fine Mie grid")

    curves = {
        "sampling_dataset": _curve(grid, *SAMPLING_INDEX),
        "sampling_mie": _mie(fine, *SAMPLING_INDEX),
        "interpolated": _curve(grid, *TARGET_INDEX),
        "explicit": _mie(grid, *TARGET_INDEX),
    }
    print("  the four curves of the left and right columns: done")
    for n_real, n_imag in CORNERS:
        curves[f"corner {n_real} {n_imag}"] = _curve(grid, n_real, n_imag)
        print(f"  grid point {n_real} + {n_imag}i: done")
    return curves


def report(curves: dict[str, xr.Dataset]) -> None:
    """State how far each approximation is from its reference."""
    print("\n=== deviation from the explicit calculation ===")
    for name in PANELS:
        error = curves["interpolated"][name] - curves["explicit"][name]
        relative = abs(error / curves["explicit"][name])
        print(
            f"{name:5s} interpolation: at most {float(relative.max()):.1%},"
            f" typically {float(relative.median()):.2%}"
        )

    print("\n=== deviation from the fine size grid ===")
    fine = curves["sampling_mie"]
    coarse = curves["sampling_dataset"]
    for name in PANELS:
        resampled = fine[name].interp(x=coarse["x"])
        relative = abs((coarse[name] - resampled) / resampled)
        print(
            f"{name:5s} size sampling: at most {float(relative.max()):.1%},"
            f" typically {float(relative.median()):.2%}"
        )


def plot(curves: dict[str, xr.Dataset], path: Path) -> None:
    """Draw the four panels in the layout of the published figure."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    limits = {"qext": (0.0, 5.0), "g": (0.5, 0.90)}

    for row, name in enumerate(PANELS):
        left, right = axes[row]
        left.plot(
            curves["sampling_mie"]["x"],
            curves["sampling_mie"][name],
            color="black",
            linewidth=0.6,
            label="Spheres at high size resolution",
        )
        left.plot(
            curves["sampling_dataset"]["x"],
            curves["sampling_dataset"][name],
            color="blue",
            linewidth=0.8,
            label="Spheres from data set",
        )

        for key in curves:
            if key.startswith("corner"):
                right.plot(
                    curves[key]["x"],
                    curves[key][name],
                    color="0.7",
                    linewidth=0.5,
                )
        right.plot(
            curves["explicit"]["x"],
            curves["explicit"][name],
            color="black",
            linewidth=0.9,
            label="Mie theory for m=1.54+0.005i",
        )
        right.plot(
            curves["interpolated"]["x"],
            curves["interpolated"][name],
            color="red",
            linewidth=0.9,
            label="Data set interpolated for m=1.54+0.005i",
        )

        for axis, title in ((left, "(a)"), (right, "(b)")):
            axis.set_ylim(*limits[name])
            axis.set_xlim(0.0, X_MAX)
            axis.grid(alpha=0.3)
            axis.text(
                0.02,
                0.93,
                title if row == 0 else {"(a)": "(c)", "(b)": "(d)"}[title],
                transform=axis.transAxes,
            )
        left.set_ylabel(PANELS[name])

    axes[0, 0].legend(fontsize=7, loc="upper right")
    axes[0, 1].legend(fontsize=7, loc="upper right")
    axes[0, 0].text(18.0, 3.4, "m=1.52+0i", fontsize=9)
    axes[0, 1].text(14.0, 3.4, "Spheres", fontsize=9)
    for axis in axes[1]:
        axis.set_xlabel("Size parameter $x_c$")

    fig.suptitle("Gasteiger and Wiegner (2018), Figure 4, recomputed")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    print(f"\nwrote {path}")


def main() -> None:
    print("Computing the sampling and interpolation curves:")
    curves = compute()
    report(curves)
    plot(curves, FIGURES / "gasteiger_fig4.png")


if __name__ == "__main__":
    main()
