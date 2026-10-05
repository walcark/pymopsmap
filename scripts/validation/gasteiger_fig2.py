"""
Reproduce Figure 2 of Gasteiger and Wiegner (2018), GMD 11:2739.

    PYMOPSMAP_DATASET_SOURCE=... \
        pixi run -e dev python -m scripts.validation.gasteiger_fig2

"Optical properties of single particles (or narrow size bins in the case of
spheres) with fixed refractive index m = 1.56 + 0.00215i as a function of size
parameter." Three panels: extinction efficiency, single scattering albedo,
asymmetry parameter, for five shapes.

The refractive index is a grid point of the data set, so nothing is
interpolated and the curves are the data set itself, read through the wrapper.
The article plots xc from 1 to 1000; the irregular shapes stop at 30.2, which
is where their own coverage ends (Table 2).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

from pymopsmap import MicroParameters
from pymopsmap.engine import run_point
from pymopsmap.engine.outputs import DEFAULT_OUTPUT
from pymopsmap.exceptions import OutsideCoverageError
from pymopsmap.psd import FixedPSD
from pymopsmap.shapes import Irregular, Shape, Sphere, Spheroid

FIGURES = Path(__file__).resolve().parents[2] / "docs" / "figures"

REFRACTIVE_INDEX = (1.56, 0.00215)
WAVELENGTH_UM = 0.5
CONCENTRATION_M3 = 1e6

# The five curves of the figure, in the colours the caption gives them.
SHAPES: dict[str, tuple[Shape, str]] = {
    "spheres": (Sphere(), "tab:blue"),
    "prolate, 1.4": (Spheroid(mode="prolate", aspect_ratio=1.4), "tab:orange"),
    "prolate, 3.0": (Spheroid(mode="prolate", aspect_ratio=3.0), "tab:green"),
    "irregular D": (Irregular(shape_id="D"), "tab:red"),
    "irregular F": (Irregular(shape_id="F"), "tab:purple"),
}

PANELS = {
    "qext": "extinction efficiency $q_{ext}$",
    "ssa": "single scattering albedo $\\omega_0$",
    "g": "asymmetry parameter $g$",
}


def size_parameters(count: int = 500) -> np.ndarray:
    """A log grid of cross-section-equivalent size parameters, 1 to 1000."""
    return np.logspace(0.0, 3.0, count)


def compute() -> xr.Dataset:
    """
    One curve per shape over the size parameter grid.

    Returns
    -------
    xr.Dataset
        Dimensions ``(shape, x)``, variables ``qext``, ``ssa`` and ``g``.
        Positions a shape does not cover hold NaN.
    """
    grid = size_parameters()
    per_shape: dict[str, xr.Dataset] = {}
    for name, (shape, _) in SHAPES.items():
        per_shape[name] = _curve(shape, grid)
        covered = int(np.isfinite(per_shape[name]["qext"]).sum())
        print(f"  {name}: {covered}/{len(grid)} points covered")
    return xr.concat(
        list(per_shape.values()),
        dim=xr.DataArray(list(per_shape), dims="shape", name="shape"),
    )


def _curve(shape: Shape, grid: np.ndarray) -> xr.Dataset:
    """Run one shape over the grid, one single particle at a time."""
    values = {name: np.full(len(grid), np.nan) for name in PANELS}
    for index, x in enumerate(grid):
        radius = x * WAVELENGTH_UM / (2.0 * np.pi)
        mode = MicroParameters(
            wavelength=[WAVELENGTH_UM],
            n_real=REFRACTIVE_INDEX[0],
            n_imag=REFRACTIVE_INDEX[1],
            shape=shape,
            psd=FixedPSD(radius=radius, n=CONCENTRATION_M3),
        )
        try:
            point = run_point([mode], DEFAULT_OUTPUT, quiet=True)
        except OutsideCoverageError:
            continue
        # The extinction efficiency is the extinction coefficient over the
        # geometric cross section the same run reports.
        values["qext"][index] = _scalar(point, "kext") / _scalar(
            point, "cross_dens"
        )
        values["ssa"][index] = _scalar(point, "ssa")
        values["g"][index] = _scalar(point, "g")

    return xr.Dataset(
        {name: (("x",), value) for name, value in values.items()},
        coords={"x": grid},
    )


def _scalar(point: xr.Dataset, name: str) -> float:
    """The single wavelength value of one output variable."""
    return float(point[name].values.ravel()[0])


def report(data: xr.Dataset) -> None:
    """Print the three statements the article makes about the curves."""
    ssa = data["ssa"]
    print(f"\nmaximum omega_0 over all shapes: {float(ssa.max()):.4f}")
    print("  article: 'maxima with values of about 0.991'")

    tail = ssa.sel(shape=["spheres", "prolate, 1.4", "prolate, 3.0"]).sel(
        x=1000.0, method="nearest"
    )
    print(
        f"omega_0 at xc = 1000, spheres and spheroids: {tail.values.round(4)}"
    )
    print("  article: 'approaches a value of about 0.551 at xc = 1000'")


def plot(data: xr.Dataset, path: Path) -> None:
    """Draw the three panels of the figure."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(7, 10), sharex=True)
    for axis, (name, label) in zip(axes, PANELS.items()):
        for shape_name, (_, colour) in SHAPES.items():
            axis.semilogx(
                data["x"],
                data[name].sel(shape=shape_name),
                color=colour,
                label=shape_name,
            )
        axis.set_ylabel(label)
        axis.grid(which="both", alpha=0.3)
    axes[0].legend(fontsize=8)
    axes[-1].set_xlabel("size parameter $x_c$")
    fig.suptitle("Gasteiger and Wiegner (2018), Figure 2, recomputed")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    print(f"\nwrote {path}")


def main() -> None:
    print("Computing the single-particle curves:")
    data = compute()
    report(data)
    plot(data, FIGURES / "gasteiger_fig2.png")


if __name__ == "__main__":
    main()
