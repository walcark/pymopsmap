"""
Reproduce Figure 6 and Table 4 of Gasteiger and Wiegner (2018), GMD 11:2739.

    PYMOPSMAP_DATASET_SOURCE=... \
        pixi run -e dev python -m scripts.validation.gasteiger_fig6

"Phase functions at lambda = 500 nm of the five COSMO-MUSCAT dust size bins
(different colors) assuming spherical particles (solid lines) and prolate
spheroids (dashed lines)."

Section 5.2 fixes everything the figure needs: the bins are "determined by the
radius limits 0.1, 0.3, 0.9, 2.6, 8, and 24 um", with "constant dv/dlnr within
each bin", a refractive index m = 1.53 + 0.0078i, the aspect ratio
distribution of Kandler et al. (2009), and volume-equivalent sizes for the
spheroids. Table 4 gives four numbers per bin and shape, which the script
prints beside its own.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

from pymopsmap import MicroParameters
from pymopsmap.engine import run_point
from pymopsmap.engine.outputs import OutputType
from pymopsmap.psd import DistrListPSD, DistrType
from pymopsmap.shapes import Sphere, SpheroidDistrFile

FIGURES = Path(__file__).resolve().parents[2] / "docs" / "figures"

WAVELENGTH_UM = 0.5
REFRACTIVE_INDEX = (1.53, 0.0078)
BIN_EDGES = [0.1, 0.3, 0.9, 2.6, 8.0, 24.0]
AR_KANDLER = str(
    Path(__file__).resolve().parents[2] / "bin/mopsmap/data/ar_kandler"
)

OUTPUTS = frozenset(
    {OutputType.INTEGRATED, OutputType.LIDAR, OutputType.PHASE_FUNCTION}
)

# Table 4, as published: omega_0, g, eta (g m-2) and Z (m2 sr-1 g-1), spheres
# then prolate spheroids.
TABLE4 = {
    ("spheres", 1): (0.9632, 0.6567, 0.2905, 4.234e-2),
    ("spheres", 2): (0.9216, 0.6866, 0.5594, 1.185e-1),
    ("spheres", 3): (0.7903, 0.8088, 2.230, 1.403e-2),
    ("spheres", 4): (0.6450, 0.8998, 6.989, 1.204e-3),
    ("spheres", 5): (0.5561, 0.9442, 22.09, 8.225e-5),
    ("spheroids", 1): (0.9628, 0.6585, 0.3000, 3.981e-2),
    ("spheroids", 2): (0.9264, 0.7111, 0.5236, 5.421e-2),
    ("spheroids", 3): (0.7934, 0.8109, 2.071, 8.901e-3),
    ("spheroids", 4): (0.6485, 0.9017, 6.633, 7.457e-4),
    ("spheroids", 5): (0.5601, 0.9419, 20.90, 8.651e-5),
}

COLOURS = ["black", "tab:blue", "tab:green", "tab:orange", "tab:red"]


def compute() -> xr.Dataset:
    """
    Every bin, both shapes.

    Returns
    -------
    xr.Dataset
        Dimensions ``(shape, bin, theta)`` for the phase function and
        ``(shape, bin)`` for the four tabulated quantities.
    """
    per_shape = []
    for shape_name in ("spheres", "spheroids"):
        per_bin = [
            _bin(shape_name, index) for index in range(1, len(BIN_EDGES))
        ]
        per_shape.append(
            xr.concat(
                per_bin,
                dim=xr.DataArray(
                    np.arange(1, len(BIN_EDGES)), dims="bin", name="bin"
                ),
            )
        )
    return xr.concat(
        per_shape,
        dim=xr.DataArray(["spheres", "spheroids"], dims="shape", name="shape"),
    )


def _bin(shape_name: str, index: int) -> xr.Dataset:
    """One bin of one shape, as the figure and the table need it."""
    spherical = shape_name == "spheres"
    mode = MicroParameters(
        wavelength=[WAVELENGTH_UM],
        n_real=REFRACTIVE_INDEX[0],
        n_imag=REFRACTIVE_INDEX[1],
        shape=(
            Sphere()
            if spherical
            else SpheroidDistrFile(distr_filename=AR_KANDLER)
        ),
        psd=DistrListPSD(
            radii=[BIN_EDGES[index - 1], BIN_EDGES[index]],
            concentrations=[1.0, 1.0],
            distr_type=DistrType.DVDLNR,
        ),
        # OPAC gives the mineral components a density of 2.6 g cm-3, which is
        # what turns an extinction coefficient into the tabulated eta.
        density=2.6,
    )
    point = run_point(
        [mode], OUTPUTS, quiet=True, size_equ="cs" if spherical else "vol"
    )
    return xr.Dataset(
        {
            "phase": (("theta",), point["phase"].values.ravel()),
            "ssa": float(point["ssa"].values.ravel()[0]),
            "g": float(point["g"].values.ravel()[0]),
            "eta": float(point["ext_to_mass"].values.ravel()[0]),
            "Z": float(point["back_to_mass"].values.ravel()[0]),
        },
        coords={"theta": point["theta"].values},
    )


def report(data: xr.Dataset) -> None:
    """Print Table 4 beside the recomputed values."""
    print(f"\n{'':22s}{'published':>12s}{'computed':>12s}{'rel. dev.':>12s}")
    for (shape_name, index), reference in TABLE4.items():
        for name, expected in zip(("ssa", "g", "eta", "Z"), reference):
            obtained = float(
                data[name].sel(shape=shape_name, bin=index).values
            )
            deviation = obtained / expected - 1.0
            label = f"{shape_name[:8]:8s} bin {index} {name:>4s}"
            print(
                f"{label:22s}{expected:12.4g}{obtained:12.4g}{deviation:11.2%}"
            )


def plot(data: xr.Dataset, path: Path) -> None:
    """Draw the ten phase functions."""
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(figsize=(8, 7))
    for index, colour in zip(data["bin"].values, COLOURS):
        for shape_name, style in (("spheres", "-"), ("spheroids", "--")):
            axis.semilogy(
                data["theta"],
                data["phase"].sel(shape=shape_name, bin=index),
                style,
                color=colour,
                label=f"Bin {index}, {shape_name}",
            )
    axis.set_xlim(0, 180)
    axis.set_xticks(np.arange(0, 181, 20))
    axis.set_xlabel("scattering angle [deg]")
    axis.set_ylabel("phase function")
    axis.grid(alpha=0.3)
    axis.legend(fontsize=8, ncol=2)
    fig.suptitle("Gasteiger and Wiegner (2018), Figure 6, recomputed")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    print(f"\nwrote {path}")


def main() -> None:
    print("Computing the five COSMO-MUSCAT dust bins:")
    data = compute()
    report(data)
    plot(data, FIGURES / "gasteiger_fig6.png")


if __name__ == "__main__":
    main()
