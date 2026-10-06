"""
Reproduce Figure 11 of Gasteiger and Wiegner (2018), GMD 11:2739.

    PYMOPSMAP_DATASET_SOURCE=... VOGEL_DATA=~/downloads/jgrd54020 \
        pixi run -e dev python -m scripts.validation.gasteiger_fig11

"Modeled wavelength-dependent optical properties for ashes from different
volcanoes."

Section 5.8 builds the ensembles one measured particle at a time: "Each single
particle is modeled as a prolate spheroid with the given size and aspect ratio,
as well as with the refractive index given for the type of ash the volcano
emits. In addition, we assume a non-absorbing fraction of X = 0.5."

The measurements are the supporting information of

    Vogel, A. et al., "Reference dataset of volcanic ash physicochemical and
    optical properties", J. Geophys. Res. Atmos. 122, 2017,
    doi:10.1002/2016JD026328

which has to be fetched from the publisher. Point ``VOGEL_DATA`` at the
directory holding its sixteen ``jgrd54020-sup-*.txt`` files.

Particles above r = 47.5 um and aspect ratios above 5 are set to those values,
as section 5.8 says, "for each volcano, less than 0.5 % of the particles was
affected by these modifications".
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import xarray as xr

from pymopsmap import MicroParameters, Specie
from pymopsmap.engine import run_point
from pymopsmap.engine.outputs import DEFAULT_OUTPUT
from pymopsmap.species import Mode

ROOT = Path(__file__).resolve().parents[2]
FIGURES = ROOT / "docs" / "figures"
REFRACTIVE_DIR = (
    ROOT / "bin/mopsmap/misc/paper_examples/sect_58_volcano_optics"
)
VOGEL_DATA = Path(
    os.getenv("VOGEL_DATA", Path.home() / "downloads" / "jgrd54020")
)

# The aspect ratio grid of the data set, which is what every measured ratio is
# spread over. MOPSMAP interpolates between its points, so a particle lands on
# the two that bracket it.
AR_GRID = np.array(
    [
        1.0,
        1.2,
        1.4,
        1.6,
        1.8,
        2.0,
        2.2,
        2.4,
        2.6,
        2.8,
        3.0,
        3.4,
        3.8,
        4.2,
        4.6,
        5.0,
    ]
)

R_MAX_UM = 47.5
AR_MAX = 5.0
NONABS_FRACTION = 0.5
WAVELENGTHS_UM = np.round(np.arange(0.3, 1.5001, 0.025), 4)

# The nine volcanoes of the figure: data file, label, ash type, and the colour
# and line style make_plot.py gives each.
VOLCANOES = [
    ("ds05", "Grimsvotn", "basalt", "#002fd5", "-"),
    ("ds06", "Mount Kelud", "basalt_andesite", "#00a7ff", "-"),
    ("ds07", "Mount Sakurajima", "andesite", "#00ce3d", "-"),
    ("ds08", "Eyjafjallajokull", "andesite", "#00ce3d", "--"),
    ("ds09", "Mount Spurr", "andesite", "#00ce3d", "-."),
    ("ds10", "Mount Redoubt", "dacite", "#ffad00", "-"),
    ("ds11", "Soufriere Hills", "dacite", "#ffad00", "--"),
    ("ds12", "Mount St. Helens", "dacite", "#ffad00", "-."),
    ("ds13", "Chaiten", "rhyolite", "#ff0000", "-"),
]


def _refractive_index(ash: str) -> tuple[list[float], list[float]]:
    """Interpolate one ash refractive index onto the wavelength grid."""
    table = np.loadtxt(REFRACTIVE_DIR / f"refrac_{ash}")
    real = np.interp(WAVELENGTHS_UM, table[:, 0], table[:, 1])
    imaginary = np.interp(WAVELENGTHS_UM, table[:, 0], table[:, 2])
    return real.tolist(), imaginary.tolist()


def _measurements(dataset: str) -> tuple[np.ndarray, np.ndarray]:
    """
    Read one volcano: a radius and an aspect ratio per measured particle.

    Column 0 is the circle-equivalent diameter, column 14 the ratio of the
    minor to the major axis, which is the reciprocal of the one MOPSMAP
    takes. The header is 27 lines and carries a non-UTF-8 volcano name.
    """
    pattern = f"*{dataset}.txt"
    path = next(iter(sorted(VOGEL_DATA.glob(pattern))), None)
    if path is None:
        raise SystemExit(
            f"no {pattern} under {VOGEL_DATA}. The supporting information of "
            "Vogel et al. (2017) has to be fetched from the publisher; set "
            "VOGEL_DATA to the directory holding it."
        )
    table = np.loadtxt(path, skiprows=27, encoding="latin-1")
    return table[:, 0] * 0.5, 1.0 / table[:, 14]


def _modes(dataset: str, ash: str) -> list[MicroParameters]:
    """One mode per measured size and aspect ratio."""
    radii, ratios = _measurements(dataset)
    real, imaginary = _refractive_index(ash)
    modes = Mode.from_particles(
        radii=radii,
        aspect_ratios=ratios,
        n_real=real,
        n_imag=imaginary,
        r_max=R_MAX_UM,
        nonabs_fraction=NONABS_FRACTION,
    )
    specie = Specie.custom(modes, name=dataset, wl=WAVELENGTHS_UM.tolist())
    return specie.at(wl=WAVELENGTHS_UM.tolist()).modes


def compute() -> xr.Dataset:
    """
    The albedo and the asymmetry parameter of the nine ashes.

    The modes are built one volcano at a time: that part is Python, holds the
    interpreter lock and costs memory, so running it nine times over would buy
    nothing. The MOPSMAP calls are what goes out in parallel, each one a
    subprocess of its own.

    The cost is the number of contributions, one NetCDF read each: nine
    thousand modes over forty-nine wavelengths, each landing between four
    refractive index grid points and split again by the non-absorbing
    fraction, is about half an hour of one core per volcano.

    Returns
    -------
    xr.Dataset
        Dimensions ``(volcano, wl)``.
    """
    prepared = []
    for dataset, label, ash, *_ in VOLCANOES:
        modes = _modes(dataset, ash)
        print(f"  {label:18s} {ash:16s} {len(modes):6d} modes")
        prepared.append(modes)

    workers = min(len(VOLCANOES), os.cpu_count() or 1)
    print(f"\n  running {workers} volcanoes at a time")
    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(_run, prepared))

    for (_, label, *_), result in zip(VOLCANOES, results):
        print(f"  {label:18s} reff = {float(result['reff'][0]):5.2f} um")
    return xr.concat(
        [result[["ssa", "g", "reff"]] for result in results],
        dim=xr.DataArray(
            [label for _, label, *_ in VOLCANOES], dims="volcano"
        ),
    )


def _run(modes: list[MicroParameters]) -> xr.Dataset:
    """One MOPSMAP call, a subprocess, so it runs beside the others."""
    return run_point(modes, DEFAULT_OUTPUT, quiet=True)


def report(data: xr.Dataset) -> None:
    """Print what section 5.8 says about the spread between the ashes."""
    print("\n=== statements of section 5.8 ===")
    at_wl = data["ssa"].sel(wl=0.55, method="nearest")
    print(
        f"omega_0 at 550 nm spans {float(at_wl.max() - at_wl.min()):.3f},"
        f" from {str(at_wl.idxmin('volcano').values)}"
        f" to {str(at_wl.idxmax('volcano').values)}"
    )
    print(
        "  article: 'up to about 0.12 with ash from Chaiten and Mt. Kelud"
        " being the least and most absorbing'"
    )

    rise = data["ssa"].isel(wl=-1) - data["ssa"].isel(wl=0)
    print(f"\nomega_0 rises by {float(rise.mean()):.3f} over the range")
    print("  article: 'typically by about 0.05'")

    print(f"variability in g: {float(data['g'].std('volcano').max()):.3f}")
    print("  article: 'less than 0.05'")

    change = data["g"].max("wl") - data["g"].min("wl")
    print(f"change in g with wavelength, at most {float(change.max()):.3f}")
    print("  article: 'less than 0.02'")

    print(f"\neffective radii: {data['reff'].isel(wl=0).values.round(1)}")
    print("  article: 'in the range from 9.5 to 21 um'")


def plot(data: xr.Dataset, path: Path) -> None:
    """Draw the two panels in the layout of the published figure."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    for _, label, _, colour, style in VOLCANOES:
        series = data.sel(volcano=label)
        axes[0].plot(
            data["wl"] * 1e3, series["ssa"], style, color=colour, label=label
        )
        axes[1].plot(data["wl"] * 1e3, series["g"], style, color=colour)

    axes[0].set_yticks(np.arange(0.72, 1.00, 0.02))
    axes[0].set_ylim(0.8, 1.027)
    axes[0].set_ylabel("single scattering albedo $\\omega_0$")
    axes[0].legend(fontsize=8, ncol=2)
    axes[1].set_yticks(np.arange(0.6, 1.00, 0.01))
    axes[1].set_ylabel("asymmetry parameter $g$")
    axes[1].set_xticks(np.arange(300, 1600, 100))
    axes[1].set_xlim(299, 1501)
    axes[1].set_xlabel("wavelength $\\lambda$ [nm]")
    for axis in axes:
        axis.grid(which="major", alpha=0.5)
        axis.grid(which="minor", alpha=0.2)

    fig.suptitle("Gasteiger and Wiegner (2018), Figure 11, recomputed")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    print(f"\nwrote {path}")


def main() -> None:
    print("Computing the nine volcanic ashes:")
    data = compute()
    report(data)
    plot(data, FIGURES / "gasteiger_fig11.png")


if __name__ == "__main__":
    main()
