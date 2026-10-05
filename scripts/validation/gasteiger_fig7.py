"""
Reproduce Figure 7 of Gasteiger and Wiegner (2018), GMD 11:2739.

    PYMOPSMAP_DATASET_SOURCE=... \
        pixi run -e dev python -m scripts.validation.gasteiger_fig7

"Optical and microphysical properties of the OPAC desert aerosol type as a
function of cutoff radius rmax", normalised on the values at rmax = 60 um.

The four modes, their concentrations and their shapes are the ones of
``bin/mopsmap/misc/paper_examples/sect_53_rmax_cutoff/calc_cutoff_optics.py``,
the script the authors ship with MOPSMAP: water soluble as spheres, the three
mineral modes as spheroids with the aspect ratio distribution of Kandler et
al. (2009). Section 5.3 states six numbers about the result, which the script
prints beside its own.
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

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "bin" / "mopsmap" / "data"
AR_KANDLER = str(DATA / "ar_kandler")

WAVELENGTHS_UM = [0.4, 0.6, 1.0]
R_MIN_UM = 0.01
REFERENCE_RMAX_UM = 60.0

# n [cm-3], r_mod [um], sigma, refractive index file, spheroidal.
MODES = [
    (6.67, 0.0212, 2.24, "refr_waso00", False),
    (0.898, 0.07, 1.95, "refr_mineral", True),
    (0.1016, 0.39, 2.00, "refr_mineral", True),
    (0.000472, 1.9, 2.15, "refr_mineral", True),
]

OUTPUTS = frozenset({OutputType.INTEGRATED})

PANELS = {
    "kext": "normalized extinction coefficient",
    "ssa": "single scat. alb. $\\omega_0$",
    "g": "asymmetry para. $g$",
}

# The cutoffs section 5.3 names: PM2.5, PM10, and the 10 um example.
MARKERS = {"PM2.5": 1.25, "PM10": 5.0, "r = 10 um": 10.0}


def cutoffs(count: int = 200) -> np.ndarray:
    """
    The cutoff grid of the reference script: 0.5 to 60 um, log spaced.

    The three cutoffs section 5.3 quotes are added to it, so the percentages
    the report prints are read at the radius the article names rather than at
    the nearest grid point.
    """
    grid = np.logspace(np.log10(0.5), np.log10(REFERENCE_RMAX_UM), count)
    return np.unique(np.concatenate([grid, list(MARKERS.values())]))


def _refractive_index(filename: str) -> tuple[list[float], list[float]]:
    """Interpolate one MOPSMAP refractive index file onto the wavelengths."""
    table = np.loadtxt(DATA / filename, comments="#")
    real = np.interp(WAVELENGTHS_UM, table[:, 0], table[:, 1])
    imaginary = np.interp(WAVELENGTHS_UM, table[:, 0], table[:, 2])
    return real.tolist(), imaginary.tolist()


def _modes(rmax: float) -> list[MicroParameters]:
    """The four desert modes, cut off at one maximum radius."""
    built = []
    for concentration, rm, sigma, filename, spheroidal in MODES:
        real, imaginary = _refractive_index(filename)
        built.append(
            MicroParameters(
                wavelength=WAVELENGTHS_UM,
                n_real=real,
                n_imag=imaginary,
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
                    rmax=rmax,
                ),
                density=2.6,
            )
        )
    return built


def compute() -> xr.Dataset:
    """
    The whole ensemble at every cutoff.

    Returns
    -------
    xr.Dataset
        Dimensions ``(rmax, wl)``, with the three plotted quantities plus the
        cross section density and the mass concentration.
    """
    grid = cutoffs()
    names = ["kext", "ssa", "g", "cross_dens", "mass_conc"]
    values = {
        name: np.full((len(grid), len(WAVELENGTHS_UM)), np.nan)
        for name in names
    }
    for index, rmax in enumerate(grid):
        point = run_point(_modes(rmax), OUTPUTS, quiet=True)
        for name in names:
            values[name][index] = point[name].values
        if index % 40 == 0:
            print(f"  rmax = {rmax:6.2f} um")

    return xr.Dataset(
        {name: (("rmax", "wl"), value) for name, value in values.items()},
        coords={"rmax": grid, "wl": WAVELENGTHS_UM},
    )


def _normalised(data: xr.Dataset, name: str) -> xr.DataArray:
    """One quantity divided by its value at the largest cutoff."""
    return data[name] / data[name].isel(rmax=-1)


def report(data: xr.Dataset) -> None:
    """Print the statements section 5.3 makes, beside the computed values."""
    mass = _normalised(data, "mass_conc").isel(wl=0)
    area = _normalised(data, "cross_dens").isel(wl=0)
    for label, cutoff in (("PM10", 5.0), ("PM2.5", 1.25)):
        print(
            f"\n{label} (rmax = {cutoff} um):"
            f" {float(mass.sel(rmax=cutoff, method='nearest')):.1%} of the"
            f" mass, {float(area.sel(rmax=cutoff, method='nearest')):.1%}"
            " of the cross section"
        )
    print(
        "  article: 59.5 % and 21.6 % of the mass, 94.4 % and 69.0 % of"
        " the cross section"
    )

    print(
        f"\nat rmax = 10 um:"
        f" {float(area.sel(rmax=10.0, method='nearest')):.1%} of the cross"
        f" section, {float(mass.sel(rmax=10.0, method='nearest')):.1%}"
        " of the mass"
    )
    print("  article: 97.8 % and 75.6 %")

    ssa = data["ssa"]
    shift = ssa.sel(rmax=1.25, method="nearest") - ssa.isel(rmax=-1)
    print(f"\nPM2.5 raises omega_0 by {shift.values.round(4)}")
    print("  article: 'about 0.035-0.071 higher'")
    shift_g = data["g"].sel(rmax=1.25, method="nearest") - data["g"].isel(
        rmax=-1
    )
    print(f"PM2.5 lowers g by {(-shift_g).values.round(4)}")
    print("  article: 'reduced by about 0.02-0.04'")


def plot(data: xr.Dataset, path: Path) -> None:
    """Draw the three panels of the figure."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    colours = ["b", "g", "r"]

    for wavelength, colour in zip(WAVELENGTHS_UM, colours):
        axes[0].semilogx(
            data["rmax"],
            _normalised(data, "kext").sel(wl=wavelength),
            color=colour,
            label=f"$\\alpha_{{ext}}$ at {wavelength * 1e3:.0f} nm",
        )
    axes[0].semilogx(
        data["rmax"],
        _normalised(data, "cross_dens").isel(wl=0),
        color="k",
        label="a - 'area'",
    )
    axes[0].semilogx(
        data["rmax"],
        _normalised(data, "mass_conc").isel(wl=0),
        color="grey",
        label="M - 'mass'",
    )
    axes[0].set_ylim(0.0, 1.0)
    axes[0].set_ylabel("normalized value")
    axes[0].legend(fontsize=8)

    for axis, name in ((axes[1], "ssa"), (axes[2], "g")):
        for wavelength, colour in zip(WAVELENGTHS_UM, colours):
            axis.semilogx(
                data["rmax"], data[name].sel(wl=wavelength), color=colour
            )
        axis.set_ylabel(PANELS[name])

    for axis in axes:
        for cutoff in MARKERS.values():
            axis.axvline(x=cutoff, color="y")
        axis.grid(which="both", alpha=0.3)
    axes[1].set_ylim(0.76, 1.0)
    axes[2].set_ylim(0.6, 0.8)
    axes[2].set_xlabel("cut-off radius $r_{max}$ [$\\mu$m]")
    axes[2].set_xlim(float(data["rmax"][0]), float(data["rmax"][-1]))

    fig.suptitle("Gasteiger and Wiegner (2018), Figure 7, recomputed")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    print(f"\nwrote {path}")


def main() -> None:
    print("Computing the OPAC desert type over the cutoff grid:")
    data = compute()
    report(data)
    plot(data, FIGURES / "gasteiger_fig7.png")


if __name__ == "__main__":
    main()
