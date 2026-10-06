"""
Draw the figures the guide is illustrated with.

    PYMOPSMAP_DATASET_SOURCE=... \
        pixi run -e dev python -m scripts.demo.gallery

Each panel answers one question a reader of docs/guide.md has just been asked
to believe. They are regenerated from the library rather than drawn by hand,
so a picture that stops matching the code is a picture that stops being
produced.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

import pymopsmap as pm
from pymopsmap.psd import LognormalPSD
from pymopsmap.shapes import Sphere
from pymopsmap.species import Mode

FIGURES = Path(__file__).resolve().parents[2] / "docs" / "figures"
WAVELENGTHS = np.round(np.linspace(0.35, 2.0, 60), 4).tolist()


def spectra(path: Path) -> None:
    """What four catalogue species look like across the solar spectrum."""
    import matplotlib.pyplot as plt

    species = [
        (pm.CAMS.SULPHATE, "sulphate"),
        (pm.CAMS.SEA_SALT, "sea salt"),
        (pm.CAMS.DUST, "dust"),
        (pm.CAMS.BLACK_CARBON, "black carbon"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for entry, label in species:
        result = pm.load(entry).compute(wl=WAVELENGTHS, rh=50, quiet=True)
        normalised = result["kext"] / result["kext"].sel(
            wl=0.55, method="nearest"
        )
        axes[0].plot(result["wl"] * 1e3, normalised, label=label)
        axes[1].plot(result["wl"] * 1e3, result["ssa"], label=label)

    axes[0].set_ylabel("extinction, normalised at 550 nm")
    axes[0].set_yscale("log")
    axes[1].set_ylabel("single scattering albedo")
    axes[1].set_ylim(0, 1.02)
    for axis in axes:
        axis.set_xlabel("wavelength [nm]")
        axis.grid(alpha=0.3)
    axes[0].legend(fontsize=9)
    fig.suptitle("Four CAMS species at 50 % relative humidity", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")


def humidity(path: Path) -> None:
    """How far three species swell, and what it does to their scattering."""
    import matplotlib.pyplot as plt

    grid = [0.0, 30.0, 50.0, 70.0, 80.0, 90.0, 95.0]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for entry, label in [
        (pm.CAMS.SEA_SALT, "sea salt"),
        (pm.CAMS.SULPHATE, "sulphate"),
        (pm.CAMS.DUST, "dust"),
    ]:
        result = pm.load(entry).compute(wl=[0.55], rh=grid, quiet=True)
        series = result.squeeze("wl", drop=True)
        axes[0].plot(
            grid, series["kext"] / series["kext"].isel(rh=0), "o-", label=label
        )
        axes[1].plot(grid, series["ssa"], "o-", label=label)

    axes[0].set_ylabel("extinction, relative to dry")
    axes[1].set_ylabel("single scattering albedo")
    for axis in axes:
        axis.set_xlabel("relative humidity [%]")
        axis.grid(alpha=0.3)
    axes[0].legend(fontsize=9)
    fig.suptitle("Hygroscopic growth at 550 nm", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")


def declared_sweep(path: Path) -> None:
    """A species whose size and absorption are both axes of its own space."""
    import matplotlib.pyplot as plt

    radii = np.logspace(np.log10(0.03), np.log10(1.0), 24)
    absorption = np.logspace(-4, -1.3, 18)
    specie = pm.Specie.custom(
        Mode(
            shape=Sphere(),
            psd=LognormalPSD(
                rm=xr.DataArray(radii, dims="rm"),
                sigma=1.6,
                n=1e9,
                rmin=0.005,
                rmax=20.0,
            ),
            n_real=1.53,
            n_imag=xr.DataArray(absorption, dims="n_imag"),
            density_dry=1.8,
        ),
        name="parameter-space",
    )
    print(f"  declared sweep: {specie.swept}")
    result = specie.compute(wl=[0.55], quiet=True).squeeze("wl", drop=True)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for axis, name, label in (
        (axes[0], "ssa", "single scattering albedo"),
        (axes[1], "g", "asymmetry parameter"),
    ):
        mesh = axis.pcolormesh(
            result["rm"],
            result["n_imag"],
            result[name].transpose("n_imag", "rm"),
            shading="nearest",
        )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel("modal radius [$\\mu$m]")
        axis.set_ylabel("imaginary refractive index")
        fig.colorbar(mesh, ax=axis, label=label)
    fig.suptitle(
        f"One species, {len(radii)} x {len(absorption)} points, "
        "declared and swept in one call",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")


def angular(path: Path) -> None:
    """The phase function, across species and across humidity."""
    import matplotlib.pyplot as plt

    outputs = frozenset({pm.OutputType.PHASE_FUNCTION})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    for entry, label in [
        (pm.CAMS.SULPHATE, "sulphate"),
        (pm.CAMS.DUST, "dust"),
        (pm.CAMS.SEA_SALT, "sea salt"),
    ]:
        point = (
            pm.load(entry)
            .compute(
                wl=[0.55], rh=50, outputs=outputs, n_angles=721, quiet=True
            )
            .squeeze("wl", drop=True)
        )
        axes[0].semilogy(point["theta"], point["phase"], label=label)
    axes[0].set_title("three species at 50 % humidity", fontsize=10)
    axes[0].legend(fontsize=9)

    humidities = [0.0, 50.0, 80.0, 95.0]
    wet = pm.load(pm.CAMS.SEA_SALT).compute(
        wl=[0.55], rh=humidities, outputs=outputs, n_angles=721, quiet=True
    )
    for value in humidities:
        axes[1].semilogy(
            wet["theta"],
            wet["phase"].sel(rh=value).squeeze("wl", drop=True),
            label=f"RH = {value:.0f} %",
        )
    axes[1].set_title("sea salt, as it takes up water", fontsize=10)
    axes[1].legend(fontsize=9)

    for axis in axes:
        axis.set_xlim(0, 180)
        axis.set_xticks(np.arange(0, 181, 30))
        axis.set_xlabel("scattering angle [deg]")
        axis.set_ylabel("phase function")
        axis.grid(alpha=0.3)
    fig.suptitle("Phase functions at 550 nm", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")


def mixture(path: Path) -> None:
    """A mixture reweighted without recomputing anything."""
    import matplotlib.pyplot as plt

    fractions = np.linspace(0.0, 1.0, 11)
    fig, axis = plt.subplots(figsize=(7, 4.2))
    albedo, asymmetry = [], []
    for dust in fractions:
        mix = pm.Mix.from_optical_depth(
            {
                pm.CAMS.DUST: max(dust, 1e-6),
                pm.CAMS.SULPHATE: max(1.0 - dust, 1e-6),
            },
            wl_ref=0.55,
            rh_ref=50,
        )
        point = mix.compute(wl=[0.55], rh=50, quiet=True).squeeze(
            "wl", drop=True
        )
        albedo.append(float(point["ssa"]))
        asymmetry.append(float(point["g"]))

    axis.plot(fractions, albedo, "o-", label="single scattering albedo")
    axis.plot(fractions, asymmetry, "s-", label="asymmetry parameter")
    axis.set_xlabel("fraction of the optical depth carried by dust")
    axis.grid(alpha=0.3)
    axis.legend(fontsize=9)
    fig.suptitle(
        "Eleven mixtures, two MOPSMAP runs: the weights apply afterwards",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")


PANELS = {
    "guide-spectra.png": spectra,
    "guide-humidity.png": humidity,
    "guide-sweep.png": declared_sweep,
    "guide-angular.png": angular,
    "guide-mixture.png": mixture,
}


def main() -> None:
    for name, draw in PANELS.items():
        draw(FIGURES / name)


if __name__ == "__main__":
    main()
