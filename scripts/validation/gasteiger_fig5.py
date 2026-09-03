"""
Reproduce Figure 5 of Gasteiger and Wiegner (2018), GMD 11:2739.

    pixi run -e dev python -m scripts.validation.gasteiger_fig5

"Properties of OPAC aerosol types as a function of relative humidity RH
calculated with the kappa parameterization (Zieger et al., 2013) implemented in
MOPSMAP": the ten OPAC climatologies, at three lidar wavelengths, over five
humidities, in four panels.

The figure is a validation target rather than a plot: what matters is whether
the pipeline lands on the published curves. Types whose components fall outside
the optical dataset at hand are reported and skipped, so a partial dataset
still produces a partial figure.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

import pymopsmap as pm
from pymopsmap.exceptions import CoverageError, DomainError
from pymopsmap.species import OPAC_COMPOSITIONS

WAVELENGTHS_UM = [0.355, 0.532, 1.064]
HUMIDITIES = [0.0, 50.0, 70.0, 80.0, 90.0]

# The four rows of the figure. The extinction is shown normalised on its dry
# value, which is what makes the ten types comparable on one panel.
PANELS = {
    "kext": "extinction coefficient, normalised on RH = 0",
    "ssa": "single scattering albedo",
    "ext_to_mass": "extinction to mass conversion factor",
    "back_to_mass": "mass to backscatter conversion factor",
}

OUTPUTS = frozenset({pm.OutputType.INTEGRATED, pm.OutputType.LIDAR})


def compute() -> xr.Dataset:
    """
    Compute every OPAC type over the humidity grid of the figure.

    Returns
    -------
    xr.Dataset
        Dimensions ``(type, rh, wl)``, with one variable per panel. Types the
        dataset cannot cover are absent.
    """
    per_type: dict[str, xr.Dataset] = {}
    for name in OPAC_COMPOSITIONS:
        try:
            result = pm.opac_mix(name).compute(
                wl=WAVELENGTHS_UM,
                rh=HUMIDITIES,
                outputs=OUTPUTS,
                quiet=True,
            )
        except (CoverageError, DomainError) as exc:
            print(f"  {name}: skipped, {exc}")
            continue
        per_type[name] = _panels(result)
        print(f"  {name}: done")

    if not per_type:
        raise SystemExit(
            "No OPAC type could be computed with the optical dataset at hand."
        )
    return xr.concat(
        list(per_type.values()),
        dim=xr.DataArray(list(per_type), dims="type", name="type"),
    )


def _panels(result: xr.Dataset) -> xr.Dataset:
    """Keep the four plotted quantities, extinction normalised on dry."""
    panels = result[list(PANELS)]
    return panels.assign(kext=panels["kext"] / panels["kext"].sel(rh=0.0))


def report(data: xr.Dataset) -> None:
    """Print the computed values, one block per panel."""
    for name, label in PANELS.items():
        print(f"\n=== {label} ===")
        table = data[name].sel(wl=0.532, method="nearest").to_pandas()
        print(table.round(4).to_string())


def plot(data: xr.Dataset, path: Path) -> None:
    """Draw the four panels, one line per aerosol type."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(
        len(PANELS), len(WAVELENGTHS_UM), figsize=(13, 12), sharex=True
    )
    for row, (name, label) in enumerate(PANELS.items()):
        for column, wavelength in enumerate(WAVELENGTHS_UM):
            axis = axes[row, column]
            for aerosol in data["type"].values:
                series = (
                    data[name]
                    .sel(type=aerosol)
                    .sel(wl=wavelength, method="nearest")
                )
                axis.plot(data["rh"], series, marker="o", label=str(aerosol))
            if row == 0:
                axis.set_title(f"{wavelength * 1e3:.0f} nm")
            if column == 0:
                axis.set_ylabel(label, fontsize=8)
            if row == len(PANELS) - 1:
                axis.set_xlabel("relative humidity [%]")
            axis.grid(alpha=0.3)
    axes[0, -1].legend(fontsize=6, loc="upper left")
    fig.suptitle(
        "Gasteiger and Wiegner (2018), Figure 5, recomputed", fontsize=11
    )
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    print(f"\nwrote {path}")


def main() -> None:
    print("Computing the OPAC types:")
    data = compute()
    report(data)
    if np.isfinite(data["kext"]).any():
        plot(data, Path("gasteiger_fig5.png"))


if __name__ == "__main__":
    main()
