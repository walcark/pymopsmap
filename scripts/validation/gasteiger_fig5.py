"""
Reproduce Figure 5 of Gasteiger and Wiegner (2018), GMD 11:2739.

    PYMOPSMAP_DATASET_SOURCE=... \
        pixi run -e dev python -m scripts.validation.gasteiger_fig5

"Properties of OPAC aerosol types as a function of relative humidity RH
calculated with the kappa parameterization (Zieger et al., 2013) implemented in
MOPSMAP": the ten OPAC climatologies, at three lidar wavelengths, over five
humidities, in four rows.

The layout follows the published one so the two can be laid side by side: same
panel grid, same axis limits and ticks, same colour per type, and the legend
split over the three panels of the top row the way the article splits it.

The five types containing soot need the extended archive of the optical data
set, soot reaching a real refractive index of 1.75 where the main archive stops
at 1.64. Without it they are reported and skipped, so a partial data set still
produces a partial figure.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

import pymopsmap as pm
from pymopsmap.exceptions import CoverageError, DomainError
from pymopsmap.species import OPAC_COMPOSITIONS

FIGURES = Path(__file__).resolve().parents[2] / "docs" / "figures"

WAVELENGTHS_UM = [0.355, 0.532, 1.064]
HUMIDITIES = [0.0, 50.0, 70.0, 80.0, 90.0]

# The wavelength every extinction coefficient is normalised on: "the extinction
# coefficient normalized to the extinction coefficient of the same aerosol type
# at RH = 0 % and lambda = 532 nm". One reference for all three columns, which
# is why the 1064 nm panel does not start at one.
REFERENCE_WL_UM = 0.532

OUTPUTS = frozenset({pm.OutputType.INTEGRATED, pm.OutputType.LIDAR})

# The four rows: variable, axis label, limits, tick step.
ROWS = [
    ("kext", "Normalized extinction", (0.0, 5.0), 1.0),
    ("ssa", "Single-scattering albedo $\\omega_0$", (0.55, 1.005), 0.1),
    (
        "ext_to_mass",
        "Extinction to mass conversion\nfactor $\\eta$ [$g\\,m^{-2}$]",
        (0.0, 2.5),
        0.5,
    ),
    (
        "back_to_mass",
        "Mass to backscatter conver-\nsion factor $Z$"
        " [$m^2\\,sr^{-1}\\,g^{-1}$]",
        (0.0, 0.12),
        0.02,
    ),
]

# The colours of the published legend, and which top-row panel carries each
# group of entries.
TYPES = {
    "urban": ("Urban", "black", 0),
    "desert": ("Desert", "orange", 0),
    "arctic": ("Arctic", "gold", 0),
    "antarctic": ("Antarctic", "purple", 0),
    "continental_clean": ("Contin. clean", "lightgreen", 1),
    "continental_average": ("Contin. average", "forestgreen", 1),
    "continental_polluted": ("Contin. polluted", "darkgreen", 1),
    "maritime_clean": ("Maritime clean", "lightskyblue", 2),
    "maritime_polluted": ("Maritime polluted", "blue", 2),
    "maritime_tropical": ("Maritime tropical", "navy", 2),
}


def compute() -> xr.Dataset:
    """
    Compute every OPAC type over the humidity grid of the figure.

    Returns
    -------
    xr.Dataset
        Dimensions ``(type, rh, wl)``, with one variable per row. Types the
        data set cannot cover are absent.
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
        per_type[name] = _rows(result)
        print(f"  {name}: done")

    if not per_type:
        raise SystemExit(
            "No OPAC type could be computed with the optical dataset at hand."
        )
    return xr.concat(
        list(per_type.values()),
        dim=xr.DataArray(list(per_type), dims="type", name="type"),
    )


def _rows(result: xr.Dataset) -> xr.Dataset:
    """Keep the four plotted quantities, extinction normalised."""
    rows = result[[name for name, *_ in ROWS]]
    reference = (
        rows["kext"].sel(rh=0.0).sel(wl=REFERENCE_WL_UM, method="nearest")
    )
    return rows.assign(kext=rows["kext"] / reference)


def report(data: xr.Dataset) -> None:
    """Print the computed values, one block per row."""
    for name, label, *_ in ROWS:
        print(f"\n=== {label.replace(chr(10), ' ')}, at 532 nm ===")
        table = data[name].sel(wl=0.532, method="nearest").to_pandas()
        print(table.round(4).to_string())
    _check_statements(data)


def _check_statements(data: xr.Dataset) -> None:
    """
    Compare against what section 5.1 says in words.

    The figure itself is a scatter of markers whose values are not tabulated,
    so these statements are the only numerical hold the article gives.
    """
    print("\n=== statements of section 5.1 ===")
    eta = data["ext_to_mass"].sel(wl=0.532, method="nearest").sel(rh=0.0)
    highest = str(eta.idxmax("type").values)
    print(f"highest eta at RH = 0: {highest} ({float(eta.max()):.3f})")
    print("  article: 'the desert aerosol type ... shows the highest values'")

    if "desert" in data["type"]:
        span = data["ssa"].sel(wl=0.532, method="nearest").sel(type="desert")
        print(
            f"desert omega_0 over the humidity grid:"
            f" {float(span.min()):.4f} to {float(span.max()):.4f}"
        )
        print("  article: 'virtually independent on the RH'")


def plot(data: xr.Dataset, path: Path) -> None:
    """Draw the figure in the layout of the published one."""
    import matplotlib.pyplot as plt

    present = [name for name in TYPES if name in set(data["type"].values)]

    fig, axes = plt.subplots(
        len(ROWS), len(WAVELENGTHS_UM), figsize=(9.5, 11), sharex=True
    )
    for row, (name, label, limits, step) in enumerate(ROWS):
        for column, wavelength in enumerate(WAVELENGTHS_UM):
            axis = axes[row, column]
            for aerosol in present:
                _, colour, _ = TYPES[aerosol]
                axis.plot(
                    data["rh"],
                    data[name]
                    .sel(type=aerosol)
                    .sel(wl=wavelength, method="nearest"),
                    marker="o",
                    markersize=4,
                    linewidth=1.2,
                    color=colour,
                    label=TYPES[aerosol][0],
                )
            axis.set_ylim(*limits)
            axis.set_yticks(
                np.round(np.arange(limits[0], limits[1] + step / 2, step), 3)
            )
            axis.set_xlim(-3.0, 93.0)
            axis.set_xticks(HUMIDITIES)
            axis.grid(alpha=0.3, linewidth=0.5)
            axis.tick_params(labelsize=8)
            if row == 0:
                axis.set_title(
                    f"$\\lambda$ = {wavelength * 1e3:.0f} nm", fontsize=10
                )
            if column == 0:
                axis.set_ylabel(label, fontsize=8)
            else:
                axis.tick_params(labelleft=False)
            if row == len(ROWS) - 1:
                axis.set_xlabel("Relative humidity $RH$ [%]", fontsize=8)

    _legends(axes[0], present)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    print(f"\nwrote {path}")


def _legends(top_row, present: list[str]) -> None:
    """Split the legend over the three top panels, as the article does."""
    from matplotlib.lines import Line2D

    for column, axis in enumerate(top_row):
        entries = [
            Line2D(
                [],
                [],
                color=TYPES[name][1],
                marker="o",
                markersize=4,
                linewidth=1.2,
                label=TYPES[name][0],
            )
            for name in present
            if TYPES[name][2] == column
        ]
        if entries:
            axis.legend(handles=entries, fontsize=7, loc="upper left")


def main() -> None:
    print("Computing the OPAC types:")
    data = compute()
    report(data)
    if np.isfinite(data["kext"]).any():
        plot(data, FIGURES / "gasteiger_fig5.png")


if __name__ == "__main__":
    main()
