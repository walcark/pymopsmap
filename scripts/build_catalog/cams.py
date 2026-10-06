"""
Build the CAMS species catalogue from the MAJA microphysical parameter table.

Run once, when the source table changes:

    pixi run -e dev python -m scripts.build_catalog.cams

Reads ``data/cams/`` and writes one file per CAMS version into
``src/pymopsmap/data/cams/``. Every convention the runtime used to fix at its
call sites is applied here instead, once:

  * wavelengths from nanometres to micrometres,
  * ``lnvar`` (the logarithm of the geometric standard deviation) to ``sigma``,
  * the imaginary refractive index to the positive MOPSMAP convention,
  * the truncation radii, previously hardcoded in the adapter, written out.

Dry densities are deliberately absent: CAMS publishes none, and inventing them
would freeze approximate values into what becomes the source of truth. MOPSMAP
therefore keeps applying its own default until a sourced table exists.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import xarray as xr

from pymopsmap.species import schema

ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIR = ROOT / "data" / "cams"
TARGET_DIR = ROOT / "src" / "pymopsmap" / "data" / "cams"

# Radius range covered by the modelled ensemble, in micrometres. Hardcoded in
# the adapter this catalogue replaces.
RMIN_UM = 0.001
RMAX_UM = 40.0

MODES = {"fine": "f", "coarse": "c"}

AMPLITUDE_NOTE = (
    "CAMS modal weights, normalised on the fine mode. They carry the ratio "
    "between modes; the absolute scale is set by the mixture."
)


def build() -> None:
    """Write one catalogue file per CAMS version."""
    source = xr.load_dataset(
        SOURCE_DIR / "cams_aer_microphysical_parameters.nc"
    )
    weights = json.loads(
        (SOURCE_DIR / "cams_aer_modes_concentrations.json").read_text()
    )

    TARGET_DIR.mkdir(parents=True, exist_ok=True)
    for version in source["cams_versions"].values:
        tree = _version_tree(source, weights, str(version))
        schema.write(tree, TARGET_DIR / f"{version}.nc")
        print(f"wrote {TARGET_DIR / f'{version}.nc'}")

    shutil.copy(
        SOURCE_DIR / "cams_aer_modes_concentrations.json",
        TARGET_DIR / "modes.json",
    )


def _version_tree(
    source: xr.Dataset, weights: dict[str, list[float]], version: str
) -> xr.DataTree:
    tree = xr.DataTree()
    tree.attrs.update(
        source="CAMS",
        version=version,
        description=source.attrs.get("description", ""),
        authors=source.attrs.get("authors", ""),
    )
    for specie in source["aerosols_species"].values:
        name = str(specie)
        selected = source.sel(aerosols_species=specie, cams_versions=version)
        if _absent_from_version(selected):
            # The source table pads every version with every specie and fills
            # the gaps with NaN. One file per version lets an absent specie
            # simply not be there, so loading it fails clearly instead of
            # yielding NaN.
            print(f"  {version}: skipping '{name}', absent from this version")
            continue
        tree[name] = xr.DataTree()
        tree[name].attrs["growth"] = schema.Growth.TABULATED.value
        for mode, suffix in MODES.items():
            amplitude = weights[name][0 if suffix == "f" else 1]
            tree[f"{name}/{mode}"] = xr.DataTree(
                _mode_dataset(selected, suffix, amplitude, f"{version}/{name}")
            )
    return tree


def _absent_from_version(selected: xr.Dataset) -> bool:
    """True when the source carries no data at all for this specie."""
    return bool(selected["mr_f"].isnull().all())


def _mode_dataset(
    selected: xr.Dataset, suffix: str, amplitude: float, where: str
) -> xr.Dataset:
    """One mode of one specie, in canonical units and conventions."""
    wl_um = selected["wavelength"].values.astype(float) * 1e-3
    rh = selected["relative_humidity"].values.astype(float)
    n_imag = _positive_imaginary_index(selected[f"mi_{suffix}"], where)

    ds = xr.Dataset(
        data_vars={
            "n_real": (("rh", "wl"), selected[f"mr_{suffix}"].values),
            "n_imag": (("rh", "wl"), n_imag),
            "rm": (("rh",), selected[f"rmodal_{suffix}"].values),
            # The source stores ln(sigma) under the name lnvar.
            "sigma": (("rh",), np.exp(selected[f"lnvar_{suffix}"].values)),
            "rmin": RMIN_UM,
            "rmax": RMAX_UM,
            "n": amplitude,
        },
        coords={"rh": rh, "wl": wl_um},
        attrs={"psd_type": "lognormal", "shape_type": "sphere"},
    )

    ds["wl"].attrs.update(units="um", long_name="wavelength")
    ds["rh"].attrs.update(units="percent", long_name="relative humidity")
    ds["n_real"].attrs.update(units="1")
    ds["n_imag"].attrs.update(units="1", sign="positive")
    ds["rm"].attrs.update(units="um", long_name="modal radius")
    ds["sigma"].attrs.update(
        units="1", long_name="geometric standard deviation"
    )
    ds["rmin"].attrs.update(units="um")
    ds["rmax"].attrs.update(units="um")
    ds["n"].attrs.update(units="m-3", note=AMPLITUDE_NOTE)
    return ds


def _positive_imaginary_index(source: xr.DataArray, where: str) -> np.ndarray:
    """
    Convert the imaginary index to the positive MOPSMAP convention.

    The source stores it negative, so the conversion is a sign flip and it
    lives here rather than at every call site. A handful of source points carry
    the opposite sign, which is a typo rather than a physical value: the
    magnitude is kept and the anomaly is reported, so a new one in a future
    source table cannot pass unnoticed.
    """
    values = source.values
    wrong_sign = values > 0
    if wrong_sign.any():
        for index in np.argwhere(wrong_sign):
            rh = source["relative_humidity"].values[index[0]]
            wl = source["wavelength"].values[index[1]]
            print(
                f"  {where}: source sign anomaly on {source.name} "
                f"at rh={rh}, wl={wl} nm: {values[tuple(index)]:.6g}. "
                "Magnitude kept, sign normalised."
            )
    return np.abs(values)


if __name__ == "__main__":
    build()
