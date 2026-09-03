"""
Build the OPAC species catalogue, kappa flavour.

Run once:

    pixi run -e dev python -m scripts.build_catalog.opac

Two sources, neither of which needs the network:

  * the dry size distributions and densities of Hess et al. (1998) Table 1c,
    transcribed below and cross-checked against the OPAC desert example of the
    MOPSMAP user guide, which pins waso, minm, miam and micm independently,
  * the refractive indices shipped with MOPSMAP in bin/mopsmap/data, on a
    common 61-point grid from 0.25 to 40 micrometres.

Only the dry state is written, with a hygroscopicity parameter, so MOPSMAP
applies the growth itself. The wet tabulated state is the other flavour, and it
needs the per-humidity GEISA files, whose host was decommissioned during the
migration of the database from IPSL to SEDOO: the pages moved, the file links
did not.

Every species is written as a sphere, which is what OPAC defines. The MOPSMAP
guide runs its mineral modes as spheroids instead, which a user expresses with
a custom species rather than a catalogue entry.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

from pymopsmap.species import schema

ROOT = Path(__file__).resolve().parents[2]
REFRACTIVE_DIR = ROOT / "bin" / "mopsmap" / "data"
TARGET_DIR = ROOT / "src" / "pymopsmap" / "data" / "opac"

# Hess, Koepke and Schult (1998), BAMS 79:831, Table 1c, as published on
# https://geisa.aeris-data.fr/opac/ :
#   sigma, modal radius [um], rmin [um], rmax [um], dry density [g cm-3]
DRY_STATE: dict[str, tuple[float, float, float, float, float]] = {
    "inso": (2.51, 0.471, 0.005, 20.0, 2.0),
    "waso": (2.24, 0.0212, 0.005, 20.0, 1.8),
    "soot": (2.00, 0.0118, 0.005, 20.0, 1.0),
    "ssam": (2.03, 0.209, 0.005, 20.0, 2.2),
    "sscm": (2.03, 1.75, 0.005, 60.0, 2.2),
    "minm": (1.95, 0.07, 0.005, 20.0, 2.6),
    "miam": (2.00, 0.39, 0.005, 20.0, 2.6),
    "micm": (2.15, 1.90, 0.005, 60.0, 2.6),
    "mitr": (2.20, 0.50, 0.02, 5.0, 2.5),
    "suso": (2.03, 0.0695, 0.005, 20.0, 1.7),
}

# Dry refractive index file shipped with MOPSMAP for each species. The sea salt
# modes share one index, as do the four mineral modes: what distinguishes them
# is their size distribution.
REFRACTIVE_FILE: dict[str, str] = {
    "inso": "refr_insoluble",
    "waso": "refr_waso00",
    "soot": "refr_soot",
    "ssam": "refr_seasalt00",
    "sscm": "refr_seasalt00",
    "minm": "refr_mineral",
    "miam": "refr_mineral",
    "micm": "refr_mineral",
    "mitr": "refr_mineral",
    "suso": "refr_suso00",
}

# Zieger et al. (2013), ACP 13:10609, Table 4, from the modified g(RH). Zero
# for the species that take up no water: MOPSMAP demands a value on every mode
# as soon as a humidity is set.
KAPPA: dict[str, float] = {
    "waso": 0.249,
    "ssam": 0.916,
    "sscm": 0.916,
    "suso": 0.804,
    "inso": 0.0,
    "soot": 0.0,
    "minm": 0.0,
    "miam": 0.0,
    "micm": 0.0,
    "mitr": 0.0,
}


def build() -> None:
    """Write the kappa flavour of the OPAC catalogue."""
    tree = xr.DataTree()
    tree.attrs.update(
        source="OPAC",
        version="kappa",
        description=(
            "OPAC aerosol components in their dry state, with the "
            "hygroscopicity parameter driving growth inside MOPSMAP."
        ),
        references=(
            "Hess, Koepke and Schult (1998), BAMS 79:831, Table 1c; "
            "Zieger et al. (2013), ACP 13:10609, Table 4; "
            "refractive indices from the MOPSMAP data directory."
        ),
    )
    for name in sorted(DRY_STATE):
        tree[name] = xr.DataTree()
        tree[name].attrs["growth"] = schema.Growth.KAPPA.value
        tree[f"{name}/only"] = xr.DataTree(_mode_dataset(name))

    TARGET_DIR.mkdir(parents=True, exist_ok=True)
    target = TARGET_DIR / "kappa.nc"
    schema.write(tree, target)
    print(f"wrote {target}")


def _mode_dataset(name: str) -> xr.Dataset:
    """The single mode of one OPAC species, dry."""
    sigma, rm, rmin, rmax, density = DRY_STATE[name]
    wl, n_real, n_imag = _read_refractive_index(REFRACTIVE_FILE[name])

    ds = xr.Dataset(
        data_vars={
            "n_real": (("wl",), n_real),
            "n_imag": (("wl",), n_imag),
            "rm": rm,
            "sigma": sigma,
            "rmin": rmin,
            "rmax": rmax,
            # Monomodal, so the amplitude is a bare reference: the mixture
            # supplies the concentration.
            "n": 1.0,
            "kappa": KAPPA[name],
            "density_dry": density,
        },
        coords={"wl": wl},
        attrs={"psd_type": "lognormal", "shape_type": "sphere"},
    )

    ds["wl"].attrs.update(units="um", long_name="wavelength")
    ds["n_real"].attrs.update(units="1")
    ds["n_imag"].attrs.update(units="1", sign="positive")
    ds["rm"].attrs.update(units="um", long_name="modal radius")
    ds["sigma"].attrs.update(
        units="1", long_name="geometric standard deviation"
    )
    ds["rmin"].attrs.update(units="um")
    ds["rmax"].attrs.update(units="um")
    ds["n"].attrs.update(
        units="m-3", note="reference amplitude; the mixture sets the scale"
    )
    ds["kappa"].attrs.update(units="1", long_name="hygroscopicity parameter")
    ds["density_dry"].attrs.update(units="g cm-3")
    return ds


def _read_refractive_index(
    filename: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a MOPSMAP refractive index file: wavelength, real, imaginary."""
    table = np.loadtxt(REFRACTIVE_DIR / filename, comments="#")
    return table[:, 0], table[:, 1], table[:, 2]


if __name__ == "__main__":
    build()
