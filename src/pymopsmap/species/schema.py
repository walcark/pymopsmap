"""The NetCDF contract of an aerosol species.

A species file is an ``xarray.DataTree``: the root carries the source and its
version, each child group is a species, and each grandchild group is one mode
of that species.

    /                    attrs: source, version, schema_rev
    +-- /sulphate        attrs: growth
        +-- /fine        attrs: psd_type, shape_type
        |     n_real, n_imag, the size distribution parameters including its
        |     amplitude, the shape parameters, density_dry, and kappa when
        |     relevant
        +-- /coarse

A species file is self-contained: it holds everything needed to compute the
species as it stands. A mixture then scales a whole species by one factor,
which preserves the ratio between its modes.

Conventions are declared, never dispatched on. Units are the one exception:
they are declared and honoured, because a finite conversion table is cheap and
lets any NetCDF reader make sense of the file. The *meaning* of a quantity is
carried by its variable name, so ``sigma`` and ``ln_sigma`` would be different
variables rather than one variable with a switch. Signs are normalised at write
time and asserted at read time, which catches a badly written file where an
attribute alone would not.
"""

from __future__ import annotations

from enum import StrEnum
from pathlib import Path
from typing import get_args

import xarray as xr
from pydantic import BaseModel

from pymopsmap.exceptions import SchemaError
from pymopsmap.psd import PSD
from pymopsmap.shapes import Shape

# Bumped on any breaking layout change. It feeds the sweep store version key.
SCHEMA_REV = 1

# Length units accepted on read, expressed in micrometres.
LENGTH_UNITS: dict[str, float] = {"um": 1.0, "nm": 1e-3, "m": 1e6}

LENGTH_VARIABLES = frozenset({"wl", "rm", "rmin", "rmax", "radius", "radii"})

CANONICAL_UNITS: dict[str, str] = {
    "wl": "um",
    "rh": "percent",
    "rm": "um",
    "rmin": "um",
    "rmax": "um",
    "radius": "um",
    "n_real": "1",
    "n_imag": "1",
    "sigma": "1",
    "kappa": "1",
    "density_dry": "g cm-3",
    "n": "m-3",
    "A": "m-3",
}

HUMIDITY_DIM = "rh"


class Growth(StrEnum):
    """How a species responds to relative humidity."""

    TABULATED = "tabulated"
    KAPPA = "kappa"
    NONE = "none"


# Size distributions MOPSMAP refuses to grow hygroscopically.
GROWTH_INCOMPATIBLE_PSD = frozenset({"mod-gamma"})


# ---------------------------------------------------------------------------
# Discriminated unions as the source of truth
# ---------------------------------------------------------------------------


def _by_discriminator(union) -> dict[str, type[BaseModel]]:
    """Map each ``type`` literal of an annotated union to its model class.

    This allows to populate the PSF_MODELS and SHAPE_MODELS list based on
    their defined child classes.
    """
    members = get_args(get_args(union)[0])
    return {get_args(m.model_fields["type"].annotation)[0]: m for m in members}


PSD_MODELS: dict[str, type[BaseModel]] = _by_discriminator(PSD)
SHAPE_MODELS: dict[str, type[BaseModel]] = _by_discriminator(Shape)


# The variable carrying the amplitude of each size distribution. Scaling a
# species multiplies it, which preserves the ratio between its modes. None
# means the amplitude lives in an external file.
AMPLITUDE_FIELD: dict[str, str | None] = {
    "fixed": "n",
    "lognormal": "n",
    "mod-gamma": "A",
    "distr-list": "concentrations",
    "bin-file": None,
}


def expected_psd_variables(psd_type: str) -> frozenset[str]:
    """Variables a mode must carry for the given size distribution type."""
    model = PSD_MODELS[psd_type]
    return frozenset(model.model_fields) - {"type"}


def expected_shape_variables(shape_type: str) -> frozenset[str]:
    """Variables a mode must carry for the given shape type."""
    model = SHAPE_MODELS[shape_type]
    return frozenset(model.model_fields) - {"type"}


# ---------------------------------------------------------------------------
# Reading and writing
# ---------------------------------------------------------------------------


def read(path: str | Path) -> xr.DataTree:
    """
    Open a species file, convert its units, and validate it.

    Parameters
    ----------
    path : str or Path
        The NetCDF file to read.

    Returns
    -------
    xr.DataTree
        The tree, with every length converted to micrometres.

    Raises
    ------
    SchemaError
        If the file does not follow the canonical schema.
    """
    tree = xr.open_datatree(path).load()
    tree = _convert_units(tree)
    validate(tree)
    return tree


def write(
    tree: xr.DataTree, path: str | Path, check_units: bool = True
) -> None:
    """
    Validate a tree and write it, stamped with the schema revision.

    Parameters
    ----------
    tree : xr.DataTree
        The species tree to write.
    path : str or Path
        Destination NetCDF file.
    check_units : bool
        Validate before writing. Set to False only to write a file that
        deliberately uses a non-canonical unit, for conversion tests.
    """
    if check_units:
        validate(tree)
    tree = tree.copy()
    tree.attrs["schema_rev"] = SCHEMA_REV
    tree.to_netcdf(path)


def _convert_units(tree: xr.DataTree) -> xr.DataTree:
    """Bring every declared length onto micrometres."""
    for node in tree.subtree:
        if not node.has_data:
            continue
        ds = node.to_dataset()
        for name in list(ds.variables):
            if name not in LENGTH_VARIABLES:
                continue
            unit = ds[name].attrs.get("units", "um")
            if unit == "um":
                continue
            if unit not in LENGTH_UNITS:
                raise SchemaError(
                    f"{node.path}: variable '{name}' declares unknown length "
                    f"unit '{unit}'. Accepted: {sorted(LENGTH_UNITS)}."
                )
            converted = ds[name] * LENGTH_UNITS[unit]
            converted.attrs = dict(ds[name].attrs) | {"units": "um"}
            ds[name] = converted
        node.dataset = ds
    return tree


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def validate(tree: xr.DataTree) -> None:
    """
    Check a species tree against the canonical schema.

    Raises
    ------
    SchemaError
        On the first inconsistency found, naming the offending group.
    """
    for specie in tree.children.values():
        growth = _read_growth(specie)
        modes = [n for n in specie.children.values() if n.has_data]
        if not modes:
            raise SchemaError(f"{specie.path}: species has no mode group.")
        for mode in modes:
            _validate_mode(mode, growth)


def _read_growth(specie: xr.DataTree) -> Growth:
    raw = specie.attrs.get("growth")
    if raw is None:
        raise SchemaError(f"{specie.path}: missing 'growth' attribute.")
    try:
        return Growth(raw)
    except ValueError:
        raise SchemaError(
            f"{specie.path}: unknown growth '{raw}'. "
            f"Expected one of {[g.value for g in Growth]}."
        ) from None


def _validate_mode(mode: xr.DataTree, growth: Growth) -> None:
    ds = mode.to_dataset()
    where = mode.path

    psd_type = _read_discriminator(ds, where, "psd_type", PSD_MODELS)
    shape_type = _read_discriminator(ds, where, "shape_type", SHAPE_MODELS)

    _check_growth_structure(ds, where, growth, psd_type)
    _check_units(ds, where)
    _check_required_variables(
        ds,
        where,
        expected_psd_variables(psd_type)
        | expected_shape_variables(shape_type),
    )
    _check_conventions(ds, where)


def _read_discriminator(
    ds: xr.Dataset, where: str, attr: str, models: dict[str, type]
) -> str:
    value = ds.attrs.get(attr)
    if value is None:
        raise SchemaError(f"{where}: missing '{attr}' attribute.")
    if value not in models:
        raise SchemaError(
            f"{where}: unknown {attr} '{value}'. "
            f"Expected one of {sorted(models)}."
        )
    return value


def _check_growth_structure(
    ds: xr.Dataset, where: str, growth: Growth, psd_type: str
) -> None:
    has_rh = HUMIDITY_DIM in ds.dims
    has_kappa = "kappa" in ds.variables

    if growth is Growth.TABULATED and not has_rh:
        raise SchemaError(
            f"{where}: growth='tabulated' requires an '{HUMIDITY_DIM}' "
            "dimension carrying the tabulated wet state."
        )
    if growth is Growth.KAPPA:
        if has_rh:
            raise SchemaError(
                f"{where}: growth='kappa' forbids an '{HUMIDITY_DIM}' "
                "dimension. The file would then hold two sources of truth."
            )
        if not has_kappa:
            raise SchemaError(
                f"{where}: growth='kappa' requires a 'kappa' variable."
            )
        if psd_type in GROWTH_INCOMPATIBLE_PSD:
            raise SchemaError(
                f"{where}: growth='kappa' cannot be combined with "
                f"psd_type='{psd_type}'; MOPSMAP does not implement it."
            )
    if growth is Growth.NONE and (has_rh or has_kappa):
        raise SchemaError(
            f"{where}: growth='none' forbids both an '{HUMIDITY_DIM}' "
            "dimension and a 'kappa' variable."
        )


def _check_units(ds: xr.Dataset, where: str) -> None:
    for name, canonical in CANONICAL_UNITS.items():
        if name not in ds.variables:
            continue
        declared = ds[name].attrs.get("units")
        if declared is None or declared == canonical:
            continue
        if name in LENGTH_VARIABLES and declared in LENGTH_UNITS:
            continue
        raise SchemaError(
            f"{where}: variable '{name}' declares units '{declared}', "
            f"expected '{canonical}'."
        )


def _check_required_variables(
    ds: xr.Dataset, where: str, required: frozenset[str]
) -> None:
    missing = sorted(required - set(ds.variables))
    if missing:
        raise SchemaError(f"{where}: missing required variable(s) {missing}.")
    for name in ("n_real", "n_imag"):
        if name not in ds.variables:
            raise SchemaError(f"{where}: missing required variable '{name}'.")


def _check_conventions(ds: xr.Dataset, where: str) -> None:
    if (ds["n_imag"] < 0).any():
        raise SchemaError(
            f"{where}: 'n_imag' holds negative values. The canonical "
            "convention is positive, matching what MOPSMAP expects."
        )
