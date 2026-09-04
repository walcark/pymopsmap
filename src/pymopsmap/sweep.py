"""
The parameter space of a computation, and the engine that walks it.

Space normalisation turns the arguments of ``compute`` into the points to
visit; xsweep then runs them, stores each result, and assembles them back onto
the dimensions that were asked for. Re-running a grid already computed calls
MOPSMAP for nothing, and an interrupted sweep resumes from what it is missing.
"""

from __future__ import annotations

import hashlib
import os
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

from pymopsmap.utils import CACHE_DIR


def distinct_values(value: Any) -> list[Any]:
    """
    The scalar values a parameter takes, whatever shape it was given in.

    A field repeats itself between pixels, so the distinct values are what
    matters to anyone asking what a request will actually compute.
    """
    if value is None:
        return [None]
    if isinstance(value, xr.DataArray):
        return [_scalar(v) for v in np.unique(np.asarray(value).ravel())]
    if _is_sequence(value):
        return [_scalar(v) for v in np.unique(np.asarray(value).ravel())]
    return [_scalar(value)]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, (list, tuple, np.ndarray))


def _scalar(value: Any) -> Any:
    return value.item() if isinstance(value, np.generic) else value


# Where the sweep store lives. One directory per store, keyed on what makes
# two computations different, so unrelated species never share an entry.
SWEEP_STORE_ENV = "PYMOPSMAP_SWEEP_STORE"


def store_root() -> Path:
    """Directory holding the sweep stores."""
    override = os.environ.get(SWEEP_STORE_ENV)
    if override:
        return Path(override)
    return CACHE_DIR / "sweeps"


def run_sweep(
    point: Callable[..., xr.Dataset],
    space: xr.Dataset,
    outputs: Iterable[str],
    version: str,
    fixed: list[str] | None = None,
    quiet: bool = False,
) -> xr.Dataset:
    """
    Walk a parameter space, one MOPSMAP run per point.

    Parameters
    ----------
    point : callable
        Computes one point. It receives the swept variables by name plus
        ``wl``, and returns a dataset over the wavelengths.
    space : xr.Dataset
        The points to visit. ``wl`` is a vector consumed whole, since MOPSMAP
        takes the entire spectral grid in one run; every other variable loops.
    fixed : list of str, optional
        Parameters of the space held constant, declared ``const`` so they add
        no dimension to the result.
    outputs : iterable of str
        The variables the point function produces.
    version : str
        Everything that makes this computation different from another one:
        the species, its source version, the schema revision, the requested
        outputs. It enters the store key, so two species cannot collide.
    quiet : bool
        Unused for now; kept so callers do not special-case it.

    Returns
    -------
    xr.Dataset
        One variable per output, carrying a dimension per swept axis.
    """
    import xsweep

    fixed = fixed or []
    # Every variable of the space is a swept parameter but the fixed ones;
    # a one-dimensional axis named after itself is promoted to a coordinate
    # by xarray, so both places have to be read.
    looped = sorted(
        {str(n) for n in space.variables if str(n) != "wl"} - set(fixed)
    )
    contract = _contract(looped, sorted(fixed), outputs)
    sweeper = xsweep.Sweeper(
        contract,
        point,
        xsweep.SweepPolicy(
            store=str(_store_path(version, space)),
            # A field repeats its values between pixels, and a MOPSMAP run
            # costs seconds: compute each distinct point once.
            dedup=True,
            # A failed point would otherwise become NaN, which is the silent
            # gap this pipeline already refuses at the MOPSMAP level.
            on_error="raise",
        ),
        version=version,
    )
    try:
        result = sweeper(space)
    except xsweep.PointFailed as failure:
        # The cause is what the caller can act on; the wrapper is bookkeeping.
        raise failure.__cause__ or failure from None
    # The status variable belongs to the store, not to the optical properties.
    return result.drop_vars("status", errors="ignore")


def _contract(
    looped: list[str], fixed: list[str], outputs: Iterable[str]
) -> str:
    """
    Build the call contract.

    Wavelength is a vector because MOPSMAP computes the whole grid in one run;
    everything else is a loop. Output dimensions beyond ``wl`` are discovered
    from the first call rather than declared, since the number of angles is a
    run parameter.
    """
    clauses = " ".join(
        filter(
            None,
            [
                f"loop({', '.join(looped)})" if looped else "",
                "vec(wl)",
                f"const({', '.join(fixed)})" if fixed else "",
            ],
        )
    )
    produced = ", ".join(f"{name}(wl)" for name in outputs)
    return f"{clauses} -> {produced}"


def _store_path(version: str, space: xr.Dataset) -> Path:
    """
    Where the results of one sweep live.

    A store holds one grid: its axes are fixed when it is created, so asking
    for a different grid is a different store rather than an extension of the
    first. The path therefore carries a fingerprint of the axes alongside the
    version, and resuming means restarting the same request, not widening it.
    """
    axes = "|".join(
        f"{name}={np.asarray(values.values).tobytes().hex()}"
        for name, values in sorted(space.coords.items())
    )
    digest = hashlib.blake2b(axes.encode(), digest_size=8).hexdigest()
    return store_root() / f"{_slug(version)}-{digest}"


def _slug(version: str) -> str:
    """A directory name that survives a version string."""
    return "".join(c if c.isalnum() or c in "-_." else "-" for c in version)


def build_space(wl: list[float], **axes: Any) -> tuple[xr.Dataset, list[str]]:
    """
    Turn the arguments of a computation into the space to sweep.

    Parameters
    ----------
    wl : list of float
        Wavelengths, always a vector: MOPSMAP takes the whole grid in one run.
    **axes : Any
        One entry per parameter. A scalar fixes it and adds no dimension. A
        sequence sweeps it under a dimension named after the parameter. A
        ``DataArray`` sweeps it under its own dimensions, which is how a field
        is computed: one value per pixel of an image, keeping its shape.

    Returns
    -------
    space : xr.Dataset
        Every parameter, swept ones carrying their own dimension and fixed
        ones stored without dimensions.
    fixed : list of str
        Names of the parameters held constant, which the contract declares
        as ``const`` so they add no dimension to the result.
    """
    variables: dict[str, Any] = {}
    fixed: list[str] = []
    for name, value in axes.items():
        if value is None:
            continue
        if isinstance(value, xr.DataArray):
            # The caller named the dimensions; keep them, whatever their rank.
            variables[name] = value
        elif _is_sequence(value):
            variables[name] = (name, np.asarray([_scalar(v) for v in value]))
        else:
            variables[name] = _scalar(value)
            fixed.append(name)
    space = xr.Dataset(
        variables, coords={"wl": ("wl", np.asarray(wl, dtype=float))}
    )
    return space, fixed
