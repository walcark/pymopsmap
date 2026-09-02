"""MOPSMAP output parsers — stdout integrated block and ASCII output files."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import xarray as xr

from pymopsmap.models.output_request import OutputRequest, OutputType

if TYPE_CHECKING:
    pass


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------


def format_mopsmap_outputs(
    out_mopsmap: dict[str, Any],
    output_types: OutputRequest | None = None,
) -> xr.Dataset:
    """
    Build a result dataset from MOPSMAP output artefacts.

    out_mopsmap keys:
      - stdout      : captured stdout string
      - output_path : Path to netcdf output (always present)
      - ascii_base  : Path prefix for ASCII output files (present when any
                       non-integrated type was requested)
    """
    from pymopsmap.models.output_request import DEFAULT_OUTPUT

    if output_types is None:
        output_types = DEFAULT_OUTPUT

    datasets: list[xr.Dataset] = []

    ascii_base: Path | None = out_mopsmap.get("ascii_base")

    if OutputType.INTEGRATED in output_types:
        # When ascii_file is active, MOPSMAP writes integrated to
        # {base}.integrated (no prefix, just columns); stdout is empty.
        integrated_file = (
            Path(str(ascii_base) + ".integrated")
            if ascii_base is not None
            else None
        )
        if integrated_file is not None and integrated_file.exists():
            datasets.append(
                _parse_integrated_file(integrated_file.read_text())
            )
        else:
            datasets.append(format_stdout(out_mopsmap["stdout"]))

    wl = None
    if datasets:
        wl = datasets[0]["wl"].values

    for otype in output_types:
        if otype == OutputType.INTEGRATED:
            continue
        ext = otype.value  # matches the file extension written by MOPSMAP
        if ascii_base is not None:
            p = Path(str(ascii_base) + f".{ext}")
            if p.exists():
                text = p.read_text()
                ds = _parse_ascii(otype, text, wl)
                if ds is not None:
                    datasets.append(ds)

    if not datasets:
        raise ValueError(
            "No MOPSMAP output could be parsed for the requested output types."
        )

    # Deduplicate variables across datasets: integrated takes precedence.
    merged_vars: set[str] = set()
    clean: list[xr.Dataset] = []
    for ds in datasets:
        overlap = set(ds.data_vars) & merged_vars
        clean.append(ds.drop_vars(list(overlap)) if overlap else ds)
        merged_vars.update(ds.data_vars)  # type: ignore[arg-type]

    return xr.merge(clean, compat="no_conflicts", join="outer")


# ---------------------------------------------------------------------------
# Integrated (stdout)
# ---------------------------------------------------------------------------


def _integrated_dataset(arr: np.ndarray) -> xr.Dataset:
    """
    Build the integrated dataset from the twelve MOPSMAP columns.

    ksca is not one of them. It is derived here, once, rather than in every
    consumer that needs a scattering weight, and it is a quantity users want
    in its own right.
    """
    kext, ssa = arr[:, 1], arr[:, 2]
    return xr.Dataset(
        data_vars={
            "kext": (("wl",), kext),
            "ssa": (("wl",), ssa),
            "ksca": (("wl",), kext * ssa),
            "g": (("wl",), arr[:, 3]),
            "reff": (("wl",), arr[:, 4]),
            "n": (("wl",), arr[:, 5]),
            "cross_dens": (("wl",), arr[:, 6]),
            "vol_dens": (("wl",), arr[:, 7]),
            "mass_conc": (("wl",), arr[:, 8]),
            "angstrom_ext": (("wl",), arr[:, 9]),
            "angstrom_sca": (("wl",), arr[:, 10]),
            "angstrom_abs": (("wl",), arr[:, 11]),
        },
        coords={"wl": arr[:, 0].astype(np.float32)},
    )


def format_stdout(stdout: str) -> xr.Dataset:
    """
    Parse MOPSMAP stdout blocks split by "integrated".
    Columns (12): wl, kext, ssa, g, reff, n, cross_dens, vol_dens,
                  mass_conc, angstrom_ext, angstrom_sca, angstrom_abs
    """
    blocks = stdout.split("integrated")
    rows: list[np.ndarray] = []
    for blk in blocks[1:]:
        lines = [
            ln
            for ln in blk.splitlines()
            if ln.strip() and not ln.strip().startswith("Warning:")
        ]
        toks = " ".join(lines).split()
        if not toks:
            continue
        rows.append(np.asarray(toks, dtype=np.float32))

    if not rows:
        raise ValueError(
            "No numeric data found in stdout (split by 'integrated')."
        )

    return _integrated_dataset(np.stack(rows, axis=0))


def _parse_integrated_file(text: str) -> xr.Dataset:
    """
    Parse {ascii_base}.integrated written when ascii_file output is active.
    Same 12 columns as stdout but without the leading 'integrated' token.
    """
    return _integrated_dataset(_read_numeric(text))


# ---------------------------------------------------------------------------
# ASCII parsers
# ---------------------------------------------------------------------------


def _parse_ascii(
    otype: OutputType, text: str, wl: np.ndarray | None
) -> xr.Dataset | None:
    if otype == OutputType.PHASE_FUNCTION:
        return parse_phase_function(text)
    if otype == OutputType.SCATTERING_MATRIX:
        return parse_scattering_matrix(text)
    if otype == OutputType.VOLUME_SCATTERING_FUNCTION:
        return parse_volume_scattering_function(text)
    if otype == OutputType.LIDAR:
        return parse_lidar(text)
    if otype == OutputType.COEFF:
        return parse_coeff(text)
    return None


def _read_numeric(text: str) -> np.ndarray:
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("!"):
            continue
        rows.append([float(v) for v in line.split()])
    if not rows:
        raise ValueError("No numeric data in ASCII output.")
    return np.array(rows, dtype=np.float32)


def parse_phase_function(text: str) -> xr.Dataset:
    """
    Parse MOPSMAP phase_function ASCII output.
    Columns: wl  theta  a1
    Returns Dataset with phase(wl, theta).
    """
    arr = _read_numeric(text)
    wl_vals = np.unique(arr[:, 0])
    theta_vals = np.unique(arr[:, 1])

    phase = arr[:, 2].reshape(len(wl_vals), len(theta_vals))
    return xr.Dataset(
        data_vars={"phase": (("wl", "theta"), phase)},
        coords={"wl": wl_vals, "theta": theta_vals},
    )


def parse_scattering_matrix(text: str) -> xr.Dataset:
    """
    Parse MOPSMAP scattering_matrix ASCII output.
    Columns: wl  theta  a1  a2  a3  a4  b1  b2
    Returns Dataset with scattering_matrix(wl, theta, element).
    """
    arr = _read_numeric(text)
    wl_vals = np.unique(arr[:, 0])
    theta_vals = np.unique(arr[:, 1])
    elements = ["a1", "a2", "a3", "a4", "b1", "b2"]

    n_wl = len(wl_vals)
    n_th = len(theta_vals)
    mat = arr[:, 2:8].reshape(n_wl, n_th, 6)

    return xr.Dataset(
        data_vars={"scattering_matrix": (("wl", "theta", "element"), mat)},
        coords={"wl": wl_vals, "theta": theta_vals, "element": elements},
    )


def parse_volume_scattering_function(text: str) -> xr.Dataset:
    """
    Parse MOPSMAP volume_scattering_function ASCII output.
    Columns: wl  theta  a1_tilde
    Returns Dataset with vol_sca_func(wl, theta).
    """
    arr = _read_numeric(text)
    wl_vals = np.unique(arr[:, 0])
    theta_vals = np.unique(arr[:, 1])

    vsf = arr[:, 2].reshape(len(wl_vals), len(theta_vals))
    return xr.Dataset(
        data_vars={"vol_sca_func": (("wl", "theta"), vsf)},
        coords={"wl": wl_vals, "theta": theta_vals},
    )


def parse_lidar(text: str) -> xr.Dataset:
    """
    Parse MOPSMAP lidar ASCII output.
    Columns (9): wl  kext  backscatter  lidar_ratio  depol_ratio
                 angstrom_back  angstrom_ext  ext_to_mass  back_to_mass
    Returns Dataset with lidar variables indexed by wl.
    """
    arr = _read_numeric(text)
    wl = arr[:, 0]
    return xr.Dataset(
        data_vars={
            "backscatter": (("wl",), arr[:, 2]),
            "lidar_ratio": (("wl",), arr[:, 3]),
            "depol_ratio": (("wl",), arr[:, 4]),
            "angstrom_back": (("wl",), arr[:, 5]),
            "angstrom_ext": (("wl",), arr[:, 6]),
            "ext_to_mass": (("wl",), arr[:, 7]),
            "back_to_mass": (("wl",), arr[:, 8]),
        },
        coords={"wl": wl},
    )


def parse_coeff(text: str) -> xr.Dataset:
    """
    Parse MOPSMAP coeff ASCII output.
    Columns (8): wl  l  alpha1  alpha2  alpha3  alpha4  beta1  beta2
    Returns Dataset with coeff(wl, l, coeff_element).
    """
    arr = _read_numeric(text)
    wl_vals = np.unique(arr[:, 0])
    l_vals = np.unique(arr[:, 1]).astype(np.int32)
    coeff_names = ["alpha1", "alpha2", "alpha3", "alpha4", "beta1", "beta2"]

    n_wl = len(wl_vals)
    n_l = len(l_vals)
    data = arr[:, 2:8].reshape(n_wl, n_l, 6)

    return xr.Dataset(
        data_vars={"coeff": (("wl", "l", "coeff_element"), data)},
        coords={"wl": wl_vals, "l": l_vals, "coeff_element": coeff_names},
    )


# ---------------------------------------------------------------------------
# Legacy helper kept for compatibility (not used in new pipeline)
# ---------------------------------------------------------------------------


def format_netcdf_file(filename: Path, wl: np.ndarray) -> xr.Dataset:
    """Parse output.nc — kept for backward compatibility."""
    xrds = xr.open_dataset(filename)
    if "nreff" in xrds.dims:
        xrds = xrds.isel(nreff=0)

    if "phase" not in xrds:
        raise KeyError("Variable 'phase' not found in netCDF file.")

    theta = np.linspace(
        0.0, 180.0, len(xrds.nthetamax.data), dtype=np.float32
    )[::-1]
    wl = np.asarray(wl, dtype=np.float32)
    phase = np.asarray(xrds["phase"].data, dtype=np.float32)

    if phase.ndim != 3:
        raise ValueError(f"'phase' must be 3D, got shape {phase.shape}.")
    if phase.shape[0] == wl.size:
        pass
    elif phase.shape[1] == wl.size:
        phase = phase.transpose(1, 0, 2)
    else:
        raise ValueError(
            f"Cannot align phase with wl: "
            f"wl.size={wl.size}, phase.shape={phase.shape}."
        )

    mueller_idx = np.arange(phase.shape[1], dtype=np.int32)
    return xr.Dataset(
        data_vars={"phase": (("wl", "mueller_idx", "theta"), phase)},
        coords={"wl": wl, "mueller_idx": mueller_idx, "theta": theta},
    )
