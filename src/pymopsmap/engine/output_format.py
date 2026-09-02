"""MOPSMAP output parsers — stdout integrated block and ASCII output files."""

from __future__ import annotations

import warnings
from collections.abc import Callable, Iterable
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


# ---------------------------------------------------------------------------
# Combination of several species into one external mixture
# ---------------------------------------------------------------------------
#
# External mixing is additive on the coefficients, so a mixture is the sum of
# its species rather than a separate MOPSMAP run. What is additive, and what
# has to be rebuilt from the additive quantities, is settled by reading how
# MOPSMAP forms each output. Line references are to bin/mopsmap/src.
#
# The rules live here, beside the parsers that produce the variables, so a new
# output cannot be added without one: combine() raises on any variable it does
# not know.

CombineRule = Callable[[str, list[xr.Dataset]], xr.DataArray]


def _total(items: Iterable[xr.DataArray]) -> xr.DataArray:
    """Sum data arrays without starting from a scalar zero."""
    total, *rest = items
    for item in rest:
        total = total + item
    return total


def _additive(name: str, parts: list[xr.Dataset]) -> xr.DataArray:
    """Densities per unit volume of air; they simply add."""
    return _total(part[name] for part in parts)


def _weighted_mean(weight: str) -> CombineRule:
    """
    Intensive quantities, averaged over what they are intensive against.

    An asymmetry parameter is the mean of cos(theta) over scattered photons,
    so it is weighted by scattering. A ratio of two additive quantities is the
    ratio of their sums, which is the same thing: for X = num/den, the sum of
    the numerators over the sum of the denominators is the den-weighted mean
    of X. That covers lidar_ratio (kext/backscatter,
    write_output_ascii.f90:305), ext_to_mass (mass/cext, line 309) and
    back_to_mass (backscatter/mass, line 310) with no constant to carry.
    """

    def rule(name: str, parts: list[xr.Dataset]) -> xr.DataArray:
        return _total(part[name] * part[weight] for part in parts) / _total(
            part[weight] for part in parts
        )

    return rule


def _single_scattering_albedo(
    name: str, parts: list[xr.Dataset]
) -> xr.DataArray:
    """A ratio of sums, never the mean of the individual albedos."""
    return _additive("ksca", parts) / _additive("kext", parts)


def _angstrom(source: Callable[[xr.Dataset], xr.DataArray]) -> CombineRule:
    """
    Rebuild a spectral log-derivative from the combined spectrum.

    The exponent of a sum is not the sum of the exponents. MOPSMAP computes it
    per interval and leaves the first wavelength undefined
    (write_output_ascii.f90:183 writes sqrt(-1), that is NaN, then line 188
    uses the previous wavelength).
    """

    def rule(name: str, parts: list[xr.Dataset]) -> xr.DataArray:
        total = _total(source(part) for part in parts)
        wl = parts[0]["wl"]
        values = np.asarray(total)
        exponent = np.full(values.shape, np.nan, dtype=float)
        # A non-absorbing mixture has kext == ksca, so the absorption
        # exponent is a 0/0: NaN is the answer, the warning is not.
        with np.errstate(divide="ignore", invalid="ignore"):
            exponent[1:] = -np.log(values[1:] / values[:-1]) / np.log(
                wl.values[1:] / wl.values[:-1]
            )
        return xr.DataArray(exponent, coords={"wl": wl}, dims=("wl",))

    return rule


def _depolarisation(name: str, parts: list[xr.Dataset]) -> xr.DataArray:
    """
    Split the backscatter into its two components, add, and re-form the ratio.

    MOPSMAP reports the total backscatter (write_output_ascii.f90:304) and
    delta = (1 - a2/a1)/(1 + a2/a1) at 180 degrees (line 306). Both a1 and a2
    are scattering-weighted sums over the modes, so the parallel and
    perpendicular components are additive even though their ratio is not.
    """
    parallel = _total(
        part["backscatter"] / (1.0 + part[name]) for part in parts
    )
    perpendicular = _total(
        part["backscatter"] * part[name] / (1.0 + part[name]) for part in parts
    )
    return perpendicular / parallel


def _effective_radius(name: str, parts: list[xr.Dataset]) -> xr.DataArray:
    """
    Rebuild the effective radius from the summed moments, for spheres only.

    MOPSMAP defines reff as r3_n_sum/r2_n_sum (mix_contributions.f90:64), two
    raw moments it does not output. What it does output carries the geometric
    constants and, for non-spherical shapes, different powers of the aspect
    ratio: cs_sum gets pi times r2_n_sum, or rat**-2 or rat**-6 times that,
    and vol_sum gets 4/3 pi times r3_n_sum, or rat**3 or rat**-6 times it
    (add_contribution.f90:269-289). The moments are therefore recoverable only
    when every mode is a sphere, where the ratio reduces to 0.75 vol/cs.
    """
    shapes = {
        shape for part in parts for shape in part.attrs.get("shape_types", [])
    }
    volume = _additive("vol_dens", parts)
    if shapes != {"sphere"}:
        warnings.warn(
            "pymopsmap: 'reff' is only defined for spherical modes when "
            f"combining species; got shapes {sorted(shapes)}. The mixed "
            "effective radius is NaN.",
            UserWarning,
            stacklevel=3,
        )
        return xr.full_like(volume, np.nan)
    return 0.75 * volume / _additive("cross_dens", parts)


COMBINE: dict[str, CombineRule] = {
    # Densities per unit volume of air
    "kext": _additive,
    "ksca": _additive,
    "backscatter": _additive,
    "n": _additive,
    "mass_conc": _additive,
    "vol_dens": _additive,
    "cross_dens": _additive,
    # Scattering-weighted
    "g": _weighted_mean("ksca"),
    "phase": _weighted_mean("ksca"),
    "scattering_matrix": _weighted_mean("ksca"),
    "vol_sca_func": _weighted_mean("ksca"),
    "coeff": _weighted_mean("ksca"),
    # Ratios of additive quantities
    "ssa": _single_scattering_albedo,
    "lidar_ratio": _weighted_mean("backscatter"),
    "ext_to_mass": _weighted_mean("kext"),
    "back_to_mass": _weighted_mean("mass_conc"),
    "depol_ratio": _depolarisation,
    # Rebuilt from the additive quantities
    "reff": _effective_radius,
    "angstrom_ext": _angstrom(lambda ds: ds["kext"]),
    "angstrom_sca": _angstrom(lambda ds: ds["ksca"]),
    "angstrom_abs": _angstrom(lambda ds: ds["kext"] - ds["ksca"]),
    "angstrom_back": _angstrom(lambda ds: ds["backscatter"]),
}


def combine(parts: list[xr.Dataset], weights: list[float]) -> xr.Dataset:
    """
    Combine per-species results into one external mixture.

    Parameters
    ----------
    parts : list of xr.Dataset
        One result per species, all over the same wavelength grid.
    weights : list of float
        A scale factor per species. Optical properties are linear in the
        number concentration, so scaling here is exact and needs no rerun.

    Returns
    -------
    xr.Dataset
        The mixture, carrying the same variables as its parts.

    Raises
    ------
    KeyError
        If a variable has no combination rule.
    """
    if len(parts) != len(weights):
        raise ValueError(
            f"parts and weights must have the same length, got "
            f"{len(parts)} and {len(weights)}."
        )

    scaled = [
        part.assign(
            {
                str(name): part[name] * weight
                for name in part.data_vars
                if COMBINE.get(str(name)) is _additive
            }
        )
        for part, weight in zip(parts, weights)
    ]

    names = sorted(str(name) for name in scaled[0].data_vars)
    unknown = [name for name in names if name not in COMBINE]
    if unknown:
        raise KeyError(
            f"No combination rule for {unknown}. Declare one in COMBINE, "
            "beside the parser that produces the variable."
        )

    out = xr.Dataset(
        {name: COMBINE[name](name, scaled) for name in names},
        coords=scaled[0].coords,
    )
    out.attrs["shape_types"] = sorted(
        {s for part in parts for s in part.attrs.get("shape_types", [])}
    )
    return out
