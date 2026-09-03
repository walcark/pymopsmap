"""Size-parameter coverage checks (Gasteiger & Wiegner 2018 limits)."""

from __future__ import annotations

import math
import warnings
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import xarray as xr

    from pymopsmap.microparams import MicroParameters

# Hard limits from Gasteiger & Wiegner (2018, GMD 11:2739–2762), Tables 1 & 2.
# x = 2π r / λ  (volume-equivalent size parameter for irregular shapes).
# Keys match Shape.type discriminator values.
_X_LIMITS: dict[str, tuple[float, float]] = {
    "sphere": (1e-6, 1005.0),
    "spheroid": (1e-6, 1005.0),
    "spheroid-lognormal": (1e-6, 1005.0),
    "spheroid-distr-file": (1e-6, 1005.0),
    "irregular": (1e-3, 30.2),
    "irregular-distr-file": (1e-3, 30.2),
    "irregular-overlay": (1e-3, 30.2),
}


def _max_radius(psd) -> float | None:
    """Return the maximum relevant radius (µm) for a coverage check."""
    from pymopsmap.psd import (
        DistrListPSD,
        FixedPSD,
        LognormalPSD,
        ModifiedGammaPSD,
    )

    if isinstance(psd, FixedPSD):
        return psd.radius
    if isinstance(psd, (LognormalPSD, ModifiedGammaPSD)):
        return psd.rmax
    if isinstance(psd, DistrListPSD):
        return float(max(psd.radii))
    return None  # FileDefinedPSD: cannot determine without reading the file


def _valid_mask(mp: MicroParameters, rh: float | None) -> np.ndarray:
    """
    Boolean mask: True where the size parameter is within dataset limits.

    The radius checked is the one MOPSMAP will work on. When growth is
    delegated to the engine it scales every radius by the growth factor
    (calc_hygroscopic_growth.f90), so a mode that fits dry can leave the
    coverage once it takes up water.
    """
    from pymopsmap.scatlib.growth import growth_factor_cubed

    x_min, x_max = _X_LIMITS.get(mp.shape.type, (0.0, float("inf")))
    r_max = _max_radius(mp.psd)

    if r_max is None:
        return np.ones(len(mp.wavelength), dtype=bool)

    if mp.kappa and rh:
        r_max *= growth_factor_cubed(mp.kappa, rh) ** (1.0 / 3.0)

    wl = np.asarray(mp.wavelength, dtype=float)
    x = 2.0 * math.pi * r_max / wl
    return (x >= x_min) & (x <= x_max)


def _clip_mp(mp: MicroParameters, mask: np.ndarray) -> MicroParameters:
    from pymopsmap.microparams import MicroParameters

    return MicroParameters(
        wavelength=[w for w, v in zip(mp.wavelength, mask) if v],
        n_real=[r for r, v in zip(mp.n_real, mask) if v],  # type: ignore[arg-type]
        n_imag=[i for i, v in zip(mp.n_imag, mask) if v],  # type: ignore[arg-type]
        shape=mp.shape,
        psd=mp.psd,
        kappa=mp.kappa,
        density=mp.density,
    )


def clip_modes_to_coverage(
    modes: MicroParameters | list[MicroParameters],
    rh: float | None = None,
) -> tuple[MicroParameters | list[MicroParameters], np.ndarray]:
    """
    Clip wavelengths outside the Gasteiger & Wiegner (2018) coverage.

    Wavelengths outside the valid size-parameter range emit a UserWarning
    and are dropped before the MOPSMAP run. The caller must reindex the
    result back to the original grid (clipped positions become NaN).

    Args:
        modes: the modes of the run.
        rh: the humidity handed to MOPSMAP, when it applies the growth itself.

    Returns:
        clipped_modes — same type as input, out-of-range wavelengths dropped.
        valid_mask    — boolean over the *original* wl axis (True = kept).
    """
    from pymopsmap.microparams import MicroParameters as MP

    single = isinstance(modes, MP)
    if single:
        assert isinstance(modes, MP)
        mp_list: list[MP] = [modes]
    else:
        assert isinstance(modes, list)
        mp_list = list(modes)

    per_mask = [_valid_mask(mp, rh) for mp in mp_list]
    combined: np.ndarray = np.ones(len(mp_list[0].wavelength), dtype=bool)
    for m in per_mask:
        combined &= m

    if combined.all():
        return modes, combined

    # Build warning message
    n_clipped = int((~combined).sum())
    wl = np.asarray(mp_list[0].wavelength, dtype=float)
    clipped_wl = wl[~combined]

    parts: list[str] = []
    for mp, mask in zip(mp_list, per_mask):
        if mask.all():
            continue
        r_max = _max_radius(mp.psd)
        x_min, x_max = _X_LIMITS.get(mp.shape.type, (0.0, float("inf")))
        if r_max is not None:
            x_clip = 2.0 * math.pi * r_max / clipped_wl
            parts.append(
                f"shape '{mp.shape.type}':"
                f" x ∈ [{x_clip.min():.3g}, {x_clip.max():.3g}]"
                f" exceeds limit [{x_min:.3g}, {x_max:.3g}]"
            )
        else:
            parts.append(
                f"shape '{mp.shape.type}': limit [{x_min:.3g}, {x_max:.3g}]"
            )

    warnings.warn(
        f"pymopsmap: {n_clipped}/{len(combined)} wavelength(s) clipped. "
        f"size-parameter coverage (Gasteiger & Wiegner 2018 Tables 1–2). "
        + "; ".join(parts)
        + ". Out-of-range positions will be NaN in the result.",
        UserWarning,
        stacklevel=4,
    )

    clipped = [_clip_mp(mp, combined) for mp in mp_list]
    result = clipped[0] if single else clipped
    return result, combined


def reindex_to_full_grid(
    ds: xr.Dataset,
    original_wl: list[float],
    valid_mask: np.ndarray,
) -> xr.Dataset:
    """
    Restore the clipped wavelengths, filling the dropped ones with NaN.

    The kept wavelengths are reassigned from the original grid before
    reindexing: MOPSMAP echoes wavelengths as float32 strings that do not
    always round-trip, so the result cannot be matched on its own values.
    Reindexing on the coordinate, rather than assigning by position, keeps
    variables that carry extra dimensions such as theta, l or element.

    Parameters
    ----------
    ds : xr.Dataset
        Result over the kept wavelengths only.
    original_wl : list[float]
        The full wavelength grid the caller requested.
    valid_mask : np.ndarray
        Boolean over ``original_wl``, True where the wavelength was kept.

    Returns
    -------
    xr.Dataset
        Result over the full grid, NaN at the clipped positions.
    """
    wl_full = np.asarray(original_wl, dtype=np.float32)
    if valid_mask.all():
        return ds
    return ds.assign_coords(wl=wl_full[valid_mask]).reindex(wl=wl_full)
