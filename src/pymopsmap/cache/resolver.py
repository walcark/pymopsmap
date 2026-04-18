"""Optical dataset resolver — maps refractive index grid to NC file paths."""

from __future__ import annotations

from itertools import product
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr

from pymopsmap.exceptions import IndexFileError

if TYPE_CHECKING:
    from pymopsmap.models import MicroParameters
    from pymopsmap.models.microparams import Shape


def _fmt_mreal(v: float) -> str:
    return f"{v:.4f}"


def _fmt_mimag(v: float) -> str:
    return f"{v:.6f}"


def _fmt_eps(v: float) -> str:
    s = f"{v:.3f}"
    if s[0] == " ":
        s = "0" + s[1:]
    return s


def _bracket(grid: np.ndarray, value: float) -> list[float]:
    """Return the 1 or 2 grid values that bracket `value`.

    When value lies exactly on a grid point, returns only that single value
    (no interpolation needed, and deduplication collapses to 1 file).
    """
    i = int(np.searchsorted(grid, value))
    # Exact match: no interpolation needed
    if i < len(grid) and np.isclose(grid[i], value, rtol=1e-9, atol=0):
        return [float(grid[i])]
    indices = set()
    if i > 0:
        indices.add(i - 1)
    if i < len(grid):
        indices.add(i)
    return sorted({grid[j] for j in indices})


class NCFileResolver:
    def __init__(self, index_path: Path):
        try:
            ds = xr.open_dataset(index_path)
        except Exception as exc:
            raise IndexFileError(
                f"Cannot open index.nc at {index_path}: {exc}"
            ) from exc
        self.avail_mreal: np.ndarray = np.sort(ds["mreal"].values.astype(float))
        self.avail_mimag: np.ndarray = np.sort(ds["mimag"].values.astype(float))
        if "eps" in ds:
            self.avail_eps: np.ndarray = np.sort(ds["eps"].values.astype(float))
        else:
            self.avail_eps = np.array([])
        ds.close()

    def resolve(self, mp: MicroParameters | list[MicroParameters]) -> list[str]:
        from pymopsmap.models import MicroParameters as MP

        modes = [mp] if isinstance(mp, MP) else mp

        files: set[str] = {"index.nc"}
        for mode in modes:
            for mr, mi in zip(mode.n_real, mode.n_imag):
                files.update(self._files_for_params(mode.shape, mr, mi))
        return sorted(files)

    def _files_for_params(self, shape: Shape, mreal: float, mimag: float) -> list[str]:
        from pymopsmap.models.microparams import (
            Irregular,
            IrregularDistrFile,
            IrregularOverlay,
            Sphere,
            Spheroid,
            SpheroidDistrFile,
            SpheroidLognormal,
        )

        mr_vals = _bracket(self.avail_mreal, mreal)
        mi_vals = _bracket(self.avail_mimag, mimag)

        if isinstance(shape, Sphere):
            return [
                f"spheres/sphere_{_fmt_mreal(mr)}_{_fmt_mimag(mi)}.nc"
                for mr, mi in product(mr_vals, mi_vals)
            ]

        if isinstance(shape, (Spheroid, SpheroidLognormal, SpheroidDistrFile)):
            eps_vals = self._eps_for_shape(shape)
            return [
                f"spheroids_merged/spheroid_merged_{_fmt_eps(eps)}_{_fmt_mreal(mr)}_{_fmt_mimag(mi)}.nc"
                for eps, mr, mi in product(eps_vals, mr_vals, mi_vals)
            ]

        if isinstance(shape, (Irregular, IrregularDistrFile, IrregularOverlay)):
            shape_id = self._irregular_id(shape)
            return [
                f"irregular/shape{shape_id}_{_fmt_mreal(mr)}_{_fmt_mimag(mi)}.nc"
                for mr, mi in product(mr_vals, mi_vals)
            ]

        return []

    def _eps_for_shape(self, shape: Shape) -> list[float]:
        from pymopsmap.models.microparams import Spheroid, SpheroidLognormal

        if len(self.avail_eps) == 0:
            return []

        if isinstance(shape, Spheroid):
            return _bracket(self.avail_eps, shape.aspect_ratio)

        if isinstance(shape, SpheroidLognormal):
            # Lognormal distribution covers a range; return all eps in range.
            ar = shape.aspect_ratio
            sigma = shape.sigma_ar
            lo = max(self.avail_eps[0], ar / (1 + 3 * sigma))
            hi = min(self.avail_eps[-1], ar * (1 + 3 * sigma))
            mask = (self.avail_eps >= lo) & (self.avail_eps <= hi)
            selected = self.avail_eps[mask].tolist()
            return selected if selected else _bracket(self.avail_eps, ar)

        # SpheroidDistrFile: can't know without reading file — use all eps
        return self.avail_eps.tolist()

    @staticmethod
    def _irregular_id(shape: Shape) -> str:
        from pymopsmap.models.microparams import (
            Irregular,
        )

        if isinstance(shape, Irregular):
            return shape.shape_id
        # For file-defined irregular, we can't know which shape_id, return a placeholder
        # that will cause a downstream error rather than silently miss a file.
        return "A"
