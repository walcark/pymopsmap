"""ResultCache: read/write OptiProps to a blake2b-keyed NetCDF store."""

from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING

from pymopsmap.utils import RESULT_CACHE_DIR
from pymopsmap.utils.caching import cache_key

if TYPE_CHECKING:
    from pymopsmap.models import (
        MicroParameters,
        OptiProps,
    )
    from pymopsmap.models.output_request import OutputRequest


class ResultCache:
    def __init__(self, root_dir: Path | None = None):
        self.root_dir = root_dir or RESULT_CACHE_DIR
        self.root_dir.mkdir(parents=True, exist_ok=True)

    def key(
        self,
        mp: MicroParameters | list[MicroParameters],
        outputs: OutputRequest,
    ) -> str:
        return cache_key(mp, extra=sorted(outputs))

    def _path(self, key: str) -> Path:
        return self.root_dir / f"{key}.nc"

    def get(self, key: str) -> OptiProps | None:
        import xarray as xr

        from pymopsmap.models import OptiProps

        p = self._path(key)
        if not p.exists() or p.stat().st_size == 0:
            return None
        try:
            ds = xr.open_dataset(p)
            return OptiProps(ds=ds)
        except Exception:
            return None

    def put(self, key: str, result: OptiProps) -> None:
        target = self._path(key)
        tmp = target.with_suffix(".tmp")
        result.ds.to_netcdf(tmp)
        os.replace(tmp, target)
