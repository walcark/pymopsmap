"""The .mopsmap accessor: what a result dataset can be turned into."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import xarray as xr

from pymopsmap.utils import get_logger

logger = get_logger(__name__)

# SMART-G reads its aerosol look-up tables under these dimension names.
SMARTG_DIMS = {"humidity": "hum", "wavelength": "wav", "mueller": "stk"}


@xr.register_dataset_accessor("mopsmap")
class MopsmapAccessor:
    """Conversions of a computed result into downstream formats."""

    def __init__(self, dataset: xr.Dataset) -> None:
        self._ds = dataset

    def to_smartg(
        self,
        path: str | Path,
        name: str,
        humidity_dim: str = "rh",
    ) -> xr.Dataset:
        """
        Write a look-up table readable by the AerOPAC class of SMART-G.

        Parameters
        ----------
        path : str or Path
            Destination NetCDF file.
        name : str
            Name stored in the table, which SMART-G uses to identify it.
        humidity_dim : str
            Which dimension of the result holds relative humidity. It is a
            parameter because the caller chooses that name when sweeping.

        Returns
        -------
        xr.Dataset
            The table that was written.
        """
        ds = self._ds
        if humidity_dim not in ds.dims:
            raise KeyError(
                f"No humidity dimension '{humidity_dim}' in the result; "
                f"it carries {sorted(map(str, ds.dims))}. Pass humidity_dim=."
            )
        if "phase" not in ds:
            raise KeyError(
                "No 'phase' in the result: request "
                "OutputType.PHASE_FUNCTION when computing it."
            )

        phase = ds["phase"]
        if "mueller_idx" not in phase.dims:
            phase = phase.expand_dims(mueller_idx=[0])

        order = (humidity_dim, "wl", "mueller_idx", "theta")
        table = xr.Dataset(
            {
                "phase": (
                    tuple(SMARTG_DIMS[k] for k in ("humidity", "wavelength"))
                    + (SMARTG_DIMS["mueller"], "theta"),
                    phase.transpose(*order).values,
                ),
                "ext": (
                    (SMARTG_DIMS["humidity"], SMARTG_DIMS["wavelength"]),
                    ds["kext"].transpose(humidity_dim, "wl").values,
                ),
                "ssa": (
                    (SMARTG_DIMS["humidity"], SMARTG_DIMS["wavelength"]),
                    ds["ssa"].transpose(humidity_dim, "wl").values,
                ),
            },
            coords={
                SMARTG_DIMS["humidity"]: ds[humidity_dim].values,
                SMARTG_DIMS["wavelength"]: ds["wl"].values,
                "theta": ds["theta"].values,
            },
            attrs={
                "name": name,
                "H_mix_min": 0,
                "H_mix_max": 99,
                "H_stra_min": 0,
                "H_stra_max": 0,
                "H_free_min": 0,
                "H_free_max": 0,
                "Z_mix": 2.0,
                "Z_free": 0.0,
                "Z_stra": 0.0,
                "date": datetime.today().strftime("%Y-%m-%d"),
                "source": "Created using MOPSMAP v1.0.",
            },
        )
        table.to_netcdf(path)
        logger.info("SMART-G look-up table written: %s", path)
        return table
