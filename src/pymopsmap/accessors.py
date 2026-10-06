"""The .mopsmap accessor: what a result dataset can be turned into."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import numpy as np
import xarray as xr

from pymopsmap.utils import get_logger

logger = get_logger(__name__)

# SMART-G reads its aerosol look-up tables under these dimension names.
SMARTG_DIMS = {"humidity": "hum", "wavelength": "wav", "mueller": "stk"}

# SMART-G stores the six independent terms of the scattering matrix in the
# order F11, F21, F33, F34, F22, F44.  MOPSMAP emits them as a1, a2, a3, a4,
# b1, b2, which is the same matrix in Mishchenko's notation: a1 = F11,
# a2 = F22, a3 = F33, a4 = F44, b1 = F12 = F21 and b2 = -F34.  Hence the
# reordering, and the sign on the last term.
SMARTG_MUELLER: tuple[tuple[int, float], ...] = (
    (0, 1.0),  # F11 <- a1
    (4, 1.0),  # F21 <- b1
    (2, 1.0),  # F33 <- a3
    (5, -1.0),  # F34 <- -b2
    (1, 1.0),  # F22 <- a2
    (3, 1.0),  # F44 <- a4
)


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
        if "scattering_matrix" not in ds and "phase" not in ds:
            raise KeyError(
                "Neither 'scattering_matrix' nor 'phase' in the result: "
                "request OutputType.SCATTERING_MATRIX when computing it."
            )

        # SMART-G 2.0 requires strictly increasing scattering angles and
        # raises otherwise, where 1.x accepted either order.  MOPSMAP
        # emits them decreasing, so the table is sorted here rather than
        # left for every reader to discover.
        if ds["theta"].size > 1 and bool(ds["theta"][0] > ds["theta"][-1]):
            ds = ds.sortby("theta")

        phase = _smartg_phase(ds, humidity_dim)

        table = xr.Dataset(
            {
                "phase": (
                    tuple(SMARTG_DIMS[k] for k in ("humidity", "wavelength"))
                    + (SMARTG_DIMS["mueller"], "theta"),
                    phase,
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
                # MOPSMAP works in micrometres, SMART-G reads these tables
                # in nanometres.
                SMARTG_DIMS["wavelength"]: ds["wl"].values * 1e3,
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


def _smartg_phase(ds: xr.Dataset, humidity_dim: str) -> np.ndarray:
    """
    Build the ``(hum, wav, stk, theta)`` phase block SMART-G expects.

    A full scattering matrix is reordered into SMART-G's six terms. A
    scalar phase function fills F11 and F22 and leaves the four polarised
    terms at zero, which is what SMART-G itself does with a table that
    carries fewer than six: usable for a scalar run, wrong for a
    polarised one.

    Parameters
    ----------
    ds : xr.Dataset
        A computed result, holding ``scattering_matrix`` or ``phase``.
    humidity_dim : str
        Which dimension of the result holds relative humidity.

    Returns
    -------
    np.ndarray
        Shape ``(n_hum, n_wl, 6, n_theta)``.
    """
    if "scattering_matrix" in ds:
        matrix = ds["scattering_matrix"].transpose(
            humidity_dim, "wl", "element", "theta"
        )
        terms = [
            sign * matrix[:, :, i, :].values for (i, sign) in SMARTG_MUELLER
        ]
        return np.stack(terms, axis=2)

    logger.warning(
        "Only a scalar phase function was computed, so the four polarised "
        "terms of the table are zero. Request "
        "OutputType.SCATTERING_MATRIX for a polarised run."
    )
    f11 = ds["phase"].transpose(humidity_dim, "wl", "theta").values
    block = np.zeros((*f11.shape[:2], 6, f11.shape[2]), dtype=f11.dtype)
    block[:, :, 0, :] = f11
    block[:, :, 4, :] = f11
    return block
