"""
The refractive index MOPSMAP will use, once it has grown the particles.

When growth is delegated to the engine, the launch file carries the dry
refractive index and MOPSMAP mixes it with water itself
(calc_hygroscopic_growth.f90). The dataset files it opens are therefore not
the ones the dry index points at, and the resolver has to predict them.

The formula is applied here for resolution only: the optical properties still
come from the Fortran, so this is not a second implementation of the physics.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from pymopsmap.utils import MOPSMAP_PATH

WATER_REFRACTIVE_INDEX_FILE = (
    Path(MOPSMAP_PATH).parent / "data" / "refr_water_segelstein"
)


def growth_factor_cubed(kappa: float, rh: float) -> float:
    """Volume growth factor of the kappa parameterisation."""
    return 1.0 + kappa * rh / (100.0 - rh)


def grown_refractive_index(
    wavelength: list[float],
    n_real: list[float],
    n_imag: list[float],
    kappa: float,
    rh: float,
) -> tuple[list[float], list[float]]:
    """
    Mix a dry refractive index with water, as MOPSMAP does.

    Parameters
    ----------
    wavelength : list of float
        Wavelengths in micrometres.
    n_real, n_imag : list of float
        The dry refractive index.
    kappa : float
        Hygroscopicity parameter.
    rh : float
        Relative humidity in percent.

    Returns
    -------
    tuple
        The real and imaginary parts after growth.
    """
    if kappa <= 0.0 or rh <= 0.0:
        return list(n_real), list(n_imag)

    gf3 = growth_factor_cubed(kappa, rh)
    water_real, water_imag = _water_refractive_index(wavelength)
    grown = [
        ((np.asarray(part) + water * (gf3 - 1.0)) / gf3).tolist()
        for part, water in ((n_real, water_real), (n_imag, water_imag))
    ]
    return grown[0], grown[1]


def _water_refractive_index(
    wavelength: list[float],
) -> tuple[np.ndarray, np.ndarray]:
    """Segelstein water, on the requested wavelengths, as MOPSMAP reads it."""
    table = np.loadtxt(WATER_REFRACTIVE_INDEX_FILE, comments="#")
    return (
        np.interp(wavelength, table[:, 0], table[:, 1]),
        np.interp(wavelength, table[:, 0], table[:, 2]),
    )
