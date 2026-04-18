"""CAMS aerosol adapter — interpolates microphysical parameters from CAMS tables."""

import json
from enum import Enum

import numpy as np
import xarray as xr
from numpy.typing import ArrayLike

from pymopsmap.classes import (
    LognormalPSD,
    MicroParameters,
    MicroParametersDispatch,
    Sphere,
)
from pymopsmap.utils import DATA_PATH, SortedPosFloat64List


class CamsAerosol(str, Enum):
    CONTINENTAL = "continen"
    SULPHATE_CAMS = "sulphate"
    SEA_SALT_CAMS = "sea_salt"
    DUST_CAMS = "dust"
    BLACK_CARBON_CAMS = "black_carbon"
    NITRATE_CAMS = "nitrate"
    AMMONIUM_CAMS = "ammonium"
    ORGANIC_MATTER_CAMS = "organic_matter"
    SECONDARY_ORGANIC = "secondary_organic"


class CamsVersion(str, Enum):
    V47_R1 = "47r1"
    V48_R1 = "48r1"
    V49_R1 = "49r1"


def read_aerosol_microphysical_parameters(
    aerosol: CamsAerosol,
    version: CamsVersion,
    wl_microns: SortedPosFloat64List,
    rh: ArrayLike,
) -> MicroParametersDispatch:
    """
    Extracts the microphysical parameters of a given CAMS aerosol, and
    then dispatch each instance
    """
    path = DATA_PATH / "cams/cams_aer_microphysical_parameters.nc"

    try:
        xrds = xr.load_dataset(path)
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Unable to open file: {path}") from e

    conc = read_aerosol_modes_concentrations(aerosol)

    aer, ver = aerosol.value, version.value

    rh = [rh] if isinstance(rh, float) else rh

    dispatch: MicroParametersDispatch = MicroParametersDispatch()
    for relhum in rh:
        xrds_sel = xrds.sel(aerosols_species=aer, cams_versions=ver)
        granulo = _read_granulometry(xrds_sel, relhum)
        refr_index = _read_refractive_index(xrds_sel, relhum, wl_microns)
        modes = []
        for mode in ["fine", "coarse"]:
            mp: MicroParameters = MicroParameters(
                wavelength=wl_microns,
                n_real=refr_index[mode][0],
                n_imag=refr_index[mode][1],
                shape=Sphere(),
                psd=LognormalPSD(
                    rm=np.round(granulo[mode][0], 6),
                    sigma=np.round(np.exp(granulo[mode][1]), 6),
                    n=conc[mode],
                    rmin=0.001,
                    rmax=40.0,
                ),
            )
            modes.append(mp)
        dispatch.append(modes=modes, params={"rh": relhum})

    return dispatch


def read_aerosol_modes_concentrations(
    aerosol: CamsAerosol,
) -> dict[str, float]:
    """
    Parse the JSON file containing modes concentration for each specie.
    Modes are organized as follow :

    {
        "specie1": [CVfine, CVcoarse],
        "specie2": ...
    }

    where (CVfine, CVcoarse) are the number of particles per m³ of the
    fine and coarse modes.

    Parameters
    ----------
    aerosol : CamsAerosol
        The CAMS aerosol of interest.

    Returns
    -------
    dict[str, float]
        A dictionnary {"fine": CVfine, "coarse": CVcoarse}.
    """
    data_path = DATA_PATH / "cams/cams_aer_modes_concentrations.json"
    try:
        with open(data_path) as f:
            data = json.load(f)
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Unable to open file: {data_path}") from e

    aer = aerosol.value
    if aer not in data:
        raise KeyError(f"Unable to find hydrophilicity for specie: {aer}")
    conc = data[aer]
    return {"fine": conc[0], "coarse": conc[1]}


def _read_granulometry(
    ds: xr.Dataset, rh: float
) -> dict[str, tuple[float, float]]:
    """
    Read the log-normal distribution arguments (rm, sigma) for the fine
    and coarse modes in an xarray dataset, for a given relative humidity.

    Parameters
    ----------
    ds : xr.Dataset
        The xarray dataset of the CAMS specie of interest.
    rh : float
        The relative humidity [%].

    Returns
    -------
    dict[str, tuple[float, float]]
        The dictionary of log-normal distribution arguments:
        {"fine": (rm_f, sigma_f), "coarse": (rm_c, sigma_c)}.
    """
    assert isinstance(rh, float)

    values = ds.interp(relative_humidity=rh)

    rm_fine = values["rmodal_f"].data
    rm_coarse = values["rmodal_c"].data
    sigma_fine = values["lnvar_f"].data
    sigma_coarse = values["lnvar_c"].data

    return {"fine": (rm_fine, sigma_fine), "coarse": (rm_coarse, sigma_coarse)}


def _read_refractive_index(
    ds: xr.Dataset, rh: float, wl_microns: ArrayLike
) -> dict[str, tuple[list[float], list[float]]]:
    """
    Read the refractive index (nr, ni) for the fine and coarse modes in
    an xarray dataset, for a given relative humidity and wavelength.
    from the microphysical parameters dataset.

    Parameters
    ----------
    ds : xr.Dataset
        The xarray dataset of the CAMS specie of interest.
    rh : float
        The relative humidity [%].
    wl_microns : ArrayLike
        The wavelengths for which to get the refractive indexes.

    Returns
    -------
    dict[str, tuple[list[float], list[float]]]
        The dictionary of refractive index real and imaginary parts.
        {"fine": (nr_arr_f, ni_arr_f), "coarse": (nr_arr_c, ni_arr_c)}.
    """
    assert isinstance(rh, float)

    wl_nm = np.array(wl_microns) * 1e3
    values = ds.interp(relative_humidity=rh, wavelength=wl_nm)

    nr_fine = list(values["mr_f"].data)
    nr_coarse = list(values["mr_c"].data)
    ni_fine = list(-values["mi_f"].data)
    ni_coarse = list(-values["mi_c"].data)

    return {"fine": (nr_fine, ni_fine), "coarse": (nr_coarse, ni_coarse)}
