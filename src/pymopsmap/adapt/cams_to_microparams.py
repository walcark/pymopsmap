"""
cams_to_microparams.py

Author  : Kévin Walcarius
Date    : 2025-01-08
Version : 1.0
License : MIT
Summary : Module to transform input CAMS formatted microparameters
          to an indexed set of MicroParams usable by Mopsmap.

Description:

The input CAMS formatted microparameters is a convention used in
the MAJA software to store the CAMS microphysical parameters.

Maja assumes (as for now) that each CAMS aerosol is a combination
of a fine and coarse mode (e.g two size distributions). Sometimes,
a CAMS aerosol only has one mode, and thus fine is equal to coarse.
Each mode is defined by a Log-Normal size distribution, whose
parameters are specified in the micrphysical parameters file.

The microphysical parameters are stored in the folder:

    pymopsmap/data/cams

The file organisation of the folder is the following:
    - a file stores the aerosols mode concentrations for each mode
    - a file stores the microphysical parameters for each mode
"""

import numpy as np
import json
import xarray as xr
from enum import Enum
from numpy.typing import ArrayLike
from pymopsmap.classes.microparams import LognormalPSD, MicroParameters, Sphere
from pymopsmap.utils import ROOT_PATH, PosFloat64List, SortedPosFloat64List

CAMS_DATA_PATH = ROOT_PATH / "data/cams"


# --------------------------------------------------------------------------
# Enumeration for the CAMS arguments
# --------------------------------------------------------------------------
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


# --------------------------------------------------------------------------
# Read microphysical parameters
# --------------------------------------------------------------------------
"""
The aerosol microphysical parameters must be a netcdf file. The netcdf
file store a xarray-like object with the given structure:

...
"""


def read_aerosol_microphysical_parameters(
    aerosol: CamsAerosol,
    version: CamsVersion,
    wl_microns: SortedPosFloat64List,
    rh: ArrayLike,
) -> tuple[list[dict[str, float]], list[list[MicroParameters]]]:
    """
    Extracts the microphysical parameters of a given CAMS aerosol. The
    parameters are the following:
        - wavelength
        - n_real (fine and coarse): the real part of the refr index
        - n_imag (fine and coarse): the imaginary part of the refr index
        - rm (fine and coarse): the modal radius of the Log-Normal size
          distribution
        - sigma (fine and coarse): the standard deviation of the log-normal
          size distribution
    """
    path = CAMS_DATA_PATH / "cams_aer_microphysical_parameters.nc"

    try:
        xrds = xr.load_dataset(path)
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Unable to open file: {path}") from e

    conc = read_aerosol_modes_concentrations(aerosol)

    aer, ver = aerosol.value, version.value

    rh = [rh] if isinstance(rh, float) else rh
    vars_li = [{"rh": relhum} for relhum in rh]

    mps = []
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
        mps.append(modes)

    return vars_li, mps


def _read_granulometry(
    ds: xr.Dataset, rh: float
) -> dict[str, tuple[float, float]]:
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
    Read the refractive index from the microphysical parameters dataset.
    The relative humidity `rh` is a float because the aerosol have a
    size distribution varying with humidity, and thus we need to perform
    a unique Mopsmap run for each humidity level.
    """
    assert isinstance(rh, float)

    wl_nm = np.array(wl_microns) * 1e3
    values = ds.interp(relative_humidity=rh, wavelength=wl_nm)

    nr_fine = list(values["mr_f"].data)
    nr_coarse = list(values["mr_c"].data)
    ni_fine = list(-values["mi_f"].data)
    ni_coarse = list(-values["mi_c"].data)

    return {"fine": (nr_fine, ni_fine), "coarse": (nr_coarse, ni_coarse)}


# --------------------------------------------------------------------------
# Read aerosol mode concentrations
# --------------------------------------------------------------------------
"""
The aerosol mode concentrations must be a JSON file with a dictionnary
of aerosol -> [fine, coarse] concentration. For instance:

{
    "continen": [1.0, 0.0],
    "black_carbon": [1.0, 0.0],
    "sulphate": [1.0, 0.0],
    "organic_matter": [1.0, 0.0],
    "nitrate": [1.0, 0.0],
    "ammonium": [1.0, 0.0],
    "secondary_organic": [1.0, 0.0],
    "sea_salt": [70.0, 3.0],
    "dust": [391.0, 8.39]
}
"""


def read_aerosol_modes_concentrations(
    aerosol: CamsAerosol,
) -> dict[str, float]:
    """
    Returns a list [CVfine, CVcoarse] where CVfine is the number of fine
    particles per m³ and CVcoarse the number of coarse particle per m³,
    for the CAMS aerosol `aerosol`.
    """
    data_path = CAMS_DATA_PATH / "cams_aer_modes_concentrations.json"
    try:
        with open(data_path, "r") as f:
            data = json.load(f)
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Unable to open file: {data_path}") from e

    aer = aerosol.value
    if aer not in data:
        raise KeyError(f"Unable to find hydrophilicity for specie: {aer}")
    conc = data[aer]
    return {"fine": conc[0], "coarse": conc[1]}


if __name__ == "__main__":
    from pymopsmap.mopsmap import compute_optical_properties
    from pymopsmap.classes import extend_optiprops

    index, mps = read_aerosol_microphysical_parameters(
        aerosol=CamsAerosol.CONTINENTAL,
        version=CamsVersion.V49_R1,
        wl_microns=np.linspace(0.490, 1.5, 100),
        rh=[0.0, 50.0, 90.0],
    )

    ops = [compute_optical_properties(mp=mp) for mp in mps]

    op_tot = extend_optiprops(index, ops)

    print(op_tot)
