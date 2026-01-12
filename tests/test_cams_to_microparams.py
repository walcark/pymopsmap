import xarray as xr
import netCDF4
import numpy as np

from pymopsmap.adapt.cams_to_microparams import (
    read_aerosol_modes_concentrations,
    CamsAerosol,
    CamsVersion,
    CAMS_DATA_PATH,
    _read_refractive_index,
)


def test_read_aerosol_mode_concentration():
    aer: CamsAerosol = CamsAerosol.AMMONIUM_CAMS

    concentration = read_aerosol_modes_concentrations(aer)
    assert isinstance(concentration, dict)
    assert concentration["fine"] >= 0
    assert concentration["coarse"] >= 0


def test_read_refractive_index():
    path = CAMS_DATA_PATH / "cams_aer_microphysical_parameters.nc"

    ds_test = xr.load_dataset(path)

    wl = np.linspace(0.550, 2.0, 10)

    for rh in [0.0, 50.0, 70.0, 95.0]:
        refr_indx = _read_refractive_index(ds_test, rh, wl)
        assert isinstance(refr_indx, dict)
        assert len(refr_indx["fine"]) == 2
        assert len(refr_indx["coarse"]) == 2
