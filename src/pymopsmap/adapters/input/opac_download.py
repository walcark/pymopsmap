"""Download OPAC microphysical parameters from the GEISHA database."""

import re
import subprocess

import xarray as xr
from tqdm import tqdm

from pymopsmap.utils import DATA_PATH

GEISHA_API_ROOT = "https://aeris-geisa.ipsl.fr/geisa_files/2011/Aerosols/OPAC/"

OPAC_SPECIES_SOLUBLE = {
    "inso": False,
    "soot": False,
    "waso": True,
    "ssam": True,
    "sscm": True,
    "suso": True,
    "minm": False,
    "miam": False,
    "micm": False,
    "mitr": False,
}

RHS = ("00", "50", "70", "80", "90", "95", "98", "99")


def parse_geisa_output(text: str) -> dict:
    """Parse a GEISA OPAC file, return PSD and refractive index table."""
    psd_keys = ("rmin", "rmax", "sigma", "rmod")
    psd_vals: list[float] = []
    wl: list[float] = []
    nr: list[float] = []
    ni: list[float] = []

    in_psd = False
    in_optics = False

    for line in text.splitlines():
        if "size distribution" in line:
            in_psd = True
            continue

        if in_psd and "optical parameters" in line:
            in_psd = False
            in_optics = True
            continue

        if in_psd and ("radius" in line or "sigma" in line or "Rmod" in line):
            m = re.findall(r"[-+]?\d*\.\d+E[+-]?\d+", line)
            if m:
                psd_vals.append(float(m[0]))

        if in_optics:
            nums = re.findall(r"[-+]?\d*\.\d+E[+-]?\d+", line)
            if len(nums) == 9:
                wl.append(float(nums[0]))
                nr.append(float(nums[7]))
                ni.append(float(nums[8]))

    return {
        **dict(zip(psd_keys, psd_vals[:4])),
        "wavelength": wl,
        "n_real": nr,
        "n_imag": ni,
    }


def fetch(specie: str) -> dict:
    """Download and parse a single OPAC GEISA species file."""
    url = f"{GEISHA_API_ROOT}{specie}"
    out = subprocess.run(
        ["curl", "-L", url],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    return parse_geisa_output(out)


def build_ds(specie: str, soluble: bool) -> xr.Dataset:
    """Build an xarray Dataset for one OPAC species over RH grid."""

    if soluble:
        datasets = []
        for rh in RHS:
            p = fetch(f"{specie}{rh}")
            ds_rh = xr.Dataset(
                data_vars={
                    "rmin": ("rh", [p["rmin"]]),
                    "rmax": ("rh", [p["rmax"]]),
                    "sigma": ("rh", [p["sigma"]]),
                    "rmod": ("rh", [p["rmod"]]),
                    "n_real": (["rh", "wavelength"], [p["n_real"]]),
                    "n_imag": (["rh", "wavelength"], [p["n_imag"]]),
                },
                coords={"rh": [int(rh)], "wavelength": p["wavelength"]},
            )
            datasets.append(ds_rh)
        ds = xr.concat(datasets, dim="rh")

    else:
        p = fetch(f"{specie}00")
        ds = xr.Dataset(
            data_vars={
                "rmin": p["rmin"],
                "rmax": p["rmax"],
                "sigma": p["sigma"],
                "rmod": p["rmod"],
                "n_real": ("wavelength", p["n_real"]),
                "n_imag": ("wavelength", p["n_imag"]),
            },
            coords={"wavelength": p["wavelength"]},
        ).expand_dims(rh=[int(rh) for rh in RHS])

    return ds


def download_opac_microparams():
    """
    Download and build OPAC microphysical parameter datasets.

    This function retrieves OPAC aerosol microphysical parameters from the
    GEISA database for all predefined aerosol species. Hygroscopic species
    are processed over a relative humidity grid (RHS), while non-hygroscopic
    species are replicated across all RH values.

    The resulting datasets contain the following variables:
        - rmin      : minimum particle radius [µm]
        - rmax      : maximum particle radius [µm]
        - sigma     : lognormal distribution width
        - rmod      : modal radius [µm]
        - n_real    : real refractive index (rh × wavelength)
        - n_imag    : imaginary refractive index (rh × wavelength),
                      negative convention

    Outputs are saved as NetCDF files, one per species.
    """
    length: int = len(OPAC_SPECIES_SOLUBLE)
    for specie, soluble in tqdm(
        OPAC_SPECIES_SOLUBLE.items(),
        total=length,
        desc="Download OPAC micro-parameters",
    ):
        ds = build_ds(specie, soluble)
        ds.to_netcdf(DATA_PATH / f"opac/{specie}.nc")


if __name__ == "__main__":
    download_opac_microparams()
