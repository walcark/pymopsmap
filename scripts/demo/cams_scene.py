"""
Optical properties of a real CAMS scene, over three successive dates.

    pixi run -e dev python -m scripts.demo.cams_scene

Downloads a window of the CAMS global reanalysis (EAC4) from the Atmosphere
Data Store, then computes the optical properties of every pixel. Two fields
vary across the scene and neither costs an extra MOPSMAP run:

  * relative humidity, which sets the state of each species,
  * the aerosol composition, whose mass concentrations weight the mixture.

EAC4 publishes mass mixing ratios, so the only conversion on the way in is by
the air density, computed from the temperature the same request brings back.

Requires ~/.cdsapirc and the dataset licence accepted once on the ADS.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import xarray as xr

import pymopsmap as pm

DATASET = "cams-global-reanalysis-eac4"
CACHE = Path("cams_scene.nc")
OUTPUT = Path("cams_scene.png")

# The North Atlantic and western Europe during storm Nelson, which drove a
# strong marine flux onto the continent at the end of March 2024. The domain
# is deliberately large: the sweep costs one run per distinct humidity, not
# one per pixel, so widening the scene is nearly free.
AREA = [62.0, -35.0, 30.0, 25.0]  # north, west, south, east
DATES = ["2024-03-27", "2024-03-28", "2024-03-29"]
TIME = "12:00"
PRESSURE_LEVEL = "1000"

# What to ask the store for, and the GRIB short name it answers with. EAC4
# carries no nitrate or ammonium: those species entered the IFS after the
# reanalysis was produced. Dust, black carbon and organic matter are left out
# for another reason: their refractive index falls outside the main archive of
# the optical dataset.
REQUEST = {
    "sulphate_aerosol_mixing_ratio": ["aermr11"],
    "sea_salt_aerosol_0.03-0.5um_mixing_ratio": ["aermr01"],
    "sea_salt_aerosol_0.5-5um_mixing_ratio": ["aermr02"],
    "sea_salt_aerosol_5-20um_mixing_ratio": ["aermr03"],
}

# The catalogue species each set of bins feeds. Sea salt is published in three
# size bins; the catalogue species already carries its own fine and coarse
# modes, so the bins are summed into one mass.
SPECIES = {
    pm.CAMS.SULPHATE: ["aermr11"],
    pm.CAMS.SEA_SALT: ["aermr01", "aermr02", "aermr03"],
}

WAVELENGTHS = [0.44, 0.55, 0.67, 0.87]
TARGET_UM = 0.55

DRY_AIR_GAS_CONSTANT = 287.058  # J kg-1 K-1
PASCALS = float(PRESSURE_LEVEL) * 100.0

# The CAMS tables stop at 95 percent, and a reanalysis reports more than that.
MAX_TABULATED_RH = 95.0


def download() -> xr.Dataset:
    """Fetch the scene, or reuse what is already on disk."""
    if CACHE.exists():
        print(f"reusing    : {CACHE}")
        return xr.open_dataset(CACHE)

    import cdsapi

    print(f"requesting : {DATASET} over {AREA} for {DATES}")
    cdsapi.Client().retrieve(
        DATASET,
        {
            "variable": [*REQUEST, "relative_humidity", "temperature"],
            "pressure_level": [PRESSURE_LEVEL],
            "date": "/".join([DATES[0], DATES[-1]]),
            "time": [TIME],
            "area": AREA,
            "data_format": "netcdf",
        },
        str(CACHE),
    )
    return xr.open_dataset(CACHE)


def air_density(scene: xr.Dataset) -> xr.DataArray:
    """Density of dry air at the level of the scene, in kg m-3."""
    return PASCALS / (DRY_AIR_GAS_CONSTANT * scene["t"].squeeze(drop=True))


def mass_concentrations(scene: xr.Dataset) -> dict:
    """
    Turn the mixing ratios into concentrations, in kg m-3.

    A mixing ratio is a mass of aerosol per mass of air, so the air density is
    what turns it into the concentration a mixture weights itself with.
    """
    density = air_density(scene)
    return {
        specie: sum(scene[name].squeeze(drop=True) for name in bins) * density
        for specie, bins in SPECIES.items()
    }


def compute(scene: xr.Dataset) -> xr.Dataset:
    """Compute the mixture over every pixel of every date."""
    humidity = scene["r"].squeeze(drop=True)
    # Rounding to the percent collapses the field onto far fewer distinct
    # values, and the tables are given at that resolution anyway.
    humidity = humidity.clip(0.0, MAX_TABULATED_RH).round()

    masses = mass_concentrations(scene)
    print(f"scene      : {dict(humidity.sizes)}")
    print(f"pixels     : {humidity.size:,}")
    distinct = len(np.unique(humidity.values))
    print(f"humidities : {distinct} distinct")
    print(
        f"runs       : {distinct * len(SPECIES)} instead of "
        f"{humidity.size * len(SPECIES):,}"
    )
    for specie, mass in masses.items():
        print(
            f"  {specie.name:<10} {float(mass.min()):.2e} to "
            f"{float(mass.max()):.2e} kg m-3"
        )

    mix = pm.Mix.from_mass(masses, rh_ref=0.0)
    started = time.monotonic()
    result = mix.compute(wl=WAVELENGTHS, rh=humidity, quiet=True)
    print(f"elapsed    : {time.monotonic() - started:.1f} s")
    return result


def plot(kext: xr.DataArray, scene: xr.Dataset, path: Path) -> None:
    """One map per date, on a shared scale."""
    import matplotlib.pyplot as plt

    steps = kext.sizes["valid_time"]
    scale = {"vmin": float(kext.min()), "vmax": float(kext.max())}

    fig, axes = plt.subplots(
        1, steps, figsize=(4.6 * steps, 4.4), constrained_layout=True
    )
    for axis, step in zip(np.atleast_1d(axes), range(steps)):
        image = kext.isel(valid_time=step)
        drawn = axis.pcolormesh(
            image["longitude"], image["latitude"], image, **scale
        )
        axis.set_title(str(image["valid_time"].values)[:10])
        axis.set_xlabel("longitude")
    np.atleast_1d(axes)[0].set_ylabel("latitude")
    fig.colorbar(drawn, ax=axes, label="kext [m-1]")
    fig.suptitle(
        f"CAMS sulphate and sea salt, extinction at "
        f"{TARGET_UM * 1e3:.0f} nm, {PRESSURE_LEVEL} hPa"
    )
    fig.savefig(path, dpi=140)
    print(f"wrote      : {path}")


def main() -> None:
    scene = download()
    result = compute(scene)
    kext = result["kext"].interp(wl=TARGET_UM)
    print(f"kext range : {float(kext.min()):.3e} to {float(kext.max()):.3e}")
    plot(kext, scene, OUTPUT)


if __name__ == "__main__":
    main()
