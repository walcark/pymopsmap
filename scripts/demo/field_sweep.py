"""
Compute optical properties over a humidity field that varies in space and time.

    pixi run -e dev python -m scripts.demo.field_sweep

A relative humidity cube of (y, x, t) is handed to ``compute`` as a
``DataArray``. Its dimensions survive into the result, so what comes back is a
cube of optical properties rather than a table to reshape by hand.

The point of the exercise is the cost: the cube holds a million values but only
a handful of distinct ones, and a MOPSMAP run costs seconds. Only the distinct
values are computed, once each, and the results are spread back over the cube.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import xarray as xr

import pymopsmap as pm

SHAPE = (100, 100, 100)  # y, x, t
HUMIDITIES = (0.0, 50.0, 90.0, 95.0)

# Enough wavelengths to interpolate between, few enough to keep the result
# small: the field is a million points, so every wavelength costs 8 MB per
# variable.
WAVELENGTHS = [0.44, 0.55, 0.67, 0.87]
TARGET_UM = 0.66

OUTPUT = Path("field_sweep.png")


def humidity_cube() -> xr.DataArray:
    """
    A moist front advancing west to east over the scene.

    Structured rather than random, so the images show something: dry ahead of
    the front, progressively moister behind it, and each pixel takes one of
    four humidities.
    """
    y, x, t = SHAPE
    columns = np.arange(x)[None, :, None]
    times = np.arange(t)[None, None, :]
    rows = np.arange(y)[:, None, None]

    # The front advances eastward with time and waves with latitude; behind
    # it the air is moist, ahead of it it stays dry.
    front = 0.85 * times + 12.0 * np.sin(2 * np.pi * rows / y)
    step = np.clip((front - columns) / (x / (2 * len(HUMIDITIES))), 0, None)
    index = np.clip(step.astype(int), 0, len(HUMIDITIES) - 1)

    return xr.DataArray(
        np.asarray(HUMIDITIES)[index],
        dims=["y", "x", "t"],
        name="rh",
    )


def compute(field: xr.DataArray) -> xr.Dataset:
    """Run the sweep and report what it actually cost."""
    specie = pm.load(pm.CAMS.SULPHATE)
    distinct = np.unique(field.values)

    print(f"field      : {dict(zip(field.dims, field.shape))}")
    print(f"points     : {field.size:,}")
    print(f"distinct   : {len(distinct)} humidities {list(distinct)}")

    started = time.monotonic()
    result = specie.compute(wl=WAVELENGTHS, rh=field, quiet=True)
    elapsed = time.monotonic() - started

    print(f"elapsed    : {elapsed:.1f} s")
    print(f"result     : {dict(result.sizes)}")
    print(f"in memory  : {result.nbytes / 1e6:.0f} MB")
    return result


# Two moments while the front is still crossing the scene; at the ends it has
# either not arrived or fully passed, and both images would be flat.
SLICES = (25, 50, 70)


def plot(kext: xr.DataArray, path: Path) -> None:
    """Two time slices of the extinction, at the interpolated wavelength."""
    import matplotlib.pyplot as plt

    images = [kext.isel(t=step) for step in SLICES]
    scale = {"vmin": float(kext.min()), "vmax": float(kext.max())}

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    for axis, image, step in zip(axes, images, SLICES):
        when = f"t = {step}"
        drawn = axis.pcolormesh(image["x"], image["y"], image, **scale)
        axis.set_title(when)
        axis.set_xlabel("x")
        axis.set_aspect("equal")
    axes[0].set_ylabel("y")
    fig.colorbar(drawn, ax=axes, label="kext [m-1]")
    fig.suptitle(
        f"CAMS sulphate extinction at {TARGET_UM * 1e3:.0f} nm, "
        "over a moving moist front"
    )
    fig.savefig(path, dpi=140)
    print(f"wrote      : {path}")


def main() -> None:
    result = compute(humidity_cube())
    kext = result["kext"].interp(wl=TARGET_UM)
    print(
        f"kext range : {float(kext.min()):.3e} to {float(kext.max()):.3e} m-1"
    )
    plot(kext, OUTPUT)


if __name__ == "__main__":
    main()
