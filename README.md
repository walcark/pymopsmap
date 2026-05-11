# PyMopsmap

<p align="center">
  <img src="https://github.com/walcark/pymopsmap/actions/workflows/ci.yml/badge.svg">
  <a href="https://pixi.sh"><img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json"></a>
  <a href="https://github.com/astral-sh/ruff"><img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json"></a>
  <a href="https://pypi.org/project/pymopsmap/"><img src="https://img.shields.io/pypi/v/pymopsmap.svg"></a>
  <a href="https://github.com/walcark/pymopsmap/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-Apache--2.0-green"></a>
  <img src="https://img.shields.io/badge/python-3.11%2B-blue">
</p>

Python wrapper for [MOPSMAP](https://mopsmap.net) — aerosol optical property computation based on Mie, T-matrix, and DDA single-particle scattering.

## Overview

PyMopsmap drives the MOPSMAP Fortran binary from Python. Given a set of aerosol microphysical parameters (shape, size distribution, refractive index), it:

1. resolves and downloads the required optical dataset files,
2. writes a MOPSMAP launch file and runs the binary,
3. parses the outputs into an `xarray`-backed `OptiProps` object,
4. caches the result on disk keyed by a blake2b hash of the inputs.

## Requirements

- Python ≥ 3.11 (managed via [pixi](https://prefix.dev/))
- MOPSMAP optical dataset — set `PYMOPSMAP_DATASET_SOURCE` to a local path or HTTP base URL

## Installation

```bash
pixi install
```

## Quick start

```python
import pymopsmap as pm

mp = pm.MicroParameters(
    wavelength=[0.44, 0.55, 0.67],
    n_real=[1.45, 1.45, 1.45],
    n_imag=[1e-3, 1e-3, 1e-3],
    shape=pm.Sphere(),
    psd=pm.LognormalPSD(rm=0.1, sigma=1.5, n=1.0, rmin=0.01, rmax=10.0),
)

op = pm.compute(mp)     # → OptiProps (xarray Dataset)
kext = pm.kext(mp)      # → DataArray indexed by wavelength
```

### Batch computation over external parameters

```python
sweep = pm.ParametricSweep()
for rh, mp_rh in zip([0, 50, 80], [mp_rh0, mp_rh50, mp_rh80]):
    sweep.add(pm.ParticleMixture([mp_rh]), {"rh": rh})

op = pm.compute(sweep)   # → OptiProps with an extra 'rh' dimension
```

### CAMS aerosol adapter

```python
from pymopsmap.adapters import cams_to_kext, CamsAerosol, CamsVersion

kext = cams_to_kext(
    aerosol=CamsAerosol.SEA_SALT_CAMS,
    version=CamsVersion.V49_R1,
    wl_microns=[0.44, 0.55, 0.67],
    rh=[0, 50, 80, 99],
)
```

### OPAC aerosol adapter

```python
from pymopsmap.adapters.input.opac import OpacMix, OpacMixName, OpacHumidityMode

mix = OpacMix(OpacMixName.CONTINENTAL_AVERAGE)

# GEISA mode: wet PSD and refractive index interpolated from GEISA tables
op = mix.compute(wavelengths=[0.44, 0.55, 0.67], rhs=[0, 50, 80])

# Kappa mode: hygroscopic growth via κ parameterisation (Zieger et al. 2013)
# with volume-weighted refractive index mixing with water
op = mix.compute(
    wavelengths=[0.44, 0.55, 0.67],
    rhs=[0, 50, 80],
    mode=OpacHumidityMode.KAPPA,
)
```

## Output types

```python
op = pm.compute(mp, output_types=frozenset({
    pm.OutputType.INTEGRATED,
    pm.OutputType.PHASE_FUNCTION,
    pm.OutputType.LIDAR,
}))
```

Available: `INTEGRATED`, `LIDAR`, `PHASE_FUNCTION`, `SCATTERING_MATRIX`, `VOLUME_SCATTERING_FUNCTION`, `COEFF`.

## Dataset cache

Files are cached under `~/.cache/pymopsmap/` by default.

```bash
export PYMOPSMAP_DATASET_SOURCE=/path/to/optical_dataset
export PYMOPSMAP_CACHE_DIR=/custom/cache/dir   # optional
```

```python
pm.cache_status(mp)   # lists cached vs missing dataset files
pm.prefetch(mp)       # download without computing
```

## Project structure

```
src/pymopsmap/
├── models/      # MicroParameters, OptiProps, OutputRequest, particle systems
├── engine/      # MOPSMAP binary interface (launch file, runner, output parser)
├── cache/       # Optical dataset files and result cache
├── adapters/
│   ├── input/   # External data formats → MicroParameters  (e.g. CAMS, OPAC)
│   └── output/  # OptiProps → external formats             (e.g. SMART-G)
└── utils/       # Logging, types, temp files, caching
```

## Development

```bash
pixi run -e dev test    # pytest + coverage
pixi run -e dev lint    # ruff
pixi run -e dev all     # fmt + lint + type-check + test
```
