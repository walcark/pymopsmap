# PyMopsmap

<p align="center">
  <img src="https://github.com/walcark/pymopsmap/actions/workflows/ci.yml/badge.svg">
  <a href="https://pixi.sh"><img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json"></a>
  <a href="https://github.com/astral-sh/ruff"><img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json"></a>
  <a href="https://pypi.org/project/pymopsmap/"><img src="https://img.shields.io/pypi/v/pymopsmap.svg"></a>
  <a href="https://github.com/walcark/pymopsmap/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-Apache--2.0-green"></a>
  <img src="https://img.shields.io/badge/python-3.11%2B-blue">
</p>

Aerosol optical properties from microphysics, on the grid you ask for.

PyMopsmap wraps [MOPSMAP](https://mopsmap.net) (Mie, T-matrix and DDA
single-particle scattering, [Gasteiger & Wiegner 2018,
GMD](https://doi.org/10.5194/gmd-11-2739-2018)) behind two calls: load an
aerosol, compute its optical properties over any set of parameters.

```python
import numpy as np
import pymopsmap as pm

sulphate = pm.load(pm.CAMS.SULPHATE)
op = sulphate.compute(rh=50, wl=np.linspace(0.4, 2.0, 100))

op.kext        # <xarray.DataArray (wl: 100)>
op.ssa
```

## Installation

```bash
pip install pymopsmap
```

The built-in aerosol catalogue (CAMS, OPAC) ships with the package. Two external
pieces are not bundled:

| Piece | How to get it |
|---|---|
| MOPSMAP binary | download from [mopsmap.net](https://mopsmap.net), place at `bin/mopsmap/mopsmap` |
| Optical dataset | set `PYMOPSMAP_DATASET_SOURCE` to a local path or HTTP base URL; files are fetched on demand into `~/.cache/pymopsmap/` |

```bash
export PYMOPSMAP_DATASET_SOURCE=https://your-server.org/mopsmap_dataset
python -c "import pymopsmap as pm; print(pm.load(pm.CAMS.SULPHATE))"
```

For development:

```bash
git clone https://github.com/walcark/pymopsmap.git
cd pymopsmap
pixi install -e dev
```

## Sweeps

`compute` returns an `xarray.Dataset` whose dimensions are the ones you asked
for. Pass a scalar and you get a scalar axis; pass a `DataArray` and you get
your own dimension name.

```python
op = sulphate.compute(rh=50, wl=wl)                       # (wl,)
op = sulphate.compute(rh=[50, 70, 90], wl=wl)             # (rh, wl)

op = sulphate.compute(
    rh=xr.DataArray([50, 70, 90], dims="rh_nominal"),
    wl=wl,
)                                                          # (rh_nominal, wl)
```

Results are stored on disk by [xsweep](https://github.com/walcark/xsweep), so
re-running a grid you have already computed calls MOPSMAP for nothing, and a
sweep interrupted halfway resumes from the points it is missing.

## Custom species

The parameter space lives on the species, not in the call. Any microphysical
parameter accepts a `DataArray`, and xarray's broadcasting rules apply: distinct
dimensions multiply, a shared dimension varies together.

```python
aer = pm.Specie.custom(
    pm.Mode(
        shape=pm.shapes.Sphere(),
        psd=pm.psd.LognormalPSD(
            rm=0.05, sigma=1.5, n=1e9, rmin=0.01, rmax=10.0
        ),
        n_real=1.45,
        n_imag=xr.DataArray([1e-4, 1e-3, 1e-2], dims="absorption"),
        density_dry=1.8,
        sweep={"rm": xr.DataArray(np.linspace(0.05, 0.5, 20), dims="rm")},
    )
)

aer.swept                        # {"rm": 20, "absorption": 3}
aer.at(wl=wl, rm=7, absorption=0)
```

Size distribution and shape parameters are swept through `sweep` rather than
inside the `psd` itself, whose fields stay typed as floats. That is deliberate:
every materialised point is validated, and a `MicroParameters` never holds an
array where a number belongs.

Put two parameters on the **same** dimension to walk a trajectory rather than a
grid:

```python
aging = np.linspace(0.05, 0.5, 20)
sweep = {
    "rm": xr.DataArray(aging, dims="aging"),
    "sigma": xr.DataArray(np.linspace(1.4, 2.1, 20), dims="aging"),
}                                # 20 points, not 400
```

A custom species saves and reloads through the same format as the built-in
catalogue:

```python
aer.to_netcdf("my_aerosol.nc")
pm.load("my_aerosol.nc").compute(wl=wl)
```

## Mixtures

`Mix` combines species as an external mixture.

```python
mix = pm.Mix({
    pm.CAMS.SULPHATE: 3.2e9,      # m-3
    pm.CAMS.SEA_SALT: 1.1e8,
})
op = mix.compute(rh=[50, 80], wl=wl)
```

Each species is computed once and the results are combined afterwards, so
changing a weight costs nothing:

```python
polluted = pm.Mix({pm.CAMS.SULPHATE: 9.6e9, pm.CAMS.SEA_SALT: 1.1e8})
polluted.compute(rh=[50, 80], wl=wl)      # no MOPSMAP run, reuses the store
```

### Weighting by mass or by optical depth

Number concentration is rarely what you have. Two other currencies express the
same weights:

```python
# CAMS gives mass mixing ratios, not number densities
mix = pm.Mix.from_mass({pm.CAMS.SULPHATE: 4.1e-9, pm.CAMS.DUST: 2.2e-8})  # kg m-3

# AERONET-style: fractions of total optical depth at a reference wavelength
mix = pm.Mix.from_optical_depth(
    {pm.CAMS.SULPHATE: 0.30, pm.CAMS.DUST: 0.70},
    wl_ref=0.56,
    rh_ref=50,
)
op = mix.compute(rh=50, wl=wl)
op.kext            # normalised extinction, sums to 1 at wl_ref
op.kext * 0.25     # for a total AOD of 0.25

mix.weights        # {CAMS.SULPHATE: 2.9e9, CAMS.DUST: 4.1e7}  in m-3
```

`rh_ref` is required, and it is the humidity **at which your fractions were
observed**, not a convention. Optical depth fractions are a derived quantity:
they depend on humidity, because species do not grow alike. Number
concentration is what stays fixed, so PyMopsmap inverts your fractions into
concentrations at `rh_ref` and carries those. Ask for another humidity
afterwards and you get what the same air mass would look like there, with
different fractions.

If your fractions and your humidity come from the same observation, `rh_ref`
equals `rh` and you type the number twice.

### OPAC climatology

```python
op = pm.opac_mix("continental_average").compute(rh=[0, 50, 80], wl=wl)
```

The ten climatologies of Hess et al. (1998) ship with the package, as does the
dry state of their ten components under `pm.OPAC`.

## Humidity

How a species responds to relative humidity is a property of the species, stated
in its data file, not a switch at call time.

| `specie.growth` | The data provides | Effect |
|---|---|---|
| `"tabulated"` | wet values on a real `rh` axis | interpolated at your `rh` |
| `"kappa"` | dry values and a hygroscopicity `kappa` | MOPSMAP grows the particles |
| `"none"` | dry only | passing `rh=` raises |

```python
sulphate.growth              # 'tabulated'
sulphate.rh_range            # (0.0, 95.0)
sulphate.compute(rh=99, wl=wl)
# DomainError: rh=99 is outside the CAMS 49r1 range [0, 95] for 'sulphate'
```

For a `"kappa"` species the file value is a default you can override:

```python
aer.compute(rh=80, wl=wl, kappa=0.3)
```

## Output types

Integrated properties come back by default. Ask for more:

```python
op = sulphate.compute(
    rh=50, wl=wl,
    outputs={pm.OutputType.PHASE_FUNCTION, pm.OutputType.LIDAR},
)
op.phase          # (wl, theta)
op.lidar_ratio    # (wl,)
```

Available: `INTEGRATED`, `LIDAR`, `PHASE_FUNCTION`, `SCATTERING_MATRIX`,
`VOLUME_SCATTERING_FUNCTION`, `COEFF`.

## Exports

The result is a plain `xarray.Dataset`, so the whole xarray API applies. Format
conversions live on a `.mopsmap` accessor:

```python
op.sel(wl=0.55, method="nearest").kext
op.mopsmap.to_smartg("lut.nc", humidity_dim="rh")
```

## Catalogue

```python
list(pm.CAMS)                # SULPHATE, SEA_SALT, DUST, BLACK_CARBON, ...
pm.CAMS.SULPHATE.versions    # ('47r1', '48r1', '49r1')

pm.load(pm.CAMS.SULPHATE, version="48r1")
pm.load(pm.OPAC.WASO)        # OPAC components keep their published codes
```

A species absent from a version is absent from its file, so loading it fails
clearly rather than yielding NaN: `secondary_organic` did not exist in CAMS
47r1.

Dataset files needed for a computation can be inspected and pre-fetched:

```python
sulphate.cache_status(wl=wl, rh=50)     # cached vs missing
sulphate.prefetch(wl=wl, rh=50)         # download without computing
```

Two errors say what a computation cannot do. `DomainError` when a wavelength or
a humidity falls outside what the species tabulates, and `CoverageError` when
the optical dataset at hand does not ship a file the refractive index needs.
The MOPSMAP dataset comes in two archives, and the main one covers refractive
indices from 1.28 to 1.64: soot and mineral dust reach past it.

## How it works

Three levels, and only the first two are public.

| Level | Object | Role |
|---|---|---|
| Description | `Specie` | the parameter space, from a NetCDF file or built in Python |
| Point | `MicroParameters` | one concrete, validated point, rendered to MOPSMAP commands |
| Combination | `Mix` | concentration-weighted species, external mixing |

A `Specie` wraps an `xarray.DataTree`: one group per mode, refractive index and
size distribution parameters as variables, the distribution and shape types as
attributes. Reading a species and writing one are the same code path, which is
why a hand-built aerosol round-trips through the catalogue format.

`compute` hands that space to [xsweep](https://github.com/walcark/xsweep),
which walks the points; each one materialises a `MicroParameters`, runs MOPSMAP
in its own directory, and the results are assembled back onto your dimensions.
A point that fails raises rather than becoming a silent NaN.

```
src/pymopsmap/
├── species/      # Specie, Mode, Mix, the NetCDF schema, the catalogue
├── shapes.py     # Sphere, Spheroid, Irregular, ...
├── psd.py        # LognormalPSD, ModifiedGammaPSD, ...
├── microparams.py# one validated point
├── engine/       # one point, one MOPSMAP run, and the output rules
├── sweep.py      # parameter space to points, and the xsweep binding
├── scatlib/      # the MOPSMAP optical dataset: resolve, download, cache
├── accessors.py  # .mopsmap accessor on the result
└── data/         # the built-in catalogue
```

## Adding a species

Write the NetCDF, load it. No code:

```python
pm.load("my_source/my_specie.nc")
```

The schema is documented in `docs/api-v2-spec.md`, section 4. To contribute a
source to the built-in catalogue, add an ingestion script under
`scripts/build_catalog/` that emits that schema.

## Roadmap

- **Growable sweeps**: a store holds one grid, so widening a request today
  recomputes it rather than extending what is already there.
- **Transparent remote dataset**: automatic download of the optical dataset
  when `PYMOPSMAP_DATASET_SOURCE` is not set, removing the manual setup step.
- **Article validation**: `scripts/validation/` reproduces Figure 5 of
  [Gasteiger & Wiegner (2018)](https://doi.org/10.5194/gmd-11-2739-2018); the
  other figures are next.
- **Tabulated OPAC**: the wet state published by GEISA, alongside the kappa
  flavour that ships today. The GEISA file host was decommissioned during the
  migration of the database, so the links on its pages no longer resolve.

## Development

```bash
pixi run -e dev test          # pytest + coverage
pixi run -e dev lint          # ruff
pixi run -e dev all           # fmt + lint + type-check + test
```
