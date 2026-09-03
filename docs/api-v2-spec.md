# PyMopsmap API v2 - Design Specification

Status: draft, pass 2 (revised after a full read of the current codebase)
Goal: one declarative aerosol description, one computation entry point, and a
public surface small enough to fit on one screen.

---

## 1. Target API

```python
import pymopsmap as pm

specie = pm.load(pm.CAMS.SULPHATE)                       # reads one NetCDF group
op = specie.compute(rh=50, wl=np.linspace(0.4, 2.0, 100))

op = specie.compute(
    rh=xr.DataArray([50, 70, 90], dims="rh_nom_dif"),
    wl=np.linspace(0.4, 2.0, 100),
)                                                        # dims (rh_nom_dif, wl)
```

Everything else (per-source adapters, sweep bookkeeping, result caching,
launch-file generation) becomes an implementation detail.

---

## 2. Problems in the current design

### 2.1 Structural

| Problem | Location |
|---|---|
| Each input adapter re-implements interpolation, sweep construction and the engine call | `adapters/input/cams.py`, `adapters/input/opac.py` |
| `ParametricSweep` can only express a list of pre-built mixtures with an attached param dict: no cartesian product, no correlated walk | `models/particle_mixture.py` |
| Every source invents its own `compute` signature | `OpacMix.compute` vs `cams_to_kext` vs `pm.compute` |
| Number concentration lives in a JSON (CAMS) or a Python dict (OPAC), decoupled from the data it applies to | `cams.py:88`, `opac.py:25` |
| `rmin`/`rmax` hardcoded for CAMS (`0.001`, `40.0`), in-file for OPAC | `cams.py:78` |
| Result cache (blake2b) duplicates what a sweep store provides, without resumability | `cache/results.py` |
| 30 names exported from `pymopsmap/__init__.py` | `__init__.py:__all__` |

### 2.2 Data-level inconsistencies

| Quantity | CAMS file | OPAC file | Fix |
|---|---|---|---|
| Wavelength | nm, `int32` (200..2350) | um | store um everywhere |
| Width | `lnvar_*` = ln(sigma), `exp()` at call site | `sigma` direct | store `sigma` |
| Imaginary index | negative (`-mi_f` at call site) | negative (`-n_imag` at call site) | store MOPSMAP positive convention |
| `n` | absent, JSON sidecar | absent, Python dict | see section 6 |
| Version | `cams_versions` coordinate | n/a | one file per version |

### 2.3 Three competing humidity mechanisms (most important design finding)

MOPSMAP has **native** hygroscopic growth: `launch_file.py` emits a global
`rH <value>`, and `commands.py` emits `mode N kappa <k>` and
`mode N density <d>`. But:

- `cams.py` interpolates the wet PSD and wet refractive index itself and passes
  `rh=None` to the engine.
- `opac.py` GEISA mode does the same.
- `opac.py` KAPPA mode applies the growth formula **and** the volume-weighted
  water mixing by hand (`_build_mode_kappa`, `_water_refrac`), explicitly
  "without needing kappa/rH commands".

Consequence: `MicroParameters.kappa` and the `rh=` argument of `pm.compute` are
dead paths for every shipped adapter, and the word `rh` means two different
things depending on the call site.

**v2 decision (settled).** There is exactly one implementation of hygroscopic
growth, and it is never a Python formula. The `growth` attribute of the stored
group says which of two physical descriptions the file carries:

| `growth` | The file contains | Who applies it |
|---|---|---|
| `"tabulated"` | a real `rh` axis, wet values tabulated | interpolation in `load` |
| `"kappa"` | no `rh` axis, dry values plus a `kappa` variable | MOPSMAP (`mode N kappa` + `rH`) |
| `"none"` | neither an `rh` axis nor a `kappa` variable | nobody, `rh=` raises |

Notes:

- The attribute names the **physics**, not the executor. Once the Python formula
  is gone there is a single possible backend for the kappa mode, so calling the
  value `"mopsmap"` would carry no information and would have to change if the
  backend ever did.
- A `kappa` value in the file and a `kappa` passed at call time are not two
  modes, only two provenances of the same number. The file variable is the
  default, `compute(kappa=...)` overrides it. One code path.
- `"none"` earns an explicit value: without it a dry-only specie reports
  "missing attribute" instead of "this specie does not model humidity".
- `growth` is checkable, and `load` asserts it: `tabulated` implies an `rh` dim;
  `kappa` implies a `kappa` variable and **no** `rh` dim (otherwise the source of
  truth is ambiguous); `none` implies neither.

`_build_mode_kappa` and `_water_refrac` are deleted. `MicroParameters.kappa` and
`MicroParameters.density` become load-bearing again, serving `growth="kappa"`.
`pm.compute(..., rh=...)` at module level is removed.

### 2.3.1 What the MOPSMAP source settles

`bin/mopsmap/src/calc_hygroscopic_growth.f90` removes every doubt about the
delegation, and adds three constraints.

**The Python port is line-for-line identical to the Fortran.** No regression
risk, the delegation is numerically exact:

| Quantity | MOPSMAP | `opac.py:_build_mode_kappa` |
|---|---|---|
| Growth factor | `(1 + kappa*rh/(100-rh))**(1/3)` | same |
| Refractive index | `(m + m_water*(gf**3 - 1)) / gf**3` | `(1-f_w)*m + f_w*m_water` with `f_w=(gf**3-1)/gf**3`, so `1-f_w = 1/gf**3`: identical |
| Density | `(density + 1.00*(gf**3 - 1)) / gf**3` | `(1-f_w)*rho_dry + f_w*1.0`: identical |

**The delegation also fixes a bug.** MOPSMAP grows both truncation bounds:

```fortran
mode(i_mode)%rmin = mode(i_mode)%rmin * growth_factor
mode(i_mode)%rmax = mode(i_mode)%rmax * growth_factor
```

`opac.py` grows only one (`rmin=... * gf`, `rmax=float(sel_dry["rmax"])`). At
high RH on a strongly hygroscopic species (`sscm`, kappa=0.916, gf**3 near 92 at
99% RH) the upper cut stays dry and truncates the wet distribution.

**Three engine constraints, to enforce at load time:**

| Constraint | Source | Consequence |
|---|---|---|
| When `rH > 0`, **every** mode of the launch file must carry a kappa, or MOPSMAP aborts | `if(mode(i_mode)%kappa.lt.0) ... stop` | `growth` must be uniform across the modes of one run, so it belongs on the specie group, not the mode. A heterogeneous `Mix` cannot go through a single launch file. Settles open question 10.4. |
| `mod_gamma` plus growth is not implemented | explicit `stop` in the same file | reject `growth="kappa"` with `psd_type="mod-gamma"` at load |
| RH is capped at 99 | `if(relative_humidity.gt.99) ... stop` | one more validation bound, engine side |

**Default density is 2.6, not 1.0.** `module_input.f90:61` declares
`real(dp) :: density = 2.6`, the mineral dust value. `cams.py` never sets a
density, so every CAMS specie (sulphate, sea salt, black carbon, ...) currently
runs at 2.6 g cm-3, making their `mass_conc` and `ext_to_mass` outputs wrong.

### 2.3.2 The `density` trap, settled

`density` meaning dry or wet depending on another attribute is the same disease
as `lnvar` versus `sigma`. The section 4.1 rule applies: the meaning goes in the
variable name.

**Store `density_dry` only. Never store a wet density.** Dry density is a
material property and is invariant; wet density is derived, by the same formula
in both cases:

| `growth` | What `load` does |
|---|---|
| `kappa` | pass `density_dry` straight through as `mode N density`; MOPSMAP wets it itself |
| `tabulated` | compute `gf**3 = (rm(rh)/rm(0))**3` from the table, apply `(rho_dry + gf**3 - 1)/gf**3`, send the result |

Three lines, moved out of `opac.py:98`, and guaranteed consistent because it is
the formula MOPSMAP would have applied. `density_dry` is required for every
group, which also fixes the CAMS default-2.6 problem above.

### 2.4 Verified bug: unvalidated interpolation domain

The CAMS `relative_humidity` grid is `0..95` step 5. The README example uses
`rh=[0, 50, 80, 99]`. `ds.interp(relative_humidity=99)` returns NaN
(no extrapolation), and the wavelength grid stops at 2350 nm with the same
behaviour. Verified:

```
rmodal_f @ rh=99         -> nan
mr_f @ rh=99, wl=550nm   -> nan
mr_f @ rh=50, wl=2500nm  -> nan
```

Downstream, the NaN splits in two:

- `rm=nan` is rejected by pydantic, but with a misleading message
  (`rm: Input should be greater than 0`) that names the wrong cause.
- `n_real=nan` / `n_imag=nan` are **accepted**: `assert_strictly_positive` tests
  `v <= 0`, and `nan <= 0` is `False`. The NaN then reaches
  `resolver._bracket(grid, nan)`, where `np.searchsorted` returns `len(grid)`,
  and the function silently returns the **last grid point**. Verified:
  `_bracket(linspace(1,2,11), nan) -> [2.0]`.

So an out-of-domain RH currently selects a wrong optical dataset file without
any warning.

**v2 rules.**
1. `load` records per-axis source bounds on the `Specie`.
2. `compute` validates every target against them **before** materialising any
   `MicroParameters`, and raises naming the axis, the requested value and the
   available range.
3. `Float64List` gains an `assert_finite` validator, so NaN can never reach the
   resolver regardless of the path taken.

### 2.5 Verified bug: coverage clipping breaks non-scalar outputs

`clip_modes_to_coverage` drops out-of-range wavelengths, then `engine/__init__.py`
reindexes with `full = np.full(len(wl_full), nan); full[valid_mask] = da.values`
and rebuilds every variable as `(("wl",), full)`. This assumes all variables are
1-D over `wl`. `phase(wl, theta)`, `scattering_matrix(wl, theta, element)` and
`coeff(wl, l, coeff_element)` are not, so clipping plus any ASCII output type
fails on shape mismatch.

**v2 rule.** Reindex with `Dataset.reindex(wl=...)` on the real wl coordinate,
not by positional assignment, so extra dims are preserved.

Note also that the clip mask depends on `rmax`, which varies with `rh` (growth)
and with any swept `rm`. The output `wl` axis is fixed across the sweep, but the
NaN pattern is per-point. This must be documented, not hidden.

### 2.6 Verified bug: all modes share one refractive index file

Most severe finding of the whole review. `utils/temp.py:get_tempfile` returns a **stable path per
filename**, held in a process-global registry:

```python
def get_tempfile(filename: str) -> Path:
    if filename not in TEMP_REGISTRY:
        TEMP_REGISTRY[filename] = get_tempdir() / filename
    return TEMP_REGISTRY[filename]
```

`write_refr_file` always asks for `"ri.txt"`. It is called once per mode, opens
the path in `"w"` mode, and returns the same string every time. So in a
multi-mode mixture each mode overwrites the previous one and every
`refrac file` command points at the last mode written. Verified:

```
mode 1 refrac file '/tmp/pymopsmap-b5e4.../ri.txt'
mode 2 refrac file '/tmp/pymopsmap-b5e4.../ri.txt'      <- same file

actual content:
0.440000 1.600000 0.050000       <- the COARSE mode index
0.550000 1.600000 0.050000
0.670000 1.600000 0.050000
```

The fine mode was built with `n_real=1.40, n_imag=1e-3`. Those values never
reach MOPSMAP.

Impact: **every CAMS specie (fine + coarse) and every OPAC mix is computed with
the last mode's refractive index applied to all modes.** Single-mode runs are
unaffected, which is why the test suite misses it.

This also blocks xsweep: two points running in parallel would share
`mopsmap.txt` and `output.nc` as well.

**v2 fix.** `engine/workspace.py`: one isolated directory per run, with per-mode
files `ri_1.txt`, `ri_2.txt`, ... `utils/temp.py` is deleted. Regression test:
assert that a two-mode launch file names two distinct `refrac file` paths.

### 2.7 Smaller findings

| Finding | Location |
|---|---|
| `cache_key` uses `pickle.dumps` of a pydantic model: not stable across Python/pydantic versions. Use canonical `model_dump_json()`. | `utils/caching.py` |
| `n_angles=2000` is hardcoded, never surfaced, and absent from the cache key | `launch_file.py:36` |
| `ResultCache.get` returns a lazily-opened `Dataset` and never closes the handle: fd leak on repeated hits | `cache/results.py:get` |
| `DATA_PATH` points outside the package (`ROOT_PATH/data`), and hatchling ships no `data/`: an installed wheel has no CAMS file | `utils/__init__.py`, `pyproject.toml` |
| `extend_optiprops` leaks the loop variable `k` into `swap_dims`, and `unstack` fills NaN on a partial grid | `models/optiprops.py` |
| `output_request_commands` is fully dead code, duplicated by `_file_suffix` | `engine/commands.py` |
| `format_netcdf_file` is dead legacy code | `engine/output_format.py` |
| `smartg.py` hardcodes `op.coord("rh")` and `phase.ndim == 4`, so it breaks as soon as the humidity dim is user-named | `adapters/output/smartg.py` |
| `OptiProps.sel(prop)` shadows `xr.Dataset.sel` with different semantics | `models/optiprops.py` |

---

## 3. Three levels

| Level | Object | Role |
|---|---|---|
| Description | `Specie` (wraps an `xr.DataTree`) | carries the parameter space; loaded from `.nc` or built in Python |
| Point | `MicroParameters` (pydantic) | one concrete point, validated, rendered to MOPSMAP commands |
| Combination | `Mix` | concentration-weighted set of species (external mixing) |

`MicroParameters` is not removed. It is demoted from public entry point to
per-point internal representation, and kept as the escape hatch for a single
computation.

---

## 4. Storage format

One NetCDF file per source *version*. Groups are species, subgroups are modes.

```
cams/49r1.nc
+-- /                    attrs: source="CAMS", version="49r1", schema_rev=1
+-- /sulphate            attrs: growth="tabulated"
    +-- /fine            attrs: psd_type="lognormal", shape_type="sphere"
    |     dims: (wl, rh)
    |     n_real(wl, rh)     units="1"
    |     n_imag(wl, rh)     units="1", sign="positive"
    |     rm(rh)             units="um"
    |     sigma(rh)          units="1"      (sigma, not ln(sigma))
    |     rmin(), rmax()     units="um"
    |     density_dry()      units="g cm-3", always dry (section 2.3.2)
    |     n()                units="m-3", required: carries the fine/coarse ratio
    |     kappa()            present iff growth="kappa"
    +-- /coarse          same layout
```

Rules:

- `attrs` carry **only** discriminator strings, reusing the existing pydantic
  `Literal` values (`psd_type`, `shape_type`). Numeric shape parameters
  (`aspect_ratio`, `zeta1`, `sigma_ar`) are variables, so they may depend on `rh`.
- `growth` per specie group, values from section 2.3.
- Version is a file, not a coordinate.
- Root `schema_rev` is an integer, bumped on any breaking layout change, and
  feeds the sweep store version string (section 7).
- `density_dry` is required for every mode and is always the dry material
  density. Wet density is derived at load, never stored (section 2.3.2).

### 4.1 Conventions: declared, never dispatched (settled)

The file states its conventions, but the reader validates them rather than
adapting to them. Conversion logic lives in the writer, once per source, not in
the reader. Three kinds of convention, three different treatments:

| Kind | Rule | Example |
|---|---|---|
| **Unit** | CF-style `units` attribute, declared **and honoured** on read (converted) | `units="um"`. Cheap, standard, a finite conversion table, and a third-party tool opening the file reads it correctly. |
| **Meaning** | Encoded in the **variable name**, never in an attribute | `sigma` and `ln_sigma` are different variables. A reader that does not know `ln_sigma` fails loudly instead of misreading it as `sigma`. |
| **Sign** | Normalised at write time; the attribute documents **and is asserted** | `n_imag` always positive, `sign="positive"`, and `load` raises if `(n_imag < 0).any()`. The assertion catches a badly written file, which an attribute alone would not. |

Canonical choices: `wl` and all radii in micrometres, `n` in m-3, `rh` in
percent, `density` in g cm-3, `n_imag` positive.

The current CAMS `lnvar_f` is the worked example: it is not a unit but a
different quantity. The fix is not an attribute switch, it is to write `sigma`
into the file and do the `exp()` in the writer. Same for `mi_f`, whose sign flip
moves from `cams.py:187` into the writer.

`load` is then a thin walk over the tree, reusing the discriminated unions
unchanged:

```python
psd = {"type": node.attrs["psd_type"], **psd_scalars}
shape = {"type": node.attrs["shape_type"], **shape_scalars}
MicroParameters(wavelength=wl, n_real=..., n_imag=..., shape=shape, psd=psd)
```

### 4.2 Coverage metadata

Two independent clipping mechanisms exist and must report through one channel:

1. **Source coverage**: the `wl` and `rh` bounds of the stored group
   (section 2.4). Violation is an error, raised before any run.
2. **Dataset coverage**: the size-parameter limits of the MOPSMAP optical
   dataset (`coverage.py`, Gasteiger & Wiegner 2018 Tables 1-2). Violation is a
   warning plus NaN, because it depends on `rmax` and cannot be known up front.

---

## 5. Custom species

The sweep dimension is carried by the `Specie`, not by `compute` kwargs. This
avoids an untypable `compute(**kwargs)` where `rm=` would be ambiguous across
modes, and inherits xarray broadcasting semantics for free.

```python
aer = pm.Specie.custom(
    pm.Mode(
        shape=pm.Sphere(),
        psd=pm.Lognormal(
            rm=xr.DataArray(np.linspace(0.05, 0.5, 20), dims="rm"),
            sigma=1.5, rmin=0.01, rmax=10.0,
        ),
        n_real=1.45,
        n_imag=xr.DataArray([1e-4, 1e-3, 1e-2], dims="k"),
    )
)
op = aer.compute(wl=np.linspace(0.4, 2.0, 100))   # dims (rm, k, wl)
```

Distinct dims multiply (cartesian product); a shared dim varies together:

```python
# rm and sigma on the SAME dim -> a 20-point trajectory, not a 20x20 grid
rm    = xr.DataArray(np.linspace(0.05, 0.5, 20), dims="aging")
sigma = xr.DataArray(np.linspace(1.4, 2.1, 20), dims="aging")
```

Implementation constraint: pydantic fields stay typed `float`. `Specie.custom`
accepts `DataArray`, stores them in the `DataTree`, and the sweep engine hands
scalars back point by point. `MicroParameters` never sees a `DataArray`.

Round trip closes the loop and is the primary format test:

```python
aer.to_netcdf("my_aerosol.nc")
pm.load("my_aerosol.nc").compute(...)   # identical result
```

A custom `Specie` built without an `rh` axis gets `growth="none"`, so passing
`rh=` raises with a message naming the specie. See open question 10.1.

---

## 6. Mixtures

`Mix` is a concentration-weighted set of species, external mixing.

```python
mix = pm.Mix({pm.CAMS.SULPHATE: 3.2e9, pm.CAMS.SEA_SALT: 1.1e8})  # m-3
op = mix.compute(rh=[50, 80], wl=np.linspace(0.4, 2.0, 100))
```

`Mix.compute` computes each `Specie` **once** and combines afterwards, rather
than flattening all modes into a single launch file. External mixing is
additive, so the result is identical to what MOPSMAP computes internally, but
any change of concentration is then free: N runs total instead of one run per
concentration combination.

Combination rules, by variable class:

| Class | Variables | Rule |
|---|---|---|
| Additive | `kext`, `ksca`, `backscatter`, `n`, `mass_conc`, `vol_dens`, `cross_dens` | direct sum |
| ksca-weighted | `g`, `phase`, `scattering_matrix`, `vol_sca_func`, `coeff` | `sum(ksca_i * X_i) / sum(ksca_i)` |
| Ratio | `ssa` | `sum(ksca_i) / sum(kext_i)` |
| Derived from moments | `reff` | recompute from summed `vol_dens` and `cross_dens`, never average |
| Spectral derivative | `angstrom_ext`, `angstrom_sca`, `angstrom_abs`, `angstrom_back`, `lidar_ratio`, `depol_ratio` | recompute from the combined spectra, never combine directly |

The last two rows are the ones that silently give wrong numbers if treated as
additive. They must be recomputed, and the combination code must fail loudly on
any output variable it does not have a rule for.

### 6.1 How the rules are wired

Not one method per output variable: that would be roughly 25 methods, and adding
a variable to `engine/outputs.py` would silently require remembering to add the
matching method in `mix.py`. Nobody remembers.

**Four rule functions plus a lookup table, declared next to the parsers that
produce the variables.**

```python
# engine/outputs.py: the parser and the rule are neighbours
COMBINE: dict[str, CombineRule] = {
    "kext":         ADDITIVE,
    "backscatter":  ADDITIVE,
    "vol_dens":     ADDITIVE,
    "cross_dens":   ADDITIVE,
    "g":            KSCA_WEIGHTED,
    "phase":        KSCA_WEIGHTED,
    "ssa":          SSA_RATIO,
    "reff":         from_moments("vol_dens", "cross_dens"),
    "angstrom_ext": respectral("kext"),
    "lidar_ratio":  reratio("kext", "backscatter"),
    ...
}
```

`mix.py` then contains **no variable name at all**: it does `COMBINE[name]` and
raises on a miss. Adding a parser without its rule fails on the first `Mix`
instead of silently producing wrong numbers on the tenth.

### 6.2 Three currencies for the same weight

Optical properties are **linear in `n`**. Scaling `n` by alpha scales `kext`,
`ksca`, `backscatter`, `mass_conc`, `vol_dens` and `cross_dens` by alpha, and
leaves `ssa`, `g`, `phase`, `reff`, the Angstrom exponents and the lidar ratio
unchanged.

So the weights condition **no** MOPSMAP run. They apply afterwards, which means
they can be **solved once the per-specie `kext` are known**. That is impossible
if all modes are flattened into a single launch file, where the concentrations
must be known before launching. This is a second, independent reason for the
compute-then-combine decision above.

| Currency | Constructor | Requires |
|---|---|---|
| Number concentration, m-3 | `pm.Mix({s: n})` | nothing |
| Mass concentration, kg m-3 | `pm.Mix.from_mass({s: c})` | `density_dry`, now mandatory in the schema |
| Optical depth fraction | `pm.Mix.from_optical_depth({s: f}, wl_ref=, rh_ref=)` | one run per specie first |

`from_mass` earns its place: CAMS provides mass mixing ratios, not number
densities, so it is probably the most useful currency for MAJA-type work.

```python
mix = pm.Mix.from_optical_depth(
    {pm.CAMS.SULPHATE: 0.30, pm.CAMS.DUST: 0.70},
    wl_ref=0.56,
    rh_ref=50,
)
op = mix.compute(rh=50, wl=wl)
op.kext          # normalised extinction, sums to 1 at wl_ref
op.kext * 0.25   # for a total AOD of 0.25
```

Normalising avoids requiring a layer thickness: `ssa`, `g` and the phase
function are scale invariant, so they are exact without it, and `kext` comes out
normalised, to be multiplied by the AOD of interest.

### 6.3 Why `rh_ref` is required and has no default

The invariant quantity is `n`, the number concentration. Humidity grows
particles, it neither creates nor destroys them.

| Quantity | Depends on `rh`? | Status |
|---|---|---|
| `n_i` (particles per m3) | no | the physical state of the air mass |
| `sigma_i(rh, wl)` (per-particle cross section) | yes | what MOPSMAP computes |
| `f_i(rh, wl)` (AOD fraction) | **yes** | **derived**: `n_i*sigma_i / sum_j(n_j*sigma_j)` |

Optical depth fractions are therefore a derived observation. Recovering the
state means inverting them:

```
n_i  proportional to  f_i / sigma_i(rh_obs, wl_obs)
```

So `rh_ref` is not a convention, it is **the humidity at which the fractions
were observed**. Given "30% of AOT at 560 nm, and 50% humidity, at 15h40", those
are one consistent pair of observations and `rh_ref = 50`.

Asking afterwards for `rh=80` then yields what **that same air mass** would look
like at 80% humidity. The fractions there are no longer 0.30 / 0.70, because
sulphate and dust do not grow alike. That is the correct behaviour.

**There is no default.** A default would silently assume the fractions hold at
some particular humidity, which is exactly the assumption not to make on the
user's behalf. When the fractions and the humidity come from the same
observation, `rh_ref` equals `rh` and the number is typed twice; that redundancy
is where a user notices the two can differ.

Forcing the fractions to stay 0.30 / 0.70 across an `rh` sweep would mean `n_i`
varying with humidity. That is not an alternative mode, it is unphysical.

`mix.weights` is inspectable so the inversion can be checked:

```python
mix.weights      # {CAMS.SULPHATE: 2.9e9, CAMS.DUST: 4.1e7}  in m-3
```

Consequence for the xsweep store: the key covers only
`(specie, rh, wl, outputs)`. Changing currency or fractions triggers **no**
recomputation.

### 6.4 Constraints

- Species must share the `wl` grid. This is already required by MOPSMAP itself:
  `wl_command` emits `wavelength from_refrac_file` for multi-wavelength runs, so
  the grid comes from the refractive-index file and is global to the launch file.
- The distribution amplitude (`n`, or `A`, or `concentrations` depending on the
  psd type) is **required per mode**, because the ratio between the modes of a
  species is part of what defines that species: a bimodal sea salt is not
  computable without knowing how much fine relative to coarse. A species file is
  therefore self-contained. What belongs to the mixture is the **absolute
  scale**: `Mix` multiplies a whole species by one factor, preserving the ratio
  between its modes. It never replaces the amplitude, since a single number
  could not express several modes.
- `OPAC_MIXES` stays a Python table (literature values, not gridded data), and
  `OpacMixName.CONTINENTAL_AVERAGE` becomes a constructor returning a `Mix`.

---

## 7. xsweep integration

| Axis | xsweep role | Reason |
|---|---|---|
| `wl` | `vec` | MOPSMAP takes the whole spectral grid in one launch file |
| `rh` | `loop` | PSD and refractive index change, one run per RH |
| swept microphysical params | `loop` | one run per point |
| `theta`, `l`, `element` | output dims | declared in the contract |

Contract shape: `"vec(wl), loop(rh) -> optics(wl)"`.

This replaces `ParametricSweep`, `extend_optiprops` **and** `cache/results.py`,
and adds resumability.

Pitfalls:

1. The xsweep memoization key covers swept points only. It captures neither the
   specie identity nor the source file version, so both must go into `version=`
   (e.g. `f"cams-49r1-sulphate-{schema_rev}"`). Otherwise two species collide in
   the same store.
2. The key must also cover the output request and `n_angles`, which the current
   cache handles for outputs and misses for `n_angles`.
3. `compute` must accept scalar / list / `ndarray` / `DataArray` and normalise to
   an `xr.Dataset` space. One `_as_space()` helper, used identically for `rh` and
   for any swept parameter.
4. The output contract must enumerate variables and their dims. The current
   merge namespace is flat and resolves collisions implicitly
   (`angstrom_ext` appears in both the integrated and the lidar blocks,
   "integrated takes precedence"). Make that explicit.

---

## 8. Public surface

Target: replace 30 exported names with 12.

| Keep | Why |
|---|---|
| `load`, `Specie`, `Mix`, `Mode` | the whole entry point |
| `CAMS`, `OPAC` | source enums |
| `Sphere`, `Spheroid`, `Irregular`, ... | irreducible, but moved under `pm.shapes` |
| `Lognormal`, `ModifiedGamma`, ... | irreducible, but moved under `pm.psd` |
| `OutputType` | output selection |

| Remove | Replacement |
|---|---|
| `pm.compute`, `pm.kext`, `pm.ssa`, `pm.phase` | `specie.compute(...)`, then attribute access on the result |
| `ParticleMixture`, `ParametricSweep` | `Specie`, `Mix` |
| `cams_to_kext`, `cams_to_optiprops`, `cams_to_smartg` | `load(CAMS.X).compute(...)` |
| `OpacMix`, `OpacHumidityMode` | `Mix`, `growth` attribute |
| `read_aerosol_microphysical_parameters` | `load` |
| `cache_status`, `prefetch` (module level) | `specie.prefetch(...)`, `specie.cache_status(...)` |
| `OptiProps` | see below |

### 8.1 Drop `OptiProps`, return `xr.Dataset`

`OptiProps` is a frozen dataclass whose only methods are `save`, `sel(prop)` and
`coord`. `sel` shadows `xr.Dataset.sel` with incompatible semantics, which is
actively confusing, and the wrapper blocks the whole xarray API.

Return a plain `xr.Dataset` and register an accessor for the pymopsmap-specific
operations:

```python
op = specie.compute(...)
op.kext                      # plain xarray
op.sel(wl=0.55, method="nearest")
op.mopsmap.to_smartg(path)   # accessor, replaces adapters/output/smartg.py
```

This alone removes `OptiProps`, `extend_optiprops`, `pm.kext`, `pm.ssa` and
`pm.phase` from the public surface.

The SMART-G accessor must take the humidity dim name as an argument
(default `"rh"`) instead of hardcoding `op.coord("rh")`, since v2 lets the user
name that dim.

---

## 9. Package layout

```
src/pymopsmap/
+-- __init__.py            # ~12 exports (section 8)
+-- exceptions.py
|
+-- species/               # LEVEL 1: the description
|   +-- __init__.py        #   load()
|   +-- specie.py          #   Specie: DataTree wrapper, growth validation,
|   |                      #   point materialisation
|   +-- mix.py             #   Mix + the combination rules of section 6
|   +-- schema.py          #   the NetCDF contract, read AND write:
|   |                      #   units, growth, schema_rev, assertions
|   +-- catalog.py         #   CAMS / OPAC -> importlib.resources paths
|
+-- shapes.py              # LEVEL 2 public: Sphere, Spheroid, Irregular, ...
+-- psd.py                 # LEVEL 2 public: Lognormal, ModifiedGamma, ...
+-- microparams.py         # LEVEL 2 internal: MicroParameters, one valid point
|
+-- engine/                # one point -> one MOPSMAP run
|   +-- __init__.py        #   run_point(modes, outputs) -> xr.Dataset
|   +-- workspace.py       #   isolated temp dir per run   <- fixes 2.6
|   +-- commands.py
|   +-- launch_file.py
|   +-- launcher.py
|   +-- outputs.py         #   OutputType + parsers (was output_format.py)
|   +-- coverage.py
|
+-- sweep.py               # xsweep binding: _as_space(), contract, version key
|
+-- scatlib/               # the MOPSMAP optical dataset (was cache/)
|   +-- __init__.py
|   +-- cache.py           #   was optical.py
|   +-- downloader.py
|   +-- resolver.py
|
+-- accessors.py           # .mopsmap accessor on xr.Dataset (to_smartg, ...)
|
+-- data/                  # PACKAGE DATA (section 10)
|   +-- cams/{47r1,48r1,49r1}.nc
|   +-- opac/opac.nc
|
+-- utils/
    +-- logging.py
    +-- types.py
```

Outside the package:

```
scripts/build_catalog/
+-- cams.py                # MAJA table -> canonical schema
+-- opac.py                # was opac_download.py, GEISA parsing
```

### 9.1 Rationale for the non-obvious choices

| Choice | Reason |
|---|---|
| `scatlib/` rather than `cache/` or `dataset/` | It is MOPSMAP's own word (`scatlib '<path>'` in the launch file). `cache/` covered two unrelated caches, the result cache dies, only this one remains. `dataset/` would collide with `xr.Dataset`. |
| `shapes.py` and `psd.py` **flat** at package root | Section 8 exposes `pm.shapes.Sphere` and `pm.psd.Lognormal`. The import path should mirror the public namespace, not add a `micro/` level the user never sees. |
| `engine/workspace.py` is new | Replaces `utils/temp.py`. One directory per run, `ri_1.txt`, `ri_2.txt`, ... per mode. Fixes section 2.6 and unblocks xsweep parallelism. |
| `species/schema.py` does read **and** write | `Specie.to_netcdf()` is public (round-trip test, section 5), so writing the canonical format belongs in the package. Only source-specific **ingestion** (GEISA regex, CAMS reshaping) goes to `scripts/`. |
| `OutputType` lives in `engine/outputs.py` | The enum value **is** the file extension MOPSMAP writes (`otype.value` in the parser). Separating them guarantees they drift. Re-exported from `__init__.py`. |
| `sweep.py` is a module, not a package | Space normalisation, contract, version key. Promote it if it passes 200 lines. |

### 9.2 What disappears

| Deleted | Replaced by |
|---|---|
| `adapters/` **entirely** | `species/catalog.py` (input), `accessors.py` (output), `scripts/build_catalog/` (ingestion) |
| `models/particle_mixture.py` | `species/specie.py` + `species/mix.py` |
| `models/optiprops.py` | bare `xr.Dataset` + `accessors.py` (section 8.1) |
| `cache/results.py` | the xsweep zarr store |
| `utils/caching.py` | the xsweep `version=` key, built in `sweep.py` |
| `utils/temp.py` | `engine/workspace.py` |
| `main.py` | nothing, a demo script calling `cams_to_kext` |
| `commands.py:output_request_commands` | already dead, duplicate of `_file_suffix` |
| `output_format.py:format_netcdf_file` | already dead |

`models/` disappears as a package: `microparams.py`, `shapes.py` and `psd.py`
move to the root, `optiprops.py` and `particle_mixture.py` die,
`output_request.py` merges into `engine/outputs.py`.

---

## 10. Data and packaging (settled)

Three data mechanisms coexist today, and only one is right:

| Data | Current mechanism | Verdict |
|---|---|---|
| MOPSMAP optical dataset (GB) | downloaded from `PYMOPSMAP_DATASET_SOURCE` into `~/.cache/pymopsmap/dataset/` | **Correct**, keep as is |
| CAMS (797 kB) | committed to `data/cams/`, read via `DATA_PATH` | absent from the wheel: `pip install` ships a library with no data |
| OPAC | downloaded at runtime from GEISA **into `DATA_PATH`** | worse: writes into the install directory, which may be read-only, and gates the first call on a third-party server |

**Decision: three tiers.**

| Tier | Content | Mechanism |
|---|---|---|
| 1 | built-in catalogue (CAMS, OPAC) | `src/pymopsmap/data/`, read via `importlib.resources` |
| 2 | user species | `pm.load("my_aerosol.nc")`, section 5 |
| 3 | MOPSMAP optical dataset | unchanged, cache plus download |

Why package data rather than download, in order of weight:

1. **`schema_rev` couples the file to the reader.** A downloaded catalogue would
   need version negotiation, both directions. Shipped together, they cannot
   disagree. This argument did not exist before the schema was versioned.
2. **Size is a non-issue.** 797 kB for CAMS, all versions; OPAC is 10 species
   times 8 RH values, tens of kB. Far under the 100 MB PyPI per-file limit.
3. **No network on the nominal path.** `pm.load(CAMS.SULPHATE)` works offline,
   in CI, behind a proxy.
4. **The escape hatch already exists.** `load()` takes a path, so extensibility
   needs no second mechanism. A remote catalogue stays addable later without an
   API change.

Consequences:

- `opac_download.py` becomes `scripts/build_catalog/opac.py`, a maintainer
  script. Its GEISA parsing is regex over text and has no business running on a
  user machine; the result is deterministic, so it runs once.
- `DATA_PATH` is deleted. Tests pass paths to `load()` instead of overriding a
  global.
- Explicit hatchling config, then verify:

```toml
[tool.hatch.build.targets.wheel]
packages = ["src/pymopsmap"]
```

```bash
python -m build && unzip -l dist/*.whl | grep '\.nc'
```

### 10.1 Repository hygiene

| Finding | Action |
|---|---|
| `data/smartg/` weighs **207 MB** and is committed. These are **generated outputs** (`{specie}_sol.nc` from `create_lut_for_smartg`), not inputs. They dominate the 130 MiB git pack. | remove from the repository and from history |
| `cams_aerosols_optical_properties.nc` (693 kB) is referenced by no code. No grep of the filename hits `src/` or `tests/`. | either it is the validation reference and deserves a test that uses it, or it is dead |
| `.gitignore` contains `data/**`, but the files are already tracked, so the rule is inert for them | clean up once `data/` moves into the package |

---

## 11. Migration order

Section 2.6 changes the priority: `engine/workspace.py` comes first. It is a
live silent-wrong-result bug on every multi-mode computation, it is independent
of the rest of the redesign, and the fix is small.

1. `engine/workspace.py`, plus a test asserting that a two-mode launch file
   names two distinct `refrac file` paths.
2. Fixes 2.4 (interpolation domain) and 2.5 (clipping reindex).
3. `species/schema.py` + `scripts/build_catalog/cams.py` + the round-trip test.
4. `species/specie.py` + `load`, keeping the current engine.
5. `sweep.py`, delete `cache/results.py` and `extend_optiprops`.
6. `species/mix.py`, rewrite OPAC as a `Mix` constructor.
7. `accessors.py`, bare `xr.Dataset`, shrink `__init__.py`.

Section 2.7 items are folded into whichever step touches their file.

---

## 12. Open questions

1. A custom `Specie` with `growth="none"` and a user passing `rh=`: raise, or
   silently ignore? (Recommendation: raise.)
2. ~~Should `growth="mopsmap"` be supported?~~ **Settled**: three values,
   `tabulated` / `kappa` / `none`, growth is always applied by MOPSMAP for the
   kappa mode, the Python formula is deleted. See section 2.3.
3. ~~Where do the species files live?~~ **Settled**: package data via
   `importlib.resources`, three tiers. See section 10.
4. ~~Does `Mix` accept species with different `growth` modes in one call?~~
   **Settled**: yes, because `Mix` computes each specie in its own run. MOPSMAP
   could not do it in a single launch file anyway (section 2.3.1).
