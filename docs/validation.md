# Validating pymopsmap against Gasteiger and Wiegner (2018)

What the article contains, what this wrapper reproduces, and what it does not.
The point is not the pictures. A published number that comes back out of the
wrapper is evidence that nothing was shifted between the launch file and the
result; a figure that matches is the same evidence, read with the eye.

Reference: Gasteiger, J. and Wiegner, M., *MOPSMAP v1.0: a versatile tool for
modeling aerosol optical properties*, Geosci. Model Dev. 11, 2739-2762, 2018.
[doi:10.5194/gmd-11-2739-2018](https://doi.org/10.5194/gmd-11-2739-2018)

## What is checked, and where

| | Target | State | Agreement | Where |
|---|---|---|---|---|
| 2 | Single particles vs size parameter | reproduced | the two values the text quotes, to 1e-3 | `scripts/validation/gasteiger_fig2.py` |
| 4 | Sampling and interpolation error | replaced | an independent Mie code, 3e-3 | `tests/integration/test_mie_reference.py` |
| 5 | OPAC types vs humidity | partial | 5 types of 10, statements of section 5.1 only | `scripts/validation/gasteiger_fig5.py` |
| 6 | Phase functions of dust size bins | reproduced | the 40 values of Table 4, to 6e-4 | `scripts/validation/gasteiger_fig6.py` |
| 7 | Desert aerosol vs cutoff radius | reproduced | the six percentages of section 5.3 | `scripts/validation/gasteiger_fig7.py` |
| T3 | One lognormal mode, two indices | reproduced | 20 values, 1.5e-3 | `tests/integration/test_gasteiger_2018.py` |
| T4 | Dust size bins at 500 nm | reproduced | 20 values, 1.5e-3; 40 in the figure script, 6e-4 | same |
| T5 | Size equivalence conventions | reproduced | 40 values, 1.5e-3 | same |
| T6 | Jacobian of a dust ensemble | reproduced | 3 values to 1.5e-3, 9 derivatives to 0.1 | same |

Figures 1, 3 and 8 are diagrams with nothing to compute. Figure 9 and section
5.6 onwards describe scripts shipped with MOPSMAP rather than the wrapper.

Running the checks needs the optical data set:

```bash
export PYMOPSMAP_DATASET_SOURCE=/path/to/mopsmap/optical_dataset
pixi run -e dev pytest tests/integration
```

## The four figures

### Figure 2, single particles against size parameter

![Figure 2 recomputed](figures/gasteiger_fig2.png)

Extinction efficiency, single scattering albedo and asymmetry parameter of one
particle at m = 1.56 + 0.00215i, for spheres, prolate spheroids of aspect
ratio 1.4 and 3.0, and irregular shapes D and F. The refractive index is a
grid point of the data set, so nothing is interpolated and the curves are the
data set read through the wrapper.

The two numbers the text quotes come back: the albedo peaks at 0.9917 against
"about 0.991", and at xc = 1000 the spheres and spheroids give 0.5516, 0.5518
and 0.5549 against "about 0.551". The irregular curves stop at xc = 30.2,
which is where Table 2 says their coverage ends.

Extinction efficiency is not a MOPSMAP output. It is `kext / cross_dens`, both
of which the integrated block gives.

### Figure 5, the OPAC types against relative humidity

![Figure 5 recomputed](figures/gasteiger_fig5.png)

Five of the ten types, for the reason given below.

### Figure 6, phase functions of five dust size bins

![Figure 6 recomputed](figures/gasteiger_fig6.png)

The five COSMO-MUSCAT bins at 500 nm, spheres against prolate spheroids with
the aspect ratio distribution of Kandler et al. (2009). Section 5.2 fixes
every input: bin edges at 0.1, 0.3, 0.9, 2.6, 8 and 24 um, constant dv/dlnr
within each bin, m = 1.53 + 0.0078i, and volume-equivalent sizes for the
spheroids. A two-point dv/dlnr table with equal values is exactly "constant
dv/dlnr", which is what `DistrListPSD` writes.

The forty numbers of Table 4 come back within 6e-4.

### Figure 7, the desert type against the cutoff radius

![Figure 7 recomputed](figures/gasteiger_fig7.png)

The authors ship the script that makes this figure, in
`bin/mopsmap/misc/paper_examples/sect_53_rmax_cutoff/`, so the four modes,
their concentrations and their shapes are known rather than inferred: water
soluble as spheres, the three mineral modes as spheroids with the Kandler
distribution.

Every percentage section 5.3 quotes comes back:

| | published | computed |
|---|---|---|
| PM10 mass | 59.5 % | 59.5 % |
| PM2.5 mass | 21.6 % | 21.6 % |
| PM10 cross section | 94.4 % | 94.4 % |
| PM2.5 cross section | 69.0 % | 69.0 % |
| cross section at rmax = 10 um | 97.8 % | 97.8 % |
| mass at rmax = 10 um | 75.6 % | 75.5 % |
| PM2.5 raises omega_0 by | 0.035 to 0.071 | 0.0352 to 0.0708 |
| PM2.5 lowers g by | 0.02 to 0.04 | 0.021 to 0.0434 |

## Why the tables carry the weight

A figure has to be read off an axis, which costs a digit or two and invites
agreement that is not there. Tables 3 to 6 print four significant digits, and
the authors also ship their own reference runs in
`bin/mopsmap/misc/paper_examples/`, whose output files carry six. Those are
what the test suite is written against.

For Table 5 the agreement is six digits on twelve quantities at once,
including the extinction coefficient, the cross section density, the mass
concentration, the lidar ratio and the linear depolarisation ratio.

## Two statements of the article that do not hold

**Table 3 gives no concentration.** Its extinction coefficients are
reproduced with N = 1000 cm-3, which is the only value that fits; the other
four rows of the table do not depend on it.

**Section 5.4 says "N0 = 10^3.66 cm-3, which results in a concentration of
N = 100 cm-3 in the range from rmin to rmax".** The two halves disagree: with
rmod = 0.5 um, sigma = 2 and rmax = 1.75 um, the lognormal of `log_distr.f90`
puts 96.46 % of its particles inside that range, so N0 = 10^3.66 would give
4410 cm-3, not 100. The in-range concentration is the half that reproduces
the table, and MOPSMAP echoes it back as n = 100.0.

Separately, the reference runs in `sect_54_size_equivalence/` were made with
the concentration in cm-3 where MOPSMAP reads m-3, so every extensive
quantity in them sits 1e6 below the published table. The test applies the
factor rather than alter the shipped expectations.

## Figure 5, the one that stays partial

Two gaps, both structural.

**Soot.** Five of the ten OPAC types contain it, and its real refractive index
reaches 1.75, outside the 1.28 to 1.64 of the main archive of the optical data
set. They need the extended archive.

**Sea salt at high humidity.** The OPAC coarse sea salt mode runs to
rmax = 60 um. At RH = 70 % the growth factor of the kappa parameterisation
takes that to 88 um, a size parameter of 1038 at 532 nm, past the 1013 the
data set covers at that refractive index. The clipping is not the wrapper
being cautious: MOPSMAP refuses the same run, here at RH = 90 %.

```
Error: Maximum specified particle size 1.259E+02 at wavelength 5.320E-01
is not covered by .../spheres/sphere_1.3600_0.000538.nc
```

(`add_contribution.f90`, lines 122-124). The article plots those points, so
its own runs used a cutoff below the OPAC value. It does not say which, so the
script leaves them out rather than guess.

What remains of Figure 5 is checked against the statements section 5.1 makes
in words, since the figure tabulates nothing.

## Six defects this campaign found

Each was a real bug, fixed in its own commit, with a test.

| Defect | Effect |
|---|---|
| An unquoted path in the launch file | a Fortran list-directed read ends on a slash, so an aspect ratio file never reached MOPSMAP |
| Escaped braces in the overlay command | `shape irregular_overlay` wrote `{self.xmin}` literally |
| The eps grid read whole | the resolver asked for all 31 aspect ratios, one of which has no file |
| `max_sizepara` at eps = 1 | it is zero, there being no spheroid file for a sphere, and it rejected every size |
| A fully clipped run | MOPSMAP got an empty wavelength file and reported a missing wavelength |
| The clipping warning | it sized the dry particle where the mask sized the grown one |

Two gaps in the wrapper were also closed: `size_equ`, without which two
columns of Table 5 and half of Table 4 cannot be expressed, and a logging
setup that forced DEBUG onto stdout and opened a file in the working
directory on import.
