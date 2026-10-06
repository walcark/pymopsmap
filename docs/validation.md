# Reproducing Gasteiger and Wiegner (2018)

> Gasteiger, J. and Wiegner, M., *MOPSMAP v1.0: a versatile tool for the
> modeling of aerosol optical properties*, Geosci. Model Dev. 11, 2739-2762,
> 2018. [doi:10.5194/gmd-11-2739-2018](https://doi.org/10.5194/gmd-11-2739-2018)

Nine figures and four tables of the article, recomputed through this wrapper.

## What is reproduced

| | Target | Agreement |
|---|---|---|
| Fig. 2 | single particles against size parameter | the two values the text quotes, to 1e-3 |
| Fig. 4 | sampling and interpolation error | spheres, against `miepython` |
| Fig. 5 | the ten OPAC types against humidity | curve by curve, read at 600 dpi |
| Fig. 6 | phase functions of five dust bins | the 40 values of Table 4, to 6e-4 |
| Fig. 7 | the desert type against the cutoff radius | the 6 percentages of section 5.3 |
| Fig. 8 | one size distribution, three equivalences | the conversion the text gives |
| Fig. 9 | dust against its index variability | the 11 values of section 5.6, to 1e-3 |
| Fig. 10 | nephelometer truncation | the 3 statements of section 5.7 |
| Fig. 11 | nine volcanic ashes | the 5 statements of section 5.8 |
| Table 3 | one lognormal mode, two indices | 20 values, to 1.5e-3 |
| Table 4 | the five COSMO-MUSCAT dust bins | 20 values, to 1.5e-3 |
| Table 5 | the three size equivalences | 40 values, six digits |
| Table 6 | the Jacobian of a dust ensemble | 3 values to 1.5e-3, 9 derivatives to 0.1 |

Figures 1 and 3 are diagrams. Tables 5 and 6 are checked against the reference
runs the authors ship in `bin/mopsmap/misc/paper_examples/`, which carry six
digits where the article prints four.

`tests/validation/test_mie_reference.py` adds the one check that does not go
through the MOPSMAP data set: the same ensembles computed with `miepython`,
agreeing to 3e-3.

## Running it

Both archives of the optical data set are needed; soot reaches a real
refractive index of 1.75, past the 1.64 the main one covers.

```bash
export PYMOPSMAP_DATASET_SOURCE=/path/to/optical_dataset
pixi run -e dev test-validation                               # the tables
pixi run -e dev python -m scripts.validation.gasteiger_fig5   # a figure
```

Each figure script prints the published values beside its own.

## The figures

### Figure 2

![Figure 2](figures/gasteiger_fig2.png)

m = 1.56 + 0.00215i, a grid point of the data set. The albedo peaks at 0.9917
against "about 0.991", and at xc = 1000 the spheres and spheroids give 0.5516,
0.5518 and 0.5549 against "about 0.551". The irregular curves stop at
xc = 30.2, where Table 2 says their coverage ends.

### Figure 4

![Figure 4](figures/gasteiger_fig4.png)

Panels (b) and (d) are complete: the four grid points around m = 1.54 + 0.005i,
the value MOPSMAP interpolates there, and `miepython` on the same size grid.
Panels (a) and (c) show the spheres only; the other pair needs a T-matrix
calculation finer than the data set holds.

### Figure 5

![Figure 5](figures/gasteiger_fig5.png)

The values are not tabulated, so the comparison is panel by panel. The
single-scattering albedo at 532 nm and RH = 0:

| type | published | computed |
|---|---|---|
| urban | 0.673 | 0.6731 |
| continental polluted | 0.785 | 0.7851 |
| continental average | 0.845 | 0.8453 |
| desert | 0.867 | 0.8672 |
| arctic | 0.749 | 0.7499 |
| maritime polluted | 0.927 | 0.9260 |
| continental clean | 0.939 | 0.9388 |

Two properties of the OPAC ensembles come from section 4.2 of the MOPSMAP user
guide rather than from OPAC: the mineral components are prolate spheroids with
the aspect ratio distribution of Kandler et al. (2009), and the coarse sea salt
mode is cut at 20 um. The catalogue this package ships follows it.

### Figure 6

![Figure 6](figures/gasteiger_fig6.png)

Section 5.2 fixes every input. The forty values of Table 4 come back within
6e-4.

### Figure 7

![Figure 7](figures/gasteiger_fig7.png)

The authors ship the script that makes this figure, so the inputs are known
rather than inferred.

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

### Figure 8

![Figure 8](figures/gasteiger_fig8.png)

No MOPSMAP run: the figure shows what `size_equ` means before any optics. The
shapes, the crossings and the three cutoffs match, and the conversion the text
gives comes back exactly, 0.757 um. The vertical scale differs by a factor
2.504, 2.513 and 2.515 on the three curves, against sqrt(2 pi) = 2.5066.

### Figure 9

![Figure 9](figures/gasteiger_fig9.png)

The red curve carries the index distribution of Kandler et al. (2011),
Appendix S1, fetched from the publisher; `KANDLER_DATA` points the script at
it. Its average comes out at 0.01746 against the 0.0175 section 5.6 gives.

| | published | computed |
|---|---|---|
| omega_0, average m_i | 0.741 | 0.7401 |
| omega_0, X = 0.5 | 0.834 | 0.8341 |
| g, measured / average / X = 0.5 | 0.744 / 0.789 / 0.749 | 0.7437 / 0.7895 / 0.7488 |
| S, measured / average / X = 0.5 | 41 / 78 / 42 sr | 40.70 / 78.35 / 41.96 sr |
| delta_l, measured / average / X = 0.5 | 0.241 / 0.212 / 0.220 | 0.2411 / 0.2114 / 0.2201 |

### Figure 10

![Figure 10](figures/gasteiger_fig10.png)

| statement of section 5.7 | computed |
|---|---|
| total scattering underestimated "by a factor of ~ 2" at rmod = 1 um | 2.16, 2.01, 1.83 |
| shape effect on Cts "less than 3 %" | 2.2 % |
| "maximum shape effect on Cbs is 7 %" | 7.0 % |

### Figure 11

![Figure 11](figures/gasteiger_fig11.png)

Each measured particle of Vogel et al. (2017) becomes a prolate spheroid of
its own size and aspect ratio, nine to thirty thousand modes per volcano. The
measurements are that article's supporting information, fetched from the
publisher; `VOGEL_DATA` points the script at it.

| statement of section 5.8 | computed |
|---|---|
| omega_0 spans "up to about 0.12" | 0.122 |
| it rises "typically by about 0.05" | 0.048 |
| variability in g "less than 0.05" | 0.014 |
| its change with wavelength "less than 0.02" | 0.018 |
| effective radii "from 9.5 to 21 um" | 9.5 to 21.0 |

Chaiten comes out the least absorbing, as the article says; the most absorbing
comes out as Grimsvotn rather than Mount Kelud, whose curves overlap at 0.807
against 0.803 at 300 nm.

## Inputs taken from elsewhere

**Table 3 gives no concentration.** Its extinction coefficients are reproduced
with N = 1000 cm-3; the other four rows do not depend on it.

**Section 5.4 gives two.** "N0 = 10^3.66 cm-3" and "a concentration of
N = 100 cm-3 in the range from rmin to rmax" differ by a factor of 44, since
the lognormal of `log_distr.f90` puts 96.46 % of its particles inside that
range. The in-range concentration reproduces the table, and MOPSMAP echoes it
back as n = 100.0.

**The reference runs in `sect_54_size_equivalence/` carry a factor 1e-6**, a
concentration in cm-3 where MOPSMAP reads m-3. The test applies the factor
rather than alter the shipped expectations.
