# Changelog

## 0.5.0

A validation campaign against the reference article, and the twelve defects it
found. Every table of Gasteiger and Wiegner (2018) is now asserted, and nine
of its figures are redrawn from this package.

### Breaking

- The OPAC catalogue follows the presets MOPSMAP ships rather than OPAC as
  published: the mineral components are prolate spheroids with the aspect
  ratio distribution of Kandler et al. (2009), per section 4.2 of the user
  guide. Its extinction-to-mass factor moves by 7 % and its mass-to-backscatter
  factor by a factor of two. The coarse sea salt mode is cut at 20 um rather
  than 60: at 60 it leaves the data set at 355 nm, where Figure 5 plots it.
- `run_point` no longer takes `size_equ`. A mode carries it, and a run whose
  modes disagree is refused.
- Logging goes through the standard library. `structlog` and the
  `PYMOPSMAP_LOG_LEVEL` and `PYMOPSMAP_LOG_FILE` variables are gone; a library
  with no handler emits nothing below WARNING, which is what it should do.
- `structlog` and `einops` are no longer dependencies, and `miepython` is a
  test dependency.

### Added

- `docs/guide.md`, illustrated, and `docs/validation.md`, which says what
  agrees with the article and what does not.
- `Mode.from_particles` turns a table of measured particles into modes,
  binning aspect ratios onto the grid the data set holds the way
  `init_shape.f90` would. `Mode.from_index_distribution` spreads a mode over a
  measured distribution of imaginary indices.
- `size_equ` and `nonabs_fraction` on a mode, both of which section 2.1 and
  section 3.1 of the article define and neither of which was reachable.
- `n_angles` on `compute`, which sets the shape of every angular output and so
  keys its store.
- Nine scripts under `scripts/validation/` reproducing figures 2, 4, 5, 6, 7,
  8, 9, 10 and 11, and `scripts/demo/gallery.py` drawing the guide.

### Fixed

- `compute` could not return any angular output: the contract declared every
  one as a single number per wavelength, so a phase function, a scattering
  matrix, a volume scattering function or an expansion failed on its first
  point.
- A swept axis named after its parameter came back as 0, 1, 2 instead of the
  values given, because xarray promotes such a variable to a coordinate.
- A species file stayed open in the xarray file cache; once that cache evicted
  it, every later read of the tree raised.
- A wavelength clipped out of coverage reset the non-absorbing fraction and the
  size equivalence of its mode.
- A spheroid of aspect ratio one is read from the sphere files, as MOPSMAP
  reads it.
- An aspect ratio distribution file is quoted in the launch file: a
  list-directed Fortran read ends the record on an unquoted slash.
- The eps grid of a spheroid distribution file is read from the file rather
  than taken whole, and the index entry at eps = 1 no longer rejects every
  size.
- A run whose every wavelength leaves the data set is refused with its reason
  rather than handed to MOPSMAP empty.
- The sweep store is keyed on the species description, so rebuilding a
  catalogue file no longer serves its old results.

### Performance

- The data set index is read once rather than once per swept point, which also
  removes a race: netCDF4 is not thread safe.
- The size limit and the required files are asked once per distinct question.
  Clipping nine thousand modes over forty-nine wavelengths drops from 562 s to
  under one.
- A refractive index file is named after its contents, so the thousands of
  modes of a measured ensemble share one.

## 0.4.0

A rewrite of the public API around a declarative description of an aerosol,
and five correctness fixes found on the way.

### Breaking

- `pm.compute`, `pm.kext`, `pm.ssa`, `pm.phase`, `ParticleMixture` and
  `ParametricSweep` are replaced by `pm.load(...).compute(...)`.
- `OptiProps` is gone: every entry point returns a plain `xarray.Dataset`, and
  the conversions this library adds live on a `.mopsmap` accessor.
- The `adapters` package is gone. `cams_to_kext`, `cams_to_optiprops`,
  `cams_to_smartg` and `OpacMix` are replaced by the catalogue and `Mix`.
- Shapes and size distributions move under `pm.shapes` and `pm.psd`. The public
  surface drops from thirty names to fifteen, pinned by a test.
- `DATA_PATH` is gone. The species catalogue ships inside the wheel.

### Added

- A canonical NetCDF schema for species, and the CAMS and OPAC catalogues
  built into it. A hand-built species and a catalogue one serialise
  identically.
- `Specie.custom` and `Mode` to describe an aerosol in Python, with swept
  parameters carried by the species rather than by the compute call.
- `Mix`, weighted by number concentration, by mass, or by fraction of the
  total optical depth.
- The scattering coefficient `ksca`, derived once at parse time.
- `DomainError` and `CoverageError`, which say what a computation cannot do
  and why.
- A swept parameter goes where its scalar would: every numeric field of a
  size distribution or a shape accepts a `DataArray`, and `compute` walks the
  dimensions the species declares. `Mode(sweep=...)` is gone, with the
  placeholder value it forced alongside it.
- Swept parameters and mixture weights accept a `DataArray`, so a scene whose
  humidity and composition vary per pixel yields a result of that shape.
  Repeated values cost one run. A mixture reports its resolved concentrations
  as a `concentration(specie, ...)` variable rather than an attribute, which an
  array could not be.
- Sweeps run through [xsweep](https://github.com/walcark/xsweep): results are
  memoised in a store, an interrupted sweep resumes, and the hand-written
  result cache is gone. Each run owns its directory, so points can run in
  parallel.

### Fixed

- Every mode of a mixture shared one refractive index file, so MOPSMAP read
  the last mode's index for all of them. Affected every CAMS species and every
  OPAC mix.
- Refractive indices were written in fixed point with six decimals, so any
  imaginary part below 5e-7 reached MOPSMAP as zero and nearby wavelengths
  collapsed onto the same value.
- Out-of-domain interpolation returned NaN silently, which then resolved to an
  arbitrary optical dataset file.
- MOPSMAP reports a missing dataset file on stdout and still exits zero, so a
  partially failed run returned a mixture of correct values and silent gaps.
- Dataset files and size-parameter coverage were resolved on the dry
  refractive index, while MOPSMAP grows the particles itself.
- Clipping a wavelength rebuilt the result as one-dimensional, which broke
  every output carrying an angle.
- The size-parameter limit was one hardcoded value per shape family, while the
  merged spheroid limit spans three orders of magnitude across the refractive
  index grid.
