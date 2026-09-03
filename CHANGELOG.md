# Changelog

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
