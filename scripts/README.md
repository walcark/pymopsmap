# Scripts

None of these is imported by the package. They produce things: the data it
ships, the figures the documentation shows, and the figures that justify
trusting it.

| | What it does | Writes to |
|---|---|---|
| `build_catalog/` | builds the species catalogue the wheel ships, from the MOPSMAP data directory and the published tables | `src/pymopsmap/data/` |
| `validation/` | redraws the figures of Gasteiger and Wiegner (2018) and prints its own numbers beside the published ones | `docs/figures/gasteiger_*.png` |
| `demo/` | draws the panels of the guide, and shows a sweep over a real CAMS scene | `docs/figures/guide-*.png` |

## Rebuilding the catalogue

The outputs are committed, so this is only needed when a source changes.

```bash
pixi run -e dev python -m scripts.build_catalog.cams
pixi run -e dev python -m scripts.build_catalog.opac
```

## Redrawing the figures

Both need the optical data set, both archives.

```bash
export PYMOPSMAP_DATASET_SOURCE=/path/to/optical_dataset
pixi run -e dev python -m scripts.validation.gasteiger_fig5
pixi run -e dev python -m scripts.demo.gallery
```

`scripts/validation/gasteiger_fig11.py` also needs the supporting information
of Vogel et al. (2017), which has to be fetched from the publisher; the script
says so and where to put it.
