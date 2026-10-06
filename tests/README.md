# Tests

Three directories, told apart by what they need and what they prove.

| | Needs | Proves |
|---|---|---|
| `unit/` | nothing but the package | that the code does what it says: validation rules, the launch file it writes, the parsers, the combination rules, the sweep plumbing. MOPSMAP is stubbed |
| `integration/` | the MOPSMAP binary and the optical data set | that the whole pipeline runs and comes back with the axes it promised |
| `validation/` | the same, plus the **extended** archive | that the numbers are right, against a published reference the package had no hand in |

```bash
pixi run test                                 # unit only, no data set needed
export PYMOPSMAP_DATASET_SOURCE=/path/to/optical_dataset
pixi run -e dev pytest tests/integration tests/validation
```

`integration/` and `validation/` skip themselves when
`PYMOPSMAP_DATASET_SOURCE` is unset, so `pixi run test` is green on a machine
that has never downloaded the 40 GB.

## What validation is held to

| File | Reference |
|---|---|
| `test_gasteiger_2018.py` | Tables 3, 4, 5 and 6 of Gasteiger and Wiegner (2018), GMD 11:2739, and the reference runs the authors ship in `bin/mopsmap/misc/paper_examples/` |
| `test_mie_reference.py` | `miepython`, an implementation of Mie theory that shares no code with MOPSMAP. The only check that does not go through the MOPSMAP data set |

`docs/validation.md` reports what agrees and what does not, figure by figure
and table by table.

## Data

`tests/data/` holds fixtures: a small dataset cache and the legacy CAMS values
the catalogue is pinned against. It is **not** the library's data. What the
package ships lives in `src/pymopsmap/data/`, and is built by
`scripts/build_catalog/`.
