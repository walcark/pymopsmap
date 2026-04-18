"""MOPSMAP computation pipeline — file resolution, execution, result caching."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pymopsmap.models import (
        MicroParameters,
        MicroParametersDispatch,
        OptiProps,
    )
    from pymopsmap.models.output_request import OutputRequest

    Modes = MicroParameters | list[MicroParameters]


def compute_optical_properties(
    dispatch: Modes | MicroParametersDispatch,
    output_types: OutputRequest | None = None,
    rh: float | None = None,
    quiet: bool = False,
) -> OptiProps:
    """
    Compute optical properties for a single mode set or a dispatch.

    Transparent pipeline:
      1. Check result cache → return hit immediately.
      2. Ensure required dataset files are present (download if missing).
      3. Write MOPSMAP launch file → run MOPSMAP → parse outputs.
      4. Store result in cache → return.
    """
    from pymopsmap.models import MicroParametersDispatch, extend_optiprops
    from pymopsmap.models.output_request import DEFAULT_OUTPUT

    if output_types is None:
        output_types = DEFAULT_OUTPUT

    if not isinstance(dispatch, MicroParametersDispatch):
        return _compute_single(dispatch, output_types=output_types, rh=rh, quiet=quiet)

    return extend_optiprops(
        index=dispatch.params,
        optiprops_li=[
            _compute_single(modes, output_types=output_types, rh=rh, quiet=quiet)
            for modes in dispatch
        ],
    )


def _compute_single(
    modes: Modes,
    output_types: OutputRequest,
    rh: float | None,
    quiet: bool,
) -> OptiProps:
    from pymopsmap.cache.downloader import DatasetDownloader
    from pymopsmap.cache.optical import OpticalDatasetCache
    from pymopsmap.cache.resolver import NCFileResolver
    from pymopsmap.cache.results import ResultCache
    from pymopsmap.utils import DATASET_CACHE_DIR

    from .launch_file import write_launching_file
    from .launcher import launch_mopsmap
    from .output_format import format_mopsmap_outputs

    result_cache = ResultCache()
    dataset_cache = OpticalDatasetCache()
    downloader = DatasetDownloader(cache=dataset_cache, quiet=quiet)

    key = result_cache.key(modes, output_types)
    cached = result_cache.get(key)
    if cached is not None:
        return cached

    # Ensure index.nc is present first
    index_path = dataset_cache.full_path("index.nc")
    if not dataset_cache.is_cached("index.nc"):
        downloader.download("index.nc")

    # Resolve required dataset files and download missing ones
    resolver = NCFileResolver(index_path)
    mp_list = [modes] if not isinstance(modes, list) else modes
    required = resolver.resolve(mp_list)
    downloader.download_missing(required)

    # Run MOPSMAP
    paths = write_launching_file(
        mp=modes,
        output_types=output_types,
        rh=rh,
        mopsmap_data_path=DATASET_CACHE_DIR,
    )
    out_mopsmap = launch_mopsmap(input_filename=paths["mopsmap"])
    out_mopsmap["ascii_base"] = paths.get("ascii_base")

    result = format_mopsmap_outputs(out_mopsmap, output_types=output_types)

    result_cache.put(key, result)
    return result


__all__ = ["compute_optical_properties"]
