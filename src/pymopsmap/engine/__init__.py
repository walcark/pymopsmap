"""MOPSMAP computation pipeline — resolution, execution, caching."""

from __future__ import annotations

from typing import TYPE_CHECKING, TypeAlias

if TYPE_CHECKING:
    from pymopsmap.models import MicroParameters, OptiProps
    from pymopsmap.models.output_request import OutputRequest

    Modes: TypeAlias = MicroParameters | list[MicroParameters]


def run_point(
    modes: Modes,
    output_types: OutputRequest,
    rh: float | None = None,
    quiet: bool = False,
) -> OptiProps:
    """
    Run MOPSMAP once, for one point of a parameter space.

    Transparent pipeline: check the result cache, ensure the required dataset
    files are present, write the launch file, run the binary, parse the
    outputs, store the result.
    """
    from pymopsmap.cache.downloader import DatasetDownloader
    from pymopsmap.cache.optical import OpticalDatasetCache
    from pymopsmap.cache.resolver import NCFileResolver
    from pymopsmap.cache.results import ResultCache
    from pymopsmap.models import OptiProps
    from pymopsmap.utils import DATASET_CACHE_DIR

    from .coverage import clip_modes_to_coverage, reindex_to_full_grid
    from .launch_file import write_launching_file
    from .launcher import launch_mopsmap
    from .output_format import format_mopsmap_outputs

    result_cache = ResultCache()
    dataset_cache = OpticalDatasetCache()
    downloader = DatasetDownloader(cache=dataset_cache, quiet=quiet)

    # Cache key is based on the original (unclipped) request so the full masked
    # result is stored and returned on subsequent calls.
    key = result_cache.key(modes, output_types, rh=rh)
    cached = result_cache.get(key)
    if cached is not None:
        return cached

    # Clip wavelengths that exceed dataset size-parameter coverage.
    # original_wl is used to reindex the result back to the full grid.
    mp_ref = modes if not isinstance(modes, list) else modes[0]
    original_wl = list(mp_ref.wavelength)
    modes_run, valid_mask = clip_modes_to_coverage(modes)

    # Ensure index.nc is present first
    index_path = dataset_cache.full_path("index.nc")
    if not dataset_cache.is_cached("index.nc"):
        downloader.download("index.nc")

    # Resolve required dataset files and download missing ones
    resolver = NCFileResolver(index_path)
    mp_list = [modes_run] if not isinstance(modes_run, list) else modes_run
    required = resolver.resolve(mp_list)
    downloader.download_missing(required)

    # Run MOPSMAP on the (possibly clipped) wavelength grid
    paths = write_launching_file(
        mp=modes_run,
        output_types=output_types,
        rh=rh,
        mopsmap_data_path=DATASET_CACHE_DIR,
    )
    out_mopsmap = launch_mopsmap(input_filename=paths["mopsmap"])
    out_mopsmap["ascii_base"] = paths.get("ascii_base")

    result = format_mopsmap_outputs(out_mopsmap, output_types=output_types)

    # Reindex to the original wavelength grid; clipped positions become NaN.
    if not valid_mask.all():
        result = OptiProps(
            ds=reindex_to_full_grid(result.ds, original_wl, valid_mask)
        )

    result_cache.put(key, result)
    return result


__all__ = ["run_point"]
