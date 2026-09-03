"""MOPSMAP computation pipeline — resolution, execution, caching."""

from __future__ import annotations

from typing import TYPE_CHECKING, TypeAlias

if TYPE_CHECKING:
    import xarray as xr

    from pymopsmap.engine.outputs import OutputRequest
    from pymopsmap.microparams import MicroParameters

    Modes: TypeAlias = MicroParameters | list[MicroParameters]


def run_point(
    modes: Modes,
    output_types: OutputRequest,
    rh: float | None = None,
    quiet: bool = False,
) -> xr.Dataset:
    """
    Run MOPSMAP once, for one point of a parameter space.

    Transparent pipeline: check the result cache, ensure the required dataset
    files are present, write the launch file, run the binary, parse the
    outputs, store the result.
    """
    from pymopsmap.exceptions import DownloadError
    from pymopsmap.scatlib.cache import OpticalDatasetCache
    from pymopsmap.scatlib.coverage import require_available
    from pymopsmap.scatlib.downloader import DatasetDownloader
    from pymopsmap.scatlib.resolver import NCFileResolver
    from pymopsmap.scatlib.results import ResultCache
    from pymopsmap.utils import DATASET_CACHE_DIR

    from .coverage import clip_modes_to_coverage, reindex_to_full_grid
    from .launch_file import write_launching_file
    from .launcher import launch_mopsmap
    from .output_format import format_mopsmap_outputs, shape_types

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
    modes_run, valid_mask = clip_modes_to_coverage(modes, rh=rh)

    # Ensure index.nc is present first
    index_path = dataset_cache.full_path("index.nc")
    if not dataset_cache.is_cached("index.nc"):
        downloader.download("index.nc")

    # Resolve required dataset files and download missing ones. A file the
    # source does not ship is a coverage gap, not a download failure.
    resolver = NCFileResolver(index_path)
    mp_list = [modes_run] if not isinstance(modes_run, list) else modes_run
    required = resolver.resolve(mp_list, rh=rh)
    try:
        downloader.download_missing(required)
    except DownloadError as exc:
        require_available(
            missing=[exc.file_path], modes=mp_list, source=exc.source
        )
        raise

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
    # Some combination rules need to know what produced the result: an
    # effective radius cannot be rebuilt from non-spherical modes.
    result.attrs["shape_types"] = shape_types(mp_list)

    # Reindex to the original wavelength grid; clipped positions become NaN.
    if not valid_mask.all():
        result = reindex_to_full_grid(result, original_wl, valid_mask)

    result_cache.put(key, result)
    return result


__all__ = ["run_point"]
