"""The data set index is read once, not once per point of a sweep."""

from __future__ import annotations

import xarray as xr


def _index(path):
    xr.Dataset(
        {"max_sizepara": (("np", "eps", "mreal", "mimag"), [[[[1005.0]]]])},
        coords={
            "np": [-1],
            "eps": [1.0],
            "mreal": [1.52],
            "mimag": [0.0],
        },
    ).to_netcdf(path)
    return path


def test_the_limits_open_the_file_once_per_path(tmp_path, monkeypatch):
    """
    netCDF4 is not thread safe, and a sweep walks its points in threads.

    Opening the index per point raced and raised "NetCDF: Not a valid ID";
    reading it once per path removes the race and the work.
    """
    from pymopsmap.scatlib import limits as limits_module

    path = _index(tmp_path / "index.nc")
    opened = []
    original = xr.open_dataset

    def counting(*args, **kwargs):
        opened.append(args[0])
        return original(*args, **kwargs)

    monkeypatch.setattr(limits_module.xr, "open_dataset", counting)
    limits_module.read_index.cache_clear()

    limits_module.SizeParameterLimits(path)
    limits_module.SizeParameterLimits(path)
    limits_module.SizeParameterLimits(path)

    assert len(opened) == 1


def test_the_resolver_opens_the_file_once_per_path(tmp_path, monkeypatch):
    from pymopsmap.scatlib import limits as limits_module
    from pymopsmap.scatlib import resolver as resolver_module

    path = _index(tmp_path / "index.nc")
    opened = []
    original = xr.open_dataset

    def counting(*args, **kwargs):
        opened.append(args[0])
        return original(*args, **kwargs)

    monkeypatch.setattr(limits_module.xr, "open_dataset", counting)
    limits_module.read_index.cache_clear()

    resolver_module.NCFileResolver(path)
    resolver_module.NCFileResolver(path)

    assert len(opened) == 1
