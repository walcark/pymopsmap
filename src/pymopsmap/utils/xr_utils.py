from typing import Iterable, TypeVar, Union, Any
from numpy.typing import NDArray
import xarray as xr
import numpy as np


T = TypeVar("T")
U = TypeVar("U", bound=np.generic)


def ensure_np_array_1d(
    data: Union[T, Iterable[T], NDArray[T]],
    dtype: type[U] | None = None,
    rounding: int | None = None
) -> NDArray[Any]:
    """
    Ensure the input is converted to a 1D Numpy array of type T.
    """
    if dtype is None:
        arr = np.asarray(data)
    else:
        arr = np.asarray(data, dtype=dtype)

    arr = np.atleast_1d(arr)

    if arr.ndim != 1:
        raise ValueError(f"Expected 1D array, got shape {arr.shape}.")

    if round is not None:
        arr = np.round(arr, rounding)

    return arr


def ensure_list_1d(
    data: Union[T, Iterable[T], NDArray[T]],
    dtype: type[U] | None = None,
    rounding: int | None = None
) -> list[Any]:
    """
    Ensure the input is converted to a python list of length >= 1.
    """
    if isinstance(data, list):
        return data

    arr = ensure_np_array_1d(data, dtype=dtype, rounding=rounding)
    return arr.tolist()


def extract_from_dataset(
    xrds: xr.Dataset,
    variable: str,
    attrs: dict[str, Union[T, Iterable[T], NDArray[T]]]
) -> xr.DataArray:
    """
    Extract a variable from a xr.Dataset.
    """
    attrs = {k: v for k, v in attrs.items()}
    xrds: xr.DataArray = xrds[variable].sel(**attrs)
    return xrds
