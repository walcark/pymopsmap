"""Mode: one aerosol mode described in Python rather than read from a file."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import xarray as xr

from pymopsmap.psd import PSD
from pymopsmap.shapes import Shape

from .schema import AMPLITUDE_FIELD


@dataclass(frozen=True)
class Mode:
    """
    One mode of a hand-built species.

    Any numeric parameter accepts a ``DataArray`` instead of a float, which
    sweeps it: distinct dimensions multiply, and two parameters sharing a
    dimension vary together. The dimensions live on the species rather than on
    the compute call, so a parameter name is never ambiguous between modes.

    Parameters
    ----------
    shape : Shape
        Particle shape.
    psd : PSD
        Size distribution, amplitude included.
    n_real, n_imag : float or DataArray
        Refractive index, on the wavelengths given to ``compute``.
    density_dry : float, optional
        Dry material density in g cm-3.
    kappa : float, optional
        Hygroscopicity. Its presence makes the species grow with humidity.
    sweep : dict, optional
        Size distribution or shape parameters to vary, as ``DataArray``.
        They are given here rather than inside the ``psd``, whose fields stay
        typed as floats so that a materialised point is always validated.
    """

    shape: Shape
    psd: PSD
    n_real: Any
    n_imag: Any
    density_dry: float | None = None
    kappa: float | None = None
    sweep: dict[str, xr.DataArray] = field(default_factory=dict)
    name: str = field(default="only")

    def to_dataset(self, wl: list[float] | None = None) -> xr.Dataset:
        """Render the mode in the canonical schema."""
        variables: dict[str, Any] = {
            "n_real": _as_variable(self.n_real),
            "n_imag": _as_variable(self.n_imag),
        }
        for name in _fields_of(self.psd):
            variables[name] = _as_variable(getattr(self.psd, name))
        for name in _fields_of(self.shape):
            variables[name] = _as_variable(getattr(self.shape, name))
        variables.update(self.sweep)
        if self.density_dry is not None:
            variables["density_dry"] = self.density_dry
        if self.kappa is not None:
            variables["kappa"] = self.kappa

        coords = {"wl": wl} if wl is not None else {}
        ds = xr.Dataset(variables, coords=coords)
        ds.attrs.update(psd_type=self.psd.type, shape_type=self.shape.type)
        _stamp_units(ds)
        return ds


def _fields_of(model) -> list[str]:
    return [name for name in type(model).model_fields if name != "type"]


def _as_variable(value: Any) -> Any:
    """Keep a DataArray as it is; anything else becomes a scalar or a list."""
    if isinstance(value, xr.DataArray):
        return value
    if isinstance(value, (list, tuple, np.ndarray)):
        return ("wl", np.asarray(value, dtype=float))
    return value


LENGTHS = {"rm", "rmin", "rmax", "radius"}


def _stamp_units(ds: xr.Dataset) -> None:
    for name in ds.variables:
        key = str(name)
        if key in LENGTHS or key == "wl":
            ds[key].attrs.update(units="um")
        elif key in AMPLITUDE_FIELD.values():
            ds[key].attrs.update(units="m-3")
        elif key == "density_dry":
            ds[key].attrs.update(units="g cm-3")
        elif key in ("n_real", "n_imag", "sigma", "kappa"):
            ds[key].attrs.update(units="1")
    if "n_imag" in ds:
        ds["n_imag"].attrs.update(sign="positive")
