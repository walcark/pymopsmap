"""Specie: a species description, and one materialised point of it."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple

import xarray as xr

from pymopsmap.models.microparams import MicroParameters
from pymopsmap.models.optiprops import OptiProps
from pymopsmap.models.output_request import DEFAULT_OUTPUT, OutputRequest
from pymopsmap.sweep import as_space
from pymopsmap.utils import check_within_grid

from .catalog import CamsSpecie, path_for
from .schema import (
    Growth,
    expected_psd_variables,
    expected_shape_variables,
    read,
)

# MOPSMAP refuses a relative humidity above this value
# (calc_hygroscopic_growth.f90).
MAX_ENGINE_RH = 99.0


class Point(NamedTuple):
    """One species evaluated at one point of its parameter space.

    Attributes
    ----------
    modes : list[MicroParameters]
        One validated set of micro-parameters per mode.
    engine_rh : float or None
        The humidity MOPSMAP must apply itself. ``None`` when the modes are
        already wet, which is the case for tabulated species: growing them a
        second time would double the effect.
    """

    modes: list[MicroParameters]
    engine_rh: float | None


@dataclass(frozen=True)
class Specie:
    """An aerosol species, backed by one group of a catalogue file."""

    name: str
    tree: xr.DataTree
    source: str
    version: str | None

    @property
    def growth(self) -> Growth:
        return Growth(self.tree.attrs["growth"])

    @property
    def modes(self) -> list[str]:
        return list(self.tree.children)

    @property
    def rh_range(self) -> tuple[float, float] | None:
        """Humidity covered by the table, or None when there is no table."""
        return self._range("rh")

    @property
    def wl_range(self) -> tuple[float, float]:
        span = self._range("wl")
        assert span is not None  # every mode carries a wavelength axis
        return span

    def _range(self, axis: str) -> tuple[float, float] | None:
        first = self.tree[self.modes[0]].to_dataset()
        if axis not in first.coords:
            return None
        values = first[axis].values
        return float(values.min()), float(values.max())

    def at(
        self,
        wl: list[float],
        rh: float | None = None,
        kappa: float | None = None,
    ) -> Point:
        """
        Materialise the species at one point of its parameter space.

        Parameters
        ----------
        wl : list of float
            Wavelengths in micrometres.
        rh : float, optional
            Relative humidity in percent. Required for a tabulated species,
            optional for a kappa one, refused for a dry one.
        kappa : float, optional
            Overrides the hygroscopicity stored in the file. Only meaningful
            for a species whose growth is driven by kappa.

        Returns
        -------
        Point
            The materialised modes, and the humidity left to the engine.
        """
        check_within_grid("wl", wl, self.wl_range, self._where)
        engine_rh = self._resolve_humidity(rh)
        selectors = {"wl": wl} | (
            {"rh": rh} if self.growth is Growth.TABULATED else {}
        )
        return Point(
            modes=[self._mode(name, selectors, kappa) for name in self.modes],
            engine_rh=engine_rh,
        )

    def compute(
        self,
        wl: list[float],
        rh: float | list[float] | None = None,
        kappa: float | None = None,
        outputs: OutputRequest = DEFAULT_OUTPUT,
        quiet: bool = False,
    ) -> OptiProps:
        """
        Compute the optical properties of the species.

        Parameters
        ----------
        wl : list of float
            Wavelengths in micrometres. MOPSMAP takes the whole grid in one
            run, so this is never a swept axis.
        rh : float or list of float, optional
            Relative humidity in percent. A list adds an ``rh`` dimension to
            the result, one MOPSMAP run per value.
        kappa : float, optional
            Overrides the hygroscopicity stored in the file.
        outputs : OutputRequest
            Which families of optical properties to compute.
        quiet : bool
            Silence the dataset download progress bars.

        Returns
        -------
        OptiProps
        """
        from pymopsmap import engine
        from pymopsmap.models import extend_optiprops

        points, dims = as_space(rh=rh)
        # Materialise every point first: an invalid request must fail before
        # any MOPSMAP run rather than halfway through a sweep.
        materialised = [
            self.at(wl=wl, rh=point["rh"], kappa=kappa) for point in points
        ]
        results = [
            engine.run_point(
                point.modes,
                output_types=outputs,
                rh=point.engine_rh,
                quiet=quiet,
            )
            for point in materialised
        ]

        if not dims:
            return results[0]
        index = [{dim: point[dim] for dim in dims} for point in points]
        return extend_optiprops(index=index, optiprops_li=results)

    @property
    def _where(self) -> str:
        version = f" {self.version}" if self.version else ""
        return f"{self.source}{version} {self.name}"

    def _resolve_humidity(self, rh: float | None) -> float | None:
        """Validate the requested humidity and say who must apply it."""
        if self.growth is Growth.NONE:
            if rh is not None:
                raise ValueError(
                    f"{self._where} does not model humidity "
                    f"(growth='none'), so rh={rh} cannot be honoured."
                )
            return None

        if self.growth is Growth.TABULATED:
            if rh is None:
                raise ValueError(
                    f"{self._where} is tabulated against humidity, so rh is "
                    "required."
                )
            span = self.rh_range
            assert span is not None  # the schema guarantees the rh axis
            check_within_grid("rh", rh, span, self._where)
            # The table already holds the wet state; the engine must not grow
            # the particles a second time.
            return None

        if rh is None:
            return None
        check_within_grid("rh", rh, (0.0, MAX_ENGINE_RH), self._where)
        return rh

    def _mode(
        self, name: str, selectors: dict, kappa: float | None
    ) -> MicroParameters:
        ds = self.tree[name].to_dataset().interp(**selectors)
        psd_type = ds.attrs["psd_type"]
        shape_type = ds.attrs["shape_type"]

        # Validated as a mapping so the discriminated unions resolve the
        # shape and psd types from their stored discriminator.
        return MicroParameters.model_validate(
            {
                "wavelength": list(selectors["wl"]),
                "n_real": ds["n_real"].values.tolist(),
                "n_imag": ds["n_imag"].values.tolist(),
                "shape": {"type": shape_type}
                | _fields(ds, expected_shape_variables(shape_type)),
                "psd": {"type": psd_type}
                | _fields(ds, expected_psd_variables(psd_type)),
                "kappa": _optional(ds, "kappa", kappa),
                "density": _optional(ds, "density_dry", None),
            }
        )


def _fields(ds: xr.Dataset, names: frozenset[str]) -> dict:
    """Read a set of variables, as scalars or as lists when they have dims."""
    return {
        name: ds[name].values.tolist() if ds[name].dims else float(ds[name])
        for name in names
    }


def _optional(
    ds: xr.Dataset, name: str, override: float | None
) -> float | None:
    if override is not None:
        return override
    return float(ds[name]) if name in ds.variables else None


def load(
    specie: CamsSpecie | str | Path, version: str | None = None
) -> Specie:
    """
    Load a species from the built-in catalogue or from a file.

    Parameters
    ----------
    specie : CamsSpecie or str or Path
        A catalogue entry, or the path to a species file holding exactly one
        species.
    version : str, optional
        Source version, for a catalogue entry only. Defaults to the most
        recent one shipped.

    Returns
    -------
    Specie
    """
    if isinstance(specie, CamsSpecie):
        path = path_for(specie, version)
        tree = read(path)
        if specie.value not in tree.children:
            raise KeyError(
                f"'{specie.value}' is not part of {tree.attrs['source']} "
                f"{tree.attrs['version']}."
            )
        name = specie.value
    else:
        tree = read(Path(specie))
        names = list(tree.children)
        if len(names) != 1:
            raise ValueError(
                f"{specie} holds several species {names}; load one from the "
                "catalogue instead of the file as a whole."
            )
        name = names[0]

    group = tree[name]
    assert isinstance(group, xr.DataTree)
    return Specie(
        name=name,
        tree=group,
        source=tree.attrs.get("source", "unknown"),
        version=tree.attrs.get("version"),
    )
