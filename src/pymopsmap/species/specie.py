"""Specie: a species description, and one materialised point of it."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, NamedTuple

import numpy as np
import xarray as xr

from pymopsmap.engine.outputs import DEFAULT_OUTPUT, OutputRequest
from pymopsmap.microparams import MicroParameters
from pymopsmap.sweep import build_space, distinct_values, run_sweep
from pymopsmap.utils import check_within_grid

from .catalog import CamsSpecie, CatalogSpecie, OpacSpecie, path_for
from .mode import Mode
from .schema import (
    AMPLITUDE_FIELD,
    HUMIDITY_DIM,
    Growth,
    expected_psd_variables,
    expected_shape_variables,
    read,
)

# MOPSMAP refuses a relative humidity above this value
# (calc_hygroscopic_growth.f90).
MAX_ENGINE_RH = 99.0


class CacheStatusReport(NamedTuple):
    """Which optical dataset files are already here, and which are not."""

    cached: list[str]
    missing: list[str]


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


# eq=False keeps identity semantics, so a Specie can key a mixture: its tree
# is not hashable and two loads of the same file are two descriptions.
@dataclass(frozen=True, eq=False)
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
    def amplitude(self) -> float:
        """
        Total stored amplitude, summed over the modes.

        It is the reference a mixture scales against: dividing by it and
        multiplying by a requested concentration preserves the ratio between
        the modes, which a plain override could not express.
        """
        total = 0.0
        for name in self.modes:
            ds = self.tree[name].to_dataset()
            field = AMPLITUDE_FIELD[ds.attrs["psd_type"]]
            if field is None:
                raise ValueError(
                    f"{self._where}: mode '{name}' keeps its amplitude in an "
                    "external file, so it cannot be scaled."
                )
            total += float(np.sum(ds[field].values))
        return total

    @property
    def rh_range(self) -> tuple[float, float] | None:
        """Humidity covered by the table, or None when there is no table."""
        return self._range("rh")

    @property
    def wl_range(self) -> tuple[float, float] | None:
        """
        Wavelengths the tabulated refractive index covers.

        None when the index is a scalar, which covers every wavelength.
        """
        return self._range("wl")

    def _range(self, axis: str) -> tuple[float, float] | None:
        first = self.tree[self.modes[0]].to_dataset()
        if axis not in first.coords:
            return None
        values = first[axis].values
        return float(values.min()), float(values.max())

    @classmethod
    def custom(
        cls,
        modes: Mode | list[Mode],
        name: str = "custom",
        wl: list[float] | None = None,
    ) -> Specie:
        """
        Build a species in Python rather than reading one from a file.

        Parameters
        ----------
        modes : Mode or list of Mode
            The modes, in order. Their names default to ``only`` for a single
            mode, and to ``mode_1``, ``mode_2``, ... for several.
        name : str
            Name of the species, used in error messages and when saving.
        wl : list of float, optional
            Wavelengths the refractive index is given on. Omit it when the
            index is a scalar.

        Returns
        -------
        Specie
        """
        listed = [modes] if isinstance(modes, Mode) else list(modes)
        tree = xr.DataTree()
        tree.attrs["source"] = "custom"
        for index, mode in enumerate(listed, start=1):
            label = mode.name if len(listed) == 1 else f"mode_{index}"
            tree[label] = xr.DataTree(mode.to_dataset(wl))
        growth = (
            Growth.KAPPA
            if any(mode.kappa is not None for mode in listed)
            else Growth.NONE
        )
        tree.attrs["growth"] = growth.value
        return cls(name=name, tree=tree, source="custom", version=None)

    @property
    def swept(self) -> dict[str, int]:
        """
        Dimensions a parameter of this species varies over, and their sizes.

        Wavelength and humidity are not among them: both are arguments of
        ``compute``, and the humidity axis of a tabulated species is the table
        it is interpolated on, not a sweep the caller asked for.
        """
        sizes: dict[str, int] = {}
        for mode in self.modes:
            for dim, size in self.tree[mode].to_dataset().sizes.items():
                if str(dim) not in ("wl", HUMIDITY_DIM):
                    sizes[str(dim)] = size
        return sizes

    def to_netcdf(self, path: str | Path) -> None:
        """
        Write the species in the catalogue format.

        A hand-built species and a catalogue one serialise identically, which
        is what makes the format one thing rather than two.
        """
        from .schema import write

        tree = xr.DataTree()
        tree.attrs["source"] = self.source
        if self.version is not None:
            tree.attrs["version"] = self.version
        tree[self.name] = self.tree.copy()
        write(tree, path)

    def at(
        self,
        wl: list[float],
        rh: float | None = None,
        kappa: float | None = None,
        **swept: int,
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
        **swept : int
            One index per swept dimension of the species, naming the point to
            materialise. See ``Specie.swept``.

        Returns
        -------
        Point
            The materialised modes, and the humidity left to the engine.
        """
        if self.wl_range is not None:
            check_within_grid("wl", wl, self.wl_range, self._where)
        engine_rh = self._resolve_humidity(rh)
        selectors = {"wl": wl} | (
            {"rh": rh} if self.growth is Growth.TABULATED else {}
        )
        return Point(
            modes=[
                self._mode(name, selectors, kappa, swept)
                for name in self.modes
            ],
            engine_rh=engine_rh,
        )

    def compute(
        self,
        wl: list[float],
        rh: float | list[float] | None = None,
        kappa: float | None = None,
        outputs: OutputRequest = DEFAULT_OUTPUT,
        n_angles: int = 2000,
        quiet: bool = False,
    ) -> xr.Dataset:
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
        xr.Dataset
            One variable per optical property, with a dimension per swept axis.
        """
        from pymopsmap import engine
        from pymopsmap.engine.outputs import shape_types, variables_for

        # A parameter the species declares as varying is an axis of the
        # sweep too. It is walked by position, since what varies is a field
        # the caller shaped, not a value the engine can be told directly.
        declared = self.swept
        space, swept, fixed = build_space(
            wl,
            rh=rh,
            **{dim: np.arange(size) for dim, size in declared.items()},
        )
        # Materialise every point first: an invalid request must fail up
        # front, with its own error, rather than midway through a sweep
        # wrapped in the engine's.
        for humidity in distinct_values(rh):
            self.at(
                wl=wl,
                rh=humidity,
                kappa=kappa,
                **dict.fromkeys(declared, 0),
            )

        def point(**at: Any) -> xr.Dataset:
            wavelengths = list(at.pop("wl"))
            indices = {dim: int(at.pop(dim)) for dim in declared}
            materialised = self.at(
                wl=wavelengths, kappa=kappa, **at, **indices
            )
            return engine.run_point(
                materialised.modes,
                output_types=outputs,
                rh=materialised.engine_rh,
                quiet=quiet,
                n_angles=n_angles,
            )

        result = run_sweep(
            point,
            space,
            outputs=variables_for(outputs),
            version=self._sweep_version(outputs, n_angles),
            swept=swept,
            fixed=fixed,
            quiet=quiet,
        )
        # The store keeps the numbers, not the attributes around them, so
        # what produced them is restated here. A mixture needs it: an
        # effective radius cannot be rebuilt from non-spherical modes.
        result = self._label_declared_axes(result, declared)
        result.attrs["shape_types"] = shape_types(
            self.at(
                wl=wl,
                rh=distinct_values(rh)[0],
                kappa=kappa,
                **dict.fromkeys(declared, 0),
            ).modes
        )
        return result

    def _label_declared_axes(
        self, result: xr.Dataset, declared: dict[str, int]
    ) -> xr.Dataset:
        """
        Put the parameter values back on the axes swept by position.

        A declared axis is walked by index, so the result would otherwise
        carry 0, 1, 2 where the caller gave radii. The values are recovered
        from whichever variable varies along that dimension alone; a
        parameter spread over several dimensions keeps its index, since no
        single value belongs to a position.
        """
        for dim in declared:
            values = self._axis_values(dim)
            if values is not None:
                result = result.assign_coords({dim: values})
        return result

    def _axis_values(self, dim: str) -> Any:
        """The values a declared dimension takes, when it names them alone."""
        for name in self.modes:
            for variable in self.tree[name].to_dataset().values():
                if variable.dims == (dim,):
                    return variable.values
        return None

    def _sweep_version(self, outputs: OutputRequest, n_angles: int) -> str:
        """
        What makes this computation different from another one.

        It keys the store, so it has to name the species and everything the
        stored numbers depend on but the swept point: leave one out and two
        species share an entry.
        """
        from .schema import SCHEMA_REV

        requested = "-".join(sorted(o.value for o in outputs))
        # The angle count sets the shape of every angular output, and the
        # contract does not name it, so two counts are two stores.
        angular = f"theta{n_angles}" if _is_angular(outputs) else ""
        return "-".join(
            part
            for part in (
                self.source,
                self.version,
                self.name,
                f"schema{SCHEMA_REV}",
                self.fingerprint(),
                requested,
                angular,
            )
            if part
        )

    def fingerprint(self) -> str:
        """
        A digest of the description the stored numbers were computed from.

        Naming the species is not enough to key a store. A catalogue file can
        be rebuilt under the same source and version, and a hand-built species
        can be edited between two calls; either way the old results are no
        longer answers to the question being asked.
        """
        digest = hashlib.blake2b(digest_size=8)
        for node in self.tree.subtree:
            digest.update(node.path.encode())
            for key, value in sorted(node.attrs.items()):
                digest.update(f"{key}={value}".encode())
            dataset = node.to_dataset()
            for name in sorted(map(str, dataset.variables)):
                values = np.asarray(dataset[name].values)
                digest.update(name.encode())
                digest.update(str(values.dtype).encode())
                digest.update(values.tobytes())
        return digest.hexdigest()

    def cache_status(
        self,
        wl: list[float],
        rh: float | list[float] | None = None,
        kappa: float | None = None,
    ) -> CacheStatusReport:
        """
        Which optical dataset files this computation needs, and which are here.

        Parameters
        ----------
        wl : list of float
            Wavelengths in micrometres.
        rh : float or list of float, optional
            Relative humidity in percent.
        kappa : float, optional
            Overrides the hygroscopicity stored in the file.

        Returns
        -------
        CacheStatusReport
        """
        from pymopsmap.scatlib import cache as cache_module

        dataset_cache = cache_module.OpticalDatasetCache()
        required = self._required_files(wl, rh, kappa)
        return CacheStatusReport(
            cached=[p for p in required if dataset_cache.is_cached(p)],
            missing=[p for p in required if not dataset_cache.is_cached(p)],
        )

    def prefetch(
        self,
        wl: list[float],
        rh: float | list[float] | None = None,
        kappa: float | None = None,
        quiet: bool = False,
    ) -> None:
        """Download every dataset file this needs, without running."""
        from pymopsmap.scatlib import cache as cache_module
        from pymopsmap.scatlib import downloader as downloader_module

        dataset_cache = cache_module.OpticalDatasetCache()
        downloader = downloader_module.DatasetDownloader(
            cache=dataset_cache, quiet=quiet
        )
        downloader.download_missing(self._required_files(wl, rh, kappa))

    def _required_files(
        self,
        wl: list[float],
        rh: float | list[float] | None,
        kappa: float | None,
    ) -> list[str]:
        """Union of the dataset files needed over every point of the sweep."""
        from pymopsmap.scatlib import cache as cache_module
        from pymopsmap.scatlib import downloader as downloader_module
        from pymopsmap.scatlib import resolver as resolver_module

        dataset_cache = cache_module.OpticalDatasetCache()
        index = dataset_cache.full_path("index.nc")
        if not dataset_cache.is_cached("index.nc"):
            downloader_module.DatasetDownloader(cache=dataset_cache).download(
                "index.nc"
            )
        resolver = resolver_module.NCFileResolver(index)

        required: set[str] = set()
        for value in distinct_values(rh):
            materialised = self.at(wl=wl, rh=value, kappa=kappa)
            required.update(
                resolver.resolve(materialised.modes, rh=materialised.engine_rh)
            )
        return sorted(required)

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
        self,
        name: str,
        selectors: dict,
        kappa: float | None,
        swept: dict[str, int] | None = None,
    ) -> MicroParameters:
        ds = self.tree[name].to_dataset()
        if swept:
            ds = ds.isel({k: v for k, v in swept.items() if k in ds.dims})
        ds = ds.interp(**{k: v for k, v in selectors.items() if k in ds.dims})
        psd_type = ds.attrs["psd_type"]
        shape_type = ds.attrs["shape_type"]

        # Validated as a mapping so the discriminated unions resolve the
        # shape and psd types from their stored discriminator.
        return MicroParameters.model_validate(
            {
                "wavelength": list(selectors["wl"]),
                "n_real": _spectrum(ds["n_real"], len(selectors["wl"])),
                "n_imag": _spectrum(ds["n_imag"], len(selectors["wl"])),
                "shape": {"type": shape_type}
                | _fields(ds, expected_shape_variables(shape_type)),
                "psd": {"type": psd_type}
                | _fields(ds, expected_psd_variables(psd_type)),
                "kappa": _optional(ds, "kappa", kappa),
                "density": _optional(ds, "density_dry", None),
                "nonabs_fraction": _optional(ds, "nonabs_fraction", None)
                or 0.0,
                "size_equ": ds.attrs.get("size_equ", "cs"),
            }
        )


def _is_angular(outputs: OutputRequest) -> bool:
    """Whether any requested output is sampled over scattering angle."""
    from pymopsmap.engine.outputs import OutputType

    return bool(
        outputs
        & {
            OutputType.PHASE_FUNCTION,
            OutputType.SCATTERING_MATRIX,
            OutputType.VOLUME_SCATTERING_FUNCTION,
        }
    )


def _spectrum(values: xr.DataArray, count: int) -> list[float]:
    """A refractive index given as a scalar covers every wavelength."""
    if values.dims:
        return values.values.tolist()
    return [float(values)] * count


def _fields(ds: xr.Dataset, names: frozenset[str]) -> dict:
    """
    Read a set of variables back into the values a model expects.

    A variable with dimensions becomes a list, a numeric scalar a float, and
    anything else is handed over as it is: a shape carries strings too, such
    as the oblate or prolate of a spheroid.
    """
    values = {}
    for name in names:
        variable = ds[name]
        if variable.dims:
            values[name] = variable.values.tolist()
        elif variable.dtype.kind in "fiu":
            values[name] = float(variable)
        else:
            values[name] = variable.item()
    return values


def _optional(
    ds: xr.Dataset, name: str, override: float | None
) -> float | None:
    if override is not None:
        return override
    return float(ds[name]) if name in ds.variables else None


def load(
    specie: CatalogSpecie | str | Path, version: str | None = None
) -> Specie:
    """
    Load a species from the built-in catalogue or from a file.

    Parameters
    ----------
    specie : CamsSpecie or OpacSpecie or str or Path
        A catalogue entry, or the path to a species file holding exactly one
        species.
    version : str, optional
        Source version, for a catalogue entry only. Defaults to the most
        recent one shipped.

    Returns
    -------
    Specie
    """
    if isinstance(specie, (CamsSpecie, OpacSpecie)):
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
