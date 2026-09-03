"""The built-in species catalogue, shipped inside the package."""

from __future__ import annotations

import importlib.resources as resources
from enum import StrEnum
from pathlib import Path

DATA_PACKAGE = "pymopsmap"

# Sorting version names is not an ordering: it happens to work for CAMS
# (47r1 < 48r1 < 49r1) and would pick 'kappa' over 'geisa' for OPAC. The
# default is therefore stated, not derived.
DEFAULT_VERSION: dict[str, str] = {"cams": "49r1", "opac": "kappa"}


class Source(StrEnum):
    """A catalogue source, one directory of version files."""

    CAMS = "cams"
    OPAC = "opac"


class CamsSpecie(StrEnum):
    """CAMS aerosol species, named as in the catalogue files."""

    AMMONIUM = "ammonium"
    BLACK_CARBON = "black_carbon"
    CONTINENTAL = "continen"
    DUST = "dust"
    NITRATE = "nitrate"
    ORGANIC_MATTER = "organic_matter"
    SEA_SALT = "sea_salt"
    SECONDARY_ORGANIC = "secondary_organic"
    SULPHATE = "sulphate"

    @property
    def source(self) -> str:
        return Source.CAMS.value

    @property
    def versions(self) -> tuple[str, ...]:
        """Every version shipped for this source."""
        return versions_for(self.source)


class OpacSpecie(StrEnum):
    """
    OPAC aerosol components, under the codes Hess et al. (1998) gave them.

    INSO insoluble, WASO water soluble, SOOT soot, SSAM and SSCM the
    accumulation and coarse sea salt modes, SUSO sulfate droplets, MINM MIAM
    MICM the nucleation, accumulation and coarse mineral modes, MITR
    transported mineral.
    """

    INSO = "inso"
    WASO = "waso"
    SOOT = "soot"
    SSAM = "ssam"
    SSCM = "sscm"
    SUSO = "suso"
    MINM = "minm"
    MIAM = "miam"
    MICM = "micm"
    MITR = "mitr"

    @property
    def source(self) -> str:
        return Source.OPAC.value

    @property
    def versions(self) -> tuple[str, ...]:
        """Every version shipped for this source."""
        return versions_for(self.source)


CatalogSpecie = CamsSpecie | OpacSpecie


def _source_dir(source: str) -> Path:
    return Path(str(resources.files(DATA_PACKAGE) / "data" / source))


def versions_for(source: str) -> tuple[str, ...]:
    """Versions shipped for a source, oldest first."""
    return tuple(sorted(p.stem for p in _source_dir(source).glob("*.nc")))


def path_for(specie: CatalogSpecie, version: str | None = None) -> Path:
    """
    Locate the catalogue file holding a species.

    Parameters
    ----------
    specie : CamsSpecie or OpacSpecie
        The species to locate.
    version : str, optional
        Source version. Defaults to the most recent one shipped.

    Raises
    ------
    ValueError
        If the requested version is not shipped.
    """
    available = versions_for(specie.source)
    if version is None:
        version = DEFAULT_VERSION[specie.source]
    if version not in available:
        raise ValueError(
            f"Unknown {specie.source} version '{version}'. "
            f"Available: {list(available)}."
        )
    return _source_dir(specie.source) / f"{version}.nc"
