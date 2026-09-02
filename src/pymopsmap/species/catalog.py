"""The built-in species catalogue, shipped inside the package."""

from __future__ import annotations

import importlib.resources as resources
from enum import StrEnum
from pathlib import Path

DATA_PACKAGE = "pymopsmap"


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
        return "cams"

    @property
    def versions(self) -> tuple[str, ...]:
        """Every version shipped for this source, oldest first."""
        return versions_for("cams")


def _source_dir(source: str) -> Path:
    return Path(str(resources.files(DATA_PACKAGE) / "data" / source))


def versions_for(source: str) -> tuple[str, ...]:
    """Versions shipped for a source, oldest first."""
    return tuple(sorted(p.stem for p in _source_dir(source).glob("*.nc")))


def path_for(specie: CamsSpecie, version: str | None = None) -> Path:
    """
    Locate the catalogue file holding a species.

    Parameters
    ----------
    specie : CamsSpecie
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
        return _source_dir(specie.source) / f"{available[-1]}.nc"
    if version not in available:
        raise ValueError(
            f"Unknown {specie.source} version '{version}'. "
            f"Available: {list(available)}."
        )
    return _source_dir(specie.source) / f"{version}.nc"
