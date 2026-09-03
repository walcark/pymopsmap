"""What the optical dataset at hand covers, and what it does not."""

from __future__ import annotations

from typing import TYPE_CHECKING

from pymopsmap.exceptions import CoverageError

if TYPE_CHECKING:
    from pymopsmap.models import MicroParameters

# The MOPSMAP dataset is published in two archives: a mandatory one covering
# refractive indices from 1.28 to 1.64, and an extended one covering the rest
# of Table 1 of Gasteiger and Wiegner (2018). Soot, mineral dust and organic
# matter reach outside the first one.
MAIN_ARCHIVE_REAL_RANGE = (1.28, 1.64)


def require_available(
    missing: list[str], modes: list[MicroParameters], source: str
) -> None:
    """
    Turn a missing dataset file into a message naming what asked for it.

    A missing file is not a download failure but a coverage gap: the
    refractive index is a legitimate grid point that this copy of the dataset
    does not ship.

    Parameters
    ----------
    missing : list of str
        Dataset files that could not be found.
    modes : list of MicroParameters
        The modes whose refractive indices required them.
    source : str
        Where the dataset was looked for.

    Raises
    ------
    CoverageError
        If anything is missing.
    """
    if not missing:
        return

    low, high = MAIN_ARCHIVE_REAL_RANGE
    outside = sorted(
        {
            round(value, 4)
            for mode in modes
            for value in mode.n_real  # type: ignore[union-attr]
            if not low <= value <= high
        }
    )

    lines = [
        f"The optical dataset at {source} does not ship "
        f"{len(missing)} required file(s), starting with {missing[0]}."
    ]
    if outside:
        lines.append(
            f"Refractive indices {outside} fall outside the range "
            f"[{low}, {high}] of the main archive; they need the extended one."
        )
    else:
        lines.append(
            "Every refractive index is within the main archive range "
            f"[{low}, {high}], so the dataset copy is incomplete."
        )
    raise CoverageError(" ".join(lines), missing=missing, source=source)
