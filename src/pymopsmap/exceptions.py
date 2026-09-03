"""Custom exceptions for pymopsmap."""


class DatasetSourceNotConfiguredError(Exception):
    """Raised when PYMOPSMAP_DATASET_SOURCE is not set."""


class DownloadError(Exception):
    def __init__(
        self, message: str, file_path: str, source: str, cause: Exception
    ):
        super().__init__(message)
        self.file_path = file_path
        self.source = source
        self.cause = cause


class MopsmapError(Exception):
    def __init__(self, message: str, returncode: int, stderr: str):
        super().__init__(message)
        self.returncode = returncode
        self.stderr = stderr


class IndexFileError(Exception):
    """Raised when index.nc is missing or unreadable."""


class DomainError(ValueError):
    """Raised when an interpolation target falls outside its source grid."""

    def __init__(
        self,
        axis: str,
        requested: float,
        available: tuple[float, float],
        context: str | None = None,
    ):
        where = f" for {context}" if context else ""
        super().__init__(
            f"{axis}={requested:g} is outside the available range "
            f"[{available[0]:g}, {available[1]:g}]{where}."
        )
        self.axis = axis
        self.requested = requested
        self.available = available
        self.context = context


class SchemaError(ValueError):
    """Raised when a species file does not follow the canonical schema."""


class CoverageError(Exception):
    """Raised when the optical dataset at hand lacks a required file."""

    def __init__(self, message: str, missing: list[str], source: str):
        super().__init__(message)
        self.missing = missing
        self.source = source
