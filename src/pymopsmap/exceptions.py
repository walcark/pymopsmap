"""Custom exceptions for pymopsmap."""


class DatasetSourceNotConfiguredError(Exception):
    """Raised when PYMOPSMAP_DATASET_SOURCE is not set and a download is needed."""


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
