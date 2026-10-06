"""The MOPSMAP optical dataset: resolve, download, cache.

Named after the scatlib command of the launch file, which is what MOPSMAP
itself calls this data.
"""

from .cache import OpticalDatasetCache
from .coverage import require_available
from .downloader import DatasetDownloader
from .resolver import NCFileResolver

__all__ = [
    "OpticalDatasetCache",
    "DatasetDownloader",
    "NCFileResolver",
    "require_available",
]
