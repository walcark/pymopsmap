"""The MOPSMAP optical dataset: resolve, download, cache.

Named after the scatlib command of the launch file, which is what MOPSMAP
itself calls this data.
"""

from .cache import OpticalDatasetCache
from .downloader import DatasetDownloader
from .resolver import NCFileResolver
from .results import ResultCache

__all__ = [
    "OpticalDatasetCache",
    "DatasetDownloader",
    "NCFileResolver",
    "ResultCache",
]
