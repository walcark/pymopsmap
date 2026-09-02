"""Optical dataset downloader — fetches NC files from a remote source."""

import shutil
import time
import urllib.request
from pathlib import Path

from tqdm import tqdm

from pymopsmap.exceptions import DatasetSourceNotConfiguredError, DownloadError
from pymopsmap.utils import DATASET_SOURCE, get_logger

from .cache import OpticalDatasetCache

logger = get_logger(__name__)

_MAX_RETRIES = 3
_RETRY_DELAYS = [1, 2, 4]


class DatasetDownloader:
    def __init__(
        self,
        cache: OpticalDatasetCache,
        source: str | Path | None = None,
        quiet: bool = False,
    ):
        self.cache = cache
        self.source = source or DATASET_SOURCE
        self.quiet = quiet

    def _resolved_source(self) -> str:
        if self.source is None:
            raise DatasetSourceNotConfiguredError(
                "PYMOPSMAP_DATASET_SOURCE is not configured. "
                "Set the env var or pass source= to DatasetDownloader."
            )
        return str(self.source)

    def download(self, relative_path: str) -> None:
        source = self._resolved_source()
        url_or_path = f"{source.rstrip('/')}/{relative_path}"
        tmp_path = self.cache.full_path(relative_path).with_suffix(".tmp")
        tmp_path.parent.mkdir(parents=True, exist_ok=True)

        last_exc: Exception = RuntimeError("no attempts made")
        for attempt in range(_MAX_RETRIES):
            try:
                if Path(source).exists() or not source.startswith("http"):
                    self._copy_local(
                        Path(url_or_path), tmp_path, relative_path
                    )
                else:
                    self._download_https(url_or_path, tmp_path, relative_path)
                self.cache.register(relative_path, tmp_path)
                return
            except (OSError, urllib.error.URLError) as exc:
                last_exc = exc
                if tmp_path.exists():
                    tmp_path.unlink(missing_ok=True)
                if attempt < _MAX_RETRIES - 1:
                    time.sleep(_RETRY_DELAYS[attempt])

        raise DownloadError(
            f"Failed to download {relative_path} after {_MAX_RETRIES} retries",
            file_path=relative_path,
            source=source,
            cause=last_exc,
        )

    def _copy_local(self, src: Path, tmp_path: Path, label: str) -> None:
        size = src.stat().st_size if src.exists() else 0
        with tqdm(
            total=size,
            unit="B",
            unit_scale=True,
            desc=label,
            disable=self.quiet,
        ) as bar:
            shutil.copy2(src, tmp_path)
            bar.update(size)

    def _download_https(self, url: str, tmp_path: Path, label: str) -> None:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req) as response:
            total = int(response.headers.get("Content-Length", 0)) or None
            with tqdm(
                total=total,
                unit="B",
                unit_scale=True,
                desc=label,
                disable=self.quiet,
            ) as bar:
                with open(tmp_path, "wb") as f:
                    while chunk := response.read(65536):
                        f.write(chunk)
                        bar.update(len(chunk))

    def download_missing(self, paths: list[str]) -> None:
        for p in paths:
            if not self.cache.is_cached(p):
                logger.debug("Downloading %s", p)
                self.download(p)
            else:
                logger.debug("Already cached: %s", p)
