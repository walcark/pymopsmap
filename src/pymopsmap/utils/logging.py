from __future__ import annotations

import logging
import logging.config
from pathlib import Path
import yaml
from typing import Optional


# Emplacement par défaut du fichier YAML
DEFAULT_LOGGING_CONFIG = Path(__file__).resolve().parent.parent.parent.parent / "config" / "logging.yaml"

_LOGGING_INITIALIZED = False


def init_logging(config_file: str | Path = DEFAULT_LOGGING_CONFIG) -> None:
    """
    Initialize Python logging from a YAML configuration file.
    Safe: calling it multiple times has no effect.
    """
    global _LOGGING_INITIALIZED

    if _LOGGING_INITIALIZED:
        return

    config_file = Path(config_file)
    if not config_file.exists():
        raise FileNotFoundError(f"Logging config file not found: {config_file}")

    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    logging.config.dictConfig(config)
    _LOGGING_INITIALIZED = True


def get_logger(
    name: str = "pymopsmap",
    config_file: Optional[str | Path] = None,
) -> logging.Logger:
    """
    Returns a logger configured according to the YAML file.
    If logging is not initialized yet, initializes it automatically.
    """

    # Initialize logging only once
    if not _LOGGING_INITIALIZED:
        init_logging(config_file or DEFAULT_LOGGING_CONFIG)

    return logging.getLogger(name)
