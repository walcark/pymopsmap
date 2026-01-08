"""
logging.py

Author  : Kévin Walcarius
Date    : 2026-01-08
Version : 1.0
License : MIT
Summary : Module used to define loggers. The configuration
          by default in defined in a YAML file. The configuration
          is ensured to by applied once and for all for each
          logger instanciated thoughrough the project.

"""

from __future__ import annotations

import logging
import logging.config
from pathlib import Path
import yaml


DEFAULT_LOGGING_CONFIG = (
    Path(__file__).resolve().parent.parent.parent.parent
    / "config"
    / "logging.yaml"
)

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
        raise FileNotFoundError(
            f"Logging config file not found: {config_file}"
        )

    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    logging.config.dictConfig(config)
    _LOGGING_INITIALIZED = True


def get_logger(name: str = "pymopsmap") -> logging.Logger:
    """
    Returns a logger configured according to the YAML file.
    If logging is not initialized yet, initializes it automatically.
    """
    if not _LOGGING_INITIALIZED:
        init_logging(DEFAULT_LOGGING_CONFIG)

    return logging.getLogger(name)
