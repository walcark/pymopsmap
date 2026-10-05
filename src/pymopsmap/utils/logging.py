"""Structured logging via structlog: colored console, optional JSON file."""

from __future__ import annotations

import logging
import os
import sys

import structlog

_INITIALIZED = False

# A library has no business flooding the console of whoever imports it, nor
# dropping a log file in their working directory: both are opt-in.
_DEFAULT_LEVEL = "WARNING"


def _configure() -> None:
    global _INITIALIZED
    if _INITIALIZED:
        return

    shared_processors: list = [  # type: ignore[type-arg]
        structlog.stdlib.add_logger_name,
        structlog.stdlib.add_log_level,
        structlog.stdlib.PositionalArgumentsFormatter(),
        structlog.processors.TimeStamper(fmt="%Y-%m-%d %H:%M:%S"),
        structlog.processors.StackInfoRenderer(),
    ]

    structlog.configure(
        processors=[
            *shared_processors,
            structlog.stdlib.ProcessorFormatter.wrap_for_formatter,
        ],
        logger_factory=structlog.stdlib.LoggerFactory(),
        wrapper_class=structlog.stdlib.BoundLogger,
        cache_logger_on_first_use=True,
    )

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(
        structlog.stdlib.ProcessorFormatter(
            processors=[
                structlog.stdlib.ProcessorFormatter.remove_processors_meta,
                structlog.dev.ConsoleRenderer(),
            ],
            foreign_pre_chain=shared_processors,
        )
    )

    root = logging.getLogger()
    root.setLevel(os.getenv("PYMOPSMAP_LOG_LEVEL", _DEFAULT_LEVEL).upper())
    root.addHandler(console_handler)

    log_file = os.getenv("PYMOPSMAP_LOG_FILE")
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(
            structlog.stdlib.ProcessorFormatter(
                processors=[
                    structlog.stdlib.ProcessorFormatter.remove_processors_meta,
                    structlog.processors.JSONRenderer(),
                ],
                foreign_pre_chain=shared_processors,
            )
        )
        root.addHandler(file_handler)

    logging.getLogger("matplotlib").setLevel(logging.WARNING)

    _INITIALIZED = True


def get_logger(name: str = "pymopsmap") -> structlog.stdlib.BoundLogger:
    _configure()
    return structlog.stdlib.get_logger(name)
