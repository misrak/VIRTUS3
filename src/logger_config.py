#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Logging configuration for VIRTUS3.

This module provides a configured logger that ensures real-time output,
which is critical when running jobs through schedulers like SLURM where
print statements may be buffered and not appear until the job completes.
"""

import logging
import sys


class ImmediateStreamHandler(logging.StreamHandler):
    """
    StreamHandler that flushes after every emit to ensure immediate output.
    This is crucial for real-time logging in batch job environments.
    """
    def emit(self, record):
        super().emit(record)
        self.flush()


def setup_logger(name, level=logging.INFO, log_file=None):
    """
    Configure and return a logger with immediate console output.

    Parameters
    ----------
    name : str
        Logger name (typically __name__ from calling module)
    level : int
        Logging level (default: logging.INFO)
    log_file : str, optional
        Path to log file. If provided, logs will also be written to this file.

    Returns
    -------
    logging.Logger
        Configured logger instance

    Examples
    --------
    >>> logger = setup_logger(__name__)
    >>> logger.info("Processing started")
    >>> logger.warning("Unmapped BAM file is empty")
    >>> logger.error("Command failed with return code 1")
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Remove any existing handlers to avoid duplicates
    logger.handlers = []

    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Console handler with immediate flush
    console_handler = ImmediateStreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # File handler (optional)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    # Prevent propagation to root logger
    logger.propagate = False

    return logger


def get_logger(name):
    """
    Get or create a logger with the given name.

    If the logger hasn't been set up yet, it will be configured with
    default settings (INFO level, console output only).

    Parameters
    ----------
    name : str
        Logger name (typically __name__ from calling module)

    Returns
    -------
    logging.Logger
        Logger instance
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        logger = setup_logger(name)
    return logger
