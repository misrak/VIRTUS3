"""
VIRTUS3: Detection of viral transcripts in single-cell RNA-seq data

This package provides a two-step pipeline combining CellRanger and Alevin
for identifying and quantifying viral transcripts in single-cell RNA sequencing data.
"""

from ._version import __version__
from .virtus3 import main, pipeline, run_command, analyze_fastq_name

__all__ = [
    '__version__',
    'main',
    'pipeline',
    'run_command',
    'analyze_fastq_name',
]

# Entry point function for console script
def cli_main():
    """Console script entry point for the virtus3 command."""
    import sys
    sys.exit(main())
