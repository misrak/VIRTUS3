#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Entry point for module execution: python -m virtus3

This allows the package to be invoked as:
    python -m virtus3 [args...]
"""

import sys
from .virtus3 import main

if __name__ == "__main__":
    sys.exit(main())
