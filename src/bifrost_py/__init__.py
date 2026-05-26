"""
Docstring needed
"""

__version__ = "1.0.0"

"""Bifrost: fiber birefringence simulation library."""
from .bifrost_py import start, info

# Import module-level magic functions
from .bifrost_py import __getattr__, __dir__

__all__ = ["start", "info"]