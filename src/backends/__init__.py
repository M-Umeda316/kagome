"""MLIP/calculator backends.

Re-exports the backend-agnostic ``Calculator`` interface and the dependency-free
``ToyCalculator`` for a stable package-level API surface. Heavy backends
(MACE/OrbMol-v2) are intentionally NOT re-exported here so that importing
``src.backends`` never pulls optional ML dependencies; import their factory
functions directly from ``src.backends.mace_backend`` / ``src.backends.orb_backend``.
"""
from src.backends.base import Calculator
from src.backends.toy import ToyCalculator

__all__ = ['Calculator', 'ToyCalculator']
