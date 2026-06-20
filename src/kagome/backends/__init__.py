"""MLIP/calculator backends.

Re-exports the backend-agnostic ``Calculator`` interface and the dependency-free
``ToyCalculator`` for a stable package-level API surface. Heavy backends
(MACE/OrbMol-v2) are intentionally NOT re-exported here so that importing
``kagome.backends`` never pulls optional ML dependencies; import their factory
functions directly from ``kagome.backends.mace_backend`` / ``kagome.backends.orb_backend``.
"""
from kagome.backends.base import Calculator
from kagome.backends.toy import ToyCalculator

__all__ = ['Calculator', 'ToyCalculator']
