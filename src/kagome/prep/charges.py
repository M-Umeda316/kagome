"""Shared partial-charge assignment for the classical prep/mixing stages.

Single implementation of the "NAGL first, Gasteiger fallback" policy (decision
D-2, specs/decisions.md 2026-06-14) used by both the initial-structure prep
(:mod:`kagome.prep.openmm_equilibrate`) and the well-mixed mixing translator
(:mod:`kagome.prep.mixing`, 2026-07-17). Extracted so the two callers cannot
drift (they previously duplicated the NAGL try/except with slightly different
availability checks and Gasteiger implementations).

All OpenFF imports are deferred to call time so this module imports cleanly in
ML-only environments.
"""
from __future__ import annotations

import logging

import numpy as np

logger = logging.getLogger(__name__)

CHARGE_METHODS = ('nagl', 'gasteiger')


def assign_charges(offmol, charge_method: str, nagl_model: str) -> str:
    """Assign partial charges to an OpenFF ``Molecule`` in place.

    ``charge_method`` is ``'nagl'`` or ``'gasteiger'`` (anything else raises
    ``ValueError``). The NAGL path falls back to Gasteiger on any failure
    (missing package, model download unavailable, unsupported chemistry) with a
    warning — prep/mixing must not hard-fail on charge niceties.

    Returns the method actually used (``'nagl'`` or ``'gasteiger'``) so callers
    can record a fallback truthfully in their metadata/cache keys.
    """
    if charge_method not in CHARGE_METHODS:
        raise ValueError(
            f'charge_method must be one of {CHARGE_METHODS}; '
            f'got {charge_method!r}'
        )
    if charge_method == 'nagl':
        try:
            from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper

            offmol.assign_partial_charges(
                nagl_model, toolkit_registry=NAGLToolkitWrapper(),
            )
            return 'nagl'
        except Exception as exc:  # noqa: BLE001 - fall back, don't fail prep
            logger.warning(
                'NAGL charge assignment failed (%s); falling back to Gasteiger.',
                exc,
            )
    _assign_gasteiger(offmol)
    return 'gasteiger'


def _assign_gasteiger(offmol) -> None:
    """Gasteiger charges via the OpenFF/RDKit wrapper, with a NaN scrub.

    Gasteiger can emit NaN for exotic valences; non-finite charges are zeroed so
    the downstream OpenMM system stays finite (warned, not fatal).
    """
    from openff.units import unit as offunit

    offmol.assign_partial_charges('gasteiger')
    charges = np.array(
        [c.m_as(offunit.elementary_charge) for c in offmol.partial_charges],
        dtype=np.float64,
    )
    if not np.isfinite(charges).all():
        logger.warning(
            'Gasteiger produced non-finite charges; zeroing the bad entries.',
        )
        charges = np.nan_to_num(charges, nan=0.0, posinf=0.0, neginf=0.0)
        offmol.partial_charges = charges * offunit.elementary_charge
