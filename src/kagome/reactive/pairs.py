"""Reaction bookkeeping wrapper around the minimal bias kernel record.

Paper: arXiv:2511.22874, Mori et al. §2.2 (identification set P, reaction
confirmation).

:class:`kagome.boost.tdbb.PairBias` is deliberately reduced to the fields the
bias kernel reads (idx_a/idx_b/is_formation/r0).  The reaction accounting a
candidate needs — which candidate a pair belongs to, whether a confirmed
formation counts toward alpha/Carothers p, and whether the pair gates the
trigger conjunction — is carried here on :class:`TrackedPair`, which *composes*
a ``PairBias`` (``self.bias``) rather than duplicating the kernel logic.

The consumers of this metadata are :class:`kagome.reactive.bonds.BondTracker`
and the selection audit log; the bias kernels only ever see ``tracked.bias``.
"""
from __future__ import annotations

from dataclasses import dataclass

from kagome.boost.tdbb import PairBias


@dataclass
class TrackedPair:
    """A biased pair plus its reaction-accounting metadata.

    ``bias`` is the minimal :class:`PairBias` the kernel consumes.  The
    ``idx_a`` / ``idx_b`` / ``is_formation`` / ``r0`` properties delegate to it
    so BondTracker (and the audit log) can read a single flat object without
    reaching through ``.bias`` — keeping their bodies unchanged from when
    ``PairBias`` still carried these fields.

    Metadata fields:

    - ``candidate_id``: which selected candidate this pair belongs to (restarts
      each cycle). ``< 0`` marks activation / legacy / pre-candidate_id paths,
      which BondTracker confirms per-pair instead of per-candidate.
    - ``counts_as_reaction``: if ``False`` a confirmed formation is excluded
      when counting reactions (alpha / Carothers p); bias and topology effects
      are unchanged. Used for nylon water-forming k-l pairs
      (specs/decisions.md 2026-07-06 A5).
    - ``is_trigger``: whether this pair belongs to the paper's identification
      set P (§2.2 step 3-4). The reaction event fires / confirms only when ALL
      trigger pairs of a candidate simultaneously satisfy their bonding
      condition. Set from ``PairSpec.score_pair`` — the bias-only nylon k-l
      water pair (``score_pair=False``) is ``is_trigger=False`` and never
      participates in the trigger conjunction.
    """
    bias: PairBias
    candidate_id: int = -1
    counts_as_reaction: bool = True
    is_trigger: bool = True

    @property
    def idx_a(self) -> int:
        return self.bias.idx_a

    @property
    def idx_b(self) -> int:
        return self.bias.idx_b

    @property
    def is_formation(self) -> bool:
        return self.bias.is_formation

    @property
    def r0(self) -> float:
        return self.bias.r0
