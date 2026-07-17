"""Unit tests for the bond-graph -> OpenFF/OpenMM mixing translator (WM-P2).

Covers, MD-free where possible (specs/decisions.md 2026-07-17 "well-mixed 測定
モード", P2):

* graph decomposition into molecules (connected components);
* radical cap-H counting (opened vinyl growth end vs. sp2 / carbonyl carbons);
* placeholder-H detection (shed initiator H, disconnected from the graph);
* cap-H injection / removal round trip (atom count, order, bond graph preserved);
* fragment parameterization cache hits on isomorphic fragments;
* a small OpenMM System build + single energy evaluation (no MD).

The OpenMM/OpenFF-dependent tests are skipped when those packages are absent so
the pure-graph tests still run in an ML-only environment. Charge assignment uses
Gasteiger (offline, deterministic) so tests do not depend on a NAGL download.
"""
from __future__ import annotations

import numpy as np
import pytest

from kagome.reactive.topology import BondTopology, apply_vinyl_addition
from kagome.prep.mixing import (
    ClassicalMix,
    FragmentParamCache,
    MixTranslatorConfig,
    build_classical_mix,
    connected_components,
    is_placeholder_h,
    radical_cap_counts,
)

from scripts._systems import (
    _INITIATOR_SMILES,
    _METHACRYLATE_SMILES,
    _MONOMER_SMILES,
    build_vinyl_copolymer_system,
    copolymer_initial_bonds,
)

def _has_openmm() -> bool:
    try:
        import openmm  # noqa: F401
        import openff.toolkit  # noqa: F401
        return True
    except Exception:
        return False


requires_openmm = pytest.mark.skipif(
    not _has_openmm(), reason='OpenMM/OpenFF not installed',
)


# ── small realistic systems ───────────────────────────────────────────────────

def _copolymer_system(monomer_specs, n_initiators, box=22.0, seed=7):
    """Build a small copolymer box and its initial bond topology.

    Returns (positions_A, species, propagation_map, topology, cell_A).
    """
    rng = np.random.default_rng(seed)
    positions, species, _template, _groups, pmap, _chain = (
        build_vinyl_copolymer_system(monomer_specs, n_initiators, box, rng)
    )
    bonds = copolymer_initial_bonds(monomer_specs, n_initiators)
    topo = BondTopology.from_bonds(bonds)
    cell = np.diag([box, box, box]).astype(float)
    return positions, species, pmap, topo, cell


def _radical_c_of_initiator(species, topology, n_initiators, initiator_smiles):
    """Global index of the (closed-shell) initiator radical carbon: the C with
    3 carbon neighbours in the first initiator block (offset 0)."""
    from scripts._systems import _find_ibn_radical_c

    return _find_ibn_radical_c(initiator_smiles)  # first initiator at offset 0


# ── graph decomposition ────────────────────────────────────────────────────────

def test_connected_components_unreacted():
    """One initiator + one MA + one MMA -> three disjoint molecules."""
    specs = [(_MONOMER_SMILES, 1), (_METHACRYLATE_SMILES, 1)]
    _pos, species, _pmap, topo, _cell = _copolymer_system(specs, 1)
    comps = connected_components(topo, len(species))
    assert len(comps) == 3
    # components are sorted by min index and internally sorted
    assert comps == sorted(comps, key=min)
    assert sum(len(c) for c in comps) == len(species)
    # atoms partition exactly
    assert sorted(a for c in comps for a in c) == list(range(len(species)))


def test_reaction_merges_components_and_makes_radical():
    """A confirmed addition joins initiator+MA, sheds a placeholder H, and opens
    the vinyl so the MA beta carbon becomes a radical needing one cap."""
    specs = [(_MONOMER_SMILES, 1), (_METHACRYLATE_SMILES, 1)]
    _pos, species, pmap, topo, _cell = _copolymer_system(specs, 1)

    radical_c = _radical_c_of_initiator(species, topo, 1, _INITIATOR_SMILES)
    # alpha carbon of the MA monomer = the vinyl_alpha_C key in pmap sitting in
    # the MA block (first monomer after the initiator).
    alpha = min(pmap)  # smallest alpha index = the MA monomer's alpha
    apply_vinyl_addition(topo, radical_c, alpha, pmap, species)

    comps = connected_components(topo, len(species))
    # initiator+MA fused (one component contains both radical_c and alpha)
    fused = [c for c in comps if radical_c in c][0]
    assert alpha in fused

    # exactly one shed placeholder H (lone, coordination 0)
    lone = [a for a in range(len(species)) if is_placeholder_h(topo, species, a)]
    assert len(lone) == 1
    assert species[lone[0]] == 'H'

    # the MA beta carbon (propagation target) is now a 1-cap radical
    beta = pmap[alpha]
    caps = radical_cap_counts(topo, species, fused)
    assert caps.get(beta) == 1


def test_radical_cap_counts_ignores_sp2_and_carbonyl():
    """Intact MA: no atom is a radical (double bonds keep order sums saturated)."""
    specs = [(_MONOMER_SMILES, 1)]
    _pos, species, _pmap, topo, _cell = _copolymer_system(specs, 0)
    caps = radical_cap_counts(topo, species, list(range(len(species))))
    assert caps == {}


def test_placeholder_h_only_for_disconnected_h():
    specs = [(_MONOMER_SMILES, 1)]
    _pos, species, _pmap, topo, _cell = _copolymer_system(specs, 0)
    # every H is bonded in an intact monomer -> no placeholders
    assert not any(
        is_placeholder_h(topo, species, a) for a in range(len(species))
    )


# ── OpenMM-dependent: build, round trip, cache, energy ─────────────────────────

@requires_openmm
def test_build_smoke_and_index_map_complete():
    specs = [(_MONOMER_SMILES, 1), (_METHACRYLATE_SMILES, 1)]
    pos, species, _pmap, topo, cell = _copolymer_system(specs, 1)
    cfg = MixTranslatorConfig(charge_method='gasteiger')
    mix = build_classical_mix(topo, species, pos, cell, cfg)

    assert isinstance(mix, ClassicalMix)
    # unreacted: no caps, no placeholders -> OMM particle count == N
    assert mix.metadata['n_cap_h'] == 0
    assert mix.metadata['n_placeholder_h'] == 0
    assert mix.system.getNumParticles() == len(species)
    # index map covers every MLIP atom exactly once
    assert len(mix.mlip_to_omm) == len(species)
    assert set(mix.mlip_to_omm) == set(range(len(species)))


@requires_openmm
def test_write_back_roundtrip_preserves_positions_and_graph():
    """After capping a radical and shedding a placeholder H, write-back restores
    every MLIP atom's coordinate and preserves the bond graph (round-trip)."""
    specs = [(_MONOMER_SMILES, 1), (_METHACRYLATE_SMILES, 1)]
    pos, species, pmap, topo, cell = _copolymer_system(specs, 1)

    radical_c = _radical_c_of_initiator(species, topo, 1, _INITIATOR_SMILES)
    alpha = min(pmap)
    apply_vinyl_addition(topo, radical_c, alpha, pmap, species)
    bonds_before = topo.bonds()

    cfg = MixTranslatorConfig(charge_method='gasteiger')
    mix = build_classical_mix(topo, species, pos, cell, cfg)

    assert mix.metadata['n_cap_h'] == 1        # MA beta radical capped
    assert mix.metadata['n_placeholder_h'] == 1  # shed initiator H
    # OMM particles = real + caps (in fragments) + placeholder (appended)
    assert mix.system.getNumParticles() == len(species) + 1

    # write-back of the (unchanged) starting coords reproduces MLIP positions
    restored = mix.write_back(mix.positions_A)
    assert restored.shape == (len(species), 3)
    np.testing.assert_allclose(restored, pos, atol=1e-9)

    # the caller's bond graph is untouched by translation (round-trip)
    assert topo.bonds() == bonds_before

    # cap H are the only omm particles mapping to None; exactly one of them
    assert sum(1 for m in mix.omm_to_mlip if m is None) == 1


@requires_openmm
def test_fragment_cache_hits_on_isomorphic_monomers():
    """Two identical free MA monomers share a graph key: charge assignment runs
    once per distinct graph (initiator + MA = 2 misses, second MA = 1 hit)."""
    specs = [(_MONOMER_SMILES, 2)]
    pos, species, _pmap, topo, cell = _copolymer_system(specs, 1)
    cache = FragmentParamCache()
    cfg = MixTranslatorConfig(charge_method='gasteiger')
    build_classical_mix(topo, species, pos, cell, cfg, cache=cache)

    assert cache.misses == 2   # initiator graph + MA graph
    assert cache.hits == 1     # the second MA reused the MA template


@requires_openmm
def test_energy_evaluation_is_finite():
    """A single energy evaluation on the built system is finite (no MD)."""
    import openmm
    from openmm import unit as ommunit

    specs = [(_MONOMER_SMILES, 1), (_METHACRYLATE_SMILES, 1)]
    pos, species, pmap, topo, cell = _copolymer_system(specs, 1)
    radical_c = _radical_c_of_initiator(species, topo, 1, _INITIATOR_SMILES)
    alpha = min(pmap)
    apply_vinyl_addition(topo, radical_c, alpha, pmap, species)

    cfg = MixTranslatorConfig(charge_method='gasteiger')
    mix = build_classical_mix(topo, species, pos, cell, cfg)

    integrator = openmm.VerletIntegrator(1.0 * ommunit.femtosecond)
    context = openmm.Context(
        mix.system, integrator, openmm.Platform.getPlatformByName('Reference'),
    )
    context.setPositions((mix.positions_A * 0.1) * ommunit.nanometer)  # Å -> nm
    energy = context.getState(getEnergy=True).getPotentialEnergy()
    assert np.isfinite(energy.value_in_unit(ommunit.kilojoule_per_mole))


@requires_openmm
def test_cache_reuse_across_two_calls():
    """Passing the same cache to a second call reuses templates (all hits)."""
    specs = [(_MONOMER_SMILES, 1)]
    pos, species, _pmap, topo, cell = _copolymer_system(specs, 1)
    cache = FragmentParamCache()
    cfg = MixTranslatorConfig(charge_method='gasteiger')

    build_classical_mix(topo, species, pos, cell, cfg, cache=cache)
    misses_first = cache.misses
    hits_before = cache.hits
    build_classical_mix(topo, species, pos, cell, cfg, cache=cache)
    # second call adds no new misses, only hits for each fragment
    assert cache.misses == misses_first
    assert cache.hits > hits_before
