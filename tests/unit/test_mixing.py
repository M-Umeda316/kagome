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
from kagome.prep.mix_md import MixMDConfig, run_mix_md
from kagome.prep.mixing import (
    ClassicalMix,
    FragmentParamCache,
    MixTranslatorConfig,
    build_classical_mix,
    connected_components,
    fragment_cache_key,
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


# ── config validation (review 2026-07-17, findings 1 & 4) ─────────────────────

def test_config_rejects_unknown_charge_method():
    with pytest.raises(ValueError, match='charge_method'):
        MixTranslatorConfig(charge_method='am1bcc')


def test_config_rejects_split_nonbonded_forces():
    """combine_nonbonded_forces=False would leave placeholder H missing from
    the split vdW force (particle-count mismatch) — rejected up front."""
    with pytest.raises(ValueError, match='combine_nonbonded_forces'):
        MixTranslatorConfig(combine_nonbonded_forces=False)


# ── cache key separation (finding 3) ──────────────────────────────────────────

def test_cache_key_separates_charge_configs():
    """Templates charged under different configs must never collide."""
    base = MixTranslatorConfig(charge_method='gasteiger')
    keys = {
        fragment_cache_key('CC', base),
        fragment_cache_key('CC', MixTranslatorConfig(charge_method='nagl')),
        fragment_cache_key('CC', MixTranslatorConfig(
            charge_method='nagl', nagl_model='other-model.pt')),
        fragment_cache_key('CC', MixTranslatorConfig(
            charge_method='gasteiger', forcefield='openff-2.1.0.offxml')),
    }
    assert len(keys) == 4          # all distinct
    # same config + same graph -> same key (cache still works)
    assert fragment_cache_key('CC', base) == fragment_cache_key(
        'CC', MixTranslatorConfig(charge_method='gasteiger'))


# ── multi-cap geometry (finding 2) ────────────────────────────────────────────

def test_multi_cap_positions_are_distinct():
    """A 2-deficit centre (carbene-like CH2) gets two NON-coincident cap H;
    coincident caps would give a zero H-C-H angle and NaN angle forces."""
    from kagome.prep.mixing import _cap_positions

    species = ['C', 'H', 'H']
    topo = BondTopology.from_bonds([(0, 1, 1.0), (0, 2, 1.0)])
    pos = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])

    caps = _cap_positions(0, 2, pos, topo, len(species))
    assert caps.shape == (2, 3)
    # both at the C-H bond length from the parent
    np.testing.assert_allclose(
        np.linalg.norm(caps - pos[0], axis=1), 1.09, atol=1e-9,
    )
    # and clearly separated from each other
    assert np.linalg.norm(caps[0] - caps[1]) > 1.0


def test_isolated_heavy_atom_caps_all_distinct():
    """An atom that lost every bond (4-deficit C) gets 4 pairwise-distinct caps."""
    from kagome.prep.mixing import _cap_positions

    topo = BondTopology()          # no bonds at all
    pos = np.zeros((1, 3))
    caps = _cap_positions(0, 4, pos, topo, 1)
    assert caps.shape == (4, 3)
    for a in range(4):
        for b in range(a + 1, 4):
            assert np.linalg.norm(caps[a] - caps[b]) > 0.5


def test_fragment_positions_spread_caps_through_build_path():
    """The full fragment path (RDKit build + positions) keeps 2 caps distinct."""
    from kagome.prep.mixing import _build_fragment_rdkit, _fragment_positions

    species = ['C', 'H', 'H']
    topo = BondTopology.from_bonds([(0, 1, 1.0), (0, 2, 1.0)])
    pos = np.array([[0.0, 0.0, 0.0], [1.09, 0.0, 0.0], [0.0, 1.09, 0.0]])

    caps = radical_cap_counts(topo, species, [0, 1, 2])
    assert caps == {0: 2}
    mol, local_to_global, cap_parents = _build_fragment_rdkit(
        topo, species, [0, 1, 2], caps,
    )
    assert cap_parents == [0, 0]
    coords = _fragment_positions(local_to_global, cap_parents, pos, topo, 3)
    cap_coords = [c for c, g in zip(coords, local_to_global) if g is None]
    assert len(cap_coords) == 2
    assert np.linalg.norm(cap_coords[0] - cap_coords[1]) > 1.0


# ── unknown-element guard (finding 5) ─────────────────────────────────────────

def test_unknown_element_deficit_raises():
    """A valence-deficient atom of an element outside NEUTRAL_VALENCE must fail
    loudly (silent skipping would corrupt the omm<->mlip map)."""
    species = ['Si', 'H', 'H']     # SiH2: deficit 2, Si not in NEUTRAL_VALENCE
    topo = BondTopology.from_bonds([(0, 1, 1.0), (0, 2, 1.0)])
    with pytest.raises(ValueError, match='NEUTRAL_VALENCE'):
        radical_cap_counts(topo, species, [0, 1, 2])


def test_unknown_element_saturated_passes():
    """A saturated unknown element is allowed through (no cap needed)."""
    species = ['Si', 'H', 'H', 'H', 'H']   # SiH4: saturated
    topo = BondTopology.from_bonds(
        [(0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0), (0, 4, 1.0)],
    )
    caps = radical_cap_counts(topo, species, list(range(5)))
    assert caps == {}


# ── write_back input validation (finding 6) ───────────────────────────────────

def test_write_back_rejects_nonfinite_input():
    """MD-divergence NaN must be reported as bad INPUT, not as a mapping bug."""
    mix = ClassicalMix(
        system=None,
        positions_A=np.zeros((2, 3)),
        box_vectors_A=np.eye(3) * 10.0,
        omm_to_mlip=[0, 1],
        mlip_to_omm={0: 0, 1: 1},
        n_mlip_atoms=2,
    )
    bad = np.zeros((2, 3))
    bad[1, 2] = np.nan
    with pytest.raises(ValueError, match='non-finite'):
        mix.write_back(bad)


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
    # metadata records the ACTUAL charge method and the minimization contract
    assert mix.metadata['charge_method'] == 'gasteiger'
    assert mix.metadata['charge_methods_used'] == ['gasteiger']
    assert mix.metadata['charge_method_requested'] == 'gasteiger'
    assert mix.metadata['nagl_fallback'] is False
    assert mix.metadata['requires_minimization'] is True


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
def test_small_box_cutoff_is_clamped_and_context_builds():
    """A box smaller than 2× the Sage cutoff must not break Context creation.

    Regression for the WM-P4 orb smoke: a dense/small mixing box (edge < 1.8 nm)
    previously failed with 'cutoff distance cannot be greater than half the
    periodic box size'. build_classical_mix now clamps the nonbonded cutoff to
    fit the box, so the Context builds and evaluates a finite energy.
    """
    import openmm
    from openmm import unit as ommunit

    specs = [(_MONOMER_SMILES, 1)]
    # box=12.5 Å -> half-box 0.625 nm < Sage's 0.9 nm default cutoff.
    pos, species, _pmap, topo, cell = _copolymer_system(specs, 1, box=12.5)

    cfg = MixTranslatorConfig(charge_method='gasteiger')
    mix = build_classical_mix(topo, species, pos, cell, cfg)

    nb = next(f for f in mix.system.getForces()
              if isinstance(f, openmm.NonbondedForce))
    cutoff_nm = nb.getCutoffDistance().value_in_unit(ommunit.nanometer)
    assert cutoff_nm <= 0.49 * (12.5 * 0.1) + 1e-9

    integrator = openmm.VerletIntegrator(1.0 * ommunit.femtosecond)
    context = openmm.Context(
        mix.system, integrator, openmm.Platform.getPlatformByName('Reference'),
    )
    context.setPositions((mix.positions_A * 0.1) * ommunit.nanometer)
    energy = context.getState(getEnergy=True).getPotentialEnergy()
    assert np.isfinite(energy.value_in_unit(ommunit.kilojoule_per_mole))


@requires_openmm
def test_element_verification_catches_reordering():
    """_verify_topology_elements: matching order passes, a species mismatch
    raises (D-4-style guard against silent atom scrambling; finding 8)."""
    from openff.toolkit import Molecule, Topology

    from kagome.prep.mixing import _verify_topology_elements

    water = Molecule.from_smiles('O')      # atom order O, H, H
    top = Topology.from_molecules([water])

    _verify_topology_elements(top, [0, 1, 2], ['O', 'H', 'H'])  # no raise
    with pytest.raises(ValueError, match='element mismatch'):
        _verify_topology_elements(top, [0, 1, 2], ['N', 'H', 'H'])
    with pytest.raises(ValueError, match='atom count'):
        _verify_topology_elements(top, [0, 1], ['O', 'H'])
    # a cap slot (None) must be H: position 0 is O -> mismatch
    with pytest.raises(ValueError, match='element mismatch'):
        _verify_topology_elements(top, [None, 1, 2], ['O', 'H', 'H'])


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


# ── WM-P4: classical mixing MD survives post-reaction injected hot contacts ────

def _dense_reacted_mix(box_A=12.0, place_box_A=30.0, seed=7, declash_injected=False):
    """Build a dense, *reacted* ClassicalMix (cap H + placeholder H injected).

    A confirmed vinyl addition (``apply_vinyl_addition``) caps the new growth-end
    radical and sheds a placeholder H that lands ~1 Å from its former parent — the
    exact injected-hot-contact geometry that crashed the WM-P4 campaign's mixing
    phase. The monomers are packed by placing them in a roomy box and then
    isotropically compressing coordinates + cell to ``box_A`` (~0.4 g/mL): the
    random placer cannot reach that density directly (min-separation), and the
    dense box is what makes the injected contacts bite under dynamics.

    ``declash_injected`` defaults to ``False`` so the warm-up regression below
    sees the ORIGINAL injected hot contacts (the geometric de-clash is exactly
    what removes them); the de-clash tests build with it True.
    """
    specs = [(_MONOMER_SMILES, 2), (_METHACRYLATE_SMILES, 2)]
    pos, species, pmap, topo, _cell = _copolymer_system(
        specs, 1, box=place_box_A, seed=seed)
    radical_c = _radical_c_of_initiator(species, topo, 1, _INITIATOR_SMILES)
    apply_vinyl_addition(topo, radical_c, min(pmap), pmap, species)

    pos = (pos % place_box_A) * (box_A / place_box_A)   # isotropic compression
    cell = np.diag([box_A, box_A, box_A]).astype(float)
    mix = build_classical_mix(
        topo, species, pos, cell,
        MixTranslatorConfig(
            charge_method='gasteiger', declash_injected=declash_injected))
    # sanity: this fixture must actually exercise the injected-atom path
    assert mix.metadata['n_cap_h'] == 1
    assert mix.metadata['n_placeholder_h'] == 1
    return mix


@requires_openmm
def test_mixing_md_survives_injected_hot_contacts_after_reaction():
    """Regression for the WM-P4 mixing NaN (specs/decisions.md 追補 2026-07-18).

    A freshly translated post-reaction system carries injected hot contacts (a
    placeholder H overlapping its former parent, a rigid cap C-H). At production
    scale (~500+ atoms, 0.5 g/mL) the ``LocalEnergyMinimizer`` meets its
    RMS-force tolerance globally while such a contact stays hot — its force is
    diluted across all atoms — so stepping 0.5 fs Langevin straight away diverges
    to ``Particle coordinate is NaN`` inside ``run_mix_md``.

    We reproduce that *un-cleared hot start* deterministically at unit scale by
    running the minimizer with a deliberately loose tolerance (a stand-in for the
    scale-dependent RMS dilution — at unit scale a 10 kJ/mol/nm minimize would
    clear the contact that production could not). Single-precision 'CPU' platform,
    matching the campaign's single-precision GPU arithmetic where the overflow
    actually turns into NaN.

    Without the soft-start warm-up (``warmup_steps=0`` — the pre-fix code path)
    the run diverges; with it (default) the run completes with finite positions.
    """
    import openmm

    mix = _dense_reacted_mix()
    common = dict(
        temperature_K=333.0, n_steps=400, timestep_fs=0.5, platform='CPU',
        minimize_tolerance_kj_mol_nm=1.0e6,   # emulate an un-cleared hot start
    )

    # (1) pre-fix behaviour: no warm-up -> the injected hot contact makes 0.5 fs
    # Langevin diverge to NaN (the exact production failure).
    with pytest.raises((openmm.OpenMMException, RuntimeError)):
        run_mix_md(mix, MixMDConfig(warmup_steps=0, **common), seed=12345)

    # (2) with the soft-start warm-up (default) the same system survives: finite
    # positions, a healthy (finite) mixing energy, and the warm-up counted apart
    # from the reported mixing steps.
    result = run_mix_md(mix, MixMDConfig(**common), seed=12345)
    assert np.isfinite(result.positions_A).all()
    assert np.isfinite(result.final_energy_kj_mol)
    assert result.n_steps == 400                 # reported mixing steps unchanged
    assert result.n_warmup_steps > 0             # warm-up ran, tracked separately
    # write-back maps cleanly (no NaN leaked into the MLIP coordinates)
    restored = mix.write_back(result.positions_A)
    assert np.isfinite(restored).all()


# ── WM-P4 Layer 1: geometric de-clash of injected atoms ───────────────────────

def test_cap_direction_maximises_clearance_away_from_a_neighbour():
    """Cap-H best-direction (pure geometry, no MD): with a real atom sitting on
    the naive vacant-valence axis, the max-clearance placement points the cap the
    OTHER way, so it never lands on the blocking neighbour."""
    from kagome.prep.mixing import _cap_positions, _vacant_axis

    # radical C at origin bonded to one neighbour at +x -> vacant axis is -x.
    species = ['C', 'C']
    topo = BondTopology.from_bonds([(0, 1, 1.0)])
    pos = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
    axis = _vacant_axis(0, pos, topo, 2)
    np.testing.assert_allclose(axis, [-1.0, 0.0, 0.0], atol=1e-9)

    # Put a blocker exactly where the naive cap (p + 1.09*axis) would go.
    blocker = pos[0] + 1.09 * axis
    pos2 = np.vstack([pos, blocker[None, :]])            # atoms: C, C, blocker
    topo2 = BondTopology.from_bonds([(0, 1, 1.0)])       # blocker unbonded
    cap = _cap_positions(0, 1, pos2, topo2, 3)[0]

    # naive placement would collide with the blocker; the de-clash placement
    # keeps well clear of it while preserving the C-H bond length.
    assert np.linalg.norm(cap - blocker) > 1.0
    np.testing.assert_allclose(np.linalg.norm(cap - pos2[0]), 1.09, atol=1e-9)


def test_legacy_cap_placement_ignores_clearance():
    """maximize_clearance=False reproduces the legacy 'opposite the bond' cap."""
    from kagome.prep.mixing import _cap_positions

    species = ['C', 'C']
    topo = BondTopology.from_bonds([(0, 1, 1.0)])
    pos = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
    cap = _cap_positions(0, 1, pos, topo, 2, maximize_clearance=False)[0]
    np.testing.assert_allclose(cap, [-1.09, 0.0, 0.0], atol=1e-9)


def test_declash_moves_only_the_injected_atom_out_of_a_hard_overlap():
    """_declash_injected relocates a coincident placeholder off its victim to a
    safe clearance, and leaves every original atom exactly where it was."""
    from kagome.prep.mixing import _box_diag, _declash_injected, _HARD_CLASH_A

    # 4 real atoms + one placeholder (index 4) placed on top of atom 0.
    pos = np.array([
        [0.0, 0.0, 0.0],
        [5.0, 0.0, 0.0],
        [0.0, 5.0, 0.0],
        [0.0, 0.0, 5.0],
        [0.0, 0.0, 0.0],        # placeholder coincident with atom 0
    ])
    box = np.diag([20.0, 20.0, 20.0]).astype(float)
    out, moved = _declash_injected(
        pos, cap_parent_omm={}, placeholder_omm=[4], box_diag=_box_diag(box))
    assert moved == 1
    # originals untouched
    np.testing.assert_array_equal(out[:4], pos[:4])
    # placeholder pushed to a real clearance from every other atom
    dists = np.linalg.norm(out[:4] - out[4], axis=1)
    assert dists.min() >= _HARD_CLASH_A


def test_declash_noop_when_no_hard_overlap():
    """An injected atom already clear of everything is left in place."""
    from kagome.prep.mixing import _box_diag, _declash_injected

    pos = np.array([
        [0.0, 0.0, 0.0],
        [5.0, 0.0, 0.0],
        [2.5, 2.5, 2.5],        # placeholder, far from both reals
    ])
    box = np.diag([20.0, 20.0, 20.0]).astype(float)
    out, moved = _declash_injected(
        pos, cap_parent_omm={}, placeholder_omm=[2], box_diag=_box_diag(box))
    assert moved == 0
    np.testing.assert_array_equal(out, pos)


@requires_openmm
def test_declash_rescues_hard_injected_overlap_end_to_end():
    """Layer 1 rescues a hard, near-singular injected overlap under real MD.

    Build a dense *reacted* copolymer (cap H + placeholder H injected) and force
    the shed placeholder H to COINCIDE with a real atom — a true LJ singularity
    that minimize + warm-up cannot escape. Only the build-time ``declash_injected``
    flag differs between the two runs; the MD config (warm-up ON) is identical, so
    the geometric de-clash is the sole cause of survival:

    * OFF -> the coincident injected atom diverges (NaN / RuntimeError);
    * ON  -> the placeholder is relocated before minimization and the run
      completes with finite positions that write back cleanly.
    """
    import openmm

    specs = [(_MONOMER_SMILES, 2), (_METHACRYLATE_SMILES, 2)]
    place_box, box_A = 30.0, 12.0
    pos, species, pmap, topo, _cell = _copolymer_system(
        specs, 1, box=place_box, seed=7)
    radical_c = _radical_c_of_initiator(species, topo, 1, _INITIATOR_SMILES)
    apply_vinyl_addition(topo, radical_c, min(pmap), pmap, species)
    pos = (pos % place_box) * (box_A / place_box)        # isotropic compression
    cell = np.diag([box_A, box_A, box_A]).astype(float)

    placeholder = [a for a in range(len(species))
                   if is_placeholder_h(topo, species, a)]
    assert len(placeholder) == 1
    g = placeholder[0]
    victim = next(a for a in range(len(species))
                  if species[a] != 'H' and a != g)
    pos[g] = pos[victim].copy()                          # force a hard overlap

    # Deliberately WEAK minimization: at unit scale the minimizer is otherwise
    # strong enough to relieve the overlap that production's 564-atom RMS-force
    # dilution could not — so we stand that dilution in with a loose tolerance
    # (same device as the WM-P4 warm-up regression) AND disable the tighter
    # second-pass minimize (warmup tol == main tol, so it is skipped) so the hard
    # overlap survives into the dynamics exactly as it did in production. The
    # soft-start warm-up stays ON — identical for both runs — so the ONLY
    # difference is the build-time de-clash.
    common = dict(
        temperature_K=333.0, n_steps=200, timestep_fs=0.5, platform='CPU',
        minimize_tolerance_kj_mol_nm=1.0e6,
        warmup_minimize_tolerance_kj_mol_nm=1.0e6,
    )

    # Layer 1 OFF: coincident injected atom -> singular -> diverges.
    mix_off = build_classical_mix(
        topo, species, pos, cell,
        MixTranslatorConfig(charge_method='gasteiger', declash_injected=False))
    assert mix_off.metadata['n_declashed'] == 0
    with pytest.raises((openmm.OpenMMException, RuntimeError)):
        run_mix_md(mix_off, MixMDConfig(**common), seed=1)

    # Layer 1 ON (default): de-clash relocates the injected atom -> finite run.
    mix_on = build_classical_mix(
        topo, species, pos, cell,
        MixTranslatorConfig(charge_method='gasteiger', declash_injected=True))
    assert mix_on.metadata['n_declashed'] >= 1
    result = run_mix_md(mix_on, MixMDConfig(**common), seed=1)
    assert np.isfinite(result.positions_A).all()
    assert np.isfinite(result.final_energy_kj_mol)
    restored = mix_on.write_back(result.positions_A)
    assert np.isfinite(restored).all()


# ── WM-P4 bridge: ungated warm-up survives orb->Sage handed-over close contacts ─

def _two_molecule_clash(gap_A, box=24.0, seed=5):
    """Two free MA monomers (UNREACTED: no injected atoms), with molecule B
    rigidly translated so its nearest atom to molecule A sits ``gap_A`` apart.

    Rigid translation keeps every intramolecular bond length — hence every Sage
    X-H constraint — intact, so this reproduces the production failure faithfully:
    orb hands over a geometry with a pair closer than the classical Sage LJ wall
    tolerates (orb's learned repulsion is softer), and the clashing atoms include
    Sage-constrained H that a geometric push-apart cannot separate. Returns
    ``(topology, species, positions_A, cell_A)``.
    """
    specs = [(_MONOMER_SMILES, 2)]
    pos, species, _pmap, topo, cell = _copolymer_system(specs, 0, box=box, seed=seed)
    comps = connected_components(topo, len(species))
    assert len(comps) == 2, 'fixture expects two disjoint monomers'
    mol_a, mol_b = comps[0], comps[1]

    pa = pos[mol_a]
    pb = pos[mol_b]
    # nearest inter-molecular pair (roomy box: ignore PBC for the setup)
    diff = pa[:, None, :] - pb[None, :, :]
    d = np.linalg.norm(diff, axis=2)
    ia, ib = np.unravel_index(int(np.argmin(d)), d.shape)
    a_glob, b_glob = mol_a[ia], mol_b[ib]
    d0 = float(d[ia, ib])
    u = (pos[b_glob] - pos[a_glob]) / d0
    # translate ALL of molecule B (rigid) so the nearest pair sits gap_A apart
    delta = (gap_A - d0) * u
    pos = pos.copy()
    for g in mol_b:
        pos[g] = pos[g] + delta
    return topo, species, pos, cell


@requires_openmm
def test_ungated_warmup_survives_noninjected_handed_over_contact():
    """Reproduces production cycle-2: a NON-injected mixing cycle whose handed-over
    coordinates contain a Sage-intolerable close contact.

    Pre-fix the warm-up + tighter 2nd-minimize were gated on ``has_injected`` so
    for a cycle with no cap/placeholder H (cycle 2 selected but confirmed nothing)
    NEITHER ran — only the loose first minimize, then 0.5 fs, which diverged. With
    the ungated bridge (default) the same system survives. The loose minimize
    tolerance stands in for production's 564-atom RMS-force dilution (the device
    used by the WM-P4 warm-up regression) so the contact reaches the dynamics.
    """
    import openmm

    topo, species, pos, cell = _two_molecule_clash(gap_A=0.9)
    mix = build_classical_mix(
        topo, species, pos, cell, MixTranslatorConfig(charge_method='gasteiger'))
    assert mix.metadata['n_cap_h'] == 0
    assert mix.metadata['n_placeholder_h'] == 0        # genuinely non-injected
    common = dict(
        temperature_K=333.0, n_steps=200, timestep_fs=0.5, platform='CPU',
        minimize_tolerance_kj_mol_nm=1.0e12,
    )

    # pre-fix behaviour (warm-up disabled, as the has_injected gate did here)
    with pytest.raises((openmm.OpenMMException, RuntimeError)):
        run_mix_md(mix, MixMDConfig(warmup_steps=0, **common), seed=3)

    # ungated bridge (default warm-up now runs every cycle) -> survives
    result = run_mix_md(mix, MixMDConfig(**common), seed=3)
    assert np.isfinite(result.positions_A).all()
    assert np.isfinite(result.final_energy_kj_mol)
    restored = mix.write_back(result.positions_A)
    assert np.isfinite(restored).all()


@requires_openmm
def test_ungated_bridge_survives_constrained_hh_contact():
    """The soft-core DECIDER: two molecules pushed to ~0.7 A closest approach with
    Sage's rigid X-H constraints intact — the case a geometric push-apart cannot
    fix (a displaced constrained H snaps back), so only bounded overdamped
    dynamics moving whole rigid groups can separate them.

    If the ungated tighter-minimize + soft-start warm-up survive this, the
    orb<->Sage wall mismatch is bridged WITHOUT a soft-core pre-relaxation and the
    graceful skip stays the only remaining net. If it diverges, a soft-core /
    force-capped pre-relaxation is required.

    EVIDENCE (this test, 2026-07-18): the bridge SURVIVES the constrained-H clash,
    so no soft-core was implemented — see specs/decisions.md 追補 2026-07-18
    (WM-P4 bridge).
    """
    import openmm

    topo, species, pos, cell = _two_molecule_clash(gap_A=0.7)
    mix = build_classical_mix(
        topo, species, pos, cell, MixTranslatorConfig(charge_method='gasteiger'))
    assert mix.metadata['n_cap_h'] == 0
    assert mix.metadata['n_placeholder_h'] == 0
    common = dict(
        temperature_K=333.0, n_steps=200, timestep_fs=0.5, platform='CPU',
        minimize_tolerance_kj_mol_nm=1.0e12,
    )
    # Without the bridge (warm-up + tighter minimize disabled) the constrained-H
    # clash is genuinely unrecoverable — proves the contact is hard and the
    # 'survives' assertion below is not vacuous.
    with pytest.raises((openmm.OpenMMException, RuntimeError)):
        run_mix_md(mix, MixMDConfig(warmup_steps=0, **common), seed=3)

    # The ungated bridge (default) rescues it: bounded overdamped dynamics move
    # the whole rigid molecules apart where a geometric push-apart could not.
    result = run_mix_md(mix, MixMDConfig(**common), seed=3)
    assert np.isfinite(result.positions_A).all()
    assert np.isfinite(result.final_energy_kj_mol)
    restored = mix.write_back(result.positions_A)
    assert np.isfinite(restored).all()
