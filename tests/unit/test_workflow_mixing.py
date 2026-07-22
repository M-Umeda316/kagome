"""Unit tests for the WM-P3 classical mixing-phase workflow integration.

Covers (specs/decisions.md 2026-07-17 "well-mixed 測定モード" (i), 追補
2026-07-18), MD-free where possible:

* MixConfig / MixMDConfig validation and n_mix_steps rounding;
* run() entry gating: mixing requires bond-topology tracking and a periodic
  cell, and the error fires before any MD;
* dataclasses.asdict round-trip of the nested mixing config (manifest
  provenance relies on it);
* one-cycle integration: mixing.jsonl audit record, mix_settle trajectory
  frames, CycleLog phase ordering, finite state, manifest extra;
* seed determinism of the full mixing pipeline (Reference platform);
* settle_steps=0 degenerate case (mixing only, no settle segment).

The OpenMM/OpenFF-dependent tests are skipped when those packages are absent
so the pure tests still run in an ML-only environment (this module never
imports openmm at module level). Charge assignment uses Gasteiger (offline,
deterministic) so tests do not depend on a NAGL download.
"""
from __future__ import annotations

import json
from dataclasses import asdict

import numpy as np
import pytest

from kagome.backends.toy import ToyCalculator
from kagome.prep.mix_md import MixMDConfig
from kagome.reactive.groups import PairSpec, ReactiveGroup, ReactionTemplate
from kagome.workflows.polymerization import (
    MixConfig,
    PolymerizationConfig,
    PolymerizationWorkflow,
    SimulationState,
    masses_from_species,
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


# ── MixConfig validation (pure) ───────────────────────────────────────────────

class TestMixConfigValidation:

    def test_mix_time_ps_nonpositive_rejected(self):
        with pytest.raises(ValueError, match='mix_time_ps'):
            MixConfig(mix_time_ps=0.0, settle_steps=0, temperature_K=300.0)
        with pytest.raises(ValueError, match='mix_time_ps'):
            MixConfig(mix_time_ps=-1.0, settle_steps=0, temperature_K=300.0)

    def test_negative_settle_steps_rejected(self):
        with pytest.raises(ValueError, match='settle_steps'):
            MixConfig(mix_time_ps=1.0, settle_steps=-1, temperature_K=300.0)

    def test_zero_settle_steps_allowed(self):
        """settle_steps=0 is a valid degenerate case (mixing only)."""
        cfg = MixConfig(mix_time_ps=1.0, settle_steps=0, temperature_K=300.0)
        assert cfg.settle_steps == 0

    def test_temperature_nonpositive_rejected(self):
        with pytest.raises(ValueError, match='temperature_K'):
            MixConfig(mix_time_ps=1.0, settle_steps=0, temperature_K=0.0)
        with pytest.raises(ValueError, match='temperature_K'):
            MixConfig(mix_time_ps=1.0, settle_steps=0, temperature_K=-10.0)

    def test_n_mix_steps_rounding(self):
        """0.05 ps at 0.5 fs/step -> exactly 100 classical steps."""
        cfg = MixConfig(mix_time_ps=0.05, settle_steps=0, temperature_K=300.0)
        assert cfg.n_mix_steps == 100

    def test_n_mix_steps_clamps_to_one(self):
        """A mixing time shorter than one timestep still runs >= 1 step."""
        cfg = MixConfig(mix_time_ps=1e-6, settle_steps=0, temperature_K=300.0)
        assert cfg.n_mix_steps == 1

    # fail-fast validation of the remaining knobs (review: a bad timestep_fs
    # silently clamped n_mix_steps to 1 and only failed deep in the first cycle).
    def test_nonpositive_timestep_rejected(self):
        with pytest.raises(ValueError, match='timestep_fs'):
            MixConfig(mix_time_ps=1.0, settle_steps=0, temperature_K=300.0,
                      timestep_fs=-1.0)
        with pytest.raises(ValueError, match='timestep_fs'):
            MixConfig(mix_time_ps=1.0, settle_steps=0, temperature_K=300.0,
                      timestep_fs=0.0)

    def test_nonpositive_friction_rejected(self):
        with pytest.raises(ValueError, match='friction_per_ps'):
            MixConfig(mix_time_ps=1.0, settle_steps=0, temperature_K=300.0,
                      friction_per_ps=0.0)
        with pytest.raises(ValueError, match='friction_per_ps'):
            MixConfig(mix_time_ps=1.0, settle_steps=0, temperature_K=300.0,
                      friction_per_ps=-2.0)

    def test_bad_platform_rejected(self):
        with pytest.raises(ValueError, match='platform'):
            MixConfig(mix_time_ps=1.0, settle_steps=0, temperature_K=300.0,
                      platform='GPU')

    def test_bad_charge_method_rejected(self):
        with pytest.raises(ValueError, match='charge_method'):
            MixConfig(mix_time_ps=1.0, settle_steps=0, temperature_K=300.0,
                      charge_method='am1bcc')

    def test_valid_optional_fields_construct(self):
        cfg = MixConfig(mix_time_ps=1.0, settle_steps=0, temperature_K=300.0,
                        timestep_fs=0.25, friction_per_ps=2.0,
                        platform='Reference', charge_method='gasteiger')
        assert cfg.platform == 'Reference'
        assert cfg.charge_method == 'gasteiger'
        assert cfg.timestep_fs == 0.25
        assert cfg.friction_per_ps == 2.0


# ── MixMDConfig validation (pure; module imports without openmm) ──────────────

class TestMixMDConfigValidation:

    def test_nonpositive_n_steps_rejected(self):
        with pytest.raises(ValueError, match='n_steps'):
            MixMDConfig(temperature_K=300.0, n_steps=0)
        with pytest.raises(ValueError, match='n_steps'):
            MixMDConfig(temperature_K=300.0, n_steps=-5)

    def test_nonpositive_temperature_rejected(self):
        with pytest.raises(ValueError, match='temperature_K'):
            MixMDConfig(temperature_K=0.0, n_steps=10)

    def test_nonpositive_timestep_rejected(self):
        with pytest.raises(ValueError, match='timestep_fs'):
            MixMDConfig(temperature_K=300.0, n_steps=10, timestep_fs=0.0)

    def test_valid_config_constructs(self):
        cfg = MixMDConfig(temperature_K=300.0, n_steps=10)
        assert cfg.timestep_fs == 0.5
        assert cfg.platform == 'CPU'


# ── run() entry gating (pure, ToyCalculator; error fires before any MD) ───────

def _tiny_workflow(mixing, initial_bonds=None):
    """2-atom A-B system: the smallest runnable workflow construction."""
    template = ReactionTemplate(
        name='simple',
        groups=['A', 'B'],
        pairs=[PairSpec('A', 'B', is_formation=True, r_min=0.5, r_max=5.0)],
    )
    groups = {
        'A': ReactiveGroup('A', [0]),
        'B': ReactiveGroup('B', [1]),
    }
    config = PolymerizationConfig(
        biased_steps=5, unbiased_steps=5, n_cycles=1, seed=3,
        mixing=mixing,
    )
    wf = PolymerizationWorkflow(
        config, ToyCalculator(), template, groups,
        initial_bonds=initial_bonds,
    )
    return wf


class TestMixingGating:

    @staticmethod
    def _mix():
        return MixConfig(mix_time_ps=1.0, settle_steps=0, temperature_K=300.0)

    def test_mixing_without_initial_bonds_rejected(self):
        """mixing needs the live bond graph: no initial_bonds -> ValueError."""
        wf = _tiny_workflow(self._mix(), initial_bonds=None)
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'C'],
            cell=np.eye(3) * 20.0,
        )
        with pytest.raises(ValueError, match='initial_bonds'):
            wf.run(state)
        # fail-fast: the error fired before any MD advanced the state
        assert state.step == 0

    def test_mixing_without_cell_rejected(self):
        """mixing translates into a periodic box: cell=None -> ValueError."""
        wf = _tiny_workflow(self._mix(), initial_bonds=[(0, 1, 1.0)])
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'C'],
            cell=None,
        )
        with pytest.raises(ValueError, match='cell'):
            wf.run(state)
        assert state.step == 0

    def test_no_mixing_runs_without_bonds_or_cell(self):
        """mixing=None (default) keeps the paper-faithful loop unchanged."""
        wf = _tiny_workflow(None)
        state = SimulationState(
            positions=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            species=['C', 'C'],
        )
        logs = wf.run(state)
        assert [l.phase for l in logs] == ['biased', 'unbiased']

    def test_no_mixing_is_bit_stable(self, tmp_path):
        """The non-mixing path is bit-stable: two runs on a fixed seed give
        numerically identical final positions AND velocities, write no
        mixing.jsonl, and never enter a 'mixing'/'mix_settle' phase.

        Pins that the mixing feature is truly inert when off (mixing=None) —
        stronger than the phase-name check above. Uses only the ToyCalculator +
        VelocityVerlet toy fixture; no orb/GPU/openmm.
        """
        def run_once(out_dir):
            wf = _tiny_workflow(None)
            state = SimulationState(
                positions=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
                velocities=np.zeros((2, 3)),
                species=['C', 'C'],
                masses=masses_from_species(['C', 'C']),
            )
            logs = wf.run(state, output_dir=out_dir)
            return state.positions.copy(), state.velocities.copy(), logs

        pos1, vel1, logs1 = run_once(tmp_path / 'run1')
        pos2, vel2, logs2 = run_once(tmp_path / 'run2')

        # bit-stable dynamics
        np.testing.assert_array_equal(pos1, pos2)
        np.testing.assert_array_equal(vel1, vel2)

        # no mixing artifacts, no mixing phases
        assert not (tmp_path / 'run1' / 'mixing.jsonl').exists()
        assert not (tmp_path / 'run2' / 'mixing.jsonl').exists()
        phases1 = [l.phase for l in logs1]
        phases2 = [l.phase for l in logs2]
        assert phases1 == phases2 == ['biased', 'unbiased']
        assert 'mixing' not in phases1 and 'mix_settle' not in phases1


# ── asdict round trip (manifest provenance) ───────────────────────────────────

class TestConfigAsdict:

    def test_nested_mixing_fields_roundtrip(self):
        cfg = PolymerizationConfig(
            mixing=MixConfig(
                mix_time_ps=2.0, settle_steps=100, temperature_K=333.0,
                charge_method='gasteiger', platform='Reference',
            ),
        )
        d = asdict(cfg)
        assert d['mixing']['mix_time_ps'] == 2.0
        assert d['mixing']['settle_steps'] == 100
        assert d['mixing']['temperature_K'] == 333.0
        assert d['mixing']['charge_method'] == 'gasteiger'
        assert d['mixing']['platform'] == 'Reference'
        assert d['mixing']['timestep_fs'] == 0.5
        # JSON-serializable as the manifest requires
        json.dumps(d['mixing'])

    def test_mixing_none_stays_none(self):
        assert asdict(PolymerizationConfig())['mixing'] is None


# ── OpenMM-dependent: one-cycle integration ───────────────────────────────────

def _mix_workflow(monomer_specs=None, settle_steps=3, platform='CPU', seed=11):
    """Tiny copolymer workflow with the mixing stage enabled.

    Returns (workflow, state). mix_time_ps=0.005 at 0.5 fs = 10 classical
    steps; Gasteiger charges (offline, deterministic).
    """
    from scripts._systems import (
        _METHACRYLATE_SMILES,
        _MONOMER_SMILES,
        build_vinyl_copolymer_system,
        copolymer_initial_bonds,
    )
    if monomer_specs is None:
        monomer_specs = [(_MONOMER_SMILES, 2), (_METHACRYLATE_SMILES, 1)]
    box = 22.0
    rng = np.random.default_rng(7)
    positions, species, template, groups, pmap, chain_c_map = (
        build_vinyl_copolymer_system(
            monomer_specs=monomer_specs, n_initiators=1,
            box_size=box, rng=rng,
        )
    )
    bonds = copolymer_initial_bonds(monomer_specs, n_initiators=1)
    config = PolymerizationConfig(
        biased_steps=5, unbiased_steps=5, n_cycles=1, seed=seed,
        save_interval=1, minimize=False, equil_steps=0,
        mixing=MixConfig(
            mix_time_ps=0.005, settle_steps=settle_steps,
            temperature_K=333.0, charge_method='gasteiger',
            platform=platform,
        ),
    )
    wf = PolymerizationWorkflow(
        config, ToyCalculator(), template, groups,
        propagation_map=pmap, chain_c_map=chain_c_map,
        initial_bonds=bonds,
    )
    state = SimulationState(
        positions=positions,
        velocities=np.zeros_like(positions),
        species=species,
        cell=np.diag([box, box, box]).astype(float),
        masses=masses_from_species(species),
    )
    return wf, state


def _trajectory_phases(path):
    """All 'phase' values found in a trajectory.jsonl (headers skipped)."""
    phases = []
    for raw in path.read_text(encoding='utf-8').splitlines():
        if not raw.strip():
            continue
        rec = json.loads(raw)
        if 'phase' in rec:
            phases.append(rec['phase'])
    return phases


@requires_openmm
def test_one_cycle_mixing_integration(tmp_path):
    """One full cycle with mixing: audit record, settle frames, log order,
    finite state, manifest provenance."""
    wf, state = _mix_workflow(settle_steps=3)
    logs = wf.run(state, output_dir=tmp_path)

    # (a) mixing.jsonl: exactly one record with the audited diagnostics
    mixing_log = tmp_path / 'mixing.jsonl'
    assert mixing_log.exists()
    records = [json.loads(l) for l in
               mixing_log.read_text(encoding='utf-8').splitlines() if l.strip()]
    assert len(records) == 1
    rec = records[0]
    assert rec['cycle'] == 0
    assert rec['mix_time_ps'] == 0.005
    assert rec['n_steps_classical'] == 10   # 0.005 ps / 0.5 fs
    assert np.isfinite(rec['rms_displacement_A'])
    assert rec['rms_displacement_A'] >= 0.0
    assert np.isfinite(rec['minimized_energy_kj_mol'])
    assert np.isfinite(rec['final_energy_kj_mol'])
    assert rec['charge_method'] == 'gasteiger'
    assert rec['nagl_fallback'] is False
    # unreacted first cycle: no radicals capped, no placeholder H shed
    assert rec['n_cap_h'] == 0
    assert rec['n_placeholder_h'] == 0
    assert rec['cache_misses'] >= 1
    assert 'seed' in rec

    # (b) trajectory: settle frames present; classical mixing frames are NOT
    # (classical diagnostics stay out of trajectory.jsonl — 追補 (h))
    phases = _trajectory_phases(tmp_path / 'trajectory.jsonl')
    assert 'mix_settle' in phases
    assert 'mixing' not in phases

    # (c) CycleLog ordering: mixing then mix_settle, after unbiased
    assert [l.phase for l in logs] == [
        'biased', 'unbiased', 'mixing', 'mix_settle',
    ]
    assert logs[2].steps == 10   # classical steps
    assert logs[3].steps == 3    # MLIP settle steps

    # (d) write-back left a physically usable state
    assert np.isfinite(state.positions).all()
    assert np.isfinite(state.velocities).all()

    # (e) manifest records the nested mixing config (provenance)
    manifest = json.loads(
        (tmp_path / 'manifest.json').read_text(encoding='utf-8'))
    assert manifest['extra']['mixing']['mix_time_ps'] == 0.005
    assert manifest['extra']['mixing']['charge_method'] == 'gasteiger'


@requires_openmm
def test_mixing_run_is_deterministic(tmp_path):
    """Same seed -> identical final positions: classical seed, MB re-draw and
    settle MD all derive from the workflow rng (追補 2026-07-18 (i)).

    Reference platform: bit-stable OpenMM dynamics regardless of threading.
    """
    from scripts._systems import _MONOMER_SMILES

    specs = [(_MONOMER_SMILES, 1)]

    def run_once(out_dir):
        wf, state = _mix_workflow(
            monomer_specs=specs, settle_steps=2, platform='Reference',
            seed=11,
        )
        wf.run(state, output_dir=out_dir)
        return state.positions.copy(), state.velocities.copy()

    pos1, vel1 = run_once(tmp_path / 'run1')
    pos2, vel2 = run_once(tmp_path / 'run2')
    np.testing.assert_allclose(pos1, pos2)
    np.testing.assert_allclose(vel1, vel2)
    # and the mixing actually moved the system (velocities re-drawn at 333 K)
    assert not np.allclose(vel1, 0.0)


@requires_openmm
def test_mixing_phase_skips_gracefully_on_divergence(tmp_path, monkeypatch):
    """Layer 2 (graceful degradation, specs/decisions.md 追補 2026-07-18 (WM-P4
    de-clash)): if the classical mixing MD diverges even after de-clash + minimize
    + warm-up, the workflow SKIPS that cycle instead of crashing the whole run.

    We force an unrecoverable mix by monkeypatching ``run_mix_md`` to raise (the
    real build runs, so this exercises the actual ``_run_mixing_phase`` skip path
    including the mixing.jsonl record). The run must complete, mark mixing skipped,
    run no ``mix_settle``, and leave a finite, usable state.
    """
    import kagome.prep.mix_md as mix_md_mod

    def _boom(mix, cfg, seed):
        raise RuntimeError('forced mixing divergence (test)')

    monkeypatch.setattr(mix_md_mod, 'run_mix_md', _boom)

    wf, state = _mix_workflow(settle_steps=3)
    logs = wf.run(state, output_dir=tmp_path)          # must NOT raise

    # (a) mixing is present but recorded as a zero-step skip; NO mix_settle
    assert [l.phase for l in logs] == ['biased', 'unbiased', 'mixing']
    assert logs[-1].steps == 0
    phases = _trajectory_phases(tmp_path / 'trajectory.jsonl')
    assert 'mix_settle' not in phases

    # (b) mixing.jsonl records the skip with a reason
    records = [json.loads(l) for l in
               (tmp_path / 'mixing.jsonl').read_text(encoding='utf-8')
               .splitlines() if l.strip()]
    assert len(records) == 1
    assert records[0]['cycle'] == 0
    assert records[0]['skipped'] is True
    assert 'skip_reason' in records[0] and records[0]['skip_reason']

    # (c) the state stayed finite and usable (pre-mixing MLIP state kept)
    assert np.isfinite(state.positions).all()
    assert np.isfinite(state.velocities).all()


@requires_openmm
def test_settle_steps_zero_skips_settle_segment(tmp_path):
    """settle_steps=0: mixing runs, but no mix_settle log or frames appear."""
    from scripts._systems import _MONOMER_SMILES

    wf, state = _mix_workflow(
        monomer_specs=[(_MONOMER_SMILES, 1)], settle_steps=0,
    )
    logs = wf.run(state, output_dir=tmp_path)

    assert [l.phase for l in logs] == ['biased', 'unbiased', 'mixing']
    phases = _trajectory_phases(tmp_path / 'trajectory.jsonl')
    assert 'mix_settle' not in phases
    assert 'mixing' not in phases
    # the mixing audit record is still written
    records = [json.loads(l) for l in
               (tmp_path / 'mixing.jsonl').read_text(encoding='utf-8')
               .splitlines() if l.strip()]
    assert len(records) == 1
    assert records[0]['cycle'] == 0
