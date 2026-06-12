# Decisions log

Use this template for each decision.

## Template
- Context:
- Paper anchor:
- Decision:
- Alternatives considered:
- Scientific risk:
- Licensing/commercial impact:
- Follow-up:

## 2026-06-11: Time variable t in f1(t) resets each biased segment
- Context: Eq. 5 defines f1(t) = min(γt, f1_max). The paper does not explicitly state whether t resets at the start of each biased segment or accumulates across the entire simulation.
- Paper anchor: Eq. 5, Section 2.1 ("monotonically from 0"), biased/unbiased alternation schedule (Fig. 1).
- Decision: t is measured from the start of each biased segment and resets to 0 when a new biased segment begins. Units of t are MD time steps within that segment.
- Alternatives considered: (a) Cumulative t across entire simulation — rejected because the biased/unbiased alternation implies periodic relaxation and restart. (b) Reset on reaction event — more complex, no explicit paper support.
- Scientific risk: Low. If the paper intends cumulative t, the boost will saturate at f1_max quickly and the practical difference is minimal for long simulations.
- Licensing/commercial impact: None.
- Follow-up: Verify against paper supplementary materials if available.

## 2026-06-11: Units convention for gamma
- Context: Eq. 5 defines f1(t) = γt. With f1 in kcal/mol and t in steps, gamma has units kcal/(mol·step).
- Paper anchor: Eq. 5, Table of hyperparameters (γ = 1.0).
- Decision: gamma is in kcal/(mol·step). With gamma=1.0 and f1_max=250 kcal/mol, the ramp saturates at step 250 of a 2000-step biased segment.
- Alternatives considered: gamma in kcal/(mol·fs), which would saturate at 250 fs = 1000 steps. The difference is quantitative, not qualitative.
- Scientific risk: Medium. Affects the effective boost ramp rate. Configurable, so can be adjusted.
- Licensing/commercial impact: None.
- Follow-up: Check if paper specifies gamma units explicitly.

## 2026-06-11: numpy as sole core dependency
- Context: The TDBB equations and selection logic are purely numerical.
- Paper anchor: N/A (implementation choice).
- Decision: Use numpy as the only runtime dependency for the core package. ASE, OpenMM, etc. are optional backend dependencies.
- Alternatives considered: Pure Python (too slow for distance matrices), scipy (overkill for initial implementation).
- Scientific risk: None.
- Licensing/commercial impact: numpy is BSD-3, fully commercial-safe.
- Follow-up: None.

## 2026-06-12: JSONL as trajectory output format
- Context: Phase 2 requires trajectory output for figure regeneration. Multiple formats considered.
- Paper anchor: Figs. 2-6 require time-series data of energy, conversion, and spatial density.
- Decision: Use JSON Lines (.jsonl) with a metadata header line. Each subsequent line is one frame.
- Alternatives considered: (a) ASE Trajectory — requires ASE dependency (LGPL, license discussion needed). (b) HDF5 — requires h5py, overkill for current scale. (c) XYZ — lacks structured metadata.
- Scientific risk: None. Format is lossless and machine-parseable.
- Licensing/commercial impact: None (json is stdlib).
- Follow-up: Consider binary format (HDF5) for production-scale runs.

## 2026-06-12: matplotlib as optional plot dependency
- Context: Figure regeneration requires a plotting library.
- Paper anchor: Figs. 2-6 reproduction.
- Decision: Add matplotlib as optional dependency under `[project.optional-dependencies] plot`. Not a core dependency.
- Alternatives considered: plotly (heavier), raw SVG (impractical for scatter/line plots).
- Scientific risk: None.
- Licensing/commercial impact: matplotlib is PSF-licensed (permissive), commercial-safe.
- Follow-up: None.

## 2026-06-12: Atomic masses optional with backward-compatible default
- Context: Phase 1 used unit masses (m=1 for all atoms). Real MD requires proper atomic masses.
- Paper anchor: Standard MD methodology; paper uses OpenMM which handles masses automatically.
- Decision: Add `masses: NDArray | None = None` to SimulationState. When None, unit masses are used (backward-compatible). When provided, forces are divided by mass before integration.
- Alternatives considered: Making masses mandatory — rejected to preserve existing tests and toy use case.
- Scientific risk: Low. Users who provide species can easily derive masses.
- Licensing/commercial impact: None.
- Follow-up: Consider making masses mandatory in a future major version.

## 2026-06-12: Integrator as swappable protocol
- Context: Paper implies NVT ensemble (constant temperature) which requires a thermostat. Phase 1 used bare velocity Verlet (NVE).
- Paper anchor: Section 2 describes NVT simulations at 300-600K.
- Decision: Extract integrator into a Protocol with implementations for VelocityVerlet and Langevin. PolymerizationWorkflow accepts an optional integrator parameter, defaulting to VelocityVerlet.
- Alternatives considered: Hard-coding Langevin as the only option — rejected because Verlet is useful for testing and debugging.
- Scientific risk: Low. Default remains Verlet for backward compatibility; users opt in to Langevin.
- Licensing/commercial impact: None.
- Follow-up: None.

## 2026-06-12: Event-based bond tracking over distance-threshold inference
- Context: Need to track which bonds form/break for conversion analysis (Eq. 11-12).
- Paper anchor: Eq. 11-12 require counting formed/broken bonds.
- Decision: Track bonds as explicit attempt/confirm events. TDBB already selects specific pairs, so we record attempts at biased-phase start and confirm outcomes after unbiased relaxation using threshold_fraction × r0.
- Alternatives considered: Post-hoc distance analysis at every frame — rejected as ambiguous (what cutoff?) and expensive.
- Scientific risk: Low. Threshold fraction (default 1.2) is configurable.
- Licensing/commercial impact: None.
- Follow-up: Tune threshold_fraction against real MLIP results.

## 2026-06-12: Reactive group atom removal after confirmed formation
- Context: Once atoms form a bond, they should not be selected as candidates again.
- Paper anchor: Implicit in polymerization chemistry — reacted sites are consumed.
- Decision: After each unbiased phase, remove confirmed-formation atoms from their reactive groups.
- Alternatives considered: (a) Keep all atoms — leads to re-boosting already-bonded atoms. (b) Reclassify into new groups — more complex, deferred.
- Scientific risk: Low for formation. Dissociation events are not handled (atoms are not re-added).
- Licensing/commercial impact: None.
- Follow-up: Add re-addition logic for dissociation if needed for curing simulations.

## 2026-06-12: MACE-MP-0 as default uMLIP backend
- Context: Phase 3 requires a real uMLIP to run paper-faithful simulations. The paper uses PFP/Matlantis (proprietary, blocked_pending_review).
- Paper anchor: Section 2, "universal machine learning interatomic potential."
- Decision: Use MACE-MP-0 (small) as the default commercial-safe uMLIP backend via ASE adapter.
- Alternatives considered: (a) PFP/Matlantis — blocked, proprietary. (b) CHGNet — Apache 2.0, viable but less mature for organic systems. (c) ANI — older, limited element coverage.
- Scientific risk: Medium. MACE-MP-0 is trained on Materials Project bulk crystals, not optimized for organic/polymer systems. Quantitative agreement with PFP-based paper results is not expected. Qualitative trends (bond formation under bias, energy landscape) should be comparable.
- Licensing/commercial impact: MACE code MIT, MACE-MP-0 weights MIT. Fully commercial-safe. ASE is LGPL-2.1 (import-only, commercial-safe). PyTorch is BSD-3.
- Follow-up: Consider MACE-OFF23 for organic chemistry (ASL license — not commercial-safe as default). Consider fine-tuning MACE on polymer-relevant data.

## 2026-06-12: Using MACE instead of paper's PFP backend
- Context: The paper uses PFP (Preferred Potential / Matlantis) which is not available under a clear commercial license.
- Paper anchor: Section 2-3, all simulations use PFP as the uMLIP.
- Decision: Acknowledge that results will differ quantitatively from the paper due to different uMLIP backend. The TDBB workflow, equations, and selection logic are backend-independent per design. Qualitative trend matching is the Phase 3 acceptance criterion.
- Alternatives considered: None that are both commercially safe and directly comparable to PFP.
- Scientific risk: High for quantitative reproduction. Low for qualitative workflow validation.
- Licensing/commercial impact: Enables open, redistributable default configuration.
- Follow-up: If PFP access is confirmed by user, add as opt-in backend behind feature flag.

## 2026-06-12: MD unit system — LAMMPS 'real' style with explicit conversion
- Context: Integrators were computing accel = F/m without unit conversion. Forces are in kcal/(mol·Å), masses in amu, time in fs. The missing factor (≈4.184e-4) caused ~2390× overestimation of force-driven acceleration, making dynamics non-physical when combining force-driven and thermostat-driven motion.
- Paper anchor: Section 2 — MD simulations with timestep 0.25 fs, physical masses, NVT ensemble.
- Decision: Introduce `src/units.py` with FORCE_CONV = 4.184e-4 [Å/fs² per kcal/(mol·Å·amu)]. All integrators multiply F/m by FORCE_CONV. Langevin c2 thermal noise scale includes FORCE_CONV. Kinetic energy conversion: KE[kcal/mol] = 0.5·m·v² / FORCE_CONV.
- Alternatives considered: (a) Convert forces to eV/Å at the backend boundary — rejected, pfpoly's internal unit is kcal/mol. (b) Dimensionless reduced units — rejected, real-unit MD required for paper reproduction.
- Scientific risk: None — this is a bug fix, not an approximation. Previous MACE E2E runs produced non-physical trajectories and must be re-run.
- Licensing/commercial impact: None.
- Follow-up: None.

## 2026-06-12: Integrator split into pre_force / post_force
- Context: The single-step integrator used the same forces for both half-kicks in velocity Verlet, giving first-order velocity accuracy. Splitting into pre_force (before force computation) and post_force (after) enables proper second-order velocity Verlet with one force evaluation per step.
- Paper anchor: Standard MD methodology; Section 2 implies production-quality NVT integration.
- Decision: Integrator Protocol provides pre_force() and post_force() methods. VelocityVerlet: pre=half-kick+drift, post=half-kick. Langevin BAOAB: pre=B-A-O-A, post=B. Workflow computes forces between the two calls.
- Alternatives considered: (a) Force-computation callback passed to integrator — rejected, mixes calculator concerns into integrator. (b) Store old/new forces in integrator state — rejected, adds hidden state.
- Scientific risk: None — strictly more correct than previous implementation.
- Licensing/commercial impact: None.
- Follow-up: None.

## 2026-06-12: Minimum image convention for bias and selection distances
- Context: TDBB bias forces, candidate selection, and bond tracking computed distances from raw coordinates, ignoring periodic boundaries. This is incorrect when atoms cross cell boundaries.
- Paper anchor: Section 2 — periodic boundary conditions used in all simulations.
- Decision: Add `src/geometry.py` with `minimum_image()` for orthorhombic cells. Thread `cell` parameter through `total_bias()`, `find_candidates()`, `score_candidates()`, and `BondTracker` methods. Cell=None (non-periodic) is the default, preserving backward compatibility.
- Alternatives considered: General triclinic MIC — deferred, all current systems use orthorhombic cells.
- Scientific risk: Low. Orthorhombic MIC is exact for the current box geometry.
- Licensing/commercial impact: None.
- Follow-up: Add triclinic MIC if needed for non-cubic cells.

## 2026-06-12: Selection pair key normalization bug fix
- Context: `find_candidates()` used `(idx_a, idx_b)` as the pair_specs key, but the recursive enumeration only checked `(prev_depth, depth)` where prev_depth < depth. If group_a appeared after group_b in the template's groups list, the key would be (high, low) and the constraint would be silently skipped.
- Paper anchor: Eq. 7 — all pair distance constraints must be enforced.
- Decision: Normalize the key to `(min(idx_a, idx_b), max(idx_a, idx_b))`. Added regression tests for reversed group order.
- Alternatives considered: Storing both orderings — rejected, normalization is simpler and sufficient.
- Scientific risk: None — pure bug fix.
- Licensing/commercial impact: None.
- Follow-up: None.

## 2026-06-12: Dissociation tracking uses r0 = λ·Σr_vdW as confirmation threshold
- Context: BondTracker.check_outcomes() uses `r > threshold_fraction · r0` to confirm a dissociation. The `r0` value is inherited from `_build_pair_biases()` which sets it to `λ · Σr_vdW` (the same r0 used for the formation potential). However, Eq. 3 (dissociation potential) does not include r0 — it is V^d(r) = f1·exp(-f2·r²), centered at r=0.
- Paper anchor: Eq. 3 (dissociation potential, no r0 term); Eq. 12 (bond breaking count).
- Decision: Use r0 = λ·Σr_vdW as the confirmation threshold for dissociation. Rationale: this is the equilibrium bond length estimate already computed for the pair, so `r > 1.3·r0` is a reasonable proxy for "bond is broken."
- Alternatives considered: (a) Use the equilibrium distance from the actual MD trajectory at the time of attempt — more accurate but requires storing the attempt distance. (b) Use covalent radii sum — not directly available in current VDW_RADII table. (c) Fixed absolute threshold (e.g. 3.0 Å) — less transferable.
- Scientific risk: Medium. The threshold affects what counts as a confirmed dissociation. Current dissociation templates are unused in practice (all active reactions are formation), so this has no immediate impact. Requires tuning when dissociation reactions are introduced.
- Licensing/commercial impact: None.
- Follow-up: If dissociation tracking is activated, validate threshold against MLIP potential energy profiles.

## 2026-06-12: OrbMol-v2 as recommended backend for organic/polymer systems
- Context: MACE-MP-0 is trained on Materials Project bulk crystals and carries medium scientific risk for organic polymer systems (documented above). The paper uses PFP/Matlantis (proprietary, blocked). A better-fit open alternative is needed.
- Paper anchor: Section 2 — uMLIP must handle organic polymer chemistry (C, H, N, O).
- Decision: Add OrbMol-v2 as an optional backend, recommended for organic/polymer runs. MACE-MP-0 remains the default for general-purpose use. OrbMol-v2 is trained on OMol25 + OPoly26 (polymer-specific DFT data, ωB97M-V/def2-TZVPD), directly relevant to this project.
- Alternatives considered: (a) MACE-OFF23 — ASL license, not commercial-safe. (b) ANI-2x — limited element coverage. (c) Fine-tuning MACE — requires curated data and GPU, higher effort.
- Scientific risk: Low-medium. OPoly26 training set includes polymer-relevant chemistry. Long-range Coulomb via PME adds physics absent in MACE-MP-0. Quantitative agreement with PFP still not expected, but qualitative trends should improve.
- Licensing/commercial impact: orb-models code Apache-2.0, OrbMol-v2 weights Apache-2.0. Fully commercial-safe. Note: nvalchemiops (for periodic PME) is blocked_pending_review; non-periodic runs or runs without long-range Coulomb do not require it.
- Follow-up: (a) Test nvalchemiops license for periodic system support. (b) Benchmark OrbMol-v2 vs MACE-MP-0 on ethylene system. (c) Windows compatibility is not guaranteed by upstream — verify before making it default anywhere.
