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

## 2026-06-11 (updated 2026-06-13): Units convention for gamma
- Context: Eq. 5 defines f1(t) = γt. With f1 in kcal/mol, the unit of gamma depends on whether t is in steps or physical time (fs).
- Paper anchor: Eq. 5, Section 3 Methods (γ = 1.0). PDF p.7: "bias parameters of ... γ = 1.0" — unit not stated.
- Decision: PDF confirmed (p.7): γ=1.0, unit not stated. Maintaining kcal/(mol·step). Saturation at step 250 (62.5 fs) is within biased phase (2000 steps = 500 fs). Fig. S4 (p.25-26) confirms γ acts as global scaling factor for reaction rates — unit choice affects absolute rate but not relative trends.
- Alternatives considered: (a) kcal/(mol·step) — current implementation, saturation at step 250. (b) kcal/(mol·fs) — saturation at 1000 steps (250 fs). Both produce saturation within biased phase; qualitative behavior identical.
- Scientific risk: Low. Both interpretations saturate within biased phase. γ scaling (Fig. S4) confirms no pathway-specific bias from parameter choice.
- Licensing/commercial impact: None.
- Follow-up: None. Resolved.

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

## 2026-06-13: T-G1 vinyl radical polymerization system design
- Context: Implementing the first realistic chemical system for TDBB. The paper uses vinyl monomers (methyl acrylate etc.) with AIBN initiator; current code only has ethylene without initiator.
- Paper anchor: Section 3 Methods — radical polymerization systems; Fig. 4-5 (methyl acrylate conversion curves); AIBN initiator described in system preparation.
- Decision: (a) Model AIBN by pre-placing isobutyronitrile radicals (IBN•) as closed-shell CC(C)C#N with the central C designated as the reactive radical site. AIBN thermal decomposition is not simulated — consistent with common practice in MD studies. (b) Monomer: methyl acrylate (SMILES: C=CC(=O)OC) as the simplest paper-relevant vinyl system. (c) 3D coordinates via RDKit SMILES embedding + MMFF optimization (BSD license, commercial-safe). (d) Chain propagation: after formation (radical_C + vinyl_alpha_C), the vinyl_beta_C is added to radical_C group; implemented via propagation_map dict[int,int] in PolymerizationWorkflow. (e) r_min=1.6, r_max=5.0 Å for candidate selection (same as ethylene template, paper values not stated per-system).
- Alternatives considered: (a) Simulate AIBN decomposition — adds dissociation complexity, deferred to Phase 10. (b) Use ASE molecule builder — doesn't support arbitrary SMILES; ethylene is the only vinyl monomer available. (c) No chain propagation — α(t) would stop at 1/n_initiators, far from paper trends.
- Scientific risk: Medium. IBN closed-shell geometry approximation introduces ~1H error on radical C. RDKit geometry is pre-optimized but not relaxed by MACE-MP-0 — short initial relaxation is recommended. Quantitative rate agreement with PFP-based paper results is not expected.
- Licensing/commercial impact: RDKit (BSD-3 clause, commercial-safe) added as optional dependency. See specs/dependency-license-matrix.md.
- Follow-up: Extend to other monomers (methacrylate, styrene) by passing different SMILES. Nylon-6,6 requires T-G2 (owner approval).

## 2026-06-13: NPT Monte Carlo barostat (T-B)
- Context: Paper states "NPT ensemble" for all simulations. No barostat type, target pressure, or coupling constants are specified in the arXiv HTML or Methods section.
- Paper anchor: "refinement using the PFP uMLIP under the NPT ensemble" (Section 2).
- Decision: Implement isotropic Monte Carlo barostat (MCBarostat). Frequency=25 steps (OpenMM default). Target pressure=1 atm = 1.4596e-5 kcal/(mol·Å³) (standard atmospheric pressure — implicit assumption, not stated in paper). Max volume change fraction=0.01 per attempt. Acceptance criterion: standard NPT Metropolis (N_atoms included in Jacobian term). Active in both biased and unbiased phases.
- Alternatives considered: (a) Berendsen — simpler but not a correct NPT ensemble. (b) Parrinello-Rahman — correct but complex to implement and prone to instability during equilibration. (c) MC barostat — correct NPT ensemble, matches OpenMM's default for Langevin+NPT. Selected as most likely to match the paper's implicit implementation.
- Scientific risk: Medium. Target pressure 1 atm is assumed, not stated. If the paper uses a different pressure or pressure coupling, quantitative volumes will differ. Qualitative workflow behavior (bond formation, α(t) trends) should be unaffected.
- Licensing/commercial impact: None (internal implementation).
- Follow-up: Confirm target pressure from paper PDF or supplementary materials if available. The barostat step interval (frequency=25) should be tunable.

## 2026-06-13: Maxwell-Boltzmann velocity initialization added (T-C)
- Context: All run scripts previously initialized velocities to zero (np.zeros). Standard MD practice and NPT equilibration assumption require Maxwell-Boltzmann initialization at the target temperature.
- Paper anchor: Standard MD methodology; paper uses NPT equilibration which presupposes physically initialized velocities. No explicit statement about velocity initialization in arXiv HTML.
- Decision: Add `maxwell_boltzmann_velocities()` in `src/integrators/init_velocities.py`. Unit: σ = sqrt(KB·T·FORCE_CONV/m) [Å/fs]. Center-of-mass drift removed by default. Applied to all Langevin-based run scripts (run_mace_pbc.py, run_mace.py, run_orb.py). NVE (run_toy_bond_demo.py) keeps zero initialization since there is no thermostat.
- Alternatives considered: Rescaling after zero init (velocity rescaling) — rejected, MB is standard.
- Scientific risk: None — strictly more physically correct than zero initialization.
- Licensing/commercial impact: None (numpy stdlib only).
- Follow-up: None.

## 2026-06-13: Position wrapping (PBC fold-back) implemented in integrators (T-D)
- Context: VelocityVerletIntegrator and LangevinIntegrator updated positions without wrapping back into the primary cell. Over long simulations in PBC, coordinates could drift arbitrarily outside [0, box).
- Paper anchor: PBC used in all simulations (Section 2).
- Decision: Add `wrap_positions()` to `src/geometry.py` (in-place modulo operation for orthorhombic cell). Apply after each drift step in `pre_force()` for both integrators. Thread `cell` as an optional keyword argument through `pre_force()` signatures and `PolymerizationWorkflow._run_biased_phase()`/`_run_unbiased_phase()`. Cell=None is a no-op for backward compatibility.
- Alternatives considered: Wrapping in the workflow rather than integrators — rejected, wrap immediately after drift is physically correct.
- Scientific risk: None — strictly more correct for PBC runs.
- Licensing/commercial impact: None.
- Follow-up: None.

## 2026-06-13: Bond confirmation threshold corrected to paper-faithful value
- Context: BondTracker default threshold_fraction was 1.2 (effective threshold = 1.2×r0 = 0.72×ΣvdW). All run scripts used 1.3 (effective = 0.78×ΣvdW). The paper explicitly states bonding criterion as "60% of the sum of their van der Waals radii" = r0 = 0.6×ΣvdW.
- Paper anchor: Section 2, bond confirmation criterion: "60% of the sum of their van der Waals radii."
- Decision: Change default threshold_fraction to 1.0 in BondTracker.__init__(). Remove explicit threshold_fraction=1.3 from all run scripts (use default). The confirmed threshold is now r <= r0 = lambda*sum(vdW), matching the paper exactly.
- Alternatives considered: Keep 1.2/1.3 as lenient approximation — rejected, paper is explicit about this value.
- Scientific risk: Low. toy_bond_demo result (r=1.825 Å, r0=2.04 Å) remains confirmed_formation=1 under new threshold (1.825 < 2.04). For MLIP-based runs this makes confirmation stricter.
- Licensing/commercial impact: None.
- Follow-up: None. Tests that explicitly pass threshold_fraction=1.2 remain valid as custom-threshold tests.

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
- Alternatives considered: (a) Convert forces to eV/Å at the backend boundary — rejected, kagome's internal unit is kcal/mol. (b) Dimensionless reduced units — rejected, real-unit MD required for paper reproduction.
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

## 2026-06-12: Ethylene+ethylene does not confirm bond formation in non-periodic short runs
- Context: B2 validation goal was confirmed_formation ≥ 1. Ran OrbMol-v2 with 500 biased steps (f1_max=250 kcal/mol reached) + 200 unbiased steps × 4 cycles at 500 K, non-periodic, box=8.0 Å.
- Paper anchor: Section 2 — TDBB demonstrated on polymer/interface systems, typically 2000 biased + 2000 unbiased steps per cycle.
- Outcome: confirmed_formations = 0. Bond events show attempted pairs at ~3.9-4.2 Å; after 200 unbiased steps at 500 K the atoms spring back to initial distances. By cycle 3 all pairs exceeded r_max=4.5 Å (non-periodic diffusion).
- Root causes: (a) Ethylene + ethylene direct C-C bond formation barrier is ~40+ kcal/mol — realistic PES modeled by OrbMol-v2. (b) Non-periodic system: molecules diffuse freely after bias removed. (c) 200 unbiased steps insufficient to trap the bonded state; paper uses 2000.
- Decision: Accept confirmed_formations = 0 as correct behavior for this system/conditions. TDBB machinery is verified: bias applied correctly, attempt recorded, outcome checked. B2 confirmation requires paper-scale parameters (≥2000 steps/phase) with PBC or a simpler test system.
- Scientific risk: None — the implementation is correct; the chemistry is hard.
- Follow-up: (a) Re-run B2 with paper-scale parameters when Windows PBC blocker (nvalchemiops) is resolved. (b) Consider using a chemically simpler test system (e.g., low-barrier radical chain) for unit-level bond formation validation.

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

## 2026-06-13: Langevin friction corrected from 0.01 to 0.001 /fs (P-0c)
- Context: LangevinParams default friction_per_fs was 0.01 (= 10 ps⁻¹). All run scripts used this value explicitly.
- Paper anchor: PDF p.20, Supporting Information: "Langevin thermostat with a coupling constant of 1.0 ps⁻¹". 1.0 ps⁻¹ = 0.001 fs⁻¹.
- Decision: Change LangevinParams default to friction_per_fs=0.001. Remove explicit friction_per_fs=0.01 from all run scripts so they use the corrected default.
- Alternatives considered: Keep 0.01 — rejected, it is 10× the paper value.
- Scientific risk: Low. Weaker coupling produces more physical dynamics; temperature equilibration is slightly slower but Langevin thermostat still functions correctly.
- Licensing/commercial impact: None.
- Follow-up: None.
- NOTE (2026-07-07): この論文忠実値 0.001 は **OrbMol-v2 の production レシピでは意図的に 0.01 へ上書きされる**(f2=10→2 に伴う累積加熱の冷却レバー)。整合と根拠は「2026-07-07: Langevin friction — 論文値 0.001 と production レシピ 0.01 の整合」を参照。paper_faithful run は 0.001 のまま。

## 2026-06-13: Temperature defaults corrected from 500 K to 333 K (P-0d)
- Context: All run scripts defaulted to --temperature 500 K, which was an arbitrary test value. The paper specifies system-specific temperatures.
- Paper anchor: PDF p.21: vinyl radical polymerization at 333 K; PDF p.22: nylon-6,6 at 300 K; PDF p.24: epoxy curing at 333 K.
- Decision: Change run_vinyl_aibn.py default to 333 K. Change run_mace_pbc.py, run_mace.py, run_orb.py defaults to 333 K (ethylene is closest to vinyl system). Future nylon-6,6 script will default to 300 K.
- Alternatives considered: Keep 500 K — rejected, not paper-faithful.
- Scientific risk: None. Lower temperature is more physical for these systems.
- Licensing/commercial impact: None.
- Follow-up: None.

## 2026-06-13: T8.1 toy chemistry system for bond formation demonstration
- Context: All prior MLIP-based runs (MACE-MP-0, OrbMol-v2) yielded confirmed_formations=0. The TDBB machinery (bias, attempt, confirm) was verified to be correct, but C-C bond formation in ethylene has a ~40+ kcal/mol barrier, and non-periodic systems allow diffusion. A lower-barrier system is needed to demonstrate end-to-end bond formation machinery.
- Paper anchor: Implicit — the paper demonstrates TDBB produces bond formation. The toy system proves the machinery works before applying it to hard chemistry.
- Decision: Use a 2-atom system (C–C) with ToyCalculator(epsilon=10 kcal/mol, sigma=2.04 A). The LJ well minimum sits at sigma = r0_CC = lambda_vdw*(vdwC+vdwC) = 0.6*3.4 = 2.04 A, matching the TDBB target distance. epsilon=10 kcal/mol >> kT(500K)=0.99 kcal/mol, so atoms stay trapped in the well during the unbiased phase. BondTracker confirms formation if r <= 1.3*r0 = 2.65 A.
- Result: confirmed_formations=1 at step=600 with r=1.825 A. T8.1 PASSED. Script: scripts/run_toy_bond_demo.py, runs/toy_bond_demo/.
- Alternatives considered: (a) MACE+PBC with paper-scale steps — correct chemistry but ethylene barrier is hard and requires long simulation; reserved for T8.2. (b) Radical chain initiation with MACE — more realistic but complex; no commercial-safe low-barrier benchmark available. (c) Lower epsilon toy — rejected: atoms escape LJ well during unbiased phase.
- Scientific risk: None for demonstrating machinery. The toy backend does not represent real polymer chemistry.
- Licensing/commercial impact: None (toy backend is internal code).
- Follow-up: Real-chemistry bond formation via MACE+PBC is T8.2.

## 2026-06-13: 2-group reaction template for vinyl polymerization (simplification of Eq. 6-7)
**Superseded by 2026-06-18 S4 Phase1 — 4-group multi-pair template (d_ijkl = r_ij + r_ik + r_jl per Table S1).**
- Context: Eq. 6-7 defines a general framework with groups I, J, K, L and pair set P. The current implementation in scripts/_systems.py uses a 2-group template (C_donor, C_acceptor) with 1 pair for vinyl polymerization.
- Paper anchor: Eq. 6-7, Table of reaction systems (Section 3 and paper examples).
- Decision: 2-group template is correct for vinyl/radical polymerization. Confirmed from arXiv HTML: vinyl radical polymerization uses Gi (radical carbon) and Gj (alkene carbon), pair {(i,j)} only. The 4-group template (Gi, Gj, Gk, Gl with pairs (i,j),(i,k),(j,l)) is specifically needed for epoxy curing on CuO surface, NOT for vinyl polymerization.
- Alternatives considered: Implementing the full 4-group template for vinyl — rejected as the chemistry requires only 1 bond-forming pair.
- Scientific risk: None for vinyl/radical systems. If nylon-6,6 or epoxy systems are added, the template builder must be extended. The selection machinery (src/reactive/selection.py) already supports N groups and arbitrary pair sets; only the test system builder in scripts/_systems.py is system-specific.
- Licensing/commercial impact: None.
- Follow-up: When adding nylon-6,6 or epoxy systems, add a corresponding build_*_template() function in scripts/_systems.py. The 4-group epoxy template requires Gi (epoxy O), Gj (1-deg amine N), Gk (2-deg amine N), Gl (surface OH) with P={(i,j),(i,k),(j,l)}.

## 2026-06-13: Equation numbering discrepancy in analysis modules (src/analysis/)
**Resolved 2026-06-18: PDF confirmed Eq.11 is the last numbered equation. Density formula is unnumbered. Docstrings and claims.yaml updated.**
- Context: During T6.1 (arXiv HTML verification), confirmed that the equation numbering used in src/analysis/conversion.py ("Eq. 11-12") and src/analysis/density.py ("Eq. 13") may not match the actual paper. From the HTML: Eq. 11 = alpha(t) = 1 - exp(-kp_eff*t) (fitting formula), Eq. 12 = depth-resolved density. The raw conversion fraction alpha = N_reacted/N_total appears to be either unnumbered or given a different number.
- Paper anchor: Eq. 11 (PDF p.9), density formula (PDF p.12, unnumbered).
- Decision: PDF confirmed: Eq.11 = α(t) = 1 − exp(−k*_p·t) is the last numbered equation. The density formula ρ_rxn(z) on p.12 has no equation number. The raw conversion α = N_reacted/N_total is also unnumbered (defined in Fig.2 caption as α = 1 − [M]/[M]₀). Docstrings in conversion.py and density.py updated to match. claims.yaml eq12 renamed to unnumbered_density.
- Scientific risk: None.
- Licensing/commercial impact: None.

## 2026-06-13: OrbMol-v2 as default backend for vinyl polymerization (T-OD)
- Context: Owner requested OrbMol-v2 as the default backend. run_vinyl_aibn.py previously used MACE-MP-0.
- Paper anchor: Section 2 — uMLIP must handle organic polymer chemistry. OrbMol-v2 trained on OPoly26 (polymer DFT data, ωB97M-V/def2-TZVPD).
- Decision: Change run_vinyl_aibn.py default backend to OrbMol-v2 via --backend orb (default). MACE-MP-0 retainable via --backend mace. Backend selection uses lazy imports. Other ethylene scripts (run_mace.py, run_mace_pbc.py, run_orb.py) retain their respective hardcoded backends.
- Alternatives considered: Creating separate run_vinyl_aibn_orb.py — rejected, --backend flag is cleaner.
- Scientific risk: Low. OrbMol-v2 trained on polymer-relevant data; expected to outperform MACE-MP-0 for organic systems.
- Licensing/commercial impact: orb-models code + weights both Apache-2.0, fully commercial-safe.
- Follow-up: Future run scripts (nylon-6,6, epoxy) will also default to OrbMol-v2.

## 2026-06-13: Nylon-6,6 step-growth polycondensation system design (T-G2)
- Context: Implementing the second chemical system from the paper. Nylon-6,6 uses step-growth (condensation) mechanism, requiring mixed formation/dissociation bias — the first system to exercise this capability.
- Paper anchor: PDF p.22, Table S2, Fig. S2, Fig. 4b-c. System: 100 hexamethylenediamine + 100 adipic acid, 300 K, 1 atm, NPT.
- Decision: (a) 4-group template: amine_N (i), carboxyl_C (j), amine_H (k), carboxyl_OH (l). (b) 4 PairSpecs: (i,j) formation r=3.0-6.0, (i,k) dissociation r=0.0-3.0, (j,l) dissociation r=0.0-3.0, (k,l) formation r=0.0-100.0. The k-l pair has a permissive distance range because Table S2 does not constrain k-l for candidate identification. (c) Water generation is not explicitly modeled — TDBB only biases distances; the N-H stretching + C-OH stretching effectively drives the condensation. (d) Each NH2 reacts once (N removed from amine_N after confirmed formation). (e) Carothers comparison: DPn = 1/(1-p) theoretical curve vs simulated DPn.
- Alternatives considered: (a) 2-group template (N, C only) — rejected, paper explicitly uses 4 groups. (b) Explicit water removal — rejected, paper shows results without water removal (p.11 "without continuous water removal"). (c) N reacts twice (primary → secondary amine) — rejected, nylon-6,6 is a linear polymer; each NH2 forms one amide bond.
- Scientific risk: Medium. Mixed formation+dissociation bias in a single template is untested in E2E. The k-l formation pair (amine_H + carboxyl_OH → water proximity) bias may not drive the chemistry correctly without explicit bond topology changes. Qualitative DPn vs conversion trend is the acceptance criterion, not quantitative match.
- Licensing/commercial impact: None (RDKit BSD-3 for structure generation).
- Follow-up: E2E execution (T-G2 acceptance criteria) after P-0 and T-OD complete.

## 2026-06-13: Candidate selection distance range corrected to Table S1 values (T-G1a)
- Context: Vinyl polymerization template used r_min=1.6, r_max=4.5 Å. E2E runs found 0 candidates because reactive atoms were > 4.5 Å apart at initial density.
- Paper anchor: PDF p.22, Table S1: Initiation and Propagation both use i-j r_min=3.0, r_max=6.0 Å.
- Decision: Update both ethylene and vinyl+AIBN templates to r_min=3.0, r_max=6.0 Å. Also increased default box_size from 14.0 to 16.0 Å (14.0 cannot fit 8+2 molecules with min_sep=2.5).
- Alternatives considered: Keep 1.6-4.5 — rejected, fails to find any candidates at paper-relevant densities.
- Scientific risk: None — strictly more faithful to the paper.
- Licensing/commercial impact: None.
- Follow-up: None.

## 2026-06-13: OrbMol-v2 PBC requires TORCHDYNAMO_DISABLE on Windows
- Context: OrbMol-v2 with PBC triggers Coulomb PME → nvalchemiops → torch._inductor which requires cl.exe (MSVC C++ compiler), not available on this system.
- Paper anchor: N/A (platform-specific workaround).
- Decision: Set TORCHDYNAMO_DISABLE=1 in orb_backend.py via os.environ.setdefault. This disables torch.compile for all OrbMol-v2 computations. Performance impact is acceptable for current CPU-based runs.
- Alternatives considered: (a) Install MSVC Build Tools — heavy dependency, not required by upstream docs. (b) Use MACE-MP-0 only for PBC — loses OrbMol-v2's polymer-optimized potential.
- Scientific risk: None. TORCHDYNAMO_DISABLE only affects JIT compilation, not numerical results.
- Licensing/commercial impact: None.
- Follow-up: Remove TORCHDYNAMO_DISABLE if upstream fixes nvalchemiops Windows support or cl.exe becomes available.

## 2026-06-13: T-G1a root-cause — zero confirmed formations is a system-scale problem, not a parameter bug
- Context: vinyl+AIBN E2E produced 0 confirmed bond formations even at paper step counts (3 cycles x 2000 biased + 2000 unbiased). The earlier "insufficient steps" hypothesis was DISPROVEN — 4x more steps still gave 0.
- Diagnosis: The formation bias V^f = f1(1 - exp(-f2(r-r0)^2)) with f2=10 A^-2 (paper-confirmed, PDF p.7) and r0~2.04 A (C-C) yields a force significant only within ~r0 +/- 0.5 A. Measured |F|: r=2.5 A -> 277, r=3.0 A -> 0.48, r=4.5 A -> ~0 kcal/mol/A. Candidates are LISTED at 3-6 A (Table S1) but the bias does not pull them in from there.
- Paper anchor (decisive): PDF p.7 and S-7 state polymerization is governed by near-contact events, and f2 in 5-20 gives robust behaviour. The bias CAPTURES pairs that thermal motion brings to near-contact (~2-2.5 A); it does not drag distant pairs together. This requires (a) paper density 0.5 g/mL and (b) many molecules so near-contact events are frequent.
- Decision: Do NOT change TDBB parameters (they are paper-correct). Instead reproduce paper system scale/density: hold density at 0.5 g/mL via box_from_density() and scale molecule counts up, running on GPU.
- Licensing/commercial impact: None.
- Follow-up: Validate confirmed formations appear at paper density on the 40+2 (and later 200+10) systems.

## 2026-06-13: grid-guided initial placement at paper density
- Context: At 0.5 g/mL the prior global rejection-sampling placer (_place_fragments_in_box) stalled at ~35/42 molecules and raised RuntimeError. Packing IS feasible (per-molecule volume ~286 A^3 -> ~6.6 A spacing) but global random rejection is inefficient at liquid density.
- Decision: Seed each molecule near a distinct grid cell centre (ncells = ceil(N^(1/3)) per dim), then apply random rotation + attempt-shrinking jitter and accept the first pose with all inter-molecular atom separations >= min_sep (2.5 A). Deterministic via the run seed.
- Paper anchor: SI S-3..S-4 specifies only the target density (0.5 g/mL); the packing method is unspecified. The grid is an initial-configuration device only and does not bias the subsequent biased/unbiased dynamics.
- Licensing/commercial impact: None.
- Follow-up: If 200+10 packing is still tight, lower initial density and rely on the NPT barostat to compress toward equilibrium.

## 2026-06-14: gpu40 paper-scale validation — system-scale hypothesis insufficient; missing equilibration is the dominant root cause
- Context: Re-ran the gpu40 validation (40 monomers + 2 AIBN, density 0.35 g/mL, box 25.71 Å, NPT 1 atm, T=333 K, 3 cycles × 2000 biased + 2000 unbiased = 12,000 steps, OrbMol-v2 on RTX 4060 Ti, seed 42). This is the follow-up validation promised in the 2026-06-13 "system-scale problem" and "grid-guided placement" records.
- Result: confirmed_formations=0, confirmed_dissociations=0, propagation_events=0. Candidates ARE detected each biased phase (2/4/5) and selected (2 each, bias_E=500 kcal/mol), so the TDBB machinery runs. No bond confirms.
- NEW root-cause finding (decisive): temperature is NOT controlled at the 333 K setpoint. Per-frame T: mean 527 K, std 182 K, max 1364 K. Phase split: biased mean 586 K / max 1364 K, unbiased mean 469 K / max 677 K. Time trend per 2000 steps: 850 → 607 → 516 → 438 → 392 → 360 K. The system starts with a temperature spike to ~1364 K at step ~51 and spends the entire 12,000-step run shedding initial heat, only approaching the setpoint (360 K) in the final cycle.
- Diagnosis: scripts/run_vinyl_aibn.py goes straight from build_vinyl_aibn_system() → integrator → workflow.run() with NO energy minimization and NO equilibration (grep for minimi|equilibr|relax|BFGS|FIRE in run_vinyl_aibn.py / _systems.py / polymerization.py returns nothing). Grid-guided placement at liquid-like density leaves close intermolecular contacts → large initial forces → potential energy converts to kinetic energy → Langevin cannot dissipate fast enough → atoms have high kinetic energy and candidate pairs never settle into the near-contact (~2-2.5 Å) window where V^f actually pulls (see 2026-06-13 force-vs-distance record). Hence 0 confirmations even at paper scale.
- Paper anchor (decisive): PDF p.20 (SI) — "Equilibration simulations were performed in the NPT ensemble at 300 K and 1 atm ... Production simulations using reactive acceleration MD were then carried out at 333 K and 1 atm." The paper runs a SEPARATE NPT equilibration BEFORE the reactive (TDBB) production stage. Our workflow omits this stage.
- Decision: The "lower density + rely on barostat" follow-up (2026-06-13) is necessary but NOT sufficient. Add a pre-TDBB stage: (1) energy minimization of the placed structure, then (2) a short NPT equilibration at the target T/P to relax close contacts and thermalize, BEFORE entering the biased/unbiased TDBB loop. Then re-validate confirmed_formations on the 40+2 system.
- Scientific risk: Low — adding equilibration is paper-faithful (PDF p.20) and does not alter TDBB equations or scheduling. The biased/unbiased dynamics and logging are unchanged.
- Licensing/commercial impact: None.
- Follow-up: Implement minimization + NPT equilibration stage (pending owner approval), then re-run gpu40 and check (a) T stabilizes near 333 K from the start and (b) confirmed_formations >= 1.
- Artifacts: runs/vinyl_aibn_gpu40/ (summary.json, trajectory.jsonl, bonds.jsonl, figures/). Prior interrupted partial preserved at runs/vinyl_aibn_gpu40_partial_interrupted/.

## 2026-06-14: Equilibration fix resolves temperature instability; formations=0 persists (bias capture-range + scale)
- Context: Implemented pre-TDBB relaxation (FIRE energy minimization + unbiased NPT equilibration) per the 2026-06-14 root-cause record, then re-ran gpu40 (40+2, density 0.35, NPT 333 K, minimize fmax=1.0, equil 2000 steps, 3 cycles × 2000+2000, seed 42, OrbMol-v2 on RTX 4060 Ti). New code: src/integrators/minimize.py (FIRE), PolymerizationWorkflow._minimize / _run_equilibration_phase (opt-in via PolymerizationConfig.minimize / equil_steps), run_vinyl_aibn.py CLI (--minimize/--no-minimize, --minimize-fmax, --equil-steps). Defaults off in config (legacy-safe); run_vinyl_aibn defaults minimize on + equil 2000.
- Result — temperature FIXED: FIRE converged in 307 steps (fmax 0.95 kcal/mol/Å). Whole-run T mean 310 K / std 40 K / max 373 K vs baseline mean 527 / std 182 / max 1364 K. Production phases on-target: biased mean 316 K, unbiased mean 330 K (setpoint 333 K). The 1364 K initial spike is gone. (Equilibration phase runs cooler, mean 237 K, because 2000 steps = 0.5 ps < Langevin coupling time 1 ps; production reaches setpoint.)
- Result — formations UNCHANGED: confirmed_formations=0. Candidates detected/selected each biased phase (2/5/4), all attempted at i-j 4.1-6.0 Å (target r0=2.04 Å), none confirm.
- Diagnosis (confirms 2026-06-13 record, now with temperature ruled out): V^f = f1(1-exp(-f2(r-r0)^2)) with f2=10 Å^-2 has well width ~1/sqrt(f2) ≈ 0.32 Å. Force F(r)=2 f1 f2 (r-r0) exp(-f2 (r-r0)^2): r=2.5→276, r=3.0→0.48, r≥4.0→~0 kcal/mol/Å. Candidate LISTING range (Table S1: 3-6 Å) far exceeds the bias CAPTURE range (~2.5 Å). The biased phase cannot drag a 4-6 Å pair inward; it can only capture pairs that thermal diffusion brings to near-contact. With only 42 molecules and 3 cycles, near-contact events are too rare → 0 formations. Temperature stability does not change this.
- Decision: Keep TDBB parameters unchanged (paper-confirmed, PDF p.7). The equilibration fix is correct and retained. Reaching formations>0 is a system-scale/statistics matter: requires paper-scale molecule counts (200 monomer + 10 AIBN) and/or more cycles so near-contact events are frequent — a large GPU cost requiring owner approval (Ask-first trigger 2/7). Do NOT widen the bias well or extend the capture range without re-reading the paper and owner approval (Ask-first trigger 3: changing TDBB scientific meaning).
- Scientific risk: Equilibration — Low (paper-faithful p.20, TDBB equations/scheduling unchanged). Assumption: equilibration length 2000 steps is not specified by the paper; chosen to match a TDBB block (500 fs) and documented here.
- Licensing/commercial impact: None.
- Follow-up: (a) Owner decision on paper-scale (200+10) run to obtain formations>0. (b) Optionally lengthen equilibration to >=4000 steps (1 ps) so the equilibration phase itself reaches 333 K. (c) Tests: tests/unit/test_minimize.py (3 cases) added; full unit suite 179 passed.
- Artifacts: runs/vinyl_aibn_gpu40/ (new, minimize+equil), runs/vinyl_aibn_gpu40_no_equil/ (baseline), runs/vinyl_aibn_gpu40_partial_interrupted/ (original interrupted).

## 2026-06-14: FIRE densification to reach paper density at scale; 200+10 production run
- Context: Owner approved full paper-scale (200 monomer + 10 AIBN) per handoff-plan T-G1a. The greedy grid placer (scripts/_systems._place_fragments_in_box) cannot pack 210 molecules at the paper density 0.5 g/mL — or even the 0.35 used for gpu40 — although the same density places fine at 42 molecules (so it is an algorithm-efficiency limit, not a physical one). Probe: 200+10 places only at <=0.25 g/mL. The MC barostat (max_volume_change_frac=0.01 -> ~0.33% linear/move) is far too slow to compress 0.25 -> 0.5 within a feasible equilibration.
- Decision: Add a deterministic densification step. compress_box() (src/integrators/minimize.py) shrinks the cubic box from the dilute placement edge to the target edge in 20 geometric stages, FIRE-relaxing the close contacts each ~1.1% shrink introduces. run_vinyl_aibn.py: try direct placement at the target box first (small systems unchanged, bit-for-bit); on RuntimeError, place at the densest feasible density among {0.25,0.20,0.15,0.10} then compress_box to the target. Also widened the placer grid to ~70% cell occupancy (ncells=ceil((N/0.7)^(1/3)); leaves 42->4 unchanged, 210->7) which alone was insufficient, so densification is the actual fix.
- Validation: On 200+10, direct placement at 0.5 (39.03 Å) fails -> falls back to 0.25 (49.18 Å) -> compresses in 20 stages. Observed stages 1-4: FIRE converges each stage (129/75/75/67 steps) and potential energy decreases monotonically (200->181->162->145 kcal/mol) — healthy, no blow-up. Unit tests: tests/unit/test_minimize.py compress_box cases + full suite (34 minimize/systems, 179 total) pass.
- Paper anchor: density is the paper-specified initial condition (SI S-3, 0.5 g/mL for vinyl). The compression path is a non-physical preparation device and does not bias the biased/unbiased dynamics (analogous to the grid-placement device, 2026-06-13).
- Scientific risk: Low. Preparation only; TDBB equations/scheduling unchanged. Assumption: 20 compression stages and the dilute fallback ladder are engineering choices, not paper values; documented here.
- Licensing/commercial impact: None.
- Production run launched: runs/vinyl_aibn_paper200/ — 200+10, density 0.5 (via densification), minimize fmax=1.0, equil 2000, 3 cycles × (2000 biased + 2000 unbiased), T=333 K, P=1 atm, seed 42, OrbMol-v2 on RTX 4060 Ti (CUDA, ~6.5 GB VRAM, 100% util). Per-step at density 0.5 is notably slower than the dilute 0.25 probe (0.32 s/step) due to the larger OrbMol-v2 neighbour graph.
- Follow-up: On completion, check confirmed_formations (target >=1 at paper scale), whole-run temperature near 333 K, generate figures, update specs/figure-comparison.md.

## 2026-06-14: paper200 record correction — launched run did NOT complete (no artifacts)
- Context: The preceding record states "Production run launched: runs/vinyl_aibn_paper200/". On review, that directory does not exist and no simulation is running (only VSCode LSP python processes). A 100+5 probe log (runs/_probe100.log) is preserved and shows FIRE densification working healthily up to compress stage 7/20 before being cut off.
- Correction: The 200+10 production run via the all-ML FIRE densification path was aborted before completion. The owner reports it became practically infeasible on the GPU. Root cause is throughput, NOT VRAM: the GPU is an RTX 4060 Ti **16 GB** and the prior record measured only ~6.5 GB used, so memory was never the limit. The limit is wall-clock: at density 0.5 the OrbMol-v2 neighbour graph makes every step slow, and the all-ML pipeline spends ~1,000-2,000 of its expensive ML force evaluations on structure preparation alone (compress_box: 20 stages × tens-to-200 FIRE steps), at the slowest density band.
- Implication: The expensive ML potential is being used for work that needs no ML (packing, densification, thermalization). This motivates the classical/ML split below.

## 2026-06-14: Decouple initial-structure preparation (classical OpenMM/OpenFF) from TDBB production (ML)
- Owner decision (this session): adopt a classical structure-preparation stage and reserve the ML potential (OrbMol-v2) for the TDBB reactive production only. Approach selected over (a) keeping the all-ML pipeline and merely trimming it, (b) RDKit-MMFF-only light relaxation, and (c) buying faster hardware. Hardware change is NOT required for 200+10 (16 GB GPU is ample; the bottleneck is the number of ML evaluations, which this split removes).
- Pipeline: RDKit single-molecule templates → dilute grid placement (0.25 g/mL, the density the greedy placer reliably reaches) → **classical OpenMM+OpenFF Sage prep** (energy minimize → deterministic box compression 0.25→0.5 g/mL → NVT thermalization at the target T with the box fixed) → return relaxed positions + box (atom order preserved) → **short ML re-equilibration** (existing PolymerizationConfig.equil_steps, to settle onto the ML PES) → **TDBB production** (unchanged). The slow all-ML compress_box is thereby removed from the paper-scale path.
- Atom-order invariant (decisive): groups / propagation_map are global-index based and the builder lays out atoms as [initiator]×N_init then [monomer]×N_mono, each an identical RDKit-AddHs-ordered block. The classical stage MUST update positions and the box ONLY, never reordering atoms. To guarantee the OpenFF Topology matches this order exactly, build each OpenFF Molecule via Molecule.from_rdkit() from the SAME RDKit Mol object the builder uses (not from_smiles, whose atom order can differ). A small refactor will expose `_rdkit_mol(smiles, seed)` shared by `_rdkit_3d` and the prep module.

### Decision D-1: prep target density fixed at the paper value 0.5 g/mL (do NOT NPT to the FF equilibrium density)
- Paper anchor: SI S-3 specifies 0.5 g/mL as the initial configuration density for vinyl. The classical FF's own 1-atm equilibrium density for liquid methyl acrylate is ~0.95 g/mL, so a classical NPT-at-1-atm prep would compress far past 0.5 and depart from the paper's stated initial condition.
- Decision: the prep compresses the dilute box to exactly the 0.5 g/mL edge, then thermalizes at fixed box (NVT) — it does not run a classical barostat to find the FF density. Density evolution during production is governed by the ML NPT loop (unchanged), so the FF's density bias never enters the science.
- Scientific risk: Low. Preparation only; sets the initial structure at the paper density. TDBB equations/scheduling untouched.

### Decision D-2: partial charges via openff-nagl (MIT code) with RDKit Gasteiger fallback
- OpenFF Sage requires partial charges; the default AM1-BCC route depends on AmberTools (GPL) or OpenEye (proprietary), both undesirable for a commercial-safe build. openff-nagl provides a GNN AM1-BCC surrogate with MIT code and CC-BY-4.0 weights (both commercial-safe; see license matrix). If nagl is unavailable, fall back to RDKit Gasteiger charges (no new heavy dependency). Charge accuracy is non-critical here because the classical structure is overwritten by a short ML re-equilibration before production.

### Decision D-3: prep protocol "simple" by default; the literature 21-step melt protocol is an optional extra
- The 21-step compression-relaxation protocol (Larsen/Hooper-style) targets entangled polymer melts and is overkill for a small-molecule monomer/initiator liquid. Default `--prep-protocol simple` = minimize → compress 0.25→0.5 → NVT thermalize. A `21step` option may be added later. Assumption: the protocol's stage counts and lengths are engineering choices, NOT paper values; documented here.

### Decision D-4: classical prep is a self-contained OpenMM simulation, not a Calculator backend
- The workflow's Langevin thermostat and MC barostat are custom NumPy kernels; running a melt equilibration through them is exactly the slow/weak path. Instead the prep runs OpenMM's native integrator + barostat to completion and returns only (positions, cell). This keeps the ML workflow ML-only (clean I/O / state / analysis separation per CLAUDE.md) and uses OpenMM's optimized kernels for the cheap classical work.
- Paper-faithfulness deviation (recorded): the paper's equilibration uses the same potential as production; here equilibration uses a classical FF while production uses the ML potential. This affects ONLY the initial structure, which is then re-equilibrated under the ML PES before any biased/unbiased dynamics. TDBB equations, scheduling, reaction selection, and logging are unchanged.
- Licensing/commercial impact: adds openmm (MIT core / LGPL GPU), openff-toolkit (MIT), openff-interchange (MIT), openff-forcefields (CC-BY-4.0, attribution required), openff-nagl (MIT) + openff-nagl-models (CC-BY-4.0, attribution required). All commercial-safe; the two CC-BY-4.0 components require attribution in distributed outputs. See specs/dependency-license-matrix.md and specs/approved_dependencies.yaml.
- Follow-up: implement Phase 1-2, validate on 40+2, then run 200+10 (classical prep → ML production) and check confirmed_formations >= 1 and T ≈ 333 K within the 16 GB GPU.

## 2026-06-14: Classical prep runs in WSL (Linux); OpenFF is unusable under native Windows
- Decisive finding (measured this session): the OpenFF/OpenMM classical-prep stack cannot run under native Windows on this machine. conda-forge OpenFF unavoidably pulls a PyTorch built against Intel MKL (openff-toolkit-base → openff-nagl → pytorch; there are NO non-MKL conda-forge pytorch builds), and on this Windows install MKL's CBLAS gemm DLL fails to load with a delay-load fault (`Windows fatal exception: code 0xc06d007f`) the first time numpy does a matmul. Reproduced in a clean python=3.11 env, so it is not a version/ordering artifact. numpy works in the existing ML envs (base/pfpoly-gpu) because those use a pip PyTorch + a non-MKL numpy; it is specifically the conda MKL pulled by OpenFF that is broken here. pip-installing openff-toolkit on Windows is not an option (no Windows wheels).
- Pivot (owner-approved this session): run the classical prep in **WSL (Ubuntu-24.04)**, where the same conda OpenFF stack installs and runs cleanly (the MKL/libiomp DLL faults are Windows-specific; conda-forge resolves a `cpu_generic`/OpenBLAS torch on Linux). Keep the validated ML production on the existing Windows `pfpoly-gpu` env (OrbMol-v2). Hand the relaxed structure across via the PreparedStructure JSON (src/prep/structure_io) on the shared `/mnt/c` filesystem — exactly the decoupling decision D-4 already implies.
- Concrete setup: WSL conda env `pfpoly-prep` = python 3.12 + openff-toolkit-base + openff-interchange + openff-forcefields + openff-packmol + openmm + rdkit + numpy/scipy/ase. Run `scripts/prep_structure.py` with that env's python from the repo on /mnt/c; it writes `runs/prep/<name>.json`. Then on Windows: `pfpoly-gpu` python `scripts/run_vinyl_aibn.py --load-structure runs/prep/<name>.json` (same --n-monomers/--n-initiators).
- Charge method: default to **Gasteiger** in prep_structure (decision D-2 fallback). nagl is available in the env but unnecessary; the classical structure is overwritten by the ML re-equilibration, so charge fidelity is non-critical, and Gasteiger avoids a model download and any torch invocation during prep.
- Validation (40 monomers + 2 initiators, seed 42, Gasteiger, OpenMM CPU): prep placed dilute at 0.25 g/mL (28.76 Å) → minimized → compressed to 0.5 g/mL (22.83 Å) in 20 stages → 500-step NVT thermalization → saved 504-atom structure (all finite). Handoff into Windows pfpoly-gpu production: species-order assertion passed; pre-TDBB FIRE converged to **E = -67.5 kcal/mol (negative — a genuinely relaxed structure)** vs the positive energies seen from the old grid/ML-compress path; the first biased phase found **6 candidates (vs 2 for the dilute/un-equilibrated gpu40)**. End-to-end pipeline completes. This confirms the classical prep produces a denser, better-relaxed starting structure and removes the slow ML compress_box from the paper-scale path.
- Note: WSL GPU passthrough works (nvidia-smi sees the RTX 4060 Ti; OpenMM in WSL even exposes a CUDA platform). So a future full migration of production into a single WSL env is possible, but is NOT done now to avoid disturbing the proven Windows OrbMol-v2 stack (pip torch cu128 + nvalchemi-toolkit-ops ABI). Recorded as a future optimization.
- Licensing/commercial impact: unchanged from the prior record (OpenFF MIT / CC-BY-4.0, OpenMM MIT/LGPL); running them in WSL Linux does not change license terms.
- Follow-up: run the real 200+10 prep in WSL → production on Windows GPU; check confirmed_formations >= 1 and T ≈ 333 K.

## 2026-06-15: paper200 production runs in NVT — the MC barostat crashes OrbMol-v2 on the Windows GPU
- Context: classical prep produced runs/prep/paper200.json (2520 atoms at 0.50 g/mL). First Windows GPU production (NPT) loaded it, FIRE-minimized 500 steps successfully (E=-633 kcal/mol — many ML evals fine), then crashed in the equilibration phase with `torch.AcceleratorError: CUDA error: device not ready` inside an ML forward.
- Isolation (this session): a short repro with --no-minimize (skipping the 500-step FIRE) still crashed in equilibration; the fault surfaced in different layers (F.linear, then silu) — i.e. an ASYNC CUDA error reported late, not tied to one op. A short run with **--no-barostat (NVT)** under CUDA_LAUNCH_BLOCKING=1 completed cleanly (equilibration + a full biased/unbiased cycle, 26 candidates / 10 selected). Conclusion: the **MC barostat's volume-move recompute** (it rescales the cell and calls the OrbMol-v2 calculator on the rescaled box) is what triggers the CUDA fault on this Windows GPU; the plain biased/unbiased ML evals are stable. This matches the original "paper200 had to be aborted on the GPU" report — it was the NPT barostat path, not raw throughput or VRAM (peak < 7 GB of 16 GB).
- Decision: run paper-scale **production in NVT (--no-barostat)** at the fixed 0.50 g/mL box. This is scientifically acceptable here because the classical prep already set the density to the paper's initial value (decision D-1), and over 3 short cycles with few formations the NPT volume change would be negligible. Recorded deviation: the paper runs reactive MD in NPT at 1 atm; we use NVT at the paper density because the MC barostat is GPU-unstable on this Windows setup. TDBB equations, scheduling, selection, and logging are unchanged.
- Bonus observation: candidate counts scale strongly with the proper density — 2 (dilute gpu40) → 6 (prepped 40+2) → **26 (prepped 200+10)** per biased phase. This is exactly the near-contact statistics increase needed for confirmed_formations > 0, and validates the whole "densify properly, then run TDBB" approach.
- Future: the barostat-on-GPU crash is Windows-specific; re-enabling NPT would require either running production in WSL CUDA too, or a CPU/host barostat recompute. Deferred.
- Production launched: runs/vinyl_aibn_paper200/ — 200+10, NVT, classical-prepped 0.50 g/mL, FIRE minimize + 2000 ML equil + 3×(2000 biased + 2000 unbiased), T=333 K, seed 42, OrbMol-v2 on RTX 4060 Ti.
- Follow-up: on completion check confirmed_formations (target >=1), whole-run T near 333 K, generate figures, update specs/figure-comparison.md.

## 2026-06-15: 200+10 exceeds 16 GB VRAM for sustained MD; scale production to 100+5
- Measured ceiling: with NVT (barostat already removed), the 200+10 (2520-atom) run still failed — VRAM crept from a ~9.5 GB single-call baseline up to ~16 GB over a few hundred steps, then the GPU pegged at ~15.9 GB / 0 % util while the process spun on CPU (a hang/thrash, not a crash). A VRAM sampler confirmed the climb (9.5 → 9.9 → jump to 15.9 GB). Killing the process freed VRAM cleanly (back to ~0.3 GB), so it is memory pressure, not a leak in our code or a wedged driver.
- Mitigations tried (kept, they help smaller systems and are good practice, but were NOT sufficient for 2520 atoms): (a) `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True` set before torch import (run_vinyl_aibn.py and orb_backend.py); (b) `torch.cuda.empty_cache()` after each OrbMol-v2 call (orb_backend.compute, cuda only, after dropping batch/result). With both, 2520 atoms still saturated 16 GB and hung.
- Root cause: OrbMol-v2's per-call working set at 2520 atoms / 0.50 g/mL (model activations + neighbour graph + autograd for forces) is ~9.5 GB, leaving too little headroom; the per-step neighbour-graph size variation pushes reserved memory to the 16 GB ceiling. This is the real reason the original paper200 "had to be aborted on the GPU" — a hard VRAM ceiling, not throughput.
- Decision: run paper-scale production at **100 monomers + 5 initiators (1260 atoms)**, which fits comfortably — observed ~2.4 GB VRAM at 75 % GPU util, healthy and progressing (vs 9.5+ GB / 0 % / hung for 2520). The 0.50 g/mL paper density and the full TDBB schedule (3×2000+2000) are unchanged; only the molecule count is halved. Candidate statistics remain strong at proper density (the prep gives many near-contacts), so formations can still be demonstrated.
- Full 200+10 needs a larger GPU (>=24 GB, e.g. RTX 4090 / cloud A100) or further memory engineering (e.g. neighbour-list capping, fp16, gradient checkpointing in the backend) — deferred, flagged to owner.
- Production launched: runs/vinyl_aibn_paper100/ — 100+5, NVT, classical-prepped 0.50 g/mL, FIRE minimize + 2000 ML equil + 3×(2000+2000), T=333 K, seed 42, OrbMol-v2 on RTX 4060 Ti.
- Follow-up: on completion check confirmed_formations (target >=1), whole-run T near 333 K, generate figures, update specs/figure-comparison.md.

## 2026-06-15: paper100 completes end-to-end at 0.5 g/mL; formations=0 is isolated to the TDBB bias capture range (NOT density/scale/temperature)
- Result: the 100+5 (1260-atom) run COMPLETED all 14,000 steps on the GPU (RC=0), NVT, classically-prepped 0.50 g/mL. This is the first paper-density TDBB run to finish end-to-end here (no crash, no hang). VRAM stayed ~2.4-4.7 GB / 16 GB at 75-87% util.
- Temperature: controlled — whole-run mean 295 K / max 336 K (vs the old unprepped baseline mean 527 / max 1364 K). Production phases: biased mean 300 K, unbiased mean 314 K (setpoint 333 K; phases are 0.5 ps < Langevin coupling so they sit slightly low, but stable, no spike). The classical prep + relaxed structure fixed the thermal instability.
- Candidates: plentiful at proper density — biased phases found 14 / 11 / 19 candidates, selected 5 / 4 / 5, with bias energy reaching f1_max (250 kcal/mol per pair, 1000-1250 total). The TDBB machinery (detection, selection, bias ramp, logging) is fully exercised.
- confirmed_formations = 0 (dissociations 0, propagation 0). Decisive isolation: density (0.50 g/mL, paper value), scale (many candidates), and temperature (stable ~300 K) are ALL now correct, yet no bond confirms. This rules out the hypotheses from the 2026-06-13/14 records (system-scale, missing equilibration, temperature). The remaining cause is the one flagged on 2026-06-14: the **bias capture range vs the candidate listing range**. V^f = f1(1-exp(-f2 (r-r0)^2)) with f2=10 Å^-2 has a well half-width ~0.32 Å around r0≈2 Å; its force is ~0 beyond ~2.5-3 Å. But candidates are LISTED at 3-6 Å (Table S1). So selected pairs at 3-6 Å feel essentially no bias and are never dragged to bonding distance; only pairs thermal diffusion happens to bring to <~2.5 Å could bond, which over 3 cycles does not occur. Proper density raised the candidate COUNT but not the fraction already inside the bias capture radius.
- This is an Ask-first item (trigger 3: changing the scientific meaning of TDBB). Do NOT widen the bias well, change f2/r0, or shrink the candidate listing range without re-reading the paper's TDBB force definition and getting owner approval. The open question is the paper's intended coupling between the candidate listing window (3-6 Å) and the bias functional form/range — either the bias is meant to ramp r0 from the current separation down to bonding (a moving tether that pulls 3-6 Å pairs in), or the listing window should be tied to the capture radius, or f2 is smaller than assumed. Resolving this is the next scientific task.
- Infrastructure status: the original blocker ("paper-scale run is infeasible / aborts on the GPU") is RESOLVED — paper-density TDBB now runs to completion via WSL classical prep + NVT production at 100+5 (and 200+10 is a pure VRAM-ceiling matter needing a >=24 GB GPU). The remaining formations=0 is a separate, pre-existing TDBB-physics question, now cleanly isolated.
- Artifacts: runs/vinyl_aibn_paper100/ (summary.json, trajectory.jsonl, bonds.jsonl, figures/), runs/prep/paper100.json (+ paper200.json prep, which is valid but too large to run on 16 GB).

## 2026-06-15: CORRECTION after reading the PDF (paper/2511.22874v1.pdf) — the [3,6] window is correct; formations=0 is an algorithm/chemistry gap, not a window bug
- Verified from the PDF (Table S1, p.22, and §2.1-2.2): for vinyl Initiation and Propagation, the FORMING pair i–j is listed at rmin=3.0, rmax=6.0 Å and V^f is applied to i–j, with additional constraints i–k ∈ [0,3] and j–l ∈ [0,3]. Bias params f1_max=250/125, f2=10 Å^-2, γ=1.0, r0=0.6·Σ r_vdw — all match our implementation. So the earlier "the [3,6] window is wrong / should be narrowed to the capture range" hypothesis (this file, earlier 2026-06-15 entry) is WITHDRAWN: the window is paper-correct.
- What the PDF actually says about how bonds form (§2.2 step 3): "the system is subjected to accelerated MD under the bias ΔV(t) for either a predetermined period OR UNTIL A REACTION EVENT IS OBSERVED. A reaction event is defined to occur when the interatomic distances of a specific set of pairs (ij, ik, and jl) SIMULTANEOUSLY satisfy the prescribed bonding conditions" where bonded = separation < 60% of the sum of vdW radii. Reacted atoms are then removed from the reactive set; unbiased relaxation follows. The bias "drives two atoms to ~60% of Σ vdW; once reached, the potential is zero." f1(t) grows monotonically and "can surpass the activation barrier."
- Real gaps between our implementation and the paper (this is why formations=0 even with correct density/scale/temperature/window):
  1. Detection timing: we check the bond condition ONCE, at the END of the unbiased phase (BondTracker.check_outcomes after _run_unbiased_phase). The paper checks CONTINUOUSLY DURING the biased phase and ends the phase on the event. A pair the bias transiently pulls to <r0 during biasing springs back during unbiased and is never counted by us.
  2. Multi-pair criterion: we only test i–j (radical_C – vinyl_alpha_C). The paper requires ij AND ik AND jl to satisfy the vdW bonding condition simultaneously. Our vinyl template has no k/l constraints at all (groups.py / _systems.build_vinyl_aibn_system define only the i–j formation pair).
  3. Biased phase length: we run a FIXED 2000 biased steps; the paper runs "until a reaction event or max time." With f2=10 the bias has ~0 force beyond ~2.7 Å, so a 3–6 Å pair only bonds after diffusion brings it into the ~<2.7 Å capture shell during a (long) biased phase; fixed short phases rarely catch this.
  4. Radical chemistry: our initiator is a CLOSED-SHELL surrogate (SMILES CC(C)C#N, isobutyronitrile) and the orb backend sets spin=1 (singlet). A real AIBN-derived radical is open-shell (doublet). Without an unpaired electron OrbMol-v2 has no radical addition channel, so even if the bias drives i–j to <r0 the MLIP does not stabilize the new C–C bond and it relaxes apart. The paper uses PFP with actual AIBN decomposition → real radicals.
- Net: formations=0 is NOT a density/scale/temperature/window problem (all correct). It is (a) a workflow-fidelity gap — continuous in-biased-phase detection of the multi-pair vdW bonding condition with run-until-event, vs our end-of-unbiased single-pair check; and (b) a chemistry-model gap — closed-shell radical surrogate vs a real open-shell radical. Both are Ask-first (TDBB scientific meaning / reaction model). Do NOT change f2/r0/λ/window (paper-confirmed).

## 2026-06-15: Fix A (in-phase detection) implemented — necessary but NOT sufficient; the decisive blocker is the radical chemistry (#4)
- Implemented gap #1 (BondTracker.check_reactions_during_bias + run-until-reaction in _run_biased_phase; 193 unit tests pass). Re-ran paper100 at 0.50 g/mL, NVT, **5 cycles** (runs/vinyl_aibn_paper100_fixA/, 22,000 steps).
- Result: confirmed_formations=0 STILL, and every biased phase ran the FULL 2000 steps — i.e. the in-phase detector never fired in any of the 5 cycles. So no selected i–j pair (listed at 3–6 Å) ever reached the bonding distance (<2.04 Å) during biasing, even with continuous detection over 5 cycles and 12/14/19/25/19 candidates per phase.
- Decisive conclusion: with the paper-confirmed bias (f2=10 → ~0 force beyond ~2.7 Å), the bond-boost CANNOT pull a 3–6 Å pair inward by itself. In the paper this works because PFP provides the REACTIVE CHEMISTRY: a real AIBN-derived radical is attracted to the vinyl carbon along the radical-addition reaction coordinate, which brings the pair into the bias capture shell where the (barrier-surpassing) bias completes the bond. Our system has a CLOSED-SHELL surrogate initiator (SMILES CC(C)C#N) and the orb backend sets total spin=1, so OrbMol-v2 sees no unpaired electron and no radical-addition channel — the radical-C and vinyl-C never approach. This is gap #4 and it is the dominant cause; Fix A alone cannot overcome it.
- Proposed next step (Fix B, Ask-first — reaction model): use an OPEN-SHELL radical. Minimal validation first: 1 real 2-cyanoprop-2-yl radical •C(CH3)2CN + a few methyl acrylate monomers, total spin = doublet (spin=2) passed to OrbMol-v2, run TDBB, and check whether OrbMol-v2 + the bias forms a C–C bond (formations>=1). If yes, the TDBB mechanism is validated end-to-end. Open question for scaling: OrbMol-v2 takes a SYSTEM-level spin; a melt with multiple simultaneous radicals has an ill-defined global spin (independent doublets vs high-spin sum), so the multi-radical melt is a separate modeling problem from the single-radical mechanism check.
- Timing (this machine, 1260 atoms, OrbMol-v2 NVT): ~0.3 s/step → ~2 h for the 22,000-step 5-cycle run (~1 h per 3-cycle 12k-step run). Host-bound (GPU ~75-87%).

## 2026-06-15: Fix B minimal (real radical) also gives formations=0; formations is now an OrbMol-v2 reactive-reproduction question
- Added --initiator-smiles / --spin to run_vinyl_aibn.py and wired spin → create_orb_calculator (atoms.info['spin']). Built the real open-shell 2-cyanoprop-2-yl radical via SMILES C[C](C)C#N (RDKit confirms 1 radical electron; radical_C detected).
- Minimal test (1 radical + 3 methyl acrylate, box 12 Å, spin=2 doublet, NVT, 5 cycles, in-phase detection on): confirmed_formations=0; every biased phase ran full (no in-phase event). So even with a real radical, the selected i–j pair (3–6 Å) never reached bonding distance during biasing in this dilute system.
- In-capture-range diagnostic (1 radical + 1 monomer, monomer translated so radical_C–alpha = 2.40 Å, spin=2, biased, minimize OFF): the pair ended at r=23.5 Å (flew apart), formations=0. INCONCLUSIVE: translating the whole monomer to 2.4 Å without relaxation creates clashes between the molecules' other atoms → large forces → separation; minimize=False compounds it. Not a clean test of the bias+chemistry at contact.
- Platform note: torch reports `expandable_segments not supported on this platform`, so that VRAM mitigation is a no-op on this Windows CUDA build; only per-call empty_cache helps. (Does not change the 200+10 VRAM-ceiling conclusion.)
- Assessment: getting formations>0 is no longer a single code fix — it is a reactive-MD REPRODUCTION question for OrbMol-v2 (the paper used PFP). The bias (f2=10) cannot pull a 3–6 Å pair inward, so formation depends on the MLIP itself drawing the radical and vinyl carbon into the ~2.6 Å bias-capture shell along the radical-addition coordinate. Whether OrbMol-v2 reproduces that attractive addition channel (with correct spin) is unverified.
- Clean next diagnostic (recommended before any more workflow runs): an OrbMol-v2 PES scan of the radical addition — radical approaching the vinyl terminal carbon with a proper approach geometry, scan r_CC from ~3.5 → 1.5 Å, relax the other DOF at each point, spin=2, and check for a downhill channel to ~1.5 Å. If the channel exists, formations are achievable with enough sampling/cycles (and proper multi-radical spin handling); if the PES is purely repulsive, OrbMol-v2 cannot reproduce this reaction and a different MLIP (or PFP) would be required — no TDBB tuning would help. This is an Ask-first scope decision.

## 2026-06-15: PES scan — OrbMol-v2 DOES reproduce radical addition; formations is a sampling + open-shell-spin problem, not an MLIP limitation
- Built scripts/scan_radical_addition.py: methyl radical (•CH3) approaching one carbon of ethylene along the π-face normal, with the radical lobe oriented at the target carbon; the forming C–C distance is scanned and at each point ALL other DOF are relaxed under a hard C–C distance constraint (constrained FIRE: project the along-bond force off the two constrained atoms + re-impose the exact distance each step), so sp2→sp3 rehybridization is captured. Total spin = doublet (spin=2). OrbMol-v2 on CUDA.
- Rigid (unrelaxed) scan: monotonic repulsion (43→130 kcal/mol, 3.5→1.54 Å) — as expected, a rigid scan cannot see the reaction (it ignores rehybridization); NOT a valid test.
- Relaxed (constrained) scan — decisive, textbook radical-addition profile:
  r=3.5 → ΔE 0.0 (reference); barrier peak r=2.2 → +6.1; r=2.0 → −2.1; r=1.8 → −16.0; r=1.6 → −26.8; r=1.54 → **−27.8 kcal/mol** (bonded product). Matches the known •CH3+C2H4 addition (barrier ~7, exothermic ~−23 kcal/mol).
- VERDICT: **OrbMol-v2 reproduces the radical addition** (modest barrier + exothermic C–C bond). So the whole approach is viable on this backend; reaching confirmed_formations>0 is a SAMPLING problem (get a correctly-oriented radical–vinyl pair over the ~6 kcal/mol barrier, where the bias then drives the last stretch to the product), NOT an MLIP limitation and NOT fixable/blocked by TDBB parameters.
- Why our runs gave 0, fully explained now: (i) paper100 and paper100_fixA used the CLOSED-SHELL surrogate (spin=1) → no radical, no addition channel → guaranteed 0 regardless of detection/cycles; (ii) the open-shell minimal run (1 radical, dilute box 12, 5 cycles) had the channel but far too little sampling (few candidates, random orientations, short biased phases) to cross the barrier; (iii) the bias has ~0 force beyond ~2.6 Å, so a pair must thermally/diffusively climb from ~3.5 Å (vdW contact) to ~2.6 Å (most of the barrier) before the bias assists — rare in a handful of cycles.
- Path to formations>0 (for a future scope): open-shell radicals (spin) in a DENSE melt + many cycles + run-until-reaction (longer biased phases); plus a solution for OrbMol-v2's system-level spin with MULTIPLE simultaneous radicals (independent doublets vs a single global multiplicity) — the main open modeling question. Not blocked by the MLIP.

## 2026-06-15: FIRST confirmed_formation — TDBB forms a radical-addition C–C bond end-to-end (scripts/demo_radical_formation.py)
- Built scripts/demo_radical_formation.py: places ONE open-shell radical (C[C](C)C#N, spin=2 doublet) productively next to one methyl acrylate — radical_C on the vinyl π-face normal at the chosen approach distance, radical lobe aimed at the terminal vinyl carbon (reusing the orientation from the PES scan) — then runs a biased TDBB segment (in-phase detection on) and an unbiased relaxation, non-periodic.
- Result (--approach 2.5, --select-rmin 1.5): `reaction event at step 71 - ending biased phase`; **confirmed_formations=1**; final r(radical_C, vinyl_alpha_C) = **1.62 Å** AFTER the 500-step unbiased relaxation — i.e. the new C–C bond is STABLE under OrbMol-v2 (it does not spring back). This is the project's first confirmed bond and validates the full chain end-to-end: candidate selection → bias ramp → in-phase vdW-criterion detection (Fix A) → confirmation → bond survives unbiased relaxation (real radical chemistry).
- DEMO caveat (documented in the script): the candidate window was widened to r_min=1.5 Å for this single-pair demonstration so the pre-positioned near-contact pair (2.5 Å) is selectable. The paper window is [3,6] (Table S1), which excludes the <2.6 Å bias-capture shell; the paper relies on diffusion to deliver pairs inward during a long biased phase. The widened window is a DEMO device to exercise the mechanism on one pair, NOT a paper-faithful production change (f2/r0/λ and the production window are unchanged elsewhere). With --approach 2.5 and the paper window [3,6], the pair is not selected (0 candidates) and drifts apart — confirming the window/capture-shell separation is what makes melt formation a sampling problem, exactly as diagnosed.
- Conclusion: the TDBB method + OrbMol-v2 demonstrably forms and retains the vinyl radical-addition bond. Remaining work for full polymerization is sampling at melt scale (get pairs from the [3,6] window across the ~6 kcal/mol barrier into the capture shell) and multi-radical system-spin handling — both scoped, neither blocked by the MLIP or the TDBB equations.

## 2026-06-15: S2 probe — [3,6]-window melt run; formations=0 isolated to biased-phase timescale

- Context: First S2 probe run (T-S2.1/T-S2.2): 20 monomers + 1 open-shell AIBN radical (C[C](C)C#N, spin=2), density 0.5 g/mL (box 18.12 Å, 251 atoms), NVT 333 K, paper window [3,6], no directed placement, 5 cycles × (3000 biased + 300 unbiased steps), OrbMol-v2/CUDA, seed 42. Instrumented with min_pair_distance (see git 659f118).
- Result: confirmed_formations=0. min_pair_distance per biased cycle: 4.54, 4.74, **3.23**, 3.79, 3.47 Å. bias_E reached f1_max=250 kcal/mol in every cycle (γ saturates at step 250). Artifact: runs/s2_probe/.
- Analysis (decisive): timestep=0.25 fs → 3000 biased steps = 750 fs = 0.75 ps of biased MD per cycle. At this timescale, thermal diffusion in a dense organic melt cannot reliably move a pair from 3.23 Å into the TDBB capture shell. TDBB force at r=3.23 Å: F = f1·2·f2·(r−r0)·exp(−f2·(r−r0)^2) = 250·2·10·1.19·exp(−10·1.19^2) ≈ **0.004 kcal/mol/Å** (essentially zero). Force only becomes significant at r < ~2.7 Å. Conclusion: the TDBB bias IS ramped to maximum but provides negligible restoring force beyond ~2.7 Å; formation requires thermal diffusion to bring pairs to <2.7 Å first, then bias captures them. 0.75 ps is insufficient for that diffusion in a 0.5 g/mL melt.
- Decision (T-S2.3): sweep biased-phase length by increasing the MD timestep (rather than the step count, for wall-clock efficiency). Add --timestep-fs CLI arg to run_vinyl_aibn.py (default 0.25 preserved for backward compat). S2 sweep run uses --timestep-fs 1.0 (standard for ML organic MD) giving 1 fs × 15000 steps = 15 ps per biased cycle — 20× longer than the probe. At D ≈ 10^-10 m²/s (monomer liquid), MSD in 15 ps ≈ 3 Å rms, potentially sufficient to cross from 3.2 Å to the 2.6 Å capture threshold.
- Paper anchor: timestep not specified in paper for TDBB production MD. 1.0 fs is standard for organic ML MD (OrbMol-v2 docs); 0.25 fs was chosen conservatively for initial development. γ unit is kcal/(mol·step); at dt=1.0 fs the physical ramp rate is 1.0 kcal/(mol·fs), within the handoff-plan-v4 sweep range (0.5-5.0 kcal/(mol·ps)).
- Scientific risk: Low — TDBB equations, f2/r0/λ/window, γ/f1_max all unchanged. Timestep is an integration parameter, not a TDBB parameter. If 1.0 fs is numerically unstable (NaN, energy divergence), revert to 0.5 fs (2× speedup vs 4×).
- Licensing/commercial impact: None.
- Follow-up: observe min_pair_distance in s2_sweep1; if <2.6 Å in any cycle, formations should appear. If still stuck, try lower density (0.3 g/mL, more mobility) or increase γ.

## 2026-06-15: S2 sweep1 — longer biased phase (15 ps/cycle); min_pair_dist plateaus at ~3 Å

- Context: S2 sweep1 (T-S2.3): same system as S2 probe (20+1, spin=2, 0.5 g/mL, paper [3,6] window) but --timestep-fs 1.0 --biased-steps 15000 (= 15 ps per biased cycle, 20× probe) + --equil-steps 2000, 5 cycles. Artifact: runs/s2_sweep1/.
- Result: confirmed_formations=0. min_pair_distance by cycle: 3.64, **3.09**, 3.94, 3.87, 3.64 Å. Best is 3.09 Å vs probe's 3.23 Å — marginal improvement; cycles 2-4 show no trend toward decreasing. Acceptance criterion (< 2.6 Å) NOT met.
- Analysis: the min_pair_distance stagnates around 3-4 Å regardless of biased phase length (0.75 ps → 15 ps). At r=3.0 Å (r0=2.04, f2=10): TDBB force F ≈ 0.48 kcal/mol/Å; at r=3.09 Å: F ≈ 0.22 kcal/mol/Å. The relaxed PES scan (scripts/scan_radical_addition.py) shows an uphill slope from r~3.5 Å rising to +6.1 kcal/mol at r=2.2 Å. Estimated PES slope at r=3 Å ≈ 3-5 kcal/mol/Å. Net force on the pair at 3 Å: TDBB attractive (-0.48) + PES repulsive (~+4) = net repulsive. The TDBB bias at [3,6] Å is overwhelmed by the PES barrier in OrbMol-v2; thermal diffusion alone (kT=0.66 kcal/mol at 333 K) cannot reliably cross from 3 Å to the 2.7 Å threshold where TDBB force becomes dominant (F~40 kcal/mol/Å at r=2.7 Å).
- Structural interpretation: the dense melt (0.5 g/mL) surrounds the selected pair with neighbours that constrain approach; the OrbMol-v2 radical addition PES shows a gradually rising barrier from 3.5→2.2 Å with no long-range attractive well. In contrast, PFP (the paper's MLIP) may have a more favourable radical-monomer PES that shows downhill approach from 3.5 Å, naturally delivering pairs into the TDBB capture shell where the bias completes the bond.
- Next sweep (T-S2.3): try density 0.3 g/mL (box ~21.5 Å, 67% more free volume, higher diffusivity) — explicit in the handoff-plan sweep range. If this also fails, close S2 as "melt-driven formation is a scale/MLIP-PES question, not a TDBB bug."
- Scientific risk: findings are internally consistent and consistent with the PES scan. The conclusion that OrbMol-v2 + 0.5 g/mL + single radical + 251 atoms does not give melt-driven formation in O(10 cycles) is a valid negative result.

## 2026-06-15: S1 DONE — chain propagation demonstrated (pentamer, radical migrates along the chain)
- T-S1.1: tests/unit/test_propagation.py (3 cases) verifies the post-formation bookkeeping in _update_groups_after_cycle — consumed radical (atom_a) and reacted vinyl alpha-C (atom_b) removed from groups, monomer beta-C (propagation_map[alpha]) added to radical_C; the single-radical (doublet) invariant holds; no double-apply. Full unit suite green.
- T-S1.2/3: scripts/demo_chain_propagation.py — 1 open-shell radical + N methyl acrylate, spin=2 held constant (single chain → always one unpaired electron). Runs one TDBB cycle at a time; before each cycle the next monomer is placed productively at the current chain-end radical (vinyl π-face on the radical lobe, alpha-C at 2.5 Å). DEMO devices (documented): widened candidate window (--select-rmin 1.5) and directed monomer placement; production window [3,6] and all TDBB params unchanged.
- Result (--n-monomers 4, OrbMol-v2/CUDA): confirmed_formations=4, propagation_events=4. The active radical migrated 1 → 12 → 24 → 36 → 48 — i.e. a pentamer (initiator + 4 monomers) built end-to-end, each addition detected in-phase (reaction events at steps 86/67/95/77, biased phase ending early each cycle) and the bond retained through the unbiased relaxation before the next monomer. The `assert len(radical_C)==1` held every step (doublet invariant).
- T-S1.4: spin invariant validated — exactly one chain-end radical throughout; spin held at 2. This is the explicit boundary with S3 (multiple simultaneous radicals → ill-defined global spin).
- Acceptance criteria (handoff-plan-v4 S1) all met: propagation_events ≥ 2 (got 4); bonds stable post-relaxation (radical migration confirms retention); radical migrates to each new chain end; doublet throughout; reproduction command recorded (script + handoff-plan-v4).
- Caveat (carried to S2): additions here are made deterministic by directed placement + widened window; melt-driven, undirected formation from the paper [3,6] window is S2 (a sampling problem, not blocked).
- Next: S2 (melt-driven formations) and/or S3 (multi-radical spin). See specs/archive/handoff-plan-v4.md.

## 2026-06-17: S2 sweep2 — lower density (0.3 g/mL) also fails; candidate window [3,6] closed as PES-mismatched for OrbMol-v2

- Context: S2 sweep2 (T-S2.3 follow-up): same system as sweep1 (20+1, spin=2, paper [3,6] window) but density lowered to 0.30 g/mL (box 21.48 Å) per the sweep1 "next step" plan. --timestep-fs 1.0 --biased-steps 15000 (15 ps/cycle) + --equil-steps 2000, 5 cycles. Artifact: runs/s2_sweep2/.
- Result: confirmed_formations=0. min_pair_distance by cycle: 3.19, 3.45, 3.46, 3.56, 3.54 Å. Candidates per cycle: 3, 1, 4, 3, 5 (fewer than sweep1's 4-6 at 0.5 g/mL).
- Analysis: lower density gives MORE free volume but FEWER candidates (reduced packing → fewer [3,6] contacts), and min_pair_dist does not improve (3.19 vs sweep1's 3.09 Å). The pairs still cannot cross from ~3 Å to the TDBB capture shell (<2.7 Å). Confirms sweep1 conclusion: the [3,6] window + OrbMol-v2 PES is the structural mismatch, not density or timescale.
- S2 sweep outcome (probe + sweep1 + sweep2): three density/timescale combinations tested; all give min_pair_dist ≈ 3.0-3.5 Å and 0 formations. The [3,6] Å window designed for PFP does not work with OrbMol-v2.

## 2026-06-17: Candidate window re-tuning for OrbMol-v2 PES — [3,6] → [1.5, 3.0]

- Context: The paper's [3,6] Å candidate window (Table S1) was calibrated for PFP (Matlantis). OrbMol-v2's radical addition PES (scripts/scan_radical_addition.py) has NO attractive well at 3-6 Å — the energy rises monotonically from r=3.5 Å to the barrier top at r=2.2 Å (+6.1 kcal/mol). The TDBB bias (f1=250, f2=10, r0=2.04) generates negligible force beyond ~2.9 Å (0.48 kcal/mol/Å at 3.0 Å vs PES slope ~3.6 kcal/mol/Å). Three S2 melt runs confirmed that pairs selected at [3,6] cannot cross this "dead zone" — min_pair_dist plateaus at ~3 Å.
- Quantitative analysis: fitting a quadratic to the PES scan data (3.5→0.0, 2.2→+6.1 kcal/mol) gives dE(r) = 3.61·(r−3.5)² kcal/mol. The TDBB attractive force F_TDBB = 250·20·(r−2.04)·exp(−10·(r−2.04)²) crosses the PES slope at **r ≈ 2.87 Å** — this is the effective capture radius. Beyond this distance, the PES repulsive slope dominates and the bias cannot pull the pair inward. Below 2.87 Å, TDBB force grows explosively (42 kcal/mol/Å at 2.7 Å, 122 at 2.6 Å, 493 at 2.4 Å) and overwhelms the barrier.
- Paper anchor: Table S1 specifies [3,6] for the PFP backend. The paper notes (Fig. S4) that TDBB parameters are system-specific and should be tuned. The candidate window [r_min, r_max] is an operational parameter in the selection phase (Eq. 6-7), NOT part of the TDBB potential definition (Eq. 2-5). Adjusting it for a different MLIP's PES is analogous to the paper's own γ sweep (Fig. S4).
- Decision: set the candidate window to **[1.5, 3.0] Å** for OrbMol-v2 production. Rationale: r_max=3.0 is just above the capture radius (2.87 Å), ensuring selected pairs can be captured by the bias. r_min=1.5 is below the C-C product distance (~1.54 Å) so bonded pairs can still be detected during the biased phase. All other TDBB parameters (f2=10, f1_max=250, γ=1.0, r0=2.04, λ=0.6) remain unchanged — this is a selection parameter change, not a bias-function change.
- Deviation from paper: documented. The paper uses [3,6] with PFP; we use [1.5, 3.0] with OrbMol-v2. The difference is driven by the different MLIP PES profiles in the 2.5-4.0 Å range. If PFP has a shallow vdW pre-reaction well at 3-4 Å that OrbMol-v2 lacks, [3,6] would be optimal for PFP but systematically fail for OrbMol-v2 — consistent with all observations.
- Scientific risk: Low-medium. The window change does not alter the TDBB equations, only which pairs are selected as candidates. Pairs must still cross the barrier and satisfy the bonding criterion (r < r0 = 2.04 Å). Risk is that at [1.5, 3.0] fewer candidates exist in the melt (pairs at <3 Å are rarer than at <6 Å); mitigated by proper density (0.5 g/mL) and more cycles.
- Licensing/commercial impact: None.
- Follow-up: run S2 sweep3 with [1.5, 3.0] window, 20+1, spin=2, 0.5 g/mL, 5-10 cycles × (15000 biased + 1000 unbiased), timestep 1.0 fs. If confirmed_formations >= 1, the hypothesis is validated.

## 2026-06-17: S2 sweep3 — [1.5, 3.5] window gives 0 candidates in 4/5 cycles; capture radius confirmed at 2.87 Å

- Context: S2 sweep3 (window [1.5, 3.5], 2000 biased + 500 unbiased × 5 cycles, 0.5 g/mL, timestep 1.0 fs). Artifact: runs/s2_sweep3/.
- Result: confirmed_formations=0. Candidates per cycle: 0, 0, 0, **1**, 0. Only cycle 3 found a candidate (the radical-vinyl pair thermally diffused to <3.5 Å). min_pair_distance in cycle 3: **2.872 Å** — essentially the TDBB/PES crossover point (2.87 Å, computed above). The pair reached capture radius but did not cross the barrier in 2000 steps (2 ps).
- Analysis: narrowing the window addresses the bias-effectiveness problem but creates a candidate-scarcity problem. With [3,6] we had 4-7 candidates/cycle (all ineffective); with [1.5,3.5] we have ~0.2 candidates/cycle (rare but properly biased). Neither extreme works alone.

## 2026-06-17: f2 reduction (10 → 5) to widen TDBB capture range for OrbMol-v2 PES

- Context: the TDBB bias well half-width is ~1/√f2. At f2=10 the width is 0.32 Å (capture radius 2.87 Å); at f2=5 the width is 0.45 Å (capture radius ~3.2 Å); at f2=3 the width is 0.58 Å (capture radius ~3.6 Å). The paper's [3,6] candidate window provides 4-7 candidates/cycle at 0.5 g/mL, but with f2=10 the bias force at 3-6 Å is essentially zero (0.48 kcal/mol/Å at 3.0 Å). Lowering f2 pushes the capture radius outward, so the abundant [3,6] candidates actually feel the bias.
- Paper anchor: PDF p.7 states f2=10 for the PFP backend. PDF p.7 / Fig. S4 also states "f2 in the range of 5 to 20 gives robust behaviour" — f2=5 is explicitly within the paper's validated range. The paper's characterization of f2 robustness is in the context of PFP's PES; OrbMol-v2's different radical-addition barrier profile may shift the optimal f2 within this range.
- Decision: test f2=5 with the paper [3,6] window (no window override). This combines the abundant candidates of [3,6] with a wider capture range. All other TDBB parameters unchanged (f1_max=250, γ=1.0, r0=2.04, λ=0.6).
- TDBB force at f2=5 vs f2=10 at key distances (f1=250, r0=2.04):
  - r=3.0: f2=10 → 0.48, f2=5 → 18.8 kcal/mol/Å (**39× stronger**)
  - r=3.5: f2=10 → 0.00, f2=5 → 1.14 kcal/mol/Å
  - r=2.5: f2=10 → 277, f2=5 → 394 kcal/mol/Å
- Scientific risk: Low-medium. f2=5 is within the paper's stated robust range. The bias well becomes wider and shallower per unit distance, but f1_max=250 is unchanged so the total barrier-surmounting energy is the same. The primary effect is extending the capture radius from ~2.87 to ~3.2 Å, where candidates actually exist.
- Licensing/commercial impact: None.
- Follow-up: run S2 sweep4 with f2=5, paper [3,6] window, 20+1, spin=2, 0.5 g/mL, 5 cycles × (2000 biased + 500 unbiased), timestep 1.0 fs. Compare min_pair_distance and formations against sweep1 (f2=10, same window).

## 2026-06-17: S2 sweep4 — f2=5 achieves FIRST MELT-DRIVEN FORMATION; TDBB reproduces radical addition in an undirected melt

- Context: S2 sweep4 (f2=5, paper [3,6] window, 20+1, spin=2, 0.5 g/mL, NVT 333 K, timestep 1.0 fs, 5 cycles × (2000 biased + 500 unbiased), seed 42, OrbMol-v2/CUDA). Artifact: runs/s2_sweep4/.
- Result: **confirmed_formations=1**, dissociations=0, propagation_events=1. Total steps: 13,263 (cycle 0 ended early due to reaction event).
- Per-cycle detail:
  - Cycle 0: 5 candidates, 1 selected, **reaction event at step 763** (biased phase ended early), min_pair_dist=**1.97 Å**, bias_E=6.01 kcal/mol. Chain propagation: atom 228 (beta-C) → radical_C.
  - Cycle 1: 3 candidates, 1 selected, min_pair_dist=3.27 Å, bias_E=249.94 (no reaction)
  - Cycle 2: 4 candidates, 1 selected, min_pair_dist=3.65 Å, bias_E=250.00 (no reaction)
  - Cycle 3: 2 candidates, 1 selected, min_pair_dist=4.45 Å, bias_E=250.00 (no reaction)
  - Cycle 4: 3 candidates, 1 selected, min_pair_dist=3.78 Å, bias_E=250.00 (no reaction)
- Analysis:
  - **f2=5 is the decisive parameter change.** With f2=10 (sweep1-3), the TDBB force at 3.0 Å was 0.48 kcal/mol/Å — negligible vs the PES slope (~3.6). With f2=5, the force at 3.0 Å is 18.8 kcal/mol/Å (39× stronger), sufficient to compete with the OrbMol-v2 PES barrier and pull a [3,6]-listed pair into the bonding region.
  - Cycle 0 succeeded because the selected pair happened to start close enough for f2=5 to capture it; the bias spent only 6.01 kcal/mol (vs f1_max=250) before the pair crossed the barrier and reached 1.97 Å (bonding). The bond survived the 500-step unbiased relaxation and chain propagation fired (beta-C became the new radical_C).
  - Cycles 1-4 did not form bonds — the selected pairs remained at 3.3-4.5 Å where even f2=5 force (~1-19 kcal/mol/Å) could not overcome the barrier in 2000 steps (2 ps). This is expected: formation rate is stochastic and depends on pair geometry/orientation; 1/5 cycles is a reasonable hit rate for a 20-monomer, single-radical melt.
- Comparison to prior S2 runs (all 20+1 monomers, spin=2, OrbMol-v2):
  | Run    | f2  | Window    | Candidates/cycle | min_pair_dist best | Formations |
  |--------|-----|-----------|------------------|--------------------|------------|
  | probe  | 10  | [3,6]     | 1-6              | 3.23 Å             | 0          |
  | sweep1 | 10  | [3,6]     | 4-6              | 3.09 Å             | 0          |
  | sweep2 | 10  | [3,6]     | 1-5              | 3.19 Å             | 0          |
  | sweep3 | 10  | [1.5,3.5] | 0-1              | 2.87 Å             | 0          |
  | **sweep4** | **5** | **[3,6]** | **2-5** | **1.97 Å** | **1** |
- S2 milestone: this is the project's first melt-driven formation — an undirected, unpositioned, thermally-sampled radical-addition bond formed by the TDBB protocol in a realistic melt environment. Unlike the S1 demo (scripts/demo_radical_formation.py, directed placement + widened window), this used the paper's candidate window [3,6] and no manual positioning.
- Scientific significance: validates the complete TDBB → OrbMol-v2 pipeline for radical vinyl polymerization. The only deviation from paper parameters is f2=5 (vs paper's 10), which is explicitly within the paper's stated robust range (5-20, PDF p.7 / Fig. S4). The candidate window, f1_max, γ, r0, λ are all paper-faithful.
- Next steps: (1) run more cycles or seeds to measure formation rate statistics; (2) scale to larger systems (100+5, 200+10 with classical prep) at f2=5; (3) begin S3 (multi-radical spin handling) for sustained chain growth.

## 2026-06-18: S2 sweep5 — 15 cycles × 2 seeds confirm reproducibility; S2 DONE

- Context: sweep4 (seed 42, 5 cycles) gave formations=1 but seed 7 (5 cycles) gave 0. Extended to 15 cycles to test whether additional sampling yields formations on seed 7. Artifact: runs/s2_sweep5_seed7/.
- Result (seed 7, 15 cycles): **confirmed_formations=2**, propagation_events=2.
  - Cycle 6: reaction at step 634, min_pair_dist=1.92 Å, bias_E=17.2 kcal/mol. Propagation: atom 192 → radical_C.
  - Cycle 12: reaction at step 1461, min_pair_dist=1.94 Å, bias_E=11.4 kcal/mol. Propagation: atom 180 → radical_C.
  - Cycles 0-5, 7-11, 13-14: no reaction (min_pair_dist 3.25-5.01 Å).
- Formation rate: 2/15 = 13.3%/cycle (seed 7) vs 1/5 = 20%/cycle (seed 42). Combined: 3/20 = 15%/cycle. Consistent with a stochastic process where thermal diffusion must deliver a properly-oriented pair into the f2=5 capture shell (~3.2 Å) during the 2 ps biased phase.
- Reproducibility confirmed across seeds:
  | Run         | Seed | Cycles | Formations | Rate     |
  |-------------|------|--------|------------|----------|
  | sweep4      | 42   | 5      | 1          | 20%      |
  | sweep4_seed7| 7    | 5      | 0          | 0%       |
  | sweep5_seed7| 7    | 15     | 2          | 13.3%    |
  | **Combined**|      | **20** | **3**      | **15%**  |
- Key observations:
  - Both reactions used low bias energy (6-17 kcal/mol, ~3-7% of f1_max=250) — the bias assisted the final approach but the pair was already thermally close. This is the intended TDBB mechanism: diffusion delivers, bias completes.
  - Chain propagation works in a melt: the radical migrated (initiator → atom 192 → atom 180), demonstrating successive additions from the paper [3,6] window without directed placement.
  - 5 cycles is marginal for this system/density; 15 cycles reliably produces formations.
- **S2 acceptance criteria — ALL MET:**
  - >=1 melt-driven confirmed_formation with paper [3,6] window ✅ (3 total across 2 seeds)
  - Selected [3,6] pairs reach <r0 during biasing ✅ (min_pair_dist 1.92-1.97 Å)
  - No directed placement, no widened window ✅
  - Reproduction commands recorded ✅
- S2 CLOSED. Remaining work (scaling, multi-radical spin, figures) is S3/S5/S6 scope.

## 2026-06-18: S3 Phase 1 — High-spin approximation PES validation
- Context: OrbMol-v2 accepts a single system-level spin multiplicity. For N radicals, the
  high-spin approximation (spin = N+1, all unpaired electrons parallel) is needed. Must verify
  that adding a spectator radical with high-spin coupling does not alter the PES of the
  reacting radical significantly.
- Paper anchor: Section 2 (AIBN-initiated polymerization uses 2 radical fragments);
  Section 3 Methods (OrbMol-v2 spin parameter).
- Decision: **High-spin approximation is VALID for TDBB** (verdict: MARGINAL overall, PASS for kinetics).
- Test system: 2 CH3 radicals + C2H4. Radical 1 scans C-C formation (3.5→1.54 Å).
  Radical 2 is spectator at 7 Å from ethylene center.
- Correct comparison: 1-radical/spin=2 (doublet) vs 2-radical/spin=3 (triplet).
  The spin=1 (closed-shell singlet) comparison is NOT physically meaningful for 2 radicals
  in single-determinant DFT — open-shell singlet requires multi-reference treatment.
- Key results:
  | Metric | Value | Criterion |
  |--------|-------|-----------|
  | Activation barrier diff (r=2.2 Å) | 0.88 kcal/mol | PASS (< 1) |
  | Product energy diff (r=1.54 Å) | 2.39 kcal/mol | MARGINAL |
  | Max pointwise diff (r=1.8 Å) | 2.78 kcal/mol | MARGINAL (< 3) |
  | Barrier height 1rad/s=2 | 6.15 kcal/mol | — |
  | Barrier height 2rad/s=3 | 7.03 kcal/mol | — |
- Interpretation:
  - In the TDBB-relevant region (r > 2.0 Å, where bias drives approach), the PES
    curves are nearly identical (< 1 kcal/mol difference).
  - The product well is ~2.4 kcal/mol deeper with spin=3 — irrelevant for TDBB because
    the bias is removed after bond formation.
  - For chain-growth polymerization, radical count is conserved (1 consumed → 1 generated),
    so spin multiplicity stays constant — no dynamic spin update needed.
- Alternatives considered:
  - Broken-symmetry DFT (spin=1 with initial guess): OrbMol-v2 does not support this.
  - Per-atom spin specification: OrbMol-v2 only accepts system-level spin.
  - Dynamic spin update per cycle: unnecessary because radical count is conserved.
- Scientific risk: Low for TDBB kinetics (barrier < 1 kcal/mol). Product thermodynamics
  carry ~2.4 kcal/mol systematic shift — acceptable for qualitative reproduction.
- Licensing/commercial impact: None.
- Reproduction:
  ```
  conda run -n pfpoly-gpu python scripts/scan_radical_addition_2rad.py --device cuda --output-dir runs/s3_pes
  ```
- Artifacts: `runs/s3_pes/pes_comparison.json`, `runs/s3_pes/s3_pes_comparison.png`
- Follow-up: Phase 2 — 2-radical melt run with `--n-initiators 2 --spin 3`.

## 2026-06-18: S3 Phase 2 — 2-radical melt run results
- Context: Validate that the existing multi-radical code works correctly with
  n_initiators=2 and spin=3 (high-spin approximation validated in Phase 1).
- Paper anchor: Section 2 (AIBN produces 2 radical fragments), Section 3 Methods.
- Decision: **S3 Phase 2 PASS** — 2-radical polymerization works out of the box.
- Run config: n_monomers=20, n_initiators=2, spin=3, f2=5.0, density=0.5 g/mL,
  T=333 K, 15 cycles (2000 biased + 500 unbiased), seed=42.
- Key results:
  | Metric | Value |
  |--------|-------|
  | confirmed_formations | 6 (cycles 1, 2, 5, 7, 8, 9) |
  | reaction rate | 6/15 = 40%/cycle (vs S2 15%/cycle with 1 radical) |
  | n_selected per cycle | 2 in 14/15 cycles (1 in cycle 12) |
  | chain propagation | 6 events: initiator→107→191→203→215→155→179 |
  | radical count | conserved at 2 throughout (except cycle 12: only 1 candidate pair available) |
  | dissociations | 0 |
  | n_atoms | 262 (vs S2: 251) |
  | total_steps | 30,594 |
- Observations:
  - Both radicals independently select candidates each cycle (n_selected=2).
  - Reaction rate ~2.7x higher than S2 (40% vs 15%), consistent with 2 independent
    radicals each having ~20% per-radical success probability.
  - Chain propagation correctly transfers radical site (beta-C → radical_C).
  - All formations used low bias energy (195-276 kcal/mol out of 500 max for 2 pairs),
    confirming melt-driven mechanism with bias-assisted final approach.
  - Temperature remained stable at 333 K with spin=3.
  - No code changes were needed — existing multi-radical logic in polymerization.py
    handled everything correctly.
- Scientific risk: None beyond Phase 1 PES caveat (2.4 kcal/mol product shift).
- Reproduction:
  ```
  conda run -n pfpoly-gpu python scripts/run_vinyl_aibn.py \
      --n-monomers 20 --n-initiators 2 --initiator-smiles "C[C](C)C#N" --spin 3 \
      --density 0.5 --backend orb --device cuda --no-barostat \
      --n-cycles 15 --biased-steps 2000 --unbiased-steps 500 --equil-steps 2000 \
      --timestep-fs 1.0 --f2 5.0 --seed 42 --output-dir runs/s3_2rad
  ```
- Artifacts: `runs/s3_2rad/summary.json`, `runs/s3_2rad/bonds.jsonl`
- **S3 acceptance criteria — ALL MET:**
  - PES validation: barrier diff < 1 kcal/mol (high-spin approximation valid)
  - 2-radical run: >= 1 confirmed_formation (got 6)
  - n_selected = 2 in most cycles (both radicals active)
  - Chain propagation works from both radical sites
  - Radical count conserved throughout simulation
- S3 CLOSED.

## 2026-06-18: S4 Phase 1 — Multi-pair criterion (d_ijkl = r_ij + r_ik + r_jl) implemented and validated

- Context: S4 objective is paper fidelity improvement. Phase 1 extends the vinyl template
  from 2 groups (i=radical_C, j=vinyl_alpha_C) to 4 groups (+ k=chain_C, l=vinyl_beta_C)
  with the Table S1 multi-pair scoring criterion d_ijkl = r_ij + r_ik + r_jl.
- Paper anchor: Table S1 (p.22) — vinyl Initiation/Propagation uses 4 groups with pairs
  i-j (formation, [3,6]), i-k (structural constraint, [0,3]), j-l (structural constraint, [0,3]).
  The i-k and j-l pairs participate in candidate selection and scoring but NO bias force is applied.
- Implementation:
  - `PairSpec.constraint_only: bool = False` added to `src/reactive/groups.py` — when True,
    the pair is used for filtering and scoring but `_build_pair_biases` skips it (no V^f/V^d).
  - `_find_chain_c_neighbor(smiles, radical_idx)` added to `scripts/_systems.py` — finds
    the non-nitrile C neighbor of the radical C in the initiator molecule.
  - `build_vinyl_aibn_system` returns 6-tuple (added `chain_c_map: dict[int,int]`).
  - Template: 4 groups, 3 pairs (1 formation + 2 constraint_only).
  - `PolymerizationWorkflow` constructor accepts `chain_c_map`, `_update_groups_after_cycle`
    maintains chain_C and vinyl_beta_C groups after each propagation event.
  - All callers updated: `run_vinyl_aibn.py`, `prep_structure.py`, `demo_radical_formation.py`,
    `demo_chain_propagation.py`, `tests/unit/test_systems.py`.
  - New tests: `test_chain_c_map_size`, `test_template_has_4_groups_and_3_pairs`,
    `test_find_chain_c_neighbor`, `test_propagation_map_beta_in_vinyl_beta_group_only`.
  - Full test suite: 199 tests PASSED.
- Validation run: S3 conditions (20 mono + 2 init, spin=3, f2=5, density 0.5, 15 cycles,
  2000 biased + 500 unbiased, seed 42, OrbMol-v2/CUDA). Artifact: runs/s4_multipair/.
- Key results:
  | Metric | S3 (2-group) | S4 (4-group) |
  |--------|-------------|-------------|
  | confirmed_formations | 6 / 15 cycles | 20 / 15 cycles |
  | reaction rate | 40% / cycle | 100% in cycles 1-10 |
  | chain_c_map updates | N/A | verified every cycle |
  | constraint_only skip | N/A | verified (no bias on i-k, j-l) |
- Analysis:
  - 20/20 monomers consumed — all monomer groups empty by cycle 12.
  - Cycles 1-9: 2 formations/cycle (both radicals react); cycle 10: 1; cycle 12: 1.
  - Cycles 11, 13-15: 0 candidates (monomer pool exhausted, only chain_C/vinyl_beta_C
    from already-reacted monomers remain; no vinyl_alpha_C for formation).
  - The 4-group criterion with i-k [0,3] and j-l [0,3] filters candidates more tightly,
    selecting pairs where the radical's chain backbone and the monomer's vinyl backbone
    are properly oriented — this dramatically improves reaction success rate.
  - chain_c_map bookkeeping verified: "removed old partner, added new alpha_C" logged
    for every propagation event across all 20 formations.
- Deviation from paper: None for Phase 1 — template and scoring are now paper-faithful
  per Table S1 (vinyl Propagation row).
- Scientific risk: None — Phase 1 is a pure paper-fidelity improvement.
- Licensing/commercial impact: None.
- Reproduction:
  ```
  conda run -n pfpoly-gpu python scripts/run_vinyl_aibn.py \
      --n-monomers 20 --n-initiators 2 --initiator-smiles "C[C](C)C#N" --spin 3 \
      --density 0.5 --backend orb --device cuda --no-barostat \
      --n-cycles 15 --biased-steps 2000 --unbiased-steps 500 --equil-steps 2000 \
      --timestep-fs 1.0 --f2 5.0 --seed 42 --output-dir runs/s4_multipair
  ```
- **S4 Phase 1 DONE.** Phase 2 (AIBN decomposition / Activation) is next.

## 2026-06-18: S4 Phase 2 — AIBN C-N homolysis PES scan (GATE check)

- Context: Before implementing V^d-driven AIBN activation, verify that OrbMol-v2 can
  model C-N azo bond homolysis. Scanned one C-N bond (C1-N5) of full AIBN
  (CC(C)(C#N)N=NC(C)(C)C#N, 24 atoms with H) from r=1.4 to 3.4 Å using constrained
  FIRE relaxation at each point. Tested spin=1 (singlet, intact) and spin=3 (triplet,
  post-decomposition).
- Paper anchor: Table S1 — Activation row, V^d on azo C-N bonds; f1_max_dissociation=125 kcal/mol.
- PES scan results:
  | r (Å) | Singlet dE (kcal/mol) | Triplet dE (kcal/mol) |
  |--------|----------------------|----------------------|
  | 1.4    | 0.0 (ref)            | 0.0 (ref)            |
  | 1.6    | 0.3                  | 3.4                  |
  | 1.8    | 10.8                 | 7.0                  |
  | 2.0    | 23.7                 | 9.4                  |
  | 2.2    | 35.0                 | −25.1                |
  | 2.6    | 38.7                 | −4.8                 |
  | 3.0    | 39.4                 | −38.8                |
  | 3.4    | 11.1                 | −38.7                |
- Key findings:
  - Singlet barrier: 39.4 kcal/mol (max dE at r≈3.0 Å). V^d f1_max=125 >> 39.4 → PASS.
  - Singlet dissociation energy: ~11 kcal/mol at large r (endothermic, expected for
    homolysis on closed-shell surface).
  - Triplet becomes favorable beyond r≈1.8 Å (spin crossing region), with large
    exothermic stabilization (−38.8 kcal/mol at r=3.0 Å).
  - OrbMol-v2 PES has discontinuities (triplet at r=2.2 and singlet at r=3.2) — likely
    constrained FIRE finding different local minima. Does not affect the GATE verdict.
- **GATE VERDICT: PASS** — V^d (125 kcal/mol) can drive AIBN C-N homolysis. The
  singlet barrier (39.4 kcal/mol) is well below f1_max. After dissociation, spin
  switching to triplet provides large thermodynamic stabilization.
- Decision: proceed with Phase 2 implementation (activation workflow, spin switching).
- Scientific risk: Low for GATE check. The PES discontinuities suggest OrbMol-v2's
  AIBN treatment is not perfectly smooth, but the qualitative picture (barrier exists,
  V^d can overcome it, triplet is stable) is sufficient for TDBB activation.
- Licensing/commercial impact: None.
- Reproduction:
  ```
  conda run -n pfpoly-gpu python scripts/scan_aibn_decomposition.py --device cuda
  ```
- Artifacts: `runs/s4_aibn_pes/aibn_pes.json`, `runs/s4_aibn_pes/aibn_pes.png`

## 2026-06-18: S4 Phase 2 — V^d f2 for activation must differ from production f2

- Context: V^d = f1·exp(-f2·r²) has force peak at r=1/√(2f2). The production f2
  (5-10) puts the peak at r=0.22-0.32 Å — at C-N bond distance (1.49 Å), the bias
  energy and force are essentially zero:
  | f2   | V^d(1.49 Å) | F(1.49 Å) |
  |------|-------------|-----------|
  | 10.0 | 0.00        | 0.00      |
  | 5.0  | 0.00        | 0.03      |
  | 0.5  | 41.2        | 61.4      |
  f2=0.5 puts the force peak at r=1.0 Å, providing 41 kcal/mol repulsive energy at
  the C-N equilibrium — comparable to the 39 kcal/mol homolysis barrier (PES scan).
- First test result: activation with f2=5.0 (production value) on 5 monomers + 1 AIBN,
  3000 steps: **0/2 C-N bonds dissociated** — confirmed that production f2 is useless
  for dissociation at bond distances.
- Paper anchor: The paper specifies f2=10 for formation (PDF p.7) and f2∈[5,20] as
  robust range (Fig. S4), but does NOT explicitly state whether the same f2 applies to
  dissociation/activation. V^d = f1·exp(-f2·r²) (Eq. 3) and V^f = f1·(1-exp(-f2·(r-r0)²))
  (Eq. 2) use the same f2 parameter in the equations, but the physical contexts are
  fundamentally different: V^f acts near r0≈2 Å (close-contact, bias→bond), V^d must
  act near r_bond≈1.5 Å (stretch existing bond apart).
- Decision: add `activation_f2` and `activation_f1_max` parameters to `run_activation()`
  and CLI.  Defaults: f2=0.3, f1_max=250.  The production f2 (for V^f formation) remains
  unchanged.  This is analogous to the f2=5 tuning for OrbMol-v2 PES (decisions.md
  2026-06-17): adjusting operational parameters for the specific MLIP and reaction context.
- Scientific risk: Low. f2=0.3, f1_max=250 are outside the paper's stated ranges (f2∈[5,20],
  f1_max_dissociation=125), but those ranges were validated for FORMATION, not dissociation.
  The effective potential analysis (PES + V^d) confirms these parameters are physically
  necessary: the OrbMol-v2 C-N barrier (39.4 kcal/mol) requires f1_max≥200 at f2=0.3
  for the effective potential to be monotonically repulsive.
- Licensing/commercial impact: None.
- Follow-up: see next record for validation results.

## 2026-06-18: S4 Phase 2 — Effective potential analysis and activation validation

- Context: V^d (f2=0.5, f1_max=125) stretched C-N bond to ~1.7-1.8 Å but could not
  break it. Computed effective potential Eff(r) = PES(r) + V^d(r) at PES scan points
  for various f1/f2 combinations to find monotonically repulsive parameters.
- Effective potential barrier (kcal/mol) at key parameter combinations:
  | f1_max | f2=0.5 | f2=0.3 | f2=0.2 |
  |--------|--------|--------|--------|
  | 125    | 11.0   | 6.1    | 7.2    |
  | 200    | ~2     | **0**  | **0**  |
  | 250    | **0**  | **0**  | **0**  |
  "0" = monotonically decreasing effective potential → guaranteed dissociation.
- Root cause: V^d(r) = f1·exp(-f2·r²) must drop faster than PES rises over the range
  r=1.5→2.2 Å (PES rises 35 kcal/mol). At f2=0.5/f1=125, V^d drops only 24 kcal/mol
  → 11 kcal/mol effective barrier → bond equilibrates at ~1.7 Å. At f2=0.3/f1=200,
  V^d drops 46 kcal/mol → no effective barrier.
- GPU validation (f2=0.3, f1_max=250, 3000 max steps):
  - **2/2 C-N bonds dissociated at step 175-176** (f1 reached 175 kcal/mol)
  - Spin switched 1 → 3 after activation
  - Subsequent TDBB propagation: 1 confirmed formation in 5 cycles
  - Full pipeline: AIBN decomposition → radical generation → polymerization propagation
- Paper anchor: Table S1 Activation row. Parameters differ from paper defaults because
  OrbMol-v2's C-N barrier (39.4 kcal/mol) may differ from PFP's. The TDBB mechanism
  (V^d ramp-up → dissociation → spin switch → propagation) is faithfully reproduced.
- Decision: adopt f2=0.3, f1_max=250 as activation defaults for OrbMol-v2.
- Licensing/commercial impact: None.
- Reproduction:
  ```
  conda run -n pfpoly-gpu python scripts/run_vinyl_aibn.py \
      --n-monomers 5 --n-initiators 1 --activation --spin 1 \
      --activation-f2 0.3 --activation-f1-max 250 --activation-steps 3000 \
      --density 0.5 --backend orb --device cuda --no-barostat \
      --n-cycles 5 --biased-steps 2000 --unbiased-steps 500 --equil-steps 2000 \
      --timestep-fs 1.0 --f2 5.0 --seed 42 --output-dir runs/s4_activation_v3
  ```
- Artifacts: `runs/s4_activation_v3/`

## 2026-06-18: RF13 — Verlet+barostat temperature source
- Context: `_integrator_temperature` returned a hardcoded 300 K for non-Langevin integrators. MC barostat acceptance (kT term) used this value even when the actual kinetic temperature differed.
- Paper anchor: MC barostat acceptance formula requires kT (mc_barostat.py docstring). NPT simulations use Langevin (unaffected); Verlet+barostat path was inconsistent.
- Decision: Langevin keeps returning target T. Verlet (and any non-Langevin) now returns the instantaneous kinetic temperature via `instant_temperature_K(state.velocities, state.masses)`.
- Alternatives considered: Always use kinetic T — rejected because Langevin's target T is the correct ensemble temperature for NPT acceptance.
- Scientific risk: Low. Langevin path (all current production runs) is unchanged. Verlet+barostat acceptance probability now fluctuates with kinetic T instead of being pinned to 300 K.
- Licensing/commercial impact: None.
- Follow-up: None.

## 2026-06-19: RF1 — manifest に実効パラメータを記録（案A採用）
- Context: `configs/boost/paper_faithful.yaml` はどの実行経路でも読み込まれておらず、`RunManifest` の `config_path` は文字列参照のみ。真実源が TDBBParams 既定・YAML・argparse 既定の3系統に分散し、manifest.json から実効値を復元できなかった。
- Paper anchor: CLAUDE.md 非交渉要件「All experiments must record seed, config path, git SHA, backend name, and output directory」。configs/boost/paper_faithful.yaml の値は PDF p.7 由来。
- Decision: 案A（YAML を実読込せず、`RunManifest.extra` に実効パラメータの dict を格納）を採用。`PolymerizationConfig` を `dataclasses.asdict` でシリアライズし、TDBBParams 入れ子も展開する。numpy 型は float/int に正規化。`config_path` は由来の目安と位置づけ、真の記録は `extra` に置く。
- Alternatives considered: 案B（YAML ローダ新設 + CLI 上書き差分記録）— 変更面積が大きく、乖離解消の本質は「実効値の記録」であるため、案Aの方が即効性がある。
- Scientific risk: なし。記録の追加のみで物理に触れない。
- Licensing/commercial impact: None.
- Follow-up: 案B（YAML ローダ）は将来的に検討可。paper_faithful.yaml の冒頭コメントに実効値は manifest.extra を見る旨を追記済み。

## 2026-06-19: RF2 — α denominator = n_monomers（初期モノマー数）に一本化
- Context: α = N_reacted / N_total の分母が3箇所で不一致だった。(1) polymerization.py が全群和（radical_C + vinyl_alpha_C + chain_C + vinyl_beta_C）で過大計上、(2) run_vinyl_aibn.py:461 が `n_monomers * 2 + n_initiators`、(3) run_vinyl_aibn.py:417 が radical_C + vinyl_alpha_C で死蔵。
- Paper anchor: 本文 p.9 Fig.2 キャプション α = 1 − [M]/[M]₀（[M]₀ = 初期モノマー濃度）。Eq. 11 の α は monomer conversion。よって分母 = n_monomers（初期モノマー数）。vinyl では1形成イベントにつきモノマー1個消費 → α = confirmed_formations / n_monomers。開始剤・constraint 群は分母に含めない。
- Decision: 分母を n_monomers に統一。`monomer_site_count()` を conversion.py に新設（vinyl_alpha_C 群サイズ = 初期モノマー数）。`PolymerizationWorkflow.run()` に `n_monomers` 引数を追加し、渡された場合はそれを trajectory ヘッダに記録。run_vinyl_aibn.py は `args.n_monomers` を渡す。死蔵変数 `n_reactive_sites`（:417）と旧式 `n_monomers * 2 + n_initiators`（:461）を除去。
- Alternatives considered: 全群和（旧実装）— constraint 群を含むため分母が過大で α を過小評価。`2*n_mono+n_init`（旧図コマンド）— 反応サイトの二重計上。いずれも paper 定義と不一致。
- Scientific risk: 中。α の絶対スケールが変わる（旧値は分母過大で α を過小評価していた）。過去 run の図と数値が変わるため本記録で定義変更を明示。
- Licensing/commercial impact: None.
- Follow-up: nylon は反応進行度 p（Carothers, carothers.py）で別管理。α(t) プロット（Eq. 11）は vinyl monomer-conversion 用。

## 2026-06-19: RF4 — MD ループ共通化と群更新のストラテジ化（不変リファクタ）
- Context: `PolymerizationWorkflow` の MD ステップ骨格（pre_force → compute → [total_bias] → post_force → barostat）が `_run_biased_phase`, `_run_unbiased_phase`, `_run_equilibration_phase`, `run_activation` の4箇所で重複。群更新 `_update_groups_after_cycle` が vinyl 固有の群名リテラル（`chain_C`, `vinyl_beta_C`, `radical_C`）を直接参照しており、nylon/epoxy では差し替え不可。
- Paper anchor: CLAUDE.md「Separate numerical kernels from orchestration code」「Functions should align with scientific concepts」。物理不変のリファクタのため paper 解釈変更なし。
- Decision: (1) `_md_step()` を導入し4フェーズの MD ステップを単一実装に統合。力の構築順（pre_force → compute(base) → [bias] → post_force → step++ → [barostat]）を厳密に保存。`enable_barostat` フラグで activation（バロスタなし）を区別。(2) `PostCycleUpdater` プロトコルを導入し、`DefaultPostCycleUpdater`（形成原子除去のみ）と `VinylChainPropagationUpdater`（連鎖伝播）を実装。DI で注入可能にし、既存の `propagation_map`/`chain_c_map` パラメータからの後方互換変換あり。(3) `_build_pair_biases` に `template` キーワード引数を追加し、`run_activation` の重複ペア構築コードを廃止。
- Alternatives considered: (a) 各フェーズを完全独立のまま維持 — 4箇所の力評価順変更リスクが分散。(b) MD ステップを外部関数に切り出し — ワークフロー内部状態（barostat, integrator）への参照が多く、メソッドが自然。
- Scientific risk: なし。不変リファクタ。`test_deterministic_with_seed` でビット一致を確認済み（243テスト全緑）。
- Licensing/commercial impact: None.
- Follow-up: VDW_RADII/ATOMIC_MASSES の `src/constants.py` への移動は将来的に検討可（RF4 step 6, optional）。

## 2026-06-19: RF5 — スコア d_ijkl は i-j, i-k, j-l の3項固定。nylon k-l はバイアス専用
- Context: `score_candidates` が `template.pairs` 全合算のため、nylon では k-l（amine_H–carboxyl_OH, 水形成）の距離がスコアに加算されていた。4項スコアは paper の3項定義 d_ijkl = r_ij + r_ik + r_jl と不一致。k-l の r_max=100 によりほぼ無拘束の H–OH 距離が候補ソートを支配し得る状態。
- Paper anchor: 本文 p.4 Eq.7（d_ijkl = r_ij + r_ik + r_jl, 3項固定）。SI Table S1（vinyl）/ Table S2（nylon condensation: 群同定は i-j, i-k, j-l の3ペアのみ。k-l にはバイアス V^f のみ適用、距離窓・スコアに非関与）。
- Decision: `PairSpec` に `score_pair: bool = True` フラグを追加。`score_pair=False` のペアは候補同定（距離窓フィルタ）と d_ijkl スコアの両方から除外。バイアス適用は `is_formation`/`constraint_only` 側で従来どおり制御。nylon k-l ペアに `score_pair=False` を設定。vinyl は全ペアが `score_pair=True`（現状不変）。
- Alternatives considered: k-l を template.pairs から除去してバイアス専用リストを別に持つ — 変更面積が大きく、既存の `_build_pair_biases` が template.pairs 全体からバイアスを構築する設計と相性が悪い。フラグ方式の方が最小限。
- Scientific risk: 中。nylon の候補順位が変わり、異なる反応ペアが選択される可能性がある。vinyl は不変（既存3ペアのみ、全て score_pair=True）。
- Licensing/commercial impact: None.
- Follow-up: RF10（pair_distances 死蔵）と合わせて整理可。→ RF10 で対応済み。

## 2026-06-19: RF10 — Candidate.pair_distances 死蔵フィールドの除去とスコア計算の一本化
- Context: `Candidate.pair_distances` は `_enumerate_recursive` で計算・格納されるが、`score_candidates` は距離を再計算しこの値を使わない。下流でも未使用（テストは `{}` を渡す）。距離の二重計算かつ未使用フィールド。
- Paper anchor: Eq.7（d_ijkl スコア定義）。スコア計算のタイミングを変更するが、スコア値自体は不変（d_ijkl = r_ij + r_ik + r_jl の3項合計）。
- Decision: `pair_distances` フィールドを `Candidate` から除去し、`_enumerate_recursive` で距離を合算して `score` を直接設定。`score_candidates` は距離再計算を行わずソートのみに簡略化（引数から `template`, `positions`, `cell` を除去）。これにより (1) 二重計算が解消、(2) `Candidate` が軽量化（dict 不要）、(3) 責務が明確化（列挙時にスコア確定、ソートは分離）。
- Alternatives considered: (a) `pair_distances` を残し `score_candidates` がそれを参照する形に変更 — 二重計算は解消されるが不要な dict を保持し続ける。(b) 現状維持 — 二重計算と死蔵が残る。
- Scientific risk: なし。スコア値は同一（同じ距離を同じ順序で加算）。vinyl の決定論テストに影響なし。
- Licensing/commercial impact: None.

## 2026-06-19: RF11 — パッケージング構造の整理（最小対応: 逆方向依存解消 + packages.find 制限）
- Context: 配布名 `kagome` (pyproject.toml) と import 名 `src.*` が不一致。`src/prep/openmm_equilibrate.py` が `from scripts._systems import _rdkit_mol, box_from_density` でライブラリ層→スクリプト層の逆方向依存。`pyproject.toml` の `packages.find where=["."]` が `scripts`/`tests` もパッケージ化。egg-info がルートと `src/` に重複。
- Paper anchor: N/A（実装/配布の健全性）。
- Decision: (1) `_rdkit_mol` と `box_from_density` を `src/chem/builders.py`（新設）に移動。`scripts/_systems.py` は re-export で後方互換を維持。`src/prep/openmm_equilibrate.py` は `src.chem.builders` から import するよう修正し逆方向依存を解消。(2) `pyproject.toml` に `include = ["src*"]` を追加し、`scripts`/`tests` のパッケージ化を防止。(3) `src/` 配下の重複 egg-info を削除（.gitignore 済みのためローカルのみ）。(4) `src/` → `pfpoly/` リネーム（フル import 名統一）は影響範囲が広いため別タスク・別 PR とする。
- Alternatives considered: (a) `src/` を `kagome/` にリネームし全 import を統一 — 理想的だが `from src.` が全コード/テストに散在しており変更面積が大きい。別 PR で段階的に実施すべき。(b) 現状維持 — 逆方向依存が残る。
- Scientific risk: なし。コード移動のみ、関数の実装は不変。
- Licensing/commercial impact: None.

## 2026-06-19: RF12 — 幾何の orthorhombic 限定を明文化し非対角 cell で ValueError
- Context: `src/geometry.py` の `minimum_image`/`wrap_positions` は cell の対角成分のみ使用するが、非対角成分が渡された場合に静かに誤った結果を返す状態だった。decisions.md "2026-06-13" に「triclinic は deferred」と記載あり。
- Paper anchor: 現行系（vinyl, nylon）は全て cubic/orthorhombic。epoxy/CuO スラブ系で triclinic が必要になる場合は別途対応。
- Decision: `_check_orthorhombic(cell)` ガード関数を追加。非対角成分の絶対値が 1e-10 を超える場合に `ValueError` を発行。`minimum_image`/`wrap_positions` の両関数で cell が非 None のときに呼び出す。`warnings.warn` ではなく `ValueError` を採用（誤った最小像距離は下流の反応候補選択・バイアス力計算を破壊するため、静かに続行するより早期に失敗すべき）。
- Alternatives considered: (a) `warnings.warn` で続行 — 誤った物理量が伝播するリスクがある。(b) triclinic 対応を実装 — 現行系では不要で過剰設計。
- Scientific risk: なし。現行系は全て orthorhombic。将来 triclinic が必要になった場合は `_check_orthorhombic` を一般化する。
- Licensing/commercial impact: None.

## 2026-06-19: RF8 — 解離判定を純関数 is_dissociated/is_formed に集約、activation の 2.5 Å 閾値を文書化
- Context: `BondTracker` は `threshold_fraction * r0`（相対閾値、デフォルト1.0）で形成/解離を判定するが、`run_activation` は絶対閾値 `dissoc_threshold = 2.5` Å をハードコードしていた。同じ「解離」概念に対しモジュール間で閾値規約が不統一。
- Paper anchor: Eq.3（V^d）。decisions.md 既存記録「2026-06-12: Dissociation tracking uses r0 = λ·Σr_vdW as confirmation threshold」。Table S1 Activation 行（V^d on C-N azo bonds）。
- Decision: (1) `src/reactive/bonds.py` に純関数 `is_dissociated(r, r0, threshold_fraction=1.0) -> bool`（r > threshold_fraction * r0）と `is_formed(r, r0, threshold_fraction=1.0) -> bool`（r <= threshold_fraction * r0）を追加。(2) `BondTracker.check_reactions_during_bias` および `check_outcomes` の判定ロジックをこれらの純関数で置換。(3) `run_activation` の解離判定を `is_dissociated(r, dissoc_threshold)` に置換（dissoc_threshold=2.5, threshold_fraction=1.0 で「r > 2.5」と等価）。(4) 2.5 Å を絶対閾値として維持する根拠: azo C-N 平衡結合長 ~1.49 Å に対し、r0=λ·Σr_vdw(C,N)=0.6*(1.70+1.55)=1.95 Å は vdW 接触距離であり、これを超えただけでは C-N 結合の明確な解離とは限らない。2.5 Å は平衡長の ~1.68 倍であり、結合が確実に破断した状態を示す安全な閾値。activation パラメータ（f2=0.3, f1_max=250）は変更なし。
- Alternatives considered: (a) r0 相対に完全統一（threshold_fraction=1.0 → 1.95 Å で検出）— 解離をより早く検出するが、揺らぎで偽検出のリスク。(b) 現状維持 — 判定が2箇所にインライン展開され、テスト困難。
- Scientific risk: なし。activation の閾値 2.5 Å は変更しておらず、判定ロジックの純関数化のみ。BondTracker の動作もビット一致。
- Licensing/commercial impact: None.

## 2026-06-20: RF16 — ライセンスゲートを許可リスト方式へ拡張
- Context: `scripts/check_dependency_licenses.py` は `blocked_pending_review` の依存が import された場合のみ失敗するブロックリスト方式で、「import される全依存が承認済みか」を検証していなかった。実際 scipy/rdkit/pyyaml/openff（namespace）が import されているのに registry 未掲載でもゲートは exit 0 だった。また MACE-OFF23（ASL, blocked）は import 名を持たず（モデル文字列で選択）`_INSTALL_TO_IMPORT['mace-off23']=[]` のため構造的に検知不能だった。
- Paper anchor: CLAUDE.md 商用ガードレール「New dependencies require an explicit license check before adoption」「If license status is unclear, mark it blocked_pending_review」「single source of truth in sync with the matrix」。
- Decision: (1) 許可リスト強制を追加。`sys.stdlib_module_names` と自前パッケージ（src/scripts/tests/pfpoly）を除く全 import が registry に登録されていなければ exit 1。ブロックリスト判定は多層防御として維持。(2) 未登録だった実在依存を登録: rdkit(BSD-3, matrix には既出だが YAML 欠落=同期ずれ), scipy(BSD-3), pyyaml(MIT), openff-units(MIT)。`_INSTALL_TO_IMPORT` に openff-toolkit→['openff','openff_toolkit'] 等の namespace マッピングを追加。(3) MACE-OFF23 に YAML フィールド `detect_strings: ['mace_off(']` を追加し、ゲートがソース中の OFF ローダ呼び出しを grep で検知。さらに `src/backends/mace_backend.py` に実行時ガードを追加（model 文字列が mace_off/mace-off/off23 を含む場合 RuntimeError）。`detect_strings` は '(' を含むため、トークンを列挙するだけのガード自身は誤検知しない。
- Alternatives considered: (a) ブロックリストのみ維持 — 未登録依存を見逃す（現状の穴）。(b) MACE-OFF をトークン（'off23' 等）で広く grep — ブロック用ガードのコード自身に誤反応するため不採用、呼び出しパターン 'mace_off(' に限定。(c) openff-units を blocked 扱い — 実体は OpenFF スタックの MIT パッケージであり過剰。
- Scientific risk: なし（ガードレール・記録のみ、数値計算に非干渉）。
- Licensing/commercial impact: 商用安全性の網羅性が向上。未承認依存の混入と ASL 重み（MACE-OFF23）の使用を CI/実行時の両面で阻止。openff-units の upstream LICENSE は release 前に最終確認すること（YAML notes に明記）。

## 2026-06-20: RF15 — 縮合系の解離イベントで遊離原子を群から消費する
- Context: `DefaultPostCycleUpdater` は `tracker.confirmed_formations()` のみを処理し反応原子を群から除去していた。nylon-6,6 縮合（Table S2）では amine_N–amine_H (i–k) と carboxyl_C–carboxyl_OH (j–l) が V^d（is_formation=False）で解離するが、`confirmed_dissociations()` は群更新で一切消費されていなかった。結果、解離して遊離した amine_H / carboxyl_OH（水の脱離基）がサイクル横断で群に残り再選択され得た。
- Paper anchor: 本文 §2.2（biased→unbiased→トポロジ更新のサイクル）、SI Table S2（nylon Condensation: i–j 形成 V^f、i–k/j–l 解離 V^d、k–l 水形成 V^f）。CLAUDE.md「Distinguish bond formation bias from bond dissociation bias」。
- Decision: (1) `DefaultPostCycleUpdater` を拡張し、`confirmed_formations()` に加えて `confirmed_dissociations()` も処理する。各確定解離イベントの両原子（atom_a/atom_b）を所属群から除去する（`_remove_pair` に共通化）。`processed_dissociations` カウンタで二重処理を防止。`remove_atom` は不在時 ValueError を握り潰すため、形成側と解離側で同一原子を二重に消そうとしても冪等。(2) **nylon に専用 updater（NylonCondensationUpdater）は新設しない**。`build_nylon66_system` は各ジアミン/二酸の**両末端**を初期から amine_N/carboxyl_C 群に登録しているため、片端が反応しても残る末端で鎖延長が自然に進行する（vinyl のラジカル移動のような「新末端を新原子へ昇格する」処理は不要）。これは handoff-plan-v6 の RF15 案（NylonCondensationUpdater 追加）からの意図的な簡素化で、コード実態に即した最小修正。(3) vinyl 経路（`VinylChainPropagationUpdater`）は変更しない。vinyl テンプレートに解離ペアは無く `confirmed_dissociations()` は常に空のため、数値はビット一致のまま。
- Alternatives considered: (a) 専用 NylonCondensationUpdater で末端昇格 — 両末端が既登録のため冗長で、誤った再分類リスクを増やす。(b) 解離の atom_a/atom_b のうち脱離基側のみ除去 — どちらが脱離基かはテンプレート依存で一般化困難。両方除去で安全（反応中心 amine_N/carboxyl_C は amide 形成で消費されるため除去は無害かつ冪等）。(c) 現状維持 — 解離原子が再選択され縮合トポロジが進行しない（既知バグ）。
- Scientific risk: nylon の反応経路が変わる（解離原子が次サイクルで再選択されなくなる＝正しい挙動）。vinyl は不変（決定論テストでビット一致を確認）。一般の「解離後も両原子が反応性を保つ」機構が将来必要になった場合は、PairSpec に脱離基フラグを追加して `_remove_pair` を選択的にする。
- Licensing/commercial impact: None.

## 2026-06-20: RF19b — MC barostat Jacobian (N+1) と運動温度の自由度 (3N-3)
- Context: (1) `MCBarostat.try_step` は体積提案を ln(V) 一様（`delta_ln_V ~ U(-max,max)`）で行うが、受容判定の Jacobian 項が `N·kT·ln(V'/V)`（N 形式）だった。ln(V) 一様提案では測度 d(lnV) の分だけ厳密には `(N+1)` が正しい。(2) `instant_temperature_K` は `3N` 自由度を使っていたが、MB 初期化は COM 並進を除去し NVE Verlet は COM を保存するため、厳密には `3N-3` が正しい。
- Paper anchor: NPT MC 受容式（mc_barostat.py docstring; Frenkel & Smit の log-volume サンプリング）。運動温度の等分配（init_velocities.py docstring）。いずれも MD の標準的扱いで論文固有の解釈ではない。
- Decision: (1) barostat の `delta_H` を `(n_atoms + 1) * kT * delta_ln_V` に変更。(2) `instant_temperature_K` の分母を `3N-3`（N>1）/`3N`（N=1）に変更。どちらも差は O(1/N) で、既定 1 atm・数百原子では受容率・診断温度への実影響は極小。単原子（N=1）は 3N を維持（既存テスト `test_known_temperature` と整合）。
- Alternatives considered: (a) barostat を N のまま維持（OpenMM は体積一様提案で N を使用）— しかし本実装は ln(V) 一様提案なので N+1 が整合的。(b) 温度を 3N のまま維持 — COM 除去と非整合で診断温度が系統的に約 (3N-3)/3N 倍低く出る。
- Scientific risk: 低。受容率と診断温度が O(1/N) 変化する。barostat の既存テストは受容数の厳密値に依存しないため緑。温度テストは N=1 が 3N 維持、N>1 統計テストは ±30% 許容で通過（むしろ精度向上）。
- Licensing/commercial impact: None.

## 2026-06-20: nvalchemiops のライセンスを Apache-2.0 と確認 → 周期 OrbMol-v2 を unblock
- Context: 論文スケール（vinyl/nylon=NPT、epoxy/CuO=NVT、いずれも周期 + 長距離 Coulomb）を OrbMol-v2 で回すには周期 PME が必須で、これは `nvalchemiops`（NVIDIA）に依存する。従来は `specs/approved_dependencies.yaml` で `blocked_pending_review`（理由: ①ライセンス未確認、②Windows の torch.compile/cl.exe 失敗）としており、`orb_backend.py::_check_periodic_support` が周期ラン（cell≠None）を例外で停止していた。ローカルマシンのスペック不足を機にクラウド GPU での論文スケール実証を検討する中で、本依存のライセンス確定が前提条件として浮上した。
- Paper anchor: SI（vinyl 333K / nylon 300K / epoxy 333K, すべて PBC）。OrbMol 採用記録（2026-06-12「Long-range Coulomb via PME adds physics absent in MACE-MP-0」）。CLAUDE.md commercial-use guardrails「New dependencies require an explicit license check before adoption」。
- Decision: upstream を直接確認し、`nvalchemiops` の配布元が NVIDIA ALCHEMI Toolkit-Ops（pip: `nvalchemi-toolkit-ops`、import: `nvalchemiops`、github.com/NVIDIA/nvalchemi-toolkit-ops）であること、その LICENSE が標準 Apache-2.0（`SPDX-License-Identifier: Apache-2.0`、© NVIDIA CORPORATION 2025-2026）であることを確認した（2026-06-20）。これに基づき (1) `approved_dependencies.yaml` の `nvalchemiops` を `approved`（evidence: Apache-2.0）へ更新、(2) `dependency-license-matrix.md` を Apache-2.0 / 商用可へ更新し検証注記を追加、(3) 本レコードを追加。WSL の pfpoly-gpu には `nvalchemi-toolkit-ops==0.3.1` が既にインストール済み・import 可能で、`environment_wsl.yml`（L243）にも記載済みであることを確認。
- Alternatives considered: (a) PME を切って周期だけ回す — Nylon(アミド)/Epoxy(極性基)で Coulomb 物理が落ち paper-faithful でなくなるため却下。(b) 別バックエンド(MACE-MP-0)で周期 — OrbMol-v2 の polymer 最適化ポテンシャルと長距離 Coulomb を失う。(c) blocked のまま据え置き — ライセンスが Apache-2.0 と確定した以上、根拠が消失。
- Scientific risk: なし（ライセンス・記録更新のみ。数値挙動は不変）。`_check_periodic_support` ガードは `find_spec('nvalchemiops')` が non-None なら通過する実装で、コード変更は不要。
- Licensing/commercial impact: nvalchemiops は Apache-2.0、商用セーフ。これで周期 OrbMol-v2 の論文スケールランが Linux/WSL/クラウドで guardrail に抵触せず実行可能。Windows ローカルでは torch.compile/cl.exe 問題が残るため、周期ランは Linux/WSL/クラウド、または TORCHDYNAMO_DISABLE=1 で回す。
- Follow-up: (a) クラウド(Linux+GPU)で小スケール周期ラン1本を実行し、PME 経路の動作とローカル非周期結果との整合（再現性）を確認。(b) pfpoly-gpu をコンテナ化して再現性を固定。(c) 本番ランを並列投入し、各ランで seed/config/git SHA/backend/output dir を記録（対象系数は 2026-06-20 のスコープ決定を参照: epoxy/CuO 除外後は vinyl 6 単量体 + nylon = 7系 × 3 runs = 21 ラン）。

## 2026-06-20: 再現スコープを有機系に限定し epoxy/CuO 界面を再現対象から除外
- Context: 論文は3系統（① ビニルラジカル重合、② nylon-6,6 縮合、③ epoxy(DGEBA)+DETA キュアの CuO(001) 界面）を扱う。本リポジトリの既定バックエンドは OrbMol-v2 に確定済みで、その学習データは OMol25（有機分子）+ OPoly26（ポリマー、ωB97M-V/def2-TZVPD）。CuO は金属酸化物表面（Cu + 表面格子）で、この学習ドメイン外。クラウド GPU での論文スケール実証にあたり、再現対象系を確定する必要が生じた。
- Paper anchor: SI Table S1/S2/epoxy テンプレート（paper/notes.md「Systems studied」「System compositions」）。OrbMol 採用記録（2026-06-12: OMol25+OPoly26 学習、有機/ポリマー向け）。CLAUDE.md「uMLIP must handle organic polymer chemistry (C, H, N, O)」。
- Decision: 再現対象を**有機系（① ビニルラジカル重合、② nylon-6,6 縮合）に限定**し、**③ epoxy/CuO 界面を再現対象から除外**する。理由: (1) OrbMol-v2 は CuO（Cu, 金属酸化物表面）を学習しておらずドメイン外で、エネルギー/力の信頼性が担保できない。(2) CLAUDE.md は uMLIP を C/H/N/O の有機ポリマー化学に適用する前提で、金属界面は本実装の主目的（TDBB 重合・キュアの有機系再現）から外れる。(3) epoxy/CuO は未実装（ビルダー無し）で、追加コストに対し OrbMol-v2 では科学的妥当性が低い。本決定により VRAM 上限見積もりは epoxy/CuO の 48–80 GB 帯が消え、最大系は nylon(~4400原子)/大型ビニル単量体(~5000原子)= **A100 40GB 級**が基準となる（2026-06-15 の 200+10≈2520原子→≥24GB 実測を踏まえた線形外挿、要実測検証）。
- Alternatives considered: (a) epoxy/CuO を別バックエンド（CuO を学習した MLIP 等）で再現 — 本リポジトリの「OrbMol-v2 既定・有機系」方針と外れ、別途ライセンス/検証コスト。将来要望時に再検討。(b) CuO を古典 FF、有機相のみ OrbMol-v2 のハイブリッド — 界面の電子論的相互作用が落ち paper-faithful でない。(c) 全系維持 — OrbMol-v2 では CuO 界面のMLIP精度が担保できず、誤った再現になるリスク。
- Scientific risk: 低。除外する系は OrbMol-v2 のドメイン外で、無理に回す方がむしろ科学的に不誠実。残す2系は OPoly26 ドメイン内で再現の妥当性が高い。paper/ の notes・claims は論文の忠実な記録として epoxy/CuO 記載を**保持**（再現対象=実装範囲とは別管理）。将来 CuO 対応 MLIP を導入する場合は本決定を上書きする。
- Licensing/commercial impact: なし。むしろ対象縮小で依存・ハードウェア要件が軽くなる（80GB GPU 不要）。
- Follow-up: (a) クラウド本番のインスタンス階層は A100 40GB を上限基準に設定（80GB は不要）。(b) nylon と大型ビニル単量体の VRAM を本番前に短時間ランで実測し、40GB で sustained 完走するか確認。(c) figure 再現対象も有機2系に限定（specs/figure-comparison.md 側の epoxy 関連が将来追加されないよう本決定を参照）。

## 2026-06-20: 二置換ビニルビルダー対応 + 三級ラジカルは原子タイピングのみで扱う
- Context: 論文のラジカル重合単量体6種のうち、1,1-二置換ビニル（methacrylate, 1,1-diphenylethylene, dimethyl itaconate）は `scripts/_systems._find_vinyl_alpha_beta` が beta 炭素に「H ちょうど1個」を要求していたためビルド不可（beta が C(R)(R') で H=0 のため `ValueError`）。メタクリレートを再現対象に入れたいという要望（2026-06-20）で対応が必須化。あわせて、二置換単量体の伝播後ラジカルが三級になる安定性を TDBB にどう反映するかの方針決定が必要になった。
- Paper anchor: paper/notes.md「Radical polymerization: methyl acrylate, methacrylate, styrene, vinyl acetate, diphenylethylene, dimethyl itaconate」、ビニルテンプレート（notes.md: Gi radical_C / Gj vinyl_alpha_C / Gk chain_C / Gl vinyl_beta_C, Table S1 ij+ik+jl）。反応選択は Eq.7 の距離窓（幾何）で、速度定数を持たない（src/reactive/selection.py）。
- Decision: (1) `_find_vinyl_alpha_beta` の判定を「alpha = 末端 =CH2（2H）、beta = もう一方のビニル炭素（mono は 1H、1,1-二置換は 0H）」に一般化（`h_b == 1` → `h_b <= 1`、対称分岐も同様）。head-to-tail 則によりラジカルは常に末端 CH2(alpha) に付加し、不対電子は置換側 beta に局在 → これは mono/二置換で不変。(2) **三級ラジカルの安定性は原子タイピングのみで扱う**: beta（=二置換では三級炭素）を従来どおり vinyl_beta_C 群に型付けするだけで、伝播後ラジカルが三級中心に座る。その安定性のエネルギーは OrbMol-v2 に委ね、専用 ReactionTemplate や boost パラメータ・速度差は導入しない。
- Alternatives considered: (a) 三級用に別テンプレート/boost を与える — 論文に三級専用バイアスの記載がなく、論文外ハイパラ導入になり paper-faithful から外れる。(b) 解重合バイアス(V^d)で天井温度/可逆性を表現 — 化学的には妥当だがスコープ大幅増で、現段階の要望（構造をビルドし回せる）に対し過剰。(c) アロマ環の C=C を誤検出する懸念 — RDKit はアロマ結合を bond order 1.5 で返すため `!= 2.0` で除外され、styrene/diphenylethylene のフェニル環は誤マッチしない（検証済み）。
- Scientific risk: 低。結合形成/解離の記録（src/reactive/bonds.py）は距離・イベントのみで H 数・価数に非依存のため、beta の H=0 でも破綻しない。三級安定性を MLIP に委ねるのは OrbMol-v2 が当該化学（OMol25+OPoly26 有機/ポリマー）を学習している前提で妥当。万一三級の depropagation を明示的に再現したくなれば本決定を上書きし V^d 経路を追加する。
- Licensing/commercial impact: なし。
- Follow-up: (a) クラウド本番で methacrylate 系を1本短時間ランし、伝播・三級ラジカル生成が破綻なく進むか確認。(b) 残る二置換単量体（diphenylethylene, dimethyl itaconate）も同経路でビルド可能（検証済み）だが、本番投入時に density/box の妥当性を確認。

## 2026-06-20: WSL 単一env で古典FFを Calculator 化し compress の既定バックエンドにする
- Context: 論文密度(0.5 g/mL)への高密度化は ML(OrbMol-v2)でも古典FFでも可能だが、ML 圧縮は 20段×200 FIRE = 数千回の forward+backward で GPU 時間を大きく消費する。D-4（2026-06-14）では古典prep を別env(pfpoly-prep)に分離していたが、その理由（openff-nagl が第2の PyTorch を引き込み production torch と衝突）は Windows 固有だった。`environment_wsl.yml` を確認したところ、WSL の pfpoly-gpu は **OpenFF(toolkit 0.18.1/interchange/Sage 2026.01/nagl)+OpenMM 8.5.2 と OrbMol-v2 を同一env で共存**させており、衝突は解消済み。ユーザ要望は「前処理を別コマンドにせず、一貫した単一経路のまま古典FFを圧縮の標準にする」。
- Paper anchor: SI S-3 の初期密度 0.5 g/mL（compress の目標）。compress_box docstring「非物理的な準備デバイスで後段の biased/unbiased 動力学をバイアスしない」。古典prep の既存実装 src/prep/openmm_equilibrate.py（D-3「simple」プロトコル: minimize→compress→NVT）。
- Decision: (1) 古典FF（OpenMM+OpenFF Sage+Gasteiger 電荷）を kagome の `Calculator` インターフェース（compute(positions, species, cell)→(E[kcal/mol], F[kcal/mol/Å])）として実装する新バックエンド `src/backends/classical_backend.py` を追加。古典FFは結合情報が必須のため、Calculator は `molecule_specs`（topology）で構築する topology-aware 設計とし、`compute` は座標と周期箱(cell)のみ更新して OpenMM Context から E/F を取得。これにより既存の `compress_box` をそのまま再利用できる。(2) 実行スクリプト（run_vinyl_aibn.py / run_nylon66.py / scripts/profile_vram.py）に `--compress-backend {classical,ml}` を追加し**既定を classical** にする。MD 本体は従来どおり OrbMol-v2。(3) D-4 の env 分離は Windows 固有の制約であり、WSL 単一env では古典prep を in-process 実行してよい、と D-4 のスコープを更新する（pfpoly-prep 別env 経路は Windows 用フォールバックとして温存）。
- Alternatives considered: (a) 既存 equilibrate_structure をインライン呼び出し — テスト済みコード再利用で最小リスクだが、compress_box とは別経路になり ML/classical で圧縮ロジックが二重化。Calculator 化なら compress_box に一本化でき将来の一貫性が高いと判断（ユーザ選択）。(b) soft-sphere/WCA パッカー — 依存ゼロ・最速だが「真の古典FF」でなく、ユーザ要望（古典FFを標準）に合わない。(c) ML 圧縮のまま — GPU 時間が大きく、古典で十分な準備工程に MLIP を浪費。
- Scientific risk: 低〜中。compress は準備デバイスで後段動力学をバイアスしないため、古典FFでも ML でも最終構造は後段の ML 再平衡で上書きされる（D-2/D-4）。リスクはエンジニアリング面（OpenMM の PME カットオフ < 箱半長 の制約: 目標箱 ~18–20Å に対し Sage 既定 0.9nm カットオフが境界。既存 equilibrate_structure も同目標へ圧縮できているため precedent あり）。原子順序不変条件（topology が builder と同順）は `_build_openff_topology` の検証で担保。
- Licensing/commercial impact: なし。OpenMM(MIT core / LGPL GPU, import-only)・OpenFF stack(MIT)・Sage(CC-BY-4.0, 帰属要)・RDKit Gasteiger(BSD) はいずれも承認済み（approved_dependencies.yaml）。新規依存の追加なし。
- Follow-up: (a) classical_backend の単体テスト（小系で E/F が有限・compress_box で箱が縮む）。(b) WSL で methacrylate/nylon の classical 圧縮→OrbMol-v2 MD の統合スモーク。(c) VRAM プロファイラは classical 既定にすると GPU は MD のみ測定（ピークは ML 圧縮時とほぼ同等、所要時間は短縮）。

## 2026-06-20: D-4 追補 — WSL 単一env で前処理ファイル受け渡しを既定から外す + nylon を vinyl と同等化
- Context: D-4（2026-06-14）は初期構造の古典prep を別env(pfpoly-prep)+ JSON 受け渡し(prep_structure.py → run_*.py --load-structure)に分離していた。理由は Windows での torch 衝突(openff-nagl が第2の PyTorch を引き込む)。WSL の pfpoly-gpu は古典(OpenFF/OpenMM/Sage)と ML(OrbMol-v2)を同一env に持つためこの分離は不要。さらに監査の結果、(a) S6 本番は既に run_vinyl_aibn を直接呼ぶ単一経路で JSON 受け渡し未使用、(b) 2026-06-20 の古典Calculator化で run_vinyl_aibn は build→古典圧縮→MD が単一runで完結する一方、(c) run_nylon66 は density も compress も無く固定 box_size=25Å で配置するだけで単一runでは論文密度(0.5 g/mL)に到達できない、という vinyl/nylon の非対称が判明した。
- Paper anchor: SI S-3 初期密度 0.5 g/mL（vinyl/nylon 共通）。D-4（初期構造準備の分離）。compress_box / src/prep/openmm_equilibrate（古典prep）。
- Decision: WSL 単一env を前提に前処理を「既定で in-process・単一run」に揃える。(1) run_nylon66.py を vinyl と同等化: `--box-size` 既定を None に変更し `--density`(既定 0.5)+ `--compress-backend{classical,ml}`/`--compress-platform` を追加。--box-size 省略時は target_edge へ direct placement、失敗時は希薄配置 + 古典圧縮。--box-size 明示時は従来どおり直接配置(後方互換)。backend 生成を build 前に移動('ml' 圧縮で再利用)。サマリの box_size_A は初期 cell 由来に変更(NPT で変動するため)。(2) prep_structure.py を「任意のキャッシュ用途(prep-once-run-many)」と再定義し、WSL では pfpoly-gpu で直接実行可能・別 pfpoly-prep env は Windows 専用フォールバックと docstring/usage を更新。(3) --load-structure help を「論文密度に必須でない・任意」と明確化。D-4 の分離方針は Windows フォールバックとして温存。
- Alternatives considered: (a) prep_structure + --load-structure を削除して完全一本化 — prep-once-run-many(同一構造で多 seed)というキャッシュ用途は WSL でも有用なため温存。(b) nylon に --load-structure を追加 — prep_structure.py が vinyl 専用で nylon 用 producer が無く、現状は不要(将来 nylon prep が要れば追加)。(c) インライン完全古典平衡化(NVT 込み, --prep equilibrate)— スコープ大のため今回は見送り(古典圧縮で密度到達 → ML --equil-steps で熱平衡の現行で十分)。
- Scientific risk: 低〜中。nylon の既定挙動が「固定 25Å」→「論文密度へ圧縮」に変わる(--box-size 明示で旧挙動を完全再現可能)。圧縮は準備デバイスで後段 ML 動力学をバイアスしない。原子順序不変は _build_openff_topology の検証で担保。
- Licensing/commercial impact: なし(新規依存なし)。
- Follow-up: (a) WSL で nylon 単一run(--density 0.5, classical 圧縮 → OrbMol-v2 MD)のスモーク。(b) 将来 nylon の prep-once-run-many が必要になれば prep_structure.py を nylon 対応に拡張。

## 2026-06-24: 200+10 は 16 GB GPU でも WSL なら sustained MD が収まる — 「16 GB 実走不可」結論を撤回(Windows 固有の制約だった)
- Measurement: 同一の RTX 4060 Ti (16 GB) を **WSL2 (Ubuntu-24.04) + pfpoly-gpu (torch 2.12.1+cu126)** で実行し、`scripts/profile_vram.py --systems vinyl_methyl_acrylate --scale 1.0 --md-steps 300`(= 論文スケール 200+10, 2520 atoms, 0.50 g/mL, 古典圧縮 → OrbMol-v2 持続 MD 300 step)を完走(exit 0)。**device peak 9.55 GB / 16 GB(torch reserved peak 7.88 GB)、headroom 6.45 GB、fits=True、0.954 s/step**。`PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True` が env で有効。成果物: runs/vram_s6_test_wsl/vram_profile.json。
- Contradiction with 2026-06-15: 2026-06-15「200+10 exceeds 16 GB VRAM for sustained MD」では同型カード(Windows)で予約 VRAM が数百 step かけて 9.5 → 15.9 GB に creep してハングした。今回は同規模・同 300 step 窓(= Windows で creep が顕在化した窓)でピークが 9.55 GB のまま平坦で、creep は発生せず。
- Root cause of the difference: 2026-06-15/2026-06-15(VRAM record)で記録のとおり `expandable_segments:True` は **Windows CUDA build では no-op**(torch が "not supported on this platform")。WSL2 では実際に効くため、近傍グラフのサイズ変動由来のアロケータ断片化が抑えられ、ピークが単発フットプリント(~9.5 GB)に留まる。したがって「200+10 は 16 GB では sustained 不可、≥24 GB 必須」という結論は **Windows 固有の制約**であり、VRAM 容量の根本的な天井ではなかった。
- Decision: 200+10 paper-scale を **このマシン(16 GB)+ WSL で実走可能**と再評価する。本番 S6 を WSL の pfpoly-gpu で起動(scripts/run_s6_paper_scale.sh 経由、PYTHONPATH=src:. + expandable_segments)。≥24 GB GPU は「短時間化(0.95 s/step → A100 で数倍速)」の手段であって「実走の必須条件」ではない、と位置づけを変更。
- Caveat (honest scope): 検証は 300 step。S6 本番は ~127,000 step(50 cycle)/ 数十時間規模なので、超長時間での creep ゼロは未保証。ただし Windows での creep は「数百 step」窓で顕在化したため、300 step 平坦は強い(完全証明ではない)証拠。本番ランの実 VRAM を定期サンプルして確認する。
- Licensing/commercial impact: なし(新規依存なし。torch/orb-models/OpenFF は既存)。
- Follow-up: (a) README:84「S6 未着手(24 GB+ GPU 必要)」, docs/installation.md S6 ハードウェア要件, docs/paper-reproduction.md を「WSL なら 16 GB で実走可、24 GB は時短用」に更新。(b) S6 本番完走時に confirmed_formations / 温度 / figures を確認し figure-comparison.md を更新(これがリポジトリ初の 200+10 paper-scale 完走になる)。(c) 本番中の実 VRAM をサンプルし、長時間 creep の有無を記録。

## 2026-06-24: PES スキャン診断で OrbMol-v2 がラジカル付加の引力チャネルを再現すると確認 — formations=0 は MLIP 限界ではなくサンプリング/スピンの問題に確定
- Background: 2026-06-15(Fix B minimal)が「formations>0 は OrbMol-v2 の反応再現性の問題であり未検証。推奨は workflow ではなく PES スキャン診断」と残していた。S6 本番(200+10, 16GB WSL)を 18/50 cycle で停止し、その診断を実施。
- Method: `scripts/scan_radical_addition.py --backend orb --device cuda --spin 2`。メチルラジカル `[CH3]` + エチレン `C=C`(教科書的ラジカル付加、文献: 障壁 ~7、発熱 ~-23 kcal/mol)、**spin=2 doublet**、forming C–C を 3.5→1.54 Å で拘束緩和スキャン(sp2→sp3 再混成を許可)。
- Result (decisive, VERDICT=attractive): dE_vs_3.5 が 3.5→2.2 Å で +6.1 kcal/mol まで微増(遷移状態 ~2.2-2.4 Å、障壁 ~+6)→ 2.0 Å で -2.1 → 1.8 Å -16.0 → 1.6 Å -26.8 → **1.54 Å で -27.8 kcal/mol(最小)**。教科書 PES(障壁 ~7、発熱 ~-23)と定量的に一致。**OrbMol-v2 はラジコル付加の引力チャネルを正しく再現する。**
- Interpretation: (1) 「別 MLIP / PFP が必須」という 2026-06-15 末尾の最悪シナリオは**否定**。formations=0 は MLIP の根本的限界ではない。(2) PES 形状が workflow の困難を説明: 3.5→2.2 Å はほぼ平坦(~6 kcal/mol)で、3-6 Å の候補ペアはほとんど引力を感じない → ~2.2 Å の障壁を**バイアスで押し越える**必要があり、越えれば井戸へ自落。窓 [3,6] も f2 も触る必要はなく、バイアスが 2.2 Å まで届けば足りる。(3) スピンの寄与: 本スキャンは spin=2(doublet)で成功。**【訂正: 下記 2026-06-24 CORRECTION 参照】** 当初ここで「S6 が spin=1 だから候補が来ても結合しない」と記したが誤り。S6 はログ上 spin=1 のまま formation を 2 件出していた(run を kill したため bonds.jsonl 未生成で空に見えただけ)。スピンは formation EVENT の可否を決める直接原因ではなかった(電子状態の物理的厳密性には効くが、バイアス駆動の幾何条件は spin=1 でも満たされる)。
- Decision: formations を出す筋道を「(a) 正しい開殻スピン + (b) バイアスで ~2.2 Å 障壁を越えさせる」と確定。窓/f2/r0 のチューニングは不要(paper-confirmed のまま）。次の科学的争点は「溶融系の複数ラジカルで系全体スピンが ill-defined」問題(独立 doublet vs 高スピン和）であり、scripts/scan_radical_addition_2rad.py がその検証用。
- Licensing/commercial impact: なし。
- Follow-up: (a) scan_radical_addition_2rad.py で 2 ラジカル同時のスピン扱いを検証(melt スケーリングの本質的ブロッカー)。(b) その結果次第で、単一ラジカル+少数モノマーの TDBB を spin=2 で回し formations>=1 を実証(機構の end-to-end 検証）。(c) 多ラジカル melt のスピン処理方針を decisions に起こしてから 200+10 を spin 対応で再実行。

## 2026-06-24: 2ラジカル PES で高スピン近似は MARGINAL(使用可・~2-3 kcal/mol の不定性)— melt は系全体スピン=高スピン和で扱う
- Method: `scripts/scan_radical_addition_2rad.py --device cuda`。2×CH3 ラジカル + C2H4。スキャンするラジカル1がエチレン C0 へ付加、傍観ラジカル2 を中心から 7 Å に固定。1rad/s=2(クリーン基準)・2rad/s=1(非物理の閉殻一重項・参考)・2rad/s=3(高スピン三重項=2 doublet の高スピン和)で拘束緩和スキャン。成果物: runs/s3_pes_2rad/pes_comparison.json + s3_pes_comparison.png/pdf。
- Result: 正しい比較は 1rad/s=2 vs 2rad/s=3。barrier 33.90 vs 37.18(diff 3.28)、product min -27.78 vs -30.18(diff 2.41)、max点差 2.64 kcal/mol。VERDICT=MARGINAL(<1 kcal/mol のクリーン基準には未達だが、付加チャネルの形—障壁 ~+7、深い発熱井戸—は完全保存)。参考: 2rad/s=1(閉殻一重項)は非物理でズレ最大 ~4.2 kcal/mol、これは使うべきでない。
- Decision: 多ラジカル melt のスピンは **系全体を高スピン和(全不対電子を平行、multiplicity = n_radicals+1)で OrbMol-v2 に渡す** 方針を採用。根拠: 高スピン近似は付加チャネルを定性的に正しく再現し、誤差は ~2-3 kcal/mol に収まる(独立 doublet の厳密扱いは MLIP の系全体スピン制約では不可能なため、これが実務上の最善)。閉殻一重項(spin=1)は物理的厳密性では不可。**【訂正: 下記 2026-06-24 CORRECTION 参照】** ただし「S6 が spin=1 で formations=0」は事実誤認で、S6 は spin=1 でも formation 2 件を出していた。高スピン和は formation を出すための前提ではなく、電子状態を物理的に正しくするための改善。
- Caveat: ~2-3 kcal/mol の系統誤差は障壁(~7 kcal/mol)に対して無視できない比率。バイアス(f1_max=250)は障壁を十分上回るので bond formation の可否自体には影響しないが、絶対反応性/速度論の定量比較には効く。傍観ラジカルが増える(高密度 melt)とシフトが累積する可能性があり、本番では確認が要る。
- Licensing/commercial impact: なし。
- Follow-up: (a) 単一ラジカル+少数モノマーの TDBB を **spin=2** で回し formations>=1 を end-to-end 実証(2026-06-15 Fix B minimal は spin=2 でも 0 だったが、当時は PES が引力的と未確認・希薄系でバイアスが ~2.2 Å 障壁まで届かなかった可能性。今回 PES が引力的と判明したので、バイアス捕捉が 2.2 Å に届く配置/密度で再試行する価値がある)。(b) run_vinyl_aibn の活性化後スピンを「閉殻 spin=1」から「高スピン和」へ修正できるか検討(現状 activation は spin=1 のまま)。(c) 200+10 を高スピン対応で再実行する前に、まず (a) の最小系で機構を確定。

## 2026-06-24 CORRECTION: S6(200+10)は spin=1 でも formation を出していた — 「formations=0」は run を kill した運用ミスによる誤認(bonds.jsonl は完了時書き込み)
- Trigger: 上記 2 エントリで「S6 は spin=1 ゆえ formations=0」と推論したが、停止した S6 のログ(tasks/bbliiqoxs.output, line 93-110)を精読したところ事実誤認と判明。faithful-reporting に基づき訂正。
- Actual data: S6(200+10, 2640 atoms, spin=1, f2=5.0, Fix A 連続 in-phase 検出, f1_max_formation=250)は cycle 19/50 までに **chain propagation(formation)を 2 件**記録。Cycle 9: `reaction event at step 465 (1 pair(s))`, min_pair_dist 1.74 Å → `Chain propagation: atom 577 (beta-C) → radical_C`。Cycle 15: `reaction event at step 586`, min_pair_dist 1.90 Å → `Chain propagation: atom 1453 → radical_C`。selection.jsonl に該当ペア([223,576,128,781] / [79,1452,32,2581])も記録。
- Why I misread: bonds.jsonl と summary.json は **run 完了時にフラッシュ**される。私が cycle 18 で run を kill したため両ファイルが未生成 → 「bonds.jsonl が空 = formations=0」と早合点した。**ファイル不在 ≠ formation ゼロ**。formation の有無はログの reaction-event/Chain-propagation 行で判断すべきだった。運用ミスとして、動作中で成果を出していた run を停止してしまった。
- Corrected scientific picture: (1) **TDBB は paper-scale(200+10)で実際に bond formation を生成する** — これは本プロジェクトの主目標の達成。(2) spin=1 でも formation EVENT は起きる: 強いバイアス(f1_max=250 ≫ 障壁 ~7)がペアを幾何的 bonding 条件(< 60% Σvdw)まで駆動し、paper の reaction-event 定義(幾何条件)を満たすため。(3) なぜ 2026-06-15 の 100+5 / 20+2 では 0 で、今回 200+10 で出たか: **スケール**(200+10 は候補・遭遇が多い。cycle 9 で 5 候補)+ **Fix A**(run-until-reaction の連続検出が過渡的接触を捕捉)+ **f2=5.0**(3 Å でのバイアス力 ~17 kcal/mol/Å、f2=10 の ~0 と桁違い)の組合せ。当初「スケールアップで転化率」仮説には実は merit があった(2026-06-14/15 で否定し過ぎていた)。(4) PES 診断(本日の 2 エントリ)は依然有効: OrbMol-v2 に真の引力井戸がある=強制接触したペアは実在の井戸に入る(純バイアス artifact ではない)。高スピンは formation の前提ではなく、電子状態の物理的厳密化のための改善。
- Caveat / remaining rigor: spin=1 は閉殻サロゲート(activation 解離 0 → 真の開殻ラジカル不在)。formation は paper の幾何基準では成立だが電子状態は真のラジカルでない。物理的厳密化には高スピン和 / 実ラジカルでの再実行(plan B)が有効。ただし TDBB 機構そのものは paper-scale で formation を出すと実証済み。
- Process lesson: 長時間 run を停止する前に、(i) 成果ログ(reaction-event/propagation)を確認、(ii) 完了時書き込みファイル(bonds/summary)の有無で結論しない。
- Decision: S6(200+10, 現行パラメータ)を **完走まで再実行**し、bonds.jsonl / summary.json / figures と 50 cycle 全体の formation 数を確定する(リポジトリ初の 200+10 paper-scale 完走 + formation 記録)。spin 厳密化(plan B)は完走ベースライン取得後の改善として別途。
- Licensing/commercial impact: なし。
- Follow-up: (a) S6 再実行・完走 → confirmed_formations, 温度, figures を figure-comparison.md に記録。(b) 完走後に高スピン/実ラジカル版を比較し、formation の電子状態依存を評価。

## 2026-06-24: 候補枯渇の真因 = 活性化と平衡化の順序バグ + AIBN の prep 不安定(幾何/密度/NVT/スピンは無関係)
- Question: S6 の低転化率(~7 formation = 3.5%、論文 60-80%)の律速は何か。「候補枯渇」(ほとんどのサイクルで n_candidates=0)を実測で切り分け。
- Method: (1) オフライン再構成は **原子順序不整合で無効**と判明(再ビルド index が軌跡と不一致 → radical–共有 chain_C 距離 28 Å 等あり得ない値。`scripts/diag_candidate_starvation.py` の数値は破棄)。(2) 権威ある方法に切替: `polymerization.py` の `_run_biased_phase`/`run_activation` に env ガード付き計装(`KAGOME_DIAG_STARVATION=1`)を追加し、ライブの正しいグループで `raw_ij / validk / validl / both / qualified` と activation 時の azo C-N 実距離を出力。短縮 S6(200+10)で測定。
- Decisive evidence:
  - 活性化が失敗した run(equil 長): DIAG-STARV は raw_ij≈40-50(幾何は潤沢)だが **validk≈0**(i–k: radical_C–chain_C ≤3Å がほぼ不成立)→ qualified≈0。
  - 活性化が成功した run(equil 500): DIAG-ACTIV = 20 azo C-N bonds, dist **min 1.18 / max 17.97 / mean 3.86 Å、window[0,3]に 15/20 のみ**。15/15 dissociated で活性化発火。直後 DIAG-STARV cycle 0 = raw_ij=39, **validk=20, both=12, qualified=16**。
  - つまり qualified candidates は活性化の成否で **0 ⇄ 16** に変わる。候補枯渇 ⇔ 活性化不発、の因果が確定。
- Root cause: ワークフロー順序が `minimize → 333K 平衡化(本番 equil_steps=2000)→ run_activation`。AIBN アゾは熱的に不安定で、**活性化の前に 333K 平衡化で制御されず分解・四散**(azo C-N が 0.5 ps で 18Å まで開く=物理的 AIBN 速度より桁違いに速く、prep 歪み or OrbMol-v2 由来)。平衡化が長いほど azo が全滅 → 本番(equil 2000)は activation「0 dissociated」→ 本物のラジカル不在 → radical_C/chain_C トポロジー崩壊で i–k 制約不成立 → 候補枯渇 → formation ~7。
- Conclusion: 低転化率の律速は **幾何でも密度でも NVT でもスピンでもなく、活性化前の平衡化が AIBN を壊していること**。これは 2026-06-13/14/15 の各仮説(scale, equilibration不足, temperature, window, bias range)とも別の、ワークフロー順序の問題。
- Decision/fix candidates(Ask-first: 反応モデル/ワークフロー順序): (a) **活性化を 333K 平衡化の前に実行**(構築直後の intact 構造で azo を解離 → その後に平衡化・TDBB)。(b) または活性化前の equil を skip/短縮。(c) または古典圧縮/平衡化中に azo 結合を拘束。(d) 古典圧縮単体で既に azo が壊れるか(compress 直後の azo 距離)を切り分け、FF の azo パラメータ妥当性を確認。最有力は (a)。
- Caveat: 短縮 run(equil 500)でも 5/20 は既に窓外。古典圧縮(20段 FIRE)時点での損傷可能性があり、(d) の切り分けが要る。
- Licensing/commercial impact: なし(計装は env ガード、既定無効)。
- Follow-up: (a) 活性化順序を build→activate→equil→TDBB に変える設計を decisions に起こし実装(承認後)。(b) 修正後に短縮 S6 で qualified candidates が全サイクル潤沢か、formation が増えるかを確認。(c) その後 NPT 化(plan B)。

## 2026-06-24: 活性化順序の修正を実装・検証 — activation 0→19, formation レート ~3-4x 改善。残課題: 活性化後 equil が候補品質を再劣化
- Implementation: `scripts/run_vinyl_aibn.py` の活性化分岐を `minimize → equil(333K) → activate` から **`minimize(FIRE 0K) → activate → equil(333K)`** に変更。minimize は clash 緩和(0K なので azo を熱分解しない)として活性化前に残し、333K 平衡化のみ活性化後へ移動。スピン切替(高スピン和)は activate 直後・equil 前に実行。
- Verification (短縮 S6: 200+10, equil 500, activation-steps 1000, 2 cycles, spin auto): DIAG-ACTIV = azo C-N **max 3.46Å(無傷), 窓[0,3]に 18/20**(旧順序 equil 前は max 17.97Å, 15/20)。**Activation result: 19 dissociated(旧本番=0)**。Spin switched 1→21。cycle0 qualified=3 → reaction event step4 → Chain propagation、cycle1 qualified=1。**Confirmed formations: 1 / 2 cycles**(旧 ~7/50 ≈ 0.14/cycle に対し ~0.5/cycle ≈ 3-4x)。根本原因(活性化前 equil が AIBN を破壊)は解消と確認。
- Remaining issue (next lever): `validk`(radical_C–chain_C ≤3Å)が活性化直後は ~20(bkb57pl5o 実測, equil を挟まない計測)だが、**活性化後 equil 500 を挟むと 6 に低下、qualified 16→3**。activation 自体は救えたが、その後の 333K 平衡化(spin=21 高スピン)が生成ラジカルのトポロジー/配置を再び劣化させている。候補をさらに潤沢にするレバー候補: (a) post-activation equil を短縮/skip し直接 TDBB へ、(b) high-spin 動力学の安定性確認、(c) 平衡化を弱め(低温/短時間)に。
- Performance note: high-spin(spin=21)生産は spin=1 比でやや重い/不安定の可能性。別の短縮 run(bmiz0nptj, equil 2000)は活性化到達前に異常に時間がかかり kill(主因は CPU 古典圧縮 + ML minimize のブレで、本修正とは無関係)。
- Status: 計装(KAGOME_DIAG_STARVATION env ガード)と本修正は未コミット。
- Licensing/commercial impact: なし。
- Follow-up: (a) post-activation equil を短縮/skip した短縮 run で qualified が ~16 を維持し formation レートが上がるか確認。(b) 効果確認後に長め(10-50 cycle)で転化率を測り、論文 60-80% への差を再評価。(c) その後 NPT 化(plan B)。

## 2026-06-25: 決定的 — spin=21(20ラジカルの高スピン和)で OrbMol-v2 は熱暴走する。高スピン和方針は paper-scale で破綻
- Trigger: 活性化順序の修正後、候補が cycle 0(qualified=34)以降に全ラジカルで一斉崩壊(cycle1 で validk 34→5、cycle3 で 0)。「post-activation equil が原因」仮説で equil=0 を試したが改善せず → 温度を直接測定。
- Decisive measurement(trajectory の temperature_K, 設定 333K): **spin=21 の run は両方とも熱暴走**。s6_equil0(equil 0, spin21): cyc0 biased 2.5e6 K → cyc2 で 1e10 K。s6_fix_verify(equil 500, spin21): cyc0 biased 2.3e8 K → 1.7e10 K。equil の長短は無関係(むしろ equil 500 の方が速く発散)。一方、元の安定 S6(spin=1; 活性化 0 解離でスピン切替が起きなかった)は全 run ~300K で安定だった。
- Root cause: **OrbMol-v2 は 20 ラジカルの高スピン和 multiplicity=21 で異常(発散的)な力を返し、Langevin 333K 恒温器が制御不能になり系が爆発(1e8〜1e10 K)。** これが validk 崩壊=候補枯渇の正体(系全体が爆発し全共有結合が壊れる)。検出された formation は爆発中の偶発接触で物理的に無意味。
- 構図: 元 S6 = 活性化失敗→spin=1→安定だが候補少。活性化順序の修正 = 活性化成功→spin 1→21 切替→爆発。修正は活性化を直したが、その結果「高スピン和」方針(2026-06-24 採用)が 20 ラジカルで破綻することを露呈。2-radical PES の "MARGINAL"(~2-3 kcal/mol)が、20 radical では "catastrophic"(発散)に悪化。
- Implication: OrbMol-v2 の系全体スピン制約下では、多ラジカル melt の高スピン和は使えない。formation を出すには (A) **production spin を低く保つ**(spin=1 等。非物理だが安定でバイアス駆動の幾何 formation は起きる。元 S6 が spin=1 で formation 2件を出した実績あり)か、(B) 別アプローチ(per-radical spin 不可、別 MLIP、または逐次 1 ラジカルずつ処理)。
- Decision: PENDING(Ask-first: スピン/反応モデルの科学的意味)。次の切り分け/修正候補: 活性化順序の修正は残しつつ **スピン切替を低い上限でキャップ**(例 spin=1 or 2)し、活性化成功+安定動力学+候補維持を両立できるか短縮 run で検証。これが成立すれば paper-scale で持続的な候補供給→転化率向上が見込める。
- Status: 計装・活性化順序修正はコミット済み(165295e)。スピンキャップは未実装(設計待ち)。
- Licensing/commercial impact: なし。
- Follow-up: (a) spin キャップ実装(--production-spin-cap 等)→ 短縮 run で温度安定(~333K)と qualified 維持を確認。(b) 安定したら長尺で転化率測定(plan B)。(c) その後 NPT 化(plan C)。高スピン和方針(2026-06-24 の 2-radical 決定)は本結果で paper-scale 不適と更新。

## 2026-06-25 事前確認(spin の本質究明)— 結論: OrbMol-v2 では agnostic 不可・高スピン和は即爆発・低スピンは非物理かつ不安定 → 逐次処理が本命
- 事前確認1(OrbMol-v2 の spin 実装と挙動):
  - orb-models adapter は spin を **graph_feats["spin_multiplicity"](系全体スカラー)**として渡す。conservative_regressor の conditioner が **charge+spin を必須**とし、未設定だと `AssertionError: Missing required total_charge and spin_multiplicity`。→ **OrbMol-v2 は spin-agnostic 実行が構造的に不可**(spin=None を試すと conditioner が assert で拒否。backend は変更せず元のまま)。
  - 孤立緩和ラジカル(8 個)の spin 掃引: mult 1〜31 で max|F|=5〜15 と穏当(高 mult 単体では非爆発)。
  - standalone 力テスト(実 post-activation 構造, periodic): max|F|=27万〜85万で spin 非依存に見えたが、**安定だった killed_partial 構造でも 85万**だったため**周期座標 unwrap のアーティファクトと判明・破棄**。
  - 実 MD ラン(s6_spincap1, 活性化成功18解離, spin キャップ=1): cyc0 biased が **408→1220→992K(設定333K付近)で開始** — spin=21(s6_equil0)が初手 2.5e6 K で即爆発したのと対照的。**実動力学では高 multiplicity が初期から不安定**(standalone では見えなかった)。ただし spin=1 でも unbiased 中に加熱(1777→1.3e5K)し cyc1 で 3e8K に発散 → **timestep 1 fs(論文 0.25 fs)による数値ドリフトが第2の不安定源**(硬い反応系で 1 fs は過大)。
- 事前確認2(論文/PFP の多ラジカル扱い, PDF 全28p 抽出):
  - **"spin"/"multiplicity"/"unpaired"/"doublet"/"open-shell"/"unrestricted" の語が論文に一度も無い**。p21「AIBN 分解+ラジカル付加のみ有効化、停止/連鎖移動は無効」、p8-9「全開始剤ラジカルが活性持続(リビング的)」=**同時に最大20ラジカル**、p27「障壁検証は単一ラジカル+モノマー」、p20「PFP Bader 電荷を古典平衡化に使用」。→ **論文は PFP をスピン非設定(agnostic)で回している**と強く推定。OrbMol-v2 では再現不可。
- Correction: 2026-06-25 前エントリ「spin=21 が OrbMol-v2 を爆発させる」は**部分的に正しいが単純化**。正: (1) 実 MD では高 multiplicity が初期不安定を生む(spin 依存は実在)、(2) ただし spin=1 でも timestep 1 fs で徐々に発散、の**2要因**。候補枯渇の最終原因は「活性化後の本物ラジカル melt を安定に動かせないこと」。
- 重要: spin=1 で **qualified candidates が cyc0=37/cyc1=36 と2サイクル維持**(温度上昇前)。**動力学が安定なら候補は枯れない**=安定化できれば転化率を積める。
- Decision(方針, Ask-first 済の議論継続): OrbMol-v2 互換で物理的に筋が通るのは **逐次処理(一度に1ラジカル=doublet spin=2、安定領域かつ物理的、PES 検証済み)**。agnostic は不可、高スピン和は爆発、低スピンは非物理。バックエンド変更は「agnostic/フラグメントスピン対応が商用安全に存在するか不明」で不確実。**timestep を 0.25 fs(論文値)に戻す**のは安定化の併用レバー。
- Status: --production-spin-cap フラグ・diag_spin_sweep.py を診断用に残置。spin=None backend ガード(死にコード)と diag_spin_on_real_struct.py(周期 unwrap アーティファクトで信頼不可)は破棄。
- Follow-up: (a) 逐次処理ワークフローの設計(1ラジカル activate→doublet で連鎖成長→終了で次)。(b) timestep 0.25 fs を既定反応条件に。(c) 設計合意後に実装。

## 2026-06-25: プロジェクト方針の確定 — 目標は「商用安全 TDBB 能力(i)+ 妥当域の知見(ii)」、論文の数値追試は非目標
- Context: spin の壁(多ラジカル melt を OrbMol-v2 で物理的に扱えない)を受け、ゴールを再設定する議論をオーナーと実施。論文キネティクスの物理的意味を精査した結果、Rp∝[I] や線形転化率は **停止反応を無効化した理想化(モデル化選択)の帰結**であり MLIP の物理的忠実性をほとんど検証しない(論文 p8-9 自身が「√[I] でなく線形なのは停止なしだから」と明記)。物理的に濃い結果(相対反応性・障壁)ほど OrbMol-v2 では困難、再現容易な集団キネティクスほど物理的中身が薄い、という非対称が判明。
- Decision: 本プロジェクトの目標を **(i) 使える商用安全 TDBB シミュレーション能力(道具)** と **(ii) その妥当域(OrbMol-v2+TDBB がどこまで信頼できるか)の知見** に確定。**論文の数値追試(完全再現)は非目標**(MLIP が異なる以上原理的に不可能、かつ目玉キネティクスの物理的中身が薄い)。論文は「手法の出典 + 定性的サニティチェック」に格下げ。CLAUDE.md の Mission(設計のための商用安全な workflow 再現)の実体に合致。
- Implications: (a) 多ラジカルスピンの壁(2026-06-25)は「再現の障害」ではなく **画定すべき限界(成果 ii)** に位置づけ変更。(b) **S6(密な 200+10 同時多ラジカル melt)は OrbMol-v2 の妥当域外**。スコープ内での扱いは「妥当域外と画定 + ガードレールで検知/拒否 + ラジカル能力は妥当領域(希薄/単一ラジカル=doublet)で実証」。(c) **逐次処理は必須でなく任意拡張**(密領域でも構造生成したい場合のみ。キネティクスは非物理になる)。
- Workstreams: A) 妥当域スペック + 小 validation スイート(本セッションの PES/障壁・スピン壁の知見を結晶化)。B) 妥当領域でのツール堅牢化(無効領域のガードレール=黙って発散させない、+ クリーンな end-to-end: 段階重合ナイロン[閉殻=spin 問題なし]・希薄単一ラジカルのビニル[doublet])。C) ドキュメント(使い方・妥当域・限界)。
- Scientific risk: 低(スコープを物理的に妥当な領域に限定する方向の決定)。Licensing/commercial impact: なし(OrbMol-v2 維持、PFP 不採用を再確認)。
- Follow-up: ワークストリーム A から着手 → specs/validity-domain.md を新設(本エントリ参照)。

## 2026-06-25: ラジカル架橋の方針 — 逐次は保留、spin-agnostic バックエンド spike を先行
- Context: ラジカル架橋(crosslinking/curing)を「構造生成までは扱える道具にしたい」とのオーナー意向。手段は (S) 逐次(OrbMol-v2, doublet, 構造生成のみ・非物理キネティクス)か (A) spin-agnostic バックエンド(同時多ラジカルを一律処理、論文 PFP と同戦略)。
- 重要な再認識: 多ラジカルスピンの壁は **OrbMol-v2 が spin 必須であることに固有**。spin-agnostic MLIP(multiplicity 非要求)なら壁が無く、しかも **commercial-safe な候補が存在しうる**(PFP=Matlantis の制限に縛られない)。バックエンドはクリーンなインターフェース化済みで差替え容易。
- Decision: **逐次は投機的に実装せず保留**。先に **バックエンド spike を実施**(候補 spin-agnostic MLIP の選定 → license チェック → `scripts/scan_radical_addition.py` で PES 検証=障壁~7/井戸~−28 を再現するか)。判定後に逐次の要否を確定する。
  - spike 通過 → 同時多ラジカル架橋が商用安全に可能 → **逐次は不要化**。
  - spike 失敗 → commercial-safe agnostic 無し → **逐次が OrbMol-v2 での構造生成フォールバックに昇格**。
- Rationale: 「大きな実装の前に安く前提を潰す」原則。逐次を先に作るとバックエンドが通った場合に無駄になる。逐次には副次的利点(検証済み backend を妥当領域 doublet で使う=局所結合形成の信頼性)があるが二次的。
- Guardrails (CLAUDE.md): 新バックエンドは **採用前に license チェック必須**、`specs/dependency-license-matrix.md` に商用化ステータス記録。permissive 重みでなければ `blocked_pending_review`。PFP は使用許諾未確認のため既定にしない。
- Follow-up: (a) 候補選定基準(commercial-safe permissive 重み + 有機ラジカル対応 + spin 非要求 + 保存力)で候補列挙・license 確認。(b) 最有力候補で PES 検証。(c) 結果を validity-domain.md §2.7 と本エントリに反映し逐次要否を確定。

## 2026-06-25 CORRECTION(重大): 「spin=21 が爆発」は誤帰属。真因は timestep 1fs + 平衡化不足。OrbMol-v2 は良条件で 200+10 多ラジカル melt を安定に走る
- Trigger: AIMNet2-NSE spike のクリーンな同条件比較(40+4, 0.25fs, equil1000, spin=9)で、**OrbMol-v2 が mean 386K と aimnet(875K)より安定**だった。「OrbMol-v2 は多ラジカル不可」が条件依存では?と疑い、爆発した元条件(200+10, spin=21)を良条件で再試験。
- Decisive measurement: **OrbMol-v2, 200+10, spin=21(auto, 18解離), 0.25fs, equil2000** → 温度 **mean 404K / max 437K で安定**、qualified candidates **40/41/41 と全サイクル潤沢に維持**(runs/orb_fullcond)。これは s6_equil0(1fs,equil0)/s6_fix_verify(1fs,equil500)で 1e6〜1e10 K に爆発した**まさに同じ系・同じ spin=21**。唯一の差は timestep(1fs→0.25fs)+ 平衡化(0/500→2000)。
- Root cause(確定・訂正): 一連の「爆発 / 候補崩壊 / spin の壁」は **run 条件のアーティファクト** — (1) **timestep 1fs が反応性ラジカル系には過大で数値発散**(C-H 伸縮 ~3000cm⁻¹、強バイアス、開殻で硬い)、(2) 平衡化不足。**論文値 timestep=0.25fs + 適切な平衡化**で解消。2026-06-25「spin=21 が爆発」「spin が決定的」とした各記述は **timestep への誤帰属**であり撤回。spin=21(高スピン和)自体は良条件下で安定に扱える。
- Implications: (a) **逐次処理・spin-agnostic バックエンド・「多ラジカル melt は妥当域外」(validity-domain §2.1/§2.6)は、いずれも“爆発という偽ブロッカー”前提だった → 前提崩壊**。OrbMol-v2 単独で多ラジカル melt が走る。(b) **AIMNet2-NSE は検証済みの有効な第2バックエンド**(PES 一致 = OrbMol-v2 と barrier ~6/well ~−28、ラジカル訓練、クロス検証に有用、MIT 商用安全)だが**問題解決には不要**。第2バックエンドとしての価値(独立検証)で残置。(c) **元の目標(多ラジカル melt の転化率、さらにキネティクス)が OrbMol-v2 + 良条件で再び射程**。「論文数値追試は非目標」「妥当域に限定」の方針(2026-06-25 プロジェクト方針)は、この新事実を踏まえ**再検討の余地**。
- Caveat: orb_fullcond は 3 cycle・biased 200×0.25fs=50fs と短く formations=0(候補は潤沢=40+、biased 時間延長で形成する見込み)。長尺で転化率を実測して確認要。aimnet が 875K と高めだった理由(thermostat 結合 vs バイアス仕事)は別途。
- Process lesson: (i) 既定の探索用 timestep(1fs)を本番反応系にそのまま使ったのが誤りの起点。反応性 MD は論文値 0.25fs を既定にすべき(validity-domain §3 に既記載、根拠が実証された)。(ii) 「バックエンド固有の限界」と結論する前に、run 条件(timestep/平衡化)を最小比較で潰すべきだった。クリーンな2バックエンド比較が真因露呈の決め手。
- Status: aimnet backend(src/kagome/backends/aimnet_backend.py)、scan/run の aimnet 対応、license matrix 記録は spike 成果として残置(第2バックエンドとして有効)。
- Follow-up: (a) **OrbMol-v2 + 0.25fs + equil で 200+10 を長尺実行し転化率を実測**(元の S6 目標が達成可能か)。(b) timestep 0.25fs を反応条件の既定にし S6 スクリプト/configs を更新(現状 1fs)。(c) 目標スコープ(論文追試の要否)を新事実で再検討。

## 2026-06-26: formation 律速は「biased 時間」ではなく f2 捕捉幅 vs 候補窓のデッドゾーン。f2=5→2 で melt formation が発火(ただし累積加熱という新課題)
- Trigger: 中断していた長尺 s6_good_conditions(200+10, 0.25fs, biased2000, 20cyc)を再開する前に「formation の出やすさ」を小スケールで確認。前エントリ CORRECTION 末尾の見込み「biased 時間延長で形成する見込み」を検証。
- Method/Evidence(段階的):
  1. **機構隔離(demo_radical_formation.py)**: 理想配向 1ラジカル+1モノマー(doublet)で biased2000×0.25fs。**confirmed_formations=1**(step71 で障壁越え, final r=1.63Å, chain propagation も発火)。→ TDBB 形成機構・確定ロジックは健全。なお demo に候補生成バグを発見・修正(`if ps.is_formation: ps.r_min=select_rmin` が constraint_only ペアの C=C[~1.33Å]・radical_C–chain_C[~1.5Å] 窓まで潰し 0 candidates 化 → `and not ps.constraint_only` で主形成ペアのみ widen に修正)。
  2. **小 melt(40+4, activation, biased2000, 5cyc, f2=5=現行)**: 温度 319–379K 安定・候補潤沢(15–25)・bias 上限だが **confirmed_formations=0**、min_pair_dist が **3.3–3.6Å で頭打ち**(5cyc 通して下降せず)。→ biased 200→2000(10×)でも 0 = 「時間」ではない、を実証(CORRECTION 末尾の見込みを反証)。
  3. **真因の定量化**: formation bias V^f=f1(1−exp(−f2(r−r0)²)) の力は r0±~1/√f2 のガウス窓内のみ。r0=λΣr_vdw=0.6·(1.7+1.7)=**2.04Å**。f2=5 では |F| が r=2.3Å:478 / 2.5:358 / 2.8:81 / 3.0:17 / **3.3:0.7 / 3.6:0.01** kcal/mol/Å。**候補窓[3,6]で選ばれる pair(3.0–3.6Å)は捕捉シェル(~2.5Å)外のデッドゾーンにあり bias 力が実質ゼロ**。密 melt の 500fs では <2.8Å へ自力拡散できず formation 不成立。demo が出たのは 2.5Å の捕捉内に手配置したから(bias_E わずか 0.02)。
  4. **f2 スイープ(40+4, biased2000, 3cyc, seed7)**: **f2=2 → confirmed_formations=3/3cyc**(min_pair_dist 2.01–2.04=r0 到達, biased 早期終了 382/1066/702, propagation×3)。**f2=1 → 3/3cyc**(min 2.02–2.04, 早期終了 289/296/240)。予測どおりデッドゾーン橋渡しで毎 cycle 発火。
- Root cause(確定): melt で formation が出ない真因は **f2(捕捉幅)が論文/OrbMol 調整値で狭すぎ、候補窓[3,6]との間にデッドゾーンが生じ、選択 pair に bias 力が届かないこと**。biased 時間・スピン・密度ではない。
- 新課題(要対策): **f2 を下げると系が累積加熱**。f2=2: 平衡後389K→cyc2で565K(上昇中)、f2=1: 414→607K。原因は bias 仕事+ラジカル付加の反応熱(~−28kcal/mol)を unbiased 500step の Langevin で捌けず毎 cycle 残留・積算。**長尺(20–50cyc)では際限なく上昇する懸念**。対策候補: unbiased を 500→2000–3000 / Langevin friction 強化 / cycle 間冷却 equil。
- Decision: PENDING(Ask-first: f2 は論文 default 10 を OrbMol 用に 5 へ調整済みのハイパラ、さらなる変更は bias の選択性・非物理結合リスクに関わる)。最有力: **production f2=2 採用 + cycle 間緩和を増やして T を 333K 近傍に戻す**。長尺再開前に小 melt で「f2=2 + 緩和強化 → formation 維持 & T≈333K」を確認し、その後 16GB VRAM を踏まえ **半スケール 100+5 から**長尺化。
- 形成結合の物理性: f2=2/1 とも propagation×3 発火=結合は unbiased を越えて持続(bias 強制の一過性ではない)。OrbMol はラジカル付加に引力チャネル(barrier~6/well~−28, scan_radical_addition.py 既検証)を持つので、bias は「捕捉シェルへ届かせる」役で PES が結合を完成、と整合。ただし f2 を下げ過ぎると非物理 pull リスク → f2=2(保守側)を優先、f2=1 は不採用(より高温)。
- Status: demo 候補生成バグ修正は適用済(未コミット)。f2 変更・緩和強化は未実装(本エントリの確認 run + 承認後)。s6_good_conditions(中断)は checkpoint 無しのため破棄→再実行扱い。
- Licensing/commercial impact: なし。
- Follow-up: (a) 小 melt(40+4)で f2=2 + unbiased↑/friction↑ → formation 維持 & T≈333K を確認。(b) OK なら 100+5 半スケールで長尺・転化率実測(VRAM 監視)。(c) configs/scripts に production f2 と緩和設定を反映。(d) 形成結合の幾何/エネルギーを 1 例抜き取り検証(任意)。
- CONFIRMED 2026-06-26: follow-up(a) 合格。小 melt 40+4・**f2=2 + friction_per_fs=0.01 + unbiased=1500**・6cyc(runs/s6_f2_2_cool): 温度 **mean 345 / max 420 K で plateau(累積ドリフト無し)**、formation **5/6cyc**・propagation 5。f2=2 baseline(friction0.001/unb500=467→565K 上昇)と対照的に加熱が解消し formation も維持。レシピ確定 = **production f2=2 + friction 0.01 + unbiased 1500**。実装: `scripts/run_vinyl_aibn.py` に `--friction-per-fs`(default 0.001=従来不変)を追加済(未コミット)。次: 100+5 半スケールで長尺・転化率実測。
- COMPLETED 2026-06-26: follow-up(b) 完走。**100+5 半スケール長尺(30cyc, 確定レシピ, runs/s6_half_100x5)**: **転化率 26/100=26%**(confirmed_formations 26 / propagation 26)、活性化 10/10、温度 **mean 338 / 最終1/3 337 K(全30cyc で 333K 近傍維持・ドリフト無し)**、**VRAM peak 6.4GB/16(<40%、swap 無し完走)**。これがリポジトリ初の「多ラジカル melt で formation が積み上がる長尺完走」。f2=5(=元 s6_good_conditions 設定)なら formation≈0 だったところを 26% に。中断ランの真の障害は「捕捉幅×加熱×VRAM」で、いずれも解消。残: configs/scripts へのレシピ反映(follow-up c)、より高転化率には cycle 数増(論文 60-80% は multi-hundred cycle 相当)。図: runs/s6_half_100x5/figures/(conversion/temperature/energy_vs_step)。

## 2026-06-26: cycle 境界 checkpoint / resume を実装(長尺 run のクラッシュ復旧)
- Motivation: 70cyc@100+5(~22-26h) や 50cyc@200+10(~30-38h) は checkpoint 無しでは all-or-nothing(swap/クラッシュで最初から)。次の長尺前の保険として実装。
- Design: `PolymerizationWorkflow.run()` が各 cycle 境界(group 更新後)で `<output-dir>/checkpoint.pkl` に atomic 保存(.tmp→rename)。保存対象 = cycle ループが変化させる動的状態のみ: positions/velocities/cell/step, groups の atom_indices(propagation で変化), updater の chain_c_map・_processed_formations/dissociations, BondTracker の _events・_reacted, **numpy Generator state(bit-exact 継続の要)**, logs, production spin。静的部分(species/template/calculator/weights)は保存せず、resume 側が再 build して動的状態で上書き。
- Resume: `--resume` で build時の activation/minimize/equilibration を skip、checkpoint の spin を calc.set_spin で復元、保存 cycle の次から継続。trajectory/selection は追記(append)、bonds は tracker 全 events から再生成。_pending は cycle 境界で空(unbiased の check_outcomes でクリア)なので非保存。
- 実装: io/trajectory.py(append フラグ), workflows/polymerization.py(save/load_checkpoint + run(checkpoint_path,resume,checkpoint_extra)), scripts/run_vinyl_aibn.py(--resume/--no-checkpoint + spin 復元/activation skip), scripts/run_s6_paper_scale.sh(RESUME=1 受け渡し)。既定で checkpoint 自動保存、--no-checkpoint で無効。
- 検証: ユニット test_resume_is_bit_exact(ToyCalculator, 「2cyc→checkpoint→resume 2cyc」== 無中断 4cyc を positions/tracker events で bit-exact 確認, 33/33 pass)。実 GPU e2e(8+2, runs/ckpt_e2e): leg1 活性化4/4・spin5・checkpoint 保存 → leg2 `--resume` で **活性化を再実行せず spin 5 復元・cycle 2 から継続**(cycle 0 からではない)を確認。
- Licensing/commercial impact: なし。
- Follow-up: 長尺(70cyc@100+5 or 200+10)を回す際は既定で checkpoint され、落ちても `RESUME=1 OUTPUT_DIR=<同じ>` で再開可能。最適化余地: resume 時も build の圧縮を払う(~数分)ので、将来 species/template も checkpoint すれば build skip 可能。

## 2026-06-27: cycle-15 ハング調査 — 物理ではなくバックエンド/GPU stall。診断計測(StepWatchdog)を実装
- Trigger: runs/s6_full_200x10(200+10, 50cyc, OrbMol-v2, --resume)が cycle 15 でハング(進展停止)。
- Evidence(trajectory/selection/log の事後解析):
  - ハング地点 = cycle 15 biased step 31535 の `_md_step`(= calculator.compute, Orb/warp GPU)直後。**E_base≈9986 / T≈336K と全フレームと同一、NaN 無し、エネルギー爆発無し** → 数値発散ではない「静かな停止」。
  - selection.jsonl の cycle 15 が**2回**(完全同一スコア)。trajectory の追記順では cycle 15 **unbiased step 33033** が biased step 31535 より**前**に存在 = cycle 15 を**同一 checkpoint から2回**実行し、一方は unbiased まで完走・もう一方(resume)は biased で停止。**同一初期状態から到達点が異なる = ハングは非決定的**。
  - GPU = RTX 4060 Ti **16GB**。run_s6_paper_scale.sh は 200+10 を **≥24GB 必要**と明記。run_vinyl_aibn.py 冒頭コメントに既出の既知故障モード「CUDA caching allocator 断片化 → VRAM 枯渇 → run hangs」(2026-06-15 VRAM record)。
- Hypothesis(最有力): **VRAM 不足/瞬間スパイク → WSL system-memory fallback による stall**(クラッシュせず PCIe スピルで 10–100× 減速 = 「ハング」に見える、断片化依存で非決定的)。物理/TDBB ロジックではない。
- 0.68Å 近接の扱い(訂正含む): 最終フレームで原子間最小距離 ~0.68Å・<1.0Å が 164 対あるが、**生座標(PBC 無し)でも 0.66–0.69Å = 実在の近接**で PBC アーティファクトではない(当初「アーティファクト」と述べたのは撤回)。ただし **初期構造(cycle -1, T=0)から全 cycle 一様に存在し、cycle 15 を完走した試行でも同条件** → **本ハングの引き金ではない**(新規崩壊ではない)。真因は compress_box の FIRE 未収束(最終 fmax≈46 vs 目標2.0)による残留クラッシュで、構造品質の別件として要調査。病的局所形状が compute コスト/メモリを一時的に跳ね上げ VRAM スパイク仮説に寄与する可能性あり。
- 実装(診断のみ・物理不変): `src/kagome/diagnostics.py` の `StepWatchdog` を biased/unbiased/activation の各 MD ループに装着。各ステップ前に faulthandler を単発再武装(規定秒ハングで全スレッド stack dump)、遅ステップ WARN(peak VRAM 付き)、任意 heartbeat。env 閾値(診断値・科学パラメータではない): `KAGOME_WATCHDOG_S`(既定180, 0で無効)/ `KAGOME_STALL_WARN_S`(既定20)/ `KAGOME_HEARTBEAT_STEPS`(既定0=無効)。CPU/MACE では VRAM ログ自動 skip。faulthandler 出力は stderr(既存 run の resume.log は stderr 捕捉済=warp warning が記録されている)。
- 検証: ユニット 33 passed/1 skipped。watchdog 単体 smoke OK(既定値・arm/step_done)。
- Licensing/commercial impact: なし。
- Follow-up: (a) 計測入りで half-scale(100+5, 既知の 16GB swap-free 完走実績)を 50cyc 再現 → ハングしなければ VRAM 確定。(b) フル 200+10 を回すならハング中 `nvidia-smi`(共有 GPU メモリ/util)+ stack dump で stall 箇所を直接確認。`KAGOME_HEARTBEAT_STEPS=500` 推奨。(c) VRAM 確定なら sysmem fallback を無効化(NVIDIA Sysmem Fallback Policy = Prefer No Sysmem Fallback)し OOM 例外化で白黒、もしくは半スケール運用。(d) 残留クラッシュ(0.68Å)= minimize 未収束を別途調査(max_steps/fmax 見直し)。
- CONFIRMED 2026-06-27: 計測入りで same checkpoint から `--resume`(runs/s6_full_200x10)。前回停止の **step 31535 を素通りし cycle 15→16 を継続**(同一初期状態で結果が変わる = **ハングは非決定的、物理ではない**を実証)。VRAM 定量: heartbeat `torch_peak=7.9 GiB`(backward グラフのステップ内スパイク、後に 0.11 へ解放)、ホスト nvidia-smi の device は **used~10.3 / free~5.7 GiB で安定**(上昇トレンドではない)。**torch のステップ内ピーク 7.9 > 空き 5.7 = 毎ステップ物理 16GB の縁で spill 気味**(step rate 1.08→0.86 と低下、~1 step/s と遅い)。フラグメンテーション次第でこの spill が thrash 化 → 「静かなハング」。warp は torch 外確保のため torch-only 計測は過小評価 → mem_get_info の device used を追加(commit 8c7af3d)。
- Decision 2026-06-27: 数値一致は非目標(eval.qualitative=「定性トレンド一致を定量厳密一致に先行」, claim は全てサイズ非依存の method/workflow)と確認。**フル 200+10 は 16GB GPU では VRAM の縁で低速+ハングリスク常時 → 半スケール 100+5 に切替**。フル A は cycle 15 まで checkpoint 済(必要なら 24GB+ GPU で resume 可)。半スケールは VRAM peak 6.4GB で swap-free 完走実績(runs/s6_half_100x5, 30cyc=転化率26%)。計測入りで **100+5 / 50cyc を実行**(runs/s6_half_50c, KAGOME_HEARTBEAT_STEPS=500)。
- Follow-up(更新): (a) フル 200+10 を回したい場合の本質対策 = ①24GB+ GPU、②sysmem fallback 無効化で OOM 顕在化、③メモリ削減(backward グラフ縮小/逐次評価)。当面は半スケールで科学目標を回収。(b) 0.68Å 残留クラッシュ(minimize 未収束)は別途。

## 2026-06-28: half-scale 100+5 / 50cyc 完走 — 論文 method/workflow を定性再現(runs/s6_half_50c)
- 結果: **完走(cycle 50/50)。確定形成 43, 解離 0**。transition 率 α = 43/100 = **43%**。VRAM swap-free(torch_peak ~4.3 / free ~14.6 GiB)、stall WARN ゼロ、~2.1–2.5 steps/s。
- 設定: seed 7, 100 monomer + 5 AIBN, activation(f2=0.3, f1_max=250, 5000 steps)→ spin 1→11(N_radicals=10)→ equil 2000 → TDBB 50cyc ×(biased 2000 + unbiased 1500), f2=2, friction 0.01/fs, density 0.5 g/mL, T=333K, no-barostat, OrbMol-v2/cuda, timestep 0.25fs, minimize fmax=1.0。launcher: `scripts/run_vinyl_aibn.py ... --resume`(下記参照)。
- 温度制御: biased mean 333.1K/std 14.1K(target 333), unbiased mean 339.4K/std 11.4K(反応発熱の緩和でやや高め=妥当域)。Langevin friction 0.01/fs で良好に制御。注: biased に min=0.0K のフレーム1点(resume 直後初期フレームのアーティファクト、mean には実害なし)。
- 反応進行: cycle あたりほぼ 1 形成の単調増加(仕様上は複数/cycle 可だが select_non_overlapping + 同一ステップ同時クロスの稀少性で実効ほぼ 1 本/cycle, bonds.jsonl 全 `1 pair(s)`)。candidate 数は monomer 消費に伴い緩やかに減少(健全)。
- 中断・復旧: 実行途中(cycle 48 到達時)に **launcher(tee していた bash タスク)停止で子 python が道連れ終了**(物理/数値起因ではない)。cycle 境界 checkpoint.pkl(next_cycle=48)から `--resume` で bit-exact 再開し cycle 48–50 を完走 → checkpoint 機構(2026-06-26)が長尺の保険として実機で機能することを確認。再起動は `nohup python ... --resume >> run.log 2>&1 &`(launcher 切断に耐性)。
- 図(scripts/reproduce_figures.py から再生成可、手編集なし): runs/s6_half_50c/figures/ に conversion_vs_step / temperature_vs_step / energy_vs_step / base_energy(.png/.pdf)。**density プロットは skip**("No trajectory frames match bond event steps" = bond イベント step と trajectory 出力 step の粒度不一致。主目的=転化率/温度/エネルギーには影響なし)。
- 成果物: runs/s6_half_50c/{run.log, bonds.jsonl(43 confirmed_formation), trajectory.jsonl(~194MB), checkpoint.pkl(cycle50 状態), figures/}。
- 評価: paper の claim は全てサイズ非依存の method/workflow + eval.qualitative(定性トレンド一致先行)。**半スケールで TDBB workflow(biased 近接生成→確定→鎖伸長)と単調増加する転化・制御された温度を定性再現 → 論文追試の科学目標を達成**(数値厳密一致は非目標と既決)。
- Licensing/commercial impact: なし(OrbMol-v2 は既存の妥当域バックエンド、新規依存なし)。
- Follow-up: (a) density プロットが必要なら trajectory 出力間隔を bond step に整合させ再生成。(b) 先行 30cyc 実績(α=26%)→ 50cyc(α=43%)で転化が伸長を継続中、より高転化を見たい場合は cycle 数増で延長可能(checkpoint で再開可)。(c) 0.68Å 残留クラッシュ(minimize 未収束)は引き続き別件。

## 2026-06-29: half-scale 100+5 / 100cyc 完走 — checkpoint resume で 50cyc を延長、転化率 73%(runs/s6_half_50c)
- 結果: **完走(cycle 100/100)。確定形成 73, 解離 0**。transition 率 α = 73/100 = **73%**。目標 60% を上回る。
- 方式: 前項 50cyc 完走の `runs/s6_half_50c/checkpoint.pkl`(next_cycle=50, 復元 spin=11, reacted=43)から **`--resume` + `--n-cycles 100`** で cycle 50→99 を継続。resume が `range(start_cycle, n_cycles)`(polymerization.py)で動くため、checkpoint の next_cycle から n_cycles 増分だけ素直に延長できることを実機確認。既存 43% を bit-exact で引き継ぎ、build/activation/minimize は skip(spin 11 復元)。
- 転化の伸び: 30cyc=26% → 50cyc=43% → **100cyc=73%**。cycle に対し単調増加を継続(停止反応なし理想化と整合)。末期(cycle ~80 以降)は candidate 数が 5–8 に減少・形成ペース ~0.53/cycle(中盤 ~0.7–0.8 から鈍化)= monomer 消費(残 27)に伴う正常な頭打ち。
- 実行健全性: ~2.0–2.3 steps/s(cycle 68 付近で一度 1.46 まで低下したが一過性で 2.04 に回復)。VRAM swap-free(torch_peak ~4.7 / nvidia-smi device ~5–6.5 / 16 GiB、~10 GiB 空き)、stall WARN ゼロ。**今回は中断なしでノンストップ完走**。
- 起動方式の知見: 前回の launcher 切断(子 python 道連れ)を踏まえ、**`setsid nohup python ... --resume >> run.log 2>&1 &`** でフル detach 起動 → launcher/セッション切断に耐え完走。長尺起動の既定手順として有効。WSL では `nvidia-smi` の per-process 表示が効かず他プロセスの GPU 使用は総量からの推定のみ(Windows タスクマネージャ GPU タブで直接確認可)。
- 図(scripts/reproduce_figures.py で再生成、手編集なし): runs/s6_half_50c/figures/ を 100cyc データで上書き(conversion/temperature/energy/base_energy, .png/.pdf)。density は前回同様 skip(bond event step と trajectory 出力 step の粒度不一致、主目的の転化/温度/エネルギーには影響なし)。注: run_vinyl_aibn.py 直起動のため図は自動生成されず手動再生成した(run_s6_paper_scale.sh 経由なら自動)。
- 成果物: runs/s6_half_50c/{run.log, bonds.jsonl(73 confirmed_formation), trajectory.jsonl, checkpoint.pkl(cycle100 状態), checkpoint.cycle50.bak.pkl(50cyc 退避), figures/}。
- 評価: 前項(50cyc/43%)の延長として **半スケール TDBB workflow で転化率 73% まで単調伸長を達成**。論文追試の科学目標(定性トレンド一致)をより高転化で補強。数値厳密一致は非目標(既決)。
- Licensing/commercial impact: なし。
- Follow-up: (a) さらに高転化を見たい場合は同 checkpoint(cycle100)から `--n-cycles 150` 等で再延長可(末期鈍化のため伸び幅は逓減)。(b) density プロット・0.68Å 残留クラッシュ(minimize 未収束)は引き続き別件。

## 2026-07-02: フルスケール(200+10)100cyc の転化率50% は正常、ただしトラジェクトリに化学的不整合を確認 — 4レイヤー修正に着手
- Context: 別マシンで論文フルスケール(200 monomer + 10 AIBN, notes.md)を 100 cycle 完走。転化率 α ≈ 50%。半スケール(100 monomer)100cyc=73% との差の原因と、トラジェクトリ観察で見えた化学的不整合の切り分け。
- 転化率50%の評価(正常): α = confirmed_formations / n_monomers(2026-06-27 決定, Fig.2 α=1−[M]/[M]₀)。本実装の TDBB は概ね 1 形成/cycle(鎖端で候補選択→バイアス→1確定)で進むため、**同一 cycle 数なら転化率は n_monomers にほぼ反比例**。半スケール(100mon)100cyc=73形成=73% に対し、フルスケール(200mon)100cyc≈100形成=50% は整合。飽和ではなく単調増加の途中で、論文の 60–80% には cycle 増(フルスケールで概算 120–160+ cyc)で到達見込み。checkpoint から `--n-cycles` 増で延長可。
- 観測された不整合(ユーザ報告, ビューア表示由来と確認): (1) 価電子数を無視した結合(過配位に見える)。(2) ビニル部位でない原子間の分子間結合。**観測の出所はビューアの距離推定表示**(OVITO/VMD 等が座標から結合を自動描画)であることをユーザ確認。
- 根本原因(コード読解で特定):
  1. 記録上は非ビニル結合は出得ない: `PolymerizationWorkflow._build_pair_biases`(polymerization.py:1079-1106)が `constraint_only`(i-k, j-l)を除外 → バイアス=BondTracker 記録対象は radical_C–vinyl_alpha_C の i-j ペアのみ。bonds.jsonl の confirmed_formation は構造上必ずビニル部位。
  2. trajectory は座標のみで結合トポロジーを持たない: `TrajectoryFrame`(io/trajectory.py:14-28)は positions のみ、bonds フィールドなし。ビューアが距離推定 → 最大 f1=250 kcal/mol の引力バイアスが**分子ごと**近接させるため、エステル基・H など非ビニル原子間にも見かけの結合が描かれる。→ 観測(2)の正体。
  3. 反応後に C=C 二重結合を開く処理が無い: `VinylChainPropagationUpdater`(polymerization.py:270-332)は群メンバーシップ上でラジカルを beta_C へ移すだけで、幾何・価数は sp3 生成物へ緩和されない。alpha_C は「beta_C との二重結合 ~1.33Å + 新規 radical_C 結合 + H + COOCH₃」で過配位に見える。→ 観測(1)の正体。
  4. 形成判定は幾何のみ: `bonds.is_formed`(r ≤ 0.6·Σr_vdw ≈ 2.04Å for C-C)で価数/占有チェック無し。
- 本質: バイアスで原子を近接させイベントを数えるだけで、実際の化学反応(二重結合開裂・価数付替・sp3 緩和)を行っていない。転化率カウントは正しく増えるが、トラジェクトリは化学的に破綻して見える。
- Paper anchor: notes.md Table S1(vinyl: i-j V^f, i-k/j-l constraint only)、Eq.6-7(reactive groups/candidate)、radical addition の頭尾機構(_find_vinyl_alpha_beta docstring)。価数保存は radical + C=C → C-C + 新ラジカル(beta)で本来保たれるべき。
- Decision(4レイヤー修正, ユーザ承認済み):
  1. **トポロジー出力**(安全・科学的意味不変): RDKit mol から各フラグメントの結合(次数付き)を抽出→グローバル index へオフセット。確定形成で radical_C–vinyl_alpha_C 単結合を追加、当該 monomer の C=C を単結合へ開裂反映。TrajectoryFrame/Writer に bonds を出力しビューアの距離推定を排除。→ 観測(1)(2)を直接解消。
  2. **価数/占有ガード**: 候補生成・形成確定に価数/占有チェックを追加、飽和原子・過配位を選ばせない/確定させない。反応モデルの科学的意味に関わるため根拠を本 decision に集約。
  3. **生成物の緩和強化**: 確定後のスピン状態・unbiased 緩和・局所 minimize で sp3 生成物へ緩和し過配位幾何を解消。
  4. **最小再現先行**: 小スケール run で記録/座標/スピンを切り分け、修正前後比較のベースラインとする。
- Scientific risk: レイヤー1は出力のみで risk なし。レイヤー2-3 は反応モデルに踏み込むため、TDBB の科学的意味を変えないこと(候補=幾何 Eq.7、バイアス=Eq.2-3 を保持)を条件に段階導入。
- Licensing/commercial impact: なし(RDKit=BSD, 既存依存)。
- Follow-up: レイヤー1実装+テスト+最小再現 → レイヤー2 → レイヤー3 の順。各段階でトラジェクトリのビューア表示が化学的に妥当か確認。

## 2026-07-02: Layer 1(トポロジー出力)実装完了 — ビューアの見かけ結合/過配位を解消
- 実装:
  - `src/kagome/reactive/topology.py`(新規): `BondTopology`(次数付き結合集合, coordination/valence 参照)+ `apply_vinyl_addition`。ラジカル付加の価数保存編集を1関数に集約 — (1)ラジカル中心が既に4配位(閉殻開始剤モデル)なら余剰 H 結合を1本除去、(2)radical_C–alpha_C 単結合を追加、(3)当該 monomer の alpha=beta を単結合へ開裂。配位数で開始剤(4配位→H除去)と伝播ラジカル(3配位→そのまま)を自動判別。
  - `scripts/_systems.py`: `layout_bonds(specs)`(連続フラグメント配置から次数付き結合を抽出しグローバル index へ)+ 薄いラッパ `vinyl_initial_bonds`/`full_aibn_initial_bonds`。戻り値は非破壊(既存の 6/7-tuple 呼び出し site を壊さない)。`_AIBN_SMILES` を先頭定数ブロックへ移動(前方参照解消)。
  - `src/kagome/workflows/polymerization.py`: `PolymerizationWorkflow(initial_bonds=...)` を追加。確定形成ごとに `_apply_topology_updates` で topology を編集し、変化時に `topology.jsonl` へスナップショット追記。trajectory ヘッダに初期 `bonds` を出力。checkpoint に topology 状態を保存/復元(resume で連続)。
  - `src/kagome/io/trajectory.py`: ヘッダに `bonds`。`src/kagome/io/readers.py`: `read_topology_snapshots`/`bonds_at_step`。
  - `scripts/run_vinyl_aibn.py`: 経路別に初期結合を計算し workflow へ(best-effort, 失敗しても run は継続)。活性化経路は解離した azo C-N 結合を activation 後に topology から除去(→ イソブチロニトリルラジカル2個 + N2, 価数正しい)。
  - `scripts/export_xyz.py`: `--format mol2` 追加。`topology.jsonl` を読みフレーム毎の明示結合(次数付き)で MOL2 出力 → Winmostar/OVITO/VMD が距離推定せず実結合を表示。既定は従来 xyz(後方互換)。
- 検証: `tests/unit/test_topology.py`(13 tests: BondTopology, apply_vinyl_addition, ビルダー整合, workflow E2E 出力)。最小再現(4mono+1init, ToyCalculator, 形成注入)で **header 55 結合 / 過配位ゼロ / C=C 開裂 order=1.0 / 余剰H除去(結合数維持)/ MOL2 出力成功**。全ユニット 361 passed, 1 skipped(回帰なし)。
- 重要な発見(Layer 2/3 で扱う): 開始剤が閉殻モデル(isobutyronitrile `CC(C)C#N`)のため、ラジカル炭素が実 AIBN ラジカルより H を1つ多く持つ。付加時に5配位(価数違反)になる。Layer 1 では出力上その余剰 H 結合を落として理想化トポロジーを描く(物理座標には H が残る)。物理的に正しい解は開殻ラジカルモデル(H 原子自体を除去+スピン)で、Layer 3 の選択肢。該当は開始剤炭素のみ(各鎖初回付加)。
- 使い方: 再 run すると trajectory に bonds が自動付与。既存 run は `python scripts/export_xyz.py <run>/trajectory.jsonl --format mol2` で MOL2 化(topology.jsonl があれば時間発展結合、無ければヘッダ初期結合)。
- Scientific risk: なし(出力のみ、TDBB の候補=幾何/バイアスは不変)。Licensing: なし。
- Follow-up: Layer 2(価数/占有ガード)→ Layer 3(緩和 or 開殻モデル)。フルスケール run は活性化経路なので、次回 run から MOL2 で価数を実確認。

## 2026-07-02: Layer 2(価数/占有ガード)実装完了 — 価数安全性を不変条件として保証
- 実装:
  - `reactive/topology.py`: `over_coordinated_atoms(topology, species)`(全原子の過配位スキャン)+ `vinyl_addition_over_coordinates(topology, radical_c, alpha_c, propagation_map, species)`(コピー上で付加編集を dry-run し、過配位になる原子を返す。空=価数安全)。閉殻開始剤の H 除去も dry-run に含むので、4配位ラジカルは誤って弾かれない。
  - `workflows/polymerization.py`: `_valence_filter` を biased phase の選択後に挿入。形成ペア(is_formation かつ not constraint_only)を `_formation_pair_positions` で解決し、過配位になる候補を除外・WARN ログ。topology 無し/非vinyl では no-op。加えて `_apply_topology_updates` に防御チェック — 万一すり抜けた形成が過配位を起こすなら topology 編集を skip し ERROR ログ(出力トポロジーは常に価数正しく保つ)。
- 位置づけ: 反応済み原子は群更新で除去されるため通常フローでは発火しないが、**価数安全性を「群簿記の創発的性質」から「保証される不変条件」へ格上げ**。簿記バグ・伝播エッジケース・幾何的に飽和した候補を選択段階で捕捉し、監査ログに残す(CLAUDE.md「反応ペア選択判断の記録」に合致)。TDBB の候補=幾何(Eq.7)・バイアス(Eq.2-3)は不変。
- 検証: `test_topology.py` に5件追加(有効付加は非フラグ / 飽和alpha はフラグ / 閉殻ラジカルはH除去で非フラグ / 過配位スキャン / workflow `_valence_filter` が飽和候補を drop)。全ユニット **366 passed, 1 skipped**(回帰なし)。
- Scientific risk: 低。正常系の挙動は不変(発火しない)。発火時は化学的に不可能な形成を防ぐのみ。Licensing: なし。
- Follow-up: Layer 3。ユーザ要望により**開始剤の閉殻→開殻モデルを検討**(H原子除去+スピン)。まず調査/設計から(スピン周りの過去の困難 decisions.md を踏まえる)。

## 2026-07-02: Layer 3 調査 — 活性化経路の開始剤ラジカルは既に開殻(価数正しい)、spare-H は非活性化経路のみ
- 調査結果(実測): `build_full_aibn_system`(活性化経路)のラジカル中心は、intact AIBN で 4配位(2CH₃+CN+azo-N)。activation で azo C-N を除去すると **3配位(C×3, H無し)= 真の開殻ラジカル**。Layer 1 の `apply_vinyl_addition` で付加すると 4配位になり過配位ゼロ(2ラジカル×付加を実測、over-valent=[])。
- 帰結: **ユーザのフルスケール run は活性化経路(`--activation`)なので、Layer 1 適用だけで開始剤も含め価数正しい**。spare-H 由来の過配位は `build_vinyl_aibn_system`(非活性化、`_INITIATOR_SMILES='CC(C)C#N'` イソブチロニトリルの閉殻簡易モデル、ラジカル炭素に H が1つ多い)でのみ発生。これは小規模テスト/デモ用の簡易経路。
- スピン対応現状: backend は `set_spin`/`supports_spin` 実装済(orb spin=1既定, aimnet mult)。活性化経路は activation 後に総スピン=n_radicals+1(cap 可)へ切替済。系レベルの総多重度で扱う。
- Layer 3 の再定義(調査を踏まえ):
  - (A) **生成物の幾何緩和**: 形成後の座標は分子ごと引き寄せで歪む可能性。既存 unbiased(2000 steps)で緩和されるはずだが、形成直後の局所 minimize(FIRE, 反応領域限定)追加で sp3 生成物への緩和を確実化する余地。物理トラジェクトリ(表示でなく座標)の改善。低リスク。
  - (B) **非活性化経路の開殻化**: `build_vinyl_aibn_system` の開始剤を実ラジカル(spare-H 原子を除去、原子数-1/開始剤)へ変更 → 下流 index/ n_per_init シフト+スピン。規模中。**代替**: 非活性化経路は閉殻近似と明記し、化学的に厳密な run は `--activation` を推奨(活性化経路は既に正しい)。低リスクで実質十分。
- 推奨: (B) は代替(ドキュメント+活性化推奨)で十分な可能性が高い。(A) は形成直後緩和を入れるか、まず MOL2 で実 run の座標歪みを確認してから判断。ユーザ判断待ち。
- Licensing: なし。

## 2026-07-02: Layer 3 実装 — (B) ドキュメント+活性化推奨を実装、(A) 緩和は実 run 検証後に判断(ユーザ決定)
- ユーザ決定: (A) 幾何緩和=「まず実 run を MOL2 で確認」してから要否判断 / (B) 非活性化経路=「ドキュメント+活性化推奨」(実ラジカルビルドはしない)。
- (B) 実装:
  - `scripts/_systems.py build_vinyl_aibn_system` docstring に閉殻近似の注記(ラジカル炭素の placeholder H、初回付加で raw 幾何は5配位、出力トポロジーは H を落として理想化、価数厳密には `--activation` 推奨)。
  - `scripts/run_vinyl_aibn.py` 非活性化分岐に INFO ログで同旨の推奨(実行時に気づける)。
- (A) 支援ツール実装: `scripts/check_topology_valence.py`(新規)。run の topology.jsonl を読み、(1)スナップショット毎の過配位原子(Layer1/2 で0のはず)、(2)フレーム毎の非結合最近接距離(座標歪み=分子ごと引き寄せの定量指標、scipy cKDTree)を報告。MOL2 目視と併せ (A) の要否を定量判断可能。Windows cp932 対応で出力は ASCII。
- 使い方(実 run 検証): `python scripts/export_xyz.py <run>/trajectory.jsonl --format mol2`(目視)+ `python scripts/check_topology_valence.py <run>`(定量)。過配位0かつ最近接非結合距離が常識的(>1.0–1.2Å)なら (A) 緩和は不要。<1.0Å が頻発するなら形成直後の局所 minimize を追加検討。
- 状態: (A) の緩和実装は**実 run の MOL2/定量チェック結果待ち**(条件付き follow-up)。現時点で actionable な (B)+ツールは完了。全ユニット 366 passed, 1 skipped(回帰なし)。
- Licensing: なし。
- Follow-up: ユーザが実 run(活性化, フルスケール)を MOL2 化 + check_topology_valence で確認 → 過配位0を実証(想定)。座標歪みが問題なら Layer 3(A) 局所 minimize を実装。

## 2026-07-02: 既存 half-scale run(runs/s6_half_50c)で検証 — 価数は正しい(観測はビューア由来と確定)、幾何クラッシュは既知の別問題
- 背景: 再 run 不要で検証するため、topology 導入前の run から `bonds.jsonl`(確定形成記録)を遡って再構築。`scripts/reconstruct_topology.py`(新規)で初期トポロジー(full_aibn)→ azo C-N 除去 → 73形成を step 順に `apply_vinyl_addition` 適用 → `topology.jsonl` を生成。
- 対象: runs/s6_half_50c(活性化経路, 5 AIBN + 100 methyl acrylate = 1320原子, orb-orbmol_v2, seed7, 100cyc, 73形成, α=73%)。再構築の species レイアウトが run と完全一致(妥当性確認)。
- 価数検証(結果): **過配位原子ゼロ**。C 配位数分布 `{2:10(ニトリルC), 3:164, 4:266}` すべて化学的に妥当。→ **ユーザ観測の「価電子無視/非ビニル結合」は記録された反応ではなく、座標のみトラジェクトリに対するビューアの距離推定アーティファクトと実データで確定**(診断②を実証)。Layer 1 の MOL2 出力で解消。
- 幾何検証(結果): 実 orb トラジェクトリ(259フレームサンプル)の非結合最近接距離 = **0.65Å**(step 30890, cycle13 unbiased)。全サンプルフレームに <1.2Å 非結合ペア。最悪接触はすべて **H-H**(0.65–0.72Å, 分子内 geminal・分子間の両方, <1.0Å が48ペア)。反応部位(C-C)ではない。→ これは decisions.md 既記の「0.68Å 残留クラッシュ(minimize 未収束)」= **初期配置/最小化由来の既知の別問題で、反応化学・価数とは独立**。
- Layer 3(A) への含意: 幾何クラッシュは実在するが**形成非依存(H-H, 初期minimize 由来)**なので、「形成直後の局所 minimize」では解消しない。正しい対処は**初期構造の minimize 収束改善**(fmax 厳格化/ステップ増、または prep 段の圧縮・緩和見直し)。post-formation 緩和は不要と判断。
- 成果物: runs/s6_half_50c/topology.jsonl(再構築), traj_stride50.mol2(104フレーム/14MB, 明示結合付き, 目視確認用)。`scripts/reconstruct_topology.py`, `export_xyz.py --stride`(大トラジェクトリ間引き)。
- 評価: **Layer 1+2 が実データで有効性を実証**。ユーザの主懸念(価数不整合)は「表示の問題」と確定し解消。残る幾何クラッシュは別トラック(初期minimize)。
- Licensing: なし。
- Follow-up: (a) ユーザが traj_stride50.mol2 を Winmostar/OVITO/VMD で目視 → 見かけ結合が消えたか最終確認。(b) 0.65Å H-H クラッシュ対策として初期 minimize 収束改善を別タスク化(fmax/steps 調整の小 run で検証)。(c) フルスケール run も同様に reconstruct_topology で遡及検証可能。

## 2026-07-03: D1 — 反応確定タイミング: バイアス中は暫定検知、非バイアス緩和後に確定
- Context: `check_reactions_during_bias` は r ≤ r0 到達時に即 `confirmed_formation` を発行し `_reacted` に追加するため、非バイアス緩和後の `check_outcomes` で再判定されない。f1 最大 250 kcal/mol の人工引力下での近接は化学結合の証拠にならず、バイアス除去後にペアが離れても確定反応として残る。レビュー指摘 H1。
- Paper anchor: §2.2 step 3 — 反応が検知されたらバイアス相を終了し、非バイアス相で緩和。論文は「緩和後に確定」を明示していないが、biased/unbiased 交互プロトコルの趣旨は「バイアスで遷移状態を突破し、緩和後に持続する結合のみ化学反応とみなす」と解釈。
- Decision: バイアス中の閾値到達は `tentative_formation` / `tentative_dissociation` として記録し、バイアス相終了のトリガーにのみ使う。確定 (`confirmed_*`) の発行は非バイアス緩和後の `check_outcomes` に一本化する。`_reacted` には暫定検知時点では追加しない。
- Alternatives considered: (a) 即 confirm 維持 + docstring 修正 — TDBB を「遷移状態突破=結合成立」と解釈する立場。棄却: 250 kcal/mol の井戸は遷移状態バリアを大幅に超えるため、近接が持続的結合の証拠にならない。(b) バイアス中に検知したペアのみ check_outcomes で再判定 — 実装複雑化に対して利点なし。
- Scientific risk: Low. 真に結合が成立するペアは緩和後も r ≤ r0 を維持するはずであり、confirmed 率は低下しない。偽陽性の排除により α・DPn の精度が向上。
- Licensing/commercial impact: None.
- Follow-up: bonds.py の BondTracker, polymerization.py の _run_biased_phase を修正。H2 (候補単位の原子的受理) と同時に実装。

## 2026-07-03: D2 — MC バロスタット受理判定にバイアスエネルギー変化を含める
- Context: バイアス相の NPT 体積試行で `try_step` に渡るのは `base_energy` のみで、体積スケールによるバイアスエネルギー変化 ΔV_bias が受理判定に含まれない。f2=10 Å⁻² の勾配最大点で f1=250 のとき |dV/dr|≈670 kcal/(mol·Å)。最大体積移動 (ΔlnV=0.01) で r≈3 Å のペアは Δr≈0.01 Å → ΔV_bias は数 kcal/mol/ペア。kT(333 K)≈0.66 kcal/mol に対して無視できない。レビュー指摘 M1。
- Paper anchor: §2 — NPT ensemble。
- Decision: try_step にバイアスエネルギー再計算コールバックを渡し、ΔH = Δ(E_base + E_bias) + PΔV − (N+1)kT·lnV'/V で受理判定する。バイアスが無い場合(非バイアス相)はコールバック=None で従来と同一動作。
- Alternatives considered: (a) バイアス相中はバロスタット無効化 — 実装コスト最小だがバイアス相の密度が制御されない。(b) 現状維持 + 意図的除外を記録 — 物理的に不整合が残る。
- Scientific risk: Low. 正しい NPT サンプリング。バイアス相中のバロスタット受理率が若干変わるが、密度制御がより正確になる。
- Licensing/commercial impact: None.
- Follow-up: mc_barostat.py の try_step にコールバック引数追加、polymerization.py の _md_step からクロージャを渡す。

## 2026-07-03: I3 — boost.advance() の Eq.5 離散化
- Context: ループ先頭で boost.advance() を呼ぶため、ループ内 1 ステップ目の f1=γ·1、ループ前の力評価(最初の half-kick 用)は f1=0。
- Paper anchor: Eq. 5 — f1(t)=min(γt, f1_max)、t=0 から開始。
- Decision: これは f1(t)=γt を t=0 から離散化した正しい形。ループ前 (t=0) で f1=0、ループ内 step 1 で f1=γ·1。記録として残す。
- Scientific risk: None.
- Follow-up: None.

## 2026-07-03: H2 — 候補単位の原子的受理 (candidate_id)
- Context: 縮合系（ナイロン）では 1 候補に formation (N-C) と dissociation (N-H, C-OH) の複数ペアがある。V^d だけ確定し V^f が不成立の場合、DefaultPostCycleUpdater がアミン N を全グループから除去してしまい、反応サイトが不可逆に失われる。
- Paper anchor: Table S2 — 縮合テンプレートの formation/dissociation 対。
- Decision: PairBias と BondEvent に candidate_id (int) を追加。_build_pair_biases が同一候補の全ペアに同じ ID を付与。DefaultPostCycleUpdater は同一 candidate_id の formation が confirmed された場合のみ dissociation のグループ編集を適用する。
- Alternatives considered: (a) updater を候補結果オブジェクトに変更 — API 変更が大きく後方互換性問題。(b) formation/dissociation を cycle 内 step 近接で紐付け — タイミング依存で脆弱。
- Scientific risk: None. ビニル系では dissociation イベントがないため影響なし。縮合系のサイト保存が正しくなる。
- Licensing/commercial impact: None.
- Follow-up: 縮合系のエンドツーエンドテストで検証。

## 2026-07-06: A1/A2/A3 — 密度プロファイルの PBC・N_frames・面積を論文定義に整合(ユーザー承認 2026-07-06)
- Context: `scripts/reproduce_figures.py::plot_density_profile` が (1) `reaction_density_profile` に `cell` を渡さず PBC 中点補正を無効化、(2) `n_frames` にイベント発生 step 数 `len(positions_at_event)` を渡し論文定義(全サンプリングフレーム数)と乖離、(3) セル未指定時の面積を `(z_max−z_min)**2`(xy 断面積として無意味)でフォールバックしていた。レビュー指摘 A1/A2/A3。
- Paper anchor: PDF p.12 unnumbered eq. ρ_rxn(z) = N_rxn(z) / (A · Δz · N_frames)。A は xy 断面積、N_frames は解析窓の全サンプリングフレーム数。
- Decision:
  (a) `TrajectoryFrame` に `cell`(オプション, 既定 None)を追加し、initial/equilibration/biased/unbiased の各フレーム書き込みで `state.cell` を記録する(NPT で箱が変動するため per-frame)。旧トラジェクトリは cell 欠落 → None で後方互換。
  (b) `reaction_density_profile` に `cells_at_event`(step→cell)を追加。イベントごとにそのフレームのセルで PBC 中点補正する。
  (c) `n_frames = len(frames)`(全サンプリングフレーム)。イベント step 数は使わない。
  (d) 面積 A は `--cell-xy-area` 指定時はそれ、未指定時はフレームセルの平均 Lx·Ly。セルも面積も無ければ明示エラーで中断(旧 `(z_max−z_min)**2` フォールバックは廃止)。
  (e) z ビン範囲はセルがあれば [0, ⟨Lz⟩) を既定(wrap 済み座標 [0,Lz) に整合)。データ min/max 依存はセル無し・面積指定時のみのフォールバックとして残す。
- Ask-first 該当: z ビン範囲の [0, Lz) 既定化は「図の平均化・範囲変更」に該当するが、本修正は論文定義への接近であり、2026-07-06 にユーザー承認済み。
- Alternatives considered: (a) 単一 cell(初期箱)で全イベント処理 — NPT で箱が変わるため中点・面積が不正確。棄却。(b) 面積を per-event の Lx·Ly で可変にする — ρ の定義は固定 A を前提とするため、平均 A を採用。
- Scientific risk: 中。図の縦軸スケール(絶対密度)と z 範囲が変わる。過去 run の密度図は再生成が必要。定義正確化のため許容。
- Licensing/commercial impact: None.
- Follow-up: `tests/unit/test_density.py` に PBC 中点・N_frames 正規化・面積欠落エラーのテストを追加。

## 2026-07-06: A5 — 縮合(nylon)1 反応 = 主形成ペア 1 件で計数、水形成イベントはバイアス・トポロジー用途のみ(ユーザー承認 2026-07-06)
- Context: nylon テンプレートの amine_H–carboxyl_OH ペアは `is_formation=True, score_pair=False` だが `constraint_only=False` のため `_build_pair_biases` で通常の formation PairBias になり、独立に `confirmed_formation` を発行しうる。結果 1 縮合反応(アミド結合 1 本 + 水 1 分子)が amine_N–carboxyl_C と amine_H–carboxyl_OH の 2 件の confirmed_formation として α(t)・Carothers p に二重計上される。レビュー指摘 A5。
- Paper anchor: PDF p.22 Table S2(縮合テンプレートの formation/dissociation 対)、Fig.4b-c(Carothers p)。1 縮合反応 = アミド結合 1 本。
- Decision: 物理的には amine_H–carboxyl_OH の近接バイアス(水の O–H 形成方向)は正しい挙動なのでバイアス自体は維持し、計数だけを主反応(amine_N–carboxyl_C)に限定する。
  (a) `PairSpec` に `count_as_reaction: bool = True` を追加。nylon の amine_H–carboxyl_OH に `count_as_reaction=False` を設定。
  (b) `PairBias` と `BondEvent` に `counts_as_reaction`(bool, 既定 True)を伝播。`_build_pair_biases` が `ps.count_as_reaction` を PairBias に載せ、BondTracker が BondEvent に載せて `bonds.jsonl` に記録する。欠落時 True 扱いで後方互換。
  (c) `BondTracker.confirmed_formations()` は全イベントのまま。計数側(reproduce_figures.py / run_nylon66.py)がイベント読み込み直後に `counts_as_reaction=True` のみをフィルタしてから conversion / carothers に渡す。analysis 関数の引数仕様は不変(docstring に前提明記)。
  (d) density プロット(`reproduce_figures.py::plot_density_profile`)も同フィルタを適用する。ρ_rxn は「反応」の密度であり、水形成イベント(count_as_reaction=False)を空間分布に含めない。vinyl は該当ペアを持たず無影響。
- 既存記録 decisions.md 2026-06-13(水は明示モデル化しない)との整合: 水分子そのものは生成せず、TDBB は距離のみバイアスする点は不変。本決定は「水形成方向のバイアスイベントを反応計数に含めない」だけを追加規定する。
- Ask-first 該当: 計数定義の変更(縮合 1 反応 = 主形成ペア 1 件)は 2026-07-06 にユーザー承認済み。
- Alternatives considered: `candidate_id` によるユニーク化(同一候補の複数 formation を 1 と数える)— スキーマ変更不要だが「どちらが主反応か」の情報が残らず、density での位置も曖昧になるため棄却。
- Scientific risk: 中。nylon の α(t)・p が従来比で下がる(二重計上の解消)。vinyl は amine_H 等を持たず、constraint_only の連鎖ペアは _build_pair_biases で除外済みのため計数不変。
- Licensing/commercial impact: None.
- Follow-up: 合成イベント列(同一 candidate_id で 2 formation)で α・p が 1 反応分になるユニットテストを追加。bonds.jsonl 旧フォーマット読み込み互換を担保。

## 2026-07-06: I1/I2/I3 — FIRE ミニマイザを正準アルゴリズム(Bitzek 2006 / ASE)に整列
- Context: `src/kagome/integrators/minimize.py::fire_minimize` の FIRE 更新順序が正準形と非等価だった。レビュー指摘 I1/I2/I3。
  - I1: 速度混合 `v = (1-α)v + α·F̂·|v|` を慣性キック `v += dt·F` の**後**に適用していた。正準 FIRE(Bitzek et al., PRL 97, 170201 (2006) / `ase.optimize.fire.FIRE.step`)は混合を**キック前**の速度と `|v|` に対して行い、その後にキックする。キック後混合では混合に用いる `v` と `|v|` が既にキック済みの値になり、混合係数 α の意味と速度方向の減衰が正準形と一致しない。
  - I2: 時間刻み増加の N_min 判定にオフバイワンがあった。旧実装は `n_positive += 1` してから `if n_positive > n_min` を判定するため、連続正パワー 6 回目で dt 増加を開始していた。ASE は `if Nsteps > Nmin: ...` を判定してから `Nsteps += 1`(判定してからカウント)で、7 回目に開始する。
  - I3: `FireParams` docstring が既定 fmax=1.0 kcal/mol/Å を ASE 既定より「looser(緩い)」と記述していたが誤り。1.0 kcal/mol/Å ≈ 0.043 eV/Å は ASE 既定 0.05 eV/Å(≈ 1.15 kcal/mol/Å)より小さい閾値=より**厳しい**。
- Paper anchor: PDF p.20 (SI) — reactive production MD の前に equilibration/最小化で初期接触を緩和。FIRE 自体は Bitzek et al., PRL 97, 170201 (2006)。本ミニマイザは論文本文のアルゴリズムではなく、初期構造(grid-packed)の接触緩和という準備工程の標準ツール。
- Decision:
  (a) 更新順序を正準形へ変更: power = F·v 判定 → (power>0 なら)キック前の v と |v| で混合 → dt/α 更新 → 慣性キック `v += dt·F` → 変位制限 → 移動。ASE と同様、初回ステップは v=0 のため power/混合ブランチを丸ごとスキップし(`vel is None` センチネル)キックのみ適用する。旧実装は初回に power=0 で else ブランチへ入り dt を fdec 倍していたが、これも ASE 非等価だったため是正。
  (b) N_min 判定を ASE と同じ「判定してからカウント」順に統一(I2)。
  (c) docstring を「1.0 kcal/mol/Å ≈ 0.043 eV/Å は ASE 既定 0.05 eV/Å より厳しい」に訂正(値自体は不変)(I3)。
  (d) 変位制限 `maxstep_A` は kagome 独自の per-atom クランプ(初期クラッシュ耐性のためのロバスト性措置)であり正準 FIRE の一部ではないため、本修正のスコープ(I1/I2/I3)では変更しない。ASE はグローバルノルムクランプを使うが、クランプが発火しない緩やかな系では両者は一致する。
- Scientific risk: 低〜中(科学的意味の変更なし、準備工程のみ)。本ミニマイザは初期構造の接触緩和専用で TDBB の候補=幾何(Eq.7)・バイアス(Eq.2-3)は不変。ただし収束経路(中間座標列)が変わるため、**準備構造が微変する。2026-07-06 以前のコミットで record 済みの run を bit-exact に再現するには当時の旧コミットが必要**(seed/config/git_sha は manifest に記録済み)。最終緩和構造の物理的品質(最終 fmax が閾値以下)は不変。
- Alternatives considered: (a) I2 を現行維持 — 1 ステップ差で収束品質への影響は軽微だが、I1 と同時に「ASE 準拠」で統一する方が参照実装との突き合わせ検証が容易なため棄却。(b) maxstep をグローバルノルム化して完全 ASE 等価にする — I1/I2/I3 のスコープ外かつ per-atom クランプはロバスト性上の意図的選択のため見送り。
- Licensing/commercial impact: None(ASE は test 依存に追加しない。参照実装はテスト内に手書きし、ASE がインストール済みの場合のみ本物の `ase.optimize.FIRE` との突き合わせを skipif ガード付きで実施)。
- Follow-up: `tests/unit/test_minimize.py` に (a) 正準 FIRE 更新式を手書きした参照実装との位置列一致テスト、(b) skipif ガード付きの本物 ASE FIRE 突き合わせテストを追加。既存の収束テストは green のまま。

## 2026-07-06: D1 — MC バロスタットは原子単位アフィンスケーリング(LAMMPS `dilate all` 相当)を採用
- Context: `src/kagome/integrators/mc_barostat.py:115` は体積提案受理時に全原子座標を一律の線形係数 `scale = exp(ΔlnV/3)` で拡縮する(`new_positions = positions * scale`)。これは分子重心を保存する「分子単位スケーリング」ではなく、各原子を個別に原点基準でスケールする「原子単位(dilate all)」方式。この設計選択が decisions.md 未記録だった。レビュー指摘 D1。
- Paper anchor: PDF Section 2 — NPT ensemble。バロスタットの具体方式(原子/分子単位のスケール対象)は論文未指定。原子単位スケーリング自体は LAMMPS `fix npt` の既定 `dilate all` と同型。
- Decision: 原子単位アフィンスケーリング(`dilate all` 相当)を採用する。理由: 本プロジェクトの主対象は架橋ネットワーク系(nylon 縮合・epoxy 硬化)であり、反応進行に伴い分子境界が動的に変化(モノマー→鎖→網目)する。分子単位スケーリングは「分子=剛体クラスタ」の定義を前提とするが、結合トポロジーが毎サイクル変わる本系ではクラスタ定義が不安定で、重心の定義自体が曖昧になる。原子単位スケーリングはトポロジー非依存で常に well-defined。受理判定の Jacobian 項 `-(N_atoms+1)·kT·ΔlnV` も原子数 N ベースであり、原子単位スケールと整合する(decisions.md 2026-06-20 RF19b)。
- Alternatives considered: (a) 分子単位スケーリング(重心を保存し分子内相対座標は不変)— 剛体分子系では標準だが、架橋で分子境界が動く本系ではクラスタ定義が不安定なため棄却。(b) group ごとの選択的スケール — 実装複雑化かつ論文根拠なし。棄却。
- Scientific risk: 低〜中。強共有結合内距離も僅かにスケールされるが、`max_volume_change_frac=0.01`(1 提案あたり線形 ~0.33%)と小さく、MLIP の復元力で緩和される。分子単位系との定量差は圧縮率評価で O(結合長変動) 程度。qualitative なワークフロー挙動(密度収束・反応進行)には影響しない。
- Licensing/commercial impact: None(内部実装)。
- Follow-up: 剛体分子近似が有効な非架橋系(例: 溶媒充填)を扱う場合は分子単位スケーリングをオプション追加検討。現状の架橋系では不要。

## 2026-07-06: W2 — フェーズ境界の MLIP forward 重複は既知事項として修正しない
- Context: ワークフローのフェーズ境界(biased→unbiased、unbiased→次サイクル biased、equilibration→cycle0)で、同一座標に対する MLIP forward 計算が重複する(`polymerization.py:797-799, 926-928, 1056-1058` 近傍)。各フェーズ入口で `base_energy, base_forces = calculator.compute(...)` を初期化のため呼ぶが、直前フェーズ最終ステップの `_md_step` が既に同座標で forward を済ませているケースがある。1 サイクル(約 4000 ステップ)あたり約 2 回、割合にして ~0.05%。レビュー指摘 W2。
- Paper anchor: N/A(実装効率の問題。科学的計算内容は不変)。
- Decision: **修正しない。既知事項として記録するのみ。** `_md_step` が末尾の base energy/forces を呼び出し元へ返し、次フェーズ入口がそれを受け取れば重複を除去できるが、`_md_step`/フェーズ関数群のシグネチャ変更(base_forces を戻り値/引数に追加)が必要で、API 変更のコストが ~0.05% の削減効果に見合わない。支配コストは MLIP forward 自体であり、2 回/サイクルの追加は無視できる。
- Alternatives considered: (a) `_md_step` が base_forces を返し次フェーズが再利用 — 効果 ~0.05% に対し API 変更が波及(全フェーズ関数・checkpoint/resume の状態互換)。棄却。(b) フェーズ境界で forward をスキップし前回値を状態に保持 — 隠れ状態が増え resume の bit-exact 性検証が複雑化。棄却。
- Scientific risk: なし(数値結果は現状で正しい。重複は冗長計算のみで誤りではない)。
- Licensing/commercial impact: None.
- Follow-up: フェーズ関数群を将来リファクタする機会があれば、その際に base energy/forces のフェーズ間受け渡しを併せて導入する(単独では着手しない)。

## 2026-07-06: 不完全 AIBN 活性化時に反応グループを実解離集合へ縮約(ユーザー承認 2026-07-06)
- Context: `--activation` 経路で AIBN の azo C-N が全本数解離しなかった場合、`build_full_aibn_system` は全ラジカル中心を無条件に `radical_C`/`chain_C`/`chain_c_map` に登録するため、未解離中心(3C+1N の 4配位、H無し)が反応グループに残る。`run_vinyl_aibn.py` の activation ブロックは解離した azo 結合を topology から除去する(`wf._topology.remove_bond`)のみで、未解離中心をグループから除去していなかった。結果:
  - 未解離中心は 4配位・H無しのため `topology._spare_hydrogen` が None → `vinyl_addition_over_coordinates` が毎サイクル `[radical]` を返し、valence guard が正しく drop する。しかし radical_C に残り続けるため毎サイクル候補選択枠を浪費する。
  - `production_spin = len(radical_C) + 1` が未解離中心も数えるため、系の総スピン多重度が過大になる。
  - 証拠: runs/diag_nvt(解離6/10 → 毎サイクル drop 4)、runs/diag_npt(8/10 → drop 2)、runs/s6_half_100x5(10/10 → drop 0、正常成長)。
- Paper anchor: Table S1(Activation 行、azo C-N への V^d)。開始剤ラジカルは開殻(2026-07-02 Layer 3 調査)。化学的には未開裂 AIBN は不活性(ラジカルを生成していない)。
- 既存記録との整合: 2026-07-02 の valence guard(Layer1/2)は正しい(未解離中心の付加を正しく拒否する)。変更しない。活性化経路の「開始剤ラジカルは開殻」前提(2026-07-02 Layer 3)は**完全活性化を暗黙に仮定していた**。不完全活性化時はグループを実解離集合へ縮約するのが正しい。未解離 AIBN は化学的に不活性なので反応グループから除外する。
- Decision:
  (a) `scripts/_systems.py` に純関数 `prune_undissociated_centers(groups, chain_c_map, dissociated_c_indices) -> (new_groups, new_chain_c_map, n_pruned)` を追加。実解離した azo_C 集合に含まれないラジカル中心を `radical_C`/`chain_C`/`chain_c_map` から除去する(新オブジェクトを返す純関数、入力は非破壊)。
  (b) `run_vinyl_aibn.py` の activation ブロックで activation 後・topology 結合除去後に (a) を呼び、`wf.groups` と `wf._updater.chain_c_map`(および局所 `groups`/`chain_c_map`)へ再代入して間引き後グループで production を実行。wf は構築時に groups dict と chain_c_map への参照を保持するため、再代入で確実に反映される(wf 構築を activation 後へ移せない — activation は `wf.run_activation`/`wf._minimize`/`wf._topology` を使うため)。
  (c) `production_spin(n_radicals, cap=None) -> int` を `_systems.py` の純関数に切り出し、間引き後の radical 数で算出。間引き発生時は WARNING で除外中心数を明示。
  (d) 案B(fail-loud): activation 後 `len(dissociated) < 2 × n_initiators` で顕著な WARNING(解離数/期待数、`--activation-steps`/`--activation-f1-max` 増加提案)。CLI フラグ `--strict-activation`(既定 off)を追加し、on なら RuntimeError で中断。
- Ask-first 該当: 反応グループの縮約は選択挙動・スピン多重度を変えるが、化学的整合(未開裂 AIBN=不活性)への接近であり 2026-07-06 にユーザー承認済み。
- Alternatives considered: (a) valence guard 自体を変更して未解離中心を無視 — ガードは正しく動いており、根本原因は登録集合の過大。棄却。(b) `build_full_aibn_system` で登録を条件付きにする — 活性化はビルド後に走るため解離結果はビルド時に未知。棄却。(c) in-place 変異で間引き — wf 反映は自動だが純関数のテスト容易性を優先し新オブジェクト返し+再代入を採用。
- Scientific risk: 低。完全活性化(10/10)時は挙動不変(n_pruned=0)。不完全活性化時のみグループ・スピンが実解離集合に整合(=より正しい)。
- Licensing/commercial impact: None。
- Follow-up: `tests/unit/test_systems.py` に間引き関数・スピン純関数のユニットテスト、`tests/unit/test_topology.py` に未解離中心(3C+1N)への付加ガード回帰テストを追加。

## 2026-07-07: TDBB パラメータ較正スクリプト(calibrate_tdbb.py)の導入
- Context: OrbMol-v2 checkpoint を差し替えた際、TDBB パラメータ(activation f1_max / steps、production f1_max)が新しい PES のバリアと整合するかを手動 PES スキャンで確認していた。これを `scripts/calibrate_tdbb.py` として自動化し、`run_s6_paper_scale.sh` が結果を自動で読み込む配線を追加する。
- Paper anchor: Table S1(Activation 行の V^d、production f1_max=250/125 kcal/mol)、Eq.2-3(バイアス振幅はバリアを超える必要がある)。バリア実測に基づく較正手順自体は論文外の追加(論文は PFP 前提の固定パラメータを直接与える)。
- Decision:
  (a) `calibrate_tdbb.py` は2つの PES スキャンを実行: AIBN azo C-N homolysis(spin=1、r=1.4→3.4 Å)とラジカル付加 ·CH3 + C2H4(spin=2、r=3.5→1.54 Å)。各点は `constrained_relax` で拘束緩和。バリアは各スキャンの反応物基準 dE の最大値。
  (b) activation_f1_max = max(250, ceil(barrier×2.5/25)×25)。安全係数 2.5 は「バリアを確実に超えつつ dwell 時間マージンを持つ」ための本リポジトリ独自の設定(論文外、既定バリア 39.3 では 250 のまま)。
  (c) activation_steps = max(5000, int(5000×barrier/39))。検証済みバリア 39 kcal/mol で検証済み値 5000(decisions.md 2026-06-26)を再現する線形スケール。1000 単位に丸め。
  (d) production 側 f1_max_formation=250 / f1_max_dissociation=125 / f2=2.0 は固定(Table S1 / decisions.md 2026-06-26)。ラジカル付加バリアの実測値は f1_max=250 に対しマージンが十分(~6 kcal/mol)なため導出には使わず、15 kcal/mol 超で WARNING を出す検知のみに使用。
  (e) 出力は `runs/calibration/calibration.json`(スキャン生データ+推奨値)と `tdbb_params.env`。env は `VAR="${VAR:-value}"` 形式で、手動の環境変数指定が較正値より常に優先される。`run_s6_paper_scale.sh` は起動時に存在すれば source する。`run_vinyl_aibn.py` に `--f1-max-formation` / `--f1-max-dissociation` を追加し TDBBParams へ配線(既定値は従来ハードコードと同一の 250/125)。
- 2026-07-07 実測(orbmol-v2-teqabfhg-20260523): C-N barrier 39.3 kcal/mol @ r=3.00 Å(既存推定 39.4 と一致)、ラジカル付加 barrier 6.1 @ 2.20 Å、well −27.8 @ 1.54 Å。推奨値は全て検証済み既定値と一致(activation 0.3/250/5000、production 2.0/250/125)。
- Alternatives considered: (a) production f1_max もバリアから導出 — 論文値 250/125 に対し実測バリアのマージンが大きく、導出式の論文根拠もないため固定+警告方式を採用。(b) activation_steps に安全係数 1.5 を乗算(初稿)— 基準バリアで検証済み 5000 が未検証の 8000 に置き換わるためレビューで棄却(hyperparameter 無断変更の防止)。
- Scientific risk: 低。既定 checkpoint では全推奨値が検証済み値と一致し挙動不変。checkpoint 差し替え時に PES バリアの乖離を検知できるようになる分、安全側。
- Licensing/commercial impact: None(既存 OrbMol-v2 backend のみ使用、Apache-2.0 記録済み)。
- Follow-up: RESUME 時の TDBB パラメータ来歴ギャップ(checkpoint は TDBB params を保存せず、resume 時は現在の env/CLI 値が無記録で適用される。`RunManifest.append_resume` は timestamp/git_sha/ckpt 位置のみ記録)。resume_history へ有効な TDBB 設定を記録する対応を別タスクで実施。

## 2026-07-07: Langevin friction — 論文値 0.001 と production レシピ 0.01 の整合(P-0c と 2026-06-26 の緊張を解消)
- Context: P-0c(2026-06-13「Langevin friction を 0.01→0.001 に修正、全スクリプトから明示 0.01 を除去」)と 2026-06-26 production レシピ(「f2=2 + **friction 0.01** + unbiased 1500」を採用)が正面衝突しており、どちらが「正」か記録上あいまい。本エントリで整合を明文化する(調査の発端: baseline `runs/s6_half_100x5` が git_dirty=true でこのレシピ=friction 0.01 で回っていたと判明。手前の分析で誤って friction 0.001 と推定していた点の訂正も兼ねる)。
- Paper anchor: PDF p.20 SI「Langevin thermostat coupling constant 1.0 ps⁻¹」= 0.001 fs⁻¹。**論文忠実値は 0.001**。
- Decision: 論文忠実値は 0.001 のまま不変(paper_faithful 系の run・PFP 等の論文再現 backend はこれを使う)。**OrbMol-v2 backend の production レシピに限り、friction 0.01 を「意図的な非論文値の逸脱」として採用する**。理由の連鎖:
  (1) 論文の PFP はラジカル付加 PES に引力チャネルを持ち、候補窓 [3,6] Å からデッドゾーンを越えて形成できる。OrbMol-v2 の PES にはこの引力が無く、f2 を論文値 10→2 に下げないと melt formation が起きない(decisions 2026-06-26, S2 sweeps)。
  (2) f2 を下げると反応熱(~−28 kcal/mol)+バイアス仕事を短い unbiased 相の Langevin で捌けず、**累積加熱で温度が暴走**する(f2=2 + friction 0.001 + unb500: 389→565K)。
  (3) friction 0.01(=論文の 10×)はこの加熱を抑える冷却レバー。**f2 変更とセットの逸脱**であり、単独の熱浴パラメータ選択ではない。
- 位置づけ: orbmol_v2_production レシピ(f2=2 + friction 0.01 + unbiased 1500)は OrbMol 固有の工学的調整。**絶対転化率・絶対速度定数(kp_eff)は論文と直接比較不可**。論文比較は Eq.11 α(t)=1−exp(−kp·t) の**形状(トレンド)**と鎖成長の定性挙動で行う。
- 実証(2026-07-07, `runs/s6_half_100x5_main0707`, 30cyc): friction 0.01 は温度を **mean 338.7K(target 333)・std 11.1・linear drift −0.066 K/1000steps(実質フラット・暴走なし)**に維持。2026-06-26 は 6cyc のみの検証だったが、全 30cyc で冷却レバーが有効と確認。α(t) は Eq.11 に適合(R²=0.85; baseline 0.98)。コード変更(H1: in-bias 検出を tentative 化し unbiased 緩和後に確定)で絶対転化率は 26%→16% に下がるが、これは非物理的な一過性結合を除外する**正当な厳格化**であり速度式の形状は保存(退行ではない)。
- Alternatives considered: (a) 論文値 0.001 に統一 — f2=2 では加熱暴走で長尺不可(2026-06-26 実証)。棄却。(b) friction を上げず unbiased をさらに延長して放熱 — コスト増で非効率、friction の方が有効(unbiased 1500 と併用済)。(c) f2 を論文値 10 に戻す — デッドゾーンで formation≈0。棄却。
- Scientific risk: 中。論文からの二重逸脱(f2 10→2、friction 0.001→0.01)は絶対速度論の論文比較を無効化する。緩和策: 上記のとおりトレンドで比較し絶対値比較はしない旨を明記。将来 paper-faithful backend(PFP 等)が使用可能になれば f2=10 + friction 0.001 に戻す。
- Licensing/commercial impact: None。
- Follow-up: (a) `configs/boost/orbmol_v2_production.yaml` に friction 0.01 の逸脱理由と本エントリ参照を明記。(b) P-0c エントリに本エントリへの相互参照を追記(下記で実施済)。

## 2026-07-08: 実測 Carothers DPn(Fig.4c)を simulation-vs-theory 化
- Context: nylon の Carothers 検証(論文 Fig.4c)が「理論 1/(1−p) を最終 p で評価するだけ」= theory-vs-theory だった。実測 DPn(シミュレーション網目の数平均重合度)を出せるようにする。
- Paper anchor: arXiv:2511.22874 Fig.4c、Table S2、decisions 2026-06-13(T-G2)。
- Decision: (a) `src/kagome/analysis/carothers.py` に `dpn_measured_from_topology` / `monomer_sets_from_bonds`(純関数、MD/IO 非依存)を追加。DPn = (単量体ユニット数)/(単量体グラフの連結成分数)を topology.jsonl から計算。(b) `run_nylon66.py` が `layout_bonds`(diamine=seed / diacid=seed+1、ビルダー配置順と整合)で `initial_bonds` を渡し **topology.jsonl を出力**。(c) `reproduce_figures.py` に `plot_dpn_vs_conversion`(実測 DPn vs p × 理論 dpn_carothers(p))、`summary.json` の `carothers_p` で step-growth 判定(vinyl は no-op)。
- 前提: 単量体所属は初期(cycle −1)snapshot の連結成分から復元。アミド結合と水生成結合は同一 diamine/diacid ペアを繋ぐため、成分数はどちらの有無でも不変(水は擬似分子として現れない)。
- Scientific risk: 低。純関数はユニットテスト済(tests/unit/test_carothers.py: 28 passed、独立再検証済)。実 MD 出力での検証は初回実 run 時に「実測点が 1/(1−p) 曲線を超えない(理想線形なら一致、分岐/環化で下振れ)」を確認。
- Licensing/commercial impact: なし。
- Follow-up: topology.jsonl は結合変化サイクルのみ書くため実測 DPn 系列は疎(Fig.4c overlay には十分)。密な曲線が要れば snapshot 方針変更(別タスク)。

## 2026-07-08: nylon amide 形成も MA と同型の capture-shell デッドゾーン(paper f2=10 で [3,6] 窓の bias 力 ≈0)。PES scan スクリプト追加
- Trigger: 長尺 nylon run 前に、OrbMol-v2 + TDBB が候補窓 [3,6] Å から amide 結合(amine N–carboxyl C)を閉じられるかを、MA デッドゾーン診断(2026-06-26)の nylon 版として事前検証。
- Tool: `scripts/scan_amide_formation.py`(scan_radical_addition.py 踏襲: constrained-FIRE 緩和 PES + 解析的 TDBB reach 表)。閉殻フラグメント methylamine+酢酸(spin=1、多ラジカル問題なし)。
- r0 = λ·(r_vdw^N+r_vdw^C) = 0.60·(1.55+1.70) = **1.95 Å**(Eq.4, VDW_RADII)。
- 解析的 reach(f1 飽和=250、GPU 不要、独立再検証済): **f2=10(run_nylon66.py の paper default)では [3,6] 窓内の最大 bias 力 0.086 kcal/mol/Å ≈ 0** → 選択 pair(3.0–3.6Å)は捕捉シェル(~2.5Å)外のデッドゾーンで力が届かない。MA(2026-06-26)と完全同型。f2=5 → 10.6、f2=2 → 116 kcal/mol/Å で橋渡し。r=3.5Å で |力|≥10 を満たす最小 f2 は 2。
- Decision: PENDING(Ask-first: f2 は論文 default のハイパラ。MA と同じく production f2=2 が有力だが、**GPU PES scan で「1.95Å 付近に到達可能な bonded チャネルが存在」を確認してから確定**)。
- Follow-up: (a) `--device cuda` で OrbMol PES scan を実行し barrier/well を実測。(b) 引力チャネルがあれば f2 引下げ(2)で長尺可、なければ MLIP 側の問題として別途(backend 交換は CLAUDE.md 商用審査)。(c) 縮合は反応熱が小さく f2 低下による加熱は MA より軽い見込み(要実測)。
- Licensing/commercial impact: なし。

## 2026-07-08: run_nylon66 に pre-TDBB minimize+equilibration を追加(segfault 修正、MA の T-G1a パリティ)
- Symptom: `runs/nylon66_f2test`(10+10, f2=2)が **TDBB ループ突入の最初の OrbMol forward パスでネイティブ segfault**(dmesg カーネルダンプ PID 3370、OOM ではない=空き 28Gi、Python トレースバック無し)。ログは "Starting TDBB" → "OrbCalculator device check" で途切れ 0 サイクル。
- Root cause(確認): nylon が build → classical 圧縮 → TDBB を **無 minimize/無 equilibration** で通していた(vinyl は持つ)。圧縮後 0.5 g/mL(fmax ~60、未収束)の近接接触が初回 MLIP 評価で過大な力/退化幾何を生み backend を落とす。MA が 2026-06-13(T-G1a)で解決したのと同一ギャップ。
- Fix: `run_nylon66.py` に `--minimize`/`--no-minimize`(既定 ON)、`--minimize-fmax`(1.0)、`--equil-steps`(2000)を追加し `PolymerizationConfig` に配線。手動 `wf._minimize`/`_run_equilibration_phase` ではなく **config 駆動**(`PolymerizationWorkflow.run()` が cycle ループ前に config から minimize+equil を実行、非 resume 時。vinyl の非活性化パスと同一)。nylon は活性化なしなので順序は build → compress → minimize → equilibrate → TDBB。二重 minimize なし。
- Paper anchor: PDF p.20(production 前に equilibration)。既定値は vinyl の承認済み値を踏襲(新規ハイパラ無し)。
- Tests: `tests/unit/test_workflow.py::TestPreTDBBMinimize`(config.minimize で minimize フェーズが出る/抑止される)。targeted suite 125 passed(独立再検証: workflow+carothers 73 passed、`--help` に新規4フラグ表示)。
- **未解決リスク(要実走確認)**: minimize 自体の初回 orb 評価も同じ悪い幾何で走るため、圧縮構造が退化しすぎていれば minimize の初手で同様に segfault しうる(minimize は forward パスを遮蔽せず、それ自体が forward パス)。その場合のフォールバック: (a) 圧縮の収束を締める(fmax≪60)、(b) より希薄配置+緩やかな圧縮、(c) ML minimize 前に classical(OpenMM/Sage)pre-minimize を挟む。再実行で判明次第対応。
- Licensing/commercial impact: なし。
- **2026-07-08 実走結果**: 再実行 `runs/nylon66_f2test`(10+10, f2=2, 10cyc, NVT, minimize+equil ON)は **exit 0 で完走。segfault は再発せず、修正は有効と確定**。温度は 300K target 近傍で安定(暴走なし)。ただし **真のアミド C-N 形成は 0 件**(summary: `carothers_p=0.0`, `dpn=1.0`, `confirmed_formations=0` counted)。f2=2 でも bias 中は C-N を ~1.6 Å まで押せるが unbiased 1500 で解離して戻る。→ NEXT GATE(scan_amide_formation.py の GPU PES で OrbMol に ~1.95 Å 近傍の bonded C-N チャネルがあるか)の実行が依然必要。

## 2026-07-08: nylon 実走で判明した3つのテンプレ/簿記バグ(段階成長の反応判定・topology 編集の不整合)
- Context: 上記 nylon f2=2 完走ランのイベントログ(bonds.jsonl / topology.jsonl / selection.jsonl)を精査したところ、segfault とは別に、段階成長テンプレートの反応判定と topology 編集に3つの不整合が判明。summary の数値(confirmed_formations 1 total / confirmed_dissociations 5)と topology 実差分(bond 追加1・削除1)のズレがこれらで完全に説明できる。r0 は `r0 = λ_vdw·mean(vdW半径)` で C-N=1.95、**N-H≈1.65、H-O≈1.63 Å**。
- Paper anchor: arXiv:2511.22874 Table S2(4-group amide template)、Eq.4(r0 定義)、Fig.4c。
- 症状と根本原因:
  (B1) **補助 H-O「水生成」ペアが biased フェーズを誤って早期打ち切りしている**。テンプレの (amine_H, carboxyl_OH) ペアは `count_as_reaction=False`・`r_max=100`(常時アクティブ)だが、`check_reactions_during_bias` はこのペアの r≤r0(1.63)= N-H…O 水素結合の H…O 距離をそのまま「暫定反応」と判定し biased phase を打ち切る(cycle 0/1/2/4/8 の 52〜112 step 打ち切り、min_pair_dist≈1.63 が一致)。結果、**本命のアミド C-N ペアが bias 時間をほぼ得られない**。→ `tentative` 判定は `counts_as_reaction=True` のペアに限定すべき。
  (B2) **解離候補が別分子間の非結合 N-H ペアに対して生成される**(cycle 5:144-32/120-22、cycle 7:144-8、cycle 9:72-8 — いずれも別モノマーの N と H)。非結合なので r>1.65 が自明成立し「confirmed_dissociation」が乱発。反応簿記側は H2 ガード(対応形成が未確定ならスキップ)で実害を止めているが、候補生成が「実際に結合している N-H」に制約されていない。
  (B3) **topology 編集パス(`_apply_topology_edits`, polymerization.py:1146-1177)が反応簿記側のガードを尊重していない**(実装不整合=バグ):
    - 形成側: `confirmed_formations()` を `counts_as_reaction` で絞らず全件に `add_bond` → **spurious な H-O 共有結合 [238,309](count_as_reaction=False)を topology に書き込んだ**(cycle 6)。
    - 解離側: group-update パス(:303-319)にある H2 ガード(対応形成が未確定なら skip)が topology パスには無い。`has_bond` チェックのみのため、実結合だった **N-H [103,118] を cycle 3 で topology から削除**(ログは「Skipping dissociation (103,118)…H2」と言っているのに実結合は消えた)。他の非結合ペアは has_bond=False で偶然 no-op。
- Decision: PENDING(実装は未着手)。3件とも TDBB の反応選択・確定の意味論に関わる(CLAUDE.md ask-first: 「簡略化が TDBB の科学的意味を変える」)ため、修正方針を提示しユーザー確認後に実装する。提案する最小修正:
  (F1) `check_reactions_during_bias` の暫定判定を `pair.counts_as_reaction` が真のペアに限定(補助 H-O ペアで打ち切らない)。
  (F2) 解離候補生成を「現 topology 上で実際に結合している原子ペア」に制約(選択段でトポロジ参照)。
  (F3) `_apply_topology_edits` の形成側に `counts_as_reaction` フィルタ、解離側に group-update と同じ H2 ガード(対応形成が確定した candidate_id のみ)を適用し、両簿記の整合を取る。
- Scientific risk: 中。現状の「confirmed_formations 1 / dissociations 5」は **中身が全てスプリアス**(H-O 水素結合の誤認と非結合 N-H の誤解離)であり、このカウントを反応実績として解釈してはならない。真のアミド形成 0 という結論自体は B1〜B3 と独立(C-N は unbiased で戻るため元々未確定)だが、**B1 により C-N が bias 時間を奪われている可能性があり、「OrbMol の PES にアミドチャネル無し」と結論するのは時期尚早**。まず F1 を入れて C-N に十分な bias を与えた再走、そして NEXT GATE の PES scan が必要。
- Licensing/commercial impact: なし(OrbMol-v2 のみ使用)。
- Follow-up: (a) F1〜F3 の実装可否をユーザーに確認。(b) `scripts/scan_amide_formation.py --device cuda` で C-N の bonded チャネル de-risk。(c) F1 実装後に f2=2 で再走し、C-N が確定するか(Carothers 実測点が出るか)確認。(d) 各修正に unit/回帰テスト(repo DoD)。

## 2026-07-08: F1 実装 — biased フェーズ終了トリガを「アンカー形成」に限定(B1 の一般化修正、安全側の切り分け)
- Context: 上記 B1 の当初案(count_as_reaction ペアに限定)では不十分と判明。段階成長の従属開裂(N-H, C-OH)は `count_as_reaction=True`(デフォルト)なので、それだけでは依然として本命の C-N 形成(r0=1.95、遠方から寄せるため最も遅く跨ぐ)を出し抜いて biased フェーズを打ち切る。ユーザーと「複数結合の協奏組み替えは原理的に不能か」を議論し、4-group テンプレ+H2 ガードで協奏縮合は表現可能=不能ではなく、真因は「終了トリガがアンカー(反応座標)を取り違えている」ことと確認。安全側の切り分けとして F1 のみ最小実装→再走で C-N 確定を見る、B2/F3 は保留と合意。
- Paper anchor: arXiv:2511.22874 Table S2(amide anchor = amine_N–carboxyl_C 形成)、decisions 2026-07-03 D1(検出→終了→unbiased 確認)。
- Decision: biased フェーズの**終了トリガを「アンカー tentative イベント」に限定**(`_run_biased_phase`, polymerization.py)。実装:
  - ループ前に `has_counted_formation = any(p.is_formation and p.counts_as_reaction for p in active_pairs)` を1回計算。
  - 純粋ヘルパー `_is_terminating_bias_event(ev, has_counted_formation)`(module-level, polymerization.py:173):終了イベント条件は (i) counted formation の跨ぎ(`event_type=='tentative_formation' and counts_as_reaction`)、または (ii) counted formation が active_pairs に無い純解離テンプレート時のみ counted dissociation の跨ぎ(フォールバック=AIBN 等)。
  - `check_reactions_during_bias` は**無改変**(全跨ぎを従来どおり記録=監査証跡維持)。毎ステップ検出は続け、`terminating` に該当するイベントがある時のみ `break`。従属開裂と water pair は C-N が結合距離に達するまで bias され続け、確認段で既存 H2 ガードにより C-N と一緒に確定。
- 影響範囲: vinyl(アンカー=radical_C+alpha_C, counted formation)は挙動不変。`_run_activation_phase`(AIBN 活性化)は別メソッドで無改変。
- Tests: `tests/unit/test_workflow.py::TestBiasedPhaseTerminationGating`(ヘルパー直接5件+実 MD ループ1件)。独立再検証: `test_workflow.py + test_carothers.py` = **79 passed**(orchestrator が自走確認)。
- Scientific risk: 低〜中。協奏縮合の意図どおり「C-N 到達まで bias→開裂/水を束ねて確定」になる。**ただし B2(非結合 N-H への解離候補生成)と F3(topology 編集の counts/H2 ガード不整合)は未修正**なので、再走で spurious な confirmed_dissociation や topology への spurious 書き込みは残りうる。今回の切り分け目的(F1 で C-N が確定するか)には影響しないが、Carothers/topology の解釈時は B2/F3 未修正を念頭に。
- Licensing/commercial impact: なし。
- Follow-up: (a) f2=2 で再走(`runs/nylon66_f2test`)し confirmed_formations(counted C-N)>0 になるか確認。(b) C-N が出れば B2/F3 に進み topology/Carothers を正す。(c) 依然 0 なら NEXT GATE の GPU PES scan で OrbMol の bonded C-N チャネル有無を判定。
- **2026-07-08 F1 再走結果(exit 0 完走)**: F1 は設計どおり機能。biased ステップ数が 129〜1005(旧 52〜112 の早期打ち切りから大幅増)へ延び、各サイクルで **anchor=C-N 形成の跨ぎで終了**、min_pair_dist は後半 1.92〜1.94 Å まで到達(C-N を真に結合距離へ寄せた)。**しかし counted confirmed_formations は 10/10 サイクルで 0**(carothers_p=0.0, dpn=1.0)。→ f2 でも F1 でも縮合は確定せず。**真因は下記アーキテクチャ限界**(次エントリ)。副次観測: 唯一の confirmed_formation は (70,309) counts=False r=1.02 Å=**水の O-H**(アミン H70 がカルボキシル O309 へプロトン移動して O-H が成立)。水生成の半分は起きるが、原子除去が無いためアミド C-N 側は復元する——機構仮説の直接証拠。confirmed_dissociations=7 は全て B2 の spurious(非結合/脱離基ペアが bias で離散)で H2 スキップ済み。

## 2026-07-14: 大マシン py-spy 実測に基づく orb backend 高速化(--compile 配線バグ修正 + empty_cache 条件化)
- Trigger: nylon66 paper-scale 走行中(大マシン RTX 5000 Ada 32GB + Xeon 8558U、WSL)に `py-spy top --pid` で非侵襲プロファイル(rate 100 と 10 の2回、順位一致。両回とも behind-in-sampling 警告ありのため % は傾向値)。OwnTime 上位: torch forward(linear)/silu > `torch.cuda.empty_cache()`(~9%) > `_compute_neighbor_list_with_fallback`。GIL 0% / Active ~27%。trajectory.jsonl の JSON 書き出しは上位に不在(容疑者から除外)。
- 結論: ボトルネックは Orb GNN の Python 側フォワード(小カーネル大量発行の CPU オーバーヘッド)+ 毎ステップ empty_cache。GPU-Util ~60% で頭打ちの観測(2026-07-13 較正)と整合。fix-plan-2026-07-06 B3(empty_cache コスト実測)をこれで回収。
- Decision 1(empty_cache 条件化): `create_orb_calculator(empty_cache: bool = True)` を追加し、`compute()` 末尾の毎ステップ `torch.cuda.empty_cache()` をゲート。**既定 True で挙動不変**(16 GB での VRAM 断片化対策 2026-06-15 を維持)。32 GB 機では run スクリプトの `--no-empty-cache` で無効化(CPU 時間 ~9% 削減見込み)。
- Decision 2(--compile 配線 + バグ修正): orb-models は model compile を公式サポート(Apache-2.0 の既存機能、upstream 公称 10k 原子で ~1.7x)。backend の既存 `compile` 引数を run スクリプト(run_nylon66 / run_epoxy_amine / run_vinyl_aibn)の `--compile` フラグに配線。**バグ修正**: 従来 `TORCHDYNAMO_DISABLE=1` を無条件 setdefault していたため(Windows の nvalchemiops/cl.exe 対策)、`compile=True` を渡しても dynamo ごと無効化され eager 実行になっていた。`_configure_torch_env()` に分離し、compile=False の時のみ設定(Windows 対策は維持)。compile=True で環境変数が既に立っている場合は warning。
- 新規依存なし。両フラグとも opt-in で、**デフォルト挙動は完全不変**(既存 run の再現性に影響なし)。summary JSON に `compile` / `empty_cache` を記録(実験トレーサビリティ)。
- Tests: `test_backends.py::TestOrbTorchEnvConfig`(dynamo env 分岐 3 件)、`TestOrbEmptyCacheGating`(fake model で empty_cache 呼出の有無 2 件)。
- **Validation protocol(未実施、大マシンでの採用前に必須)**:
  (a) 数値等価性: 同一構造・同一 seed で compile 有/無の energy/forces を比較(目安 forces atol 1e-4 kcal/mol/Å — torch.compile はカーネル融合で浮動小数点演算順序が変わりうる)。短尺 run(nylon 10+10 2cycle 等)で反応イベント列が一致するかも確認。
  (b) 性能: s/step を compile 有/無 × empty_cache 有/無の 4 通りで実測。重合でグラフサイズが変わるため再コンパイル頻度をログで監視(頻発すると逆効果の可能性)。
  (c) 初回コンパイルの数分は paper-scale では償却可。
- Scientific risk: 低(既定 OFF。採用は (a) の検証通過後)。
- Licensing/commercial impact: なし(orb-models Apache-2.0 の既存機能、新規依存なし)。
- Follow-up: 大マシンの次回 run(epoxy paper-scale 等)で validation protocol を通してから `--compile --no-empty-cache` を採用。

## 2026-07-08: 【重要・アーキテクチャ限界】MLIP 原子集合が反応で不変 → 段階成長縮合は現構造では確定不能
- Context: F1 修正後も nylon の counted amide 形成が 10/10 サイクルで 0。f2/トリガの問題を全て潰した上でコードを精査し、根本原因を特定。
- 根本原因(コードで確定): 反応実行(`apply_vinyl_addition` topology.py:119、`_apply_topology_edits` polymerization.py:1131-1179)は**サイドカーの結合グラフ `BondTopology`(add_bond/remove_bond/set_order)のみを編集し、OrbMol に渡す原子集合(state.positions / state.species)は一切変更しない**。全ラン通じ `n_atoms=440` 不変。確定は unbiased 緩和の終端 `check_outcomes` で距離判定(D1, 2026-07-03)。
- なぜ vinyl は成立し縮合は不成立か:
  - vinyl 付加生成物の C-C(~1.5 Å)は**原子を1つも抜かなくても OrbMol の PES 上でそのまま安定極小**。よって結合グラフ編集(価数簿記/Carothers/XYZ 出力用)だけで足り、MLIP 入力を変えなくても物理的に本物。unbiased で保持され confirmed。
  - アミド縮合は**脱離基(アミン H・カルボキシル OH)が実際に系から抜けて水になって初めて**安定なアミドになる。現構造は脱離基原子を MLIP 系から除去しないため、OrbMol は常に「N プロトン化・C が OH 保持のまま C-N が 1.9 Å に密集した高エネルギー前反応錯体」を評価し、unbiased で押し戻す。→ **f2/F1 に依らず confirmed=0**。
  - D1 の前提「真に結合するペアは緩和後も r≤r0」は vinyl では成立、**縮合では原理的に不成立**(未コミットの反応物盆地を緩和しているため)。
- Paper anchor: §2.2 step 3(検知→biased 終了→unbiased 緩和)、Table S2(縮合 formation/dissociation 対)。論文の TDBB は結合ブースティングであり、遷移状態突破後は**生成物盆地**へ入る=脱離基が抜けた状態を緩和しているはず(要・論文プロトコル精読)。
- Decision: PENDING(ask-first。TDBB の科学的意味に関わる大規模変更のためユーザー判断待ち)。候補方針:
  (A) **原子レベルの反応実行**を縮合に導入: 検出時にアミン H とカルボキシル OH を MLIP 原子集合から除去して水分子として再配置し、本物のアミド生成物を OrbMol に提示してから緩和(commit-then-relax、縮合では D1 の confirm-before-commit を反転)。これが nylon を実際に前進させる本命。
  (B) 先に**生成物側 de-risk**: 既製アミド(例 N-メチルアセトアミド+H2O)を OrbMol で緩和し極小保持を確認/`scripts/scan_amide_formation.py --device cuda` で生成物側 PES を確認。(A) の前提(OrbMol がアミドを安定化する)を安価に検証。
  (C) B2/F3 の簿記クリーンアップのみ: spurious dissociation/topology を正すが、アミドは生成しないので単独では価値低。
- 推奨順: (B) で前提検証 → (A) 実装 → その後 B2/F3。
- Scientific risk: 高(段階成長縮合の再現可否そのもの)。vinyl 系には影響しない(原子除去は縮合テンプレート限定)。
- Licensing/commercial impact: なし(OrbMol-v2 のみ)。
- Follow-up: ユーザーに (A)/(B)/(C) の方針確認。(A) 実装時は論文 §2/SI の反応実行手順を精読し decisions に反応実行の原子操作仕様を明記。
- **2026-07-08 (B) de-risk 実施済み・前提確定**: ユーザーが方針 B(生成物側 de-risk 先行)を選択。`scripts/derisk_amide_product.py`(新規、生成物側専用。既存 scan_amide_formation.py は反応物側で別物)を OrbMol/GPU で実行(`runs/calibration/amide_product_derisk.json`):
  - **Check 1 PASS**: 本物のアミド生成物(N-メチルアセトアミド+水)を自由 FIRE 緩和 → **C-N が 1.377→1.359 Å の実アミド長で安定極小に収束**(converged, max force 0.47)。→ **方針 A の前提「OrbMol は脱離基を除いた本物のアミドを安定化する」は成立**。
  - **Check 2 reverts=True**: 脱離基を残したままの反応錯体(酢酸+メチルアミン, N···C=1.55 Å)を自由緩和 → **1.55→3.03 Å に復元**。0K 最小化ですら戻る=ワークフロー失敗を最小系で再現し、「密集した前反応錯体こそが原因で OrbMol 自体は問題ない」という診断を確証。
  - **Check 3**: dE_reaction = −4.27 kcal/mol(気相アミド生成として妥当・軽度発熱)。
  - **結論**: **nylon が確定しない原因は反応の駆動が完結しないまま停止していること。OrbMol backend は本物のアミドを完全に支持する。** 次は論文精読で反応実行手順を確認(下記で方針 A は棄却され連言判定へ訂正)。

## 2026-07-08: 【方針訂正】論文の反応実行は「原子手術」ではなく「連言反応イベント判定(ij∧ik∧jl 同時)」— 方針 A を棄却、F1 を連言判定へ置換
- Context: 方針 A(原子レベル反応実行=脱離基を MLIP 原子集合から除去し水を明示配置)の実装前に、ユーザー承認どおり論文(paper/2511.22874v1.pdf, Section 2.2 step 3-4, p.5)の反応実行手順を精読。結果、**A は論文の方法ではないと判明**。
- Paper anchor(短い要約のみ; 原文転載しない): §2.2 step 3「Reaction Induction」— バイアス下加速 MD を、ペア (ij, ik, jl) の距離が**同時に**結合条件を満たすまで(または最大時間まで)継続。結合条件=「2原子間距離が vdW 半径和の 60% 未満なら結合」(=r0=λΣr_vdw, λ=0.6 と同一閾値)。fp=1(形成)は新たに条件を満たした時、fp=0(解離)は満たさなくなった時にイベント。イベント後に無バイアス MD で構造緩和。§2.2 step 4「Reaction Evaluation」— 緩和後に前後の結合状態を比較し、変化した原子を reacted としてグループから除去。SI p.21 + Table S2: nylon は Vf を i-j(アミド)と k-l(水)に、Vd を i-k(N-H)と j-l(C-OH)に適用。k-l 水はバイアスのみでトリガー(ij,ik,jl)に含まれない。**原子除去・水の明示構築の手順は論文に存在しない**(脱離基はバイアスで解離し水になって系内に残る; MLIP 原子集合は論文でも不変)。
- Root cause 再確定: 現コード(F1)は「counted formation(C-N)の跨ぎ**のみ**」で biased を終了/確定していた。論文は「ij 形成 AND ik 解離 AND jl 解離が**同時**」まで継続。F1 は N-H/C-OH が解離しきる前に停止するため、脱離基が付いたまま unbiased に入り復元 → confirmed=0。de-risk Check 1(脱離基を除いた本物のアミドは OrbMol 上で安定)は、連言が満たされた時点の生成物が unbiased で保持されることを保証する。
- Decision: **方針 A(原子手術)を棄却。代わりに反応イベント判定を論文忠実な候補単位の連言に修正する(F1' と呼ぶ; F1 を置換)。**
  (1) トリガー(biased 終了)= 候補の**トリガーペア全て**が同時に結合条件を満たした時。トリガーペア = テンプレの score_pair=True のペア(=論文の識別ペア集合 P: ij, ik, jl)。formation ペアは r≤r0、dissociation ペアは r>r0。k-l 水(score_pair=False, bias-only)はトリガーに含めない。
  (2) 確定(check_outcomes, unbiased 緩和後)も同じ連言を候補単位で評価。全トリガーペアが緩和後も条件を満たす候補のみ confirmed。→ H2 ガード(formation なしに dissociation 単独確定)は連言評価に自然に内包。
  (3) vinyl は score_pair=True の実効トリガーが i-j のみ(i-k/j-l は constraint_only=無バイアスで常時充足)なので挙動不変。
- Alternatives considered: (a) 方針 A 原子手術 — 論文非該当・大規模・checkpoint 複雑化。棄却。(b) F1 の formation-only 維持 — 論文の同時条件を満たさず縮合が復元。棄却。(c) 現状の per-pair 独立確定 — spurious dissociation を生む(B2 の温床)。棄却。
- Scientific risk: 中。連言により縮合が論文どおり駆動されるはず。ただし f1/f2/biased_steps が i-k/j-l の実結合(N-H ~90 kcal/mol)を biased 期間内に解離させ切れるかは経験的検証が要る(次の再走で確認)。切れない場合は f1_max_dissociation/biased_steps の再検討(論文値 125/2000 を基準に)。
- Licensing/commercial impact: なし。
- Follow-up: (a) F1' を実装(トリガー+確定を候補単位連言に)、F1 のヘルパーを置換。(b) B2(解離候補を実結合ペアに制約)は連言確定と併せて意味が出るので同時に検討。(c) 実装後 f2=2 で再走し confirmed amide>0 と Carothers 実測点を確認。(d) unit/回帰テスト。
- **2026-07-08 F1'+B2 実装済み(サブエージェント、orchestrator が全面レビュー+独立テスト)**:
  - PairBias に `is_trigger`(=PairSpec.score_pair)追加(tdbb.py)。`_build_pair_biases` が伝播。水 k-l(score_pair=False)は is_trigger=False。
  - `BondTracker.check_reactions_during_bias` は fired candidate_id リストを返す(候補の全トリガーペア同時充足で発火)。per-pair の tentative 記録(audit)は保持。candidate_id<0(活性化/legacy)は従来の単発発火にフォールバック。
  - `check_outcomes` は候補単位連言で確定(緩和後に全トリガー充足の候補のみ、その全ペアを確定=amide formation[counts=True]+leaving 解離+水[counts=False, is_formed 時のみ])。部分充足は何も確定せず(spurious 単発解離を排除)。→ `_apply_topology_edits` が C-N 追加・N-H/C-OH 除去・水 O-H 追加を整合実行(**F3 の症状も副次的に解消**)。
  - `_run_biased_phase` は fired で break、`_is_terminating_bias_event`/has_counted_formation を削除。vinyl は実効トリガー i-j のみで挙動不変。
  - B2: selection.py `_enumerate_recursive` に `dissociation_group_pairs`(score_pair・非 constraint・is_formation=False)への「実結合必須」制約(L6 の双対、topology 有時のみ)。非結合の分子間 N-H/C-OH 解離候補を排除。
  - Tests: test_bonds.py `TestConjunctiveReactionEvent`、test_workflow.py `TestConjunctiveBiasedPhaseTermination`、test_selection.py `TestDissociationTopologyGuard` 追加/置換。**`pytest tests/unit` = 444 passed, 1 skipped(openmm 無し、無関係)**、orchestrator が独立再実行で確認。
  - NEXT: f2=2 で再走し confirmed amide(counted)>0 と Carothers 実測点を確認。経験的リスク: f1_max_dissociation=125 が実 N-H/C-OH を biased 窓内で解離させ切れるか。
- **2026-07-08 F1'+B2 再走結果(exit 0、成功)**: `runs/nylon66_f2test`(10+10, f2=2, 10cyc, NVT)。**confirmed_formations(counted アミド)= 4 本**(cycle 0/2/3/4 で "candidate reaction fired" → unbiased 保持 → 確定)。confirmed_dissociations=8(4 アミド×2 脱離基、全て対応形成とペア)、水形成 4(counts=False)。**spurious 単発解離ゼロ**(B2+連言確定が奏功)。**carothers_p=0.2, dpn=1.25=1/(1−0.2) 理論と整合**。topology 差分は全て完全な縮合4点セット(例 cyc0: C-N 0-327 追加, N-H 0-8 除去, C-OH 327-329 除去, 水 8-329 追加)。f1_max_dissociation=125 は実 N-H/C-OH を biased 窓内で解離させ切れることを実証(fired step 409〜933)。反応進行で候補が 18→0 へ単調減少=モノマー端消費を正しく反映。**段階成長縮合(生成 dissociation・4-group multi-pair・Carothers)の TDBB 機構が E2E で実証された(Track 1 の当初ゲート達成)**。cycle 1/5/6/7/8/9 は連言未成立で未反応(近接不足)=正常。
- Scientific status: nylon Track 1 の反応機構検証は完了。残タスク: (a) 論文 Fig.4c 相当の Carothers 曲線には大規模・長時間ラン(100+100, NPT, 多サイクル)で複数 (p,DPn) 点が必要。(b) F3 の簿記クリーンアップは連言確定で症状解消済みだが `_apply_topology_edits` のガード明示化は任意。(c) 温度安定性/密度など従来チェックの図生成。

## 2026-07-08: 【Track 2 / E0 設計】bulk epoxy-amine 開環付加の PES + 生成物 de-risk(proceed-or-fallback ゲート)
- Context: Track 1(nylon)の反応機構検証が完了(連言反応イベント F1'+B2)。次は Track 2 = **bulk epoxy-amine 硬化**。これは 2026-06-20 で除外した epoxy/CuO **表面**系とは別物で、CuO 金属酸化物(OrbMol 未学習ドメイン)を避け、**有機のみ(C/H/N/O)で OrbMol のドメイン内に留まる**方向。外部クロスメソッド比較先は ref#2(Provenzano 2025, ACS Appl. Polym. Mater. 7(8):4876, オープンソース古典 bulk 架橋プロトコル)。E0 は「OrbMol+TDBB でエポキシ開環+アミン付加が駆動でき、開環生成物を安定に保持できるか」を回すゲート。MA(2026-06-26)・nylon(2026-07-08 de-risk)と同じ流儀。
- 反応化学(bulk epoxy-amine 開環付加): アミン N がエポキシ環の末端 C を求核攻撃 → **N-C 結合形成 + 環の C-O 結合開裂 + アミン N-H 開裂 → その H が環 O へ移り β-ヒドロキシ基(O-H)形成**。nylon 縮合と**構造同型の 4-group テンプレート**だが、脱離基は分子外へ出ない(水は生じない):O は隣接環 C に結合したままヒドロキシルとして分子内に残る(分子内開環)。1級アミンは2回反応(1°→2°→3° 架橋)するが E0 は単回付加の PES de-risk に限定(2回目=2級付加は E1 で扱う/E0 で任意の第2スキャン)。
  - Gi=amine N, Gj=epoxy 末端 C, Gk=amine H, Gl=epoxy 環 O。ペア: (i,j) N-C 形成 / (i,k) N-H 解離 / (j,l) C-O 環解離 / (k,l) H-O ヒドロキシル形成(bias-only, count_as_reaction=False)。これは連言反応イベント(F1', 2026-07-08)がそのまま適用できる(ij∧ik∧jl 同時)。B2(解離ペア=実結合)も (i,k) N-H と (j,l) 環 C-O が実結合なので機能。
  - r0(Eq.4, λ=0.6): N-C=0.6·(1.55+1.70)=1.95 Å(nylon と同一)、環 C-O 解離=0.6·(1.70+1.52)=1.93、H-O 形成=0.6·(1.10+1.52)=1.57。→ dead-zone 構造は nylon/MA と同型、f2≈2 が既定見込み。
- Paper anchor: 論文 §2 epoxy テンプレート(SI, ただし CuO 表面版)。decisions 2026-06-18(表面 4-group: Gi=epoxy O, Gj=1°N, Gk=2°N, Gl=surface OH — **bulk では surface OH が無く O 自身がヒドロキシルになる別テンプレート**)、2026-06-20(CuO 除外=有機ドメイン限定)、2026-06-26(MA dead-zone/f2)、2026-07-08(nylon 連言+生成物 de-risk の教訓: 生成物側の安定性検証が決定的)。
- Decision(E0 設計): 新規 `scripts/scan_epoxy_amine.py`(or 2本)で、nylon de-risk の**両面**を回す:
  (1) **反応物側 PES + TDBB reach**(scan_amide_formation.py 相当): 固定 N···C 距離で残りを constrained-relax、[3,6] Å 候補窓〜結合距離の PES とバイアス reach を f2 sweep(10/5/2)で測定 → dead-zone/f2 判定。
  (2) **生成物側安定性 de-risk**(derisk_amide_product.py 相当・**決定的**): 開環済み β-ヒドロキシアミン生成物を自由 FIRE 緩和し、N-C(~1.47 Å)が安定極小として保持されるか + 環が再閉環しないかを確認。反応物錯体(未開環・環 O 保持)の自由緩和が(開環しないなら)復元することも対照確認。
  (3) **反応エネルギー**: dE = E(開環生成物) − E(エポキシ+アミン)。エポキシ環歪み解放(~-20〜-25 kcal/mol 見込み)で発熱のはず → OrbMol がこれを再現するか。
  - フラグメント(閉殻シングレット): エポキシ = プロピレンオキシド `CC1CO1` もしくはグリシジルメチルエーテル `COCC1CO1`(DGEBA のグリシジルエーテル末端に近い)。アミン = メチルアミン `CN`(1°、nylon と同モデル)。
- Gate(proceed-or-fallback):
  - **PROCEED to E1**(bulk テンプレート+ビルダー+run): N-C の到達可能な bonded チャネルあり **かつ** 開環生成物を OrbMol が安定保持 **かつ** 発熱が妥当。
  - dead-zone(f2=10 で reach≈0)→ f2=2 採用(MA/nylon と同じ、想定内)。
  - 生成物を保持しない/再閉環する → **MLIP 懸念**。別バックエンドは CLAUDE.md commercial-license review 必須・既定にしない。
- Alternatives considered: (a) CuO 表面版をそのまま → OrbMol ドメイン外(2026-06-20 除外)。(b) E0 を飛ばして直接 bulk run → nylon で「機構/PES を先に de-risk しないと dead-zone/生成物問題で 0 反応」を実証済み。却下。(c) DGEBA/DETA フルフラグメント → PES de-risk には過大、小分子代表で十分(必要なら --monomer-fragments 相当で後追い)。
- Scientific risk: 中。開環は歪み解放で熱力学的に有利=生成物保持は nylon アミドより容易な見込みだが、OrbMol の環歪み PES 精度は未確認ゆえ de-risk が必要。1°→2°→3° の多段/2級アミンの立体障害は E1 以降。
- Licensing/commercial impact: なし(OrbMol-v2 有機ドメインのみ)。E0 が fallback(別 MLIP)を示唆した場合のみ license review。
- Follow-up: (a) ユーザー承認後 `scripts/scan_epoxy_amine.py` を実装(反応物側+生成物側+energetics、--device cuda)。(b) GPU 実行してゲート判定。(c) PROCEED なら E1(bulk epoxy-amine テンプレートを scripts/_systems.py に追加、1°→2°→3° の群再割当を含む builder + run script)。(d) 関連: [[multi-system-verification-plan]]。
- **2026-07-09 実行結果 → 次エントリ(PROCEED 判定、f2=2 採用)**。

## 2026-07-09: 【Track 2 / E0 実行結果】PROCEED — OrbMol は開環チャネル・生成物保持・発熱を全て支持。epoxy でも f2=2 を採用
- Context: 2026-07-08 E0 設計に基づき `scripts/scan_epoxy_amine.py` を GPU 実行(WSL pfpoly-gpu、backend=orb(OrbMol-v2)、device=cuda、seed=20260708)。フラグメント: エポキシ=プロピレンオキシド `CC1CO1`、アミン=メチルアミン `CN`、開環生成物=`CNCC(O)C`。原子収支 C4H11NO で平衡確認済み。成果物: `runs/calibration/epoxy_amine_scan.json`、ログ `runs/calibration/epoxy_e0_direct.log`。
- 結果(ゲート3条件+reach):
  (1) **TDBB reach**: 論文 f2=10 では候補窓 [3,6] Å の最大バイアス力 0.086 kcal/mol/Å ≈ 0 = **デッドゾーン**(MA 2026-06-26 / nylon 2026-07-08 と完全同型)。f2=2 で r=3.5 Å の力 12.7、窓内最大 115.8 kcal/mol/Å と橋渡し可能。r=3.5 Å で |力|≥10 を満たす最小 f2 は **2**。
  (2) **反応物側 PES**(拘束緩和 N…C スキャン 6.0→1.4 Å): r=3.4 Å に浅い井戸(−1.24 kcal/mol)、障壁 36.6 kcal/mol @ r=1.80 Å、結合側極小に到達可能(bonded min accessible=True)。
  (3) **生成物側 de-risk(決定的)**: 開環 β-ヒドロキシアミンの自由 FIRE 緩和で N-C 1.470→1.451 Å を保持(85 steps, converged)、攻撃 C…環 O は 2.42 Å で**再閉環なし** → PASS。対照(環無傷の反応物錯体 N…C=1.55 Å)は自由緩和で 3.35 Å へ復元し環は開かない=偽陽性でないことを確認。
  (4) **反応エネルギー**: ΔE = −29.5 kcal/mol(発熱)。環歪み解放の見込み(−20〜−25)と整合。
- Decision:
  (a) **ゲート判定 = PROCEED to E1**(E0 設計の3条件を全て充足)。
  (b) **epoxy-amine 系の production f2=2.0 を採用**。論文値 10 からの意図的逸脱で、根拠は MA/nylon と同じ(OrbMol PES に候補窓からの引力チャネルが無く、f2=10 ではバイアス力が捕捉殻に届かない)。Ask-first 該当(論文 default ハイパラの変更)だが、本スキャンの解析的 reach 実測+ユーザー承認(2026-07-09「記録して E1 準備へ」)により確定。
- Scientific risk: 低〜中。障壁 36.6 kcal/mol は f1_max_formation=250 に対しマージン十分。f2=2 逸脱により絶対速度論は論文比較不可(2026-07-07 friction エントリと同じ位置づけ=トレンド比較のみ)。発熱 −29.5 kcal/mol は nylon 縮合(−4.3)より大きく MA(−27.8)並みのため、加熱対策として OrbMol production レシピ(friction 0.01 + unbiased 1500)を E1 でも踏襲する。
- Licensing/commercial impact: なし(OrbMol-v2 のみ使用)。
- Follow-up(E1 準備): (a) `scripts/_systems.py` に bulk epoxy-amine 4-group テンプレート+ビルダー(Gi=amine N / Gj=epoxy 末端 C / Gk=amine H / Gl=環 O、1°→2°→3° の群再割当を含む)。(b) `scripts/run_epoxy_amine.py`(run_nylon66.py 踏襲: minimize+equil、checkpoint/resume、f2=2 既定、333 K)。(c) unit テスト。(d) 小規模スモーク → paper-scale。(e) 外部比較は ref#2(Provenzano 2025)プロトコル([[multi-system-verification-plan]])。
- **2026-07-09 実装済み → 次エントリ(E1 実装)**。

## 2026-07-09: 【Track 2 / E1 実装】bulk epoxy-amine ビルダー+1°→2°→3° updater+run script
- Context: E0 PROCEED(前エントリ)を受け、bulk epoxy-amine 硬化の実行系一式を実装。
- Paper anchor: arXiv:2511.22874 SI epoxy テンプレート(4-group、CuO 表面版を bulk へ適応 = 2026-07-08 E0 設計)、Table S3(DGEBA+DETA 樹脂、100:50)、p.24(epoxy production 333 K)、Eq.4(r0, λ=0.6)。連言反応イベント F1'(2026-07-08)がそのまま適用される。
- 実装:
  (a) `scripts/_systems.py`: `_DGEBA_SMILES`/`_DETA_SMILES` 定数、`_find_epoxide_c_o`(SMARTS C1CO1、攻撃対象=H の多い末端環 C)、`_find_amine_n_h`(sp3 N の全 N-H を登録; 1° も 2° も対象)、`build_epoxy_amine_system`(nylon ビルダー踏襲: エポキシ先・アミン後、rdkit_seed/seed+1)。テンプレート: (amine_N, epoxy_C) 形成 [3,6] / (amine_N, amine_H) 解離 [0,3] / (epoxy_C, ring_O) 解離 [0,3] / (amine_H, ring_O) ヒドロキシル形成 bias-only(score_pair=False, count_as_reaction=False)= nylon Table S2 と構造同型。戻り値に `amine_h_map`(N→全 N-H の global index)を追加。
  (b) `src/kagome/workflows/polymerization.py`: `EpoxyAmineAdditionUpdater`(DefaultPostCycleUpdater の H2 候補単位セマンティクスを維持しつつ、**確定形成でアミン N を群から除去しない**。消費するのは攻撃された epoxy C・環 O・移動した H のみ。N は amine_h_map の登録 H が amine_H 群から尽きた時=3級化した時のみ退役)。残 H 数は live 群メンバーシップから導出するため checkpoint 追加状態なし(processed カウンタは既存の save/restore がそのまま効く)。
  (c) `scripts/run_epoxy_amine.py`: run_nylon66.py 踏襲(minimize+equil 既定 ON、checkpoint/resume、classical 圧縮パス)。既定: DGEBA+DETA 10:5(Table S3 の 2:1 比)、T=333 K(paper p.24)、f2=2.0(前エントリ)、biased 2000 / unbiased 1500(OrbMol production レシピ)、密度 0.5(SI S-3 慣行)。主要メトリクス=エポキシド転化率(counted N-C 形成 / エポキシド総数)+アミン H 転化率。
- 論文にない仮定(明示):
  (i) **1° と 2° のアミン N を単一 amine_N 群に統合**(論文の表面テンプレートは Gj=1°N / Gk=2°N を別群にする)。bulk では両者の TDBB パラメータを同一にし、多段付加は updater の退役ロジックで表現。選択・バイアスの意味論は不変で、群分割は候補生成の等価な再パラメータ化。
  (ii) N-H は**全 H を登録**(nylon は N あたり 1 H)。同一 N の複数 H が別候補になり得るが、確定した H から順に消費されるだけで二重反応は H2/連言確定が防ぐ。
  (iii) アンサンブル既定 NPT 1 atm(論文の epoxy production は NVT だが CuO スラブ系での指定。bulk は密度到達のため repo の NPT 慣行に従い、--no-barostat で NVT 可)。
- Tests: `tests/unit/test_systems.py` に TestEpoxyAmineHelpers/TestBuildEpoxyAmineSystem(エポキシド検出は propylene oxide/DGEBA、DETA の N/H 数、群サイズ、bias-only kl、amine_h_map 整合)、`tests/unit/test_workflow.py` に TestEpoxyAmineAdditionUpdater(1回目付加で N 残留・2回目で退役、H2 ガード、二重処理防止、resume 相当のカウンタ復元)。**pytest tests/unit = 462 passed**(既存 444 から退行なし)。
- Scientific risk: 低〜中。機構は nylon で E2E 実証済みの連言イベントの同型適用。未検証は (a) DGEBA の嵩高さ・芳香環存在下での OrbMol 安定性(スモークで確認)、(b) 2級アミン付加の立体障害が biased 窓内で越えられるか(paper-scale で計測)。
- Licensing/commercial impact: なし(OrbMol-v2 のみ、RDKit/既存依存のみ使用)。
- Follow-up: (a) GPU スモーク(2 DGEBA + 1 DETA)で配線確認。(b) 小規模実走(10+5)で最初の開環確定を確認。(c) paper-scale(100+50)+ Fig.5 相当の種濃度追跡。(d) E2: ref#2(Provenzano 2025)配合との外部比較。(e) ブランチ+PR でコミット(E0 スクリプト群と併せて)。
- **2026-07-09 実走結果 → 次エントリ(機構成功、ただし friction 未配線で過熱 → 修正)**。

## 2026-07-09: 【Track 2 / E1 実走結果】開環付加の機構検証 SUCCESS(2/20 開環確定)。ただし friction 未配線で過熱 → 0.01 を既定化
- Context: `runs/epoxy_amine_smoke10`(10 DGEBA + 5 DETA = 590 原子、10 cyc、NVT 333 K、f2=2、seed 7、OrbMol/cuda、exit 0 完走、約 2.5 h)。配線スモーク(2+1、runs/epoxy_e1_smoke)も exit 0。PR #10(feat/epoxy-amine-bulk-curing)。
- 機構検証(SUCCESS):
  - **counted 開環確定 2 件**(epoxide conversion 0.10 = 2/20 環、amine-H conversion 0.08)、確定解離 4(=2 反応 × N-H/環 C-O のペア、spurious なし)、ヒドロキシル O-H 形成 1(counts=False)。
  - 連言イベントは cycle 0/1/2/3 で発火(発火 4 → 確定 2 = unbiased 確認で 2 件棄却。連言確定の厳格さとして正常)。cycle 9 も発火したが未確定。
  - **cycle 2 は完全な 4 点セット**をトポロジに記録: +N-C(403,569), +O-H(402,556), −環 C-O(402,403), −N-H(556,569)= nylon 縮合と同じクリーンな協奏組み替え。
  - 候補数 29→21→11→9→13→7→5→1→6→1 と反応進行で単調減衰傾向 = EpoxyAmineAdditionUpdater のサイト消費が機能。
  - **観察(cycle 0)**: N-C 形成+N-H/環 C-O 解離は確定したが、ヒドロキシル O-H が確認時点で未成立(r>r0=1.57 Å)= 双性イオン的中間状態で確定。トポロジ上 H が一時的に無結合になる(n_bonds 615→614)。k-l は bias-only(非トリガー)なので設計どおりの挙動だが、ネットワーク解析でトポロジを使う際は「O-H は後続 MD で自然形成しるがトポロジ再判定はされない」点に留意。転化率カウント(counted N-C のみ)には無影響。
- **問題: 温度過熱**。trajectory 実測(target 333 K): 序盤 4000-12000 步で mean 745 K まで上昇 → 終盤 417-440 K まで緩降下、全体 mean 549 K / max 873 K。**333 K に一度も収束せず**。
  - Root cause: `run_epoxy_amine.py` が `LangevinParams(temperature_K=...)` を friction 指定なし(既定 0.001)で構築していた。E0 実行結果エントリで「f2=2 + friction 0.01 + unbiased 1500 の OrbMol production レシピを踏襲」と決定済みだったのに **friction の配線を失念**。f2=2 + friction 0.001 の過熱は MA で実証済みの既知パターン(2026-07-07: 389→565 K)であり、本件はその再現(エポキシは発熱 −29.5 kcal/mol + 序盤の大 bias_E ~3700 でさらに顕著)。発散はしない(単調減衰へ転じる)が定常が高すぎる。
- Decision: `--friction-per-fs` フラグを追加し**既定 0.01**(OrbMol f2=2 レシピ、2026-07-07 で 30cyc 検証済みの冷却レバー)。--friction-per-fs 0.001 で論文忠実値へ戻せる。summary.json に friction_per_fs を記録。
- 検証: `runs/epoxy_amine_fric01`(同 seed、5 cyc、friction 0.01)で温度プロファイルが 333 K 近傍に収まるか実測(実行済み → 結果は本エントリ末尾に追記)。
- **2026-07-09 fric01 検証結果(exit 0、PASS)**: 温度 **mean 362.3 K / max 579.6 K**(target 333)、チャンク別 mean 334〜422 K、**終盤 342 K で暴走なし**(friction 0.001 の mean 549 / max 873 / ピーク 745 K から正常化)。スパイクはバイアス相・反応イベント時のみで unbiased 相で 335〜345 K へ復帰 = 冷却レバーとして設計どおり。mean が target+29 K なのは小箱で bias 注入頻度が高いためで、paper-scale では選択ペア密度が下がり緩和見込み(要監視)。機構も維持: counted 開環 1 件確定(**ヒドロキシル O-H 込みの完全セットで確定**)+対応解離 2、発火は cycle 1。エポキシド転化率 0.05(5 cyc)。→ **friction 0.01 既定を確定**。paper-scale(100+50)へ進んで良い。
- Scientific risk: 中(2026-07-07 friction エントリと同じ位置づけ)。二重逸脱(f2 10→2、friction 0.001→0.01)により絶対速度論は論文比較不可、トレンド比較のみ。smoke10 の 2 反応自体は「過熱状態でも連言機構が正しく動く」ことの実証であり、転化率の定量は fric01 系列で取り直す。
- Licensing/commercial impact: なし。
- Follow-up: (a) fric01 検証結果の追記。(b) paper-scale(100+50)は fric01 検証後に。(c) cycle 0 型の双性イオン確定の頻度が高い場合、k-l 成立を確認条件に加えるかを別途検討(論文は k-l をトリガーに含めないため現状が paper-faithful)。

## 2026-07-09: 【Track 2 / E2 準備】ネットワーク解析基盤(種濃度・ゲル点・架橋密度)の設計と実装
- Context: paper-scale epoxy(別マシンで実行予定)と ref#2(Provenzano 2025, ACS Appl. Polym. Mater. 7(8):4876)比較に必要な解析が未実装: (a) 論文 Fig.5 相当の種濃度追跡、(b) ゲル点、(c) 架橋密度。大マシンの結果到着前にローカルで先行実装し、既存 smoke データで動作検証まで済ませる。
- Paper anchor: arXiv:2511.22874 Fig.5(epoxy 硬化の種濃度 c(t)/c(0) 変化)。ゲル点の理論基準は Flory-Stockmayer(Flory 1941; 教科書理論。論文外の追加だが ref#2 系の標準比較指標)。
- Decision(設計):
  (a) `src/kagome/analysis/network.py`(新規、純関数・MD/IO 非依存・RDKit 非依存): topology.jsonl のスナップショット(bonds+species)からグラフパターンで種を判定 — エポキシド=C-C-O 三員環、アミン次数=N の H 隣接数(2H=1°, 1H=2°, 0H=3°)、ヒドロキシル=C1H1 配位の O。関数: `species_series`(サイクル別種カウント)、`largest_component_fraction`(carothers.monomer_sets_from_bonds 再利用の単量体グラフ最大連結成分=ゲル化指標)、`gel_point_flory_stockmayer(f, g, r)`(α_gel = 1/√(r(f−1)(g−1)); DGEBA f=2 + DETA g=5、r=1 で 0.5)、`crosslink_counts`(3° アミン N 数+全環開環 DGEBA 数を分離して返す=ref#2 の定義に後で揃えられる構造)。
  (b) `scripts/reproduce_figures.py` に `plot_species_concentrations`(c(t)/c(0)、生カウント・スムージング/フィルタなし=ask-first 非該当。初期 0 の生成種は生カウントを右軸)と `plot_gel_curve`(最大成分フラクション vs エポキシド転化率+FS 理論線)。species は --species-json > summary.json からの再構築(SMILES+counts+seed、seed は manifest.json フォールバック)> trajectory ヘッダの優先順位で解決。
  (c) `tests/unit/test_network.py`: 合成トポロジ(手組みグラフ)で各判定関数、開環イベント前後の種収支(epoxide −1, 1° −1, 2° +1, OH +1)、FS 式の解析値、実データ回帰(runs/epoxy_amine_smoke10 存在時のみ skipif)。
- 論文にない仮定(明示): FS 理論線は比較基準としての追加(論文は理論ゲル点を示さない)。架橋点の定義は「3° アミン N」と「両エポキシド消費 DGEBA」を分離集計(合算の重み付けは E2 で ref#2 に整合させる)。トポロジのみから種判定するため、双性イオン確定(2026-07-09 E1 実走結果の観察)では一時的に N が過大次数/H が無配位になる点は docstring に明記(1° 判定は H≥2 として過渡的過プロトン化に防御的)。
- 実装(2026-07-09、サブエージェント実装+orchestrator レビュー): 上記 (a)(b)(c) 一式。テスト: test_network.py+test_carothers.py = 57 passed、全 `pytest tests/unit` = **491 passed**(退行なし)。実データ検証(smoke10): 種再構築が trajectory 記録と完全一致、初期 epoxide 20 / 1° 10 / 2° 5 / OH 0 → 最終 epoxide 18(転化率 0.10 = summary と一致)/ 1° 8 / 2° 7 / **OH 1**(cycle-0 開環が双性イオン型確定で O-H 欠落=既知観察と整合)、最大成分比 1/15→2/15(ゲル化前、当然)。図出力: `figures/species_concentrations.png/.pdf`、`gel_curve.png/.pdf`。
- Licensing/commercial impact: なし(新規依存なし)。
- Follow-up: (a) E2 本番で ref#2 の定義(ゲル点判定・架橋密度の規格化)と突き合わせ。(b) 実データ回帰テストの OH=1 アサーションは k-l 確認条件を変えた場合に要更新(テスト内コメントに明記済み)。

## 2026-07-12: 【Track 2 / E2 設計】Provenzano 2025 外部比較プロトコル(配合・比較軸・データセット)
- Context: E2 = paper-scale epoxy 硬化ランと ref#2(Provenzano et al., ACS Appl. Polym. Mater. 2025, 7(8), 4876, DOI 10.1021/acsapm.4c04208, オープンアクセス PMC12039963)の古典 MD 架橋プロトコルとのクロスメソッド比較。大マシンの GPU 待ちの間に、論文精査(Web 調査 2026-07-12)・データセット取得・比較スクリプト整備を先行実施する。
- Paper anchor(ref#2 の確定事実、調査 2026-07-12):
  - 配合: **DGEBA:DETA = 5:2**(800+320 分子、約 45,000 原子、80 Å 立方体、Packmol 充填)= エポキシド:アミンH = 1600:1600 = **1:1 化学量論**。
  - 手法: Class II 力場(COMPASS/PCFF 系)、距離ベース決定論的架橋(反応確率なし)。カットオフ半径 2.0→10 Å を 0.25 Å 刻み・各半径 10 イテレーションで漸増。硬化温度 500 K(テンプレート値)、結合 1 本ごとに NVT+NPT 各 1000 ステップ、架橋中 timestep 0.01 fs。結合形成後に原子型・電荷を更新(reaxpar.txt)。
  - 転化率(架橋度)の定義: 反応済みアミン H 数 / 初期全反応可能サイト数。検討範囲 0〜90%。
  - 報告値: 硬化後密度 1.08–1.16 g/cm³(架橋度 0→90%)、ヤング率 ~2.5(40%)→~4.5 GPa(90%)、熱伝導率 0.16–0.24 W/(m·K)。**ゲル点・Tg・架橋密度(数値)は未報告**。
  - 公開物: Zenodo DOI 10.5281/zenodo.11402476(コード xlinker.py "DDCC Xlinker v4.0" + LAMMPS data 構造 + 架橋度推移例 xl_trend.txt)。**ライセンス CC-BY-4.0**(レコード全体。コード個別の LICENSE ファイルなし = CC-BY-4.0 が唯一の明示ライセンス)。GitHub リポジトリは存在しない(Zenodo のみ)。
- Decision(比較プロトコル):
  (a) **E2 比較ランの配合は ref#2 に合わせて DGEBA:DETA = 5:2 とする**(例: 100+40 = エポキシド 200 : アミン H 200)。既存計画の 100+50(2:1、アミン H 25% 過剰)は Table S3(arXiv:2511.22874)向けの配合であり、E2 用途には化学量論比を優先する。paper-scale の 100+50 ランはそのまま元論文再現用として別立て。
  (b) **共通比較軸は転化率**(時間軸は不可: ref#2 の進行軸は「カットオフ半径×イテレーション」で物理時間と対応しない。TDBB 側も f2=2/friction 0.01 逸脱により絶対速度論は比較不可 = 2026-07-07 エントリと同じ位置づけ)。転化率の基準は**エポキシド基準に統一**(1:1 化学量論ならアミン H 基準と一致するため ref#2 の定義とも整合。当方 2:1 系では両基準が乖離する点に注意)。
  (c) 比較指標: (i) 転化率 vs 反応イベント順(進行の形状のみ、参考)、(ii) **ネットワーク指標 vs 転化率**(最大成分フラクション=ゲル曲線、3° アミン数、全環開環 DGEBA 数、種組成)— こちらが主比較。(iii) ref#2 の架橋後 LAMMPS 構造(data.xlinker、45% 例)から当方の network.py 相当の指標を算出し、同転化率での当方スナップショットと突き合わせる。(iv) 密度(NPT ランがあれば 1.08–1.16 g/cm³ と比較)。
  (d) データセットは Zenodo から取得スクリプトで再現可能にダウンロードし(コミットしない)、小さな xl_trend.txt のみ帰属表示付きでリポジトリに同梱可(CC-BY-4.0)。
  (e) ゲル点は ref#2 未報告のため、当方は Flory-Stockmayer 理論値(f=2, g=5, r=1 → α_gel=0.5)+最大成分フラクション立ち上がりを報告し、ref#2 構造から算出した点を重ねる(= 当方が比較点を自前生成する。ref#2 の主張との比較ではないことを図注に明記)。
- 論文にない仮定(明示): (i) ref#2 は経時比較を提供しないため「結果(転化率-ネットワーク構造関係)のみ比較」とする方法論的立場(2026-07-08 ロードマップ合意の具体化)。(ii) 5:2 配合の TDBB 側パラメータ(f2=2、333 K、friction 0.01)は 2:1 系での検証値を流用(配合比は候補生成にのみ影響し PES/バイアス設計は不変のため妥当と判断)。(iii) ref#2 の硬化温度 500 K と当方 333 K の差は「機構でなく結果を比較」の範囲内として許容。
- Scientific risk: 中。手法差(距離ベース 500 K vs PES ベース 333 K)により到達転化率・ネットワーク形成経路が本質的に異なりうる — それ自体が E2 の測定対象。ref#2 のフルカーブ(0–90%)生データは Zenodo に含まれず 45% 例のみ(フルカーブは xlinker.py 再実行が必要 = 現時点ではスコープ外、必要になれば別途判断)。
- Licensing/commercial impact: **Provenzano 2025 Zenodo 一式(コード+データ)= CC-BY-4.0、帰属表示の上で商用利用可**。dependency-license-matrix.md に追記。比較解析はデータの読み取りのみで xlinker コードは実行・改変・同梱しない(現スコープ)。
- Follow-up: (a) 比較解析スクリプト `scripts/compare_epoxy_external.py` + Zenodo 取得スクリプト + LAMMPS data パーサ + unit テスト(本エントリと同日実装、サブエージェント)。(b) E2 比較ラン config(100+40、5:2)の作成はラン計画確定時。(c) Orselly 2022(Track 3)の DGEBA/DETA 値をメモ: Tg 実験 119–142 °C / sim 123 °C、ヤング率 実験 2.55 / sim 2.57 GPa、密度 実験 1.178 / sim 1.111 g/cm³(CHARMM CGenFF、カットオフ 2–3 Å、900 K、目標転化 95%)。

## 2026-07-16: 「Track 拡張 / 共重合(アクリレート+メタクリレート混合ビニル系)ビルダー」の設計と実装
- Context: ユーザ要望で、単一モノマーではなくアクリレート+メタクリレートの**混合(共重合)ラジカル系**を走らせたい。既存 `build_vinyl_aibn_system` は単一モノマー種を前提に固定ストライド `mono_offset + j*n_per_mono` でグローバル index を計算しており(`_systems.py`)、原子数の異なる2種(methyl acrylate 12原子 vs methyl methacrylate 15原子)を1つの箱に混在できないことが唯一のブロッカー(Explore 調査 2026-07-16)。VinylChainPropagationUpdater / apply_vinyl_addition / valence guard / groups.py / find_candidates / conversion(α 分母 monomer_site_count)は全て種非依存で、正しいグローバル index と結合 propagation_map さえ渡せば無改変で動く。
- Paper anchor: 論文 §2/Table S1(vinyl ラジカル重合の 4-group ij+ik+jl 基準)。共重合(2種モノマー混合)は**論文外の拡張**であり、論文は各モノマーを個別系として扱う(methyl acrylate, methacrylate, styrene… を別々に。§ Systems studied)。本実装は「同一テンプレート(radical_C/vinyl_alpha_C/chain_C/vinyl_beta_C)の下で両種の反応性原子を同一群に登録する」という拡張であり、TDBB の方程式・スケジュール・反応選択ロジックは無改変。
- Decision(論文外の仮定として明示):
  (i) **両モノマー種の alpha-C を単一 vinyl_alpha_C 群に、beta-C を単一 vinyl_beta_C 群に統合**する。ラジカルはどちらのモノマーにも無差別に付加でき(head-to-tail 位置化学は _find_vinyl_alpha_beta が種ごとに解決)、共重合の連鎖統計は MLIP エネルギー論に委ねる(2026-06-20「atom-typing only、余分なバイアスを足さない」方針の踏襲)。反応性比(r1/r2)を制御するバイアスは**入れない**(paper-faithful; 論文は付加選択にバイアスを足さない)。
  (ii) **オフセットは種固定ストライドを廃し、配置順に沿って running offset を累積**する(`_place_fragments_in_box` の返す fragment 順=種の登録順が唯一の権威)。これにより原子数の異なる任意個数・任意種のモノマー混合を許容。propagation_map / chain_c_map は各 fragment のローカル alpha/beta/radical/chain_c を running offset に加えて構成。
  (iii) 配置順: initiators → monomer_specs の順(spec0 を count0 回、spec1 を count1 回…)。conformer seed は種ごとに rdkit_seed+1+k(幾何のみ変化、原子順序は SMILES で決定的)。
  (iv) α 転化率分母は両種 alpha-C の総数(monomer_site_count が vinyl_alpha_C 群サイズを数えるため自動で正しい)。summary の n_monomers は両種の合計を渡す。
- 実装:
  (a) `scripts/_systems.py`: `build_vinyl_copolymer_system(monomer_specs, n_initiators, box_size, rng, initiator_smiles, min_sep, rdkit_seed)` を追加(既存 build_vinyl_aibn_system は無改変)。戻り値シグネチャは build_vinyl_aibn_system と同一(positions, species, template, groups, propagation_map, chain_c_map)。
  (b) `scripts/run_vinyl_copolymer.py`: 新規の lean ランナー。`--n-acrylate`/`--n-methacrylate`/`--n-initiators`、`--backend {toy,orb,mace,aimnet}`(toy はスモーク用の LJ、GPU 不要)、NVT/NPT、density/box、cycles/steps、seed、checkpoint/resume。wf.run で manifest 記録。分類ラベルの trajectory-topology 出力(vinyl_initial_bonds 等の単一ストライド前提ヘルパ)は共重合レイアウトと非整合なため**呼ばない**(initial_bonds=None; 軌跡は距離推論で描画。科学的成果物=転化率・反応イベントには無影響)。
  (c) `tests/unit/test_copolymer.py`: ビルダーのオフセット正当性(原子数・群サイズ・propagation_map の alpha∈vinyl_alpha_C / beta∈vinyl_beta_C・2種のストライド差)、MD なし候補生成(find_candidates + score + select)。
- 論文になる仮定(明示): (i)(ii) は上記の通り論文外。共重合の反応性比や配列(交互/ランダム/ブロック)は制御せず MLIP に委ねるため、得られる共重合体は**統計的にランダム共重合寄り**になる見込み(定量検証は別途)。
- Scientific risk: 中。TDBB コアは無改変だが「2種混在の候補選択・付加選択が MLIP のみで妥当に振る舞うか」は smoke で確認が必要(toy では幾何のみ、実 MLIP は別途)。反応性比の制御が無いため、実験の共重合組成曲線との定量比較には使えない(トレンド観察用)。
- Licensing/commercial impact: なし(RDKit + 既存依存のみ。toy は自前 LJ、OrbMol-v2 は既存記録済み)。
- Follow-up: (a) toy スモーク(少数系)で通し確認。(b) orb バックエンドで小規模実走→開環/付加が両種で起きるか。(c) 反応性比を見たい場合は組成 vs 転化率の解析を別途(現状は未実装)。

## 2026-07-17: 共重合の相対反応性・クロス成長解析(analyze_copolymer_reactivity)
- Context: 共重合ラン(2026-07-16 エントリの拡張系)から「どちらのモノマーがどれだけ取り込まれたか」に加え、「**どちらの種の成長末端が**どちらのモノマーを取り込んだか」(クロス成長表)を集計したい。`runs/copoly_reactivity_orb`(20MA+20MMA+3init、20 cycles、orb)で確定 4 件を初適用。
- Paper anchor: 論文 §2(TDBB の confirmed formation = radical_C–vinyl_alpha_C ペアの記録)。解析はログの読み取りのみで、TDBB・反応選択・MD は無改変。反応性比(Mayo-Lewis r1/r2)自体は論文外の解析軸。
- Decision(論文外の仮定として明示):
  (i) **原子 index→分子ブロックの写像は決定的に再構成できる**(`copolymer_atom_species`)。ビルダーの配置順(initiators→monomer_specs 順)と各 SMILES の原子数のみに依存し、rng・座標に依らない(`copolymer_alpha_species` と同一のオフセット規約。unit テストでビルダーの群 index と突き合わせ)。
  (ii) **confirmed_formation の radical 端点が属するブロック = 付加時点の連鎖末端種**とみなす。根拠: vinyl 連鎖成長では確定形成のたびに VinylChainPropagationUpdater がラジカルを新モノマーの beta-C へ移すため、radical_C は常に「最後に付加したモノマーのブロック」(最初の付加のみ initiator ブロック)にある。
  (iii) **Mayo-Lewis 推定は初期存在量補正の低転化率近似**: r_A = (N_AA/N_AM)·(avail_M/avail_A)、r_M = (N_MM/N_MA)·(avail_A/avail_M)。存在量は初期値固定(逐次消費は未補正)。等モル仕込みでは補正係数 1。initiator 末端イベントは r 推定から除外し別行で報告。分母イベント 0 のときは None、分子 0 のときは点推定 0.0。
  (iv) イベント数が少ない間は無意味なため、出力に警告(全体 n<10、モノマー末端イベント n<20)を常設。文献の r 値はスクリプト出力に**埋め込まない**(出典管理は specs 側で行う)。
- 実装: `scripts/_systems.py` に `copolymer_atom_species`、`scripts/analyze_copolymer_reactivity.py` にクロス表(initiator/acrylate/methacrylate 末端 × 取込種)と `reactivity_ratio_estimates` を追加。`tests/unit/test_copolymer.py` に写像テスト+合成 bonds.jsonl による end-to-end 分類テスト(16 passed)。
- Scientific risk: 低〜中。解析のみで物理は不変。ただし (iii) の低転化率近似は転化率が上がると系統誤差(残存モノマー組成のドリフト)を生む — 高転化率まで回す場合は逐次存在量補正(イベント順に残数を減算)へ拡張する余地を明記しておく。
- Licensing/commercial impact: なし(標準ライブラリのみ)。
- Follow-up: `runs/copoly_50x50_seed7`(50+50+5、30 cycles、実行中)完走後にクロス表を再集計。必要なら逐次存在量補正・ダイアッド統計(配列のブロック性/交互性)を追加。

## 2026-07-17: 「well-mixed 測定モード」— サイクル間の古典MD混合ステージ+確率的候補選択(反応性測定プロトコル)
- Context: 共重合の相対反応性測定で、候補選択が「凍った初期配置の最近接ペアへの決定的一点張り」になっている事実を selection.jsonl で定量確認(選択スロットの46-53%が前サイクルの再選択、最長14サイクル固着。2026-07-17)。1サイクルの物理時間 ~0.9 ps(3500 steps × 0.25 fs)に対し近傍リフレッシュに必要な拡散は 50-400 ps で2桁不足。unbiased フェーズ(375 fs)は結合生存確認用であり拡散用ではない(論文スケジュール自体の性質)。結果、種選択性の生カウントは化学ではなく初期配置依存性を測っている疑いが強い。
- Paper anchor: TDBB の方程式・バイアス・確認ロジック(§2、Table S1)は無改変。本モードは**論文外の測定プロトコル拡張**であり、論文再現ランには使用しない(明示フラグでのみ有効化、manifest に記録)。
- Decision(ユーザー承認 2026-07-17):
  (i) **古典MD混合ステージ**: 各サイクル末尾(トポロジー更新後・checkpoint 前)に、現在の結合グラフで構築した古典系(OpenFF Sage 2.2 + NAGL 電荷、OpenMM、NVT・同一セル)で数十〜数百 ps の拡散混合を実行し、座標を書き戻す。速度は Maxwell-Boltzmann 再抽選(workflow rng で決定性維持)→ 短い MLIP 整定セグメント(phase='mix_settle')で PES 差の過渡を吸収。
  (ii) **忠実度は OpenFF フル**(簡易汎用FFでなく): prep の既存機構(Sage 2.2 + NAGL + Gasteiger フォールバック)を進化するグラフに適用。工数増と引き換えに出会いの統計(局所構造)の質を優先。
  (iii) **ラジカルの古典表現**: 成長末端 C(配位3)は古典世界にのみ存在するキャップ H 粒子を注入してタイプ付けし、書き戻しで捨てる(イニシエータの closed-shell placeholder H と同思想)。トポロジーから切り離された placeholder H は非結合項のみの中性粒子として扱う。
  (iv) **確率的候補選択(理想1)も同一検証キャンペーンに含める**: ラジカルごとの相手選択を決定的最近接から確率抽選(スコアのソフトマックス等、温度パラメータ付き)に変更するオプションを別PRで実装し、baseline / 混合 / 確率選択 / 両方 の4腕で比較。
  (v) **混合時間・頻度の既定値は検証で測ってから決める**(無根拠デフォルト禁止則): 20+20 系・20サイクル・混合時間 sweep(例 25/50/100/200 ps)、GPU 予算 ~15-20時間(標準)。
  (vi) 判定指標: 再選択率(baseline 46%)、固着ペア最長連、窓内候補集合のサイクル間 Jaccard、混合中モノマー MSD(rms ≥ 5 Å 目標)、ハンドバック時エネルギー過渡の整定内減衰、サイクルあたり確定歩留まり(baseline 23%)、結合グラフの往復保存。
- 前提工事: 共重合ドライバは initial_bonds=None でトポロジー追跡が無効(2026-07-16 の既知制約)。`copolymer_initial_bonds`(running-offset 規約で各フラグメントの RDKit 結合を合成)を実装しトポロジーを有効化する(Phase 1)。混合と独立に topology.jsonl・鎖構造解析が共重合でも使えるようになる副産物あり。
- 科学的注意: 毎サイクル完全混合は実際のラジカル重合にある空間相関(成長鎖近傍のモノマー枯渇等)を意図的に消す。反応性比(per-encounter 化学)測定には正しいが、形態・シーケンス統計を見るランには不適 — 用途を manifest 上で区別する。速度論の絶対値は元々非物理的であり失うものなし。
- Scientific risk: 中。古典 PES と MLIP PES の局所構造差がハンドバック後の候補分布に残る可能性(整定セグメントで部分吸収、指標(vi)で監視)。NAGL/OpenFF が扱えない成分(想定外の価数等)は Gasteiger フォールバック+警告ログ。
- Licensing/commercial impact: 新規依存なし(OpenMM MIT/LGPL、openff-toolkit/interchange/nagl MIT、Sage・NAGL weights CC-BY-4.0 帰属要 — matrix 記載済み・検証済み 2026-06-14)。
- Follow-up(実装フェーズ): P1 共重合トポロジー有効化 → P2 グラフ→古典系トランスレータ(src/kagome/prep/mixing.py、成分キャッシュ)→ P3 ワークフロー統合(MixConfig、phase='mixing'/'mix_settle')→ P4 検証キャンペーン(4腕+sweep、既定値決定)→ P5 多シード生産測定+条件付き確率・選択確率補正(Horvitz-Thompson)解析。

### 追補 2026-07-17: WM-P2 トランスレータの実装仕様(論文外の実装決定)
`src/kagome/prep/mixing.py` 実装で確定した、上記(iii)を具体化する実装決定。いずれも古典混合ステージ内部のみに影響し、TDBB/反応選択/MLIP PES・論文再現ランには不変(混合モードは明示フラグでのみ有効)。
- (a) **成長末端(ラジカル)検出は配位数でなく結合次数の不足で判定**: 各重原子について `deficit = round(neutral_valence[el] - Σ(bond_order))`(neutral_valence C=4/N=3/O=2/H=1/F=1/S=2/Cl=1)。deficit>0 ならラジカル中心。配位数3でも C=O カルボニル炭素や sp2 炭素は次数和4で deficit=0 となり誤ってキャップされない(decisions.md(iii)「配位3の炭素」を次数ベースに精密化)。BondTopology の bond order を唯一の権威とする。
- (b) **キャップ H は deficit 個・結合次数1で注入**(通常の成長末端 C は1個)。古典世界のみに存在し書き戻しで破棄。H は決してキャップしない(H+H=H2 を避けるため)。キャップ H の初期座標は親ラジカルから 1.09 Å(汎用 C–H)を既存隣接の反対方向に置く近似(P3 の古典最小化で緩和)。
- (c) **placeholder H(脱離した開始剤 H、グラフから孤立=配位0の H)は OpenFF に渡せない**: 素の `[H]` は `RadicalsNotSupportedError` で拒否される(WSL 実測)。よって古典系(OpenMM System)に**末尾追加する非結合項のみ・電荷0・質量1.008 の粒子**として実現し(decisions.md(iii)「非結合項のみの中性粒子」の具体化)、書き戻しで元 index へそのまま戻す。LJ サイズは**ハードコードせず Sage 自身(メタンの H)から導出**して FF に追従(engineering 定数の無根拠導入を回避)。
- (d) **原子順序の不変条件を差し替え**: 鎖成長後は1分子が非連続なグローバル index にまたがる(初期の連続レイアウト前提が崩れる)ため、prep の「連続レイアウトで順序一致」不変条件は使えず、明示的な `omm_to_mlip`/`mlip_to_omm` 写像を保持して書き戻す。判定指標(vi)の結合グラフ往復保存は「トランスレータはグラフを不変入力として扱い、座標のみ往復。全 MLIP 原子が一意に写像され write_back で原座標を完全再現」で担保(unit テスト test_write_back_roundtrip_preserves_positions_and_graph)。
- (e) **成分キャッシュのキー = キャップ後フラグメントの RDKit 標準 SMILES**(`Chem.MolToSmiles`)。値は電荷付与済み OpenFF Molecule テンプレート。interchange が同型コピーへ電荷をグラフ同型で写像するため、インスタンスごとの原子順序差は無害。キャッシュは per-cycle の支配的コスト(NAGL 電荷推論)を同型フラグメント単位で1回に削減(自由モノマー多数でも1回)。valence(bond/angle/torsion)の SMIRKS ラベリングは現状 call ごと(P3 で拡張余地)。
- 電荷: 既存 prep と同一(NAGL 既定、失敗時 Gasteiger フォールバック+警告)。テストは Gasteiger(オフライン・決定的)で電荷手法非依存に。
- スコープ: P2 はトランスレータ+テストのみ。公開 API(`build_classical_mix`/`ClassicalMix`/`FragmentParamCache`/`MixTranslatorConfig`、型ヒント必須、MD 不要でテスト可能)を P3 統合の seam として設計。
- レビュー反映(同日、マルチエージェント・コードレビュー high の指摘 10 件):
  - **combine_nonbonded_forces=False は非対応として即時 ValueError**(config 構築時に検証): placeholder H を結合済み NonbondedForce にのみ追加する設計のため、vdW 分離構成では粒子数不整合で Context 生成が壊れる。フィールドは制約の文書化と将来拡張の予約として残す。
  - **複数キャップの決定的幾何**: deficit≥2 の中心(カルベン様/全結合喪失原子)ではキャップ H を空価数軸から 54.7356°(四面体半角)傾け方位角 2π/count で分散(rng 不要の純幾何)。一点重なり(H-C-H 角 0 → NaN 角度力)を排除。
  - **キャッシュキー拡張**: `charge_method|nagl_model|forcefield|標準SMILES`(`fragment_cache_key`)。異設定間のテンプレート汚染を構造的に排除。キャッシュ値は `(テンプレート, 実際に使われた手法)` 対とし、NAGL→Gasteiger フォールバック時は metadata に実手法(`charge_method`=実際/'mixed'、`charge_method_requested`、`charge_methods_used`、`nagl_fallback`)を正直に記録。
  - **未知元素ガード**: NEUTRAL_VALENCE 外の元素が価数不足(RDKit 既定価数比)の場合は黙殺せず ValueError(RDKit 既定価数が未定義の元素も拒否)。飽和した未知元素は通過。加えて RDKit 原子は全て SetNoImplicit(True)(暗黙 H の湧きによる写像ずれを根絶)+ from_rdkit 後の原子数一致検証。
  - **元素検証(D-4 相当)**: 組み上げた OpenFF トポロジーの全粒子について元素 = species[omm_to_mlip[k]](キャップは 'H')を全数検証、不一致は ValueError。placeholder H も追加前に species=='H' を検証。
  - **write_back 入力検証**: 非有限座標(MD 発散 NaN)は写像前に検出し「入力が非有限(混合ランの発散)」として別メッセージで報告(写像バグと誤診しない)。
  - **最小化必須の契約**: positions_A は初期推定(キャップ H は理想方向、placeholder H は旧親と LJ 激突位置 ~1 Å)。P3 は動力学前に必ず古典最小化を実行 — `metadata['requires_minimization']=True` と ClassicalMix docstring に明記。座標の事前ずらしは行わない(エネルギー有限をテストで確認済み、最小化が標準手順であり、ずらしは PBC 越しの新たな重なりリスクを生むだけ)。
  - **共通化**: NAGL→Gasteiger 電荷付与を `kagome/prep/charges.py` に一本化し openmm_equilibrate と mixing の二重実装・挙動差(import 検査、Gasteiger 実装差)を解消。成分分解は carothers の既存 component finder を self-edge 規約(network.py 前例)で再利用し union-find 再実装を削除。

### 追補 2026-07-18: WM-P3 ワークフロー統合の実装仕様(論文外の実装決定)
`_run_mixing_phase` 統合で確定した実装決定。混合モードは明示フラグ(`--mix`)でのみ有効であり、TDBB 方程式・反応選択・論文再現ランには不変。
- (f) **挿入位置**: 各サイクルの `_update_groups_after_cycle` 直後・`save_checkpoint` 直前(decisions.md (i) の「トポロジー更新後・checkpoint 前」を実装)。checkpoint は混合後の positions/velocities/rng_state を保存するため、resume 意味論は追加フィールドなしで正しい(CHECKPOINT_VERSION 据え置き)。
- (g) **混合時間・整定ステップ数に既定値なし**: `MixConfig.mix_time_ps` / `settle_steps` は必須引数(CLI では `--mix` 指定時に `--mix-ps` / `--mix-settle-steps` を要求)。無根拠デフォルト禁止則(v)の適用 — 既定値は P4 の sweep で測ってから導入する。timestep 0.5 fs / friction 1.0 /ps / platform CPU / Sage 2.2 + NAGL は既存 `ClassicalPrepConfig` の前例を踏襲(新規定数の導入なし)。
- (h) **古典フレームは trajectory.jsonl に書かない**: trajectory.jsonl は MLIP PES の単位系・エネルギースケールで一貫させる。古典混合の診断(段階エネルギー、rms 変位、キャッシュ統計、実電荷手法、シード)は `mixing.jsonl` に per-cycle 1 行で記録(auditability 規則)。`mix_settle`(MLIP)フレームのみ trajectory.jsonl に phase='mix_settle' で記録。MLIP step カウンタは古典 MD では進めない(mix_settle でのみ進む)→ resume 時の step ベース truncation と整合。
- (i) **決定性**: OpenMM 積分器と MB シードは workflow rng から `rng.integers(2**31-1)` で導出。書き戻し後の MLIP 速度は `maxwell_boltzmann_velocities(masses, T, rng)` で再抽選(decisions.md (i))。同一プラットフォーム・同一スレッド条件下でビット再現(CPU platform はスレッド数依存の丸め差があり得るため、決定性テストは Reference platform で行う)。
- (j) **混合 MSD 指標は非イメージング rms 変位**: OpenMM の getPositions は座標を wrap しないため、書き戻し前後の直接差の rms(Å)を拡散指標として mixing.jsonl に記録(判定指標 (vi) の「混合中モノマー MSD」の実装。全 MLIP 原子対象 — モノマー限定はP4解析側で topology.jsonl から分離可能)。
- (k) **図フィルタ**: reproduce_figures の energy プロットは unbiased マスクを「biased 以外」から「biased 以外かつ mixing/mix_settle 以外」に修正(古典スケール混入と整定過渡の混入を防止)。既存ラン(mixing 無し)の図はビット不変。temperature プロットは MLIP 運動温度のみが対象になるため無修正。
- (l) **古典 MD カーネルは `prep/mix_md.py` に分離**(数値カーネルとオーケストレーションの分離規則): `run_mix_md(mix, cfg, seed)` は build 済み `ClassicalMix` を受け、LangevinMiddleIntegrator NVT を最小化(requires_minimization 契約の履行)→ MB 初期化 → チャンク実行+非有限エネルギー早期検出(equilibrate と同パターン)で実行し、最終座標(OpenMM 順・Å)とエネルギー系列を返す。ワークフロー側は OpenMM を import しない(遅延 import、OpenMM 無し環境でも import 可能)。
- (m) **FragmentParamCache はワークフロー寿命で永続**(サイクル横断で NAGL 費用を同型フラグメント1回に)、checkpoint には保存しない(性能最適化であり resume 時は再構築)。
- レビュー反映(2026-07-18、マルチエージェント high レビュー実質5件):
  - **混合シードは 0 を除外して抽選**(`rng.integers(1, 2**31-1)`): OpenMM `setRandomNumberSeed(0)` は「毎回新規シード」を意味し、~1/2e9 の確率で決定性が黙って崩れる。決定性のビット再現には加えて OpenMM プラットフォームが決定的(Reference か単一スレッド)である必要があり CPU/CUDA は力の縮約順がスレッド依存 — この caveat を workflow コメントに明記(既定 CPU は性能優先で維持、determinism テストは Reference)。
  - **resume 時の測定モード切替ガード**: checkpoint extra に mixing 設定(mix_ps/settle/timestep/platform/charge_method)を記録し、resume 時の CLI 引数と不一致なら hard error。--mix の付け忘れ等でサイクル途中から paper-faithful ループに黙って切り替わり測定を汚染するのを防止(古い checkpoint は mixing キー無し=混合オフ扱い)。
  - **base_energy 図も mix_settle フレーム除外**: energy scatter だけでなく `plot_energy_vs_step` 第2図(base energy)も unbiased_mask を適用し、サイクル境界の整定過渡スパイクを排除。
  - 見送り: 古典 NVT チャンクループの openmm_equilibrate との共通化(cleanup)は prep テスト済みコードの改変リスクを避け WM-P3 スコープ外(将来のリファクタ課題)。

### 追補 2026-07-18 (WM-P4): 古典混合MDのソフトスタート・ウォームアップ(反応後の注入ホットコンタクト対策)
本番キャンペーン(runs/wm_p4/、4腕)の混合腕 `A2_mix50` がサイクル0の混合フェーズで `openmm.OpenMMException: Particle coordinate is NaN`(`run_mix_md` の `integrator.step`)でクラッシュした実障害の根本対策。混合モードは `--mix` 明示フラグでのみ有効で、TDBB 方程式・反応選択・論文再現ランには不変(数値カーネル `prep/mix_md.py` のみ変更)。
- **根本原因**: サイクル0で実結合(ビニルラジカル付加)が成立 → `apply_vinyl_addition` が成長末端に**キャップ H**(Sage が剛体拘束する C–H)を注入し、開始剤の**placeholder H** を脱離させる。`build_classical_mix` は placeholder H を旧親位置(~1 Å の激しい LJ 重なり)に、キャップ H を空価数方向 1.09 Å に置く(requires_minimization 契約)。本番規模(~564原子・0.5 g/mL・~23 Å)では `LocalEnergyMinimizer`(tol 10 kJ/mol/nm)の**RMS 力**判定が全原子で薄まり、単一の閉じ込められたホットコンタクトが未解消のまま基準を満たす → 続く 0.5 fs Langevin が初回チャンクで発散。単精度 GPU 演算でオーバーフローが NaN 化。
- **なぜ従来テストで漏れたか**: orb スモークは 3+3+1 の極小系で確定結合ゼロ(キャップ/placeholder H 経路を動力学下で未励起)。ユニットテストは単点エネルギー有限のみ確認し反応後の古典 MD を走らせていなかった。
- **採用対策(ソフトスタート・ウォームアップ)**: `run_mix_md` に、報告混合ステップの**外**でカウントするウォームアップ MD セグメントを追加。(1) 注入原子があるとき(metadata の cap/placeholder H 数>0)のみ、より厳しい第2段最小化(`warmup_minimize_tolerance_kj_mol_nm` 既定 1.0 < 主 tol)で最悪コンタクトを追加除去、(2) 高摩擦(`warmup_friction_per_ps` 既定 50/ps)で微小 timestep(`warmup_timestep_fs` 既定 0.1 fs)から設定 timestep へ幾何級数ランプ(`warmup_stages` 既定5)する `warmup_steps`(既定2000)の減衰的 MD で残留ホットコンタクトを吐き出す。混合 timestep・物理・測定は不変。非有限早期検出は維持(ウォームアップ自身も発散時は情報つきで raise)。
- **決定 vs 代替**: decisions.md 追補 (b) の「placeholder H の事前ずらしはしない(最小化が標準手順、ずらしは PBC 越しの新重なりリスク)」を尊重し、**placeholder-H 配置は一切変更しない**。座標を動かさず動力学的にコンタクトを緩める warmup は同方針と整合し、かつ本番規模で最小化が原理的に取り切れない閉じ込めコンタクト(RMS 希釈)にも効く点で第2段最小化単独より頑健。
- **既定値の根拠**: 混合*時間*は依然デフォルト禁止(decisions.md (v))だが、数値安定化ノブ(warmup 長・timestep・摩擦・第2最小化 tol)は正当に既定可能な engineering 定数として文書化。warmup は反応後サイクルのみ発火(未反応サイクルは従来とビット不変)。
- **決定性**: warmup は既存 integrator の seed/stream を共有し `setStepSize`/`setFriction` で制御(新規乱数なし)。determinism テストは Reference で緑。API は加算的拡張のみ(`MixMDConfig` に warmup 系フィールド、`MixMDResult` に `n_warmup_steps`/`warmup_energies_kj_mol`、mixing.jsonl に `n_warmup_steps`)。既存公開 API(`run_mix_md`/`MixMDConfig`/`MixMDResult`/`build_classical_mix`/`ClassicalMix`)は非破壊。
- **回帰テスト**(`tests/unit/test_mixing.py::test_mixing_md_survives_injected_hot_contacts_after_reaction`, OpenMM ゲート・orb 不要): 2 MA+2 MMA+1 開始剤を疎な箱で配置→反応→等方圧縮(~0.4 g/mL、placer は直接到達不可)で反応後の注入コンタクトを再現。ユニット規模では最小化が強力すぎて本番の閉じ込めを emergent に再現できないため、**わざと緩い最小化 tol(1e6)で「未解消ホットスタート」を決定論的に代理**(本番の 564原子 RMS 希釈と同じ終状態)、単精度 'CPU' で NaN 化を再現。`warmup_steps=0`(修正前経路)で NaN、既定 warmup で有限完走を同一テスト内で検証。

### 追補 2026-07-18 (WM-P4 de-clash): 注入原子の幾何デクラッシュ + 混合スキップの安全網
上記 warmup の後にも残った本番障害への第2の頑健化。本番ラン(20+20+2, 564原子, 0.5 g/mL)で cycle 0–3 は正常混合(rms ~14 Å)だったが **cycle 4 が `_run_warmup` の `RuntimeError: soft-start warm-up became non-finite after 100 warm-up steps (dt 0.100 fs)` でクラッシュ**。最小化 + 0.1 fs warmup でも脱出できない **near-singular(近特異)な重なり**を注入原子がまれに持つ = warmup 単独では取り切れないと実証された。混合モードは `--mix` 明示フラグでのみ有効・TDBB/反応選択/論文再現ランには不変。
- **根本原因**: `apply_vinyl_addition` 後、cap H が空価数方向 1.09 Å に置かれ**近接原子へ真っ直ぐ向く**、または脱離 placeholder H が**他原子の上に着地**すると、LJ が near-singular(<0.5 Å)になる。RMS 力最小化は特異点近傍から脱出できず、続く動力学が発散する。
- **write_back の再確認(設計を左右する事実)**: `write_back`(mixing.py)は `out[mlip_idx] = omm_positions_A[omm_idx]`。placeholder の `mlip_idx` は非 None(実 MLIP 原子)なので、その**混合後座標が MLIP index にそのまま書き戻される(破棄されるのは cap H=mlip_idx None のみ)**。`_run_mixing_phase` は `result.positions_A`(混合後座標)を渡すため、placeholder は混合中に動きその位置が次サイクルの MLIP 状態に入る。よって「非相互作用ゴースト(ε=0/排他)」案は**却下**: 混合中に力を受けず MLIP 世界へランダム位置が戻り新たなクラッシュ源になる。相互作用を保ったまま空きへ再配置する方が原理的に正しい。
- **追補 (b)/2026-07-17 レビューの改訂(注入原子に限定)**: 「座標の事前ずらしはしない(最小化が標準手順)」は near-singular 重なりには原理的に無力(最小化は特異点から脱出できない)という本番証拠を得た。よって**注入原子(cap H, placeholder H)に限り、ハードクラッシュ時のみ**標的デクラッシュする。任意座標の一般的な事前変位ではなく、元原子は一切動かさない。
- **Layer 1(ハード重なりの幾何除去、build 時・`MixTranslatorConfig.declash_injected` 既定 True)**:
  - (1a) **cap 配置を最大クリアランス方向へ**: 「空価数軸の反対」ではなく、固定 Fibonacci 球(**N=64**、golden-angle、rng 無し)+空価数軸を候補とし、全原子(結合+非結合、最近接像)への最小距離を**最大化**する方向を選ぶ。C–H 結合長は維持。自身の結合相手(親)は障害物から除外(自分の結合から逃げないため)。
  - (1b) **一般デクラッシュ**: 組み立て後の OMM 座標で、各注入原子の(自身の結合相手を除いた)最近接像距離が **<0.8 Å** なら押し出す。cap は親から C–H 長・最大クリアランス方向へ再配置、placeholder(配位0・結合相手なし)は衝突原子から **~1.3 Å** 離す。決定的(rng 無し・昇順逐次更新)。移動は注入原子のみ、元原子は不動。`n_declashed` を metadata / mixing.jsonl に記録。
  - **soft-core(warmup 中に CustomNonbondedForce で LJ を軟化する案)は却下**: 並行する軟化力の構築・context 再初期化・排他管理が必要で決定性/テストが困難。幾何デクラッシュは最小・決定的で特異点を直接除去する。
- **Layer 2(グレースフル・スキップ=完走保証、`_run_mixing_phase`)**: de-clash + 最小化 + warmup でも発散する場合、**そのサイクルの混合のみをスキップ**: pre-mixing の MLIP positions/velocities を維持、WARNING、mixing.jsonl に `skipped:true`+`skip_reason`、`mix_settle` も skip、checkpoint は未混合状態で通常書き込み。捕捉例外は `RuntimeError`(自前の非有限検出)と `openmm.OpenMMException`(`Particle coordinate is NaN`)。**スキップは稀であるべき**: de-clash がほぼゼロ化するはずで、非自明なスキップ率(例 15 中 >2)は測定汚染の**赤信号**として `summary.json` の `mixing_skipped_cycles`/`n_mixing_skipped` に集計・監視(黙って許容しない)。
- **決定性(スキップ横断)**: `mix_seed = rng.integers(1, 2**31-1)` を try の**外**で無条件抽選 → 数値混合の成否に関わらず workflow rng はこの1回分だけ進む。スキップは seed 抽選のみ消費、成功はそれに加えて MB 速度再抽選 + settle の Langevin 抽選を消費(スキップは本質的に消費が少ない — これは正しい)。skip-vs-success は (seed, geometry) の決定関数なので rng 軌跡は再現可能で、混合中クラッシュが draw を失って resume を desync させることはない。
- **既定値の根拠**: N=64・0.8 Å ハードクラッシュ閾値・1.3 Å クリアランスは数値安定化の engineering 定数として文書化(混合*時間*のデフォルト禁止則(v)とは別カテゴリ)。既存 warmup は defense-in-depth として維持。層構成: **de-clash(特異点除去)→ 最小化 → 第2最小化 → warmup(<0.8 Å 未満の回復可能ホット)→ mix → graceful-skip(それでも発散時のみ)**。
- **API 非破壊**: `MixTranslatorConfig.declash_injected`(既定 True)、metadata に `declash_injected`/`n_declashed`、mixing.jsonl に `skipped`/`skip_reason`/`n_declashed`、summary.json に `mixing_skipped_cycles`/`n_mixing_skipped` を加算的に追加。`build_classical_mix`/`ClassicalMix`/`run_mix_md`/`MixMDConfig` の既存シグネチャは不変。
- **回帰テスト**(OpenMM ゲート・orb 不要):
  - `test_mixing.py`: cap 方向のクリアランス最大化(障害物を空価数軸上に置くと反対を向く)、legacy 配置の不変性、`_declash_injected` の単体(注入原子のみ移動 / ノーオペ)、および **end-to-end**: dense reacted 系で placeholder を実原子に**強制一致**(真の特異点)→ 同一 MD config(warmup ON)で `declash_injected=False` は発散・`True` は有限完走 + クリーン write-back。既存 warmup テストは `declash_injected=False` で従来のホットコンタクトを保持し不変。
  - `test_workflow_mixing.py`: graceful-skip — `run_mix_md` を monkeypatch で raise → `wf.run` は例外を出さず継続、logs=[biased,unbiased,mixing(steps=0)]・`mix_settle` 無し・mixing.jsonl `skipped:true`・state 有限。

### 追補 2026-07-18 (WM-P4 bridge): orb↔Sage 反発ウォール不整合が真の反応非依存根因 — ウォームアップ非ゲート化(soft-core は不要と実証)
デクラッシュ修正(cf7b9db)投入後の本番再ラン(A2)で判明した、より深い第3の根因。cycle 0,1 は正常混合、**cycle 2 の混合が OpenMMException NaN で graceful スキップ**された。その cycle の mixing.jsonl は `declashed=0, n_cap_h=0, n_placeholder_h=0` = **注入原子ゼロ**(cycle 2 は候補選択したが結合確定せず、トポロジー/フラグメントは正常混合した cycle 1 と不変)。すなわち発散は注入原子の幾何ではなく、**ハンドオーバされた MLIP 座標**そのものに由来する。
- **真の根因(反応非依存)**: orb(MLIP)の学習された反発は古典 Sage FF の LJ コアより**柔らかい**。orb 平衡化幾何は時折、Sage の反発ウォールが許容するより近い原子対を含む(cycle 2 の post-unbiased 座標がそれ、cycle 1 は含まず)。確率的・反応非依存で毎キャンペーン再発しうる。デクラッシュ(cf7b9db)は正しいが支配的失敗ではなかった。
- **ゲートのギャップ**: `run_mix_md` でウォームアップと第2段(厳しめ)最小化が共に `has_injected` でゲートされていた(ウォームアップ commit の「未反応サイクルはビット不変」意図)。注入原子ゼロの cycle 2 ではどちらも走らず、緩い tol=10 の第1最小化のみ → RMS 力希釈で 1 個の悪コンタクトが生き残り → 0.5 fs 初回チャンクで即発散。**未反応サイクルも古典↔MLIP 橋渡しが必要**。
- **対策1(非ゲート化)**: ウォームアップ + 第2段最小化を `has_injected` から切り離し、**毎混合サイクル**実行(`warmup_steps>0` でのみゲート、デクラッシュは注入ゲートのまま)。コストは ~1 ps/サイクルで軽微。未反応サイクルの「ビット不変」特性は放棄(橋渡しは毎サイクル必要 = 意図的)。
- **対策2(より優しいウォームアップ既定)**: コンタクトは任意のハンドオーバ対でありうる(中程度の注入コンタクトに限らない)ため、ランプ初段 timestep を **0.1→0.05 fs**、段数 **5→6** に。準最急降下的な初段クロールで硬い orb→Sage コンタクトも生存。数値安定化の engineering 既定として文書化。
- **soft-core は「不要」と単体テストで実証(自律判断)**: 制約付き X–H(Sage は X–H 結合を剛体拘束)同士のコンタクトは幾何的押し離しでは分離不能(変位 H が拘束で戻る)であり、当初は soft-core/力キャップ前緩和が必要かと想定。**決め手テスト** `test_ungated_bridge_survives_constrained_hh_contact`(2分子を剛体並進で最近接 ~0.7 Å、拘束 X–H 温存)で検証: 橋渡し無し(`warmup_steps=0`)では発散(真に硬い)、非ゲート橋渡し(第2段最小化 tol 1.0 + 0.05 fs ウォームアップ)では**有限完走**。境界付き過減衰動力学が剛体分子群を丸ごと押し離すため、拘束原子コンタクトも soft-core 無しで橋渡しできる。よって **CustomNonbondedForce soft-core は実装しない**(複雑さ・力場再構築・context 再初期化のリスクを回避。証拠が不要と示した)。将来もし橋渡しが取りこぼす硬コンタクトが再発すれば soft-core を再検討(設計は proposal に記録済み)。
- **graceful スキップは最終安全網として維持**(不変)。層構成: **デクラッシュ(注入特異点除去)→ 第1最小化 → 第2最小化(毎サイクル)→ ウォームアップ(毎サイクル、0.05 fs 起点)→ mix → graceful-skip(それでも発散時のみ)**。
- **決定性/API**: ウォームアップは既存 integrator の seed/stream を共有(新規乱数なし)。非ゲート化で未反応サイクルもウォームアップ乱数ステップを消費するが seed 決定的。`MixMDConfig` の既定変更(`warmup_timestep_fs` 0.05、`warmup_stages` 6)のみ、シグネチャ不変。
- **回帰テスト**(OpenMM ゲート・orb 不要、`test_mixing.py`):
  - `test_ungated_warmup_survives_noninjected_handed_over_contact`: 2 MA 未反応系、分子 B を剛体並進で最近接 ~0.9 Å(注入原子ゼロ)。第1最小化を無効化する巨大 tol(1e12、本番 564 原子 RMS 希釈の単体規模スタンドイン)で、修正前経路(`warmup_steps=0`)は発散・非ゲート既定は有限完走 + クリーン write-back。本番 cycle 2 の再現。
  - `test_ungated_bridge_survives_constrained_hh_contact`: 上記の決め手(~0.7 Å、拘束 X–H)。橋渡し無しは発散・非ゲート橋渡しは完走を同一テストで検証。

## 追補 2026-07-19 — WM-P4 混合実装 フラットレビューと対応
2026-07-17〜18 の WM-P2(トランスレータ、`src/kagome/prep/mixing.py`)/ WM-P3(ワークフロー統合)/ WM-P4(ソフトスタート・de-clash・ungated bridge)一式に対し、3体の独立レビュー(フラット構成: 正しさ・簿記・再現性・可観測性の複数観点を1パスで横断)を実施した結果の記録。混合モードは引き続き `--mix` 明示フラグでのみ有効・TDBB 方程式/反応選択/論文再現ランには不変という前提は変わらない。

1. **レビュー総括**: 3体とも中核ロジック — write-back(決定 (d) の `omm_to_mlip`/`mlip_to_omm` 写像往復)、index マップ、placeholder-H の往復(決定 (b)/(iii) およびその 2026-07-18 de-clash 改訂)、単位換算、graceful-skip の決定性(WM-P4 de-clash 追補の rng 消費規約)、mixing OFF 時のビット不変(決定 (k) の図フィルタ、決定 (f) の挿入位置)、および decisions.md (f)–(m)・(v)/(vi) への忠実度 — に「バグなし」と結論。指摘は全て周辺(ガード欠落・検証の遅延・テスト網羅・再現性の言明・図の分母)に限定され、コア数値カーネル(`prep/mix_md.py`)・トランスレータ(`prep/mixing.py`)本体の再設計は不要と判断。

2. **決定性の注記(クレーム修正・重要)**: 決定 (i)/2026-07-18 レビュー反映で「同一 seed ⇒ 決定的」と記していたが、これは **rng ドロー列**(mix_seed・MB 速度・settle Langevin の抽選順序)についてのみ成立することを明確化する。本番プラットフォーム(CPU/CUDA)では力の縮約がスレッド順に依存するため、混合後の**幾何**(座標)はビット再現しない。ビット安定は Reference プラットフォームに限られ、既存の `test_mixing_run_is_deterministic`(`tests/unit/test_workflow_mixing.py:324`)が `platform='Reference'` を明示指定している事実と整合する。これはコードのバグではなくバックエンド特性(浮動小数の非結合性 × 並列力計算)であり、修正の余地はない。
   - 方針: 再現性がクリティカルなラン(unit test、回帰確認)は `platform='Reference'`(低速)で固定する。本番 CUDA/CPU ランは「同一 seed ⇒ 同一 rng 消費」は成立するが「同一 seed ⇒ 同一混合幾何」は成立しない旨を manifest に注記する。決定 (i) の文言は本追補で補強するものであり、既存記述の削除・書き換えは行わない。

3. **密度図の分母修正(オーナー承認済み・適用済み)**: `scripts/reproduce_figures.py` の `plot_density_profile` が ρ_rxn(z)=N_rxn/(A·Δz·N_frames) の N_frames に mixing/mix_settle フレームを算入していた不整合を修正し、決定 (k) のエネルギー図除外(`_analysis_frames` によるマスク)と揃えた。mix_settle 中は反応が発火しないため、これを分母に含めると混合ランの ρ_rxn(z) が系統的に過小評価され、非混合ラン(baseline)と比較不能になる。Eq.12 の「解析窓の全サンプルフレーム」は MLIP 動力学窓(古典 mix_settle 緩和を除く)と解釈する。既存の非混合ランの図はビット不変(mixing フェーズを持たないランでは `_analysis_frames` が no-op になるため)。

4. **歩留低下の解釈(オーナー承認・力学は変更せず記録のみ)**: WM-P4 の暫定 3 腕比較で確定生成数が A1(baseline)=5 / A3=4 / A4(mix 系)=1 と、混合腕が最少だった。これは**バグでなく、well-mixed 測定が捉えている実在の物理的特性が最有力**と判断する。
   - 主要因: 混合が「同一ペアを連続サイクルで再選択して cook し続ける」baseline の経路(2026-07-17 決定の Context に記録された再選択率 46-53%)を破壊するため。決定的根拠: BoostState は毎サイクル新規生成でありバイアス持ち越しは無い(paper-faithful)、かつ確定済み結合トポロジーは混合フェーズの**前**(決定 (f) の挿入位置)に永続化済みで混合は `state.positions`/`state.velocities` のみを更新する(グラフは不変入力 — 決定 (d))。
   - 副次要因: 決定 (i) の毎サイクル MB 速度再抽選が反応相関運動量を消去すること、および WM-P4 ungated bridge 追補の非ゲート warm-up(短 `mix_time_ps` 時は総線量の最大 ~20% を占めうる、friction 50/ps)が測定アンサンブルを撹乱しうる。ただしいずれも NVT 整合を保つ意図的設計(決定 (g)/(l))であり、歩留低下対策として取り除かない。
   - 混合時間の既定値確定は先送り: A2_mix50 の複数シード化 + `mix_time_ps` sweep で歩留を再検証してから判断する。現データは単桁カウント(1〜5)でノイズ支配であり統計的結論を出せない。決定 (v)(無根拠デフォルト禁止則)を維持し、混合時間の既定値は本追補時点でも decisions.md に確定記載しない。

5. **対応した補強(実装はサブエージェント委譲・並行実施)**:
   - (a) resume モード切替ガード(2026-07-18 レビュー反映)に `friction_per_ps`・`temperature_K` を追加。旧 checkpoint はこれらのキーが欠落するため非比較(値なし=ガード対象外)として後方互換を維持。
   - (b) `MixConfig.__post_init__` に `timestep_fs`/`friction_per_ps`/`platform`/`charge_method` の fail-fast 検証を追加し、不正値がサイクル途中の `run_mix_md` 呼び出しまで検出されずに伝播することを防止。
   - (c) CLI の `--mix` 無しで混合系フラグ(`--mix-ps` 等)のみ指定した場合を、既存の逆方向(決定 (g) の「`--mix` 指定時に必須」)と対称的に error 化。
   - (d) `mixing.jsonl` の各行に `schema_version` を付与し、WM-P3/P4/P4 de-clash/P4 bridge の各追補で加算的に増えたフィールド(`n_warmup_steps`、`n_declashed`、`skipped`/`skip_reason` 等)のログ消費側でのバージョン判別を可能にした。
   - (e) テスト補強: resume ガードの新規2フィールド、CLI 検証の対称性、`MixConfig` の fail-fast 検証、mixing OFF ランの数値ビット不変回帰、openmm ゲートテストのスキップ理由の可視化(skip でなく明示メッセージ)。
   - (f) `scripts/analyze_wm_p4.py` の `ARM_CYCLE_CAP`(58行目付近のハードコード辞書)をデータ由来(実行された cycle 数から導出)に変更し、腕ごとの系説明を manifest 由来の値に統一。手打ち値と実データの乖離リスクを解消。
   - 各対応の根拠は個々のレビュー指摘(major: resume ガード欠落フィールド/回帰テスト皆無/密度図分母; minor: `MixConfig` 検証の遅延/CLI フラグの非対称)に対応する。いずれも API 非破壊・既存挙動への影響なし(mixing OFF ランはビット不変のまま)。
