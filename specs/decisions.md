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
- Context: Eq. 6-7 defines a general framework with groups I, J, K, L and pair set P. The current implementation in scripts/_systems.py uses a 2-group template (C_donor, C_acceptor) with 1 pair for vinyl polymerization.
- Paper anchor: Eq. 6-7, Table of reaction systems (Section 3 and paper examples).
- Decision: 2-group template is correct for vinyl/radical polymerization. Confirmed from arXiv HTML: vinyl radical polymerization uses Gi (radical carbon) and Gj (alkene carbon), pair {(i,j)} only. The 4-group template (Gi, Gj, Gk, Gl with pairs (i,j),(i,k),(j,l)) is specifically needed for epoxy curing on CuO surface, NOT for vinyl polymerization.
- Alternatives considered: Implementing the full 4-group template for vinyl — rejected as the chemistry requires only 1 bond-forming pair.
- Scientific risk: None for vinyl/radical systems. If nylon-6,6 or epoxy systems are added, the template builder must be extended. The selection machinery (src/reactive/selection.py) already supports N groups and arbitrary pair sets; only the test system builder in scripts/_systems.py is system-specific.
- Licensing/commercial impact: None.
- Follow-up: When adding nylon-6,6 or epoxy systems, add a corresponding build_*_template() function in scripts/_systems.py. The 4-group epoxy template requires Gi (epoxy O), Gj (1-deg amine N), Gk (2-deg amine N), Gl (surface OH) with P={(i,j),(i,k),(j,l)}.

## 2026-06-13: Equation numbering discrepancy in analysis modules (src/analysis/)
- Context: During T6.1 (arXiv HTML verification), confirmed that the equation numbering used in src/analysis/conversion.py ("Eq. 11-12") and src/analysis/density.py ("Eq. 13") may not match the actual paper. From the HTML: Eq. 11 = alpha(t) = 1 - exp(-kp_eff*t) (fitting formula), Eq. 12 = depth-resolved density. The raw conversion fraction alpha = N_reacted/N_total appears to be either unnumbered or given a different number.
- Paper anchor: Eq. 11-12 (per HTML numbering; PDF needed for confirmation).
- Decision: Do not change docstrings or code until the PDF confirms equation numbers. The functional implementations are correct regardless of numbering. Flag the discrepancy in claims.yaml and paper/notes.md.
- Alternatives considered: Renaming based on HTML — rejected because HTML rendering may be incomplete for math-heavy sections.
- Scientific risk: None for computation. Risk is in documentation confusion only.
- Licensing/commercial impact: None.
- Follow-up: Verify equation numbers from PDF and update docstrings in src/analysis/conversion.py and src/analysis/density.py accordingly.

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
