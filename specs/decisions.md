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
