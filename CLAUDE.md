# CLAUDE.md

## Mission
This repository exists to reproduce the method in:
- "Ready-to-Use Polymerization Simulations Combining Universal Machine Learning Interatomic Potential with Time-Dependent Bond Boosting for Polymer and Interface Design"

The goal is to build a reproducible, commercial-safe implementation of the paper's workflow for polymerization and curing simulations.

## Operating mode
Use **paper-faithful where specified, backend-agnostic where licensing is unclear**.

That means:
- Reproduce the TDBB workflow, equations, scheduling, and reaction-selection logic as described in the paper.
- Keep the machine-learning interatomic potential behind a clean backend interface.
- Do not make any proprietary or license-unclear backend the default.
- Prefer commercial-use-permitted software and clearly documented licenses.

## Non-negotiables
- Every implementation change must reference at least one paper artifact: claim, equation, figure, table, or method paragraph.
- Any assumption not explicitly supported by the paper must be written to `specs/decisions.md` before implementation.
- Figures must be reproducible from scripts; no manual editing.
- All experiments must record seed, config path, git SHA, backend name, and output directory.
- Keep TDBB implementation independent from any specific model provider.
- Never paste large copyrighted passages from the paper into repo files; store only short notes and structured summaries.

## Commercial-use guardrails
- New dependencies require an explicit license check before adoption.
- New model backends require a documented commercialization status in `specs/dependency-license-matrix.md`.
- If software is permissive but model weights are restricted, treat the backend as restricted.
- If license status is unclear, mark it `blocked_pending_review`.
- Do not assume that use in a paper implies permission for commercial use.
- Do not make PFP/Matlantis-style backends the default unless usage rights are confirmed by the user.

## Repository map
- `paper/`: structured notes from the paper
- `specs/`: requirements, decisions, acceptance criteria, experiment matrix, dependency licenses
- `configs/`: all reproducible experiment configurations
- `src/boost/`: time-dependent bond boosting implementation
- `src/reactive/`: reactive group definitions and candidate selection
- `src/workflows/`: polymerization and curing loops
- `src/backends/`: MLIP/calculator adapters
- `scripts/`: entry points, checks, reporting, figure generation
- `.claude/agents/`: specialized sub-agents
- `.claude/skills/`: reusable workflows
- `.claude/hooks/`: local safety and quality hooks

## Default workflow
1. Read the paper and update `paper/claims.yaml`.
2. Translate a claim or equation into a task in `specs/tasks.md` or `specs/decisions.md`.
3. Implement the smallest viable slice.
4. Add or update unit and integration tests.
5. Run a minimal end-to-end reproduction.
6. Save artifacts in `runs/`.
7. Regenerate figures from scripts.
8. Update docs.

## Coding rules
- Python 3.11+
- Type hints required for public interfaces
- Configs should be structured with `dataclass`, `pydantic`, or typed dictionaries
- Functions should align with scientific concepts, not premature abstraction layers
- Use deterministic seeds everywhere possible
- Separate numerical kernels from orchestration code
- Prefer NumPy, PyTorch, OpenMM, ASE-style abstractions only when justified
- Keep I/O, simulation state, and analysis cleanly separated

## Scientific implementation rules
- Treat TDBB equations as the most paper-critical component.
- Preserve units explicitly in code and config comments.
- Log both biased and unbiased segments.
- Distinguish bond formation bias from bond dissociation bias in APIs.
- Record all reactive-pair selection decisions for debugging and auditability.
- Make reaction candidate generation testable without running MD.

## Ask-first triggers
Stop and ask before proceeding if:
- The paper supports multiple plausible mathematical interpretations.
- A backend with unclear commercial rights is about to be added.
- A simplification changes the scientific meaning of TDBB.
- A figure is being reproduced with altered smoothing, filtering, or averaging.
- A default hyperparameter is introduced without a paper citation or explicit user approval.

## Output expectations
When working in this repo, prefer outputs in this order:
1. Updated spec or decision record
2. Minimal code change
3. Tests
4. Reproduction command
5. Artifact description

## Definition of done
A task is done only when:
- The code is implemented
- Tests pass
- Assumptions are documented
- Reproduction commands are written
- Outputs can be traced to configs and seeds
