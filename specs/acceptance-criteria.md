# Acceptance criteria

## Phase 1: Framework correctness
- TDBB equations are implemented and unit-tested against paper-derived expectations.
- Formation and dissociation biases can be evaluated independently.
- Candidate selection enforces non-overlapping atom usage.
- The workflow alternates biased and unbiased segments.

## Phase 2: Reproduction readiness
- One toy polymerization example runs end-to-end.
- The run emits machine-readable metadata.
- A figure script can regenerate at least one publication-style plot from raw outputs.

## Phase 3: Scientific reproduction
- Paper-faithful config reproduces the expected qualitative trend.
- Deviations from the paper are explained in `specs/decisions.md`.
- Backend choice and commercial-use status are documented.
