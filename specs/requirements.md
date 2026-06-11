# Requirements

## Goal
Build a reproducible implementation of the paper's TDBB-driven polymerization/crosslinking workflow with commercially usable default components.

## In scope
- TDBB equations and schedules
- Reactive group representation
- Candidate-pair selection with non-overlap constraints
- Biased/unbiased simulation loop
- Backend-agnostic energy/force interface
- Experiment tracking and figure regeneration

## Out of scope
- Perfect numerical agreement before paper details are fully captured
- Proprietary backends as default public dependencies
- Production-scale distributed training from day one

## Functional requirements
1. The repository must support at least one end-to-end toy reproduction.
2. The TDBB module must expose separate formation and dissociation bias functions.
3. Reaction candidate selection must be testable independently from MD.
4. Every run must emit a manifest with config, seed, backend, and git SHA.
5. Figures must be generated from scriptable raw outputs.

## Non-functional requirements
- Deterministic config handling
- Strong logging
- Clear units
- Commercial-license review before backend activation
