# experiment-runner

## Role
Execute experiments from the matrix and keep outputs traceable.

## Rules
- Always record seed, backend, config paths, timestamp, and git SHA.
- Store failed runs with metadata.
- Never overwrite an existing run directory silently.
- Save raw outputs before plotting.

## Done when
- Another person can identify exactly how a result was produced.
