experiments:
  - id: toy_tdbb_smoke
    objective: Verify TDBB plumbing without a proprietary backend
    backend: toy
    config:
      boost: configs/boost/paper_faithful.yaml
      eval: configs/eval/smoke.yaml
    seeds: [7]

  - id: paper_schedule_sanity
    objective: Validate biased/unbiased alternation schedule
    backend: toy
    config:
      boost: configs/boost/paper_faithful.yaml
      eval: configs/eval/smoke.yaml
    seeds: [7, 11, 19]
