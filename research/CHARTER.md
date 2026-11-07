# ResSim Research Programme — Charter

**Mandate (user, C19):** look for breakthroughs anywhere: new capabilities, faster algorithms, new findings.
The assistant both plans and executes the work.

## Method (one cycle)
1. **Gap scan:** literature search; record each gap with a citation in GAP_MAP.md.
2. **Hypotheses:** write each one so it can be falsified. Rank by novelty × impact × testability ÷ cost. Record in HYPOTHESES.md.
3. **Pre-register:** before running anything, write down the success and kill criteria plus the verification tests.
4. **Experiment:** use the smallest model that can decide the question, as reproducible scripts in exp/ with logs and a JSON results file.
5. **Verify independently:** do not reuse the same derivation. Use a second formulation, an analytic limit, FD vs adjoint, convexity or KKT checks, and multiple starting points.
6. **Rethink:** state what is known, what is new, what is assumed, and what could break it. Then promote the result: to the engine/UI, to a paper draft, or kill it.

## Breakthrough criteria (all must hold)
- **(a) Novelty:** not found in a targeted literature search. The search is logged, and the claim is always worded "not found", never "does not exist".
- **(b) Correctness:** passes at least 2 independent checks.
- **(c) Impact:** gives a new capability, ≥ 2× speed, or changes a decision a practitioner would make.
- **(d) Reproducibility:** one script, under 10 minutes on 2 cores and 2 GB.

## Resource limits
2 CPU cores, 2 GB RAM, no GPU, no commercial simulators. Every experiment must fit within these limits, or be designed to run later on the user's machine.

## Honesty rules
- Model results are not field results; every claim states its assumptions.
- Negative results are recorded too.
