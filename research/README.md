# ResSim research programme

**Index:**
- `CHARTER.md`: method, breakthrough criteria, limits.
- `GAP_MAP.md`: literature gaps.
- `HYPOTHESES.md`: ranked hypotheses.
- `FINDINGS_LOG.md`: **all findings, publication-grade**.
- `CYCLE*_REPORT.md`: one report per research cycle.
- `LIT_WATCH.md`: recent-paper watch. Run with `python3 lit/watch.py`, which uses OpenAlex.
- `papers/P1_convex_seismicity/`: manuscript draft and publication plan.
- `exp/`: reproducible scripts. Each prints its numbers and writes JSON and PNG output.

**Reproduce everything (≈ 3 minutes, Python with numpy, scipy and matplotlib):**
```
cd exp
python3 h1_seis_opt.py h1_results.json
python3 h1_verify.py
python3 h1_dose.py
python3 h1_fig.py
python3 h1b_poro.py
python3 h1b_extra.py
python3 h1b_fig.py
```

**Cycle 4 (C22):** `python3 h1f_theta_law.py` (≈ 40–60 min on 2 cores: 30 SLSQP optima on refined grids) and `python3 h1g_kernel.py` (seconds). Results: `CYCLE4_REPORT.md`, F12–F13.

**Cycle 5 (C23):** `python3 exp/h6_aspo.py` (seconds; data in `data/aspo/`). Pre-registration `CYCLE5_PLAN.md`, results `CYCLE5_REPORT.md`, F14–F16 (Äspö HRL second site: volume-driven AE, no resolvable rate effect, cyclic reduction not significant). Engine: `seis::fitStageKernels`, rate-and-state fault cycles `engine/include/ressim/rsf.h` (G3).
