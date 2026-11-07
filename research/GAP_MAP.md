# Gap map (Cycle 1 scan, 2026-10)

| Area | State of the art (evidence) | Gap? | Verdict |
|---|---|---|---|
| Nonlinear solvers (FIM) | OBL/RAJ linearisation (Computation 2025, 10.3390/computation13080191); localized Newton ≥13× fewer computations (arXiv 2008.01539); adaptive ASPEN↔Newton (JutulDarcy 2025); nonlinear elimination in fractures, −15% iterations (arXiv 2504.18900); ML initial guess −14%; NCPR; nonlinear multigrid | Generic tricks already published | **Low novelty**. Engineering value only (H3). |
| HF tip on coarse grids | ILSA + universal asymptote accurate on coarse meshes; LEFM needs ~10× more elements (Peirce–Detournay CMAME 2008; Dontsov–Peirce 2016; PyFrac); layered-stress asymptote (arXiv 2305.13714) | Sub-grid asymptote **inside EDFM growth** is rarely done | Medium (H2) |
| DDM + EDFM coupling | Published (JRMGE 2026, S1674775526001204) | Not new by itself | Implement as engineering only |
| Injection schedule vs induced seismicity | Contradictory claims: constant best (Boyet, Phil Trans A 2024); gradual step-up best (Alghannam & Juanes, Nat Commun 2020); cyclic best (PSU 2024); cyclic worse under volume control (2024); cycling lowers max M0 (Geng GRL 2026). Segall & Lu (JGR 2015): "more work required to determine optimal tapering" | **No proven global-optimum theory; no explanation for the contradictions** | **High (H1, Cycle 1)** |
| Seismicity-constrained optimisation | LP on Coulomb rate / seismogenic index (Hill et al. JGR 2024); stochastic gradients with 880 controls (Frontiers 2022, 10.3389/feart.2022.1053951); feedback control (arXiv 2310.02700); integral form of N(t) (Heimisson & Segall JGR 2018) | Nobody has stated that the Dieterich + linear-poroelastic problem is **convex** (targeted search, not found) | H1 |
| Stress inversion | Bayesian, from breakouts and DFIT | Integration with the simulator is missing | H4 (later) |

## Cycle 9 update (C47, 2026-11-02)

* **Closed:** whether injection-design convexity survives multiphase CO₂ flow with a coupled fault. It does
  (F25). The reason is that the rate → stress nonlinearity is small *and mostly favourable* compared with the Jensen
  margin (F26). Simulator-driven SCP reaches the global optimum in about 11 runs (F27).
* **Closed (negative):** a Dieterich surrogate as a design proxy for the seismic moment of a rate-and-state fault (F28).
* **Citation fixed:** Candela et al. 2022 (Dieterich-constrained injection optimisation, Oklahoma) was already in this map
  as "stochastic gradients, Frontiers 2022" (Cycle 1), but P1 did not cite it by name. It is now cited. The paper makes no
  structural or convexity claim, so the P1 novelty stands.
* **New gaps:**
  * two-way coupled faults (slip-dependent permeability), where NL could become harmful (H10e);
  * gas-compressibility-aware sub-cell storage for EDFM leak-off, since SGL over-corrects in gas at t90 (F30, H2e);
  * moment-aware schedule design (H11).

## Cycle 12 update (C50, 2026-11-04)

- **Closed:** F35 single-geometry gap (H10g, F37); explicit-lag concern (F38); Äspö small-catalogue gap (continuous catalogue is public, F39/F40).
- **Prior art re-checked for the manuscripts:** Hill et al. 2024 JGR (LP on stress-rate constraints, Raton Basin), Candela et al. 2022 (Dieterich optimisation, Oklahoma), Boyet et al. 2024 (Basel, permeability enhancement, constant injection best), Zheng/Jha/Jafarpour SPE-215155 (ensemble control CO₂), Niemz 2020 / 2021, Kwiatek 2018, Zang 2017 / 2018 (Rock Mech Rock Eng), Norbeck & Horne 2015 and Li 2020 (EFM leak-off; weakens P4 novelty). None tests convexity or gives a breakdown threshold; none does a formal Äspö kernel comparison.
- **Still open for P2:** 3-D multi-fault geometry; implicit two-way coupling; hysteretic permeability law; heterogeneous pre-stress.
- **Still open for P4:** KGD/Carter benchmark, fine-grid reference, comparison with Vinsome–Westerveld.
- **UHS / cyclic storage seismicity:** forward HM studies exist (e.g. 3-D FE fault-slip studies of seasonal H₂ cycles; EGU26 Hontomin analogue with RSF rates). No convex/optimal cycle design with a certificate found (one search, 2026-11-04) → candidate field (NEW_FIELDS.md).
