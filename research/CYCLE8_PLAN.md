# Research Cycle 8 — plan (pre-registration)

Written before the analyses below were run. H7b had finished when this file was written; its design and verdict
thresholds are fixed in the header of `exp/h7b_rsf_ensemble.cpp`.

## H2c — sub-grid leak-off beyond the H2b deck family
* Decks: `exp/h2c_mkdeck.py`, made from `examples/HFGROW_INJECTION.DATA` like H2b (dx 10 / 5 / 2.5 m, SUBGRID and SUBGRID LEAK).
  * **CP:** a corner-point grid (COORD / ZCORN) dipping 2° along y, water injection.
  * **GAS:** dry-gas injection at 3000 sm³/d into the oil zone (GAS phase, PVDG / SGOF, GOC above the reservoir).
* Measure: t60 / t90 are the times at which a wing first reaches 60 / 90 m, the H2b definition (it reproduces the
  H2b value t60 = 25.49 d at 10 m).
* Verdict per variant:
  * SUPPORTED if the SUBGRID LEAK t60 spread over the three dx is < 1 d **and** smaller than the spread with
    SUBGRID only.
  * REFUTED if the SUBGRID LEAK spread is ≥ 2 d.
  * Anything in between is INCONCLUSIVE.

## H8b — Utah FORGE 16A(78)-32, April 2022 (third site for P1)
* Data:
  * Pumping: GDR 1379 (DOI 10.15121/1871203), 1 s treating pressure and slurry rate for stages 1–3.
  * Catalogue: GES v2 (GDR 1429, DOI 10.15121/1908927), status 2 events only.
* **Clock alignment (data issue, decided before any kernel fit).** With the pump clock read as MDT (UTC−6), every
  stage shows hundreds of events 20–50 min *before* pumping begins, and there is no injection then. The
  catalogue is in UTC. A pump clock set to UTC−5 (CDT, the service company's clock) removes the problem.
  * Primary analysis: shift the pump times by −60 min.
  * Robustness: the verdict must hold for shifts of −50 and −70 min.
* Completeness: M ≥ −1.5. This is the largest of the three per-stage max-curvature Mc values (−1.8 / −1.7 / −1.5).
* Windows: 5-min bins from pump start (shifted) to the last pumping + 3 h, or to the end of the catalogue for that
  stage.
* Kernels: Poisson counts per bin, one rate factor per stage.
  * **KV:** volume pumped in the bin, with a 1 % floor.
  * **KE:** hydraulic energy Σ q·p in the bin, with a 1 % floor.
  * **KVτ:** volume with exponential memory, ∫ q(s) e^{−(t−s)/τ} ds; τ is profiled (1 extra parameter).
  * **KD:** Dieterich rate with S(t) = surface treating pressure; Aσ is profiled (1 extra parameter).
  * **KC:** constant.
* Verdict (AIC summed over the 3 stages, on the primary shift):
  * SUPPORTED (volume-type) if the best of {KV, KE, KVτ} beats KD by ΔAIC ≥ 2.
  * REFUTED if KD beats all of them by ≥ 2.
  * Otherwise INCONCLUSIVE.
* Caveat recorded in advance: surface treating pressure includes friction that depends on the rate. The
  downhole pressure is not in the data set.
