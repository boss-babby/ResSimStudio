# What drives acoustic emission during mine-scale hydraulic fracturing? A stage-level kernel test on the Äspö HRL experiments with the continuous catalogue and corrected event times

*Manuscript draft v1.0, 2026-11-04. Target: Geophysical Journal International (alternatives: Solid Earth; Seismological Research Letters). Authors: [to be decided]. Generative-AI assistance was used for code, analysis and drafting; disclose per journal policy. AI tools are not authors.*

> **Status (honest, for the authors).** The analysis is complete and runs in about 30 s (`research/exp/h6c_aspo.py`; figures
> `scripts/mkfigs.py`). The pre-registered predictions (CYCLE12_PLAN, H6c) were all supported on the primary catalogue. The
> paper is a quantitative reanalysis of public data. Its contribution is the formal model comparison, the correction of a timing
> defect and a calibrated statement of what the data can and cannot show about cyclic injection. It must be framed as building on
> Zang et al. (2017), Kwiatek et al. (2018) and Niemz et al. (2020), not as overturning them. Before submission, contact the data
> authors (GFZ) about the timing defect (§2.2) and offer them the opportunity to comment.

## Abstract

The 2015 Äspö Hard Rock Laboratory experiments created six hydraulic fractures from one borehole at 410 m depth using
conventional, progressive-cyclic and pulse injection. They are a key data set for the claim that cyclic "soft" injection reduces
induced seismicity. We reanalysed the public hydraulic records together with the continuous acoustic-emission (AE) catalogue of
Niemz et al. (2020), which has 19,756 events, at the level of 31 injection stages. Five kernels for the expected AE count of a
stage were compared, fitted on the two conventional tests in the same rock:
* hydraulic energy;
* injected volume;
* a Dieterich rate-and-state response to the interval pressure;
* stage duration;
* a constant.

In the magnitude-complete catalogue (M_AE ≥ 2.5, 2,322 events), hydraulic energy is preferred (ΔAIC = 53 over volume, 545
over Dieterich, 736 over duration). The ranking is unchanged for completeness thresholds from 2.15 to 3.0, for the full catalogue
and for the 196-event triggered catalogue. A stage bootstrap never gives the Dieterich kernel an advantage. At fixed volume and
fixed cumulative volume, higher-rate stages produce more AE (exponent 0.4–0.9). Only 14 % of stage AE occur after shut-in,
against 50 % predicted by the Dieterich kernel, although the interval pressure persists for minutes.

The progressive-cyclic test HF3 produced 0.25–0.52 times the energy-normalised AE rate of the conventional tests in the same rock,
consistently across catalogues. However, with three or four conventional tests the between-test scatter is too large for the
reduction to be significant at the 5 % level, except in one of eight catalogue/scatter combinations (one-sided p = 0.04–0.31).

We also report a defect in the public catalogue files: the seconds field of `origin_time` contains the day of the month. Correcting
it does not change the conclusions of a first analysis that used the defective times, but it shifts individual numbers.

## Key points

1. Stage-level AE at Äspö follow hydraulic energy, not a Dieterich response to pressure. The result is robust to the
   completeness magnitude, the catalogue and a stage bootstrap.
2. The AE stop with the flow: the post-shut-in fraction is 0.14, against 0.50 predicted by Dieterich. Rate matters at fixed
   volume.
3. Cyclic injection (HF3) shows a consistent 2–4× reduction that the available number of tests cannot establish statistically.
   The public catalogue times contain a defect, which we document and correct.

## 1. Introduction

Cyclic or "soft" injection has been proposed as a way to stimulate permeability with fewer and smaller induced events (Zang et
al., 2017, 2018). The Äspö HRL experiments of 2015 are one of very few in-situ tests in which the injection protocol was varied
under otherwise similar conditions. Six hydraulic fractures (HF1–HF6) were propagated from a horizontal borehole at 410 m depth
in three rock types, with conventional continuous injection, progressive cyclic injection (HF3) and pulse hammer injection (HF5).

Earlier studies reached the following conclusions:
* Zang et al. (2017): fluid-injection protocol and monitoring set-up.
* Kwiatek et al. (2018): from a triggered catalogue of about 200 high-quality events, seismic energy correlates with hydraulic
  energy.
* Niemz et al. (2020): from a continuous, full-waveform-based catalogue of about 20,000 events, HF3 has a lower seismicity rate
  per hydraulic energy, a larger b-value and a smaller maximum magnitude than the conventional tests.
* Niemz et al. (2021): reduced fracture aftergrowth in HF3, from tilt data.

Two questions remain open, and both matter for how the experiments are used in induced-seismicity models.
1. **Which physical driver best predicts the AE count of an injection stage?** Rate-and-state seismicity models (Dieterich, 1994),
   as used in many forecasts of fluid-induced seismicity, are driven by stress (here pore pressure). Fracture-mechanics arguments
   instead point to injected volume (McGarr, 2014) or hydraulic energy. A formal comparison on the same data has not been made.
2. **How certain is the HF3 reduction?** It comes from one cyclic test compared with a few conventional tests, which differ in
   rock type and duration. The statistical uncertainty that the between-test scatter implies has not been quantified.

We address both questions with a stage-level Poisson model. The hypotheses were stated before the continuous catalogue was
analysed (Section 2.4).

## 2. Data and methods

### 2.1 Data

The hydraulic data (flow rate and interval pressure) and both AE catalogues are from the public archive of Zang et al. (2023),
described by Zang et al. (2024):
* the triggered catalogue of Kwiatek et al. (2018), 196 events with moment magnitude;
* the continuous catalogue of Niemz et al. (2020), 19,756 events with AE magnitude M_AE for 4,302 of them.

Hydraulic records were reduced to 1-s means. The continuous catalogue contains no events on 10 June 2015, the day of the pulse
test HF5. HF5 is therefore excluded from all fits (Figure 1).

### 2.2 A timing defect in the public catalogue files

In both catalogue files the seconds field of `origin_time` equals the day of the month. Values are confined to {03, 04, 09, 11,
12}, which are the test days, in 196/196 and 19,756/19,756 events. The event name of the continuous catalogue
(YYYYMMDDHHMMSS.µs) is monotonic in time and self-consistent, whereas `origin_time` is not monotonic. We therefore used the event
name as the origin time. For the triggered catalogue we used the name of the matched continuous event; 4 of 196 names could not be
parsed and were dropped.

A first version of this analysis (Section 4.4) used the defective times. Event times were then wrong by up to ±60 s, which
affects the assignment of events near the stage boundaries and the split between co- and post-injection events.

### 2.3 Stages and kernels

A stage is a continuous period with flow > 0.3 L/min, with gaps < 5 s bridged, a duration ≥ 10 s and a median interval pressure
≥ 1 MPa. The last condition excludes packer inflation and line filling. Each stage window extends to the next stage or to 600 s
after shut-in, whichever comes first. This gives 31 stages, 25 of them outside HF5, and assigns 2,315 of the 2,322 complete
events.

The expected AE count of stage i in test j is μ_ij = a_j · K_i, with a per-test rate factor a_j and one of the following kernels:
* **K2, hydraulic energy:** Σ q·p over the stage.
* **K3, injected volume:** Σ q.
* **K1, Dieterich:** the Dieterich (1994) rate response to the interval pressure, integrated over the stage window. Aσ is fitted.
* **K5, duration.**
* **K0, a constant.**
* **K4, an energy-Kaiser kernel:** energy above the previous maximum.

Exploratory multiplicative extensions test whether the rate matters at fixed volume:
* **K6:** V·q̄^β.
* **K7:** V·(1 + V_cum)^γ, where V_cum is the cumulative volume of earlier stages. This tests the confound that later stages
  have both higher rates and larger pre-existing fractures.
* **K8:** V·q̄^β·(1 + V_cum)^γ.

The kernels are fitted by Poisson maximum likelihood on the two conventional tests in the same rock (HF1, HF2) and compared by
AIC. Rate factors for all tests are then estimated with the preferred kernel.

**HF3 test.** We compare the HF3 factor with the mean of HF1 and HF2. The comparison uses a predictive t-test: the log factor of a
new conventional test is assumed to scatter like the conventional tests (HF1, HF2, HF4 and HF6, or without HF4), with df = n − 1.

**Post-shut-in fraction.** For HF1 and HF2 we compare the observed post-shut-in fraction with the Dieterich prediction.

### 2.4 Pre-registration

The analysis on the corrected and continuous catalogues was pre-registered (H6c, CYCLE12_PLAN in the archive) after the timing
defect was found, and before any stage-level count on the continuous catalogue was computed. The primary catalogue is (iii),
M_AE ≥ 2.5, the completeness magnitude of Niemz et al. (2020). These plan files are internal project records: they were not deposited with an external registry, so their dates are not independently time-stamped. They will be archived unchanged with the code so that readers can compare plans and results. Predictions:
* **P-a:** K1 worse than the better of K2 and K3 by ΔAIC ≥ 10.
* **P-b:** K2 and K3 better than K5 and K0 by ΔAIC ≥ 10.
* **P-c:** HF3 ratio < 1 but not significant at 5 %.
* **P-d:** the observed post-shut-in fraction is below the Dieterich prediction.

The robustness checks were also declared in advance:
* a stage bootstrap within tests (1000 resamples);
* leave-one-test-out;
* completeness thresholds 2.15 and 3.0.

## 3. Results

### 3.1 Kernel comparison (Figure 2, Table 1)

| catalogue | n | K2 energy | K3 volume | K1 Dieterich | K5 duration | K0 const. | K6 V·q̄^β | K7 V·(1+V_cum)^γ |
|---|---|---|---|---|---|---|---|---|
| triggered, times corrected | 192 | **77.7** | 83.2 | 134.5 | 159.4 | 177.6 | 71.5 | 72.2 |
| continuous M_AE ≥ 3.0 | 192 | **102.2** | 113.7 | 159.3 | 225.0 | 226.3 | 53.0 | 62.1 |
| **continuous M_AE ≥ 2.5 (primary)** | 2,322 | **246.4** | 299.4 | 791.4 | 982.4 | 1198.9 | 207.0 | 228.8 |
| continuous M_AE ≥ 2.15 | 4,146 | **234.9** | 266.4 | 899.7 | 981.4 | 1407.4 | 254.4 | 263.2 |
| continuous, all events | 19,756 | **1296.6** | 1308.7 | 3441.5 | 2400.4 | 4859.4 | 1014.2 | 976.2 |

*Table 1. AIC of the stage-count kernels (training tests HF1 and HF2; lower is better). Bold: best pre-registered kernel.*

Hydraulic energy is the best pre-registered kernel in every catalogue. In the primary catalogue the Dieterich kernel is worse by
ΔAIC = 545, and the stage bootstrap puts the 95 % interval of ΔAIC(K1 − best of K2/K3) at 70–981. All 1000 resamples exceed 10.
The Dieterich kernel is even worse than stage duration in the full catalogue. Energy beats volume in the primary catalogue
(ΔAIC = 53), with a bootstrap 95 % interval of ΔAIC(K2 − K3) from −84 to +2. Energy is preferred in most resamples, but the
interval just reaches +2, so we report the energy preference as strong but not decisive.

All four pre-registered predictions are supported on the primary catalogue. P-a and P-b are also supported on the full catalogue.

Figure 3 shows the stage-level fit. The energy kernel follows the observed counts across two orders of magnitude. The Dieterich
kernel compresses them: it predicts similar counts for stages of very different volume.

### 3.2 The rate matters at fixed volume

In the primary catalogue, K6 gives a rate exponent β = 0.45–0.65 (95 % profile interval). After controlling for cumulative volume
(K8), β = 0.4–0.9, and the cumulative-volume exponent alone (K7) is γ = 0.25–0.30. Higher-rate stages therefore produce more AE
per injected volume, even after the fracture-growth confound is accounted for. This is consistent with the energy kernel, which
weights volume by pressure, and with the energy-rate correlation of Kwiatek et al. (2018).

The exponent depends on the completeness threshold: β = 0.10–0.25 for M_AE ≥ 2.15, and 1.3–2.3 for M_AE ≥ 3.0. In the full
catalogue the sign reverses. We attribute the reversal to detection saturation during high-rate bursts and therefore base
conclusions only on complete subsets.

### 3.3 The AE stop with the flow (Figure 5)

In HF1 and HF2, 210 of 1,465 complete events (14.3 %) occurred after shut-in, against 50.0 % predicted by the Dieterich kernel
(binomial z = −27). The deficit appears in every catalogue version (14–29 % observed, 40–66 % predicted). The interval pressure
decays over minutes after shut-in (Figure 1), so a pressure-driven model must produce a long tail. The observed AE instead follow
the flow. Niemz et al. (2021) analysed fracture growth at Äspö jointly from seismometer-derived tilt and AE; a direct comparison of our stage-level post-shut-in fractions with their growth estimates is left to future work. The deficit indicates that AE are produced
mainly at the propagating fracture tip, not by pressure diffusion into a stressed volume.

### 3.4 Cyclic injection: a consistent but statistically unresolved reduction (Figure 4)

The energy-normalised rate factor of HF3, relative to the mean of HF1 and HF2 in the same rock, is:
* 0.25 in the triggered catalogue (times corrected);
* 0.47 in the full continuous catalogue;
* 0.41 at M_AE ≥ 2.5;
* 0.52 at M_AE ≥ 2.15.

The direction agrees with Niemz et al. (2020).

The significance of the reduction depends on how the between-test scatter is estimated:

| catalogue | scatter from HF1, HF2, HF4, HF6 | scatter from HF1, HF2, HF6 |
|---|---|---|
| triggered, corrected | σ = 1.15, t = −0.99, p = 0.20 | σ = 0.67, t = −1.72, p = 0.11 |
| continuous, all | σ = 0.32, t = −1.92, p = 0.08 | σ = 0.28, t = −2.19, p = 0.08 |
| continuous M_AE ≥ 2.5 | σ = 1.35, t = −0.55, p = 0.31 | σ = 0.23, t = −3.23, **p = 0.04** |
| continuous M_AE ≥ 2.15 | σ = 0.44, t = −1.20, p = 0.16 | σ = 0.23, t = −2.33, p = 0.07 |

*Table 2. One-sided predictive t-test of the HF3 rate factor against the HF1/HF2 mean (σ = between-test log scatter).*

HF4 was in a different rock and was stopped early. It has very few complete events (7) and dominates the scatter when it is
included. Leave-one-test-out changes the HF3 ratio by at most a factor of 1.3.

Our reading: every catalogue version shows a 2–4× reduction for HF3, but with only three or four reference tests, the data cannot
exclude that HF3 is simply a low-activity test. A definitive statement needs more replicate tests per protocol, ideally in the
same rock. Alternatively, a model of between-test variability could be calibrated on other experiments.

### 3.5 Effect of the timing correction

Using the defective times, the earlier version of this analysis found AIC values of 80.4 (energy), 86.3 (volume), 137.2
(Dieterich), 158 (duration) and 265.7 (Kaiser), with an HF3 ratio of 0.25 (t = −1.03) and a post-shut-in fraction of 0.25
against 0.48. With corrected times the same catalogue gives 77.7, 83.2, 134.5, 159.4 and 322.9, an HF3 ratio of 0.247
(t = −0.99) and 0.29 against 0.49. The qualitative conclusions are unchanged. Analyses at finer time resolution (inter-event
times, co-/post-injection splits within seconds) would be more strongly affected.

## 4. Discussion

### 4.1 Implications for seismicity models of stimulation

Dieterich-type models have been applied to fluid-induced seismicity at scales from the laboratory to the reservoir. At Äspö,
however, they fail at the stage level for two reasons. They predict the wrong distribution of counts across stages, and they
predict a long post-shut-in tail that does not occur.

The energy kernel's success, together with the positive rate exponent at fixed volume, suggests that AE production is controlled
by the work done on the propagating fracture. In this regime, schedule-optimisation frameworks based on Dieterich convexity (e.g.
our companion paper) do not apply: hydraulic-fracture AE are outside their scope. Conversely, shear reactivation of pre-existing
faults by pore-pressure diffusion, the setting of most reservoir-scale induced seismicity, is not tested by these data.

### 4.2 Implications for cyclic soft stimulation

Our results neither confirm nor refute the benefit of cyclic injection. They quantify the evidence: a factor 2–4 reduction,
consistent across catalogues, with a one-sided p-value between 0.04 and 0.31 depending on the analysis choices. The magnitude
statistics of HF3 (b-value, maximum magnitude) reported by Niemz et al. (2020) are a separate line of evidence that our count-based
test does not address.

A future in-situ campaign should include at least three replicate tests per protocol in the same rock to resolve a factor-2 effect
against the observed between-test scatter.

### 4.3 Limitations

1. **Stage definition.** Stages are defined with fixed thresholds; we did not explore alternative thresholds systematically.
2. **Pressure measurement.** The Dieterich kernel uses the interval pressure at the borehole, not the pressure on the fracture
   faces. A fracture-scale pressure model could improve K1, but it would add parameters.
3. **Completeness.** It varies with noise during injection. We addressed this with three thresholds, and the conclusions on
   kernels hold for all of them.
4. **Rock types.** HF4 and HF6 are in different rock types from HF1–HF3, so the rate factors combine protocol and rock effects.
5. **Stage independence.** The Poisson model treats stages as independent given the kernel. Overdispersion would widen intervals
   but would not change the ranking, which is supported by the bootstrap.

## 5. Conclusions

1. At the stage level, AE during the Äspö hydraulic-fracturing experiments are best predicted by hydraulic energy. They are
   predicted poorly by a Dieterich response to pressure (ΔAIC ≥ 545 in the complete catalogue).
2. AE stop with the flow (14 % after shut-in, against 50 % predicted) and increase with injection rate at fixed volume.
3. Cyclic injection (HF3) shows a consistent 2–4× lower energy-normalised AE rate. With the available tests this reduction is
   not statistically robust.
4. The public catalogue `origin_time` field contains a seconds/day defect. Users should take times from the event names.

## Open Research

Data: Zang et al. (2023), GFZ Data Services, https://doi.org/10.5880/GFZ.2.6.2023.004 (CC BY 4.0). The scripts reproducing all
results and figures (`h6_aspo.py`, `h6c_aspo.py`, `mkfigs.py`) and the pre-registration file will be archived on Zenodo [DOI to
be added].

## References

Dieterich, J. (1994). A constitutive law for rate of earthquake production and its application to earthquake clustering. *Journal of Geophysical Research*, 99(B2), 2601–2618. https://doi.org/10.1029/93JB02581

Kwiatek, G., Martínez-Garzón, P., Plenkers, K., Leonhardt, M., Zang, A., von Specht, S., Dresen, G., & Bohnhoff, M. (2018). Insights into complex subdecimeter fracturing processes occurring during a water injection experiment at depth in Äspö Hard Rock Laboratory, Sweden. *Journal of Geophysical Research: Solid Earth*, 123, 6616–6635. https://doi.org/10.1029/2017JB014715

McGarr, A. (2014). Maximum magnitude earthquakes induced by fluid injection. *Journal of Geophysical Research: Solid Earth*, 119, 1008–1019. https://doi.org/10.1002/2013JB010597

Niemz, P., Cesca, S., Heimann, S., Grigoli, F., von Specht, S., Hammer, C., Zang, A., & Dahm, T. (2020). Full-waveform-based characterization of acoustic emission activity in a mine-scale experiment: a comparison of conventional and advanced hydraulic fracturing schemes. *Geophysical Journal International*, 222(1), 189–206. https://doi.org/10.1093/gji/ggaa127

Niemz, P., Dahm, T., Milkereit, C., Cesca, S., Petersen, G., & Zang, A. (2021). Insights into hydraulic fracture growth gained from a joint analysis of seismometer-derived tilt signals and acoustic emissions. *Journal of Geophysical Research: Solid Earth*, 126. https://doi.org/10.1029/2021JB023057

Zang, A., et al. (2017). Hydraulic fracture monitoring in hard rock at 410 m depth with an advanced fluid-injection protocol and extensive sensor array. *Geophysical Journal International*, 208, 790–813. https://doi.org/10.1093/gji/ggw430

Zang, A., Zimmermann, G., Hofmann, H., Stephansson, O., Min, K.-B., & Kim, K. Y. (2018). How to reduce fluid-injection-induced seismicity. *Rock Mechanics and Rock Engineering*, 52, 475–493. https://doi.org/10.1007/s00603-018-1467-4

Zang, A., Niemz, P., von Specht, S., Zimmermann, G., Milkereit, C., Plenkers, K., & Klee, G. (2023). Comprehensive data set of in-situ hydraulic stimulation experiments for geothermal purposes at the Äspö Hard Rock Laboratory (Sweden). GFZ Data Services. https://doi.org/10.5880/GFZ.2.6.2023.004

Zang, A., Niemz, P., von Specht, S., Zimmermann, G., Milkereit, C., Plenkers, K., & Klee, G. (2024). Comprehensive data set of in situ hydraulic stimulation experiments for geothermal purposes at the Äspö Hard Rock Laboratory (Sweden). *Earth System Science Data*, 16, 295–310. https://doi.org/10.5194/essd-16-295-2024

## Figures

![Figure 1](figs/fig1_tests.pdf)

**Figure 1.** Hydraulic records (flow q, blue; interval pressure p, grey) and AE rate per minute (grey: all events of the
continuous catalogue; red: M_AE ≥ 2.5) for the six tests (first record of each test; HF6 continued on the next day). The
continuous catalogue has no events on the day of HF5.

![Figure 2](figs/fig2_dAIC.pdf)

**Figure 2.** ΔAIC of each pre-registered kernel relative to the best one, for five catalogue versions (log scale; the dashed line
marks ΔAIC = 10). Hydraulic energy is the best kernel in all versions.

![Figure 3](figs/fig3_obs_exp.pdf)

**Figure 3.** Observed against expected stage counts (M_AE ≥ 2.5) for the energy kernel (left) and the Dieterich pressure kernel
(right), with per-test rate factors. The grey band is ±2√μ. Zero counts are plotted at 0.3.

![Figure 4](figs/fig4_factors.pdf)

**Figure 4.** Energy-normalised rate factors per test relative to the mean of HF1 and HF2, for four catalogue versions. The
progressive-cyclic test HF3 is shaded.

![Figure 5](figs/fig5_postshutin.pdf)

**Figure 5.** Fraction of stage AE after shut-in in HF1 and HF2: observed (95 % binomial interval) against the prediction of the
fitted Dieterich kernel.
