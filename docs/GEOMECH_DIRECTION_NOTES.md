# Deep geomechanics / fracture direction: assessment notes (after C18)

Purpose: judge the external suggestions (EFEM or ghost-penalty XFEM/CutFEM; VCCT instead of the interaction integral), list better
options, state the mathematical limits honestly, and propose how ResSim could grow from a reservoir simulator into an integrated
geological / geothermal analysis tool. These are working notes for a research round. They are not a literature review yet:
references are given by author and topic and must be checked before they are quoted.

## 1. Verdict on the suggestions

| Suggestion | Verdict | Why |
|---|---|---|
| EFEM instead of XFEM (phase 1-2) | **Agree, strongly.** | EFEM uses the same cut-cell geometry as our EDFM. The displacement jump lives in each element and is condensed out statically, so there are no new global unknowns and the small-cut conditioning problem is much milder. LLNL's GEOS pairs EDFM flow with EFEM mechanics in this way (Cusini, White, Castelletto, Settgast). Limits: one constant jump per element means no tip singularity, an aperture that is piecewise constant along the fracture, and possible stress locking for kinked or curved fractures. EFEM is right for opening and slip of existing fractures and faults at reservoir scale, but not for accurate propagation. |
| Ghost-penalty XFEM (CutFEM) | Correct remedy, but **only if we really use XFEM**. | The ghost penalty (Burman et al.) does fix the ill-conditioning from small cuts. The price is face-coupling terms (face loops on hex and corner-point grids) and penalty parameters to tune. Cheaper alternatives: aggregated FEM (Badia, Verdugo et al.), merging small cut cells, or a robust block preconditioner. With EFEM the problem largely goes away. |
| VCCT instead of the interaction integral (phase 3) | **The suggestion is backwards for embedded cracks.** | VCCT needs nodes on the crack front and elements aligned with the crack. That is natural with cohesive / interface-element meshes, but awkward when the crack cuts elements arbitrarily (XFEM/EFEM). The domain interaction integral is the mesh-independent choice in 3-D. **More important:** on reservoir grids (5-50 m cells) *any* stress-intensity extraction is crude. Hydraulic fractures are mostly viscosity- or leak-off-dominated, not toughness-dominated, so the tip is governed by the asymptote w ~ x^(2/3) and not by K_I. |

**Better replacement for phase 3:** front tracking from the tip asymptote. Take the universal tip asymptote (Dontsov and Peirce),
which spans the toughness, viscosity and leak-off regimes, and invert it in the cells next to the front. This gives the distance
from each cell to the front and the local front velocity (implicit level set algorithm, ILSA, Peirce and Detournay). The method
was designed for coarse grids, so growth becomes nearly grid-independent. Our HFGROW criterion K_I = p_net·sqrt(pi·min(L, H/2))
would become a special case (the toughness limit).

## 2. A better overall architecture: multi-fidelity with one shared fracture object

| Level | Scale | Method | Gives | Main limitation |
|---|---|---|---|---|
| A: field / reservoir | 100 m - 10 km | EDFM flow + **EFEM** fracture/fault mechanics on the existing hex FE grid; frictional contact (augmented Lagrangian / semismooth Newton); thermal stress; rate-and-state faults (quasi-static + radiation damping) | aperture from mechanics, slip, dilation, permeability evolution, induced-seismicity proxies, compaction | fault nucleation length is not resolved on coarse cells (see §4) |
| B: hydraulic fracture | 1 m - 1 km | **planar-3D displacement-discontinuity BEM** (only the fracture is meshed; elasticity is exact in a homogeneous medium) + lubrication FV + ILSA / universal tip asymptote + leak-off coupled to the reservoir through EDFM; H-matrix compression | high-quality HF geometry (height growth, multi-cluster stress shadow, proppant later) | planar or mildly non-planar; layered moduli only approximately |
| C: near-wellbore / lab | mm - 10 m | **phase-field fracture** (variational, Bourdin-Francfort-Marigo; hydraulic versions by Wheeler-Wick, Miehe) on a fine FE mesh | crack initiation, branching, merging, curved paths, interaction with natural fractures, core / lab test analysis | cost (mesh h << length scale l), dependence on l, approximate crack opening |

The key design choice is **one fracture object shared by all three levels** (geometry, aperture, permeability, history). Level B
grows it and writes the geometry into EDFM. Level A carries the slip and aperture of all fractures and faults. Level C runs as a
sub-model with boundary stresses taken from A. Commercial tools usually have only one of these.

## 3. "COMSOL / MOOSE / TOUGH quality" for EGS: what it means in practice

Quality has to be shown by benchmarks, not claimed:

- **Coupled THM:**
  - DOE Geothermal Technologies Office code comparison problems (White et al., Geothermics 2018);
  - Mandel and Cryer (have);
  - Terzaghi (have);
  - Lauwerier (have).
- **Fracture mechanics:**
  - Sneddon crack opening;
  - a sliding crack with friction (analytical);
  - KGD / PKN / radial HF in all regimes (Savitski and Detournay; Detournay's review);
  - multi-cluster interaction (have a PKN version).
- **Induced seismicity:**
  - SCEC SEAS benchmarks for rate-and-state (BP1-QD and later).
- **Field data:**
  - Utah FORGE stimulations (pressure records, microseismic clouds);
  - Bedretto;
  - Soultz.

Physics still missing for EGS:

- thermal stress in fracture opening and slip;
- rate-and-state fault slip with seismic-moment output (Gutenberg-Richter statistics from the cloud);
- shear-stimulation permeability on a DFN coupled to mechanics (we have the Barton-Bandis / Coulomb proxies);
- silica / calcite chemistry in fractures;
- traffic-light / adaptive injection control.

## 4. Mathematical limitations to state up front

1. **Tip singularities.** In LEFM stress ~ r^(-1/2). In a viscosity-dominated HF, w ~ x^(2/3) and the pressure has a log
   singularity. Lubrication is degenerate (conductivity ~ w^3), and the tip is a moving boundary that may have a fluid lag. Coarse
   grids cannot resolve the tip, so we need asymptotes or enrichment, never plain refinement.
2. **Strain-softening plasticity loses ellipticity**, so localisation becomes mesh-dependent. The C16 Drucker-Prager softening is
   affected unless it is regularised (viscoplastic, nonlocal / gradient, or fracture-energy scaling of the softening modulus).
   **Audit item for the next round.**
3. **Rate-and-state faults** need the nucleation length h* ≈ G·D_c / ((b − a)·σ'_n) to be resolved by several elements. Typical
   values are metres to tens of metres; reservoir cells are often coarser. On coarse grids only Coulomb / quasi-static slip is
   meaningful, so event magnitudes need sub-grid or BEM fault models.
4. **Coulomb contact** is non-smooth (complementarity) and non-associative. The Jacobian is non-symmetric and solutions may be
   non-unique. This needs semismooth Newton or an augmented Lagrangian, plus a careful time-step choice.
5. **Fixed-stress split** convergence is proven for linear poroelasticity (Mikelić and Wheeler). With contact, plasticity or
   fracture growth it can stall, so we need Anderson acceleration or a monolithic fallback.
6. **DFN to continuum upscaling** is valid only above a representative volume, which fractured granite may not have at cell scale.
7. **Phase-field fracture** only converges to sharp cracks as l → 0 with h << l. The nucleation stress depends on l unless a
   strength-consistent model is used (Kumar, Bourdin, Francfort, Lopez-Pamies).
8. **Inverse problems** (stress state, fracture network from microseismic data) are ill-posed. Report uncertainty, not single
   answers.

## 5. Candidate novel contributions (honest; none is guaranteed to be a breakthrough)

- **Adjoint-based injection design with seismic-risk constraints:**
  - gradients of the induced seismicity measures (Coulomb stress rate, rate-and-state seismicity rate) and of heat recovery with
    respect to the injection schedule, through the coupled THM model;
  - we already have the adjoint framework for all three engines, and very few tools offer this.
- **A sub-grid tip model inside EDFM:**
  - the tip-asymptote front location within a host cell, so that reservoir-scale growth is grid-independent;
  - this would close the gap between reservoir simulators and dedicated HF simulators.
- **One fracture object across EDFM, EFEM, BEM and phase field** in one desktop tool with a UI: an integration result rather than
  a new algorithm, but rare.
- **Bayesian stress-state inversion** (DFIT, breakouts, microseismic) feeding all of the above with uncertainty.

## 6. Proposed process

1. **Research round** (mostly reading and small prototypes, little production code). Output:
   - a literature review with a method-comparison table;
   - decisions on the architecture in §2;
   - the benchmark list from §3, each with an analytic or published reference;
   - three proof-of-concept checks:
     - 2-D DDM KGD vs the analytical regimes;
     - EFEM sliding crack vs the analytical solution;
     - a mesh-dependence check of the current DP softening, with a regularisation fix.
2. **Implementation rounds:**
   - EFEM + contact (A);
   - planar-3D DDM HF with ILSA (B);
   - rate-and-state faults + seismicity outputs (A);
   - phase-field sub-model (C);
   - adjoint injection design.
   Every round is verified against §3.

## 7. C20 revision after implementing G1 (critical update)

Section 1 recommended element-constant EFEM. **The implementation and its tests revised that recommendation** (research finding F9):

* With one constant jump per element (statically condensed), opening and slip are accurate only on meshes aligned with the fracture.
  For a 30° crack the slip was over-predicted by 81 % / 36 % / 20 % at a/h = 3.5 / 7 / 14 (and still 25 % → 20 % at a/h = 12 → 19 for a longer crack), the Sneddon volume error *grew* from 10 % to 12 % with refinement, a spurious tangential jump of about 25 % of the opening did not decay, and a locked frictional case failed to converge.
  This agrees with the ISRM 2015 EFEM-vs-XFEM comparison, Cusini et al. 2021 and the GGGR 2024 fault-slip benchmark.
* Adopted instead: **continuous jumps with nodal Heaviside enrichment** (the XFEM jump without tip functions; the tip edge carries zero jump).
  This costs 2 extra global unknowns per enriched node (banded, still cheap in 2-D). Errors are 1–2 % (opening) and 5 % (slip) at a/h = 12, independent of the orientation.
* The constant-jump model stays in the code (`jumpModel = 1`) for comparison and for very coarse reservoir-scale grids where speed matters more than accuracy.
* Implication for 3-D (G1b): use the same continuous enrichment on the Q8 mesh. The fracture–EDFM cut geometry is still shared, which was the main argument for EFEM.


## 8. C25 decision: the non-XFEM path (implemented)

Chosen alternative to XFEM: **mesh = flow grid (isoparametric hexahedra for corner point) + EDFM fracture pieces with an
LEFM tip criterion against the FE closure stress + Coulomb/Dieterich/RSF fault layer driven by the FE stress change.**
Every component is verifiable on its own (patch test, brick equivalence, KGD/PKN growth, closed-form thermo-elastic
unclamping), it runs in black oil, compositional/CO2 and thermal, and it works on corner-point grids where XFEM cuts
would create slivers. What it gives up (documented in ENGINES_AND_MATH §6.19): slip discontinuities in the FE mesh,
curved fracture paths, FE apertures. Next steps if needed: implicit fixed-stress for BO/comp, aperture from the FE
normal stress, phase field (G4) as the curved-path option inside the run.
