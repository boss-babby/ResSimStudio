# ResSim Studio: engines, mathematics and solution methods

This is the technical reference for every simulation engine and numerical tool in ResSim Studio. For each one it lists
the equations, the primary unknowns, the discretisation, the nonlinear and linear solution methods, and the
correlations and models it uses. It is kept up to date round by round, alongside `ROADMAP.md`, which covers features and
status. Source paths are relative to `ressim/`.

Contents

1. [Common framework](#1-common-framework)
2. [Grid and discretisation](#2-grid-and-discretisation)
3. [Black-oil engine](#3-black-oil-engine-simulatorcpp)
4. [Compositional EOS engine](#4-compositional-eos-engine-compsimcpp-eoscpp)
5. [CO2 storage engine (CO2STORE)](#5-co2-storage-engine-co2brinecpp--compsimcpp)
6. [Thermal engine (STARS-like)](#6-thermal-engine-thermsimcpp-steamcpp)
7. [Linear algebra](#7-linear-algebra-linsolvecpp-amgcpp-bsrcpp)
8. [Nonlinear solution and time stepping](#8-nonlinear-solution-and-time-stepping)
9. [Wells and field management](#9-wells-and-field-management)
10. [EOR and PVT laboratory tools](#10-eor-and-pvt-laboratory-tools-eorcpp-eoscpp)
11. [Streamlines](#11-streamlines-streamlinecpp)
12. [Studies: sensitivity, uncertainty, history matching](#12-studies-sensitivity-uncertainty-history-matching-studycpp)
13. [External engines (OPM Flow, ECLIPSE)](#13-external-engines)
14. [Verification](#14-verification)
15. [Open items and candidate advanced methods](#15-open-items-and-candidate-advanced-methods)

---

## 1. Common framework

| Piece | File | What it is |
|---|---|---|
| Forward-mode AD `AD<N>` / `Dual<N>` | `ad.h`, `nonlinear.h` | Dual numbers with a fixed number of derivatives. The black-oil residual is written once for a generic scalar, and the exact Jacobian comes out of the same code (the ForwardDiff idiom from Julia). |
| `jacobian()`, `derivative()` | `nonlinear.h` | Value plus exact Jacobian of a small vector function. |
| `localNewton()` | `nonlinear.h` | Damped Newton for small dense systems: exact AD Jacobian, pivoted LU, Armijo backtracking, and box bounds that keep iterates physical. Used for local constitutive solves. |
| `OscillationDetector` | `nonlinear.h` | OPM-style detection of oscillating or stagnating Newton residuals. It triggers update relaxation (ω = 1 → 0.5). |
| `PidStepController` | `nonlinear.h` | PID time-step control on a relative change measure (Valli / Gustafsson type). |
| BSR block matrices | `bsr.h` | nb × nb dense blocks per connection. This is the native layout of every FIM Jacobian: black oil nb = 2–3, polymer 4, compositional nc + 1(+1), thermal nc + 3. |

Units: the engines use SI internally (Pa, K, m³, mol, J, s). Deck units (FIELD / METRIC) are converted at import and
export.

---

## 2. Grid and discretisation

- **Geometry:** Cartesian (DX/DY/DZ/TOPS) and corner-point (COORD/ZCORN) with ACTNUM, MINPV, MULTPV, NTG, PINCH,
  faults/MULTFLT, NNC, MULTREGT, NEWTRAN.
- **Transmissibility:** two-point flux approximation (TPFA), harmonic average of the half-transmissibilities
  T_i = K_i A_i·n_i / |d_i|² (corner-point: face-area vectors and centre-to-face vectors, the ECLIPSE/OPM convention),
  T = (T_a⁻¹ + T_b⁻¹)⁻¹ · multipliers.
- **Connection graph:** the solver is connection-based, so faults, NNCs, LGRs and dual-porosity links are simply extra
  connections.
- **LGR (CARFIN/REFINE):** sub-cells of the host hexahedron. Fine–coarse connections redistribute the host connection,
  which keeps fault throw and MULTFLT.
- **Dual porosity / permeability:** NZ doubled (matrix first). Matrix–fracture transfer T_mf = CDARCY·K·V·σ (SIGMA,
  SIGMAV, LTOSIGMA, MULTMF), and NODPPM for DUALPERM.
- **Phase upwinding:** single-point upstream weighting per phase, decided by the sign of the phase potential difference
  ΔΦ_α = p_a − p_b + Pc terms − ρ_α g (D_a − D_b), with ρ_α the saturation-weighted average.

---

## 3. Black-oil engine (`simulator.cpp`)

**Equations** (per cell, surface-volume conservation):

  ∂/∂t [φ (S_o b_o + S_g b_g R_v)] + ∇·(b_o λ_o v_o + R_v b_g λ_g v_g) = q_o
  ∂/∂t [φ (S_g b_g + S_o b_o R_s)] + ∇·(b_g λ_g v_g + R_s b_o λ_o v_o) = q_g
  ∂/∂t [φ S_w b_w] + ∇·(b_w λ_w v_w) = q_w

Here b = 1/B and λ = k_r/μ. Darcy flux: v_α = −T (ΔΦ_α).

- **Primary variables:** p_o, S_w, and a third variable chosen by a three-state switch per cell: S_g (three-phase), R_s
  (undersaturated oil) or R_v (undersaturated wet gas). The switch is applied at each Newton update.
- **Jacobian:** exact, by forward AD (`AD<N>`).
- **PVT:** PVTO/PVDO/PVTG/PVDG/PVTW/PVCDO tables. PVTG uses "right-extreme" interpolation that keeps the distance to
  the dew point (OPM rule). DRSDT (ALL/FREE) and VAPPARS (S_o/S_o,max)^v.
- **Saturation functions:** SWOF/SGOF and the SWFN/SGFN/SOF3 family. Three-phase oil kr: ECLIPSE default
  (saturation-weighted), STONE1 (+STONE1EX) and STONE2. End-point scaling: two-point and three-point (SCALECRS) for kr
  and Pc. SWATINIT Pc scaling.
- **Hysteresis:** Carlson scanning curves for the non-wetting phases (gas, and oil in oil–water). EHYSTR Killough is
  approximated by Carlson. Land trapping for the compositional core.
- **Equilibration:** EQUIL with hydrostatic integration of the phase pressures, plus RSVD/RVVD, THPRES and multiple
  EQLNUM regions.
- **Aquifers:** Fetkovich (analytic pseudo-steady state) and Carter–Tracy (van Everdingen–Hurst dimensionless pressure
  table).
- **Polymer (`polymer.cpp`):** OPM model. Todd–Longstaff mixing, Langmuir-type adsorption table, RRF permeability
  reduction, dead pore volume, shear thinning (PLYSHLOG).

---

## 4. Compositional EOS engine (`compsim.cpp`, `eos.cpp`)

**Formulation:** ECLIPSE-300 FIM / Acs–Coats–Watts volume balance.

- **Unknowns:** p, N_1..N_nc (+ N_w).
- **Component balances:** N_i − N_i^n + Δt (Σ_faces F_i + q_i) = 0
- **Volume balance:** PV(p) − V_L − V_V (− V_W) = 0
- **Flux:** F_i = T Σ_α [k_rα/μ_α · ξ_α · x_iα]_up ΔΦ_α

**Thermodynamics (`eos.cpp`):**

- **Equations of state:** Peng–Robinson (1976), PR78 m(ω) for ω > 0.49, and SRK. Van der Waals mixing with BIC and a
  Peneloux volume shift.
- **Fugacity coefficients:** from the reduced residual Helmholtz energy F(T, V, n) (Michelsen & Mollerup 2007). Exact
  analytic derivatives with respect to n, p and T.
- **Stability:** Michelsen tangent-plane distance with vapour- and liquid-like trial phases. Successive substitution,
  then a second-order method in α = 2√W.
- **Flash:**
  1. Rachford–Rice with a safeguarded Newton (Leibovici–Neoschil / Nielsen–Lia bounds).
  2. Successive substitution with GDEM acceleration.
  3. Full Newton on the vapour mole numbers (quadratic convergence).
  4. Warm start from the previous step's K-values; the stability test is skipped while a cell stays two-phase.
- **Absent components (z_i = 0):** trace partition K = φ_L/φ_V, giving an identity row in the flash Newton.
- **Derivatives of the flash:** implicit-function theorem (`flashSensitivities`), so the FIM Jacobian is analytic.
- **Phase properties:** Lohrenz–Bray–Clark viscosity (analytic derivatives), Li critical-temperature rule for labelling
  single-phase fluids, saturation pressure by Newton with a stability-test bisection fallback near critical.
- **Separators:** E300 FIELDSEP semantics (the last stage is the stock tank), WELLSTRE/WINJGAS streams, SEPCOND.
- **Water vaporisation (C7):** VAPWAT, written as `-- RESSIM: WATERVAP` because E300 EOS mode has no VAPWAT. The
  free-water Raoult law gives y_w = P_sat(T)/p (IAPWS-IF97). The gas phase takes y_w/(1 − y_w) moles of water per mole of
  hydrocarbon vapour at the gas molar volume, and the aqueous phase gives them up; when the free water runs out
  (dry-out), all of it is vaporised. Derivatives are exact forward-AD (dy_w/dp = −P_sat/p²), and water flows in the gas
  phase.
- **Water:** PVTW-style aqueous phase. CO2SOL: CO2 dissolves in the brine through f_CO2(EOS) = x_aq·H(p, T, m), with H
  calibrated to Spycher–Pruess / Duan–Sun. SOLUBILI tables (Rs_w(p)) are an alternative. Garcia apparent molar volume
  and Islam–Carlson viscosity.
- **Speed (C6):** in-place LU, OpenMP evalCell, AMG hierarchy reuse, inexact Newton (linear tolerance 1e-5),
  Aziz–Settari target-change time-step control, early divergence cut.

---

## 5. CO2 storage engine (`co2brine.cpp` + `compsim.cpp`)

- **Components and phases:** CO2 and H2O, with NaCl as a fixed molality. Phases are the aqueous phase and a CO2-rich
  phase.
- **Mutual solubility:**
  - Spycher, Pruess & Ennis-King (2003): Redlich–Kwong with a(T) = 7.54e7 − 4.13e4·T, b = 27.80 / 18.18 cm³/mol, true
    equilibrium constants and Poynting corrections.
  - Salting-out: Spycher & Pruess (2005) with the Duan–Sun (2003) activity coefficient, solved by a fixed point on the
    salt mole fraction.
- **CO2 phase:** RK molar volume and Fenghour–Wakeham–Vesovic viscosity.
- **Brine:** Batzle–Wang density and viscosity. Dissolved CO2 by the Garcia (2001) apparent molar volume.
- **Trapping mechanisms (C5):**
  - Structural trapping.
  - Residual trapping by Land trapping with hysteresis (Carlson / Killough).
  - Solubility trapping with molecular diffusion.
  - Salt dry-out (precipitation, with porosity and permeability reduction).
  - Mineral trapping by transition-state-theory kinetics: r = A k₂₅ exp[−E/R (1/T − 1/298.15)] (1 − Ω), with Palandri &
    Kharaka parameters.
- **Inventory:** six classes. Every molecule is counted exactly once.

---

## 6. Thermal engine (`thermsim.cpp`, `steam.cpp`)

STARS-style steam and hot-water processes: steamflood, cyclic steam stimulation (CSS) and SAGD.

**Unknowns per cell:** p, N_1..N_nc (hydrocarbon component moles), N_w (water moles) and U (internal energy of fluids +
rock, J), so nb = nc + 3.

**Equations:**

- **Components:** N_k − N_k^n + Δt (Σ F_k + q_k) = 0
- **Energy:** U − U^n + Δt (Σ [F_E + F_cond] + q_E + Q_loss) = 0
- **Volume:** PV(p, T) − V_W − V_O − V_G = 0, with PV = PV₀ (1 + c_r Δp + c_T ΔT)
- **Fluxes:**
  - Component flux: F_k = T Σ_α [k_r/μ · ξ · x_k]_up ΔΦ_α
  - Convective energy flux: F_E = T Σ_α [k_r/μ · ξ · h]_up ΔΦ_α
  - Conduction: F_cond = (A/d) λ_eff (T_a − T_b), with a harmonic average of λ.

**Why U instead of T:** with energy as the unknown, the state stays unique for univariant systems (pure water / steam,
or steam with a dead oil). There, T = T_sat(p) and the steam fraction is set by the energy. A T-based formulation needs
variable switching instead.

**(p, N, U) flash:**

- An energy search for E(T) = U over PT flashes, using an Illinois-modified regula falsi on a bracket (all quantities
  monotone in T).
- Univariant case: if U lies inside the latent-heat jump [E(T_s − ε), E(T_s + ε)], the phase split is interpolated in
  energy between the two flashed end states, with T = T_s.
- **New this round, near-univariant blend:** a trace of volatile hydrocarbon makes E(T) extremely steep but still
  continuous. If the temperature bracket collapses (|ΔT| < 1e-11·T) before the energy matches, the split is
  interpolated in energy between the bracket end states, using the true end-point energies rather than the
  Illinois-modified ones. The resulting state is smooth in U and p.
  - This fixed the steamflood failure. The finite-difference check had shown a wrong d(volume)/dp in the injector cell,
    because an energy error of slope × ΔT_bracket was larger than the local finite-difference step.

**Phase behaviour (free-water three-phase):**

- Phases are aqueous (W), oil (O) and gas (G = hydrocarbon vapour + steam).
- Water / steam equilibrium: y_w p = P_sat(T) (IAPWS-IF97).
- Hydrocarbons: K = (A + B/p + C·p) exp(−D/(T − E)) (STARS / E300 KVCR, Crookston form). Wilson K-values are the
  fallback.
- With K'_i = K_i/(1 − y_w), the oil/gas split is a standard Rachford–Rice problem. Without an aqueous phase, water is a
  gas-only component.

**Properties:**

- **Water and steam:** IAPWS-IF97 regions 1, 2 and 4 (v, h, u, c_p, P_sat, T_sat), with smooth metastable extrapolation
  for derivatives near saturation. Viscosity from IAPWS-2008, or a WATVISCT table. Steam viscosity from IAPWS or E300
  THSVC (A + B·T + C·p^D).
- **Hydrocarbon liquids:** ρ = ρ_ref exp(c_p Δp − c_T1 ΔT − c_T2/2 Δ(T²)) (DREF, CREF, THERMEX1); ideal mixing.
- **Hydrocarbon gas:** Z = ZFACTOR.
- **Enthalpies:** liquid h = SPECHA·ΔT + SPECHB/2·ΔT². Gas h = liquid h + HEATVAP·(1 − T/Tc)^HEATVAPE (Watson form), or
  SPECHG/H.
- **Viscosity:** OILVISCT per component with log-linear mole-fraction mixing; GASVISCT or GASVISCF (μ = A·T^B).
- **Rock:** H_r = HEATCR·ΔT + HEATCRT·ΔT²/2 per rock volume. Conductivity THCONR·(1 − THCONSF·S_g), or the
  THCROCK/THC* phase average.
- **Heat loss to the over- and underburden:** Vinsome–Westerveld semi-analytic profile per boundary face,
  θ(z, t) = (θ_i + p z + q z²) e^(−z/d), with the moments (I, θ_N) carried per face. It is fully implicit in the cell
  temperature.

**Discretisation and solver:**

- Local cell properties use forward differences of the local (p, N, U) state, with regime-consistent one-sided
  differences. Face, well and heat-loss terms use the chain rule.
- Block system solved by GMRES + BlockCPR.
- Convergence: mass 1e-6, energy 1e-3 K-equivalent, volume 1e-6.
- Time step: at most 14 Newton iterations, cut ×0.25. Growth follows the Aziz–Settari rule (1 + ω) η / (Δ + ω η), ω = 1/2,
  with targets ΔT 30 K, ΔS 0.2 and Δp 10 %, at most ×2, and at most ×1.25 for three steps after a cut (C8). Restart with a small step (0.05 d) after a well control change.
- **Robustness additions (C7).** Together these took the 5-year steamflood from a failure, or 15 minutes with 776
  cuts, down to about 80 s with 64 cuts.
  1. **Path-independent (p, N, U) flash.** The regula falsi bracket shrinks to a few ulps before blending, so the state no
     longer depends on the starting temperature. This removed residual noise of about 3e-6 relative in volume.
  2. **Appleyard-type local chopping in (p, N, U).** After the global update, every cell is re-flashed (properties only).
     If S_w or S_g would jump by more than 0.25, that cell's update is shortened (factor 0.25/ΔS, at least 0.3, up to 3
     times).
  3. **Oscillation / stagnation relaxation.** A two-cycle or stalled residual reduces ω (1 → 0.6 → 0.36 → 0.3, OPM
     style). ω grows back (÷0.6) once the residual drops by more than (1 − ω/2) per iteration (C8; the old
     "halves again" test could never pass at ω = 0.3).
  4. **Relaxed local tolerances from iteration 6** (mass 1e-5, energy 1e-2 K, volume 2e-5), accepted only when the global
     balances are tight: |Σ R_mass| / Σ N < 1e-9 and |Σ R_E| / Σ C < 1e-6 K. This follows OPM's CNV / MB split.
  5. **Well connection hysteresis.** A flowing connection shuts only when its drawdown reverses by more than 0.02 bar.
     This stops the on/off ping-pong of near-balanced connections; in CSS this ping-pong had made the cycle-3 restart
     fail.

**Wells:**

- Producers on BHP / oil / water / liquid control. Injectors on rate (CWE) with a BHP limit, or on BHP.
- **Producer steam trap (C9, STARS OPERATE STEAMTRAP).** The produced fluid is kept ΔT_sub below saturation. This is a
  dynamic lower BHP limit bhp ≥ P_sat(T_hot + ΔT_sub), where T_hot is the temperature of the hottest connection cell.
  - It combines with the rate / BHP controls; the most restrictive constraint is active.
  - It is implicit: d bhp/dX = P_sat′(T) · dT/dX of that cell enters the well-coupling Jacobian blocks.
  - The subcool is held to 1e-5 K when active, and the Jacobian check of the producer columns gives 8e-4.
  - Deck form: `-- RESSIM: STEAMTRAP well dT` in SCHEDULE (from that report step on). The SAGD example uses 10 °C.
- Injected enthalpy follows WINJTEMP semantics: quality with T, or quality with P (saturated); single-phase (T, P); or a
  given enthalpy.
- Mixture-density wellbore head, lagged.

**Deck:** E300 THERMAL keywords, two-way (import ↔ export). See `therm_deck.h` for the list.

### 6.1 C8 additions: reactions, solids, combustion, UCG (THMC "C")

**Component types.** Each component is a K-value fluid (KVCR / dead oil), a non-condensable gas (always in the gas
phase: O2, N2, CO2, CO, H2, CH4) or an immobile **solid** (coke, coal, char). Solids carry moles N_s per cell, take pore
space with molar volume MW/ρ_s, have their own heat capacity (SPECHS-type c_s) and never flow.

**Pore volume.** PV(p, T) = V_b φ_ref(p, T); the fluid pore volume is PV − V_s with V_s = Σ N_s MW_s / ρ_s. The volume
equation becomes PV − V_s − V_W − V_O − V_G = 0, so deposition (coke) or consumption (coal, char) of solids changes the
space available to fluids directly and implicitly.

**Reactions (STARS form).** For reaction r with reactant / product coefficients ν⁻, ν⁺ (nc + 1 entries, water last):

- rate (mol/m³ bulk/s): r = A exp(−E_a / RT) Π_k C_k^{o_k}, C_k = moles of k in its reacting phase per m³ bulk;
- **O2 partial-pressure kinetics:** for the flagged component the factor is its gas partial pressure y_k p (kPa) instead
  of a concentration (STARS O2PP);
- optional **minimum temperature** (rate 0 below), **reversible** form r (1 − Q/K_eq), ln K_eq = A_eq + B_eq/T (used
  for the water-gas shift);
- sources: component k gets V_b (ν⁺_k − ν⁻_k) r, energy gets V_b ΔH r (ΔH > 0 exothermic). The internal energy U
  remains a **sensible** energy, so the energy balance reads ΔU = H_in − H_out − Q_loss + ∫Q_react dt (closed to
  < 1e-7 in the tests).
- All rates are fully implicit (derivatives through the local (p, N, U) state), evaluated in the property vector
  Y at offset 4 + 3·stride. Mass is conserved by construction (Σ MW ν⁻ = Σ MW ν⁺ checked in the UI).

**Capillary pressure** in the thermal flux: Φ_W = p − P_cow(S_w), Φ_G = p + P_cog(S_g) from SWOF / SGOF.

**Conduction** uses the harmonic face average λ_f = d / (d_a/λ_a + d_b/λ_b) (E300 practice; arithmetic averaging over-
conducts across contrasts such as a coal seam against shale).

**Permeability alteration** from the fluid porosity φ_f = (PV − V_s)/V_b (solids and mechanics):

- Kozeny–Carman: k/k₀ = (φ/φ₀)³ ((1 − φ₀)/(1 − φ))²;
- exponential: k/k₀ = exp(C (φ − φ₀));
- clamped to [1/k_max, k_max]; applied to face transmissibilities (harmonic of the two cell multipliers) and well
  indices after each accepted step (explicit in time, like STARS PERMCK lagging).

**Gas injectors** inject a WELLSTRE stream (WINJGAS 'w' STREAM 's'): air, enriched air, O2/steam for UCG, CO2. Solid
fractions in a stream are dropped. The injected enthalpy is the ideal-gas enthalpy of the stream at WINJTEMP T.

**Examples:**

| Example | Physics | Result (2 cores) |
|---|---|---|
| `combustion` | Laboratory combustion tube 1.8 m, adiabatic, heavy oil; cracking to coke + heavy / light / coke burning (O2 partial pressure), air 30 sm³/m²/h, igniter at 400 °C | front at 1.2 m after 12 h (≈ 0.1 m/h), peak 630 °C, burned zone S_g = 1, steam plateau ≈ 150 °C ahead, oil bank; 282 steps, 8 s |
| `ucg` | Underground coal gasification, 25×1×3 seam with a bottom channel; pyrolysis, char combustion / partial oxidation, steam gasification, Boudouard, reversible water-gas shift, CO / H2 / CH4 burning; O2-enriched air | ignition, cavity growth (φ_f 0.05 → 0.6, k × 3800, Kozeny–Carman), syngas CO 0.7 / H2 0.3 / CH4 0.24 mol/s, cavity T ≈ 1350–1700 °C; 6 s |
| `sagd_thm` | SAGD with coupled geomechanics (§6.2) | surface heave 3.6 mm, reservoir top 19 mm after 91 d; 1 cm / 4 cm after 1.5 y |

The kinetic parameters of the examples are illustrative (literature-type values tuned for robust demonstrations), not a
calibrated reaction model for a specific oil or coal.

### 6.2 Geomechanics (THM) (`geomech.cpp`)

**Model.** Quasi-static linear thermo-poroelasticity, small strain, tension positive:
div σ = 0, σ = D ε − (α Δp + 3 K_dr β ΔT) m, ε = sym ∇u, m = (1 1 1 0 0 0)ᵀ,
with drained isotropic stiffness D(E, ν), Biot α, linear thermal expansion β, K_dr = E / (3 (1 − 2ν)).

**Discretisation.** Trilinear 8-node hexahedra with 2×2×2 Gauss integration on the Cartesian flow grid, padded with
non-reservoir over-, under- and side-burden elements (own E, ν, β). Rollers on the sides and base, free (or fixed) top.
Corner-point grids run on a box of the mean column sizes with an (i, j, k) cell map. The SPD system is solved by PCG with
symmetric Gauss–Seidel, warm-started from the last displacement (relative tolerance 1e-8).

**Coupling: fixed-stress split** (Kim, Tchelepi & Juanes 2011; unconditionally stable). Lagrangian porosity

φ = φ₀ + α ε_v + Δp / N − β_φ ΔT, 1/N = (α − φ₀)/K_s, β_φ = 3β (α − φ₀).

At fixed mean total stress the flow sees dφ = c_p dp + c_T dT with c_p = α²/K_dr + 1/N and c_T = 3αβ − β_φ. Each time
step: flow solve with the fixed-stress PV law → mechanics → re-centre the PV law on the mechanics porosity
(PV = V_b [φ_mech + c_p (p − p_k) + c_T (T − T_k)]) → repeat until max |φ_mech − φ_flow| < 1e-4 (≤ 4 iterations).
One-way coupling (mechanics after the step) is also available.

**Initial stress and failure.** σ_v0 = ∫ ρ_b g dz (gradient input), σ'_h0 = k₀ σ'_v0. Reported per cell: displacement
(heave), volumetric strain, effective mean / vertical / horizontal stress, Mohr–Coulomb shear index
τ_m / (c cos φ_f + σ'_m sin φ_f) and tensile index −σ'_3 / T₀ (capped at 10). Maps in Results: UPLIFT, DISPX, STRAINV,
SIGEFFM/V/H, FAILURE, POROM, POROF, PERMMULT, SOLID_<name>; summary FUPL, FUPR, FMXF, FRHR, FRHT.

**Boundary conditions (C9, BCMECH-like).** The six faces of the mechanics grid are set independently:
- 0 W, 1 E, 2 S, 3 N, 4 top, 5 base;
- each face is a roller (normal displacement fixed), free, or fixed (all components);
- a free face may carry a normal stress increment s (compression positive), e.g. a surface load or a tectonic stress
  change. It is applied as the traction −s n with consistent nodal loads s·A/4 per node of each bilinear face element;
- a rigid-body guard adds a roller when a direction has no constraint.

The deck form is `-- RESSIM: MECHBC face type stress`.

**Verification.** Uniaxial-strain column under uniform Δp and ΔT: heave H (α Δp + 3K β ΔT)/M, M = K + 4G/3. Porosity
and c_p are reproduced to 4e-10 (exact for linear fields on hex8). A 1 MPa surface load gives settlement s(H − z)/M,
Δσ'_v = s and Δσ'_h = s ν/(1 − ν), all to 1e-8.

### 6.3 C8 solver review (independent self-review findings and fixes)

| Finding | Fix |
|---|---|
| Arithmetic conduction average | harmonic (E300) |
| No capillary pressure in thermal | SWOF / SGOF P_c in the phase potentials |
| Porosity changes did not affect the rock-fluid volume (fixed PV, unused φ) | PV law with c_r, c_T, solids and fixed-stress re-centring |
| Recovery test of the relaxation (r₀ < 0.5 r₁) can never pass at ω = 0.3 (best reduction 1 − ω = 0.7) → linear convergence, cuts | test r₀ < (1 − ω/2) r₁ (SAGD 91 d: cuts 66 → 30, 25 → 19 s) |
| Time-step control η/Δ allowed ×2 growth straight after a cut (cut-double-cut cycles) | Aziz–Settari (1 + ω) η / (Δ + ω η), ω = 1/2; growth ≤ 1.25 for 3 steps after a cut |
| Global Newton damping (one difficult cell throttled the whole update) | per-cell damping of Δp and ΔU |
| Rate/BHP control could flip-flop between Newton iterations | at most two switches per time step, then the BHP limit is kept |
| IF97 above 1073 K undefined (combustion / UCG temperatures) | region 5 (to 2273 K) with a smooth blend |
| Injector wellbore density assumed liquid water | gas injectors: ideal-gas density of the stream |
| Rock energy with fixed grain volume | kept deliberately: grain mass is conserved, solids carry their own energy |

### 6.4 C10: unified multi-component framework (solubility, PR gas, trapping, inventory)

The goal is that a thermal / reactive case uses **the same phase-behaviour and trapping model** as the compositional
and CO2STORE cores. No separate module is involved: CO2, H2, CH4 and the other components are ordinary components with
library data (`gasdata.cpp`). That is how one run can produce UCG syngas and store CO2 at the same time.

**Component data.** The `gasdata` library holds 14 entries: H2, CO, CO2, CH4, C2H6, C3H8, N2, O2, H2S, Ar, He, NH3,
SO2 and H2O. For each it stores:
- critical constants (Poling–Prausnitz–O'Connell);
- IAPWS G7-04 Henry coefficients;
- the partial molar volume at infinite dilution v∞ (Garcia 2001 fit for CO2);
- the Sechenov coefficient k_s;
- the ideal-gas c_p;
- the element formula.

Names are looked up case-insensitively, with aliases (C1, METHANE, …).

**Aqueous solubility (Henry + Krichevsky–Kasarnovsky + Sechenov).** For a soluble component i:

  y_i φ_i p = x_aq,i · k_H,i(T) · exp(v∞_i (p − p_sat(T)) / RT) · exp(k_s,i m_NaCl)

- k_H follows IAPWS G7-04: ln(k_H/p_sat) = A/T_r + B τ^0.355/T_r + C T_r^−0.41 e^τ, with T_r = T/647.096 and τ = 1 − T_r.
  The Wagner–Pruss p_sat is used inside the formula.
- A van't Hoff form k_H = k_H298 exp(−ΔT (1/T − 1/298.15)) is available as an alternative.
- The analytic d ln k_H/dT gives the dissolution enthalpy h_sol = −R T² d ln k_H/dT.
- Verification: the G7-04 check values are reproduced to 5e-5 in ln k_H. CO2 in water at 50 °C and 8.4 / 15.4 / 24.4 MPa
  with the PR gas gives x = 0.0179 / 0.0211 / 0.0228. That is −2.5 / −5.4 / −4.6 % against Spycher–Pruess (2005).

**Three-phase flash: aqueous / oil / gas, convex formulation.**
- Every component has a fugacity constant per phase:
  - H_iG = φ_i p;
  - H_iO = K_i H_iG (KVCR K-values);
  - H_iA as above.
  - Water uses H_wA = p_sat (Raoult; K_w p above T_c) and H_wG = p.
- Phase j holds x_ij = f_i / H_ij. With affinities a_ij = min_k H_ik / H_ij, the phase fractions β minimise the convex
  function (Michelsen 1994 type)

  Q(β) = Σ_j β_j − Σ_i z_i ln(Σ_j β_j a_ij),  β ≥ 0

  The KKT conditions are exactly the phase-stability conditions: Σ_i x_ij = 1 for present phases and ≤ 1 for absent ones.
- Solver: a projected damped Newton method with an Armijo line search, Levenberg regularisation, and an active-set
  update. A bound variable whose Newton step points outwards leaves the free set. When the gas is non-ideal, an outer
  successive substitution on φ_i(y) runs to 1e-13.
- Absent phases report their trial compositions (useful for the 3D maps).
- The general flash is used only when a component is soluble. Otherwise the legacy free-water Rachford–Rice flash runs.
  The test suite checks that both give identical answers (to 1e-13) when nothing dissolves, including supercritical
  water.

**Peng–Robinson gas (`GASEOS PR`).**
- The components with critical constants form a PR sub-mixture. Constants come from the library by name or from
  `CRIT` lines; the deck TCRIT is never used.
- This sub-mixture provides φ_i (outer loop above), the molar volume (replacing its Z-factor volume), the residual
  enthalpy, and the LBC viscosity.
- Steam and non-PR components keep their own models (Amagat volume, mole-averaged viscosity).

**Aqueous phase with solutes.** Ideal dilute mixing:
- v = x_w v_w + Σ x_i v∞_i;
- h = x_w h_w + Σ x_i (h_i^ig + h_sol,i);
- u = h − p v;
- Islam–Carlson (2012) CO2 viscosity factor (1 + 4.65 x_CO2^1.0134).

The surface water volume counts only the water moles. Dissolved gas is reported as surface gas.

**Reactions with dissolved species.** `PHASE 0` on a reactant makes its concentration the aqueous one, n_W x_aq,i / V_b.
This is how mineral carbonation by dissolved CO2 works (CarbFix type). Example: CaO + CO2(aq) → CaCO3 in `ucg_co2`.

**Gas hysteresis (`HYSTGAS`).** This is the same model as the compositional core:
- Land trapping: S_gt = S_gc + (S_g,max − S_gc) / (1 + C (S_g,max − S_gc)).
- Carlson free-gas k_rg on the imbibition scanning curve.
- S_g,max is the historical maximum, lagged one step.
- Maps: SGMAX and SGTRAP.

**Storage inventory of the tracked component (`TRACK`, CO2 by default).** It uses the same classes as the CO2 cores,
and every molecule is counted once:
- dissolved in water;
- dissolved in oil;
- residual (gas up to S_gt);
- structural (mobile gas whose upstream outflow is below `structTol` × its volume per year);
- mobile;
- mineral, Σ ΔN_solid · trackEq · M (solids with a `FORMULA ... trackEq`, e.g. CALCITE trackEq 1).

The closed-cell test closes this inventory to round-off. The UI shows it in the "<name> storage" tab, the same tab as
CO2STORE runs.

**Element balance.** `FORMULA` lines (or library formulas) let the Reactions editor check the element balance of every
user reaction. Parentheses and fractional pseudo-formulas such as coal C0.85H0.7O0.1 are supported.

**Results.** Thermal runs write the following maps per report step, so the 3D viewer can animate how every component
moves (same arrays as the compositional core):
- ZMF*k* (overall, without water and solids);
- XMF*k* (oil);
- YMF*k* (gas);
- AMF*k* (aqueous);
- SGTRAP and SGMAX;
- SOLID_*name*.

The component selector of the viewer now uses the run's component list.

**Deck lines (SI).** These are carried in ECLIPSE comments. Long lines continue on `-- RESSIM+` lines; earlier exports
silently truncated reaction lines longer than 128 columns.
- `AQSOL name LIB | AQSOL name 1 A B C vInf ksalt | AQSOL name 2 kH298 dT 0 vInf ksalt`
- `GASEOS PR|IDEAL`
- `CRIT name Tc Pc ω`
- `SALINITY m`
- `HYSTGAS sgtMax|OFF`
- `TRACK name [structTol]`
- `FORMULA name formula [trackEq]`

**Supercritical water (C11).** The C10 limitation (water gas-only above T_c with a jump of K_w at T_c) is removed:

- **IF97 region 3** (ρ, T) is implemented in `steam.cpp`. Below T_c the density root is bracketed with the Wagner–Pruss
  auxiliary saturated densities ρ'(T), ρ''(T), so the root never jumps to the wrong branch near the critical point. At
  ≥ 623.15 K, `steam::liquid` / `steam::vapour` use region 3 between the B23 line and saturation (check values of the
  standard, worst relative error 6.9e-8; isobars 15–40 MPa have no jumps in h, ρ).
- **K_w and H_wA** continue C1 through T_c: the vapour-pressure branch is extended above T_c with a smooth quadratic, so
  the liquid / gas split of water is a continuous function of T at every p.
- **Water in a gas mixture near the critical point.** Ideal mixing at the partial pressure p_w puts the gas-water state
  onto a dense, liquid-like branch just above T_c, which makes U(T) of the cell non-monotone (two flash roots, so Newton
  oscillates). Water's enthalpy and internal energy in the gas are therefore evaluated at
  p_e = min(p_w, p_3(T, 0.75 ρ_c)) (`steam::vapourMixCap`), and u_w = h_w(T, p_e) − z(T, p_e) R T. The volume still uses
  the true p_w. In a U(T) scan of water + 0–30 % CO2 at 18–40 MPa, the negative slopes drop from up to 1.4e4 J/K to
  ≤ 150 J/K, and pure water is monotone everywhere.
- **Injector at a BHP limit.** Injection-only components can no longer flow back (backflow was clamped at N = 0, which
  stalled the time step). Both control modes use the same C1 no-backflow ramp, scaled by each well's connection width.
  The ucg_co2 variant at 800 sm3/d and 50 bar (which failed at day 62.6 in C10) now runs to the end.

**Known limitation (near-critical CO2–water cavity) — addressed in C12 by the fold construction (§6.6).** The ucg_co2 variant
at 150 sm3/d CO2 with a 120 bar BHP failed at day 244.3 in C11. The cause is one cell at p ≈ 22.86 MPa that sits on the 646–648 K fold of the mixture energy, where the
remaining non-monotonicity comes from ideal mixing of real-gas components. A proper fix needs a real-mixture EOS for
water + gases (PR with water, or a GERG-type model), or an equal-area construction across the fold in the energy flash
(treating the fold as a phase jump). This is planned for the next round. The 80 sm3/d and 800 sm3/d variants run.

### 6.5 C10 research notes: what EGS, H2 gasification and nanoparticle processes need from the engine

Each line below is a requirement for the engine, with its status.

- **UCG + CCS (brown / black / blue H2).**
  - Colours:
    - brown = lignite;
    - black = bituminous coal (≈ 19–20 t CO2 per t H2 unabated);
    - blue = with CCS (1–4 t CO2/t H2).
  - The engine needs the following, all **available in C10**:
    - (a) gasification kinetics (C8);
    - (b) CO2 retained in the cavity, in solution, as residual gas and as minerals;
    - (c) per-component product accounting (C9 Components tab: H2/CO, LHV);
    - (d) element balance of user chemistry.
  - Literature: CO2 storage in UCG cavities happens mainly through adsorption on char and coal, plus matrix–cleat
    diffusion. Overburden fractures are the leakage risk.
  - **Next:** a competitive extended-Langmuir sorption component on coal/char (CO2 VL 11.24 m³/t, PL 5.76 MPa; CH4 VL
    8.77 m³/t, PL 6.69 MPa; TOUGH2 EOS7C-ECBM). It is planned as an equilibrium "sorbed" state per solid with its own
    accumulation term.
- **ISCG / in-situ H2 from oil (STARS-type):** coke gasification + water-gas shift (~34 mol% H2 at 800 °C in the
  literature). This is expressible today with the C8 reaction framework. Use `TRACK H2` to get an H2 inventory.
- **EGS (enhanced geothermal).** Needs:
  - single-phase water / steam heat extraction with conduction to the rock (available);
  - THM with thermal contraction (available, §6.2);
  - discrete fractures — EDFM is the industry standard, MINC for natural fractures. **Not yet available:** dual
    porosity exists in the black-oil core only;
  - hydro-shear stimulation (Barton–Bandis closure, Coulomb slip from cooling and pore pressure) — **planned** on top of
    the geomechanics failure index;
  - Tracer and silica chemistry can be expressed with the reaction framework.
  - Reference codes: FEHM, TOUGH-FLAC, MOOSE; Utah FORGE data.
  - Gridding decision: EDFM fracture cells (non-neighbour connections already exist) and radial / LGR refinement
    around wells are the most valuable additions. Unstructured PEBI grids are lower priority.
- **Nanoparticles (EOR, fines, wettability).** Needs:
  - an advection–dispersion transport component (a thermal component with no gas/oil K-value is already transported);
  - Langmuir-type retention with a critical-velocity release term;
  - pore-throat entrapment;
  - K = K0 (φ/φ0)³ from the retained volume (`PERMLAW` Kozeny–Carman exists);
  - wettability alteration through k_r interpolation by adsorbed concentration — **planned** (same structure as the
    polymer adsorption of the black-oil core).
- **Advanced CO2.**
  - CPG (CO2-plume geothermal) = CO2STORE + THERMAL: this works in the thermal core today, with PR CO2, solubility and
    hysteresis.
  - CarbFix-type mineralisation (> 95 % within 2 years): use an aqueous-CO2 reaction with a `FORMULA` carbonate
    product.
  - Palandri–Kharaka rate laws (log k25, E_a) can be entered as reactions (calcite −5.81 / 23.5 kJ, albite −12.56 /
    69.8 kJ, anorthite −9.12 / 17.8 kJ, kaolinite −13.18 / 22.2 kJ).

### 6.6 C12: Langmuir sorption, embedded discrete fractures, the near-critical fold, Newton fixes

**Competitive Langmuir sorption (`SORB`, `thermal::addSorption`).** Each sorbing component gets a *site* solid and a
*sorbed* solid. A reversible kinetic reaction connects them:

- r = A·[c·[S] − [Ads]/b(T)], with b(T) = b_ref·exp(ΔH/R·(1/T − 1/T_ref)) and A = b_ref/(2τ), so τ is the relaxation time
  at c = 1/b.
- c is the gas partial pressure (GAS) or the aqueous concentration (AQ); AUTO picks GAS when the component is a
  non-condensable gas.
- Components that share a site name compete for the same sites (extended Langmuir).
- An optional sorbed-phase density occupies pore volume. EQ/NOEQ selects whether the initial loading is in equilibrium
  with the initial gas.
- Because sorption is a reaction, it reuses the existing mass and energy balance, Jacobian, inventory (the sorbed state
  counts as trapped) and element balance unchanged.
- Deck line: `SORB comp capacity b tRef dH tau [AUTO|GAS|AQ] [adsRho] [site] [EQ|NOEQ]`.
- Example `ecbm`: CO2 displaces CH4 on coal; sorbed CH4 drops from 1.23e5 to 4.17e4 kmol.

**EDFM (`FRACTURE`, `FRACSTRESS`, `thermal::addFractures`).** Planar vertical fractures are cut into segments by the host
cells of a Cartesian box. Each segment becomes an extra cell with:

- pore volume w·L·H (aperture w);
- a matrix–fracture connection T_mf = 2·k_m·A/⟨d⟩;
- fracture–fracture connections along the trace (cubic law k_f = w²/12) and at intersections (star–delta);
- well–fracture connections WI_f = 2π·k_f·w₀ / max(ln(r_e/r_w), 1), with r_e = 0.14·√(L² + H²), for every perforated
  host cell (`connectWellsToFractures`).

Fracture mechanics is explicit, updated once per accepted step:

- **Barton–Bandis closure:** w = w_open − V_m·σ′/(K_n0·V_m + σ′), with effective normal stress
  σ′ = σ_n + E·α/(1 − ν)·ΔT − p. Cooling unloads the fracture.
- **Coulomb hydro-shear:** when τ₀ > μ·σ′, dilation (τ₀ − μσ′)/k_s·tan ψ is added. It is irreversible.
- Stresses come from `FRACSTRESS` (S_Hmax / S_hmin gradients, azimuth, E, ν, α_T), resolved on each fracture plane.
- The aperture change is capped at **25 % per step**. Without the cap, the explicit update on a fast cold front made Newton
  fail. Transmissibilities rescale with w³ (cubic law) and pore volumes with w.

Results and UI:
- Results maps FRACAP / FRACP / FRACT per fracture cell.
- An editor in the Thermal → Fractures section.
- Example `egs_edfm`: 21 × 11 × 3 box at 200 °C, three fractures; 145 steps, 9.8 s.

**Fold construction (near-critical CO2–water).** This closes the C11 known limitation (§6.4). Ideal mixing of real-gas
components leaves a narrow fold of the cell energy E(T) near water's critical point, where dE/dT < 0 between a local
maximum E1 at T1 and a local minimum E2 at T2. On the fold:

- For U in [E2, E1], three temperatures solve E(T) = U, and the root that was found depended on the starting value.
- The fold is now treated like a first-order transition. The cell is a two-state mixture of the outer-branch roots P
  (E = E2, T < T1) and Q (E = E1, T > T2), blended with θ = (U − E2)/(E1 − E2).
- T(U), the phase split and the volumes are therefore single valued and continuous.
- Fold detection uses a *fixed* 0.5 K sample grid on 605–690 K, so it does not flicker under Jacobian perturbations. It is
  active only for 0.98–1.35 p_c and 615–680 K.
- An earlier version (v1) searched a window around the current root. It flickered and failed at day 235.
- `RESSIM_NOFOLD=1` disables the construction for comparison.
- Result: the ucg_co2 150 sm3/d / 120 bar variant, which failed at day 244.3 in C11, now runs through day 270 (430 steps, 2606 Newton iterations; checked to that time).

**Newton fixes found by the C12 analytic cases.**

- **Iteration-0 criterion.** The loose absolute tolerances (mass 1e-6, energy 1e-3) let slowly diffusing cells skip
  whole steps. With small Δt, the per-step change of a far cell stayed below the tolerance and was never applied, so
  pressure and heat fronts lagged: a bias that depended on time-step size and was found with Lauwerier and Theis.
  Iteration 0 now accepts only a state converged to a 1e4× tighter level.
- **Temperature chop.** The Appleyard-type local chop also limits |ΔT| ≤ 80 K per update, with the same property-only
  flash check. This targets cells with almost no rock heat capacity, such as EDFM fracture cells. The full regression
  showed no side effects.

**Injector start-up (observed in SPE4, not changed).** Injectivity follows ECLIPSE / CMG: total block mobility × injected
molar density. A cold, viscous well block therefore holds a steam injector at its BHP limit until the block heats up,
about 0.5 d in SPE4 1A. Chevron's convention in SPE4 (k_r = 1, injected-fluid viscosity) injects at full rate at once.

---

### 6.7 C14: MINC, hydraulic fracture propagation, restart

**MINC (multiple interacting continua, Pruess & Narasimhan 1985).** With `DUALPORO` and `MINC n [dims] [f_1..f_n]`
every matrix block is split into n nested shells (shell 1 next to the fracture). For fracture sets in `dims` directions
the matrix block is a slab / square / cube of size L; the shell boundaries are the iso-distance surfaces at the volume
fractions V_j (default V_j ∝ r^j, r = 1 + 3/n, outer shell thinnest because the early transient is steepest there).
Transmissibility between consecutive shells (and fracture → shell 1) is T = k_m A_j / d_j with the interface area A_j of
the nested surface and the distance d_j between the shell centroids (fracture → shell 1: half the outer-shell thickness).
The connections are ordinary 'M' grid connections, so every black-oil feature (wells, PVT, hysteresis, restart, adjoint
excluded) applies; the Warren–Root σ term is not used when shells exist. Verification `-minc`: 1-D slab drainage after
a fracture pressure drop, recovery F(t) vs the analytic series F = 1 − Σ 8/((2m+1)²π²) exp(−(2m+1)²π² t_D):
max |F − F_exact| 0.013 (10 shells), 0.026 (4), 0.373 (Warren–Root, pseudo-steady).

**Hydraulic fracture propagation (`hydrofrac.cpp`).** Planar fracture, two symmetric wings, Newtonian fluid:

- mass: ∂A/∂t + ∂q/∂x + q_L = 0, Carter leak-off q_L = 2 h C_L / √(t − τ(x)) (τ = exposure time);
- PKN (height H, L ≫ H): A = (π/4) H w, q = −(π H w³ / 64 μ) ∂p/∂x, p_net = E′ w / (2H), E′ = E/(1 − ν²);
- KGD (plane strain): A = w per unit height, q = −(w³ / 12 μ) ∂p/∂x, p_i = Σ_j C_ij w_j with constant
  displacement-discontinuity elements C_ij = E′/(4π) [1/(x_i − x_j + a) − 1/(x_i − x_j − a)] + mirror wing;
- discretisation: fixed uniform grid, backward Euler, Newton (tridiagonal for PKN, dense LU for KGD), w ≥ 0 as a
  complementarity condition (a cell that runs dry leaks less than Carter; the realised leak-off = injected − stored, so
  the volume balance is exact). The zero-toughness front is the degenerate-diffusion front (w³ face average);
- toughness: the next cell may open only if K_I ≥ K_IC; PKN K_I = p_net √(πH/2); KGD displacement correlation on the
  tip element K_I = E′ w √(2π) / (8 √r) / 1.253 (the factor calibrates constant DD elements: a uniformly loaded crack
  gives 1.238–1.253 for 5–200 elements); with toughness at most one cell opens per step and the step shrinks when the
  front is step-limited;
- front position from the near-tip asymptote (w³ ∝ r PKN, w^{3/2} ∝ r KGD viscous, w² ∝ r LEFM) extrapolated from the
  two cells behind the partially filled tip cell.

Verification `-hf`: PKN vs Nordgren (1972) L = 0.68 (E′q³/(2μH⁴))^{1/5} t^{4/5}: −1.7 % (L), +1.7 % (w0); KGD vs
Adachi–Detournay M-vertex L = 0.6152 (E′Q³/12μ)^{1/6} t^{2/3}: +5.5 % (100 cells) → +4.0 % (200), w0 −2.7 %;
KGD K-vertex L = 0.9324 (E′Qt/K′)^{2/3}: −0.4 %; PKN high leak-off limit L = q√t/(πHC_L): +2.8 % (efficiency 0.03).
EGS stimulation sizing (designer): one stage per plane, centred on the injector, aperture = retained fraction × mean
width; stages are independent (no stress shadowing - phase 2).

**Restart.** `saveState` / `loadState` serialise everything that evolves (stateIO); a restart replays the schedule
silently up to the restart step (well / group controls, events) and then loads the state, so the continuation is
bitwise identical to the uninterrupted run. Lessons: never re-store the "old" accumulation after loading (it is part of
the state); the derived arrays (T, p, S) of the thermal core belong to the state because they seed the next flash
temperature guess and the mechanics load (refreshing them from X changed SAGD-THM results by 1.4e-2).

### 6.8 C15: hydraulic fractures in black-oil / compositional runs, adsorbed gas, stress-dependent permeability, fault & fracture stability

**Hydraulic fractures outside thermal (`grid.cpp` HF block, `simulator.cpp`, `compsim.cpp`).** A fracture is a vertical
planar rectangle `-- RESSIM: HFRAC name x0 y0 x1 y1 zTop zBot width perm poro [gamma] [satnum]` embedded in any host grid
(Cartesian or corner-point) with the EDFM construction of §6.6: the trace is clipped by every host column, each piece
becomes an extra cell (pore volume φ_f w L h, permeability k_f or the cubic law w²/12 when perm = 0) appended after the
matrix / LGR / MINC cells. Connections: matrix–fracture T = 2 k_m A / ⟨d⟩ (⟨d⟩ = average normal distance of the host to
the plane), fracture–fracture along the trace and vertically (harmonic k_f w h / ½L), intersections of different
fractures inside the same host (Moinfar et al. 2014), well–fracture for every completed host cell the fracture
crosses (the perforation cluster). The cells are ordinary 'M' cells, so black-oil, dual porosity (fractures connect to
the fracture continuum), compositional, restart and result files apply. ECLIPSE / OPM export keeps the
HFRAC lines as comments (ResSim-only); `ecl_write` omits the extra cells with a warning.

**Multistage designer (`hydrofrac.cpp` StageLayout).** Stages along a horizontal well (spacing, toe offset, strike
⟂ S_h by default), each half-length / height from the C14 PKN / KGD solver or entered, propped width
w_p = C_p / ((1 − φ_p) ρ_p) from the proppant concentration; the stages are written as HFRAC lines.

**Adsorbed gas (`-- RESSIM: LANGMUIR V_L P_L ρ_grain`).** Gas content V = V_L p / (P_L + p) (sm³ per tonne of rock),
stored mass per cell (1 − φ) V_b ρ_grain V(p) ρ_g,sc, added to the gas accumulation (fully implicit, derivative
V_L P_L / (P_L + p)²). Reported separately in place (free vs adsorbed).

**Stress-dependent permeability (`-- RESSIM: PERMSTRESS γ` and per-fracture γ).** k = k_0 exp(−γ Δσ′) with
Δσ′ = −Δp (constant total stress): T_ab ← T_ab · exp(−(γ_a (p_a,i − p_a) + γ_b (p_b,i − p_b)) / 2), differentiated
in the Jacobian; closes propped fractures and stimulated matrix during depletion (Raghavan & Chin 2004).

**Fault & fracture stability (`faultstab.cpp`, Tools → Geomechanics).** Stress tensor from S_v, S_H (azimuth), S_h and
p (gradients or the run), resolved on every plane (faults, HFRAC planes, user planes):
σ_n = n·σn, τ = |σn − σ_n n|, σ_n′ = σ_n − α p; slip tendency T_s = τ/σ_n′ (Morris 1996), dilation tendency
T_d = (σ_1 − σ_n)/(σ_1 − σ_3) (Ferrill 1999), CFF = τ − μσ_n′ − c, critical pressure p_c = (σ_n − (τ − c)/μ)/α and the
margin Δp_c. Stress path of a laterally extensive reservoir dS_h = dS_H = α(1 − 2ν)/(1 − ν) dp (Segall 1989, Engelder &
Fischer 1994; can be switched off for localised injection). Views: Mohr circles (initial / current) with the Coulomb
envelope, evolution of T_s and Δp_c along the pressure history of the run, lower-hemisphere stereonet of T_s with the
planes' poles, per-cell CFF / T_s map from the run pressures.

**Fracture movement.** Aperture: Barton–Bandis closure Δv_n = σ_n′ V_m / (K_n0 V_m + σ_n′) plus shear dilation
min(δ tanψ, w_dil,max) → w = w_res + V_m − Δv_n + w_dil, permeability multiplier (w/w_0)³. Slip: quasi-dynamic spring
slider with rate-and-state friction (Dieterich–Ruina aging law) and radiation damping η = G/(2v_s) (Rice 1993):
τ_0 + Δτ(t) − kδ − ηV = σ_n′(t)[f_0 + a ln(V/V_0) + b ln(V_0θ/D_c)], dθ/dt = 1 − Vθ/D_c, k = 7πG/(16r);
unstable (stick-slip) when k < k_c = σ_n′(b − a)/D_c; events (V > 1e-3 m/s) report slip, M_0 = Gπr²δ,
M_w = (2/3)(log10 M_0 − 9.1). Integrated with Dormand–Prince 5(4), error control on δ/D_c and ln θ.

Verification: `-stab` (optimal plane 60.5° for μ 0.6: T_s = 0.494, Δp_c 2.45 MPa vs analytic; weakening patch
k = 1.37e8 < k_c = 6.94e8 → one event at 132.9 d, slip 0.010 m, M_w 2.87; strengthening (a > b): aseismic creep);
`-hfbo` (EDFM fracture vs Gringarten infinite-conductivity solution max error 3.7 %; corner-point = Cartesian to
5e-9; 6 stages raise oil 2137 → 3688 m³; PERMSTRESS γ 0.004 /bar 1002 → 752 m³). Example `SHALE_HF.DATA`
(10-stage well, 5673 matrix + 1380 fracture cells, free 1.08e8 + adsorbed 6.31e7 sm³ gas in place).

### 6.9 C16: plasticity and damage, HPC (mixed precision, GPU, MPI), thermal CO2 fixes

**Elasto-plastic geomechanics** (`geomech.cpp`). Drucker–Prager yield f = √J2 + α I1 − k matched to Mohr–Coulomb
(`dpMatch` 0 outer cone α = 2 sin φ/(√3(3 − sin φ)), 1 inner cone (3 + sin φ), 2 plane strain α = tan φ/√(9 + 12 tan² φ),
k = 3c/√(9 + 12 tan² φ)). Non-associated flow with dilation angle ψ (α_ψ from ψ). Implicit return mapping per Gauss point:
elastic trial → Δλ = f / (G + K α α_ψ + H(…)); apex return when the trial lies beyond the cone; the residual cone
(c → c_res) is returned to when softening crosses it. Cohesion c(ε̄_p) = max(c_res, c0 + H ε̄_p), with H per unit equivalent plastic
strain. Global equilibrium by initial-stiffness iterations with **Aitken Δ² relaxation** (plain iteration saw-toothed under softening),
committed / trial plastic state (commit after the converged flow step, so a cut time step restores the state).
Damage permeability k/k0 = min(kmax, exp(c_k ε_p,vol)) lagged per step into the transmissibilities of the thermal core.
Boundary conditions per face: fixed, traction, roller and **rigid plate** (one common normal displacement, load- or displacement-
controlled; plate force reported). Deck: `-- RESSIM: GEOPLAST model match ψ H c_res c_k kmax [iter tol]`,
`-- RESSIM: MECHBC face type stress [disp dispControl]`.

**Mixed precision.** `BILU0f` stores the block ILU(0) factors in float (50 % memory, same pattern); the Krylov vectors and
the outer solve stay double. A failed float factorisation (tiny pivot) falls back to double automatically.

**GPU solver.** Own D3D11 device (the render context is not thread safe) + mutex; compute kernels SpMV (BSR), block-Jacobi
apply, axpbypcz and two-dot reductions; a float BiCGSTAB inner solve with a block-Jacobi preconditioner is refined in double
on the CPU (`refineWithDevice`). Triangular BILU0 solves are sequential, therefore the device uses block-Jacobi. Fallback to the CPU
path on any device error.

**MPI.** `mpi_rt` loads msmpi.dll / libmpi.so.12 at run time (one ABI shim, no link dependency). `ResSimCLI --mpi-batch jobdir`
distributes case files round-robin over the ranks (job.txt, progress.txt, cancel.txt, done.txt, case_XXXX.rsres). Studies use it
for case-parallel runs; one simulation is not domain-decomposed.

**Thermal CO2 fixes (found by the Pruess EGS verification).** (1) With `gasPR` the injected enthalpy and molar density of the gas
stream now include the Peng–Robinson residual terms (before: ideal-gas enthalpy — cold dense CO2 entered with ≈ 12 kJ/mol too much
energy, i.e. effectively warm). (2) LBC over-predicts dense CO2 viscosity (2.7e-4 vs 1.35e-4 Pa s at 20 °C / 500 bar); CO2-rich gas
now uses μ = μ_LBC (μ_Fenghour / μ_LBC,CO2)^x_CO2 at the PR molar density (exact for pure CO2, reverts to LBC for other gases).

**Produced steam.** `WellReport::qsteam` — water leaving a producer in the vapour phase (CWE): per connection the water inflow is
split by the aqueous / vapour molar mobilities, q_steam = q_water · Σ_j w_j a_G,j / Σ_j w_j (a_W,j + a_G,j), w_j = WI_j (p_j − p_wj).

### 6.10 C17: implicit steam limit, MINC everywhere, conduction geometry, dispersion, multi-cluster fracturing

**Implicit steam-rate constraint** (`thermsim.cpp`, `STEAMTRAP 'well' subcool steamMax`). C16 applied the producer steam limit
lagged (rate cut after the step) and ≈ 33 % of the SPE4 reports exceeded it. Now the limit is a well control inside Newton: per
connection the steam rate measure A_s,j = WI_j · (M_w/ρ_w,CWE) · λ_v,w,j (vapour-phase water molar mobility, with its derivatives)
gives the BHP that produces exactly q_s,max, b_S = (Σ A_s,j p_j − q_s,max)/Σ A_s,j. When b_S is above the BHP of the primary control
(liquid rate or BHP limit), the well switches to "rate = q_s,max on the steam measure" with the full Jacobian of A_s; a 0.02 bar
hysteresis band and a two-switch rule stop flip-flopping. SPE4 2A: steam ≤ 1.02 × limit in 100 % of the post-breakthrough reports.

**MINC in the thermal and compositional engines.** Both deck builders used the unexpanded per-cell input arrays (`c.grid`), so MINC
shells and LGR cells (global index ≥ grid size) read defaults. They now map every simulator cell to its host input cell
exactly as `Simulator::setup` does: the deck builders now work on a copy whose grid arrays are expanded
(`expandLgrArrays` with the LGR-local PORO/PERM/…, `expandMincArrays`, `applyDpnum`, `prepareHydraulicFractures`), and the arrays
that are not expanded (enumerated initial state, HEATCR/THCONR/TEMPI/XMF) map each cell to its host (`hostQ`: MINC shell →
`mincHost`, LGR cell → `lgrHost`, fracture cell → its host). Review finds: in the compositional engine LGR cells had porosity 0.2,
no diffusion and SATNUM 1, and in both engines the Peaceman index of a completion inside an LGR read PERMX past the end of the
array; DPNUM = 0 folding was ignored by the compositional diffusion. Check `-dlgr` (3 × 3 LGR around a producer, compositional
and thermal decks). Matrix–fracture / shell diffusion now uses the shape-factor geometry with the
matrix porosity, T_d = (A/d)_geo · φ_m (before: no diffusion across 'M' connections).

**Conduction / diffusion geometry** (`grid.cpp`, `Connection::geo`). The thermal conduction and the compositional diffusion
transmissibilities were derived as T / k_harmonic, which is wrong when the two half-distances differ and the permeabilities
differ (layered rock). Every geometric connection now carries the pure half-cell geometry
geo = mult / (1/g_a + 1/g_b), g = A / d_half (the same half-cell vectors as the transmissibility, times MULTX/MULTFLT/MULTREGT,
without NTG and permeability); conduction G = geo · λ, diffusion T_d = geo · NTG_h · φ_h. MINC shells carry the
proximity-function geometry A_j / Δy_j. LGR / NNC / PINCH / HF connections fall back to T / k_h. Check: two layers of 1 m and 10 m give
G = 100/5.5 = 18.18 m exactly (before 2-cell permeability-weighted value). Dual-porosity 'M' conduction had been halved — fixed.

**Mechanical dispersion** (`compsim.cpp`, OPM `DISPERC` GRID array, per-cell dispersivity α). For each connection and each phase
present on both sides, with the TPFA molar phase flux F_α = T λ_α ΔΦ_α:

  J_i,α = (α_f / d) · |F_α| · (x_i,a − x_i,b),   α_f = 2 α_a α_b/(α_a + α_b),  d = centre distance,

i.e. the isotropic tensor α_L = α_T: α |v| ξ A/d Δx with the face-normal Darcy velocity (|v| A ξ = |F|). Full Jacobian including
∂|F|/∂(p, N). Verification (`-cdisp`): 1D miscible C1/C3 displacement, the front variance grows by 2 α x_f on top of the upwind
numerical dispersion — α recovered 0.96 (α = 1 m) and 1.86 (α = 2 m; the inlet boundary truncates the larger spreading);
DISPERC export → import round trip and α/d on the connections exact. UI: Trapping/diffusion tab, "Mechanical dispersion".

**Multi-cluster stages: stress shadowing + limited entry** (`hydrofrac.cpp`, `hf::simulateClusters`). Sequential stages in pure PKN
only raise the BHP without changing geometry, so the shadow is modelled where it matters — simultaneous clusters sharing the stage
rate. Cluster k is a PKN fracture (zero toughness, Carter leak-off) fed by q_k from one wellbore pressure:

  p_w = S_hmin + Δσ_k + Δσ_shadow,k + p_net,k(0) + K_perf q_k²,  Σ_k q_k = Q   (fully implicit with all fractures)

  Δσ_shadow,k = Σ_{j≠k} p̄_j [1 − ξ³/(1 + ξ²)^{3/2}],  ξ = |s_k − s_j| / c_j,  c_j = min(H, 2L_j)/2   (Sneddon, normal through the centre)

  K_perf = 8 ρ / (π² n² d⁴ C_d²)   (limited entry).

Checks (`-hf`): 1 cluster = Nordgren within 1.2 %; 4 clusters 2000 m apart = Nordgren(Q/4) within 0.1 %; 5 clusters at 10 m: outer
clusters 239.5 m / interior 95.5–109.5 m, rate shares 0.387/0.061/0.104 (shadow 1.0–2.1 MPa); limited entry (4 × 10 mm perforations)
evens the lengths to 163.5–169.5 m. Mass balance 1e-16. UI: Hydraulic fractures → designer → "Multi-cluster stage" (cluster table,
rate shares, half-length bars, "use the mean cluster half-length" for the stages).

**Figure export.** Log axes (decades with 2–9 minor ticks) in the vector PDF/SVG export and the y-label shrink-to-fit.

### 6.11 C18: compositional / thermal adjoint, in-run fracture growth, single-run MPI domain decomposition, GPU ILU

**Compositional and thermal discrete adjoint** (`stepadjoint.h/.cpp`, generic step-adjoint driver used by both engines). The black-oil adjoint of C9/C11 is
extended to the other two FIM engines with the same recipe: the forward run records every converged step (state, dt, lagged data),
the backward sweep solves Jᵀ λ_n = −(∂F/∂x_n)ᵀ − (∂R_{n+1}/∂x_n)ᵀ λ_{n+1} with the transposed step Jacobian, and the gradient is
dF/dθ = Σ_n λ_nᵀ ∂R_n/∂θ (+ the explicit ∂F/∂θ). The compositional residual is AD-exact (the flash-derivatives of the converged
EOS state enter through the implicit-function theorem); the thermal residual contains finite-difference Jacobian blocks (steam
tables), so its gradient carries a 1e-4 relative error. Checks (`-cadj`): compositional porosity gradient +347.1291 vs central FD
+347.1287 (7.5e-6); thermal 33 steps, F = 160.000, 1.1e-4 – 4.0e-4 relative error (tolerance 2e-3). Frozen in the adjoint: rock
hysteresis, group allocation, the lagged steam-trap switch (identical to the black-oil limitations).

**In-run hydraulic-fracture growth** (`simulator.cpp` `setupHfGrowth / growFractures`, deck `-- RESSIM: HFGROW name K_IC σ_min
[dσ/dz] [L0] [closedMult]`). The HFRAC trace is the maximum extent; every EDFM piece of the trace exists from the start (fixed
unknown layout, so the Jacobian pattern never changes). Pieces farther than L0 from the injection point (the pieces in the hosts of
the completed wells, else the trace centre) are *closed*: their fracture–fracture transmissibilities are scaled by k_host,n / k_f ·
closedMult, i.e. they conduct like the host rock; the host–fracture transmissibility (already a matrix value) is kept. After each
converged time step, for every fracture and wing, the outermost open column (all layers at that distance along the trace) is the
tip, with the highest fracture pressure p_tip of that column and its depth z_tip:

  p_net = p_tip − (σ_min + dσ/dz · (z_tip − z_mid)),  K_I = p_net √(π c),  c = min(L_wing, H/2)

(c = L: KGD / Griffith crack while the wing is shorter than half the height, c = H/2: PKN-type constant-height regime, which is
what `hydrofrac.h` uses for K_I = p_net √(πH/2)). When K_I ≥ K_IC the next column opens (one column per wing per step), its
transmissibilities are restored, the event (t, wing, L, K_I, p_net) is logged and stored, and the next step may not grow
(dt_next ≤ dt) while the fracture moves. The newly opened sliver fills from the fracture in the following implicit steps, which
gives the fluid-volume coupling (pressure drops after an opening, rises again under injection). State `RESSIM-BO-STATE-2` stores
the open flags and the event list, so restarts continue a growing fracture exactly; snapshots carry `HFOPEN` (1 = open piece).
Check (`-hfg`, 200 m trace, 25 m initial wings, 2 sm³/d water into 0.1 mD rock, K_IC 1 MPa√m, H = 10 m):

| case | result |
|---|---|
| K_IC unreachable (1e6) | wings stay 25 / 25 m, BHP 392.6 bar at 86 d |
| σ_min 320 bar | symmetric growth, 75 m at 6 d, full 100 m at 13 d, final BHP 342.2 bar (fully open trace 342.1) |
| σ_min 340 bar | 75 m at 30 d; injection pressure while propagating 335.4 – 344.2 bar, i.e. pinned at σ_min + K_IC/√(πH/2) = 342.5 bar ± one opening transient |
| restart in mid-growth (55 m wing) | identical wing lengths, events and BHP (0 difference) |

Scope: black-oil engine (the compositional engine keeps the full trace open, the checker says so); the growth path is a straight
vertical plane (no turning — direction comes from the input trace, e.g. along S_Hmax); adjoint gradients are refused for growing
fractures (the transmissibilities change with time). Opening is quantised by the EDFM piece length (= host cell size along the
trace); refine the grid along the trace for smoother growth curves.

**Single-run MPI domain decomposition** (`mpi_dd.h/.cpp`, `bsr.h` `setDistributedSolver`). The Newton linear solve of the black-oil
engine is distributed: rank 0 assembles the Jacobian (as before), the cell rows are split into contiguous blocks (one per rank), the
sparsity pattern is scattered once and the values per solve, and BiCGSTAB runs over all ranks with halo exchange of the block
vectors (Sendrecv in rotating-offset rounds; zero-length rounds are kept so partners never wait for a message that is not sent) and
Allreduce dot products. Preconditioner: **CPR** two-stage — stage 1 a global AMG V-cycle on the quasi-IMPES pressure system,
gathered to rank 0; stage 2 block-Jacobi BILU0 on each rank's block. Plain block-Jacobi BILU0 was too weak (SPE9LIKE 7307 linear
iterations vs serial CPR 849). With CPR: SPE9LIKE 900 d on 2 ranks 1176 linear iterations, 9.94 s (serial 9.27 s on this 2-core
box, OMP = 1), results identical (FOPT 2.64888e7). `-mpidd` checks: one rank = serial BILU0 iterations; 2 ranks CPR 11 vs BJ 29
iterations, error 1e-8; demo run FOPT identical. Run: `mpiexec -n N ResSimCLI X.DATA --mpi-dd -o X.rsres` (ranks > 0 serve the
solver, rank 0 simulates and writes the results), or in the app *Settings → Parallel → MPI → Single runs: domain decomposition*.
Limits: the assembly and the AMG pressure stage are on rank 0 (the solve is distributed, the assembly is not), black-oil only,
ranks ≤ physical cores (MPICH busy-polls; oversubscription is ~500× slower).

**MPI launcher fixes for Intel MPI on Windows.** `mpiexec` of the Intel Fortran runtime fails with `HYD_spawn … error 5 Access is
denied` when the libfabric provider directory is not on PATH, and when its proxies cannot break away from the job object of the
calling program. The launcher now prepends `…\opt\mpi\libfabric\bin` (oneAPI and Fortran-runtime layouts) to PATH, the app's
job object allows breakaway (`JOB_OBJECT_LIMIT_BREAKAWAY_OK`), `RESSIM_MPI_LIB` points ResSimCLI at the launcher's own MPI DLL
(impi.dll for Intel, msmpi.dll for MS-MPI), and `--mpi-info` / *Test MPI* print a diagnosis (MPI family, DLL, error-5 checklist:
Unblock / antivirus / Intel < 2021.6 → `-localonly` or MS-MPI). Not yet confirmed on the user's machine.

**GPU: multicolour BILU0 on the device, real timing** (`gpu_compute.cpp`, multicolour BILU0 in `bsr.h/.cpp`). Graph colouring of the block matrix
(greedy, ≤ 8–40 colours) makes ILU(0) parallel inside a colour: the factorisation is done on the CPU in colour order, the
triangular solves are one compute dispatch per colour on the GPU (forward then backward). The device timer separates *device busy*
(GPU timestamp queries) from wall time and counts readbacks. CPU reference (`-mp`, 120×100, 3×3 blocks): BILU0 20 iterations, MC-BILU0
33 (the price of the colouring), float storage + double refinement 5e-11. On wine/llvmpipe (no real GPU here) 108000 cells: CPU
BILU0 0.148 s, GPU Jacobi 0.652 s, GPU MC-BILU0 0.545 s (busy 0.382 s) — software rasteriser numbers; the *Benchmark* button
(Settings → HPC) and `ResSimStudio --gpu-bench` produce the same table on the user's GPU. A residual of exactly 0 at iteration 0
with b ≠ 0 (drivers that drop UAV writes) is now a device error, not convergence.

**Results file v2.** `.rsres` RunStats are size-prefixed (C17 files stored the raw struct; new GPU/DD counters would have shifted
everything after it); v1 files are read with the C17 layout. Engine restart states use the same rule (`RESSIM-BO-STATE-2`).

### 6.12 C20: seismicity-aware injection design, comp/thermal DD, HFGROW adjoint, P3D, embedded fracture contact (G1)

**Seismicity-aware injection designer (`seismicity.cpp`, research finding F1/F8).**
Dieterich (1994) seismicity rate on each fault patch $k$ under a Coulomb stress history $S_k(t)$:
$$R_k(t)=\frac{r\,e^{S_k(t)/A\bar\sigma}}{1+\frac{1}{t_a}\int_0^t e^{S_k(\tau)/A\bar\sigma}\,d\tau},\qquad N_k=\int_0^T R_k\,dt = r\,t_a\ln\Big(1+\tfrac1{t_a}\textstyle\int_0^T e^{S_k/A\bar\sigma}dt\Big).$$
With $S_k=\sum_w G_{kw} * q_w$ linear in the rates (poroelastic Green's functions), $\sum_k N_k$ is a log-sum-exp of affine functions → **convex** in the rate history (F1).
The designer computes $G$ from the black-oil engine (unit-rate steps per injector and period, $\Delta \mathrm{CFS}=\mu\,\beta\,\Delta p$ at the patch cells), minimises $\sum N_k$ subject to the target volume and rate bounds (projected Newton / SQP on the convex programme), then **re-simulates** the design.
Superposition alone was 17.8 % off (step-sequence dependence + amplification $e^{\Delta S/A\bar\sigma}$), so a sequential convex correction adds the simulated-minus-predicted stress of the verified schedule as an affine offset and re-optimises (F8): model error 17.8 % → 0.15 %, design −18.0 % events vs constant rate. The rate integral is checked against the ODE form to 3.3e-13.

**Compositional / thermal MPI domain decomposition.** The C18 distributed CPR-BiCGSTAB is now hooked into the compositional and thermal Newton loops (blocks up to `kMaxBlock` = 16; larger blocks fall back to serial with a log line). `ResSimCLI X.DATA --mpi-dd` under mpiexec: COMP and CO2STORE identical to serial, STEAMFLOOD +0.025 % in cumulative quantities at 2 ranks.

**HFGROW adjoint.** The growth state `hfOpen_` is captured per step (`AdjRec::hfo`) and restored with the fracture transmissibilities during the backward sweep, i.e. the discrete adjoint of the *frozen topology* (growth events are not differentiated: they are discrete). Check (`-hfadj`): $dF/dk$ adjoint 429.371 vs central FD 429.406 (8.1e-5), $dF/d\phi$ 1.8e-3; the FD step must stay inside one growth branch (h = 2e-4).

**P3D height growth (`hydrofrac.cpp`).** Symmetric three-layer section (pay height $h$, stress contrast $\Delta\sigma$ in the bounding layers). Each cross-section grows in height $H\ge h$ until $K_I=K_{IC}$ at both vertical tips:
$$p_{net}=\frac{K_{IC}}{\sqrt{\pi H/2}}+\frac{2\Delta\sigma}{\pi}\arccos\frac{h}{H},\qquad A(p)=\frac{\pi}{2E'}\Big[K_{IC}\sqrt{\tfrac{2}{\pi}}H^{3/2}+\frac{2\Delta\sigma}{\pi}h\sqrt{H^2-h^2}\Big],$$
contained sections ($p_{net}<K_{IC}/\sqrt{\pi h/2}$) use $p=2E'A/(\pi h^2)$. Lateral flow is PKN-type in the area form $\partial_t A+\partial_x q+2H_\ell c_L/\sqrt{t-\tau}=0$ with the section conductance $q=-\frac{A^3}{\pi^2H^2\mu}\,\partial_x p$ (equal to the PKN $\frac{\pi H w_0^3}{64\mu}$ for an elliptic section, $A=\pi H w_0/4$) and analytic $dA/dp$ (a numerical slope across the containment kink stalls the solver). Height cap and a runaway warning when $p_{net}\to\Delta\sigma$. Checks (`-p3d`): section inverse 5e-15; $\Delta\sigma$ = 1 GPa reproduces PKN (1.2e-5); $\Delta\sigma$ 10 / 6 MPa give L 1020 / 832 m with H 44.6 / 87.7 m; no-leak-off volume balance exact.

**Embedded fracture contact, G1 (`efem.cpp`).** 2-D plane strain, bilinear quads, fractures cut the mesh arbitrarily. Displacement
$$\mathbf u=\sum_a N_a\mathbf u_a+\sum_{a\in\mathcal E}N_a\,(H(\mathbf x)-H(\mathbf x_a))\,\mathbf b_a,\qquad [\![\mathbf u]\!]_\Gamma=\sum_{a\in\mathcal E}N_a\mathbf b_a ,$$
$H=1$ on the $+\mathbf n$ side; $\mathcal E$ = nodes of cut elements whose support does not contain a tip (zero jump on the tip edge). Cut elements are split along the fracture line, triangulated and integrated with a degree-2 rule. Local jump $\mathbf w=R^T[\![\mathbf u]\!]$ ($R=[\mathbf n\ \mathbf m]$); penalty contact $t_n=K w_n$ for $w_n\le0$ ($K$ = penalty·E/h), Coulomb return mapping $|t_t|\le c+\mu|t_n|$ with plastic-slip history per Gauss point; fluid pressure $p$ acts as $-p\,\mathbf e_n$. Interface terms $\int_\Gamma N_a R(\mathbf t-p\mathbf e_n)\,d\Gamma$ and $\int_\Gamma N_aN_bRTR^T d\Gamma$; full Newton on the banded system. Convergence is tested on the residual of the return-mapped tractions, and $w_n=0$ counts as closed so the first iterate uses the stick tangent.
The classic element-constant (SKON) EFEM is kept as `jumpModel = 1` for comparison: it over-predicts slip by 20–25 % on meshes oblique to the fracture and does not converge with refinement (F9, confirms ISRM 2015 / Cusini 2021 / GGGR 2024).
Checks (`-efem`, a/h = 12, 100×100): Sneddon pressurised crack at 0° and 30°: volume −0.4 % / −1.7 %, max opening +1.8 % / +1.0 %, profile error 1.9 % / 1.1 % of $w_0$; inclined frictional crack (30°, μ 0.3, σ_yy −10 MPa) vs the analytic elliptical slip $b(s)=\frac{4(1-\nu^2)}{E}(\tau-\mu\sigma_n)\sqrt{a^2-s^2}$: slip integral −5.0 %, max slip −1.9 %; μ 0.8: locked (|slip| 4e-8 m vs 2e-3 m frictionless).
Next (G-plan): EFEM ↔ EDFM pressure coupling (done in C21, §6.13), rate-and-state on the embedded faults and phase-field growth (done in C23, §6.15–6.16); 3-D G1b on hexahedra and DDM hydraulic fracture deferred to C24+; adjoint injection design with contact moved to research (H7).

### 6.13 C21: parallel-performance fixes, sub-grid fracture tip (H2), flow ↔ contact coupling (G1-c)

**Why MPI DD / GPU gave no speed-up on a 4-core laptop (H3, `research/CYCLE3_REPORT.md`).** Measured on CH15 (54k cells, 2 cores):
1 thread 48.6 s, OpenMP 2 threads 31.4 s (1.55×), old DD 2 ranks 48.4 s (no gain), C21 hybrid DD 36.0 s (1.37×). Causes and fixes:
(i) blocking MPI receives busy-poll a core each, and libgomp's default active wait spins idle OpenMP threads, so ranks × threads oversubscribe the
cores and Task Manager shows ~100 % CPU with no progress. The serving ranks now *politely wait* (`MPI_Iprobe` + sleep with a wake tag; `RESSIM_MPI_SPIN=1`
restores spinning), and `--mpi-dd` with threads prints a hint to set `OMP_WAIT_POLICY=PASSIVE` (libgomp reads it before `main`). (ii) the default
thread count was the number of logical CPUs: hyper-threads add memory-bandwidth contention to bandwidth-bound sparse kernels. Default is now one thread
per **physical** core (`physicalCores()`; a VM that reports every vCPU as core 0 falls back to the logical count when that is ≥ 2). (iii) the CPR
preconditioner was rebuilt on every Newton iteration; it is now rebuilt lazily (pressure AMG kept while the iteration count stays low, never reused
after a pattern change). (iv) DD distributes only the linear solve and costs 40 % more iterations (block-Jacobi coupling across ranks), so on ONE machine
it can never beat shared-memory OpenMP; it pays off across machines. (v) the GPU path has no device AMG: the pressure stage goes back to the host every
iteration and the D3D11 ILU is copy-based, so a laptop iGPU is slower than the CPU. Level-scheduled BILU0 (`RESSIM_LEVEL_ILU=1`) is bitwise identical
to the serial sweep but slower (22.2 vs 19.3 s apply: rows of a level are scattered in memory); a level-ordered contiguous copy is the next step.

**Sub-grid hydraulic-fracture tip, `HFGROW … SUBGRID` (research hypothesis H2).** The default in-run growth opens whole piece columns once per step
(length and time quantised by the grid). With `SUBGRID` each wing has a continuous length $L$; the tip piece has open fraction
$\phi=\mathrm{clamp}\big((L-s_{in})/\Delta s,0,1\big)$. Because $k_f=w^2/12$ is many orders above the host permeability, any $\phi>0$ makes the piece
hydraulically part of the fracture, so continuity has to come from the **leak-off area**: fracture–fracture transmissibility full, host–fracture
transmissibility $\times\phi$ (a φ-linear fracture conductance acts as a step and re-quantised the growth). Inside every time step the tip lengths solve
$K_I(L,p(L))=K_{IC}$ (bracket expansion, then safeguarded regula falsi, each trial re-solving the step; $|r|\le1\,\%$ or a 0.05 m bracket), with
$K_I=p_{net}\sqrt{\pi\min(L,H/2)}$ and $p_{net}$ the open-length-weighted mean fracture pressure minus $\sigma_{min}(z)$. The tip state is saved in
restarts (state v3). Results (HF_dx decks, 200 m trace, time for the + wing to reach 60 / 90 m):

| scheme | dx 10 m | dx 5 m | dx 2.5 m | cost at dx 10 |
|---|---|---|---|---|
| discrete (default) | 22.0 / 56.0 d | 24.0 / 56.0 d | 25.0 / 58.5 d | 75 steps, 1.6 s |
| SUBGRID | 29.5 / 70.2 d | 28.0 / 73.6 d | 26.7 / 60.9 d | 864 steps, 38 s |

The two schemes converge to the same limit **from opposite sides** (the discrete scheme opens a whole column too early; the continuous tip stalls each
time it enters a fresh host cell, where the leak-off per area is ~30× that of the pre-pressurised cells behind it). H2 is therefore **refuted** as
stated: the dominant quantisation error is host-cell leak-off, not the tip position. Use: SUBGRID is an opt-in research option (smooth $L(t)$, $K_I$ within
4 %); running both schemes brackets the grid-converged answer. A sub-grid leak-off model (Carter-type fresh-rock correction of the tip host cell) is the
follow-up.

**Flow ↔ contact-mechanics coupling, G1-c (`geocouple.cpp`).** Staggered and lagged by one report step: after each step the fracture-piece pressures of
one layer give each fracture a pressure profile $p(s)$ (EFEM `Fracture::profS/profP`; $s$ = distance from the midpoint), the EFEM contact problem of §6.12 is
solved over the grid's x–y extent with the slip history carried over, and each piece gets the hydraulic aperture
$$a=\begin{cases}a_{res}+V_m+w_n & w_n>0\ \text{(open)}\\[2pt] a_{res}+V_m-\dfrac{V_m\,\sigma_n'}{K_{n0}V_m+\sigma_n'} & \text{closed, } \sigma_n'=-t_n\end{cases}\;+\;\min\!\big(|w_t|\tan\psi,\ a_{dil}\big)$$
(Barton–Bandis hyperbolic closure, shear dilation), and the conductance multiplier $m=(a/a_{ref})^3$ (cubic law, $a_{ref}$ = HFRAC width), clamped to
$[10^{-6},10^6]$ and applied to fracture–fracture transmissibilities as the harmonic mean $2m_am_b/(m_a+m_b)$ (`Simulator::setHfConductance`).
The map-view (plane-strain) opening is scaled by $\min(1,H/L_f)$, the PKN/KGD width ratio $(2\Delta pH/E')/(4\Delta p\,a/E')$; without it a 200 m,
10 m high fracture opened 0.7 m instead of ~9 mm. The two crack-end cut elements are excluded from the piece interpolation (penalty tractions at the
tip are unreliable; within ~2 elements of a tip $\sigma_n'$ is still under-predicted, i.e. the aperture is conservative-high). Checks (`-geoc`):
Barton–Bandis limits exact; closed crack $t_n=\sigma_{yy}+p$: −3.998 vs −4 MPa; a constant profile reproduces the uniform-pressure solution exactly;
a 37 → 31 MPa profile opens only the high-pressure end (0.62 mm vs 0); coupled injection with a closed fracture ($\sigma_n'\approx6$ MPa, 0.1 mm HFRAC):
BHP 372.1 vs 345.9 bar uncoupled; pressurised fracture ($\sigma_{yy}=-31.5$ MPa): opens near the injector (8.0 mm at the well, 2.1 mm at the tips; hand
PKN estimate 9.3 mm) and lowers the BHP 414.0 → 406.6 bar. Limitations: total-stress elasticity with pressure on the fracture faces only (no poroelastic
matrix stress), explicit one-step lag (no iteration), one 2-D layer drives all layers, black-oil only. Use: Tools > Geomechanics (fracture contact)
> *Coupled flow run*, or `ResSimCLI X.DATA --geocouple sxx,syy[,E_GPa,nu,friction,layer]`.

### 6.14 C22: optimal-schedule structure and closed-form design potential (research H1f, F12), exact plot export

**Model.** Dieterich seismicity in the linear regime ($t\ll t_a$): expected excess events $I=\int_0^{\infty}(e^{\Pi s(t)}-1)\,dt$, with the
Coulomb change normalised by its constant-rate value at shut-in, $\Pi=\Delta S_{max}/A\bar\sigma$. A first-order stress response is
$\theta\,\dot s=u-s$ (time in units of the injection window $T_i$), with $u=q/\bar q\in[0,u_{max}]$ and $\int_0^1u\,dt=1$. After shut-in,
$s=s(1)e^{-(t-1)/\theta}$, so the tail is $\theta\,\mathrm{Tail}(\Pi s(1))$ with $\mathrm{Tail}(P)=\int_0^\infty(e^{Pe^{-x}}-1)dx=\mathrm{Ei}(P)-\gamma-\ln P$.

**Structure (Pontryagin).** $H$ is linear in $u$, so the optimum is made of bang arcs and singular arcs. On a singular arc
$\partial H/\partial u\equiv0\Rightarrow\dot\lambda=0\Rightarrow e^{\Pi s}$ is constant, which gives a **constant-stress plateau** $s^*$ with $u=s^*$.
The end condition is $\lambda(1)=\theta\,\partial_s\mathrm{Tail}$, which is smaller than the plateau multiplier, so the optimum ends with a
**maximum-rate burst** of length $t_b$. The structure is therefore charging arc ($u_{max}$, until $s=s^*$) → plateau → burst. The three-arc
cost is evaluated with Simpson's rule on each arc plus the analytic tail, and minimised over $(s^*,t_b)$: $s^*$ follows from the volume
constraint, and $t_b$ by golden section on $[0,\min(0.5,3\theta)]$. It reproduces the discretised optimum (SLSQP, 310–500 variables) to ≤ 0.001 %
for Π 1–8, θ 0.01–0.05 and $u_{max}$ 3 and 10. For θ → 0,

$$G=1-\frac{I_{opt}}{I_{const}}=\theta\,\kappa_2+O(\theta^2),\qquad
\kappa_2=\frac{e^\Pi\,[\Pi\,m(u)+R(\Pi,u)-\mathrm{Ein}(\Pi)+B^*]}{e^\Pi-1},$$

with $m(u)=(u-1)\ln\frac{u}{u-1}$, $R=\int_0^{\ln\frac{u}{u-1}}\big(1-e^{\Pi(u(1-e^{-x})-1)}\big)dx$, $\mathrm{Ein}(\Pi)=\gamma+\ln\Pi+E_1(\Pi)$, and $B^*$ the
optimised burst term. Without $B^*$ this is the pre-registered (refuted) law.

**Validity.** The closed form is exact for θ ≤ 0.1, or for $u_{max}\le3$ and θ ≤ 0.4. Beyond that, a coasting arc ($u=0$) appears before the
burst, and the closed form becomes a lower bound (flag `GainEstimate::exact`; the value is clamped at ≥ 0 because constant rate is feasible).
For diffusive responses (radial $\ln t$, linear $\sqrt t$) the same structure appears in weakened form, and the gain is ¼–¾ of the
first-order value (F13).

**Engine.** `seis::scheduleGain(Π, θ, u_max)` implements the closed form. `designInjection` computes the following on the governing patch
(the one with the largest weighted constant-rate excess):
- Π = max ΔCFS / Aσ̄ of the constant schedule;
- θ = (time to 63 % of the shut-in stress) / injection window;
- `farField` when the stress keeps rising more than 5 % after shut-in (closed form not applicable).

`DesignSetup::lastRateFactor` bounds the last control interval at a multiple of the mean rate. This is a box bound, so the problem stays
convex; a penalty on the shut-in stress would make it a difference of convex functions. The burst raises $s(1)$ by 26–58 %, and the bound
trades a little count for a lower shut-in stress. Tests: `-c22`.

**Exact plot export.** ImPlot and ImGui draw calls are captured through hooks (series, geometry, text and colormap scale) into a `Snapshot`:
filled triangles, polylines, text with font size and clip, in screen pixels. The `Snapshot` is replayed on the vector canvases
(PDF base-14 / WinAnsi, SVG) or rasterised (WIC PNG at 1–4× pixel scale). The page is 0.5755 pt per pixel. In print mode, neutral greys are
mapped through the theme's print-colour map (dark background → white, text → black, mid-grey grid kept as a light grey); saturated series
colours and colormaps are kept.

### 6.15 C23 G3: rate-and-state friction on an embedded fault (`rsf.cpp`)

Quasi-dynamic earthquake-cycle model on one straight fault (2-D plane strain, mode II, elastic full space), n equal elements of size dx.

* Friction (regularised Dieterich–Ruina, Rice et al. 2001): τ_f = σ′ a asinh( V/(2V₀) · exp(ψ/a) ); aging law
  dψ/dt = (bV₀/D_c)(exp((f₀−ψ)/b) − V/V₀) − α σ̇′/σ′ (Linker–Dieterich normal-stress term); steady state f_ss = f₀ + (a−b) ln(V/V₀).
* Elasticity: τ_i = τ⁰_i + Σ_j K_ij(δ_j − V_pl t) + Δτ_i(t) with the constant-slip edge-dislocation kernel
  K_ij = G/(2π(1−ν)) [1/(x_i−ξ⁺_j) − 1/(x_i−ξ⁻_j)] (back-slip loading at V_pl).
* Quasi-dynamic balance τ_i − ηV_i = τ_f(V_i, ψ_i, σ′_i), η = G/(2c_s) (radiation damping, Rice 1993), solved per element by a bracketed
  monotone Newton iteration; (δ, ψ) integrated by Dormand–Prince RK4(5) with error control on ψ and slip (rtol).
* σ′ = σ_n − p(t), floored at 0.1 MPa. Pore pressure: Theis line source with superposed shut-in, a table from the G1-c coupled run
  (`fromCoupled`), or a callback. Initial σ_n, τ either uniform or resolved from the EFEM contact solve along one fracture (`fromEfem`;
  τ capped just below the steady-state strength).
* Events: max V > 1 mm/s opens, below closes; M₀ = G Σ Δδ dx · W with W = rupture length (square rupture), Mw = ⅔(log₁₀M₀ − 9.1).
* Resolution: nucleation length h* = 2 G′ D_c b / (π σ′ (b−a)²), G′ = G/(1−ν) (Rubin–Ampuero 2005, a/b > 0.5 branch), cohesive zone L_b = G′ D_c/(σ′b);
  the UI warns when L_b/dx < 3.

Verification (`-c23`): one element = spring-slider, stick–slip only for k < k_c = σ′(b−a)/D_c (33 events at k/k_c = 0.5, none at 2);
finite fault (2 km, VW patch 0.2–0.8, 128 elements): h* 306 m, L_b 53 m, periodic events T = 5.13 yr, slip = V_pl T within 10 %,
Mw 4.43, Δτ 3.3 MPa; rtol halving changes T by 0.025 %; steady sliding |V/V_pl − 1| = 6e-6; injection 50 m from the patch advances the
first event from 10.0 yr to 0.002 yr; a velocity-strengthening fault only creeps; `fromEfem` stresses match the contact solve (31.73 vs
31.74 MPa). Limits: scalar slip on one planar fault, no fault-fault interaction, 2-D (W = L assumption for Mw), no fault opening.

### 6.16 C23 G4: phase-field fracture (`phasefield.cpp`, experimental)

Variational fracture with AT2 regularisation and pressure work (Bourdin, Chukwudozie & Yoshioka 2012):
Π(u,v) = ∫(v²+k)ψ(ε) − p∫u·∇v + G_c∫[(1−v)²/(4ℓ) + ℓ|∇v|²]. Hybrid formulation (Ambati 2015): isotropic degradation in momentum,
Amor volumetric–deviatoric split ψ⁺ in the history field H = max ψ⁺ (Miehe 2010); staggered, both sub-problems linear:

* momentum ∫(v²+k) ε(δu):C:ε(u) = ∫p ∇v·δu; phase ∫(2H + G_c/(2ℓ)) v δv + 2ℓG_c ∇v·∇δv = ∫G_c/(2ℓ) δv + p u·∇δv, then v ← min(v, v_prev).
* Mesh: rectilinear, uniform h in a core box (3a × 0.5a), geometric grading 1.12 to a 20a plate with fixed boundary.
* **Q4 + Wilson incompatible modes**, statically condensed per element (history strains include the modes). Plain Q4 shear-locks on the
  elongated elements where graded rows meet the core: frozen-profile Sneddon ratio 0.97 → 0.62 on refinement (dead end, kept in notes).
* **Initial crack = three node rows (|y| ≤ h) with v = 0.** With one row, the elements next to the crack line keep g = v² ≈ 0.002–0.05 at
  their Gauss points and resist opening the more the finer the mesh (volume ∝ (ℓ/a)^0.6 — found with the frozen AT2 profile, `PF_FROZEN=1`).
* Linear systems: BiCGSTAB + smoothed-aggregation AMG (u ≈ 150–170 iterations at 15k dofs, v ≈ 18), ILU0 fallback; `PF_DEBUG=1` prints them.
* Growth detection: pressure control is unstable after p_c, so p_c is bracketed by the load step in which the half-length (v < 0.5 on the
  crack line) jumps by more than max(4h, 0.1a). Sneddon ratio V/(2πpa²/E′) is taken at the first load step (AT2 has no elastic threshold;
  V/p creeps up as the tips soften).

Verification (`-c23`, a = 0.1 m, E 10 GPa, ν 0.25, G_c 100 J/m², ℓ = 2h):

| h (m) | ℓ/a | nodes | Sneddon ratio | p_c / p_Griffith (step midpoint) | time (2 cores) |
|---|---|---|---|---|---|
| 0.01 | 0.2 | 7 345 | 1.801 | 0.875 | 45 s |
| 0.005 | 0.1 | 16 095 | 1.380 | 0.925 | 206 s |
| l → 0 (linear extrapolation) | 0 | | **0.960** | **0.975** | |

p_Griffith = √(E′G_c/(πa)) with the physical G_c: at fixed ℓ the measured p_c hardly depends on h (ℓ = 0.01: 1.757 vs 1.783 MPa for
h = 0.005 / 0.01), so the Bourdin discretisation factor G_c,eff = G_c(1 + 3h/(8ℓ)) is reported but not applied. Honest reading: both
quantities carry the known O(ℓ/a) error of diffuse models; at ℓ/a = 0.1 the volume is ~40 % high. Limits: 2-D, horizontal initial crack,
pressure control (growth after p_c not followed), no fluid flow / leak-off inside the crack, no contact for closing cracks.

### 6.17 C23: stage-kernel calibration (`seis::fitStageKernels`) and native publication export

**Kernel calibration** (research Cycle 5, F14). Stage-wise event counts N_s ~ Poisson(λ_t K_s) with one productivity λ_t per test
(profiled out in closed form: λ_t = ΣN/ΣK over the test) and kernels K ∈ {1, duration, V, E = ∫pq dt, V q^β, V(1+V_cum)^γ,
V q^β (1+V_cum)^γ} (V_cum = volume injected earlier in the same test); exponents on a profile-likelihood grid in [−2, 3]; AIC ranking;
95 % profile-likelihood interval of β (Δ log L = 1.92). (The Dieterich kernel of the research script, research/exp/h6_aspo.py, needs
the pressure history and is not part of the in-app fit.) Äspö HF1 + HF2 (12 stages): best V(1+V_cum)^0.45 (AIC 71.5), joint
β ∈ [−1.25, 0.75] → no resolvable rate effect; corr(log q, log(1+V_cum)) = 0.91 explains why V q^β alone looks significant.

**Publication export.** The figure model (`figure.h`) gained native layers — Heatmap (cell values + colour map + colour bar), Shape
(polylines, rectangles, quads, markers for fractures / wells / cell outlines), ColorBar and text — rendered in z-order by the same vector
back end as the curves (PDF / SVG / PNG). The UI plots draw through `pk::` helpers (`plotkit.h`) that record these layers while they
draw, so the journal-style re-draw contains everything the screen shows; *Exact* mode (C22) still copies the screen draw lists.

### 6.18 C24: fault mechanics and induced seismicity inside the main run (`faultmech.cpp`), in-run G1-c

**Purpose.** Up to C23 the fault-stability, seismicity, rate-and-state and contact tools were separate *Tools* panels that
re-used a finished run. C24 makes them part of the model: every `FAULTS` face with a `FAULTMECH` row is a mechanical object
that is updated after each report step of the normal Run, can feed back into the flow (two-way), and is drawn in the main
3D view and the Results tab *Faults & seismicity*.

**Deck (ResSim-only `-- RESSIM:` comment lines, so ECLIPSE / OPM ignore them):**
```
-- RESSIM: STRESSFIELD svGrad shRatio sHRatio azimuthDeg stressPath sv0
-- RESSIM: FAULTMECH name mu cohesion aSigma ta_yr r0 bValue mMin permFactor
-- RESSIM: FAULTRSF  name layer a b Dc_mm
-- RESSIM: GEOCOUPLE sxx syy E_GPa nu friction layer
```

**Coupling loop (sequential, one report step lag).**
1. *Setup:* the fault faces are taken from the built grid (`FAULTS` boxes → the cell face polygon from the corner points,
   unit normal n, area A, the active cells a / b on both sides and the connections crossing the face).
2. *t = 0 and after every report step:* face pressure p_f = ½(p_a + p_b) from the flow solution.
3. *Stress:* total stress in the grid frame, S = S_H e₁e₁ᵀ + S_h e₂e₂ᵀ + S_v e_ze_zᵀ with S_v = S_v0 + g_v z,
   S_h = r_h S_v + A_p Δp, S_H = r_H S_v + A_p Δp (A_p = stress path, poroelastic horizontal change per unit Δp, Δp = p_f − p_f(0)).
   Traction t = S n, σ_n = n·t, τ = |t − σ_n n|.
4. *Coulomb:* CFF = τ − μ(σ_n − p_f) − c; ΔCFS = CFF − CFF(0); slip tendency τ / σ′_n. For the vertical strike-slip
   test geometry ∂CFF/∂p = μ(1 − A_p) exactly (τ does not change when both horizontal stresses shift equally).
5. *Seismicity (Dieterich 1994):* per face, state γ with dγ/dt = (1 − γ Ṡ)/(Aσ), R = r / (γ Ṡ_r), Ṡ_r = Aσ / t_a.
   Within a report step Ṡ = Ṡ_r + ΔCFF/Δt is constant, so γ is advanced with the exact exponential solution
   γ(t+h) = (γ − 1/Ṡ) e^{−Ṡh/Aσ} + 1/Ṡ; the expected count ∫R dt uses Simpson on 16 sub-steps.
   Background r = r₀ · A (events ≥ M_min per km² per year). Verified to 1 % against the closed form
   N(t) = r k [t + τ ln((1 + (k−1)e^{−t/τ})/k)], k = Ṡ/Ṡ_r, τ = Aσ/Ṡ (`-c24`; the test gives 1.32734 vs 1.32734).
6. *Catalogue:* Poisson (normal for > 50) number of events per face and step, uniform time within the step and position on
   the face polygon, magnitudes from a truncated Gutenberg–Richter law, M ≤ M_cap with M_cap = full-fault rupture
   (M₀ = (16/7) Δσ (A/π)^{3/2}, Δσ = 3 MPa).
7. *Two-way coupling:* a face with CFF ≥ 0 gets the transmissibility factor F^x, x = clamp(CFF / (0.05 σ′_n), 0, 1),
   irreversible (max over time); `Simulator::setConnectionMultipliers` multiplies the original transmissibility of the
   crossing connections before the next report step (fault-valve behaviour; F = 1 is one-way monitoring).
8. *After the run (FAULTRSF):* the faces of one layer form a trace; the quasi-dynamic rate-and-state model of §6.15
   (aging law, radiation damping) is driven by the simulated σ_n(t), p(t), τ(t) per element (interpolated from the faces)
   and adds rupture events (kind 1). rsf.h assumes a square rupture; C24 scales M₀ by min(1, H/L) with H the fault height.
9. *G1-c (GEOCOUPLE):* `geocouple::coupleOnce` (§6.14) after every report step of the same run, so hydraulic-fracture
   apertures / conductance follow the contact solution during the main simulation.

**Test case (`examples/FAULT_INJECTION.DATA`, `-c24`).** 30×20×5 cells of 50×50×10 m, fault F1 (100 faces, MULTFLT 0.05)
150 m east of a water injector (800 sm³/d, 280 bar cap, shut in at 500 d). Stress: Sv 0.24 bar/m, Shmin/Sv 0.65,
SHmax/Sv 1.10 at 60°, μ 0.6 → initial CFF −4.8 bar at 1980 m (hand calculation reproduced to 1e-3 bar).
Results: max ΔCFS 14.6 bar (two-way) vs 11.2 bar (one-way), 43 % of F1 reactivated at 100 d, 100 % at 200 d;
two-way coupling cuts the pressure jump across F1 from 47.6 to 43.7 bar at its peak and to ≈ 5 bar by shut-in, and raises
the expected events from 352 to 479 (the conducting fault carries pressure to the far side). Run time ≈ 1 s.

**Limits (honest).** One report step lag; the stress is analytical (linear poroelastic stress path, no 3-D
elastic stress transfer between faces, no thermal stress); Dieterich population model with a constant Aσ per face; the
catalogue is a stochastic realisation of the expected count (no ETAS triggering); FAULTRSF is 2-D plane strain on one
layer; black-oil engine only (compositional / thermal runs ignore FAULTMECH). The G4 phase-field model stays a sub-scale
tool (it is a 2-D lab-scale solver; coupling it into field runs would need volume-controlled growth and a flow model
inside the crack, ROADMAP G4-next).

### 6.19 C25: finite-element geomechanics for every engine on Cartesian and corner-point grids — the non-XFEM path (`geoseq.cpp`, `geomech.cpp`)

**Why not XFEM.** XFEM (or EFEM-3D) needs a fracture surface cutting the mechanics mesh, enrichment of the cut elements, a
contact/return-mapping kernel on the enriched dofs and a 3-D Newton with AMG — not verifiable on our hardware at field
sizes (§6.12, ROADMAP G1b) and fragile on corner-point meshes (sliver cuts). The C25 path keeps the three things a reservoir
engineer actually needs — (i) the full stress state on the real (dipping, faulted) grid, (ii) fracture growth against that
stress, (iii) fault stresses from that stress — with components that are each verifiable:

1. **Mechanics mesh = the flow grid.** Cartesian tensor grids: the brick elements of `geomech::Solver::setup`. Corner-point
   grids (COORD/ZCORN, dipping TOPS, irregular DX): trilinear **isoparametric hexahedra** (`Solver::setupHex`), node
   (i,j,k) = mean of the corner points of the up-to-8 cells sharing it. At a fault pillar the node is the mean of both
   sides, i.e. the throw is smeared over one element: the mechanics mesh is conforming, the flow grid keeps its NNCs.
   Over/under/side-burden padding extrudes the edge nodes. 2×2×2 Gauss quadrature; B-matrix from J⁻¹∂N/∂ξ; the brick and the
   isoparametric element give identical results on a box (test: |Δu| 2e-17, |Δσ'| 4e-8 Pa) and the patch test (uniform Δp,
   uniaxial strain ε_v = αΔp/M, Δσ'_zz = −αΔp) is exact on a distorted mesh (rel. error 4e-9).
2. **Sequential coupling (`geoseq::Coupler`) for black oil and compositional/CO2.** After every report step:
   Δp per cell → FE solve (elastic or Drucker–Prager) → σ (full tensor, compression +), ε_v, u, failure index.
   Two-way (GEOMECH coupling 1, explicit, one report-step lag):
   PV ← PV_flow · m, m = (φ_FE/φ₀)·exp(−c_r(p−p₀)) (the FE porosity replaces the flow rock compressibility; the engines
   keep moles/masses of the old level, so the update is mass-conservative — test: oil balance 9e-7, comp moles unchanged
   to 1e-12); k-multiplier k/k₀ = exp(b ε_p)·exp(−c_σ(σ'_m − σ'_m0)) (damage × stress law, GEOMECH item 25 `permStress`)
   → matrix transmissibility × harmonic mean of the two cells. Black oil composes geo × fault multipliers in one
   `setConnectionMultipliers` call; compositional/thermal map by cell pairs (`setConnectionMultipliersByCells`).
   Thermal decks use the THM mechanics of the thermal engine (fixed-stress, §6.x) on the same hex mesh.
3. **Fracture propagation on the grid (HFGROW … FESTRESS).** EDFM fracture pieces along the trace open column by column
   when K_I = p_net·√(π·min(L, H/2)) ≥ K_IC (KGD → PKN). With FESTRESS the closure stress of a piece is the **smallest
   principal FE total stress at its host cell** and the tip column is driven by its *largest* net pressure
   max_k (p_k − σ_min,k) (identical to the old rule for a linear σ(z) parallel to the fluid gradient). Black oil grows
   inside the time step; the compositional engine grows between report steps (`geoseq::HfSequential`, same rule).
   Poro-elastic stressing is therefore included automatically: in the CO2 example (2 mD) the tip σ_min rises with the
   pressure (p_max ≈ 211 bar under a 230 bar cap), the fracture grows one column and arrests — the analytic gradient
   would have predicted unlimited growth. In GEOMECH_CP the stress path matters decisively: with Biot α = 1 and
   k₀ = 0.45 the FE σ_min rose by ≈ 0.62Δp and the best net pressure stayed at −34 bar (no growth); the example therefore
   uses α = 0.7 and k₀ = 0.38 (initial σ_min 320.9 bar, unchanged by α), giving 7 growth events one-way and 6 two-way.
4. **Fault stresses from the FE stress (STRESSFIELD … FE).** σ_face = σ_in-situ(analytic) + (σ_FE(t) − σ_FE(t₀)) projected
   on the face normal; captures stress arching, stressing of the far side and thermal stress. Thermal decks with the
   analytic source can add ΔS_h = ΔS_H = −Eα ΔT/(1−ν) (`TH E α ν`; test: 50 K cooling, E 30 GPa, α 8e-6, ν 0.25 →
   Δσ'_n −160 bar, ΔCFF +96 bar with μ 0.6, exact).

**Limitations (stated honestly).** Explicit one-step lag for BO/comp (stable for the report-step sizes used, but the pore
volume jumps once per report step; thermal is iterative fixed-stress; C26: iterative option for black oil, §6.20); fault throws are smeared in the mechanics (no
slip discontinuity — fault slip is represented by the Coulomb/Dieterich/RSF layer, not by the FE mesh); fracture paths
are prescribed traces (no curving), widths are EDFM widths, not FE apertures (C26: DDM apertures with DDMWIDTH, §6.20); growth is column-wise (sub-grid tip option
for BO only).

**Compositional robustness fix (C25).** A tiny EDFM cell (0.03 m³) at the CO2 solubility limit flips phase state every
Newton iteration (2-cycle). Two standard remedies were added: damping of the update when a period-2 residual is detected
(halved per cycle, ≥ 1/16) and an OPM-style relaxed acceptance when the cells above the strict tolerance hold < 1e-4 of the
pore volume, are within 1e-2 of their content and the global material balance is strict (< 1e-12).

**Examples / tests (`-c25`).** GEOMECH_CP.DATA (black oil, corner point with 8 m fault throw; E 10 GPa, ν 0.25, α 0.7, k₀ 0.38), CO2_GEOMECH.DATA (CO2STORE,
tight 2 mD, fracture growth + fault), EGS_FAULT.DATA (thermal doublet, THM + fault + RSF).

### 6.20 C26: boundary-element (DDM) mechanics of faults and hydraulic fractures, iterative FE coupling, phase field with volume control (`ddm.cpp`, `hf_ddm.cpp`, `faultmech.cpp`, `geoseq.cpp`, `phasefield.cpp`)

**Kernels.** Displacement-discontinuity method with triangular dislocation elements (TDE): the closed-form full-space and
half-space solutions of Nikkhoo & Walter (2015, GJI 201:1119), taken from the MIT-licensed cutde translation
(`engine/src/tde_kernels.inc`, rendered from the cutde templates for double precision, namespace `ressim::ddm::tde`).
Every quadrilateral fault face or fracture piece is split into two triangles. The influence of element j on element i is
the traction t_i = σ(x_i; Δu_j = e_k) · n_i. **Important pitfall (C26):** TDE stresses evaluated *in the element plane* (the
quad centroid lies on the diagonal of its own triangles; grid neighbours lie on the extended diagonals) are finite but
wrong (no NaN, so a NaN retry never fires). The collocation point is therefore evaluated at x_i ± δ n_i (δ = 10⁻³ element
size) and averaged, which gives the principal value. Checks: kernel vs finite-difference displacement gradient 1.0e-9 (full and
half space); free-surface traction of the half-space kernel ~1e-17; opening self-term K_ii = −1.2004 G/h for a square element
(theory 1.2008), shear self-term −1.0504 G/h (theory).

**Fault slip (FAULTSLIP).** `-- RESSIM: FAULTSLIP E_GPa nu dMu [slipRef_mm] [cutoff] [FULLSPACE|HALFSPACE]`. Quasi-static
frictional equilibrium on the FAULTS faces after every report step: the Coulomb excess of the C24 stress model (STRESSFIELD,
pore pressure of the run) is relaxed by slip, τ_i + Σ_j K_ij s_j ≤ μ_i σ'_n,i, solved by a projected Gauss–Seidel / active-set
iteration on the dense DDM matrix (faces ≤ 3000; neighbour cut-off radius in element sizes). Slip-weakening μ → μ − dMu
after slipRef produces cascades. Outputs per face: slip, Δτ from slip, seismic moment M₀ = G Σ A s, moment magnitude, and
the slip-corrected Dieterich rate. Verification (`-c26`): penny crack 448 elements, volume within 5.0 % and centre opening
within 1.7 % of Sneddon; 80 × 10 m rectangle centre opening 1.971 vs 1.875 mm; FAULT_INJECTION: 100 faces slip, max
3.40 mm, M₀ 8.63e11 N·m (= G Σ A s exactly), expected events 352 → 192 (slip relaxes the stress that the
Dieterich model would otherwise accumulate); slipped faces end at their strength within 0.5 bar; faces not yet slipped
gain up to 0.16 bar of shear stress (stress transfer).

**Hydraulic-fracture apertures (HFGROW … DDMWIDTH [E nu]).** The open pieces of an in-run fracture form a DDM mesh; the
net pressure p − σ_min of every piece is the opening load, and the opening w solves K w = −(p − σ_min) (closed pieces
w = 0, active set). The fracture conductance follows from w: multiplier (w_eff / w_in)³ (cubic law) or linear. HFGROW_INJECTION:
21/21 pieces open, max aperture 0.48 mm, conductance × 0.0042; the rate-controlled run changes only in the local
pressure (max 343.29 vs 343.22 bar) because the injected volume is fixed.

**Iterative two-way FE coupling (GEOMECH coupling = 2).** The explicit scheme (coupling 1) applies the pore-volume
multipliers m(p) of the FE state of report step n during step n+1. The iterative scheme repeats each report step from its
saved start state (`saveState` / `loadState`) until max|m_new − m_used| ≤ tol (GEOMECH maxIter, tol): old-time masses keep the
start multipliers and the new time level uses the latest m, so mass is conserved in every pass and the accepted step is a
fixed point of flow and mechanics. `loadState` now keeps the summary and snapshots up to the restored time when the
simulator has already passed it (a rewind); a fresh simulator (restart) is cleared as before. Measured (`-c26`,
GEOMECH_CP with E = 1.5 GPa, flow time steps ≤ 5 d in every run): with 150 d report steps the explicit scheme is 0.239 bar
from the 5 d iterative reference in average pressure, the iterative scheme 0.095 bar (2.3 passes per step, 2.3× the cost of
explicit). **Honest note:** with the report steps of the example decks (5–50 d) explicit and iterative differ by < 0.01 bar
— the explicit lag only matters for long report steps or very soft rock. Black oil (GUI and CLI); the compositional runner
keeps the explicit scheme and says so in the log.

**Phase field with volume control, any orientation and in-situ stress (C24-b, `phasefield.cpp`).** The crack faces carry the
traction (pI − S)·n, S the in-situ stress (compression positive): the scalar pressure in Bourdin's terms becomes the tensor
T = pI − S in the momentum load ∫∇v·T·δu and in the phase-field load u·T·∇δv. This yields the net pressure p − σ_n and the
shear traction τ on the faces while the far field stays unloaded (u is the perturbation of the in-situ state). Volume control:
for a fixed phase field the momentum problem is linear, u = p u₁ + u_S, so each staggered iteration sets p = (V_target −
V(u_S)) / V(u₁) exactly; growth after the peak pressure is stable and followed. The initial crack may have any angle (the
uniform core box is rotated; a three-node band |n| ≤ h(|cos θ| + |sin θ|)). Griffith reference with stress:
p_c = σ_n + √(E′G_c/(πa) − τ²). Solvers (C26): both systems are made symmetric (Dirichlet rows and columns), solved by PCG with
AMG (hierarchy built once per load step and refreshed in the staggered iterations; component-decoupled operator for the
aggregation) — 3.7× faster than BiCGSTAB with a fresh AMG per solve. Verification (`-c26`, a = 0.1 m, h = 0.01 m, l = 2h):
horizontal crack, no stress: peak pressure 0.85 × Griffith (the known AT2 size effect at l/a = 0.2), then the pressure falls
while the crack grows, p√(πa)/K_IC ≈ 1.05 at a ≈ 2a₀ (toughness regime); 30° crack under isotropic 5 MPa: net peak ratio
0.861 vs 0.855 horizontal (orientation and stress invariance); 45° crack under S_xx = 6, S_yy = 5 MPa (σ_n 5.5, τ 0.5 MPa):
the crack turns towards the maximum compressive stress (tip direction 45° → 16.9°); its net peak pressure is 0.64 × the straight-crack mixed-mode Griffith value (the kinked path is not the straight extension the formula assumes). **Limits:** 2-D plane strain, rectilinear
mesh (mild mesh bias for inclined cracks), no fluid flow inside the crack (uniform pressure), not coupled into the field run —
the phase field remains the near-wellbore / laboratory tool; in-run fracture paths stay the LEFM grid-piece growth (§6.11/§6.19). *(C28: in-run coupling by volume added, §6.21.)*

### 6.21 C28: phase field inside the field run, distributed assembly, contiguous ILU(0), comparison items

**Phase-field fracture growth in the run (`HFGROW … PHASEFIELD`, `simulator.cpp` PfLink + `phasefield.cpp`).**

*Coupling.* Each PHASEFIELD fracture owns one persistent volume-controlled phase-field model of its plane (§6.20). It runs in a
worker thread as a coroutine: `phasefield::run` asks for the next target volume through `Input::nextVolume` and reports every
converged load step through `Input::onStep`. Displacement, phase field and history H = max ψ⁺ stay in memory between report
steps. At the end of every report step (`pfReportStep`) the simulator computes the crack volume per unit height
  V = η · V_inj / H.
Here V_inj is the cumulative water injected by the wells completed in the host cells of the fracture's pieces (surface volume ≈
reservoir volume; stb → m³), H is the fracture height (zBot − zTop, else the height of the injection-column pieces) and η is
the fluid efficiency. The simulator then posts V and waits for the step. The phase-field half-length L_PF (radial v = 0.5
measure) opens every piece with |s| ≤ L_PF on both wings (growth only, HfEvent, `applyHfTrans`). For that fracture this
replaces the LEFM piece test of §6.11.

*Net-pressure formulation.* The in-plane closure stress σ is isotropic for a straight crack. By superposition, the problem
with far-field compression σ and face pressure p equals the unstressed problem with face pressure p − σ. The model therefore
runs at S = 0 and reports p = σ + p_net. This avoids letting the O(h) error of the compressed regularised crack (σ ≈ 100 p_net)
swamp the net pressure.

*Mesh.* a₀ = max(initLen, half-height of the injection piece); L_max = trace extent. h = min(a₀/8, L_max/60), floored at
L_max/400, with ℓ = 2h. Core box: 1.08 L_max × max(0.5 a₀, 6h). Growth stops at 0.93 of the core, which is reported as
"maximum extent".

*Material.* E′ = E/(1 − ν²); E comes from the PHASEFIELD item, else DDMWIDTH E, else the GEOMECH E when enabled, else
20 GPa. G_c = K_IC²/E′.

*Verification (`-c28pf`).* 2 m³/d water into a 0.1 mD layer, η = 0.002, K_IC 1 MPa√m, σ 320 bar:
  * V = 0.0352 m²/m after 88 d gives L_PF = 32.1 m, against the toughness-regime Griffith–Sneddon length
    a = (E′V/(2√π K_IC))^{2/3} = 35.5 m (ratio 0.90, the AT2 l/a effect of §6.16).
  * Crack pressure: 321.1 bar = σ + K_IC/√(πa).
  * Monotone growth; wings opened 25 → 35 m; 6 s for 16 report steps.

*Limits (stated in the UI).*
  * One-way coupling by volume: the flow model's fracture pressure (392.7 bar in the test) is reported but does not drive the
    phase field, so η sets the length.
  * There is no fluid-pressure gradient inside the phase-field crack: it carries a uniform pressure, and the pieces carry
    the flow.
  * The geometry is 2-D plane strain (KGD-type). The state is not saved in restart files: after a restart the model starts
    from the initial crack and integrates the injection after the restart, while the open pieces are kept.
  * Copies of a Simulator share the link object.

**Sub-grid leak-off (research H2b, `HFGROW … SUBGRID LEAK`, `setHfTips`).**
  * *Model.* For a SUBGRID tip piece with open fraction 0 < φ < 1, the exposed part of its host cell is treated as its own
    storage: p_sub = p_rest + (p_cell − p_rest)/φ, where p_rest is the host pressure when the tip first entered the cell.
  * *Leak-off.* The host–fracture transmissibility becomes T₀·φ·m with m = clamp((p_f − p_sub)/(p_f − p_cell), 10⁻³, 1),
    lagged to the start-of-step pressures. m → 1 as φ → 1.
  * *Cycle 7 result (amended in Cycle 8, F24).* t60 = 25.49 / 26.34 / 26.00 d at 10 / 5 / 2.5 m vs 33.58 / 28.77 / 26.97 d without the correction (C28 engine). With the C29 engine (well-start repair) the same decks give 25.31 / 25.55 / 26.00 vs 28.67 / 27.42 / 26.84 d, and the dipping corner-point variant 23.00 / 22.80 / 23.65 vs 25.84 / 24.53 / 24.00 d
    (F21).

**Distributed assembly (C21-3, `Simulator::setDistributedAssembly`, `mpi_dd.cpp` `solveCollective`).**
  * *Partition.* Rank r owns the cell slab [a₀, a₁). It evaluates the properties and assembles only those Jacobian rows
    (connections with at least one owned cell; wells on the rank of their first completion).
  * *Lockstep.* The Newton state is replicated. The residual norms and the max-change/CFL quantities are reduced with
    `sumN` / `maxN`, so every rank takes identical decisions (convergence, cuts, well switching).
  * *Solve.* CPR-BiCGSTAB over the ranks with rank-local rows, no matrix scatter. The update is all-gathered.
  * *Results.* Bitwise identical across ranks and equal to serial (FOPT 5,478,750.61). On SPE9-like data (9,000 cells,
    2 ranks): 5.82 s vs 7.32 s with rank-0 assembly and 7.47 s serial. Linear iterations 1,175 vs 848 serial, from the
    block-Jacobi coupling of the distributed preconditioner.

**Level-ordered contiguous BILU0 (C23-b, `bsr.cpp` `buildContiguous`).**
  * *Renumbering.* The forward level sets of the block-ILU(0) graph give a topological order π. Rows are renumbered by π
    and the factor is stored contiguously per level (the value map keeps the original matrix positions). The factor sweeps
    and both triangular solves then stream through memory with one OpenMP loop per level.
  * *Exactness.* For a structurally symmetric pattern, ILU(0) in any topological order of the lower-triangular dependency
    graph gives the same factors (each row needs only its final predecessor rows). The factor is therefore bitwise equal
    to the serial one. The structural-symmetry check falls back to the plain path otherwise.
  * *Default.* On for n ≥ 20,000, ≥ 128 rows per level and parallel execution. `RESSIM_LEVEL_ILU` = 2 forces it on, 0 or 1
    turns it off.
  * *Measurement (54k-cell block Laplacian, nb 3, 2 threads).* 19 iterations in all variants; apply time 2.18 ms serial,
    2.39 ms level-scheduled, 2.09 ms contiguous (464 levels, 116 rows/level). The gain is meant for ≥ 4 cores.

**FOSM prediction uncertainty and data worth (`study.cpp`, iTOUGH2).**
  * *Prediction variance.* After the match, C_p = s₀²(JᵀWJ)⁻¹ (§12.0). For each prediction z with sensitivity row g
    (central differences of the same runs): σ_z² = g C_p gᵀ, with the 95 % interval z ± 1.96 σ_z.
  * *Data worth.* Recompute C_p with the rows of one observed data set removed, or with a candidate data set added (its
    sensitivities from the runs, weight from the stated error), and report Δσ_z/σ_z.
  * *Monotonicity.* Removing rows can only enlarge (JᵀWJ)⁻¹ in the Loewner order, and adding rows can only shrink it.
    Removing the only data set makes the prediction unbounded. These properties are tested.
  * *Test numbers.* FOPT 8.337e6 ± 888 (independent linear check 911); removing the injector BHP +44 %; adding WGOR −44 %.

**Tracers (C28 additions, `tracer.cpp`).**
  * *Partitioning tracer* with coefficient K into the partner phase. The partner carries K·c·b_carrier/b_partner per surface
    volume, so the storage is (φ b_w S_w + K φ b_o S_o)c and the retardation is R = 1 + K S_o/S_w (PITT, Cooke 1971):
    S_or = (t_p/t_w − 1)/(t_p/t_w − 1 + K).
  * *Dispersion / diffusion.* Conductance per connection α_L|F|/d + (φ S b)_avg D_m A/d, implicit, so it remains an M-matrix.
  * *Decay chains (Bateman).* The parent's decayed amount D_p of the step feeds the daughter with the in-step survival factor
    g = λ₁(e^{−λ₁Δt} − e^{−λ₂Δt})/((λ₂ − λ₁)(1 − e^{−λ₁Δt})), which is exact for any step size. A stable end member conserves
    the total.
  * *Test numbers.* t50 ratio 1.504 (theory 1.5), S_or 0.2013 (true 0.2); chain 0.015625 / 0.3515625 / 0.6328125 exact.

**Geostatistics (`geostat.cpp`, GSLIB conventions).**
  * *Variograms.* Spherical, exponential and Gaussian with nugget and geometric anisotropy (major/minor ranges at an azimuth,
    vertical range).
  * *Kriging.* Ordinary kriging (Lagrange multiplier) or simple kriging in a moving neighbourhood.
  * *SGS.* Normal-score transform, random path, simple kriging on data plus previously simulated nodes, back-transform by
    the data quantiles with linear tails inside the data range.
  * *Test numbers.* Kriging is exact at the data; SGS on 80×80 gives mean −0.115, variance 0.961 and variogram error 0.042.

**Gaussian-process proxy (`fitProxyGP`).**
  * *Kernel.* k(u,v) = s_f² exp(−½Σ((u_k − v_k)/ℓ_k)²) on the unit cube, plus noise s_n². Constant mean; ARD length scales
    and s_f² by maximum likelihood (Cholesky, multi-start on log-parameters).
  * *Prediction.* k*ᵀα, variance s_f² − k*ᵀK⁻¹k*. Q² from closed-form leave-one-out residuals α_i/[K⁻¹]_ii.
  * *Test numbers.* Q² 1.000 vs 0.960 quadratic; 99 % coverage of the 95 % band.

**Production analysis (`prodanalysis.cpp`).**
  * *Decline.* Arps q = q_i(1 + bD_i(t − t₀))^{−1/b}, fitted by log-rate least squares from the peak, with closed-form
    cumulative and EUR.
  * *IPR.* Composite Vogel: J above p_b; below, q = q_b + (J p_b/1.8)[1 − 0.2 p/p_b − 0.8 (p/p_b)²].
  * *VLP.* Homogeneous no-slip mixture on the case PVT, hydrostatic plus Darcy–Weisbach friction (Haaland), marched in
    segments; the operating point is the root of IPR − VLP.
  * *VFPPROD.* Generated from the same model on the flow × THP × WCT × GOR nodes (exact at the nodes).
  * *Test numbers.* q_i 1006 / D_i 0.00208 / b 0.65 (true 1000 / 0.002 / 0.6); VLP 3053 vs Colebrook 3070.5 psi.
  * *Not modelled.* No hold-up or flow-pattern correlation, no surface-network solver.

**FAULTSLIP in the thermal and compositional engines (`coupled_deck.cpp`).** The headless coupled runner applies the same
per-report-step fault mechanics as the black-oil driver:
  * *Thermal.* Thermo-elastic horizontal stress change Δσ_h = E α_T ΔT/(1 − ν) (STRESSFIELD … TH E alpha nu, §6.19) on the
    faces, or the THM mechanics tensor with STRESSFIELD FE.
  * *Compositional.* Total stress from the FE geomechanics.

It is used by the GUI, ResSimCLI, `runSingle` (studies) and the MPI batch KPIs. EGS test: dCFS 95.3 bar, slip 8.18 mm,
61 expected events; CO2: uplift 3.72 mm.

### 6.22 C29: salt precipitation, wellbore friction / ICD, surface network, two-way phase field, device CPR-AMG

**Salt precipitation and dry-out (`-- RESSIM: SALT c0 [cSol [rhoSalt [phiR [Gamma [T_C [cInj]]]]]]`, `salt.cpp`).** This is
TOUGH3 EWASG physics in the black-oil engine. It is operator-split after every converged flow step:
1. *Dissolved NaCl* is the internal water tracer NACL (kg per m³ of water). It uses the implicit upwind transport of §9.1 with
   the converged water fluxes; the injected concentration is cInj (or a WTRACER NACL value).
2. *Solubility equilibrium per cell.*
   * Any excess over cSol (default 317 kg/m³) precipitates as halite.
   * Solid redissolves into undersaturated brine.
   * The solid saturation is S_s = m_solid / (ρ_halite PV₀).
3. *Dry-out.*
   * Gas that enters a water-bearing cell picks up water vapour to saturation (local equilibrium).
   * The humidity y of the gas leaving dried cells comes from inflow mixing. It is solved by Gauss–Seidel sweeps over the
     converged gas fluxes; injected gas is dry.
   * Evaporation rate: E = ω Σ_in q_g,res (1 − y_up), with ω = p_sat(T) M_w / (R T ρ_w) (IAPWS saturation line).
   * E is the water sink of the next flow step, applied as E·S_w/(S_w + 0.005). This is implicit in S_w, so no negative
     saturation can occur.
   * The salt stays behind, so the brine concentrates and precipitates.
4. *Rock property changes.*
   * Porosity: φ(1 − S_s), through the pore-volume reference.
   * Permeability: the Verma & Pruess (1988) tubes-in-series factor
     k/k₀ = θ²(1 − Γ + Γ/w²)/(1 − Γ + Γ[θ/(θ + w − 1)]²), with θ = (1 − S_s − φ_r)/(1 − φ_r) and w = 1 + (1/Γ)/(1/φ_r − 1).
     It is applied to connection transmissibilities (harmonic mean of the cell factors) and to well connection factors.
5. *Restart.* The solid and evaporation state are saved with the restart (state v5). Snapshots carry SSALT and KSALT.

*Verification (`-c29salt`), 20 × 1 × 1 column at 150 °C with dry gas injection:*
* S_w in the injection cell drops from 0.200 to 6.8e-4; without the salt model it stays at 0.200.
* S_s = 0.0227 and k factor 0.800, which equals the closed-form Verma–Pruess value to 1e-12.
* Salt mass is conserved (in place + solid + produced) to 1e-6.
* An initial brine of 400 kg/m³ is capped at 317, conserving mass to 1e-9.
* Restart matches to 1.8e-7.

*Limitations:*
* The vapour is not a component: it leaves with the gas and does not condense.
* Halite is the only salt.
* The PV and transmissibility factors are relative to the start of the run. A two-way GEOMECH or FAULTMECH permeability
  coupling of the same cells, or HF transmissibilities, override them (last writer wins).

**Wellbore friction and inflow-control devices (`-- RESSIM: WELLBORE well ID rough [icdK [rhoCal]]`, FlexWell-lite,
`simulator.cpp wellboreFriction`).**
* *Tubing path.* The tubing is the polyline datum → connection 0 → 1 → … in COMPDAT order (heel to toe). Each segment carries
  the reservoir rates of all connections beyond it (toe side), at their mobility-weighted density.
* *Friction.* Darcy–Weisbach with 64/Re (laminar) or Haaland (turbulent). The connection pressure shifts by the accumulated
  drop (producers +, injectors −). The drop is explicit and under-relaxed ½ with the previous step, like the hydrostatic heads.
* *ICD.* Δp = K ρ_mix/ρ_cal q². It is solved per connection from the local balance q = a(D₀ − Δp), where a is the reservoir
  productivity and D₀ the drawdown at the lagged BHP. The positive root u = 2D₀/(1 + √(1 + 4 c a² D₀)) is exact for the
  connection.
  * The explicit form Δp = Kρq²(lagged q) is unstable when the gain 2Kρq·a exceeds ~1. Relaxation does not cure this, so
    it is not used.
* *Verification (`-c29wb`).*
  * ID 0.1 ft: 42.2 psi of friction at the deepest connection; FOPT 4.409e6 → 4.387e6.
  * ICD with K = 2e-5 psi/(rb/d)²: drop 851.5 psi at 8274 rb/d, i.e. K ρ_mix/ρ_cal q² with ρ_mix/ρ_cal = 0.62.
* *Not modelled.* Transient wellbore storage, multiphase slip / hold-up, and annulus or concentric strings (full FlexWell
  / T2Well).

**Surface network (ECLIPSE NODEPROP / BRANPROP, `network.cpp`).**
* *Input.* The keywords stay verbatim, so they round-trip exactly and OPM Flow runs them too.
* *Node assignment.* Nodes are groups. A well belongs to the first network node found walking up its group chain.
* *Balance.* Explicit, once per time step (like NETBALAN with a zero interval).
  * Node rates are the latest converged rates.
  * Pressures go top-down from the fixed-pressure terminal node: p_down = VFP_branch(q, p_up, WCT, GOR, ALQ). VFP 0 means no
    loss; 9999 means the branch is removed.
  * Each producer gets its node pressure as its THP limit. updateThp converts it with the well VFP table (bisection on
    IPR = VFP), or with the mixture head to the datum when there is no table.
* *Verification (`-c29net`).*
  * Node G pressure = 600 + 0.05 q exactly (C29 explicit); since C30 within the implicit balance tolerance 0.2 %
    (767.99 vs 766.86 psi at q 3337 stb/d).
  * The producer is THP-limited at the node pressure, with BHP above its own 500 psi limit (C29 1858 psi, C30 1305.5 psi).
  * Less free gas near the well gives more late oil: FOPT 5.65e6 (C30 5.644e6) vs 4.65e6. The same trend appears with fixed BHP limits
    of 1000 / 1858 / 2500 psi: 4.56 / 4.80 / 5.25e6.
* *Limitations (C29).* One balance per step with lagged node rates; no compressors or chokes beyond VFP / ALQ. **The lagged
  balance oscillated on gassy decks (NETWORK.DATA) - replaced in C30 by an implicit balance, see §6.23.**

**Two-way phase-field coupling (`HFGROW … PHASEFIELD eta TWOWAY`).**
* *One-way baseline (C28, §6.21).* The crack volume comes from a fixed fluid efficiency, V = η V_inj/H.
* *TWOWAY volume.* From the second phase-field step on, the volume follows the flow model:
  * V = max(V_last, min(V_cap, C · p_net,flow)).
  * C = V/p_net is the compliance of the last phase-field step.
  * p_net,flow is the flow model's pressure in the open pieces minus closure.
  * The mass cap is V_cap = V_last + ΔV_inj/H: the crack cannot store more than the fluid injected in the step.
* *Sub-stepping.* The volume increment is applied in sub-steps (×1.5 at most each, at most 16) so that the staggered
  iterations follow a fast-growing crack.
* *Feedback to the flow model.*
  * Each piece's conductance multiplier is (w(x)/w_ref)³, from the Sneddon aperture w(x) = 4p_net/E′ √(a² − x²), clamped to
    [1e-3, 1e3].
  * This uses `setHfConductance`, shared with G1-c and DDMWIDTH: last call wins.
  * The crack never shrinks.
* *Verification (`-c29pf`).* 0.6 m³/d into the C28 test layer:
  * Two-way half-length 59.8 m vs 23.3 m one-way.
  * |p_crack − p_flow| at the last step is 1.26 bar vs 8.02 bar one-way.
  * 16 aperture updates.
* *Meaning.* η = 0.002 one-way under-predicts the crack whenever the flow model holds the fracture above the LEFM
  pressure. Two-way removes the assumed η.

**CPR-AMG on the device (`DeviceCprAmg`, `bsr.cpp`; D3D11 path `gpu_compute.cpp` mode 2, UI: HPC → GPU preconditioner →
CPR-AMG).**
* *Stage 1.* Quasi-IMPES weights W restrict the residual to pressure; E prolongs back.
* *Pressure V-cycle.* Smoothed-aggregation hierarchy in float32 CSR (`DevCsr`), damped-Jacobi smoothing, and an explicit
  dense inverse on the coarsest level (≤ ~150 unknowns, so it is one SpMV on the device).
* *Stage 2.* Block-Jacobi (inverted 3×3 diagonal blocks).
* *Accuracy and fallback.* Double-precision iterative refinement around the float solve, as for float BILU0 (C16). The CPU
  emulation `CpuFloatBlockSolver::cprAmg` runs the same algorithm and backs the tests. If the device solve fails, the
  solver falls back to the CPU.
* *Verification (`-c29gpu`).*
  * 120 × 100 block system: 23 iterations vs 43 with the device MC-BILU0, error 6.6e-11.
  * Demo Newton loop (6 steps): FOPT identical to the CPU CPR, **921 vs 5118 linear iterations** for the device MC-BILU0,
    no fallbacks.

**Well-solution repair (Newton update).** A rate-controlled well whose BHP leaves the flowing window gets zero rate and
zero derivative. For an injector the window is BHP above the lowest connection pressure minus head; for a producer it is
below the highest. With zero derivative the next update explodes and the iteration flip-flops. This happened when dry gas
was injected into a tiny EDFM fracture cell. Such a well is now put back inside the window: half-way between the edge and
the previous BHP. Water HF runs also need fewer cuts (H2c CP deck: 16 vs 144). The edge is the lowest connection value for
injectors and the highest for producers (the first connection that can flow); taking the maximum for injectors broke three
phase-field tests. Side effect: time-step paths change, reference values shift by ~2 %, and the old sub-grid-tip (SG)
grid study changes a lot (t60 at 10 m 33.58 → 28.67 d). F21 is amended accordingly (research Cycle 8, F24).

**Dynamic grid refinement (DynaGrid) — deferred, with the reason.** The engine's grid is structured (Cartesian /
corner-point plus LGR / EDFM cells) with a static active set. Run-time amalgamation or refinement would rebuild the active
set, connection list, Jacobian pattern, well connections, AMG hierarchy, restart layout and the cell maps of every
coupled module (tracers, salt, EDFM, geomechanics, fault mechanics). A static LGR around fronts (C14) and the
sub-grid HF tip models (H2 / H2b) cover the main use cases. Design note: ROADMAP C29.

### 6.23 C30: implicit network balance, fracture / phase-field field output, geomechanics post-processing, colour maps

**Why the C29 network run oscillated (user report, NETWORK.DATA).** FOPR jumped 3493 → 49 → 2584 → 20 stb/d between
report steps from ~3.5 y and the wells "died" at ~6.5 y (FOPT 7.12e6). There were two causes, both fixed:

1. *Wrong ratios at zero rate.* updateThp evaluates the tubing curve BHP = VFP(q, THP, WCT, GOR). When the trial BHP
   reached the shut-in pressure, q = 0 and the code set GOR = WCT = 0. A gas-free column is the heaviest lift curve, so
   exactly at q → 0 the VFP BHP exceeded the reservoir pressure: no intersection, and the well was declared unable to
   flow. At a vanishing drawdown the ratios are now their limits, the rates at BHP − 10⁻³ p.
2. *Explicit balance with lagged rates.* Node pressures came from the rates of the previous step,
   p_node = G(q_{n−1}), and the wells then flowed at p_node.
   * Linearised, the loop is δp_{n+1} = −(dG/dq)(dq/dp) δp_n with a gain ≈ 0.74 here, rising to > 1 as the GOR
     increases.
   * The result is a period-2 limit cycle: a high-rate step gives a high branch loss, which gives a low rate, then a low
     loss, then a high rate.
   * The C30 balance is implicit. After a converged step the node pressures G(p) are recomputed from the rates the
     wells delivered at p. The step is re-solved until |G(p) − p| ≤ 0.2 % per node (at most 12 trials).
   * G is decreasing in p, so the residual r(p) = G(p) − p has one root per node. The update is regula falsi with the
     Illinois modification, and the first trial is the fixed-point value.
   * `RESSIM_NET_EXPLICIT=1` restores the explicit balance for diagnostics.

*Result (NETWORK.DATA).*
* Both wells hold 3000 stb/d until ~730 d.
* After that the decline is smooth: FOPR 3968 / 2790 / 1040 / 56 stb/d at 3.5 / 5 / 7 / 10 y, node PLAT 453 → 302 psi.
* FOPT 1.0544e7 (explicit with fixed ratios 1.0536e7) vs 1.397e7 without the network; 86 steps, 2.5 s.
* `-c30net` checks no rate reversal after the plateau, a smallest step ratio > 0.5, and wells alive at 7 y. The
  measured values are 0 reversals and a smallest ratio of 0.90.

**Fracture field output (`hf_field.cpp`, `-- RESSIM: HFFIELD ON|OFF [E nu] [radius] [budget]`).** One-way and display only:
it does not call `setHfConductance`, so DDMWIDTH / PF TWOWAY keep the feedback.
* *Opening.* After every report step every piece gets p_net = p_piece − σ_c, where σ_c is the FESTRESS FE minimum
  stress or σ_min + ∇σ (z − z_mid). The open pieces solve K_o w = −p_net with contact (w ≥ 0). K_o is the dense
  opening matrix of the rectangular pieces, each built from two Nikkhoo–Walter triangular dislocations. The domain is a
  half space when the fractures are shallower than twice their extent, otherwise a full space.
* *Matrix field.* Matrix cells within radius × the largest fracture extent are evaluated nearest first, up to the
  budget (cells × triangles, default 6 M).
  * Quantities: stress change, displacement and volumetric strain from the superposition of all open triangles.
  * Cells whose centre lies in a fracture plane are averaged at ±¼ cell width along the normal, because the TDE field is
    singular in the plane.
* *Arrays.*
  * HFWIDTH (mm), HFPNET, HFCOND.
  * HFDSXX / YY / ZZ / XY, HFDSNN (normal-stress change on planes parallel to the nearest fracture: stress shadow),
    HFDSM, HFDTAU.
  * HFDISPX / Y / Z, HFDISP (mm), HFSTRV (µε).
* *Auto mode.* On for HFGROW fractures with K_IC > 0 that are not phase field, and for DDMWIDTH fractures. ON forces
  every fracture, OFF disables.
* *Elastic constants.* E, ν default to HFFIELD → DDMWIDTH → PF → GEOMECH E[0] → 20 GPa / 0.25.
* *Cost (HFGROW_INJECTION).* 21 pieces and 705 cells: 0.015 s per update.
* *Verification (`-c30hff`).* Probes just off the widest open piece:
  * σ_nn = −3.221e5 Pa vs p_net 3.22e5 Pa (traction condition);
  * displacement jump 0.4715 mm vs DDM opening 0.4735 mm;
  * the deck round trip and the arrays are checked as well.

**Phase-field field output.** `phasefield::Input::onField` delivers each completed load step's field on its
rectilinear mesh, in the local crack frame:
* node phase field v, displacement u;
* element stress change (degraded by v² + k);
* the crack-opening profile w(x) = ∫ u·∇v dy.

The Simulator keeps one field per PfRecord. Beyond 48 fields every second older one is dropped, and the latest is
always kept.

*Arrays on the flow grid.* PFDAMAGE (1 − v), PFDSXX / YY / XY, PFDSNN, PFDISP (mm, intact nodes v > 0.5 only; the
displacement of broken nodes is not physical) and PFWIDTH (mm, COD at the piece). Phase-field fractures are excluded
from the automatic DDM field, because a DDM solve under the flow model's net pressure would contradict the phase-field
opening (7.6 mm / 59 bar vs 0.68 mm). `HFFIELD ON` forces them in.

*Verification (`-c30pf`).*
* ∫w dx / V = 1.0000 on every step (consistency of the output).
* The first-step centre opening is 1.20 × Sneddon 4 p a / E′, the known regularisation widening at ℓ / a = 0.2.

*HFGROW_PHASEFIELD.* PFDISP max 0.32 mm, PFWIDTH max 0.68 mm, PFDS ±0.95 bar (Sneddon p_net ≈ 1.2 bar at a = 32 m).

**Geomechanics post-processing for any run (`geoseq::postProcess`, `geoseq::mechMaps`).**
* *Load.* The snapshots of a finished run load the FE mesh of the grid: Cartesian bricks or corner-point isoparametric
  hexahedra, with the padding and boundary conditions of `geomech::Input`.
  * Sources: ResSim black oil / compositional / thermal, or snapshots read from an ECLIPSE / OPM UNRST.
  * Each snapshot loads its pressure change from the first snapshot, plus a TEMP change when present (thermal stress).
* *Solve.* One FE solve per snapshot in time order, path-dependent with plasticity.
* *Arrays (the same as the in-run coupling).* UPLIFT, DISPX / Y / MAG (m), STRAINV, SIGEFFM / V / H, SIGMIN, SQRTJ2
  (pressure units, compression +), FAILURE (Mohr–Coulomb index), POROM, KDAMAGE, and EPSPLAS / EPSPV with plasticity.
* *Use.* It runs on a worker thread from Results → Fields & sections → Geomechanics (post), or from the result viewer's
  *Geomechanics* button for imported cases.
* *Verification (`-c30post`).*
  * Equals the in-run one-way coupler stepping with the report steps: 4.7e-10 m of 1.1e-2 m.
  * The first snapshot gives zero displacement.
  * Twice the pressure change gives twice the uplift (elastic).
  * The same result after writing the run as EGRID / UNRST and reading it back: difference 0.

**Colour maps (`app/src/colormaps.cpp`).**
* 15 maps: Turbo, Viridis, Jet, CoolWarm, Grey, Plasma, Inferno, Magma, Cividis, Rainbow, RdBu, Spectral, Seismic,
  Terrain, Hot. Each can be reversed.
* Used in the 3-D viewers, the section maps, the fracture panels and the phase-field plots.
* `-c30cmap` checks: rgb in [0, 1], reverse = map(1 − t), continuity (largest step 0.024 between 513 samples),
  monotonic luminance of the perceptual maps, and a light centre for the diverging maps.

### 6.24 C31: user-defined quantities (UDQ), UDA targets and ACTIONX in the black-oil engine; well section view (`udq.cpp`, `sim_udq.cpp`)

Until C30, UDQ and ACTIONX were kept verbatim for OPM Flow / ECLIPSE but ignored by the ResSim engine, so decks with run-time
logic gave different answers in ResSim. C31 evaluates them.

**UDQ (`udq::Set`).**
* Records: `DEFINE name expr /`, `ASSIGN name [selector] value /`, `UPDATE name ON|OFF|NEXT /`, `UNITS name 'text' /`.
  Variable types follow the first letters: FU (field scalar), WU (one value per well), GU (one value per group).
* Expressions are parsed by recursive descent into a tree with this precedence (low → high):
  comparisons `== != < > <= >=` (results 1 / 0), `+ −`, `* /`, unary minus, `^` (right associative), then primaries.
  Primaries are numbers, summary quantities (F*; W* with an optional `'pattern'` or well name; G* with a group name),
  other UDQs, TIME / DAYS / DAY / MNTH / YEAR, parenthesised expressions and functions.
* Values are scalars or vectors over wells / groups. Binary operations broadcast a scalar onto a vector. Undefined
  entries (a well outside the pattern, division by zero, log of a non-positive number) are NaN and propagate.
* Functions:
  * elementwise: ABS EXP LN LOG NINT SQRT DEF UNDEF IDV;
  * ranking: SORTA SORTD (rank 1..n of the defined entries);
  * set → scalar: SUM, AVEA (arithmetic), AVEG (geometric), AVEH (harmonic), MAX, MIN, PROD, NORM1, NORM2, NORMI.
    Set functions skip undefined entries.
* *When it is evaluated.* All DEFINEs active at the report step are evaluated in deck order after every accepted time
  step. UPDATE OFF freezes a value; NEXT evaluates once, then freezes. At the start of a report step only the ASSIGN
  records of that step are applied.
  * The reason: a restarted run has no previous-step summary, so re-evaluating DEFINEs there would read zeros and give
    a different UDA value than the continuous run (found by the restart test).
* *Field quantities.* The F / W / G vectors a UDQ can read are the summary vectors of the run. C31 adds the reservoir
  voidage rates FVPR / FVIR (Σ q_res over producers / injectors, reservoir barrels or rm³ per day). These are needed for
  voidage-replacement logic and are computed from the well reservoir rates. ECLIPSE has FVPR / FVIR as well.

**UDA (user-defined arguments).** A UDQ name in a numeric item of WCONPROD (ORAT WRAT GRAT LRAT RESV BHP THP), WCONINJE
(RATE RESV BHP THP) or WELTARG becomes the target.
* The importer stores (target, name) pairs in `WellControl::uda` and writes the names back on export.
* `validateCase` accepts a UDA in place of a missing numeric target.
* After each accepted step the target is set to the UDQ value: the well's own entry for a WU quantity, the scalar for
  an FU quantity. It applies from the **next** time step. The step being solved never changes its own target, so
  Newton still sees a fixed well equation.
* A NaN value keeps the previous target.
* Under a group control the value replaces only the well's own limit (`own`), not the group-assigned share.
* *Verification.* The injector rate equals FUINJ of the previous step within the well-equation tolerance: largest
  relative difference 4.1e-6. The same-step FUINJ differs by a factor of 21 in the step where the target jumps, which
  proves the lag.

**ACTIONX.**
* *Parsing.* The header gives name, maximum runs and minimum interval (days). The condition lines have the form
  `lhs op rhs [AND|OR] /`; the trailing connector may also come without the slash. `(` `)` group terms. The body runs to
  ENDACTIO.
* *Precedence.* AND binds tighter than OR.
* *Well sets.* A comparison on a W quantity is true when at least one matching well satisfies it. That set of wells is
  carried along: AND takes the intersection and OR the union of the well sets. `'?'` in the body expands to one record
  per triggering well.
* *Other operands.* Months (JAN..DEC), DAY / MNTH / YEAR, TIME, numbers, summary quantities and UDQs.
* *Timing (OPM Flow practice).* Actions are evaluated at the end of every report step in which they are active and
  allowed (runs < max and the interval has elapsed). The body is applied through `applyScheduleText` and takes effect
  from the next report step.
  * Supported body keywords: WELOPEN, WCONPROD / WCONINJE, WELTARG, GCONPROD / GCONINJE, COMPDAT on existing
    connections, UDQ.
  * New connections are not added at run time; the engine warns instead.
* *Restart.* The state file (version 6) stores run counts, last-run times, the applied bodies (replayed on load) and
  the UDQ values.
* *Import.* UDQ / ACTIONX are now information messages ("evaluated by the ResSim engine"). PYACTION stays a warning.

**Example `ACTIONX_UDQ.DATA`** (11 × 11 × 3, SPE1 fluids, 4 producers + 1 injector).
* UDQs: FUWCT = FWPR / (FWPR + FOPR); injector RATE = UDA FUINJ = 1.05 FVPR + 500; producer ORAT = UDA FUPLIM
  (ASSIGN 1500).
* Action WCUT: `WWCT 'P*' > 0.8` → WELOPEN '?' SHUT, up to 4 runs.
* Action PRESSURE: `FPR > 5300 AND DAY >= 1` → ASSIGN FUPLIM 2200 and WELTARG I1 BHP 5800.
* Result: PRESSURE fires at 330 d. WCUT shuts P1 at 870 d and P2–P4 together at 1110 d (one run, three wells through
  `'?'`). FOPT is 8.316e6 stb.
* A restart after step 19 (after the actions fired) gives an identical FOPT.

**Verification (`-c31expr`, `-c31deck`).**
* Expression checks: precedence, broadcasting, NaN propagation, set functions, SORTA / SORTD, UPDATE NEXT / ON, save /
  load, ACTIONX parsing, well-set AND / OR, `'?'` expansion and malformed conditions. Example values: FUAVG = 458.258
  (AVEG), FUMX = 8.36364, NEXT/ON sequence 0 10 10 30 40.
* Deck checks:
  * UDA import / export round trip;
  * FUWCT equals FWPR / FLPR at every step (difference 0 over 43 points) and appears in the ECLIPSE summary output;
  * UDA lag;
  * run limits;
  * ASSIGN inside an action;
  * bitwise restart;
  * the import report classifies UDQ / ACTIONX as information.

**Well section view (UI, `ui_fieldviz.cpp`).**
* *Curtain.* The plan polyline is the deduplicated well path (wellhead, then the connection-cell centres). A vertical
  well gets an east–west line through the head. The polyline is extended at both ends by m × the cell size and sampled
  every Δx/4.
* *Cells.* Each sample takes the nearest column centre (the corner average of the first layer with corner data).
  Samples farther than 1.5 cell sizes from every column are gaps. Runs of the same column form segments [s₀, s₁]. Each
  active cell k of a segment is drawn as a quad between the averaged top and bottom corner depths. The segment list is
  cached per well / margin / grid.
* *Trajectory.* Path points and connection pieces are projected onto the curtain coordinate (nearest point of the
  polyline).
* *Log tracks.* The selected arrays in the connection cells, plotted against TVD (cell-centre depth) or MD (cumulative
  3-D distance from the wellhead along the connection sequence), for this step and the first step.

### 6.25 C32: multisegment wells, run-time ACTIONX connections and the geomodelling workflow (ResGeo) (`msw.cpp`, `sim_msw.cpp`, `sim_udq.cpp`, `geomodel_*.cpp`, `seismic.cpp`)

**Multisegment wells (`msw.h`).** WELSEGS / COMPSEGS / WSEGVALV / WSEGSICD / WSEGAICD stay verbatim in the deck (exact
round trip; OPM Flow runs the same deck). The engine parses the definitions in force at every report step.
* *Segment tree.* Segment 1 is the top node at the BHP reference; every other segment has an outlet. Node depth D and
  length L along the tubing come from WELSEGS (INC increments are converted to absolute values; ABS ranges are
  interpolated between the outlet and the range end).
* *Connections → segments.* COMPSEGS gives each connection the node of its branch nearest to the centre of the
  perforated interval (the ECLIPSE rule). Unlisted connections take the nearest node by depth (with a warning).
* *Segment flow and mixture.* q_s = Σ (reservoir inflow of all connections in the subtree of s). Density and viscosity
  are inflow-weighted; with no flow, the mobility-weighted cell values are used.
* *Pressure.* p_s = p_out + ρ_mix g (D_s − D_out) + sgn(q) (Δp_f + Δp_dev).
  * Friction is Darcy–Weisbach, Δp_f = f ρ v² |L_s − L_out| / (2 d), with f = 64/Re for Re < 2300 and
    Haaland otherwise: 1/√f = −1.8 log₁₀[(ε/d / 3.7)^1.11 + 6.9/Re]. WELSEGS item 6 'H--' switches friction off.
  * WSEGVALV: Δp = ρ v² / (2 C_v²), with v = q / A_c.
  * WSEGSICD (OPM form): Δp = K (ρ/ρ_cal)^¾ (μ/μ_cal)^¼ q².
  * WSEGAICD: Δp = a (ρ/ρ_cal) (μ_cal/μ)^y q^x, with x and y from items 13 and 14. Before C32 the parser read items 11
    and 12; this is now fixed and checked by a test.
* *Device coupling.* The subtree inflow responds linearly, q = A (D₀ − Δp), where A is the subtree productivity and D₀
  the drawdown without the device. Δp solves Δp − c (A (D₀ − Δp))^x = 0 on [0, D₀] by bisection. The device enters the
  well equation as a conductance factor (D₀ − Δp)/D₀ on all connections of the subtree, multiplied along the path to
  the top.
  * A lagged pressure shift of the size of the drawdown was unstable: it blocked or flooded connections whenever the
    BHP moved.
* *Relaxation and shut valves.* The friction and device parts are under-relaxed by ½ against the previous step; this
  state is part of the restart. A SHUT device closes every connection in its subtree.
* *Not modelled:* acceleration pressure drop, drift-flux slip, the SICD/AICD emulsion and critical-water terms,
  WSEGLABY / WSEGPULL, segment storage.
* *Verification (`-c32msw`, 15 tests):*
  * the friction, hydrostatic and valve drops equal hand calculations;
  * a wide frictionless tubing reproduces the standard well;
  * SICDs equalise inflow (the coefficient of variation drops ≥ 30 %) and delay breakthrough (60 → 135 d);
  * a SHUT valve isolates its subtree;
  * restart is bitwise identical;
  * WSEGAICD items 13 and 14 are read correctly.

**ACTIONX new connections (`Simulator::rebuildActionConnections`).** When an action body contains COMPDAT for
completions that have no connection yet, `refreshConnectionWindows` rebuilds the connection set of every well, using
the same code as setup (Peaceman factors; LGR, dual-porosity and EDFM targets; report-step windows).
* The dynamic well state is kept: control, mode, BHP, group and economic state.
* Per-connection state (q_w, polymer factors, friction part) follows its cell. New connections start from the cell
  values, and the hydrostatic heads are re-initialised from the cell densities.
* The Jacobian pattern is rebuilt (the cell rows gain well columns), and the AMG hierarchy is marked stale.
* Restart: the action replay at load time rebuilds the same set before the saved state is read again.
* Not supported at run time: new wells from WELSPECS, runs with multisegment wells, and domain-decomposed runs. These
  give a warning; OPM Flow applies them.
* *Verification (`-c32actx`):* ACTIONX_UDQ with P2 completed in layer 1, deepened to layers 1–3 by
  `ACTIONX DEEPEN: FOPT > 1E6`.
  * The action fires at 210 d. P2's BHP at its rate target rises from 4951 to 5045 psia.
  * The run equals the deck with the same COMPDAT written after report step 6 to a relative 5e-12 (P2 cumulative and
    FOPT).
  * Restart after the action is identical.

**Geomodelling workflow (`geomodel.h`, the Geomodel workspace / ResGeo).** This is a Petrel/RMS-style chain in the
project length unit, with depths positive down.
1. **Data.**
   * LAS 2.0 logs (wrapped or unwrapped, NULL values).
   * Well heads; deviation surveys (Petrel well-path files, or MD/INC/AZI columns via minimum curvature: dogleg β,
     ratio factor RF = (2/β) tan(β/2)); well tops.
   * IRAP classic and ZMAP+ surfaces; XYZ points gridded by ordinary kriging.
   * Fault sticks.
   * SEG-Y rev 1, IBM or IEEE floats. Inline/crossline and CDP X/Y come from the trace headers, with the coordinate
     scalar applied; the plan geometry is fitted from the first and last traces.
2. **Seismic interpretation.**
   * *Seeded auto-tracking* is a breadth-first flood over traces. At each neighbour it takes the extremum of the chosen
     event (peak, trough, or + → − zero crossing) within ±window of the parent time, and accepts it when
     |Δt| ≤ maxJump. For peaks and troughs it also requires an amplitude ≥ 35 % of the parent's; without this test the
     tracker crossed faults on noise.
   * Use one seed per fault block.
   * *Velocity model:* v(z) = v₀ + k (z − z_d), so t(z) = (2000/k) ln(1 + k (z − z_d)/v₀) ms TWT and
     z(t) = z_d + (v₀/k)(e^{k t/2000} − 1).
   * *Attributes* at depth points: amplitude, RMS or mean |a| in ±window.
3. **Structure and grid.**
   * *Lattice:* vertical pillars on a regular, optionally rotated plan lattice.
   * *Fault blocks:* each column gets a block from its side of every fault trace at the horizon depth; the trace comes
     from interpolating the sticks at that depth.
   * *Corner depths:* each column samples its corner depths inside its own block. A split node steps away from the
     pillar into the block, in steps of 0.5·k·a with a = 0.75·max(dx, dy), until it is on the block's side of the
     trace at its own estimated depth (dipping faults) and outside the map's ramp cell. It then extrapolates linearly
     back to the pillar. The ZCORN split nodes therefore carry the throw.
   * *Fault faces:* faces between blocks become FAULTS / MULTFLT records.
   * *Zones:* layered proportionally, parallel to the top, or parallel to the base. Cells thinner than minThickness are
     inactive.
   * *Well tie:* residuals at the tops are kriged (Gaussian variogram, range tieRange) and added to the horizon.
   * *Geometry:* bulk volume from the pillar quadrilateral (shoelace area) × thickness.
4. **Well blocking.** Each curve is averaged over the trajectory length inside every cell: arithmetic, geometric (PERM*),
   harmonic, or most-of (FACIES / LITH*).
5. **Facies.**
   * SIS: indicator kriging per code on a random path, conditioned to the blocked cells.
   * TGS: one Gaussian field, with thresholds from the cumulative proportions.
   * Object channels: sinuous ribbons y = A sin(2π x/λ) of width w and thickness h, inserted until the target
     proportion is reached; the wells are then forced to their facies.
6. **Petrophysics.**
   * Normal-score SGS (or kriging) per facies.
   * Optional collocated co-simulation with a secondary variable (Markov-1): y = ρ s + √(1−ρ²) r, where r is simulated
     with the residual data.
   * Poro-perm fit per facies: log₁₀ k = a + b φ, plus a N(0, σ) scatter.
7. **Volumetrics.**
   * STOIIP = Σ V_b · NTG · φ · (1 − S_w) / B_o over the cells above the OWC (GIIP above the GOC).
   * Optional Brooks–Corey saturation-height function: S_w = S_wirr + (1 − S_wirr)(P_d / P_c)^λ with
     P_c = Δρ g h above the contact.
   * Reported per zone, facies and fault block, with P90 / P50 / P10 over realisations (P90 is the 10th percentile).
8. **Upscaling.**
   * Coarse blocks are fx × fy columns by fz layers within each zone; faults stay only on coarse boundaries.
   * φ_c = V_p / (NTG_c · V_b,c), so the pore volume is conserved exactly. NTG is bulk-weighted.
   * k is arithmetic, harmonic, geometric, or flow-based: a local TPFA solve with a unit pressure drop and no-flow
     sides, solved by PCG.
9. **Export.**
   * GRDECL: SPECGRID / COORD / ZCORN / ACTNUM / properties / FAULTS / MULTFLT.
   * `toCase` builds the simulation grid and properties directly. Wells are sampled along their trajectories every
     min(dx, dy)/10; the dominant direction sets the completion direction. Names are truncated to 8 characters, and
     no controls are set: the UI adds LRAT producers.
* *Verification (`-c32geo`, 28 tests; demo project = faulted anticline with a 30 m throw, 6 wells including a
  horizontal one, and a SEG-Y volume):*
  * format round trips; ZMAP exact; minimum curvature exact on a circular arc;
  * auto-tracking picks 1271/1271 traces, depth-converted RMS 2.4 m against the map, whose own error is 4 m;
  * mean grid throw 27.5 m (true 30); well-tie RMS 2.32 → 0.04 m;
  * facies proportions within ±0.1; the co-simulation correlation rises from 0.05 to 0.56;
  * STOIIP equals the hand sum; the saturation-height function lowers it; realisation percentiles are ordered;
  * flow-based upscaling is exact for layers (k_x = arithmetic 505, k_z = harmonic 19.802);
  * coarse pore volume equals fine pore volume (1.10522e8 m³);
  * the geomodel runs in the black-oil engine through `toCase`.
* *ResSimCLI `--geo-demo [dir] [--field]`* writes the demo project, runs the whole chain headless and writes
  GEOMODEL_FINE.GRDECL and GEOMODEL_COARSE.GRDECL.

### 6.26 C33: scripting API, report generator, deck regression, opm-tests fixes and solver threading (`api.cpp`, `capi.cpp`, `report.cpp`, `regress.cpp`, `simulator.cpp`, `sim_udq.cpp`, `deck.cpp`, `linsolve.cpp`, `bsr.cpp`)

**Scripting API.**
* `api::Session` (C++) wraps one run: `open(deck)` / `setKeyword("KEY=VALUE")` / `run()`, or step-wise `begin()` →
  `step()` → `applySchedule(text)` → `step()` …, plus summary vectors, per-cell arrays of every report step, in-place
  volumes and the engine log.
* `applySchedule(text)` parses SCHEDULE keywords with the deck importer (`applyScheduleText`) and applies them from the
  current report step: controls, WELOPEN, WEFAC, WLIST, UDQ, new wells (WELSPECS + COMPDAT). This is the same path
  ACTIONX bodies use. A step-wise run reproduces the whole run exactly (bit-identical FOPT in `-c33`).
* `capi.h` exposes it as a plain C API (`rs_open`, `rs_begin`, `rs_step`, `rs_apply_schedule`, `rs_summary_get`,
  `rs_cell_array`, `rs_report`, `rs_compare_json`, …): libressim.so on Linux, ressim.dll on Windows.
* The Python package `python/ressim` is a ctypes wrapper of the C API. It needs no compiler and no pybind11, and uses
  numpy arrays when numpy is installed. It also covers ECLIPSE/OPM summary reading (`open_summary`), summary
  comparison (`compare`) and the report generator (`Session.report`).
* Examples: quickstart, parameter sweep, comparison with OPM Flow, and closed-loop steering (shut a producer when the
  water cut exceeds a limit).

**Report generator** (`report.cpp`; no Office, zlib or other external library).
* Writes DOCX (WordprocessingML with native, editable charts in `word/charts/chartN.xml`) and PPTX (one slide per
  section; tables split so that at most two text lines sit above a table; image slides).
* Also writes self-contained HTML (SVG charts, data-URI PNG maps) and Markdown.
* ZIP container with CRC-32. Raw DEFLATE uses LZ77 (32 KB window, hash chains) plus fixed Huffman codes; PNG is
  zlib-wrapped with Adler-32. The inflater is used for tests.
* Content: model description, initial in-place volumes, field and per-well charts, pressure / saturation layer maps of
  the first and last report step, engine messages, and an optional reference-run overlay with a relRMS deviation table.

**Deck regression** (`regress.cpp`, `ResSimCLI --regress manifest`, app *Run → Deck regression*).
* The manifest lists cases (deck + reference SMSPEC), keys and tolerances.
* Per vector: relRMS = RMS(ref − run) / max|ref|, with the run interpolated linearly in time onto the reference times.
  Vectors that are zero in both runs are ignored; a case where every common vector is zero passes.
* A case passes when every vector is within tol (default 5 %).
* Wall time is compared with the history (> 1.10 × the previous run is flagged). Outputs: HTML, DOCX, JUnit XML and
  history.tsv.
* `tools/fetch_opm_tests.py` downloads OPM/opm-tests (202 cases with OPM Flow references) and writes the manifest.
  Decks with keywords the engine does not model are listed as SKIP with the reason. The `nightly_regression` scripts
  run the whole chain.

**Engine fixes found by the regression.**
1. **WEFAC.** Well equations use flowing rates q. The reservoir source, group / field rates and cumulatives use
   e·q (ECLIPSE / OPM).
   * SPE9_CP_SHORT: WWPR 35 % → 2 %.
2. **Run-time wells and retyping.** WELSPECS / COMPDAT inside ACTIONX bodies or scripts create wells during the run.
   The well type follows the last control keyword (producer ↔ injector) instead of being fixed at import.
   * The import-time rule "periods of the other type are shut" now applies only to thermal / compositional decks.
3. **WLIST and full glob matching** (`*`, `?`, escaped `\*`).
   * List names start with `*`. A pattern that matches list names gives the union of those lists; otherwise it is a
     well glob.
   * ACTIONX_WELL_TEMPL-01…06 now pass.
4. **Stopped wells.** An open well with a zero rate target (e.g. WCONHIST with zero rates) has no inflow and no
   cross-flow. It now reports the zero-inflow BHP instead of 0, as OPM / ECLIPSE do:
   p_bhp = Σ_k WI_k (p_k − Δp_k^head) / Σ_k WI_k.
   * 0_BASE_MODEL2: WBHP:INJ1 relRMS 12.1 % → 1.4 %, WBHP:INJ2 13.3 % → < 1 %.
5. **MULTX- / MULTY- / MULTZ-.** These multiply the face to the minus neighbour. They are folded into MULTX / MULTY /
   MULTZ of that neighbour, which is exact for neighbour connections; fault NNCs use the low-index cell's multiplier.
   * 2_MULTXYZ_MODEL2: FWCT 25.7 % → 1.9 %.
6. **Solver threading.**
   * `setNumThreads(0)` ignored OMP_NUM_THREADS.
   * The BILU0 level-scheduled triangular solves ran in parallel on levels of only a few rows: barrier cost dominated
     on 2 cores.
   * Parallel gates are now n > 20 000 and n / levels ≥ 128 (BILU0), and kParMin 20 000 (vector kernels).
   * Results: ACTIONX_WCONHIST 48.4 s → 8.6 s; 0_BASE_MODEL2 320 s → 19.7 s. Found with a static gprof build, because
     dynamic gprof hides libgomp time.

**Verification (`-c33`, 17 checks; full suite 848 checks pass, MPI `-mpidd` pass).**
* CRC-32 check value; DEFLATE round trips (text, random, long runs, empty); PNG IDAT / Adler-32.
* DOCX / PPTX structure and size.
* Step-wise = whole run; WOPT = FOPT; mid-run WELOPEN SHUT.
* Manifest parsing; self-reference PASS; a +10 % FOPR deviation flagged (relRMS 6.4 %); writeOutputs.
* WEFAC rate and water balance (in-place change = FWIT − FWPT within 1 %); WLIST pattern; a run-time well P2 that has its
  own summary vectors.

### 6.27 C34: group guide rates, aquifers, directional relative permeability, rock compaction and other opm-tests items (`simulator.cpp`, `grid.cpp`, `deck.cpp`)

C34 worked through the opm-tests families that still failed after C33. Each item was checked against the OPM Flow
reference of a specific deck. The OPM source code was read where the manual leaves the behaviour open.

**Guide rates (GUIDERAT, GCONPROD items 9/10, GRCTRL decks).**
* The well guide rate is computed from the well potentials with the OPM formula
  GR = pot^A / (B + C·R1^D + E·R2^F).
  * OIL and RES phases: R1 = WOR, R2 = GOR.
  * LIQ: water and gas fractions of liquid.
  * GAS: WGR and OGR.
  * The guide rate is converted to the controlled phase with the current rate ratios.
  * Minimum interval, the allow-increase flag and damping are honoured.
* A group target is shared hierarchically. A subgroup with its own guide rate takes its share first, and its wells split
  that share by their own guide rates. The subgroup guide rate can be:
  * a value;
  * FORM, i.e. GUIDERAT applied to the summed potentials;
  * POTN, the potential;
  * otherwise, the sum of its wells.
* The group's own target and its other limits compete. The most violated constraint controls; for example, LRAT with an
  ORAT limit stays on LRAT.

**VREP / REIN coupled to the current production (NUPCOL-like).**
* OPM updates group targets inside the first NUPCOL Newton iterations. The lagged target in ResSim injected nothing in
  the first step and under-replaced voidage afterwards.
* After a converged step, VREP / REIN injection targets are now recomputed from that step's production. If a target moved
  by more than 1 %, the step is re-solved (at most two passes). Producer allocations are kept.
* Effect: KRNUM-02Z relRMS fell from 36 % to 8.9 %.

**Aquifers.**
* Numerical aquifers (AQUNUM / AQUCON):
  * The aquifer cells form a chain with NNC transmissibility 1/(1/T1 + 1/T2).
  * Reservoir faces connect to the first aquifer cell with T = Tc·Tcell/(Tc + Tcell)·mult.
  * The aquifer's initial pressure comes from AQUNUM.
* Carter–Tracy aquifers:
  * use the AQUTAB influence table (AQUCT item 11, piecewise linear);
  * AQUANCON item 9 replaces the face area.
* Results:
  * 2D_OW_CTAQUIFER: 67 % → 1.0 %.
  * 3D_2AQU_NUM: 5.0 % (pass).
  * 3D_1AQU_3CELLS: 24 % → 16 %.

**Directional relative permeability (KRNUMX/Y/Z, SATOPTS DIRECT).**
* The upstream mobility of each X / Y / Z connection is evaluated with the region of that direction. Wells use SATNUM.
* KRNUM-02X and KRNUM-02Y: 47 % → 4.3 % (pass). SPE1CASE2_KRNUM: pass.

**WELPI.** The connection factors present at the time of the keyword are scaled by PI_target / PI_current.

**Rock compaction (ROCKCOMP, ROCKOPTS, ROCKTAB, ROCK2D / ROCK2DTR / ROCKWNOD, OVERBURD, ROCKNUM).**
* Pore volume = ROCK compressibility term × table multiplier. ROCK cr still applies, as in OPM.
* Effective pressure:
  * Pe = p_o − overburden(depth), or p − p_initial with STORE;
  * IRREVERS uses min(p, p_min).
* p_min is +∞ until the first converged step, as OPM's minRefPressure.
* ROCKTAB is interpolated linearly with linear extrapolation.
* Water compaction (ROCK2D) uses bilinear interpolation in (Pe, max(Sw, Sw_max) − Sw_init), following OPM
  TabulatedTwoDFunction with the vertical policy.
* The transmissibility multiplier is applied to the upstream mobility of the cell. This covers connections and wells;
  OPM averages the two cells instead.
* Region: ROCKNUM if given, else the ROCKOPTS array (default PVTNUM).
* HYSTER is approximated by IRREVERS, with a warning.
* SPE1CASE2_ROCK2DTR: 1.5 % (pass).

**End-point arrays without ENDSCALE.**
* OPM Flow applies SWL / SWCR / … even without ENDSCALE in RUNSPEC (opm-tests model2 6_ENDSCALE). ResSim now does the
  same, with a warning, and the export adds ENDSCALE.
* Model2 chain: 6_ENDSCALE 16.5 % → 9.4 %; 7_HYSTERESIS 18 % → 3.8 % (pass); 8_MULTIPLY_TRANXYZ 15 % → 6.5 %.

**EDITNNC / EDITNNCR.**
* The transmissibility of a non-neighbour connection between two cells is multiplied (EDITNNC) or replaced (EDITNNCR).
  This applies to fault juxtapositions, pinch-outs and explicit NNCs; Cartesian neighbours are skipped, as in OPM.
* 9_EDITNNC_MODEL2: 34 % → 12 %.

**OPERATE / OPERATER.**
* These are folded into the arrays like EQUALS.
* Supported operations: MULTA, POLY, SLOG, LOG10, LOGE, INV, MULTX, ADDX, COPY, MAXLIM, MINLIM, MULTP, ABS and MULTIPLY.
* OPERATER selects cells by OPERNUM; OPERATE uses a box.

**WECON workovers.**
* WECON item 7:
  * CON closes the worst-offending completion (COMPLUMP lumps, else single connections);
  * +CON and PLUG also close every completion below it.
  * The well stays open while it has open connections.
* Other WECON items:
  * item 10 POTN applies the rate limits to potentials;
  * items 11/12 give the secondary water-cut limit and its workover.
* Closed connections persist through restarts.
* 3D_WECON: 72 % → 19 %.

**Restart state v7.** Adds WEFAC (WellRates.eff), guide rates, WELPI-scaled connection factors, econ-closed connections
and the rock-compaction history.

**Regression after C34.** The fetched opm-tests set grew to 307 cases (from 202). Against OPM Flow at 5 % relRMS:
**102 PASS, 147 FAIL, 5 ERROR, 53 SKIP**.
* New passes include:
  * all four GRCTRL MSW variants;
  * KRNUM-02X / 02Y and SPE1CASE2_KRNUM;
  * 2D_OW_CTAQUIFER, 3D_2AQU_NUM, AQUNUM-02, AQUCT-03 and AQUFET-01;
  * 7_HYSTERESIS_MODEL2;
  * 0B_ROCKTAB_MODEL6 and SPE1CASE2_ROCK2DTR;
  * NNC_AND_EDITNNC.
* WECON_PLUSCON_COMPLUMP now runs; it failed at t = 0 before, because OPERATER was missing.
* Still open:
  * GRCTRL STW variants (15–22 %, the water rate of PROD3);
  * ROCKTAB / ROCKTABH / ROCK_REGIONS_DEFAULT (25–43 %);
  * AQUFLUX;
  * KRNUM-03 (gas, 7–17 %);
  * the 9_* group / injection family (19–35 %);
  * 0B_WELPI (13 %);
  * 3D_WECON (19 %) and 3D_WECON_9 (56 %);
  * UDT-1D (convergence at 32.7 d);
  * polymer decks without PLYROCK.
* *Update C44:* 9_EDITNNC_MODEL2 now passes (3.0 %, then 2.0 %). 0B_WELPI 13.4 %, GRCTRL STW 9.9-14 % material. The model2
  well-PI and BHP-limit offsets traced back to WPIMULT applied from time zero and to FBHPDEF being ignored (§9.w).

**Verification (`-c34`, 7 checks plus the end-point scaling group; full suite 855 checks pass, MPI `-mpidd` pass, MinGW
0 warnings, test_api.py pass).**
* OPERATE / OPERATER change exactly the selected cells.
* EDITNNC multiplies only the NNC.
* A unit ROCKTAB reproduces the base run, and a varying table changes it.
* WECON CON closes completions.
* VREP injects in the first step: FWIR 2546 vs FLPR 2525.

### 6.28 C35: chemical EOR in the black-oil engine — surfactant, alkaline and foam (`chem.cpp`, `simulator.cpp`, `deck.cpp`, `ecl.cpp`)

**Scope.** The ECLIPSE options SURFACT (SURFACTW), ALKALINE and FOAM run in the own engine. They follow the ECLIPSE
Technical Description. Where OPM Flow has the same option, its definitions are used as well (OPM implements only the
tabular foam model). The keywords stay verbatim in `Case::extraKeywords`, so export / re-import is lossless.
* `SURFNUM` (REGIONS) and `SFOAM` (SOLUTION) are first-class arrays.
* Reference decks: `decks/chem/ALKALINE_ASP.DATA` (alkaline–surfactant–polymer, 1-D, FIELD) and
  `decks/chem/FOAM.DATA` (SPE1-type 10×10×3 gas injection with foam).

**Transport (one equation per chemical, carried by its phase: water for surfactant and alkaline, gas or water for foam
per FOAMOPTS).** For cell *i* with carrier pore volume $V_i = \phi S_\alpha b_\alpha V$ (surface volume of the carrier):

$$\frac{V_i^{n+1}c_i^{n+1} + A_i(c^{n+1}) - V_i^n c_i^n - A_i^n}{\Delta t} + \sum_j F_{ij}\,c_{up} + q_i c_{w}
 + \lambda_i V_i^{n+1} c_i^{n+1} = 0$$

* $A_i = (V_b-V_p)\,\rho_r\,\Gamma(c)\,m_{ads}$ is the adsorbed amount. $\Gamma$ comes from SURFADS / ALKADS / FOAMADS
  with the rock density from SURFROCK / ALKROCK / FOAMROCK. Index 2 means irreversible: $\Gamma = \max(\Gamma(c), \Gamma_{max})$.
* $m_{ads}$ = ALSURFAD(c_alk) for surfactant, 1 otherwise. ALPOLADS multiplies the polymer adsorption the same way.
* $F_{ij}$ is the converged upwind surface flux of the carrier phase.
* $q_i$ uses the well model's own connection rates: the water rates stored in the assembly, and gas-injection rates
  scaled to the well's gas rate. Injected mass therefore equals rate × WSURFACT / WALKALIN / WFOAM concentration exactly.
* Decay applies to foam only: $\lambda = \ln 2 / T_{1/2,w}(S_w) + \ln 2 / T_{1/2,o}(S_o)$ from FOAMDCYW / FOAMDCYO,
  integrated as $\mathrm{expm1}(\lambda\Delta t)/\Delta t$.
* Solution method: Newton on the concentrations with ILU0-BiCGSTAB, falling back to GMRES. Updates are chopped
  (Appleyard). Order: alkaline → surfactant (sees ALSURFAD) → foam.
* Coupling: the flow equations are solved fully implicitly. The chemical effects below are re-evaluated after each
  transport, and the outer loop repeats the step until they stop changing (the same loop as polymer).

**Surfactant effects.**
* Water viscosity: $\mu_w(c) = $ SURFVISC(c). The water mobility is multiplied by SURFVISC(0)/SURFVISC(c).
* Interfacial tension: $\sigma = $ SURFST(c$_s$) · ALSURFST(c$_a$). Pcow is scaled by $\sigma/$SURFST(0).
* Capillary number: $N_c = |K\nabla\Phi_w| / \sigma$ in SI.
  * Conversions: field $K\nabla\Phi$ in mD·psi/ft × 9.869233e-16·6894.757/0.3048; metric mD·bar/m × 9.869233e-16·1e5.
  * IFT: lb/in × 175.1268 → N/m.
  * Per cell, $|K\nabla\Phi_w|$ is the Euclidean norm of the X / Y / Z averages of the face values
    $T_f|\Delta\Phi_w|/(c_{Darcy}A_f)$.
* Miscibility $F = $ SURFCAPD(log$_{10}N_c$), clipped to [0, 1].
  * The critical end points blend: $S_{wcr} = F S^m_{wcr} + (1-F)S^i_{wcr}$, and the same for $S_{owcr}$.
  * $k_{rw}$ and $k_{ro}$ blend: $k_r = F\,k_r^{m}(S;\,EP) + (1-F)\,k_r^{i}(S;\,EP)$, with the immiscible table from
    SATNUM and the miscible table from SURFNUM, both evaluated on the blended end points. $k_{rg}$ is unchanged.
  * Without SURFNUM there is no miscible blend; a note is written.

**Foam effects.** Gas mobility is multiplied by $1-(1-M_c)(1-M_p)(1-M_s)$:
* $M_c$ = FOAMMOB(c$_f$) per PVTNUM;
* $M_p$ = FOAMMOBP(p$_o$), the pressure weakening;
* $M_s$ = FOAMMOBS(u$_g$), the shear weakening with the gas Darcy velocity.

FOAMOPTS FUNC is read as the tabular model (with a note). Deviation from OPM: OPM reduces the mobility of the foam
transport phase and ignores FOAMMOBP / FOAMMOBS and decay. We reduce the gas mobility and apply the full ECLIPSE
product.

**Reporting.**
* Component vectors follow the tracers in `WellRates::comp` / `SummaryPoint::fcomp*`, in the order ALKALINE, SURFACT, FOAM.
* ECL summary vectors: FTPR / FTPT / FTIR / FTIT + SUR / ALK / FOA, FSPR / FSPT / FSIR / FSIT, and WTPR / WTIR per well.
* Snapshot fields: SURFACT, ALKALINE, FOAM, SURFCNM (log$_{10}N_c$), SURFMISC (F), SURFSTM ($\sigma/\sigma_0$) and FOAMMRF
  (gas mobility factor).
* `Simulator::chemInPlace(i)` gives the amount in solution + adsorbed; `chemDecayed(i)` gives the cumulative decay.
* Restart state v8 stores the concentrations, adsorption history, multipliers and decay totals.

**Verification (`-c35`, 15 checks; full suite 878 checks pass after C35-3).**
* The ASP deck imports with 0 unsupported keywords.
* Injected masses equal rate × concentration × time to 1e-6: 1.0e6 lb alkaline and 6.0e5 lb surfactant.
* Mass balances close (injected = produced + in solution + adsorbed) to < 1e-8.
* IFT ratio reaches 7e-5 and F reaches 1.
* Surfactant raises FOPT from 26 384 to 28 174 stb (28 090 before the chemSolve line search, which converges the transport more tightly).
* SURFACT with zero surfactant / alkaline and SURFCAPD = 0 reproduces the deck without the option to 1e-7.
* Export → import gives identical results, and a restart after step 2 gives identical results.
* FOAM deck:
  * mass balance with decay closes to < 1e-15;
  * the gas mobility factor falls to 0.23;
  * the maximum injector BHP rises from 8289 to 9317 psia;
  * FGPT drops from 7.30e8 to 5.32e8 Mscf.
* **opm-tests SPE1FOAM against OPM Flow** (`decks/chem/SPE1FOAM.DATA` with the reference in `decks/chem/opm-ref`):
  relRMS is FOPR 0.07 %, FGOR 0.25 % and injector BHP 0.33 % (0.24 / 0.76 / 0.50 % before the chemSolve line search).
  * Before the rock-volume fix below, FGOR was off by 6.7 %.

**Unit fix found by this comparison (also affects polymer).** The adsorbing rock volume is $(1-\phi)V_b$, as in OPM.
* Before, it was computed as $V_b - V_p$, with $V_b$ in ft³ and $V_p$ in rb. In FIELD units that is
  $(1-\phi/5.61)V_b$.
* Rock densities follow the OPM keyword dimensions:
  * PLYROCK / FOAMROCK are densities (lb/ft³, kg/m³);
  * SURFROCK / ALKROCK are mass per reservoir volume (lb/rb, kg/rm³).
* Effect on the FIELD polymer example: it moves from 1.2 % to 1.5 % of its OPM value.
* The porosity handed to the sequential geomechanics coupling had the same rb/ft³ mix-up in FIELD units
  (`coupled_bo.cpp`, `coupled_deck.cpp`, `geoseq.cpp`); it now uses the true porosity.

**chemSolve robustness (C35-3).** Plain Newton on the chemical concentrations can cycle when the adsorption isotherm
saturates (FOAMADS / SURFADS reach their plateau at small $c$) in cells with a tiny carrier volume. In
SPE1CASE2_SOLVENT_FOAM, cell 1 flipped between $c = 0$ and $0.125$ with a flat residual. Each iteration now:
* evaluates the scaled residual $r = \max_i |R_i| / s_i$;
* takes the Newton step clamped to $\pm 0.5\,s_i$;
* runs a backtracking line search, halving $\lambda$ up to 12 times until $r(c + \lambda\,\delta c) < r(c)$;
* stops after at most 60 iterations.

Reducing the step limit alone did not help.

### 6.29 C35: SOLVENT, the four-component black-oil solvent model (`solvent.cpp`, `simulator.cpp`, `chem.cpp`)

**Scope.**
* Keywords: RUNSPEC `SOLVENT`; PROPS `PVDS` (solvent $B_s$ and $\mu_s$ against pressure), `SDENSITY` and `SSFN`;
  SOLUTION `SSOL`; SCHEDULE `WSOLVENT`, the solvent fraction $f_S$ of a gas injector's rate.
* Everything stays verbatim except SSOL, which is a grid array, so the deck round trip is exact.
* Combines with FOAM. `FOAMOPTS SOLVENT` makes solvent the foam transport phase.
* Follows OPM Flow `blackoilsolventmodules.hh`, immiscible part.
* Not modelled:
  * the miscible model (MISC, PMISC, MSFN, SORWMIS, SGCWMIS, SOF2, TLMIXPAR / TLPMIXPA, Todd–Longstaff). These keywords
    are kept, a note is written, and the immiscible model runs. OPM Flow runs them via Run > Engine;
  * solvent with CO2SOL / H2SOL (error);
  * gas–water systems without oil (the engine needs oil);
  * the MPI domain-decomposed run.

**Formulation: sequential-implicit composition inside a fully implicit flow solve.**
* **Flow unknowns.** The fully implicit equations are unchanged: water, oil and the total gas-like component G+S (free
  gas + dissolved gas + solvent, in surface volumes).
  * The gas-phase variable is $S_{gT} = S_g + S_s$.
  * As in OPM, the saturation functions, end-point scaling and hysteresis see $S_{gT}$.
* **Composition.** Each cell carries $z = N_s/(N_s+N_g)$, frozen during a flow solve.
* **Flash in `evalCell`, splitting $S_{gT}$.** With free gas present ($R_s = R_{s,sat}$):
  $$S_s = \frac{z\,(S_{gT} b_g + S_o b_o R_s)}{b_s(1-z) + z\,b_g}, \qquad S_g = S_{gT} - S_s$$
  When the gas phase would be pure solvent (undersaturated oil): $S_s = S_{gT}$, $S_g = 0$ and
  $R_s = S_s b_s (1-z)/(z S_o b_o)$.
* **Relative permeability.** With $k_{rgT} = k_{rg}(S_{gT})$:
  * gas: $k_{rg} = k_{rgT}\cdot \mathrm{SSFN}_g(S_g/S_{gT})$;
  * solvent: $k_{rs} = k_{rgT}\cdot \mathrm{SSFN}_s(S_s/S_{gT})$.
* **Solvent flux.** Gas-phase pressure, density $\rho_s = $ SDENSITY $\cdot b_s$, and an upstream mobility
  $k_{rs} b_s/\mu_s$. It is added to the G+S equation in both assembly paths, CSR and block.
* **Wells.**
  * Producers take $w_i\,\Delta p\cdot\lambda_s$ of solvent.
  * Injector mobility includes $\lambda_s$. A gas injector with WSOLVENT injects $f_S q_g$ solvent and $(1-f_S)q_g$ gas.
  * Reporting follows OPM: WGIR is the total gas-like injection, WNIR the solvent part. WGPR excludes produced solvent,
    and so does FGIP.
* **Composition equation.** Solved after each flow solve, implicit in $z$, with the flow state and potentials frozen:
  $$\frac{z_i N_{T,i} - N_{s,i}^n}{\Delta t} + \sum_{j\,\text{out}} T_{ij}\,\lambda_s(z_i)\,\Delta\Phi_{s,ij}
     - \sum_{j\,\text{in}} T_{ij}\,\lambda_s(z_j)\,\Delta\Phi_{s,ij} + q^{prod}_{s,i}(z_i) - f_S q^{inj}_{g,i} = 0$$
  * $N_T$ is the converged G+S amount.
  * Newton uses ILU0-BiCGSTAB with GMRES as fallback.
  * Converged when $|R|\,\Delta t / N_T \le 10^{-9}$; $\Delta z$ is clamped to 0.2; at most 40 iterations.
  * Cells without free gas (state 0) that receive solvent are switched to state 1.
* **Outer coupling.** The loop is shared with polymer and the chemicals. Order: solvent transport → re-evaluate → chemical
  transport. It repeats the flow solve until $z$ and the chemical effects stop changing.
  * If the solvent first appears in a cell during an outer pass, the flow is re-solved before the chemicals are transported.
* **Accepted state.** At step acceptance the cells are re-evaluated with the final $z$ (the last transport follows the
  last flow solve).
  * The old-time solvent amount stored is exactly the transported $z N_T$. Solvent mass is therefore conserved to
  round-off (8e-16 on SPE1CASE2_SOLVENT).
  * The flow accumulations use the consistent final split.
  * Restart state v9 stores $z$, the old-time solvent amounts and the cumulative totals, so a restart reproduces the
    continuous run bit for bit.
* **Initialisation.** SSOL gives $z = n_s/(n_s+n_g)$ from the initial solvent and gas amounts.

**Foam with solvent.** With FOAMOPTS SOLVENT, the foam is carried in the solvent: carrier volume $V = z N_T$ for
mass consistency with the transport. The foam factor multiplies the solvent mobility.

**Reporting.**
* Summary: FNIR / FNPR / FNIT / FNPT / FNIP, and WNIR / WNPR / WNIT / WNPT per well, written to ECL summaries.
* Snapshot: SSOLVENT and SOLVFRAC ($z$). SGAS is the free gas only, so SOIL = 1 − SWAT − SGAS − SSOLVENT.
* `Simulator::solventInPlace()` gives the solvent in place.

**Verification (`-solv`, 8 checks).**
* The three opm-tests SPE1CASE2 solvent decks import with 0 unsupported keywords.
* SPE1CASE2_SOLVENT (10 years of solvent injection) runs to the end. Max SSOLVENT is 0.65.
* Solvent balance: injected = produced + in place, to 8e-16.
* WSOLVENT 0 reproduces the run without the option exactly: FOPT 52 388 518.65, FGPT 387 430 661.7.
* Injecting solvent instead of gas lowers the produced hydrocarbon gas.
* **opm-tests SPE1CASE2_SOLVENT_FOAM against OPM Flow:** relRMS FGOR 0.38 %, injector BHP 0.88 %, producer BHP 0.57 %.
  At 10 years, WGPR:PROD is 25 234 vs 25 240 Mscf/d.
* Export → import gives identical FPR / FNIP.
* A restart after 30 report steps is identical: FPR 4947.624351 in both.
* Deck and reference files: `decks/solvent/` (README).

### 6.30 C35: gas-lift optimisation — LIFTOPT, WLIFTOPT, GLIFTOPT (`gaslift.cpp`, `simulator.cpp`)

**Scope.**
* Keywords (SCHEDULE, kept verbatim, exact round trip):
  * `LIFTOPT`: increment, minimum economic gradient, minimum interval between optimisations, optimise in all Newton
    iterations;
  * `WLIFTOPT`: well / pattern, YES / NO, maximum lift gas rate, weight, minimum lift gas rate, incremental-gas weight,
    allocate extra gas;
  * `GLIFTOPT`: group, maximum lift gas supply, maximum total gas (produced + lift).
* The ALQ of a well enters its VFPPROD lift curve $p_{bhp} = \mathrm{VFP}(q, p_{thp}, \mathrm{WFR}, \mathrm{GFR}, \mathrm{ALQ})$.
  * For ALQ type GRAT (or blank) the ALQ is the lift gas rate.
  * For IGLR / TGLR the rate is converted: IGLR $= q_{lift}/q_l$, TGLR $= (q_g + q_{lift})/q_l$.
* Fixed ALQ from WCONPROD item 12 / WELTARG LIFT is honoured as before.
* Follows the OPM Flow algorithm (`GasLiftSingleWell`, `GasLiftStage2`, `GasLiftGroupInfo`) in its objective and
  constraints. The search differs, see below.

**Lift performance of a well (frozen reservoir state).** At the start of every time step, for each THP-controlled
producer with WLIFTOPT YES:
* candidate ALQ values $a_j = j\cdot\Delta a$ (the increment) up to $a_{max}$:
  * $a_{max}$ = WLIFTOPT item 3, otherwise the largest value of the table's ALQ axis;
  * the last point is $a_{max}$ itself, at most 400 points;
* the operating point of the inflow (connection rates at the current cell state) and the lift curve at the THP limit
  (`thpOperatingPoint`) at each $a_j$, with BHP $\ge$ the BHP limit, gives $q_o(a_j), q_w(a_j), q_g(a_j)$;
* the well's own rate limits (ORAT / WRAT / GRAT / LRAT) cap the rates proportionally.

**Allocation.** A greedy allocation in whole increments. Each pass chooses the well and block $k > j$ (current
index $j$) with the largest weighted gradient
$$G = \frac{w_f\, s\, (q_o(a_k) - q_o(a_j))}{(a_k - a_j) + w_g\,(q_g(a_k) - q_g(a_j))}$$
($w_f$ = WLIFTOPT weight, $w_g$ = incremental-gas weight, OPM `calcEcoGradient`) and accepts it while $G >$ the
minimum economic gradient.
* **Look-ahead over blocks.** $k$ runs over all larger indices, the upper concave hull of the lift curve. A dead well
  that needs several increments before it flows therefore gets them as one block. This replaces OPM's separate
  "enough to flow" stage.
* **Minimum lift gas.** A positive WLIFTOPT minimum starts the well at that rate; a negative one starts it at the first
  flowing rate.
* **Constraints along the well's group chain** (efficiency-weighted sums, the wells not optimised at their current
  rates):
  * GLIFTOPT maximum lift gas supply and maximum total gas: a block that would exceed them is rejected (larger blocks
    too);
  * GCONPROD ORAT / WRAT / GRAT / LRAT targets and limits: the part of the oil gain that the group limit would cut back
    does not count, through
    $s = \mathrm{clamp}\big((Q_{cap}-Q)/(\varepsilon\,\Delta Q), 0, 1\big)$ (OPM `updateRatesToGroupLimits`).
* **Interval and fixed wells.** LIFTOPT item 3 (minimum interval) is honoured. WLIFTOPT NO with item 3 fixes the ALQ at
  that rate (NO 0.0 = no lift gas); NO without item 3 keeps the current ALQ.
* **After the allocation.** The THP-derived BHP limits are recomputed (`updateThp`), then the group targets are
  allocated on the potentials with the new lift gas.

**Differences from OPM.**
* OPM adjusts the ALQ one increment at a time per well inside the Newton iterations (stage 1), then redistributes
  between wells from the lowest decremental to the highest incremental gradient (stage 2), and finally removes surplus
  increments.
* ResSim solves the same discrete problem once per time step on the frozen lift curves: a greedy concave-hull
  allocation, which is optimal for concave curves under one shared limit.
* LIFTOPT item 4 (all iterations) is accepted but not used, since the allocation is per time step.

**Restart, reporting, API.**
* Restart state v10 stores the ALQ in force, the time of the last optimisation, the lift gas total and the run count.
* State v10 also stores the BHP limit of THP-controlled wells and the group-allocated rate targets carried over from
  the last time step. Before C35 the replay of the schedule reset them to the deck values, so a restart of any THP case
  differed slightly from the continuous run; it is now bitwise identical.
* Summary: FGLIR / FGLIT (efficiency weighted) and WGLIR / WGLIT, written to ECL summaries.
* `Simulator::wellAlq(i)` and `gasLiftRuns()` give the ALQ of a well and the number of optimisations.

**Verification (`-glift`, 11 checks; decks/gaslift, opm-tests).**
* GASLIFT-01 (9×9×2, six THP producers that do not flow without lift gas, WEFAC 0.8, VREP injection, FIELD ORAT
  10 000):
  * import: 7 gas-lift keywords, 0 unsupported;
  * FOPT 9.70e6 sm3 with lift gas vs 0.54e6 without, against about 10.25e6 in the commercial simulator (opm-tests
    plots, −5.3 %);
  * allocations are whole increments of 4e4 or the table maximum 4e5;
  * WLIFTOPT NO 0.0 gives zero lift gas from 2021-08-01 to 2021-10-01;
  * GLIFTOPT PLAT-1 9e5 holds: at most 8.96e5 after 2022-08-01;
  * an economic gradient of 10 allocates nothing;
  * a restart after 12 report steps is identical (FOPT 9 703 718.515, FGLIT 1 378 784 000);
  * export → import is identical.
* GASLIFT-04 (MODEL05 three-phase PVT, VFP with ALQ type GRAT), against the commercial simulator:
  * FOPT 1.2517e7 vs about 1.255e7 (−0.3 %);
  * FGPT 4.415e8 vs about 4.42e8 (−0.1 %);
  * FPR 119.0 vs about 116.5 bar;
  * the ORAT 10 000 plateau is held by lift gas.
* Known difference: in GASLIFT-01 our wells water out and die from about day 800 (ECL about day 1350). The cause is the
  late-time water distribution of this coarse 9×9×2 model, not the allocation; the same decks agree to 0.3 % in
  GASLIFT-04.
* OPM Flow publishes no reference runs for these decks. Its old plotted run does not lift the wells (FOPT about
  0.9e6).

### 6.31 C35: surface network and facilities — OPM semantics, satellites, consumption, sales, automatic chokes (`network.cpp`, `simulator.cpp`)

The C29 extended network (NODEPROP / BRANPROP, implicit balance) was rewritten to follow OPM Flow exactly. The checks
were opm-tests `network/` and `model2/`, with the OPM sources (`NetworkKeywordHandlers.cpp`,
`BlackoilWellModelNetworkPressureComputation.hpp`) as the specification.

**Topology**
* BRANPROP: VFP table 0 **removes** the branch; 9999 is a flow line with **no pressure loss** (p_child = p_parent).
  Item 4 is the ALQ value; item 5 (ALQ_SURFACE_DENSITY) is read.
* Standard network (GRUPNET), built from the group parent *when the keyword is read*:
  * terminal pressure ≥ 0 → fixed-pressure node, linked to the parent by a 9999 branch;
  * VFP ≤ 0 → branch dropped;
  * otherwise a branch to the parent with that VFP / ALQ;
  * item 6 (FLO) adds the lift gas to the branch flow.

  Standard and extended networks are not mixed. When GRUPTREE re-parents a group, its uptree branch moves to the new
  parent and keeps its VFP / ALQ (`Schedule::addGroupToGroup`). This is how NETWORK-01-REROUTE(_STD) works.
* NEFAC f: a node passes f × its flow up the tree (the efficiency of the node's group).

**Flows.** The inflow at a leaf node is the efficiency-weighted rate of the whole group subtree:

  q_leaf = Σ_w eff_w q_w + [Σ_w eff_w ALQ_w] + q_sat

The bracketed lift-gas term is used only when the branch adds lift gas. q_sat is the GSATPROD satellite rate
(oil, water, gas and its lift gas). Flows are accumulated from the leaves to the root. At each node the gas is first
reduced by the net consumption, q_g ← max(q_g − (C − I), 0), and then q_up += NEFAC × q_down.

**Pressures**, from the root down:
* a terminal node has its NODEPROP / GRUPNET pressure;
* a 9999 branch gives p = p_up;
* an automatic choke (NODEPROP item 3 = YES, together with GCONPROD for the choke group) gives p = THP of the group;
* any other branch gives p = VFP(q_o, q_w, q_g, p_up, ALQ).

A node pressure ≤ 1 atm is invalid (the wells are shut, as in OPM). Wells on THP control use the pressure of their
group node as the THP limit. The network and the well rates are balanced implicitly inside the step (C29 loop).

**Facilities**
* **GSATPROD** (satellite group): fixed oil / water / gas / RV rates and lift gas with no wells. They are added to
  FOPR / FWPR / FGPR and to the group totals, and they enter the network at the satellite's node.
* **GCONSUMP** consumption C and import I per group:
  * the net rate C − I lowers the gas sent to injection (REIN) and to sales;
  * item 4 = node: a named node, `'NO-NODES'` (no network effect), or defaulted (the nearest network node at or
    above the group, used only if it is a leaf);
  * FGPR is **not** reduced (as in OPM).
* **Gas balance:** FGSR = FGPR − FGIR − FGCR + FGIMR. FGCR / FGCT, FGIMR / FGIMT and FGSR / FGST are reported (also
  per group) and stored in restart v10.
* **GCONSALE** target S: the gas available to group g is

    A_g = Σ_sat q_g,sat + Σ_{w∈g} eff_w q_g,w − Σ_{h⊆g} (C_h − I_h).

  Gas re-injection becomes T = REIN · A_g, capped at A_g − S, so the sales target is met first. Items 3–5 (maximum /
  minimum sales rate and the procedure) are read and written; the procedure is not acted on.
* **Automatic choke (AS_CHOKE):** `chokeGroup` finds the network node whose group has a GCONPROD target. Each
  outer iteration moves the node pressure (bisection on the pressure, for the oil, water or gas target) until the
  group rate equals the target; the choke is then the drop p_node − p_parent (dpChoke, restart v10). If the target
  cannot be reached even with the choke fully open, the choke stays open.
* **GCONPROD exceed procedures:** items 11–13 (water / gas / liquid procedures; defaulted to item 7) are read and
  written. **NONE** means the limit is only reported, not enforced (OPM `groupLimitAction`). A RATE procedure caps
  the subordinate wells.

**THP operating point.** Network pressures change the THP limit of a well between Newton iterations. With low node
pressures, the VFP curve and the inflow curve can then cross twice, or touch at the lift limit. The operating point
uses the lowest-BHP stable root of a 48-sample scan. If the scan finds no sign change, a golden-section search finds
the minimum of BHP_VFP(q) − BHP_IPR(q) in a narrowed window. If the minimum is positive, the well is dead at that
THP.

**Results vs OPM Flow (relRMS, `-facil`):**

| Deck | Vector errors |
|---|---|
| NETWORK-01 and NETWORK-01_STANDARD | FOPR 0.14 % |
| NETWORK-01-REROUTE / REROUTE_STD | 0.28 / 0.30 % |
| NETWORK-GSATPROD-GCONSUMP | FOPR 0.61 %, FGCR 0.00 %, FGSR 0.67 % |
| model2 9_3C_GINJ_GAS_GCONSUMP_STW | FGSR 0.40 % (750 000 sales from 2019-06-01), FGPR 1.91 %, FPR 0.08 % |

NETWORK_MODEL5_STDW_AUTOCHK has no OPM reference run:
* B1 holds ORAT 6010 against a target of 6000, with a choke drop of 4.53 bar at the first report step;
* from day 60 the network cannot deliver the target and the choke opens.

**Limitations**
* Near the lift limit (an AUTOCHK variant with ORAT 4000), the rate is bistable. Ours stays below the target.
* The NETWORK-01-WTEST difference (9 %) comes from the WECON / WTEST threshold timing, not from the network.
* There is no NETBALAN tolerance control, no compressors (NETCOMPA) and no separator nodes. The GCONSALE procedures
  and the min / max sales limits are not acted on. GASLIFT-05..15 (gas lift in the network / MSW) run through OPM.
* The opm-tests decks are not in the repository (about 12 MB of includes). Fetch them with
  `python tools/fetch_opm_tests.py --dirs network,model2` or point `RESSIM_OPMT` at a copy; otherwise the reference
  tests are skipped.

### 6.32 C36: 3-D seismic — out-of-core SEG-Y, attributes, flattening, 3-D scene (`seismic3d.cpp`, `geoscene3d.cpp`)

**Out-of-core SEG-Y (`SegyVolume`).** Opening a file scans only the 240-byte trace headers: inline, crossline,
CDP X / Y (with the coordinate scalar). The lattice comes from the minimum inline / crossline and the most frequent
increments, giving (i, j) per trace. The plan geometry `x = x0 + i·ix + j·jx`, `y = y0 + i·iy + j·jy` is the
least-squares fit over all traces; the RMS residual is reported as `fitError`. Traces, inline / crossline / time
slices and sub-cubes are read from the file through a small LRU trace cache. Supported sample formats: 1 (IBM),
2 (int32), 3 (int16), 5 (IEEE) and 8 (int8), big-endian or byte-swapped. Missing traces read as zeros. `proxy(si, sj)`
builds a decimated in-memory cube in one sequential pass: a stride per plan axis, and box-averaging in time when the
sample count is reduced. This proxy is what the volume renderer uploads.

**Complex-trace attributes** (Taner, Koehler and Sheriff 1979). The analytic signal is `u + i·H[u]`, where the Hilbert
transform is computed by FFT on a power-of-two zero pad.
* Envelope: `A = sqrt(u² + H²)`.
* Phase: `φ = atan2(H, u)`.
* Instantaneous frequency: `f = Δφ / (2π·dt)`. Δφ is the wrapped phase difference of consecutive samples, i.e. the
  argument of `z_{t+1}·conj(z_t)`. A central difference of the unwrapped phase is biased, so it is not used.
* Cosine of phase; sweetness `A / sqrt(f̄)`; RMS amplitude over a window.
* Relative acoustic impedance: the running integral of the trace with its mean removed (low-cut).

Multi-trace attributes:
* Semblance coherence (Marfurt et al. 1998): `Σ_t (Σ_k u_k)² / (K·Σ_t Σ_k u_k²)` over a 3 × 3 trace aperture and a
  vertical window.
* Variance: `1 − coherence` normalised.

Whole volumes are streamed inline by inline into a new IEEE SEG-Y with the same headers.

**Flattening.** Each trace of a section is shifted so that the horizon sample `s_h(tr)` moves to the reference
`s_ref`, the mean horizon time, using linear interpolation. Traces with an undefined horizon stay unshifted.

**3-D scene.** Model coordinates are x east, y north and z the display depth (positive down), relative to a rounded
local origin to keep float precision with UTM coordinates. The vertical domain is either true depth or two-way time.
The time ↔ depth conversion uses the linear velocity law v = v0 + k·z:
* `TWT = 2000·ln(1 + k·z/v0)/k` ms;
* inverse `z = v0·(e^{k·TWT/2000} − 1)/k`;
* `TWT = 2000·z/v0` when k = 0.

The shaders apply the same law per sample, so time slices and sections are curved sheets in the depth domain.

The volume is ray-marched in a pixel shader over a 3-D float texture of the proxy cube. The opacity transfer function
is `α = opacity·smoothstep(threshold, 1, |a|/clip)`, composited front to back and stopped at the opaque-scene depth.
Horizons carry GPU contour lines from `fwidth` of the depth; derivatives are taken outside any varying branch.
Picking intersects the pixel ray analytically with the slice planes / sheets, and with horizons by marching plus
bisection on `z_ray(t) − z_h(x(t), y(t))`.

### 6.33 C36: geomodelling round 2 — structural editing, fault-parallel pillars, MPS, uncertainty loop (`geomodel_c36.cpp`, `geomodel_grid.cpp`)

**Horizon editing.** The brush weight is `w = exp(−3 r²/R²)` for r < R. Modes:
* raise / lower: `z += w·dz`;
* smooth: `z += w·s·(z̄₃ₓ₃ − z)`;
* flatten: `z += w·(z_target − z)`.

Whole-surface tools:
* Smoothing: a 3 × 3 binomial filter over the defined nodes. Interior nodes of a plane are unchanged, to 1e-13.
* Hole filling: Gauss–Seidel on the discrete Laplace equation over the undefined nodes; exact for planar data.
* Minimum thickness: top-down `z_h ← max(z_h, z_{h−1} + t_min)`, which removes crossing horizons.
* Conformable (isochore) horizons.

**Fault-parallel (sloped) pillars.** With `GridSpec.pillars = 1`, every pillar takes the plan slope of the nearest
fault plane. The slope is `s_f = (trace(z_mid + 50) − trace(z_mid − 50))/100` at the nearest trace segment, tapered
with distance: `s = Σ_f s_f·exp(−(d_f/L)²)` (L = taper distance, default 4 cells), normalised when Σw > 1.

The pillar position is `p(z) = p(z_mid) + s·(z − z_mid)`. A node depth must satisfy `z = h(p(z))`, which is solved by
fixed-point iteration (8 iterations; contraction factor |s·∇h| < 1). COORD stores the sloped end points. All geometry
(cell centres, `cellOf`, upscaling) is derived from COORD, so no extra state is needed. `cellOf` inverts the bilinear
pillar offsets by a fixed point (≤ 6 iterations).

Verification on the demo:
* pillar slope at the fault 0.48 vs the fault dip 0.50;
* 6e-22 at 2.8 km from the fault;
* horizons met exactly at the sloped pillar lines outside the 4-cell fault band;
* bulk volume −0.16 % vs vertical pillars;
* `cellOf(centre)` = cell for 7200 / 7200 cells.

**Direct Sampling MPS** (Mariethoz, Renard and Straubhaar 2010). The random path visits the uninformed cells. For each
one:
1. The n ≤ `maxNbr` nearest informed cells (hard data and earlier nodes) within the search radius give lag vectors.
   These are rotated into the training-image frame by the TI azimuth and scaled by the TI pixel size.
2. TI positions are scanned at random, up to `maxScan` of the TI, with mismatch distance `D = #{k : TI(x + h_k) ≠ v_k}`.
   Out-of-TI lags count as mismatches, and the scan aborts early once D reaches the current best.
3. The scan stops at the first candidate with `D ≤ ⌊t·n⌋`; otherwise the best candidate is kept.

**Servosystem.** The best candidate is tracked per centre category k. The category minimising
`D_k − λ·n·(p_k^target − p_k^current)` is chosen (λ = `servo`, default 0.5; target = the TI marginal or user
proportions). Without it, an incomplete scan favours the background category: the channel fraction was 0.17–0.23 vs
0.29 in the TI. With λ = 0.5 it is 0.285–0.292, unconditional and conditional.

Training images: GSLIB / SGeMS ASCII, or a generated sinuous-channel TI. Verification: hard data honoured exactly.
Indicator continuity `P(I(x + 12) = 1 | I(x) = 1)` along / across the channels:

| Model | Along | Across |
|---|---|---|
| TI | 0.44 | 0.26 |
| MPS | 0.38 | 0.25 |
| SIS with the same proportions | 0.30 | 0.32 (isotropic) |

**Uncertainty loop.** One realisation = structural field → grid → blocking → facies (SIS / TGS / channels / MPS) →
porosity → volumes. Inputs:
* triangular distributions for OWC, GOC, porosity and NTG multipliers, Sw / Swirr and Bo, with inverse CDF
  `F⁻¹(u) = a + sqrt(u(b − a)(c − a))` for u < F_c (and the mirror expression above);
* plain Monte Carlo or Latin hypercube sampling (one sample per stratum, independent permutations per input);
* new facies / porosity seeds per realisation.

**Structural field.** `δz(x) = σ·sqrt(2/M)·Σ_m cos(k_m·x + φ_m)`, with M = 256 and k_m ~ N(0, (√6/r)² I). This is
the spectral representation of the Gaussian covariance `exp(−3h²/r²)`. Simple kriging of the field values at the
(deduplicated) well tops is subtracted, so δz = 0 at the wells and the variance recovers away from them. Measured:
|δz| < 0.005 at 18 tops; sd 11.6 for target 10 beyond one range.

The same field is added to every horizon, keeping the zone thicknesses. Without a structural perturbation, the grid
and blocking are built once and shared by all realisations (OpenMP over realisations).

**Outputs.**
* P90 / P50 / P10 (P90 = the low case, 90 % exceedance) and the realisations closest to them, which are rebuilt
  exactly from their stored sample (STOIIP reproduced to 1e-12).
* A tornado: each input at its P10 / P90 quantile with everything else at the mode, a ±1.2816σ structural bar and a
  stochastic P90–P10 bar, sorted by swing.

Verified unbiased: the structure-only mean equals the base case within 0.3 %. With an OWC-only Latin hypercube, the
mean matches the linearised expectation.

### 6.34 C36: multisegment wells — labyrinth ICD, acceleration, drift flux, BHP reference depth (`msw.cpp`, `sim_msw.cpp`)

**WSEGLABY (device 4).** ECLIPSE does not publish the formula, and OPM Flow parses the keyword but does not simulate it.
So the model is rebuilt from the keyword items as a standard loss-coefficient device:
`Δp = ½ρv² [γ_in + γ_out + (n_ch − 1)·γ_lab + f(Re, R/D)·L_ch/D]`, with
* `v = q / (N_config·A_ch)`, `L_ch = L1 + (n_ch − 1)·L2`;
* f from Haaland (laminar 64/Re).

The defaults are the keyword's metric values (A 6e-5 m², L1 0.1863 m, L2 0.2832 m, D 0.006316 m, R 1e-5 m, γ 0.35 / 0.5 / 0.5).
They are converted for FIELD decks. A zero area or diameter is treated as SHUT.
The device drop is balanced with the subtree inflow exactly like the other devices (bisection on `m − Δp(A·(D0 − m)) = 0`).

**Acceleration (WELSEGS item 6 'HFA').** The axial momentum balance of a segment gives
`p_seg − p_out += (M_seg − Σ M_child) / A_seg`.
* Homogeneous flow: `M = ρq²/A`.
* Drift flux: `M = (ρ_g j_g²/α + ρ_l j_l²/(1 − α))·A`.
* Perforation inflow enters radially and carries no axial momentum.

The sign does not depend on the flow direction: injection decelerating into the reservoir gives a pressure recovery.

**Drift flux (WELSEGS item 7 'DF', WSEGDFPA).** This is the gas–liquid model of Shi et al. (2005, SPE J. 10(1)):
* `v_g = C0·j + v_d`;
* `C0 = A / (1 + (A − 1)γ²)`, `γ = clamp((β − B)/(1 − B))`, `β = max(α, F_v·α·|j| / v_sgf)`;
* `v_d = (1 − αC0)·C0·K(α)·V_c·m(θ) / (αC0·√(ρ_g/ρ_l) + 1 − αC0)`;
* `V_c = (σgΔρ/ρ_l²)^¼`;
* K = 1.53/C0 for α ≤ a1, Ku for α ≥ a2, and linear in between;
* Ku from the Wallis–Makkenchery bond-number fit (C_K = 142, C_w = 0.008; tends to 3.2 for large pipes);
* `v_sgf = Ku·√(ρ_l/ρ_g)·V_c`;
* inclination factor Hasan–Kabir `m = √cos θ·(1 + sin θ)^1.2` (WSEGDFMD default 'H-K'), negative for downflow and zero horizontally.

Defaults are A 1.2, B 0.3, F_v 1, a1 0.2, a2 0.4, σ = 0.03 N/m. The holdup solves `α(C0·j + v_d) − j_g = 0` by bisection on [0, 1].
It sets the hydrostatic density `α ρ_g + (1 − α) ρ_l` of co-current flow toward the outlet.
Friction stays homogeneous. Counter-current flow and injection are treated as homogeneous.

The simulator passes, per connection, the free-gas fraction of the reservoir inflow (mobility / B_g share) and the gas / liquid densities.

**Limits:**
* oil and water travel as one liquid;
* phase densities are the inflow (reservoir-condition) values, so there is no flash or expansion up the tubing;
* there is no segment storage.

**BHP reference depth (bug fix).** For a multisegment well, the BHP reference depth is now the top-segment depth (WELSEGS item 2).
This overrides WELSPECS item 5, as in OPM (`Well::updateSegments → updateRefDepth(depthTopSegment())`) and ECLIPSE.
Before the fix, a WELSPECS depth of 0 was read as "defaulted" and set to the top connection. For tubing modelled up to surface,
that gave a negative top-node pressure and an intermittently dead well.

**Verification (`-c36msw`, 9 checks):**
* Labyrinth Δp equals the hand formula (17187.5 Pa), and the in-well device balance agrees to 1e-6.
* Acceleration matches the hand momentum balance (0.069493 bar at the toe, 0.061078 bar for the contraction).
* Drift-flux limits:
  * bubbly upflow: α 0.00617 < λ 0.00990, C0 1.2, v_d 0.2044 ≈ 1.53·V_c;
  * horizontal: v_d = 0;
  * downflow: α 0.187 > λ 0.167;
  * annular-like: α 0.905 vs λ 0.909, C0 → 1.
* Vertical 1000 m column: holdup 0.165 vs λ 0.25 gives a head of 67.2 vs 61.3 bar.
* MSW_HORIZONTAL with WSEGLABY ICDs: inflow coefficient of variation 0.566 vs 1.219 for the open tubing; FWPT −2.1 %; FOPT +5.2 %.
* HFA on MSW_HORIZONTAL: max acceleration 1.95 psi (≈ ρv²); FOPT −0.02 %.
* SPE1 with the producer's 8400 ft tubing to surface (0.5 ft ID, top node at the 1000 psi limit):
  * DF holdup 0.864 vs no-slip 0.872;
  * final bottom-node pressure 2283 vs 2270 psi;
  * FOPT −0.06 %.
  * The effect is small because the flow is high-gas-fraction and fast (C0 → 1).

### 6.35 C36: network balancing and automatic compressors — NETBALAN, NETCOMPA, COMPOFF; gas-lift NONE caps (`network.cpp`, `gaslift.cpp`)

**NETBALAN** follows the OPM `Network::Balance` semantics.
* Item 1, interval: 0 balances every time step; > 0 balances at the first step of each interval, with node pressures and THP limits held in between; < 0 (NUPCOL) is treated as every step.
* Item 2, node-pressure tolerance (default 0.1 bar, converted for FIELD): replaces the internal relative tolerance of the implicit Illinois balance.
* Item 3, maximum iterations.
* Item 7, MAX_BALANCE_ERROR: an unconverged balance above it cuts the time step.
* Items 4, 5, 6 and 8 are read but not used.

The last balance time is part of restart state v11.

**NETCOMPA** (an ECLIPSE keyword; OPM parses it but does not simulate it). Items:
* 1 / 2: the inlet / outlet node, i.e. the branch;
* 3: group checked (default: the inlet node);
* 4: phase of the target (GAS / OIL / WAT / LIQ);
* 5: VFP table used while on (0 keeps the branch table);
* 6: ALQ;
* 7: fuel-gas rate;
* 8: extraction group;
* 9–11: ON-OFF or MULTI with n levels and the ALQ of level 1;
* 12: switching sequence.

At the start of each time step, before the balance, the compressors of the lowest sequence number whose group produces below 99.9 % of its GCONPROD target in that phase go up one level.
* MULTI: ALQ is linear from level 1 to item 6, and fuel = item 7 × level / n.
* Fuel gas is taken out of the network flow at the extraction group's node (like GCONSUMP) and reported in FGCR.
* COMPOFF switches all compressors off at its report step; they come back on while the target is still missed.
* Compressor levels are part of restart state v11.

**Gas-lift group caps with procedure NONE (bug fix).** A GCONPROD phase limit now caps the lift-gas allocation only if it is the control mode or its exceed procedure (item 7 / 11–13) is not NONE.
This follows OPM `GroupKeywordHandlers` (production_controls) and `GasLiftGroupInfo` (`has_control`).
GASLIFT-01 with PLAT-2 LRAT 3000 and RATE vs NONE: max FGLIR 1.504e6 vs 1.696e6, FGLIT 1.379e9 vs 1.400e9 (`-glift`).

**Verification (`-c36net`):** NETWORK.DATA with GCONPROD PLAT ORAT 6000 and a compressor table (suction 50 psi plus the same line loss). Results:
* The compressor switches on once, at 730.5 d, when PLAT drops below 6000 stb/d. FOPT is 1.1632e7 vs 1.0544e7 without it (+10.3 %).
* MULTI 3 levels with 900 Mscf/d fuel: FGCR steps 300 / 600 / 900.
* COMPOFF switches the compressor off and it switches back on.
* NETBALAN 365 d: the PLAT pressure changes 9 times in 85 steps (83 with every-step balancing).
* NETBALAN 0.01 psi / 20 iterations: FOPT within 0.01 % of the default tolerance.
* Restart with the compressor on is bit-identical (state v11).

Node pressures are now kept across report steps. This holds them between interval balances, and it is also a better first guess for the implicit balance (`-c30net` and `-facil` unchanged).


### 6.36 C37: multisegment wells — tubing flash, oil–water slip, segment storage, WSEGFLIM, WSEGPULL (`msw.cpp`, `sim_msw.cpp`)

**Tubing flash (on by default, `-- RESSIM: MSWOPTS flash storage`).** Before C37 the segments carried the reservoir
volumes of the inflow, so dissolved gas stayed in solution all the way up the tubing. Now each segment sums the
**surface** rates of its subtree (qo, qw, qg incl. dissolved / vaporised parts) and re-splits them at the segment
pressure with the PVT tables of the well's PVT region:
* `rs = min(qg/qo, Rs_sat(p))`, `rv = min(qo/qg, Rv_sat(p))`, and `rv = 0` if `rs·rv ≥ 1`;
* free oil `o' = (qo − rv·qg)/(1 − rs·rv)`, free gas `g' = (qg − rs·qo)/(1 − rs·rv)`;
* local rates `qo_loc = o'·Bo(p, rs)`, `qg_loc = g'·Bg(p)`, `qw_loc = qw·Bw(p)`;
* densities `ρo = (ρos + rs·ρgs)/Bo`, `ρg = (ρgs + rv·ρos)/Bg`, and viscosities from the same tables.

The segment solve iterates (up to 6 passes) until the pressures and the phase split agree. On SPE1 with 0.5 ft tubing to
surface the gas fraction at the top is 0.683 instead of ~0. The column drop is 1591.7 psi with flash and 2333.6 psi
without. FOPT goes down (4.925e7 vs 4.957e7 stb) because the expanding gas adds friction in this narrow tubing. That is
physics, not an error: in a wider tubing the lighter head wins.

**Oil–water slip (WSEGDFPA items 7–11).** Shi et al. (2005) give the oil-in-water holdup used inside the liquid phase:
* profile parameter `C0 = A/(1 + (A − 1)γ²)`, with `γ = clamp((αo − B1)/(B2 − B1), 0, 1)` (defaults A 1.2, B1 0.4, B2 0.7);
* drift velocity `vd = 1.53·Vc·(1 − αo)^N·m(θ)·f_ow`, where `Vc = (σ g Δρ/ρw²)^¼`, σ = 0.025 N/m and N defaults to 2;
* holdup `αo = jo/(C0·j + vd)`.

The result feeds the liquid density of the gas–liquid drift model (`ρ_L = αo ρo + (1 − αo) ρw`). It is active with the
DF model when item 7 > 0. In upflow the lighter oil overtakes the water: α_o is 0.0265 against a no-slip value of 0.0385.
The test water column has 0.154 oil holdup (no-slip 0.3) and a head density of 978 kg/m³ (no-slip 939).

**Segment storage (transient wellbore).** Each segment keeps its committed head density and gas holdup in the well state,
and these are saved in restarts. In a step of length Δt the content moves toward the steady value with the swept-volume
weight `w = |q|Δt / (V + |q|Δt)`. V is WELSEGS item 11, or else A·L. A shut segment keeps its content (a static column
after a shut-in). The scheme is an explicit mixing-tank update of the head term only. Friction and the devices stay
steady, so this is a storage model, not a full transient wellbore (no pressure waves, no slugging). On SPE1 the
cumulative effect is about 4e-9 relative. It matters for shut-ins / start-ups and large segment volumes.

**WSEGFLIM — segment flow limit (assumed semantics).** The ECLIPSE manual text was not available, so the semantics are
an assumption, built from the keyword items (OPM parses but does not simulate it). On segment range SEG1–SEG2 the inflow
of each limited phase (1 OIL, 2 WAT, 3 GAS, 4 LIQ) is held to FLOW_LIMIT. A smooth conductance factor multiplies the
connection factors of the segment:
* `f(Q0) = (1 + (Q0/Q_lim)^a)^(−1/a)`, with a ≥ 1 (default 8; item A);
* f is the product over the two limit pairs, with a floor of 1e-6;
* Q0 is the **unlimited** potential of the segment.

Q0 is explicit, one step lagged. The fixed point `Q0 = q_now/f(Q0_used)` is extrapolated one step with the clamped
growth `clamp(Q0_new/Q0_used, 0.8, 1.25)`, and the potential is kept in the well state for restarts. On SPE1 with OIL
8000 stb/d the median FOPR after 60 d is 8014.9 (0.2 %). Transients reach +10 % (max 8800) when the potential grows by
~20 %/step. A fully implicit limiter would need coupling into the well equations; that is not done. The limiter only
binds while the well is BHP-limited, because under rate control the BHP compensates.

**WSEGPULL — downhole pump / pressure boost (assumed semantics).** On SEG1–SEG2 a production pump adds
`Δp_boost = P·max(0, DP − G·max(0, qo_sub − Qo))`, where P is the on/off factor (item P), DP the boost (item DP), Qo the
oil rate (item OIL_RATE_LIMIT) above which the boost decays with GRADIENT_FACTOR G, and qo_sub the surface oil rate of the
subtree. It acts in the segment pressure balance, toward the top. On SPE1 a boost of 800 psi lowers the bottom node from
2237 to 1525 psi and raises FOPT from 4.925e7 to 5.057e7.

Tests: `-c37msw` 13/13. They cover the parser, the OW holdup limits, the limiter function, the flash / OW column /
pump + limiter / storage unit solves, SPE1 flash on / off, WSEGFLIM, WSEGPULL, storage on / off, restart from step 5
(bitwise) and the MSWOPTS export round trip. `-c36msw` now compares HO and DF with flash and storage off (C36 meaning).

### 6.37 C37: PYACTION — Python actions in the ResSim engine (`pyaction.cpp`, `sim_udq.cpp`)

`PYACTION` (OPM Flow keyword: record 1 NAME + RUN_COUNT SINGLE / UNLIMITED / FIRST_TRUE, record 2 the script file,
relative to the deck) is evaluated by the ResSim black-oil engine.

**Process model.** There is no embedded interpreter, so ResSim needs no Python headers or DLL at build time. One
persistent child process runs a generated bootstrap (`python -u ressim_pyaction_host_<pid>.py`). It is started with
fork + exec on POSIX and CreateProcess on Windows, and talks over pipes in line-based JSON. The interpreter comes from
`Case::pythonExe` (app setting), then `RESSIM_PYTHON`, then `python3` / `python` / `py -3`. A ping handshake returns the
version. Without Python the run continues with a warning and the actions are disabled.

**opm_embedded shim.** The bootstrap registers a module `opm_embedded` with:
* `current_summary_state`: `ss['FOPR']`, `in`, `get`, `well_var`, `group_var`, `has_well_var`, `has_group_var`,
  `wells(pattern)`, `groups()`, `elapsed()` in seconds, and `update`;
* `current_schedule`: `open_well` / `shut_well` / `stop_well` → WELOPEN, `insert_keywords(text)`, `well_names`,
  `get_well(name).status()`, `.start`;
* `current_report_step`, a minimal `current_ecl_state` (title, units), and `OpmLog`.

Summary values are fetched **on demand**: each lookup is a query message answered by the engine from the UDQ / ACTIONX
source (all F*, W*, G* vectors the UDQ engine knows, plus UDQ values), and cached for the call. Scripts run in one
persistent namespace per action:
* new style: the whole script is executed at each call, so `if 'setup_done' not in locals():` works;
* legacy style: the module body runs once, then `run(ecl_state, schedule, report_step, summary_state,
  actionx_callback)` is called. The bool return drives FIRST_TRUE, and `actionx_callback(name, wells)` applies the body
  of the deck's ACTIONX `name`, with `?` expanded to the wells.

`print()` is captured into the run log. A Python exception sends its traceback to the log, and that action is not called
again.

**Semantics.**
* When it runs: at the end of every report step, after the ACTIONX blocks. The keywords the scripts return are applied
  through the ACTIONX path (`actionApply`), so they act from the next report step. They are stored in `actStep_` /
  `actText_` and replayed on restart.
* Run counts: SINGLE = one call; FIRST_TRUE = until run() returns True; UNLIMITED = every call.
* Restart state (v12): call counts and finished flags are saved. **Not** restored: the Python globals.
* Differences to OPM Flow: OPM calls the scripts every time step (ResSim: report step); ecl_state is minimal; schedule
  editing goes through deck text (WELOPEN / insert_keywords), not OPM's object API.

Tests: `-c37py` 10/10:
* JSON round trip, record parsing, and host protocol (queries, globals, print, errors, clean stop);
* SPE1 with PYACTION insert_keywords is bitwise identical to the equivalent ACTIONX (FOPT 25 990 000.4351 both);
* the legacy run() + actionx_callback form and FIRST_TRUE;
* a broken script: logged, run unchanged;
* restart from step 20 equals the continuous run; export verbatim.

Example: `decks/pyaction/SPE1_PYACTION.DATA` + `spe1_rate_cut.py`.

### 6.38 C37: WCYCLE well cycling, GRIDUNIT, WTEST 'P' closures (`simulator.cpp`, `deck.cpp`)

**WCYCLE** `well on off start maxDt YES|NO` (days) is kept as a verbatim extra keyword (exact round trip). The records
in force are re-derived at each report step from every WCYCLE up to that step, including those added by ACTIONX or
PYACTION. A later record for the same name replaces the earlier one. The semantics follow OPM Flow's
`WCYCLE::nextTimeStep / wellStatus / efficiencyScale`:

* Each well has an *open time* t_o and a *close time* t_c.
  * t_o is set when the schedule opens the well: a status change, a new well, an injector/producer switch, or **any**
    well-level `WELOPEN OPEN`, even for a well that is already open (OPM `REQUEST_OPEN_WELL`).
  * A WTEST re-opening also sets t_o.
  * Both times are cleared while the schedule keeps the well shut.
* **Status at the start of every time step:**
  * open while t < t_o + on; then shut with t_c = t;
  * shut while t < t_c + off; then open again with t_o = t.
  * A cycle only switches wells that the schedule keeps open.
* **Step size:**
  * With `YES`, the step is cut so that it ends exactly on t_o + on or t_c + off.
  * The first step after a cycle re-opening is at most `maxDt`.
* **Start-up ramp:** for t − t_o < start after a *cycle* re-opening (not after a schedule opening), the well efficiency
  factor is f = f_WEFAC · min(1, (t − t_o + Δt)/start). It scales the reservoir source, the group rates and the cumulatives.
* **Restart:** t_o, t_c and the schedule-opening flag are saved in restart state v13.

**WTEST 'P'**: an open well that cannot flow at its limits is closed at the end of the step and re-tested every WTEST
interval. Examples are an injector whose BHP limit is below the reservoir pressure, and a producer whose BHP limit is
above it. This mirrors OPM, where non-operable wells are closed with reason PHYSICAL. Wells without WTEST 'P' keep the
earlier behaviour: they stay open at zero rate.

**GRIDUNIT** `METRES|FEET|CM [MAP]` scales COORD / ZCORN into the deck unit system on import. The factors are
toM = 1, 0.3048 and 0.01; the stored factor is toM in METRIC decks and toM/0.3048 in FIELD decks.

**opm-tests effect** (`-c37opm`):
* wcycle: 9 of 9 cases pass. Before C37 all 9 failed with FGIR 35–71 %.
* gridunit: the CM / FEET variants dropped from 80–99 % to the same 8.3–8.4 % as the METRES base.
  * That residual is time-step discretisation of the decline after the BHP switch. OPM's step 12.77 → 15 d puts the
    well on BHP control for the whole step, while ResSim's steps land nearer the switch at about 13.6 d.
  * The pressures before the switch agree to 0.01 bar.

### 6.39 C37: large-model performance (1 M-cell round)

**Benchmark generator:** `tools/gen_bigcase.py NX NY NZ out.DATA [--months M]`. It writes a synthetic SPE1-PVT deck with
layered sinusoidal log-permeability, 4 LRAT producers and 1 water injector. Measured on 2 cores and 2 GB RAM with
6 months of simulation:

| cells | wall | peak RSS | steps / Newton / linear | assemble / precond / solve (s) |
|---|---|---|---|---|
| 9 k | 2.8 s | 45 MB | 13 / 46 / 301 | 0.6 / 0.8 / 1.3 |
| 25 k | 5.9 s | 106 MB | 12 / 33 / 323 | 1.3 / 0.9 / 3.5 |
| 50 k | 11.9 s | 206 MB | 10 / 32 / 318 | 2.4 / 1.8 / 7.3 |
| 100 k (before C37) | 45.3 s | 454 MB | 12 / 47 / 624 | 7.0 / 6.7 / 30.6 |
| 100 k | 42.4 s | 405 MB | 12 / 47 / 624 | 6.7 / 5.0 / 29.7 |
| 100 k, 1 thread | 65.8 s | 395 MB | same | 15.0 / 6.5 / 42.9 |
| 200 k | 108.8 s | 801 MB | 12 / 62 / 817 | 17.2 / 13.1 / 76.5 |

**C37 changes.** The results are bit-identical: FOPT 1.82625e6 at 100 k cells.
* **BILU0:** in contiguous level-scheduled mode the factor no longer copies the matrix values into `LU_` (46 MB at 100 k).
  `LU_` keeps only the pattern.
* **WellSchur:** an FNV hash of the block pattern and of the well couplings lets the Schur matrix S be reused between
  Newton iterations. Only its values are refilled, and `finalizeDiag` is skipped. S.col is trimmed after a rebuild.
* **Grid::build:** trims the connection vector's growth slack (32 MB at 100 k) before the adjacency is built.

**Profile (gprof, 1 thread, 100 k cells):**
* BILU0 apply: 36 % of the time.
* Block SpMV: 16 %. Each SpMV takes about 8 ms per call, which is memory-bandwidth bound.
* Assembly: 17 %.
* AMG Gauss-Seidel: 7 %.

The solve time per linear iteration is linear in n, at 0.47 µs per cell (4.4 ms at 9 k cells, 93.6 ms at 200 k).
The number of linear iterations per Newton step rises slowly, from 6.5 to 13.

**Float32 BILU0 in CPR** (`-- RESSIM: TUNING … mixedPrecision=1`) does **not** pay off here. At 100 k cells it took
46.9 s against 42.4 s with 627 iterations, because the float factor has no level-parallel apply. It stays off by default.

**Extrapolation to 1 M cells (same physics, 6 months):**
* Memory: about 4.0 GB (about 4.0 kB per cell; the scalar CSR copy is 0.8 kB per cell of it).
* Time on 2 cores: about 80 Newton and 1,400 linear iterations, giving roughly 11 min of solve and 14 min in total.
* On an 8-core workstation with about 3× the memory bandwidth: about 5 min.
* The 1 M-cell run itself cannot be done on this 2 GB build machine.

**Next steps:**
* ~~Drop the scalar CSR Jacobian when the block path is active (−0.8 kB per cell).~~ Done in C45 (block assembly, §9.v).
* Add a level-parallel float BILU0 and a fused SpMV + dot. C45 adds the lean float32 factorisation (§9.v); a fused SpMV + dot is still open.

### 6.40 C38: production-data import and history generation (`proddata.cpp`)

**Input.** A table with a header row and one row per well and date (CSV, TSV, semicolon or whitespace). Columns are
recognised by name (case-insensitive aliases, unit suffixes allowed: `OIL`, `ORAT`, `WOPRH`, `Oil (sm3)`, `GAS_KSM3`,
`BHP_bar`, `DAYS`, `HOURS`, …) or mapped by hand in the UI. Dates: ISO, `DD-MON-YYYY`, `D MON YYYY`, `DD/MM/YYYY`
(or month first), dotted, `YYYY-MM`, `MON-YYYY`, ISO with a time part, or YEAR / MONTH / DAY columns.

**Periods.** The sorted unique dates `d_0 < d_1 < … < d_{n-1}` define the history periods. With "dates open the period"
(default), period i is `[d_i, d_{i+1})`; the last one is closed one calendar month later for monthly data (median date
spacing 28–31 d, month-end dates stay month-end), otherwise by the median spacing. With "dates close the period",
period i is `[d_{i-1}, d_i)`.

**Rates.** For volumes `V` over a period of length `L` and uptime `u ≤ L` (DAYS, or HOURS/24):

* calendar-average rate (no uptime, or WEFAC off): `q = V / L`, `WEFAC = 1`;
* with uptime and WEFAC on: flowing rate `q = V / u` and `WEFAC = u / L`.

Both give the same cumulative `q · WEFAC · L = V`, which the engine honours (C33 WEFAC: rates, cumulatives and the
reservoir source use factor × rate). Unit factors scale the liquid and gas columns.

**QC** (severity info / warning / error):

* unparseable date or missing well → error, row dropped (a second header row of units is recognised and skipped);
* negative value → 0 (warning); non-numeric cell → missing; BHP / THP ≤ 0 → missing;
* duplicate well and date → volumes summed (rates averaged), warning;
* uptime longer than the period → clamped; volumes with zero days on → calendar rate;
* producing and injecting in one period → injection used; water and gas injection together → water used;
* jumps of the water-oil ratio or GOR by more than a factor (default 10) against the previous producing period →
  outlier warning;
* well names longer than 8 characters (ECLIPSE limit), wells that never flow, periods without data after the first
  record (shut there), gaps of more than 3× the typical spacing;
* at merge time: wells not defined in the deck (WELSPECS) and deck wells without data.

**Keywords per period.** Producers get `WCONHIST 'W' OPEN <ctrl> qo qw qg 1* 1* thp bhp /` (control RESV by default, or
ORAT / WRAT / GRAT / LRAT / BHP; measured THP / BHP as observed values), shut periods `SHUT`; injectors
`WCONINJH 'W' WATER|GAS OPEN q bhp thp /`; `WEFAC` only where the factor changes against the well's previous open
period.

**Merge into a deck (text level).** The deck's SCHEDULE is split into time blocks at its DATES and TSTEP (TSTEP values
with `n*v` repeats; DATES with an optional `HH:MM:SS`). A history period starting at the same time as a block is appended
to that block (after its WELSPECS / COMPDAT); otherwise a new block is inserted. The schedule is extended to the end of
the last period. Boundaries are written as DATES for whole days (consecutive empty boundaries grouped into one DATES
keyword) and as TSTEP otherwise; END and everything after it is kept. Optionally, the WCONPROD / WCONHIST / WCONINJE /
WCONINJH records of the imported wells are removed, so the history alone controls them (wildcards are left alone).
The merge refuses decks whose dates do not increase, warns about SCHEDULE INCLUDE files (their dates are invisible to
it) and uses the case START when the text has none. A stand-alone include (`DATES` + keywords, no DATES at START) and an
observed-data CSV (calendar-average rates at the period ends, for history matching) can be written instead.

**Verification (`-c38prod`).** SPE1 with a 4-month table under ORAT control, uptime 14 of 28 days in February
(WEFAC 0.5) and a gas injector under WCONINJH: FOPT(120 d) = 218 999.96 against 219 000 in the table, and FGIT = 4.5·10⁶
exactly. Mid-month dates insert new DATES, and the schedule is extended to 104 d with FOPT = 90 000 as in the table.

### 6.41 C38: float32 factors in the level-parallel BILU0, NONNC, UDT, Newton minimum

**Float factors (`bsr.cpp`).** The contiguous (level-ordered) BILU0 path of C28 / C37 factorises in double precision.
With the float option, the factors are then copied to a float32 array and the double array is freed until the next
factorisation. The triangular sweeps read float blocks and accumulate in double (`z_i -= Σ float(L_ik) · z_k`), so only
the factor entries are rounded (relative 6·10⁻⁸). The steady-state factor storage halves, and so does the memory traffic
of the sweeps. Peak memory during the factorisation is unchanged (the double array exists while factorising).
`TUNING … mixedPrecision = 1` (`-- RESSIM: TUNING 0.3 0 2 1 1`) now selects this path for the CPR second stage and the
stand-alone BILU0. The serial BILU0f (C16) is only used when the contiguous path is off (small or thin models).

* Unit test (`-c38ilu`, 7 200 cells × 3, 7-point): factors 3.43 MB → 1.72 MB, apply 0.58 → 0.50 ms, |Δz|/|z| = 3.4·10⁻⁸.
* 100 k-cell benchmark (gen_bigcase, 6 months, 2 threads): 58.0 s → 50.8 s (−12 %), identical time steps, Newton (47)
  and linear (624) iterations, identical FOPT.

**NONNC (`grid.cpp`).** Removes the X / Y connections whose two cells are not Cartesian neighbours (fault-throw
connections across displaced layers). Kept: Cartesian neighbours, Z / PINCH connections, explicit NNC, LGR and
dual-porosity connections. On FAULTED.DATA, 176 of 10 916 connections are dropped. opm-tests editnnc/NONNC:
13.7 % → 3.05 % (pass).

**UDT (`udq.cpp`).** A UDT table `TU_NAME` (1-D) is used as `TU_NAME[expr]` in UDQ expressions, elementwise for well /
group sets. Interpolation: `LC` linear with constant extrapolation, `LL` linear with linear extrapolation, `NV` nearest
point. Union operators `UADD USUB UMUL UDIV UPOW UMIN UMAX UEQ UNE UGT UGE ULT ULE` have the lowest precedence; the result
is defined where either operand is (an undefined operand gives the other one).

The evaluation order now follows OPM (`UDQConfig`): the due ASSIGNs first, then the latest DEFINE of every name – well,
group, then field quantities – in the order of the name's first appearance. So `ASSIGN X 0 / DEFINE X0 X / DEFINE X (…)
UMIN X0` evaluates X before X0, sees an undefined X0, and yields the table value (opm-tests UDT-1D: 180, then 153.667,
127.006, …; ours 180, 153.675, 126.957).

**Newton minimum and TSFMAX.** At least one Newton update per time step (OPM `NewtonMinIterations`). Without it, a
near-steady state whose residual is below the tolerances at iteration 0 kept the previous well rates (UDT-1D: rates at
the BHP limit stopped declining at 0.016 m³/d instead of decaying to 3·10⁻⁴). Opt-out: `RESSIM_NEWTON_MIN0=1`. TUNING
TSFMAX now defaults to 3 (ECLIPSE / OPM) when TUNING is given; before, ResSim kept its own 2. With both, the UDT-1D time
steps coincide with OPM's (0.1, 0.4, 1, 2, 3 … d).

**EDIT-section TRANX / TRANY / TRANZ (`deck.cpp`, `grid.cpp`).** `MULTIPLY 'TRANX' f box` in EDIT used to be folded
into MULTX. Now it is stored separately (`GridInput::tranMultX/Y/Z`) and multiplies only the Cartesian-neighbour
connection (i,j,k)–(i+1,j,k): for corner-point lateral faces only when k₁ = k₂. Fault-juxtaposition connections across
the same face keep their transmissibility, because OPM writes `TRANX[c1]` only for `gc2 − gc1 = 1`
(`Transmissibility::createTransmissibilityArrays_`). MULTX still applies to every connection across the face. The
exporter writes an EDIT section `MULTIPLY` with runs along I. model2 8_MULTIPLY_TRANXYZ: 5.94 → 5.03 %. Unit test
`-c38tran` on FAULTED: 3 800 neighbours halved, 96 fault NNCs unchanged, round trip exact.

EDITNNC keeps the ECLIPSE rule (non-neighbours only). OPM would also scale a Cartesian neighbour, but every entry in
model2 9_EDITNNC is a fault juxtaposition, so the result is the same either way. Its remaining 12.5 % (late WWPR) is
the ~4 % base-model offset seen in 0_BASE_MODEL2, amplified at water breakthrough; it is not EDITNNC (374 / 374
entries matched).

**SWATINIT below SWL (`simulator.cpp`).** As in OPM `EclMaterialLawManager::applySwatinit`: S_w = max(SWATINIT, SWL_cell),
SWL unchanged, P_cow scaled by P_c,eq / P_cow(S_w). Before, the cell started below connate. model2 5_SWATINIT:
4.61 → 1.65 %.

**Production data: end dates.** With "the date marks the end", the date is inclusive: 31 JAN closes January, and the
period is [1 JAN, 1 FEB 00:00). Monthly end-of-month data therefore land on the deck's monthly report dates. The merge
warns when the largest producing GOR is implausible for the deck units (> 100 Mscf/stb FIELD, > 18 000 sm³/sm³ METRIC),
the typical scf-vs-Mscf mistake. End-to-end on SPE1 (3 months of table data): ORAT history reproduces FOPT
620 000 / 1 160 000 / 1 760 000 stb exactly, injection FGIT exact; RESV history gives +1.5 … +6.5 % oil (voidage
matching with the model's PVT, expected).

**opm-tests effect (model2 + editnnc + udt, 59 cases, full key list).** C34: 12 pass / 41 fail / 3 error / 3 skip.
C38 (Newton minimum, TSFMAX, NONNC, UDT): 13 pass / 43 fail / 0 error / 3 skip; three cases improved by > 20 %, none got
worse by > 25 %. The three former UDT errors now run; their residual difference is at noise-level rates.

### 6.42 C40: geomodel data analysis — log calculator, regression, variogram fitting, well tie (`geo_analysis.cpp`)

**Log calculator.** A recursive-descent parser builds an expression tree once; the tree is evaluated per sample (vectors of
equal length; curves, `MD` and `TVD = at(md).z` as variables). Grammar, lowest precedence first:
`or := and {'||' and}`, `and := cmp {'&&' cmp}`, `cmp := add {(<|<=|>|>=|==|!=) add}`, `add := mul {(+|-) mul}`,
`mul := unary {(*|/) unary}`, `unary := (-|+|!) unary | pow`, `pow := atom ['^' unary]` (right-associative, so `-2^2 = -4` and
`2^3^2 = 512`). Null samples (LAS null value) become NaN and propagate through arithmetic, functions, comparisons and `if`;
`&&` / `||` return a definite result only when one operand decides it (`0 && NaN = 0`, `1 || NaN = 1`); division by zero, `sqrt` of
a negative number and `ln` / `log10` of a non-positive number give NaN. Results are written back with the log's null value.

**Regression (cross-plots).** Ordinary least squares on finite pairs, y = a + b x, or log10 y = a + b x with y ≤ 0 excluded;
R² = 1 − SSE / SST. The log-linear form is the poro-perm transform of the Petrophysics step (`PoroPerm{a, b, r2, n}`).

**Variogram fitting.** For γ(h) = c0 + c · f(h / a) with the unit structures of `geostat` (spherical; exponential
1 − e^{−3h/a} and Gaussian 1 − e^{−3h²/a²} with practical range a), the fit minimises Σ w_k (γ̂_k − γ(h_k))². For a fixed range
the problem is linear in (c0, c): the 2 × 2 normal equations are solved and projected onto c0, c ≥ 0 (c0 = 0 → c from the
one-parameter solution; c = 0 → c0 = weighted mean). The range is found by an 80-point logarithmic scan over
[0.05 h_max, 2 h_max] followed by a golden-section refinement in ±15 %. Weights start as N_k / h_k² and are then iterated three
times as Cressie's N_k / γ(h_k)² (Cressie 1985, *Math. Geol.* 17:563), which favours short lags with many pairs. Verification: the
three models are recovered to < 1e-7 from exact points; an AR(1) field with practical range 40 gives a fitted range of 46 (sampling
error of 30 × 400 samples). Applying a fit to a geomodel step stores c0 / (c0 + c) and c / (c0 + c), because the SIS / SGS
simulations work with unit total sill.

**Well tie.** With slowness s(z) from DT (µs/ft or µs/m; automatic: median < 300 → µs/ft) converted to s per project length unit,
two-way time is integrated from the first valid sample, t(z_0) = t_vm(z_0) of the current velocity model and
t_i = t_{i−1} + 2000 · (z_i − z_{i−1}) · ½ (s_i + s_{i−1}) ms (z = TVD from the trajectory, non-increasing depths skipped).
Slowness and density are resampled linearly onto a regular time axis (default 1 ms); AI = ρ / s (g/cc · m/s);
r_k = (AI_{k+1} − AI_k) / (AI_{k+1} + AI_k); the synthetic is r * w with the Ricker wavelet w(t) = (1 − 2π²f²t²) e^{−π²f²t²}
truncated at ±1.5 / f. The seismic trace is taken at the nearest trace to the well position at mid-log MD. The bulk shift Δ
maximises the Pearson correlation between the synthetic, linearly interpolated at trace times t − Δ, and the trace over their
overlap (scan in 0.5 ms steps within ±Δ_max). The tied time–depth pairs (z_i, t_i + Δ) calibrate the linear-gradient model
v(z) = v0 + k (z − z_d), t(z) = (2000 / k) ln(1 + k (z − z_d) / v0): for k ∈ [0, 2] in 61 steps, v0 is found by golden section (the
misfit is unimodal in v0), and the pair with the smallest RMS time misfit is returned. Over a short logged interval v0 and k
trade off (only the interval velocity is constrained), so the calibrated model is meant for the depth range of the logs.
Verification on the demo (its DT / RHOB carry the same impedance steps as the synthetic SEG-Y and the background velocity
1800 + 0.6 z m/s): correlation 0.95 at all six wells with |Δ| ≤ 2 ms; starting from a wrong model (k = 0.65, ~30 ms early) the
shift is found (+30 ms) and the calibrated model predicts the true TWT within 0.3 ms.

**Geomodel parameters in Studies (`study.cpp`, `geostat::applyPropModel`).** A case stores property models
m = (property, method, γ(h; c0, c, a_max, a_min, a_v, θ), seed s, hard data {x_i, y_i, z_i, v_i}, optional a + b·φ transform).
For a run with parameter values v, `applyParams` first builds the modified recipes — *GeoRange*: (a_max, a_min, a_v) ← v·(a_max,
a_min, a_v); *GeoSeed*: s ← max(1, round(v)) — and re-runs kriging / SGS (§6.21 *Geostatistics*, the same code as the Model panel) for
every model whose recipe changed; then PERMX = PERMY = clamp(10^(a + b φ), 10⁻⁴, 10⁵) mD and PERMZ = (kv/kh)·PERMX on porosity
models with a transform; only then the ordinary parameters (multipliers, contacts, …) act on the new fields. Unchanged
recipes leave the case bit-identical, so the base run of a study equals the base case. The seed is a discrete variable: in
sensitivity it gives realisation spread, in ensemble / ES-MDA runs it samples geological uncertainty; gradient-based history
matching should keep it fixed (its finite-difference derivative is meaningless). Verification (`-c40geo`, SPE1-like grid with
5 hard data): equal seeds give identical fields, other seeds differ, the data cells keep their values, PERM matches the
transform to 10⁻⁶; averaged over 6 seeds the lag-1 correlation along I rises from 0.05 (ranges × 0.25) to 0.72 (× 4); the
models survive the project text round trip exactly.

## 7. Linear algebra (`linsolve.cpp`, `amg.cpp`, `bsr.cpp`)

- **Krylov:** restarted flexible GMRES and BiCGSTAB on an abstract `LinearOperator` (CSR, BSR, Schur-reduced).
  OpenMP-parallel kernels.
- **WellSchur:** exact Schur complement of the BHP unknowns, S = A_cc − A_cw D_ww⁻¹ A_wc. The BHPs are recovered
  afterwards.
- **BlockCPR (two-stage CPR):**
  1. Decoupling: true-IMPES (accumulation-based weights, the default) or quasi-IMPES.
  2. AMG V-cycle on the pressure system.
  3. BILU0 on the full block system.
- **AMG:** plain (unsmoothed) aggregation by default (OPM Flow style, robust on anisotropic thin layers), with smoothed
  aggregation optional. Galerkin coarse operators and an exact coarsest LU. `refresh()` reuses the aggregates and
  transfer operators across Newton iterations.
- **BILU0:** block ILU(0) with exactly inverted diagonal blocks (pivoted LU). It is robust for cells with a vanishing
  phase, where a scalar ILU0 breaks down.
- **C12 profiling (Egg model, 18,553 cells, 2 threads).** `rundeck` prints the Newton-loop profile, and `RESSIM_DUMP` writes
  a Jacobian that the stand-alone benchmark can read:

  | Component | Time |
  |---|---|
  | Total run | 29.1 s |
  | Assembly | 6.0 s |
  | CPR setup | 6.7 s |
  | Krylov | 16.1 s (3,850 BiCGSTAB iterations, ≈ 4.2 ms each) |
  | Block SpMV | 0.48 ms |
  | AMG V-cycle | 0.68 ms |
  | BILU0 apply | 0.56 ms |
  | Whole CPR apply | 1.8 ms |

  - BiCGSTAB costs two SpMVs and two preconditioner applications per iteration. FGMRES needs 14 iterations instead of 11
    for the same time (one application per iteration, plus orthogonalisation), so the solver is bound by memory bandwidth,
    not by the Krylov method.
  - C12 change: the second CPR stage forms r − A·z₁ from only the pressure column of each block, since z₁ = (x_p, 0, …).
    This costs 1/n_b of a block SpMV and is exact. Egg: 29.1 → 28.1 s with identical results.
  - Next steps:
    - a single-precision (float) copy of the BILU0 / AMG data, which halves memory traffic;
    - multicolour Gauss–Seidel and BILU0, so the triangular sweeps can be threaded;
    - assembly reuse of the flux stencils.
  - JutulDarcy reports 8.6 s for Egg on its authors' workstation. A like-for-like ratio needs both codes on the same
    machine; Julia / JutulDarcy is not installed here.
- **C26 AMG hierarchy reuse across time steps** (`TUNING` field 7 `amgReuse`, default 10; 0 = rebuild every time step). Within a
  Newton loop the aggregates and transfer operators were already reused (`AMG::refresh`, values only). C26 also keeps them at
  the first Newton iteration of the next time steps, for up to N steps **within one report step** (a fresh build at every
  report step keeps restarts and iterative-coupling rewinds bitwise identical to the continuous run — reuse across report
  steps broke the 1e-10 restart checks). A quality monitor records the BiCGSTAB iterations of the last fresh build; a reused
  solve that fails or needs more than 1.5 × base + 3 iterations marks the hierarchy stale and the next step rebuilds it.
  Egg (CPR forced): 29.65 → 28.78 s (−3 %; −9 % when reuse also crossed report steps), identical results (FOPT 558,594);
  SPE9-like: no gain. **Honest result:** a small win — the bulk of the CPR set-up is the BILU0 factorisation, not the AMG
  aggregation.

---

## 8. Nonlinear solution and time stepping

- Fully implicit Newton–Raphson in all engines.
- **Black oil:**
  - Appleyard-type chopping of the saturation and pressure updates (dS_max, dp_max,rel 0.3).
  - Three-state variable switching.
  - Oscillation detection with relaxation.
  - In-Newton well control switching.
  - Relaxed CNV / MB tolerances on the last iterations (OPM practice).
- **Compositional:** inexact Newton (Eisenstat–Walker-style fixed tolerance 1e-5), early divergence cut, Aziz–Settari
  target-change control.
- **Black-oil target-change control (C27):** `-- RESSIM: TUNING` items 8 (dsTarget) and 9 (dpTarget). After an accepted
  step with largest changes ΔS and Δp, the next step is dt_next = dt · min_k (1+ω) η_k / (Δ_k + ω η_k) with ω = 1
  (Aziz & Settari NORM rule). It is capped by the growth factor and the Newton-count rule. *Verified (-c27):* small case
  with dS 0.02: 228 vs 53 steps; FOPT within 0.9 % of the uncontrolled run (truncation error of the coarse steps); last
  ΔS bounded.
- **Thermal:** see section 6.
- **Local solves:** safeguarded Newton (`localNewton`), regula falsi / Illinois (energy flash), and Rachford–Rice with
  bounded Newton.

---

## 9. Wells and field management

- **Standard well model:** connection rate q = WI · λ · (Φ_cell − p_wb,conn). Peaceman WI with directional
  completions; wellbore head lagged one step.
- **Controls:** WCONPROD, WCONINJE, WCONHIST, WCONINJH, WELTARG, WELOPEN. RESV at field-average conditions. A zero rate
  stops the well.
- **Groups:** GRUPTREE, GCONPROD, GCONINJE, GCONPRI (as limits). Shares are allocated by potentials with water filling.
  GUIDERAT formulas exist for compositional groups.
- **VFP / THP:** 5-D multilinear interpolation of VFPPROD/VFPINJ. The THP limit is the IPR–VFP intersection found by
  bisection.
- **Economic limits and well tests (C27):** WECON and WTEST stay verbatim keywords (exact round trip) and are read at
  their report step.
  - After every accepted step, each open producer with a limit is checked against: min oil / gas / liquid rate,
    max water cut, max GOR, max WGR.
  - Workover WELL shuts the well. NONE only reports the violation. Item 8 YES ends the run.
  - WTEST with reason 'E' re-opens the well every *interval* days, at most *n* times. The well is shut again if the
    limit is still violated.
  - The well state (economic flag, next test time, test count, stop flag) is part of the restart state (v4).
  - *Verified (-c27):* shut-in at the limit; NONE run identical; 3 re-openings and 4 shut-ins with WTEST; restart
    bitwise; deck round trip.

### 9.1 Passive tracers (C27, `tracer.cpp`)

**Input.** TRACER (name, phase WAT / OIL / GAS), WTRACER (well pattern, tracer, injected concentration) and TBLKF<name>
(per-cell initial concentration) stay verbatim, so OPM Flow runs the same deck.
`-- RESSIM: TRACERPROP name t½ k_d c₀` adds:
- first-order decay λ = ln 2 / t½;
- linear equilibrium adsorption: adsorbed amount per pore volume = k_d · c;
- a uniform initial value c₀.

**Equation.** The concentration c is the tracer amount per surface volume of the carrier phase α. After every converged
flow step, one linear equation per tracer is solved on the converged state (sequential transport, as in the ECLIPSE /
OPM tracer models):

  [C_i^{n+1} c_i^{n+1} − C_i^n c_i^n] / Δt + Σ_j F_ij^+ c_i^{n+1} − Σ_j F_ij^− c_j^{n+1} + q_i^prod c_i^{n+1} − q_i^inj c^inj
      + (e^{λΔt} − 1)/Δt · C_i^{n+1} c_i^{n+1} = 0,

where:
- C = PV b_α S_α + k_d PV is the capacity (free phase + adsorbed);
- F_ij = λ_α,up T_ij ΔΦ_α,ij is the upwind surface flux of phase α (threshold pressures included);
- the well rates are the converged connection rates of phase α.

The decay factor is the integrating-factor form of λC c. It is exact for decay alone, whatever λΔt is (short-lived
radionuclide tracers), and it only adds to the diagonal.

**Properties.** The matrix is an M-matrix (upwind, implicit), so 0 ≤ c ≤ max(c_inj, c₀) and no CFL limit applies.
- The equation is solved with ILU0-BiCGSTAB (GMRES fallback) on the polymer-transport pattern.
- Concentrations below 1e-30 of the scale are set to zero.
- The free-phase fluxes are antisymmetric, so the discrete balance is exact up to the solver tolerance.
- Well tracer rates use the component-vector convention (+ produced / − injected, amount/s). Field and well
  cumulatives come from the component accumulation.
- Snapshots carry `TR_<name>`.
- Restart state: concentrations and the decayed amount. The capacity is recomputed from the loaded state.

**Not modelled:** partitioning between phases (SOR tracer tests), dispersion / diffusion (beyond the upwind numerical
dispersion), parent → daughter chains (TOUGH3 EOS7R).

**Verified (-c27):**
- Balance error 6e-13 relative; max c − 1 = 7e-9 (solver tolerance); the flow solution is unchanged.
- k_d 0.5 delays breakthrough (1141 vs 1096 d).
- Decay without flow: M/M₀ = exp(−λt) to 1e-10 with steps of 0.25 half-lives.
- Restart bitwise; deck round trip.

---

### 9.x C41 additions: cross-flow, rock-compaction hysteresis, connection economics, facies, fault networks

**Cross-flow (WELSPECS item 10, default YES).** *Producer*: a connection with p_cell < p_conn injects the wellbore mixture
(surface fractions x_c of the well's gross production at the last assembly): q_t = λ_t WI (p_cell − p_conn), λ_t = Σ_p k_rp/(μ_p B_p),
surface component rates q_c = x_c q_t / V with V = x_w/b_w + (x_o − R_v x_g)/(d b_o) + (x_g − R_s x_o)/(d b_g), d = 1 − R_s R_v
(OPM StandardWell::computePerfRate). *Injector*: a connection with p_cell > p_conn produces at the cell mobilities,
q_p = WI λ_p b_p (p_cell − p_conn) (+ R_s, R_v terms); the rate control counts the net injected-phase surface rate. The cross-flow
terms carry no adjoint derivatives yet, and multisegment wells keep the previous no-cross-flow behaviour by default (C45: the
per-segment mixture exists as an experimental option, §9.v). The mixture is part of the restart state (version 14).

**ROCKCOMP HYSTER / ROCKTABH.** Each ROCKTABH region table is a set of elastic curves; the first row of every curve lies on the
deflation curve, which is used for p ≤ p_min (lowest pressure reached, updated after each converged step). For p > p_min the
multiplier is OPM's guided bilinear interpolation (UniformXTabulated2DFunction, LeftExtreme): with curve start pressures x_i ≤ p_min < x_{i+1},
α = (p_min − x_i)/(x_{i+1} − x_i), Δ = x_{i+1} − x_i, m = (1 − α) f_i(p − αΔ) + α f_{i+1}(p + (1 − α)Δ), continuous with the deflation
curve at p = p_min; linear extrapolation beyond the table ends. OVERBURD and ROCKOPTS STORE shift p and p_min alike.

**CECON.** After the WECON checks of every converged step, each open completion (COMPLUMP lump or single connection) of a producer
with limits is tested on its summed surface rates: water cut q_w/(q_o + q_w) (1 for q_o < 0), GOR q_g/q_o (∞ without oil),
WGR q_w/q_g (∞ without gas). CON closes the completion, +CON also every completion after it in the connection order, WELL shuts the
well; a well without open connections is shut (WTEST 'E' re-tests it). Stopped wells are skipped unless item 10 = YES. C42 adds
two items that OPM ignores (see §9.y).

**Facies by truncated Gaussian simulation.** One standard-normal SGS field Y (variogram of the recipe) is truncated at
t_k = Φ⁻¹(Σ_{j≤k} p_j); facies f = k where t_{k−1} < Y ≤ t_k. Hard facies codes at wells are conditioned at the normal score of the
middle of their class, Φ⁻¹((P_{k−1} + P_k)/2). Each cell gets φ_f, k_f (PERMY = PERMX, PERMZ = (k_v/k_h) k_f) and optionally
SATNUM = f. Proportions are honoured in expectation (8-realisation mean 0.533 / 0.288 / 0.179 for targets 0.5 / 0.3 / 0.2).

**Fault networks.** Junctions: a fault tip extended by the snap tolerance that meets another fault's trace (Y) or two traces that
cross (X). Truncation: on nine depth levels between top and base the minor trace is cut at its intersection with the major and a
terminal stick is inserted on the major trace. Throw profile: horizon depths sampled on both sides of the trace; hanging wall on
the dip side; throw = z_HW − z_FW; negative values beyond 5 % of the maximum are flagged as reverse; the profile also reports whether
the throw tapers towards the fault tips.

**Regression metric (material worst).** For liquid-phase W/F/G {O, W, L}P{R, T} vectors the deviation is scaled by
max(own scale, max over time of the reference oil + water series of the same entity); all-zero vectors are skipped. It is reported
next to the strict worst deviation and does not change the pass / fail status.

### 9.y C42 additions: geomodel import, 2-D seismic, calculators, CECON extensions

**Delimited-table sniffing (`geo_import.cpp`).**
- Delimiter: tab, ';', ',' and '|' are tried in that order on the first 30 body lines. The first one whose most frequent non-zero
  count of occurrences outside quotes covers at least 80 % of the lines, with at most 20 % of lines lacking it, wins. Otherwise
  runs of blanks are used.
- Decimal comma (only with a delimiter other than ','): cells of the form `d,d` are counted against numeric cells with a '.';
  a decimal comma is assumed when the first count is non-zero and at least as large.
- Header: row 0 is a header when some column is non-numeric in row 0 and numeric in row 1, or when rows 0 and 1 are both fully
  non-numeric and row 2 has numbers. In that case row 1 is a units row; it is stored and shown, not interpreted.
- Petrel `BEGIN HEADER … END HEADER` blocks give the column names. `# KEY: value` lines go to the metadata (well name).

**Column roles.** Each column has a role and a factor. The guess is in two steps:
1. Normalised header names (upper case, alphanumerics only) are matched against synonym lists (e.g. EASTING / XCOORD → X,
   KB / RKB → KB, DEPT / MD / MEASDEPTH → MD).
2. Columns that are still unassigned follow the canonical order of the kind (name, X, Y, KB, TD for heads; MD, INC, AZI for surveys).

The factor is chosen per column in the wizard: 1, ft → m (0.3048) or m → ft. `validateMap` reports what a builder still needs.

**Surveys.**
- MD / INC / AZI: minimum curvature from the head and KB.
- Absolute X / Y / TVDSS (or Z elevation, or TVD below KB): the stations are taken as given. When there is no MD column,
  MD is the cumulative 3-D arc length.
- Offsets DX / DY: added to the head.

Inclination and azimuth for display come from the station differences: inc = atan2(√(Δx² + Δy²), Δz), azi = atan2(Δx, Δy).

**Well matching.** Keys are compared case-insensitively. If no key matches exactly, the longest project well name that
prefixes the file name at a separator ('_', '-', '.', blank) is used. So `A1_survey.txt` → A1 and never A10.

**2-D lines.**
- *Geometry*: 2-D or 3-D is decided by the number of distinct inline / crossline values. More than one of each, forming a lattice,
  means 3-D (at most 4 traces per inline × crossline product and more than half of the traces on distinct pairs); otherwise the file is read as a line. The plan position comes from the X / Y bytes scaled by the
  coordinate scalar (byte 71; negative = divide), with source X / Y as the fallback.
- *Tracking*: from every seed the pick moves trace by trace to the extremum (peak, trough or zero crossing) within ±window of the
  previous pick, refined to sub-sample time by a parabola through the extremum sample and its two neighbours. Tracking stops when the
  jump exceeds maxJump or the event amplitude falls below 5 % of the seed amplitude.
- *Intersections*: segment–segment tests between consecutive traces of every pair of lines. The mis-tie is the difference
  between the picks of the two lines at the crossing traces.
- *Gridding*: picks are converted to depth through the velocity model and kriged (or inverse-distance weighted) onto the horizon
  lattice, or onto a padded bounding lattice when no horizon exists.

**Surface calculator.** Every surface is sampled bilinearly at the nodes of the target lattice (NaN outside its extent). The
expression is evaluated node by node with the log-calculator parser: operators, comparisons, `if`, min / max, `isnan`, with NaN
propagation. Variables are the surface names in upper case with non-alphanumerics replaced by '_' (a leading digit gets a '_'
prefix), plus X and Y.

**Property calculator.** The same parser runs over the cells. Variables: PORO, PERMX, PERMZ, NTG, FACIES, ZONE, BLOCK, DEPTH (cell
centre), TOP, THICK, BULK, I, J, K (1-based), X, Y and the calculated properties.

Writing the result:
- A new property is stored as is.
- Overwriting changes only the cells with a finite result.
- Permeability overwrites keep kv/kh: k_z ← k_z · k_x,new / k_x,old.
- Facies results are rounded to integer codes.

**Statistics.** `describe()` gives n, nulls, min, max, mean, the sample standard deviation (n − 1) and percentiles P10 / P50 / P90.
The percentiles use linear interpolation between order statistics at rank p(n − 1). These are percentiles, not the reserve
convention, so P90 is the high value.

**CECON items 11–12 and PLUG (beyond OPM).**
- A producing completion whose summed oil (gas) surface rate is below item 11 (12) violates its limit. Only that completion is
  closed (CON), whatever item 9 says.
- PLUG is a plug-back by depth. A water violation (WCT, WGR) closes the offending completion and every completion at the same or
  greater depth. A gas violation (GOR) closes it and every completion at the same or shallower depth.
- The WPLUG lengths are not used.
- This is a ResSim convention and has not been verified against ECLIPSE.

### 9.z C43 additions: geomodel tools (`geo_tools.cpp`) and well performance events

**Mis-tie correction.**
- For every intersection of two lines a and b, and every horizon picked on both, the misfit is d = t_a − t_b.
- The shifts s minimise Σ (d + s_a − s_b)². One line of each connected component of the intersection graph is held at s = 0.
  The normal equations are solved by Gaussian elimination with partial pivoting.
- Lines without intersections get component −1 and no shift.
- Applying the result adds s to the line's start time t0 and to its picks.
- Test: four pairs, rms misfit 6.671 → 0.000 ms, shifts 0 / −8 / +5 / 0 ms.

**Polygons.**
- Point-in-polygon uses the even–odd crossing rule. Area uses the shoelace formula; length is the sum of the segments.
- Clipping sets the excluded lattice nodes to NaN.
- The volumetrics mask is a per-cell 0 / 1 array: a cell counts when its centre is inside. Contacts per fault block override
  the global contacts.

**Layer-cake velocity model.**
- Interval velocity of zone k at well w: V = 2000 Δz / ΔTWT (m/s, TWT in ms), from the tops matched by name.
- Interface modes: mean of the wells, an inverse-distance map (12 nearest wells, `gridPoints` without a variogram), or a
  user value.
- Depth: z_k = z_{k−1} + V_k (t_k − t_{k−1}) / 2000. The residual at each well is corrected with an IDW map. Interfaces are
  forced non-decreasing.
- Test: V_int 2000 / 3030 m/s; rms tie 1.837 → 0.000 m.

**Checkshots (V0 + kZ).**
- Model: t = (2000/k) ln(1 + k h / V0), with t in ms and h the depth below datum; k → 0 gives t = 2000 h / V0.
- For each k, V0 is found by bisection on the mean residual (t is monotone in V0); k ∈ [0, 3] s⁻¹ by golden-section search
  on the rms.
- Test: V0 = 1800, k = 0.6 recovered exactly.

**Seismic calculator.**
- The expression is evaluated with the C42 calculator engine (`evalExpression`, geo_analysis) sample by sample, one inline
  at a time.
- Variables: AMP, T (ms), I, J. Any attribute whose name appears (ENVELOPE, PHASE, FREQ, RMS, …) is computed for that inline
  with the given window.
- Crop / decimate keeps the index ranges with a step in each direction.

**Well performance events WPWE0–7 (`Simulator::pweTrack`).** These follow OPM `data::WellPerformanceEvents` (opm-simulators
PR 7395). The tracker compares two snapshots, taken at the start of the first step and at the end of every accepted time step.
Each snapshot holds the set of open connection cells and whether the well is an injector.

| Event | Rule |
|---|---|
| WPWE1 | \|open_now \ open_before\|, only while the well is open with deck status OPEN |
| WPWE2 | \|open_before \ open_now\| minus the cells closed only by the reach of '+CON' (or ResSim PLUG) |
| WPWE3 | 1 for a reach past the offender, or when some connection closed and every open connection belongs to the shallowest completion (or none is open) |
| WPWE4 / 7 | simulator decisions (economic limits, failed WTEST) flagged at the closure site |
| WPWE5 / 6 | change of the injector flag between the snapshots |

Economic closures run after the step is accepted, so a retried step cannot report events twice or lose them. A connection
opened and closed again within one step leaves no WPWE1 / WPWE2 entry; this is a property of the snapshot comparison, as in OPM.

WTEST 'E' re-opens the connections closed by economic workovers when the well has none open. Their original COMPDAT end steps
are kept.

Test: '+CON' on the top completion gives WPWE2 = 1, WPWE3 = 1 and WPWE7 = 1 in the same step. CON on the top gives WPWE2 = 1 and
WPWE3 = 0. CON on the bottom gives WPWE3 = 1. With all connections limited and WTEST 'E' every 30 d, each closure gives
WPWE2 = 2 and each retest WPWE1 = 2.

WPWE0 is always 0 (no drilling queue). C44 verifies the tracker against the OPM reference output of the opm-tests `wpwe/` decks
and adds WPWE4 (WELSPECS item 9 = STOP); see §9.w.

### 9.w C44 additions: OPM reference checks, well-control fixes, hysteresis, restart

**WPWE against OPM Flow.** `tools/wpwe_compare.cpp` and the test `-c44` run the opm-tests `wpwe/` decks (`decks/opm_wpwe/`,
ODbL) and compare every WPWE*n* event series with the OPM reference summary at the reference report steps. Five decks are
identical (WORKOVER-SHUT-STOP, COMPLUMP-CONNECTION-COUNT, DECK-CONNECTION-CHANGES, DECK-STATUS-CONVERSION,
MIXED-CLOSURE-CAUSES). WTEST-REOPEN differs only in the re-open timing (accepted). The rules that were needed:

| Rule | Behaviour |
|---|---|
| WELSPECS item 9 | STOP: a simulator closure (economic WELL workover, failed test) stops the well -> WPWE4; SHUT -> WPWE7 |
| No open connection left | the well is shut (WPWE7) even with item 9 = STOP, as in OPM |
| WELOPEN on one layer | a multi-layer COMPDAT record is split into one record per layer; the other layers keep their state |
| Topmost connection | the first connection in deck order (`cComp`); splits insert the new record after the old one |
| Late well | a well defined later takes its baseline when it appears (no WPWE1 for the initial connections) |

**Injector cross-flow recirculation.** With cross-flow allowed (WELSPECS item 10 YES), fluid can flow from the reservoir into
an injector through connections where the reservoir pressure exceeds the wellbore pressure. Before C44 that fluid was counted
as surface production of the injector. In ECLIPSE and OPM it mixes in the wellbore and leaves again through the injecting
connections. ResSim now records, per injector and Newton update, the back-produced surface volumes $q^{x}_{\alpha}$ of the
non-injected phases and the total injected volume $Q_{inj}$. Each injecting connection $c$ with rate $q_c$ then carries the extra
sink $-q^{x}_{\alpha} q_c / Q_{inj}$ of each such phase (lagged by one assembly, constant in the Jacobian). Nothing is reported
as production when $Q_{inj} > 0$. Check: GUIDERAT_POTENTIAL_FALLBACK FGPR 4.0e5 = OPM at all report steps (was 1.33e6),
WBHP:INJ2 within 0.3 %.

**First guide rate.** OPM uses the well potential when a well has no guide-rate value yet (`GuideRate::get`, no entry in
`values`). ResSim did the same only for flowing wells. Now a well with no value and no flow also returns its potential. Check:
GUIDERAT_DEFAULT_WGRUPCON day 1, PROD1 1.88e5 vs OPM 1.64e5 (was 3.4e5); days 4-28 within 1-2 %.

**End-point scaling: no SWCR >= SWL clamp.** OPM `EclEpsScalingPoints::extractScaled` takes the per-cell end points as given.
ResSim used to raise SWCR to SWL and SGCR to SGL; the clamp is removed (only SWU > SWCR and the oil residuals below
1 - SWL are kept for a valid table).

**EHYSTR models 1 and 3.** For the wetting phases (water in oil-water, oil in gas-oil) models 0 and 2 use the drainage curve
and models 1 and 3 use the imbibition curve (OPM `EclHysteresisTwoPhaseLaw::twoPhaseSatKrw`). ResSim now evaluates krw and
krog on the imbibition end points for models 1 / 3; the non-wetting phases keep the Carlson scanning curves. Test: with
ISWCR = 0.4 the qfs example produces 41 % less water under model 1 than under model 0. The model2 9_* water cut at early
times now matches OPM (9_2B FWCT at day 1: 0.002869 vs 0.002867).

**WPIMULT per connection and step.** ECLIPSE and OPM multiply the connection factors of the connections that exist when the
keyword is read. The factor is cumulative, connections defined later do not inherit it, and a COMPDAT redefinition resets
it. ResSim used one well multiplier for the whole run, so a WPIMULT in 1998 also reduced the injectivity of connections opened
in 2001, and repeated WPIMULT records compounded from time zero. Now `Completion::mult` holds the cumulative factor of each
open period. WPIMULT at step $n$ ends the open period of each matching connection at $n$ and starts a new period with
mult × factor. Connections opened in the same step are multiplied in place. I / J / K select connections (one layer of a
multi-layer record is split off); completion numbers are kept for OPM. The export writes COMPDAT and then one WPIMULT record
per layer at the step where each period starts. `WellInput::piMult` stays as a whole-run multiplier for the UI and studies.

**WELPI for injectors.** The injector PI is $\sum_k WI_k\,\lambda_{t,k}\,b_{\alpha}$, where $\alpha$ is the phase injected under the
controls of the WELPI step. WELPI is now applied *after* that step's controls; before C44 it ran before them, so a WELPI in the
step where a water injector turns to gas used $b_w$ instead of $b_g$. Check: 0B_WELPI_MODEL2 INJ2 (WELPI 1e4 in the switch step)
WBHP 315.7 vs OPM 315.8 bar (was 266), and 410.8 vs 411.0 at 2.28e6 sm3/d. Still open in this deck: PROD3 BHP about 2.5 % high
after its WELPI, and a period of lower PROD2 / field oil rate around day 350.

**Export of PI multipliers.** COMPDAT is followed by one WPIMULT record per open period that starts at that step. The record
uses K = 0 (the whole column) unless the well has another connection range in that column; then one record per layer is
written. A whole-run well multiplier (`piMult`, set in the UI or by studies) is folded into these factors, so OPM applies it
too. A `-- RESSIM: PIMULT well f` comment lets ResSim restore it on import.

**FBHPDEF.** The default producer / injector BHP limits (items 1, 2) apply to the well controls that follow when their BHP item
is defaulted. This includes WCONHIST and WCONINJH, whose limits are otherwise 1 atm and 100000 psia. The keyword is kept in the
deck and the controls are written with the BHP defaulted, so the round trip is unchanged. Together with WPIMULT this fixes
0C_BASE_FBHPDEF: INJ1 reached 5764 bar before (OPM 500); now WBHP:INJ1 is 439.9 vs 443.5 and FOPR 8446 vs 8436 at day 823.

**Restart state v15.** The black-oil restart blob now also stores, per well: the connections closed by economic workovers
(used by WTEST 'E' to re-open them), the WPWE baseline (open connection cells, injector flag, last events) and the injector
recirculation totals. Test: a restart taken while all connections are closed by CON gives the same re-opening, WPWE1 / WPWE2
series and FOPT as the continuous run. Older blobs (v1-v14) still load.

**Regression against opm-tests (63 cases: krnum, network, model2).** Changes from reg2 to reg3 (recirculation + first guide
rate): no case got worse; GUIDERAT_DEFAULT_WGRUPCON material metric 15.5 % -> 3.2 %, 9_3C 20.1 -> 4.95 %, 9_3D 32.5 -> 8.5 %,
9_3F 25 -> 14.7 %, 0A1 / 0A2 GRCTRL STW 10.4 / 14.4 -> 9.9 / 14.1 %. GUIDERAT_POTENTIAL_FALLBACK still fails per well because OPM
alternates the group allocation between wells from step to step; the field rates are identical. Still open: KRNUM-03Z
(earlier gas breakthrough in OPM, vertical flow), NETWORK-01-WTEST WBHP:PROD2, 3_A_MPI_MULTFLT_SCHED, 6_ENDSCALE (4.7 %
material), CO2STORE, and the late-time FOPR drift of the 9_* GCONPROD decks.

After WPIMULT / FBHPDEF / WELPI (reg5): **18 pass** (was 17). 0C_BASE_FBHPDEF passes (29.3 % -> 4.2 %), 0B_WELPI_MODEL2
13.4 % -> 11.7 % (material 8.8 %), and no other case changed by more than 0.2 %. Report: `docs/regression_C44_opm.html`.

**Verification C44.** `-c44`: 20 checks (WPWE vs OPM ×5, STOP, WELOPEN layer, fault cut ×3, restart v15, EHYSTR 1, injector
recirculation vs OPM, WPIMULT / FBHPDEF import + round trip). Full suite: **1073 checks pass, 0 fail**.

### 9.v C45 additions: memory of large models, multisegment cross-flow (experimental), panel translation

**Block assembly (no scalar Jacobian).** With the block solver on (default), the cell rows of the Jacobian are no longer
assembled into the scalar CSR `J_` and copied to BSR. The assembly writes every nv×nv cell block straight into the
`WellSchur` block matrix (`schur_.bindCells`). `J_` keeps only the well rows and the cell→well columns. These are
the parts the Schur complement needs. The connection blocks are found through a per-cell sorted neighbour list
(`blkStart_`, `selfPos_`, `adjPos_`) built once in `buildPattern`. Paths that still need the full scalar matrix switch
back to it:
* `RESSIM_DUMP` / `RESSIM_DUMP_FAIL` (matrix dumps);
* the adjoint (`ensureScalarJacobian()` rebuilds the pattern with `forceScalarJ_`, since the adjoint transposes `J_`);
* `RESSIM_SCALAR_JAC=1`, the switch for checking the old path.

The output is bit-identical. The arithmetic per entry is the same and only the storage location changes.

**Lean linear algebra** (`-- RESSIM: TUNING` item 10 `leanMemory` = −1 auto / 0 off / 1 on; `SOLVER` fourth token; UI
*Solver → Low-memory linear solver*; `RESSIM_LEAN=0|1` overrides).
* *BILU0:* the second stage factorises each row in double scratch and stores it straight into float32 (`bilu0RowF`).
  No double copy of the factors exists at any time. The Krylov vectors, the residual and the SpMV stay double.
* *Rescue:* if a lean solve stalls (relative residual > 1e-3), a double factorisation of the same matrix takes over for
  that solve. This rescue is shared with mixed precision.
* *Auto mode:* switches on from 500 000 active cells. Below that, results are unchanged.
* *Always on, any size:* AMG takes over the fine pressure matrix (`setup(CSR&&)`) instead of copying it, and the CPR
  accumulation scratch `accD_` is released after the pressure matrix has been formed.
* *Solvent:* solvent face arrays (`cfS_`, `cfDirS_`) exist only when SOLVENT is active.
* *Grid connections:* `Connection` is packed from 64 to 56 bytes.

Measured peak RSS on 2 threads (`tools/gen_bigcase.py`, 6 months; FOPT/FGPT/FPR identical to the C44 run):

| cells | C44 | C45 default | C45 lean |
|---|---|---|---|
| 25 k | 111 MB | 94 MB | — |
| 100 k | 416 MB | 322 MB | 294 MB |
| 200 k | 801 MB | — | 584 MB |

The 200 k lean run is 2.88 kB per cell, which extrapolates to about 2.75 GiB at 1 M cells. That is below the 3 GB target
and below the 4.0 GB estimated in §6.39. The 1 M-cell run itself still cannot be done on the 2 GB build machine.
In the 200 k lean run (9 steps, 30 Newton, 333 linear iterations), the float factors converged without the rescue.

**Producer cross-flow in multisegment wells (experimental, off by default).** Switch it on with
`-- RESSIM: MSWOPTS flash storage 1` (third token) or `RESSIM_MSW_XFLOW=1`.
* *Injected fluid:* a cross-flowing connection (p_cell < p_conn) injects the fluid of its own segment, following OPM
  MultisegmentWell. That is the surface-fraction mixture of the gross inflow of the segment's subtree (`mswMixture`,
  `xmixConn`; restart state version 16), instead of the well-level mixture used by standard wells (§9.x).
* *Storage sweep:* the injected reservoir volume is carried down from the top through every segment on the connection's
  path (`SolveOpts::xfSweep`), so the segment-storage content (§6.36) follows the fluid that actually passes.
* *First fill:* the tubing's first content is the top-segment mixture, not the steady value of the potential rates.
* *Not counted as segment flow:* cross-flowing connections still count as zero in the segment flow and pressure-drop
  solution (`mswHeads` uses |q_res| of producing connections).

*Why it stays off.* The connection heads come from the explicit start-of-step solution of the segment network. With
cross-flow on, MSW-2D-VERT-02 gets worse: FGOR 3.6 % becomes 23.8 % against OPM. The mechanism is a flip-flop between steps:
1. A light head makes the lower connection produce oil.
2. Under ORAT control, the BHP rises and the upper (gas) connection stops.
3. The tubing fills with liquid, and the next step's heavy head makes the lower connection inject about 10⁶ sm³/d of gas.

Two attempts did not remove the flip-flop: under-relaxing the top-segment density (58 %) and the first-fill and sweep
fixes alone (23.8 %). The cure is an implicit coupling of the segment pressures and fractions, as in OPM. No regression case
improved materially with cross-flow (0A1–0A4 unchanged; 9_3B/C/D within ±0.4 %; BASE_MSW_WSEGAICD 14.4 → 14.9 %).
The default therefore stays the C44 behaviour: back-flowing MSW connections are shut for the step.
*C49 update:* the flip-flop is cured by the converged connection heads of §9.r (FGOR 23.8 % → 1.87 % with cross-flow).
Cross-flow stays off by default because no OPM case improves with it (§9.r).

**Translation of panels (i18n).** The UI string table now covers the Model / Results panels, property and table widgets
and the help tooltips: 2 011 keys. `lang/es.txt` has 1 660 of the 1 697 new keys translated; the other 37 stay English
on purpose (units, keyword names, formats). Dynamic strings go through `TRS` (cached by content). Labels translated
inside ImGui use `TRL` (`translated###English`), so widget IDs, saved layouts and selections do not depend on the language.
Plot titles are translated with `TRL`, so the plot IDs stay English. Popup and window names, plot legends (series IDs) and the
built-in keyword documentation stay English. `tools/i18n_extract.py --update lang/xx.txt`
appends new keys to a language file.

### 9.u C46 additions: well PI, completion numbers, SCHEDULE fault multipliers

**WELPI.** The PI that is reported and rescaled is the OPM `computeConnLevelProdInd` value. For each open connection

  PI_p = CF_k · λ_p / B_p.

A producer with both oil and gas adds the cross terms

  PI_o += Rv · PI_g,  PI_g += Rs · PI_o,

both computed from the values before the addition. An injector uses its injected phase.

At a WELPI event the factor is f = PI_target / PI_current. It multiplies every connection of the well defined up to that
step (`cWelpi`), including shut ones, which matches OPM `prepareWellPIScaling`. A COMPDAT record for an existing cell
after the event starts a new record without the factor (`welpiSince`). ACTIONX rebuilds of the connection list carry the
factor over to the rebuilt records.

**Completion numbers.** COMPLUMP numbers apply where given. Elsewhere the number is the order of first definition, one
per layer (OPM `WellConnections::addConnection`). WELOPEN and WPIMULT items 6–7 select by these numbers. A record
selected only in part is split. The split carries any later closing step (`untilStep`) into the new period.

**SCHEDULE MULTFLT.** Each event (step s, fault pattern, factor m) multiplies the current transmissibility of every
connection across the fault faces from step s onward, NNCs of displaced columns included:

  T_c ← m · T_c.

Repeated factors are cumulative (OPM `applyModifierDeck`). Restarts re-apply every event with s ≤ the restart step.

**Open finding (6_ENDSCALE).** Initial SWAT matches the OPM restart file cell by cell, to 2·10⁻¹¹. Late water fronts in
cells with SWL below the table connate saturation move more slowly than in the reference. Evaluating krow without end-point
scaling reproduces that deck (worst 5.5 → 2.7 %), but it makes 7_HYSTERESIS, 8_MULTIPLY and 9_EDITNNC worse. Two-point
krow scaling is mathematically the same as OPM's current `EclEpsScalingPoints` and `EclEpsTwoPhaseLaw`, so the formula
was not changed. The cause is still open.

### 9.t C47 additions: low-pressure gas tables, OPM reference drift, fault history export

**Low-pressure padding of gas PVT tables (OPM `LowPressureTablePadding`, opm-common PR 3779).** Tables are evaluated
with linear extrapolation in b = 1/B. A dry-gas table that starts well above 1 bar therefore predicts b < 0 near
atmospheric pressure. Such pressures occur in well potentials at a BHP of 1 atm. ResSim now pads PVDG and PVTG tables
the way OPM Flow does. Let rows 0 and 1 have pressures p₀ and p₁, bᵢ = 1/Bᵢ and (b/μ)ᵢ. Then:

  B_lim = max(10·B₀, 1),  p_lim = p where the line through (p₀, b₀), (p₁, b₁) gives b = 1/B_lim,
  Rv_lim = Rv(p_lim),  μ_lim = 1 / (B_lim · (b/μ)(p_lim))   (both linear in p through rows 0 and 1).

Padding is needed if p₀ > 1 bar and b(1 bar) < 1/B_lim. Otherwise it is needed if p₀ > 0 and b(0) < 0; in that case the
limit row is (0, 1.1·B₀, Rv₀, μ₀). If p_lim < p₀, these rows are prepended:
* (1 bar, Rv_lim, 1.1·B_lim, μ_lim), only when p_lim > 1 bar;
* (p_lim, Rv_lim, B_lim, μ_lim).

For PVTG, each padding row is a saturated record without undersaturated rows. The existing record-completion logic then
gives it the undersaturated slope of the next record. The pressure 1 bar is 14.5038 psia for FIELD decks.
`FluidModel::gasPaddingRows()` reports the number of rows added. Example: the opm-tests model2 table starts at 50 bar
(B = 0.024958, 0.017643 at 70 bar). There b(1 bar) ≈ −0.62, so rows are added at 1 bar (B = 1.1) and at 2.95 bar (B = 1).

The padding changes nothing on 9_3A, 9_4B, 9_1A, 6_ENDSCALE and 0A1_STW: the results are identical with and without it,
because no cell or well reaches those pressures. It is still needed to keep potentials and low-pressure producers
physical, and to stay consistent with OPM ≥ 2025.

**OPM reference drift (C47 study).** The opm-tests reference summaries are regenerated automatically by OPM pull
requests. We compared each case's own reference history (`tools/dev/refdrift*.py/.cpp`, 18 cases) against the time-aligned
difference (`tools/dev/smrydiff.py`; OPM UNSMRY rows are ministeps, so rows are aligned by TIME). Findings:
* **Stable references** (OPM's own versions differ by ≤ 2.5 %): 6_ENDSCALE, 9_1A, 0A1, 0B_WELPI, ACTIONX_M1, ACTIONX_WEFAC,
  7_HYSTERESIS, WSEGSICD / AICD. 9_4B differs by ≤ 4.8 %. A ResSim deviation above the tolerance on these cases is a ResSim
  issue.
* **Moving targets** (OPM's own versions differ by more than the 5 % tolerance): NETWORK-01-WTEST 40.8 %, T1L 35.1 %,
  UDQ_WCONPROD 30.2 %, UDQ_ACTIONX 53.3 %, 9_3A 12–16 %.
* **Bisection** to the first reference commit that changed:
  * T1L changed with opm-simulators PR 6311 (group targets kept during local well switching). ResSim passes the earlier
    reference (3.2 %).
  * 9_3A changed with opm-common PR 3779 (the padding above). ResSim reproduces the 2024-09 reference to 5.4 %.
* **9_3A is reference-unstable.** The 2024-09 and 2025-02 OPM references differ by 45 % in FWCT, 43 % in WWPR:PROD3 and 35 % in
  WGPR:PROD1, all at t = 151 d, and by 31 % in WGIR:INJ1 at t = 136 d. They also take different numbers of ministeps (31 vs
  26). A change that only adds table rows at ≤ 3 bar therefore flips a transient. Deviations on this case measure
  transient timing, not physics, and are not pursued.

**Wells completed only by ACTIONX (opm-tests ACTIONX_GCONINJE / ACTIONX_WCONINJE).** Suppose the WELSPECS reference depth
is defaulted and a well has no connection at setup. Its reference depth was then 0, and `rebuildActionConnections` kept that
value after the action added connections. The wellbore head bhp → connection became ρ_w·g·z ≈ 261 bar at 2667 m. The water
injector could not lower its rate below 21 000 sm³/d even at the BHP floor, so the rate equation stalled and the run stopped
with "time step below minimum" at t = 578 d. The rebuilt well now takes the defaulted reference depth from its connections
(OPM `Well::updateRefDepth`). In addition, the ACTIONX month token `JLY` (the ECLIPSE spelling of July, already accepted in
DATES) is now recognised in conditions. Results: both cases go from ERROR to FAIL. ACTIONX_WCONINJE (14.3 %) differs only by
one reference sample, 6895 bar at the report step where the action opens WI01 (OPM reports the BHP limit before the well
flows). The rest of its WBHP / WWIR history matches, including the 3500 → 4000 sm³/d switch of ACT-04. ACTIONX_GCONINJE
(10.3 %) has the WGPR:OP02 offset common to this model family (also in ACTIONX_WEFAC / GCONPROD). Debug aid:
`RESSIM_DEBUG=1 RESSIM_DEBUG_WELLS=1` prints the open flag, mode, target, rate, stall counter and connection pressures of
every well at each Newton iteration.

**Per-face fault history (`--fault-csv FILE`, environment `RESSIM_FAULT_CSV`).** After a FAULTMECH run,
`faultmech::writeFaultCsv` writes one row per report step and fault face with the columns
`t, face, fault, k (1-based), x, y, depth, area, aSig, cff, dcfs, rate, cum`. Time is in days. Stresses, coordinates and
areas are in the case units: `Output::paPerUnit` and `mPerUnit` convert them to Pa and m. `rate` is the seismicity rate
relative to background, and `cum` is the cumulative expected number of events. Cycle 9 of the research track uses this file. CO2STORE decks no longer
print the "not modelled" warning: they run on the compositional CO2STORE core (§5).

### 9.r C49 additions: converged multisegment heads, fault-permeability ramp, translated window titles

**Converged connection heads of multisegment wells (Picard coupling).** Up to C48 the connection heads h (pressure of each
connection relative to the BHP) came from one segment-network solve with the start-of-step state; they were then held fixed
during the Newton iterations. Now, after the step converges, the segment network is solved again from the *end-of-step*
state, G(h). If any head moved by more than 0.02 atm, the step is re-solved (warm start, as in the VREP loop) with new
heads. This repeats for at most 8 passes.

* *Why not plain substitution.* On MSW-2D-VERT-02 the map G is a near step: the 512 m tubing segment above the lower
  connection is gas-filled while that connection cross-flows and liquid-filled as soon as it produces (26 bar). Plain
  substitution h ← G(h) and Aitken / Irons–Tuck relaxation both ended in a 2-cycle.
* *What is used.* G is decreasing in each head (a heavier head gives less inflow and a lighter column). Each head is
  therefore bracketed separately: regula falsi with the Illinois modification on r_j = G_j(h) − h_j, the same scheme as the
  network balance (§6.23, C30). The first trial is the substitution value.
* *Convergence.* A step counts as converged when every residual is below 0.02 atm, or when every remaining head is bracketed
  within 0.02 atm (a jump of G). Otherwise, after 8 passes the step keeps the heads it was solved with. On MSW-2D-VERT-02
  with cross-flow, 40 of 58 steps end at the pass limit with a bracket that is still open; the result is nevertheless close
  to the run without cross-flow (FGPT −0.19 %).
* *Limits, stated plainly.* This is an outer (Picard) coupling, not OPM's fully implicit segment equations: the Jacobian has
  no ∂h/∂x terms, and each pass costs one more Newton solve. MSW decks run about 1.1–2.6× longer (the 24 OPM MSW cases:
  390 s → 416 s in total). Adjoint runs keep the explicit heads, because the adjoint differentiates the start-of-step head
  update.
* *Option.* `-- RESSIM: MSWOPTS flash storage xflow heads`: the fourth token is 1 (converged), 0 (explicit, as before C49)
  or −1 (default = converged). `RESSIM_MSW_HEADS=0/1` overrides the deck. The UI is in *Model › Wells › Multisegment well*.
  The run log ends with one statistics line (steps checked, extra solves, steps left at the pass limit).
  `RESSIM_DEBUG_MSWHEADS=1` prints the residuals of every pass and the segment state.

*Results.* A/B test on the 24 OPM MSW reference cases, worst relative RMS: no case is worse.
* TEST_WSEGSICD: 93.4 % → 34.0 %.
* BASE_MSW_WSEGAICD: 24.4 % → 24.1 %.
* 0A4_MSW: 4.26 % → 4.20 %.
* All others are within ±0.05 %.

MSW-2D-VERT-02 FGOR against OPM:

| Cross-flow | Explicit heads | Converged heads |
|---|---|---|
| off | 1.29 % | 1.29 % |
| on | 23.8 % (flip-flop, FOPT +3.9 %) | 1.87 % (FOPT exact) |

With converged heads, cross-flow changes no OPM case by more than 0.1 %, so it stays off by default.

**FAULTMECH ramp width (item 10).** The two-way transmissibility factor of a reactivated face is permFactor^x with
x = clamp(CFF / (w σ'n), 0, 1). w was fixed at 0.05; it is now item 10 (default 0.05, written only when it differs).
Within the ramp, ln(multiplier) = ln(permFactor) · CFF / (w σ'n). A wider ramp is therefore *the same shape with a smaller
slope*, not a smoother shape. Research Cycle 11 (H10f) uses this. The fault CSV (`--fault-csv`) has two more columns:
`sne` (σ'n) and `mult` (the transmissibility multiplier in force).

**Window titles and dock tabs translated.** A display-only hook in Dear ImGui (`GImGuiTitleTranslate`, used in
RenderWindowTitleBarContents, the docked-window tab label, TabItemCalcSize(window) and TabBarGetTabName) translates the
visible title. The window name, hence the ImGui ID, `ressim_layout.ini`, DockBuilder and popup IDs, stays English.
18 new Spanish titles.

**Stopped wells report a BHP.** A well with STOP status (WCONPROD / WCONINJE item 3 or WELOPEN) has no surface flow. It
now reports the zero-inflow BHP as WBHP: p_bhp = Σ_k WI_k (p_k − Δh_k) / Σ_k WI_k over the open connections, the same
formula C33 uses for open wells with a zero rate target. SHUT wells keep WBHP = 0. OPM and ECLIPSE do the same.
opm-tests SPE1CASE1_SHUT_CONTROLS (2026) checks this: WBHP of S-ORAT / S-GRUP / FS-ORAT / FS-GRUP went from 99 / 99 / 70 /
70 % to 0.005–0.06 %. A STOP well is still modelled without cross-flow between its connections; OPM lets a stopped well
cross-flow. The small FS-* residual comes from that.

**Per-well summary of compositional / CO2STORE / thermal runs via the API.** `api::runCase`, which the regression runner,
studies and the C API use, converted the compositional and thermal report steps into field totals only. The per-well rates,
BHP, FWCT and FGOR were lost: an OPM CO2STORE comparison showed WBHP:INJ = 0. They are now kept, together with the well
names. CO2STORE and CO2STORE_DIFFUSIVE now pass (worst FPR 0.59 %).

**Clean skip for SOLVENT + CO2SOL / H2SOL.** OPM's model takes the solvent properties from CO2 / H2 tables instead of PVDS.
These decks are now marked "not supported by the ResSim engine" (regression SKIP) instead of failing in the run. Boundary
conditions (BCCON / BCPROP) are still not modelled; they are reported as warnings at import (CO2STORE_DRSDTCON,
CO2STORE_GW_DIRICHLET fail for that reason).

**Full OPM regression (209 cases).** The C48 rerun of all cases in `tools/fetch_opm_tests.py` was completed: 102 PASS /
77 FAIL / 30 SKIP at the 5 % relative-RMS criterion (`docs/regression_209_full.txt`). The criterion was not loosened. The
largest groups of failures:
* the 9_* network / group model family (WWPR of a low-rate producer, 21–30 %);
* KRNUM-02Z / 03*;
* the tracer cases with water-tracer coupling;
* MULTREGT-01 / MULTFLT-01;
* boundary-condition decks.

**Known limitation kept: ACTIONX_GCONPROD.** After the ACT-03 control jump at 1097 d, ResSim's 31-day Newton step fails and
is cut (31 → 9.3 d; then 10.85, 10.85, 14, 14, 15.5, 15.5 d), while OPM converges monthly steps. The deck's compounding UDA
(FU_FGPR = 0.9 · FGPR per step) makes the target depend on the number of steps, which explains the remaining 17 %
(F34). This is documented, not limited: forcing report-sized steps would hide a Newton failure.

### 9.s C48 additions: well reference depth, group UDA, gas-aware sub-grid leak-off, comp / thermal snapshots

**Defaulted BHP reference depth.** If WELSPECS item 5 is defaulted, the reference depth is now the shallowest cell-centre
depth of the connections the well had at its *first* definition step, **SHUT connections included**. This is OPM
`Well::updateRefDepth`: `ref_depth = connections.get(0).depth()`, set once. ResSim used to take the shallowest *open*
connection.

In the ACTIONX model family, OP02 defines K3–4 and K7–10 SHUT and later opens COMPLUMP 2 (K7–10) with WELOPEN. The reported
WBHP was therefore about 3 bar high, and the 200-bar limit was reached 180 d late. Results after the fix:
* ACTIONX_WEFAC: 9.46 % → 1.8 % (pass).
* ACTIONX_GCONINJE: 10.3 % → 1.04 % (pass).

Connections added later by COMPDAT do not move the reference depth. Wells that ACTIONX completes at run time keep the C47
rule (their first run-time connections).

**Group UDA.** GCONPROD items 3–6 and 14 (ORAT WRAT GRAT LRAT RESV) and GCONINJE items 4–7 (RATE RESV REIN VREP) may hold an
`FU_` / `GU_` UDQ name. The name is stored in `GroupProdControl::uda` / `GroupInjControl::uda`. `Simulator::applyUda()` sets
the value every time step before the group targets are allocated; an undefined UDQ keeps the previous value. The keyword is
written back verbatim.

ACTIONX_GCONPROD (ACT-02: GRAT = FU_FGPR = 0.9 FGPR) is now exact up to 1066 d (FGPR 0.00 %); it was 69.8 % → 17.1 % overall.
The remaining difference comes from the deck design. The UDA compounds once per *time step*, and after the ACT-03 control jump
ResSim takes ministeps while OPM keeps monthly steps. That makes the result step-size dependent, and it is accepted.

**SUBGRID LEAKC (H2e, research F32).** Option of HFGROW SUBGRID (`... SUBGRID LEAKC`); the default is off. The leaked fluid in
the exposed sub-cell is stored with the sub-cell's own total compressibility:

    p_sub = p_res + (p_c − p_res) · (c_t,cell / c_t,sub) / φ,     S_g,sub = S_g0 + ΔS_g / φ

Here c_t comes from the cell phase compressibilities (`cellCt()`), and S_g0 is the gas saturation when the tip enters the
cell (`hfSg0_`). For water, c_t,sub = c_t,cell, so the result is identical to SUBGRID LEAK. For dry gas at 600 sm³/d it
removes the t90 over-correction: spread 0.81 d instead of 2.19 d over dx 10 / 5 / 2.5 m.

**Comp / thermal 3-D snapshots (`cell_snapshot.h`).** `comp::cellSnapshot` and `thermal::cellSnapshot` build the
per-report-step cell maps (PRESSURE, saturations, XCO2, SGTRAP, ZMFk/XMFk/YMFk, TEMP/DTEMP, SSTEAM, THM maps …). The GUI runner
and `ResSimCLI -o out.rsres` share them, so a headless compositional, CO2STORE or thermal run now writes the same 3-D results
as the GUI. `coupleddeck::runComp` / `runThermal` take an `onState` callback.

**CLI `--csv` well columns.** For black-oil runs, `--csv` also writes WOPR / WWPR / WGPR / WWIR / WGIR / WBHP:<well>.

**Closed as investigated, unresolved.**
* KRNUM-03X/Y/Z: the directional relperm and the connection-mobility logic match OPM's FlowProblem / WellInterface code.
* 0A1 GRCTRL STW: ResSim STW is within 0.9–2.9 % of OPM's *MSW* reference, while OPM's own STW and MSW versions differ by 17 %
  on WWPR:PROD3. Our standard-well head algorithm matches StandardWellConnections line by line, so the 0A1–0A4 STW failures are
  recorded as an OPM standard-well treatment difference.

**Regression (C48, model2 + actionx + krnum, 73 cases).** 31 pass, 40 fail, 2 skip
(`docs/regression_C48_model2_actionx_krnum.html`). Against C46 there are no regressions. ACTIONX_WEFAC and ACTIONX_GCONINJE now
pass, and ACTIONX_GCONPROD improved from 69.8 % to 17.1 %.
The full 209-case rerun was stopped after 115 cases for the release: 48 pass, 52 fail, 15 skip
(`docs/regression_C48_full_partial.txt`). The remaining 94 cases move to C49.

## 10. EOR and PVT laboratory tools (`eor.cpp`, `eos.cpp`)

- **Swelling test / P-X diagram:** saturation pressure versus mol % of injection gas. The first-contact MMP is the
  cricondenbar of the pseudo-binary.
- **Multi-contact MMP, mixing-cell method (Ahmadi & Johns 2011):** key tie-line lengths are extrapolated to zero.
- **Slim-tube MMP (C7):**
  - Each displacement runs a 1-D compositional slim tube on the FIM engine.
  - Recovery at 1.2 PV injected is computed at several pressures.
  - The MMP is taken at the breakover of the recovery-versus-pressure curve (two-line intersection).
- **MMP correlations:** Yellig–Metcalfe and Cronquist.
- **WAG schedule generator.**

---

## 11. Streamlines (`streamline.cpp`)

- Pollock semi-analytical tracing in each cell's reference unit cube, using the total face fluxes.
- Works on corner-point grids, faults and NNCs.
- Outputs time-of-flight, drainage and sweep regions.

---

## 12. Studies: sensitivity, uncertainty, history matching (`study.cpp`)

- **Sensitivity:** tornado (one at a time) and Latin-hypercube Monte Carlo.
- **Proxy models:** design of experiments with a quadratic response surface.
- **History matching:**
  - Levenberg–Marquardt on normalised parameters, with a finite-difference Jacobian evaluated in a parallel batch.
  - ES-MDA (Emerick & Reynolds 2013): ensemble smoother with multiple data assimilation and inflation factors
    Σ α_i⁻¹ = 1.
  - **C9: analytic gradients** (`adjoint.cpp`, `study.cpp`). See §12.1.

### 12.0 Inversion statistics (C27, iTOUGH2 style, `hmStatistics`)

**Setup.** At the matched parameters p̂ (normalised, n of them), with m weighted residuals r = W^{1/2}(d_obs − d(p̂)) and
the Jacobian J from n + 1 runs (forward differences, run in parallel):
- error variance s₀² = rᵀr / (m − n);
- covariance C_pp = s₀² (JᵀJ)⁻¹ through the symmetric eigen-decomposition (Jacobi; condition number reported);
- standard deviation σ_k = √C_kk, mapped to physical units (log-scaled parameters as log10);
- 95 % interval p̂_k ± 1.96 σ_k (normal approximation, m ≫ n);
- correlation R_kl = C_kl / (σ_k σ_l);
- composite scaled sensitivity CSS_k = √(Σ_i (J_ik p_k)² / m).

**Residuals.** Mean, standard deviation, skewness and kurtosis are reported, plus:
- leverage h_i = diag(J(JᵀJ)⁻¹Jᵀ), with Σ h_i = n;
- studentised residuals W_i = r_i / (s₀ √(1 − h_i)); |W| > 1.96 = outlier;
- the share of each observation group (well × vector) in the objective.

**Limits.** The analysis is linear (FOSM) and valid near the optimum. Strong correlations (|R| → 1) mean the data cannot
separate the parameters; the C27 test finds krw–μw +0.998 (mobility-ratio non-uniqueness).

**Verified (-c27):** symmetric covariance with unit-diagonal correlation; 95 % intervals contain the truth (1 % noise);
Σ shares = 100 %; Σ h = 3.0000.

**Coupled responses in studies (C27).** `runSingle` runs the coupled black-oil driver (`coupled_bo.cpp`) when the case has
fault mechanics / FE geomechanics / DDM. The responses expected events, max Mw, max ΔCFS, max uplift and max failure are
read from the run KPIs, so tornado / MC / HM / statistics work on them. FAULT_SLIP gives 18.61 expected events (C26 CLI:
identical ΔCFS 5.296). Tornado ×0.5 / ×1.5 injection gives 16.69 / 18.61: the ×1.5 case equals the base because the
injector is BHP-limited.

**Van Genuchten–Mualem tables (C27, `tools.cpp`).** Se = (S − S_r)/(1 − S_r − S_nr);
k_rw = k_rw,max √Se [1 − (1 − Se^{1/m})^m]²; non-wetting phase either Parker (1 − Se)^{1/2}(1 − Se^{1/m})^{2m} or the TOUGH
Corey form; Pc = P₀ (Se^{−1/m} − 1)^{1−m}, capped. Written as SWOF or SGOF rows. *Verified:* max |k_rw − formula| 3e-17,
k_rw(0.5; m 0.5) = 0.0127.

### 12.1 Discrete adjoint and forward sensitivities (C9)

The misfit F(θ) = Σ r_i²/n_obs, with r_i = (y_i − d_i)/σ_i, is differentiated exactly through the *discrete* black-oil
model. Each step n solves R_n(x_n, x_{n−1}, d_n; θ) = 0, where x_n holds the cell primaries and the well BHPs. The
connection heads d_n are explicit: computeHeads runs at the start of the step from the previous state, so
d_n = H(x_{n−1}, bhp_{n−1}, d_{n−1}; θ). They are treated as an augmented lagged state with their own multipliers μ_n.

- **Recording:** the forward run stores the converged x_n, dt_n and the inputs of the head update. A *replay* model then
  rebuilds J_n = ∂R_n/∂x_n at the recorded states. The replay deviation ‖R_n(recorded)‖ is reported and is 0 in the
  tests.
- **Adjoint (reverse) recursion**, for n = N … 1:
  - J_nᵀ λ_n = −∂F/∂x_n + D_accᵀ λ_{n+1}/PV + H_{n+1,x}ᵀ μ_{n+1}
  - μ_n = H_{n+1,d}ᵀ μ_{n+1} − R_{n,d}ᵀ λ_n − ∂F/∂d_n
- **Gradient:**
  dF/dθ = Σ_obs ∂F/∂y · ∂y/∂θ|_explicit + Σ_n λ_nᵀ ∂R_n/∂θ − Σ_n μ_nᵀ ∂H_n/∂θ.
  - ∂R/∂θ and ∂H/∂θ are central differences of the residual / head update of perturbed replay models at the *recorded*
    state (no perturbed simulations). This keeps every parameter type generic: KX/KY/KZ multipliers, PORO, rel-perm
    end points, viscosities, well PI, aquifer, …
  - Cost: one forward run + one backward pass + 2 × n_param residual evaluations per step. The parameter contractions
    run in parallel threads.
- **Forward (tangent) mode** gives all sensitivities dy/dθ for the Jacobian of Levenberg–Marquardt (“LM analytic”):
  - sd_n = H_θ + H_x s_{n−1} + H_d sd_{n−1}
  - J_n s_n = −R_θ + D_acc s_{n−1}/PV − R_d sd_n
  - The parameter columns are solved in parallel.
- **Observation operator:** WBHP, rates (WOPR/WWPR/WGPR/WWIR, F*), FPR, WWCT/WGOR. Rates depend on the cell
  state, the BHP and the head of each connection (p_conn = bhp + d_k). The head term ∂y/∂d is what makes rate
  observations of rate-controlled wells exactly consistent with the well equation.
- **Linear solves:** ILU(0)-GMRES at 1e-12 relative residual, with BiCGSTAB refinement as a fallback.
- **Frozen lagged data:** group / THP targets, the RESV conversion, hysteresis history, DRSDT limits and aquifer state
  keep their recorded values (frozen-coefficient adjoint). Control and phase switching events make F piecewise
  smooth. The gradient is exact on the current branch. Finite differences of *full runs* across a switching event
  differ by O(1e-3).
- **Verification (`-adj`):**

  | Check | Result |
  |---|---|
  | adjoint vs forward mode (duality) | 9e-9 |
  | adjoint vs central FD of full runs, 40-step horizon without switching events | 6e-7 |
  | same, full 240-step horizon (crosses control switches) | 1.4e-3 |

  Per-step BHP sensitivities match FD to 1e-6 – 1e-9 away from the events.
- **Optimiser:** projected L-BFGS (memory 6) with Armijo backtracking in the normalised box [0,1]ⁿ. Each iteration costs
  one recorded forward run + one backward pass, independent of the number of parameters. On the test case it reaches
  the truth (misfit 14.7 → 0.022) in 10 simulations, against 17 for LM with the analytic Jacobian.
- **UI:** Studies → Levenberg–Marquardt → *Optimiser / gradient* (FD LM, analytic LM, adjoint L-BFGS). *Check gradient*
  shows a table of adjoint, forward and FD values with relative errors.
- **C12: analytic field gradients (`AdjFieldReq`).** For per-cell PERM / PORO multipliers, per-face transmissibility and
  per-connection WI, the term λ_nᵀ ∂R_n/∂ln θ is accumulated analytically during the reverse pass. There are no replay
  differences, so the cost does not depend on the number of parameters:
  - Flux is linear in T_f, and T_f is the harmonic average of the two half-transmissibilities, so
    ∂T_f/∂ln k_a = T_f²/t_a.
  - The connection flux is linear in WI ∝ k.
  - The accumulation is linear in the pore volume.
  - Results (`-adjf`):

    | Check | Result |
    |---|---|
    | Adjoint vs central FD of full runs | 4.3e-6 |
    | Gradient time | 0.25 s (forward run 0.70 s) |
    | Per-cell permeability HM with projected L-BFGS (UI: Studies → field HM) | misfit 42.11 → 0.091 in 15 iterations (4 s) |

  - Not available with LGR or dual porosity; those grids are reported as such.
  - Still open for full JutulDarcy parity:
    - a transposed block-CPR solve (now ILU(0)-GMRES);
    - adjoints for the compositional and thermal cores.

### 12.2 Multi-component product analysis (C9, `components.cpp`)

This applies to every multi-component run (EOS compositional, CO2STORE, thermal / combustion / UCG), not only thermal.

- **Rates:** the engines provide per-well component molar rates. The report-step average comes from the cumulative
  produced and injected moles (Σ rate·Δt reproduces the simulator's produced moles exactly, verified to 0). Field rates
  and cumulatives are built from these.
- **Species recognition** from CNAMES: H2, CO, CO2, C1/CH4 … C6, N2, O2, H2S, H2O/WATER/STEAM, NH3, AR.
- **Heating values:** standard 25 °C heats of combustion (HHV/LHV, kJ/mol): H2 285.83/241.82, CO 282.98, CH4
  890.36/802.34, … For C5+ pseudo components: LHV 44.0 MJ/kg, HHV 47.3 MJ/kg.
- **Stream analysis**, for the field or any well:
  - dry and gas mole fractions;
  - H2 (kmol/d, t/d, sm³/d);
  - syngas CO + H2 and the H2/CO ratio;
  - LHV/HHV power (MW) and gas heating value (MJ/sm³ at 15 °C, 1 atm);
  - CO2 and O2 slip (produced/injected);
  - per-well cumulative breakdown.
- **Results:** *Components* tab with component rates, dry composition, cumulatives, hydrogen/syngas and heating-value
  plots, plus the component table at the time cursor.

---

## 13. External engines

- **OPM Flow:** native Windows flow.exe. The app exports the deck, follows the live log, and imports EGRID / INIT /
  UNRST / SMSPEC / UNSMRY.
- **ECLIPSE:** runs through eclrun (or eclipse/e300/frontsim directly). THERMAL and COMPS decks select e300
  automatically.

---

## 14. Verification

- Test suite (`tests/run_tests.cpp`): see ROADMAP for the current count (C9 groups `-adj`, `-mc`, plus `-thm` / `-thermal`), passing native and on Windows. Highlights:
  - Analytic Buckley–Leverett.
  - Material-balance closure.
  - IF97 reference values.
  - Vinsome–Westerveld against the analytic half-space solution (within 5 %).
  - Energy balance of thermal floods < 1e-5.
  - Flash finite-difference checks (5e-7).
  - CO2SOL Jacobian check (1.4e-8).
- **Comparisons:** OPM Flow on the examples (relRMS 0.1–2 %), Norne (FPR within 0.3 % of ECLIPSE), and SPE1, SPE3, SPE5
  and SPE9.
- **Thermal Jacobian:** `Simulator::jacobianCheck(dt[, cell])` compares the assembled Jacobian against central
  differences.
- **C10 (`-sol`):**
  - IAPWS G7-04 check values;
  - general flash ≡ legacy flash (1e-13);
  - Spycher–Pruess CO2 solubility (≤ 5.4 %);
  - closed-cell dissolution: exact mass / energy, inventory closure, Jacobian 2.6e-7;
  - Land hysteresis in 1D;
  - UCG + CO2 round trip and run (energy balance 4e-7).
- **C11 (`-scw`):**
  - IF97 region 3 check values (6.9e-8) and jump-free isobars;
  - water flash continuity at T_c (|Δn_W| 7.8e-9);
  - a supercritical box heated from below (p_inj error −152 Pa, no cuts);
  - injector Jacobian at the rate and BHP limits (1.0e-3 / 9.4e-4);
  - 1D conduction against the erf solution (0.015 K of a 20 K step);
  - the **Theis** well test against E1 (worst 1.16 %, Δp(100 m, 1 d) = 2.079 bar).

  The planned analytic cases Lauwerier and Terzaghi were added in C12 (below). Mandel is still planned.
- **C9:** adjoint / forward-sensitivity gradient checks (§12.1), component mole conservation of the product analysis
  (§12.2), the steam-trap constraint (subcool held to 1e-5 K, Jacobian of the producer cell 8e-4) and oedometric
  surface-load settlement (5e-10).
- **C12** (groups `-ver`, `-sorb`, `-edfm`, `-spe4`, `-spe11`, `-adjf`):

  | Case | Reference | Result |
  |---|---|---|
  | **Lauwerier** (hot water in a thin layer, conductive cap / base rock) | analytic Laplace solution | \|ΔT_D\| < 0.05 along the layer |
  | **Terzaghi** 1-D consolidation (THM fixed-stress) | series solution | error 0.001 in the degree of consolidation |
  | **Langmuir** equilibrium and kinetics | isotherm; relaxation time τ | `-sorb` |
  | **ECBM** CO2 / CH4 competition | extended Langmuir (VL, PL of CO2 / CH4) | CH4 desorbs; mass closure |
  | **EDFM** single fracture | cubic law | 0.07 % |
  | **EDFM** influx to a fracture | analytic linear influx | 2.1 % |
  | Hydro-shear | Coulomb + dilation | 100 → 287 µm, irreversible |
  | EGS EDFM deck | — | export / import round trip and run: 129 fracture cells, 6 well–fracture connections, injected / produced 4.05e7 / 4.04e7 kmol |
  | **SPE4 problem 1A** (cyclic steam, 13 × 4 r-z, three 365-d cycles) | Aziz et al. (1987), text values | see below |
  | **SPE11A-like** slice (2.8 × 1.2 m, 20 °C, 1.1 bar; SPE11A facies 1/3/4/5 data, wells and schedule; layered geometry instead of the real facies map) | Carroll et al. (1991) solubility | see below |
  | **Egg model** (EclGym deck, 3600 d) | — | runs, see below |

  - **SPE4 1A details:**
    - OOIP 79,880 sm³ (paper 79,227–80,572) and OWIP 66,175 sm³ (paper 65,693–66,986).
    - Injected heat in cycle 1: 3,430 GJ against the paper's 3,594 GJ (IF97 datum). The 4.6 % shortfall is the BHP-limited
      start-up of the injector (§6.6); after 0.5 d the rate is full.
    - Heat loss after cycle 1: 616 GJ (paper 491–703 GJ).
    - Cumulative oil after three cycles 18.3 MSTB, water 26.8 MSTB. The paper only plots these (Figs. 2–3); the PDF text
      layer gives ≈ 18–21 MSTB, which is consistent but should not be read as a tabulated reference.
    - Run: 1,361 steps, 10,186 Newton iterations, 146 cuts, 18 s. That is many Newton iterations per step: the
      steam-front Newton is a speed target.
  - **SPE11A-like details:**
    - CO2 solubility at the top conditions is 0.0400 mol/kg (Henry / Carroll 0.0408).
    - The 4.590 g injection protocol is delivered exactly. CO2 mass balance closes to 3e-10.
    - Free gas stays below the seals (0 cm³ above the upper seal).
    - Dissolved fraction at 5 h: 76 % on 10 cm cells, 68 % on 5 cm cells. Dissolution depends on the grid: coarser cells mix more.
      At 5 d all free gas has dissolved.
    - The top boundary is a constant-pressure boundary: a BHP producer plus a BHP water injector. A producer alone closes
      the domain, and dissolution then drains the pressure below 1 bar.
    - A comparison with the SPE11 participants' maps needs the real facies geometry (the CSP `.geo` file); this is planned.
  - **Egg model details:**
    - FOPT 558,594 m³; FWPT 1.724e6 m³; FPR 399.3 bar; recovery ≈ 66 % of the ≈ 845k m³ OOIP.
    - Injected 2.29e6 m³ against 2.283e6 m³ produced (reservoir volumes).
    - 131 steps, 450 Newton iterations, 0 cuts.
    - Jansen et al. (2014) publish no numeric table, so this is a run check plus the §7 speed profile.
  - Injector drawdown-ramp FD check: moved to day 10, on a copy of the state. At day 30 the BHP-mode state sits within
    ~3 Pa of the limit, inside the FD step, and central differences straddle the C1 kink. Sweep: 2.3e-3 at ≥ 10 Pa from the
    limit, against a spurious 6.3e-2 at 3 Pa.

- **C13** (groups `-grin`, plus new checks in `-adj` and `-edfm`):
  - **Gringarten (1975) single-fracture heat extraction.** A vertical fracture of length L = 200 m and height H = 20 m in
    impermeable rock (λ = 2.5 W/m/K) carries ṁ = 0.5 kg/s of 320 K water into 470 K rock. The analytic outlet temperature
    is T_D = (T0 − T_out)/(T0 − T_inj) = erfc( λ·A / (ṁ·c_w·√(κ t)) ), with A = L·H (one face; both faces lose heat symmetrically) and
    κ = λ/(ρc)_r.

    | t (yr) | 0.25 | 0.5 | 1 | 2 | 5 | max \|ΔT_D\| |
    |---|---|---|---|---|---|---|
    | analytic | 0.020 | 0.099 | 0.243 | 0.409 | 0.602 | – |
    | EDFM, refined matrix (0.25 m first cell, ratio 1.5) | 0.046 | 0.130 | 0.264 | 0.419 | 0.603 | 0.031 |
    | EDFM, uniform 25 m matrix | 0.262 | 0.293 | 0.351 | 0.449 | 0.630 | 0.242 |

    The conduction boundary layer √(κt) is about 3 m after 0.25 yr. A 25 m matrix cell spreads the fracture's heat sink over the
    whole cell at once, so it over-predicts early cooling. The geometric refinement resolves the layer. Smaller time steps
    change nothing (the error is spatial). What remains is first-order upwinding along the plane, and it falls with
    the number of segments. The EGS designer therefore refines around every fracture plane by default (2 m first cell, ratio 1.8, a compromise for field-size models). The multistage
    example then gives production temperatures of 162.7 / 121.9 / 75.4 °C at 91 d / 1 yr / 10 yr. The uniform grid
    gave 128.4 / 118.5 / 78.9.
  - **Optimiser records:** L-BFGS stores every accepted iterate (all parameter values, the gradient, its run, cumulative
    runs and wall time; \|g\| 53 → 4.8e-6 on the `-adj` case). The field HM stores the ln k / ln φ snapshots per iterate
    (the last snapshot equals the final multipliers to 5.7e-8) and the sensitivity maps at the prior and at the optimum.
  - **EGS designer:** transverse fractures connect both horizontal wells. The produced temperature falls monotonically
    between T_inj and T_res. Net thermal power is positive and declining (3 stages, 1 yr: 131.4 → 96.1 °C,
    10.35 → 6.08 MW).

- **C16 verifications**
  - **Plasticity (`-gmp`)**: platen compression reaches the Mohr–Coulomb strength 45.0672 MPa (error 1.2e-8); with softening the
    peak is 36.82 MPa (inner cone), falling to the residual 35.23 → 39.87 MPa plateau; Mandel pressure error ≤ 0.0021 (Mandel–Cryer
    rise reproduced), Cryer sphere peak 1.0850 vs 1.0859; SAGD THM 91 d: ε_p 8.72e-4, k/k0 1.299, energy balance 3.8e-9.
  - **Mixed precision / GPU / MPI (`-mp`)**: BILU0 and BILU0f both 18 iterations, refined residual 2.2e-13; demo FOPT 5.47875e6
    identical for double, mixed and GPU paths (53 device solves, no fallback, software path in the sandbox); MPI batch FOPT
    2.73938e6 equal to the direct run (6 decks over 3 ranks in 8 s).
  - **SPE4 problem 2A (`-spe4b`)**: steamflood of the problem-1 oil in the 1/8 inverted nine-spot, the specified 9 × 5 × 4 point-
    distributed grid along the diagonal (25 points per layer, area = pattern/8 exactly), Peaceman r_e = 0.2 Δ, well fractions 1/8,
    1/4, 1/8, Yanosik–McCracken nine-point stencil. Time to 10 STB/d steam: **near 370 d** (participants 257–539), **far 1,645 d**
    (1,398–2,245). The five-point stencil gives the reversed order (far 735 d before near 970 d) — the grid-orientation effect the
    paper discusses. OOIP 251,792 STB per pattern; Np 188 MSTB in 10 yr; 1,168 steps, 13 s. The 10 STB/d steam limit is a lagged
    back-pressure control on 5-day reports (steam ≤ 1.5 × limit in ≈ 67 % of the post-breakthrough reports) — not an implicit constraint.
  - **SPE11A official geometry (`-spe11g`)**: the CSP `spe11a.geo` rasterised (tools/spe11_raster.py, point-in-polygon of the cell
    centres; maps at 56 × 24, 140 × 60, 280 × 120), facies 1–6 of Table 1 with Brooks–Corey Pc capped at 9.5e4 Pa, facies 7
    inactive, both wells and the protocol, Dirichlet top. 56 × 24: injected 4.590 g, mass balance 5e-8; max mobile CO2 in Box A at
    18,000 s (FluidFlower 14,800; forecasts 17,700–18,000); Box A at 72 h dissolved 2.31 g (exp 3.10, forecasts 2.11–3.60), seal
    0.14 g (0.28–0.57), mobile 0 (0–0.87); Box B dissolved 0.03 g (0.30–0.78) — the 5 cm grid dissolves too fast and the plume barely
    reaches Box B. **C17: with the specified dispersivity 1e-2 m (DISPERC model, §6.10)** — max mobile in Box A at 17,400 s,
    Box A dissolved 2.30 g, seal 0.22 g (closer to 0.28–0.57), Box B dissolved 0.008 g, mass balance 7e-9 (1539 steps, 150 s).
    `SPE11_DISP=0` reproduces the C16 numbers. Finer maps are provided; the 140 × 60 run did not finish within the sandbox session.
  - **Pruess (2006) CO2-EGS five-spot (`-egs`)**: 1/8 element with point-distributed control volumes (21 × 21 diagonal points),
    hand-built MINC (cubic blocks D = 50 m, fracture 2 % with φ 0.5, five matrix shells by the proximity function
    1 − (1 − 2x/D)³), rock and operating data of Table 3, wells Peaceman r_w = 0.1 m (not given in the paper). CO2 / water: initial
    heat **1.46** (paper ≈ 1.5), 25-yr mean 1.77, initial mass flow **3.46** (paper ≈ 3.7). Absolute rates (CO2 34.6 MW / 110 kg/s
    vs ≈ 75 MW / 280 kg/s) depend on the unspecified well representation (r_w = 3 m gives +45 %) and are not claimed. Grid check:
    the 11 × 11 element gives ratios 1.25 / 2.97 (coarse cold zone around the injector).
  - **SPE3 compositional**: still deferred — the published EOS characterisation is not available to us; the black-oil SPE3 case
    remains (FOPR within 0.09 %).
  - **C17 checks**: `-tminc` thermal and `-cminc` compositional MINC (8 shells within 0.006 of the fine 1D model, Warren–Root
    0.14), layered conduction geometry exact, `-cdisp` dispersion (α recovered within 4–7 %, DISPERC round trip), `-dlgr` LGR
    completions in compositional / thermal decks, `-hf` multi-cluster stress shadowing (Nordgren limits 1.2 % / 0.1 %), `-spe4b`
    implicit steam limit (100 % of reports within 1.02 × limit). Native suite 407 PASS, Windows (wine) 406 PASS.
  - **C18 checks**: `-cadj` compositional (7.5e-6 vs FD) and thermal (≤ 4e-4) adjoint gradients; `-mp` multicolour BILU0 and
    float refinement; `-mpidd` distributed CPR solve (1 / 2 ranks = serial solution, demo FOPT identical); `-hfg` in-run LEFM
    fracture growth (initial opening, symmetric propagation, stress dependence, pressure pinned at σ_min + K_IC/√(πH/2), exact
    restart, HFOPEN); results-file v1/v2 round trip. Native suite 433 PASS, Windows (wine) 432 PASS.
  - **C32 checks** (§6.25): `-c32msw` multisegment wells (friction / hydrostatic / valve = hand values, frictionless
    tubing = standard well, SICD equalisation and later breakthrough, SHUT subtree, exact restart, WSEGAICD items);
    `-c32geo` geomodelling (exchange formats, SEG-Y IBM / IEEE, minimum curvature, auto-tracking RMS 2.4 m, fault throw
    27.5 / 30 m, well tie 2.32 → 0.04 m, facies / SGS / co-simulation, volumetrics = hand sum, flow-based upscaling exact
    for layers, PV conserved, geomodel → engine run); `-c32actx` ACTIONX COMPDAT at run time = the same COMPDAT in the
    schedule (5e-12), exact restart.
  - **C40 checks** (§6.42): `-c40geo` – log calculator expressions / null propagation, log-linear regression (R² = 1 on
    exact data), variogram fit (relative error 1.7·10⁻⁸ on exact model values; AR(1) field → range 46 vs 40), well tie on the
    demo (correlation ≥ 0.95, |Δ| ≤ 2 ms) and with a wrong velocity model (+30 ms found, TWT error 0.3 ms), geomodel Study
    parameters (seed reproducibility, range → correlation 0.05 → 0.72, PROPMODEL round trip).

---

## 15. Open items and candidate advanced methods

Ideas that could improve robustness and speed further. The user may also suggest more.

| Method | Where | Expected gain |
|---|---|---|
| ~~UDQ / UDA / ACTIONX in the own engine~~ | field management | **done in C31 (§6.24)**; **multisegment wells and connections added by an action at run time done in C32 (§6.25)**; open: PYACTION, wells created by an action, MSW acceleration / drift-flux terms and ACTIONX connections in MSW / DD runs |
| ~~Geomodelling (data → seismic → structure → grid → facies → petrophysics → volumetrics → upscaling)~~ | geomodel / ResGeo | **foundation done in C32 (§6.25)**; next: stratigraphic (non-vertical) pillars and faults as pillars, horizon interpretation editing, well correlation flattening, MPS / object libraries, cloud transform, uncertainty workflows across structure + properties, seismic inversion attributes |
| Analytic local derivatives (AD) for the thermal (p, N, U) flash via the implicit-function theorem | thermal | Jacobian exact everywhere; no finite-difference noise near univariant states |
| Nonlinear preconditioning: ASPEN / multiplicative Schwarz on the hard (phase-changing) cells | thermal, compositional | 2–5× fewer global Newton iterations on steam fronts |
| Anderson acceleration of the energy / flash fixed points | thermal, CO2 brine | fewer local iterations |
| Trust-region Newton on fractional flow (Jenny–Tchelepi / Wang–Tchelepi) | black oil, thermal | unconditional Newton convergence for large Δt |
| Sequential fully implicit (SFI) with an outer loop | all | memory and speed for large models |
| Multiscale finite-volume pressure (MsFV / MsRSB) as the CPR first stage | all | million-cell pressure solves |
| GPU BSR SpMV and BILU0 (D3D11 compute / CUDA) | linear algebra | block-Jacobi BiCGSTAB **C16**; **multicolour BILU0 on the device + device/wall timing C18 (§6.11)**; **CPR-AMG on the device C29 (§6.22)**; next: CUDA/HIP backend, numbers from a real discrete GPU |
| ~~Adjoint gradients~~ | history matching | **done in C9 (§12.1)**; field gradients **C12**; **compositional and thermal adjoint C18 (§6.11)**; open: transposed CPR solve, adjoint through growing fractures |
| ~~Geomechanics (poro-elastic, fixed-stress split)~~ | coupled | **done in C8 (§6.2)**; fracture dilation permeability done in C12 (EDFM); **C16: Drucker–Prager plasticity with softening and damage permeability (§6.9)**; next: AMG for elasticity |
| Multigrid (AMG) preconditioner for the elasticity system | geomechanics | 5–10× fewer CG iterations than SGS on large boxes |
| ~~Continuous supercritical water~~ | thermal | **done in C11** (IF97 region 3, C1 K_w, isochore cap); next: real-mixture EOS for water + gases, or an equal-area fold construction in the energy flash (near-critical CO2–water cavities) |
| ~~Adjoint with analytic ∂R/∂θ~~ (per-face T, WI, PERM / PORO fields) | history matching | **done in C12 (§12.1)**; the transposed CPR solve is still open |
| ~~Competitive extended-Langmuir sorption~~ | THMC, CO2 in UCG / ECBM | **done in C12 (§6.6)** |
| From the C27 CMG / TOUGH3 comparison: ~~EWASG salt~~, ~~wellbore friction / ICD~~, ~~surface network~~ **done in C29 (§6.22)** (explicit / lagged forms); implicit network balance **done in C30 (§6.23)**; open: transient discretised wellbore (full FlexWell / T2Well), chokes / compressors, dynamic grid refinement (DynaGrid, deferred: design note in ROADMAP) | black oil, wells | see `COMPARISON_CMG_TOUGH3.md` §4 |
| ~~EDFM fracture cells + Barton–Bandis / Coulomb hydro-shear~~ | EGS | **done in C12 (§6.6)**; next: implicit (in-Newton) aperture, inclined fractures, corner-point hosts, propagation |
| Mixed-precision CPR (float BILU0 / AMG data), multicolour smoothers | linear algebra | float BILU0 **done in C16**; next: float AMG hierarchy, multicolour smoothers |
| ~~Mandel consolidation check~~ | geomechanics verification | **done in C16** (Mandel + Cryer) |
| ~~MINC (multiple interacting continua)~~ | dual porosity, EGS | **done in C14 (§6.7)** for the black-oil core; **thermal and compositional in C17 (§6.10)** |
| Fracture propagation phase 2 | EGS, unconventional | PKN / KGD sizing **C14**; multi-cluster stress shadow **C17**; **in-run LEFM growth of EDFM planes C18 (§6.11)**; P3D height growth **C20 (§6.12)**; embedded contact mechanics G1 **C20**; next: EFEM↔EDFM coupling, 3-D, DDM HF, curving paths |
| Elastoplastic geomechanics + continuum damage | geomechanics | Mohr–Coulomb / Drucker–Prager return mapping on the existing hexahedral FE grid, scalar damage D with E(1 − D) and k = k0 exp(β D); chosen over XFEM (no enrichment bookkeeping, works for corner-point hosts); checks: Mandel, elastoplastic wellbore (Salençon) |
| Single-run MPI domain decomposition, phase 2 | HPC | **phase 1 C18**: distributed CPR linear solve (assembly on rank 0); next: distributed assembly and AMG, METIS-style partitions, compositional / thermal |
| XFEM / embedded-discontinuity geomechanics on the EDFM planes | geomechanics | Heaviside-enriched hexahedra cut by the HFRAC planes: aperture from the displacement jump, slip with Coulomb friction, tip stress-intensity from the interaction integral; replaces the empirical aperture laws (phased plan in ROADMAP C18) |
| Higher-order (TVD) upwinding along EDFM planes | EGS | removes the remaining ≈ 0.03 Gringarten error without more segments |
| Real SPE11A geometry (CSP .geo) and SPE11B/C | CO2 verification | SPE11A geometry **done in C16** (coarse); isotropic dispersion (DISPERC) **C17**; next: full tensor (α_L ≠ α_T), 1 cm runs, SPE11B/C |
| Nanoparticle transport with retention / release and wettability alteration | EOR | nanofluid and fines-migration studies |
| Detailed coal / oil kinetics libraries (Belgrave, Kapadia UCG) with calibration to tube runs | THMC | quantitative combustion and UCG forecasts |
