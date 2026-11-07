# C35-3 SOLVENT design note (developer)

## Why the solvent cannot be transported sequentially like the chemicals

Solvent occupies pore volume. Suppose the FIM solves p, Sw and X with the solvent saturation Ss lagged. A Gauss-Seidel
analysis at an injector (incompressible limit) gives

  e^{k+1} = -(λs / λ(o+w)) e^k

so the outer loop diverges whenever the solvent is more mobile than the liquids, which is the normal case. The pressure
solve must therefore see the solvent volume.

## Chosen formulation: FIM on the total gas-like component with a lagged composition, plus composition transport

* **FIM equations:** water, oil, and G+S (the gas component plus the solvent, both in surface volume).
* **FIM unknowns:** p, Sw, and X = SgT (the total gas-phase saturation Sg + Ss) in state 1, or Rs in state 0. State 2
  (no oil) keeps X = Rv.
* **Lagged composition:** z = Ns / (Ns + Ng), the solvent fraction of the cell's G+S amount. Ng includes the
  dissolved gas.
* **Flash in evalCell (AD), state 1:**
  * Free-gas branch: Rs = Rs_sat(p), So = 1 - Sw - SgT, and
    Ss = z (SgT bg + So bo Rs) / (bs (1 - z) + z bg).
    The branch is valid if Ss <= SgT; then Sg = SgT - Ss.
  * Otherwise, pure-solvent gas phase: Sg = 0, Ss = SgT, and Rs = Ss bs (1 - z) / (z So bo(p, Rs)), found by fixed
    point. Rs is capped at Rs_sat, where the branches meet continuously.
  * With z = 0 the flash reproduces the no-solvent model exactly.
* **Flash, state 2:** Ss = z SgT bg / (bs (1 - z) + z bg) with SgT = 1 - Sw.
* **Saturation functions use SgT:** krg(SgT), kro(Sw, SgT), pcog(SgT).
  * Immiscible: mob_g = bg/μg · krg(SgT) · SSFN_g(Sg/SgT), and mob_s = bs/μs · krg(SgT) · SSFN_s(Ss/SgT).
  * Miscible (MISC): M = MISC(Fs) · PMISC(p).
    * sor = M·SORWMIS(Sw) + (1 - M)·sogcr, and sgc = M·SGCWMIS(Sw) + (1 - M)·sgcr.
    * F = (SgT - sgc) / (So + SgT - sor - sgc).
    * mkrgt = MSFN_krsg(F)·SOF2(Sn) and mkro = MSFN_kro(F)·SOF2(Sn).
    * kro = (1 - M) kro + M mkro, and krgT = (1 - M) krgT + M mkrgt, then the SSFN split.
    * Todd-Longstaff effective viscosities and densities with TLMIXPAR (× TLPMIXPA(p)), as in OPM
      blackoilsolventmodules.hh, effectiveProperties().
* **Accumulation:** acc[G+S] = pv (bg Sg + bo Rs So + bs Ss).
* **Flux:** the gas equation also gets the solvent flux mob_s ΔΦ_s, where ΔΦ_s uses the gas-phase pressure and the
  solvent density. It is added in both flux assembly paths (simulator.cpp ~1938 and ~2015).
* **Wells:**
  * Producer: the G+S equation outflow includes wi·mob_s·dd. Reported gas excludes the solvent, which is reported
    separately (FNPR / WNPR).
  * Gas injector with WSOLVENT fraction f: the total surface rate q drives the G+S equation. The reservoir volume uses
    1/b_mix = (1 - f)/bg + f/bs. The solvent part f·q goes to the transport.
* **Composition transport after the FIM converges:** keep p, Sw and the FIM state X fixed. Unknown z per cell (Newton,
  ILU0-BiCGSTAB).
  * Solvent mass pv·Ss(z)·bs.
  * Upwind fluxes T·mob_s(z_up)·ΔΦ_s with frozen potentials.
  * Well rates as above.
  * Then the outer loop repeats the FIM until z stops changing (the same loop as polymer / chem). This is the standard
    sequential-fully-implicit compositional split (Møyner & Lie). The composition coupling is weak, so it converges.
* **State changes in the transport:** a state-0 cell receiving solvent switches to state 1 with SgT = Ns / (pv bs).
  In the FIM update, SgT < 0 with z ≈ 0 switches to state 0.
* **Reporting:**
  * SGAS = Sg, SSOLVENT = Ss, SOIL = 1 - Sw - SgT.
  * FNPR / FNPT / FNIR / FNIT, WNPR / WNIR, FNIP.
* **Restart:** state v9 adds z.
* **Deck keywords:**
  * RUNSPEC SOLVENT and MISCIBLE.
  * PROPS PVDS, SDENSITY, SSFN, MISC, PMISC, MSFN, SORWMIS, SGCWMIS, TLMIXPAR, TLPMIXPA.
  * REGIONS MISCNUM.
  * SOLUTION SSOL.
  * SCHEDULE WSOLVENT.
  * FOAMOPTS SOLVENT makes the solvent the foam carrier.
  * Not supported (note): CO2SOL / H2SOL.
* **Not supported with solvent (error / note):** GPU assembly, MPI DD, adjoint.
* **Validation targets (opm-tests, with OPM references):**
  * SPE1CASE2_SOLVENT_FOAM: immiscible solvent plus foam with the solvent carrier.
  * SPE5CASE_SOLVENT_EHYSTR2_1: miscible.
  * GRUP_SOLVENT: CO2SOL, so out of scope.
