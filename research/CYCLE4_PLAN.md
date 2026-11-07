# Research Cycle 4 — pre-registration

Date: 2026-10-20 (written before any Cycle-4 experiment). Charter: `CHARTER.md`. Inputs: `CYCLE3_REPORT.md` (candidates F10, F11).

## H1f — a closed form for the schedule gain G(Π, θ) (F10 → theory)

Model (as Cycles 1–3): Dieterich (1994) rate-and-state seismicity, linear (background-dominated, t_a ≫ T) regime, so the expected
excess count is I = ∫ (e^{Π s(t)} − 1) dt with s the stress normalised by the constant-rate maximum, Π = ΔS_max/(Aσ). Single
exponential rate → stress kernel: T ṡ = −s + u, u = q/q̄ ∈ [0, u_max], ∫_0^{T_i} u dt = T_i, θ = T/T_i. Observation window: injection
plus the full decay tail.

Derivation made before the numerics (Pontryagin): H is linear in u, so the optimum is bang–singular: a maximum-rate ramp, then a
singular arc with ∂H/∂u ≡ 0 ⇒ λ̇ = 0 ⇒ e^{Πs} constant ⇒ **constant stress plateau** with u = s*. An early stop is never optimal
(its first-order cost t_f (Π e^Π − e^Π + 1) > 0). Expanding both schedules to first order in θ gives

    G(Π, θ, u_max) = θ · κ(Π, u_max) + O(θ²),
    κ = e^Π [ Π m(u) + R(Π, u) − Ein(Π) ] / (e^Π − 1),
    m(u) = (u − 1) ln(u/(u − 1)),   R = ∫_0^{ln(u/(u−1))} (1 − e^{Π(u(1 − e^{−x}) − 1)}) dx,
    Ein(Π) = γ + ln Π + E1(Π)  (entire exponential integral);   u → ∞: κ = e^Π (Π − Ein(Π)) / (e^Π − 1).

Pre-registered tests (numerical optimum = discretised convex programme, exact exponential stepping, ≥ 400 steps):

| Test | Pass criterion |
|---|---|
| F1 structure | the numerical optimum is ramp-at-u_max + constant-stress plateau: stress on the plateau varies < 2 % of its mean, no early stop; Π ∈ {1, 3, 6}, θ = 0.05 |
| F2 small-θ law | numerical G within ± 15 % of θ κ at θ ≤ 0.02 for Π ∈ {1, 3, 6}, u_max ∈ {3, 10} |
| F3 interior maximum | G(θ) has an interior maximum for every Π tested (window with tail) — report θ*(Π) and G_max(Π); descriptive |
| F4 Otaniemi consistency | the law with the Otaniemi Π = 3.06 and u_max reproduces the small-θ slope of E8 within a factor 2; if not, identify the missing ingredient |

A FAIL of F1 or F2 refutes the derivation; F3/F4 are descriptive.

## H2b — sub-grid leak-off for HF growth (F11 follow-up)

Pre-registered pass (from CYCLE3_REPORT): spread of t60 over dx 10 / 5 / 2.5 m below 1 d with a Carter-type fresh-rock leak-off for the
tip host cell. Run only if the engine change is small and the runs fit this cycle; otherwise deferred with reason.

## H5 — coupled G1-c in the seismicity designer

Deferred unless H1f/H2b finish early (engineering-heavy; no new hypothesis yet).

## Second field site for F10 (Basel 2006)

Open data search first (injection log and catalogue licences); if no openly licensed hourly injection log is found, F10 stays CANDIDATE
and the closed form (H1f) is the cycle's contribution.
