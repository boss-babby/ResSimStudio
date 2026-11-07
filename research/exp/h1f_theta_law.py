#!/usr/bin/env python3
"""Research Cycle 4, H1f: closed form of the schedule gain G(Pi, theta) for one exponential stress kernel.

Model (CYCLE4_PLAN.md): excess count I = int (exp(Pi s/s_c) - 1) dt (Dieterich, background-dominated regime),
T ds/dt = -s + u, u in [0, u_max] during injection [0, Ti], int u dt = Ti, u = 0 afterwards (tail fully counted).
s_c = max stress of the constant-rate schedule, theta = T/Ti.  Ti = 1.

Outputs exp/h1f_results.json and exp/h1f_fig.png.
"""
import json, math, sys
import numpy as np
from scipy.optimize import minimize, brentq
from scipy.integrate import quad
from scipy.special import exp1

EULER = 0.5772156649015329


def Ein(x):
    return EULER + math.log(x) + exp1(x)


# ------------------------------------------------------------------ closed forms
def kappa(Pi, u):
    """first-order coefficient: G = theta * kappa + O(theta^2)"""
    if math.isinf(u):
        m, R = 1.0, 0.0
    else:
        xr = math.log(u / (u - 1.0))
        m = (u - 1.0) * xr
        R = quad(lambda x: 1.0 - math.exp(Pi * (u * (1.0 - math.exp(-x)) - 1.0)), 0.0, xr)[0]
    return math.exp(Pi) * (Pi * m + R - Ein(Pi)) / (math.exp(Pi) - 1.0)


def tail_int(Pi_eff, T):
    """int_0^inf (exp(Pi_eff e^{-t/T}) - 1) dt = T (Ei(Pi_eff) - gamma - ln Pi_eff)"""
    if Pi_eff <= 0:
        return 0.0
    from scipy.special import expi
    return T * (expi(Pi_eff) - EULER - math.log(Pi_eff))


def cost_constant(Pi, theta):
    T = theta
    sc = 1.0 - math.exp(-1.0 / theta)
    f = lambda t: math.exp(Pi * (1.0 - math.exp(-t / T)) / sc) - 1.0
    I_inj = quad(f, 0.0, 1.0, limit=400, points=[min(1.0, 5 * T)])[0]
    return I_inj + tail_int(Pi, T)     # s(Ti)/sc = 1


def cost_ramp_plateau(Pi, theta, u):
    """exact cost of the Pontryagin candidate (ramp at u until s = s*, then s = s* to Ti, then decay), any theta"""
    T = theta
    sc = 1.0 - math.exp(-1.0 / theta)

    def vol_err(ss):
        tr = T * math.log(u / (u - ss))
        return u * tr + ss * (1.0 - tr) - 1.0

    hi = min(u * (1 - 1e-12), 50.0)
    # feasible plateau needs t_r <= Ti
    smax_ramp = u * (1.0 - math.exp(-1.0 / T))
    hi = min(hi, smax_ramp)
    if vol_err(hi) < 0:
        return float("nan"), float("nan"), float("nan")
    ss = brentq(vol_err, 1e-12, hi, xtol=1e-14)
    tr = T * math.log(u / (u - ss))
    ramp = quad(lambda t: math.exp(Pi * u * (1.0 - math.exp(-t / T)) / sc) - 1.0, 0.0, tr, limit=200)[0]
    plat = (1.0 - tr) * (math.exp(Pi * ss / sc) - 1.0)
    return ramp + plat + tail_int(Pi * ss / sc, T), ss / sc, tr


def tailP(P):
    from scipy.special import expi
    return expi(P) - EULER - math.log(P) if P > 0 else 0.0


def burst_B(Pi, u):
    """first-order benefit of a terminal max-rate burst (units of T e^Pi); returns (B*, beta*)"""
    from scipy.optimize import minimize_scalar
    if math.isinf(u):   # impulse of normalised size w = v/T
        f = lambda w: -(Pi * w - math.exp(-Pi) * (tailP(Pi * (1 + w)) - tailP(Pi)))
        r = minimize_scalar(f, bounds=(0.0, 5.0), method="bounded", options={"xatol": 1e-10})
        return -r.fun, r.x

    def B(beta):
        if beta <= 0:
            return 0.0
        se = u - (u - 1.0) * math.exp(-beta)
        burst = quad(lambda x: math.exp(Pi * (u - (u - 1.0) * math.exp(-x)) - Pi) - 1.0, 0.0, beta)[0]
        return Pi * beta * (u - 1.0) - burst - math.exp(-Pi) * (tailP(Pi * se) - tailP(Pi))
    r = minimize_scalar(lambda b: -B(b), bounds=(0.0, 3.0), method="bounded", options={"xatol": 1e-10})
    return -r.fun, r.x


def kappa2(Pi, u):
    """revised first-order coefficient (bang-singular-bang: ramp, plateau, terminal burst)"""
    k1 = kappa(Pi, u) * (math.exp(Pi) - 1.0) / math.exp(Pi)
    return math.exp(Pi) * (k1 + burst_B(Pi, u)[0]) / (math.exp(Pi) - 1.0)


def cost_three_arc(Pi, theta, u):
    """exact cost of ramp (u) + plateau (s*) + terminal burst (u, length tb) at any theta, minimised over tb"""
    from scipy.optimize import minimize_scalar
    T = theta
    sc = 1.0 - math.exp(-1.0 / theta)

    def cost(tb):
        def vol_err(ss):
            tr = T * math.log(u / (u - ss))
            return u * tr + ss * (1.0 - tr - tb) + u * tb - 1.0
        hi = min(u * (1 - 1e-12), u * (1.0 - math.exp(-(1.0 - tb) / T)))
        if vol_err(1e-12) > 0 or vol_err(hi) < 0:
            return float("inf"), 0, 0
        ss = brentq(vol_err, 1e-12, hi, xtol=1e-14)
        tr = T * math.log(u / (u - ss))
        if tr + tb > 1.0:
            return float("inf"), 0, 0
        ramp = quad(lambda t: math.exp(Pi * u * (1.0 - math.exp(-t / T)) / sc) - 1.0, 0.0, tr, limit=200)[0]
        plat = (1.0 - tr - tb) * (math.exp(Pi * ss / sc) - 1.0)
        burst = quad(lambda t: math.exp(Pi * (u - (u - ss) * math.exp(-t / T)) / sc) - 1.0, 0.0, tb, limit=200)[0] if tb > 0 else 0.0
        se = u - (u - ss) * math.exp(-tb / T)
        return ramp + plat + burst + tail_int(Pi * se / sc, T), ss / sc, se / sc
    r = minimize_scalar(lambda b: cost(b)[0], bounds=(0.0, min(0.5, 3 * T)), method="bounded", options={"xatol": 1e-9 * max(T, 1e-3)})
    c0 = cost(0.0)
    best = (r.x, *cost(r.x)) if cost(r.x)[0] <= c0[0] else (0.0, *c0)
    return best[1], best[2], best[3], best[0]    # cost, s*/sc, s_end/sc, tb


# ------------------------------------------------------------------ discretised convex optimum
def grid(theta):
    """injection grid refined where the optimum has structure: the start ramp (~0.1 T) and the terminal burst (~0.05 T)"""
    T = theta
    pts = [0.0]
    def seg(a, b, n):
        a, b = max(a, pts[-1]), min(b, 1.0)
        if b > a + 1e-15:
            pts.extend(np.linspace(a, b, max(1, n) + 1)[1:].tolist())
    seg(0.0, min(T, 0.3), 60)
    seg(min(T, 0.3), min(12 * T, 0.6), 60)
    seg(min(12 * T, 0.6), 1.0 - min(2 * T, 0.3), 100)
    seg(1.0 - min(2 * T, 0.3), 1.0 - min(0.4 * T, 0.1), 30)
    seg(1.0 - min(0.4 * T, 0.1), 1.0, 60)
    tin = np.array(pts)
    tail = 1.0 + T * np.concatenate([np.linspace(0, 1, 41)[1:], np.linspace(1, 6, 61)[1:], np.linspace(6, 40, 51)[1:]])
    return tin, np.concatenate([tin, tail])


class Problem:
    def __init__(self, Pi, theta, u_max):
        self.Pi, self.T, self.u = Pi, theta, u_max
        self.tin, self.t = grid(theta)
        self.dt = np.diff(self.t)
        self.nin = len(self.tin) - 1
        self.a = np.exp(-self.dt / self.T)
        self.sc = 1.0 - math.exp(-1.0 / theta)
        # within-step quadrature: 3-point Gauss on each step using exact s(t) = u + (s0 - u) e^{-tau/T}
        g = np.array([-math.sqrt(0.6), 0.0, math.sqrt(0.6)])
        self.gx = (g + 1) / 2
        self.gw = np.array([5, 8, 5]) / 18.0
        self.e = np.exp(-np.outer(self.dt, self.gx) / self.T)          # (nsteps, 3)

    def forward(self, uin):
        n = len(self.dt)
        uu = np.zeros(n); uu[:self.nin] = uin
        s0 = np.zeros(n + 1)
        for k in range(n):
            s0[k + 1] = uu[k] + (s0[k] - uu[k]) * self.a[k]
        return uu, s0

    def cost_grad(self, uin):
        Pi, sc = self.Pi, self.sc
        uu, s0 = self.forward(uin)
        sq = uu[:, None] + (s0[:-1, None] - uu[:, None]) * self.e           # stress at quadrature points
        ex = np.exp(Pi * sq / sc)
        J = float(np.sum(self.dt[:, None] * self.gw * (ex - 1.0)))
        # adjoint: dJ/ds0[k] and dJ/du[k]
        dq = self.dt[:, None] * self.gw * ex * Pi / sc                       # dJ/dsq
        n = len(self.dt)
        lam = 0.0                                                            # dJ/ds0[k+1] accumulated
        gu = np.zeros(n)
        for k in range(n - 1, -1, -1):
            # s0[k+1] = u + (s0[k] - u) a ; sq[k] = u + (s0[k] - u) e
            ds0k = float(np.sum(dq[k] * self.e[k])) + lam * self.a[k]
            gu[k] = float(np.sum(dq[k] * (1.0 - self.e[k]))) + lam * (1.0 - self.a[k])
            lam = ds0k
        return J, gu[:self.nin]

    def solve(self, x0=None):
        dtin = self.dt[:self.nin]
        if x0 is None:
            x0 = np.ones(self.nin)
        cons = [{"type": "eq", "fun": lambda x: float(np.dot(x, dtin)) - 1.0, "jac": lambda x: dtin}]
        J0 = self.cost_grad(np.ones(self.nin))[0]
        f = lambda x: tuple(v / J0 if i == 0 else v / J0 for i, v in enumerate(self.cost_grad(x)))
        r = minimize(f, x0, jac=True, method="SLSQP", bounds=[(0.0, self.u)] * self.nin, constraints=cons,
                     options={"maxiter": 1500, "ftol": 1e-12})
        return r.x, r.fun * J0, r, J0


def numeric_gain(Pi, theta, u):
    P = Problem(Pi, theta, u)
    x, Jopt, r, Jc = P.solve()
    return 1.0 - Jopt / Jc, x, P, r


def main():
    out = {"model": "I = int (exp(Pi s/s_c) - 1) dt; T s' = -s + u; tail counted", "tests": {}}
    # ---- F1 structure at theta = 0.05
    f1 = []
    for Pi in (1.0, 3.0, 6.0):
        G, x, P, r = numeric_gain(Pi, 0.05, 10.0)
        uu, s0 = P.forward(x)
        s = s0[1:P.nin + 1] / P.sc
        ramp_end = int(np.argmax(x < 0.98 * P.u))
        burst_start = P.nin - int(np.argmax(x[::-1] < 0.98 * P.u))
        plat = s[ramp_end + 3:max(ramp_end + 4, burst_start - 3)]
        var = float((plat.max() - plat.min()) / plat.mean())
        early_stop = bool(np.any(x[int(0.9 * P.nin):] < 0.5 * plat.mean()))
        Jc = cost_constant(Pi, 0.05)
        J2 = cost_ramp_plateau(Pi, 0.05, 10.0)[0]
        J3, ss3, se3, tb3 = cost_three_arc(Pi, 0.05, 10.0)
        tail_burst_vol = float(np.dot(np.maximum(x[burst_start - 2:] - plat.mean(), 0), P.dt[burst_start - 2:P.nin]))
        f1.append(dict(Pi=Pi, G_numeric=G, G_two_arc=1 - J2 / Jc, G_three_arc=1 - J3 / Jc, plateau_var=var, early_stop=early_stop,
                       terminal_burst_volume=tail_burst_vol, tb_three_arc=tb3, s_end_three_arc=se3, s_plateau_three_arc=ss3, ok=bool(r.success)))
        print(f"F1 Pi {Pi}: G numeric {100*G:.3f}% | 2-arc (pre-reg) {100*(1-J2/Jc):.3f}% | 3-arc {100*(1-J3/Jc):.3f}% (tb {tb3:.4f}, s* {ss3:.4f}, "
              f"s_end {se3:.3f}) | plateau var {100*var:.2f}% | early stop {early_stop} | terminal burst vol {tail_burst_vol:.4f}", flush=True)
    out["tests"]["F1"] = dict(rows=f1, passed_as_registered=all(r["plateau_var"] < 0.02 and not r["early_stop"] for r in f1),
                              note="structure found: bang-singular-bang (terminal max-rate burst), not the pre-registered ramp+plateau")
    # ---- F2: pre-registered law (kappa) and revised law (kappa2) vs the numerical optimum
    f2 = []
    for Pi, tag in ((1.0, "reg"), (3.0, "reg"), (6.0, "reg"), (2.0, "fresh"), (4.5, "fresh"), (8.0, "fresh")):
        for u in (3.0, 10.0):
            for th in ((0.01, 0.02) if tag == "reg" else (0.015,)):
                Jc = cost_constant(Pi, th)
                G3 = 1 - cost_three_arc(Pi, th, u)[0] / Jc
                Gn = numeric_gain(Pi, th, u)[0]
                l1, l2 = th * kappa(Pi, u), th * kappa2(Pi, u)
                f2.append(dict(Pi=Pi, u=u, theta=th, set=tag, G_numeric=Gn, G_three_arc=G3, law_prereg=l1, law_revised=l2,
                               ratio_prereg=Gn / l1, ratio_revised=Gn / l2))
                print(f"F2 [{tag}] Pi {Pi} u {u:4.0f} th {th:.3f}: numeric {100*Gn:.3f}%  3-arc {100*G3:.3f}%  law(pre-reg) {100*l1:.3f}% "
                      f"(ratio {Gn/l1:.3f})  law(revised) {100*l2:.3f}% (ratio {Gn/l2:.3f})", flush=True)
    out["tests"]["F2"] = dict(rows=f2, passed_prereg=all(abs(r["ratio_prereg"] - 1) <= 0.15 for r in f2 if r["set"] == "reg"),
                              passed_revised_fresh=all(abs(r["ratio_revised"] - 1) <= 0.15 for r in f2 if r["set"] == "fresh"))
    # ---- F3 interior maximum (3-arc closed form, spot-checked numerically at theta*)
    f3 = []
    ths = np.exp(np.linspace(math.log(0.003), math.log(5.0), 24))
    for Pi in (1.0, 3.0, 6.0, 10.0):
        for u in (3.0, 10.0):
            Gs, G2s = [], []
            for th in ths:
                Jc = cost_constant(Pi, th)
                Gs.append(1 - cost_three_arc(Pi, th, u)[0] / Jc)
                J2 = cost_ramp_plateau(Pi, th, u)[0]
                G2s.append(1 - J2 / Jc if not math.isnan(J2) else float("nan"))
            Gs = np.array(Gs)
            k = int(np.nanargmax(Gs))
            Gn = numeric_gain(Pi, float(ths[k]), u)[0]
            f3.append(dict(Pi=Pi, u=u, theta=ths.tolist(), G=Gs.tolist(), G_no_burst=G2s, theta_star=float(ths[k]), G_max=float(Gs[k]),
                           G_numeric_at_star=Gn))
            print(f"F3 Pi {Pi:4.1f} u {u:4.0f}: theta* {ths[k]:.3f}  G_max {100*Gs[k]:.2f}% (numeric {100*Gn:.2f}%)  "
                  f"no-burst G at theta* {100*G2s[k]:.2f}%", flush=True)
    out["tests"]["F3"] = dict(rows=f3, passed=all(ths[0] * 1.01 < r["theta_star"] < ths[-1] * 0.99 for r in f3))
    # ---- F4 Otaniemi
    Pi_o = 3.06
    E8 = [(0.0086, 0.006), (0.026, 0.020), (0.086, 0.055)]
    slope_E8 = float(np.mean([g / t for t, g in E8]))
    f4 = dict(Pi=Pi_o, slope_E8=slope_E8, kappa={}, kappa2={})
    for u in (2.0, 3.0, 4.6, 10.0, math.inf):
        f4["kappa"][str(u)] = kappa(Pi_o, u)
        f4["kappa2"][str(u)] = kappa2(Pi_o, u)
    print(f"F4 Otaniemi: E8 slope {slope_E8:.3f}; kappa (no burst) " + ", ".join(f"u {k}: {v:.3f}" for k, v in f4["kappa"].items())
          + "; kappa2 " + ", ".join(f"u {k}: {v:.3f}" for k, v in f4["kappa2"].items()), flush=True)
    out["tests"]["F4"] = f4
    json.dump(out, open("h1f_results.json", "w"), indent=1)
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(10.5, 4))
        for r in f3:
            if r["u"] == 10.0:
                l, = ax[0].semilogx(r["theta"], 100 * np.array(r["G"]), label=f"Pi = {r['Pi']:g}")
                ax[0].semilogx(r["theta"], 100 * np.array(r["G_no_burst"]), "--", color=l.get_color(), lw=0.8)
                tt = ths[ths < 0.05]
                ax[0].semilogx(tt, 100 * tt * kappa2(r["Pi"], 10.0), "k:", lw=0.8)
        ax[0].set_xlabel("theta = T_k / T_inj"); ax[0].set_ylabel("gain vs constant rate (%)"); ax[0].legend(fontsize=8)
        ax[0].set_title("G(Pi, theta), u_max = 10\n(dashed: no shut-in burst; dotted: theta*kappa2)", fontsize=9)
        G, x, P, r = numeric_gain(3.0, 0.05, 10.0)
        uu, s0 = P.forward(x)
        ax[1].plot(P.tin[1:], x, label="optimal rate u", lw=1); ax[1].plot(P.tin[1:], s0[1:P.nin + 1] / P.sc, label="stress s/s_c")
        uc, sc0 = P.forward(np.ones(P.nin)); ax[1].plot(P.tin[1:], sc0[1:P.nin + 1] / P.sc, "--", label="stress, constant rate")
        ax[1].set_ylim(0, 3); ax[1].set_xlabel("t / T_inj"); ax[1].legend(fontsize=8); ax[1].set_title("Pi = 3, theta = 0.05: bang-singular-bang", fontsize=9)
        fig.tight_layout(); fig.savefig("h1f_fig.png", dpi=120)
    except Exception as e:
        print("figure skipped:", e)
    print("F1 (as registered)", "PASS" if out["tests"]["F1"]["passed_as_registered"] else "FAIL")
    print("F2 pre-registered law", "PASS" if out["tests"]["F2"]["passed_prereg"] else "FAIL",
          "| revised law on fresh Pi", "PASS" if out["tests"]["F2"]["passed_revised_fresh"] else "FAIL")
    print("F3", "PASS" if out["tests"]["F3"]["passed"] else "FAIL")


if __name__ == "__main__":
    main()
