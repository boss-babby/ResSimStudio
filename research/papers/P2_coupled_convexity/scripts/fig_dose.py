#!/usr/bin/env python3
"""P2 Fig. 4 + Table 2: convexity-violation frequency vs reached transmissibility multiplier m1 (G1 = H10e/H10f, G2 = H10g/E12).
Binomial GLM  logit p = a + b ln m1 (+ c [G2]) ; likelihood-ratio test for the geometry term; parametric bootstrap for m1(p).
Data are copied from the result logs (research/exp/h10e.log, h10f_*.log, h10g.log, e12.log); confirmatory flag per point.
"""
import csv, json, math, os
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(HERE, "..", "figs"); os.makedirs(FIG, exist_ok=True)
# (geometry, fluid, ramp w, m1, violations, M, status, source)
D = [("G1", "CO2", 0.05, 2.92, 9, 40, "confirmatory", "H10e NS"),
     ("G1", "water", 0.05, 4.53, 19, 40, "confirmatory", "H10e WNS"),
     ("G1", "CO2", 0.08, 2.13, 3, 40, "exploratory", "H10f NSR08 s13"),
     ("G1", "water", 0.08, 2.92, 12, 40, "exploratory", "H10f WNSR08 s13"),
     ("G1", "CO2", 0.08, 2.13, 5, 40, "exploratory", "E12a NSR08 s29"),
     ("G1", "water", 0.08, 2.92, 14, 40, "exploratory", "E12a WNSR08 s29"),
     ("G1", "CO2", 0.15, 1.58, 1, 40, "exploratory", "H10f NSR15"),
     ("G1", "water", 0.15, 1.93, 0, 40, "exploratory", "H10f WNSR15"),
     ("G1", "CO2", 0.5, 1.17, 0, 40, "confirmatory", "H10f-a NSR"),
     ("G1", "water", 0.5, 1.26, 0, 40, "confirmatory", "H10f-a WNSR"),
     ("G2", "CO2", 0.27425, 1.31, 0, 40, "confirmatory", "H10g LOW"),
     ("G2", "CO2", 0.03761, 3.60, 21, 40, "confirmatory", "H10g HIGH"),
     ("G2", "water", 0.27425, 1.35, 0, 40, "confirmatory", "H10g LOW"),
     ("G2", "water", 0.04438, 3.61, 16, 40, "confirmatory", "H10g HIGH"),
     ("G2", "CO2", 0.08606, 2.02, 6, 40, "exploratory", "E12b s29"),
     ("G2", "CO2", 0.05689, 2.60, 13, 40, "exploratory", "E12b s29"),
     ("G2", "water", 0.10155, 2.01, 1, 40, "exploratory", "E12b s29"),
     ("G2", "water", 0.06713, 2.61, 11, 40, "exploratory", "E12b s29")]

lm = np.array([math.log(r[3]) for r in D]); k = np.array([r[4] for r in D], float); n = np.array([r[5] for r in D], float)
g2 = np.array([1.0 if r[0] == "G2" else 0.0 for r in D]); wat = np.array([1.0 if r[1] == "water" else 0.0 for r in D])


def nll(th, X):
    eta = X @ th; return float(np.sum(n * np.logaddexp(0, eta) - k * eta))


def fit(X):
    r = minimize(nll, np.zeros(X.shape[1]), args=(X,), method="BFGS"); return r.x, r.fun


X0 = np.column_stack([np.ones_like(lm), lm]); th0, f0 = fit(X0)
th1, f1 = fit(np.column_stack([X0, g2])); th2, f2 = fit(np.column_stack([X0, wat]))
p_geo = chi2.sf(2 * (f0 - f1), 1); p_fluid = chi2.sf(2 * (f0 - f2), 1)
inv = lambda th, p: math.exp((math.log(p / (1 - p)) - th[0]) / th[1])
rng = np.random.default_rng(7); boot = []
pf = 1 / (1 + np.exp(-(X0 @ th0)))
for _ in range(2000):
    kb = rng.binomial(n.astype(int), pf).astype(float)
    r = minimize(lambda th: float(np.sum(n * np.logaddexp(0, X0 @ th) - kb * (X0 @ th))), th0, method="BFGS")
    if r.x[1] > 0:
        boot.append([inv(r.x, p) for p in (0.05, 0.15, 0.5)])
boot = np.array(boot); ci = np.percentile(boot, [2.5, 97.5], axis=0)
res = {"a": th0[0], "b": th0[1], "m1_at_p": {p: inv(th0, p) for p in (0.05, 0.15, 0.5)},
       "m1_ci95": {p: list(ci[:, i]) for i, p in enumerate((0.05, 0.15, 0.5))},
       "LRT_geometry": {"coef": th1[2], "p": p_geo}, "LRT_fluid": {"coef": th2[2], "p": p_fluid},
       }
# saturated-model deviance for goodness of fit
eps = 1e-12; ph = k / n
sat = -np.sum(k * np.log(np.clip(ph, eps, 1)) + (n - k) * np.log(np.clip(1 - ph, eps, 1)))
res["deviance"] = 2 * (f0 - sat); res["dof"] = len(D) - 2; res["gof_p"] = float(chi2.sf(res["deviance"], len(D) - 2))
# overdispersion: quasi-binomial dispersion phi; widen the bootstrap CI on the log scale by sqrt(phi)
phi = max(1.0, res["deviance"] / res["dof"]); res["phi"] = phi
res["m1_ci95_quasi"] = {p: [math.exp(math.log(res["m1_at_p"][p]) - math.sqrt(phi) * (math.log(res["m1_at_p"][p]) - math.log(lo))),
                            math.exp(math.log(res["m1_at_p"][p]) + math.sqrt(phi) * (math.log(hi) - math.log(res["m1_at_p"][p])))]
                        for p, (lo, hi) in res["m1_ci95"].items()}
json.dump(res, open(os.path.join(FIG, "fig4_dose_response.json"), "w"), indent=1, default=float)
with open(os.path.join(FIG, "fig4_dose_response.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["geometry", "fluid", "ramp_w", "m1", "violations", "M", "status", "source"]); w.writerows(D)

plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans", "axes.linewidth": 0.8})
fig, ax = plt.subplots(figsize=(3.5, 2.8))
mm = np.exp(np.linspace(math.log(1.05), math.log(5), 200))
ax.plot(mm, 1 / (1 + np.exp(-(th0[0] + th0[1] * np.log(mm)))), color="0.25", lw=1.2, label="logistic fit (pooled)")
sty = {("G1", "CO2"): ("o", "#1f77b4"), ("G1", "water"): ("s", "#1f77b4"), ("G2", "CO2"): ("o", "#d62728"), ("G2", "water"): ("s", "#d62728")}
for r in D:
    mk, c = sty[(r[0], r[1])]
    ax.plot(r[3], r[4] / r[5], mk, ms=5, mfc=c if r[6] == "confirmatory" else "white", mec=c, mew=1.1)
ax.axhspan(0, 0.05, color="0.92", zorder=0); ax.axhline(0.15, color="0.6", lw=0.6, ls="--")
ax.set_xscale("log"); ax.set_xticks([1.2, 1.5, 2, 3, 4]); ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xlabel(r"reached transmissibility multiplier $m_1$ (u = 1)"); ax.set_ylabel("fraction of pairs violating convexity")
for lab, mk, c in (("G1  x = 150 m", "o", "#1f77b4"), ("G2  x = 220 m", "o", "#d62728")):
    ax.plot([], [], mk, color=c, label=lab)
ax.plot([], [], "o", mfc="white", mec="0.3", label="open: exploratory"); ax.plot([], [], "s", color="0.3", label="square: water")
ax.legend(fontsize=6.5, frameon=False, loc="upper left"); ax.set_ylim(-0.02, 0.6)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig4_dose_response.pdf")); fig.savefig(os.path.join(FIG, "fig4_dose_response.png"), dpi=200)
print(json.dumps(res, indent=1, default=float))
