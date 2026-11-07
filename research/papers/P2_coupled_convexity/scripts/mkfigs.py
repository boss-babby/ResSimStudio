#!/usr/bin/env python3
"""P2 figures 1, 2, 3, 5 from the archived experiment outputs in research/exp (no simulator needed).
Fig. 4 (dose-response) is scripts/fig_dose.py.  Output: ../figs/figN_*.pdf + .png (+ .csv of plotted data)."""
import csv, json, math, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.normpath(os.path.join(HERE, "..", "..", "..", "exp")); FIG = os.path.join(HERE, "..", "figs")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 8.5, "font.family": "DejaVu Sans", "axes.linewidth": 0.8, "legend.frameon": False})
J = lambda f: json.load(open(os.path.join(EXP, f)))


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".pdf")); fig.savefig(os.path.join(FIG, name + ".png"), dpi=200); plt.close(fig)


# ---------------- Fig. 1: model and two-way coupling law ----------------
fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.5), gridspec_kw={"width_ratios": [1.6, 1]})
a.add_patch(Rectangle((0, 0), 400, 50, fc="#f2efe6", ec="k", lw=0.8))
for xf, c, lab, ha, dx in ((150, "#1f77b4", "G1  x = 150 m", "right", -4), (220, "#d62728", "G2  x = 220 m", "left", 4)):
    a.plot([xf, xf], [0, 50], color=c, lw=2.2); a.text(xf + dx, 6, lab, color=c, fontsize=7.5, ha=ha)
a.plot([5, 5], [0, 35], color="k", lw=0.8); a.plot([5, 5], [35, 50], color="k", lw=3.5); a.text(9, 30, "injector, perforated\nlayers 8-10", fontsize=7)
a.add_patch(Rectangle((390, 0), 10, 50, fc="#d9d9d9", ec="none")); a.text(352, 22, "pore-volume\nboundary\n(MULTPV 10$^4$)", fontsize=6.5, ha="center")
a.annotate("", xy=(120, 43), xytext=(12, 43), arrowprops=dict(arrowstyle="->", lw=1.0))
a.text(22, 47.5, "CO$_2$ (CO2STORE) or water", fontsize=7)
a.set_xlim(-5, 405); a.set_ylim(52, -2); a.set_xlabel("x (m)"); a.set_ylabel("depth below 1000 m (m)")
a.set_title("(a) 40 x 1 x 10 section, k = 50 / 5 mD, sealing fault MULTFLT 0.02", fontsize=8, loc="left")
cff = np.linspace(-1, 3, 300); sn = 95.0
for w, ls in ((0.5, "-"), (0.15, "--"), (0.08, "-."), (0.05, ":")):
    x = np.clip(cff / (w * sn), 0, 1); b.plot(cff, 50 ** x, ls, color="k", lw=1, label="w = %g" % w)
b.axhline(2, color="0.6", lw=0.6); b.set_yscale("log"); b.set_xlabel("CFF on face (bar)"); b.set_ylabel("T multiplier  m = P$^{x}$")
b.set_title("(b) two-way law, P = 50, x = CFF/(w$\\sigma'_n$)", fontsize=8, loc="left"); b.legend(fontsize=6.5)
fig.tight_layout(); save(fig, "fig1_model")

# ---------------- Fig. 2: midpoint gap distributions per regime ----------------
groups = [("passive fault (one-way)", [("h10_results.json", k) for k in ("N", "F", "W", "WF")] +
           [("h10c_results.json", k) for k in ("NN", "WNN", "NH", "WNH")] + [("h10e_ctrl_results.json", k) for k in ("NSC", "WNSC")]),
          ("two-way, smooth law", [("h10f_a_results.json", k) for k in ("NSR", "WNSR")] +
           [("h10g_results.json", k) for k in ("X:22:0:32.68:0.27425:50", "X:22:1:27.06:0.27425:50")]),
          ("two-way, steep law", [("h10e_results.json", k) for k in ("NS", "NHS", "WNS", "WNHS")] +
           [("h10g_results.json", k) for k in ("X:22:0:32.68:0.03761:50", "X:22:1:27.06:0.04438:50")])]
lab = {"X:22:0:32.68:0.27425:50": "G2N-LOW", "X:22:1:27.06:0.27425:50": "G2W-LOW",
       "X:22:0:32.68:0.03761:50": "G2N-HIGH", "X:22:1:27.06:0.04438:50": "G2W-HIGH"}
fig, ax = plt.subplots(figsize=(7.2, 2.9)); xs, names, rows = 0, [], []
cols = ["#4d7fb3", "#6aa84f", "#cc4125"]
for gi, (gname, items) in enumerate(groups):
    x0 = xs
    for f, k in items:
        c = J(f)["cases"][k]; d = 100 * np.array([p["delta"] for p in c["pairs"]])
        jit = (np.random.default_rng(int(10 * xs)).random(len(d)) - 0.5) * 0.5
        ax.scatter(xs + jit, d, s=5, color=cols[gi], alpha=0.7, lw=0)
        ax.plot([xs - 0.35, xs + 0.35], [100 * c["eps"]] * 2, color="k", lw=1)
        ax.text(xs, 23, "%d/%d" % (c["violations"], len(d)), ha="center", fontsize=6)
        names.append(lab.get(k, k)); rows.append([gname, lab.get(k, k), c["violations"], len(d), c["eps"], float(np.median(d))]); xs += 1
    ax.text((x0 + xs - 1) / 2, -68, gname, ha="center", fontsize=7.5, color=cols[gi]); xs += 0.6
ax.axhline(0, color="0.5", lw=0.6)
ax.set_xticks([i + 0.6 * (i >= 10) + 0.6 * (i >= 14) for i in range(len(names))]); ax.set_xticklabels(names, rotation=60, fontsize=6.5)
ax.set_ylabel(r"midpoint gap $\delta$ (%)"); ax.set_ylim(-72, 27)
ax.text(0.995, 0.03, "bar = tolerance $\\varepsilon$; violation if $\\delta > \\varepsilon$; numbers = violations / pairs",
        transform=ax.transAxes, ha="right", fontsize=6.5)
fig.tight_layout(); save(fig, "fig2_midpoint_gaps")
with open(os.path.join(FIG, "fig2_midpoint_gaps.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["regime", "case", "violations", "pairs", "eps", "median_delta_pct"]); w.writerows(rows)

# ---------------- Fig. 3: rate-and-state objectives are not convex ----------------
RX = re.compile(r"PAIR\s+\d+ .* N (\S+) (\S+) mid (\S+) \| M0 (\S+) (\S+) mid (\S+) \| D (\S+) (\S+) mid (\S+)")
def pairs(fn):
    out = []
    for line in open(os.path.join(EXP, fn)):
        m = RX.search(line)
        if m:
            v = [float(t) for t in m.groups()]; out.append(v)
    return np.array(out)
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.4)); rows = []
for ax, (fn, title) in zip(axs, (("h7_200.log", "(a) single RSF fault, 200 pairs"), ("h7b.log", "(b) 16-member RSF ensemble, 60 pairs"))) :
    P = pairs(fn)
    for j, (nm, c, mk) in enumerate((("RSF event count", "#1f77b4", "o"), ("RSF moment", "#d62728", "^"), ("Dieterich (same pressure)", "0.3", "s"))):
        avg = (P[:, 3 * j] + P[:, 3 * j + 1]) / 2; r = P[:, 3 * j + 2] / avg
        ax.scatter(np.full(len(r), j) + (np.random.default_rng(j).random(len(r)) - 0.5) * 0.5, r, s=5, color=c, marker=mk, alpha=0.6, lw=0)
        rows.append([fn, nm, int(np.sum(P[:, 3 * j + 2] > avg * (1 + 1e-9))), len(r), float(np.max(r))])
        # moment is printed with 3 significant figures in the logs; annotate the full-precision counts of F20 / F22
        off = {("h7_200.log", 1): 82, ("h7b.log", 1): 19}.get((fn, j), rows[-1][2]); rows[-1].append(off)
        ax.text(j, 1.42, "%d/%d" % (off, len(r)), ha="center", fontsize=6.5)
    ax.axhline(1, color="k", lw=0.7); ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["count", "moment", "Dieterich"])
    ax.set_ylim(0.3, 1.5); ax.set_title(title, fontsize=7.5, loc="left")
axs[0].set_ylabel("mid / average of endpoints")
axs[0].text(0.02, 0.03, "> 1 : convexity violated", transform=axs[0].transAxes, fontsize=6.5)
# panel c: optimisation for the non-convex moment
h11 = [l for l in open(os.path.join(EXP, "h11.log"))]
tr = [float(re.search(r"RSF (\S+)", l).group(1)) for l in h11 if l.startswith("TRAIN") or l.startswith("GP-EI")]
best = np.minimum.accumulate(tr); ax = axs[2]
ax.plot(np.arange(1, len(tr) + 1), np.array(tr) / 1e15, ".", color="0.6", ms=3, label="evaluated")
ax.plot(np.arange(1, len(tr) + 1), best / 1e15, "-", color="#d62728", lw=1.2, label="GP-EI best so far")
ax.axhline(2.35775, color="k", lw=0.8, ls="--", label="constant rate"); ax.axhline(2.21503, color="#1f77b4", lw=0.8, label="ES median best (180 evals)")
ax.axhline(4.0522, color="0.4", lw=0.8, ls=":", label="Dieterich-surrogate optimum")
ax.axvline(12.5, color="0.8", lw=0.6); ax.set_ylim(2.1, 6.2); ax.set_xlabel("RSF evaluation"); ax.set_ylabel("mean moment (10$^{15}$ N m)")
ax.set_title("(c) designing for moment (E = 16)", fontsize=7.5, loc="left"); ax.legend(fontsize=5.5, loc="upper right")
fig.tight_layout(); save(fig, "fig3_rsf_nonconvex")
with open(os.path.join(FIG, "fig3_rsf_nonconvex.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["log", "objective", "strict_violations_from_log", "pairs", "max_ratio", "violations_full_precision"]); w.writerows(rows)

# ---------------- Fig. 5: SCP certified design on the convex regime ----------------
d = J("h10b_results.json"); fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.4))
for ax, k, t in ((axs[0], "N", "(a) CO$_2$, fault in plume (x = 150 m)"), (axs[1], "F", "(b) CO$_2$, pressure-only fault (x = 300 m)")):
    h = d[k]["primary"]["hist"]; runs = [r[3] for r in h]; val = [r[2] for r in h]
    ax.plot(runs, val, "o-", color="#d62728", ms=3, lw=1, label="SCP iterate")
    ax.axhline(d[k]["grid_best"], color="k", ls="--", lw=0.8, label="best of %d-point grid" % d[k]["grid_n"])
    ax.set_xlabel("simulator runs"); ax.set_ylabel("expected events N"); ax.set_title(t, fontsize=7.5, loc="left")
    ax.set_ylim(min(val) * 0.97, min(max(val), 2 * d[k]["grid_best"])); ax.legend(fontsize=6.5)
    ax.text(0.97, 0.55, "final %.4f\n(%.2f %% vs grid)" % (d[k]["primary"]["N"], 100 * d[k]["primary"]["gap"]), transform=ax.transAxes, ha="right", fontsize=7)
fig.tight_layout(); save(fig, "fig5_scp")
print("figs written to", os.path.normpath(FIG))
