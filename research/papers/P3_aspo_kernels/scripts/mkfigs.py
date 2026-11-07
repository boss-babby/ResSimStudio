#!/usr/bin/env python3
"""P3 figures from the GFZ Aspo data (research/data/aspo) and exp/h6c_results.json.
fig1 hydraulics + AE rate per test; fig2 dAIC per kernel and catalogue; fig3 observed vs expected stage counts;
fig4 per-test rate factors; fig5 post-shut-in fraction.  Output ../figs/*.pdf|png|csv"""
import csv, json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); EXP = os.path.normpath(os.path.join(HERE, "..", "..", "..", "exp"))
FIG = os.path.join(HERE, "..", "figs"); os.makedirs(FIG, exist_ok=True); sys.path.insert(0, EXP)
import h6_aspo as h6, h6c_aspo as h6c
plt.rcParams.update({"font.size": 8.5, "font.family": "DejaVu Sans", "axes.linewidth": 0.8, "legend.frameon": False})
R = json.load(open(os.path.join(EXP, "h6c_results.json")))
def save(fig, n): fig.savefig(os.path.join(FIG, n + ".pdf")); fig.savefig(os.path.join(FIG, n + ".png"), dpi=200); plt.close(fig)
SCH = {"HF1": "conventional, ÄG", "HF2": "conventional, ÄG", "HF3": "progressive cyclic, ÄG", "HF4": "conventional, fgDG",
       "HF5": "pulse hammer, fgDG", "HF6": "conventional, fgG"}

# ---- fig1: hydraulics + AE (continuous, MAE >= 2.5) per test (first record part)
ev = h6c.cat_continuous(2.5); evall = h6c.cat_continuous()
fig = plt.figure(figsize=(7.2, 8.2)); gs = fig.add_gridspec(9, 2, height_ratios=[2, 1, 0.55] * 3, hspace=0.12, wspace=0.42); rows = []
for j, tn in enumerate(h6.TESTS):
    r0, c0 = 3 * (j // 2), j % 2
    ax = fig.add_subplot(gs[r0, c0]); axe = fig.add_subplot(gs[r0 + 1, c0], sharex=ax)
    d = h6.load_test(tn)[0]; t, q, p = d[:, 0], d[:, 1], d[:, 2]; t0 = t[0]; tm = (t - t0) / 60
    ax.plot(tm, np.maximum(q, 0), color="#1f77b4", lw=0.7); ax.set_ylabel("q (L/min)", color="#1f77b4", fontsize=7)
    a2 = ax.twinx(); a2.plot(tm, p, color="0.4", lw=0.6); a2.set_ylabel("p (MPa)", color="0.4", fontsize=7)
    plt.setp(ax.get_xticklabels(), visible=False)
    sel = (ev[:, 0] >= t[0]) & (ev[:, 0] <= t[-1]); sela = (evall[:, 0] >= t[0]) & (evall[:, 0] <= t[-1])
    bw = 1.0 if tm[-1] > 20 else 0.25; bins = np.arange(0, tm[-1] + bw, bw)
    h, _ = np.histogram((ev[sel, 0] - t0) / 60, bins); ha, _ = np.histogram((evall[sela, 0] - t0) / 60, bins)
    axe.bar(bins[:-1], ha / bw, width=bw, color="0.75", align="edge", label="all AE")
    axe.bar(bins[:-1], h / bw, width=bw, color="#d62728", align="edge", label="M$_{AE}$≥2.5")
    axe.set_yscale("symlog", linthresh=1); axe.set_ylabel("AE/min", fontsize=7); axe.set_xlabel("minutes", fontsize=7)
    note = "  (not covered by the catalogue)" if sela.sum() == 0 else ""
    ax.set_title("%s  %s  %d / %d AE%s" % (tn, SCH[tn], sel.sum(), sela.sum(), note), fontsize=7.3, loc="left")
    if j == 1: axe.legend(fontsize=6, loc="upper left")
    if sela.sum() == 0:
        axe.set_yscale("linear"); axe.set_ylim(0, 1); axe.set_yticks([]); axe.text(0.5, 0.5, "no catalogue events on this day", transform=axe.transAxes, ha="center", va="center", fontsize=7)
    rows.append([tn, SCH[tn], int(sel.sum()), int(sela.sum())])
fig.subplots_adjust(left=0.08, right=0.93, top=0.975, bottom=0.03); save(fig, "fig1_tests")

# ---- fig2: dAIC relative to the best pre-registered kernel, per catalogue
cats = [("i_triggered_fixed", "triggered (196), times fixed"), ("x_continuous_M3.0", "continuous M$_{AE}$≥3.0"),
        ("iii_continuous_M2.5", "continuous M$_{AE}$≥2.5 (primary)"), ("x_continuous_M2.15", "continuous M$_{AE}$≥2.15"), ("ii_continuous_all", "continuous, all")]
K = [("K2", "energy Σqp"), ("K3", "volume"), ("K1", "Dieterich(p)"), ("K5_duration", "duration"), ("K0_const", "constant")]
fig, ax = plt.subplots(figsize=(7.2, 2.7)); w = 0.16; out = []
for j, (c, lab) in enumerate(cats):
    a = R[c]["aic"]; b = min(a[k] for k, _ in K)
    vals = [a[k] - b for k, _ in K]; out.append([c] + vals)
    ax.bar(np.arange(len(K)) + (j - 2) * w, np.maximum(vals, 0.3), w, label=lab + "  (n = %d)" % R[c]["n_events"])
ax.set_yscale("log"); ax.set_xticks(range(len(K))); ax.set_xticklabels([l for _, l in K]); ax.axhline(10, color="k", lw=0.6, ls="--")
ax.text(0, 0.5, "best kernel\nin every\ncatalogue", ha="center", fontsize=7)
ax.set_ylabel("ΔAIC vs best kernel"); ax.legend(fontsize=6.3, ncol=2, loc="upper left"); ax.set_ylim(0.3, 3e4)
fig.tight_layout(); save(fig, "fig2_dAIC")
with open(os.path.join(FIG, "fig2_dAIC.csv"), "w", newline="") as f:
    wr = csv.writer(f); wr.writerow(["catalogue"] + [k for k, _ in K]); wr.writerows(out)

# ---- fig3: observed vs expected per stage (primary catalogue), K2 vs K1, all tests (per-test factors)
h6.events = lambda: ev; _, S, raw = h6.build()
n = np.array([s["co"] + s["post"] for s in S]); fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True); rows3 = []
for ax, (KK, lab, kw) in zip(axs, (("K2", "energy kernel", {}), ("K1", "Dieterich pressure kernel", {"Asig": R["iii_continuous_M2.5"]["K1_Asig"]}))):
    sh = h6.kernel_shapes(S, raw, KK, kw.get("Asig"))[0]; fac = h6.poisson_fit(S, sh, h6.TESTS)[1]
    mu = np.array([fac[s["test"]] * sh[i] for i, s in enumerate(S)])
    for tn, c in zip(h6.TESTS, plt.cm.tab10.colors):
        m = np.array([s["test"] == tn for s in S]); ax.loglog(np.maximum(mu[m], 0.3), np.maximum(n[m], 0.3), "o", ms=4, color=c, label=tn)
    lim = [0.3, 600]; ax.plot(lim, lim, "k:", lw=0.8)
    x = np.geomspace(0.3, 600, 100); ax.fill_between(x, np.maximum(x - 2 * np.sqrt(x), 0.3), x + 2 * np.sqrt(x), color="0.9", zorder=0)
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_xlabel("expected stage count (%s)" % lab); ax.set_title(lab, fontsize=8, loc="left")
    rows3 += [[KK, s["test"], s["part"], s["k"], int(n[i]), float(mu[i])] for i, s in enumerate(S)]
axs[0].set_ylabel("observed stage count (M$_{AE}$≥2.5)"); axs[0].legend(fontsize=6.5, loc="upper left")
fig.tight_layout(); save(fig, "fig3_obs_exp")
with open(os.path.join(FIG, "fig3_obs_exp.csv"), "w", newline="") as f:
    wr = csv.writer(f); wr.writerow(["kernel", "test", "part", "stage", "observed", "expected"]); wr.writerows(rows3)

# ---- fig4: per-test rate factors (energy kernel), normalised by the HF1/HF2 mean
fig, ax = plt.subplots(figsize=(3.5, 2.8)); tests = ["HF1", "HF2", "HF3", "HF4", "HF6"]; mk = {"i_triggered_fixed": "s", "ii_continuous_all": "^", "iii_continuous_M2.5": "o", "x_continuous_M2.15": "v"}
for j, c in enumerate(mk):
    f = R[c]["robust"]["factors"]; ag = math.exp(np.mean(np.log([f["HF1"], f["HF2"]])))
    ax.plot(np.arange(len(tests)) + (j - 1.5) * 0.12, [max(f[t], 1e-3) / ag for t in tests], mk[c], ms=5, mfc="white" if c != "iii_continuous_M2.5" else "k", color="k", label=c.split("_", 1)[1].replace("_", " "))
ax.axhline(1, color="0.5", lw=0.6); ax.axvspan(1.6, 2.4, color="#fde0dd", zorder=0)
ax.set_yscale("log"); ax.set_xticks(range(len(tests))); ax.set_xticklabels(tests); ax.set_ylabel("rate factor / mean(HF1, HF2)")
ax.legend(fontsize=6, loc="lower left"); ax.text(2, 0.075, "cyclic", ha="center", fontsize=7, color="#b2182b")
fig.tight_layout(); save(fig, "fig4_factors")

# ---- fig5: post-shut-in fraction observed vs Dieterich prediction (HF1+HF2)
fig, ax = plt.subplots(figsize=(3.5, 2.6)); labs = []
for j, (c, lab) in enumerate(cats):
    h = R[c]["H1h"]; o, pr, nn = h["observed"], h["dieterich"], h["n"]; e = 1.96 * math.sqrt(o * (1 - o) / nn)
    ax.errorbar(j, o, yerr=e, fmt="o", color="k", ms=4); ax.plot(j, pr, "x", color="#d62728", ms=7, mew=1.5); labs.append(lab.split(" (")[0].replace("continuous ", "cont. "))
ax.plot([], [], "o", color="k", label="observed (95 % binomial)"); ax.plot([], [], "x", color="#d62728", label="Dieterich prediction")
ax.set_xticks(range(len(cats))); ax.set_xticklabels(labs, rotation=30, ha="right", fontsize=6.5); ax.set_ylabel("post-shut-in fraction"); ax.set_ylim(0, 0.8)
ax.legend(fontsize=6.5, loc="upper left"); fig.tight_layout(); save(fig, "fig5_postshutin")
with open(os.path.join(FIG, "table1_tests.csv"), "w", newline="") as f:
    wr = csv.writer(f); wr.writerow(["test", "scheme", "AE_M2.5_first_record", "AE_all_first_record"]); wr.writerows(rows)
print("ok")
