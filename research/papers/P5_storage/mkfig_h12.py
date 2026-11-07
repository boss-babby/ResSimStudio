#!/usr/bin/env python3
"""P5 candidate figure from exp/h12_results.json + exp/h12c_results.json (H12a pre-registered, H12b pre-registered, E13 exploratory)."""
import json, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
E = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "exp")
a = json.load(open(os.path.join(E, "h12_results.json"))); c = json.load(open(os.path.join(E, "h12c_results.json")))
fig, ax = plt.subplots(1, 3, figsize=(12.5, 3.8))
# (a) H12a deltas
for k, (case, col) in enumerate((("SW", "#1f77b4"), ("SG", "#d62728"))):
    d = sorted(100 * p["delta"] for p in a["H12a"][case]["pairs"])
    ax[0].plot(d, range(1, 41), "o", ms=3.5, color=col, label="%s: %d/40 violations" % ("water" if case == "SW" else "depleted gas", a["H12a"][case]["violations"]))
ax[0].axvline(0, color="k", lw=0.8); ax[0].set_xlabel("midpoint gap δ (%)  (> 0 = non-convex)"); ax[0].set_ylabel("pair rank")
ax[0].set_title("(a) Signed storage schedules are convex", fontsize=10); ax[0].legend(fontsize=8, loc="upper left")
# (b) schedules at Asig = 10
L = c["levels"]["10.0"]; t = [45 * i for i in range(9)]
for key, lab, col, ls in (("uD", "Dieterich optimum", "#2ca02c", "-"), ("uK_lex", "Kaiser optimum (lexicographic)", "#9467bd", "--"), ("uR", "robust (min-max)", "k", ":")):
    u = L[key]; ax[1].step(t, u + [u[-1]], where="post", color=col, ls=ls, lw=1.6, label=lab)
ax[1].step(t, [1, 1, 1, 1, -1, -1, -1, -1, -1], where="post", color="0.6", lw=1, label="constant")
ax[1].axhline(0, color="k", lw=0.5); ax[1].set_xlabel("time (d)"); ax[1].set_ylabel("rate multiplier (+ inject, − withdraw)")
ax[1].set_title("(b) Optimal cycles, Π ≈ 0.4", fontsize=10); ax[1].legend(fontsize=7.5)
# (c) cross-costs
labs, xs = [], []
import numpy as np
w = 0.2
for j, (lv, pi) in enumerate((("1.0", "Π ≈ 4"), ("10.0", "Π ≈ 0.4"))):
    L = c["levels"][lv]
    vals_D = [L["base_ND_rel"], 0.0, L["XD_lex"], L["robust_ND_rel"]]
    vals_K = [L["base_NK_rel"], L["XK"], 0.0, L["robust_NK_rel"]]
    x = np.arange(4) + 5 * j
    ax[2].bar(x - w / 2, [100 * v for v in vals_D], w, color="#2ca02c", label="cost under Dieterich" if j == 0 else None)
    ax[2].bar(x + w / 2, [100 * v for v in vals_K], w, color="#9467bd", label="cost under Kaiser" if j == 0 else None)
    labs += ["const.", "Dieterich\nopt.", "Kaiser\nopt.", "robust"]; xs += list(x)
    ax[2].text(1.5 + 5 * j, 44, pi, ha="center", fontsize=9)
ax[2].set_xticks(xs); ax[2].set_xticklabels(labs, fontsize=7); ax[2].set_ylabel("excess over model optimum (%)"); ax[2].set_ylim(0, 48)
ax[2].set_title("(c) Memory-model regret", fontsize=10); ax[2].legend(fontsize=8, loc="center right")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "figs", "fig_h12_storage." + ext), dpi=200)
