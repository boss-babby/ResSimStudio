"""Writes the figure specs of P1 (figs/*.spec) from the experiment outputs and renders them with ResSim Studio's figure
engine (tools/paper_figs.cpp -> app/src/figure.cpp). Usage: python3 mkfigs.py /path/to/paper_figs
Inputs: fig_synth.json (fig_data.py), ../../exp/h1e_curves.json, ../../exp/h6_results.json; fig 5 comes from
`paper_figs faultrun examples/FAULT_INJECTION.DATA figs/fig5_coupled_fault` (ResSim engine run)."""
import json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.join(HERE, "..", "..", "exp")
FIG = os.path.join(HERE, "figs"); os.makedirs(FIG, exist_ok=True)
EXE = sys.argv[1] if len(sys.argv) > 1 else "/tmp/sd/pfig"
COL = dict(constant="#000000", step="#D55E00", taper="#009E73", opt="#0072B2", grey="#8C8C8C", pink="#CC79A7", sky="#56B4E9", yellow="#E69F00")


class Spec:
    def __init__(s, rows, cols, w, h, legend="auto"): s.L = [f"FIGURE {rows} {cols} {w} {h} {legend}"]
    def panel(s, title, xl, yl, logx=False, logy=False, limits=None, vlines=()):
        s.L.append(f"PANEL {title}|{xl}|{yl}")
        if logx: s.L.append("LOGX")
        if logy: s.L.append("LOGY")
        if limits: s.L.append("LIMITS %g %g %g %g" % tuple(limits))
        for v in vlines: s.L.append(f"VLINE {v:g}")
    def series(s, label, x, y, color, width=1.2, kind="line", dash=0):
        s.L.append(f"SERIES {label}|{color}|{width}|{kind}|{dash}")
        s.L += [f"{a:.9g} {b:.9g}" for a, b in zip(x, y)]
    def text(s, x, y, t): s.L.append(f"TEXT {x:g} {y:g}|{t}")
    def render(s, name):
        p = os.path.join(FIG, name + ".spec"); open(p, "w").write("\n".join(s.L) + "\n")
        subprocess.run([EXE, "render", p, os.path.join(FIG, name)], check=True)
        os.remove(p)


# ---------------------------------------------------------------- Fig 1: counting-window artefact (synthetic, T1)
d = json.load(open(os.path.join(HERE, "fig_synth.json")))
Ti, Th = d["Ti_d"], d["Th_d"]
names = [("constant", "constant", COL["constant"], 0), ("step-up(4)", "step-up (4 steps)", COL["step"], 1),
         ("taper-last30%", "taper (last 30 %)", COL["taper"], 3), ("optimal", "global optimum (Thm 1)", COL["opt"], 0)]
f = Spec(2, 2, 170, 120)
f.panel("(a) Injection schedules, equal volume", "time (days)", "rate (m3/day)", limits=(0, Ti, 0, None) if False else None)
for k, lab, c, dash in names:
    tc = d["t_ctrl_d"]; q = d["sched"][k]["q"]
    xs, ys = [], []
    dt = tc[1] - tc[0]
    for t, v in zip(tc, q): xs += [t - dt / 2, t + dt / 2]; ys += [v, v]
    f.series(lab, xs, ys, c, 1.6 if k == "optimal" else 1.2, "line", dash)
f.panel("(b) Coulomb stress change on the patch", "time (days)", "dS (MPa)", vlines=[Ti])
for k, lab, c, dash in names: f.series("#" + lab, d["sched"][k]["t"], d["sched"][k]["dS"], c, 1.6 if k == "optimal" else 1.2, "line", dash)
f.panel("(c) Cumulative excess events", "time (days)", "N/r - t (days of background)", vlines=[Ti])
for k, lab, c, dash in names: f.series("#" + lab, d["sched"][k]["t"], d["sched"][k]["E"], c, 1.6 if k == "optimal" else 1.2, "line", dash)
f.panel("(d) Ranking depends on the counting window", "schedule", "events / constant-rate events")
r = d["rank"]; c0 = r["constant"]
order = ["constant", "step-up(4)", "taper-last30%", "optimal"]
f.series("counted to shut-in", [i + 0.8 for i in range(4)], [r[k][0] / c0[0] for k in order], COL["sky"], 1, "bars:0.36")
f.series("counted incl. post-injection tail", [i + 1.2 for i in range(4)], [r[k][1] / c0[1] for k in order], COL["opt"], 1, "bars:0.36")
for i, lab in enumerate(["const", "step-up", "taper", "optimum"]): f.text(i + 0.72, 0.08, lab)
f.render("fig1_window_artefact")

# ---------------------------------------------------------------- Fig 2: convexity check + theta law
f = Spec(1, 2, 170, 62)
f.panel("(a) E along random segments q(l) = (1-l) qA + l qB", "l", "excess events (days)")
for i, s in enumerate(d["segments"]):
    c = [COL["opt"], COL["step"], COL["taper"]][i]
    f.series(f"segment {i + 1}" if i == 0 else f"#segment {i + 1}", s["lam"], s["E"], c, 1.4)
    f.series("chord" if i == 0 else "#chord", s["lam"], s["chord"], c, 0.9, "line", 1)
th = [0.0009, 0.0086, 0.026, 0.086, 0.26, 0.86, 2.57]
f.panel("(b) Otaniemi: value of design vs theta = T_k / T_inj", "theta", "change vs constant rate (%)", logx=True)
f.series("global optimum", th, [0.0, -0.6, -2.0, -5.5, -9.2, -7.0, -3.9], COL["opt"], 1.6, "markers")
f.series("operator schedule", th, [91, 51, 28, 8, 0.2, -0.3, -0.3], COL["step"], 1.4, "markers", 1)
f.render("fig2_convexity_theta")

# ---------------------------------------------------------------- Fig 3: Otaniemi 2018 field replication
cpath = os.path.join(EXP, "h1e_curves.json")
if os.path.exists(cpath):
    c = json.load(open(cpath)); tsi, tho = c["t_shutin_d"], c["t_holdout_d"]
    f = Spec(2, 2, 170, 120)
    f.panel("(a) OTN-3 injection rate", "days since 2018-06-04", "rate (l/min)", vlines=[tsi])
    f.series("actual", c["t10_d"], c["q_lpm"], COL["grey"], 0.6)
    xs, ys = [], []
    tc = c["tc_d"]; dt = tc[1] - tc[0] if len(tc) > 1 else 1
    for t, v in zip(tc, c["qopt_lpm"]): xs += [t, t + dt]; ys += [v, v]
    f.series("optimal, same volume", xs, ys, COL["opt"], 1.4)
    f.panel("(b) Wellhead pressure: linear kernel", "days since 2018-06-04", "WHP (bar)", vlines=[tho, tsi])
    f.series("observed", c["whp_t_d"], c["whp"], COL["grey"], 0.5, "scatter")
    f.series("kernel fit (all)", c["t10_d"][:len(c["whp_fit"])], c["whp_fit"], COL["opt"], 1.2)
    f.series("fit on first 60 %", c["t10_d"][:len(c["whp_fit_tr"])], c["whp_fit_tr"], COL["step"], 1.0, "line", 1)
    f.panel("(c) Hourly event rate (M_L >= -0.95)", "days since 2018-06-04", "events / h", logy=True, vlines=[tsi])
    obs = [(t, v) for t, v in zip(c["th_d"], c["counts"]) if v > 0]
    f.series("observed", [a for a, _ in obs], [b for _, b in obs], COL["grey"], 0.5, "scatter")
    f.series("Dieterich fit", c["th_d"], [max(v, 1e-3) for v in c["lam"]], COL["opt"], 1.2)
    f.series("hold-out forecast", c["th_d"], [max(v, 1e-3) for v in c["lam_ho"]], COL["taper"], 1.0, "line", 1)
    f.panel("(d) Cumulative events", "days since 2018-06-04", "events", vlines=[tho, tsi])
    f.series("observed", c["th_d"], c["cum_obs"], COL["constant"], 1.4)
    f.series("Dieterich fit", c["th_d"], c["cum_mod"], COL["opt"], 1.2, "line", 1)
    f.series("hold-out forecast", c["th_d"], c["cum_ho"], COL["taper"], 1.0, "line", 2)
    f.render("fig3_otaniemi")
else:
    print("h1e_curves.json missing - run research/exp/h1e_otaniemi.py first")

# ---------------------------------------------------------------- Fig 4 (v0.5): Otaniemi 2020 (OTN-2) + theta-law at two sites
c2p = os.path.join(EXP, "h8_curves.json")
if os.path.exists(c2p):
    c = json.load(open(c2p)); tsi = c["t_shutin_d"]
    e18 = json.load(open(os.path.join(EXP, "h1e_results.json")))["E8_exploratory"]
    e20 = json.load(open(os.path.join(EXP, "h8_results.json")))["E8_exploratory"]
    f = Spec(2, 2, 170, 120)
    f.panel("(a) OTN-2 (2020) injection rate", "days since start", "rate (l/min)", vlines=[tsi])
    f.series("actual", c["t10_d"], c["q_lpm"], COL["grey"], 0.6)
    xs, ys = [], []
    tc = c["tc_d"]; dt = tc[1] - tc[0] if len(tc) > 1 else 1
    for t, v in zip(tc, c["qopt_lpm"]): xs += [t, t + dt]; ys += [v, v]
    f.series("optimal, same volume", xs, ys, COL["opt"], 1.4)
    f.panel("(b) OTN-2 cumulative events (M >= -1.45)", "days since start", "events", vlines=[tsi])
    f.series("observed", c["th_d"], c["cum_obs"], COL["constant"], 1.4)
    f.series("Dieterich fit", c["th_d"], c["cum_mod"], COL["opt"], 1.2, "line", 1)
    def keyed(e):
        ks = sorted(e, key=lambda k: e[k]["Tk_over_Tinj"]); return ks
    f.panel("(c) Gain of the optimum over constant rate", "theta = T_k / T_inj", "gain (%)", logx=True)
    k18, k20 = keyed(e18), keyed(e20)
    f.series("Otaniemi 2018", [e18[k]["Tk_over_Tinj"] for k in k18], [100 * e18[k]["gain_opt_vs_const"] for k in k18], COL["step"], 1.4)
    f.series("Otaniemi 2020", [e20[k]["Tk_over_Tinj"] for k in k20], [100 * e20[k]["gain_opt_vs_const"] for k in k20], COL["opt"], 1.4)
    f.panel("(d) Cost of the operator schedule vs constant rate", "theta = T_k / T_inj", "extra events (%)", logx=True)
    f.series("Otaniemi 2018", [e18[k]["Tk_over_Tinj"] for k in k18], [100 * e18[k]["actual_vs_const"] for k in k18], COL["step"], 1.4)
    f.series("Otaniemi 2020", [e20[k]["Tk_over_Tinj"] for k in k20], [100 * e20[k]["actual_vs_const"] for k in k20], COL["opt"], 1.4)
    f.render("fig4_otaniemi2020_theta")
else:
    print("h8_curves.json missing - run research/exp/h8 pipeline first")

# ---------------------------------------------------------------- Fig 5: Aspo HRL scope test
h = json.load(open(os.path.join(EXP, "h6_results.json")))
f = Spec(1, 2, 170, 64)
f.panel("(a) Aspo HRL: AE per injection stage vs volume", "injected volume per stage (l)", "AE events (co + post)", logx=True)
tests = sorted(set(s["test"] for s in h["stages"]))
cols = [COL["opt"], COL["step"], COL["taper"], COL["pink"], COL["yellow"], COL["sky"]]
for i, t in enumerate(tests):
    S = [s for s in h["stages"] if s["test"] == t and s["V"] > 0]
    f.series(t, [s["V"] for s in S], [s["co"] + s["post"] for s in S], cols[i % len(cols)], 1.0, "scatter")
M = h["H6a"]["models"]
lab = [("K7_V_cumV", "V (1+Vcum)^g"), ("K6_V_qbeta", "V q^b"), ("K8_V_q_cumV", "V q^b (1+Vcum)^g"), ("K2", "energy"), ("K3", "volume"),
       ("K1", "Dieterich (pressure)"), ("K5_duration", "duration"), ("K0_const", "constant")]
f.panel("(b) Stage-count models: AIC relative to best", "model", "dAIC")
best = min(M[k]["aic"] for k, _ in lab)
f.series("dAIC", [i + 1 for i in range(len(lab))], [M[k]["aic"] - best for k, _ in lab], COL["opt"], 1, "bars:0.6")

f.render("fig5_aspo")
print("AIC legend:", [(i + 1, l, round(M[k]["aic"] - best, 2)) for i, (k, l) in enumerate(lab)])
