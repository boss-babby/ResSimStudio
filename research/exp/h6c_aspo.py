#!/usr/bin/env python3
"""H6c (CYCLE12_PLAN.md, pre-registered): Cycle 5 Aspo pipeline (h6_aspo.py, unchanged) on corrected event times and the
continuous Niemz et al. (2020) catalogue.  usage: python3 h6c_aspo.py  -> h6c.log (stdout), h6c_results.json
Catalogues: i = 196-event triggered catalogue with corrected (event-name) times; ii = continuous, all; iii = continuous MAE >= 2.5;
exploratory: MAE >= 2.15, >= 3.0.  The GFZ 'origin_time' seconds field equals the day of month (data defect) -> event_name is used."""
import contextlib, datetime as dt, io, json, math, os, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h6_aspo as h6

HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "..", "data", "aspo")
ts = lambda s: dt.datetime.strptime(s, "%Y%m%d%H%M%S.%f").replace(tzinfo=dt.timezone.utc).timestamp()


def cat_triggered_fixed():
    out, bad = [], 0
    for line in list(open(os.path.join(DATA, "aspo_catalog_kwiatek2018niemz2020_match.csv")))[1:]:
        s = line.rstrip("\n").split(";")
        try: out.append((ts(s[5].strip()), float(s[4])))
        except ValueError: bad += 1
    return np.array(sorted(out)), bad


def cat_continuous(mmin=None):
    out = []
    for line in list(open(os.path.join(DATA, "aspo_catalog_niemz2020.csv")))[1:]:
        s = line.rstrip("\n").split(";")
        m = float(s[5]) if s[5] not in ("None", "") else float("nan")
        if mmin is None or (m == m and m >= mmin): out.append((ts(s[0]), m))
    return np.array(sorted(out))


def run_pipeline(ev):
    """h6.main() with the event list replaced; its outputs go to a temp dir. Returns (results dict, stdout)."""
    h6.events = lambda: ev
    tmp = tempfile.mkdtemp(); old = h6.HERE; h6.HERE = tmp; buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf): h6.main()
        res = json.load(open(os.path.join(tmp, "h6_results.json")))
    finally:
        h6.HERE = old
        for f in os.listdir(tmp): os.remove(os.path.join(tmp, f))
        os.rmdir(tmp)
    return res, buf.getvalue()


def robustness(ev, K1_Asig, kernel_h6b, rng):
    """stage bootstrap of dAIC(K1 - min(K2, K3)) and dAIC(K2 - K3) on HF1+HF2; leave-one-test-out HF3 ratio."""
    h6.events = lambda: ev
    _, S, raw = h6.build(); train = ["HF1", "HF2"]
    sh = {K: h6.kernel_shapes(S, raw, K)[0] for K in ("K2", "K3")}; sh["K1"] = h6.kernel_shapes(S, raw, "K1", K1_Asig)[0]
    idx = {t: [i for i, s in enumerate(S) if s["test"] == t] for t in train}
    d1, d23 = [], []
    for _ in range(1000):
        pick = [i for t in train for i in rng.choice(idx[t], len(idx[t]), replace=True)]
        Sb = [S[i] for i in pick]; a = {}
        for K in ("K1", "K2", "K3"):
            ll, _ = h6.poisson_fit(Sb, sh[K][pick], train); a[K] = 2 * (len(train) + (K == "K1")) - 2 * ll
        d1.append(a["K1"] - min(a["K2"], a["K3"])); d23.append(a["K2"] - a["K3"])
    # leave-one-test-out on the conventional tests for the HF3 ratio (kernel used by H6b)
    shK = sh[kernel_h6b] if kernel_h6b in sh else h6.kernel_shapes(S, raw, kernel_h6b)[0]
    fac = h6.poisson_fit(S, shK, h6.TESTS)[1]; conv = ["HF1", "HF2", "HF4", "HF6"]
    loto = {}
    for drop in conv:
        ref = [c for c in conv if c != drop and c in ("HF1", "HF2", "HF4", "HF6")]
        ag = [c for c in ("HF1", "HF2") if c != drop] or ["HF4", "HF6"]
        mu = np.mean(np.log([max(fac[c], 1e-12) for c in ag])); loto[drop] = math.exp(math.log(max(fac["HF3"], 1e-12)) - mu)
    pct = lambda v: [float(np.percentile(v, 2.5)), float(np.median(v)), float(np.percentile(v, 97.5))]
    return {"dAIC_K1_minus_best_K2K3": pct(d1), "frac_boot_K1_worse_by_10": float(np.mean(np.array(d1) >= 10)),
            "dAIC_K2_minus_K3": pct(d23), "loto_HF3_ratio": loto, "factors": fac}


if __name__ == "__main__":
    rng = np.random.default_rng(12)
    trig, bad = cat_triggered_fixed()
    cats = {"i_triggered_fixed": trig, "ii_continuous_all": cat_continuous(), "iii_continuous_M2.5": cat_continuous(2.5),
            "x_continuous_M2.15": cat_continuous(2.15), "x_continuous_M3.0": cat_continuous(3.0)}
    print("dropped unparseable names (i):", bad)
    allres = {"dropped_i": bad}
    for name, ev in cats.items():
        res, txt = run_pipeline(ev)
        M = res["H6a"]["models"]; aic = {k: v["aic"] for k, v in M.items()}
        best23 = min(aic["K2"], aic["K3"])
        pa = aic["K1"] - best23 >= 10
        pb = all(min(aic["K0_const"], aic["K5_duration"]) - aic[k] >= 10 for k in ("K2", "K3"))
        b = res["H6b"]; h = res["H1h"]
        pc = (b["ratio_HF3"] < 1) and not b["significant"]
        pdd = h["observed"] < h["dieterich"]
        rob = robustness(ev, M["K1"]["Asig"], b["kernel"], rng)
        r = {"n_events": int(len(ev)), "assigned": res["assigned"], "aic": aic, "K1_Asig": M["K1"]["Asig"],
             "K7_gamma": M["K7_V_cumV"]["beta"], "K7_ci": M["K7_V_cumV"]["beta_ci95"], "K6_beta_ci": M["K6_V_qbeta"]["beta_ci95"],
             "K8_beta_ci": M["K8_V_q_cumV"]["beta_ci95"],
             "H6b": {k: b[k] for k in ("kernel", "ratio_HF3", "t", "lower5_ratio", "significant", "sigma")},
             "H1h": h, "P_a": pa, "P_b": pb, "P_c": pc, "P_d": pdd, "robust": rob}
        allres[name] = r
        print("\n== %s  n=%d  assigned %s" % (name, len(ev), res["assigned"]))
        print("  AIC " + "  ".join("%s %.1f" % (k, v) for k, v in sorted(aic.items(), key=lambda x: x[1])))
        print("  K1 Asig %.3f MPa | K7 gamma %.2f CI %s | K6 beta CI %s | K8 beta CI %s" % (M["K1"]["Asig"], r["K7_gamma"], r["K7_ci"], r["K6_beta_ci"], r["K8_beta_ci"]))
        print("  H6b kernel %s HF3/AG %.3f t %.2f 5%%-bound %.3f significant %s sigma %.2f" % (b["kernel"], b["ratio_HF3"], b["t"], b["lower5_ratio"], b["significant"], b["sigma"]))
        print("  H1h post fraction obs %.3f (n %d) vs Dieterich %.3f" % (h["observed"], h["n"], h["dieterich"]))
        print("  boot dAIC(K1-best23) %s  P(>=10) %.3f | dAIC(K2-K3) %s | LOTO HF3 ratio %s" % (
            ["%.1f" % x for x in rob["dAIC_K1_minus_best_K2K3"]], rob["frac_boot_K1_worse_by_10"], ["%.1f" % x for x in rob["dAIC_K2_minus_K3"]],
            {k: round(v, 3) for k, v in rob["loto_HF3_ratio"].items()}))
        if name.startswith(("ii", "iii")):
            print("  PREDICTIONS P-a %s P-b %s" % (pa, pb) + (" P-c %s P-d %s" % (pc, pdd) if name.startswith("iii") else ""))
    json.dump(allres, open(os.path.join(HERE, "h6c_results.json"), "w"), indent=1, default=float)
