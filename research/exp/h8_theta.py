#!/usr/bin/env python3
"""Cycle 6 (POST HOC, not pre-registered): theta-law forecast of the 2020 OTN-2 gain from the 2018 theta-scan.
Forecast = stress-weighted mean of G_2018(theta_j) over the 2020 kernel components (weights a_j = steady-state
pressure share, theta_j = T_j / T_inj,2020), G_2018 interpolated in log(theta) from h1e_results.json E8."""
import json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
r18 = json.load(open(os.path.join(HERE, "h1e_results.json"))); r20 = json.load(open(os.path.join(HERE, "h8_results.json")))
e8 = r18["E8_exploratory"]
th18 = np.array([v["Tk_over_Tinj"] for v in e8.values()]); g18 = np.array([v["gain_opt_vs_const"] for v in e8.values()])
k = r20["E1"]["kernel_all"]; a = np.array(k["a"]); T = np.array(k["T"]); Tinj = r20["data"]["shut_in_h"]
th = T / Tinj
G = np.where(th < th18[0], 0.0, np.interp(np.log(np.maximum(th, 1e-12)), np.log(th18), g18, right=g18[-1]))
fc = float((a * G).sum() / a.sum())
real = r20["E4"]["gain_vs_constant"]
th_eff = float(math.exp((a[a > 0] * np.log(th[a > 0])).sum() / a[a > 0].sum()))
out = dict(theta_components=th.tolist(), weights=(a / a.sum()).tolist(), G2018_at_theta=G.tolist(), forecast_gain=fc,
           realised_gain=real, ratio=real / fc if fc > 0 else None, theta_eff_logmean=th_eff,
           Pi_2018=r18["E4"]["Pi"], Pi_2020=r20["E4"]["Pi"], Pi_law_2020=r20["E4"]["Pi_law_gain_pct"] / 100,
           note="post hoc; Pi differs (3.06 vs 3.43), no Pi correction applied")
json.dump(out, open(os.path.join(HERE, "h8_theta.json"), "w"), indent=1)
print(f"theta-law forecast {100 * fc:.1f}%  realised {100 * real:.1f}%  Pi-law {r20['E4']['Pi_law_gain_pct']:.1f}%  theta_eff {th_eff:.3f}")
