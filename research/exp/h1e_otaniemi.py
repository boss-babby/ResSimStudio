#!/usr/bin/env python3
"""H1e: field replication of the convex seismicity-design theory on Otaniemi 2018 (pre-registered in CYCLE3_PLAN.md).

Pipeline
  1. hourly / 10-min series from data/otaniemi (exp/h1e_fetch.py)
  2. E1  rate -> wellhead-pressure kernel  WHP = p0 + c1 q + c2 q^2 + sum_j a_j y_j,  y_j = exp-filter(q, T_j), a_j >= 0 (NNLS)
         reservoir pressure p_res = sum_j a_j y_j (instantaneous friction / skin terms excluded)
  3. E2  Dieterich (1994) rate with stress proxy S = p_res filtered by a diffusion delay tau_d:
         R = r exp(S/P*) / (1 + I/t_a),  I = int exp(S/P*) dt;  Poisson likelihood of hourly counts M_L >= M_c
  4. E3  hold-out (fit before stage 4) + seismogenic-index baseline
  5. E4-E6 convex retrospective optimum (6-h controls, same volume, 0 <= q <= q_max, WHP <= observed max), window artefact,
         robustness; H1d: b-value, P(M_L >= 2), maximum volume under the actual hazard
Writes exp/h1e_results.json and exp/h1e_fig.png.  Runtime ~ 2-4 min on 2 cores.
"""
import csv, json, math, os, sys, time
import numpy as np
from scipy.optimize import nnls, minimize
from scipy.signal import lfilter

HERE = os.path.dirname(os.path.abspath(__file__))
SITE = os.environ.get("SITE", "otn3")             # otn3 = Otaniemi 2018 (Cycle 3), otn2 = Otaniemi 2020 (Cycle 6, H8')
DATA = os.path.join(HERE, "..", "data", "otaniemi" if SITE == "otn3" else "otaniemi2020")
TAG = "h1e" if SITE == "otn3" else "h8"
T0_LABEL = "2018-06-04" if SITE == "otn3" else "2020-05-05"
t_start = time.time()
OUT = {}

# ------------------------------------------------------------------ data
H = np.array([[float(r["t_h"]), float(r["rate_m3min"]), float(r["whp_bar"])] for r in csv.DictReader(open(os.path.join(DATA, "hydraulic_10min.csv")))])
C = np.array([[float(r["t_h"]), float(r["ML"])] for r in csv.DictReader(open(os.path.join(DATA, "catalogue.csv")))])
dt10 = 1.0 / 6.0
T_END = 2682.0 if SITE == "otn3" else float(math.ceil(C[:, 0].max()))   # end of catalogue (h)
n10 = int(math.ceil(T_END / dt10))
q10 = np.zeros(n10); whp10 = np.full(n10, np.nan)
k = np.clip((H[:, 0] / dt10).astype(int), 0, n10 - 1)
q10[k] = np.maximum(H[:, 1], 0.0); whp10[k] = H[:, 2]
T_SI = 1168.0 if SITE == "otn3" else float(math.ceil(H[H[:, 1] > 0.01, 0].max() * 6) / 6)   # shut-in (h after T0 UTC)
q10[int(T_SI / dt10) + 1:] = 0.0
V_TOT = q10.sum() * 10.0                          # m3 (q in m3/min, 10-min bins)
Q_MAX = float(np.percentile(H[:, 1], 99.9))
WHP_MAX = float(np.nanmax(whp10))
STAGE4 = 657.0 if SITE == "otn3" else 0.6 * T_SI  # otn3: start of stage 4 (2018-07-01 06:55 UTC); otn2: pre-registered 60 %
OUT["data"] = dict(volume_m3=V_TOT, q_max_m3min=Q_MAX, whp_max_bar=WHP_MAX, shut_in_h=T_SI, end_h=T_END, events=int(len(C)))

# ------------------------------------------------------------------ magnitude of completeness and b (Aki-Utsu)
ml = C[:, 1]
bins = np.arange(-3.5, 2.0, 0.1)
hcnt, _ = np.histogram(ml, bins)
mc_maxc = bins[np.argmax(hcnt)] + 0.05
MC = round(mc_maxc + 0.2, 2)
sel = ml >= MC
b_val = math.log10(math.e) / (ml[sel].mean() - MC)
b_err = b_val / math.sqrt(sel.sum())
OUT["catalogue"] = dict(Mc_maxcurv=mc_maxc, Mc=MC, n_above=int(sel.sum()), b=b_val, b_se=b_err, ML_max=float(ml.max()))
print(f"Mc {MC:.2f} (max curvature {mc_maxc:.2f} + 0.2), N(M>=Mc) {sel.sum()}, b {b_val:.3f} +- {b_err:.3f}, ML max {ml.max():.2f}")

# hourly counts
NH = int(math.ceil(T_END))
counts = np.bincount(np.clip(C[sel, 0].astype(int), 0, NH - 1), minlength=NH).astype(float)

# ------------------------------------------------------------------ E1 kernel
TJ = np.array([0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000.0])


def expfilt(q, T, dt):
    """y' = (q - y)/T, exact for piecewise-constant q (recursive filter); columns of a 2-D q are filtered independently."""
    a = math.exp(-dt / T)
    return lfilter([1 - a], [1, -a], q, axis=0)


Y10 = np.stack([expfilt(q10, T, dt10) for T in TJ], axis=1)          # n10 x J
ok = ~np.isnan(whp10)
ok &= np.arange(n10) * dt10 < T_SI + 2
A_full = np.column_stack([np.ones(n10), -np.ones(n10), q10, q10 ** 2, Y10])


def fit_kernel(mask):
    coef, _ = nnls(A_full[mask], whp10[mask])
    return coef


def p_res_from(coef, Y):
    return Y @ coef[4:]


tr = ok & (np.arange(n10) * dt10 < 0.6 * T_SI)
te = ok & (np.arange(n10) * dt10 >= 0.6 * T_SI)
c_tr = fit_kernel(tr)
pred = A_full @ c_tr
rng = np.nanmax(whp10[ok]) - np.nanmin(whp10[ok])
rmse_te = float(np.sqrt(np.mean((pred[te] - whp10[te]) ** 2)))
rmse_tr = float(np.sqrt(np.mean((pred[tr] - whp10[tr]) ** 2)))
c_all = fit_kernel(ok)
rmse_all = float(np.sqrt(np.mean(((A_full @ c_all)[ok] - whp10[ok]) ** 2)))
OUT["E1"] = dict(rmse_train_bar=rmse_tr, rmse_test_bar=rmse_te, whp_range_bar=float(rng), rel_test=rmse_te / rng, rmse_all_bar=rmse_all,
                 kernel_all=dict(p0=float(c_all[0] - c_all[1]), c1=float(c_all[2]), c2=float(c_all[3]), a=[float(x) for x in c_all[4:]], T=[float(x) for x in TJ]),
                 verdict="PASS" if rmse_te / rng <= 0.10 else ("WEAK" if rmse_te / rng <= 0.20 else "FAIL (non-stationary)"))
print(f"E1 kernel: RMSE train {rmse_tr:.1f} bar, test {rmse_te:.1f} bar = {100 * rmse_te / rng:.1f}% of range {rng:.0f} bar -> {OUT['E1']['verdict']}")
print("   kernel a_j (bar per m3/min):", np.round(c_all[4:], 1), " c1", round(c_all[2], 1), " c2", round(c_all[3], 1))

# hourly reservoir pressure operator: S_h = M q_ctrl (built from unit pulses through the same filters)
def hourly_from10(x10):
    n = NH * 6
    x = np.zeros(n) if x10.ndim == 1 else np.zeros((n,) + x10.shape[1:])
    x[:min(n, len(x10))] = x10[:n]
    return x.reshape((NH, 6) + x10.shape[1:]).mean(axis=1)


# ------------------------------------------------------------------ E2 Dieterich fit
def lnR_series(s, ln_r, ln_ta):
    """s = S/P* (hourly); returns ln R (per hour) and ln I (cumulative), exact-ish with hourly trapezoid in logs."""
    lnI = np.logaddexp.accumulate(s)              # ln sum_{j<=i} exp(s_j) * 1 h
    lnR = ln_r + s - np.logaddexp(0.0, lnI - ln_ta)
    return lnR, lnI


def delay(p_h, tau):
    if tau <= 1e-6:
        return p_h
    return expfilt(p_h, tau, 1.0)


def nll(theta, p_h, idx):
    ln_r, ln_P, ln_ta, ln_tau = theta
    P = math.exp(ln_P)
    s = delay(p_h, math.exp(ln_tau)) / P
    lnR, _ = lnR_series(s, ln_r, ln_ta)
    lam = np.exp(np.clip(lnR[idx], -50, 50))
    return float(np.sum(lam - counts[idx] * np.log(lam + 1e-300)))


BOUNDS = [(-40, 5), (math.log(0.5), math.log(3000)), (math.log(1), math.log(1e6)), (math.log(0.01), math.log(3000))]


def fit_dieterich(p_h, idx, starts=None):
    best = None
    starts = starts or [(-12, math.log(15), math.log(300), math.log(5)), (-20, math.log(8), math.log(1000), math.log(30)),
                        (-8, math.log(30), math.log(100), math.log(1)), (-15, math.log(12), math.log(3000), math.log(100))]
    for x0 in starts:
        r = minimize(nll, np.array(x0, float), args=(p_h, idx), method="Nelder-Mead", bounds=BOUNDS,
                     options=dict(maxiter=4000, maxfev=6000, xatol=1e-4, fatol=1e-3))
        if best is None or r.fun < best.fun:
            best = r
    return best


p_h_all = hourly_from10(p_res_from(c_all, Y10))
idx_all = np.arange(NH)
fit = fit_dieterich(p_h_all, idx_all)
ln_r, ln_P, ln_ta, ln_tau = fit.x
P_star, t_a, tau_d = math.exp(ln_P), math.exp(ln_ta), math.exp(ln_tau)
s_all = delay(p_h_all, tau_d) / P_star
lnR, lnI = lnR_series(s_all, ln_r, ln_ta)
lam = np.exp(lnR)
cum_obs, cum_mod = np.cumsum(counts), np.cumsum(lam)
corr = float(np.corrcoef(cum_obs, cum_mod)[0, 1])
i_si = int(T_SI)
tail_obs, tail_mod = float(counts[i_si:].sum()), float(lam[i_si:].sum())
OUT["E2"] = dict(r_per_h=math.exp(ln_r), P_star_bar=P_star, t_a_h=t_a, tau_d_h=tau_d, nll=float(fit.fun), corr_cum=corr,
                 N_obs=float(cum_obs[-1]), N_mod=float(cum_mod[-1]), tail_obs=tail_obs, tail_mod=tail_mod,
                 verdict="PASS" if corr >= 0.98 and abs(tail_mod / max(tail_obs, 1) - 1) <= 0.30 else
                 ("WEAK" if corr >= 0.95 and 0.5 <= tail_mod / max(tail_obs, 1) <= 2 else "FAIL"))
print(f"E2 fit: P* {P_star:.2f} bar, t_a {t_a:.0f} h, tau_d {tau_d:.1f} h, r {math.exp(ln_r):.3g}/h; corr(cum) {corr:.4f}; "
      f"N {cum_obs[-1]:.0f} vs {cum_mod[-1]:.0f}; tail {tail_obs:.0f} vs {tail_mod:.0f} -> {OUT['E2']['verdict']}")

# ------------------------------------------------------------------ E3 hold-out
tr10 = ok & (np.arange(n10) * dt10 < STAGE4)
c_ho = fit_kernel(tr10)
p_h_ho = hourly_from10(p_res_from(c_ho, Y10))
idx_tr = np.arange(int(STAGE4))
fit_ho = fit_dieterich(p_h_ho, idx_tr)
s_ho = delay(p_h_ho, math.exp(fit_ho.x[3])) / math.exp(fit_ho.x[1])
lam_ho = np.exp(lnR_series(s_ho, fit_ho.x[0], fit_ho.x[2])[0])
N_si_obs, N_end_obs = float(counts[:i_si].sum()), float(counts.sum())
N_si_ho, N_end_ho = float(lam_ho[:i_si].sum()), float(lam_ho.sum())
Vcum_h = np.cumsum(hourly_from10(q10) * 60.0)
sigma_si = counts[:int(STAGE4)].sum() / Vcum_h[int(STAGE4) - 1]
N_si_SI = float(sigma_si * Vcum_h[i_si - 1])
err = lambda a, b: a / b - 1
OUT["E3"] = dict(N_shutin_obs=N_si_obs, N_shutin_forecast=N_si_ho, N_shutin_SI_baseline=N_si_SI, N_end_obs=N_end_obs, N_end_forecast=N_end_ho,
                 err_shutin=err(N_si_ho, N_si_obs), err_end=err(N_end_ho, N_end_obs), err_SI=err(N_si_SI, N_si_obs),
                 params=dict(P_star_bar=math.exp(fit_ho.x[1]), t_a_h=math.exp(fit_ho.x[2]), tau_d_h=math.exp(fit_ho.x[3]), r_per_h=math.exp(fit_ho.x[0])))
OUT["E3"]["verdict"] = "PASS" if (abs(OUT["E3"]["err_shutin"]) <= 0.3 and abs(OUT["E3"]["err_end"]) <= 0.3 and abs(OUT["E3"]["err_shutin"]) < abs(OUT["E3"]["err_SI"])) else "FAIL"
print(f"E3 hold-out: shut-in {N_si_ho:.0f} vs obs {N_si_obs:.0f} ({100 * OUT['E3']['err_shutin']:+.0f}%), end {N_end_ho:.0f} vs {N_end_obs:.0f} "
      f"({100 * OUT['E3']['err_end']:+.0f}%), SI baseline at shut-in {N_si_SI:.0f} ({100 * OUT['E3']['err_SI']:+.0f}%) -> {OUT['E3']['verdict']}")

# ------------------------------------------------------------------ E9 (Cycle 6, pre-registered): cross-stimulation transfer
if SITE == "otn2":
    ref = json.load(open(os.path.join(HERE, "h1e_results.json")))["E2"]
    P9, ta9, tau9 = ref["P_star_bar"], ref["t_a_h"], ref["tau_d_h"]
    s9 = delay(p_h_all, tau9) / P9
    shape = np.exp(lnR_series(s9, 0.0, math.log(ta9))[0])
    r9 = counts.sum() / shape.sum()                  # Poisson MLE of the productivity alone
    lam9 = r9 * shape
    corr9 = float(np.corrcoef(np.cumsum(counts), np.cumsum(lam9))[0, 1])
    e_si9 = float(lam9[:i_si].sum() / max(counts[:i_si].sum(), 1) - 1)
    OUT["E9"] = dict(P_star_bar_2018=P9, t_a_h_2018=ta9, tau_d_h_2018=tau9, r_per_h=float(r9), corr_cum=corr9, err_shutin=e_si9,
                     tail_obs=float(counts[i_si:].sum()), tail_mod=float(lam9[i_si:].sum()),
                     verdict="PASS" if corr9 >= 0.95 and abs(e_si9) <= 0.30 else ("WEAK" if corr9 >= 0.95 else "FAIL"))
    print(f"E9 transfer (2018 P*, t_a, tau_d; r refit): corr(cum) {corr9:.4f}, shut-in error {100 * e_si9:+.0f}%, tail "
          f"{counts[i_si:].sum():.0f} vs {lam9[i_si:].sum():.0f} -> {OUT['E9']['verdict']}")

# ------------------------------------------------------------------ E4-E6 design
DTC = 6.0                                          # control interval (h)
NC = int(math.ceil(T_SI / DTC))
# unit pulse responses (hourly) of p_res and of the WHP for each control interval, via 10-min filters
U10 = np.zeros((n10, NC))
for j in range(NC):
    U10[int(j * DTC / dt10):int(min((j + 1) * DTC, T_SI) / dt10), j] = 1.0
YU = [expfilt(U10, T, dt10) for T in TJ]          # each n10 x NC


def operators(coef, tau):
    Pres10 = sum(coef[4 + jj] * YU[jj] for jj in range(len(TJ)))       # n10 x NC  (bar per m3/min)
    Ph = hourly_from10(Pres10)
    Sh = delay(Ph, tau) if tau > 1e-6 else Ph
    return Sh, Pres10


def events_and_grad(x, Sh, P, ln_r, ln_ta, nh):
    s = (Sh[:nh] @ x) / P
    m = s.max()
    w = np.exp(s - m)
    I = w.sum() * math.exp(m)
    ta, r = math.exp(ln_ta), math.exp(ln_r)
    N = r * ta * math.log1p(I / ta)
    g = r * ta / (ta + I) * math.exp(m) * (Sh[:nh].T @ w) / P
    return N, g


def actual_ctrl():
    qa = hourly_from10(q10)
    return np.array([qa[int(j * DTC):int(min((j + 1) * DTC, T_SI))].mean() for j in range(NC)])


x_act = actual_ctrl()
dur = np.array([min((j + 1) * DTC, T_SI) - j * DTC for j in range(NC)]) * 60.0     # minutes per control
x_act *= V_TOT / (x_act @ dur)                    # exact volume
x_const = np.full(NC, V_TOT / dur.sum())


def design(params, coef, nh, x0s, whp_cap=True, ub=None):
    ln_r_, ln_P_, ln_ta_, ln_tau_ = params
    Sh, Pres10 = operators(coef, math.exp(ln_tau_))
    P = math.exp(ln_P_)
    sc = events_and_grad(x_const, Sh, P, ln_r_, ln_ta_, nh)[0]
    f = lambda x: events_and_grad(x, Sh, P, ln_r_, ln_ta_, nh)[0] / sc
    gf = lambda x: events_and_grad(x, Sh, P, ln_r_, ln_ta_, nh)[1] / sc
    cons = [dict(type="eq", fun=lambda x: (x @ dur - V_TOT) / V_TOT, jac=lambda x: dur / V_TOT)]
    if whp_cap:
        sub = slice(0, int(T_SI / dt10), 3)        # every 30 min
        Psub = Pres10[sub]
        Qsub = U10[sub]
        p0 = coef[0] - coef[1]
        cons.append(dict(type="ineq", fun=lambda x: (WHP_MAX - (p0 + coef[2] * (Qsub @ x) + coef[3] * (Qsub @ x) ** 2 + Psub @ x)) / WHP_MAX,
                         jac=lambda x: -((coef[2] + 2 * coef[3] * (Qsub @ x))[:, None] * Qsub + Psub) / WHP_MAX))
    best = None
    for x0 in x0s:
        bnds = [(0, Q_MAX if ub is None else ub[j]) for j in range(NC)]
        x0 = np.minimum(x0, [b[1] for b in bnds])
        r = minimize(f, x0, jac=gf, method="SLSQP", bounds=bnds, constraints=cons, options=dict(maxiter=500, ftol=1e-12))
        if best is None or r.fun < best.fun:
            best = r
    return best.x, best, Sh, P


params_full = (ln_r, ln_P, ln_ta, ln_tau)
nh_full, nh_si = NH, i_si


def N_of(x, params, coef, nh):
    Sh, _ = operators(coef, math.exp(params[3]))
    return events_and_grad(x, Sh, math.exp(params[1]), params[0], params[2], nh)[0]


x_opt, res_opt, Sh_full, _ = design(params_full, c_all, nh_full, [x_const, x_act])
x_opt_nocap, _, _, _ = design(params_full, c_all, nh_full, [x_const, x_act], whp_cap=False)
N_act, N_cst, N_opt = (N_of(x, params_full, c_all, nh_full) for x in (x_act, x_const, x_opt))
N_opt_nocap = N_of(x_opt_nocap, params_full, c_all, nh_full)
Pi = float((Sh_full @ x_const).max() / P_star)
# Pi-law (F3 table, Cycle 1): gain vs Pi, log-log interpolation
PI_T = np.array([0.19, 0.76, 1.53, 3.82, 7.64]); G_T = np.array([0.6, 2, 4.7, 18.4, 46.7])
gain_law = float(np.exp(np.interp(math.log(Pi), np.log(PI_T), np.log(G_T)))) if Pi <= PI_T[-1] else float("nan")
OUT["E4"] = dict(N_actual=N_act, N_constant=N_cst, N_optimal=N_opt, N_optimal_no_whp_cap=N_opt_nocap,
                 gain_vs_constant=1 - N_opt / N_cst, gain_vs_actual=1 - N_opt / N_act, actual_vs_constant=N_act / N_cst - 1,
                 Pi=Pi, Pi_law_gain_pct=gain_law, slsqp_success=bool(res_opt.success), slsqp_msg=str(res_opt.message))
print(f"E4 design (full window): actual {N_act:.0f}, constant {N_cst:.0f}, optimal {N_opt:.0f} (-{100 * (1 - N_opt / N_cst):.1f}% vs constant, "
      f"-{100 * (1 - N_opt / N_act):.1f}% vs actual; no WHP cap {N_opt_nocap:.0f}); Pi {Pi:.2f}, Pi-law gain {gain_law:.1f}%")

# E4b fairness: stage changes need pauses (stage 1->2, 2->3, 3->4 gaps from the stimulation log, local = UTC+3)
GAPS = [(11 * 24 + 21.56 - 3, 12 * 24 + 20.85 - 3), (19 * 24 + 0.05 - 3, 20 * 24 + 8.55 - 3), (27 * 24 + 7.77 - 3, 27 * 24 + 9.92 - 3)]
if SITE == "otn2":                                 # operator pauses >= 3 h taken from the pumping record (exploratory)
    act_ = np.where(q10[:int(T_SI / dt10)] > 0.01)[0] * dt10
    GAPS = [(float(a_), float(b_)) for a_, b_ in zip(act_[:-1], act_[1:]) if b_ - a_ >= 3.0]
    OUT["E4_pauses_h"] = GAPS
ub_gap = np.full(NC, Q_MAX)
for a_, b_ in GAPS:
    for j in range(NC):
        if j * DTC < b_ and (j + 1) * DTC > a_:
            ub_gap[j] = 0.0
x_cst_gap = np.where(ub_gap > 0, 1.0, 0.0); x_cst_gap *= V_TOT / (x_cst_gap @ dur)
x_opt_gap, _, _, _ = design(params_full, c_all, nh_full, [x_cst_gap], ub=ub_gap)
N_cst_gap, N_opt_gap = N_of(x_cst_gap, params_full, c_all, nh_full), N_of(x_opt_gap, params_full, c_all, nh_full)
OUT["E4"]["with_stage_pauses"] = dict(N_constant=N_cst_gap, N_optimal=N_opt_gap, actual_vs_constant=N_act / N_cst_gap - 1,
                                       gain_vs_actual=1 - N_opt_gap / N_act, gain_vs_constant=1 - N_opt_gap / N_cst_gap)
print(f"E4b with stage-change pauses: constant {N_cst_gap:.0f}, optimal {N_opt_gap:.0f}; actual is {100 * (N_act / N_cst_gap - 1):+.1f}% vs constant")

# E7 (EXPLORATORY, added after seeing E4; not pre-registered): why is the realised gain far below the Pi law?
# Hypothesis: long kernel memory (T >= window) makes the stress a function of cumulative volume (dose invariance),
# which no schedule can change; only the fast part of the stress is controllable. Test at FIXED Pi.
E7 = {}
for tag, keep in (("full", lambda T: True), ("fast_only(T<=30h)", lambda T: T <= 30), ("slow_only(T>=100h)", lambda T: T >= 100)):
    cmod = c_all.copy()
    for jj, T in enumerate(TJ):
        if not keep(T):
            cmod[4 + jj] = 0.0
    Sh_m, _ = operators(cmod, tau_d)
    scale = Pi * P_star / max((Sh_m @ x_const).max(), 1e-12)        # same Pi as the fitted model
    cmod[4:] *= scale
    xo, _, _, _ = design(params_full, cmod, nh_full, [x_const], whp_cap=False)
    Nc, No = N_of(x_const, params_full, cmod, nh_full), N_of(xo, params_full, cmod, nh_full)
    Na = N_of(x_act, params_full, cmod, nh_full)
    E7[tag] = dict(gain_opt_vs_const=1 - No / Nc, actual_vs_const=Na / Nc - 1)
    print(f"E7 {tag:20s}: Pi {Pi:.2f}: optimal vs constant -{100 * (1 - No / Nc):.1f}%, actual vs constant {100 * (Na / Nc - 1):+.1f}%")
OUT["E7_exploratory"] = E7

# E8 (EXPLORATORY, after E7 refuted the memory hypothesis): single-exponential kernels at FIXED Pi; the kernel time
# scale T_k relative to the injection duration decides how much an optimum can gain over constant rate, while Pi
# sets the penalty for rate fluctuations (actual cyclic schedule).
E8 = {}
for jj, T in enumerate(TJ):
    if T < 1:
        continue
    cmod = np.zeros_like(c_all); cmod[4 + jj] = 1.0
    params0 = (ln_r, ln_P, ln_ta, math.log(0.01))          # no extra delay
    Sh_m, _ = operators(cmod, 0.0)
    cmod[4:] *= Pi * P_star / max((Sh_m @ x_const).max(), 1e-12)
    xo, _, _, _ = design(params0, cmod, nh_full, [x_const], whp_cap=False)
    Nc, No, Na = (N_of(x, params0, cmod, nh_full) for x in (x_const, xo, x_act))
    E8[f"{T:g}"] = dict(Tk_over_Tinj=T / T_SI, gain_opt_vs_const=1 - No / Nc, actual_vs_const=Na / Nc - 1)
    print(f"E8 T_k {T:6g} h (T_k/T_inj {T / T_SI:7.4f}): optimal vs constant -{100 * (1 - No / Nc):5.1f}%, actual vs constant {100 * (Na / Nc - 1):+6.1f}%")
OUT["E8_exploratory"] = E8

# KKT / multistart check
x_rand = np.random.default_rng(1).uniform(0, Q_MAX, NC); x_rand *= V_TOT / (x_rand @ dur); x_rand = np.minimum(x_rand, Q_MAX)
x_ms, res_ms, _, _ = design(params_full, c_all, nh_full, [x_rand])
OUT["E4"]["multistart_rel_diff"] = abs(N_of(x_ms, params_full, c_all, nh_full) / N_opt - 1)

# E5 window artefact
x_opt_si, _, _, _ = design(params_full, c_all, nh_si, [x_const, x_act])
row = {}
for name, x in (("actual", x_act), ("constant", x_const), ("opt_full", x_opt), ("opt_shutin", x_opt_si)):
    row[name] = dict(N_shutin=N_of(x, params_full, c_all, nh_si), N_full=N_of(x, params_full, c_all, nh_full))
rank_si = sorted(row, key=lambda k: row[k]["N_shutin"]); rank_full = sorted(row, key=lambda k: row[k]["N_full"])
OUT["E5"] = dict(table=row, rank_shutin=rank_si, rank_full=rank_full, flip=rank_si != rank_full)
print("E5 window:", {k: (round(v['N_shutin']), round(v['N_full'])) for k, v in row.items()}, " rank@shut-in", rank_si, " rank@full", rank_full)

# E6 robustness: design with hold-out parameters, evaluate with full parameters
params_ho = tuple(fit_ho.x)
x_opt_ho, _, _, _ = design(params_ho, c_ho, nh_full, [x_const, x_act])
N_opt_ho_eval = N_of(x_opt_ho, params_full, c_all, nh_full)
kept = (N_cst - N_opt_ho_eval) / max(N_cst - N_opt, 1e-12)
OUT["E6"] = dict(N_design_holdout_eval_full=N_opt_ho_eval, gain_kept_fraction=kept, verdict="PASS" if kept >= 0.5 else ("WEAK" if kept > 0 else "FAIL"))
print(f"E6 robustness: hold-out design evaluated with full params {N_opt_ho_eval:.0f} -> keeps {100 * kept:.0f}% of the gain -> {OUT['E6']['verdict']}")

# ------------------------------------------------------------------ H1d D3 magnitudes
def p_exceed(N, m):
    return 1 - math.exp(-N * 10 ** (-b_val * (m - MC)))


OUT["D3"] = {nm: dict(N=Nv, P_ML_ge_2=p_exceed(Nv, 2.0), P_ML_ge_1_5=p_exceed(Nv, 1.5)) for nm, Nv in (("actual", N_act), ("constant", N_cst), ("optimal", N_opt))}
OUT["D3"]["observed_N"] = N_end_obs
OUT["D3"]["observed_P_ML_ge_2_from_catalogue_N"] = p_exceed(N_end_obs, 2.0)
# G-R roll-off check: b vs Mc and expected vs observed count above M_L 1.5 / 2.0 (unbounded G-R)
bvm = {}
for mc_ in (-0.95, -0.5, 0.0, 0.5, 1.0):
    m_ = ml[ml >= mc_]
    bvm[str(mc_)] = dict(n=int(len(m_)), b=float(math.log10(math.e) / (m_.mean() - mc_)) if len(m_) > 10 else None)
OUT["D3"]["b_vs_Mc"] = bvm
exp2 = N_end_obs * 10 ** (-b_val * (2.0 - MC)); exp15 = N_end_obs * 10 ** (-b_val * (1.5 - MC))
OUT["D3"]["GR_expected_ge2"] = exp2; OUT["D3"]["observed_ge2"] = int((ml >= 2.0).sum())
OUT["D3"]["GR_expected_ge1_5"] = exp15; OUT["D3"]["observed_ge1_5"] = int((ml >= 1.5).sum())
OUT["D3"]["poisson_p_zero_ge2"] = math.exp(-exp2)
print("D3 b vs Mc:", {k: (v['n'], round(v['b'], 2) if v['b'] else None) for k, v in bvm.items()},
      f"; unbounded G-R expects {exp2:.1f} events >= 2.0 (observed {OUT['D3']['observed_ge2']}, P(0) {math.exp(-exp2):.3f}), "
      f"{exp15:.1f} >= 1.5 (observed {OUT['D3']['observed_ge1_5']})")
# maximum volume with the hazard of the actual schedule (convex programme: max sum q s.t. N(q) <= N_act, bounds, WHP cap)
Sh, Pres10 = operators(c_all, tau_d)
fN = lambda x: events_and_grad(x, Sh, P_star, ln_r, ln_ta, nh_full)
cons = [dict(type="ineq", fun=lambda x: (N_act - fN(x)[0]) / N_act, jac=lambda x: -fN(x)[1] / N_act)]
r = minimize(lambda x: -(x @ dur) / V_TOT, x_opt, jac=lambda x: -dur / V_TOT, method="SLSQP", bounds=[(0, Q_MAX)] * NC, constraints=cons,
             options=dict(maxiter=500, ftol=1e-12))
OUT["D3"]["max_volume_same_hazard_m3"] = float(r.x @ dur)
OUT["D3"]["max_volume_ratio"] = float(r.x @ dur / V_TOT)
print(f"D3: b {b_val:.2f}; P(ML>=2): actual {OUT['D3']['actual']['P_ML_ge_2']:.3f}, constant {OUT['D3']['constant']['P_ML_ge_2']:.3f}, "
      f"optimal {OUT['D3']['optimal']['P_ML_ge_2']:.3f}; max volume at the actual hazard {r.x @ dur:.0f} m3 ({r.x @ dur / V_TOT:.2f}x)")

# ------------------------------------------------------------------ figure
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    th = np.arange(NH) + 0.5
    fig, ax = plt.subplots(4, 1, figsize=(10, 11), sharex=True)
    t10 = np.arange(n10) * dt10
    ax[0].plot(t10 / 24, q10 * 1000, lw=0.6, label="actual"); ax[0].set_ylabel("rate (l/min)")
    xc = np.repeat(x_opt, 1); tt = np.arange(NC) * DTC / 24
    ax[0].step(tt, x_opt * 1000, where="post", lw=1.2, color="C3", label="optimal (same volume)"); ax[0].legend()
    ax[1].plot(t10[ok] / 24, whp10[ok], ".", ms=1, label="WHP observed")
    ax[1].plot(t10 / 24, A_full @ c_all, lw=0.8, label="kernel fit"); ax[1].plot(t10 / 24, A_full @ c_tr, lw=0.8, ls="--", label="kernel fit on first 60%")
    ax[1].axvline(0.6 * T_SI / 24, color="k", lw=0.5); ax[1].set_ylabel("WHP (bar)"); ax[1].legend(fontsize=8)
    ax[2].bar(th / 24, counts, width=1 / 24, color="0.6", label=f"observed M_L >= {MC}")
    ax[2].plot(th / 24, lam, color="C3", lw=1, label="Dieterich fit"); ax[2].plot(th / 24, lam_ho, color="C2", lw=1, ls="--", label="hold-out forecast")
    ax[2].set_ylabel("events / h"); ax[2].set_yscale("log"); ax[2].set_ylim(0.1, None); ax[2].legend(fontsize=8)
    ax[3].plot(th / 24, cum_obs, "k", label="observed"); ax[3].plot(th / 24, cum_mod, "C3", label="fit"); ax[3].plot(th / 24, np.cumsum(lam_ho), "C2--", label="hold-out")
    for a in ax: a.axvline(T_SI / 24, color="k", ls=":", lw=0.8)
    ax[3].set_ylabel("cumulative events"); ax[3].set_xlabel(f"days since {T0_LABEL} 00:00 UTC"); ax[3].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, f"{TAG}_fig.png"), dpi=110)
except Exception as e:
    print("figure skipped:", e)

# curves for the paper figures (rendered by ResSim Studio's figure engine, tools/paper_figs.cpp)
try:
    th = np.arange(NH) + 0.5; t10 = np.arange(n10) * dt10; st = max(1, n10 // 3000)
    CUR = dict(th_d=(th / 24).tolist(), counts=np.asarray(counts, float).tolist(), lam=np.asarray(lam, float).tolist(), lam_ho=np.asarray(lam_ho, float).tolist(),
               cum_obs=np.asarray(cum_obs, float).tolist(), cum_mod=np.asarray(cum_mod, float).tolist(), cum_ho=np.cumsum(lam_ho).tolist(),
               t10_d=(t10[::st] / 24).tolist(), q_lpm=(np.asarray(q10)[::st] * 1000).tolist(), tc_d=(np.arange(NC) * DTC / 24).tolist(), qopt_lpm=(np.asarray(x_opt) * 1000).tolist(),
               whp_t_d=(t10[ok][::st] / 24).tolist(), whp=np.asarray(whp10)[ok][::st].tolist(), whp_fit=(A_full @ c_all)[::st].tolist(), whp_fit_tr=(A_full @ c_tr)[::st].tolist(),
               t_shutin_d=float(T_SI / 24), t_holdout_d=float(0.6 * T_SI / 24))
    json.dump(CUR, open(os.path.join(HERE, f"{TAG}_curves.json"), "w"))
except Exception as e:
    print("curves skipped:", e)
OUT["runtime_s"] = time.time() - t_start
json.dump(OUT, open(os.path.join(HERE, f"{TAG}_results.json"), "w"), indent=1, default=float)
print(f"done in {OUT['runtime_s']:.0f} s")
