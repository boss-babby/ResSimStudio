// C28 production analysis: decline curves, IPR / VLP nodal analysis, VFPPROD generation (see prodanalysis.h)
#include "ressim/prodanalysis.h"
#include "ressim/fluid.h"
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace ressim {
namespace prod {

// ============================================================================ Arps decline
namespace {
// log-rate residual sum for given b, Di; qi by the closed-form least squares in log space
double arpsSse(const std::vector<double>& t, const std::vector<double>& lq, double b, double Di, double* lnQi) {
    const size_t n = t.size();
    double s = 0;
    std::vector<double> g(n);
    for (size_t i = 0; i < n; ++i) { g[i] = b > 1e-9 ? std::log1p(b * Di * t[i]) / b : Di * t[i]; s += lq[i] + g[i]; }
    const double a = s / double(n);
    double e = 0; for (size_t i = 0; i < n; ++i) { const double r = lq[i] - (a - g[i]); e += r * r; }
    if (lnQi) *lnQi = a;
    return e;
}
}  // namespace

ArpsFit fitArps(const std::vector<double>& tIn, const std::vector<double>& qIn, double bMax, double fixB) {
    ArpsFit F;
    // from the peak rate on, positive rates only
    size_t pk = 0; double qm = -1;
    for (size_t i = 0; i < qIn.size() && i < tIn.size(); ++i) if (qIn[i] > qm) { qm = qIn[i]; pk = i; }
    std::vector<double> t, lq, q;
    for (size_t i = pk; i < qIn.size() && i < tIn.size(); ++i) if (qIn[i] > 0 && std::isfinite(qIn[i])) { t.push_back(tIn[i]); q.push_back(qIn[i]); }
    if (t.size() < 3) { F.error = "at least 3 positive rates after the peak are needed"; return F; }
    F.t0 = t.front();
    for (auto& v : t) v -= F.t0;
    for (double v : q) lq.push_back(std::log(v));
    const double span = std::max(1e-6, t.back());
    double best = 1e300, bb = 0, bd = 0;
    const int nb = fixB >= 0 ? 1 : 41;
    for (int ib = 0; ib < nb; ++ib) {
        const double b = fixB >= 0 ? fixB : bMax * ib / (nb - 1);
        // golden-section on log Di over [1e-6 / span, 1e3 / span]
        double lo = std::log(1e-6 / span), hi = std::log(1e3 / span);
        // coarse scan first (the SSE can be multi-modal in Di for b > 1)
        double bl = lo, bf = 1e300;
        for (int k = 0; k <= 60; ++k) { const double x = lo + (hi - lo) * k / 60.0; const double f = arpsSse(t, lq, b, std::exp(x), nullptr); if (f < bf) { bf = f; bl = x; } }
        double a = std::max(lo, bl - (hi - lo) / 60.0), c = std::min(hi, bl + (hi - lo) / 60.0);
        const double gr = 0.6180339887498949;
        double x1 = c - gr * (c - a), x2 = a + gr * (c - a);
        double f1 = arpsSse(t, lq, b, std::exp(x1), nullptr), f2 = arpsSse(t, lq, b, std::exp(x2), nullptr);
        for (int it = 0; it < 80; ++it) {
            if (f1 < f2) { c = x2; x2 = x1; f2 = f1; x1 = c - gr * (c - a); f1 = arpsSse(t, lq, b, std::exp(x1), nullptr); }
            else { a = x1; x1 = x2; f1 = f2; x2 = a + gr * (c - a); f2 = arpsSse(t, lq, b, std::exp(x2), nullptr); }
        }
        const double x = f1 < f2 ? x1 : x2, f = std::min(f1, f2);
        if (f < best) { best = f; bb = b; bd = std::exp(x); }
    }
    double lnQi = 0; arpsSse(t, lq, bb, bd, &lnQi);
    F.qi = std::exp(lnQi); F.Di = bd; F.b = bb; F.points = int(t.size());
    F.rmseLog = std::sqrt(best / double(t.size()));
    double mean = 0; for (double v : q) mean += v; mean /= double(q.size());
    double sst = 0, sse = 0;
    for (size_t i = 0; i < q.size(); ++i) { const double m = arpsRate(F, t[i] + F.t0); sst += (q[i] - mean) * (q[i] - mean); sse += (q[i] - m) * (q[i] - m); }
    F.r2 = sst > 0 ? 1 - sse / sst : 1.0;
    F.ok = true;
    return F;
}

double arpsRate(const ArpsFit& f, double t) {
    const double dt = std::max(0.0, t - f.t0);
    if (f.b < 1e-9) return f.qi * std::exp(-f.Di * dt);
    return f.qi / std::pow(1 + f.b * f.Di * dt, 1.0 / f.b);
}
double arpsCum(const ArpsFit& f, double t) {
    const double dt = std::max(0.0, t - f.t0);
    if (f.Di <= 0) return f.qi * dt;
    const double q = arpsRate(f, t);
    if (f.b < 1e-9) return (f.qi - q) / f.Di;
    if (std::fabs(f.b - 1) < 1e-9) return f.qi / f.Di * std::log(f.qi / q);
    return std::pow(f.qi, f.b) / ((1 - f.b) * f.Di) * (std::pow(f.qi, 1 - f.b) - std::pow(q, 1 - f.b));
}
double arpsEUR(const ArpsFit& f, double qEcon, double tMax, double* tEnd) {
    double te = tMax;
    if (qEcon > 0 && qEcon < f.qi && f.Di > 0) {
        const double r = f.qi / qEcon;
        const double tq = f.b < 1e-9 ? std::log(r) / f.Di : (std::pow(r, f.b) - 1) / (f.b * f.Di);
        te = std::min(te, tq);
    } else if (qEcon >= f.qi) te = 0;
    if (tEnd) *tEnd = te;
    return arpsCum(f, f.t0 + te);
}

// ============================================================================ IPR
double iprRate(const Ipr& r, double pwf) {
    pwf = std::clamp(pwf, 0.0, r.pr);
    if (r.J <= 0 || r.pr <= 0) return 0.0;
    if (r.pb <= 0 || pwf >= r.pb) return r.J * (r.pr - pwf);
    if (r.pr <= r.pb) {   // saturated reservoir: Vogel with qmax = J pr / 1.8
        const double qmax = r.J * r.pr / 1.8, x = pwf / r.pr;
        return qmax * (1 - 0.2 * x - 0.8 * x * x);
    }
    const double qb = r.J * (r.pr - r.pb), qv = r.J * r.pb / 1.8, x = pwf / r.pb;
    return qb + qv * (1 - 0.2 * x - 0.8 * x * x);
}
double iprAof(const Ipr& r) { return iprRate(r, 0.0); }

// ============================================================================ VLP (homogeneous no-slip)
namespace {
struct UnitsSI { double p, len, qLiq, rs, rho, bgToRm3; };
UnitsSI unitsOf(const Case& c) {
    if (c.units == UnitSystem::Field) return {6894.757, 0.3048, 0.158987, 178.1076, 16.01846, 0.0056146};   // psi, ft, stb, Mscf/stb, lb/ft3, rb/Mscf
    return {1e5, 1.0, 1.0, 1.0, 1.0, 1.0};                                                                  // bar, m, sm3, sm3/sm3, kg/m3, rm3/sm3
}
}  // namespace

double vlpBhp(const Case& c, const Tubing& tb, double qOil, double thp, double wct, double gor) {
    FluidModel fl; fl.setup(c, 0);
    const UnitsSI u = unitsOf(c);
    const double D = tb.diameter * u.len, eps = tb.roughness * u.len, A = 0.25 * 3.141592653589793 * D * D;
    const double L = std::max(tb.lengthMD, tb.depthTVD) * u.len, H = tb.depthTVD * u.len;
    const double cosT = L > 0 ? H / L : 1.0;
    const int ns = std::max(4, tb.segments);
    const double dL = L / ns;
    wct = std::clamp(wct, 0.0, 0.999);
    const double qo = std::max(0.0, qOil) * u.qLiq / 86400.0;            // sm3/s
    const double qw = qo * wct / (1 - wct);
    const double R = std::max(0.0, gor) * u.rs;                           // sm3/sm3
    const DensitySC& rs = fl.rhoSC();
    const double ro = rs.oil * u.rho, rw = rs.water * u.rho, rg = rs.gas * u.rho;
    const double mdot = qo * ro + qw * rw + qo * R * rg;                   // kg/s
    double p = thp * u.p;
    for (int s = 0; s < ns; ++s) {
        // two-pass (predictor-corrector) gradient per segment
        double grad = 0;
        for (int pass = 0; pass < 2; ++pass) {
            const double pm = p + 0.5 * grad * dL;
            const double pc = std::max(pm / u.p, 1e-3);                  // case units for the PVT
            double rsat = fl.live() ? fl.rsSat(pc) : 0.0;
            rsat = std::min(rsat * u.rs, R);
            const double bo = fl.boS(pc, rsat / u.rs), bw = fl.bwS(pc);   // 1/B
            const double muo = bo / std::max(1e-30, fl.boMu(AD3(pc), AD3(rsat / u.rs)).v), muw = bw / std::max(1e-30, fl.bwMu(AD3(pc)).v);
            const double bg = fl.bgS(pc), mug = bg / std::max(1e-30, fl.bgMu(AD3(pc)).v);
            const double Qo = bo > 0 ? qo / bo : 0.0, Qw = bw > 0 ? qw / bw : 0.0;
            const double qgFree = qo * std::max(0.0, R - rsat);
            const double Qg = bg > 0 ? qgFree / (bg / u.bgToRm3) : 0.0;    // 1/Bg in case units -> rm3 / sm3
            const double Qt = Qo + Qw + Qg;
            if (Qt <= 0) { grad = 9.80665 * (ro * 0.8 + rw * 0.2) * cosT; break; }
            const double rho = mdot / Qt;
            const double mu = ((Qo * muo + Qw * muw + Qg * mug) / Qt) * 1e-3;   // cP -> Pa s
            const double v = Qt / A;
            const double Re = rho * v * D / std::max(mu, 1e-9);
            double f;
            if (Re < 1e-6) f = 0;
            else if (Re < 2300) f = 64.0 / Re;
            else { const double x = -1.8 * std::log10(std::pow(eps / D / 3.7, 1.11) + 6.9 / Re); f = 1.0 / (x * x); }
            grad = rho * 9.80665 * cosT + f * rho * v * v / (2 * D);
        }
        p += grad * dL;
    }
    return p / u.p;
}

OperatingPoint operatingPoint(const Case& c, const Ipr& ipr, const Tubing& tb, double thp, double wct, double gor) {
    OperatingPoint o;
    const double aof = iprAof(ipr);
    if (aof <= 0) { o.note = "no inflow (PI or reservoir pressure is zero)"; return o; }
    // f(q) = VLP(q) - IPR^-1(q); IPR^-1 by bisection on pwf
    auto pwfOf = [&](double q) {
        double lo = 0, hi = ipr.pr;
        for (int it = 0; it < 80; ++it) { const double m = 0.5 * (lo + hi); if (iprRate(ipr, m) > q) lo = m; else hi = m; }
        return 0.5 * (lo + hi);
    };
    auto f = [&](double q) { return vlpBhp(c, tb, q, thp, wct, gor) - pwfOf(q); };
    const int nScan = 40;
    double qPrev = aof * 1e-4, fPrev = f(qPrev);
    if (fPrev > 0) { o.note = "the well does not flow at this THP (VLP above the IPR at every rate)"; return o; }
    // the stable operating point is the rightmost crossing (VLP rising through the IPR)
    double qa = -1, qb = -1;
    for (int k = 1; k <= nScan; ++k) {
        const double q = aof * (1e-4 + (1 - 1e-4) * k / nScan);
        const double fq = f(q);
        if (fPrev <= 0 && fq > 0) { qa = qPrev; qb = q; }
        qPrev = q; fPrev = fq;
    }
    if (qa < 0) { o.ok = true; o.q = aof; o.pwf = 0; o.note = "VLP below the IPR up to the open-flow potential"; return o; }
    for (int it = 0; it < 60; ++it) { const double m = 0.5 * (qa + qb); if (f(m) <= 0) qa = m; else qb = m; }
    o.q = 0.5 * (qa + qb); o.pwf = pwfOf(o.q); o.ok = true;
    return o;
}

VfpTable makeVfpProd(const Case& c, const Tubing& tb, int id, const std::vector<double>& rates, const std::vector<double>& thps,
                     const std::vector<double>& wcts, const std::vector<double>& gors) {
    VfpTable v; v.id = id; v.inj = false; v.datum = tb.depthTVD;
    v.flo = "OIL"; v.wfr = "WCT"; v.gfr = "GOR";
    v.floV = rates; v.thpV = thps; v.wfrV = wcts; v.gfrV = gors; v.alqV = {0.0};
    const size_t nf = rates.size(), nt = thps.size(), nw = wcts.size(), ng = gors.size();
    v.bhp.assign(nf * nt * nw * ng, 0.0);
    for (size_t it = 0; it < nt; ++it)
        for (size_t iw = 0; iw < nw; ++iw)
            for (size_t ig = 0; ig < ng; ++ig)
                for (size_t f = 0; f < nf; ++f)
                    v.bhp[((it * nw + iw) * ng + ig) * nf + f] = vlpBhp(c, tb, rates[f], thps[it], wcts[iw], gors[ig]);
    return v;
}

Ipr iprFromResults(const Case& c, const std::vector<SummaryPoint>& summary, int well, std::string* note) {
    Ipr r;
    FluidModel fl; fl.setup(c, 0);
    if (summary.empty() || well < 0) { if (note) *note = "no results"; return r; }
    // last point with the well producing
    for (size_t k = summary.size(); k-- > 0;) {
        const SummaryPoint& S = summary[k];
        if (size_t(well) >= S.wells.size()) continue;
        const WellRates& W = S.wells[size_t(well)];
        if (!W.open || W.oil <= 0 || W.bhp <= 0) continue;
        r.pr = S.fpr;
        const double dd = S.fpr - W.bhp;
        if (dd <= 0) continue;
        double rsOil = 0;
        if (fl.live() && W.oil > 0) rsOil = W.gas / W.oil;
        r.pb = fl.live() ? fl.pbOfRs(rsOil) : 0.0;
        // invert the composite Vogel for J from the observed point (exact above pb, consistent below)
        Ipr t = r; t.J = 1.0;
        const double q1 = iprRate(t, W.bhp);
        r.J = q1 > 0 ? W.oil / q1 : 0.0;
        if (note) { char b[160]; std::snprintf(b, sizeof b, "from t = %.1f d: q_o %.4g at BHP %.4g, average pressure %.4g", S.time, W.oil, W.bhp, S.fpr); *note = b; }
        return r;
    }
    if (note) *note = "the well never produced oil";
    return r;
}

}  // namespace prod
}  // namespace ressim
