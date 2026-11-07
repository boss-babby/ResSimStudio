// Compositional FIM simulator core (milestones C2 / C3) - see compsim.h for the formulation.
#include "ressim/binio.h"
#include "ressim/compsim.h"
#include "ressim/steam.h"
#include "ressim/bsr.h"
#include "ressim/stepadjoint.h"
#include <map>
#include <memory>
#include <chrono>
#include "ressim/linsolve.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <chrono>
#include <cstdlib>
#include <stdexcept>

namespace ressim {
namespace comp {

namespace {
constexpr int MAXC = 12;            // components
constexpr int KD = MAXC + 2;        // derivative slots: p, N_1..N_nc, N_w
constexpr double MWW = 18.01528e-3; // kg/mol
constexpr double GRAV = 9.80665;

// forward-mode dual number with a fixed number of derivative slots (p, N_1..N_nc, N_w)
struct DV {
    double v = 0;
    double d[KD] = {};
    DV() = default;
    explicit DV(double x) : v(x) {}
};
inline DV operator+(const DV& a, const DV& b) { DV r(a.v + b.v); for (int k = 0; k < KD; ++k) r.d[k] = a.d[k] + b.d[k]; return r; }
inline DV operator-(const DV& a, const DV& b) { DV r(a.v - b.v); for (int k = 0; k < KD; ++k) r.d[k] = a.d[k] - b.d[k]; return r; }
inline DV operator*(const DV& a, const DV& b) { DV r(a.v * b.v); for (int k = 0; k < KD; ++k) r.d[k] = a.d[k] * b.v + a.v * b.d[k]; return r; }
inline DV operator/(const DV& a, const DV& b) {
    DV r(a.v / b.v); const double ib = 1.0 / b.v;
    for (int k = 0; k < KD; ++k) r.d[k] = (a.d[k] - r.v * b.d[k]) * ib;
    return r;
}
inline DV operator*(double s, const DV& a) { DV r(s * a.v); for (int k = 0; k < KD; ++k) r.d[k] = s * a.d[k]; return r; }
inline DV operator*(const DV& a, double s) { return s * a; }
inline DV operator/(double s, const DV& b) { return DV(s) / b; }
inline DV operator+(double s, const DV& a) { DV r = a; r.v += s; return r; }
inline DV operator-(double s, const DV& a) { DV r(s - a.v); for (int k = 0; k < KD; ++k) r.d[k] = -a.d[k]; return r; }
inline DV operator-(const DV& a) { return 0.0 - a; }
inline DV powDV(const DV& a, double e) {   // a^e for a >= 0
    if (a.v <= 0) return DV(0.0);
    DV r(std::pow(a.v, e)); const double f = e * std::pow(a.v, e - 1.0);
    for (int k = 0; k < KD; ++k) r.d[k] = f * a.d[k];
    return r;
}
inline DV sqrtDV(const DV& a) {
    if (a.v <= 0) return DV(0.0);
    DV r(std::sqrt(a.v)); const double f = 0.5 / r.v;
    for (int k = 0; k < KD; ++k) r.d[k] = f * a.d[k];
    return r;
}
inline DV clamp01(const DV& x) { if (x.v <= 0) return DV(0.0); if (x.v >= 1) return DV(1.0); return x; }
// piecewise-linear table lookup (constant extrapolation), with the segment slope as derivative
inline DV tab(const std::vector<double>& xs, const std::vector<double>& ys, const DV& x) {
    const size_t n = xs.size();
    if (n == 0 || ys.size() != n) return DV(0.0);
    if (x.v <= xs[0]) return DV(ys[0]);
    if (x.v >= xs[n - 1]) return DV(ys[n - 1]);
    const size_t i = size_t(std::upper_bound(xs.begin(), xs.end(), x.v) - xs.begin());
    const double s = (ys[i] - ys[i - 1]) / std::max(1e-300, xs[i] - xs[i - 1]);
    DV r(ys[i - 1] + s * (x.v - xs[i - 1]));
    for (int k = 0; k < KD; ++k) r.d[k] = s * x.d[k];
    return r;
}

struct Phase {
    bool on = false;
    DV vol;                 // phase volume in the cell (m3)
    DV mob;                 // molar mobility kr / mu * xi  (mol / (m3 Pa s))
    DV vmob;                // volumetric mobility kr / mu (1 / (Pa s))
    DV rho;                 // mass density (kg/m3)
    DV vm;                  // molar volume (m3/mol)
    DV pc;                  // phase pressure - oil pressure (Pa): -Pcow (water), +Pcgo (gas), 0 (oil); set for absent phases too
    DV x[MAXC + 1];         // composition (slot nc = water: 1 in the aqueous phase, 0 in the hydrocarbon phases)
};
struct CellEval { Phase ph[3]; DV vtot; };   // 0 = liquid (oil), 1 = vapour (gas), 2 = water

// PVTW water: molar volume (m3/mol) and its pressure derivative
inline void waterVm(const CompCase& c, double p, double& vm, double& dvm) {
    const double X = c.cw * (p - c.pwRef), den = 1 + X + 0.5 * X * X;
    const double bw = c.bwRef / den, dbw = -c.bwRef * c.cw * (1 + X) / (den * den);
    vm = MWW * bw / c.rhoWsc; dvm = MWW * dbw / c.rhoWsc;
}

// ---- CO2SOL helpers: ln H(p) (piecewise linear in ln p) with d lnH / dp
inline double solLnH(const CompCase& c, double p, double* dp = nullptr) {
    const auto& X = c.solLnP; const auto& Y = c.solLnH;
    if (X.size() < 2) { if (dp) *dp = 0; return 0; }
    const double lp = std::log(std::max(p, 1.0));
    size_t k = 0;
    if (lp <= X.front()) k = 0; else if (lp >= X.back()) k = X.size() - 2;
    else k = size_t(std::upper_bound(X.begin(), X.end(), lp) - X.begin()) - 1;
    const double s = (Y[k + 1] - Y[k]) / (X[k + 1] - X[k]);
    if (dp) *dp = s / std::max(p, 1.0);
    return Y[k] + s * (lp - X[k]);
}
constexpr double MWCO2 = 44.0095e-3;
constexpr double VPHI_CO2 = 3.7e-5;   // apparent molar volume of dissolved CO2, m3/mol (Garcia 2001, ~37 cm3/mol)
// HC-side CO2 fugacity for HC mole numbers n (ln f, Pa) and d ln f / d n_CO2 (value-only flash; used by the s solve)
inline double hcLnFugCo2(const CompCase& c, const eos::TParams& tp, int nc, double p, const double* n, std::vector<double>& K, double* dlnf_dnc) {
    const int ic = c.co2Comp;
    double Nt = 0; for (int i = 0; i < nc; ++i) Nt += std::max(n[i], 0.0);
    Nt = std::max(Nt, 1e-300);
    std::vector<double> z(static_cast<size_t>(nc));
    for (int i = 0; i < nc; ++i) z[size_t(i)] = std::max(n[i], 0.0) / Nt;
    if (z[size_t(ic)] <= 0) { if (dlnf_dnc) *dlnf_dnc = 1e300; return -1e300; }
    eos::FlashOptions fo; if (!K.empty()) fo.Kinit = &K;
    const eos::FlashResult f = eos::flash(c.fluid, tp, p, z, fo);
    eos::EvalOptions o; o.dn = true; o.dP = false;
    if (f.phases == 2) {
        K = f.K;
        // d ln f / d n_c at fixed p: two-phase -> through the flash sensitivities of the vapour composition
        eos::FlashSens S; eos::flashSensitivities(c.fluid, tp, p, z, f, S, false);
        const eos::PhaseEval& V = f.V;
        const double yc = f.y[size_t(ic)];
        double d = 0;
        for (int k = 0; k < nc; ++k) {
            const double a = (k == ic ? 1.0 / yc : 0.0) + (V.dlnphi_dn.size() == size_t(nc) * nc ? V.dlnphi_dn[size_t(ic) * nc + k] : 0.0);
            d += a * S.dy_dz[size_t(k) * nc + ic] / Nt;
        }
        if (dlnf_dnc) *dlnf_dnc = d;
        return std::log(std::max(yc, 1e-300)) + V.lnphi[size_t(ic)] + std::log(p);
    }
    K.clear();
    const eos::PhaseEval e = eos::evalPhase(c.fluid, tp, p, z, o);
    if (dlnf_dnc) *dlnf_dnc = 1.0 / std::max(n[ic], 1e-300) - 1.0 / Nt + e.dlnphi_dn[size_t(ic) * nc + ic] / Nt;
    return std::log(z[size_t(ic)]) + e.lnphi[size_t(ic)] + std::log(p);
}
// moles of CO2 dissolved in the aqueous phase: g(s) = ln f_HC(N - s e_c) - ln(s / (N_w + s)) - ln H(p) = 0, safeguarded
// Newton on s in (0, N_c); s0 = warm start (< 0: none)
inline double solveAqueousCo2(const CompCase& c, const eos::TParams& tp, int nc, double p, const double* N, std::vector<double> K, double s0) {
    const int ic = c.co2Comp;
    const double Nc = std::max(N[ic], 0.0), Nw = std::max(N[nc], 0.0);
    if (Nc <= 1e-300 || Nw <= 1e-300) return 0.0;
    const double lnH = solLnH(c, p);
    std::vector<double> n(N, N + nc);
    auto g = [&](double s, double* dg) {
        n[size_t(ic)] = Nc - s;
        double dl = 0; const double lf = hcLnFugCo2(c, tp, nc, p, n.data(), K, &dl);
        if (dg) *dg = -dl - 1.0 / s + 1.0 / (Nw + s);
        return lf - std::log(s / (Nw + s)) - lnH;
    };
    double lo = 0, hi = Nc;   // g(lo+) = +inf, g(hi-) = -inf
    double s = (s0 > 0 && s0 < Nc) ? s0 : std::min(0.5 * Nc, Nw * 0.01);
    if (!(s > 0)) s = 0.5 * Nc;
    for (int it = 0; it < 60; ++it) {
        double dg = 0; const double v = g(s, &dg);
        if (v > 0) lo = s; else hi = s;
        if (std::fabs(v) < 1e-14) break;
        double sn = (dg < 0 && std::isfinite(dg)) ? s - v / dg : -1;
        if (!(sn > lo && sn < hi)) sn = (lo > 0 ? std::sqrt(lo * hi) : 0.5 * (lo + hi));
        if (lo <= 0 && !(sn > 0)) sn = 0.5 * hi;
        if (std::fabs(sn - s) < 1e-15 * Nc) { s = sn; break; }
        s = sn;
    }
    return s;
}

// ---- saturation functions (tables or Corey) of one cell: SATNUM table t, IMBNUM imbibition table imb, hysteresis state
struct SatFn {
    const CompCase& c;
    const SatTables& t;
    const SatTables* imb = nullptr;     // imbibition table (IMBNUM), nullptr = none
    double shy = 0;                     // historical maximum gas saturation (lagged)
    double swco() const { return t.sw.empty() ? c.swc : t.sw.front(); }
    DV krw(const DV& Sw) const {
        if (!t.sw.empty()) return tab(t.sw, t.krw, Sw);
        return c.krw0 * powDV(clamp01((1.0 / std::max(1e-12, 1 - c.swc - c.sorw)) * (Sw - DV(c.swc))), c.nw);
    }
    DV pcow(const DV& Sw) const { return t.pcow.empty() || t.sw.empty() ? DV(0.0) : tab(t.sw, t.pcow, Sw); }
    DV krgDrain(const DV& Sg) const {
        if (!t.sg.empty()) return tab(t.sg, t.krg, Sg);
        return c.krg0 * powDV(clamp01((1.0 / std::max(1e-12, 1 - c.swc - c.sgc)) * (Sg - DV(c.sgc))), c.ng);
    }
    static double critical(const std::vector<double>& s, const std::vector<double>& kr) {   // largest S with kr = 0
        double v = s.empty() ? 0.0 : s.front();
        for (size_t i = 0; i < s.size() && i < kr.size(); ++i) { if (kr[i] <= 0) v = s[i]; else break; }
        return v;
    }
    double sgcrD() const { return t.sg.empty() ? c.sgc : critical(t.sg, t.krg); }
    double sgMaxD() const { return t.sg.empty() ? 1.0 - c.swc : t.sg.back(); }
    double sgcrI() const {
        if (c.trap.sgtMax > 0) return std::max(c.trap.sgtMax, sgcrD());
        if (imb && imb != &t && !imb->sg.empty()) return critical(imb->sg, imb->krg);
        return sgcrD();
    }
    bool hystActive() const { return c.hyst.on && sgcrI() > sgcrD() + 1e-9; }
    // Land / Killough trapped gas saturation of the scanning curve that starts at shy
    double trapped(double sh) const {
        const double d = sgcrD(), i = sgcrI(), mx = sgMaxD();
        if (!c.hyst.on || i <= d + 1e-9 || sh <= d) return std::min(sh, d);
        const double C = 1.0 / (i - d) - 1.0 / std::max(1e-9, mx - d);
        return d + (sh - d) / (1.0 + C * (sh - d));
    }
    DV krgImb(const DV& S) const {
        if (imb && imb != &t && !imb->sg.empty() && c.trap.sgtMax <= 0) return tab(imb->sg, imb->krg, S);
        // no table: drainage curve rescaled to start at the imbibition critical saturation
        const double d = sgcrD(), i = sgcrI(), mx = sgMaxD();
        if (S.v <= i) return DV(0.0);
        return krgDrain(DV(d) + ((mx - d) / std::max(1e-9, mx - i)) * (S - DV(i)));
    }
    DV krg(const DV& Sg) const {
        if (!hystActive() || Sg.v >= shy || shy <= sgcrD() + 1e-9) return krgDrain(Sg);
        const double d = sgcrD(), i = sgcrI(), mx = sgMaxD();
        const double C = 1.0 / (i - d) - 1.0 / std::max(1e-9, mx - d);
        const double sgt = d + (shy - d) / (1.0 + C * (shy - d));
        if (Sg.v <= sgt) return DV(0.0);
        if (c.hyst.model == 1) {   // Killough: imbibition curve mapped onto [sgt, shy], scaled to meet the drainage curve at shy
            const DV sn = DV(i) + ((mx - i) / std::max(1e-12, shy - sgt)) * (Sg - DV(sgt));
            const double kImax = krgImb(DV(mx)).v, kD = krgDrain(DV(shy)).v;
            const double sc = kImax > 0 ? kD / kImax : 0.0;
            return sc * krgImb(sn);
        }
        // Carlson / Land: krg = krg_drainage(free gas saturation), S_gf = (a + sqrt(a^2 + 4 a / C)) / 2, a = Sg - Sgt
        const DV a = Sg - DV(sgt);
        const DV sgf = 0.5 * (a + sqrtDV(a * a + (4.0 / C) * a));
        return krgDrain(DV(d) + sgf);
    }
    DV pcgo(const DV& Sg) const { return t.pcog.empty() || t.sg.empty() ? DV(0.0) : tab(t.sg, t.pcog, Sg); }
    DV krow(const DV& So) const {
        if (!t.so.empty()) return tab(t.so, t.krow, So);
        return c.kro0 * powDV(clamp01((1.0 / std::max(1e-12, 1 - c.swc - c.sorw)) * (So - DV(c.sorw))), c.now);
    }
    DV krog(const DV& So) const {
        if (!t.so.empty()) return tab(t.so, t.krog, So);
        return c.kro0 * powDV(clamp01((1.0 / std::max(1e-12, 1 - c.swc - c.sorg)) * (So - DV(c.sorg + c.swc))), c.no);
    }
    // three-phase oil relative permeability
    DV kro(const DV& So, const DV& Sw, const DV& Sg, bool water) const {
        const double sc = swco();
        const DV kw = krow(So), kg = krog(So);
        if (c.kroModel == 1) {   // Stone II (normalised by krocw = krow(Swc))
            const double krocw = std::max(1e-12, krow(DV(1.0 - sc)).v);
            const DV krwv = water ? krw(Sw) : DV(0.0), krgv = krg(Sg);
            DV r = krocw * (((1.0 / krocw) * kw + krwv) * ((1.0 / krocw) * kg + krgv) - krwv - krgv);
            if (r.v < 0) r = DV(0.0);
            return r;
        }
        DV dw = Sw - DV(sc); if (dw.v < 0) dw = DV(0.0);   // ECLIPSE default model
        const DV den = Sg + dw;
        return den.v > 1e-12 ? (Sg * kg + dw * kw) / den : kw;
    }
};
inline const SatTables& tableOf(const CompCase& c, int k) { return k > 0 && k <= int(c.satRegions.size()) ? c.satRegions[size_t(k - 1)] : c.sat; }
inline SatFn satFn(const CompCase& c, int q, const std::vector<double>* shy) {
    const int r = c.satnum.empty() || q < 0 ? 0 : c.satnum[size_t(q)];
    const int ri = c.imbnum.empty() || q < 0 ? r : c.imbnum[size_t(q)];
    const SatTables& t = tableOf(c, r);
    return SatFn{c, t, ri != r ? &tableOf(c, ri) : nullptr, shy && q >= 0 && size_t(q) < shy->size() ? (*shy)[size_t(q)] : 0.0};
}

}  // namespace

double gasRelPerm(const CompCase& c, int table, int imbTable, double sg, double sgMax, double* trapped) {
    const SatTables& t = tableOf(c, table);
    const SatFn sf{c, t, imbTable >= 0 && imbTable != table ? &tableOf(c, imbTable) : nullptr, sgMax};
    if (trapped) *trapped = sf.trapped(std::max(sgMax, sg));
    return sf.krg(DV(sg)).v;
}

namespace {
// ---- cell evaluation: properties and derivatives with respect to (p, N_1..N_nc, N_w) of one cell
struct EvalCtx {
    const CompCase* c; const eos::TParams* tp; int nc;
    CompStats* st;
    const co2::Fluid* cf;
    const std::vector<double>* shy = nullptr;   // hysteresis state (historical max Sg per cell)
    std::vector<double>* sAq = nullptr;          // CO2SOL warm start: aqueous fraction of the cell's CO2 (per cell)
};

// ---- CO2STORE cell evaluation: binary CO2 / H2O flash (co2brine), aqueous in slot 0, CO2-rich phase in slot 1
void evalCellCo2(const EvalCtx& E, int q, double p, const double* N, int& nph, CellEval& out) {
    const CompCase& c = *E.c;
    const co2::Fluid& F = *E.cf;
    const double n0 = std::max(N[0], 0.0), n1 = std::max(N[1], 0.0);
    const double Nt = std::max(n0 + n1, 1e-300), z = n0 / Nt;
    const co2::FlashResult f = co2::flash(F, p, z);
    E.st->flashes++;
    nph = f.phases;
    out = CellEval();
    DV Ntd(Nt); Ntd.d[1] = Ntd.d[2] = 1.0;
    DV Z(z); Z.d[1] = (1.0 - z) / Nt; Z.d[2] = -z / Nt;
    auto aq = [&](Phase& P, const DV& X, const DV& moles) {
        P.on = true;
        double dvp = 0, dvx = 0;
        const double v = co2::aqueousMolarVolume(F, p, X.v, &dvp, &dvx);
        DV vm(v); for (int k = 0; k < KD; ++k) vm.d[k] = (k == 0 ? dvp : 0.0) + dvx * X.d[k];
        P.vm = vm; P.vol = moles * vm;
        const DV xi = DV(1.0) / vm;
        P.x[0] = X; P.x[1] = 1.0 - X;
        P.rho = (F.mwBrinePerWater * (1.0 - X) + co2::MW_CO2 * X) * xi;
        P.vmob = DV(1.0 / F.muBrine);
        P.mob = (1.0 / F.muBrine) * xi;
    };
    auto gas = [&](Phase& P, const DV& Y, const DV& moles) {
        P.on = true;
        double dv = 0;
        const double v = co2::gasMolarVolume(F, p, &dv);
        DV vm(v); vm.d[0] = dv;
        P.vm = vm; P.vol = moles * vm;
        const DV xi = DV(1.0) / vm;
        P.x[0] = Y; P.x[1] = 1.0 - Y;
        P.rho = (co2::MW_CO2 * Y + co2::MW_H2O * (1.0 - Y)) * xi;
        double dmu = 0;
        const double mu = co2::co2Viscosity(F.T, co2::MW_CO2 / v, &dmu);
        DV mud(mu); mud.d[0] = dmu * (-co2::MW_CO2 / (v * v)) * dv;
        P.vmob = DV(1.0) / mud;
        P.mob = xi / mud;
    };
    if (f.phases == 2) {
        DV X(f.x); X.d[0] = f.dx_dP;
        DV Y(f.y); Y.d[0] = f.dy_dP;
        const DV beta = (Z - X) / (Y - X);
        aq(out.ph[0], X, (1.0 - beta) * Ntd);
        gas(out.ph[1], Y, beta * Ntd);
    } else if (f.gas) gas(out.ph[1], Z, Ntd);
    else aq(out.ph[0], Z, Ntd);
    out.vtot = DV(0.0);
    for (int a = 0; a < 2; ++a) if (out.ph[a].on) out.vtot = out.vtot + out.ph[a].vol;
    const DV Sw = out.ph[0].on ? out.ph[0].vol / out.vtot : DV(0.0);
    const DV Sg = out.ph[1].on ? out.ph[1].vol / out.vtot : DV(0.0);
    const SatFn sf = satFn(c, q, E.shy);
    out.ph[1].pc = sf.pcow(Sw) + sf.pcgo(Sg);
    if (out.ph[0].on) { const DV kr = sf.krw(Sw); out.ph[0].mob = kr * out.ph[0].mob; out.ph[0].vmob = kr * out.ph[0].vmob; }
    if (out.ph[1].on) { const DV kr = sf.krg(Sg); out.ph[1].mob = kr * out.ph[1].mob; out.ph[1].vmob = kr * out.ph[1].vmob; }
}
void evalCell(const EvalCtx& E, int q, double p, const double* N, std::vector<double>& K, int& nph, CellEval& out) {
    if (E.c->co2store) { K.clear(); evalCellCo2(E, q, p, N, nph, out); return; }
    const int nc = E.nc;
    const eos::Model& m = E.c->fluid;
    const CompCase& c = *E.c;
    // CO2SOL: split the CO2 moles into aqueous s and hydrocarbon N_c - s; the hydrocarbon phases are evaluated for the
    // hydrocarbon moles and their derivatives transformed with ds/d(p, N) (implicit function theorem) below
    const bool sol = c.co2sol && c.water && c.co2Comp >= 0 && c.co2Comp < nc;
    const int ic = sol ? c.co2Comp : -1;
    double sAq = 0, Nloc[MAXC + 1];
    const double* const Nfull = N;
    if (sol) {
        const double warm = (E.sAq && size_t(q) < E.sAq->size()) ? (*E.sAq)[size_t(q)] * std::max(N[ic], 0.0) : -1.0;
        sAq = solveAqueousCo2(c, *E.tp, nc, p, N, K, warm);
        for (int i = 0; i <= nc; ++i) Nloc[i] = N[i];
        Nloc[ic] = std::max(N[ic] - sAq, 0.0);
        if (E.sAq && size_t(q) < E.sAq->size()) (*E.sAq)[size_t(q)] = N[ic] > 0 ? sAq / N[ic] : 0.0;
        N = Nloc;
    }
    DV lnfD;   // ln f_CO2 of the hydrocarbon phases (CO2SOL), derivatives w.r.t. (p, hydrocarbon moles)
    auto lnFug = [&](const eos::PhaseEval& e0, const DV* X, const std::vector<double>& xv) {
        eos::EvalOptions o; o.dn = true; o.dP = true;
        const eos::PhaseEval e = (e0.dlnphi_dn.size() == size_t(nc) * nc && e0.dlnphi_dP.size() == size_t(nc)) ? e0 : eos::evalPhase(m, *E.tp, p, xv, o);
        DV r(std::log(std::max(X[ic].v, 1e-300)) + e.lnphi[size_t(ic)] + std::log(p));
        for (int k = 0; k < KD; ++k) {
            double s = X[ic].d[k] / std::max(X[ic].v, 1e-300);
            for (int j = 0; j < nc; ++j) s += e.dlnphi_dn[size_t(ic) * nc + j] * X[j].d[k];
            if (k == 0) s += e.dlnphi_dP[size_t(ic)] + 1.0 / p;
            r.d[k] = s;
        }
        return r;
    };
    double Nt = 0; for (int i = 0; i < nc; ++i) Nt += std::max(N[i], 0.0);
    Nt = std::max(Nt, 1e-300);
    std::vector<double> z(static_cast<size_t>(nc));
    for (int i = 0; i < nc; ++i) z[size_t(i)] = std::max(N[i], 0.0) / Nt;
    eos::FlashOptions fo;
    if (!K.empty()) fo.Kinit = &K;
    const eos::FlashResult f = eos::flash(m, *E.tp, p, z, fo);
    E.st->flashes++;
    if (K.empty() || f.phases == 1) E.st->stabilityTests++;
    DV Ntd(Nt); for (int j = 0; j < nc; ++j) Ntd.d[1 + j] = 1.0;
    out = CellEval();
    auto fillPhase = [&](Phase& P, const eos::PhaseEval& e, const DV* X, const DV& moles) {
        P.on = true;
        DV v(e.vm);
        for (int k = 0; k < KD; ++k) { double s = k == 0 ? e.dvm_dP : 0.0; for (int i = 0; i < nc; ++i) s += e.dvm_dx[size_t(i)] * X[i].d[k]; v.d[k] = s; }
        P.vol = moles * v; P.vm = v;
        const DV xi = DV(1.0) / v;
        DV mw(0.0);
        for (int i = 0; i < nc; ++i) { P.x[i] = X[i]; mw = mw + m.comp[size_t(i)].mw * X[i]; }
        P.rho = mw * xi;
        std::vector<double> xv(static_cast<size_t>(nc)), dmx(static_cast<size_t>(nc));
        for (int i = 0; i < nc; ++i) xv[size_t(i)] = X[i].v;
        double dmxi = 0;
        const double mu = eos::lbcViscosity(m, xv, c.T, xi.v, dmx.data(), &dmxi);
        DV mud(mu);
        for (int k = 0; k < KD; ++k) { double s = dmxi * xi.d[k]; for (int i = 0; i < nc; ++i) s += dmx[size_t(i)] * X[i].d[k]; mud.d[k] = s; }
        P.vmob = DV(1.0) / mud;   // kr applied below
        P.mob = xi / mud;
    };
    if (f.phases == 2) {
        nph = 2; K = f.K;
        eos::FlashSens S;
        eos::flashSensitivities(m, *E.tp, p, z, f, S, false);
        DV beta(f.beta); beta.d[0] = S.dbeta_dP;
        DV X[MAXC], Y[MAXC];
        for (int j = 0; j < nc; ++j) beta.d[1 + j] = S.dbeta_dz[size_t(j)] / Nt;
        for (int i = 0; i < nc; ++i) {
            X[i] = DV(f.x[size_t(i)]); Y[i] = DV(f.y[size_t(i)]);
            X[i].d[0] = S.dx_dP[size_t(i)]; Y[i].d[0] = S.dy_dP[size_t(i)];
            for (int j = 0; j < nc; ++j) { X[i].d[1 + j] = S.dx_dz[size_t(i) * nc + j] / Nt; Y[i].d[1 + j] = S.dy_dz[size_t(i) * nc + j] / Nt; }
        }
        fillPhase(out.ph[0], f.L, X, (1.0 - beta) * Ntd);
        fillPhase(out.ph[1], f.V, Y, beta * Ntd);
        if (sol) lnfD = lnFug(f.V, Y, f.y);
    } else {
        nph = 1; K.clear();
        DV Z[MAXC];
        for (int i = 0; i < nc; ++i) { Z[i] = DV(z[size_t(i)]); for (int j = 0; j < nc; ++j) Z[i].d[1 + j] = ((i == j ? 1.0 : 0.0) - z[size_t(i)]) / Nt; }
        const int a = f.vapourLabel ? 1 : 0;
        fillPhase(out.ph[a], f.vapourLabel ? f.V : f.L, Z, Ntd);
        if (sol) lnfD = lnFug(f.vapourLabel ? f.V : f.L, Z, z);
    }
    double dsd[KD] = {};   // ds / d(p, N_1..N_nc, N_w)
    if (sol && sAq > 0) {
        const double Nw = std::max(Nfull[nc], 0.0);
        double dlnH = 0; solLnH(c, p, &dlnH);
        const double gs = -lnfD.d[1 + ic] - 1.0 / sAq + 1.0 / (Nw + sAq);
        if (gs < 0 && std::isfinite(gs)) {
            dsd[0] = -(lnfD.d[0] - dlnH) / gs;
            for (int j = 0; j < nc; ++j) dsd[1 + j] = -lnfD.d[1 + j] / gs;
            dsd[1 + nc] = -(1.0 / (Nw + sAq)) / gs;
        }
        // hydrocarbon-phase quantities: Q(p, N_hc(p, N)), N_hc,c = N_c - s
        auto tr = [&](DV& Q) { const double qc = Q.d[1 + ic]; for (int k = 0; k < KD; ++k) Q.d[k] -= qc * dsd[k]; };
        for (int a = 0; a < 2; ++a) {
            Phase& P = out.ph[a];
            if (!P.on) continue;
            tr(P.vol); tr(P.mob); tr(P.vmob); tr(P.rho); tr(P.vm); tr(P.pc);
            for (int i = 0; i < nc; ++i) tr(P.x[i]);
        }
    }
    if (c.water && sol) {   // aqueous phase with dissolved CO2: moles N_w + s, x_c = s / (N_w + s)
        Phase& W = out.ph[2];
        W.on = true;
        double vw, dvw; waterVm(c, p, vw, dvw);
        const double Nw = std::max(Nfull[nc], 0.0);
        DV sD(sAq); for (int k = 0; k < KD; ++k) sD.d[k] = dsd[k];
        DV nwd(Nw); nwd.d[1 + nc] = 1.0;
        const DV na = nwd + sD;
        const DV xc = na.v > 0 ? sD / na : DV(0.0);
        DV vwd(vw); vwd.d[0] = dvw;
        const DV v = (1.0 - xc) * vwd + VPHI_CO2 * xc;
        W.vol = na * v; W.vm = v;
        W.rho = ((1.0 - xc) * MWW + MWCO2 * xc) / v;
        const DV mu = c.muw * (1.0 + 4.65 * powDV(xc, 1.0134));
        W.vmob = DV(1.0) / mu;
        W.mob = DV(1.0) / (mu * v);
        W.x[nc] = 1.0 - xc; W.x[ic] = xc;
    } else if (c.water) {   // aqueous phase: depends on p and N_w only
        Phase& W = out.ph[2];
        W.on = true;
        double vw, dvw; waterVm(c, p, vw, dvw);
        const double Nw = std::max(N[nc], 0.0);
        DV v(vw); v.d[0] = dvw;
        DV nwd(Nw); nwd.d[1 + nc] = 1.0;
        W.vol = nwd * v; W.vm = v;
        W.rho = DV(MWW) / v;
        W.vmob = DV(1.0 / c.muw);
        W.mob = DV(1.0 / c.muw) / v;
        W.x[nc] = DV(1.0);
    }
    if (c.waterVap && c.water && out.ph[1].on && out.ph[2].on) {
        // water vaporisation: free-water Raoult law y_w = Psat(T) / p; the water vapour takes the gas molar volume
        Phase& V = out.ph[1]; Phase& W = out.ph[2];
        const double ps = steam::psat(c.T);
        DV yw(ps / p); yw.d[0] = -ps / (p * p);
        if (yw.v < 0.95) {
            const DV nV = V.vol / V.vm, nAq = W.vol / W.vm;
            DV nwv = yw / (1.0 - yw) * nV;
            if (nwv.v > 0.999 * nAq.v) nwv = 0.999 * nAq;   // dry-out: (almost) all the water is in the vapour
            const DV nT = nV + nwv, fh = nV / nT;
            V.rho = (V.rho * V.vm * nV + MWW * nwv) / (nT * V.vm);
            for (int i = 0; i < nc; ++i) V.x[i] = V.x[i] * fh;
            V.x[nc] = nwv / nT;
            V.vol = nT * V.vm;
            W.vol = W.vol - nwv * W.vm;
        }
    }
    out.vtot = DV(0.0);
    for (int a = 0; a < 3; ++a) if (out.ph[a].on) out.vtot = out.vtot + out.ph[a].vol;
    auto applyKr = [&](Phase& P, const DV& kr) { P.mob = kr * P.mob; P.vmob = kr * P.vmob; };
    if (!c.water && c.sat.empty() && c.satRegions.empty() && c.kroModel == 0 && !c.hyst.on) {
        // two hydrocarbon phases, Corey curves on the hydrocarbon pore space (immobile connate water)
        const double span = std::max(1e-12, 1.0 - c.sgc - c.sorg);
        for (int a = 0; a < 2; ++a) {
            Phase& P = out.ph[a];
            if (!P.on) continue;
            const DV S = out.ph[1 - a].on ? P.vol / out.vtot : DV(1.0);
            DV sn = (1.0 / span) * (S - DV(a == 0 ? c.sorg : c.sgc));
            if (sn.v >= 1.0) sn = DV(1.0);
            applyKr(P, (a == 0 ? c.kro0 : c.krg0) * powDV(sn, a == 0 ? c.no : c.ng));
        }
        return;
    }
    // saturations on the total pore space (connate water fills 1 - fpv when water is immobile)
    const SatFn sf = satFn(c, q, E.shy);
    DV Sw, So, Sg;
    if (c.water) {
        Sw = out.ph[2].vol / out.vtot;
        So = out.ph[0].on ? out.ph[0].vol / out.vtot : DV(0.0);
        Sg = out.ph[1].on ? out.ph[1].vol / out.vtot : DV(0.0);
    } else {
        const double fh = 1.0 - c.swc;
        Sw = DV(c.swc);
        So = out.ph[0].on ? fh * (out.ph[0].vol / out.vtot) : DV(0.0);
        Sg = out.ph[1].on ? fh * (out.ph[1].vol / out.vtot) : DV(0.0);
    }
    out.ph[1].pc = sf.pcgo(Sg);
    if (c.water) {
        applyKr(out.ph[2], sf.krw(Sw));
        out.ph[2].pc = -sf.pcow(Sw);
    }
    if (out.ph[1].on) applyKr(out.ph[1], sf.krg(Sg));
    if (out.ph[0].on) applyKr(out.ph[0], sf.kro(So, Sw, Sg, c.water));
}

// ---- assembly
struct Sys {
    BSR J;
    std::vector<double> R, acc;
    std::vector<int> kab, kba;   // block index of (a, b) and (b, a) per connection
};
int findBlk(const BSR& J, int r, int c) { for (int k = J.rowPtr[size_t(r)]; k < J.rowPtr[size_t(r) + 1]; ++k) if (J.col[size_t(k)] == c) return k; return -1; }
void buildPattern(int n, int nb, const std::vector<CompSimulator::Conn>& conns, const std::vector<CompWell>& wells, Sys& s) {
    std::vector<std::vector<int>> nbr(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) nbr[size_t(i)].push_back(i);
    for (const CompWell& w : wells)   // the eliminated well equation couples all perforated cells of a well
        if (w.open)
            for (int a : w.cells) for (int b : w.cells) if (a != b) nbr[size_t(a)].push_back(b);
    for (const auto& cn : conns) { nbr[size_t(cn.a)].push_back(cn.b); nbr[size_t(cn.b)].push_back(cn.a); }
    s.J.n = n; s.J.nb = nb; s.J.rowPtr.assign(size_t(n) + 1, 0); s.J.col.clear();
    for (int i = 0; i < n; ++i) {
        auto& v = nbr[size_t(i)]; std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end());
        s.J.rowPtr[size_t(i) + 1] = s.J.rowPtr[size_t(i)] + int(v.size());
        s.J.col.insert(s.J.col.end(), v.begin(), v.end());
    }
    s.J.val.assign(size_t(s.J.rowPtr[size_t(n)]) * nb * nb, 0.0);
    s.J.finalizeDiag();
    s.kab.resize(conns.size()); s.kba.resize(conns.size());
    for (size_t e = 0; e < conns.size(); ++e) { s.kab[e] = findBlk(s.J, conns[e].a, conns[e].b); s.kba[e] = findBlk(s.J, conns[e].b, conns[e].a); }
}

struct AsmIn {
    const CompCase& c; const eos::TParams& tp; int nc;
    const std::vector<double>& pv0; const std::vector<CompSimulator::Conn>& conns;
    const std::vector<CompSimulator::WellAux>& aux;
    const co2::Fluid& cf;
    const std::vector<double>* shy = nullptr;       // hysteresis state
    const std::vector<double>* kMult = nullptr;     // permeability factor per cell (well indices)
    const std::vector<double>* pvRef = nullptr;     // reference pore volume (mineral reactive area)
    const std::vector<double>* minReact = nullptr;  // remaining reactive silicate (mol) per cell and reaction
    const std::vector<double>* minProd = nullptr;   // reacted moles per cell and reaction
    std::vector<double>* minRate = nullptr;         // out: reaction rates (mol/s) of the last assembly
    std::vector<double>* sAq = nullptr;             // CO2SOL warm start per cell
};
constexpr double RGAS = 8.314462618;

// injected-fluid molar density xi(p) and its derivative
void injXi(const AsmIn& in, const CompWell& w, double p, double& xi, double& dxi) {
    if (w.water) { double vw, dvw; waterVm(in.c, p, vw, dvw); xi = 1.0 / vw; dxi = -dvw / (vw * vw); return; }
    if (in.c.co2store) {   // injected stream: CO2-rich (z_CO2 >= 0.5) or brine
        double v, dv;
        if (!w.zinj.empty() && w.zinj[0] >= 0.5) v = co2::gasMolarVolume(in.cf, p, &dv);
        else { double dx; v = co2::aqueousMolarVolume(in.cf, p, 0.0, &dv, &dx); }
        xi = 1.0 / v; dxi = -dv / (v * v); return;
    }
    eos::EvalOptions o; o.dn = false; o.dP = true;
    const eos::PhaseEval e = eos::evalPhase(in.c.fluid, in.tp, p, w.zinj, o);
    xi = 1.0 / e.vm; dxi = -e.dvm_dP / (e.vm * e.vm);
}

// assemble R (n * nv) and the block Jacobian; acc = accumulation blocks (for the true-IMPES weights)
void assemble(const AsmIn& in, double dt, const std::vector<double>& p, const std::vector<double>& N, const std::vector<double>& NOld,
              std::vector<std::vector<double>>& K, std::vector<int>& nph, CompStats& st, Sys& s, std::vector<CompWellReport>* rep) {
    const CompCase& c = in.c;
    const int nc = in.nc, ncw = nc + (c.water ? 1 : 0), n = c.size(), nv = ncw + 1, b2 = nv * nv;
    EvalCtx E{&c, &in.tp, nc, &st, &in.cf, in.shy, in.sAq};
    const int icSol = (c.co2sol && c.water) ? c.co2Comp : -1;
    auto inPh = [&](int i, int al) { return al == 2 ? (i == nc || i == icSol) : (i != nc || (al == 1 && c.waterVap)); };
    std::vector<CellEval> ce(static_cast<size_t>(n));
    // cell evaluations (flash + phase properties + derivatives) are independent: thread-parallel, per-thread counters
    {
        int fl = 0, stb = 0;
#pragma omp parallel if (n > 64) reduction(+ : fl, stb)
        {
            CompStats ls; EvalCtx El = E; El.st = &ls;
#pragma omp for schedule(dynamic, 16)
            for (int q = 0; q < n; ++q) evalCell(El, q, p[size_t(q)], &N[size_t(q) * ncw], K[size_t(q)], nph[size_t(q)], ce[size_t(q)]);
            fl += ls.flashes; stb += ls.stabilityTests;
        }
        st.flashes += fl; st.stabilityTests += stb;
    }
    s.R.assign(size_t(n) * nv, 0.0);
    std::fill(s.J.val.begin(), s.J.val.end(), 0.0);
    s.acc.assign(size_t(n) * b2, 0.0);
    // accumulation + volume balance
    for (int q = 0; q < n; ++q) {
        double* R = &s.R[size_t(q) * nv];
        double* D = s.J.blk(s.J.diag[size_t(q)]);
        double* A = &s.acc[size_t(q) * b2];
        for (int i = 0; i < ncw; ++i) {
            R[i] = N[size_t(q) * ncw + i] - NOld[size_t(q) * ncw + i];
            D[i * nv + 1 + i] += 1.0; A[i * nv + 1 + i] = 1.0;
        }
        const double fpv = (c.water || c.co2store) ? 1.0 : 1.0 - c.swc;   // fluid-filled fraction modelled by the volume balance
        const double pvh = in.pv0[size_t(q)] * (1 + c.cr * (p[size_t(q)] - c.pref)) * fpv;
        const DV& vh = ce[size_t(q)].vtot;
        R[ncw] = pvh - vh.v;
        for (int k = 0; k < nv; ++k) {
            const double d = (k == 0 ? in.pv0[size_t(q)] * c.cr * fpv : 0.0) - vh.d[k];
            D[ncw * nv + k] += d; A[ncw * nv + k] = d;
        }
    }
    // inter-cell fluxes: phase potentials with capillary pressure and gravity, phase-wise upwinding
    for (size_t e = 0; e < in.conns.size(); ++e) {
        const int a = in.conns[e].a, bq = in.conns[e].b;
        const CellEval& A = ce[size_t(a)]; const CellEval& B = ce[size_t(bq)];
        const double gdd = c.gravity ? GRAV * in.conns[e].dD : 0.0;
        for (int al = 0; al < 3; ++al) {
            const Phase& PA = A.ph[al]; const Phase& PB = B.ph[al];
            if (!PA.on && !PB.on) continue;
            const double wa = PA.on ? (PB.on ? 0.5 : 1.0) : 0.0, wb = PB.on ? (PA.on ? 0.5 : 1.0) : 0.0;
            const double rho = wa * (PA.on ? PA.rho.v : 0.0) + wb * (PB.on ? PB.rho.v : 0.0);
            const double dphi = (p[size_t(a)] + PA.pc.v) - (p[size_t(bq)] + PB.pc.v) + rho * gdd;
            double dA[KD], dB[KD];
            for (int k = 0; k < nv; ++k) {
                dA[k] = (k == 0 ? 1.0 : 0.0) + PA.pc.d[k] + (PA.on ? wa * PA.rho.d[k] * gdd : 0.0);
                dB[k] = (k == 0 ? -1.0 : 0.0) - PB.pc.d[k] + (PB.on ? wb * PB.rho.d[k] * gdd : 0.0);
            }
            const bool upA = dphi >= 0;
            const Phase& U = upA ? PA : PB;
            if (!U.on) continue;
            const double T = in.conns[e].T;
            double* Raa = s.J.blk(s.J.diag[size_t(a)]); double* Rab = s.J.blk(s.kab[e]);
            double* Rbb = s.J.blk(s.J.diag[size_t(bq)]); double* Rba = s.J.blk(s.kba[e]);
            for (int i = 0; i < ncw; ++i) {
                if (!inPh(i, al)) continue;   // water only in the aqueous phase, hydrocarbons in oil / gas (+ CO2 in water: CO2SOL)
                const DV mx = U.mob * U.x[i];
                const double F = T * mx.v * dphi;
                s.R[size_t(a) * nv + i] += dt * F;
                s.R[size_t(bq) * nv + i] -= dt * F;
                for (int k = 0; k < nv; ++k) {
                    const double dFa = T * (mx.v * dA[k] + (upA ? mx.d[k] * dphi : 0.0));
                    const double dFb = T * (mx.v * dB[k] + (upA ? 0.0 : mx.d[k] * dphi));
                    Raa[i * nv + k] += dt * dFa; Rab[i * nv + k] += dt * dFb;
                    Rba[i * nv + k] -= dt * dFa; Rbb[i * nv + k] -= dt * dFb;
                }
            }
            // C17 mechanical dispersion (DISPERC, isotropic alpha_L = alpha_T): J_i = (alpha / d) |F_alpha| (x_i,a - x_i,b),
            // F_alpha = T lambda_alpha dPhi the molar phase flux, i.e. alpha |v| xi A / d (x_a - x_b) with the TPFA face velocity
            const double Tdp = in.conns[e].Tdisp;
            if (Tdp > 0 && PA.on && PB.on) {
                const double Fm = T * U.mob.v * dphi, aF = std::fabs(Fm), sg = Fm >= 0 ? 1.0 : -1.0;
                if (aF > 0) {
                    double dFmA[KD], dFmB[KD];
                    for (int k = 0; k < nv; ++k) {
                        dFmA[k] = T * (U.mob.v * dA[k] + (upA ? U.mob.d[k] * dphi : 0.0));
                        dFmB[k] = T * (U.mob.v * dB[k] + (upA ? 0.0 : U.mob.d[k] * dphi));
                    }
                    for (int i = 0; i < ncw; ++i) {
                        if (!inPh(i, al)) continue;
                        const double dx = PA.x[i].v - PB.x[i].v;
                        const double F = Tdp * aF * dx;
                        s.R[size_t(a) * nv + i] += dt * F;
                        s.R[size_t(bq) * nv + i] -= dt * F;
                        for (int k = 0; k < nv; ++k) {
                            const double dFa = Tdp * (sg * dFmA[k] * dx + aF * PA.x[i].d[k]);
                            const double dFb = Tdp * (sg * dFmB[k] * dx - aF * PB.x[i].d[k]);
                            Raa[i * nv + k] += dt * dFa; Rab[i * nv + k] += dt * dFb;
                            Rba[i * nv + k] -= dt * dFa; Rbb[i * nv + k] -= dt * dFb;
                        }
                    }
                }
            }
        }
    }
    // molecular diffusion (DIFFUSE): per phase present on both sides, J_i = Td D_i <S xi> (x_i,a - x_i,b), Td = phi A / d,
    // arithmetic mean of S xi (x Millington-Quirk phi^(1/3) S^(7/3) with tortuosity 1)
    const bool diff = c.co2store ? c.trap.diffusion : (!c.diffOil.empty() || !c.diffGas.empty());
    if (diff) {
        auto coef = [&](int al, int i) -> double {
            if (c.co2store) return al == 0 ? (i == 0 ? c.trap.dCo2Water : c.trap.dH2oWater) : (i == 0 ? c.trap.dCo2Gas : c.trap.dH2oGas);
            const std::vector<double>& D = al == 0 ? c.diffOil : c.diffGas;
            return size_t(i) < D.size() ? D[size_t(i)] : 0.0;
        };
        const bool mq = c.trap.tortuosity == 1;
        auto weight = [&](const CellEval& C, int q, int al) {
            const Phase& P = C.ph[al];
            DV S = P.vol / C.vtot;
            if (!c.co2store && !c.water) S = (1.0 - c.swc) * S;
            if (mq) { const double phi = c.cellPoro.empty() ? 0.2 : c.cellPoro[size_t(q)]; return std::cbrt(phi) * powDV(S, 10.0 / 3.0) / P.vm; }
            return S / P.vm;
        };
        for (size_t e = 0; e < in.conns.size(); ++e) {
            const double Td = in.conns[e].Td;
            if (Td <= 0) continue;
            const int a = in.conns[e].a, bq = in.conns[e].b;
            const CellEval& A = ce[size_t(a)]; const CellEval& B = ce[size_t(bq)];
            double* Raa = s.J.blk(s.J.diag[size_t(a)]); double* Rab = s.J.blk(s.kab[e]);
            double* Rbb = s.J.blk(s.J.diag[size_t(bq)]); double* Rba = s.J.blk(s.kba[e]);
            for (int al = 0; al < 2; ++al) {
                const Phase& PA = A.ph[al]; const Phase& PB = B.ph[al];
                if (!PA.on || !PB.on) continue;
                const DV wA = weight(A, a, al), wB = weight(B, bq, al);
                const double wm = 0.5 * (wA.v + wB.v);
                for (int i = 0; i < nc; ++i) {
                    const double D = coef(al, i);
                    if (D <= 0) continue;
                    const double dx = PA.x[i].v - PB.x[i].v;
                    const double F = Td * D * wm * dx;
                    s.R[size_t(a) * nv + i] += dt * F;
                    s.R[size_t(bq) * nv + i] -= dt * F;
                    for (int k = 0; k < nv; ++k) {
                        const double dFa = Td * D * (0.5 * wA.d[k] * dx + wm * PA.x[i].d[k]);
                        const double dFb = Td * D * (0.5 * wB.d[k] * dx - wm * PB.x[i].d[k]);
                        Raa[i * nv + k] += dt * dFa; Rab[i * nv + k] += dt * dFb;
                        Rba[i * nv + k] -= dt * dFa; Rbb[i * nv + k] -= dt * dFb;
                    }
                }
            }
        }
    }
    // mineral trapping (CO2STORE): kinetic carbonation, r = A k(T) Sw (1 - Q/K), Q/K = mEq / m_CO2(aq), limited by the
    // remaining reactive silicate (precipitation) or the carbonate formed so far (dissolution)
    if (c.co2store && !c.trap.minerals.empty() && in.minRate && in.minReact && in.minProd && in.pvRef) {
        const size_t nm = c.trap.minerals.size();
        for (int q = 0; q < n; ++q) {
            for (size_t r = 0; r < nm; ++r) (*in.minRate)[size_t(q) * nm + r] = 0.0;
            const CellEval& C = ce[size_t(q)];
            if (!C.ph[0].on) continue;
            const DV X = C.ph[0].x[0];
            const DV m = X / ((1.0 - X) * MWW);
            const DV Sw = C.ph[0].vol / C.vtot;
            double* R = &s.R[size_t(q) * nv];
            double* D = s.J.blk(s.J.diag[size_t(q)]);
            for (size_t r = 0; r < nm; ++r) {
                const co2::MineralReaction& M = c.trap.minerals[r];
                const double kT = std::pow(10.0, M.logk25) * std::exp(-M.Ea / RGAS * (1.0 / c.T - 1.0 / 298.15));
                const double Ak = M.area * (*in.pvRef)[size_t(q)] * kT;
                DV qk = m.v > 1e-12 ? M.mEq / m : DV(10.0);
                if (qk.v > 10.0) qk = DV(10.0);
                DV rate = Ak * (Sw * (1.0 - qk));
                const double inv = rate.v > 0 ? (*in.minReact)[size_t(q) * nm + r] : (*in.minProd)[size_t(q) * nm + r];
                const double cap = std::max(0.0, inv) / dt;
                if (std::fabs(rate.v) > cap) rate = DV(rate.v > 0 ? cap : -cap);
                (*in.minRate)[size_t(q) * nm + r] = rate.v;
                R[0] += dt * M.nuCO2 * rate.v; R[1] += dt * M.nuH2O * rate.v;
                for (int k = 0; k < nv; ++k) { D[k] += dt * M.nuCO2 * rate.d[k]; D[nv + k] += dt * M.nuH2O * rate.d[k]; }
            }
        }
    }
    // wells. Every flowing term t (connection k, phase alpha for producers / total for injectors) has
    //   component outflow  q_i,t = B_i,t (Phi_t - bhp),   Phi_t = p_k + Pc_alpha,k - H_k   (H: lagged wellbore head)
    // producers: B_i,t = WI mob_alpha x_i,alpha; injectors: B_i,t = z_inj,i WI xi_inj(p_k) sum_alpha kr/mu (q < 0).
    // Rate control: sum_t Lw_t (Phi_t - bhp) = +-target, Lw = rate-weighted mobility -> bhp eliminated analytically,
    //   d bhp / d u_j = sum_{t in j} (dLw_t (Phi_t - bhp) + Lw_t dPhi_t) / sum Lw.
    if (rep) rep->resize(c.wells.size());
    const double vwsc = MWW / c.rhoWsc;
    for (size_t widx = 0; widx < c.wells.size(); ++widx) {
        const CompWell& w = c.wells[widx];
        CompWellReport* wr = rep ? &(*rep)[widx] : nullptr;
        if (wr) { wr->molar.assign(size_t(ncw), 0.0); wr->resv = 0; wr->bhp = 0; wr->onBhp = false; }
        if (!w.open || w.cells.empty()) continue;
        if (w.injector && w.water && !c.water) throw std::runtime_error("compositional core: water injector needs the mobile water phase");
        const CompSimulator::WellAux& ax = in.aux[widx];
        struct Term { size_t k; int q, al; DV B[MAXC + 1]; DV lw; DV lwc[7]; DV phi; DV vol; };
        std::vector<Term> terms;
        std::vector<double> zi(size_t(ncw), 0.0);
        if (w.injector) {
            if (w.water) zi[size_t(nc)] = 1.0;
            else { double zs = 0; for (int i = 0; i < nc; ++i) zs += w.zinj[size_t(i)]; for (int i = 0; i < nc; ++i) zi[size_t(i)] = w.zinj[size_t(i)] / zs; }
        }
        for (size_t k = 0; k < w.cells.size(); ++k) {
            const int q = w.cells[k];
            const double wi = (k < w.wi.size() ? w.wi[k] : 1.0) * (in.kMult && !in.kMult->empty() ? (*in.kMult)[size_t(q)] : 1.0);
            const double H = k < ax.H.size() ? ax.H[k] : 0.0;
            const CellEval& C = ce[size_t(q)];
            if (w.injector) {
                Term t; t.k = k; t.q = q; t.al = -1;
                DV lt(0.0); for (int al = 0; al < 3; ++al) if (C.ph[al].on) lt = lt + C.ph[al].vmob;
                double xi, dxi; injXi(in, w, p[size_t(q)], xi, dxi);
                DV X(xi); X.d[0] = dxi;
                t.vol = wi * lt;
                t.lw = t.vol * X;
                for (int i = 0; i < ncw; ++i) t.B[i] = zi[size_t(i)] * t.lw;
                t.phi = DV(p[size_t(q)] - H); t.phi.d[0] = 1.0;
                terms.push_back(t);
                continue;
            }
            for (int al = 0; al < 3; ++al) {
                const Phase& P = C.ph[al];
                if (!P.on) continue;
                Term t; t.k = k; t.q = q; t.al = al;
                for (int i = 0; i < ncw; ++i) t.B[i] = !inPh(i, al) ? DV(0.0) : wi * (P.mob * P.x[i]);
                const DV m = wi * P.mob;
                // surface yields of the phase: sum_i x_i,alpha * (per-component separator yield), lagged yields
                auto yield = [&](const std::vector<double>& f) { DV y(0.0); if (al < 2 && f.size() >= size_t(nc)) for (int i = 0; i < nc; ++i) y = y + f[size_t(i)] * P.x[i]; return y; };
                auto lwFor = [&](int ctl) -> DV {
                    switch (ctl) {
                        case 1: return al == 2 ? DV(0.0) : m;
                        case 2: return wi * P.vmob;
                        case 3: return al == 2 ? DV(0.0) : m * yield(ax.fo);
                        case 4: return al == 2 ? DV(0.0) : m * yield(ax.fg);
                        case 5: return c.co2store ? m * yield(ax.fw) : al == 2 ? vwsc * m : DV(0.0);
                        case 6: return c.co2store ? m * yield(ax.fw) : al == 2 ? vwsc * m : m * yield(ax.fo);
                        default: return DV(0.0);
                    }
                };
                t.lwc[w.control >= 0 && w.control <= 6 ? w.control : 0] = lwFor(w.control);
                for (const auto& l : w.limits) if (l.first >= 1 && l.first <= 6) t.lwc[l.first] = lwFor(l.first);
                t.lw = lwFor(w.control);
                t.vol = wi * P.vmob;
                t.phi = DV(p[size_t(q)] + P.pc.v - H);
                for (int kk = 0; kk < nv; ++kk) t.phi.d[kk] = (kk == 0 ? 1.0 : 0.0) + P.pc.d[kk];
                terms.push_back(t);
            }
        }
        const double sgn = w.injector ? -1.0 : 1.0;    // flow direction: producer Phi > bhp, injector bhp > Phi
        bool rateMode = w.injector ? !w.injBhpControl : w.control != 0;
        const double target = w.injector ? w.rate : w.target;
        // all potentials relative to the first term (avoids cancellation for large well indices: Phi - bhp << p)
        const double ref = terms.empty() ? 0.0 : terms[0].phi.v;
        std::vector<double> ph(terms.size());
        for (size_t t = 0; t < terms.size(); ++t) ph[t] = terms[t].phi.v - ref;
        double b = (w.injector ? w.bhpMax : w.bhp) - ref;   // bhp - ref
        std::vector<char> open(terms.size(), 1);
        bool shut = terms.empty();
        int activeCtl = w.injector ? -1 : w.control;   // rate constraint in force (producers: primary or a secondary limit)
        auto lwOf = [&](const Term& T, int ctl) -> const DV& { return ctl < 0 ? T.lw : T.lwc[ctl]; };
        if (rateMode && !shut) {
            bool zero = target <= 0;
            if (!w.injector) for (const auto& l : w.limits) if (l.first >= 1 && l.first <= 6 && l.second >= 0 && l.second <= 0) zero = true;
            if (zero) { rateMode = false; shut = true; }   // zero rate target / limit stops the well
            else {
                // every rate constraint gives the BHP that meets it exactly; producers run at the highest of them (the most
                // restrictive), so that no limit is exceeded
                std::vector<std::pair<int, double>> cands{{activeCtl, target}};
                if (!w.injector) for (const auto& l : w.limits) if (l.first >= 1 && l.first <= 6 && l.second > 0 && l.first != w.control) cands.push_back(l);
                bool found = false; double bBest = 0; int ctlBest = activeCtl; std::vector<char> openBest = open;
                for (const auto& cd : cands) {
                    std::vector<char> op(terms.size(), 1);
                    double bc = 0; bool okc = true;
                    for (int pass = 0; pass < 8; ++pass) {   // open set: terms flowing in the right direction
                        double sl = 0, slp = 0;
                        for (size_t t = 0; t < terms.size(); ++t) if (op[t]) { const double l = lwOf(terms[t], cd.first).v; sl += l; slp += l * ph[t]; }
                        if (sl <= 0) { okc = false; break; }
                        bc = (slp - sgn * cd.second) / sl;
                        bool changed = false;
                        for (size_t t = 0; t < terms.size(); ++t) { const char o = sgn * (ph[t] - bc) > 0; if (o != op[t]) { op[t] = o; changed = true; } }
                        if (!changed) break;
                    }
                    if (!okc) continue;
                    if (!found || (w.injector ? bc < bBest : bc > bBest)) { found = true; bBest = bc; ctlBest = cd.first; openBest = op; }
                }
                if (!found) rateMode = false;
                else { b = bBest; activeCtl = ctlBest; open = openBest; }
                // BHP limit: producer lower limit, injector upper limit (<= 0: none)
                if (!rateMode) b = (w.injector ? w.bhpMax : w.bhp) - ref;
                else if (!w.injector && b + ref < w.bhp) { rateMode = false; b = w.bhp - ref; }
                else if (w.injector && w.bhpMax > 0 && b + ref > w.bhpMax) { rateMode = false; b = w.bhpMax - ref; }
                if (!rateMode && w.injector && w.bhpMax <= 0) shut = true;
            }
        }
        if (!rateMode) for (size_t t = 0; t < terms.size(); ++t) open[t] = !shut && sgn * (ph[t] - b) > 0;
        if (wr) { wr->bhp = shut ? 0.0 : b + ref; wr->onBhp = !rateMode && !shut; wr->control = rateMode ? activeCtl : 0; }
        // d bhp / d (vars of each connection)
        const size_t nk = w.cells.size();
        std::vector<double> dbhp(nk * size_t(nv), 0.0);
        if (rateMode) {
            double sl = 0; for (size_t t = 0; t < terms.size(); ++t) if (open[t]) sl += lwOf(terms[t], activeCtl).v;
            for (size_t t = 0; t < terms.size(); ++t) {
                if (!open[t]) continue;
                const Term& T = terms[t];
                const DV& L = lwOf(T, activeCtl);
                for (int kk = 0; kk < nv; ++kk) dbhp[T.k * nv + kk] += (L.d[kk] * (ph[t] - b) + L.v * T.phi.d[kk]) / sl;
            }
        }
        for (size_t t = 0; t < terms.size(); ++t) {
            if (!open[t]) continue;
            const Term& T = terms[t];
            const double dp = ph[t] - b;
            double* R = &s.R[size_t(T.q) * nv];
            double* D = s.J.blk(s.J.diag[size_t(T.q)]);
            if (wr) wr->resv += T.vol.v * dp;
            for (int i = 0; i < ncw; ++i) {
                const DV& Bi = T.B[i];
                if (Bi.v == 0.0) { bool any = false; for (int kk = 0; kk < nv; ++kk) any = any || Bi.d[kk] != 0.0; if (!any) continue; }
                const double qi = Bi.v * dp;
                R[i] += dt * qi;
                if (wr) wr->molar[size_t(i)] += qi;
                for (int kk = 0; kk < nv; ++kk) D[i * nv + kk] += dt * (Bi.d[kk] * dp + Bi.v * T.phi.d[kk]);
                if (rateMode)
                    for (size_t j = 0; j < nk; ++j) {
                        const int qj = w.cells[j];
                        double* Bj = qj == T.q ? D : s.J.blk(findBlk(s.J, T.q, qj));
                        for (int kk = 0; kk < nv; ++kk) Bj[i * nv + kk] -= dt * Bi.v * dbhp[j * nv + kk];
                    }
            }
        }
    }
}

// phase volumes (liquid, vapour, water) of a cell at its current state; optional phase compositions
void phaseVolumes(const CompCase& c, const eos::TParams& tp, int nc, int ncw, double p, const double* N, const std::vector<double>& K, double v[3],
                  const co2::Fluid* cf = nullptr, double* xaq = nullptr, double* nwVap = nullptr) {
    if (nwVap) *nwVap = 0;
    if (c.co2store && cf) {
        const double Nt = std::max(N[0], 0.0) + std::max(N[1], 0.0), z = Nt > 0 ? std::max(N[0], 0.0) / Nt : 0.0;
        const co2::FlashResult f = co2::flash(*cf, p, z);
        v[0] = (1 - f.beta) * Nt * co2::aqueousMolarVolume(*cf, p, f.phases == 2 ? f.x : z);
        v[1] = f.beta * Nt * co2::gasMolarVolume(*cf, p);
        v[2] = 0;
        if (xaq) *xaq = f.phases == 2 ? f.x : (f.gas ? f.x : z);
        return;
    }
    double sA = 0, Nl[MAXC + 1];
    const double Nw0 = ncw > nc ? std::max(N[nc], 0.0) : 0.0;
    if (c.co2sol && ncw > nc && c.co2Comp >= 0 && c.co2Comp < nc) {   // CO2SOL: aqueous CO2 first
        sA = solveAqueousCo2(c, tp, nc, p, N, K, -1.0);
        for (int i = 0; i <= nc; ++i) Nl[i] = N[i];
        Nl[c.co2Comp] = std::max(0.0, N[c.co2Comp] - sA);
        N = Nl;
    }
    double Nt = 0; for (int i = 0; i < nc; ++i) Nt += std::max(N[i], 0.0);
    std::vector<double> z(static_cast<size_t>(nc));
    for (int i = 0; i < nc; ++i) z[size_t(i)] = std::max(N[i], 0.0) / std::max(Nt, 1e-300);
    eos::FlashOptions fo; if (!K.empty()) fo.Kinit = &K;
    const eos::FlashResult f = eos::flash(c.fluid, tp, p, z, fo);
    v[0] = v[1] = v[2] = 0;
    if (xaq) *xaq = Nw0 + sA > 0 ? sA / (Nw0 + sA) : 0.0;
    if (f.phases == 2) { v[0] = (1 - f.beta) * f.L.vm * Nt; v[1] = f.beta * f.V.vm * Nt; }
    else if (f.vapourLabel) v[1] = f.V.vm * Nt; else v[0] = f.L.vm * Nt;
    if (ncw > nc) { double vw, dvw; waterVm(c, p, vw, dvw); v[2] = Nw0 * vw + sA * VPHI_CO2;
        if (c.waterVap && v[1] > 0 && Nw0 > 0) {   // water vaporisation (as evalCell)
            const double yw = steam::psat(c.T) / p;
            if (yw < 0.95) {
                const double vmV = f.V.vm, nV = v[1] / vmV;
                const double nwv = std::min(yw / (1 - yw) * nV, 0.999 * Nw0);
                v[1] += nwv * vmV; v[2] -= nwv * vw;
                if (nwVap) *nwVap = nwv;
            }
        }
    }
}
}  // namespace

// ------------------------------------------------------------------------------------------------- setup
void CompSimulator::setup(const CompCase& c) {
    c_ = c;
    n_ = c.size(); nc_ = c.fluid.nc(); ncw_ = nc_ + (c.water ? 1 : 0); nv_ = ncw_ + 1;
    if (nc_ < 2 || nc_ > MAXC) throw std::runtime_error("compositional core: 2..12 components supported");
    if (nv_ > kMaxBlock) throw std::runtime_error("compositional core: block size exceeds kMaxBlock");
    tp_ = eos::atTemperature(c.fluid, c.T);
    if (c.co2store) {
        if (nc_ != 2 || c.water) throw std::runtime_error("CO2STORE: the model needs exactly the components CO2 and H2O (no separate water phase)");
        cf_ = co2::makeFluid(c.T, c.salinity);
        co2Vg_ = co2::gasSurfaceVolume(c.stdT, c.stdP);
        co2Vw_ = co2::brineSurfaceVolume(cf_, c.stdT, c.stdP);
    }
    c_.solLnP.clear(); c_.solLnH.clear();
    if (c.co2sol) {
        if (!c.water) throw std::runtime_error("CO2SOL needs the aqueous phase (WATER)");
        if (c_.co2Comp < 0)
            for (int i = 0; i < nc_; ++i) {
                std::string u = c.fluid.comp[size_t(i)].name;
                for (char& ch : u) ch = char(std::toupper(static_cast<unsigned char>(ch)));
                if (u == "CO2") c_.co2Comp = i;
            }
        if (c_.co2Comp < 0 || c_.co2Comp >= nc_) throw std::runtime_error("CO2SOL: no component named CO2 in the EOS model");
        cf_ = co2::makeFluid(c.T, c.solSalinity);
        // H(p) = f_CO2^pure(p) (this EOS) / x_sat(p) (Spycher-Pruess / Duan-Sun at the salinity): consistent saturation
        std::vector<double> e(size_t(nc_), 0.0); e[size_t(c_.co2Comp)] = 1.0;
        eos::EvalOptions o; o.dn = false; o.dP = false;
        const int np = 140; const double l0 = std::log(5e4), l1 = std::log(2e8);
        for (int k = 0; k < np; ++k) {
            const double lp = l0 + (l1 - l0) * k / (np - 1), P = std::exp(lp);
            double x = 0, y = 0;
            if (c.solTabP.size() >= 2) {   // SOLUBILI table: linear in pressure, clamped at the ends
                const auto& TP = c.solTabP; const auto& TX = c.solTabX;
                size_t j = 1; while (j + 1 < TP.size() && TP[j] < P) ++j;
                const double w = std::clamp((P - TP[j - 1]) / std::max(1e-9, TP[j] - TP[j - 1]), 0.0, 1.0);
                x = TX[j - 1] + w * (TX[j] - TX[j - 1]);
                if (P > TP.back()) x = TX.back();
            } else co2::solubility(cf_, P, x, y);
            x = std::max(x, 1e-12);
            const eos::PhaseEval pe = eos::evalPhase(c.fluid, tp_, P, e, o);
            c_.solLnP.push_back(lp); c_.solLnH.push_back(pe.lnphi[size_t(c_.co2Comp)] + lp - std::log(x));
        }
    }
    auto at = [&](const std::vector<double>& v, int i) { return v.size() == 1 ? v[0] : v[size_t(i)]; };
    pv0_.assign(size_t(n_), 0); depth_.assign(size_t(n_), 0);
    conns_.clear();
    if (!c.cellPV.empty()) {   // general cell list (deck grid)
        pv0_ = c.cellPV;
        depth_ = c.cellDepth.size() == c.cellPV.size() ? c.cellDepth : std::vector<double>(c.cellPV.size(), 0.0);
        for (const auto& cn : c.conns)
            if (cn.T > 0 && cn.a != cn.b && cn.a >= 0 && cn.b >= 0 && cn.a < n_ && cn.b < n_)
                conns_.push_back({cn.a, cn.b, cn.T, depth_[size_t(cn.b)] - depth_[size_t(cn.a)], cn.Td, cn.Tdisp});
    } else {
        std::vector<double> top(size_t(c.nx) * c.ny, c.top);
        for (int k = 0; k < c.nz; ++k)
            for (int j = 0; j < c.ny; ++j)
                for (int i = 0; i < c.nx; ++i) {
                    const int q = i + c.nx * (j + c.ny * k);
                    const double dzq = at(c.dz, q);
                    depth_[size_t(q)] = top[size_t(i + c.nx * j)] + 0.5 * dzq;
                    top[size_t(i + c.nx * j)] += dzq;
                    pv0_[size_t(q)] = at(c.dx, q) * at(c.dy, q) * dzq * at(c.poro, q);
                }
        // TPFA transmissibilities (harmonic averages of half-cell conductances), m3
        auto add = [&](int a, int b, char dir) {
            double ga, gb, ka, kb;   // geometric half-cell factors A / (d/2) and permeabilities
            if (dir == 'x') { ga = at(c.dy, a) * at(c.dz, a) / (0.5 * at(c.dx, a)); gb = at(c.dy, b) * at(c.dz, b) / (0.5 * at(c.dx, b)); ka = at(c.permx, a); kb = at(c.permx, b); }
            else if (dir == 'y') { ga = at(c.dx, a) * at(c.dz, a) / (0.5 * at(c.dy, a)); gb = at(c.dx, b) * at(c.dz, b) / (0.5 * at(c.dy, b)); ka = at(c.permy, a); kb = at(c.permy, b); }
            else { ga = at(c.dx, a) * at(c.dy, a) / (0.5 * at(c.dz, a)); gb = at(c.dx, b) * at(c.dy, b) / (0.5 * at(c.dz, b)); ka = at(c.permz, a); kb = at(c.permz, b); }
            const double ha = ka * ga, hb = kb * gb, da = at(c.poro, a) * ga, db = at(c.poro, b) * gb;
            const double dist = dir == 'x' ? 0.5 * (at(c.dx, a) + at(c.dx, b)) : dir == 'y' ? 0.5 * (at(c.dy, a) + at(c.dy, b)) : 0.5 * (at(c.dz, a) + at(c.dz, b));
            if (ha > 0 && hb > 0) conns_.push_back({a, b, ha * hb / (ha + hb), depth_[size_t(b)] - depth_[size_t(a)], da > 0 && db > 0 ? da * db / (da + db) : 0.0,
                                                    c.dispersivity > 0 && dist > 0 ? c.dispersivity / dist : 0.0});
        };
        for (int k = 0; k < c.nz; ++k)
            for (int j = 0; j < c.ny; ++j)
                for (int i = 0; i < c.nx; ++i) {
                    const int q = i + c.nx * (j + c.ny * k);
                    if (i + 1 < c.nx) add(q, q + 1, 'x');
                    if (j + 1 < c.ny) add(q, q + c.nx, 'y');
                    if (k + 1 < c.nz) add(q, q + c.nx * c.ny, 'z');
                }
    }
    T0_.resize(conns_.size()); for (size_t e = 0; e < conns_.size(); ++e) T0_[e] = conns_[e].T;
    // separator trains, ECLIPSE 300 FIELDSEP / SEPCOND semantics: the LAST stage is the stock tank (its liquid is the
    // stock-tank oil, at the stage conditions); no stage given = a single flash at standard conditions.
    // 0 = field separator, k = SEPCOND train k
    sepTp_.clear(); sepSt_.clear();
    for (size_t tr = 0; tr <= c.sepTrains.size(); ++tr) {
        std::vector<std::pair<double, double>> stages = tr == 0 ? c.separator : c.sepTrains[tr - 1];
        if (stages.empty()) stages.push_back({c.stdT, c.stdP});
        std::vector<eos::TParams> v;
        for (const auto& st : stages) v.push_back(eos::atTemperature(c.fluid, st.first));
        sepTp_.push_back(v); sepSt_.push_back(stages);
    }
    // initial state: flash (p, z) and fill the hydrocarbon pore volume
    p_.assign(size_t(n_), 0); N_.assign(size_t(n_) * ncw_, 0); K_.assign(size_t(n_), {}); nph_.assign(size_t(n_), 1);
    for (int q = 0; q < n_; ++q) {
        p_[size_t(q)] = at(c.p0, q);
        const std::vector<double>& z = c.z0.size() == 1 ? c.z0[0] : c.z0[size_t(q)];
        if (c.co2store) {   // overall composition fills the whole pore volume
            const double zs = std::max(1e-300, z[0] + z[1]), zc = std::max(1e-10, z[0] / zs);
            double v[3]; const double Nn[2] = {zc, 1 - zc};
            phaseVolumes(c, tp_, nc_, ncw_, p_[size_t(q)], Nn, K_[size_t(q)], v, &cf_);
            const double pvt = pv0_[size_t(q)] * (1 + c.cr * (p_[size_t(q)] - c.pref));
            const double Nt = pvt / std::max(1e-300, v[0] + v[1]);
            N_[size_t(q) * ncw_] = Nt * zc; N_[size_t(q) * ncw_ + 1] = Nt * (1 - zc);
            continue;
        }
        const eos::FlashResult f = eos::flash(c.fluid, tp_, p_[size_t(q)], z);
        double vmix;
        if (f.phases == 2) { vmix = (1 - f.beta) * f.L.vm + f.beta * f.V.vm; K_[size_t(q)] = f.K; nph_[size_t(q)] = 2; }
        else vmix = f.vapourLabel ? f.V.vm : f.L.vm;
        const double pvt = pv0_[size_t(q)] * (1 + c.cr * (p_[size_t(q)] - c.pref));
        double sw = c.water ? (c.sw0.empty() ? c.swc : at(c.sw0, q)) : c.swc;
        sw = std::min(sw, 1.0 - 1e-6);   // keep a trace of hydrocarbon in the water zone (the cell flash needs moles)
        const double Nt = pvt * (1 - sw) / vmix;
        double zs = 0; for (double v : z) zs += v;
        for (int i = 0; i < nc_; ++i) N_[size_t(q) * ncw_ + i] = std::max(Nt * z[size_t(i)] / zs, 1e-14 * Nt);
        if (c.water) { double vw, dvw; waterVm(c, p_[size_t(q)], vw, dvw); N_[size_t(q) * ncw_ + nc_] = pvt * sw / vw; }
        if (c_.co2sol && c.water) {   // CO2SOL: brine in equilibrium with the hydrocarbon CO2 fugacity
            const int ic = c_.co2Comp;
            std::vector<double> Kq = K_[size_t(q)];
            double dl = 0; const double lf = hcLnFugCo2(c_, tp_, nc_, p_[size_t(q)], &N_[size_t(q) * ncw_], Kq, &dl);
            const double x = std::min(0.2, std::exp(std::min(0.0, lf - solLnH(c_, p_[size_t(q)]))));
            if (x > 0 && z[size_t(ic)] > 0) {
                double vw, dvw; waterVm(c, p_[size_t(q)], vw, dvw);
                const double na = pvt * sw / ((1 - x) * vw + x * VPHI_CO2);
                N_[size_t(q) * ncw_ + nc_] = na * (1 - x);
                N_[size_t(q) * ncw_ + ic] += na * x;
            }
        }
    }
    pOld_ = p_; NOld_ = N_;
    sAq_.assign(size_t(n_), -1.0);
    cumInj_.assign(size_t(ncw_), 0); cumProd_.assign(size_t(ncw_), 0);
    rep_.assign(c.wells.size(), CompWellReport());
    t_ = 0; st_ = CompStats();
    // trapping / interaction state
    pvMult_.assign(size_t(n_), 1.0); kMult_.assign(size_t(n_), 1.0); pvEff_ = pv0_;
    sgHist_.assign(size_t(n_), 0.0);
    if (c.hyst.on) sgHist_ = gasSaturation();
    const size_t nm = c.co2store ? c.trap.minerals.size() : 0;
    minProd_.assign(size_t(n_) * nm, 0.0); minRate_.assign(size_t(n_) * nm, 0.0); minReact_.assign(size_t(n_) * nm, 0.0);
    for (int q = 0; q < n_; ++q)
        for (size_t r = 0; r < nm; ++r) {
            const co2::MineralReaction& M = c.trap.minerals[r];
            minReact_[size_t(q) * nm + r] = std::max(0.0, M.reactant0) * pv0_[size_t(q)] / std::max(1e-12, M.vmReactant);
        }
    saltDiss_.assign(size_t(n_), 0.0); saltSolid_.assign(size_t(n_), 0.0);
    if (c.co2store && c.trap.saltPrecip)
        for (int q = 0; q < n_; ++q) {
            double v[3], xa = 0; phaseVolumes(c_, tp_, nc_, ncw_, p_[size_t(q)], &N_[size_t(q) * ncw_], K_[size_t(q)], v, &cf_, &xa);
            const double nAq = v[0] > 0 ? v[0] / co2::aqueousMolarVolume(cf_, p_[size_t(q)], xa) : 0.0;
            saltDiss_[size_t(q)] = c.salinity * MWW * nAq * (1 - xa);
        }
    beginStep();
}

// ------------------------------------------------------------------------------------------------- separator
// single-phase separator stream: vapour by Li's label OR by its compressibility factor. At separator / stock-tank
// pressures Li's mixture critical temperature mislabels sub-critical-temperature gases (pure CO2 at 60 F and 1 atm,
// T < Tc but far below its vapour pressure) as liquid; a low-pressure root with Z > 0.5 is a gas.
static bool sepVapour(const eos::FlashResult& f) {
    if (f.vapourLabel) return true;
    const eos::PhaseEval& P = f.L.vm > 0 ? f.L : f.V;
    return P.Z > 0.5;
}
double CompSimulator::gasSurfaceVolumePerMole() const { return eos::Rgas * c_.stdT / c_.stdP; }
double CompSimulator::waterSurfaceVolumePerMole() const { return c_.co2store ? co2Vw_ : MWW / c_.rhoWsc; }

void CompSimulator::surfaceVolumes(const std::vector<double>& z, double& vo, double& vg, int train) const {
    const size_t tr = train > 0 && size_t(train) < sepTp_.size() ? size_t(train) : 0;
    const std::vector<eos::TParams>& sepTp = sepTp_[tr];
    const auto& stages = sepSt_[tr];
    vo = vg = 0;
    double zs = 0; for (int i = 0; i < nc_; ++i) zs += std::max(z[size_t(i)], 0.0);
    if (zs <= 0) return;
    std::vector<double> x(static_cast<size_t>(nc_));
    for (int i = 0; i < nc_; ++i) x[size_t(i)] = std::max(z[size_t(i)], 0.0) / zs;
    double nL = 1, gas = 0;
    const size_t ns = sepTp.size();
    for (size_t s = 0; s < ns && nL > 0; ++s) {
        const double P = stages[s].second;
        const eos::FlashResult f = eos::flash(c_.fluid, sepTp[s], P, x);
        if (f.phases == 2) { gas += nL * f.beta; nL *= 1 - f.beta; x = f.x; }
        else if (sepVapour(f)) { gas += nL; nL = 0; }
    }
    if (nL > 0) {
        eos::EvalOptions o; o.dn = false; o.dP = false; o.root = eos::Root::Liquid;
        vo = nL * eos::evalPhase(c_.fluid, sepTp.back(), stages.back().second, x, o).vm;
    }
    vg = gas * gasSurfaceVolumePerMole();
}

void CompSimulator::componentYields(const std::vector<double>& z, std::vector<double>& fo, std::vector<double>& fg, int train) const {
    const size_t tr = train > 0 && size_t(train) < sepTp_.size() ? size_t(train) : 0;
    const std::vector<eos::TParams>& sepTp = sepTp_[tr];
    const auto& stages = sepSt_[tr];
    fo.assign(size_t(nc_), 0.0); fg.assign(size_t(nc_), 0.0);
    double zs = 0; for (int i = 0; i < nc_; ++i) zs += std::max(z[size_t(i)], 0.0);
    if (zs <= 0) return;
    // feed with a trace of every component (yields defined for components absent from the stream)
    std::vector<double> x(static_cast<size_t>(nc_));
    for (int i = 0; i < nc_; ++i) x[size_t(i)] = std::max(z[size_t(i)], 0.0) / zs + 1e-8;
    double s = 0; for (double v : x) s += v; for (double& v : x) v /= s;
    const std::vector<double> feed = x;
    double nL = 1;
    const size_t ns = sepTp.size();
    for (size_t st = 0; st < ns && nL > 0; ++st) {
        const double P = stages[st].second;
        const eos::FlashResult f = eos::flash(c_.fluid, sepTp[st], P, x);
        if (f.phases == 2) { nL *= 1 - f.beta; x = f.x; }
        else if (sepVapour(f)) nL = 0;
    }
    double vmL = 0;
    if (nL > 0) {
        eos::EvalOptions o; o.dn = false; o.dP = false; o.root = eos::Root::Liquid;
        vmL = eos::evalPhase(c_.fluid, sepTp.back(), stages.back().second, x, o).vm;
    }
    const double vg = gasSurfaceVolumePerMole();
    for (int i = 0; i < nc_; ++i) {
        const double li = nL > 0 ? nL * x[size_t(i)] : 0.0;           // moles of i in the stock-tank oil per mole of feed
        const double fl = std::clamp(li / feed[size_t(i)], 0.0, 1.0);   // fraction of component i reporting to the oil
        fo[size_t(i)] = fl * vmL;   // sum_i feed_i fo_i = nL vmL
        fg[size_t(i)] = (1.0 - fl) * vg;
    }
}

// lagged well quantities of the coming time step: wellbore density (connection heads) and separator factors
void CompSimulator::beginStep() {
    aux_.resize(c_.wells.size());
    if (rep_.size() != c_.wells.size()) rep_.resize(c_.wells.size());
    EvalCtx E{&c_, &tp_, nc_, &st_, &cf_, &sgHist_};
    for (size_t w = 0; w < c_.wells.size(); ++w) {
        const CompWell& W = c_.wells[w];
        WellAux& A = aux_[w];
        const size_t nk = W.cells.size();
        A.H.assign(nk, 0.0);
        if (!W.open || nk == 0) continue;
        double pm = 0; for (int q : W.cells) pm += p_[size_t(q)] / double(nk);
        double rho = 0;
        std::vector<double> zs(size_t(nc_), 0.0);
        if (W.injector) {
            if (W.water) { double vw, dvw; waterVm(c_, pm, vw, dvw); rho = MWW / vw; }
            else if (c_.co2store) rho = (!W.zinj.empty() && W.zinj[0] >= 0.5) ? co2::co2Density(cf_, pm) : co2::aqueousDensity(cf_, pm, 0.0);
            else { eos::EvalOptions o; o.dn = false; o.dP = false; rho = eos::evalPhase(c_.fluid, tp_, pm, W.zinj, o).rho; }
        } else {
            double sm = 0, smr = 0;
            for (size_t k = 0; k < nk; ++k) {
                const int q = W.cells[k];
                std::vector<double> Kq = K_[size_t(q)]; int nphq = 0;
                CellEval C; evalCell(E, q, p_[size_t(q)], &N_[size_t(q) * ncw_], Kq, nphq, C);
                const double wi = (k < W.wi.size() ? W.wi[k] : 1.0) * (kMult_.empty() ? 1.0 : kMult_[size_t(q)]);
                for (int al = 0; al < 3; ++al) {
                    const Phase& P = C.ph[al];
                    if (!P.on) continue;
                    sm += P.vmob.v; smr += P.vmob.v * P.rho.v;
                    if (al < 2) for (int i = 0; i < nc_; ++i) zs[size_t(i)] += wi * P.mob.v * P.x[i].v;
                }
                if (k == 0 && sm <= 0) { for (int al = 0; al < 2; ++al) if (C.ph[al].on) { smr = C.ph[al].rho.v; sm = 1; break; } }
            }
            rho = sm > 0 ? smr / sm : 0.0;
            double zt = 0; for (double v : zs) zt += v;
            if (zt <= 0) {   // immobile: use the in-place composition of the first connection
                const int q = W.cells[0];
                for (int i = 0; i < nc_; ++i) zs[size_t(i)] = N_[size_t(q) * ncw_ + i];
            }
            if (c_.co2store) { A.fo = {0.0, 0.0}; A.fg = {co2Vg_, 0.0}; A.fw = {0.0, co2Vw_}; }
            else {
                bool needY = W.control >= 3;
                for (const auto& l : W.limits) needY = needY || l.first >= 3;
                if (needY) componentYields(zs, A.fo, A.fg, W.sep);
            }
        }
        if (c_.gravity) {
            const double dref = W.refDepth >= 0 ? W.refDepth : depth_[size_t(W.cells[0])];
            for (size_t k = 0; k < nk; ++k) A.H[k] = rho * GRAV * (depth_[size_t(W.cells[k])] - dref);
        }
    }
}

// ------------------------------------------------------------------------------------------------- public API
void CompSimulator::residualAndJacobian(double dt, std::vector<double>& R, std::vector<double>* Jd) {
    Sys s; buildPattern(n_, nv_, conns_, c_.wells, s);
    if (aux_.size() != c_.wells.size()) beginStep();
    AsmIn in{c_, tp_, nc_, pvEff_, conns_, aux_, cf_, &sgHist_, &kMult_, &pv0_, &minReact_, &minProd_, &minRate_, &sAq_};
    assemble(in, dt, p_, N_, NOld_, K_, nph_, st_, s, &rep_);
    R = s.R;
    if (Jd) {
        const int m = n_ * nv_;
        Jd->assign(size_t(m) * m, 0.0);
        for (int r = 0; r < n_; ++r)
            for (int k = s.J.rowPtr[size_t(r)]; k < s.J.rowPtr[size_t(r) + 1]; ++k) {
                const int cc = s.J.col[size_t(k)];
                const double* B = s.J.blk(k);
                for (int i = 0; i < nv_; ++i) for (int j = 0; j < nv_; ++j) (*Jd)[size_t(r * nv_ + i) * m + size_t(cc * nv_ + j)] = B[i * nv_ + j];
            }
    }
}

// ------------------------------------------------------------------------------------------------- group control
void CompSimulator::surfaceRates(size_t w, const CompWellReport& X, double& qo, double& qg, double& qw, double& qinj) const {
    qo = qg = qw = qinj = 0;
    if (X.molar.empty()) return;
    std::vector<double> zp(size_t(nc_), 0.0); double np = 0;
    for (int i = 0; i < nc_; ++i) { const double v = X.molar[size_t(i)]; if (v > 0) { zp[size_t(i)] = v; np += v; } else qinj -= v; }
    if (c_.co2store) {
        const double c0 = X.molar[0], w0 = X.molar[1];
        if (c0 > 0) qg = c0 * co2Vg_;
        if (w0 > 0) qw = w0 * co2Vw_;
        qinj = std::max(0.0, -c0) + std::max(0.0, -w0);
        return;
    }
    if (np > 0) { double vo, vg; surfaceVolumes(zp, vo, vg, w < c_.wells.size() ? c_.wells[w].sep : 0); qo = vo * np; qg = vg * np; }
    if (c_.water) { const double v = X.molar[size_t(nc_)]; if (v > 0) qw = v * waterSurfaceVolumePerMole(); else qinj -= v; }
}

double CompSimulator::groupQuantity(size_t w, const CompWellReport& X, const CompGroup& G) const {
    double qo, qg, qw, qi; surfaceRates(w, X, qo, qg, qw, qi);
    if (G.injector) return qi;
    switch (G.control) { case 2: return std::max(0.0, X.resv); case 3: return qo; case 4: return qg; case 5: return qw; case 6: return qo + qw; default: return qo; }
}

void CompSimulator::allocateGroups(const std::vector<char>& fixed, const std::vector<double>& fixRate) {
    for (const CompGroup& G : groups_) {
        std::vector<int> act;
        for (int w : G.wells) if (w >= 0 && size_t(w) < c_.wells.size() && c_.wells[size_t(w)].open && !c_.wells[size_t(w)].cells.empty()) act.push_back(w);
        if (act.empty()) continue;
        double rem = G.target;
        std::vector<int> free;
        for (int w : act) { if (fixed[size_t(w)]) rem -= fixRate[size_t(w)]; else free.push_back(w); }
        rem = std::max(0.0, rem);
        // guide rates: productivity x drawdown to the BHP limit x last yield of the target quantity per mole
        std::vector<double> g(free.size(), 0.0), f(free.size(), -1.0);
        double fsum = 0; int fn = 0;
        for (size_t k = 0; k < free.size(); ++k) {
            const CompWell& W = c_.wells[size_t(free[k])];
            double dd = 0;
            for (size_t j = 0; j < W.cells.size(); ++j) {
                const int q = W.cells[j];
                const double wi = (j < W.wi.size() ? W.wi[j] : 1.0) * kMult_[size_t(q)];
                dd += G.injector ? wi * std::max(0.0, (W.bhpMax > 0 ? W.bhpMax : 1.2 * p_[size_t(q)]) - p_[size_t(q)])
                                 : wi * std::max(0.0, p_[size_t(q)] - W.bhp);
            }
            g[k] = dd;
            const size_t w = size_t(free[k]);
            if (!G.injector && w < rep_.size() && !rep_[w].molar.empty()) {
                double nprod = 0; for (int i = 0; i < ncw_; ++i) nprod += std::max(0.0, rep_[w].molar[size_t(i)]);
                if (nprod > 0) { f[k] = groupQuantity(w, rep_[w], G) / nprod; fsum += f[k]; ++fn; }
            }
        }
        const double fmean = fn > 0 ? fsum / fn : 1.0;
        double gs = 0;
        for (size_t k = 0; k < free.size(); ++k) { if (!G.injector) g[k] *= f[k] >= 0 ? f[k] : fmean; gs += g[k]; }
        // GUIDERAT formula (producers with previous rates): potential of the nominated phase = current rate scaled to
        // the drawdown at the BHP limit; the guide rate is converted to the group's control phase by the well's ratio
        if (!G.injector && G.grPhase >= 0) {
            std::vector<double> gr(free.size(), -1.0);
            bool all = true;
            for (size_t k = 0; k < free.size(); ++k) {
                const size_t w = size_t(free[k]);
                if (w >= rep_.size() || rep_[w].molar.empty()) { all = false; break; }
                const CompWellReport& X = rep_[w];
                const CompWell& W = c_.wells[w];
                const double qo = X.qo, qg = X.qg, qw = X.qw;
                double qN = 0, r1 = 0, r2 = 0;
                switch (G.grPhase) {
                    case 1: qN = qo + qw; r1 = qN > 0 ? qw / qN : 0; r2 = qN > 0 ? qg / qN : 0; break;
                    case 2: qN = qg; r1 = qg > 0 ? qw / qg : 0; r2 = qg > 0 ? qo / qg : 0; break;
                    default: qN = qo; r1 = qo > 0 ? qw / qo : 0; r2 = qo > 0 ? qg / qo : 0; break;
                }
                if (!(qN > 0)) { all = false; break; }
                double ddLim = 0, ddCur = 0;
                for (size_t j = 0; j < W.cells.size(); ++j) {
                    const int q = W.cells[j];
                    const double wi = (j < W.wi.size() ? W.wi[j] : 1.0) * kMult_[size_t(q)];
                    ddLim += wi * std::max(0.0, p_[size_t(q)] - W.bhp); ddCur += wi * std::max(0.0, p_[size_t(q)] - X.bhp);
                }
                const double pot = qN * (ddCur > 0 ? std::max(1.0, ddLim / ddCur) : 1.0) * G.potScale;
                const double den = G.grB + G.grC * std::pow(std::max(0.0, r1 * G.r1Scale), G.grD) + G.grE * std::pow(std::max(0.0, r2 * G.r2Scale), G.grF);
                const double GR = std::pow(pot, G.grA) / std::max(1e-12, den);
                double nprod = 0; for (int i = 0; i < ncw_; ++i) nprod += std::max(0.0, X.molar[size_t(i)]);
                const double fq = nprod > 0 ? groupQuantity(w, X, G) / nprod : 0.0;     // control quantity per mole
                const double fn2 = nprod > 0 ? qN / nprod : 0.0;                         // nominated phase per mole
                gr[k] = fn2 > 0 ? GR * fq / fn2 : 0.0;
            }
            if (all) { gs = 0; for (size_t k = 0; k < free.size(); ++k) { g[k] = gr[k]; gs += g[k]; } }
        }
        for (size_t k = 0; k < free.size(); ++k) {
            CompWell& W = c_.wells[size_t(free[k])];
            const double share = rem * (gs > 0 ? g[k] / gs : 1.0 / double(free.size()));
            if (G.injector) { W.rate = share; W.injBhpControl = false; }
            else { W.control = G.control; W.target = share; }
        }
        for (int w : act) if (fixed[size_t(w)]) {
            CompWell& W = c_.wells[size_t(w)];
            if (G.injector) { W.rate = fixRate[size_t(w)]; W.injBhpControl = false; } else { W.control = G.control; W.target = fixRate[size_t(w)]; }
        }
    }
}

bool CompSimulator::advanceTo(double tEnd, double dtMax) {
    static const bool dbgEnv = std::getenv("RESSIM_COMP_DEBUG") != nullptr;
    dbg_ = dbgEnv;
    Sys s; buildPattern(n_, nv_, conns_, c_.wells, s);
    if (rep_.size() != c_.wells.size()) rep_.resize(c_.wells.size());
    const int m = n_ * nv_;
    AsmIn in{c_, tp_, nc_, pvEff_, conns_, aux_, cf_, &sgHist_, &kMult_, &pv0_, &minReact_, &minProd_, &minRate_, &sAq_};
    std::vector<CompWellReport> wr;
    std::vector<char> gFixed(c_.wells.size(), 0); std::vector<double> gFixRate(c_.wells.size(), 0.0);
    int gpass = 0;
    BlockCPR cpr;   // re-set up every Newton iteration; the AMG hierarchy (coarsening) is reused within a time step
    while (t_ < tEnd * (1 - 1e-12)) {
        double dt = std::min({dtLast_, dtMax, tEnd - t_});
        if (tEnd - t_ - dt < 1e-3 * dt) dt = tEnd - t_;
        if (forcedIdx_ < forced_.size()) { dt = std::min(forced_[forcedIdx_], tEnd - t_); if (tEnd - t_ - dt < 1e-6 * dt) dt = tEnd - t_; }
        pOld_ = p_; NOld_ = N_;
        updateRockFactors();
        if (!groups_.empty()) allocateGroups(gFixed, gFixRate);
        beginStep();
        const std::vector<std::vector<double>> Kold = K_;
        bool ok = false;
        int it = 0;
        double eHist[32] = {0};
        double thCycle = 1.0;   // C25 damping of a phase-state 2-cycle (see below)
        for (; it < (thCycle < 1.0 ? 25 : 15); ++it) {
            wr = rep_;
            auto tA = std::chrono::steady_clock::now();
            assemble(in, dt, p_, N_, NOld_, K_, nph_, st_, s, &wr);
            tAsm_ += std::chrono::duration<double>(std::chrono::steady_clock::now() - tA).count();
            // convergence: component residuals relative to the cell's moles, volume residual relative to its pore volume
            // + global material balance per component (sum of the residuals = in-place change - net inflow), E300 TOLMB style
            double eN = 0, eV = 0, eMB = 0;
            std::vector<double> sumR(static_cast<size_t>(ncw_), 0.0), sumN(static_cast<size_t>(ncw_), 0.0);
            double pvTot = 0, pvBad = 0, eBad = 0;   // C25: pore volume of the cells violating the strict tolerance
            for (int q = 0; q < n_; ++q) {
                double Nt = 0; for (int i = 0; i < ncw_; ++i) Nt += N_[size_t(q) * ncw_ + i];
                double eq = 0;
                for (int i = 0; i < ncw_; ++i) {
                    eq = std::max(eq, std::fabs(s.R[size_t(q) * nv_ + i]) / Nt);
                    sumR[size_t(i)] += s.R[size_t(q) * nv_ + i]; sumN[size_t(i)] += Nt;
                }
                const double ev = std::fabs(s.R[size_t(q) * nv_ + ncw_]) / pvEff_[size_t(q)];
                eN = std::max(eN, eq); eV = std::max(eV, ev);
                pvTot += pvEff_[size_t(q)];
                if (eq >= tolN_ || ev >= tolV_) { pvBad += pvEff_[size_t(q)]; eBad = std::max(eBad, std::max(eq, ev)); }
            }
            for (int i = 0; i < ncw_; ++i) eMB = std::max(eMB, std::fabs(sumR[size_t(i)]) / sumN[size_t(i)]);
            if (!std::isfinite(eN) || !std::isfinite(eV)) break;
            if (eN < tolN_ && eV < tolV_ && eMB < 1e-12) { ok = true; break; }
            // relaxed local tolerance after the strict iterations (OPM tolerance_cnv_relaxed / ECLIPSE XXXCNV style):
            // a 2-cycle in one cell (upwind / phase-state switching) with residual ~1e-6 of its moles is accepted as
            // long as the global material balance is strictly satisfied
            if (relaxedOk_ && it >= 4 && eN < 1e-5 && eV < 1e-5 && eMB < 1e-12) { ok = true; ++st_.relaxed; break; }
            // C25 (OPM relaxed_max_pv_fraction style): the cells still above the strict tolerance hold < 1e-4 of the pore
            // volume (tiny fracture / sliver cells at a phase boundary) and are within 1e-2 of their content, the global
            // material balance is strict -> accept. Without it one 0.03 m3 EDFM cell flipping phase state cuts the step to 1e-6 d.
            if (relaxedOk_ && it >= 3 && eMB < 1e-12 && eBad < 1e-2 && pvBad < 1e-4 * pvTot) { ok = true; ++st_.relaxed; break; }
            // early chop: Newton stagnating (residual not reduced 3x over the last 3 iterations while still far from
            // convergence) -> cutting now saves the remaining iterations of a step that would fail anyway
            eHist[it] = std::max(eN, eV);
            if (dbg_ && std::getenv("RESSIM_COMP_DEBUG")[0] == '2') {
                int qN = 0; double mN = 0;
                for (int q = 0; q < n_; ++q) { double Nt = 0; for (int i = 0; i < ncw_; ++i) Nt += N_[size_t(q) * ncw_ + i];
                    for (int i = 0; i < ncw_; ++i) if (std::fabs(s.R[size_t(q) * nv_ + i]) / Nt > mN) { mN = std::fabs(s.R[size_t(q) * nv_ + i]) / Nt; qN = q; } }
                std::fprintf(stderr, "    it %d eN %.3e (cell %d nph %d) eV %.3e eMB %.3e\n", it, eN, qN, nph_[size_t(qN)], eV, eMB);
            }
            // C25: phase-state 2-cycle (one cell flipping between one and two phases at the solubility limit - typically a
            // tiny-volume fracture / sliver cell next to a CO2 front): the iterates jump across the kink of the flash, the
            // residual repeats with period 2. Damping the update (halved per detected cycle, >= 1/16) moves the iterate into
            // the kink and converges linearly instead of cutting the time step to nothing.
            if (it >= 3 && std::fabs(eHist[it] - eHist[it - 2]) < 1e-2 * eHist[it] && std::fabs(eHist[it - 1] - eHist[it - 3]) < 1e-2 * eHist[it - 1] &&
                std::fabs(eHist[it] - eHist[it - 1]) > 0.1 * eHist[it] && thCycle > 1.0 / 16) { thCycle *= 0.5; ++st_.relaxed; }
            if (it >= 7 && eHist[it] > 1e-5 && eHist[it] > 0.3 * eHist[it - 3] && thCycle >= 1.0) break;
            auto tS = std::chrono::steady_clock::now();
            cpr.setup(s.J, s.acc.data(), it > 0);
            tSetup_ += std::chrono::duration<double>(std::chrono::steady_clock::now() - tS).count();
            tS = std::chrono::steady_clock::now();
            std::vector<double> rhs(static_cast<size_t>(m)), dx(static_cast<size_t>(m), 0.0);
            for (int r = 0; r < m; ++r) rhs[size_t(r)] = -s.R[size_t(r)];
            BsrOperator op(s.J);
            // inexact Newton: the linear tolerance (relative to the current residual) only needs to be well below the
            // Newton contraction; 1e-5 keeps quadratic-like convergence with far fewer Krylov iterations than 1e-8
            // C20: MPI domain-decomposed solve when a distributed solver is installed (ResSimCLI --mpi-dd); the block
            // system has pressure as the first unknown of every block, as the DD solver's CPR stage requires. Serial
            // GMRES + CPR on failure.
            LinResult lr;
            bool ddDone = false;
            if (distributedSolver()) {
                int itsD = 0; std::vector<double> xd;
                if (distributedSolver()->solve(s.J, rhs, xd, linTol_, 300, itsD) && xd.size() == dx.size()) {
                    dx.swap(xd); lr.converged = true; lr.iterations = itsD; lr.relResidual = linTol_; ddDone = true;
                }
            }
            if (!ddDone) lr = gmres(op, cpr, rhs, dx, linTol_, 300, 60);
            st_.linear += lr.iterations;
            tSolve_ += std::chrono::duration<double>(std::chrono::steady_clock::now() - tS).count();
            if (!lr.converged && lr.relResidual > 1e-3) break;
            // damped update: pressure change limited to 20 % per iteration, moles kept positive (component-wise)
            double th = thCycle;
            for (int q = 0; q < n_; ++q) { const double d = std::fabs(dx[size_t(q) * nv_]); if (d > 0.2 * p_[size_t(q)]) th = std::min(th, 0.2 * p_[size_t(q)] / d); }
            for (int q = 0; q < n_; ++q) {
                p_[size_t(q)] += th * dx[size_t(q) * nv_];
                double Nt = 0; for (int i = 0; i < ncw_; ++i) Nt += N_[size_t(q) * ncw_ + i];
                for (int i = 0; i < ncw_; ++i) {
                    double& v = N_[size_t(q) * ncw_ + i];
                    v = std::max(v + th * dx[size_t(q) * nv_ + 1 + i], 1e-14 * Nt);
                }
            }
            st_.newton++;
        }
        if (!ok) {
            if (dbg_) std::fprintf(stderr, "  cut at t %.4g d dt %.4g d after %d its\n", t_ / 86400, dt / 86400, it);
            p_ = pOld_; N_ = NOld_; K_ = Kold;
            dtLast_ = dt * 0.4; st_.cuts++; cutLast_ = true;
            if (dtLast_ < 1e-3 || forcedIdx_ < forced_.size()) return false;
            continue;
        }
        if (!groups_.empty() && gpass < 4) {   // group targets: members on their BHP limit are fixed, the rest re-distributed
            bool redo = false;
            for (const CompGroup& G : groups_) {
                double Q = 0; bool anyFree = false; std::vector<int> nowFixed;
                for (int w : G.wells) {
                    if (w < 0 || size_t(w) >= wr.size() || !c_.wells[size_t(w)].open) continue;
                    const double q = groupQuantity(size_t(w), wr[size_t(w)], G);
                    Q += q;
                    if (gFixed[size_t(w)]) continue;
                    if (wr[size_t(w)].onBhp) nowFixed.push_back(w); else anyFree = true;
                }
                if (Q < G.target * (1 - 1e-4) && !nowFixed.empty() && anyFree) {
                    for (int w : nowFixed) { gFixed[size_t(w)] = 1; gFixRate[size_t(w)] = groupQuantity(size_t(w), wr[size_t(w)], G); }
                    redo = true;
                }
            }
            if (redo) { p_ = pOld_; N_ = NOld_; K_ = Kold; ++gpass; continue; }
        }
        gpass = 0; std::fill(gFixed.begin(), gFixed.end(), 0);
        // well reports: component rates, surface rates (separator flash of the produced stream), cumulatives
        const double vgsc = gasSurfaceVolumePerMole(), vwsc = waterSurfaceVolumePerMole();
        for (size_t w = 0; w < wr.size(); ++w) {
            CompWellReport& R = rep_[w];
            const CompWellReport& X = wr[w];
            R.bhp = X.bhp; R.onBhp = X.onBhp; R.molar = X.molar; R.resv = X.resv;
            R.qo = R.qg = R.qw = R.qgi = R.qwi = 0;
            if (X.molar.empty()) continue;
            std::vector<double> zp(size_t(nc_), 0.0); double np = 0, ni = 0;
            for (int i = 0; i < nc_; ++i) { const double v = X.molar[size_t(i)]; if (v > 0) { zp[size_t(i)] = v; np += v; } else ni -= v; }
            if (c_.co2store) {   // CO2 -> surface gas, H2O -> brine
                const double c0 = X.molar[0], w0 = X.molar[1];
                if (c0 > 0) R.qg = c0 * co2Vg_; else R.qgi = -c0 * co2Vg_;
                if (w0 > 0) R.qw = w0 * co2Vw_; else R.qwi = -w0 * co2Vw_;
                ni = 0;
            } else if (np > 0) { double vo, vg; surfaceVolumes(zp, vo, vg, w < c_.wells.size() ? c_.wells[w].sep : 0); R.qo = vo * np; R.qg = vg * np; }
            if (!c_.co2store) R.qgi = ni * vgsc;
            if (c_.water) { const double v = X.molar[size_t(nc_)]; if (v > 0) R.qw = v * vwsc; else R.qwi = -v * vwsc; }
            R.cumO += R.qo * dt; R.cumG += R.qg * dt; R.cumW += R.qw * dt; R.cumGI += R.qgi * dt; R.cumWI += R.qwi * dt;
            R.cumMolP.resize(size_t(ncw_), 0.0); R.cumMolI.resize(size_t(ncw_), 0.0);
            for (int i = 0; i < ncw_; ++i) {
                const double v = X.molar[size_t(i)];
                if (v > 0) { cumProd_[size_t(i)] += v * dt; R.cumMolP[size_t(i)] += v * dt; }
                else { cumInj_[size_t(i)] -= v * dt; R.cumMolI[size_t(i)] -= v * dt; }
            }
        }
        t_ += dt; st_.steps++;
        if (forcedIdx_ < forced_.size()) ++forcedIdx_;
        if (adjOn_) adjRecord(dt);
        if (dbg_) std::fprintf(stderr, "step %d t %.4g d dt %.4g d it %d cuts %d\n", st_.steps, t_ / 86400, dt / 86400, it, st_.cuts);
        // next step: Newton-count factor combined with target-change control (Aziz-Settari, omega = 1) on the cell
        // pressure and overall-composition changes; growth limited after a chop
        {
            double dpm = 0, dzm = 0;
            for (int q = 0; q < n_; ++q) {
                dpm = std::max(dpm, std::fabs(p_[size_t(q)] - pOld_[size_t(q)]));
                double Nt = 0, No = 0;
                for (int i = 0; i < ncw_; ++i) { Nt += N_[size_t(q) * ncw_ + i]; No += NOld_[size_t(q) * ncw_ + i]; }
                if (Nt > 0 && No > 0)
                    for (int i = 0; i < ncw_; ++i)
                        dzm = std::max(dzm, std::fabs(N_[size_t(q) * ncw_ + i] / Nt - NOld_[size_t(q) * ncw_ + i] / No));
            }
            double f = it <= 4 ? 2.0 : it <= 8 ? 1.3 : 1.0;
            const double dpT = dpTarget_, dzT = dzTarget_;
            f = std::min(f, 2.0 * dpT / (dpm + dpT));
            f = std::min(f, 2.0 * dzT / (dzm + dzT));
            if (cutLast_) f = std::min(f, 1.25);
            f = std::max(f, 0.5);
            dtLast_ = dt * f;
            cutLast_ = false;
        }
        endStep(dt);
    }
    return true;
}

std::vector<double> CompSimulator::gasSaturation() const {
    std::vector<double> sg(static_cast<size_t>(n_), 0.0);
    for (int q = 0; q < n_; ++q) {
        double v[3]; phaseVolumes(c_, tp_, nc_, ncw_, p_[size_t(q)], &N_[size_t(q) * ncw_], K_[size_t(q)], v, &cf_);
        const double f = (c_.water || c_.co2store) ? 1.0 : 1.0 - c_.swc;
        sg[size_t(q)] = v[1] / (v[0] + v[1] + v[2]) * f;
    }
    return sg;
}

std::vector<double> CompSimulator::waterSaturation() const {
    std::vector<double> sw(static_cast<size_t>(n_), c_.swc);
    if (c_.co2store) {
        for (int q = 0; q < n_; ++q) {
            double v[3]; phaseVolumes(c_, tp_, nc_, ncw_, p_[size_t(q)], &N_[size_t(q) * ncw_], K_[size_t(q)], v, &cf_);
            sw[size_t(q)] = v[0] / std::max(1e-300, v[0] + v[1]);
        }
        return sw;
    }
    if (!c_.water) return sw;
    for (int q = 0; q < n_; ++q) {
        double v[3]; phaseVolumes(c_, tp_, nc_, ncw_, p_[size_t(q)], &N_[size_t(q) * ncw_], K_[size_t(q)], v, &cf_);
        sw[size_t(q)] = v[2] / (v[0] + v[1] + v[2]);
    }
    return sw;
}

double CompSimulator::vaporisedWater() const {
    if (!c_.waterVap || !c_.water) return 0.0;
    double s = 0;
    for (int q = 0; q < n_; ++q) {
        double v[3], nw = 0;
        phaseVolumes(c_, tp_, nc_, ncw_, p_[size_t(q)], &N_[size_t(q) * ncw_], K_[size_t(q)], v, &cf_, nullptr, &nw);
        s += nw;
    }
    return s;
}

double CompSimulator::componentInPlace(int i) const {
    double s = 0; for (int q = 0; q < n_; ++q) s += N_[size_t(q) * ncw_ + i];
    return s;
}

double CompSimulator::hcpvPressure() const {
    double sv = 0, svp = 0;
    for (int q = 0; q < n_; ++q) {
        double v[3]; phaseVolumes(c_, tp_, nc_, ncw_, p_[size_t(q)], &N_[size_t(q) * ncw_], K_[size_t(q)], v, &cf_);
        const double h = v[0] + v[1];
        sv += h; svp += h * p_[size_t(q)];
    }
    return sv > 0 ? svp / sv : 0.0;
}

void CompSimulator::co2Inventory(double& mobile, double& immobile, double& dissolved) const {
    const Co2Inventory I = co2Trapping();
    mobile = I.mobile + I.structural; immobile = I.residual; dissolved = I.dissolved;
}

// ------------------------------------------------------------------------------------------------- trapping (C5)
void CompSimulator::updateRockFactors() {
    const bool salt = c_.co2store && c_.trap.saltPrecip;
    const size_t nm = c_.co2store ? c_.trap.minerals.size() : 0;
    if (!salt && nm == 0) { if (!geoPv_.empty() || !geoK_.empty() || !connX_.empty()) composeFactors(); return; }
    for (int q = 0; q < n_; ++q) {
        double frac = 0;
        if (salt) frac += saltSolid_[size_t(q)] * co2::MW_NACL / c_.trap.rhoSalt / pv0_[size_t(q)];
        for (size_t r = 0; r < nm; ++r) frac += minProd_[size_t(q) * nm + r] * c_.trap.minerals[r].dVsolid / pv0_[size_t(q)];
        const double phiR = std::clamp(1.0 - frac, 0.01, 2.0);
        pvMult_[size_t(q)] = phiR;
        kMult_[size_t(q)] = phiR >= 1.0 ? 1.0 : std::max(1e-8, co2::permFactor(c_.trap, phiR));
        pvEff_[size_t(q)] = pv0_[size_t(q)] * phiR;
    }
    composeFactors();
}

// C25: chemical (salt / mineral) x external (geomechanics) factors -> effective pore volumes and transmissibilities
void CompSimulator::composeFactors() {
    for (int q = 0; q < n_; ++q) {
        const double g = geoPv_.empty() ? 1.0 : geoPv_[size_t(q)];
        pvEff_[size_t(q)] = pv0_[size_t(q)] * (pvMult_.empty() ? 1.0 : pvMult_[size_t(q)]) * g;
    }
    for (size_t e = 0; e < conns_.size() && e < T0_.size(); ++e) {
        const int a = conns_[e].a, b = conns_[e].b;
        const double ka = (kMult_.empty() ? 1.0 : kMult_[size_t(a)]) * (geoK_.empty() ? 1.0 : geoK_[size_t(a)]);
        const double kb = (kMult_.empty() ? 1.0 : kMult_[size_t(b)]) * (geoK_.empty() ? 1.0 : geoK_[size_t(b)]);
        conns_[e].T = T0_[e] * (ka + kb > 0 ? 2 * ka * kb / (ka + kb) : 0.0) * (connX_.empty() ? 1.0 : connX_[e]);
    }
}

void CompSimulator::setGeoFactors(const std::vector<double>& pvFactor, const std::vector<double>& permFactor) {
    auto fix = [&](const std::vector<double>& in, std::vector<double>& out) {
        if (int(in.size()) != n_) { out.clear(); return; }
        out = in; for (double& v : out) v = std::isfinite(v) ? std::clamp(v, 1e-4, 1e4) : 1.0;
    };
    fix(pvFactor, geoPv_); fix(permFactor, geoK_);
    composeFactors();
}

void CompSimulator::setConnectionMultipliersByCells(const std::vector<std::pair<int, int>>& cells, const std::vector<double>& mult) {
    std::map<std::pair<int, int>, size_t> idx;
    for (size_t e = 0; e < conns_.size(); ++e) idx[{std::min(conns_[e].a, conns_[e].b), std::max(conns_[e].a, conns_[e].b)}] = e;
    connX_.assign(conns_.size(), 1.0);
    for (size_t q = 0; q < cells.size() && q < mult.size(); ++q) {
        auto it = idx.find({std::min(cells[q].first, cells[q].second), std::max(cells[q].first, cells[q].second)});
        if (it != idx.end()) connX_[it->second] = std::isfinite(mult[q]) ? std::clamp(mult[q], 1e-8, 1e8) : 1.0;
    }
    composeFactors();
}

void CompSimulator::endStep(double dt) {
    if (c_.hyst.on) {   // scanning curves start from the historical maximum gas saturation
        const std::vector<double> sg = gasSaturation();
        for (int q = 0; q < n_; ++q) sgHist_[size_t(q)] = std::max(sgHist_[size_t(q)], sg[size_t(q)]);
    }
    const size_t nm = c_.co2store ? c_.trap.minerals.size() : 0;
    for (size_t i = 0; i < minRate_.size() && nm; ++i) {
        minProd_[i] += minRate_[i] * dt;
        minReact_[i] = std::max(0.0, minReact_[i] - minRate_[i] * dt);
    }
    if (!(c_.co2store && c_.trap.saltPrecip)) return;
    // ---- dissolved salt: implicit upwind transport with the converged brine (aqueous H2O) fluxes, then precipitation /
    // re-dissolution against the solubility limit (water vaporised into the CO2 leaves its salt behind)
    CompStats st; EvalCtx E{&c_, &tp_, nc_, &st, &cf_, &sgHist_};
    std::vector<CellEval> ce(static_cast<size_t>(n_));
    std::vector<double> W(static_cast<size_t>(n_), 0.0);
    for (int q = 0; q < n_; ++q) {
        std::vector<double> Kq; int np = 0;
        evalCell(E, q, p_[size_t(q)], &N_[size_t(q) * ncw_], Kq, np, ce[size_t(q)]);
        const Phase& A = ce[size_t(q)].ph[0];
        W[size_t(q)] = A.on ? A.vol.v / A.vm.v * A.x[1].v : 0.0;
    }
    std::vector<double> out(static_cast<size_t>(n_), 0.0), src(static_cast<size_t>(n_), 0.0);
    std::vector<std::vector<std::pair<int, double>>> inflow(static_cast<size_t>(n_));
    for (const Conn& cn : conns_) {
        const Phase& PA = ce[size_t(cn.a)].ph[0]; const Phase& PB = ce[size_t(cn.b)].ph[0];
        if (!PA.on && !PB.on) continue;
        const double wa = PA.on ? (PB.on ? 0.5 : 1.0) : 0.0, wb = PB.on ? (PA.on ? 0.5 : 1.0) : 0.0;
        const double rho = wa * (PA.on ? PA.rho.v : 0.0) + wb * (PB.on ? PB.rho.v : 0.0);
        const double dphi = (p_[size_t(cn.a)] + PA.pc.v) - (p_[size_t(cn.b)] + PB.pc.v) + (c_.gravity ? rho * GRAV * cn.dD : 0.0);
        const bool upA = dphi >= 0;
        const Phase& U = upA ? PA : PB;
        if (!U.on) continue;
        const double F = cn.T * U.mob.v * U.x[1].v * std::fabs(dphi);   // mol H2O / s in the brine
        const int u = upA ? cn.a : cn.b, d = upA ? cn.b : cn.a;
        out[size_t(u)] += F;
        inflow[size_t(d)].push_back({u, F});
    }
    for (size_t w = 0; w < c_.wells.size() && w < rep_.size() && w < aux_.size(); ++w) {
        const CompWell& Wl = c_.wells[w];
        if (!Wl.open || Wl.cells.empty()) continue;
        const double bhp = rep_[w].bhp;
        if (!Wl.injector) {
            for (size_t k = 0; k < Wl.cells.size(); ++k) {
                const int q = Wl.cells[k];
                const Phase& P = ce[size_t(q)].ph[0];
                if (!P.on) continue;
                const double wi = (k < Wl.wi.size() ? Wl.wi[k] : 1.0) * kMult_[size_t(q)];
                const double dp = p_[size_t(q)] + P.pc.v - (k < aux_[w].H.size() ? aux_[w].H[k] : 0.0) - bhp;
                if (dp > 0) out[size_t(q)] += wi * P.mob.v * P.x[1].v * dp;
            }
        } else if (!Wl.zinj.empty() && Wl.zinj[0] < 0.5 && rep_[w].molar.size() > 1) {   // brine injector: salt at the model salinity
            const double h2o = std::max(0.0, -rep_[w].molar[1]);
            double sw = 0; std::vector<double> lt(Wl.cells.size(), 0.0);
            for (size_t k = 0; k < Wl.cells.size(); ++k) {
                const CellEval& C = ce[size_t(Wl.cells[k])];
                for (int al = 0; al < 2; ++al) if (C.ph[al].on) lt[k] += C.ph[al].vmob.v;
                lt[k] *= k < Wl.wi.size() ? Wl.wi[k] : 1.0; sw += lt[k];
            }
            for (size_t k = 0; k < Wl.cells.size(); ++k)
                src[size_t(Wl.cells[k])] += (sw > 0 ? lt[k] / sw : 1.0 / double(Wl.cells.size())) * h2o * c_.salinity * MWW;
        }
    }
    const std::vector<double> m0 = saltDiss_;
    std::vector<double> m = m0;
    const double tiny = 1e-30;
    for (int sweep = 0; sweep < 4000; ++sweep) {
        double ch = 0, sc = 0;
        for (int dir = 0; dir < 2; ++dir)
            for (int ii = 0; ii < n_; ++ii) {
                const int q = dir == 0 ? ii : n_ - 1 - ii;
                double rhs = m0[size_t(q)] + dt * src[size_t(q)];
                for (const auto& f : inflow[size_t(q)]) rhs += dt * f.second * m[size_t(f.first)] / std::max(W[size_t(f.first)], tiny);
                const double diag = 1.0 + (W[size_t(q)] > tiny ? dt * out[size_t(q)] / W[size_t(q)] : 0.0);
                const double v = rhs / diag;
                ch = std::max(ch, std::fabs(v - m[size_t(q)])); sc = std::max(sc, std::fabs(v));
                m[size_t(q)] = v;
            }
        if (ch <= 1e-12 * std::max(sc, 1e-30)) break;
    }
    for (int q = 0; q < n_; ++q) {
        const double tot = std::max(0.0, m[size_t(q)]) + saltSolid_[size_t(q)];
        const double cap = c_.trap.saltSolubility * W[size_t(q)] * MWW;
        saltDiss_[size_t(q)] = std::min(tot, cap);
        saltSolid_[size_t(q)] = tot - saltDiss_[size_t(q)];
    }
}

Co2Inventory CompSimulator::co2Trapping() const {
    Co2Inventory I;
    if (!c_.co2store) return I;
    CompStats st; EvalCtx E{&c_, &tp_, nc_, &st, &cf_, &sgHist_};
    std::vector<CellEval> ce(static_cast<size_t>(n_));
    for (int q = 0; q < n_; ++q) { std::vector<double> Kq; int np = 0; evalCell(E, q, p_[size_t(q)], &N_[size_t(q) * ncw_], Kq, np, ce[size_t(q)]); }
    // gas-phase outflow per cell (reservoir m3/s): mobile CO2 that does not leave its cell is held by a seal
    std::vector<double> gout(static_cast<size_t>(n_), 0.0);
    for (const Conn& cn : conns_) {
        const Phase& PA = ce[size_t(cn.a)].ph[1]; const Phase& PB = ce[size_t(cn.b)].ph[1];
        if (!PA.on && !PB.on) continue;
        const double wa = PA.on ? (PB.on ? 0.5 : 1.0) : 0.0, wb = PB.on ? (PA.on ? 0.5 : 1.0) : 0.0;
        const double rho = wa * (PA.on ? PA.rho.v : 0.0) + wb * (PB.on ? PB.rho.v : 0.0);
        const double dphi = (p_[size_t(cn.a)] + PA.pc.v) - (p_[size_t(cn.b)] + PB.pc.v) + (c_.gravity ? rho * GRAV * cn.dD : 0.0);
        const bool upA = dphi >= 0;
        const Phase& U = upA ? PA : PB;
        if (!U.on) continue;
        gout[size_t(upA ? cn.a : cn.b)] += cn.T * U.vmob.v * std::fabs(dphi);
    }
    const double year = 365.25 * 86400.0;
    for (int q = 0; q < n_; ++q) {
        const CellEval& C = ce[size_t(q)];
        if (C.ph[0].on) I.dissolved += C.ph[0].vol.v / C.ph[0].vm.v * C.ph[0].x[0].v * co2::MW_CO2;
        if (!C.ph[1].on) continue;
        const Phase& G = C.ph[1];
        const double mG = G.vol.v / G.vm.v * G.x[0].v * co2::MW_CO2;
        const double Sg = G.vol.v / C.vtot.v;
        const SatFn sf = satFn(c_, q, &sgHist_);
        const double sgt = sf.trapped(std::max(sgHist_[size_t(q)], Sg));
        const double fr = Sg > 0 ? std::min(1.0, sgt / Sg) : 1.0;
        I.residual += fr * mG;
        if (sf.krg(DV(Sg)).v <= 0) I.immobileKr += mG;
        const double mm = (1 - fr) * mG;
        if (mm <= 0) continue;
        const double vMob = (1 - fr) * G.vol.v;
        if (gout[size_t(q)] * year < c_.trap.structTol * vMob) I.structural += mm; else I.mobile += mm;
    }
    const size_t nm = c_.trap.minerals.size();
    for (int q = 0; q < n_; ++q)
        for (size_t r = 0; r < nm; ++r) I.mineral += minProd_[size_t(q) * nm + r] * c_.trap.minerals[r].nuCO2 * co2::MW_CO2;
    for (double v : saltSolid_) I.saltSolid += v * co2::MW_NACL;
    return I;
}

std::vector<double> CompSimulator::trappedGasSaturation() const {
    std::vector<double> t(static_cast<size_t>(n_), 0.0);
    const std::vector<double> sg = gasSaturation();
    for (int q = 0; q < n_; ++q) t[size_t(q)] = satFn(c_, q, &sgHist_).trapped(std::max(sgHist_[size_t(q)], sg[size_t(q)]));
    return t;
}

std::vector<double> CompSimulator::saltSaturation() const {
    std::vector<double> v(static_cast<size_t>(n_), 0.0);
    for (int q = 0; q < n_ && size_t(q) < saltSolid_.size(); ++q) v[size_t(q)] = saltSolid_[size_t(q)] * co2::MW_NACL / c_.trap.rhoSalt / pv0_[size_t(q)];
    return v;
}

void CompSimulator::saltInPlace(double& dissolved, double& solid) const {
    dissolved = 0; solid = 0;
    for (double v : saltDiss_) dissolved += v;
    for (double v : saltSolid_) solid += v;
}

std::vector<double> CompSimulator::mineralCO2() const {
    std::vector<double> v(static_cast<size_t>(n_), 0.0);
    const size_t nm = c_.co2store ? c_.trap.minerals.size() : 0;
    for (int q = 0; q < n_; ++q)
        for (size_t r = 0; r < nm; ++r) v[size_t(q)] += minProd_[size_t(q) * nm + r] * c_.trap.minerals[r].nuCO2 * co2::MW_CO2 / pv0_[size_t(q)];
    return v;
}

void CompSimulator::compositions(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z) const {
    x.assign(size_t(n_) * nc_, 0.0); y.assign(size_t(n_) * nc_, 0.0); z.assign(size_t(n_) * nc_, 0.0);
    for (int q = 0; q < n_; ++q) {
        const double* N = &N_[size_t(q) * ncw_];
        double Nt = 0; for (int i = 0; i < nc_; ++i) Nt += std::max(N[i], 0.0);
        std::vector<double> zz(static_cast<size_t>(nc_));
        for (int i = 0; i < nc_; ++i) { zz[size_t(i)] = std::max(N[i], 0.0) / std::max(Nt, 1e-300); z[size_t(q) * nc_ + i] = zz[size_t(i)]; }
        if (c_.co2store) {
            const co2::FlashResult f = co2::flash(cf_, p_[size_t(q)], zz[0]);
            const bool aq = f.phases == 2 || !f.gas, gs = f.phases == 2 || f.gas;
            if (aq) { const double xa = f.phases == 2 ? f.x : zz[0]; x[size_t(q) * nc_] = xa; x[size_t(q) * nc_ + 1] = 1 - xa; }
            if (gs) { const double yg = f.phases == 2 ? f.y : zz[0]; y[size_t(q) * nc_] = yg; y[size_t(q) * nc_ + 1] = 1 - yg; }
            continue;
        }
        eos::FlashOptions fo; if (!K_[size_t(q)].empty()) fo.Kinit = &K_[size_t(q)];
        const eos::FlashResult f = eos::flash(c_.fluid, tp_, p_[size_t(q)], zz, fo);
        for (int i = 0; i < nc_; ++i) {
            if (f.phases == 2) { x[size_t(q) * nc_ + i] = f.x[size_t(i)]; y[size_t(q) * nc_ + i] = f.y[size_t(i)]; }
            else if (f.vapourLabel) y[size_t(q) * nc_ + i] = zz[size_t(i)];
            else x[size_t(q) * nc_ + i] = zz[size_t(i)];
        }
    }
}

double CompSimulator::aqueousCo2Moles(int q) const {
    if (!c_.co2sol || !c_.water || c_.co2Comp < 0) return 0.0;
    const double Nc = N_[size_t(q) * ncw_ + c_.co2Comp];
    if (size_t(q) < sAq_.size() && sAq_[size_t(q)] >= 0) return sAq_[size_t(q)] * std::max(Nc, 0.0);   // last assembly (converged state)
    return solveAqueousCo2(c_, tp_, nc_, p_[size_t(q)], &N_[size_t(q) * ncw_], K_[size_t(q)], -1.0);
}

void CompSimulator::co2SolSplit(double& hcKg, double& aqKg) const {
    hcKg = aqKg = 0;
    if (!c_.co2sol || !c_.water || c_.co2Comp < 0) return;
    for (int q = 0; q < n_; ++q) {
        const double s = aqueousCo2Moles(q), Nc = std::max(0.0, N_[size_t(q) * ncw_ + c_.co2Comp]);
        aqKg += s * MWCO2; hcKg += std::max(0.0, Nc - s) * MWCO2;
    }
}

std::vector<double> CompSimulator::dissolvedCO2() const {
    std::vector<double> x(static_cast<size_t>(n_), 0.0);
    if (c_.co2sol && c_.water && c_.co2Comp >= 0) {   // CO2SOL: mole fraction of CO2 in the aqueous phase
        for (int q = 0; q < n_; ++q) {
            const double s = aqueousCo2Moles(q), Nw = std::max(0.0, N_[size_t(q) * ncw_ + nc_]);
            x[size_t(q)] = Nw + s > 0 ? s / (Nw + s) : 0.0;
        }
        return x;
    }
    if (!c_.co2store) return x;
    for (int q = 0; q < n_; ++q) {
        double v[3], xa = 0; phaseVolumes(c_, tp_, nc_, ncw_, p_[size_t(q)], &N_[size_t(q) * ncw_], K_[size_t(q)], v, &cf_, &xa);
        x[size_t(q)] = v[0] > 0 ? xa : 0.0;
    }
    return x;
}

void CompSimulator::inPlaceSurface(double& oil, double& gas, double& wat) const {
    oil = gas = wat = 0;
    if (c_.co2store) {
        for (int q = 0; q < n_; ++q) { gas += N_[size_t(q) * ncw_] * co2Vg_; wat += N_[size_t(q) * ncw_ + 1] * co2Vw_; }
        return;
    }
    std::vector<double> z(static_cast<size_t>(nc_));
    for (int q = 0; q < n_; ++q) {
        double Nt = 0;
        for (int i = 0; i < nc_; ++i) { z[size_t(i)] = N_[size_t(q) * ncw_ + i]; Nt += z[size_t(i)]; }
        double vo, vg; surfaceVolumes(z, vo, vg);
        oil += vo * Nt; gas += vg * Nt;
        if (c_.water) wat += N_[size_t(q) * ncw_ + nc_] * waterSurfaceVolumePerMole();
    }
}


// ---------------------------------------------------------------- C14 restart state
template <class Ar> void CompSimulator::stateIO(Ar& ar) {
    int n = n_, ncw = ncw_; uint64_t nw = c_.wells.size(), ncn = conns_.size();
    const int n0 = n, ncw0 = ncw; const uint64_t nw0 = nw, ncn0 = ncn;
    ar(n, ncw, nw, ncn);
    if (Ar::reading && (n != n0 || ncw != ncw0 || nw != nw0 || ncn != ncn0))
        throw std::runtime_error("compositional restart state does not match the model (cells " + std::to_string(n) + " / " + std::to_string(n0) +
                                 ", components " + std::to_string(ncw) + " / " + std::to_string(ncw0) + ", wells " + std::to_string(nw) + " / " + std::to_string(nw0) + ")");
    ar(t_, dtLast_, cutLast_, p_, N_, pOld_, NOld_, K_, nph_, cumInj_, cumProd_, st_, sAq_);
    ar(sgHist_, saltDiss_, saltSolid_, minProd_, minReact_, minRate_, pvMult_, kMult_, pvEff_);
    for (Conn& e : conns_) ar(e.T, e.Td);
    if (Ar::reading) { aux_.resize(size_t(nw)); rep_.resize(size_t(nw)); }
    for (WellAux& a : aux_) ar(a.H, a.fo, a.fg, a.fw);
    for (CompWellReport& r : rep_) {
        ar(r.bhp, r.onBhp, r.control, r.qo, r.qg, r.qw, r.qgi, r.qwi, r.resv, r.molar);
        ar(r.cumO, r.cumG, r.cumW, r.cumGI, r.cumWI, r.cumMolP, r.cumMolI);
    }
}

std::string CompSimulator::saveState() const {
    StateOut o;
    o.w.str("RESSIM-COMP-STATE-1");
    const_cast<CompSimulator*>(this)->stateIO(o);
    return o.w.buf;
}

void CompSimulator::loadState(const std::string& blob) {
    StateIn in(blob);
    if (in.r.str() != "RESSIM-COMP-STATE-1") throw std::runtime_error("not a compositional restart state");
    stateIO(in);
}


// ================================================================================================ C18 discrete adjoint
// Forward model per accepted step n (exactly the residual of advanceTo / assemble):
//     R_n = [N_n - N_{n-1}] + dt_n F(p_n, N_n; theta)  (component rows),   PV(p_n) - V_fluid(p_n, N_n)  (volume row)
// with the rate-controlled BHPs eliminated (exact Schur complement), so x_n = (p, N) per cell and
//     dR_n/dx_{n-1} = -I on the component rows (moles are the primary unknowns) - no other coupling.
// J_n is the analytic (AD) Jacobian of the Newton solver re-evaluated at the recorded converged state. dR_n/dtheta and
// dy/dtheta|_x are central differences of two perturbed models evaluated at the SAME recorded states (no re-simulation);
// the initial-state dependence (pore-volume multipliers change the initial moles) enters through dR_1/dx_0.
// Observations y are well rates / BHP from the well model at x_n: dy/dx_n by central differences over the perforated
// cells (the only cells a well report depends on). Lagged well data (wellbore heads, separator yields of rate targets),
// group allocation, hysteresis history and rock factors are frozen at their recorded values (frozen-coefficient
// discrete adjoint, as for the black-oil engine).
namespace {
void mergeWells(std::vector<CompWell>& dst, const std::vector<CompWell>& src) {
    bool same = dst.size() == src.size();
    for (size_t w = 0; same && w < src.size(); ++w) same = dst[w].cells == src[w].cells && dst[w].wi.size() == src[w].wi.size();
    if (!same) { dst = src; return; }
    for (size_t w = 0; w < src.size(); ++w) { std::vector<double> wi = dst[w].wi; dst[w] = src[w]; dst[w].wi = wi; }
}
std::vector<int> regionOrAll(const std::vector<int>& cells, int n) {
    if (!cells.empty()) return cells;
    std::vector<int> a(static_cast<size_t>(n)); for (int i = 0; i < n; ++i) a[size_t(i)] = i; return a;
}
void expandTo(std::vector<double>& v, int n) { if (v.size() == 1 && n > 1) v.assign(size_t(n), v[0]); }
}  // namespace

CompAdjParam compParamPermMult(const std::string& name, const std::vector<int>& cells) {
    CompAdjParam P; P.name = name;
    P.apply = [cells](CompCase& c, double m) {
        const int n = c.size();
        const std::vector<int> reg = regionOrAll(cells, n);
        if (c.cellPV.empty()) {
            expandTo(c.permx, n); expandTo(c.permy, n); expandTo(c.permz, n);
            for (int q : reg) if (q >= 0 && q < n) { c.permx[size_t(q)] *= m; c.permy[size_t(q)] *= m; c.permz[size_t(q)] *= m; }
        } else {
            std::vector<char> in(size_t(n), 0); for (int q : reg) if (q >= 0 && q < n) in[size_t(q)] = 1;
            for (auto& cn : c.conns) if ((cn.a >= 0 && cn.a < n && in[size_t(cn.a)]) || (cn.b >= 0 && cn.b < n && in[size_t(cn.b)])) cn.T *= m;
        }
    };
    return P;
}
CompAdjParam compParamPoroMult(const std::string& name, const std::vector<int>& cells) {
    CompAdjParam P; P.name = name;
    P.apply = [cells](CompCase& c, double m) {
        const int n = c.size();
        const std::vector<int> reg = regionOrAll(cells, n);
        if (c.cellPV.empty()) { expandTo(c.poro, n); for (int q : reg) if (q >= 0 && q < n) c.poro[size_t(q)] *= m; }
        else for (int q : reg) if (q >= 0 && q < n) c.cellPV[size_t(q)] *= m;
    };
    return P;
}
CompAdjParam compParamWellIndexMult(const std::string& name, int well) {
    CompAdjParam P; P.name = name;
    P.apply = [well](CompCase& c, double m) { if (well >= 0 && size_t(well) < c.wells.size()) for (double& v : c.wells[size_t(well)].wi) v *= m; };
    return P;
}

void CompSimulator::setAdjointRecording(bool on) {
    adjOn_ = on; adjRec_.clear();
    if (on) { if (aux_.size() != c_.wells.size()) beginStep(); adjRecord(0.0); }
}
std::vector<double> CompSimulator::stepSizes() const {
    std::vector<double> d; for (size_t k = 1; k < adjRec_.size(); ++k) d.push_back(adjRec_[k].dt); return d;
}
void CompSimulator::adjRecord(double dt) {
    AdjRec r; r.t = t_; r.dt = dt; r.p = p_; r.N = N_; r.K = K_; r.nph = nph_; r.aux = aux_; r.wells = c_.wells;
    r.pvEff = pvEff_; r.kMult = kMult_;
    if (c_.hyst.on) r.sgHist = sgHist_;
    bool changed = false;
    for (size_t e = 0; e < conns_.size() && e < T0_.size(); ++e) if (conns_[e].T != T0_[e]) { changed = true; break; }
    if (changed) { r.connT.resize(conns_.size()); for (size_t e = 0; e < conns_.size(); ++e) r.connT[e] = conns_[e].T; }
    adjRec_.push_back(std::move(r));
}
void CompSimulator::adjEval(const AdjRec& prev, const AdjRec& cur, const CompSimulator& base, std::vector<double>& R, BSR* J,
                            std::vector<CompWellReport>* wr) {
    mergeWells(c_.wells, cur.wells);
    // lagged well data (wellbore heads, separator yields) re-derived from the start-of-step state of THIS model, so
    // their dependence on x_{n-1} and theta is part of the discrete adjoint (identical to the recorded data for the base)
    p_ = prev.p; N_ = prev.N; K_ = prev.K; kMult_ = cur.kMult;
    if (c_.hyst.on && !cur.sgHist.empty()) sgHist_ = cur.sgHist;
    beginStep();
    p_ = cur.p; N_ = cur.N; pOld_ = prev.p; NOld_ = prev.N; K_ = cur.K; nph_ = cur.nph;
    pvEff_.resize(size_t(n_));
    for (int q = 0; q < n_; ++q) {
        const double b0 = base.pv0_[size_t(q)];
        pvEff_[size_t(q)] = pv0_[size_t(q)] * (b0 > 0 && size_t(q) < cur.pvEff.size() ? cur.pvEff[size_t(q)] / b0 : 1.0);
    }
    kMult_ = cur.kMult;
    if (c_.hyst.on && !cur.sgHist.empty()) sgHist_ = cur.sgHist;
    for (size_t e = 0; e < conns_.size(); ++e)
        conns_[e].T = !cur.connT.empty() && e < base.T0_.size() && base.T0_[e] > 0 ? T0_[e] * cur.connT[e] / base.T0_[e] : T0_[e];
    Sys s; buildPattern(n_, nv_, conns_, c_.wells, s);
    AsmIn in{c_, tp_, nc_, pvEff_, conns_, aux_, cf_, &sgHist_, &kMult_, &pv0_, &minReact_, &minProd_, &minRate_, &sAq_};
    CompStats tmp; std::vector<CompWellReport> w2;
    assemble(in, cur.dt, p_, N_, NOld_, K_, nph_, tmp, s, wr ? wr : &w2);
    R = s.R;
    if (J) *J = s.J;
}
double CompSimulator::adjObsValue(const CompAdjObs& o, const CompWellReport& X) const {
    if (o.qty == 3) return X.bhp;
    if (X.molar.empty()) return 0.0;
    if (o.qty >= 10) { const size_t i = size_t(o.qty - 10); return i < X.molar.size() ? X.molar[i] : 0.0; }
    double qo = 0, qg = 0, qw = 0, qgi = 0, qwi = 0;
    const double vgsc = gasSurfaceVolumePerMole(), vwsc = waterSurfaceVolumePerMole();
    std::vector<double> zp(size_t(nc_), 0.0); double np = 0, ni = 0;
    for (int i = 0; i < nc_; ++i) { const double v = X.molar[size_t(i)]; if (v > 0) { zp[size_t(i)] = v; np += v; } else ni -= v; }
    if (c_.co2store) {
        const double c0 = X.molar[0], w0 = X.molar[1];
        if (c0 > 0) qg = c0 * co2Vg_; else qgi = -c0 * co2Vg_;
        if (w0 > 0) qw = w0 * co2Vw_; else qwi = -w0 * co2Vw_;
    } else {
        if (np > 0) { double vo, vg; surfaceVolumes(zp, vo, vg, size_t(o.well) < c_.wells.size() ? c_.wells[size_t(o.well)].sep : 0); qo = vo * np; qg = vg * np; }
        qgi = ni * vgsc;
        if (c_.water && X.molar.size() > size_t(nc_)) { const double v = X.molar[size_t(nc_)]; if (v > 0) qw = v * vwsc; else qwi = -v * vwsc; }
    }
    switch (o.qty) { case 0: return qo; case 1: return qg; case 2: return qw; case 4: return qgi; case 5: return qwi; default: return 0.0; }
}

namespace {
// observation -> (step a, weight wa, step b, weight wb), linear in time between recorded step ends
struct ObsMap { int ka = 1, kb = 1; double wa = 1, wb = 0; };
template <class Rec> ObsMap mapObs(const std::vector<Rec>& rec, double t) {
    ObsMap m; const int N = int(rec.size()) - 1;
    if (N < 1) return m;
    if (t <= rec[1].t) { m.ka = m.kb = 1; return m; }
    if (t >= rec[size_t(N)].t) { m.ka = m.kb = N; return m; }
    int k = 1; while (k < N && rec[size_t(k + 1)].t < t) ++k;
    const double a = (t - rec[size_t(k)].t) / std::max(1e-300, rec[size_t(k + 1)].t - rec[size_t(k)].t);
    m.ka = k; m.kb = k + 1; m.wa = 1 - a; m.wb = a;
    return m;
}
}  // namespace

double CompSimulator::adjointMisfit(const std::vector<CompAdjObs>& obs) {
    if (adjRec_.size() < 2) return 0.0;
    const std::string saved = saveState();
    const std::vector<CompWell> wSaved = c_.wells; const std::vector<double> pvS = pvEff_; const std::vector<Conn> cS = conns_;
    std::map<int, std::vector<CompWellReport>> rep;
    std::vector<double> R;
    double F = 0;
    for (const auto& o : obs) {
        const ObsMap m = mapObs(adjRec_, o.t);
        double y = 0;
        for (int side = 0; side < 2; ++side) {
            const int k = side ? m.kb : m.ka; const double w = side ? m.wb : m.wa;
            if (w == 0) continue;
            if (!rep.count(k)) adjEval(adjRec_[size_t(k - 1)], adjRec_[size_t(k)], *this, R, nullptr, &rep[k]);
            y += w * (size_t(o.well) < rep[k].size() ? adjObsValue(o, rep[k][size_t(o.well)]) : 0.0);
        }
        const double r = (y - o.value) / o.sigma; F += r * r;
    }
    c_.wells = wSaved; pvEff_ = pvS; conns_ = cS; loadState(saved);
    return F;
}

CompAdjResult CompSimulator::adjointGradient(const std::vector<CompAdjObs>& obs, const std::vector<CompAdjParam>& params) {
    CompAdjResult res;
    const auto t0 = std::chrono::steady_clock::now();
    const int N = int(adjRec_.size()) - 1, np = int(params.size()), m = n_ * nv_;
    if (N < 1) { res.error = "no recorded steps (setAdjointRecording before the run)"; return res; }
    const std::string saved = saveState();
    const std::vector<CompWell> wSaved = c_.wells; const std::vector<double> pvS = pvEff_; const std::vector<Conn> cS = conns_;
    try {
        // observations: values, residuals, dF/dy
        std::vector<ObsMap> om(obs.size());
        std::map<int, std::vector<CompWellReport>> rep;
        std::vector<double> R;
        std::vector<double> dFdy(obs.size(), 0.0);
        for (size_t j = 0; j < obs.size(); ++j) {
            om[j] = mapObs(adjRec_, obs[j].t);
            double y = 0;
            for (int side = 0; side < 2; ++side) {
                const int k = side ? om[j].kb : om[j].ka; const double w = side ? om[j].wb : om[j].wa;
                if (w == 0) continue;
                if (!rep.count(k)) adjEval(adjRec_[size_t(k - 1)], adjRec_[size_t(k)], *this, R, nullptr, &rep[k]);
                y += w * (size_t(obs[j].well) < rep[k].size() ? adjObsValue(obs[j], rep[k][size_t(obs[j].well)]) : 0.0);
            }
            const double r = (y - obs[j].value) / obs[j].sigma;
            res.F += r * r; dFdy[j] = 2 * r / obs[j].sigma;
        }
        // per step: the observations it carries (index, weight)
        std::vector<std::vector<std::pair<int, double>>> atStep(size_t(N) + 1);
        for (size_t j = 0; j < obs.size(); ++j) {
            if (om[j].wa != 0) atStep[size_t(om[j].ka)].push_back({int(j), om[j].wa});
            if (om[j].wb != 0) atStep[size_t(om[j].kb)].push_back({int(j), om[j].wb});
        }
        // perturbed models (setup of the modified case: transmissibilities, pore volumes, initial moles, well indices)
        CompCase base = c_; base.wells = adjRec_[0].wells;
        std::vector<std::unique_ptr<CompSimulator>> sp(static_cast<size_t>(np)), sm(static_cast<size_t>(np));
        std::vector<std::vector<double>> dN0(static_cast<size_t>(np)), dP0(static_cast<size_t>(np));
        for (int j = 0; j < np; ++j) {
            const double h = params[size_t(j)].h;
            CompCase cp = base, cm = base;
            params[size_t(j)].apply(cp, 1 + h); params[size_t(j)].apply(cm, 1 - h);
            sp[size_t(j)].reset(new CompSimulator()); sm[size_t(j)].reset(new CompSimulator());
            sp[size_t(j)]->setup(cp); sm[size_t(j)]->setup(cm);
            dN0[size_t(j)].resize(N_.size());
            for (size_t i = 0; i < N_.size(); ++i) dN0[size_t(j)][i] = (sp[size_t(j)]->N_[i] - sm[size_t(j)]->N_[i]) / (2 * h);
            dP0[size_t(j)].resize(p_.size());
            for (size_t i = 0; i < p_.size(); ++i) dP0[size_t(j)][i] = (sp[size_t(j)]->p_[i] - sm[size_t(j)]->p_[i]) / (2 * h);
        }
        res.grad.assign(size_t(np), 0.0);
        // explicit dependence of the observations on theta at fixed states
        for (int k = 1; k <= N; ++k) {
            if (atStep[size_t(k)].empty()) continue;
            for (int j = 0; j < np; ++j) {
                std::vector<CompWellReport> rp, rm;
                sp[size_t(j)]->adjEval(adjRec_[size_t(k - 1)], adjRec_[size_t(k)], *this, R, nullptr, &rp);
                sm[size_t(j)]->adjEval(adjRec_[size_t(k - 1)], adjRec_[size_t(k)], *this, R, nullptr, &rm);
                for (const auto& ow : atStep[size_t(k)]) {
                    const CompAdjObs& o = obs[size_t(ow.first)];
                    if (size_t(o.well) >= rp.size()) continue;
                    const double dy = (sp[size_t(j)]->adjObsValue(o, rp[size_t(o.well)]) - sm[size_t(j)]->adjObsValue(o, rm[size_t(o.well)])) / (2 * params[size_t(j)].h);
                    res.grad[size_t(j)] += dFdy[size_t(ow.first)] * ow.second * dy;
                }
            }
        }
        StepAdjointProblem P;
        P.steps = N;
        P.jacobian = [&](int n, BSR& J) { adjEval(adjRec_[size_t(n - 1)], adjRec_[size_t(n)], *this, R, &J, nullptr); };
        // (dR_n/dx_{n-1})^T lam: -I on the component rows; perforated cells additionally through the lagged well data
        // (heads with gravity, separator yields of surface-rate controls) -> central differences of R_n in x_{n-1}
        bool auxActive = c_.gravity;
        for (const auto& rc : adjRec_) for (const CompWell& W : rc.wells) {
            bool y = W.control >= 3; for (const auto& l : W.limits) y = y || l.first >= 3;
            if (y && !W.injector) auxActive = true;
        }
        // observations of step nNext depend on x_{nNext-1} through the lagged well data (rate-control yields, heads):
        // out += sum dF/dy (dy_nNext / dx_{nNext-1}), central differences over the perforated cells
        auto objPrev = [&](int nNext, std::vector<double>& out) {
            if (!auxActive || nNext < 1 || nNext > N) return;
            const auto& list = atStep[size_t(nNext)];
            if (list.empty()) return;
            std::vector<char> seen(static_cast<size_t>(n_), 0); std::vector<int> cells;
            for (const CompWell& W : adjRec_[size_t(nNext)].wells) if (W.open) for (int q : W.cells) if (q >= 0 && q < n_ && !seen[size_t(q)]) { seen[size_t(q)] = 1; cells.push_back(q); }
            AdjRec pp = adjRec_[size_t(nNext - 1)];
            std::vector<CompWellReport> rp, rm;
            for (int q : cells) {
                double Nt = 0; for (int i = 0; i < ncw_; ++i) Nt += pp.N[size_t(q) * ncw_ + i];
                for (int v = 0; v < nv_; ++v) {
                    double& x = v == 0 ? pp.p[size_t(q)] : pp.N[size_t(q) * ncw_ + size_t(v - 1)];
                    const double x0 = x, d = v == 0 ? 1e-6 * std::max(std::fabs(x0), 1e5) : 1e-6 * std::max(Nt, 1e-30);
                    // only the lagged data may depend on x_{n-1} here: the accumulation term does not enter the well reports
                    x = x0 + d; adjEval(pp, adjRec_[size_t(nNext)], *this, R, nullptr, &rp);
                    x = x0 - d; adjEval(pp, adjRec_[size_t(nNext)], *this, R, nullptr, &rm);
                    x = x0;
                    double g = 0;
                    for (const auto& ow : list) {
                        const CompAdjObs& o = obs[size_t(ow.first)];
                        if (size_t(o.well) >= rp.size()) continue;
                        g += dFdy[size_t(ow.first)] * ow.second * (adjObsValue(o, rp[size_t(o.well)]) - adjObsValue(o, rm[size_t(o.well)])) / (2 * d);
                    }
                    out[size_t(q) * nv_ + size_t(v)] += g;
                }
            }
        };
        P.addObjectiveGradient = [&](int n, std::vector<double>& rhs) {
            if (n < N) objPrev(n + 1, rhs);
            const auto& list = atStep[size_t(n)];
            if (list.empty()) return;
            std::vector<char> seen(static_cast<size_t>(n_), 0); std::vector<int> cells;
            for (const auto& ow : list) {
                const int w = obs[size_t(ow.first)].well;
                if (w < 0 || size_t(w) >= adjRec_[size_t(n)].wells.size()) continue;
                for (int q : adjRec_[size_t(n)].wells[size_t(w)].cells) if (q >= 0 && q < n_ && !seen[size_t(q)]) { seen[size_t(q)] = 1; cells.push_back(q); }
            }
            AdjRec tmp = adjRec_[size_t(n)];
            std::vector<CompWellReport> rp, rm;
            for (int q : cells) {
                double Nt = 0; for (int i = 0; i < ncw_; ++i) Nt += tmp.N[size_t(q) * ncw_ + i];
                for (int v = 0; v < nv_; ++v) {
                    double& x = v == 0 ? tmp.p[size_t(q)] : tmp.N[size_t(q) * ncw_ + size_t(v - 1)];
                    const double x0 = x, d = v == 0 ? 1e-6 * std::max(std::fabs(x0), 1e5) : 1e-6 * std::max(Nt, 1e-30);
                    x = x0 + d; adjEval(adjRec_[size_t(n - 1)], tmp, *this, R, nullptr, &rp);
                    x = x0 - d; adjEval(adjRec_[size_t(n - 1)], tmp, *this, R, nullptr, &rm);
                    x = x0;
                    double g = 0;
                    for (const auto& ow : list) {
                        const CompAdjObs& o = obs[size_t(ow.first)];
                        if (size_t(o.well) >= rp.size()) continue;
                        g += dFdy[size_t(ow.first)] * ow.second * (adjObsValue(o, rp[size_t(o.well)]) - adjObsValue(o, rm[size_t(o.well)])) / (2 * d);
                    }
                    rhs[size_t(q) * nv_ + size_t(v)] += g;
                }
            }
        };
        auto prevCouplingT = [&](int n, const std::vector<double>& lam, std::vector<double>& rhs) {
            std::vector<char> perf(static_cast<size_t>(n_), 0);
            if (auxActive) for (const CompWell& W : adjRec_[size_t(n)].wells) if (W.open) for (int q : W.cells) if (q >= 0 && q < n_) perf[size_t(q)] = 1;
            for (int q = 0; q < n_; ++q) if (!perf[size_t(q)]) for (int i = 0; i < ncw_; ++i) rhs[size_t(q) * nv_ + 1 + i] -= lam[size_t(q) * nv_ + i];
            if (!auxActive) return;
            AdjRec pp = adjRec_[size_t(n - 1)];
            std::vector<double> Rp, Rm;
            for (int q = 0; q < n_; ++q) {
                if (!perf[size_t(q)]) continue;
                double Nt = 0; for (int i = 0; i < ncw_; ++i) Nt += pp.N[size_t(q) * ncw_ + i];
                for (int v = 0; v < nv_; ++v) {
                    double& x = v == 0 ? pp.p[size_t(q)] : pp.N[size_t(q) * ncw_ + size_t(v - 1)];
                    const double x0 = x, d = v == 0 ? 1e-6 * std::max(std::fabs(x0), 1e5) : 1e-6 * std::max(Nt, 1e-30);
                    x = x0 + d; adjEval(pp, adjRec_[size_t(n)], *this, Rp, nullptr, nullptr);
                    x = x0 - d; adjEval(pp, adjRec_[size_t(n)], *this, Rm, nullptr, nullptr);
                    x = x0;
                    double sdot = 0; for (int r = 0; r < m; ++r) sdot += lam[size_t(r)] * (Rp[size_t(r)] - Rm[size_t(r)]);
                    rhs[size_t(q) * nv_ + size_t(v)] += sdot / (2 * d);
                }
            }
        };
        P.addPrevCoupling = prevCouplingT;
        P.addParamTerm = [&](int n, const std::vector<double>& lam, std::vector<double>& grad) {
            std::vector<double> Rp, Rm;
            for (int j = 0; j < np; ++j) {
                sp[size_t(j)]->adjEval(adjRec_[size_t(n - 1)], adjRec_[size_t(n)], *this, Rp, nullptr, nullptr);
                sm[size_t(j)]->adjEval(adjRec_[size_t(n - 1)], adjRec_[size_t(n)], *this, Rm, nullptr, nullptr);
                double s = 0; for (int r = 0; r < m; ++r) s += lam[size_t(r)] * (Rp[size_t(r)] - Rm[size_t(r)]);
                grad[size_t(j)] += s / (2 * params[size_t(j)].h);
            }
        };
        P.addInitialTerm = [&](const std::vector<double>& lam1, std::vector<double>& grad) {
            std::vector<double> g0(static_cast<size_t>(m), 0.0);
            prevCouplingT(1, lam1, g0);   // (dR_1/dx_0)^T lam_1
            objPrev(1, g0);                // + dF/dx_0 through the lagged well data of step 1
            for (int j = 0; j < np; ++j) {
                double s = 0;
                for (int q = 0; q < n_; ++q) {
                    s += g0[size_t(q) * nv_] * dP0[size_t(j)][size_t(q)];
                    for (int i = 0; i < ncw_; ++i) s += g0[size_t(q) * nv_ + 1 + i] * dN0[size_t(j)][size_t(q) * ncw_ + i];
                }
                grad[size_t(j)] += s;
            }
        };
        const StepAdjointStats st = solveStepAdjoint(P, res.grad);
        res.ok = st.ok; res.error = st.error; res.linear = st.linear; res.steps = N;
    } catch (const std::exception& e) { res.ok = false; res.error = e.what(); }
    c_.wells = wSaved; pvEff_ = pvS; conns_ = cS; loadState(saved);
    res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return res;
}

}  // namespace comp
}  // namespace ressim
