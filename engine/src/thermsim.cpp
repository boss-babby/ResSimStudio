// Thermal (steam / hot water) simulator core - see thermsim.h for the formulation.
#include "ressim/thermsim.h"
#include "ressim/binio.h"
#include "ressim/steam.h"
#include "ressim/bsr.h"
#include "ressim/linsolve.h"
#include "ressim/eos.h"
#include "ressim/gasdata.h"
#include "ressim/co2brine.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <set>
#include <map>
#include <memory>
#include <chrono>
#include <stdexcept>
#include "ressim/stepadjoint.h"

namespace ressim {
namespace thermal {

namespace {
constexpr double kConRamp = 2e3;   // Pa, maximum width of the C1 drawdown ramp of injector connections
constexpr double kRgas = 8.314462618;
constexpr double kG = 9.80665;
constexpr int W = 0, O = 1, G = 2;   // phase slots

double interp(const std::vector<double>& x, const std::vector<double>& y, double v) {
    if (x.empty()) return 0;
    if (v <= x.front()) return y.front();
    if (v >= x.back()) return y.back();
    const size_t k = size_t(std::upper_bound(x.begin(), x.end(), v) - x.begin());
    const double t = (v - x[k - 1]) / (x[k] - x[k - 1]);
    return y[k - 1] + t * (y[k] - y[k - 1]);
}

double kValue(const Component& q, double p, double T) {
    if (q.type == 1) return 1e8;          // non-condensable: gas phase only
    if (!q.volatileK()) return 0.0;
    const double k = (q.kv1 + q.kv2 / p + q.kv3 * p) * std::exp(q.kv4 / (T - q.kv5));
    return std::max(k, 0.0);
}
double liqMolarVolume(const Case& c, const Component& q, double p, double T) {
    const double dT = T - c.tRef;
    const double den = std::max(1e-6, (1 - q.cpL * (p - c.pRef)) * (1 + q.ct1 * dT + q.ct2 * dT * dT));
    return q.mw * den / q.rhoRef;
}
double liqViscosity(const Component& q, double T) {
    if (!q.viscT.empty()) {
        if (q.viscT.size() == 1) return q.viscV[0];
        // log-linear in mu vs T
        std::vector<double> lv(q.viscV.size());
        for (size_t i = 0; i < lv.size(); ++i) lv[i] = std::log(std::max(q.viscV[i], 1e-12));
        return std::exp(interp(q.viscT, lv, T));
    }
    return q.avisc * std::exp(q.bvisc / T);
}
double gasViscosity(const Component& q, double T) { return q.avg * std::pow(T, q.bvg); }
double hvap(const Component& q, double T) {
    if (q.tcrit > 0) return T < q.tcrit ? q.hvapA * std::pow(1 - T / q.tcrit, q.hvapB) : 0.0;
    return q.hvapA;
}
// E300 component enthalpies relative to the reference temperature (J/mol)
bool hasOilCp(const Component& q) { return q.cpl1 != 0 || q.cpl2 != 0; }
bool hasGasCp(const Component& q) { return q.cpg1 != 0 || q.cpg2 != 0; }
double hGas(const Case& c, const Component& q, double T, double p);
double hLiq(const Case& c, const Component& q, double T) {
    const double d = T - c.tRef;
    if (q.tcrit > 0 && T >= q.tcrit && (hasGasCp(q) || hasOilCp(q))) return hGas(c, q, T, c.stdP);
    if (hasOilCp(q)) return q.cpl1 * d + 0.5 * q.cpl2 * d * d;
    if (hasGasCp(q)) return q.hvaps + q.cpg1 * d + 0.5 * q.cpg2 * d * d - hvap(q, T);
    return 0.0;
}
double hGas(const Case& c, const Component& q, double T, double p) {
    const double d = T - c.tRef;
    if (hasGasCp(q)) return q.hvaps + q.cpg1 * d + 0.5 * q.cpg2 * d * d - (hasOilCp(q) ? q.z1 * (p - c.stdP) : 0.0);
    if (hasOilCp(q)) return q.cpl1 * d + 0.5 * q.cpl2 * d * d + hvap(q, T);
    return 0.0;
}
double gasMolarVolume(const Component& q, double p, double T) { return std::max(q.z0 * kRgas * T / p - q.z1, 1e-12); }
double waterViscosity(const Case& c, double T, double rho) {
    if (!c.watViscT.empty()) {
        if (c.watViscT.size() == 1) return c.watViscV[0];
        std::vector<double> lv(c.watViscV.size());
        for (size_t i = 0; i < lv.size(); ++i) lv[i] = std::log(std::max(c.watViscV[i], 1e-12));
        return std::exp(interp(c.watViscT, lv, T));
    }
    return steam::viscosity(T, rho);
}
double steamViscosity(const Case& c, double T, double p, double rho) {
    if (c.steamViscIapws) return steam::viscosity(T, rho);
    const double mu = c.steamVisc[0] + c.steamVisc[1] * (T - 273.15) + c.steamVisc[2] * std::pow(std::max(p, 0.0) * 1e-6, c.steamVisc[3]);
    return std::max(mu, 1e-4) * 1e-3;
}

// standard Rachford-Rice: sum z (K - 1) / (1 + b (K - 1)) = 0; returns b clipped to [0, 1]
double rachfordRice(const std::vector<double>& z, const std::vector<double>& K) {
    double f0 = 0, s1 = 0;
    bool inf1 = false;
    for (size_t i = 0; i < z.size(); ++i) {
        if (z[i] <= 0) continue;
        f0 += z[i] * (K[i] - 1);
        if (K[i] <= 0) inf1 = true; else s1 += z[i] / K[i];
    }
    if (f0 <= 0) return 0.0;              // bubble point not reached: liquid
    if (!inf1 && s1 <= 1) return 1.0;     // dew point not reached: vapour
    double lo = 0, hi = 1;   // f(0) > 0 > f(1-): the root lies inside (0, 1), poles of K < 1 are beyond 1
    auto f = [&](double b) { double s = 0; for (size_t i = 0; i < z.size(); ++i) if (z[i] > 0) s += z[i] * (K[i] - 1) / (1 + b * (K[i] - 1)); return s; };
    double b = 0.5 * (lo + hi);
    for (int it = 0; it < 200; ++it) {
        const double v = f(b);
        if (v > 0) lo = b; else hi = b;
        double d = 0; for (size_t i = 0; i < z.size(); ++i) if (z[i] > 0) { const double t = 1 + b * (K[i] - 1); d -= z[i] * (K[i] - 1) * (K[i] - 1) / (t * t); }
        double bn = d < 0 ? b - v / d : 0.5 * (lo + hi);
        if (!(bn > lo && bn < hi)) bn = 0.5 * (lo + hi);
        if (std::fabs(bn - b) < 1e-15 || hi - lo < 1e-15) { b = bn; break; }
        b = bn;
    }
    return std::clamp(b, 0.0, 1.0);
}
} // namespace

// Water "vapour pressure" used by the flash K_w = p_v / p (Raoult) and H_wA = p_v. Below Tc: IF97 saturation line.
// Above Tc a C1 extension p_v = p_c exp(a (1 - Tc/T) + b (T - Tc)^2), a matched to d ln psat / d ln T at Tc and b
// = 2e-3 K^-2 so that a separate liquid-like water branch survives only ~10-20 K above Tc (IF97 region-1 liquid
// properties are extrapolated there) and vanishes continuously. Replaces the old jump K_w = psat/p -> 1e3 at Tc,
// whose discontinuity made near-critical cells (UCG cavity walls, deep SCW) oscillate between two flash branches.
static double waterVapourPressure(double T) {
    const double Tc = steam::kTc;
    if (T < Tc) return steam::psat(T);
    static const double pc0 = steam::psat(Tc - 1e-7);
    static const double a = [] {
        const double T1 = steam::kTc - 0.01, T2 = steam::kTc - 0.11;
        return (std::log(steam::psat(T1)) - std::log(steam::psat(T2))) / (0.1) * steam::kTc;   // d ln p / d ln T at Tc
    }();
    // d/dT [a (1 - Tc/T)] = a Tc / T^2 -> a / Tc at Tc, equal to d ln psat / dT there
    const double dT = T - Tc;
    return pc0 * std::exp(std::min(a * (1.0 - Tc / T) + 2e-3 * dT * dT, 60.0));
}

FlashState flashPT(const Case& c, double p, double T, const std::vector<double>& zin) {
    if (c.aqueousSolutes()) return flashGeneral(c, p, T, zin);
    const int nc = int(c.comps.size());
    FlashState f; f.x.assign(size_t(nc), 0.0); f.y.assign(size_t(nc) + 1, 0.0); f.a.assign(size_t(nc) + 1, 0.0);
    double tot = 0; for (double v : zin) tot += std::max(v, 0.0);
    if (!(tot > 0)) return f;
    std::vector<double> z(static_cast<size_t>(nc) + 1);
    for (int i = 0; i <= nc; ++i) z[size_t(i)] = std::max(zin[size_t(i)], 0.0) / tot;
    double zh = 0; for (int i = 0; i < nc; ++i) zh += z[size_t(i)];
    const double zw = z[size_t(nc)];
    const double Kw = std::min(waterVapourPressure(T) / p, 1e3);
    std::vector<double> K(static_cast<size_t>(nc));
    for (int i = 0; i < nc; ++i) K[size_t(i)] = kValue(c.comps[size_t(i)], p, T);
    auto finish = [&]() { f.nW *= tot; f.nO *= tot; f.nG *= tot; if (f.nW > 0) f.a[size_t(nc)] = 1.0; return f; };
    if (zh < 1e-14) {
        if (Kw < 1) f.nW = 1; else { f.nG = 1; f.y[size_t(nc)] = 1; }
        if (zh > 0) { f.nO = zh; f.nW = std::max(0.0, f.nW - zh); for (int i = 0; i < nc; ++i) f.x[size_t(i)] = z[size_t(i)] / zh; }
        return finish();
    }
    if (Kw < 1) {   // aqueous phase present (if enough water)
        std::vector<double> zp(static_cast<size_t>(nc)), Kp(static_cast<size_t>(nc));
        for (int i = 0; i < nc; ++i) { zp[size_t(i)] = z[size_t(i)] / zh; Kp[size_t(i)] = K[size_t(i)] / (1 - Kw); }
        const double b = rachfordRice(zp, Kp);
        const double V = b * zh / (1 - Kw), L = zh * (1 - b), Wq = zw - V * Kw;
        if (Wq >= 0) {
            f.nW = Wq; f.nO = L; f.nG = V;
            for (int i = 0; i < nc; ++i) {
                const double xi = zp[size_t(i)] / (1 + b * (Kp[size_t(i)] - 1));
                if (b < 1) f.x[size_t(i)] = xi;
                f.y[size_t(i)] = b < 1 ? K[size_t(i)] * xi : zp[size_t(i)] * (1 - Kw);
            }
            if (b >= 1) {   // oil composition for reporting only
                double s = 0; for (int i = 0; i < nc; ++i) s += K[size_t(i)] > 0 ? f.y[size_t(i)] / K[size_t(i)] : 0;
                for (int i = 0; i < nc; ++i) f.x[size_t(i)] = (s > 0 && K[size_t(i)] > 0) ? f.y[size_t(i)] / K[size_t(i)] / s : zp[size_t(i)];
            }
            f.y[size_t(nc)] = Kw;
            return finish();
        }
    }
    // no aqueous phase: water is a gas-only component
    //   g(V) = sum z_i (K_i - 1) / (1 + V (K_i - 1)) + z_w / V = 0, decreasing on (0, 1]
    auto g = [&](double V) {
        double s = zw > 0 ? zw / V : 0.0;
        for (int i = 0; i < nc; ++i) if (z[size_t(i)] > 0) s += z[size_t(i)] * (K[size_t(i)] - 1) / (1 + V * (K[size_t(i)] - 1));
        return s;
    };
    double V;
    double g1 = zw;
    bool inf1 = false;
    for (int i = 0; i < nc; ++i) if (z[size_t(i)] > 0) { if (K[size_t(i)] <= 0) inf1 = true; else g1 += z[size_t(i)] * (K[size_t(i)] - 1) / K[size_t(i)]; }
    if (!inf1 && g1 >= 0) V = 1.0;
    else if (zw <= 0 && g(0.0) <= 0) V = 0.0;
    else {
        double lo = 0, hi = 1;
        V = 0.5;
        for (int it = 0; it < 200; ++it) {
            const double v = g(V);
            if (v > 0) lo = V; else hi = V;
            double d = zw > 0 ? -zw / (V * V) : 0.0;
            for (int i = 0; i < nc; ++i) if (z[size_t(i)] > 0) { const double t = 1 + V * (K[size_t(i)] - 1); d -= z[size_t(i)] * (K[size_t(i)] - 1) * (K[size_t(i)] - 1) / (t * t); }
            double vn = d < 0 ? V - v / d : 0.5 * (lo + hi);
            if (!(vn > lo && vn < hi)) vn = 0.5 * (lo + hi);
            if (std::fabs(vn - V) < 1e-15 || hi - lo < 1e-15) { V = vn; break; }
            V = vn;
        }
    }
    f.nG = V; f.nO = 1 - V; f.nW = 0;
    for (int i = 0; i < nc; ++i) {
        const double xi = z[size_t(i)] / (1 + V * (K[size_t(i)] - 1));
        f.x[size_t(i)] = V < 1 ? xi : 0.0;
        f.y[size_t(i)] = V > 0 ? (V < 1 ? K[size_t(i)] * xi : z[size_t(i)]) : 0.0;
    }
    if (V > 0) f.y[size_t(nc)] = zw / V;
    if (V >= 1) { double s = 0; for (int i = 0; i < nc; ++i) s += z[size_t(i)]; for (int i = 0; i < nc; ++i) f.x[size_t(i)] = s > 0 ? z[size_t(i)] / s : 0; }
    return finish();
}

// ------------------------------------------------------------------------------------------------ C10: gas EOS, solutes
namespace {
// Peng-Robinson model of the non-water gas components with critical constants (Component::tcG > 0), cached per thread
// (keyed on the component data, so copies of a Case share the cache)
struct GasEos {
    eos::Model m;
    std::vector<int> idx;               // thermal component -> PR index (-1: not in the PR mixture)
    std::vector<char> isCO2;            // Islam-Carlson brine viscosity correction
    std::vector<double> key;
};
const GasEos& gasCache(const Case& c) {
    thread_local GasEos cache;
    thread_local bool valid = false;
    std::vector<double> key;
    key.reserve(c.comps.size() * 5 + 1);
    key.push_back(c.gasPR ? 1.0 : 0.0);
    for (const Component& q : c.comps) { key.push_back(q.tcG); key.push_back(q.pcG); key.push_back(q.omegaG); key.push_back(q.mw); key.push_back(double(q.type)); }
    if (valid && key == cache.key) return cache;
    cache = GasEos(); cache.key = key; valid = true;
    cache.idx.assign(c.comps.size(), -1); cache.isCO2.assign(c.comps.size(), 0);
    for (size_t i = 0; i < c.comps.size(); ++i) {
        const gasdata::GasInfo* gi = gasdata::find(c.comps[i].name);
        cache.isCO2[i] = gi && std::string(gi->name) == "CO2";
    }
    if (c.gasPR)
        for (size_t i = 0; i < c.comps.size(); ++i) {
            const Component& q = c.comps[i];
            if (q.solid() || !(q.tcG > 0) || !(q.pcG > 0)) continue;
            eos::Component e; e.name = q.name; e.tc = q.tcG; e.pc = q.pcG; e.omega = q.omegaG; e.mw = q.mw;
            e.zc = 0.3074; e.vc = 0.3074 * kRgas * q.tcG / q.pcG;   // LBC critical volume estimate (Zc of PR)
            cache.idx[i] = int(cache.m.comp.size());
            cache.m.comp.push_back(e);
        }
    return cache;
}
const GasEos* gasEos(const Case& c) { const GasEos& g = gasCache(c); return g.m.comp.empty() ? nullptr : &g; }
// PR state of the dry gas subset: fugacity coefficients (per thermal component, 1 for non-PR), molar volume, residual
// enthalpy and LBC viscosity of the PR subset; s = mole fraction of the PR subset in the gas
struct GasPR { double s = 0, vm = 0, hres = 0, mu = 0; std::vector<double> lnphi; };
bool gasPRState(const Case& c, const GasEos& ge, double p, double T, const std::vector<double>& y, GasPR& g, bool props) {
    const int nc = int(c.comps.size());
    g.lnphi.assign(size_t(nc), 0.0);
    std::vector<double> n(ge.m.comp.size(), 0.0);
    double s = 0;
    for (int i = 0; i < nc; ++i) if (ge.idx[size_t(i)] >= 0) { const double v = std::max(y[size_t(i)], 0.0); n[size_t(ge.idx[size_t(i)])] = v; s += v; }
    g.s = s;
    if (!(s > 1e-300)) {   // trace limit: fugacity coefficients at infinite dilution in ... (use the equimolar PR mixture)
        for (double& v : n) v = 1.0;
    }
    const eos::TParams tp = eos::atTemperature(ge.m, T);
    eos::EvalOptions o; o.dn = false; o.dP = false; o.dT = false; o.hres = props; o.root = eos::Root::Vapour;
    const eos::PhaseEval e = eos::evalPhase(ge.m, tp, p, n, o);
    for (int i = 0; i < nc; ++i) if (ge.idx[size_t(i)] >= 0) g.lnphi[size_t(i)] = e.lnphi[size_t(ge.idx[size_t(i)])];
    g.vm = e.vm; g.hres = e.hres;
    if (props) {
        std::vector<double> x(n.size());
        double t = 0; for (double v : n) t += v;
        for (size_t k = 0; k < n.size(); ++k) x[k] = n[k] / t;
        g.mu = eos::lbcViscosity(ge.m, x, T, 1.0 / std::max(e.vm, 1e-12));
        // C16: LBC over-predicts dense CO2 viscosity (x2 at 20 C / 500 bar); CO2-rich gas uses the Fenghour-Wakeham-Vesovic
        // (1998) reference correlation at the same molar density, blended as mu_LBC (mu_F / mu_LBC,CO2)^x_CO2
        for (int i = 0; i < nc; ++i) {
            const int k = ge.idx[size_t(i)];
            if (k < 0 || !ge.isCO2[size_t(i)] || !(x[size_t(k)] > 0)) continue;
            std::vector<double> xp(n.size(), 0.0); xp[size_t(k)] = 1.0;
            const double rhoM = 1.0 / std::max(e.vm, 1e-12);
            const double muL = eos::lbcViscosity(ge.m, xp, T, rhoM), muF = co2::co2Viscosity(T, rhoM * ge.m.comp[size_t(k)].mw);
            if (muL > 0 && muF > 0 && std::isfinite(muF)) g.mu *= std::pow(muF / muL, x[size_t(k)]);
            break;
        }
    }
    return std::isfinite(g.vm) && g.vm > 0;
}
// Henry constant (Pa) of a soluble component at T (saturation curve), its d ln kH / dT and the partial molar volume
double henryK(const Component& q, double T) {
    if (q.aqModel == 1) return gasdata::henryIapws(q.hA, q.hB, q.hC, T);
    if (q.aqModel == 2) return q.kh298 * std::exp(-q.hdT * (1.0 / T - 1.0 / 298.15));
    return 0.0;
}
double henryDlnDT(const Component& q, double T) {
    if (q.aqModel == 1) return gasdata::henryIapwsDlnDT(q.hA, q.hB, q.hC, T);
    if (q.aqModel == 2) return q.hdT / (T * T);
    return 0.0;
}
double vInfOf(const Component& q, double T) { return q.vInf < 0 ? gasdata::vInfCO2(T) : q.vInf; }
// aqueous fugacity constant H_iA = kH exp(vInf (p - psat) / RT) exp(ksalt m) (Pa)
double aqueousH(const Case& c, const Component& q, double p, double T) {
    const double ps = steam::psat(std::min(T, steam::kTc - 1e-6));
    return henryK(q, T) * std::exp(vInfOf(q, T) * (p - ps) / (kRgas * T) + q.ksalt * c.salinity);
}
}  // namespace

FlashState flashGeneral(const Case& c, double p, double T, const std::vector<double>& zin) {
    const int nc = int(c.comps.size()), nw = nc;
    FlashState f; f.x.assign(size_t(nc), 0.0); f.y.assign(size_t(nc) + 1, 0.0); f.a.assign(size_t(nc) + 1, 0.0);
    auto zOf = [&](int i) { return (i < nc && c.comps[size_t(i)].solid()) ? 0.0 : std::max(zin[size_t(i)], 0.0); };   // solids: not in the flash
    double tot = 0; for (int i = 0; i <= nc; ++i) tot += zOf(i);
    if (!(tot > 0)) return f;
    std::vector<double> z(static_cast<size_t>(nc) + 1);
    for (int i = 0; i <= nc; ++i) z[size_t(i)] = zOf(i) / tot;
    const GasEos* ge = gasEos(c);
    const double inf = std::numeric_limits<double>::infinity();
    // fugacity constants H_ij (x_ij = f_i / H_ij), phases j = W (aqueous), O, G; gas reference H_iG = phi_i p
    std::vector<double> lnphi(static_cast<size_t>(nc), 0.0), K(static_cast<size_t>(nc)), HA(static_cast<size_t>(nc), inf);
    for (int i = 0; i < nc; ++i) {
        const Component& q = c.comps[size_t(i)];
        K[size_t(i)] = kValue(q, p, T);
        if (q.soluble()) HA[size_t(i)] = aqueousH(c, q, p, T);
    }
    const double HwA = std::min(waterVapourPressure(T), 1e3 * p);    // Raoult water, C1 across Tc (K_w capped at 1e3)
    std::vector<double> a(static_cast<size_t>(nc + 1) * 3, 0.0);      // affinities a[i*3 + j]
    double beta[3] = {1.0 / 3, 1.0 / 3, 1.0 / 3};
    std::vector<double> E(static_cast<size_t>(nc) + 1);
    auto Qf = [&](const double* b, bool& ok) {
        double q = b[0] + b[1] + b[2];
        ok = true;
        for (int i = 0; i <= nc; ++i) {
            if (z[size_t(i)] <= 0) continue;
            const double e = b[0] * a[size_t(i) * 3] + b[1] * a[size_t(i) * 3 + 1] + b[2] * a[size_t(i) * 3 + 2];
            if (!(e > 0)) { ok = false; return inf; }
            q -= z[size_t(i)] * std::log(e);
        }
        return q;
    };
    std::vector<double> yv(static_cast<size_t>(nc) + 1, 0.0);
    for (int outer = 0; outer < 40; ++outer) {
        // affinities: a_ij = min_k H_ik / H_ij (infinite H -> 0, zero H -> the only phase)
        for (int i = 0; i <= nc; ++i) {
            double H[3];
            if (i == nw) { H[0] = HwA; H[1] = inf; H[2] = p; }
            else {
                const Component& q = c.comps[size_t(i)];
                const double hg = std::exp(lnphi[size_t(i)]) * p;
                if (q.solid()) { H[0] = H[1] = H[2] = inf; }
                else { H[2] = hg; H[1] = K[size_t(i)] * hg; H[0] = HA[size_t(i)]; }
            }
            double hmin = std::min(H[0], std::min(H[1], H[2]));
            for (int j = 0; j < 3; ++j) {
                double v;
                if (hmin <= 0) v = H[j] <= 0 ? 1.0 : 0.0;
                else v = std::isinf(H[j]) ? 0.0 : hmin / H[j];
                a[size_t(i) * 3 + size_t(j)] = v;
            }
        }
        // projected (damped) Newton on the convex Q(beta), beta >= 0
        bool okq; double Q = Qf(beta, okq);
        if (!okq) { beta[0] = beta[1] = beta[2] = 1.0 / 3; Q = Qf(beta, okq); }
        for (int it = 0; it < 200; ++it) {
            double g[3] = {1, 1, 1}, Hm[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
            for (int i = 0; i <= nc; ++i) {
                if (z[size_t(i)] <= 0) continue;
                const double* ai = &a[size_t(i) * 3];
                const double e = beta[0] * ai[0] + beta[1] * ai[1] + beta[2] * ai[2];
                E[size_t(i)] = e;
                for (int j = 0; j < 3; ++j) {
                    g[j] -= z[size_t(i)] * ai[j] / e;
                    for (int k = 0; k < 3; ++k) Hm[j][k] += z[size_t(i)] * ai[j] * ai[k] / (e * e);
                }
            }
            bool fr[3]; double gmax = 0;
            for (int j = 0; j < 3; ++j) {
                fr[j] = beta[j] > 0 || g[j] < 0;
                if (fr[j]) gmax = std::max(gmax, std::fabs(g[j]));
            }
            if (gmax < 1e-14) break;
            // solve (H_FF + mu I) d = -g_F; variables at their bound whose Newton step points outwards leave the free set
            double d[3] = {0, 0, 0};
            bool solved = false;
            for (int pass = 0; pass < 3 && !solved; ++pass) {
            int idx[3], nf = 0; for (int j = 0; j < 3; ++j) if (fr[j]) idx[nf++] = j;
            if (nf == 0) break;
            double tr = 0; for (int k = 0; k < nf; ++k) tr += Hm[idx[k]][idx[k]];
            std::fill(d, d + 3, 0.0);
            for (double mu = 0; !solved && mu < 1e3 * (tr + 1); mu = mu == 0 ? 1e-13 * (tr + 1e-300) : mu * 100) {
                double A[3][4];
                for (int r = 0; r < nf; ++r) { for (int q = 0; q < nf; ++q) A[r][q] = Hm[idx[r]][idx[q]] + (r == q ? mu : 0.0); A[r][nf] = -g[idx[r]]; }
                bool sing = false;
                for (int col = 0; col < nf && !sing; ++col) {
                    int pv = col; for (int r = col + 1; r < nf; ++r) if (std::fabs(A[r][col]) > std::fabs(A[pv][col])) pv = r;
                    if (!(std::fabs(A[pv][col]) > 1e-300)) { sing = true; break; }
                    for (int q = 0; q <= nf; ++q) std::swap(A[col][q], A[pv][q]);
                    for (int r = 0; r < nf; ++r) if (r != col) { const double m = A[r][col] / A[col][col]; for (int q = col; q <= nf; ++q) A[r][q] -= m * A[col][q]; }
                }
                if (sing) continue;
                for (int r = 0; r < nf; ++r) d[idx[r]] = A[r][nf] / A[r][r];
                double gd = 0; for (int j = 0; j < 3; ++j) gd += g[j] * d[j];
                if (gd < 0 && std::isfinite(gd)) solved = true;
            }
            if (solved) {
                bool drop = false;
                for (int j = 0; j < 3; ++j) if (fr[j] && beta[j] <= 0 && d[j] < 0) { fr[j] = false; drop = true; }
                if (drop) solved = false;
            }
            }
            if (!solved) { for (int j = 0; j < 3; ++j) d[j] = (beta[j] > 0 || g[j] < 0) ? -g[j] : 0.0; }   // projected steepest descent
            // largest step keeping beta >= 0
            double amax = 1.0; int hit = -1;
            for (int j = 0; j < 3; ++j) if (d[j] < 0) { const double t = -beta[j] / d[j]; if (t < amax) { amax = t; hit = j; } }
            double gd = 0; for (int j = 0; j < 3; ++j) gd += g[j] * d[j];
            double al = amax, bn[3]; bool accepted = false;
            for (int ls = 0; ls < 60; ++ls) {
                for (int j = 0; j < 3; ++j) bn[j] = std::max(0.0, beta[j] + al * d[j]);
                if (al == amax && hit >= 0) bn[hit] = 0.0;
                bool ok; const double Qn = Qf(bn, ok);
                if (ok && Qn <= Q + 1e-4 * al * gd + 1e-15 * std::fabs(Q)) { Q = Qn; accepted = true; break; }
                al *= 0.5;
            }
            if (!accepted) break;
            double ch = 0; for (int j = 0; j < 3; ++j) { ch = std::max(ch, std::fabs(bn[j] - beta[j])); beta[j] = bn[j]; }
            if (ch < 1e-16) break;
        }
        for (int i = 0; i <= nc; ++i) {
            const double* ai = &a[size_t(i) * 3];
            E[size_t(i)] = beta[0] * ai[0] + beta[1] * ai[1] + beta[2] * ai[2];
        }
        if (!ge) break;   // ideal gas: no fugacity iteration
        // gas composition (trial composition when the gas is absent) and new fugacity coefficients
        double sy = 0;
        for (int i = 0; i <= nc; ++i) { yv[size_t(i)] = z[size_t(i)] > 0 && E[size_t(i)] > 0 ? z[size_t(i)] * a[size_t(i) * 3 + 2] / E[size_t(i)] : 0.0; sy += yv[size_t(i)]; }
        if (sy > 0) for (double& v : yv) v /= sy;
        GasPR gp;
        if (!gasPRState(c, *ge, p, T, yv, gp, false)) break;
        double dl = 0;
        for (int i = 0; i < nc; ++i) { dl = std::max(dl, std::fabs(gp.lnphi[size_t(i)] - lnphi[size_t(i)])); lnphi[size_t(i)] = gp.lnphi[size_t(i)]; }
        if (dl < 1e-13) break;
    }
    // phase amounts and compositions
    double sx = 0, sa = 0, sy = 0;
    for (int i = 0; i <= nc; ++i) {
        if (z[size_t(i)] <= 0 || !(E[size_t(i)] > 0)) continue;
        const double ci = z[size_t(i)] / E[size_t(i)];
        if (i < nc) { f.x[size_t(i)] = ci * a[size_t(i) * 3 + 1]; sx += f.x[size_t(i)]; }
        f.y[size_t(i)] = ci * a[size_t(i) * 3 + 2]; sy += f.y[size_t(i)];
        f.a[size_t(i)] = ci * a[size_t(i) * 3]; sa += f.a[size_t(i)];
    }
    // present phases sum to 1 at convergence; normalise (absent phases: trial compositions for reporting)
    if (sx > 0) for (double& v : f.x) v /= sx;
    if (sy > 0) for (double& v : f.y) v /= sy;
    if (sa > 0) for (double& v : f.a) v /= sa;
    f.nW = beta[0] * tot; f.nO = beta[1] * tot; f.nG = beta[2] * tot;
    if (!(f.nW > 1e-300 * tot)) { f.nW = 0; std::fill(f.a.begin(), f.a.end(), 0.0); }
    if (!(f.nO > 1e-300 * tot)) f.nO = 0;
    if (!(f.nG > 1e-300 * tot)) { f.nG = 0; std::fill(f.y.begin(), f.y.end(), 0.0); }
    return f;
}

// ------------------------------------------------------------------------------------------------ phase properties
namespace {
struct Phase {
    bool on = false;
    double n = 0, v = 0, rho = 0, mu = 1, h = 0, u = 0;   // moles, m3/mol, kg/m3, Pa s, J/mol, J/mol
};
struct PhaseSet { Phase ph[3]; double T = 0; FlashState f; };

void phaseProps(const Case& c, double p, double T, const FlashState& f, PhaseSet& S) {
    const int nc = int(c.comps.size());
    S.T = T;
    Phase& w = S.ph[W]; Phase& o = S.ph[O]; Phase& g = S.ph[G];
    w = Phase(); o = Phase(); g = Phase();
    if (f.nW > 0) {
        const steam::Props L = steam::liquid(T, p);
        w.on = true; w.n = f.nW; w.v = L.v * steam::kMw; w.rho = 1.0 / L.v; w.h = L.h * steam::kMw; w.u = L.u * steam::kMw;
        w.mu = waterViscosity(c, T, w.rho);
        const double aw = f.a.size() > size_t(nc) ? f.a[size_t(nc)] : 1.0;
        if (aw < 1.0 && aw > 0) {   // dissolved gases: ideal mixing at infinite dilution (partial molar volumes / enthalpies)
            const GasEos& gc = gasCache(c);
            double v = aw * w.v, m = aw * steam::kMw, h = aw * w.h, u = aw * w.u, xco2 = 0;
            for (int i = 0; i < nc; ++i) {
                const double xa = f.a[size_t(i)];
                if (xa <= 0) continue;
                const Component& q = c.comps[size_t(i)];
                const double vi = vInfOf(q, T);
                const double hi = hGas(c, q, T, p) - kRgas * T * T * henryDlnDT(q, T);   // ideal gas + dissolution enthalpy
                v += xa * vi; m += xa * q.mw; h += xa * hi; u += xa * (hi - p * vi);
                if (gc.isCO2[size_t(i)]) xco2 += xa;
            }
            w.v = v; w.rho = m / v; w.h = h; w.u = u;
            if (xco2 > 0) w.mu *= 1 + 4.65 * std::pow(xco2, 1.0134);   // Islam & Carlson (2012) CO2-brine viscosity
        }
    }
    if (f.nO > 0) {
        o.on = true; o.n = f.nO;
        double v = 0, m = 0, lmu = 0, h = 0;
        for (int i = 0; i < nc; ++i) {
            const Component& q = c.comps[size_t(i)];
            const double x = f.x[size_t(i)];
            if (x <= 0) continue;
            v += x * liqMolarVolume(c, q, p, T); m += x * q.mw;
            lmu += x * std::log(liqViscosity(q, T)); h += x * hLiq(c, q, T);
        }
        o.v = std::max(v, 1e-12); o.rho = m / o.v; o.mu = std::exp(lmu); o.h = h; o.u = h - p * o.v;
    }
    if (f.nG > 0) {
        g.on = true; g.n = f.nG;
        const double yw = f.y[size_t(nc)];
        const double pw = yw * p;
        double zsum = 0, m = 0, h = 0, num = 0, den = 0;
        for (int i = 0; i < nc; ++i) {
            const Component& q = c.comps[size_t(i)];
            const double y = f.y[size_t(i)];
            if (y <= 0) continue;
            zsum += y * gasMolarVolume(q, p, T); m += y * q.mw; h += y * hGas(c, q, T, p);
            num += y * gasViscosity(q, T); den += y;   // E300: mole-fraction average
        }
        const GasEos* ge = c.gasPR ? gasEos(c) : nullptr;
        if (ge) {   // Peng-Robinson for the PR subset: replace its Z-factor volume / viscosity, add its residual enthalpy
            GasPR gp;
            if (gasPRState(c, *ge, p, T, f.y, gp, true) && gp.s > 0) {
                double vOld = 0, numOld = 0;
                for (int i = 0; i < nc; ++i) {
                    if (ge->idx[size_t(i)] < 0) continue;
                    const double y = f.y[size_t(i)];
                    if (y <= 0) continue;
                    vOld += y * gasMolarVolume(c.comps[size_t(i)], p, T); numOld += y * gasViscosity(c.comps[size_t(i)], T);
                }
                zsum += gp.s * gp.vm - vOld;
                h += gp.s * gp.hres;
                num += gp.s * gp.mu - numOld;
            }
        }
        double uCorr = 0;
        if (yw > 0) {
            // water in the gas: ideal mixing at its partial pressure; enthalpy / internal energy evaluated at
            // min(p_w, p on the 0.75 rho_c isochore) - keeps U(T) monotone near the water critical point (see
            // steam::vapourMixCap); the volume uses the true p_w
            const double pe = std::min(pw, steam::vapourMixCap(T));
            const steam::Props V = steam::vapour(T, pe);
            const double zw = steam::vapourZ(T, pw), zwe = pe < pw ? steam::vapourZ(T, pe) : zw;
            uCorr = yw * (zw - zwe) * kRgas * T;
            zsum += yw * zw * kRgas * T / p; m += yw * steam::kMw; h += yw * V.h * steam::kMw;
            const double rhoS = pw > 0 ? pw * steam::kMw / (zw * kRgas * T) : 0.0;
            num += yw * steamViscosity(c, T, p, rhoS); den += yw;
        }
        g.v = std::max(zsum, 1e-12); g.rho = m / g.v; g.h = h; g.u = h - p * g.v + uCorr;
        g.mu = den > 0 ? num / den : 1e-5;
    }
}
double fluidEnergy(const PhaseSet& S) {
    double e = 0;
    for (const Phase& q : S.ph) if (q.on) e += q.n * q.u;
    return e;
}
} // namespace

// ------------------------------------------------------------------------------------------------ local cell state
struct Simulator::Local {
    double T = 0, S[3] = {0, 0, 0};
    bool on[3] = {false, false, false};
    double rho[3] = {0, 0, 0}, xi[3] = {0, 0, 0}, mu[3] = {1, 1, 1}, kr[3] = {0, 0, 0};
    std::vector<double> x, y, a;
    double n[3] = {0, 0, 0};            // phase moles
    double pv = 0;
    int regime = 0;   // phase-presence bits (W 1, O 2, G 4) + 8 on the saturation line
};

// output layout: 0 T, 1 PV - V_fluid (m3), 2 lambda_eff, 3 total mobility sum kr/mu, then per phase (W, O, G):
// rho, e = kr/mu xi h, c_k = kr/mu xi x_k (k = 0..nc, water last)
bool Simulator::evalCell(int q, const double* X, double* Y, Local* L, double Tguess) const {
    const Case& c = c_;
    const int nc = nc_;
    const double p = X[0];
    if (!(p > 0) || !std::isfinite(p)) return false;
    std::vector<double> z(static_cast<size_t>(nc) + 1);
    double Nt = 0, Nall = 0;
    for (int k = 0; k <= nc; ++k) {
        const double v = std::max(X[1 + k], 0.0);
        Nall += v;
        if (k < nc && c.comps[size_t(k)].solid()) { z[size_t(k)] = 0; continue; }   // solids: not in the flash
        z[size_t(k)] = v; Nt += v;
    }
    if (!(Nt > 0)) return false;
    const double U = X[nb_ - 1];
    const double Cr = crock_[size_t(q)];
    auto solidEnergy = [&](double T) {
        if (!hasSolid_) return 0.0;
        double e = 0;
        for (int k = 0; k < nc; ++k) if (c.comps[size_t(k)].solid()) e += std::max(X[1 + k], 0.0) * hLiq(c, c.comps[size_t(k)], T);
        return e;
    };
    auto energyAt = [&](double T, PhaseSet& S) {
        S.f = flashPT(c, p, T, z);
        phaseProps(c, p, T, S.f, S);
        return fluidEnergy(S) + rockEnergy(q, T) + solidEnergy(T);
    };
    const double Cest = Cr + Nall * 80.0;
    const double tol = 1e-13 * (std::fabs(U) + Cest * 300.0);
    // phase split interpolated in energy between two flashed states (an energy jump that the temperature cannot resolve:
    // univariant water / steam, or near-univariant with a trace of volatile hydrocarbon in the gas)
    auto blend = [&](const PhaseSet& Sl, const PhaseSet& Sh, double th, double Tb, PhaseSet& Sout) {
        th = std::clamp(th, 0.0, 1.0);
        FlashState f;
        f.nW = (1 - th) * Sl.f.nW + th * Sh.f.nW;
        f.nO = (1 - th) * Sl.f.nO + th * Sh.f.nO;
        f.nG = (1 - th) * Sl.f.nG + th * Sh.f.nG;
        f.x.assign(size_t(nc), 0.0); f.y.assign(size_t(nc) + 1, 0.0);
        for (int i = 0; i < nc; ++i)
            if (f.nO > 0) f.x[size_t(i)] = ((1 - th) * Sl.f.nO * Sl.f.x[size_t(i)] + th * Sh.f.nO * Sh.f.x[size_t(i)]) / f.nO;
        for (int i = 0; i <= nc; ++i)
            if (f.nG > 0) f.y[size_t(i)] = ((1 - th) * Sl.f.nG * Sl.f.y[size_t(i)] + th * Sh.f.nG * Sh.f.y[size_t(i)]) / f.nG;
        f.a.assign(size_t(nc) + 1, 0.0);
        for (int i = 0; i <= nc; ++i)
            if (f.nW > 0) f.a[size_t(i)] = ((1 - th) * Sl.f.nW * (Sl.f.a.empty() ? 0.0 : Sl.f.a[size_t(i)]) + th * Sh.f.nW * (Sh.f.a.empty() ? 0.0 : Sh.f.a[size_t(i)])) / f.nW;
        Sout.f = f;
        phaseProps(c, p, Tb, Sout.f, Sout);
    };
    PhaseSet S;
    bool solved = false, univariant = false;
    // univariant check: an energy jump at Tsat(p) (water / steam with non-volatile or no hydrocarbons)
    if (z[size_t(nc)] > 1e-12 && p < 0.999 * steam::kPc) {
        const double Ts = steam::tsat(p);
        const double eps = 1e-7 * Ts;
        PhaseSet Sl, Sh;
        const double El = energyAt(Ts - eps, Sl), Eh = energyAt(Ts + eps, Sh);
        if (Eh - El > 1e-6 * Cest && U >= El && U <= Eh) {
            blend(Sl, Sh, (U - El) / (Eh - El), Ts, S);
            solved = true; univariant = true;
            if (L) L->regime |= 8;
        }
    }
    if (!solved) {
        double T = std::clamp(Tguess > 0 ? Tguess : 350.0, 260.0, 2400.0);
        PhaseSet St;
        double E = energyAt(T, St);
        if (std::fabs(E - U) <= tol) { S = St; solved = true; }
        else {
            double lo, hi, Elo, Ehi;
            double step = 2.0;
            if (E < U) { lo = T; Elo = E; hi = T; Ehi = E;
                while (Ehi < U) { lo = hi; Elo = Ehi; hi = std::min(hi + step, 2500.0); step *= 2; PhaseSet Sx; Ehi = energyAt(hi, Sx); if (hi >= 2500.0) break; }
            } else { hi = T; Ehi = E; lo = T; Elo = E;
                while (Elo > U) { hi = lo; Ehi = Elo; lo = std::max(lo - step, 250.0); step *= 2; PhaseSet Sx; Elo = energyAt(lo, Sx); if (lo <= 250.0) break; }
            }
            if (!(Elo <= U && U <= Ehi)) return false;
            // Illinois regula falsi on E(T) = U; true end-point energies / states are kept for the final blend
            int side = 0;
            double Tm = T, Em = E;
            double tElo = Elo, tEhi = Ehi;
            PhaseSet Slo, Shi;
            bool haveLo = false, haveHi = false, conv = false;
            for (int it = 0; it < 400; ++it) {
                Tm = (Ehi - Elo) > 0 ? lo + (U - Elo) * (hi - lo) / (Ehi - Elo) : 0.5 * (lo + hi);
                if (!(Tm > lo && Tm < hi)) Tm = 0.5 * (lo + hi);
                Em = energyAt(Tm, St);
                if (std::fabs(Em - U) <= tol) { conv = true; break; }
                if (Em < U) { lo = Tm; Elo = Em; tElo = Em; Slo = St; haveLo = true; if (side == -1) Ehi = U + 0.5 * (Ehi - U); side = -1; }
                else { hi = Tm; Ehi = Em; tEhi = Em; Shi = St; haveHi = true; if (side == 1) Elo = U - 0.5 * (U - Elo); side = 1; }
                // bracket down to a few ulps: the blend then interpolates between (almost) identical end states, so the result
                // is the same as an exact hit of E(T) = U (path independent - no noise from the starting temperature)
                if (hi - lo <= 8.0 * std::numeric_limits<double>::epsilon() * Tm) break;
            }
            if (!conv && haveLo && haveHi && tEhi > tElo) {
                // bracket collapsed on a (near-)jump of E(T): interpolate the phase split in energy (smooth in U)
                blend(Slo, Shi, (U - tElo) / (tEhi - tElo), 0.5 * (lo + hi), S);
                if (L) L->regime |= 16;
            } else S = St;
            solved = true;
        }
    }
    // C12 fold construction (near-critical water + gases). Ideal mixing of the real-gas components leaves a narrow fold of
    // the mixture energy E(T) near the water critical point (dE/dT < 0 between a local maximum E1 at T1 and a local
    // minimum E2 at T2 > T1): for U in [E2, E1] three temperatures satisfy E(T) = U and the root found depended on the
    // starting temperature. The fold is replaced by a two-state mixture, as for a first-order transition: the end states
    // are the outer-branch roots P (E(T_P) = E2, T_P < T1) and Q (E(T_Q) = E1, T_Q > T2), blended with
    // theta = (U - E2) / (E1 - E2). T(U), the phase split and the volumes are then single valued and continuous in U,
    // independent of the guess, and the end states move smoothly with p and z (E1, E2 are extremal values).
    static const bool noFold = std::getenv("RESSIM_NOFOLD") != nullptr;
    if (!noFold && solved && !univariant && z[size_t(nc)] > 1e-12 && Nt - z[size_t(nc)] > 1e-9 * Nt &&
        p > 0.98 * steam::kPc && p < 1.35 * steam::kPc && S.T > 615.0 && S.T < 680.0) {
        // fixed sample grid (independent of the root found): the fold detection does not flicker under perturbations
        const double a0 = 605.0, a1 = 690.0;
        const int M = 170;
        std::vector<double> Ts(size_t(M) + 1), Es(size_t(M) + 1);
        PhaseSet Sx;
        for (int i = 0; i <= M; ++i) { Ts[size_t(i)] = a0 + (a1 - a0) * i / M; Es[size_t(i)] = energyAt(Ts[size_t(i)], Sx); }
        int i1 = -1, i2 = -1;
        for (int i = 1; i < M && i1 < 0; ++i) if (Es[size_t(i)] > Es[size_t(i) - 1] && Es[size_t(i)] >= Es[size_t(i) + 1]) i1 = i;
        if (i1 > 0) for (int i = i1 + 1; i < M && i2 < 0; ++i) if (Es[size_t(i)] < Es[size_t(i) - 1] && Es[size_t(i)] <= Es[size_t(i) + 1]) i2 = i;
        if (i1 > 0 && i2 > 0) {
            auto golden = [&](double lo, double hi, bool mx, double& Tx) {
                const double g = 0.5 * (std::sqrt(5.0) - 1.0);
                double x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo);
                PhaseSet St2;
                double f1 = energyAt(x1, St2), f2 = energyAt(x2, St2);
                for (int it = 0; it < 40 && hi - lo > 1e-7; ++it) {
                    if (mx ? f1 > f2 : f1 < f2) { hi = x2; x2 = x1; f2 = f1; x1 = hi - g * (hi - lo); f1 = energyAt(x1, St2); }
                    else { lo = x1; x1 = x2; f1 = f2; x2 = lo + g * (hi - lo); f2 = energyAt(x2, St2); }
                }
                Tx = 0.5 * (lo + hi);
                return energyAt(Tx, St2);
            };
            double T1, T2;
            const double E1 = golden(Ts[size_t(i1) - 1], Ts[size_t(i1) + 1], true, T1);
            const double E2 = golden(Ts[size_t(i2) - 1], Ts[size_t(i2) + 1], false, T2);
            if (E1 - E2 > 1e-12 * Cest && U >= E2 && U <= E1) {
                // monotone root on [lo, hi] (E(lo) <= target <= E(hi)) by Illinois regula falsi
                auto root = [&](double lo, double hi, double target, PhaseSet& Sr) {
                    double flo = energyAt(lo, Sr) - target, fhi = energyAt(hi, Sr) - target;
                    if (flo > 0 || fhi < 0) return false;
                    int side = 0;
                    for (int it = 0; it < 200; ++it) {
                        double m = fhi - flo > 0 ? lo - flo * (hi - lo) / (fhi - flo) : 0.5 * (lo + hi);
                        if (!(m > lo && m < hi)) m = 0.5 * (lo + hi);
                        const double fm = energyAt(m, Sr) - target;
                        if (std::fabs(fm) <= tol || hi - lo <= 8.0 * std::numeric_limits<double>::epsilon() * m) return true;
                        if (fm < 0) { lo = m; flo = fm; if (side == -1) fhi *= 0.5; side = -1; }
                        else { hi = m; fhi = fm; if (side == 1) flo *= 0.5; side = 1; }
                    }
                    return true;
                };
                int jP = -1, kQ = -1;
                for (int j = i1 - 1; j >= 0 && jP < 0; --j) if (Es[size_t(j)] <= E2) jP = j;
                for (int k = i2 + 1; k <= M && kQ < 0; ++k) if (Es[size_t(k)] >= E1) kQ = k;
                PhaseSet SP, SQ;
                if (jP >= 0 && kQ >= 0 && root(Ts[size_t(jP)], T1, E2, SP) && root(T2, Ts[size_t(kQ)], E1, SQ)) {
                    const double th = (U - E2) / (E1 - E2);
                    blend(SP, SQ, th, (1 - th) * SP.T + th * SQ.T, S);
                    if (L) L->regime |= 32;
                }
            }
        }
    }
    const double T = S.T;
    // volumes, saturations
    const double pv = poreVolume(q, p, T);
    const double Vs = hasSolid_ ? solidVolume(q, X, p, T) : 0.0;
    double V[3], Vf = 0;
    for (int a = 0; a < 3; ++a) { V[a] = S.ph[a].on ? S.ph[a].n * S.ph[a].v : 0.0; Vf += V[a]; }
    double Sat[3];
    for (int a = 0; a < 3; ++a) Sat[a] = Vf > 0 ? V[a] / Vf : 0.0;
    // relative permeabilities (Corey or tables; ECLIPSE default three-phase oil)
    // ECLIPSE default three-phase oil model: krow(So) from the oil-water table at Sw' = 1 - So, krog(So) from the gas-oil
    // table at Sg' = 1 - So - Swc (liquid saturation So + Swc), kro = (Sg krog + (Sw - Swc) krow) / (Sg + Sw - Swc)
    const double sw = Sat[W], sg = Sat[G], so = std::max(0.0, 1 - sw - sg);
    double krw, krow, krg, krog;
    if (!c.swT.empty()) { krw = interp(c.swT, c.krwT, sw); krow = interp(c.swT, c.krowT, 1 - so); }
    else {
        const double den = std::max(1e-12, 1 - c.swc - c.sorw);
        krw = c.krw0 * std::pow(std::clamp((sw - c.swc) / den, 0.0, 1.0), c.nw);
        krow = c.kro0 * std::pow(std::clamp((so - c.sorw) / den, 0.0, 1.0), c.now);
    }
    if (!c.slT.empty()) { krg = interp(c.slT, c.krgT, 1 - sg); krog = interp(c.slT, c.krogT, so + c.swc); }
    else {
        const double dg = std::max(1e-12, 1 - c.swc - c.sorg - c.sgc), dl = std::max(1e-12, 1 - c.swc - c.sorg);
        krg = c.krg0 * std::pow(std::clamp((sg - c.sgc) / dg, 0.0, 1.0), c.ng);
        krog = c.kro0 * std::pow(std::clamp((so - c.sorg) / dl, 0.0, 1.0), c.nog);
    }
    if (c.hyst.on && size_t(q) < sgHist_.size() && sg < sgHist_[size_t(q)]) {
        // imbibition scanning curve (Land trapping, Carlson free-gas krg) from the historical maximum
        const double shy = sgHist_[size_t(q)];
        const double sgt = trappedGas(q, shy);
        if (sgt > 0 && shy > sgt) {
            double d = c.sgc;
            if (!c.slT.empty()) { d = 0; for (size_t k = 0; k < c.slT.size() && k < c.krgT.size(); ++k) if (c.krgT[k] <= 0) d = std::max(d, 1 - c.slT[k]); }
            const double i = std::max(c.hyst.sgtMax, d), mx = 1 - c.swc;
            if (i > d + 1e-9) {
                const double C = 1.0 / (i - d) - 1.0 / std::max(1e-9, mx - d);
                if (sg <= sgt) krg = 0;
                else {
                    const double aa = sg - sgt, sgf = 0.5 * (aa + std::sqrt(aa * aa + 4.0 / C * aa)), se = d + sgf;
                    if (!c.slT.empty()) krg = interp(c.slT, c.krgT, 1 - se);
                    else krg = c.krg0 * std::pow(std::clamp((se - c.sgc) / std::max(1e-12, 1 - c.swc - c.sorg - c.sgc), 0.0, 1.0), c.ng);
                }
            }
        }
    }
    double kro;
    const double swe = std::max(0.0, sw - c.swc);
    if (sg + swe > 1e-12) kro = (sg * krog + swe * krow) / (sg + swe);
    else kro = std::min(krow, krog);
    double kr[3] = {S.ph[W].on ? krw : 0.0, S.ph[O].on ? kro : 0.0, S.ph[G].on ? krg : 0.0};
    // capillary pressure: phase pressure offsets (p_w = p - Pcow, p_g = p + Pcgo)
    double pcw = 0, pcg = 0;
    if (hasPc_) {
        if (!c.pcowT.empty() && c.pcowT.size() == c.swT.size()) pcw = interp(c.swT, c.pcowT, sw);
        else if (c.pcw0 != 0) pcw = c.pcw0 * (1 - std::clamp((sw - c.swc) / std::max(1e-12, 1 - c.swc - c.sorw), 0.0, 1.0));
        if (!c.pcgoT.empty() && c.pcgoT.size() == c.slT.size()) pcg = interp(c.slT, c.pcgoT, 1 - sg);
        else if (c.pcg0 != 0) pcg = c.pcg0 * std::clamp((sg - c.sgc) / std::max(1e-12, 1 - c.swc - c.sorg - c.sgc), 0.0, 1.0);
    }
    const double pcOff[3] = {-pcw, 0.0, pcg};
    Y[0] = T;
    Y[1] = pv - Vs - Vf;
    Y[2] = std::max(0.0, thc_[size_t(q)] * (1 - thcSf_[size_t(q)] * Sat[G]));
    double lt = 0;
    for (int a = 0; a < 3; ++a) if (S.ph[a].on) lt += kr[a] / S.ph[a].mu;
    Y[3] = lt;
    for (int a = 0; a < 3; ++a) {
        double* y = Y + 4 + a * stride_;
        const Phase& P = S.ph[a];
        y[0] = P.on ? P.rho : 0.0;
        y[1] = pcOff[a];
        const double mob = P.on ? kr[a] / P.mu / P.v : 0.0;   // molar mobility
        y[2] = mob * P.h;
        for (int k = 0; k <= nc; ++k) {
            double frac = 0;
            if (P.on) {
                if (a == W) frac = S.f.a[size_t(k)];
                else if (a == O) frac = k < nc ? S.f.x[size_t(k)] : 0.0;
                else frac = S.f.y[size_t(k)];
            }
            y[3 + k] = mob * frac;
        }
    }
    // reactions: component sources (mol/s) and heat release (W)
    if (rbase_ > 0) {
        double* yr = Y + rbase_;
        std::fill(yr, yr + nc + 2, 0.0);
        const double Vb = bulk_[size_t(q)];
        auto conc = [&](const Reaction& rx, int i) -> double {
            if (i < nc && c.comps[size_t(i)].solid()) return std::max(X[1 + i], 0.0) / Vb;
            if (rx.o2pp == i) return S.ph[G].on ? std::max(S.f.y[size_t(i)], 0.0) * p * 1e-3 : 0.0;
            int ph = rx.phase.empty() ? -1 : rx.phase[size_t(i)];
            if (ph < 0) ph = i == nc ? W : (c.comps[size_t(i)].type == 1 ? G : O);
            if (!S.ph[ph].on) {
                if (i == nc && ph == W && S.ph[G].on) ph = G; else return 0.0;   // water reacts as steam without an aqueous phase
            }
            double xf = ph == W ? S.f.a[size_t(i)] : ph == O ? (i < nc ? S.f.x[size_t(i)] : 0.0) : S.f.y[size_t(i)];
            return S.ph[ph].n * std::max(xf, 0.0) / Vb;
        };
        for (const Reaction& rx : c.reactions) {
            if (T < rx.tMin) continue;
            const double k = rx.A * std::exp(-rx.Ea / (kRgas * T));
            double f = 1, b = 1;
            for (int i = 0; i <= nc; ++i) {
                const double n1 = rx.order.empty() ? rx.sto[size_t(i)] : rx.order[size_t(i)];
                if (n1 > 0) f *= std::pow(conc(rx, i), n1);
                if (rx.reversible && rx.prod[size_t(i)] > 0) b *= std::pow(conc(rx, i), rx.prod[size_t(i)]);
            }
            double r = k * f;
            if (rx.reversible) r = k * (f - b / std::exp(rx.keqA + rx.keqB / T));
            r *= Vb;
            if (r == 0 || !std::isfinite(r)) continue;
            for (int i = 0; i <= nc; ++i) yr[i] += (rx.prod[size_t(i)] - rx.sto[size_t(i)]) * r;
            yr[nc + 1] += rx.dH * r;
        }
    }
    if (L) {
        L->T = T; L->pv = pv;
        for (int a = 0; a < 3; ++a) {
            L->S[a] = Sat[a]; L->on[a] = S.ph[a].on; L->rho[a] = S.ph[a].rho; L->mu[a] = S.ph[a].mu; L->kr[a] = kr[a];
            L->xi[a] = S.ph[a].on ? 1.0 / S.ph[a].v : 0.0;
            if (S.ph[a].on) L->regime |= 1 << a;
        }
        L->x = S.f.x; L->y = S.f.y; L->a = S.f.a;
        for (int a = 0; a < 3; ++a) L->n[a] = S.ph[a].on ? S.ph[a].n : 0.0;
    }
    return true;
}

// ------------------------------------------------------------------------------------------------ setup
bool Simulator::setup(const Case& cin, std::string* error) {
    auto fail = [&](const std::string& m) { if (error) *error = m; return false; };
    c_ = cin;
    Case& c = c_;
    nc_ = int(c.comps.size());
    sgHist_.clear();
    if (c.hyst.on && !(c.hyst.sgtMax > 0 && c.hyst.sgtMax < 1)) return fail("gas hysteresis: the maximum trapped gas saturation must lie in (0, 1)");
    for (const Component& q : c.comps) {
        if (q.aqModel == 2 && !(q.kh298 > 0)) return fail("component " + q.name + ": van't Hoff Henry constant must be > 0");
        if (c.gasPR && q.tcG > 0 && !(q.pcG > 0)) return fail("component " + q.name + ": Peng-Robinson needs a critical pressure");
    }
    nb_ = nc_ + 3;
    if (nb_ > kMaxBlock) return fail("too many components for the thermal core (max " + std::to_string(kMaxBlock - 3) + ")");
    stride_ = nc_ + 4;
    hasSolid_ = false;
    for (const Component& q : c.comps) if (q.solid()) hasSolid_ = true;
    for (const Reaction& rx : c.reactions) {
        if (rx.sto.size() != size_t(nc_) + 1 || rx.prod.size() != size_t(nc_) + 1) return fail("reaction " + rx.name + ": stoichiometry needs nc + 1 entries (water last)");
        if (!rx.order.empty() && rx.order.size() != size_t(nc_) + 1) return fail("reaction " + rx.name + ": orders need nc + 1 entries");
        if (!rx.phase.empty() && rx.phase.size() != size_t(nc_) + 1) return fail("reaction " + rx.name + ": phases need nc + 1 entries");
        if (rx.o2pp > nc_) return fail("reaction " + rx.name + ": O2PP component out of range");
    }
    hasPc_ = c.pcw0 != 0 || c.pcg0 != 0 || !c.pcowT.empty() || !c.pcgoT.empty();
    rbase_ = c.reactions.empty() ? 0 : 4 + 3 * stride_;
    m_ = 4 + 3 * stride_ + (c.reactions.empty() ? 0 : nc_ + 2);
    n_ = c.size();
    if (n_ <= 0) return fail("empty grid");
    const int n = n_;
    auto at = [](const std::vector<double>& v, int i, double d) { return v.empty() ? d : v.size() == 1 ? v[0] : v[size_t(i)]; };
    pv0_.assign(size_t(n), 0); depth_.assign(size_t(n), 0); bulk_.assign(size_t(n), 0);
    conns_.clear();
    if (!c.cellPV.empty()) {
        for (int i = 0; i < n; ++i) {
            pv0_[size_t(i)] = c.cellPV[size_t(i)];
            depth_[size_t(i)] = at(c.cellDepth, i, 0);
            bulk_[size_t(i)] = c.cellBulk.empty() ? pv0_[size_t(i)] / 0.3 : c.cellBulk[size_t(i)];
        }
        conns_ = c.conns;
    } else {
        const int nx = c.nx, ny = c.ny, nz = c.nz;
        auto idx = [&](int i, int j, int k) { return (k * ny + j) * nx + i; };
        std::vector<double> DX(static_cast<size_t>(n)), DY(static_cast<size_t>(n)), DZ(static_cast<size_t>(n));
        for (int k = 0; k < nz; ++k) {
            double ztop = c.top;
            for (int kk = 0; kk < k; ++kk) ztop += at(c.dz, idx(0, 0, kk), 1.0);
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    const int q = idx(i, j, k);
                    DX[size_t(q)] = at(c.dx, q, 1.0); DY[size_t(q)] = at(c.dy, q, 1.0); DZ[size_t(q)] = at(c.dz, q, 1.0);
                    bulk_[size_t(q)] = DX[size_t(q)] * DY[size_t(q)] * DZ[size_t(q)];
                    pv0_[size_t(q)] = bulk_[size_t(q)] * at(c.poro, q, 0.3);
                    depth_[size_t(q)] = ztop + 0.5 * DZ[size_t(q)];
                }
        }
        auto addConn = [&](int a, int b, int dir) {
            double ka, kb, la, lb, A;
            if (dir == 0) { ka = at(c.permx, a, 1e-13); kb = at(c.permx, b, 1e-13); la = DX[size_t(a)]; lb = DX[size_t(b)]; A = DY[size_t(a)] * DZ[size_t(a)]; }
            else if (dir == 1) { ka = at(c.permy, a, at(c.permx, a, 1e-13)); kb = at(c.permy, b, at(c.permx, b, 1e-13)); la = DY[size_t(a)]; lb = DY[size_t(b)]; A = DX[size_t(a)] * DZ[size_t(a)]; }
            else { ka = at(c.permz, a, 0.1 * at(c.permx, a, 1e-13)); kb = at(c.permz, b, 0.1 * at(c.permx, b, 1e-13)); la = DZ[size_t(a)]; lb = DZ[size_t(b)]; A = DX[size_t(a)] * DY[size_t(a)]; }
            const double T = (ka > 0 && kb > 0) ? A / (0.5 * la / ka + 0.5 * lb / kb) : 0.0;
            conns_.push_back({a, b, T, A / (0.5 * la + 0.5 * lb)});
        };
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    const int q = idx(i, j, k);
                    if (i + 1 < nx) addConn(q, idx(i + 1, j, k), 0);
                    if (j + 1 < ny) addConn(q, idx(i, j + 1, k), 1);
                    if (k + 1 < nz) addConn(q, idx(i, j, k + 1), 2);
                }
        if (c.heatLoss && c.lossFaces.empty())
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    const int t = idx(i, j, 0), b = idx(i, j, nz - 1);
                    c.lossFaces.push_back({t, DX[size_t(t)] * DY[size_t(t)]});
                    c.lossFaces.push_back({b, DX[size_t(b)] * DY[size_t(b)]});
                }
    }
    if (!c.gravity) std::fill(depth_.begin(), depth_.end(), 0.0);
    crock_.assign(size_t(n), 0);
    crock1_.assign(size_t(n), 0); thc_.assign(size_t(n), 0); thcSf_.assign(size_t(n), 0);
    for (int i = 0; i < n; ++i) {
        const double vr = std::max(0.0, bulk_[size_t(i)] - pv0_[size_t(i)]);
        crock_[size_t(i)] = vr * at(c.rockCvCell, i, c.rockCv); crock1_[size_t(i)] = vr * at(c.rockCvTCell, i, c.rockCvT);
        const double phi = bulk_[size_t(i)] > 0 ? std::min(1.0, pv0_[size_t(i)] / bulk_[size_t(i)]) : 0.3;
        const double kl = 0.5 * (c.lamO + c.lamW);
        const double th = c.thconr.empty() ? phi * kl + (1 - phi) * c.lamRock : at(c.thconr, i, 0.0);
        thc_[size_t(i)] = th;
        const double sf = at(c.thconsfCell, i, c.thconsf);
        thcSf_[size_t(i)] = sf >= 0 ? sf : (th > 0 ? phi * (kl - c.lamG) / th : 0.0);
    }
    // pore volume law: rock compressibility / thermal expansion (replaced by the fixed-stress form with geomechanics)
    pvA_ = pv0_; pvCp_.assign(size_t(n), 0); pvCt_.assign(size_t(n), 0); pvP_.assign(size_t(n), c.pref); pvT_.assign(size_t(n), c.tRef);
    for (int i = 0; i < n; ++i) { pvCp_[size_t(i)] = pv0_[size_t(i)] * c.cr; pvCt_[size_t(i)] = pv0_[size_t(i)] * c.ct; }
    // initial state: moles from saturations at (p0, T0), U from the energy at T0, p adjusted to the volume balance
    nb_ = nc_ + 3;
    X_.assign(size_t(n) * nb_, 0.0);
    p_.assign(size_t(n), 0); T_.assign(size_t(n), 0); sw_.assign(size_t(n), 0); so_.assign(size_t(n), 0); sg_.assign(size_t(n), 0);
    for (int q = 0; q < n; ++q) {
        double p = at(c.p0, q, 1e7);
        const double T = at(c.T0, q, 300.0), sw = at(c.sw0, q, c.swc), sg = at(c.sg0, q, 0.0);
        std::vector<double> x = c.x0.empty() ? std::vector<double>(size_t(nc_), nc_ ? 1.0 / nc_ : 0.0) : (c.x0.size() == 1 ? c.x0[0] : c.x0[size_t(q)]);
        std::vector<double> y = c.y0.empty() ? x : (c.y0.size() == 1 ? c.y0[0] : c.y0[size_t(q)]);
        if (int(x.size()) != nc_ || int(y.size()) != nc_) return fail("initial composition size differs from the number of components");
        double* X = &X_[size_t(q) * nb_];
        // solids (mol per m3 of pore volume) occupy part of the pore space
        if (hasSolid_ && !c.solid0.empty()) {
            const std::vector<double>& s0 = c.solid0.size() == 1 ? c.solid0[0] : c.solid0[size_t(q)];
            if (int(s0.size()) != nc_) return fail("initial solid concentrations need nc entries");
            for (int i = 0; i < nc_; ++i) if (c.comps[size_t(i)].solid()) X[1 + i] = std::max(0.0, s0[size_t(i)]) * poreVolume(q, p, T);
        }
        const double pv = poreVolume(q, p, T) - (hasSolid_ ? solidVolume(q, X, p, T) : 0.0);
        if (!(pv > 0)) return fail("solids fill the pore volume of cell " + std::to_string(q));
        const double vw = steam::liquid(T, p).v * steam::kMw;
        for (int i = 0; i < nc_; ++i) if (c.comps[size_t(i)].solid()) { x[size_t(i)] = 0; y[size_t(i)] = 0; }
        double vo = 0; for (int i = 0; i < nc_; ++i) vo += x[size_t(i)] * liqMolarVolume(c, c.comps[size_t(i)], p, T);
        double vg = 0; for (int i = 0; i < nc_; ++i) vg += y[size_t(i)] * gasMolarVolume(c.comps[size_t(i)], p, T);
        if (const GasEos* ge = c.gasPR ? gasEos(c) : nullptr) {   // PR subset volume
            std::vector<double> yy(size_t(nc_) + 1, 0.0); for (int i = 0; i < nc_; ++i) yy[size_t(i)] = y[size_t(i)];
            GasPR gp;
            if (gasPRState(c, *ge, p, T, yy, gp, false) && gp.s > 0) {
                double vOld = 0; for (int i = 0; i < nc_; ++i) if (ge->idx[size_t(i)] >= 0) vOld += y[size_t(i)] * gasMolarVolume(c.comps[size_t(i)], p, T);
                vg += gp.s * gp.vm - vOld;
            }
        }
        if (!(vg > 0)) vg = kRgas * T / p;
        const double so = std::max(0.0, 1 - sw - sg);
        const double nO = nc_ ? so * pv / std::max(vo, 1e-12) : 0.0, nG = nc_ ? sg * pv / vg : 0.0;
        for (int i = 0; i < nc_; ++i) if (!c.comps[size_t(i)].solid()) X[1 + i] = nO * x[size_t(i)] + nG * y[size_t(i)];
        X[1 + nc_] = sw * pv / vw + (nc_ ? 0.0 : (1 - sw) * pv / vw);
        // energy at T0
        std::vector<double> z(X + 1, X + 2 + nc_);
        double eSol = 0;
        for (int i = 0; i < nc_; ++i) if (c.comps[size_t(i)].solid()) { eSol += z[size_t(i)] * hLiq(c, c.comps[size_t(i)], T); z[size_t(i)] = 0; }
        PhaseSet S;
        auto E = [&](double pp) { S.f = flashPT(c, pp, T, z); phaseProps(c, pp, T, S.f, S); return fluidEnergy(S) + rockEnergy(q, T) + eSol; };
        // pressure from the volume balance at fixed moles and T (secant)
        auto vres = [&](double pp) {
            S.f = flashPT(c, pp, T, z); phaseProps(c, pp, T, S.f, S);
            double Vf = 0; for (const Phase& ph : S.ph) if (ph.on) Vf += ph.n * ph.v;
            return poreVolume(q, pp, T) - (hasSolid_ ? solidVolume(q, X, pp, T) : 0.0) - Vf;
        };
        double p1 = p, f1 = vres(p1), p2 = p * 1.001, f2 = vres(p2);
        for (int it = 0; it < 40 && std::fabs(f2) > 1e-12 * pv; ++it) {
            if (f2 == f1) break;
            double pn = p2 - f2 * (p2 - p1) / (f2 - f1);
            pn = std::clamp(pn, 0.5 * p2, 2.0 * p2);
            p1 = p2; f1 = f2; p2 = pn; f2 = vres(p2);
        }
        if (std::isfinite(p2) && p2 > 0 && std::fabs(p2 / p - 1) < 0.5) p = p2;
        X[0] = p;
        X[nb_ - 1] = E(p);
        p_[size_t(q)] = p; T_[size_t(q)] = T;
    }
    Xn_ = X_;
    Y_.assign(size_t(n) * m_, 0.0); dY_.assign(size_t(n) * m_ * nb_, 0.0);
    loss_.assign(c.lossFaces.size(), LossState());
    for (size_t f = 0; f < c.lossFaces.size(); ++f) loss_[f].Ti = c.obT > 0 ? c.obT : T_[size_t(c.lossFaces[f].cell)];
    rep_.assign(c.wells.size(), WellReport());
    J_ = BSR();
    surfYield_.assign(c.wells.size(), std::vector<double>());
    wbRho_.assign(c.wells.size(), -1.0);
    for (size_t w = 0; w < c.wells.size(); ++w) {
        const Well& wl = c.wells[w];
        if (wl.cells.size() != wl.wi.size()) return fail("well " + wl.name + ": connection and WI counts differ");
        for (int q : wl.cells) if (q < 0 || q >= n) return fail("well " + wl.name + ": connection outside the grid");
        if (wl.injector && int(wl.zinj.size()) != nc_ + 1) return fail("well " + wl.name + ": injected composition needs nc + 1 entries (water last)");
        rep_[w].molar.assign(size_t(nc_) + 1, 0.0);
    }
    t_ = 0; st_ = Stats(); cumLoss_ = 0; cumHin_ = 0;
    if (!evalAll()) return fail("initial state evaluation failed");
    for (int q = 0; q < n; ++q) T_[size_t(q)] = Y_[size_t(q) * m_];
    p0_ = p_; T0init_ = T_;
    cumReact_ = 0; reactRate_ = 0;
    sgHist_ = sg_;
    solidInit_.assign(size_t(n) * size_t(nc_), 0.0);
    for (int q = 0; q < n; ++q) for (int i = 0; i < nc_; ++i) if (c.comps[size_t(i)].solid()) solidInit_[size_t(q) * nc_ + i] = X_[size_t(q) * nb_ + 1 + i];
    if (c.track >= nc_) return fail("tracked component out of range");
    // permeability alteration: base transmissibilities and initial fluid porosity
    connT0_.resize(conns_.size());
    for (size_t e = 0; e < conns_.size(); ++e) connT0_[e] = conns_[e].T;
    permMult_.assign(size_t(n), 1.0);
    phiF0_ = fluidPorosity();
    // geomechanics: fixed-stress pore volume law centred at the initial state
    geo_ = geomech::Solver();
    gmap_.clear(); gn_ = 0;
    if (c.geo.enabled) {
        std::vector<double> phi = porosity();
        std::string err;
        if (c.cellPV.empty()) {
            gn_ = n; gmap_.resize(size_t(n));
            for (int q = 0; q < n; ++q) gmap_[size_t(q)] = q;
            if (!geo_.setup(c.geo, c.nx, c.ny, c.nz, c.dx, c.dy, c.dz, c.top, p_, phi, &err)) return fail(err);
        } else {
            gn_ = c.geoNx * c.geoNy * c.geoNz;
            if (gn_ <= 0 || int(c.geoCell.size()) != n) return fail("geomechanics: the mechanics box (geoNx / geoCell) is missing");
            gmap_ = c.geoCell;
            double pm = 0; for (double v : p_) pm += v; pm /= std::max(1, n);
            std::vector<double> pg(size_t(gn_), pm), fg(size_t(gn_), 0.2);
            for (int q = 0; q < n; ++q) {
                const int g = gmap_[size_t(q)];
                if (g == -1) continue;   // unmapped (EDFM fracture cell)
                if (g < 0 || g >= gn_) return fail("geomechanics: cell map outside the mechanics box");
                pg[size_t(g)] = p_[size_t(q)]; fg[size_t(g)] = phi[size_t(q)];
            }
            const bool ok = c.geoNodes.empty() ? geo_.setup(c.geo, c.geoNx, c.geoNy, c.geoNz, c.geoDx, c.geoDy, c.geoDz, c.geoTop, pg, fg, &err)
                                               : geo_.setupHex(c.geo, c.geoNx, c.geoNy, c.geoNz, c.geoNodes, pg, fg, &err);
            if (!ok) return fail(err);
        }
        for (int q = 0; q < n; ++q) {
            const double Vb = bulk_[size_t(q)];
            const int g = gmap_[size_t(q)];
            if (g < 0) continue;
            pvA_[size_t(q)] = poreVolume(q, p_[size_t(q)], T_[size_t(q)]);
            pvP_[size_t(q)] = p_[size_t(q)]; pvT_[size_t(q)] = T_[size_t(q)];
            pvCp_[size_t(q)] = Vb * geo_.cp(g); pvCt_[size_t(q)] = Vb * geo_.cT(g);
        }
        if (!evalAll()) return fail("initial state evaluation failed (geomechanics)");
    }
    // C12: EDFM fracture apertures (Barton-Bandis / Coulomb) and their pore-volume linearisation
    fracW_.clear();
    for (const auto& f : c.fracCells) {
        if (f.cell < 0 || f.cell >= n) return fail("fracture cell index out of range");
        fracW_.push_back(f.w0);
    }
    for (int e : c.fracFF) if (e < 0 || e >= int(conns_.size())) return fail("fracture connection index out of range");
    if (!c.fracCells.empty()) { updateFractures(false); if (!evalAll()) return fail("initial state evaluation failed (fractures)"); }
    return true;
}

double Simulator::solidVolume(int q, const double* X, double p, double T) const {
    (void)q;
    double v = 0;
    for (int k = 0; k < nc_; ++k) {
        const Component& s = c_.comps[size_t(k)];
        if (s.solid()) v += std::max(X[1 + k], 0.0) * liqMolarVolume(c_, s, p, T);
    }
    return v;
}
std::vector<double> Simulator::porosity() const {
    std::vector<double> o(static_cast<size_t>(n_), 0.0);
    for (int q = 0; q < n_; ++q) o[size_t(q)] = bulk_[size_t(q)] > 0 ? poreVolume(q, p_[size_t(q)], T_[size_t(q)]) / bulk_[size_t(q)] : 0.0;
    return o;
}
std::vector<double> Simulator::fluidPorosity() const {
    std::vector<double> o = porosity();
    if (hasSolid_)
        for (int q = 0; q < n_; ++q)
            if (bulk_[size_t(q)] > 0) o[size_t(q)] -= solidVolume(q, &X_[size_t(q) * nb_], p_[size_t(q)], T_[size_t(q)]) / bulk_[size_t(q)];
    return o;
}
void Simulator::updatePermeability() {
    const Case& c = c_;
    // C16: damage permeability of yielded rock (geomechanics plasticity) multiplies the porosity law
    const bool damage = geo_.ready() && ((c.geo.plasticity == 1 && c.geo.permPlastic > 0) || c.geo.permStress > 0);   // C25: + stress-dependent k
    if (c.permLaw == 0 && !damage) return;
    const std::vector<double> phi = fluidPorosity();
    for (int q = 0; q < n_; ++q) {
        const double f0 = std::max(phiF0_[size_t(q)], 1e-6), f = std::clamp(phi[size_t(q)], 1e-6, 0.9999);
        double m = 1;
        if (c.permLaw == 1) m = std::pow(f / f0, 3) * std::pow((1 - f0) / (1 - f), 2);
        else if (c.permLaw != 0) m = std::exp(c.permC * (f - f0));
        double mx = c.permLaw ? c.permMultMax : 1.0;
        if (damage && gmap_[size_t(q)] >= 0) { m *= geo_.permMultiplier()[size_t(gmap_[size_t(q)])]; mx = std::max(mx, c.geo.permPlasticMax); }
        permMult_[size_t(q)] = std::clamp(m, 1.0 / std::max(mx, 1.0), std::max(mx, 1.0));
    }
    for (size_t e = 0; e < conns_.size(); ++e) {
        const double ma = permMult_[size_t(conns_[e].a)], mb = permMult_[size_t(conns_[e].b)];
        conns_[e].T = connT0_[e] * 2 * ma * mb / (ma + mb) * (e < connX_.size() ? connX_[e] : 1.0);
    }
}

void Simulator::setConnectionMultipliersByCells(const std::vector<std::pair<int, int>>& cells, const std::vector<double>& mult) {
    std::map<std::pair<int, int>, size_t> idx;
    for (size_t e = 0; e < conns_.size(); ++e) idx[{std::min(conns_[e].a, conns_[e].b), std::max(conns_[e].a, conns_[e].b)}] = e;
    std::vector<double> old = connX_;
    connX_.assign(conns_.size(), 1.0);
    for (size_t q = 0; q < cells.size() && q < mult.size(); ++q) {
        auto it = idx.find({std::min(cells[q].first, cells[q].second), std::max(cells[q].first, cells[q].second)});
        if (it != idx.end()) connX_[it->second] = std::isfinite(mult[q]) ? std::clamp(mult[q], 1e-8, 1e8) : 1.0;
    }
    const Case& c = c_;
    const bool law = c.permLaw != 0 || (geo_.ready() && ((c.geo.plasticity == 1 && c.geo.permPlastic > 0) || c.geo.permStress > 0));
    if (law) { updatePermeability(); return; }
    for (size_t e = 0; e < conns_.size(); ++e) {
        const double was = e < old.size() ? old[e] : 1.0;
        if (connX_[e] != 1.0 || was != 1.0) conns_[e].T = connT0_[e] * connX_[e];
    }
}

// C12: fracture apertures from the effective normal stress (Barton-Bandis), thermoelastic unloading and Coulomb shear
// dilation (irreversible, updated at step ends); cubic-law transmissibilities, implicit (linearised) fracture pore volume
void Simulator::updateFractures(bool stepEnd) {
    const Case& c = c_;
    std::vector<double> ratio(static_cast<size_t>(n_), 1.0);
    for (size_t f = 0; f < c.fracCells.size(); ++f) {
        Case::FracCell& F = c_.fracCells[f];
        const int q = F.cell;
        const double p = p_[size_t(q)], T = T_[size_t(q)];
        const double sEff = F.sn0 + F.thermStress * (T - F.T0) - p;
        double wbb = F.wOpen, dwds = 0;
        if (F.kn0 > 0 && F.vm > 0) {
            if (sEff > 0) {
                const double d = F.kn0 * F.vm + sEff;
                wbb = F.wOpen - F.vm * sEff / d;
                dwds = -F.vm * F.kn0 * F.vm / (d * d);
            } else { wbb = F.wOpen; dwds = 0; }
        }
        if (stepEnd && F.ks > 0) {
            const double excess = F.tau0 - F.mu * std::max(0.0, sEff);
            if (excess > 0) F.wShear = std::max(F.wShear, excess / F.ks * F.tanPsi);
        }
        double w = std::max(wbb + F.wShear, 1e-3 * F.w0);
        // explicit (step-lagged) aperture: at most 25 % change per time step, so that the pore-volume jump stays a small
        // perturbation of the next Newton solve (large thermo-elastic unloading follows over a few steps)
        if (stepEnd && f < fracW_.size() && fracW_[f] > 0) w = std::min(std::max(w, fracW_[f] / 1.25), fracW_[f] * 1.25);
        if (f < fracW_.size()) fracW_[f] = w;
        const double pv0 = pv0_[size_t(q)], r = w / F.w0;
        pvA_[size_t(q)] = pv0 * r; pvP_[size_t(q)] = p; pvT_[size_t(q)] = T;
        pvCp_[size_t(q)] = pv0 / F.w0 * (-dwds);                  // ds'/dp = -1
        pvCt_[size_t(q)] = pv0 / F.w0 * dwds * F.thermStress;     // ds'/dT = dsigma_n/dT
        ratio[size_t(q)] = r * r * r;
    }
    for (int e : c.fracFF) {
        const int a = conns_[size_t(e)].a, b = conns_[size_t(e)].b;
        double base = 1;
        if (c.permLaw != 0 || (geo_.ready() && ((c.geo.plasticity == 1 && c.geo.permPlastic > 0) || c.geo.permStress > 0))) { const double ma = permMult_[size_t(a)], mb = permMult_[size_t(b)]; base = 2 * ma * mb / (ma + mb); }
        const double ra = ratio[size_t(a)], rb = ratio[size_t(b)];
        conns_[size_t(e)].T = connT0_[size_t(e)] * base * 2 * ra * rb / (ra + rb);
    }
}

std::vector<double> Simulator::fractureApertures() const { return fracW_; }

bool Simulator::evalAll() {
    const int n = n_, nb = nb_, m = m_;
    bool ok = true;
    std::vector<Local> Ls(static_cast<size_t>(n));
#pragma omp parallel for schedule(dynamic, 16) if (useParallel(n > 64))
    for (int q = 0; q < n; ++q) {
        const double* X = &X_[size_t(q) * nb];
        double* Y = &Y_[size_t(q) * m];
        const double Tg = T_[size_t(q)];
        if (!evalCell(q, X, Y, &Ls[size_t(q)], Tg)) { ok = false; continue; }
        double Xp[kMaxBlock], Yp[256];
        std::vector<double> Yv;
        double* yp = Yp;
        if (m > 256) { Yv.resize(size_t(m)); yp = Yv.data(); }
        double Nt = 0; for (int k = 0; k <= nc_; ++k) Nt += std::max(X[1 + k], 0.0);
        const double Cest = crock_[size_t(q)] + Nt * 80.0;
        for (int j = 0; j < nb; ++j) {
            std::memcpy(Xp, X, sizeof(double) * size_t(nb));
            double h;
            if (j == 0) h = 1e-7 * std::max(X[0], 1e5);
            else if (j == nb - 1) h = 1e-4 * Cest;
            else h = 1e-7 * Nt + 1e-30;
            // one-sided differences that stay in the regime (phase set, saturation line) of the base state when possible
            Xp[j] += h;
            Local Lp;
            const int reg0 = Ls[size_t(q)].regime;
            bool fwd = evalCell(q, Xp, yp, &Lp, Ls[size_t(q)].T);
            if (!fwd || Lp.regime != reg0) {
                double Xb[kMaxBlock]; std::memcpy(Xb, X, sizeof(double) * size_t(nb)); Xb[j] = X[j] - h;
                Local Lb;
                double Ybk[256]; std::vector<double> Ybv; double* yb = Ybk;
                if (m > 256) { Ybv.resize(size_t(m)); yb = Ybv.data(); }
                const bool bwd = (j == 0 || Xb[j] >= 0 || j == nb - 1) && evalCell(q, Xb, yb, &Lb, Ls[size_t(q)].T);
                if (bwd && (Lb.regime == reg0 || !fwd)) { std::memcpy(yp, yb, sizeof(double) * size_t(m)); h = -h; }
                else if (!fwd) { ok = false; break; }
            }
            double* D = &dY_[size_t(q) * m * nb];
            for (int r = 0; r < m; ++r) D[r * nb + j] = (yp[r] - Y[r]) / h;
        }
    }
    if (!ok) return false;
    for (int q = 0; q < n; ++q) {
        const Local& L = Ls[size_t(q)];
        p_[size_t(q)] = X_[size_t(q) * nb];
        T_[size_t(q)] = L.T; sw_[size_t(q)] = L.S[W]; so_[size_t(q)] = L.S[O]; sg_[size_t(q)] = L.S[G];
    }
    st_.flashes += n * (nb + 1);
    return true;
}

// ------------------------------------------------------------------------------------------------ surface volumes
void Simulator::surfaceVolumes(const std::vector<double>& molar, double& vo, double& vw, double& vg) const {
    vo = vw = vg = 0;
    double tot = 0; for (double v : molar) tot += std::max(v, 0.0);
    if (!(tot > 0)) return;
    const FlashState f = flashPT(c_, c_.stdP, c_.stdT, molar);
    PhaseSet S; phaseProps(c_, c_.stdP, c_.stdT, f, S);
    if (S.ph[O].on) vo = S.ph[O].n * S.ph[O].v;
    if (S.ph[W].on) {
        const double aw = f.a.size() > size_t(nc_) ? f.a[size_t(nc_)] : 1.0;
        if (aw < 1.0) {   // gas dissolved in the stock-tank water is reported as surface gas (ideal at standard conditions)
            vw = S.ph[W].n * aw * steam::liquid(c_.stdT, c_.stdP).v * steam::kMw;
            vg += S.ph[W].n * (1 - aw) * kRgas * c_.stdT / c_.stdP;
        } else vw = S.ph[W].n * S.ph[W].v;
    }
    if (S.ph[G].on) {
        // steam in the stock-tank gas condenses: counted as water
        const double yw = f.y[size_t(nc_)];
        vg += S.ph[G].n * (1 - yw) * kRgas * c_.stdT / c_.stdP;
        vw += S.ph[G].n * yw * steam::liquid(c_.stdT, c_.stdP).v * steam::kMw;
    }
}

double Simulator::fluidMoles(int k) const {
    double s = 0; for (int q = 0; q < n_; ++q) s += X_[size_t(q) * nb_ + 1 + k]; return s;
}
double Simulator::energyInPlace() const {
    double s = 0; for (int q = 0; q < n_; ++q) s += X_[size_t(q) * nb_ + nb_ - 1]; return s;
}
double Simulator::stockTankOilInPlace() const {
    double s = 0;
    for (int q = 0; q < n_; ++q) {
        std::vector<double> z(static_cast<size_t>(nc_) + 1, 0.0);
        for (int k = 0; k < nc_; ++k) z[size_t(k)] = c_.comps[size_t(k)].solid() ? 0.0 : X_[size_t(q) * nb_ + 1 + k];
        double vo, vw, vg; surfaceVolumes(z, vo, vw, vg); s += vo;
    }
    return s;
}

// ------------------------------------------------------------------------------------------------ Newton
namespace {
int findBlock(const BSR& J, int r, int c) {
    for (int k = J.rowPtr[size_t(r)]; k < J.rowPtr[size_t(r) + 1]; ++k) if (J.col[size_t(k)] == c) return k;
    return -1;
}
}

bool Simulator::newtonStep(double dt, bool& converged, double& resNorm) {
    // (called with evalAll() done for the current iterate)
    const Case& c = c_;
    const int n = n_, nb = nb_, m = m_, nc = nc_, nk = nc + 1;
    const int eE = nk, eV = nk + 1;   // equation rows: components 0..nc, energy, volume
    const int vU = nb - 1;
    BSR& J = J_;
    if (J.n != n || J.nb != nb) {
        std::vector<std::set<int>> nbr(static_cast<size_t>(n));
        for (int q = 0; q < n; ++q) nbr[size_t(q)].insert(q);
        for (const auto& e : conns_) { nbr[size_t(e.a)].insert(e.b); nbr[size_t(e.b)].insert(e.a); }
        for (const Well& w : c.wells) for (int a : w.cells) for (int b : w.cells) nbr[size_t(a)].insert(b);
        J = BSR(); J.n = n; J.nb = nb; J.rowPtr.assign(size_t(n) + 1, 0);
        for (int q = 0; q < n; ++q) { J.rowPtr[size_t(q) + 1] = J.rowPtr[size_t(q)] + int(nbr[size_t(q)].size()); for (int b : nbr[size_t(q)]) J.col.push_back(b); }
        J.val.assign(size_t(J.rowPtr[size_t(n)]) * nb * nb, 0.0);
        J.finalizeDiag();
    }
    std::fill(J.val.begin(), J.val.end(), 0.0);
    std::vector<double> R(static_cast<size_t>(n) * nb, 0.0), acc(size_t(n) * nb * nb, 0.0);
    auto Yv = [&](int q, int r) { return Y_[size_t(q) * m + r]; };
    auto dYv = [&](int q, int r) { return &dY_[size_t(q) * m * nb + size_t(r) * nb]; };
    auto blk = [&](int r, int col) { return J.blk(findBlock(J, r, col)); };
    // accumulation + volume balance
    std::vector<double> Ntot(static_cast<size_t>(n)), Cc(static_cast<size_t>(n)), sE(static_cast<size_t>(n)), sV(static_cast<size_t>(n));
    for (int q = 0; q < n; ++q) {
        const double* X = &X_[size_t(q) * nb]; const double* Xo = &Xn_[size_t(q) * nb];
        double Nt = 0; for (int k = 0; k < nk; ++k) Nt += std::max(X[1 + k], 0.0);
        Ntot[size_t(q)] = Nt; Cc[size_t(q)] = crock_[size_t(q)] + Nt * 80.0;
        sE[size_t(q)] = Nt / Cc[size_t(q)];           // energy rows scaled to mol-equivalents of 1 K
        sV[size_t(q)] = Nt / std::max(pv0_[size_t(q)], 1e-30);   // volume rows scaled to moles
        double* Rq = &R[size_t(q) * nb];
        double* D = J.blk(J.diag[size_t(q)]);
        double* A = &acc[size_t(q) * nb * nb];
        for (int k = 0; k < nk; ++k) { Rq[k] = X[1 + k] - Xo[1 + k]; D[k * nb + 1 + k] += 1; A[k * nb + 1 + k] = 1; }
        Rq[eE] = X[vU] - Xo[vU]; D[eE * nb + vU] += 1; A[eE * nb + vU] = 1;
        Rq[eV] = Yv(q, 1);
        const double* dv = dYv(q, 1);
        for (int j = 0; j < nb; ++j) { D[eV * nb + j] += dv[j]; A[eV * nb + j] = dv[j]; }
    }
    // reactions: local sources (component balances and energy)
    if (rbase_ > 0)
        for (int q = 0; q < n; ++q) {
            double* Rq = &R[size_t(q) * nb];
            double* D = J.blk(J.diag[size_t(q)]);
            for (int r = 0; r <= nk; ++r) {
                const int row = r < nk ? r : eE;
                Rq[row] -= dt * Yv(q, rbase_ + r);
                const double* d = dYv(q, rbase_ + r);
                for (int j = 0; j < nb; ++j) D[row * nb + j] -= dt * d[j];
            }
        }
    // inter-cell flow: convection (phase-wise upwinding) + conduction
    const int base = 4, stride = stride_;
    for (const auto& e : conns_) {
        const int a = e.a, b = e.b;
        double Fa[kMaxBlock * kMaxBlock], Fb[kMaxBlock * kMaxBlock];   // d(flux_eq)/dX_a, dX_b (rows: nk comps + energy)
        double F[kMaxBlock];
        std::fill(F, F + nk + 1, 0.0);
        std::fill(Fa, Fa + (nk + 1) * nb, 0.0); std::fill(Fb, Fb + (nk + 1) * nb, 0.0);
        const double dD = depth_[size_t(a)] - depth_[size_t(b)];
        for (int ph = 0; ph < 3; ++ph) {
            const int r0 = base + ph * stride;
            const double ra = Yv(a, r0), rb = Yv(b, r0);
            const bool ona = ra > 0, onb = rb > 0;
            if (!ona && !onb) continue;
            double rho, wa, wb;
            if (ona && onb) { rho = 0.5 * (ra + rb); wa = wb = 0.5; } else if (ona) { rho = ra; wa = 1; wb = 0; } else { rho = rb; wa = 0; wb = 1; }
            const double pca = Yv(a, r0 + 1), pcb = Yv(b, r0 + 1);
            const double dPhi = (X_[size_t(a) * nb] + pca - X_[size_t(b) * nb] - pcb) - rho * kG * dD;
            const int up = dPhi >= 0 ? a : b;
            // d dPhi / dX
            double dPa[kMaxBlock], dPb[kMaxBlock];
            const double* dra = dYv(a, r0); const double* drb = dYv(b, r0);
            for (int j = 0; j < nb; ++j) { dPa[j] = -wa * dra[j] * kG * dD; dPb[j] = -wb * drb[j] * kG * dD; }
            dPa[0] += 1; dPb[0] -= 1;
            if (hasPc_) { const double* dpa = dYv(a, r0 + 1); const double* dpb = dYv(b, r0 + 1); for (int j = 0; j < nb; ++j) { dPa[j] += dpa[j]; dPb[j] -= dpb[j]; } }
            // energy (row nk) and components (rows 0..nc)
            for (int r = 0; r <= nk; ++r) {
                const int yr = r < nk ? r0 + 3 + r : r0 + 2;
                const double mob = Yv(up, yr);
                const double* dm = dYv(up, yr);
                const double f = e.T * mob * dPhi;
                F[r] += f;
                for (int j = 0; j < nb; ++j) {
                    Fa[r * nb + j] += e.T * mob * dPa[j];
                    Fb[r * nb + j] += e.T * mob * dPb[j];
                }
                double* Fu = up == a ? Fa : Fb;
                for (int j = 0; j < nb; ++j) Fu[r * nb + j] += e.T * dm[j] * dPhi;
            }
        }
        // conduction
        {
            // harmonic averaging of the cell conductivities (series resistances of the two half cells; the arithmetic
            // mean over-estimates the heat flow across conductivity contrasts - shale / sand, burnt cavity / coal)
            const double la = Yv(a, 2), lb = Yv(b, 2), ls = la + lb;
            const double lam = ls > 0 ? 2 * la * lb / ls : 0.0;
            const double gla = ls > 0 ? 2 * lb * lb / (ls * ls) : 0.0, glb = ls > 0 ? 2 * la * la / (ls * ls) : 0.0;
            const double Ta = Yv(a, 0), Tb = Yv(b, 0);
            const double f = e.G * lam * (Ta - Tb);
            F[nk] += f;
            const double* dla = dYv(a, 2); const double* dlb = dYv(b, 2);
            const double* dTa = dYv(a, 0); const double* dTb = dYv(b, 0);
            for (int j = 0; j < nb; ++j) {
                Fa[nk * nb + j] += e.G * (gla * dla[j] * (Ta - Tb) + lam * dTa[j]);
                Fb[nk * nb + j] += e.G * (glb * dlb[j] * (Ta - Tb) - lam * dTb[j]);
            }
        }
        double* Raa = blk(a, a); double* Rab = blk(a, b); double* Rba = blk(b, a); double* Rbb = blk(b, b);
        for (int r = 0; r <= nk; ++r) {
            R[size_t(a) * nb + r] += dt * F[r];
            R[size_t(b) * nb + r] -= dt * F[r];
            for (int j = 0; j < nb; ++j) {
                Raa[r * nb + j] += dt * Fa[r * nb + j]; Rab[r * nb + j] += dt * Fb[r * nb + j];
                Rba[r * nb + j] -= dt * Fa[r * nb + j]; Rbb[r * nb + j] -= dt * Fb[r * nb + j];
            }
        }
    }
    // heat loss (Vinsome-Westerveld): Q = a theta + b, theta = T - Ti
    if (!loss_.empty()) {
        const double kap = c.obK / c.obCv;
        const double tEnd = t_ + dt;
        const double d = 0.5 * std::sqrt(kap * tEnd), kdt = kap * dt, Dd = 3 * d * d + kdt;
        for (size_t f = 0; f < loss_.size(); ++f) {
            const int q = c.lossFaces[f].cell;
            const double kA = c.obK * c.lossFaces[f].area;
            const double aC = kA * (1 / d - (kdt / d - d * d * d / kdt) / Dd);
            const double bC = -kA * (loss_[f].I + d * d * d * loss_[f].thetaN / kdt) / Dd;
            const double th = Yv(q, 0) - loss_[f].Ti;
            R[size_t(q) * nb + eE] += dt * (aC * th + bC);
            const double* dT = dYv(q, 0);
            double* D = J.blk(J.diag[size_t(q)]);
            for (int j = 0; j < nb; ++j) D[eE * nb + j] += dt * aC * dT[j];
        }
    }
    // wells
    for (size_t wi = 0; wi < c.wells.size(); ++wi) {
        const Well& w = c.wells[wi];
        WellReport& rp = rep_[wi];
        std::fill(rp.molar.begin(), rp.molar.end(), 0.0); rp.hRate = 0; rp.tFluid = 0;
        double tfNum = 0, tfDen = 0;
        if (!w.open || w.cells.empty()) continue;
        const size_t ncon = w.cells.size();
        const double dref = w.refDepth >= 0 ? w.refDepth : depth_[size_t(w.cells[0])];
        const double rwb = wbRho_[wi] > 0 ? wbRho_[wi] : 800.0;
        std::vector<double> dh(ncon);
        for (size_t j = 0; j < ncon; ++j) dh[j] = rwb * kG * (depth_[size_t(w.cells[j])] - dref);
        // per connection: total coefficient A_j (rate measure per unit drawdown) and its derivative
        std::vector<double> Aj(ncon, 0.0), dAj(ncon * nb, 0.0);
        double hinj = 0, xinj = 1, Tinj = 0;
        if (w.injector) {
            // injected enthalpy (J/mol) and molar density (mol/m3) at the wellbore pressure (lagged BHP)
            const double pinj = rp.bhp > 0 ? rp.bhp : X_[size_t(w.cells[0]) * nb] + 1e6;
            const double zw = w.zinj[size_t(nc)];
            double hw = 0, vw = 0, hh = 0, vh = 0, mwInj = zw * steam::kMw;
            // E300 WINJTEMP: the enthalpy is fixed by (quality, T) or (quality, P) or (P, T) or H; the molar density for the
            // inflow uses the wellbore pressure
            double Ti = w.tinj > 0 ? w.tinj : 300.0, Pi = w.pinj > 0 ? w.pinj : pinj;
            if (w.quality >= 0) {
                if (w.tinj > 0) { Ti = std::min(w.tinj, steam::kTc - 1e-3); Pi = steam::psat(Ti); }
                else Ti = steam::tsat(std::min(Pi, 0.999 * steam::kPc));
            }
            if (zw > 0) {
                if (w.quality >= 0) {
                    const double Ps = steam::psat(Ti);
                    const steam::Props L = steam::liquid(Ti, Ps), V = steam::vapour(Ti, Ps);
                    hw = steam::kMw * ((1 - w.quality) * L.h + w.quality * V.h);
                    const steam::Props Lw = steam::liquid(Ti, std::max(pinj, Ps)), Vw = steam::vapour(Ti, std::min(pinj, Ps));
                    vw = steam::kMw * ((1 - w.quality) * Lw.v + w.quality * Vw.v);
                } else {
                    const steam::Props P = Ti < steam::tsat(std::min(Pi, 0.999 * steam::kPc)) ? steam::liquid(Ti, Pi) : steam::vapour(Ti, Pi);
                    hw = steam::kMw * P.h; vw = steam::kMw * P.v;
                }
            }
            for (int i = 0; i < nc; ++i) {
                const double zi = w.zinj[size_t(i)];
                if (zi <= 0) continue;
                const Component& qc = c.comps[size_t(i)];
                mwInj += zi * qc.mw;
                if (qc.volatileK()) { hh += zi * hGas(c, qc, Ti, Pi); vh += zi * gasMolarVolume(qc, pinj, Ti); }
                else { hh += zi * hLiq(c, qc, Ti); vh += zi * liqMolarVolume(c, qc, pinj, Ti); }
            }
            if (const GasEos* ge = c.gasPR ? gasEos(c) : nullptr) {
                // C16 fix: Peng-Robinson subset of the injected gas - its molar volume and residual enthalpy (consistent with
                // the in-place energy of phaseProps; before, dense cold CO2 entered with ideal-gas enthalpy, i.e. too warm)
                std::vector<double> yg(static_cast<size_t>(nc), 0.0);
                double vOld = 0;
                for (int i = 0; i < nc; ++i) {
                    const Component& qc = c.comps[size_t(i)];
                    if (ge->idx[size_t(i)] < 0 || !qc.volatileK() || w.zinj[size_t(i)] <= 0) continue;
                    yg[size_t(i)] = w.zinj[size_t(i)]; vOld += w.zinj[size_t(i)] * gasMolarVolume(qc, pinj, Ti);
                }
                GasPR gp;
                if (vOld > 0 && gasPRState(c, *ge, pinj, Ti, yg, gp, true) && gp.s > 0) { vh += gp.s * gp.vm - vOld; hh += gp.s * gp.hres; }
            }
            hinj = w.hinj > 0 ? w.hinj * mwInj : zw * hw + hh;
            Tinj = Ti;
            xinj = 1.0 / std::max(zw * vw + vh, 1e-12);   // mol/m3
            for (size_t j = 0; j < ncon; ++j) {
                const int q = w.cells[j];
                const double wij = w.wi[j] * permMult_[size_t(q)];
                Aj[j] = wij * Yv(q, 3) * xinj;
                const double* d = dYv(q, 3);
                for (int k = 0; k < nb; ++k) dAj[j * nb + k] = wij * d[k] * xinj;
            }
        } else {
            // weights of the controlled surface quantity per mole of component (lagged surface yields)
            std::vector<double> s(static_cast<size_t>(nk), 0.0);
            if (w.control > 0) {
                if (surfYield_[wi].empty()) {
                    std::vector<double> z(&X_[size_t(w.cells[0]) * nb + 1], &X_[size_t(w.cells[0]) * nb + 1 + nk]);
                    std::vector<double> yv(static_cast<size_t>(2 * nk), 0.0);
                    double vo, vww, vg; surfaceVolumes(z, vo, vww, vg);
                    double zh = 0; for (int k = 0; k < nc; ++k) zh += z[size_t(k)];
                    for (int k = 0; k < nc; ++k) yv[size_t(k)] = zh > 0 ? vo / zh : 0;
                    yv[size_t(nk + nc)] = steam::liquid(c.stdT, c.stdP).v * steam::kMw;
                    surfYield_[wi] = yv;
                }
                const auto& yv = surfYield_[wi];
                for (int k = 0; k < nk; ++k) {
                    const double o = yv[size_t(k)], wv = yv[size_t(nk + k)];
                    s[size_t(k)] = w.control == 1 ? o : w.control == 2 ? wv : o + wv;
                }
            }
            for (size_t j = 0; j < ncon; ++j) {
                const int q = w.cells[j];
                double a = 0; std::vector<double> da(static_cast<size_t>(nb), 0.0);
                for (int ph = 0; ph < 3; ++ph)
                    for (int k = 0; k < nk; ++k) {
                        const int yr = base + ph * stride + 3 + k;
                        const double sk = w.control > 0 ? s[size_t(k)] : 1.0;
                        a += sk * Yv(q, yr);
                        const double* d = dYv(q, yr);
                        for (int l = 0; l < nb; ++l) da[size_t(l)] += sk * d[l];
                    }
                Aj[j] = w.wi[j] * permMult_[size_t(q)] * a;
                for (int l = 0; l < nb; ++l) dAj[j * nb + l] = w.wi[j] * permMult_[size_t(q)] * da[size_t(l)];
            }
        }
        // BHP: fixed, or eliminated from the rate target with the active-connection set
        const double sign = w.injector ? 1.0 : -1.0;   // injector flow ~ (bhp - p), producer ~ (p - bhp)
        bool rateCtl = w.injector ? !w.bhpControl : w.control > 0;
        double target = w.injector ? w.rate : w.target;
        double bhpLim = w.injector ? (w.bhpControl ? w.bhpMax : (w.bhpMax > 0 ? w.bhpMax : 1e12)) : w.bhp;
        rp.steamLimited = false;
        if (!w.injector && w.steamMax > 0 && ncon > 0) {
            // C17 steam-rate limit: rate measure of the produced steam (CWE m3 per unit drawdown) per connection
            std::vector<double> As(ncon, 0.0), dAs(ncon * nb, 0.0);
            const int yr = base + 2 * stride + 3 + nc;   // vapour-phase water molar mobility
            const double cwe = steam::kMw / 1000.0;
            double sa = 0, sp = 0, sa0 = 0, sp0 = 0;
            for (size_t j = 0; j < ncon; ++j) {
                const int q = w.cells[j];
                const double wij = w.wi[j] * permMult_[size_t(q)];
                As[j] = wij * cwe * Yv(q, yr);
                const double* d = dYv(q, yr);
                for (int k = 0; k < nb; ++k) dAs[j * nb + k] = wij * cwe * d[k];
                const double pj = X_[size_t(q) * nb] - dh[j];
                sa += As[j]; sp += As[j] * pj; sa0 += Aj[j]; sp0 += Aj[j] * pj;
            }
            if (sa > 0) {
                const double bS = (sp - w.steamMax) / sa;                      // BHP that produces exactly steamMax
                const double bP = rateCtl && sa0 > 0 ? (sp0 - target) / sa0 : bhpLim;   // BHP of the primary control
                const double bPrim = std::max(bP, bhpLim);
                if (steamMode_.size() != c.wells.size()) { steamMode_.assign(c.wells.size(), -1); steamSw_.assign(c.wells.size(), 0); }
                int& sm = steamMode_[wi]; int& ss = steamSw_[wi];
                const double hyst = sm == 1 ? -2e3 : 2e3;                      // 0.02 bar band against flip-flop
                int want = bS > bPrim + hyst ? 1 : 0;
                if (sm >= 0 && want != sm && ss >= 2) want = 1;                // after two switches stay on the conservative side
                if (!assembleOnly_) { if (sm >= 0 && want != sm) ++ss; sm = want; }
                if (want == 1) {
                    rateCtl = true; target = w.steamMax; rp.steamLimited = true;
                    Aj = As; dAj = dAs;
                }
            }
        }
        // steam trap: dynamic BHP floor Psat(T_hot + subcool) with its derivative d Psat / dT (central difference)
        int trapJ = -1; double dTrap = 0;
        if (!w.injector && w.subcool > 0 && ncon > 0) {
            size_t jc = 0; double Tm = -1;
            for (size_t j = 0; j < ncon; ++j) { const double Tj = Yv(w.cells[j], 0); if (Tj > Tm) { Tm = Tj; jc = j; } }
            const double Tt = Tm + w.subcool, Tcap = steam::kTc - 0.5;
            const double pt = steam::psat(std::min(Tt, Tcap));
            if (pt > bhpLim) {
                bhpLim = pt; trapJ = int(jc);
                if (Tt < Tcap) { const double e = 0.01; dTrap = (steam::psat(Tt + e) - steam::psat(Tt - e)) / (2 * e); }
            }
        }
        rp.trapActive = false;
        // injector connections: no back-flow (a connection must not produce the injection stream - that removed injected
        // components the cell does not hold and stalled Newton for injectors pinned at their BHP limit), flow through the C1
        // drawdown ramp q = A g(x): g = 0 (x <= 0), x^2 / 2d (0 < x < d), x - d/2 (x >= d), in rate and BHP mode alike (the
        // same law on both sides of a rate / BHP switch). Width d = 5% of the drawdown carrying the target, lagged per
        // attempt, capped at kConRamp. Producers keep the linear law on the active-connection set.
        const bool injRamp = w.injector;
        double rampD = kConRamp;
        if (injRamp) {
            if (wellRamp_.size() < c.wells.size()) wellRamp_.assign(c.wells.size(), -1.0);
            double& rd = wellRamp_[wi];
            if (!assembleOnly_ && (newtonIt_ == 0 || rd <= 0)) {
                double sa = 0; for (size_t j = 0; j < ncon; ++j) sa += Aj[j];
                rd = (target > 0 && sa > 0) ? std::clamp(0.05 * target / sa, 1e-6, kConRamp) : std::min(kConRamp, 1e-4 * std::max(bhpLim, 1e5));
            }
            if (rd > 0) rampD = rd;
        }
        auto ramp = [&](double x, double& gp) {
            if (!injRamp) { gp = 1.0; return x; }
            if (x <= 0) { gp = 0.0; return 0.0; }
            if (x < rampD) { gp = x / rampD; return 0.5 * x * x / rampD; }
            gp = 1.0; return x - 0.5 * rampD;
        };
        double bhp = bhpLim;
        std::vector<char> act(ncon, 1);
        // connection on / off with a small hysteresis band (0.02 bar): a connection that flowed in the previous iteration
        // only shuts when its drawdown reverses clearly - stops on/off ping-pong of near-balanced connections in Newton
        if (wellAct_.size() < c.wells.size()) wellAct_.resize(c.wells.size());
        std::vector<char>& prevAct = wellAct_[wi];
        if (prevAct.size() != ncon) prevAct.assign(ncon, 1);
        auto flowing = [&](size_t j, double dpj) { return dpj > 0 || (prevAct[j] && dpj > -2e3); };
        act = prevAct;
        bool onBhp = !rateCtl;
        if (rateCtl && injRamp) {
            // injector rate target: sum_j A_j g(bhp - p_j + dh_j) = q (monotone in bhp) - safeguarded Newton in a bracket
            std::fill(act.begin(), act.end(), char(1));
            double sa = 0, sp = 0, pmin = 1e300, pmax = -1e300;
            for (size_t j = 0; j < ncon; ++j) {
                const double pj = X_[size_t(w.cells[j]) * nb] - dh[j];
                sa += Aj[j]; sp += Aj[j] * pj; pmin = std::min(pmin, pj); pmax = std::max(pmax, pj);
            }
            if (!(sa > 0) || !(target > 0)) { bhp = bhpLim; onBhp = true; }
            else {
                double lo = pmin, hi = pmax + target / sa + rampD;
                double b = std::clamp((target + sp) / sa + 0.5 * rampD, lo, hi);
                for (int it = 0; it < 100; ++it) {
                    double F = -target, dF = 0;
                    for (size_t j = 0; j < ncon; ++j) { double gp; F += Aj[j] * ramp(b - (X_[size_t(w.cells[j]) * nb] - dh[j]), gp); dF += Aj[j] * gp; }
                    if (F > 0) hi = b; else lo = b;
                    double bn = dF > 0 ? b - F / dF : 0.5 * (lo + hi);
                    if (!(bn > lo && bn < hi)) bn = 0.5 * (lo + hi);
                    const bool done = std::fabs(bn - b) <= 1e-13 * std::max(1.0, std::fabs(b)) || hi - lo <= 1e-13 * std::max(1.0, std::fabs(hi));
                    b = bn;
                    if (done) break;
                }
                bhp = b;
                if (bhp > bhpLim) { onBhp = true; bhp = bhpLim; }
            }
        } else if (rateCtl) {
            for (int pass = 0; pass < 8; ++pass) {
                double sa = 0, sp = 0;
                for (size_t j = 0; j < ncon; ++j) if (act[j]) { sa += Aj[j]; sp += Aj[j] * (X_[size_t(w.cells[j]) * nb] - dh[j]); }
                if (!(sa > 0)) { bhp = bhpLim; onBhp = true; break; }
                bhp = w.injector ? (target + sp) / sa : (sp - target) / sa;
                bool changed = false;
                for (size_t j = 0; j < ncon; ++j) {
                    const double dpj = sign * (bhp - (X_[size_t(w.cells[j]) * nb] - dh[j]));
                    const char nw = flowing(j, dpj) ? 1 : 0;
                    if (nw != act[j]) { act[j] = nw; changed = true; }
                }
                if (!changed) break;
            }
            if (!onBhp && ((w.injector && bhp > bhpLim) || (!w.injector && bhp < bhpLim))) { onBhp = true; bhp = bhpLim; }
        }
        // rate <-> BHP switching inside the Newton loop: after two switches in a time step the well stays on its BHP limit
        // (no rate / BHP flip-flop between iterations; the conservative side of the constraint)
        if (rateCtl && wellMode_.size() == c.wells.size()) {
            int want = onBhp ? 1 : 0;
            int& md = wellMode_[wi]; int& ns = wellSw_[wi];
            if (md >= 0 && want != md && ns >= 2) want = 1;
            if (!assembleOnly_) { if (md >= 0 && want != md) ++ns; md = want; }
            if (want == 1 && !onBhp) { onBhp = true; bhp = bhpLim; }
        }
        if (onBhp && !injRamp) for (size_t j = 0; j < ncon; ++j) act[j] = flowing(j, sign * (bhp - (X_[size_t(w.cells[j]) * nb] - dh[j]))) ? 1 : 0;
        if (injRamp) std::fill(act.begin(), act.end(), char(1));
        if (!assembleOnly_) prevAct = act;
        // d bhp / dX_l (rate control)
        std::vector<double> dB(ncon * nb, 0.0);
        const bool trapOn = onBhp && trapJ >= 0;   // BHP held at the steam-trap floor: bhp = Psat(T(X_hot) + subcool)
        if (trapOn) {
            const double* dT = dYv(w.cells[size_t(trapJ)], 0);
            for (int k = 0; k < nb; ++k) dB[size_t(trapJ) * nb + k] = dTrap * dT[k];
            rp.trapActive = true;
        }
        const bool bhpDep = !onBhp || trapOn;
        if (!onBhp && injRamp) {   // implicit function F(bhp, X) = 0: d bhp / dX = -F_X / F_bhp
            double sg = 0; for (size_t j = 0; j < ncon; ++j) { double gp; ramp(bhp - (X_[size_t(w.cells[j]) * nb] - dh[j]), gp); sg += Aj[j] * gp; }
            if (sg > 0)
                for (size_t l = 0; l < ncon; ++l) {
                    double gp; const double gx = ramp(bhp - (X_[size_t(w.cells[l]) * nb] - dh[l]), gp);
                    for (int k = 0; k < nb; ++k) dB[l * nb + k] = ((k == 0 ? Aj[l] * gp : 0.0) - dAj[l * nb + k] * gx) / sg;
                }
        } else if (!onBhp) {
            double sa = 0; for (size_t j = 0; j < ncon; ++j) if (act[j]) sa += Aj[j];
            for (size_t l = 0; l < ncon; ++l) {
                if (!act[l]) continue;
                const double dpl = X_[size_t(w.cells[l]) * nb] - dh[l] - bhp;
                for (int k = 0; k < nb; ++k) dB[l * nb + k] = (dAj[l * nb + k] * dpl + (k == 0 ? Aj[l] : 0.0)) / sa;
            }
        }
        rp.bhp = bhp; rp.onBhp = onBhp;
        if (onBhp) rp.steamLimited = false;
        // connection flows
        for (size_t j = 0; j < ncon; ++j) {
            if (!act[j]) continue;
            const int q = w.cells[j];
            const double dp = X_[size_t(q) * nb] - dh[j] - bhp;   // > 0 producer, < 0 injector
            double* Rq = &R[size_t(q) * nb];
            if (w.injector) {
                // q_j = A_j g(bhp - p + dh) (mol/s into the cell), A_j = WI lambda_t xi_inj
                double gp; const double gx = ramp(-dp, gp);
                const double qj = Aj[j] * gx;
                for (int k = 0; k < nk; ++k) Rq[k] -= dt * w.zinj[size_t(k)] * qj;
                Rq[eE] -= dt * hinj * qj;
                for (int k = 0; k < nk; ++k) rp.molar[size_t(k)] -= w.zinj[size_t(k)] * qj;
                rp.hRate -= hinj * qj;
                tfNum += qj * Tinj; tfDen += qj;
                // derivatives: dq_j/dX_q (own) and dq_j/dX_l via bhp
                double dq[kMaxBlock];
                for (int k = 0; k < nb; ++k) dq[k] = dAj[j * nb + k] * gx - (k == 0 ? Aj[j] * gp : 0.0);
                double* D = J.blk(J.diag[size_t(q)]);
                for (int r = 0; r < nk; ++r) for (int k = 0; k < nb; ++k) D[r * nb + k] -= dt * w.zinj[size_t(r)] * dq[k];
                for (int k = 0; k < nb; ++k) D[eE * nb + k] -= dt * hinj * dq[k];
                if (bhpDep)
                    for (size_t l = 0; l < ncon; ++l) {
                        if (!act[l] && int(l) != trapJ) continue;
                        double* B = blk(q, w.cells[l]);
                        for (int k = 0; k < nb; ++k) {
                            const double g = Aj[j] * gp * dB[l * nb + k];
                            for (int r = 0; r < nk; ++r) B[r * nb + k] -= dt * w.zinj[size_t(r)] * g;
                            B[eE * nb + k] -= dt * hinj * g;
                        }
                    }
            } else {
                // producer: q_jk = WI c_k(q) g(dp), energy WI e(q) g(dp)
                double gp; const double gx = ramp(dp, gp);
                for (int r = 0; r <= nk; ++r) {
                    double mob = 0; double dm[kMaxBlock]; std::fill(dm, dm + nb, 0.0);
                    for (int ph = 0; ph < 3; ++ph) {
                        const int yr = r < nk ? base + ph * stride + 3 + r : base + ph * stride + 2;
                        mob += Yv(q, yr);
                        const double* d = dYv(q, yr);
                        for (int k = 0; k < nb; ++k) dm[k] += d[k];
                    }
                    const double wij = w.wi[j] * permMult_[size_t(q)];
                    const double f = wij * mob * gx;
                    const int row = r < nk ? r : eE;
                    Rq[row] += dt * f;
                    if (r < nk) { rp.molar[size_t(r)] += f; tfNum += f * T_[size_t(q)]; tfDen += f; } else rp.hRate += f;
                    double* D = J.blk(J.diag[size_t(q)]);
                    for (int k = 0; k < nb; ++k) D[row * nb + k] += dt * wij * (dm[k] * gx + (k == 0 ? mob * gp : 0.0));
                    if (bhpDep)
                        for (size_t l = 0; l < ncon; ++l) {
                            if (!act[l] && int(l) != trapJ) continue;
                            double* B = blk(q, w.cells[l]);
                            for (int k = 0; k < nb; ++k) B[row * nb + k] -= dt * wij * mob * gp * dB[l * nb + k];
                        }
                }
            }
        }
        rp.tFluid = tfDen > 0 ? tfNum / tfDen : 0.0;
    }
    if (assembleOnly_) { lastR_ = R; lastJ_ = J.val; return true; }
    // scaling of the energy / volume rows, convergence measures
    double mMax = 0, eMax = 0, vMax = 0, eArg2 = 0;
    int vArg = -1, eArg = -1;
    for (int q = 0; q < n; ++q) {
        double* Rq = &R[size_t(q) * nb];
        for (int k = 0; k < nk; ++k) { mMax = std::max(mMax, std::fabs(Rq[k]) / std::max(Ntot[size_t(q)], 1e-30)); }
        eMax = std::max(eMax, std::fabs(Rq[eE]) / Cc[size_t(q)]);
        { const double vv = std::fabs(Rq[eV]) / std::max(pv0_[size_t(q)], 1e-30); if (vv > vMax) { vMax = vv; vArg = q; } }
        { const double ee = std::fabs(Rq[eE]) / Cc[size_t(q)]; if (ee > eArg2) { eArg2 = ee; eArg = q; } }
        Rq[eE] *= sE[size_t(q)]; Rq[eV] *= sV[size_t(q)];
        for (int k = J.rowPtr[size_t(q)]; k < J.rowPtr[size_t(q) + 1]; ++k) {
            double* B = J.blk(k);
            for (int j = 0; j < nb; ++j) { B[eE * nb + j] *= sE[size_t(q)]; B[eV * nb + j] *= sV[size_t(q)]; }
        }
        double* A = &acc[size_t(q) * nb * nb];
        for (int j = 0; j < nb; ++j) { A[eE * nb + j] *= sE[size_t(q)]; A[eV * nb + j] *= sV[size_t(q)]; }
    }
    resNorm = std::max({mMax / 1e-6, eMax / 1e-3, vMax / 1e-6});
    converged = mMax < tolM_ && eMax < tolE_ && vMax < tolV_;
    // C12: at least one Newton update per step unless the start state is already converged to a much tighter level (the
    // loose tolerances let slowly diffusing cells skip whole steps: with small time steps the per-step change of a far
    // cell stays below 1e-6 N and was never applied - the pressure / heat front then lagged, a time-step-size bias)
    if (newtonIt_ == 0) converged = mMax < 1e-10 && eMax < 1e-7 && vMax < 1e-10;
    if (!converged && newtonIt_ >= 6 && tolM_ >= 1e-6) {
        // relaxed local (CNV-type) tolerances on the late iterations, provided the global balances are tight (OPM practice:
        // a few cells oscillating on a phase-state boundary must not force a time-step cut)
        double sM = 0, sN = 0, sE2 = 0, sC = 0;
        for (int k = 0; k < nk; ++k) { double a = 0; for (int q = 0; q < n; ++q) a += R[size_t(q) * nb + k]; sM = std::max(sM, std::fabs(a)); }
        for (int q = 0; q < n; ++q) { sN += Ntot[size_t(q)]; sE2 += R[size_t(q) * nb + eE] / sE[size_t(q)]; sC += Cc[size_t(q)]; }
        const double mb = sM / std::max(sN, 1e-30), eb = std::fabs(sE2) / std::max(sC, 1e-30);
        converged = mMax < 1e-5 && eMax < 1e-2 && vMax < 2e-5 && mb < 1e-9 && eb < 1e-6;
    }
    resHist_.push_back(resNorm);
    {   // oscillation / stagnation of the residual: relax the update (omega 1 -> 0.5 -> 0.3)
        const size_t h = resHist_.size();
        if (h >= 3) {
            const double r0 = resHist_[h - 1], r2 = resHist_[h - 3], r1 = resHist_[h - 2];
            const bool osc = std::fabs(r0 - r2) < 0.2 * r0 && std::fabs(r1 - r0) > 0.05 * r0 && r0 > 0.5 * r1;   // 2-cycle
            const bool stag = h >= 4 && r0 > 0.85 * std::min(r1, r2) && r1 > 0.85 * r2;   // no progress over two iterations
            if (osc || stag) relax_ = std::max(0.3, relax_ * 0.6);
            // converging again: back towards full steps. A damped Newton step can at best reduce the residual by (1 - omega),
            // so the recovery test scales with omega (a fixed 0.5 would lock omega = 0.3 in: linear convergence -> cuts)
            else if (r0 < (1.0 - 0.5 * relax_) * r1) relax_ = std::min(1.0, relax_ / 0.6);
        }
    }
    static const bool dbg = std::getenv("THERM_DEBUG") != nullptr;
    if (dbg) std::fprintf(stderr, "    newton dt %.4g d: mass %.2e energy %.2e (cell %d T %.3f sw %.4f sg %.4f) vol %.2e (cell %d T %.3f sw %.4f sg %.4f p %.6f)\n", dt / 86400, mMax, eMax,
                          eArg, eArg >= 0 ? T_[size_t(eArg)] : 0.0, eArg >= 0 ? sw_[size_t(eArg)] : 0.0, eArg >= 0 ? sg_[size_t(eArg)] : 0.0, vMax,
                          vArg, vArg >= 0 ? T_[size_t(vArg)] : 0.0, vArg >= 0 ? sw_[size_t(vArg)] : 0.0, vArg >= 0 ? sg_[size_t(vArg)] : 0.0, vArg >= 0 ? p_[size_t(vArg)] / 1e6 : 0.0);
    if (converged) return true;
    // solve J dx = -R
    std::vector<double> rhs(R.size()), dx(R.size(), 0.0);
    for (size_t i = 0; i < R.size(); ++i) rhs[i] = -R[i];
    LinResult lr;
    bool ddDone = false;
    if (distributedSolver()) {   // C20: MPI domain-decomposed solve (ResSimCLI --mpi-dd); pressure is the first unknown of each block
        int itsD = 0; std::vector<double> xd;
        if (distributedSolver()->solve(J, rhs, xd, 1e-6, 400, itsD) && xd.size() == dx.size()) {
            dx.swap(xd); lr.converged = true; lr.iterations = itsD; lr.relResidual = 1e-6; ddDone = true;
        }
    }
    if (!ddDone) {
        BlockCPR cpr;
        cpr.setup(J, acc.data(), false);
        BsrOperator op(J);
        lr = gmres(op, cpr, rhs, dx, 1e-6, 400, 80);
    }
    st_.linear += lr.iterations;
    if (dbg) std::fprintf(stderr, "      linear its %d converged %d\n", lr.iterations, int(lr.converged));
    for (double v : dx) if (!std::isfinite(v)) return false;
    // damped update, cell by cell (local chopping as in commercial FIM codes - a single difficult cell, e.g. the steam
    // injection block, no longer throttles the whole field's Newton update): |dp| <= 0.3 p, |dU| <= 50 K equivalent,
    // moles >= 0; the global relaxation factor (oscillation / stagnation) multiplies the local factors
    const std::vector<double> Xold = X_;
    for (int q = 0; q < n; ++q) {
        double om = 1.0;
        {
            const double p = X_[size_t(q) * nb], dp = dx[size_t(q) * nb];
            if (std::fabs(dp) > 0.3 * p) om = std::min(om, 0.3 * p / std::fabs(dp));
            const double dU = dx[size_t(q) * nb + vU];
            if (std::fabs(dU) > 50.0 * Cc[size_t(q)]) om = std::min(om, 50.0 * Cc[size_t(q)] / std::fabs(dU));
        }
        om *= relax_;
        double* X = &X_[size_t(q) * nb];
        X[0] += om * dx[size_t(q) * nb];
        for (int k = 0; k < nk; ++k) X[1 + k] = std::max(0.0, X[1 + k] + om * dx[size_t(q) * nb + 1 + k]);
        X[vU] += om * dx[size_t(q) * nb + vU];
    }
    // Appleyard-type local chopping in (p, N, U): a cell whose saturations would jump by more than 0.25 in one update
    // (steam appearing / collapsing on the energy update) takes a shortened update; checked with a property-only flash
#pragma omp parallel for schedule(dynamic, 16) if (useParallel(n > 64))
    for (int q = 0; q < n; ++q) {
        double* X = &X_[size_t(q) * nb];
        const double* X0 = &Xold[size_t(q) * nb];
        double yb[256]; std::vector<double> yv; double* y = yb;
        if (m_ > 256) { yv.resize(size_t(m_)); y = yv.data(); }
        for (int cut = 0; cut < 3; ++cut) {
            Local L;
            if (!evalCell(q, X, y, &L, T_[size_t(q)])) { for (int j = 0; j < nb; ++j) X[j] = X0[j] + 0.5 * (X[j] - X0[j]); continue; }
            const double dS = std::max(std::fabs(L.S[0] - sw_[size_t(q)]), std::fabs(L.S[2] - sg_[size_t(q)]));
            // C12: temperature swings as well (cells without rock heat capacity, e.g. EDFM fracture cells, where the
            // (N, U) update of a cold front overshoots T by hundreds of K)
            const double dT = std::fabs(L.T - T_[size_t(q)]);
            if (dS <= 0.25 && dT <= 80.0) break;
            const double f = std::max(0.3, std::min(dS > 0.25 ? 0.25 / dS : 1.0, dT > 80.0 ? 80.0 / dT : 1.0));
            for (int j = 0; j < nb; ++j) X[j] = X0[j] + f * (X[j] - X0[j]);
        }
    }
    return true;
}

void Simulator::finishStep(double dt) {
    const Case& c = c_;
    // gas hysteresis: historical maximum gas saturation (lagged into the next step, as ECLIPSE / the compositional core)
    if (c.hyst.on) for (int q = 0; q < n_ && size_t(q) < sgHist_.size(); ++q) sgHist_[size_t(q)] = std::max(sgHist_[size_t(q)], sg_[size_t(q)]);
    // reaction heat
    reactRate_ = 0;
    if (rbase_ > 0) for (int q = 0; q < n_; ++q) reactRate_ += Y_[size_t(q) * m_ + rbase_ + nc_ + 1];
    cumReact_ += reactRate_ * dt;
    // heat-loss states and cumulative loss
    if (!loss_.empty()) {
        const double kap = c.obK / c.obCv;
        const double tEnd = t_ + dt;
        const double d = 0.5 * std::sqrt(kap * tEnd), kdt = kap * dt, Dd = 3 * d * d + kdt;
        for (size_t f = 0; f < loss_.size(); ++f) {
            LossState& L = loss_[f];
            const double th = T_[size_t(c.lossFaces[f].cell)] - L.Ti;
            const double p = (kdt * th / d + L.I - d * d * d * (th - L.thetaN) / kdt) / Dd;
            const double qq = 0.5 * ((th - L.thetaN) / kdt + 2 * p / d - th / (d * d));
            const double Q = c.obK * c.lossFaces[f].area * (th / d - p);
            cumLoss_ += Q * dt;
            L.I = th * d + p * d * d + 2 * qq * d * d * d;
            L.thetaN = th;
        }
    }
    // wells: surface rates, cumulatives, lagged yields and wellbore densities
    const int nk = nc_ + 1;
    for (size_t wi = 0; wi < c.wells.size(); ++wi) {
        const Well& w = c.wells[wi];
        WellReport& r = rep_[wi];
        r.qo = r.qw = r.qg = r.qinj = r.qcwe = r.qsteam = 0;
        if (!w.open) continue;
        r.cumMolP.resize(size_t(nk), 0.0); r.cumMolI.resize(size_t(nk), 0.0);
        for (int k = 0; k < nk && k < int(r.molar.size()); ++k) {
            const double v = r.molar[size_t(k)];
            if (v > 0) r.cumMolP[size_t(k)] += v * dt; else r.cumMolI[size_t(k)] -= v * dt;
        }
        if (w.injector) {
            double qt = 0; for (int k = 0; k < nk; ++k) qt -= r.molar[size_t(k)];
            r.qinj = qt;
            r.qcwe = -r.molar[size_t(nc_)] * steam::kMw / 1000.0;   // cold-water equivalent (m3/s)
            r.cumCwe += r.qcwe * dt;
        } else {
            double vo, vw, vg; surfaceVolumes(r.molar, vo, vw, vg);
            r.qo = vo; r.qw = vw; r.qg = vg;
            r.cumO += vo * dt; r.cumW += vw * dt; r.cumG += vg * dt;
            {   // C16: produced steam (CWE) - the water inflow of every connection split by its aqueous / vapour molar mobility
                // (shared drawdown p - p_w; the phase pressures differ only by Pc), the SPE4 "steam rate" of a producer
                const int baseY = 4;
                const double dref = w.refDepth >= 0 ? w.refDepth : depth_[size_t(w.cells.empty() ? 0 : w.cells[0])];
                const double rwb = wbRho_[wi] > 0 ? wbRho_[wi] : 800.0;
                double stm = 0, wat = 0;
                for (size_t j = 0; j < w.cells.size(); ++j) {
                    const int q = w.cells[j];
                    const double dp = std::max(0.0, p_[size_t(q)] - (r.bhp + rwb * kG * (depth_[size_t(q)] - dref)));
                    const double aW = Y_[size_t(q) * m_ + baseY + 0 * stride_ + 3 + nc_], aG = Y_[size_t(q) * m_ + baseY + 2 * stride_ + 3 + nc_];
                    const double f = w.wi[j] * permMult_[size_t(q)] * dp;
                    stm += f * aG; wat += f * (aW + aG);
                }
                const double nw = std::max(0.0, r.molar[size_t(nc_)]);
                r.qsteam = wat > 0 ? nw * stm / wat * steam::kMw / 1000.0 : 0.0;
                r.cumSteam += r.qsteam * dt;
            }
            // lagged surface yields per mole of component
            std::vector<double> yv(static_cast<size_t>(2 * nk), 0.0);
            double zh = 0; for (int k = 0; k < nc_; ++k) zh += std::max(r.molar[size_t(k)], 0.0);
            double nwp = std::max(r.molar[size_t(nc_)], 0.0);
            if (zh > 0) for (int k = 0; k < nc_; ++k) yv[size_t(k)] = vo / zh;
            yv[size_t(nk + nc_)] = nwp > 0 ? vw / nwp : steam::liquid(c.stdT, c.stdP).v * steam::kMw;
            if (zh > 0 || nwp > 0) surfYield_[wi] = yv;
        }
        r.cumH += r.hRate * dt;
        cumHin_ -= r.hRate * dt;
        // wellbore density: mobility-weighted mass density of the first connection cell (producers) / injected fluid
        const int q = w.cells.empty() ? 0 : w.cells[0];
        double num = 0, den = 0;
        const int base = 4, stride = stride_;
        for (int ph = 0; ph < 3; ++ph) {
            double mob = 0; for (int k = 0; k < nk; ++k) mob += Y_[size_t(q) * m_ + base + ph * stride + 3 + k];
            num += mob * Y_[size_t(q) * m_ + base + ph * stride]; den += mob;
        }
        if (!w.injector && den > 0) wbRho_[wi] = num / den;
        if (w.injector) {
            // injected-fluid density at the wellbore pressure (IF97; was a fixed 1.1e-3 / 0.1 m3/kg mixture)
            const double pw = std::clamp(r.bhp > 0 ? r.bhp : p_[size_t(q)], 1e4, 0.999 * steam::kPc);
            if (w.quality >= 0) {
                const double Ts = w.tinj > 0 ? std::min(w.tinj, steam::kTc - 1e-3) : steam::tsat(pw);
                const double Ps = steam::psat(Ts);
                const double vl = steam::liquid(Ts, std::max(pw, Ps)).v, vv = steam::vapour(Ts, std::min(pw, Ps)).v;
                wbRho_[wi] = 1.0 / ((1 - w.quality) * vl + w.quality * vv);
            } else {
                const double Ti = w.tinj > 0 ? w.tinj : 300.0;
                wbRho_[wi] = Ti < steam::tsat(pw) ? 1.0 / steam::liquid(Ti, pw).v : 1.0 / steam::vapour(Ti, pw).v;
                if (w.zinj[size_t(nc_)] < 0.5) wbRho_[wi] = std::min(wbRho_[wi], pw * 0.029 / (kRgas * Ti));   // gas (air) injector
            }
        }
    }
}

bool Simulator::solveStep(double dt, int& its) {
    bool conv = false; double res = 0;
    relax_ = 1.0; resHist_.clear();
    wellMode_.assign(c_.wells.size(), -1); wellSw_.assign(c_.wells.size(), 0);
    steamMode_.assign(c_.wells.size(), -1); steamSw_.assign(c_.wells.size(), 0);
    for (its = 0; its < 14; ++its) {
        newtonIt_ = its;
        if (!evalAll()) return false;
        if (!newtonStep(dt, conv, res)) return false;
        ++st_.newton;
        if (conv) return true;
    }
    return false;
}

// mechanics at the current flow state; the pore-volume law is re-centred on the mechanics porosity (fixed-stress
// split: PV = V_b [phi_mech + c_p (p - p_k) + c_T (T - T_k)]). 'change' = max |phi_mech - phi_flow|.
bool Simulator::mechanicsUpdate(bool iterate, double& change) {
    (void)iterate;
    const int n = n_;
    std::vector<double> dp(static_cast<size_t>(gn_), 0.0), dT(static_cast<size_t>(gn_), 0.0);
    for (int q = 0; q < n; ++q) {
        if (gmap_[size_t(q)] < 0) continue;
        const size_t g = size_t(gmap_[size_t(q)]);
        dp[g] = p_[size_t(q)] - p0_[size_t(q)]; dT[g] = T_[size_t(q)] - T0init_[size_t(q)];
    }
    if (!geo_.solve(dp, dT)) return false;
    ++st_.mechSolves; st_.mechCg += geo_.lastCg();
    change = 0;
    const std::vector<double>& phi = geo_.porosity();
    for (int q = 0; q < n; ++q) {
        if (gmap_[size_t(q)] < 0) continue;
        const double Vb = bulk_[size_t(q)];
        const double pf = poreVolume(q, p_[size_t(q)], T_[size_t(q)]) / Vb;
        const double pm = phi[size_t(gmap_[size_t(q)])];
        change = std::max(change, std::fabs(pm - pf));
        pvA_[size_t(q)] = Vb * pm; pvP_[size_t(q)] = p_[size_t(q)]; pvT_[size_t(q)] = T_[size_t(q)];
    }
    return true;
}

bool Simulator::advanceTo(double tEnd, double dtMax) {
    const int n = n_;
    const bool mech = geo_.ready();
    while (t_ < tEnd * (1 - 1e-12)) {
        double dt = std::min({dt_, dtMax, tEnd - t_});
        if (tEnd - t_ - dt < 1e-3 * dt) dt = tEnd - t_;
        const bool forced = forcedIdx_ < forced_.size();
        if (forced) { dt = std::min(forced_[forcedIdx_], tEnd - t_); if (tEnd - t_ - dt < 1e-6 * dt) dt = tEnd - t_; }
        Xn_ = X_;
        const std::vector<double> T0 = T_, S0 = sw_, G0 = sg_, P0 = p_;
        const std::vector<WellReport> rep0 = rep_;
        const std::vector<double> pvA0 = pvA_, pvP0 = pvP_, pvT0 = pvT_;
        bool ok = false; int its = 0, cutAttempts = 0;
        for (int attempt = 0; attempt < 14 && !ok; ++attempt) {
            ok = solveStep(dt, its);
            if (ok && mech && c_.geo.coupling == 1) {
                // iterative fixed-stress coupling: mechanics -> porosity -> flow until the porosity change is small
                for (int k = 0; k < std::max(1, c_.geo.maxIter); ++k) {
                    double ch = 0;
                    if (!mechanicsUpdate(true, ch)) { ok = false; break; }
                    ++st_.mechIters;
                    if (ch < c_.geo.tol) break;
                    int its2 = 0;
                    if (!solveStep(dt, its2)) { ok = false; break; }
                    its = std::max(its, its2);
                }
            }
            if (ok) break;
            // cut
            X_ = Xn_; T_ = T0; rep_ = rep0; pvA_ = pvA0; pvP_ = pvP0; pvT_ = pvT0;
            evalAll();
            dt *= 0.25; ++st_.cuts; ++cutAttempts;
            if (dt < 1.0 || forced) return false;
        }
        if (!ok) return false;
        if (forced) ++forcedIdx_;
        if (adjOn_) adjRecord(t_, dt);
        finishStep(dt);
        if (mech && c_.geo.coupling == 0) { double ch = 0; mechanicsUpdate(false, ch); }
        if (mech) geo_.commit();   // C16: accept the plastic state of this step
        updatePermeability();
        if (!c_.fracCells.empty()) updateFractures(true);
        t_ += dt; dtLast_ = dt; ++st_.steps;
        // next step size from the changes (temperature 30 K, saturation 0.2, pressure 10 %) and the Newton count
        double dT = 0, dS = 0, dP = 0;
        for (int q = 0; q < n; ++q) {
            dT = std::max(dT, std::fabs(T_[size_t(q)] - T0[size_t(q)]));
            dS = std::max({dS, std::fabs(sw_[size_t(q)] - S0[size_t(q)]), std::fabs(sg_[size_t(q)] - G0[size_t(q)])});
            dP = std::max(dP, std::fabs(p_[size_t(q)] - P0[size_t(q)]) / std::max(P0[size_t(q)], 1e5));
        }
        // Aziz-Settari change control  f = (1 + w) eta / (d + w eta)  (w = 1/2, smoother than eta / d); growth after a time
        // step cut is limited to 1.25 for the next three steps (stops the cut -> double -> cut cycle on steam fronts)
        auto ctl = [](double eta, double d) { const double w = 0.5; return (1 + w) * eta / (d + w * eta); };
        double f = 2.0;
        if (dT > 0) f = std::min(f, ctl(30.0, dT));
        if (dS > 0) f = std::min(f, ctl(0.2, dS));
        if (dP > 0) f = std::min(f, ctl(0.1, dP));
        if (its > 8) f = std::min(f, 0.7);
        if (cutAttempts > 0) sinceCut_ = 0; else ++sinceCut_;
        if (sinceCut_ < 3) f = std::min(f, 1.25);
        f = std::max(f, 0.3);
        dt_ = std::max(1.0, dt * f);
    }
    return true;
}

double Simulator::jacobianCheck(double dt, int* worstRow, int* worstCol, int onlyCell) {
    const int n = n_, nb = nb_;
    const std::vector<double> X0 = X_;
    const std::vector<WellReport> rep0 = rep_;
    bool conv; double res;
    assembleOnly_ = true;
    if (!evalAll()) { assembleOnly_ = false; return -1; }
    newtonStep(dt, conv, res);
    const std::vector<double> R0 = lastR_, J0 = lastJ_;
    const BSR& J = J_;
    std::vector<double> rowNorm(R0.size(), 0.0);
    for (int q = 0; q < n; ++q)
        for (int k = J.rowPtr[size_t(q)]; k < J.rowPtr[size_t(q) + 1]; ++k)
            for (int r = 0; r < nb; ++r)
                for (int j = 0; j < nb; ++j) rowNorm[size_t(q) * nb + r] = std::max(rowNorm[size_t(q) * nb + r], std::fabs(J0[size_t(k) * nb * nb + r * nb + j]) * (j == 0 ? X0[size_t(J.col[size_t(k)]) * nb] : 1.0));
    double worst = 0;
    for (int qc = 0; qc < n; ++qc)
        for (int j = 0; j < nb; ++j) {
            if (onlyCell >= 0 && qc != onlyCell) break;
            double Nt = 0; for (int k = 1; k < nb - 1; ++k) Nt += std::max(X0[size_t(qc) * nb + k], 0.0);
            const double h = j == 0 ? 1e-6 * X0[size_t(qc) * nb] : (j == nb - 1 ? 1e-5 * (crock_[size_t(qc)] + 80 * Nt) : 1e-6 * Nt);
            std::vector<double> Rp, Rm;
            X_ = X0; X_[size_t(qc) * nb + j] += h; rep_ = rep0;
            if (!evalAll()) continue;
            newtonStep(dt, conv, res); Rp = lastR_;
            X_ = X0; X_[size_t(qc) * nb + j] -= h; rep_ = rep0;
            if (!evalAll()) continue;
            newtonStep(dt, conv, res); Rm = lastR_;
            for (int q = 0; q < n; ++q) {
                int kb = -1;
                for (int k = J.rowPtr[size_t(q)]; k < J.rowPtr[size_t(q) + 1]; ++k) if (J.col[size_t(k)] == qc) kb = k;
                for (int r = 0; r < nb; ++r) {
                    const double fd = (Rp[size_t(q) * nb + r] - Rm[size_t(q) * nb + r]) / (2 * h);
                    const double an = kb >= 0 ? J0[size_t(kb) * nb * nb + r * nb + j] : 0.0;
                    const double scale = (j == 0 ? X0[size_t(qc) * nb] : 1.0);
                    const double e = std::fabs(fd - an) * scale / std::max(rowNorm[size_t(q) * nb + r], 1e-300);
                    if (e > worst) { worst = e; if (worstRow) *worstRow = q * nb + r; if (worstCol) *worstCol = qc * nb + j;
                        if (std::getenv("THERM_DEBUG")) std::fprintf(stderr, "  jac q %d r %d qc %d j %d fd %.6e an %.6e (fwd %.6e) rownorm %.3e\n", q, r, qc, j, fd, an,
                            (Rp[size_t(q) * nb + r] - R0[size_t(q) * nb + r]) / h, rowNorm[size_t(q) * nb + r]); }
                }
            }
        }
    X_ = X0; rep_ = rep0; assembleOnly_ = false;
    evalAll();
    return worst;
}

// ------------------------------------------------------------------------------------------------ C10: compositions, trapping
double Simulator::trappedGas(int q, double sh) const {
    const Case& c = c_;
    double d = c.sgc;
    if (!c.slT.empty()) { d = 0; for (size_t k = 0; k < c.slT.size() && k < c.krgT.size(); ++k) if (c.krgT[k] <= 0) d = std::max(d, 1 - c.slT[k]); }
    (void)q;
    if (!c.hyst.on) return std::min(sh, d);
    const double i = std::max(c.hyst.sgtMax, d), mx = 1 - c.swc;
    if (i <= d + 1e-9 || sh <= d) return std::min(sh, d);
    const double C = 1.0 / (i - d) - 1.0 / std::max(1e-9, mx - d);
    return d + (sh - d) / (1.0 + C * (sh - d));
}

std::vector<double> Simulator::trappedGasSaturation() const {
    std::vector<double> t(static_cast<size_t>(n_), 0.0);
    for (int q = 0; q < n_; ++q) {
        const double sh = q < int(sgHist_.size()) ? std::max(sgHist_[size_t(q)], sg_[size_t(q)]) : sg_[size_t(q)];
        t[size_t(q)] = std::min(sg_[size_t(q)], trappedGas(q, sh));
    }
    return t;
}

void Simulator::phaseCompositions(std::vector<double>& x, std::vector<double>& y, std::vector<double>& a) const {
    const size_t nk = size_t(nc_) + 1;
    x.assign(size_t(n_) * nk, 0.0); y.assign(size_t(n_) * nk, 0.0); a.assign(size_t(n_) * nk, 0.0);
    std::vector<double> Yb(static_cast<size_t>(m_));
    for (int q = 0; q < n_; ++q) {
        Local L;
        if (!evalCell(q, &X_[size_t(q) * nb_], Yb.data(), &L, T_[size_t(q)])) continue;
        for (size_t k = 0; k < nk; ++k) {
            if (L.on[O] && k < L.x.size() && k < size_t(nc_)) x[size_t(q) * nk + k] = L.x[k];
            if (L.on[G] && k < L.y.size()) y[size_t(q) * nk + k] = L.y[k];
            if (L.on[W] && k < L.a.size()) a[size_t(q) * nk + k] = L.a[k];
        }
    }
}

Inventory Simulator::trackedInventory() const {
    Inventory I;
    const Case& c = c_;
    const int k = c.track;
    if (k < 0 || k >= nc_) return I;
    const double mw = c.comps[size_t(k)].mw;
    std::vector<Local> Ls(static_cast<size_t>(n_));
    std::vector<double> Yb(static_cast<size_t>(m_));
    std::vector<char> ok(static_cast<size_t>(n_), 0);
    for (int q = 0; q < n_; ++q) ok[size_t(q)] = evalCell(q, &X_[size_t(q) * nb_], Yb.data(), &Ls[size_t(q)], T_[size_t(q)]) ? 1 : 0;
    // gas outflow per cell (reservoir m3/s, upstream weighting as the flux assembly)
    std::vector<double> gout(static_cast<size_t>(n_), 0.0);
    for (const Case::Conn& e : conns_) {
        const int a = e.a, b = e.b;
        if (!ok[size_t(a)] || !ok[size_t(b)]) continue;
        const Local& A = Ls[size_t(a)]; const Local& B = Ls[size_t(b)];
        if (!A.on[G] && !B.on[G]) continue;
        const double rho = A.on[G] && B.on[G] ? 0.5 * (A.rho[G] + B.rho[G]) : A.on[G] ? A.rho[G] : B.rho[G];
        const double dPhi = (p_[size_t(a)] - p_[size_t(b)]) - rho * kG * (depth_[size_t(a)] - depth_[size_t(b)]);
        const int up = dPhi >= 0 ? a : b;
        const Local& U = Ls[size_t(up)];
        if (!U.on[G]) continue;
        gout[size_t(up)] += e.T * U.kr[G] / U.mu[G] * std::fabs(dPhi);
    }
    const double year = 365.25 * 86400.0;
    for (int q = 0; q < n_; ++q) {
        if (!ok[size_t(q)]) continue;
        const Local& L = Ls[size_t(q)];
        if (L.on[W] && size_t(k) < L.a.size()) I.dissolvedAq += L.n[W] * L.a[size_t(k)] * mw;
        if (L.on[O] && size_t(k) < L.x.size()) I.dissolvedOil += L.n[O] * L.x[size_t(k)] * mw;
        if (!L.on[G] || size_t(k) >= L.y.size()) continue;
        const double mG = L.n[G] * L.y[size_t(k)] * mw;
        const double Sg = L.S[G];
        const double sh = q < int(sgHist_.size()) ? std::max(sgHist_[size_t(q)], Sg) : Sg;
        const double sgt = trappedGas(q, sh);
        const double fr = Sg > 0 ? std::min(1.0, sgt / Sg) : 1.0;
        I.residual += fr * mG;
        if (L.kr[G] <= 0) I.immobileKr += mG;
        const double mm = (1 - fr) * mG;
        if (mm <= 0) continue;
        const double vMob = (1 - fr) * Sg * L.pv;
        if (gout[size_t(q)] * year < c.structTol * vMob) I.structural += mm; else I.mobile += mm;
    }
    for (int q = 0; q < n_; ++q)
        for (int i = 0; i < nc_; ++i) {
            const Component& s = c.comps[size_t(i)];
            if (!s.solid() || s.trackEq == 0) continue;
            I.mineral += (X_[size_t(q) * nb_ + 1 + i] - solidInit_[size_t(q) * nc_ + i]) * s.trackEq * mw;
        }
    return I;
}

// ------------------------------------------------------------------------------------------------ C12 EDFM
bool addFractures(Case& c, const std::vector<Fracture>& fr, const FractureStress& st, std::string* error) {
    auto fail = [&](const std::string& m) { if (error) *error = "EDFM: " + m; return false; };
    if (fr.empty()) return true;
    const int nx = c.nx, ny = c.ny, nz = c.nz, nbox = nx * ny * nz;
    if (nbox <= 0) return fail("empty grid");
    // general form already built (deck model of an all-active Cartesian grid: cell q = box cell q, per-cell dx/dy/dz given)
    const bool general = !c.cellPV.empty();
    if (general && (int(c.cellPV.size()) != nbox || int(c.dx.size()) != nbox || int(c.dy.size()) != nbox || int(c.dz.size()) != nbox))
        return fail("the case uses general connections that are not an all-active Cartesian box");
    auto at = [](const std::vector<double>& v, int i, double d) { return v.empty() ? d : v.size() == 1 ? v[0] : v[size_t(i)]; };
    auto idx = [&](int i, int j, int k) { return (k * ny + j) * nx + i; };
    // ---- the box in general-connection form (same geometry / transmissibilities as Simulator::setup)
    std::vector<double> DX(static_cast<size_t>(nbox)), DY(static_cast<size_t>(nbox)), DZ(static_cast<size_t>(nbox));
    std::vector<double> ztopK(static_cast<size_t>(nz) + 1, c.top);
    for (int k = 0; k < nz; ++k) ztopK[size_t(k) + 1] = ztopK[size_t(k)] + at(c.dz, idx(0, 0, k), 1.0);
    if (general) {
        for (int q = 0; q < nbox; ++q) { DX[size_t(q)] = c.dx[size_t(q)]; DY[size_t(q)] = c.dy[size_t(q)]; DZ[size_t(q)] = c.dz[size_t(q)]; }
        for (int k = 0; k < nz; ++k) {
            const int q = idx(0, 0, k);
            if (q < int(c.cellDepth.size())) ztopK[size_t(k)] = c.cellDepth[size_t(q)] - 0.5 * DZ[size_t(q)];
            ztopK[size_t(k) + 1] = ztopK[size_t(k)] + DZ[size_t(q)];
        }
    }
    if (!general) c.cellPV.assign(size_t(nbox), 0), c.cellBulk.assign(size_t(nbox), 0), c.cellDepth.assign(size_t(nbox), 0);
    for (int k = 0; k < nz && !general; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const int q = idx(i, j, k);
                DX[size_t(q)] = at(c.dx, q, 1.0); DY[size_t(q)] = at(c.dy, q, 1.0); DZ[size_t(q)] = at(c.dz, q, 1.0);
                c.cellBulk[size_t(q)] = DX[size_t(q)] * DY[size_t(q)] * DZ[size_t(q)];
                c.cellPV[size_t(q)] = c.cellBulk[size_t(q)] * at(c.poro, q, 0.3);
                c.cellDepth[size_t(q)] = ztopK[size_t(k)] + 0.5 * DZ[size_t(q)];
            }
    if (!general) c.conns.clear();
    auto addConn = [&](int a, int b, int dir) {
        double ka, kb, la, lb, A;
        if (dir == 0) { ka = at(c.permx, a, 1e-13); kb = at(c.permx, b, 1e-13); la = DX[size_t(a)]; lb = DX[size_t(b)]; A = DY[size_t(a)] * DZ[size_t(a)]; }
        else if (dir == 1) { ka = at(c.permy, a, at(c.permx, a, 1e-13)); kb = at(c.permy, b, at(c.permx, b, 1e-13)); la = DY[size_t(a)]; lb = DY[size_t(b)]; A = DX[size_t(a)] * DZ[size_t(a)]; }
        else { ka = at(c.permz, a, 0.1 * at(c.permx, a, 1e-13)); kb = at(c.permz, b, 0.1 * at(c.permx, b, 1e-13)); la = DZ[size_t(a)]; lb = DZ[size_t(b)]; A = DX[size_t(a)] * DY[size_t(a)]; }
        const double T = (ka > 0 && kb > 0) ? A / (0.5 * la / ka + 0.5 * lb / kb) : 0.0;
        c.conns.push_back({a, b, T, A / (0.5 * la + 0.5 * lb)});
    };
    for (int k = 0; k < nz && !general; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const int q = idx(i, j, k);
                if (i + 1 < nx) addConn(q, idx(i + 1, j, k), 0);
                if (j + 1 < ny) addConn(q, idx(i, j + 1, k), 1);
                if (k + 1 < nz) addConn(q, idx(i, j, k + 1), 2);
            }
    if (c.heatLoss && c.lossFaces.empty() && !general)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const int t = idx(i, j, 0), b = idx(i, j, nz - 1);
                c.lossFaces.push_back({t, DX[size_t(t)] * DY[size_t(t)]});
                c.lossFaces.push_back({b, DX[size_t(b)] * DY[size_t(b)]});
            }
    if (c.geo.enabled && c.geoCell.empty()) {
        c.geoNx = nx; c.geoNy = ny; c.geoNz = nz; c.geoDx = c.dx; c.geoDy = c.dy; c.geoDz = c.dz; c.geoTop = c.top;
        c.geoCell.resize(size_t(nbox));
        for (int q = 0; q < nbox; ++q) c.geoCell[size_t(q)] = q;
    }
    // column faces (rows of the first layer / column)
    std::vector<double> xf(static_cast<size_t>(nx) + 1, 0.0), yf(static_cast<size_t>(ny) + 1, 0.0);
    for (int i = 0; i < nx; ++i) xf[size_t(i) + 1] = xf[size_t(i)] + DX[size_t(idx(i, 0, 0))];
    for (int j = 0; j < ny; ++j) yf[size_t(j) + 1] = yf[size_t(j)] + DY[size_t(idx(0, j, 0))];
    // ---- fracture pieces
    struct Piece { int frac, i, j, k, cell, host; double ta, tb, h, w, kf; };
    std::vector<Piece> pcs;
    std::vector<int> hostOf;
    const double azi = st.azimuthDeg * M_PI / 180.0, ax = std::cos(azi), ay = std::sin(azi);
    for (int f = 0; f < int(fr.size()); ++f) {
        const Fracture& F = fr[size_t(f)];
        const double Lx = F.x1 - F.x0, Ly = F.y1 - F.y0, L = std::hypot(Lx, Ly);
        if (!(L > 0)) return fail(F.name + ": zero-length trace");
        if (!(F.aperture > 0)) return fail(F.name + ": aperture must be > 0");
        const double tx = Lx / L, ty = Ly / L, nxv = -ty, nyv = tx;
        const double zt = F.zTop < 0 ? ztopK[0] : F.zTop, zb = F.zBot < 0 ? ztopK[size_t(nz)] : F.zBot;
        const double kf = F.perm > 0 ? F.perm : F.aperture * F.aperture / 12.0;
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                // Liang-Barsky clip of the trace against the column rectangle
                double t0 = 0, t1 = L;
                const double pp[4] = {-tx, tx, -ty, ty};
                const double qq[4] = {F.x0 - xf[size_t(i)], xf[size_t(i) + 1] - F.x0, F.y0 - yf[size_t(j)], yf[size_t(j) + 1] - F.y0};
                bool in = true;
                for (int e = 0; e < 4 && in; ++e) {
                    if (std::fabs(pp[e]) < 1e-14) { if (qq[e] < 0) in = false; continue; }
                    const double r = qq[e] / pp[e];
                    if (pp[e] < 0) t0 = std::max(t0, r); else t1 = std::min(t1, r);
                }
                if (!in || t1 - t0 <= 1e-9 * L) continue;
                // mean distance of the host column from the fracture plane (20 x 20 midpoint rule)
                double dsum = 0;
                const double hx = (xf[size_t(i) + 1] - xf[size_t(i)]) / 20, hy = (yf[size_t(j) + 1] - yf[size_t(j)]) / 20;
                for (int a = 0; a < 20; ++a)
                    for (int b = 0; b < 20; ++b) {
                        const double x = xf[size_t(i)] + (a + 0.5) * hx, y = yf[size_t(j)] + (b + 0.5) * hy;
                        dsum += std::fabs((x - F.x0) * nxv + (y - F.y0) * nyv);
                    }
                const double dAvg = std::max(dsum / 400, 1e-6);
                for (int k = 0; k < nz; ++k) {
                    const double h = std::min(zb, ztopK[size_t(k) + 1]) - std::max(zt, ztopK[size_t(k)]);
                    if (h <= 1e-9) continue;
                    const int host = idx(i, j, k);
                    const int cell = int(c.cellPV.size());
                    const double A = (t1 - t0) * h;
                    c.cellPV.push_back(F.aperture * A); c.cellBulk.push_back(F.aperture * A);
                    c.cellDepth.push_back(std::max(zt, ztopK[size_t(k)]) + 0.5 * h);
                    hostOf.push_back(host);
                    const double kn = at(c.permx, host, 1e-13) * nxv * nxv + at(c.permy, host, at(c.permx, host, 1e-13)) * nyv * nyv;
                    c.conns.push_back({host, cell, 2 * A * kn / dAvg, 2 * A / dAvg});
                    pcs.push_back({f, i, j, k, cell, host, t0, t1, h, F.aperture, kf});
                    // mechanics of the plane
                    Case::FracCell fc; fc.cell = cell; fc.frac = f; fc.host = host; fc.w0 = F.aperture; fc.segLen = t1 - t0; fc.segH = h; fc.kf = kf;
                    fc.xa = F.x0 + tx * t0; fc.ya = F.y0 + ty * t0; fc.xb = F.x0 + tx * t1; fc.yb = F.y0 + ty * t1;
                    fc.zt = std::max(zt, ztopK[size_t(k)]); fc.zb = fc.zt + h;
                    const double depth = c.cellDepth.back();
                    const double SH = st.sHmaxGrad * depth, Sh = st.sHminGrad * depth;
                    const double na = nxv * ax + nyv * ay;
                    fc.sn0 = SH * na * na + Sh * (1 - na * na);
                    fc.tau0 = std::fabs(SH - Sh) * std::fabs(na) * std::sqrt(std::max(0.0, 1 - na * na));
                    fc.p0 = at(c.p0, host, 1e7); fc.T0 = at(c.T0, host, 300.0);
                    if (F.mech) {
                        fc.kn0 = F.kn0; fc.vm = F.vm; fc.mu = F.mu; fc.ks = F.ks; fc.tanPsi = std::tan(F.dilationDeg * M_PI / 180.0);
                        fc.thermStress = st.alphaT > 0 ? st.E * st.alphaT / (1 - st.nu) : 0.0;
                        const double s0 = std::max(0.0, fc.sn0 - fc.p0);
                        fc.wOpen = F.aperture + (fc.kn0 > 0 && fc.vm > 0 ? fc.vm * s0 / (fc.kn0 * fc.vm + s0) : 0.0);
                    } else { fc.kn0 = 0; fc.vm = 0; fc.ks = 0; fc.wOpen = F.aperture; }
                    c.fracCells.push_back(fc);
                }
            }
    }
    if (pcs.empty()) return fail("no fracture intersects the grid");
    // ---- fracture-fracture connections: along each trace (same layer), vertically (same column), at intersections
    auto ffConn = [&](int a, int b, double ta, double tb, double g) {
        if (!(ta > 0) || !(tb > 0)) return;
        c.fracFF.push_back(int(c.conns.size()));
        c.conns.push_back({a, b, ta * tb / (ta + tb), g});
    };
    for (size_t u = 0; u < pcs.size(); ++u)
        for (size_t v = u + 1; v < pcs.size(); ++v) {
            const Piece& A = pcs[u]; const Piece& B = pcs[v];
            const double La = A.tb - A.ta, Lb = B.tb - B.ta;
            if (A.frac == B.frac) {
                if (A.k == B.k && (std::fabs(A.tb - B.ta) < 1e-9 * (La + Lb) || std::fabs(B.tb - A.ta) < 1e-9 * (La + Lb))) {
                    const double h = std::min(A.h, B.h);
                    ffConn(A.cell, B.cell, A.kf * A.w * h / (0.5 * La), B.kf * B.w * h / (0.5 * Lb), A.w * h / (0.5 * (La + Lb)));
                } else if (A.i == B.i && A.j == B.j && std::abs(A.k - B.k) == 1) {
                    ffConn(A.cell, B.cell, A.kf * A.w * La / (0.5 * A.h), B.kf * B.w * Lb / (0.5 * B.h), A.w * La / (0.5 * (A.h + B.h)));
                }
            } else if (A.i == B.i && A.j == B.j && A.k == B.k) {
                // intersection of the two traces inside this column
                const Fracture& Fa = fr[size_t(A.frac)]; const Fracture& Fb = fr[size_t(B.frac)];
                const double ax1 = Fa.x1 - Fa.x0, ay1 = Fa.y1 - Fa.y0, bx1 = Fb.x1 - Fb.x0, by1 = Fb.y1 - Fb.y0;
                const double den = ax1 * by1 - ay1 * bx1;
                if (std::fabs(den) < 1e-14) continue;
                const double sa = ((Fb.x0 - Fa.x0) * by1 - (Fb.y0 - Fa.y0) * bx1) / den;   // fraction along A's trace
                const double sb = ((Fb.x0 - Fa.x0) * ay1 - (Fb.y0 - Fa.y0) * ax1) / den;
                const double LA = std::hypot(ax1, ay1), LB = std::hypot(bx1, by1);
                const double tA = sa * LA, tB = sb * LB;
                if (tA < A.ta || tA > A.tb || tB < B.ta || tB > B.tb) continue;
                // mean distance of each piece from the intersection line (Moinfar et al. 2014): (a^2 + b^2) / (2 (a + b))
                const double a1 = tA - A.ta, a2 = A.tb - tA, b1 = tB - B.ta, b2 = B.tb - tB;
                const double dA = std::max(1e-6, (a1 * a1 + a2 * a2) / (2 * (a1 + a2))), dB = std::max(1e-6, (b1 * b1 + b2 * b2) / (2 * (b1 + b2)));
                const double h = std::min(A.h, B.h);
                ffConn(A.cell, B.cell, A.kf * A.w * h / dA, B.kf * B.w * h / dB, h * 0.5 * (A.w + B.w) / (0.5 * (dA + dB)));
            }
        }
    // ---- per-cell arrays: fracture cells take their host's values
    auto ext = [&](std::vector<double>& v) {
        if (int(v.size()) != nbox) return;
        for (int h : hostOf) v.push_back(v[size_t(h)]);
    };
    auto extV = [&](std::vector<std::vector<double>>& v) {
        if (int(v.size()) != nbox) return;
        for (int h : hostOf) v.push_back(v[size_t(h)]);
    };
    ext(c.p0); ext(c.T0); ext(c.sw0); ext(c.sg0); ext(c.thconr); ext(c.thconsfCell); ext(c.rockCvCell); ext(c.rockCvTCell);
    ext(c.poro); ext(c.permx); ext(c.permy); ext(c.permz);
    extV(c.x0); extV(c.y0); extV(c.solid0);
    if (c.geo.enabled) for (size_t f = 0; f < hostOf.size(); ++f) c.geoCell.push_back(-1);
    connectWellsToFractures(c, c.wells);
    return true;
}

int connectWellsToFractures(const Case& c, std::vector<Well>& wells, double rw) {
    int added = 0;
    for (Well& w : wells) {
        const size_t nOrig = w.cells.size();
        for (size_t k = 0; k < nOrig; ++k) {
            const int host = w.cells[k];
            for (const Case::FracCell& fc : c.fracCells) {
                if (fc.host != host) continue;
                if (std::find(w.cells.begin(), w.cells.end(), fc.cell) != w.cells.end()) continue;
                const double re = 0.14 * std::sqrt(fc.segLen * fc.segLen + fc.segH * fc.segH);
                const double wi = 2 * M_PI * fc.kf * fc.w0 / std::max(std::log(std::max(re, 1e-6) / rw), 1.0);
                w.cells.push_back(fc.cell); w.wi.push_back(wi); ++added;
            }
        }
    }
    return added;
}

// ------------------------------------------------------------------------------------------------ C12 sorption
bool addSorption(Case& c, const Sorption& s, std::string* error) {
    auto fail = [&](const std::string& m) { if (error) *error = "SORB " + s.comp + ": " + m; return false; };
    auto upc = [](std::string v) { for (char& ch : v) ch = char(std::toupper(static_cast<unsigned char>(ch))); return v; };
    int ic = -1;
    for (size_t i = 0; i < c.comps.size(); ++i) if (upc(c.comps[i].name) == upc(s.comp)) ic = int(i);
    if (ic < 0) return fail("no such component");
    if (c.comps[size_t(ic)].solid()) return fail("a solid cannot sorb");
    if (!(s.capacity > 0) || !(s.bRef > 0)) return fail("capacity and b must be > 0");
    // insert a component at position 'at' (before water) in every per-component array
    auto insertComp = [&](const Component& k) {
        const size_t at = c.comps.size();
        c.comps.push_back(k);
        for (auto& r : c.reactions) {
            r.sto.insert(r.sto.begin() + long(std::min(at, r.sto.size())), 0.0);
            r.prod.insert(r.prod.begin() + long(std::min(at, r.prod.size())), 0.0);
            if (!r.order.empty()) r.order.insert(r.order.begin() + long(std::min(at, r.order.size())), 0.0);
            if (!r.phase.empty()) r.phase.insert(r.phase.begin() + long(std::min(at, r.phase.size())), -1);
            if (r.o2pp >= int(at)) ++r.o2pp;
        }
        for (auto& v : c.solid0) if (v.size() >= at) v.insert(v.begin() + long(at), 0.0);
        for (auto& v : c.x0) if (v.size() >= at) v.insert(v.begin() + long(at), 0.0);
        for (auto& v : c.y0) if (v.size() >= at) v.insert(v.begin() + long(at), 0.0);
        for (auto& w : c.wells) if (w.zinj.size() >= at) w.zinj.insert(w.zinj.begin() + long(at), 0.0);
        return int(at);
    };
    const Component K = c.comps[size_t(ic)];   // copy: insertComp reallocates
    int is = -1;
    for (size_t i = 0; i < c.comps.size(); ++i) if (c.comps[i].solid() && upc(c.comps[i].name) == upc(s.site)) is = int(i);
    const bool newSite = is < 0;
    if (newSite) {
        Component st; st.name = upc(s.site); st.type = 2; st.dead = true; st.mw = 1e-6; st.rhoRef = 1e9;
        is = insertComp(st);
    }
    Component ad; ad.name = upc(K.name) + "_ADS"; ad.type = 2; ad.dead = true;
    ad.mw = c.comps[size_t(ic)].mw + c.comps[size_t(is)].mw;
    ad.rhoRef = s.adsRho > 0 ? s.adsRho : 1e9;
    ad.cpl1 = c.comps[size_t(ic)].cpg1 > 0 ? c.comps[size_t(ic)].cpg1 : c.comps[size_t(ic)].cpl1;
    ad.formula = c.comps[size_t(ic)].formula;
    const int ia = insertComp(ad);
    const size_t nc = c.comps.size();
    Reaction r; r.name = "SORB_" + upc(K.name);
    r.sto.assign(nc + 1, 0.0); r.prod.assign(nc + 1, 0.0);
    r.sto[size_t(ic)] = 1; r.sto[size_t(is)] = 1; r.prod[size_t(ia)] = 1;
    const Component& Kc = c.comps[size_t(ic)];
    const int ph = s.phase >= 0 ? s.phase : (Kc.volatileK() ? 2 : 0);
    if (ph == 2) r.o2pp = ic;
    else { r.phase.assign(nc + 1, -1); r.phase[size_t(ic)] = ph; }
    r.reversible = true;
    r.keqA = std::log(s.bRef) - s.dH / (kRgas * s.tRef);   // ln b(T) = ln bRef + dH/R (1/T - 1/tRef)
    r.keqB = s.dH / kRgas;
    r.A = s.bRef / (2.0 * std::max(s.tau, 1e-6));
    r.dH = s.dH;
    c.reactions.push_back(r);
    // initial sites (mol per m3 of pore volume) and the equilibrium loading with the initial gas
    const int n = c.cellPV.empty() ? c.nx * c.ny * c.nz : int(c.cellPV.size());
    auto phiOf = [&](int q) {
        if (!c.cellPV.empty() && !c.cellBulk.empty()) return c.cellPV[size_t(q)] / std::max(1e-30, c.cellBulk[size_t(q)]);
        return c.poro.size() == 1 ? c.poro[0] : c.poro[size_t(q)];
    };
    if (c.solid0.empty()) c.solid0.assign(1, std::vector<double>(nc, 0.0));
    if (c.solid0.size() == 1 && n > 1) c.solid0.assign(size_t(n), c.solid0[0]);
    for (auto& v : c.solid0) v.resize(nc, 0.0);
    for (int q = 0; q < int(c.solid0.size()); ++q) {
        const double phi = std::max(1e-6, phiOf(q));
        std::vector<double>& v = c.solid0[size_t(q)];
        double load = 0;
        if (s.initEquilibrium && ph == 2) {
            const double sg = c.sg0.empty() ? 0.0 : (c.sg0.size() == 1 ? c.sg0[0] : c.sg0[size_t(q)]);
            const double p = c.p0.empty() ? 1e5 : (c.p0.size() == 1 ? c.p0[0] : c.p0[size_t(q)]);
            const double T = c.T0.empty() ? s.tRef : (c.T0.size() == 1 ? c.T0[0] : c.T0[size_t(q)]);
            const auto& yv = c.y0.empty() ? c.x0 : c.y0;
            const double y = yv.empty() ? 0.0 : (yv.size() == 1 ? yv[0] : yv[size_t(q)]).size() > size_t(ic) ? (yv.size() == 1 ? yv[0] : yv[size_t(q)])[size_t(ic)] : 0.0;
            if (sg > 0 && y > 0) {
                const double bT = std::exp(r.keqA + r.keqB / T), pp = y * p * 1e-3;
                const double freeSites = newSite ? s.capacity : v[size_t(is)] * phi;   // shared site: what is left
                load = freeSites * bT * pp / (1 + bT * pp);
            }
        }
        if (newSite) v[size_t(is)] = (s.capacity - load) / phi;
        else v[size_t(is)] -= load / phi;
        v[size_t(ia)] = load / phi;
    }
    return true;
}


// ---------------------------------------------------------------- C14 restart state
template <class Ar> void Simulator::stateIO(Ar& ar) {
    int n = n_, nc = nc_; uint64_t nw = c_.wells.size(), nf = c_.fracCells.size(), ncn = conns_.size();
    const int n0 = n, nc0 = nc; const uint64_t nw0 = nw, nf0 = nf, ncn0 = ncn;
    ar(n, nc, nw, nf, ncn);
    if (Ar::reading && (n != n0 || nc != nc0 || nw != nw0 || nf != nf0 || ncn != ncn0))
        throw std::runtime_error("thermal restart state does not match the model (cells " + std::to_string(n) + " / " + std::to_string(n0) +
                                 ", components " + std::to_string(nc) + " / " + std::to_string(nc0) + ", wells " + std::to_string(nw) + " / " +
                                 std::to_string(nw0) + ", fracture cells " + std::to_string(nf) + " / " + std::to_string(nf0) + ")");
    ar(t_, dt_, dtLast_, cumLoss_, cumHin_, cumReact_, reactRate_, sinceCut_);
    ar(X_, T_, p_, sw_, so_, sg_, loss_, sgHist_, permMult_);
    ar(pvA_, pvCp_, pvP_, pvCt_, pvT_, fracW_, surfYield_, wbRho_, wellRamp_, wellAct_, st_);
    for (Case::FracCell& F : c_.fracCells) ar(F.wShear);
    for (Case::Conn& e : conns_) ar(e.T);
    for (WellReport& r : rep_) {
        ar(r.bhp, r.onBhp, r.trapActive, r.qo, r.qw, r.qg, r.qinj, r.qcwe, r.hRate, r.tFluid, r.molar);
        ar(r.cumO, r.cumW, r.cumG, r.cumCwe, r.cumH, r.cumMolP, r.cumMolI);
    }
    uint8_t mech = geo_.ready() ? 1 : 0; const uint8_t mech0 = mech;
    ar(mech);
    if (Ar::reading && mech != mech0) throw std::runtime_error("thermal restart: geomechanics on / off differs from the model");
    if (mech) geo_.stateIO(ar);
}

std::string Simulator::saveState() const {
    StateOut o;
    o.w.str("RESSIM-THERMAL-STATE-1");
    const_cast<Simulator*>(this)->stateIO(o);
    return o.w.buf;
}

void Simulator::loadState(const std::string& blob) {
    StateIn in(blob);
    if (in.r.str() != "RESSIM-THERMAL-STATE-1") throw std::runtime_error("not a thermal restart state");
    if (rep_.size() != c_.wells.size()) rep_.resize(c_.wells.size());
    stateIO(in);
    Xn_ = X_;
    // evalAll() refreshes the derived arrays from X; a continuous run keeps those of its last Newton evaluation (the
    // final update may follow it), which also seed the next flash (temperature guess) and the mechanics - keep the saved
    // ones so the continuation is bitwise identical (Y / dY are rebuilt at the next Newton iteration anyway)
    const std::vector<double> p = p_, T = T_, sw = sw_, so = so_, sg = sg_;
    const Stats st = st_;
    evalAll();
    p_ = p; T_ = T; sw_ = sw; so_ = so; sg_ = sg; st_ = st;
}


// ================================================================================================ C18 discrete adjoint
// Step residual (newtonStep, unscaled rows): component and energy rows  X_k - X_k^{n-1} + dt (...)  with the energy U a
// primary unknown, volume row PV(p, T) - sum V_alpha. Hence dR_n/dx_{n-1} = -I on the component and energy rows. J_n is
// the simulator's own Jacobian at the recorded converged state: chain rule over local cell properties whose derivatives
// are forward differences of the (p, N, U) flash - so the thermal adjoint is exact up to that differencing error (~1e-6
// relative per entry). dR/dtheta and dy/dtheta|_x: central differences of perturbed models at the recorded states; the
// initial-state dependence (pore volume -> moles / energy) through dR_1/dx_0. Lagged data (heat-loss memory, wellbore
// densities, surface yields, control modes and connection activity) are frozen at the recorded values.
namespace {
void mergeThermWells(std::vector<Well>& dst, const std::vector<Well>& src) {
    bool same = dst.size() == src.size();
    for (size_t w = 0; same && w < src.size(); ++w) same = dst[w].cells == src[w].cells && dst[w].wi.size() == src[w].wi.size();
    if (!same) { dst = src; return; }
    for (size_t w = 0; w < src.size(); ++w) { std::vector<double> wi = dst[w].wi; dst[w] = src[w]; dst[w].wi = wi; }
}
std::vector<int> regionOrAllT(const std::vector<int>& cells, int n) {
    if (!cells.empty()) return cells;
    std::vector<int> a(static_cast<size_t>(n)); for (int i = 0; i < n; ++i) a[size_t(i)] = i; return a;
}
void expandToT(std::vector<double>& v, int n) { if (v.size() == 1 && n > 1) v.assign(size_t(n), v[0]); }
struct TObsMap { int ka = 1, kb = 1; double wa = 1, wb = 0; };
}  // namespace

AdjParam paramPermMult(const std::string& name, const std::vector<int>& cells) {
    AdjParam P; P.name = name;
    P.apply = [cells](Case& c, double m) {
        const int n = c.size(); const std::vector<int> reg = regionOrAllT(cells, n);
        if (c.cellPV.empty()) {
            expandToT(c.permx, n); expandToT(c.permy, n); expandToT(c.permz, n);
            for (int q : reg) if (q >= 0 && q < n) { c.permx[size_t(q)] *= m; c.permy[size_t(q)] *= m; c.permz[size_t(q)] *= m; }
        } else {
            std::vector<char> in(static_cast<size_t>(n), 0); for (int q : reg) if (q >= 0 && q < n) in[size_t(q)] = 1;
            for (auto& cn : c.conns) if ((cn.a >= 0 && cn.a < n && in[size_t(cn.a)]) || (cn.b >= 0 && cn.b < n && in[size_t(cn.b)])) cn.T *= m;
        }
    };
    return P;
}
AdjParam paramPoroMult(const std::string& name, const std::vector<int>& cells) {
    AdjParam P; P.name = name;
    P.apply = [cells](Case& c, double m) {
        const int n = c.size(); const std::vector<int> reg = regionOrAllT(cells, n);
        if (c.cellPV.empty()) { expandToT(c.poro, n); for (int q : reg) if (q >= 0 && q < n) c.poro[size_t(q)] *= m; }
        else for (int q : reg) if (q >= 0 && q < n) c.cellPV[size_t(q)] *= m;
    };
    return P;
}
AdjParam paramWellIndexMult(const std::string& name, int well) {
    AdjParam P; P.name = name;
    P.apply = [well](Case& c, double m) { if (well >= 0 && size_t(well) < c.wells.size()) for (double& v : c.wells[size_t(well)].wi) v *= m; };
    return P;
}
AdjParam paramConductivityMult(const std::string& name) {
    AdjParam P; P.name = name;
    P.apply = [](Case& c, double m) { c.lamRock *= m; c.lamW *= m; c.lamO *= m; c.lamG *= m; for (double& v : c.thconr) v *= m; };
    return P;
}

void Simulator::setAdjointRecording(bool on) {
    adjOn_ = on; adjRec_.clear();
    if (on) { AdjRec r; r.t0 = r.t = t_; r.X = X_; r.T = T_; adjRec_.push_back(std::move(r)); }
}
std::vector<double> Simulator::stepSizes() const {
    std::vector<double> d; for (size_t k = 1; k < adjRec_.size(); ++k) d.push_back(adjRec_[k].dt); return d;
}
void Simulator::adjRecord(double t0, double dt) {
    AdjRec r; r.t0 = t0; r.t = t0 + dt; r.dt = dt;
    r.X = X_; r.T = T_; r.sgHist = sgHist_; r.permMult = permMult_;
    r.pvA = pvA_; r.pvCp = pvCp_; r.pvP = pvP_; r.pvCt = pvCt_; r.pvT = pvT_;
    r.connT.resize(conns_.size()); for (size_t e = 0; e < conns_.size(); ++e) r.connT[e] = conns_[e].T;
    r.fracW = fracW_; r.wbRho = wbRho_; r.wellRamp = wellRamp_; r.loss = loss_; r.surfYield = surfYield_; r.wellAct = wellAct_;
    r.wellMode = wellMode_; r.steamMode = steamMode_; r.rep = rep_; r.wells = c_.wells;
    adjRec_.push_back(std::move(r));
}
void Simulator::adjEval(const AdjRec& prev, const AdjRec& cur, const Simulator& base, std::vector<double>& R, BSR* J) {
    X_ = cur.X; Xn_ = prev.X; T_ = cur.T; t_ = cur.t0;
    sgHist_ = cur.sgHist; permMult_ = cur.permMult;
    pvP_ = cur.pvP; pvT_ = cur.pvT;
    pvA_.resize(size_t(n_)); pvCp_.resize(size_t(n_)); pvCt_.resize(size_t(n_));
    for (int q = 0; q < n_; ++q) {   // the recorded pore-volume law scaled to this model's reference pore volume
        const double b0 = base.pv0_[size_t(q)], f = b0 > 0 ? pv0_[size_t(q)] / b0 : 1.0;
        pvA_[size_t(q)] = cur.pvA[size_t(q)] * f; pvCp_[size_t(q)] = cur.pvCp[size_t(q)] * f; pvCt_[size_t(q)] = cur.pvCt[size_t(q)] * f;
    }
    for (size_t e = 0; e < conns_.size(); ++e)
        conns_[e].T = e < base.connT0_.size() && base.connT0_[e] > 0 ? cur.connT[e] * connT0_[e] / base.connT0_[e] : connT0_[e];
    fracW_ = cur.fracW; wbRho_ = cur.wbRho; wellRamp_ = cur.wellRamp; loss_ = cur.loss; surfYield_ = cur.surfYield; wellAct_ = cur.wellAct;
    wellMode_ = cur.wellMode; steamMode_ = cur.steamMode; rep_ = cur.rep;
    mergeThermWells(c_.wells, cur.wells);
    assembleOnly_ = true;
    bool conv = false; double res = 0;
    const bool ok = evalAll() && newtonStep(cur.dt, conv, res);
    assembleOnly_ = false;
    if (!ok) throw std::runtime_error("thermal adjoint: property evaluation failed at a recorded state");
    R = lastR_;
    if (J) *J = J_;
}
double Simulator::adjObsValue(const AdjObs& o, const WellReport& r) const {
    switch (o.qty) {
        case 3: return r.bhp;
        case 4: return r.hRate;
        case 6: return r.tFluid;
        default: break;
    }
    if (o.qty >= 10) { const size_t k = size_t(o.qty - 10); return k < r.molar.size() ? r.molar[k] : 0.0; }
    if (r.molar.empty() || (size_t(o.well) < c_.wells.size() && c_.wells[size_t(o.well)].injector)) return 0.0;
    double vo, vw, vg; surfaceVolumes(r.molar, vo, vw, vg);
    return o.qty == 0 ? vo : o.qty == 1 ? vg : o.qty == 2 ? vw : 0.0;
}

double Simulator::adjointMisfit(const std::vector<AdjObs>& obs) {
    AdjResult r = adjointGradient(obs, {});
    return r.F;
}

AdjResult Simulator::adjointGradient(const std::vector<AdjObs>& obs, const std::vector<AdjParam>& params) {
    AdjResult res;
    const auto tw0 = std::chrono::steady_clock::now();
    const int N = int(adjRec_.size()) - 1, np = int(params.size()), nb = nb_, nk = nc_ + 1, eE = nk, vU = nb - 1, m = n_ * nb;
    if (N < 1) { res.error = "no recorded steps (setAdjointRecording before the run)"; return res; }
    const std::string saved = saveState();
    const std::vector<Well> wSaved = c_.wells;
    const std::vector<double> XnS = Xn_;
    try {
        auto mapObs = [&](double t) {
            TObsMap mm;
            if (t <= adjRec_[1].t) return mm;
            if (t >= adjRec_[size_t(N)].t) { mm.ka = mm.kb = N; return mm; }
            int k = 1; while (k < N && adjRec_[size_t(k + 1)].t < t) ++k;
            const double a = (t - adjRec_[size_t(k)].t) / std::max(1e-300, adjRec_[size_t(k + 1)].t - adjRec_[size_t(k)].t);
            mm.ka = k; mm.kb = k + 1; mm.wa = 1 - a; mm.wb = a; return mm;
        };
        std::vector<double> R;
        std::map<int, std::vector<WellReport>> repAt;
        auto reportsAt = [&](int k) -> const std::vector<WellReport>& {
            auto it = repAt.find(k);
            if (it != repAt.end()) return it->second;
            adjEval(adjRec_[size_t(k - 1)], adjRec_[size_t(k)], *this, R, nullptr);
            return repAt[k] = rep_;
        };
        std::vector<double> dFdy(obs.size(), 0.0);
        std::vector<std::vector<std::pair<int, double>>> atStep(size_t(N) + 1);
        for (size_t j = 0; j < obs.size(); ++j) {
            const TObsMap om = mapObs(obs[j].t);
            double y = 0;
            for (int side = 0; side < 2; ++side) {
                const int k = side ? om.kb : om.ka; const double w = side ? om.wb : om.wa;
                if (w == 0) continue;
                const auto& rp = reportsAt(k);
                y += w * (size_t(obs[j].well) < rp.size() ? adjObsValue(obs[j], rp[size_t(obs[j].well)]) : 0.0);
                atStep[size_t(k)].push_back({int(j), w});
            }
            const double r = (y - obs[j].value) / obs[j].sigma;
            res.F += r * r; dFdy[j] = 2 * r / obs[j].sigma;
        }
        res.grad.assign(size_t(np), 0.0);
        if (np > 0) {
            Case base = c_; base.wells = adjRec_[1].wells;
            std::vector<std::unique_ptr<Simulator>> sp(static_cast<size_t>(np)), sm(static_cast<size_t>(np));
            std::vector<std::vector<double>> dX0(static_cast<size_t>(np));
            for (int j = 0; j < np; ++j) {
                const double h = params[size_t(j)].h;
                Case cp = base, cm = base;
                params[size_t(j)].apply(cp, 1 + h); params[size_t(j)].apply(cm, 1 - h);
                sp[size_t(j)].reset(new Simulator()); sm[size_t(j)].reset(new Simulator());
                std::string e1, e2;
                if (!sp[size_t(j)]->setup(cp, &e1) || !sm[size_t(j)]->setup(cm, &e2)) throw std::runtime_error("thermal adjoint: perturbed model setup failed: " + e1 + e2);
                dX0[size_t(j)].resize(X_.size());
                for (size_t i = 0; i < X_.size(); ++i) dX0[size_t(j)][i] = (sp[size_t(j)]->X_[i] - sm[size_t(j)]->X_[i]) / (2 * h);
            }
            // explicit dependence of the observations on theta
            for (int k = 1; k <= N; ++k) {
                if (atStep[size_t(k)].empty()) continue;
                for (int j = 0; j < np; ++j) {
                    sp[size_t(j)]->adjEval(adjRec_[size_t(k - 1)], adjRec_[size_t(k)], *this, R, nullptr);
                    sm[size_t(j)]->adjEval(adjRec_[size_t(k - 1)], adjRec_[size_t(k)], *this, R, nullptr);
                    for (const auto& ow : atStep[size_t(k)]) {
                        const AdjObs& o = obs[size_t(ow.first)];
                        if (size_t(o.well) >= sp[size_t(j)]->rep_.size()) continue;
                        const double dy = (sp[size_t(j)]->adjObsValue(o, sp[size_t(j)]->rep_[size_t(o.well)]) -
                                           sm[size_t(j)]->adjObsValue(o, sm[size_t(j)]->rep_[size_t(o.well)])) / (2 * params[size_t(j)].h);
                        res.grad[size_t(j)] += dFdy[size_t(ow.first)] * ow.second * dy;
                    }
                }
            }
            StepAdjointProblem P;
            P.steps = N;
            P.jacobian = [&](int n, BSR& J) { adjEval(adjRec_[size_t(n - 1)], adjRec_[size_t(n)], *this, R, &J); };
            P.addObjectiveGradient = [&](int n, std::vector<double>& rhs) {
                const auto& list = atStep[size_t(n)];
                if (list.empty()) return;
                std::vector<char> seen(static_cast<size_t>(n_), 0); std::vector<int> cells;
                for (const auto& ow : list) {
                    const int w = obs[size_t(ow.first)].well;
                    if (w < 0 || size_t(w) >= adjRec_[size_t(n)].wells.size()) continue;
                    for (int q : adjRec_[size_t(n)].wells[size_t(w)].cells) if (q >= 0 && q < n_ && !seen[size_t(q)]) { seen[size_t(q)] = 1; cells.push_back(q); }
                }
                AdjRec tmp = adjRec_[size_t(n)];
                for (int q : cells) {
                    double Nt = 0; for (int k = 1; k < nb - 1; ++k) Nt += std::max(tmp.X[size_t(q) * nb + k], 0.0);
                    for (int v = 0; v < nb; ++v) {
                        double& x = tmp.X[size_t(q) * nb + v];
                        const double x0 = x;
                        const double d = v == 0 ? 1e-6 * std::max(std::fabs(x0), 1e5) : v == vU ? 1e-6 * (crock_[size_t(q)] + 80 * Nt) : 1e-6 * std::max(Nt, 1e-30);
                        x = x0 + d; adjEval(adjRec_[size_t(n - 1)], tmp, *this, R, nullptr); const std::vector<WellReport> rp = rep_;
                        x = x0 - d; adjEval(adjRec_[size_t(n - 1)], tmp, *this, R, nullptr); const std::vector<WellReport> rm = rep_;
                        x = x0;
                        double g = 0;
                        for (const auto& ow : list) {
                            const AdjObs& o = obs[size_t(ow.first)];
                            if (size_t(o.well) >= rp.size()) continue;
                            g += dFdy[size_t(ow.first)] * ow.second * (adjObsValue(o, rp[size_t(o.well)]) - adjObsValue(o, rm[size_t(o.well)])) / (2 * d);
                        }
                        rhs[size_t(q) * nb + size_t(v)] += g;
                    }
                }
            };
            P.addPrevCoupling = [&](int, const std::vector<double>& lam, std::vector<double>& rhs) {
                for (int q = 0; q < n_; ++q) {
                    for (int k = 0; k < nk; ++k) rhs[size_t(q) * nb + 1 + k] -= lam[size_t(q) * nb + k];
                    rhs[size_t(q) * nb + vU] -= lam[size_t(q) * nb + eE];
                }
            };
            P.addParamTerm = [&](int n, const std::vector<double>& lam, std::vector<double>& grad) {
                std::vector<double> Rp, Rm;
                for (int j = 0; j < np; ++j) {
                    sp[size_t(j)]->adjEval(adjRec_[size_t(n - 1)], adjRec_[size_t(n)], *this, Rp, nullptr);
                    sm[size_t(j)]->adjEval(adjRec_[size_t(n - 1)], adjRec_[size_t(n)], *this, Rm, nullptr);
                    double s = 0; for (int r = 0; r < m; ++r) s += lam[size_t(r)] * (Rp[size_t(r)] - Rm[size_t(r)]);
                    grad[size_t(j)] += s / (2 * params[size_t(j)].h);
                }
            };
            P.addInitialTerm = [&](const std::vector<double>& lam1, std::vector<double>& grad) {
                for (int j = 0; j < np; ++j) {
                    double s = 0;
                    for (int q = 0; q < n_; ++q) {
                        for (int k = 0; k < nk; ++k) s -= lam1[size_t(q) * nb + k] * dX0[size_t(j)][size_t(q) * nb + 1 + k];
                        s -= lam1[size_t(q) * nb + eE] * dX0[size_t(j)][size_t(q) * nb + vU];
                    }
                    grad[size_t(j)] += s;
                }
            };
            const StepAdjointStats st = solveStepAdjoint(P, res.grad);
            res.ok = st.ok; res.error = st.error; res.linear = st.linear;
        } else res.ok = true;
        res.steps = N;
    } catch (const std::exception& e) { res.ok = false; res.error = e.what(); }
    c_.wells = wSaved; loadState(saved); Xn_ = XnS;
    evalAll();
    res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - tw0).count();
    return res;
}

} // namespace thermal
} // namespace ressim
