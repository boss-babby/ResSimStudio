// C21 geomechanics G1-c: staggered EDFM flow <-> EFEM contact coupling (see geocouple.h)
#include "ressim/geocouple.h"
#include "ressim/simulator.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>

namespace ressim {
namespace geocouple {

double aperture(const BartonBandis& bb, double sigmaNeff, double opening, double slip) {
    const double Vm = std::max(bb.Vm, 0.0), Kn = std::max(bb.Kn0, 1e-30);
    double a = bb.aRes + Vm;
    if (opening > 0) a += opening;
    else if (sigmaNeff > 0) a -= Vm * sigmaNeff / (Kn * Vm + sigmaNeff);
    const double psi = bb.dilationDeg * 3.14159265358979323846 / 180.0;
    if (psi > 0 && slip > 0) a += std::min(std::fabs(slip) * std::tan(psi), std::max(bb.maxDilation, 0.0));
    return std::max(a, 1e-12);
}

namespace {
struct Geo {
    double toM = 1, toPa = 1e5, Lx = 0, Ly = 0;
    std::vector<double> s;            // per piece (m)
    std::vector<int> frac, layer;
};

bool geometry(const Simulator& sim, Geo& G, std::string& err) {
    const Case& c = sim.caseData();
    const GridInput& g = c.grid;
    const Units U = Units::of(c.units);
    G.toPa = 101325.0 / U.atm;
    G.toM = c.units == UnitSystem::Field ? 0.3048 : 1.0;
    const size_t ncell = size_t(g.nx) * g.ny * g.nz;
    auto dxAt = [&](const std::vector<double>& v, int i, int j, bool isX) -> double {
        if (v.size() == ncell) return v[size_t(g.index(i, j, 0))];
        if (isX && v.size() == size_t(g.nx)) return v[size_t(i)];
        if (!isX && v.size() == size_t(g.ny)) return v[size_t(j)];
        if (v.size() == 1) return v[0];
        return -1;
    };
    G.Lx = 0; G.Ly = 0;
    for (int i = 0; i < g.nx; ++i) { const double d = dxAt(g.dx, i, 0, true); if (d <= 0) { err = "G1-c coupling needs a Cartesian DX / DY grid"; return false; } G.Lx += d; }
    for (int j = 0; j < g.ny; ++j) { const double d = dxAt(g.dy, 0, j, false); if (d <= 0) { err = "G1-c coupling needs a Cartesian DX / DY grid"; return false; } G.Ly += d; }
    G.Lx *= G.toM; G.Ly *= G.toM;
    const int nH = sim.numHydraulicFractureCells();
    if (nH == 0 || int(g.hfPieces.size()) != nH) { err = "no hydraulic / embedded fracture pieces (HFRAC) in this case"; return false; }
    G.s.assign(size_t(nH), 0.0); G.frac.assign(size_t(nH), 0); G.layer.assign(size_t(nH), 0);
    for (int k = 0; k < nH; ++k) {
        const auto& P = g.hfPieces[size_t(k)];
        const auto& F = g.hfracs[size_t(P.frac)];
        const double L = std::hypot(F.x1 - F.x0, F.y1 - F.y0);
        const double mx = (F.x1 - F.x0) / std::max(L, 1e-12), my = (F.y1 - F.y0) / std::max(L, 1e-12);
        const double xc = 0.5 * (P.xa + P.xb) - 0.5 * (F.x0 + F.x1), yc = 0.5 * (P.ya + P.yb) - 0.5 * (F.y0 + F.y1);
        G.s[size_t(k)] = (xc * mx + yc * my) * G.toM;
        G.frac[size_t(k)] = P.frac;
        G.layer[size_t(k)] = P.host / std::max(1, g.nx * g.ny);
    }
    return true;
}

// linear interpolation of y(x) on sorted x (clamped)
double interp(const std::vector<std::pair<double, double>>& xy, double x) {
    if (xy.empty()) return 0;
    if (x <= xy.front().first) return xy.front().second;
    if (x >= xy.back().first) return xy.back().second;
    const auto it = std::lower_bound(xy.begin(), xy.end(), std::make_pair(x, -1e300));
    const auto& b = *it; const auto& a = *(it - 1);
    const double w = (x - a.first) / std::max(b.first - a.first, 1e-300);
    return a.second + w * (b.second - a.second);
}
}  // namespace

bool coupleOnce(Simulator& sim, const Input& in, efem::Result& prev, StepOut& st, Output* geo, std::string& err) {
    Geo G;
    if (!geometry(sim, G, err)) return false;
    const Case& c = sim.caseData();
    const GridInput& g = c.grid;
    const int nH = int(G.s.size()), nf = int(g.hfracs.size());
    if (geo && geo->s.empty()) { geo->s = G.s; geo->frac = G.frac; }
    efem::Input E;
    E.Lx = G.Lx; E.Ly = G.Ly; E.nx = std::max(4, in.nx); E.ny = std::max(4, in.ny);
    E.E = in.E; E.nu = in.nu; E.sxx = in.sxx; E.syy = in.syy; E.penalty = in.penalty;
    // an axis-parallel fracture lying on an element edge gives degenerate cuts: refine that direction by one element
    // (deterministic, so the slip history of the previous step stays compatible)
    for (int pass = 0; pass < 8; ++pass) {
        bool bump = false;
        const double hx = E.Lx / E.nx, hy = E.Ly / E.ny;
        for (const auto& F : g.hfracs) {
            const double y0 = F.y0 * G.toM, y1 = F.y1 * G.toM, x0 = F.x0 * G.toM, x1 = F.x1 * G.toM;
            if (std::fabs(y1 - y0) < 1e-9 * E.Ly) { const double r = std::fmod(y0, hy) / hy; if (std::min(r, 1 - r) < 0.02) { ++E.ny; bump = true; } }
            if (std::fabs(x1 - x0) < 1e-9 * E.Lx) { const double r = std::fmod(x0, hx) / hx; if (std::min(r, 1 - r) < 0.02) { ++E.nx; bump = true; } }
        }
        if (!bump) break;
    }
    st.p.assign(size_t(nH), 0.0);
    for (int k = 0; k < nH; ++k) st.p[size_t(k)] = sim.hfCellPressure(k) * G.toPa;
    for (int f = 0; f < nf; ++f) {
        const auto& F = g.hfracs[size_t(f)];
        efem::Fracture R;
        R.x0 = F.x0 * G.toM; R.y0 = F.y0 * G.toM; R.x1 = F.x1 * G.toM; R.y1 = F.y1 * G.toM;
        R.friction = in.friction; R.cohesion = in.cohesion;
        std::vector<std::pair<double, double>> prof;
        for (int k = 0; k < nH; ++k) if (G.frac[size_t(k)] == f && G.layer[size_t(k)] == in.layer) prof.push_back({G.s[size_t(k)], st.p[size_t(k)]});
        if (prof.empty()) for (int k = 0; k < nH; ++k) if (G.frac[size_t(k)] == f) prof.push_back({G.s[size_t(k)], st.p[size_t(k)]});
        std::sort(prof.begin(), prof.end());
        for (const auto& q : prof) { R.profS.push_back(q.first); R.profP.push_back(q.second); }
        if (R.profS.size() == 1) { R.pressure = R.profP[0]; R.profS.clear(); R.profP.clear(); }
        E.fracs.push_back(R);
    }
    efem::Result res;
    const bool havePrev = !prev.cut.empty();
    if (!efem::solve(E, res, havePrev ? &prev : nullptr) || !res.ok) { err = "contact solve failed: " + res.error; return false; }
    st.newton = res.newton;
    // per fracture: cut-element profiles of opening, slip, normal traction
    using Prof = std::vector<std::vector<std::pair<double, double>>>;
    Prof wn = Prof(size_t(nf)), wt = Prof(size_t(nf)), tn = Prof(size_t(nf));
    // the crack-end elements (tip singularity, partially cut) and slivers carry unreliable penalty tractions: drop them
    std::vector<double> sMin(size_t(nf), 1e300), sMax(size_t(nf), -1e300), lenMax(size_t(nf), 0.0);
    std::vector<int> ncut(size_t(nf), 0);
    for (const auto& C : res.cut) if (C.frac >= 0 && C.frac < nf) {
        sMin[size_t(C.frac)] = std::min(sMin[size_t(C.frac)], C.s); sMax[size_t(C.frac)] = std::max(sMax[size_t(C.frac)], C.s);
        lenMax[size_t(C.frac)] = std::max(lenMax[size_t(C.frac)], C.len); ++ncut[size_t(C.frac)];
    }
    for (const auto& C : res.cut) {
        if (C.frac < 0 || C.frac >= nf) continue;
        if (ncut[size_t(C.frac)] >= 4 && (C.s <= sMin[size_t(C.frac)] || C.s >= sMax[size_t(C.frac)] || C.len < 0.25 * lenMax[size_t(C.frac)])) continue;
        wn[size_t(C.frac)].push_back({C.s, C.wn});
        wt[size_t(C.frac)].push_back({C.s, std::max(std::fabs(C.wtp), C.state == 2 ? std::fabs(C.wt) : 0.0)});
        tn[size_t(C.frac)].push_back({C.s, C.tn});
    }
    // finite-height (PKN) correction of the map-view opening
    std::vector<double> wfac(size_t(nf), 1.0);
    if (in.height != 0)
        for (int f = 0; f < nf; ++f) {
            const auto& F = g.hfracs[size_t(f)];
            const double Lf = std::hypot(F.x1 - F.x0, F.y1 - F.y0) * G.toM;
            double H = in.height;
            if (H < 0) {
                std::vector<double> hl(size_t(std::max(1, g.nz)), 0.0);
                for (int k = 0; k < nH; ++k) if (G.frac[size_t(k)] == f) {
                    const auto& P = g.hfPieces[size_t(k)];
                    const size_t l = size_t(std::clamp(G.layer[size_t(k)], 0, std::max(0, g.nz - 1)));
                    hl[l] = std::max(hl[l], std::fabs(P.zb - P.zt) * G.toM);
                }
                H = std::accumulate(hl.begin(), hl.end(), 0.0);
            }
            if (H > 0 && Lf > 0) wfac[size_t(f)] = std::min(1.0, H / Lf);
        }
    for (int f = 0; f < nf; ++f) { std::sort(wn[size_t(f)].begin(), wn[size_t(f)].end()); std::sort(wt[size_t(f)].begin(), wt[size_t(f)].end()); std::sort(tn[size_t(f)].begin(), tn[size_t(f)].end()); }
    st.wn.assign(size_t(nH), 0.0); st.wt.assign(size_t(nH), 0.0); st.tn.assign(size_t(nH), 0.0); st.a.assign(size_t(nH), 0.0); st.mult.assign(size_t(nH), 1.0);
    st.maxOpening = 0; st.maxSlip = 0; st.minMult = 1e300; st.maxMult = 0; st.openPieces = 0; st.slipPieces = 0;
    for (int k = 0; k < nH; ++k) {
        const int f = G.frac[size_t(k)];
        if (wn[size_t(f)].empty()) continue;
        const double s = G.s[size_t(k)];
        const double w = interp(wn[size_t(f)], s) * wfac[size_t(f)], sl = interp(wt[size_t(f)], s), t = interp(tn[size_t(f)], s);
        const double a = aperture(in.bb, std::max(0.0, -t), std::max(0.0, w), sl);
        const double aRef = in.aRef > 0 ? in.aRef : std::max(g.hfracs[size_t(f)].width * G.toM, 1e-9);
        const double m = std::clamp(std::pow(a / aRef, 3.0), 1e-6, 1e6);
        st.wn[size_t(k)] = w; st.wt[size_t(k)] = sl; st.tn[size_t(k)] = t; st.a[size_t(k)] = a; st.mult[size_t(k)] = m;
        st.maxOpening = std::max(st.maxOpening, w); st.maxSlip = std::max(st.maxSlip, sl);
        st.minMult = std::min(st.minMult, m); st.maxMult = std::max(st.maxMult, m);
        if (w > 0) ++st.openPieces;
        if (sl > 0) ++st.slipPieces;
    }
    if (st.minMult > st.maxMult) st.minMult = st.maxMult = 1;
    sim.setHfConductance(st.mult);
    st.t = sim.time();
    prev = std::move(res);
    return true;
}

bool runCoupled(Simulator& sim, const Input& in, Output& out, std::atomic<bool>* cancel, const std::function<void(const std::string&)>& log) {
    out = Output();
    efem::Result prev;
    auto t0 = std::chrono::steady_clock::now();
    double efemSec = 0;
    auto once = [&]() -> bool {
        StepOut st; std::string err;
        const auto a = std::chrono::steady_clock::now();
        if (!coupleOnce(sim, in, prev, st, &out, err)) { out.error = err; return false; }
        efemSec += std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
        if (log) {
            char b[256];
            std::snprintf(b, sizeof b, "  G1-c t = %.4g d: %d open / %d slipping pieces, max opening %.3g mm, max slip %.3g mm, conductance x %.3g .. %.3g",
                          st.t, st.openPieces, st.slipPieces, st.maxOpening * 1e3, st.maxSlip * 1e3, st.minMult, st.maxMult);
            log(b);
        }
        out.steps.push_back(std::move(st));
        return true;
    };
    try {
        if (!once()) return false;
        while (sim.runReportStep()) {
            if (cancel && cancel->load()) { out.error = "cancelled"; return false; }
            if (!once()) return false;
        }
    } catch (const std::exception& e) { out.error = e.what(); return false; }
    out.efemSeconds = efemSec;
    (void)t0;
    out.ok = true;
    return true;
}

}  // namespace geocouple
}  // namespace ressim
