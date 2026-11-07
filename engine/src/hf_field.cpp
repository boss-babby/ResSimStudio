// C30: one-way stress / opening field of the in-run hydraulic fractures - see hf_field.h
#include "ressim/hf_field.h"
#include "ressim/case.h"
#include "ressim/simulator.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>

namespace ressim {
namespace hffield {

namespace {
inline double dot3(const ddm::V3& a, const ddm::V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
// distance (case units) from a point to the vertical rectangle of a piece; inside: true when the projection falls on it
double pieceDistance(const GridInput::HfPiece& P, double x, double y, double z, double& normalDist, bool& inside) {
    const double tx = P.xb - P.xa, ty = P.yb - P.ya, L2 = tx * tx + ty * ty;
    double t = L2 > 0 ? ((x - P.xa) * tx + (y - P.ya) * ty) / L2 : 0.0;
    const double tc = std::clamp(t, 0.0, 1.0);
    const double px = P.xa + tc * tx - x, py = P.ya + tc * ty - y;
    const double dz = std::max({0.0, P.zt - z, z - P.zb});
    const double L = std::sqrt(L2);
    normalDist = L > 0 ? std::fabs((x - P.xa) * ty - (y - P.ya) * tx) / L : std::hypot(px, py);
    inside = t >= -1e-9 && t <= 1 + 1e-9 && dz == 0;
    return std::sqrt(px * px + py * py + dz * dz);
}
}  // namespace

bool Field::setup(const Simulator& sim, std::string& note) {
    active_ = have_ = false; st_ = Stats{}; pieces_.clear(); el_.clear(); slipDir_.clear(); Ko_.clear(); zMid_.clear(); obs_.clear();
    const Case& c = sim.caseData();
    const Grid& g = sim.grid();
    const auto& H = c.grid.hfracs;
    const int nH = g.numHfCells();
    if (c.hfField.on == 0) { note = "HFFIELD OFF"; return false; }
    if (nH == 0 || c.grid.hfPieces.size() != size_t(nH)) { note = "no fracture pieces"; return false; }
    // eligible fractures: a closure stress is defined (HFGROW / DDMWIDTH, or HFFIELD ON with sigmaMin > 0)
    std::vector<char> elig(H.size(), 0);
    bool any = false;
    for (size_t f = 0; f < H.size(); ++f) {
        const auto& F = H[f];
        // phase-field fractures carry their own opening / stress (PFWIDTH / PFDS*): the DDM field from the flow pressure is a
        // different crack model, included only on request (HFFIELD ON)
        elig[f] = ((F.kic > 0 && !F.phaseField) || F.ddmWidth || (c.hfField.on == 1 && (F.sigma > 0 || F.phaseField))) ? 1 : 0;
        any = any || elig[f];
    }
    if (!any) { note = "no fracture with a closure stress (HFGROW / DDMWIDTH)"; return false; }
    const bool fld = c.units == UnitSystem::Field;
    Mlen_ = fld ? 0.3048 : 1.0; Pa_ = fld ? 6894.757 : 1e5;
    double E = 20e9, nu = 0.25;
    if (c.hfField.E > 0) { E = c.hfField.E * 1e9; nu = c.hfField.nu; }
    else {
        bool got = false;
        for (const auto& F : H) if (F.ddmWidth && F.ddmE > 0) { E = F.ddmE * 1e9; nu = F.ddmNu; got = true; break; }
        if (!got) for (const auto& F : H) if (F.phaseField && F.pfE > 0) { E = F.pfE * 1e9; nu = F.pfNu; got = true; break; }
        if (!got && c.geoSeq.enabled && !c.geoSeq.E.empty() && c.geoSeq.E[0] > 0) { E = c.geoSeq.E[0]; nu = c.geoSeq.nu.empty() ? 0.25 : c.geoSeq.nu[0]; }
    }
    st_.E = E; st_.nu = nu; nu_ = nu; G_ = E / (2 * (1 + nu));
    // elements, extent, depth
    double ext = 0, Hmax = 0, zMin = 1e300;
    std::vector<double> sMin(H.size(), 1e300), sMax(H.size(), -1e300);
    for (int k = 0; k < nH; ++k) {
        const auto& P = c.grid.hfPieces[size_t(k)];
        if (!elig[size_t(P.frac)]) continue;
        const ddm::V3 a = {P.xa * Mlen_, P.ya * Mlen_, -P.zt * Mlen_}, b = {P.xb * Mlen_, P.yb * Mlen_, -P.zt * Mlen_};
        const ddm::V3 cc = {P.xb * Mlen_, P.yb * Mlen_, -P.zb * Mlen_}, d = {P.xa * Mlen_, P.ya * Mlen_, -P.zb * Mlen_};
        ddm::Element e = ddm::quadElement(a, b, cc, d);
        if (e.tris.size() != 2 || !(e.area > 0)) continue;
        pieces_.push_back(k);
        zMid_.push_back(0.5 * (P.zt + P.zb));
        el_.push_back(e);
        zMin = std::min(zMin, P.zt * Mlen_);
        Hmax = std::max(Hmax, std::fabs(P.zb - P.zt) * Mlen_);
    }
    for (size_t f = 0; f < H.size(); ++f) if (elig[f]) {
        const auto& F = H[f];
        ext = std::max(ext, 0.5 * std::hypot(F.x1 - F.x0, F.y1 - F.y0) * Mlen_);
        if (F.zTop >= 0 && F.zBot > F.zTop) Hmax = std::max(Hmax, (F.zBot - F.zTop) * Mlen_);
    }
    const int n = int(pieces_.size());
    st_.pieces = n;
    if (n == 0) { note = "no fracture elements"; return false; }
    if (n > 2500) { note = "more than 2500 fracture pieces - field output skipped (use HFFIELD OFF to silence)"; return false; }
    const double extent = std::max(ext, 0.5 * Hmax);
    st_.halfSpace = zMin < 2 * std::max(extent, Hmax);
    const auto t0 = std::chrono::steady_clock::now();
    Ko_ = ddm::buildOpening(el_, G_, nu_, st_.halfSpace);
    {   // opening sense of buildOpening (it flips K when the +n discontinuity tensions its own element)
        double S[6];
        sense_ = 1;
        if (ddm::elementStress(el_[0].c, el_[0], el_[0].n, G_, nu_, st_.halfSpace, S)) {
            const ddm::V3& m = el_[0].n;
            const double tn = m[0] * (S[0] * m[0] + S[3] * m[1] + S[4] * m[2]) + m[1] * (S[3] * m[0] + S[1] * m[1] + S[5] * m[2]) + m[2] * (S[4] * m[0] + S[5] * m[1] + S[2] * m[2]);
            if (tn > 0) sense_ = -1;
        }
    }
    slipDir_.resize(size_t(n));
    for (int q = 0; q < n; ++q)
        for (int t = 0; t < 2; ++t) {
            ddm::V3 nn, s, d; ddm::triFrame(el_[size_t(q)].tris[size_t(t)], nn, s, d);
            const ddm::V3& b = el_[size_t(q)].n;
            slipDir_[size_t(q)][size_t(t)] = {sense_ * dot3(b, s), sense_ * dot3(b, d), sense_ * dot3(b, nn)};
        }
    // observation cells: active matrix cells within the radius, nearest first, capped by the budget
    const double R = c.hfField.radius * std::max(extent, 1.0);
    st_.radiusM = R;
    struct Cand { double d; Obs o; };
    std::vector<Cand> cand;
    const int nc = g.numCells();
    if (g.corners.size() >= size_t(nc)) {
        double bx0 = 1e300, bx1 = -1e300, by0 = 1e300, by1 = -1e300, bz0 = 1e300, bz1 = -1e300;
        for (int k : pieces_) {
            const auto& P = c.grid.hfPieces[size_t(k)];
            bx0 = std::min({bx0, P.xa, P.xb}); bx1 = std::max({bx1, P.xa, P.xb}); by0 = std::min({by0, P.ya, P.yb}); by1 = std::max({by1, P.ya, P.yb});
            bz0 = std::min(bz0, P.zt); bz1 = std::max(bz1, P.zb);
        }
        const double Rc = R / Mlen_;
        for (int gq = 0; gq < nc; ++gq) {
            const int a = g.activeOf[size_t(gq)];
            if (a < 0) continue;
            const auto& C = g.corners[size_t(gq)];
            double x = 0, y = 0, z = 0;
            for (int v = 0; v < 8; ++v) { x += C[size_t(v)].x; y += C[size_t(v)].y; z += C[size_t(v)].z; }
            x /= 8; y /= 8; z /= 8;
            if (x < bx0 - Rc || x > bx1 + Rc || y < by0 - Rc || y > by1 + Rc || z < bz0 - Rc || z > bz1 + Rc) continue;
            double best = 1e300, bestN = 0; int bq = -1; bool bIn = false;
            for (int q = 0; q < n; ++q) {
                double dn; bool in;
                const double dd = pieceDistance(c.grid.hfPieces[size_t(pieces_[size_t(q)])], x, y, z, dn, in);
                if (dd < best) { best = dd; bq = q; bestN = dn; bIn = in; }
            }
            if (bq < 0 || best > Rc) continue;
            Cand cd; cd.d = best; cd.o.act = a; cd.o.n = el_[size_t(bq)].n;
            cd.o.x = {x * Mlen_, y * Mlen_, -z * Mlen_};
            const double cw = std::max(1e-6, std::min(gq < int(g.cellDx.size()) ? g.cellDx[size_t(gq)] : 1.0, gq < int(g.cellDy.size()) ? g.cellDy[size_t(gq)] : 1.0));
            if (bIn && bestN < 0.05 * cw) {   // centre in the fracture plane: average the two faces at +-cw / 4
                const double off = 0.25 * cw * Mlen_;
                cd.o.twin = true;
                cd.o.x2 = {cd.o.x[0] - off * cd.o.n[0], cd.o.x[1] - off * cd.o.n[1], cd.o.x[2] - off * cd.o.n[2]};
                cd.o.x = {cd.o.x[0] + off * cd.o.n[0], cd.o.x[1] + off * cd.o.n[1], cd.o.x[2] + off * cd.o.n[2]};
            }
            cand.push_back(cd);
        }
    }
    std::sort(cand.begin(), cand.end(), [](const Cand& A, const Cand& B) { return A.d < B.d; });
    st_.cellsInRadius = int(cand.size());
    const double perCell = 2.0 * 2.0 * double(n) * (st_.halfSpace ? 7.0 : 1.0);   // triangles x (stress + displacement), half space ~7x
    const size_t cap = size_t(std::max(1.0, c.hfField.budget * 1e6 / perCell));
    double twinCost = 0;
    for (const auto& cd : cand) {
        twinCost += cd.o.twin ? 2 : 1;
        if (twinCost > double(cap)) break;
        obs_.push_back(cd.o);
    }
    st_.cells = int(obs_.size());
    st_.setupSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (st_.cells < st_.cellsInRadius)
        st_.note = "stress field limited to the " + std::to_string(st_.cells) + " nearest of " + std::to_string(st_.cellsInRadius) +
                   " cells in the radius (raise the HFFIELD budget)";
    w_.assign(size_t(n), 0.0); wPrev_.assign(size_t(n), -1.0);
    wAll_.assign(size_t(nH), 0.0); pnAll_.assign(size_t(nH), 0.0);
    active_ = true;
    note = st_.note;
    return true;
}

void Field::evalPoint(const ddm::V3& x, double S[6], double u[3]) const {
    for (int i = 0; i < 6; ++i) S[i] = 0;
    u[0] = u[1] = u[2] = 0;
    ddm::V3 o = x;
    if (st_.halfSpace && o[2] > 0) o[2] = 0;
    for (size_t q = 0; q < el_.size(); ++q) {
        const double w = w_[q];
        if (!(w > 0)) continue;
        const ddm::Element& e = el_[q];
        // near the element plane (inside 1 % of its size): the robust +-delta average of ddm::elementStress
        const ddm::V3 r = {o[0] - e.c[0], o[1] - e.c[1], o[2] - e.c[2]};
        const double dn = std::fabs(dot3(r, e.n)), rr = std::sqrt(dot3(r, r));
        if (dn < 1e-2 * e.size && rr < 1.5 * e.size) {
            double St[6];
            const ddm::V3 b = {sense_ * w * e.n[0], sense_ * w * e.n[1], sense_ * w * e.n[2]};
            if (ddm::elementStress(o, e, b, G_, nu_, st_.halfSpace, St)) for (int i = 0; i < 6; ++i) S[i] += St[i];
            for (int t = 0; t < 2; ++t) {
                const auto& d = slipDir_[q][size_t(t)];
                double ut[3];
                if (ddm::triDisp(o, e.tris[size_t(t)], {w * d[0], w * d[1], w * d[2]}, nu_, st_.halfSpace, ut)) for (int i = 0; i < 3; ++i) u[i] += ut[i];
            }
            continue;
        }
        for (int t = 0; t < 2; ++t) {
            const auto& d = slipDir_[q][size_t(t)];
            const ddm::V3 sl = {w * d[0], w * d[1], w * d[2]};
            double St[6], ut[3];
            if (ddm::triStress(o, e.tris[size_t(t)], sl, G_, nu_, st_.halfSpace, St)) for (int i = 0; i < 6; ++i) S[i] += St[i];
            if (ddm::triDisp(o, e.tris[size_t(t)], sl, nu_, st_.halfSpace, ut)) for (int i = 0; i < 3; ++i) u[i] += ut[i];
        }
    }
}

void Field::fieldAt(const double x[3], double S[6], double u[3]) const { evalPoint({x[0], x[1], x[2]}, S, u); }

void Field::addArrays(const Simulator& sim, Snapshot& s) {
    if (!active_) return;
    const auto t0 = std::chrono::steady_clock::now();
    const Case& c = sim.caseData();
    const Grid& g = sim.grid();
    const int n = int(pieces_.size()), nH = g.numHfCells();
    const size_t na = s.pressure.size();
    const auto& sMin = sim.minStress();
    std::vector<int> idx;
    std::vector<double> pn(size_t(n), 0.0);
    for (int q = 0; q < n; ++q) {
        const int k = pieces_[size_t(q)];
        const auto& P = c.grid.hfPieces[size_t(k)];
        const auto& F = c.grid.hfracs[size_t(P.frac)];
        double sig = 0; bool feOk = false;
        if (F.feStress && !sMin.empty()) {
            const int host = g.hfCells[size_t(k)].host;
            const int a = host >= 0 && size_t(host) < g.activeOf.size() ? g.activeOf[size_t(host)] : -1;
            if (a >= 0 && size_t(a) < sMin.size() && std::isfinite(sMin[size_t(a)]) && sMin[size_t(a)] > 0) { sig = sMin[size_t(a)]; feOk = true; }
        }
        if (!feOk) {
            const double zm = 0.5 * (F.zTop + F.zBot) > 0 ? 0.5 * (F.zTop + F.zBot) : zMid_[size_t(q)];
            sig = F.sigma + F.sigmaGrad * (zMid_[size_t(q)] - zm);
        }
        const double pnet = sim.hfCellPressure(k) - sig;
        pnAll_[size_t(k)] = pnet;
        if (sim.hfCellOpenFraction(k) > 0) { idx.push_back(q); pn[size_t(q)] = pnet * Pa_; }
    }
    const int m = int(idx.size());
    std::fill(w_.begin(), w_.end(), 0.0);
    if (m > 0) {
        std::vector<double> K(static_cast<size_t>(m) * static_cast<size_t>(m)), p(static_cast<size_t>(m)), ws;
        for (int r = 0; r < m; ++r) {
            p[size_t(r)] = pn[size_t(idx[size_t(r)])];
            for (int q = 0; q < m; ++q) K[size_t(r) * size_t(m) + size_t(q)] = Ko_[size_t(idx[size_t(r)]) * size_t(n) + size_t(idx[size_t(q)])];
        }
        ddm::solveOpening(K, m, p, ws, 1.0);
        for (int r = 0; r < m; ++r) w_[size_t(idx[size_t(r)])] = std::max(0.0, ws[size_t(r)]);
    }
    st_.open = 0; st_.maxWidth = 0;
    for (int q = 0; q < n; ++q) {
        wAll_[size_t(pieces_[size_t(q)])] = w_[size_t(q)];
        if (w_[size_t(q)] > 0) { ++st_.open; st_.maxWidth = std::max(st_.maxWidth, w_[size_t(q)]); }
    }
    // fracture-cell and host arrays
    std::vector<float> W(na, 0.0f), PN(na, 0.0f), CM(na, 0.0f);
    const auto& cond = sim.hfConductance();
    for (int q = 0; q < n; ++q) {
        const int k = pieces_[size_t(q)];
        const int a = g.activeOf[size_t(g.hfCells[size_t(k)].global)];
        const int h = g.hfCells[size_t(k)].host;
        const int ah = h >= 0 && size_t(h) < g.activeOf.size() ? g.activeOf[size_t(h)] : -1;
        const float wmm = float(w_[size_t(q)] * 1e3), p = float(pnAll_[size_t(k)]);
        if (a >= 0 && size_t(a) < na) { W[size_t(a)] = wmm; PN[size_t(a)] = p; CM[size_t(a)] = float(cond.size() == size_t(nH) ? cond[size_t(k)] : 1.0); }
        if (ah >= 0 && size_t(ah) < na) {
            W[size_t(ah)] = std::max(W[size_t(ah)], wmm);
            if (w_[size_t(q)] > 0) PN[size_t(ah)] = PN[size_t(ah)] == 0.0f ? p : std::max(PN[size_t(ah)], p);
        }
    }
    s.extra.emplace_back("HFWIDTH", std::move(W));
    s.extra.emplace_back("HFPNET", std::move(PN));
    s.extra.emplace_back("HFCOND", std::move(CM));
    // matrix field (cached while the openings do not change)
    static const char* names[11] = {"HFDSXX", "HFDSYY", "HFDSZZ", "HFDSXY", "HFDSNN", "HFDSM", "HFDTAU", "HFDISPX", "HFDISPY", "HFDISPZ", "HFDISP"};
    double dw = 0;
    for (int q = 0; q < n; ++q) dw = std::max(dw, std::fabs(w_[size_t(q)] - wPrev_[size_t(q)]));
    if (!have_ || dw > 1e-9) {
        cache_.assign(12, std::vector<float>(na, 0.0f));
        const int no = int(obs_.size());
        double shadow = 0;
        const double Ev = 2 * G_ * (1 + nu_);
#pragma omp parallel for schedule(dynamic, 8) reduction(max : shadow)
        for (int i = 0; i < no; ++i) {
            const Obs& o = obs_[size_t(i)];
            double S[6], u[3];
            evalPoint(o.x, S, u);
            if (o.twin) {
                double S2[6], u2[3];
                evalPoint(o.x2, S2, u2);
                for (int j = 0; j < 6; ++j) S[j] = 0.5 * (S[j] + S2[j]);
                for (int j = 0; j < 3; ++j) u[j] = 0.5 * (u[j] + u2[j]);
            }
            const ddm::V3& nn = o.n;
            const double snn = nn[0] * (S[0] * nn[0] + S[3] * nn[1] + S[4] * nn[2]) + nn[1] * (S[3] * nn[0] + S[1] * nn[1] + S[5] * nn[2]) +
                               nn[2] * (S[4] * nn[0] + S[5] * nn[1] + S[2] * nn[2]);
            const double tr = S[0] + S[1] + S[2];
            const size_t a = size_t(o.act);
            if (a >= na) continue;
            cache_[0][a] = float(-S[0] / Pa_); cache_[1][a] = float(-S[1] / Pa_); cache_[2][a] = float(-S[2] / Pa_); cache_[3][a] = float(-S[3] / Pa_);
            cache_[4][a] = float(-snn / Pa_); cache_[5][a] = float(-tr / 3 / Pa_);
            cache_[6][a] = float(std::sqrt(0.25 * (S[0] - S[1]) * (S[0] - S[1]) + S[3] * S[3]) / Pa_);
            cache_[7][a] = float(u[0] * 1e3); cache_[8][a] = float(u[1] * 1e3); cache_[9][a] = float(u[2] * 1e3);
            cache_[10][a] = float(std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) * 1e3);
            cache_[11][a] = float((1 - 2 * nu_) / Ev * tr * 1e6);
            shadow = std::max(shadow, -snn / Pa_);
        }
        st_.maxShadow = shadow;
        st_.evaluated = no;
        wPrev_ = w_; have_ = true;
    }
    for (int j = 0; j < 11; ++j) s.extra.emplace_back(names[j], cache_[size_t(j)]);
    s.extra.emplace_back("HFSTRV", cache_[11]);
    ++st_.updates;
    st_.lastSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

}  // namespace hffield
}  // namespace ressim
