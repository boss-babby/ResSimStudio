// C26: DDM apertures of the in-run hydraulic fractures - see hf_ddm.h
#include "ressim/hf_ddm.h"
#include "ressim/case.h"
#include "ressim/ddm.h"
#include "ressim/simulator.h"
#include <algorithm>
#include <chrono>
#include <cmath>

namespace ressim {
namespace hfddm {

bool Aperture::setup(const Simulator& sim, std::string& err) {
    active_ = false; out_ = Output{}; pieces_.clear(); Ko_.clear(); wIn_.clear(); zMid_.clear(); cubic_.clear();
    const Case& c = sim.caseData();
    const Grid& g = sim.grid();
    const auto& H = c.grid.hfracs;
    bool any = false;
    for (const auto& F : H) if (F.ddmWidth) any = true;
    if (!any) { err = "no HFGROW ... DDMWIDTH fracture"; return false; }
    const int nH = g.numHfCells();
    if (nH == 0 || c.grid.hfPieces.size() != size_t(nH)) { err = "the fracture pieces are not built (HFRAC outside the grid?)"; return false; }
    const bool fld = c.units == UnitSystem::Field;
    Mlen_ = fld ? 0.3048 : 1.0; Pa_ = fld ? 6894.757 : 1e5;
    // elastic moduli: DDMWIDTH E nu, else the GEOMECH moduli, else 20 GPa / 0.25 (first DDMWIDTH fracture decides)
    double E = 20e9, nu = 0.25;
    for (const auto& F : H) if (F.ddmWidth) {
        if (F.ddmE > 0) { E = F.ddmE * 1e9; nu = F.ddmNu; }
        else if (c.geoSeq.enabled && !c.geoSeq.E.empty() && c.geoSeq.E[0] > 0) { E = c.geoSeq.E[0]; nu = c.geoSeq.nu.empty() ? 0.25 : c.geoSeq.nu[0]; }
        break;
    }
    const double G = E / (2 * (1 + nu));
    std::vector<ddm::Element> el;
    halfSpace_ = true;
    for (int k = 0; k < nH; ++k) {
        const auto& P = c.grid.hfPieces[size_t(k)];
        const auto& F = H[size_t(P.frac)];
        if (!F.ddmWidth) continue;
        if (P.zt < 0) halfSpace_ = false;
        const ddm::V3 a = {P.xa * Mlen_, P.ya * Mlen_, -P.zt * Mlen_}, b = {P.xb * Mlen_, P.yb * Mlen_, -P.zt * Mlen_};
        const ddm::V3 cc = {P.xb * Mlen_, P.yb * Mlen_, -P.zb * Mlen_}, d = {P.xa * Mlen_, P.ya * Mlen_, -P.zb * Mlen_};
        ddm::Element e = ddm::quadElement(a, b, cc, d);
        if (e.tris.empty() || !(e.area > 0)) continue;
        pieces_.push_back(k);
        el.push_back(e);
        wIn_.push_back(std::max(1e-6, F.width * Mlen_));
        zMid_.push_back(0.5 * (P.zt + P.zb));
        cubic_.push_back(F.perm <= 0 ? 1 : 0);
    }
    if (pieces_.empty()) { err = "DDMWIDTH fractures have no pieces"; return false; }
    if (pieces_.size() > 3000) { err = "more than 3000 fracture pieces - DDM apertures skipped"; return false; }
    const auto t0 = std::chrono::steady_clock::now();
    Ko_ = ddm::buildOpening(el, G, nu, halfSpace_);
    out_.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    out_.elements = int(pieces_.size());
    out_.width.assign(size_t(nH), 0.0); out_.mult.assign(size_t(nH), 1.0); out_.pnet.assign(size_t(nH), 0.0);
    active_ = true;
    return true;
}

void Aperture::update(Simulator& sim) {
    if (!active_) return;
    const Case& c = sim.caseData();
    const Grid& g = sim.grid();
    const int n = int(pieces_.size()), nH = g.numHfCells();
    const auto& sMin = sim.minStress();
    std::vector<double> pn(size_t(n), 0.0);
    std::vector<int> idx;   // open pieces (the growth state decides which pieces exist as fracture)
    for (int q = 0; q < n; ++q) {
        const int k = pieces_[size_t(q)];
        const auto& P = c.grid.hfPieces[size_t(k)];
        const auto& F = c.grid.hfracs[size_t(P.frac)];
        double sig = 0;
        bool feOk = false;
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
        out_.pnet[size_t(k)] = pnet;
        if (sim.hfCellOpenFraction(k) > 0) { idx.push_back(q); pn[size_t(q)] = pnet * Pa_; }
    }
    const int m = int(idx.size());
    std::vector<double> w(size_t(n), 0.0);
    if (m > 0) {
        std::vector<double> K(static_cast<size_t>(m) * static_cast<size_t>(m)), p(static_cast<size_t>(m)), ws;
        for (int r = 0; r < m; ++r) {
            p[size_t(r)] = pn[size_t(idx[size_t(r)])];
            for (int s = 0; s < m; ++s) K[size_t(r) * size_t(m) + size_t(s)] = Ko_[size_t(idx[size_t(r)]) * size_t(n) + size_t(idx[size_t(s)])];
        }
        ddm::solveOpening(K, m, p, ws, 1.0);
        for (int r = 0; r < m; ++r) w[size_t(idx[size_t(r)])] = ws[size_t(r)];
    }
    std::vector<double> mult(size_t(nH), 1.0);
    const std::vector<double>& prev = sim.hfConductance();
    if (prev.size() == size_t(nH)) mult = prev;   // keep multipliers of other fractures (G1-c)
    out_.maxWidth = 0; out_.volume = 0; out_.open = 0;
    for (int q = 0; q < n; ++q) {
        const int k = pieces_[size_t(q)];
        const auto& P = c.grid.hfPieces[size_t(k)];
        const double area = std::hypot(P.xb - P.xa, P.yb - P.ya) * std::fabs(P.zb - P.zt) * Mlen_ * Mlen_;
        out_.width[size_t(k)] = w[size_t(q)];
        if (w[size_t(q)] > 0) { ++out_.open; out_.volume += w[size_t(q)] * area; out_.maxWidth = std::max(out_.maxWidth, w[size_t(q)]); }
        const double we = std::max(w[size_t(q)], 0.1 * wIn_[size_t(q)]);
        const double r = we / wIn_[size_t(q)];
        mult[size_t(k)] = sim.hfCellOpenFraction(k) > 0 ? (cubic_[size_t(q)] ? r * r * r : r) : 1.0;
        out_.mult[size_t(k)] = mult[size_t(k)];
    }
    sim.setHfConductance(mult);
    ++out_.updates;
    Output::Step st; st.t = sim.time(); st.maxWidth = out_.maxWidth; st.volume = out_.volume; st.open = out_.open;
    st.minMult = 1e300; st.maxMult = -1e300; st.maxPnet = -1e300;
    for (int q = 0; q < n; ++q) {
        const int k = pieces_[size_t(q)];
        st.minMult = std::min(st.minMult, out_.mult[size_t(k)]); st.maxMult = std::max(st.maxMult, out_.mult[size_t(k)]);
        st.maxPnet = std::max(st.maxPnet, out_.pnet[size_t(k)]);
    }
    if (n == 0) { st.minMult = st.maxMult = 1; st.maxPnet = 0; }
    st.width.assign(out_.width.begin(), out_.width.end());
    out_.steps.push_back(std::move(st));
}

}  // namespace hfddm
}  // namespace ressim
