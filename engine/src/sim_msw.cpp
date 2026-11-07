// C32 multisegment wells in the black-oil Simulator: setup per report step, segment pressure solve inside computeHeads,
// SHUT devices, per-step records (msw.h for the model).
#include "ressim/simulator.h"
#include "ressim/msw.h"
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace ressim {

void Simulator::mswSetup(int step) {
    bool any = false;
    for (const auto& e : case_.extraKeywords) if (e.name == "WELSEGS") { any = true; break; }
    if (!any) { if (!msw_.empty()) { msw_.clear(); mswOf_.clear(); mswSegOf_.clear(); mswLast_.clear(); } return; }
    // signature of the definitions in force (setup only when it changes)
    std::string sig;
    for (const auto& e : case_.extraKeywords)
        if (e.section == DeckSection::Schedule && e.step <= step &&
            (e.name == "WELSEGS" || e.name == "COMPSEGS" || e.name == "WSEGVALV" || e.name == "WSEGSICD" || e.name == "WSEGAICD" ||
             e.name == "WSEGLABY" || e.name == "WSEGDFPA"))
            sig += e.name + std::to_string(e.step) + e.body;
    if (sig == mswSig_ && mswOf_.size() == wells_.size()) return;
    mswSig_ = sig;
    std::string err;
    std::vector<msw::Well> W = msw::parse(case_.extraKeywords, step, case_.units == UnitSystem::Field, err);
    if (!err.empty()) log_("Warning: multisegment wells: " + err);
    // keep the relaxation state of wells whose segment count is unchanged
    msw_ = W;
    mswOf_.assign(wells_.size(), -1);
    mswSegOf_.assign(wells_.size(), {});
    mswLast_.assign(msw_.size(), msw::Solution{});
    const int nx = grid_.nx, ny = grid_.ny;
    for (size_t m = 0; m < msw_.size(); ++m) {
        const msw::Well& M = msw_[m];
        int wi = -1;
        for (size_t i = 0; i < wells_.size(); ++i) if (case_.wells[size_t(wells_[i].input)].name == M.name) wi = int(i);
        if (wi < 0) { log_("Warning: WELSEGS for unknown well " + M.name); continue; }
        WellRT& w = wells_[size_t(wi)];
        if (case_.wells[size_t(w.input)].lgr >= 0) { log_("Warning: WELSEGS " + M.name + ": wells in local grids stay standard wells"); continue; }
        mswOf_[size_t(wi)] = int(m);
        // C36: the BHP reference depth of a multisegment well is the top-segment node depth (WELSEGS item 2), overriding
        // WELSPECS item 5 - OPM Well::updateSegments -> updateRefDepth(segments.depthTopSegment()), ECLIPSE likewise
        if (std::fabs(w.refDepth - M.topDepth) > 1e-9) {
            log_("Multisegment well " + M.name + ": BHP reference depth set to the top segment depth " + std::to_string(M.topDepth) + " (WELSEGS item 2)");
            w.refDepth = M.topDepth;
        }
        auto& so = mswSegOf_[size_t(wi)];
        so.assign(w.cell.size(), -1);
        int unlisted = 0;
        for (size_t k = 0; k < w.cell.size(); ++k) {
            const int gidx = grid_.globalOf[size_t(w.cell[k])];
            const int I = gidx % nx, J = (gidx / nx) % ny, K = gidx / (nx * ny);
            for (const auto& c : M.conns) if (c.i == I && c.j == J && c.k == K) { so[k] = M.indexOf(c.seg); break; }
            if (so[k] < 0) {   // not in COMPSEGS: nearest node by depth
                double best = 1e300;
                for (size_t s = 1; s < M.segs.size(); ++s) if (std::fabs(M.segs[s].depth - w.cdepth[k]) < best) { best = std::fabs(M.segs[s].depth - w.cdepth[k]); so[k] = int(s); }
                if (so[k] < 0) so[k] = 0;
                ++unlisted;
            }
        }
        if (unlisted) log_("Warning: WELSEGS " + M.name + ": " + std::to_string(unlisted) + " connection(s) not in COMPSEGS - attached to the segment of nearest depth");
        if (!M.warn.empty()) log_("Warning: WELSEGS " + M.name + ": " + M.warn);
        int ndev = 0, nshut = 0;
        for (const auto& g : M.segs) { ndev += g.device > 0; nshut += g.shut; }
        char b[200];
        std::snprintf(b, sizeof b, "Report step %d: multisegment well %s - %d segments, %d connections, %d devices (%d shut)",
                      step + 1, M.name.c_str(), int(M.segs.size()), int(M.conns.size()), ndev, nshut);
        log_(b);
        {
            const size_t ns = M.segs.size(), nc = w.cell.size();
            if (w.dpF.size() != ns + nc && w.dpF.size() != 2 * ns + nc && w.dpF.size() != 4 * ns + nc) w.dpF.clear();
        }   // state: friction per segment + device factor per connection
    }
}

// connections whose segment (or an ancestor) has a SHUT device are closed (called after updateConnections)
static bool shutPath(const msw::Well& M, int s) {
    for (int guard = 0; s >= 0 && guard < 100000; ++guard) {
        if (M.segs[size_t(s)].shut) return true;
        if (M.segs[size_t(s)].id == 1) break;
        s = M.indexOf(M.segs[size_t(s)].outlet);
    }
    return false;
}

void Simulator::mswHeads(WellRT& w, const std::vector<double>& qres, const std::vector<double>& rho, const std::vector<double>& mu,
                         const std::vector<double>& pCell, const std::vector<double>& aRes,
                         const std::vector<double>* gasFrac, const std::vector<double>* rhoGas, const std::vector<double>* rhoLiq,
                         const std::vector<double>* qSurf) {
    const int wi = int(&w - wells_.data());
    const int m = mswOfWell(wi);
    if (m < 0) return;
    const msw::Well& M = msw_[size_t(m)];
    const auto& so = mswSegOf_[size_t(wi)];
    const size_t n = w.cell.size(), ns = M.segs.size();
    const bool prod = w.type == WellType::Producer;
    msw::Inflow in; in.resize(ns);
    std::vector<double> rw(ns, 0), mw(ns, 0), cw(ns, 0), gw(ns, 0), lw(ns, 0);
    const bool ph = gasFrac && rhoGas && rhoLiq && gasFrac->size() == n && rhoGas->size() == n && rhoLiq->size() == n;
    const bool haveQs = qSurf && qSurf->size() == 3 * n;   // C37 tubing flash
    if (haveQs) in.qs.assign(3 * ns, 0.0);
    for (size_t k = 0; k < n; ++k) {
        const int s = k < so.size() ? so[k] : -1;
        if (s < 0) continue;
        if (w.wi[k] <= 0) { if (in.rho[size_t(s)] <= 0) { in.rho[size_t(s)] = rho[k]; in.mu[size_t(s)] = mu[k]; } continue; }
        // qres: connection reservoir rate magnitude; + toward the top (production), injection flows down (negative)
        const double q = prod ? std::fabs(qres[k]) : -std::fabs(qres[k]);
        const double aq = std::fabs(q);
        in.q[size_t(s)] += q;
        if (haveQs) for (int c = 0; c < 3; ++c) in.qs[3 * size_t(s) + size_t(c)] += (prod ? 1.0 : -1.0) * std::fabs((*qSurf)[3 * k + size_t(c)]);
        rw[size_t(s)] += aq * rho[k]; mw[size_t(s)] += aq * mu[k]; cw[size_t(s)] += aq;
        if (ph && prod) {   // C36: free gas / liquid split of the inflow (producers; injection stays homogeneous)
            const double qgk = (*gasFrac)[k] * q;
            in.qg[size_t(s)] += qgk;
            in.rhoG[size_t(s)] += std::fabs(qgk) * (*rhoGas)[k]; gw[size_t(s)] += std::fabs(qgk);
            in.rhoL[size_t(s)] += std::fabs(q - qgk) * (*rhoLiq)[k]; lw[size_t(s)] += std::fabs(q - qgk);
        }
        if (in.rho[size_t(s)] <= 0) { in.rho[size_t(s)] = rho[k]; in.mu[size_t(s)] = mu[k]; }
        // reservoir pressure referred to the segment node depth with the connection-cell density
        const double a = k < aRes.size() ? aRes[k] : 0.0;
        in.a[size_t(s)] += a;
        in.pres[size_t(s)] += a * (pCell[k] - rho[k] * U_.gravity * (w.cdepth[k] - M.segs[size_t(s)].depth));
    }
    for (size_t s = 0; s < ns; ++s) {
        if (cw[s] > 0) { in.rho[s] = rw[s] / cw[s]; in.mu[s] = mw[s] / cw[s]; }
        if (in.a[s] > 0) in.pres[s] /= in.a[s];
        in.rhoG[s] = gw[s] > 0 ? in.rhoG[s] / gw[s] : 0.0;
        in.rhoL[s] = lw[s] > 0 ? in.rhoL[s] / lw[s] : 0.0;
    }
    const bool field = case_.units == UnitSystem::Field;
    msw::Units u;
    u.g = U_.gravity;
    u.toM = field ? 0.3048 : 1.0;
    u.toRho = field ? 16.01846337 : 1.0;
    u.toQ = (field ? 0.158987294928 : 1.0) / 86400.0;
    u.toPa = 101325.0 / U_.atm;
    // C37 tubing flash: black-oil PVT of the well's first connection cell at the segment pressure
    if (haveQs && O_ + W_ + G_ > 0) {
        const FluidModel* Fp = &fl(w.cell[0]);
        const bool liveO = live_, vapG = vap_, hO = O_, hW = W_, hG = G_;
        in.flash = [Fp, liveO, vapG, hO, hW, hG](double p, const double* qq) {
            const FluidModel& F = *Fp;
            msw::Inflow::Flash r;
            const DensitySC& r0 = F.rhoSC();
            const double gf = F.gasFactor();
            const double qo = hO ? qq[0] : 0.0, qw = hW ? qq[1] : 0.0, qg = hG ? qq[2] : 0.0;
            double rs = 0, rv = 0;
            if (hO && hG) {
                if (liveO && qo > 1e-12) rs = std::min(qg / qo, F.rsSat(p));
                if (vapG && qg > 1e-12) rv = std::min(qo / qg, F.rvSat(p));
                if (rs * rv >= 1.0 - 1e-9) rv = 0;
            }
            const double d = 1.0 - rs * rv;
            const double fo = hO ? std::max(0.0, (qo - rv * qg) / d) : 0.0;   // stock-tank oil in the oil phase
            const double fg = hG ? std::max(0.0, (qg - rs * qo) / d) : 0.0;   // surface gas in the gas phase
            const AD3 P(p);
            if (fo > 0) {
                const double bo = std::max(F.boS(p, rs), 1e-30);
                r.qo = fo / bo; r.rhoO = (r0.oil + rs * r0.gas * gf) * bo;
                r.muO = F.bo(P, AD3(rs)).v / std::max(F.boMu(P, AD3(rs)).v, 1e-30);
            }
            if (qw > 0) {
                const double bw = std::max(F.bwS(p), 1e-30);
                r.qw = qw / bw; r.rhoW = r0.water * bw;
                r.muW = F.bw(P).v / std::max(F.bwMu(P).v, 1e-30);
            }
            if (fg > 0) {
                const double bg = std::max(vapG ? F.bg(P, AD3(rv)).v : F.bgS(p), 1e-30);
                r.qg = fg / bg; r.rhoG = (r0.gas * gf + rv * r0.oil) * bg;
                r.muG = vapG ? bg / std::max(F.bgMu(P, AD3(rv)).v, 1e-30) : bg / std::max(F.bgMu(P).v, 1e-30);
            }
            return r;
        };
        in.flashAtSegment = case_.mswOpt.flash;
    }
    // top node pressure: BHP at the well's reference depth, head to the top-segment depth with its density
    const double rTop = !in.rho.empty() && in.rho[0] > 0 ? in.rho[0] : w.headRho;
    const double pTop = w.bhp + rTop * U_.gravity * (M.topDepth - w.refDepth);
    const bool haveState = w.dpF.size() >= ns + n;
    // C37 committed segment storage (head density, gas holdup) behind the device factors: [nonHydro ns][factor n][rho ns][holdup ns]
    // C37 well state layout: [nonHydro ns][connection factor n][WSEGFLIM potential ns] + committed storage [rho ns][holdup ns]
    std::vector<double> stoC, flimP;
    if (w.dpF.size() == 2 * ns + n || w.dpF.size() == 4 * ns + n) flimP.assign(w.dpF.begin() + long(ns + n), w.dpF.begin() + long(2 * ns + n));
    if (w.dpF.size() == 4 * ns + n) stoC.assign(w.dpF.begin() + long(2 * ns + n), w.dpF.end());
    msw::SolveOpts so2;
    so2.dt = headDt_;
    if (case_.mswOpt.storage && !stoC.empty()) so2.sto = &stoC;
    if (!flimP.empty()) so2.flimPrev = &flimP;
    // C45: tubing fluid injected at cross-flowing connections travels down from the top through every segment on the path
    std::vector<double> xfS;
    if (w.xfRes.size() == n) {
        xfS.assign(ns, 0.0);
        for (size_t k = 0; k < n; ++k) {
            const int s0 = k < so.size() ? so[k] : -1;
            if (s0 < 0 || !(w.xfRes[k] > 0)) continue;
            for (int s = s0, guard = 0; s > 0 && guard < int(ns); ++guard) {
                xfS[size_t(s)] += w.xfRes[k];
                s = M.indexOf(M.segs[size_t(s)].outlet);
            }
        }
        so2.xfSweep = &xfS;
    }
    msw::Solution S = msw::solve(M, pTop, in, u, haveState ? std::vector<double>(w.dpF.begin(), w.dpF.begin() + long(ns)) : std::vector<double>(ns, 0.0), so2);   // first solve: relaxed from zero (potential rates overstate the friction)
    // restart-safe state: relaxed friction per segment, then the device factor per connection (applied to wiBase)
    w.dpF = S.nonHydro;
    w.dpF.resize(ns + n, 1.0);
    w.dpF.insert(w.dpF.end(), S.flimQ0.begin(), S.flimQ0.end());
    if (!stoC.empty()) w.dpF.insert(w.dpF.end(), stoC.begin(), stoC.end());
    for (size_t k = 0; k < n; ++k) {
        const int s = k < so.size() ? so[k] : -1;
        double f = s >= 0 ? S.devFac[size_t(s)] : 1.0;
        // (no relaxation of the factor: it acts as an implicit conductance in the next solve)
        w.dpF[ns + k] = f;
        if (k < w.wi.size() && w.wi[k] > 0 && k < w.wiBase.size()) w.wi[k] = w.wiBase[k] * f;
    }
    w.dpConn.assign(n, 0.0);
    w.dpFric.assign(n, 0.0);
    for (size_t k = 0; k < n; ++k) {
        const int s = k < so.size() ? so[k] : -1;
        if (s < 0) continue;
        const double rc = S.rhoHead[size_t(s)] > 0 ? S.rhoHead[size_t(s)] : S.rhoSeg[size_t(s)] > 0 ? S.rhoSeg[size_t(s)] : rho[k];
        w.dpConn[k] = S.p[size_t(s)] + rc * U_.gravity * (w.cdepth[k] - M.segs[size_t(s)].depth) - w.bhp;
        w.dpFric[k] = S.nonHydro[size_t(s)];
    }
    w.headRho = rTop;
    if (std::getenv("RESSIM_DEBUG_MSWHEADS")) {
        for (size_t s = 0; s < ns; ++s) std::fprintf(stderr, "   seg%zu p %.2f rhoH %.1f rhoS %.1f fr %.3g acc %.3g dev %.3g q %.4g in.q %.4g\n", s, S.p[s], S.rhoHead[s], S.rhoSeg[s], S.dFric[s], S.dAcc[s], S.dDev[s], S.qSeg[s], in.q[s]);
        for (size_t k = 0; k < n; ++k) std::fprintf(stderr, "   conn%zu qres %.4g pcell %.2f\n", k, qres[k], pCell[k]);
    }
    mswLast_[size_t(m)] = std::move(S);
}

// C45 producer cross-flow in a multisegment well: a connection that takes fluid from the tubing receives the fluid of its
// segment - the gross inflow of the segment's subtree (fluid flows toward the outlet); a subtree without inflow takes the
// mixture of the first segment toward the top that has one, and the well mixture (xmix) as the last resort
void Simulator::mswMixture(WellRT& w, const std::vector<double>& gIn) {
    const int wi = int(&w - wells_.data());
    const int m = mswOfWell(wi);
    const size_t n = w.cell.size();
    if (m < 0 || gIn.size() != 3 * n) return;
    const msw::Well& M = msw_[size_t(m)];
    const auto& so = mswSegOf_[size_t(wi)];
    const size_t ns = M.segs.size();
    std::vector<double> sub(3 * ns, 0.0);
    for (size_t k = 0; k < n; ++k) {
        const int s = k < so.size() ? so[k] : -1;
        if (s < 0) continue;
        for (int c = 0; c < 3; ++c) sub[3 * size_t(s) + size_t(c)] += std::max(0.0, gIn[3 * k + size_t(c)]);
    }
    const std::vector<std::vector<int>> ch = msw::children(M);
    const std::vector<int> ord = msw::topoOrder(M);
    std::vector<int> parent(ns, -1);
    for (size_t s = 0; s < ns && s < ch.size(); ++s) for (int c : ch[s]) if (c >= 0 && size_t(c) < ns) parent[size_t(c)] = int(s);
    for (size_t r = ord.size(); r-- > 0;) {                // children before parents
        const int s = ord[r], p = s >= 0 && size_t(s) < ns ? parent[size_t(s)] : -1;
        if (p >= 0) for (int c = 0; c < 3; ++c) sub[3 * size_t(p) + size_t(c)] += sub[3 * size_t(s) + size_t(c)];
    }
    w.xmixConn.assign(3 * n, 0.0);
    for (size_t k = 0; k < n; ++k) {
        int s = k < so.size() ? so[k] : -1;
        int guard = 0;
        while (s >= 0 && guard++ < int(ns) + 1) {
            const double t = sub[3 * size_t(s)] + sub[3 * size_t(s) + 1] + sub[3 * size_t(s) + 2];
            if (t > 0) { for (int c = 0; c < 3; ++c) w.xmixConn[3 * k + size_t(c)] = sub[3 * size_t(s) + size_t(c)] / t; break; }
            s = parent[size_t(s)];
        }
        if (w.xmixConn[3 * k] + w.xmixConn[3 * k + 1] + w.xmixConn[3 * k + 2] <= 0)
            for (int c = 0; c < 3; ++c) w.xmixConn[3 * k + size_t(c)] = w.xmix[c];
    }
}

void Simulator::mswRecord() {
    for (size_t m = 0; m < msw_.size() && m < mswLast_.size(); ++m) {
        const msw::Solution& S = mswLast_[m];
        if (S.p.empty()) continue;
        MswRecord r; r.t = time_; r.msw = int(m);
        r.p = S.p; r.q = S.qSeg; r.rho = S.rhoSeg; r.dFric = S.dFric; r.dDev = S.dDev; r.dAcc = S.dAcc; r.holdup = S.holdup; r.noSlip = S.noSlip;
        r.oilHold = S.oilHold; r.dPump = S.dPump; r.flimFac = S.flimFac; r.rhoHead = S.rhoHead;
        mswRec_.push_back(std::move(r));
        // C37: commit the segment content of the accepted step (storage state for the next step)
        for (size_t wi = 0; wi < mswOf_.size() && wi < wells_.size(); ++wi) {
            if (mswOf_[wi] != int(m)) continue;
            WellRT& w = wells_[wi];
            const size_t ns = S.p.size(), nc = w.cell.size();
            if (w.dpF.size() < 2 * ns + nc || S.rhoHead.size() != ns) continue;
            const bool first = w.dpF.size() == 2 * ns + nc;   // no committed content yet (well opened in this step)
            w.dpF.resize(2 * ns + nc);
            if (first && mswXflow() && S.rhoHead[0] > 0) {
                // C45: the first content of the tubing is the produced / injected well mixture of the top segment, not the
                // steady value of the lagged start-of-step rates: those are potential rates (every open connection flowing)
                // and left stagnant segments - below a cross-flowing or closed connection - holding a liquid column for the
                // rest of the run (MSW-2D-VERT-02: +25 bar at the lower connection). Flowing segments are swept to their own
                // content within a step, as before.
                w.dpF.insert(w.dpF.end(), ns, S.rhoHead[0]);
                w.dpF.insert(w.dpF.end(), ns, S.holdup[0]);
            } else {
                w.dpF.insert(w.dpF.end(), S.rhoHead.begin(), S.rhoHead.end());
                w.dpF.insert(w.dpF.end(), S.holdup.begin(), S.holdup.end());
            }
        }
    }
}

// exposed for updateConnections
bool mswConnectionShut(const msw::Well& M, int seg) { return seg >= 0 && shutPath(M, seg); }

}  // namespace ressim
