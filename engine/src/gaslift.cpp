// ResSim - C35 gas-lift optimisation (ECLIPSE / OPM Flow LIFTOPT, WLIFTOPT, GLIFTOPT).
//
// Input (kept verbatim, exact deck round trip): LIFTOPT (increment, minimum economic gradient, minimum interval between
// optimisations, optimise in all Newton iterations), WLIFTOPT (well, YES / NO, maximum lift gas rate (default: largest
// ALQ of the well's VFPPROD table), weighting factor, minimum lift gas rate (< 0: enough to flow), incremental gas
// weighting factor), GLIFTOPT (group, maximum lift gas supply, maximum total gas = produced + lift gas).
//
// Method (ENGINES_AND_MATH 6.30):
//   * lift performance per well at the start of a time step (frozen reservoir state, as the OPM GasLiftSingleWell
//     stage): for ALQ_j = j * increment (j = 0 .. up to the maximum) the operating point of the inflow and the VFP lift
//     curve at the THP limit (Simulator::thpOperatingPoint), capped by the BHP limit and the well's own rate limits
//     -> q_o(ALQ_j), q_w, q_g;
//   * allocation in whole increments to the well with the largest weighted incremental gradient
//         G = w_f dq_o / (dALQ + w_g dq_g)                                         (OPM calcEcoGradient)
//     while G > the minimum economic gradient, subject to the well maximum, the group lift gas supply and total gas
//     limits (GLIFTOPT) of every group above the well, and the group production limits / targets (GCONPROD): the part
//     of an increment's oil that a group limit would cut back does not count (OPM updateRatesToGroupLimits);
//   * the gradient of a well looks ahead over several increments (concave hull of the lift curve), so a well that only
//     starts to flow after a few increments of lift gas (dead well) is allocated as one block - this replaces the
//     separate 'enough to flow' rule of OPM; a negative WLIFTOPT minimum still forces that block;
//   * the result sets the ALQ (lift gas rate) used by the VFP lookups of the time step; the group targets are then
//     allocated on the potentials with the new lift gas.
#include "ressim/simulator.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

namespace ressim {
namespace {
std::string glUp(std::string s) { for (char& c : s) c = char(std::toupper(static_cast<unsigned char>(c))); return s; }
bool glMatch(const std::string& pat, const std::string& name) {   // '*' / '?' anywhere
    const std::string a = glUp(pat), b = glUp(name);
    size_t i = 0, j = 0, star = std::string::npos, mark = 0;
    while (j < b.size()) {
        if (i < a.size() && (a[i] == '?' || a[i] == b[j])) { ++i; ++j; }
        else if (i < a.size() && a[i] == '*') { star = i++; mark = j; }
        else if (star != std::string::npos) { i = star + 1; j = ++mark; }
        else return false;
    }
    while (i < a.size() && a[i] == '*') ++i;
    return i == a.size();
}
bool given(const std::vector<std::string>& r, size_t i) { return i < r.size() && !r[i].empty() && r[i].find('*') == std::string::npos; }
char glBuf[512];
template <class... A> std::string gfmt(const char* f, A... a) { std::snprintf(glBuf, sizeof glBuf, f, a...); return glBuf; }
}  // namespace

void Simulator::gasLiftControls(int step) {
    if (step == 0 || glw_.size() != wells_.size()) {
        if (step == 0) { gloInc_ = gloEco_ = gloInterval_ = gloCum_ = 0; gloLast_ = -1e30; gloRuns_ = 0; glg_.clear(); for (auto& w : wells_) { w.alq = 0; w.alqLift = false; } }
        glw_.assign(wells_.size(), GasLiftWell());
    }
    if (step >= 0)
        for (const ExtraKeyword& k : case_.extraKeywords) {
            if (k.section != DeckSection::Schedule || k.step != step) continue;
            if (k.name == "LIFTOPT") {
                for (const auto& r : verbatimRecords(k.body)) {
                    if (r.empty()) continue;
                    gloInc_ = verbatimNum(r, 0, 0);
                    gloEco_ = verbatimNum(r, 1, 0);
                    gloInterval_ = verbatimNum(r, 2, 0);
                    log_(gloInc_ > 0 ? gfmt("  gas-lift optimisation: increment %g, minimum economic gradient %g, interval %g d", gloInc_, gloEco_, gloInterval_)
                                     : std::string("  gas-lift optimisation switched off (LIFTOPT increment <= 0)"));
                    break;
                }
            } else if (k.name == "WLIFTOPT") {
                for (const auto& r : verbatimRecords(k.body)) {
                    if (r.empty() || r[0].empty()) continue;
                    std::string pat = r[0];
                    if (pat.size() >= 2 && pat.front() == '\'') pat = pat.substr(1, pat.size() - 2);
                    int hits = 0;
                    for (size_t i = 0; i < wells_.size(); ++i) {
                        if (!glMatch(pat, case_.wells[size_t(wells_[i].input)].name)) continue;
                        GasLiftWell& G = glw_[i];
                        G.on = verbatimStr(r, 1, "YES") == "YES";
                        G.maxRate = given(r, 2) ? verbatimNum(r, 2, -1) : -1;
                        G.weight = given(r, 3) ? verbatimNum(r, 3, 1) : 1.0;
                        if (!(G.weight > 0)) G.weight = 1.0;
                        G.minRate = given(r, 4) ? verbatimNum(r, 4, 0) : 0.0;
                        G.gasW = given(r, 5) ? verbatimNum(r, 5, 0) : 0.0;
                        // NO: the lift gas rate stays fixed at item 3 (if given), else at its current value
                        G.fixed = G.on ? -1 : (given(r, 2) ? std::max(0.0, G.maxRate) : (wells_[i].alqLift ? wells_[i].alq : -1));
                        ++hits;
                    }
                    if (!hits) log_("Warning: WLIFTOPT - no well matches '" + pat + "'");
                }
            } else if (k.name == "GLIFTOPT") {
                for (const auto& r : verbatimRecords(k.body)) {
                    if (r.empty() || r[0].empty()) continue;
                    std::string g = r[0];
                    if (g.size() >= 2 && g.front() == '\'') g = g.substr(1, g.size() - 2);
                    GasLiftGroup& G = glg_[glUp(g)];
                    G.maxLift = given(r, 1) ? verbatimNum(r, 1, -1) : -1;
                    G.maxTotal = given(r, 2) ? verbatimNum(r, 2, -1) : -1;
                    if (G.maxLift < 0 && G.maxTotal < 0) glg_.erase(glUp(g));
                }
            }
        }
    // ALQ in force of the wells that are not optimised: WLIFTOPT NO value, else WCONPROD item 12 / WELTARG LIFT
    for (size_t i = 0; i < wells_.size(); ++i) {
        WellRT& w = wells_[i];
        const GasLiftWell& G = glw_[i];
        if (G.on && w.type == WellType::Producer) { if (!w.alqLift) { w.alq = 0; w.alqLift = true; } continue; }
        if (G.fixed >= 0) { w.alq = G.fixed; w.alqLift = true; }
        else { w.alq = std::max(0.0, w.own.alq); w.alqLift = false; }
    }
}

double Simulator::liftGasRate(const WellRT& w) const {
    if (w.type != WellType::Producer || !w.open || !(w.alq > 0)) return 0.0;
    if (w.alqLift) return w.alq;
    // a fixed WCONPROD / WELTARG ALQ is a lift gas rate only for GRAT (or blank) ALQ tables
    const VfpTable* t = wellVfp(w);
    if (!t) return 0.0;
    std::string ty = t->alq;
    ty.erase(std::remove(ty.begin(), ty.end(), '\''), ty.end());
    return (ty.empty() || ty == "GRAT") ? w.alq : 0.0;
}

bool Simulator::gasLiftOptimize() {
    if (!(gloInc_ > 0) || glw_.size() != wells_.size()) return false;
    bool anyOn = false;
    for (size_t i = 0; i < wells_.size(); ++i) if (glw_[i].on && wells_[i].type == WellType::Producer && wells_[i].open) anyOn = true;
    if (!anyOn) return false;
    if (gloLast_ > -1e29 && time_ < gloLast_ + gloInterval_ - 1e-9) return false;
    gloLast_ = time_;
    ++gloRuns_;
    const double inc = gloInc_;
    struct Cand { size_t w; std::vector<double> a, qo, qw, qg; int j = 0; };
    std::vector<Cand> cs;
    std::vector<char> isCand(wells_.size(), 0);
    for (size_t i = 0; i < wells_.size(); ++i) {
        WellRT& w = wells_[i];
        const GasLiftWell& G = glw_[i];
        if (!G.on || w.type != WellType::Producer || !w.open) continue;
        const std::string& nm = case_.wells[size_t(w.input)].name;
        const VfpTable* t = wellVfp(w);
        if (!t || !(std::max(w.own.thp, w.thpNet) > 0)) {
            if (gloNotes_.insert(nm).second) log_("NOTE: gas lift - well " + nm + " has no THP limit / VFPPROD table: lift gas has no effect, not optimised");
            w.alq = 0; continue;
        }
        double amax = G.maxRate;
        if (amax < 0) {   // default: the largest ALQ of the table (GRAT tables; other ALQ types need item 3)
            std::string ty = t->alq;
            ty.erase(std::remove(ty.begin(), ty.end(), '\''), ty.end());
            amax = (ty.empty() || ty == "GRAT") && !t->alqV.empty() ? t->alqV.back() : 0.0;
            if (!(amax > 0) && gloNotes_.insert(nm + "#max").second)
                log_("NOTE: gas lift - well " + nm + ": no maximum lift gas rate (WLIFTOPT item 3) and no GRAT ALQ axis in VFP table - not optimised");
        }
        Cand c; c.w = i;
        for (int j = 0; j < 400; ++j) {
            double a = j * inc;
            if (a > amax + 1e-9 * inc) { if (c.a.empty() || c.a.back() < amax - 1e-9 * inc) a = amax; else break; }
            c.a.push_back(a);
            if (a >= amax) break;
        }
        if (c.a.empty()) c.a.push_back(0.0);
        for (double a : c.a) {
            double b = 0, q[3] = {0, 0, 0};
            if (!thpOperatingPoint(w, a, true, b, q)) q[0] = q[1] = q[2] = 0;
            // the well's own rate limits (WCONPROD / WELTARG) cap the deliverable rates
            double f = 1.0;
            if (w.own.orat > 0 && q[PO] > 0) f = std::min(f, w.own.orat / q[PO]);
            if (w.own.wrat > 0 && q[PW] > 0) f = std::min(f, w.own.wrat / q[PW]);
            if (w.own.grat > 0 && q[PG] > 0) f = std::min(f, w.own.grat / q[PG]);
            if (w.own.lrat > 0 && q[PO] + q[PW] > 0) f = std::min(f, w.own.lrat / (q[PO] + q[PW]));
            c.qo.push_back(f * q[PO]); c.qw.push_back(f * q[PW]); c.qg.push_back(f * q[PG]);
        }
        // minimum lift gas: > 0 at least that rate (unless the well already meets a rate limit with less), < 0 enough to flow
        int j0 = 0;
        if (G.minRate > 0) { while (j0 + 1 < int(c.a.size()) && c.a[size_t(j0)] < G.minRate - 1e-9 * inc) ++j0; }
        else if (G.minRate < 0) { while (j0 + 1 < int(c.a.size()) && !(c.qo[size_t(j0)] > 0)) ++j0; if (!(c.qo[size_t(j0)] > 0)) j0 = 0; }
        for (int j = 0; j < j0; ++j)
            if (c.qo[size_t(j)] > 0 && c.qo[size_t(j)] >= c.qo[size_t(j0)] * (1 - 1e-9)) { j0 = j; break; }
        c.j = j0;
        cs.push_back(std::move(c));
        isCand[i] = 1;
    }
    if (cs.empty()) return false;
    // group sums (efficiency weighted): oil, water, gas, lift gas; base = the wells that are not optimised (current rates)
    struct GS { double o = 0, w = 0, g = 0, a = 0; };
    std::map<std::string, GS> gs;
    auto chainOf = [&](size_t i) -> const std::vector<std::string>& {
        static const std::vector<std::string> none;
        return i < wellChain_.size() ? wellChain_[i] : none;
    };
    for (size_t i = 0; i < wells_.size(); ++i) {
        const WellRT& w = wells_[i];
        if (isCand[i] || w.type != WellType::Producer || !w.open) continue;
        for (const std::string& g : chainOf(i)) {
            GS& s = gs[glUp(g)];
            s.o += w.eff * w.rates.oil; s.w += w.eff * w.rates.water; s.g += w.eff * w.rates.gas; s.a += w.eff * liftGasRate(w);
        }
    }
    for (const Cand& c : cs) {
        const WellRT& w = wells_[c.w];
        const size_t j = size_t(c.j);
        for (const std::string& g : chainOf(c.w)) {
            GS& s = gs[glUp(g)];
            s.o += w.eff * c.qo[j]; s.w += w.eff * c.qw[j]; s.g += w.eff * c.qg[j]; s.a += w.eff * c.a[j];
        }
    }
    // group production caps (GCONPROD targets and limits) by upper-case name
    std::map<std::string, const GroupProdControl*> caps;
    for (const auto& [n, g] : gprod_) caps[glUp(n)] = &g;
    for (int iter = 0; iter < 100000; ++iter) {
        double best = -1e300; int bc = -1, bk = 0;
        for (int ci = 0; ci < int(cs.size()); ++ci) {
            const Cand& c = cs[size_t(ci)];
            const WellRT& w = wells_[c.w];
            const GasLiftWell& G = glw_[c.w];
            const size_t j = size_t(c.j);
            for (size_t k = j + 1; k < c.a.size(); ++k) {
                const double dA = c.a[k] - c.a[j], dqo = c.qo[k] - c.qo[j], dqw = c.qw[k] - c.qw[j], dqg = c.qg[k] - c.qg[j];
                bool ok = true; double sc = 1.0;
                for (const std::string& g0 : chainOf(c.w)) {
                    const std::string g = glUp(g0);
                    const GS& s = gs[g];
                    auto gl = glg_.find(g);
                    if (gl != glg_.end()) {
                        if (gl->second.maxLift >= 0 && s.a + w.eff * dA > gl->second.maxLift * (1 + 1e-9) + 1e-9) ok = false;
                        if (gl->second.maxTotal >= 0 && s.g + s.a + w.eff * (dqg + dA) > gl->second.maxTotal * (1 + 1e-9) + 1e-9) ok = false;
                    }
                    auto cp = caps.find(g);
                    if (cp != caps.end()) {
                        auto cut = [&](double cap, double Q, double dQ) {
                            if (!(cap > 0) || !(dQ > 0)) return;
                            sc = std::min(sc, std::clamp((cap - Q) / (w.eff * dQ), 0.0, 1.0));
                        };
                        // OPM (GroupKeywordHandlers GCONPROD, GasLiftGroupInfo): a phase limit constrains the lift
                        // allocation only when it is the group's control mode or its exceed procedure is not NONE
                        const GroupProdControl& P = *cp->second;
                        auto active = [&](Control m) { return P.mode == m || glUp(P.actionFor(m)) != "NONE"; };
                        if (active(Control::ORAT)) cut(P.orat, s.o, dqo);
                        if (active(Control::WRAT)) cut(P.wrat, s.w, dqw);
                        if (active(Control::GRAT)) cut(P.grat, s.g, dqg);
                        if (active(Control::LRAT)) cut(P.lrat, s.o + s.w, dqo + dqw);
                    }
                }
                if (!ok) break;   // larger blocks need even more lift gas
                const double den = dA + G.gasW * dqg;
                if (!(den > 0)) continue;
                const double grad = G.weight * sc * dqo / den;
                if (sc * dqo > 0 && grad > best + 1e-12 * std::fabs(best)) { best = grad; bc = ci; bk = int(k); }
            }
        }
        if (bc < 0 || !(best > gloEco_)) break;
        Cand& c = cs[size_t(bc)];
        const WellRT& w = wells_[c.w];
        const size_t j = size_t(c.j), k = size_t(bk);
        for (const std::string& g : chainOf(c.w)) {
            GS& s = gs[glUp(g)];
            s.o += w.eff * (c.qo[k] - c.qo[j]); s.w += w.eff * (c.qw[k] - c.qw[j]); s.g += w.eff * (c.qg[k] - c.qg[j]); s.a += w.eff * (c.a[k] - c.a[j]);
        }
        c.j = bk;
    }
    bool changed = false;
    const bool dbg = std::getenv("RESSIM_DEBUG") != nullptr;
    for (const Cand& c : cs) {
        WellRT& w = wells_[c.w];
        const double a = c.a[size_t(c.j)];
        if (std::fabs(a - w.alq) > 1e-9 * std::max(1.0, inc)) changed = true;
        w.alq = a; w.alqLift = true;
        if (dbg) std::fprintf(stderr, "  [gas lift t %.2f %s alq %.6g qo %.6g (alq 0: %.6g)]\n", time_, case_.wells[size_t(w.input)].name.c_str(),
                              a, c.qo[size_t(c.j)], c.qo[0]);
    }
    return changed;
}

}  // namespace ressim
