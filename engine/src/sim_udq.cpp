// ResSim - C31 UDQ / UDA / ACTIONX in the black-oil simulator (see udq.h for the deck syntax and semantics)
#include "ressim/simulator.h"
#include "ressim/udq.h"
#include "ressim/io.h"
#include "ressim/deck_parser.h"
#include "ressim/pyaction.h"
#include <filesystem>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <algorithm>
#include <limits>

namespace ressim {

namespace {
const double NaN = std::numeric_limits<double>::quiet_NaN();
// proleptic Gregorian calendar (H. Hinnant)
long daysFromCivilU(int y, int m, int d) {
    y -= m <= 2;
    const long era = (y >= 0 ? y : y - 399) / 400;
    const unsigned yoe = unsigned(y - era * 400);
    const unsigned doy = unsigned((153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1);
    const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return era * 146097 + long(doe) - 719468;
}
void civilFromDaysU(long z, int& y, int& m, int& d) {
    z += 719468;
    const long era = (z >= 0 ? z : z - 146096) / 146097;
    const unsigned doe = unsigned(z - era * 146097);
    const unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    const unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    const unsigned mp = (5 * doy + 2) / 153;
    d = int(doy - (153 * mp + 2) / 5 + 1);
    m = int(mp < 10 ? mp + 3 : mp - 9);
    y = int(yoe) + int(era) * 400 + (m <= 2);
}
std::string g(double v) { char b[40]; std::snprintf(b, sizeof(b), "%.6g", v); return b; }
}  // namespace

int Simulator::numActions() const { return actions_ ? int(actions_->size()) : 0; }

void Simulator::udqSetup() {
    udq_.reset(); actions_ = std::make_shared<std::vector<udq::Action>>(); actEvents_.clear(); actStep_.clear(); actText_.clear();
    baseSchedule_.clear(); baseWells_.clear(); udqGroups_.clear();
    wCum_.assign(wells_.size() * 5, 0.0);
    auto set = std::make_shared<udq::Set>();
    bool anyUdq = false, anyUda = false;
    for (const ExtraKeyword& e : case_.extraKeywords) {
        if (e.name == "UDQ") {
            std::string err;
            if (!set->add(e.body, e.step, err)) log_("Warning: UDQ (report step " + std::to_string(e.step + 1) + "): " + err);
            anyUdq = true;
        } else if (e.name == "UDT") {   // C38 user-defined tables for the UDQ expressions
            std::string err;
            if (!set->addTables(e.body, e.step, err)) log_("Warning: UDT (report step " + std::to_string(e.step + 1) + "): " + err);
        } else if (e.name == "ACTIONX") {
            udq::Action a; std::string err;
            if (!udq::parseAction(e.body, e.step, a, err)) { log_("Warning: " + err + " - the action is not evaluated by the ResSim engine"); continue; }
            actions_->push_back(std::move(a));
        }
    }
    for (const auto& st : case_.schedule) {
        for (const auto& wc : st.controls) if (!wc.uda.empty()) anyUda = true;
        for (const auto& g : st.gprod) if (!g.uda.empty()) anyUda = true;   // C48 group UDA
        for (const auto& g : st.ginj) if (!g.uda.empty()) anyUda = true;
    }
    // the UDQ set is also needed for actions (conditions on UDQs, UDQ records inside an action body)
    if (anyUdq || anyUda || !actions_->empty()) udq_ = set;
    if (!actions_->empty()) { baseSchedule_ = case_.schedule; baseWells_ = case_.wells; }
    // groups: every group of a well chain and of GRUPTREE (FIELD last)
    for (const auto& w : case_.wells) {
        std::string gname = w.group;
        for (int d = 0; d < 32 && !gname.empty() && gname != "FIELD"; ++d) {
            if (std::find(udqGroups_.begin(), udqGroups_.end(), gname) == udqGroups_.end()) udqGroups_.push_back(gname);
            std::string parent = "FIELD";
            for (const auto& [ch, pa] : case_.groupTree) if (ch == gname) parent = pa;
            gname = parent;
        }
    }
    for (const auto& [ch, pa] : case_.groupTree) {
        if (ch != "FIELD" && std::find(udqGroups_.begin(), udqGroups_.end(), ch) == udqGroups_.end()) udqGroups_.push_back(ch);
        if (pa != "FIELD" && std::find(udqGroups_.begin(), udqGroups_.end(), pa) == udqGroups_.end()) udqGroups_.push_back(pa);
    }
    udqGroups_.push_back("FIELD");
    if (udq_) {
        const auto nm = udq_->names();
        std::string l;
        for (const auto& n : nm) l += (l.empty() ? "" : " ") + n;
        log_("UDQ: " + std::to_string(nm.size()) + " quantities" + (l.empty() ? std::string() : " (" + l + ")") +
             (anyUda ? ", UDA targets in well controls" : "") +
             (actions_->empty() ? std::string() : ", " + std::to_string(actions_->size()) + " ACTIONX blocks evaluated at the end of every report step"));
    }
}

// FIELD / group / well quantities of the current state for UDQ and ACTIONX expressions (s must outlive src)
void Simulator::udqSource(const SummaryPoint& s, udq::Source& src) const {
    for (const auto& w : wells_) src.wells.push_back(case_.wells[size_t(w.input)].name);
    src.groups = udqGroups_;
    src.days = time_;
    {
        int y, m, d;
        civilFromDaysU(daysFromCivilU(case_.startYear, case_.startMonth, case_.startDay) + long(std::floor(time_ + 1e-9)), y, m, d);
        src.year = y; src.month = m; src.day = d;
    }
    const SummaryPoint* sp = &s;
    auto ratio = [](double a, double b) { return b > 0 ? a / b : 0.0; };
    auto wellQ = [this, ratio](const std::string& kw, size_t i) -> double {
        const WellRT& w = wells_[i];
        const WellRates& r = w.rates;
        const double* cu = i * 5 + 4 < wCum_.size() ? &wCum_[i * 5] : nullptr;
        if (kw == "WOPR") return r.oil;
        if (kw == "WWPR") return r.water;
        if (kw == "WGPR") return r.gas;
        if (kw == "WLPR") return r.oil + r.water;
        if (kw == "WWIR") return r.waterInj;
        if (kw == "WGIR") return r.gasInj;
        if (kw == "WVPR" || kw == "WVIR") return r.resv;
        if (kw == "WBHP") return w.bhp;
        if (kw == "WTHP") return w.thpNet > 0 ? w.thpNet : (w.ctrl.thp > 0 ? w.ctrl.thp : NaN);
        if (kw == "WWCT") return ratio(r.water, r.oil + r.water);
        if (kw == "WGOR") return ratio(r.gas, r.oil);
        if (kw == "WGLR") return ratio(r.gas, r.oil + r.water);
        if (kw == "WWGR") return ratio(r.water, r.gas);
        if (kw == "WOGR") return ratio(r.oil, r.gas);
        if (kw == "WSTAT") return w.open ? 1.0 : 0.0;
        if (kw.size() == 5 && kw.compare(0, 4, "WPWE") == 0 && kw[4] >= '0' && kw[4] <= '7') return double(w.pweLast[kw[4] - '0']);   // C43
        if (kw == "WOPRH") return r.oilH >= 0 ? r.oilH : NaN;
        if (kw == "WWPRH") return r.waterH >= 0 ? r.waterH : NaN;
        if (kw == "WGPRH") return r.gasH >= 0 ? r.gasH : NaN;
        if (cu) {
            if (kw == "WOPT") return cu[0];
            if (kw == "WWPT") return cu[1];
            if (kw == "WGPT") return cu[2];
            if (kw == "WLPT") return cu[0] + cu[1];
            if (kw == "WWIT") return cu[3];
            if (kw == "WGIT") return cu[4];
        }
        return NaN;
    };
    src.well = wellQ;
    src.field = [this, sp, ratio](const std::string& kw) -> double {
        const SummaryPoint& s = *sp;
        if (kw == "FVPR" || kw == "FVIR") {   // reservoir voidage rates of the producers / injectors
            double t = 0;
            for (size_t i = 0; i < wells_.size() && i < s.wells.size(); ++i)
                if ((wells_[i].type == WellType::Producer) == (kw == "FVPR")) t += std::fabs(s.wells[i].resv);
            return t;
        }
        if (kw == "FOPR") return s.fopr;
        if (kw == "FWPR") return s.fwpr;
        if (kw == "FGPR") return s.fgpr;
        if (kw == "FLPR") return s.fopr + s.fwpr;
        if (kw == "FWIR") return s.fwir;
        if (kw == "FGIR") return s.fgir;
        if (kw == "FOPT") return s.fopt;
        if (kw == "FWPT") return s.fwpt;
        if (kw == "FGPT") return s.fgpt;
        if (kw == "FLPT") return s.fopt + s.fwpt;
        if (kw == "FWIT") return s.fwit;
        if (kw == "FGIT") return s.fgit;
        if (kw == "FPR") return s.fpr;
        if (kw == "FPRP") return s.fprp;
        if (kw == "FWCT") return s.fwct;
        if (kw == "FGOR") return s.fgor;
        if (kw == "FGLR") return ratio(s.fgpr, s.fopr + s.fwpr);
        if (kw == "FOIP") return s.oip;
        if (kw == "FWIP") return s.wip;
        if (kw == "FGIP") return s.gip;
        if (kw == "FAQR") return s.faqrS;
        if (kw == "FAQT") return s.faqtS;
        return NaN;
    };
    // group = sum over the wells whose group chain contains it (FIELD: all wells)
    auto chainHas = [this](size_t i, const std::string& grp) {
        if (grp == "FIELD") return true;
        std::string gname = case_.wells[size_t(wells_[i].input)].group;
        for (int d = 0; d < 32 && !gname.empty() && gname != "FIELD"; ++d) {
            if (gname == grp) return true;
            std::string parent = "FIELD";
            for (const auto& [ch, pa] : case_.groupTree) if (ch == gname) parent = pa;
            gname = parent;
        }
        return false;
    };
    src.group = [this, wellQ, chainHas, ratio](const std::string& kw, size_t gi) -> double {
        if (gi >= udqGroups_.size() || kw.size() < 2) return NaN;
        const std::string& grp = udqGroups_[gi];
        auto sum = [&](const std::string& wkw) {
            double t = 0;
            for (size_t i = 0; i < wells_.size(); ++i) if (chainHas(i, grp)) { const double v = wellQ(wkw, i); if (std::isfinite(v)) t += v; }
            return t;
        };
        if (kw == "GWCT") return ratio(sum("WWPR"), sum("WLPR"));
        if (kw == "GGOR") return ratio(sum("WGPR"), sum("WOPR"));
        if (kw == "GGLR") return ratio(sum("WGPR"), sum("WLPR"));
        static const char* sums[] = {"OPR", "WPR", "GPR", "LPR", "WIR", "GIR", "VPR", "VIR", "OPT", "WPT", "GPT", "LPT", "WIT", "GIT", "OPRH", "WPRH", "GPRH"};
        const std::string tail = kw.substr(1);
        for (const char* t : sums) if (tail == t) return sum("W" + tail);
        return NaN;
    };
}

void Simulator::udqEvaluate(const SummaryPoint* sp, SummaryPoint* out, bool assignsOnly) {
    if (!udq_) return;
    static const SummaryPoint zero;
    const SummaryPoint& s = sp ? *sp : (summary_.empty() ? zero : summary_.back());
    udq::Source src;
    udqSource(s, src);
    udq_->evaluate(src, step_, assignsOnly);
    if (out) {
        out->udq.clear();
        for (const std::string& n : udq_->names()) {
            const udq::Value* v = udq_->get(n);
            if (!v) continue;
            if (v->kind == udq::Value::Scalar) out->udq.push_back({n, v->s});
            else {
                const auto& names = v->kind == udq::Value::Well ? src.wells : src.groups;
                for (size_t i = 0; i < v->v.size() && i < names.size(); ++i) out->udq.push_back({n + ":" + names[i], v->v[i]});
            }
        }
    }
}

void Simulator::applyUda() {
    if (!udq_) return;
    // C48 group UDA (GCONPROD / GCONINJE rate items given by an FU_ / GU_ name); undefined: the previous value stays
    auto gval = [&](const std::string& name, const std::string& group) {
        const auto it = std::find(udqGroups_.begin(), udqGroups_.end(), group);
        return udq_->groupValue(name, it == udqGroups_.end() ? size_t(-1) : size_t(it - udqGroups_.begin()));
    };
    for (auto& [gn, g] : gprod_)
        for (const auto& [target, name] : g.uda) {
            const double v = gval(name, gn);
            if (!std::isfinite(v)) continue;
            if (target == "ORAT") g.orat = v; else if (target == "WRAT") g.wrat = v; else if (target == "GRAT") g.grat = v;
            else if (target == "LRAT") g.lrat = v; else if (target == "RESV") g.resv = v;
        }
    for (auto& m : ginj_)
        for (auto& [gn, g] : m)
            for (const auto& [target, name] : g.uda) {
                const double v = gval(name, gn);
                if (!std::isfinite(v)) continue;
                if (target == "RATE") g.rate = v; else if (target == "RESV") g.resv = v; else if (target == "REIN") g.rein = v;
                else if (target == "VREP") g.vrep = v;
            }
    for (size_t i = 0; i < wells_.size(); ++i) {
        WellRT& w = wells_[i];
        if (w.own.uda.empty()) continue;
        for (const auto& [target, name] : w.own.uda) {
            const double v = udq_->wellValue(name, i);
            if (!std::isfinite(v)) continue;            // undefined: the previous target stays in force
            auto put = [&](WellControl& c) {
                if (target == "ORAT") c.orat = v; else if (target == "WRAT") c.wrat = v; else if (target == "GRAT") c.grat = v;
                else if (target == "LRAT") c.lrat = v; else if (target == "RESV") c.resv = v; else if (target == "BHP") c.bhp = v;
                else if (target == "THP") c.thp = v; else if (target == "RATE") c.rate = v;
            };
            put(w.own);
            if (!w.grpCtl) put(w.ctrl);
            if (target == "BHP") { w.bhpDeck = v; if (w.mode == Control::BHP && !w.grpCtl) w.bhp = v; }
        }
    }
}

// C32: ACTIONX COMPDAT adds completions to case_.wells during the run. The connection set of every well is rebuilt
// from the (modified) well input exactly as at setup - Peaceman factors, LGR / dual-porosity / EDFM targets, report-step
// windows - while the dynamic state of the wells (controls, BHP, mode, group / economic state) is kept; per-connection
// state follows its cell (new connections start from the cell values) and the Jacobian pattern is rebuilt (the cell rows
// gain the well columns of the new connections). Wells added by WELSPECS inside an action are not created at run time.
bool Simulator::rebuildActionConnections(bool quiet) {
    if (!msw_.empty()) {
        if (!quiet) log_("  Warning: ACTIONX adds connections in a run with multisegment wells - not supported by the ResSim engine yet (ignored; OPM Flow applies it)");
        return false;
    }
    if (distributedSolver() != nullptr || daOn_) {
        if (!quiet) log_("  Warning: ACTIONX adds connections in a domain-decomposed run - ignored by the ResSim engine (OPM Flow applies it)");
        return false;
    }
    const size_t nOld = wells_.size();
    std::vector<WellRT> old = std::move(wells_);
    try { buildWells(); }
    catch (const std::exception& e) {
        wells_ = std::move(old);
        if (!quiet) log_(std::string("  Warning: ACTIONX connections not added: ") + e.what());
        // restore the cell -> well lists of the old connection set
        cellWells_.assign(size_t(grid_.numActive()), {});
        for (size_t wi = 0; wi < wells_.size(); ++wi)
            for (size_t c = 0; c < wells_[wi].cell.size(); ++c) {
                auto& cw = cellWells_[size_t(wells_[wi].cell[c])];
                if (std::find(cw.begin(), cw.end(), int(wi)) == cw.end()) cw.push_back(int(wi));
            }
        return false;
    }
    if (wells_.size() > nOld) {
        // C33: wells created at run time (WELSPECS + COMPDAT in an ACTIONX body or a script) join the run shut; their
        // control records (applied from the next report step) open them. Per-well run-time arrays grow with them.
        std::string names;
        for (size_t i = nOld; i < wells_.size(); ++i) names += (names.empty() ? "" : ", ") + case_.wells[size_t(wells_[i].input)].name;
        if (!quiet) log_("  ACTIONX: " + std::to_string(wells_.size() - nOld) + " new well(s) created: " + names);
        wellNames_.clear(); for (const auto& w : case_.wells) wellNames_.push_back(w.name);
        if (!wCum_.empty()) wCum_.resize(wells_.size() * 5, 0.0);
        if (!pcInj_.empty()) pcInj_.resize(wells_.size(), 0.0);
        for (auto& t : trInj_) t.resize(wells_.size(), 0.0);
        if (!econLim_.empty()) econLim_.resize(wells_.size());
        if (!wtest_.empty()) wtest_.resize(wells_.size());
        wellChain_.clear();                               // group membership is rebuilt at the next report step
    }
    int added = 0;
    for (size_t wi = 0; wi < nOld && wi < wells_.size(); ++wi) {
        WellRT& fresh = wells_[wi];
        WellRT st = std::move(old[wi]);
        added += int(fresh.cell.size()) - int(st.cell.size());
        // per-connection dynamic values follow their cell (first match)
        auto follow = [&](const std::vector<double>& v, double def) {
            std::vector<double> r(fresh.cell.size(), def);
            if (v.size() != st.cell.size()) return r;
            for (size_t c = 0; c < fresh.cell.size(); ++c)
                for (size_t o = 0; o < st.cell.size(); ++o) if (st.cell[o] == fresh.cell[c]) { r[c] = v[o]; break; }
            return r;
        };
        if (!st.qwConn.empty()) st.qwConn = follow(st.qwConn, 0.0);
        if (!st.wmulConn.empty()) st.wmulConn = follow(st.wmulConn, 1.0);
        if (!st.zwConn.empty()) st.zwConn = follow(st.zwConn, 1.0);
        if (!st.dpFric.empty()) st.dpFric = follow(st.dpFric, 0.0);
        st.dpConn.clear();                              // re-initialised from the cell densities at the next assembly
        st.cell = fresh.cell; st.wi = fresh.wi; st.cdepth = fresh.cdepth; st.wiBase = fresh.wiBase; st.dlwi = fresh.dlwi;
        // C47: a defaulted WELSPECS reference depth follows the connections (OPM Well::updateRefDepth). A well whose
        // connections all arrive at run time (opm-tests ACTIONX_GCONINJE: injectors opened by ACTIONX) had the setup
        // value 0 -> a ~260 bar wellbore head and a stalled rate control ("time step below minimum")
        if (!(case_.wells[size_t(st.input)].refDepth > 0)) st.refDepth = fresh.refDepth;
        // C46: the run-time WELPI factor follows its connection record (same cell and window start); a new period of a cell with
        // a different WPIMULT factor (a WPIMULT split) keeps the cell's factor, a COMPDAT redefinition starts unscaled (OPM)
        {
            std::vector<double> fw(fresh.cell.size(), 1.0);
            for (size_t c = 0; c < fresh.cell.size(); ++c) {
                int same = -1, cellHit = -1;
                for (size_t o = 0; o < st.cell.size() && o < st.cWelpi.size(); ++o) {
                    if (st.cell[o] != fresh.cell[c]) continue;
                    if (st.cFrom[o] == fresh.cFrom[c]) { same = int(o); break; }
                    if (st.cFrom[o] <= fresh.cFrom[c] && (cellHit < 0 || st.cFrom[o] > st.cFrom[size_t(cellHit)])) cellHit = int(o);
                }
                if (same >= 0) fw[c] = st.cWelpi[size_t(same)];
                else if (cellHit >= 0 && st.cMult[size_t(cellHit)] != fresh.cMult[c]) fw[c] = st.cWelpi[size_t(cellHit)];
                fresh.wiBase[c] *= fw[c]; fresh.wi[c] *= fw[c];
            }
            st.cWelpi = fw;
        }
        st.cFrom = fresh.cFrom; st.cUntil = fresh.cUntil; st.cComp = fresh.cComp; st.cpos = fresh.cpos;
        st.cDef = fresh.cDef; st.cMult = fresh.cMult;
        st.rw = fresh.rw; st.re = fresh.re; st.plen = fresh.plen; st.cellWellIdx = fresh.cellWellIdx; st.rowCells = fresh.rowCells;
        fresh = std::move(st);
    }
    buildPattern();
    amgStale_ = true;                                   // the pressure hierarchy belongs to the old pattern
    if (!quiet) log_("  ACTIONX: connection set rebuilt (" + std::to_string(added) + " connection(s) added), Jacobian pattern updated");
    return true;
}

void Simulator::refreshConnectionWindows(bool quiet) {
    // C32: completions without a connection that are open now or later -> rebuild the connection set once
    bool missing = false;
    for (const WellRT& w : wells_) {
        const WellInput& in = case_.wells[size_t(w.input)];
        std::vector<char> used(in.completions.size(), 0);
        for (size_t c = 0; c < w.cComp.size(); ++c) if (size_t(w.cComp[c]) < used.size()) used[size_t(w.cComp[c])] = 1;
        for (size_t k = 0; k < in.completions.size(); ++k)
            if (!used[k] && in.completions[k].open && (in.completions[k].untilStep < 0 || in.completions[k].untilStep > step_)) missing = true;
    }
    if (missing || case_.wells.size() > wells_.size()) rebuildActionConnections(quiet);   // C33: also new wells
    for (WellRT& w : wells_) {
        const WellInput& in = case_.wells[size_t(w.input)];
        std::vector<char> used(in.completions.size(), 0);
        for (size_t c = 0; c < w.cell.size() && c < w.cComp.size(); ++c) {
            const size_t k = size_t(w.cComp[c]);
            if (k >= in.completions.size()) continue;
            used[k] = 1;
            const Completion& cp = in.completions[k];
            w.cFrom[c] = std::max(0, cp.fromStep);
            w.cUntil[c] = cp.open ? cp.untilStep : w.cFrom[c];   // a completion that never opens: empty window
            if (c < w.cMult.size() && w.cMult[c] > 0 && cp.mult != w.cMult[c]) {   // C46: WPIMULT changed this record in place
                const double r = cp.mult / w.cMult[c];
                w.wiBase[c] *= r; if (w.wi[c] > 0) w.wi[c] *= r;
                w.cMult[c] = cp.mult;
            }
        }
        if (quiet) continue;
        for (size_t k = 0; k < in.completions.size(); ++k)
            if (!used[k] && in.completions[k].open && (in.completions[k].untilStep < 0 || in.completions[k].untilStep > step_))
                log_("  Warning: ACTIONX opens a connection of " + in.name + " in an inactive or unknown cell (or the connection set could not be rebuilt) - ignored");
    }
}

void Simulator::actionApply(int step, const std::string& body, bool quiet) {
    // UDQ records inside the action body extend the quantity set from this step on
    const deck::RawDeck D = deck::parseDeckText("SCHEDULE\n" + body, "", "ACTIONX");
    std::string other;
    for (const auto& k : D.keywords) {
        if (k.name == "UDQ") {
            std::string err;
            if (udq_ && !udq_->add(deck::keywordBodyText(k), step, err) && !quiet) log_("  Warning: UDQ in ACTIONX: " + err);
        }
    }
    std::vector<std::string> msgs;
    applyScheduleText(case_, step, body, &msgs);
    if (!quiet) for (const auto& m : msgs) log_("  ACTIONX: " + m);
    refreshConnectionWindows(quiet);
}

std::vector<std::string> Simulator::applySchedule(const std::string& text) {
    std::vector<std::string> msgs;
    if (step_ >= int(case_.schedule.size())) { msgs.push_back("the schedule is finished - nothing applied"); return msgs; }
    const deck::RawDeck D = deck::parseDeckText("SCHEDULE\n" + text, "", "script");
    for (const auto& k : D.keywords) {
        if (k.name != "UDQ") continue;
        std::string err;
        if (udq_ && !udq_->add(deck::keywordBodyText(k), step_, err)) msgs.push_back("UDQ: " + err);
        else if (!udq_) msgs.push_back("UDQ: the deck defines no UDQ set - record ignored");
    }
    if (baseSchedule_.empty()) { baseSchedule_ = case_.schedule; baseWells_ = case_.wells; }   // restart replay base
    const size_t nW0 = case_.wells.size();
    applyScheduleText(case_, step_, text, &msgs);
    refreshConnectionWindows(true);
    for (size_t i = nW0; i < case_.wells.size(); ++i)
        msgs.push_back("new well " + case_.wells[i].name + (i < wells_.size() ? " created" : " defined but not added to the run (multisegment / distributed run)"));
    actStep_.push_back(step_); actText_.push_back(text);
    log_("Report step " + std::to_string(step_) + " (" + g(time_) + " days): script schedule input applied");
    return msgs;
}

void Simulator::actionsEndOfStep() {
    if (!actions_ || actions_->empty() || !udq_) return;
    const int evalStep = step_ - 1;                    // the report step that just ended
    if (step_ >= int(case_.schedule.size())) return;    // nothing left to act on
    // the conditions see the state of the last accepted time step (the UDQs were evaluated with it)
    static const SummaryPoint zero;
    const SummaryPoint& s = summary_.empty() ? zero : summary_.back();
    udq::Source src;
    udqSource(s, src);
    for (udq::Action& a : *actions_) {
        if (a.fromStep > evalStep) continue;
        if (a.runs >= a.maxRuns) continue;
        if (a.runs > 0 && time_ - a.lastRun < a.minWaitDays - 1e-9) continue;
        const udq::Trigger t = udq::evaluateAction(a, src, *udq_);
        if (!t.fired) continue;
        ++a.runs; a.lastRun = time_;
        std::vector<std::string> wn;
        for (size_t i : t.wells) wn.push_back(src.wells[i]);
        const std::string body = udq::expandBody(a.body, wn);
        std::string wl;
        for (const auto& w : wn) wl += (wl.empty() ? "" : " ") + w;
        log_("Report step " + std::to_string(evalStep + 1) + " (" + g(time_) + " days): ACTIONX " + a.name + " fires (run " +
             std::to_string(a.runs) + "/" + std::to_string(a.maxRuns) + ")" + (wl.empty() ? std::string() : " - wells " + wl));
        actionApply(step_, body, false);
        actStep_.push_back(step_); actText_.push_back(body);
        ActionEvent ev; ev.t = time_; ev.step = evalStep; ev.name = a.name; ev.wells = wl;
        actEvents_.push_back(ev);
    }
}

// ======================================================================== C37 PYACTION
int Simulator::numPyActions() const { return pyActions_ ? int(pyActions_->size()) : 0; }
std::string Simulator::pythonVersion() const { return pyHost_ ? pyHost_->version() : std::string(); }

void Simulator::pyActionSetup() {
    pyActions_.reset(); pyHost_.reset(); pyFailed_ = false;
    std::vector<std::pair<int, std::string>> bodies;
    for (const ExtraKeyword& e : case_.extraKeywords) if (e.name == "PYACTION") bodies.emplace_back(e.step, e.body);
    if (bodies.empty()) return;
    std::string err;
    auto acts = std::make_shared<std::vector<pyact::Action>>(pyact::parseActions(bodies, err));
    if (!err.empty()) log_("Warning: " + err);
    if (acts->empty()) return;
    for (auto& a : *acts) {
        std::filesystem::path fp = std::filesystem::u8path(a.file);
        if (fp.is_relative() && !case_.deckDir.empty()) fp = std::filesystem::u8path(case_.deckDir) / fp;
        a.file = fp.u8string();
        std::error_code ec;
        if (!std::filesystem::exists(fp, ec)) log_("Warning: PYACTION " + a.name + ": script " + a.file + " not found");
    }
    pyActions_ = acts;
    if (baseSchedule_.empty()) { baseSchedule_ = case_.schedule; baseWells_ = case_.wells; }   // restart replay base
    log_("PYACTION: " + std::to_string(acts->size()) + " Python action(s), run at the end of every report step (the Python process starts at the first call)");
}

void Simulator::pyActionsEndOfStep() {
    if (!pyActions_ || pyActions_->empty() || pyFailed_) return;
    const int evalStep = step_ - 1;
    if (step_ >= int(case_.schedule.size())) return;    // nothing left to act on
    bool any = false;
    for (const auto& a : *pyActions_) if (!a.done && a.fromStep <= evalStep) any = true;
    if (!any) return;
    if (!pyHost_) pyHost_ = std::make_shared<pyact::Host>();
    if (!pyHost_->running()) {
        std::string err;
        if (!pyHost_->start(case_.pythonExe, err)) {
            log_("Warning: PYACTION disabled - " + err);
            pyFailed_ = true;
            return;
        }
        log_("PYACTION: Python " + pyHost_->version() + " started");
    }
    static const SummaryPoint zero;
    const SummaryPoint& s = summary_.empty() ? zero : summary_.back();
    udq::Source src;
    udqSource(s, src);
    auto idxOf = [](const std::vector<std::string>& v, const std::string& n) -> long {
        for (size_t i = 0; i < v.size(); ++i) if (v[i] == n) return long(i);
        return -1;
    };
    // values requested by the scripts: [kind F|W|G, vector, entity]; UDQ names are served from the UDQ set
    auto query = [&](const pyact::Json& q) -> pyact::Json {
        if (q.t != pyact::Json::Arr || q.a.size() < 2) return pyact::Json();
        const std::string kind = q.a[0].s, name = q.a[1].s, ent = q.a.size() > 2 ? q.a[2].s : std::string();
        double v = NaN;
        try {
            if (kind == "F") {
                if (udq_ && udq_->has(name)) v = udq_->scalar(name);
                else if (src.field) v = src.field(name);
            } else if (kind == "W") {
                const long i = idxOf(src.wells, ent);
                if (i >= 0) {
                    if (udq_ && udq_->has(name)) v = udq_->wellValue(name, size_t(i));
                    else if (src.well) v = src.well(name, size_t(i));
                }
            } else if (kind == "G") {
                const long i = idxOf(src.groups, ent);
                if (i >= 0) {
                    if (udq_ && udq_->has(name)) v = udq_->groupValue(name, size_t(i));
                    else if (src.group) v = src.group(name, size_t(i));
                }
            }
        } catch (...) { v = NaN; }
        return std::isfinite(v) ? pyact::Json::num(v) : pyact::Json();
    };
    pyact::Json wl = pyact::Json::arr(), gl = pyact::Json::arr(), ws = pyact::Json::obj();
    for (size_t i = 0; i < src.wells.size(); ++i) {
        wl.a.push_back(pyact::Json::str(src.wells[i]));
        ws.set(src.wells[i], pyact::Json::str(i < wells_.size() && wells_[i].open ? "OPEN" : "SHUT"));
    }
    for (const auto& gname : src.groups) gl.a.push_back(pyact::Json::str(gname));
    pyact::Json start = pyact::Json::arr();
    start.a = {pyact::Json::num(case_.startYear), pyact::Json::num(case_.startMonth), pyact::Json::num(case_.startDay)};
    for (pyact::Action& a : *pyActions_) {
        if (a.done || a.fromStep > evalStep) continue;
        pyact::Json req = pyact::Json::obj();
        req.set("cmd", pyact::Json::str("run")).set("name", pyact::Json::str(a.name)).set("file", pyact::Json::str(a.file))
           .set("step", pyact::Json::num(evalStep + 1)).set("elapsed", pyact::Json::num(time_ * 86400.0)).set("days", pyact::Json::num(time_))
           .set("start", start).set("wells", wl).set("groups", gl).set("well_status", ws)
           .set("units", pyact::Json::str(case_.units == UnitSystem::Field ? "FIELD" : "METRIC")).set("title", pyact::Json::str(case_.title));
        pyact::Json resp; std::string err;
        const std::string head = "Report step " + std::to_string(evalStep + 1) + " (" + g(time_) + " days): PYACTION " + a.name;
        if (!pyHost_->call(req, query, resp, err)) {
            log_("Warning: " + head + " - " + err + "; PYACTION disabled for the rest of the run");
            pyFailed_ = true;
            return;
        }
        ++a.calls;
        if (const pyact::Json* lg = resp.get("log"))
            for (const auto& m : lg->a) if (m.t == pyact::Json::Arr && m.a.size() >= 2) log_("  PYACTION " + a.name + " [" + m.a[0].s + "] " + m.a[1].s);
        const pyact::Json* ok = resp.get("ok");
        if (!ok || !ok->b) {
            const pyact::Json* e = resp.get("error");
            std::string msg = e ? e->s : std::string("unknown error");
            while (!msg.empty() && (msg.back() == '\n' || msg.back() == '\r')) msg.pop_back();
            log_("Warning: " + head + " - Python error:\n" + msg);
            a.done = true;                                // a failing script is not called again
            continue;
        }
        const pyact::Json* res = resp.get("result");
        const bool legacyTrue = res && res->t == pyact::Json::Bool && res->b;
        if (a.runCount == "SINGLE") a.done = true;
        else if (a.runCount == "FIRST_TRUE" && legacyTrue) a.done = true;
        // ACTIONX callbacks (legacy run(..., actionx_callback)): the body of the named action with its well list
        std::string applied;
        if (const pyact::Json* ax = resp.get("actionx")) {
            for (const auto& c : ax->a) {
                if (c.t != pyact::Json::Arr || c.a.empty()) continue;
                const std::string an = c.a[0].s;
                std::vector<std::string> wn;
                if (c.a.size() > 1) for (const auto& w : c.a[1].a) wn.push_back(w.s);
                const udq::Action* hit = nullptr;
                if (actions_) for (const auto& x : *actions_) if (x.name == an) hit = &x;
                if (!hit) { log_("Warning: " + head + " - actionx_callback: no ACTIONX named " + an); continue; }
                applied += udq::expandBody(hit->body, wn);
                if (!applied.empty() && applied.back() != '\n') applied += '\n';
            }
        }
        if (const pyact::Json* kw = resp.get("keywords")) {
            applied += kw->s;
            if (!applied.empty() && applied.back() != '\n') applied += '\n';
        }
        bool blank = true;
        for (char c : applied) if (!std::isspace(static_cast<unsigned char>(c))) blank = false;
        if (blank) continue;
        log_(head + " applies schedule input (call " + std::to_string(a.calls) + ")");
        actionApply(step_, applied, false);
        actStep_.push_back(step_); actText_.push_back(applied);
        ActionEvent ev; ev.t = time_; ev.step = evalStep; ev.name = "PY:" + a.name; ev.notes = applied.substr(0, 200);
        actEvents_.push_back(ev);
    }
}

double Simulator::groupRateTarget(const std::string& g, const char* q) const {
    const auto it = gprod_.find(g);
    if (it == gprod_.end()) return -1;
    const std::string t(q);
    const GroupProdControl& c = it->second;
    return t == "ORAT" ? c.orat : t == "WRAT" ? c.wrat : t == "GRAT" ? c.grat : t == "LRAT" ? c.lrat : t == "RESV" ? c.resv : -1;
}

}  // namespace ressim
