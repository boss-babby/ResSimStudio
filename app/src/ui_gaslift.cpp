// C35 gas-lift optimisation editor (Wells page): LIFTOPT, WLIFTOPT, GLIFTOPT records per report step, the readiness of
// the producers (THP control + VFPPROD table with an ALQ axis), and the lift gas allocation of the last run.
// The keywords are kept verbatim in Case::extraKeywords (what the engine reads and the deck writer exports); a keyword
// is only re-written when it is edited here, so imported decks keep their layout and comments otherwise.
#include "app.h"
#include "widgets.h"
#include "implot.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <sstream>
#include <string>

using namespace ressim;

namespace {
const double kUnset = std::nan("");
bool unset(double v) { return std::isnan(v); }

// records of a keyword body: '/'-terminated, '--' comments removed, N*v expanded, '1*' -> "*"
std::vector<std::vector<std::string>> glRecords(const std::string& body) {
    std::vector<std::vector<std::string>> out;
    std::vector<std::string> cur;
    std::istringstream is(body);
    std::string line;
    while (std::getline(is, line)) {
        const size_t cm = line.find("--");
        if (cm != std::string::npos) line.erase(cm);
        std::string tok;
        bool q = false;
        auto flush = [&] {
            if (tok.empty()) return;
            const size_t st = tok.find('*');
            if (st != std::string::npos && tok[0] != '\'') {
                const int n = st > 0 ? std::atoi(tok.substr(0, st).c_str()) : 1;
                const std::string v = tok.substr(st + 1);
                for (int k = 0; k < std::max(1, n); ++k) cur.push_back(v.empty() ? "*" : v);
            } else cur.push_back(tok);
            tok.clear();
        };
        for (char ch : line) {
            if (ch == '\'') { q = !q; tok += ch; continue; }
            if (!q && ch == '/') { flush(); out.push_back(cur); cur.clear(); continue; }
            if (!q && (ch == ' ' || ch == '\t' || ch == '\r')) { flush(); continue; }
            tok += ch;
        }
        flush();
    }
    return out;
}
std::string unq(std::string s) { if (s.size() >= 2 && s.front() == '\'') s = s.substr(1, s.size() - 2); return s; }
double numOr(const std::vector<std::string>& r, size_t i, double d) {
    if (i >= r.size() || r[i] == "*" || r[i].empty()) return d;
    char* e = nullptr; const double v = std::strtod(r[i].c_str(), &e);
    return (e && *e == 0) ? v : d;
}
bool yesOr(const std::vector<std::string>& r, size_t i, bool d) {
    if (i >= r.size() || r[i] == "*") return d;
    std::string s = unq(r[i]);
    for (char& c : s) c = char(std::toupper(static_cast<unsigned char>(c)));
    return s.rfind("Y", 0) == 0;
}
std::string numTxt(double v) { if (unset(v)) return "1*"; char b[32]; std::snprintf(b, sizeof b, "%.6g", v); return b; }

struct LiftRow { int step = 0; double inc = 0, eco = 0, interval = 0; bool allIts = true; };
struct WellRow { int step = 0; std::string well; bool opt = true; double maxRate = kUnset, weight = kUnset, minRate = kUnset, gasW = kUnset; bool extra = false; };
struct GroupRow { int step = 0; std::string group; double maxLift = kUnset, maxTotal = kUnset; };

template <class R, class F> std::vector<R> readKw(const Case& c, const char* name, F parse) {
    std::vector<R> v;
    for (const auto& k : c.extraKeywords)
        if (k.name == name)
            for (const auto& r : glRecords(k.body)) { if (r.empty()) continue; R x; x.step = k.step; if (parse(r, x)) v.push_back(x); }
    std::stable_sort(v.begin(), v.end(), [](const R& a, const R& b) { return a.step < b.step; });
    return v;
}
template <class R, class F> void writeKw(Case& c, const char* name, const std::vector<R>& v, F line) {
    c.extraKeywords.erase(std::remove_if(c.extraKeywords.begin(), c.extraKeywords.end(), [&](const ExtraKeyword& k) { return k.name == name; }),
                          c.extraKeywords.end());
    std::map<int, std::string> body;
    for (const R& x : v) body[x.step] += line(x) + " /\n";
    for (const auto& [st, b] : body) {
        ExtraKeyword k; k.name = name; k.section = DeckSection::Schedule; k.step = st; k.beforeControls = false;
        k.body = b + (std::string(name) == "LIFTOPT" ? "" : "/\n");
        k.note = "gas lift (C35 editor)";
        c.extraKeywords.push_back(k);
    }
}
}  // namespace

void App::editGasLift() {
    Case& C = case_;
    bool ch = false;
    const bool fld = C.units == UnitSystem::Field;
    const char* gu = fld ? "Mscf/d" : "sm3/d";
    ui::hint("Gas-lift optimisation (ECLIPSE / OPM LIFTOPT, WLIFTOPT, GLIFTOPT) in the own engine: at every time step the lift "
             "gas is allocated in whole increments to the THP-controlled producer with the largest weighted incremental "
             "gradient  w * dQoil / (dALQ + wg * dQgas)  while it exceeds the minimum economic gradient, within the well "
             "maximum, the group lift gas supply / total gas limits (GLIFTOPT) and the group production targets (GCONPROD). "
             "The ALQ enters the VFPPROD lift curves (ENGINES_AND_MATH 6.30). OPM Flow and ECLIPSE read the same keywords.");
    // ---- readiness of the producers: THP limit + VFP table with an ALQ axis
    {
        std::map<std::string, const WellControl*> last;   // latest control record per well
        for (const auto& S : C.schedule) for (const auto& wc : S.controls) last[wc.well] = &wc;
        int ready = 0, prod = 0;
        std::string missing;
        for (const auto& w : C.wells) {
            if (w.type != WellType::Producer) continue;
            ++prod;
            auto it = last.find(w.name);
            const WellControl* wc = it == last.end() ? nullptr : it->second;
            const VfpTable* t = nullptr;
            if (wc && wc->vfp > 0) for (const auto& v : C.vfp) if (!v.inj && v.id == wc->vfp) t = &v;
            if (wc && wc->thp > 0 && t && t->alqV.size() > 1) ++ready;
            else missing += (missing.empty() ? "" : ", ") + w.name;
        }
        if (prod) {
            ImGui::Text("Producers ready for gas lift: %d of %d", ready, prod);
            ImGui::SameLine(); ui::help("A producer needs a THP limit (WCONPROD item 10) and a VFPPROD table (item 11) whose ALQ axis "
                                        "has more than one value; the ALQ type should be GRAT (lift gas rate) or blank. The maximum "
                                        "lift gas rate defaults to the largest ALQ of the table.");
            if (!missing.empty()) ImGui::TextColored(col4(Pal::Warn), "not ready: %s", missing.c_str());
        }
        std::set<std::string> types;
        for (const auto& v : C.vfp) if (!v.inj && v.alqV.size() > 1) { std::string a = unq(v.alq); types.insert(a.empty() ? "blank" : a); }
        if (!types.empty()) {
            std::string s; for (auto& t : types) s += (s.empty() ? "" : ", ") + t;
            ImGui::TextDisabled("VFPPROD ALQ types in the case: %s", s.c_str());
        }
    }
    // the start day of each report step (for the step columns)
    std::vector<double> day(C.schedule.size() + 1, 0.0);
    for (size_t i = 0; i < C.schedule.size(); ++i) day[i + 1] = day[i] + C.schedule[i].length;
    const int nst = std::max(0, int(C.schedule.size()) - 1);
    auto dayOf = [day](int st) { return day[size_t(std::clamp(st, 0, int(day.size()) - 1))]; };

    // ---- LIFTOPT
    ImGui::SeparatorText("LIFTOPT - optimisation settings");
    {
        std::vector<LiftRow> rows = readKw<LiftRow>(C, "LIFTOPT", [](const std::vector<std::string>& r, LiftRow& x) {
            x.inc = numOr(r, 0, 0); x.eco = numOr(r, 1, 0); x.interval = numOr(r, 2, 0); x.allIts = yesOr(r, 3, true); return true; });
        ui::RecordTable<LiftRow> t;
        t.integer("Step", [](LiftRow& x) -> int& { return x.step; }, 1)
         .computed("Day", [dayOf](const LiftRow& x) { return dayOf(x.step); })
         .num(fld ? "Increment (Mscf/d, <= 0 off)" : "Increment (sm3/d, <= 0 off)", [](LiftRow& x) -> double& { return x.inc; })
         .num("Min economic gradient (oil / gas)", [](LiftRow& x) -> double& { return x.eco; })
         .num("Min interval (days)", [](LiftRow& x) -> double& { return x.interval; })
         .flag("All Newton iterations", [](LiftRow& x) -> bool& { return x.allIts; });
        t.onChange = [nst](LiftRow& x) { x.step = std::clamp(x.step, 0, nst); x.eco = std::max(0.0, x.eco); x.interval = std::max(0.0, x.interval); };
        if (t.draw("liftopt", rows, ImGui::GetFontSize() * (3.0f + float(std::min<size_t>(rows.size(), 6))))) {
            writeKw(C, "LIFTOPT", rows, [](const LiftRow& x) {
                return "  " + numTxt(x.inc) + "  " + numTxt(x.eco) + "  " + numTxt(x.interval) + "  " + (x.allIts ? "YES" : "NO"); });
            ch = true;
        }
        if (rows.empty() && ui::accentButton("Switch gas-lift optimisation on")) {
            double amax = 0;
            for (const auto& v : C.vfp) if (!v.inj && !v.alqV.empty()) amax = std::max(amax, v.alqV.back());
            LiftRow x; x.inc = amax > 0 ? amax / 10 : (fld ? 1000 : 3e4);
            rows.push_back(x);
            writeKw(C, "LIFTOPT", rows, [](const LiftRow& r) { return "  " + numTxt(r.inc) + "  " + numTxt(r.eco) + "  " + numTxt(r.interval) + "  YES"; });
            if (readKw<WellRow>(C, "WLIFTOPT", [](const std::vector<std::string>&, WellRow&) { return true; }).empty()) {
                std::vector<WellRow> wr;
                for (const auto& w : C.wells) if (w.type == WellType::Producer) { WellRow r; r.well = w.name; wr.push_back(r); }
                writeKw(C, "WLIFTOPT", wr, [](const WellRow& r) { return " '" + r.well + "'  YES"; });
            }
            notify("Gas lift: LIFTOPT + WLIFTOPT for all producers added - producers need THP control and a VFPPROD table with an ALQ axis");
            ch = true;
        }
    }
    // ---- WLIFTOPT
    ImGui::SeparatorText("WLIFTOPT - wells");
    {
        std::vector<WellRow> rows = readKw<WellRow>(C, "WLIFTOPT", [](const std::vector<std::string>& r, WellRow& x) {
            x.well = unq(r[0]); x.opt = yesOr(r, 1, true);
            x.maxRate = numOr(r, 2, kUnset); x.weight = numOr(r, 3, kUnset); x.minRate = numOr(r, 4, kUnset); x.gasW = numOr(r, 5, kUnset);
            x.extra = yesOr(r, 6, false); return !x.well.empty(); });
        ui::RecordTable<WellRow> t;
        auto isU = [](double& (*f)(WellRow&)) { return [f](WellRow& x) { return unset(f(x)); }; };
        auto mk = [](double& (*f)(WellRow&)) { return [f](WellRow& x) { f(x) = kUnset; }; };
        double& (*fMax)(WellRow&) = [](WellRow& x) -> double& { return x.maxRate; };
        double& (*fW)(WellRow&) = [](WellRow& x) -> double& { return x.weight; };
        double& (*fMin)(WellRow&) = [](WellRow& x) -> double& { return x.minRate; };
        double& (*fG)(WellRow&) = [](WellRow& x) -> double& { return x.gasW; };
        static std::string hMax, hMin;
        hMax = std::string("Max lift gas (") + gu + ")"; hMin = std::string("Min lift gas (") + gu + ", < 0 flow)";
        t.integer("Step", [](WellRow& x) -> int& { return x.step; }, 1)
         .text("Well / pattern", [](WellRow& x) -> std::string& { return x.well; })
         .flag("Optimise", [](WellRow& x) -> bool& { return x.opt; })
         .optNum(hMax.c_str(), fMax, isU(fMax), mk(fMax), "table max")
         .optNum("Weight", fW, isU(fW), mk(fW), "1")
         .optNum(hMin.c_str(), fMin, isU(fMin), mk(fMin), "0")
         .optNum("Gas weight", fG, isU(fG), mk(fG), "0")
         .flag("Extra lift gas", [](WellRow& x) -> bool& { return x.extra; });
        t.onChange = [nst](WellRow& x) { x.step = std::clamp(x.step, 0, nst); if (!unset(x.maxRate)) x.maxRate = std::max(0.0, x.maxRate); };
        t.newRow = [&](int at) {
            WellRow x = rows.empty() ? WellRow{} : rows[size_t(std::clamp(at, 0, int(rows.size()) - 1))];
            if (x.well.empty()) for (const auto& w : C.wells) if (w.type == WellType::Producer) { x.well = w.name; break; }
            return x;
        };
        if (t.draw("wliftopt", rows, ImGui::GetFontSize() * (3.0f + float(std::min<size_t>(rows.size(), 8))))) {
            writeKw(C, "WLIFTOPT", rows, [](const WellRow& x) {
                std::string s = " '" + x.well + "'  " + (x.opt ? "YES" : "NO") + "  " + numTxt(x.maxRate) + "  " + numTxt(x.weight) + "  " +
                                numTxt(x.minRate) + "  " + numTxt(x.gasW);
                if (x.extra) s += "  YES";
                return s; });
            ch = true;
        }
        ui::hint("Optimise = NO with a maximum fixes the lift gas at that rate (NO 0.0 switches lift gas off); NO without a maximum "
                 "keeps the current rate. Weight scales the oil gain of the well (priority). 'Well / pattern' accepts OP* style "
                 "patterns.");
    }
    // ---- GLIFTOPT
    ImGui::SeparatorText("GLIFTOPT - group limits");
    {
        std::vector<GroupRow> rows = readKw<GroupRow>(C, "GLIFTOPT", [](const std::vector<std::string>& r, GroupRow& x) {
            x.group = unq(r[0]); x.maxLift = numOr(r, 1, kUnset); x.maxTotal = numOr(r, 2, kUnset); return !x.group.empty(); });
        std::vector<std::string> groups = {"FIELD"};
        for (const auto& w : C.wells) if (!w.group.empty() && std::find(groups.begin(), groups.end(), w.group) == groups.end()) groups.push_back(w.group);
        for (const auto& k : C.extraKeywords)
            if (k.name == "GRUPTREE")
                for (const auto& r : glRecords(k.body))
                    for (size_t i = 0; i < r.size() && i < 2; ++i) {
                        const std::string g = unq(r[i]);
                        if (!g.empty() && g != "*" && std::find(groups.begin(), groups.end(), g) == groups.end()) groups.push_back(g);
                    }
        static std::string hL, hT;
        hL = std::string("Max lift gas supply (") + gu + ")"; hT = std::string("Max total gas, produced + lift (") + gu + ")";
        ui::RecordTable<GroupRow> t;
        t.integer("Step", [](GroupRow& x) -> int& { return x.step; }, 1)
         .choiceStr("Group", [](GroupRow& x) -> std::string& { return x.group; }, groups)
         .optNum(hL.c_str(), [](GroupRow& x) -> double& { return x.maxLift; }, [](GroupRow& x) { return unset(x.maxLift); },
                 [](GroupRow& x) { x.maxLift = kUnset; }, "no limit")
         .optNum(hT.c_str(), [](GroupRow& x) -> double& { return x.maxTotal; }, [](GroupRow& x) { return unset(x.maxTotal); },
                 [](GroupRow& x) { x.maxTotal = kUnset; }, "no limit");
        t.onChange = [nst](GroupRow& x) { x.step = std::clamp(x.step, 0, nst); };
        t.newRow = [&](int at) {
            GroupRow x = rows.empty() ? GroupRow{} : rows[size_t(std::clamp(at, 0, int(rows.size()) - 1))];
            if (x.group.empty()) x.group = groups.size() > 1 ? groups[1] : "FIELD";
            return x;
        };
        if (t.draw("gliftopt", rows, ImGui::GetFontSize() * (3.0f + float(std::min<size_t>(rows.size(), 5))))) {
            writeKw(C, "GLIFTOPT", rows, [](const GroupRow& x) { return " '" + x.group + "'  " + numTxt(x.maxLift) + "  " + numTxt(x.maxTotal); });
            ch = true;
        }
        if (rows.empty() && ImGui::SmallButton("+ Group lift gas limit")) {
            GroupRow x; x.group = groups.size() > 1 ? groups[1] : "FIELD"; rows.push_back(x);
            writeKw(C, "GLIFTOPT", rows, [](const GroupRow& r) { return " '" + r.group + "'  " + numTxt(r.maxLift) + "  " + numTxt(r.maxTotal); });
            ch = true;
        }
    }
    // ---- allocation of the last run
    {
        std::vector<double> tt, fg;
        std::vector<std::vector<double>> wl;
        std::vector<std::string> wn;
        double tot = 0;
        {
            std::lock_guard<std::mutex> l(runner_.mutex());
            const RunResults& R = runner_.results();
            if (R.valid && R.summary.size() > 1) {
                std::vector<size_t> idx;
                for (size_t w = 0; w < R.wellNames.size(); ++w)
                    for (const auto& sp : R.summary) if (w < sp.wells.size() && sp.wells[w].lift > 0) { idx.push_back(w); break; }
                for (size_t w : idx) { wn.push_back(R.wellNames[w]); wl.emplace_back(); }
                for (const auto& sp : R.summary) {
                    tt.push_back(sp.time); fg.push_back(sp.fglir); tot = sp.fglit;
                    for (size_t j = 0; j < idx.size(); ++j) wl[j].push_back(idx[j] < sp.wells.size() ? sp.wells[idx[j]].lift : 0.0);
                }
            }
        }
        if (!fg.empty() && *std::max_element(fg.begin(), fg.end()) > 0) {
            ImGui::SeparatorText("Lift gas allocation (last run)");
            ImGui::Text("Total lift gas injected (FGLIT): %.4g %s", tot, fld ? "Mscf" : "sm3");
            if (ImPlot::BeginPlot("Lift gas", ImVec2(-1, ImGui::GetFontSize() * 15))) {
                ImPlot::SetupAxes("days", gu, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlotSpec a; a.LineColor = col4(Pal::Accent); a.LineWeight = 2.5f; a.Flags = ImPlotStairsFlags_PreStep;   // the value of t[i] is the rate over (t[i-1], t[i]]
                ImPlot::PlotStairs("FGLIR (field, efficiency weighted)", tt.data(), fg.data(), int(tt.size()), a);
                for (size_t j = 0; j < wn.size(); ++j) {
                    ImPlotSpec b; b.LineWeight = 1.2f; b.Flags = ImPlotStairsFlags_PreStep;
                    ImPlot::PlotStairs(("WGLIR " + wn[j]).c_str(), tt.data(), wl[j].data(), int(tt.size()), b);
                }
                ImPlot::EndPlot();
            }
            ui::hint("Summary vectors FGLIR / FGLIT and WGLIR / WGLIT per well (Results > custom plots, ECLIPSE summary export).");
        }
    }
    if (ch) markDirty();
}
