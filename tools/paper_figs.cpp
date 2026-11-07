// C24: publication figures for research papers, rendered by ResSim Studio's own figure engine (app/src/figure.cpp, the
// same code as File > Export figure > Publication in the app).
//   paper_figs render <spec.txt> <outBase>      figure spec (written by research/papers/*/mkfigs.py) -> outBase.pdf / .svg / .csv
//   paper_figs faultrun <deck.DATA> <outBase>   runs the coupled fault-mechanics case in the ResSim engine (two-way and
//                                                one-way) and renders the results figure + outBase_numbers.txt
// Build (native, after tools/nb.sh):
//   g++ -std=c++17 -O2 -fopenmp -Iengine/include -Iapp/src tools/paper_figs.cpp app/src/figure.cpp $(ls build/native/*.o) -o /tmp/sd/pfig -lpthread -ldl
//
// Spec format (one item per line, '|' separates text fields):
//   FIGURE rows cols widthMM heightMM [legend: auto|tr|tl|br|bl|below|none]
//   PANEL title|xlabel|ylabel          LOGX  LOGY  LIMITS x0 x1 y0 y1  VLINE x  TEXT x y|text
//   SERIES label|color|width|kind|dash  (color = #rrggbb or auto:N; kind = line|scatter|markers|bars:W) then "x y" lines
#include "figure.h"
#include "ressim/case.h"
#include "ressim/io.h"
#include "ressim/simulator.h"
#include "ressim/faultmech.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace ressim;

static std::vector<std::string> split(const std::string& s, char d) {
    std::vector<std::string> o; std::string c; std::istringstream is(s);
    while (std::getline(is, c, d)) o.push_back(c);
    return o;
}
static uint32_t color(const std::string& c) {
    if (c.rfind("auto:", 0) == 0) return fig::distinctColor(std::stoi(c.substr(5)), 0);
    if (!c.empty() && c[0] == '#' && c.size() >= 7) {
        const int r = std::stoi(c.substr(1, 2), nullptr, 16), g = std::stoi(c.substr(3, 2), nullptr, 16), b = std::stoi(c.substr(5, 2), nullptr, 16);
        return fig::rgba(r, g, b);
    }
    return fig::rgba(0, 0, 0);
}
static fig::Style journalStyle(double wmm, double hmm, const std::string& leg) {
    fig::Style st;
    st.widthMM = wmm; st.heightMM = hmm; st.font = 0; st.fontPt = 8; st.lineWidthPt = 1.0;
    st.grid = true; st.minorTicks = true; st.titles = true; st.box = true; st.useLimits = true; st.layers = true;
    st.legend = leg == "tr" ? fig::Legend::InsideTopRight : leg == "tl" ? fig::Legend::InsideTopLeft : leg == "br" ? fig::Legend::InsideBottomRight
              : leg == "bl" ? fig::Legend::InsideBottomLeft : leg == "below" ? fig::Legend::Below : leg == "none" ? fig::Legend::None : fig::Legend::Auto;
    return st;
}
static bool writeAll(const fig::Figure& F, const fig::Style& st, const std::string& base) {
    std::string err;
    bool ok = fig::writePDF(F, st, base + ".pdf", &err);
    ok = fig::writeSVG(F, st, base + ".svg", &err) && ok;
    ok = fig::writeCSV(F, base + ".csv", &err) && ok;
    if (!ok) std::fprintf(stderr, "write %s: %s\n", base.c_str(), err.c_str());
    else std::printf("wrote %s.pdf / .svg / .csv (%zu panels)\n", base.c_str(), F.panels.size());
    return ok;
}

static int render(const std::string& spec, const std::string& base) {
    std::ifstream in(spec);
    if (!in) { std::fprintf(stderr, "cannot open %s\n", spec.c_str()); return 1; }
    fig::Figure F; F.name = base;
    double wmm = 160, hmm = 100; std::string leg = "auto";
    fig::Panel* P = nullptr; fig::Series* S = nullptr;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (line.rfind("FIGURE", 0) == 0) {
            std::istringstream is(line.substr(6)); is >> F.rows >> F.cols >> wmm >> hmm; if (!(is >> leg)) leg = "auto";
            continue;
        }
        if (line.rfind("PANEL ", 0) == 0) {
            const auto f = split(line.substr(6), '|');
            F.panels.emplace_back(); P = &F.panels.back(); S = nullptr;
            P->title = f.size() > 0 ? f[0] : ""; P->xlabel = f.size() > 1 ? f[1] : ""; P->ylabel = f.size() > 2 ? f[2] : "";
            continue;
        }
        if (!P) continue;
        if (line == "LOGX") { P->logX = true; continue; }
        if (line == "LOGY") { P->logY = true; continue; }
        if (line.rfind("LIMITS", 0) == 0) { std::istringstream is(line.substr(6)); is >> P->x0 >> P->x1 >> P->y0 >> P->y1; P->hasLimits = true; continue; }
        if (line.rfind("VLINE", 0) == 0) { P->vlines.push_back(std::stod(line.substr(5))); continue; }
        if (line.rfind("TEXT ", 0) == 0) {
            const auto f = split(line.substr(5), '|');
            std::istringstream is(f[0]); fig::Shape sh; sh.kind = fig::Shape::Text; double x = 0, y = 0; is >> x >> y;
            sh.x = {x}; sh.y = {y}; sh.text = f.size() > 1 ? f[1] : ""; sh.halign = 0; sh.keepColor = true; sh.stroke = fig::rgba(40, 40, 40);
            P->shapes.push_back(sh);
            continue;
        }
        if (line.rfind("SERIES ", 0) == 0) {
            const auto f = split(line.substr(7), '|');
            P->series.emplace_back(); S = &P->series.back();
            S->label = f.size() > 0 ? f[0] : "";
            S->legend = !S->label.empty() && S->label[0] != '#';
            S->color = color(f.size() > 1 ? f[1] : "auto:0");
            S->width = f.size() > 2 ? std::stof(f[2]) : 1.0f;
            const std::string kind = f.size() > 3 ? f[3] : "line";
            S->scatter = kind == "scatter"; S->markers = kind == "markers";
            if (kind.rfind("bars:", 0) == 0) { S->barWidth = std::stod(kind.substr(5)); S->barLine = fig::rgba(30, 30, 30); }
            S->dash = f.size() > 4 ? std::stoi(f[4]) : 0;
            continue;
        }
        if (S) { std::istringstream is(line); double x, y; if (is >> x >> y) { S->x.push_back(x); S->y.push_back(y); } }
    }
    return writeAll(F, journalStyle(wmm, hmm, leg), base) ? 0 : 1;
}

// ---- the coupled fault-mechanics run of the main model (C24) as a figure of the paper
struct RunOut { faultmech::Output o; std::vector<double> t, dpAcross, injBhp; double sec = 0; };
static bool runCase(const Case& c, RunOut& R) {
    Simulator sim; sim.setup(c); sim.initialize();
    faultmech::Monitor M; std::string err;
    if (!M.setup(c, sim.grid(), err)) { std::fprintf(stderr, "faultmech: %s\n", err.c_str()); return false; }
    const Grid& G = sim.grid();
    const int aW = G.activeOf[size_t(14 + G.nx * (9 + G.ny * 2))], aE = G.activeOf[size_t(15 + G.nx * (9 + G.ny * 2))];
    auto rec = [&]() {
        const auto& P = sim.snapshots().back().pressure;
        R.t.push_back(sim.time()); R.dpAcross.push_back(P[size_t(aW)] - P[size_t(aE)]);
    };
    M.step(sim.time(), sim.snapshots().back().pressure); rec();
    while (sim.runReportStep()) { M.step(sim.time(), sim.snapshots().back().pressure); M.applyCoupling(sim); rec(); }
    M.finish(nullptr, [](const std::string& s) { std::printf("%s\n", s.c_str()); });
    R.o = M.out();
    return true;
}
static int faultrun(const std::string& deck, const std::string& base) {
    Case c = importEclipseDeck(deck);
    RunOut two, one;
    if (!runCase(c, two)) return 1;
    Case c1 = c; for (auto& r : c1.faultMech.faults) { r.permFactor = 1; r.rsfLayer = 0; }
    if (!runCase(c1, one)) return 1;
    const auto& S2 = two.o.series[0]; const auto& S1 = one.o.series[0];
    const double tShut = 500;
    fig::Figure F; F.rows = 2; F.cols = 2;
    const uint32_t cTwo = fig::rgba(0, 114, 178), cOne = fig::rgba(213, 94, 0), cGrey = fig::rgba(110, 110, 110);
    auto ser = [](const std::string& l, const std::vector<double>& x, const std::vector<double>& y, uint32_t col, float w, int dash = 0, bool sc = false) {
        fig::Series s; s.label = l; s.x = x; s.y = y; s.color = col; s.width = w; s.dash = dash; s.scatter = sc; return s;
    };
    {   fig::Panel P; P.title = "(a) Coulomb stress change on fault F1"; P.xlabel = "time (days)"; P.ylabel = "max dCFS (bar)";
        P.series.push_back(ser("two-way coupled", S2.t, S2.maxDcfs, cTwo, 1.6f));
        P.series.push_back(ser("one-way", S1.t, S1.maxDcfs, cOne, 1.4f, 1));
        P.vlines.push_back(tShut); F.panels.push_back(P); }
    {   fig::Panel P; P.title = "(b) Reactivated fault area (CFF >= 0)"; P.xlabel = "time (days)"; P.ylabel = "area (%)";
        std::vector<double> a2(S2.reactFrac), a1(S1.reactFrac); for (double& v : a2) v *= 100; for (double& v : a1) v *= 100;
        P.series.push_back(ser("two-way coupled", S2.t, a2, cTwo, 1.6f));
        P.series.push_back(ser("one-way", S1.t, a1, cOne, 1.4f, 1));
        P.vlines.push_back(tShut); F.panels.push_back(P); }
    {   fig::Panel P; P.title = "(c) Induced events (M >= -1)"; P.xlabel = "time (days)"; P.ylabel = "cumulative events";
        P.series.push_back(ser("expected, two-way", S2.t, S2.expected, cTwo, 1.6f));
        P.series.push_back(ser("catalogue, two-way", S2.t, S2.events, cTwo, 1.0f, 2));
        P.series.push_back(ser("expected, one-way", S1.t, S1.expected, cOne, 1.4f, 1));
        P.vlines.push_back(tShut); F.panels.push_back(P); }
    {   fig::Panel P; P.title = "(d) Pressure jump across F1 (injector row)"; P.xlabel = "time (days)"; P.ylabel = "p_W - p_E (bar)";
        P.series.push_back(ser("two-way coupled", two.t, two.dpAcross, cTwo, 1.6f));
        P.series.push_back(ser("one-way", one.t, one.dpAcross, cOne, 1.4f, 1));
        P.vlines.push_back(tShut); F.panels.push_back(P); }
    (void)cGrey;
    if (!writeAll(F, journalStyle(170, 115, "auto"), base)) return 1;
    // magnitude-time figure of the stochastic catalogue
    fig::Figure G; G.rows = 1; G.cols = 2;
    {   fig::Panel P; P.title = "(a) Event magnitudes, two-way run"; P.xlabel = "time (days)"; P.ylabel = "Mw";
        std::vector<double> t, m, tr, mr;
        for (const auto& e : two.o.events) { (e.kind == 0 ? t : tr).push_back(e.t); (e.kind == 0 ? m : mr).push_back(e.mw); }
        P.series.push_back(ser("Dieterich + G-R sample", t, m, cTwo, 1.0f, 0, true));
        if (!tr.empty()) P.series.push_back(ser("rate-and-state rupture", tr, mr, fig::rgba(204, 121, 167), 1.0f, 0, true));
        P.vlines.push_back(tShut); G.panels.push_back(P); }
    {   fig::Panel P; P.title = "(b) Magnitude-frequency"; P.xlabel = "Mw"; P.ylabel = "N (>= Mw)"; P.logY = true;
        std::vector<double> m; for (const auto& e : two.o.events) if (e.kind == 0) m.push_back(e.mw);
        std::sort(m.begin(), m.end());
        std::vector<double> x, y; for (size_t i = 0; i < m.size(); ++i) { x.push_back(m[i]); y.push_back(double(m.size() - i)); }
        P.series.push_back(ser("catalogue", x, y, cTwo, 1.0f, 0, true));
        if (!m.empty()) {
            const double b = c.faultMech.faults[0].bValue, mm = c.faultMech.faults[0].mMin;
            P.series.push_back(ser("input b = 1", {mm, m.back()}, {double(m.size()), double(m.size()) * std::pow(10.0, -b * (m.back() - mm))}, cOne, 1.4f, 1));
        }
        G.panels.push_back(P); }
    if (!writeAll(G, journalStyle(170, 62, "auto"), base + "_mw")) return 1;
    // numbers for the text
    std::ofstream o(base + "_numbers.txt");
    auto mx = [](const std::vector<double>& v) { return v.empty() ? 0.0 : *std::max_element(v.begin(), v.end()); };
    double dpi2 = 0, dpi1 = 0; for (size_t i = 0; i < two.t.size(); ++i) if (two.t[i] <= tShut) dpi2 = std::max(dpi2, two.dpAcross[i]);
    for (size_t i = 0; i < one.t.size(); ++i) if (one.t[i] <= tShut) dpi1 = std::max(dpi1, one.dpAcross[i]);
    size_t k100 = 0; for (size_t i = 0; i < S2.t.size(); ++i) if (std::fabs(S2.t[i] - 100) < std::fabs(S2.t[k100] - 100)) k100 = i;
    size_t k100b = 0; for (size_t i = 0; i < S1.t.size(); ++i) if (std::fabs(S1.t[i] - 100) < std::fabs(S1.t[k100b] - 100)) k100b = i;
    size_t kS = 0; for (size_t i = 0; i < S2.t.size(); ++i) if (S2.t[i] <= tShut + 1e-9) kS = i;
    int nD = 0, nR = 0; double mwMax = -9; for (const auto& e : two.o.events) { (e.kind == 0 ? nD : nR)++; mwMax = std::max(mwMax, double(e.mw)); }
    char b[1024];
    std::snprintf(b, sizeof b,
                  "faces %zu\nmax_dCFS_two %.3f\nmax_dCFS_one %.3f\nreact_100d_two %.1f\nreact_100d_one %.1f\nexpected_shutin_two %.2f\nexpected_end_two %.2f\nexpected_end_one %.2f\n"
                  "post_shutin_fraction_two %.3f\ncatalogue_dieterich %d\ncatalogue_rsf %d\nmw_max %.2f\nmw_cap %.2f\ndp_across_two %.2f\ndp_across_one %.2f\ntwo_way_updates %ld\n"
                  "critical_t0 %d\n",
                  two.o.faces.size(), mx(S2.maxDcfs), mx(S1.maxDcfs), 100 * S2.reactFrac[k100], 100 * S1.reactFrac[k100b], S2.expected[kS], S2.expected.back(), S1.expected.back(),
                  S2.expected.back() > 0 ? 1 - S2.expected[kS] / S2.expected.back() : 0.0, nD, nR, mwMax, two.o.mMax.empty() ? 0.0 : two.o.mMax[0], dpi2, dpi1, two.o.twoWayUpdates,
                  two.o.criticalAtStart);
    o << b;
    for (const auto& r : two.o.rsf) { std::snprintf(b, sizeof b, "rsf layer %d faces %d hstar %.1f Lb %.1f dx %.2f events %d err %s\n", r.layer, r.faces, r.hStar, r.Lb, r.dx, r.events, r.error.c_str()); o << b; }
    std::printf("%s", b);
    return 0;
}

int main(int argc, char** argv) {
    if (argc >= 4 && std::string(argv[1]) == "render") return render(argv[2], argv[3]);
    if (argc >= 4 && std::string(argv[1]) == "faultrun") return faultrun(argv[2], argv[3]);
    std::fprintf(stderr, "usage: paper_figs render <spec> <outBase> | faultrun <deck> <outBase>\n");
    return 2;
}
