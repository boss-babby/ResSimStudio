// C48 tests: SUBGRID LEAKC (gas-aware sub-cell storage), defaulted BHP reference depth with SHUT connections (OPM
// Well::updateRefDepth), group UDA (GCONPROD rate given by a UDQ name), comp / thermal cell snapshots (CLI -o)
#include "ressim/cell_snapshot.h"
#include "ressim/comp_deck.h"
#include "ressim/therm_deck.h"
#include <filesystem>
#include <fstream>

static std::string c48Replace(std::string s, const std::string& a, const std::string& b, bool* ok = nullptr) {
    const size_t q = s.find(a);
    if (q == std::string::npos) { if (ok) *ok = false; return s; }
    s.replace(q, a.size(), b);
    return s;
}
static std::string c48Write(const std::string& name, const std::string& text) {
    namespace fs = std::filesystem;
    const std::string dir = (fs::temp_directory_path() / "ressim_c48").string();
    fs::create_directories(dir);
    const std::string p = dir + "/" + name;
    std::ofstream(p, std::ios::binary) << text;
    return p;
}

static void testC48Leakc() {
    std::printf("== C48 SUBGRID LEAKC: parse, export, re-import (H2e gas-aware sub-cell storage)\n");
    const std::string path = findExample("HFGROW_INJECTION.DATA");
    if (path.empty()) { std::printf("    (examples/HFGROW_INJECTION.DATA not found - skipped)\n"); return; }
    bool ok = true;
    const std::string d = c48Replace(c32ReadText(path), "-- RESSIM: HFGROW F1 1.0 340 0 20\n", "-- RESSIM: HFGROW F1 1.0 340 0 20 SUBGRID LEAKC\n", &ok);
    check(ok, "C48: HFGROW_INJECTION layout (HFGROW line)");
    if (!ok) return;
    try {
        const Case c = importEclipseDeck(c48Write("LEAKC.DATA", d), nullptr, nullptr);
        check(!c.grid.hfracs.empty() && c.grid.hfracs[0].subLeak == 2, "C48: SUBGRID LEAKC -> subLeak = 2");
        const std::string out = c48Write("LEAKC_OUT.DATA", "");
        exportEclipseDeck(c, out);
        const std::string txt = c32ReadText(out);
        check(txt.find("SUBGRID LEAKC") != std::string::npos, "C48: the exporter writes SUBGRID LEAKC");
        const Case r = importEclipseDeck(out, nullptr, nullptr);
        check(!r.grid.hfracs.empty() && r.grid.hfracs[0].subLeak == 2, "C48: LEAKC survives the round trip");
        // plain SUBGRID LEAK stays 1
        const Case p = importEclipseDeck(c48Write("LEAK.DATA", c48Replace(d, "SUBGRID LEAKC", "SUBGRID LEAK")), nullptr, nullptr);
        check(!p.grid.hfracs.empty() && p.grid.hfracs[0].subLeak == 1, "C48: SUBGRID LEAK stays the H2b scheme (subLeak = 1)");
    } catch (const std::exception& e) { check(false, (std::string("C48: LEAKC deck: ") + e.what()).c_str()); }
}

static void testC48RefDepthShut() {
    std::printf("== C48 defaulted BHP reference depth counts connections defined SHUT (OPM Well::updateRefDepth)\n");
    const std::string path = findExample("SPE1.DATA");
    if (path.empty()) { std::printf("    (examples/SPE1.DATA not found - skipped)\n"); return; }
    bool ok = true;
    std::string base = c48Replace(c32ReadText(path), " 'PROD' 'PROD' 10 10 8400 OIL /\n", " 'PROD' 'PROD' 10 10 1* OIL /\n", &ok);
    const std::string shut = c48Replace(base, " 'PROD' 10 10 3 3 OPEN", " 'PROD' 10 10 1 1 SHUT 1* 1* 0.5 1* 0 1* Z /\n 'PROD' 10 10 3 3 OPEN", &ok);
    check(ok, "C48: SPE1 layout (WELSPECS / COMPDAT)");
    if (!ok) return;
    try {
        auto refOf = [](const std::string& p) {
            Case c = importEclipseDeck(p, nullptr, nullptr);
            auto s = std::make_unique<Simulator>();
            s->setup(c);
            int w = -1;
            for (size_t i = 0; i < c.wells.size(); ++i) if (c.wells[i].name == "PROD") w = int(i);
            return s->wellRefDepth(w);
        };
        const double rOpen = refOf(c48Write("REF_OPEN.DATA", base)), rShut = refOf(c48Write("REF_SHUT.DATA", shut));
        std::printf("    reference depth: open K3 only %.2f, with SHUT K1 %.2f\n", rOpen, rShut);
        check(rOpen > 8300 && rShut > 8300, "C48: reference depths are cell-centre depths");
        check(rShut < rOpen - 10, "C48: a SHUT top connection sets the reference depth (shallower than the open K3)");
    } catch (const std::exception& e) { check(false, (std::string("C48: reference depth: ") + e.what()).c_str()); }
}

static void testC48GroupUda() {
    std::printf("== C48 group UDA: GCONPROD rate item given by a UDQ name follows the UDQ during the run\n");
    const std::string path = findExample("SPE1.DATA");
    if (path.empty()) { std::printf("    (examples/SPE1.DATA not found - skipped)\n"); return; }
    bool ok = true;
    const std::string anchor = "WCONPROD\n 'PROD' OPEN ORAT 20000 1* 1* 1* 1* 1000 /\n/\n";
    const std::string d = c48Replace(c32ReadText(path), anchor,
                                     anchor + "\nUDQ\nASSIGN FU_CAP 5000 /\n/\n\nGCONPROD\n 'FIELD' ORAT FU_CAP /\n/\n", &ok);
    check(ok, "C48: SPE1 layout (WCONPROD)");
    if (!ok) return;
    try {
        Case c = importEclipseDeck(c48Write("GUDA.DATA", d), nullptr, nullptr);
        bool parsed = false;
        for (const auto& st : c.schedule)
            for (const auto& g : st.gprod)
                for (const auto& u : g.uda) if (g.group == "FIELD" && u.first == "ORAT" && u.second == "FU_CAP") parsed = true;
        check(parsed, "C48: GCONPROD item 3 'FU_CAP' recorded as a group UDA");
        auto s = std::make_unique<Simulator>();
        s->setup(c); s->initialize();
        for (int k = 0; k < 2 && s->runReportStep(); ++k) {}
        const double t = s->groupRateTarget("FIELD", "ORAT");
        std::printf("    FIELD ORAT target after 2 report steps: %.6g\n", t);
        check(std::fabs(t - 5000) < 1e-6, "C48: the group target takes the UDQ value (5000)");
        double qo = 0;
        for (const auto& r : s->connectionRates()) qo += r.oil;
        std::printf("    field oil rate %.6g (well limit 20000)\n", qo);
        check(qo > 0 && qo < 5000 * 1.02, "C48: the field oil rate is held at the UDA group target");
    } catch (const std::exception& e) { check(false, (std::string("C48: group UDA: ") + e.what()).c_str()); }
}

static void testC48CellSnapshots() {
    std::printf("== C48 comp / thermal cell snapshots (shared by the GUI runner and ResSimCLI -o)\n");
    const std::string pc = findExample("CO2STORE.DATA");
    if (!pc.empty()) {
        try {
            const Case c = importEclipseDeck(pc, nullptr, nullptr);
            comp::DeckModel m; std::string err;
            if (!comp::deckModel(c, m, &err)) { check(false, ("C48: CO2STORE model: " + err).c_str()); }
            else {
                std::vector<comp::CompSummary> o; std::vector<Snapshot> snaps;
                const double pu = c.units == UnitSystem::Field ? 6894.757293168 : 1e5;
                comp::runDeckModel(m, o, [&](const comp::CompSummary& S, comp::CompSimulator& sim) {
                    snaps.push_back(comp::cellSnapshot(sim, S.time, pu)); return true; }, &err, 2);
                check(!snaps.empty(), "C48: CO2STORE snapshots written per report step");
                if (!snaps.empty()) {
                    const Snapshot& s = snaps.back();
                    const size_t n = s.pressure.size();
                    check(n > 0 && s.sg.size() == n && s.sw.size() == n && s.xco2.size() == n, "C48: CO2STORE PRESSURE / SWAT / SGAS / XCO2 per active cell");
                    double pmin = 1e30, pmax = -1e30;
                    for (float v : s.pressure) { pmin = std::min(pmin, double(v)); pmax = std::max(pmax, double(v)); }
                    check(pmin > 1 && pmax < 2000, "C48: CO2STORE pressure map in bar");
                    check(!s.field("SGTRAP").empty() || !s.field("ZMF1").empty(), "C48: CO2STORE named arrays (SGTRAP / ZMFk) present");
                }
            }
        } catch (const std::exception& e) { check(false, (std::string("C48: CO2STORE snapshots: ") + e.what()).c_str()); }
    }
    const std::string pt = findExample("STEAMFLOOD.DATA");
    if (!pt.empty()) {
        try {
            const Case c = importEclipseDeck(pt, nullptr, nullptr);
            thermal::DeckModel m; std::string err;
            if (!thermal::deckModel(c, m, &err)) { check(false, ("C48: STEAMFLOOD model: " + err).c_str()); }
            else {
                std::vector<thermal::Summary> o; std::vector<Snapshot> snaps; std::vector<float> t0;
                const size_t nAct = m.cellGlobal.size();
                thermal::runDeckModel(m, o, [&](const thermal::Summary& S, thermal::Simulator& sim) {
                    snaps.push_back(thermal::cellSnapshot(sim, c.thermal, S.time, nAct ? nAct : size_t(sim.size()), c.units == UnitSystem::Field, t0));
                    return true; }, &err, 2);
                check(!snaps.empty(), "C48: STEAMFLOOD snapshots written per report step");
                if (!snaps.empty()) {
                    const Snapshot& s = snaps.back();
                    check(s.pressure.size() == nAct && !s.field("TEMP").empty() && s.field("TEMP").size() == nAct,
                          "C48: STEAMFLOOD PRESSURE / TEMP per active cell");
                    check(!s.field("DTEMP").empty(), "C48: STEAMFLOOD temperature change map (DTEMP) present");
                }
            }
        } catch (const std::exception& e) { check(false, (std::string("C48: STEAMFLOOD snapshots: ") + e.what()).c_str()); }
    }
}
