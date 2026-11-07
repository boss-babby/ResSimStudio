// C47 tests: OPM low-pressure padding of gas PVT tables, CO2STORE import messages, per-face fault history CSV
#include "ressim/fluid.h"
#include "ressim/faultmech.h"
#include <cmath>
#include <cstdio>
#include <fstream>

static void testC47GasPadding() {
    std::printf("== C47 OPM low-pressure padding of PVDG / PVTG (opm-common PR 3779): limiting row + 1 bar row\n");
    Case c = makeDemoCase();
    c.units = UnitSystem::Metric;
    c.phases.vapoil = false;
    // model2 PVDG (PVT-DRY-GAS-LIVE-OIL.INC): b = 1/Bg extrapolated from 50 / 70 bar is negative at 1 bar -> padding
    c.props.pvdg = {{50, 0.024958, 0.01440}, {70, 0.017643, 0.01490}, {90, 0.013615, 0.01544}, {110, 0.011081, 0.01605}};
    FluidModel f; f.setup(c);
    check(f.gasPaddingRows() == 2, "C47: model2 PVDG gets two padding rows (limiting pressure + 1 bar)");
    // OPM: B_lim = max(10 Bg0, 1) = 1; p_lim where b = 1 on the line through rows 0 / 1; 1 bar row with 1.1 B_lim
    const double b0 = 1 / 0.024958, b1 = 1 / 0.017643, pl = 50 + (1.0 - b0) / (b1 - b0) * 20;
    check(std::fabs(f.bgS(pl) - 1.0) < 1e-9, "C47: b = 1 / B_lim = 1 at the limiting pressure");
    check(std::fabs(f.bgS(1.0) - 1 / 1.1) < 1e-9, "C47: b = 1 / (1.1 B_lim) at 1 bar");
    check(f.bgS(1.01325) > 0 && f.bgS(0.5) > 0, "C47: positive gas b at 1 atm and below (no negative densities in potentials)");
    check(std::fabs(f.bgS(80) - 0.5 * (1 / 0.017643 + 1 / 0.013615)) < 1e-9, "C47: table interior unchanged");
    // a table that already reaches low pressure is not padded
    Case d = c; d.props.pvdg = {{1, 0.9, 0.012}, {50, 0.024958, 0.01440}, {70, 0.017643, 0.01490}};
    FluidModel g; g.setup(d);
    check(g.gasPaddingRows() == 0, "C47: no padding when the table starts at 1 bar");
    // zero-pressure branch: b(1 bar) above 1/B_lim but b(0) < 0 -> one row (0, 1.1 Bg0, mu0)
    Case e = c; e.props.pvdg = {{20, 0.06, 0.013}, {40, 0.028, 0.0137}};   // b 16.7 -> 35.7: b(1) = -1.4 < 0.1 -> first branch
    FluidModel h; h.setup(e);
    check(h.gasPaddingRows() >= 1, "C47: steep low-pressure table padded");
}

static void testC47Co2storeImport() {
    std::printf("== C47 CO2STORE decks: no false 'not modelled' warning, validation passes (ResSim CO2STORE core)\n");
    std::vector<std::string> w;
    try {
        const Case c = importEclipseDeck(findDeck("examples/CO2STORE.DATA"), &w);
        bool falseWarn = false;
        for (const auto& s : w) if (s.find("CO2 storage (CO2STORE) is not modelled") != std::string::npos || s.find("cannot model this deck") != std::string::npos) falseWarn = true;
        check(!falseWarn, "C47: CO2STORE import carries no 'not modelled' / 'cannot model' warning");
        check(validateCase(c).empty(), "C47: CO2STORE case validates");
        check(c.comp.enabled && c.comp.co2store, "C47: CO2STORE runs on the compositional CO2STORE core");
    } catch (const std::exception& ex) { check(false, (std::string("C47: CO2STORE import: ") + ex.what()).c_str()); }
}

static void testC47FaultCsv() {
    std::printf("== C47 per-face fault history CSV (--fault-csv / RESSIM_FAULT_CSV)\n");
    faultmech::Output o;
    o.names = {"F1"};
    faultmech::Face F; F.fault = 0; F.k = 2; F.area = 500; F.aSig = 1; F.c = {150.f, 50.f, 1012.5f};
    o.faces = {F, F};
    for (int s = 0; s < 3; ++s) {
        faultmech::StepOut S; S.t = 10.0 * s;
        S.cff = {-36.f + s, -35.f + s}; S.dcfs = {float(s), float(0.5 * s)}; S.rate = {1.f, 1.f}; S.cum = {0.01f * s, 0.02f * s};
        o.steps.push_back(S);
    }
    const std::string path = "c47_fault.csv";
    check(faultmech::writeFaultCsv(o, path), "C47: fault CSV written");
    std::ifstream in(path); std::string line; int n = 0; std::string header, last;
    while (std::getline(in, line)) { if (n == 0) header = line; last = line; ++n; }
    check(header == "t,face,fault,k,x,y,depth,area,aSig,cff,dcfs,rate,cum,sne,mult", "C47: fault CSV header (C49: + sne, mult)");
    check(n == 1 + 3 * 2, "C47: one row per report step and face");
    check(last.rfind("20,1,F1,3,", 0) == 0, "C47: last row = last step, face 1, fault name, 1-based layer");
    std::remove(path.c_str());
}

static void testC47ActionWellsJly() {
    std::printf("== C47 ACTIONX: month JLY; a well whose connections all arrive by ACTIONX gets its defaulted reference depth\n");
    {
        using namespace ressim::udq;
        Source s; s.field = [](const std::string&) { return (double)NAN; };
        Set set; Action a; std::string m;
        const bool p = parseAction("J /\nMNTH >= JLY AND YEAR >= 2020 /\n/\nWELOPEN\n'W' OPEN /\n/\n", 0, a, m);
        s.day = 1; s.year = 2020;
        s.month = 6; const bool jun = evaluateAction(a, s, set).fired;
        s.month = 7; const bool jul = evaluateAction(a, s, set).fired;
        s.month = 8; const bool aug = evaluateAction(a, s, set).fired;
        check(p && !jun && jul && aug, "C47: 'MNTH >= JLY' = July (ECLIPSE spelling; opm-tests ACTIONX_WCONINJE ACT-04)");
    }
    const std::string path = findExample("ACTIONX_UDQ.DATA");
    if (path.empty()) { std::printf("    (examples/ACTIONX_UDQ.DATA not found - skipped)\n"); return; }
    namespace fs = std::filesystem;
    std::string d = c32ReadText(path);
    auto replaceOnce = [](std::string s, const std::string& a, const std::string& b) { const size_t q = s.find(a); if (q != std::string::npos) s.replace(q, a.size(), b); return s; };
    const std::string pAll = " 'P*' 2* 1 3 OPEN 1* 1* 0.5 /\n";
    if (d.find(pAll) == std::string::npos) { check(false, "C47: ACTIONX_UDQ layout (P* COMPDAT)"); return; }
    // P2 has no connection at all until the action adds layers 2-3
    d = replaceOnce(d, pAll, " 'P1' 2* 1 3 OPEN 1* 1* 0.5 /\n 'P3' 2* 1 3 OPEN 1* 1* 0.5 /\n 'P4' 2* 1 3 OPEN 1* 1* 0.5 /\n");
    const std::string action = "ACTIONX\n NEWP2 1 /\n FOPT > 1E6 /\n/\nCOMPDAT\n 'P2' 2* 2 3 OPEN 1* 1* 0.5 /\n/\nWELOPEN\n 'P2' OPEN /\n/\nENDACTIO\n\n";
    const std::string tsBlock = "TSTEP\n 10 20 30 30 24*60\n/\n";
    if (d.find(tsBlock) == std::string::npos) { check(false, "C47: ACTIONX_UDQ layout (TSTEP)"); return; }
    d = replaceOnce(d, tsBlock, action + tsBlock);
    const std::string dir = (fs::temp_directory_path() / "ressim_c47_actx").string();
    fs::create_directories(dir);
    const std::string p = dir + "/NEWP2.DATA";
    std::ofstream(p, std::ios::binary) << d;
    try {
        Case c = importEclipseDeck(p, nullptr, nullptr);
        auto s = std::make_unique<Simulator>();
        std::vector<std::string> logs;
        s->log = [&logs](const std::string& l) { logs.push_back(l); };
        s->setup(c); s->initialize();
        // check right after the action (fires at the end of report step 7, 270 d): later actions / limits may close P2
        std::vector<Simulator::ConnRate> cr;
        for (int st = 0; s->runReportStep(); ++st) if (st == 8) cr = s->connectionRates();
        if (std::getenv("C47_VERBOSE")) for (const auto& l : logs) if (l.find("ACTION") != std::string::npos || l.find("P2") != std::string::npos) std::printf("    | %s\n", l.c_str());
        double qo = 0, worstHead = 0; int n = 0;
        if (std::getenv("C47_VERBOSE")) for (const auto& r : cr) std::printf("    conn %s k %d qo %.4g pb %.5g pc %.5g\n", r.well.c_str(), r.k, r.oil, r.pBlock, r.pConn);
        for (const auto& r : cr)
            if (r.well == "P2") { ++n; qo += r.oil; worstHead = std::max(worstHead, std::fabs(r.pConn - r.pBlock) / std::max(1.0, r.pBlock)); }
        check(n == 2, "C47: the action added both P2 connections");
        check(qo > 0, "C47: P2 produces through the run-time connections");
        check(worstHead < 0.3, "C47: P2 connection pressures near the block pressures (reference depth = top connection, not 0)");
    } catch (const std::exception& e) { check(false, (std::string("C47: NEWP2 run: ") + e.what()).c_str()); }
    std::error_code ec; fs::remove_all(dir, ec);
}
