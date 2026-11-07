// C46 tests: geomodel project file (.rgeo) round trip and object-tree folders
#include "ressim/geo_project.h"
#include <cmath>
#include <cstdio>

static void testC46GeoProject() {
    using namespace ressim::geo;
    std::printf("== C46 geomodel project file (.rgeo): exact round trip of inputs, interpretation, folders and settings\n");
    Project P; P.field = true;
    Surface s; s.name = "Top Reservoir"; s.nx = 3; s.ny = 2; s.x0 = 451000.125; s.y0 = 6780000.5; s.dx = 25; s.dy = 37.5; s.rotDeg = 12.25;
    s.z = {1500.1, 1501.2, std::nan(""), 1503.3, 1.0 / 3.0, 1505.5};
    P.surfaces.push_back(s);
    GeoWell w; w.name = "A-1 H"; w.x = 451100; w.y = 6780050; w.kb = 24.5; w.trajSrc = 1;
    w.traj = {{0, 0, 0, 451100, 6780050, -24.5}, {1200, 45.5, 130.25, 451300.5, 6779900.25, 1100.75}};
    WellLog lg; lg.well = w.name; lg.file = "C:/data/a 1.las"; lg.nullValue = -999.25; lg.md = {1000, 1000.5, 1001};
    LogCurve c; c.mnem = "PHIE"; c.unit = "v/v"; c.descr = "effective porosity"; c.v = {0.21, -999.25, 0.1 + 0.2}; lg.curves.push_back(c);
    LogCurve f; f.mnem = "FACIES"; f.discrete = true; f.v = {1, 2, 2}; lg.curves.push_back(f);
    w.logs.push_back(lg); w.tops = {{"Top Reservoir", 1010.25}, {"Base", 1100}};
    P.wells.push_back(w);
    Polygon pg; pg.name = "lease outline"; pg.closed = true; pg.pts = {{451000, 6780000}, {452000, 6780000}, {452000, 6781000}};
    P.polygons.push_back(pg);
    Fault fl; fl.name = "F1"; fl.throwHint = 15; fl.mult = 0.01; fl.sticks = {{{451500, 6780000, 1400}, {451520, 6780010, 1700}}, {{451600, 6780500, 1390}}};
    P.faults.push_back(fl);
    SeismicLine ln; ln.name = "IL 1001"; ln.path = "/data/line 1001.sgy"; ln.ns = 4; ln.format = 1; ln.t0 = 0; ln.dt = 4;
    ln.cdp = {10, 11}; ln.x = {451000, 451025}; ln.y = {6780000, 6780000}; ln.a.assign(8, 1.0f);
    P.lines.push_back(ln);
    LineHorizon lh; lh.name = "H1"; lh.t = {{1200.5, std::nan("")}}; P.lineHorizons.push_back(lh);
    P.cubePath = "/data/full stack.segy";
    P.folders.addFolder(ObjKind::Well, "Producers 2024");
    P.folders.assign(ObjKind::Well, "A-1 H", "Producers 2024");
    P.folders.assign(ObjKind::Surface, "Top Reservoir", "Interpreted");
    P.settings = {{"ws.order", "A-1 H|B-2"}, {"template", "Porosity|PHI|v/v|1|0|1|0|0|0.35"}, {"empty", ""}};

    const std::string t1 = projectToText(P);
    Project Q; std::string err;
    const bool ok = projectFromText(t1, Q, err);
    check(ok, ("project text parses back" + (err.empty() ? std::string() : ": " + err)).c_str());
    const std::string t2 = projectToText(Q);
    check(t1 == t2, "write -> read -> write is byte-identical");
    bool same = Q.field && Q.surfaces.size() == 1 && Q.wells.size() == 1 && Q.polygons.size() == 1 && Q.faults.size() == 1 && Q.lines.size() == 1;
    if (same) {
        const Surface& a = Q.surfaces[0];
        same = a.name == s.name && a.nx == 3 && a.ny == 2 && a.x0 == s.x0 && a.rotDeg == s.rotDeg && std::isnan(a.z[2]) && a.z[4] == s.z[4];
        const GeoWell& b = Q.wells[0];
        same = same && b.name == w.name && b.kb == w.kb && b.trajSrc == 1 && b.traj.size() == 2 && b.traj[1].azi == 130.25 && b.traj[1].tvd == 1100.75 &&
               b.logs.size() == 1 && b.logs[0].file == lg.file && b.logs[0].curves.size() == 2 && b.logs[0].curves[0].v[2] == 0.1 + 0.2 &&
               b.logs[0].curves[0].descr == c.descr && b.logs[0].curves[1].discrete && b.tops.size() == 2 && b.tops[0].surface == "Top Reservoir";
        same = same && Q.polygons[0].name == pg.name && Q.polygons[0].closed && Q.polygons[0].pts.size() == 3;
        same = same && Q.faults[0].sticks.size() == 2 && Q.faults[0].sticks[0][1][2] == 1700 && Q.faults[0].mult == 0.01 && Q.faults[0].throwHint == 15;
        same = same && Q.lines[0].path == ln.path && Q.lines[0].cdp == ln.cdp && Q.lines[0].x == ln.x && Q.lines[0].ns == 4 && Q.lines[0].a.empty();
        same = same && Q.lineHorizons.size() == 1 && Q.lineHorizons[0].t[0][0] == 1200.5 && std::isnan(Q.lineHorizons[0].t[0][1]);
        same = same && Q.cubePath == P.cubePath && Q.settings.size() == 3 && Q.settings[1].second == P.settings[1].second && Q.settings[2].second.empty();
    }
    check(same, "names with spaces, exact doubles (%.17g), NaN nodes / picks, logs, tops, surveys, sticks, line geometry (no amplitudes), settings");
    check(Q.folders.folderOf(ObjKind::Well, "A-1 H") == "Producers 2024" && Q.folders.folderOf(ObjKind::Surface, "Top Reservoir") == "Interpreted" &&
              Q.folders.foldersOf(ObjKind::Surface).size() == 1 && Q.folders.folderOf(ObjKind::Polygon, "lease outline").empty(),
          "folder membership survives (assign creates the folder; other kinds stay at the top level)");
    // folder operations
    ObjFolders F = Q.folders;
    F.renameObject(ObjKind::Well, "A-1 H", "A-1 HST");
    F.renameFolder(ObjKind::Well, "Producers 2024", "Producers");
    const bool r1 = F.folderOf(ObjKind::Well, "A-1 HST") == "Producers" && F.folderOf(ObjKind::Well, "A-1 H").empty();
    F.addFolder(ObjKind::Well, "Injectors");
    const bool dupe = !F.addFolder(ObjKind::Well, "Injectors") && !F.addFolder(ObjKind::Well, "");
    F.assign(ObjKind::Well, "A-1 HST", "Injectors");
    const bool r2 = F.folderOf(ObjKind::Well, "A-1 HST") == "Injectors";
    F.removeFolder(ObjKind::Well, "Injectors");
    const bool r3 = F.folderOf(ObjKind::Well, "A-1 HST").empty() && F.foldersOf(ObjKind::Well).size() == 1;
    F.assign(ObjKind::Fault, "F1", "North"); F.forget(ObjKind::Fault, "F1");
    const bool r4 = F.folderOf(ObjKind::Fault, "F1").empty() && F.foldersOf(ObjKind::Fault).size() == 1;   // folder stays, member gone
    check(r1 && dupe && r2 && r3 && r4, "folders: rename object / folder, duplicate refused, move, remove folder (members to the top level), forget");
    // file save / load (atomic tmp + rename) and malformed input
    const std::string path = "c46_project_test.rgeo";
    std::string e2; Project R;
    const bool io = saveProject(path, P, &e2) && loadProject(path, R, e2) && projectToText(R) == t1;
    std::remove(path.c_str());
    check(io, "saveProject / loadProject file round trip");
    Project bad; std::string e3;
    const bool rej = !projectFromText("NOT A PROJECT\n", bad, e3) && !e3.empty() &&
                     !projectFromText("RESGEO-PROJECT 1\nSURFACE x 2 2 0 0 1 1 0\n1 2\nEND\n", bad, e3);
    check(rej, ("wrong header / truncated surface rejected with a message (" + e3 + ")").c_str());
}

// C46 completion numbers (COMPLUMP / definition order) in WELOPEN and WPIMULT, WELPI on connections defined SHUT, the
// surface-volume PI per phase (WPIO / WPIW / WPIG / WPIL / WPI summary vectors; OPM: gas incl. Rs x oil, oil incl. Rv x gas)
static void testC46CompletionNumbers() {
    std::printf("== C46 completion numbers (WELOPEN / WPIMULT items 6-7, COMPLUMP), WELPI on SHUT-defined connections, WPI vectors\n");
    std::string d = eclipseDeckToString(smallCase());
    const size_t p0 = d.find("WCONPROD"), p1 = d.find("END", p0);
    if (p0 == std::string::npos || p1 == std::string::npos) { check(false, "C46: deck layout"); return; }
    d = d.substr(0, p0) +
        "COMPDAT\n 'PROD' 9 10 1 1 SHUT 1* 1* 0.5 1* 0 1* Z /\n/\n\nCOMPLUMP\n 'PROD' 9 10 1 1 7 /\n/\n\n"
        "WPIMULT\n 'PROD' 0.5 1* 1* 1* 7 7 /\n/\n\n"
        "WCONPROD\n 'PROD' OPEN ORAT 6000 /\n/\n\nWCONINJE\n 'INJ' WATER OPEN RATE 8000 /\n/\n\nTSTEP\n 91.3125 /\n\n"
        "WELPI\n 'PROD' 30 /\n/\n\nTSTEP\n 91.3125 /\n\n"
        "WELOPEN\n 'PROD' OPEN 1* 1* 1* 7 7 /\n/\n\nTSTEP\n 91.3125 /\n\n"
        "WELOPEN\n 'PROD' SHUT 1* 1* 1* 2 2 /\n/\n\nTSTEP\n 91.3125 /\n\nEND\n";
    auto rec = [](const Case& c, int i, int j, int k) -> const Completion* {
        const Completion* best = nullptr;
        for (const auto& w : c.wells) if (w.name == "PROD")
            for (const auto& q : w.completions) if (q.i == i && q.j == j && k >= q.k1 && k <= q.k2 && (!best || q.fromStep > best->fromStep)) best = &q;
        return best;
    };
    try {
        const Case a = eclipseDeckFromString(d, "");
        const Completion* n3 = rec(a, 8, 9, 0);
        const Completion* l2 = rec(a, 9, 9, 1);
        const Completion* l1 = rec(a, 9, 9, 0);
        std::printf("    lumped connection 7: open %d from step %d (defined at %d), WPIMULT %.3g; layer 2 shut at step %d; layer 1 until %d\n",
                    n3 ? int(n3->open) : -1, n3 ? n3->fromStep : -9, n3 ? n3->defStep : -9, n3 ? n3->mult : -1, l2 ? l2->untilStep : -9, l1 ? l1->untilStep : -9);
        check(n3 && n3->open && n3->fromStep == 2 && n3->defStep == 0 && std::fabs(n3->mult - 0.5) < 1e-12,
              "C46: COMPLUMP number 7 - WPIMULT scales the SHUT connection, WELOPEN items 6-7 open it at step 2 (definition step kept)");
        check(l2 && l2->k1 == 1 && l2->k2 == 1 && l2->untilStep == 3 && l1 && l1->k1 == 0 && l1->k2 == 0 && l1->untilStep < 0,
              "C46: default completion numbers = definition order per layer - WELOPEN SHUT 2 2 splits the record and shuts layer 2 only");
        const Case b = eclipseDeckFromString(eclipseDeckToString(a), "");
        const Completion* m3 = rec(b, 8, 9, 0);
        check(m3 && m3->open && m3->fromStep == 2 && std::fabs(m3->mult - 0.5) < 1e-12, "C46: deck round trip keeps the opening step and the WPIMULT factor");
        // run: WELPI at step 1 -> WPIO = 30 at the start of step 1; opening the lumped connection raises the PI; WPIW / WPIG reported
        double logPi = -1, logTarget = -1;
        Simulator s;
        s.log = [&](const std::string& m) { const size_t q = m.find("WELPI PROD: PI "); if (q != std::string::npos) std::sscanf(m.c_str() + q + 15, "%lf -> %lf", &logPi, &logTarget); };
        s.setup(a); s.initialize();
        while (s.runReportStep()) {}
        int pw = -1; for (size_t w = 0; w < a.wells.size(); ++w) if (a.wells[w].name == "PROD") pw = int(w);
        double piEnd0 = -1, piStart1 = -1, piEnd1 = -1, piStart2 = -1, piEnd2 = -1, piStart3 = -1, pref = -1, pg = -1;
        for (const SummaryPoint& p : s.summary()) {
            const double v = p.wells[size_t(pw)].pi[0];
            if (p.time <= 91.3125 + 1e-9) piEnd0 = v;
            if (p.time > 91.3125 + 1e-9 && p.time <= 182.625 + 1e-9) { if (piStart1 < 0) { piStart1 = v; pref = p.wells[size_t(pw)].piPref; pg = p.wells[size_t(pw)].pi[2]; } piEnd1 = v; }
            if (p.time > 182.625 + 1e-9 && p.time <= 273.9375 + 1e-9) { if (piStart2 < 0) piStart2 = v; piEnd2 = v; }
            if (p.time > 273.9375 + 1e-9 && piStart3 < 0) piStart3 = v;
        }
        std::printf("    WPIO end of step 0 %.6g, WELPI log PI %.6g -> %.4g; step 1 %.4g .. %.4g (WPI %.4g, WPIG %.4g), step 2 %.4g .. %.4g (connection 7 opened), step 3 %.4g (layer 2 shut)\n",
                    piEnd0, logPi, logTarget, piStart1, piEnd1, pref, pg, piStart2, piEnd2, piStart3);
        check(std::fabs(logPi - piEnd0) < 1e-6 * piEnd0 && std::fabs(logTarget - 30) < 1e-12 && piStart1 > 0.8 * 30 && piStart1 <= 30 * 1.001 && std::fabs(pref - piStart1) < 1e-12,
              "C46: WELPI uses the reported oil PI (WPIO = WPI for an oil producer) and rescales it to the target");
        check(piStart2 > 1.2 * piEnd1 && piStart3 < piEnd2, "C46: the PI follows the connections opened / shut by completion number");
    } catch (const std::exception& e) { std::printf("    %s\n", e.what()); check(false, "C46: completion-number deck parses"); }
}

static void testC46SchedMultflt() {
    std::printf("== C46 MULTFLT in SCHEDULE: cumulative fault multipliers applied from their report step (OPM semantics)\n");
    Case c0 = smallCase();
    c0.schedule.resize(4);
    FaultInput F; F.name = "F1"; FaultFace ff; ff.i1 = ff.i2 = 4; ff.j1 = 0; ff.j2 = 9; ff.k1 = 0; ff.k2 = 1; ff.face = 'X'; F.faces.push_back(ff);
    c0.grid.faults.push_back(F);
    std::string d = eclipseDeckToString(c0);
    // the exporter writes the 4 equal steps as one TSTEP 4*dt record -> split it: 2 steps, MULTFLT, 2 steps (from step 2)
    const size_t p = d.find("\nTSTEP", d.find("\nSCHEDULE")), q = p == std::string::npos ? p : d.find("/\n", p);
    if (q == std::string::npos) { check(false, "C46: MULTFLT deck layout"); return; }
    d.replace(p + 1, q + 2 - (p + 1), "TSTEP\n 2*91.3125 /\n\nMULTFLT\n 'F1' 0.001 /\n/\n\nMULTFLT\n 'F*' 0.5 /\n/\n\nTSTEP\n 2*91.3125 /\n");
    try {
        const Case a = eclipseDeckFromString(d, "");
        double gm = -1; for (const auto& f : a.grid.faults) if (f.name == "F1") gm = f.mult;
        std::printf("    %zu SCHEDULE MULTFLT events (first at step %d, x%g / x%g); GRID multiplier %g\n", a.schedMultflt.size(),
                    a.schedMultflt.empty() ? -1 : a.schedMultflt[0].step, a.schedMultflt.empty() ? 0 : a.schedMultflt[0].mult,
                    a.schedMultflt.size() > 1 ? a.schedMultflt[1].mult : 0, gm);
        check(a.schedMultflt.size() == 2 && a.schedMultflt[0].step == 2 && a.schedMultflt[1].step == 2 && std::fabs(gm - 1) < 1e-12,
              "C46: SCHEDULE MULTFLT parsed as time-dependent events (not folded into the GRID multiplier)");
        const Case b = eclipseDeckFromString(eclipseDeckToString(a), "");
        check(b.schedMultflt.size() == 2 && b.schedMultflt[0].step == 2 && std::fabs(b.schedMultflt[0].mult - 0.001) < 1e-15,
              "C46: deck round trip keeps the SCHEDULE MULTFLT at its report step");
        // run: no change before step 2, then x0.001 x 0.5 on the 20 connections across the fault (10 rows x 2 layers)
        auto run = [](const Case& cc, std::vector<std::string>& logs) {
            Simulator s; s.log = [&](const std::string& m) { if (m.find("MULTFLT (SCHEDULE)") != std::string::npos) logs.push_back(m); };
            s.setup(cc); s.initialize();
            std::vector<double> bhp;
            int pi = -1; for (size_t w = 0; w < cc.wells.size(); ++w) if (cc.wells[w].type != WellType::Producer) pi = int(w);
            while (s.runReportStep()) bhp.push_back(s.summary().back().wells[size_t(pi)].bhp);
            return bhp;
        };
        std::vector<std::string> la, l0;
        const std::vector<double> ba = run(a, la);
        Case n = a; n.schedMultflt.clear();
        const std::vector<double> b0 = run(n, l0);
        for (const auto& m : la) std::printf("    %s\n", m.c_str());
        std::printf("    injector BHP end of steps 1 / 2 / 3: %.6g / %.6g / %.6g (no MULTFLT: %.6g / %.6g / %.6g)\n",
                    ba.size() > 2 ? ba[1] : 0, ba.size() > 2 ? ba[2] : 0, ba.size() > 3 ? ba[3] : 0, b0.size() > 2 ? b0[1] : 0, b0.size() > 2 ? b0[2] : 0, b0.size() > 3 ? b0[3] : 0);
        check(la.size() == 2 && la[0].find("on 20 connections") != std::string::npos, "C46: MULTFLT events act on the 20 connections across the fault");
        check(ba.size() >= 4 && b0.size() >= 4 && std::fabs(ba[1] - b0[1]) < 1e-9 * std::fabs(b0[1]) && ba[3] > b0[3] + 1e-3,
              "C46: identical before step 2, the sealed fault raises the injector pressure afterwards");
    } catch (const std::exception& e) { std::printf("    %s\n", e.what()); check(false, "C46: SCHEDULE MULTFLT deck parses"); }
}
