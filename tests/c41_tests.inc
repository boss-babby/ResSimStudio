// C41 tests: workflow recorder format / runner, facies (truncated Gaussian) property models in Studies, fault network
#include "ressim/workflow.h"
#include "ressim/geostat.h"
#include "ressim/geomodel.h"

static void testC41Workflow() {
    using namespace ressim::wf;
    std::printf("== C41 workflow recorder: format, validation, runner\n");
    const std::vector<Step> w = {{Kind::Note, "demo"}, {Kind::Command, "Run simulation"}, {Kind::Wait, ""}, {Kind::ForVariants, ""},
                                 {Kind::Command, "B"}, {Kind::End, ""}, {Kind::Command, "Keyword: SWOF"}};
    std::vector<Step> r; std::string err;
    const bool rt = parse(serialize(w), r, err) && r.size() == w.size();
    bool same = rt;
    for (size_t i = 0; same && i < w.size(); ++i) same = r[i].kind == w[i].kind && r[i].arg == w[i].arg;
    std::vector<Step> x;
    const bool e1 = !parse("cmd\tA\n", x, err) && err.find("header") != std::string::npos;
    const bool e2 = !parse("RESSIM_WORKFLOW 1\nfor_variants\nfor_variants\nend\n", x, err) && err.find("nested") != std::string::npos;
    const bool e3 = !parse("RESSIM_WORKFLOW 1\nfor_variants\ncmd\tA\n", x, err) && err.find("not closed") != std::string::npos;
    const bool e4 = !parse("RESSIM_WORKFLOW 1\nfrobnicate\tA\n", x, err) && err.find("unknown step") != std::string::npos;
    const bool crlf = parse("\xEF\xBB\xBFRESSIM_WORKFLOW 1\r\n# comment\r\ncmd\tA\r\n\r\nwait\r\n", x, err) && x.size() == 2 && x[0].arg == "A";
    check(same && e1 && e2 && e3 && e4 && crlf, "C41: workflow file - round trip, BOM / CRLF / comments, header / nesting / unclosed loop / unknown step errors");

    // runner: a mock registry, a simulation that is busy for 4 ticks after "Run simulation", 3 variants
    std::vector<std::string> ran; int busyTicks = 0, loaded = -1, ticks = 0;
    std::vector<int> loads;
    auto exec = [&](const std::string& l) {
        if (l != "Run simulation" && l != "B" && l != "Keyword: SWOF") return false;
        ran.push_back(l + (loaded >= 0 ? "@" + std::to_string(loaded) : ""));
        if (l == "Run simulation") busyTicks = 4;
        return true;
    };
    auto busy = [&] { return busyTicks > 0; };
    Runner R;
    const bool started = R.start(w);
    while (R.active() && ticks < 200) {
        R.tick(exec, busy, [] { return 3; }, [&](int i) { loaded = i; loads.push_back(i); return true; });
        if (busyTicks > 0) --busyTicks;
        ++ticks;
    }
    const std::vector<std::string> want = {"Run simulation", "B@0", "B@1", "B@2", "Keyword: SWOF@2"};
    std::printf("    runner: %zu commands in %d ticks, variants loaded %zu, error '%s'\n", ran.size(), ticks, loads.size(), R.error.c_str());
    // the wait step must not finish while the run is busy: the command after it (B@0) comes >= 4 ticks after the run
    Runner R2; std::vector<std::string> ran2;
    R2.start({{Kind::Command, "B"}, {Kind::Command, "missing"}, {Kind::Command, "B"}});
    for (int t = 0; t < 20 && R2.active(); ++t) R2.tick([&](const std::string& l) { if (l != "B") return false; ran2.push_back(l); return true; }, [] { return false; }, nullptr, nullptr);
    Runner R3; std::vector<std::string> ran3;
    R3.start({{Kind::ForVariants, ""}, {Kind::Command, "B"}, {Kind::End, ""}, {Kind::Command, "B"}});
    for (int t = 0; t < 20 && R3.active(); ++t) R3.tick([&](const std::string& l) { ran3.push_back(l); return true; }, [] { return false; }, [] { return 0; }, [](int) { return true; });
    check(started && ran == want && R.error.empty() && loads == std::vector<int>({0, 1, 2}) && ticks >= 8,
          "C41: workflow runner - waits for the simulation, loops over every variant, runs the rest with the last variant loaded");
    check(ran2.size() == 1 && R2.error.find("missing") != std::string::npos && ran3.size() == 1,
          "C41: workflow runner - stops at an unknown command with an error; a variant loop without variants is skipped");
}

static void testC41Facies() {
    using namespace ressim;
    std::printf("== C41 facies property model (truncated Gaussian) and Studies re-simulation\n");
    Case base = makeDemoCase();
    Grid G; G.build(base);
    PropModelInput m;
    m.property = 3; m.model = 1; m.nugget = 0; m.sill = 1; m.rangeMax = m.rangeMin = 2000; m.rangeVert = 40; m.seed = 5;
    m.facProp = {0.5, 0.3, 0.2}; m.facPoro = {0.25, 0.15, 0.05}; m.facPerm = {500, 50, 1}; m.facSatnum = true; m.kvkh = 0.2; m.source = "test";
    const int qs[4] = {0, 9, 57, 140};
    for (int i = 0; i < 4; ++i) { const int q = qs[i]; if (q < int(G.center.size())) m.data.push_back({G.center[size_t(q)].x, G.center[size_t(q)].y, G.center[size_t(q)].z, double(1 + i % 3)}); }
    base.propModels = {m};
    Case a = base;
    const geostat::ApplyReport rep = geostat::applyFaciesModel(a, m);
    std::printf("    %s\n", rep.message.c_str());
    const size_t n = a.grid.satnum.size();
    bool consistent = rep.ok && n == size_t(a.grid.nx * a.grid.ny * a.grid.nz);
    for (size_t q = 0; consistent && q < n; ++q) {
        const int f = a.grid.satnum[q] - 1;
        consistent = f >= 0 && f < 3 && a.grid.poro[q] == m.facPoro[size_t(f)] && a.grid.permx[q] == m.facPerm[size_t(f)] && a.grid.permy[q] == a.grid.permx[q] &&
                     std::fabs(a.grid.permz[q] - 0.2 * a.grid.permx[q]) < 1e-12;
    }
    bool hard = true;
    for (int i = 0; i < 4; ++i) if (qs[i] < int(n)) hard = hard && a.grid.satnum[size_t(qs[i])] == 1 + i % 3;
    // proportions over several realisations (unconditional part dominates)
    double pr[3] = {0, 0, 0};
    const int nr = 8;
    for (int s = 0; s < nr; ++s) {
        Case c = base; PropModelInput mm = m; mm.seed = unsigned(100 + s); mm.data.clear();
        geostat::applyFaciesModel(c, mm);
        for (int f : c.grid.satnum) pr[f - 1] += 1.0 / double(c.grid.satnum.size() * nr);
    }
    std::printf("    mean proportions over %d realisations: %.3f %.3f %.3f (target 0.5 0.3 0.2)\n", nr, pr[0], pr[1], pr[2]);
    const bool prop = std::fabs(pr[0] - 0.5) < 0.1 && std::fabs(pr[1] - 0.3) < 0.1 && std::fabs(pr[2] - 0.2) < 0.1;
    check(consistent && hard && prop, "C41: truncated Gaussian facies - per-facies PORO / PERM / SATNUM, hard data honoured, proportions near the targets");
    // Studies: GeoSeed re-simulates the facies model; persistence through the project text
    StudyParam ps; ps.kind = ParamKind::GeoSeed; ps.target = -1; ps.base = 5; ps.low = 1; ps.high = 100; ps.log = false;
    const Case s7 = applyParams(base, {ps}, {7}), s7b = applyParams(base, {ps}, {7}), s8 = applyParams(base, {ps}, {8});
    const bool studies = s7.grid.satnum == s7b.grid.satnum && s7.grid.satnum != s8.grid.satnum && s7.grid.permx != s8.grid.permx;
    const Case rt = projectFromText(projectText(base));
    const bool persist = rt.propModels.size() == 1 && rt.propModels[0].property == 3 && rt.propModels[0].facSatnum && rt.propModels[0].facProp.size() == 3 &&
                         std::fabs(rt.propModels[0].facPerm[1] - 50) < 1e-9 && std::fabs(rt.propModels[0].facPoro[2] - 0.05) < 1e-12 && rt.propModels[0].data.size() == m.data.size();
    const bool label = paramLabel(ps, base).find("seed") != std::string::npos;
    Case bad = base; PropModelInput one = m; one.facProp = {1.0};
    const bool rejects = !geostat::applyFaciesModel(bad, one).ok;
    check(studies && persist && label && rejects, "C41: facies model in Studies - GeoSeed reproducible / different per seed, PROPFAC persisted, < 2 facies rejected");
}

static void testC41FaultNetwork() {
    using namespace ressim::geo;
    std::printf("== C41 fault network: junctions, truncation, closed fault blocks, throw profile\n");
    auto vstick = [](double x, double y, double dxdz = 0) { return std::vector<P3>{{x, y, 900.0}, {x + dxdz * 400, y, 1300.0}}; };
    Fault F1; F1.name = "F1"; for (double y = 0; y <= 2000; y += 500) F1.sticks.push_back(vstick(1000, y));
    Fault F2; F2.name = "F2"; for (double x : {1150.0, 1400.0, 1700.0, 2000.0}) F2.sticks.push_back(vstick(x, 1050));        // ends 150 short of F1
    Fault F3; F3.name = "F3"; for (double x : {700.0, 850.0, 1200.0, 1600.0, 2000.0}) F3.sticks.push_back(vstick(x, 450));   // crosses F1
    const double tol = 250;
    const std::vector<Fault> fl = {F1, F2, F3};
    const auto J = detectJunctions(fl, 1050, tol);
    const FaultJunction* jy = nullptr; const FaultJunction* jx = nullptr;
    for (const auto& j : J) { if (j.minor == 1 && j.major == 0) jy = &j; if (j.type == 'X' && ((j.minor == 2 && j.major == 0) || (j.minor == 0 && j.major == 2))) jx = &j; }
    std::printf("    %zu junctions; Y gap %.1f at (%.0f, %.0f); X at (%.0f, %.0f)\n", J.size(), jy ? jy->gap : -1, jy ? jy->x : 0, jy ? jy->y : 0, jx ? jx->x : 0, jx ? jx->y : 0);
    const bool det = jy && jy->type == 'Y' && std::fabs(jy->gap - 150) < 1e-6 && std::fabs(jy->x - 1000) < 1e-6 && jx && !jx->truncate && std::fabs(jx->x - 1000) < 1e-6 && std::fabs(jx->y - 450) < 1e-6;
    check(det, "C41: fault junctions - Y junction with its gap (minor tip short of the major), X crossing found and not truncated by default");
    // truncation: F2 snapped onto F1; F3 cut at F1 (east side kept, the side with most sticks)
    std::vector<FaultJunction> rules = J;
    for (auto& j : rules) if (j.type == 'X') { j.truncate = true; j.major = 0; j.minor = 2; }
    FaultNetReport rep;
    const auto T = truncateFaults(fl, rules, tol, &rep);
    double minX2 = 1e300, minX3 = 1e300;
    for (double z : {900.0, 1100.0, 1300.0}) {
        for (const auto& p : T[1].traceAt(z)) minX2 = std::min(minX2, p[0]);
        for (const auto& p : T[2].traceAt(z)) minX3 = std::min(minX3, p[0]);
    }
    std::printf("    truncated %d, sticks removed %d; F2 west end %.6f, F3 west end %.6f\n", rep.truncated, rep.sticksRemoved, minX2, minX3);
    const bool trunc = rep.truncated == 2 && rep.sticksRemoved == 2 && std::fabs(minX2 - 1000) < 1e-6 && std::fabs(minX3 - 1000) < 1e-6 && T[0].sticks.size() == F1.sticks.size();
    // closed fault blocks: without the network F2 leaves the NE block open (2 blocks); with it, 3 (F3 not used here)
    Surface h1, h2;
    for (Surface* h : {&h1, &h2}) { h->nx = h->ny = 21; h->x0 = h->y0 = 0; h->dx = h->dy = 100; h->z.assign(441, h == &h1 ? 1000.0 : 1100.0); }
    h1.name = "TOP"; h2.name = "BASE";
    GridSpec gs; gs.nx = gs.ny = 20; gs.dx = gs.dy = 100;
    std::string e1, e2;
    const GeoGrid g0 = buildGrid(gs, {h1, h2}, {F1, F2}, {}, e1);
    const GeoGrid g1 = buildGrid(gs, {h1, h2}, {T[0], T[1]}, {}, e2);
    std::printf("    fault blocks: %d without the network, %d with it\n", g0.nBlocks, g1.nBlocks);
    check(trunc && e1.empty() && e2.empty() && g0.nBlocks == 2 && g1.nBlocks == 3,
          "C41: fault truncation - minor faults end exactly on the major at every depth, the Y-junction gap closes a fault block");
    // throw profile: F dips east (hanging wall east), the horizon is 20 deeper east of the fault
    Fault Fd; Fd.name = "Fd"; for (double y = 0; y <= 2000; y += 500) Fd.sticks.push_back(vstick(1000, y, 0.1));
    Surface hs = h1; for (int j = 0; j < 21; ++j) for (int i = 0; i < 21; ++i) hs.z[size_t(i + 21 * j)] = i * 100.0 > 1010 ? 1020.0 : 1000.0;
    const ThrowProfile P = faultThrowProfile(Fd, hs, 150, 11);
    bool all20 = !P.throwV.empty();
    for (double t : P.throwV) all20 = all20 && std::fabs(t - 20) < 1e-9;
    Surface hr = hs; for (double& z : hr.z) z = 2020 - z;   // reverse throw
    const ThrowProfile Pr = faultThrowProfile(Fd, hr, 150, 11);
    std::printf("    throw profile: %zu samples, max %.3f, reverses %d, tips taper %d; reversed horizon max %.3f first %.3f\n", P.throwV.size(), P.maxThrow, int(P.reverses), int(P.tipsTaper),
                Pr.maxThrow, Pr.throwV.empty() ? 0.0 : Pr.throwV[0]);
    check(all20 && P.throwV.size() == 11 && !P.reverses && !P.tipsTaper && !Pr.throwV.empty() && std::fabs(Pr.throwV[0] + 20) < 1e-9,
          "C41: horizon-fault throw profile - normal throw positive on the hanging wall (dip side), constant throw flagged as not tapering, reverse throw negative");
}
