// =====================================================================================================================
// C25 tests: isoparametric FE geomechanics on corner-point grids, sequential geomechanics in black oil / compositional,
// FE closure stress for HF growth, fault mechanics with FE stress in all three engines, thermo-elastic fault stress
// =====================================================================================================================
#include "ressim/geoseq.h"
#include "ressim/comp_deck.h"
#include "ressim/therm_deck.h"

static geomech::Input c25Input() {
    geomech::Input in; in.enabled = true; in.coupling = 1;
    in.E = {10e9}; in.nu = {0.25}; in.biot = {1.0}; in.alphaT = {1e-5}; in.Ks = 0;
    in.sigmaVGrad = 23e3; in.k0 = 0.7;
    return in;
}

static void testC25HexPatch() {
    std::printf("== C25 isoparametric hexahedra: brick equivalence and patch test on a distorted mesh\n");
    const int nx = 4, ny = 3, nz = 5;
    const double dx = 20, dy = 30, dz = 8, top = 1500;
    const size_t n = size_t(nx * ny * nz), nn = size_t((nx + 1) * (ny + 1) * (nz + 1));
    auto nid = [&](int i, int j, int k) { return 3 * ((size_t(k) * (ny + 1) + size_t(j)) * (nx + 1) + size_t(i)); };
    std::vector<double> X(3 * nn);
    for (int k = 0; k <= nz; ++k) for (int j = 0; j <= ny; ++j) for (int i = 0; i <= nx; ++i) { const size_t q = nid(i, j, k); X[q] = i * dx; X[q + 1] = j * dy; X[q + 2] = top + k * dz; }
    std::vector<double> p0(n, 150e5), phi(n, 0.2), dp(n), dT(n, 0.0);
    for (size_t c = 0; c < n; ++c) dp[c] = 1e6 * (1 + 0.3 * std::sin(0.7 * double(c)));
    for (int pad = 0; pad < 2; ++pad) {
        geomech::Input in = c25Input();
        if (pad) { in.overLayers = 2; in.underLayers = 1; in.sidePad = 1; }
        geomech::Solver A, B; std::string e1, e2;
        const bool okA = A.setup(in, nx, ny, nz, {dx}, {dy}, {dz}, top, p0, phi, &e1);
        const bool okB = B.setupHex(in, nx, ny, nz, X, p0, phi, &e2);
        check(okA && okB && B.generalHex(), "brick and isoparametric setups");
        if (!okA || !okB) { std::printf("    %s %s\n", e1.c_str(), e2.c_str()); return; }
        A.solve(dp, dT); B.solve(dp, dT);
        double du = 0, ds = 0, umax = 0, smax = 0;
        for (size_t c = 0; c < n; ++c) {
            du = std::max(du, std::fabs(A.dispZ()[c] - B.dispZ()[c])); umax = std::max(umax, std::fabs(A.dispZ()[c]));
            for (int r = 0; r < 6; ++r) { ds = std::max(ds, std::fabs(A.sigEffFull()[6 * c + size_t(r)] - B.sigEffFull()[6 * c + size_t(r)])); smax = std::max(smax, std::fabs(A.sigEffFull()[6 * c + size_t(r)])); }
        }
        std::printf("    %s: max |du_z| %.3g of %.3g m, max |d sigma'| %.3g of %.3g Pa\n", pad ? "padded" : "no padding", du, umax, ds, smax);
        check(du <= 1e-6 * umax + 1e-12 && ds <= 1e-6 * smax, pad ? "isoparametric = brick elements on a Cartesian box (with padding)" : "isoparametric = brick elements on a Cartesian box");
    }
    // patch test: interior nodes moved by up to 25 % of the spacing, uniform dp -> uniaxial strain, exact for any distortion
    std::vector<double> Y = X;
    unsigned s = 12345u;
    auto rnd = [&]() { s = s * 1103515245u + 12345u; return double((s >> 8) & 0xFFFF) / 65535.0 - 0.5; };
    for (int k = 1; k < nz; ++k) for (int j = 1; j < ny; ++j) for (int i = 1; i < nx; ++i) {
        const size_t q = nid(i, j, k); Y[q] += 0.5 * dx * rnd(); Y[q + 1] += 0.5 * dy * rnd(); Y[q + 2] += 0.5 * dz * rnd();
    }
    geomech::Input in = c25Input();
    geomech::Solver D; std::string e;
    const bool ok = D.setupHex(in, nx, ny, nz, Y, p0, phi, &e);
    check(ok, "distorted isoparametric mesh accepted");
    if (!ok) { std::printf("    %s\n", e.c_str()); return; }
    std::vector<double> sz0(n); for (size_t c = 0; c < n; ++c) sz0[c] = D.sigEffFull()[6 * c + 2];
    const double dpu = 1e6;
    D.solve(std::vector<double>(n, dpu), dT);
    const double E = 10e9, nu = 0.25, M = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
    const double evExact = dpu / M;
    double eErr = 0, sErr = 0;
    for (size_t c = 0; c < n; ++c) {
        eErr = std::max(eErr, std::fabs(D.volStrain()[c] - evExact) / evExact);
        sErr = std::max(sErr, std::fabs((D.sigEffFull()[6 * c + 2] - sz0[c]) + dpu) / dpu);
    }
    std::printf("    distorted mesh, uniform dp 10 bar: eps_v exact %.6g, max rel error %.2g; d sigma'_zz = -dp max rel error %.2g; sigma_min cell 0 %.4g bar\n",
                evExact, eErr, sErr, D.minTotalStress(0) / 1e5);
    check(eErr < 1e-6 && sErr < 1e-6, "patch test: uniaxial-strain solution exact on the distorted mesh (eps_v = alpha dp / M, d sigma'_v = -alpha dp)");
}

static void testC25CornerPointMesh() {
    std::printf("== C25 mechanics mesh of a corner-point grid (GEOMECH_CP.DATA)\n");
    const std::string path = findExample("GEOMECH_CP.DATA");
    if (path.empty()) { std::printf("    (examples/GEOMECH_CP.DATA not found - skipped)\n"); return; }
    const Case c = importEclipseDeck(path);
    check(c.grid.cornerPoint() && c.geoSeq.enabled && c.geoSeq.coupling == 1 && std::fabs(c.geoSeq.permStress - 2e-8) < 1e-20 && c.geoSeq.overLayers == 2,
          "GEOMECH (25 items incl. permStress) read into Case::geoSeq of a black-oil deck");
    check(!c.grid.hfracs.empty() && c.grid.hfracs[0].feStress == 1 && c.faultMech.source == 1, "HFGROW ... FESTRESS and STRESSFIELD ... FE parsed");
    Grid g; g.build(c);
    check(!geoseq::tensorBox(g), "corner-point grid detected as non-tensor (isoparametric mesh)");
    const std::vector<double> X = geoseq::gridNodes(g, 1.0);
    check(X.size() == 3 * size_t(17 * 11 * 5), "node array (nx+1)(ny+1)(nz+1) x 3");
    // the fault pillar: nodes there are the average of the up- and down-thrown corners (throw smeared over one column)
    const size_t qf = 3 * size_t(8);   // node (8, 0, 0)
    std::printf("    node at the fault pillar: x %.1f depth %.3f m (up-thrown 2016.0, down-thrown 2024.0)\n", X[qf], X[qf + 2]);
    check(std::fabs(X[qf + 2] - 2020.0) < 1e-3, "fault-pillar node = mean of the two sides");
    std::vector<double> p(size_t(g.numActive()), 200e5), phi(size_t(g.numActive()), 0.18);
    geoseq::Coupler C; std::string err;
    const bool ok = C.setup(c.geoSeq, g, g.globalOf, 1.0, p, phi, err);
    check(ok && C.hexMesh(), "sequential coupler on the corner-point grid (isoparametric hexahedra)");
    if (!ok) { std::printf("    %s\n", err.c_str()); return; }
    C.step(p);
    const int a = 0;
    double S[6]; C.totalStress(a, S);
    const double z = g.depth[0];
    std::printf("    cell 0 (z %.1f m): total sigma_xx %.2f, sigma_zz %.2f bar (Sv grad %.0f Pa/m -> %.2f), sigma_min %.2f bar\n",
                z, S[0] / 1e5, S[2] / 1e5, c.geoSeq.sigmaVGrad, c.geoSeq.sigmaVGrad * z / 1e5, C.minStress()[0] / 1e5);
    check(std::fabs(S[2] - c.geoSeq.sigmaVGrad * z) < 0.02 * c.geoSeq.sigmaVGrad * z, "initial vertical total stress = overburden gradient x depth (2 %)");
    check(C.minStress()[0] > 0 && C.minStress()[0] < S[2], "sigma_min is the horizontal stress (k0 < 1)");
    // pressurisation -> expansion, uplift
    std::vector<double> p2 = p; for (double& v : p2) v += 20e5;
    C.step(p2);
    const auto pv = C.pvFactor(p2, 5e-5 / 1e5);
    double pmin = 1e9, pmax = 0; for (double v : pv) { pmin = std::min(pmin, v); pmax = std::max(pmax, v); }
    std::printf("    +20 bar: surface uplift %.3f mm, PV factor %.5f .. %.5f, max failure %.3f\n", C.solver().surfaceUplift() * 1e3, pmin, pmax, C.maxFailure());
    check(C.solver().surfaceUplift() > 0 && pmin > 1.0 - 1e-3 && pmax < 1.05, "pressurisation lifts the surface; FE / flow PV factor near 1 (consistent compressibilities)");
}

static void testC25BlackOilCoupled() {
    std::printf("== C25 black oil + corner point: two-way FE geomechanics, FESTRESS fracture growth, fault mechanics with FE stress\n");
    const std::string path = findExample("GEOMECH_CP.DATA");
    if (path.empty()) { std::printf("    (examples/GEOMECH_CP.DATA not found - skipped)\n"); return; }
    struct R { double pAvg = 0, oip0 = 0, oip1 = 0, fopt = 0, fwit = 0; size_t hfEv = 0; double maxSmin = 0, minPv = 1, maxPv = 1, maxDcfs = 0; bool fe = false; double sec = 0; };
    auto run = [&](int coupling, bool fe) {
        R r;
        Case c = importEclipseDeck(path);
        c.geoSeq.coupling = coupling;
        if (!fe) { c.grid.hfracs[0].feStress = 0; c.faultMech.source = 0; }
        const auto t0 = std::chrono::steady_clock::now();
        Simulator sim; sim.setup(c); sim.initialize();
        const Grid& G = sim.grid();
        const double paU = 1e5;
        geoseq::Coupler C; std::string err;
        std::vector<double> p(size_t(G.numActive())), phi(size_t(G.numActive()), 0.18);
        auto pres = [&]() { const auto& P = sim.snapshots().back().pressure; for (size_t a = 0; a < p.size(); ++a) p[a] = double(P[a]) * paU; };
        pres();
        if (!C.setup(c.geoSeq, G, G.globalOf, 1.0, p, phi, err)) { std::printf("    setup: %s\n", err.c_str()); return r; }
        faultmech::Monitor M;
        if (!M.setup(c, G, err)) { std::printf("    faultmech: %s\n", err.c_str()); return r; }
        if (fe) M.setStressProvider([&C](int cell, double S[6]) { return C.totalStress(cell, S); });
        auto couple = [&]() {
            pres();
            C.step(p);
            if (C.twoWay()) {
                const auto m = C.pvFactor(p, c.props.rock.cr / paU);
                for (double v : m) { r.minPv = std::min(r.minPv, v); r.maxPv = std::max(r.maxPv, v); }
                sim.setPoreVolumeMultipliers(m);
            }
            if (fe) { auto sm = C.minStress(); for (double& v : sm) v /= paU; sim.setMinStress(sm); }
            M.step(sim.time(), sim.snapshots().back().pressure);
            // composed transmissibility: fault x geomechanics (matrix connections)
            std::map<int, double> T;
            if (C.twoWay()) {
                std::vector<std::pair<int, int>> pr; std::vector<int> idx;
                for (size_t q = 0; q < G.conns.size(); ++q) {
                    const auto& cn = G.conns[q];
                    if (G.globalOf[size_t(cn.a)] >= G.numCells() || G.globalOf[size_t(cn.b)] >= G.numCells()) continue;
                    pr.push_back({cn.a, cn.b}); idx.push_back(int(q));
                }
                const auto f = C.connFactor(pr);
                for (size_t q = 0; q < f.size(); ++q) if (f[q] != 1.0) T[idx[q]] = f[q];
            }
            for (const auto& e : M.connMultipliers()) { auto it = T.find(e.first); if (it == T.end()) T[e.first] = e.second; else it->second *= e.second; }
            if (!T.empty()) { std::vector<int> cs; std::vector<double> ms; for (const auto& e : T) { cs.push_back(e.first); ms.push_back(e.second); } sim.setConnectionMultipliers(cs, ms); }
        };
        couple();
        r.oip0 = sim.summary().empty() ? 0.0 : sim.summary().front().oip;
        while (sim.runReportStep()) couple();
        const auto& S = sim.summary().back();
        r.oip1 = S.oip; r.fopt = S.fopt; r.fwit = S.fwit;
        double pv = 0, pp = 0; for (int a = 0; a < G.numActive(); ++a) { pv += G.poreVolRef[size_t(a)]; pp += sim.snapshots().back().pressure[size_t(a)] * G.poreVolRef[size_t(a)]; }
        r.pAvg = pp / pv;
        r.hfEv = sim.hfGrowthEvents().size();
        for (double v : sim.minStress()) r.maxSmin = std::max(r.maxSmin, v);
        for (double v : M.out().series[0].maxDcfs) r.maxDcfs = std::max(r.maxDcfs, v);
        r.fe = M.out().feStress;
        r.sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        return r;
    };
    const R one = run(0, true), two = run(1, true), ana = run(1, false);
    std::printf("    one-way: p_avg %.3f bar, %zu growth events, max dCFS %.3f bar (%.1f s)\n", one.pAvg, one.hfEv, one.maxDcfs, one.sec);
    std::printf("    two-way: p_avg %.3f bar, PV factor %.5f .. %.5f, %zu growth events, max sigma_min %.2f bar, max dCFS %.3f bar (%.1f s)\n",
                two.pAvg, two.minPv, two.maxPv, two.hfEv, two.maxSmin, two.maxDcfs, two.sec);
    std::printf("    analytic stress (HFGROW gradient, STRESSFIELD path): %zu growth events, max dCFS %.3f bar\n", ana.hfEv, ana.maxDcfs);
    const double mb = std::fabs((two.oip0 - two.oip1) - two.fopt) / std::max(1.0, two.fopt);
    std::printf("    two-way oil balance: OIP0 - OIP %.6g vs FOPT %.6g sm3 (rel %.2g)\n", two.oip0 - two.oip1, two.fopt, mb);
    check(one.fe && two.fe, "fault mechanics takes the FE stress change (STRESSFIELD ... FE)");
    check(two.maxPv > 1.0 + 1e-5 && std::fabs(two.pAvg - one.pAvg) > 1e-3, "two-way coupling changes the pore volumes and the pressure");
    check(mb < 5e-3, "mass conservation with the geomechanical pore-volume update (oil balance)");
    check(two.maxSmin > 0 && two.hfEv > 0, "HFGROW grows against the FE minimum stress (FESTRESS)");
    check(std::fabs(two.maxDcfs - ana.maxDcfs) > 1e-3, "FE stress change differs from the analytic stress path (stress arching / shear)");
}

static void testC25Comp() {
    std::printf("== C25 CO2 storage (compositional): sequential FE geomechanics, fault mechanics, report-step fracture growth\n");
    const std::string path = findExample("CO2_GEOMECH.DATA");
    if (path.empty()) { std::printf("    (examples/CO2_GEOMECH.DATA not found - skipped)\n"); return; }
    const Case c = importEclipseDeck(path);
    check(c.comp.enabled && c.geoSeq.enabled && c.faultMech.on() && !c.grid.hfracs.empty(), "CO2STORE deck with GEOMECH / FAULTMECH / HFGROW");
    comp::DeckModel m; std::string err;
    if (!comp::deckModel(c, m, &err)) { std::printf("    deckModel: %s\n", err.c_str()); check(false, "compositional model"); return; }
    Case cx = c; expandLgrArrays(cx); expandMincArrays(cx); applyDpnum(cx); prepareHydraulicFractures(cx);
    Grid gx; gx.build(cx);
    check(size_t(gx.numActive()) == m.base.cellPV.size(), "simulator cells = active cells of the expanded grid");
    geoseq::Coupler C; faultmech::Monitor M; geoseq::HfSequential H;
    bool okC = false, okM = false, okH = false, conserved = true, pvSet = false;
    double maxPv = 1, cff0 = 0, cffEnd = 0; int steps = 0;
    const double pu = 1e5;
    std::vector<comp::CompSummary> out;
    const bool ran = comp::runDeckModel(m, out, [&](const comp::CompSummary& S, comp::CompSimulator& sim) {
        const size_t ns = sim.pressure().size();
        if (steps == 0) {
            std::vector<double> phi(ns, 0.3);
            okC = C.setup(c.geoSeq, gx, gx.globalOf, 1.0, sim.pressure(), phi, err);
            okM = M.setup(c, gx, err);
            if (okM && okC) M.setStressProvider([&C](int cell, double S6[6]) { return C.totalStress(cell, S6); });
            okH = H.setup(cx, gx, err);
        }
        if (okC) {
            C.step(sim.pressure());
            double n0 = 0; for (double v : sim.moles()) n0 += v;
            sim.setGeoFactors(C.pvFactor(sim.pressure(), m.base.cr), C.permFactor());
            double n1 = 0; for (double v : sim.moles()) n1 += v;
            conserved = conserved && std::fabs(n1 - n0) <= 1e-12 * std::fabs(n0);
            pvSet = pvSet || !sim.geoPvFactor().empty();
            for (double v : sim.geoPvFactor()) maxPv = std::max(maxPv, v);
        }
        std::vector<float> pc(ns); for (size_t q = 0; q < ns; ++q) pc[q] = float(sim.pressure()[q] / pu);
        if (okM) {
            M.step(S.time, pc);
            const auto& st = M.out().steps.back();
            double mx = -1e300; for (float v : st.cff) mx = std::max(mx, double(v));
            if (steps == 0) cff0 = mx;
        cffEnd = mx;
        }
        if (okH) {
            std::vector<double> pcd(pc.begin(), pc.end());
            std::vector<double> sm = C.minStress(); for (double& v : sm) v /= pu;
            if (steps > 0) H.grow(S.time, pcd, sm);
            std::map<std::pair<int, int>, double> T; H.connMultipliers(T);
            for (const auto& e : M.connMultipliers()) { const auto& cn = gx.conns[size_t(e.first)]; T[{std::min(cn.a, cn.b), std::max(cn.a, cn.b)}] = e.second; }
            std::vector<std::pair<int, int>> pr; std::vector<double> mv; for (const auto& e : T) { pr.push_back(e.first); mv.push_back(e.second); }
            if (!pr.empty()) sim.setConnectionMultipliersByCells(pr, mv);
        }
        ++steps;
        return steps <= 8;
    }, &err, 8);
    std::printf("    %d report steps: PV factor max %.5f, fault faces %zu (FE %d), max CFF %.3f -> %.3f bar, fracture pieces %d open of %d, %zu growth events\n",
                steps, maxPv, okM ? M.out().faces.size() : size_t(0), okM ? int(M.out().feStress) : 0, cff0, cffEnd, okH ? H.openPieces() : 0, okH ? H.pieces() : 0, okH ? H.events().size() : size_t(0));
    if (!ran) std::printf("    run: %s\n", err.c_str());
    check(ran && okC && pvSet && maxPv > 1.0, "compositional run with the FE pore-volume feedback (injection dilates the pores)");
    check(conserved, "moles unchanged by the geomechanical pore-volume update (mass conservative)");
    check(okM && M.out().feStress && M.out().engine == "CO2 storage" && std::isfinite(cffEnd) && cffEnd > cff0, "fault mechanics in the CO2 storage run with the FE stress (Coulomb stress rises with injection)");
    check(okH && H.pieces() > 0 && H.openPieces() > 0, "report-step LEFM fracture model of the compositional engine set up");
}

static void testC25Thermal() {
    std::printf("== C25 EGS (thermal): thermo-elastic fault stress and THM mechanics stress on the fault\n");
    const std::string path = findExample("EGS_FAULT.DATA");
    if (path.empty()) { std::printf("    (examples/EGS_FAULT.DATA not found - skipped)\n"); return; }
    Case c = importEclipseDeck(path);
    check(c.thermal.enabled && c.thermal.geo.enabled && c.faultMech.on() && c.faultMech.source == 1 && std::fabs(c.faultMech.thE - 30) < 1e-12,
          "thermal deck with GEOMECH and STRESSFIELD ... FE TH");
    Grid g; g.build(c);
    {   // analytic source: uniform 50 K cooling at constant pressure -> d sigma_n = E alpha / (1 - nu) dT, tau unchanged
        Case ca = c; ca.faultMech.source = 0;
        faultmech::Monitor M; std::string err;
        const bool ok = M.setup(ca, g, err);
        check(ok && M.out().thermalStress, "analytic source with the thermo-elastic term");
        if (!ok) { std::printf("    %s\n", err.c_str()); return; }
        std::vector<float> p(size_t(g.numActive()), 300.0f), T0(size_t(g.numActive()), 473.15f), T1(size_t(g.numActive()), 423.15f);
        M.step(0, p, &T0); M.step(30, p, &T1);
        const auto& s0 = M.out().steps[0]; const auto& s1 = M.out().steps[1];
        const double coef = 30e9 * 8e-6 / 0.75 / 1e5;
        double dsn = s1.sne[0] - s0.sne[0], dcff = s1.cff[0] - s0.cff[0];
        std::printf("    50 K cooling: d sigma'_n %.4f bar (exact %.4f), dCFF %.4f bar (exact %.4f)\n", dsn, -coef * 50, dcff, 0.6 * coef * 50);
        check(std::fabs(dsn + coef * 50) < 1e-2 && std::fabs(dcff - 0.6 * coef * 50) < 1e-2, "thermo-elastic unclamping of the fault = E alpha dT / (1 - nu)");
    }
    // coupled thermal run, FE (THM) stress
    thermal::DeckModel m; std::string err;
    if (!thermal::deckModel(c, m, &err)) { std::printf("    deckModel: %s\n", err.c_str()); check(false, "thermal model"); return; }
    faultmech::Monitor M;
    const bool okM = M.setup(c, g, err);
    const thermal::Simulator* simT = nullptr;
    if (okM) M.setStressProvider([&simT](int cell, double S6[6]) {
        const geomech::Solver* gm = simT ? simT->mechanics() : nullptr;
        if (!gm) return false;
        const auto& mp = simT->mechCellMap();
        if (cell < 0 || size_t(cell) >= mp.size() || mp[size_t(cell)] < 0) return false;
        gm->totalStress(mp[size_t(cell)], S6); return true;
    });
    int steps = 0; double cff0 = 0, cffEnd = 0, dTmin = 0;
    std::vector<thermal::Summary> out;
    const int na = g.numActive();
    const bool ran = okM && thermal::runDeckModel(m, out, [&](const thermal::Summary& S, thermal::Simulator& sim) {
        simT = &sim;
        std::vector<float> p(static_cast<size_t>(na)), T(static_cast<size_t>(na));
        for (int q = 0; q < na; ++q) { p[size_t(q)] = float(sim.pressure()[size_t(q)] / 1e5); T[size_t(q)] = float(sim.temperature()[size_t(q)]); }
        M.step(S.time, p, &T);
        const auto& st = M.out().steps.back();
        double mx = -1e300; for (float v : st.cff) mx = std::max(mx, double(v));
        if (steps == 0) cff0 = mx;
        cffEnd = mx;
        dTmin = std::min(dTmin, double(*std::min_element(T.begin(), T.end())) - 473.15);
        const auto mult = M.connMultipliers();
        if (!mult.empty()) {
            std::vector<std::pair<int, int>> pr; std::vector<double> mv;
            for (const auto& e : mult) { const auto& cn = g.conns[size_t(e.first)]; pr.push_back({cn.a, cn.b}); mv.push_back(e.second); }
            sim.setConnectionMultipliersByCells(pr, mv);
        }
        ++steps;
        return true;
    }, &err, 8);
    simT = nullptr;
    std::printf("    %d report steps: max CFF %.3f -> %.3f bar (FE %d), %zu faces, expected events %.4g\n", steps, cff0, cffEnd, okM ? int(M.out().feStress) : 0,
                okM ? M.out().faces.size() : size_t(0), okM && !M.out().series.empty() ? M.out().series[0].expected.back() : 0.0);
    if (!ran) std::printf("    run: %s\n", err.c_str());
    check(ran && M.out().feStress && M.out().engine == "thermal" && std::isfinite(cffEnd), "fault mechanics in the thermal run with the THM mechanics stress");
    check(cffEnd > cff0, "injection (pressure + cooling) moves the fault towards failure");
    const Case rt = projectFromText(projectText(c));
    check(rt.thermal.geo.enabled && rt.faultMech.source == 1 && std::fabs(rt.faultMech.thAlpha - 8e-6) < 1e-15, "thermal GEOMECH / STRESSFIELD FE TH survive the deck round trip");
}

static void testC25RoundTrip() {
    std::printf("== C25 deck round trip of GEOMECH (black oil), HFGROW FESTRESS, STRESSFIELD FE\n");
    const std::string path = findExample("GEOMECH_CP.DATA");
    if (path.empty()) { std::printf("    (examples/GEOMECH_CP.DATA not found - skipped)\n"); return; }
    const Case c = importEclipseDeck(path);
    const Case rt = projectFromText(projectText(c));
    check(rt.geoSeq.enabled && rt.geoSeq.coupling == 1 && std::fabs(rt.geoSeq.permStress - 2e-8) < 1e-20 && rt.geoSeq.sidePad == 1 &&
          std::fabs(rt.geoSeq.E[0] - 1e10) < 1, "GEOMECH of a black-oil deck survives the round trip");
    check(!rt.grid.hfracs.empty() && rt.grid.hfracs[0].feStress == 1 && rt.faultMech.source == 1, "FESTRESS / FE tokens survive the round trip");
}
