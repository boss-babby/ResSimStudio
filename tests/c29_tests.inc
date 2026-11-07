// =====================================================================================================================
// C29 tests (remaining items): EWASG-type salt precipitation / dry-out, wellbore friction + ICD, surface network,
// two-way phase-field coupling, device CPR-AMG (float32 emulation of the GPU algorithm)
// =====================================================================================================================

static double c29VermaPruess(double Ss, double phiR, double G) {
    const double phiRel = 1 - Ss;
    if (phiRel <= phiR) return 0;
    const double th = (phiRel - phiR) / (1 - phiR), w = 1 + (1 / G) / (1 / phiR - 1), r = th / (th + w - 1);
    return th * th * (1 - G + G / (w * w)) / (1 - G + G * r * r);
}

// 1-D 20-cell gas-injection column (20 ft cells) for the dry-out tests
static Case c29SaltColumn(double c0, double tempC, bool gasInj) {
    Case c = makeDemoCase();
    GeomInput gi; gi.nx = 20; gi.ny = 1; gi.nz = 1; gi.dx = gi.dy = 20; gi.dz = {20};
    generateGeometry(gi, c.grid);
    const size_t n = size_t(c.grid.size());
    c.grid.poro.assign(n, 0.2); c.grid.permx.assign(n, 100); c.grid.permy.clear(); c.grid.permz.assign(n, 10); c.grid.ntg.clear();
    for (auto& w : c.wells) {
        if (w.type == WellType::Producer) { w.i = 19; w.j = 0; } else { w.i = 0; w.j = 0; if (gasInj) w.type = WellType::GasInjector; }
        w.completions.clear();
        Completion cp; cp.i = w.i; cp.j = w.j; cp.k1 = 0; cp.k2 = 0; w.completions.push_back(cp);
    }
    c.schedule.resize(10);
    for (auto& s : c.schedule) { s.length = 10; s.controls.clear(); }
    WellControl wp; wp.well = "PROD"; wp.mode = Control::ORAT; wp.orat = 300; wp.bhp = 2000;
    WellControl wi; wi.well = "INJ"; wi.mode = Control::RATE; wi.rate = gasInj ? 2000 : 300; wi.bhp = 9000;
    c.schedule[0].controls = {wp, wi};
    c.salt.on = true; c.salt.c0 = c0; c.salt.tempC = tempC;
    return c;
}

static void testC29Salt() {
    std::printf("== C29 salt precipitation / dry-out (TOUGH EWASG-type, black oil)\n");
    const double vf = 0.158987294928;
    // 1) initial supersaturated brine precipitates at t = 0; mass conserved; Verma-Pruess factors
    {
        Case lo = c29SaltColumn(100, 60, false), hi = c29SaltColumn(400, 60, false);
        Simulator a; a.log = [](const std::string&) {}; a.setup(lo); a.initialize();
        Simulator b; b.log = [](const std::string&) {}; b.setup(hi); b.initialize();
        size_t tA = 99, tB = 99;
        for (size_t t = 0; t < a.tracers().size(); ++t) if (a.tracers()[t].name == "NACL") tA = t;
        for (size_t t = 0; t < b.tracers().size(); ++t) if (b.tracers()[t].name == "NACL") tB = t;
        const bool have = tA != 99 && tB != 99;
        const double cap0 = have ? a.tracerInPlace(tA) / 100.0 : 0;   // surface water volume (stb)
        const double total = 400 * cap0 * vf, after = have ? b.tracerInPlace(tB) * vf + b.saltSolidMass() : 0;
        double kerr = 0, smax = 0;
        for (size_t q = 0; q < b.saltSolid().size(); ++q) {
            smax = std::max(smax, b.saltSolid()[q]);
            kerr = std::max(kerr, std::fabs(b.saltPermFactor()[q] - c29VermaPruess(b.saltSolid()[q], 0.8, 0.8)));
        }
        double cmax = 0; if (have) for (int q = 0; q < 20; ++q) cmax = std::max(cmax, b.tracerConc()[tB * 20 + size_t(q)]);
        std::printf("    c0 400 kg/m3 > solubility 317: solid %.6g kg (S_s max %.4g), dissolved max %.6g kg/m3, salt total %.9g / %.9g kg, k-factor error %.1e\n",
                    b.saltSolidMass(), smax, cmax, after, total, kerr);
        check(have && a.saltSolidMass() == 0 && b.saltSolidMass() > 0, "brine above the solubility precipitates at t = 0 (below: none)");
        check(have && std::fabs(cmax - 317) < 1e-9, "dissolved NaCl capped at the solubility");
        check(have && std::fabs(after / total - 1) < 1e-9, "salt mass (dissolved + solid) conserved by the equilibrium step (1e-9)");
        check(kerr < 1e-12 && smax > 0, "permeability factor = Verma-Pruess tubes-in-series formula");
    }
    // 2) dry-out by dry-gas injection at 150 C
    {
        Case c = c29SaltColumn(250, 150, true);
        Case n = c; n.salt.on = false;
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        size_t tr = 99;
        for (size_t t = 0; t < s.tracers().size(); ++t) if (s.tracers()[t].name == "NACL") tr = t;
        const double m0 = tr != 99 ? s.tracerInPlace(tr) * vf : 0;
        const double sw0 = s.snapshots().front().sw[0];
        while (s.runReportStep()) {}
        Simulator r; r.log = [](const std::string&) {}; r.setup(n); r.initialize();
        while (r.runReportStep()) {}
        const double m1 = tr != 99 ? s.tracerInPlace(tr) * vf + s.saltSolidMass() : 0;
        const double mp = s.summary().back().fwpt * vf * 250.0;   // produced brine: connate water at c0 (the dry-out front stays near the injector)
        const double sw1 = s.snapshots().back().sw[0], swN = r.snapshots().back().sw[0];
        const double S0 = s.saltSolid()[0], k0 = s.saltPermFactor()[0];
        const double gi = s.summary().back().fgit, giN = r.summary().back().fgit;
        // complete dry-out: all the salt of the injection cell becomes solid, S_s = c0 V_w / (rho_s PV) at the initial state
        std::printf("    injector cell: S_w %.4f -> %.3g (no salt model %.4f), solid S_s %.4g, k factor %.4g | evaporated %.5g stb | salt %.9g -> %.9g kg"
                    " + produced %.6g | gas injected %.6g vs %.6g MMscf\n", sw0, sw1, swN, S0, k0, s.saltEvaporated(), m0, m1, mp, gi / 1e3, giN / 1e3);
        check(tr != 99 && s.saltEvaporated() > 0 && sw1 < 0.25 * sw0 && swN > 0.5 * sw0, "dry gas evaporates the water of the injection cell (local equilibrium sink)");
        check(S0 > 0 && k0 < 1 && std::fabs(k0 - c29VermaPruess(S0, 0.8, 0.8)) < 1e-12, "dry-out precipitates salt and reduces the permeability (Verma-Pruess)");
        check(m0 > 0 && std::fabs((m1 + mp) / m0 - 1) < 1e-6, "salt mass conserved through transport + evaporation + precipitation (in place + solid + produced, 1e-6)");
        // restart round trip of the salt state
        bool rok = false;
        try {
            Simulator a; a.log = [](const std::string&) {}; a.setup(c); a.initialize();
            for (int k = 0; k < 4; ++k) a.runReportStep();
            const std::string st = a.saveState();
            Simulator b; b.log = [](const std::string&) {}; b.setup(c); b.initialize(); b.loadState(st);
            while (a.runReportStep()) {}
            while (b.runReportStep()) {}
            const double d = std::fabs(a.saltSolidMass() - b.saltSolidMass()) / std::max(1e-30, a.saltSolidMass());
            const double dp = std::fabs(a.summary().back().fpr - b.summary().back().fpr) / a.summary().back().fpr;
            std::printf("    restart after step 4: solid %.9g vs %.9g kg, FPR rel diff %.1e\n", a.saltSolidMass(), b.saltSolidMass(), dp);
            rok = d < 1e-6 && dp < 1e-6;
        } catch (const std::exception& e) { std::printf("    restart exception: %s\n", e.what()); }
        check(rok, "salt state (solid, evaporation sink) survives save / load (state v5)");
    }
    // 3) deck round trip
    {
        Case c = c29SaltColumn(250, 150, true);
        c.salt.cSol = 300; c.salt.phiR = 0.85; c.salt.gamma = 0.7; c.salt.cInj = 12;
        Case::WellboreSpec B; B.well = "PROD"; B.id = 0.25; B.rough = 1e-4; B.icdK = 0.003; B.rhoCal = 62.4;
        c.wellbores.push_back(B);
        const std::string path = "/tmp/c29_salt.DATA";
        exportEclipseDeck(c, path);
        const Case d = importEclipseDeck(path);
        const bool ok = d.salt.on && d.salt.c0 == 250 && d.salt.cSol == 300 && d.salt.phiR == 0.85 && d.salt.gamma == 0.7 && d.salt.tempC == 150 &&
                        d.salt.cInj == 12 && d.wellbores.size() == 1 && d.wellbores[0].well == "PROD" && d.wellbores[0].id == 0.25 &&
                        d.wellbores[0].icdK == 0.003 && d.wellbores[0].rhoCal == 62.4;
        check(ok, "RESSIM SALT / WELLBORE deck round trip");
    }
}

static void testC29Wellbore() {
    std::printf("== C29 wellbore friction + inflow-control devices (FlexWell-lite)\n");
    auto run = [](const Case& c, double& fopt, double& dpw, double& resv, double& bhp) {
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        while (s.runReportStep()) {}
        fopt = s.summary().back().fopt;
        size_t ip = 0; for (size_t w = 0; w < c.wells.size(); ++w) if (c.wells[w].name == "PROD") ip = w;
        dpw = s.summary().back().wells[ip].dpWellbore; resv = s.summary().back().wells[ip].resv; bhp = s.summary().back().wells[ip].bhp;
    };
    Case base = smallCase();
    base.schedule.resize(8);
    for (auto& st : base.schedule) for (auto& wc : st.controls) if (wc.well == "PROD") { wc.mode = Control::BHP; wc.bhp = 3000; wc.orat = 20000; }
    double f0, d0, q0, b0, f1, d1, q1, b1, f2, d2, q2, b2;
    run(base, f0, d0, q0, b0);
    Case fr = base; { Case::WellboreSpec B; B.well = "PROD"; B.id = 0.1; B.rough = 5e-5; fr.wellbores.push_back(B); }
    run(fr, f1, d1, q1, b1);
    // single connection at the datum (no friction length): ICD drop only
    Case ic = base;
    for (auto& w : ic.wells) if (w.name == "PROD") for (auto& cp : w.completions) cp.k2 = cp.k1;
    double fi0, di0, qi0, bi0;
    run(ic, fi0, di0, qi0, bi0);
    const double K = 2e-5;
    { Case::WellboreSpec B; B.well = "PROD"; B.id = 2.0; B.icdK = K; B.rhoCal = 62.4; ic.wellbores.push_back(B); }
    run(ic, f2, d2, q2, b2);
    const double ratio = d2 / (K * q2 * q2);   // = rho_mix / rho_cal (oil-dominated inflow: ~0.6 - 1.0)
    std::printf("    BHP 3000 psi, 8 steps: FOPT %.6g (no wellbore) vs %.6g (ID 0.1 ft: friction %.4g psi at the deepest connection)\n", f0, f1, d1);
    std::printf("    single connection + ICD K %.1e psi/(rb/d)^2: FOPT %.6g -> %.6g, ICD drop %.4g psi at q %.5g rb/d -> rho_mix / rho_cal %.4f\n", K, fi0, f2, d2, q2, ratio);
    check(d0 == 0 && d1 > 1 && f1 < f0, "tubing friction raises the connection pressures of a producer on BHP control: less oil");
    check(d2 > 0 && f2 < fi0 && ratio > 0.5 && ratio < 1.1, "ICD drop = K rho_mix / rho_cal q^2 (single connection, oil density ratio)");
}

static void testC29Network() {
    std::printf("== C29 surface network (NODEPROP / BRANPROP extended network)\n");
    Case base = smallCase();
    base.schedule.resize(10);
    for (auto& st : base.schedule) for (auto& wc : st.controls) if (wc.well == "PROD") { wc.orat = 9000; wc.bhp = 500; }
    // synthetic pipeline table 7: p_in = p_out + 0.05 q_oil (exact under multilinear interpolation)
    VfpTable T; T.id = 7; T.flo = "OIL"; T.wfr = "WCT"; T.gfr = "GOR";
    T.floV = {0, 50000}; T.thpV = {0, 10000}; T.wfrV = {0}; T.gfrV = {0}; T.alqV = {0};
    for (double thp : T.thpV) for (double flo : T.floV) T.bhp.push_back(thp + 0.05 * flo);
    Case net = base;
    net.vfp.push_back(T);
    std::string grp = "G";
    for (const auto& w : net.wells) if (w.name == "PROD") grp = w.group;
    ExtraKeyword np; np.name = "NODEPROP"; np.section = DeckSection::Schedule; np.step = 0; np.body = " 'FIELD' 600 /\n/\n";
    ExtraKeyword bp; bp.name = "BRANPROP"; bp.section = DeckSection::Schedule; bp.step = 0; bp.body = " '" + grp + "' 'FIELD' 7 /\n/\n";
    net.extraKeywords.push_back(np); net.extraKeywords.push_back(bp);
    auto go = [](const Case& c, double& fopt, std::vector<Simulator::NetNode>& nodes, double& thp, double& bhp, double& qo, double& qFirst) {
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        while (s.runReportStep()) {}
        fopt = s.summary().back().fopt; nodes = s.networkNodes();
        size_t ip = 0; for (size_t w = 0; w < c.wells.size(); ++w) if (c.wells[w].name == "PROD") ip = w;
        thp = s.summary().back().wells[ip].thpNet; bhp = s.summary().back().wells[ip].bhp; qo = s.summary().back().wells[ip].oil;
        qFirst = s.summary().front().wells[ip].oil; if (s.summary().size() > 1) qFirst = s.summary()[1].wells[ip].oil;
    };
    double f0, f1, t0, t1, b0, b1, q0, q1, e0, e1;
    std::vector<Simulator::NetNode> n0, n1;
    go(base, f0, n0, t0, b0, q0, e0);
    go(net, f1, n1, t1, b1, q1, e1);
    if (std::getenv("C29V")) {
        for (double bl : {1000.0, 1858.0, 2500.0}) {
            Case b2 = base; for (auto& st : b2.schedule) for (auto& wc : st.controls) if (wc.well == "PROD") wc.bhp = bl;
            double f2, t2, bb2, q2, e2; std::vector<Simulator::NetNode> n2; go(b2, f2, n2, t2, bb2, q2, e2);
            std::printf("      | no network, BHP limit %.0f: FOPT %.6g, final %.5g stb/d\n", bl, f2, q2);
        }
    }
    double pG = -1, qG = 0, pF = -1;
    for (const auto& N : n1) { if (N.name == grp) { pG = N.p; qG = N.qo; } if (N.name == "FIELD") pF = N.p; }
    std::printf("    no network: FOPT %.6g, PROD %.5g stb/d at BHP %.5g | network: FOPT %.6g, node %s %.6g psi (= 600 + 0.05 x %.6g), PROD THP limit %.6g, %.5g stb/d at BHP %.5g | first-step oil rate %.5g vs %.5g\n",
                f0, q0, b0, f1, grp.c_str(), pG, qG, t1, q1, b1, e0, e1);
    check(n0.empty() && n1.size() == 2 && pF == 600, "network nodes from NODEPROP / BRANPROP (fixed terminal pressure)");
    // C30: the balance is implicit - the node pressure applied to the wells agrees with VFP(q) of the rates they deliver at it
    // within the balance tolerance (0.2 %); C29 (explicit) computed it from the lagged rates exactly
    check(std::fabs(pG - (600 + 0.05 * qG)) < 2.5e-3 * pG, "branch pressure = VFP pipeline table of the delivered rate (p_down = VFP(q, p_up), implicit balance tolerance)");
    check(std::fabs(t1 - pG) < 1e-9 * pG && b1 > pG && b1 > 500 + 1 && q1 > q0, "producer constrained by the node pressure (THP limit -> BHP above its own 500 psi limit; less free gas, more late oil)");
}

static void testC29PfTwoWay() {
    std::printf("== C29 two-way phase-field coupling (HFGROW ... PHASEFIELD eta TWOWAY)\n");
    const char* rate = std::getenv("C29RATE") ? std::getenv("C29RATE") : "0.6";
    auto deck = [&](const std::string& grow) {
        std::string d = hfDeck(true, false, grow);
        const std::string a = "WCONPROD\n P OPEN ORAT 5 4* 10 /";
        const size_t at = d.find(a);
        if (at != std::string::npos) d.replace(at, a.size(), std::string("WCONINJE\n P WATER OPEN RATE ") + rate + " 1* 700 /");
        const size_t q = d.find(" 0.2 0 1 0\n");
        if (q != std::string::npos) d.replace(q, 11, " 0.2 0 1 0\n 0.3 0.02 0.7 0\n");
        const size_t ts = d.find("TSTEP\n");
        if (ts != std::string::npos) d.replace(ts, d.find("END") - ts, "TSTEP\n 4*1 6*4 6*10 /\n");
        return d;
    };
    bool ok = false; std::string err; std::vector<Simulator::PfRecord> rec; int twoLog = 0; double wl = 0;
    bool rt = false;
    const auto t0 = std::chrono::steady_clock::now();
    try {
        Case c = eclipseDeckFromString(deck("-- RESSIM: HFGROW F1 1.0 320 0 20 1 PHASEFIELD 0.002 TWOWAY\n"), "");
        rt = !c.grid.hfracs.empty() && c.grid.hfracs[0].pfTwoWay == 1 && projectText(c).find("TWOWAY") != std::string::npos;
        Simulator s; s.log = [&](const std::string& m) { if (m.find("(two-way)") != std::string::npos) ++twoLog; if (std::getenv("C29V") && m.find("phase field") != std::string::npos) std::printf("      | %s\n", m.c_str()); };
        s.setup(c); s.initialize();
        while (s.runReportStep()) {}
        rec = s.phaseFieldRecords(); wl = s.hfWingLength(0, 1);
        ok = true;
    } catch (const std::exception& e) { err = e.what(); }
    // the same deck one-way (eta volume): crack pressure vs flow-model pressure at the end
    double oneA = 0, oneDp = -1;
    try {
        Case c = eclipseDeckFromString(deck("-- RESSIM: HFGROW F1 1.0 320 0 20 1 PHASEFIELD 0.002\n"), "");
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        while (s.runReportStep()) {}
        for (const auto& r : s.phaseFieldRecords()) if (r.active && r.volume > 0) { oneA = r.halfLength; oneDp = std::fabs(r.pressure * 1e-5 - r.pFlow); }
    } catch (const std::exception&) {}
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (!ok) std::printf("    exception: %s\n", err.c_str());
    int two = 0, act = 0; bool mono = true, vmono = true; double lastA = 0, lastV = 0, maxDev = 0, lastDp = -1;
    for (const auto& r : rec) {
        if (!r.active || r.volume <= 0) continue;
        ++act; if (r.note == "two-way") ++two;
        mono = mono && r.halfLength >= lastA - 1e-9; lastA = r.halfLength;
        vmono = vmono && r.volume >= lastV * (1 - 1e-6); lastV = r.volume;
        if (r.note == "two-way" && r.pressure > 0) { maxDev = std::max(maxDev, std::fabs(r.pressure - r.pFlow * 1e5) / r.pressure); lastDp = std::fabs(r.pressure * 1e-5 - r.pFlow); }
    }
    std::printf("    rate 0.6 m3/d, %.1f s (both runs): %d active phase-field steps (%d pressure-driven), final half-length %.2f m vs one-way (eta 0.002) %.2f m;"
                " |p_crack - p_flow| at the last step %.2f bar vs one-way %.2f bar; wing %.1f m, %d aperture updates\n",
                sec, act, two, lastA, oneA, lastDp, oneDp, wl, twoLog);
    check(rt, "TWOWAY flag parsed and written back");
    check(ok && two >= 3 && twoLog >= 3, "pressure-driven crack volume after the first step + conductance feedback (Sneddon aperture)");
    check(ok && mono && vmono, "irreversible crack: volume and half-length never decrease");
    check(lastDp >= 0 && oneDp > 0 && lastDp < 2.0 && lastDp < 0.25 * oneDp, "two-way closes the crack / flow pressure gap (< 2 bar, < 1/4 of one-way)");
}

static void testC29DeviceCprAmg() {
    std::printf("== C29 device CPR-AMG (float32 emulation of the GPU V-cycle)\n");
    const BSR A = makeBlockLaplace(120, 100, 3, 13);
    const int N = A.size();
    std::vector<double> xt(static_cast<size_t>(N)), b(static_cast<size_t>(N));
    for (int e = 0; e < N; ++e) xt[size_t(e)] = std::cos(0.003 * e) + 0.2 * (e % 3);
    A.spmv(xt.data(), b.data());
    DeviceCprAmg D;
    const bool built = D.build(A);
    // coarsest inverse check
    double cinv = 1;
    if (built && !D.coarseInv.rp.empty()) {
        const DevCsr& C = D.lv.back().A;
        std::vector<float> v(size_t(C.rows)), Av(size_t(C.rows)), w(size_t(C.rows));
        for (int i = 0; i < C.rows; ++i) v[size_t(i)] = float(std::sin(0.37 * i) + 1.5);
        C.spmv(v.data(), Av.data()); D.coarseInv.spmv(Av.data(), w.data());
        double e = 0, nv = 0; for (int i = 0; i < C.rows; ++i) { e = std::max(e, double(std::fabs(w[size_t(i)] - v[size_t(i)]))); nv = std::max(nv, double(std::fabs(v[size_t(i)]))); }
        cinv = e / nv;
    }
    CpuFloatBlockSolver bj; bj.multicolour = true; CpuFloatBlockSolver cp; cp.cprAmg = true;
    std::vector<double> x1, x2;
    auto t0 = std::chrono::steady_clock::now();
    const LinResult r1 = refineWithDevice(bj, A, b, x1, 1e-10, 8, 3000);
    auto t1 = std::chrono::steady_clock::now();
    const LinResult r2 = refineWithDevice(cp, A, b, x2, 1e-10, 8, 3000);
    auto t2 = std::chrono::steady_clock::now();
    double e1 = 0, e2 = 0, xn = 0;
    for (int e = 0; e < N; ++e) { e1 = std::max(e1, std::fabs(x1[size_t(e)] - xt[size_t(e)])); e2 = std::max(e2, std::fabs(x2[size_t(e)] - xt[size_t(e)])); xn = std::max(xn, std::fabs(xt[size_t(e)])); }
    std::printf("    120x100 nb3 block system: %zu AMG levels (coarsest %d, dense inverse error %.1e) | MC-BILU0 emulation %d its %.0f ms err %.1e | CPR-AMG emulation %d its %.0f ms err %.1e\n",
                D.lv.size(), D.lv.empty() ? 0 : D.lv.back().A.rows, cinv, r1.iterations, 1e3 * std::chrono::duration<double>(t1 - t0).count(), e1 / xn,
                r2.iterations, 1e3 * std::chrono::duration<double>(t2 - t1).count(), e2 / xn);
    check(built && D.lv.size() >= 2 && D.W.rows == A.n && D.E.rows == N && cinv < 1e-3, "device hierarchy: W / E restriction, AMG levels, explicit coarsest inverse");
    check(r2.converged && e2 / xn < 1e-8, "float32 CPR-AMG device algorithm + double refinement reaches double accuracy");
    // inside the Newton loop: CPU CPR, device MC-BILU0 emulation, device CPR-AMG emulation
    Case c = makeDemoCase();
    c.schedule.resize(6);
    c.tuning.blockSolver = true;
    double fopt[3]; RunStats st[3];
    CpuFloatBlockSolver mcd; mcd.multicolour = true;
    for (int k = 0; k < 3; ++k) {
        Case ck = c; ck.tuning.gpu = (k > 0);
        setGpuBlockSolver(k == 1 ? &mcd : k == 2 ? &cp : nullptr);
        Simulator sim; sim.log = [](const std::string&) {}; sim.setup(ck); sim.initialize();
        while (sim.runReportStep()) {}
        fopt[k] = sim.summary().back().fopt; st[k] = sim.stats();
    }
    setGpuBlockSolver(nullptr);
    std::printf("    demo 6 steps: FOPT %.8g (CPU CPR) / %.8g (device MC-BILU0) / %.8g (device CPR-AMG); linear its %d / %d / %d; device solves %d / %d, fallbacks %d / %d\n",
                fopt[0], fopt[1], fopt[2], st[0].linear, st[1].linear, st[2].linear, st[1].gpuSolves, st[2].gpuSolves, st[1].gpuFallbacks, st[2].gpuFallbacks);
    check(std::fabs(fopt[2] / fopt[0] - 1) < 1e-6 && st[2].gpuSolves > 0 && st[2].gpuFallbacks == 0, "Newton loop with the device CPR-AMG: same result, no fallbacks");
    check(st[2].linear * 3 < st[1].linear, "device CPR-AMG needs > 3x fewer linear iterations than the device MC-BILU0 in the Newton loop");
}
