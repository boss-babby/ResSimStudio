// ===================================================================================== C18 tests (included by run_tests.cpp)

// compositional discrete adjoint vs finite differences of re-simulations on the same step sequence
static void testCompAdjoint() {
    using namespace ressim::eos;
    using namespace ressim::comp;
    std::printf("Compositional discrete adjoint (C18)\n");
    const std::vector<double> zoil = {0.5, 0.03, 0.07, 0.2, 0.15, 0.05}, zgas = {0.77, 0.20, 0.03, 0, 0, 0};
    const double T = fahrenheitToK(160);
    for (int variant = 0; variant < 2; ++variant) {
    // variant 0: one layer, no gravity, BHP producer; variant 1: two layers with gravity (lagged wellbore heads), a
    // BHP-limited injector and a stock-tank-oil-rate producer (lagged separator yields in the rate constraint) - the lagged
    // well data are differentiated as functions of the previous state, so both are exact discrete adjoints
    const int nx = 6, ny = 4;
    const int nz = variant == 1 ? 2 : 1, nn = nx * ny * nz;
    CompCase c; c.nx = nx; c.ny = ny; c.nz = nz; c.dx = {30}; c.dy = {30}; c.dz = {8}; c.poro = {0.2};
    c.permx.resize(size_t(nn)); c.permy.resize(size_t(nn)); c.permz = {3e-14};
    for (int q = 0; q < nn; ++q) { c.permx[size_t(q)] = c.permy[size_t(q)] = 1e-13 * (1.0 + 0.6 * std::sin(1.7 * q)); }
    c.top = 2500; c.gravity = variant == 1; c.fluid = spe5Fluid(); c.T = T; c.cr = 5e-10; c.pref = 2e7;
    c.p0 = {psiaToPa(2600)}; c.z0 = {zoil}; c.separator = {{300.0, 7e5}, {288.71, 101325.0}};
    CompWell inj; inj.name = "I"; inj.injector = true; inj.cells = variant == 1 ? std::vector<int>{0, nx * ny} : std::vector<int>{0}; inj.wi = variant == 1 ? std::vector<double>{3e-12, 3e-12} : std::vector<double>{3e-12}; inj.zinj = zgas; inj.rate = 60.0; inj.bhpMax = psiaToPa(4500);
    CompWell prd; prd.name = "P"; prd.cells = variant == 1 ? std::vector<int>{nx * ny - 1, nn - 1} : std::vector<int>{nx * ny - 1}; prd.wi = variant == 1 ? std::vector<double>{3e-12, 3e-12} : std::vector<double>{3e-12}; prd.bhp = psiaToPa(1800); prd.control = variant == 1 ? 3 : 0; prd.target = 40.0 / 86400.0;
    c.wells = {inj, prd};
    std::vector<int> left, right;
    for (int q = 0; q < nn; ++q) ((q % nx) < nx / 2 ? left : right).push_back(q);
    std::vector<CompAdjParam> params = {compParamPermMult("perm left", left), compParamPermMult("perm right", right),
                                        compParamPoroMult("poro left", left), compParamWellIndexMult("WI P", 1)};
    const double tEnd = 120 * 86400.0;
    auto run = [&](const CompCase& cc, const std::vector<double>* forced, CompSimulator& sim) {
        sim.setup(cc);
        sim.setNewtonTolerance(1e-11, 1e-11, false, 1e-10);
        if (forced) sim.setForcedSteps(*forced);
        sim.setAdjointRecording(true);
        bool ok = true;
        for (int r = 1; r <= 12 && ok; ++r) ok = sim.advanceTo(r * 10 * 86400.0, 8 * 86400.0);
        return ok;
    };
    CompSimulator base; const bool ok0 = run(c, nullptr, base);
    // synthetic "history": the base response scaled (a misfit with non-zero gradient)
    std::vector<CompAdjObs> obs;
    {
        CompSimulator probe; std::vector<double> dts = base.stepSizes(); run(c, &dts, probe);
        for (int r = 1; r <= 12; ++r)
            for (int qty : {0, 1, 3, 10}) {
                CompAdjObs o; o.well = qty == 3 ? 0 : 1; o.qty = qty; o.t = r * 10 * 86400.0;
                o.value = 0; o.sigma = 1;
                obs.push_back(o);
            }
        // value = 0.9 x simulated, sigma = 5 % of it
        std::vector<CompAdjObs> one(1);
        for (auto& o : obs) { one[0] = o; one[0].value = 0; one[0].sigma = 1; const double y = std::sqrt(probe.adjointMisfit(one)); o.value = 0.9 * y; o.sigma = std::max(1e-30, 0.05 * y); }
    }
    if (variant == 1) { const auto& W = base.wellReports(); std::printf("    end: inj BHP %.1f psia onBhp %d, prod BHP %.1f psia onBhp %d qo %.3f sm3/d qg %.1f\n", W[0].bhp / 6894.757, int(W[0].onBhp), W[1].bhp / 6894.757, int(W[1].onBhp), W[1].qo * 86400, W[1].qg * 86400); }
    const double F0 = base.adjointMisfit(obs);
    const CompAdjResult g = base.adjointGradient(obs, params);
    std::printf("  base run %s: %d steps, F = %.4f; adjoint %s in %.2f s (%d linear its)\n", ok0 ? "ok" : "FAILED", base.adjointSteps(), F0,
                g.ok ? "ok" : g.error.c_str(), g.seconds, g.linear);
    const std::vector<double> dts = base.stepSizes();
    double worst = 0, gmax = 0;
    for (double v : g.grad) gmax = std::max(gmax, std::fabs(v));
    for (size_t j = 0; j < params.size(); ++j) {
        const double h = 1e-3;
        CompCase cp = c, cm = c; params[j].apply(cp, 1 + h); params[j].apply(cm, 1 - h);
        CompSimulator sp, sm; const bool okp = run(cp, &dts, sp), okm = run(cm, &dts, sm);
        const double fd = (sp.adjointMisfit(obs) - sm.adjointMisfit(obs)) / (2 * h);
        const double rel = std::fabs(g.grad[j] - fd) / std::max(std::fabs(fd), 1e-6 * gmax);   // (an insensitive parameter: relative to the largest gradient)
        worst = std::max(worst, rel);
        std::printf("    dF/d(%s): adjoint %+.6e  FD %+.6e  rel.err %.2e%s\n", params[j].name.c_str(), g.grad[j], fd, rel, okp && okm ? "" : " (FD run failed)");
    }
    if (variant == 0) check(ok0 && g.ok && F0 > 0 && worst < 1e-4, "compositional adjoint gradient (perm / poro / WI multipliers) matches finite differences of re-simulations");
    else check(ok0 && g.ok && F0 > 0 && worst < 1e-4, "compositional adjoint with gravity heads, BHP-limited injector, ORAT producer (lagged well data differentiated) matches finite differences");
    }
}

// thermal discrete adjoint (hot-water flood of a heavy oil, 2-D) vs finite differences of re-simulations
static void testThermAdjoint() {
    std::printf("Thermal discrete adjoint (C18)\n");
    using ressim::thermal::Case; using ressim::thermal::Simulator; using ressim::thermal::Well; using ressim::thermal::AdjObs;
    using ressim::thermal::AdjParam; using ressim::thermal::AdjResult; using ressim::thermal::paramPermMult; using ressim::thermal::paramPoroMult;
    using ressim::thermal::paramWellIndexMult; using ressim::thermal::paramConductivityMult;
    const int nx = 5, ny = 3, nn = nx * ny;
    Case c; c.nx = nx; c.ny = ny; c.nz = 1; c.dx = {10}; c.dy = {10}; c.dz = {5}; c.poro = {0.3}; c.gravity = false; c.cr = 5e-10;
    c.permx.resize(size_t(nn)); c.permy.resize(size_t(nn)); c.permz = {1e-12};
    for (int q = 0; q < nn; ++q) c.permx[size_t(q)] = c.permy[size_t(q)] = 3e-12 * (1.0 + 0.5 * std::sin(1.3 * q));
    c.comps = {thermHeavy(1.0)}; c.p0 = {3e6}; c.T0 = {300}; c.sw0 = {0.3}; c.x0 = {{1.0}};
    Well inj; inj.name = "INJ"; inj.injector = true; inj.cells = {0}; inj.wi = {5e-12}; inj.zinj = {0, 1}; inj.bhpMax = 9e6;
    inj.rate = 20.0 / 86400 * 1000 / ressim::steam::kMw; inj.tinj = 450;   // 20 m3/d of hot water (single phase)
    Well pr; pr.name = "PROD"; pr.cells = {nn - 1}; pr.wi = {5e-12}; pr.bhp = 2e6;
    c.wells = {inj, pr};
    std::vector<int> left, right;
    for (int q = 0; q < nn; ++q) ((q % nx) < 2 ? left : right).push_back(q);
    std::vector<AdjParam> params = {paramPermMult("perm left", left), paramPermMult("perm right", right), paramPoroMult("poro left", left),
                                    paramWellIndexMult("WI PROD", 1), paramConductivityMult("conductivity")};
    auto run1 = [&](const Case& cc, const std::vector<double>* forced, Simulator& sim, double tol) {
        std::string err;
        if (!sim.setup(cc, &err)) return false;
        sim.setNewtonTolerance(tol, tol * 1e3, tol);
        if (forced) sim.setForcedSteps(*forced);
        sim.setAdjointRecording(true);
        bool ok = true;
        for (int r = 1; r <= 8 && ok; ++r) ok = sim.advanceTo(r * 15 * 86400.0, 10 * 86400.0);
        return ok;
    };
    // forced-step FD runs cannot cut the step: when the very tight tolerance is not reachable in floating point (seen with
    // the Windows build: different libm / FMA rounding), retry the run with a 10x / 100x looser tolerance
    auto run = [&](const Case& cc, const std::vector<double>* forced, Simulator& sim) {
        for (double tol : {1e-11, 1e-10, 1e-9}) {
            if (run1(cc, forced, sim, tol)) return true;
            sim = Simulator();
        }
        return false;
    };
    Simulator base; const bool ok0 = run(c, nullptr, base);
    std::vector<AdjObs> obs;
    for (int r = 1; r <= 8; ++r)
        for (int qty : {0, 2, 4, 6, 3}) { AdjObs o; o.well = qty == 3 ? 0 : 1; o.qty = qty; o.t = r * 15 * 86400.0; obs.push_back(o); }
    {   // synthetic history: 0.9 x simulated, sigma 5 %
        std::vector<AdjObs> one(1);
        for (auto& o : obs) { one[0] = o; one[0].value = 0; one[0].sigma = 1; const double y = std::sqrt(base.adjointMisfit(one)); o.value = 0.9 * y; o.sigma = std::max(1e-30, 0.05 * y); }
    }
    const double F0 = base.adjointMisfit(obs);
    const AdjResult g = base.adjointGradient(obs, params);
    const auto& W = base.wellReports();
    std::printf("  base run %s: %d steps, producer T %.1f K, qo %.3f m3/d, qw %.3f m3/d; F = %.3f; adjoint %s in %.2f s (%d linear its)\n",
                ok0 ? "ok" : "FAILED", base.adjointSteps(), W[1].tFluid, W[1].qo * 86400, W[1].qw * 86400, F0, g.ok ? "ok" : g.error.c_str(), g.seconds, g.linear);
    const std::vector<double> dts = base.stepSizes();
    double worst = 0, gmax = 0;
    for (double v : g.grad) gmax = std::max(gmax, std::fabs(v));
    for (size_t j = 0; j < params.size(); ++j) {
        const double h = 1e-3;
        Case cp = c, cm = c; params[j].apply(cp, 1 + h); params[j].apply(cm, 1 - h);
        Simulator sp, sm; const bool okp = run(cp, &dts, sp), okm = run(cm, &dts, sm);
        const double fd = (sp.adjointMisfit(obs) - sm.adjointMisfit(obs)) / (2 * h);
        const double rel = std::fabs(g.grad[j] - fd) / std::max(std::fabs(fd), 1e-6 * gmax);
        worst = std::max(worst, rel);
        std::printf("    dF/d(%s): adjoint %+.6e  FD %+.6e  rel.err %.2e%s\n", params[j].name.c_str(), g.grad[j], fd, rel, okp && okm ? "" : " (FD run failed)");
    }
    check(ok0 && g.ok && F0 > 0 && worst < 2e-3, "thermal adjoint gradient (perm / poro / WI / conductivity) matches finite differences within 2e-3");
}

// C18: multicolour block ILU(0) (CPU OpenMP version and the float32 emulation of the GPU kernels)
static void testMulticolourIlu() {
    std::printf("== C18 multicolour BILU0 (parallel ILU: CPU colours + GPU kernel emulation)\n");
    const BSR A = makeBlockLaplace(120, 100, 3, 11);
    // irregular variant: add a few long-range couplings (wells / NNCs) -> more colours
    BSR B = A;
    {
        std::vector<std::vector<std::pair<int, std::vector<double>>>> rows(static_cast<size_t>(B.n));
        const int b2 = 9;
        for (int i = 0; i < B.n; ++i)
            for (int p = B.rowPtr[size_t(i)]; p < B.rowPtr[size_t(i) + 1]; ++p)
                rows[size_t(i)].push_back({B.col[size_t(p)], std::vector<double>(B.val.begin() + size_t(p) * b2, B.val.begin() + size_t(p + 1) * b2)});
        auto add = [&](int i, int j, double t) {
            for (auto* r : {&rows[size_t(i)], &rows[size_t(j)]}) {
                const int o = (r == &rows[size_t(i)]) ? j : i, s = (r == &rows[size_t(i)]) ? i : j;
                std::vector<double> off(9, 0.0); for (int e = 0; e < 3; ++e) off[size_t(e) * 4] = -t;
                r->push_back({o, off});
                for (auto& q : *r) if (q.first == s) for (int e = 0; e < 3; ++e) q.second[size_t(e) * 4] += 1.02 * t;
            }
        };
        for (int w = 0; w < 6; ++w) for (int k = 1; k <= 8; ++k) add(1000 * w + 37, 1000 * w + 37 + 117 * k, 0.5);
        B.col.clear(); B.val.clear(); B.rowPtr.assign(1, 0);
        for (auto& r : rows) {
            std::sort(r.begin(), r.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
            for (auto& q : r) { B.col.push_back(q.first); B.val.insert(B.val.end(), q.second.begin(), q.second.end()); }
            B.rowPtr.push_back(int(B.col.size()));
        }
        B.finalizeDiag();
    }
    for (int v = 0; v < 2; ++v) {
        const BSR& M = v ? B : A;
        const int N = M.size();
        std::vector<double> xt(static_cast<size_t>(N)), b(static_cast<size_t>(N));
        for (int e = 0; e < N; ++e) xt[size_t(e)] = std::sin(0.013 * e) + 0.1 * (e % 3);
        M.spmv(xt.data(), b.data());
        BsrOperator op(M);
        BILU0 P0; P0.factor(M);
        MCBILU0 P1; P1.factor(M);
        // colouring validity: no entry couples two cells of one colour
        const Colouring& C = P1.colouring();
        bool valid = true;
        for (int i = 0; i < M.n; ++i) for (int p = M.rowPtr[size_t(i)]; p < M.rowPtr[size_t(i) + 1]; ++p)
            if (M.col[size_t(p)] != i && C.colour[size_t(i)] == C.colour[size_t(M.col[size_t(p)])]) valid = false;
        // equivalence: MC-BILU0 = natural BILU0 of the permuted matrix
        const BSR Pm = permuteBsr(M, C);
        BILU0 Pp; Pp.factor(Pm);
        std::vector<double> r(static_cast<size_t>(N)), z1(static_cast<size_t>(N)), rp(static_cast<size_t>(N)), zp(static_cast<size_t>(N));
        for (int e = 0; e < N; ++e) r[size_t(e)] = std::cos(0.7 * e);
        P1.apply(r.data(), z1.data());
        for (int k = 0; k < M.n; ++k) for (int e = 0; e < 3; ++e) rp[size_t(k) * 3 + e] = r[size_t(C.perm[size_t(k)]) * 3 + e];
        Pp.apply(rp.data(), zp.data());
        double de = 0, zn = 0;
        for (int k = 0; k < M.n; ++k) for (int e = 0; e < 3; ++e) {
            de = std::max(de, std::fabs(zp[size_t(k) * 3 + e] - z1[size_t(C.perm[size_t(k)]) * 3 + e])); zn = std::max(zn, std::fabs(zp[size_t(k) * 3 + e]));
        }
        std::vector<double> x0, x1;
        auto t0 = std::chrono::steady_clock::now();
        const LinResult r0 = bicgstab(op, P0, b, x0, 1e-9, 3000);
        auto t1 = std::chrono::steady_clock::now();
        const LinResult r1 = bicgstab(op, P1, b, x1, 1e-9, 3000);
        auto t2 = std::chrono::steady_clock::now();
        double e1 = 0, xn = 0; for (int e = 0; e < N; ++e) { e1 = std::max(e1, std::fabs(x1[size_t(e)] - xt[size_t(e)])); xn = std::max(xn, std::fabs(xt[size_t(e)])); }
        CpuFloatBlockSolver dev; dev.multicolour = true;
        std::vector<double> xr;
        const LinResult rr = refineWithDevice(dev, M, b, xr, 1e-10, 8, 3000);
        double er = 0; for (int e = 0; e < N; ++e) er = std::max(er, std::fabs(xr[size_t(e)] - xt[size_t(e)]));
        std::printf("  %s: %d colours (sizes", v ? "Laplace + 48 long couplings" : "5-point 120x100 nb3", C.colours);
        for (int c = 0; c < std::min(C.colours, 6); ++c) std::printf(" %d", C.start[size_t(c) + 1] - C.start[size_t(c)]);
        std::printf(") | BILU0 %d its %.1f ms, MC-BILU0 %d its %.1f ms, err %.1e | equivalence %.1e | float32 device emulation + refinement: err %.1e (%d its)\n",
                    r0.iterations, 1e3 * std::chrono::duration<double>(t1 - t0).count(), r1.iterations, 1e3 * std::chrono::duration<double>(t2 - t1).count(),
                    e1 / xn, de / zn, er / xn, rr.iterations);
        check(valid && (v || C.colours == 2), v ? "greedy colouring of an irregular block graph is a valid distance-1 colouring"
                                                : "5-point stencil gets the red-black (2-colour) ordering");
        check(de / zn < 1e-12, "MC-BILU0 (per-colour parallel factor + sweeps) = natural BILU0 of P A P^T");
        check(r1.converged && e1 / xn < 1e-6 && r1.iterations <= 3 * r0.iterations + 5, "MC-BILU0 BiCGSTAB converges (iteration penalty of the colour ordering bounded)");
        check(rr.converged && er / xn < 1e-8, "float32 MC-BILU0 device emulation + double refinement reaches double accuracy");
    }
}

// C18: MPI domain-decomposed linear solve. In the normal suite (one rank) it checks the partition and the algorithm;
// under "mpiexec -n P run_tests -mpidd" rank 0 drives a real multi-rank solve + a full simulation while the other ranks
// serve (serveDD).
static void testMpiDD(bool underMpi) {
    std::printf("== C18 MPI domain decomposition (distributed BiCGSTAB, halo exchange, block-Jacobi BILU0)\n");
    const BSR A = makeBlockLaplace(120, 100, 3, 13);
    {
        const auto st = mpi::partitionRows(A, 4);
        bool mono = st.front() == 0 && st.back() == A.n;
        int mx = 0, mn = 1 << 30;
        for (int q = 0; q < 4; ++q) {
            mono = mono && st[size_t(q) + 1] >= st[size_t(q)];
            const int nz = A.rowPtr[size_t(st[size_t(q) + 1])] - A.rowPtr[size_t(st[size_t(q)])];
            mx = std::max(mx, nz); mn = std::min(mn, nz);
        }
        check(mono && mx < 1.02 * mn + 20, "nnz-balanced contiguous partition (4 parts)");
    }
    const int N = A.size();
    std::vector<double> xt(static_cast<size_t>(N)), b(static_cast<size_t>(N));
    for (int e = 0; e < N; ++e) xt[size_t(e)] = std::sin(0.002 * e) + 0.3 * (e % 3);
    A.spmv(xt.data(), b.data());
    mpi::DDSolver dd;
    std::vector<double> x;
    int its = 0, itsBJ = 0;
    dd.cpr = false;
    std::vector<double> xbj;
    const bool okBJ = dd.solve(A, b, xbj, 1e-10, 3000, itsBJ);
    dd.cpr = true;
    const auto t0 = std::chrono::steady_clock::now();
    const bool ok = dd.solve(A, b, x, 1e-10, 3000, its);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    double e1 = 0, xn = 0; for (int e = 0; e < N; ++e) { e1 = std::max(e1, std::fabs(x[size_t(e)] - xt[size_t(e)])); xn = std::max(xn, std::fabs(xt[size_t(e)])); }
    BILU0 M; M.factor(A); BsrOperator op(A); std::vector<double> xs;
    const LinResult rs = bicgstab(op, M, b, xs, 1e-10, 3000);
    std::printf("  %d ranks: CPR (global AMG pressure + block-Jacobi BILU0) %d its, %.3f s, err %.1e | block-Jacobi BILU0 only %d its | serial BILU0: %d its | halo %lld doubles, %.3f s\n",
                dd.stats.ranks, its, sec, e1 / xn, itsBJ, rs.iterations, dd.stats.haloDoubles, dd.stats.haloSec);
    check(okBJ && ok && its < itsBJ, "CPR two-stage preconditioner needs fewer iterations than block-Jacobi BILU0");
    check(ok && e1 / xn < 1e-7, underMpi ? "distributed solve over all ranks reaches the serial solution" : "single-rank DD solve = serial solution");
    if (dd.stats.ranks == 1) check(std::abs(itsBJ - rs.iterations) <= 1, "one rank: block-Jacobi = global BILU0 (same iterations)");
    // second solve with the same pattern: no pattern re-scatter
    std::vector<double> b2(b); for (auto& v : b2) v *= 2;
    dd.solve(A, b2, x, 1e-10, 3000, its);
    check(dd.stats.patternSetups == 1, "pattern scattered once, values per solve");
    // full black-oil run with the distributed solver hooked into the Newton loop
    Case c = makeDemoCase();
    c.schedule.resize(6);
    double fo[2]; int nw[2]; RunStats st[2];
    for (int k = 0; k < 2; ++k) {
        setDistributedSolver(k ? &dd : nullptr);
        Simulator sim; sim.setup(c); sim.initialize();
        while (sim.runReportStep()) {}
        fo[k] = sim.summary().back().fopt; nw[k] = sim.stats().newton; st[k] = sim.stats();
    }
    setDistributedSolver(nullptr);
    std::printf("  demo 6 steps: FOPT serial %.8g / DD %.8g, Newton %d / %d, DD solves %d (fallbacks %d), %.3f s in DD solves\n", fo[0], fo[1], nw[0], nw[1],
                st[1].ddSolves, st[1].ddFallbacks, st[1].ddSec);
    check(st[1].ddSolves > 0 && st[1].ddFallbacks == 0 && std::fabs(fo[1] - fo[0]) < 1e-5 * fo[0], "simulation with the domain-decomposed solver = serial run");
    if (underMpi) dd.shutdown();
}

// ---- C18 in-run hydraulic-fracture growth (LEFM tip criterion inside the flow run)
static void testHfGrowth() {
    std::printf("== C18 in-run hydraulic-fracture growth (K_I = p_net sqrt(pi min(L, H/2)) >= K_IC opens the next piece)\n");
    // water injection into the tight-oil EDFM case: the 200 m trace is the maximum extent, 25 m wings open at the start
    auto deck = [](const std::string& grow, const std::string& rate) {
        std::string d = hfDeck(true, false, grow);
        const std::string a = "WCONPROD\n P OPEN ORAT 5 4* 10 /";
        const size_t at = d.find(a);
        if (at != std::string::npos) d.replace(at, a.size(), "WCONINJE\n P WATER OPEN RATE " + rate + " 1* 700 /");
        const size_t q = d.find(" 0.2 0 1 0\n");   // slightly mobile connate water (krw > 0 at the initial saturation) for the injector
        if (q != std::string::npos) d.replace(q, 11, " 0.2 0 1 0\n 0.3 0.02 0.7 0\n");
        return d;
    };
    struct R { bool ok = false; double bhp = 0, wl[2] = {0, 0}, wi = 0, t75 = 1e9, pMin = 1e9, pMax = -1e9; int ev = 0, hf = 0; std::string err; };
    auto run = [](const std::string& d, Simulator* keep = nullptr) {
        R r;
        try {
            Case c = eclipseDeckFromString(d, "");
            Simulator local; Simulator& s = keep ? *keep : local;
            s.log = [](const std::string&) {};
            s.setup(c); s.initialize();
            while (s.runReportStep()) {
                // injection pressure while the fracture is still growing (between the first and the last opening)
                const bool growing = !s.hfGrowthEvents().empty() && s.hfWingLength(0, 1) < 99.9;
                if (growing && !s.summary().back().wells.empty()) {
                    r.pMin = std::min(r.pMin, s.summary().back().wells[0].bhp); r.pMax = std::max(r.pMax, s.summary().back().wells[0].bhp);
                }
            }
            for (const auto& e : s.hfGrowthEvents()) if (e.side > 0 && e.length >= 75 - 1e-9 && e.t < r.t75) r.t75 = e.t;
            r.bhp = s.summary().back().wells.empty() ? 0 : s.summary().back().wells[0].bhp;
            r.wl[0] = s.hfWingLength(0, -1); r.wl[1] = s.hfWingLength(0, 1);
            r.ev = int(s.hfGrowthEvents().size()); r.hf = s.numHydraulicFractureCells();
            r.wi = s.summary().back().fwit;
            r.ok = true;
        } catch (const std::exception& e) { r.err = e.what(); std::printf("    exception: %s\n", e.what()); }
        return r;
    };
    const std::string g1 = "-- RESSIM: HFGROW F1 1.0 320 0 20\n", gHi = "-- RESSIM: HFGROW F1 1.0 340 0 20\n", gClosed = "-- RESSIM: HFGROW F1 1e6 320 0 20\n";
    const R full = run(deck("", "2")), grow = run(deck(g1, "2")), hi = run(deck(gHi, "2")), closed = run(deck(gClosed, "2"));
    const double pCrit = 1.0e6 / std::sqrt(M_PI * 5.0) / 1e5;   // K_IC / sqrt(pi H/2) in bar (H = 10 m: PKN regime)
    std::printf("    final BHP: fully open trace %.2f bar | growth sigma 320: wings %.1f / %.1f m, %d events, 75 m at %.2f d, BHP %.2f | sigma 340: wings %.1f / %.1f m,"
                " 75 m at %.2f d, BHP while growing %.2f - %.2f (sigma + p_crit = %.2f) | no growth (K_IC 1e6): wings %.1f / %.1f m, BHP %.2f\n",
                full.bhp, grow.wl[0], grow.wl[1], grow.ev, grow.t75, grow.bhp, hi.wl[0], hi.wl[1], hi.t75, hi.pMin, hi.pMax, 340 + pCrit, closed.wl[0], closed.wl[1], closed.bhp);
    check(grow.ok && closed.ok && std::fabs(closed.wl[0] - 25) < 1.5 && std::fabs(closed.wl[1] - 25) < 1.5 && closed.ev == 0,
          "HF growth: initial opening = pieces within initLen of the well (25 m wings), none grow when K_IC is unreachable");
    check(grow.ok && grow.ev > 0 && grow.wl[0] > 35 && grow.wl[1] > 35, "HF growth: pressurised fracture propagates (K_I >= K_IC)");
    check(grow.ok && std::fabs(grow.wl[0] - grow.wl[1]) <= 10.5, "HF growth: symmetric wings around a centred injector (within one piece)");
    check(hi.ok && hi.t75 > 2 * grow.t75 && hi.t75 < 1e8, "HF growth: higher closure stress -> slower propagation");
    check(hi.ok && hi.pMin > 340 - 2 * pCrit && hi.pMax < 340 + 2 * pCrit,
          "HF growth: injection pressure pinned within +-2 K_IC / sqrt(pi H / 2) of the closure stress while the fracture propagates");
    check(full.ok && grow.bhp <= closed.bhp + 1e-6 && grow.bhp >= full.bhp - 1e-6, "HF growth: injection pressure between the fully open and the non-growing fracture");
    // deck round trip and restart of the growth state
    try {
        Case c = eclipseDeckFromString(deck(g1, "2"), "");
        const Case c2 = eclipseDeckFromString(eclipseDeckToString(c), "");
        const auto& F = c2.grid.hfracs.at(0);
        check(std::fabs(F.kic - 1.0) < 1e-12 && std::fabs(F.sigma - 320) < 1e-9 && std::fabs(F.initLen - 20) < 1e-9 && std::fabs(F.closedMult - 1) < 1e-12,
              "HFGROW: deck export / import round trip");
        c = eclipseDeckFromString(deck(gHi, "2"), "");   // restart in the middle of the propagation (sigma 340: 14 - 66 d)
        Simulator a; a.log = [](const std::string&) {};
        a.setup(c); a.initialize();
        const int half = int(c.schedule.size()) * 3 / 5;
        for (int k = 0; k < half; ++k) a.runReportStep();
        const double wMid = a.hfWingLength(0, 1);
        const std::string blob = a.saveState();
        while (a.runReportStep()) {}
        Simulator b; b.log = [](const std::string&) {};
        b.setup(c); b.initialize(); b.loadState(blob);
        while (b.runReportStep()) {}
        const double d = std::fabs(a.summary().back().wells[0].bhp - b.summary().back().wells[0].bhp);
        std::printf("    restart at step %d (wing %.0f m): wing lengths %.1f / %.1f vs %.1f / %.1f, BHP difference %.2e bar\n", half, wMid, a.hfWingLength(0, -1), a.hfWingLength(0, 1),
                    b.hfWingLength(0, -1), b.hfWingLength(0, 1), d);
        check(wMid > 25 && wMid < 100 && d < 1e-6 && a.hfWingLength(0, 1) == b.hfWingLength(0, 1) && a.hfGrowthEvents().size() == b.hfGrowthEvents().size(), "HF growth: restart reproduces the grown fracture and the injection pressure");
        bool snap = false;
        for (const auto& e : a.snapshots().back().extra) if (e.first == "HFOPEN") snap = true;
        check(snap, "HF growth: HFOPEN snapshot array for the 3-D view");
    } catch (const std::exception& e) { std::printf("    exception: %s\n", e.what()); check(false, "HFGROW round trip / restart"); }
}
