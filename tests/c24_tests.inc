// =====================================================================================================================
// C24 tests: fault mechanics + induced seismicity inside the flow run (faultmech.h), in-run G1-c
// =====================================================================================================================
#include "ressim/faultmech.h"

static std::string findExample(const char* name) {
    for (const char* pre : {"examples/", "../examples/", "../../examples/"}) {
        const std::string p = std::string(pre) + name;
        if (std::ifstream(p).good()) return p;
    }
    return "";
}

static void testFaultMechStress() {
    std::printf("== C24 fault mechanics: geometry, stress resolution, Dieterich rate vs closed form\n");
    const std::string path = findExample("FAULT_INJECTION.DATA");
    if (path.empty()) { std::printf("    (examples/FAULT_INJECTION.DATA not found - skipped)\n"); return; }
    const Case c = importEclipseDeck(path);
    Grid g; g.build(c);
    faultmech::Monitor M; std::string err;
    const bool ok = M.setup(c, g, err);
    if (!ok) std::printf("    setup error: %s\n", err.c_str());
    check(ok, "FAULTMECH resolves the FAULTS faces");
    if (!ok) return;
    const auto& F = M.out().faces;
    size_t withConn = 0; for (const auto& f : F) withConn += f.conns.empty() ? 0 : 1;
    std::printf("    %zu faces, %zu with connections across, normal of face 0 (%.3f %.3f %.3f), area %.1f m2\n", F.size(), withConn, F[0].n[0], F[0].n[1], F[0].n[2], F[0].area);
    check(F.size() == 100 && withConn == F.size(), "20 x 5 faces of fault F1, each with its cross-fault connection");
    check(std::fabs(std::fabs(F[0].n[0]) - 1) < 1e-6 && std::fabs(F[0].area - 500) < 1e-6, "x-face normal and area 50 x 10 m");
    // t0 with a uniform 200 bar: hand calculation at the face depth
    std::vector<float> p(size_t(g.numActive()), 200.0f);
    M.step(0, p);
    const auto& f0 = M.out().faces[0];
    const double z = f0.c[2], Sv = 0.24 * z, SH = 1.10 * Sv, Sh = 0.65 * Sv, th = 60 * M_PI / 180;
    const double snH = SH * std::cos(th) * std::cos(th) + Sh * std::sin(th) * std::sin(th), tauH = (SH - Sh) * std::sin(th) * std::cos(th);
    const double cffH = tauH - 0.6 * (snH - 200);
    std::printf("    z %.1f m: sigma_n %.3f (hand %.3f) bar, tau %.3f (hand %.3f), CFF %.3f (hand %.3f)\n", z, f0.sn0, snH, f0.tau0, tauH, M.out().steps[0].cff[0], cffH);
    check(std::fabs(f0.sn0 - snH) < 1e-3 && std::fabs(f0.tau0 - tauH) < 1e-3 && std::fabs(M.out().steps[0].cff[0] - cffH) < 1e-3, "stress tensor resolved on the face = hand calculation");
    // constant pressurisation rate -> constant Coulomb stressing rate: Dieterich closed form N(t)
    const double rate = 0.05;   // bar/day
    const double days = 400;
    for (int s = 1; s <= 80; ++s) { std::fill(p.begin(), p.end(), float(200 + rate * days * s / 80)); M.step(days * s / 80, p); }
    const auto& S1 = M.out().steps[1];
    const double dCdp = S1.dcfs[0] / (rate * days / 80);
    const double aSig = f0.aSig, ta = 100, Sr = aSig / ta, Sdot = Sr + dCdp * rate * 365.25;   // per year
    const double r = c.faultMech.faults[0].r0 * f0.area * 1e-6, k = Sdot / Sr, tau = aSig / Sdot, t = days / 365.25;
    const double Nexact = r * k * (t + tau * std::log((1 + (k - 1) * std::exp(-t / tau)) / k));
    const double Nnum = M.out().steps.back().cum[0];
    std::printf("    dCFF/dp %.3f, A sigma %.3f bar, S_dot / S_r %.0f: expected events on face 0 %.6g (closed form %.6g)\n", dCdp, aSig, k, Nnum, Nexact);
    check(std::fabs(dCdp - (0.6 - 0.5 * 0.6)) < 1e-3, "dCFF/dp = mu (1 - A) for this stress path (tau independent of the isotropic horizontal change)");
    check(std::fabs(Nnum / Nexact - 1) < 0.01, "Dieterich expected count within 1 % of the closed form for a constant stressing rate");
}

static double runFaultCase(const Case& c, faultmech::Output& out, double& dpAcross) {
    Simulator sim; sim.setup(c); sim.initialize();
    faultmech::Monitor M; std::string err;
    if (!M.setup(c, sim.grid(), err)) { std::printf("    setup: %s\n", err.c_str()); return -1; }
    M.step(sim.time(), sim.snapshots().back().pressure);
    const auto t0 = std::chrono::steady_clock::now();
    double dpMaxInj = 0;
    const Grid& G = sim.grid();
    const int aW = G.activeOf[size_t(14 + G.nx * (9 + G.ny * 2))], aE = G.activeOf[size_t(15 + G.nx * (9 + G.ny * 2))];
    while (sim.runReportStep()) {
        const auto& P = sim.snapshots().back().pressure;
        M.step(sim.time(), P);
        M.applyCoupling(sim);
        if (sim.time() <= 500 + 1e-9) dpMaxInj = std::max(dpMaxInj, double(P[size_t(aW)] - P[size_t(aE)]));
    }
    M.finish(nullptr, [](const std::string& s) { std::printf("%s\n", s.c_str()); });
    out = M.out();
    dpAcross = dpMaxInj;
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

static void testFaultMechRun() {
    std::printf("== C24 fault mechanics in the flow run: injection next to a critically stressed fault (two-way vs one-way)\n");
    const std::string path = findExample("FAULT_INJECTION.DATA");
    if (path.empty()) { std::printf("    (examples/FAULT_INJECTION.DATA not found - skipped)\n"); return; }
    Case c = importEclipseDeck(path);
    faultmech::Output two, one; double dp2 = 0, dp1 = 0;
    const double s2 = runFaultCase(c, two, dp2);
    Case c1 = c; for (auto& r : c1.faultMech.faults) { r.permFactor = 1; r.rsfLayer = 0; }
    const double s1 = runFaultCase(c1, one, dp1);
    check(s2 >= 0 && s1 >= 0 && two.steps.size() > 10, "coupled runs complete");
    if (two.steps.size() < 2 || one.steps.size() < 2) return;
    const auto& R2 = two.series[0]; const auto& R1 = one.series[0];
    size_t iMax = 0; for (size_t i = 0; i < R2.t.size(); ++i) if (R2.reactFrac[i] > R2.reactFrac[iMax]) iMax = i;
    int nD = 0, nR = 0; double mwMax = -9;
    for (const auto& e : two.events) { (e.kind == 0 ? nD : nR)++; mwMax = std::max(mwMax, double(e.mw)); }
    std::printf("    two-way: max dCFS %.2f bar, reactivated %.0f %% of F1 at t = %.0f d, expected events %.3g (one-way %.3g), catalogue %d Dieterich + %d RSF, Mw max %.2f, M_cap %.2f, %ld coupling updates (%.1f s)\n",
                *std::max_element(R2.maxDcfs.begin(), R2.maxDcfs.end()), 100 * R2.reactFrac[iMax], R2.t[iMax], R2.expected.back(), R1.expected.back(), nD, nR, mwMax, two.mMax[0], two.twoWayUpdates, s2);
    std::printf("    reactivated fraction of F1 over time:");
    for (size_t i = 0; i < R2.t.size(); i += 2) std::printf(" %.0fd:%.0f%%", R2.t[i], 100 * R2.reactFrac[i]);
    std::printf("\n    pressure jump across F1 at the injector row during injection: two-way %.2f bar, one-way %.2f bar; critical faces at t0: %d\n", dp2, dp1, two.criticalAtStart);
    const double bg = c.faultMech.faults[0].r0 * 20 * 5 * 500 * 1e-6 * (R2.t.back() / 365.25);
    check(*std::max_element(R2.maxDcfs.begin(), R2.maxDcfs.end()) > 5 && R2.reactFrac[iMax] > 0, "injection raises the Coulomb stress and reactivates part of the fault");
    check(R2.expected.back() > 10 * bg, "expected induced events >> background over the same time");
    check(nD > 0 && !two.events.empty(), "stochastic catalogue realisation is not empty");
    check(dp2 < dp1, "two-way coupling: reactivated faces conduct, the pressure jump across the fault shrinks");
    check(R2.expected.back() <= R1.expected.back() * 1.5, "two-way coupling relieves the injection-side pressure (no runaway in the count)");
    check(!two.rsf.empty() && two.rsf[0].error.empty() && two.rsf[0].Lb / std::max(1e-9, two.rsf[0].dx) > 1, "rate-and-state trace on layer 3 runs on the simulated forcing");
    bool mono = true; for (size_t i = 1; i < R2.expected.size(); ++i) mono = mono && R2.expected[i] >= R2.expected[i - 1] - 1e-12;
    check(mono, "expected event count is non-decreasing in time");
    // deck round trip of the new keywords
    const Case rt = projectFromText(projectText(c));
    check(rt.faultMech.on() && rt.faultMech.faults.size() == 1 && rt.faultMech.faults[0].permFactor == 20 && rt.faultMech.faults[0].rsfLayer == 3 &&
          std::fabs(rt.faultMech.sHRatio - 1.10) < 1e-12, "STRESSFIELD / FAULTMECH / FAULTRSF survive the deck round trip");
}
