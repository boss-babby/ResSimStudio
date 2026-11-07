// =====================================================================================================================
// C23 tests: G3 rate-and-state friction on embedded faults (quasi-dynamic cycles, injection, EFEM / G1-c coupling)
// =====================================================================================================================
#include "ressim/rsf.h"
#include "ressim/efem.h"
#include "ressim/geocouple.h"

static void testRsfSpringSlider() {
    std::printf("== C23 G3 rate-and-state: one element = spring-slider, stick-slip below the critical stiffness only\n");
    const double yr = 3.15576e7;
    for (double kr : {0.5, 2.0}) {
        rsf::Input in; in.n = 1;
        // self stiffness G / (pi (1 - nu) L / 2) = kr * k_c with k_c = (b - a) sigma / Dc = 2.5e8 Pa/m
        in.x1 = in.G / (M_PI * (1 - in.nu) * 0.5 * kr * 2.5e8);
        in.vInit = 1.01e-9; in.tEnd = 60 * yr;
        rsf::Output o; const bool ok = rsf::run(in, o);
        std::printf("    k/k_c %.1f (k %.3g, k_c %.3g Pa/m): %zu events, %ld steps (%ld rejected), %.2f s\n", kr, o.kSelf, o.kc, o.events.size(), o.steps, o.rejected, o.seconds);
        check(ok && std::fabs(o.kSelf / o.kc - kr) < 1e-6, "element self-stiffness from the dislocation kernel equals the design value");
        if (kr < 1) {
            check(o.events.size() >= 10, "k < k_c: periodic stick-slip events");
            if (o.events.size() >= 6) {
                const double T1 = o.events[4].tStart - o.events[3].tStart, T2 = o.events[5].tStart - o.events[4].tStart;
                double slipCycle = 0;   // seismic + aseismic slip per cycle must equal Vpl T (slip budget)
                std::printf("    recurrence %.4f / %.4f yr, coseismic slip %.4f m, Vpl T %.4f m, Mw %.2f\n", T1 / yr, T2 / yr, o.events[4].slipMax, in.Vpl * T2, o.events[4].Mw);
                slipCycle = in.Vpl * T2;
                check(std::fabs(T1 - T2) < 0.005 * T2, "limit cycle: constant recurrence interval");
                check(o.events[4].slipMax < slipCycle && o.events[4].slipMax > 0.5 * slipCycle, "coseismic slip is the dominant part of the Vpl T slip budget");
            }
        } else check(ok && o.events.empty(), "k > k_c: stable sliding, no events");
    }
}

static void testRsfFaultCycle() {
    std::printf("== C23 G3 rate-and-state: finite fault with a velocity-weakening patch (back-slip loading, radiation damping)\n");
    const double yr = 3.15576e7;
    rsf::Input in; in.n = 128; in.tEnd = 60 * yr;
    rsf::Output o; const bool ok = rsf::run(in, o);
    std::printf("    h* %.0f m, Lb %.1f m, dx %.1f m, %zu events in 60 yr, %ld steps, %.2f s\n", o.hStar, o.Lb, o.dx, o.events.size(), o.steps, o.seconds);
    check(ok && o.Lb / o.dx > 3, "cohesive zone resolved (Lb / dx > 3)");
    check(o.events.size() >= 8, "repeating events on a patch of ~4 nucleation lengths");
    if (o.events.size() >= 6) {
        const auto& E = o.events;
        const size_t k = E.size() - 1;
        const double T1 = E[k - 1].tStart - E[k - 2].tStart, T2 = E[k].tStart - E[k - 1].tStart;
        std::printf("    last events: T %.3f / %.3f yr, slip %.3f m, rupture %.0f m, stress drop %.2f MPa, Mw %.2f\n", T1 / yr, T2 / yr, E[k].slipMax, E[k].length, E[k].stressDrop * 1e-6, E[k].Mw);
        check(std::fabs(T1 - T2) < 0.01 * T2, "periodic system-size events");
        check(E[k].slipMax > 0.9 * in.Vpl * T2 && E[k].slipMax < 1.1 * in.Vpl * T2, "slip budget: coseismic slip at the patch = Vpl x recurrence within 10 %");
        check(E[k].stressDrop > 0.5e6 && E[k].stressDrop < 20e6 && E[k].Mw > 3 && E[k].Mw < 6, "stress drop and magnitude in the physical range");
        check(E[k].vMax > 0.05 && E[k].vMax < 10, "peak slip rate limited by radiation damping (quasi-dynamic, ~m/s)");
        // time-step convergence
        rsf::Input i2 = in; i2.rtol = 1e-8; i2.tEnd = 25 * yr;
        rsf::Output o2; rsf::run(i2, o2);
        if (o2.events.size() >= 2 && E.size() >= 2) {
            const double Ta = E[1].tStart - E[0].tStart, Tb = o2.events[1].tStart - o2.events[0].tStart;
            std::printf("    rtol 1e-7 vs 1e-8: first recurrence %.4f vs %.4f yr\n", Ta / yr, Tb / yr);
            check(std::fabs(Ta - Tb) < 0.005 * Tb, "tolerance-converged recurrence (0.5 %)");
        } else check(false, "rtol 1e-8 run produced events");
    }
    // steady sliding: velocity-strengthening fault started at Vpl stays at Vpl
    rsf::Input vs = in; vs.aVW = vs.aVS = 0.02; vs.tEnd = 5 * yr;
    rsf::Output ov; rsf::run(vs, ov);
    double dv = 0; for (double v : ov.finalV) dv = std::max(dv, std::fabs(v / vs.Vpl - 1));
    std::printf("    steady sliding check: max |V/Vpl - 1| = %.2e after 5 yr, %zu events\n", dv, ov.events.size());
    check(ov.ok && ov.events.empty() && dv < 1e-5, "steady state at the loading rate is preserved (kernel rows of a uniform slip sum consistently)");
}

static void testRsfInjection() {
    std::printf("== C23 G3 rate-and-state: injection-induced nucleation (2-D Theis pore pressure on the fault)\n");
    const double yr = 3.15576e7, day = 86400;
    rsf::Input base; base.n = 96; base.x1 = 1500; base.Vpl = 1e-10; base.tEnd = 40 * yr;
    base.vwFrom = 0.25; base.vwTo = 0.75;
    rsf::Output o0; rsf::run(base, o0);
    rsf::Input inj = base; inj.inj.on = true; inj.inj.at = 0.5; inj.inj.offset = 20; inj.inj.pStar = 1.5e6; inj.inj.diffusivity = 0.05;
    inj.inj.tStart = 0; inj.inj.tStop = 200 * day;
    rsf::Output o1; rsf::run(inj, o1);
    const double t0 = o0.events.empty() ? 1e300 : o0.events[0].tStart, t1 = o1.events.empty() ? 1e300 : o1.events[0].tStart;
    std::printf("    first event: tectonic only %.3f yr, with 200 d injection %.4f yr (%.0f d), Mw %.2f, max dp %.2f MPa\n", t0 / yr, t1 / yr, t1 / day,
                o1.events.empty() ? 0.0 : o1.events[0].Mw, o1.pMax.empty() ? 0.0 : *std::max_element(o1.pMax.begin(), o1.pMax.end()) * 1e-6);
    check(o1.ok && t1 < 0.1 * t0, "pore pressure on a velocity-weakening patch advances the first event (clock advance)");
    // velocity-strengthening fault: aseismic slip, no events
    rsf::Input sv = inj; sv.aVW = sv.aVS = 0.022; sv.tEnd = 2 * yr;
    rsf::Output o2; rsf::run(sv, o2);
    double smax = 0; for (double s : o2.finalSlip) smax = std::max(smax, s);
    double vm = 0; for (double v : o2.vmax) vm = std::max(vm, v);
    std::printf("    velocity strengthening + injection: %zu events, max slip %.3g mm, max V %.2e m/s\n", o2.events.size(), smax * 1e3, vm);
    check(o2.ok && o2.events.empty() && smax > 10 * sv.Vpl * sv.tEnd && vm < 1e-3, "velocity-strengthening fault responds to injection by aseismic creep only");
}

static void testRsfEmbedded() {
    std::printf("== C23 G3 rate-and-state on an EFEM-embedded fault: initial stresses from the contact solve + G1-c pressure table\n");
    efem::Input ei; ei.Lx = ei.Ly = 10; ei.nx = ei.ny = 60; ei.E = 20e9; ei.nu = 0.25; ei.sxx = -20e6; ei.syy = -40e6;
    efem::Fracture F; const double th = 40 * M_PI / 180;   // fault dip from the x axis
    F.x0 = 5 - 2 * std::cos(th); F.y0 = 5 - 2 * std::sin(th); F.x1 = 5 + 2 * std::cos(th); F.y1 = 5 + 2 * std::sin(th); F.friction = 0.9;
    ei.fracs.push_back(F);
    efem::Result er; const bool eok = efem::solve(ei, er);
    rsf::Input in; in.n = 40;
    std::string err;
    const bool ok = eok && rsf::fromEfem(ei, er, 0, in, err);
    // analytic resolved stresses for the locked fault (friction 0.9 > tau/sigma): n = (-sin, cos)
    const double sn = -(ei.sxx * std::sin(th) * std::sin(th) + ei.syy * std::cos(th) * std::cos(th));
    const double tau = std::fabs((ei.syy - ei.sxx) * std::sin(th) * std::cos(th));
    const double snC = in.sigmaNProf.empty() ? 0 : in.sigmaNProf[20], taC = in.tau0Prof.empty() ? 0 : in.tau0Prof[20];
    std::printf("    centre: sigma_n %.2f MPa (analytic %.2f), tau %.2f MPa (analytic %.2f)%s\n", snC * 1e-6, sn * 1e-6, taC * 1e-6, tau * 1e-6, err.empty() ? "" : (" " + err).c_str());
    check(ok && std::fabs(snC - sn) < 0.05 * sn && std::fabs(taC - tau) < 0.05 * tau, "fromEfem: fault stresses from the embedded contact solve match the resolved far field (5 %)");
    // G1-c pressure table -> forcing
    geocouple::Output co; co.ok = true;
    for (int k = 0; k < 5; ++k) { co.s.push_back(-2 + k); co.frac.push_back(0); }
    for (int st = 0; st < 3; ++st) { geocouple::StepOut so; so.t = st * 10.0; for (int k = 0; k < 5; ++k) so.p.push_back(10e6 + st * 1e6 * (k == 2 ? 2 : 1)); co.steps.push_back(so); }
    const bool ok2 = rsf::fromCoupled(co, 0, in, err);
    std::printf("    coupled table: %zu times, dp centre at the last step %.2f MPa, end %.2f MPa\n", in.forceT.size(), in.forceP.empty() ? 0.0 : in.forceP.back()[20] * 1e-6, in.forceP.empty() ? 0.0 : in.forceP.back()[0] * 1e-6);
    check(ok2 && in.forceT.size() == 3 && std::fabs(in.forceT[2] - 20 * 86400.0) < 1e-6 && std::fabs(in.forceP[2][0] - 2e6) < 1e3 && in.forceP[2][20] > 2e6 && in.forceP[0][20] == 0,
          "fromCoupled: piece pressures relative to the first step, mapped along the fault, days -> s");
    in.tEnd = 30 * 86400.0; in.Vpl = 0;
    rsf::Output o; const bool ok3 = rsf::run(in, o);
    std::printf("    embedded-fault run: %ld steps, %zu events, max slip %.3g mm\n", o.steps, o.events.size(), o.finalSlip.empty() ? 0.0 : *std::max_element(o.finalSlip.begin(), o.finalSlip.end()) * 1e3);
    check(ok3, "rate-and-state run on the embedded fault with the coupled pressure table completes");
}

static void testKernelFitAspo() {
    std::printf("== C23 research Cycle 5: Aspo stage-count kernel fit (seis::fitStageKernels)\n");
    // stage table of HF1 / HF2 (research/exp/h6_aspo.py; V litres, E kJ, mean rate L/min, counts co+post)
    struct R { int test; double V, E, q; int n; };
    const R rows[] = {{0, 0.77, 7.9, 0.96, 0}, {0, 2.66, 24.3, 1.29, 3}, {0, 4.40, 35.8, 1.58, 10}, {0, 5.14, 42.2, 2.41, 12}, {0, 4.67, 41.1, 4.12, 3}, {0, 9.25, 81.5, 4.44, 20},
                      {1, 2.96, 26.7, 1.67, 4}, {1, 5.03, 44.4, 1.63, 6}, {1, 1.22, 10.1, 2.36, 0}, {1, 5.19, 47.8, 2.45, 17}, {1, 4.84, 46.8, 4.61, 15}, {1, 9.42, 101.5, 4.75, 50}};
    std::vector<seis::StageObs> S;
    for (const R& r : rows) { seis::StageObs o; o.test = r.test; o.volume = r.V; o.energy = r.E; o.rate = r.q; o.count = r.n; S.push_back(o); }
    const seis::KernelFit f = seis::fitStageKernels(S);
    for (const auto& m : f.models) std::printf("    %-34s logL %8.2f  AIC %7.2f%s\n", m.name.c_str(), m.logL, m.aic, m.exponent != 0 ? ("  exponent " + std::to_string(m.exponent)).c_str() : "");
    std::printf("    best: %s; rate exponent at fixed volume %.2f [%.2f, %.2f]; verdict: %s\n", f.best.c_str(), f.rateExponent, f.rateLo, f.rateHi, f.verdict.c_str());
    check(!f.models.empty() && f.best.find("volume") != std::string::npos, "volume-type kernels beat the constant / duration kernels on stage counts");
    check(f.rateLo <= f.rateExponent && f.rateExponent <= f.rateHi, "profile-likelihood interval brackets the rate exponent");
    check(f.rateLo < 0 && f.rateHi > 0 && f.cumExponent > 0, "Aspo (research F14): no resolvable rate effect once the growing fracture (cumulative volume) is accounted for");
}

#include "ressim/phasefield.h"
static void testPhaseFieldSneddonGriffith() {
    std::printf("== C23 G4: phase-field (AT2) pressurised crack - Sneddon volume and Griffith critical pressure, l -> 0\n");
    // Bourdin-Chukwudozie-Yoshioka benchmark: a = 0.1 m, E 10 GPa, nu 0.25, Gc 100 J/m2, l = 2h. Two resolutions (l/a = 0.2, 0.1);
    // diffuse models converge like O(l/a), so the check is on the l -> 0 extrapolation (measured C23: 0.96 and 0.975)
    double sr[2] = {0, 0}, pr[2] = {0, 0}; bool ok = true;
    const double hs[2] = {0.01, 0.005};
    for (int i = 0; i < 2; ++i) {
        phasefield::Input in; in.h = hs[i];
        phasefield::Output o;
        const bool r = phasefield::run(in, o);
        ok = ok && r && o.pcNumeric > 0;
        sr[i] = o.sneddonRatio; pr[i] = o.pcNumeric / o.pcTheory;
        std::printf("    h %.3f l %.3f (%d nodes, %.0f s): Sneddon ratio %.3f, p_c %.3f MPa in [%.3f, %.3f] vs Griffith %.3f (ratio %.3f)\n",
                    o.h, o.ell, o.nodes, o.seconds, sr[i], o.pcNumeric * 1e-6, o.pcBracketLo * 1e-6, o.pcBracketHi * 1e-6, o.pcTheory * 1e-6, pr[i]);
    }
    const double sx = 2 * sr[1] - sr[0], px = 2 * pr[1] - pr[0];
    std::printf("    l -> 0 extrapolation: Sneddon %.3f, p_c %.3f\n", sx, px);
    check(ok, "phase-field runs complete and detect unstable growth");
    check(std::fabs(sr[1] - 1) < std::fabs(sr[0] - 1) && std::fabs(pr[1] - 1) < std::fabs(pr[0] - 1), "volume and critical pressure converge towards Sneddon / Griffith as l, h -> 0");
    check(std::fabs(sx - 1) < 0.10, "extrapolated crack-opening volume within 10 % of Sneddon (2 pi p a^2 / E')");
    check(std::fabs(px - 1) < 0.08, "extrapolated critical pressure within 8 % of Griffith sqrt(E' Gc / (pi a))");
    check(pr[1] > 0.8 && pr[1] < 1.0, "at l/a = 0.1 the diffuse crack fails 0-20 % below Griffith (known AT2 size effect)");
}
