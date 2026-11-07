// =====================================================================================================================
// C22 tests: research Cycle 4 / H1f integration (closed-form design potential, shut-in-burst control)
// =====================================================================================================================

static void testScheduleGain() {
    std::printf("== C22 H1f closed-form schedule gain (3-arc optimum incl. terminal burst) vs research/exp/h1f numerics\n");
    struct Ref { double Pi, th, u, g3, g2; };
    // numbers from research/exp/h1f.log (SLSQP optimum == 3-arc closed form to 1e-5; 2-arc = no terminal burst)
    const Ref refs[] = {{1, 0.05, 10, 0.02176, 0.01237}, {3, 0.05, 10, 0.07946, 0.05392}, {6, 0.05, 10, 0.18135, 0.14058},
                        {3, 0.02, 10, 0.03554, -1}, {1, 0.01, 3, 0.00418, -1}, {6, 0.02, 3, 0.06999, -1}};
    for (const Ref& r : refs) {
        const seis::GainEstimate g = seis::scheduleGain(r.Pi, r.th, r.u);
        std::printf("    Pi %g theta %.3g u %g: gain %.4f %% (ref %.4f %%), no-burst %.4f %% (ref %s), burst vol %.4f, s* %.4f, s_end %.4f\n", r.Pi, r.th, r.u,
                    100 * g.gain, 100 * r.g3, 100 * g.gainNoBurst, r.g2 > 0 ? std::to_string(100 * r.g2).c_str() : "-", g.burstFraction, g.sPlateau, g.sEnd);
        check(std::fabs(g.gain - r.g3) < 0.02 * r.g3, "closed-form 3-arc gain matches the numerical optimum within 2 %");
        if (r.g2 > 0) check(std::fabs(g.gainNoBurst - r.g2) < 0.02 * r.g2, "no-burst (2-arc) gain matches the research script within 2 %");
        check(g.gain >= g.gainNoBurst, "allowing the shut-in burst never lowers the attainable gain");
    }
    const seis::GainEstimate a = seis::scheduleGain(3, 0.005, 10), b = seis::scheduleGain(3, 0.01, 10);
    check(a.gain > 0 && std::fabs(b.gain / a.gain - 2) < 0.1, "small-theta regime: gain is linear in theta (theta law)");
    const seis::GainEstimate big = seis::scheduleGain(3, 3.0, 10);
    check(big.gain >= 0 && !big.exact && a.exact, "outside the validated range the reference is flagged (lower bound) and never negative");
}

static void testSeismicityLastInterval() {
    std::printf("== C22 seismicity designer: design potential + last-interval rate bound (no shut-in burst)\n");
    const char* deckText = R"(RUNSPEC
DIMENS
 41 41 1 /
OIL
WATER
METRIC
TABDIMS
 1 1 20 20 /
WELLDIMS
 1 1 1 1 /
START
 1 JAN 2025 /
GRID
DX
 1681*20 /
DY
 1681*20 /
DZ
 1681*10 /
TOPS
 1681*2000 /
PORO
 1681*0.1 /
PERMX
 1681*10 /
PERMY
 1681*10 /
PERMZ
 1681*1 /
PROPS
PVCDO
 200 1.0 1e-4 1.0 0 /
PVTW
 200 1.0 4e-5 0.5 0 /
DENSITY
 800 1000 1 /
ROCK
 200 4e-5 /
SWOF
 0.2 0 1 0
 0.3 0.02 0.7 0
 1.0 1 0 0 /
SOLUTION
EQUIL
 2005 200 1000 0 /
SCHEDULE
WELSPECS
 I G 21 21 2005 WATER /
/
COMPDAT
 I 21 21 1 1 OPEN 1* 1* 0.2 /
/
WCONINJE
 I WATER OPEN RATE 100 1* 900 /
/
TSTEP
 10*1 /
END
)";
    ressim::Case cs = ressim::eclipseDeckFromString(deckText, "");
    seis::DesignSetup s;
    s.wells = {"I"};
    seis::Patch p; p.name = "F1"; p.i = 26; p.j = 20; p.k = 0; p.aSigma = 0.1; p.taYears = 50; p.rate = 2;
    s.patches = {p};
    s.nCtrl = 6; s.injDays = 6; s.postDays = 12; s.subSteps = 2; s.postSteps = 8; s.totalVolume = 600; s.correctionIters = 0;
    const seis::DesignOutput d0 = seis::designInjection(cs, s);
    s.lastRateFactor = 1.0;
    const seis::DesignOutput d1 = seis::designInjection(cs, s);
    for (const auto& l : d0.logLines) if (l.find("Design potential") != std::string::npos || l.find("Optimum") != std::string::npos) std::printf("      %s\n", l.c_str());
    check(d0.ok && d1.ok, "both designs run");
    if (!(d0.ok && d1.ok)) return;
    std::printf("    free q/mean:"); for (double q : d0.qOpt[0]) std::printf(" %.3f", q / 100.0); std::printf("  (max rate factor %.3g)\n", s.maxRateFactor);
    std::printf("    free: Pi %.3g theta %.3g endsWithBurst %d, q_last/mean %.3f | bounded: q_last/mean %.3f, events %.5g vs %.5g\n", d0.Pi, d0.theta,
                int(d0.endsWithBurst), d0.qOpt[0].back() / (600.0 / 6), d1.qOpt[0].back() / (600.0 / 6), d1.eventsOpt, d0.eventsOpt);
    check(d0.Pi > 0 && d0.theta > 0 && (d0.farField || d0.gainTheory > 0), "design potential (Pi, theta, regime, first-order reference gain) is reported");
    check(d0.qOpt[0].front() > d0.qOpt[0][1], "optimum opens with the charging burst (first arc of the H1f bang-singular-bang structure)");
    check(d1.qOpt[0].back() <= (600.0 / 6) * (1 + 1e-6), "last-interval bound is honoured");
    check(!d1.endsWithBurst, "bounded design does not end with a shut-in burst");
    double vol = 0; for (double v : d1.qOpt[0]) vol += v * s.injDays / s.nCtrl;
    check(std::fabs(vol - 600) < 1e-6 * 600, "bounded design still injects the requested volume");
    check(d1.eventsOpt >= d0.eventsOpt * (1 - 1e-9), "extra constraint cannot beat the free optimum (model prediction)");
    check(d1.eventsOpt <= d1.eventsConst * (1 + 1e-9), "bounded design predicted no worse than constant rate");
}

static void testSeismicityNearField() {
    std::printf("== C22 seismicity designer, near-field patch: closed-form potential vs designed gain, shut-in burst control\n");
    const char* deckText = R"(RUNSPEC
DIMENS
 41 41 1 /
OIL
WATER
METRIC
TABDIMS
 1 1 20 20 /
WELLDIMS
 1 1 1 1 /
START
 1 JAN 2025 /
GRID
DX
 1681*20 /
DY
 1681*20 /
DZ
 1681*10 /
TOPS
 1681*2000 /
PORO
 1681*0.1 /
PERMX
 1681*10 /
PERMY
 1681*10 /
PERMZ
 1681*1 /
PROPS
PVCDO
 200 1.0 1e-4 1.0 0 /
PVTW
 200 1.0 4e-5 0.5 0 /
DENSITY
 800 1000 1 /
ROCK
 200 4e-5 /
SWOF
 0.2 0 1 0
 0.3 0.02 0.7 0
 1.0 1 0 0 /
SOLUTION
EQUIL
 2005 200 1000 0 /
SCHEDULE
WELSPECS
 I G 21 21 2005 WATER /
/
COMPDAT
 I 21 21 1 1 OPEN 1* 1* 0.2 /
/
WCONINJE
 I WATER OPEN RATE 100 1* 900 /
/
TSTEP
 10*1 /
END
)";
    ressim::Case cs = ressim::eclipseDeckFromString(deckText, "");
    seis::DesignSetup s;
    s.wells = {"I"};
    seis::Patch p; p.name = "NF"; p.i = 22; p.j = 21; p.k = 0; p.aSigma = 0.1; p.taYears = 50; p.rate = 2;
    s.patches = {p};
    s.nCtrl = 12; s.injDays = 6; s.postDays = 3; s.subSteps = 2; s.postSteps = 6; s.totalVolume = 600; s.correctionIters = 0;
    const seis::DesignOutput d0 = seis::designInjection(cs, s);
    s.lastRateFactor = 0.5;
    const seis::DesignOutput d1 = seis::designInjection(cs, s);
    for (const auto& l : d0.logLines) if (l.find("Design potential") != std::string::npos || l.find("Optimum") != std::string::npos) std::printf("      %s\n", l.c_str());
    check(d0.ok && d1.ok, "both designs run");
    if (!(d0.ok && d1.ok)) return;
    std::printf("    free q/mean:"); for (double q : d0.qOpt[0]) std::printf(" %.3f", q / 100.0);
    std::printf("\n    bounded q/mean:"); for (double q : d1.qOpt[0]) std::printf(" %.3f", q / 100.0);
    const double gD = 1 - d0.eventsOpt / d0.eventsConst, gB = 1 - d1.eventsOpt / d1.eventsConst;
    std::printf("\n    Pi %.3g theta %.3g farField %d burst %d | designed gain %.2f %% (closed form %.2f %%), no-burst %.2f %% (closed form %.2f %%)\n", d0.Pi,
                d0.theta, int(d0.farField), int(d0.endsWithBurst), 100 * gD, 100 * d0.gainTheory, 100 * gB, 100 * d0.gainTheoryNoBurst);
    check(!d0.farField && d0.gainTheory > 0, "near-field patch: first-order reference potential reported");
    check(gD <= d0.gainTheory + 1e-9, "designed gain (diffusive engine kernel) stays below the first-order reference (H1g: kernel shape)");
    const size_t n = d0.qOpt[0].size();
    check(d0.qOpt[0].front() > 1.5 * 100 && d0.qOpt[0][n - 1] > d0.qOpt[0][n - 2], "bang-singular-bang shape: charging burst first, terminal uptick last");
    check(!d1.endsWithBurst && d1.qOpt[0].back() <= 50 * (1 + 1e-6), "last-interval bound (0.5 x mean) is honoured");
    check(d1.eventsOpt >= d0.eventsOpt * (1 - 1e-9) && d1.eventsOpt <= d1.eventsConst * (1 + 1e-9), "bounded optimum lies between free optimum and constant rate");
}
