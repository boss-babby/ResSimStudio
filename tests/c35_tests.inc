// C35 tests: chemical EOR in the black-oil engine - surfactant (SURFACT), alkaline (ALKALINE) and foam (FOAM), with the
// user's ECLIPSE example decks decks/chem/ALKALINE_ASP.DATA (alkaline-surfactant-polymer) and decks/chem/FOAM.DATA.

namespace c35 {
static std::string readText(const std::string& path) {
    std::ifstream f(path); std::stringstream ss; ss << f.rdbuf(); return ss.str();
}
static std::string sub(std::string t, const std::string& a, const std::string& b, int* n = nullptr) {
    int k = 0;
    for (size_t p = t.find(a); p != std::string::npos; p = t.find(a, p + b.size())) { t.replace(p, a.size(), b); ++k; }
    if (n) *n = k;
    return t;
}
struct Run { bool ok = false; std::vector<ressim::SummaryPoint> sum; std::vector<double> inPlace, decayed; std::vector<std::string> names;
             double maxMisc = 0, minSt = 1, minFoam = 1; std::string err; };
static Run run(const std::string& deck, int stopAfter = -1, std::string* state = nullptr) {
    Run R;
    try {
        const ressim::Case c = ressim::eclipseDeckFromString(deck);
        ressim::Simulator sim;
        sim.setup(c); sim.initialize();
        int k = 0;
        while (sim.runReportStep()) { if (stopAfter >= 0 && ++k >= stopAfter) { if (state) *state = sim.saveState(); break; } }
        R.sum = sim.summary();
        for (size_t i = 0; i < sim.numChem(); ++i) { R.inPlace.push_back(sim.chemInPlace(i)); R.decayed.push_back(sim.chemDecayed(i)); R.names.push_back(sim.chemName(i)); }
        for (const auto& s : sim.snapshots())
            for (const auto& [n, v] : s.extra) {
                if (n == "SURFMISC") for (float x : v) R.maxMisc = std::max(R.maxMisc, double(x));
                if (n == "SURFSTM") for (float x : v) R.minSt = std::min(R.minSt, double(x));
                if (n == "FOAMMRF") for (float x : v) R.minFoam = std::min(R.minFoam, double(x));
            }
        R.ok = true;
    } catch (const std::exception& e) { R.err = e.what(); std::printf("    run failed: %s\n", e.what()); }
    return R;
}
// injected - produced - in place - decayed, relative to the injected amount
static double balance(const Run& R, size_t i, double* inj = nullptr) {
    if (R.sum.empty() || i >= R.inPlace.size()) return 1e9;
    const auto& p = R.sum.back();
    const size_t k = p.fcompIT.size() - R.inPlace.size() + i;   // chemicals follow the tracers
    const double I = k < p.fcompIT.size() ? p.fcompIT[k] : 0.0, P = k < p.fcompPT.size() ? p.fcompPT[k] : 0.0;
    if (inj) *inj = I;
    return std::fabs(I - P - R.inPlace[i] - R.decayed[i]) / std::max(I, 1e-30);
}
static double well(const ressim::SummaryPoint& p, size_t w, double ressim::WellRates::*m) { return w < p.wells.size() ? p.wells[w].*m : 0.0; }
}  // namespace c35

static void testC35Chem() {
    std::printf("== C35 chemical EOR: surfactant / alkaline / foam (ASP and FOAM example decks)\n");
    const std::string fa = findDeck("decks/chem/ALKALINE_ASP.DATA"), ff = findDeck("decks/chem/FOAM.DATA");
    if (fa.empty() || ff.empty()) { std::printf("  [SKIP] decks/chem not found\n"); return; }
    const std::string asp = c35::readText(fa), foam = c35::readText(ff);
    // ---- import: no engine-unsupported features, SURFNUM read, keywords kept verbatim
    {
        bool ok = false;
        try {
            const ressim::Case c = ressim::eclipseDeckFromString(asp);
            int nChem = 0;
            for (const auto& k : c.extraKeywords) if (k.name == "SURFST" || k.name == "SURFCAPD" || k.name == "ALSURFST" || k.name == "WSURFACT" || k.name == "WALKALIN") ++nChem;
            std::printf("    ASP import: %zu engine-unsupported, SURFNUM %zu values, %d chemical keywords kept\n", c.engineUnsupported.size(), c.grid.surfnum.size(), nChem);
            for (auto& u : c.engineUnsupported) std::printf("      %s\n", u.c_str());
            ok = c.engineUnsupported.empty() && c.grid.surfnum.size() == 40 && c.grid.surfnum[0] == 2 && nChem >= 9;
        } catch (const std::exception& e) { std::printf("    %s\n", e.what()); }
        check(ok, "C35: ALKALINE_ASP deck imports for the own engine (SURFACT / ALKALINE / POLYMER), SURFNUM and chemical keywords kept");
    }
    // ---- ASP run: mass balance of both chemicals, low IFT -> miscibility, more oil than without surfactant
    const c35::Run A = c35::run(asp);
    check(A.ok && A.names.size() == 2, "C35: ASP deck runs to the end with alkaline + surfactant components");
    if (A.ok && A.names.size() == 2) {
        double ia = 0, is = 0;
        const double ba = c35::balance(A, 0, &ia), bs = c35::balance(A, 1, &is);
        std::printf("    %s injected %.6g (balance error %.2e), %s injected %.6g (balance error %.2e); max miscibility %.3f, min IFT ratio %.3g\n",
                    A.names[0].c_str(), ia, ba, A.names[1].c_str(), is, bs, A.maxMisc, A.minSt);
        // 200 stb/d x 50 lb/stb x 100 d alkaline, 200 x 30 x 100 surfactant
        check(std::fabs(ia - 1.0e6) < 1e-6 * 1.0e6 && std::fabs(is - 6.0e5) < 1e-6 * 6.0e5, "C35: injected chemical masses = rate x concentration x time");
        check(ba < 1e-6 && bs < 1e-6, "C35: alkaline and surfactant mass balances (injected = produced + in solution + adsorbed) close");
        check(A.maxMisc > 0.05 && A.minSt < 1e-2, "C35: surfactant lowers the IFT (ALSURFST x SURFST) and the capillary number reaches the SURFCAPD miscible range");
        int n = 0;
        const std::string noSurf = c35::sub(asp, "WSURFACT\n 'I'  30.0 /", "WSURFACT\n 'I'  0.0 /", &n);
        const c35::Run B = c35::run(noSurf);
        if (B.ok) {
            const double fa1 = A.sum.back().fopt, fb1 = B.sum.back().fopt;
            std::printf("    FOPT with surfactant %.6g, without %.6g (%d WSURFACT records zeroed)\n", fa1, fb1, n);
            check(n == 2 && fa1 > fb1 * 1.01, "C35: surfactant (low IFT -> miscible curves) recovers more oil than the same flood without it");
        } else check(false, "C35: ASP deck without surfactant runs");
    }
    // ---- consistency: SURFCAPD = 0 everywhere with no surfactant injected equals the deck without the SURFACT option
    {
        std::string d0 = c35::sub(asp, "WSURFACT\n 'I'  30.0 /", "WSURFACT\n 'I'  0.0 /");
        d0 = c35::sub(d0, "  -2    1.0\n  10    1.0 /", "  -2    0.0\n  10    0.0 /");
        d0 = c35::sub(d0, "WALKALIN\n 'I' 50.0 /", "WALKALIN\n 'I' 0.0 /");   // ALSURFST would otherwise scale Pcow even without surfactant
        std::string d1 = d0;
        // remove the surfactant option and its data from d1
        d1 = c35::sub(d1, "\nSURFACT\n", "\n");
        for (const char* kw : {"SURFVISC", "SURFADS", "SURFST", "SURFCAPD", "SURFROCK", "SURFNUM"}) {
            const size_t p = d1.find(std::string("\n") + kw + "\n");
            if (p == std::string::npos) continue;
            size_t e = p + 1;
            const int tabs = std::string(kw) == "SURFNUM" ? 1 : 2;
            for (int t = 0; t < tabs; ++t) { e = d1.find('/', e); if (e == std::string::npos) break; ++e; }
            if (e != std::string::npos) d1.erase(p + 1, e - p - 1);
        }
        // drop every WSURFACT keyword: its records up to the terminating line that holds only "/"
        for (size_t p = d1.find("\nWSURFACT"); p != std::string::npos; p = d1.find("\nWSURFACT")) {
            size_t e = d1.find('\n', p + 1);
            while (e != std::string::npos) {
                const size_t nx = d1.find('\n', e + 1);
                std::string ln = d1.substr(e + 1, (nx == std::string::npos ? d1.size() : nx) - e - 1);
                ln.erase(0, ln.find_first_not_of(" \t"));
                e = nx;
                if (ln.rfind("/", 0) == 0) break;
            }
            d1.erase(p, (e == std::string::npos ? d1.size() : e) - p);
        }
        if (d1.find("TSTEP\n  50 /") == std::string::npos) std::printf("    warning: deck edit removed the last TSTEP\n");
        const c35::Run X = c35::run(d0), Y = c35::run(d1);
        bool ok = X.ok && Y.ok && X.names.size() == 2 && Y.names.size() == 1;
        if (ok) {
            const double a = X.sum.back().fopt, b = Y.sum.back().fopt;
            std::printf("    zero-surfactant SURFACT run FOPT %.10g, run without SURFACT %.10g\n", a, b);
            ok = std::fabs(a - b) <= 1e-7 * std::max(1.0, b);
        }
        check(ok, "C35: SURFACT with zero surfactant / alkaline and SURFCAPD = 0 reproduces the run without the option");
    }
    // ---- round trip: export -> import keeps the chemical model and the results
    {
        bool ok = false;
        try {
            const ressim::Case c = ressim::eclipseDeckFromString(asp);
            const std::string back = ressim::eclipseDeckToString(c);
            const c35::Run B = c35::run(back);
            if (B.ok && A.ok) {
                std::printf("    exported deck: FOPT %.10g vs %.10g; SURFNUM %s\n", B.sum.back().fopt, A.sum.back().fopt, back.find("SURFNUM") != std::string::npos ? "written" : "missing");
                ok = std::fabs(B.sum.back().fopt - A.sum.back().fopt) <= 1e-9 * A.sum.back().fopt && back.find("SURFNUM") != std::string::npos &&
                     back.find("WSURFACT") != std::string::npos && back.find("ALSURFAD") != std::string::npos;
            }
        } catch (const std::exception& e) { std::printf("    %s\n", e.what()); }
        check(ok, "C35: ASP deck round trip (export -> import) keeps SURFNUM, the chemical keywords and identical results");
    }
    // ---- restart (state version 8): stop after 2 report steps, continue from the saved state
    {
        std::string st;
        const c35::Run H = c35::run(asp, 2, &st);
        bool ok = false;
        if (H.ok && !st.empty() && A.ok) {
            try {
                const ressim::Case c = ressim::eclipseDeckFromString(asp);
                ressim::Simulator sim; sim.setup(c); sim.initialize();
                sim.runReportStep(); sim.runReportStep();
                sim.loadState(st);
                while (sim.runReportStep()) {}
                const double a = sim.summary().back().fopt, b = A.sum.back().fopt;
                std::printf("    restart after step 2: FOPT %.10g vs continuous %.10g\n", a, b);
                ok = std::fabs(a - b) <= 1e-8 * b;
            } catch (const std::exception& e) { std::printf("    %s\n", e.what()); }
        }
        check(ok, "C35: restart from a saved state with chemical components continues the run exactly");
    }
    // ---- FOAM deck: gas-phase foam with adsorption, decay and pressure / shear weakening
    {
        bool imp = false;
        try { const ressim::Case c = ressim::eclipseDeckFromString(foam); imp = c.engineUnsupported.empty(); for (auto& u : c.engineUnsupported) std::printf("      %s\n", u.c_str()); }
        catch (const std::exception& e) { std::printf("    %s\n", e.what()); }
        check(imp, "C35: FOAM example deck imports for the own engine");
        const c35::Run F = c35::run(foam);
        check(F.ok && F.names.size() == 1, "C35: FOAM deck runs to the end");
        if (F.ok && F.names.size() == 1) {
            double inj = 0;
            const double b = c35::balance(F, 0, &inj);
            std::printf("    foam injected %.6g, decayed %.6g, in place %.6g, balance error %.2e, min gas mobility factor %.3g\n", inj, F.decayed[0], F.inPlace[0], b, F.minFoam);
            check(inj > 0 && b < 1e-6 && F.decayed[0] > 0, "C35: foam mass balance with adsorption and first-order decay (FOAMDCYW / FOAMDCYO) closes");
            check(F.minFoam < 0.5, "C35: foam reduces the gas mobility (FOAMMOB with FOAMMOBP / FOAMMOBS weakening)");
            const c35::Run G = c35::run(c35::sub(foam, "'INJECTOR'  1.1 /", "'INJECTOR'  0.0 /"));
            if (G.ok) {
                // injector is rate controlled: a lower gas mobility needs a higher injection pressure
                double bf = 0, bg = 0;
                for (const auto& p : F.sum) bf = std::max(bf, c35::well(p, 1, &ressim::WellRates::bhp));
                for (const auto& p : G.sum) bg = std::max(bg, c35::well(p, 1, &ressim::WellRates::bhp));
                const double gf = F.sum.back().fgpt, gg = G.sum.back().fgpt;
                std::printf("    max injector BHP with foam %.6g, without %.6g; FGPT %.6g vs %.6g\n", bf, bg, gf, gg);
                check(bf > bg && gf < gg, "C35: foam raises the gas injection pressure and delays / reduces the produced gas");
            } else check(false, "C35: FOAM deck without foam runs");
        }
    }
    // ---- OPM Flow reference: opm-tests spe1_foam (tabular foam model, adsorption) - relative RMS over the report times
    {
        const std::string fd = findDeck("decks/chem/SPE1FOAM.DATA"), fr = findDeck("decks/chem/opm-ref/SPE1FOAM.SMSPEC");
        if (fd.empty() || fr.empty()) std::printf("  [SKIP] decks/chem/SPE1FOAM not found\n");
        else {
            ressim::ecl::Summary ref;
            std::string err;
            const c35::Run R = c35::run(c35::readText(fd));
            bool ok = R.ok && ressim::ecl::readSummary(fr, ref, &err);
            if (ok) {
                auto rel = [&](const char* key, std::function<double(const ressim::SummaryPoint&)> f) {
                    const int v = ref.find(key);
                    if (v < 0) return 1e9;
                    double se = 0, sr = 0;
                    for (const auto& p : R.sum) {
                        if (p.time <= 0) continue;
                        const double r = ref.at(v, p.time);
                        if (!std::isfinite(r)) continue;
                        se += (f(p) - r) * (f(p) - r); sr += r * r;
                    }
                    return sr > 0 ? std::sqrt(se / sr) : 1e9;
                };
                const double eo = rel("FOPR", [](const ressim::SummaryPoint& p) { return p.fopr; });
                const double eg = rel("FGOR", [](const ressim::SummaryPoint& p) { return p.fopr > 0 ? p.fgpr / p.fopr : 0.0; });
                const double eb = rel("WBHP:INJ", [](const ressim::SummaryPoint& p) { return c35::well(p, 1, &ressim::WellRates::bhp); });
                std::printf("    SPE1FOAM vs OPM Flow: relRMS FOPR %.2f %%, FGOR %.2f %%, WBHP:INJ %.2f %%\n", 100 * eo, 100 * eg, 100 * eb);
                ok = eo < 0.02 && eg < 0.02 && eb < 0.02;
            } else std::printf("    %s\n", err.c_str());
            check(ok, "C35: opm-tests SPE1FOAM matches the OPM Flow reference (FOPR, FGOR, injector BHP within 2 % relRMS)");
        }
    }
}

// ---- C35-3 SOLVENT: four-component solvent model (opm-tests spe1_solvent decks, OPM Flow reference)
namespace c35s {
struct Run { bool ok = false; std::vector<ressim::SummaryPoint> sum; double nip = 0, ssMax = 0; std::string state; };
static Run run(const std::string& deck, int stopAfter = -1) {
    Run R;
    try {
        const ressim::Case c = ressim::eclipseDeckFromString(deck);
        ressim::Simulator sim; sim.setup(c); sim.initialize();
        int k = 0;
        while (sim.runReportStep()) if (stopAfter >= 0 && ++k >= stopAfter) { R.state = sim.saveState(); break; }
        R.sum = sim.summary(); R.nip = sim.solventInPlace();
        for (const auto& s : sim.snapshots()) for (float x : s.field("SSOLVENT")) R.ssMax = std::max(R.ssMax, double(x));
        R.ok = !R.sum.empty();
    } catch (const std::exception& e) { std::printf("    run failed: %s\n", e.what()); }
    return R;
}
}  // namespace c35s

static void testC35Solvent() {
    std::printf("== C35 SOLVENT: four-component solvent model (spe1_solvent decks)\n");
    const std::string ff = findDeck("decks/solvent/SPE1CASE2_SOLVENT_FOAM.DATA"), fs = findDeck("decks/solvent/SPE1CASE2_SOLVENT.DATA"),
                      fr = findDeck("decks/solvent/opm-ref/SPE1CASE2_SOLVENT_FOAM.SMSPEC");
    if (ff.empty() || fs.empty() || fr.empty()) { std::printf("  [SKIP] decks/solvent not found\n"); return; }
    const std::string dF = c35::readText(ff), dS = c35::readText(fs);
    {
        bool ok = false;
        try {
            const ressim::Case c = ressim::eclipseDeckFromString(dS);
            int n = 0;
            for (const auto& k : c.extraKeywords) if (k.name == "SOLVENT" || k.name == "PVDS" || k.name == "SSFN" || k.name == "SDENSITY" || k.name == "WSOLVENT") ++n;
            for (auto& u : c.engineUnsupported) std::printf("      %s\n", u.c_str());
            std::printf("    import: %zu engine-unsupported, %d solvent keywords kept\n", c.engineUnsupported.size(), n);
            ok = c.engineUnsupported.empty() && n >= 5;
        } catch (const std::exception& e) { std::printf("    %s\n", e.what()); }
        check(ok, "C35: SOLVENT decks import for the own engine (SOLVENT, PVDS, SDENSITY, SSFN, WSOLVENT kept verbatim)");
    }
    // ---- immiscible solvent flood: solvent balance, solvent saturation
    const c35s::Run S = c35s::run(dS);
    check(S.ok, "C35: SPE1CASE2_SOLVENT runs to the end");
    if (S.ok) {
        const auto& p = S.sum.back();
        const double err = std::fabs(p.fnit - p.fnpt - S.nip) / std::max(p.fnit, 1e-30);
        std::printf("    solvent injected %.6g, produced %.6g, in place %.6g (balance error %.2e), max SSOLVENT %.3f, FOPT %.6g\n", p.fnit, p.fnpt, S.nip, err, S.ssMax, p.fopt);
        check(p.fnit > 0 && err < 1e-3 && S.ssMax > 0.05, "C35: solvent mass balance (injected = produced + in place) closes; the solvent occupies pore space (SSOLVENT)");
        // WSOLVENT 0: pure gas injection must equal the run without the SOLVENT option
        int n1 = 0, n2 = 0;
        const std::string g0 = c35::sub(dS, "'INJ'  1.0/", "'INJ'  0.0/", &n1);
        const std::string g1 = c35::sub(g0, "\nSOLVENT\n", "\n", &n2);
        const c35s::Run A = c35s::run(g0), B = c35s::run(g1);
        bool ok = A.ok && B.ok && n1 == 1 && n2 == 1;
        if (ok) {
            const double a = A.sum.back().fopt, b = B.sum.back().fopt, ga = A.sum.back().fgpt, gb = B.sum.back().fgpt;
            std::printf("    WSOLVENT 0: FOPT %.10g / FGPT %.10g; without SOLVENT: FOPT %.10g / FGPT %.10g\n", a, ga, b, gb);
            ok = std::fabs(a - b) <= 1e-6 * b && std::fabs(ga - gb) <= 1e-6 * gb;
        }
        check(ok, "C35: SOLVENT with no solvent injected reproduces the run without the option (z = 0 limit)");
        check(A.ok && p.fgpt < A.sum.back().fgpt, "C35: injecting solvent instead of gas lowers the produced (hydrocarbon) gas");
    }
    // ---- OPM Flow reference: SPE1CASE2_SOLVENT_FOAM (immiscible solvent carrying foam, FOAMOPTS SOLVENT)
    {
        ressim::ecl::Summary ref; std::string err;
        const c35s::Run F = c35s::run(dF);
        bool ok = F.ok && ressim::ecl::readSummary(fr, ref, &err);
        if (ok) {
            auto rel = [&](const char* key, std::function<double(const ressim::SummaryPoint&)> f) {
                const int v = ref.find(key);
                if (v < 0) return 1e9;
                double se = 0, sr = 0;
                for (const auto& q : F.sum) {
                    if (q.time <= 0) continue;
                    const double r = ref.at(v, q.time);
                    if (!std::isfinite(r)) continue;
                    se += (f(q) - r) * (f(q) - r); sr += r * r;
                }
                return sr > 0 ? std::sqrt(se / sr) : 1e9;
            };
            const double eg = rel("FGOR", [](const ressim::SummaryPoint& q) { return q.fopr > 0 ? q.fgpr / q.fopr : 0.0; });
            const double ei = rel("WBHP:INJ", [](const ressim::SummaryPoint& q) { return c35::well(q, 1, &ressim::WellRates::bhp); });
            const double ep = rel("WBHP:PROD", [](const ressim::SummaryPoint& q) { return c35::well(q, 0, &ressim::WellRates::bhp); });
            std::printf("    SPE1CASE2_SOLVENT_FOAM vs OPM Flow: relRMS FGOR %.2f %%, WBHP:INJ %.2f %%, WBHP:PROD %.2f %%\n", 100 * eg, 100 * ei, 100 * ep);
            ok = eg < 0.02 && ei < 0.02 && ep < 0.02;
        } else std::printf("    %s\n", err.c_str());
        check(ok, "C35: opm-tests SPE1CASE2_SOLVENT_FOAM matches the OPM Flow reference (FGOR, both BHPs within 2 % relRMS)");
    }
    // ---- round trip and restart (state version 9)
    if (S.ok) {
        bool ok = false;
        try {
            const std::string back = ressim::eclipseDeckToString(ressim::eclipseDeckFromString(dS));
            const c35s::Run B = c35s::run(back);
            ok = B.ok && std::fabs(B.sum.back().fpr - S.sum.back().fpr) <= 1e-9 * S.sum.back().fpr && std::fabs(B.nip - S.nip) <= 1e-9 * S.nip &&
                 std::fabs(B.sum.back().fopt - S.sum.back().fopt) <= 1e-9 * S.sum.back().fopt && back.find("WSOLVENT") != std::string::npos;
            std::printf("    exported deck: FPR %.10g vs %.10g, FNIP %.10g vs %.10g\n", B.ok ? B.sum.back().fpr : 0.0, S.sum.back().fpr, B.ok ? B.nip : 0.0, S.nip);
        } catch (const std::exception& e) { std::printf("    %s\n", e.what()); }
        check(ok, "C35: SOLVENT deck round trip (export -> import) gives identical results");
        const c35s::Run H = c35s::run(dS, 30);
        ok = false;
        if (!H.state.empty()) {
            try {
                ressim::Simulator sim; sim.setup(ressim::eclipseDeckFromString(dS)); sim.initialize();
                sim.loadState(H.state);
                while (sim.runReportStep()) {}
                const double a = sim.summary().back().fpr, b = S.sum.back().fpr, na = sim.solventInPlace(), nb = S.nip;
                const double ia = sim.summary().back().fnit, ib = S.sum.back().fnit;
                std::printf("    restart after 30 report steps: FPR %.10g vs %.10g, FNIP %.10g vs %.10g, FNIT %.8g vs %.8g\n", a, b, na, nb, ia, ib);
                ok = std::fabs(a - b) <= 1e-7 * b && std::fabs(na - nb) <= 1e-6 * nb && std::fabs(ia - ib) <= 1e-9 * ib;
            } catch (const std::exception& e) { std::printf("    %s\n", e.what()); }
        }
        check(ok, "C35: restart (state version 9) with solvent continues the run exactly");
    }
}

// ---------------------------------------------------------------- C35-4 gas-lift optimisation (LIFTOPT / WLIFTOPT / GLIFTOPT)
namespace c35g {
struct Run { bool ok = false; std::vector<ressim::SummaryPoint> sum; std::vector<std::string> wells; std::string state; };
static Run run(const std::string& deck, const std::string& dir = "", int stopAfter = -1) {
    Run R;
    try {
        const ressim::Case c = ressim::eclipseDeckFromString(deck, dir);
        for (const auto& w : c.wells) R.wells.push_back(w.name);
        ressim::Simulator sim; sim.setup(c); sim.initialize();
        int k = 0;
        while (sim.runReportStep()) if (stopAfter >= 0 && ++k >= stopAfter) { R.state = sim.saveState(); break; }
        R.sum = sim.summary();
        R.ok = !R.sum.empty();
    } catch (const std::exception& e) { std::printf("    run failed: %s\n", e.what()); }
    return R;
}
}  // namespace c35g

static void testC35GasLift() {
    std::printf("== C35 gas-lift optimisation (opm-tests gaslift: GASLIFT-01, GASLIFT-04 vs the commercial simulator)\n");
    const std::string f1 = findDeck("decks/gaslift/GASLIFT-01.DATA"), f4 = findDeck("decks/gaslift/GASLIFT-04.DATA");
    if (f1.empty() || f4.empty()) { std::printf("  [SKIP] decks/gaslift not found\n"); return; }
    const std::string d1 = c35::readText(f1), d4 = c35::readText(f4);
    const std::string dir4 = f4.substr(0, f4.find_last_of("/\\") + 1);
    {
        bool ok = false;
        try {
            const ressim::Case c = ressim::eclipseDeckFromString(d1);
            int n = 0, bad = 0;
            for (const auto& k : c.extraKeywords) if (k.name == "LIFTOPT" || k.name == "WLIFTOPT" || k.name == "GLIFTOPT") ++n;
            for (const auto& u : c.engineUnsupported) if (u.find("LIFTOPT") != std::string::npos) { ++bad; std::printf("      %s\n", u.c_str()); }
            std::printf("    import: %d gas-lift keywords kept, %d reported unsupported\n", n, bad);
            ok = n >= 7 && bad == 0;
        } catch (const std::exception& e) { std::printf("    %s\n", e.what()); }
        check(ok, "C35: LIFTOPT / WLIFTOPT / GLIFTOPT import as engine keywords (kept verbatim)");
    }
    const c35g::Run A = c35g::run(d1);
    int nOff = 0;
    const std::string d1off = c35::sub(d1, "LIFTOPT\n         4.0E4", "LIFTOPT\n         0.0  ", &nOff);
    const c35g::Run B = c35g::run(d1off);
    check(A.ok && B.ok && nOff == 2, "C35: GASLIFT-01 runs with and without gas-lift optimisation");
    {   // C36: a GCONPROD limit with procedure NONE (not the control mode) does not cap the lift allocation (OPM has_control)
        int nn = 0;
        const c35g::Run N = c35g::run(c35::sub(d1, "3E3    RATE", "3E3    NONE", &nn));
        double mA = 0, mN = 0; for (const auto& p : A.sum) mA = std::max(mA, p.fglir); for (const auto& p : N.sum) mN = std::max(mN, p.fglir);
        std::printf("    PLAT-2 LRAT 3000 with procedure RATE: max FGLIR %.4g, FGLIT %.4g; NONE: %.4g, %.4g\n", mA, A.ok ? A.sum.back().fglit : 0.0, mN, N.ok ? N.sum.back().fglit : 0.0);
        check(nn == 1 && N.ok && A.ok && mN > 1.05 * mA && N.sum.back().fglit > A.sum.back().fglit,
              "C36: GCONPROD limits with exceed procedure NONE do not cap the gas-lift allocation (OPM GCONPROD / GasLiftGroupInfo semantics)");
    }
    if (A.ok && B.ok) {
        const double a = A.sum.back().fopt, b = B.sum.back().fopt, ecl = 10.25e6;   // ECL: opm-tests gaslift/plots (digitised)
        std::printf("    GASLIFT-01 FOPT with gas lift %.4g, without %.4g sm3 (commercial simulator %.4g, OPM Flow reference plot ~0.9e6); FGLIT %.4g sm3\n",
                    a, b, ecl, A.sum.back().fglit);
        check(a > 5 * b && std::fabs(a - ecl) < 0.10 * ecl,
              "C35: lift gas makes the dead wells flow - GASLIFT-01 FOPT within 10 % of the commercial simulator (OPM Flow: wells not flowing)");
        // allocation rules: whole increments (or the well maximum = largest table ALQ), lift off while WLIFTOPT NO 0,
        // GLIFTOPT PLAT-1 9E5 from 2022-08-01, LIFTOPT increment 0 = no lift gas
        bool inc = true, off = true, grp = true, any = false;
        double worstGrp = 0;
        for (const auto& p : A.sum) {
            double plat1 = 0;
            for (size_t w = 0; w < p.wells.size() && w < A.wells.size(); ++w) {
                const double L = p.wells[w].lift;
                if (L > 0) any = true;
                const double r = std::fmod(L, 4.0e4);
                if (L > 0 && std::fabs(L - 4.0e5) > 1e-6 && std::min(r, 4.0e4 - r) > 1e-6) inc = false;
                if (A.wells[w].rfind("OP-A", 0) == 0 || A.wells[w].rfind("OP-B", 0) == 0) plat1 += p.wells[w].eff * L;
            }
            if (p.time > 213.5 && p.time <= 274.5 && p.fglir > 0) off = false;
            if (p.time > 578.5) { worstGrp = std::max(worstGrp, plat1); if (plat1 > 9.0e5 * (1 + 1e-9)) grp = false; }
        }
        std::printf("    largest PLAT-1 lift gas after the GLIFTOPT limit: %.6g (limit 9e5); FGLIT without optimisation %.4g\n", worstGrp, B.sum.back().fglit);
        check(any && inc, "C35: lift gas is allocated in whole LIFTOPT increments (or the well maximum)");
        check(off, "C35: WLIFTOPT NO 0.0 fixes the lift gas at 0 (2021-08-01 to 2021-10-01)");
        check(grp && worstGrp > 0, "C35: GLIFTOPT maximum lift gas supply of group PLAT-1 is respected");
        check(B.sum.back().fglit == 0, "C35: LIFTOPT increment 0 switches the optimisation off (no lift gas)");
        // economic gradient: a gradient above every incremental oil / gas ratio allocates nothing
        int nE = 0;
        const std::string dE = c35::sub(d1, "4.0E4  0.0    0.0", "4.0E4  10.0   0.0", &nE);
        const c35g::Run E = c35g::run(dE);
        bool okE = E.ok && nE == 1;
        if (okE) {
            double before = 0;   // lift gas before the second LIFTOPT (2022-01-01, day 366; START 2020-12-31)
            for (const auto& p : E.sum) if (p.time < 365.5) before = std::max(before, p.fglir);
            std::printf("    minimum economic gradient 10 sm3/sm3: largest FGLIR before 2022-01 %.4g\n", before);
            okE = before == 0;
        }
        check(okE, "C35: the minimum economic gradient stops the allocation (LIFTOPT item 2)");
        // restart (state version 10) and export -> import
        const c35g::Run H = c35g::run(d1, "", 12);
        bool okR = false;
        if (!H.state.empty()) {
            try {
                ressim::Simulator sim; sim.setup(ressim::eclipseDeckFromString(d1)); sim.initialize();
                sim.loadState(H.state);
                while (sim.runReportStep()) {}
                const auto& q = sim.summary().back();
                std::printf("    restart after 12 report steps: FOPT %.10g vs %.10g, FGLIT %.10g vs %.10g\n", q.fopt, A.sum.back().fopt, q.fglit, A.sum.back().fglit);
                okR = std::fabs(q.fopt - A.sum.back().fopt) <= 1e-7 * A.sum.back().fopt && std::fabs(q.fglit - A.sum.back().fglit) <= 1e-7 * A.sum.back().fglit;
            } catch (const std::exception& e) { std::printf("    %s\n", e.what()); }
        }
        check(okR, "C35: restart (state version 10) with gas lift continues the run exactly");
        bool okX = false;
        try {
            const std::string back = ressim::eclipseDeckToString(ressim::eclipseDeckFromString(d1));
            const c35g::Run X = c35g::run(back);
            okX = X.ok && std::fabs(X.sum.back().fopt - A.sum.back().fopt) <= 1e-9 * A.sum.back().fopt && back.find("WLIFTOPT") != std::string::npos &&
                  back.find("GLIFTOPT") != std::string::npos;
        } catch (const std::exception& e) { std::printf("    %s\n", e.what()); }
        check(okX, "C35: gas-lift deck round trip (export -> import) gives identical results");
    }
    // ---- GASLIFT-04: three-phase MODEL05 PVT, ALQ (GRAT) lift curves - commercial simulator results (opm-tests plots)
    {
        const c35g::Run G = c35g::run(d4, dir4);
        bool ok = G.ok;
        if (ok) {
            const auto& p = G.sum.back();
            const double eo = 12.55e6, eg = 4.42e8, ep = 116.5;   // ECL FOPT, FGPT, FPR at 2026-01-01 (digitised, +-1 %)
            std::printf("    GASLIFT-04 FOPT %.4g (ECL %.4g, %+.1f %%), FGPT %.4g (ECL %.4g, %+.1f %%), FPR %.1f (ECL %.1f), FGLIT %.4g\n", p.fopt, eo,
                        100 * (p.fopt / eo - 1), p.fgpt, eg, 100 * (p.fgpt / eg - 1), p.fpr, ep, p.fglit);
            double plateau = 1e30;   // field ORAT 10 000 while the lift gas is on (first 200 days)
            for (const auto& q : G.sum) if (q.time > 20 && q.time < 200) plateau = std::min(plateau, q.fopr);
            std::printf("    lowest FOPR on the 10 000 sm3/d plateau (day 20-200): %.1f\n", plateau);
            ok = std::fabs(p.fopt / eo - 1) < 0.03 && std::fabs(p.fgpt / eg - 1) < 0.03 && std::fabs(p.fpr / ep - 1) < 0.04 && plateau > 9500;
        }
        check(ok, "C35: GASLIFT-04 matches the commercial simulator (FOPT, FGPT within 3 %, FPR within 4 %, ORAT plateau held by lift gas)");
    }
}

// ============================================================================ C35-4b surface network / facilities
namespace c35f {
struct Run { bool ok = false; std::vector<ressim::SummaryPoint> sum; std::vector<std::string> wells; ressim::UnitSystem units = ressim::UnitSystem::Metric;
             std::vector<ressim::Simulator::NetNode> nodes; std::string state; };
static Run run(const std::string& deck, const std::string& dir = "", int stopAfter = -1) {
    Run R;
    try {
        const ressim::Case c = ressim::eclipseDeckFromString(deck, dir);
        for (const auto& w : c.wells) R.wells.push_back(w.name);
        R.units = c.units;
        ressim::Simulator sim; sim.setup(c); sim.initialize();
        int k = 0;
        while (sim.runReportStep()) if (stopAfter >= 0 && ++k >= stopAfter) { R.state = sim.saveState(); break; }
        R.sum = sim.summary(); R.nodes = sim.networkNodes();
        R.ok = !R.sum.empty();
    } catch (const std::exception& e) { std::printf("    run failed: %s\n", e.what()); }
    return R;
}
static double nodeP(const Run& R, const std::string& n) { for (const auto& N : R.nodes) if (N.name == n) return N.p; return -1; }
}  // namespace c35f

static void testC35Facilities() {
    std::printf("== C35 surface network & facilities (BRANPROP 0 / 9999, GRUPNET, GSATPROD, GCONSUMP, GCONSALE, automatic choke)\n");
    const std::string fn = findDeck("examples/NETWORK.DATA");
    if (fn.empty()) { std::printf("  [SKIP] examples/NETWORK.DATA not found\n"); return; }
    std::string base = c35::readText(fn);
    int n0 = 0;
    base = c35::sub(base, " 40*91.3125", " 4*91.3125", &n0);
    check(n0 == 1, "C35 facilities: NETWORK example shortened to 4 report steps");
    auto withNet = [&](const std::string& net, const std::string& extra) {
        int n = 0;
        std::string d = c35::sub(base, "NODEPROP\n 'FIELD' 300 /\n 'PLAT' /\n/\n\nBRANPROP\n 'PLAT' 'FIELD' 2 /\n/\n", net + extra, &n);
        return n == 1 ? d : std::string();
    };
    const std::string ext = "NODEPROP\n 'FIELD' 300 /\n 'PLAT' /\n/\n\nBRANPROP\n 'PLAT' 'FIELD' 2 /\n/\n";
    const std::string satCons = "GSATPROD\n 'SAT' 400 100 800 1* 0 /\n/\nGCONSUMP\n 'PLAT' 300 50 /\n/\n";
    const c35f::Run B = c35f::run(withNet(ext, ""));
    const c35f::Run S = c35f::run(withNet(ext, satCons));
    {   // import: the facility keywords are engine keywords
        bool ok = false;
        try {
            const ressim::Case c = ressim::eclipseDeckFromString(withNet(ext, satCons + "GCONSALE\n 'FIELD' 100 200 50 RATE /\n/\nNEFAC\n 'PLAT' 0.9 /\n/\n"));
            int n = 0, bad = 0;
            for (const auto& k : c.extraKeywords) if (k.name == "GSATPROD" || k.name == "GCONSUMP" || k.name == "GCONSALE" || k.name == "NEFAC") ++n;
            for (const auto& u : c.engineUnsupported) if (u.find("GSATPROD") != std::string::npos || u.find("GCONSUMP") != std::string::npos || u.find("GCONSALE") != std::string::npos) ++bad;
            std::printf("    import: %d facility keywords kept, %d reported unsupported\n", n, bad);
            ok = n == 4 && bad == 0;
        } catch (const std::exception& e) { std::printf("    %s\n", e.what()); }
        check(ok, "C35: GSATPROD / GCONSUMP / GCONSALE / NEFAC import as engine keywords (kept verbatim)");
    }
    check(B.ok && S.ok, "C35: network example runs with satellite production and gas consumption");
    if (B.ok && S.ok) {
        double worst = 0, worstC = 0;
        for (const auto& p : S.sum) {
            if (p.time <= 0) continue;   // initial point (before the schedule)
            double o = 0, w = 0; for (const auto& r : p.wells) { o += r.eff * r.oil; w += r.eff * r.water; }
            worst = std::max({worst, std::fabs(p.fopr - o - 400), std::fabs(p.fwpr - w - 100)});
            worstC = std::max({worstC, std::fabs(p.fgcr - 300), std::fabs(p.fgimr - 50), std::fabs(p.fgsr - (p.fgpr - p.fgir - p.fgcr + p.fgimr))});
        }
        const auto& L = S.sum.back();
        std::printf("    satellite bookkeeping |FOPR - wells - 400| max %.2e; FGCR %.6g FGIMR %.6g FGSR %.6g; FGCT %.8g vs %.8g\n", worst, L.fgcr, L.fgimr, L.fgsr, L.fgct, 300 * L.time);
        check(worst < 1e-6 && worstC < 1e-6 && std::fabs(L.fgct - 300 * (L.time - S.sum.front().time + S.sum.front().dt)) < 1e-6 * L.fgct,
              "C35: GSATPROD rates in FOPR / FWPR, GCONSUMP -> FGCR / FGIMR / FGSR = FGPR - FGIR - FGCR + FGIMR, totals");
        // satellite flow through the PLAT node raises its pressure (more flow through the branch)
        const double pB = c35f::nodeP(B, "PLAT"), pS = c35f::nodeP(S, "PLAT");
        std::printf("    PLAT node pressure: base %.4f, with satellite + consumption %.4f\n", pB, pS);
        check(pB > 0 && pS > 0 && std::fabs(pS - pB) > 1e-3, "C35: satellite production / consumption change the network node flow and pressure");
    }
    {   // GCONSUMP 'NO-NODES' leaves the network alone
        const c35f::Run N = c35f::run(withNet(ext, "GCONSUMP\n 'PLAT' 3000 /\n/\n")), M = c35f::run(withNet(ext, "GCONSUMP\n 'PLAT' 3000 1* 'NO-NODES' /\n/\n"));
        bool ok = N.ok && M.ok && B.ok;
        if (ok) {
            const double a = c35f::nodeP(B, "PLAT"), b = c35f::nodeP(M, "PLAT"), c = c35f::nodeP(N, "PLAT");
            std::printf("    PLAT pressure: no consumption %.5f, 'NO-NODES' %.5f, consumption at the node %.5f\n", a, b, c);
            ok = std::fabs(a - b) < 1e-9 * a && std::fabs(a - c) > 1e-4;
        }
        check(ok, "C35: GCONSUMP at the source node reduces the branch gas flow; 'NO-NODES' has no network effect (OPM)");
    }
    {   // BRANPROP 9999 = no-loss flow line, 0 = branch removed (OPM / ECLIPSE)
        const c35f::Run N = c35f::run(withNet("NODEPROP\n 'FIELD' 300 /\n 'PLAT' /\n/\n\nBRANPROP\n 'PLAT' 'FIELD' 9999 /\n/\n", ""));
        const c35f::Run D = c35f::run(withNet(ext, "") .size() ? c35::sub(withNet(ext, ""), "TSTEP\n 4*91.3125", "TSTEP\n 91.3125 /\nBRANPROP\n 'PLAT' 'FIELD' 0 /\n/\nTSTEP\n 3*91.3125") : std::string());
        bool ok = N.ok && D.ok;
        if (ok) {
            const double p = c35f::nodeP(N, "PLAT");
            double thpEnd = 0; for (const auto& r : D.sum.back().wells) thpEnd = std::max(thpEnd, r.thpNet);
            std::printf("    9999: PLAT %.6f (FIELD 300); after BRANPROP 0: largest network THP %.3g, PLAT parent '%s'\n", p, thpEnd, D.nodes.empty() ? "-" : D.nodes.back().parent.c_str());
            ok = std::fabs(p - 300) < 1e-9 && thpEnd == 0;
        }
        check(ok, "C35: BRANPROP VFP 9999 = no pressure loss, VFP 0 removes the branch (wells leave the network)");
    }
    {   // GRUPNET standard network == the same extended network
        const c35f::Run G = c35f::run(withNet("GRUPNET\n 'FIELD' 300 /\n 'PLAT' 1* 2 /\n/\n", ""));
        bool ok = G.ok && B.ok;
        if (ok) { std::printf("    GRUPNET FOPT %.10g vs NODEPROP / BRANPROP %.10g\n", G.sum.back().fopt, B.sum.back().fopt); ok = std::fabs(G.sum.back().fopt - B.sum.back().fopt) <= 1e-9 * B.sum.back().fopt; }
        check(ok, "C35: GRUPNET standard network gives the same run as the equivalent extended network");
    }
    {   // GCONPROD procedure NONE: a limit that is only reported (OPM / ECLIPSE); RATE enforces it
        const c35f::Run N = c35f::run(withNet(ext, "GCONPROD\n 'PLAT' NONE 1000 4* /\n/\n")), R = c35f::run(withNet(ext, "GCONPROD\n 'PLAT' NONE 1000 3* RATE /\n/\n"));
        bool ok = N.ok && R.ok && B.ok;
        if (ok) {
            double maxR = 0; for (const auto& p : R.sum) maxR = std::max(maxR, p.fopr);
            std::printf("    FOPT base %.8g, ORAT limit 1000 NONE %.8g, RATE %.8g (max FOPR %.6g)\n", B.sum.back().fopt, N.sum.back().fopt, R.sum.back().fopt, maxR);
            ok = std::fabs(N.sum.back().fopt - B.sum.back().fopt) <= 1e-9 * B.sum.back().fopt && maxR <= 1000 * 1.002 && B.sum.back().fopr > 1000;
        }
        check(ok, "C35: GCONPROD exceed procedure NONE only reports the limit, RATE enforces it");
        try {
            const ressim::Case c = ressim::eclipseDeckFromString(withNet(ext, "GCONPROD\n 'PLAT' NONE 1000 1* 1* 1* NONE 1* 1* 1* RATE 1* WELL /\n/\n"));
            bool got = false;
            for (const auto& st : c.schedule) for (const auto& g : st.gprod) if (g.group == "PLAT") got = g.actW == "RATE" && g.actG.empty() && g.actL == "WELL" && g.actionFor(ressim::Control::ORAT) == "NONE";
            check(got, "C35: GCONPROD items 11-13 (water / gas / liquid procedures, default item 7) are read");
        } catch (const std::exception& e) { std::printf("    %s\n", e.what()); check(false, "C35: GCONPROD items 11-13 are read"); }
    }
    {   // restart with satellite / consumption totals (state version 10)
        const c35f::Run H = c35f::run(withNet(ext, satCons), "", 2);
        bool ok = false;
        if (H.ok && !H.state.empty() && S.ok) {
            try {
                const ressim::Case c = ressim::eclipseDeckFromString(withNet(ext, satCons));
                ressim::Simulator sim; sim.setup(c); sim.initialize();
                sim.runReportStep(); sim.runReportStep();
                sim.loadState(H.state);
                while (sim.runReportStep()) {}
                const auto& a = sim.summary().back(); const auto& b = S.sum.back();
                std::printf("    restart after step 2: FOPT %.10g vs %.10g, FGCT %.10g vs %.10g, FGST %.10g vs %.10g\n", a.fopt, b.fopt, a.fgct, b.fgct, a.fgst, b.fgst);
                ok = std::fabs(a.fopt - b.fopt) <= 1e-8 * b.fopt && std::fabs(a.fgct - b.fgct) <= 1e-8 * b.fgct && std::fabs(a.fgst - b.fgst) <= 1e-8 * std::fabs(b.fgst);
            } catch (const std::exception& e) { std::printf("    %s\n", e.what()); }
        }
        check(ok, "C35: restart continues satellite / consumption / sales totals exactly");
    }
    // ---- opm-tests decks vs OPM Flow (python tools/fetch_opm_tests.py --dirs network,model2; or RESSIM_OPMT=<root>)
    const char* env = std::getenv("RESSIM_OPMT");
    const std::string root = env ? env : "opm-tests";
    auto exists = [](const std::string& p) { std::ifstream f(p); return f.good(); };
    struct Ref { const char* dir; const char* deck; const char* keys; double tol; };
    const Ref refs[] = {
        {"network", "NETWORK-01", "FOPR,FGPR,WBHP:PROD1,WBHP:PROD2", 1.0},
        {"network", "NETWORK-01_STANDARD", "FOPR,FGPR,WBHP:PROD1,WBHP:PROD2", 1.0},
        {"network", "NETWORK-01-REROUTE", "FOPR,FGPR,WBHP:PROD1,WBHP:PROD2", 1.0},
        {"network", "NETWORK-01-REROUTE_STD", "FOPR,FGPR,WBHP:PROD1,WBHP:PROD2", 1.0},
        {"network", "NETWORK-GSATPROD-GCONSUMP", "FOPR,FGPR,FGCR,FGSR,WBHP:PROD1", 1.5},
        {"model2", "9_3C_GINJ_GAS_GCONSUMP_STW", "FGPR,FGSR,FPR", 2.5},
    };
    for (const Ref& r : refs) {
        const std::string dir = root + "/" + r.dir + "/", deck = dir + r.deck + ".DATA", ref = dir + "opm-simulation-reference/flow/" + r.deck + ".SMSPEC";
        if (!exists(deck) || !exists(ref)) { std::printf("  [SKIP] %s (fetch opm-tests: python tools/fetch_opm_tests.py --dirs network,model2)\n", r.deck); continue; }
        const c35f::Run R = c35f::run(c35::readText(deck), dir);
        bool ok = R.ok;
        double worst = 0;
        if (ok) {
            ressim::ecl::Summary A; std::string err;
            ok = ressim::ecl::readSummary(ref, A, &err);
            if (ok) {
                const ressim::ecl::Summary Bs = ressim::ecl::summaryFromRun(R.sum, R.wells, R.units, {});
                std::vector<std::string> keys; { std::string k; for (const char* p = r.keys;; ++p) { if (*p == ',' || !*p) { keys.push_back(k); k.clear(); if (!*p) break; } else k += *p; } }
                const auto D = ressim::ecl::compareSummaries(A, Bs, keys);
                std::printf("    %-28s", r.deck);
                for (const auto& d : D) { std::printf(" %s %.2f%%", d.key.c_str(), 100 * d.relRms); worst = std::max(worst, 100 * d.relRms); }
                std::printf("\n");
                ok = D.size() == keys.size();
            }
        }
        char msg[160]; std::snprintf(msg, sizeof msg, "C35: %s vs OPM Flow - relRMS of the main vectors <= %.1f %%", r.deck, r.tol);
        check(ok && worst <= r.tol, msg);
    }
    {   // automatic choke (no OPM reference): the manifold node B1 holds the group ORAT target 6000 at the first report step
        const std::string dir = root + "/network/", deck = dir + "NETWORK_MODEL5_STDW_AUTOCHK.DATA";
        if (!exists(deck)) std::printf("  [SKIP] NETWORK_MODEL5_STDW_AUTOCHK (fetch opm-tests network)\n");
        else {
            std::string d = c35::readText(deck);
            const size_t i = d.find(" 1 FEB 2020 /");
            bool ok = i != std::string::npos;
            if (ok) {
                d = d.substr(0, i) + " 1 FEB 2020 /\n/\nEND\n";
                const c35f::Run R = c35f::run(d, dir);
                ok = R.ok;
                if (ok) {
                    double qb = 0;
                    const auto& p = R.sum.back();
                    for (size_t w = 0; w < p.wells.size() && w < R.wells.size(); ++w) if (R.wells[w].rfind("B-", 0) == 0) qb += p.wells[w].eff * p.wells[w].oil;
                    const double pB1 = c35f::nodeP(R, "B1"), pM = c35f::nodeP(R, "M5S");
                    std::printf("    AUTOCHK: B1 oil rate %.1f (target 6000), node B1 %.3f vs M5S %.3f (choke dp %.3f)\n", qb, pB1, pM, pB1 - pM);
                    ok = std::fabs(qb - 6000) < 0.02 * 6000 && pB1 > pM + 0.1;
                }
            }
            check(ok, "C35: automatic choke (NODEPROP AS_CHOKE) raises the manifold pressure until the group ORAT target is met");
        }
    }
}
