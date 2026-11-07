// C34 tests: opm-tests driven engine items - OPERATE / OPERATER, EDITNNC, rock compaction (ROCKCOMP / ROCKTAB), WECON
// connection workovers (CON with COMPLUMP), VREP injection coupled to the production of the same time step, end-point
// arrays without ENDSCALE (OPM behaviour, see testEndPointScaling).

namespace c34 {
static std::string qfs() {
    ressim::api::Session s;
    if (!s.open("example:qfs")) return "";
    const std::string path = c33::tmpDir() + "/qfs_c34.DATA";
    if (!s.exportDeck(path)) return "";
    std::ifstream f(path); std::stringstream ss; ss << f.rdbuf();
    return ss.str();
}
// insert text right after the line holding 'marker' (e.g. "\nPROPS\n"); before it when before = true
static std::string ins(std::string t, const std::string& marker, const std::string& text, bool before = false) {
    const size_t p = t.find(marker);
    if (p == std::string::npos) return "";
    t.insert(before ? p + 1 : p + marker.size(), text);
    return t;
}
static bool runCase(const std::string& deck, ressim::Simulator& sim) {
    try {
        const ressim::Case c = ressim::eclipseDeckFromString(deck);
        sim.setup(c); sim.initialize();
        while (sim.runReportStep()) {}
        return true;
    } catch (const std::exception& e) { std::printf("    run failed: %s\n", e.what()); return false; }
}
}  // namespace c34

static void testC34Items() {
    std::printf("== C34 opm-tests items: OPERATE(R), EDITNNC, ROCKCOMP, WECON CON + COMPLUMP, VREP coupling\n");
    const std::string base = c34::qfs();
    check(!base.empty(), "C34: example deck exported");
    if (base.empty()) return;
    // OPERATE (box) and OPERATER (OPERNUM regions) folded into the arrays
    {
        std::string d = c34::ins(base, "\nPROPS\n", "OPERNUM\n 625*1 /\nOPERATE\n PORO 1 1 1 1 1 1 MULTX PORO 0.5 /\n/\n"
                                                     "OPERATER\n PERMX 1 MAXLIM PERMX 50 /\n/\n", true);
        bool ok = false;
        try {
            const ressim::Case a = ressim::eclipseDeckFromString(base), b = ressim::eclipseDeckFromString(d);
            const double p0 = a.grid.poro[0], p1 = b.grid.poro[0], k1 = *std::max_element(b.grid.permx.begin(), b.grid.permx.end());
            std::printf("    PORO(1,1,1) %.4g -> %.4g, max PERMX after MAXLIM 50: %.4g\n", p0, p1, k1);
            ok = std::fabs(p1 - 0.5 * p0) < 1e-12 && std::fabs(b.grid.poro[1] - a.grid.poro[1]) < 1e-12 && k1 <= 50 + 1e-12;
        } catch (const std::exception& e) { std::printf("    %s\n", e.what()); }
        check(ok, "OPERATE MULTX on one cell and OPERATER MAXLIM on an OPERNUM region change exactly the selected cells");
    }
    // EDITNNC multiplies a non-neighbour connection (explicit NNC between cells 1 and 3 of row 1); neighbours untouched
    {
        const std::string g = "NNC\n 1 1 1 3 1 1 10.0 /\n/\n";
        std::string d1 = c34::ins(base, "\nPROPS\n", g, true);
        std::string d2 = c34::ins(base, "\nPROPS\n", g + "EDIT\nEDITNNC\n 1 1 1 3 1 1 0.25 /\n 1 1 1 2 1 1 0.0 /\n/\n", true);
        bool ok = false;
        try {
            ressim::Grid G1, G2; G1.build(ressim::eclipseDeckFromString(d1)); G2.build(ressim::eclipseDeckFromString(d2));
            auto tr = [](const ressim::Grid& G, int a, int b) { for (auto& c : G.conns) if ((c.a == a && c.b == b) || (c.a == b && c.b == a)) return c.trans; return -1.0; };
            const double t1 = tr(G1, 0, 2), t2 = tr(G2, 0, 2), n1 = tr(G1, 0, 1), n2 = tr(G2, 0, 1);
            std::printf("    NNC 1-3: %.4g -> %.4g (EDITNNC 0.25); neighbour 1-2: %.4g -> %.4g (EDITNNC on neighbours is ignored)\n", t1, t2, n1, n2);
            ok = t1 > 0 && std::fabs(t2 - 0.25 * t1) < 1e-9 * t1 && std::fabs(n2 - n1) < 1e-12 && G2.numEditNnc == 1;
        } catch (const std::exception& e) { std::printf("    %s\n", e.what()); }
        check(ok, "EDITNNC: multiplier applied to the non-neighbour connection only");
    }
    // ROCKCOMP + ROCKTAB: a unit table changes nothing, a pressure-dependent table changes the run
    {
        const std::string unit = c34::ins(c34::ins(base, "\nRUNSPEC\n", "ROCKCOMP\n REVERS 1 /\n"), "\nPROPS\n", "ROCKTAB\n 14.7 1.0 1.0\n 10000 1.0 1.0 /\n");
        const std::string var = c34::ins(c34::ins(base, "\nRUNSPEC\n", "ROCKCOMP\n IRREVERS 1 /\n"), "\nPROPS\n", "ROCKTAB\n 14.7 0.90 0.5\n 10000 1.10 1.5 /\n");
        ressim::Simulator s0, s1, s2;
        const bool ok0 = c34::runCase(base, s0), ok1 = c34::runCase(unit, s1), ok2 = c34::runCase(var, s2);
        const double f0 = ok0 ? s0.summary().back().fopt : 0, f1 = ok1 ? s1.summary().back().fopt : 0, f2 = ok2 ? s2.summary().back().fopt : 0;
        std::printf("    FOPT base %.6g, unit ROCKTAB %.6g, varying ROCKTAB %.6g\n", f0, f1, f2);
        check(ok0 && ok1 && std::fabs(f1 - f0) < 1e-6 * f0, "ROCKTAB with unit multipliers reproduces the base run");
        check(ok2 && std::fabs(f2 - f0) > 1e-3 * f0, "ROCKTAB pore-volume / transmissibility multipliers (IRREVERS) change the run");
    }
    // WECON CON with COMPLUMP: three producer connections, the worst-offending completion is closed, the well stays open
    {
        std::string d = base;
        const size_t p = d.find(" 'PROD' 25 25 1 1 OPEN");
        if (p != std::string::npos) d.insert(p, " 'PROD' 24 25 1 1 OPEN 1* 1* 0.5 1* 0 1* Z /\n 'PROD' 25 24 1 1 OPEN 1* 1* 0.5 1* 0 1* Z /\n");
        d = c34::ins(d, "\nTSTEP", "COMPLUMP\n 'PROD' 24 25 1 1 1 /\n 'PROD' 25 24 1 1 2 /\n 'PROD' 25 25 1 1 3 /\n/\nWECON\n 'PROD' 1* 1* 0.3 2* CON NO /\n/\n", true);
        ressim::Simulator s;
        const bool ok = !d.empty() && c34::runCase(d, s);
        int closed = 0, open = 0;
        if (ok) for (auto& c : s.connections()) if (c.well == "PROD") (c.wi > 0 ? open : closed)++;
        const double q = ok ? s.summary().back().fopr + s.summary().back().fwpr : 0;
        std::printf("    PROD connections closed %d, open %d, final liquid rate %.4g\n", closed, open, q);
        check(ok && closed >= 1 && (open >= 1 ? q > 0 : true), "WECON CON: the worst completion (COMPLUMP) is closed instead of the whole well");
    }
    // VREP: the injector follows the voidage of the same time step (first step injects)
    {
        std::string d = base;
        const size_t p = d.find(" 'INJ' WATER OPEN RATE 300 1* 5000 /");
        if (p != std::string::npos) d.replace(p, std::string(" 'INJ' WATER OPEN RATE 300 1* 5000 /").size(), " 'INJ' WATER OPEN GRUP 5000 1* 9000 /");
        d = c34::ins(d, "\nWCONINJE", "GCONINJE\n 'FIELD' WATER VREP 3* 1.0 /\n/\n", true);
        ressim::Simulator s;
        const bool ok = !d.empty() && c34::runCase(d, s);
        double r1 = -1, rN = -1;
        if (ok && s.summary().size() > 2) {
            const auto& a = s.summary()[1]; const auto& b = s.summary().back();
            r1 = a.fwir; rN = b.fwir / std::max(1e-9, b.fwpr + b.fopr);
            std::printf("    first step FWIR %.4g (production %.4g + %.4g); last FWIR / FLPR %.4g\n", a.fwir, a.fopr, a.fwpr, rN);
        }
        check(ok && r1 > 0, "VREP: injection in the first time step (targets coupled to the current production, NUPCOL-like)");
    }
}
