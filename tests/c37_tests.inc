// C37 tests: MSW tubing flash, oil-water slip, segment storage, WSEGFLIM, WSEGPULL; PYACTION
#include "ressim/pyaction.h"
namespace c37m {
// SPE1 vertical MSW to surface (c36m::spe1Vertical) with extra SCHEDULE text inserted before WCONPROD / RESSIM lines appended
std::string spe1Msw(const std::string& model, const std::string& extraSched, const std::string& ressimLines = "") {
    std::string t = c36m::spe1Vertical(model);
    if (t.empty()) return t;
    if (!extraSched.empty()) {
        const size_t p = t.find("WCONPROD");
        if (p == std::string::npos) return "";
        t.insert(p, extraSched + "\n");
    }
    if (!ressimLines.empty()) {
        const size_t p = t.find("SCHEDULE");
        if (p == std::string::npos) return "";
        t.insert(p, ressimLines + "\n");
    }
    return t;
}
struct Run { double fopt = 0, maxFopr = 0, medFopr = 0, pBot = 0, colDp = 0; std::vector<ressim::Simulator::MswRecord> rec; std::string log; };
Run run(const std::string& text, int stopAfter = -1, std::string* blob = nullptr) {
    Run r;
    std::vector<ressim::DeckMessage> rep;
    const ressim::Case c = ressim::eclipseDeckFromString(text, "", nullptr, &rep);
    ressim::Simulator s; s.log = [&r](const std::string& m) { r.log += m + "\n"; }; s.setup(c); s.initialize();
    while ((stopAfter < 0 || s.reportStep() < stopAfter) && s.runReportStep()) {}
    if (blob) *blob = s.saveState();
    if (s.summary().empty()) return r;
    r.fopt = s.summary().back().fopt;
    for (const auto& p : s.summary()) if (p.time > 60) r.maxFopr = std::max(r.maxFopr, p.fopr);   // (after the first lagged steps)
    r.rec = s.mswRecords();
    int nc = 0;
    for (const auto& q : r.rec) if (q.t > 30 && q.t < 300 && q.p.size() > 1) { r.colDp += q.p.back() - q.p[0]; ++nc; }
    if (nc) r.colDp /= nc;
    std::vector<double> rates;
    for (const auto& p : s.summary()) if (p.time > 60) rates.push_back(p.fopr);
    if (!rates.empty()) { std::sort(rates.begin(), rates.end()); r.medFopr = rates[rates.size() / 2]; }
    if (!r.rec.empty()) { r.pBot = r.rec.back().p.back(); }
    return r;
}
}  // namespace c37m

static void testC37Msw() {
    std::printf("== C37 multisegment wells: tubing flash, oil-water slip, segment storage, WSEGFLIM, WSEGPULL\n");
    using namespace ressim::msw;
    const double PI = 3.14159265358979323846;
    ressim::msw::Units u; u.g = 9.80665e-5; u.toM = 1; u.toRho = 1; u.toQ = 1.0 / 86400.0; u.toPa = 1e5;
    // ---- parser
    {
        std::vector<ExtraKeyword> kw(4);
        kw[0].name = "WELSEGS"; kw[0].body = " 'W' 1000 0 1* INC 'HF-' 'DF' /\n 2 3 1 1 100 50 0.1 1E-4 1* 0.75 /\n/\n";
        kw[1].name = "WSEGFLIM"; kw[1].body = " 'W' 3 3 OIL 500 12 WAT 200 1* /\n/\n";
        kw[2].name = "WSEGPULL"; kw[2].body = " 'W' 2 2 0.9 25 300 0.05 /\n/\n";
        kw[3].name = "WSEGDFPA"; kw[3].body = " 1* 1* 1* 1* 1* 1* 0.8 1.1 0.3 0.8 1.5 /\n";
        std::string err;
        const auto W = parse(kw, 0, false, err);
        bool ok = W.size() == 1 && err.empty();
        if (ok) {
            const Segment& g3 = W[0].segs[size_t(W[0].indexOf(3))];
            const Segment& g2 = W[0].segs[size_t(W[0].indexOf(2))];
            ok = std::fabs(g3.volume - 0.75) < 1e-12 && g3.flimPh[0] == 1 && std::fabs(g3.flimQ[0] - 500) < 1e-12 && std::fabs(g3.flimA[0] - 12) < 1e-12 &&
                 g3.flimPh[1] == 2 && std::fabs(g3.flimQ[1] - 200) < 1e-12 && std::fabs(g3.flimA[1] - 8) < 1e-12 &&
                 g2.pump && std::fabs(g2.pumpP - 0.9) < 1e-12 && std::fabs(g2.pumpDp - 25) < 1e-12 && std::fabs(g2.pumpQo - 300) < 1e-12 && std::fabs(g2.pumpG - 0.05) < 1e-12 &&
                 std::fabs(W[0].df.owFactor - 0.8) < 1e-12 && std::fabs(W[0].df.owA - 1.1) < 1e-12 && std::fabs(W[0].df.owB1 - 0.3) < 1e-12 &&
                 std::fabs(W[0].df.owB2 - 0.8) < 1e-12 && std::fabs(W[0].df.owN - 1.5) < 1e-12 && std::fabs(W[0].df.A - 1.2) < 1e-12;
        }
        check(ok, "parser: WELSEGS item 10 volume, WSEGFLIM (2 limits, default exponent), WSEGPULL items, WSEGDFPA oil-water items 7-11");
    }
    // ---- oil-water holdup limits
    {
        DriftParams P; P.owFactor = 1.0;
        const double ro = 750, rw = 1020;
        const double Vc = std::pow(P.owSigma * 9.80665 * (rw - ro) / (rw * rw), 0.25);
        double C0 = 0, vd = 0;
        const double a1 = owHoldup(0.02, 0.5, ro, rw, 0.0, P, &C0, &vd);   // little oil, upward: oil rises faster
        const bool b1 = a1 < 0.02 / 0.52 && std::fabs(a1 * (C0 * 0.52 + vd) - 0.02) < 1e-9 && std::fabs(C0 - 1.2) < 1e-9 && vd > 0.9 * 1.53 * Vc;
        const double a2 = owHoldup(0.3, 0.3, ro, rw, PI / 2, P, &C0, &vd);   // horizontal: no drift
        const bool b2 = std::fabs(vd) < 1e-12 && std::fabs(a2 * C0 * 0.6 - 0.3) < 1e-9;
        const double a3 = owHoldup(0.1, 0.5, ro, rw, PI, P, &C0, &vd);       // downward: oil lags
        const bool b3 = vd < 0 && a3 > 0.1 / (C0 * 0.6);
        DriftParams P0 = P; P0.owFactor = 1.0;
        const double a4 = owHoldup(0.2, 0.2, 1000, 1000, 0.0, P0, &C0, &vd);  // equal densities: no drift
        const bool b4 = std::fabs(vd) < 1e-12;
        std::printf("    oil-water: bubbly up alpha_o %.5f (no-slip %.5f, v_d %.4f, 1.53 V_c %.4f), horizontal %.4f, downward %.4f (no-slip %.4f)\n",
                    a1, 0.02 / 0.52, vd, 1.53 * Vc, a2, a3, 0.1 / 0.6);
        check(b1 && b2 && b3 && b4 && a4 > 0, "oil-water drift flux (Shi et al. 2005): oil rises in upflow, no drift horizontally / at equal densities, oil lags in downflow");
    }
    // ---- flow limiter
    {
        const double f1 = flowLimitFactor(10, 1000, 8), f2 = flowLimitFactor(5000, 1000, 8), f3 = flowLimitFactor(5000, 1000, 50), f4 = flowLimitFactor(1000, 1000, 8);
        bool mono = true; double last = 0;
        for (double q0 = 1; q0 < 1e5; q0 *= 1.3) { const double q = q0 * flowLimitFactor(q0, 1000, 8); mono = mono && q >= last - 1e-9 && q <= 1000 * (1 + 1e-12); last = q; }
        std::printf("    limiter: f(0.01 Qlim) %.6f, f(5 Qlim) %.5f (a 8) %.5f (a 50), f(Qlim) %.5f\n", f1, f2, f3, f4);
        check(f1 > 0.999999 && std::fabs(f2 * 5 - 1) < 1e-3 && std::fabs(f3 * 5 - 1) < 1e-9 && std::fabs(f4 - std::pow(2.0, -1.0 / 8)) < 1e-12 && mono,
              "WSEGFLIM limiter: f = (1 + (Q0/Qlim)^a)^(-1/a) - no effect far below the limit, limited rate monotone and <= the limit");
    }
    // ---- tubing flash in the segment solve (synthetic PVT: Rs = 0.1 p, b_o = 1, gas b = p / 100, rho_g = 0.8 p / 100)
    {
        Well w; w.name = "F";
        Segment t; t.id = 1; t.depth = 0; t.diam = 0.1; w.segs.push_back(t);
        for (int i = 2; i <= 11; ++i) { Segment s; s.id = i; s.outlet = i - 1; s.length = 100.0 * (i - 1); s.depth = 100.0 * (i - 1); s.diam = 0.1; w.segs.push_back(s); }
        Inflow in; in.resize(11);
        in.q[10] = 300; in.rho[10] = 800; in.mu[10] = 1; in.a[10] = 10; in.pres[10] = 250;
        in.qs.assign(33, 0.0); in.qs[30] = 300; in.qs[32] = 300 * 20.0;   // oil 300 sm3/d, gas 6000 sm3/d (Rs 20 -> bubble point 200 bar)
        in.flash = [](double p, const double* q) {
            Inflow::Flash F;
            const double rs = std::min(q[2] / std::max(q[0], 1e-12), 0.1 * p);
            const double fg = std::max(0.0, q[2] - rs * q[0]), bg = p / 100.0;
            F.qo = q[0]; F.rhoO = 800 - 2 * (20 - rs); F.muO = 1;
            F.qg = fg > 0 ? fg / bg : 0.0; F.rhoG = 0.8 * bg; F.muG = 0.02;
            return F;
        };
        in.flashAtSegment = true;
        const Solution S = solve(w, 150.0, in, u);
        Inflow inR = in; inR.flashAtSegment = false;   // fluid at the reservoir pressure (250 bar: all gas dissolved)
        const Solution SR = solve(w, 150.0, inR, u);
        // free gas fraction rises toward the top, none at the bottom where p > 150 bar
        const double fTop = S.noSlip[1], fBot = S.noSlip[10];
        std::printf("    flash: top-node gas fraction %.3f, bottom %.3f; p(toe) %.2f bar vs %.2f bar at reservoir conditions; q top %.1f vs %.1f rm3/d\n",
                    fTop, fBot, S.p[10], SR.p[10], S.qSeg[1], SR.qSeg[1]);
        check(fTop > 0.2 && fBot < 0.5 * fTop && SR.noSlip[1] < 1e-12 && S.p[10] < SR.p[10] - 5 && S.qSeg[1] > 1.5 * SR.qSeg[1] && S.qgLoc[1] > 0,
              "tubing flash: gas released up the tubing below the bubble point -> lighter column, larger local rate");
    }
    // ---- oil-water slip in a vertical column (flash provides the oil / water split)
    {
        Well w; w.name = "OW"; w.driftFlux = true; w.df.owFactor = 1.0;
        Segment t; t.id = 1; t.depth = 0; t.diam = 0.15; w.segs.push_back(t);
        Segment s; s.id = 2; s.outlet = 1; s.length = 1000; s.depth = 1000; s.diam = 0.15; w.segs.push_back(s);
        Inflow in; in.resize(2); in.q[1] = 200; in.rho[1] = 0.3 * 750 + 0.7 * 1020; in.mu[1] = 1; in.qs = {0, 0, 0, 60, 140, 0};
        in.flash = [](double, const double* q) { Inflow::Flash F; F.qo = q[0]; F.qw = q[1]; F.rhoO = 750; F.rhoW = 1020; return F; };
        const Solution S = solve(w, 50.0, in, u);
        Well w0 = w; w0.df.owFactor = -1;
        const Solution S0 = solve(w0, 50.0, in, u);
        std::printf("    oil-water column: oil holdup %.4f vs no-slip %.4f, head density %.1f vs %.1f kg/m3\n", S.oilHold[1], S0.oilHold[1], S.rhoHead[1], S0.rhoHead[1]);
        check(S.oilHold[1] < 0.3 - 0.01 && std::fabs(S0.oilHold[1] - 0.3) < 1e-12 && S.rhoHead[1] > S0.rhoHead[1] + 2 &&
              std::fabs(S.rhoHead[1] - (S.oilHold[1] * 750 + (1 - S.oilHold[1]) * 1020)) < 1e-9,
              "oil-water slip (WSEGDFPA item 7 > 0): faster oil leaves more water in an upward column -> heavier head");
    }
    // ---- pump and flow limiter in the solve
    {
        Well w; w.name = "P";
        Segment t; t.id = 1; t.depth = 1000; t.diam = 0.1; w.segs.push_back(t);
        Segment s; s.id = 2; s.outlet = 1; s.length = 100; s.depth = 1000; s.diam = 0.1; s.pump = true; s.pumpP = 0.8; s.pumpDp = 20; s.pumpQo = 100; s.pumpG = 0.01; w.segs.push_back(s);
        Segment d; d.id = 3; d.outlet = 2; d.length = 101; d.depth = 1000; d.diam = 0.1; d.flimPh[0] = 1; d.flimQ[0] = 150; d.flimA[0] = 8; w.segs.push_back(d);
        Inflow in; in.resize(3); in.q[2] = 400; in.rho[2] = 800; in.mu[2] = 1; in.a[2] = 5; in.pres[2] = 200;
        in.qs = {0, 0, 0, 0, 0, 0, 300, 50, 0};
        in.flash = [](double, const double* q) { Inflow::Flash F; F.qo = q[0] * 1.2; F.qw = q[1]; F.rhoO = 780; F.rhoW = 1000; return F; };
        Well w0 = w; w0.segs[1].pump = false; w0.segs[2].flimPh[0] = 0;
        const Solution S = solve(w, 100.0, in, u), S0 = solve(w0, 100.0, in, u);
        const double boost = 0.8 * (20 - 0.01 * (300 - 100));
        const double D0 = 200 - S.p[2];
        const double qres0 = 5 * D0, qo0 = 300.0 / 400.0 * qres0, f = flowLimitFactor(qo0, 150, 8);
        std::printf("    pump boost %.4f bar (hand %.4f); limiter factor %.5f (hand %.5f, potential oil %.1f -> %.1f sm3/d)\n", -S.dPump[1], boost, S.flimFac[2], f, qo0, qo0 * f);
        check(std::fabs(-S.dPump[1] - boost) < 1e-9 && std::fabs(S.p[1] - S0.p[1] + boost) < 1e-9 && std::fabs(S.flimFac[2] - f) < 1e-12 &&
              std::fabs(S.devFac[2] - f) < 1e-12 && qo0 * f <= 150 && S0.devFac[2] == 1.0,
              "WSEGPULL boost P (DP - G (q_o - q_lim)) lowers the upstream pressure; WSEGFLIM factor on the subtree conductance");
    }
    // ---- segment storage
    {
        Well w; w.name = "S";
        Segment t; t.id = 1; t.depth = 0; t.diam = 0.1; w.segs.push_back(t);
        Segment s; s.id = 2; s.outlet = 1; s.length = 500; s.depth = 500; s.diam = 0.1; w.segs.push_back(s);
        Inflow in; in.resize(2); in.q[1] = 50; in.rho[1] = 700; in.mu[1] = 1;
        std::vector<double> sto = {0.0, 1000.0, 0.0, 0.0};   // previous content: water column
        SolveOpts o; o.dt = 0.25; o.sto = &sto;
        const Solution S = solve(w, 10.0, in, u, {}, o);
        const double V = 0.25 * PI * 0.01 * 500, sw = 50.0 * 0.25, wgt = sw / (V + sw);
        const double hand = 1000 + wgt * (700 - 1000);
        Inflow i0 = in; i0.q[1] = 0; i0.rho[1] = 700;
        const Solution Sz = solve(w, 10.0, i0, u, {}, o);
        std::printf("    storage: V %.3f m3, swept %.2f m3 -> head density %.2f (hand %.2f); shut segment keeps %.1f kg/m3\n", V, sw, S.rhoHead[1], hand, Sz.rhoHead[1]);
        check(std::fabs(S.rhoHead[1] - hand) < 1e-9 && std::fabs(Sz.rhoHead[1] - 1000) < 1e-9, "segment storage: content moves toward the inflow with w = q dt / (V + q dt); a shut segment keeps its column");
    }
    // ---- deck runs: SPE1 vertical MSW to surface
    {
        const std::string base = c37m::spe1Msw("HO", "");
        if (!base.empty()) {
            const std::string off = c37m::spe1Msw("HO", "", "-- RESSIM: MSWOPTS 0 1");
            const c37m::Run RF = c37m::run(base), RO = c37m::run(off);
            double topF = 0, topO = 0;
            for (const auto& r : RF.rec) if (r.t > 2 && r.t < 400 && r.noSlip.size() > 2) topF = std::max(topF, r.noSlip[1]);
            for (const auto& r : RO.rec) if (r.t > 2 && r.t < 400 && r.noSlip.size() > 2) topO = std::max(topO, r.noSlip[1]);
            std::printf("    SPE1 tubing to surface: flash FOPT %.6g (column drop %.1f psi at 30-300 d, top gas fraction %.3f) | no flash FOPT %.6g (%.1f psi, %.3f)\n",
                        RF.fopt, RF.colDp, topF, RO.fopt, RO.colDp, topO);
            check(topF > 0.2 && topO < 0.05 && RF.colDp < RO.colDp - 200,
                  "SPE1 tubing flash: gas released below the bubble point lightens the tubing column (bottom-to-top drop)");
            // WSEGFLIM oil limit on the bottom segment
            const std::string lim = c37m::spe1Msw("HO", "WSEGFLIM\n 'PROD' 17 17 OIL 8000 20 /\n/\n");
            const c37m::Run RL = c37m::run(lim);
            std::printf("    WSEGFLIM OIL 8000: FOPR after 60 d median %.1f max %.1f (unlimited max %.1f), FOPT %.6g vs %.6g\n", RL.medFopr, RL.maxFopr, RF.maxFopr, RL.fopt, RF.fopt);
            check(RF.maxFopr > 9000 && std::fabs(RL.medFopr - 8000) < 0.02 * 8000 && RL.maxFopr < 8000 * 1.12 && RL.fopt < RF.fopt,
                  "WSEGFLIM on SPE1: the flow-limiting valve holds the oil rate at the 8000 stb/d limit (median within 2 %, lagged transients < 12 %)");
            // WSEGPULL pump near the bottom
            const std::string pmp = c37m::spe1Msw("HO", "WSEGPULL\n 'PROD' 16 16 1.0 800 /\n/\n");
            const c37m::Run RP = c37m::run(pmp);
            std::printf("    WSEGPULL 800 psi: FOPT %.6g vs %.6g, bottom node %.1f vs %.1f psi\n", RP.fopt, RF.fopt, RP.pBot, RF.pBot);
            check(RP.fopt > RF.fopt * 1.005 && RP.pBot < RF.pBot - 300, "WSEGPULL on SPE1: the pump boost lowers the bottom-hole pressure -> more oil");
            // storage on / off and restart continuity
            const std::string noSto = c37m::spe1Msw("HO", "", "-- RESSIM: MSWOPTS 1 0");
            const c37m::Run RS = c37m::run(noSto);
            std::printf("    storage on FOPT %.8g, off %.8g (rel. diff %.2e)\n", RF.fopt, RS.fopt, std::fabs(RF.fopt - RS.fopt) / RF.fopt);
            check(std::fabs(RF.fopt - RS.fopt) < 5e-3 * RF.fopt, "segment storage: negligible for a flowing well at daily steps (< 0.5 % FOPT)");
            std::string blob;
            c37m::run(base, 5, &blob);
            std::vector<ressim::DeckMessage> rep;
            const ressim::Case c = ressim::eclipseDeckFromString(base, "", nullptr, &rep);
            ressim::Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
            const bool loaded = !blob.empty();
            if (loaded) s.loadState(blob);
            while (s.runReportStep()) {}
            const double fr = s.summary().empty() ? 0.0 : s.summary().back().fopt;
            std::printf("    restart from report step 5: FOPT %.10g vs continuous %.10g\n", fr, RF.fopt);
            check(loaded && std::fabs(fr - RF.fopt) < 1e-6 * RF.fopt, "restart: segment storage state carried in the well state (continuous = restarted)");
            // export round trip of the options and verbatim keywords
            const ressim::Case co = ressim::eclipseDeckFromString(c37m::spe1Msw("HO", "WSEGFLIM\n 'PROD' 17 17 OIL 8000 20 /\n/\n", "-- RESSIM: MSWOPTS 0 0"), "", nullptr, &rep);
            const std::string txt = ressim::eclipseDeckToString(co);
            const ressim::Case c2 = ressim::eclipseDeckFromString(txt, "", nullptr, &rep);
            check(!co.mswOpt.flash && !co.mswOpt.storage && txt.find("-- RESSIM: MSWOPTS 0 0") != std::string::npos && !c2.mswOpt.flash && !c2.mswOpt.storage &&
                  txt.find("WSEGFLIM") != std::string::npos, "export: MSWOPTS written (non-default) and re-imported, WSEGFLIM kept verbatim");
        }
    }
}

// ============================================================================ C37-3: PYACTION
namespace c37p {
namespace fs = std::filesystem;
const char* kCut = "WCONPROD\n 'PROD' OPEN ORAT 5000 1* 1* 1* 1* 1000 /\n/\n";
// SPE1 with extra SCHEDULE text before its TSTEP block
std::string spe1With(const std::string& extra) {
    std::ifstream f(findExample("SPE1.DATA"));
    std::stringstream ss; ss << f.rdbuf();
    std::string t = ss.str();
    const size_t p = t.find("TSTEP", t.find("WCONINJE"));
    if (p == std::string::npos) return "";
    t.insert(p, extra + "\n");
    return t;
}
struct PRun { double fopt = 0, foprEnd = 0; std::string log; size_t events = 0; std::string blob; int calls = 0; };
PRun run(const std::string& text, const std::string& dir, int stopAfter = -1, const std::string* restart = nullptr) {
    PRun r;
    std::vector<ressim::DeckMessage> rep;
    const ressim::Case c = ressim::eclipseDeckFromString(text, dir, nullptr, &rep);
    ressim::Simulator s; s.log = [&r](const std::string& m) { r.log += m + "\n"; }; s.setup(c); s.initialize();
    if (restart) s.loadState(*restart);
    while ((stopAfter < 0 || s.reportStep() < stopAfter) && s.runReportStep()) {}
    if (stopAfter >= 0) r.blob = s.saveState();
    if (s.summary().empty()) return r;
    r.fopt = s.summary().back().fopt; r.foprEnd = s.summary().back().fopr;
    r.events = s.actionEvents().size();
    return r;
}
size_t count(const std::string& h, const std::string& n) { size_t k = 0, p = 0; while ((p = h.find(n, p)) != std::string::npos) { ++k; p += n.size(); } return k; }
}  // namespace c37p

static void testC37Py() {
    std::printf("== C37 PYACTION: JSON protocol, opm_embedded shim, schedule / ACTIONX callback, restart\n");
    using namespace ressim::pyact;
    namespace fs = std::filesystem;
    // ---- JSON
    {
        Json o = Json::obj();
        Json a = Json::arr(); a.a = {Json::num(1.5), Json::boolean(true), Json(), Json::str("x\"y\\z\n\t\x01")};
        o.set("a", a).set("u", Json::str("\xC3\xA9")).set("n", Json::num(-2.5e-12)).set("nan", Json::num(std::nan("")));
        const std::string d = dump(o);
        Json b; std::string err;
        const bool ok = parse(d, b, err);
        Json u; const bool ok2 = parse("{\"s\":\"\\u00e9\\ud83d\\ude00\", \"v\": [1e3, -0.5, NaN]}", u, err);
        check(ok && b.get("a") && b.get("a")->a.size() == 4 && b.get("a")->a[3].s == "x\"y\\z\n\t\x01" && b.get("u")->s == "\xC3\xA9" &&
              b.get("n")->n == -2.5e-12 && b.get("nan")->t == Json::Null && ok2 && u.get("s")->s == "\xC3\xA9\xF0\x9F\x98\x80" &&
              u.get("v")->a[0].n == 1000 && u.get("v")->a[2].t == Json::Null,
              "JSON writer / parser round trip (escapes, UTF-8, surrogate pairs, NaN -> null)");
    }
    // ---- PYACTION records
    {
        std::string err;
        auto A = parseActions({{0, " 'CTRL' 'UNLIMITED' /\n 'ctrl.py' /\n"}, {2, " WELLCONTROLLER FIRST_TRUE / -- comment\n 'dir/wc.py' /\n"},
                               {3, " 'ONE' /\n 'one.py' /\n"}, {4, " 'BAD' 'SOMETIMES' /\n 'b.py' /\n"}}, err);
        check(A.size() == 4 && A[0].runCount == "UNLIMITED" && A[0].file == "ctrl.py" && A[1].name == "WELLCONTROLLER" && A[1].runCount == "FIRST_TRUE" &&
              A[1].file == "dir/wc.py" && A[1].fromStep == 2 && A[2].runCount == "SINGLE" && A[3].runCount == "SINGLE" && err.find("SOMETIMES") != std::string::npos,
              "PYACTION records: name / run count (SINGLE default, unknown -> SINGLE + warning), script file, report step");
    }
    // ---- the Python child
    Host probe;
    std::string perr;
    if (!probe.start("", perr)) {
        std::printf("    (no Python 3 interpreter: %s) - PYACTION run tests skipped\n", perr.c_str());
        check(perr.find("Python") != std::string::npos, "PYACTION without Python: clear error, run tests skipped");
        return;
    }
    std::printf("    Python %s\n", probe.version().c_str());
    const fs::path dir = fs::temp_directory_path() / "ressim_c37py";
    fs::create_directories(dir);
    auto write = [&](const std::string& name, const std::string& body) { std::ofstream f(dir / name); f << body; };
    {   // direct protocol call with queries
        write("probe.py", "import opm_embedded\nss = opm_embedded.current_summary_state\nsch = opm_embedded.current_schedule\n"
                          "if 'n' not in locals():\n    n = 0\nn += 1\nprint('call', n, ss['FOPR'], ss.well_var('P1', 'WBHP'), 'FXYZ' in ss, ss.wells('P*'))\n"
                          "sch.shut_well('P1')\nif n == 2:\n    raise ValueError('boom')\n");
        Json req = Json::obj();
        Json wl = Json::arr(); wl.a = {Json::str("P1"), Json::str("I1")};
        req.set("cmd", Json::str("run")).set("name", Json::str("PROBE")).set("file", Json::str((dir / "probe.py").string())).set("step", Json::num(3)).set("wells", wl);
        int nq = 0;
        auto q = [&](const Json& x) -> Json {
            ++nq;
            if (x.a[1].s == "FOPR") return Json::num(1234.5);
            if (x.a[1].s == "WBHP" && x.a[2].s == "P1") return Json::num(250);
            return Json();
        };
        Json r1, r2; std::string e1, e2;
        const bool c1 = probe.call(req, q, r1, e1), c2 = probe.call(req, q, r2, e2);
        const std::string l1 = (r1.get("log") && !r1.get("log")->a.empty()) ? r1.get("log")->a[0].a[1].s : "";
        std::printf("    call 1: '%s', keywords %s| call 2 error: %s\n", l1.c_str(), r1.get("keywords") ? r1.get("keywords")->s.substr(0, 30).c_str() : "-",
                    (r2.get("error") ? r2.get("error")->s.substr(r2.get("error")->s.rfind("ValueError")).c_str() : "-"));
        check(c1 && c2 && r1.get("ok")->b && l1 == "call 1 1234.5 250 False ['P1']" && r1.get("keywords")->s.find("WELOPEN\n 'P1' SHUT /") == 0 &&
              !r2.get("ok")->b && r2.get("error")->s.find("boom") != std::string::npos && nq == 6,
              "Python host: persistent globals, summary queries answered by the engine, print() -> log, schedule shim, traceback on errors");
        probe.stop();
        check(!probe.running(), "Python host stops cleanly");
    }
    // ---- SPE1: PYACTION (new-style script) vs the equivalent ACTIONX
    write("cut.py", std::string("import opm_embedded\nss = opm_embedded.current_summary_state\nsch = opm_embedded.current_schedule\n"
                                "if 'setup_done' not in locals():\n    setup_done = True\n    done = False\n"
                                "if not done and ss['FOPT'] > 1e7:\n    done = True\n"
                                "    opm_embedded.OpmLog.info('rate cut at %.1f days, WBHP %.1f' % (ss.elapsed() / 86400, ss.well_var('PROD', 'WBHP')))\n"
                                "    sch.insert_keywords(\"\"\"") + c37p::kCut + "\"\"\")\n");
    write("legacy.py", "state = {'n': 0}\n"
                       "def run(ecl_state, schedule, report_step, summary_state, actionx_callback):\n"
                       "    state['n'] += 1\n    print('legacy call', state['n'], 'step', report_step)\n"
                       "    if summary_state.well_var('PROD', 'WOPT') > 1e7:\n        actionx_callback('CUTX', ['PROD'])\n        return True\n    return False\n");
    write("broken.py", "import opm_embedded\nraise RuntimeError('broken script')\n");
    const std::string base = c37p::spe1With("");
    const std::string tAx = c37p::spe1With(std::string("ACTIONX\n 'CUTA' 1 /\n FOPT > 1E7 /\n/\n") + c37p::kCut + "ENDACTIO\n");
    const std::string tPy = c37p::spe1With("PYACTION\n 'CUT' 'UNLIMITED' /\n 'cut.py' /\n");
    const std::string tLeg = c37p::spe1With(std::string("ACTIONX\n 'CUTX' 1 /\n FOPT < 0 /\n/\n") + c37p::kCut + "ENDACTIO\n\nPYACTION\n LEG FIRST_TRUE /\n 'legacy.py' /\n");
    const std::string tBad = c37p::spe1With("PYACTION\n 'BAD' 'UNLIMITED' /\n 'broken.py' /\n");
    if (base.empty()) { check(false, "SPE1 deck found"); return; }
    const c37p::PRun R0 = c37p::run(base, dir.string()), RA = c37p::run(tAx, dir.string()), RP = c37p::run(tPy, dir.string()),
                     RL = c37p::run(tLeg, dir.string()), RB = c37p::run(tBad, dir.string());
    std::printf("    FOPT base %.6g | ACTIONX %.6g | PYACTION %.6g | legacy run() + callback %.6g | broken script %.6g\n", R0.fopt, RA.fopt, RP.fopt, RL.fopt, RB.fopt);
    std::printf("    ACTIONX %.4f vs PYACTION %.4f, applies %zu\n", RA.fopt, RP.fopt, c37p::count(RP.log, "PYACTION CUT applies"));
    std::printf("    final FOPR: base %.1f, PYACTION %.1f; legacy calls %zu, applies %zu\n", R0.foprEnd, RP.foprEnd, c37p::count(RL.log, "legacy call"),
                c37p::count(RL.log, "PYACTION LEG applies"));
    {
        const size_t a = RP.log.find("rate cut at");
        std::printf("    log: %s\n", a == std::string::npos ? "-" : RP.log.substr(a, RP.log.find('\n', a) - a).c_str());
    }
    check(RA.fopt < 0.9 * R0.fopt && std::fabs(RP.fopt - RA.fopt) < 1e-9 * RA.fopt && RP.foprEnd <= 5000 * (1 + 1e-4) &&
          c37p::count(RP.log, "PYACTION CUT applies") == 1 && RP.log.find("rate cut at") != std::string::npos,
          "SPE1: PYACTION insert_keywords acts exactly like the equivalent ACTIONX (same report-step boundary, bitwise FOPT)");
    check(std::fabs(RL.fopt - RA.fopt) < 1e-9 * RA.fopt && c37p::count(RL.log, "PYACTION LEG applies") == 1 && RL.log.find("legacy call 1 step 1") != std::string::npos,
          "legacy run(): actionx_callback applies the named ACTIONX body (module body executed once, run() every call)");
    {
        const size_t nCalls = c37p::count(RL.log, "legacy call");
        const size_t lastCall = RL.log.rfind("legacy call"), apply = RL.log.find("PYACTION LEG applies");
        check(nCalls >= 2 && apply != std::string::npos && lastCall < apply, "FIRST_TRUE: the script is not called again after it returned True");
    }
    check(std::fabs(RB.fopt - R0.fopt) < 1e-12 * R0.fopt && RB.log.find("broken script") != std::string::npos && c37p::count(RB.log, "Python error") == 1,
          "a failing script: traceback in the log, run continues unchanged, not called again");
    // ---- restart across the PYACTION change
    {
        const c37p::PRun h = c37p::run(tPy, dir.string(), 20);
        const c37p::PRun r = c37p::run(tPy, dir.string(), -1, &h.blob);
        std::printf("    restart from report step 20: FOPT %.2f vs continuous %.2f\n", r.fopt, RP.fopt);
        check(!h.blob.empty() && std::fabs(r.fopt - RP.fopt) < 1e-9 * RP.fopt, "restart: schedule changes made by PYACTION are replayed (state v12)");
    }
    // ---- export keeps PYACTION verbatim, import report says it is evaluated
    {
        std::vector<ressim::DeckMessage> rep;
        const ressim::Case c = ressim::eclipseDeckFromString(tPy, dir.string(), nullptr, &rep);
        const std::string out = ressim::eclipseDeckToString(c);
        bool full = false;
        for (const auto& m : rep) if (m.text.find("PYACTION (Python actions): run by the ResSim engine") != std::string::npos) full = true;
        check(out.find("PYACTION") != std::string::npos && out.find("'cut.py'") != std::string::npos && full && !c.deckDir.empty(),
              "PYACTION kept verbatim on export, import report marks it as evaluated, deck directory recorded");
    }
    std::error_code ec;
    fs::remove_all(dir, ec);
}

// ============================================================================ C37-5: opm-tests triage (GRIDUNIT, WCYCLE)
static void testC37Opm() {
    std::printf("== C37 opm-tests items: GRIDUNIT, WCYCLE\n");
    // ---- GRIDUNIT: COORD / ZCORN given in CM / FEET in a METRIC deck give the same grid and the same run
    {
        std::ifstream f(findExample("FAULTED.DATA"));
        std::stringstream ss; ss << f.rdbuf();
        const std::string base = ss.str();
        std::vector<ressim::DeckMessage> rep;
        const ressim::Case c0 = ressim::eclipseDeckFromString(base, "", nullptr, &rep);
        auto variant = [&](const char* unit, double fac) {
            ressim::Case c = c0;
            for (double& x : c.grid.coord) x *= fac;
            for (double& x : c.grid.zcorn) x *= fac;
            std::string t = ressim::eclipseDeckToString(c);
            const size_t g = t.find("\nGRID");
            if (g == std::string::npos) return std::string();
            const size_t e = t.find('\n', g + 1);
            t.insert(e + 1, std::string("GRIDUNIT\n '") + unit + "' /\n\n");
            return t;
        };
        double worst = 0; bool ok = !c0.grid.coord.empty();
        for (int v = 0; v < 2 && ok; ++v) {
            const std::string t = v == 0 ? variant("CM", 100.0) : variant("FEET", 1.0 / 0.3048);
            ok = !t.empty();
            if (!ok) break;
            std::vector<ressim::DeckMessage> r2;
            const ressim::Case c = ressim::eclipseDeckFromString(t, "", nullptr, &r2);
            ok = c.grid.coord.size() == c0.grid.coord.size() && c.grid.zcorn.size() == c0.grid.zcorn.size();
            for (size_t i = 0; ok && i < c.grid.coord.size(); ++i) worst = std::max(worst, std::fabs(c.grid.coord[i] - c0.grid.coord[i]) / (1 + std::fabs(c0.grid.coord[i])));
            for (size_t i = 0; ok && i < c.grid.zcorn.size(); ++i) worst = std::max(worst, std::fabs(c.grid.zcorn[i] - c0.grid.zcorn[i]) / (1 + std::fabs(c0.grid.zcorn[i])));
        }
        std::printf("    GRIDUNIT CM / FEET -> metres: worst relative geometry difference %.3g\n", worst);
        check(ok && worst < 1e-6, "GRIDUNIT CM and FEET: COORD / ZCORN converted to the deck unit system");
        // GRIDUNIT METRES in a METRIC deck is a no-op; the keyword is reported as modelled
        const std::string t = variant("METRES", 1.0);
        std::vector<ressim::DeckMessage> r3;
        const ressim::Case c = ressim::eclipseDeckFromString(t, "", nullptr, &r3);
        bool same = c.grid.zcorn.size() == c0.grid.zcorn.size();
        for (size_t i = 0; same && i < c.grid.zcorn.size(); ++i) same = std::fabs(c.grid.zcorn[i] - c0.grid.zcorn[i]) < 1e-6 * (1 + std::fabs(c0.grid.zcorn[i]));
        bool ign = false;
        for (const auto& m : r3) if (m.text.find("GRIDUNIT") != std::string::npos && m.text.find("ignored") != std::string::npos) ign = true;
        check(same && !ign, "GRIDUNIT METRES: unchanged geometry, keyword not reported as ignored");
    }
    // ---- WCYCLE on the SPE1 producer: on 20 d, off 40 d, 3 d start-up ramp, 1 d first step, steps land on the switches
    {
        struct WRun { std::vector<ressim::SummaryPoint> sum; std::string log; std::string blob; double fopt = 0; };
        auto run = [](const std::string& text, int stopAfter, const std::string* restart) {
            WRun r;
            std::vector<ressim::DeckMessage> rep;
            const ressim::Case c = ressim::eclipseDeckFromString(text, "", nullptr, &rep);
            ressim::Simulator s; s.log = [&r](const std::string& m) { r.log += m + "\n"; }; s.setup(c); s.initialize();
            if (restart) s.loadState(*restart);
            while ((stopAfter < 0 || s.reportStep() < stopAfter) && s.runReportStep()) {}
            if (stopAfter >= 0) r.blob = s.saveState();
            r.sum = s.summary();
            if (!r.sum.empty()) r.fopt = r.sum.back().fopt;
            return r;
        };
        const std::string yes = c37p::spe1With("WCYCLE\n 'PROD' 20 40 3 1 'YES' /\n/\n");
        const std::string no = c37p::spe1With("WCYCLE\n 'PROD' 20 40 3 1 'NO' /\n/\n");
        const WRun base = run(c37p::spe1With(""), -1, nullptr);
        const WRun a = run(yes, -1, nullptr);
        auto at = [](const WRun& r, double t) { for (const auto& p : r.sum) if (std::fabs(p.time - t) < 1e-6) return &p; return (const ressim::SummaryPoint*)nullptr; };
        bool land = at(a, 20) && at(a, 60) && at(a, 120) && at(a, 180);
        bool offZero = true, onFlow = true; int nOff = 0;
        for (const auto& p : a.sum) {
            const double ph = std::fmod(p.time, 60.0);
            if (p.time > 1 && ph > 1e-6 && ph < 20 - 1e-6 && p.time < 600) onFlow = onFlow && p.fopr > 0;
            if (ph > 20 + 1e-6 && p.time < 600) { offZero = offZero && p.fopr == 0; ++nOff; }
        }
        std::printf("    WCYCLE YES: %zu steps (base %zu), switch points 20/60/120/180 d %s, %d off-period points, FOPT %.6g (base %.6g)\n",
                    a.sum.size(), base.sum.size(), land ? "present" : "MISSING", nOff, a.fopt, base.fopt);
        check(land, "WCYCLE CONTROLLED_TIMESTEP YES: time steps land on the switch times");
        check(offZero && nOff > 0 && onFlow, "WCYCLE: no production in the off periods, production in the on periods");
        // first step after a cycle re-opening <= MAX_TIMESTEP, and its rate is ramped by (dt / START_TIME)
        const ressim::SummaryPoint* p61 = at(a, 61);
        const ressim::SummaryPoint* p64 = nullptr;
        for (const auto& p : a.sum) if (p.time >= 63 - 1e-6 && p.time < 80) { p64 = &p; break; }
        const double ratio = p61 && p64 && p64->fopr > 0 ? p61->fopr / p64->fopr : -1;
        std::printf("    re-opening at 60 d: first step ends at 61 d %s, FOPR(61) / FOPR(>=63) = %.3f (ramp 1/3 expected)\n", p61 ? "yes" : "no", ratio);
        check(p61 && ratio > 0.2 && ratio < 0.5, "WCYCLE: first step after a re-opening <= MAX_TIMESTEP with the start-up efficiency ramp");
        check(a.fopt < 0.6 * base.fopt && a.fopt > 0.2 * base.fopt, "WCYCLE: cumulative production about the on-fraction (1/3) of the base run");
        const WRun b = run(no, -1, nullptr);
        bool noLand = !at(b, 20);
        std::printf("    WCYCLE NO: 20 d point %s, FOPT %.6g\n", noLand ? "absent (switch at a natural step)" : "present", b.fopt);
        check(noLand && b.fopt > 0, "WCYCLE CONTROLLED_TIMESTEP NO: switches at the natural time steps");
        // restart in the middle of an off period is exact
        const WRun h = run(yes, 3, nullptr);
        const WRun r = run(yes, -1, &h.blob);
        std::printf("    restart after report step 3: FOPT %.10g vs continuous %.10g\n", r.fopt, a.fopt);
        check(!h.blob.empty() && std::fabs(r.fopt - a.fopt) <= 1e-9 * std::max(1.0, a.fopt), "WCYCLE: restart reproduces the continuous run");
        // round trip: WCYCLE is written back
        std::vector<ressim::DeckMessage> rep;
        const ressim::Case c = ressim::eclipseDeckFromString(yes, "", nullptr, &rep);
        const std::string out = ressim::eclipseDeckToString(c);
        check(out.find("WCYCLE") != std::string::npos && out.find("'PROD'") != std::string::npos, "WCYCLE: kept in the exported deck");
    }
}
