// C32 tests: multisegment wells (msw.h) - parser, segment solve vs hand formulas, device balance, deck runs
// (MSW_HORIZONTAL.DATA: standard-well limit, friction heel-toe effect, SICD inflow equalisation, SHUT valves, restart)

static void testC32MswUnit() {
    std::printf("== C32 multisegment wells: parser and segment solve\n");
    using namespace ressim::msw;
    std::vector<ExtraKeyword> kw(3);
    kw[0].name = "WELSEGS"; kw[0].step = 0;
    kw[0].body = " 'W' 1000 0 1* INC 'HFA' /\n 2 5 1 1 100 50 0.2 1E-4 /\n 6 7 2 3 2*10 0.15 1E-4 /\n/\n";
    kw[1].name = "COMPSEGS"; kw[1].step = 0;
    kw[1].body = " 'W' /\n 1 1 1 1 160 180 /\n 1 1 2 1 390 410 /\n 2 1 1 2 215 225 /\n/\n";
    kw[2].name = "WSEGVALV"; kw[2].step = 0;
    kw[2].body = " 'W' 7 0.7 0.001 /\n/\n";
    std::string err;
    auto W = parse(kw, 0, true, err);
    bool ok = W.size() == 1 && W[0].segs.size() == 7 && err.empty();
    if (ok) {
        const Well& w = W[0];
        const Segment& s5 = w.segs[size_t(w.indexOf(5))];
        const Segment& s7 = w.segs[size_t(w.indexOf(7))];
        ok = std::fabs(s5.length - 400) < 1e-9 && std::fabs(s5.depth - 1200) < 1e-9 && s7.outlet == 6 && w.segs[size_t(w.indexOf(6))].outlet == 3 &&
             std::fabs(s7.length - 220) < 1e-9 && std::fabs(s7.depth - 1120) < 1e-9 && s7.device == 1 && w.conns.size() == 3 &&
             w.conns[0].seg == 3 && w.conns[1].seg == 5 && w.conns[2].seg == 7;
        std::printf("    seg5 L %.1f z %.1f, seg7 L %.1f z %.1f outlet %d, conn segs %d %d %d %s\n", s5.length, s5.depth, s7.length, s7.depth, s7.outlet,
                    w.conns[0].seg, w.conns[1].seg, w.conns[2].seg, err.c_str());
    }
    check(ok, "WELSEGS INC ranges / branches, COMPSEGS allocation to the nearest node of the branch, WSEGVALV on a segment");
    {   // ABS geometry
        std::vector<ExtraKeyword> k2(1); k2[0].name = "WELSEGS";
        k2[0].body = " 'A' 2000 100 1* ABS /\n 2 5 1 1 500 2400 0.1 1E-5 /\n/\n";
        std::string e2; auto A = parse(k2, 0, false, e2);
        const bool ok2 = A.size() == 1 && A[0].segs.size() == 5 && std::fabs(A[0].segs[2].length - 300) < 1e-9 && std::fabs(A[0].segs[2].depth - 2200) < 1e-9 &&
                         std::fabs(A[0].segs[4].length - 500) < 1e-9;
        check(ok2, "WELSEGS ABS: intermediate nodes interpolated between the outlet and the range end");
    }
    // segment solve: one horizontal water segment, METRIC (bar, m, kg/m3, rm3/d)
    {
        Well w; w.name = "H";
        Segment t; t.id = 1; t.depth = 1000; w.segs.push_back(t);
        Segment s; s.id = 2; s.outlet = 1; s.length = 1000; s.depth = 1000; s.diam = 0.1; s.rough = 1e-4; w.segs.push_back(s);
        Segment v; v.id = 3; v.outlet = 2; v.length = 1000; v.depth = 1100; v.diam = 0.1; v.rough = 1e-4; w.segs.push_back(v);
        Inflow in; in.resize(3);
        in.q[2] = 1000; in.rho[2] = 1000; in.mu[2] = 1.0;   // 1000 rm3/d of water enters at the bottom node
        ressim::msw::Units u; u.g = 9.80665e-5; u.toM = 1; u.toRho = 1; u.toQ = 1.0 / 86400.0; u.toPa = 1e5;
        const Solution S = solve(w, 100.0, in, u);
        const double qs = 1000.0 / 86400.0, A = 0.25 * 3.14159265358979323846 * 0.01, vel = qs / A, Re = 1000 * vel * 0.1 / 1e-3;
        const double t1 = -1.8 * std::log10(std::pow(1e-3 / 3.7, 1.11) + 6.9 / Re), f = 1 / (t1 * t1);
        const double dpF = f * 1000 * vel * vel * 1000 / (2 * 0.1) / 1e5;
        const double dpH = 1000 * 9.80665e-5 * 100;
        std::printf("    friction %.6f bar (hand %.6f), hydrostatic %.6f bar (hand %.6f), p bottom %.6f\n", S.dFric[1], dpF, S.dHydro[2], dpH, S.p[2]);
        check(std::fabs(S.dFric[1] - dpF) < 1e-9 * dpF && std::fabs(S.dHydro[2] - dpH) < 1e-12 && S.dFric[2] == 0 &&
              std::fabs(S.p[2] - (100 + dpF + dpH)) < 1e-9, "segment solve: Haaland friction over the segment length + rho g dz (vertical segment, zero length change)");
        // laminar
        Inflow in2 = in; in2.mu[2] = 500;
        const Solution S2 = solve(w, 100.0, in2, u);
        const double Re2 = 1000 * vel * 0.1 / 0.5, dp2 = 64 / Re2 * 1000 * vel * vel * 1000 / (2 * 0.1) / 1e5;
        check(std::fabs(S2.dFric[1] - dp2) < 1e-9 * dp2, "laminar segment: f = 64 / Re (Hagen-Poiseuille)");
        // valve device balance: dp = c q^2 with q = A (D0 - dp)
        Well wv = w; wv.segs[2].device = 1; wv.segs[2].cv = 0.7; wv.segs[2].area = 2e-4;
        Inflow iv = in; iv.a[2] = 50; iv.pres[2] = 150;   // 50 rm3/d/bar
        const Solution SV = solve(wv, 100.0, iv, u);
        const double c = 1000.0 / (2 * 0.49 * 4e-8) / (86400.0 * 86400.0) / 1e5;
        const double pin = 100 + SV.dFric[1] + SV.dHydro[2] + SV.dFric[2];
        const double qd = 50 * (150 - pin - SV.dDev[2]);
        std::printf("    valve: dp %.5f bar, c q^2 %.5f (q %.2f rm3/d through the device)\n", SV.dDev[2], c * qd * qd, qd);
        check(SV.dDev[2] > 0 && std::fabs(SV.dDev[2] - c * qd * qd) < 1e-6 * SV.dDev[2], "WSEGVALV: device drop balanced with the subtree inflow (dp = rho v^2 / 2 Cv^2 at q = A (D0 - dp))");
    }
}

namespace {
ressim::Case mswVariant(const ressim::Case& base, int variant) {   // 0 standard well, 1 tubing only, 2 tubing + SICD, 3 wide tubing no devices, 4 SICD + heel valves shut
    ressim::Case c = base;
    std::vector<ressim::ExtraKeyword> keep;
    for (auto& e : c.extraKeywords) {
        const bool mswKw = e.name == "WELSEGS" || e.name == "COMPSEGS" || e.name == "WSEGSICD" || e.name == "WSEGVALV";
        if (variant == 0 && mswKw) continue;
        if ((variant == 1 || variant == 3) && e.name == "WSEGSICD") continue;
        if (variant == 3 && e.name == "WELSEGS") { size_t p; while ((p = e.body.find("0.2034")) != std::string::npos) e.body.replace(p, 6, "3.0000"); }
        keep.push_back(e);
    }
    if (variant == 4) {
        ressim::ExtraKeyword v; v.name = "WSEGVALV"; v.section = ressim::DeckSection::Schedule; v.step = 0;
        v.body = " 'P1' 20 0.7 0.001 4* SHUT /\n 'P1' 21 0.7 0.001 4* SHUT /\n/\n";
        keep.push_back(v);
    }
    c.extraKeywords = keep;
    return c;
}
struct MswRun { double fopt = 0, fwpt = 0, wct2 = 0, bt = -1; std::vector<double> devQ; std::vector<ressim::Simulator::MswRecord> rec; int nseg = 0; double maxFric = 0; std::string blob; };
MswRun runMsw(const ressim::Case& c, int stopAfter = -1) {
    MswRun r;
    ressim::Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
    while ((stopAfter < 0 || s.reportStep() < stopAfter) && s.runReportStep()) {}
    if (stopAfter >= 0) { r.blob = s.saveState(); return r; }
    const auto& S = s.summary();
    r.fopt = S.back().fopt; r.fwpt = S.back().fwpt;
    for (const auto& p : S) { if (p.time <= 730.5) r.wct2 = p.fwct; if (r.bt < 0 && p.fwct > 0.1) r.bt = p.time; }
    r.rec = s.mswRecords();
    if (!s.mswWells().empty()) r.nseg = int(s.mswWells()[0].segs.size());
    // device-segment flows (segments 20..37) at the first record after 40 days
    for (const auto& rc : r.rec)
        if (rc.t > 40) {
            for (size_t i = 19; i < rc.q.size() && i < 37; ++i) r.devQ.push_back(rc.q[i]);
            break;
        }
    for (const auto& rc : r.rec) for (double d : rc.dFric) r.maxFric = std::max(r.maxFric, std::fabs(d));
    return r;
}
double cv(const std::vector<double>& v) {
    if (v.empty()) return 0;
    double m = 0, s = 0;
    for (double x : v) m += x;
    m /= double(v.size());
    for (double x : v) s += (x - m) * (x - m);
    return std::sqrt(s / double(v.size())) / std::max(1e-30, std::fabs(m));
}
}  // namespace

static void testC32MswDeck() {
    std::printf("== C32 multisegment wells in the black-oil engine (MSW_HORIZONTAL.DATA)\n");
    const std::string path = findExample("MSW_HORIZONTAL.DATA");
    if (path.empty()) { std::printf("    (examples/MSW_HORIZONTAL.DATA not found - skipped)\n"); return; }
    std::vector<DeckMessage> rep;
    const Case base = importEclipseDeck(path, nullptr, &rep);
    int warnMsw = 0, infoMsw = 0;
    for (const auto& m : rep) {
        const bool k = m.keyword == "WELSEGS" || m.keyword == "COMPSEGS" || m.keyword == "WSEGSICD" || m.text.find("multisegment") != std::string::npos;
        if (k && m.level >= 1) ++warnMsw;
        if (k && m.level == 0) ++infoMsw;
    }
    check(warnMsw == 0 && infoMsw >= 1, "import report: WELSEGS / COMPSEGS / WSEGSICD are evaluated (information, not warnings)");
    {
        const std::string txt = eclipseDeckToString(base);
        check(txt.find("WELSEGS") != std::string::npos && txt.find("WSEGSICD") != std::string::npos && txt.find(" 2 19 1 1 150 0 0.2034") != std::string::npos,
              "export keeps the WELSEGS / COMPSEGS / WSEGSICD text verbatim");
    }
    const MswRun R0 = runMsw(mswVariant(base, 0)), R1 = runMsw(mswVariant(base, 1)), R2 = runMsw(base), R3 = runMsw(mswVariant(base, 3)), R4 = runMsw(mswVariant(base, 4));
    std::printf("    standard well : FOPT %.6g  FWPT %.6g  WCT(2y) %.3f\n", R0.fopt, R0.fwpt, R0.wct2);
    std::printf("    MSW 3 ft tube : FOPT %.6g  FWPT %.6g  WCT(2y) %.3f  max friction %.3g psi\n", R3.fopt, R3.fwpt, R3.wct2, R3.maxFric);
    std::printf("    MSW tubing    : FOPT %.6g  FWPT %.6g  WCT(2y) %.3f  max friction %.3g psi, inflow CV %.3f\n", R1.fopt, R1.fwpt, R1.wct2, R1.maxFric, cv(R1.devQ));
    std::printf("    MSW + SICD    : FOPT %.6g  FWPT %.6g  WCT(2y) %.3f  inflow CV %.3f\n", R2.fopt, R2.fwpt, R2.wct2, cv(R2.devQ));
    std::printf("    SICD + 2 shut : FOPT %.6g  device flows of the shut heel segments %.3g %.3g\n", R4.fopt, R4.devQ.size() > 1 ? R4.devQ[0] : -1.0, R4.devQ.size() > 1 ? R4.devQ[1] : -1.0);
    check(R2.nseg == 37 && !R2.rec.empty() && R2.devQ.size() == 18, "37 segments set up (18 tubing + 18 ICD device segments), segment records per step");
    check(std::fabs(R3.fopt - R0.fopt) < 5e-3 * R0.fopt && R3.maxFric < 1.0,
          "wide frictionless tubing (3 ft, no devices) reproduces the standard well (FOPT within 0.5 %)");
    // friction: heel / toe inflow ratio larger than with the wide tubing
    auto ht = [](const std::vector<double>& q) { return q.size() >= 18 ? (q[0] + q[1] + q[2]) / std::max(1e-12, q[15] + q[16] + q[17]) : 0.0; };
    const MswRun& RW = R3;
    std::vector<double> qWide;
    for (const auto& rc : RW.rec) if (rc.t > 40) { for (size_t i = 19; i < rc.q.size() && i < 37; ++i) qWide.push_back(rc.q[i]); break; }
    std::printf("    device flows (heel..toe) SICD:"); for (double q : R2.devQ) std::printf(" %.0f", q); std::printf("\n    tubing only:"); for (double q : R1.devQ) std::printf(" %.0f", q); std::printf("\n");
    std::printf("    heel / toe inflow: wide tubing %.3f, 2.441 in tubing %.3f, with SICD %.3f\n", ht(qWide), ht(R1.devQ), ht(R2.devQ));
    check(R1.maxFric > 5 && ht(R1.devQ) > ht(qWide) * 1.02, "tubing friction (2.441 in) raises the heel share of the inflow (heel / toe ratio up)");
    check(cv(R2.devQ) < 0.7 * cv(R1.devQ), "spiral ICDs equalise the inflow along the well (coefficient of variation down >= 30 %)");
    std::printf("    water breakthrough (WCT > 0.1): tubing %.0f d, SICD %.0f d\n", R1.bt, R2.bt);
    check(R2.bt > R1.bt && R2.fwpt < R1.fwpt && R2.fopt > R1.fopt, "ICDs delay the water breakthrough, produce less water and more oil than the open tubing");
    check(R4.devQ.size() >= 2 && std::fabs(R4.devQ[0]) < 1e-9 && std::fabs(R4.devQ[1]) < 1e-9 && R4.fopt > 0, "WSEGVALV SHUT closes the connections of its subtree (zero flow through the two heel devices)");
    // restart (relaxation state per segment in the v5 wellbore fields)
    {
        const MswRun a = runMsw(base, 7);
        Simulator b; b.log = [](const std::string&) {}; b.setup(base); b.initialize(); b.loadState(a.blob);
        while (b.runReportStep()) {}
        const double f1 = b.summary().back().fopt;
        std::printf("    restart after report step 7: FOPT %.10g vs continuous %.10g\n", f1, R2.fopt);
        check(std::fabs(f1 - R2.fopt) <= 1e-9 * R2.fopt, "restart of a multisegment well continues identically");
    }
}

// ================================================================ C32 geomodelling workflow (geomodel.h)
#include "ressim/geomodel.h"

static void testC32GeoFormats() {
    std::printf("== C32 geomodel: exchange formats (IRAP / ZMAP / XYZ / LAS / SEG-Y), minimum curvature, gridding\n");
    using namespace ressim::geo;
    namespace fs = std::filesystem;
    const std::string dir = (fs::temp_directory_path() / "ressim_c32_fmt").string();
    fs::create_directories(dir);
    {   // IRAP round trip (rotated lattice)
        Surface s; s.name = "T"; s.nx = 7; s.ny = 5; s.x0 = 1000; s.y0 = 2000; s.dx = 25; s.dy = 40; s.rotDeg = 30;
        for (int j = 0; j < s.ny; ++j) for (int i = 0; i < s.nx; ++i) s.z.push_back(1500 + 3 * i - 2 * j + 0.125 * i * j);
        s.z[3] = std::numeric_limits<double>::quiet_NaN();
        const bool w = writeIrap(dir + "/t.irap", s);
        std::string e; const Surface r = readIrap(dir + "/t.irap", e);
        double md = 0; bool nanOk = r.valid() && std::isnan(r.z[3]);
        if (r.valid()) for (size_t k = 0; k < s.z.size(); ++k) if (k != 3) md = std::max(md, std::fabs(r.z[k] - s.z[k]));
        double x, y; s.nodeXY(4, 2, x, y);
        const double a = r.valid() ? r.at(x, y) : -1;
        std::printf("    IRAP: %dx%d rot %.1f, max node diff %.2e, at(node 4,2) %.4f (expected %.4f) %s\n", r.nx, r.ny, r.rotDeg, md, a, s.node(4, 2), e.c_str());
        check(w && r.nx == 7 && r.ny == 5 && std::fabs(r.rotDeg - 30) < 1e-9 && md < 1e-3 && nanOk && std::fabs(a - s.node(4, 2)) < 1e-3,
              "IRAP classic: write / read round trip with rotation and undefined nodes, bilinear sampling");
    }
    {   // ZMAP+ (hand-written, 3 x 4 nodes, column-major from the top-left)
        FILE* f = std::fopen((dir + "/t.zmap").c_str(), "w");
        std::fprintf(f, "! test\n@T, GRID, 4\n 15, 1E+30, , 4, 1\n 4, 3, 0.0, 200.0, 0.0, 300.0\n 0.0, 0.0, 0.0\n@\n");
        // columns x = 0, 100, 200; each column from ymax (300) to ymin (0): z = 1000 + x + 2 y
        for (int i = 0; i < 3; ++i) for (int j = 3; j >= 0; --j) std::fprintf(f, " %.1f", 1000.0 + 100 * i + 2 * 100.0 * j);
        std::fprintf(f, "\n");
        std::fclose(f);
        std::string e; const Surface z = readZmap(dir + "/t.zmap", e);
        const double v = z.valid() ? z.at(150, 250) : -1;
        std::printf("    ZMAP: %dx%d, z(150,250) = %.3f (expected %.3f) %s\n", z.nx, z.ny, v, 1000 + 150 + 2 * 250.0, e.c_str());
        check(z.nx == 3 && z.ny == 4 && std::fabs(v - 1650) < 1e-6, "ZMAP+ grid: header, column-major node order from the north-west corner");
    }
    {   // XYZ points of a plane -> kriging on a lattice reproduces the plane
        FILE* f = std::fopen((dir + "/p.xyz").c_str(), "w");
        std::mt19937 r(5); std::uniform_real_distribution<double> U(0, 1000);
        for (int k = 0; k < 60; ++k) { const double x = U(r), y = U(r); std::fprintf(f, "%.3f %.3f %.4f\n", x, y, 2000 + 0.05 * x - 0.02 * y); }
        std::fclose(f);
        std::string e; const auto pts = readXyz(dir + "/p.xyz", e);
        Surface t; t.nx = 11; t.ny = 11; t.dx = t.dy = 100; t.z.assign(121, 0);
        ressim::geostat::Variogram vg; vg.model = ressim::geostat::VarioModel::Gaussian; vg.rangeMax = vg.rangeMin = 2000; vg.nugget = 0; vg.sill = 1;
        const Surface g = gridPoints(pts, t, vg, 16);
        double err = 0; for (int j = 2; j < 9; ++j) for (int i = 2; i < 9; ++i) err = std::max(err, std::fabs(g.node(i, j) - (2000 + 5.0 * i - 2.0 * j)));
        std::printf("    XYZ: %zu points, kriged plane max error inside %.3f\n", pts.size(), err);
        check(pts.size() == 60 && err < 0.5, "scattered XYZ points: ordinary kriging onto a lattice reproduces a planar surface");
    }
    {   // LAS 2.0 round trip
        WellLog L; L.well = "T1";
        LogCurve a{"GR", "GAPI", "gamma", {}, false}, b{"PHIE", "V/V", "phi", {}, false};
        for (int k = 0; k < 50; ++k) { L.md.push_back(1000 + 0.5 * k); a.v.push_back(50 + k); b.v.push_back(k == 7 ? std::numeric_limits<double>::quiet_NaN() : 0.2 + 0.001 * k); }
        L.curves = {a, b};
        const bool w = writeLas(dir + "/t.las", L);
        std::string e; const WellLog r = readLas(dir + "/t.las", e);
        const LogCurve* c = r.curve("phie");
        const double v = c ? r.valueAt(*c, 1010.25) : -1;
        std::printf("    LAS: well '%s', %zu samples, %zu curves, PHIE(1010.25) %.5f, null -> NaN %d %s\n", r.well.c_str(), r.md.size(), r.curves.size(), v, c ? int(std::isnan(c->v[7])) : -1, e.c_str());
        check(w && r.md.size() == 50 && r.curves.size() == 2 && c && std::isnan(c->v[7]) && std::fabs(v - 0.2205) < 1e-4,
              "LAS 2.0: write / read round trip, NULL value -> undefined, curve lookup and depth interpolation");
    }
    {   // SEG-Y round trip (IEEE) and IBM float decoding
        SeismicCube c; c.ni = 4; c.nj = 3; c.ns = 20; c.dt = 2; c.t0 = 100; c.il0 = 10; c.xl0 = 50; c.dil = 1; c.dxl = 1;
        c.x0 = 500000; c.y0 = 6000000; c.ix = 12.5; c.iy = 0; c.jx = 0; c.jy = 25;
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 3; ++j) for (int s = 0; s < 20; ++s) c.a.push_back(float(std::sin(0.3 * s + i) * (j + 1)));
        const bool w = writeSegy(dir + "/t.sgy", c);
        std::string e; const SeismicCube r = readSegy(dir + "/t.sgy", e);
        double md = 0; if (r.a.size() == c.a.size()) for (size_t k = 0; k < c.a.size(); ++k) md = std::max(md, double(std::fabs(r.a[k] - c.a[k])));
        double x, y; r.xy(3, 2, x, y);
        std::printf("    SEG-Y: %dx%dx%d dt %.1f t0 %.0f il0 %d xl0 %d, max sample diff %.1e, trace (3,2) at %.2f %.2f %s\n", r.ni, r.nj, r.ns, r.dt, r.t0, r.il0, r.xl0, md, x, y, e.c_str());
        check(w && r.ni == 4 && r.nj == 3 && r.ns == 20 && std::fabs(r.dt - 2) < 1e-9 && std::fabs(r.t0 - 100) < 1e-9 && r.il0 == 10 && r.xl0 == 50 && md == 0 &&
              std::fabs(x - 500037.5) < 0.01 && std::fabs(y - 6000050) < 0.01, "SEG-Y rev1: IEEE samples, inline / crossline, CDP coordinates with the scalar, sample interval, delay");
        // IBM: patch the format code to 1 and the first sample to 0xC276A000 (-118.625)
        std::fstream f(dir + "/t.sgy", std::ios::in | std::ios::out | std::ios::binary);
        const unsigned char fmt[2] = {0, 1}, ibm[4] = {0xC2, 0x76, 0xA0, 0x00};
        f.seekp(3200 + 24); f.write(reinterpret_cast<const char*>(fmt), 2);
        f.seekp(3600 + 240); f.write(reinterpret_cast<const char*>(ibm), 4);
        f.close();
        const SeismicCube b = readSegy(dir + "/t.sgy", e);
        std::printf("    IBM float 0xC276A000 -> %.4f\n", b.a.empty() ? 0.0 : double(b.a[0]));
        check(b.format == 1 && !b.a.empty() && std::fabs(b.a[0] + 118.625) < 1e-4, "SEG-Y IBM System/360 float decoding (format 1)");
    }
    {   // minimum curvature: constant build 0 -> 90 deg over a radius R gives TVD R and departure R
        GeoWell w; w.name = "MC"; w.x = 0; w.y = 0; w.kb = 0;
        const double R = 300, pi = 3.14159265358979323846;
        for (int k = 0; k <= 30; ++k) { Station s; s.md = 1000 + k * (R * pi / 2) / 30; s.inc = 3.0 * k; s.azi = 90; w.traj.push_back(s); }
        Station s0; s0.md = 0; s0.inc = 0; s0.azi = 90; w.traj.insert(w.traj.begin(), s0);
        minimumCurvature(w);
        const Station& e = w.traj.back();
        std::printf("    minimum curvature: end TVD %.4f (expected %.4f), east %.4f (expected %.4f), north %.2e\n", e.tvd, 1000 + R, e.x, R, e.y);
        check(std::fabs(e.tvd - (1000 + R)) < 1e-6 && std::fabs(e.x - R) < 1e-6 && std::fabs(e.y) < 1e-6, "minimum curvature: exact on a circular build section (TVD and departure = radius)");
        const P3 p = w.at(1000 + R * pi / 4);
        check(std::fabs(p[2] - (1000 + R * std::sin(pi / 4))) < 0.5 && std::fabs(w.mdAtTvd(1000 + R * std::sin(pi / 4)) - (1000 + R * pi / 4)) < 1.0,
              "trajectory position / MD at TVD inside a build section");
    }
    std::error_code ec; fs::remove_all(dir, ec);
}

namespace {
struct GeoProj {
    std::vector<ressim::geo::Surface> hz;
    std::vector<ressim::geo::Fault> faults;
    std::vector<ressim::geo::GeoWell> wells;
    ressim::geo::SeismicCube cube;
    std::string err;
};
GeoProj loadDemo(const std::string& dir, bool field) {
    using namespace ressim::geo;
    GeoProj P;
    std::string e;
    P.err.clear();
    if (writeDemoProject(dir, field, e).empty() || !e.empty()) { P.err = e; return P; }
    for (const char* h : {"TopRes", "MidRes", "BaseRes"}) { Surface s = readIrap(dir + "/" + h + ".irap", e); s.name = h; P.hz.push_back(s); }
    P.faults = readFaultSticks(dir + "/F1_sticks.txt", e);
    P.wells = readWellHeads(dir + "/well_heads.txt", e);
    for (auto& w : P.wells) {
        if (std::filesystem::exists(dir + "/" + w.name + "_survey.txt")) readDeviation(dir + "/" + w.name + "_survey.txt", w, e);
        w.logs.push_back(readLas(dir + "/" + w.name + ".las", e));
    }
    for (const auto& t : readTops(dir + "/well_tops.txt", e)) for (auto& w : P.wells) if (w.name == t.first) w.tops.push_back(t.second);
    P.cube = readSegy(dir + "/seismic.segy", e);
    if (!e.empty()) P.err = e;
    return P;
}
}  // namespace

static void testC32GeoWorkflow() {
    std::printf("== C32 geomodel: demo project workflow (seismic -> structure -> grid -> facies -> petrophysics -> volumes -> upscaling -> simulation)\n");
    using namespace ressim::geo;
    namespace fs = std::filesystem;
    const std::string dir = (fs::temp_directory_path() / "ressim_c32_demo").string();
    GeoProj P = loadDemo(dir, false);
    size_t ntops = 0; for (auto& w : P.wells) ntops += w.tops.size();
    double w6reach = 0;
    for (auto& w : P.wells) if (w.name == "W6" && !w.traj.empty()) w6reach = std::hypot(w.traj.back().x - w.x, w.traj.back().y - w.y);
    std::printf("    loaded: %zu horizons, %zu fault(s) with %zu sticks, %zu wells, %zu tops, W6 reach %.0f m, seismic %dx%dx%d %s\n", P.hz.size(), P.faults.size(),
                P.faults.empty() ? size_t(0) : P.faults[0].sticks.size(), P.wells.size(), ntops, w6reach, P.cube.ni, P.cube.nj, P.cube.ns, P.err.c_str());
    check(P.err.empty() && P.hz.size() == 3 && P.hz[0].valid() && P.faults.size() == 1 && P.faults[0].sticks.size() == 5 && P.wells.size() == 6 && ntops == 18 &&
              w6reach > 500 && P.cube.ni == 41 && P.cube.nj == 31 && P.cube.ns == 176,
          "demo project: IRAP horizons, fault sticks, well heads + survey + LAS + tops, SEG-Y volume all read back");
    if (!P.err.empty() || P.hz.size() < 3 || P.wells.size() < 6) return;
    // ---- seismic interpretation: auto-track the top reservoir peak from one seed, depth-convert, compare with the map
    VelocityModel vm; vm.v0 = 1800; vm.k = 0.6;
    const double zt = 2003.0;
    check(std::fabs(timeToDepth(vm, depthToTime(vm, zt)) - zt) < 1e-9 && std::fabs(depthToTime(vm, 0)) < 1e-12, "velocity model v0 + k z: time <-> depth round trip");
    {
        double si, sj; P.cube.ij(1200, 1500, si, sj);
        const double z0 = P.hz[0].at(1200, 1500);
        double si2, sj2; P.cube.ij(3400, 1500, si2, sj2);   // one seed per fault block
        const TimeHorizon th = autoTrack(P.cube, {{si, sj, depthToTime(vm, z0)}, {si2, sj2, depthToTime(vm, P.hz[0].at(3400, 1500))}}, 0, 12, 12);
        int picked = 0; for (double t : th.t) picked += std::isfinite(t);
        const Surface d = horizonToDepth(P.cube, th, vm, "TopRes_seis");
        double s2 = 0; int n = 0, nf = 0;
        for (int j = 0; j < d.ny; ++j) for (int i = 0; i < d.nx; ++i) {
            double x, y; d.nodeXY(i, j, x, y);
            const double a = d.node(i, j), b = P.hz[0].at(x, y);
            if (!std::isfinite(a) || !std::isfinite(b)) continue;
            if (std::fabs(x - 2600 - 0.1 * (y - 1500)) < 150) { ++nf; continue; }   // fault zone excluded
            s2 += (a - b) * (a - b); ++n;
        }
        const double rms = n ? std::sqrt(s2 / n) : 1e9;
        std::printf("    auto-track: %d / %d traces picked, depth-converted vs interpreted map RMS %.2f m over %d nodes\n", picked, P.cube.ni * P.cube.nj, rms, n);
        check(picked > 0.9 * P.cube.ni * P.cube.nj && rms < 6.0, "seeded peak auto-tracking (a seed per fault block) + v0 + kz depth conversion agrees with the depth map (map error 4 m)");
        const std::vector<P3> pts = {{1200, 1500, z0}, {1200, 1500, z0 + 40}};
        const auto at = extractAttribute(P.cube, vm, pts, 0, 0);
        std::printf("    amplitude at the top / 40 m below: %.3f / %.3f\n", at[0], at[1]);
        check(at.size() == 2 && at[0] > 0.6 && std::fabs(at[1]) < at[0], "attribute extraction at depth (amplitude at the top-reservoir peak)");
    }
    // ---- structural grid
    GridSpec gs; gs.x0 = 0; gs.y0 = 0; gs.dx = 100; gs.dy = 100; gs.nx = 40; gs.ny = 30;
    Zone z1; z1.name = "Upper"; z1.top = 0; z1.base = 1; z1.nLayers = 8;
    Zone z2; z2.name = "Lower"; z2.top = 1; z2.base = 2; z2.nLayers = 8;
    gs.zones = {z1, z2};
    gs.tieToWells = false;
    std::string e;
    const GeoGrid g0 = buildGrid(gs, P.hz, P.faults, P.wells, e);
    gs.tieToWells = true; gs.tieRange = 1500;
    GeoGrid g = buildGrid(gs, P.hz, P.faults, P.wells, e);
    const int nx = g.grid.nx, ny = g.grid.ny, nz = g.grid.nz;
    // throw across the fault faces (top of layer 1)
    double thr = 0; int nt = 0;
    for (const auto& F : g.grid.faults) for (const auto& f : F.faces) for (int j = f.j1; j <= f.j2; ++j) for (int i = f.i1; i <= f.i2; ++i) {
        if (f.face != 'X' || i + 1 >= nx) continue;
        auto zc = [&](int ii, int a, int b) { return g.grid.zcorn[size_t((0 * 2 * ny + 2 * j + b) * 2 * nx + 2 * ii + a)]; };
        thr += 0.5 * ((zc(i + 1, 0, 0) - zc(i, 1, 0)) + (zc(i + 1, 0, 1) - zc(i, 1, 1))); ++nt;
    }
    thr = nt ? thr / nt : 0;
    size_t nfaces = 0; for (const auto& F : g.grid.faults) nfaces += F.faces.size();
    double vol = 0; for (double v : g.bulk) vol += v;
    std::printf("    grid %dx%dx%d, %d fault blocks, %zu fault face records, mean throw across X faces %.1f m, bulk %.4g m3, %s\n", nx, ny, nz, g.nBlocks, nfaces, thr, vol, e.c_str());
    check(e.empty() && nz == 16 && g.nBlocks == 2 && nfaces > 0 && thr > 20 && thr < 40, "pillar grid: fault blocks from the stick traces, split ZCORN nodes carry the throw, FAULTS faces");
    // well tie: grid top at the well tops before / after
    double r0 = 0, r1 = 0; int nw = 0;
    for (const auto& w : P.wells) for (const auto& t : w.tops) {
        if (t.surface != "TopRes") continue;
        const P3 p = w.at(t.md);
        double fi, fj; g.columnOf(p[0], p[1], fi, fj);
        const int i = int(fi), j = int(fj);
        if (i < 0 || j < 0 || i >= nx || j >= ny) continue;
        const double u = fi - i, v = fj - j;
        auto top = [&](const GeoGrid& G) {
            auto zc = [&](int a, int b) { return G.grid.zcorn[size_t((0 * 2 * ny + 2 * j + b) * 2 * nx + 2 * i + a)]; };
            return (1 - u) * (1 - v) * zc(0, 0) + u * (1 - v) * zc(1, 0) + (1 - u) * v * zc(0, 1) + u * v * zc(1, 1);
        };
        r0 += std::pow(top(g0) - p[2], 2); r1 += std::pow(top(g) - p[2], 2); ++nw;
    }
    r0 = nw ? std::sqrt(r0 / nw) : 0; r1 = nw ? std::sqrt(r1 / nw) : 0;
    std::printf("    well tie (TopRes, %d wells): RMS misfit %.2f m before, %.2f m after (%zu report lines)\n", nw, r0, r1, g.tieReport.size());
    check(nw >= 5 && r0 > 1.0 && r1 < 0.35 * r0, "well tie: kriged residual correction moves the horizons onto the well tops");
    bool okGeom = true;
    for (size_t c = 0; c < g.bulk.size(); ++c) if (g.bot[c] < g.top[c] - 1e-9) okGeom = false;
    check(okGeom && g.layerZone.size() == 16 && g.layerZone[7] == 0 && g.layerZone[8] == 1, "zones: proportional layering, no inverted cells, layer -> zone map");
    // ---- blocking
    const auto B = blockWells(g, P.wells);
    std::map<std::string, int> per; double phiS = 0; int nS = 0;
    for (const auto& b : B) { per[b.well]++; auto f = b.v.find("FACIES"); auto p = b.v.find("PHIE"); if (f != b.v.end() && p != b.v.end() && f->second == 1) { phiS += p->second; ++nS; } }
    std::printf("    blocked cells: %zu (W1 %d, W6 %d), sand cells %d with mean PHIE %.3f\n", B.size(), per["W1"], per["W6"], nS, nS ? phiS / nS : 0.0);
    check(per.size() == 6 && per["W1"] >= 14 && per["W6"] > per["W1"] && nS > 10 && phiS / std::max(1, nS) > 0.15,
          "well-log blocking along vertical and horizontal trajectories (arithmetic PHIE, most-of FACIES)");
    // ---- facies
    FaciesSpec fs0; fs0.method = 0; fs0.vg.rangeMax = fs0.vg.rangeMin = 800; fs0.vg.rangeVert = 6; fs0.seed = 3;
    const FaciesResult F0 = modelFacies(g, B, fs0);
    FaciesSpec fs1 = fs0; fs1.method = 1;
    const FaciesResult F1 = modelFacies(g, B, fs1);
    FaciesSpec fs2 = fs0; fs2.method = 2; fs2.codes = {0, 1}; fs2.proportions = {0.5, 0.5}; fs2.chWidth = 400; fs2.chThick = 8; fs2.chAzimuth = 30;
    const FaciesResult F2 = modelFacies(g, B, fs2);
    int nData = 0; std::vector<double> pd(2, 0);
    for (const auto& b : B) { auto f = b.v.find("FACIES"); if (f != b.v.end()) { ++nData; pd[size_t(f->second > 0.5)] += 1; } }
    const double pSand = pd[1] / std::max(1, nData);
    std::printf("    facies data %d (sand %.2f): SIS sand %.3f honoured %d, TGS sand %.3f honoured %d, channels sand %.3f (mismatch %d forced)\n", nData, pSand,
                F0.proportions.size() > 1 ? F0.proportions[1] : -1, F0.conditioned, F1.proportions.size() > 1 ? F1.proportions[1] : -1, F1.conditioned,
                F2.proportions.size() > 1 ? F2.proportions[1] : -1, F2.mismatches);
    check(F0.conditioned == nData && std::fabs(F0.proportions[1] - pSand) < 0.1, "SIS facies: honours every blocked cell, reproduces the data proportion (+-0.1)");
    check(F1.conditioned >= nData - 1 && std::fabs(F1.proportions[1] - pSand) < 0.12, "truncated Gaussian facies: thresholds from the proportions, data honoured");
    check(F2.conditioned == nData && std::fabs(F2.proportions[1] - 0.5) < 0.08, "object-based channels: target channel proportion, wells forced to their facies");
    // ---- petrophysics: SGS per facies, then collocated co-simulation with the extracted seismic amplitude
    PetroSpec ps; ps.curve = "PHIE"; ps.vg.rangeMax = ps.vg.rangeMin = 1000; ps.vg.rangeVert = 8; ps.minV = 0; ps.maxV = 0.4; ps.seed = 11;
    const PetroResult phi = modelProperty(g, B, &F0.facies, ps);
    double hon = 0;
    for (const auto& b : B) { auto p = b.v.find("PHIE"); if (p != b.v.end()) hon = std::max(hon, std::fabs(phi.v[size_t(b.cell)] - p->second)); }
    std::printf("    PHIE SGS per facies: mean %.4f sd %.4f, data %d, max misfit at the blocked cells %.2e (%s)\n", phi.mean, phi.sd, phi.data, hon, phi.note.c_str());
    check(phi.data > 50 && hon < 2e-3 && phi.mean > 0.08 && phi.mean < 0.25, "porosity by SGS per facies: blocked values honoured, plausible distribution");
    std::vector<P3> topPts; for (const auto& c : g.center) topPts.push_back({c[0], c[1], c[2]});
    // secondary: amplitude of the top-reservoir peak at each column (porosity-sensitive)
    std::vector<double> sec(g.bulk.size(), 0.0);
    {
        std::vector<P3> colTop;
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) { const size_t c = size_t(i + nx * j); colTop.push_back({g.center[c][0], g.center[c][1], g.top[c]}); }
        const auto am = extractAttribute(P.cube, vm, colTop, 1, 12);
        for (size_t c = 0; c < sec.size(); ++c) sec[c] = am[c % size_t(nx * ny)];
    }
    PetroSpec pc = ps; pc.perFacies = false; pc.secondary = sec; pc.corrSecondary = 0.7;
    PetroSpec pn = ps; pn.perFacies = false;
    const PetroResult phC = modelProperty(g, B, nullptr, pc), phN = modelProperty(g, B, nullptr, pn);
    auto corr = [&](const std::vector<double>& a) {
        double ma = 0, mb = 0, sab = 0, saa = 0, sbb = 0; int n = 0;
        for (size_t c = 0; c < a.size(); ++c) if (g.grid.actnum[c]) { ma += a[c]; mb += sec[c]; ++n; }
        ma /= n; mb /= n;
        for (size_t c = 0; c < a.size(); ++c) if (g.grid.actnum[c]) { sab += (a[c] - ma) * (sec[c] - mb); saa += (a[c] - ma) * (a[c] - ma); sbb += (sec[c] - mb) * (sec[c] - mb); }
        return sab / std::sqrt(std::max(1e-30, saa * sbb));
    };
    const double cC = corr(phC.v), cN = corr(phN.v);
    std::printf("    correlation with the seismic amplitude: SGS %.3f, collocated co-simulation (rho 0.7) %.3f\n", cN, cC);
    check(cC > cN + 0.2 && cC > 0.35, "collocated co-simulation (Markov-1 residual form) imprints the seismic secondary");
    // ---- poro-perm
    const auto fits = fitPoroPerm(B, "PHIE", "PERM", &F0.facies, F0.codes);
    std::printf("    poro-perm per facies:"); for (size_t q = 0; q < fits.size(); ++q) std::printf(" [code %d: log k = %.2f + %.2f phi, r2 %.2f, n %d]", F0.codes[q], fits[q].a, fits[q].b, fits[q].r2, fits[q].n); std::printf("\n");
    bool fitOk = false; for (size_t q = 0; q < fits.size(); ++q) if (F0.codes[q] == 1 && fits[q].b > 8 && fits[q].b < 25 && fits[q].r2 > 0.3) fitOk = true;
    check(fitOk, "poro-perm transform fitted on the blocked sand (log10 k slope 13 in the logs; mixed-cell blocking - arithmetic phi, geometric k - steepens it)");
    const auto kx = permFromPoro(phi.v, fits, &F0.facies, F0.codes, 0.2, 5);
    std::vector<double> kz(kx.size()); for (size_t c = 0; c < kx.size(); ++c) kz[c] = 0.1 * kx[c];
    std::vector<double> ntg(g.bulk.size(), 1.0);
    // ---- volumetrics
    VolSpec vs; vs.contacts.owc = 2080; vs.contacts.goc = -1e30; vs.satHeight = false; vs.swConst = 0.25; vs.bo = 1.25;
    const VolResult V = volumetrics(g, phi.v, ntg, &F0.facies, vs);
    double hand = 0; for (size_t c = 0; c < g.bulk.size(); ++c) if (g.grid.actnum[c] && g.center[c][2] < 2080) hand += g.bulk[c] * phi.v[c] * 0.75 / 1.25;
    double zsum = 0; for (const auto& r : V.byZone) zsum += r.stoiip;
    std::printf("    STOIIP %.4g sm3 (hand %.4g), zones %zu sum %.4g, blocks %zu, facies %zu\n", V.total.stoiip, hand, V.byZone.size(), zsum, V.byBlock.size(), V.byFacies.size());
    check(std::fabs(V.total.stoiip - hand) < 1e-6 * hand && std::fabs(zsum - hand) < 1e-6 * hand && V.byZone.size() == 2 && V.byBlock.size() == 2,
          "volumetrics: STOIIP = sum bulk NTG phi (1 - Sw) / Bo above the OWC, split by zone / fault block / facies");
    VolSpec vh = vs; vh.satHeight = true; vh.sh.swirr = 0.15; vh.sh.pd = 0.05; vh.sh.drhoOW = 0.0225; vh.sh.lambda = 2;   // bar, bar/m
    const VolResult VH = volumetrics(g, phi.v, ntg, nullptr, vh);
    VolSpec v15 = vs; v15.swConst = 0.15;
    const VolResult V15 = volumetrics(g, phi.v, ntg, nullptr, v15);
    std::vector<double> real;
    for (unsigned s = 1; s <= 5; ++s) { PetroSpec q = ps; q.seed = 100 + s; real.push_back(volumetrics(g, modelProperty(g, B, &F0.facies, q).v, ntg, nullptr, vs).total.stoiip); }
    const Percentiles pc5 = percentiles(real);
    std::printf("    saturation-height STOIIP %.4g (constant Sw = Swirr 0.15: %.4g); 5 realisations P90 %.4g P50 %.4g P10 %.4g\n", VH.total.stoiip, V15.total.stoiip, pc5.p90, pc5.p50, pc5.p10);
    check(VH.total.stoiip < V15.total.stoiip && VH.total.stoiip > 0.7 * V15.total.stoiip && pc5.p90 <= pc5.p50 && pc5.p50 <= pc5.p10 && pc5.p10 > pc5.p90,
          "saturation-height (Brooks-Corey above the free-water level) lowers the volume near the contact; realisation percentiles ordered");
    // ---- upscaling
    {
        std::vector<double> kl(2 * 2 * 4); for (int l = 0; l < 4; ++l) for (int q = 0; q < 4; ++q) kl[size_t(q + 4 * l)] = l % 2 ? 10.0 : 1000.0;
        const double kxE = flowBasedPerm(2, 2, 4, kl, 10, 10, 1, 0), kzE = flowBasedPerm(2, 2, 4, kl, 10, 10, 1, 2);
        std::printf("    flow-based box (layers 1000 / 10 mD): kx %.4f (arithmetic 505), kz %.4f (harmonic %.4f)\n", kxE, kzE, 4.0 / (2 / 1000.0 + 2 / 10.0));
        check(std::fabs(kxE - 505) < 1e-6 * 505 && std::fabs(kzE - 4.0 / (2 / 1000.0 + 2 / 10.0)) < 1e-6 * kzE, "flow-based upscaling is exact for layered media (arithmetic along, harmonic across the layers)");
    }
    UpscaleSpec us; us.fx = 2; us.fy = 2; us.fz = 3; us.permMethod = 3;
    const Upscaled Up = upscale(g, phi.v, ntg, kx, kx, kz, us);
    double pvF = 0, pvC = 0;
    for (size_t c = 0; c < g.bulk.size(); ++c) if (g.grid.actnum[c]) pvF += g.bulk[c] * phi.v[c];
    for (size_t c = 0; c < Up.g.bulk.size(); ++c) if (Up.g.grid.actnum[c]) pvC += Up.g.bulk[c] * Up.poro[c] * Up.ntg[c];
    std::printf("    upscaled %dx%dx%d (zone-respecting layers %zu), PV fine %.6g coarse %.6g, faults %zu\n", Up.g.grid.nx, Up.g.grid.ny, Up.g.grid.nz, Up.g.layerZone.size(), pvF, pvC, Up.g.grid.faults.size());
    check(Up.g.grid.nx == 20 && Up.g.grid.ny == 15 && Up.g.grid.nz == 6 && std::fabs(pvC - pvF) < 1e-9 * pvF && !Up.g.grid.faults.empty(),
          "upscaling: zone-respecting coarsening (8 -> 3+3 layers per zone), pore volume conserved, fault faces kept on block boundaries");
    // ---- export + simulation
    const std::string gp = dir + "/model.grdecl";
    const bool wg = writeGrdecl(gp, Up.g, {{"PORO", &Up.poro}, {"PERMX", &Up.kx}, {"PERMY", &Up.ky}, {"PERMZ", &Up.kz}, {"NTG", &Up.ntg}});
    std::string txt; { std::ifstream f(gp); std::stringstream ss; ss << f.rdbuf(); txt = ss.str(); }
    check(wg && txt.find("SPECGRID") != std::string::npos && txt.find("ZCORN") != std::string::npos && txt.find("FAULTS") != std::string::npos && txt.find("MULTFLT") != std::string::npos,
          "GRDECL export: SPECGRID / COORD / ZCORN / ACTNUM / properties / FAULTS / MULTFLT");
    const std::string ex = findExample("CHANNELS.DATA");
    if (!ex.empty()) {
        ressim::Case c = ressim::importEclipseDeck(ex, nullptr, nullptr);
        c.wells.clear();
        toCase(c, Up.g, Up.poro, Up.ntg, Up.kx, Up.ky, Up.kz, P.wells);
        c.equil.datumDepth = 2050; c.equil.datumPressure = 210; c.equil.owc = 2080; c.equil.goc = 1800; c.equil.rsvdDepth.clear(); c.equil.rsvdRs.clear();
        c.schedule.resize(4);
        for (auto& s : c.schedule) { s.length = 60; s.controls.clear(); s.gprod.clear(); s.ginj.clear(); }
        for (const auto& w : c.wells) { ressim::WellControl wc; wc.well = w.name; wc.mode = ressim::Control::LRAT; wc.lrat = 150; wc.bhp = 120; c.schedule[0].controls.push_back(wc); }
        c.extraKeywords.clear();
        int ncomp = 0; for (const auto& w : c.wells) ncomp += int(w.completions.size());
        ressim::Simulator s; s.log = [](const std::string&) {};
        bool ok = true;
        try { s.setup(c); s.initialize(); while (s.runReportStep()) {} } catch (const std::exception& ex2) { std::printf("    simulation error: %s\n", ex2.what()); ok = false; }
        const double fopt = ok && !s.summary().empty() ? s.summary().back().fopt : 0;
        std::printf("    simulation case: %zu wells, %d connections, FOPT after 240 d %.4g sm3 (STOIIP %.4g)\n", c.wells.size(), ncomp, fopt, V.total.stoiip);
        check(ok && c.wells.size() == 6 && ncomp >= 6 * 3 && fopt > 0 && fopt < V.total.stoiip, "geomodel -> simulation case: corner-point grid + properties + wells along the trajectories run in the engine");
    }
    std::error_code ec2; fs::remove_all(dir, ec2);
}

// ---------------------------------------------------------------------------------------------------------------------
// C32: ACTIONX COMPDAT adds connections during the run (connection set + Jacobian pattern rebuilt, state kept)
namespace {
std::string c32ReadText(const std::string& p) { std::ifstream f(p, std::ios::binary); std::stringstream ss; ss << f.rdbuf(); return ss.str(); }
double c32WellCum(const std::vector<SummaryPoint>& S, size_t w) {
    double cum = 0, t0 = 0;
    for (const auto& p : S) { if (w < p.wells.size()) cum += p.wells[w].oil * (p.time - t0); t0 = p.time; }
    return cum;
}
}  // namespace

static void testC32ActionConnections() {
    std::printf("== C32 ACTIONX COMPDAT: new connections at run time (ACTIONX_UDQ.DATA, P2 deepened when FOPT > 1e6 stb)\n");
    const std::string path = findExample("ACTIONX_UDQ.DATA");
    if (path.empty()) { std::printf("    (examples/ACTIONX_UDQ.DATA not found - skipped)\n"); return; }
    namespace fs = std::filesystem;
    const std::string base = c32ReadText(path);
    auto replaceOnce = [](std::string s, const std::string& a, const std::string& b) { const size_t p = s.find(a); if (p != std::string::npos) s.replace(p, a.size(), b); return s; };
    // P2 completed in layer 1 only; the other producers keep layers 1-3
    std::string d = replaceOnce(base, " 'P*' 2* 1 3 OPEN 1* 1* 0.5 /\n", " 'P*' 2* 1 3 OPEN 1* 1* 0.5 /\n 'P2' 2* 2 3 SHUT 1* 1* 0.5 /\n");
    const std::string action = "ACTIONX\n DEEPEN 1 /\n FOPT > 1E6 /\n/\nCOMPDAT\n 'P2' 2* 2 3 OPEN 1* 1* 0.5 /\n/\nENDACTIO\n\n";
    // the 'SHUT' layers 2-3 above would only be a window change; for a true new connection P2 must not have them at all
    d = replaceOnce(d, " 'P2' 2* 2 3 SHUT 1* 1* 0.5 /\n", "");
    d = replaceOnce(d, " 'P*' 2* 1 3 OPEN 1* 1* 0.5 /\n", " 'P1' 2* 1 3 OPEN 1* 1* 0.5 /\n 'P3' 2* 1 3 OPEN 1* 1* 0.5 /\n 'P4' 2* 1 3 OPEN 1* 1* 0.5 /\n 'P2' 2* 1 1 OPEN 1* 1* 0.5 /\n");
    const std::string tsBlock = "TSTEP\n 10 20 30 30 24*60\n/\n";
    const std::string dA = replaceOnce(d, tsBlock, action + tsBlock);
    const std::string dir = (fs::temp_directory_path() / "ressim_c32_actx").string();
    fs::create_directories(dir);
    auto runDeck = [&](const std::string& txt, const std::string& name, std::vector<std::string>* logs, Case* cOut) {
        const std::string p = dir + "/" + name + ".DATA";
        std::ofstream(p, std::ios::binary) << txt;
        Case c = importEclipseDeck(p, nullptr, nullptr);
        if (cOut) *cOut = c;
        auto s = std::make_unique<Simulator>();
        s->log = [logs](const std::string& l) { if (logs) logs->push_back(l); };
        s->setup(c); s->initialize();
        while (s->runReportStep()) {}
        return s;
    };
    std::vector<std::string> logA;
    Case cA;
    auto A = runDeck(dA, "ACTX_A", &logA, &cA);
    const auto C = runDeck(d, "ACTX_C", nullptr, nullptr);   // no action
    size_t iP2 = 0;
    for (size_t w = 0; w < A->caseData().wells.size(); ++w) if (A->caseData().wells[w].name == "P2") iP2 = w;
    int evStep = -1; double evT = 0;
    for (const auto& e : A->actionEvents()) if (e.name == "DEEPEN") { evStep = e.step; evT = e.t; }
    bool rebuilt = false, warned = false;
    for (const auto& l : logA) { if (l.find("connection set rebuilt (2 connection(s) added)") != std::string::npos) rebuilt = true; if (l.find("cannot add connections") != std::string::npos) warned = true; }
    const double cumA = c32WellCum(A->summary(), iP2), cumC = c32WellCum(C->summary(), iP2);
    std::printf("    DEEPEN fired at %.1f d (end of report step %d); P2 cumulative oil %.6g (action) vs %.6g (layer 1 only), FOPT %.6g vs %.6g\n", evT, evStep + 1,
                cumA, cumC, A->summary().back().fopt, C->summary().back().fopt);
    check(evStep >= 0 && rebuilt && !warned, "ACTIONX COMPDAT: the connection set of P2 is rebuilt at run time (2 new connections, no 'cannot add' warning)");
    // right after the action P2 (rate controlled in both runs) needs less drawdown with three connections: higher BHP
    double bhpA = 0, bhpC = 0, oA = 0, oC = 0;
    for (const auto& p : A->summary()) if (p.time > evT + 1e-6) { bhpA = p.wells[iP2].bhp; oA = p.wells[iP2].oil; break; }
    for (const auto& p : C->summary()) if (p.time > evT + 1e-6) { bhpC = p.wells[iP2].bhp; oC = p.wells[iP2].oil; break; }
    std::printf("    first step after the action: P2 BHP %.1f vs %.1f psia (oil %.1f vs %.1f stb/d); later the deeper layers water out sooner, "
                "so the cumulative is lower\n", bhpA, bhpC, oA, oC);
    check(evStep >= 0 && (bhpA > bhpC + 5 || oA > 1.02 * oC), "the new connections take effect at once: P2 inflow performance improves (higher BHP at the rate target)");
    if (evStep < 0) return;
    // reference: the same COMPDAT written into the schedule at the report step the action acted on
    std::vector<double> ts{10, 20, 30, 30}; for (int k = 0; k < 24; ++k) ts.push_back(60);
    std::string tsA = "TSTEP\n", tsB = "TSTEP\n";
    for (int k = 0; k < int(ts.size()); ++k) (k <= evStep ? tsA : tsB) += " " + std::to_string(int(ts[size_t(k)]));
    const std::string ref = tsA + "\n/\nCOMPDAT\n 'P2' 2* 2 3 OPEN 1* 1* 0.5 /\n/\n" + tsB + "\n/\n";
    const std::string dB = replaceOnce(d, tsBlock, ref);
    const auto B = runDeck(dB, "ACTX_B", nullptr, nullptr);
    const double cumB = c32WellCum(B->summary(), iP2);
    const double fA = A->summary().back().fopt, fB = B->summary().back().fopt;
    std::printf("    reference deck (COMPDAT after report step %d): P2 cumulative %.6g (rel. diff %.2e), FOPT %.8g vs %.8g (rel. diff %.2e)\n", evStep + 1, cumB,
                std::fabs(cumA - cumB) / cumB, fB, fA, std::fabs(fA - fB) / fB);
    check(std::fabs(cumA - cumB) < 1e-6 * cumB && std::fabs(fA - fB) < 1e-6 * fB,
          "run-time connections = the same COMPDAT written into the schedule (P2 cumulative and FOPT agree to 1e-6)");
    // restart after the action: the replay rebuilds the connections, the saved well state matches
    const int stepR = evStep + 3;
    if (stepR < int(cA.schedule.size())) {
        Simulator a; a.log = [](const std::string&) {}; a.setup(cA); a.initialize();
        while (a.reportStep() < stepR && a.runReportStep()) {}
        const std::string blob = a.saveState();
        Simulator b; b.log = [](const std::string&) {}; b.setup(cA); b.initialize(); b.loadState(blob);
        while (b.runReportStep()) {}
        const double f1 = b.summary().back().fopt;
        std::printf("    restart after report step %d: FOPT %.10g vs continuous %.10g\n", stepR, f1, fA);
        check(std::fabs(f1 - fA) <= 1e-9 * fA, "restart after a run-time COMPDAT continues identically (connections rebuilt by the action replay)");
    }
}

// WSEGAICD item positions (OPM keyword JSON: 13 FLOW_RATE_EXPONENT, 14 VISC_EXPONENT) and the generator's record layout
static void testC32MswAicdItems() {
    std::printf("== C32 WSEGAICD items 13 / 14 (flow-rate and viscosity exponents)\n");
    std::vector<ressim::ExtraKeyword> kws(3);
    kws[0].name = "WELSEGS";  kws[0].body = " 'P1' 8030 0 1* INC 'HF-' /\n 2 3 1 1 150 0 0.2034 1.0E-4 /\n 4 4 2 2 0.1 0 0.2034 1.0E-4 /\n 5 5 3 3 0.1 0 0.2034 1.0E-4 /\n/\n";
    kws[1].name = "COMPSEGS"; kws[1].body = " 'P1' /\n 2 2 2 2 0.0 0.1 /\n 3 2 2 3 0.0 0.1 /\n/\n";
    kws[2].name = "WSEGAICD"; kws[2].body = " 'P1' 4 5 0.0008 1* 62.4 1.45 5* 3.05 0.67 /\n/\n";
    std::string err;
    const auto W = ressim::msw::parse(kws, 0, true, err);
    double x = 0, y = 0, st = 0; int dev = 0;
    if (!W.empty()) for (const auto& sg : W[0].segs) if (sg.device == 3) { ++dev; x = sg.x; y = sg.y; st = sg.strength; }
    std::printf("    %d AICD segments: strength %.4g, x %.3g, y %.3g %s\n", dev, st, x, y, err.c_str());
    check(dev == 2 && std::fabs(x - 3.05) < 1e-12 && std::fabs(y - 0.67) < 1e-12 && std::fabs(st - 0.0008) < 1e-15, "WSEGAICD: exponents read from items 13 / 14 (defaults '5*' skip items 8-12)");
}
