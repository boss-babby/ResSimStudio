// C42 tests: geomodel import with column assignment (tables, guesses, typed builders, merging, manual wells),
// 2-D SEG-Y lines (read / write, 2-D vs 3-D detection, tracking, intersections, gridding)
#include "ressim/geo_import.h"
#include <cstdio>
#include <fstream>

static void testC42Import() {
    using namespace ressim::geo;
    std::printf("== C42 import: delimited tables, column guesses, builders\n");
    // (1) odd header names, semicolon + decimal comma, units row, quoted names with blanks
    const std::string heads =
        "# exported by some tool\n"
        "Bohrung;Rechtswert (m);Hochwert (m);Elevation;Total Depth\n"
        "-;m;m;m;m\n"
        "\"A 1\";450000,5;6780000,25;25,0;2100\n"
        "B-2;450800;6780400;27,5;-999.25\n"
        "C-3;451300,0;6781000,0;24;2300\n";
    TableFormat F; std::string err;
    TextTable T = parseTextTable(heads, F, err);
    ColumnMap M = guessColumns(T, F, ImportKind::WellHeads);
    const bool sniff = err.empty() && T.delimUsed == ';' && T.decimalCommaUsed && T.hasHeader && T.ncol == 5 && T.rows.size() == 3 && T.meta.count("UNITS");
    // "Bohrung" / "Rechtswert" / "Hochwert" (German) are not synonyms: Elevation -> KB and Total Depth -> TD by name, then the
    // positional fallback takes the first text column as the well name and the first free numeric columns as X / Y
    const bool guess = M.role[0] == ColRole::WellName && M.role[1] == ColRole::X && M.role[2] == ColRole::Y && M.role[3] == ColRole::KB &&
                       M.role[4] == ColRole::TD && validateMap(M, ImportKind::WellHeads).empty();
    // the user overrides a role (what the mapping dialog does): swap X / Y and back
    M.role[1] = ColRole::Y; M.role[2] = ColRole::X;
    const bool swapped = buildImport(T, F, M, ImportKind::WellHeads, ImportOptions()).wells[0].y == 450000.5;
    M.role[1] = ColRole::X; M.role[2] = ColRole::Y;
    ImportResult R = buildImport(T, F, M, ImportKind::WellHeads, ImportOptions());
    const bool built = R.wells.size() == 3 && R.wells[0].name == "A 1" && std::fabs(R.wells[0].x - 450000.5) < 1e-9 && std::fabs(R.wells[1].kb - 27.5) < 1e-9 &&
                       std::fabs(R.wells[0].td - 2100) < 1e-9 && R.wells[1].td == 0;
    check(sniff, "C42: table sniffing - ';' delimiter, decimal comma, header + units row, quoted names");
    check(guess && swapped && built, "C42: unknown header names -> synonym + positional guess, manual override, 3 heads (KB, TD, null TD)");

    // (2) synonyms: Easting / Northing / KB elev / UWI, comma separated, unit factor ft -> m
    TextTable T2 = parseTextTable("UWI,Easting,Northing,KB_Elev_ft\nW1,1000,2000,100\nW2,1500,2500,90\n", F, err);
    ColumnMap M2 = guessColumns(T2, F, ImportKind::WellHeads);
    M2.scale[3] = 0.3048;
    ImportResult R2 = buildImport(T2, F, M2, ImportKind::WellHeads, ImportOptions());
    check(M2.role[0] == ColRole::WellName && M2.role[1] == ColRole::X && M2.role[2] == ColRole::Y && M2.role[3] == ColRole::KB && R2.wells.size() == 2 &&
              std::fabs(R2.wells[0].kb - 30.48) < 1e-9, "C42: synonym guesses (UWI / Easting / Northing / KB_Elev) and a per-column unit factor");

    // (3) Petrel header block (well tops with Z elevation and MD)
    const std::string ptops =
        "# Petrel well tops\nVERSION 2\nBEGIN HEADER\nX\nY\nZ\nMD\nWell\nSurface\nEND HEADER\n"
        "450000 6780000 -1975 2000 \"A 1\" \"Top Res\"\n450000 6780000 -2075 2100 \"A 1\" Base\n";
    TextTable T3 = parseTextTable(ptops, F, err);
    ColumnMap M3 = guessColumns(T3, F, ImportKind::Tops);
    ImportResult R3 = buildImport(T3, F, M3, ImportKind::Tops, ImportOptions());
    check(T3.hasHeader && T3.ncol == 6 && M3.find(ColRole::WellName) == 4 && M3.find(ColRole::Surface) == 5 && M3.find(ColRole::MD) == 3 && R3.tops.size() == 2 &&
              R3.tops[0].surface == "Top Res" && R3.tops[0].hasMd && std::fabs(R3.tops[0].tvdss - 1975) < 1e-9, "C42: Petrel BEGIN/END HEADER tops (quoted names, Z elevation -> TVDSS)");

    // (4) one deviation file with several wells (tab separated), sorted by MD per well
    const std::string dev = "WELL\tMD\tINCL\tAZIM\nA 1\t0\t0\t0\nA 1\t1000\t0\t0\nA 1\t2000\t60\t90\nB-2\t0\t0\t0\nB-2\t1500\t30\t180\nA 1\t1500\t30\t90\n";
    TextTable T4 = parseTextTable(dev, F, err);
    ColumnMap M4 = guessColumns(T4, F, ImportKind::Deviation);
    ImportResult R4 = buildImport(T4, F, M4, ImportKind::Deviation, ImportOptions());
    check(T4.delimUsed == '\t' && R4.trajMode == 0 && R4.wells.size() == 2 && R4.wells[0].traj.size() == 4 && R4.wells[0].traj[2].md == 1500 && R4.wells[1].traj.size() == 2,
          "C42: multi-well deviation file (tab), grouped per well and sorted by MD");

    // (5) headerless survey X Y TVDSS for a default well, points with set names, fault sticks without stick ids
    TextTable T5 = parseTextTable("450000 6780000 -25\n450000 6780000 1000\n450100 6780000 1100\n", F, err);
    ColumnMap M5; M5.role = {ColRole::X, ColRole::Y, ColRole::TVDSS}; M5.scale = {1, 1, 1};
    ImportOptions O5; O5.defaultWell = "D-4";
    ImportResult R5 = buildImport(T5, F, M5, ImportKind::Deviation, O5);
    TextTable T6 = parseTextTable("x y z name\n1 2 -100 TopA\n2 3 -110 TopA\n5 5 -200 TopB\n", F, err);
    ColumnMap M6 = guessColumns(T6, F, ImportKind::Points);
    ImportResult R6 = buildImport(T6, F, M6, ImportKind::Points, ImportOptions());
    TextTable T7 = parseTextTable("100 0 1000 F1\n100 0 1500 F1\n200 0 1000 F1\n210 0 1600 F1\n", F, err);
    ColumnMap M7 = guessColumns(T7, F, ImportKind::FaultSticks);
    ImportResult R7 = buildImport(T7, F, M7, ImportKind::FaultSticks, ImportOptions());
    check(R5.trajMode == 1 && R5.wells.size() == 1 && R5.wells[0].name == "D-4" && R5.wells[0].trajSrc == 1 &&
              R6.points.size() == 2 && R6.points[0].second.size() == 2 && std::fabs(R6.points[0].second[0][2] - 100) < 1e-9 &&
              R7.faults.size() == 1 && R7.faults[0].name == "F1" && R7.faults[0].sticks.size() == 2,
          "C42: headerless X/Y/TVDSS survey, named point sets (Z elevation -> depth), fault sticks split at depth resets");
    check(!validateMap(ColumnMap{{ColRole::MD, ColRole::Inc}, {1, 1}}, ImportKind::Deviation).empty() &&
              !validateMap(ColumnMap{{ColRole::X, ColRole::X, ColRole::Z}, {1, 1, 1}}, ImportKind::Points).empty(),
          "C42: map validation (incomplete MD/INC/AZI, a role assigned twice)");
}

static void testC42Wells() {
    using namespace ressim::geo;
    std::printf("== C42 wells: merging, name matching, manual wells, XYZ / offset surveys\n");
    std::vector<GeoWell> P;
    TableFormat F; std::string err;
    ImportResult H = buildImport(parseTextTable("name x y kb\nA1 1000 2000 30\nA10 3000 2000 25\n", F, err), F,
                                 guessColumns(parseTextTable("name x y kb\nA1 1000 2000 30\n", F, err), F, ImportKind::WellHeads), ImportKind::WellHeads, ImportOptions());
    MergeStats s1 = mergeHeads(P, H);
    const bool heads = s1.added == 2 && P.size() == 2;
    // file-name matching: "A1_survey" -> A1 (not A10), "a10-dev" -> A10, unknown -> -1
    const bool match = matchWell(P, "A1_survey") == 0 && matchWell(P, "a10-dev") == 1 && matchWell(P, "A100") == -1 && matchWell(P, "a1") == 0;
    // MD/INC/AZI survey for A1 (from a file without a well column -> defaultWell from the file name)
    TextTable T = parseTextTable("MD INC AZI\n0 0 0\n1000 0 0\n1500 30 90\n", F, err);
    ImportOptions O; O.defaultWell = "A1_survey";
    ImportResult D = buildImport(T, F, guessColumns(T, F, ImportKind::Deviation), ImportKind::Deviation, O);
    MergeStats s2 = mergeSurveys(P, D);
    const P3 e = P[0].at(1500);
    const bool surv = s2.updated == 1 && P.size() == 2 && P[0].traj.size() == 3 && e[0] > 1000 + 100 && std::fabs(e[1] - 2000) < 1e-6 && std::fabs(P[0].at(0)[2] + 30) < 1e-9;
    // the head moves afterwards: the MD/INC/AZI trajectory follows
    ImportResult H2; GeoWell h; h.name = "a1"; h.x = 1100; h.y = 2000; h.kb = 40; H2.wells = {h};
    mergeHeads(P, H2);
    const bool follow = std::fabs(P[0].at(1500)[0] - (e[0] + 100)) < 1e-6 && std::fabs(P[0].at(0)[2] + 40) < 1e-9;
    // manual well: vertical to TD, then an XYZ survey typed in a table (no MD -> arc length) gives inclinations
    GeoWell m; m.name = "NEW-1"; m.x = 500; m.y = 500; m.kb = 20;
    setVertical(m, 2000);
    const bool vert = m.traj.size() == 2 && std::fabs(m.at(2000)[2] - 1980) < 1e-9;
    std::vector<Station> st(3);
    st[0].md = -1; st[0].x = 500; st[0].y = 500; st[0].tvd = 980;
    st[1].md = -1; st[1].x = 500; st[1].y = 500; st[1].tvd = 1480;
    st[2].md = -1; st[2].x = 1000; st[2].y = 500; st[2].tvd = 1480;
    const bool xyz = trajFromXyz(m, st, err) && std::fabs(m.traj[0].md - 1000) < 1e-9 && std::fabs(m.traj[2].md - 2000) < 1e-9 &&
                     std::fabs(m.traj[2].inc - 90) < 1e-6 && std::fabs(m.traj[2].azi - 90) < 1e-6 && m.trajSrc == 1;
    // DX / DY offset survey + TVD below KB for a well without a head in the project -> created at (0,0) with a warning
    TextTable T3 = parseTextTable("well,md,tvd,ns,ew\nZ9,0,0,0,0\nZ9,1000,1000,0,0\nZ9,1200,1190,10,15\n", F, err);
    ColumnMap M3 = guessColumns(T3, F, ImportKind::Deviation);
    ImportResult D3 = buildImport(T3, F, M3, ImportKind::Deviation, ImportOptions());
    MergeStats s3 = mergeSurveys(P, D3);
    const bool offs = M3.find(ColRole::DX) == 4 && M3.find(ColRole::DY) == 3 && D3.trajMode == 2 && D3.tvdBelowKb && s3.added == 1 && !s3.warnings.empty() &&
                      std::fabs(P.back().at(1200)[0] - 15) < 0.5 && std::fabs(P.back().at(1200)[1] - 10) < 0.5;
    // tops by TVDSS on the deviated well, replacement of an existing top, skipped unknown well
    std::vector<TopRow> tops(3);
    tops[0].well = "A1"; tops[0].surface = "Top"; tops[0].tvdss = 900; tops[0].hasTvd = true;
    tops[1].well = "A1"; tops[1].surface = "TOP"; tops[1].md = 950; tops[1].hasMd = true;
    tops[2].well = "XX"; tops[2].surface = "Top"; tops[2].md = 1; tops[2].hasMd = true;
    MergeStats s4 = mergeTops(P, tops, false);
    const bool tp = P[0].tops.size() == 1 && std::fabs(P[0].tops[0].md - 950) < 1e-9 && s4.warnings.size() == 1;
    MergeStats s5 = mergeTops(P, {tops[0]}, false);
    const bool tp2 = std::fabs(P[0].tops[0].md - 940) < 0.5 && s5.updated == 1;   // TVDSS 900 + KB 40 on the vertical part
    check(heads && match, "C42: head merge, file-name / case-insensitive well matching (A1_survey -> A1, not A10)");
    check(surv && follow, "C42: survey merged by file name; the trajectory follows a corrected head (x and KB)");
    check(vert && xyz, "C42: manual well - vertical to TD, typed XYZ survey with MD from arc length and inclinations 0 / 90 deg");
    check(offs, "C42: NS / EW offset survey with TVD below KB creates a missing well (warning) and positions it");
    check(tp && tp2, "C42: tops by MD / TVDSS, case-insensitive replacement, unknown wells reported");
    check(writeDeviationTable(P[0]).find("MD\tINC\tAZI") == 0, "C42: deviation table export header");
}

static void testC42Seismic2D() {
    using namespace ressim::geo;
    std::printf("== C42 2-D seismic: SEG-Y lines, detection, tracking, intersections, gridding\n");
    // two crossing lines with a dipping reflector (Ricker-like wavelet), t(x) = 800 + 0.05 (x - 1000)
    auto makeLine = [](const char* name, double xa, double ya, double xb, double yb, int n) {
        SeismicLine L; L.name = name; L.ns = 300; L.dt = 4; L.t0 = 0;
        for (int t = 0; t < n; ++t) {
            const double f = double(t) / (n - 1), x = xa + f * (xb - xa), y = ya + f * (yb - ya);
            L.x.push_back(x); L.y.push_back(y); L.cdp.push_back(100 + t);
            const double tr = 800 + 0.05 * (x - 1000) + 0.02 * (y - 1000);
            for (int s = 0; s < L.ns; ++s) { const double u = (s * L.dt - tr) / 12.0; L.a.push_back(float((1 - 2 * u * u) * std::exp(-u * u))); }
        }
        return L;
    };
    const SeismicLine A = makeLine("L1", 0, 1000, 2000, 1000, 101), B = makeLine("L2", 1000, 0, 1000, 2000, 81);
    const std::string pa = "/tmp/c42_l1.sgy", pb = "/tmp/c42_l2.sgy";
    const bool wr = writeSegy2D(pa, A) && writeSegy2D(pb, B);
    std::string err, info;
    SeismicLine A2 = readSegy2D(pa, SegyByteMap(), err);
    const bool rd = err.empty() && A2.ntr() == 101 && A2.ns == 300 && std::fabs(A2.x[50] - 1000) < 0.01 && A2.cdp[0] == 100 && std::fabs(A2.at(50, 200) - A.at(50, 200)) < 1e-6;
    const int k2 = segyGeometryKind(pa, SegyByteMap(), info, err);
    // a 3-D cube written by writeSegy must be detected as 3-D
    SeismicCube C; C.ni = 6; C.nj = 5; C.ns = 20; C.a.assign(size_t(6 * 5 * 20), 0.5f);
    const bool w3 = writeSegy("/tmp/c42_c3.sgy", C);
    const int k3 = segyGeometryKind("/tmp/c42_c3.sgy", SegyByteMap(), info, err);
    // wrong coordinate byte location: X = Y = 0 everywhere -> a clear message
    SegyByteMap bad; bad.xByte = 41; bad.yByte = 45; bad.fallbackXY = false;
    std::string e2; readSegy2D(pa, bad, e2);
    check(wr && rd && k2 == 2 && w3 && k3 == 3 && e2.find("byte locations") != std::string::npos,
          "C42: 2-D SEG-Y write / read (CDP, X/Y, samples), 2-D vs 3-D detection, wrong byte location reported");
    // tracking a peak from one seed along each line
    std::vector<SeismicLine> lines = {A2, readSegy2D(pb, SegyByteMap(), err)};
    LineHorizon H; H.name = "H1";
    double maxErr = 0; int picked = 0;
    for (size_t l = 0; l < lines.size(); ++l) {
        const auto& L = lines[l];
        const int s0 = L.ntr() / 2;
        const double tSeed = 800 + 0.05 * (L.x[size_t(s0)] - 1000) + 0.02 * (L.y[size_t(s0)] - 1000) + 6;   // 6 ms off: snapped
        H.t.push_back(trackLine(L, {{s0, tSeed}}, 0, 12, 8));
        for (int t = 0; t < L.ntr(); ++t) {
            const double tt = H.t.back()[size_t(t)];
            if (!std::isfinite(tt)) continue;
            ++picked; maxErr = std::max(maxErr, std::fabs(tt - (800 + 0.05 * (L.x[size_t(t)] - 1000) + 0.02 * (L.y[size_t(t)] - 1000))));
        }
    }
    std::printf("   tracked %d / %d traces, max |error| %.3f ms\n", picked, lines[0].ntr() + lines[1].ntr(), maxErr);
    check(picked == lines[0].ntr() + lines[1].ntr() && maxErr < 0.5, "C42: 1-D seeded tracking along both lines (sub-sample, all traces, < 0.5 ms)");
    const auto ties = lineIntersections(lines);
    double mis = 1e9;
    if (ties.size() == 1) mis = std::fabs(H.t[0][size_t(ties[0].ta)] - H.t[1][size_t(ties[0].tb)]);
    check(ties.size() == 1 && std::fabs(ties[0].x - 1000) < 1e-6 && std::fabs(ties[0].y - 1000) < 1e-6 && ties[0].ta == 50 && ties[0].tb == 40 && mis < 0.5,
          "C42: line intersection at (1000, 1000), traces 50 / 40, mis-tie < 0.5 ms");
    const auto ip = interpolatePicks(10, {{2, 100}, {6, 140}, {9, std::nan("")}});
    check(!std::isfinite(ip[1]) && std::fabs(ip[4] - 120) < 1e-9 && std::fabs(ip[6] - 140) < 1e-9 && !std::isfinite(ip[7]), "C42: manual picks interpolated between clicked traces");
    // picks -> depth points -> grid on an automatic lattice
    VelocityModel v; v.v0 = 2000; v.k = 0;
    const auto pts = lineHorizonPoints(lines, H, v, 2);
    Surface lat = latticeFor(pts, 21);
    ressim::geostat::Variogram vg; vg.rangeMax = vg.rangeMin = 3000;
    Surface s = gridPoints(pts, lat, vg);
    const double zc = s.at(1000, 1000), zt = timeToDepth(v, 800);
    std::printf("   %zu points, lattice %d x %d, depth at the tie %.2f (exact %.2f)\n", pts.size(), lat.nx, lat.ny, zc, zt);
    check(lat.valid() && pts.size() > 80 && std::fabs(zc - zt) < 2.0, "C42: line picks -> depth points -> kriged surface on an automatic lattice");
    std::remove(pa.c_str()); std::remove(pb.c_str()); std::remove("/tmp/c42_c3.sgy");
}

static void testC42Calc() {
    std::printf("== C42 surface / property calculator, statistics\n");
    using namespace ressim::geo;
    Surface A; A.name = "Top Ness"; A.nx = 5; A.ny = 4; A.x0 = 0; A.y0 = 0; A.dx = 100; A.dy = 100;
    A.z.assign(20, 0.0);
    for (int j = 0; j < 4; ++j) for (int i = 0; i < 5; ++i) A.z[size_t(i + 5 * j)] = 2000 + 10 * i;
    Surface B = A; B.name = "2 base"; B.dx = 50; B.dy = 50; B.nx = 9; B.ny = 7;   // finer lattice, smaller extent in y? covers 0..400 x 0..300
    B.z.assign(63, 0.0);
    for (int j = 0; j < 7; ++j) for (int i = 0; i < 9; ++i) B.z[size_t(i + 9 * j)] = 2100 + 0.1 * (50.0 * i);   // 2100 + 0.1 x
    Surface out; std::string err;
    const bool ok = surfaceCalculator({A, B}, 0, "_2_BASE - TOP_NESS", out, err);
    // thickness at node (i,j) = 2100 + 0.1*100 i - (2000 + 10 i) = 100
    bool thick = ok && out.valid() && out.nx == 5;
    for (double z : out.z) thick = thick && std::fabs(z - 100.0) < 1e-9;
    Surface o2;
    const bool ok2 = surfaceCalculator({A, B}, 0, "if(X > 150, TOP_NESS + 5, nan)", o2, err);
    const bool cond = ok2 && std::isnan(o2.z[0]) && std::isnan(o2.z[1]) && std::fabs(o2.z[2] - 2025) < 1e-9;
    std::string e3;
    Surface o3;
    const bool bad = !surfaceCalculator({A, B}, 0, "TOP_NESS + NOPE", o3, e3) && !e3.empty();
    check(exprName("Top Ness") == "TOP_NESS" && exprName("2 base") == "_2_BASE" && thick && cond && bad,
          "C42: surface calculator - resampling to the target lattice (isochore = 100), if()/X, nan masking, unknown-name error");
    std::vector<double> v = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, NAN};
    const Stats s = describe(v);
    check(s.n == 10 && s.nNull == 1 && s.min == 1 && s.max == 10 && std::fabs(s.mean - 5.5) < 1e-12 && std::fabs(s.p50 - 5.5) < 1e-12 &&
          std::fabs(s.p10 - 1.9) < 1e-12 && std::fabs(s.p90 - 9.1) < 1e-12 && std::fabs(s.sd - 3.0276503540974917) < 1e-9,
          "C42: describe() - n / nulls / min / max / mean / sd / P10 / P50 / P90");
    // property calculator path = evalExpression over named grid vectors
    std::vector<double> poro = {0.1, 0.2, 0.3}, fac = {0, 1, 1}, r;
    std::map<std::string, const std::vector<double>*> vars = {{"PORO", &poro}, {"FACIES", &fac}};
    const bool pc = evalExpression("if(FACIES == 1, PORO * 0.9, 0)", vars, 3, r, err) && r[0] == 0 &&
                    std::fabs(r[1] - 0.18) < 1e-12 && std::fabs(r[2] - 0.27) < 1e-12;
    check(pc, "C42: property calculator expression - facies-conditional porosity edit");
}

// C42 CECON extensions: minimum completion rates (items 11 / 12) and depth-ordered PLUG
static void testC42Cecon() {
    std::printf("== C42 CECON minimum rates and PLUG\n");
    auto run = [](const std::string& body, std::string& log, RunStats& st, double& fopt, double& qo) {
        Case c = smallCase();
        ExtraKeyword k; k.name = "CECON"; k.section = DeckSection::Schedule; k.step = 0; k.body = body;
        c.extraKeywords.push_back(k);
        Simulator s; log.clear(); s.log = [&](const std::string& m) { log += m + "\n"; }; s.setup(c); s.initialize();
        while (s.runReportStep()) {}
        st = s.stats(); fopt = s.summary().back().fopt; qo = s.summary().back().fopr;
        size_t p = 0; int shown = 0;
        while ((p = log.find("CECON", p)) != std::string::npos && shown < 3) {
            const size_t a = log.rfind('\n', p), b = log.find('\n', p);
            std::printf("    | %s\n", log.substr(a == std::string::npos ? 0 : a + 1, b - (a == std::string::npos ? 0 : a + 1)).c_str());
            p = b; ++shown;
        }
    };
    const RunResult r0 = runSingle(smallCase());
    std::string L1, L2, L3; RunStats s1, s2, s3; double f1, q1, f2, q2, f3, q3;
    // item 11 on the lower connection: always below a huge minimum -> that connection closes, workover WELL is not applied
    run(" 'PROD' 0 0 2 2 1* 1* 1* 'WELL' 'NO' 1e9 /\n", L1, s1, f1, q1);
    std::printf("    min oil rate on K=2: shut-ins %d, FOPT %.6g (reference %.6g), final oil rate %.4g\n", s1.econShut, f1, r0.summary.back().fopt, q1);
    check(s1.econShut == 0 && q1 > 0 && L1.find("oil rate") != std::string::npos && L1.find("workover CON closes 1") != std::string::npos,
          "C42: CECON item 11 - a completion below the minimum oil rate is closed (CON), the well keeps producing");
    // water-cut limit on the upper connection: PLUG closes it and every deeper one -> no connection left -> well shut;
    // CON closes only the upper one
    run(" 'PROD' 0 0 1 1 1e-6 1* 1* 'PLUG' 'NO' /\n", L2, s2, f2, q2);
    run(" 'PROD' 0 0 1 1 1e-6 1* 1* 'CON' 'NO' /\n", L3, s3, f3, q3);
    std::printf("    WCT limit on K=1: PLUG shut-ins %d (oil %.4g), CON shut-ins %d (oil %.4g)\n", s2.econShut, q2, s3.econShut, q3);
    const bool fired = L3.find("CECON limit") != std::string::npos;
    check(fired && s2.econShut == 1 && q2 < 0.05 * q3 && L2.find("workover PLUG closes 2") != std::string::npos && s3.econShut == 0 && q3 > 0,
          "C42: CECON PLUG - a water violation closes the completion and all deeper ones (CON only the offender)");
}
