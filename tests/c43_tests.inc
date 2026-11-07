// C43 tests: geomodel tools (mis-tie correction, polygons, polygon / per-segment volumetrics, layer-cake velocity
// model, checkshot fit, seismic calculator and cropping)
#include "ressim/geo_tools.h"
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <cstring>

static void testC43Mistie() {
    std::printf("== C43 2-D mis-tie correction\n");
    using namespace ressim::geo;
    auto mkLine = [](const std::string& nm, double x0, double y0, double x1, double y1, int n) {
        SeismicLine L; L.name = nm; L.ns = 10; L.dt = 4;
        for (int i = 0; i < n; ++i) { L.x.push_back(x0 + (x1 - x0) * i / (n - 1)); L.y.push_back(y0 + (y1 - y0) * i / (n - 1)); L.cdp.push_back(i + 1); }
        L.a.assign(size_t(n) * 10, 0.0f);
        return L;
    };
    std::vector<SeismicLine> lines = {mkLine("L0", 0, 500, 1000, 500, 101), mkLine("L1", 300, 0, 300, 1000, 101),
                                      mkLine("L2", 700, 0, 700, 1000, 101), mkLine("L3", 2000, 0, 2000, 1000, 51)};
    const double err[4] = {0, 8, -5, 3};
    std::vector<LineHorizon> hz(2);
    for (int h = 0; h < 2; ++h) {
        hz[size_t(h)].name = h ? "B" : "A";
        for (size_t l = 0; l < lines.size(); ++l) {
            std::vector<double> t;
            for (int k = 0; k < lines[l].ntr(); ++k) t.push_back(1000 + 200 * h + 0.05 * lines[l].x[size_t(k)] + 0.02 * lines[l].y[size_t(k)] + err[l]);
            hz[size_t(h)].t.push_back(t);
        }
    }
    const auto ties = lineIntersections(lines);
    const MistieFit F = fitMisties(lines, ties, hz, 0);
    std::printf("    ties %zu, pairs %d, rms %.3f -> %.3f ms, shifts %.3f %.3f %.3f %.3f\n", ties.size(), F.pairs, F.rmsBefore, F.rmsAfter,
                F.shift[0], F.shift[1], F.shift[2], F.shift[3]);
    // the picks of the crossing lines are taken at the nearest trace -> sub-trace position error of a few hundredths of a ms
    check(ties.size() == 2 && F.pairs == 4 && std::fabs(F.shift[0]) < 1e-12 && std::fabs(F.shift[1] + 8) < 0.2 && std::fabs(F.shift[2] - 5) < 0.2 &&
          F.shift[3] == 0 && F.component[3] == -1 && F.rmsAfter < 0.1 * F.rmsBefore,
          "C43: least-squares bulk shift per line removes the mis-ties (reference line fixed, unconnected line untouched)");
    std::vector<SeismicLine> l2 = lines; std::vector<LineHorizon> h2 = hz;
    applyLineShifts(l2, h2, F.shift);
    check(std::fabs(l2[1].t0 - lines[1].t0 - F.shift[1]) < 1e-12 && std::fabs(h2[0].t[1][5] - hz[0].t[1][5] - F.shift[1]) < 1e-12,
          "C43: shifts applied to the line data (t0) and to every pick");
}

static void testC43Polygons() {
    std::printf("== C43 polygons and polygon / per-segment volumetrics\n");
    using namespace ressim::geo;
    const std::string txt = "# Lease A\n0 0\n100 0\n100 100\n0 100\n0 0\n999 999\n# Path\n0 0 0\n50 50 0\n999 999\n"
                            "10 10 0 P1\n20 10 0 P1\n20 20 0 P1\n50 50 0 P2\n60 50 0 P2\n60 60 0 P2\n";
    const auto P = parsePolygons(txt);
    const bool parsed = P.size() == 4 && P[0].name == "Lease A" && P[0].closed && P[0].pts.size() == 4 && !P[1].closed && P[1].name == "Path" &&
                        P[2].name == "P1" && P[3].name == "P2" && P[3].pts.size() == 3;
    const bool geom = std::fabs(polygonArea(P[0]) - 10000) < 1e-9 && insidePolygon(P[0], 50, 50) && !insidePolygon(P[0], 150, 50) &&
                      !insidePolygon(P[1], 10, 10) && std::fabs(polylineLength(P[1]) - std::sqrt(5000.0)) < 1e-9;
    writePolygons("/tmp/c43_poly.txt", P);
    std::string e;
    const auto Q = readPolygons("/tmp/c43_poly.txt", e);
    const bool rt = Q.size() == 4 && Q[0].pts.size() == 4 && Q[0].name == "Lease A" && !Q[1].closed;
    Surface s; s.nx = 11; s.ny = 11; s.dx = s.dy = 20; s.z.assign(121, 1.0);
    const Surface c = clipSurface(s, P[0], true);
    int def = 0; for (double z : c.z) def += std::isfinite(z) ? 1 : 0;   // nodes 0..100 in x and y (6 x 6, edges included / excluded by the rule)
    check(parsed && geom && rt && def >= 25 && def <= 36, "C43: polygon parsing (names, 999 separators, id column), area / inside / length, file round trip, surface clipping");
    // grid with two flat surfaces: polygon mask and per-block contacts in the volumetrics
    Surface top; top.name = "T"; top.nx = 21; top.ny = 21; top.dx = top.dy = 100; top.z.assign(441, 1000.0);
    Surface bot = top; bot.name = "B"; for (double& z : bot.z) z = 1100.0;
    GridSpec gs; gs.dx = gs.dy = 100; gs.nx = 20; gs.ny = 20; gs.tieToWells = false;
    Zone z; z.name = "Z"; z.top = 0; z.base = 1; z.nLayers = 10; gs.zones = {z};
    const GeoGrid g0 = buildGrid(gs, {top, bot}, {}, {}, e);
    std::vector<double> poro(g0.bulk.size(), 0.2), ntg(g0.bulk.size(), 1.0);
    VolSpec vs; vs.satHeight = false; vs.swConst = 0.0; vs.bo = 1.0; vs.contacts.owc = 1050;
    const VolResult all = volumetrics(g0, poro, ntg, nullptr, vs);
    Polygon half; half.pts = {{0, 0}, {1000, 0}, {1000, 2000}, {0, 2000}};
    VolSpec vm = vs; vm.mask = polygonMask(g0, half);
    const VolResult V = volumetrics(g0, poro, ntg, nullptr, vm);
    GeoGrid g1 = g0;
    g1.block.assign(size_t(gs.nx * gs.ny), 0);
    for (int j = 0; j < gs.ny; ++j) for (int i = 10; i < gs.nx; ++i) g1.block[size_t(i + gs.nx * j)] = 1;
    g1.nBlocks = 2;
    VolSpec vb = vs; vb.blockContacts.resize(2); vb.blockContacts[0].owc = 1050; vb.blockContacts[1].owc = 1100;
    const VolResult VB = volumetrics(g1, poro, ntg, nullptr, vb);
    std::printf("    STOIIP all %.6g, west half %.6g, per-block contacts %.6g (expected 1.5 x all)\n", all.total.stoiip, V.total.stoiip, VB.total.stoiip);
    check(e.empty() && all.total.stoiip > 0 && std::fabs(V.total.stoiip / all.total.stoiip - 0.5) < 1e-9 && std::fabs(VB.total.stoiip / all.total.stoiip - 1.5) < 1e-9,
          "C43: volumetrics inside a polygon (half the area) and with contacts per fault block (deeper OWC in block 2)");
}

static void testC43LayerCake() {
    std::printf("== C43 layer-cake velocity model, checkshots\n");
    using namespace ressim::geo;
    Surface t1; t1.name = "T1"; t1.nx = 41; t1.ny = 41; t1.dx = t1.dy = 50; t1.z.resize(41 * 41);
    for (int j = 0; j < 41; ++j) for (int i = 0; i < 41; ++i) t1.z[size_t(i + 41 * j)] = 1000 + 0.05 * (50.0 * i);   // ms
    Surface t2 = t1; t2.name = "T2"; for (double& t : t2.z) t += 100;
    // truth: layer 1 at 2000 m/s from the datum 0, layer 2 at 3000 m/s -> z1 = t1, z2 = z1 + 150
    std::vector<GeoWell> wells;
    const double wx[4] = {300, 1500, 1800, 600}, wy[4] = {400, 600, 1700, 1500};
    for (int w = 0; w < 4; ++w) {
        GeoWell W; W.name = "W" + std::to_string(w + 1); W.x = wx[w]; W.y = wy[w]; W.kb = 0;
        const double z1 = 1000 + 0.05 * wx[w];
        W.tops = {{"Top1", z1}, {"Top2", z1 + 150 + (w == 2 ? 6.0 : 0.0)}};   // W3: base 6 m deeper than the model
        wells.push_back(W);
    }
    LayerCakeInput in; in.twt = {t1, t2}; in.topName = {"Top1", "Top2"}; in.mode = {0, 0}; in.residual = false;
    const LayerCakeResult A = layerCakeDepth(in, wells);
    in.residual = true;
    const LayerCakeResult B = layerCakeDepth(in, wells);
    const double z2 = B.depth[1].at(1800, 1700), z2far = B.depth[1].at(1500, 600);
    std::printf("    Vint %.1f / %.1f m/s (wells %d / %d), rms at the tops %.3f -> %.3f m, base at W3 %.2f (top %.2f), at W2 %.2f\n", A.vMean[0], A.vMean[1],
                A.nWells[0], A.nWells[1], B.rmsBefore, B.rmsAfter, z2, 1000 + 0.05 * 1800 + 156, z2far);
    check(std::fabs(A.vMean[0] - 2000) < 1e-6 && std::fabs(A.vMean[1] - 3030.0) < 1e-6 && A.nWells[1] == 4 && B.rmsAfter < 0.02 * B.rmsBefore &&
          std::fabs(z2 - (1000 + 90 + 156)) < 0.05 && std::fabs(z2far - (1000 + 75 + 150)) < 0.05,
          "C43: interval velocities from time surfaces + tops (W3 +6 m raises the layer-2 mean by 30 m/s), residual correction ties the tops");
    // checkshots from v = 1800 + 0.6 z
    std::vector<double> zz, tt;
    for (double z = 100; z <= 3000; z += 100) { zz.push_back(z); tt.push_back(2000.0 / 0.6 * std::log(1 + 0.6 * z / 1800.0)); }
    double rms = 0;
    const VelocityModel vm = fitCheckshots(zz, tt, 0, &rms);
    std::printf("    checkshots: v0 %.2f k %.4f rms %.4f ms\n", vm.v0, vm.k, rms);
    check(std::fabs(vm.v0 - 1800) < 2 && std::fabs(vm.k - 0.6) < 0.002 && rms < 0.05, "C43: checkshot fit recovers v0 + k z (1800, 0.6)");
}

static void testC43Seismic() {
    std::printf("== C43 seismic calculator and cropping\n");
    using namespace ressim::geo;
    SeismicCube c; c.name = "S"; c.ni = 6; c.nj = 5; c.ns = 40; c.dt = 4; c.t0 = 0; c.x0 = 100; c.y0 = 200; c.ix = 25; c.jy = 25;
    c.a.resize(size_t(c.ni * c.nj * c.ns));
    for (int i = 0; i < c.ni; ++i) for (int j = 0; j < c.nj; ++j) for (int s = 0; s < c.ns; ++s)
        c.a[(size_t(i) * c.nj + j) * c.ns + s] = float(std::sin(0.3 * s) * (1 + i));
    SeismicCube o; std::string e;
    const bool ok = seismicCalculator(c, "if(T > 80, AMP * 2, 0) + I", 5, o, e);
    const bool val = ok && std::fabs(o.at(2, 1, 30) - (2 * c.at(2, 1, 30) + 2)) < 1e-5 && std::fabs(o.at(3, 0, 10) - 3) < 1e-6;
    SeismicCube env;
    const bool ok2 = seismicCalculator(c, "ENVELOPE - abs(AMP)", 5, env, e);
    double mn = 1e30; for (int s = 5; s < 35; ++s) mn = std::min(mn, double(env.at(1, 1, s)));
    std::string e3; SeismicCube bad;
    const bool rej = !seismicCalculator(c, "AMP + XYZ", 5, bad, e3) && !e3.empty();
    check(val && ok2 && mn > -0.05 && rej, "C43: seismic calculator - AMP / T / I and an on-demand attribute (envelope >= |amplitude|), unknown names rejected");
    const SeismicCube k = cropCube(c, 1, 4, 1, 3, 10, 29, 1, 2, 2);
    double x, y; k.xy(0, 0, x, y);
    check(k.ni == 4 && k.nj == 2 && k.ns == 10 && k.t0 == 40 && k.dt == 8 && std::fabs(x - 125) < 1e-9 && std::fabs(y - 225) < 1e-9 &&
          k.at(1, 1, 2) == c.at(2, 3, 14) && k.il0 == 2 && k.dxl == 2,
          "C43: sub-volume cropping with decimation keeps geometry, header numbering and samples");
}

// C43 P10: geomodel C API (rs_geo_*, used by python/ressim/geo.py)
#if defined(_WIN32) && !defined(RESSIM_CAPI_BUILD)
#define RESSIM_CAPI_BUILD   // static link of the engine objects into the test exe: no dllimport thunks
#endif
#include "ressim/capi.h"
static void testC43GeoCapi() {
    std::printf("== C43 geomodel C API (Python / MATLAB / Julia scripting)\n");
    rs_geo* g = rs_geo_new();
    const std::string dir = "/tmp/ressim_c43_capi";
    const int nf = rs_geo_demo(g, dir.c_str(), 0);
    const int ns = rs_geo_num_surfaces(g), nw = rs_geo_num_wells(g);
    std::printf("    demo: %d files, %d surfaces, %d wells, seismic %s\n", nf, ns, nw, rs_geo_seismic_json(g) ? rs_geo_seismic_json(g) : "-");
    check(nf > 0 && ns == 3 && nw == 6 && rs_geo_seismic_json(g) != nullptr, "C43 capi: demo project loaded (3 surfaces, 6 wells, SEG-Y)");
    // surface values / calculator / sampling
    const int n0 = rs_geo_surface_values(g, 0, nullptr, 0);
    std::vector<double> z(static_cast<size_t>(n0));
    rs_geo_surface_values(g, 0, z.data(), n0);
    const int ic = rs_geo_surface_calc(g, "TopRes + 25", 0, "TopShift");
    std::vector<double> z2(static_cast<size_t>(n0));
    rs_geo_surface_values(g, ic, z2.data(), n0);
    double dmax = 0; int nd = 0;
    for (int k = 0; k < n0; ++k) if (std::isfinite(z[size_t(k)])) { dmax = std::max(dmax, std::fabs(z2[size_t(k)] - z[size_t(k)] - 25)); ++nd; }
    check(ic == 3 && nd > 0 && dmax < 1e-9, "C43 capi: surface calculator result added (TopRes + 25)");
    // well JSON + position + log
    const std::string wj = rs_geo_well_json(g, 0);
    double md[2] = {0, 500}, xyz[6];
    const int np = rs_geo_well_position(g, 0, md, 2, xyz);
    const int nl = rs_geo_well_log(g, 0, "PHIE", nullptr, nullptr, 0);
    std::printf("    well 0: %s ... ; log PHIE %d samples; tvd(500) %.2f\n", wj.substr(0, 60).c_str(), nl, xyz[5]);
    check(wj.find("\"tops\":[{") != std::string::npos && np == 2 && nl > 10 && xyz[5] > 0, "C43 capi: wells with tops, positions and logs");
    // polygon + clip
    char sj[512]; std::snprintf(sj, sizeof sj, "%s", rs_geo_surface_json(g, 0));
    double x0, y0, dx, dy; int nx, ny;
    std::sscanf(std::strstr(sj, "\"nx\":"), "\"nx\":%d,\"ny\":%d,\"x0\":%lf,\"y0\":%lf,\"dx\":%lf,\"dy\":%lf", &nx, &ny, &x0, &y0, &dx, &dy);
    const double xm = x0 + 0.5 * (nx - 1) * dx;
    const double pts[8] = {x0 - 1, y0 - 1, xm, y0 - 1, xm, y0 + ny * dy, x0 - 1, y0 + ny * dy};
    const int ip = rs_geo_add_polygon(g, "west", 1, pts, 4);
    const int left = rs_geo_clip_surface(g, ic, ip, 1);
    std::printf("    clip: %d of %d nodes kept\n", left, nd);
    check(ip == 0 && left > 0.4 * nd && left < 0.6 * nd, "C43 capi: clip a surface by a polygon");
    check(rs_geo_write_polygons(g, (dir + "/poly.txt").c_str()) == 1 && rs_geo_load_polygons(g, (dir + "/poly.txt").c_str()) == 1 &&
          rs_geo_num_polygons(g) == 2, "C43 capi: polygons written and read back");
    // seismic calculator + crop
    float tr[2000];
    const int nsmp = rs_geo_seismic_trace(g, 3, 4, tr, 2000);
    const float a0 = tr[nsmp / 2];
    rs_geo_seismic_calc(g, "AMP * 2", 6);
    rs_geo_seismic_trace(g, 3, 4, tr, 2000);
    const int ntr = rs_geo_seismic_crop(g, 0, 9, 0, 9, 0, -1 + nsmp, 2, 2, 1);
    check(nsmp > 0 && std::fabs(tr[nsmp / 2] - 2 * a0) < 1e-5f * (1 + std::fabs(a0)) && ntr == 25, "C43 capi: seismic calculator and crop");
    check(rs_geo_surface_calc(g, "NOPE + 1", 0, "x") == -1 && std::strlen(rs_last_error()) > 0, "C43 capi: errors reported through rs_last_error");
    rs_geo_free(g);
}

// C43 well performance events WPWE0-7 (OPM data::WellPerformanceEvents): CECON workovers on the two-connection producer
// of smallCase (K = 1 top, K = 2 bottom), and the WTEST 'E' retest of a well whose connections were all closed
static void testC43Wpwe() {
    std::printf("== C43 WPWE well performance events\n");
    struct Out { int sum[8] = {}; double tFirst[8]; int maxv[8] = {}; std::string log; };
    auto run = [](const std::string& cecon, const std::string& wtest) {
        Out o; for (double& t : o.tFirst) t = -1;
        Case c = smallCase();
        ExtraKeyword k; k.name = "CECON"; k.section = DeckSection::Schedule; k.step = 0; k.body = cecon; c.extraKeywords.push_back(k);
        if (!wtest.empty()) { ExtraKeyword t; t.name = "WTEST"; t.section = DeckSection::Schedule; t.step = 0; t.body = wtest; c.extraKeywords.push_back(t); }
        Simulator s; s.log = [&](const std::string& m) { o.log += m + "\n"; }; s.setup(c); s.initialize();
        while (s.runReportStep()) {}
        int pw = -1;
        for (size_t w = 0; w < c.wells.size(); ++w) if (c.wells[w].name == "PROD") pw = int(w);
        for (const SummaryPoint& p : s.summary())
            if (pw >= 0 && size_t(pw) < p.wells.size())
                for (int e = 0; e < 8; ++e) { const int v = p.wells[size_t(pw)].pwe[e]; o.sum[e] += v; o.maxv[e] = std::max(o.maxv[e], v); if (v && o.tFirst[e] < 0) o.tFirst[e] = p.time; }
        std::printf("    sums WPWE0..7: %d %d %d %d %d %d %d %d  (first WPWE2 at t = %.4g d)\n", o.sum[0], o.sum[1], o.sum[2], o.sum[3], o.sum[4], o.sum[5], o.sum[6], o.sum[7], o.tFirst[2]);
        return o;
    };
    // '+CON' on the top completion: offender counted (WPWE2 = 1), the swept bottom one not; reach -> closed to bottom; no
    // connection left -> the simulator shuts the well (WPWE7)
    const Out a = run(" 'PROD' 0 0 1 1 1e-6 1* 1* '+CON' 'NO' /\n", "");
    check(a.sum[2] == 1 && a.sum[3] == 1 && a.sum[7] == 1 && a.tFirst[2] == a.tFirst[3] && a.tFirst[3] == a.tFirst[7] && a.sum[1] == 0 && a.sum[0] == 0,
          "C43: WPWE '+CON' - only the offender counts as closed (WPWE2 = 1), closed to bottom (WPWE3) and shut (WPWE7) in the same step");
    // CON on the top completion: the bottom one keeps flowing -> not closed to bottom, no status event
    const Out b = run(" 'PROD' 0 0 1 1 1e-6 1* 1* 'CON' 'NO' /\n", "");
    check(b.sum[2] == 1 && b.sum[3] == 0 && b.sum[7] == 0, "C43: WPWE CON on the top completion - WPWE2 = 1, not closed to bottom, the well stays open");
    // CON on the bottom completion: the well flows over its topmost completion alone -> closed to bottom
    const Out d = run(" 'PROD' 0 0 2 2 1e-6 1* 1* 'CON' 'NO' /\n", "");
    check(d.sum[2] == 1 && d.sum[3] == 1 && d.sum[7] == 0, "C43: WPWE CON on the bottom completion - flowing over the topmost completion alone sets WPWE3");
    // every connection limited (K1 = K2 = 0: all), CON -> both close in one step (WPWE2 = 2), shut; WTEST 'E' every 30 d
    // re-opens them (WPWE1 = 2 at each retest while the limit closes them again)
    const Out e = run(" 'PROD' 0 0 0 0 1e-6 1* 1* 'CON' 'NO' /\n", " 'PROD' 30 'E' /\n");
    const bool reopened = e.log.find("WTEST re-opens 2 connection(s)") != std::string::npos;
    check(e.maxv[2] == 2 && e.sum[7] >= 1 && reopened && e.maxv[1] == 2 && e.tFirst[1] > e.tFirst[2],
          "C43: WPWE + WTEST 'E' - all connections closed (WPWE2 = 2, WPWE7), the retest re-opens them (WPWE1 = 2)");
    // no limits -> no events at all (the initial configuration is not an event)
    const Out n = run(" 'PROD' 0 0 1 1 0.9999 1* 1* 'CON' 'NO' /\n", "");
    int tot = 0; for (int x : n.sum) tot += x;
    check(tot == 0, "C43: WPWE - no events without simulator / schedule changes (initial configuration excluded)");
}
