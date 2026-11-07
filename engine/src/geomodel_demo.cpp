// C32 geomodelling: synthetic but complete demo project (faulted anticline) in industry exchange formats
#include "ressim/geomodel.h"
#include <cmath>
#include <cstdio>
#include <random>
#include <filesystem>

namespace ressim {
namespace geo {

namespace {
// truth model in metres (the files are scaled to feet for a field project)
const double kThrow = 30.0;
double faultX(double y, double z) { return 2600.0 + 0.1 * (y - 1500.0) + 0.5 * (z - 2050.0); }   // fault plane dipping east (hanging wall east)
bool hangingWall(double x, double y, double z) { return x > faultX(y, z); }
double topTrue(double x, double y) {
    double z = 2000.0 + 3e-5 * (x - 1800.0) * (x - 1800.0) + 5e-5 * (y - 1500.0) * (y - 1500.0);
    if (hangingWall(x, y, z)) z += kThrow;
    return z;
}
double midTrue(double x, double y) { return topTrue(x, y) + 35.0 + 5.0 * std::sin(x / 700.0); }
double baseTrue(double x, double y) { return midTrue(x, y) + 40.0; }
double hTrue(int h, double x, double y) { return h == 0 ? topTrue(x, y) : h == 1 ? midTrue(x, y) : baseTrue(x, y); }
double mapError(double x, double y) { return 4.0 * std::sin(x / 900.0) * std::cos(y / 800.0); }   // depth-conversion error of the interpreted maps
double phiTrend(double x, double y) { return 0.20 + 0.04 * std::sin(x / 1100.0 + 0.5) * std::cos(y / 900.0); }
const char* kHz[3] = {"TopRes", "MidRes", "BaseRes"};
}  // namespace

std::vector<std::string> writeDemoProject(const std::string& dir, bool field, std::string& err) {
    std::vector<std::string> files;
    { std::error_code ec; std::filesystem::create_directories(dir, ec); }
    const double L = field ? 3.28084 : 1.0;
    const std::string sep = (!dir.empty() && dir.back() != '/' && dir.back() != '\\') ? "/" : "";
    auto P = [&](const std::string& f) { files.push_back(dir + sep + f); return dir + sep + f; };
    std::mt19937 rng(20261023);
    std::normal_distribution<double> N(0, 1);
    std::uniform_real_distribution<double> U(0, 1);
    // ---- horizons (interpreted maps = truth + depth-conversion error), IRAP classic
    for (int h = 0; h < 3; ++h) {
        Surface s; s.name = kHz[h]; s.nx = 81; s.ny = 61; s.x0 = 0; s.y0 = 0; s.dx = 50 * L; s.dy = 50 * L;
        s.z.resize(size_t(s.nx * s.ny));
        for (int j = 0; j < s.ny; ++j) for (int i = 0; i < s.nx; ++i) { const double x = i * 50.0, y = j * 50.0; s.z[size_t(i + s.nx * j)] = (hTrue(h, x, y) + mapError(x, y)) * L; }
        if (!writeIrap(P(std::string(kHz[h]) + ".irap"), s)) { err = "cannot write " + files.back(); return files; }
    }
    // ---- fault sticks
    {
        FILE* f = std::fopen(P("F1_sticks.txt").c_str(), "w");
        if (!f) { err = "cannot write " + files.back(); return files; }
        std::fprintf(f, "# x y z stick fault  (fault F1, dipping east, throw ~%g)\n", kThrow * L);
        for (int s = 0; s < 5; ++s) {
            const double y = -100.0 + s * 800.0;
            for (double z : {1850.0, 2250.0}) std::fprintf(f, "%.2f %.2f %.2f %d F1\n", faultX(y, z) * L, y * L, z * L, s + 1);
        }
        std::fclose(f);
    }
    // ---- wells: heads, one deviated survey, logs, tops
    struct WD { const char* name; double x, y; bool dev; };
    const WD wd[6] = {{"W1", 1200, 1500, false}, {"W2", 2000, 800, false}, {"W3", 2100, 2300, false}, {"W4", 3200, 1500, false}, {"W5", 1500, 2600, false}, {"W6", 2300, 1400, true}};
    const double kb = 25.0;
    FILE* fh = std::fopen(P("well_heads.txt").c_str(), "w");
    FILE* ft = std::fopen(P("well_tops.txt").c_str(), "w");
    if (!fh || !ft) { err = "cannot write well files"; if (fh) std::fclose(fh); if (ft) std::fclose(ft); return files; }
    std::fprintf(fh, "# name x y kb\n");
    std::fprintf(ft, "# well surface md\n");
    for (const auto& w : wd) {
        GeoWell g; g.name = w.name; g.x = w.x; g.y = w.y; g.kb = kb;
        std::fprintf(fh, "%s %.2f %.2f %.2f\n", w.name, w.x * L, w.y * L, kb * L);
        double td;
        if (w.dev) {
            FILE* fd = std::fopen(P(std::string(w.name) + "_survey.txt").c_str(), "w");
            if (!fd) { err = "cannot write survey"; std::fclose(fh); std::fclose(ft); return files; }
            std::fprintf(fd, "# WELL NAME: %s\n# HEAD X: %.2f\n# HEAD Y: %.2f\n# KB: %.2f\nMD INC AZI\n", w.name, w.x * L, w.y * L, kb * L);
            std::vector<std::array<double, 3>> st;
            st.push_back({0, 0, 300}); st.push_back({1800, 0, 300});
            for (double md = 1830; md <= 2240; md += 30) st.push_back({md, std::min(88.0, (md - 1800) / 30.0 * 6.0), 300});
            st.push_back({2900, 88, 300});
            for (const auto& s : st) {
                std::fprintf(fd, "%.2f %.3f %.3f\n", s[0] * L, s[1], s[2]);
                Station q; q.md = s[0]; q.inc = s[1]; q.azi = s[2]; g.traj.push_back(q);
            }
            std::fclose(fd);
            minimumCurvature(g);
            td = 2900;
        } else td = baseTrue(w.x, w.y) + 60 + kb;
        // logs: Markov sand / shale beds inside the reservoir zones
        WellLog lg; lg.well = w.name;
        LogCurve gr{"GR", "GAPI", "gamma ray", {}, false}, phi{"PHIE", "V/V", "effective porosity", {}, false}, perm{"PERM", "MD", "permeability", {}, false},
            fac{"FACIES", "", "0 shale 1 sand", {}, true};
        // C40: sonic + density for the well tie - background v(z) = 1800 + 0.6 z m/s (the velocity of the seismic model),
        // impedance steps at the three horizons matching the reflection coefficients of the synthetic SEG-Y, small
        // deterministic bed-scale variation (no extra random draws: the other demo outputs stay identical)
        LogCurve dtc{"DT", "US/F", "compressional slowness", {}, false}, rhob{"RHOB", "G/C3", "bulk density", {}, false};
        const double top = topTrue(w.x, w.y) + kb - 100;
        bool sand = false;
        bool topSeen[3] = {false, false, false};
        for (double md = std::max(0.0, std::floor(top)); md <= td; md += 0.5) {
            const P3 p = g.at(md);
            const double z = p[2];
            int zone = -1;
            for (int h = 0; h < 3; ++h)
                if (!topSeen[h] && z >= hTrue(h, p[0], p[1])) { topSeen[h] = true; std::fprintf(ft, "%s %s %.2f\n", w.name, kHz[h], md * L); }
            if (z >= topTrue(p[0], p[1]) && z < midTrue(p[0], p[1])) zone = 0;
            else if (z >= midTrue(p[0], p[1]) && z < baseTrue(p[0], p[1])) zone = 1;
            const double ns = zone == 0 ? 0.62 : 0.45;           // net-to-gross of the zone
            if (zone >= 0) {                                       // bed switching, mean bed ~ 3 m
                const double pSwitch = 0.5 / 3.0;
                if (U(rng) < pSwitch) sand = U(rng) < ns;
            } else sand = false;
            double ph, k, g0;
            if (sand) { ph = std::max(0.08, std::min(0.32, phiTrend(p[0], p[1]) + 0.02 * N(rng))); k = std::pow(10.0, -0.6 + 13.0 * ph + 0.25 * N(rng)); g0 = 38 + 6 * N(rng); }
            else { ph = std::max(0.01, 0.05 + 0.012 * N(rng)); k = std::pow(10.0, -2 + 0.3 * N(rng)); g0 = 108 + 8 * N(rng); }
            lg.md.push_back(md * L);
            gr.v.push_back(g0); phi.v.push_back(ph); perm.v.push_back(k); fac.v.push_back(sand ? 1 : 0);
            {
                double m = 1.0;   // impedance multiplier below the crossed horizons
                const double rcH[3] = {0.08 * (1.0 + 4.0 * (phiTrend(p[0], p[1]) - 0.2)), 0.08 * -0.45, 0.08 * 0.6};
                for (int h = 0; h < 3; ++h) if (topSeen[h]) m *= (1 + rcH[h]) / (1 - rcH[h]);
                const double bed = 1.0 + (sand ? -0.006 : 0.006) + 0.003 * std::sin(md * 1.7);
                const double v = (1800.0 + 0.6 * z) * std::sqrt(m) * bed, rho = 2.30 * std::sqrt(m) / std::sqrt(bed);
                dtc.v.push_back(304800.0 / v); rhob.v.push_back(rho);
            }
        }
        lg.curves = {gr, phi, perm, fac, dtc, rhob};
        if (!writeLas(P(std::string(w.name) + ".las"), lg)) { err = "cannot write " + files.back(); std::fclose(fh); std::fclose(ft); return files; }
    }
    std::fclose(fh); std::fclose(ft);
    // ---- SEG-Y: convolutional synthetic (Ricker 25 Hz), top-reservoir amplitude rises with porosity
    {
        SeismicCube c; c.name = "demo"; c.ni = 41; c.nj = 31; c.dt = 4; c.t0 = 1500; c.ns = 176;
        c.il0 = 1001; c.xl0 = 2001; c.x0 = 0; c.y0 = 0; c.ix = 100 * L; c.iy = 0; c.jx = 0; c.jy = 100 * L;
        c.textHeader = "C 1 ResGeo demo survey - synthetic faulted anticline (convolutional model, Ricker 25 Hz)";
        c.textHeader.resize(3200, ' ');
        c.a.assign(size_t(c.ni) * size_t(c.nj) * size_t(c.ns), 0.0f);
        VelocityModel vm; vm.v0 = 1800; vm.k = 0.6; vm.datum = 0;
        const double f0 = 25.0, pi2 = 9.8696044;
        for (int i = 0; i < c.ni; ++i) for (int j = 0; j < c.nj; ++j) {
            const double x = i * 100.0, y = j * 100.0;
            const double amp[3] = {1.0 + 4.0 * (phiTrend(x, y) - 0.2), -0.45, 0.6};
            for (int h = 0; h < 3; ++h) {
                const double tH = depthToTime(vm, hTrue(h, x, y));
                for (int s = 0; s < c.ns; ++s) {
                    const double tt = (c.t0 + s * c.dt - tH) / 1000.0, a = pi2 * f0 * f0 * tt * tt;
                    c.a[(size_t(i) * size_t(c.nj) + size_t(j)) * size_t(c.ns) + size_t(s)] += float(amp[h] * (1 - 2 * a) * std::exp(-a));
                }
            }
            for (int s = 0; s < c.ns; ++s) c.a[(size_t(i) * size_t(c.nj) + size_t(j)) * size_t(c.ns) + size_t(s)] += float(0.03 * N(rng));
        }
        if (!writeSegy(P("seismic.segy"), c)) { err = "cannot write " + files.back(); return files; }
    }
    // ---- README with the suggested workflow settings
    {
        FILE* f = std::fopen(P("README.txt").c_str(), "w");
        if (f) {
            const char* u = field ? "ft" : "m";
            std::fprintf(f,
                "ResGeo demo project (%s units) - synthetic faulted anticline\n\n"
                "well_heads.txt      name x y kb\n"
                "W1..W6.las          LAS 2.0: GR, PHIE, PERM (mD), FACIES (0 shale / 1 sand), DT (us/ft), RHOB (g/cc), index MD\n"
                "W6_survey.txt       deviation survey MD INC AZI (minimum curvature; horizontal in the lower zone)\n"
                "well_tops.txt       well surface md  (TopRes / MidRes / BaseRes)\n"
                "TopRes/MidRes/BaseRes.irap   IRAP classic depth maps (interpreted: +-4 %s depth-conversion error vs the wells)\n"
                "F1_sticks.txt       fault sticks of F1 (dips east, throw %.0f %s on the east side)\n"
                "seismic.segy        SEG-Y rev1 IEEE, inlines 1001-1041 (x), crosslines 2001-2031 (y), %.0f %s bins, 1500-2200 ms\n\n"
                "Suggested settings: velocity v(z) = %.0f + 0.6 z (%s/s), grid origin (0,0), 40 x 30 cells of %.0f %s,\n"
                "zones TopRes-MidRes (8 layers), MidRes-BaseRes (8 layers); the top-reservoir amplitude correlates with porosity.\n",
                field ? "FIELD" : "METRIC", u, kThrow * L, u, 100 * L, u, 1800 * L, u, 100 * L, u);
            std::fclose(f);
        }
    }
    return files;
}


std::string runDemoWorkflow(const std::string& dir, bool field, std::string& err) {
    err.clear();
    if (writeDemoProject(dir, field, err).empty() || !err.empty()) return std::string();
    const double L = field ? 3.28084 : 1.0;
    std::string e;
    std::vector<Surface> hz;
    for (const char* h : {"TopRes", "MidRes", "BaseRes"}) { Surface s = readIrap(dir + "/" + h + ".irap", e); s.name = h; hz.push_back(s); }
    const auto faults = readFaultSticks(dir + "/F1_sticks.txt", e);
    auto wells = readWellHeads(dir + "/well_heads.txt", e);
    for (auto& w : wells) {
        if (std::filesystem::exists(dir + "/" + w.name + "_survey.txt")) readDeviation(dir + "/" + w.name + "_survey.txt", w, e);
        w.logs.push_back(readLas(dir + "/" + w.name + ".las", e));
    }
    for (const auto& t : readTops(dir + "/well_tops.txt", e)) for (auto& w : wells) if (w.name == t.first) w.tops.push_back(t.second);
    if (!e.empty()) { err = e; return std::string(); }
    GridSpec gs; gs.dx = gs.dy = 100 * L; gs.nx = 40; gs.ny = 30;
    Zone z1; z1.name = "Upper"; z1.top = 0; z1.base = 1; z1.nLayers = 8;
    Zone z2; z2.name = "Lower"; z2.top = 1; z2.base = 2; z2.nLayers = 8;
    gs.zones = {z1, z2}; gs.tieRange = 1500 * L;
    const GeoGrid g = buildGrid(gs, hz, faults, wells, err);
    if (!err.empty()) return std::string();
    const auto B = blockWells(g, wells);
    FaciesSpec fs; fs.method = 0; fs.vg.rangeMax = fs.vg.rangeMin = 800 * L; fs.vg.rangeVert = 6 * L; fs.seed = 3;
    const FaciesResult F = modelFacies(g, B, fs);
    PetroSpec ps; ps.vg.rangeMax = ps.vg.rangeMin = 1000 * L; ps.vg.rangeVert = 8 * L; ps.maxV = 0.4; ps.seed = 11;
    const PetroResult phi = modelProperty(g, B, &F.facies, ps);
    const auto fits = fitPoroPerm(B, "PHIE", "PERM", &F.facies, F.codes);
    const auto kx = permFromPoro(phi.v, fits, &F.facies, F.codes, 0.2, 5);
    std::vector<double> kz(kx.size()), ntg(kx.size(), 1.0);
    for (size_t c = 0; c < kx.size(); ++c) kz[c] = 0.1 * kx[c];
    VolSpec vs; vs.contacts.owc = 2080 * L; vs.satHeight = false; vs.swConst = 0.25; vs.bo = 1.25; vs.field = field;
    const VolResult V = volumetrics(g, phi.v, ntg, &F.facies, vs);
    std::vector<double> real;
    for (unsigned r = 1; r <= 5; ++r) { PetroSpec q = ps; q.seed = 100 + r; real.push_back(volumetrics(g, modelProperty(g, B, &F.facies, q).v, ntg, nullptr, vs).total.stoiip); }
    const Percentiles P = percentiles(real);
    UpscaleSpec us; us.fx = 2; us.fy = 2; us.fz = 3; us.permMethod = 3;
    const Upscaled U = upscale(g, phi.v, ntg, kx, kx, kz, us);
    const bool w1 = writeGrdecl(dir + "/GEOMODEL_FINE.GRDECL", g, {{"PORO", &phi.v}, {"NTG", &ntg}, {"PERMX", &kx}, {"PERMY", &kx}, {"PERMZ", &kz}});
    const bool w2 = writeGrdecl(dir + "/GEOMODEL_COARSE.GRDECL", U.g, {{"PORO", &U.poro}, {"NTG", &U.ntg}, {"PERMX", &U.kx}, {"PERMY", &U.ky}, {"PERMZ", &U.kz}});
    if (!w1 || !w2) { err = "cannot write the GRDECL files in " + dir; return std::string(); }
    int act = 0; for (int a : g.grid.actnum) act += a != 0;
    char b[512];
    std::string rep;
    std::snprintf(b, sizeof b, "ResGeo demo workflow (%s units) in %s\n", field ? "FIELD" : "METRIC", dir.c_str()); rep += b;
    std::snprintf(b, sizeof b, "  data      : 3 horizons, %zu fault(s), %zu wells\n", faults.size(), wells.size()); rep += b;
    std::snprintf(b, sizeof b, "  grid      : %d x %d x %d, %d active, %d fault blocks\n", g.grid.nx, g.grid.ny, g.grid.nz, act, g.nBlocks); rep += b;
    for (const auto& t : g.tieReport) rep += "  well tie  : " + t + "\n";
    std::snprintf(b, sizeof b, "  blocking  : %zu cells\n", B.size()); rep += b;
    for (size_t k = 0; k < F.codes.size(); ++k) { std::snprintf(b, sizeof b, "  facies    : code %d  %.1f %%\n", F.codes[k], 100 * F.proportions[k]); rep += b; }
    std::snprintf(b, sizeof b, "  porosity  : mean %.4f sd %.4f (%d data)\n", phi.mean, phi.sd, phi.data); rep += b;
    for (size_t q = 0; q < fits.size(); ++q) { std::snprintf(b, sizeof b, "  poro-perm : facies %d  log10 k = %.3f + %.3f phi (r2 %.2f)\n", q < F.codes.size() ? F.codes[q] : -1, fits[q].a, fits[q].b, fits[q].r2); rep += b; }
    std::snprintf(b, sizeof b, "  STOIIP    : %.5g %s  (5 realisations: P90 %.5g  P50 %.5g  P10 %.5g)\n", V.total.stoiip, field ? "stb" : "sm3", P.p90, P.p50, P.p10); rep += b;
    std::snprintf(b, sizeof b, "  upscaled  : %d x %d x %d (flow-based permeability)\n", U.g.grid.nx, U.g.grid.ny, U.g.grid.nz); rep += b;
    rep += "  written   : GEOMODEL_FINE.GRDECL, GEOMODEL_COARSE.GRDECL (INCLUDE them in the GRID section)\n";
    return rep;
}

}  // namespace geo
}  // namespace ressim
