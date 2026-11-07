// C40 tests: geomodel data analysis - log calculator expressions, regression, variogram fitting, well tie
#include "ressim/geo_analysis.h"
#include "ressim/study.h"
#include "ressim/io.h"

static void testC40GeoAnalysis() {
    namespace fs = std::filesystem;
    using namespace ressim::geo;
    std::printf("== C40 data analysis: log calculator, regression, variogram fit, well tie\n");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    // ---- expressions: precedence, functions, NaN propagation, errors
    {
        const std::vector<double> A = {1, 2, nan, 4}, B = {2, 0, 1, -1};
        const std::map<std::string, const std::vector<double>*> vars = {{"A", &A}, {"b", &B}};
        std::vector<double> o; std::string err;
        auto ev = [&](const char* e) { o.clear(); return evalExpression(e, vars, 4, o, err); };
        auto same = [](const std::vector<double>& x, const std::vector<double>& y) {
            if (x.size() != y.size()) return false;
            for (size_t i = 0; i < x.size(); ++i) if (std::isnan(x[i]) != std::isnan(y[i]) || (!std::isnan(x[i]) && std::fabs(x[i] - y[i]) > 1e-12)) return false;
            return true;
        };
        bool ok = ev("a + B * 2") && same(o, {5, 2, nan, 2});
        ok = ok && ev("if(A > 1, A / B, -1)") && same(o, {-1, nan, nan, -4});          // division by zero -> NaN, NaN condition -> NaN
        ok = ok && ev("clamp(pow(A, 2), 2, 10)") && same(o, {2, 4, nan, 10});
        ok = ok && ev("log10(100) + abs(-2)^2") && same(o, {6, 6, 6, 6});
        ok = ok && ev("-2^2 + 2^3^2") && same(o, {508, 508, 508, 508});                // unary minus below ^, ^ right-associative
        ok = ok && ev("1 < 2 && 3 >= 3 || 0") && same(o, {1, 1, 1, 1});
        ok = ok && ev("isnan(A) + !(B != 0)") && same(o, {0, 1, 1, 0});
        const bool e1 = !ev("A +") && err.find("position") != std::string::npos;
        const bool e2 = !ev("C * 2") && err.find("unknown curve 'C'") != std::string::npos;
        const bool e3 = !ev("clamp(A, 1)") && err.find("3 argument") != std::string::npos;
        std::printf("    expressions %s, errors: '%s'\n", ok ? "OK" : "WRONG", err.c_str());
        check(ok && e1 && e2 && e3, "C40: log-calculator expressions - precedence, right-associative ^, functions, NaN propagation, syntax / name / arity errors");
        // curve added to a log: nulls in, nulls out, case-insensitive names, MD / TVD variables
        WellLog lg; lg.md = {1000, 1001, 1002}; lg.nullValue = -999.25;
        lg.curves.push_back(LogCurve{"GR", "GAPI", "", {30, -999.25, 120}, false});
        GeoWell w; w.name = "T"; w.kb = 0;
        const bool a1 = addCalculatedCurve(lg, &w, "vsh", "V/V", "clamp((gr - 30) / (120 - 30), 0, 1) + 0 * TVD", err);
        const LogCurve* v = lg.curve("VSH");
        const bool a2 = a1 && v && v->v[0] == 0 && v->v[1] == -999.25 && std::fabs(v->v[2] - 1) < 1e-12 && v->descr.rfind("= ", 0) == 0;
        const bool a3 = addCalculatedCurve(lg, &w, "VSH", "V/V", "MD - 1000", err) && lg.curves.size() == 2 && lg.curve("VSH")->v[2] == 2;   // replaces
        check(a2 && a3, "C40: calculated curve - null samples stay null, new curve replaces one with the same name, MD / TVD usable");
    }
    // ---- regression: exact log-linear poro-perm
    {
        std::vector<double> x, y;
        for (int i = 0; i <= 25; ++i) { x.push_back(0.05 + 0.01 * i); y.push_back(std::pow(10.0, -0.6 + 13.0 * x.back())); }
        y.push_back(-5); x.push_back(0.2);   // non-positive y skipped in log mode
        const Regression r = ressim::geo::regress(x, y, true), l = ressim::geo::regress({1, 2, 3, 4}, {3, 5, 7, 9}, false);
        std::printf("    regression: log10 k = %.4f + %.4f phi (R2 %.6f, n %d); linear a %.3f b %.3f\n", r.a, r.b, r.r2, r.n, l.a, l.b);
        check(std::fabs(r.a + 0.6) < 1e-9 && std::fabs(r.b - 13) < 1e-9 && r.r2 > 0.999999 && r.n == 26 && std::fabs(l.a - 1) < 1e-12 && std::fabs(l.b - 2) < 1e-12,
              "C40: cross-plot regression - log-linear poro-perm recovered exactly, non-positive y excluded");
    }
    // ---- variogram fitting: recover known models from noise-free points; then from a sampled 1-D Gaussian field
    {
        bool ok = true; double worst = 0;
        const ressim::geostat::VarioModel models[3] = {ressim::geostat::VarioModel::Spherical, ressim::geostat::VarioModel::Exponential, ressim::geostat::VarioModel::Gaussian};
        for (auto m : models) {
            ressim::geostat::ExpVario ev;
            for (int k = 1; k <= 24; ++k) { ev.lag.push_back(5.0 * k); ev.gamma.push_back(0.1 + 0.9 * structure(m, 5.0 * k / 50.0)); ev.pairs.push_back(100 + 10 * k); }
            const VarioFit f = fitVariogram(ev, m);
            const double e = std::max({std::fabs(f.nugget - 0.1), std::fabs(f.sill - 0.9), std::fabs(f.range - 50) / 50});
            worst = std::max(worst, e);
            ok = ok && f.ok && e < 0.01;
            const VarioFit g = fitVariogram(ev, m, 0.1);   // fixed nugget
            ok = ok && std::fabs(g.nugget - 0.1) < 1e-12 && std::fabs(g.range - 50) < 0.5;
        }
        // field: AR(1) sequence = exponential covariance with practical range 3 / -ln(rho) samples
        std::mt19937 rng(11); std::normal_distribution<double> N(0, 1);
        const double rho = std::exp(-3.0 / 40.0);   // practical range 40
        std::vector<ressim::geostat::Point> pts; double v = 0;
        for (int r = 0; r < 30; ++r) { v = N(rng); for (int i = 0; i < 400; ++i) { v = rho * v + std::sqrt(1 - rho * rho) * N(rng); pts.push_back({1000.0 * r, 0, double(i), v}); } }
        const auto ev = ressim::geostat::experimental(pts, 2.0, 40, 1);
        const VarioFit fe = fitVariogram(ev, ressim::geostat::VarioModel::Exponential);
        std::printf("    variogram fit: worst relative error %.2e on noise-free points; AR(1) field -> exponential nugget %.3f sill %.3f range %.1f (true 0 / 1 / 40)\n",
                    worst, fe.nugget, fe.sill, fe.range);
        check(ok, "C40: variogram fit - nugget, sill and range of spherical / exponential / Gaussian models recovered (< 1 %), fixed nugget honoured");
        check(fe.ok && fe.nugget < 0.1 && std::fabs(fe.nugget + fe.sill - 1) < 0.15 && fe.range > 28 && fe.range < 55,
              "C40: variogram fit on a simulated exponential (AR(1)) field recovers the range within sampling error");
    }
    // ---- well tie on the demo project: sonic + density built consistently with the synthetic SEG-Y
    {
        const std::string dir = (fs::temp_directory_path() / "ressim_c40_tie").string();
        std::string err;
        writeDemoProject(dir, false, err);
        std::string e;
        auto wells = readWellHeads(dir + "/well_heads.txt", e);
        for (auto& w : wells) {
            if (fs::exists(dir + "/" + w.name + "_survey.txt")) readDeviation(dir + "/" + w.name + "_survey.txt", w, e);
            w.logs.push_back(readLas(dir + "/" + w.name + ".las", e));
        }
        for (const auto& t : readTops(dir + "/well_tops.txt", e)) for (auto& w : wells) if (w.name == t.first) w.tops.push_back(t.second);
        const SeismicCube cube = readSegy(dir + "/seismic.segy", e);
        VelocityModel vm; vm.v0 = 1800; vm.k = 0.6; vm.datum = 0;
        double minCorr = 1, maxShift = 0, maxRms = 0; int n = 0;
        for (const auto& w : wells) {
            if (w.logs.empty() || !w.logs[0].curve("DT")) continue;
            TieInput in; in.freq = 25;
            const TieResult r = wellTie(w, w.logs[0], cube, vm, in);
            if (!r.ok) { std::printf("    %s: %s\n", w.name.c_str(), r.msg.c_str()); minCorr = -1; continue; }
            minCorr = std::min(minCorr, r.corr); maxShift = std::max(maxShift, std::fabs(r.shift)); maxRms = std::max(maxRms, r.vmRmsMs); ++n;
        }
        // wrong starting velocity model: the bulk shift absorbs it and the calibrated model predicts the true times
        VelocityModel bad; bad.v0 = 1800; bad.k = 0.65; bad.datum = 0;   // ~30 ms early at the reservoir
        const GeoWell& w1 = wells.empty() ? GeoWell() : wells[0];
        TieResult rb = w1.logs.empty() ? TieResult() : wellTie(w1, w1.logs[0], cube, bad, TieInput());
        double dtTrue = 1e9;
        if (rb.ok && !rb.tdZ.empty()) { const double z = rb.tdZ[rb.tdZ.size() / 2]; dtTrue = std::fabs(depthToTime(rb.vmFit, z) - depthToTime(vm, z)); }
        std::printf("    well tie: %d wells, min correlation %.3f, max |shift| %.1f ms, velocity-fit misfit %.2f ms; wrong v(z): shift %+.1f ms, corr %.3f, calibrated time error %.2f ms\n",
                    n, minCorr, maxShift, maxRms, rb.shift, rb.corr, dtTrue);
        check(n == 6 && minCorr > 0.8 && maxShift < 4 && maxRms < 2, "C40: well tie on the demo - synthetic (sonic x density -> RC * Ricker 25 Hz) correlates > 0.8 with the SEG-Y at every well, |shift| < 4 ms");
        check(rb.ok && std::fabs(rb.shift) > 20 && rb.corr > 0.8 && dtTrue < 4, "C40: well tie with a wrong velocity model - bulk shift found, calibrated v(z) predicts the true TWT (< 4 ms)");
        std::error_code ec; fs::remove_all(dir, ec);
    }
}

// C40: geomodel parameters in Studies - the property model is re-simulated with other seeds / variogram ranges
static void testC40StudyGeoParams() {
    using namespace ressim;
    std::printf("== C40 geomodel parameters in Studies (GeoSeed, GeoRange, property-model persistence)\n");
    Case base = makeDemoCase();
    const GridInput& g0 = base.grid;
    PropModelInput m;
    m.property = 0; m.method = 1; m.model = 1; m.nugget = 0; m.sill = 1; m.rangeMax = m.rangeMin = 1500; m.rangeVert = 30; m.seed = 3;
    m.minV = 0.02; m.maxV = 0.4; m.permFromPoro = true; m.ppA = -0.6; m.ppB = 13; m.kvkh = 0.1; m.source = "test";
    // hard data: a few porosity values at cell centres of the first layer
    const double dx = g0.dx.empty() ? 1000 : g0.dx[0], dy = g0.dy.empty() ? 1000 : g0.dy[0];
    Grid G; G.build(base);
    for (int q : {0, 7, 33, 120, 255}) if (q < int(G.center.size())) m.data.push_back({G.center[size_t(q)].x, G.center[size_t(q)].y, G.center[size_t(q)].z, 0.12 + 0.0005 * q});
    base.propModels = {m};
    StudyParam ps; ps.kind = ParamKind::GeoSeed; ps.target = -1; ps.base = 3; ps.low = 1; ps.high = 100; ps.log = false;
    StudyParam pr; pr.kind = ParamKind::GeoRange; pr.target = 0; pr.base = 1; pr.low = 0.25; pr.high = 4;
    const std::vector<StudyParam> P = {ps, pr};
    const Case a = applyParams(base, P, {7, 1}), b = applyParams(base, P, {7, 1}), c = applyParams(base, P, {8, 1}), u = applyParams(base, P, {3, 1});
    const bool repro = a.grid.poro == b.grid.poro, differs = a.grid.poro != c.grid.poro, untouched = u.grid.poro == base.grid.poro;
    // honoured data, poro-perm transform
    bool hard = true, perm = true;
    for (const auto& d : m.data) {
        double best = 1e300; size_t qb = 0;
        for (size_t q = 0; q < G.center.size(); ++q) { const double e = std::hypot(G.center[q].x - d[0], G.center[q].y - d[1]) + std::fabs(G.center[q].z - d[2]); if (e < best) { best = e; qb = q; } }
        hard = hard && std::fabs(a.grid.poro[qb] - d[3]) < 1e-9;
    }
    for (size_t q = 0; q < a.grid.poro.size(); ++q) perm = perm && std::fabs(a.grid.permx[q] - std::clamp(std::pow(10.0, -0.6 + 13 * a.grid.poro[q]), 1e-4, 1e5)) < 1e-6 * a.grid.permx[q] && std::fabs(a.grid.permz[q] - 0.1 * a.grid.permx[q]) < 1e-9 * a.grid.permx[q];
    // longer range -> smoother field: lag-1 correlation along I, averaged over 6 seeds
    auto lag1 = [&](const Case& cc) {
        const auto& p = cc.grid.poro; const int nx = cc.grid.nx, ny = cc.grid.ny, nz = cc.grid.nz;
        double m0 = 0; for (double v : p) m0 += v; m0 /= double(p.size());
        double num = 0, den = 0;
        for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i + 1 < nx; ++i) { const size_t q = size_t(i + nx * (j + ny * k)); num += (p[q] - m0) * (p[q + 1] - m0); }
        for (double v : p) den += (v - m0) * (v - m0);
        return den > 0 ? num / den : 0.0;
    };
    double rShort = 0, rLong = 0;
    for (int sd = 10; sd < 16; ++sd) { rShort += lag1(applyParams(base, P, {double(sd), 0.25})) / 6; rLong += lag1(applyParams(base, P, {double(sd), 4.0})) / 6; }
    // persistence through the project text (deck with RESSIM lines)
    Case rt = projectFromText(projectText(base));
    const bool persist = rt.propModels.size() == 1 && rt.propModels[0].data.size() == m.data.size() && rt.propModels[0].seed == 3 && rt.propModels[0].permFromPoro &&
                         std::fabs(rt.propModels[0].ppB - 13) < 1e-12 && std::fabs(rt.propModels[0].rangeVert - 30) < 1e-9 && rt.propModels[0].source == "test" &&
                         std::fabs(rt.propModels[0].data[2][3] - m.data[2][3]) < 1e-12;
    std::printf("    seed 7 twice identical %d, seed 8 differs %d, base seed leaves PORO unchanged %d, hard data honoured %d, PERM = 10^(a+b phi) %d; "
                "lag-1 correlation range x0.25 %.3f vs x4 %.3f; project round trip %d (dx %.0f dy %.0f)\n", repro, differs, untouched, hard, perm, rShort, rLong, persist, dx, dy);
    check(repro && differs && untouched && hard && perm, "C40: GeoSeed re-simulates the case's property model - reproducible per seed, conditioned to the hard data, PERM from the poro-perm transform");
    check(rLong > rShort + 0.15, "C40: GeoRange multiplies the variogram ranges - longer ranges give a smoother (more correlated) porosity field");
    check(persist, "C40: property models persist with the project (RESSIM PROPMODEL / PROPDATA lines)");
}
