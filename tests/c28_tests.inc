// =====================================================================================================================
// C28 tests (comparison gaps + completion round): FOSM prediction uncertainty + data worth, partitioning tracers and
// decay chains, geostatistics (kriging / SGS) and Gaussian-process proxies, decline-curve / nodal analysis, fault slip in
// the thermal and compositional engines, MPI study responses, research cycle 7 drivers
// =====================================================================================================================

static void testC28Fosm() {
    std::printf("== C28 FOSM prediction uncertainty + data worth (iTOUGH2)\n");
    Case base = smallCase();
    std::vector<StudyParam> P(1);
    P[0].kind = ParamKind::PermXY; P[0].base = 1; P[0].low = 0.3; P[0].high = 3; P[0].log = true;
    const std::vector<double> truthV = {1.5};
    const Case truthCase = applyParams(base, P, truthV);
    const RunResult truth = runSingle(truthCase);
    base.observations = syntheticObservations(truth, truthCase, {ObsQty::WBHP}, 0.01, 7);
    // keep only the producers' WBHP so that candidate data (WWCT, FPR ...) exist
    StudyCallbacks cb;
    HMStatsOptions o; o.predict = {Response::FOPT, Response::FWPT};
    const HMStatistics S = hmStatistics(base, P, truthV, o, cb);
    if (!S.ok) { std::printf("    error: %s\n", S.error.c_str()); check(false, "FOSM statistics"); return; }
    // independent derivative of FOPT with respect to log10 k (central difference, different step)
    auto fopt = [&](double v) { const Case c = applyParams(base, P, {v}); return responseValue(runSingle(c), c, Response::FOPT, o.eco, {}, {}); };
    const double d = 0.01, x0 = std::log10(truthV[0]);
    const double g = (fopt(std::pow(10.0, x0 + d)) - fopt(std::pow(10.0, x0 - d))) / (2 * d);
    const double ref = std::fabs(g) * S.stdDev[0];
    std::printf("    FOPT %.6g +- %.4g (FOSM)   |dFOPT/dlog10k| sigma_x = %.4g   FWPT %.6g +- %.4g\n", S.pred[0].value, S.pred[0].stdDev, ref, S.pred[1].value, S.pred[1].stdDev);
    check(S.pred.size() == 2 && std::isfinite(S.pred[0].stdDev) && S.pred[0].stdDev > 0, "FOSM prediction standard deviation finite and positive");
    check(std::fabs(S.pred[0].stdDev - ref) < 0.05 * ref, "one parameter: sigma_z = |dz/dx| sigma_x (independent central difference, 5 %)");
    check(std::fabs(S.pred[0].hi - S.pred[0].lo - 2 * 1.96 * S.pred[0].stdDev) < 1e-9 * std::max(1.0, S.pred[0].value), "95 % prediction interval");
    bool monoE = true, monoC = true; int nE = 0, nC = 0; bool removeAllUnbounded = false;
    for (const auto& w : S.worth) {
        std::printf("    %-14s %-9s %3d pts  dsd FOPT %+.3g %%  FWPT %+.3g %%\n", w.set.c_str(), w.existing ? "observed" : "candidate", w.points, w.dStdPct[0], w.dStdPct[1]);
        if (w.existing) { ++nE; monoE = monoE && w.dStdPct[0] >= -1e-9; if (!std::isfinite(w.dStdPct[0])) removeAllUnbounded = true; }
        else { ++nC; monoC = monoC && w.dStdPct[0] <= 1e-9; }
    }
    check(nE >= 1 && nC >= 1, "data worth lists observed and candidate data sets");
    check(monoE && monoC, "removing data never reduces, adding data never increases the prediction uncertainty (Loewner order)");
    // with a single observed data set, removing it leaves the parameter undetermined
    check(nE != 1 || removeAllUnbounded, "removing the only data set makes the prediction unbounded");
}

// 1-D water-swept core at residual oil (S_or = 0.2 immobile): PITT, dispersion, decay-chain tests
static Case c28CoreCase(const std::vector<Case::TracerProp>& props, const std::vector<std::string>& tracers, bool inject, double days = 900, int steps = 90) {
    Case c = makeDemoCase();
    GeomInput gi; gi.nx = 60; gi.ny = 1; gi.nz = 1; gi.dx = 10; gi.dy = 10; gi.dz = {10};
    generateGeometry(gi, c.grid);
    const size_t n = size_t(c.grid.size());
    c.grid.poro.assign(n, 0.2); c.grid.permx.assign(n, 500); c.grid.permy.clear(); c.grid.permz.assign(n, 50); c.grid.ntg.clear();
    c.init.pressure.assign(n, 4800); c.init.swat.assign(n, 0.8); c.init.rs.assign(n, 1.0);
    for (auto& w : c.wells) {
        if (w.type == WellType::Producer) { w.i = 59; w.j = 0; } else { w.i = 0; w.j = 0; }
        for (auto& cp : w.completions) { cp.i = w.i; cp.j = w.j; cp.k1 = 0; cp.k2 = 0; }
    }
    ReportStep st = c.schedule.front();
    for (auto& wc : st.controls) {
        if (wc.well == "PROD") { wc.mode = Control::BHP; wc.bhp = 4000; wc.orat = 1e6; wc.status = inject ? WellStatus::Open : WellStatus::Shut; }
        if (wc.well == "INJ") { wc.mode = Control::RATE; wc.rate = 7; wc.bhp = 1e5; wc.status = inject ? WellStatus::Open : WellStatus::Shut; }
    }
    st.length = days / steps;
    c.schedule.assign(size_t(steps), st);
    for (size_t k = 1; k < c.schedule.size(); ++k) c.schedule[k].controls.clear();
    ExtraKeyword t; t.name = "TRACER"; t.section = DeckSection::Props;
    for (const auto& nm : tracers) t.body += " '" + nm + "' 'WAT' /\n";
    t.body += "/\n";
    c.extraKeywords.push_back(t);
    if (inject) {
        ExtraKeyword w; w.name = "WTRACER"; w.section = DeckSection::Schedule; w.step = 0;
        for (const auto& nm : tracers) w.body += " 'INJ' '" + nm + "' 1.0 /\n";
        w.body += "/\n";
        c.extraKeywords.push_back(w);
    }
    c.tracerProps = props;
    return c;
}

static void testC28Tracers() {
    std::printf("== C28 tracers: partitioning (PITT), dispersion / diffusion, decay chains (Bateman)\n");
    struct Run { std::vector<double> t; std::vector<std::vector<double>> cProd; std::vector<double> inPlace0, inPlace, prodCum, injCum, dec; double soAvg = 0; };
    auto run = [](const Case& c) {
        Run R;
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        const size_t nt = s.tracers().size();
        for (size_t k = 0; k < nt; ++k) R.inPlace0.push_back(s.tracerInPlace(k));
        R.cProd.assign(nt, {});
        while (s.runReportStep()) {
            const auto& P = s.summary().back();
            R.t.push_back(P.time);
            int pw = -1; for (size_t w = 0; w < c.wells.size(); ++w) if (c.wells[w].type == WellType::Producer) pw = int(w);
            for (size_t k = 0; k < nt; ++k) {
                const double qw = pw >= 0 && size_t(pw) < P.wells.size() ? P.wells[size_t(pw)].water : 0.0;
                const double rt = pw >= 0 && size_t(pw) < P.wells.size() && k < P.wells[size_t(pw)].comp.size() ? P.wells[size_t(pw)].comp[k] * 86400.0 : 0.0;
                R.cProd[k].push_back(qw > 0 ? rt / qw : 0.0);
            }
        }
        const auto& E = s.summary().back();
        for (size_t k = 0; k < nt; ++k) {
            R.inPlace.push_back(s.tracerInPlace(k)); R.dec.push_back(s.tracerDecayed(k));
            R.prodCum.push_back(k < E.fcompPT.size() ? E.fcompPT[k] : 0.0); R.injCum.push_back(k < E.fcompIT.size() ? E.fcompIT[k] : 0.0);
        }
        return R;
    };
    auto tAt = [](const Run& R, size_t k, double level) {   // first time the produced concentration reaches 'level' (linear interpolation)
        for (size_t i = 1; i < R.t.size(); ++i)
            if (R.cProd[k][i] >= level && R.cProd[k][i - 1] < level)
                return R.t[i - 1] + (level - R.cProd[k][i - 1]) / (R.cProd[k][i] - R.cProd[k][i - 1]) * (R.t[i] - R.t[i - 1]);
        return -1.0;
    };
    // 1) partitioning inter-well tracer test: retardation R = 1 + K S_or / S_w = 1 + 2 * 0.2 / 0.8 = 1.5
    {
        Case::TracerProp P; P.name = "TP"; P.kPart = 2.0;
        const Run R = run(c28CoreCase({P}, {"TW", "TP"}, true, 900, 180));
        const double tw = tAt(R, 0, 0.5), tp = tAt(R, 1, 0.5);
        const double ratio = tp / tw, Sor = 0.2, Kp = 2.0, sorEst = (ratio - 1) / (ratio - 1 + Kp);
        double bal[2];
        for (int k = 0; k < 2; ++k) bal[k] = (R.inPlace[size_t(k)] + R.prodCum[size_t(k)] + R.dec[size_t(k)] - R.injCum[size_t(k)] - R.inPlace0[size_t(k)]) / std::max(1e-30, R.injCum[size_t(k)]);
        std::printf("    t50: conservative %.4g d, partitioning (K = 2) %.4g d, ratio %.4f (theory 1.5) -> S_or estimate %.4f (true 0.2); balance %.1e / %.1e\n",
                    tw, tp, ratio, sorEst, bal[0], bal[1]);
        check(tw > 0 && tp > 0 && std::fabs(ratio - 1.5) < 0.03, "PITT: partitioning tracer retarded by 1 + K S_or / S_w (3 %)");
        check(std::fabs(sorEst - Sor) < 0.01, "residual oil saturation recovered from the tracer separation (Cooke 1971) within 0.01");
        check(std::fabs(bal[0]) < 1e-8 && std::fabs(bal[1]) < 1e-8, "mass balance of the conservative and the partitioning tracer (1e-8)");
    }
    // 2) dispersion: alphaL = 30 ft spreads the front - variance increase 2 alpha t50^2 / L (advection-dispersion, 30 %)
    {
        Case::TracerProp P; P.name = "TD"; P.alphaL = 30;
        const Run R = run(c28CoreCase({P}, {"TW", "TD"}, true, 900, 180));
        const double L = 590.0;   // injector to producer centre distance (ft)
        auto sig = [&](size_t k) { const double a = tAt(R, k, 0.1587), b = tAt(R, k, 0.8413); return 0.5 * (b - a); };
        const double s0 = sig(0), s1 = sig(1), t50 = tAt(R, 0, 0.5);
        const double dvar = s1 * s1 - s0 * s0, theory = 2 * 30.0 * t50 * t50 / L;
        const double bal = (R.inPlace[1] + R.prodCum[1] + R.dec[1] - R.injCum[1] - R.inPlace0[1]) / R.injCum[1];
        std::printf("    front spread sigma_t: %.4g d (numerical only) -> %.4g d with alphaL 30 ft; variance increase %.4g d^2 vs 2 alpha t50^2 / L = %.4g; balance %.1e\n",
                    s0, s1, dvar, theory, bal);
        check(s1 > s0 && std::fabs(dvar / theory - 1) < 0.3, "longitudinal dispersion widens the breakthrough by 2 alphaL t^2 / L (30 %)");
        check(std::fabs(bal) < 1e-8, "dispersion keeps the mass balance (1e-8)");
    }
    // 3) molecular diffusion in a closed core: left half at c = 1 equilibrates to 0.5, mass conserved
    {
        Case::TracerProp P; P.name = "TM"; P.dm = 50.0;   // ft2/day (large: equilibrates within the run)
        Case c = c28CoreCase({P}, {"TM"}, false, 3000, 30);
        ExtraKeyword k; k.name = "TBLKFTM"; k.section = DeckSection::Solution;
        k.body = " 30*1.0 30*0.0 /\n";
        c.extraKeywords.push_back(k);
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        const double m0 = s.tracerInPlace(0);
        while (s.runReportStep()) {}
        double mn = 1e30, mx = -1e30; for (double v : s.tracerConc()) { mn = std::min(mn, v); mx = std::max(mx, v); }
        std::printf("    diffusion (closed core, D = 50 ft2/d, 3000 d): c in [%.5f, %.5f], mass %.10g -> %.10g\n", mn, mx, m0, s.tracerInPlace(0));
        check(m0 > 0 && std::fabs(s.tracerInPlace(0) / m0 - 1) < 1e-9, "diffusion conserves the tracer (closed system)");
        check(std::fabs(mn - 0.5) < 0.02 && std::fabs(mx - 0.5) < 0.02, "diffusion equilibrates the two halves to c = 0.5");
    }
    // 4) decay chain parent -> daughter without flow: exact Bateman solution for any step
    {
        Case::TracerProp A; A.name = "PA"; A.halfLife = 100; A.init = 1.0;
        Case::TracerProp B; B.name = "DA"; B.halfLife = 300; B.parent = "PA";
        Case::TracerProp C; C.name = "ST"; C.parent = "DA";                    // stable grand-daughter
        const Run R = run(c28CoreCase({A, B, C}, {"PA", "DA", "ST"}, false, 600, 6));
        const double l1 = std::log(2.0) / 100, l2 = std::log(2.0) / 300, t = R.t.back(), N0 = R.inPlace0[0];
        const double n1 = N0 * std::exp(-l1 * t), n2 = N0 * l1 / (l2 - l1) * (std::exp(-l1 * t) - std::exp(-l2 * t));
        const double n3 = N0 - n1 - n2;
        std::printf("    chain PA (100 d) -> DA (300 d) -> ST (stable), 100-day steps: N1 %.10f (%.10f) N2 %.10f (%.10f) N3 %.10f (%.10f) of N0\n",
                    R.inPlace[0] / N0, n1 / N0, R.inPlace[1] / N0, n2 / N0, R.inPlace[2] / N0, n3 / N0);
        check(std::fabs(R.inPlace[0] / n1 - 1) < 1e-9 && std::fabs(R.inPlace[1] / n2 - 1) < 1e-9, "Bateman parent / daughter amounts exact for 1-half-life steps (1e-9)");
        check(std::fabs((R.inPlace[0] + R.inPlace[1] + R.inPlace[2]) / N0 - 1) < 1e-9, "chain with a stable end member conserves the total amount");
    }
    // 5) deck round trip of the new TRACERPROP items
    {
        Case::TracerProp P; P.name = "TP"; P.kPart = 2.5; P.alphaL = 12; P.dm = 0.3; P.halfLife = 50;
        Case::TracerProp Q; Q.name = "TQ"; Q.parent = "TP";
        const Case c = c28CoreCase({P, Q}, {"TP", "TQ"}, true, 300, 10);
        const std::string path = "/tmp/c28_tracer.DATA";
        exportEclipseDeck(c, path);
        const Case d = importEclipseDeck(path);
        bool ok = d.tracerProps.size() == 2;
        if (ok) { const auto& a = d.tracerProps[0]; const auto& b = d.tracerProps[1];
                  ok = a.kPart == 2.5 && a.alphaL == 12 && a.dm == 0.3 && a.halfLife == 50 && a.parent.empty() && b.parent == "TP" && b.kPart == 0; }
        check(ok, "RESSIM TRACERPROP K / parent / alphaL / Dm deck round trip");
    }
}

#include "ressim/geostat.h"
#include "ressim/prodanalysis.h"
#include "ressim/coupled_deck.h"

static void testC28Geostat() {
    std::printf("== C28 geostatistics: kriging exactness, SGS variogram / histogram reproduction, case application\n");
    using namespace geostat;
    Options o; o.vg.model = VarioModel::Spherical; o.vg.nugget = 0; o.vg.sill = 1; o.vg.rangeMax = o.vg.rangeMin = 10; o.vg.rangeVert = 1;
    o.maxNbr = 24; o.seed = 11;
    std::vector<Xyz> T; const int N = 80;
    for (int j = 0; j < N; ++j) for (int i = 0; i < N; ++i) T.push_back({double(i), double(j), 0.0});
    // kriging: exact at the data, variance 0 there, sill far away; ordinary-kriging weights sum to 1 (mean reproduced)
    std::vector<Point> D = {{10, 10, 0, 1.0}, {30, 12, 0, 3.0}, {55, 60, 0, 2.0}, {70, 20, 0, 5.0}, {20, 70, 0, 4.0}};
    const KrigeOut K = krige(T, D, o);
    double errData = 0, varData = 0;
    for (const auto& p : D) { const size_t q = size_t(p.y) * N + size_t(p.x); errData = std::max(errData, std::fabs(K.est[q] - p.v)); varData = std::max(varData, K.var[q]); }
    const size_t far = size_t(79) * N + 79; double mean = 0; for (const auto& p : D) mean += p.v; mean /= double(D.size());
    std::printf("    kriging: max error at data %.2e, max variance at data %.2e, far-field estimate %.4f (data mean %.4f), variance %.4f\n",
                errData, varData, K.est[far], mean, K.var[far]);
    check(errData < 1e-9 && varData < 1e-9, "kriging is an exact interpolator with zero variance at the data");
    check(K.var[far] > 0.99 && std::fabs(K.est[far] - mean) < 1e-6, "beyond the search: the ordinary-kriging mean and the sill variance");
    // unconditional SGS (Gaussian, normalScore off): experimental variogram reproduces the model
    Options u = o; u.normalScore = false; u.skMean = 0;
    const std::vector<double> z = sgs(T, {}, u);
    std::vector<Point> P; for (size_t q = 0; q < T.size(); ++q) P.push_back({T[q][0], T[q][1], 0, z[q]});
    const ExpVario E = experimental(P, 2.0, 10, 0, 0, 3000000);
    double zm = 0, zv = 0; for (double v : z) zm += v; zm /= double(z.size()); for (double v : z) zv += (v - zm) * (v - zm); zv /= double(z.size());
    double relErr = 0; int nl = 0; std::string line;
    for (size_t k = 0; k < E.lag.size(); ++k) {
        const double g = gamma(o.vg, E.lag[k], 0, 0);
        char b[48]; std::snprintf(b, sizeof b, " %.0f:%.2f/%.2f", E.lag[k], E.gamma[k], g); line += b;
        if (E.lag[k] <= 10) { relErr += std::fabs(E.gamma[k] / zv - g) / g; ++nl; }   // normalised by the realisation variance (ergodic fluctuation)
    }
    relErr /= std::max(1, nl);
    std::printf("    SGS 80x80 (range 10): mean %.3f, variance %.3f, lag:gamma/model%s; mean rel. error up to the range %.3f\n", zm, zv, line.c_str(), relErr);
    check(std::fabs(zm) < 0.25 && std::fabs(zv - 1) < 0.25, "SGS realisation: mean 0 / variance 1 within the ergodic fluctuation");
    check(relErr < 0.15, "SGS reproduces the variogram model up to the range (15 %)");
    // conditional SGS with normal score: data honoured exactly, values inside the data range
    const std::vector<double> zc = sgs(T, D, o);
    double errC = 0, mn = 1e30, mx = -1e30;
    for (const auto& p : D) errC = std::max(errC, std::fabs(zc[size_t(p.y) * N + size_t(p.x)] - p.v));
    for (double v : zc) { mn = std::min(mn, v); mx = std::max(mx, v); }
    std::printf("    conditional SGS: max data error %.2e, values in [%.3f, %.3f] (data 1 .. 5)\n", errC, mn, mx);
    check(errC < 1e-12 && mn >= 1.0 - 1e-12 && mx <= 5.0 + 1e-12, "conditional SGS honours the data and the normal-score back transform stays in the data range");
    // application to a case: porosity from the well cells; the well cells keep their values
    Case c = smallCase();
    for (size_t q = 0; q < c.grid.poro.size(); ++q) c.grid.poro[q] = 0.1 + 0.2 * double(q % 7) / 6.0;
    std::vector<double> before;
    for (const auto& w : c.wells) for (const auto& cp : w.completions) for (int k = cp.k1; k <= cp.k2; ++k) before.push_back(c.grid.poro[size_t(cp.i + c.grid.nx * (cp.j + c.grid.ny * k))]);
    ApplyInput ai; ai.property = 0; ai.method = 1; ai.opt.vg.rangeMax = ai.opt.vg.rangeMin = 1500; ai.opt.vg.rangeVert = 60; ai.opt.seed = 3;
    const ApplyReport R = applyToCase(c, ai);
    std::vector<double> after;
    for (const auto& w : c.wells) for (const auto& cp : w.completions) for (int k = cp.k1; k <= cp.k2; ++k) after.push_back(c.grid.poro[size_t(cp.i + c.grid.nx * (cp.j + c.grid.ny * k))]);
    std::printf("    case: %s\n", R.message.c_str());
    check(R.ok && R.cells == 200 && R.dataPoints == int(before.size()) && before == after, "SGS porosity on the case honours the well cells");
    ApplyInput ak = ai; ak.property = 1; ak.method = 0;
    const ApplyReport R2 = applyToCase(c, ak);
    check(R2.ok && c.grid.permz.size() == c.grid.permx.size() && std::fabs(c.grid.permz[5] / c.grid.permx[5] - 0.1) < 1e-12, "kriged log-permeability with PERMZ = kv/kh PERMX");
}

static void testC28GpProxy() {
    std::printf("== C28 Gaussian-process (kriging) proxy vs. quadratic response surface\n");
    std::vector<StudyParam> P(3);
    for (auto& p : P) { p.kind = ParamKind::PermXY; p.base = 1; p.low = 0; p.high = 1; p.log = false; p.active = true; }
    auto f = [](const std::vector<double>& v) { return std::sin(4 * v[0]) + 2 * v[1] * v[1] * v[0] + 0.3 * v[2]; };
    const auto X = designPoints(P, DesignKind::LatinHypercube, 40, 5);
    std::vector<double> y; for (const auto& v : X) y.push_back(f(v));
    const ProxyModel gp = fitProxyGP(P, X, y), pq = fitProxy(P, X, y, true);
    std::mt19937 rng(9); std::uniform_real_distribution<double> U(0, 1);
    double eg = 0, ep = 0, cover = 0; const int nt = 300;
    for (int i = 0; i < nt; ++i) {
        std::vector<double> v = {U(rng), U(rng), U(rng)};
        const double t = f(v), g = proxyPredict(gp, P, v), q = proxyPredict(pq, P, v), sd = proxyStd(gp, P, v);
        eg += (g - t) * (g - t); ep += (q - t) * (q - t);
        if (std::fabs(g - t) <= 1.96 * sd + 1e-9) cover += 1;
    }
    eg = std::sqrt(eg / nt); ep = std::sqrt(ep / nt); cover /= nt;
    const double sdTrain = proxyStd(gp, P, X[0]);
    std::printf("    GP: R2 %.5f Q2 %.4f, test RMSE %.4g (95 %% band covers %.0f %%), sd at a training point %.2e | quadratic: R2 %.4f Q2 %.4f, test RMSE %.4g\n",
                gp.r2, gp.q2, eg, 100 * cover, sdTrain, pq.r2, pq.q2, ep);
    check(gp.ok && gp.kind == 1 && gp.samples == 40 && gp.gpEll.size() == 3, "GP proxy fitted (ARD length scales by maximum likelihood)");
    check(gp.q2 > 0.97 && gp.q2 > pq.q2, "GP leave-one-out Q2 > 0.97 and better than the quadratic surface");
    check(eg < 0.6 * ep, "GP test error well below the quadratic surface (non-polynomial response)");
    check(cover > 0.85 && sdTrain < 0.05, "GP predictive band: >= 85 % coverage, near-zero sd at the training points");
}

static void testC28ProdAnalysis() {
    std::printf("== C28 production analysis: Arps decline, IPR / VLP nodal analysis, VFPPROD from the tubing model\n");
    using namespace prod;
    // Arps hyperbolic with 2 % noise, starting after a 60-day build-up
    std::mt19937 rng(4); std::normal_distribution<double> N(0, 0.02);
    std::vector<double> t, q;
    for (int i = 0; i <= 120; ++i) {
        const double ti = 10.0 * i;
        const double qi = ti < 60 ? 300 + 10 * ti : 1000.0 / std::pow(1 + 0.6 * 0.002 * (ti - 60), 1 / 0.6);
        t.push_back(ti); q.push_back(qi * (1 + (ti < 60 ? 0 : N(rng))));
    }
    const ArpsFit F = fitArps(t, q);
    double tEnd = 0; const double eur = arpsEUR(F, 50, 1e5, &tEnd);
    const ArpsFit Ft{true, "", 1000, 0.002, 0.6, 60, 1, 0, 0};
    double tEndT = 0; const double eurT = arpsEUR(Ft, 50, 1e5, &tEndT);
    std::printf("    hyperbolic: qi %.1f (1000) Di %.5f (0.002) b %.3f (0.6) R2 %.4f, %d points from t0 = %.0f d; EUR to 50/d %.4g (true %.4g) at %.0f d\n",
                F.qi, F.Di, F.b, F.r2, F.points, F.t0, eur, eurT, tEnd);
    check(F.ok && std::fabs(F.b - 0.6) < 0.12 && std::fabs(F.Di / 0.002 - 1) < 0.15 && std::fabs(F.qi / 1000 - 1) < 0.05 && F.t0 == 60, "Arps fit recovers qi, Di, b from noisy data (fit from the peak)");
    check(std::fabs(eur / eurT - 1) < 0.08, "EUR to the economic limit within 8 %");
    // closed-form checks: cumulative = integral of the rate (trapezoid on a fine grid)
    double cum = 0; for (int i = 0; i < 20000; ++i) { const double a = 60 + 0.05 * i, b = a + 0.05; cum += 0.025 * (arpsRate(Ft, a) + arpsRate(Ft, b)); }
    check(std::fabs(cum / arpsCum(Ft, 1060) - 1) < 1e-6, "Arps cumulative = integral of the rate");
    std::vector<double> qe; for (double ti : t) qe.push_back(800 * std::exp(-0.001 * ti));
    const ArpsFit Fe = fitArps(t, qe);
    check(Fe.ok && Fe.b < 0.03 && std::fabs(Fe.Di / 0.001 - 1) < 0.02, "exponential decline identified (b ~ 0)");
    // IPR: continuity at the bubble point, AOF of the composite Vogel
    const Ipr I{3000, 2000, 2.0};
    const double qb = iprRate(I, 2000), qb2 = iprRate(I, 1999.999), aof = iprAof(I);
    check(std::fabs(qb - 2000) < 1e-9 && std::fabs(qb2 - qb) < 0.01 && std::fabs(aof - (2000 + 2.0 * 2000 / 1.8)) < 1e-9, "composite Vogel IPR: J above pb, continuous, AOF = qb + J pb / 1.8");
    // VLP: water column + Darcy-Weisbach friction (independent Colebrook calculation in SI units)
    const Case c = makeDemoCase();
    Tubing tb; tb.depthTVD = tb.lengthMD = 6561.68; tb.diameter = 0.25; tb.roughness = 5e-5;
    const double thp = 200;
    const double qw = 1999.0, wct = qw / (qw + 1.0);
    const double bhp = vlpBhp(c, tb, 1.0, thp, wct, 0.0);
    {
        const double Bw = 1.029 * (1 - 3.13e-6 * ((thp + bhp) / 2 - 4014.7)), rho = 64.79 * 16.01846 / Bw, mu = 0.31e-3;
        const double D = 0.25 * 0.3048, A = 0.25 * 3.141592653589793 * D * D, v = qw * 0.158987 * Bw / 86400 / A, Re = rho * v * D / mu;
        double fD = 0.02; for (int it = 0; it < 50; ++it) fD = std::pow(-2 * std::log10(5e-5 * 0.3048 / D / 3.7 + 2.51 / (Re * std::sqrt(fD))), -2);
        const double dp = (rho * 9.80665 + fD * rho * v * v / (2 * D)) * 2000.0 / 6894.757;
        std::printf("    VLP water 1999 stb/d in 3-in tubing, 2000 m: BHP %.2f psi vs Colebrook %.2f psi (Re %.3g, f %.4f)\n", bhp, thp + dp, Re, fD);
        check(std::fabs(bhp - (thp + dp)) < 0.02 * dp, "VLP = hydrostatic + Darcy-Weisbach friction (Haaland vs Colebrook, 2 %)");
    }
    // operating point and generated VFP table
    const Ipr I2{4500, 0, 1.5};
    const OperatingPoint op = operatingPoint(c, I2, tb, 300, 0.3, 1.0);
    const double vlpAt = op.ok ? vlpBhp(c, tb, op.q, 300, 0.3, 1.0) : 0;
    std::printf("    operating point: q %.1f stb/d at pwf %.1f psi (VLP %.1f, IPR %.1f) %s\n", op.q, op.pwf, vlpAt, op.ok ? iprRate(I2, op.pwf) : 0.0, op.note.c_str());
    check(op.ok && std::fabs(vlpAt - op.pwf) < 0.5 && std::fabs(iprRate(I2, op.pwf) - op.q) < 0.5, "operating point: IPR and VLP intersect");
    const VfpTable V = makeVfpProd(c, tb, 3, {100, 500, 1500, 3000}, {150, 300}, {0, 0.5}, {0.8, 2.0});
    const double e1 = V.eval(500, 300, 0.5, 2.0, 0), d1 = vlpBhp(c, tb, 500, 300, 0.5, 2.0);
    // between the nodes the table is multilinear: the centre of a cell is the mean of its 16 corners
    double avg = 0;
    for (double r : {500.0, 1500.0}) for (double th : {150.0, 300.0}) for (double w : {0.0, 0.5}) for (double g : {0.8, 2.0}) avg += vlpBhp(c, tb, r, th, w, g) / 16;
    const double e2 = V.eval(1000, 225, 0.25, 1.4, 0), d2 = vlpBhp(c, tb, 1000, 225, 0.25, 1.4);
    std::printf("    VFPPROD %zu points: node %.2f vs %.2f; cell centre %.2f = corner mean %.2f (direct model %.2f psi: the GOR axis is strongly non-linear - use more GOR nodes)\n",
                V.bhp.size(), e1, d1, e2, avg, d2);
    check(V.bhp.size() == 32 && std::fabs(e1 - d1) < 1e-9 && std::fabs(e2 - avg) < 1e-9, "VFPPROD table: exact at the nodes, multilinear between (table layout = VfpTable::eval)");
}

static void testC28FaultSlipEngines() {
    std::printf("== C28 FAULTSLIP in the thermal (EGS) and compositional (CO2) engines - headless coupled runner\n");
    {
        const std::string path = findExample("EGS_FAULT.DATA");
        if (path.empty()) { std::printf("    (EGS_FAULT.DATA not found - skipped)\n"); }
        else {
            Case c = importEclipseDeck(path);
            c.faultMech.slipOn = true; c.faultMech.slipE = 30; c.faultMech.slipNu = 0.25; c.faultMech.slipDMu = 0.05;
            auto runT = [](const Case& cc, coupledbo::Output& co) {
                thermal::DeckModel m; std::string err;
                if (!thermal::deckModel(cc, m, &err)) return false;
                std::vector<thermal::Summary> o;
                return coupleddeck::runThermal(cc, m, o, co, {}, {}, &err, 8);
            };
            coupledbo::Output a, b;
            Case cn = c; cn.faultMech.thE = 0; cn.faultMech.source = 0;   // no thermo-elastic term, analytic poro path
            Case ct = c; ct.faultMech.source = 0;                          // analytic path + thermo-elastic term
            const bool okA = runT(cn, a), okB = runT(ct, b);
            const double sA = coupledbo::kpiMax(a.kpi, "_max_slip_m"), sB = coupledbo::kpiMax(b.kpi, "_max_slip_m");
            const double dA = coupledbo::kpiMax(a.kpi, "_max_dcfs"), dB = coupledbo::kpiMax(b.kpi, "_max_dcfs");
            std::printf("    EGS 8 steps: max dCFS %.4g -> %.4g bar with the thermo-elastic stress; max slip %.4g -> %.4g mm\n", dA, dB, sA * 1e3, sB * 1e3);
            for (const auto& l : b.report) std::printf("   %s\n", l.c_str());
            check(okA && okB && std::isfinite(sA) && std::isfinite(sB), "thermal run with FAULTSLIP: slip KPIs from the headless runner");
            check(dB > dA && sB >= sA, "cooling unclamps the fault: larger dCFS and at least as much slip");
        }
    }
    {
        const std::string path = findExample("CO2_GEOMECH.DATA");
        if (path.empty()) { std::printf("    (CO2_GEOMECH.DATA not found - skipped)\n"); return; }
        Case c = importEclipseDeck(path);
        c.faultMech.slipOn = true; c.faultMech.slipE = 20; c.faultMech.slipNu = 0.25;
        comp::DeckModel m; std::string err;
        const bool okM = comp::deckModel(c, m, &err);
        coupledbo::Output co; std::vector<comp::CompSummary> o;
        const bool ok = okM && coupleddeck::runComp(c, m, o, co, {}, {}, &err, 8);
        const double s = coupledbo::kpiMax(co.kpi, "_max_slip_m"), up = coupledbo::kpiMax(co.kpi, "surface_uplift_m"), ev = coupledbo::kpiSum(co.kpi, "_expected_events");
        std::printf("    CO2 8 steps: ok %d, %zu KPIs: max slip %.4g mm, uplift %.4g mm, expected events %.4g\n", ok, co.kpi.size(), s * 1e3, up * 1e3, ev);
        check(ok && std::isfinite(s) && std::isfinite(up) && std::isfinite(ev), "compositional run with FE geomechanics + FAULTSLIP: slip / uplift / seismicity KPIs");
        // studies on compositional cases now run through the same runner
        Case cs = c; cs.schedule.resize(std::min<size_t>(cs.schedule.size(), 4));
        const RunResult r = runSingle(cs);
        check(r.ok && !r.summary.empty() && !r.kpi.empty(), "runSingle runs compositional cases (studies / MPI batch) with the coupled KPIs");
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// C28 (C21-3) distributed assembly. In the normal suite two "virtual ranks" are threads running the same Simulator in
// lockstep with shared-memory reductions and a gathered serial solve (tests the restricted evaluation / assembly, the
// reduced norms and the lockstep decisions without MPI); under "mpiexec -n P run_tests -mpida" the real collective solver
// (DDSolver::solveCollective) is used on every rank.
#include <thread>
#include <mutex>
#include <condition_variable>
namespace c28da {
struct Barrier {
    std::mutex m; std::condition_variable cv; int n = 2, arrived = 0; long gen = 0;
    void wait() {
        std::unique_lock<std::mutex> lk(m);
        const long g = gen;
        if (++arrived == n) { arrived = 0; ++gen; cv.notify_all(); }
        else cv.wait(lk, [&] { return gen != g; });
    }
};
struct Shared {
    Barrier bar;
    std::vector<double>* slot[2] = {nullptr, nullptr};
    const BSR* S[2] = {nullptr, nullptr}; const std::vector<double>* b[2] = {nullptr, nullptr};
    std::vector<double> x; bool ok = false; int its = 0;
    std::vector<int> starts;
};
}  // namespace c28da

static void testC28DistAssembly() {
    std::printf("== C28 distributed assembly (C21-3): each rank evaluates / assembles its own rows, replicated Newton state\n");
    Case c = makeDemoCase();
    c.schedule.resize(6);
    c.tuning.threads = 1;
    double fo0 = 0, fw0 = 0; int nw0 = 0;
    {
        Simulator sim; sim.log = [](const std::string&) {}; sim.setup(c); sim.initialize();
        while (sim.runReportStep()) {}
        fo0 = sim.summary().back().fopt; fw0 = sim.summary().back().fwpt; nw0 = sim.stats().newton;
    }
    c28da::Shared sh;
    double fo[2] = {0, 0}, fw[2] = {0, 0}; int nw[2] = {0, 0}, solves[2] = {0, 0}; int range[2] = {0, 0}; bool thrown = false;
    auto worker = [&](int k) {
        try {
            Simulator sim; sim.log = [](const std::string&) {}; sim.setup(c);
            DistAssembly d;
            d.starts = sim.cellPartition(2);
            d.a0 = d.starts[size_t(k)]; d.a1 = d.starts[size_t(k) + 1];
            range[k] = d.a1 - d.a0;
            if (k == 0) sh.starts = d.starts;
            auto reduce = [&sh, k](double* v, int n, bool mx) {
                std::vector<double> mine(v, v + n);
                sh.slot[k] = &mine;
                sh.bar.wait();
                for (int i = 0; i < n; ++i) {
                    const double a = (*sh.slot[0])[size_t(i)], b = (*sh.slot[1])[size_t(i)];
                    v[i] = mx ? std::max(a, b) : a + b;
                }
                sh.bar.wait();
            };
            d.sumN = [reduce](double* v, int n) { reduce(v, n, false); };
            d.maxN = [reduce](double* v, int n) { reduce(v, n, true); };
            d.solve = [&sh, k](const BSR& S, const std::vector<double>& b, std::vector<double>& x, double tol, int maxIt, int& its) {
                sh.S[k] = &S; sh.b[k] = &b;
                sh.bar.wait();
                if (k == 0) {   // gather the owned rows of both "ranks" and solve serially
                    BSR F = *sh.S[0];
                    std::vector<double> bf = *sh.b[0];
                    const size_t b2 = size_t(F.nb) * F.nb;
                    const int r0 = sh.starts[1], r1 = sh.starts[2];
                    const BSR& S1 = *sh.S[1];
                    std::copy(S1.val.begin() + size_t(S1.rowPtr[size_t(r0)]) * b2, S1.val.begin() + size_t(S1.rowPtr[size_t(r1)]) * b2,
                              F.val.begin() + size_t(F.rowPtr[size_t(r0)]) * b2);
                    std::copy(sh.b[1]->begin() + size_t(r0) * F.nb, sh.b[1]->begin() + size_t(r1) * F.nb, bf.begin() + size_t(r0) * F.nb);
                    BILU0 M; M.factor(F); BsrOperator op(F);
                    sh.x.clear();
                    const LinResult lr = bicgstab(op, M, bf, sh.x, tol, std::max(maxIt, 400));
                    sh.ok = lr.converged; sh.its = lr.iterations;
                }
                sh.bar.wait();
                x = sh.x; its = sh.its;
                const bool ok = sh.ok;
                sh.bar.wait();
                return ok;
            };
            sim.setDistributedAssembly(d);
            sim.initialize();
            while (sim.runReportStep()) {}
            fo[k] = sim.summary().back().fopt; fw[k] = sim.summary().back().fwpt; nw[k] = sim.stats().newton; solves[k] = sim.stats().ddSolves;
        } catch (...) { thrown = true; }
    };
    std::thread t1(worker, 1);
    worker(0);
    t1.join();
    std::printf("  demo 6 steps: FOPT serial %.9g | rank 0 %.9g | rank 1 %.9g; Newton %d / %d / %d; rows %d + %d; collective solves %d\n",
                fo0, fo[0], fo[1], nw0, nw[0], nw[1], range[0], range[1], solves[0]);
    check(!thrown && range[0] > 0 && range[1] > 0, "two-rank partition of the cells, both ranks own rows");
    check(fo[0] == fo[1] && fw[0] == fw[1] && nw[0] == nw[1], "lockstep: both ranks take identical decisions (bitwise equal results)");
    check(solves[0] > 0 && std::fabs(fo[0] - fo0) < 1e-5 * fo0 && std::fabs(fw[0] - fw0) <= 1e-4 * std::max(1.0, fw0),
          "distributed-assembly run = serial run (FOPT / FWPT)");
}

// mpiexec -n P run_tests -mpida : the real collective solve over MPI
static void testC28DistAssemblyMpi() {
    const int me = mpi::rank(), P = mpi::size();
    if (me == 0) std::printf("== C28 distributed assembly over MPI (%d ranks): own rows, no matrix scatter, all-gathered update\n", P);
    Case c = makeDemoCase();
    c.schedule.resize(6);
    c.tuning.threads = 1;
    double fo0 = 0; int nw0 = 0;
    {
        Simulator sim; sim.log = [](const std::string&) {}; sim.setup(c); sim.initialize();
        while (sim.runReportStep()) {}
        fo0 = sim.summary().back().fopt; nw0 = sim.stats().newton;
    }
    mpi::DDSolver dd;
    Simulator sim; sim.log = [](const std::string&) {}; sim.setup(c);
    DistAssembly d;
    d.starts = sim.cellPartition(P);
    d.a0 = d.starts[size_t(me)]; d.a1 = d.starts[size_t(me) + 1];
    d.sumN = [](double* v, int n) { mpi::allreduceSumN(v, n); };
    d.maxN = [](double* v, int n) { for (int i = 0; i < n; ++i) v[i] = mpi::allreduceMax(v[i]); };
    const std::vector<int> st = d.starts;
    d.solve = [&dd, st](const BSR& S, const std::vector<double>& b, std::vector<double>& x, double tol, int maxIt, int& its) {
        return dd.solveCollective(S, b, x, st, tol, maxIt, its);
    };
    sim.setDistributedAssembly(d);
    const auto t0 = std::chrono::steady_clock::now();
    sim.initialize();
    while (sim.runReportStep()) {}
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const double fo = sim.summary().back().fopt;
    double foMax = mpi::allreduceMax(fo), foMin = -mpi::allreduceMax(-fo);
    if (me == 0) {
        std::printf("  FOPT serial %.9g / distributed %.9g (spread over ranks %.2e), Newton %d / %d, %d collective solves (%d its), %.3f s, pattern set-ups %d\n",
                    fo0, fo, foMax - foMin, nw0, sim.stats().newton, dd.stats.solves, dd.stats.iterations, sec, dd.stats.patternSetups);
        check(foMax == foMin, "all ranks hold the identical replicated state");
        check(dd.stats.solves > 0 && std::fabs(fo - fo0) < 1e-5 * fo0, "MPI distributed-assembly run = serial run");
    }
}

// C28 (C23-b) level-ordered contiguous BILU0: same ILU(0) (same iterations, same solution), timing of the sweeps
static void c28SetEnv(const char* k, const char* v) {
#ifdef _WIN32
    _putenv_s(k, v);
#else
    setenv(k, v, 1);
#endif
}
static void testC28LevelIlu() {
    std::printf("== C28 level-ordered contiguous BILU0 (C23-b): renumbering by forward level keeps ILU(0) exact\n");
    const BSR A = makeBlockLaplace(240, 225, 3, 7);
    const int N = A.size();
    std::vector<double> xt(static_cast<size_t>(N)), b(static_cast<size_t>(N));
    for (int e = 0; e < N; ++e) xt[size_t(e)] = std::cos(0.001 * e) + 0.2 * (e % 3);
    A.spmv(xt.data(), b.data());
    const char* modes[3] = {"0", "1", "2"};
    const char* names[3] = {"serial sweep", "C21 level-scheduled (scattered rows)", "C28 level-ordered contiguous"};
    LinResult res[3]; double tApply[3] = {0, 0, 0}, err[3] = {0, 0, 0}; bool contig = false; int levels = 0;
    std::vector<double> x0;
    for (int m = 0; m < 3; ++m) {
        c28SetEnv("RESSIM_LEVEL_ILU", modes[m]);
        BILU0 M; M.factor(A);
        if (m == 2) { contig = M.contiguousLevels(); levels = M.forwardLevels(); }
        BsrOperator op(A); std::vector<double> x;
        res[m] = bicgstab(op, M, b, x, 1e-10, 2000);
        std::vector<double> z(static_cast<size_t>(N));
        const auto t0 = std::chrono::steady_clock::now();
        for (int r = 0; r < 40; ++r) M.apply(b.data(), z.data());
        tApply[m] = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / 40.0;
        if (m == 0) x0 = z;
        for (int e = 0; e < N; ++e) err[m] = std::max(err[m], std::fabs(z[size_t(e)] - x0[size_t(e)]) / (1e-30 + std::fabs(x0[size_t(e)]) + 1e-12));
    }
    c28SetEnv("RESSIM_LEVEL_ILU", "0");
    for (int m = 0; m < 3; ++m)
        std::printf("  %-40s %4d its, apply %.2f ms, max rel. difference of M^-1 b vs serial %.1e\n", names[m], res[m].iterations, 1e3 * tApply[m], err[m]);
    std::printf("  %d forward levels, %d rows per level on average\n", levels, levels > 0 ? A.n / levels : 0);
    check(contig, "contiguous level storage built (structurally symmetric pattern, topological renumbering)");
    check(res[2].converged && std::abs(res[2].iterations - res[0].iterations) <= 1 && err[2] < 1e-9,
          "renumbered ILU(0) = original ILU(0): same iterations, same preconditioned vector (round-off only)");
}

// C28 phase field in the field run (C26-4 open item): HFGROW ... PHASEFIELD eta - a persistent volume-controlled phase-field
// model of the fracture plane advanced once per report step by the injected volume, its half-length opens the EDFM pieces
static void testC28PhaseFieldInRun() {
    std::printf("== C28 phase-field fracture growth inside the flow run (HFGROW ... PHASEFIELD)\n");
    auto deck = [](const std::string& grow) {
        std::string d = hfDeck(true, false, grow);
        const std::string a = "WCONPROD\n P OPEN ORAT 5 4* 10 /";
        const size_t at = d.find(a);
        if (at != std::string::npos) d.replace(at, a.size(), "WCONINJE\n P WATER OPEN RATE 2 1* 700 /");
        const size_t q = d.find(" 0.2 0 1 0\n");
        if (q != std::string::npos) d.replace(q, 11, " 0.2 0 1 0\n 0.3 0.02 0.7 0\n");
        const size_t ts = d.find("TSTEP\n");
        if (ts != std::string::npos) d.replace(ts, d.find("END") - ts, "TSTEP\n 4*1 6*4 6*10 /\n");
        return d;
    };
    const double eta = 0.002;
    char g[160]; std::snprintf(g, sizeof g, "-- RESSIM: HFGROW F1 1.0 320 0 20 1 PHASEFIELD %.6g\n", eta);
    bool ok = false; std::string err;
    std::vector<Simulator::PfRecord> rec; double wl[2] = {0, 0}, wit = 0; int ev = 0;
    std::string roundTrip;
    const auto t0 = std::chrono::steady_clock::now();
    try {
        Case c = eclipseDeckFromString(deck(g), "");
        roundTrip = projectText(c);
        Simulator s; s.log = [](const std::string&) {};
        s.setup(c); s.initialize();
        while (s.runReportStep()) {}
        rec = s.phaseFieldRecords(); wl[0] = s.hfWingLength(0, -1); wl[1] = s.hfWingLength(0, 1);
        wit = s.summary().back().fwit; ev = int(s.hfGrowthEvents().size());
        ok = true;
    } catch (const std::exception& e) { err = e.what(); }
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (!ok) std::printf("    exception: %s\n", err.c_str());
    // toughness-regime plane-strain crack: V = 2 sqrt(pi) K a^1.5 / E'  ->  a(V) = (E' V / (2 sqrt(pi) K))^(2/3)
    const double Ep = 20e9 / (1 - 0.25 * 0.25), K = 1e6;
    bool mono = !rec.empty(), pOk = true; double lastA = 0, ratio = 0; int act = 0;
    for (const auto& r : rec) {
        if (!r.active || r.volume <= 0) continue;
        ++act;
        mono = mono && r.halfLength >= lastA - 1e-9; lastA = r.halfLength;
        pOk = pOk && r.pressure > 320e5;
    }
    if (act > 0) {
        const auto& r = rec.back();
        const double aT = std::pow(Ep * r.volume / (2 * std::sqrt(M_PI) * K), 2.0 / 3.0);
        ratio = r.halfLength / std::max(aT, 20.0);
        std::printf("    %d report steps, %.1f s: injected %.4g m3, V = eta Vinj / H = %.4g m2/m, phase-field half-length %.2f m (toughness a(V) %.2f m),"
                    " crack pressure %.4g bar (closure 320, flow model %.4g bar), wings %.1f / %.1f m, %d growth events\n",
                    int(rec.size()), sec, wit, r.volume, r.halfLength, aT, r.pressure / 1e5, r.pFlow, wl[0], wl[1], ev);
    }
    check(ok && act >= 3 && roundTrip.find("PHASEFIELD") != std::string::npos, "deck PHASEFIELD option parsed / written, phase-field steps inside the run");
    check(ok && mono && pOk && lastA > 20.0, "crack grows monotonically (irreversible phase field) under the injected volume, pressure above the closure stress");
    check(ok && ev > 0 && wl[0] > 20 && std::fabs(wl[0] - wl[1]) <= 10.5 && std::fabs(wl[1] - lastA) <= 10.5,
          "phase-field half-length opens the EDFM pieces of both wings (within one piece)");
    {   // C28 H2b: SUBGRID LEAK flag (sub-grid leak-off) parsed and written
        bool rt = false;
        try {
            Case c2 = eclipseDeckFromString(deck("-- RESSIM: HFGROW F1 1.0 320 0 20 1 SUBGRID LEAK\n"), "");
            rt = !c2.grid.hfracs.empty() && c2.grid.hfracs[0].subgrid == 1 && c2.grid.hfracs[0].subLeak == 1 &&
                 projectText(c2).find("SUBGRID LEAK") != std::string::npos;
        } catch (...) {}
        check(rt, "HFGROW SUBGRID LEAK (H2b sub-grid leak-off) deck round trip");
    }
    check(ok && ratio > 0.75 && ratio < 1.25, "phase-field length within the O(l/a) band (+-25 %) of the toughness-regime Sneddon-Griffith length");
}
