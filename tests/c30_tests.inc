// =====================================================================================================================
// C30 tests: implicit surface-network balance (no limit cycle), fracture field output (HFFIELD: DDM opening, traction
// boundary condition, displacement jump), phase-field field output (crack-opening integral = crack volume), FE
// geomechanics post-processing of finished results (= in-run one-way coupling, elastic linearity, ECLIPSE binary round
// trip) and the colour maps of the viewers
// =====================================================================================================================
#include "ressim/hf_field.h"
#include "ressim/phasefield.h"
#include "../app/src/colormaps.cpp"

static void testC30Network() {
    std::printf("== C30 surface network: implicit balance of the node pressures (NETWORK.DATA)\n");
    const std::string path = findExample("NETWORK.DATA");
    if (path.empty()) { std::printf("    (examples/NETWORK.DATA not found - skipped)\n"); return; }
    Case c = importEclipseDeck(path);
    Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
    const auto t0 = std::chrono::steady_clock::now();
    while (s.runReportStep()) {}
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const auto& S = s.summary();
    // reversals of the field oil rate after the plateau (the C29 explicit balance flip-flopped 3493 -> 49 -> 2584 -> 20 stb/d)
    int rev = 0, sgnPrev = 0; double minRatio = 1e30, q2557 = 0;
    for (size_t k = 1; k < S.size(); ++k) {
        const double a = S[k - 1].fopr, b = S[k].fopr, d = b - a;
        if (S[k].time > 760 && std::fabs(d) > 0.02 * std::max(a, b)) {
            const int sg = d > 0 ? 1 : -1;
            if (sgnPrev != 0 && sg != sgnPrev) ++rev;
            sgnPrev = sg;
        }
        if (S[k].time > 760 && S[k].time < 3000 && a > 0) minRatio = std::min(minRatio, b / a);
        if (S[k].time <= 2557.5) q2557 = b;
    }
    std::printf("    %zu steps, %.1f s: FOPT %.6g stb, FOPR at 2557 d %.5g stb/d, rate reversals after the plateau %d, smallest step ratio %.3f, %zu nodes\n",
                S.size(), sec, S.back().fopt, q2557, rev, minRatio, s.networkNodes().size());
    check(!s.networkNodes().empty(), "network read (NODEPROP / BRANPROP)");
    check(rev <= 1 && minRatio > 0.5, "no limit cycle: the field rate declines monotonically after the plateau (no dead / revived wells)");
    check(q2557 > 300 && S.back().fopt > 1.0e7, "the wells keep producing through the network (the C29 explicit balance killed them at ~6.5 y, FOPT 7.1e6)");
}

static void testC30HfField() {
    std::printf("== C30 fracture field output (HFFIELD): DDM opening, traction condition, displacement jump, arrays\n");
    const std::string path = findExample("HFGROW_INJECTION.DATA");
    if (path.empty()) { std::printf("    (examples/HFGROW_INJECTION.DATA not found - skipped)\n"); return; }
    Case c = importEclipseDeck(path);
    {   // deck round trip of the HFFIELD keyword
        Case d = c; d.hfField.on = 1; d.hfField.E = 15; d.hfField.nu = 0.21; d.hfField.radius = 1.5; d.hfField.budget = 3;
        const std::string txt = projectText(d);
        bool rt = false;
        try { const Case e = projectFromText(txt); rt = e.hfField.on == 1 && e.hfField.E == 15 && e.hfField.nu == 0.21 && e.hfField.radius == 1.5 && e.hfField.budget == 3; } catch (...) {}
        check(txt.find("HFFIELD ON") != std::string::npos && rt, "HFFIELD ON E nu radius budget round-trips through the deck");
    }
    Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
    while (s.runReportStep()) {}
    const hffield::Stats* st = s.hfFieldStats();
    check(st != nullptr && st->updates > 0 && st->open > 0 && st->maxWidth > 0, "field output active for a growing fracture (auto)");
    if (!st || st->open == 0) return;
    const Snapshot& sn = s.snapshots().back();
    const size_t na = size_t(s.grid().numActive());
    bool arrays = true;
    for (const char* nm : {"HFWIDTH", "HFPNET", "HFDSNN", "HFDSXX", "HFDISP", "HFSTRV"}) arrays = arrays && sn.field(nm).size() == na;
    check(arrays, "HFWIDTH / HFPNET / HFDSNN / HFDSXX / HFDISP / HFSTRV arrays on the active cells");
    // probe at the widest open piece: sigma_nn just off the face = -p_net (tension +), displacement jump = opening
    const Grid& G = s.grid();
    const auto& W = s.caseData().grid.hfPieces;
    const double lm = c.units == UnitSystem::Field ? 0.3048 : 1.0, pa = c.units == UnitSystem::Field ? 6894.757 : 1e5;
    std::vector<double> width, pnet;
    {   // per fracture cell from the arrays (HFWIDTH mm, HFPNET case units)
        const auto& A = sn.field("HFWIDTH"); const auto& P = sn.field("HFPNET");
        for (int k = 0; k < G.numHfCells(); ++k) {
            const int a = G.activeOf[size_t(G.hfCells[size_t(k)].global)];
            width.push_back(a >= 0 ? double(A[size_t(a)]) * 1e-3 : 0); pnet.push_back(a >= 0 ? double(P[size_t(a)]) : 0);
        }
    }
    size_t kb = 0;
    for (size_t k = 0; k < width.size(); ++k) if (width[k] > width[kb]) kb = k;
    if (kb >= W.size()) { check(false, "fracture pieces available"); return; }
    const auto& Pc = W[kb];
    const double L = std::hypot(Pc.xb - Pc.xa, Pc.yb - Pc.ya) * lm, H = std::fabs(Pc.zb - Pc.zt) * lm;
    const double nx = -(Pc.yb - Pc.ya) * lm / L, ny = (Pc.xb - Pc.xa) * lm / L, dl = 0.01 * std::min(L, H);
    // a quarter along the piece (the centre lies on the diagonal shared by the two triangles)
    const double xc = (0.375 * Pc.xa + 0.625 * Pc.xb) * lm, yc = (0.375 * Pc.ya + 0.625 * Pc.yb) * lm, zc = -0.5 * (Pc.zt + Pc.zb) * lm;
    double sn2 = 0, du = 0;
    for (int side : {-1, 1}) {
        const double x[3] = {xc + side * dl * nx, yc + side * dl * ny, zc};
        double Sx[6], u[3];
        if (!s.hfFieldAt(x, Sx, u)) { check(false, "field probe"); return; }
        sn2 += 0.5 * (Sx[0] * nx * nx + Sx[1] * ny * ny + 2 * Sx[3] * nx * ny);
        du += side * (u[0] * nx + u[1] * ny);
    }
    const double pn = pnet[kb] * pa;
    std::printf("    %d pieces, %d open, %d cells, %s space, E %.3g GPa | widest piece %zu: w %.4g mm, p_net %.4g Pa, sigma_nn off the face %.4g Pa, displacement jump %.4g mm\n",
                st->pieces, st->open, st->cells, st->halfSpace ? "half" : "full", st->E * 1e-9, kb, width[kb] * 1e3, pn, sn2, du * 1e3);
    check(pn > 0 && std::fabs(sn2 + pn) < 0.15 * pn, "traction condition: normal-stress change on the open face = -p_net (15 %)");
    check(std::fabs(du - width[kb]) < 0.15 * width[kb], "displacement jump across the face = DDM opening (15 %)");
}

static void testC30PhaseFieldOutput() {
    std::printf("== C30 phase-field field output: crack-opening integral, damage, stress change\n");
    phasefield::Input in; in.h = 0.01; in.control = 1; in.steps = 10;
    {
        const double Ep = in.E / (1 - in.nu * in.nu), pg = std::sqrt(Ep * in.Gc / (M_PI * in.crackA));
        in.vMax = 3 * 2 * M_PI * pg * in.crackA * in.crackA / Ep;
    }
    std::vector<phasefield::Field> F;
    in.onField = [&F](const phasefield::Field& f) { F.push_back(f); };
    phasefield::Output o; phasefield::run(in, o);
    check(o.ok && !F.empty() && F.size() == o.steps.size(), "one field per completed load step");
    if (F.empty()) return;
    double worst = 0, last = 0;
    for (const auto& f : F) {
        double I = 0;
        for (size_t q = 0; q + 1 < f.codX.size(); ++q) I += 0.5 * (f.cod[q] + f.cod[q + 1]) * (f.codX[q + 1] - f.codX[q]);
        last = I / std::max(1e-30, f.volume);
        worst = std::max(worst, std::fabs(last - 1));
    }
    double sned = 0;
    {   // independent check: the first (pre-growth) step against the Sneddon opening w(0) = 4 p a / E'
        const auto& f0 = F.front();
        const double Ep = in.E / (1 - in.nu * in.nu);
        double wc = 0;
        size_t qc = 0; for (size_t q = 0; q < f0.codX.size(); ++q) if (std::fabs(f0.codX[q]) < std::fabs(f0.codX[qc])) qc = q;
        if (!f0.cod.empty()) wc = f0.cod[qc];
        sned = wc / (4 * f0.p * std::max(1e-9, f0.half) / Ep);
    }
    const auto& f = F.back();
    const size_t nx1 = f.x.size(), ny1 = f.y.size();
    float vmin = 1; for (float v : f.v) vmin = std::min(vmin, v);
    double smax = 0; for (float v : f.syy) smax = std::max(smax, std::fabs(double(v)));
    std::printf("    %zu steps, %zu x %zu nodes: int cod dx / V = %.4f (worst step %.3g), half-length %.4g m, min phase field %.3g, max |d sigma_yy| %.3g MPa\n",
                F.size(), nx1, ny1, last, worst, f.half, vmin, smax * 1e-6);
    check(f.v.size() == nx1 * ny1 && f.sxx.size() == (nx1 - 1) * (ny1 - 1) && f.ux.size() == nx1 * ny1, "field array sizes (nodes / elements)");
    std::printf("    first step: centre opening / Sneddon 4 p a / E' = %.3f\n", sned);
    check(worst < 0.05, "crack-opening profile integrates to the crack volume (5 %, every step; consistency of the output)");
    check(sned > 0.85 && sned < 1.3, "first-step opening at the centre within -15 / +30 % of Sneddon (regularised crack is wider)");
    check(vmin < 0.05 && smax > 0, "damage and stress change present");
}

static void testC30PostGeomech() {
    std::printf("== C30 FE geomechanics as post-processing: = in-run one-way coupling, elastic linearity, ECLIPSE binary round trip\n");
    Case c = smallCase();
    c.schedule.resize(8);
    geomech::Input in; in.enabled = true; in.coupling = 0;
    in.E = {10e9}; in.nu = {0.25}; in.biot = {1.0}; in.alphaT = {1e-5}; in.Ks = 0; in.sigmaVGrad = 23e3; in.k0 = 0.7;
    in.overLayers = 2; in.underLayers = 1; in.sidePad = 1;
    Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
    const Grid& G = s.grid();
    const size_t na = size_t(G.numActive());
    const bool fu = c.units == UnitSystem::Field;
    const double paU = fu ? 6894.757 : 1e5, lm = fu ? 0.3048 : 1.0;
    // in-run one-way coupler stepping with the report steps (reference)
    geoseq::Coupler C; std::string err;
    std::vector<double> p(na), phi(na, 0.2);
    auto pres = [&]() { const auto& P = s.snapshots().back().pressure; for (size_t a = 0; a < na; ++a) p[a] = double(P[a]) * paU; };
    pres();
    for (size_t a = 0; a < na; ++a) { const int g = G.globalOf[a]; if (G.bulkVol[size_t(g)] > 0) phi[a] = std::clamp(G.poreVolRef[a] / G.bulkVol[size_t(g)], 0.01, 0.6); }
    if (!C.setup(in, G, G.globalOf, lm, p, phi, err)) { check(false, "coupler setup"); std::printf("    %s\n", err.c_str()); return; }
    std::vector<double> upRef;
    while (s.runReportStep()) { pres(); C.step(p); }
    for (size_t a = 0; a < na; ++a) { const int b = C.boxOf(int(a)); upRef.push_back(b >= 0 ? C.solver().dispZ()[size_t(b)] : 0.0); }
    std::vector<Snapshot> snaps = s.snapshots();
    const geoseq::PostResult R = geoseq::postProcess(in, G, fu, snaps);
    bool all = R.error.empty() && R.steps == int(snaps.size());
    for (const auto& sn : snaps) all = all && sn.field("UPLIFT").size() == na && sn.field("SIGEFFM").size() == na && sn.field("FAILURE").size() == na;
    double d = 0, m = 0, u0 = 0;
    const auto& U = snaps.back().field("UPLIFT");
    for (size_t a = 0; a < na && U.size() == na; ++a) { d = std::max(d, std::fabs(double(U[a]) - (-upRef[a]))); m = std::max(m, std::fabs(upRef[a])); }
    for (float v : snaps.front().field("UPLIFT")) u0 = std::max(u0, std::fabs(double(v)));
    std::printf("    %d snapshots in %.2f s: max uplift %.4g m, subsidence %.4g m, failure %.3g | |post - in-run| %.3g of %.3g m, first step %.2g m\n",
                R.steps, R.seconds, R.maxUplift, R.maxSubsidence, R.maxFailure, d, m, u0);
    check(all, "UPLIFT / SIGEFFM / FAILURE ... on every snapshot");
    // the vertical displacement sign convention differs (dispZ + down, UPLIFT + up): compare either sign, the test reports it
    double dAlt = 0; for (size_t a = 0; a < na && U.size() == na; ++a) dAlt = std::max(dAlt, std::fabs(double(U[a]) - upRef[a]));
    std::printf("    UPLIFT vs solver dispZ: %s\n", d < dAlt ? "UPLIFT = -dispZ (dispZ positive down)" : "UPLIFT = dispZ");
    check(m > 0 && std::min(d, dAlt) <= 1e-4 * m + 1e-9, "post-processing = in-run one-way coupling (same mesh, same steps)");
    check(R.maxSubsidence > 0 || R.maxUplift > 0, "uplift / subsidence statistics of the reservoir cells");
    check(u0 < 1e-12, "the first snapshot is the reference state (zero displacement)");
    {   // linear elasticity: doubling the pressure change doubles the displacement
        std::vector<Snapshot> two = s.snapshots();
        for (auto& sn : two) for (size_t a = 0; a < na; ++a) sn.pressure[a] = float(2.0 * double(sn.pressure[a]) - double(two.front().pressure[a]));
        for (size_t a = 0; a < na; ++a) two.front().pressure[a] = s.snapshots().front().pressure[a];
        geoseq::postProcess(in, G, fu, two);
        const auto& U2 = two.back().field("UPLIFT");
        double e = 0; for (size_t a = 0; a < na && U2.size() == na; ++a) e = std::max(e, std::fabs(double(U2[a]) - 2 * double(U[a])));
        check(U2.size() == na && e <= 2e-3 * std::max(1e-12, 2 * m), "elastic: twice the pressure change -> twice the uplift");
    }
    {   // ECLIPSE binary round trip (the same path an imported OPM / ECLIPSE case takes)
        const std::string base = "/tmp/ressim_c30_post";
        std::vector<std::string> names; for (const auto& w : c.wells) names.push_back(w.name);
        std::string e2; bool ok = ecl::writeResultCase(base, c, G, s.summary(), s.snapshots(), names, &e2);
        ecl::ResultCase rc;
        ok = ok && ecl::openResultCase(base + ".EGRID", rc, &e2) && rc.grid;
        double dr = 1e30; size_t nr = 0;
        if (ok) {
            std::vector<Snapshot> rs = ecl::snapshotsFromRestart(rc, -1, 0);
            nr = rs.size();
            const geoseq::PostResult R2 = geoseq::postProcess(in, *rc.grid, rc.unitSystem == 2, rs);
            if (R2.error.empty() && !rs.empty() && rs.back().field("UPLIFT").size() == na) {
                dr = 0; const auto& U3 = rs.back().field("UPLIFT");
                for (size_t a = 0; a < na; ++a) dr = std::max(dr, std::fabs(double(U3[a]) - double(U[a])));
            }
        }
        std::printf("    ECLIPSE round trip: %zu restart steps, units %d, |post(restart) - post(run)| %.3g m %s\n", nr, rc.unitSystem, dr, e2.c_str());
        check(ok && nr == snaps.size() && dr <= 1e-3 * m + 1e-9, "post-processing of the written EGRID / UNRST case = post-processing of the run");
        for (const char* ext : {".EGRID", ".INIT", ".SMSPEC", ".UNSMRY", ".UNRST"}) std::remove((base + ext).c_str());
    }
}

static void testC30Colormaps() {
    std::printf("== C30 colour maps of the viewers\n");
    const int n = int(Colormap::Count);
    bool range = true, rev = true, names = true; double jumpMax = 0; std::string worst;
    std::set<std::string> seen;
    for (int m = 0; m < n; ++m) {
        const Colormap cm = Colormap(m);
        names = names && colormapName(cm) && seen.insert(colormapName(cm)).second;
        float prev[3]; colormapRGB(cm, 0, prev);
        for (int i = 0; i <= 512; ++i) {
            const float t = float(i) / 512;
            float a[3], b[3], r[3];
            colormapRGB(cm, t, a); colormapRGB(cm, 1 - t, b); colormapRGB(cm, t, r, true);
            for (int k = 0; k < 3; ++k) {
                range = range && a[k] >= 0 && a[k] <= 1;
                rev = rev && std::fabs(r[k] - b[k]) < 1e-6f;
                const double j = std::fabs(a[k] - prev[k]);
                if (j > jumpMax) { jumpMax = j; worst = colormapName(cm); }
                prev[k] = a[k];
            }
        }
    }
    auto lum = [](Colormap cm, float t) { float c[3]; colormapRGB(cm, t, c); return 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2]; };
    bool mono = true;
    for (Colormap cm : {Colormap::Viridis, Colormap::Plasma, Colormap::Inferno, Colormap::Magma, Colormap::Cividis, Colormap::Grey})
        for (int i = 1; i <= 64; ++i) mono = mono && lum(cm, i / 64.f) > lum(cm, (i - 1) / 64.f) - 1e-4;
    bool div = true;
    for (int m = 0; m < n; ++m) if (colormapDiverging(Colormap(m))) div = div && lum(Colormap(m), 0.5f) > std::min(lum(Colormap(m), 0.f), lum(Colormap(m), 1.f));
    std::printf("    %d maps, largest step between 513 samples %.3f (%s)\n", n, jumpMax, worst.c_str());
    check(n >= 15 && names, "15 named colour maps");
    check(range && rev, "rgb in [0, 1]; reversed map = map(1 - t)");
    check(jumpMax < 0.05, "continuous (no jumps between neighbouring samples)");
    check(mono, "perceptual maps (viridis, plasma, inferno, magma, cividis, grey): monotonic luminance");
    check(div, "diverging maps are lighter at the centre than at the darker end");
}
