// C36 tests: out-of-core SEG-Y, seismic attributes, flattening (seismic3d.h)
#include "ressim/seismic3d.h"

static void testC36Seismic() {
    namespace fs = std::filesystem;
    using namespace ressim::geo;
    std::printf("== C36 seismic: out-of-core SEG-Y, attributes, flattening\n");
    const std::string dir = (fs::temp_directory_path() / "ressim_c36_seis").string();
    fs::create_directories(dir);
    const int NI = 40, NJ = 30, NS = 200;
    std::string e;
    const bool w5 = writeSyntheticSegy(dir + "/s5.sgy", NI, NJ, NS, 4.0, 5, false, e);
    const SeismicCube C = readSegy(dir + "/s5.sgy", e);
    SegyVolume V;
    const bool ok5 = w5 && V.open(dir + "/s5.sgy", e);
    {
        std::printf("    volume %dx%dx%d il0 %d dil %d xl0 %d, origin %.2f %.2f, inline step (%.3f, %.3f), xline step (%.3f, %.3f), fit error %.2e, |a|max %.3f rms %.3f\n",
                    V.ni, V.nj, V.ns, V.il0, V.dil, V.xl0, V.x0, V.y0, V.ix, V.iy, V.jx, V.jy, V.fitError, V.ampAbsMax, V.ampRms);
        check(ok5 && V.ni == NI && V.nj == NJ && V.ns == NS && V.il0 == 100 && V.dil == 2 && V.xl0 == 300 && V.dxl == 1 && std::fabs(V.dt - 4) < 1e-12 &&
              std::fabs(V.ix - 20) < 1e-3 && std::fabs(V.iy - 15) < 1e-3 && std::fabs(V.jx + 15) < 1e-3 && std::fabs(V.jy - 20) < 1e-3 && V.fitError < 0.01 && V.missing == 0,
              "C36: out-of-core SEG-Y index - lattice, inline step 2, rotated plan geometry by least squares over all traces");
    }
    if (ok5 && C.ni == NI) {
        double d = 0;
        const std::vector<float> il = V.inlineSlice(17), xl = V.xlineSlice(11), ts = V.timeSlice(93);
        for (int j = 0; j < NJ; ++j) for (int s = 0; s < NS; ++s) d = std::max(d, double(std::fabs(il[size_t(j * NS + s)] - C.at(17, j, s))));
        for (int i = 0; i < NI; ++i) for (int s = 0; s < NS; ++s) d = std::max(d, double(std::fabs(xl[size_t(i * NS + s)] - C.at(i, 11, s))));
        for (int i = 0; i < NI; ++i) for (int j = 0; j < NJ; ++j) d = std::max(d, double(std::fabs(ts[size_t(i * NJ + j)] - C.at(i, j, 93))));
        std::vector<float> tr(NS); V.trace(5, 7, tr.data());
        for (int s = 0; s < NS; ++s) d = std::max(d, double(std::fabs(tr[size_t(s)] - C.at(5, 7, s))));
        const SeismicCube S = V.subCube(3, 9, 4, 20, 50, 120);
        for (int i = 0; i < S.ni; ++i) for (int j = 0; j < S.nj; ++j) for (int s = 0; s < S.ns; ++s) d = std::max(d, double(std::fabs(S.at(i, j, s) - C.at(i + 3, j + 4, s + 50))));
        std::printf("    inline / crossline / time slice / trace / sub-cube vs the in-memory reader: max diff %.1e; sub-cube %dx%dx%d t0 %.0f il0 %d\n", d, S.ni, S.nj, S.ns, S.t0, S.il0);
        check(d == 0 && S.ni == 7 && S.nj == 17 && S.ns == 71 && std::fabs(S.t0 - 200) < 1e-9 && S.il0 == 106,
              "C36: slices, traces and sub-cubes read on demand equal the in-memory SEG-Y");
        const SeismicCube P = V.proxy(14, 50);
        double dp = 0;
        for (int a = 0; a < P.ni; ++a) for (int b = 0; b < P.nj; ++b) for (int s = 0; s < P.ns; ++s) {
            double m = 0; for (int q = 0; q < 4; ++q) m += C.at(a * 3, b * 3, s * 4 + q);
            dp = std::max(dp, std::fabs(m / 4 - P.at(a, b, s)));
        }
        std::printf("    proxy %dx%dx%d (stride 3 / 3 / 4), dt %.0f, max diff vs box average %.1e\n", P.ni, P.nj, P.ns, P.dt, dp);
        check(P.ni == 14 && P.nj == 10 && P.ns == 50 && std::fabs(P.dt - 16) < 1e-12 && dp < 1e-5, "C36: decimated proxy cube (plan stride + time box average) for volume rendering");
    }
    {   // other sample formats and byte orders read as the same signal
        struct F { int fmt; bool le; double minCorr; };
        const F fl[] = {{1, false, 0.99999}, {2, false, 0.99999}, {3, false, 0.9999}, {8, false, 0.99}, {5, true, 1.0 - 1e-12}, {3, true, 0.9999}};
        bool all = true;
        for (const F& f : fl) {
            const std::string p = dir + "/f" + std::to_string(f.fmt) + (f.le ? "le" : "be") + ".sgy";
            SegyVolume X;
            const bool ok = writeSyntheticSegy(p, 12, 10, 120, 4.0, f.fmt, f.le, e) && X.open(p, e);
            double sab = 0, saa = 0, sbb = 0;
            if (ok) {
                const SeismicCube A = X.toCube();
                SegyVolume R; std::string e2; writeSyntheticSegy(dir + "/ref.sgy", 12, 10, 120, 4.0, 5, false, e2); R.open(dir + "/ref.sgy", e2);
                const SeismicCube B = R.toCube();
                for (size_t k = 0; k < A.a.size() && k < B.a.size(); ++k) { sab += double(A.a[k]) * B.a[k]; saa += double(A.a[k]) * A.a[k]; sbb += double(B.a[k]) * B.a[k]; }
            }
            const double corr = saa > 0 && sbb > 0 ? sab / std::sqrt(saa * sbb) : 0;
            std::printf("    format %d %s: %s, little-endian detected %d, correlation with IEEE %.7f\n", f.fmt, f.le ? "LE" : "BE", ok ? "read" : e.c_str(), int(X.littleEndian), corr);
            all = all && ok && X.littleEndian == f.le && corr >= f.minCorr;
        }
        check(all, "C36: SEG-Y sample formats 1 / 2 / 3 / 5 / 8 and byte-swapped (little-endian) files");
    }
    {   // a trace missing from the file: index marks it, slices read zeros there
        std::ifstream in(dir + "/s5.sgy", std::ios::binary);
        std::vector<char> b((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        const size_t tl = 240 + 4 * NS, k = 3600 + tl * size_t(4 * NJ + 6);   // trace (4, 6)
        std::vector<char> c2(b.begin(), b.begin() + long(k));
        c2.insert(c2.end(), b.begin() + long(k + tl), b.end());
        { std::ofstream o(dir + "/miss.sgy", std::ios::binary); o.write(c2.data(), long(c2.size())); }
        SegyVolume M;
        const bool ok = M.open(dir + "/miss.sgy", e);
        std::vector<float> t(NS, 1.0f);
        const bool got = ok && M.trace(4, 6, t.data());
        double mx = 0; for (float v : t) mx = std::max(mx, double(std::fabs(v)));
        std::printf("    missing trace: %lld missing, trace (4,6) present %d, max |a| %.1e\n", ok ? M.missing : -1LL, int(got), mx);
        check(ok && M.missing == 1 && !got && mx == 0 && M.ni == NI && M.nj == NJ, "C36: missing traces are detected and read as zeros");
    }
    {   // complex trace: a 30 Hz cosine with a Gaussian envelope
        const int n = 500; const double dt = 2.0;
        std::vector<float> x(n), env(n), fr(n), ph(n);
        for (int i = 0; i < n; ++i) { const double t = i * dt / 1000.0, g = std::exp(-std::pow((i - 250) / 120.0, 2)); x[size_t(i)] = float(g * std::cos(2 * 3.14159265358979 * 30 * t)); }
        traceAttribute(SeisAttr::Envelope, x.data(), n, dt, 3, env.data());
        traceAttribute(SeisAttr::Frequency, x.data(), n, dt, 3, fr.data());
        traceAttribute(SeisAttr::Phase, x.data(), n, dt, 3, ph.data());
        double eErr = 0, fErr = 0;
        for (int i = 150; i < 350; ++i) { eErr = std::max(eErr, std::fabs(env[size_t(i)] - std::exp(-std::pow((i - 250) / 120.0, 2)))); fErr = std::max(fErr, std::fabs(fr[size_t(i)] - 30.0)); }
        std::printf("    30 Hz Gaussian-windowed cosine: envelope error %.2e, instantaneous frequency error %.3f Hz, phase(250) %.2f deg\n", eErr, fErr, ph[250]);
        check(eErr < 0.01 && fErr < 0.3, "C36: Hilbert transform - envelope and instantaneous frequency of an analytic test signal");
    }
    if (ok5 && C.ni == NI) {
        // coherence: identical traces -> 1; on the synthetic cube the fault (i - j = NI/4) is a low-coherence lineament
        std::vector<float> a(100); for (int s = 0; s < 100; ++s) a[size_t(s)] = float(std::sin(0.2 * s));
        const double s1 = semblance({a.data(), a.data(), a.data()}, 10, 60);
        const SeismicCube Co = cubeAttribute(C, SeisAttr::Coherence, 4);
        double onF = 0, offF = 0; int nOn = 0, nOff = 0;
        for (int i = 1; i < NI - 1; ++i) for (int j = 1; j < NJ - 1; ++j) {
            const double dd = i - j - NI / 4.0;
            double m = 0; for (int s = 30; s < NS - 30; ++s) m += Co.at(i, j, s); m /= (NS - 60);
            if (std::fabs(dd) <= 1.0) { onF += m; ++nOn; } else if (std::fabs(dd) >= 4.0) { offF += m; ++nOff; }
        }
        onF /= std::max(1, nOn); offF /= std::max(1, nOff);
        std::printf("    semblance of identical traces %.6f; mean coherence on the fault %.3f vs away %.3f\n", s1, onF, offF);
        check(std::fabs(s1 - 1) < 1e-12 && onF < offF - 0.05, "C36: semblance coherence (3 x 3 traces) shows the fault as a low-coherence lineament");
        // streamed attribute volume == in-memory attribute cube
        double dmax = 0;
        bool okv = true;
        for (SeisAttr at : {SeisAttr::Envelope, SeisAttr::Coherence, SeisAttr::Sweetness}) {
            okv = okv && volumeAttribute(V, at, 4, dir + "/attr.sgy", e);
            const SeismicCube A = readSegy(dir + "/attr.sgy", e), B = cubeAttribute(C, at, 4);
            okv = okv && A.a.size() == B.a.size();
            for (size_t k = 0; k < A.a.size() && k < B.a.size(); ++k) dmax = std::max(dmax, double(std::fabs(A.a[k] - B.a[k])));
        }
        std::printf("    streamed attribute volumes (envelope, coherence, sweetness) vs in-memory: max diff %.1e\n", dmax);
        check(okv && dmax < 1e-5, "C36: attribute volumes streamed inline by inline to SEG-Y equal the in-memory attribute cube");
    }
    {   // flattening: traces shifted by an integer number of samples become identical after flattening
        const int ntr = 6, ns = 80;
        std::vector<float> sec(size_t(ntr * ns), 0.0f);
        std::vector<double> hz(ntr);
        for (int t = 0; t < ntr; ++t) { hz[size_t(t)] = 30 + 2 * t; for (int s = 0; s < ns; ++s) sec[size_t(t * ns + s)] = float(std::exp(-std::pow((s - hz[size_t(t)]) / 3.0, 2))); }
        const std::vector<float> f = flattenSection(sec, ntr, ns, hz, 40);
        double d = 0; for (int t = 1; t < ntr; ++t) for (int s = 0; s < 60; ++s) d = std::max(d, double(std::fabs(f[size_t(t * ns + s)] - f[size_t(s)])));
        int pk = 0; for (int s = 0; s < ns; ++s) if (f[size_t(s)] > f[size_t(pk)]) pk = s;
        std::printf("    flattening: max trace difference %.1e, reflector at sample %d (reference 40)\n", d, pk);
        check(d < 1e-6 && pk == 40, "C36: section flattened on a horizon (reflector flat at the reference time)");
    }
    {   // larger file: open / slice timings with a small trace cache (out-of-core behaviour)
        const auto t0 = std::chrono::steady_clock::now();
        const bool ok = writeSyntheticSegy(dir + "/big.sgy", 160, 120, 400, 4.0, 1, false, e);
        const auto t1 = std::chrono::steady_clock::now();
        SegyVolume B; B.cacheTraces = 256;
        const bool op = ok && B.open(dir + "/big.sgy", e);
        const auto t2 = std::chrono::steady_clock::now();
        size_t n = 0;
        if (op) { n += B.inlineSlice(80).size(); n += B.xlineSlice(60).size(); n += B.timeSlice(200).size(); }
        const auto t3 = std::chrono::steady_clock::now();
        auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
        std::printf("    %.1f MB IBM file: write %.0f ms, header index %.0f ms, inline + crossline + time slice %.0f ms (%zu samples)\n",
                    double(B.fileBytes) / 1048576.0, ms(t0, t1), ms(t1, t2), ms(t2, t3), n);
        check(op && n == size_t(120 * 400 + 160 * 400 + 160 * 120), "C36: out-of-core browsing of a 31 MB IBM-float volume (index + three slices)");
        fs::remove(dir + "/big.sgy");
    }
    std::error_code ec; fs::remove_all(dir, ec);
}

// C36 geomodelling round 2: surface editing, fault-parallel sloped pillars, Direct Sampling MPS, uncertainty loop
static void testC36Geo() {
    namespace fs = std::filesystem;
    using namespace ressim::geo;
    std::printf("== C36 geomodel: surface editing, sloped pillars, MPS, uncertainty loop\n");
    // ---- surface editing on a plane z = 1000 + 0.01 x (rotated lattice)
    {
        Surface s; s.name = "P"; s.nx = 81; s.ny = 61; s.dx = 25; s.dy = 25; s.x0 = 5000; s.y0 = 8000; s.rotDeg = 30;
        s.z.resize(size_t(s.nx) * size_t(s.ny));
        for (int j = 0; j < s.ny; ++j) for (int i = 0; i < s.nx; ++i) { double x, y; s.nodeXY(i, j, x, y); s.z[size_t(i + s.nx * j)] = 1000 + 0.01 * (x - 5000); }
        const Surface s0 = s;
        double cx, cy; s.nodeXY(40, 30, cx, cy);
        const int n = brushSurface(s, cx, cy, 200, 15, 0);
        double maxOut = 0;
        for (int j = 0; j < s.ny; ++j) for (int i = 0; i < s.nx; ++i) {
            double x, y; s.nodeXY(i, j, x, y);
            if (std::hypot(x - cx, y - cy) > 200) maxOut = std::max(maxOut, std::fabs(s.node(i, j) - s0.node(i, j)));
        }
        const double dC = s.node(40, 30) - s0.node(40, 30), dR = s.node(44, 30) - s0.node(44, 30);   // r = 100: w = exp(-0.75)
        Surface f = s; brushSurface(f, cx, cy, 150, 1012.5, 2);
        double fi, fj; surfaceIJ(s, cx, cy, fi, fj);
        std::printf("    brush: %d nodes, centre +%.3f, r=100 +%.3f (expected %.3f), outside max %.1e; flatten centre %.3f; node (%.3f, %.3f)\n", n, dC, dR, 15 * std::exp(-0.75), maxOut, f.node(40, 30), fi, fj);
        check(n > 100 && std::fabs(dC - 15) < 1e-9 && std::fabs(dR - 15 * std::exp(-0.75)) < 1e-9 && maxOut == 0 && std::fabs(f.node(40, 30) - 1012.5) < 1e-9 &&
              std::fabs(fi - 40) < 1e-9 && std::fabs(fj - 30) < 1e-9, "C36: surface brush - Gaussian taper exp(-3 r^2/R^2), raise / flatten, rotation-aware node indices, nothing outside R");
        // roughness reduction + plane preserved by smoothing; holes in-painted on the plane
        Surface r = s0; std::mt19937 rng(5); std::normal_distribution<double> N(0, 2);
        for (double& z : r.z) z += N(rng);
        auto rough = [](const Surface& a) { double q = 0; int m = 0; for (int j = 1; j + 1 < a.ny; ++j) for (int i = 1; i + 1 < a.nx; ++i) { const double l = a.node(i + 1, j) + a.node(i - 1, j) + a.node(i, j + 1) + a.node(i, j - 1) - 4 * a.node(i, j); q += l * l; ++m; } return std::sqrt(q / m); };
        const double r0 = rough(r); smoothSurface(r, 3); const double r1 = rough(r);
        Surface p = s0; smoothSurface(p, 2);
        double dPlane = 0; for (int j = 2; j + 2 < p.ny; ++j) for (int i = 2; i + 2 < p.nx; ++i) dPlane = std::max(dPlane, std::fabs(p.node(i, j) - s0.node(i, j)));
        Surface h = s0; for (int j = 20; j < 35; ++j) for (int i = 30; i < 50; ++i) h.z[size_t(i + h.nx * j)] = std::nan("");
        const int nh = fillSurfaceHoles(h);
        double dHole = 0; for (int j = 20; j < 35; ++j) for (int i = 30; i < 50; ++i) dHole = std::max(dHole, std::fabs(h.node(i, j) - s0.node(i, j)));
        std::vector<Surface> two = {s0, isochoreSurface(s0, -5, "X")};
        const int moved = enforceThickness(two, 2.0);
        double minT = 1e9; for (size_t q = 0; q < two[1].z.size(); ++q) minT = std::min(minT, two[1].z[q] - two[0].z[q]);
        std::printf("    smoothing roughness %.3f -> %.3f, plane change %.1e; %d holes filled, max error %.2e; crossing horizon: %d nodes moved, min thickness %.3f\n", r0, r1, dPlane, nh, dHole, moved, minT);
        check(r1 < 0.2 * r0 && dPlane < 1e-9 && nh == 300 && dHole < 1e-3 && moved == int(two[1].z.size()) && std::fabs(minT - 2.0) < 1e-9,
              "C36: smoothing keeps planes, Laplace hole filling is exact on a plane, enforced minimum thickness removes crossings");
    }
    // ---- sloped pillars on the demo project (no well tie: the horizons are honoured exactly at the pillar lines)
    const std::string dir = (fs::temp_directory_path() / "ressim_c36_geo").string();
    std::string e;
    writeDemoProject(dir, false, e);
    std::vector<Surface> hz;
    for (const char* n : {"TopRes", "MidRes", "BaseRes"}) { Surface s = readIrap(dir + "/" + n + ".irap", e); s.name = n; hz.push_back(s); }
    const auto faults = readFaultSticks(dir + "/F1_sticks.txt", e);
    auto wells = readWellHeads(dir + "/well_heads.txt", e);
    for (auto& w : wells) {
        if (fs::exists(dir + "/" + w.name + "_survey.txt")) readDeviation(dir + "/" + w.name + "_survey.txt", w, e);
        w.logs.push_back(readLas(dir + "/" + w.name + ".las", e));
    }
    for (const auto& t : readTops(dir + "/well_tops.txt", e)) for (auto& w : wells) if (w.name == t.first) w.tops.push_back(t.second);
    GridSpec gs; gs.dx = gs.dy = 100; gs.nx = 40; gs.ny = 30; gs.x0 = hz[0].x0; gs.y0 = hz[0].y0; gs.rotDeg = hz[0].rotDeg;
    Zone z1; z1.name = "Upper"; z1.top = 0; z1.base = 1; z1.nLayers = 8;
    Zone z2; z2.name = "Lower"; z2.top = 1; z2.base = 2; z2.nLayers = 8;
    gs.zones = {z1, z2}; gs.tieRange = 1500; gs.tieToWells = false;
    std::string e0, e1;
    const GeoGrid gv = buildGrid(gs, hz, faults, wells, e0);
    GridSpec gp = gs; gp.pillars = 1;
    const GeoGrid gq = buildGrid(gp, hz, faults, wells, e1);
    {
        // fault dip from the sticks (plan offset per unit depth)
        const auto a = faults[0].traceAt(2100), b = faults[0].traceAt(2200);
        double fx = 0, fy = 0; for (size_t k = 0; k < a.size(); ++k) { fx += (b[k][0] - a[k][0]) / 100; fy += (b[k][1] - a[k][1]) / 100; }
        fx /= double(a.size()); fy /= double(a.size());
        const double fSlope = std::hypot(fx, fy);
        const int nx = gq.grid.nx, ny = gq.grid.ny, nz = gq.grid.nz;
        // slope of the pillar nearest the fault trace / far away
        double nearS = 0, farS = 1e9, nearD = 1e300, farD = 0;
        const auto tr = faults[0].traceAt(2150);
        for (int jj = 0; jj <= ny; ++jj) for (int ii = 0; ii <= nx; ++ii) {
            const double* p = &gq.grid.coord[size_t(6 * (ii + (nx + 1) * jj))];
            const double sl = std::hypot(p[3] - p[0], p[4] - p[1]) / (p[5] - p[2]);
            double d = 1e300; for (const auto& q : tr) d = std::min(d, std::hypot(q[0] - p[0], q[1] - p[1]));
            if (d < nearD) { nearD = d; nearS = sl; }
            if (d > farD) { farD = d; farS = sl; }
        }
        // conformity: top horizon depth at the sloped pillar position (nodes away from the fault, single block)
        double conf = 0; int nconf = 0;
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
            const double z = gq.grid.zcorn[size_t((0 * 2 * ny + 2 * j) * 2 * nx + 2 * i)];
            double x, y; gq.pillarXY(i, j, z, x, y);
            double d = 1e300; for (const auto& q : faults[0].traceAt(z)) d = std::min(d, std::hypot(q[0] - x, q[1] - y));
            if (d < 500) continue;   // the fault band (horizons re-extrapolated from each side) is 4 cells wide
            conf = std::max(conf, std::fabs(z - hz[0].at(x, y))); ++nconf;
        }
        double bv = 0, bq = 0; int av = 0, aq = 0;
        for (size_t c = 0; c < gv.bulk.size(); ++c) { if (gv.grid.actnum[c]) { bv += gv.bulk[c]; ++av; } if (gq.grid.actnum[c]) { bq += gq.bulk[c]; ++aq; } }
        int same = 0, tot = 0;
        for (int k = 0; k < nz; k += 3) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
            const size_t c = size_t(i + nx * (j + ny * k));
            if (!gq.grid.actnum[c] || gq.bulk[c] <= 0) continue;
            ++tot; same += gq.cellOf(gq.center[c][0], gq.center[c][1], gq.center[c][2]) == int(c);
        }
        std::printf("    fault dip slope %.4f; pillar at %.0f m from the trace slope %.4f, farthest (%.0f m) %.2e; max pillar slope %.4f (vertical grid %.1e)\n", fSlope, nearD, nearS, farD, farS, gq.maxPillarSlope(), gv.maxPillarSlope());
        std::printf("    top horizon at the sloped pillar lines: max error %.3f (%d nodes); bulk vertical %.6g sloped %.6g (%+.2f %%), active %d / %d; cellOf(centre) %d / %d\n",
                    conf, nconf, bv, bq, 100 * (bq / bv - 1), av, aq, same, tot);
        check(e1.empty() && fSlope > 0.05 && std::fabs(nearS - fSlope) < 0.25 * fSlope && farS < 1e-3 * fSlope && gv.maxPillarSlope() == 0 && conf < 0.05 && nconf > 300,
              "C36: fault-parallel pillars - dip of the fault plane at the fault, vertical far away, horizons met at the sloped pillar lines");
        check(std::fabs(bq / bv - 1) < 0.03 && same >= int(0.99 * tot), "C36: sloped-pillar grid - bulk volume within 3 % of the vertical-pillar grid, cellOf inverts the sloped geometry");
    }
    // ---- MPS: channel TI -> 2-D grid, conditioned, anisotropic continuity reproduced (SIS for contrast)
    {
        const TrainingImage ti = channelTrainingImage(200, 200, 1, 6, 80, 12, 0.28, 3);
        Surface t0; t0.nx = 3; t0.ny = 3; t0.dx = 3000; t0.dy = 3000; t0.z.assign(9, 1000.0); t0.name = "T";
        Surface t1 = isochoreSurface(t0, 10, "B");
        GridSpec g2; g2.nx = 80; g2.ny = 80; g2.dx = g2.dy = 25; Zone z; z.top = 0; z.base = 1; z.nLayers = 1; g2.zones = {z}; g2.tieToWells = false;
        std::string e2;
        const GeoGrid G = buildGrid(g2, {t0, t1}, {}, {}, e2);
        std::vector<Blocked> data;
        std::mt19937 rng(9);
        for (int q = 0; q < 60; ++q) {   // "truth" = TI window [50,130)^2 (grid i -> TI x for azimuth 90)
            Blocked b; b.cell = int(rng() % 6400u); b.v["FACIES"] = ti.at(50 + b.cell % 80, 50 + b.cell / 80, 0); data.push_back(b);
        }
        FaciesSpec fsp; fsp.codes = {0, 1}; fsp.curve = "FACIES";
        MpsSpec m; m.azimuthDeg = 90; m.tiCell = 25; m.seed = 21;   // TI +x -> map east = grid i
        const auto t1c = std::chrono::steady_clock::now();
        const FaciesResult R = modelFaciesMps(G, data, fsp, ti, m);
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1c).count();
        double tiFrac = 0; for (int v : ti.v) tiFrac += v; tiFrac /= double(ti.v.size());
        // continuity: unconditional MPS vs unconditional SIS with the same proportions
        MpsSpec mu = m; mu.seed = 22;
        const FaciesResult U = modelFaciesMps(G, {}, fsp, ti, mu);
        std::vector<int> cat(U.facies.begin(), U.facies.end());
        const double alongM = indicatorContinuity(cat, 80, 80, 1, 12, 0, 1), acrossM = indicatorContinuity(cat, 80, 80, 1, 0, 12, 1);
        const double alongT = indicatorContinuity(ti.v, 200, 200, 1, 12, 0, 1), acrossT = indicatorContinuity(ti.v, 200, 200, 1, 0, 12, 1);
        FaciesSpec sis = fsp; sis.method = 0; sis.vg.rangeMax = sis.vg.rangeMin = 200; sis.vg.rangeVert = 10; sis.proportions = {1 - tiFrac, tiFrac};
        const FaciesResult Rs = modelFacies(G, {}, sis);
        std::vector<int> cs(Rs.facies.begin(), Rs.facies.end());
        const double alongS = indicatorContinuity(cs, 80, 80, 1, 12, 0, 1), acrossS = indicatorContinuity(cs, 80, 80, 1, 0, 12, 1);
        const std::string tiPath = dir + "/ti.gslib";
        std::string e3; writeTrainingImage(tiPath, ti); const TrainingImage back = readTrainingImage(tiPath, e3);
        std::printf("    MPS %.2f s: conditional channel fraction %.3f, unconditional %.3f (TI %.3f, servosystem); %d conditioned, %d mismatches; TI file round trip %s\n",
                    sec, R.proportions.size() > 1 ? R.proportions[1] : -1.0, U.proportions.size() > 1 ? U.proportions[1] : -1.0, tiFrac, R.conditioned, R.mismatches, back.v == ti.v ? "exact" : "DIFFERS");
        std::printf("    P(I(x+12)|I(x)) along / across the channels: MPS %.2f / %.2f, TI %.2f / %.2f, SIS %.2f / %.2f\n", alongM, acrossM, alongT, acrossT, alongS, acrossS);
        check(R.mismatches == 0 && R.conditioned > 50 && R.proportions.size() == 2 && std::fabs(R.proportions[1] - tiFrac) < 0.05 && U.proportions.size() == 2 &&
              std::fabs(U.proportions[1] - tiFrac) < 0.05 && back.v == ti.v,
              "C36: Direct Sampling MPS honours the well data and the training-image proportions (servosystem); GSLIB training image round trip");
        check(alongM > 1.5 * acrossM && alongT > 1.5 * acrossT && alongM / acrossM > 1.3 * (alongS / std::max(1e-9, acrossS)),
              "C36: MPS reproduces the channel continuity of the training image (along >> across), unlike two-point SIS");
    }
    // ---- uncertainty loop on the demo project
    {
        UncInputs in;
        in.gs = gs; in.gs.tieToWells = true; in.hz = hz; in.faults = faults; in.wells = wells;
        in.fs.method = 0; in.fs.vg.rangeMax = in.fs.vg.rangeMin = 800; in.fs.vg.rangeVert = 6; in.fs.seed = 3;
        in.ps.vg.rangeMax = in.ps.vg.rangeMin = 1000; in.ps.vg.rangeVert = 8; in.ps.maxV = 0.4; in.ps.seed = 11;
        in.vs.contacts.owc = 2080; in.vs.satHeight = false; in.vs.swConst = 0.25; in.vs.bo = 1.25;
        UncSpec u; u.nReal = 12; u.zSd = 8; u.zRange = 1500; u.owc = Tri{2060, 2080, 2095}; u.poroMult = Tri{0.9, 1.0, 1.1}; u.tornadoStochastic = 4;
        const auto t0 = std::chrono::steady_clock::now();
        const UncResult R = uncertaintyLoop(in, u);
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("    %s", R.report.c_str());
        std::printf("    %.1f s; P50 realisation #%d\n", sec, R.iP50);
        bool owcBar = false; for (const auto& t : R.tornado) if (t.name == "OWC depth") owcBar = t.low < R.base && R.base < t.high;
        // deterministic rebuild of the P50 realisation
        std::string er;
        double again = 0;
        if (R.iP50 >= 0) again = realise(in, R.real[size_t(R.iP50)].s, er).vol.total.stoiip;
        check(R.err.empty() && R.stoiip.p90 < R.stoiip.p50 && R.stoiip.p50 < R.stoiip.p10 && R.tornado.size() >= 4 && owcBar && R.tornado[0].high - R.tornado[0].low >= R.tornado.back().high - R.tornado.back().low,
              "C36: uncertainty loop - P90 < P50 < P10, tornado sorted by swing with the OWC bar around the base case");
        {   // Latin hypercube: the triangular CDF of the OWC samples puts exactly one sample in each of the nReal strata
            std::vector<int> hit(size_t(u.nReal), 0);
            const double A = u.owc.lo, B = u.owc.hi, C = u.owc.mode;
            for (const auto& x : R.real) {
                const double o = x.s.owc, F = o < C ? (o - A) * (o - A) / ((B - A) * (C - A)) : 1 - (B - o) * (B - o) / ((B - A) * (B - C));
                ++hit[size_t(std::clamp(int(F * u.nReal), 0, u.nReal - 1))];
            }
            bool lhs = true; for (int h : hit) lhs = lhs && h == 1;
            check(lhs, "C36: uncertainty loop - Latin hypercube sampling, one OWC sample per probability stratum");
        }
        check(R.iP50 >= 0 && er.empty() && std::fabs(again / R.real[size_t(R.iP50)].stoiip - 1) < 1e-12, "C36: a realisation (structure field, contacts, seeds) is rebuilt exactly from its sample");
        // structural field: zero at the well tops, ~ sd elsewhere
        std::vector<std::array<double, 2>> pts; for (const auto& w : wells) for (const auto& t : w.tops) { const P3 q = w.at(t.md); pts.push_back({q[0], q[1]}); }
        const auto P = perturbHorizons(hz, 10, 1500, 77, pts);
        double atW = 0; for (const auto& p : pts) atW = std::max(atW, std::fabs(P[0].at(p[0], p[1]) - hz[0].at(p[0], p[1])));
        double s2 = 0; int ns = 0;
        for (int j = 0; j < hz[0].ny; ++j) for (int i = 0; i < hz[0].nx; ++i) {
            double x, y; hz[0].nodeXY(i, j, x, y);
            double d = 1e300; for (const auto& p : pts) d = std::min(d, std::hypot(p[0] - x, p[1] - y));
            if (d < 1500) continue;
            const double v = P[0].node(i, j) - hz[0].node(i, j); s2 += v * v; ++ns;
        }
        const double sdFar = ns ? std::sqrt(s2 / ns) : 0;
        std::printf("    structural field: max |dz| at the %zu well tops %.2e (bilinear), sd far from wells %.2f (target 10, %d nodes)\n", pts.size(), atW, sdFar, ns);
        {   // MPS facies inside the loop (fs.method 3 + training image)
            UncInputs im = in; im.fs.method = 3; im.fs.codes = {0, 1, 2};
            auto T = std::make_shared<TrainingImage>(channelTrainingImage(150, 150, 1, 5, 60, 10, 0.3, 4));
            im.ti = T; im.mps.tiCell = 100; im.mps.azimuthDeg = 30;
            UncSpec um; um.nReal = 3; um.tornadoStochastic = 0; um.owc = Tri{2070, 2080, 2090};
            const UncResult Rm = uncertaintyLoop(im, um);
            UncInputs bad = im; bad.ti.reset();
            const UncResult Rb = uncertaintyLoop(bad, um);
            std::printf("    MPS facies in the loop: err '%s', STOIIP P50 %.4g; without a training image: '%s'\n", Rm.err.c_str(), Rm.stoiip.p50, Rb.err.c_str());
            check(Rm.err.empty() && Rm.stoiip.p50 > 0 && !Rb.err.empty(), "C36: uncertainty loop with MPS facies (training image shared by all realisations)");
        }
        check(atW < 0.3 && sdFar > 4 && sdFar < 16, "C36: structural random field is tied to the wells (zero at the tops) with the requested sd away from them");
    }
}

// ============================================================================ C36-5: MSW labyrinth / acceleration / drift flux
namespace c36m {
std::string spe1Vertical(const std::string& flowModel) {   // SPE1 producer with its full tubing to surface as a multisegment well
    std::ifstream f(findExample("SPE1.DATA"));
    std::stringstream ss; ss << f.rdbuf();
    std::string t = ss.str();
    size_t p = t.find(" 'PROD' 'PROD' 10 10 8400 OIL /");
    if (p == std::string::npos) return "";
    t.replace(p, std::string(" 'PROD' 'PROD' 10 10 8400 OIL /").size(), " 'PROD' 'PROD' 10 10 0 OIL /");
    p = t.find("WCONPROD");
    if (p == std::string::npos) return "";
    t.insert(p, "WELSEGS\n 'PROD' 0 0 1* ABS 'HFA' '" + flowModel + "' /\n 2 17 1 1 8400 8400 0.5 1E-4 /\n/\n\nCOMPSEGS\n 'PROD' /\n 10 10 3 1 8375 8425 /\n/\n\n");
    return t;
}
}  // namespace c36m

static void testC36Msw() {
    std::printf("== C36 multisegment wells: WSEGLABY, acceleration (HFA), drift flux (DF)\n");
    using namespace ressim::msw;
    const double PI = 3.14159265358979323846;
    ressim::msw::Units u; u.g = 9.80665e-5; u.toM = 1; u.toRho = 1; u.toQ = 1.0 / 86400.0; u.toPa = 1e5;
    // ---- parser
    {
        std::vector<ExtraKeyword> kw(3);
        kw[0].name = "WELSEGS"; kw[0].body = " 'W' 1000 0 1* INC 'HFA' 'DF' /\n 2 3 1 1 100 50 0.1 1E-4 /\n 4 4 2 3 0.1 0 0.1 1E-4 /\n/\n";
        kw[1].name = "WSEGLABY"; kw[1].body = " 'W' 4 4 2 3 /\n/\n";
        kw[2].name = "WSEGDFPA"; kw[2].body = " 1.1 0.15 0.35 1.25 0.25 1.0 /\n";
        std::string err;
        const auto W = parse(kw, 0, true, err);
        bool ok = W.size() == 1 && err.empty();
        if (ok) {
            const Segment& g = W[0].segs[size_t(W[0].indexOf(4))];
            ok = W[0].accelOn && W[0].driftFlux && g.device == 4 && g.labN == 2 && g.labCh == 3 &&
                 std::fabs(g.labA - 6e-5 / (0.3048 * 0.3048)) < 1e-12 && std::fabs(g.labD - 0.006316 / 0.3048) < 1e-12 &&
                 std::fabs(W[0].df.A - 1.25) < 1e-12 && std::fabs(W[0].df.a2 - 0.35) < 1e-12 && std::fabs(W[0].df.vdFactor - 1.1) < 1e-12;
        }
        check(ok, "parser: WELSEGS 'HFA' 'DF', WSEGLABY defaults (metric, converted for FIELD), WSEGDFPA parameters");
    }
    // ---- labyrinth device: hand formula and the device balance
    {
        Segment g; g.device = 4; g.labN = 2; g.labCh = 3; g.labA = 6e-5; g.labL1 = 0.1863; g.labL2 = 0.2832; g.labD = 0.006316; g.labR = 1e-5;
        const double q = 3e-4, rho = 900, mu = 2e-3;
        const double v = q / (2 * 6e-5), Re = rho * v * 0.006316 / mu;
        const double t1 = -1.8 * std::log10(std::pow(1e-5 / 0.006316 / 3.7, 1.11) + 6.9 / Re), f = 1 / (t1 * t1);
        const double hand = 0.5 * rho * v * v * (0.35 + 0.5 + 2 * 0.5 + f * (0.1863 + 2 * 0.2832) / 0.006316);
        const double dp = labyrinthDp(g, q, rho, mu, 1.0);
        // in a well: tubing node 2 at the top depth, device segment 3 feeding it
        Well w; w.name = "L";
        Segment t; t.id = 1; t.depth = 1000; w.segs.push_back(t);
        Segment s; s.id = 2; s.outlet = 1; s.length = 100; s.depth = 1000; w.segs.push_back(s);
        Segment d = g; d.id = 3; d.outlet = 2; d.branch = 2; d.length = 100.1; d.depth = 1000; w.segs.push_back(d);
        Inflow in; in.resize(3); in.q[2] = 30; in.rho[2] = 900; in.mu[2] = 2; in.a[2] = 5; in.pres[2] = 120;
        const Solution S = solve(w, 100.0, in, u);
        const double pin = 100 + S.dFric[1] + S.dFric[2] + S.dHydro[2];
        const double qd = 5 * (120 - pin - S.dDev[2]);
        const double dq = labyrinthDp(g, qd / 86400.0, 900, 2e-3, 1.0) / 1e5;
        std::printf("    labyrinth dp %.6g Pa (hand %.6g); in-well device drop %.4f bar = dp(q %.3f rm3/d) %.4f bar\n", dp, hand, S.dDev[2], qd, dq);
        check(std::fabs(dp - hand) < 1e-9 * hand && S.dDev[2] > 0 && std::fabs(S.dDev[2] - dq) < 1e-6 * S.dDev[2],
              "WSEGLABY: rho v^2 / 2 (inlet + outlet + bends + Haaland L / D), balanced with the subtree inflow");
    }
    // ---- acceleration: momentum flux balance
    {
        Well w; w.name = "A"; w.accelOn = true;
        Segment t; t.id = 1; t.depth = 1000; t.diam = 0.1; w.segs.push_back(t);
        Segment s2; s2.id = 2; s2.outlet = 1; s2.length = 100; s2.depth = 1000; s2.diam = 0.08; w.segs.push_back(s2);
        Segment s3; s3.id = 3; s3.outlet = 2; s3.length = 200; s3.depth = 1000; s3.diam = 0.1; w.segs.push_back(s3);
        Inflow in; in.resize(3); in.q[2] = 2000; in.rho[2] = 800; in.mu[2] = 1;
        const Solution S = solve(w, 100.0, in, u);
        const double q = 2000 / 86400.0, A2 = 0.25 * PI * 0.0064, A3 = 0.25 * PI * 0.01;
        const double h3 = 800 * q * q / (A3 * A3) / 1e5, h2 = 800 * q * q * (1 / A2 - 1 / A3) / A2 / 1e5;
        Well w0 = w; w0.accelOn = false;
        const Solution S0 = solve(w0, 100.0, in, u);
        std::printf("    acceleration: seg3 %.6f bar (hand %.6f), seg2 (contraction) %.6f bar (hand %.6f); p toe %.5f vs HF- %.5f\n", S.dAcc[2], h3, S.dAcc[1], h2, S.p[2], S0.p[2]);
        check(std::fabs(S.dAcc[2] - h3) < 1e-9 * h3 && std::fabs(S.dAcc[1] - h2) < 1e-9 * h2 && std::fabs(S.p[2] - S0.p[2] - h2 - h3) < 1e-9,
              "HFA acceleration: (M_seg - sum M_child) / A_seg, M = rho q^2 / A (perforation inflow without axial momentum)");
    }
    // ---- drift-flux holdup: limits
    {
        DriftParams P;
        const double rl = 800, rg = 100, D = 0.1;
        const double Vc = std::pow(P.sigma * 9.80665 * (rl - rg) / (rl * rl), 0.25);
        double C0 = 0, vd = 0;
        // bubbly upward: alpha = jg / (C0 j + v_d), v_d -> 1.53 V_c (1 - alpha C0) ... at small alpha
        const double a1 = driftHoldup(0.005, 0.5, rg, rl, D, 0.0, P, &C0, &vd);
        const double lam1 = 0.005 / 0.505;
        const bool b1 = std::fabs(a1 * (C0 * 0.505 + vd) - 0.005) < 1e-9 && a1 < lam1 && std::fabs(C0 - 1.2) < 1e-9 && vd > 0.9 * 1.53 * Vc && vd <= 1.53 * Vc;
        std::printf("    bubbly up: alpha %.5f (no-slip %.5f), C0 %.3f, v_d %.4f m/s (1.53 V_c %.4f)\n", a1, lam1, C0, vd, 1.53 * Vc);
        // horizontal: no drift, alpha = jg / (C0 j)
        const double a2 = driftHoldup(0.3, 0.3, rg, rl, D, PI / 2, P, &C0, &vd);
        const bool b2 = std::fabs(vd) < 1e-12 && std::fabs(a2 * C0 * 0.6 - 0.3) < 1e-9;
        // downward co-current: gas slower -> alpha above the no-slip fraction
        const double a3 = driftHoldup(0.1, 0.5, rg, rl, D, PI, P, &C0, &vd);
        const bool b3 = a3 > 0.1 / 0.6 && vd < 0;
        // high mixture velocity: C0 -> 1, v_d small -> alpha close to the no-slip fraction
        const double a4 = driftHoldup(20, 2, rg, rl, D, 0.0, P, &C0, &vd);
        const bool b4 = std::fabs(a4 - 20.0 / 22.0) < 0.02 && C0 < 1.05;
        std::printf("    horizontal alpha %.4f (C0 %.3f), downward %.4f (no-slip %.4f), annular-like %.4f (no-slip %.4f, C0 %.3f)\n", a2, C0, a3, 0.1 / 0.6, a4, 20.0 / 22.0, C0);
        check(b1 && b2 && b3 && b4, "drift flux (Shi et al. 2005): bubble-rise limit, no drift horizontally, gas lags in downflow, C0 -> 1 at high velocity");
    }
    // ---- drift flux in the segment solve: heavier column than homogeneous flow
    {
        Well w; w.name = "V"; w.driftFlux = true;
        Segment t; t.id = 1; t.depth = 0; t.diam = 0.1; w.segs.push_back(t);
        Segment s; s.id = 2; s.outlet = 1; s.length = 1000; s.depth = 1000; s.diam = 0.1; w.segs.push_back(s);
        Inflow in; in.resize(2); in.q[1] = 400; in.qg[1] = 100; in.rho[1] = 0.75 * 800 + 0.25 * 100; in.mu[1] = 1; in.rhoG[1] = 100; in.rhoL[1] = 800;
        const Solution S = solve(w, 50.0, in, u);
        Well wh = w; wh.driftFlux = false;
        const Solution SH = solve(wh, 50.0, in, u);
        std::printf("    vertical 1000 m: holdup %.4f vs no-slip %.4f, head density %.1f vs %.1f kg/m3, head %.3f vs %.3f bar\n", S.holdup[1], S.noSlip[1], S.rhoHead[1], SH.rhoHead[1], S.dHydro[1], SH.dHydro[1]);
        check(S.holdup[1] < S.noSlip[1] && std::fabs(SH.holdup[1] - 0.25) < 1e-12 && S.dHydro[1] > SH.dHydro[1] * 1.01 &&
              std::fabs(S.rhoHead[1] - (S.holdup[1] * 100 + (1 - S.holdup[1]) * 800)) < 1e-9,
              "DF segment solve: gas slip lowers the in-situ gas fraction and raises the hydrostatic head of an upward column");
    }
    // ---- deck runs: labyrinth ICDs on MSW_HORIZONTAL, HFA on/off
    {
        const std::string path = findExample("MSW_HORIZONTAL.DATA");
        if (!path.empty()) {
            std::vector<DeckMessage> rep;
            const Case base = importEclipseDeck(path, nullptr, &rep);
            Case lab = mswVariant(base, 1), acc = base;
            ExtraKeyword L; L.name = "WSEGLABY"; L.section = DeckSection::Schedule; L.step = 0; L.body = " 'P1' 20 37 1 2 /\n/\n";
            lab.extraKeywords.push_back(L);
            for (auto& e : acc.extraKeywords) if (e.name == "WELSEGS") { size_t p = e.body.find("'HF-'"); if (p != std::string::npos) e.body.replace(p, 5, "'HFA'"); }
            const MswRun R1 = runMsw(mswVariant(base, 1)), RL = runMsw(lab), R2 = runMsw(base), RA = runMsw(acc);
            double maxAcc = 0, maxDevL = 0;
            for (const auto& rc : RA.rec) if (rc.t > 2) for (double d : rc.dAcc) maxAcc = std::max(maxAcc, std::fabs(d));   // (day 1: potential rates)
            for (const auto& rc : RL.rec) for (double d : rc.dDev) maxDevL = std::max(maxDevL, std::fabs(d));
            std::printf("    tubing only FOPT %.6g FWPT %.6g CV %.3f | labyrinth FOPT %.6g FWPT %.6g CV %.3f max device %.1f psi | SICD FOPT %.6g\n",
                        R1.fopt, R1.fwpt, cv(R1.devQ), RL.fopt, RL.fwpt, cv(RL.devQ), maxDevL, R2.fopt);
            std::printf("    HFA: FOPT %.6g vs HF- %.6g, max acceleration %.4f psi\n", RA.fopt, R2.fopt, maxAcc);
            check(RL.devQ.size() == 18 && maxDevL > 1 && cv(RL.devQ) < 0.8 * cv(R1.devQ) && RL.fwpt < R1.fwpt,
                  "WSEGLABY on MSW_HORIZONTAL: labyrinth ICDs equalise the inflow and cut the water production vs the open tubing");
            check(maxAcc > 0.1 && maxAcc < 3.0 && std::fabs(RA.fopt - R2.fopt) < 2e-3 * R2.fopt,
                  "HFA on MSW_HORIZONTAL: acceleration drops of the order of rho v^2 (< 3 psi for ~7000 rb/d in 2.441 in tubing), FOPT within 0.2 %");
            const std::string txt = eclipseDeckToString(lab);
            check(txt.find("WSEGLABY") != std::string::npos && txt.find(" 'P1' 20 37 1 2 /") != std::string::npos, "WSEGLABY kept verbatim on export");
        }
    }
    // ---- deck runs: SPE1 producer with the full tubing to surface, HO vs DF
    {
        const std::string tHO = c36m::spe1Vertical("HO"), tDF = c36m::spe1Vertical("DF");
        if (!tHO.empty()) {
            auto run = [](const std::string& t, double& fopt, double& fgpt, double& hold, double& ns, double& pBot) {
                std::vector<DeckMessage> rep;
                Case c = eclipseDeckFromString(t, "", nullptr, &rep);
                c.mswOpt.flash = false;   // C36 comparison without the C37 tubing flash (gas stays dissolved; flash: -c37msw)
                c.mswOpt.storage = false; // ... and without segment storage (the recorded holdup is the steady one, HO = no-slip)
                ressim::Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
                while (s.runReportStep()) {}
                fopt = s.summary().back().fopt; fgpt = s.summary().back().fgpt;
                hold = ns = 0;
                for (const auto& r : s.mswRecords()) if (r.holdup.size() > 8 && r.noSlip[8] > ns) { ns = r.noSlip[8]; hold = r.holdup[8]; }
                pBot = s.mswRecords().empty() ? 0.0 : s.mswRecords().back().p.back();
                return s.mswRecords().size();
            };
            double oH, gH, hH, nH, oD, gD, hD, nD, pH, pD;
            const size_t nh = run(tHO, oH, gH, hH, nH, pH), nd = run(tDF, oD, gD, hD, nD, pD);
            std::printf("    SPE1 vertical MSW to surface: HO FOPT %.6g FGPT %.6g | DF FOPT %.6g FGPT %.6g; mid-tubing gas fraction at peak: HO %.3f, DF holdup %.3f (no-slip %.3f)\n",
                        oH, gH, oD, gD, hH, hD, nD);
            std::printf("    final bottom-node pressure HO %.1f psi, DF %.1f psi (top node at the 1000 psi BHP limit)\n", pH, pD);
            check(nh > 0 && nd > 0 && nD > 0.05 && hD < nD && std::fabs(hH - nH) < 1e-12 && pD > pH && oD < oH,
                  "SPE1 with tubing to surface: drift flux keeps more liquid in the column (holdup < no-slip) -> heavier head, less oil than homogeneous flow");
        }
    }
}

// ============================================================================ C36-6: NETBALAN, NETCOMPA / COMPOFF
namespace c36n {
std::string networkDeck(const std::string& extraSchedule, bool compTable) {
    std::ifstream f(findExample("NETWORK.DATA"));
    std::stringstream ss; ss << f.rdbuf();
    std::string t = ss.str();
    if (compTable) {   // VFPPROD 3: the pipeline of table 2 behind a compressor - suction pressure 50 psi + the same line loss
        const size_t a = t.find("VFPPROD\n 2 0"), b = t.find("TUNING", a);
        if (a == std::string::npos || b == std::string::npos) return "";
        std::istringstream is(t.substr(a, b - a));
        std::string line, out;
        int n = 0;
        const double pout[4] = {100, 300, 600, 1000};
        while (std::getline(is, line)) {
            ++n;
            if (n == 2) { out += " 3 0 'OIL' 'WCT' 'GOR' 'THP' 1* 'FIELD' 'BHP' /\n"; continue; }
            std::istringstream ls(line);
            std::vector<std::string> tok; std::string x;
            while (ls >> x) tok.push_back(x);
            if (n >= 8 && tok.size() == 9) {
                const int it = std::atoi(tok[0].c_str());
                char b2[200];
                std::snprintf(b2, sizeof b2, " %s %s %s %s", tok[0].c_str(), tok[1].c_str(), tok[2].c_str(), tok[3].c_str());
                out += b2;
                for (int k = 4; k < 8; ++k) { std::snprintf(b2, sizeof b2, " %.2f", 50 + std::atof(tok[size_t(k)].c_str()) - pout[it - 1]); out += b2; }
                out += " /\n";
            } else out += line + "\n";
        }
        t.insert(b, out);
    }
    const size_t p = t.find("TSTEP\n 40*91.3125");
    if (p == std::string::npos) return "";
    t.insert(p, "GCONPROD\n 'PLAT' ORAT 6000 /\n/\n\n" + extraSchedule);
    return t;
}
struct NetRun { double fopt = 0, fgcrMax = 0; std::vector<double> fgcr, t, pPlat; std::string log; std::string blob; bool ok = false; };
NetRun run(const std::string& deck, int stopAfter = -1, const std::string& blob = "") {
    NetRun r;
    if (deck.empty()) return r;
    std::vector<DeckMessage> rep;
    const Case c = eclipseDeckFromString(deck, "", nullptr, &rep);
    ressim::Simulator s; s.log = [&](const std::string& l) { r.log += l + "\n"; }; s.setup(c); s.initialize();
    if (!blob.empty()) s.loadState(blob);
    while ((stopAfter < 0 || s.reportStep() < stopAfter) && s.runReportStep()) {}
    if (stopAfter >= 0) { r.blob = s.saveState(); r.ok = true; return r; }
    for (const auto& p : s.summary()) { r.fgcr.push_back(p.fgcr); r.t.push_back(p.time); r.fgcrMax = std::max(r.fgcrMax, p.fgcr); }
    r.fopt = s.summary().back().fopt;
    int ip = -1;
    for (size_t i = 0; i < s.networkNodes().size(); ++i) if (s.networkNodes()[i].name == "PLAT") ip = int(i);
    for (const auto& h : s.networkHistory()) if (ip >= 0 && size_t(ip) < h.p.size()) r.pPlat.push_back(h.p[size_t(ip)]);
    r.ok = true;
    return r;
}
}  // namespace c36n

static void testC36Network() {
    std::printf("== C36 surface network: NETBALAN, automatic compressors (NETCOMPA / COMPOFF)\n");
    if (findExample("NETWORK.DATA").empty()) { std::printf("    (examples/NETWORK.DATA not found - skipped)\n"); return; }
    const std::string comp = "NETCOMPA\n 'PLAT' 'FIELD' 'PLAT' 'OIL' 3 /\n/\n\n";
    const c36n::NetRun B = c36n::run(c36n::networkDeck("", true)), C = c36n::run(c36n::networkDeck(comp, true));
    const size_t on = C.log.find("compressor on branch PLAT -> FIELD switched on");
    std::printf("    GCONPROD PLAT ORAT 6000: FOPT without compressor %.6g, with NETCOMPA %.6g (%+.1f %%)\n", B.fopt, C.fopt, 100 * (C.fopt / B.fopt - 1));
    if (on != std::string::npos) std::printf("    %s\n", C.log.substr(C.log.rfind('\n', on) + 1, C.log.find('\n', on) - C.log.rfind('\n', on) - 1).c_str());
    check(B.ok && C.ok && on != std::string::npos && C.fopt > 1.03 * B.fopt && C.log.find("switched on", on + 48) == std::string::npos,
          "NETCOMPA: the compressor switches on once when the group misses its oil target; the boosted branch table raises the recovery");
    {   // import keeps the keywords; information level
        std::vector<DeckMessage> rep;
        const Case c = eclipseDeckFromString(c36n::networkDeck(comp + "NETBALAN\n 0 0.5 8 /\n\n", true), "", nullptr, &rep);
        int nc = 0, nb = 0, warn = 0;
        for (const auto& e : c.extraKeywords) { nc += e.name == "NETCOMPA"; nb += e.name == "NETBALAN"; }
        for (const auto& m : rep) if ((m.keyword == "NETCOMPA" || m.keyword == "NETBALAN") && m.level >= 2) ++warn;
        const std::string txt = eclipseDeckToString(c);
        check(nc == 1 && nb == 1 && warn == 0 && txt.find("NETCOMPA") != std::string::npos && txt.find(" 'PLAT' 'FIELD' 'PLAT' 'OIL' 3 /") != std::string::npos,
              "NETCOMPA / NETBALAN imported without errors and exported verbatim");
    }
    // multi-level compressor with fuel gas: one level per time step, FGCR = fuel x level / levels
    {
        const c36n::NetRun M = c36n::run(c36n::networkDeck("NETCOMPA\n 'PLAT' 'FIELD' 'PLAT' 'OIL' 3 0 900 'PLAT' 'MULTI' 3 0 /\n/\n\n", true));
        std::vector<double> lv;
        for (double g : M.fgcr) if (g > 0 && (lv.empty() || std::fabs(g - lv.back()) > 1e-6)) lv.push_back(g);
        std::printf("    multi-level: FGCR steps"); for (double g : lv) std::printf(" %.1f", g); std::printf(" Mscf/d, FOPT %.6g\n", M.fopt);
        check(lv.size() >= 3 && std::fabs(lv[0] - 300) < 1e-6 && std::fabs(lv[1] - 600) < 1e-6 && std::fabs(lv[2] - 900) < 1e-6 && M.fgcrMax <= 900 + 1e-6 &&
              M.log.find("up to level 3 of 3") != std::string::npos,
              "NETCOMPA MULTI: levels 1..3 on successive time steps, fuel gas 900 x level / 3 reported in FGCR");
    }
    // COMPOFF: switched off at the given report step (and on again when the target is still missed)
    {
        std::string d = c36n::networkDeck(comp, true);
        const size_t p = d.find("TSTEP\n 40*91.3125");
        d.replace(p, std::string("TSTEP\n 40*91.3125\n/").size(), "TSTEP\n 20*91.3125\n/\n\nCOMPOFF\n\nTSTEP\n 20*91.3125\n/");
        const c36n::NetRun O = c36n::run(d);
        const size_t a = O.log.find("COMPOFF: 1 automatic compressor(s) switched off");
        check(a != std::string::npos && O.log.find("switched on", a) != std::string::npos, "COMPOFF switches the compressors off; NETCOMPA switches them on again while the target is missed");
    }
    // NETBALAN: interval 365 d holds the node pressure between balances; tolerance / iterations
    {
        const c36n::NetRun Y = c36n::run(c36n::networkDeck("NETBALAN\n 365 /\n\n", false));
        const c36n::NetRun T = c36n::run(c36n::networkDeck("NETBALAN\n 0 0.01 20 /\n\n", false));
        const c36n::NetRun D = c36n::run(c36n::networkDeck("", false));
        int changes = 0;
        for (size_t i = 1; i < Y.pPlat.size(); ++i) changes += std::fabs(Y.pPlat[i] - Y.pPlat[i - 1]) > 1e-9;
        int changesD = 0;
        for (size_t i = 1; i < D.pPlat.size(); ++i) changesD += std::fabs(D.pPlat[i] - D.pPlat[i - 1]) > 1e-9;
        std::printf("    NETBALAN 365 d: PLAT pressure changed %d times over %zu steps (every step: %d); FOPT %.6g | tight 0.01 psi: %.6g | default: %.6g\n",
                    changes, Y.pPlat.size(), changesD, Y.fopt, T.fopt, D.fopt);
        check(Y.ok && changes <= 11 && changesD > 2 * changes && std::fabs(T.fopt - D.fopt) < 5e-3 * D.fopt,
              "NETBALAN: interval balancing (node pressures held for 365 d), node tolerance / max iterations used by the implicit balance");
    }
    // restart with the compressor on (state v11: level + last balance time)
    {
        const std::string d = c36n::networkDeck(comp, true);
        const c36n::NetRun A = c36n::run(d, 20);
        const c36n::NetRun R = c36n::run(d, -1, A.blob);
        std::printf("    restart after report step 20 (compressor on): FOPT %.3f vs continuous %.3f\n", R.fopt, C.fopt);
        {   // C36 state v11 (network); later rounds bump the version - any state >= 11 carries the network fields
            const size_t q = A.blob.find("RESSIM-BO-STATE-");
            const int ver = q == std::string::npos ? 0 : std::atoi(A.blob.c_str() + q + 16);
            check(ver >= 11, "restart state version >= 11 (network state)");
        }
        check(std::fabs(R.fopt - C.fopt) < 1e-6 * C.fopt, "restart with an automatic compressor on reproduces the continuous run");
    }
}
