// C36 headless harness of the Geomodel 3-D scene builders and the seismic interpretation UI (ui_geo3d.cpp).
// The D3D11 renderer is replaced by a recording mock (geoscene3d.h in this directory); ImGui / ImPlot run without a
// backend so the full UI code path (geomodel3D, geoSeismicSection) executes for every attribute and pick mode.
// Build + run: bash tests/geo3d_headless/run.sh
#include "ui_geo3d.cpp"
#include "ressim/geomodel.h"
#include "ressim/seismic3d.h"
#include <filesystem>
#include <iostream>
#include <thread>

int g_plotPalette = 0;
namespace fig {
uint32_t distinctColor(int i, int) {
    static const uint32_t c[] = {0xFF1F77B4u, 0xFFFF7F0Eu, 0xFF2CA02Cu, 0xFFD62728u, 0xFF9467BDu, 0xFF8C564Bu, 0xFFE377C2u, 0xFF7F7F7Fu};
    return c[unsigned(i) % 8u];
}
}  // namespace fig
bool fileDialog(bool, const char*, const char*, std::string&) { return false; }
std::vector<double> geoGridProperty(const App::GeoState& S, int which) {
    std::vector<double> v(S.g.bulk.size(), std::nan(""));
    for (size_t c = 0; c < v.size(); ++c) {
        if (!S.g.grid.actnum.empty() && !S.g.grid.actnum[c]) continue;
        v[c] = which == 0 ? S.g.top[c] : S.g.bot[c] - S.g.top[c];
    }
    return v;
}
const char* kGeoGridProps[7] = {"Top depth", "Facies", "Porosity", "log10 permeability X", "Fault block", "Zone", "Cell thickness"};

static int fails = 0, checks = 0;
static void check(const char* name, bool ok, const std::string& detail = "") {
    ++checks; if (!ok) ++fails;
    std::printf("  [%s] %s %s\n", ok ? "PASS" : "FAIL", name, detail.c_str());
}
static std::string fmt(const char* f, double a, double b = 0, double c = 0) { char s[256]; std::snprintf(s, sizeof s, f, a, b, c); return s; }

static void loadDemo(App::GeoState& S) {
    const std::string dir = (std::filesystem::temp_directory_path() / "ResGeo_h3d").string();
    std::string e;
    writeDemoProject(dir, false, e);
    for (const char* h : {"TopRes", "MidRes", "BaseRes"}) { Surface s = readIrap(dir + "/" + h + ".irap", e); s.name = h; S.hz.push_back(s); }
    S.faults = readFaultSticks(dir + "/F1_sticks.txt", e);
    S.wells = readWellHeads(dir + "/well_heads.txt", e);
    for (auto& w : S.wells) {
        if (std::filesystem::exists(dir + "/" + w.name + "_survey.txt")) readDeviation(dir + "/" + w.name + "_survey.txt", w, e);
        w.logs.push_back(readLas(dir + "/" + w.name + ".las", e));
    }
    for (const auto& t : readTops(dir + "/well_tops.txt", e)) for (auto& w : S.wells) if (w.name == t.first) w.tops.push_back(t.second);
    S.cube = readSegy(dir + "/seismic.segy", e); S.hasCube = S.cube.ni > 0;
    S.vm.v0 = 1800; S.vm.k = 0.6; S.vm.datum = 0;
    S.gs = GridSpec(); S.gs.dx = S.gs.dy = 100; S.gs.nx = 40; S.gs.ny = 30;
    Zone z1; z1.name = "Upper"; z1.top = 0; z1.base = 1; z1.nLayers = 8;
    Zone z2; z2.name = "Lower"; z2.top = 1; z2.base = 2; z2.nLayers = 8;
    S.gs.zones = {z1, z2}; S.gs.tieRange = 1500;
    const Surface& h = S.hz[0];
    S.gs.x0 = h.x0; S.gs.y0 = h.y0; S.gs.rotDeg = h.rotDeg;
    S.g = buildGrid(S.gs, S.hz, S.faults, S.wells, e); S.hasGrid = S.g.grid.nx > 0;
}

static bool runFrame(App& app, App::GeoState& S, int view) {
    ImGuiIO& io = ImGui::GetIO();
    io.DisplaySize = ImVec2(1600, 1000); io.DeltaTime = 1.0f / 60;
    ImGui::NewFrame();
    ImGui::SetNextWindowPos(ImVec2(0, 0)); ImGui::SetNextWindowSize(ImVec2(1600, 1000));
    ImGui::Begin("h");
    if (view == 1) app.geoSeismicSection(S, ImVec2(1200, 800)); else app.geomodel3D(S, ImVec2(1200, 800));
    ImGui::End();
    ImGui::Render();
    return ImGui::GetDrawData() != nullptr;
}

int main() {
    std::printf("C36 geomodel 3-D scene - headless checks\n");
    App::GeoState S; loadDemo(S);
    check("demo project + grid", S.hasCube && S.hz.size() == 3 && S.hasGrid && !S.faults.empty(),
          fmt("cube %.0f x %.0f, grid cells %.0f", S.cube.ni, S.cube.nj, double(S.g.bulk.size())));
    Geo3DView V; V.sc.init();
    SeisGeom g; seisGeom(S, g);
    V.ox = std::round(g.x0); V.oy = std::round(g.y0); V.originSet = true;
    seisAttrRange(SeisAttr::Amplitude, S.cube.a, 0.99, V.rLo, V.rHi);
    S.il = g.ni / 2; S.xl = g.nj / 2; V.showTS = true;
    float lo[3], hi[3]; sceneBounds(S, V, lo, hi);
    {
        double x, y; g.xy(g.ni - 1, g.nj - 1, x, y);
        const bool inX = x - V.ox >= lo[0] - 1e-3 && x - V.ox <= hi[0] + 1e-3 && y - V.oy >= lo[1] - 1e-3 && y - V.oy <= hi[1] + 1e-3;
        const double zt = timeToDepth(S.vm, g.tEnd());
        check("scene bounds contain the seismic lattice", inX && hi[2] >= zt - 1e-3 && lo[2] <= timeToDepth(S.vm, g.t0) + 1e-3, fmt("z %.1f .. %.1f", lo[2], hi[2]));
        check("local origin keeps coordinates small", std::fabs(lo[0]) < 1e5 && std::fabs(hi[1]) < 1e5);
    }
    // ---- slices
    buildSlices(S, V, g);
    {
        const auto& sl = V.sc.slices;
        check("inline + crossline + time slice built", sl.size() == 3, fmt("%.0f slices", double(sl.size())));
        if (sl.size() == 3) {
            const auto& A = sl[0];
            double maxd = 0;
            for (int s = 0; s < g.ns; s += 7) for (int j = 0; j < g.nj; j += 3) maxd = std::max(maxd, double(std::fabs(A.tex[size_t(s) * size_t(g.nj) + size_t(j)] - S.cube.at(S.il, j, s))));
            double x0, y0; g.xy(S.il, 0, x0, y0);
            const double dz0 = std::fabs(A.pos[2] - timeToDepth(S.vm, g.t0)), dzN = std::fabs(A.pos[size_t(A.nu) * size_t(A.nv - 1) * 3 + 2] - timeToDepth(S.vm, g.tEnd()));
            check("inline texture = cube samples (transposed)", A.tw == g.nj && A.th == g.ns && maxd == 0, fmt("max diff %.3g", maxd));
            check("inline slice geometry (plan + depth conversion)", std::fabs(A.pos[0] - (x0 - V.ox)) < 1e-2 && std::fabs(A.pos[1] - (y0 - V.oy)) < 1e-2 && dz0 < 1e-2 && dzN < 1e-2,
                  fmt("dz top %.3g bottom %.3g", dz0, dzN));
            const auto& B = sl[1];
            double maxb = 0;
            for (int s = 0; s < g.ns; s += 5) for (int i = 0; i < g.ni; i += 2) maxb = std::max(maxb, double(std::fabs(B.tex[size_t(s) * size_t(g.ni) + size_t(i)] - S.cube.at(i, S.xl, s))));
            check("crossline texture = cube samples", B.tw == g.ni && maxb == 0, fmt("max diff %.3g", maxb));
            const auto& T = sl[2];
            double maxt = 0;
            for (int i = 0; i < g.ni; ++i) for (int j = 0; j < g.nj; ++j) maxt = std::max(maxt, double(std::fabs(T.tex[size_t(j) * size_t(g.ni) + size_t(i)] - S.cube.at(i, j, V.curTS))));
            check("time slice texture = cube samples", T.tw == g.ni && T.th == g.nj && maxt == 0, fmt("sample %.0f, max diff %.3g", V.curTS, maxt));
        }
    }
    // ---- picking: rays against the slices
    {
        const double tm = g.t0 + 0.43 * (g.tEnd() - g.t0);
        double x, y; g.xy(S.il - 5, 7.25, x, y);
        double o[3] = {x - V.ox, y - V.oy, zOfT(S, V, tm)};
        const double L = std::hypot(g.ix, g.iy); double d[3] = {g.ix / L, g.iy / L, 0};
        int kind = -1; double fi = 0, fj = 0, t = 0, z = 0;
        const bool hit = rayHitSlices(S, V, g, o, d, kind, fi, fj, t, z);
        check("ray along the inline direction hits the inline slice", hit && kind == 0 && std::fabs(fi - S.il) < 1e-6 && std::fabs(fj - 7.25) < 1e-6 && std::fabs(t - tm) < 1e-6,
              fmt("fi %.6f fj %.6f t err %.2g", fi, fj, t - tm));
        g.xy(S.il + 3.3, S.xl + 2.2, x, y);
        double o2[3] = {x - V.ox, y - V.oy, lo[2] - 100}, d2[3] = {0, 0, 1};
        const bool hit2 = rayHitSlices(S, V, g, o2, d2, kind, fi, fj, t, z);
        const double tts = g.t0 + V.curTS * g.dt;
        check("vertical ray hits the time slice", hit2 && kind == 2 && std::fabs(fi - (S.il + 3.3)) < 1e-6 && std::fabs(t - tts) < 1e-6, fmt("fi %.6f t %.4f", fi, t));
        // seed + fault stick from the picks
        App::GeoState P; P.cube = S.cube; P.hasCube = true; P.vm = S.vm;
        addSeedAt(P, g, S.il, 7.25, tm);
        double xs, ys, xc = 0, yc = 0; g.xy(S.il, 7.25, xs, ys);
        if (!P.seeds.empty()) P.cube.xy(P.seeds[0][0], P.seeds[0][1], xc, yc);
        check("pick -> tracking seed (cube indices)", P.seeds.size() == 1 && std::hypot(xc - xs, yc - ys) < 1e-6 && P.seeds[0][2] == tm);
        addStickPoint(P, g, S.il, 5, tm + 40); addStickPoint(P, g, S.il, 6, tm - 30); addStickPoint(P, g, S.il, 5.5, tm);
        finishStick(P);
        const bool sorted = P.draft.sticks.size() == 1 && P.draft.sticks[0].size() == 3 && P.draft.sticks[0][0][2] < P.draft.sticks[0][1][2] && P.draft.sticks[0][1][2] < P.draft.sticks[0][2][2];
        check("fault stick digitising (sorted by depth, depth-converted)", sorted && std::fabs(P.draft.sticks[0][0][2] - timeToDepth(S.vm, tm - 30)) < 1e-9);
    }
    // ---- meshes: horizons + grid
    V.showGrid = true; V.gridProp = 0;
    buildMeshes(S, V);
    {
        const auto& M = V.sc.meshes;
        check("3 horizon meshes + grid mesh", M.size() == 4, fmt("%.0f meshes", double(M.size())));
        if (M.size() == 4) {
            const auto& H = M[0];
            size_t n = H.pos.size() / 3; double maxd = 0; bool idxOk = true;
            for (size_t k = 0; k < n; ++k) if (std::isfinite(H.val[k])) maxd = std::max(maxd, double(std::fabs(H.pos[3 * k + 2] - H.val[k])));
            for (unsigned q : H.idx) idxOk = idxOk && q < n;
            check("horizon mesh: depth z, valid triangles", maxd < 1e-2 && idxOk && H.idx.size() > 3 * 100, fmt("%.0f triangles, z err %.3g", double(H.idx.size() / 3), maxd));
            const auto& G = M[3];
            const size_t faces = G.idx.size() / 6;
            const int nx = S.g.grid.nx, ny = S.g.grid.ny, nz = S.g.grid.nz;
            int act = 0; for (int a : S.g.grid.actnum) act += a != 0;
            double zmin = 1e300, zmax = -1e300;   // corner depths of the active cells
            for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
                if (!S.g.grid.actnum[size_t(i + nx * (j + ny * k))]) continue;
                for (int f = 0; f < 2; ++f) for (int bb = 0; bb < 2; ++bb) for (int aa = 0; aa < 2; ++aa) {
                    const double z = S.g.grid.zcorn[size_t(((2 * k + f) * 2 * ny + 2 * j + bb) * 2 * nx + 2 * i + aa)]; zmin = std::min(zmin, z); zmax = std::max(zmax, z); }
            }
            bool inside = true;
            for (size_t k = 0; k < G.pos.size() / 3; ++k) inside = inside && std::isfinite(G.pos[3 * k]) && G.pos[3 * k + 2] >= zmin - 1 && G.pos[3 * k + 2] <= zmax + 1;
            check("grid shell: exposed faces, positions inside the grid depth range", faces >= size_t(2 * nx * ny) && inside && G.uv.size() == G.pos.size() / 3 * 2,
                  fmt("%.0f faces (%.0f active cells)", double(faces), act) + fmt(", corner z %.1f .. %.1f", zmin, zmax));
            // box filter: one layer -> exactly the top + bottom + side faces of that layer
            V.gk[0] = V.gk[1] = 3; buildMeshes(S, V);
            const size_t f1 = V.sc.meshes.back().idx.size() / 6;
            int actK = 0; for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) actK += S.g.grid.actnum[size_t(i + nx * (j + ny * 3))] != 0;
            check("K filter: one layer shows its top and bottom faces", f1 >= size_t(2 * actK), fmt("%.0f faces, %.0f active cells in layer 3", double(f1), actK));
            V.gk[0] = 0; V.gk[1] = 1 << 20;
        }
    }
    // ---- faults, wells, labels
    buildSolidsLines(S, V, lo, hi);
    {
        size_t wellsWithPath = 0, tops = 0; for (const auto& w : S.wells) { wellsWithPath += w.totalMd() > 0 || !w.logs.empty(); tops += w.tops.size(); }
        int wl = 0, tl = 0; for (const auto& L : V.labels) { if (!L.small) ++wl; else ++tl; }
        check("fault surface + well tubes", V.sc.solids.size() == 2, fmt("%.0f solids", double(V.sc.solids.size())));
        check("well name labels + top labels", wl == int(wellsWithPath) && tl >= int(tops), fmt("%.0f names, %.0f small labels", wl, tl));
        if (V.sc.solids.size() == 2) {
            const auto& W = V.sc.solids[1];
            bool ok = W.pos.size() == W.nrm.size() && W.col.size() / 4 == W.pos.size() / 3;
            for (unsigned q : W.idx) ok = ok && q < W.pos.size() / 3;
            check("well tube buffers consistent", ok, fmt("%.0f vertices", double(W.pos.size() / 3)));
        }
        V.wellLog = 0; buildSolidsLines(S, V, lo, hi);
        check("well tubes coloured by a log", V.sc.solids.size() == 2);
    }
    // ---- volume
    V.showVol = true; buildVolume(S, V);
    check("volume (proxy) desc", V.sc.hasVol && V.sc.vol.ni == S.cube.ni && V.sc.vol.nj == S.cube.nj && V.sc.vol.ns == S.cube.ns && std::fabs(V.sc.vol.x0 - (S.cube.x0 - V.ox)) < 1e-9 && V.sc.vol.zMode == 1,
          fmt("%.0f x %.0f x %.0f", V.sc.vol.ni, V.sc.vol.nj, V.sc.vol.ns));
    {
        double maxd = 0;
        for (int i = 0; i < V.sc.vol.ni; i += 3) for (int j = 0; j < V.sc.vol.nj; j += 3) for (int s = 0; s < V.sc.vol.ns; s += 11)
            maxd = std::max(maxd, double(std::fabs(V.sc.vol.data[(size_t(i) * size_t(V.sc.vol.nj) + size_t(j)) * size_t(V.sc.vol.ns) + size_t(s)] - S.cube.at(i, j, s))));
        check("volume samples = cube", maxd == 0);
    }
    // ---- time domain
    V.domainTime = true; buildSlices(S, V, g);
    check("time domain: slice z = TWT x v0/2000", !V.sc.slices.empty() && std::fabs(V.sc.slices[0].pos[2] - g.t0 * S.vm.v0 / 2000) < 1e-3 && V.sc.slices[0].nv == 2);
    V.domainTime = false;
    // ---- out-of-core source == in-memory source
    {
        const std::string p = (std::filesystem::temp_directory_path() / "h3d_ooc.sgy").string();
        std::string e;
        writeSyntheticSegy(p, 48, 40, 150, 2.0, 1, false, e);
        App::GeoState O;
        O.segy = std::make_shared<SegyVolume>();
        const bool op = O.segy->open(p, e);
        O.cube = O.segy->toCube(); O.hasCube = true;
        App::GeoState M; M.cube = O.cube; M.hasCube = true;
        double md = 0;
        for (int k = 0; k < 3; ++k) {
            const auto a = seisSlice(O, k, 9), b = seisSlice(M, k, 9);
            if (a.size() != b.size()) md = 1e30;
            else for (size_t q = 0; q < a.size(); ++q) md = std::max(md, double(std::fabs(a[q] - b[q])));
        }
        const auto ca = seisAttrSection(O, 0, 20, SeisAttr::Coherence), cb = seisAttrSection(M, 0, 20, SeisAttr::Coherence);
        double mc = ca.size() == cb.size() ? 0 : 1e30; for (size_t q = 0; q < ca.size() && q < cb.size(); ++q) mc = std::max(mc, double(std::fabs(ca[q] - cb[q])));
        check("out-of-core slices = in-memory slices (inline, crossline, time)", op && md == 0, fmt("max diff %.3g", md));
        check("out-of-core coherence section = in-memory", mc < 1e-6, fmt("max diff %.3g", mc));
        // proxy: seeds picked on the full-resolution source land on the same plan point of the proxy
        O.cube = O.segy->proxy(16, 64);
        SeisGeom og; seisGeom(O, og);
        addSeedAt(O, og, 30, 21, 100);
        double xs, ys, xc = 0, yc = 0; og.xy(30, 21, xs, ys); if (!O.seeds.empty()) O.cube.xy(O.seeds[0][0], O.seeds[0][1], xc, yc);
        check("seed on a full-res section -> proxy cube indices", O.seeds.size() == 1 && std::hypot(xc - xs, yc - ys) < 1e-6, fmt("proxy %.0f x %.0f", O.cube.ni, O.cube.nj));
        std::filesystem::remove(p);
    }
    // ---- C36-4 horizon brush: lattice mapping, ray / horizon intersection, brush + flatten
    {
        const std::vector<Surface> hz0 = S.hz;
        const Surface& h = S.hz[0];
        double xn, yn, xl, yl, fi, fj; h.nodeXY(7, 5, xn, yn); geoLatticeXY(h, 7, 5, xl, yl);
        double xf, yf; geoLatticeXY(h, 7.3, 5.6, xf, yf); surfaceIJ(h, xf, yf, fi, fj);
        check("lattice -> world (rotation aware) inverts surfaceIJ", std::hypot(xn - xl, yn - yl) < 1e-6 && std::fabs(fi - 7.3) < 1e-9 && std::fabs(fj - 5.6) < 1e-9, fmt("node err %.2e, fi %.6f fj %.6f", std::hypot(xn - xl, yn - yl), fi, fj));
        Geo3DView V2; V2.ox = std::round(h.x0); V2.oy = std::round(h.y0);
        double xw, yw; geoLatticeXY(h, h.nx * 0.45, h.ny * 0.55, xw, yw);
        const double zt = h.at(xw, yw);
        double o[3] = {xw - V2.ox, yw - V2.oy, zt - 700}, d[3] = {0, 0, 1}, X, Y, D;
        const bool hv = rayHitHorizon(S, V2, 0, o, d, X, Y, D);
        double dd[3] = {0.35, -0.2, 1}; const double nn = std::sqrt(dd[0] * dd[0] + dd[1] * dd[1] + 1); for (double& c : dd) c /= nn;
        double o2[3] = {xw - V2.ox - 300 * dd[0], yw - V2.oy - 300 * dd[1], zt - 300 * dd[2]}, X2, Y2, D2;
        const bool ho = rayHitHorizon(S, V2, 0, o2, dd, X2, Y2, D2);
        const double tt = (X2 - V2.ox - o2[0]) / dd[0], zr = o2[2] + tt * dd[2], yr = o2[1] + tt * dd[1] + V2.oy;
        check("vertical / oblique ray hits the horizon surface", hv && std::fabs(D - zt) < 1e-6 && std::hypot(X - xw, Y - yw) < 1e-6 && ho && std::fabs(zr - D2) < 1e-3 && std::fabs(yr - Y2) < 1e-6,
              fmt("vertical dz %.2e, oblique dz %.2e", std::fabs(D - zt), std::fabs(zr - D2)));
        double om[3] = {o[0] + 1e6, o[1], o[2]};
        check("ray outside the horizon misses", !rayHitHorizon(S, V2, 0, om, d, X, Y, D));
        S.editHz = 0; S.brushMode = 0; S.brushAmt = 10; S.brushR = 300; S.hzUndo.clear();
        const int rev0 = S.hzRev;
        const int n = geoApplyBrush(S, xw, yw);
        const double dz = S.hz[0].at(xw, yw) - zt;
        S.brushMode = 2; S.flatDepth = 0;
        double xw2, yw2; geoLatticeXY(h, h.nx * 0.3, h.ny * 0.3, xw2, yw2);
        const double z2 = S.hz[0].at(xw2, yw2);
        geoApplyBrush(S, xw2, yw2);
        double xw3, yw3; geoLatticeXY(h, h.nx * 0.3 + 1, h.ny * 0.3, xw3, yw3);
        geoApplyBrush(S, xw3, yw3);
        check("brush raise / lower + flatten (depth of the first click), undo entries, cache revision", n > 10 && dz > 5 && dz <= 10 + 1e-9 && std::fabs(S.flatDepth - z2) < 1e-9 &&
              S.hzUndo.size() == 3 && S.hzRev == rev0 + 3 && S.hzUndo[0][0].z == hz0[0].z,
              fmt("%.0f nodes, dz %.3f, flatten depth err %.2e", n, dz, std::fabs(S.flatDepth - z2)));
        S.brushMode = 0;
        S.hz = hz0; S.hzUndo.clear();
    }
    // ---- headless UI frames (ImGui + ImPlot, no renderer)
    {
        ImGui::CreateContext(); ImPlot::CreateContext();
        ImGuiIO& io = ImGui::GetIO();
        io.IniFilename = nullptr;
        unsigned char* px; int tw, th; io.Fonts->GetTexDataAsRGBA32(&px, &tw, &th);
        App app;
        S.v3.reset();
        bool ok = true;
        for (int a = 0; a < int(SeisAttr::Count); ++a) {
            S.secAttr = a;
            for (int pm = 0; pm < 4; ++pm) { S.pickMode = pm; S.secOrient = a % 2; S.flatten = a % 3 == 0; ok = ok && runFrame(app, S, 1); }
            ok = ok && runFrame(app, S, 4);
        }
        // the attribute cube job: wait and run frames until it lands
        S.secAttr = int(SeisAttr::Envelope); S.v3->showTS = true; S.v3->showVol = true;
        for (int f = 0; f < 400 && ok; ++f) { ok = runFrame(app, S, 4); if (S.v3->attrCube && S.v3->attrKey == S.v3->attrJobKey) break; std::this_thread::sleep_for(std::chrono::milliseconds(20)); }
        ok = ok && runFrame(app, S, 4);
        check("UI frames: section (10 attributes x 4 pick modes) + 3-D scene", ok);
        check("attribute cube computed in the background and shown", S.v3->attrCube && S.v3->sc.hasVol && S.v3->tsFromAttr, fmt("rebuilds: slices %.0f, volume %.0f", S.v3->sc.nSlices, S.v3->sc.nVol));
        const int before = S.v3->sc.nSlices; runFrame(app, S, 4); runFrame(app, S, 4);
        check("no rebuild when nothing changed", S.v3->sc.nSlices == before);
        S.il = std::min(S.il + 1, S.cube.ni - 1); runFrame(app, S, 4);
        check("moving the inline rebuilds the slices", S.v3->sc.nSlices == before + 1);
        S.pickMode = 3; runFrame(app, S, 4);
        const int m0 = S.v3->sc.nMeshes; runFrame(app, S, 4);
        const bool same = S.v3->sc.nMeshes == m0;
        ++S.hzRev; runFrame(app, S, 4);
        check("edit-horizon mode frames; a horizon edit (revision) rebuilds the meshes only then", same && S.v3->sc.nMeshes == m0 + 1, fmt("meshes %.0f -> %.0f", m0, S.v3->sc.nMeshes));
        S.pickMode = 0;
        ImPlot::DestroyContext(); ImGui::DestroyContext();
    }
    std::printf("C36 geo3d headless: %d / %d PASS\n", checks - fails, checks);
    return fails ? 1 : 0;
}
