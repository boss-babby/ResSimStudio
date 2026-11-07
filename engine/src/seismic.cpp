// C32 geomodelling: SEG-Y rev 1 read / write, seeded horizon auto-tracking, velocity model, attribute extraction
#include "ressim/geomodel.h"
#include "ressim/geo_import.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <deque>
#include <fstream>
#include <limits>
#include <map>

namespace ressim {
namespace geo {

namespace {
const double kNaN = std::numeric_limits<double>::quiet_NaN();
uint16_t be16(const unsigned char* p) { return uint16_t((p[0] << 8) | p[1]); }
int32_t be32(const unsigned char* p) { return int32_t((uint32_t(p[0]) << 24) | (uint32_t(p[1]) << 16) | (uint32_t(p[2]) << 8) | uint32_t(p[3])); }
void put16(unsigned char* p, int v) { p[0] = (unsigned char)((v >> 8) & 255); p[1] = (unsigned char)(v & 255); }
void put32(unsigned char* p, int32_t v) { const uint32_t u = uint32_t(v); p[0] = (unsigned char)(u >> 24); p[1] = (unsigned char)(u >> 16); p[2] = (unsigned char)(u >> 8); p[3] = (unsigned char)u; }
float ibm2ieee(const unsigned char* p) {
    const uint32_t u = (uint32_t(p[0]) << 24) | (uint32_t(p[1]) << 16) | (uint32_t(p[2]) << 8) | uint32_t(p[3]);
    if ((u & 0x7fffffff) == 0) return 0.0f;
    const int sign = (u >> 31) ? -1 : 1, expo = int((u >> 24) & 0x7f) - 64;
    const double frac = double(u & 0x00ffffff) / 16777216.0;
    return float(sign * frac * std::pow(16.0, expo));
}
float ieeeBE(const unsigned char* p) {
    const uint32_t u = (uint32_t(p[0]) << 24) | (uint32_t(p[1]) << 16) | (uint32_t(p[2]) << 8) | uint32_t(p[3]);
    float f; std::memcpy(&f, &u, 4); return f;
}
// EBCDIC -> ASCII for the printable range of the textual header (others become spaces)
char ebcdic(unsigned char c) {
    static const char* tbl =
        "                                                                "      // 0x00-0x3f
        "           .<(+|&         !$*); -/        ,%_>?         `:#@'=\""   // 0x40-0x7f
        " abcdefghi       jklmnopqr       ~stuvwxyz                      "    // 0x80-0xbf
        "{ABCDEFGHI      }JKLMNOPQR      \\ STUVWXYZ      0123456789      ";   // 0xc0-0xff
    return tbl[c];
}
}  // namespace

bool SeismicCube::ij(double x, double y, double& i, double& j) const {
    const double det = ix * jy - iy * jx;
    if (std::fabs(det) < 1e-12) return false;
    const double u = x - x0, v = y - y0;
    i = (u * jy - v * jx) / det; j = (-u * iy + v * ix) / det;
    return i >= -1e-9 && j >= -1e-9 && i <= ni - 1 + 1e-9 && j <= nj - 1 + 1e-9;
}
float SeismicCube::sampleT(double i, double j, double t) const {
    if (ni < 1 || nj < 1 || ns < 2) return float(kNaN);
    const double s = (t - t0) / dt;
    if (i < 0 || j < 0 || s < 0 || i > ni - 1 || j > nj - 1 || s > ns - 1) return float(kNaN);
    const int i0 = std::min(int(i), std::max(ni - 2, 0)), j0 = std::min(int(j), std::max(nj - 2, 0)), s0 = std::min(int(s), ns - 2);
    const double fi = ni > 1 ? i - i0 : 0, fj = nj > 1 ? j - j0 : 0, fs = s - s0;
    auto A = [&](int a, int b, int c) { return double(at(std::min(a, ni - 1), std::min(b, nj - 1), c)); };
    double v = 0;
    for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int c = 0; c < 2; ++c)
        v += (a ? fi : 1 - fi) * (b ? fj : 1 - fj) * (c ? fs : 1 - fs) * A(i0 + a, j0 + b, s0 + c);
    return float(v);
}

SeismicCube readSegy(const std::string& path, std::string& err) { return readSegy(path, SegyByteMap(), err); }
SeismicCube readSegy(const std::string& path, const SegyByteMap& bm, std::string& err) {
    SeismicCube C;
    auto B = [](int byte1) { return size_t(std::max(1, std::min(237, byte1)) - 1); };
    std::ifstream f(path, std::ios::binary);
    if (!f) { err = "cannot open " + path; return C; }
    f.seekg(0, std::ios::end);
    const long long size = f.tellg();
    f.seekg(0);
    if (size < 3600 + 240) { err = path + ": too small for SEG-Y"; return C; }
    std::vector<unsigned char> th(3200), bh(400);
    f.read(reinterpret_cast<char*>(th.data()), 3200);
    f.read(reinterpret_cast<char*>(bh.data()), 400);
    bool ascii = true;
    for (int k = 0; k < 80; ++k) if (th[size_t(k)] >= 0x80) ascii = false;
    C.textHeader.resize(3200);
    for (int k = 0; k < 3200; ++k) C.textHeader[size_t(k)] = ascii ? char(th[size_t(k)]) : ebcdic(th[size_t(k)]);
    const int dtus = be16(&bh[16]), ns = be16(&bh[20]), fmt = be16(&bh[24]);
    int ext = be16(&bh[304]);   // rev 1: number of extended textual headers
    if (ext < 0 || ext > 100) ext = 0;
    if (fmt != 1 && fmt != 5) { err = path + ": sample format " + std::to_string(fmt) + " not supported (IBM float 1 / IEEE float 5)"; return C; }
    if (ns <= 0) { err = path + ": zero samples per trace"; return C; }
    C.ns = ns; C.dt = dtus / 1000.0; C.format = fmt;
    const long long start = 3600 + 3200LL * ext, tlen = 240 + 4LL * ns;
    const long long ntr = (size - start) / tlen;
    struct Tr { int il, xl; double x, y; long long off; };
    std::vector<Tr> tr; tr.reserve(size_t(ntr));
    std::vector<unsigned char> h(240);
    for (long long t = 0; t < ntr; ++t) {
        f.seekg(start + t * tlen);
        f.read(reinterpret_cast<char*>(h.data()), 240);
        int il = be32(&h[B(bm.inlineByte)]), xl = be32(&h[B(bm.xlineByte)]);
        if (il == 0 && xl == 0) { il = be32(&h[8]); xl = be32(&h[20]); }   // fallback: field record / CDP number
        int sc = int(int16_t(be16(&h[B(bm.scalarByte)])));
        double scale = sc == 0 ? 1.0 : sc > 0 ? double(sc) : 1.0 / double(-sc);
        double x = be32(&h[B(bm.xByte)]) * scale, y = be32(&h[B(bm.yByte)]) * scale;
        if (bm.fallbackXY && x == 0 && y == 0) { x = be32(&h[72]) * scale; y = be32(&h[76]) * scale; }
        if (t == 0) { const int d0 = int(int16_t(be16(&h[108]))); C.t0 = d0; }
        tr.push_back({il, xl, x, y, start + t * tlen + 240});
    }
    if (tr.empty()) { err = path + ": no traces"; return C; }
    std::vector<int> ils, xls;
    for (auto& t : tr) { ils.push_back(t.il); xls.push_back(t.xl); }
    std::sort(ils.begin(), ils.end()); ils.erase(std::unique(ils.begin(), ils.end()), ils.end());
    std::sort(xls.begin(), xls.end()); xls.erase(std::unique(xls.begin(), xls.end()), xls.end());
    C.il0 = ils.front(); C.xl0 = xls.front();
    C.dil = ils.size() > 1 ? ils[1] - ils[0] : 1; C.dxl = xls.size() > 1 ? xls[1] - xls[0] : 1;
    C.ni = (ils.back() - ils.front()) / std::max(1, C.dil) + 1; C.nj = (xls.back() - xls.front()) / std::max(1, C.dxl) + 1;
    if (double(C.ni) * C.nj * C.ns > 4e8) { err = path + ": volume too large for this machine"; return SeismicCube(); }
    C.a.assign(size_t(C.ni) * size_t(C.nj) * size_t(C.ns), 0.0f);
    std::vector<unsigned char> buf(size_t(4 * ns));
    std::map<std::pair<int, int>, std::pair<double, double>> pos;
    for (auto& t : tr) {
        const int i = (t.il - C.il0) / std::max(1, C.dil), j = (t.xl - C.xl0) / std::max(1, C.dxl);
        if (i < 0 || j < 0 || i >= C.ni || j >= C.nj) continue;
        f.seekg(t.off);
        f.read(reinterpret_cast<char*>(buf.data()), 4 * ns);
        for (int s = 0; s < ns; ++s) C.a[(size_t(i) * size_t(C.nj) + size_t(j)) * size_t(ns) + size_t(s)] = fmt == 1 ? ibm2ieee(&buf[size_t(4 * s)]) : ieeeBE(&buf[size_t(4 * s)]);
        pos[{i, j}] = {t.x, t.y};
    }
    // plan geometry from three corner traces (least squares over all would be more robust; corners are exact for regular surveys)
    auto P = [&](int i, int j, double& x, double& y) { auto it = pos.find({i, j}); if (it == pos.end()) return false; x = it->second.first; y = it->second.second; return true; };
    double xa, ya, xb, yb, xc, yc;
    if (P(0, 0, xa, ya)) {
        C.x0 = xa; C.y0 = ya;
        if (C.ni > 1 && P(C.ni - 1, 0, xb, yb)) { C.ix = (xb - xa) / (C.ni - 1); C.iy = (yb - ya) / (C.ni - 1); }
        if (C.nj > 1 && P(0, C.nj - 1, xc, yc)) { C.jx = (xc - xa) / (C.nj - 1); C.jy = (yc - ya) / (C.nj - 1); }
    }
    const size_t sl = path.find_last_of("/\\");
    C.name = path.substr(sl == std::string::npos ? 0 : sl + 1);
    return C;
}

bool writeSegy(const std::string& path, const SeismicCube& C) {
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    std::vector<unsigned char> th(3200, ' ');
    std::string txt = C.textHeader.empty() ? std::string("C 1 ResSim / ResGeo SEG-Y rev1, IEEE float") : C.textHeader;
    for (size_t k = 0; k < 3200 && k < txt.size(); ++k) th[k] = (unsigned char)txt[k];
    std::fwrite(th.data(), 1, 3200, f);
    std::vector<unsigned char> bh(400, 0);
    put16(&bh[16], int(std::lround(C.dt * 1000))); put16(&bh[20], C.ns); put16(&bh[24], 5); put16(&bh[28], 1);
    put16(&bh[300], 0x0100); put16(&bh[302], 1);
    std::fwrite(bh.data(), 1, 400, f);
    std::vector<unsigned char> h(240), d(size_t(4 * C.ns));
    int seq = 0;
    for (int i = 0; i < C.ni; ++i)
        for (int j = 0; j < C.nj; ++j) {
            std::fill(h.begin(), h.end(), 0);
            ++seq;
            put32(&h[0], seq); put32(&h[4], seq);
            put32(&h[8], C.il0 + i * C.dil); put32(&h[20], C.xl0 + j * C.dxl);
            put16(&h[70], -100);
            double x, y; C.xy(i, j, x, y);
            put32(&h[72], int32_t(std::lround(x * 100))); put32(&h[76], int32_t(std::lround(y * 100)));
            put32(&h[180], int32_t(std::lround(x * 100))); put32(&h[184], int32_t(std::lround(y * 100)));
            put32(&h[188], C.il0 + i * C.dil); put32(&h[192], C.xl0 + j * C.dxl);
            put16(&h[108], int(std::lround(C.t0)));
            put16(&h[114], C.ns); put16(&h[116], int(std::lround(C.dt * 1000)));
            std::fwrite(h.data(), 1, 240, f);
            for (int s = 0; s < C.ns; ++s) {
                const float v = C.at(i, j, s); uint32_t u; std::memcpy(&u, &v, 4);
                d[size_t(4 * s)] = (unsigned char)(u >> 24); d[size_t(4 * s + 1)] = (unsigned char)(u >> 16); d[size_t(4 * s + 2)] = (unsigned char)(u >> 8); d[size_t(4 * s + 3)] = (unsigned char)u;
            }
            std::fwrite(d.data(), 1, d.size(), f);
        }
    std::fclose(f);
    return true;
}

// seeded flood-fill tracking: from each picked trace, neighbours are picked at the extremum within +-window of the
// neighbour's predicted time (parabolic sub-sample refinement); a pick is rejected if it jumps more than maxJump
TimeHorizon autoTrack(const SeismicCube& c, const std::vector<std::array<double, 3>>& seeds, int mode, double window, double maxJump) {
    TimeHorizon H; H.ni = c.ni; H.nj = c.nj; H.t.assign(size_t(c.ni) * size_t(c.nj), kNaN);
    if (c.ns < 3) return H;
    auto pick = [&](int i, int j, double tPred) -> double {
        const int s0 = std::max(1, int(std::floor((tPred - window - c.t0) / c.dt))), s1 = std::min(c.ns - 2, int(std::ceil((tPred + window - c.t0) / c.dt)));
        double best = -1e300; int bs = -1;
        for (int s = s0; s <= s1; ++s) {
            const double a = c.at(i, j, s - 1), b = c.at(i, j, s), d = c.at(i, j, s + 1);
            if (mode == 2) { if (b >= 0 && d < 0) { const double t = c.t0 + c.dt * (s + b / (b - d)); const double sc = -std::fabs(t - tPred); if (sc > best) { best = sc; bs = s; } } continue; }
            const double v = mode == 0 ? b : -b;
            if ((mode == 0 ? (b >= a && b >= d) : (b <= a && b <= d)) && v > best) { best = v; bs = s; }
        }
        if (bs < 0) return kNaN;
        if (mode == 2) {
            const double b = c.at(i, j, bs), d = c.at(i, j, bs + 1);
            return c.t0 + c.dt * (bs + b / (b - d));
        }
        const double a = c.at(i, j, bs - 1), b = c.at(i, j, bs), d = c.at(i, j, bs + 1), den = a - 2 * b + d;
        const double off = std::fabs(den) > 1e-20 ? 0.5 * (a - d) / den : 0.0;
        return c.t0 + c.dt * (bs + std::max(-0.5, std::min(0.5, off)));
    };
    // amplitude of a pick (trilinear at the picked time); a neighbour pick weaker than 35 % of its parent is rejected
    // (no event there - noise or a side lobe), which stops the tracker at faults and terminations instead of jumping
    std::vector<double> amp(H.t.size(), 0.0);
    auto ampAt = [&](int i, int j, double t) { const double v = c.sampleT(i, j, t); return std::isfinite(v) ? (mode == 1 ? -v : v) : 0.0; };
    std::deque<std::pair<int, int>> q;
    for (const auto& s : seeds) {
        const int i = int(std::lround(s[0])), j = int(std::lround(s[1]));
        if (i < 0 || j < 0 || i >= c.ni || j >= c.nj) continue;
        const double t = pick(i, j, s[2]);
        if (std::isfinite(t)) { H.t[size_t(i * c.nj + j)] = t; amp[size_t(i * c.nj + j)] = ampAt(i, j, t); q.push_back({i, j}); }
    }
    const int di[4] = {1, -1, 0, 0}, dj[4] = {0, 0, 1, -1};
    while (!q.empty()) {
        auto [i, j] = q.front(); q.pop_front();
        const double t0 = H.t[size_t(i * c.nj + j)];
        for (int d = 0; d < 4; ++d) {
            const int a = i + di[d], b = j + dj[d];
            if (a < 0 || b < 0 || a >= c.ni || b >= c.nj || std::isfinite(H.t[size_t(a * c.nj + b)])) continue;
            const double t = pick(a, b, t0);
            if (!std::isfinite(t) || std::fabs(t - t0) > maxJump) continue;
            const double am = ampAt(a, b, t);
            if (mode != 2 && am < 0.35 * amp[size_t(i * c.nj + j)]) continue;
            H.t[size_t(a * c.nj + b)] = t; amp[size_t(a * c.nj + b)] = am; q.push_back({a, b});
        }
    }
    return H;
}

// v(z) = v0 + k (z - datum): one-way time t(z) = ln(1 + k (z - datum) / v0) / k; TWT ms = 2000 t
double depthToTime(const VelocityModel& v, double z) {
    const double dz = z - v.datum;
    const double t = std::fabs(v.k) > 1e-12 ? std::log(std::max(1e-12, 1 + v.k * dz / v.v0)) / v.k : dz / v.v0;
    return 2000.0 * t;
}
double timeToDepth(const VelocityModel& v, double twt) {
    const double t = twt / 2000.0;
    return v.datum + (std::fabs(v.k) > 1e-12 ? v.v0 * (std::exp(v.k * t) - 1) / v.k : v.v0 * t);
}
Surface horizonToDepth(const SeismicCube& c, const TimeHorizon& h, const VelocityModel& v, const std::string& name) {
    Surface s; s.name = name; s.nx = c.ni; s.ny = c.nj;
    s.x0 = c.x0; s.y0 = c.y0;
    s.dx = std::hypot(c.ix, c.iy); s.dy = std::hypot(c.jx, c.jy);
    s.rotDeg = std::atan2(c.iy, c.ix) * 180.0 / 3.14159265358979323846;
    s.z.assign(size_t(s.nx) * size_t(s.ny), kNaN);
    for (int i = 0; i < c.ni; ++i) for (int j = 0; j < c.nj; ++j) {
        const double t = h.t[size_t(i * c.nj + j)];
        if (std::isfinite(t)) s.z[size_t(i + s.nx * j)] = timeToDepth(v, t);
    }
    // crossline axis left-handed (j to the right of i): mirror into the Surface frame
    const double cr = c.ix * c.jy - c.iy * c.jx;
    if (cr < 0) {   // j points clockwise from i: flip j
        Surface f = s; f.x0 = c.x0 + (c.nj - 1) * c.jx; f.y0 = c.y0 + (c.nj - 1) * c.jy;
        for (int i = 0; i < s.nx; ++i) for (int j = 0; j < s.ny; ++j) f.z[size_t(i + s.nx * j)] = s.z[size_t(i + s.nx * (s.ny - 1 - j))];
        return f;
    }
    return s;
}
std::vector<double> extractAttribute(const SeismicCube& c, const VelocityModel& v, const std::vector<P3>& pts, int kind, double win) {
    std::vector<double> out(pts.size(), kNaN);
    for (size_t k = 0; k < pts.size(); ++k) {
        double i, j;
        if (!c.ij(pts[k][0], pts[k][1], i, j)) continue;
        const double t = depthToTime(v, pts[k][2]);
        if (kind == 0) { out[k] = c.sampleT(i, j, t); continue; }
        double s = 0; int n = 0;
        for (double tt = t - win; tt <= t + win + 1e-9; tt += c.dt) {
            const double a = c.sampleT(i, j, tt);
            if (!std::isfinite(a)) continue;
            s += kind == 1 ? a * a : std::fabs(a); ++n;
        }
        if (n) out[k] = kind == 1 ? std::sqrt(s / n) : s / n;
    }
    return out;
}

}  // namespace geo
}  // namespace ressim
