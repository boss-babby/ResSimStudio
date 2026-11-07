// C32 geomodelling: surfaces (IRAP classic, ZMAP+, XYZ), LAS 2.0, well heads / deviation / tops, fault sticks
#include "ressim/geomodel.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>

namespace ressim {
namespace geo {

namespace {
const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kDeg = 3.14159265358979323846 / 180.0;
std::string trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
    return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
}
std::string upper(std::string s) { for (char& c : s) c = char(std::toupper(static_cast<unsigned char>(c))); return s; }
std::vector<std::string> split(const std::string& s) {   // whitespace / comma / semicolon / tab, quotes kept together
    std::vector<std::string> t; std::string cur; bool q = false;
    for (char c : s) {
        if (c == '"' || c == '\'') { q = !q; continue; }
        if (!q && (c == ' ' || c == '\t' || c == ',' || c == ';' || c == '\r' || c == '\n')) { if (!cur.empty()) { t.push_back(cur); cur.clear(); } continue; }
        cur += c;
    }
    if (!cur.empty()) t.push_back(cur);
    return t;
}
bool isNum(const std::string& s, double& v) {
    if (s.empty()) return false;
    char* e = nullptr; v = std::strtod(s.c_str(), &e);
    return e && *e == 0;
}
bool readAll(const std::string& path, std::string& txt) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    std::ostringstream o; o << f.rdbuf(); txt = o.str();
    return true;
}
}  // namespace

// ================================================================ surfaces
void Surface::nodeXY(int i, int j, double& x, double& y) const {
    const double c = std::cos(rotDeg * kDeg), s = std::sin(rotDeg * kDeg);
    const double u = i * dx, v = j * dy;
    x = x0 + c * u - s * v; y = y0 + s * u + c * v;
}
double Surface::at(double x, double y) const {
    if (!valid()) return kNaN;
    const double c = std::cos(rotDeg * kDeg), s = std::sin(rotDeg * kDeg);
    const double u = (c * (x - x0) + s * (y - y0)) / dx, v = (-s * (x - x0) + c * (y - y0)) / dy;
    if (u < -1e-9 || v < -1e-9 || u > nx - 1 + 1e-9 || v > ny - 1 + 1e-9) return kNaN;
    const int i = std::min(std::max(int(std::floor(u)), 0), nx - 2), j = std::min(std::max(int(std::floor(v)), 0), ny - 2);
    const double fu = u - i, fv = v - j;
    const double z00 = node(i, j), z10 = node(i + 1, j), z01 = node(i, j + 1), z11 = node(i + 1, j + 1);
    // undefined nodes: weights renormalised over the defined ones (edges of a faulted / clipped surface)
    double w[4] = {(1 - fu) * (1 - fv), fu * (1 - fv), (1 - fu) * fv, fu * fv}, zz[4] = {z00, z10, z01, z11}, sw = 0, sz = 0;
    for (int k = 0; k < 4; ++k) if (std::isfinite(zz[k])) { sw += w[k]; sz += w[k] * zz[k]; }
    return sw > 0.25 ? sz / sw : kNaN;
}
double Surface::atBlock(double x, double y, double tx, double ty) const {
    for (int k = 0; k <= 8; ++k) {
        const double f = k / 8.0;
        const double z = at(x + f * (tx - x), y + f * (ty - y));
        if (std::isfinite(z)) return z;
    }
    return kNaN;
}
void Surface::range(double& lo, double& hi) const {
    lo = 1e300; hi = -1e300;
    for (double v : z) if (std::isfinite(v)) { lo = std::min(lo, v); hi = std::max(hi, v); }
}

// IRAP classic ASCII:  -996 ny dx dy / xmin xmax ymin ymax / nx rot xrot yrot / 0 0 0 0 0 0 0 / values (9999900 undefined)
Surface readIrap(const std::string& path, std::string& err) {
    Surface s;
    std::string txt;
    if (!readAll(path, txt)) { err = "cannot open " + path; return s; }
    std::istringstream is(txt);
    std::vector<double> v; std::string tok;
    while (is >> tok) { double d; if (isNum(tok, d)) v.push_back(d); }
    if (v.size() < 19 || std::fabs(v[0] + 996) > 0.5) { err = path + ": not an IRAP classic ASCII grid (-996 header)"; return s; }
    s.ny = int(v[1]); s.dx = v[2]; s.dy = v[3];
    const double xmin = v[4], ymin = v[6];
    s.nx = int(v[8]); s.rotDeg = v[9];
    const double xr = v[10], yr = v[11];
    s.x0 = std::isfinite(xr) && (xr != 0 || yr != 0) && s.rotDeg != 0 ? xr : xmin;
    s.y0 = std::isfinite(yr) && (xr != 0 || yr != 0) && s.rotDeg != 0 ? yr : ymin;
    const size_t n = size_t(s.nx) * size_t(s.ny);
    if (s.nx < 2 || s.ny < 2 || v.size() < 19 + n) { err = path + ": IRAP grid has " + std::to_string(v.size() - 19) + " values, expected " + std::to_string(n); s = Surface(); return s; }
    s.z.resize(n);
    for (size_t k = 0; k < n; ++k) { const double d = v[19 + k]; s.z[k] = d >= 9999900 - 0.5 ? kNaN : d; }
    const size_t sl = path.find_last_of("/\\");
    s.name = path.substr(sl == std::string::npos ? 0 : sl + 1);
    const size_t dt = s.name.find_last_of('.'); if (dt != std::string::npos) s.name = s.name.substr(0, dt);
    return s;
}
bool writeIrap(const std::string& path, const Surface& s) {
    FILE* f = std::fopen(path.c_str(), "w");
    if (!f) return false;
    double xmax = s.x0 + (s.nx - 1) * s.dx, ymax = s.y0 + (s.ny - 1) * s.dy;
    std::fprintf(f, "-996 %d %.6f %.6f\n%.6f %.6f %.6f %.6f\n%d %.6f %.6f %.6f\n0 0 0 0 0 0 0\n", s.ny, s.dx, s.dy, s.x0, xmax, s.y0, ymax, s.nx,
                 s.rotDeg, s.x0, s.y0);
    for (size_t k = 0; k < s.z.size(); ++k) std::fprintf(f, "%.4f%s", std::isfinite(s.z[k]) ? s.z[k] : 9999900.0, (k % 6 == 5 || k + 1 == s.z.size()) ? "\n" : " ");
    std::fclose(f);
    return true;
}
// ZMAP+: '!' comments, '@name, GRID, nPerLine' / 'width, null, nullstr, decimals, start' / 'rows, cols, xmin, xmax, ymin, ymax'
// / '0.0, 0.0, 0.0' / '@' then values column by column from the top (north) row down
Surface readZmap(const std::string& path, std::string& err) {
    Surface s;
    std::string txt;
    if (!readAll(path, txt)) { err = "cannot open " + path; return s; }
    std::istringstream is(txt);
    std::string line; int stage = 0; std::vector<double> hdr; double nul = 1e30; std::vector<double> vals;
    while (std::getline(is, line)) {
        const std::string t = trim(line);
        if (t.empty() || t[0] == '!') continue;
        if (t[0] == '@') {
            if (stage == 0) { stage = 1; auto f = split(t.substr(1)); if (!f.empty()) s.name = trim(f[0]); continue; }
            stage = 2; continue;
        }
        auto f = split(t);
        if (stage == 1) {   // comma-separated header: empty fields (e.g. the null-value text) keep their position
            std::vector<std::string> fl;
            if (t.find(',') != std::string::npos) {
                std::istringstream cs(t);
                std::string field;
                while (std::getline(cs, field, ',')) fl.push_back(trim(field));
                while (!fl.empty() && fl.back().empty()) fl.pop_back();
            }
            else fl = f;
            for (auto& x : fl) { double d; if (isNum(x, d)) hdr.push_back(d); else hdr.push_back(kNaN); }
        }
        else if (stage == 2) for (auto& x : f) { double d; if (isNum(x, d)) vals.push_back(d); }
    }
    // hdr: width null nullstr decimals start | rows cols xmin xmax ymin ymax | 0 0 0
    std::vector<double> h; for (double d : hdr) h.push_back(d);
    if (h.size() < 11) { err = path + ": incomplete ZMAP+ header"; return s; }
    nul = std::isfinite(h[1]) ? h[1] : (std::isfinite(h[2]) ? h[2] : 1e30);
    const int rows = int(h[5]), cols = int(h[6]);
    const double xmin = h[7], xmax = h[8], ymin = h[9], ymax = h[10];
    if (rows < 2 || cols < 2 || vals.size() < size_t(rows) * size_t(cols)) { err = path + ": ZMAP+ value count mismatch"; return s; }
    s.nx = cols; s.ny = rows; s.x0 = xmin; s.y0 = ymin; s.dx = (xmax - xmin) / (cols - 1); s.dy = (ymax - ymin) / (rows - 1);
    s.z.assign(size_t(rows) * size_t(cols), kNaN);
    size_t k = 0;
    for (int c = 0; c < cols; ++c)
        for (int r = 0; r < rows; ++r) {
            const double d = vals[k++];
            s.z[size_t(c + cols * (rows - 1 - r))] = std::fabs(d - nul) < 1e-6 * (1 + std::fabs(nul)) ? kNaN : d;
        }
    return s;
}
std::vector<P3> readXyz(const std::string& path, std::string& err) {
    std::vector<P3> p;
    std::ifstream f(path);
    if (!f) { err = "cannot open " + path; return p; }
    std::string line;
    while (std::getline(f, line)) {
        const auto t = split(line);
        if (t.size() < 3) continue;
        double x, y, z;
        if (isNum(t[0], x) && isNum(t[1], y) && isNum(t[2], z) && z < 9999899) p.push_back({x, y, std::fabs(z)});
    }
    if (p.empty()) err = path + ": no x y z points";
    return p;
}
Surface gridPoints(const std::vector<P3>& pts, const Surface& tmpl, const geostat::Variogram& vg, int maxNbr) {
    Surface s = tmpl; s.z.assign(size_t(s.nx) * size_t(s.ny), kNaN);
    if (pts.empty()) return s;
    std::vector<geostat::Xyz> tg; tg.reserve(s.z.size());
    for (int j = 0; j < s.ny; ++j) for (int i = 0; i < s.nx; ++i) { double x, y; s.nodeXY(i, j, x, y); tg.push_back({x, y, 0.0}); }
    if (vg.rangeMax > 0) {
        std::vector<geostat::Point> d; d.reserve(pts.size());
        for (const auto& p : pts) d.push_back({p[0], p[1], 0.0, p[2]});
        geostat::Options o; o.vg = vg; o.vg.rangeVert = 1e9; o.maxNbr = maxNbr; o.searchRanges = 3.0;
        const auto K = geostat::krige(tg, d, o);
        for (size_t k = 0; k < tg.size() && k < K.est.size(); ++k) s.z[k] = K.est[k];
    } else {
        for (size_t k = 0; k < tg.size(); ++k) {
            double sw = 0, sz = 0;
            for (const auto& p : pts) {
                const double d2 = (p[0] - tg[k][0]) * (p[0] - tg[k][0]) + (p[1] - tg[k][1]) * (p[1] - tg[k][1]);
                if (d2 < 1e-12) { sw = 1; sz = p[2]; break; }
                sw += 1 / d2; sz += p[2] / d2;
            }
            s.z[k] = sz / sw;
        }
    }
    return s;
}

// ================================================================ LAS 2.0
const LogCurve* WellLog::curve(const std::string& mnem) const {
    const std::string m = upper(mnem);
    for (const auto& c : curves) if (upper(c.mnem) == m) return &c;
    return nullptr;
}
double WellLog::valueAt(const LogCurve& c, double d) const {
    if (md.size() < 2 || c.v.size() != md.size()) return kNaN;
    const bool asc = md.back() > md.front();
    auto it = asc ? std::lower_bound(md.begin(), md.end(), d) : std::lower_bound(md.begin(), md.end(), d, std::greater<double>());
    if (it == md.end()) return std::fabs(d - md.back()) < 1e-9 ? c.v.back() : kNaN;
    const size_t k = size_t(it - md.begin());
    if (k == 0) return std::fabs(d - md[0]) < 1e-9 ? c.v[0] : kNaN;
    const double a = c.v[k - 1], b = c.v[k];
    if (!std::isfinite(a) || !std::isfinite(b)) return std::isfinite(a) ? a : b;
    if (c.discrete) return std::fabs(d - md[k - 1]) <= std::fabs(md[k] - d) ? a : b;
    const double f = (d - md[k - 1]) / (md[k] - md[k - 1]);
    return a + f * (b - a);
}
WellLog readLas(const std::string& path, std::string& err) {
    WellLog L; L.file = path;
    std::ifstream f(path);
    if (!f) { err = "cannot open " + path; return L; }
    std::string line, sec; bool wrap = false;
    std::vector<LogCurve> all;
    std::vector<double> data;
    while (std::getline(f, line)) {
        const std::string t = trim(line);
        if (t.empty() || t[0] == '#') continue;
        if (t[0] == '~') { sec = upper(t.substr(1, 1)); continue; }
        if (sec == "A") { for (const auto& x : split(t)) { double d; if (isNum(x, d)) data.push_back(d); } continue; }
        // header line: MNEM.UNIT  DATA : DESCRIPTION
        const size_t dot = t.find('.');
        if (dot == std::string::npos) continue;
        const std::string mn = trim(t.substr(0, dot));
        const size_t sp = t.find_first_of(" \t", dot);
        const std::string unit = sp == std::string::npos ? t.substr(dot + 1) : t.substr(dot + 1, sp - dot - 1);
        const size_t colon = t.rfind(':');
        const std::string val = trim(t.substr(sp == std::string::npos ? t.size() : sp, colon == std::string::npos || colon < sp ? std::string::npos : colon - sp));
        const std::string desc = colon == std::string::npos ? std::string() : trim(t.substr(colon + 1));
        const std::string M = upper(mn);
        if (sec == "V" && M == "WRAP") wrap = upper(val).rfind("YES", 0) == 0;
        if (sec == "W") {
            if (M == "NULL") { double d; if (isNum(val, d)) L.nullValue = d; }
            if (M == "WELL") L.well = val.empty() ? desc : val;
        }
        if (sec == "C") { LogCurve c; c.mnem = mn; c.unit = unit; c.descr = desc; all.push_back(c); }
    }
    (void)wrap;   // wrapped and unwrapped data are read the same way (values in curve order)
    if (all.empty()) { err = path + ": no ~C curve section"; return L; }
    const size_t nc = all.size(), nr = data.size() / nc;
    if (nr == 0) { err = path + ": no ~A data"; return L; }
    for (size_t c = 0; c < nc; ++c) all[c].v.resize(nr);
    for (size_t r = 0; r < nr; ++r)
        for (size_t c = 0; c < nc; ++c) {
            const double d = data[r * nc + c];
            all[c].v[r] = std::fabs(d - L.nullValue) < 1e-6 * (1 + std::fabs(L.nullValue)) ? kNaN : d;
        }
    L.md = all[0].v;
    for (size_t c = 1; c < nc; ++c) {
        LogCurve& cv = all[c];
        const std::string M = upper(cv.mnem);
        cv.discrete = M.find("FAC") != std::string::npos || M.find("LITH") != std::string::npos || M == "ZONE";
        L.curves.push_back(cv);
    }
    if (L.well.empty()) { const size_t sl = path.find_last_of("/\\"); L.well = path.substr(sl == std::string::npos ? 0 : sl + 1); const size_t d = L.well.find('.'); if (d != std::string::npos) L.well = L.well.substr(0, d); }
    return L;
}
bool writeLas(const std::string& path, const WellLog& L) {
    FILE* f = std::fopen(path.c_str(), "w");
    if (!f || L.md.empty()) { if (f) std::fclose(f); return false; }
    std::fprintf(f, "~VERSION INFORMATION\n VERS.   2.0 : CWLS LOG ASCII STANDARD - VERSION 2.0\n WRAP.   NO  : ONE LINE PER DEPTH STEP\n");
    std::fprintf(f, "~WELL INFORMATION\n STRT.M  %.4f : START DEPTH\n STOP.M  %.4f : STOP DEPTH\n STEP.M  %.4f : STEP\n NULL.   %.2f : NULL VALUE\n WELL.   %s : WELL\n",
                 L.md.front(), L.md.back(), L.md.size() > 1 ? L.md[1] - L.md[0] : 0.0, L.nullValue, L.well.c_str());
    std::fprintf(f, "~CURVE INFORMATION\n DEPT.M : MEASURED DEPTH\n");
    for (const auto& c : L.curves) std::fprintf(f, " %s.%s : %s\n", c.mnem.c_str(), c.unit.c_str(), c.descr.c_str());
    std::fprintf(f, "~A\n");
    for (size_t r = 0; r < L.md.size(); ++r) {
        std::fprintf(f, "%10.3f", L.md[r]);
        for (const auto& c : L.curves) std::fprintf(f, " %12.5f", r < c.v.size() && std::isfinite(c.v[r]) ? c.v[r] : L.nullValue);
        std::fprintf(f, "\n");
    }
    std::fclose(f);
    return true;
}

// ================================================================ wells
void minimumCurvature(GeoWell& w) {
    auto& T = w.traj;
    if (T.empty()) return;
    std::sort(T.begin(), T.end(), [](const Station& a, const Station& b) { return a.md < b.md; });
    if (T[0].md > 0) T.insert(T.begin(), Station{0, T[0].inc * 0, T[0].azi, 0, 0, 0});
    T[0].x = w.x; T[0].y = w.y; T[0].tvd = -w.kb + T[0].md;   // TVDSS: depth below datum = MD - KB for the vertical top part
    for (size_t k = 1; k < T.size(); ++k) {
        const Station &a = T[k - 1];
        Station& b = T[k];
        const double i1 = a.inc * kDeg, i2 = b.inc * kDeg, a1 = a.azi * kDeg, a2 = b.azi * kDeg, dmd = b.md - a.md;
        double cosd = std::cos(i2 - i1) - std::sin(i1) * std::sin(i2) * (1 - std::cos(a2 - a1));
        cosd = std::max(-1.0, std::min(1.0, cosd));
        const double dl = std::acos(cosd);
        const double rf = dl > 1e-9 ? 2 / dl * std::tan(dl / 2) : 1.0;
        b.x = a.x + dmd / 2 * (std::sin(i1) * std::sin(a1) + std::sin(i2) * std::sin(a2)) * rf;   // east
        b.y = a.y + dmd / 2 * (std::sin(i1) * std::cos(a1) + std::sin(i2) * std::cos(a2)) * rf;   // north
        b.tvd = a.tvd + dmd / 2 * (std::cos(i1) + std::cos(i2)) * rf;
    }
}
P3 GeoWell::at(double md) const {
    if (traj.size() < 2) return {x, y, md - kb};
    if (md <= traj.front().md) return {traj.front().x, traj.front().y, traj.front().tvd - (traj.front().md - md)};
    for (size_t k = 1; k < traj.size(); ++k)
        if (md <= traj[k].md) {
            const Station &a = traj[k - 1], &b = traj[k];
            const double f = (md - a.md) / std::max(1e-12, b.md - a.md);
            return {a.x + f * (b.x - a.x), a.y + f * (b.y - a.y), a.tvd + f * (b.tvd - a.tvd)};
        }
    const Station& b = traj.back();
    const double ex = md - b.md, i = b.inc * kDeg, az = b.azi * kDeg;   // straight extrapolation
    return {b.x + ex * std::sin(i) * std::sin(az), b.y + ex * std::sin(i) * std::cos(az), b.tvd + ex * std::cos(i)};
}
double GeoWell::totalMd() const {
    double m = traj.empty() ? 0 : traj.back().md;
    for (const auto& l : logs) if (!l.md.empty()) m = std::max(m, std::max(l.md.front(), l.md.back()));
    return m;
}
double GeoWell::mdAtTvd(double tvd) const {
    const double L = std::max(totalMd(), 1.0);
    double prev = at(0)[2];
    const int n = 4000;
    for (int k = 1; k <= n; ++k) {
        const double m = L * k / n, z = at(m)[2];
        if ((prev - tvd) * (z - tvd) <= 0 && z != prev) return L * (k - 1) / n + (tvd - prev) / (z - prev) * (L / n);
        prev = z;
    }
    return kNaN;
}
std::vector<GeoWell> readWellHeads(const std::string& path, std::string& err) {
    std::vector<GeoWell> w;
    std::ifstream f(path);
    if (!f) { err = "cannot open " + path; return w; }
    std::string line;
    while (std::getline(f, line)) {
        const std::string t = trim(line);
        if (t.empty() || t[0] == '#' || t[0] == '!') continue;
        const auto x = split(t);
        double a, b, c = 0;
        if (x.size() >= 3 && isNum(x[1], a) && isNum(x[2], b)) {
            GeoWell g; g.name = x[0]; g.x = a; g.y = b;
            if (x.size() >= 4 && isNum(x[3], c)) g.kb = c;
            w.push_back(g);
        }
    }
    if (w.empty()) err = path + ": no 'name x y [kb]' rows";
    return w;
}
bool readDeviation(const std::string& path, GeoWell& w, std::string& err) {
    std::ifstream f(path);
    if (!f) { err = "cannot open " + path; return false; }
    std::string line; std::vector<std::string> cols;
    std::vector<std::vector<double>> rows;
    while (std::getline(f, line)) {
        const std::string t = trim(line);
        if (t.empty()) continue;
        if (t[0] == '#') {
            const std::string U = upper(t);
            const size_t c = t.find(':');
            if (U.find("WELL NAME") != std::string::npos && c != std::string::npos) w.name = trim(t.substr(c + 1));
            double d;
            if (U.find("HEAD X") != std::string::npos && c != std::string::npos && isNum(trim(t.substr(c + 1)), d)) w.x = d;
            if (U.find("HEAD Y") != std::string::npos && c != std::string::npos && isNum(trim(t.substr(c + 1)), d)) w.y = d;
            if (U.find("KB") != std::string::npos && c != std::string::npos && isNum(trim(t.substr(c + 1)), d)) w.kb = d;
            continue;
        }
        if (t.rfind("====", 0) == 0) continue;
        const auto x = split(t);
        double d;
        if (!x.empty() && !isNum(x[0], d)) { cols.clear(); for (auto& s : x) cols.push_back(upper(s)); continue; }
        std::vector<double> r; for (auto& s : x) { if (isNum(s, d)) r.push_back(d); }
        if (r.size() >= 3) rows.push_back(r);
    }
    if (rows.empty()) { err = path + ": no survey rows"; return false; }
    auto col = [&](const char* n) { for (size_t i = 0; i < cols.size(); ++i) if (cols[i] == n) return int(i); return -1; };
    int iMd = col("MD"), iInc = col("INCL"), iAz = col("AZIM");
    if (iInc < 0) iInc = col("INC");
    if (iAz < 0) iAz = col("AZI");
    if (iMd < 0) { iMd = 0; iInc = 1; iAz = 2; }
    if (iInc < 0 || iAz < 0) {   // Petrel file without angles: use X Y TVD
        const int iX = col("X"), iY = col("Y"), iT = col("TVD");
        if (iX < 0 || iY < 0 || iT < 0) { err = path + ": need MD INCL AZIM or MD X Y TVD columns"; return false; }
        w.traj.clear();
        for (auto& r : rows) { Station s; s.md = r[size_t(iMd)]; s.x = r[size_t(iX)]; s.y = r[size_t(iY)]; s.tvd = r[size_t(iT)] - w.kb; w.traj.push_back(s); }
        if (!w.traj.empty()) { w.x = w.traj[0].x; w.y = w.traj[0].y; }
        return true;
    }
    w.traj.clear();
    for (auto& r : rows) {
        if (size_t(std::max(iMd, std::max(iInc, iAz))) >= r.size()) continue;
        Station s; s.md = r[size_t(iMd)]; s.inc = r[size_t(iInc)]; s.azi = r[size_t(iAz)]; w.traj.push_back(s);
    }
    minimumCurvature(w);
    return true;
}
std::vector<std::pair<std::string, Top>> readTops(const std::string& path, std::string& err) {
    std::vector<std::pair<std::string, Top>> out;
    std::ifstream f(path);
    if (!f) { err = "cannot open " + path; return out; }
    std::string line;
    while (std::getline(f, line)) {
        const std::string t = trim(line);
        if (t.empty() || t[0] == '#') continue;
        const auto x = split(t);
        double d;
        if (x.size() >= 3 && isNum(x[2], d)) { Top tp; tp.surface = x[1]; tp.md = d; out.push_back({x[0], tp}); }
    }
    if (out.empty()) err = path + ": no 'well surface md' rows";
    return out;
}

// ================================================================ faults
std::vector<std::array<double, 2>> Fault::traceAt(double depth) const {
    std::vector<std::array<double, 2>> tr;
    for (const auto& s : sticks) {
        if (s.empty()) continue;
        if (s.size() == 1) { tr.push_back({s[0][0], s[0][1]}); continue; }
        // interpolate / extrapolate along the stick at the depth
        size_t k = 1;
        while (k + 1 < s.size() && s[k][2] < depth) ++k;
        const P3 &a = s[k - 1], &b = s[k];
        const double f = std::fabs(b[2] - a[2]) > 1e-9 ? (depth - a[2]) / (b[2] - a[2]) : 0.0;
        tr.push_back({a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1])});
    }
    return tr;
}
std::vector<Fault> readFaultSticks(const std::string& path, std::string& err) {
    std::vector<Fault> out;
    std::ifstream f(path);
    if (!f) { err = "cannot open " + path; return out; }
    std::map<std::string, std::map<int, std::vector<P3>>> by;
    std::vector<std::string> order;
    std::string line; int autoStick = 0; bool newStick = true;
    const size_t sl = path.find_last_of("/\\");
    std::string defName = path.substr(sl == std::string::npos ? 0 : sl + 1);
    if (defName.find('.') != std::string::npos) defName = defName.substr(0, defName.find('.'));
    while (std::getline(f, line)) {
        const std::string t = trim(line);
        if (t.empty()) { newStick = true; continue; }
        if (t[0] == '#' || t[0] == '!') continue;
        const auto x = split(t);
        double a, b, c, d;
        if (x.size() < 3 || !isNum(x[0], a) || !isNum(x[1], b) || !isNum(x[2], c)) continue;
        if (std::fabs(a + 999) < 1e-6 || std::fabs(c + 999) < 1e-6 || a >= 999999) { newStick = true; continue; }
        int stick; std::string fault = defName;
        if (x.size() >= 4 && isNum(x[3], d)) { stick = int(d); if (x.size() >= 5) fault = x[4]; }
        else { if (newStick) ++autoStick; stick = autoStick; }
        newStick = false;
        if (!by.count(fault)) order.push_back(fault);
        by[fault][stick].push_back({a, b, std::fabs(c)});
    }
    for (const auto& n : order) {
        Fault F; F.name = n;
        for (auto& [id, pts] : by[n]) { std::sort(pts.begin(), pts.end(), [](const P3& p, const P3& q) { return p[2] < q[2]; }); F.sticks.push_back(pts); }
        out.push_back(F);
    }
    if (out.empty()) err = path + ": no fault sticks";
    return out;
}

}  // namespace geo
}  // namespace ressim
