// Publication-quality figure export: see figure.h.
#include "figure.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>

namespace fig {

// ===================================================================== palettes
const std::vector<Palette>& palettes() {
    static const std::vector<Palette> p = {
        {"Tableau 10", {rgba(31, 119, 180), rgba(255, 127, 14), rgba(44, 160, 44), rgba(214, 39, 40), rgba(148, 103, 189),
                        rgba(140, 86, 75), rgba(227, 119, 194), rgba(127, 127, 127), rgba(188, 189, 34), rgba(23, 190, 207)}},
        {"Okabe-Ito (colour-blind safe)", {rgba(0, 114, 178), rgba(230, 159, 0), rgba(0, 158, 115), rgba(213, 94, 0), rgba(86, 180, 233),
                                           rgba(204, 121, 167), rgba(240, 228, 66), rgba(0, 0, 0)}},
        {"Paul Tol bright", {rgba(68, 119, 170), rgba(238, 102, 119), rgba(34, 136, 51), rgba(204, 187, 68), rgba(102, 204, 238),
                             rgba(170, 51, 119), rgba(187, 187, 187)}},
        {"Paul Tol muted", {rgba(51, 34, 136), rgba(136, 204, 238), rgba(68, 170, 153), rgba(17, 119, 51), rgba(153, 153, 51),
                            rgba(221, 204, 119), rgba(204, 102, 119), rgba(136, 34, 85), rgba(170, 68, 153)}},
        {"Dark 2", {rgba(27, 158, 119), rgba(217, 95, 2), rgba(117, 112, 179), rgba(231, 41, 138), rgba(102, 166, 30), rgba(230, 171, 2),
                    rgba(166, 118, 29), rgba(102, 102, 102)}},
        {"Grayscale", {rgba(0, 0, 0), rgba(90, 90, 90), rgba(150, 150, 150), rgba(40, 40, 40), rgba(120, 120, 120), rgba(180, 180, 180)}},
    };
    return p;
}

static uint32_t hsv(double h, double s, double v) {
    h = h - std::floor(h);
    const double i = std::floor(h * 6), f = h * 6 - i, p = v * (1 - s), q = v * (1 - f * s), t = v * (1 - (1 - f) * s);
    double r, g, b;
    switch (int(i) % 6) {
    case 0: r = v; g = t; b = p; break;
    case 1: r = q; g = v; b = p; break;
    case 2: r = p; g = v; b = t; break;
    case 3: r = p; g = q; b = v; break;
    case 4: r = t; g = p; b = v; break;
    default: r = v; g = p; b = q; break;
    }
    return rgba(int(r * 255 + 0.5), int(g * 255 + 0.5), int(b * 255 + 0.5));
}

uint32_t distinctColor(int i, int palette) {
    const auto& P = palettes();
    const auto& c = P[size_t(std::max(0, std::min(palette, int(P.size()) - 1)))].colors;
    if (i >= 0 && i < int(c.size())) return c[size_t(i)];
    const int k = std::max(0, i - int(c.size()));
    // golden-angle hues, alternating value levels so neighbours differ in lightness too
    return hsv(0.07 + k * 0.618033988749895, 0.70 - 0.15 * (k % 2), 0.92 - 0.22 * ((k / 2) % 2));
}

// ===================================================================== font metrics
static const short kHelv[95] = {278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278, 333, 278, 278, 556, 556, 556, 556,
                                556, 556, 556, 556, 556, 556, 278, 278, 584, 584, 584, 556, 1015, 667, 667, 722, 722, 667, 611, 778,
                                722, 278, 500, 667, 556, 833, 722, 778, 667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, 278,
                                278, 278, 469, 556, 333, 556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500, 222, 833, 556, 556,
                                556, 556, 333, 500, 278, 556, 500, 722, 500, 500, 500, 334, 260, 334, 584};
static const short kTimes[95] = {250, 333, 408, 500, 500, 833, 778, 180, 333, 333, 500, 564, 250, 333, 250, 278, 500, 500, 500, 500,
                                 500, 500, 500, 500, 500, 500, 278, 278, 564, 564, 564, 444, 921, 722, 667, 667, 722, 611, 556, 722,
                                 722, 333, 389, 722, 611, 889, 722, 722, 556, 722, 667, 556, 611, 722, 722, 944, 722, 722, 611, 333,
                                 278, 333, 469, 500, 333, 444, 500, 444, 500, 444, 333, 500, 500, 278, 278, 500, 278, 778, 500, 500,
                                 500, 500, 333, 389, 278, 500, 500, 722, 500, 500, 444, 480, 200, 480, 541};

namespace {
struct Seg { std::string text; bool sup; };
std::vector<Seg> segments(const std::string& s) {
    std::vector<Seg> out;
    size_t i = 0;
    std::string cur;
    while (i < s.size()) {
        if (s[i] == '^' && i + 1 < s.size() && s[i + 1] == '{') {
            const size_t e = s.find('}', i + 2);
            if (e != std::string::npos) {
                if (!cur.empty()) { out.push_back({cur, false}); cur.clear(); }
                out.push_back({s.substr(i + 2, e - i - 2), true});
                i = e + 1;
                continue;
            }
        }
        cur += s[i++];
    }
    if (!cur.empty()) out.push_back({cur, false});
    return out;
}
// UTF-8 -> code points
std::vector<unsigned> decode(const std::string& s) {
    std::vector<unsigned> cp;
    for (size_t i = 0; i < s.size();) {
        const unsigned char c = (unsigned char)s[i];
        unsigned u = c; int n = 0;
        if (c >= 0xF0) { u = c & 7; n = 3; } else if (c >= 0xE0) { u = c & 15; n = 2; } else if (c >= 0xC0) { u = c & 31; n = 1; }
        ++i;
        for (int k = 0; k < n && i < s.size(); ++k, ++i) u = (u << 6) | ((unsigned char)s[i] & 63);
        cp.push_back(u);
    }
    return cp;
}
double rawWidth(const std::string& s, int font) {
    const short* w = font == 1 ? kTimes : kHelv;
    double t = 0;
    for (unsigned u : decode(s)) t += (u >= 32 && u < 127) ? w[u - 32] : (u == 0xD7 ? (font ? 564 : 584) : (font ? 500 : 556));
    return t / 1000.0;
}
const double kSup = 0.68;   // superscript size factor
const double kRise = 0.38;  // superscript baseline rise (fraction of font size)
} // namespace

double textWidthPt(const std::string& s, double size, int font) {
    double w = 0;
    for (const auto& g : segments(s)) w += rawWidth(g.text, font) * size * (g.sup ? kSup : 1.0);
    return w;
}
double Canvas::textWidth(const std::string& s, double sizePt) const { return textWidthPt(s, sizePt, font); }
void Canvas::fillTriangles(const double* xy, int nTri, uint32_t color) {
    // fallback: outline each triangle (canvases override this with a real fill)
    for (int t = 0; t < nTri; ++t) {
        const double* q = xy + 6 * t;
        const double pts[8] = {q[0], q[1], q[2], q[3], q[4], q[5], q[0], q[1]};
        polyline(pts, 4, color, 0.1, 0);
    }
}
void Canvas::textFit(double x, double y, const std::string& s, double sizePt, double, bool vertical, uint32_t color) {
    text(x, y, s, sizePt, 0, vertical, color);
}
uint32_t printColor(uint32_t c, const Snapshot* sn) {
    if (sn)
        for (size_t i = 0; i < sn->printFrom.size() && i < sn->printTo.size(); ++i)
            if (sn->printFrom[i] == c) return sn->printTo[i];
    const int r = int(c & 255), g = int((c >> 8) & 255), b = int((c >> 16) & 255);
    const int mx = std::max(r, std::max(g, b)), mn = std::min(r, std::min(g, b));
    if (mx - mn > 40) return c;   // a real colour (series, colormap) stays
    if (mn > 35 && mx < 215) return c;   // mid greys may be data (grey colour-map cells): kept
    return (c & 0xFF000000u) | uint32_t(255 - r) | (uint32_t(255 - g) << 8) | (uint32_t(255 - b) << 16);
}

// ===================================================================== ticks
std::vector<double> niceTicks(double a, double b, int target, double* stepOut) {
    std::vector<double> t;
    if (!(b > a) || !std::isfinite(a) || !std::isfinite(b)) { if (stepOut) *stepOut = 1; return t; }
    const double raw = (b - a) / std::max(1, target);
    const double mag = std::pow(10.0, std::floor(std::log10(raw)));
    const double r = raw / mag;
    const double step = (r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10) * mag;
    const double first = std::ceil(a / step - 1e-9) * step;
    for (double v = first; v <= b + step * 1e-9; v += step) t.push_back(std::fabs(v) < step * 1e-9 ? 0.0 : v);
    if (stepOut) *stepOut = step;
    return t;
}
int axisExponent(double a, double b) {
    const double m = std::max(std::fabs(a), std::fabs(b));
    if (!(m > 0) || !std::isfinite(m)) return 0;
    if (m >= 1e5 || m < 1e-3) return int(std::floor(std::log10(m) / 3.0)) * 3;
    return 0;
}
std::string formatTick(double v, double step, int e) {
    const double sc = std::pow(10.0, -e);
    const double sv = v * sc, ss = step * sc;
    int dec = 0;
    if (ss > 0) dec = std::max(0, -int(std::floor(std::log10(ss) + 1e-9)));
    // step like 2.5 needs one more decimal
    if (ss > 0 && std::fabs(ss * std::pow(10.0, dec) - std::round(ss * std::pow(10.0, dec))) > 1e-6) ++dec;
    dec = std::min(dec, 8);
    char b[64];
    std::snprintf(b, sizeof(b), "%.*f", dec, sv);
    std::string s = b;
    if (s == "-0" || s.find_first_not_of("-0.") == std::string::npos) s = "0";
    // typographic minus
    if (!s.empty() && s[0] == '-') s = "\xE2\x88\x92" + s.substr(1);
    return s;
}

// ===================================================================== layout & render
namespace {

struct Range { double lo = 0, hi = 1; };
void dataRange(const Panel& p, int axis, bool x, Range& r, bool& any, bool logAx = false) {
    any = false;
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    auto take = [&](double d) { if (std::isfinite(d) && (!logAx || d > 0)) { lo = std::min(lo, d); hi = std::max(hi, d); any = true; } };
    for (const auto& s : p.series) {
        if (!x && s.axis != axis) continue;
        const auto& v = x ? s.x : s.y;
        for (double d : v) take(d);
        if (!x) for (double d : s.ylo) take(d);
        if (s.barWidth > 0) {
            if (x) for (double d : s.x) { take(d - 0.5 * s.barWidth); take(d + 0.5 * s.barWidth); }
            else take(0.0);
        }
    }
    if (axis == 0) {   // C23 layers live on the primary axes
        for (const auto& h : p.heatmaps) { take(x ? h.x0 : h.y0); take(x ? h.x1 : h.y1); }
        for (const auto& sh : p.shapes)
            if (sh.kind == Shape::Polyline || sh.kind == Shape::Polygon) for (double d : (x ? sh.x : sh.y)) take(d);
    }
    if (!any) { r = {0, 1}; return; }
    if (logAx) {   // transformed (log10) range, padded to whole decades when it spans several
        double a = std::log10(lo), b = std::log10(hi);
        if (b - a < 1e-9) { a -= 0.5; b += 0.5; }
        else if (b - a >= 2) { a = std::floor(a + 1e-9); b = std::ceil(b - 1e-9); }
        else { const double pad = (b - a) * 0.05; a -= pad; b += pad; }
        r = {a, b};
        return;
    }
    if (hi - lo < 1e-12 * std::max(1.0, std::fabs(hi))) {
        const double d = std::fabs(hi) > 0 ? std::fabs(hi) * 0.1 : 1.0;
        lo -= d; hi += d;
    } else if (!x && p.heatmaps.empty()) {
        const double pad = (hi - lo) * 0.05;
        // keep zero-based quantities anchored at zero
        lo = (lo >= 0 && lo - pad < 0) ? 0.0 : lo - pad;
        hi += pad;
    }
    r = {lo, hi};
}

struct Axis {
    bool log = false;   // C17: r, ticks in log10 space
    bool inv = false;   // C23: inverted axis (values grow downwards / to the left)
    double tf(double v) const { return log ? (v > 0 ? std::log10(v) : std::numeric_limits<double>::quiet_NaN()) : v; }
    Range r;
    std::vector<double> ticks;
    double step = 1;
    int exp = 0;
    std::vector<std::string> labels;
    double labelW = 0;
};
std::string logLabel(double t) {   // tick label of a log axis at log10 value t
    const double n = std::round(t);
    if (std::fabs(t - n) < 1e-9) {
        if (n >= -2 && n <= 3) { char b[32]; std::snprintf(b, sizeof b, "%g", std::pow(10.0, n)); return b; }
        return "10^{" + std::string(n < 0 ? "\xE2\x88\x92" : "") + std::to_string(int(std::fabs(n))) + "}";
    }
    char b[32]; std::snprintf(b, sizeof b, "%g", std::pow(10.0, t));
    return b;
}
void makeAxis(Axis& a, int target, double fs, int font) {
    if (a.log) {
        a.ticks.clear(); a.exp = 0; a.step = 1; a.labels.clear(); a.labelW = 0;
        const int d0 = int(std::ceil(a.r.lo - 1e-9)), d1 = int(std::floor(a.r.hi + 1e-9));
        const int nd = d1 - d0 + 1;
        const int every = std::max(1, (nd + std::max(1, target) - 1) / std::max(1, target));
        for (int d = d0; d <= d1; ++d) if ((d - d0) % every == 0) a.ticks.push_back(double(d));
        if (nd < 3)   // short range: 1-2-5 sequence
            for (int d = d0 - 1; d <= d1; ++d)
                for (double m : {2.0, 5.0}) { const double t = d + std::log10(m); if (t >= a.r.lo - 1e-9 && t <= a.r.hi + 1e-9) a.ticks.push_back(t); }
        std::sort(a.ticks.begin(), a.ticks.end());
        for (double t : a.ticks) { a.labels.push_back(logLabel(t)); a.labelW = std::max(a.labelW, textWidthPt(a.labels.back(), fs, font)); }
        return;
    }
    a.ticks = niceTicks(a.r.lo, a.r.hi, target, &a.step);
    double m = 0;
    for (double t : a.ticks) m = std::max(m, std::fabs(t));
    a.exp = axisExponent(0, m);
    a.labels.clear();
    a.labelW = 0;
    for (double t : a.ticks) {
        a.labels.push_back(formatTick(t, a.step, a.exp));
        a.labelW = std::max(a.labelW, textWidthPt(a.labels.back(), fs, font));
    }
}
std::string expLabel(int e) { return e ? "\xC3\x97" "10^{" + std::to_string(e) + "}" : std::string(); }

struct LegendItem { std::string label; uint32_t color; double width; bool scatter, markers; int dash; int swatch = 0; int marker = 0; uint32_t stroke = 0; };
// C23: light neutral colours made for the dark theme (white / light grey lines and text) print dark on paper
uint32_t ink(uint32_t c, bool keep = false) {
    if (keep) return c;
    const int r = int(c & 255), g = int(c >> 8 & 255), b = int(c >> 16 & 255);
    const int mx = std::max(r, std::max(g, b)), mn = std::min(r, std::min(g, b));
    if (mx - mn < 40 && (r * 299 + g * 587 + b * 114) / 1000 > 150) {
        const uint32_t a = std::max<uint32_t>(c >> 24, 90u);
        return (rgba(45, 45, 45) & 0x00FFFFFFu) | (a << 24);
    }
    return c;
}
std::vector<LegendItem> legendItems(const Panel& p) {
    std::vector<LegendItem> v;
    std::set<std::string> seen;
    for (const auto& s : p.series) {
        if (!s.legend || s.label.empty() || s.label.compare(0, 2, "##") == 0) continue;
        if (!seen.insert(s.label).second) continue;
        LegendItem it{s.label, ink(s.color), s.width, s.scatter, s.markers, s.dash};
        if (!s.ylo.empty() || s.barWidth > 0) { it.swatch = 1; it.stroke = s.barLine; }
        v.push_back(it);
    }
    for (const auto& sh : p.shapes) {
        if (sh.label.empty() || sh.label.compare(0, 2, "##") == 0 || !seen.insert(sh.label).second) continue;
        LegendItem it{sh.label, ink(sh.kind == Shape::Polygon || sh.kind == Shape::Marker ? (sh.fill ? sh.fill : sh.stroke) : sh.stroke, sh.keepColor),
                      std::max(0.4, sh.widthPx / 1.6), false, false, sh.dash};
        it.swatch = sh.kind == Shape::Polygon ? 1 : sh.kind == Shape::Marker ? 2 : 0;
        it.marker = sh.marker;
        it.stroke = ink(sh.stroke, sh.keepColor);
        v.push_back(it);
    }
    return v;
}
// C23 marker shapes (radius r in pt)
void drawMarker(Canvas& c, double x, double y, double r, int shape, uint32_t fill, uint32_t stroke, double lw) {
    const uint32_t sc = stroke ? stroke : fill;
    switch (shape) {
    case 1: if (fill) c.fillRect(x - r, y - r, 2 * r, 2 * r, fill); c.strokeRect(x - r, y - r, 2 * r, 2 * r, sc, lw); break;
    case 2: case 3: {
        double t[6];
        double pts[10];
        int np;
        if (shape == 2) { const double q[8] = {x, y - r * 1.25, x - r, y, x, y + r * 1.25, x + r, y}; std::copy(q, q + 8, pts); np = 4; }
        else { const double q[6] = {x, y - r * 1.2, x - r * 1.1, y + r * 0.8, x + r * 1.1, y + r * 0.8}; std::copy(q, q + 6, pts); np = 3; }
        if (fill) {
            for (int k = 1; k + 1 < np; ++k) {   // fan, counter-clockwise on the page (y down)
                t[0] = pts[0]; t[1] = pts[1]; t[2] = pts[2 * k]; t[3] = pts[2 * k + 1]; t[4] = pts[2 * k + 2]; t[5] = pts[2 * k + 3];
                if ((t[2] - t[0]) * (t[5] - t[1]) - (t[3] - t[1]) * (t[4] - t[0]) < 0) { std::swap(t[2], t[4]); std::swap(t[3], t[5]); }
                c.fillTriangles(t, 1, fill);
            }
        }
        pts[2 * np] = pts[0]; pts[2 * np + 1] = pts[1];
        c.polyline(pts, np + 1, sc, lw, 0);
        break;
    }
    case 4: { const double a[4] = {x - r, y - r, x + r, y + r}, b[4] = {x - r, y + r, x + r, y - r}; c.polyline(a, 2, sc, lw, 0); c.polyline(b, 2, sc, lw, 0); break; }
    case 5: { const double a[4] = {x - r, y, x + r, y}, b[4] = {x, y - r, x, y + r}; c.polyline(a, 2, sc, lw, 0); c.polyline(b, 2, sc, lw, 0); break; }
    default: c.circle(x, y, r, fill, sc, lw); break;
    }
}
double dashLen(int dash, double w, double* gap) {
    switch (dash) {
    case 1: *gap = 2.2 * std::max(w, 0.6); return 4.0 * std::max(w, 0.6);
    case 2: *gap = 1.6 * std::max(w, 0.6); return 0.8 * std::max(w, 0.6);
    case 3: *gap = 1.8 * std::max(w, 0.6); return 5.0 * std::max(w, 0.6);
    default: *gap = 0; return 0;
    }
}

struct Ctx {
    const Style& st;
    Canvas& c;
    double fs, tfs, lineH, tick, pad;
    int font;
};

void drawLegendBox(Ctx& X, const std::vector<LegendItem>& items, double x, double y, double maxW, bool box, double* outW = nullptr, double* outH = nullptr,
                   bool measureOnly = false, int columns = 1) {
    const double fs = X.tfs, sample = fs * 2.2, gap = fs * 0.45, rowH = fs * 1.25, padIn = fs * 0.45;
    // column widths
    double itemW = 0;
    for (const auto& it : items) itemW = std::max(itemW, sample + gap + textWidthPt(it.label, fs, X.font));
    columns = std::max(1, columns);
    const int rows = (int(items.size()) + columns - 1) / columns;
    const double colGap = fs * 1.2;
    const double w = std::min(maxW, padIn * 2 + columns * itemW + (columns - 1) * colGap), h = padIn * 2 + rows * rowH - (rowH - fs) * 0.3;
    if (outW) *outW = w;
    if (outH) *outH = h;
    if (measureOnly) return;
    if (box) {
        X.c.fillRect(x, y, w, h, rgba(255, 255, 255, 230));
        X.c.strokeRect(x, y, w, h, rgba(160, 160, 160), 0.4);
    }
    for (size_t i = 0; i < items.size(); ++i) {
        const int col = int(i) / rows, row = int(i) % rows;
        const double ix = x + padIn + col * (itemW + colGap), iy = y + padIn + row * rowH + fs * 0.5;
        const auto& it = items[i];
        const double lw = X.st.lineWidthPt * it.width;
        if (it.swatch == 1) {   // filled area / bars
            X.c.fillRect(ix, iy - fs * 0.38, sample, fs * 0.76, it.color);
            if (it.stroke) X.c.strokeRect(ix, iy - fs * 0.38, sample, fs * 0.76, it.stroke, 0.4);
        } else if (it.swatch == 2) {
            drawMarker(X.c, ix + sample * 0.5, iy, fs * 0.32, it.marker, it.color, it.stroke, 0.5);
        } else {
            if (!it.scatter) {
                const double pts[4] = {ix, iy, ix + sample, iy};
                X.c.polyline(pts, 2, it.color, lw, it.dash);
            }
            if (it.scatter || it.markers) X.c.circle(ix + sample * 0.5, iy, std::max(1.4, lw * 1.3), it.color, it.color, 0.3);
        }
        X.c.text(ix + sample + gap, iy + fs * 0.34, it.label, fs, 0, false, rgba(0, 0, 0));
    }
}

struct PanelLayout {
    Axis x, y, y2, cb; double l = 0, r = 0, t = 0, b = 0; bool hasY2 = false;
    ColorBar bar; bool hasCb = false; double cbW = 0;   // C23 colour bar (explicit or from the first heatmap)
};

void computeAxes(Ctx& X, const Panel& p, PanelLayout& L, double cellW, double cellH) {
    Range xr, yr, y2r; bool any;
    L.x.log = p.logX; L.y.log = p.logY; L.y2.log = p.logY2;
    dataRange(p, 0, true, xr, any, p.logX);
    dataRange(p, 0, false, yr, any, p.logY);
    bool any2 = false;
    dataRange(p, 1, false, y2r, any2, p.logY2);
    L.hasY2 = any2 || (p.hasY2 && p.hasLimits && X.st.useLimits && any2);
    if (X.st.useLimits && p.hasLimits) {   // visible app ranges (log axes: only a positive range)
        auto lim = [](bool lg, double a, double b, Range& r) {
            if (!(b > a)) return;
            if (!lg) { r = {a, b}; return; }
            if (a > 0) r = {std::log10(a), std::log10(b)};
        };
        lim(p.logX, p.x0, p.x1, xr);
        lim(p.logY, p.y0, p.y1, yr);
        if (p.hasY2) lim(p.logY2, p.y20, p.y21, y2r);
    }
    L.x.inv = p.invertX; L.y.inv = p.invertY; L.y2.inv = p.invertY;
    if (p.equal && !p.logX && !p.logY && cellW > 0 && cellH > 0) {   // C23: equal scaling - widen the short axis around its centre
        const double dx = xr.hi - xr.lo, dy = yr.hi - yr.lo, ratio = cellW / cellH;
        if (dx > 0 && dy > 0) {
            if (dx / dy < ratio) { const double c0 = 0.5 * (xr.lo + xr.hi), h = 0.5 * dy * ratio; xr = {c0 - h, c0 + h}; }
            else { const double c0 = 0.5 * (yr.lo + yr.hi), h = 0.5 * dx / ratio; yr = {c0 - h, c0 + h}; }
        }
    }
    L.x.r = xr; L.y.r = yr; L.y2.r = y2r;
    // colour bar
    L.hasCb = false;
    if (X.st.layers) {
        if (p.cbar.on) { L.bar = p.cbar; L.hasCb = true; }
        else if (!p.heatmaps.empty()) {
            const Heatmap& h = p.heatmaps.front();
            L.bar.on = true; L.bar.vmin = h.vmin; L.bar.vmax = h.vmax; L.bar.lut = h.lut;
            L.bar.label = h.label.compare(0, 2, "##") == 0 ? std::string() : h.label.substr(0, h.label.find("##"));
            L.hasCb = true;
        }
    }
    if (L.hasCb) {
        L.cb.r = {std::min(L.bar.vmin, L.bar.vmax), std::max(L.bar.vmin, L.bar.vmax)};
        if (!(L.cb.r.hi > L.cb.r.lo)) L.cb.r.hi = L.cb.r.lo + 1;
        makeAxis(L.cb, std::max(2, std::min(7, int(cellH / (X.tfs * 2.4)))), X.tfs, X.font);
        L.cbW = X.fs * 0.7 + X.fs * 0.9 + X.pad + L.cb.labelW + (L.bar.label.empty() ? 0 : X.lineH + X.fs * 0.2) + (L.cb.exp ? 0 : 0);
    }
    // cellW / cellH: plot frame size (estimate on the first pass). Fewer ticks until the labels do not collide.
    int tx = std::max(2, std::min(9, int(cellW / (X.fs * 4.5)))), ty = std::max(2, std::min(8, int(cellH / (X.tfs * 2.4))));
    for (;;) {
        makeAxis(L.x, tx, X.tfs, X.font);
        const double need = double(L.x.ticks.size()) * (L.x.labelW + X.tfs * 1.2);
        if (need <= cellW || tx <= 2) break;
        --tx;
    }
    for (;;) {
        makeAxis(L.y, ty, X.tfs, X.font);
        if (double(L.y.ticks.size()) * X.tfs * 1.6 <= cellH || ty <= 2) break;
        --ty;
    }
    if (L.hasY2) makeAxis(L.y2, ty, X.tfs, X.font);
}

void drawPanel(Ctx& X, const Panel& p, const PanelLayout& L, double px, double py, double pw, double ph, bool legendInside) {
    Canvas& c = X.c;
    const Style& st = X.st;
    const uint32_t black = rgba(0, 0, 0), gridC = rgba(222, 222, 222);
    // xT / yT: transformed (axis-space) coordinates - ticks; mx / my: data values (log axes transform, <= 0 -> NaN)
    auto xT = [&](double t) { const double f = (t - L.x.r.lo) / (L.x.r.hi - L.x.r.lo); return px + (L.x.inv ? 1 - f : f) * pw; };
    auto yT = [&](const Axis& a, double t) { const double f = (t - a.r.lo) / (a.r.hi - a.r.lo); return py + ph - (a.inv ? 1 - f : f) * ph; };
    auto mx = [&](double v) { return xT(L.x.tf(v)); };
    auto my = [&](const Axis& a, double v) { return yT(a, a.tf(v)); };
    const bool axes = !p.noAxes;
    // grid
    if (st.grid && axes) {
        for (double t : L.x.ticks) { const double xx = xT(t), pts[4] = {xx, py, xx, py + ph}; c.polyline(pts, 2, gridC, 0.35, 0); }
        for (double t : L.y.ticks) { const double yy = yT(L.y, t), pts[4] = {px, yy, px + pw, yy}; c.polyline(pts, 2, gridC, 0.35, 0); }
    }
    // series
    c.pushClip(px, py, pw, ph);
    // C23 native layers below the series: heatmaps (merged runs of equal colour), cell edges, filled polygons
    const bool layers = st.layers;
    if (layers) {
        for (const Heatmap& h : p.heatmaps) {
            if (h.rows <= 0 || h.cols <= 0 || h.v.size() < size_t(h.rows) * size_t(h.cols) || L.x.log || L.y.log) continue;
            const double cw = (h.x1 - h.x0) / h.cols, chh = (h.y1 - h.y0) / h.rows;
            for (int r = 0; r < h.rows; ++r) {
                const double ya = my(L.y, h.y1 - r * chh), yb = my(L.y, h.y1 - (r + 1) * chh);
                int c0 = 0;
                while (c0 < h.cols) {
                    auto val = [&](int cc) { return h.colMajor ? h.v[size_t(cc) * h.rows + r] : h.v[size_t(r) * h.cols + cc]; };
                    const uint32_t col = h.colorAt(val(c0));
                    int c1 = c0 + 1;
                    while (c1 < h.cols && h.colorAt(val(c1)) == col) ++c1;
                    const double xa = mx(h.x0 + c0 * cw), xb = mx(h.x0 + c1 * cw);
                    // a hair of overlap hides the seams between adjacent rectangles in PDF viewers
                    const double ov = 0.15;
                    c.fillRect(std::min(xa, xb) - ov, std::min(ya, yb) - ov, std::fabs(xb - xa) + 2 * ov, std::fabs(yb - ya) + 2 * ov, col);
                    c0 = c1;
                }
            }
            if (st.cellEdges || p.cellEdges) {
                const uint32_t ec = rgba(40, 40, 40, 110);
                for (int k = 0; k <= h.cols; ++k) { const double xx = mx(h.x0 + k * cw), a[4] = {xx, my(L.y, h.y0), xx, my(L.y, h.y1)}; c.polyline(a, 2, ec, 0.25, 0); }
                for (int k = 0; k <= h.rows; ++k) { const double yy = my(L.y, h.y1 - k * chh), a[4] = {mx(h.x0), yy, mx(h.x1), yy}; c.polyline(a, 2, ec, 0.25, 0); }
            }
        }
        std::vector<double> poly, tri;
        for (const Shape& sh : p.shapes) {
            if (sh.kind != Shape::Polygon || sh.x.size() < 3 || sh.y.size() < sh.x.size()) continue;
            poly.clear();
            for (size_t i = 0; i < sh.x.size(); ++i) { poly.push_back(mx(sh.x[i])); poly.push_back(my(L.y, sh.y[i])); }
            if (sh.fill >> 24) {
                tri.clear();
                for (size_t k = 1; k + 1 < sh.x.size(); ++k) {
                    double t[6] = {poly[0], poly[1], poly[2 * k], poly[2 * k + 1], poly[2 * k + 2], poly[2 * k + 3]};
                    if ((t[2] - t[0]) * (t[5] - t[1]) - (t[3] - t[1]) * (t[4] - t[0]) < 0) { std::swap(t[2], t[4]); std::swap(t[3], t[5]); }
                    tri.insert(tri.end(), t, t + 6);
                }
                if (!tri.empty()) c.fillTriangles(tri.data(), int(tri.size() / 6), sh.fill);
            }
            if (sh.stroke >> 24 && sh.widthPx > 0) {
                poly.push_back(poly[0]); poly.push_back(poly[1]);
                c.polyline(poly.data(), int(poly.size() / 2), ink(sh.stroke, sh.keepColor), std::max(0.25, sh.widthPx * 0.5), sh.dash);
            }
        }
    }
    for (double v : p.vlines) {
        const double tv = L.x.tf(v);
        if (!std::isfinite(tv) || tv < L.x.r.lo || tv > L.x.r.hi) continue;
        const double xx = xT(tv), pts[4] = {xx, py, xx, py + ph};
        c.polyline(pts, 2, rgba(120, 120, 120), 0.5, 1);
    }
    std::vector<double> buf;
    for (const auto& s : p.series) {
        const Axis& ya = (s.axis == 1 && L.hasY2) ? L.y2 : L.y;
        const double lw = st.lineWidthPt * s.width;
        const size_t n = std::min(s.x.size(), s.y.size());
        const uint32_t scol = ink(s.color);
        if (layers && !s.ylo.empty()) {   // C23: filled band between ylo and y
            const size_t m = std::min(n, s.ylo.size());
            std::vector<double> tri;
            for (size_t i = 0; i + 1 < m; ++i) {
                const double xa = mx(s.x[i]), xb = mx(s.x[i + 1]);
                const double a0 = my(ya, s.ylo[i]), a1 = my(ya, s.y[i]), b0 = my(ya, s.ylo[i + 1]), b1 = my(ya, s.y[i + 1]);
                if (!std::isfinite(xa + xb + a0 + a1 + b0 + b1)) continue;
                double t1[6] = {xa, a0, xb, b0, xb, b1}, t2[6] = {xa, a0, xb, b1, xa, a1};
                for (double* t : {t1, t2}) {
                    if ((t[2] - t[0]) * (t[5] - t[1]) - (t[3] - t[1]) * (t[4] - t[0]) < 0) { std::swap(t[2], t[4]); std::swap(t[3], t[5]); }
                    tri.insert(tri.end(), t, t + 6);
                }
            }
            if (!tri.empty()) c.fillTriangles(tri.data(), int(tri.size() / 6), s.color);
            continue;
        }
        if (layers && s.barWidth > 0) {   // C23: bars from 0
            for (size_t i = 0; i < n; ++i) {
                const double xa = mx(s.x[i] - 0.5 * s.barWidth), xb = mx(s.x[i] + 0.5 * s.barWidth), y0 = my(ya, 0.0), y1 = my(ya, s.y[i]);
                if (!std::isfinite(xa + xb + y0 + y1)) continue;
                c.fillRect(std::min(xa, xb), std::min(y0, y1), std::fabs(xb - xa), std::fabs(y1 - y0), s.color);
                if (s.barLine >> 24) c.strokeRect(std::min(xa, xb), std::min(y0, y1), std::fabs(xb - xa), std::fabs(y1 - y0), ink(s.barLine), 0.4);
            }
            continue;
        }
        if (!s.scatter) {
            buf.clear();
            auto flush = [&]() { if (buf.size() >= 4) c.polyline(buf.data(), int(buf.size() / 2), scol, lw, s.dash); buf.clear(); };
            for (size_t i = 0; i < n; ++i) {
                const double xx = mx(s.x[i]), yy = my(ya, s.y[i]);
                if (!std::isfinite(xx) || !std::isfinite(yy)) { flush(); continue; }
                // drop points that coincide on the output (keeps PDFs of long runs small)
                if (buf.size() >= 2 && std::fabs(buf[buf.size() - 2] - xx) < 0.05 && std::fabs(buf.back() - yy) < 0.05) continue;
                buf.push_back(xx); buf.push_back(yy);
            }
            flush();
        }
        if (s.scatter || s.markers) {
            const double r = std::max(1.2, lw * 1.25);
            for (size_t i = 0; i < n; ++i)
                if (std::isfinite(mx(s.x[i])) && std::isfinite(my(ya, s.y[i]))) c.circle(mx(s.x[i]), my(ya, s.y[i]), r, scol, scol, 0.3);
        }
    }
    // C23 native layers above the series: lines, markers, text
    if (layers) {
        std::vector<double> pl;
        for (const Shape& sh : p.shapes) {
            if (sh.kind == Shape::Polyline && sh.x.size() >= 2) {
                pl.clear();
                for (size_t i = 0; i < sh.x.size() && i < sh.y.size(); ++i) { pl.push_back(mx(sh.x[i])); pl.push_back(my(L.y, sh.y[i])); }
                c.polyline(pl.data(), int(pl.size() / 2), ink(sh.stroke, sh.keepColor), std::max(0.25, sh.widthPx * 0.5), sh.dash);
            } else if (sh.kind == Shape::Marker && !sh.x.empty() && !sh.y.empty()) {
                const double xx = mx(sh.x[0]), yy = my(L.y, sh.y[0]);
                if (std::isfinite(xx) && std::isfinite(yy))
                    drawMarker(c, xx, yy, std::max(1.0, sh.sizePx * 0.55), sh.marker, sh.fill, sh.stroke ? ink(sh.stroke, sh.keepColor) : 0, std::max(0.3, sh.widthPx * 0.4));
            }
        }
    }
    c.popClip();
    if (layers) {   // text annotations may extend beyond the frame a little (labels of wells at the edge)
        for (const Shape& sh : p.shapes) {
            if (sh.kind != Shape::Text || sh.text.empty() || sh.x.empty() || sh.y.empty()) continue;
            const double xx = mx(sh.x[0]) + sh.offX * 0.75, yy = my(L.y, sh.y[0]) + sh.offY * 0.75;
            if (!std::isfinite(xx) || !std::isfinite(yy) || xx < px - X.fs || xx > px + pw + X.fs || yy < py - X.fs || yy > py + ph + X.fs) continue;
            c.text(xx, sh.vertical ? yy : yy + X.tfs * 0.34, sh.text, X.tfs, sh.vertical ? 1 : sh.halign, sh.vertical, ink(sh.stroke, sh.keepColor));
        }
    }
    if (axes) {
    // frame + ticks
    const double fw = 0.6;
    if (st.box) c.strokeRect(px, py, pw, ph, black, fw);
    else {
        const double a[6] = {px, py, px, py + ph, px + pw, py + ph};
        c.polyline(a, 3, black, fw, 0);
        if (L.hasY2) { const double b[4] = {px + pw, py, px + pw, py + ph}; c.polyline(b, 2, black, fw, 0); }
    }
    const double tl = X.tick, mtl = X.tick * 0.55;
    auto minorTicks = [&](const Axis& a, double lo, double hi, auto&& draw) {
        if (!st.minorTicks || a.ticks.empty()) return;
        if (a.log) {   // 2..9 x 10^n
            for (int d = int(std::floor(lo)) - 1; d <= int(std::ceil(hi)); ++d)
                for (int m = 2; m <= 9; ++m) { const double t = d + std::log10(double(m)); if (t >= lo - 1e-12 && t <= hi + 1e-12) draw(t); }
            return;
        }
        const double m = a.step / std::pow(10.0, std::floor(std::log10(a.step) + 1e-9));
        const int nm = std::fabs(m - 2) < 1e-6 ? 4 : 5;
        const double ms = a.step / nm;
        for (double v = a.ticks.front() - a.step; v <= hi + 1e-12; v += ms) if (v >= lo - 1e-12) draw(v);
    };
    for (size_t i = 0; i < L.x.ticks.size(); ++i) {
        const double xx = xT(L.x.ticks[i]);
        const double a[4] = {xx, py + ph, xx, py + ph - tl};
        c.polyline(a, 2, black, fw, 0);
        if (st.box) { const double b[4] = {xx, py, xx, py + tl}; c.polyline(b, 2, black, fw, 0); }
        c.text(xx, py + ph + X.pad + X.tfs * 0.78, L.x.labels[i], X.tfs, 1, false, black);
    }
    minorTicks(L.x, L.x.r.lo, L.x.r.hi, [&](double v) {
        const double xx = xT(v), a[4] = {xx, py + ph, xx, py + ph - mtl};
        c.polyline(a, 2, black, fw * 0.8, 0);
        if (st.box) { const double b[4] = {xx, py, xx, py + mtl}; c.polyline(b, 2, black, fw * 0.8, 0); }
    });
    for (size_t i = 0; i < L.y.ticks.size(); ++i) {
        const double yy = yT(L.y, L.y.ticks[i]);
        const double a[4] = {px, yy, px + tl, yy};
        c.polyline(a, 2, black, fw, 0);
        if (st.box && !L.hasY2) { const double b[4] = {px + pw, yy, px + pw - tl, yy}; c.polyline(b, 2, black, fw, 0); }
        c.text(px - X.pad, yy + X.tfs * 0.34, L.y.labels[i], X.tfs, 2, false, black);
    }
    minorTicks(L.y, L.y.r.lo, L.y.r.hi, [&](double v) {
        const double yy = yT(L.y, v), a[4] = {px, yy, px + mtl, yy};
        c.polyline(a, 2, black, fw * 0.8, 0);
        if (st.box && !L.hasY2) { const double b[4] = {px + pw, yy, px + pw - mtl, yy}; c.polyline(b, 2, black, fw * 0.8, 0); }
    });
    if (L.hasY2) {
        for (size_t i = 0; i < L.y2.ticks.size(); ++i) {
            const double yy = yT(L.y2, L.y2.ticks[i]);
            const double a[4] = {px + pw, yy, px + pw - tl, yy};
            c.polyline(a, 2, black, fw, 0);
            c.text(px + pw + X.pad, yy + X.tfs * 0.34, L.y2.labels[i], X.tfs, 0, false, black);
        }
    }
    }
    // C23 colour bar right of the frame (after the right axis labels)
    if (L.hasCb && !L.bar.lut.empty()) {
        const double bx = px + pw + (L.hasY2 ? L.y2.labelW + X.pad * 2 + (p.y2label.empty() ? 0 : X.lineH + X.fs * 0.1) : 0) + X.fs * 0.7;
        const double bw = X.fs * 0.9;
        const int nseg = std::min<int>(256, int(L.bar.lut.size()));
        for (int k = 0; k < nseg; ++k) {
            const double f0 = double(k) / nseg, f1 = double(k + 1) / nseg;
            const size_t li = std::min(L.bar.lut.size() - 1, size_t((k + 0.5) / nseg * L.bar.lut.size()));
            const double ya = py + ph - f1 * ph, yb = py + ph - f0 * ph;
            c.fillRect(bx, ya - 0.1, bw, (yb - ya) + 0.2, L.bar.lut[L.bar.vmin <= L.bar.vmax ? li : L.bar.lut.size() - 1 - li] | 0xFF000000u);
        }
        c.strokeRect(bx, py, bw, ph, black, 0.5);
        for (size_t i = 0; i < L.cb.ticks.size(); ++i) {
            const double yy = yT(L.cb, L.cb.ticks[i]);
            const double a[4] = {bx + bw, yy, bx + bw - X.tick * 0.8, yy};
            c.polyline(a, 2, black, 0.5, 0);
            c.text(bx + bw + X.pad, yy + X.tfs * 0.34, L.cb.labels[i], X.tfs, 0, false, black);
        }
        if (L.cb.exp) c.text(bx, py - X.pad * 0.8, expLabel(L.cb.exp), X.tfs, 0, false, black);
        if (!L.bar.label.empty()) c.text(bx + bw + X.pad * 2 + L.cb.labelW + X.fs * 0.95, py + ph * 0.5, L.bar.label, X.fs, 1, true, black);
    }
    if (axes) {
    // axis factor labels (x10^n) at the top of the y axes / end of the x axis
    if (L.y.exp) c.text(px, py - X.pad * 0.8, expLabel(L.y.exp), X.tfs, 0, false, black);
    if (L.hasY2 && L.y2.exp) c.text(px + pw, py - X.pad * 0.8, expLabel(L.y2.exp), X.tfs, 2, false, black);
    if (L.x.exp) c.text(px + pw, py + ph + X.pad + X.tfs * 0.78 + X.lineH, expLabel(L.x.exp), X.tfs, 2, false, black);
    // labels
    // C17: axis labels longer than the frame side shrink (down to 65 %) and are then shortened with an ellipsis - never clipped
    auto fitLabel = [&](const std::string& lab, double room, double& size) {
        size = X.fs;
        const double w = textWidthPt(lab, X.fs, X.font);
        if (w <= room || w <= 0) return lab;
        size = std::max(0.65 * X.fs, X.fs * room / w);
        std::string t = lab;
        while (t.size() > 2 && textWidthPt(t + "...", size, X.font) > room) {
            t.pop_back();
            while (!t.empty() && (static_cast<unsigned char>(t.back()) & 0xC0) == 0x80) t.pop_back();   // whole UTF-8 characters
            if (!t.empty() && static_cast<unsigned char>(t.back()) >= 0xC0) t.pop_back();
        }
        const size_t sp = t.rfind("^{");   // never leave an open superscript
        if (sp != std::string::npos && t.find('}', sp) == std::string::npos) t.resize(sp);
        return t.size() == lab.size() ? t : t + "...";
    };
    double lfs = X.fs;
    if (!p.xlabel.empty()) { const std::string t = fitLabel(p.xlabel, pw, lfs); c.text(px + pw * 0.5, py + ph + X.pad + X.tfs * 0.78 + X.lineH * 1.05 + X.fs * 0.2, t, lfs, 1, false, black); }
    if (!p.ylabel.empty()) { const std::string t = fitLabel(p.ylabel, ph, lfs); c.text(px - X.pad * 2 - L.y.labelW - X.fs * 0.25, py + ph * 0.5, t, lfs, 1, true, black); }
    if (L.hasY2 && !p.y2label.empty()) { const std::string t = fitLabel(p.y2label, ph, lfs); c.text(px + pw + X.pad * 2 + L.y2.labelW + X.fs * 0.95, py + ph * 0.5, t, lfs, 1, true, black); }
    }
    if (st.titles && !p.title.empty()) c.text(px + pw * 0.5, py - X.pad - (L.y.exp || (L.hasY2 && L.y2.exp) ? X.lineH : 0), p.title, X.fs, 1, false, black);
    // legend inside
    if (legendInside) {
        const auto items = legendItems(p);
        if (!items.empty()) {
            double w, h;
            drawLegendBox(X, items, 0, 0, pw * 0.9, true, &w, &h, true);
            Legend where = st.legend;
            if (where == Legend::Auto || where == Legend::OutsideRight || where == Legend::Below) {
                // corner with the fewest data points underneath
                int best = 0; long bestN = -1;
                for (int k = 0; k < 4; ++k) {
                    const double bx = (k == 0 || k == 2) ? px + pw - w : px, by = (k < 2) ? py : py + ph - h;
                    long cnt = 0;
                    for (const auto& s : p.series) {
                        const Axis& ya = (s.axis == 1 && L.hasY2) ? L.y2 : L.y;
                        const size_t n = std::min(s.x.size(), s.y.size());
                        for (size_t i = 0; i < n; ++i) {
                            const double xx = mx(s.x[i]), yy = my(ya, s.y[i]);
                            if (xx >= bx - 2 && xx <= bx + w + 2 && yy >= by - 2 && yy <= by + h + 2) ++cnt;
                        }
                    }
                    if (bestN < 0 || cnt < bestN) { bestN = cnt; best = k; }
                }
                where = best == 0 ? Legend::InsideTopRight : best == 1 ? Legend::InsideTopLeft : best == 2 ? Legend::InsideBottomRight : Legend::InsideBottomLeft;
            }
            const double m = X.fs * 0.4;
            double lx = px + pw - w - m, ly = py + m;
            if (where == Legend::InsideTopLeft) lx = px + m;
            if (where == Legend::InsideBottomRight) ly = py + ph - h - m;
            if (where == Legend::InsideBottomLeft) { lx = px + m; ly = py + ph - h - m; }
            drawLegendBox(X, items, lx, ly, pw * 0.9, true);
        }
    }
}

} // namespace

Figure recolored(const Figure& f, int palette) {
    // same original colour -> same new colour (keeps the grouping of the app plots: per case / per vector)
    Figure g = f;
    std::vector<uint32_t> seen;
    for (auto& p : g.panels)
        for (auto& s : p.series) {
            const uint32_t key = s.color | 0xFF000000u;
            size_t k = 0;
            while (k < seen.size() && seen[k] != key) ++k;
            if (k == seen.size()) seen.push_back(key);
            s.color = (distinctColor(int(k), palette) & 0x00FFFFFFu) | (s.color & 0xFF000000u);   // bands keep their transparency
        }
    return g;
}

void renderSnapshot(const Snapshot& sn, const Style& st, Canvas& c) {
    const double W = st.widthMM / 25.4 * 72.0, H = st.heightMM / 25.4 * 72.0;
    c.font = st.font;
    c.begin(W, H, st.transparent);
    const double sw = sn.x1 - sn.x0, sh = sn.y1 - sn.y0;
    if (!(sw > 0 && sh > 0)) { c.end(); return; }
    const double k = std::min(W / sw, H / sh), ox = (W - sw * k) * 0.5, oy = (H - sh * k) * 0.5;
    auto col = [&](uint32_t u) { return st.printColors ? printColor(u, &sn) : u; };
    if (!st.transparent) c.fillRect(ox, oy, sw * k, sh * k, col(sn.bg) | 0xFF000000u);
    c.pushClip(ox, oy, sw * k, sh * k);
    int cur = -2;
    bool open = false;
    auto setClip = [&](int ci) {
        if (ci == cur) return;
        if (open) c.popClip();
        open = false;
        cur = ci;
        if (ci >= 0 && size_t(4 * ci + 3) < sn.clips.size()) {
            const double* q = &sn.clips[size_t(4 * ci)];
            c.pushClip(ox + (q[0] - sn.x0) * k, oy + (q[1] - sn.y0) * k, (q[2] - q[0]) * k, (q[3] - q[1]) * k);
            open = true;
        }
    };
    auto drawText = [&](const SnapText& t) {
        setClip(t.clip);
        c.textFit(ox + (t.x - sn.x0) * k, oy + (t.y - sn.y0) * k, t.s, t.size * k, t.width * k, t.vertical, col(t.color));
    };
    std::vector<double> buf;
    size_t ti = 0;
    for (size_t r = 0; r < sn.runs.size(); ++r) {
        while (ti < sn.texts.size() && sn.texts[ti].order <= int(r)) drawText(sn.texts[ti++]);
        const Snapshot::Run& R = sn.runs[r];
        if (R.ntri <= 0 || size_t(6 * (R.tri0 + R.ntri)) > sn.xy.size()) continue;
        setClip(R.clip);
        buf.resize(size_t(6 * R.ntri));
        const double* q = &sn.xy[size_t(6 * R.tri0)];
        for (int i = 0; i < 3 * R.ntri; ++i) {
            buf[size_t(2 * i)] = ox + (q[2 * i] - sn.x0) * k;
            buf[size_t(2 * i + 1)] = oy + (q[2 * i + 1] - sn.y0) * k;
        }
        c.fillTriangles(buf.data(), R.ntri, col(R.color));
    }
    while (ti < sn.texts.size()) drawText(sn.texts[ti++]);
    setClip(-2);
    c.popClip();
    c.end();
}

void render(const Figure& fIn, const Style& st, Canvas& c) {
    if (st.exact && !fIn.snap.empty()) { renderSnapshot(fIn.snap, st, c); return; }
    if (st.palette >= 0) {
        Style s2 = st;
        s2.palette = -1;
        render(recolored(fIn, st.palette), s2, c);
        return;
    }
    const Figure& f = fIn;
    const double W = st.widthMM / 25.4 * 72.0, H = st.heightMM / 25.4 * 72.0;
    c.font = st.font;
    c.begin(W, H, st.transparent);
    if (f.panels.empty()) { c.end(); return; }
    Ctx X{st, c, st.fontPt, st.fontPt * 0.9, st.fontPt * 1.15, st.fontPt * 0.38, st.fontPt * 0.35, st.font};
    const int n = int(f.panels.size());
    int cols = std::max(1, std::min(f.cols, n)), rows = std::max(1, (n + cols - 1) / cols);
    if (f.rows * f.cols >= n && f.rows > 0) rows = std::max(1, std::min(f.rows, (n + cols - 1) / cols));
    // Auto legend: if the per-panel legends would cover too much of small panels, collect them below the figure
    Style stLocal = st;
    if (st.legend == Legend::Auto) {
        const double cw = (W - X.fs) / cols, ch = (H - X.fs) / rows;
        const double pwE = cw * 0.72, phE = ch * 0.62;     // rough plot frame size
        for (const auto& p : f.panels) {
            const auto items = legendItems(p);
            if (items.empty()) continue;
            double w = 0, h = 0;
            drawLegendBox(X, items, 0, 0, pwE * 0.9, true, &w, &h, true);
            if (w * h > 0.2 * pwE * phE || w > 0.55 * pwE || h > 0.55 * phE) { stLocal.legend = Legend::Below; break; }
        }
    }
    const Style& stL = stLocal;
    // shared legend: same entries in every panel and more than one panel
    bool shared = false;
    std::vector<LegendItem> sharedItems;
    if (stL.legend != Legend::None) {
        sharedItems = legendItems(f.panels[0]);
        shared = n > 1 && !sharedItems.empty();
        for (int i = 1; i < n && shared; ++i) {
            const auto o = legendItems(f.panels[size_t(i)]);
            if (o.size() != sharedItems.size()) shared = false;
            for (size_t k = 0; shared && k < o.size(); ++k) if (o[k].label != sharedItems[k].label) shared = false;
        }
        if (stL.legend == Legend::OutsideRight || stL.legend == Legend::Below) {
            // outside legends always collect everything
            if (!shared) {
                std::set<std::string> seen;
                sharedItems.clear();
                for (const auto& p : f.panels) for (const auto& it : legendItems(p)) if (seen.insert(it.label).second) sharedItems.push_back(it);
            }
            shared = !sharedItems.empty();
        }
    }
    const double margin = X.fs * 0.5;
    double legW = 0, legH = 0;
    const bool legRight = shared && stL.legend == Legend::OutsideRight;
    int legCols = 1;
    if (shared && !legRight) {
        // fit as many columns as possible below the panels
        double itemW = 0;
        for (const auto& it : sharedItems) itemW = std::max(itemW, X.tfs * 2.2 + X.tfs * 0.45 + textWidthPt(it.label, X.tfs, X.font));
        legCols = std::max(1, std::min(int(sharedItems.size()), int((W - 2 * margin) / (itemW + X.tfs * 1.2))));
        drawLegendBox(X, sharedItems, 0, 0, W - 2 * margin, false, &legW, &legH, true, legCols);
    } else if (legRight) {
        drawLegendBox(X, sharedItems, 0, 0, W * 0.4, true, &legW, &legH, true, 1);
    }
    const double areaW = W - 2 * margin - (legRight ? legW + X.fs * 0.6 : 0);
    const double areaH = H - 2 * margin - (shared && !legRight ? legH + X.fs * 0.3 : 0);
    const double cellW = areaW / cols, cellH = areaH / rows;
    // axes and common margins (aligned plot frames)
    std::vector<PanelLayout> L(static_cast<size_t>(n));
    double ml = 0, mr = 0, mt = 0, mb = 0;
    // pass 1: margins from a frame estimate; pass 2: ticks for the real frame size, final margins
    for (int pass = 0; pass < 2; ++pass) {
        const double fw = pass == 0 ? cellW * 0.75 : std::max(10.0, cellW - ml - mr);
        const double fh = pass == 0 ? cellH * 0.65 : std::max(10.0, cellH - mt - mb);
        ml = mr = mt = mb = 0;
        for (int i = 0; i < n; ++i) {
            const Panel& p = f.panels[size_t(i)];
            computeAxes(X, p, L[size_t(i)], fw, fh);
            const PanelLayout& q = L[size_t(i)];
            double l = q.y.labelW + X.pad * 2 + (p.ylabel.empty() ? 0 : X.lineH + X.fs * 0.1);
            double r = q.hasY2 ? q.y2.labelW + X.pad * 2 + (p.y2label.empty() ? 0 : X.lineH + X.fs * 0.1) : X.fs * 0.9;
            if (q.hasCb) r = (q.hasY2 ? r : 0) + q.cbW + X.fs * 0.3;
            double t = (st.titles && !p.title.empty() ? X.lineH * 1.1 : 0) + (q.y.exp || (q.hasY2 && q.y2.exp) ? X.lineH : X.fs * 0.5);
            double b = X.pad + X.tfs * 1.05 + (p.xlabel.empty() && !q.x.exp ? 0 : X.lineH * 1.05) + X.fs * 0.15;
            if (p.noAxes) { l = X.fs * 0.4; b = X.fs * 0.4; if (!q.hasCb) r = X.fs * 0.4; }
            ml = std::max(ml, l); mr = std::max(mr, r); mt = std::max(mt, t); mb = std::max(mb, b);
        }
    }
    for (int i = 0; i < n; ++i) {
        const int row = i / cols, col = i % cols;
        const double cx = margin + col * cellW, cy = margin + row * cellH;
        const double px = cx + ml, py = cy + mt, pw = std::max(10.0, cellW - ml - mr), ph = std::max(10.0, cellH - mt - mb);
        drawPanel(X, f.panels[size_t(i)], L[size_t(i)], px, py, pw, ph, !shared && stL.legend != Legend::None);
    }
    if (shared) {
        if (legRight) drawLegendBox(X, sharedItems, W - margin - legW, margin + mt, legW, true, nullptr, nullptr, false, 1);
        else drawLegendBox(X, sharedItems, (W - legW) * 0.5, H - margin - legH, legW, false, nullptr, nullptr, false, legCols);
    }
    c.end();
}

// ===================================================================== SVG
namespace {
std::string hex(uint32_t c) {
    char b[8];
    std::snprintf(b, sizeof(b), "#%02x%02x%02x", c & 255, (c >> 8) & 255, (c >> 16) & 255);
    return b;
}
double alphaOf(uint32_t c) { return ((c >> 24) & 255) / 255.0; }
std::string xmlEsc(const std::string& s) {
    std::string o;
    for (char ch : s) {
        switch (ch) {
        case '&': o += "&amp;"; break;
        case '<': o += "&lt;"; break;
        case '>': o += "&gt;"; break;
        case '"': o += "&quot;"; break;
        default: o += ch;
        }
    }
    return o;
}
std::string num(double v) {
    char b[32];
    std::snprintf(b, sizeof(b), "%.2f", v);
    std::string s = b;
    while (s.size() > 1 && s.find('.') != std::string::npos && (s.back() == '0' || s.back() == '.')) { const char ch = s.back(); s.pop_back(); if (ch == '.') break; }
    if (s == "-0") s = "0";
    return s;
}

class SvgCanvas : public Canvas {
public:
    std::ostringstream o;
    int clipId = 0, depth = 0;
    void begin(double w, double h, bool transparent) override {
        o << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
          << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"" << num(w) << "pt\" height=\"" << num(h) << "pt\" viewBox=\"0 0 "
          << num(w) << " " << num(h) << "\">\n<!-- ResSim Studio figure -->\n"
          << "<g font-family=\"" << (font == 1 ? "'Times New Roman', Times, serif" : "Helvetica, Arial, sans-serif") << "\">\n";
        if (!transparent) o << "<rect x=\"0\" y=\"0\" width=\"" << num(w) << "\" height=\"" << num(h) << "\" fill=\"#ffffff\"/>\n";
    }
    void end() override { while (depth-- > 0) o << "</g>\n"; o << "</g>\n</svg>\n"; }
    void stroke(uint32_t color, double w, int dash) {
        o << " fill=\"none\" stroke=\"" << hex(color) << "\" stroke-width=\"" << num(w) << "\" stroke-linecap=\"round\" stroke-linejoin=\"round\"";
        if (alphaOf(color) < 0.999) o << " stroke-opacity=\"" << num(alphaOf(color)) << "\"";
        double gap = 0;
        const double d = dashLen(dash, w, &gap);
        if (d > 0) {
            if (dash == 3) o << " stroke-dasharray=\"" << num(d) << "," << num(gap) << "," << num(w * 0.8) << "," << num(gap) << "\"";
            else o << " stroke-dasharray=\"" << num(d) << "," << num(gap) << "\"";
        }
    }
    void polyline(const double* xy, int n, uint32_t color, double w, int dash) override {
        if (n < 2) return;
        o << "<polyline points=\"";
        for (int i = 0; i < n; ++i) o << (i ? " " : "") << num(xy[2 * i]) << "," << num(xy[2 * i + 1]);
        o << "\"";
        stroke(color, w, dash);
        o << "/>\n";
    }
    void fillRect(double x, double y, double w, double h, uint32_t color) override {
        o << "<rect x=\"" << num(x) << "\" y=\"" << num(y) << "\" width=\"" << num(w) << "\" height=\"" << num(h) << "\" fill=\"" << hex(color) << "\"";
        if (alphaOf(color) < 0.999) o << " fill-opacity=\"" << num(alphaOf(color)) << "\"";
        o << "/>\n";
    }
    void strokeRect(double x, double y, double w, double h, uint32_t color, double sw) override {
        o << "<rect x=\"" << num(x) << "\" y=\"" << num(y) << "\" width=\"" << num(w) << "\" height=\"" << num(h) << "\"";
        stroke(color, sw, 0);
        o << "/>\n";
    }
    void circle(double x, double y, double r, uint32_t fill, uint32_t strokeC, double sw) override {
        o << "<circle cx=\"" << num(x) << "\" cy=\"" << num(y) << "\" r=\"" << num(r) << "\" fill=\"" << hex(fill) << "\" stroke=\"" << hex(strokeC)
          << "\" stroke-width=\"" << num(sw) << "\"/>\n";
    }
    void text(double x, double y, const std::string& s, double size, int ha, bool vertical, uint32_t color) override {
        const double w = textWidth(s, size);
        double dx = ha == 1 ? -w * 0.5 : ha == 2 ? -w : 0;
        o << "<text";
        if (vertical) o << " transform=\"translate(" << num(x) << "," << num(y) << ") rotate(-90)\" x=\"" << num(dx) << "\" y=\"0\"";
        else o << " x=\"" << num(x + dx) << "\" y=\"" << num(y) << "\"";
        o << " font-size=\"" << num(size) << "\" fill=\"" << hex(color) << "\">";
        for (const auto& g : segments(s)) {
            if (g.sup) o << "<tspan font-size=\"" << num(size * kSup) << "\" dy=\"" << num(-size * kRise) << "\">" << xmlEsc(g.text) << "</tspan><tspan dy=\"" << num(size * kRise) << "\">&#8203;</tspan>";
            else o << xmlEsc(g.text);
        }
        o << "</text>\n";
    }
    void fillTriangles(const double* xy, int nTri, uint32_t color) override {
        if (nTri <= 0 || alphaOf(color) <= 0) return;
        o << "<path d=\"";
        for (int t = 0; t < nTri; ++t) {
            const double* q = xy + 6 * t;
            o << "M" << num(q[0]) << " " << num(q[1]) << "L" << num(q[2]) << " " << num(q[3]) << "L" << num(q[4]) << " " << num(q[5]) << "Z";
        }
        o << "\" fill=\"" << hex(color) << "\"";
        if (alphaOf(color) < 0.999) o << " fill-opacity=\"" << num(alphaOf(color)) << "\"";
        o << "/>\n";
    }
    void textFit(double x, double y, const std::string& s, double size, double widthPt, bool vertical, uint32_t color) override {
        if (s.empty()) return;
        o << "<text";
        if (vertical) o << " transform=\"translate(" << num(x) << "," << num(y) << ") rotate(-90)\" x=\"0\" y=\"0\"";
        else o << " x=\"" << num(x) << "\" y=\"" << num(y) << "\"";
        o << " font-size=\"" << num(size) << "\" fill=\"" << hex(color) << "\"";
        if (alphaOf(color) < 0.999) o << " fill-opacity=\"" << num(alphaOf(color)) << "\"";
        if (widthPt > 0) o << " textLength=\"" << num(widthPt) << "\" lengthAdjust=\"spacingAndGlyphs\"";
        o << " xml:space=\"preserve\">" << xmlEsc(s) << "</text>\n";
    }
    void pushClip(double x, double y, double w, double h) override {
        ++clipId;
        o << "<clipPath id=\"c" << clipId << "\"><rect x=\"" << num(x) << "\" y=\"" << num(y) << "\" width=\"" << num(w) << "\" height=\"" << num(h)
          << "\"/></clipPath>\n<g clip-path=\"url(#c" << clipId << ")\">\n";
        ++depth;
    }
    void popClip() override { if (depth > 0) { o << "</g>\n"; --depth; } }
};

// ===================================================================== PDF
std::string winAnsi(const std::string& s) {
    std::string o;
    for (unsigned u : decode(s)) {
        if (u < 128) o += char(u);
        else if (u == 0x2212 || u == 0x2013) o += char(0x96);   // minus -> en dash
        else if (u == 0x2014) o += char(0x97);
        else if (u == 0x2022) o += char(0x95);
        else if (u == 0x2026) o += char(0x85);
        else if (u == 0x200B) continue;
        else if (u < 256) o += char(u);
        else o += '?';
    }
    return o;
}
std::string pdfStr(const std::string& s) {
    std::string o = "(";
    for (char ch : winAnsi(s)) {
        if (ch == '(' || ch == ')' || ch == '\\') { o += '\\'; o += ch; }
        else o += ch;
    }
    return o + ")";
}

class PdfCanvas : public Canvas {
public:
    std::ostringstream o;
    double W = 0, H = 0;
    void begin(double w, double h, bool transparent) override {
        W = w; H = h;
        o << "1 J 1 j\n";
        if (!transparent) o << "1 1 1 rg 0 0 " << num(w) << " " << num(h) << " re f\n";
    }
    void end() override {}
    void rgb(uint32_t c, bool strokeOp) {
        o << num((c & 255) / 255.0) << " " << num(((c >> 8) & 255) / 255.0) << " " << num(((c >> 16) & 255) / 255.0) << (strokeOp ? " RG\n" : " rg\n");
    }
    void style(uint32_t color, double w, int dash) {
        rgb(color, true);
        o << num(w) << " w\n";
        double gap = 0;
        const double d = dashLen(dash, w, &gap);
        if (d > 0) {
            if (dash == 3) o << "[" << num(d) << " " << num(gap) << " " << num(w * 0.8) << " " << num(gap) << "] 0 d\n";
            else o << "[" << num(d) << " " << num(gap) << "] 0 d\n";
        } else o << "[] 0 d\n";
    }
    void polyline(const double* xy, int n, uint32_t color, double w, int dash) override {
        if (n < 2 || alphaOf(color) <= 0) return;
        style(color, w, dash);
        for (int i = 0; i < n; ++i) o << num(xy[2 * i]) << " " << num(H - xy[2 * i + 1]) << (i ? " l\n" : " m\n");
        o << "S\n";
    }
    void fillRect(double x, double y, double w, double h, uint32_t color) override {
        if (alphaOf(color) <= 0) return;
        rgb(color, false);
        o << num(x) << " " << num(H - y - h) << " " << num(w) << " " << num(h) << " re f\n";
    }
    void strokeRect(double x, double y, double w, double h, uint32_t color, double sw) override {
        style(color, sw, 0);
        o << num(x) << " " << num(H - y - h) << " " << num(w) << " " << num(h) << " re S\n";
    }
    void circle(double x, double y, double r, uint32_t fill, uint32_t strokeC, double sw) override {
        const double k = 0.5523 * r, Y = H - y;
        rgb(fill, false);
        rgb(strokeC, true);
        o << num(sw) << " w [] 0 d\n";
        o << num(x + r) << " " << num(Y) << " m\n"
          << num(x + r) << " " << num(Y + k) << " " << num(x + k) << " " << num(Y + r) << " " << num(x) << " " << num(Y + r) << " c\n"
          << num(x - k) << " " << num(Y + r) << " " << num(x - r) << " " << num(Y + k) << " " << num(x - r) << " " << num(Y) << " c\n"
          << num(x - r) << " " << num(Y - k) << " " << num(x - k) << " " << num(Y - r) << " " << num(x) << " " << num(Y - r) << " c\n"
          << num(x + k) << " " << num(Y - r) << " " << num(x + r) << " " << num(Y - k) << " " << num(x + r) << " " << num(Y) << " c\nb\n";
    }
    void text(double x, double y, const std::string& s, double size, int ha, bool vertical, uint32_t color) override {
        const double w = textWidth(s, size);
        const double d = ha == 1 ? -w * 0.5 : ha == 2 ? -w : 0;
        rgb(color, false);
        o << "BT\n";
        if (vertical) o << "0 1 -1 0 " << num(x) << " " << num(H - y + d) << " Tm\n";
        else o << "1 0 0 1 " << num(x + d) << " " << num(H - y) << " Tm\n";
        for (const auto& g : segments(s)) {
            if (g.sup) o << "/F1 " << num(size * kSup) << " Tf " << num(size * kRise) << " Ts " << pdfStr(g.text) << " Tj 0 Ts\n";
            else o << "/F1 " << num(size) << " Tf " << pdfStr(g.text) << " Tj\n";
        }
        o << "ET\n";
    }
    std::vector<int> alphas;   // ExtGState /Ann (fill + stroke alpha nn %) used by the content
    bool alphaOn(uint32_t c) {   // opens "q /Ann gs" for a translucent colour; returns true if Q is needed
        const int a = int(std::lround(alphaOf(c) * 100));
        if (a >= 100) return false;
        if (std::find(alphas.begin(), alphas.end(), a) == alphas.end()) alphas.push_back(a);
        o << "q /A" << a << " gs\n";
        return true;
    }
    void fillTriangles(const double* xy, int nTri, uint32_t color) override {
        if (nTri <= 0 || alphaOf(color) <= 0) return;
        const bool q = alphaOn(color);
        rgb(color, false);
        for (int t = 0; t < nTri; ++t) {
            const double* p = xy + 6 * t;
            o << num(p[0]) << " " << num(H - p[1]) << " m " << num(p[2]) << " " << num(H - p[3]) << " l " << num(p[4]) << " " << num(H - p[5]) << " l h\n";
        }
        o << "f\n";
        if (q) o << "Q\n";
    }
    void textFit(double x, double y, const std::string& s, double size, double widthPt, bool vertical, uint32_t color) override {
        if (s.empty() || alphaOf(color) <= 0) return;
        const double w = textWidth(s, size);
        const double tz = (w > 0 && widthPt > 0) ? std::clamp(100.0 * widthPt / w, 40.0, 250.0) : 100.0;
        const bool q = alphaOn(color);
        rgb(color, false);
        o << "BT\n";
        if (vertical) o << "0 1 -1 0 " << num(x) << " " << num(H - y) << " Tm\n";
        else o << "1 0 0 1 " << num(x) << " " << num(H - y) << " Tm\n";
        o << "/F1 " << num(size) << " Tf " << num(tz) << " Tz " << pdfStr(s) << " Tj 100 Tz\nET\n";
        if (q) o << "Q\n";
    }
    void pushClip(double x, double y, double w, double h) override {
        o << "q " << num(x) << " " << num(H - y - h) << " " << num(w) << " " << num(h) << " re W n\n";
    }
    void popClip() override { o << "Q\n"; }
};
} // namespace

std::string svgString(const Figure& f, const Style& s) {
    SvgCanvas c;
    render(f, s, c);
    return c.o.str();
}

bool writeSVG(const Figure& f, const Style& s, const std::string& path, std::string* err) {
    std::ofstream out(path, std::ios::binary);
    if (!out) { if (err) *err = "cannot write " + path; return false; }
    out << svgString(f, s);
    return bool(out);
}

bool writePDF(const Figure& f, const Style& s, const std::string& path, std::string* err) {
    PdfCanvas c;
    render(f, s, c);
    const std::string content = c.o.str();
    std::ostringstream pdf;
    std::vector<long> off;
    pdf << "%PDF-1.4\n%\xE2\xE3\xCF\xD3\n";
    auto obj = [&](const std::string& body) { off.push_back(long(pdf.tellp())); pdf << off.size() << " 0 obj\n" << body << "\nendobj\n"; };
    char date[32];
    const std::time_t t = std::time(nullptr);
    std::strftime(date, sizeof(date), "D:%Y%m%d%H%M%S", std::gmtime(&t));
    obj("<< /Type /Catalog /Pages 2 0 R >>");
    obj("<< /Type /Pages /Kids [3 0 R] /Count 1 >>");
    std::string gs;
    for (int a : c.alphas) gs += " /A" + std::to_string(a) + " << /Type /ExtGState /ca " + num(a / 100.0) + " /CA " + num(a / 100.0) + " >>";
    if (!gs.empty()) gs = " /ExtGState <<" + gs + " >>";
    obj("<< /Type /Page /Parent 2 0 R /MediaBox [0 0 " + num(c.W) + " " + num(c.H) + "] /Resources << /Font << /F1 5 0 R >>" + gs + " >> /Contents 4 0 R >>");
    obj("<< /Length " + std::to_string(content.size()) + " >>\nstream\n" + content + "endstream");
    obj(std::string("<< /Type /Font /Subtype /Type1 /BaseFont /") + (s.font == 1 ? "Times-Roman" : "Helvetica") + " /Encoding /WinAnsiEncoding >>");
    obj(std::string("<< /Producer (ResSim Studio) /Creator (ResSim Studio figure export) /Title ") + pdfStr(f.name.empty() ? "Figure" : f.name) +
        " /CreationDate (" + date + "Z) >>");
    const long xref = long(pdf.tellp());
    pdf << "xref\n0 " << off.size() + 1 << "\n0000000000 65535 f \n";
    for (long v : off) { char b[24]; std::snprintf(b, sizeof(b), "%010ld 00000 n \n", v); pdf << b; }
    pdf << "trailer\n<< /Size " << off.size() + 1 << " /Root 1 0 R /Info " << off.size() << " 0 R >>\nstartxref\n" << xref << "\n%%EOF\n";
    std::ofstream out(path, std::ios::binary);
    if (!out) { if (err) *err = "cannot write " + path; return false; }
    const std::string b = pdf.str();
    out.write(b.data(), std::streamsize(b.size()));
    return bool(out);
}

bool writeCSV(const Figure& f, const std::string& path, std::string* err) {
    std::ofstream out(path);
    if (!out) { if (err) *err = "cannot write " + path; return false; }
    // one block per panel: x,y pairs of every series side by side
    for (const auto& p : f.panels) {
        out << "# " << (p.title.empty() ? p.ylabel : p.title) << "\n";
        size_t rows = 0;
        bool first = true;
        for (const auto& s : p.series) {
            out << (first ? "" : ",") << "\"" << p.xlabel << " (" << s.label << ")\",\"" << s.label << "\"";
            first = false;
            rows = std::max(rows, std::min(s.x.size(), s.y.size()));
        }
        out << "\n";
        for (size_t r = 0; r < rows; ++r) {
            first = true;
            for (const auto& s : p.series) {
                out << (first ? "" : ",");
                first = false;
                if (r < s.x.size() && r < s.y.size()) {
                    char b[64];
                    std::snprintf(b, sizeof(b), "%.10g,%.10g", s.x[r], s.y[r]);
                    out << b;
                } else out << ",";
            }
            out << "\n";
        }
        out << "\n";
    }
    return bool(out);
}

// ===================================================================== recorder
Recorder& recorder() { static Recorder r; return r; }
void Recorder::beginFigure(const std::string& id, int rows, int cols) {
    cur_ = Figure();
    cur_.name = id;
    cur_.rows = rows; cur_.cols = cols;
    active_ = true;
    panelOpen_ = false;
}
void Recorder::beginPanel(const std::string& title, const std::string& xl, const std::string& yl, const std::string& y2l) {
    if (!active_) return;
    panel_ = Panel();
    panel_.title = title.compare(0, 2, "##") == 0 ? std::string() : title;
    panel_.xlabel = xl; panel_.ylabel = yl; panel_.y2label = y2l;
    panelOpen_ = true;
}
void Recorder::series(const std::string& label, const double* x, const double* y, int n, uint32_t color, float width, bool scatter, int axis,
                      bool markers, bool legend) {
    if (!active_ || !panelOpen_) return;
    Series s;
    const size_t hh = label.find("##");
    s.label = hh == std::string::npos ? label : label.substr(0, hh);
    s.x.assign(x, x + std::max(0, n));
    s.y.assign(y, y + std::max(0, n));
    s.color = color | 0xFF000000u;
    s.width = width;
    s.scatter = scatter;
    s.axis = axis;
    s.markers = markers;
    s.legend = legend && !s.label.empty();
    panel_.series.push_back(std::move(s));
}
Panel* Recorder::openPanel() { return active_ && panelOpen_ ? &panel_ : nullptr; }
Panel* Recorder::lastPanel() { return active_ && !panelOpen_ && !cur_.panels.empty() ? &cur_.panels.back() : nullptr; }
uint32_t Heatmap::colorAt(double val) const {
    if (lut.empty()) return 0xFF808080u;
    double t = vmax != vmin ? (val - vmin) / (vmax - vmin) : 0.0;
    if (!std::isfinite(t)) t = 0;
    t = std::min(1.0, std::max(0.0, t));
    return lut[std::min(lut.size() - 1, size_t(t * double(lut.size() - 1) + 0.5))] | 0xFF000000u;
}
Series* Recorder::openSeries() { return active_ && panelOpen_ && !panel_.series.empty() ? &panel_.series.back() : nullptr; }
void Recorder::vline(double x) { if (active_ && panelOpen_) panel_.vlines.push_back(x); }
void Recorder::limits(double x0, double x1, double y0, double y1, bool hasY2, double y20, double y21) {
    if (!active_ || !panelOpen_) return;
    panel_.hasLimits = true;
    panel_.x0 = x0; panel_.x1 = x1; panel_.y0 = y0; panel_.y1 = y1;
    panel_.hasY2 = hasY2; panel_.y20 = y20; panel_.y21 = y21;
}
void Recorder::axisScale(bool logX, bool logY, bool logY2) {
    if (!active_ || !panelOpen_) return;
    panel_.logX = logX; panel_.logY = logY; panel_.logY2 = logY2;
}
void Recorder::endPanel() {
    if (!active_ || !panelOpen_) return;
    cur_.panels.push_back(std::move(panel_));
    panelOpen_ = false;
}
void Recorder::endFigure() {
    if (!active_) return;
    if (panelOpen_) endPanel();
    figs_[cur_.name] = std::move(cur_);
    active_ = false;
}
const Figure* Recorder::get(const std::string& id) const {
    auto it = figs_.find(id);
    return it == figs_.end() ? nullptr : &it->second;
}

} // namespace fig
