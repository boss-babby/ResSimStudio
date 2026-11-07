#include "plotkit.h"
#include "implot_internal.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cctype>
#include <cstring>
#include <string>
#include <vector>

int g_plotPalette = 0;
float g_plotLineScale = 1.0f;

namespace {
int g_axis = 0;
bool g_hasY2 = false;
}

PlotGrid fitGrid(int n, ImVec2 avail, float minW, float minH, float aspect, int forceCols) {
    PlotGrid g;
    n = std::max(1, n);
    float best = -1;
    for (int cols = 1; cols <= n; ++cols) {
        if (forceCols > 0 && cols != forceCols) continue;
        const int rows = (n + cols - 1) / cols;
        const float cw = avail.x / cols, ch = avail.y / rows;
        // effective plot size: limited by the width (at the preferred aspect) and the height
        float score = std::min(cw / aspect, ch);
        if (cw < minW) score *= 0.25f * cw / minW;
        if (ch < minH) score *= 0.5f + 0.5f * ch / minH;
        score *= 1.0f - 0.04f * float(rows * cols - n);   // few empty cells
        if (score > best) { best = score; g.cols = cols; g.rows = rows; }
    }
    g.size = avail;
    if (avail.y / g.rows < minH) { g.size.y = g.rows * minH; g.scroll = true; }
    return g;
}

int formatSI(double v, char* buf, int size, void*) {
    const double a = std::fabs(v);
    if (a >= 1e9) return std::snprintf(buf, size_t(size), "%gG", v / 1e9);
    if (a >= 1e6) return std::snprintf(buf, size_t(size), "%gM", v / 1e6);
    if (a >= 1e4) return std::snprintf(buf, size_t(size), "%gk", v / 1e3);
    return std::snprintf(buf, size_t(size), "%g", v);
}
void setupSI(ImAxis axis) { ImPlot::SetupAxisFormat(axis, formatSI, nullptr); }

void recPanel(const char* title, const char* xl, const char* yl, const char* y2l) {
    g_axis = 0;
    g_hasY2 = false;
    fig::recorder().beginPanel(title ? title : "", xl ? xl : "", yl ? yl : "", y2l ? y2l : "");
    ImPlot::GetCaptureHooks().RecordAll = fig::recorder().openPanel() != nullptr;   // C23: heatmaps / text of this plot are recorded
}
void recEnd() {
    ImPlot::GetCaptureHooks().RecordAll = false;
    if (!fig::recorder().recording()) return;
    if (fig::Panel* P = fig::recorder().openPanel())
        if (const ImPlotPlot* pl = ImPlot::GetCurrentPlot()) {   // C23: equal scaling / inverted axes
            P->equal = ImHasFlag(pl->Flags, ImPlotFlags_Equal);
            P->invertX = pl->Axes[ImAxis_X1].IsInverted();
            P->invertY = pl->Axes[ImAxis_Y1].IsInverted();
        }
    const ImPlotRect r = ImPlot::GetPlotLimits(ImAxis_X1, ImAxis_Y1);
    ImPlotRect r2 = r;
    if (g_hasY2) r2 = ImPlot::GetPlotLimits(ImAxis_X1, ImAxis_Y2);
    fig::recorder().limits(r.X.Min, r.X.Max, r.Y.Min, r.Y.Max, g_hasY2, r2.Y.Min, r2.Y.Max);
    if (const ImPlotPlot* pl = ImPlot::GetCurrentPlot())   // C17: log axes carry over to the exported figure
        fig::recorder().axisScale(pl->Axes[ImAxis_X1].Scale == ImPlotScale_Log10, pl->Axes[ImAxis_Y1].Scale == ImPlotScale_Log10,
                                  pl->Axes[ImAxis_Y2].Scale == ImPlotScale_Log10);
    fig::recorder().endPanel();
}
void recSeries(const char* label, const double* x, const double* y, int n, ImU32 col, float weight, bool scatter, bool markers) {
    if (!fig::recorder().recording()) return;
    fig::recorder().series(label ? label : "", x, y, n, uint32_t(col), weight / 1.6f, scatter, g_axis, markers);
}
void recVLine(double x) { fig::recorder().vline(x); }
void useY2(bool on) {
    ImPlot::SetAxes(ImAxis_X1, on ? ImAxis_Y2 : ImAxis_Y1);
    g_axis = on ? 1 : 0;
    if (on) g_hasY2 = true;
}

// ---------------------------------------------------------------- single-plot capture (C14)
namespace {
struct PText { fig::SnapText t; ImVec4 clip; int idx; };
struct Capture {
    fig::Panel panel;
    std::vector<std::pair<double, uint32_t>> hlines;   // horizontal lines (y, colour), drawn over the x range at the end
    std::string title;
    int unnamed = 0;
    ImGuiID plotId = 0;
    std::function<void(const fig::Figure&)> done;
    // C22 exact snapshots
    fig::Snapshot snap, dsnap;      // single captured plot / requested recorded figure (dashboard)
    std::vector<PText> texts;       // strings of the current geometry region (text hook)
    int lastTarget = 0;             // 1 snap, 2 dsnap: where a colour scale right after the plot goes
    bool pendingFig = false;        // data capture finished - delivered at the next frame (a colour scale may follow)
    fig::Figure fig;
    std::string dashId;
    int dashState = 0;              // 1 requested, 2 capturing this frame
    std::function<void(const std::string&, const fig::Snapshot&)> dashDone;
};
Capture& cap() { static Capture c; return c; }
std::string cleanLabel(const char* s) {
    std::string l = s ? s : "";
    const size_t h = l.find("##");
    if (h != std::string::npos) l.resize(h);
    return l;
}
void hookBegin(ImGuiID id, const char* title, void*) {
    Capture& c = cap();
    c.panel = fig::Panel(); c.hlines.clear(); c.unnamed = 0;
    c.plotId = id;
    c.snap = fig::Snapshot();
    c.title = cleanLabel(title);
    c.panel.title = c.title;
}
void hookSeries(const char* label, const double* xs, const double* ys, int n, ImU32 col, float weight, int kind, int yaxis, void*) {
    Capture& c = cap();
    if ((uint32_t(col) & 0xFF000000u) == 0 || (kind == 1 && weight <= 0)) return;   // invisible helper items (fit ranges)
    if (kind == 5) { for (int i = 0; i < n; ++i) c.panel.vlines.push_back(xs[i]); return; }
    if (kind == 6) { for (int i = 0; i < n; ++i) c.hlines.push_back({ys[i], uint32_t(col)}); return; }
    fig::Series s;
    s.label = cleanLabel(label);
    if (s.label.empty()) { s.legend = false; s.label = "series " + std::to_string(++c.unnamed); }
    s.x.assign(xs, xs + n); s.y.assign(ys, ys + n);
    s.color = uint32_t(col) | 0xFF000000u;
    s.width = std::max(0.4f, weight / 1.6f);
    s.scatter = kind == 1;
    s.markers = kind == 3 || kind == 7;
    s.axis = yaxis > 0 ? 1 : 0;
    if (kind == 4) {   // shaded band: its two bounding curves, one legend entry
        for (const auto& o : c.panel.series) if (o.label == s.label && o.color == s.color) s.legend = false;
        s.dash = 2;
    }
    c.panel.series.push_back(std::move(s));
}
void hookAxis(int axis, const char* label, double mn, double mx, void*) {
    fig::Panel& p = cap().panel;
    const std::string l = cleanLabel(label);
    const ImPlotPlot* pl = ImPlot::GetCurrentPlot();
    const bool lg = pl && axis >= 0 && axis < ImAxis_COUNT && pl->Axes[axis].Scale == ImPlotScale_Log10;
    if (axis == ImAxis_X1) p.logX = lg; else if (axis == ImAxis_Y1) p.logY = lg; else if (axis == ImAxis_Y2) p.logY2 = lg;
    if (pl) {   // C23
        p.equal = ImHasFlag(pl->Flags, ImPlotFlags_Equal);
        p.noAxes = ImHasFlag(pl->Axes[ImAxis_X1].Flags, ImPlotAxisFlags_NoTickLabels) && ImHasFlag(pl->Axes[ImAxis_Y1].Flags, ImPlotAxisFlags_NoTickLabels);
        if (axis == ImAxis_X1) p.invertX = pl->Axes[axis].IsInverted();
        if (axis == ImAxis_Y1) p.invertY = pl->Axes[axis].IsInverted();
    }
    if (axis == ImAxis_X1) { p.xlabel = l; p.x0 = mn; p.x1 = mx; p.hasLimits = true; }
    else if (axis == ImAxis_Y1) { p.ylabel = l; p.y0 = mn; p.y1 = mx; }
    else if (axis == ImAxis_Y2) { p.y2label = l; p.y20 = mn; p.y21 = mx; p.hasY2 = true; }
}
void hookEnd(void*) {
    Capture& c = cap();
    for (const auto& h : c.hlines) {
        fig::Series s; s.label = "y = " + std::to_string(h.first); s.legend = false; s.color = h.second | 0xFF000000u; s.dash = 1;
        s.x = {c.panel.x0, c.panel.x1}; s.y = {h.first, h.first};
        c.panel.series.push_back(std::move(s));
    }
    fig::Figure F;
    F.name = c.title.empty() ? std::string("plot") : c.title;
    F.rows = F.cols = 1;
    F.panels.push_back(c.panel);
    c.fig = std::move(F);
    c.pendingFig = true;   // handed over by plotCaptureNewFrame() together with the exact snapshot
}

// ---- C22 exact snapshot: geometry of the draw list -> fig::Snapshot
std::string stripIcons(const char* b, const char* e) {
    std::string o;
    for (const char* p = b; p < e;) {
        unsigned int c = 0;
        const int n = ImTextCharFromUtf8(&c, p, e);
        if (n <= 0) break;
        const bool icon = (c >= 0xE000 && c <= 0xF8FF) || c >= 0xF0000 || c == 0xFFFD;
        if (!icon) o.append(p, size_t(n));
        p += n;
    }
    return o;
}
void hookText(ImDrawList* dl, ImFont* font, float size, const ImVec2& pos, ImU32 col, const char* b, const char* e, const ImVec4& clip, bool vertical) {
    if (!font || !b) return;
    if (!e) e = b + strlen(b);
    ImFontBaked* fb = font->GetFontBaked(size);
    const float asc = fb && fb->Size > 0 ? fb->Ascent * (size / fb->Size) : size * 0.8f;
    int line = 0;
    for (const char* l = b; l <= e;) {
        const char* nl = (const char*)memchr(l, '\n', size_t(e - l));
        const char* le = nl ? nl : e;
        std::string str = stripIcons(l, le);
        if (str.find_first_not_of(' ') != std::string::npos) {
            PText t;
            t.t.s = str;
            t.t.size = size;
            t.t.width = font->CalcTextSizeA(size, FLT_MAX, 0, str.c_str(), str.c_str() + str.size()).x;
            t.t.color = uint32_t(col);
            t.t.vertical = vertical;
            if (!vertical) { t.t.x = pos.x; t.t.y = pos.y + asc + line * size; }
            else { t.t.x = pos.x + asc + line * size; t.t.y = pos.y; }
            t.clip = clip;
            t.idx = dl->IdxBuffer.Size;
            cap().texts.push_back(std::move(t));
        }
        if (!nl) break;
        l = nl + 1;
        ++line;
    }
}
struct SnapBuilder {
    fig::Snapshot& S;
    bool brk = true;
    int clipIndex(ImVec4 c, const ImVec4& fr) {
        c.x = std::max(c.x, fr.x); c.y = std::max(c.y, fr.y); c.z = std::min(c.z, fr.z); c.w = std::min(c.w, fr.w);
        const size_t n = S.clips.size() / 4;
        for (size_t i = n; i-- > 0 && i + 8 > n;)   // recent clips first
            if (S.clips[4 * i] == c.x && S.clips[4 * i + 1] == c.y && S.clips[4 * i + 2] == c.z && S.clips[4 * i + 3] == c.w) return int(i);
        S.clips.insert(S.clips.end(), {double(c.x), double(c.y), double(c.z), double(c.w)});
        return int(n);
    }
    void tri(double x0, double y0, double x1, double y1, double x2, double y2, uint32_t col, int clip) {
        const double cr = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);
        if (std::fabs(cr) < 1e-7) return;
        if (cr < 0) { std::swap(x1, x2); std::swap(y1, y2); }
        const int nt = int(S.xy.size() / 6);
        S.xy.insert(S.xy.end(), {x0, y0, x1, y1, x2, y2});
        if (!brk && !S.runs.empty()) {
            fig::Snapshot::Run& r = S.runs.back();
            if (r.color == col && r.clip == clip && r.tri0 + r.ntri == nt) { ++r.ntri; return; }
        }
        S.runs.push_back({col, clip, nt, 1});
        brk = false;
    }
    static uint32_t mix(uint32_t a, uint32_t b) {
        uint32_t o = 0;
        for (int k = 0; k < 32; k += 8) o |= (((a >> k & 255) + (b >> k & 255) + 1) / 2) << k;
        return o;
    }
    static int diff(uint32_t a, uint32_t b) {
        int d = 0;
        for (int k = 0; k < 32; k += 8) d = std::max(d, std::abs(int(a >> k & 255) - int(b >> k & 255)));
        return d;
    }
    // per-vertex colours (gradients, e.g. colour scales): subdivide until the colours are close
    void triCol(const ImVec2& a, const ImVec2& b, const ImVec2& c, uint32_t ca, uint32_t cb, uint32_t cc, int clip, int depth) {
        if (depth >= 5 || (diff(ca, cb) <= 6 && diff(cb, cc) <= 6 && diff(ca, cc) <= 6)) {
            tri(a.x, a.y, b.x, b.y, c.x, c.y, mix(mix(ca, cb), cc), clip);
            return;
        }
        const ImVec2 ab((a.x + b.x) * 0.5f, (a.y + b.y) * 0.5f), bc((b.x + c.x) * 0.5f, (b.y + c.y) * 0.5f), ca2((c.x + a.x) * 0.5f, (c.y + a.y) * 0.5f);
        const uint32_t cab = mix(ca, cb), cbc = mix(cb, cc), cca = mix(cc, ca);
        triCol(a, ab, ca2, ca, cab, cca, clip, depth + 1);
        triCol(ab, b, bc, cab, cb, cbc, clip, depth + 1);
        triCol(ca2, bc, c, cca, cbc, cc, clip, depth + 1);
        triCol(ab, bc, ca2, cab, cbc, cca, clip, depth + 1);
    }
};
void addGeometry(fig::Snapshot& S, ImDrawList* dl, int idx0, int idx1, const ImVec4& fr, int kind) {
    Capture& c = cap();
    const bool first = S.empty() && S.plots == 0;
    if (first) {
        S.x0 = fr.x; S.y0 = fr.y; S.x1 = fr.z; S.y1 = fr.w;
        S.bg = ImGui::GetColorU32(ImGuiCol_WindowBg) | 0xFF000000u;
        // theme colours -> print colours (first match wins: backgrounds, text, then grid / ticks / borders)
        auto map = [&](ImU32 from, ImU32 to) {
            for (uint32_t f : S.printFrom) if (f == uint32_t(from)) return false;
            S.printFrom.push_back(uint32_t(from));
            S.printTo.push_back(uint32_t((to & 0x00FFFFFFu) | (from & 0xFF000000u)));
            return true;
        };
        const ImU32 white = IM_COL32(255, 255, 255, 255), ink = IM_COL32(25, 25, 25, 255);
        map(S.bg, white);
        map(ImGui::GetColorU32(ImGuiCol_WindowBg), white);
        map(ImGui::GetColorU32(ImGuiCol_ChildBg), white);
        map(ImPlot::GetStyleColorU32(ImPlotCol_FrameBg), white);
        map(ImPlot::GetStyleColorU32(ImPlotCol_PlotBg), white);
        map(ImPlot::GetStyleColorU32(ImPlotCol_LegendBg), white);
        for (ImPlotCol k : {ImPlotCol_TitleText, ImPlotCol_InlayText, ImPlotCol_AxisText, ImPlotCol_LegendText}) map(ImPlot::GetStyleColorU32(k), ink);
        map(ImGui::GetColorU32(ImGuiCol_Text), ink);
        const ImU32 grid = ImPlot::GetStyleColorU32(ImPlotCol_AxisGrid);
        if (map(grid, IM_COL32(185, 185, 185, 255))) S.printTo.back() |= 0xFF000000u;   // grid: opaque light grey on paper
        map(ImPlot::GetStyleColorU32(ImPlotCol_AxisTick), IM_COL32(70, 70, 70, 255));
        map(ImPlot::GetStyleColorU32(ImPlotCol_PlotBorder), IM_COL32(60, 60, 60, 255));
        map(ImPlot::GetStyleColorU32(ImPlotCol_LegendBorder), IM_COL32(120, 120, 120, 255));
        map(ImGui::GetColorU32(ImGuiCol_Border), IM_COL32(120, 120, 120, 255));
    } else {
        S.x0 = std::min(S.x0, double(fr.x)); S.y0 = std::min(S.y0, double(fr.y));
        S.x1 = std::max(S.x1, double(fr.z)); S.y1 = std::max(S.y1, double(fr.w));
    }
    if (kind == 0) ++S.plots;
    SnapBuilder B{S};
    const ImVec2 white = dl->_Data->TexUvWhitePixel;
    auto isWhite = [&](const ImDrawVert& v) { return v.uv.x == white.x && v.uv.y == white.y; };
    std::vector<PText>& T = c.texts;
    size_t ti = 0;
    auto emitText = [&](const PText& p) {
        const ImVec4& cl = p.clip;
        if (cl.z <= cl.x || cl.w <= cl.y) return;
        fig::SnapText t = p.t;
        t.clip = B.clipIndex(cl, fr);
        t.order = int(S.runs.size());
        S.texts.push_back(std::move(t));
        B.brk = true;
    };
    for (const ImDrawCmd& cmd : dl->CmdBuffer) {
        if (cmd.UserCallback) continue;
        const int a = std::max(int(cmd.IdxOffset), idx0), b = std::min(int(cmd.IdxOffset + cmd.ElemCount), idx1);
        if (a >= b) continue;
        const ImVec4 cl = cmd.ClipRect;
        if (cl.z <= fr.x || cl.x >= fr.z || cl.w <= fr.y || cl.y >= fr.w) continue;
        const int ci = B.clipIndex(cl, fr);
        for (int k = a; k + 2 < b; k += 3) {
            while (ti < T.size() && T[ti].idx <= k) emitText(T[ti++]);
            const unsigned i0 = cmd.VtxOffset + dl->IdxBuffer[k], i1 = cmd.VtxOffset + dl->IdxBuffer[k + 1], i2 = cmd.VtxOffset + dl->IdxBuffer[k + 2];
            if (int(std::max(i0, std::max(i1, i2))) >= dl->VtxBuffer.Size) continue;
            const ImDrawVert &v0 = dl->VtxBuffer[int(i0)], &v1 = dl->VtxBuffer[int(i1)], &v2 = dl->VtxBuffer[int(i2)];
            if (!isWhite(v0) || !isWhite(v1) || !isWhite(v2)) continue;          // glyphs / images (text comes from the hook)
            if (((v0.col | v1.col | v2.col) & IM_COL32_A_MASK) == 0) continue;
            const float bx0 = std::min(v0.pos.x, std::min(v1.pos.x, v2.pos.x)), bx1 = std::max(v0.pos.x, std::max(v1.pos.x, v2.pos.x));
            const float by0 = std::min(v0.pos.y, std::min(v1.pos.y, v2.pos.y)), by1 = std::max(v0.pos.y, std::max(v1.pos.y, v2.pos.y));
            if (bx1 < cl.x || bx0 > cl.z || by1 < cl.y || by0 > cl.w) continue;  // fully clipped
            if (v0.col == v1.col && v1.col == v2.col) B.tri(v0.pos.x, v0.pos.y, v1.pos.x, v1.pos.y, v2.pos.x, v2.pos.y, v0.col, ci);
            else B.triCol(v0.pos, v1.pos, v2.pos, v0.col, v1.col, v2.col, ci, 0);
        }
    }
    while (ti < T.size()) emitText(T[ti++]);
    T.clear();
}
void hookGeometry(ImDrawList* dl, int idx0, int idx1, const ImVec4& frame, int kind, void*) {
    Capture& c = cap();
    const ImPlot::ImPlotCaptureHooks& H = ImPlot::GetCaptureHooks();
    int target = 0;
    if (kind == 1) target = c.lastTarget;
    else if (H.Active) target = 1;
    else if (c.dashState == 2 && fig::recorder().recording() && fig::recorder().currentId() == c.dashId) target = 2;
    if (kind == 0) c.lastTarget = target;
    if (target == 0) { c.texts.clear(); return; }
    addGeometry(target == 1 ? c.snap : c.dsnap, dl, idx0, idx1, frame, kind);
}
// ---- C23 native layers (heatmaps, bands, bars, text, colour scales) for the re-drawn publication layout
void forTargets(const std::function<void(fig::Panel&)>& f, bool recorderToo = true) {
    if (ImPlot::GetCaptureHooks().Active) f(cap().panel);
    if (recorderToo) if (fig::Panel* rp = fig::recorder().openPanel()) f(*rp);
}
void hookHeatmap(const char* label, const double* v, int rows, int cols, bool colMajor, double smin, double smax, double x0, double y0, double x1,
                 double y1, const ImU32* lut, int nlut, void*) {
    fig::Heatmap h;
    h.label = label ? label : "";
    h.rows = rows; h.cols = cols; h.colMajor = colMajor;
    h.v.assign(v, v + size_t(rows) * size_t(cols));
    h.vmin = smin; h.vmax = smax; h.x0 = x0; h.y0 = y0; h.x1 = x1; h.y1 = y1;
    h.lut.assign(lut, lut + nlut);
    forTargets([&](fig::Panel& p) { p.heatmaps.push_back(h); });
}
void hookBand(const char* label, const double* xs, const double* a, const double* b, int n, ImU32 col, void*) {
    fig::Series s;
    s.label = cleanLabel(label);
    if (s.label.empty()) s.legend = false;
    s.x.assign(xs, xs + n); s.ylo.assign(a, a + n); s.y.assign(b, b + n);
    s.color = uint32_t(col);
    const ImPlotPlot* pl = ImPlot::GetCurrentPlot();
    s.axis = pl && pl->CurrentY == ImAxis_Y2 ? 1 : 0;
    forTargets([&](fig::Panel& p) { p.series.push_back(s); }, false);   // dashboards record their bands explicitly
}
void hookBars(const char* label, const double* xs, const double* ys, int n, double width, ImU32 fill, ImU32 line, bool, void*) {
    fig::Series s;
    s.label = cleanLabel(label);
    if (s.label.empty()) s.legend = false;
    s.x.assign(xs, xs + n); s.y.assign(ys, ys + n);
    s.color = uint32_t(fill) | 0xFF000000u; s.barWidth = width; s.barLine = uint32_t(line);
    forTargets([&](fig::Panel& p) { p.series.push_back(s); }, false);
}
void hookTextItem(const char* text, double x, double y, float ox, float oy, ImU32 col, bool vertical, void*) {
    if (!text || !*text) return;
    fig::Shape t;
    t.kind = fig::Shape::Text; t.text = cleanLabel(text); t.x = {x}; t.y = {y}; t.offX = ox; t.offY = oy; t.stroke = uint32_t(col);
    t.vertical = vertical; t.halign = 1;
    if (t.text.empty()) return;
    forTargets([&](fig::Panel& p) { p.shapes.push_back(t); });
}
void hookColorScale(const char* label, double mn, double mx, const ImU32* lut, int nlut, void*) {
    Capture& c = cap();
    fig::ColorBar b;
    b.on = true; b.label = cleanLabel(label); b.vmin = mn; b.vmax = mx; b.lut.assign(lut, lut + nlut);
    if (c.pendingFig && c.lastTarget == 1 && !c.fig.panels.empty()) c.fig.panels.back().cbar = b;   // scale drawn right after the captured plot
    if (fig::Panel* rp = fig::recorder().lastPanel()) if (!rp->heatmaps.empty() && !rp->cbar.on) rp->cbar = b;
}
void hookMenu(ImGuiID id, const char*, void*) {
    if (ImGui::MenuItem("Export this plot (PDF / SVG / PNG / CSV)...")) ImPlot::GetCaptureHooks().Target = id;
}
} // namespace

void installPlotCapture(std::function<void(const fig::Figure&)> onCaptured) {
    cap().done = std::move(onCaptured);
    ImPlot::ImPlotCaptureHooks& H = ImPlot::GetCaptureHooks();
    H.Begin = hookBegin; H.Series = hookSeries; H.Axis = hookAxis; H.End = hookEnd; H.ContextMenu = hookMenu;
    H.Geometry = hookGeometry;
    H.Heatmap = hookHeatmap; H.Band = hookBand; H.Bars = hookBars; H.Text = hookTextItem; H.ColorScale = hookColorScale;
    H.HoverButton = true;
    GImDrawListTextHook = hookText;
}

void plotCaptureNewFrame() {
    Capture& c = cap();
    ImPlot::ImPlotCaptureHooks& H = ImPlot::GetCaptureHooks();
    c.lastTarget = 0;
    c.texts.clear();
    GImDrawListTextHookList = nullptr;
    if (c.pendingFig) {
        c.pendingFig = false;
        c.fig.snap = std::move(c.snap);
        c.snap = fig::Snapshot();
        if (c.done) c.done(c.fig);
    }
    if (c.dashState == 2) {
        H.GeomAll = false;
        c.dashState = 0;
        fig::Snapshot s = std::move(c.dsnap);
        c.dsnap = fig::Snapshot();
        if (c.dashDone) c.dashDone(c.dashId, s);
    } else if (c.dashState == 1) {
        c.dashState = 2;
        c.dsnap = fig::Snapshot();
        H.GeomAll = true;
    }
}
void requestFigureSnapshot(const std::string& figId, std::function<void(const std::string&, const fig::Snapshot&)> onDone) {
    Capture& c = cap();
    c.dashId = figId;
    c.dashDone = std::move(onDone);
    c.dashState = 1;
}
ImGuiID lastCapturedPlot() { return cap().plotId; }
void recapturePlot(ImGuiID id) { if (id) ImPlot::GetCaptureHooks().Target = id; }

// ---------------------------------------------------------------- C23 native layer drawing helpers
namespace pk {
namespace {
void add(fig::Shape&& s) {
    if (ImPlot::GetCaptureHooks().Active) cap().panel.shapes.push_back(s);
    if (fig::Panel* rp = fig::recorder().openPanel()) rp->shapes.push_back(std::move(s));
}
ImVec2 P(double x, double y) { return ImPlot::PlotToPixels(x, y); }
}
bool recording() { return ImPlot::GetCaptureHooks().Active || fig::recorder().openPanel() != nullptr; }
void line(double x0, double y0, double x1, double y1, ImU32 col, float th, const char* legend, int dash) {
    ImPlot::GetPlotDrawList()->AddLine(P(x0, y0), P(x1, y1), col, th);
    if (!recording()) return;
    fig::Shape s; s.kind = fig::Shape::Polyline; s.x = {x0, x1}; s.y = {y0, y1}; s.stroke = uint32_t(col); s.widthPx = th; s.dash = dash;
    if (legend) s.label = legend;
    add(std::move(s));
}
void polyline(const double* x, const double* y, int n, ImU32 col, float th, const char* legend) {
    if (n < 2) return;
    std::vector<ImVec2> pts(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) pts[size_t(i)] = P(x[i], y[i]);
    ImPlot::GetPlotDrawList()->AddPolyline(pts.data(), n, col, 0, th);
    if (!recording()) return;
    fig::Shape s; s.kind = fig::Shape::Polyline; s.x.assign(x, x + n); s.y.assign(y, y + n); s.stroke = uint32_t(col); s.widthPx = th;
    if (legend) s.label = legend;
    add(std::move(s));
}
void rect(double x0, double y0, double x1, double y1, ImU32 fill, ImU32 stroke, float th, const char* legend) {
    const ImVec2 a = P(x0, y0), b = P(x1, y1);
    const ImVec2 lo(std::min(a.x, b.x), std::min(a.y, b.y)), hi(std::max(a.x, b.x), std::max(a.y, b.y));
    ImDrawList* dl = ImPlot::GetPlotDrawList();
    if (fill & IM_COL32_A_MASK) dl->AddRectFilled(lo, hi, fill);
    if (stroke & IM_COL32_A_MASK) dl->AddRect(lo, hi, stroke, 0, 0, th);
    if (!recording()) return;
    fig::Shape s; s.kind = fig::Shape::Polygon; s.x = {x0, x1, x1, x0}; s.y = {y0, y0, y1, y1};
    s.fill = uint32_t(fill); s.stroke = uint32_t(stroke); s.widthPx = (stroke & IM_COL32_A_MASK) ? th : 0;
    if (legend) s.label = legend;
    add(std::move(s));
}
void quad(const double* x4, const double* y4, ImU32 fill, ImU32 stroke, float th) {
    const ImVec2 p0 = P(x4[0], y4[0]), p1 = P(x4[1], y4[1]), p2 = P(x4[2], y4[2]), p3 = P(x4[3], y4[3]);
    ImDrawList* dl = ImPlot::GetPlotDrawList();
    if (fill & IM_COL32_A_MASK) dl->AddQuadFilled(p0, p1, p2, p3, fill);
    if (stroke & IM_COL32_A_MASK) dl->AddQuad(p0, p1, p2, p3, stroke, th);
    if (!recording()) return;
    fig::Shape s; s.kind = fig::Shape::Polygon; s.x.assign(x4, x4 + 4); s.y.assign(y4, y4 + 4);
    s.fill = uint32_t(fill); s.stroke = uint32_t(stroke); s.widthPx = (stroke & IM_COL32_A_MASK) ? th : 0;
    add(std::move(s));
}
void marker(double x, double y, float r, ImU32 fill, ImU32 stroke, int shape, float th, const char* legend) {
    const ImVec2 c = P(x, y);
    ImDrawList* dl = ImPlot::GetPlotDrawList();
    if (shape == 1) {
        if (fill & IM_COL32_A_MASK) dl->AddRectFilled(ImVec2(c.x - r, c.y - r), ImVec2(c.x + r, c.y + r), fill);
        if (stroke & IM_COL32_A_MASK) dl->AddRect(ImVec2(c.x - r, c.y - r), ImVec2(c.x + r, c.y + r), stroke, 0, 0, th);
    } else {
        if (fill & IM_COL32_A_MASK) dl->AddCircleFilled(c, r, fill);
        if (stroke & IM_COL32_A_MASK) dl->AddCircle(c, r, stroke, 0, th);
    }
    if (!recording()) return;
    fig::Shape s; s.kind = fig::Shape::Marker; s.x = {x}; s.y = {y}; s.fill = uint32_t(fill); s.stroke = uint32_t(stroke);
    s.sizePx = r; s.marker = shape; s.widthPx = th; s.keepColor = true;
    if (legend) s.label = legend;
    add(std::move(s));
}
void text(double x, double y, const char* str, ImU32 col, ImVec2 off, int halign) {
    if (!str || !*str) return;
    const ImVec2 sz = ImGui::CalcTextSize(str);
    ImVec2 p = P(x, y);
    p.x += off.x - (halign == 1 ? sz.x * 0.5f : halign == 2 ? sz.x : 0.0f);
    p.y += off.y - sz.y * 0.5f;
    ImPlot::GetPlotDrawList()->AddText(p, col, str);
    if (!recording()) return;
    fig::Shape s; s.kind = fig::Shape::Text; s.x = {x}; s.y = {y}; s.text = str; s.stroke = uint32_t(col); s.offX = off.x; s.offY = off.y;
    s.halign = halign;
    add(std::move(s));
}
void colorScale(const char* label, double lo, double hi, const ImU32* lut, int n) {
    if (!recording() || n <= 0) return;
    fig::ColorBar b; b.on = true; b.label = label ? label : ""; b.vmin = lo; b.vmax = hi; b.lut.assign(lut, lut + n);
    if (ImPlot::GetCaptureHooks().Active) cap().panel.cbar = b;
    if (fig::Panel* rp = fig::recorder().openPanel()) rp->cbar = b;
}
void cellEdges(bool on) {
    if (ImPlot::GetCaptureHooks().Active) cap().panel.cellEdges = on;
    if (fig::Panel* rp = fig::recorder().openPanel()) rp->cellEdges = on;
}
} // namespace pk

// ============================================================================ C40 linked time cursor
namespace {
double g_tcDays = NAN, g_tcEpoch = 0, g_tcDragged = NAN;
bool g_tcOn = false;
bool g_tcDragFlag = false;

void timeCursorPreEnd(void*) {
    if (!g_tcOn || !std::isfinite(g_tcDays)) return;
    if (fig::recorder().recording() || ImPlot::GetCaptureHooks().Active) return;   // never in exported figures
    ImPlotPlot* pl = ImPlot::GetCurrentPlot();
    if (!pl) return;
    const ImPlotAxis& ax = pl->Axes[ImAxis_X1];
    double f = 0, off = 0;
    if (ax.Scale == ImPlotScale_Time) { f = 86400.0; off = g_tcEpoch; if (g_tcEpoch <= 0) return; }
    else {
        if (!ax.HasLabel()) return;
        std::string l = pl->GetAxisLabel(ax);
        for (char& c : l) c = char(std::tolower((unsigned char)c));
        if (l.find("distance") != std::string::npos || l.find("lag") != std::string::npos || l.find("min)") != std::string::npos) return;
        if (l.find("day") != std::string::npos) f = 1.0;
        else if (l.find("year") != std::string::npos) f = 1.0 / 365.25;
        else if (l.find("month") != std::string::npos) f = 1.0 / 30.4375;
        else return;
    }
    double x = off + g_tcDays * f;
    const ImPlotRange& r = ax.Range;
    if (x < std::min(r.Min, r.Max) || x > std::max(r.Min, r.Max)) return;
    ImVec4 col = ImGui::GetStyleColorVec4(ImGuiCol_CheckMark);
    col.w = 0.85f;
    const int id = 0x7C40 + int(pl->ID & 0x3FF);
    if (ImPlot::DragLineX(id, &x, col, 1.5f, ImPlotDragToolFlags_NoFit)) { g_tcDragged = (x - off) / f; g_tcDragFlag = true; }
    const float px = ImPlot::PlotToPixels(x, 0).x;
    if (ImGui::IsMouseHoveringRect(ImVec2(px - 4, pl->PlotRect.Min.y), ImVec2(px + 4, pl->PlotRect.Max.y)))
        ImGui::SetTooltip("Report time %.1f days - drag to move all views", g_tcDays);
}
}  // namespace

void timeCursorSet(double days, double startEpochSec, bool enabled) {
    g_tcDays = days; g_tcEpoch = startEpochSec; g_tcOn = enabled;
    ImPlot::GetCaptureHooks().PreEnd = timeCursorPreEnd;
}
bool timeCursorDragged(double* days) {
    if (!g_tcDragFlag) return false;
    g_tcDragFlag = false;
    if (days) *days = g_tcDragged;
    return std::isfinite(g_tcDragged);
}
