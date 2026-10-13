#include "plotkit.h"
#include <algorithm>
#include <cmath>
#include <cstdio>

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
}
void recEnd() {
    if (!fig::recorder().recording()) return;
    const ImPlotRect r = ImPlot::GetPlotLimits(ImAxis_X1, ImAxis_Y1);
    ImPlotRect r2 = r;
    if (g_hasY2) r2 = ImPlot::GetPlotLimits(ImAxis_X1, ImAxis_Y2);
    fig::recorder().limits(r.X.Min, r.X.Max, r.Y.Min, r.Y.Max, g_hasY2, r2.Y.Min, r2.Y.Max);
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
