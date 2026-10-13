// Shared helpers of the result plots: distinct series colours (selectable palette), plot grids that fit the
// available space, compact SI tick labels, and recording of the drawn series for publication export (figure.h).
#pragma once
#include "imgui.h"
#include "implot.h"
#include "figure.h"
#include <string>

extern int g_plotPalette;      // fig::palettes() index used for curve colours
extern float g_plotLineScale;  // line width multiplier of the app plots

inline ImU32 seriesColor(int i) { return ImU32(fig::distinctColor(i, g_plotPalette)); }

struct PlotGrid { int rows = 1, cols = 1; ImVec2 size; bool scroll = false; };
// rows x cols for n plots in 'avail': plots as large and close to 'aspect' (w/h) as possible, at least min size;
// size = total subplot size (taller than avail -> scroll == true, draw inside a scrolling child)
PlotGrid fitGrid(int n, ImVec2 avail, float minW, float minH, float aspect = 1.7f, int forceCols = 0);

// ImPlot tick formatter: 1.25k, 3.4M, 2G ...
int formatSI(double value, char* buf, int size, void* user);
void setupSI(ImAxis axis);

// recording (no-ops when fig::recorder() is not recording)
void recPanel(const char* title, const char* xlabel, const char* ylabel, const char* y2label = "");
void recEnd();                 // call before ImPlot::EndPlot(): stores visible limits
void recSeries(const char* label, const double* x, const double* y, int n, ImU32 col, float weight, bool scatter = false, bool markers = false);
void recVLine(double x);
void useY2(bool on);           // ImPlot::SetAxes(X1, Y2 / Y1) + recorder axis
