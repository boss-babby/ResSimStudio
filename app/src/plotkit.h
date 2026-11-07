// Shared helpers of the result plots: distinct series colours (selectable palette), plot grids that fit the
// available space, compact SI tick labels, and recording of the drawn series for publication export (figure.h).
#pragma once
#include "imgui.h"
#include "implot.h"
#include "figure.h"
#include <functional>
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

// Single-plot export for every ImPlot plot of the app (C14): right-click > "Export this plot..." or the export button
// in the corner of a hovered plot captures that plot's visible series, axis labels and ranges in the next frame
// (ImPlot capture hooks) and hands them over as a one-panel figure.
void installPlotCapture(std::function<void(const fig::Figure&)> onCaptured);
// C22 exact ("as on screen") export: the capture also records the plot's draw geometry (grid, heatmaps, cell
// backgrounds, fills, custom draw-list markers, legend, labels; a colour scale drawn right after it) as fig::Snapshot.
// Call plotCaptureNewFrame() once per frame before the UI; captured figures are delivered from there.
void plotCaptureNewFrame();
// Snapshot of all plots recorded under a fig::recorder() figure id during the next frame (dashboards).
void requestFigureSnapshot(const std::string& figId, std::function<void(const std::string&, const fig::Snapshot&)> onDone);
ImGuiID lastCapturedPlot();       // ImPlot id of the last single-plot capture
void recapturePlot(ImGuiID id);   // capture that plot again in the next frame

// C40 linked time cursor: a vertical line at the current report time on every plot whose x axis is time (label with
// "day" / "year" / "month", or a date axis). Dragging it moves the report step of all views. The App sets the cursor
// before drawing (days since START; NaN hides it) and reads a drag back after the frame.
void timeCursorSet(double days, double startEpochSec, bool enabled);
bool timeCursorDragged(double* days);   // true once after the user dragged the cursor to 'days'

// C23 native publication layers: draw on the current plot (plot coordinates, ImPlot draw list - callers push the plot clip
// rect as before) AND record the primitive natively, so the re-drawn (journal-style) export shows cell grids, fractures,
// wells, outlines and labels with publication fonts and line widths. Light neutral colours print dark on paper.
namespace pk {
bool recording();   // a capture or a recorded figure panel is active (skip building shapes otherwise)
void line(double x0, double y0, double x1, double y1, ImU32 col, float thickPx = 1.0f, const char* legend = nullptr, int dash = 0);
void polyline(const double* x, const double* y, int n, ImU32 col, float thickPx = 1.0f, const char* legend = nullptr);
void rect(double x0, double y0, double x1, double y1, ImU32 fill, ImU32 stroke = 0, float thickPx = 1.0f, const char* legend = nullptr);
void quad(const double* x4, const double* y4, ImU32 fill, ImU32 stroke = 0, float thickPx = 1.0f);
void marker(double x, double y, float radiusPx, ImU32 fill, ImU32 stroke = 0, int shape = 0, float thickPx = 1.0f, const char* legend = nullptr);
// text centred vertically on (x, y) + pixel offset; halign 0 left, 1 centre, 2 right
void text(double x, double y, const char* s, ImU32 col, ImVec2 pixOffset = ImVec2(0, 0), int halign = 0);
void cellEdges(bool on);   // outline the cells of the heatmaps of the current plot in the export
// colour scale of custom-drawn cell maps (export only; the screen keeps its own widget): lut sampled over [lo, hi]
void colorScale(const char* label, double lo, double hi, const ImU32* lut, int n);
}
