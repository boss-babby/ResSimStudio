// C22: headless check of the exact ("as on screen") plot export: real Dear ImGui + ImPlot frames (no GPU), the
// capture hooks of plotkit.cpp, and the SVG / PDF writers of figure.cpp. Draws a plot with a heatmap, a cell grid
// and fracture markers on the plot draw list, text, a shaded band, a line and a colour scale, captures it and
// writes exact_test.svg / exact_test.pdf (+ the re-drawn layout for comparison).
//   g++ -std=c++17 -O1 -Ithird_party/imgui -Ithird_party/implot -Iapp/src tools/exact_export_test.cpp \
//       third_party/imgui/imgui*.cpp third_party/implot/implot*.cpp app/src/plotkit.cpp app/src/figure.cpp -o /tmp/sd/exact
#include "imgui.h"
#include "implot.h"
#include "implot_internal.h"
#include "plotkit.h"
#include "figure.h"
#include <cmath>
#include <cstdio>
#include <vector>

static fig::Figure g_fig;
static bool g_got = false;

static void drawPlot(bool target) {
    ImGui::SetNextWindowPos(ImVec2(0, 0));
    ImGui::SetNextWindowSize(ImVec2(900, 600));
    ImGui::Begin("Board", nullptr, ImGuiWindowFlags_NoDecoration);
    if (target) ImPlot::GetCaptureHooks().Target = ImGui::GetID("Fracture positions");
    static std::vector<double> heat;
    const int nx = 12, ny = 8;
    if (heat.empty()) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) heat.push_back(std::sin(i * 0.5) * std::cos(j * 0.6));
    if (ImPlot::BeginPlot("Fracture positions", ImVec2(700, 500))) {
        ImPlot::SetupAxes("x (m)", "y (m)");
        ImPlot::SetupAxesLimits(0, 120, 0, 80, ImGuiCond_Always);
        ImPlot::PlotHeatmap("##perm", heat.data(), ny, nx, -1, 1, nullptr, ImPlotPoint(0, 0), ImPlotPoint(120, 80));
        // C23: custom geometry through the pk:: helpers (drawn + recorded natively)
        for (int i = 0; i <= nx; ++i) pk::line(i * 10.0, 0, i * 10.0, 80, IM_COL32(255, 255, 255, 90));
        for (int j = 0; j <= ny; ++j) pk::line(0, j * 10.0, 120, j * 10.0, IM_COL32(255, 255, 255, 90));
        for (int f = 0; f < 4; ++f) {   // fracture markers
            const double x = 25 + f * 22;
            pk::line(x, 20, x, 60, IM_COL32(255, 60, 40, 255), 3.0f, "fracture");
            pk::marker(x, 40, 5, IM_COL32(255, 220, 60, 255), IM_COL32(0, 0, 0, 255), 0, 1.0f, f == 0 ? "perforation" : nullptr);
        }
        pk::rect(80, 5, 115, 15, IM_COL32(48, 50, 56, 255), IM_COL32(255, 210, 90, 255), 1.5f, "inactive");
        pk::text(80, 15, "LGR-1", IM_COL32(255, 210, 90, 255), ImVec2(2, 8), 0);
        { double bx[3] = {20, 50, 90}, by[3] = {70, 74, 72}; ImPlot::PlotBars("bars", bx, by, 3, 8.0); }
        double wx[2] = {5, 115}, wy[2] = {40, 40};
        ImPlot::PlotLine("well", wx, wy, 2);
        ImPlot::PlotText("H-1", 5, 44);
        double t[3] = {10, 60, 110}, lo[3] = {5, 8, 6}, hi[3] = {12, 15, 10};
        ImPlotSpec sp; sp.FillAlpha = 0.3f;
        ImPlot::PlotShaded("band", t, lo, hi, 3, sp);
        ImPlot::EndPlot();
    }
    ImGui::SameLine();
    ImPlot::ColormapScale("k (mD)", -1, 1, ImVec2(80, 500));
    ImGui::End();
}

static void frame(bool target) {
    plotCaptureNewFrame();
    ImGui::NewFrame();
    drawPlot(target);
    ImGui::Render();
}

int main() {
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.DisplaySize = ImVec2(1000, 700);
    io.DeltaTime = 1.0f / 60;
    io.Fonts->AddFontDefault();
    unsigned char* px = nullptr; int w = 0, h = 0;
    io.Fonts->GetTexDataAsRGBA32(&px, &w, &h);
    installPlotCapture([](const fig::Figure& f) { g_fig = f; g_got = true; });
    frame(false);
    frame(true);    // capture frame
    frame(false);   // delivery at the start of this frame
    if (!g_got) { std::printf("FAIL: nothing captured\n"); return 1; }
    const fig::Snapshot& S = g_fig.snap;
    int ntri = 0; for (const auto& r : S.runs) ntri += r.ntri;
    std::printf("panels %zu series %zu | snapshot %d plots, %.0f x %.0f px, %zu runs, %d triangles, %zu texts, %zu clips\n",
                g_fig.panels.size(), g_fig.panels.empty() ? size_t(0) : g_fig.panels[0].series.size(), S.plots, S.x1 - S.x0, S.y1 - S.y0,
                S.runs.size(), ntri, S.texts.size(), S.clips.size() / 4);
    for (const auto& t : S.texts) std::printf("  text '%s' at %.0f,%.0f size %.1f w %.1f %s\n", t.s.c_str(), t.x, t.y, t.size, t.width, t.vertical ? "V" : "");
    fig::Style st;
    st.exact = true;
    st.widthMM = 160;
    st.heightMM = st.widthMM * S.aspect();
    std::string err;
    bool ok = fig::writeSVG(g_fig, st, "exact_test.svg", &err) && fig::writePDF(g_fig, st, "exact_test.pdf", &err);
    st.printColors = true;
    ok = ok && fig::writePDF(g_fig, st, "exact_test_print.pdf", &err);
    st.exact = false; st.printColors = false; st.heightMM = 100;
    ok = ok && fig::writePDF(g_fig, st, "redrawn_test.pdf", &err);
    st.cellEdges = true;
    ok = ok && fig::writeSVG(g_fig, st, "redrawn_test.svg", &err);
    const fig::Panel& P0 = g_fig.panels.at(0);
    size_t nBand = 0, nBar = 0;
    for (const auto& se : P0.series) { nBand += !se.ylo.empty(); nBar += se.barWidth > 0; }
    std::printf("native layers: heatmaps %zu, shapes %zu, colour bar %d (%s), bands %zu, bars %zu\n", P0.heatmaps.size(), P0.shapes.size(), int(P0.cbar.on),
                P0.cbar.label.c_str(), nBand, nBar);
    const bool layersOk = P0.heatmaps.size() == 1 && P0.shapes.size() >= 30 && P0.cbar.on && nBand == 1 && nBar == 1;
    std::printf("layers %s\n", layersOk ? "OK" : "MISSING");
    ok = ok && layersOk;
    std::printf("%s %s\n", ok ? "WROTE" : "WRITE FAILED", err.c_str());
    const bool pass = ok && S.plots == 1 && ntri > 200 && S.texts.size() >= 8;
    std::printf("%s\n", pass ? "PASS" : "FAIL");
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    return pass ? 0 : 1;
}
