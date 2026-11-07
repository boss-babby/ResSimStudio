// Workspaces: switch (Model / Deck Editor / Results), the full-window Results workspace, and the publication
// figure export window (journal size presets, typography, palettes, editable labels, live WYSIWYG preview,
// vector PDF / SVG, high-resolution PNG and CSV data).
#include "app.h"
#include "i18n.h"
#include "icons.h"
#include <cfloat>
#include "figure.h"
#include "figure_png.h"
#include "plotkit.h"
#include "theme.h"
#include "widgets.h"
#include "imgui.h"
#include "imgui_internal.h"
#include <algorithm>
#include <cfloat>
#include <cmath>

namespace {

ImVec4 colA(ImU32 c, float a) { ImVec4 v = col4(c); v.w = a; return v; }

// ------------------------------------------------------------------ figure preview on an ImDrawList
class DrawListCanvas : public fig::Canvas {
public:
    DrawListCanvas(ImDrawList* dl, ImVec2 origin, float scale) : dl_(dl), o_(origin), s_(scale) {}
    void begin(double w, double h, bool transparent) override {
        const ImVec2 a = P(0, 0), b = P(w, h);
        if (transparent) {   // checkerboard
            const float c = 8;
            dl_->AddRectFilled(a, b, IM_COL32(255, 255, 255, 255));
            for (float y = a.y; y < b.y; y += c)
                for (float x = a.x + (int((y - a.y) / c) % 2) * c; x < b.x; x += 2 * c)
                    dl_->AddRectFilled(ImVec2(x, y), ImVec2(std::min(x + c, b.x), std::min(y + c, b.y)), IM_COL32(230, 230, 230, 255));
        } else {
            dl_->AddRectFilled(a, b, IM_COL32(255, 255, 255, 255));
        }
    }
    void end() override {}
    void polyline(const double* xy, int n, uint32_t color, double w, int dash) override {
        if (n < 2) return;
        const float th = std::max(0.6f, float(w) * s_);
        if (dash == 0) {
            pts_.resize(size_t(n));
            for (int i = 0; i < n; ++i) pts_[size_t(i)] = P(xy[2 * i], xy[2 * i + 1]);
            dl_->AddPolyline(pts_.data(), n, color, ImDrawFlags_None, th);
            return;
        }
        // dash patterns in units of the line width (as in the vector writers)
        static const double pat[4][4] = {{0, 0, 0, 0}, {6, 3.5, 6, 3.5}, {1, 2.5, 1, 2.5}, {6, 2.5, 1, 2.5}};
        const double lw = std::max(0.5, w);
        int k = 0;
        double left = pat[dash][0] * lw;
        bool on = true;
        for (int i = 1; i < n; ++i) {
            double x0 = xy[2 * i - 2], y0 = xy[2 * i - 1];
            const double x1 = xy[2 * i], y1 = xy[2 * i + 1];
            double seg = std::hypot(x1 - x0, y1 - y0);
            while (seg > 1e-9) {
                const double step = std::min(seg, left);
                const double t = step / seg;
                const double xm = x0 + (x1 - x0) * t, ym = y0 + (y1 - y0) * t;
                if (on) dl_->AddLine(P(x0, y0), P(xm, ym), color, th);
                x0 = xm; y0 = ym; seg -= step; left -= step;
                if (left <= 1e-9) { k = (k + 1) % 4; left = pat[dash][k] * lw; on = !on; }
            }
        }
    }
    void fillRect(double x, double y, double w, double h, uint32_t c) override { dl_->AddRectFilled(P(x, y), P(x + w, y + h), c); }
    void strokeRect(double x, double y, double w, double h, uint32_t c, double lw) override {
        dl_->AddRect(P(x, y), P(x + w, y + h), c, 0, 0, std::max(0.6f, float(lw) * s_));
    }
    void circle(double x, double y, double r, uint32_t fill, uint32_t stroke, double lw) override {
        const float rr = std::max(1.0f, float(r) * s_);
        if (fill >> 24) dl_->AddCircleFilled(P(x, y), rr, fill);
        if (stroke >> 24) dl_->AddCircle(P(x, y), rr, stroke, 0, std::max(0.6f, float(lw) * s_));
    }
    void text(double x, double y, const std::string& sIn, double sizePt, int halign, bool vertical, uint32_t color) override {
        // superscripts "^{n}" are shown inline as "^n" in the preview
        std::string s;
        for (size_t i = 0; i < sIn.size(); ++i) {
            if (sIn[i] == '^' && i + 1 < sIn.size() && sIn[i + 1] == '{') { s += '^'; ++i; continue; }
            if (sIn[i] == '}') continue;
            if (sIn.compare(i, 3, "\xE2\x88\x92") == 0) { s += '-'; i += 2; continue; }   // U+2212 is not in the UI font
            s += sIn[i];
        }
        ImFont* f = g_fonts.ui ? g_fonts.ui : ImGui::GetFont();
        float px = std::max(1.0f, float(sizePt) * s_);
        float wpx = f->CalcTextSizeA(px, FLT_MAX, 0, s.c_str()).x;
        // match the width of the export fonts (Helvetica / Times metrics) so the preview shows the real layout
        const float want = float(fig::textWidthPt(sIn, sizePt, font)) * s_;
        if (wpx > want * 1.02f && wpx > 0) { px *= want / wpx; wpx = want; }
        const float asc = f->GetFontBaked(px)->Ascent;
        const ImVec2 at = P(x, y);
        const float off = halign == 1 ? -wpx * 0.5f : halign == 2 ? -wpx : 0.0f;
        if (!vertical) {
            dl_->AddText(f, px, ImVec2(at.x + off, at.y - asc), color, s.c_str());
            return;
        }
        // draw horizontally, then rotate the new vertices 90 deg CCW about the anchor
        const int v0 = dl_->VtxBuffer.Size;
        dl_->AddText(f, px, ImVec2(at.x + off, at.y - asc), color, s.c_str());
        for (int i = v0; i < dl_->VtxBuffer.Size; ++i) {
            ImVec2& p = dl_->VtxBuffer[i].pos;
            const float dx = p.x - at.x, dy = p.y - at.y;
            p = ImVec2(at.x + dy, at.y - dx);
        }
    }
    void fillTriangles(const double* xy, int nTri, uint32_t color) override {
        if (nTri <= 0 || (color >> 24) == 0) return;
        const ImVec2 uv = dl_->_Data->TexUvWhitePixel;
        dl_->PrimReserve(3 * nTri, 3 * nTri);
        for (int t = 0; t < nTri; ++t) {
            const ImDrawIdx base = ImDrawIdx(dl_->_VtxCurrentIdx);
            for (int k = 0; k < 3; ++k) dl_->PrimWriteVtx(P(xy[6 * t + 2 * k], xy[6 * t + 2 * k + 1]), uv, color);
            dl_->PrimWriteIdx(base); dl_->PrimWriteIdx(ImDrawIdx(base + 1)); dl_->PrimWriteIdx(ImDrawIdx(base + 2));
        }
    }
    void textFit(double x, double y, const std::string& s, double sizePt, double widthPt, bool vertical, uint32_t color) override {
        ImFont* f = g_fonts.ui ? g_fonts.ui : ImGui::GetFont();
        const float px = std::max(1.0f, float(sizePt) * s_);
        const float wpx = f->CalcTextSizeA(px, FLT_MAX, 0, s.c_str()).x, want = float(widthPt) * s_;
        const float sx = (wpx > 0 && want > 0) ? std::clamp(want / wpx, 0.5f, 2.0f) : 1.0f;
        const float asc = f->GetFontBaked(px)->Ascent;
        const ImVec2 at = P(x, y);
        const int v0 = dl_->VtxBuffer.Size;
        dl_->AddText(f, px, ImVec2(at.x, at.y - asc), color, s.c_str());
        for (int i = v0; i < dl_->VtxBuffer.Size; ++i) {
            ImVec2& p = dl_->VtxBuffer[i].pos;
            const float dx = (p.x - at.x) * sx, dy = p.y - at.y;
            p = vertical ? ImVec2(at.x + dy, at.y - dx) : ImVec2(at.x + dx, at.y + dy);
        }
    }
    void pushClip(double x, double y, double w, double h) override { dl_->PushClipRect(P(x, y), P(x + w, y + h), true); }
    void popClip() override { dl_->PopClipRect(); }

private:
    ImVec2 P(double x, double y) const { return ImVec2(o_.x + float(x) * s_, o_.y + float(y) * s_); }
    ImDrawList* dl_;
    ImVec2 o_;
    float s_;
    std::vector<ImVec2> pts_;
};

struct SizePreset { const char* name; double w, h, font, line; };
const SizePreset kPresets[] = {
    {"Custom", 0, 0, 0, 0},
    {"Journal - single column (85 x 65 mm)", 85, 65, 8, 1.0},
    {"Journal - 1.5 column (120 x 85 mm)", 120, 85, 8, 1.0},
    {"Journal - double column (174 x 110 mm)", 174, 110, 9, 1.1},
    {"Thesis / report (150 x 100 mm)", 150, 100, 10, 1.2},
    {"Slide 16:9 (254 x 143 mm)", 254, 143, 16, 2.0},
    {"Poster panel (300 x 200 mm)", 300, 200, 20, 2.4},
};
constexpr int kNumPresets = int(sizeof(kPresets) / sizeof(kPresets[0]));

ImVec4 toV4(uint32_t c) { return ImGui::ColorConvertU32ToFloat4(c); }

} // namespace

// ============================================================================ workspace switch
void App::workspaceSwitch() {
    struct Item { const char* label; Workspace w; const char* tip; };
    static const Item items[] = {
        {TRN(ICON_BOX " Model"), Workspace::Model, TRN("Model workspace: inputs, 3D viewer, run monitor (Ctrl+1)")},
        {TRN(ICON_FILE_CODE " Deck Editor"), Workspace::DeckEditor, TRN("Deck Editor - writing mode: intelligent ECLIPSE / OPM deck editing (Ctrl+2)")},
        {TRN(ICON_CHART_LINE " Results"), Workspace::Results, TRN("Results workspace: run dashboard, cases, vectors, compare, 3D, figure export (Ctrl+3)")},
        {TRN(ICON_FLASK_CONICAL " Studies"), Workspace::Studies, TRN("Studies workspace: sensitivity, uncertainty, history matching, proxy models (Ctrl+4)")},
        {TRN(ICON_MOUNTAIN " Geomodel"), Workspace::Geomodel, TRN("Geomodel workspace (ResGeo): wells / surfaces / faults / SEG-Y -> structural grid -> facies -> petrophysics -> volumetrics -> upscaling -> simulation case (Ctrl+5)")},
    };
    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(1, ImGui::GetStyle().ItemSpacing.y));
    for (int i = 0; i < 5; ++i) {
        if (geoMode_ && items[i].w != Workspace::Geomodel) continue;   // ResGeo shows the geomodelling workspace only
        const bool on = workspace_ == items[i].w;
        if (i) ImGui::SameLine();
        if (on) {
            ImGui::PushStyleColor(ImGuiCol_Button, colA(Pal::Accent, 0.85f));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colA(Pal::Accent, 0.95f));
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1, 1, 1, 1));
        }
        if (ImGui::Button(TR(items[i].label))) workspace_ = items[i].w;
        if (on) ImGui::PopStyleColor(3);
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("%s", TR(items[i].tip));
    }
    ImGui::PopStyleVar();
}

// ============================================================================ Results workspace
void App::resultsWorkspace() {
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(vp->WorkPos);
    ImGui::SetNextWindowSize(vp->WorkSize);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0);
    ImGui::Begin("##results", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoSavedSettings |
                                           ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoMove);
    ImGui::PopStyleVar(2);
    const float fs = ImGui::GetFontSize();

    // ---- header: sections + plot style
    static const char* tabs[] = {TRN("Run dashboard"), TRN("Cases, vectors & compare")};
    for (int i = 0; i < 2; ++i) {
        if (i) ImGui::SameLine();
        const bool on = resultsTab_ == i;
        if (on) ImGui::PushStyleColor(ImGuiCol_Button, colA(Pal::Accent, 0.55f));
        if (ImGui::Button(TRL(tabs[i]))) resultsTab_ = i;
        if (on) ImGui::PopStyleColor();
    }
    ImGui::SameLine(0, fs * 1.5f);
    ImGui::TextDisabled("%s", TR("Curve palette"));
    ImGui::SameLine();
    ImGui::SetNextItemWidth(fs * 9);
    const auto& P = fig::palettes();
    if (ImGui::BeginCombo("##pal", P[size_t(std::clamp(g_plotPalette, 0, int(P.size()) - 1))].name)) {
        for (int i = 0; i < int(P.size()); ++i) {
            ImGui::PushID(i);
            if (ImGui::Selectable("##p", i == g_plotPalette, 0, ImVec2(fs * 16, 0))) g_plotPalette = i;
            ImGui::SameLine(0, 0);
            ImDrawList* dl = ImGui::GetWindowDrawList();
            ImVec2 p = ImGui::GetCursorScreenPos();
            for (int k = 0; k < 8; ++k)
                dl->AddRectFilled(ImVec2(p.x + k * fs * 0.8f, p.y + fs * 0.15f), ImVec2(p.x + k * fs * 0.8f + fs * 0.7f, p.y + fs * 0.95f), fig::distinctColor(k, i));
            ImGui::Dummy(ImVec2(fs * 6.6f, fs));
            ImGui::SameLine();
            ImGui::TextUnformatted(P[size_t(i)].name);
            ImGui::PopID();
        }
        ImGui::EndCombo();
    }
    ImGui::SameLine();
    ImGui::TextDisabled("%s", TR("Line width"));
    ImGui::SameLine();
    ImGui::SetNextItemWidth(fs * 6);
    ImGui::SliderFloat("##lw", &g_plotLineScale, 0.5f, 3.0f, "%.1fx");
    if (resultsTab_ == 1) {
        ImGui::SameLine(0, fs * 1.5f);
        if (ImGui::Button("+ Current run")) rvAddRun();
        ImGui::SameLine();
        if (ImGui::Button("Open results...")) {
            std::string p;
            if (fileDialog(false, "ECLIPSE / OPM Flow results\0*.SMSPEC;*.FSMSPEC;*.EGRID;*.FEGRID;*.UNRST;*.FUNRST;*.INIT;*.FINIT;*.UNSMRY;*.FUNSMRY;*.DATA;*.data\0All files (*.*)\0*.*\0", "SMSPEC", p))
                rvOpen(p);
        }
    }
    ImGui::Separator();
    ImGui::BeginChild("##resbody", ImVec2(0, 0), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
    if (resultsTab_ == 0) plotsBody();
    else rvBody();
    ImGui::EndChild();
    ImGui::End();
}

// ============================================================================ figure export
void App::exportButton(const char* figId) {
    ImGui::PushID(figId);
    if (ImGui::SmallButton("Export figure...")) {
        if (const fig::Figure* f = fig::recorder().get(figId)) {
            figx_.fig = *f;
            figx_.id = figId;
            figx_.open = true;
            figx_.cols = 0;
            figx_.plotId = 0;
            figx_.snapWait = false;
            figx_.snapTried = false;
            if (figx_.preset > 0) {
                const SizePreset& p = kPresets[figx_.preset];
                figx_.style.widthMM = p.w; figx_.style.heightMM = p.h; figx_.style.fontPt = p.font; figx_.style.lineWidthPt = p.line;
            }
            ImGui::SetWindowFocus("Export figure");
        } else {
            notify("Nothing to export yet - the plot has not been drawn");
        }
    }
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort))
        ImGui::SetTooltip("Publication-quality export of these plots: vector PDF / SVG, high-resolution PNG, CSV data.\n"
                          "Journal sizes, fonts, line widths, palettes and editable labels with a live preview.");
    ImGui::PopID();
}

void App::figureExportWindow() {
    if (!figx_.open) return;
    const float fs = ImGui::GetFontSize();
    ImGui::SetNextWindowSize(ImVec2(fs * 66, fs * 40), ImGuiCond_FirstUseEver);
    if (figDeferredFocus_) { ImGui::SetNextWindowFocus(); figDeferredFocus_ = false; }
    if (!ImGui::Begin("Export figure", &figx_.open, ImGuiWindowFlags_NoDocking)) { ImGui::End(); return; }
    fig::Style& st = figx_.style;
    fig::Figure& F0 = figx_.fig;

    // working copy with the layout override
    fig::Figure F;   // (without the snapshot - exact mode draws F0 directly, no per-frame copy of the geometry)
    F.name = F0.name; F.panels = F0.panels; F.rows = F0.rows; F.cols = F0.cols;
    const int n = int(F.panels.size());
    if (figx_.cols > 0 && n > 0) { F.cols = std::min(figx_.cols, n); F.rows = (n + F.cols - 1) / F.cols; }

    // C22: exact snapshot of a recorded (dashboard) figure - requested on demand, arrives one frame later
    auto requestSnap = [this]() {
        if (figx_.id.empty()) return;
        figx_.snapWait = true;
        requestFigureSnapshot(figx_.id, [this](const std::string& id, const fig::Snapshot& sn) {
            if (id != figx_.id) return;
            figx_.snapWait = false;
            figx_.snapTried = true;
            figx_.fig.snap = sn;
        });
    };
    if (st.exact && F0.snap.empty() && !figx_.id.empty() && !figx_.snapWait && !figx_.snapTried) requestSnap();
    const bool exact = st.exact && !F0.snap.empty();
    if (exact) {   // the page has the on-screen aspect ratio
        st.heightMM = std::clamp(st.widthMM * F0.snap.aspect(), 10.0, 900.0);
        const double pxW = F0.snap.x1 - F0.snap.x0;
        if (figx_.pxScale > 0 && pxW > 0) st.dpi = std::clamp(int(std::lround(figx_.pxScale * pxW / (st.widthMM / 25.4))), 72, 2400);
    }

    ImGui::BeginChild("##figset", ImVec2(fs * 21, 0), ImGuiChildFlags_Borders | ImGuiChildFlags_ResizeX);
    ImGui::PushItemWidth(-fs * 6.5f);
    ImGui::SeparatorText("Mode");
    {
        int mode = st.exact ? 0 : 1;
        if (ImGui::RadioButton("Publication (journal style)", &mode, 1)) st.exact = false;
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort))
            ImGui::SetTooltip("Re-drawn natively with journal fonts, line widths, palettes, ticks and legends - including heatmaps / cell maps\n"
                              "with a colour bar, model-cell outlines, filled bands, bars, fractures, wells, outlines and labels (C23).\n"
                              "Light lines and text made for the dark screen theme print dark on paper.");
        if (ImGui::RadioButton("Exact (as on screen)", &mode, 0)) { st.exact = true; figx_.snapTried = false; }
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort))
            ImGui::SetTooltip("A copy of the plot exactly as the board shows it (screen fonts and theme): grid, cell backgrounds, markers,\n"
                              "legend, axes and colour scale as vector geometry (PDF / SVG) or a sharp raster (PNG).");
        if (st.exact) {
            if (figx_.snapWait) ImGui::TextColored(ImVec4(0.9f, 0.75f, 0.3f, 1), "Capturing the plots from the screen...");
            else if (F0.snap.empty()) {
                ImGui::TextWrapped("No on-screen snapshot yet: the plots must be visible. Show them and press Capture.");
            } else {
                ImGui::TextDisabled("%d plot(s), %.0f x %.0f px, %d text items", std::max(1, F0.snap.plots), F0.snap.x1 - F0.snap.x0,
                                    F0.snap.y1 - F0.snap.y0, int(F0.snap.texts.size()));
            }
            if (ImGui::SmallButton(F0.snap.empty() ? "Capture from screen" : "Re-capture from screen")) {
                if (!figx_.id.empty()) { figx_.snapTried = false; requestSnap(); }
                else if (figx_.plotId) recapturePlot(figx_.plotId);
            }
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("Take the plots again as they are now on the board (zoom, colours, markers)");
            ImGui::Checkbox("Print colours (light background)", &st.printColors);
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort))
                ImGui::SetTooltip("Inverts the neutral greys of the dark theme (background, text, grid) - series and colour maps keep their colours");
        }
    }
    ImGui::SeparatorText("Size");
    if (ImGui::BeginCombo("Preset", kPresets[figx_.preset].name)) {
        for (int i = 0; i < kNumPresets; ++i)
            if (ImGui::Selectable(kPresets[i].name, i == figx_.preset)) {
                figx_.preset = i;
                if (i > 0) { st.widthMM = kPresets[i].w; st.heightMM = kPresets[i].h; st.fontPt = kPresets[i].font; st.lineWidthPt = kPresets[i].line; }
            }
        ImGui::EndCombo();
    }
    if (ImGui::InputDouble("Width (mm)", &st.widthMM, 1, 10, "%.0f")) figx_.preset = 0;
    ImGui::BeginDisabled(exact);
    if (ImGui::InputDouble("Height (mm)", &st.heightMM, 1, 10, "%.0f")) figx_.preset = 0;
    ImGui::EndDisabled();
    if (exact && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) ImGui::SetTooltip("Exact mode keeps the on-screen aspect ratio");
    st.widthMM = std::clamp(st.widthMM, 30.0, 600.0);
    st.heightMM = std::clamp(st.heightMM, 20.0, 600.0);
    int cols = figx_.cols;
    if (n > 1 && ImGui::SliderInt("Columns", &cols, 0, std::min(n, 6), cols == 0 ? "as in app" : "%d")) figx_.cols = cols;

    ImGui::BeginDisabled(exact);   // the publication-layout settings do not apply to an exact snapshot
    ImGui::SeparatorText("Typography & lines");
    static const char* fonts[] = {"Helvetica / Arial (sans)", "Times (serif)"};
    ImGui::Combo("Font", &st.font, fonts, 2);
    float fpt = float(st.fontPt), lpt = float(st.lineWidthPt);
    if (ImGui::SliderFloat("Font (pt)", &fpt, 5, 28, "%.1f")) { st.fontPt = fpt; figx_.preset = 0; }
    if (ImGui::SliderFloat("Line (pt)", &lpt, 0.3f, 4, "%.2f")) { st.lineWidthPt = lpt; figx_.preset = 0; }
    const auto& P = fig::palettes();
    std::string palName = st.palette < 0 ? "Same as app plots" : P[size_t(std::clamp(st.palette, 0, int(P.size()) - 1))].name;
    if (ImGui::BeginCombo("Colours", palName.c_str())) {
        if (ImGui::Selectable("Same as app plots", st.palette < 0)) st.palette = -1;
        for (int i = 0; i < int(P.size()); ++i) if (ImGui::Selectable(P[size_t(i)].name, st.palette == i)) st.palette = i;
        ImGui::EndCombo();
    }
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("Okabe-Ito and Tol palettes are colour-blind safe; 'Grayscale' for print.");

    ImGui::SeparatorText("Axes & legend");
    static const char* legends[] = {"Auto", "Inside top-right", "Inside top-left", "Inside bottom-right", "Inside bottom-left", "Outside right", "Below", "None"};
    int lg = int(st.legend);
    if (ImGui::Combo("Legend", &lg, legends, 8)) st.legend = fig::Legend(lg);
    ImGui::Checkbox("Grid", &st.grid);
    ImGui::SameLine(); ImGui::Checkbox("Minor ticks", &st.minorTicks);
    ImGui::Checkbox("Box", &st.box);
    ImGui::SameLine(); ImGui::Checkbox("Panel titles", &st.titles);
    ImGui::Checkbox("Use visible axis ranges", &st.useLimits);
    ImGui::Checkbox("Maps, cells & markers", &st.layers);
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort))
        ImGui::SetTooltip("Heatmaps / cell maps with a colour bar, filled bands, bars, fractures, wells, outlines and text labels");
    ImGui::SameLine(); ImGui::Checkbox("Cell outlines", &st.cellEdges);
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("Thin outlines around every heatmap cell (the model grid)");
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("On: the zoom of the app plots. Off: fit all data with nice round limits.");

    ImGui::EndDisabled();
    ImGui::SeparatorText("Raster (PNG)");
    static const int dpis[] = {150, 300, 600, 1200};
    if (exact) {
        static const char* sc[] = {"Use dpi below", "1x screen pixels", "2x screen pixels", "3x", "4x"};
        ImGui::Combo("Pixel scale", &figx_.pxScale, sc, 5);
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("1x = the same number of pixels as on screen; 2x-4x = sharper copies of the same picture");
    }
    int di = 2;
    for (int i = 0; i < 4; ++i) if (st.dpi == dpis[i]) di = i;
    static const char* dpiNames[] = {"150 dpi (screen)", "300 dpi (print)", "600 dpi (line art)", "1200 dpi"};
    ImGui::BeginDisabled(exact && figx_.pxScale > 0);
    if (ImGui::Combo("Resolution", &di, dpiNames, 4)) st.dpi = dpis[di];
    ImGui::EndDisabled();
    ImGui::Checkbox("Transparent background", &st.transparent);
    ImGui::TextDisabled("%d x %d px", int(st.widthMM / 25.4 * st.dpi + 0.5), int(st.heightMM / 25.4 * st.dpi + 0.5));
    ImGui::PopItemWidth();

    ImGui::SeparatorText("Labels");
    if (exact) ImGui::TextDisabled("(label edits apply to the re-drawn mode)");
    for (int i = 0; i < int(F0.panels.size()); ++i) {
        fig::Panel& p = F0.panels[size_t(i)];
        ImGui::PushID(i);
        const std::string h = (p.title.empty() ? std::string("Panel ") + std::to_string(i + 1) : p.title) + "###panel";
        if (ImGui::TreeNode(h.c_str())) {
            ImGui::PushItemWidth(-fs * 4);
            ui::inputString("Title", p.title);
            ui::inputString("X axis", p.xlabel);
            ui::inputString("Y axis", p.ylabel);
            if (p.hasY2 || !p.y2label.empty()) ui::inputString("Y2 axis", p.y2label);
            ImGui::PopItemWidth();
            ImGui::TextDisabled("Series");
            for (int k = 0; k < int(p.series.size()); ++k) {
                fig::Series& s = p.series[size_t(k)];
                ImGui::PushID(k);
                ImVec4 c = toV4(s.color);
                if (ImGui::ColorEdit3("##c", &c.x, ImGuiColorEditFlags_NoInputs)) s.color = ImGui::ColorConvertFloat4ToU32(c);
                ImGui::SameLine();
                static const char* dashes[] = {"solid", "dashed", "dotted", "dash-dot"};
                ImGui::SetNextItemWidth(fs * 4.5f);
                ImGui::Combo("##d", &s.dash, dashes, 4);
                ImGui::SameLine();
                ImGui::Checkbox("##l", &s.legend);
                if (ImGui::IsItemHovered()) ImGui::SetTooltip("Show in legend");
                ImGui::SameLine();
                ImGui::SetNextItemWidth(-1);
                ui::inputString("##lab", s.label);
                ImGui::PopID();
            }
            ImGui::TreePop();
        }
        ImGui::PopID();
    }
    ImGui::BeginDisabled(figx_.id.empty());
    if (ImGui::SmallButton("Reload from plots")) {
        if (const fig::Figure* f = fig::recorder().get(figx_.id)) { figx_.fig = *f; figx_.snapTried = false; }
    }
    ImGui::EndDisabled();
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("Take the current state of the app plots again (discards label edits)");
    ImGui::EndChild();

    // ---- preview + export buttons
    ImGui::SameLine();
    ImGui::BeginGroup();
    auto save = [&](int kind) {
        static const char* filters[] = {"PDF (vector)\0*.pdf\0", "SVG (vector)\0*.svg\0", "PNG image\0*.png\0", "CSV data\0*.csv\0"};
        static const char* exts[] = {"pdf", "svg", "png", "csv"};
        std::string base = figx_.lastDir.empty() ? std::string() : figx_.lastDir + "\\";
        std::string name = figx_.id;
        std::replace(name.begin(), name.end(), '.', '_');
        std::string path = base + name + "." + exts[kind];
        if (!fileDialog(true, filters[kind], exts[kind], path)) return;
        const size_t sl = path.find_last_of("\\/");
        if (sl != std::string::npos) figx_.lastDir = path.substr(0, sl);
        std::string err;
        bool ok = false;
        const fig::Figure& FX = exact ? F0 : F;
        if (kind == 0) ok = fig::writePDF(FX, st, path, &err);
        else if (kind == 1) ok = fig::writeSVG(FX, st, path, &err);
        else if (kind == 2) ok = fig::writePNG(FX, st, path, &err);
        else ok = fig::writeCSV(F, path, &err);
        notify(ok ? "Figure exported: " + path : "Export failed: " + err);
    };
    if (ui::accentButton("PDF (vector)...")) save(0);
    ImGui::SameLine();
    if (ImGui::Button("SVG (vector)...")) save(1);
    ImGui::SameLine();
    if (ImGui::Button("PNG...")) save(2);
    ImGui::SameLine();
    if (ImGui::Button("CSV data...")) save(3);
    ImGui::SameLine();
    if (exact) ImGui::TextDisabled("  %.0f x %.0f mm, exact, PNG %d x %d px", st.widthMM, st.heightMM, int(st.widthMM / 25.4 * st.dpi + 0.5), int(st.heightMM / 25.4 * st.dpi + 0.5));
    else ImGui::TextDisabled("  %.0f x %.0f mm, %d panel(s)", st.widthMM, st.heightMM, n);

    const ImVec2 av = ImGui::GetContentRegionAvail();
    ImGui::BeginChild("##figprev", av, ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar);
    {
        ImDrawList* dl = ImGui::GetWindowDrawList();
        const ImVec2 p0 = ImGui::GetCursorScreenPos();
        dl->AddRectFilled(p0, ImVec2(p0.x + av.x, p0.y + av.y), IM_COL32(60, 63, 70, 255));
        const double wPt = st.widthMM / 25.4 * 72, hPt = st.heightMM / 25.4 * 72;
        const float pad = fs;
        const float sc = float(std::min((av.x - 2 * pad) / wPt, (av.y - 2 * pad) / hPt));
        if (sc > 0) {
            const ImVec2 o(p0.x + (av.x - float(wPt) * sc) * 0.5f, p0.y + (av.y - float(hPt) * sc) * 0.5f);
            dl->AddRectFilled(ImVec2(o.x + 4, o.y + 4), ImVec2(o.x + float(wPt) * sc + 4, o.y + float(hPt) * sc + 4), IM_COL32(0, 0, 0, 90));
            dl->PushClipRect(o, ImVec2(o.x + float(wPt) * sc, o.y + float(hPt) * sc), true);
            DrawListCanvas cv(dl, o, sc);
            fig::render(exact ? F0 : F, st, cv);
            dl->PopClipRect();
        }
        ImGui::Dummy(av);
    }
    ImGui::EndChild();
    ImGui::EndGroup();
    ImGui::End();
}
