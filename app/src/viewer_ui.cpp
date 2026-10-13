// Shared 3D viewer UI for the Model viewer and the Results viewer: one toolbar (colour map, range, display,
// filters / slices, transparency, well focus, streamlines, camera presets, PNG snapshot) and one canvas
// (orbit / pan / zoom, cell probe, well picking, context menu, cell value history).
#include "app.h"
#include "widgets.h"
#include "theme.h"
#include "icons.h"
#include "plotkit.h"
#include "implot.h"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>

using namespace ressim;

namespace {

std::array<float, 3> cellCentre(const Grid& g, int q) {
    std::array<float, 3> c{0, 0, 0};
    for (const auto& p : g.corners[size_t(q)]) { c[0] += p.x * 0.125f; c[1] += p.y * 0.125f; c[2] += p.z * 0.125f; }
    return c;
}
// centre of the cell face normal to 'axis' (0 = I, 1 = J, 2 = K), on the + side if plus
std::array<float, 3> faceCentre(const Grid& g, int q, int axis, bool plus) {
    std::array<float, 3> c{0, 0, 0};
    const auto& cr = g.corners[size_t(q)];
    for (int k = 0; k < 8; ++k)
        if (((k >> axis) & 1) == (plus ? 1 : 0)) { c[0] += cr[size_t(k)].x * 0.25f; c[1] += cr[size_t(k)].y * 0.25f; c[2] += cr[size_t(k)].z * 0.25f; }
    return c;
}

// flow layout: continue on the same row if the item fits, else wrap
void flow(float w) {
    const float right = ImGui::GetWindowPos().x + ImGui::GetWindowContentRegionMax().x;
    ImGui::SameLine();
    if (ImGui::GetCursorScreenPos().x + w > right) ImGui::NewLine();
}
float btnW(const char* label) { return ImGui::CalcTextSize(label, nullptr, true).x + ImGui::GetStyle().FramePadding.x * 2; }
void tip(const char* t) { if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("%s", t); }

} // namespace

// ----------------------------------------------------------------------------- wells: real trajectories
std::vector<WellGlyph> App::wellGlyphs(const Grid& g, const Case& c, int step) const {
    std::vector<WellGlyph> out;
    const float h = g.bbMax.z - g.bbMin.z;
    const float zHead = g.bbMin.z - std::max(0.25f * h, 1.0f);
    auto valid = [&](int i, int j, int k) { return i >= 0 && i < g.nx && j >= 0 && j < g.ny && k >= 0 && k < g.nz; };
    for (const auto& w : c.wells) {
        WellGlyph gl;
        gl.name = w.name;
        gl.injector = w.type != WellType::Producer;
        if (step >= 0) {   // status from the latest control record up to this step
            bool any = false;
            for (int s = 0; s <= step && s < int(c.schedule.size()); ++s)
                for (const auto& wc : c.schedule[size_t(s)].controls)
                    if (wc.well == w.name) { gl.open = wc.status == WellStatus::Open; any = true; }
            if (!any) gl.open = step < 0;
        }
        struct Conn { int q; char dir; bool open; };
        std::vector<Conn> conns;
        for (const auto& cp : w.completions) {
            const bool active = cp.open && (step < 0 || (cp.fromStep <= step && (cp.untilStep < 0 || step < cp.untilStep)));
            for (int k = cp.k1; k <= cp.k2; ++k)
                if (valid(cp.i, cp.j, k)) conns.push_back({cp.i + g.nx * (cp.j + g.ny * k), cp.dir, active});
        }
        int hi = w.i, hj = w.j;
        if (!valid(hi, hj, 0)) {
            if (conns.empty()) continue;
            int k; g.ijk(conns[0].q, hi, hj, k);
        }
        const auto top = faceCentre(g, hi + g.nx * hj, 2, false);
        gl.path.push_back({top[0], top[1], zHead});
        if (conns.empty()) gl.path.push_back({top[0], top[1], g.bbMax.z});
        else {
            if (conns.front().dir != 'X' && conns.front().dir != 'Y') gl.path.push_back(faceCentre(g, conns.front().q, 2, false));
            for (const auto& cn : conns) gl.path.push_back(cellCentre(g, cn.q));
            if (conns.back().dir != 'X' && conns.back().dir != 'Y') gl.path.push_back(faceCentre(g, conns.back().q, 2, true));
            for (const auto& cn : conns) {
                const int axis = cn.dir == 'X' ? 0 : cn.dir == 'Y' ? 1 : 2;
                WellGlyph::Perf pf;
                pf.a = faceCentre(g, cn.q, axis, false);
                pf.b = faceCentre(g, cn.q, axis, true);
                pf.open = cn.open;
                gl.perfs.push_back(pf);
            }
        }
        out.push_back(std::move(gl));
    }
    return out;
}

// ----------------------------------------------------------------------------- streamline vertices
std::vector<LineVertex> streamlineVertices(const StreamlineResult& S, int colorBy, Colormap cm, const std::string& focus, bool focusOnly) {
    std::vector<LineVertex> v;
    size_t np = 0;
    for (auto& L : S.lines) np += L.pts.size();
    v.reserve(np * 2);
    const float lt = std::log10(float(S.maxTof) + 1.0f);
    for (const Streamline& L : S.lines) {
        const bool mine = !focus.empty() && ((L.source >= 0 && L.source < int(S.sources.size()) && S.sources[size_t(L.source)] == focus) ||
                                             (L.sink >= 0 && L.sink < int(S.sinks.size()) && S.sinks[size_t(L.sink)] == focus));
        if (!focus.empty() && focusOnly && !mine) continue;
        const float alpha = focus.empty() || mine ? 0.95f : 0.18f;
        float rgb[3] = {1, 1, 1};
        if (colorBy == 1 || colorBy == 2) {
            const int id = colorBy == 1 ? L.source : L.sink;
            const ImVec4 c = id < 0 ? ImVec4(0.6f, 0.6f, 0.6f, 1) : ImGui::ColorConvertU32ToFloat4(seriesColor(id + (colorBy == 2 ? 3 : 0)));
            rgb[0] = c.x; rgb[1] = c.y; rgb[2] = c.z;
        }
        for (size_t k = 1; k < L.pts.size(); ++k)
            for (int e = 0; e < 2; ++e) {
                const size_t q = k - 1 + size_t(e);
                if (colorBy == 0) colormapRGB(cm, lt > 0 ? std::log10(L.tof[q] + 1.0f) / lt : 0.0f, rgb);
                v.push_back({L.pts[q].x, L.pts[q].y, L.pts[q].z, rgb[0], rgb[1], rgb[2], alpha});
            }
    }
    return v;
}

// ----------------------------------------------------------------------------- toolbar
void App::viewerToolbar(ViewerBinding& b, ViewerUi& u) {
    Viewport3D& V = *b.view;
    ViewSettings& vs = V.settings();
    const float fs = ImGui::GetFontSize();
    ImGui::PushID(b.id.c_str());
    const int nx = b.grid ? b.grid->nx : 1, ny = b.grid ? b.grid->ny : 1, nz = b.grid ? b.grid->nz : 1;

    // colour map + range
    static const char* cmaps[] = {"Turbo", "Viridis", "Jet", "Cool-warm", "Grey"};
    flow(fs * 7.5f);
    int cm = int(vs.cmap);
    ImGui::SetNextItemWidth(fs * 7.5f);
    if (ImGui::Combo("##cmap", &cm, cmaps, 5)) vs.cmap = Colormap(cm);
    tip("Colour map");
    flow(btnW("Auto range") + ImGui::GetFrameHeight());
    ImGui::Checkbox("Auto range", b.autoRange);
    if (!*b.autoRange) {
        flow(fs * 12);
        ImGui::SetNextItemWidth(fs * 12);
        ImGui::DragFloatRange2("##range", b.rmin, b.rmax, std::max(1e-4f, (*b.rmax - *b.rmin) * 0.005f), 0, 0, "%.4g", "%.4g");
    }
    // opacity directly in the toolbar (CMG-style transparency)
    flow(fs * 7);
    ImGui::SetNextItemWidth(fs * 7);
    ImGui::SliderFloat("##opacity", &vs.opacity, 0.05f, 1.0f, ICON_EYE " %.2f");
    tip("Cell opacity: lower it to see wells, perforations and streamlines inside the grid");

    // display
    flow(btnW(ICON_LAYERS " Display"));
    if (ImGui::Button(ICON_LAYERS " Display")) ImGui::OpenPopup("display");
    if (ImGui::BeginPopup("display")) {
        ImGui::SeparatorText("Geometry");
        ImGui::SetNextItemWidth(fs * 10);
        ImGui::SliderFloat("Vertical exaggeration", &vs.zScale, 1.0f, 50.0f, "%.1fx", ImGuiSliderFlags_Logarithmic);
        ImGui::Checkbox("Cell edges", &vs.edges); ImGui::SameLine(); ImGui::Checkbox("Bounding box", &vs.bbox);
        ImGui::SameLine(); ImGui::Checkbox("Axes", &vs.axes);
        ImGui::SetNextItemWidth(fs * 10);
        ImGui::SliderFloat("Cell opacity", &vs.opacity, 0.05f, 1.0f, "%.2f");
        ImGui::SeparatorText("Colours");
        ImGui::Checkbox("Legend", &vs.legend); ImGui::SameLine(); ImGui::Checkbox("Logarithmic scale", &vs.logScale);
        ImGui::SetNextItemWidth(fs * 10);
        ImGui::SliderInt("Discrete colour bands", &vs.bands, 0, 24, vs.bands < 2 ? "continuous" : "%d");
        ImGui::SeparatorText("Wells");
        ImGui::Checkbox("Show wells", &vs.wells); ImGui::SameLine(); ImGui::Checkbox("Labels", &vs.labels);
        ImGui::SameLine(); ImGui::Checkbox("Draw on top", &vs.wellsXray);
        ImGui::SetNextItemWidth(fs * 10);
        ImGui::SliderFloat("Well radius", &vs.wellRadius, 0.3f, 4.0f, "%.1fx", ImGuiSliderFlags_Logarithmic);
        ImGui::TextDisabled("Wells follow their completion path (deviated / horizontal);\nperforations: green producer, blue injector, grey shut.");
        ImGui::EndPopup();
    }

    // dual porosity: matrix and fracture cells share the same geometry -> show one medium at a time (matrix first)
    const int nzm = b.grid ? b.grid->nzMatrix : 0;
    if (nzm > 0 && u.dualGrid != b.grid) {
        u.dualGrid = b.grid;
        if (vs.kMin <= 0 && vs.kMax >= nz - 1) { vs.kMin = 0; vs.kMax = nzm - 1; }
    } else if (nzm == 0) u.dualGrid = nullptr;
    if (nzm > 0) {
        flow(ImGui::GetFontSize() * 6.5f);
        int med = (vs.kMin >= nzm) ? 1 : (vs.kMax < nzm ? 0 : 2);
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6.5f);
        static const char* meds[] = {"Matrix", "Fracture", "Both"};
        if (ImGui::Combo("##medium", &med, meds, 3)) {
            if (med == 0) { vs.kMin = 0; vs.kMax = nzm - 1; }
            else if (med == 1) { vs.kMin = nzm; vs.kMax = nz - 1; }
            else { vs.kMin = 0; vs.kMax = nz - 1; }
            u.sliceMode = 0;
        }
        tip("Dual porosity: matrix cells are layers 1..NZ/2, fracture cells NZ/2+1..NZ of the same blocks");
    }

    // filters and slices
    flow(btnW(ICON_FILTER " Filter"));
    const bool filtered = vs.iMin > 0 || vs.jMin > 0 || vs.kMin > 0 || vs.iMax < nx - 1 || vs.jMax < ny - 1 || vs.kMax < nz - 1 || vs.filterValue;
    if (filtered) { ImVec4 ac = col4(Pal::Accent); ac.w = 0.55f; ImGui::PushStyleColor(ImGuiCol_Button, ac); }
    if (ImGui::Button(ICON_FILTER " Filter")) ImGui::OpenPopup("filter");
    if (filtered) ImGui::PopStyleColor();
    tip("IJK ranges, single I / J / K slices and value filter");
    if (ImGui::BeginPopup("filter")) {
        ImGui::SeparatorText("Single slice");
        static const char* sm[] = {"None", "I slice", "J slice", "K layer"};
        ImGui::SetNextItemWidth(fs * 7);
        bool sch = ImGui::Combo("##slicemode", &u.sliceMode, sm, 4);
        if (u.sliceMode > 0) {
            const int n = u.sliceMode == 1 ? nx : u.sliceMode == 2 ? ny : nz;
            ImGui::SameLine();
            ImGui::SetNextItemWidth(fs * 10);
            int s1 = std::clamp(u.slice, 0, n - 1) + 1;
            if (ImGui::SliderInt("##slice", &s1, 1, n)) sch = true;
            u.slice = s1 - 1;
            ImGui::SameLine();
            if (ImGui::ArrowButton("##sdn", ImGuiDir_Left)) { u.slice = std::max(0, u.slice - 1); sch = true; }
            ImGui::SameLine();
            if (ImGui::ArrowButton("##sup", ImGuiDir_Right)) { u.slice = std::min(n - 1, u.slice + 1); sch = true; }
        }
        if (sch) {
            vs.iMin = vs.jMin = vs.kMin = 0; vs.iMax = vs.jMax = vs.kMax = 1 << 30;
            if (u.sliceMode == 1) vs.iMin = vs.iMax = u.slice;
            if (u.sliceMode == 2) vs.jMin = vs.jMax = u.slice;
            if (u.sliceMode == 3) vs.kMin = vs.kMax = u.slice;
        }
        ImGui::SeparatorText("Index ranges (1-based, like ECLIPSE)");
        int i0 = std::max(0, vs.iMin) + 1, i1 = std::min(nx - 1, vs.iMax) + 1;
        int j0 = std::max(0, vs.jMin) + 1, j1 = std::min(ny - 1, vs.jMax) + 1;
        int k0 = std::max(0, vs.kMin) + 1, k1 = std::min(nz - 1, vs.kMax) + 1;
        bool ch = false;
        ImGui::SetNextItemWidth(fs * 12); ch |= ImGui::DragIntRange2("I range", &i0, &i1, 0.2f, 1, nx);
        ImGui::SetNextItemWidth(fs * 12); ch |= ImGui::DragIntRange2("J range", &j0, &j1, 0.2f, 1, ny);
        ImGui::SetNextItemWidth(fs * 12); ch |= ImGui::DragIntRange2("K range", &k0, &k1, 0.1f, 1, nz);
        if (ch) { vs.iMin = i0 - 1; vs.iMax = i1 - 1; vs.jMin = j0 - 1; vs.jMax = j1 - 1; vs.kMin = k0 - 1; vs.kMax = k1 - 1; u.sliceMode = 0; }
        ImGui::SeparatorText("Value filter");
        ImGui::Checkbox("Only cells with values in", &vs.filterValue);
        if (vs.filterValue) {
            ImGui::SetNextItemWidth(fs * 12);
            ImGui::DragFloatRange2("##vf", &vs.filterMin, &vs.filterMax, std::max(1e-4f, (*b.rmax - *b.rmin) * 0.005f), 0, 0, "%.4g", "%.4g");
            if (ImGui::SmallButton("Use colour range")) { vs.filterMin = *b.rmin; vs.filterMax = *b.rmax; }
        }
        ImGui::Separator();
        if (ImGui::Button(ICON_ROTATE_CCW " Reset filters")) { vs.iMin = vs.jMin = vs.kMin = 0; vs.iMax = vs.jMax = vs.kMax = 1 << 30; vs.filterValue = false; u.sliceMode = 0; }
        ImGui::EndPopup();
    }

    // well focus
    const auto& W = V.wells();
    if (!W.empty()) {
        flow(fs * 9);
        ImGui::SetNextItemWidth(fs * 9);
        const std::string cur = vs.focusWell >= 0 && vs.focusWell < int(W.size()) ? W[size_t(vs.focusWell)].name : std::string("All wells");
        if (ImGui::BeginCombo("##focus", (std::string(ICON_LOCATE_FIXED " ") + cur).c_str(), ImGuiComboFlags_HeightLarge)) {
            if (ImGui::Selectable("All wells (no focus)", vs.focusWell < 0)) vs.focusWell = -1;
            for (int i = 0; i < int(W.size()); ++i) {
                ImGui::PushID(i);
                ImGui::TextColored(col4(W[size_t(i)].injector ? IM_COL32(80, 160, 250, 255) : IM_COL32(60, 200, 110, 255)), W[size_t(i)].injector ? ICON_DROPLETS : ICON_CIRCLE_DOT);
                ImGui::SameLine();
                if (ImGui::Selectable(W[size_t(i)].name.c_str(), vs.focusWell == i)) { vs.focusWell = i; V.focusWellCamera(i); }
                ImGui::PopID();
            }
            ImGui::EndCombo();
        }
        tip("Focus a well: the camera centres on it, other wells are dimmed and (optionally) only its streamlines are shown.\nTip: click a well in the view.");
    }

    // streamlines
    if (b.sl) {
        SLState& sl = *b.sl;
        flow(btnW(ICON_ROUTE " Streamlines") + ImGui::GetFrameHeight() * 2);
        ImGui::BeginDisabled(!b.slAvailable && !b.slPrepare);
        if (ImGui::Checkbox(ICON_ROUTE " Streamlines", &sl.show) && sl.show && !b.slAvailable && b.slPrepare) b.slPrepare();
        ImGui::EndDisabled();
        tip(b.slAvailable ? "Pollock streamlines from injectors (and aquifers) to producers at the selected step" : b.slHint.c_str());
        ImGui::SameLine();
        if (ImGui::Button(ICON_SETTINGS "##sl")) ImGui::OpenPopup("slpop");
        tip("Streamline options and well allocation");
        if (ImGui::BeginPopup("slpop")) {
            if (!b.slAvailable) {
                ImGui::TextColored(col4(Pal::Warn), "%s", b.slHint.empty() ? "Streamlines are not available for these data." : b.slHint.c_str());
                if (b.slPrepare && ImGui::Button("Prepare streamlines")) b.slPrepare();
                ImGui::Separator();
            }
            ImGui::SetNextItemWidth(fs * 10);
            ImGui::SliderInt("Number of lines", &sl.lines, 20, 3000, "%d", ImGuiSliderFlags_Logarithmic);
            static const char* cb[] = {"Time of flight", "Injector (source)", "Producer (sink)"};
            ImGui::SetNextItemWidth(fs * 10);
            ImGui::Combo("Colour by", &sl.colorBy, cb, 3);
            ImGui::Checkbox("Hide cells", &sl.hideCells); ImGui::SameLine(); ImGui::Checkbox("X-ray (draw through cells)", &sl.xray);
            ImGui::Checkbox("With a focused well: only its streamlines", &sl.focusOnly);
            if (sl.res) {
                const StreamlineResult& S = *sl.res;
                ImGui::Text("%d lines (%d stalled), max TOF %.0f days, traced in %.1f ms", int(S.lines.size()), S.stalled, S.maxTof, sl.ms);
                if (!sl.source.empty()) ImGui::TextDisabled("%s", sl.source.c_str());
                if (!S.sources.empty() && !S.sinks.empty()) {
                    ImGui::SeparatorText("Well allocation (% of each injector's flow)");
                    const int ns = int(S.sinks.size());
                    if (ImGui::BeginTable("alloc", ns + 2, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_ScrollX, ImVec2(std::min(fs * 60, fs * 6 * (ns + 2)), 0))) {
                        ImGui::TableSetupColumn("from \\ to");
                        for (auto& n : S.sinks) ImGui::TableSetupColumn(n.c_str());
                        ImGui::TableSetupColumn("not reaching");
                        ImGui::TableHeadersRow();
                        for (size_t i = 0; i < S.sources.size(); ++i) {
                            double tot = 0; for (double x : S.pairFlux[i]) tot += x;
                            ImGui::TableNextRow();
                            ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(S.sources[i].c_str());
                            for (int k = 0; k <= ns; ++k) {
                                ImGui::TableSetColumnIndex(k + 1);
                                const double f = tot > 0 ? 100.0 * S.pairFlux[i][size_t(k)] / tot : 0.0;
                                if (f > 0.05) ImGui::TextColored(col4(f > 50 ? Pal::Ok : f > 10 ? Pal::Warn : Pal::Muted), "%.1f%%", f);
                                else ImGui::TextDisabled("-");
                            }
                        }
                        ImGui::EndTable();
                    }
                }
            }
            ImGui::TextDisabled("Cell properties TOF / DRAINAGE / SWEEP: Model viewer, Dynamic group.");
            ImGui::EndPopup();
        }
    }

    // camera presets + snapshot
    flow(ImGui::GetFrameHeight() * 6 + 20);
    if (ImGui::Button(ICON_SQUARE "##top")) V.viewPreset(0);
    tip("Top view (map)");
    ImGui::SameLine(0, 2);
    if (ImGui::Button(ICON_PANEL_BOTTOM "##front")) V.viewPreset(1);
    tip("Front view (looking north)");
    ImGui::SameLine(0, 2);
    if (ImGui::Button(ICON_PANEL_RIGHT "##side")) V.viewPreset(2);
    tip("Side view (looking west)");
    ImGui::SameLine(0, 2);
    if (ImGui::Button(ICON_BOX "##iso")) V.viewPreset(3);
    tip("3D view");
    ImGui::SameLine(0, 2);
    if (ImGui::Button(ICON_MAXIMIZE "##fit")) { V.resetCamera(); vs.focusWell = -1; }
    tip("Fit the whole model (clears the well focus)");
    ImGui::SameLine(0, 2);
    if (ImGui::Button(ICON_CAMERA "##png")) ImGui::OpenPopup("snap");
    tip("Save the view as a PNG image");
    if (ImGui::BeginPopup("snap")) {
        static int scale = 2;
        ImGui::TextUnformatted("Save view as PNG");
        ImGui::RadioButton("1x (screen size)", &scale, 1); ImGui::SameLine();
        ImGui::RadioButton("2x", &scale, 2); ImGui::SameLine();
        ImGui::RadioButton("3x (print)", &scale, 3);
        ImGui::TextDisabled("%d x %d pixels", V.width() * scale, V.height() * scale);
        if (ImGui::Button(ICON_DOWNLOAD " Save...")) {
            std::string path = safeFileName(b.pngName.empty() ? std::string("view") : b.pngName) + ".png";
            if (fileDialog(true, "PNG image (*.png)\0*.png\0", "png", path)) {
                std::string err;
                if (V.saveImage(path, scale, &err)) notify("View saved: " + path);
                else notify("Could not save the image: " + err);
            }
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }
    ImGui::PopID();
}

// ----------------------------------------------------------------------------- canvas
void App::viewerCanvas(ViewerBinding& b, ViewerUi& u, const std::function<void(int cell)>& tooltip) {
    Viewport3D& V = *b.view;
    ViewSettings& vs = V.settings();
    ImGui::PushID(b.id.c_str());
    const ImVec2 avail = ImGui::GetContentRegionAvail();
    const ImVec2 p0 = ImGui::GetCursorScreenPos();
    const int w = std::max(32, int(avail.x)), h = std::max(32, int(avail.y));
    ID3D11ShaderResourceView* srv = V.render(w, h);
    ImGui::InvisibleButton("##viewport", ImVec2(float(w), float(h)),
                           ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight | ImGuiButtonFlags_MouseButtonMiddle);
    ImGui::GetWindowDrawList()->AddImage((ImTextureID)(intptr_t)srv, p0, ImVec2(p0.x + w, p0.y + h));
    ImGuiIO& io = ImGui::GetIO();
    const float mx = io.MousePos.x - p0.x, my = io.MousePos.y - p0.y;
    if (ImGui::IsItemActive()) {
        if (ImGui::IsMouseDragging(ImGuiMouseButton_Left, 2.0f) && !io.KeyShift) V.orbit(io.MouseDelta.x, io.MouseDelta.y);
        if (ImGui::IsMouseDragging(ImGuiMouseButton_Right, 2.0f) || ImGui::IsMouseDragging(ImGuiMouseButton_Middle, 0.0f) ||
            (ImGui::IsMouseDragging(ImGuiMouseButton_Left, 0.0f) && io.KeyShift))
            V.pan(io.MouseDelta.x, io.MouseDelta.y);
    }
    const Grid* g = b.grid;
    if (ImGui::IsItemHovered()) {
        if (io.MouseWheel != 0) V.zoom(io.MouseWheel);
        if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) { V.resetCamera(); vs.focusWell = -1; }
        const int well = V.pickWell(mx, my);
        const bool click = ImGui::IsMouseReleased(ImGuiMouseButton_Left) && io.MouseDragMaxDistanceSqr[0] < 16.0f;
        if (click && well >= 0) {
            if (vs.focusWell == well) vs.focusWell = -1;
            else { vs.focusWell = well; V.focusWellCamera(well); }
        }
        const int cell = well >= 0 ? -1 : V.pick(mx, my);
        if (ImGui::IsMouseReleased(ImGuiMouseButton_Right) && io.MouseDragMaxDistanceSqr[1] < 16.0f) {
            u.ctxCell = cell; u.ctxWell = well;
            ImGui::OpenPopup("##vctx");
        }
        if (!ImGui::IsItemActive()) {
            if (well >= 0 && well < int(V.wells().size())) {
                const WellGlyph& wg = V.wells()[size_t(well)];
                ImGui::BeginTooltip();
                ImGui::TextColored(col4(wg.injector ? IM_COL32(90, 165, 250, 255) : IM_COL32(70, 205, 120, 255)), "%s %s", wg.injector ? ICON_DROPLETS : ICON_CIRCLE_DOT, wg.name.c_str());
                ImGui::TextDisabled("%s, %zu connections%s", wg.injector ? "injector" : "producer", wg.perfs.size(), wg.open ? "" : " (shut)");
                ImGui::TextDisabled(vs.focusWell == well ? "click: clear focus" : "click: focus this well");
                ImGui::EndTooltip();
            } else if (cell >= 0 && g && cell < g->numActive() && tooltip) {
                ImGui::BeginTooltip();
                tooltip(cell);
                if (b.cellHistory) ImGui::TextDisabled("right-click: value history, slices");
                ImGui::EndTooltip();
            }
        }
    }
    if (ImGui::BeginPopup("##vctx")) {
        if (u.ctxCell >= 0 && g && u.ctxCell < g->numActive()) {
            int i, j, k; g->ijk(g->globalOf[size_t(u.ctxCell)], i, j, k);
            ImGui::TextDisabled("Cell (%d, %d, %d)   %s = %.5g", i + 1, j + 1, k + 1, b.propName.c_str(), V.valueOf(u.ctxCell));
            if (b.cellHistory && ImGui::MenuItem(ICON_CHART_LINE "  Plot value history of this cell")) {
                u.histT.clear(); u.histV.clear();
                if (b.cellHistory(u.ctxCell, u.histT, u.histV)) {
                    u.histCell = u.ctxCell; u.histOpen = true;
                    char t[96]; std::snprintf(t, sizeof(t), "%s in cell (%d, %d, %d)", b.propName.c_str(), i + 1, j + 1, k + 1);
                    u.histProp = t;
                } else notify("No time series available for this property");
            }
            if (ImGui::BeginMenu(ICON_SCAN_LINE "  Slice through this cell")) {
                if (ImGui::MenuItem("I slice")) { vs.iMin = vs.iMax = i; vs.jMin = vs.kMin = 0; vs.jMax = vs.kMax = 1 << 30; u.sliceMode = 1; u.slice = i; }
                if (ImGui::MenuItem("J slice")) { vs.jMin = vs.jMax = j; vs.iMin = vs.kMin = 0; vs.iMax = vs.kMax = 1 << 30; u.sliceMode = 2; u.slice = j; }
                if (ImGui::MenuItem("K layer")) { vs.kMin = vs.kMax = k; vs.iMin = vs.jMin = 0; vs.iMax = vs.jMax = 1 << 30; u.sliceMode = 3; u.slice = k; }
                ImGui::EndMenu();
            }
            if (ImGui::MenuItem(ICON_COPY "  Copy cell info")) {
                char t[160]; std::snprintf(t, sizeof(t), "%d\t%d\t%d\t%s\t%.8g", i + 1, j + 1, k + 1, b.propName.c_str(), V.valueOf(u.ctxCell));
                ImGui::SetClipboardText(t);
            }
            ImGui::Separator();
        }
        if (u.ctxWell >= 0 && u.ctxWell < int(V.wells().size())) {
            const std::string n = V.wells()[size_t(u.ctxWell)].name;
            if (ImGui::MenuItem((std::string(ICON_LOCATE_FIXED "  Focus ") + n).c_str())) { vs.focusWell = u.ctxWell; V.focusWellCamera(u.ctxWell); }
        }
        if (vs.focusWell >= 0 && ImGui::MenuItem(ICON_X "  Clear well focus")) vs.focusWell = -1;
        if (ImGui::MenuItem(ICON_FILTER "  Reset filters / slices")) { vs.iMin = vs.jMin = vs.kMin = 0; vs.iMax = vs.jMax = vs.kMax = 1 << 30; vs.filterValue = false; u.sliceMode = 0; }
        if (ImGui::MenuItem(ICON_MAXIMIZE "  Reset camera")) V.resetCamera();
        ImGui::EndPopup();
    }
    ImGui::GetWindowDrawList()->AddText(ImVec2(p0.x + 10, p0.y + h - ImGui::GetFontSize() - 8), IM_COL32(150, 156, 170, 150),
                                        "Left drag: rotate   Right/Shift drag: pan   Wheel: zoom   Click a well: focus   Right-click: menu");
    // cell history window
    if (u.histOpen) {
        ImGui::SetNextWindowSize(ImVec2(560, 340), ImGuiCond_FirstUseEver);
        const std::string title = std::string(ICON_CHART_LINE " Cell history - ") + b.id + "###hist" + b.id;
        if (ImGui::Begin(title.c_str(), &u.histOpen, ImGuiWindowFlags_NoDocking)) {
            ImGui::TextUnformatted(u.histProp.c_str());
            if (ImPlot::BeginPlot("##cellhist", ImVec2(-1, -1))) {
                ImPlot::SetupAxes("Time (days)", b.propName.c_str(), ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                std::vector<double> t, v;
                for (size_t q = 0; q < u.histT.size() && q < u.histV.size(); ++q) if (std::isfinite(u.histV[q])) { t.push_back(u.histT[q]); v.push_back(u.histV[q]); }
                if (!t.empty()) {
                    ImPlotSpec sp; sp.LineColor = col4(Pal::Accent); sp.LineWeight = 2.0f; sp.Marker = ImPlotMarker_Circle; sp.MarkerSize = 3;
                    ImPlot::PlotLine(b.propName.c_str(), t.data(), v.data(), int(t.size()), sp);
                }
                ImPlot::EndPlot();
            }
        }
        ImGui::End();
    }
    ImGui::PopID();
}
