// C46 Geomodel project file (.rgeo): save / open of the workspace inputs and interpretation (engine geo_project.h),
// plus the application settings that belong to a project (templates, well-section layout, SEG-Y byte map, velocity)
#include "i18n.h"
#include "geo_state.h"
#include "ui_geoimport.h"
#include "widgets.h"
#include "icons.h"
#include "platform.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <sstream>

namespace {
std::vector<std::string> splitBar(const std::string& s) {
    std::vector<std::string> out; std::string cur;
    for (char c : s) { if (c == '|') { out.push_back(cur); cur.clear(); } else cur += c; }
    out.push_back(cur);
    return out;
}
std::string num(double v) { char b[40]; std::snprintf(b, sizeof b, "%.17g", v); return b; }
double toD(const std::string& s, double d = 0) { char* e = nullptr; const double v = std::strtod(s.c_str(), &e); return e == s.c_str() ? d : v; }
std::string joinNames(const std::vector<std::string>& v) { std::string o; for (const auto& s : v) { if (!o.empty()) o += '|'; o += s; } return o; }
}  // namespace

ressim::geo::Project geoToProject(const App::GeoState& S) {
    ressim::geo::Project P;
    P.field = S.field;
    P.surfaces = S.hz; P.wells = S.wells; P.polygons = S.polys; P.faults = S.faults;
    P.lines.reserve(S.lines.size());
    for (const auto& l : S.lines) {   // geometry only: amplitudes are re-read from the SEG-Y file
        SeismicLine g; g.name = l.name; g.path = l.path; g.ns = l.ns; g.format = l.format; g.t0 = l.t0; g.dt = l.dt; g.cdp = l.cdp; g.x = l.x; g.y = l.y;
        P.lines.push_back(std::move(g));
    }
    P.lineHorizons = S.lineHz;
    P.cubePath = S.hasCube ? S.cubePath : std::string();
    P.folders = S.folders;
    auto set = [&](const std::string& k, const std::string& v) { P.settings.push_back({k, v}); };
    for (const auto& t : S.templates)
        set("template", t.name + "|" + t.match + "|" + t.unit + "|" + std::to_string(t.cmap) + "|" + (t.reverse ? "1" : "0") + "|" + (t.fixed ? "1" : "0") + "|" +
                            (t.logScale ? "1" : "0") + "|" + num(t.lo) + "|" + num(t.hi));
    for (const auto& t : S.wsTracks) set("ws.track", t.curve + "|" + (t.logX ? "1" : "0") + "|" + num(t.width) + "|" + (t.useTpl ? "1" : "0"));
    set("ws.hidden", joinNames(S.wsHidden));
    set("ws.order", joinNames(S.wsOrder));
    set("ws.flatten", S.wsFlat);
    const SegyByteMap& m = S.segyMap;
    set("segy.bytes", std::to_string(m.inlineByte) + "|" + std::to_string(m.xlineByte) + "|" + std::to_string(m.cdpByte) + "|" + std::to_string(m.xByte) + "|" +
                          std::to_string(m.yByte) + "|" + std::to_string(m.scalarByte) + "|" + (m.fallbackXY ? "1" : "0"));
    std::string on; for (char c : S.lineOn) on += c ? '1' : '0';
    set("lines.on", on);
    set("velocity", num(S.vm.v0) + "|" + num(S.vm.k) + "|" + num(S.vm.datum));
    set("horizon.name", S.thName);
    return P;
}

// replaces the project inputs of S; derived models are cleared (rebuilt from the Grid step on)
void geoFromProject(App::GeoState& S, ressim::geo::Project&& P, const std::string& path) {
    S.field = P.field;
    S.hz = std::move(P.surfaces); S.wells = std::move(P.wells); S.polys = std::move(P.polygons); S.faults = std::move(P.faults);
    S.folders = std::move(P.folders);
    S.templates.clear(); S.wsTracks.clear(); S.wsHidden.clear(); S.wsOrder.clear(); S.wsFlat.clear();
    std::vector<char> lineOn;
    for (const auto& [k, v] : P.settings) {
        const auto f = splitBar(v);
        if (k == "template" && f.size() >= 9) {
            App::GeoState::PropTemplate t; t.name = f[0]; t.match = f[1]; t.unit = f[2]; t.cmap = std::atoi(f[3].c_str());
            t.reverse = f[4] == "1"; t.fixed = f[5] == "1"; t.logScale = f[6] == "1"; t.lo = toD(f[7]); t.hi = toD(f[8], 1);
            S.templates.push_back(t); S.templatesInit = true;
        } else if (k == "ws.track" && f.size() >= 4) {
            App::GeoState::WsTrack t; t.curve = f[0]; t.logX = f[1] == "1"; t.width = float(toD(f[2], 1)); t.useTpl = f[3] == "1";
            S.wsTracks.push_back(t);
        } else if (k == "ws.hidden") { for (const auto& s : f) if (!s.empty()) S.wsHidden.push_back(s); }
        else if (k == "ws.order") { for (const auto& s : f) if (!s.empty()) S.wsOrder.push_back(s); }
        else if (k == "ws.flatten") S.wsFlat = v;
        else if (k == "segy.bytes" && f.size() >= 7) {
            SegyByteMap& m = S.segyMap;
            m.inlineByte = std::atoi(f[0].c_str()); m.xlineByte = std::atoi(f[1].c_str()); m.cdpByte = std::atoi(f[2].c_str());
            m.xByte = std::atoi(f[3].c_str()); m.yByte = std::atoi(f[4].c_str()); m.scalarByte = std::atoi(f[5].c_str()); m.fallbackXY = f[6] == "1";
        } else if (k == "lines.on") { for (char c : v) lineOn.push_back(c == '1'); }
        else if (k == "velocity" && f.size() >= 3) { S.vm.v0 = toD(f[0], S.vm.v0); S.vm.k = toD(f[1], S.vm.k); S.vm.datum = toD(f[2], S.vm.datum); }
        else if (k == "horizon.name") S.thName = v;
    }
    // 2-D lines: amplitudes from the referenced SEG-Y files (geometry and picks are in the project either way)
    int reread = 0, missing = 0;
    S.lines = std::move(P.lines);
    for (auto& l : S.lines) {
        std::error_code ec;
        if (!l.path.empty() && std::filesystem::exists(l.path, ec)) {
            std::string e; SeismicLine r = readSegy2D(l.path, S.segyMap, e);
            if (r.ntr() == l.ntr() && r.ns > 0) { l.a = std::move(r.a); l.ns = r.ns; l.textHeader = r.textHeader; ++reread; continue; }
        }
        l.a.assign(size_t(std::max(0, l.ntr())) * size_t(std::max(0, l.ns)), 0.0f);   // geometry only (blank section)
        ++missing;
    }
    S.lineHz = std::move(P.lineHorizons);
    for (auto& h : S.lineHz) {
        h.t.resize(S.lines.size());
        for (size_t l = 0; l < S.lines.size(); ++l) h.t[l].resize(size_t(S.lines[l].ntr()), std::nan(""));
    }
    if (S.lineHz.empty() && !S.lines.empty()) { LineHorizon h; h.name = "H1"; h.t.resize(S.lines.size()); for (size_t l = 0; l < S.lines.size(); ++l) h.t[l].assign(size_t(S.lines[l].ntr()), std::nan("")); S.lineHz.push_back(h); }
    S.lineOn = lineOn; S.lineOn.resize(S.lines.size(), 1);
    S.curLine = S.lines.empty() ? -1 : 0; S.curLineHz = 0;
    // 3-D volume
    S.hasCube = false; S.segy.reset(); S.cube = SeismicCube(); S.cubePath.clear();
    std::string cubeMsg;
    if (!P.cubePath.empty()) {
        std::error_code ec;
        if (!std::filesystem::exists(P.cubePath, ec)) cubeMsg = "; 3-D volume not found: " + P.cubePath;
        else {
            const auto bytes = std::filesystem::file_size(P.cubePath, ec);
            std::string e;
            if (!ec && bytes > 400ull * 1048576ull) {
                auto v = std::make_shared<SegyVolume>();
                const SegyByteMap& bm = S.segyMap;
                v->byteIl = bm.inlineByte; v->byteXl = bm.xlineByte; v->byteX = bm.xByte; v->byteY = bm.yByte; v->byteScalar = bm.scalarByte;
                if (v->open(P.cubePath, e)) { S.segy = v; S.cube = v->proxy(160, 512); S.hasCube = true; S.il = v->ni / 2; S.xl = v->nj / 2; }
            } else {
                S.cube = readSegy(P.cubePath, S.segyMap, e);
                S.hasCube = e.empty() && S.cube.ni > 0;
            }
            if (S.hasCube) S.cubePath = P.cubePath; else cubeMsg = "; 3-D volume: " + e;
        }
    }
    // derived data and caches
    S.hasGrid = S.hasFacies = S.hasPetro = S.hasVol = S.hasUp = S.hasTh = false;
    S.blocked.clear(); S.seeds.clear(); S.stick.clear(); S.draft = Fault(); S.hzUndo.clear();
    S.mapSurface = 0; S.curPoly = -1; S.volPoly = -1; S.isecPath = -1;
    ++S.hzRev; ++S.lineRev; ++S.cubeRev;
    S.projectPath = path;
    S.dir = std::filesystem::path(path).parent_path().string();
    char b[512];
    std::snprintf(b, sizeof b, "project %s: %zu surfaces, %zu wells, %zu polygons, %zu faults, %zu 2-D lines (%d re-read%s)%s",
                  path.c_str(), S.hz.size(), S.wells.size(), S.polys.size(), S.faults.size(), S.lines.size(), reread,
                  missing ? (", " + std::to_string(missing) + " without the SEG-Y file").c_str() : "", cubeMsg.c_str());
    S.msg(b);
}

bool geoSaveProject(App::GeoState& S, bool saveAs) {
    std::string p = S.projectPath.empty() ? std::string("project.rgeo") : S.projectPath;
    if ((saveAs || S.projectPath.empty()) && !plat::fileDialog(true, "ResGeo project (*.rgeo)\0*.rgeo\0All files (*.*)\0*.*\0", "rgeo", p)) return false;
    std::string e;
    if (!ressim::geo::saveProject(p, geoToProject(S), &e)) { S.msg("Error: " + e); return false; }
    S.projectPath = p;
    S.msg("project saved: " + p + " (grids and property models are rebuilt from the inputs)");
    return true;
}
bool geoOpenProject(App::GeoState& S, const std::string& path) {
    ressim::geo::Project P; std::string e;
    if (!ressim::geo::loadProject(path, P, e)) { S.msg("Error: " + path + ": " + e); return false; }
    geoFromProject(S, std::move(P), path);
    return true;
}

// Data step: project file buttons
void geoProjectBar(App::GeoState& S) {
    if (ImGui::Button(TR(ICON_FOLDER_OPEN " Open project..."))) {
        std::string p;
        if (plat::fileDialog(false, "ResGeo project (*.rgeo)\0*.rgeo\0All files (*.*)\0*.*\0", "rgeo", p)) geoOpenProject(S, p);
    }
    ImGui::SameLine();
    if (ImGui::Button(TR(ICON_SAVE " Save project"))) geoSaveProject(S, false);
    ImGui::SameLine();
    if (ImGui::Button(TR("Save as..."))) geoSaveProject(S, true);
    if (!S.projectPath.empty()) { ImGui::SameLine(); ImGui::TextDisabled("%s", std::filesystem::path(S.projectPath).filename().string().c_str()); }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", TR("Surfaces, wells (survey, logs, tops), polygons, faults, 2-D line geometry and picks, the 3-D SEG-Y "
                                                       "reference, tree folders, templates and the well-section layout. Grids and property models are "
                                                       "rebuilt from these inputs."));
}
