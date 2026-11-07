// Headless harness stand-in for app/src/app.h (C36): only what ui_geo3d.cpp / geo_state.h need.
#pragma once
#include "imgui.h"
#include <functional>
#include <memory>
#include <string>
bool fileDialog(bool save, const char* filter, const char* defExt, std::string& path);
class App {
public:
    struct GeoState;
    void geomodel3D(GeoState& S, ImVec2 size);
    void geoSeismicSection(GeoState& S, ImVec2 size);
    void notify(const std::string&) {}
};
