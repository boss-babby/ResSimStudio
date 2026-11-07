#pragma once
// C46 geomodel project file (.rgeo): the interpretation and input data of the Geomodel workspace / ResGeo in one
// versioned text file - surfaces, wells (head, survey, logs, tops), polygons, faults (sticks), 2-D seismic lines
// (geometry + file reference; amplitudes are re-read from the SEG-Y file on open), line horizons, the 3-D SEG-Y
// reference, user folders of the object tree and free key / value settings of the application (templates, views).
// Grids and property models are derived data: they are rebuilt from these inputs (Grid step), not stored.
#include <string>
#include <utility>
#include <vector>
#include "ressim/geomodel.h"
#include "ressim/geo_import.h"
#include "ressim/geo_tools.h"

namespace ressim {
namespace geo {

// object kinds of the tree folders
enum class ObjKind : int { Well = 0, Surface = 1, Polygon = 2, Fault = 3, Line = 4 };
constexpr int kObjKinds = 5;
const char* objKindName(ObjKind k);   // "well", "surface", ...

struct ObjFolders {
    std::vector<std::pair<int, std::string>> folders;                // (kind, folder name), in creation order
    std::vector<std::pair<int, std::pair<std::string, std::string>>> member;   // (kind, (object name, folder))
    // folder of an object ("" = at the root of its kind)
    std::string folderOf(ObjKind k, const std::string& obj) const;
    void assign(ObjKind k, const std::string& obj, const std::string& folder);   // "" removes the assignment
    bool addFolder(ObjKind k, const std::string& name);                         // false: exists / empty
    void renameFolder(ObjKind k, const std::string& from, const std::string& to);
    void removeFolder(ObjKind k, const std::string& name);                      // members go back to the root
    void renameObject(ObjKind k, const std::string& from, const std::string& to);
    void forget(ObjKind k, const std::string& obj);                             // object deleted
    std::vector<std::string> foldersOf(ObjKind k) const;
};

struct Project {
    bool field = false;
    std::vector<Surface> surfaces;
    std::vector<GeoWell> wells;
    std::vector<Polygon> polygons;
    std::vector<Fault> faults;
    std::vector<SeismicLine> lines;       // amplitudes 'a' are not written; empty after load unless re-read
    std::vector<LineHorizon> lineHorizons;
    std::string cubePath;                 // 3-D SEG-Y volume ("" = none)
    ObjFolders folders;
    std::vector<std::pair<std::string, std::string>> settings;   // application key / value pairs (repeatable keys)
};

bool saveProject(const std::string& path, const Project& p, std::string* err = nullptr);
bool loadProject(const std::string& path, Project& p, std::string& err);
// in-memory variants (tests, C API)
std::string projectToText(const Project& p);
bool projectFromText(const std::string& text, Project& p, std::string& err);

}  // namespace geo
}  // namespace ressim
