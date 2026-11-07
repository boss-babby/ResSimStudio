// C42 geomodel data import (Petrel-style "import with column assignment"):
//   * generic delimited text tables - delimiter / header / comment / quote / decimal-comma sniffing, Petrel
//     "BEGIN HEADER ... END HEADER" blocks, "# key: value" metadata lines, skip rows, preview of the first N rows
//   * a column map (role + unit factor per column) with synonym-based auto-assignment that the user can override
//   * typed builders for well heads, deviation surveys (one or many wells per file; MD/INC/AZI, X/Y/TVD(SS)/Z or
//     DX/DY offsets), well tops (MD or TVDSS), scattered points (one or many sets) and fault sticks
//   * well matching / merging (case-insensitive, file-name prefix) for batch imports of many files
//   * 2-D SEG-Y lines (configurable trace-header byte locations, formats 1/2/3/5/8, either byte order), 1-D horizon
//     tracking along a line, line intersections (mis-tie analysis) and conversion of line picks to depth points
#pragma once
#include "ressim/geomodel.h"
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace ressim {
namespace geo {

// ---------------------------------------------------------------- delimited text tables
struct TableFormat {
    char delim = 0;                      // 0 = auto, ',' ';' '\t' '|' or ' ' (runs of blanks)
    int skipLines = 0;                   // raw lines dropped before anything else
    int headerRow = -2;                  // -2 auto, -1 no header, k = k-th non-comment line (after skipLines) is the header
    std::string comment = "#!";          // first-character comment leaders ("--" also accepted when '-' is listed)
    bool quotes = true;                  // "..." / '...' fields keep blanks and delimiters
    int decimalComma = -1;               // -1 auto, 0 '.', 1 ',' as decimal separator (only with a non-comma delimiter)
    std::vector<std::string> nulls = {"-999.25", "-999", "-9999", "9999900", "NaN", "nan", "NA", "N/A", ""};
    int maxRows = 0;                     // 0 = all rows (preview: e.g. 200)
};
struct TextTable {
    std::vector<std::string> header;     // one name per column (header row, Petrel header block or C1..Cn)
    std::vector<std::vector<std::string>> rows;   // data cells (missing cells = "")
    std::map<std::string, std::string> meta;      // "# WELL NAME: A-1" style lines (upper-case key)
    std::vector<std::string> rawHead;    // first raw lines of the file (preview)
    int ncol = 0;
    long long totalRows = 0;             // data rows in the file (also when maxRows truncated 'rows')
    bool hasHeader = false;
    char delimUsed = ' ';
    bool decimalCommaUsed = false;
};
TextTable parseTextTable(const std::string& text, const TableFormat& fmt, std::string& err);
TextTable readTextTable(const std::string& path, const TableFormat& fmt, std::string& err);
// numeric cell (honours decimal comma / null values); false = text or null
bool tableNumber(const TextTable& t, const TableFormat& fmt, const std::string& cell, double& v);

// ---------------------------------------------------------------- column roles / map
enum class ColRole {
    Ignore, WellName, X, Y, KB, TD, MD, Inc, Azi, TVD, TVDSS, Z, DX, DY, Surface, StickId, FaultName, Value, Count
};
const char* colRoleName(ColRole r);
enum class ImportKind { WellHeads, Deviation, Tops, Points, FaultSticks, Count };
const char* importKindName(ImportKind k);
std::vector<ColRole> rolesFor(ImportKind k);          // assignable roles (Ignore first)
struct ColumnMap {
    std::vector<ColRole> role;           // per column
    std::vector<double> scale;           // per column factor applied to numeric values (unit conversion)
    int find(ColRole r) const { for (size_t i = 0; i < role.size(); ++i) if (role[i] == r) return int(i); return -1; }
};
std::string normColName(const std::string& s);           // upper case, letters and digits only
ColumnMap guessColumns(const TextTable& t, const TableFormat& fmt, ImportKind k);
// "" when the map can build the kind, else the reason ("needs X and Y columns", ...)
std::string validateMap(const ColumnMap& m, ImportKind k);

// ---------------------------------------------------------------- typed builders
struct ImportOptions {
    std::string defaultWell;             // deviation / tops / LAS without a well column (e.g. from the file name)
    std::string defaultSet;              // point set / fault name without a name column
    int zSign = 0;                       // Z column of points / sticks / tops: 0 auto (mostly negative = elevation), 1 elevation (depth = -Z), 2 depth
};
struct TopRow { std::string well, surface; double md = 0, tvdss = 0; bool hasMd = false, hasTvd = false, tvdBelowKb = false; };
struct ImportResult {
    std::vector<GeoWell> wells;          // heads (x, y, kb) or survey wells (traj; see trajMode)
    int trajMode = 0;                    // 0 MD/INC/AZI (minimum curvature), 1 absolute X/Y/TVDSS, 2 DX/DY + TVDSS offsets
    bool tvdBelowKb = false;             // survey depth column was TVD below the KB (converted with the project KB on merge)
    std::vector<double> td;              // per well (heads with a TD column), NaN = none
    std::vector<TopRow> tops;
    std::vector<std::pair<std::string, std::vector<P3>>> points;   // named point sets (depth positive down)
    std::vector<Fault> faults;
    int rowsUsed = 0, rowsSkipped = 0;
    std::vector<std::string> warnings;
};
ImportResult buildImport(const TextTable& t, const TableFormat& fmt, const ColumnMap& m, ImportKind k, const ImportOptions& o);

// ---------------------------------------------------------------- well matching / merging
std::string wellKey(const std::string& name);            // upper case, blanks / '_' / '-' kept, trimmed
// project well with the same key, else the longest project well name that prefixes 'nameOrFile' at a separator; -1 none
int matchWell(const std::vector<GeoWell>& wells, const std::string& nameOrFile);
struct MergeStats { int added = 0, updated = 0; std::vector<std::string> warnings; };
MergeStats mergeHeads(std::vector<GeoWell>& proj, const ImportResult& r);
// survey wells into the project: creates missing wells (head from absolute surveys, else 0,0 with a warning),
// recomputes the trajectory from the project head / KB
MergeStats mergeSurveys(std::vector<GeoWell>& proj, const ImportResult& r);
MergeStats mergeTops(std::vector<GeoWell>& proj, const std::vector<TopRow>& tops, bool createMissing);
// trajectory editing helpers (well manager)
void setVertical(GeoWell& w, double td);                 // two stations, MD 0 .. td, vertical
bool trajFromXyz(GeoWell& w, std::vector<Station> st, std::string& err);   // st: md (or <0 = from arc length), x, y, tvd(ss)
void inclinationsFromXyz(GeoWell& w);                    // fills inc / azi from x / y / tvd (display / export)
std::string writeDeviationTable(const GeoWell& w);       // "MD INC AZI TVDSS X Y" text (clipboard / file)

// ---------------------------------------------------------------- 2-D seismic lines
struct SegyByteMap {                     // 1-based trace-header byte positions (SEG-Y rev 1 defaults)
    int inlineByte = 189, xlineByte = 193, cdpByte = 21, xByte = 181, yByte = 185, scalarByte = 71;
    bool fallbackXY = true;              // X/Y zero -> source X/Y (73 / 77)
};
struct SeismicLine {
    std::string name, path, textHeader;
    int ns = 0, format = 5;
    double t0 = 0, dt = 4;
    std::vector<int> cdp;
    std::vector<double> x, y;            // plan position per trace
    std::vector<float> a;                // [trace][sample]
    int ntr() const { return int(x.size()); }
    float at(int tr, int s) const { return a[size_t(tr) * size_t(ns) + size_t(s)]; }
    double length() const;
    std::vector<double> distance() const;   // cumulative distance per trace
};
// 2 = 2-D line, 3 = 3-D volume, 0 = unreadable (err); 'info' describes the trace geometry found
int segyGeometryKind(const std::string& path, const SegyByteMap& m, std::string& info, std::string& err);
SeismicLine readSegy2D(const std::string& path, const SegyByteMap& m, std::string& err);
// 3-D post-stack volume with user byte locations (the two-argument readSegy uses the defaults)
SeismicCube readSegy(const std::string& path, const SegyByteMap& m, std::string& err);
bool writeSegy2D(const std::string& path, const SeismicLine& l);   // IEEE float, CDP + X/Y headers (tests, export)
struct LineHorizon {                     // picks of one horizon on every loaded line
    std::string name;
    std::vector<std::vector<double>> t;  // [line][trace] TWT ms, NaN = not picked
};
// 1-D seeded tracking along a line (mode 0 peak, 1 trough, 2 zero crossing + to -): from every seed outwards, the
// extremum within +-window ms of the previous pick; stops where the jump exceeds maxJump ms or the event fades
std::vector<double> trackLine(const SeismicLine& l, const std::vector<std::pair<int, double>>& seeds, int mode, double window, double maxJump);
// manual picks: linear interpolation of the TWT between clicked traces (outside the first / last: NaN)
std::vector<double> interpolatePicks(int ntr, std::vector<std::pair<int, double>> pts);
struct LineTie { int a = -1, b = -1, ta = 0, tb = 0; double x = 0, y = 0; };   // trace ta of line a crosses trace tb of b
std::vector<LineTie> lineIntersections(const std::vector<SeismicLine>& lines);
// depth points of a line horizon (every 'step'-th trace) through the velocity model
std::vector<P3> lineHorizonPoints(const std::vector<SeismicLine>& lines, const LineHorizon& h, const VelocityModel& v, int step = 1);
// bounding lattice for gridding scattered data when no horizon lattice exists (n x n nodes, padded 5 %)
Surface latticeFor(const std::vector<P3>& pts, int n);

}  // namespace geo
}  // namespace ressim
