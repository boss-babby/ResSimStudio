// Reader for ECLIPSE / OPM Flow result files (binary big-endian and formatted):
//   EGRID / FEGRID     corner-point grid (GRIDHEAD, COORD, ZCORN, ACTNUM)
//   INIT  / FINIT      static cell properties (PORV, PERMX, PORO, ... active-cell ordering)
//   UNRST / FUNRST     unified restart (SEQNUM + INTEHEAD/DOUBHEAD + PRESSURE, SWAT, SGAS, RS, ...)
//   X0001 / F0001 ...  non-unified restart
//   SMSPEC / FSMSPEC   summary specification;  UNSMRY / FUNSMRY or S0001 / A0001 ... summary data
// Files are indexed on open; array data is read on demand (large restart files stay on disk).
#pragma once
#include "ressim/case.h"
#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace ressim {
class Grid;
struct SummaryPoint;
struct Snapshot;

namespace ecl {

enum class Type : uint8_t { INTE, REAL, DOUB, LOGI, CHAR, MESS, Other };

struct ArrayInfo {
    std::string name;     // trimmed keyword (8 chars max)
    Type type = Type::Other;
    std::string typeName; // raw 4-char type ("INTE", "C016", ...)
    int64_t count = 0;
    int strLen = 8;       // CHAR / C0nn element length
    int64_t offset = 0;   // binary: file offset of the first data record; formatted: offset of the first data line
};

// One ECLIPSE file, indexed.
class File {
public:
    bool open(const std::string& path, std::string* err = nullptr);
    bool isOpen() const { return !path_.empty(); }
    const std::string& path() const { return path_; }
    bool formatted() const { return formatted_; }
    const std::vector<ArrayInfo>& arrays() const { return idx_; }
    // index of the n-th occurrence of an array (optionally only after array 'from'), -1 if absent
    int find(const std::string& name, int occurrence = 0, int from = 0) const;
    std::vector<double> numeric(int i) const;        // INTE/REAL/DOUB/LOGI as double
    std::vector<float> floats(int i) const;
    std::vector<int> ints(int i) const;
    std::vector<std::string> strings(int i) const;   // CHAR/C0nn, trimmed
    std::vector<double> numeric(const std::string& name, int occurrence = 0) const { int i = find(name, occurrence); return i < 0 ? std::vector<double>() : numeric(i); }
    std::vector<int> ints(const std::string& name, int occurrence = 0) const { int i = find(name, occurrence); return i < 0 ? std::vector<int>() : ints(i); }
private:
    template <class F> void readValues(int i, F&& put) const;   // put(index, double, const std::string*)
    std::string path_;
    bool formatted_ = false;
    std::vector<ArrayInfo> idx_;
    mutable std::ifstream in_;
};

// ---------------------------------------------------------------- grid
struct EGrid {
    int nx = 0, ny = 0, nz = 0;
    std::vector<double> coord, zcorn;
    std::vector<int> actnum;      // size nx*ny*nz (empty = all active)
    std::string mapUnits, gridUnits;
    int numActive() const;
};
bool readEGrid(const std::string& path, EGrid& g, std::string* err = nullptr);
// Grid usable by the 3D viewer: corner geometry + active map taken exactly from ACTNUM (no physics set-up).
std::shared_ptr<Grid> displayGrid(const EGrid& g);

// ---------------------------------------------------------------- summary
struct SummaryVar {
    std::string keyword, wgname, unit;
    int num = 0;
    std::string key;              // "FOPR", "WOPR:PROD", "BPR:10,10,3", "RPR:2", "CWIR:INJ:1,1,1"
    char category = 'M';          // F W G R B C A S T(time) N(numerics) M(misc)
};
struct Summary {
    std::vector<SummaryVar> vars;
    std::vector<double> time;                    // days since start, per ministep
    std::vector<int> reportIndex;                // ministep index of each report step (end of step)
    std::vector<std::vector<float>> values;      // [var][ministep]
    int startDay = 1, startMonth = 1, startYear = 2000;
    int nx = 0, ny = 0, nz = 0;
    std::string unitSystem;                      // FIELD / METRIC / LAB (from units of TIME/FPR when known)
    int find(const std::string& key) const;      // exact key, -1 if missing
    std::vector<double> series(int var) const;
    double at(int var, double t) const;          // linear interpolation in time (NaN if no data)
};
// smspec: path of the .SMSPEC/.FSMSPEC (data file(s) located automatically)
bool readSummary(const std::string& smspec, Summary& s, std::string* err = nullptr);
// Write ECLIPSE binary summary files BASE.SMSPEC + BASE.UNSMRY (big-endian Fortran records, the format read by
// ResInsight, OPM tools, resdata/ecl2df). TIME is written first; report steps follow s.reportIndex.
bool writeSummary(const Summary& s, const std::string& base, std::string* err = nullptr);

// ResSim run -> same vector naming (FOPR, FWCT, WBHP:NAME, ...)
Summary summaryFromRun(const std::vector<SummaryPoint>& pts, const std::vector<std::string>& wellNames, UnitSystem units);

// external summary -> ResSim run representation (one point per ministep; wells matched by name, missing = shut)
std::vector<SummaryPoint> runFromSummary(const Summary& s, const std::vector<std::string>& wellNames);

// ---------------------------------------------------------------- restart
struct RestartStep {
    int report = 0;               // SEQNUM / file number
    double days = 0;
    int day = 0, month = 0, year = 0;
    int file = 0;                 // index into Restart::files
    int first = 0, last = 0;      // array index range [first, last) in that file
};
class Restart {
public:
    bool open(const std::vector<std::string>& files, std::string* err = nullptr);   // one UNRST or many X/F files
    const std::vector<RestartStep>& steps() const { return steps_; }
    std::vector<std::string> arrayNames(int step) const;          // per-cell arrays (count == numActive)
    std::vector<float> load(int step, const std::string& name) const;
    void setNumActive(int na) { na_ = na; }
private:
    std::vector<std::unique_ptr<File>> files_;
    std::vector<RestartStep> steps_;
    int na_ = 0;
};

// ---------------------------------------------------------------- whole case
struct CaseFiles {
    std::string base;             // directory + case name without extension
    std::string egrid, init, smspec;
    std::vector<std::string> restart;   // UNRST or X0000... files
};
CaseFiles locateCase(const std::string& anyFile);   // from CASE.DATA / CASE.EGRID / CASE.SMSPEC / ...

struct ResultCase {
    std::string name, path;
    CaseFiles files;
    bool hasGrid = false, hasInit = false, hasSummary = false, hasRestart = false;
    EGrid egrid;
    std::shared_ptr<Grid> grid;           // display grid (may be null when no EGRID)
    File init;
    std::vector<std::string> initArrays;  // per active cell static arrays in INIT
    Summary summary;
    Restart restart;
    int unitSystem = 0;                   // 1 METRIC 2 FIELD 3 LAB (INTEHEAD[2]), 0 unknown
    std::vector<std::string> warnings;
    std::vector<float> initArray(const std::string& name) const;   // per active cell (PORV reduced to active)
};
bool openResultCase(const std::string& anyFile, ResultCase& rc, std::string* err = nullptr);
// restart steps as ResSim snapshots (PRESSURE, SWAT, SGAS, RS; SOIL = 1 - SWAT - SGAS), per active cell of rc.grid
std::vector<Snapshot> snapshotsFromRestart(const ResultCase& rc, int maxSteps = -1);

// ---------------------------------------------------------------- comparison
struct VectorDiff {
    std::string key;
    double refFinal = 0, cmpFinal = 0;
    double maxAbs = 0, rms = 0;           // over the reference report times
    double relRms = 0;                    // rms / max(|ref|) (0..)
    int samples = 0;
};
// compares the vectors present in both cases at the report times of 'ref' (or all its ministeps)
std::vector<VectorDiff> compareSummaries(const Summary& ref, const Summary& cmp, const std::vector<std::string>& keys = {});

} // namespace ecl
} // namespace ressim
