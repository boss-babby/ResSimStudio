// C14: RunResults <-> .rsres, ECLIPSE result export (see results_io.h)
#include "results_io.h"
#include "ressim/binio.h"
#include "ressim/ecl.h"
#include "ressim/grid.h"
#include "ressim/io.h"
#include <chrono>
#include <ctime>
#include <stdexcept>

namespace {
using namespace ressim;

// C13 geothermal / EDFM analysis: one BLOB "geothermal"
template <class Ar> void geoIO(Ar& ar, RunResults& r) {
    ar(r.fracNames, r.heatInPlace, r.tInjRef, r.cooledFrac, r.wellT, r.wellMW);
    uint64_t n = r.fracSegs.size(); ar(n);
    if (Ar::reading) r.fracSegs.resize(size_t(n));
    for (auto& g : r.fracSegs) ar(g.frac, g.host, g.xa, g.ya, g.xb, g.yb, g.zt, g.zb);
    n = r.fracSnap.size(); ar(n);
    if (Ar::reading) r.fracSnap.resize(size_t(n));
    for (auto& s : r.fracSnap) ar(s.ap, s.p, s.T, s.cond);
    n = r.wellFracLinks.size(); ar(n);
    if (Ar::reading) r.wellFracLinks.resize(size_t(n));
    for (auto& l : r.wellFracLinks) ar(l.first, l.second);
}

// C24 coupled fault mechanics: BLOB "faultmech"
template <class Ar> void fmIO(Ar& ar, faultmech::Output& o) {
    uint32_t ver = 2; ar(ver);   // 2: C26 FAULTSLIP slip arrays
    ar(o.ok, o.error, o.names, o.mMax, o.paPerUnit, o.mPerUnit, o.criticalAtStart);
    int64_t tw = o.twoWayUpdates; ar(tw); o.twoWayUpdates = long(tw);
    uint64_t n = o.faces.size(); ar(n);
    if (Ar::reading) o.faces.resize(size_t(n));
    for (auto& f : o.faces) {
        ar(f.fault, f.g, f.a, f.b, f.k, f.area, f.conns, f.p0, f.sn0, f.tau0, f.cff0, f.aSig);
        for (int i = 0; i < 3; ++i) ar(f.c[size_t(i)], f.n[size_t(i)]);
        for (auto& q : f.quad) for (int i = 0; i < 3; ++i) ar(q[size_t(i)]);
    }
    n = o.steps.size(); ar(n);
    if (Ar::reading) o.steps.resize(size_t(n));
    for (auto& s : o.steps) ar(s.t, s.cff, s.dcfs, s.st, s.sne, s.rate, s.cum, s.mult);
    n = o.series.size(); ar(n);
    if (Ar::reading) o.series.resize(size_t(n));
    for (auto& s : o.series) ar(s.t, s.maxDcfs, s.maxCff, s.maxSt, s.reactFrac, s.expected, s.events, s.maxMw);
    n = o.events.size(); ar(n);
    if (Ar::reading) o.events.resize(size_t(n));
    for (auto& e : o.events) ar(e.t, e.x, e.y, e.z, e.mw, e.fault, e.kind);
    n = o.rsf.size(); ar(n);
    if (Ar::reading) o.rsf.resize(size_t(n));
    for (auto& r : o.rsf) ar(r.fault, r.layer, r.faces, r.hStar, r.Lb, r.dx, r.seconds, r.events, r.error);
    if (ver >= 2) {
        ar(o.slip, o.halfSpace, o.shearModulus, o.ddmSeconds, o.maxSlipIter, o.slipNote);
        int64_t dp = o.ddmPairs, dd = o.ddmDropped; ar(dp, dd); o.ddmPairs = long(dp); o.ddmDropped = long(dd);
        for (auto& s : o.steps) ar(s.slip, s.dtauSlip, s.moment, s.momentStep, s.slipping);
        for (auto& s : o.series) ar(s.maxSlip, s.moment);
        for (auto& e : o.events) ar(e.faces);
    }
}
// C26 DDM apertures of the in-run hydraulic fractures: BLOB "hfddm"
template <class Ar> void hdIO(Ar& ar, hfddm::Output& o) {
    uint32_t ver = 1; ar(ver);
    ar(o.width, o.mult, o.pnet, o.maxWidth, o.volume, o.open, o.elements, o.updates, o.seconds);
    uint64_t n = o.steps.size(); ar(n);
    if (Ar::reading) o.steps.resize(size_t(n));
    for (auto& s : o.steps) ar(s.t, s.maxWidth, s.volume, s.minMult, s.maxMult, s.maxPnet, s.open, s.width);
}
// C28 in-run phase field: BLOB "hfpf"
template <class Ar> void pfIO(Ar& ar, std::vector<Simulator::PfRecord>& v) {
    uint32_t ver = 1; ar(ver);
    uint64_t n = v.size(); ar(n);
    if (Ar::reading) v.resize(size_t(n));
    for (auto& r : v) { int act = r.active ? 1 : 0; ar(r.t, r.frac, r.volume, r.pressure, r.halfLength, r.pFlow, r.sigma, r.stagger, act, r.note); r.active = act != 0; }
}
// C24 in-run G1-c contact coupling: BLOB "geocouple"
template <class Ar> void gcIO(Ar& ar, geocouple::Output& o) {
    uint32_t ver = 1; ar(ver);
    ar(o.ok, o.error, o.s, o.frac, o.efemSeconds);
    uint64_t n = o.steps.size(); ar(n);
    if (Ar::reading) o.steps.resize(size_t(n));
    for (auto& s : o.steps) ar(s.t, s.p, s.wn, s.wt, s.tn, s.a, s.mult, s.maxOpening, s.maxSlip, s.minMult, s.maxMult, s.openPieces, s.slipPieces, s.newton);
}

std::string nowIso() {
    const std::time_t t = std::time(nullptr);
    char b[32]; std::tm tmv{};
#ifdef _WIN32
    localtime_s(&tmv, &t);
#else
    localtime_r(&t, &tmv);
#endif
    std::strftime(b, sizeof b, "%Y-%m-%d %H:%M", &tmv);
    return b;
}
} // namespace

ressim::ResultsFile resultsToFile(const RunResults& r) {
    ResultsFile f;
    f.engine = r.engine; f.caseText = projectText(r.input); f.title = r.input.title; f.created = nowIso();
    f.summary = r.summary; f.snaps = r.snaps; f.wellNames = r.wellNames;
    for (auto t : r.wellTypes) f.wellTypes.push_back(int(t));
    f.components = r.components; f.trackName = r.trackName; f.initial = r.initial; f.stats = r.stats; f.totalTime = r.totalTime;
    if (!r.fracSegs.empty() || !r.wellT.empty() || !r.fracNames.empty()) {
        StateOut o; geoIO(o, const_cast<RunResults&>(r));
        f.blobs.emplace_back("geothermal", std::move(o.w.buf));
    }
    if (!r.opmDeck.empty()) f.blobs.emplace_back("opmDeck", r.opmDeck);
    if (r.faultMech) { StateOut o; faultmech::Output c = *r.faultMech; fmIO(o, c); f.blobs.emplace_back("faultmech", std::move(o.w.buf)); }
    if (r.hfDdm) { StateOut o; hfddm::Output c = *r.hfDdm; hdIO(o, c); f.blobs.emplace_back("hfddm", std::move(o.w.buf)); }
    if (r.hfPf) { StateOut o; std::vector<Simulator::PfRecord> c = *r.hfPf; pfIO(o, c); f.blobs.emplace_back("hfpf", std::move(o.w.buf)); }
    if (r.geoCouple) { StateOut o; geocouple::Output c = *r.geoCouple; gcIO(o, c); f.blobs.emplace_back("geocouple", std::move(o.w.buf)); }
    f.restarts = r.restarts;
    return f;
}

RunResults resultsFromFile(const ressim::ResultsFile& f) {
    RunResults r;
    r.engine = f.engine;
    if (!f.caseText.empty()) r.input = projectFromText(f.caseText);
    r.summary = f.summary; r.snaps = f.snaps; r.wellNames = f.wellNames;
    for (int t : f.wellTypes) r.wellTypes.push_back(WellType(t));
    r.components = f.components; r.trackName = f.trackName; r.initial = f.initial; r.stats = f.stats; r.totalTime = f.totalTime;
    if (const std::string* g = f.blob("geothermal")) { StateIn in(*g); geoIO(in, r); }
    if (const std::string* d = f.blob("opmDeck")) r.opmDeck = *d;
    if (const std::string* b = f.blob("faultmech")) { StateIn in(*b); auto o = std::make_shared<faultmech::Output>(); fmIO(in, *o); r.faultMech = o; }
    if (const std::string* b = f.blob("hfddm")) { StateIn in(*b); auto o = std::make_shared<hfddm::Output>(); hdIO(in, *o); r.hfDdm = o; }
    if (const std::string* b = f.blob("hfpf")) { StateIn in(*b); auto o = std::make_shared<std::vector<Simulator::PfRecord>>(); pfIO(in, *o); r.hfPf = o; }
    if (const std::string* b = f.blob("geocouple")) { StateIn in(*b); auto o = std::make_shared<geocouple::Output>(); gcIO(in, *o); r.geoCouple = o; }
    r.restarts = f.restarts;
    auto grid = std::make_shared<Grid>();
    grid->build(r.input);
    if (!r.snaps.empty() && !r.snaps[0].pressure.empty() && r.snaps[0].pressure.size() != size_t(grid->numActive()))
        throw std::runtime_error("the stored case does not reproduce the grid of the results (" + std::to_string(grid->numActive()) + " active cells, results have " +
                                 std::to_string(r.snaps[0].pressure.size()) + ")");
    r.grid = grid;
    r.valid = !r.summary.empty();
    return r;
}

void saveRunResults(const RunResults& r, const std::string& path) { saveResultsFile(resultsToFile(r), path); }
RunResults loadRunResults(const std::string& path) { return resultsFromFile(loadResultsFile(path)); }

std::string exportRunEclipse(const RunResults& r, const std::string& base, std::vector<std::string>* files, std::vector<std::string>* warnings) {
    if (!r.grid) return "no grid in the results";
    std::string err; ecl::WriteReport rep;
    if (!ecl::writeResultCase(base, r.input, *r.grid, r.summary, r.snaps, r.wellNames, &err, &rep)) return err.empty() ? "write failed" : err;
    if (files) *files = rep.files;
    if (warnings) *warnings = rep.warnings;
    return "";
}
