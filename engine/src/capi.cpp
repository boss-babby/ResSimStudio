// C33 C API (capi.h) on top of api::Session
#define RESSIM_CAPI_BUILD
#include "ressim/capi.h"
#include "ressim/api.h"
#include "ressim/examples.h"
#include "ressim/report.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <sstream>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace ressim;

struct rs_session {
    std::unique_ptr<api::Session> s;         // null for summary-only sessions
    ecl::Summary ext;                        // rs_open_summary
    bool summaryOnly = false;
    std::string buf;                         // returned strings
    const ecl::Summary& summary() { return summaryOnly ? ext : s->summary(); }
};

namespace {
thread_local std::string g_err;
int fail(const std::string& e) { g_err = e; return -1; }
std::string jsonEsc(const std::string& s) {
    std::string o;
    for (unsigned char c : s) {
        if (c == '"' || c == '\\') { o += '\\'; o += char(c); }
        else if (c == '\n') o += "\\n";
        else if (c < 32) o += ' ';
        else o += char(c);
    }
    return o;
}
std::string jnum(double v) {
    if (!std::isfinite(v)) return "null";
    char b[40]; std::snprintf(b, sizeof b, "%.12g", v); return b;
}
#define NEED_SESSION(s) if (!(s)) return fail("null session"); if ((s)->summaryOnly) return fail("summary-only session")
}  // namespace

void rs_capi_set_error(const std::string& e) { g_err = e; }   // C43: used by capi_geo.cpp

extern "C" {

const char* rs_version(void) { return "ResSim C API 1.0 (engine C33)"; }
const char* rs_last_error(void) { return g_err.c_str(); }
void rs_set_threads(int n) {
#ifdef _OPENMP
    if (n > 0) omp_set_num_threads(n);
#else
    (void)n;
#endif
}

rs_session* rs_open(const char* path) {
    if (!path) { g_err = "null path"; return nullptr; }
    auto r = std::make_unique<rs_session>();
    r->s = std::make_unique<api::Session>();
    if (!r->s->open(path)) { g_err = r->s->error(); return nullptr; }
    return r.release();
}
rs_session* rs_open_text(const char* text, const char* dir) {
    if (!text) { g_err = "null deck text"; return nullptr; }
    auto r = std::make_unique<rs_session>();
    r->s = std::make_unique<api::Session>();
    if (!r->s->openText(text, dir ? dir : "")) { g_err = r->s->error(); return nullptr; }
    return r.release();
}
rs_session* rs_open_summary(const char* smspec) {
    if (!smspec) { g_err = "null path"; return nullptr; }
    auto r = std::make_unique<rs_session>();
    r->summaryOnly = true;
    std::string err;
    if (!ecl::readSummary(smspec, r->ext, &err)) { g_err = err.empty() ? std::string("cannot read ") + smspec : err; return nullptr; }
    return r.release();
}
void rs_close(rs_session* s) { delete s; }

int rs_set_keyword(rs_session* s, const char* kv) {
    NEED_SESSION(s);
    if (!kv) return fail("null argument");
    return s->s->setKeyword(kv) ? 0 : fail(s->s->error());
}
int rs_num_warnings(rs_session* s) { NEED_SESSION(s); return int(s->s->warnings().size()); }
const char* rs_warning(rs_session* s, int i) {
    if (!s || s->summaryOnly || i < 0 || i >= int(s->s->warnings().size())) return "";
    s->buf = s->s->warnings()[size_t(i)];
    return s->buf.c_str();
}
const char* rs_info_json(rs_session* s) {
    if (!s) { g_err = "null session"; return ""; }
    std::ostringstream o;
    if (s->summaryOnly) {
        const ecl::Summary& S = s->ext;
        o << "{\"summary_only\": true, \"vectors\": " << S.vars.size() << ", \"ministeps\": " << S.time.size() << ", \"units\": \"" << jsonEsc(S.unitSystem)
          << "\", \"dims\": [" << S.nx << ", " << S.ny << ", " << S.nz << "], \"start\": [" << S.startYear << ", " << S.startMonth << ", " << S.startDay << "]}";
    } else {
        const Case& c = s->s->caseData();
        double tEnd = 0;
        std::ostringstream steps;
        for (size_t i = 0; i < c.schedule.size(); ++i) { tEnd += c.schedule[i].length; steps << (i ? ", " : "") << jnum(tEnd); }
        o << "{\"title\": \"" << jsonEsc(c.title) << "\", \"units\": \"" << (c.units == UnitSystem::Field ? "FIELD" : "METRIC") << "\", \"engine\": \""
          << jsonEsc(api::engineName(c)) << "\", \"dims\": [" << c.grid.nx << ", " << c.grid.ny << ", " << c.grid.nz << "], \"corner_point\": "
          << (c.grid.zcorn.empty() ? "false" : "true") << ", \"start\": [" << c.startYear << ", " << c.startMonth << ", " << c.startDay << "], \"report_times\": [" << steps.str()
          << "], \"wells\": [";
        for (size_t i = 0; i < c.wells.size(); ++i) {
            const auto& w = c.wells[i];
            o << (i ? ", " : "") << "{\"name\": \"" << jsonEsc(w.name) << "\", \"type\": \""
              << (w.type == WellType::Producer ? "producer" : w.type == WellType::WaterInjector ? "water_injector" : "gas_injector") << "\", \"connections\": " << w.completions.size() << "}";
        }
        o << "], \"faults\": " << c.grid.faults.size() << ", \"unsupported\": [";
        for (size_t i = 0; i < c.engineUnsupported.size(); ++i) o << (i ? ", " : "") << "\"" << jsonEsc(c.engineUnsupported[i]) << "\"";
        o << "], \"warnings\": " << s->s->warnings().size() << "}";
    }
    s->buf = o.str();
    return s->buf.c_str();
}
int rs_list_examples(char* buf, int len) {
    std::string all;
    for (const auto& e : listExamples()) { all += e.id; all += '\n'; }
    if (buf && len > 0) { const size_t n = std::min(all.size(), size_t(len - 1)); std::memcpy(buf, all.data(), n); buf[n] = 0; }
    return int(all.size()) + 1;
}

int rs_run(rs_session* s) { NEED_SESSION(s); return s->s->run() ? 0 : fail(s->s->error()); }
int rs_begin(rs_session* s) { NEED_SESSION(s); return s->s->begin() ? 0 : fail(s->s->error()); }
int rs_step(rs_session* s) {
    NEED_SESSION(s);
    const int r = s->s->step();
    if (r < 0) return fail(s->s->error());
    return r;
}
int rs_apply_schedule(rs_session* s, const char* text) {
    NEED_SESSION(s);
    if (!text) return fail("null text");
    std::vector<std::string> msgs;
    if (!s->s->applySchedule(text, &msgs)) return fail(s->s->error());
    if (!msgs.empty()) { g_err.clear(); for (const auto& m : msgs) g_err += m + "\n"; }
    return int(msgs.size());
}
double rs_time(rs_session* s) { if (!s || s->summaryOnly) return 0.0; return s->s->time(); }
int rs_report_step(rs_session* s) { NEED_SESSION(s); return s->s->reportStep(); }
int rs_num_report_steps(rs_session* s) { NEED_SESSION(s); return s->s->numReportSteps(); }
const char* rs_run_json(rs_session* s) {
    if (!s || s->summaryOnly) { g_err = "no run"; return ""; }
    const api::RunOutput& r = s->s->output();
    std::ostringstream o;
    o << "{\"ok\": " << (r.ok ? "true" : "false") << ", \"error\": \"" << jsonEsc(r.error) << "\", \"engine\": \"" << jsonEsc(r.engine) << "\", \"seconds\": " << jnum(r.seconds)
      << ", \"timesteps\": " << r.stats.timesteps << ", \"newton\": " << r.stats.newton << ", \"linear\": " << r.stats.linear << ", \"cuts\": " << r.stats.cuts
      << ", \"snapshots\": " << r.snaps.size() << ", \"initial_in_place\": {\"oil\": " << jnum(r.initial.oil) << ", \"water\": " << jnum(r.initial.water)
      << ", \"gas_free\": " << jnum(r.initial.gasFree) << ", \"gas_dissolved\": " << jnum(r.initial.gasDissolved) << "}, \"kpi\": {";
    for (size_t i = 0; i < r.kpi.size(); ++i) o << (i ? ", " : "") << "\"" << jsonEsc(r.kpi[i].first) << "\": " << jnum(r.kpi[i].second);
    o << "}, \"log\": [";
    const size_t from = r.log.size() > 200 ? r.log.size() - 200 : 0;
    for (size_t i = from; i < r.log.size(); ++i) o << (i > from ? ", " : "") << "\"" << jsonEsc(r.log[i]) << "\"";
    o << "]}";
    s->buf = o.str();
    return s->buf.c_str();
}

int rs_summary_count(rs_session* s) { if (!s) return fail("null session"); return int(s->summary().vars.size()); }
const char* rs_summary_key(rs_session* s, int i) {
    if (!s) return "";
    const auto& S = s->summary();
    if (i < 0 || i >= int(S.vars.size())) return "";
    s->buf = S.vars[size_t(i)].key;
    return s->buf.c_str();
}
const char* rs_summary_unit(rs_session* s, int i) {
    if (!s) return "";
    const auto& S = s->summary();
    if (i < 0 || i >= int(S.vars.size())) return "";
    s->buf = S.vars[size_t(i)].unit;
    return s->buf.c_str();
}
int rs_summary_length(rs_session* s) {
    if (!s) return fail("null session");
    const auto& S = s->summary();
    return S.values.empty() ? 0 : int(S.values.front().size());
}
int rs_summary_get(rs_session* s, const char* key, double* buf, int len) {
    if (!s || !key) return fail("null argument");
    const auto& S = s->summary();
    const int i = S.find(key);
    std::vector<double> v;
    if (i >= 0) v = S.series(i);
    else if (!std::strcmp(key, "TIME")) v = S.time;
    else return fail(std::string("no summary vector ") + key);
    if (buf) for (int k = 0; k < len && k < int(v.size()); ++k) buf[k] = v[size_t(k)];
    return int(v.size());
}

int rs_grid_dims(rs_session* s, int* nx, int* ny, int* nz, int* na) {
    NEED_SESSION(s);
    const Case& c = s->s->caseData();
    if (nx) *nx = c.grid.nx;
    if (ny) *ny = c.grid.ny;
    if (nz) *nz = c.grid.nz;
    if (na) *na = s->s->output().globalOf.empty() ? -1 : int(s->s->output().globalOf.size());
    return 0;
}
int rs_static_array(rs_session* s, const char* name, double* buf, int len) {
    NEED_SESSION(s);
    if (!name) return fail("null name");
    const auto v = s->s->staticArray(name);
    if (v.empty()) return fail(std::string("no static array ") + name);
    if (buf) for (int k = 0; k < len && k < int(v.size()); ++k) buf[k] = v[size_t(k)];
    return int(v.size());
}
int rs_active_index(rs_session* s, int* buf, int len) {
    NEED_SESSION(s);
    const auto& g = s->s->output().globalOf;
    if (g.empty()) return fail("run or begin first (the active cells are known after setup)");
    if (buf) for (int k = 0; k < len && k < int(g.size()); ++k) buf[k] = g[size_t(k)];
    return int(g.size());
}
int rs_num_snapshots(rs_session* s) { NEED_SESSION(s); return int(s->s->output().snaps.size()); }
double rs_snapshot_time(rs_session* s, int i) {
    if (!s || s->summaryOnly) return std::nan("");
    const auto& v = s->s->output().snaps;
    if (i < 0) i += int(v.size());
    if (i < 0 || i >= int(v.size())) return std::nan("");
    return v[size_t(i)].time;
}
const char* rs_cell_array_names(rs_session* s, int snap) {
    if (!s || s->summaryOnly) return "";
    s->buf.clear();
    for (const auto& n : s->s->cellArrayNames(snap)) { s->buf += n; s->buf += '\n'; }
    return s->buf.c_str();
}
int rs_cell_array(rs_session* s, int snap, const char* name, float* buf, int len) {
    NEED_SESSION(s);
    if (!name) return fail("null name");
    const auto v = s->s->cellArray(snap, name);
    if (v.empty()) return fail(std::string("no cell array ") + name + " in snapshot " + std::to_string(snap));
    if (buf) for (int k = 0; k < len && k < int(v.size()); ++k) buf[k] = v[size_t(k)];
    return int(v.size());
}

int rs_write_eclipse(rs_session* s, const char* base) {
    NEED_SESSION(s);
    if (!base) return fail("null path");
    std::string err;
    return s->s->writeEclipse(base, &err) ? 0 : fail(err);
}
int rs_export_deck(rs_session* s, const char* path) {
    NEED_SESSION(s);
    if (!path) return fail("null path");
    std::string err;
    return s->s->exportDeck(path, &err) ? 0 : fail(err);
}
int rs_report(rs_session* s, const char* path, const char* options) {
    NEED_SESSION(s);
    if (!path) return fail("null path");
    if (!s->s->hasResults()) return fail("run the case first");
    const report::Options opt = report::parseOptions(options ? options : "");
    const report::Doc d = report::simulationReport(s->s->caseData(), s->s->output(), s->s->summary(), opt);
    std::string err;
    return report::writeDoc(path, d, err) ? 0 : fail(err);
}

const char* rs_compare_json(rs_session* ref, rs_session* cmp, const char* keys) {
    if (!ref || !cmp) { g_err = "null session"; return ""; }
    std::vector<std::string> k;
    if (keys && *keys) {
        std::string t = keys; size_t a = 0;
        while (a <= t.size()) { size_t e = t.find(',', a); if (e == std::string::npos) e = t.size(); if (e > a) k.push_back(t.substr(a, e - a)); a = e + 1; }
    }
    const auto d = ecl::compareSummaries(ref->summary(), cmp->summary(), k);
    std::ostringstream o;
    o << "[";
    for (size_t i = 0; i < d.size(); ++i)
        o << (i ? ", " : "") << "{\"key\": \"" << jsonEsc(d[i].key) << "\", \"ref_final\": " << jnum(d[i].refFinal) << ", \"cmp_final\": " << jnum(d[i].cmpFinal)
          << ", \"max_abs\": " << jnum(d[i].maxAbs) << ", \"rms\": " << jnum(d[i].rms) << ", \"rel_rms\": " << jnum(d[i].relRms) << ", \"samples\": " << d[i].samples << "}";
    o << "]";
    cmp->buf = o.str();
    return cmp->buf.c_str();
}

}  // extern "C"
