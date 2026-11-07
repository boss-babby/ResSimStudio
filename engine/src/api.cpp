// C33 scripting / automation layer (api.h)
#include "ressim/api.h"
#include "ressim/coupled_bo.h"
#include "ressim/examples.h"
#include "ressim/io.h"
#include "ressim/study.h"
#include <cctype>
#include <chrono>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace ressim {
namespace api {

namespace {
std::string upper(std::string s) { for (auto& ch : s) ch = char(std::toupper(static_cast<unsigned char>(ch))); return s; }
std::string dirOf(const std::string& p) {
    const size_t k = p.find_last_of("/\\");
    return k == std::string::npos ? std::string() : p.substr(0, k);
}
std::string extOf(const std::string& p) {
    const size_t k = p.find_last_of('.');
    const size_t s = p.find_last_of("/\\");
    if (k == std::string::npos || (s != std::string::npos && k < s)) return "";
    return upper(p.substr(k + 1));
}
}  // namespace

bool applySet(std::string& text, const std::string& kv, std::string& err) {
    const size_t eq = kv.find('=');
    if (eq == std::string::npos || eq == 0) { err = "set needs KEY=VALUE: " + kv; return false; }
    std::string key = upper(kv.substr(0, eq)), val = kv.substr(eq + 1);
    while (!key.empty() && key.back() == ' ') key.pop_back();
    std::vector<std::string> L;
    { size_t p = 0; while (p <= text.size()) { size_t e = text.find('\n', p); if (e == std::string::npos) e = text.size(); L.push_back(text.substr(p, e - p)); p = e + 1; } }
    auto firstTok = [](const std::string& l) {
        size_t a = l.find_first_not_of(" \t"); if (a == std::string::npos) return std::string();
        size_t b = l.find_first_of(" \t\r/", a);
        return upper(l.substr(a, b == std::string::npos ? std::string::npos : b - a));
    };
    auto join = [&]() { text.clear(); for (size_t i = 0; i < L.size(); ++i) { text += L[i]; if (i + 1 < L.size()) text += '\n'; } };
    if (key.rfind("RESSIM:", 0) == 0) {
        const std::string k = key.substr(7);
        for (auto& l : L) {
            if (l.rfind("-- RESSIM:", 0) != 0) continue;
            if (firstTok(l.substr(10)) == k) { l = "-- RESSIM: " + k + " " + val; join(); return true; }
        }
        size_t at = L.size();
        for (size_t i = 0; i < L.size(); ++i) if (L[i].rfind("-- RESSIM:", 0) == 0) at = i + 1;
        if (at == L.size()) for (size_t i = 0; i < L.size(); ++i) if (firstTok(L[i]) == "SCHEDULE") { at = i; break; }
        L.insert(L.begin() + long(at), "-- RESSIM: " + k + " " + val);
        join();
        return true;
    }
    for (size_t i = 0; i < L.size(); ++i) {
        if (L[i].rfind("--", 0) == 0 || firstTok(L[i]) != key) continue;
        size_t j = i + 1;
        for (; j < L.size(); ++j) {
            std::string l = L[j]; const size_t c = l.find("--"); if (c != std::string::npos) l = l.substr(0, c);
            if (l.find('/') != std::string::npos) break;
        }
        if (j >= L.size()) { err = "set: keyword " + key + " has no terminating '/'"; return false; }
        L.erase(L.begin() + long(i) + 1, L.begin() + long(j) + 1);
        L.insert(L.begin() + long(i) + 1, " " + val + " /");
        join();
        return true;
    }
    err = "set: keyword " + key + " not found in the deck";
    return false;
}

std::string engineName(const Case& c) {
    if (c.thermal.enabled) return "ResSim thermal";
    if (c.comp.enabled) return c.comp.co2store ? "ResSim CO2STORE" : "ResSim compositional";
    if (coupledbo::needed(c)) return "ResSim coupled black oil";
    return "ResSim black oil";
}

RunOutput runCase(const Case& c, const std::function<void(double)>& progress, std::atomic<bool>* cancel) {
    RunOutput r;
    r.engine = engineName(c);
    for (const auto& w : c.wells) r.wellNames.push_back(w.name);
    const auto t0 = std::chrono::steady_clock::now();
    const std::string bad = validateCase(c);
    if (!bad.empty()) { r.error = bad; return r; }
    if (c.thermal.enabled || c.comp.enabled) {
        RunResult rr = runSingle(c, cancel);
        r.ok = rr.ok; r.error = rr.error; r.summary = std::move(rr.summary); r.stats = rr.stats; r.kpi = rr.kpi; r.wellNames = rr.wellNames;
    } else {
        try {
            Simulator s;
            s.cancel = cancel;
            s.log = [&r](const std::string& m) { if (r.log.size() < 20000) r.log.push_back(m); };
            s.setup(c);
            s.initialize();
            r.initial = s.inPlace();
            if (coupledbo::needed(c)) {
                coupledbo::Output co;
                coupledbo::run(s, c, co, {}, cancel, progress);
                r.kpi = co.kpi;
                for (const auto& l : co.report) r.log.push_back(l);
            } else {
                while (!s.finished()) {
                    if (!s.runReportStep()) break;
                    if (progress) progress(s.time());
                }
            }
            r.summary = s.summary();
            r.snaps = s.snapshots();
            r.stats = s.stats();
            r.globalOf = s.grid().globalOf;
            r.depth = s.grid().depth;
            r.porv = s.grid().poreVolRef;
            r.wellNames = s.wellNames();
            r.ok = s.finished() && !(cancel && cancel->load());
            if (!r.ok) r.error = "cancelled / stopped early";
        } catch (const std::exception& e) { r.ok = false; r.error = e.what(); }
    }
    r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (r.ok && r.summary.empty()) { r.ok = false; r.error = "run produced no time steps"; }
    return r;
}

ecl::Summary toSummary(const RunOutput& r, const Case& c) {
    ecl::Summary s = ecl::summaryFromRun(r.summary, r.wellNames, c.units, ecl::componentNames(c));
    s.startDay = c.startDay; s.startMonth = c.startMonth; s.startYear = c.startYear;
    return s;
}

// ---------------------------------------------------------------- Session
Session::Session() = default;
Session::~Session() = default;

bool Session::open(const std::string& path) {
    error_.clear(); warnings_.clear(); out_ = RunOutput(); sim_.reset(); sumValid_ = false;
    try {
        if (path.rfind("example:", 0) == 0) {
            c_ = makeExample(path.substr(8));
            text_ = eclipseDeckToString(c_); baseDir_.clear(); fromProject_ = false;
            return true;
        }
        std::ifstream f(path, std::ios::binary);
        if (!f) { error_ = "cannot open " + path; return false; }
        std::stringstream ss; ss << f.rdbuf();
        text_ = ss.str(); baseDir_ = dirOf(path);
        fromProject_ = extOf(path) == "RSIM";
        if (fromProject_) { c_ = projectFromText(text_, baseDir_); return true; }
        std::vector<DeckMessage> rep;
        c_ = eclipseDeckFromString(text_, baseDir_, &warnings_, &rep);
        if (c_.title.empty() || c_.title == "Untitled") {
            std::string n = path.substr(path.find_last_of("/\\") == std::string::npos ? 0 : path.find_last_of("/\\") + 1);
            if (n.find('.') != std::string::npos) n = n.substr(0, n.find('.'));
            c_.title = n;
        }
        return true;
    } catch (const std::exception& e) { error_ = e.what(); return false; }
}

bool Session::openText(const std::string& deckText, const std::string& baseDir) {
    error_.clear(); warnings_.clear(); out_ = RunOutput(); sim_.reset(); sumValid_ = false;
    try {
        text_ = deckText; baseDir_ = baseDir; fromProject_ = false;
        std::vector<DeckMessage> rep;
        c_ = eclipseDeckFromString(text_, baseDir_, &warnings_, &rep);
        return true;
    } catch (const std::exception& e) { error_ = e.what(); return false; }
}

bool Session::setKeyword(const std::string& kv) {
    std::string t = text_, err;
    if (!applySet(t, kv, err)) { error_ = err; return false; }
    try {
        Case nc;
        std::vector<std::string> w;
        if (fromProject_) nc = projectFromText(t, baseDir_);
        else { std::vector<DeckMessage> rep; nc = eclipseDeckFromString(t, baseDir_, &w, &rep); }
        nc.title = c_.title;
        c_ = std::move(nc); text_ = t; warnings_ = w;
        out_ = RunOutput(); sim_.reset(); sumValid_ = false;
        return true;
    } catch (const std::exception& e) { error_ = e.what(); return false; }
}

bool Session::run(const std::function<void(double)>& progress, std::atomic<bool>* cancel) {
    sim_.reset(); sumValid_ = false;
    out_ = runCase(c_, progress, cancel);
    if (!out_.ok) error_ = out_.error;
    return out_.ok;
}

bool Session::begin() {
    sumValid_ = false; out_ = RunOutput(); out_.engine = engineName(c_);
    if (c_.thermal.enabled || c_.comp.enabled) { error_ = "step-wise runs are available for the black-oil engine (use run() for " + out_.engine + ")"; return false; }
    const std::string bad = validateCase(c_);
    if (!bad.empty()) { error_ = bad; return false; }
    try {
        sim_ = std::make_unique<Simulator>();
        RunOutput* o = &out_;
        sim_->log = [o](const std::string& m) { if (o->log.size() < 20000) o->log.push_back(m); };
        sim_->setup(c_);
        sim_->initialize();
        out_.initial = sim_->inPlace();
        out_.wellNames = sim_->wellNames();
        out_.globalOf = sim_->grid().globalOf; out_.depth = sim_->grid().depth; out_.porv = sim_->grid().poreVolRef;
        t0_ = 0;
        finishStepping();
        return true;
    } catch (const std::exception& e) { error_ = e.what(); sim_.reset(); return false; }
}

void Session::finishStepping() {
    if (!sim_) return;
    out_.summary = sim_->summary();
    out_.wellNames = sim_->wellNames();      // C33: wells created at run time (ACTIONX / script WELSPECS)
    out_.snaps = sim_->snapshots();
    out_.stats = sim_->stats();
    sumValid_ = false;
}

int Session::step() {
    if (!sim_) { error_ = "call begin() first"; return -1; }
    if (sim_->finished()) return 0;
    const auto t0 = std::chrono::steady_clock::now();
    try {
        const bool more = sim_->runReportStep();
        out_.seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        finishStepping();
        out_.ok = sim_->finished();
        if (!more && !sim_->finished()) { error_ = "the run stopped early"; return -1; }
        return sim_->finished() ? 0 : 1;
    } catch (const std::exception& e) { error_ = e.what(); return -1; }
}

bool Session::applySchedule(const std::string& text, std::vector<std::string>* messages) {
    if (!sim_) { error_ = "call begin() first"; return false; }
    if (sim_->finished()) { error_ = "the schedule is finished"; return false; }
    try {
        std::vector<std::string> m = sim_->applySchedule(text);
        if (messages) *messages = m;
        return true;
    } catch (const std::exception& e) { error_ = e.what(); return false; }
}

double Session::time() const {
    if (sim_) return sim_->time();
    return out_.summary.empty() ? 0.0 : out_.summary.back().time;
}
int Session::reportStep() const { return sim_ ? sim_->reportStep() : (hasResults() ? numReportSteps() : 0); }

const ecl::Summary& Session::summary() {
    if (!sumValid_) { sum_ = toSummary(out_, c_); sumValid_ = true; }
    return sum_;
}

std::vector<double> Session::staticArray(const std::string& nameIn) const {
    const std::string n = upper(nameIn);
    const GridInput& g = c_.grid;
    auto d = [](const std::vector<double>& v) { return v; };
    auto i2d = [](const std::vector<int>& v) { return std::vector<double>(v.begin(), v.end()); };
    if (n == "PORO") return d(g.poro);
    if (n == "PERMX") return d(g.permx);
    if (n == "PERMY") return d(g.permy);
    if (n == "PERMZ") return d(g.permz);
    if (n == "NTG") return d(g.ntg);
    if (n == "DX") return d(g.dx);
    if (n == "DY") return d(g.dy);
    if (n == "DZ") return d(g.dz);
    if (n == "TOPS") return d(g.tops);
    if (n == "ACTNUM") return i2d(g.actnum);
    if (n == "SATNUM") return i2d(g.satnum);
    if (n == "PVTNUM") return i2d(g.pvtnum);
    if (n == "EQLNUM") return i2d(g.eqlnum);
    if (n == "FIPNUM") return i2d(g.fipnum);
    return {};
}

std::vector<float> Session::cellArray(int snap, const std::string& name) const {
    if (snap < 0) snap += int(out_.snaps.size());
    if (snap < 0 || snap >= int(out_.snaps.size())) return {};
    return out_.snaps[size_t(snap)].field(upper(name));
}

std::vector<std::string> Session::cellArrayNames(int snap) const {
    std::vector<std::string> r;
    if (snap < 0) snap += int(out_.snaps.size());
    if (snap < 0 || snap >= int(out_.snaps.size())) return r;
    const Snapshot& s = out_.snaps[size_t(snap)];
    const char* base[] = {"PRESSURE", "SWAT", "SOIL", "SGAS", "RS", "RV", "POLYMER", "PADS", "XCO2"};
    for (const char* b : base) if (!s.field(b).empty()) r.push_back(b);
    for (const auto& e : s.extra) r.push_back(e.first);
    return r;
}

bool Session::writeEclipse(const std::string& base, std::string* err) const {
    if (out_.summary.empty()) { if (err) *err = "no results"; return false; }
    try {
        Grid grid;
        if (sim_) grid = sim_->grid();
        else {
            Simulator s; s.setup(c_); grid = s.grid();
        }
        return ecl::writeResultCase(base, c_, grid, out_.summary, out_.snaps, out_.wellNames, err);
    } catch (const std::exception& e) { if (err) *err = e.what(); return false; }
}

bool Session::exportDeck(const std::string& path, std::string* err) const {
    try { exportEclipseDeck(c_, path); return true; }
    catch (const std::exception& e) { if (err) *err = e.what(); return false; }
}

}  // namespace api
}  // namespace ressim
