// C33 deck regression harness (regress.h)
#include "ressim/regress.h"
#include "ressim/api.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <set>
#include <thread>
#include <atomic>
#include <sstream>

namespace ressim {
namespace regress {

namespace {
std::string trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}
std::vector<std::string> splitCsv(const std::string& s) {
    std::vector<std::string> r;
    size_t a = 0;
    while (a <= s.size()) {
        size_t e = s.find(',', a);
        if (e == std::string::npos) e = s.size();
        const std::string t = trim(s.substr(a, e - a));
        if (!t.empty()) r.push_back(t);
        a = e + 1;
    }
    return r;
}
std::string join(const std::string& dir, const std::string& p) {
    if (p.empty() || dir.empty() || p[0] == '/' || p[0] == '\\' || (p.size() > 1 && p[1] == ':')) return p;
    return dir + "/" + p;
}
std::string fmtv(double v, int prec = 4) {
    if (!std::isfinite(v)) return "-";
    char b[40];
    std::snprintf(b, sizeof b, "%.*g", prec, v);
    return b;
}
std::string xesc(const std::string& s) {
    std::string o;
    for (char c : s) {
        if (c == '&') o += "&amp;"; else if (c == '<') o += "&lt;"; else if (c == '>') o += "&gt;"; else if (c == '"') o += "&quot;";
        else o += c;
    }
    return o;
}
std::string nowStamp() {
    std::time_t t = std::time(nullptr);
    char b[32];
    std::strftime(b, sizeof b, "%Y-%m-%dT%H:%M", std::localtime(&t));
    return b;
}
}  // namespace

bool readManifest(const std::string& path, Manifest& m, std::string& err) {
    std::ifstream f(path);
    if (!f) { err = "cannot open manifest " + path; return false; }
    const size_t k = path.find_last_of("/\\");
    m.dir = k == std::string::npos ? "." : path.substr(0, k);
    std::string line;
    int ln = 0;
    while (std::getline(f, line)) {
        ++ln;
        const size_t h = line.find('#');
        if (h != std::string::npos) line = line.substr(0, h);
        line = trim(line);
        if (line.empty()) continue;
        std::istringstream is(line);
        std::string cmd;
        is >> cmd;
        if (cmd == "tol") { is >> m.tol; }
        else if (cmd == "keys") { std::string rest; std::getline(is, rest); m.keys = splitCsv(rest); }
        else if (cmd == "case") {
            CaseSpec c;
            is >> c.name;
            std::string tok;
            while (is >> tok) {
                const size_t q = tok.find('=');
                if (q == std::string::npos) { err = path + ":" + std::to_string(ln) + ": expected key=value, got '" + tok + "'"; return false; }
                const std::string key = tok.substr(0, q), val = tok.substr(q + 1);
                if (key == "deck") c.deck = join(m.dir, val);
                else if (key == "ref") c.reference = join(m.dir, val);
                else if (key == "tol") c.tol = std::atof(val.c_str());
                else if (key == "keys") c.keys = splitCsv(val);
                else if (key == "group") c.group = val;
                else { err = path + ":" + std::to_string(ln) + ": unknown item '" + key + "'"; return false; }
            }
            if (c.name.empty() || c.deck.empty()) { err = path + ":" + std::to_string(ln) + ": case needs a name and deck="; return false; }
            m.cases.push_back(c);
        } else if (cmd == "skip") {
            CaseSpec c;
            is >> c.name;
            std::getline(is, c.skipReason);
            c.skipReason = trim(c.skipReason);
            if (c.skipReason.empty()) c.skipReason = "listed as skip";
            m.cases.push_back(c);
        } else { err = path + ":" + std::to_string(ln) + ": unknown directive '" + cmd + "'"; return false; }
    }
    return true;
}

void compare(CaseResult& r, const std::vector<std::string>& keysIn, double tol) {
    r.vectors.clear();
    r.worst = 0; r.worstKey.clear();
    std::vector<std::string> keys;
    for (const auto& k : keysIn) {
        if (k.size() > 2 && k.substr(k.size() - 2) == ":*") {
            const std::string kw = k.substr(0, k.size() - 2);
            for (const auto& v : r.ref.vars)
                if (v.keyword == kw && !v.wgname.empty() && r.run.find(v.key) >= 0 && std::find(keys.begin(), keys.end(), v.key) == keys.end()) keys.push_back(v.key);
        } else if (r.ref.find(k) >= 0 && r.run.find(k) >= 0) keys.push_back(k);
    }
    const auto diffs = ecl::compareSummaries(r.ref, r.run, keys);
    bool pass = true;
    for (const auto& d : diffs) {
        VectorResult v;
        v.key = d.key; v.relRms = d.relRms; v.refFinal = d.refFinal; v.cmpFinal = d.cmpFinal;
        // vectors that are zero in both runs carry no information
        const int ri = r.ref.find(d.key);
        double mx = 0;
        if (ri >= 0) for (float x : r.ref.values[size_t(ri)]) mx = std::max(mx, double(std::fabs(x)));
        if (mx < 1e-12 && std::fabs(d.cmpFinal) < 1e-9) continue;
        v.pass = std::isfinite(d.relRms) && d.relRms <= tol;
        pass = pass && v.pass;
        if (!std::isfinite(d.relRms) || d.relRms > r.worst) { r.worst = std::isfinite(d.relRms) ? d.relRms : 1e9; r.worstKey = d.key; }
        r.vectors.push_back(v);
    }
    // C41 "material" deviation (reported next to the strict metric, does not change the status): liquid-phase rate /
    // total vectors of a well (or the field) are also normalised by that well's liquid scale (max of oil + water), so a
    // water rate of a few sm3/d next to thousands of sm3/d of oil is judged against the well's production - the same
    // idea as the absolute + relative tolerances of the OPM / ECLIPSE summary comparison tools
    double matWorst = 0; std::string matKey;
    {
        auto seriesMax = [&](const std::string& k) {
            const int i = r.ref.find(k);
            double mx = 0;
            if (i >= 0) for (float x : r.ref.values[size_t(i)]) mx = std::max(mx, double(std::fabs(x)));
            return mx;
        };
        auto sumMax = [&](const std::string& k1, const std::string& k2) {
            const int i = r.ref.find(k1), j = r.ref.find(k2);
            double mx = 0;
            if (i >= 0 && j >= 0)
                for (size_t t = 0; t < r.ref.values[size_t(i)].size() && t < r.ref.values[size_t(j)].size(); ++t)
                    mx = std::max(mx, double(std::fabs(r.ref.values[size_t(i)][t]) + std::fabs(r.ref.values[size_t(j)][t])));
            else mx = std::max(seriesMax(k1), seriesMax(k2));
            return mx;
        };
        for (const auto& d : diffs) {
            const size_t c = d.key.find(':');
            const std::string kw = d.key.substr(0, c), wg = c == std::string::npos ? std::string() : d.key.substr(c);
            double rel = d.relRms;
            if (kw.size() == 4 && (kw[0] == 'W' || kw[0] == 'F' || kw[0] == 'G') && (kw[1] == 'O' || kw[1] == 'W' || kw[1] == 'L') && (kw[3] == 'R' || kw[3] == 'T') && kw[2] == 'P') {
                const std::string a = std::string(1, kw[0]) + "OP" + kw[3] + wg, b = std::string(1, kw[0]) + "WP" + kw[3] + wg;
                const double liq = sumMax(a, b), own = d.relRms > 0 ? d.rms / d.relRms : 0.0;
                const double sc = std::max(own, liq);
                rel = sc > 1e-4 ? d.rms / sc : 0.0;
            }
            if (!std::isfinite(rel)) rel = 1e9;
            const int ri = r.ref.find(d.key);
            double mx = 0;
            if (ri >= 0) for (float x : r.ref.values[size_t(ri)]) mx = std::max(mx, double(std::fabs(x)));
            if (mx < 1e-12 && std::fabs(d.cmpFinal) < 1e-9) continue;
            if (rel > matWorst) { matWorst = rel; matKey = d.key; }
        }
    }
    std::sort(r.vectors.begin(), r.vectors.end(), [](const VectorResult& a, const VectorResult& b) { return a.relRms > b.relRms; });
    if (r.vectors.empty()) {
        if (!diffs.empty()) { r.status = "PASS"; r.message = std::to_string(diffs.size()) + " common vectors, all zero in both runs"; return; }
        r.status = "SKIP"; r.message = "no common summary vectors with the reference (grid / INIT-only test or empty SUMMARY section)"; return;
    }
    r.status = pass ? "PASS" : "FAIL";
    if (!pass) {
        int nf = 0;
        for (const auto& v : r.vectors) nf += v.pass ? 0 : 1;
        r.message = std::to_string(nf) + " of " + std::to_string(r.vectors.size()) + " vectors above " + fmtv(100 * tol, 3) + " %";
    }
    if (!matKey.empty()) r.message += (r.message.empty() ? "" : "; ") + std::string("material worst ") + matKey + " " + fmtv(100 * matWorst, 3) + " %";
}

std::vector<CaseResult> runAll(const Manifest& m, const Options& o) {
    std::vector<CaseResult> out;
    const std::string hist = o.history.empty() ? o.outDir + "/history.tsv" : o.history;
    const auto prev = readHistory(hist);
    auto log = [&](const std::string& s) { if (o.log) o.log(s); };
    for (const CaseSpec& c : m.cases) {
        if (!o.filter.empty() && c.name.find(o.filter) == std::string::npos) continue;
        if (o.cancel && o.cancel->load()) { log("cancelled by the user - remaining cases not run"); break; }
        CaseResult r;
        r.name = c.name; r.group = c.group; r.deck = c.deck;
        r.tol = o.tolOverride > 0 ? o.tolOverride : c.tol > 0 ? c.tol : m.tol;
        if (auto it = prev.find(c.name); it != prev.end()) r.prevSeconds = it->second;
        if (!c.skipReason.empty()) { r.status = "SKIP"; r.message = c.skipReason; out.push_back(std::move(r)); log("SKIP  " + c.name + "  " + c.skipReason); continue; }
        if (!c.reference.empty()) {
            std::string err;
            if (!ecl::readSummary(c.reference, r.ref, &err)) { r.status = "ERROR"; r.message = "reference: " + (err.empty() ? c.reference : err); out.push_back(std::move(r)); log("ERROR " + c.name + "  " + out.back().message); continue; }
        }
        api::Session s;
        if (!s.open(c.deck)) { r.status = "ERROR"; r.message = "import: " + s.error(); out.push_back(std::move(r)); log("ERROR " + c.name + "  " + out.back().message); continue; }
        if (!s.caseData().engineUnsupported.empty()) {
            r.status = "SKIP";
            r.message = "not supported by the ResSim engine: " + s.caseData().engineUnsupported.front();
            r.engine = api::engineName(s.caseData());
            out.push_back(std::move(r));
            log("SKIP  " + c.name + "  " + out.back().message);
            continue;
        }
        const auto t0 = std::chrono::steady_clock::now();
        std::atomic<bool> cancel{false}, done{false}, byUser{false};
        std::thread watchdog;
        if (o.maxSeconds > 0 || o.cancel)
            watchdog = std::thread([&]() {
                while (!done.load()) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(100));
                    if (o.cancel && o.cancel->load()) { byUser = true; cancel = true; break; }
                    if (o.maxSeconds > 0 && std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() > o.maxSeconds) { cancel = true; break; }
                }
            });
        const bool ok = s.run({}, &cancel);
        done = true;
        if (watchdog.joinable()) watchdog.join();
        if (byUser.load()) {
            r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            r.status = "ERROR"; r.message = "cancelled by the user"; out.push_back(std::move(r));
            log("ERROR " + c.name + "  cancelled - remaining cases not run");
            break;
        }
        if (cancel.load()) { r.seconds = o.maxSeconds; r.status = "ERROR"; r.message = "timeout after " + fmtv(o.maxSeconds, 4) + " s"; out.push_back(std::move(r)); log("ERROR " + c.name + "  timeout"); continue; }
        r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        r.engine = s.output().engine;
        r.steps = s.output().stats.timesteps; r.newton = s.output().stats.newton;
        if (!ok) { r.status = "ERROR"; r.message = "run: " + s.error(); out.push_back(std::move(r)); log("ERROR " + c.name + "  " + out.back().message); continue; }
        r.run = s.summary();
        if (c.reference.empty()) { r.status = "PASS"; r.message = "ran (no reference)"; }
        else compare(r, c.keys.empty() ? m.keys : c.keys, r.tol);
        r.slow = r.prevSeconds > 0.5 && r.seconds > o.slowFactor * r.prevSeconds;
        log(r.status + std::string(r.status.size() < 5 ? 5 - r.status.size() : 0, ' ') + " " + c.name + "  worst " + r.worstKey + " " + fmtv(100 * r.worst, 3) + " %  " +
            fmtv(r.seconds, 3) + " s" + (r.slow ? "  SLOWER than " + fmtv(r.prevSeconds, 3) + " s" : "") +
            (r.message.find("material worst") != std::string::npos ? "  [" + r.message.substr(r.message.find("material worst")) + "]" : ""));
        out.push_back(std::move(r));
    }
    return out;
}

report::Doc reportDoc(const std::vector<CaseResult>& R, const Options& o) {
    report::Doc d;
    d.title = "Deck regression report";
    d.subtitle = "ResSim engine vs reference summaries (OPM Flow / ECLIPSE)";
    int np = 0, nf = 0, ne = 0, ns = 0, nslow = 0;
    for (const auto& r : R) { np += r.status == "PASS"; nf += r.status == "FAIL"; ne += r.status == "ERROR"; ns += r.status == "SKIP"; nslow += r.slow; }
    d.h(1, "Summary");
    d.p(std::to_string(R.size()) + " cases: " + std::to_string(np) + " pass, " + std::to_string(nf) + " fail, " + std::to_string(ne) + " error, " + std::to_string(ns) +
        " skipped; " + std::to_string(nslow) + " slower than " + fmtv(100 * (o.slowFactor - 1), 3) + " % over the previous run.");
    d.p("Criterion: for every compared vector the RMS deviation over the reference report times divided by the largest |reference value| must not exceed the case tolerance. "
        "Vectors that are zero in both runs are not counted.");
    report::Table t;
    t.header = {"Case", "Status", "Worst vector", "rel. RMS (%)", "Tol (%)", "Vectors", "Time (s)", "Previous (s)", "Note"};
    t.caption = "All cases (red rows: FAIL / ERROR or slower)";
    for (size_t i = 0; i < R.size(); ++i) {
        const auto& r = R[i];
        t.rows.push_back({r.name, r.status, r.worstKey, r.status == "PASS" || r.status == "FAIL" ? fmtv(100 * r.worst, 3) : "", fmtv(100 * r.tol, 3),
                          std::to_string(r.vectors.size()), r.seconds > 0 ? fmtv(r.seconds, 3) : "", r.prevSeconds > 0 ? fmtv(r.prevSeconds, 3) : "", r.message});
        if (r.status == "FAIL" || r.status == "ERROR" || r.slow) t.highlightRows.push_back(int(i));
    }
    d.table(t);
    // detail: failures first, then passes, up to chartCases
    std::vector<const CaseResult*> order;
    for (const auto& r : R) if (r.status == "FAIL") order.push_back(&r);
    for (const auto& r : R) if (r.status == "PASS" && !r.vectors.empty()) order.push_back(&r);
    int shown = 0;
    if (!order.empty()) d.h(1, "Case details");
    for (const CaseResult* r : order) {
        if (shown++ >= o.chartCases) break;
        d.h(2, r->name + " (" + r->status + ")");
        d.p("Deck " + r->deck + "; engine " + r->engine + "; " + std::to_string(r->steps) + " time steps, " + std::to_string(r->newton) + " Newton iterations, " + fmtv(r->seconds, 3) + " s.");
        report::Table v; v.header = {"Vector", "Reference final", "ResSim final", "rel. RMS (%)", "Pass"};
        for (size_t i = 0; i < r->vectors.size() && i < 12; ++i) {
            const auto& x = r->vectors[i];
            v.rows.push_back({x.key, fmtv(x.refFinal, 6), fmtv(x.cmpFinal, 6), fmtv(100 * x.relRms, 3), x.pass ? "yes" : "NO"});
            if (!x.pass) v.highlightRows.push_back(int(i));
        }
        d.table(v);
        // overlay chart of the (up to) two worst field vectors and the worst well vector
        std::vector<std::string> keys;
        for (const auto& x : r->vectors) if (x.key.find(':') == std::string::npos && keys.size() < 2) keys.push_back(x.key);
        for (const auto& x : r->vectors) if (x.key.find(':') != std::string::npos) { keys.push_back(x.key); break; }
        for (const auto& k : keys) {
            report::Chart ch; ch.title = r->name + ": " + k; ch.xLabel = "Time (days)"; ch.yLabel = k;
            const int a = r->ref.find(k), b = r->run.find(k);
            if (a < 0 || b < 0) continue;
            report::Series s1; s1.name = "reference"; s1.x = r->ref.time; s1.y = r->ref.series(a); s1.rgb = 0x404040; s1.dashed = true;
            report::Series s2; s2.name = "ResSim"; s2.x = r->run.time; s2.y = r->run.series(b); s2.rgb = 0xD62728;
            ch.series = {s1, s2};
            d.chart(ch);
        }
    }
    return d;
}

bool writeJUnit(const std::string& path, const std::vector<CaseResult>& R, std::string& err) {
    std::ofstream f(path);
    if (!f) { err = "cannot write " + path; return false; }
    int nf = 0, ne = 0, ns = 0;
    double tt = 0;
    for (const auto& r : R) { nf += r.status == "FAIL"; ne += r.status == "ERROR"; ns += r.status == "SKIP"; tt += r.seconds; }
    f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<testsuite name=\"ressim-deck-regression\" tests=\"" << R.size() << "\" failures=\"" << nf << "\" errors=\"" << ne
      << "\" skipped=\"" << ns << "\" time=\"" << tt << "\">\n";
    for (const auto& r : R) {
        f << "  <testcase classname=\"" << xesc(r.group.empty() ? "decks" : r.group) << "\" name=\"" << xesc(r.name) << "\" time=\"" << r.seconds << "\">";
        if (r.status == "FAIL") f << "<failure message=\"" << xesc(r.message) << "\">worst " << xesc(r.worstKey) << " " << 100 * r.worst << " %</failure>";
        else if (r.status == "ERROR") f << "<error message=\"" << xesc(r.message) << "\"/>";
        else if (r.status == "SKIP") f << "<skipped message=\"" << xesc(r.message) << "\"/>";
        f << "</testcase>\n";
    }
    f << "</testsuite>\n";
    return bool(f);
}

std::map<std::string, double> readHistory(const std::string& path) {
    std::map<std::string, double> m;
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        std::istringstream is(line);
        std::string date, name, status;
        double sec = 0;
        if (std::getline(is, date, '\t') && std::getline(is, name, '\t') && std::getline(is, status, '\t') && (is >> sec))
            if (status == "PASS" || status == "FAIL") m[name] = sec;   // the latest run that completed
    }
    return m;
}

bool appendHistory(const std::string& path, const std::vector<CaseResult>& R) {
    std::ofstream f(path, std::ios::app);
    if (!f) return false;
    const std::string d = nowStamp();
    for (const auto& r : R) f << d << '\t' << r.name << '\t' << r.status << '\t' << r.seconds << '\t' << r.worst << '\n';
    return bool(f);
}

int writeOutputs(const std::vector<CaseResult>& R, const Options& o, std::string& err) {
    { std::error_code ec; std::filesystem::create_directories(o.outDir, ec); }
    const report::Doc d = reportDoc(R, o);
    if (!report::writeHtml(o.outDir + "/report.html", d, err)) return -1;
    if (!report::writeDocx(o.outDir + "/report.docx", d, err)) return -1;
    if (!writeJUnit(o.outDir + "/junit.xml", R, err)) return -1;
    appendHistory(o.history.empty() ? o.outDir + "/history.tsv" : o.history, R);
    int bad = 0;
    for (const auto& r : R) bad += r.status == "FAIL" || r.status == "ERROR";
    return bad;
}

}  // namespace regress
}  // namespace ressim
