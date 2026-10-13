#include "ecl_engine.h"
#include "opm_engine.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>

#ifdef _WIN32
#include <windows.h>
#endif

namespace fs = std::filesystem;

namespace {
bool isFile(const fs::path& p) { std::error_code ec; return fs::is_regular_file(p, ec); }
bool isDir(const fs::path& p) { std::error_code ec; return fs::is_directory(p, ec); }
fs::path u8(const std::string& s) { return fs::u8path(s); }
std::string lower(std::string s) { std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return char(std::tolower(c)); }); return s; }
std::string envStr(const char* name) {
#ifdef _WIN32
    wchar_t buf[32768];
    std::wstring wn(name, name + std::strlen(name));
    DWORD n = GetEnvironmentVariableW(wn.c_str(), buf, 32768);
    if (n == 0 || n >= 32768) return {};
    return fs::path(std::wstring(buf, n)).u8string();
#else
    const char* v = std::getenv(name); return v ? v : "";
#endif
}

#ifdef _WIN32
const char* kLauncher[] = {"eclrun.exe"};
const char* kSimExe[] = {"eclipse.exe", "e300.exe", "frontsim.exe"};
const char* kPlatforms[] = {"pc_x86_64", "win_x86_64", "pc"};
const char kPathSep = ';';
#else
const char* kLauncher[] = {"eclrun"};
const char* kSimExe[] = {"eclipse", "e300", "frontsim"};
const char* kPlatforms[] = {"linux_x86_64", "linux"};
const char kPathSep = ':';
#endif

// "2023.1" > "2019.2" > "2010.1"; non-numeric names sort last
bool versionNewer(const std::string& a, const std::string& b) {
    double va = std::atof(a.c_str()), vb = std::atof(b.c_str());
    if (va != vb) return va > vb;
    return a > b;
}
bool looksLikeVersion(const std::string& s) { return s.size() >= 4 && std::isdigit((unsigned char)s[0]) && s.find('.') != std::string::npos; }

// launcher or simulator inside an ECLIPSE root (C:\ecl): macros\eclrun.exe first, else the newest
// <version>\bin\<platform>\eclipse.exe
std::string findInRoot(const fs::path& d) {
    if (!isDir(d)) return {};
    for (const char* n : kLauncher) {
        if (isFile(d / n)) return (d / n).u8string();
        if (isFile(d / "macros" / n)) return (d / "macros" / n).u8string();
    }
    for (const char* n : kSimExe) {
        if (isFile(d / n)) return (d / n).u8string();
        for (const char* p : kPlatforms) if (isFile(d / "bin" / p / n)) return (d / "bin" / p / n).u8string();
    }
    std::vector<std::string> vers;
    std::error_code ec;
    for (auto& e : fs::directory_iterator(d, ec)) {
        if (!e.is_directory(ec)) continue;
        const std::string nm = e.path().filename().u8string();
        if (looksLikeVersion(nm)) vers.push_back(nm);
    }
    std::sort(vers.begin(), vers.end(), versionNewer);
    for (auto& v : vers)
        for (const char* p : kPlatforms) {
            const fs::path f = d / v / "bin" / p / kSimExe[0];
            if (isFile(f)) return f.u8string();
        }
    return {};
}

// one or two levels below Program Files\Schlumberger: ECLIPSE / Simulation Launcher / "ECLIPSE 2023.1" style folders
std::string findBelow(const fs::path& d, int depth) {
    std::string r = findInRoot(d);
    if (!r.empty() || depth <= 0 || !isDir(d)) return r;
    std::error_code ec;
    std::vector<fs::path> subs;
    for (auto& e : fs::directory_iterator(d, ec)) {
        if (!e.is_directory(ec)) continue;
        const std::string nm = lower(e.path().filename().u8string());
        if (nm.find("ecl") != std::string::npos || nm.find("simulat") != std::string::npos || nm.find("launcher") != std::string::npos ||
            looksLikeVersion(nm))
            subs.push_back(e.path());
    }
    std::sort(subs.begin(), subs.end(), [](const fs::path& a, const fs::path& b) { return a.filename().u8string() > b.filename().u8string(); });
    for (auto& s : subs) { r = findBelow(s, depth - 1); if (!r.empty()) return r; }
    return {};
}
} // namespace

bool isEclrunLauncher(const std::string& exe) {
    const std::string n = lower(u8(exe).filename().u8string());
    return n.rfind("eclrun", 0) == 0;
}

std::string detectEclipseExecutable(const std::string& hint, std::vector<std::string>* searched) {
    auto note = [&](const std::string& s) { if (searched) searched->push_back(s); };
    if (!hint.empty()) {
        note(hint);
        if (isFile(u8(hint))) return hint;
        std::string r = findBelow(u8(hint), 2);
        if (!r.empty()) return r;
    }
    for (const char* ev : {"ECLIPSE_EXE", "ECLRUN", "ECLPATH"}) {
        const std::string v = envStr(ev);
        if (v.empty()) continue;
        note(std::string("%") + ev + "%");
        if (isFile(u8(v))) return v;
        std::string r = findBelow(u8(v), 1);
        if (!r.empty()) return r;
    }
    note("PATH");
    {
        const std::string path = envStr("PATH");
        std::string sim;
        size_t a = 0;
        while (a <= path.size()) {
            size_t b = path.find(kPathSep, a);
            if (b == std::string::npos) b = path.size();
            std::string dir = path.substr(a, b - a);
            if (dir.size() >= 2 && dir.front() == '"' && dir.back() == '"') dir = dir.substr(1, dir.size() - 2);
            if (!dir.empty()) {
                for (const char* n : kLauncher) if (isFile(u8(dir) / n)) return (u8(dir) / n).u8string();
                if (sim.empty()) for (const char* n : kSimExe) if (isFile(u8(dir) / n)) { sim = (u8(dir) / n).u8string(); break; }
            }
            a = b + 1;
        }
        if (!sim.empty()) return sim;
    }
    std::vector<fs::path> roots;
#ifdef _WIN32
    for (const char* drv : {"C:\\", "D:\\", "E:\\", "F:\\"}) roots.push_back(fs::path(drv) / "ecl");
    for (const char* ev : {"ProgramW6432", "ProgramFiles", "ProgramFiles(x86)"}) {
        const std::string v = envStr(ev);
        if (v.empty()) continue;
        roots.push_back(u8(v) / "Schlumberger");
        roots.push_back(u8(v) / "SLB");
        roots.push_back(u8(v) / "ecl");
    }
    roots.push_back("C:\\Schlumberger");
    roots.push_back("C:\\SLB");
#else
    roots = {"/ecl", "/opt/ecl", "/usr/local/ecl", "/opt/slb", "/opt/schlumberger"};
#endif
    for (auto& r : roots) {
        if (!isDir(r)) continue;
        note(r.u8string());
        std::string f = findBelow(r, 2);
        if (!f.empty()) return f;
    }
    return {};
}

std::vector<std::string> eclipseVersions(const std::string& exe, const std::string& program) {
    std::vector<std::string> out;
    if (exe.empty() || !isFile(u8(exe))) return out;
    if (isEclrunLauncher(exe)) {
        ChildProcess p;
        std::string err, text;
        std::error_code ec;
        if (p.start(exe, {"--report-versions", program}, fs::temp_directory_path(ec).u8string(), &err)) {
            p.pump([&](const std::string& l) { if (text.size() < 4000) text += l + " "; });
            std::istringstream ss(text);
            std::string v;
            while (ss >> v) if (looksLikeVersion(v)) out.push_back(v);
        }
    } else {
        // <root>\<version>\bin\<platform>\eclipse.exe
        fs::path p = u8(exe).parent_path();
        for (int up = 0; up < 3 && !p.empty(); ++up, p = p.parent_path()) {
            const std::string nm = p.filename().u8string();
            if (looksLikeVersion(nm)) { out.push_back(nm); break; }
        }
    }
    std::sort(out.begin(), out.end(), versionNewer);
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

// ------------------------------------------------------------------------------------------------ progress parsing
bool EclProgress::parse(const std::string& line) {
    const char* p = line.c_str();
    while (*p == ' ' || *p == '\t') ++p;
    // continuation lines of an @--ERROR message start with '@' as well
    if (collecting_ && *p == '@' && std::strncmp(p, "@--", 3) != 0) {
        if (++collecting_ > 4) return false;     // first 3 lines of each message (the rest is boilerplate)
        std::string t = p + 1;
        while (!t.empty() && t[0] == ' ') t.erase(0, 1);
        while (!t.empty() && t.back() == ' ') t.pop_back();
        if (!t.empty() && lastError.size() < 600) lastError += (lastError.empty() || lastError.back() == ':' ? " " : " ") + t;
        return false;
    }
    collecting_ = 0;
    if (std::strncmp(p, "@--", 3) == 0) {
        // "@--ERROR", "@--  ERROR" (2023 releases pad the message class), "@--BUG", "@--PROBLEM"
        const char* k = p + 3;
        while (*k == ' ') ++k;
        const bool err = std::strncmp(k, "ERROR", 5) == 0, bug = std::strncmp(k, "BUG", 3) == 0;
        if (std::strncmp(k, "PROBLEM", 7) == 0) ++problems;
        if (err || bug) {
            err ? ++errors : ++bugs;
            if (errors + bugs <= 4 && lastError.size() < 600) {
                if (!lastError.empty()) lastError += " | ";
                lastError += err ? "ERROR:" : "BUG:";
                collecting_ = 1;
            }
            return true;
        }
        return false;
    }
    if (const char* t = std::strstr(p, "TIME=")) {
        // E100: " STEP    5 TIME=    100.00  DAYS (  +30.0  DAYS REPT  3 ITS) (1-APR-2015)"
        double d = 0;
        if (std::sscanf(t, "TIME= %lf", &d) == 1) {
            if (std::strncmp(p, "STEP", 4) == 0) {
                ++timesteps;
                if (const char* its = std::strstr(p, " ITS)")) {
                    const char* q = its;
                    while (q > p && q[-1] == ' ') --q;
                    while (q > p && std::isdigit((unsigned char)q[-1])) --q;
                    newton += std::atoi(q);
                }
            }
            if (d >= day) day = d;
            return true;
        }
    }
    int v = 0;
    if (std::sscanf(p, "Errors %d", &v) == 1) { reportedErrors = v; return true; }
    const std::string l = lower(line);
    if (l.find("licen") != std::string::npos &&
        (l.find("fail") != std::string::npos || l.find("denied") != std::string::npos || l.find("cannot") != std::string::npos ||
         l.find("unable") != std::string::npos || l.find("not available") != std::string::npos || l.find("expired") != std::string::npos ||
         l.find("error") != std::string::npos)) {
        licenseProblem = true;
        std::string t = p;
        if (lastError.find(t) == std::string::npos && lastError.size() < 600) lastError += (lastError.empty() ? "" : " | ") + t;
        return true;
    }
    return false;
}

std::string eclipsePrtErrors(const std::string& prtPath, int maxMessages) {
    std::ifstream f(u8(prtPath));
    if (!f) return {};
    EclProgress p;
    std::string line;
    int n = 0;
    while (std::getline(f, line) && n++ < 2000000) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        p.parse(line);
        if (p.errors + p.bugs > maxMessages) break;
    }
    return p.lastError;
}
