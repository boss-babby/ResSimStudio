// C46 geomodel project file (.rgeo) - see geo_project.h
//
// Format (text, one record per line; names are the rest of the line, numbers in %.17g, "nan" for undefined values):
//   RESGEO-PROJECT 1
//   UNITS METRIC|FIELD
//   SET <key> <value...>
//   CUBE <path>
//   FOLDER <kind> <name>                     kind = well / surface / polygon / fault / line
//   MEMBER <kind> <folder>\t<object>
//   SURFACE <name>      / nx ny x0 y0 dx dy rot / nx*ny values
//   WELL <name>         / x y kb trajSrc td nStations nTops nLogs / stations (md inc azi x y tvd) /
//                         TOP <md> <surface> ... / LOG <file> / null nMd nCurves / md values /
//                         CURVE <discrete> <mnem>\t<unit>\t<descr> / values
//   POLYGON <name>      / closed n / x y pairs
//   FAULT <name>        / throwHint mult nSticks / per stick: n, then x y z triples
//   LINE <name>         / PATH <path> / ns format t0 dt ntr / cdp x y per trace
//   LINEHZ <name>       / nLines / per line: ntr, then TWT values
//   END
#include "ressim/geo_project.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace ressim {
namespace geo {

static const char* kKind[kObjKinds] = {"well", "surface", "polygon", "fault", "line"};
const char* objKindName(ObjKind k) { const int i = int(k); return i >= 0 && i < kObjKinds ? kKind[i] : "?"; }
static int kindOf(const std::string& s) { for (int i = 0; i < kObjKinds; ++i) if (s == kKind[i]) return i; return -1; }

// ---------------------------------------------------------------- folders
std::string ObjFolders::folderOf(ObjKind k, const std::string& obj) const {
    for (const auto& m : member) if (m.first == int(k) && m.second.first == obj) return m.second.second;
    return std::string();
}
void ObjFolders::assign(ObjKind k, const std::string& obj, const std::string& folder) {
    forget(k, obj);
    if (folder.empty()) return;
    addFolder(k, folder);
    member.push_back({int(k), {obj, folder}});
}
bool ObjFolders::addFolder(ObjKind k, const std::string& name) {
    if (name.empty()) return false;
    for (const auto& f : folders) if (f.first == int(k) && f.second == name) return false;
    folders.push_back({int(k), name});
    return true;
}
void ObjFolders::renameFolder(ObjKind k, const std::string& from, const std::string& to) {
    if (to.empty() || from == to) return;
    bool exists = false;
    for (const auto& f : folders) if (f.first == int(k) && f.second == to) exists = true;
    for (auto& f : folders) if (f.first == int(k) && f.second == from) f.second = to;
    if (exists) {   // merged into an existing folder: drop the duplicate entry
        bool seen = false;
        folders.erase(std::remove_if(folders.begin(), folders.end(), [&](const std::pair<int, std::string>& f) {
            if (f.first != int(k) || f.second != to) return false;
            if (!seen) { seen = true; return false; }
            return true;
        }), folders.end());
    }
    for (auto& m : member) if (m.first == int(k) && m.second.second == from) m.second.second = to;
}
void ObjFolders::removeFolder(ObjKind k, const std::string& name) {
    folders.erase(std::remove_if(folders.begin(), folders.end(), [&](const std::pair<int, std::string>& f) { return f.first == int(k) && f.second == name; }), folders.end());
    member.erase(std::remove_if(member.begin(), member.end(), [&](const auto& m) { return m.first == int(k) && m.second.second == name; }), member.end());
}
void ObjFolders::renameObject(ObjKind k, const std::string& from, const std::string& to) {
    for (auto& m : member) if (m.first == int(k) && m.second.first == from) m.second.first = to;
}
void ObjFolders::forget(ObjKind k, const std::string& obj) {
    member.erase(std::remove_if(member.begin(), member.end(), [&](const auto& m) { return m.first == int(k) && m.second.first == obj; }), member.end());
}
std::vector<std::string> ObjFolders::foldersOf(ObjKind k) const {
    std::vector<std::string> out;
    for (const auto& f : folders) if (f.first == int(k)) out.push_back(f.second);
    return out;
}

// ---------------------------------------------------------------- writer
namespace {
std::string clean(const std::string& s) {   // names / paths on one line
    std::string o = s;
    for (char& c : o) if (c == '\n' || c == '\r' || c == '\t') c = ' ';
    return o;
}
struct W {
    std::ostringstream o;
    int col = 0;
    void num(double v) {
        char b[40];
        if (!std::isfinite(v)) std::snprintf(b, sizeof b, "nan");
        else std::snprintf(b, sizeof b, "%.17g", v);
        o << (col ? " " : "") << b;
        if (++col == 8) { o << '\n'; col = 0; }
    }
    void endNums() { if (col) { o << '\n'; col = 0; } }
    template <class T> void nums(const std::vector<T>& v) { for (const T& x : v) num(double(x)); endNums(); }
};
}  // namespace

std::string projectToText(const Project& p) {
    W w;
    auto& o = w.o;
    o << "RESGEO-PROJECT 1\n";
    o << "UNITS " << (p.field ? "FIELD" : "METRIC") << "\n";
    for (const auto& s : p.settings) {
        std::string k = clean(s.first.empty() ? "-" : s.first);
        for (char& c : k) if (c == ' ') c = '_';
        o << "SET " << k << ' ' << clean(s.second) << "\n";
    }
    if (!p.cubePath.empty()) o << "CUBE " << clean(p.cubePath) << "\n";
    for (const auto& f : p.folders.folders) o << "FOLDER " << objKindName(ObjKind(f.first)) << ' ' << clean(f.second) << "\n";
    for (const auto& m : p.folders.member) o << "MEMBER " << objKindName(ObjKind(m.first)) << ' ' << clean(m.second.second) << '\t' << clean(m.second.first) << "\n";
    for (const auto& s : p.surfaces) {
        o << "SURFACE " << clean(s.name) << "\n";
        w.num(s.nx); w.num(s.ny); w.num(s.x0); w.num(s.y0); w.num(s.dx); w.num(s.dy); w.num(s.rotDeg); w.endNums();
        w.nums(s.z);
    }
    for (const auto& g : p.wells) {
        o << "WELL " << clean(g.name) << "\n";
        w.num(g.x); w.num(g.y); w.num(g.kb); w.num(g.trajSrc); w.num(g.td); w.num(double(g.traj.size())); w.num(double(g.tops.size())); w.num(double(g.logs.size())); w.endNums();
        for (const auto& st : g.traj) { w.num(st.md); w.num(st.inc); w.num(st.azi); w.num(st.x); w.num(st.y); w.num(st.tvd); w.endNums(); }
        for (const auto& t : g.tops) { char b[40]; std::snprintf(b, sizeof b, "%.17g", t.md); o << "TOP " << b << ' ' << clean(t.surface) << "\n"; }
        for (const auto& l : g.logs) {
            o << "LOG " << clean(l.file) << "\n";
            w.num(l.nullValue); w.num(double(l.md.size())); w.num(double(l.curves.size())); w.endNums();
            w.nums(l.md);
            for (const auto& c : l.curves) {
                o << "CURVE " << (c.discrete ? 1 : 0) << ' ' << clean(c.mnem) << '\t' << clean(c.unit) << '\t' << clean(c.descr) << "\n";
                w.num(double(c.v.size())); w.endNums();
                w.nums(c.v);
            }
        }
    }
    for (const auto& pg : p.polygons) {
        o << "POLYGON " << clean(pg.name) << "\n";
        w.num(pg.closed ? 1 : 0); w.num(double(pg.pts.size())); w.endNums();
        for (const auto& q : pg.pts) { w.num(q[0]); w.num(q[1]); }
        w.endNums();
    }
    for (const auto& f : p.faults) {
        o << "FAULT " << clean(f.name) << "\n";
        w.num(f.throwHint); w.num(f.mult); w.num(double(f.sticks.size())); w.endNums();
        for (const auto& s : f.sticks) {
            w.num(double(s.size())); w.endNums();
            for (const auto& q : s) { w.num(q[0]); w.num(q[1]); w.num(q[2]); }
            w.endNums();
        }
    }
    for (const auto& l : p.lines) {
        o << "LINE " << clean(l.name) << "\n";
        o << "PATH " << clean(l.path) << "\n";
        w.num(l.ns); w.num(l.format); w.num(l.t0); w.num(l.dt); w.num(double(l.x.size())); w.endNums();
        for (size_t t = 0; t < l.x.size(); ++t) { w.num(t < l.cdp.size() ? l.cdp[t] : int(t + 1)); w.num(l.x[t]); w.num(t < l.y.size() ? l.y[t] : 0.0); }
        w.endNums();
    }
    for (const auto& h : p.lineHorizons) {
        o << "LINEHZ " << clean(h.name) << "\n";
        w.num(double(h.t.size())); w.endNums();
        for (const auto& v : h.t) { w.num(double(v.size())); w.endNums(); w.nums(v); }
    }
    o << "END\n";
    return o.str();
}

bool saveProject(const std::string& path, const Project& p, std::string* err) {
    const std::string tmp = path + ".tmp";
    {
        std::ofstream f(tmp, std::ios::binary);
        if (!f) { if (err) *err = "cannot write " + path; return false; }
        f << projectToText(p);
        if (!f) { if (err) *err = "write error on " + path; return false; }
    }
    std::remove(path.c_str());   // (Windows rename does not replace)
    if (std::rename(tmp.c_str(), path.c_str()) != 0) { if (err) *err = "cannot replace " + path; return false; }
    return true;
}

// ---------------------------------------------------------------- reader
namespace {
struct R {
    std::vector<std::string> lines;
    size_t li = 0;     // next line
    std::istringstream cur;
    bool fail = false;
    std::string why;
    bool next(std::string& key, std::string& rest) {   // next record line (skips number lines are not expected here)
        while (li < lines.size()) {
            std::string s = lines[li++];
            if (!s.empty() && s.back() == '\r') s.pop_back();
            if (s.empty() || s[0] == '#') continue;
            const size_t sp = s.find(' ');
            key = s.substr(0, sp);
            rest = sp == std::string::npos ? std::string() : s.substr(sp + 1);
            return true;
        }
        return false;
    }
    double num() {   // next number, across lines
        std::string tok;
        while (!(cur >> tok)) {
            if (li >= lines.size()) { fail = true; why = "unexpected end of file"; return 0; }
            cur.clear(); cur.str(lines[li++]);
        }
        char* e = nullptr;
        const double v = std::strtod(tok.c_str(), &e);
        if (e == tok.c_str()) { fail = true; why = "number expected, found '" + tok + "' (line " + std::to_string(li) + ")"; return 0; }
        return v;
    }
    size_t count(size_t maxN = 200000000) {
        const double v = num();
        if (!(v >= 0) || v > double(maxN)) { fail = true; why = "bad count (line " + std::to_string(li) + ")"; return 0; }
        return size_t(v);
    }
    void endNums() { cur.clear(); cur.str(std::string()); }   // numbers of a block end at a record line
};
}  // namespace

bool projectFromText(const std::string& text, Project& p, std::string& err) {
    p = Project{};
    R r;
    { std::istringstream is(text); std::string s; while (std::getline(is, s)) r.lines.push_back(s); }
    std::string key, rest;
    if (!r.next(key, rest) || key != "RESGEO-PROJECT") { err = "not a ResGeo project file"; return false; }
    if (std::atoi(rest.c_str()) > 1) { err = "project file version " + rest + " is newer than this program (1)"; return false; }
    bool ended = false;
    while (!r.fail && r.next(key, rest)) {
        if (key == "END") { ended = true; break; }
        else if (key == "UNITS") p.field = rest == "FIELD";
        else if (key == "SET") {
            const size_t sp = rest.find(' ');
            p.settings.push_back({rest.substr(0, sp), sp == std::string::npos ? std::string() : rest.substr(sp + 1)});
        } else if (key == "CUBE") p.cubePath = rest;
        else if (key == "FOLDER") {
            const size_t sp = rest.find(' ');
            const int k = kindOf(rest.substr(0, sp));
            if (k >= 0 && sp != std::string::npos) p.folders.addFolder(ObjKind(k), rest.substr(sp + 1));
        } else if (key == "MEMBER") {
            const size_t sp = rest.find(' '), tb = rest.find('\t');
            const int k = kindOf(rest.substr(0, sp));
            if (k >= 0 && sp != std::string::npos && tb != std::string::npos && tb > sp)
                p.folders.assign(ObjKind(k), rest.substr(tb + 1), rest.substr(sp + 1, tb - sp - 1));
        } else if (key == "SURFACE") {
            Surface s; s.name = rest;
            s.nx = int(r.num()); s.ny = int(r.num()); s.x0 = r.num(); s.y0 = r.num(); s.dx = r.num(); s.dy = r.num(); s.rotDeg = r.num();
            if (s.nx < 0 || s.ny < 0 || double(s.nx) * double(s.ny) > 2e8) { r.fail = true; r.why = "bad surface size"; break; }
            s.z.resize(size_t(s.nx) * size_t(s.ny));
            for (double& z : s.z) z = r.num();
            for (double& z : s.z) if (!std::isfinite(z)) z = std::nan("");
            r.endNums();
            p.surfaces.push_back(std::move(s));
        } else if (key == "WELL") {
            GeoWell g; g.name = rest;
            g.x = r.num(); g.y = r.num(); g.kb = r.num(); g.trajSrc = int(r.num()); g.td = r.num();
            const size_t ns = r.count(), nt = r.count(), nl = r.count();
            g.traj.resize(ns);
            for (auto& st : g.traj) { st.md = r.num(); st.inc = r.num(); st.azi = r.num(); st.x = r.num(); st.y = r.num(); st.tvd = r.num(); }
            r.endNums();
            for (size_t t = 0; t < nt && !r.fail; ++t) {
                if (!r.next(key, rest) || key != "TOP") { r.fail = true; r.why = "TOP expected in well " + g.name; break; }
                const size_t sp = rest.find(' ');
                Top tp; tp.md = std::strtod(rest.c_str(), nullptr); tp.surface = sp == std::string::npos ? std::string() : rest.substr(sp + 1);
                g.tops.push_back(tp);
            }
            for (size_t l = 0; l < nl && !r.fail; ++l) {
                if (!r.next(key, rest) || key != "LOG") { r.fail = true; r.why = "LOG expected in well " + g.name; break; }
                WellLog L; L.well = g.name; L.file = rest;
                L.nullValue = r.num();
                const size_t nmd = r.count(), nc = r.count();
                L.md.resize(nmd);
                for (double& v : L.md) v = r.num();
                r.endNums();
                for (size_t c = 0; c < nc && !r.fail; ++c) {
                    if (!r.next(key, rest) || key != "CURVE") { r.fail = true; r.why = "CURVE expected in well " + g.name; break; }
                    LogCurve C;
                    C.discrete = !rest.empty() && rest[0] == '1';
                    const std::string f = rest.size() > 2 ? rest.substr(2) : std::string();
                    const size_t t1 = f.find('\t'), t2 = t1 == std::string::npos ? std::string::npos : f.find('\t', t1 + 1);
                    C.mnem = f.substr(0, t1);
                    if (t1 != std::string::npos) C.unit = f.substr(t1 + 1, t2 == std::string::npos ? std::string::npos : t2 - t1 - 1);
                    if (t2 != std::string::npos) C.descr = f.substr(t2 + 1);
                    C.v.resize(r.count());
                    for (double& v : C.v) v = r.num();
                    r.endNums();
                    L.curves.push_back(std::move(C));
                }
                g.logs.push_back(std::move(L));
            }
            p.wells.push_back(std::move(g));
        } else if (key == "POLYGON") {
            Polygon pg; pg.name = rest;
            pg.closed = r.num() != 0;
            pg.pts.resize(r.count());
            for (auto& q : pg.pts) { q[0] = r.num(); q[1] = r.num(); }
            r.endNums();
            p.polygons.push_back(std::move(pg));
        } else if (key == "FAULT") {
            Fault f; f.name = rest;
            f.throwHint = r.num(); f.mult = r.num();
            f.sticks.resize(r.count(1000000));
            for (auto& s : f.sticks) { s.resize(r.count()); for (auto& q : s) { q[0] = r.num(); q[1] = r.num(); q[2] = r.num(); } }
            r.endNums();
            p.faults.push_back(std::move(f));
        } else if (key == "LINE") {
            SeismicLine l; l.name = rest;
            if (!r.next(key, rest) || key != "PATH") { r.fail = true; r.why = "PATH expected in line " + l.name; break; }
            l.path = rest;
            l.ns = int(r.num()); l.format = int(r.num()); l.t0 = r.num(); l.dt = r.num();
            const size_t nt = r.count();
            l.cdp.resize(nt); l.x.resize(nt); l.y.resize(nt);
            for (size_t t = 0; t < nt; ++t) { l.cdp[t] = int(r.num()); l.x[t] = r.num(); l.y[t] = r.num(); }
            r.endNums();
            p.lines.push_back(std::move(l));
        } else if (key == "LINEHZ") {
            LineHorizon h; h.name = rest;
            h.t.resize(r.count(100000));
            for (auto& v : h.t) { v.resize(r.count()); for (double& x : v) { x = r.num(); if (!std::isfinite(x)) x = std::nan(""); } }
            r.endNums();
            p.lineHorizons.push_back(std::move(h));
        } else {
            // unknown record of a newer minor version: skipped (numbers that follow it would be read as records and
            // ignored because they are not keywords)
        }
    }
    if (r.fail) { err = r.why; return false; }
    if (!ended) { err = "file is truncated (no END record)"; return false; }
    return true;
}

bool loadProject(const std::string& path, Project& p, std::string& err) {
    std::ifstream f(path, std::ios::binary);
    if (!f) { err = "cannot open " + path; return false; }
    std::stringstream ss; ss << f.rdbuf();
    return projectFromText(ss.str(), p, err);
}

}  // namespace geo
}  // namespace ressim
