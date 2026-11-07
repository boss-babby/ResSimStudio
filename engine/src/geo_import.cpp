// C42 geomodel data import: delimited tables with column assignment, typed builders, well merging, 2-D SEG-Y lines
// (see geo_import.h)
#include "ressim/geo_import.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>

namespace ressim {
namespace geo {

namespace {
const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kPi = 3.14159265358979323846;
std::string trimS(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
    return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
}
std::string upperS(std::string s) { for (char& c : s) c = char(std::toupper(static_cast<unsigned char>(c))); return s; }
std::string stripQuotes(const std::string& s) {
    if (s.size() >= 2 && ((s.front() == '"' && s.back() == '"') || (s.front() == '\'' && s.back() == '\''))) return s.substr(1, s.size() - 2);
    return s;
}
std::vector<std::string> splitCells(const std::string& line, char delim, bool quotes) {
    std::vector<std::string> out; std::string cur; char q = 0; bool any = false;
    for (char c : line) {
        if (quotes && (c == '"' || c == '\'') && (q == 0 || q == c)) {
            if (q == 0 && trimS(cur).empty()) { q = c; any = true; continue; }
            if (q == c) { q = 0; continue; }
        }
        if (!q) {
            if (delim == ' ' ? (c == ' ' || c == '\t') : c == delim) {
                if (delim == ' ') { if (!cur.empty() || any) { out.push_back(cur); cur.clear(); any = false; } }
                else { out.push_back(trimS(cur)); cur.clear(); any = false; }
                continue;
            }
        }
        if (c == '\r' || c == '\n') continue;
        cur += c;
    }
    if (delim == ' ') { if (!cur.empty() || any) out.push_back(cur); }
    else out.push_back(trimS(cur));
    return out;
}
int countOutsideQuotes(const std::string& line, char d) {
    int n = 0; char q = 0;
    for (char c : line) {
        if (c == '"' || c == '\'') { if (q == 0) q = c; else if (q == c) q = 0; continue; }
        if (!q && c == d) ++n;
    }
    return n;
}
bool plainNumber(const std::string& s0, bool decComma, double& v) {
    std::string s = trimS(s0);
    if (s.empty()) return false;
    if (decComma) for (char& c : s) if (c == ',') c = '.';
    char* e = nullptr; v = std::strtod(s.c_str(), &e);
    return e && *e == 0 && std::isfinite(v);
}
bool isCommentLine(const std::string& t, const std::string& lead) {
    if (t.empty()) return false;
    if (t.size() >= 2 && t[0] == '-' && t[1] == '-' && lead.find('-') != std::string::npos) return true;
    return t[0] != '-' && lead.find(t[0]) != std::string::npos;
}
}  // namespace

// ================================================================ tables
TextTable parseTextTable(const std::string& text, const TableFormat& fmt, std::string& err) {
    TextTable T;
    std::vector<std::string> lines;
    {
        std::string cur;
        for (char c : text) { if (c == '\n') { lines.push_back(cur); cur.clear(); } else cur += c; }
        if (!cur.empty()) lines.push_back(cur);
    }
    for (size_t k = 0; k < lines.size() && k < 30; ++k) T.rawHead.push_back(lines[k]);
    // Petrel header block / meta / comments / separators
    std::vector<std::string> body;
    bool inHdr = false, petrelHdr = false, hasPetrel = false;
    for (const auto& l : lines) if (upperS(trimS(l)) == "BEGIN HEADER") hasPetrel = true;
    std::vector<std::string> petrelNames;
    for (size_t k = size_t(std::max(0, fmt.skipLines)); k < lines.size(); ++k) {
        const std::string t = trimS(lines[k]);
        if (t.empty()) continue;
        const std::string U = upperS(t);
        if (U == "BEGIN HEADER") { inHdr = true; petrelHdr = true; continue; }
        if (U == "END HEADER") { inHdr = false; continue; }
        if (inHdr) { petrelNames.push_back(stripQuotes(t)); continue; }
        if (isCommentLine(t, fmt.comment)) {
            const size_t c = t.find(':');
            if (c != std::string::npos) {
                std::string key = trimS(t.substr(t.find_first_not_of("#!-") == std::string::npos ? 0 : t.find_first_not_of("#!-"), std::string::npos));
                const size_t c2 = key.find(':');
                if (c2 != std::string::npos) T.meta[upperS(trimS(key.substr(0, c2)))] = trimS(key.substr(c2 + 1));
            }
            continue;
        }
        if (t.rfind("====", 0) == 0 || t.rfind("----", 0) == 0) continue;
        if (hasPetrel && body.empty() && U.rfind("VERSION", 0) == 0) continue;
        body.push_back(lines[k]);
    }
    if (body.empty()) { err = "no data lines"; return T; }
    // delimiter
    char d = fmt.delim;
    if (d == 0) {
        d = ' ';
        const size_t nS = std::min<size_t>(body.size(), 30);
        for (char c : {'\t', ';', ',', '|'}) {
            std::map<int, int> hist; int zero = 0;
            for (size_t k = 0; k < nS; ++k) { const int n = countOutsideQuotes(body[k], c); if (n == 0) ++zero; ++hist[n]; }
            int mode = 0, best = 0; for (auto& h : hist) if (h.second > best && h.first > 0) { best = h.second; mode = h.first; }
            if (mode > 0 && best >= int(0.8 * double(nS)) && zero <= int(0.2 * double(nS))) { d = c; break; }
        }
    }
    T.delimUsed = d;
    std::vector<std::vector<std::string>> cells;
    for (const auto& l : body) cells.push_back(splitCells(l, d, fmt.quotes));
    for (auto& r : cells) for (auto& c : r) c = trimS(c);
    // decimal comma
    bool dc = false;
    if (d != ',') {
        if (fmt.decimalComma == 1) dc = true;
        else if (fmt.decimalComma < 0) {
            int hits = 0, dots = 0;
            for (size_t k = 0; k < cells.size() && k < 40; ++k)
                for (auto& c : cells[k]) {
                    const size_t p = c.find(',');
                    if (p != std::string::npos && p > 0 && p + 1 < c.size() && std::isdigit((unsigned char)c[p - 1]) && std::isdigit((unsigned char)c[p + 1]) && c.find(',', p + 1) == std::string::npos) ++hits;
                    if (c.find('.') != std::string::npos) { double v; if (plainNumber(c, false, v)) ++dots; }
                }
            dc = hits > 0 && hits >= dots;
        }
    }
    T.decimalCommaUsed = dc;
    auto numericCount = [&](const std::vector<std::string>& r) { int n = 0; double v; for (auto& c : r) if (plainNumber(c, dc, v)) ++n; return n; };
    // header
    size_t first = 0;
    if (petrelHdr && !petrelNames.empty()) { T.header = petrelNames; T.hasHeader = true; }
    else if (fmt.headerRow >= 0) {
        const size_t h = size_t(fmt.headerRow);
        if (h < cells.size()) { T.header = cells[h]; T.hasHeader = true; first = h + 1; }
    } else if (fmt.headerRow == -2 && !cells.empty()) {
        const auto& r0 = cells[0];
        bool hdr = false;
        if (cells.size() == 1) hdr = numericCount(r0) == 0;
        else {
            const auto& r1 = cells[1];
            double v;
            int flips = 0;
            for (size_t c = 0; c < r0.size() && c < r1.size(); ++c) if (!plainNumber(r0[c], dc, v) && plainNumber(r1[c], dc, v)) ++flips;
            hdr = flips > 0 || (numericCount(r0) == 0 && numericCount(r1) == 0 && cells.size() > 2 && numericCount(cells[2]) > 0);
        }
        if (hdr) {
            T.header = r0; T.hasHeader = true; first = 1;
            // units row ("m  m  deg") right below the header
            if (cells.size() > 2 && numericCount(cells[1]) == 0 && numericCount(cells[2]) > 0) {
                std::string u; for (auto& c : cells[1]) u += (u.empty() ? "" : " ") + c;
                T.meta["UNITS"] = u; first = 2;
            }
        }
    }
    int nc = int(T.header.size());
    for (size_t k = first; k < cells.size(); ++k) nc = std::max(nc, int(cells[k].size()));
    T.ncol = nc;
    for (int c = int(T.header.size()); c < nc; ++c) T.header.push_back("C" + std::to_string(c + 1));
    for (auto& h : T.header) h = stripQuotes(trimS(h));
    for (size_t k = first; k < cells.size(); ++k) {
        ++T.totalRows;
        if (fmt.maxRows > 0 && int(T.rows.size()) >= fmt.maxRows) continue;
        auto r = cells[k]; r.resize(size_t(nc));
        T.rows.push_back(std::move(r));
    }
    if (T.totalRows == 0) err = "no data rows below the header";
    return T;
}
TextTable readTextTable(const std::string& path, const TableFormat& fmt, std::string& err) {
    std::ifstream f(path, std::ios::binary);
    if (!f) { err = "cannot open " + path; return TextTable(); }
    std::ostringstream o; o << f.rdbuf();
    std::string e;
    TextTable t = parseTextTable(o.str(), fmt, e);
    if (!e.empty()) err = path + ": " + e;
    return t;
}
bool tableNumber(const TextTable& t, const TableFormat& fmt, const std::string& cell, double& v) {
    const std::string c = trimS(cell);
    for (const auto& n : fmt.nulls) if (c == n) return false;
    if (!plainNumber(c, t.decimalCommaUsed, v)) return false;
    for (const auto& n : fmt.nulls) { double nv; if (!n.empty() && plainNumber(n, false, nv) && std::fabs(nv - v) < 1e-9 * std::max(1.0, std::fabs(nv))) return false; }
    return true;
}

// ================================================================ roles
const char* colRoleName(ColRole r) {
    static const char* n[] = {"(ignore)", "Well name", "X (east)", "Y (north)", "KB / datum elevation", "TD (total depth MD)", "MD",
                              "Inclination", "Azimuth", "TVD (below KB)", "TVDSS (depth, +down)", "Z (elevation, -down)", "DX (east offset)", "DY (north offset)",
                              "Surface / marker", "Stick id", "Fault name", "Value"};
    const int i = int(r);
    return i >= 0 && i < int(ColRole::Count) ? n[i] : "?";
}
const char* importKindName(ImportKind k) {
    static const char* n[] = {"Well heads", "Deviation surveys", "Well tops", "Points (XYZ)", "Fault sticks"};
    const int i = int(k);
    return i >= 0 && i < int(ImportKind::Count) ? n[i] : "?";
}
std::vector<ColRole> rolesFor(ImportKind k) {
    using R = ColRole;
    switch (k) {
        case ImportKind::WellHeads: return {R::Ignore, R::WellName, R::X, R::Y, R::KB, R::TD};
        case ImportKind::Deviation: return {R::Ignore, R::WellName, R::MD, R::Inc, R::Azi, R::X, R::Y, R::TVD, R::TVDSS, R::Z, R::DX, R::DY};
        case ImportKind::Tops: return {R::Ignore, R::WellName, R::Surface, R::MD, R::TVD, R::TVDSS, R::Z};
        case ImportKind::Points: return {R::Ignore, R::X, R::Y, R::Z, R::TVDSS, R::Surface, R::Value};
        case ImportKind::FaultSticks: return {R::Ignore, R::X, R::Y, R::Z, R::TVDSS, R::StickId, R::FaultName};
        default: return {R::Ignore};
    }
}
std::string normColName(const std::string& s0) {
    std::string s, o;
    int depth = 0;
    for (char c : s0) { if (c == '(' || c == '[' || c == '{') { ++depth; continue; } if (c == ')' || c == ']' || c == '}') { depth = std::max(0, depth - 1); continue; } if (!depth) s += c; }
    if (trimS(s).empty()) s = s0;
    for (char c : s) if (std::isalnum((unsigned char)c)) o += char(std::toupper((unsigned char)c));
    return o;
}
namespace {
struct Syn { ColRole r; std::vector<const char*> exact; std::vector<const char*> contains; };
std::vector<Syn> synonyms(ImportKind k) {
    using R = ColRole;
    std::vector<Syn> s = {
        {R::WellName, {"WELL", "WELLNAME", "NAME", "UWI", "WELLID", "BOREHOLE", "WELLBORE", "API", "WELLIDENTIFIER", "BH", "HOLE", "WELLNO"}, {"WELLNAME", "WELLBORE", "BOREHOLE", "UWI"}},
        {R::X, {"X", "EASTING", "EAST", "XCOORD", "XCOORDINATE", "SURFACEX", "WELLHEADX", "HEADX", "XLOC", "UTMX", "XUTM", "E", "XM", "XLOCATION", "TOPHOLEX", "SURFX", "LOCX"}, {"EASTING"}},
        {R::Y, {"Y", "NORTHING", "NORTH", "YCOORD", "YCOORDINATE", "SURFACEY", "WELLHEADY", "HEADY", "YLOC", "UTMY", "YUTM", "N", "YM", "YLOCATION", "TOPHOLEY", "SURFY", "LOCY"}, {"NORTHING"}},
        {R::KB, {"KB", "KBELEV", "KBELEVATION", "ELEV", "ELEVATION", "RKB", "RTE", "DF", "KELLYBUSHING", "DATUM", "WELLDATUM", "ZDATUM", "REFELEV", "REFERENCEELEVATION", "DATUMELEVATION", "DFE", "GL", "WELLDATUMVALUE"}, {"KELLY", "KBELEV", "DATUMELEV"}},
        {R::TD, {"TD", "TOTALDEPTH", "TDMD", "BOTTOMMD", "MAXMD", "TOTALMD", "TDKB"}, {"TOTALDEPTH"}},
        {R::MD, {"MD", "MEASUREDDEPTH", "DEPTH", "DEPT", "MDKB", "MDRKB", "MDM", "MDFT", "MEASDEPTH", "MDEPTH"}, {"MEASURED"}},
        {R::Inc, {"INC", "INCL", "INCLINATION", "DEVI", "DEVIATION", "ANGLE", "INCLIN", "DEV", "INCDEG", "HOLEANGLE"}, {"INCLIN"}},
        {R::Azi, {"AZI", "AZIM", "AZIMUTH", "AZ", "AZIMUTHGN", "AZIMUTHTN", "HAZI", "DIRECTION", "AZIDEG", "AZIMGN", "AZIMTN"}, {"AZIM"}},
        {R::TVD, {"TVD", "TVDKB", "TVDRKB", "TRUEVERTICALDEPTH", "TVDM"}, {}},
        {R::TVDSS, {"TVDSS", "SSTVD", "SUBSEA", "TVDMSL", "SSDEPTH", "SUBSEADEPTH", "DEPTHSS", "TVDSD"}, {"TVDSS", "SUBSEA"}},
        {R::Z, {"Z", "ZM", "ZCOORD", "ELEVATIONZ", "ZVALUE"}, {}},
        {R::DX, {"DX", "DELTAX", "EW", "EASTWEST", "EOFFSET", "EASTOFFSET", "XOFFSET", "DEAST", "EWOFFSET"}, {}},
        {R::DY, {"DY", "DELTAY", "NS", "NORTHSOUTH", "NOFFSET", "NORTHOFFSET", "YOFFSET", "DNORTH", "NSOFFSET"}, {}},
        {R::Surface, {"SURFACE", "HORIZON", "TOP", "MARKER", "FORMATION", "PICK", "SURFACENAME", "ZONE", "HORIZONNAME", "TOPNAME", "MARKERNAME", "SET", "LAYER"}, {"SURFACE", "HORIZON", "MARKER", "FORMATION"}},
        {R::StickId, {"STICK", "STICKID", "STICKNO", "STICKNUMBER", "SEGMENT", "SEGMENTID", "ID"}, {"STICK"}},
        {R::FaultName, {"FAULT", "FAULTNAME", "FAULTID", "NAME"}, {"FAULT"}},
        {R::Value, {"VALUE", "ATTRIBUTE", "V", "PROPERTY", "AMPLITUDE"}, {}},
    };
    if (k == ImportKind::WellHeads) for (auto& x : s) if (x.r == ColRole::KB) { x.exact.push_back("Z"); x.exact.push_back("ZM"); }
    if (k == ImportKind::Points || k == ImportKind::FaultSticks) for (auto& x : s) if (x.r == ColRole::TVDSS) { x.exact.push_back("DEPTH"); x.exact.push_back("TVD"); }
    if (k == ImportKind::Points) for (auto& x : s) if (x.r == ColRole::Surface) { x.exact.push_back("NAME"); x.exact.push_back("SETNAME"); x.exact.push_back("LABEL"); }
    return s;
}
std::vector<char> numericColumns(const TextTable& t, const TableFormat& fmt) {
    std::vector<char> num(size_t(t.ncol), 0);
    for (int c = 0; c < t.ncol; ++c) {
        int n = 0, nn = 0; double v;
        for (size_t k = 0; k < t.rows.size() && k < 60; ++k) {
            const std::string& s = t.rows[k][size_t(c)];
            bool isNull = false; for (auto& x : fmt.nulls) if (trimS(s) == x) isNull = true;
            if (isNull) continue;
            ++n; if (plainNumber(s, t.decimalCommaUsed, v)) ++nn;
        }
        num[size_t(c)] = n > 0 && nn >= int(0.8 * n) ? 1 : 0;
    }
    return num;
}
}  // namespace

ColumnMap guessColumns(const TextTable& t, const TableFormat& fmt, ImportKind k) {
    ColumnMap m; m.role.assign(size_t(t.ncol), ColRole::Ignore); m.scale.assign(size_t(t.ncol), 1.0);
    const auto allowed = rolesFor(k);
    auto ok = [&](ColRole r) { return std::find(allowed.begin(), allowed.end(), r) != allowed.end(); };
    const auto syn = synonyms(k);
    const auto num = numericColumns(t, fmt);
    std::set<ColRole> used;
    auto assign = [&](int c, ColRole r) {
        if (m.role[size_t(c)] != ColRole::Ignore || used.count(r) || !ok(r)) return;
        const bool wantText = r == ColRole::WellName || r == ColRole::Surface || r == ColRole::FaultName;
        if (wantText == bool(num[size_t(c)]) && !t.rows.empty() && r != ColRole::StickId && r != ColRole::WellName) return;   // type clash (well names may be numeric)
        m.role[size_t(c)] = r; used.insert(r);
    };
    if (t.hasHeader) {
        for (int pass = 0; pass < 2; ++pass)
            for (int c = 0; c < t.ncol; ++c) {
                const std::string n = normColName(t.header[size_t(c)]);
                if (n.empty()) continue;
                for (const auto& s : syn) {
                    bool hit = false;
                    if (pass == 0) { for (auto e : s.exact) if (n == e) hit = true; }
                    else for (auto e : s.contains) if (n.find(e) != std::string::npos) hit = true;
                    if (hit) { assign(c, s.r); if (m.role[size_t(c)] != ColRole::Ignore) break; }
                }
            }
    }
    // positional fallback for missing required roles (headerless files, or names nobody can match)
    auto firstFree = [&](bool numeric) { for (int c = 0; c < t.ncol; ++c) if (m.role[size_t(c)] == ColRole::Ignore && bool(num[size_t(c)]) == numeric) return c; return -1; };
    auto need = [&](ColRole r, bool numeric) { if (used.count(r)) return; const int c = firstFree(numeric); if (c >= 0) { m.role[size_t(c)] = r; used.insert(r); } };
    const bool fallback = !t.hasHeader || !validateMap(m, k).empty();
    if (fallback) {
        using R = ColRole;
        switch (k) {
            case ImportKind::WellHeads: need(R::WellName, false); need(R::X, true); need(R::Y, true); if (!t.hasHeader) { need(R::KB, true); need(R::TD, true); } break;
            case ImportKind::Deviation:
                if (!used.count(R::X) && !used.count(R::DX)) { if (!t.hasHeader) need(R::WellName, false); need(R::MD, true); need(R::Inc, true); need(R::Azi, true); }
                break;
            case ImportKind::Tops: need(R::WellName, false); need(R::Surface, false); if (!used.count(R::TVDSS) && !used.count(R::Z) && !used.count(R::TVD)) need(R::MD, true); break;
            case ImportKind::Points: need(R::X, true); need(R::Y, true); if (!used.count(R::TVDSS)) need(R::Z, true); if (!t.hasHeader) need(R::Surface, false); break;
            case ImportKind::FaultSticks: need(R::X, true); need(R::Y, true); if (!used.count(R::TVDSS)) need(R::Z, true); if (!t.hasHeader) { need(R::StickId, true); need(R::FaultName, false); } break;
            default: break;
        }
    }
    return m;
}
std::string validateMap(const ColumnMap& m, ImportKind k) {
    using R = ColRole;
    auto has = [&](R r) { return m.find(r) >= 0; };
    std::set<R> seen;
    for (auto r : m.role) { if (r == R::Ignore) continue; if (seen.count(r)) return std::string("column role '") + colRoleName(r) + "' is assigned twice"; seen.insert(r); }
    const bool depth = has(R::TVD) || has(R::TVDSS) || has(R::Z);
    switch (k) {
        case ImportKind::WellHeads: if (!has(R::WellName)) return "assign a Well name column"; if (!has(R::X) || !has(R::Y)) return "assign the X and Y columns"; return "";
        case ImportKind::Deviation:
            if (has(R::Inc) || has(R::Azi)) { if (!(has(R::MD) && has(R::Inc) && has(R::Azi))) return "MD / Inclination / Azimuth surveys need all three columns"; return ""; }
            if (has(R::X) || has(R::Y)) { if (!(has(R::X) && has(R::Y) && depth)) return "X / Y surveys need X, Y and a depth column (TVD, TVDSS or Z)"; return ""; }
            if (has(R::DX) || has(R::DY)) { if (!(has(R::DX) && has(R::DY) && depth)) return "offset surveys need DX, DY and a depth column (TVD, TVDSS or Z)"; return ""; }
            return "assign MD + Inclination + Azimuth, or X + Y + depth, or DX + DY + depth";
        case ImportKind::Tops: if (!has(R::Surface)) return "assign a Surface / marker column"; if (!has(R::MD) && !depth) return "assign an MD or a depth (TVD / TVDSS / Z) column"; return "";
        case ImportKind::Points: case ImportKind::FaultSticks: if (!has(R::X) || !has(R::Y)) return "assign the X and Y columns"; if (!has(R::Z) && !has(R::TVDSS)) return "assign a Z or TVDSS column"; return "";
        default: return "unknown import kind";
    }
}

// ================================================================ builders
ImportResult buildImport(const TextTable& t, const TableFormat& fmt, const ColumnMap& m, ImportKind k, const ImportOptions& o) {
    using R = ColRole;
    ImportResult res;
    const std::string bad = validateMap(m, k);
    if (!bad.empty()) { res.warnings.push_back(bad); return res; }
    auto num = [&](const std::vector<std::string>& row, R r) {
        const int c = m.find(r); if (c < 0 || size_t(c) >= row.size()) return kNaN;
        double v; if (!tableNumber(t, fmt, row[size_t(c)], v)) return kNaN;
        return v * (size_t(c) < m.scale.size() ? m.scale[size_t(c)] : 1.0);
    };
    auto txt = [&](const std::vector<std::string>& row, R r) { const int c = m.find(r); return c < 0 || size_t(c) >= row.size() ? std::string() : stripQuotes(trimS(row[size_t(c)])); };
    // Z sign: auto = elevation when most values are negative
    double zSign = 1;
    if (m.find(R::Z) >= 0) {
        if (o.zSign == 1) zSign = -1;
        else if (o.zSign == 2) zSign = 1;
        else { int neg = 0, n = 0; for (auto& r : t.rows) { const double z = num(r, R::Z); if (std::isfinite(z)) { ++n; if (z < 0) ++neg; } } zSign = neg * 2 > n ? -1 : 1; }
    }
    auto depthOf = [&](const std::vector<std::string>& row, bool& belowKb) {
        belowKb = false;
        double d = num(row, R::TVDSS); if (std::isfinite(d)) return d;
        d = num(row, R::Z); if (std::isfinite(d)) return zSign * d;
        d = num(row, R::TVD); if (std::isfinite(d)) { belowKb = true; return d; }
        return kNaN;
    };
    if (k == ImportKind::WellHeads) {
        std::map<std::string, size_t> at;
        for (const auto& row : t.rows) {
            const std::string nm = txt(row, R::WellName); const double x = num(row, R::X), y = num(row, R::Y);
            if (nm.empty() || !std::isfinite(x) || !std::isfinite(y)) { ++res.rowsSkipped; continue; }
            GeoWell w; w.name = nm; w.x = x; w.y = y;
            const double kb = num(row, R::KB); if (std::isfinite(kb)) w.kb = kb;
            const double td = num(row, R::TD); if (std::isfinite(td)) w.td = td;
            const std::string key = wellKey(nm);
            if (at.count(key)) { res.warnings.push_back("duplicate well '" + nm + "' - the last row wins"); res.wells[at[key]] = w; }
            else { at[key] = res.wells.size(); res.wells.push_back(w); }
            ++res.rowsUsed;
        }
        res.td.clear(); for (auto& w : res.wells) res.td.push_back(w.td > 0 ? w.td : kNaN);
    } else if (k == ImportKind::Deviation) {
        const bool angles = m.find(R::Inc) >= 0, absXY = m.find(R::X) >= 0;
        res.trajMode = angles ? 0 : absXY ? 1 : 2;
        std::map<std::string, size_t> at;
        bool kbFlagSet = false;
        for (const auto& row : t.rows) {
            std::string nm = txt(row, R::WellName);
            if (nm.empty()) nm = o.defaultWell.empty() ? std::string("WELL") : o.defaultWell;
            Station s; bool ok = true;
            s.md = num(row, R::MD);
            if (res.trajMode == 0) { s.inc = num(row, R::Inc); s.azi = num(row, R::Azi); ok = std::isfinite(s.md) && std::isfinite(s.inc) && std::isfinite(s.azi); }
            else {
                bool kb = false; s.tvd = depthOf(row, kb);
                if (!kbFlagSet) { res.tvdBelowKb = kb; kbFlagSet = true; }
                s.x = num(row, res.trajMode == 1 ? R::X : R::DX); s.y = num(row, res.trajMode == 1 ? R::Y : R::DY);
                ok = std::isfinite(s.tvd) && std::isfinite(s.x) && std::isfinite(s.y);
                if (!std::isfinite(s.md)) s.md = -1;   // from arc length on merge
            }
            if (!ok) { ++res.rowsSkipped; continue; }
            const std::string key = wellKey(nm);
            if (!at.count(key)) { at[key] = res.wells.size(); GeoWell w; w.name = nm; w.trajSrc = res.trajMode == 1 ? 1 : 0; res.wells.push_back(w); }
            res.wells[at[key]].traj.push_back(s);
            ++res.rowsUsed;
        }
        for (auto& w : res.wells) {
            if (res.trajMode == 0) {
                std::stable_sort(w.traj.begin(), w.traj.end(), [](const Station& a, const Station& b) { return a.md < b.md; });
                for (auto& s : w.traj) { if (s.inc < 0 || s.inc > 180) res.warnings.push_back(w.name + ": inclination outside 0..180 deg"); s.azi = std::fmod(std::fmod(s.azi, 360.0) + 360.0, 360.0); }
            }
        }
    } else if (k == ImportKind::Tops) {
        for (const auto& row : t.rows) {
            TopRow r; r.well = txt(row, R::WellName); if (r.well.empty()) r.well = o.defaultWell;
            r.surface = txt(row, R::Surface);
            r.md = num(row, R::MD); r.hasMd = std::isfinite(r.md);
            bool kb = false; r.tvdss = depthOf(row, kb); r.hasTvd = std::isfinite(r.tvdss); r.tvdBelowKb = kb;
            if (r.well.empty() || r.surface.empty() || (!r.hasMd && !r.hasTvd)) { ++res.rowsSkipped; continue; }
            res.tops.push_back(r); ++res.rowsUsed;
        }
        if (res.tops.empty() && m.find(R::WellName) < 0 && o.defaultWell.empty()) res.warnings.push_back("no Well name column and no target well");
    } else if (k == ImportKind::Points) {
        std::map<std::string, size_t> at;
        for (const auto& row : t.rows) {
            const double x = num(row, R::X), y = num(row, R::Y); bool kb = false; const double z = depthOf(row, kb);
            if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) { ++res.rowsSkipped; continue; }
            std::string nm = txt(row, R::Surface); if (nm.empty()) nm = o.defaultSet.empty() ? std::string("points") : o.defaultSet;
            if (!at.count(nm)) { at[nm] = res.points.size(); res.points.push_back({nm, {}}); }
            res.points[at[nm]].second.push_back({x, y, z});
            ++res.rowsUsed;
        }
    } else if (k == ImportKind::FaultSticks) {
        std::map<std::string, size_t> fi;
        std::map<std::string, std::map<std::string, size_t>> si;   // fault -> stick key -> index
        std::string lastFault; double lastZ = kNaN; int autoStick = 0;
        const bool hasId = m.find(R::StickId) >= 0;
        for (const auto& row : t.rows) {
            const double x = num(row, R::X), y = num(row, R::Y); bool kb = false; const double z = depthOf(row, kb);
            if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) { ++res.rowsSkipped; lastZ = kNaN; ++autoStick; continue; }
            std::string fn = txt(row, R::FaultName); if (fn.empty()) fn = o.defaultSet.empty() ? std::string("fault") : o.defaultSet;
            std::string sk;
            if (hasId) sk = txt(row, R::StickId);
            else { if (fn != lastFault || (std::isfinite(lastZ) && z < lastZ - 1e-9)) ++autoStick; sk = std::to_string(autoStick); }
            lastFault = fn; lastZ = z;
            if (!fi.count(fn)) { fi[fn] = res.faults.size(); Fault f; f.name = fn; res.faults.push_back(f); }
            Fault& F = res.faults[fi[fn]];
            auto& smap = si[fn];
            if (!smap.count(sk)) { smap[sk] = F.sticks.size(); F.sticks.emplace_back(); }
            F.sticks[smap[sk]].push_back({x, y, z});
            ++res.rowsUsed;
        }
        for (auto& f : res.faults) for (auto& s : f.sticks) std::stable_sort(s.begin(), s.end(), [](const P3& a, const P3& b) { return a[2] < b[2]; });
    }
    if (res.rowsUsed == 0) res.warnings.push_back("no usable rows (check the column assignment, the null value and the decimal separator)");
    return res;
}

// ================================================================ wells
std::string wellKey(const std::string& name) { return upperS(trimS(stripQuotes(trimS(name)))); }
int matchWell(const std::vector<GeoWell>& wells, const std::string& nameOrFile) {
    const std::string k = wellKey(nameOrFile);
    for (size_t i = 0; i < wells.size(); ++i) if (wellKey(wells[i].name) == k) return int(i);
    int best = -1; size_t bl = 0;
    for (size_t i = 0; i < wells.size(); ++i) {
        const std::string p = wellKey(wells[i].name);
        if (p.empty() || p.size() >= k.size() || k.compare(0, p.size(), p) != 0) continue;
        const char nx = k[p.size()];
        if (std::isalnum((unsigned char)nx)) continue;   // "A1" must not match "A10_survey"
        if (p.size() > bl) { bl = p.size(); best = int(i); }
    }
    return best;
}
static int exactWell(const std::vector<GeoWell>& wells, const std::string& n) {
    const std::string k = wellKey(n);
    for (size_t i = 0; i < wells.size(); ++i) if (wellKey(wells[i].name) == k) return int(i);
    return -1;
}
static void recompute(GeoWell& w) {
    if (w.traj.empty()) return;
    if (w.trajSrc == 0) minimumCurvature(w);
}
MergeStats mergeHeads(std::vector<GeoWell>& proj, const ImportResult& r) {
    MergeStats st;
    for (const auto& h : r.wells) {
        const int i = exactWell(proj, h.name);
        if (i < 0) { GeoWell w; w.name = h.name; w.x = h.x; w.y = h.y; w.kb = h.kb; w.td = h.td; if (w.td > 0) setVertical(w, w.td); proj.push_back(w); ++st.added; continue; }
        GeoWell& w = proj[size_t(i)];
        const double dx = h.x - w.x, dy = h.y - w.y, dkb = h.kb - w.kb;
        w.x = h.x; w.y = h.y; w.kb = h.kb; if (h.td > 0) w.td = h.td;
        if (w.trajSrc == 0) recompute(w);
        else {
            for (auto& s : w.traj) { s.x += dx; s.y += dy; }
            (void)dkb;   // absolute surveys keep their TVDSS
            if (std::fabs(dx) + std::fabs(dy) > 1e-9) st.warnings.push_back(w.name + ": X/Y survey shifted with the new head");
        }
        ++st.updated;
    }
    return st;
}
void inclinationsFromXyz(GeoWell& w) {
    auto& T = w.traj;
    const size_t n = T.size();
    for (size_t k = 0; k < n; ++k) {
        const size_t a = k ? k - 1 : 0, b = std::min(n - 1, k + 1);
        if (a == b) { T[k].inc = 0; T[k].azi = 0; continue; }
        const double dx = T[b].x - T[a].x, dy = T[b].y - T[a].y, dz = T[b].tvd - T[a].tvd;
        const double h = std::sqrt(dx * dx + dy * dy), L = std::sqrt(h * h + dz * dz);
        if (L < 1e-12) { T[k].inc = k ? T[k - 1].inc : 0; T[k].azi = k ? T[k - 1].azi : 0; continue; }
        T[k].inc = std::atan2(h, dz) * 180.0 / kPi;
        T[k].azi = h > 1e-9 ? std::fmod(std::atan2(dx, dy) * 180.0 / kPi + 360.0, 360.0) : (k ? T[k - 1].azi : 0);
    }
}
bool trajFromXyz(GeoWell& w, std::vector<Station> st, std::string& err) {
    if (st.size() < 2) { err = w.name + ": a survey needs at least two stations"; return false; }
    bool haveMd = true; for (auto& s : st) if (!(s.md >= 0)) haveMd = false;
    if (haveMd) std::stable_sort(st.begin(), st.end(), [](const Station& a, const Station& b) { return a.md < b.md; });
    else {
        double md = st[0].tvd + w.kb;   // vertical from the datum to the first station
        if (md < 0) md = 0;
        st[0].md = md;
        for (size_t k = 1; k < st.size(); ++k) {
            md += std::sqrt(std::pow(st[k].x - st[k - 1].x, 2) + std::pow(st[k].y - st[k - 1].y, 2) + std::pow(st[k].tvd - st[k - 1].tvd, 2));
            st[k].md = md;
        }
    }
    for (size_t k = 1; k < st.size(); ++k) if (st[k].md < st[k - 1].md) { err = w.name + ": MD decreases at station " + std::to_string(k + 1); return false; }
    w.traj = st; w.trajSrc = 1;
    inclinationsFromXyz(w);
    return true;
}
void setVertical(GeoWell& w, double td) {
    w.traj.clear(); w.trajSrc = 0; w.td = td;
    if (!(td > 0)) return;
    Station a; a.md = 0; Station b; b.md = td;
    w.traj = {a, b};
    minimumCurvature(w);
}
MergeStats mergeSurveys(std::vector<GeoWell>& proj, const ImportResult& r) {
    MergeStats st;
    for (const auto& in : r.wells) {
        int i = matchWell(proj, in.name);
        if (i < 0) {
            GeoWell w; w.name = in.name;
            if (r.trajMode == 1 && !in.traj.empty()) { w.x = in.traj[0].x; w.y = in.traj[0].y; if (in.traj[0].md == 0 && !r.tvdBelowKb) w.kb = -in.traj[0].tvd; }
            else st.warnings.push_back(in.name + ": no well head in the project - placed at (0, 0); import or type the head");
            proj.push_back(w); i = int(proj.size()) - 1; ++st.added;
        } else ++st.updated;
        GeoWell& w = proj[size_t(i)];
        if (r.trajMode == 0) {
            w.traj = in.traj; w.trajSrc = 0;
            minimumCurvature(w);
        } else {
            std::vector<Station> s = in.traj;
            for (auto& q : s) {
                if (r.tvdBelowKb) q.tvd -= w.kb;
                if (r.trajMode == 2) { q.x += w.x; q.y += w.y; }
            }
            std::string e;
            // positions kept exactly as surveyed; a later head correction shifts them (mergeHeads, trajSrc 1)
            if (!trajFromXyz(w, s, e)) { st.warnings.push_back(e); continue; }
        }
        if (!w.traj.empty()) w.td = std::max(w.td, w.traj.back().md);
    }
    return st;
}
MergeStats mergeTops(std::vector<GeoWell>& proj, const std::vector<TopRow>& tops, bool createMissing) {
    MergeStats st;
    std::set<std::string> missing;
    for (const auto& t : tops) {
        int i = matchWell(proj, t.well);
        if (i < 0) {
            if (!createMissing) { missing.insert(t.well); continue; }
            GeoWell w; w.name = t.well; proj.push_back(w); i = int(proj.size()) - 1; ++st.added;
            st.warnings.push_back(t.well + ": created without a head (0, 0)");
        }
        GeoWell& w = proj[size_t(i)];
        double md = t.md;
        if (!t.hasMd) {
            const double tvdss = t.tvdBelowKb ? t.tvdss - w.kb : t.tvdss;
            md = w.traj.size() < 2 ? tvdss + w.kb : w.mdAtTvd(tvdss);
            if (!std::isfinite(md)) { st.warnings.push_back(w.name + " / " + t.surface + ": depth " + std::to_string(tvdss) + " not reached by the trajectory"); continue; }
        }
        bool rep = false;
        for (auto& x : w.tops) if (upperS(x.surface) == upperS(t.surface)) { x.md = md; rep = true; }
        if (!rep) { Top tp; tp.surface = t.surface; tp.md = md; w.tops.push_back(tp); }
        std::stable_sort(w.tops.begin(), w.tops.end(), [](const Top& a, const Top& b) { return a.md < b.md; });
        ++st.updated;
    }
    for (const auto& m : missing) st.warnings.push_back("well '" + m + "' not in the project - its tops were skipped");
    return st;
}
std::string writeDeviationTable(const GeoWell& w) {
    std::ostringstream o;
    o << "MD\tINC\tAZI\tTVDSS\tX\tY\n";
    char b[200];
    for (const auto& s : w.traj) { std::snprintf(b, sizeof b, "%.3f\t%.4f\t%.4f\t%.3f\t%.3f\t%.3f\n", s.md, s.inc, s.azi, s.tvd, s.x, s.y); o << b; }
    return o.str();
}

// ================================================================ 2-D SEG-Y
namespace {
int rd16(const unsigned char* p, bool le) { return le ? int(int16_t(uint16_t(p[0] | (p[1] << 8)))) : int(int16_t(uint16_t((p[0] << 8) | p[1]))); }
int32_t rd32(const unsigned char* p, bool le) {
    const uint32_t u = le ? (uint32_t(p[0]) | (uint32_t(p[1]) << 8) | (uint32_t(p[2]) << 16) | (uint32_t(p[3]) << 24))
                          : ((uint32_t(p[0]) << 24) | (uint32_t(p[1]) << 16) | (uint32_t(p[2]) << 8) | uint32_t(p[3]));
    return int32_t(u);
}
float sampleOf(const unsigned char* p, int fmt, bool le) {
    switch (fmt) {
        case 1: {
            const uint32_t u = uint32_t(rd32(p, le));
            if ((u & 0x7fffffff) == 0) return 0.0f;
            const int sign = (u >> 31) ? -1 : 1, expo = int((u >> 24) & 0x7f) - 64;
            return float(sign * (double(u & 0x00ffffff) / 16777216.0) * std::pow(16.0, expo));
        }
        case 2: return float(rd32(p, le));
        case 3: return float(rd16(p, le));
        case 5: { const uint32_t u = uint32_t(rd32(p, le)); float f; std::memcpy(&f, &u, 4); return f; }
        case 8: return float(int8_t(p[0]));
        default: return 0.0f;
    }
}
int bytesOfFmt(int f) { return f == 3 ? 2 : f == 8 ? 1 : 4; }
bool validFmt(int f) { return f == 1 || f == 2 || f == 3 || f == 5 || f == 8; }
char ebcdicC(unsigned char c) {
    static const char* tbl =
        "                                                                "
        "           .<(+|&         !$*); -/        ,%_>?         `:#@'=\""
        " abcdefghi       jklmnopqr       ~stuvwxyz                      "
        "{ABCDEFGHI      }JKLMNOPQR      \\ STUVWXYZ      0123456789      ";
    return tbl[c];
}
struct SegyScan {
    bool le = false; int fmt = 5, ns = 0, bps = 4; double dt = 4, t0 = 0; long long start = 3600, tlen = 0, ntr = 0;
    std::string text;
    struct Tr { int il, xl, cdp; double x, y; };
    std::vector<Tr> tr;
};
bool scanSegy(const std::string& path, const SegyByteMap& m, SegyScan& S, std::string& err, bool headersOnly) {
    std::ifstream f(path, std::ios::binary);
    if (!f) { err = "cannot open " + path; return false; }
    f.seekg(0, std::ios::end); const long long size = f.tellg(); f.seekg(0);
    if (size < 3600 + 240) { err = path + ": too small for SEG-Y"; return false; }
    std::vector<unsigned char> th(3200), bh(400);
    f.read(reinterpret_cast<char*>(th.data()), 3200); f.read(reinterpret_cast<char*>(bh.data()), 400);
    bool ascii = true; for (int k = 0; k < 80; ++k) if (th[size_t(k)] >= 0x80) ascii = false;
    S.text.resize(3200); for (int k = 0; k < 3200; ++k) S.text[size_t(k)] = ascii ? char(th[size_t(k)]) : ebcdicC(th[size_t(k)]);
    int fmt = rd16(&bh[24], false); S.le = false;
    if (!validFmt(fmt) && validFmt(rd16(&bh[24], true))) { S.le = true; fmt = rd16(&bh[24], true); }
    if (!validFmt(fmt)) { err = path + ": sample format " + std::to_string(fmt) + " not supported (1, 2, 3, 5, 8)"; return false; }
    S.fmt = fmt; S.bps = bytesOfFmt(fmt);
    int ext = rd16(&bh[304], S.le); if (ext < 0 || ext > 100) ext = 0;
    S.start = 3600 + 3200LL * ext;
    std::vector<unsigned char> h(240);
    f.seekg(S.start); f.read(reinterpret_cast<char*>(h.data()), 240);
    if (!f) { err = path + ": no traces"; return false; }
    S.ns = rd16(&bh[20], S.le); if (S.ns <= 0) S.ns = rd16(&h[114], S.le);
    S.dt = rd16(&bh[16], S.le) / 1000.0; if (S.dt <= 0) S.dt = rd16(&h[116], S.le) / 1000.0; if (S.dt <= 0) S.dt = 4;
    if (S.ns <= 0) { err = path + ": zero samples per trace"; return false; }
    S.t0 = rd16(&h[108], S.le);
    S.tlen = 240 + (long long)S.bps * S.ns;
    S.ntr = (size - S.start) / S.tlen;
    if (S.ntr < 1) { err = path + ": no traces"; return false; }
    auto at = [](int byte1) { return size_t(std::clamp(byte1, 1, 237) - 1); };
    for (long long t = 0; t < S.ntr; ++t) {
        f.seekg(S.start + t * S.tlen); f.read(reinterpret_cast<char*>(h.data()), 240);
        if (!f) break;
        SegyScan::Tr r;
        r.il = rd32(&h[at(m.inlineByte)], S.le); r.xl = rd32(&h[at(m.xlineByte)], S.le); r.cdp = rd32(&h[at(m.cdpByte)], S.le);
        const int sc = rd16(&h[at(m.scalarByte)], S.le);
        const double scale = sc == 0 ? 1.0 : sc > 0 ? double(sc) : 1.0 / double(-sc);
        r.x = rd32(&h[at(m.xByte)], S.le) * scale; r.y = rd32(&h[at(m.yByte)], S.le) * scale;
        if (m.fallbackXY && r.x == 0 && r.y == 0) { r.x = rd32(&h[72], S.le) * scale; r.y = rd32(&h[76], S.le) * scale; }
        S.tr.push_back(r);
        if (headersOnly && S.tr.size() >= 20000) break;
    }
    return !S.tr.empty();
}
}  // namespace

double SeismicLine::length() const { const auto d = distance(); return d.empty() ? 0 : d.back(); }
std::vector<double> SeismicLine::distance() const {
    std::vector<double> d(x.size(), 0.0);
    for (size_t k = 1; k < x.size(); ++k) d[k] = d[k - 1] + std::hypot(x[k] - x[k - 1], y[k] - y[k - 1]);
    return d;
}
int segyGeometryKind(const std::string& path, const SegyByteMap& m, std::string& info, std::string& err) {
    SegyScan S;
    if (!scanSegy(path, m, S, err, true)) return 0;
    std::set<int> il, xl; std::set<std::pair<int, int>> pr;
    for (auto& t : S.tr) { il.insert(t.il); xl.insert(t.xl); pr.insert({t.il, t.xl}); }
    const double n = double(S.tr.size());
    char b[300];
    std::snprintf(b, sizeof b, "%lld traces x %d samples (dt %.3g ms, format %d%s); %zu inline / %zu crossline values in the scanned headers",
                  S.ntr, S.ns, S.dt, S.fmt, S.le ? " LE" : "", il.size(), xl.size());
    info = b;
    const bool lattice = il.size() > 1 && xl.size() > 1 && double(il.size()) * double(xl.size()) <= 4.0 * n && pr.size() > 0.5 * n;
    return lattice ? 3 : 2;
}
SeismicLine readSegy2D(const std::string& path, const SegyByteMap& m, std::string& err) {
    SeismicLine L;
    SegyScan S;
    if (!scanSegy(path, m, S, err, false)) return L;
    if (double(S.tr.size()) * S.ns > 3e8) { err = path + ": line too large for this machine"; return L; }
    std::ifstream f(path, std::ios::binary);
    L.path = path;
    const size_t sl = path.find_last_of("/\\");
    L.name = path.substr(sl == std::string::npos ? 0 : sl + 1);
    if (L.name.find('.') != std::string::npos) L.name = L.name.substr(0, L.name.rfind('.'));
    L.textHeader = S.text; L.ns = S.ns; L.dt = S.dt; L.t0 = S.t0; L.format = S.fmt;
    L.a.resize(S.tr.size() * size_t(S.ns));
    std::vector<unsigned char> buf(size_t(S.bps) * size_t(S.ns));
    for (size_t t = 0; t < S.tr.size(); ++t) {
        f.seekg(S.start + (long long)t * S.tlen + 240);
        f.read(reinterpret_cast<char*>(buf.data()), std::streamsize(buf.size()));
        for (int s = 0; s < S.ns; ++s) L.a[t * size_t(S.ns) + size_t(s)] = sampleOf(&buf[size_t(s) * size_t(S.bps)], S.fmt, S.le);
        L.x.push_back(S.tr[t].x); L.y.push_back(S.tr[t].y);
        L.cdp.push_back(S.tr[t].cdp != 0 ? S.tr[t].cdp : S.tr[t].xl != 0 ? S.tr[t].xl : int(t + 1));
    }
    bool allZero = true; for (size_t t = 0; t < L.x.size(); ++t) if (L.x[t] != 0 || L.y[t] != 0) allZero = false;
    if (allZero) err = path + ": every trace has X = Y = 0 - set the coordinate byte locations (e.g. CDP X 181 / Y 185 or source X 73 / Y 77)";
    return L;
}
bool writeSegy2D(const std::string& path, const SeismicLine& l) {
    std::ofstream f(path, std::ios::binary);
    if (!f) return false;
    std::string th(3200, ' ');
    const std::string c1 = "C 1 RESSIM 2-D LINE " + l.name;
    std::memcpy(&th[0], c1.data(), std::min<size_t>(c1.size(), 80));
    f.write(th.data(), 3200);
    unsigned char bh[400] = {0};
    auto w16 = [](unsigned char* p, int v) { p[0] = (unsigned char)((v >> 8) & 255); p[1] = (unsigned char)(v & 255); };
    auto w32 = [](unsigned char* p, int32_t v) { const uint32_t u = uint32_t(v); p[0] = (unsigned char)(u >> 24); p[1] = (unsigned char)(u >> 16); p[2] = (unsigned char)(u >> 8); p[3] = (unsigned char)u; };
    w16(&bh[16], int(std::lround(l.dt * 1000))); w16(&bh[20], l.ns); w16(&bh[24], 5); w16(&bh[300], 0x0100);
    f.write(reinterpret_cast<const char*>(bh), 400);
    for (int t = 0; t < l.ntr(); ++t) {
        unsigned char h[240] = {0};
        w32(&h[0], t + 1); w32(&h[20], t < int(l.cdp.size()) ? l.cdp[size_t(t)] : t + 1);
        w16(&h[70], -100);
        w32(&h[180], int32_t(std::lround(l.x[size_t(t)] * 100))); w32(&h[184], int32_t(std::lround(l.y[size_t(t)] * 100)));
        w32(&h[188], 1); w32(&h[192], t < int(l.cdp.size()) ? l.cdp[size_t(t)] : t + 1);
        w16(&h[108], int(std::lround(l.t0))); w16(&h[114], l.ns); w16(&h[116], int(std::lround(l.dt * 1000)));
        f.write(reinterpret_cast<const char*>(h), 240);
        std::vector<unsigned char> s(size_t(4 * l.ns));
        for (int k = 0; k < l.ns; ++k) { uint32_t u; const float v = l.at(t, k); std::memcpy(&u, &v, 4); w32(&s[size_t(4 * k)], int32_t(u)); }
        f.write(reinterpret_cast<const char*>(s.data()), std::streamsize(s.size()));
    }
    return bool(f);
}

std::vector<double> trackLine(const SeismicLine& l, const std::vector<std::pair<int, double>>& seeds, int mode, double window, double maxJump) {
    const int n = l.ntr();
    std::vector<double> t(size_t(std::max(0, n)), kNaN);
    if (n < 1 || l.ns < 3) return t;
    // event time near tPrev on trace tr (sub-sample, parabolic); NaN if none; 'amp' = event amplitude
    auto pick = [&](int tr, double tPrev, double& amp) {
        const int s0 = std::max(1, int(std::floor((tPrev - window - l.t0) / l.dt))), s1 = std::min(l.ns - 2, int(std::ceil((tPrev + window - l.t0) / l.dt)));
        int best = -1; double bv = 0, bd = 1e300;
        for (int s = s0; s <= s1; ++s) {
            const double a = l.at(tr, s), am = l.at(tr, s - 1), ap = l.at(tr, s + 1);
            if (mode == 2) {
                if (a > 0 && ap <= 0) { const double d = std::fabs(l.t0 + s * l.dt - tPrev); if (d < bd) { bd = d; best = s; } }
                continue;
            }
            const bool ext = mode == 0 ? (a >= am && a >= ap && a > 0) : (a <= am && a <= ap && a < 0);
            if (!ext) continue;
            const double d = std::fabs(l.t0 + s * l.dt - tPrev);
            if (best < 0 || d < bd - 1e-9 || (std::fabs(d - bd) < 1e-9 && std::fabs(a) > std::fabs(bv))) { best = s; bd = d; bv = a; }
        }
        if (best < 0) { amp = 0; return kNaN; }
        const double a = l.at(tr, best), am = l.at(tr, best - 1), ap = l.at(tr, best + 1);
        double off = 0;
        if (mode == 2) { off = std::fabs(a - ap) > 1e-30 ? a / (a - ap) : 0; amp = std::fabs(am) + std::fabs(ap); }
        else { const double den = am - 2 * a + ap; off = std::fabs(den) > 1e-30 ? std::clamp(0.5 * (am - ap) / den, -0.5, 0.5) : 0; amp = std::fabs(a); }
        return l.t0 + (best + off) * l.dt;
    };
    for (const auto& sd : seeds) {
        const int tr0 = std::clamp(sd.first, 0, n - 1);
        if (std::isfinite(t[size_t(tr0)])) continue;
        double a0; const double ts = pick(tr0, sd.second, a0);
        if (!std::isfinite(ts)) continue;
        t[size_t(tr0)] = ts;
        for (int dir = -1; dir <= 1; dir += 2) {
            double prev = ts;
            for (int tr = tr0 + dir; tr >= 0 && tr < n; tr += dir) {
                if (std::isfinite(t[size_t(tr)])) break;
                double a; const double tt = pick(tr, prev, a);
                if (!std::isfinite(tt) || std::fabs(tt - prev) > maxJump || a < 0.05 * a0) break;
                t[size_t(tr)] = tt; prev = tt;
            }
        }
    }
    return t;
}
std::vector<double> interpolatePicks(int ntr, std::vector<std::pair<int, double>> pts) {
    std::vector<double> t(size_t(std::max(0, ntr)), kNaN);
    std::sort(pts.begin(), pts.end());
    std::vector<std::pair<int, double>> u;
    for (auto& p : pts) { if (p.first < 0 || p.first >= ntr || !std::isfinite(p.second)) continue; if (!u.empty() && u.back().first == p.first) u.back() = p; else u.push_back(p); }
    if (u.size() == 1) t[size_t(u[0].first)] = u[0].second;
    for (size_t k = 1; k < u.size(); ++k)
        for (int tr = u[k - 1].first; tr <= u[k].first; ++tr) {
            const double f = u[k].first > u[k - 1].first ? double(tr - u[k - 1].first) / double(u[k].first - u[k - 1].first) : 0.0;
            t[size_t(tr)] = u[k - 1].second + f * (u[k].second - u[k - 1].second);
        }
    return t;
}
std::vector<LineTie> lineIntersections(const std::vector<SeismicLine>& lines) {
    std::vector<LineTie> out;
    for (size_t a = 0; a < lines.size(); ++a)
        for (size_t b = a + 1; b < lines.size(); ++b) {
            const SeismicLine &A = lines[a], &B = lines[b];
            const int na = A.ntr(), nb = B.ntr();
            if (na < 2 || nb < 2) continue;
            const int blk = 32;
            auto box = [](const SeismicLine& L, int i0, int i1, double bb[4]) {
                bb[0] = bb[2] = 1e300; bb[1] = bb[3] = -1e300;
                for (int i = i0; i <= i1; ++i) { bb[0] = std::min(bb[0], L.x[size_t(i)]); bb[1] = std::max(bb[1], L.x[size_t(i)]); bb[2] = std::min(bb[2], L.y[size_t(i)]); bb[3] = std::max(bb[3], L.y[size_t(i)]); }
            };
            for (int i0 = 0; i0 < na - 1; i0 += blk) {
                const int i1 = std::min(na - 1, i0 + blk);
                double ba[4]; box(A, i0, i1, ba);
                for (int j0 = 0; j0 < nb - 1; j0 += blk) {
                    const int j1 = std::min(nb - 1, j0 + blk);
                    double bbx[4]; box(B, j0, j1, bbx);
                    if (ba[1] < bbx[0] || bbx[1] < ba[0] || ba[3] < bbx[2] || bbx[3] < ba[2]) continue;
                    for (int i = i0; i < i1; ++i)
                        for (int j = j0; j < j1; ++j) {
                            const double px = A.x[size_t(i)], py = A.y[size_t(i)], rx = A.x[size_t(i + 1)] - px, ry = A.y[size_t(i + 1)] - py;
                            const double qx = B.x[size_t(j)], qy = B.y[size_t(j)], sx = B.x[size_t(j + 1)] - qx, sy = B.y[size_t(j + 1)] - qy;
                            const double den = rx * sy - ry * sx;
                            if (std::fabs(den) < 1e-12) continue;
                            const double u = ((qx - px) * sy - (qy - py) * sx) / den, v = ((qx - px) * ry - (qy - py) * rx) / den;
                            if (u < 0 || u >= 1 || v < 0 || v >= 1) continue;
                            LineTie T; T.a = int(a); T.b = int(b); T.ta = u < 0.5 ? i : i + 1; T.tb = v < 0.5 ? j : j + 1;
                            T.x = px + u * rx; T.y = py + u * ry;
                            out.push_back(T);
                        }
                }
            }
        }
    return out;
}
std::vector<P3> lineHorizonPoints(const std::vector<SeismicLine>& lines, const LineHorizon& h, const VelocityModel& v, int step) {
    std::vector<P3> p;
    step = std::max(1, step);
    for (size_t l = 0; l < lines.size() && l < h.t.size(); ++l)
        for (int tr = 0; tr < lines[l].ntr() && tr < int(h.t[l].size()); tr += step) {
            const double t = h.t[l][size_t(tr)];
            if (std::isfinite(t)) p.push_back({lines[l].x[size_t(tr)], lines[l].y[size_t(tr)], timeToDepth(v, t)});
        }
    return p;
}
Surface latticeFor(const std::vector<P3>& pts, int n) {
    Surface s;
    if (pts.empty() || n < 2) return s;
    double x0 = 1e300, x1 = -1e300, y0 = 1e300, y1 = -1e300;
    for (auto& p : pts) { x0 = std::min(x0, p[0]); x1 = std::max(x1, p[0]); y0 = std::min(y0, p[1]); y1 = std::max(y1, p[1]); }
    double wx = std::max(x1 - x0, 1.0), wy = std::max(y1 - y0, 1.0);
    x0 -= 0.05 * wx; y0 -= 0.05 * wy; wx *= 1.1; wy *= 1.1;
    s.nx = n; s.ny = std::max(2, int(std::lround(n * wy / wx)));
    s.x0 = x0; s.y0 = y0; s.dx = wx / (s.nx - 1); s.dy = wy / (s.ny - 1);
    s.z.assign(size_t(s.nx) * size_t(s.ny), kNaN);
    return s;
}

}  // namespace geo
}  // namespace ressim
