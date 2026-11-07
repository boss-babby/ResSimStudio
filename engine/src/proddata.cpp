// C38 production-data import, QC and history SCHEDULE generation (see proddata.h)
#include "ressim/proddata.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <set>
#include <sstream>

namespace ressim {
namespace pdata {

namespace {
const double NaN = std::numeric_limits<double>::quiet_NaN();
const char* kMon[12] = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};

std::string upper(std::string s) { for (auto& c : s) c = char(std::toupper((unsigned char)c)); return s; }
std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace((unsigned char)s[a])) ++a;
    while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
}
std::string unquote(std::string s) {
    s = trim(s);
    if (s.size() >= 2 && (s.front() == '"' || s.front() == '\'') && s.back() == s.front()) s = s.substr(1, s.size() - 2);
    return trim(s);
}
std::vector<std::string> splitDelim(const std::string& line, char d) {
    std::vector<std::string> out; std::string cur; char q = 0;
    for (char c : line) {
        if (q) { if (c == q) q = 0; cur += c; continue; }
        if (c == '"' || c == '\'') { q = c; cur += c; continue; }
        if (c == d) { out.push_back(cur); cur.clear(); continue; }
        cur += c;
    }
    out.push_back(cur);
    return out;
}
std::vector<std::string> splitWs(const std::string& line) {
    std::vector<std::string> out; std::string cur; char q = 0;
    for (char c : line) {
        if (q) { if (c == q) q = 0; cur += c; continue; }
        if (c == '"' || c == '\'') { q = c; cur += c; continue; }
        if (std::isspace((unsigned char)c)) { if (!cur.empty()) out.push_back(cur); cur.clear(); continue; }
        cur += c;
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}
bool parseNum(const std::string& s0, double& v) {
    std::string s = unquote(s0);
    if (s.empty()) return false;
    const std::string u = upper(s);
    if (u == "NA" || u == "N/A" || u == "NAN" || u == "-" || u == "NULL" || u == "1*") return false;
    std::string t; for (char c : s) if (c != ' ' && c != '_') t += c;   // 1 234 / 1_234
    char* e = nullptr; v = std::strtod(t.c_str(), &e);
    return e && *e == 0 && std::isfinite(v);
}
int monthOf(const std::string& s) {
    const std::string u = upper(unquote(s)).substr(0, 3);
    for (int i = 0; i < 12; ++i) if (u == kMon[i]) return i + 1;
    if (u == "JLY") return 7;
    return 0;
}
int monthLen(int y, int m) {
    static const int L[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    const bool leap = (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
    return m == 2 && leap ? 29 : L[m - 1];
}
long addMonth(long day) {
    int y, m, d; civilFromDays(day, y, m, d);
    const bool lastDay = d == monthLen(y, m);
    if (++m > 12) { m = 1; ++y; }
    d = lastDay ? monthLen(y, m) : std::min(d, monthLen(y, m));
    return daysFromCivil(y, m, d);
}
long subMonth(long day) {
    int y, m, d; civilFromDays(day, y, m, d);
    const bool lastDay = d == monthLen(y, m);
    if (--m < 1) { m = 12; --y; }
    d = lastDay ? monthLen(y, m) : std::min(d, monthLen(y, m));
    return daysFromCivil(y, m, d);
}
std::string fmt(double v) { char b[32]; std::snprintf(b, sizeof(b), "%.7g", v); return b; }
std::string isoDate(long day) { int y, m, d; civilFromDays(day, y, m, d); char b[16]; std::snprintf(b, sizeof(b), "%04d-%02d-%02d", y, m, d); return b; }
std::string quoteWell(const std::string& w) { return "'" + w + "'"; }
std::string stripComment(const std::string& line) {
    char q = 0;
    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (q) { if (c == q) q = 0; continue; }
        if (c == '\'' || c == '"') { q = c; continue; }
        if (c == '-' && i + 1 < line.size() && line[i + 1] == '-') return line.substr(0, i);
    }
    return line;
}
bool isKeyword(const std::string& stripped, std::string& name) {
    const auto t = splitWs(stripped);
    if (t.size() != 1) return false;
    const std::string& s = t[0];
    if (s.size() > 8 || !std::isupper((unsigned char)s[0])) return false;
    for (char c : s) if (!(std::isupper((unsigned char)c) || std::isdigit((unsigned char)c) || c == '_' || c == '-' || c == '+')) return false;
    name = s;
    return true;
}
}  // namespace

Record::Record() { for (double& x : v) x = NaN; }

int Table::nErrors() const { int n = 0; for (const auto& i : issues) n += i.severity == 2; return n; }
int Table::nWarnings() const { int n = 0; for (const auto& i : issues) n += i.severity == 1; return n; }

long daysFromCivil(int y, int m, int d) {
    y -= m <= 2;
    const long era = (y >= 0 ? y : y - 399) / 400;
    const unsigned yoe = unsigned(y - era * 400);
    const unsigned doy = unsigned((153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1);
    const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return era * 146097 + long(doe) - 719468;
}
void civilFromDays(long z, int& y, int& m, int& d) {
    z += 719468;
    const long era = (z >= 0 ? z : z - 146096) / 146097;
    const unsigned doe = unsigned(z - era * 146097);
    const unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    const long yy = long(yoe) + era * 400;
    const unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    const unsigned mp = (5 * doy + 2) / 153;
    d = int(doy - (153 * mp + 2) / 5 + 1);
    m = int(mp < 10 ? mp + 3 : mp - 9);
    y = int(yy + (m <= 2));
}
std::string dateString(long day) { int y, m, d; civilFromDays(day, y, m, d); return std::to_string(d) + " '" + kMon[m - 1] + "' " + std::to_string(y); }

bool parseDate(const std::string& s0, bool monthFirst, long& day) {
    std::string s = unquote(s0);
    if (s.empty()) return false;
    // drop a time part (2020-01-31T00:00:00 / 2020-01-31 00:00)
    if (const size_t p = s.find('T'); p != std::string::npos && p >= 8 && std::isdigit((unsigned char)s[p - 1])) s = s.substr(0, p);
    std::vector<std::string> f; std::string cur;
    for (char c : s) {
        if (c == '-' || c == '/' || c == '.' || c == ' ') { if (!cur.empty()) f.push_back(cur); cur.clear(); }
        else cur += c;
    }
    if (!cur.empty()) f.push_back(cur);
    if (f.size() >= 4 && f[3].find(':') != std::string::npos) f.resize(3);
    auto isInt = [](const std::string& x) { if (x.empty()) return false; for (char c : x) if (!std::isdigit((unsigned char)c)) return false; return true; };
    int y = 0, m = 0, d = 0;
    if (f.size() == 2) {                                     // 2020-01 / 01/2020 / JAN-2020
        if (isInt(f[0]) && f[0].size() == 4 && isInt(f[1])) { y = std::atoi(f[0].c_str()); m = std::atoi(f[1].c_str()); d = 1; }
        else if (isInt(f[1]) && f[1].size() == 4) { y = std::atoi(f[1].c_str()); m = isInt(f[0]) ? std::atoi(f[0].c_str()) : monthOf(f[0]); d = 1; }
        else return false;
    } else if (f.size() == 3) {
        if (isInt(f[0]) && f[0].size() == 4) {               // ISO y m d
            y = std::atoi(f[0].c_str()); m = isInt(f[1]) ? std::atoi(f[1].c_str()) : monthOf(f[1]); d = std::atoi(f[2].c_str());
        } else if (!isInt(f[1])) {                           // 01-JAN-2020 / 1 JAN 2020
            d = std::atoi(f[0].c_str()); m = monthOf(f[1]); y = std::atoi(f[2].c_str());
        } else if (!isInt(f[0])) {                           // JAN 1 2020
            m = monthOf(f[0]); d = std::atoi(f[1].c_str()); y = std::atoi(f[2].c_str());
        } else {
            int a = std::atoi(f[0].c_str()), b = std::atoi(f[1].c_str());
            y = std::atoi(f[2].c_str());
            if (a > 12) { d = a; m = b; } else if (b > 12) { m = a; d = b; }
            else if (monthFirst) { m = a; d = b; } else { d = a; m = b; }
        }
        if (y < 100 && y >= 0 && f[f.size() - 1].size() <= 2) y += y < 50 ? 2000 : 1900;   // 2-digit years
    } else return false;
    if (y < 1800 || y > 2300 || m < 1 || m > 12 || d < 1 || d > monthLen(y, m)) return false;
    day = daysFromCivil(y, m, d);
    return true;
}

const char* colName(int role) {
    static const char* n[C_NCOL] = {"WELL", "DATE", "YEAR", "MONTH", "DAY", "OIL", "WATER", "GAS", "WATER INJ", "GAS INJ", "BHP",
                                    "THP", "DAYS ON", "HOURS ON", "STATIC P"};
    return role >= 0 && role < C_NCOL ? n[role] : "(ignored)";
}

int columnRole(const std::string& h0) {
    std::string h; for (char c : upper(unquote(h0))) if (std::isalnum((unsigned char)c)) h += c;
    // strip a unit suffix in brackets was removed with the punctuation: "OIL (SM3)" -> "OILSM3" - match by prefix lists
    struct A { int role; const char* names[16]; };
    static const A al[] = {
        {C_WELL, {"WELL", "WELLNAME", "WELLID", "UWI", "NAME", "WELLBORE", nullptr}},
        {C_DATE, {"DATE", "PRODDATE", "PRODUCTIONDATE", "MONTHDATE", "PERIOD", "TIME", "DATETIME", nullptr}},
        {C_YEAR, {"YEAR", "YR", nullptr}},
        {C_MONTH, {"MONTH", "MON", nullptr}},
        {C_DAY, {"DAY", nullptr}},
        {C_OIL, {"OIL", "ORAT", "WOPR", "WOPRH", "QO", "OILRATE", "OILVOL", "OILVOLUME", "OILPROD", "OILPRODUCTION", "WOPT", "WOPTH", "BOPD", "STBD", nullptr}},
        {C_WATER, {"WATER", "WRAT", "WWPR", "WWPRH", "QW", "WATERRATE", "WATERVOL", "WATERVOLUME", "WATERPROD", "WATERPRODUCTION", "BWPD", nullptr}},
        {C_GAS, {"GAS", "GRAT", "WGPR", "WGPRH", "QG", "GASRATE", "GASVOL", "GASVOLUME", "GASPROD", "GASPRODUCTION", "MSCFD", nullptr}},
        {C_WINJ, {"WINJ", "WATERINJ", "WATERINJECTION", "WWIR", "WWIRH", "WI", "INJWATER", "WATINJ", nullptr}},
        {C_GINJ, {"GINJ", "GASINJ", "GASINJECTION", "WGIR", "WGIRH", "GI", "INJGAS", nullptr}},
        {C_BHP, {"BHP", "WBHP", "WBHPH", "PWF", "FBHP", nullptr}},
        {C_THP, {"THP", "WTHP", "WTHPH", "FTHP", "WHP", nullptr}},
        {C_DAYS, {"DAYS", "DAYSON", "UPTIME", "UPTIMEDAYS", "ONSTREAM", "ONSTREAMDAYS", "PRODDAYS", "PRODUCINGDAYS", "INJDAYS", nullptr}},
        {C_HOURS, {"HOURS", "HOURSON", "UPTIMEHOURS", "ONSTREAMHOURS", "PRODHOURS", nullptr}},
        {C_PRES, {"PRES", "WBP", "SIBHP", "STATICPRESSURE", "RESPRESSURE", "PR", "WBP9", "WBP5", nullptr}},
    };
    for (const auto& a : al) for (int i = 0; a.names[i]; ++i) if (h == a.names[i]) return a.role;
    // with a unit suffix: OILSM3, GASKSM3, WATERM3D, BHPBAR ...
    for (const auto& a : al) {
        if (a.role == C_DAY || a.role == C_MONTH || a.role == C_YEAR || a.role == C_DATE) continue;
        for (int i = 0; a.names[i]; ++i) {
            const std::string n = a.names[i];
            if (n.size() >= 3 && h.size() > n.size() && h.compare(0, n.size(), n) == 0) {
                const std::string rest = h.substr(n.size());
                static const char* units[] = {"SM3", "KSM3", "MSM3", "M3", "M3D", "SM3D", "STB", "BBL", "MSTB", "MSCF", "MMSCF", "SCF",
                                              "BAR", "BARSA", "BARA", "PSI", "PSIA", "KPA", "MPA", "D", "PERDAY", "DAY", "RATE", "VOL", "VOLUME", nullptr};
                for (int u = 0; units[u]; ++u) if (rest == units[u]) return a.role;
            }
        }
    }
    return C_NONE;
}

bool parseTable(const std::string& text, const Options& o, Table& t) {
    t = Table{};
    auto issue = [&](int sev, int row, const std::string& well, long day, const std::string& msg) {
        Issue i; i.severity = sev; i.row = row; i.well = well; i.date = day ? isoDate(day) : std::string(); i.msg = msg; t.issues.push_back(i);
    };
    std::vector<std::pair<int, std::string>> lines;
    {
        std::istringstream in(text); std::string l; int n = 0;
        while (std::getline(in, l)) {
            ++n;
            if (!l.empty() && l.back() == '\r') l.pop_back();
            if (n == 1 && l.size() >= 3 && (unsigned char)l[0] == 0xEF && (unsigned char)l[1] == 0xBB && (unsigned char)l[2] == 0xBF) l = l.substr(3);
            const std::string tl = trim(l);
            if (tl.empty() || tl[0] == '#' || tl.compare(0, 2, "--") == 0) continue;
            lines.push_back({n, l});
        }
    }
    if (lines.empty()) { issue(2, 0, "", 0, "the table is empty"); return false; }
    const std::string& head = lines[0].second;
    char delim = 0;
    const size_t nc = std::count(head.begin(), head.end(), ','), nt = std::count(head.begin(), head.end(), '\t'), ns = std::count(head.begin(), head.end(), ';');
    if (nt > 0 && nt >= nc) delim = '\t'; else if (nc > 0 && nc >= ns) delim = ','; else if (ns > 0) delim = ';';
    auto fields = [&](const std::string& l) { return delim ? splitDelim(l, delim) : splitWs(l); };
    t.columns = fields(head);
    for (auto& c : t.columns) c = unquote(c);
    int idx[C_NCOL]; std::fill(idx, idx + C_NCOL, -1);
    for (size_t c = 0; c < t.columns.size(); ++c) {
        int r = columnRole(t.columns[c]);
        if (c < o.roleOverride.size() && o.roleOverride[c] != -2) r = o.roleOverride[c] >= C_NONE && o.roleOverride[c] < C_NCOL ? o.roleOverride[c] : C_NONE;
        if (r != C_NONE && idx[r] >= 0) { issue(1, lines[0].first, "", 0, "column '" + t.columns[c] + "' repeats " + colName(r) + " - ignored"); r = C_NONE; }
        if (r != C_NONE) idx[r] = int(c);
        else if (!t.columns[c].empty()) issue(0, lines[0].first, "", 0, "column '" + t.columns[c] + "' not recognised - ignored");
        t.roles.push_back(r);
    }
    if (idx[C_WELL] < 0) { issue(2, lines[0].first, "", 0, "no WELL column in the header (WELL, WELLNAME, UWI, NAME)"); return false; }
    if (idx[C_DATE] < 0 && (idx[C_YEAR] < 0 || idx[C_MONTH] < 0)) { issue(2, lines[0].first, "", 0, "no DATE column (or YEAR + MONTH columns) in the header"); return false; }
    const bool hasRate = idx[C_OIL] >= 0 || idx[C_WATER] >= 0 || idx[C_GAS] >= 0 || idx[C_WINJ] >= 0 || idx[C_GINJ] >= 0;
    if (!hasRate) { issue(2, lines[0].first, "", 0, "no rate / volume column (OIL, WATER, GAS, WINJ, GINJ)"); return false; }
    const bool hasUp = idx[C_DAYS] >= 0 || idx[C_HOURS] >= 0;
    // ---- rows
    std::map<std::pair<std::string, long>, size_t> seen;
    std::map<std::pair<std::string, long>, int> dupCount;
    std::set<std::string> longNames;
    for (size_t li = 1; li < lines.size(); ++li) {
        const int row = lines[li].first;
        auto f = fields(lines[li].second);
        auto cell = [&](int role) -> std::string { const int c = idx[role]; return c >= 0 && c < int(f.size()) ? f[size_t(c)] : std::string(); };
        Record r; r.row = row; r.well = unquote(cell(C_WELL));
        long day = 0; bool okDate;
        if (idx[C_DATE] >= 0) okDate = parseDate(cell(C_DATE), o.monthFirst, day);
        else {
            double yv, mv, dv = 1;
            okDate = parseNum(cell(C_YEAR), yv) && (parseNum(cell(C_MONTH), mv) || (mv = monthOf(cell(C_MONTH))) > 0);
            if (okDate && idx[C_DAY] >= 0 && !parseNum(cell(C_DAY), dv)) dv = 1;
            okDate = okDate && yv > 1800 && mv >= 1 && mv <= 12 && dv >= 1 && dv <= monthLen(int(yv), int(mv));
            if (okDate) day = daysFromCivil(int(yv), int(mv), int(dv));
        }
        if (r.well.empty() || !okDate) {
            if (li == 1 && !r.well.empty()) { issue(0, row, "", 0, "second header row (units?) skipped"); continue; }
            issue(2, row, r.well, 0, r.well.empty() ? "no well name - row dropped" : "date '" + trim(cell(C_DATE)) + "' not understood - row dropped");
            continue;
        }
        r.day = day;
        if (r.well.size() > 8 && longNames.insert(r.well).second) issue(1, row, r.well, 0, "well name longer than 8 characters (ECLIPSE limit)");
        for (int role = C_OIL; role < C_NCOL; ++role) {
            if (idx[role] < 0) continue;
            const std::string s = cell(role);
            double v;
            if (!parseNum(s, v)) { if (!trim(s).empty() && upper(unquote(s)) != "NA" && unquote(s) != "-" && upper(unquote(s)) != "1*") issue(1, row, r.well, day, std::string(colName(role)) + " value '" + trim(s) + "' not a number - missing"); continue; }
            if (role == C_BHP || role == C_THP || role == C_PRES) { if (v <= 0) continue; }
            else if (v < 0) { issue(1, row, r.well, day, std::string("negative ") + colName(role) + " " + fmt(v) + " set to 0"); v = 0; }
            if (role == C_OIL || role == C_WATER || role == C_WINJ) v *= o.liqFactor;
            if (role == C_GAS || role == C_GINJ) v *= o.gasFactor;
            r.v[role] = v;
        }
        if (std::isnan(r.v[C_DAYS]) && !std::isnan(r.v[C_HOURS])) r.v[C_DAYS] = r.v[C_HOURS] / 24.0;
        const auto key = std::make_pair(r.well, day);
        auto it = seen.find(key);
        if (it == seen.end()) { seen[key] = t.recs.size(); t.recs.push_back(r); continue; }
        Record& a = t.recs[it->second];
        const int k = ++dupCount[key];
        for (int role = C_OIL; role < C_NCOL; ++role) {
            if (std::isnan(r.v[role])) continue;
            if (std::isnan(a.v[role])) { a.v[role] = r.v[role]; continue; }
            const bool add = (o.volumes && role >= C_OIL && role <= C_GINJ) || role == C_DAYS || role == C_HOURS;
            a.v[role] = add ? a.v[role] + r.v[role] : (a.v[role] * k + r.v[role]) / (k + 1);
        }
        issue(1, row, r.well, day, std::string("duplicate well / date - ") + (o.volumes ? "volumes summed" : "rates averaged") + " with line " + std::to_string(a.row));
    }
    if (t.recs.empty()) { issue(2, 0, "", 0, "no valid data rows"); return false; }
    std::stable_sort(t.recs.begin(), t.recs.end(), [](const Record& a, const Record& b) { return a.well != b.well ? a.well < b.well : a.day < b.day; });
    {
        std::set<std::string> ws; std::vector<std::pair<int, std::string>> first;
        for (const auto& r : t.recs) if (ws.insert(r.well).second) first.push_back({r.row, r.well});
        int minRow = 1 << 30;
        std::map<std::string, int> fr; for (const auto& r : t.recs) { auto& x = fr[r.well]; x = x ? std::min(x, r.row) : r.row; minRow = std::min(minRow, r.row); }
        for (auto& p : first) p.first = fr[p.second];
        std::sort(first.begin(), first.end());
        for (auto& p : first) t.wells.push_back(p.second);
    }
    for (const auto& r : t.recs) t.dates.push_back(r.day);
    std::sort(t.dates.begin(), t.dates.end()); t.dates.erase(std::unique(t.dates.begin(), t.dates.end()), t.dates.end());
    // ---- periods
    const int nd = int(t.dates.size());
    std::vector<long> diffs; for (int i = 1; i < nd; ++i) diffs.push_back(t.dates[size_t(i)] - t.dates[size_t(i) - 1]);
    long med = 30;
    if (!diffs.empty()) { auto d2 = diffs; std::nth_element(d2.begin(), d2.begin() + d2.size() / 2, d2.end()); med = std::max(1L, d2[d2.size() / 2]); }
    const bool monthly = med >= 28 && med <= 31;
    t.periods.resize(size_t(nd));
    for (int i = 0; i < nd; ++i) {
        Period& p = t.periods[size_t(i)];
        if (!o.datesAreEnds) { p.d0 = t.dates[size_t(i)]; p.d1 = i + 1 < nd ? t.dates[size_t(i) + 1] : (monthly ? addMonth(p.d0) : p.d0 + med); }
        else {   // C38: an end date is inclusive (31 JAN closes January: the period ends at 1 FEB 00:00, the schedule time)
            p.d1 = t.dates[size_t(i)] + 1;
            p.d0 = i > 0 ? t.dates[size_t(i) - 1] + 1 : (monthly ? subMonth(t.dates[size_t(i)]) + 1 : p.d1 - med);
        }
    }
    if (nd > 2) {
        for (int i = 1; i < nd; ++i) if (diffs[size_t(i) - 1] > 3 * med) issue(0, 0, "", t.dates[size_t(i) - 1], "gap of " + std::to_string(diffs[size_t(i) - 1]) + " days to the next date (typical " + std::to_string(med) + ") - one long period");
    }
    auto periodOf = [&](long day) { return int(std::lower_bound(t.dates.begin(), t.dates.end(), day) - t.dates.begin()); };
    // ---- per well and period
    for (const auto& w : t.wells) {
        std::vector<const Record*> byP(size_t(nd), nullptr);
        for (const auto& r : t.recs) if (r.well == w) byP[size_t(periodOf(r.day))] = &r;
        int pFirst = -1; for (int p = 0; p < nd; ++p) if (byP[size_t(p)]) { pFirst = p; break; }
        bool lastInj = false; int nMissing = 0;
        double prevWor = -1, prevGor = -1;
        for (int p = pFirst; p < nd; ++p) {
            const Record* r = byP[size_t(p)];
            WellPeriod x; x.well = w; x.period = p;
            const double L = double(t.periods[size_t(p)].d1 - t.periods[size_t(p)].d0);
            if (!r) {
                if (!o.shutMissing) continue;
                ++nMissing; x.shut = true; x.inj = lastInj; t.wp.push_back(x); continue;
            }
            auto val = [&](int role) { return std::isnan(r->v[role]) ? 0.0 : r->v[role]; };
            double up = hasUp && !std::isnan(r->v[C_DAYS]) ? r->v[C_DAYS] : L;
            if (up > L + 1e-6) { issue(1, r->row, w, r->day, "uptime " + fmt(up) + " d longer than the period (" + fmt(L) + " d) - clamped"); up = L; }
            const double qo = val(C_OIL), qw = val(C_WATER), qg = val(C_GAS), wi = val(C_WINJ), gi = val(C_GINJ);
            const bool prod = qo > 0 || qw > 0 || qg > 0, inj = wi > 0 || gi > 0;
            if (prod && inj) issue(1, r->row, w, r->day, "produces and injects in the same period - injection used");
            if (wi > 0 && gi > 0) issue(1, r->row, w, r->day, "water and gas injection in the same period - water used (alternate WAG periods for both)");
            bool useEf = o.useWefac && hasUp && !std::isnan(r->v[C_DAYS]);
            if (useEf && up <= 1e-9 && (prod || inj)) { issue(1, r->row, w, r->day, "volumes with 0 days on - calendar rate used"); useEf = false; }
            const double div = o.volumes ? (useEf ? up : L) : 1.0;
            x.wefac = useEf ? std::max(0.0, std::min(1.0, up / L)) : 1.0;
            if (!o.volumes && useEf && x.wefac > 0) {}   // rates given = flowing rates
            x.orat = qo / div; x.wrat = qw / div; x.grat = qg / div; x.winj = wi / div; x.ginj = gi / div;
            x.bhp = val(C_BHP); x.thp = val(C_THP); x.pres = val(C_PRES);
            x.inj = inj; x.shut = !(prod || inj) || (useEf && up <= 1e-9);
            if (x.inj && wi > 0) x.ginj = 0;
            lastInj = x.inj;
            if (!x.inj && qo > 0) {
                const double wor = qw / qo, gor = qg / qo;
                auto jump = [&](double a, double b) { return a > 0 && b > 0 && (b > o.outlierFactor * a || a > o.outlierFactor * b); };
                if (prevWor > 0.01 && wor > 0.01 && jump(prevWor, wor)) issue(1, r->row, w, r->day, "water-oil ratio jumps " + fmt(prevWor) + " -> " + fmt(wor) + " (outlier?)");
                if (prevGor > 0 && gor > 0 && jump(prevGor, gor)) issue(1, r->row, w, r->day, "GOR jumps " + fmt(prevGor) + " -> " + fmt(gor) + " (outlier?)");
                prevWor = wor; prevGor = gor;
            }
            t.wp.push_back(x);
        }
        if (nMissing > 0) issue(0, 0, w, 0, std::to_string(nMissing) + " period(s) without data after the first record - shut there");
        bool any = false; for (const auto& x : t.wp) if (x.well == w && !x.shut) any = true;
        if (!any) issue(1, 0, w, 0, "never produces or injects in the table");
    }
    issue(0, 0, "", 0, std::to_string(t.recs.size()) + " rows, " + std::to_string(t.wells.size()) + " wells, " + std::to_string(nd) + " periods " +
          isoDate(t.periods.front().d0) + " .. " + isoDate(t.periods.back().d1) + (monthly ? " (monthly)" : "") +
          (o.volumes ? ", volumes -> rates" : ", rates") + (hasUp && o.useWefac ? " with uptime (WEFAC)" : ""));
    return true;
}

std::string periodKeywords(const Table& t, const Options& o, int period) {
    std::ostringstream h, ij, ef;
    // WEFAC only when it changes against the previous period of the well (default 1)
    std::map<std::string, std::pair<int, double>> prevEf;   // well -> (period, wefac) of its last open period before this one
    for (const auto& y : t.wp)
        if (y.period < period && !y.shut) { auto& e = prevEf[y.well]; if (e.first <= y.period) e = {y.period, y.wefac}; }
    for (const auto& x : t.wp) {
        if (x.period != period) continue;
        if (!x.inj) {
            if (x.shut) h << "  " << quoteWell(x.well) << " SHUT " << o.prodControl << " 0 0 0 /\n";
            else {
                h << "  " << quoteWell(x.well) << " OPEN " << o.prodControl << " " << fmt(x.orat) << " " << fmt(x.wrat) << " " << fmt(x.grat);
                const bool thp = o.writeThp && x.thp > 0, bhp = o.writeBhp && x.bhp > 0;
                if (thp || bhp) h << " 1* 1* " << (thp ? fmt(x.thp) : std::string("1*")) << " " << (bhp ? fmt(x.bhp) : std::string("1*"));
                h << " /\n";
            }
        } else {
            const bool gas = x.ginj > 0 && x.winj <= 0;
            if (x.shut) ij << "  " << quoteWell(x.well) << " WATER SHUT 0 /\n";
            else {
                ij << "  " << quoteWell(x.well) << (gas ? " GAS" : " WATER") << " OPEN " << fmt(gas ? x.ginj : x.winj);
                const bool thp = o.writeThp && x.thp > 0, bhp = o.writeBhp && x.bhp > 0;
                if (thp || bhp) ij << " " << (bhp ? fmt(x.bhp) : std::string("1*")) << " " << (thp ? fmt(x.thp) : std::string("1*"));
                ij << " /\n";
            }
        }
        double prev = 1;
        if (auto it = prevEf.find(x.well); it != prevEf.end()) prev = it->second.second;
        if (!x.shut && std::fabs(x.wefac - prev) > 1e-4) ef << "  " << quoteWell(x.well) << " " << fmt(x.wefac) << " /\n";
    }
    std::string s;
    if (!h.str().empty()) s += "WCONHIST\n" + h.str() + "/\n";
    if (!ij.str().empty()) s += "WCONINJH\n" + ij.str() + "/\n";
    if (!ef.str().empty()) s += "WEFAC\n" + ef.str() + "/\n";
    return s;
}

namespace {
bool findStart(const std::vector<std::string>& L, int& y, int& m, int& d) {
    for (size_t i = 0; i < L.size(); ++i) {
        std::string kw;
        if (!isKeyword(trim(stripComment(L[i])), kw) || kw != "START") continue;
        for (size_t j = i + 1; j < L.size() && j < i + 10; ++j) {
            const std::string s = trim(stripComment(L[j]));
            if (s.empty()) continue;
            auto tk = splitWs(s.substr(0, s.find('/')));
            if (tk.size() < 3) return false;
            d = std::atoi(unquote(tk[0]).c_str()); m = monthOf(tk[1]); y = std::atoi(unquote(tk[2]).c_str());
            return m > 0 && y > 0 && d > 0;
        }
    }
    return false;
}
}  // namespace

std::string mergeIntoDeck(const std::string& deck, const Table& t, const Options& o, std::vector<Issue>& issues,
                          const std::vector<std::string>& deckWells) {
    auto issue = [&](int sev, const std::string& well, const std::string& msg) { Issue i; i.severity = sev; i.well = well; i.msg = msg; issues.push_back(i); };
    std::vector<std::string> L;
    { std::istringstream in(deck); std::string l; while (std::getline(in, l)) { if (!l.empty() && l.back() == '\r') l.pop_back(); L.push_back(l); } }
    int sy = 0, sm = 0, sd = 0;
    if (!findStart(L, sy, sm, sd)) {
        if (o.startYear <= 0) { issue(2, "", "no START date in the deck - nothing merged"); return deck; }
        sy = o.startYear; sm = o.startMonth; sd = o.startDay;
        issue(0, "", "START not found in the deck text - " + std::to_string(sd) + "/" + std::to_string(sm) + "/" + std::to_string(sy) + " used");
    }
    const long day0 = daysFromCivil(sy, sm, sd);
    if (t.periods.empty()) { issue(2, "", "no history periods"); return deck; }
    std::set<std::string> imported(t.wells.begin(), t.wells.end());
    if (!deckWells.empty()) {
        std::set<std::string> dw(deckWells.begin(), deckWells.end());
        for (const auto& w : t.wells) if (!dw.count(w)) issue(1, w, "well not defined in the deck (WELSPECS) - its history keywords would fail; add the well first");
        for (const auto& w : deckWells) if (!imported.count(w)) issue(0, w, "deck well without data in the table - its deck controls are kept");
    }
    {   // C38 unit plausibility against the deck unit system: FIELD gas is Mscf (a table in scf gives GOR x1000)
        bool field = false, metric = false;
        for (const auto& l : L) { std::string kw; if (isKeyword(trim(stripComment(l)), kw)) { if (kw == "FIELD") field = true; else if (kw == "METRIC") metric = true; else if (kw == "GRID") break; } }
        double gorMax = 0; std::string gw;
        for (const auto& w : t.wp) if (!w.inj && w.orat > 1e-9 && w.grat / w.orat > gorMax) { gorMax = w.grat / w.orat; gw = w.well; }
        const double lim = field ? 100.0 : metric ? 18000.0 : 0.0;   // Mscf/stb resp. sm3/sm3 - beyond gas condensates
        if (lim > 0 && gorMax > lim)
            issue(1, gw, "GOR " + fmt(gorMax) + (field ? " Mscf/stb" : " sm3/sm3") + " is implausible for the " + (field ? "FIELD" : "METRIC") +
                         " deck - check the gas unit" + (field ? " (FIELD rates are Mscf: scf data need gas factor 0.001)" : " (thousand sm3 data need gas factor 1000)"));
    }
    // ---- locate SCHEDULE
    int sch = -1;
    for (size_t i = 0; i < L.size(); ++i) { std::string kw; if (isKeyword(trim(stripComment(L[i])), kw) && kw == "SCHEDULE") { sch = int(i); break; } }
    std::string prefix, tail;
    struct Block { double t = 0; std::string text; };
    std::vector<Block> blocks(1);
    size_t i0;
    if (sch < 0) {
        issue(1, "", "no SCHEDULE section - one is appended");
        size_t endAt = L.size();
        for (size_t i = 0; i < L.size(); ++i) { std::string kw; if (isKeyword(trim(stripComment(L[i])), kw) && kw == "END") { endAt = i; break; } }
        for (size_t i = 0; i < endAt; ++i) prefix += L[i] + "\n";
        prefix += "\nSCHEDULE\n";
        for (size_t i = endAt; i < L.size(); ++i) tail += L[i] + "\n";
        i0 = L.size();
    } else {
        for (int i = 0; i <= sch; ++i) prefix += L[size_t(i)] + "\n";
        i0 = size_t(sch) + 1;
    }
    enum Mode { NORMAL, DATES, TSTEP, CTRL } mode = NORMAL;
    static const std::set<std::string> ctrlKw = {"WCONPROD", "WCONHIST", "WCONINJE", "WCONINJH"};
    std::string ctrlName; std::vector<std::string> ctrlLines;
    int nDropped = 0; bool sawInclude = false, unsorted = false;
    auto flushCtrl = [&]() {
        std::vector<std::string> kept; int nrec = 0;
        for (const auto& l : ctrlLines) {
            const std::string s = trim(stripComment(l));
            if (s.empty()) { kept.push_back(l); continue; }
            const auto tk = splitWs(s);
            if (!tk.empty() && imported.count(unquote(tk[0]))) { ++nDropped; continue; }
            kept.push_back(l); ++nrec;
        }
        if (nrec > 0) { blocks.back().text += ctrlName + "\n"; for (const auto& l : kept) blocks.back().text += l + "\n"; blocks.back().text += "/\n"; }
        ctrlLines.clear();
    };
    auto addBlock = [&](double tt) {
        if (tt < blocks.back().t - 1e-9) unsorted = true;
        Block b; b.t = tt; blocks.push_back(b);
    };
    for (size_t i = i0; i < L.size(); ++i) {
        const std::string& line = L[i];
        const std::string s = trim(stripComment(line));
        std::string kw;
        const bool isKw = isKeyword(s, kw);
        if (mode == CTRL) {
            if (s == "/") { flushCtrl(); mode = NORMAL; continue; }
            if (!isKw) { ctrlLines.push_back(line); continue; }
            flushCtrl(); mode = NORMAL;   // keyword without terminator: fall through
        }
        if (mode == DATES) {
            if (s.empty()) continue;
            if (s == "/") { mode = NORMAL; continue; }
            if (!isKw) {
                auto tk = splitWs(s.substr(0, s.find('/')));
                if (tk.size() >= 3) {
                    const int d = std::atoi(unquote(tk[0]).c_str()), m = monthOf(tk[1]), y = std::atoi(unquote(tk[2]).c_str());
                    double secs = 0;
                    if (tk.size() >= 4) { int hh = 0, mm = 0, ss = 0; if (std::sscanf(unquote(tk[3]).c_str(), "%d:%d:%d", &hh, &mm, &ss) >= 2) secs = hh * 3600.0 + mm * 60.0 + ss; }
                    if (m > 0) addBlock(double(daysFromCivil(y, m, d) - day0) + secs / 86400.0);
                }
                continue;
            }
            mode = NORMAL;
        }
        if (mode == TSTEP) {
            if (!isKw) {
                const size_t sl = s.find('/');
                for (const auto& tok : splitWs(s.substr(0, sl))) {
                    const size_t st = tok.find('*');
                    int rep = 1; double v = 0;
                    if (st != std::string::npos) { rep = std::atoi(tok.substr(0, st).c_str()); v = std::atof(tok.substr(st + 1).c_str()); }
                    else v = std::atof(tok.c_str());
                    for (int r = 0; r < rep; ++r) addBlock(blocks.back().t + v);
                }
                if (sl != std::string::npos) mode = NORMAL;
                continue;
            }
            mode = NORMAL;
        }
        if (isKw) {
            if (kw == "DATES") { mode = DATES; continue; }
            if (kw == "TSTEP") { mode = TSTEP; continue; }
            if (kw == "END") { for (size_t j = i; j < L.size(); ++j) tail += L[j] + "\n"; break; }
            if (kw == "INCLUDE") sawInclude = true;
            if (o.replaceControls && ctrlKw.count(kw)) { mode = CTRL; ctrlName = kw; ctrlLines.clear(); continue; }
        }
        blocks.back().text += line + "\n";
    }
    if (mode == CTRL) flushCtrl();
    if (unsorted) { issue(2, "", "the deck's DATES / TSTEP are not increasing - nothing merged"); return deck; }
    if (sawInclude) issue(1, "", "the SCHEDULE has INCLUDE files - DATES / TSTEP inside them are not seen by the merge");
    if (nDropped > 0) issue(0, "", std::to_string(nDropped) + " existing WCONPROD / WCONHIST / WCONINJE / WCONINJH record(s) of imported wells removed");
    // ---- history events
    int nIns = 0, nApp = 0;
    for (int p = 0; p < int(t.periods.size()); ++p) {
        double tp = double(t.periods[size_t(p)].d0 - day0);
        if (tp < 0) {
            if (t.periods[size_t(p)].d1 - day0 <= 0) { issue(1, "", "period " + isoDate(t.periods[size_t(p)].d0) + " ends before START - skipped"); continue; }
            issue(1, "", "period " + isoDate(t.periods[size_t(p)].d0) + " starts before START - applied at START");
            tp = 0;
        }
        const std::string kws = periodKeywords(t, o, p);
        if (kws.empty()) continue;
        const std::string text = "-- history " + isoDate(t.periods[size_t(p)].d0) + " .. " + isoDate(t.periods[size_t(p)].d1) + " (production-data import)\n" + kws;
        size_t at = blocks.size();
        for (size_t b = 0; b < blocks.size(); ++b) if (std::fabs(blocks[b].t - tp) < 1e-6) { at = b; break; }
        if (at < blocks.size()) { blocks[at].text += text; ++nApp; continue; }
        size_t pos = 0; while (pos < blocks.size() && blocks[pos].t < tp) ++pos;
        Block nb; nb.t = tp; nb.text = text; blocks.insert(blocks.begin() + long(pos), nb); ++nIns;
    }
    const double tEnd = double(t.periods.back().d1 - day0);
    if (tEnd > blocks.back().t + 1e-6) { Block nb; nb.t = tEnd; blocks.push_back(nb); }
    issue(0, "", "history merged: " + std::to_string(nApp) + " period(s) at existing report dates, " + std::to_string(nIns) + " new report date(s); schedule ends " +
          isoDate(day0 + long(std::floor(blocks.back().t + 1e-9))));
    // ---- emit
    std::string out = prefix + blocks[0].text;
    std::vector<long> pendingDates;
    auto flushDates = [&]() {
        if (pendingDates.empty()) return;
        out += "DATES\n"; for (long d : pendingDates) out += "  " + dateString(d) + " /\n"; out += "/\n";
        pendingDates.clear();
    };
    for (size_t b = 1; b < blocks.size(); ++b) {
        const double tt = blocks[b].t;
        const bool whole = std::fabs(tt - std::round(tt)) < 1e-6;
        if (whole) pendingDates.push_back(day0 + long(std::llround(tt)));
        else { flushDates(); out += "TSTEP\n  " + fmt(tt - blocks[b - 1].t) + " /\n"; }
        if (trim(blocks[b].text).empty()) continue;
        flushDates();
        out += blocks[b].text;
    }
    flushDates();
    out += tail;
    return out;
}

std::string historyInclude(const Table& t, const Options& o) {
    std::string s = "-- production history generated by ResSim Studio (production-data import)\n";
    for (int p = 0; p < int(t.periods.size()); ++p) {
        const std::string k = periodKeywords(t, o, p);
        const long d0 = t.periods[size_t(p)].d0;
        const bool atStart = o.startYear > 0 && d0 == daysFromCivil(o.startYear, o.startMonth, o.startDay);
        if (!atStart) s += "DATES\n  " + dateString(d0) + " /\n/\n";
        s += k;
    }
    s += "DATES\n  " + dateString(t.periods.back().d1) + " /\n/\n";
    return s;
}

std::string observedCsv(const Table& t, const Options& o) {
    std::string s = "WELL,DATE,DAYS,WOPRH,WWPRH,WGPRH,WWIRH,WGIRH,WBHPH\n";
    const long day0 = o.startYear > 0 ? daysFromCivil(o.startYear, o.startMonth, o.startDay) : t.periods.front().d0;
    for (const auto& x : t.wp) {
        const Period& p = t.periods[size_t(x.period)];
        const double f = x.shut ? 0 : x.wefac;
        s += x.well + "," + isoDate(p.d1) + "," + fmt(double(p.d1 - day0)) + "," + fmt(x.orat * f) + "," + fmt(x.wrat * f) + "," + fmt(x.grat * f) + "," +
             fmt(x.winj * f) + "," + fmt(x.ginj * f) + "," + (x.bhp > 0 ? fmt(x.bhp) : std::string()) + "\n";
    }
    return s;
}

}  // namespace pdata
}  // namespace ressim
