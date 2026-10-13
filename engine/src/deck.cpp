// ECLIPSE 100 deck import and export (two-way).
// Import = generic deck reader (deck_parser.cpp, OPM keyword database) + the semantic mapping below.
#include "ressim/io.h"
#include "ressim/grid.h"
#include "ressim/simulator.h"
#include "ressim/deck_parser.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cctype>
#include <map>
#include <set>
#include <algorithm>
#include <limits>
#include <functional>
#include <memory>
#include <cstdio>
#include <numeric>

namespace ressim {

namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();

std::string num(double v) {
    // shortest representation that reads back to exactly the same double (never merges distinct table values)
    char b[40];
    for (int p = 10; p <= 17; ++p) {
        std::snprintf(b, sizeof(b), "%.*g", p, v);
        if (std::strtod(b, nullptr) == v) break;
    }
    return b;
}
std::string fmtNum(double v) { return num(v); }

const std::set<std::string> kDimsKeywords = {"WELLDIMS", "TABDIMS", "EQLDIMS", "REGDIMS", "AQUDIMS", "FAULTDIM"};
const char* const kDimsNote = "minimum dimensions from the imported deck (written as max(this, model needs))";
bool isDimsMinimum(const ExtraKeyword& e) { return e.section == DeckSection::Runspec && kDimsKeywords.count(e.name) && e.note == kDimsNote; }

// items of a one-record keyword body ("130 36 15 84 /", "1* 2*5 /"): defaulted items are 0
std::vector<long> dimsItems(const std::string& body) {
    std::vector<long> v;
    std::istringstream ss(body);
    std::string t;
    while (ss >> t) {
        if (t.rfind("--", 0) == 0) { std::string rest; std::getline(ss, rest); continue; }
        if (t == "/" ) break;
        bool end = false;
        if (t.back() == '/') { t.pop_back(); end = true; if (t.empty()) break; }
        const size_t star = t.find('*');
        if (star != std::string::npos) {
            const long n = star == 0 ? 1 : std::max(1L, std::atol(t.substr(0, star).c_str()));
            const long val = star + 1 < t.size() ? std::atol(t.substr(star + 1).c_str()) : 0;
            for (long i = 0; i < n && v.size() < 64; ++i) v.push_back(val);
        } else v.push_back(std::atol(t.c_str()));
        if (end) break;
    }
    return v;
}

// ECLIPSE reads at most 132 characters per line (the rest is silently dropped -> "VALUES NOT DEFINED IN BOX"): arrays are
// wrapped by width (<= 80 columns) as well as by count.
template <class T>
void writeArray(std::ostream& o, const std::vector<T>& v, int perLine = 8) {
    size_t i = 0; int col = 0; size_t width = 0;
    while (i < v.size()) {
        size_t j = i + 1;
        while (j < v.size() && v[j] == v[i]) ++j;
        const size_t cnt = j - i;
        const std::string item = cnt > 1 ? std::to_string(cnt) + "*" + num(double(v[i])) : num(double(v[i]));
        if (col > 0 && (col % perLine == 0 || width + 1 + item.size() > 80)) { o << "\n"; col = 0; width = 0; }
        o << ' ' << item;
        width += 1 + item.size();
        ++col;
        i = j;
    }
    if (col > 0) o << "\n";
}

// Safety net for the whole deck: lines longer than 'maxLen' are broken at blanks outside quotes (a record may span lines);
// over-long comment lines are shortened. Keeps decks readable by ECLIPSE (132-column limit) and OPM Flow.
std::string wrapLongLines(const std::string& text, size_t maxLen = 128) {
    std::string out;
    out.reserve(text.size() + text.size() / 64);
    size_t a = 0;
    while (a < text.size()) {
        size_t b = text.find('\n', a);
        if (b == std::string::npos) b = text.size();
        std::string line = text.substr(a, b - a);
        a = b + 1;
        bool cr = !line.empty() && line.back() == '\r';
        if (cr) line.pop_back();
        while (line.size() > maxLen) {
            size_t lead = line.find_first_not_of(" \t");
            if (lead != std::string::npos && line.compare(lead, 2, "--") == 0) { line.resize(maxLen); break; }
            // last blank outside quotes before maxLen, and before any trailing "--" comment
            size_t cut = std::string::npos;
            bool q = false;
            for (size_t k = 0; k < maxLen && k < line.size(); ++k) {
                const char ch = line[k];
                if (ch == '\'') q = !q;
                else if (!q && ch == '-' && k + 1 < line.size() && line[k + 1] == '-') break;
                else if (!q && (ch == ' ' || ch == '\t') && k > 0) cut = k;
            }
            if (cut == std::string::npos || cut == 0) break;    // nothing safe to split (one huge token)
            out.append(line, 0, cut);
            out += cr ? "\r\n" : "\n";
            size_t rest = line.find_first_not_of(" \t", cut);
            line = rest == std::string::npos ? std::string() : " " + line.substr(rest);
        }
        out += line;
        if (b < text.size() || !line.empty()) out += cr ? "\r\n" : "\n";
    }
    if (!text.empty() && text.back() != '\n' && !out.empty() && out.back() == '\n') out.pop_back();
    return out;
}

std::string upper(std::string s) { for (char& ch : s) ch = char(std::toupper((unsigned char)ch)); return s; }

bool matchName(const std::string& pat, const std::string& name) {
    if (!pat.empty() && pat.back() == '*') return name.compare(0, pat.size() - 1, pat, 0, pat.size() - 1) == 0;
    return pat == name;
}

const char* kMonths[12] = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};
long daysFromCivil(int y, int m, int d) {           // Howard Hinnant's algorithm
    y -= m <= 2;
    const long era = (y >= 0 ? y : y - 399) / 400;
    const unsigned yoe = unsigned(y - era * 400);
    const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
    const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return era * 146097 + long(doe) - 719468;
}
int monthFrom(const std::string& s) {
    const std::string u = upper(s);
    if (u == "JLY") return 7;
    for (int m = 0; m < 12; ++m) if (u == kMonths[m]) return m + 1;
    try { const int v = std::stoi(u); if (v >= 1 && v <= 12) return v; } catch (...) {}
    throw std::runtime_error("Unknown month '" + s + "'");
}

double lin(const std::vector<double>& x, const std::vector<double>& y, double v) {
    if (x.empty()) return 0;
    if (v <= x.front()) return y.front();
    if (v >= x.back()) return y.back();
    const size_t i = size_t(std::upper_bound(x.begin(), x.end(), v) - x.begin()) - 1;
    const double dx = x[i + 1] - x[i];
    return dx > 0 ? y[i] + (v - x[i]) * (y[i + 1] - y[i]) / dx : y[i];
}

// ======================================================================== keyword support table
using deck::Support;
struct SupportInfo { Support s; const char* note; };

const std::map<std::string, SupportInfo>& supportTable() {
    static const std::map<std::string, SupportInfo> t = [] {
        std::map<std::string, SupportInfo> m;
        auto add = [&](Support s, const char* note, const char* list) {
            std::istringstream is(list); std::string k;
            while (is >> k) m[k] = SupportInfo{s, note};
        };
        add(Support::Full, "", "RUNSPEC GRID EDIT PROPS REGIONS SOLUTION SUMMARY SCHEDULE END TITLE DIMENS OIL WATER GAS DISGAS "
            "FIELD METRIC START TABDIMS EQLDIMS REGDIMS INCLUDE PATHS SKIP SKIP100 SKIP300 ENDSKIP "
            "DX DY DZ DXV DYV DZV TOPS PORO PERMX PERMY PERMZ NTG ACTNUM MULTX MULTY MULTZ MULTPV COORD ZCORN SPECGRID "
            "BOX ENDBOX EQUALS MULTIPLY ADD COPY MINVALUE MAXVALUE EQUALREG MULTIREG ADDREG COPYREG MULTNUM FLUXNUM OPERNUM "
            "FAULTS MULTFLT NNC MINPV MINPORV PVTW ROCK DENSITY GRAVITY PVDO PVDG PVTO PVCDO SWOF SGOF SWFN SGFN SOF3 SOF2 "
            "SLGOF SATNUM PVTNUM EQLNUM FIPNUM EQUIL RSVD PBVD PRESSURE SWAT SGAS RS AQUFETP AQUCT AQUANCON "
            "SWATINIT THPRES VAPOIL PVTG RVVD RV WELSPECS COMPDAT WCONPROD WCONHIST WCONINJE WCONINJH WELOPEN WPIMULT WELTARG TSTEP DATES TUNING ROCKOPTS NEWTRAN OLDTRAN "
            "ENDSCALE SWL SWCR SWU SGL SGCR SGU SOWCR SOGCR "
            "POLYMER PLYVISC PLYADS PLYROCK PLYMAX PLMIXPAR WPOLYMER SPOLY PLYSHLOG PLMIXNUM STONE STONE1 STONE2 STONE1EX DRSDT "
            "DUALPORO DUALPERM SIGMA SIGMAV LTOSIGMA LX LY LZ MULTMF NODPPM DPGRID PINCH");
        add(Support::Partial, "SHRATE (shear rate based polymer shear model) is not modelled - the velocity form of PLYSHLOG is used", "SHRATE");
        add(Support::Partial, "applied to neighbour connections only (not to NNCs)", "MULTREGT");
        add(Support::Partial, "multisegment data ignored - treated as a standard well", "WELSEGS COMPSEGS");
        add(Support::Partial, "group control is not simulated", "GRUPTREE");
        add(Support::Ignored, "output / dimension keyword - no effect on the results",
            "UNIFOUT UNIFIN FMTOUT FMTIN MULTOUT MULTIN INIT NOECHO ECHO NOWARN WARN MESSAGES MONITOR NOMONITO NOSIM "
            "GRIDFILE MAPAXES MAPUNITS GRIDUNIT NOINSPEC NORSSPEC RPTGRID RPTPROPS RPTREGS RPTSOL RPTSCHED "
            "RPTRST RPTRUNSP RPTSMRY RPTONLY RPTONLYO OUTSOL WELLDIMS AQUDIMS NSTACK FAULTDIM NUPCOL VFPPDIMS VFPIDIMS "
            "SMRYDIMS WSEGDIMS UDQDIMS UDADIMS ACTDIMS NETWORK COORDSYS NONNC OPTIONS OPTIONS3 SAVE EXTRAPMS DEBUG DEBUG3 "
            "GRIDOPTS SATOPTS EQLOPTS LIVEOIL FULLIMP IMPES CPR NOHYST ZIPPY2 FORMFEED ALL MSUMLINS MSUMNEWT SEPARATE "
            "EXCEL RUNSUM NARROW PERFORMA DATE NEXTSTEP NEXT RPTPRINT PRINTREP DIMPES DIFFDP UNCODHMD MEMORY NOINSPEC "
            "INSPEC RSSPEC ECHO PETOPTS PARALLEL WRFT WRFTPLT COMPORD TRACERS GDORIENT NOGGF EHYSTR UDQPARAM "
            "EXIT SMRYDIMS");
        return m;
    }();
    return t;
}

// summary vectors the exporter always writes (not kept as extra keywords)
static const std::set<std::string>& exporterSummary() {
    static const std::set<std::string> s = {"FOPR", "FWPR", "FGPR", "FLPR", "FWIR", "FGIR", "FOPT", "FWPT", "FGPT", "FWIT", "FGIT",
        "FPR", "FPRP", "FWCT", "FGOR", "FOIP", "FWIP", "FGIP", "FAQR", "FAQT", "FCIR", "FCPR", "FCIT", "FCPT", "FCIP", "FCAD", "WCIR", "WCPR", "WCIT", "WCPT", "WBHP", "WOPR", "WWPR", "WGPR", "WLPR", "WWIR", "WGIR",
        "WWCT", "WGOR", "WOPT", "WWPT", "WGPT", "WWIT", "WGIT"};
    return s;
}

// ======================================================================== deck -> Case
using deck::Keyword;
using deck::Record;

struct Importer {
    const deck::RawDeck& D;
    std::vector<DeckMessage> msgs;
    Case c;
    const Keyword* K = nullptr;
    std::string section = "RUNSPEC";
    int ntpvt = 1, ntsfun = 1, ntequl = 1;
    int n = 0, nl = 0;
    struct Box { int i1, i2, j1, j2, k1, k2; } box{0, 0, 0, 0, 0, 0};
    bool boxActive = false;
    std::map<std::string, std::vector<double>> arr;
    std::set<std::string> usedArrays;
    std::map<int, std::vector<AquiferConnection>> aqConn;
    std::map<int, double> aqVisc;
    std::vector<WellControl> pending;
    std::vector<WellPolymer> pendingPoly;
    double pendingDrsdt = -2; bool pendingDrsdtFree = false;
    std::map<std::string, WellControl> lastCtrl;
    std::map<std::string, bool> lastInj, pendingInj;     // control came from an injector keyword
    std::set<std::pair<size_t, std::string>> injCtrl;    // (report step, well) of injector controls
    long dayNow = 0;
    double daysElapsed = 0;
    bool anyStep = false;
    // tables (per region)
    std::vector<PvtRegion> pvt;
    std::vector<char> hasPvtw, hasRock, hasDensity, hasPvdo, hasPvdg, hasPvto;
    std::vector<std::vector<double>> swof, sgof, swfn, sgfn, sof3, sof2, slgof;   // [table] -> flat rows
    std::vector<EquilInput> equil;
    std::vector<std::vector<double>> rsvd, pbvd;
    std::vector<int> equilItem7, equilItem8;
    std::vector<std::vector<double>> rvvd;
    bool sawEquil = false;
    struct MultRegT { int r1, r2; double m; std::string dirs, set; };
    std::vector<MultRegT> multregt;
    // grouped messages
    struct Group { int level; std::string keyword, file; int line; std::vector<std::string> items; int count = 0; };
    std::map<std::string, Group> groups;
    std::vector<std::string> groupOrder;

    explicit Importer(const deck::RawDeck& d) : D(d) {
        c = Case();
        c.title = "Imported deck";
        c.units = UnitSystem::Metric;            // ECLIPSE default
        c.phases.water = c.phases.oil = c.phases.gas = c.phases.disgas = false;
        c.props.pvdo.clear(); c.props.pvdg.clear(); c.props.pvto.clear(); c.props.swof.clear(); c.props.sgof.clear();
        c.wells.clear(); c.schedule.clear(); c.aquifers.clear();
        dayNow = daysFromCivil(c.startYear, c.startMonth, c.startDay);
    }

    // ------------------------------------------------------------------ messages
    std::string fileOf(const Keyword* k) const {
        if (!k || k->file < 0 || k->file >= int(D.files.size())) return D.files.empty() ? "deck" : D.files[0];
        return D.files[size_t(k->file)];
    }
    std::string where(int line = 0) const { return K ? D.where(K->file, line > 0 ? line : K->line) : std::string("deck"); }
    void msg(int level, const std::string& text, int line = 0) {
        DeckMessage m; m.level = level; m.keyword = K ? K->name : "Deck"; m.file = fileOf(K);
        m.line = line > 0 ? line : (K ? K->line : 0); m.text = text;
        msgs.push_back(m);
    }
    void warn(const std::string& text, int line = 0) { msg(1, text, line); }
    void info(const std::string& text, int line = 0) { msg(0, text, line); }
    // one message per (keyword, text); the items (well names, ...) are listed
    // a single SCHEDULE record ResSim cannot represent: kept verbatim as a one-record extra keyword at the current step
    void keepRecord(const Keyword& k, const Record& r, const std::string& note) {
        Keyword one; one.name = k.name; one.records.push_back(r);
        ExtraKeyword e; e.name = k.name; e.section = DeckSection::Schedule; e.step = int(c.schedule.size()); e.note = note;
        e.beforeControls = false;     // replaces / complements a control record
        e.body = deck::keywordBodyText(one);
        if (!c.extraKeywords.empty()) {    // merge with the previous one-record keyword of the same name and step
            ExtraKeyword& p = c.extraKeywords.back();
            const std::string tail = "/\n/\n";
            if (p.name == e.name && p.step == e.step && p.note == e.note && p.body.size() >= tail.size() &&
                p.body.compare(p.body.size() - tail.size(), tail.size(), tail) == 0) {
                p.body.erase(p.body.size() - 2); p.body += e.body; return;
            }
        }
        c.extraKeywords.push_back(std::move(e));
    }
    void grouped(int level, const std::string& text, const std::string& item = "") {
        const std::string key = (K ? K->name : "Deck") + "\x1f" + text;
        auto it = groups.find(key);
        if (it == groups.end()) {
            Group g; g.level = level; g.keyword = K ? K->name : "Deck"; g.file = fileOf(K); g.line = K ? K->line : 0;
            it = groups.emplace(key, g).first; groupOrder.push_back(key);
        }
        Group& g = it->second;
        ++g.count;
        if (!item.empty() && std::find(g.items.begin(), g.items.end(), item) == g.items.end()) g.items.push_back(item);
    }
    [[noreturn]] void fail(const std::string& text, int line = 0) const {
        throw std::runtime_error((K ? K->name : std::string("Deck")) + " (" + where(line) + "): " + text);
    }

    // ------------------------------------------------------------------ item access (UDA-tolerant)
    double num(const Record& r, size_t i, double def) {
        if (!r.has(i)) return def;
        if (r.isNumber(i)) return r.num(i, def);
        grouped(1, "non-numeric value in a numeric item (UDA/UDQ?) - default used", r.items[i].s);
        return def;
    }
    int inum(const Record& r, size_t i, int def) { return int(std::lround(num(r, i, double(def)))); }
    static std::string str(const Record& r, size_t i, const std::string& def = "") { return r.str(i, def); }
    std::vector<double> numbers(const Record& r) {
        std::vector<double> v(r.items.size());
        for (size_t i = 0; i < r.items.size(); ++i) {
            if (r.items[i].def) { v[i] = NaN; continue; }
            if (!r.isNumber(i)) fail("expected a number but found '" + r.items[i].s + "'", r.line);
            v[i] = r.num(i, NaN);
        }
        return v;
    }

    // ------------------------------------------------------------------ grid arrays
    Box fullBox() const { return {0, c.grid.nx - 1, 0, c.grid.ny - 1, 0, c.grid.nz - 1}; }
    Box curBox() const { return boxActive ? box : fullBox(); }
    void needDims() const { if (n == 0) fail("DIMENS must come before grid data"); }
    static double arrayDefault(const std::string& name) {
        if (name == "NTG" || name == "MULTPV" || name == "ACTNUM" || name.compare(0, 4, "MULT") == 0) return 1.0;
        if (name.size() > 3 && name.compare(name.size() - 3, 3, "NUM") == 0) return 1.0;
        return NaN;
    }
    std::vector<double>& array(const std::string& name) {
        needDims();
        auto& v = arr[name];
        if (v.empty()) v.assign(size_t(n), arrayDefault(name));
        return v;
    }
    Box readBox(const Record& r, size_t first, const Box& prev) {
        Box b = prev;
        const Box f = fullBox();
        auto get = [&](size_t i, int def) { return r.has(i) ? inum(r, i, 0) - 1 : def; };
        b.i1 = get(first, b.i1); b.i2 = get(first + 1, b.i2); b.j1 = get(first + 2, b.j1);
        b.j2 = get(first + 3, b.j2); b.k1 = get(first + 4, b.k1); b.k2 = get(first + 5, b.k2);
        if (b.i1 < 0 || b.i2 > f.i2 || b.j1 < 0 || b.j2 > f.j2 || b.k1 < 0 || b.k2 > f.k2 || b.i1 > b.i2 || b.j1 > b.j2 || b.k1 > b.k2)
            fail(fmtBox(b) + " is outside the " + std::to_string(c.grid.nx) + "x" + std::to_string(c.grid.ny) + "x" +
                 std::to_string(c.grid.nz) + " grid", r.line);
        return b;
    }
    static std::string fmtBox(const Box& b) {
        return "box " + std::to_string(b.i1 + 1) + "-" + std::to_string(b.i2 + 1) + ", " + std::to_string(b.j1 + 1) + "-" +
               std::to_string(b.j2 + 1) + ", " + std::to_string(b.k1 + 1) + "-" + std::to_string(b.k2 + 1);
    }
    template <class F> void forBox(const Box& b, F f) {
        for (int k = b.k1; k <= b.k2; ++k)
            for (int j = b.j1; j <= b.j2; ++j)
                for (int i = b.i1; i <= b.i2; ++i) f(c.grid.index(i, j, k));
    }

    void dataArray(const Keyword& k) {
        const std::string& name = k.name;
        const std::vector<double>& v = k.values;
        if (name == "MAPAXES") return;
        if (section == "SCHEDULE") { grouped(1, "grid property changes in SCHEDULE are not supported - ignored"); return; }
        needDims();
        if (name == "DXV" || name == "DYV" || name == "DZV") {
            const size_t want = size_t(name == "DXV" ? c.grid.nx : name == "DYV" ? c.grid.ny : c.grid.nz);
            if (v.size() != want) fail("expected " + std::to_string(want) + " values, found " + std::to_string(v.size()));
            auto& a = array(name == "DXV" ? "DX" : name == "DYV" ? "DY" : "DZ");
            for (int kk = 0; kk < c.grid.nz; ++kk)
                for (int j = 0; j < c.grid.ny; ++j)
                    for (int i = 0; i < c.grid.nx; ++i)
                        a[size_t(c.grid.index(i, j, kk))] = v[size_t(name == "DXV" ? i : name == "DYV" ? j : kk)];
            return;
        }
        if (name == "COORD") {
            const size_t want = size_t(6) * size_t(c.grid.nx + 1) * size_t(c.grid.ny + 1);
            if (v.size() != want) fail("COORD needs 6*(NX+1)*(NY+1) = " + std::to_string(want) + " values, found " + std::to_string(v.size()));
            c.grid.coord = v; return;
        }
        if (name == "ZCORN") {
            if (v.size() != size_t(8) * size_t(n)) fail("ZCORN needs 8*NX*NY*NZ = " + std::to_string(size_t(8) * n) + " values, found " + std::to_string(v.size()));
            c.grid.zcorn = v; return;
        }
        const Box b = curBox();
        auto& a = array(name);
        if (name == "TOPS" && !boxActive && v.size() == size_t(nl)) { for (int q = 0; q < nl; ++q) a[size_t(q)] = v[size_t(q)]; return; }
        const size_t bs = size_t(b.i2 - b.i1 + 1) * size_t(b.j2 - b.j1 + 1) * size_t(b.k2 - b.k1 + 1);
        if (v.size() < bs) warn("expected " + std::to_string(bs) + " values, found " + std::to_string(v.size()) + " - the remaining cells keep their previous / default value");
        else if (v.size() > bs) warn("expected " + std::to_string(bs) + " values, found " + std::to_string(v.size()) + " - extra values ignored");
        size_t m = 0;
        forBox(b, [&](int q) { if (m < v.size()) { if (!std::isnan(v[m])) a[size_t(q)] = v[m]; } ++m; });
    }

    // EQUALS / MULTIPLY / ADD / COPY / MINVALUE / MAXVALUE
    void operations(const Keyword& k) {
        needDims();
        const std::string& op = k.name;
        Box prev = curBox();
        for (const Record& r : k.records) {
            if (r.items.empty()) continue;
            if (op == "COPY") {
                const std::string src = upper(str(r, 0)), dst = upper(str(r, 1));
                const Box b = readBox(r, 2, prev); prev = b;
                if (!arr.count(src) && std::isnan(arrayDefault(src))) { warn("COPY source " + src + " is not defined - ignored", r.line); continue; }
                auto& s = array(src); auto& d = array(dst);
                forBox(b, [&](int q) { d[size_t(q)] = s[size_t(q)]; });
                continue;
            }
            std::string name = upper(str(r, 0));
            const double v = num(r, 1, NaN);
            if (std::isnan(v)) { warn("value missing for " + name + " - record ignored", r.line); continue; }
            const Box b = readBox(r, 2, prev); prev = b;
            if (name == "TRANX" || name == "TRANY" || name == "TRANZ") {
                // EDIT-section transmissibility edits: a multiplication equals a transmissibility multiplier
                if (op != "MULTIPLY") { grouped(1, op + " on " + name + " is not supported (only MULTIPLY) - ignored"); continue; }
                name = std::string("MULT") + name[4];
                if (!arr.count("TRAN-EDIT")) { arr["TRAN-EDIT"]; info("MULTIPLY on TRANX/TRANY/TRANZ applied as MULTX/MULTY/MULTZ"); }
            }
            if (op != "EQUALS" && !arr.count(name) && std::isnan(arrayDefault(name))) { grouped(1, op + " on an array before it is defined - ignored", name); continue; }
            auto& a = array(name);
            forBox(b, [&](int q) {
                double& x = a[size_t(q)];
                if (op == "EQUALS") x = v;
                else if (op == "MULTIPLY") x *= v;
                else if (op == "ADD") x += v;
                else if (op == "MINVALUE") x = std::max(x, v);
                else x = std::min(x, v);
            });
        }
    }

    // EQUALREG / MULTIREG / ADDREG / COPYREG
    void regionOperations(const Keyword& k) {
        needDims();
        const std::string& op = k.name;
        for (const Record& r : k.records) {
            if (r.items.empty()) continue;
            const bool copy = op == "COPYREG";
            const std::string name = upper(str(r, 0));
            const std::string src = copy ? name : "";
            const std::string dst = copy ? upper(str(r, 1)) : name;
            const double v = copy ? 0 : num(r, 1, NaN);
            const int reg = inum(r, 2, -1);
            const std::string set = upper(str(r, 3, "M"));
            const std::string regArr = set.empty() || set[0] == 'M' ? "MULTNUM" : set[0] == 'F' ? "FLUXNUM" : "OPERNUM";
            if (!copy && std::isnan(v)) { warn("value missing - record ignored", r.line); continue; }
            if (reg < 0) { warn("region number missing - record ignored", r.line); continue; }
            const std::vector<double>& R = array(regArr);
            if (copy && !arr.count(src)) { warn("COPYREG source " + src + " is not defined - ignored", r.line); continue; }
            auto& a = array(dst);
            const std::vector<double>* s = copy ? &array(src) : nullptr;
            for (int q = 0; q < n; ++q) {
                if (int(std::lround(R[size_t(q)])) != reg) continue;
                double& x = a[size_t(q)];
                if (copy) x = (*s)[size_t(q)];
                else if (op == "EQUALREG") x = v;
                else if (op == "MULTIREG") x *= v;
                else x += v;
            }
        }
    }

    // ------------------------------------------------------------------ tables
    // fill defaulted entries of columns 1.. by linear interpolation in column 0
    void fillDefaults(std::vector<double>& t, int ncol) {
        const size_t rows = t.size() / size_t(ncol);
        for (size_t r = 0; r < rows; ++r) if (std::isnan(t[r * ncol])) fail("the first column of a table cannot be defaulted");
        for (int col = 1; col < ncol; ++col)
            for (size_t r = 0; r < rows; ++r) {
                if (!std::isnan(t[r * ncol + col])) continue;
                long lo = long(r) - 1, hi = long(r) + 1;
                while (lo >= 0 && std::isnan(t[size_t(lo) * ncol + col])) --lo;
                while (hi < long(rows) && std::isnan(t[size_t(hi) * ncol + col])) ++hi;
                if (lo < 0 && hi >= long(rows)) fail("a table column is entirely defaulted");
                if (lo < 0) t[r * ncol + col] = t[size_t(hi) * ncol + col];
                else if (hi >= long(rows)) t[r * ncol + col] = t[size_t(lo) * ncol + col];
                else {
                    const double x0 = t[size_t(lo) * ncol], x1 = t[size_t(hi) * ncol], x = t[r * ncol];
                    const double y0 = t[size_t(lo) * ncol + col], y1 = t[size_t(hi) * ncol + col];
                    t[r * ncol + col] = x1 != x0 ? y0 + (x - x0) * (y1 - y0) / (x1 - x0) : y0;
                }
            }
    }
    // FromKw table keywords: one record per table; an empty record repeats the previous table
    std::vector<std::vector<double>> tables(const Keyword& k, int ncol, int ntab) {
        std::vector<std::vector<double>> out;
        for (const Record& r : k.records) {
            std::vector<double> v = numbers(r);
            if (v.empty()) {
                if (out.empty()) fail("the first table is empty");
                out.push_back(out.back()); continue;
            }
            if (v.size() % size_t(ncol)) fail("table " + std::to_string(out.size() + 1) + ": " + std::to_string(v.size()) +
                                              " values is not a multiple of " + std::to_string(ncol) + " columns", r.line);
            fillDefaults(v, ncol);
            out.push_back(std::move(v));
        }
        if (int(out.size()) < ntab && !out.empty()) {
            warn(std::to_string(out.size()) + " table(s) given but " + std::to_string(ntab) + " expected - the last one is repeated");
            while (int(out.size()) < ntab) out.push_back(out.back());
        }
        return out;
    }
    void ensurePvt(int nt) {
        if (int(pvt.size()) >= nt) return;
        pvt.resize(size_t(nt));
        for (auto* v : {&hasPvtw, &hasRock, &hasDensity, &hasPvdo, &hasPvdg, &hasPvto}) v->resize(size_t(nt), 0);
    }
    DensitySC defaultDensity() const {
        DensitySC d;
        if (c.units == UnitSystem::Field) { d.oil = 37.457; d.water = 62.366; d.gas = 0.062428; }
        else { d.oil = 600; d.water = 999.014; d.gas = 1.0; }
        return d;
    }
    // one record per PVT region; an empty record repeats the previous region
    template <class F> void perRegion(const Keyword& k, int nreg, F f) {
        ensurePvt(nreg);
        int t = 0;
        const Record* last = nullptr;
        for (const Record& r : k.records) {
            if (t >= nreg) break;
            const Record* use = r.items.empty() ? last : &r;
            if (!use) fail("the first record is empty");
            f(t, *use); last = use; ++t;
        }
        for (; t < nreg && last; ++t) f(t, *last);
    }

    void props(const Keyword& k) {
        const std::string& kw = k.name;
        if (kw == "PVTW") {
            perRegion(k, ntpvt, [&](int t, const Record& r) {
                auto& w = pvt[size_t(t)].pvtw;
                w.pref = num(r, 0, w.pref); w.bw = num(r, 1, 1.0); w.cw = num(r, 2, c.units == UnitSystem::Field ? 3e-6 : 4e-5);
                w.visc = num(r, 3, 0.5); w.cv = num(r, 4, 0); hasPvtw[size_t(t)] = 1;
            });
        } else if (kw == "ROCK") {
            if (k.records.size() > size_t(ntpvt)) warn("ROCK given per SATNUM/ROCKNUM region (ROCKOPTS) - the PVT-region mapping is used");
            perRegion(k, ntpvt, [&](int t, const Record& r) {
                auto& rk = pvt[size_t(t)].rock;
                rk.pref = num(r, 0, c.units == UnitSystem::Field ? 14.7 : 1.01325); rk.cr = num(r, 1, 0); hasRock[size_t(t)] = 1;
            });
        } else if (kw == "DENSITY") {
            const DensitySC def = defaultDensity();
            perRegion(k, ntpvt, [&](int t, const Record& r) {
                auto& d = pvt[size_t(t)].density;
                d.oil = num(r, 0, def.oil); d.water = num(r, 1, def.water); d.gas = num(r, 2, def.gas); hasDensity[size_t(t)] = 1;
            });
        } else if (kw == "GRAVITY") {
            const DensitySC def = defaultDensity();
            const double air = c.units == UnitSystem::Field ? 0.0763 : 1.2232;
            perRegion(k, ntpvt, [&](int t, const Record& r) {
                auto& d = pvt[size_t(t)].density;
                const double api = num(r, 0, 45.5), sgw = num(r, 1, 1.0), sgg = num(r, 2, 0.7773);
                d.oil = 141.5 / (api + 131.5) * def.water; d.water = sgw * def.water; d.gas = sgg * air; hasDensity[size_t(t)] = 1;
            });
            info("converted to DENSITY");
        } else if (kw == "PVCDO") {
            perRegion(k, ntpvt, [&](int t, const Record& r) {
                const double pref = num(r, 0, NaN), bo = num(r, 1, NaN), co = num(r, 2, 0), mu = num(r, 3, NaN), cv = num(r, 4, 0);
                if (std::isnan(pref) || std::isnan(bo) || std::isnan(mu)) fail("Pref, Bo and viscosity are required", r.line);
                auto& tab = pvt[size_t(t)].pvdo; tab.clear();
                for (int m = 0; m <= 12; ++m) {
                    const double p = pref * (0.05 + 0.25 * m);
                    tab.push_back({p, bo * std::exp(-co * (p - pref)), mu * std::exp(cv * (p - pref))});
                }
                hasPvdo[size_t(t)] = 1;
            });
            info("converted to a 13-row PVDO table per region");
        } else if (kw == "PVDO" || kw == "PVDG") {
            auto tabs = tables(k, 3, ntpvt);
            ensurePvt(int(tabs.size()));
            for (size_t t = 0; t < tabs.size() && t < size_t(ntpvt); ++t) {
                auto& tab = kw == "PVDO" ? pvt[t].pvdo : pvt[t].pvdg; tab.clear();
                for (size_t i = 0; i < tabs[t].size(); i += 3) tab.push_back({tabs[t][i], tabs[t][i + 1], tabs[t][i + 2]});
                (kw == "PVDO" ? hasPvdo : hasPvdg)[t] = 1;
            }
        } else if (kw == "PVTO" || kw == "PVTG") {
            // Tables: records, each table closed by a terminator record
            std::vector<std::vector<const Record*>> tabs(1);
            for (const Record& r : k.records) {
                if (r.terminator) { tabs.emplace_back(); continue; }
                if (!r.items.empty()) tabs.back().push_back(&r);
            }
            while (!tabs.empty() && tabs.back().empty() && int(tabs.size()) > ntpvt) tabs.pop_back();
            ensurePvt(ntpvt);
            for (size_t t = 0; t < size_t(ntpvt); ++t) {
                const bool have = t < tabs.size() && !tabs[t].empty();
                if (!have) {
                    if (t == 0) fail("the first table is empty");
                    if (kw == "PVTO") { pvt[t].pvto = pvt[t - 1].pvto; hasPvto[t] = 1; }
                    else { pvt[t].pvdg = pvt[t - 1].pvdg; pvt[t].pvtg = pvt[t - 1].pvtg; hasPvdg[t] = 1; }
                    continue;
                }
                if (kw == "PVTO") {
                    auto& out = pvt[t].pvto; out.clear();
                    for (const Record* r : tabs[t]) {
                        std::vector<double> v = numbers(*r);
                        if (v.size() < 4 || (v.size() - 1) % 3) fail("each record needs Rs followed by (P Bo mu) triplets", r->line);
                        for (double x : v) if (std::isnan(x)) fail("defaulted PVTO values are not supported", r->line);
                        PVTORecord rec; rec.rs = v[0]; rec.pb = v[1]; rec.bo = v[2]; rec.mu = v[3];
                        for (size_t i = 4; i + 2 < v.size(); i += 3) rec.under.push_back({v[i], v[i + 1], v[i + 2]});
                        out.push_back(rec);
                    }
                    hasPvto[t] = 1;
                } else {
                    auto& out = pvt[t].pvdg; out.clear();
                    auto& wet = pvt[t].pvtg; wet.clear();
                    for (const Record* r : tabs[t]) {
                        std::vector<double> v = numbers(*r);
                        if (v.size() < 4 || (v.size() - 1) % 3) fail("each record needs Pg followed by (Rv Bg mu) triplets", r->line);
                        for (double x : v) if (std::isnan(x)) fail("defaulted PVTG values are not supported", r->line);
                        out.push_back({v[0], v[2], v[3]});
                        PVTGRecord rec; rec.pg = v[0];
                        for (size_t i = 1; i + 2 < v.size(); i += 3) rec.rows.push_back({v[i], v[i + 1], v[i + 2]});
                        wet.push_back(rec);
                    }
                    hasPvdg[t] = 1;
                }
            }
            if (int(tabs.size()) > ntpvt && !(tabs.size() == size_t(ntpvt) + 1 && tabs.back().empty()))
                warn(std::to_string(tabs.size()) + " tables found but TABDIMS NTPVT = " + std::to_string(ntpvt));
        } else if (kw == "PLYVISC" || kw == "PLYADS") {
            auto t = tables(k, 2, kw == "PLYVISC" ? ntpvt : ntsfun);
            for (auto& v : t) for (size_t i = 2; i < v.size(); i += 2) if (!(v[i] > v[i - 2])) fail("the concentrations must increase");
            (kw == "PLYVISC" ? c.polymer.plyvisc : c.polymer.plyads) = t;
        } else if (kw == "STONE1") {
            c.props.kroModel = 1;
        } else if (kw == "STONE2" || kw == "STONE") {
            c.props.kroModel = 2;
        } else if (kw == "STONE1EX") {
            if (!k.records.empty()) c.props.stone1Exp = num(k.records[0], 0, 1.0);
            if (c.props.stone1Exp <= 0) fail("the exponent must be positive");
        } else if (kw == "PLYROCK") {
            c.polymer.plyrock.clear();
            const Record* last = nullptr;
            for (const Record& r : k.records) {
                if (int(c.polymer.plyrock.size()) >= ntsfun) break;
                const Record* u = r.items.empty() ? last : &r;
                if (!u) fail("the first record is empty");
                PolymerRock R;
                R.deadPoreVolume = num(*u, 0, 0); R.rrf = num(*u, 1, 1); R.rockDensity = num(*u, 2, 0);
                R.adsIndex = inum(*u, 3, 1) == 2 ? 2 : 1; R.maxAdsorption = num(*u, 4, 0);
                if (R.rrf < 1) fail("the residual resistance factor must be >= 1");
                c.polymer.plyrock.push_back(R); last = u;
            }
            while (last && int(c.polymer.plyrock.size()) < ntsfun) c.polymer.plyrock.push_back(c.polymer.plyrock.back());
        } else if (kw == "PLYMAX") {
            c.polymer.cmax.clear(); c.polymer.cmaxSalt.clear();
            for (const Record& r : k.records) {
                if (r.items.empty()) { if (c.polymer.cmax.empty()) continue; c.polymer.cmax.push_back(c.polymer.cmax.back()); c.polymer.cmaxSalt.push_back(c.polymer.cmaxSalt.back()); continue; }
                c.polymer.cmax.push_back(num(r, 0, 0)); c.polymer.cmaxSalt.push_back(num(r, 1, 0));
            }
            if (c.polymer.cmax.empty()) { c.polymer.cmax = {0.0}; c.polymer.cmaxSalt = {0.0}; }
        } else if (kw == "PLMIXPAR") {
            c.polymer.mixing.clear();
            for (const Record& r : k.records) {
                if (r.items.empty()) { if (!c.polymer.mixing.empty()) c.polymer.mixing.push_back(c.polymer.mixing.back()); continue; }
                c.polymer.mixing.push_back(num(r, 0, 1));
            }
            if (c.polymer.mixing.empty()) c.polymer.mixing = {1.0};
        } else if (kw == "PLYSHLOG") {
            c.polymer.plyshlog.clear();
            std::vector<const Record*> recs;
            for (const Record& r : k.records) recs.push_back(&r);
            for (size_t i = 0; i + 1 < recs.size(); i += 2) {
                PolymerInput::Shear S;
                S.refConc = num(*recs[i], 0, 0); S.refSalt = num(*recs[i], 1, -1); S.refTemp = num(*recs[i], 2, -1);
                S.rows = numbers(*recs[i + 1]);
                if (S.rows.size() < 4 || S.rows.size() % 2) fail("PLYSHLOG: the table needs at least 2 rows of (velocity, multiplier)", recs[i + 1]->line);
                for (size_t q = 2; q < S.rows.size(); q += 2) if (!(S.rows[q] > S.rows[q - 2])) fail("PLYSHLOG: the water velocities must increase");
                for (size_t q = 0; q < S.rows.size(); q += 2) if (!(S.rows[q] > 0)) fail("PLYSHLOG: the water velocities must be positive");
                if (S.refSalt >= 0 || S.refTemp >= 0) warn("PLYSHLOG: reference salinity / temperature are not used (no BRINE / TEMP)");
                c.polymer.plyshlog.push_back(S);
            }
            if (c.polymer.plyshlog.empty()) fail("PLYSHLOG: no table");
            while (int(c.polymer.plyshlog.size()) < ntpvt) c.polymer.plyshlog.push_back(c.polymer.plyshlog.back());
        } else if (kw == "SWOF" || kw == "SGOF" || kw == "SLGOF") {
            auto& dst = kw == "SWOF" ? swof : kw == "SGOF" ? sgof : slgof;
            dst = tables(k, 4, ntsfun);
        } else if (kw == "SWFN" || kw == "SGFN" || kw == "SOF3") {
            auto& dst = kw == "SWFN" ? swfn : kw == "SGFN" ? sgfn : sof3;
            dst = tables(k, 3, ntsfun);
        } else if (kw == "SOF2") {
            sof2 = tables(k, 2, ntsfun);
        }
    }

    // ------------------------------------------------------------------ SOLUTION
    void solution(const Keyword& k) {
        const std::string& kw = k.name;
        if (kw == "EQUIL") {
            equil.clear(); equilItem7.clear(); equilItem8.clear();
            for (const Record& r : k.records) {
                if (int(equil.size()) >= ntequl) break;
                if (r.items.empty()) { if (equil.empty()) fail("the first record is empty"); equil.push_back(equil.back()); equilItem7.push_back(equilItem7.back()); equilItem8.push_back(equilItem8.back()); continue; }
                EquilInput e;
                e.datumDepth = num(r, 0, NaN); e.datumPressure = num(r, 1, NaN);
                if (std::isnan(e.datumDepth) || std::isnan(e.datumPressure)) fail("datum depth and pressure are required", r.line);
                e.owc = num(r, 2, e.datumDepth); e.pcowOWC = num(r, 3, 0);
                e.goc = num(r, 4, e.datumDepth); e.pcogGOC = num(r, 5, 0);
                e.nacc = inum(r, 8, 0);
                equil.push_back(e); equilItem7.push_back(inum(r, 6, -99)); equilItem8.push_back(inum(r, 7, -99));   // -99: defaulted
            }
            sawEquil = !equil.empty();
        } else if (kw == "THPRES") {
            c.thpres.clear();
            for (const Record& r : k.records) {
                if (r.items.empty()) continue;
                ThresholdPressure t; t.r1 = inum(r, 0, 0); t.r2 = inum(r, 1, 0); t.value = num(r, 2, -1);
                if (t.r1 < 1 || t.r2 < 1 || t.r1 == t.r2) { warn("THPRES record with invalid regions ignored"); continue; }
                c.thpres.push_back(t);
            }
        } else if (kw == "RSVD" || kw == "PBVD") {
            auto tabs = tables(k, 2, ntequl);
            (kw == "RSVD" ? rsvd : pbvd) = tabs;
        } else if (kw == "RVVD") {
            rvvd = tables(k, 2, ntequl);
        } else if (kw == "AQUFETP" || kw == "AQUCT") {
            for (const Record& r : k.records) {
                if (r.items.empty()) continue;
                AquiferInput a;
                a.id = inum(r, 0, 1); a.depth = num(r, 1, 0); a.p0 = num(r, 2, -1);
                if (kw == "AQUFETP") {
                    a.type = AquiferType::Fetkovich;
                    a.volume = num(r, 3, NaN); a.ct = num(r, 4, NaN); a.pi = num(r, 5, NaN);
                    if (std::isnan(a.volume) || std::isnan(a.ct) || std::isnan(a.pi)) fail("volume, compressibility and PI are required", r.line);
                } else {
                    a.type = AquiferType::CarterTracy;
                    a.perm = num(r, 3, NaN); a.poro = num(r, 4, 1.0); a.ct = num(r, 5, NaN); a.rInner = num(r, 6, NaN);
                    a.thickness = num(r, 7, NaN); a.angle = num(r, 8, 360);
                    if (std::isnan(a.perm) || std::isnan(a.ct) || std::isnan(a.rInner) || std::isnan(a.thickness))
                        fail("k, ct, r0 and h are required", r.line);
                    if (r.has(10) && inum(r, 10, 1) != 1) warn("only the built-in infinite-acting influence table (AQUTAB 1) is supported", r.line);
                }
                bool dup = false;
                for (auto& x : c.aquifers) if (x.id == a.id) { x = a; dup = true; }
                if (dup) warn("aquifer " + std::to_string(a.id) + " redefined", r.line); else c.aquifers.push_back(a);
            }
        } else if (kw == "AQUANCON") {
            for (const Record& r : k.records) {
                if (r.items.empty()) continue;
                AquiferConnection a;
                const int id = inum(r, 0, 1);
                a.i1 = inum(r, 1, 1) - 1; a.i2 = inum(r, 2, 1) - 1; a.j1 = inum(r, 3, 1) - 1; a.j2 = inum(r, 4, 1) - 1;
                a.k1 = inum(r, 5, 1) - 1; a.k2 = inum(r, 6, 1) - 1;
                try { a.face = faceOf(str(r, 7, "I+"), a.minus); } catch (const std::exception& e) { fail(e.what(), r.line); }
                if (r.has(8)) grouped(1, "explicit influx coefficients ignored (face areas are used)");
                a.mult = num(r, 9, 1.0);
                if (upper(str(r, 10, "NO")) == "YES") grouped(1, "connections to internal faces are not supported");
                aqConn[id].push_back(a);
            }
        }
    }

    // ------------------------------------------------------------------ SCHEDULE
    WellInput* findWell(const std::string& name) {
        for (auto& w : c.wells) if (w.name == name) return &w;
        return nullptr;
    }
    void addControl(const WellControl& wcIn, int inj = -1) {   // inj: 1 injector, 0 producer, -1 same as before
        WellControl wc = wcIn;
        // ECLIPSE / OPM Flow: a well without any connection at this time ignores an OPEN request and remains shut
        // (connections that exist but are SHUT do not count - they may be opened by WELOPEN in the same step)
        if (wc.status == WellStatus::Open) {
            if (const WellInput* w = findWell(wc.well)) {
                const int cur = anyStep ? int(c.schedule.size()) : 0;
                bool any = false;
                for (auto& q : w->completions)
                    if ((q.open && q.fromStep <= cur && (q.untilStep < 0 || q.untilStep > cur)) || (!q.open && q.fromStep <= cur)) any = true;
                if (!any) { wc.status = WellStatus::Shut; grouped(1, "well has no connections when its control is set - it remains SHUT (as in ECLIPSE)", wc.well); }
            }
        }
        lastCtrl[wc.well] = wc;
        const bool isInj = inj < 0 ? lastInj[wc.well] : inj == 1;
        lastInj[wc.well] = isInj; pendingInj[wc.well] = isInj;
        for (auto& p : pending) if (p.well == wc.well) { p = wc; return; }
        pending.push_back(wc);
    }
    bool ctrlSeen = false;            // a native well-control keyword was read in the current report step
    void addStep(double len) {
        ctrlSeen = false;
        if (!(len > 0)) fail("time step lengths must be positive");
        ReportStep s; s.length = len; s.controls = pending; pending.clear();
        s.wpolymer = pendingPoly; pendingPoly.clear();
        s.drsdt = pendingDrsdt; s.drsdtFree = pendingDrsdtFree; pendingDrsdt = -2; pendingDrsdtFree = false;
        for (const auto& wc : s.controls) if (pendingInj[wc.well]) injCtrl.insert({c.schedule.size(), wc.well});
        pendingInj.clear();
        c.schedule.push_back(s);
        anyStep = true;
    }
    static WellStatus statusOf(const std::string& s) {
        const std::string u = upper(s);
        if (u == "SHUT") return WellStatus::Shut;
        if (u == "STOP") return WellStatus::Stop;
        return WellStatus::Open;              // OPEN / AUTO
    }
    static char faceOf(const std::string& s, bool& minus) {
        std::string u = upper(s);
        minus = u.size() > 1 && u[1] == '-';
        char f = u.empty() ? 'X' : u[0];
        if (f == 'I') f = 'X'; else if (f == 'J') f = 'Y'; else if (f == 'K') f = 'Z';
        if (f != 'X' && f != 'Y' && f != 'Z') throw std::runtime_error("unknown face '" + s + "'");
        return f;
    }
    std::vector<WellInput*> wellsMatching(const std::string& pat) {
        std::vector<WellInput*> ws;
        for (auto& w : c.wells) if (matchName(pat, w.name)) ws.push_back(&w);
        return ws;
    }

    void schedule(const Keyword& k) {
        if (k.name == "WCONPROD" || k.name == "WCONHIST" || k.name == "WCONINJE" || k.name == "WCONINJH" || k.name == "WELTARG" || k.name == "WELOPEN")
            ctrlSeen = true;
        const std::string& kw = k.name;
        if (kw == "WELSPECS") {
            for (const Record& r : k.records) {
                if (r.items.empty()) continue;
                const std::string name = str(r, 0);
                if (name.empty()) fail("well name missing", r.line);
                WellInput* w = findWell(name);
                const bool isNew = !w;
                if (!w) { c.wells.push_back(WellInput{}); w = &c.wells.back(); w->completions.clear(); w->name = name; }
                w->group = str(r, 1, "G");
                w->i = inum(r, 2, 1) - 1; w->j = inum(r, 3, 1) - 1;
                w->refDepth = num(r, 4, -1);
                const std::string ph = upper(str(r, 5, "OIL"));
                if (isNew) w->type = ph == "WATER" || ph == "WAT" ? WellType::WaterInjector : ph == "GAS" ? WellType::GasInjector : WellType::Producer;
                {   // items the ResSim engine does not use, kept for OPM Flow
                    std::string tail; size_t last = 0;
                    for (size_t i = 6; i < r.items.size(); ++i) if (!r.items[i].def) last = i + 1;
                    for (size_t i = 6; i < last; ++i) {
                        const auto& it = r.items[i];
                        tail += (tail.empty() ? "" : " ") + (it.def ? std::string("1*") : (it.quoted ? "'" + it.s + "'" : it.s));
                    }
                    w->specsTail = tail;
                }
                if (anyStep && isNew) grouped(0, "wells defined after the first report step exist from time zero (shut until opened)", name);
            }
        } else if (kw == "COMPDAT") {
            for (const Record& r : k.records) {
                if (r.items.empty()) continue;
                const std::string name = str(r, 0);
                std::vector<WellInput*> ws = wellsMatching(name);
                if (ws.empty()) { warn("well '" + name + "' is not defined in WELSPECS - record ignored", r.line); continue; }
                for (WellInput* w : ws) {
                    Completion cp;
                    cp.i = (r.has(1) && inum(r, 1, 0) > 0) ? inum(r, 1, 0) - 1 : w->i;
                    cp.j = (r.has(2) && inum(r, 2, 0) > 0) ? inum(r, 2, 0) - 1 : w->j;
                    cp.k1 = inum(r, 3, 1) - 1; cp.k2 = inum(r, 4, cp.k1 + 1) - 1;
                    cp.open = statusOf(str(r, 5, "OPEN")) == WellStatus::Open;
                    cp.cf = num(r, 7, -1); cp.diameter = num(r, 8, c.units == UnitSystem::Field ? 0.5 : 0.1524);
                    cp.kh = num(r, 9, -1); cp.skin = num(r, 10, 0);
                    const std::string d = upper(str(r, 12, "Z"));
                    cp.dir = d.empty() ? 'Z' : d[0];
                    cp.satTable = inum(r, 6, 0); cp.dfactor = num(r, 11, -1); cp.r0 = num(r, 13, -1);
                    if (cp.i < 0 || cp.i >= c.grid.nx || cp.j < 0 || cp.j >= c.grid.ny || cp.k1 < 0 || cp.k2 >= c.grid.nz || cp.k1 > cp.k2) {
                        warn("connection of " + w->name + " outside the grid - ignored", r.line); continue;
                    }
                    // connections specified in later report steps open (or shut) at that step
                    const int cur = anyStep ? int(c.schedule.size()) : 0;
                    bool replaced = false;
                    for (auto& e : w->completions) {
                        if (!(e.i == cp.i && e.j == cp.j && e.k1 == cp.k1 && e.k2 == cp.k2)) continue;
                        replaced = true;
                        if (cur == 0) { e = cp; continue; }
                        if (!cp.open) { if (e.open && e.untilStep < 0 && cur > e.fromStep) e.untilStep = cur; else if (e.fromStep >= cur) e.open = false; }
                        else if (!e.open || e.untilStep >= 0) {
                            const Completion old = e; e = cp; e.fromStep = (!old.open) ? cur : old.fromStep; e.untilStep = -1;
                            if (old.open && old.untilStep >= 0) grouped(1, "connection re-opened after being shut: treated as open throughout", w->name);
                        } else { const int f = e.fromStep; e = cp; e.fromStep = f; }     // new data for an open connection
                    }
                    if (!replaced) { cp.fromStep = cur; w->completions.push_back(cp); }   // SHUT ones too (WELOPEN may open them)
                }
            }
        } else if (kw == "WCONPROD" || kw == "WCONHIST") {
            const bool hist = kw == "WCONHIST";
            for (const Record& r : k.records) {
                if (r.items.empty()) continue;
                WellControl wc; wc.well = str(r, 0); wc.status = statusOf(str(r, 1, "OPEN")); wc.history = hist;
                const std::string m = upper(str(r, 2, hist ? "RESV" : "NONE"));
                wc.orat = num(r, 3, -1); wc.wrat = num(r, 4, -1); wc.grat = num(r, 5, -1);
                if (hist) { wc.bhp = num(r, 9, -1); wc.lrat = std::max(0.0, wc.orat) + std::max(0.0, wc.wrat);
                            wc.vfp = inum(r, 6, 0); wc.alq = num(r, 7, -1); wc.thp = num(r, 8, -1); }
                else { wc.lrat = num(r, 6, -1); wc.resv = num(r, 7, -1); wc.bhp = num(r, 8, -1);
                       wc.thp = num(r, 9, -1); wc.vfp = inum(r, 10, 0); wc.alq = num(r, 11, -1); }
                if (m == "ORAT") wc.mode = Control::ORAT; else if (m == "WRAT") wc.mode = Control::WRAT;
                else if (m == "GRAT") wc.mode = Control::GRAT; else if (m == "LRAT") wc.mode = Control::LRAT;
                else if (m == "BHP") wc.mode = Control::BHP;
                else if (m == "RESV") {
                    wc.mode = Control::RESV;
                    if (hist) wc.resv = -1;   // computed by the simulator from the observed rates
                } else {
                    // GRUP / THP / NONE / CRAT: pick the most restrictive explicit target we have (the deck mode is kept for export)
                    wc.deckMode = m;
                    wc.mode = wc.orat > 0 ? Control::ORAT : wc.lrat > 0 ? Control::LRAT : wc.wrat > 0 ? Control::WRAT :
                              wc.grat > 0 ? Control::GRAT : Control::BHP;
                    grouped(1, "control mode " + m + " is not modelled by the ResSim engine (kept for OPM Flow) - ResSim uses " +
                               std::string(wc.mode == Control::BHP ? "BHP" : wc.mode == Control::ORAT ? "ORAT" : wc.mode == Control::LRAT ? "LRAT" :
                                           wc.mode == Control::WRAT ? "WRAT" : "GRAT"), wc.well);
                }
                if (wc.mode == Control::BHP && wc.bhp <= 0 && wc.status == WellStatus::Open) {
                    wc.bhp = c.units == UnitSystem::Field ? 14.7 : 1.01325;
                    grouped(1, "no BHP limit given - 1 atm used", wc.well);
                }
                const std::vector<WellInput*> ws = wellsMatching(wc.well);
                if (ws.empty()) { warn("well '" + wc.well + "' is not defined - record ignored", r.line); continue; }
                for (WellInput* w : ws) {
                    w->type = WellType::Producer;
                    WellControl x = wc; x.well = w->name; addControl(x, 0);
                }
            }
        } else if (kw == "WCONINJE" || kw == "WCONINJH") {
            const bool hist = kw == "WCONINJH";
            for (const Record& r : k.records) {
                if (r.items.empty()) continue;
                WellControl wc; wc.well = str(r, 0);
                std::string type = upper(str(r, 1, "WATER"));
                if (type == "WAT") type = "WATER";
                if (type != "WATER" && type != "GAS") { keepRecord(k, r, "injected fluid"); grouped(1, "injection of " + type + " is not modelled by the ResSim engine - record kept for OPM Flow", wc.well); continue; }
                wc.status = statusOf(str(r, 2, "OPEN"));
                if (hist) {
                    wc.mode = upper(str(r, 11, "RATE")) == "BHP" ? Control::BHP : Control::RATE;
                    wc.rate = num(r, 3, 0); wc.bhp = num(r, 4, -1);
                } else {
                    const std::string m = upper(str(r, 3, "RATE"));
                    wc.rate = num(r, 4, -1); wc.bhp = num(r, 6, -1);
                    wc.resv = num(r, 5, -1); wc.thp = num(r, 7, -1); wc.vfp = inum(r, 8, 0);
                    if (m == "BHP") wc.mode = Control::BHP;
                    else if (m == "RATE") wc.mode = Control::RATE;
                    else if (m == "RESV") { grouped(1, "RESV injection: ResSim uses it as a surface RATE (kept as RESV for OPM Flow)", wc.well); wc.deckMode = m;
                                            wc.mode = Control::RATE; if (wc.rate <= 0) { wc.rate = wc.resv > 0 ? wc.resv : 0; wc.deckMode = "RESV*"; } }
                    else { grouped(1, "injection control " + m + " is not modelled by the ResSim engine (kept for OPM Flow) - ResSim uses " + (wc.rate > 0 ? std::string("RATE") : std::string("BHP")), wc.well);
                           wc.deckMode = m; wc.mode = wc.rate > 0 ? Control::RATE : Control::BHP; }
                }
                if (wc.bhp <= 0) wc.bhp = c.units == UnitSystem::Field ? 1e5 : 6895;   // ECLIPSE default BHP limit (100000 psia)
                const std::vector<WellInput*> ws = wellsMatching(wc.well);
                if (ws.empty()) { warn("well '" + wc.well + "' is not defined - record ignored", r.line); continue; }
                for (WellInput* w : ws) {
                    if (w->type == WellType::Producer) w->type = type == "GAS" ? WellType::GasInjector : WellType::WaterInjector;
                    WellControl x = wc; x.well = w->name; x.phase = type == "GAS" ? 2 : 1; addControl(x, 1);
                }
            }
        } else if (kw == "WELTARG") {
            for (const Record& r : k.records) {
                if (r.items.empty()) continue;
                const std::string pat = str(r, 0), what = upper(str(r, 1));
                const double v = num(r, 2, NaN);
                if (std::isnan(v)) continue;
                static const std::set<std::string> native = {"ORAT", "WRAT", "GRAT", "LRAT", "RESV", "BHP", "RATE", "THP", "VFP", "LIFT"};
                if (!native.count(what)) { keepRecord(k, r, "WELTARG target"); grouped(1, "WELTARG " + what + " is not modelled by the ResSim engine - kept for OPM Flow", pat); continue; }
                for (WellInput* w : wellsMatching(pat)) {
                    auto it = lastCtrl.find(w->name);
                    if (it == lastCtrl.end()) { grouped(1, "no previous control for the well - ignored", w->name); continue; }
                    WellControl wc = it->second;
                    if (what == "ORAT") wc.orat = v; else if (what == "WRAT") wc.wrat = v; else if (what == "GRAT") wc.grat = v;
                    else if (what == "LRAT") wc.lrat = v; else if (what == "RESV") wc.resv = v; else if (what == "BHP") wc.bhp = v;
                    else if (what == "RATE") wc.rate = v;
                    else if (what == "THP") wc.thp = v; else if (what == "VFP") wc.vfp = int(v + 0.5); else if (what == "LIFT") wc.alq = v;
                    addControl(wc);
                }
            }
        } else if (kw == "DRSDT") {
            if (!k.records.empty() && !k.records[0].items.empty()) {
                const Record& r = k.records[0];
                const double v = num(r, 0, 0);
                if (v < 0) fail("DRSDT must be >= 0");
                pendingDrsdt = v >= 1e19 ? -1 : v;
                std::string opt = r.has(1) ? str(r, 1) : "ALL";
                for (auto& ch : opt) ch = char(std::toupper((unsigned char)ch));
                pendingDrsdtFree = opt == "FREE";
            }
        } else if (kw == "WPOLYMER") {
            for (const Record& r : k.records) {
                if (r.items.empty()) continue;
                WellPolymer wp; wp.well = str(r, 0); wp.conc = num(r, 1, 0);
                if (r.has(3) || r.has(4)) grouped(1, "WPOLYMER group concentrations are not modelled - the well concentration is used", wp.well);
                if (num(r, 2, 0) != 0) grouped(1, "WPOLYMER salt concentration needs BRINE - ignored", wp.well);
                pendingPoly.push_back(wp);
            }
        } else if (kw == "WPIMULT") {
            for (const Record& r : k.records) {
                if (r.items.empty()) continue;
                const double m = num(r, 1, 1.0);
                if (r.has(2) || r.has(3) || r.has(4)) grouped(1, "multipliers on individual connections applied to the whole well");
                for (WellInput* w : wellsMatching(str(r, 0))) w->piMult *= m;
            }
        } else if (kw == "WELOPEN") {
            for (const Record& r : k.records) {
                if (r.items.empty()) continue;
                const WellStatus st = statusOf(str(r, 1, "OPEN"));
                if (r.has(2) || r.has(3) || r.has(4) || r.has(5) || r.has(6)) {
                    // connection-level: I J K (0 / default = any); completion numbers (COMPLUMP) are kept verbatim
                    const int I = inum(r, 2, 0), J = inum(r, 3, 0), K = inum(r, 4, 0), C1 = inum(r, 5, 0), C2 = inum(r, 6, 0);
                    const int cur = anyStep ? int(c.schedule.size()) : 0;
                    bool representable = C1 <= 0 && C2 <= 0;
                    std::vector<Completion*> hit;
                    for (WellInput* w : wellsMatching(str(r, 0)))
                        for (auto& q : w->completions) {
                            if (q.fromStep > cur) continue;
                            if ((I > 0 && q.i + 1 != I) || (J > 0 && q.j + 1 != J) || (K > 0 && (K < q.k1 + 1 || K > q.k2 + 1))) continue;
                            if (K > 0 && q.k1 != q.k2) representable = false;      // one layer of a multi-layer connection
                            hit.push_back(&q);
                        }
                    if (!representable) { keepRecord(k, r, "connection-level WELOPEN"); grouped(1, "WELOPEN on completion numbers / part of a connection range: not modelled by the ResSim engine - kept for OPM Flow", str(r, 0)); continue; }
                    for (Completion* q : hit) {
                        if (st == WellStatus::Open) {
                            if (!q->open) { q->open = true; q->fromStep = cur; q->untilStep = -1; }
                            else if (q->untilStep >= 0 && q->untilStep <= cur) { q->untilStep = -1; grouped(1, "connection re-opened after being shut: treated as open throughout", str(r, 0)); }
                        } else if (q->open && q->untilStep < 0) {
                            if (q->fromStep >= cur) q->open = false; else q->untilStep = cur;
                        }
                    }
                    continue;
                }
                for (WellInput* w : wellsMatching(str(r, 0))) {
                    auto it = lastCtrl.find(w->name);
                    WellControl wc;
                    if (it != lastCtrl.end()) wc = it->second;
                    else if (st == WellStatus::Open) { grouped(1, "no previous control - ignored", w->name); continue; }
                    else { wc.well = w->name; wc.mode = Control::BHP; }
                    wc.status = st;
                    addControl(wc);
                }
            }
        } else if (kw == "TSTEP") {
            for (const Record& r : k.records)
                for (double v : numbers(r)) { if (std::isnan(v)) fail("defaulted step length"); addStep(v); daysElapsed += v; }
        } else if (kw == "DATES") {
            for (const Record& r : k.records) {
                if (r.items.empty()) continue;
                long d = 0;
                try { d = daysFromCivil(inum(r, 2, 0), monthFrom(str(r, 1, "JAN")), inum(r, 0, 1)); }
                catch (const std::exception& e) { fail(e.what(), r.line); }
                double t = 0;
                if (r.has(3)) {   // HH:MM:SS
                    int hh = 0, mm = 0, ss = 0; std::sscanf(r.items[3].s.c_str(), "%d:%d:%d", &hh, &mm, &ss);
                    t = (hh + mm / 60.0 + ss / 3600.0) / 24.0;
                }
                const double target = double(d - dayNow) + t;
                const double len = target - daysElapsed;
                if (std::fabs(len) < 1e-9) { info("repeated date " + str(r, 0) + " " + str(r, 1) + " " + str(r, 2) + " ignored", r.line); continue; }
                if (len < 0) fail("dates must increase (" + str(r, 0) + " " + str(r, 1) + " " + str(r, 2) + ")", r.line);
                addStep(len); daysElapsed = target;
            }
        } else if (kw == "TUNING") {
            static const Record empty;
            const Record& r1 = k.records.size() > 0 ? k.records[0] : empty;
            const Record& r2 = k.records.size() > 1 ? k.records[1] : empty;
            const Record& r3 = k.records.size() > 2 ? k.records[2] : empty;
            Tuning& t = c.tuning;
            // Only the time-step controls are taken over. Like OPM Flow, the ECLIPSE linear-solver limits (LITMAX),
            // the minimum step (TSMINZ, 0.1 day by default) and DDSLIM (1e6) are not meaningful for this solver.
            t.initialDt = num(r1, 0, t.initialDt); t.maxDt = num(r1, 1, t.maxDt);
            t.growth = std::max(1.2, num(r1, 4, t.growth)); t.cutFactor = std::min(0.5, std::max(0.1, num(r1, 6, t.cutFactor)));
            t.cnvTol = std::max(1e-5, num(r2, 1, t.cnvTol)); t.mbTol = std::max(1e-9, num(r2, 2, t.mbTol));
            t.maxNewton = std::max(t.maxNewton, inum(r3, 0, t.maxNewton));
            if (num(r3, 7, 1e6) < 1.0) t.dsMax = num(r3, 7, t.dsMax);
            if (anyStep) grouped(0, "TUNING after the first report step applies to the whole run");
        }
    }

    // ------------------------------------------------------------------ RESSIM extension lines
    void ressim(const Record& r) {
        const std::string key = upper(str(r, 0));
        if (key == "TUNING") { c.tuning.dpMaxRel = num(r, 1, c.tuning.dpMaxRel); c.tuning.threads = inum(r, 2, 0); c.tuning.linearSolver = inum(r, 3, 2); c.tuning.cprWeights = inum(r, 4, 1) ? 1 : 0; }
        else if (key == "AQVISC") aqVisc[inum(r, 1, 0)] = num(r, 2, -1);
        else if (key == "PARAM") {
            StudyParam p; p.kind = ParamKind(std::clamp(inum(r, 1, 0), 0, int(ParamKind::ProdRate))); p.target = inum(r, 2, 0);
            p.i1 = inum(r, 3, -1); p.i2 = inum(r, 4, -1); p.j1 = inum(r, 5, -1); p.j2 = inum(r, 6, -1);
            p.k1 = inum(r, 7, -1); p.k2 = inum(r, 8, -1);
            p.base = num(r, 9, 1); p.low = num(r, 10, 0.5); p.high = num(r, 11, 2);
            p.log = inum(r, 12, 1) != 0; p.active = inum(r, 13, 1) != 0; p.label = str(r, 14);
            c.studyParams.push_back(p);
        } else if (key == "OBS") {
            Observation ob; const std::string q = upper(str(r, 1)); bool f = false;
            for (int i = 0; i <= int(ObsQty::FGOR); ++i) if (q == toString(ObsQty(i))) { ob.qty = ObsQty(i); f = true; }
            if (!f) { warn("unknown observation quantity " + q); return; }
            ob.well = str(r, 2, "-"); if (ob.well == "-") ob.well.clear();
            ob.time = num(r, 3, 0); ob.value = num(r, 4, 0); ob.sigma = num(r, 5, -1); ob.weight = num(r, 6, 1);
            c.observations.push_back(ob);
        } else if (key == "INITRS") {
            // RESSIM: INITRS <region> <Rs>  (constant initial Rs - not expressible in ECLIPSE without RSVD)
            const int reg = inum(r, 1, 1) - 1; const double rs = num(r, 2, -1);
            if (reg >= 0 && reg < int(equil.size())) equil[size_t(reg)].initialRs = rs;
        }
    }

    // ------------------------------------------------------------------ main loop
    void run() {
        static const std::set<std::string> sections = {"RUNSPEC", "GRID", "EDIT", "PROPS", "REGIONS", "SOLUTION", "SUMMARY", "SCHEDULE"};
        static const std::set<std::string> gridOps = {"EQUALS", "MULTIPLY", "ADD", "COPY", "MINVALUE", "MAXVALUE"};
        static const std::set<std::string> regOps = {"EQUALREG", "MULTIREG", "ADDREG", "COPYREG"};
        static const std::set<std::string> propKws = {"PVTW", "ROCK", "DENSITY", "GRAVITY", "PVCDO", "PVDO", "PVDG", "PVTO", "PVTG",
                                                      "PLYVISC", "PLYADS", "PLYROCK", "PLYMAX", "PLMIXPAR", "PLYSHLOG", "STONE", "STONE1", "STONE2", "STONE1EX",
                                                      "SWOF", "SGOF", "SLGOF", "SWFN", "SGFN", "SOF3", "SOF2"};
        static const std::set<std::string> solKws = {"EQUIL", "RSVD", "PBVD", "AQUFETP", "AQUCT", "AQUANCON", "THPRES", "RVVD"};
        static const std::set<std::string> schKws = {"WELSPECS", "COMPDAT", "WCONPROD", "WCONHIST", "WCONINJE", "WCONINJH", "WELTARG",
                                                     "WPIMULT", "WELOPEN", "TSTEP", "DATES", "TUNING", "WPOLYMER", "DRSDT"};
        std::map<std::string, std::pair<int, const Keyword*>> unsupported;   // name -> (count, first)
        // keywords written by the exporter itself (not kept as extra keywords, to avoid duplicates)
        static const std::set<std::string> exporterOwned = {"UNIFOUT", "INIT", "RPTRST", "WELLDIMS", "AQUDIMS", "FAULTDIM", "EQLOPTS",
                                                            "NOECHO", "ECHO", "TABDIMS", "EQLDIMS", "REGDIMS", "START", "TITLE"};
        auto keep = [&](const Keyword& k, const char* note) {
            static const std::map<std::string, DeckSection> secOf = {
                {"RUNSPEC", DeckSection::Runspec}, {"GRID", DeckSection::Grid}, {"EDIT", DeckSection::Edit}, {"PROPS", DeckSection::Props},
                {"REGIONS", DeckSection::Regions}, {"SOLUTION", DeckSection::Solution}, {"SUMMARY", DeckSection::Summary},
                {"SCHEDULE", DeckSection::Schedule}};
            if (exporterOwned.count(k.name)) {
                // dimension keywords: the source values are kept as MINIMA (groups in GRUPTREE, VFP / MSW sizes, ... are not
                // all visible to the exporter); it writes max(source, what the model needs) item by item
                if (!kDimsKeywords.count(k.name)) return;
                ExtraKeyword e; e.name = k.name; e.note = kDimsNote; e.section = DeckSection::Runspec;
                e.body = deck::keywordBodyText(k);
                c.extraKeywords.push_back(std::move(e));
                return;
            }
            ExtraKeyword e; e.name = k.name; e.note = note;
            auto it = secOf.find(section); e.section = it == secOf.end() ? DeckSection::Runspec : it->second;
            // a keyword in a section where it is not valid (e.g. summary vectors without a SUMMARY header) is written to
            // its own section - OPM Flow / ECLIPSE reject it otherwise
            if (const deck::KeywordSpec* sp = deck::findKeyword(k.name)) {
                const unsigned m = sp->sections;
                if (m && !(m & (1u << int(e.section)))) {
                    int to = (m & (1u << int(DeckSection::Summary))) ? int(DeckSection::Summary) : -1;
                    for (int b = 0; b < 8 && to < 0; ++b) if (m & (1u << b)) to = b;
                    static const char* nm[] = {"RUNSPEC", "GRID", "EDIT", "PROPS", "REGIONS", "SOLUTION", "SUMMARY", "SCHEDULE"};
                    grouped(1, std::string("not valid in ") + nm[int(e.section)] + " - written to the " + nm[to] + " section", k.name);
                    e.section = DeckSection(to);
                }
            }
            e.step = e.section == DeckSection::Schedule ? int(c.schedule.size()) : 0;
            e.beforeControls = !ctrlSeen;
            e.body = deck::keywordBodyText(k);
            c.extraKeywords.push_back(std::move(e));
        };
        bool ended = false;
        for (const Keyword& k : D.keywords) {
            K = &k;
            const std::string& kw = k.name;
            try {
                if (ended) break;
                if (sections.count(kw)) { section = kw; continue; }
                if (kw == "END") { ended = true; continue; }
                if (kw == "__RESSIM") { for (const Record& r : k.records) ressim(r); continue; }
                if (section == "SUMMARY") { if (!exporterSummary().count(kw)) keep(k, "summary request"); continue; }   // output requests (OPM Flow writes them)
                if (!k.known) { keep(k, "unknown keyword"); unsupported[kw].first++; if (!unsupported[kw].second) unsupported[kw].second = &k; continue; }
                // ---------------- RUNSPEC
                if (kw == "OIL") { c.phases.oil = true; continue; }
                if (kw == "WATER") { c.phases.water = true; continue; }
                if (kw == "GAS") { c.phases.gas = true; continue; }
                if (kw == "DISGAS") { c.phases.disgas = true; continue; }
                if (kw == "VAPOIL") { c.phases.vapoil = true; continue; }
                if (kw == "POLYMER") { c.phases.polymer = true; continue; }
                if (kw == "ENDSCALE") {
                    c.grid.endscale = true;
                    if (!k.records.empty() && !k.records[0].items.empty()) {
                        const Record& r = k.records[0];
                        auto opt = [&](size_t i) { return i < r.items.size() && !r.items[i].def ? upper(r.items[i].s) : std::string(); };
                        if (opt(0) == "DIRECT") warn("ENDSCALE: directional end-point scaling is not modelled by the ResSim engine - isotropic scaling is used");
                        if (opt(1) == "IRREVERS") warn("ENDSCALE: irreversible scaling is not modelled by the ResSim engine");
                    }
                    continue;
                }
                if (kw == "FIELD") { c.units = UnitSystem::Field; continue; }
                if (kw == "METRIC") { c.units = UnitSystem::Metric; continue; }
                if (kw == "LAB" || kw == "PVT-M") fail("unit system " + kw + " is not supported (use FIELD or METRIC)");
                if (kw == "TITLE") { std::string t = k.raw; while (!t.empty() && std::isspace((unsigned char)t.back())) t.pop_back();
                                     size_t b = 0; while (b < t.size() && std::isspace((unsigned char)t[b])) ++b; c.title = t.substr(b); continue; }
                if (kw == "DIMENS" || kw == "SPECGRID") {
                    if (k.records.empty()) fail("record missing");
                    const Record& r = k.records.front();
                    const int nx = inum(r, 0, 0), ny = inum(r, 1, 0), nz = inum(r, 2, 0);
                    if (n != 0 && (nx != c.grid.nx || ny != c.grid.ny || nz != c.grid.nz)) fail("dimensions differ from DIMENS");
                    if (nx < 1 || ny < 1 || nz < 1) fail("dimensions must be positive");
                    if (kw == "SPECGRID" && inum(r, 3, 1) > 1) warn("multiple reservoirs (NUMRES > 1) are treated as one");
                    c.grid.nx = nx; c.grid.ny = ny; c.grid.nz = nz;
                    n = nx * ny * nz; nl = nx * ny;
                    continue;
                }
                if (kw == "TABDIMS") { if (!k.records.empty()) { ntsfun = std::max(1, inum(k.records[0], 0, 1)); ntpvt = std::max(1, inum(k.records[0], 1, 1)); } keep(k, ""); continue; }
                if (kw == "EQLDIMS") { if (!k.records.empty()) ntequl = std::max(1, inum(k.records[0], 0, 1)); continue; }
                if (kw == "START") {
                    if (k.records.empty()) continue;
                    const Record& r = k.records.front();
                    c.startDay = inum(r, 0, 1);
                    try { c.startMonth = monthFrom(str(r, 1, "JAN")); } catch (const std::exception& e) { fail(e.what()); }
                    c.startYear = inum(r, 2, 1983);
                    dayNow = daysFromCivil(c.startYear, c.startMonth, c.startDay);
                    continue;
                }
                // ---------------- unsupported input that changes the model fundamentally
                if (kw == "GDFILE" || kw == "IMPORT") fail("binary grid files are not supported - convert the grid to COORD/ZCORN text");
                if (kw == "RESTART") { keep(k, "restart"); warn("restart runs are not supported by the ResSim engine (it starts at time zero) - kept for OPM Flow"); continue; }
                if (kw == "LGR" || kw == "CARFIN" || kw == "RADFIN" || kw == "RADFIN4") { keep(k, "local grid refinement"); grouped(1, "local grid refinements are not modelled by the ResSim engine - kept for OPM Flow"); continue; }
                if (kw == "DUALPORO") { c.grid.dual = std::max(c.grid.dual, 1); continue; }
                if (kw == "DUALPERM") { c.grid.dual = 2; continue; }
                if (kw == "NODPPM") { c.grid.nodppm = true; continue; }
                if (kw == "DPGRID") continue;   // fracture cells always take the matrix geometry
                if (kw == "SIGMA") { if (!k.records.empty()) c.grid.sigma = num(k.records[0], 0, 0.0); continue; }
                if (kw == "LTOSIGMA") {
                    c.grid.ltosigma = true;
                    if (!k.records.empty()) for (int d = 0; d < 3; ++d) c.grid.ltoF[d] = num(k.records[0], d, 4.0);
                    continue;
                }
                if (kw == "GRAVDR" || kw == "GRAVDRM" || kw == "GRAVDRB" || kw == "DZMTRX") {
                    if (kw != "DZMTRX") c.grid.gravdr = true;
                    keep(k, "dual porosity gravity drainage");
                    grouped(1, "dual porosity gravity drainage (GRAVDR/DZMTRX) is not modelled by the ResSim engine - kept for OPM Flow / ECLIPSE");
                    continue;
                }
                if (kw == "RADIAL") fail("radial grids are not supported");
                // ---------------- GRID / EDIT / REGIONS / SOLUTION arrays
                if (k.isArray) {
                    dataArray(k);
                    continue;
                }
                if (kw == "BOX") { needDims(); boxActive = false; if (!k.records.empty()) { box = readBox(k.records[0], 0, fullBox()); boxActive = true; } continue; }
                if (kw == "ENDBOX") { boxActive = false; continue; }
                if (gridOps.count(kw)) { operations(k); continue; }
                if (regOps.count(kw)) { regionOperations(k); continue; }
                if (kw == "OPERATE") { keep(k, "array arithmetic"); warn("OPERATE (array arithmetic) is not applied by the ResSim engine - kept for OPM Flow"); continue; }
                if (kw == "NEWTRAN") { c.grid.newtran = true; continue; }
                if (kw == "OLDTRAN") { c.grid.newtran = false; continue; }
                if (kw == "MINPV" || kw == "MINPORV") { if (!k.records.empty()) c.grid.minpv = num(k.records[0], 0, 1e-6); continue; }
                if (kw == "PINCH") {
                    auto& P = c.grid.pinch; P.on = true;
                    if (!k.records.empty()) {
                        const Record& r = k.records[0];
                        P.threshold = num(r, 0, 0.001);
                        P.gap = upper(str(r, 1, "GAP")) != "NOGAP";
                        P.maxGap = num(r, 2, 1e20);
                        P.topbot = upper(str(r, 3, "TOPBOT")) != "ALL";
                        P.multzAll = upper(str(r, 4, "TOP")) == "ALL";
                    }
                    if (!P.topbot) warn("PINCH item 4 ALL (transmissibility through the gap cells) - the ResSim engine uses the TOPBOT form");
                    continue;
                }
                if (kw == "FAULTS") {
                    needDims();
                    for (const Record& r : k.records) {
                        if (r.items.empty()) continue;
                        const std::string name = str(r, 0, "F");
                        FaultFace ff;
                        ff.i1 = inum(r, 1, 1) - 1; ff.i2 = inum(r, 2, 1) - 1; ff.j1 = inum(r, 3, 1) - 1; ff.j2 = inum(r, 4, 1) - 1;
                        ff.k1 = inum(r, 5, 1) - 1; ff.k2 = inum(r, 6, 1) - 1;
                        try { ff.face = faceOf(str(r, 7, "X"), ff.minus); } catch (const std::exception& e) { fail(e.what(), r.line); }
                        FaultInput* f = nullptr;
                        for (auto& x : c.grid.faults) if (x.name == name) f = &x;
                        if (!f) { c.grid.faults.push_back(FaultInput{name, 1.0, {}}); f = &c.grid.faults.back(); }
                        f->faces.push_back(ff);
                    }
                    continue;
                }
                if (kw == "MULTFLT") {
                    if (section == "SCHEDULE") { grouped(1, "fault multipliers in SCHEDULE are applied from time zero"); }
                    for (const Record& r : k.records) {
                        if (r.items.empty()) continue;
                        const std::string name = str(r, 0); const double m = num(r, 1, 1.0);
                        bool any = false;
                        for (auto& f : c.grid.faults) if (matchName(name, f.name)) { f.mult *= m; any = true; }
                        if (!any) grouped(1, "fault not defined in FAULTS - ignored", name);
                    }
                    continue;
                }
                if (kw == "MULTREGT") {
                    if (section == "SCHEDULE") { grouped(1, "MULTREGT in SCHEDULE is not supported - ignored"); continue; }
                    for (const Record& r : k.records) {
                        if (r.items.empty()) continue;
                        MultRegT m; m.r1 = inum(r, 0, -1); m.r2 = inum(r, 1, -1); m.m = num(r, 2, 1.0);
                        m.dirs = upper(str(r, 3, "XYZ")); const std::string set = upper(str(r, 5, "M"));
                        m.set = set.empty() || set[0] == 'M' ? "MULTNUM" : set[0] == 'F' ? "FLUXNUM" : "OPERNUM";
                        multregt.push_back(m);
                    }
                    continue;
                }
                if (kw == "NNC") {
                    for (const Record& r : k.records) {
                        if (r.items.empty()) continue;
                        NNCInput e;
                        e.i1 = inum(r, 0, 1) - 1; e.j1 = inum(r, 1, 1) - 1; e.k1 = inum(r, 2, 1) - 1;
                        e.i2 = inum(r, 3, 1) - 1; e.j2 = inum(r, 4, 1) - 1; e.k2 = inum(r, 5, 1) - 1; e.trans = num(r, 6, 0);
                        c.grid.nnc.push_back(e);
                    }
                    continue;
                }
                if (propKws.count(kw)) { props(k); continue; }
                if (solKws.count(kw)) { solution(k); continue; }
                if (schKws.count(kw)) { schedule(k); continue; }
                if (kw == "ACTIONX" || kw == "PYACTION") { keep(k, "action block"); grouped(1, "action blocks are not evaluated by the ResSim engine - kept for OPM Flow"); continue; }
                if (kw == "WELSEGS" || kw == "COMPSEGS") { keep(k, "multisegment well"); grouped(1, "multisegment well data: standard well in the ResSim engine - kept for OPM Flow"); continue; }
                // ---------------- everything else: known keyword without a mapping
                auto st = supportTable().find(kw);
                if (st != supportTable().end()) {
                    if (st->second.s == Support::Full) continue;
                    if (st->second.s == Support::Ignored) { keep(k, "output / option keyword"); continue; }
                    if (st->second.s == Support::Partial) { keep(k, st->second.note); grouped(1, st->second.note + std::string(" - kept for OPM Flow")); continue; }
                }
                keep(k, "not modelled by the ResSim engine");
                unsupported[kw].first++;
                if (!unsupported[kw].second) unsupported[kw].second = &k;
            } catch (const std::runtime_error& e) {
                const std::string w = e.what();
                if (w.find(" (") != std::string::npos && w.compare(0, kw.size(), kw) == 0) throw;
                fail(w);
            }
        }
        K = nullptr;
        for (auto& [name, cf] : unsupported) {
            K = cf.second;
            std::string t = cf.second->known ? "not modelled by the ResSim engine - kept verbatim for OPM Flow (Advanced keywords)"
                                             : "unknown keyword - kept verbatim (Advanced keywords)";
            if (cf.first > 1) t += " (" + std::to_string(cf.first) + " occurrences)";
            warn(t);
        }
        K = nullptr;
        if (!pending.empty()) warn("well controls after the last TSTEP/DATES are ignored");
    }

    // ------------------------------------------------------------------ saturation table conversion (family II -> I)
    static std::vector<double> col(const std::vector<double>& t, int ncol, int c) {
        std::vector<double> v; for (size_t i = size_t(c); i < t.size(); i += size_t(ncol)) v.push_back(t[i]); return v;
    }
    static std::vector<double> sortedUnique(std::vector<double> v, double lo, double hi) {
        std::vector<double> o;
        for (double x : v) if (x >= lo - 1e-12 && x <= hi + 1e-12) o.push_back(std::min(hi, std::max(lo, x)));
        std::sort(o.begin(), o.end());
        std::vector<double> u;
        for (double x : o) if (u.empty() || x - u.back() > 1e-9) u.push_back(x);
        return u;
    }
    SatRegion satTables(int t) {
        SatRegion R;
        auto rows = [](const std::vector<double>& v, std::vector<SatRow>& out) {
            out.clear();
            for (size_t i = 0; i + 4 <= v.size(); i += 4) out.push_back({v[i], v[i + 1], v[i + 2], v[i + 3]});
        };
        const size_t T = size_t(t);
        if (T < swof.size()) rows(swof[T], R.swof);
        if (T < sgof.size()) rows(sgof[T], R.sgof);
        if (R.sgof.empty() && T < slgof.size()) {
            const auto& v = slgof[T];
            for (size_t i = v.size(); i >= 4; i -= 4) R.sgof.push_back({1.0 - v[i - 4], v[i - 3], v[i - 2], v[i - 1]});
        }
        const bool fam2 = T < swfn.size() || T < sgfn.size();
        if (fam2) {
            auto kroOf = [&](double so, bool gasOil) -> double {
                if (T < sof3.size()) { const auto x = col(sof3[T], 3, 0); return lin(x, col(sof3[T], 3, gasOil ? 2 : 1), so); }
                if (T < sof2.size()) { const auto x = col(sof2[T], 2, 0); return lin(x, col(sof2[T], 2, 1), so); }
                return NaN;
            };
            const std::vector<double>* soNodes = T < sof3.size() ? &sof3[T] : T < sof2.size() ? &sof2[T] : nullptr;
            const int soCols = T < sof3.size() ? 3 : 2;
            double swco = 0;
            if (T < swfn.size() && R.swof.empty()) {
                const auto sw = col(swfn[T], 3, 0), krw = col(swfn[T], 3, 1), pc = col(swfn[T], 3, 2);
                swco = sw.front();
                std::vector<double> nodes = sw;
                if (soNodes) for (double so : col(*soNodes, soCols, 0)) nodes.push_back(1.0 - so);
                for (double s : sortedUnique(nodes, sw.front(), sw.back())) {
                    const double kro = kroOf(1.0 - s, false);
                    if (std::isnan(kro)) fail("SWFN needs SOF3 or SOF2 for the oil relative permeability");
                    R.swof.push_back({s, lin(sw, krw, s), kro, lin(sw, pc, s)});
                }
            } else if (!R.swof.empty()) swco = R.swof.front().s;
            if (T < sgfn.size() && R.sgof.empty()) {
                const auto sg = col(sgfn[T], 3, 0), krg = col(sgfn[T], 3, 1), pc = col(sgfn[T], 3, 2);
                std::vector<double> nodes = sg;
                if (soNodes) for (double so : col(*soNodes, soCols, 0)) nodes.push_back(1.0 - swco - so);
                for (double s : sortedUnique(nodes, sg.front(), std::min(sg.back(), 1.0 - swco))) {
                    const double kro = kroOf(1.0 - swco - s, true);
                    if (std::isnan(kro)) fail("SGFN needs SOF3 or SOF2 for the oil relative permeability");
                    R.sgof.push_back({s, lin(sg, krg, s), kro, lin(sg, pc, s)});
                }
            }
        }
        return R;
    }

    // ------------------------------------------------------------------ assemble the Case
    void finish() {
        K = nullptr;
        if (n == 0) throw std::runtime_error("Deck: DIMENS is missing");
        if (!c.phases.oil && !c.phases.water) throw std::runtime_error("Deck: no OIL or WATER phase keyword in RUNSPEC");
        GridInput& g = c.grid;

        // ---- geometry
        const bool cp = !g.zcorn.empty() || !g.coord.empty();
        auto has = [&](const char* name) { return arr.count(name) > 0; };
        auto take = [&](const char* name, std::vector<double>& dst) {
            usedArrays.insert(name);
            auto it = arr.find(name);
            if (it == arr.end()) { dst.clear(); return false; }
            dst = it->second; return true;
        };
        auto undefined = [&](const std::vector<double>& v) { size_t m = 0; for (double x : v) if (std::isnan(x)) ++m; return m; };
        usedArrays.insert({"DX", "DY", "DZ", "TOPS"});
        if (cp) {
            if (g.zcorn.empty() || g.coord.empty()) throw std::runtime_error("Deck: corner-point grids need both COORD and ZCORN");
            if (has("DX") || has("TOPS")) msg(0, "COORD/ZCORN given - DX/DY/DZ/TOPS ignored");
            g.dx.clear(); g.dy.clear(); g.dz.clear(); g.tops.clear();
        } else {
            for (const char* nm : {"DX", "DY", "DZ"}) {
                std::vector<double>& dst = nm[1] == 'X' ? g.dx : nm[1] == 'Y' ? g.dy : g.dz;
                if (!take(nm, dst)) throw std::runtime_error(std::string("Deck: ") + nm + " is missing (or COORD/ZCORN)");
                if (undefined(dst)) throw std::runtime_error(std::string("Deck: ") + nm + " is not defined for every cell");
            }
            auto it = arr.find("TOPS");
            if (it == arr.end()) throw std::runtime_error("Deck: TOPS is missing");
            std::vector<double>& t = it->second;
            for (int q = 0; q < nl; ++q) if (std::isnan(t[size_t(q)])) throw std::runtime_error("Deck: TOPS is missing in the top layer");
            bool onlyTop = true;
            for (int q = nl; q < n; ++q) if (!std::isnan(t[size_t(q)])) onlyTop = false;
            if (onlyTop) g.tops.assign(t.begin(), t.begin() + nl);
            else {
                for (int q = nl; q < n; ++q) if (std::isnan(t[size_t(q)])) t[size_t(q)] = t[size_t(q - nl)] + g.dz[size_t(q - nl)];
                g.tops = t;
            }
        }

        // ---- ACTNUM, PORO, PERM
        std::vector<double> act; take("ACTNUM", act);
        if (!take("PORO", g.poro)) throw std::runtime_error("Deck: PORO is missing");
        if (!take("PERMX", g.permx)) throw std::runtime_error("Deck: PERMX is missing");
        take("PERMY", g.permy); take("PERMZ", g.permz);
        if (g.permy.empty()) msg(0, "PERMY not given - PERMX is used");
        if (g.permz.empty()) msg(0, "PERMZ not given - PERMX is used");
        size_t madeInactive = 0;
        for (int q = 0; q < n; ++q) {
            const size_t u = size_t(q);
            const bool a = act.empty() || act[u] != 0;
            bool bad = std::isnan(g.poro[u]) || std::isnan(g.permx[u]) || (!g.permy.empty() && std::isnan(g.permy[u])) ||
                       (!g.permz.empty() && std::isnan(g.permz[u]));
            if (bad) {
                if (a) ++madeInactive;
                if (act.empty()) act.assign(size_t(n), 1.0);
                act[u] = 0;
                if (std::isnan(g.poro[u])) g.poro[u] = 0;
                if (std::isnan(g.permx[u])) g.permx[u] = 0;
                if (!g.permy.empty() && std::isnan(g.permy[u])) g.permy[u] = 0;
                if (!g.permz.empty() && std::isnan(g.permz[u])) g.permz[u] = 0;
            }
        }
        if (madeInactive) msg(1, std::to_string(madeInactive) + " cells without PORO/PERM values were made inactive");
        g.actnum.clear(); for (double v : act) g.actnum.push_back(v != 0 ? 1 : 0);
        take("NTG", g.ntg); take("MULTX", g.multx); take("MULTY", g.multy); take("MULTZ", g.multz); take("MULTPV", g.multpv);
        for (auto* v : {&g.ntg, &g.multx, &g.multy, &g.multz, &g.multpv}) for (double& x : *v) if (std::isnan(x)) x = 1.0;

        // ---- dual porosity
        take("SIGMAV", g.sigmav); take("LX", g.lx); take("LY", g.ly); take("LZ", g.lz); take("MULTMF", g.multmf);
        take("DZMTRXV", g.dzmtrxv);
        if (g.dual) {
            if (g.nz % 2) throw std::runtime_error("Deck: DUALPORO/DUALPERM needs an even NZ (matrix layers + fracture layers)");
            for (double& x : g.multmf) if (std::isnan(x)) x = 1.0;
            for (double& x : g.sigmav) if (std::isnan(x)) x = g.sigma;
            for (auto* v : {&g.lx, &g.ly, &g.lz, &g.dzmtrxv}) for (double& x : *v) if (std::isnan(x)) x = 0.0;
            const bool anySig = g.sigma > 0 || !g.sigmav.empty() || (g.ltosigma && !(g.lx.empty() && g.ly.empty() && g.lz.empty()));
            if (!anySig) msg(1, "dual porosity without SIGMA / SIGMAV / LTOSIGMA - matrix and fracture are not coupled");
            msg(0, std::string(g.dual == 1 ? "dual porosity (DUALPORO)" : "dual porosity / dual permeability (DUALPERM)") +
                   ": " + std::to_string(g.nz / 2) + " matrix + " + std::to_string(g.nz / 2) + " fracture layers, " +
                   (g.ltosigma ? "sigma from LX/LY/LZ (LTOSIGMA)" : !g.sigmav.empty() ? "SIGMAV per cell" : "SIGMA " + [&] { char b[32]; std::snprintf(b, sizeof b, "%g", g.sigma); return std::string(b); }()) +
                   (g.nodppm ? ", NODPPM" : ", fracture k x fracture porosity"));
            if (arr.count("DPNUM")) msg(1, "DPNUM (single porosity regions) is not modelled - every block is dual porosity");
        } else if (g.sigma > 0 || !g.sigmav.empty()) msg(1, "SIGMA/SIGMAV given without DUALPORO - ignored");

        // ---- MULTREGT -> MULTX/Y/Z
        if (!multregt.empty()) {
            for (auto* v : {&g.multx, &g.multy, &g.multz}) if (v->empty()) v->assign(size_t(n), 1.0);
            for (const MultRegT& m : multregt) {
                usedArrays.insert(m.set);
                const std::vector<double>& R = array(m.set);
                size_t touched = 0;
                auto apply = [&](int q1, int q2, std::vector<double>& mult) {
                    const int a = int(std::lround(R[size_t(q1)])), b = int(std::lround(R[size_t(q2)]));
                    if (a == b) return;
                    const bool match = ((m.r1 < 0 || a == m.r1) && (m.r2 < 0 || b == m.r2)) || ((m.r1 < 0 || b == m.r1) && (m.r2 < 0 || a == m.r2));
                    if (match) { mult[size_t(q1)] *= m.m; ++touched; }
                };
                for (int k = 0; k < g.nz; ++k)
                    for (int j = 0; j < g.ny; ++j)
                        for (int i = 0; i < g.nx; ++i) {
                            const int q = g.index(i, j, k);
                            if (m.dirs.find('X') != std::string::npos && i + 1 < g.nx) apply(q, g.index(i + 1, j, k), g.multx);
                            if (m.dirs.find('Y') != std::string::npos && j + 1 < g.ny) apply(q, g.index(i, j + 1, k), g.multy);
                            if (m.dirs.find('Z') != std::string::npos && k + 1 < g.nz) apply(q, g.index(i, j, k + 1), g.multz);
                        }
                msg(0, "MULTREGT " + std::to_string(m.r1) + "-" + std::to_string(m.r2) + " x" + fmtNum(m.m) + ": " +
                           std::to_string(touched) + " cell faces");
            }
            if (!g.faults.empty() || cp) msg(0, "MULTREGT is applied to neighbour connections only (not to fault NNCs)");
        }

        // ---- regions
        auto intArray = [&](const char* name, std::vector<int>& dst, int nmax, const char* what) {
            usedArrays.insert(name);
            dst.clear();
            auto it = arr.find(name);
            if (it == arr.end()) return;
            dst.resize(size_t(n));
            int mx = 1; bool all1 = true;
            for (int q = 0; q < n; ++q) {
                const double v = it->second[size_t(q)];
                dst[size_t(q)] = std::isnan(v) ? 1 : int(std::lround(v));
                if (g.actnum.empty() || g.actnum[size_t(q)]) { mx = std::max(mx, dst[size_t(q)]); if (dst[size_t(q)] != 1) all1 = false; }
                if (dst[size_t(q)] < 1) dst[size_t(q)] = 1;
            }
            if (all1) { dst.clear(); return; }
            if (nmax > 0 && mx > nmax)
                throw std::runtime_error(std::string("Deck: ") + name + " uses region " + std::to_string(mx) + " but only " +
                                         std::to_string(nmax) + " " + what + " are defined");
        };
        // PVT regions
        ensurePvt(ntpvt);
        const DensitySC defRho = defaultDensity();
        for (int t = 0; t < ntpvt; ++t) {
            const size_t u = size_t(t);
            if (!hasDensity[u]) pvt[u].density = t > 0 && hasDensity[u - 1] ? pvt[u - 1].density : defRho;
            if (!hasPvtw[u] && c.phases.water) {
                if (t > 0) pvt[u].pvtw = pvt[u - 1].pvtw;
                else { msg(1, "PVTW is missing - default water properties used"); pvt[u].pvtw.bw = 1.0; }
            }
            if (!hasRock[u]) { if (t > 0) pvt[u].rock = pvt[u - 1].rock; else msg(0, "ROCK is missing - incompressible rock"); }
        }
        intArray("PVTNUM", g.pvtnum, ntpvt, "PVT tables");
        intArray("SATNUM", g.satnum, ntsfun, "saturation tables");
        intArray("EQLNUM", g.eqlnum, ntequl, "equilibration regions");
        intArray("FIPNUM", g.fipnum, 0, "");
        int npvtUsed = 1, nsatUsed = 1, neqlUsed = 1;
        for (int v : g.pvtnum) npvtUsed = std::max(npvtUsed, v);
        for (int v : g.satnum) nsatUsed = std::max(nsatUsed, v);
        for (int v : g.eqlnum) neqlUsed = std::max(neqlUsed, v);
        // dissolved gas sanity
        if (c.phases.disgas && !(c.phases.oil && c.phases.gas)) { c.phases.disgas = false; msg(1, "DISGAS needs OIL and GAS - ignored"); }
        if (c.phases.oil && c.phases.gas && !c.phases.disgas && hasPvto[0]) {
            msg(1, "PVTO given without DISGAS - live oil enabled"); c.phases.disgas = true;
        }
        if (c.phases.disgas && !hasPvto[0] && hasPvdo[0]) { c.phases.disgas = false; msg(1, "DISGAS without PVTO - dead oil (PVDO) is used"); }
        if (c.phases.vapoil && !(c.phases.oil && c.phases.gas)) { c.phases.vapoil = false; msg(1, "VAPOIL needs OIL and GAS - ignored"); }
        if (c.phases.vapoil && pvt[0].pvtg.empty()) { c.phases.vapoil = false; msg(1, "VAPOIL without PVTG - dry gas (PVDG) is used"); }
        if (!c.phases.vapoil && !pvt[0].pvtg.empty()) msg(1, "PVTG without VAPOIL - vaporised oil is off, the saturated gas FVF / viscosity are used");
        for (int t = 0; t < npvtUsed; ++t) {
            const size_t u = size_t(t);
            const std::string reg = npvtUsed > 1 ? " for PVT region " + std::to_string(t + 1) : "";
            if (c.phases.oil && !c.phases.disgas && pvt[u].pvdo.empty()) throw std::runtime_error("Deck: PVDO (or PVCDO) is missing" + reg);
            if (c.phases.disgas && pvt[u].pvto.empty()) throw std::runtime_error("Deck: PVTO is missing" + reg);
            if (c.phases.gas && pvt[u].pvdg.empty()) throw std::runtime_error("Deck: PVDG (or PVTG) is missing" + reg);
            if (c.phases.vapoil && pvt[u].pvtg.size() < 2) throw std::runtime_error("Deck: PVTG with at least 2 pressures is needed (VAPOIL)" + reg);
        }
        PropsInput& P = c.props;
        const PvtRegion& p0 = pvt[0];
        P.pvtw = p0.pvtw; P.rock = p0.rock; P.density = p0.density; P.pvdo = p0.pvdo; P.pvdg = p0.pvdg; P.pvto = p0.pvto; P.pvtg = p0.pvtg;
        P.pvtExtra.assign(pvt.begin() + 1, pvt.begin() + npvtUsed);
        // saturation regions
        std::vector<SatRegion> sats;
        for (int t = 0; t < nsatUsed; ++t) {
            SatRegion R = satTables(t);
            const std::string reg = nsatUsed > 1 ? " for saturation region " + std::to_string(t + 1) : "";
            if (c.phases.water && c.phases.oil && R.swof.size() < 2) throw std::runtime_error("Deck: SWOF (or SWFN + SOF3/SOF2) is missing" + reg);
            if (c.phases.gas && c.phases.oil && R.sgof.size() < 2) throw std::runtime_error("Deck: SGOF (or SGFN/SLGOF) is missing" + reg);
            sats.push_back(R);
        }
        if (!swfn.empty() || !sgfn.empty() || !slgof.empty()) msg(0, "saturation functions converted to SWOF/SGOF form");
        P.swof = sats[0].swof; P.sgof = sats[0].sgof;
        P.satExtra.assign(sats.begin() + 1, sats.end());
        if (ntpvt > npvtUsed) msg(0, std::to_string(ntpvt - npvtUsed) + " PVT table(s) not referenced by PVTNUM - dropped");
        if (ntsfun > nsatUsed) msg(0, std::to_string(ntsfun - nsatUsed) + " saturation table(s) not referenced by SATNUM - dropped");

        // ---- initial conditions
        InitialState& I = c.init;
        const bool enumerated = arr.count("PRESSURE") > 0;
        usedArrays.insert({"PRESSURE", "SWAT", "SGAS", "RS", "RV"});
        if (sawEquil) {
            if (enumerated) msg(1, "both EQUIL and PRESSURE given - equilibration (EQUIL) is used");
            while (int(equil.size()) < neqlUsed) { equil.push_back(equil.back()); equilItem7.push_back(equilItem7.back()); equilItem8.push_back(equilItem8.back()); }
            for (size_t e = 0; e < equil.size(); ++e) {
                EquilInput& E = equil[e];
                if (c.phases.vapoil) {
                    // item 8 <= 0 -> Rv saturated at the GOC pressure; > 0 -> RVVD table (a defaulted item with RVVD = table)
                    std::vector<double> dz, rv;
                    if (e < rvvd.size()) { const auto& t = rvvd[e]; for (size_t i = 0; i + 1 < t.size(); i += 2) { dz.push_back(t[i]); rv.push_back(t[i + 1]); } }
                    const int it8 = equilItem8[e];
                    const bool useTable = !dz.empty() && (it8 > 0 || it8 == -99);
                    if (!dz.empty() && !useTable) msg(0, "EQUIL item 8 = " + std::to_string(it8) + " for region " + std::to_string(e + 1) + ": RVVD ignored, Rv saturated at the GOC");
                    if (useTable) { E.rvvdDepth = dz; E.rvvdRv = rv; E.initialRv = -1; }
                }
                if (!c.phases.disgas) continue;
                std::vector<double> dz, rs;
                if (e < rsvd.size()) { const auto& t = rsvd[e]; for (size_t i = 0; i + 1 < t.size(); i += 2) { dz.push_back(t[i]); rs.push_back(t[i + 1]); } }
                else if (e < pbvd.size()) {
                    // PBVD -> RSVD through the saturated PVTO curve of PVT region 1
                    std::vector<double> pb, rsat;
                    for (const auto& r : P.pvto) { pb.push_back(r.pb); rsat.push_back(r.rs); }
                    const auto& t = pbvd[e];
                    for (size_t i = 0; i + 1 < t.size(); i += 2) { dz.push_back(t[i]); rs.push_back(lin(pb, rsat, t[i + 1])); }
                    msg(0, "PBVD converted to RSVD with the saturated PVTO curve");
                }
                // ECLIPSE: item 7 <= 0 -> Rs saturated at the GOC (tables ignored); > 0 -> RSVD/PBVD table. A defaulted
                // item 7 with tables present is taken as "use the table" (common in hand-edited decks).
                const int it7 = equilItem7[e];
                const bool useTable = !dz.empty() && (it7 > 0 || (it7 == -99 && rsvd.size() + pbvd.size() > 0));
                if (!dz.empty() && !useTable && it7 != -99) msg(0, "EQUIL item 7 = " + std::to_string(it7) + " for region " + std::to_string(e + 1) + ": RSVD/PBVD ignored, Rs saturated at the GOC");
                if (useTable) { E.rsvdDepth = dz; E.rsvdRs = rs; E.initialRs = -1; }
            }
            c.equil = equil[0];
            c.equilExtra.assign(equil.begin() + 1, equil.begin() + neqlUsed);
            if (int(equil.size()) > neqlUsed) msg(0, std::to_string(equil.size() - size_t(neqlUsed)) + " EQUIL region(s) not referenced by EQLNUM - dropped");
        } else if (enumerated) {
            I.pressure = arr["PRESSURE"];
            if (undefined(I.pressure)) {
                size_t m = 0;
                for (int q = 0; q < n; ++q) if (std::isnan(I.pressure[size_t(q)]) && (g.actnum.empty() || g.actnum[size_t(q)])) ++m;
                if (m) throw std::runtime_error("Deck: PRESSURE is not defined for " + std::to_string(m) + " active cells");
                for (double& x : I.pressure) if (std::isnan(x)) x = 0;
            }
            auto opt = [&](const char* name, std::vector<double>& dst, double def) {
                auto it = arr.find(name);
                if (it == arr.end()) { dst.clear(); return; }
                dst = it->second; for (double& x : dst) if (std::isnan(x)) x = def;
            };
            opt("SWAT", I.swat, 0.0); opt("SGAS", I.sgas, 0.0); opt("RS", I.rs, 0.0);
            if (c.phases.vapoil) opt("RV", I.rv, 0.0); else I.rv.clear();
            if (c.phases.water && c.phases.oil && I.swat.empty()) msg(1, "SWAT missing - connate water saturation used");
            if (c.phases.disgas && I.rs.empty()) msg(1, "RS missing - saturated oil assumed");
            msg(0, "enumerated initial conditions (PRESSURE/SWAT/SGAS/RS)");
            // keep a sensible EQUIL for the UI
            double sp = 0, sd = 0; int cnt = 0;
            for (int q = 0; q < n && cnt < 1; ++q) if (I.pressure[size_t(q)] > 0) { sp = I.pressure[size_t(q)]; sd = 0; ++cnt; }
            (void)sd; c.equil.datumPressure = sp;
        } else {
            throw std::runtime_error("Deck: no initial conditions (EQUIL or PRESSURE/SWAT) in SOLUTION");
        }
        if (arr.count("SWATINIT")) {
            usedArrays.insert("SWATINIT");
            if (c.init.active()) msg(1, "SWATINIT ignored - the initial state is enumerated (PRESSURE / SWAT)");
            else if (!(c.phases.water && c.phases.oil)) msg(1, "SWATINIT ignored - needs water and oil");
            else {
                c.grid.swatinit = arr.at("SWATINIT");
                for (double& x : c.grid.swatinit) if (std::isnan(x)) x = -1;
                msg(0, "SWATINIT: water-oil capillary pressure scaled per cell to honour the initial water saturation");
            }
        }

        if (arr.count("PLMIXNUM")) {
            usedArrays.insert("PLMIXNUM");
            if (c.phases.polymer) {
                c.grid.plmixnum.clear();
                for (double x : arr.at("PLMIXNUM")) c.grid.plmixnum.push_back(std::isnan(x) ? 1 : int(x + 0.5));
                int mx = 1; for (int v : c.grid.plmixnum) mx = std::max(mx, v);
                if (mx > c.polymer.numMixRegions()) msg(2, "PLMIXNUM refers to region " + std::to_string(mx) + " but PLYMAX / PLMIXPAR have " + std::to_string(c.polymer.numMixRegions()) + " record(s)");
            }
        }
        if (arr.count("SPOLY")) {
            usedArrays.insert("SPOLY");
            if (!c.phases.polymer) msg(1, "SPOLY given without POLYMER - ignored");
            else { c.grid.spoly = arr.at("SPOLY"); for (double& x : c.grid.spoly) if (std::isnan(x)) x = 0; }
        }
        if (c.phases.polymer) {
            if (!c.phases.water) { c.phases.polymer = false; msg(1, "POLYMER needs WATER - ignored"); }
            else {
                if (c.polymer.plyvisc.empty() || c.polymer.plyads.empty() || c.polymer.plyrock.empty() || !(c.polymer.cmax[0] > 0))
                    msg(2, "POLYMER needs PLYVISC, PLYADS, PLYROCK and PLYMAX");
                msg(0, "polymer flooding: Todd-Longstaff viscosity mixing, adsorption, RRF and dead pore volume modelled");
            }
        }
        // ---- saturation end-point scaling
        {
            static const char* epsNames[] = {"SWL", "SWCR", "SWU", "SGL", "SGCR", "SGU", "SOWCR", "SOGCR"};
            std::vector<double>* dst[] = {&c.grid.swl, &c.grid.swcr, &c.grid.swu, &c.grid.sgl, &c.grid.sgcr, &c.grid.sgu, &c.grid.sowcr, &c.grid.sogcr};
            std::string found;
            for (int i = 0; i < 8; ++i) {
                if (!arr.count(epsNames[i])) continue;
                usedArrays.insert(epsNames[i]);
                found += std::string(found.empty() ? "" : ", ") + epsNames[i];
                if (!c.grid.endscale) continue;
                *dst[i] = arr.at(epsNames[i]);
                for (double& x : *dst[i]) if (std::isnan(x)) x = -1;
            }
            if (!found.empty() && !c.grid.endscale) msg(1, found + " given without ENDSCALE in RUNSPEC - end-point scaling is off (as in ECLIPSE / OPM Flow)");
            else if (!found.empty()) msg(0, "end-point scaling (two-point, kr and Pc): " + found);
            if (arr.count("SCALECRS")) msg(1, "SCALECRS: three-point scaling is not modelled - two-point scaling is used");
        }

        // ---- aquifers
        for (auto& [id, cs] : aqConn) {
            bool found = false;
            for (auto& a : c.aquifers) if (a.id == id) { a.conns = cs; found = true; }
            if (!found) msg(1, "AQUANCON: aquifer " + std::to_string(id) + " is not defined (AQUFETP/AQUCT) - connections ignored");
        }
        std::vector<AquiferInput> keep;
        for (auto& a : c.aquifers) {
            if (aqVisc.count(a.id)) a.visc = aqVisc[a.id];
            if (a.conns.empty()) msg(1, "aquifer " + std::to_string(a.id) + " has no AQUANCON connections - ignored");
            else keep.push_back(a);
        }
        c.aquifers = keep;

        // ---- arrays that were read but have no effect
        static const std::set<std::string> silentArrays = {"MAPAXES", "TRAN-EDIT"};
        std::vector<std::string> unusedNonTrivial;
        for (auto& [name, v] : arr) {
            if (usedArrays.count(name) || silentArrays.count(name) || name.compare(0, 3, "FIP") == 0) continue;
            if (name == "MULTNUM" || name == "FLUXNUM" || name == "OPERNUM") continue;   // used by *REG operations
            const double d = arrayDefault(name);
            bool trivial = !std::isnan(d);
            if (trivial) for (double x : v) if (!std::isnan(x) && x != d) { trivial = false; break; }
            if (!trivial) unusedNonTrivial.push_back(name);
        }
        if (!unusedNonTrivial.empty()) {
            std::string list;
            for (auto& s : unusedNonTrivial) list += (list.empty() ? "" : ", ") + s;
            msg(1, "arrays not used by the simulator: " + list);
        }

        if (c.schedule.empty()) throw std::runtime_error("Deck: SCHEDULE has no TSTEP or DATES");
        // a well has one type for the whole run: controls of the other type (e.g. a producer that is converted to
        // an injector later) are turned into shut-ins
        {
            std::map<std::string, int> converted;
            for (size_t si = 0; si < c.schedule.size(); ++si)
                for (auto& wc : c.schedule[si].controls) {
                    const WellInput* w = nullptr;
                    for (auto& x : c.wells) if (x.name == wc.well) w = &x;
                    if (!w) continue;
                    const bool prodWell = w->type == WellType::Producer;
                    const bool prodCtrl = !injCtrl.count({si, wc.well});
                    if (prodWell == prodCtrl) continue;
                    const bool zero = prodCtrl ? (wc.mode != Control::BHP && wc.orat <= 0 && wc.wrat <= 0 && wc.grat <= 0 && wc.lrat <= 0) : wc.rate <= 0;
                    if (!zero && wc.status == WellStatus::Open) ++converted[wc.well];
                    wc.status = WellStatus::Shut;
                    if (!prodWell) { wc.mode = Control::RATE; wc.rate = 0; } else { wc.mode = Control::BHP; }
                }
            for (auto& [name, cnt] : converted)
                msg(1, "well " + name + " changes between producer and injector - the periods as " +
                       std::string("the other type are simulated as shut (") + std::to_string(cnt) + " control records)");
        }
        const std::string err = validateCase(c);
        if (!err.empty()) msg(1, "validation: " + err);

        // flush grouped messages
        for (const std::string& key : groupOrder) {
            const Group& gr = groups[key];
            DeckMessage m; m.level = gr.level; m.keyword = gr.keyword; m.file = gr.file; m.line = gr.line;
            m.text = key.substr(key.find('\x1f') + 1);
            if (!gr.items.empty()) {
                std::string l; size_t shown = 0;
                for (const auto& it : gr.items) { if (shown == 8) break; l += (l.empty() ? "" : ", ") + it; ++shown; }
                if (gr.items.size() > shown) l += ", ... (" + std::to_string(gr.items.size()) + " in total)";
                m.text += ": " + l;
            } else if (gr.count > 1) m.text += " (" + std::to_string(gr.count) + " times)";
            msgs.push_back(m);
        }
        groups.clear(); groupOrder.clear();
    }
};

void deliver(const deck::RawDeck& D, std::vector<DeckMessage>& own, std::vector<std::string>* warnings, std::vector<DeckMessage>* report) {
    std::vector<DeckMessage> all;
    for (const auto& m : D.messages) {
        DeckMessage d; d.level = m.level == deck::Level::Info ? 0 : m.level == deck::Level::Warning ? 1 : 2;
        d.keyword = m.keyword; d.file = (m.file >= 0 && m.file < int(D.files.size())) ? D.files[size_t(m.file)] : "";
        d.line = m.line; d.text = m.text;
        all.push_back(d);
    }
    all.insert(all.end(), own.begin(), own.end());
    if (warnings)
        for (const auto& m : all) {
            if (m.level < 1) continue;
            std::string f = m.file; const size_t s = f.find_last_of("/\\"); if (s != std::string::npos) f = f.substr(s + 1);
            warnings->push_back(m.keyword + (f.empty() || m.line <= 0 ? "" : " (" + f + ":" + std::to_string(m.line) + ")") + ": " + m.text);
        }
    if (report) report->insert(report->end(), all.begin(), all.end());
}

Case importRaw(const deck::RawDeck& D, std::vector<std::string>* warnings, std::vector<DeckMessage>* report) {
    Importer im(D);
    try {
        im.run();
        im.finish();
    } catch (const std::exception& e) {
        DeckMessage m; m.level = 2; m.keyword = "Deck"; m.text = e.what();
        im.msgs.push_back(m);
        deliver(D, im.msgs, nullptr, report);
        throw;
    }
    deliver(D, im.msgs, warnings, report);
    return im.c;
}

}  // namespace

// ======================================================================== public: keyword support
namespace deck {
Support ressimSupport(const std::string& keyword, std::string* note) {
    auto it = supportTable().find(keyword);
    if (it != supportTable().end()) { if (note) *note = it->second.note; return it->second.s; }
    const KeywordSpec* s = findKeyword(keyword);
    if (s && s->sections == 64) { if (note) *note = "summary output request (all vectors are computed internally)"; return Support::Ignored; }
    if (note) *note = s ? "not simulated by the ResSim engine - kept on import and written back for OPM Flow" : "not an OPM Flow / ECLIPSE keyword";
    return s ? Support::Unsupported : Support::Unknown;
}
const char* toString(Support s) {
    switch (s) {
        case Support::Full: return "supported";
        case Support::Partial: return "partly supported";
        case Support::Ignored: return "no effect";
        case Support::Unsupported: return "not supported";
        default: return "unknown";
    }
}
}  // namespace deck

// ======================================================================== public: import
Case eclipseDeckFromString(const std::string& text, const std::string& baseDir, std::vector<std::string>* warnings,
                           std::vector<DeckMessage>* report) {
    std::string dir = baseDir;
    if (!dir.empty() && dir.back() != '/' && dir.back() != '\\') dir += '/';
    const deck::RawDeck D = deck::parseDeckText(text, dir, "deck");
    return importRaw(D, warnings, report);
}

Case importEclipseDeck(const std::string& path, std::vector<std::string>* warnings, std::vector<DeckMessage>* report) {
    const deck::RawDeck D = deck::parseDeckFile(path);
    Case c = importRaw(D, warnings, report);
    return c;
}

// ======================================================================== public: export
std::string eclipseDeckToString(const Case& c) {
    std::ostringstream o;
    // extra (verbatim) keywords of a section / of a report step, in their original order
    auto extraText = [&](const ExtraKeyword& e) {
        std::string b = e.body;
        if (!b.empty() && b.back() != '\n') b += '\n';
        return e.name + "\n" + b + "\n";
    };
    auto writeExtra = [&](DeckSection sec, const char* header = nullptr, const std::set<std::string>* skip = nullptr) {
        bool any = false;
        for (auto& e : c.extraKeywords) {
            if (e.section != sec || (skip && skip->count(e.name)) || isDimsMinimum(e)) continue;
            if (!any && header) o << header << "\n\n";
            any = true;
            o << extraText(e);
        }
    };
    // before: keywords that preceded the well controls of the step in the source deck (VFP tables, groups, ...)
    auto extraAtStep = [&](int s, bool last, int before = -1) {
        std::string out;
        for (auto& e : c.extraKeywords)
            if (e.section == DeckSection::Schedule && (e.step == s || (last && e.step > s)) && (before < 0 || e.beforeControls == (before == 1)))
                out += extraText(e);
        return out;
    };
    const GridInput& g = c.grid;
    const PropsInput& p = c.props;
    const bool live = c.phases.gas && c.phases.oil && c.phases.disgas;
    const bool cp = g.cornerPoint();
    o << "-- ECLIPSE 100 deck generated by ResSim Studio\n";
    o << "-- Lines starting with '-- RESSIM:' carry settings ECLIPSE has no keyword for (ignored by ECLIPSE).\n\n";
    o << "RUNSPEC\n\nTITLE\n" << (c.title.empty() ? std::string("Untitled") : c.title) << "\n\n";
    o << "DIMENS\n " << g.nx << ' ' << g.ny << ' ' << g.nz << " /\n\n";
    if (c.phases.oil) o << "OIL\n";
    if (c.phases.water) o << "WATER\n";
    if (c.phases.gas) o << "GAS\n";
    if (live) o << "DISGAS\n";
    const bool vap = c.phases.gas && c.phases.oil && c.phases.vapoil;
    if (vap) o << "VAPOIL\n";
    if (c.grid.dual == 1) o << "DUALPORO\n";
    if (c.grid.dual == 2) o << "DUALPERM\n";
    const bool polyOn = c.phases.polymer && c.phases.water;
    if (polyOn) o << "POLYMER\n";
    if (g.hasEndPointScaling()) o << "ENDSCALE\n NODIR REVERS /\n";
    o << "\n" << (c.units == UnitSystem::Field ? "FIELD" : "METRIC") << "\n\n";
    const int nsat = p.numSatRegions(), npvt = p.numPvtRegions(), neql = 1 + int(c.equilExtra.size());
    size_t rows = 20, pres = 20, rsvdRows = 20;
    for (int r = 0; r < nsat; ++r) { const SatRegion S = p.satRegion(r); rows = std::max({rows, S.swof.size(), S.sgof.size()}); }
    if (polyOn) { for (auto& t : c.polymer.plyvisc) rows = std::max(rows, t.size() / 2); for (auto& t : c.polymer.plyads) rows = std::max(rows, t.size() / 2); }
    for (int r = 0; r < npvt; ++r) {
        const PvtRegion R = p.pvtRegion(r);
        rows = std::max({rows, R.pvdo.size(), R.pvdg.size(), R.pvto.size()});
        for (auto& x : R.pvto) pres = std::max(pres, x.under.size() + 1);
        if (vap) { pres = std::max(pres, R.pvtg.size()); for (auto& x : R.pvtg) rows = std::max(rows, x.rows.size()); }
    }
    for (int e = 0; e < neql; ++e) rsvdRows = std::max(rsvdRows, (e == 0 ? c.equil : c.equilExtra[size_t(e - 1)]).rsvdDepth.size());
    for (int e = 0; e < neql; ++e) rsvdRows = std::max(rsvdRows, (e == 0 ? c.equil : c.equilExtra[size_t(e - 1)]).rvvdDepth.size());
    // dimension keyword = max(model needs, imported minimum) item by item; 0 = defaulted (1*)
    // items that COUNT the tables written (TABDIMS NTSFUN/NTPVT, EQLDIMS NTEQUL) must match the model exactly
    auto writeDims = [&](const char* name, std::vector<long> need, bool always) {
        bool have = always;
        const size_t exact = std::strcmp(name, "TABDIMS") == 0 ? 2 : std::strcmp(name, "EQLDIMS") == 0 ? 1 : 0;
        for (auto& e : c.extraKeywords)
            if (isDimsMinimum(e) && e.name == name) {
                const std::vector<long> m = dimsItems(e.body);
                if (m.size() > need.size()) need.resize(m.size(), 0);
                for (size_t i = exact; i < m.size(); ++i) need[i] = std::max(need[i], m[i]);
                have = true;
            }
        if (!have) return;
        while (!need.empty() && need.back() <= 0) need.pop_back();
        o << name << "\n";
        for (long x : need) { o << ' '; if (x > 0) o << x; else o << "1*"; }
        o << " /\n\n";
    };
    writeDims("TABDIMS", {nsat, npvt, long(rows), long(pres)}, true);
    {
        std::vector<long> rd(10, 0);
        if (polyOn && c.polymer.numMixRegions() > 1) rd[9] = c.polymer.numMixRegions();
        writeDims("REGDIMS", rd, polyOn && c.polymer.numMixRegions() > 1);
    }
    writeDims("EQLDIMS", {neql, 0, long(rsvdRows)}, true);
    if (!c.thpres.empty()) o << "EQLOPTS\n 'THPRES' /\n\n";
    size_t maxConn = 1; for (auto& w : c.wells) { size_t m = 0; for (auto& q : w.completions) m += size_t(q.k2 - q.k1 + 1); maxConn = std::max(maxConn, m); }
    std::map<std::string, size_t> grp;   // wells per group (groups as written in WELSPECS)
    for (auto& w : c.wells) ++grp[w.group.empty() ? std::string("G") : w.group];
    size_t maxPerGroup = 1; for (auto& kv : grp) maxPerGroup = std::max(maxPerGroup, kv.second);
    // item 4: max children of a group - wells of a group, or the groups directly under FIELD
    // groups also come from GRUPTREE records kept verbatim
    std::set<std::string> groups;
    for (auto& kv : grp) groups.insert(kv.first);
    for (auto& e : c.extraKeywords)
        if (e.name == "GRUPTREE") {
            std::istringstream ss(e.body); std::string line;
            while (std::getline(ss, line)) {
                std::vector<std::string> t; std::string cur; bool q = false;
                for (char ch : line) {
                    if (ch == '\'') { q = !q; continue; }
                    if (!q && (ch == ' ' || ch == '\t' || ch == '/')) { if (!cur.empty()) t.push_back(cur); cur.clear(); if (ch == '/') break; continue; }
                    cur += ch;
                }
                if (!cur.empty()) t.push_back(cur);
                if (t.size() >= 2 && t[0].rfind("--", 0) != 0) { groups.insert(upper(t[0])); if (upper(t[1]) != "FIELD") groups.insert(upper(t[1])); }
            }
        }
    groups.erase("FIELD");
    writeDims("WELLDIMS", {long(std::max<size_t>(1, c.wells.size())), long(maxConn), long(std::max<size_t>(1, groups.size()) + 1),
                           long(std::max(maxPerGroup, groups.size()))}, true);
    {
        size_t nc = 1;
        for (auto& a : c.aquifers) { size_t m = 0; for (auto& q : a.conns) m += size_t(q.i2 - q.i1 + 1) * (q.j2 - q.j1 + 1) * (q.k2 - q.k1 + 1); nc = std::max(nc, m); }
        writeDims("AQUDIMS", {0, 0, 0, 0, long(c.aquifers.size()), long(nc)}, !c.aquifers.empty());
    }
    { size_t nf = 0; for (auto& f : g.faults) nf += f.faces.size(); writeDims("FAULTDIM", {long(std::max<size_t>(1, nf))}, !g.faults.empty()); }
    o << "START\n " << c.startDay << " '" << kMonths[std::max(0, std::min(11, c.startMonth - 1))] << "' " << c.startYear << " /\n\n";
    o << "UNIFOUT\n\n";

    // ---------------- GRID
    writeExtra(DeckSection::Runspec);
    o << "GRID\n\n";
    auto arr = [&](const char* name, const auto& v) {
        if (v.empty()) return;
        o << name << "\n"; writeArray(o, v); o << "/\n\n";
    };
    if (cp) {
        o << "SPECGRID\n " << g.nx << ' ' << g.ny << ' ' << g.nz << " 1 F /\n\n";
        o << "COORD\n";
        for (size_t i = 0; i < g.coord.size(); i += 6)
            o << ' ' << num(g.coord[i]) << ' ' << num(g.coord[i + 1]) << ' ' << num(g.coord[i + 2]) << "  " << num(g.coord[i + 3])
              << ' ' << num(g.coord[i + 4]) << ' ' << num(g.coord[i + 5]) << "\n";
        o << "/\n\n";
        arr("ZCORN", g.zcorn);
    } else {
        arr("DX", g.dx); arr("DY", g.dy); arr("DZ", g.dz); arr("TOPS", g.tops);
    }
    arr("PORO", g.poro); arr("PERMX", g.permx); arr("PERMY", g.permy.empty() ? g.permx : g.permy);
    arr("PERMZ", g.permz.empty() ? g.permx : g.permz); arr("NTG", g.ntg); arr("ACTNUM", g.actnum);
    arr("MULTX", g.multx); arr("MULTY", g.multy); arr("MULTZ", g.multz); arr("MULTPV", g.multpv);
    if (g.minpv > 0) o << "MINPV\n " << num(g.minpv) << " /\n\n";
    // OLDTRAN is the ECLIPSE default for block-centred grids (and OPM Flow rejects the keyword), so only NEWTRAN is written
    if (!g.cornerPoint() && g.newtran) o << "NEWTRAN\n\n";
    if (g.pinch.on)
        o << "PINCH\n " << num(g.pinch.threshold) << (g.pinch.gap ? " GAP " : " NOGAP ") << (g.pinch.maxGap < 1e19 ? num(g.pinch.maxGap) : std::string("1*"))
          << (g.pinch.topbot ? " TOPBOT " : " ALL ") << (g.pinch.multzAll ? "ALL" : "TOP") << " /\n\n";
    if (g.dual) {
        if (g.sigma > 0) o << "SIGMA\n " << num(g.sigma) << " /\n\n";
        if (!g.sigmav.empty()) arr("SIGMAV", g.sigmav);
        if (g.ltosigma) o << "LTOSIGMA\n " << num(g.ltoF[0]) << ' ' << num(g.ltoF[1]) << ' ' << num(g.ltoF[2]) << " /\n\n";
        if (!g.lx.empty()) arr("LX", g.lx);
        if (!g.ly.empty()) arr("LY", g.ly);
        if (!g.lz.empty()) arr("LZ", g.lz);
        if (!g.multmf.empty()) arr("MULTMF", g.multmf);
        if (!g.dzmtrxv.empty()) arr("DZMTRXV", g.dzmtrxv);
        if (g.nodppm) o << "NODPPM\n\n";
    }
    if (!g.faults.empty()) {
        o << "FAULTS\n";
        for (auto& f : g.faults)
            for (auto& ff : f.faces)
                o << " '" << f.name << "' " << ff.i1 + 1 << ' ' << ff.i2 + 1 << ' ' << ff.j1 + 1 << ' ' << ff.j2 + 1 << ' ' << ff.k1 + 1 << ' '
                  << ff.k2 + 1 << " '" << ff.face << (ff.minus ? "-" : "") << "' /\n";
        o << "/\n\nMULTFLT\n";
        for (auto& f : g.faults) o << " '" << f.name << "' " << num(f.mult) << " /\n";
        o << "/\n\n";
    }
    if (!g.nnc.empty()) {
        o << "NNC\n";
        for (auto& e : g.nnc)
            o << ' ' << e.i1 + 1 << ' ' << e.j1 + 1 << ' ' << e.k1 + 1 << "  " << e.i2 + 1 << ' ' << e.j2 + 1 << ' ' << e.k2 + 1 << "  " << num(e.trans) << " /\n";
        o << "/\n\n";
    }
    // aquifer datum depths (needed by ECLIPSE): area-weighted connection depth when not given
    std::vector<double> aqDepth;
    for (auto& a : c.aquifers) aqDepth.push_back(a.depth);
    if (!c.aquifers.empty()) {
        bool need = false; for (double d : aqDepth) if (!(d > 0)) need = true;
        if (need) {
            try {
                Grid G; G.build(c);
                for (size_t ia = 0; ia < c.aquifers.size(); ++ia) {
                    if (aqDepth[ia] > 0) continue;
                    double sd = 0, sa = 0;
                    for (auto& q : c.aquifers[ia].conns)
                        for (int k = q.k1; k <= q.k2; ++k) for (int j = q.j1; j <= q.j2; ++j) for (int i = q.i1; i <= q.i2; ++i) {
                            const int gi = g.index(i, j, k); const int a = G.activeOf[gi];
                            if (a < 0) continue;
                            const double A = G.faceArea(gi, q.face, q.minus) * q.mult;
                            sd += A * G.depth[a]; sa += A;
                        }
                    if (sa > 0) aqDepth[ia] = sd / sa;
                }
            } catch (...) {}
        }
        o << "AQUANCON\n";
        for (auto& a : c.aquifers)
            for (auto& q : a.conns) {
                const char f = q.face == 'X' ? 'I' : q.face == 'Y' ? 'J' : 'K';
                o << ' ' << a.id << ' ' << q.i1 + 1 << ' ' << q.i2 + 1 << ' ' << q.j1 + 1 << ' ' << q.j2 + 1 << ' ' << q.k1 + 1 << ' ' << q.k2 + 1
                  << " '" << f << (q.minus ? '-' : '+') << "' 1* " << num(q.mult) << " 'NO' /\n";
            }
        o << "/\n\n";
    }
    o << "INIT\n\n";

    // ---------------- PROPS
    writeExtra(DeckSection::Grid);
    writeExtra(DeckSection::Edit, "EDIT");
    o << "PROPS\n\n";
    std::vector<PvtRegion> PR; for (int r = 0; r < npvt; ++r) PR.push_back(p.pvtRegion(r));
    std::vector<SatRegion> SR; for (int r = 0; r < nsat; ++r) SR.push_back(p.satRegion(r));
    o << "PVTW\n";
    for (auto& R : PR) o << ' ' << num(R.pvtw.pref) << ' ' << num(R.pvtw.bw) << ' ' << num(R.pvtw.cw) << ' ' << num(R.pvtw.visc) << ' ' << num(R.pvtw.cv) << " /\n";
    o << "\nROCK\n";
    for (auto& R : PR) o << ' ' << num(R.rock.pref) << ' ' << num(R.rock.cr) << " /\n";
    o << "\nDENSITY\n";
    for (auto& R : PR) o << ' ' << num(R.density.oil) << ' ' << num(R.density.water) << ' ' << num(R.density.gas) << " /\n";
    o << "\n";
    auto table3 = [&](const char* kw, auto get) {          // one PVDO/PVDG table per PVT region
        // OPM Flow / ECLIPSE need strictly increasing pressure and strictly decreasing B: repair flat / reversed rows
        // (OPM: PVDG Bg strictly, PVDO Bo non-strictly)
        const bool strict = kw[3] == 'G';
        o << kw << "\n";
        for (auto& R : PR) {
            double lp = -1e300, lb = 1e300; int fixed = 0;
            std::string rows;
            for (auto& r : get(R)) {
                if (r.p <= lp) { ++fixed; continue; }
                double B = r.B;
                if (strict ? B >= lb : B > lb) { B = strict ? lb * (1.0 - 1e-6) : lb; ++fixed; }
                rows += " " + num(r.p) + ' ' + num(B) + ' ' + num(r.mu) + "\n";
                lp = r.p; lb = B;
            }
            if (fixed) o << "-- " << fixed << " row(s) adjusted by ResSim: pressure must increase and B strictly decrease\n";
            o << rows << "/\n";
        }
        o << "\n";
    };
    if (c.phases.oil) {
        if (live) {
            o << "PVTO\n";
            for (size_t t = 0; t < PR.size(); ++t) {
                if (t) o << "-- PVT region " << t + 1 << "\n";
                for (auto& r : PR[t].pvto) {
                    o << " " << num(r.rs) << "  " << num(r.pb) << ' ' << num(r.bo) << ' ' << num(r.mu);
                    for (auto& u : r.under) o << "\n        " << num(u.p) << ' ' << num(u.B) << ' ' << num(u.mu);
                    o << " /\n";
                }
                o << "/\n";
            }
            o << "\n";
        } else table3("PVDO", [](const PvtRegion& R) -> const std::vector<PVDRow>& { return R.pvdo; });
    }
    if (vap) {
        o << "PVTG\n";
        for (size_t t = 0; t < PR.size(); ++t) {
            if (t) o << "-- PVT region " << t + 1 << "\n";
            for (auto& r : PR[t].pvtg) {
                o << " " << num(r.pg);
                for (size_t j = 0; j < r.rows.size(); ++j)
                    o << (j ? "\n        " : "  ") << num(r.rows[j].rv) << ' ' << num(r.rows[j].B) << ' ' << num(r.rows[j].mu);
                o << " /\n";
            }
            o << "/\n";
        }
        o << "\n";
    } else if (c.phases.gas) table3("PVDG", [](const PvtRegion& R) -> const std::vector<PVDRow>& { return R.pvdg; });
    auto satTable = [&](const char* kw, bool gas) {
        o << kw << "\n";
        for (auto& S : SR) { for (auto& r : gas ? S.sgof : S.swof) o << " " << num(r.s) << ' ' << num(r.krw) << ' ' << num(r.kro) << ' ' << num(r.pc) << "\n"; o << "/\n"; }
        o << "\n";
    };
    if (c.phases.water && c.phases.oil) satTable("SWOF", false);
    if (c.phases.gas && c.phases.oil) satTable("SGOF", true);
    if (c.phases.water && c.phases.oil && c.phases.gas && c.props.kroModel == 1) {
        o << "STONE1\n\n";
        if (c.props.stone1Exp != 1.0) { o << "STONE1EX\n"; for (size_t r = 0; r < SR.size(); ++r) o << " " << num(c.props.stone1Exp) << " /\n"; o << "\n"; }
    } else if (c.phases.water && c.phases.oil && c.phases.gas && c.props.kroModel == 2) o << "STONE2\n\n";

    if (polyOn) {
        const PolymerInput& P = c.polymer;
        auto ptab = [&](const char* kw, const std::vector<std::vector<double>>& T) {
            if (T.empty()) return;
            o << kw << "\n";
            for (auto& t : T) { for (size_t i = 0; i + 1 < t.size(); i += 2) o << " " << num(t[i]) << ' ' << num(t[i + 1]) << "\n"; o << "/\n"; }
            o << "\n";
        };
        ptab("PLYVISC", P.plyvisc);
        ptab("PLYADS", P.plyads);
        if (!P.plyrock.empty()) {
            o << "PLYROCK\n";
            for (auto& R : P.plyrock) o << " " << num(R.deadPoreVolume) << ' ' << num(R.rrf) << ' ' << num(R.rockDensity) << ' ' << R.adsIndex << ' ' << num(R.maxAdsorption) << " /\n";
            o << "\n";
        }
        o << "PLYMAX\n";
        for (size_t r = 0; r < P.cmax.size(); ++r) o << " " << num(P.cmax[r]) << ' ' << num(r < P.cmaxSalt.size() ? P.cmaxSalt[r] : 0.0) << " /\n";
        o << "\nPLMIXPAR\n";
        for (double m : P.mixing) o << " " << num(m) << " /\n";
        o << "\n";
        if (!P.plyshlog.empty()) {
            o << "PLYSHLOG\n";
            for (auto& S : P.plyshlog) {
                o << " " << num(S.refConc);
                if (S.refSalt >= 0 || S.refTemp >= 0) o << ' ' << (S.refSalt >= 0 ? num(S.refSalt) : std::string("1*")) << ' ' << (S.refTemp >= 0 ? num(S.refTemp) : std::string("1*"));
                o << " /\n";
                for (size_t i = 0; i + 1 < S.rows.size(); i += 2) o << " " << num(S.rows[i]) << ' ' << num(S.rows[i + 1]) << "\n";
                o << "/\n";
            }
            o << "\n";
        }
    }
    if (!g.swatinit.empty() && !c.init.active()) arr("SWATINIT", g.swatinit);
    if (g.hasEndPointScaling()) {
        // cells without a value (< 0) get the table end point of their SATNUM region
        std::vector<SatFunctions> sf(size_t(p.numSatRegions()));
        bool sfOk = true;
        try { for (int r = 0; r < int(sf.size()); ++r) sf[size_t(r)].setup(c, r); } catch (...) { sfOk = false; }
        auto epsArr = [&](const char* name, const std::vector<double>& v, double SatFunctions::EndPoints::*m) {
            if (v.empty()) return;
            std::vector<double> w = v;
            for (size_t q = 0; q < w.size(); ++q) {
                if (w[q] >= 0) continue;
                const int r = q < g.satnum.size() ? std::clamp(g.satnum[q] - 1, 0, int(sf.size()) - 1) : 0;
                w[q] = sfOk ? sf[size_t(r)].tableEndPoints().*m : 0.0;
            }
            arr(name, w);
        };
        using EP = SatFunctions::EndPoints;
        epsArr("SWL", g.swl, &EP::swl); epsArr("SWCR", g.swcr, &EP::swcr); epsArr("SWU", g.swu, &EP::swu);
        epsArr("SGL", g.sgl, &EP::sgl); epsArr("SGCR", g.sgcr, &EP::sgcr); epsArr("SGU", g.sgu, &EP::sgu);
        epsArr("SOWCR", g.sowcr, &EP::sowcr); epsArr("SOGCR", g.sogcr, &EP::sogcr);
    }

    writeExtra(DeckSection::Props);
    // ---------------- REGIONS
    const bool mixOut = polyOn && !g.plmixnum.empty();
    if (!g.satnum.empty() || !g.pvtnum.empty() || !g.eqlnum.empty() || !g.fipnum.empty() || mixOut) {
        o << "REGIONS\n\n";
        arr("SATNUM", g.satnum); arr("PVTNUM", g.pvtnum); arr("EQLNUM", g.eqlnum); arr("FIPNUM", g.fipnum);
        if (mixOut) arr("PLMIXNUM", g.plmixnum);
        writeExtra(DeckSection::Regions);
    } else writeExtra(DeckSection::Regions, "REGIONS");

    // ---------------- SOLUTION
    o << "SOLUTION\n\n";
    if (c.init.active()) {
        // enumerated initial state
        arr("PRESSURE", c.init.pressure);
        if (c.phases.water && c.phases.oil) arr("SWAT", c.init.swat);
        if (c.phases.gas) arr("SGAS", c.init.sgas);
        if (live) arr("RS", c.init.rs);
        if (vap && !c.init.rv.empty()) arr("RV", c.init.rv);
    } else {
        std::vector<const EquilInput*> EQ{&c.equil};
        for (auto& x : c.equilExtra) EQ.push_back(&x);
        // Rs initialisation per region: RSVD table / constant Rs (as a two-row RSVD table) / saturated at the GOC (item 7 = 0)
        // OPM Flow / ECLIPSE need the datum at the GOC when no RSVD table is given, so the ResSim default
        // ("Rs saturated at the GOC pressure") is written as a constant RSVD table computed by the equilibration.
        std::vector<double> satRs(EQ.size(), -1.0);
        std::vector<double> satRv(EQ.size(), -1.0);
        bool needSat = false, needSatRv = false;
        for (auto* e : EQ) if (e->rsvdDepth.empty() && e->initialRs < 0) needSat = true;
        for (auto* e : EQ) if (e->rvvdDepth.empty() && e->initialRv < 0) needSatRv = true;
        if ((live && needSat) || (vap && needSatRv)) {
            try {
                Simulator s; s.log = [](const std::string&) {};
                s.setup(c); s.initialize();
                const auto& r = s.initialRsByRegion();
                for (size_t e = 0; e < EQ.size() && e < r.size(); ++e) satRs[e] = r[e];
                const auto& v = s.initialRvByRegion();
                for (size_t e = 0; e < EQ.size() && e < v.size(); ++e) satRv[e] = v[e];
            } catch (...) {}
        }
        if (live && needSat) {
            for (size_t e = 0; e < EQ.size(); ++e)
                if (satRs[e] < 0 && EQ[e]->rsvdDepth.empty() && EQ[e]->initialRs < 0) {
                    // fallback: saturated Rs at the datum pressure from the PVTO curve of PVT region 1
                    double best = 0, bd = 1e300;
                    for (auto& r : p.pvto) { const double d = std::fabs(r.pb - EQ[e]->datumPressure); if (d < bd) { bd = d; best = r.rs; } }
                    satRs[e] = best;
                }
        }
        if (vap && needSatRv) {
            for (size_t e = 0; e < EQ.size(); ++e)
                if (satRv[e] < 0 && EQ[e]->rvvdDepth.empty() && EQ[e]->initialRv < 0 && !p.pvtg.empty()) {
                    // fallback: saturated Rv at the datum pressure (linear in the PVTG saturated curve)
                    const auto& T = p.pvtg; const double pd = EQ[e]->datumPressure;
                    double rv = T.front().rows.front().rv;
                    for (size_t k = 1; k < T.size(); ++k) {
                        rv = T[k].rows.front().rv;
                        if (pd <= T[k].pg) { const double w = (pd - T[k - 1].pg) / (T[k].pg - T[k - 1].pg);
                                             rv = T[k - 1].rows.front().rv + std::clamp(w, 0.0, 1.0) * (rv - T[k - 1].rows.front().rv); break; }
                    }
                    satRv[e] = rv;
                }
        }
        auto rsOf = [&](size_t e) { return EQ[e]->initialRs >= 0 ? EQ[e]->initialRs : satRs[e]; };
        bool anyTable = live;
        o << "EQUIL\n";
        for (size_t ei = 0; ei < EQ.size(); ++ei) {
            const EquilInput* e = EQ[ei];
            const bool tab = anyTable && (!e->rsvdDepth.empty() || rsOf(ei) >= 0);
            o << ' ' << num(e->datumDepth) << ' ' << num(e->datumPressure) << ' ' << num(e->owc) << ' ' << num(e->pcowOWC)
              << ' ' << num(e->goc) << ' ' << num(e->pcogGOC) << ' ' << (anyTable ? (tab ? "1" : "0") : "1*") << (vap ? " 1 " : " 1* ") << e->nacc << " /\n";
        }
        o << "\n";
        if (anyTable) {
            o << "RSVD\n";
            for (size_t ei = 0; ei < EQ.size(); ++ei) {
                const EquilInput* e = EQ[ei];
                if (!e->rsvdDepth.empty()) {
                    for (size_t i = 0; i < e->rsvdDepth.size() && i < e->rsvdRs.size(); ++i) o << ' ' << num(e->rsvdDepth[i]) << ' ' << num(e->rsvdRs[i]) << "\n";
                } else {
                    const double rs = std::max(0.0, rsOf(ei));
                    o << ' ' << num(std::min(e->goc, e->datumDepth) - 1000) << ' ' << num(rs) << "\n " << num(std::max(e->owc, e->datumDepth) + 1000) << ' ' << num(rs) << "\n";
                }
                o << "/\n";
            }
            o << "\n";
        }
        if (vap) {
            // Rv per region: RVVD table, or the constant / saturated-at-GOC value as a two-row table
            o << "RVVD\n";
            for (size_t ei = 0; ei < EQ.size(); ++ei) {
                const EquilInput* e = EQ[ei];
                if (!e->rvvdDepth.empty()) {
                    for (size_t i = 0; i < e->rvvdDepth.size() && i < e->rvvdRv.size(); ++i) o << ' ' << num(e->rvvdDepth[i]) << ' ' << num(e->rvvdRv[i]) << "\n";
                } else {
                    const double rv = std::max(0.0, e->initialRv >= 0 ? e->initialRv : satRv[ei]);
                    o << ' ' << num(std::min(e->goc, e->datumDepth) - 1000) << ' ' << num(rv) << "\n " << num(std::max(e->owc, e->datumDepth) + 1000) << ' ' << num(rv) << "\n";
                }
                o << "/\n";
            }
            o << "\n";
        }
    }
    if (!c.thpres.empty()) {
        o << "THPRES\n";
        for (auto& t : c.thpres) o << ' ' << t.r1 << ' ' << t.r2 << ' ' << (t.value >= 0 ? num(t.value) : std::string("1*")) << " /\n";
        o << "/\n\n";
    }
    bool anyFet = false, anyCT = false;
    for (auto& a : c.aquifers) (a.type == AquiferType::Fetkovich ? anyFet : anyCT) = true;
    auto p0s = [&](double v) { return v > 0 ? num(v) : std::string("1*"); };
    if (anyFet) {
        o << "AQUFETP\n";
        for (size_t i = 0; i < c.aquifers.size(); ++i) {
            const auto& a = c.aquifers[i];
            if (a.type != AquiferType::Fetkovich) continue;
            o << ' ' << a.id << ' ' << (aqDepth[i] > 0 ? num(aqDepth[i]) : std::string("1*")) << ' ' << p0s(a.p0) << ' ' << num(a.volume)
              << ' ' << num(a.ct) << ' ' << num(a.pi) << " 1 /\n";
        }
        o << "/\n\n";
    }
    if (anyCT) {
        o << "AQUCT\n";
        for (size_t i = 0; i < c.aquifers.size(); ++i) {
            const auto& a = c.aquifers[i];
            if (a.type != AquiferType::CarterTracy) continue;
            o << ' ' << a.id << ' ' << (aqDepth[i] > 0 ? num(aqDepth[i]) : std::string("1*")) << ' ' << p0s(a.p0) << ' ' << num(a.perm)
              << ' ' << num(a.poro) << ' ' << num(a.ct) << ' ' << num(a.rInner) << ' ' << num(a.thickness) << ' ' << num(a.angle) << " 1 1 /\n";
        }
        o << "/\n\n";
    }
    for (auto& a : c.aquifers) if (a.visc > 0) o << "-- RESSIM: AQVISC " << a.id << ' ' << num(a.visc) << "\n";

    if (polyOn && !g.spoly.empty()) arr("SPOLY", g.spoly);
    writeExtra(DeckSection::Solution);
    // ---------------- SUMMARY
    const std::set<std::string>& ownSummary = exporterSummary();
    o << "SUMMARY\n\nFOPR\nFWPR\nFGPR\nFLPR\nFWIR\nFGIR\nFOPT\nFWPT\nFGPT\nFWIT\nFGIT\nFPR\nFPRP\nFWCT\nFGOR\nFOIP\nFWIP\nFGIP\n";
    if (!c.aquifers.empty()) o << "FAQR\nFAQT\n";
    if (polyOn) o << "FCIR\nFCPR\nFCIT\nFCPT\nFCIP\nFCAD\nWCIR\n/\nWCPR\n/\nWCIT\n/\nWCPT\n/\n";
    o << "WBHP\n/\nWOPR\n/\nWWPR\n/\nWGPR\n/\nWLPR\n/\nWWIR\n/\nWGIR\n/\nWWCT\n/\nWGOR\n/\n"
         "WOPT\n/\nWWPT\n/\nWGPT\n/\nWWIT\n/\nWGIT\n/\n\n";

    writeExtra(DeckSection::Summary, nullptr, &ownSummary);
    // ---------------- SCHEDULE
    o << "SCHEDULE\n\n";
    const Tuning& t = c.tuning;
    o << "TUNING\n " << num(t.initialDt) << ' ' << num(t.maxDt) << ' ' << num(t.minDt) << " 1* " << num(t.growth) << " 1* " << num(t.cutFactor) << " /\n"
      << " 1* " << num(t.cnvTol) << ' ' << num(t.mbTol) << ' ' << num(t.linTol) << " /\n"
      << ' ' << t.maxNewton << " 1* " << t.maxLinear << " 1* 1* 1* 1* " << num(t.dsMax) << " /\n";
    o << "-- RESSIM: TUNING " << num(t.dpMaxRel) << ' ' << t.threads << ' ' << t.linearSolver << ' ' << t.cprWeights << "\n\n";
    o << "RPTRST\n 'BASIC=2' /\n\n";
    // COMPDAT for the connections that open (or shut) at report step s
    auto compRecord = [&](const WellInput& w, const Completion& q, bool open) {
        std::ostringstream r;
        r << " '" << w.name << "' " << q.i + 1 << ' ' << q.j + 1 << ' ' << q.k1 + 1 << ' ' << q.k2 + 1 << ' ' << (open ? "OPEN" : "SHUT") << ' ' << (q.satTable > 0 ? std::to_string(q.satTable) : std::string("1*")) << ' ';
        if (q.cf > 0) r << num(q.cf); else r << "1*";
        r << ' ' << num(q.diameter) << ' ';
        if (q.kh > 0) r << num(q.kh); else r << "1*";
        r << ' ' << num(q.skin) << ' ' << (q.dfactor >= 0 ? num(q.dfactor) : std::string("1*")) << ' ' << char(std::toupper(q.dir));
        if (q.r0 > 0) r << ' ' << num(q.r0);
        r << " /\n";
        return r.str();
    };
    auto compEvents = [&](size_t s) {
        std::string out;
        for (auto& w : c.wells)
            for (auto& q : w.completions) {
                if (s == 0 ? q.fromStep <= 0 : q.fromStep == int(s)) out += compRecord(w, q, q.open);
                if (s > 0 && q.open && q.untilStep == int(s)) out += compRecord(w, q, false);
            }
        return out;
    };
    auto writeCompdat = [&](size_t s) { const std::string e = compEvents(s); if (!e.empty()) o << "COMPDAT\n" << e << "/\n\n"; };
    if (!c.wells.empty()) {
        o << "WELSPECS\n";
        for (auto& w : c.wells) {
            const char* ph = w.type == WellType::Producer ? "OIL" : w.type == WellType::WaterInjector ? "WATER" : "GAS";
            o << " '" << w.name << "' '" << (w.group.empty() ? std::string("G") : w.group) << "' " << w.i + 1 << ' ' << w.j + 1 << ' ';
            if (w.refDepth > 0) o << num(w.refDepth); else o << "1*";
            o << " " << ph;
            if (!w.specsTail.empty()) o << ' ' << w.specsTail;
            o << " /\n";
        }
        o << "/\n\n";
        writeCompdat(0);
        bool anyPI = false; for (auto& w : c.wells) if (w.piMult != 1.0) anyPI = true;
        if (anyPI) {
            o << "WPIMULT\n";
            for (auto& w : c.wells) if (w.piMult != 1.0) o << " '" << w.name << "' " << num(w.piMult) << " /\n";
            o << "/\n\n";
        }
    }
    if (!c.observations.empty()) {
        o << "-- Observed data for history matching (ResSim extension)\n";
        for (auto& ob : c.observations)
            o << "-- RESSIM: OBS " << toString(ob.qty) << ' ' << (ob.well.empty() ? std::string("-") : ob.well) << ' ' << num(ob.time) << ' '
              << num(ob.value) << ' ' << num(ob.sigma) << (ob.weight != 1 ? " " + num(ob.weight) : std::string()) << "\n";
        o << "\n";
    }
    if (!c.studyParams.empty()) {
        o << "-- Sensitivity / history-matching parameters (ResSim extension)\n";
        for (auto& p : c.studyParams)
            o << "-- RESSIM: PARAM " << int(p.kind) << ' ' << p.target << ' ' << p.i1 << ' ' << p.i2 << ' ' << p.j1 << ' ' << p.j2 << ' ' << p.k1 << ' '
              << p.k2 << ' ' << num(p.base) << ' ' << num(p.low) << ' ' << num(p.high) << ' ' << int(p.log) << ' ' << int(p.active) << " '" << p.label
              << "'\n";
        o << "\n";
    }
    auto findWell = [&](const std::string& nm) -> const WellInput* { for (auto& w : c.wells) if (w.name == nm) return &w; return nullptr; };
    auto dflt = [&](double v) { return v < 0 ? std::string("1*") : num(v); };
    auto statusTag = [](WellStatus s) { return s == WellStatus::Open ? "OPEN" : s == WellStatus::Shut ? "SHUT" : "STOP"; };
    // open connections of a well during report step s
    auto openConns = [&](const WellInput& w, int s) {
        int n = 0;
        for (auto& q : w.completions) if (q.open && q.fromStep <= s && (q.untilStep < 0 || q.untilStep > s)) ++n;
        return n;
    };
    std::map<std::string, WellControl> inForce;      // control record in force per well
    size_t s = 0;
    while (s < c.schedule.size()) {
        std::ostringstream prod, inj, hist;
        if (s > 0) writeCompdat(s);
        // ResSim opens a well with an OPEN control as soon as its first connection opens; ECLIPSE / OPM need the
        // control record again at that step (a control given to a well without connections leaves it shut)
        std::vector<WellControl> ctrls = c.schedule[s].controls;
        for (auto& w : c.wells) {
            if (s == 0 || openConns(w, int(s)) == 0 || openConns(w, int(s) - 1) > 0) continue;
            auto it = inForce.find(w.name);
            if (it == inForce.end() || it->second.status != WellStatus::Open) continue;
            bool given = false;
            for (auto& wc : ctrls) if (matchName(wc.well, w.name)) given = true;
            if (!given) ctrls.push_back(it->second);
        }
        for (auto& wc : ctrls) for (auto& w : c.wells) if (matchName(wc.well, w.name)) { WellControl x = wc; x.well = w.name; inForce[w.name] = x; }
        for (auto& wc : ctrls) {
            const WellInput* w = findWell(wc.well);
            bool isProd = w ? w->type == WellType::Producer : wc.mode != Control::RATE;
            if (!w) for (auto& x : c.wells) if (matchName(wc.well, x.name)) { isProd = x.type == WellType::Producer; break; }
            const char* st = statusTag(wc.status);
            std::string mode = wc.deckMode.empty() ? std::string(toString(wc.mode)) : wc.deckMode;
            // a defaulted control (typically 'W' SHUT 15*) was read as NONE: NONE is not a valid WCONPROD/WCONHIST mode in
            // ECLIPSE and OPM rejects it ("no value has been specified for this control") - write it defaulted again
            if (mode == "NONE") mode = "1*";
            const std::string vfp = wc.vfp > 0 ? std::to_string(wc.vfp) : std::string("1*");
            if (wc.history && isProd)
                hist << " '" << wc.well << "' " << st << ' ' << mode << ' ' << dflt(wc.orat) << ' ' << dflt(wc.wrat) << ' '
                     << dflt(wc.grat) << ' ' << vfp << ' ' << dflt(wc.alq) << ' ' << dflt(wc.thp) << ' ' << dflt(wc.bhp) << " /\n";
            else if (isProd)
                prod << " '" << wc.well << "' " << st << ' ' << mode << ' ' << dflt(wc.orat) << ' ' << dflt(wc.wrat) << ' '
                     << dflt(wc.grat) << ' ' << dflt(wc.lrat) << ' ' << dflt(wc.resv) << ' ' << dflt(wc.bhp) << ' ' << dflt(wc.thp) << ' ' << vfp << ' '
                     << dflt(wc.alq) << " /\n";
            else {
                bool gas = wc.phase ? wc.phase == 2 : (w && w->type == WellType::GasInjector);
                if (!w && !wc.phase) for (auto& x : c.wells) if (matchName(wc.well, x.name)) { gas = x.type == WellType::GasInjector; break; }
                const bool resvStar = wc.deckMode == "RESV*";     // RESV injector without a surface-rate limit
                const std::string im = resvStar ? std::string("RESV") : mode;
                inj << " '" << wc.well << "' " << (gas ? "GAS" : "WATER") << ' ' << st << ' '
                    << im << ' ' << (resvStar ? std::string("1*") : dflt(wc.rate)) << ' ' << dflt(wc.resv) << ' ' << dflt(wc.bhp) << ' '
                    << dflt(wc.thp) << ' ' << vfp << " /\n";
            }
        }
        o << extraAtStep(int(s), false, 1);
        auto trimDefaults = [](std::string s) {          // " 1* 1* /" at the end of a record == " /"
            for (size_t p; (p = s.find(" 1* /\n")) != std::string::npos;) s.erase(p, 3);
            return s;
        };
        if (!prod.str().empty()) o << "WCONPROD\n" << trimDefaults(prod.str()) << "/\n\n";
        if (!hist.str().empty()) o << "WCONHIST\n" << hist.str() << "/\n\n";
        if (!inj.str().empty()) o << "WCONINJE\n" << trimDefaults(inj.str()) << "/\n\n";
        if (c.schedule[s].drsdt > -1.5 && c.phases.disgas) {
            o << "DRSDT\n " << (c.schedule[s].drsdt < 0 ? std::string("1e20") : num(c.schedule[s].drsdt)) << (c.schedule[s].drsdtFree ? " 'FREE'" : "") << " /\n\n";
        }
        if (polyOn && !c.schedule[s].wpolymer.empty()) {
            o << "WPOLYMER\n";
            for (auto& wp : c.schedule[s].wpolymer) o << " '" << wp.well << "' " << num(wp.conc) << " 0 /\n";
            o << "/\n\n";
        }
        o << extraAtStep(int(s), false, 0);
        std::vector<double> lens{c.schedule[s].length};
        size_t e2 = s + 1;
        while (e2 < c.schedule.size() && c.schedule[e2].controls.empty() && c.schedule[e2].wpolymer.empty() && c.schedule[e2].drsdt < -1.5 && compEvents(e2).empty() && extraAtStep(int(e2), false).empty())
            lens.push_back(c.schedule[e2++].length);
        o << "TSTEP\n"; writeArray(o, lens); o << "/\n\n";
        s = e2;
    }
    o << extraAtStep(int(c.schedule.size()), true);    // after the last report step
    o << "END\n";
    return wrapLongLines(o.str());
}

void exportEclipseDeck(const Case& c, const std::string& path) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot write " + path);
    f << eclipseDeckToString(c);
}

}  // namespace ressim
