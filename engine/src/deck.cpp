// ECLIPSE 100 deck import and export (two-way).
// Import = generic deck reader (deck_parser.cpp, OPM keyword database) + the semantic mapping below.
#include "ressim/eos.h"
#include "ressim/io.h"
#include "ressim/udq.h"
#include "ressim/therm_deck.h"
#include "ressim/grid.h"
#include "ressim/simulator.h"
#include "ressim/deck_parser.h"
#include "ressim/examples.h"
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
#include <filesystem>

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
            if (lead == 0 && (line.compare(0, 10, "-- RESSIM:") == 0 || line.compare(0, 10, "-- RESSIM+") == 0)) {
                // ResSim extension lines carry data: continue on '-- RESSIM+' lines (joined again by the parser)
                size_t cut = line.rfind(' ', maxLen);
                if (cut == std::string::npos || cut <= 11) { line.resize(maxLen); break; }
                out.append(line, 0, cut);
                out += cr ? "\r\n" : "\n";
                size_t rest = line.find_first_not_of(' ', cut);
                line = rest == std::string::npos ? std::string() : "-- RESSIM+ " + line.substr(rest);
                if (line.empty()) break;
                continue;
            }
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

// ECLIPSE / OPM name patterns: '*' = any run of characters, '?' = one character (shell-style, anywhere in the pattern)
bool matchName(const std::string& pat, const std::string& name) {
    if (pat.find_first_of("*?") == std::string::npos) return pat == name;
    size_t p = 0, n = 0, star = std::string::npos, mark = 0;
    while (n < name.size()) {
        if (p < pat.size() && (pat[p] == '?' || pat[p] == name[n])) { ++p; ++n; }
        else if (p < pat.size() && pat[p] == '*') { star = p++; mark = n; }
        else if (star != std::string::npos) { p = star + 1; n = ++mark; }
        else return false;
    }
    while (p < pat.size() && pat[p] == '*') ++p;
    return p == pat.size();
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
            "SWATINIT THPRES VAPOIL PVTG RVVD RV WELSPECS COMPDAT WCONPROD WCONHIST WCONINJE WCONINJH WELOPEN WELOPENL WPIMULT WELTARG TSTEP DATES TUNING ROCKOPTS NEWTRAN OLDTRAN "
            "ENDSCALE SWL SWCR SWU SGL SGCR SGU SOWCR SOGCR "
            "POLYMER PLYVISC PLYADS PLYROCK PLYMAX PLMIXPAR WPOLYMER SPOLY PLYSHLOG PLMIXNUM STONE STONE1 STONE2 STONE1EX DRSDT SATOPTS EHYSTR SCALECRS VAPPARS IMBNUM ISWL ISWCR ISWU ISGL ISGCR ISGU ISOWCR ISOGCR NOHYST "
            "DUALPORO DUALPERM SIGMA SIGMAV LTOSIGMA LX LY LZ MULTMF NODPPM DPGRID PINCH GRAVDR DZMTRX DZMTRXV DPNUM");
        add(Support::Partial, "C35 SOLVENT: simulated by the ResSim black-oil engine (four-component solvent model, immiscible SSFN relative permeabilities; the miscible MISC / Todd-Longstaff data are kept for OPM Flow)", "SOLVENT PVDS SDENSITY SSFN WSOLVENT SSOL MISC PMISC MSFN SORWMIS SGCWMIS TLMIXPAR TLPMIXPA");
        add(Support::Partial, "SHRATE (shear rate based polymer shear model) is not modelled - the velocity form of PLYSHLOG is used", "SHRATE");
        add(Support::Partial, "C35 chemical EOR: simulated by the ResSim black-oil engine (sequential-implicit transport coupled to the flow; surfactant viscosity / IFT / capillary-number miscibility / adsorption, alkaline IFT and adsorption multipliers, foam gas mobility reduction with pressure / shear / decay)", "SURFACT SURFACTW ALKALINE FOAM SURFVISC SURFST SURFCAPD SURFADS SURFROCK SURFNUM WSURFACT ALKADS ALKROCK ALSURFST ALSURFAD ALPOLADS WALKALIN FOAMOPTS FOAMMOB FOAMMOBP FOAMMOBS FOAMADS FOAMROCK FOAMDCYW FOAMDCYO WFOAM SFOAM");
        add(Support::Partial, "passive tracer: simulated by the ResSim black-oil engine (C27/C28, upwind implicit transport; decay, adsorption, oil-water partitioning, decay chains, dispersion via RESSIM TRACERPROP)", "TRACER WTRACER");
        add(Support::Partial, "SCHEDULE-section MULTREGT is not applied (GRID/EDIT records are fully modelled)", "MULTREGT");
        add(Support::Partial, "C31: evaluated by the ResSim engine - UDQ DEFINE / ASSIGN / UPDATE / UNITS with field / well / group quantities, UDA targets in WCONPROD / WCONINJE / WELTARG; ACTIONX conditions (AND / OR / parentheses, well sets, '?') at the end of each report step - new connections from an action are not added at run time", "UDQ ACTIONX");
        add(Support::Partial, "group control: targets allocated by well potentials; priority / sales options are not modelled", "GRUPTREE GCONPROD GCONINJE GCONPRI VFPPROD VFPINJ");
        add(Support::Partial, "surface network (extended BRANPROP / NODEPROP or standard GRUPNET): balanced by the ResSim engine every time step (C29 / C35: node pressures from the branch VFP tables with branch ALQ, no-loss 9999 branches, NEFAC, lift gas added to node flows, GCONSUMP gas taken at nodes, satellite flows -> producer THP limits, automatic chokes; C36: NETBALAN interval / node tolerance / max iterations / max balance error, NETCOMPA automatic compressors switched on by group targets (VFP / ALQ / levels / fuel gas), COMPOFF)", "BRANPROP NODEPROP NETBALAN GRUPNET NEFAC NETCOMPA COMPOFF");
        add(Support::Full, "C35 facilities: satellite group production (added to group / field rates, fixed in group control) and gas consumption / import (FGCR / FGIMR / FGSR)", "GSATPROD GCONSUMP GCONSALE");
        add(Support::Full, "C37: grid length unit of COORD / ZCORN (METRES, FEET, CM) - converted to the deck unit system on import", "GRIDUNIT");
        add(Support::Ignored, "output / dimension keyword - no effect on the results",
            "UNIFOUT UNIFIN FMTOUT FMTIN MULTOUT MULTIN INIT NOECHO ECHO NOWARN WARN MESSAGES MONITOR NOMONITO NOSIM "
            "GRIDFILE MAPAXES MAPUNITS NOINSPEC NORSSPEC RPTGRID RPTPROPS RPTREGS RPTSOL RPTSCHED "
            "RPTRST RPTRUNSP RPTSMRY RPTONLY RPTONLYO OUTSOL WELLDIMS AQUDIMS NSTACK FAULTDIM NUPCOL VFPPDIMS VFPIDIMS "
            "SMRYDIMS WSEGDIMS UDQDIMS UDADIMS UDTDIMS ACTDIMS NETWORK COORDSYS OPTIONS OPTIONS3 SAVE EXTRAPMS DEBUG DEBUG3 "
            "GRIDOPTS EQLOPTS LIVEOIL FULLIMP IMPES CPR ZIPPY2 FORMFEED ALL MSUMLINS MSUMNEWT SEPARATE "
            "EXCEL RUNSUM NARROW PERFORMA DATE NEXTSTEP NEXT RPTPRINT PRINTREP DIMPES DIFFDP UNCODHMD MEMORY NOINSPEC "
            "INSPEC RSSPEC ECHO PETOPTS PARALLEL WRFT WRFTPLT COMPORD TRACERS GDORIENT NOGGF UDQPARAM "
            "EXIT SMRYDIMS");
        return m;
    }();
    return t;
}

// E300 compositional keywords used by the ResSim compositional core (kept verbatim for export / OPM / ECLIPSE)
static const std::set<std::string>& compositionalKeywords() {
    static const std::set<std::string> s = {"COMPS", "NCOMPS", "EOS", "PRCORR", "CNAMES", "TCRIT", "PCRIT", "VCRIT", "ZCRIT", "ACF", "MW", "BIC",
        "SSHIFT", "PARACHOR", "VCRITVIS", "ZCRITVIS", "LBCCOEF", "OMEGAA", "OMEGAB", "RTEMP", "RTEMPA", "TEMPVD", "ZMFVD", "ZI", "STCOND",
        "WELLSTRE", "WINJGAS", "FIELDSEP", "NCOMPS", "SALINITY", "RTEMPVD",
        "DIFFUSE", "DIFFCWAT", "DIFFCGAS", "DIFFC", "PRECSALT", "VAPWAT", "SALTSOL", "PERMFACT", "SEPCOND", "WSEPCOND", "CO2SOL"};
    return s;
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
    // local grid refinement context (CARFIN / REFINE ... ENDFIN): grid arrays and BOX act on the local cells
    int lgrCtx = -1;
    double dzmtrx = 0;                 // DZMTRX (scalar matrix block height)
    std::vector<double> dpnumVals;     // DPNUM as read
    struct SavedGlobal { std::map<std::string, std::vector<double>> arr; int n = 0, nl = 0, nx = 0, ny = 0, nz = 0; Box box{0, 0, 0, 0, 0, 0}; bool boxActive = false; } savedG;
    std::vector<std::map<std::string, std::vector<double>>> lgrArr;
    std::map<int, std::vector<AquiferConnection>> aqConn;
    std::map<int, int> aqTabOf;                          // C34 AQUCT influence table number per aquifer
    bool satDirect = false;                              // C34 SATOPTS DIRECT
    double gridUnitFac = 1.0;                            // C37 GRIDUNIT: COORD / ZCORN length unit -> deck length unit
    std::vector<std::vector<double>> aqutab;             // C34 AQUTAB tables 2, 3, ... (tD, pD pairs)
    std::map<int, double> aqVisc;
    std::vector<WellControl> pending;
    std::vector<WellPolymer> pendingPoly;
    double pendingDrsdt = -2; bool pendingDrsdtFree = false;
    double pendingVap1 = -1, pendingVap2 = -1;
    bool pendingVapAfterDrsdt = false;
    std::vector<GroupProdControl> pendingGprod;
    GuideRateSpec pendingGuide;
    std::vector<std::pair<std::string, double>> pendingWelpi;   // C34
    std::map<std::string, double> ressimPiMult;                   // C44 "-- RESSIM: PIMULT well f": whole-run PI multiplier baked into WPIMULT
    std::vector<GroupInjControl> pendingGinj;
    bool ehystr = false, noHyst = false;
    std::map<std::string, WellControl> lastCtrl;
    std::vector<std::pair<std::string, std::array<int, 5>>> complump;   // C46 COMPLUMP records (well pattern; I J K1 K2 N)
    std::map<std::string, bool> lastInj, pendingInj;     // control came from an injector keyword
    std::set<std::pair<size_t, std::string>> injCtrl;    // (report step, well) of injector controls
    long dayNow = 0;
    double daysElapsed = 0;
    bool anyStep = false;
    // tables (per region)
    std::vector<PvtRegion> pvt;
    std::vector<char> hasPvtw, hasRock, hasDensity, hasPvdo, hasPvdg, hasPvto;
    std::vector<std::vector<double>> swof, sgof, swfn, sgfn, sof3, sof2, slgof, sgwfn;   // [table] -> flat rows
    std::vector<EquilInput> equil;
    std::vector<std::vector<double>> rsvd, pbvd;
    std::vector<int> equilItem7, equilItem8;
    std::vector<std::vector<double>> rvvd;
    bool sawEquil = false;
    struct MultRegT { int r1, r2; double m; std::string dirs, set, nncMode; };
    std::vector<MultRegT> multregt;
    // grouped messages
    struct Group { int level; std::string keyword, file; int line; std::vector<std::string> items; int count = 0; };
    std::map<std::string, Group> groups;
    std::vector<std::string> groupOrder;

    bool thermalDeck = false;   // THERMAL present: thermalFromDeck builds the model after the black-oil import
    explicit Importer(const deck::RawDeck& d) : D(d) {
        for (const auto& k : d.keywords) if (k.name == "THERMAL") thermalDeck = true;
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
    // a feature the ResSim engine cannot model: the import continues (placeholders) and the case is only run by
    // OPM Flow / ECLIPSE on the original deck
    bool co2storeDeck() const { for (const auto& k : D.keywords) if (k.name == "CO2STORE") return true; return false; }
    void unsup(const std::string& reason) {
        for (auto& r : c.engineUnsupported) if (r == reason) return;
        c.engineUnsupported.push_back(reason);
        if (thermalDeck && (reason.rfind("compositional model", 0) == 0 || reason.rfind("thermal model", 0) == 0)) return;
        if (reason.rfind("compositional model", 0) == 0) {
            msg(0, "compositional deck (COMPS): runs on the ResSim EOS compositional core (engine ResSim) or with OPM Flow / ECLIPSE E300 (original deck)");
            return;
        }
        if (reason.rfind("CO2 storage", 0) == 0) {   // C47: the CO2STORE core runs it (compositionalFromDeck removes the marker)
            msg(0, "CO2 storage deck (CO2STORE): runs on the ResSim CO2STORE core (engine ResSim) or with OPM Flow (original deck)");
            return;
        }
        if (reason.rfind("gas-water system", 0) == 0 && co2storeDeck()) return;
        msg(1, "ResSim engine: " + reason + " - the deck is imported for viewing / editing and runs with OPM Flow or ECLIPSE (original deck)");
    }
    [[noreturn]] void fail(const std::string& text, int line = 0) const {
        throw std::runtime_error((K ? K->name : std::string("Deck")) + " (" + where(line) + "): " + text);
    }

    // ------------------------------------------------------------------ item access (UDA-tolerant)
    // a numeric item that may hold a UDQ name (UDA, C31): the name is recorded in wc.uda and 'def' returned
    double unum(const Record& r, size_t i, double def, WellControl& wc, const char* target) {
        if (r.has(i) && !r.isNumber(i) && !r.items[i].def) {
            const std::string u = upper(r.items[i].s);
            if (u.size() > 2 && (u[0] == 'F' || u[0] == 'W' || u[0] == 'G') && u[1] == 'U') {
                for (auto& x : wc.uda) if (x.first == target) { x.second = u; return def; }
                wc.uda.push_back({target, u});
                grouped(0, "UDA: the target follows a user-defined quantity (UDQ) during the run", wc.well + " " + target + "=" + u);
                return def;
            }
        }
        return num(r, i, def);
    }
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
        if (v.empty()) {
            if (lgrCtx >= 0) v = inheritFromHost(name);
            else v.assign(size_t(n), arrayDefault(name));
        }
        return v;
    }
    // ---- local grid refinements
    // a local array first touched inside CARFIN / REFINE starts from the host cells' current values (NaN = not given)
    std::vector<double> inheritFromHost(const std::string& name) {
        std::vector<double> v(size_t(n), NaN);
        const LgrInput& L = c.grid.lgrs[size_t(lgrCtx)];
        auto it = savedG.arr.find(name);
        if (it == savedG.arr.end() || it->second.size() != size_t(savedG.n)) return v;
        std::vector<int> hi, hj, hk; std::vector<double> a, b;
        L.split(0, hi, a, b); L.split(1, hj, a, b); L.split(2, hk, a, b);
        for (int k = 0; k < L.nz; ++k)
            for (int j = 0; j < L.ny; ++j)
                for (int i = 0; i < L.nx; ++i)
                    v[size_t(L.index(i, j, k))] = it->second[size_t(hi[size_t(i)] + savedG.nx * (hj[size_t(j)] + savedG.ny * hk[size_t(k)]))];
        return v;
    }
    void enterLgr(int l) {
        if (lgrCtx >= 0) leaveLgr();
        const LgrInput& L = c.grid.lgrs[size_t(l)];
        savedG.arr = std::move(arr); savedG.n = n; savedG.nl = nl; savedG.nx = c.grid.nx; savedG.ny = c.grid.ny; savedG.nz = c.grid.nz;
        savedG.box = box; savedG.boxActive = boxActive;
        arr = std::move(lgrArr[size_t(l)]);
        n = L.size(); nl = L.nx * L.ny; c.grid.nx = L.nx; c.grid.ny = L.ny; c.grid.nz = L.nz; boxActive = false;
        lgrCtx = l;
    }
    void leaveLgr() {
        if (lgrCtx < 0) return;
        lgrArr[size_t(lgrCtx)] = std::move(arr);
        arr = std::move(savedG.arr); n = savedG.n; nl = savedG.nl; c.grid.nx = savedG.nx; c.grid.ny = savedG.ny; c.grid.nz = savedG.nz;
        box = savedG.box; boxActive = savedG.boxActive;
        lgrCtx = -1;
    }
    void carfin(const Keyword& k) {
        needDims();
        leaveLgr();
        if (k.records.empty()) fail("record missing");
        const Record& r = k.records[0];
        LgrInput L;
        L.name = str(r, 0);
        if (L.name.empty()) fail("LGR name missing");
        L.i1 = inum(r, 1, 1) - 1; L.i2 = inum(r, 2, 1) - 1; L.j1 = inum(r, 3, 1) - 1; L.j2 = inum(r, 4, 1) - 1;
        L.k1 = inum(r, 5, 1) - 1; L.k2 = inum(r, 6, 1) - 1;
        L.nx = inum(r, 7, L.i2 - L.i1 + 1); L.ny = inum(r, 8, L.j2 - L.j1 + 1); L.nz = inum(r, 9, L.k2 - L.k1 + 1);
        L.nwmax = inum(r, 10, 0);
        L.parent = upper(str(r, 11, "GLOBAL"));
        if (L.parent.empty()) L.parent = "GLOBAL";
        if (L.i1 < 0 || L.j1 < 0 || L.k1 < 0 || L.i2 >= c.grid.nx || L.j2 >= c.grid.ny || L.k2 >= c.grid.nz || L.i1 > L.i2 || L.j1 > L.j2 || L.k1 > L.k2)
            fail("host box of LGR " + L.name + " is outside the grid");
        if (c.grid.lgrIndex(L.name) >= 0) fail("duplicate LGR name " + L.name);
        if (L.parent != "GLOBAL") unsup("nested local grid refinement (CARFIN " + L.name + " in " + L.parent + ")");
        c.grid.lgrs.push_back(L);
        lgrArr.emplace_back();
        enterLgr(int(c.grid.lgrs.size()) - 1);
    }
    void finArray(const Keyword& k) {
        LgrInput& L = c.grid.lgrs[size_t(lgrCtx)];
        const char d = k.name[1];
        if (k.name[0] == 'N') {
            std::vector<int>& dst = d == 'X' ? L.nxfin : d == 'Y' ? L.nyfin : L.nzfin;
            dst.clear(); for (double v : k.values) dst.push_back(std::isnan(v) ? 1 : int(std::lround(v)));
        } else {
            std::vector<double>& dst = d == 'X' ? L.hxfin : d == 'Y' ? L.hyfin : L.hzfin;
            dst.clear(); for (double v : k.values) dst.push_back(std::isnan(v) ? 1.0 : v);
        }
        std::string why;
        if (!L.valid(&why) && why.find("FIN") != std::string::npos && why.find(k.name.substr(0, 1) + d) != std::string::npos) warn(why);
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
        if (lgrCtx >= 0 && (name == "DX" || name == "DY" || name == "DZ" || name == "TOPS" || name == "DXV" || name == "DYV" || name == "DZV" ||
                            name == "COORD" || name == "ZCORN")) {
            grouped(1, "geometry inside a local grid refinement is not modelled (local cells subdivide their host cell) - ignored", name);
            return;
        }
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
            c.grid.coord = v;
            if (gridUnitFac != 1.0) for (double& x : c.grid.coord) x *= gridUnitFac;   // C37 GRIDUNIT
            return;
        }
        if (name == "ZCORN") {
            if (v.size() != size_t(8) * size_t(n)) fail("ZCORN needs 8*NX*NY*NZ = " + std::to_string(size_t(8) * n) + " values, found " + std::to_string(v.size()));
            c.grid.zcorn = v;
            if (gridUnitFac != 1.0) for (double& x : c.grid.zcorn) x *= gridUnitFac;   // C37 GRIDUNIT
            return;
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
                if (!arr.count(src) && !(lgrCtx >= 0 && savedG.arr.count(src)) && std::isnan(arrayDefault(src))) { warn("COPY source " + src + " is not defined - ignored", r.line); continue; }
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
                // C38: kept apart from MULTX - the TRANX edit does not touch fault NNCs across the face
                name = std::string("MULTTR") + name[4];
            }
            if (op != "EQUALS" && !arr.count(name) && !(lgrCtx >= 0 && savedG.arr.count(name)) && std::isnan(arrayDefault(name))) { grouped(1, op + " on an array before it is defined - ignored", name); continue; }
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

    // C34 OPERATE (box) / OPERATER (OPERNUM region): target = f(source, alpha, beta) - folded into the arrays like EQUALS
    // (ECLIPSE / OPM Operate: MULTA a*X+b, POLY Y+a*X^b, SLOG 10^(a+b*X), LOG10, LOGE, INV, MULTX a*X, ADDX X+a, COPY,
    // MAXLIM min(a, X), MINLIM max(a, X), MULTP a*X^b, ABS, MULTIPLY Y*X)
    void operate(const Keyword& k) {
        needDims();
        const bool region = k.name == "OPERATER";
        Box prev = curBox();
        for (const Record& r : k.records) {
            if (r.items.empty()) continue;
            const std::string dst = upper(str(r, 0));
            const size_t o = region ? 2 : 7;
            const std::string op = upper(str(r, o)), src = upper(str(r, o + 1, dst));
            const double alpha = num(r, o + 2, 0.0), beta = num(r, o + 3, 0.0);
            std::vector<char> sel(size_t(n), 0);
            if (region) {
                const int reg = inum(r, 1, -1);
                if (reg < 0) { warn("OPERATER region number missing - record ignored", r.line); continue; }
                if (!arr.count("OPERNUM")) { warn("OPERATER needs OPERNUM - record ignored", r.line); continue; }
                const std::vector<double>& R = array("OPERNUM");
                for (int q = 0; q < n; ++q) sel[size_t(q)] = int(std::lround(R[size_t(q)])) == reg;
            } else {
                const Box b = readBox(r, 1, prev); prev = b;
                forBox(b, [&](int q) { sel[size_t(q)] = 1; });
            }
            if (!arr.count(src) && std::isnan(arrayDefault(src))) { warn(k.name + " source " + src + " is not defined - ignored", r.line); continue; }
            static const std::set<std::string> ops = {"MULTA", "POLY", "SLOG", "LOG10", "LOGE", "INV", "MULTX", "ADDX", "COPY", "MAXLIM", "MINLIM", "MULTP", "ABS", "MULTIPLY"};
            if (!ops.count(op)) { warn(k.name + ": operation " + op + " is not known - record ignored", r.line); continue; }
            const std::vector<double> s = array(src);
            auto& a = array(dst);
            for (int q = 0; q < n; ++q) {
                if (!sel[size_t(q)]) continue;
                const double X = s[size_t(q)];
                double& Y = a[size_t(q)];
                if (std::isnan(X)) continue;
                if (op == "MULTA") Y = alpha * X + beta;
                else if (op == "POLY") Y = Y + alpha * std::pow(X, beta);
                else if (op == "SLOG") Y = std::pow(10.0, alpha + beta * X);
                else if (op == "LOG10") Y = X > 0 ? std::log10(X) : Y;
                else if (op == "LOGE") Y = X > 0 ? std::log(X) : Y;
                else if (op == "INV") Y = X != 0 ? 1.0 / X : Y;
                else if (op == "MULTX") Y = alpha * X;
                else if (op == "ADDX") Y = X + alpha;
                else if (op == "COPY") Y = X;
                else if (op == "MAXLIM") Y = std::min(alpha, X);
                else if (op == "MINLIM") Y = std::max(alpha, X);
                else if (op == "MULTP") Y = alpha * std::pow(X, beta);
                else if (op == "ABS") Y = std::fabs(X);
                else if (op == "MULTIPLY") Y = Y * X;
            }
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
        } else if (kw == "SGWFN") {   // gas-water: Sg krg krw Pcgw
            sgwfn = tables(k, 4, ntsfun);
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
                    if (r.has(10) && inum(r, 10, 1) > 1) aqTabOf[a.id] = inum(r, 10, 1);   // C34 AQUTAB influence table
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
                a.coeff = num(r, 8, -1);   // C34: explicit influx coefficient replaces the face area (OPM)
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
        if (wc.status != WellStatus::Open) wc.reopen = false;
        lastCtrl[wc.well] = wc; lastCtrl[wc.well].reopen = false;   // C37: the WELOPEN event belongs to this step only
        const bool isInj = inj < 0 ? lastInj[wc.well] : inj == 1;
        lastInj[wc.well] = isInj; pendingInj[wc.well] = isInj;
        for (auto& p : pending) if (p.well == wc.well) { const bool ro = wc.status == WellStatus::Open && (p.reopen || wc.reopen); p = wc; p.reopen = ro; return; }
        pending.push_back(wc);
    }
    double fbhpProd = -1, fbhpInj = -1;   // C44 FBHPDEF default BHP limits in force (deck units; < 0 = ECLIPSE defaults)
    bool ctrlSeen = false;            // a native well-control keyword was read in the current report step
    void addStep(double len) {
        ctrlSeen = false;
        if (!(len > 0)) fail("time step lengths must be positive");
        ReportStep s; s.length = len; s.controls = pending; pending.clear();
        s.wpolymer = pendingPoly; pendingPoly.clear();
        s.drsdt = pendingDrsdt; s.drsdtFree = pendingDrsdtFree; pendingDrsdt = -2; pendingDrsdtFree = false;
        s.vap1 = pendingVap1; s.vap2 = pendingVap2; pendingVap1 = pendingVap2 = -1;
        s.vapAfterDrsdt = pendingVapAfterDrsdt; pendingVapAfterDrsdt = false;
        s.gprod = std::move(pendingGprod); s.ginj = std::move(pendingGinj); pendingGprod.clear(); pendingGinj.clear();
        s.guiderat = pendingGuide; pendingGuide = GuideRateSpec();
        s.welpi = std::move(pendingWelpi); pendingWelpi.clear();
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
    // C33 WLIST well lists ('*NAME'): the list state follows the schedule (NEW / ADD / DEL / MOV); a pattern starting
    // with '*' that names one or more lists (wildcards allowed, '*PROD?') expands to the wells of those lists
    std::map<std::string, std::vector<std::string>> wlists;
    // C46: a WELPI for well 'name' at a report step in [from, to) (its runtime scaling lives on the connection records)
    bool welpiSince(const std::string& name, int from, int to) const {
        for (int t = std::max(0, from); t < to && t < int(c.schedule.size()); ++t)
            for (const auto& wp : c.schedule[size_t(t)].welpi) if (matchName(wp.first, name)) return true;
        return false;
    }
    // C46 completion numbers (OPM / ECLIPSE): a connection's number is its 1-based position in the well's connection list
    // (order of first definition, one per layer) unless a COMPLUMP record covers it (records of the parser and the kept
    // COMPLUMP keywords up to step 'cur', so ACTIONX bodies see the deck's lumps); C1 / C2 <= 0 = no bound
    std::set<std::array<int, 4>> compNumberSel(const WellInput& w, int cur, int C1, int C2) const {
        std::vector<std::array<int, 4>> order;   // i j k lgr
        for (const Completion& q : w.completions) {
            if ((q.defStep >= 0 ? std::min(q.defStep, q.fromStep) : q.fromStep) > cur) continue;
            for (int kk = q.k1; kk <= q.k2; ++kk) {
                const std::array<int, 4> c4{q.i, q.j, kk, q.lgr};
                if (std::find(order.begin(), order.end(), c4) == order.end()) order.push_back(c4);
            }
        }
        std::vector<std::pair<std::string, std::array<int, 5>>> lumps = complump;
        for (const ExtraKeyword& e : c.extraKeywords) {
            if (e.name != "COMPLUMP" || e.step > cur) continue;
            for (const auto& rr : verbatimRecords(e.body)) {
                if (rr.size() < 6) continue;
                auto iv = [&](size_t i) { return rr[i].empty() || rr[i] == "*" || rr[i] == "1*" ? 0 : std::atoi(rr[i].c_str()); };
                std::string nm = rr[0]; if (nm.size() >= 2 && nm.front() == '\'' && nm.back() == '\'') nm = nm.substr(1, nm.size() - 2);
                lumps.push_back({nm, {iv(1), iv(2), iv(3), iv(4), iv(5)}});
            }
        }
        std::set<std::array<int, 4>> sel;
        for (size_t n = 0; n < order.size(); ++n) {
            int num = int(n) + 1;
            for (const auto& [pat, L] : lumps)
                if (matchName(pat, w.name) && (L[0] <= 0 || L[0] == order[n][0] + 1) && (L[1] <= 0 || L[1] == order[n][1] + 1) &&
                    (L[2] <= 0 || order[n][2] + 1 >= L[2]) && (L[3] <= 0 || order[n][2] + 1 <= L[3]))
                    num = L[4];
            if ((C1 <= 0 || num >= C1) && (C2 <= 0 || num <= C2)) sel.insert(order[n]);
        }
        return sel;
    }
    std::vector<WellInput*> wellsMatching(const std::string& pat) {
        std::vector<WellInput*> ws;
        if (pat.size() > 1 && pat[0] == '\\') {   // '\*RO*': escaped leading '*' = a well-name pattern, not a well list
            for (auto& w : c.wells) if (matchName(pat.substr(1), w.name)) ws.push_back(&w);
            return ws;
        }
        if (pat.size() > 1 && pat[0] == '*') {
            std::vector<std::string> names; bool anyList = false;
            for (const auto& [ln, members] : wlists)
                if (matchName(pat, ln)) { anyList = true; for (const auto& m : members) if (std::find(names.begin(), names.end(), m) == names.end()) names.push_back(m); }
            if (anyList) {
                for (const auto& n : names) for (auto& w : c.wells) if (w.name == n) ws.push_back(&w);
                return ws;
            }
        }
        for (auto& w : c.wells) if (matchName(pat, w.name)) ws.push_back(&w);
        return ws;
    }
    void wlistRecords(const Keyword& k) {
        for (const Record& r : k.records) {
            if (r.items.empty()) continue;
            std::string ln = str(r, 0); if (ln.empty() || ln[0] != '*') ln = "*" + ln;
            const std::string act = upper(str(r, 1, "NEW"));
            std::vector<std::string> ws;
            for (size_t i = 2; i < r.items.size(); ++i) {
                const std::string n = str(r, i);
                if (n.empty()) continue;
                bool found = false;
                for (WellInput* w : wellsMatching(n)) { ws.push_back(w->name); found = true; }
                if (!found) ws.push_back(n);   // a well defined later keeps its name
            }
            auto& L = wlists[ln];
            if (act == "NEW") L = ws;
            else if (act == "ADD") { for (auto& n : ws) if (std::find(L.begin(), L.end(), n) == L.end()) L.push_back(n); }
            else if (act == "DEL") { for (auto& n : ws) L.erase(std::remove(L.begin(), L.end(), n), L.end()); }
            else if (act == "MOV") {
                for (auto& [o, M] : wlists) if (o != ln) for (auto& n : ws) M.erase(std::remove(M.begin(), M.end(), n), M.end());
                for (auto& n : ws) if (std::find(L.begin(), L.end(), n) == L.end()) L.push_back(n);
            } else warn("WLIST: unknown action '" + act + "' - record ignored", r.line);
        }
    }

    void schedule(const Keyword& k) {
        if (k.name == "WCONPROD" || k.name == "WCONHIST" || k.name == "WCONINJE" || k.name == "WCONINJH" || k.name == "WELTARG" || k.name == "WELOPEN" || k.name == "WELOPENL")
            ctrlSeen = true;
        const std::string& kw = k.name;
        if (kw == "WELSPECS" || kw == "WELSPECL") {
            const bool local = kw == "WELSPECL";
            for (const Record& r0 : k.records) {
                if (r0.items.empty()) continue;
                int lgrIdx = -1;
                Record r = r0;
                if (local) {   // WELSPECL: item 3 is the local grid; i, j are local indices
                    const std::string ln = str(r0, 2);
                    lgrIdx = c.grid.lgrIndex(ln);
                    if (lgrIdx < 0) { warn("WELSPECL: local grid '" + ln + "' is not defined by CARFIN - record ignored", r0.line); continue; }
                    if (r.items.size() > 2) r.items.erase(r.items.begin() + 2);
                }
                const std::string name = str(r, 0);
                if (name.empty()) fail("well name missing", r.line);
                WellInput* w = findWell(name);
                const bool isNew = !w;
                if (!w) { c.wells.push_back(WellInput{}); w = &c.wells.back(); w->completions.clear(); w->name = name; }
                w->group = str(r, 1, "G");
                w->i = inum(r, 2, 1) - 1; w->j = inum(r, 3, 1) - 1;
                w->lgr = lgrIdx;
                w->refDepth = num(r, 4, -1);
                const std::string ph = upper(str(r, 5, "OIL"));
                w->prefPhase = ph == "WAT" ? std::string("WATER") : ph;
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
        } else if (kw == "COMPDAT" || kw == "COMPDATL") {
            const bool local = kw == "COMPDATL";
            for (const Record& r0 : k.records) {
                if (r0.items.empty()) continue;
                int lgrIdx = -1;
                Record r = r0;
                if (local) {   // COMPDATL: item 2 is the local grid; i, j, k are local indices
                    const std::string ln = str(r0, 1);
                    lgrIdx = c.grid.lgrIndex(ln);
                    if (lgrIdx < 0) { warn("COMPDATL: local grid '" + ln + "' is not defined by CARFIN - record ignored", r0.line); continue; }
                    if (r.items.size() > 1) r.items.erase(r.items.begin() + 1);
                }
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
                    cp.lgr = lgrIdx;
                    if (!local && w->lgr >= 0) grouped(1, "COMPDAT for a well declared with WELSPECL - the indices are taken as global (use COMPDATL)", w->name);
                    if (!local && (!r.has(1) || !r.has(2)) && w->lgr >= 0) { warn("COMPDAT: defaulted I/J of a local-grid well - record ignored (use COMPDATL)", r.line); continue; }
                    const int NX = lgrIdx >= 0 ? c.grid.lgrs[size_t(lgrIdx)].nx : c.grid.nx, NY = lgrIdx >= 0 ? c.grid.lgrs[size_t(lgrIdx)].ny : c.grid.ny,
                              NZ = lgrIdx >= 0 ? c.grid.lgrs[size_t(lgrIdx)].nz : c.grid.nz;
                    if (cp.i < 0 || cp.i >= NX || cp.j < 0 || cp.j >= NY || cp.k1 < 0 || cp.k2 >= NZ || cp.k1 > cp.k2) {
                        warn("connection of " + w->name + " outside the " + (local ? "local grid" : std::string("grid")) + " - ignored", r.line); continue;
                    }
                    // connections specified in later report steps open (or shut) at that step
                    const int cur = anyStep ? int(c.schedule.size()) : 0;
                    // the latest interval of this connection (a connection shut and re-opened has one entry per open period)
                    Completion* e = nullptr;
                    for (auto& x : w->completions)
                        if (x.i == cp.i && x.j == cp.j && x.k1 == cp.k1 && x.k2 == cp.k2 && x.lgr == cp.lgr && (!e || x.fromStep >= e->fromStep)) e = &x;
                    if (!e) { cp.fromStep = cur; w->completions.push_back(cp); continue; }   // SHUT ones too (WELOPEN may open them)
                    if (cur == 0) { *e = cp; continue; }
                    if (!cp.open) { if (e->open && e->untilStep < 0 && cur > e->fromStep) e->untilStep = cur; else if (e->fromStep >= cur) e->open = false; }
                    else if (!e->open) { const Completion old = *e; *e = cp; e->fromStep = cur; e->untilStep = -1; (void)old; }
                    else if (e->untilStep >= 0) {
                        const bool same = e->cf == cp.cf && e->kh == cp.kh && e->skin == cp.skin && e->diameter == cp.diameter && e->satTable == cp.satTable &&
                                          e->dfactor == cp.dfactor && e->r0 == cp.r0 && e->dir == cp.dir && e->mult == cp.mult;
                        if (e->untilStep >= cur && (same || e->fromStep >= cur)) { const int f = e->fromStep; *e = cp; e->fromStep = f; e->untilStep = -1; }   // shut and re-opened in the same step
                        else if (e->untilStep >= cur) { e->untilStep = cur; cp.fromStep = cur; cp.untilStep = -1; w->completions.push_back(cp); }   // ... with new data
                        else { cp.fromStep = cur; cp.untilStep = -1; w->completions.push_back(cp); }       // a new open period
                    } else if (e->fromStep >= cur) { const int f = e->fromStep; *e = cp; e->fromStep = f; }   // same step: replace
                    else if (e->cf != cp.cf || e->kh != cp.kh || e->skin != cp.skin || e->diameter != cp.diameter || e->satTable != cp.satTable ||
                             e->dfactor != cp.dfactor || e->r0 != cp.r0 || e->dir != cp.dir || e->mult != cp.mult ||
                             welpiSince(w->name, e->defStep >= 0 ? e->defStep : e->fromStep, cur)) {   // C44: COMPDAT resets WPIMULT; C46: and WELPI (OPM: a fresh connection)
                        // new data for an open connection: the old values stay in force until this step (not retroactive)
                        e->untilStep = cur; cp.fromStep = cur; cp.untilStep = -1; w->completions.push_back(cp);
                    }
                }
            }
        } else if (kw == "WCONPROD" || kw == "WCONHIST") {
            const bool hist = kw == "WCONHIST";
            for (const Record& r : k.records) {
                if (r.items.empty()) continue;
                WellControl wc; wc.well = str(r, 0); wc.status = statusOf(str(r, 1, "OPEN")); wc.history = hist;
                const std::string m = upper(str(r, 2, hist ? "RESV" : "NONE"));
                if (hist) { wc.orat = num(r, 3, -1); wc.wrat = num(r, 4, -1); wc.grat = num(r, 5, -1); }
                if (hist) { wc.bhpHist = num(r, 9, -1); wc.bhp = m == "BHP" ? wc.bhpHist : -1; wc.lrat = std::max(0.0, wc.orat) + std::max(0.0, wc.wrat);
                            wc.vfp = inum(r, 6, 0); wc.alq = num(r, 7, -1); wc.thp = num(r, 8, -1); }
                else { wc.orat = unum(r, 3, -1, wc, "ORAT"); wc.wrat = unum(r, 4, -1, wc, "WRAT"); wc.grat = unum(r, 5, -1, wc, "GRAT");
                       wc.lrat = unum(r, 6, -1, wc, "LRAT"); wc.resv = unum(r, 7, -1, wc, "RESV"); wc.bhp = unum(r, 8, -1, wc, "BHP");
                       wc.thp = unum(r, 9, -1, wc, "THP"); wc.vfp = inum(r, 10, 0); wc.alq = num(r, 11, -1); }
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
                    if (m == "THP" && wc.thp > 0 && wc.vfp > 0)
                        grouped(0, "THP control: BHP limit from the VFP table (rate limits stay in force)", wc.well);
                    else if (m == "GRUP")
                        grouped(0, "group control: the well shares its group's target (own rates are limits)", wc.well);
                    else if (wc.status == WellStatus::Open)
                        grouped(1, "control mode " + m + " is not modelled by the ResSim engine (kept for OPM Flow) - ResSim uses " +
                                   std::string(wc.mode == Control::BHP ? "BHP" : wc.mode == Control::ORAT ? "ORAT" : wc.mode == Control::LRAT ? "LRAT" :
                                               wc.mode == Control::WRAT ? "WRAT" : "GRAT"), wc.well);
                }
                if (wc.bhp <= 0 && fbhpProd > 0 && !wc.udaOf("BHP")) { wc.bhp = fbhpProd; wc.bhpDefault = true; }   // C44 FBHPDEF
                if (wc.mode == Control::BHP && wc.bhp <= 0 && wc.status == WellStatus::Open && !wc.udaOf("BHP")) {
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
                    wc.rate = num(r, 3, 0); wc.bhpHist = num(r, 4, -1); wc.bhp = wc.mode == Control::BHP ? wc.bhpHist : -1;
                } else {
                    const std::string m = upper(str(r, 3, "RATE"));
                    wc.rate = unum(r, 4, -1, wc, "RATE"); wc.bhp = unum(r, 6, -1, wc, "BHP");
                    wc.resv = unum(r, 5, -1, wc, "RESV"); wc.thp = unum(r, 7, -1, wc, "THP"); wc.vfp = inum(r, 8, 0);
                    if (m == "BHP") wc.mode = Control::BHP;
                    else if (m == "RATE") wc.mode = Control::RATE;
                    else if (m == "RESV") { grouped(1, "RESV injection: ResSim uses it as a surface RATE (kept as RESV for OPM Flow)", wc.well); wc.deckMode = m;
                                            wc.mode = Control::RATE; if (wc.rate <= 0) { wc.rate = wc.resv > 0 ? wc.resv : 0; wc.deckMode = "RESV*"; } }
                    else {
                        if (m == "THP" && wc.thp > 0 && wc.vfp > 0) grouped(0, "THP control: BHP limit from the VFP table (rate limit stays in force)", wc.well);
                        else if (m == "GRUP") grouped(0, "group control: the injector shares its group's target", wc.well);
                        else grouped(1, "injection control " + m + " is not modelled by the ResSim engine (kept for OPM Flow) - ResSim uses " + (wc.rate > 0 ? std::string("RATE") : std::string("BHP")), wc.well);
                        wc.deckMode = m; wc.mode = wc.rate > 0 ? Control::RATE : Control::BHP;
                    }
                }
                if (wc.bhp <= 0 && fbhpInj > 0 && !wc.udaOf("BHP")) { wc.bhp = fbhpInj; wc.bhpDefault = true; }   // C44 FBHPDEF
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
                WellControl udaProbe; udaProbe.well = pat;
                static const std::set<std::string> udaTargets = {"ORAT", "WRAT", "GRAT", "LRAT", "RESV", "BHP", "RATE", "THP"};
                const double v = udaTargets.count(what) ? unum(r, 2, NaN, udaProbe, what.c_str()) : num(r, 2, NaN);
                if (std::isnan(v) && udaProbe.uda.empty()) continue;
                static const std::set<std::string> native = {"ORAT", "WRAT", "GRAT", "LRAT", "RESV", "BHP", "RATE", "THP", "VFP", "LIFT"};
                if (!native.count(what)) { keepRecord(k, r, "WELTARG target"); grouped(1, "WELTARG " + what + " is not modelled by the ResSim engine - kept for OPM Flow", pat); continue; }
                for (WellInput* w : wellsMatching(pat)) {
                    auto it = lastCtrl.find(w->name);
                    if (it == lastCtrl.end()) { grouped(1, "no previous control for the well - ignored", w->name); continue; }
                    WellControl wc = it->second;
                    wc.uda.erase(std::remove_if(wc.uda.begin(), wc.uda.end(), [&](const std::pair<std::string, std::string>& u) { return u.first == what; }), wc.uda.end());
                    if (!udaProbe.uda.empty()) { wc.uda.push_back(udaProbe.uda[0]); addControl(wc); continue; }
                    if (what == "ORAT") wc.orat = v; else if (what == "WRAT") wc.wrat = v; else if (what == "GRAT") wc.grat = v;
                    else if (what == "LRAT") wc.lrat = v; else if (what == "RESV") wc.resv = v; else if (what == "BHP") { wc.bhp = v; wc.bhpDefault = false; }
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
                pendingVapAfterDrsdt = false;
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
            // C44: ECLIPSE / OPM semantics - the factor multiplies the connection factors of the connections that exist (open) at
            // this report step, cumulatively; connections defined later, or redefined by COMPDAT, do not inherit it. A connection
            // open before this step gets a new open period from this step (the earlier period keeps its factor).
            const int cur = anyStep ? int(c.schedule.size()) : 0;
            for (const Record& r : k.records) {
                if (r.items.empty()) continue;
                const double m = num(r, 1, 1.0);
                const int I = inum(r, 2, 0), J = inum(r, 3, 0), K = inum(r, 4, 0), C1 = inum(r, 5, 0), C2 = inum(r, 6, 0);
                const bool byComp = C1 > 0 || C2 > 0;   // C46: completion numbers (COMPLUMP / definition order), one layer at a time
                for (WellInput* w : wellsMatching(str(r, 0))) {
                    const std::set<std::array<int, 4>> sel = byComp ? compNumberSel(*w, cur, C1, C2) : std::set<std::array<int, 4>>();
                    auto inSel = [&](const Completion& q, int kk) { return sel.count(std::array<int, 4>{q.i, q.j, kk, q.lgr}) > 0; };
                    auto matchIJK = [&](const Completion& q) {
                        if (q.lgr >= 0 && (I > 0 || J > 0 || K > 0)) return false;
                        return !((I > 0 && q.i + 1 != I) || (J > 0 && q.j + 1 != J) || (K > 0 && (K < q.k1 + 1 || K > q.k2 + 1)));
                    };
                    auto match = [&](const Completion& q) { return matchIJK(q) && (!byComp || (q.k1 == q.k2 && inSel(q, q.k1))); };
                    if (K > 0 || byComp)   // one layer of a multi-layer record: split it into one record per layer (as WELOPEN)
                        for (size_t qi = 0; qi < w->completions.size(); ++qi) {
                            Completion& q = w->completions[qi];
                            if (q.k1 == q.k2 || q.fromStep > cur || !matchIJK(q)) continue;
                            if (byComp) { bool any = false; for (int kk = q.k1; kk <= q.k2; ++kk) any = any || inSel(q, kk); if (!any) continue; }
                            const Completion base = q;
                            std::vector<Completion> parts;
                            for (int kk = base.k1; kk <= base.k2; ++kk) { Completion n = base; n.k1 = n.k2 = kk; parts.push_back(n); }
                            w->completions[qi] = parts[0];
                            w->completions.insert(w->completions.begin() + std::ptrdiff_t(qi) + 1, parts.begin() + 1, parts.end());
                            qi += parts.size() - 1;
                        }
                    std::vector<size_t> hit;
                    for (size_t qi = 0; qi < w->completions.size(); ++qi) {
                        const Completion& q = w->completions[qi];
                        // C46: connections that exist SHUT at this step are scaled too (OPM) - also those a later WELOPEN opens
                        // (already parsed when an ACTIONX body runs: defined at defStep, open from a future fromStep)
                        const bool existsShut = (!q.open && q.fromStep <= cur) || (q.open && q.defStep >= 0 && q.defStep <= cur && q.fromStep > cur);
                        if (existsShut && match(q)) { w->completions[qi].mult *= m; continue; }
                        if (q.fromStep > cur || !q.open || (q.untilStep >= 0 && q.untilStep <= cur) || !match(q)) continue;
                        hit.push_back(qi);
                    }
                    for (auto it = hit.rbegin(); it != hit.rend(); ++it) {   // descending: inserts keep the lower indices valid
                        Completion& q = w->completions[*it];
                        if (q.fromStep >= cur) { q.mult *= m; continue; }
                        // C46: the new period keeps a closure already scheduled later (an ACTIONX body runs after the whole deck
                        // was parsed: a later WELOPEN SHUT / COMPDAT SHUT of this connection must survive the split)
                        Completion n = q; n.fromStep = cur; n.untilStep = q.untilStep > cur ? q.untilStep : -1; n.mult *= m;
                        q.untilStep = cur;
                        w->completions.insert(w->completions.begin() + std::ptrdiff_t(*it) + 1, n);   // deck order kept (topmost connection)
                    }
                }
            }
        } else if (kw == "WELOPEN" || kw == "WELOPENL") {
            const bool localOpen = kw == "WELOPENL";
            for (const Record& r1 : k.records) {
                if (r1.items.empty()) continue;
                Record r = r1;
                int openLgr = -1;
                if (localOpen) {   // WELOPENL: item 2 is the local grid; I J K are local indices of that LGR
                    openLgr = c.grid.lgrIndex(str(r1, 1));
                    if (openLgr < 0) { warn("WELOPENL: local grid '" + str(r1, 1) + "' is not defined by CARFIN - record ignored", r1.line); continue; }
                    if (r.items.size() > 1) r.items.erase(r.items.begin() + 1);
                }
                const WellStatus st = statusOf(str(r, 1, "OPEN"));
                if (r.has(2) || r.has(3) || r.has(4) || r.has(5) || r.has(6)) {
                    // connection-level: I J K (0 / default = any); completion numbers (COMPLUMP) are kept verbatim
                    const int I = inum(r, 2, 0), J = inum(r, 3, 0), K = inum(r, 4, 0), C1 = inum(r, 5, 0), C2 = inum(r, 6, 0);
                    const int cur = anyStep ? int(c.schedule.size()) : 0;
                    // C46 completion numbers (OPM / ECLIPSE): a connection's number is its 1-based position in the well's
                    // connection list (order of first definition, one per layer) unless a COMPLUMP record covers it;
                    // C1 / C2 defaulted = no bound. The selected connections are opened / shut one layer at a time.
                    const bool byComp = C1 > 0 || C2 > 0;
                    std::map<const WellInput*, std::set<std::array<int, 4>>> compSel;
                    if (byComp) for (WellInput* w : wellsMatching(str(r, 0))) compSel[w] = compNumberSel(*w, cur, C1, C2);
                    auto selected = [&](const WellInput* w, const Completion& q, int kk) {
                        auto it = compSel.find(w);
                        return it != compSel.end() && it->second.count(std::array<int, 4>{q.i, q.j, kk, q.lgr}) > 0;
                    };
                    bool representable = true;
                    // C44: one layer of a multi-layer COMPDAT record -> split the record into one record per layer (the COMPDAT
                    // connection factor / Kh items are per connection, so the split changes nothing else; order k1..k2 kept)
                    if (K > 0 || byComp)
                        for (WellInput* w : wellsMatching(str(r, 0)))
                            for (size_t qi = 0; qi < w->completions.size(); ++qi) {
                                Completion& q = w->completions[qi];
                                if (q.k1 == q.k2 || q.fromStep > cur) continue;
                                if (K > 0 && (K < q.k1 + 1 || K > q.k2 + 1)) continue;
                                if (byComp) { bool any = false; for (int kk = q.k1; kk <= q.k2; ++kk) any = any || selected(w, q, kk); if (!any) continue; }
                                if (localOpen && q.lgr != openLgr) continue;
                                if ((I > 0 && q.i + 1 != I) || (J > 0 && q.j + 1 != J)) continue;
                                const Completion base = q;
                                std::vector<Completion> parts;
                                for (int kk = base.k1; kk <= base.k2; ++kk) { Completion n = base; n.k1 = n.k2 = kk; parts.push_back(n); }
                                w->completions[qi] = parts[0];
                                w->completions.insert(w->completions.begin() + std::ptrdiff_t(qi) + 1, parts.begin() + 1, parts.end());
                                qi += parts.size() - 1;
                            }
                    std::vector<std::pair<WellInput*, size_t>> hit;   // latest interval of each matching connection
                    for (WellInput* w : wellsMatching(str(r, 0)))
                        for (size_t qi = 0; qi < w->completions.size(); ++qi) {
                            const Completion& q = w->completions[qi];
                            if (q.fromStep > cur) continue;
                            if (localOpen && q.lgr != openLgr) continue;
                            if ((I > 0 && q.i + 1 != I) || (J > 0 && q.j + 1 != J) || (K > 0 && (K < q.k1 + 1 || K > q.k2 + 1))) continue;
                            if (byComp && (q.k1 != q.k2 || !selected(w, q, q.k1))) continue;
                            bool later = false;
                            for (const Completion& x : w->completions)
                                if (&x != &q && x.i == q.i && x.j == q.j && x.k1 == q.k1 && x.k2 == q.k2 && x.lgr == q.lgr && x.fromStep > q.fromStep && x.fromStep <= cur) later = true;
                            if (later) continue;
                            if (K > 0 && q.k1 != q.k2) representable = false;      // one layer of a multi-layer connection
                            hit.push_back({w, qi});
                        }
                    if (!representable) { keepRecord(k, r1, "connection-level WELOPEN"); grouped(1, "WELOPEN on part of a connection range: not modelled by the ResSim engine - kept for OPM Flow", str(r, 0)); continue; }
                    for (auto& [w, qi] : hit) {
                        Completion* q = &w->completions[qi];
                        if (st == WellStatus::Open) {
                            if (!q->open) { if (q->defStep < 0 && q->fromStep < cur) q->defStep = std::max(0, q->fromStep); q->open = true; q->fromStep = cur; q->untilStep = -1; }
                            else if (q->untilStep >= 0 && q->untilStep == cur) q->untilStep = -1;                 // shut and re-opened in the same step
                            else if (q->untilStep >= 0 && q->untilStep < cur) { Completion n2 = *q; n2.fromStep = cur; n2.untilStep = -1; w->completions.push_back(n2); }   // new open period
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
                    wc.reopen = st == WellStatus::Open;   // C37 WCYCLE: WELOPEN OPEN restarts the on-period
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
            t.growth = std::max(1.2, num(r1, 4, 3.0));   // C38: TSFMAX defaults to 3 when TUNING is given (ECLIPSE / OPM) t.cutFactor = std::min(0.5, std::max(0.1, num(r1, 6, t.cutFactor)));
            t.cnvTol = std::max(1e-5, num(r2, 1, t.cnvTol)); t.mbTol = std::max(1e-9, num(r2, 2, t.mbTol));
            t.maxNewton = std::max(t.maxNewton, inum(r3, 0, t.maxNewton));
            if (num(r3, 7, 1e6) < 1.0) t.dsMax = num(r3, 7, t.dsMax);
            if (anyStep) grouped(0, "TUNING after the first report step applies to the whole run");
        }
    }

    // ------------------------------------------------------------------ RESSIM extension lines
    void ressim(const Record& r) {
        const std::string key = upper(str(r, 0));
        if (key == "TUNING") { c.tuning.dpMaxRel = num(r, 1, c.tuning.dpMaxRel); c.tuning.threads = inum(r, 2, 0); c.tuning.linearSolver = inum(r, 3, 2); c.tuning.cprWeights = inum(r, 4, 1) ? 1 : 0;
            c.tuning.mixedPrecision = inum(r, 5, 0) != 0; c.tuning.gpu = inum(r, 6, 0); c.tuning.amgReuse = std::max(0, inum(r, 7, 10));
            c.tuning.dsTarget = std::max(0.0, num(r, 8, 0.0)); c.tuning.dpTarget = std::max(0.0, num(r, 9, 0.0));
            c.tuning.leanMemory = std::clamp(inum(r, 10, -1), -1, 1); }
        else if (key == "AQVISC") aqVisc[inum(r, 1, 0)] = num(r, 2, -1);
        else if (key == "PARAM") {
            StudyParam p; p.kind = ParamKind(std::clamp(inum(r, 1, 0), 0, kNumParamKinds - 1)); p.target = inum(r, 2, 0);
            p.i1 = inum(r, 3, -1); p.i2 = inum(r, 4, -1); p.j1 = inum(r, 5, -1); p.j2 = inum(r, 6, -1);
            p.k1 = inum(r, 7, -1); p.k2 = inum(r, 8, -1);
            p.base = num(r, 9, 1); p.low = num(r, 10, 0.5); p.high = num(r, 11, 2);
            p.log = inum(r, 12, 1) != 0; p.active = inum(r, 13, 1) != 0; p.label = str(r, 14);
            c.studyParams.push_back(p);
        } else if (key == "PROPMODEL") {   // C40 property model: prop method model nugget sill rmax rmin rvert az seed maxnbr fromwells minV maxV pp a b kvkh 'source'
            PropModelInput m;
            m.property = std::clamp(inum(r, 1, 0), 0, 3); m.method = inum(r, 2, 1); m.model = std::clamp(inum(r, 3, 0), 0, 2);
            m.nugget = num(r, 4, 0); m.sill = num(r, 5, 1); m.rangeMax = num(r, 6, 1000); m.rangeMin = num(r, 7, 1000); m.rangeVert = num(r, 8, 20);
            m.azimuthDeg = num(r, 9, 0); m.seed = unsigned(std::max(1, inum(r, 10, 1))); m.maxNbr = std::max(1, inum(r, 11, 16));
            m.fromWells = inum(r, 12, 0) != 0; m.minV = num(r, 13, -1e300); m.maxV = num(r, 14, 1e300);
            m.permFromPoro = inum(r, 15, 0) != 0; m.ppA = num(r, 16, 0); m.ppB = num(r, 17, 0); m.kvkh = num(r, 18, 0.1); m.source = str(r, 19);
            c.propModels.push_back(m);
        } else if (key == "PROPFAC") {     // C41 facies of the last PROPMODEL: satnum(0/1) then proportion poro perm per facies
            if (!c.propModels.empty()) {
                PropModelInput& m = c.propModels.back();
                m.facSatnum = inum(r, 1, 0) != 0; m.facProp.clear(); m.facPoro.clear(); m.facPerm.clear();
                for (size_t q = 2; q + 2 < r.items.size(); q += 3) { m.facProp.push_back(num(r, q, 0)); m.facPoro.push_back(num(r, q + 1, -1)); m.facPerm.push_back(num(r, q + 2, -1)); }
            }
        } else if (key == "PROPDATA") {    // x y z value of the last PROPMODEL (repeated: several points per line)
            if (!c.propModels.empty())
                for (size_t q = 1; q + 3 < r.items.size(); q += 4)
                    c.propModels.back().data.push_back({num(r, q, 0), num(r, q + 1, 0), num(r, q + 2, 0), num(r, q + 3, 0)});
        } else if (key == "PIMULT") {     // C44: whole-run PI multiplier; the WPIMULT records of the deck already contain it
            const double v = num(r, 2, 1.0);
            if (v > 0) ressimPiMult[str(r, 1)] = v;
        } else if (key == "OBS") {
            Observation ob; const std::string q = upper(str(r, 1)); bool f = false;
            for (int i = 0; i <= int(ObsQty::FGOR); ++i) if (q == toString(ObsQty(i))) { ob.qty = ObsQty(i); f = true; }
            if (!f) { warn("unknown observation quantity " + q); return; }
            ob.well = str(r, 2, "-"); if (ob.well == "-") ob.well.clear();
            ob.time = num(r, 3, 0); ob.value = num(r, 4, 0); ob.sigma = num(r, 5, -1); ob.weight = num(r, 6, 1);
            c.observations.push_back(ob);
        } else if (key == "HFRAC") {
            // RESSIM: HFRAC name x0 y0 x1 y1 zTop zBot width perm poro [gamma] [satnum]   (C15 hydraulic fracture, EDFM)
            GridInput::HFrac F;
            F.name = str(r, 1, "HF" + std::to_string(c.grid.hfracs.size() + 1));
            F.x0 = num(r, 2, 0); F.y0 = num(r, 3, 0); F.x1 = num(r, 4, 0); F.y1 = num(r, 5, 0);
            F.zTop = num(r, 6, -1); F.zBot = num(r, 7, -1); F.width = num(r, 8, F.width); F.perm = num(r, 9, 0);
            F.poro = num(r, 10, F.poro); F.gamma = num(r, 11, 0); F.satnum = inum(r, 12, 0);
            if (!(F.width > 0) || (F.x0 == F.x1 && F.y0 == F.y1)) { warn("HFRAC " + F.name + ": zero length or width - ignored"); return; }
            c.grid.hfracs.push_back(F);
        } else if (key == "HFGROW") {
            // RESSIM: HFGROW name KIC sigmaMin [sigmaGrad] [initLen] [closedMult]   (C18 in-run LEFM growth of HFRAC 'name')
            const std::string nm = str(r, 1, "");
            GridInput::HFrac* F = nullptr;
            for (auto& h : c.grid.hfracs) if (h.name == nm) F = &h;
            if (!F) { warn("HFGROW " + nm + ": no HFRAC of that name before this line - ignored"); return; }
            F->kic = std::max(0.0, num(r, 2, 1.0)); F->sigma = num(r, 3, 0); F->sigmaGrad = num(r, 4, 0); F->initLen = std::max(0.0, num(r, 5, 0));
            F->closedMult = std::clamp(num(r, 6, 1.0), 1e-12, 1e6);
            for (size_t q = 2; q < r.items.size(); ++q) {   // C21 'SUBGRID' flag anywhere after the numbers
                std::string t = str(r, int(q), "");
                for (auto& ch : t) ch = char(std::toupper((unsigned char)ch));
                if (t == "SUBGRID") F->subgrid = 1;
                if (t == "LEAK") F->subLeak = std::max(F->subLeak, 1);
                if (t == "LEAKC") F->subLeak = 2;   // C48 H2e: gas-aware sub-cell storage (compressibility ratio)
                if (t == "FESTRESS") F->feStress = 1;   // C25: closure stress from the FE geomechanics (GEOMECH)
                if (t == "PHASEFIELD") {                // C28: in-run phase-field growth [eta [h_m [E_GPa nu]]]
                    F->phaseField = 1; F->subgrid = 0;
                    if (r.has(q + 1) && r.isNumber(q + 1)) F->pfEta = std::clamp(r.num(q + 1, 0.5), 1e-4, 1.0);
                    if (r.has(q + 2) && r.isNumber(q + 2)) F->pfH = std::max(0.0, r.num(q + 2, 0));
                    if (r.has(q + 3) && r.isNumber(q + 3)) F->pfE = std::max(0.0, r.num(q + 3, 0));
                    if (r.has(q + 4) && r.isNumber(q + 4)) F->pfNu = std::clamp(r.num(q + 4, 0.25), 0.0, 0.49);
                }
                if (t == "TWOWAY") F->pfTwoWay = 1;     // C29 two-way phase-field coupling (pressure-driven volume + aperture)
                if (t == "DDMWIDTH") {                  // C26: mechanical aperture from the 3-D DDM opening solve [E_GPa nu]
                    F->ddmWidth = 1;
                    if (r.has(q + 1) && r.isNumber(q + 1)) F->ddmE = std::max(0.0, r.num(q + 1, 0));
                    if (r.has(q + 2) && r.isNumber(q + 2)) F->ddmNu = std::clamp(r.num(q + 2, 0.25), 0.0, 0.49);
                }
            }
            if (F->phaseField) F->subgrid = 0;   // the phase-field model sets the wing length itself
        } else if (key == "STRESSFIELD") {
            // RESSIM: STRESSFIELD svGrad shminRatio shmaxRatio azimuthDeg [stressPath] [sv0]   (C24, total stress, case units)
            // C25: ... [FE] [TH thE_GPa thAlpha thNu]  (FE: stress tensor of the finite-element geomechanics; TH: thermo-elastic
            // horizontal stress of thermal runs with the analytic field)
            auto& F = c.faultMech;
            auto nq = [&](size_t i, double d) { return r.has(i) && r.isNumber(i) ? r.num(i, d) : d; };
            size_t nNum = r.items.size();
            for (size_t q = 1; q < r.items.size(); ++q) {
                std::string t = str(r, q, ""); for (auto& ch : t) ch = char(std::toupper((unsigned char)ch));
                if (t == "FE" || t == "TH") { nNum = std::min(nNum, q); }
                if (t == "FE") F.source = 1;
                if (t == "TH") { F.thE = std::max(0.0, nq(q + 1, 0)); F.thAlpha = std::max(0.0, nq(q + 2, 0)); F.thNu = std::clamp(nq(q + 3, 0.25), 0.0, 0.49); }
            }
            auto nb = [&](size_t i, double d) { return i < nNum ? num(r, i, d) : d; };
            F.svGrad = std::max(0.0, nb(1, 0)); F.shRatio = nb(2, F.shRatio); F.sHRatio = nb(3, F.sHRatio);
            F.azimuth = nb(4, 0); F.stressPath = nb(5, F.stressPath); F.sv0 = nb(6, 0);
        } else if (key == "FAULTMECH") {
            // RESSIM: FAULTMECH name friction cohesion Asigma ta_years r0 bvalue Mmin [permFactor [ramp]]   (C24; C49 ramp width
            // of the two-way factor as a fraction of sigma'_n, default 0.05)
            FaultMechInput::Row R; R.name = str(r, 1, "");
            if (R.name.empty()) { warn("FAULTMECH without a fault name - ignored"); return; }
            R.mu = num(r, 2, R.mu); R.coh = num(r, 3, 0); R.aSig = std::max(0.0, num(r, 4, 0)); R.ta = std::max(1e-6, num(r, 5, R.ta));
            R.r0 = std::max(0.0, num(r, 6, R.r0)); R.bValue = std::max(0.1, num(r, 7, R.bValue)); R.mMin = num(r, 8, R.mMin);
            R.permFactor = std::max(1e-6, num(r, 9, 1.0));
            R.permRamp = std::max(1e-4, num(r, 10, 0.05));
            for (auto& q : c.faultMech.faults) if (q.name == R.name) { const int L = q.rsfLayer; const double a = q.rsfA, b = q.rsfB, d = q.rsfDcMm; q = R; q.rsfLayer = L; q.rsfA = a; q.rsfB = b; q.rsfDcMm = d; return; }
            c.faultMech.faults.push_back(R);
        } else if (key == "TRACERPROP") {
            // RESSIM: TRACERPROP name halfLife_days kd [init] [K] [parent] [alphaL] [Dm]   (C27 decay / adsorption / init, C28 partitioning / chain / dispersion)
            Case::TracerProp T; T.name = upper(str(r, 1, ""));
            if (T.name.empty()) { grouped(1, "RESSIM TRACERPROP needs a tracer name", ""); return; }
            T.halfLife = std::max(0.0, num(r, 2, 0.0)); T.kd = std::max(0.0, num(r, 3, 0.0)); T.init = std::max(0.0, num(r, 4, 0.0));
            // C28: [K_partition] [parent | 1*] [alphaL] [Dm]
            T.kPart = std::max(0.0, num(r, 5, 0.0));
            { std::string pa = upper(str(r, 6, "")); if (!pa.empty() && pa != "1*" && pa != "*" && pa != "NONE" && pa != T.name) T.parent = pa; }
            T.alphaL = std::max(0.0, num(r, 7, 0.0)); T.dm = std::max(0.0, num(r, 8, 0.0));
            for (auto& q : c.tracerProps) if (q.name == T.name) { q = T; return; }
            c.tracerProps.push_back(T);
        } else if (key == "WELLBORE") {
            // RESSIM: WELLBORE well ID roughness [icdK [rhoCal]]   (C29 wellbore friction + inflow-control devices)
            Case::WellboreSpec B; B.well = upper(str(r, 1, ""));
            if (B.well.empty()) { grouped(1, "RESSIM WELLBORE needs a well name", ""); return; }
            B.id = std::max(1e-4, num(r, 2, c.units == UnitSystem::Field ? 0.333 : 0.1));
            B.rough = std::max(0.0, num(r, 3, c.units == UnitSystem::Field ? 5e-5 : 1.5e-5));
            B.icdK = std::max(0.0, num(r, 4, 0.0)); B.rhoCal = std::max(0.0, num(r, 5, 0.0));
            for (auto& q : c.wellbores) if (q.well == B.well) { q = B; return; }
            c.wellbores.push_back(B);
        } else if (key == "MSWOPTS") {
            // RESSIM: MSWOPTS flash storage [xflow [heads]]   (C37 multisegment wells: tubing flash, segment storage; C45 cross-flow;
            // 1 / 0; C49 heads: 1 converged, 0 explicit, -1 / default converged)
            c.mswOpt.flash = inum(r, 1, 1) != 0; c.mswOpt.storage = inum(r, 2, 1) != 0; c.mswOpt.xflow = inum(r, 3, 0) != 0;
            c.mswOpt.heads = std::max(-1, std::min(1, inum(r, 4, -1)));
        } else if (key == "SALT") {
            // RESSIM: SALT c0 [cSol [rhoSalt [phiR [Gamma [T_C [cInj]]]]]]   (C29 EWASG-type salt precipitation / dry-out)
            auto& S = c.salt; S.on = true;
            S.c0 = std::max(0.0, num(r, 1, 0.0)); S.cSol = std::max(1e-6, num(r, 2, S.cSol)); S.rhoSalt = std::max(1.0, num(r, 3, S.rhoSalt));
            S.phiR = std::clamp(num(r, 4, S.phiR), 0.01, 0.99); S.gamma = std::clamp(num(r, 5, S.gamma), 0.01, 0.99);
            S.tempC = std::clamp(num(r, 6, S.tempC), 0.0, 370.0); S.cInj = std::max(0.0, num(r, 7, 0.0));
        } else if (key == "HFFIELD") {
            // RESSIM: HFFIELD ON|OFF [E_GPa nu] [radius] [budget]   (C30 fracture stress / opening field output)
            auto& H = c.hfField;
            const std::string m = upper(str(r, 1, "ON"));
            H.on = (m == "OFF" || m == "0" || m == "NO") ? 0 : 1;
            auto nq = [&](size_t i, double d) { return r.has(i) && r.isNumber(i) ? r.num(i, d) : d; };
            H.E = std::max(0.0, nq(2, 0.0)); H.nu = std::clamp(nq(3, 0.25), 0.0, 0.49);
            H.radius = std::clamp(nq(4, 1.0), 0.05, 20.0); H.budget = std::clamp(nq(5, 6.0), 0.01, 1000.0);
        } else if (key == "FAULTSLIP") {
            // RESSIM: FAULTSLIP E_GPa nu dMu [slipRef_mm] [cutoff] [FULLSPACE]   (C26 DDM stress transfer + quasi-static slip)
            auto& F = c.faultMech;
            F.slipOn = true;
            auto nq = [&](size_t i, double d) { return r.has(i) && r.isNumber(i) ? r.num(i, d) : d; };
            F.slipE = std::max(1e-3, nq(1, F.slipE)); F.slipNu = std::clamp(nq(2, F.slipNu), 0.0, 0.49); F.slipDMu = std::max(0.0, nq(3, 0.0));
            F.slipRefMm = std::max(1e-6, nq(4, F.slipRefMm)); F.slipCutoff = nq(5, F.slipCutoff);
            for (size_t q = 1; q < r.items.size(); ++q) {
                std::string t = str(r, q, ""); for (auto& ch : t) ch = char(std::toupper((unsigned char)ch));
                if (t == "FULLSPACE") F.slipHalfSpace = 0;
                if (t == "HALFSPACE") F.slipHalfSpace = 1;
            }
        } else if (key == "FAULTRSF") {
            // RESSIM: FAULTRSF name layer a b Dc_mm   (C24 quasi-dynamic rate-and-state trace of FAULTMECH 'name')
            const std::string nm = str(r, 1, "");
            FaultMechInput::Row* R = nullptr;
            for (auto& q : c.faultMech.faults) if (q.name == nm) R = &q;
            if (!R) { warn("FAULTRSF " + nm + ": no FAULTMECH of that name before this line - ignored"); return; }
            R->rsfLayer = std::max(0, inum(r, 2, 1)); R->rsfA = num(r, 3, R->rsfA); R->rsfB = num(r, 4, R->rsfB); R->rsfDcMm = std::max(1e-6, num(r, 5, R->rsfDcMm));
        } else if (key == "GEOCOUPLE") {
            // RESSIM: GEOCOUPLE sxx syy E_GPa nu friction layer   (C24 G1-c in the main run; stresses compression +, case units)
            auto& G = c.geoCouple; G.on = true;
            G.sxx = num(r, 1, 0); G.syy = num(r, 2, 0); G.eGPa = std::max(1e-3, num(r, 3, G.eGPa)); G.nu = std::clamp(num(r, 4, G.nu), 0.0, 0.49);
            G.friction = std::max(0.0, num(r, 5, G.friction)); G.layer = std::max(1, inum(r, 6, 1));
        } else if (parseGeoRecord(key, r, c.geoSeq)) {
            // RESSIM: GEOMECH / MECHBC / GEOPLAST of a black-oil / compositional deck (C25 sequential FE geomechanics);
            // a thermal deck reads them into thermal.geo (thermalFromDeck) - the copy here is not exported for thermal cases
        } else if (key == "LANGMUIR") {
            // RESSIM: LANGMUIR VL PL rhoGrain   (C15 adsorbed gas: VL sm3/t | scf/ton, PL pressure units, grain density g/cm3)
            c.grid.langVL = std::max(0.0, num(r, 1, 0)); c.grid.langPL = std::max(0.0, num(r, 2, 0)); c.grid.langRho = std::max(0.0, num(r, 3, 2.6));
        } else if (key == "PERMSTRESS") {
            // RESSIM: PERMSTRESS gamma   (C15 matrix permeability k = k0 exp(-gamma (p_i - p)), 1 / pressure units)
            c.grid.permGamma = std::max(0.0, num(r, 1, 0));
        } else if (key == "MINC") {
            // RESSIM: MINC <shells> [<fracture directions 1|2|3>] [<volume fraction outer> ... <inner>]
            c.grid.minc = std::max(0, inum(r, 1, 0)); c.grid.mincDims = std::clamp(inum(r, 2, 3), 1, 3);
            c.grid.mincVol.clear();
            for (int j = 0; j < c.grid.minc && int(r.items.size()) > 3 + j; ++j) c.grid.mincVol.push_back(num(r, 3 + j, 0));
            if (!c.grid.mincVol.empty() && int(c.grid.mincVol.size()) != c.grid.minc) { warn("MINC: one volume fraction per shell - geometric default used"); c.grid.mincVol.clear(); }
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
                                                      "SWOF", "SGOF", "SLGOF", "SWFN", "SGFN", "SOF3", "SOF2", "SGWFN"};
        static const std::set<std::string> solKws = {"EQUIL", "RSVD", "PBVD", "AQUFETP", "AQUCT", "AQUANCON", "THPRES", "RVVD"};
        static const std::set<std::string> schKws = {"WELSPECS", "WELSPECL", "COMPDAT", "COMPDATL", "WCONPROD", "WCONHIST", "WCONINJE", "WCONINJH", "WELTARG",
                                                     "WPIMULT", "WELOPEN", "WELOPENL", "TSTEP", "DATES", "TUNING", "WPOLYMER", "DRSDT"};
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
                if (sections.count(kw)) { leaveLgr(); section = kw; continue; }
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
                if (kw == "NONNC") { c.grid.nonnc = true; grouped(0, "NONNC: no non-neighbour connections across fault throws (C38)"); continue; }
                if (kw == "GRIDUNIT") {   // C37: unit of the corner-point geometry (exported grids often come in METRES / CM)
                    std::string u = k.records.empty() ? std::string("METRES") : upper(str(k.records[0], 0));
                    if (u.empty()) u = "METRES";
                    const double toM = u.rfind("FEET", 0) == 0 ? 0.3048 : u.rfind("CM", 0) == 0 ? 0.01 : u.rfind("MET", 0) == 0 ? 1.0 : -1.0;
                    if (toM < 0) { warn("GRIDUNIT: unknown unit '" + u + "' - METRES assumed"); gridUnitFac = c.units == UnitSystem::Field ? 1.0 / 0.3048 : 1.0; continue; }
                    gridUnitFac = c.units == UnitSystem::Field ? toM / 0.3048 : toM;
                    if (std::fabs(gridUnitFac - 1.0) > 1e-12)
                        grouped(0, "GRIDUNIT: COORD / ZCORN given in " + u + " - converted to the deck unit system (factor " + std::to_string(gridUnitFac) + ")");
                    continue;
                }
                if (kw == "POLYMER") { c.phases.polymer = true; continue; }
                // ---- hysteresis / three-point scaling / vaporisation limits (modelled; written back by the exporter)
                if (kw == "SATOPTS") {
                    bool other = false;
                    for (const Record& r : k.records) for (const auto& it : r.items) {
                        const std::string o = upper(it.s);
                        if (o == "HYSTER") c.props.hyst.on = true;
                        else if (o == "DIRECT") satDirect = true;   // C34 directional relative permeability (KRNUMX/Y/Z)
                        else if (o == "IRREVERS") { satDirect = true; grouped(1, "SATOPTS IRREVERS: the X- / Y- / Z- tables are taken equal to X / Y / Z"); }
                        else if (!o.empty()) other = true;
                    }
                    if (other) { warn("SATOPTS: only HYSTER and DIRECT are modelled by the ResSim engine (SURFTENS is not)"); keep(k, "option keyword"); }
                    continue;
                }
                if (kw == "NOHYST") { c.props.hyst.on = false; noHyst = true; continue; }
                if (kw == "EHYSTR") {
                    if (!k.records.empty()) {
                        const Record& r = k.records[0];
                        c.props.hyst.curv = num(r, 0, 0.1);
                        c.props.hyst.model = inum(r, 1, 0);
                        c.props.hyst.killough = num(r, 2, 0.1);
                        c.props.hyst.what = r.has(4) ? upper(str(r, 4)) : "BOTH";
                        if (c.props.hyst.model < 0 || c.props.hyst.model > 3)
                            warn("EHYSTR model " + std::to_string(c.props.hyst.model) + " is not modelled - Carlson (model 0) is used");
                        if (c.props.hyst.what != "KR") warn("EHYSTR: capillary pressure hysteresis is not modelled - relative permeability hysteresis only");
                    }
                    ehystr = true;
                    continue;
                }
                if (kw == "SCALECRS") {
                    c.grid.scalecrs = true;
                    if (!k.records.empty() && k.records[0].has(0)) { const std::string o = upper(str(k.records[0], 0)); c.grid.scalecrs = !(o == "NO" || o == "N"); }
                    continue;
                }
                if (kw == "VAPPARS") {
                    if (!k.records.empty()) {
                        const double v1 = std::max(0.0, num(k.records[0], 0, 0)), v2 = std::max(0.0, num(k.records[0], 1, 0));
                        if (section == "SCHEDULE") { pendingVap1 = v1; pendingVap2 = v2; pendingVapAfterDrsdt = pendingDrsdt > -1.5; }
                        else { c.props.vap1 = v1; c.props.vap2 = v2; }
                    }
                    continue;
                }
                if (section == "RUNSPEC" && kw == "SOLVENT") {
                    keep(k, "C35 solvent option (ResSim black-oil engine)");
                    grouped(0, "SOLVENT: simulated by the ResSim black-oil engine (C35 four-component solvent model: G+S fully implicit, solvent fraction transported, immiscible SSFN relative permeabilities) and kept verbatim for OPM Flow / ECLIPSE");
                    continue;
                }
                if (section == "RUNSPEC" && (kw == "SURFACT" || kw == "SURFACTW" || kw == "ALKALINE" || kw == "FOAM")) {
                    keep(k, "C35 chemical EOR option (ResSim black-oil engine)");
                    grouped(0, "chemical EOR (SURFACT / ALKALINE / FOAM): simulated by the ResSim black-oil engine (C35 sequential-implicit chemical transport coupled to the flow) and kept verbatim for OPM Flow / ECLIPSE");
                    continue;
                }
                if (section == "RUNSPEC") {       // process options that change the physics: the own engine would give wrong answers
                    static const std::map<std::string, const char*> kPhysics = {
                        {"BRINE", "brine tracking (BRINE)"},
                        {"THERMAL", "thermal model (THERMAL)"}, {"TEMP", "temperature option (TEMP)"}, {"CO2STORE", "CO2 storage (CO2STORE)"},
                        {"H2STORE", "hydrogen storage (H2STORE)"}, {"MICP", "MICP"}, {"BIOFILM", "biofilm"}, {"COMPS", "compositional model (COMPS)"},
                        {"GASWAT", "gas-water compositional (GASWAT)"}, {"MECH", "geomechanics (MECH)"}, {"COAL", "coal-bed methane (COAL)"},
                        {"API", "API tracking (API)"}, {"MULTIN", "multiple realisations (MULTIN)"}};
                    auto it = kPhysics.find(kw);
                    if (it != kPhysics.end()) { unsup(std::string(it->second) + " is not modelled"); keep(k, "physics option"); continue; }
                }
                if (kw == "ENDSCALE") {
                    c.grid.endscale = true;
                    if (!k.records.empty() && !k.records[0].items.empty()) {
                        const Record& r = k.records[0];
                        auto opt = [&](size_t i) { return i < r.items.size() && !r.items[i].def ? upper(r.items[i].s) : std::string(); };
                        if (opt(0) == "DIRECT") warn("ENDSCALE: directional end-point scaling is not modelled by the ResSim engine - isotropic scaling is used");
                        if (opt(1) == "IRREVERS") warn("ENDSCALE: irreversible scaling is not modelled by the ResSim engine");
                        // table counts (NTENDP / NSENDP for ENPTVD / ENKRVD) or non-default options: keep the record
                        // verbatim so the exported deck sizes the pass-through depth tables like the original
                        bool extra = false;
                        for (size_t i = 0; i < r.items.size(); ++i) if (!r.items[i].def && !(i == 0 && opt(0) == "NODIR") && !(i == 1 && opt(1) == "REVERS")) extra = true;
                        if (extra) keep(k, "ENDSCALE table dimensions");
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
                // ---- local grid refinements (Cartesian CARFIN, solved in place)
                if (kw == "LGR") continue;    // RUNSPEC sizes (MAXLGR, MAXCLS, ...): written by the exporter
                if (kw == "RADFIN" || kw == "RADFIN4" || kw == "AMALGAM") {
                    keep(k, "local grid refinement");
                    unsup(std::string(kw == "AMALGAM" ? "LGR amalgamation (AMALGAM)" : "radial local grid refinement (") + (kw == "AMALGAM" ? "" : kw + ")"));
                    continue;
                }
                if (kw == "CARFIN") { carfin(k); continue; }
                if (kw == "REFINE") {
                    leaveLgr();
                    const std::string nm = k.records.empty() ? std::string() : str(k.records[0], 0);
                    const int l = c.grid.lgrIndex(nm);
                    if (l < 0) fail("REFINE: local grid '" + nm + "' is not defined by CARFIN");
                    enterLgr(l);
                    continue;
                }
                if (kw == "ENDFIN") { if (lgrCtx < 0) warn("ENDFIN without CARFIN / REFINE - ignored"); leaveLgr(); continue; }
                if (kw == "NXFIN" || kw == "NYFIN" || kw == "NZFIN" || kw == "HXFIN" || kw == "HYFIN" || kw == "HZFIN") {
                    if (lgrCtx < 0) { warn(kw + " outside CARFIN ... ENDFIN - ignored"); continue; }
                    finArray(k);
                    continue;
                }
                if (kw == "LGRLOCK" || kw == "LGRFREE" || kw == "LGROFF" || kw == "LGRON" || kw == "LGRCOPY") {
                    keep(k, "local grid refinement option");
                    grouped(0, "LGR solution options: the ResSim engine always solves local grids in place, fully implicitly with the global grid");
                    continue;
                }
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
                if (kw == "GRAVDR") { c.grid.gravdr = true; continue; }
                if (kw == "DZMTRX") {   // one matrix block height for every cell
                    double v = 0;
                    if (!k.values.empty()) v = k.values[0];
                    else if (!k.records.empty()) v = num(k.records[0], 0, 0.0);
                    if (std::isfinite(v) && v > 0) dzmtrx = v;
                    continue;
                }
                if (kw == "GRAVDRM" || kw == "GRAVDRB") {
                    keep(k, "dual porosity gravity drainage");
                    grouped(1, "alternative dual porosity gravity drainage models (GRAVDRM / GRAVDRB) are not modelled by the ResSim engine "
                               "(GRAVDR is) - kept for OPM Flow / ECLIPSE");
                    continue;
                }
                if (kw == "DPNUM") {    // one value per matrix cell (NX*NY*NZ/2); a full-grid array is also accepted
                    needDims();
                    dpnumVals = k.values;
                    continue;
                }
                if (kw == "RADIAL") fail("radial grids are not supported");
                // ---------------- GRID / EDIT / REGIONS / SOLUTION arrays
                if (kw == "XMF" || kw == "YMF") { keep(k, "initial phase composition (thermal / compositional)"); continue; }
                if (k.isArray) {
                    dataArray(k);
                    continue;
                }
                if (kw == "BOX") { needDims(); boxActive = false; if (!k.records.empty()) { box = readBox(k.records[0], 0, fullBox()); boxActive = true; } continue; }
                if (kw == "ENDBOX") { boxActive = false; continue; }
                if (gridOps.count(kw)) { operations(k); continue; }
                if (regOps.count(kw)) { regionOperations(k); continue; }
                if (kw == "OPERATE" || kw == "OPERATER") { operate(k); continue; }   // C34
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
                    if (section == "SCHEDULE") {   // C46: time-dependent, cumulative (applied by the engine from this report step)
                        for (const Record& r : k.records)
                            if (!r.items.empty()) c.schedMultflt.push_back({int(c.schedule.size()), str(r, 0), num(r, 1, 1.0)});
                        keep(k, "fault multipliers (SCHEDULE)");
                        grouped(0, "MULTFLT in SCHEDULE: cumulative fault transmissibility multipliers applied by the ResSim engine from this report step on (C46) and kept verbatim");
                        continue;
                    }
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
                        m.nncMode = upper(str(r, 4, "ALL"));
                        m.set = set.empty() || set[0] == 'M' ? "MULTNUM" : set[0] == 'F' ? "FLUXNUM" : "OPERNUM";
                        multregt.push_back(m);
                    }
                    continue;
                }
                if (kw == "ROCKCOMP" || kw == "ROCKOPTS" || kw == "ROCKTAB" || kw == "ROCKTABH" || kw == "ROCK2D" || kw == "ROCK2DTR" || kw == "ROCKWNOD" || kw == "OVERBURD") {
                    // C34 rock compaction: modelled by the ResSim black-oil engine; the keywords are kept verbatim
                    RockCompInput& R = c.props.rockComp;
                    if (kw == "ROCKCOMP" && !k.records.empty()) {
                        const Record& r = k.records[0];
                        const std::string h = upper(str(r, 0, "REVERS"));
                        R.irreversible = h == "IRREVERS" || h == "HYSTER";
                        R.hyster = h == "HYSTER";   // C41: ROCKTABH elastic reload curves (OPM FlowProblem::rockCompPoroMultiplier)
                        R.water = upper(str(r, 2, "NO")) == "YES";
                        R.on = h != "NO" || R.water;
                    } else if (kw == "ROCKOPTS" && !k.records.empty()) {
                        const Record& r = k.records[0];
                        R.store = upper(str(r, 1, "NOSTORE")) == "STORE";
                        const std::string t = upper(str(r, 2, "PVTNUM"));
                        if (t == "SATNUM" || t == "ROCKNUM" || t == "PVTNUM") R.regionArray = t;
                    } else if (kw == "ROCKTAB") {
                        for (const Record& r : k.records) {
                            std::vector<double> v = numbers(r);
                            RockCompInput::Tab1 t;
                            if (v.empty() && !R.tab.empty()) t = R.tab.back();
                            const size_t nc = v.size() % 5 == 0 && v.size() % 3 != 0 ? 5 : 3;   // DIRECT: tm in X, Y, Z (X used)
                            for (size_t m = 0; m + nc <= v.size(); m += nc) { t.p.push_back(v[m]); t.pv.push_back(v[m + 1]); t.tm.push_back(v[m + 2]); }
                            R.tab.push_back(t);
                        }
                    } else if (kw == "ROCKTABH") {   // C41: one record per elastic curve, a lone '/' ends a region's table
                        R.tabH.assign(1, {});
                        for (const Record& r : k.records) {
                            if (r.terminator) { if (!R.tabH.back().empty()) R.tabH.emplace_back(); continue; }
                            std::vector<double> v = numbers(r);
                            const size_t nc = v.size() % 5 == 0 && v.size() % 3 != 0 ? 5 : 3;
                            RockCompInput::Tab1 t;
                            for (size_t m = 0; m + nc <= v.size(); m += nc) { t.p.push_back(v[m]); t.pv.push_back(v[m + 1]); t.tm.push_back(v[m + 2]); }
                            if (!t.p.empty()) R.tabH.back().push_back(t);
                        }
                        if (R.tabH.back().empty()) R.tabH.pop_back();
                        R.tab.clear();
                        for (const auto& curves : R.tabH) {
                            RockCompInput::Tab1 d;
                            for (const auto& cv : curves) { d.p.push_back(cv.p[0]); d.pv.push_back(cv.pv[0]); d.tm.push_back(cv.tm[0]); }
                            R.tab.push_back(d);
                        }
                    } else if (kw == "ROCK2D" || kw == "ROCK2DTR") {
                        size_t reg = 0; bool started = false;
                        for (const Record& r : k.records) {
                            std::vector<double> v = numbers(r);
                            if (v.empty()) { if (started) { ++reg; started = false; } continue; }
                            if (R.tab2.size() <= reg) R.tab2.resize(reg + 1);
                            RockCompInput::Tab2& t = R.tab2[reg];
                            std::vector<double> vals(v.begin() + 1, v.end());
                            if (kw == "ROCK2D") { t.p.push_back(v[0]); t.pv.push_back(vals); }
                            else t.tm.push_back(vals);
                            started = true;
                        }
                    } else if (kw == "ROCKWNOD") {
                        for (const Record& r : k.records) { std::vector<double> v = numbers(r); if (!v.empty()) R.swNodes.push_back(v); }
                    } else if (kw == "OVERBURD") {
                        for (const Record& r : k.records) { std::vector<double> v = numbers(r); if (v.empty() && !R.overburd.empty()) v = R.overburd.back(); R.overburd.push_back(v); }
                    }
                    keep(k, "rock compaction: modelled by the ResSim engine and kept verbatim");
                    continue;
                }
                if (kw == "AQUTAB") {   // C34 Carter-Tracy influence tables (table 1 = built-in)
                    for (const Record& r : k.records) {
                        std::vector<double> v = numbers(r);
                        if (v.empty() && !aqutab.empty()) v = aqutab.back();
                        aqutab.push_back(v);
                    }
                    keep(k, "Carter-Tracy influence table: modelled by the ResSim engine and kept verbatim");
                    continue;
                }
                if (kw == "AQUNUM" || kw == "AQUCON") {   // C34 numerical aquifers (ResSim engine; kept verbatim)
                    for (const Record& r : k.records) {
                        if (r.items.empty()) continue;
                        if (kw == "AQUNUM") {
                            NumAquCell a; a.id = inum(r, 0, 1); a.i = inum(r, 1, 1) - 1; a.j = inum(r, 2, 1) - 1; a.k = inum(r, 3, 1) - 1;
                            a.area = num(r, 4, 0); a.length = num(r, 5, 0); a.poro = num(r, 6, 0); a.perm = num(r, 7, 0); a.depth = num(r, 8, 0);
                            a.p0 = num(r, 9, -1); a.pvt = inum(r, 10, 0); a.sat = inum(r, 11, 0);
                            c.grid.numAquCells.push_back(a);
                        } else {
                            NumAquCon a; a.id = inum(r, 0, 1);
                            a.i1 = inum(r, 1, 1) - 1; a.i2 = inum(r, 2, 1) - 1; a.j1 = inum(r, 3, 1) - 1; a.j2 = inum(r, 4, 1) - 1; a.k1 = inum(r, 5, 1) - 1; a.k2 = inum(r, 6, 1) - 1;
                            std::string f = upper(str(r, 7, "I+"));
                            for (auto& ch : f) if (ch == '\'' || ch == ' ') ch = 0;
                            f = f.c_str();
                            a.face = !f.empty() && (f[0] == 'J' || f[0] == 'Y') ? 'Y' : !f.empty() && (f[0] == 'K' || f[0] == 'Z') ? 'Z' : 'X';
                            a.minus = f.find('-') != std::string::npos;
                            a.mult = num(r, 8, 1.0); a.transOpt = inum(r, 9, 0);
                            const std::string ai = upper(str(r, 10, "NO")); a.allowInternal = !ai.empty() && (ai[0] == 'Y' || ai[0] == 'T' || ai == "1");
                            c.grid.numAquCons.push_back(a);
                        }
                    }
                    keep(k, "numerical aquifer: modelled by the ResSim engine and kept verbatim");
                    continue;
                }
                if (kw == "EDITNNC" || kw == "EDITNNCR") {   // C34
                    for (const Record& r : k.records) {
                        if (r.items.empty()) continue;
                        NNCInput e;
                        e.i1 = inum(r, 0, 1) - 1; e.j1 = inum(r, 1, 1) - 1; e.k1 = inum(r, 2, 1) - 1;
                        e.i2 = inum(r, 3, 1) - 1; e.j2 = inum(r, 4, 1) - 1; e.k2 = inum(r, 5, 1) - 1; e.trans = num(r, 6, kw == "EDITNNC" ? 1.0 : 0.0);
                        (kw == "EDITNNC" ? c.grid.editnnc : c.grid.editnncr).push_back(e);
                    }
                    keep(k, kw == "EDITNNC" ? "NNC transmissibility multipliers: modelled by the ResSim engine and kept verbatim"
                                            : "NNC transmissibility replacement: modelled by the ResSim engine and kept verbatim");
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
                if (kw == "GUIDERAT") {   // C34: used by the black-oil group allocation and the compositional core; kept verbatim
                    if (!k.records.empty()) {
                        const Record& r = k.records[0];
                        GuideRateSpec g; g.set = true;
                        g.minInterval = num(r, 0, 0); g.phase = upper(str(r, 1, "NONE"));
                        g.A = num(r, 2, 0); g.B = num(r, 3, 0); g.C = num(r, 4, 0); g.D = num(r, 5, 0); g.E = num(r, 6, 0); g.F = num(r, 7, 0);
                        g.allowIncrease = upper(str(r, 8, "YES")) != "NO"; g.damping = num(r, 9, 1.0);
                        pendingGuide = g;
                    }
                    keep(k, "guide-rate formula: modelled by the ResSim engine (group target shares) and kept verbatim");
                    continue;
                }
                if (kw == "GRUPTREE" || kw == "GCONPROD" || kw == "GCONINJE" || kw == "GCONPRI" || kw == "VFPPROD" || kw == "VFPINJ") {
                    // simulated by the ResSim engine; the keyword is also kept verbatim so the exported deck is unchanged
                    auto ctl = [&](const std::string& m) {
                        return m == "ORAT" ? Control::ORAT : m == "WRAT" ? Control::WRAT : m == "GRAT" ? Control::GRAT :
                               m == "LRAT" ? Control::LRAT : m == "RESV" ? Control::RESV : Control::BHP;
                    };
                    // C48: a group rate item may hold a UDQ name (UDA, FU_ / GU_): recorded, the value follows the UDQ in the run
                    auto gnum = [&](const Record& r, size_t i, std::vector<std::pair<std::string, std::string>>& uda, const char* target,
                                    const std::string& gname) {
                        if (r.has(i) && !r.isNumber(i) && !r.items[i].def) {
                            const std::string u = upper(r.items[i].s);
                            if (u.size() > 2 && (u[0] == 'F' || u[0] == 'G') && u[1] == 'U') {
                                uda.push_back({target, u});
                                grouped(0, "UDA: the group target follows a user-defined quantity (UDQ) during the run", gname + " " + target + "=" + u);
                                return -1.0;
                            }
                        }
                        return num(r, i, -1);
                    };
                    if (kw == "GRUPTREE") {
                        for (const Record& r : k.records) if (r.has(0)) c.groupTree.push_back({str(r, 0), str(r, 1, "FIELD")});
                    } else if (kw == "GCONPROD") {
                        for (const Record& r : k.records) {
                            if (!r.has(0)) continue;
                            GroupProdControl g; g.group = str(r, 0);
                            const std::string m = upper(str(r, 1, "NONE"));
                            g.mode = ctl(m);
                            if (m != "NONE" && m != "FLD" && g.mode == Control::BHP) grouped(1, "GCONPROD control mode " + m + " is not modelled - limits only", g.group);
                            g.orat = gnum(r, 2, g.uda, "ORAT", g.group); g.wrat = gnum(r, 3, g.uda, "WRAT", g.group);
                            g.grat = gnum(r, 4, g.uda, "GRAT", g.group); g.lrat = gnum(r, 5, g.uda, "LRAT", g.group);
                            g.action = upper(str(r, 6, "NONE")); g.resv = gnum(r, 13, g.uda, "RESV", g.group);
                            g.actW = upper(str(r, 10, "")); g.actG = upper(str(r, 11, "")); g.actL = upper(str(r, 12, ""));   // C35 per-phase procedures
                            g.guideRate = num(r, 8, -1); g.guideDef = upper(str(r, 9, ""));   // C34 group guide rate
                            for (const std::string* a : {&g.action, &g.actW, &g.actG, &g.actL})
                                if (!a->empty() && *a != "NONE" && *a != "RATE") grouped(1, "GCONPROD procedure " + *a + " on exceeding a limit is not modelled (RATE is used)", g.group);
                            pendingGprod.push_back(g);
                        }
                    } else if (kw == "GCONPRI") {
                        for (const Record& r : k.records) {
                            if (!r.has(0)) continue;
                            GroupProdControl g; g.group = str(r, 0); g.priority = true; g.action = "RATE";
                            g.orat = num(r, 1, -1); g.wrat = num(r, 3, -1); g.grat = num(r, 5, -1); g.lrat = num(r, 7, -1); g.resv = num(r, 9, -1);
                            pendingGprod.push_back(g);
                        }
                        grouped(1, "GCONPRI: priority ranking is approximated by well potentials (group limits honoured)");
                    } else if (kw == "GCONINJE") {
                        for (const Record& r : k.records) {
                            if (!r.has(0)) continue;
                            GroupInjControl g; g.group = str(r, 0);
                            const std::string ph = upper(str(r, 1, "WATER"));
                            if (ph != "WATER" && ph != "WAT" && ph != "GAS") { grouped(1, "GCONINJE phase " + ph + " is not modelled", g.group); continue; }
                            g.phase = ph == "GAS" ? 2 : 1;
                            g.mode = upper(str(r, 2, "NONE"));
                            g.rate = gnum(r, 3, g.uda, "RATE", g.group); g.resv = gnum(r, 4, g.uda, "RESV", g.group);
                            g.rein = gnum(r, 5, g.uda, "REIN", g.group); g.vrep = gnum(r, 6, g.uda, "VREP", g.group);
                            g.reinGroup = str(r, 10, ""); g.vrepGroup = str(r, 11, "");
                            if (g.mode == "FLD") g.mode = "NONE";
                            pendingGinj.push_back(g);
                        }
                    } else {   // VFPPROD / VFPINJ
                        const bool inj = kw == "VFPINJ";
                        if (k.records.size() < (inj ? 3u : 6u)) fail("incomplete table");
                        VfpTable t; t.inj = inj;
                        const Record& h = k.records[0];
                        t.id = inum(h, 0, 0); t.datum = num(h, 1, 0);
                        t.flo = upper(str(h, 2, inj ? "WAT" : "OIL"));
                        size_t ri = 1;
                        auto axis = [&](std::vector<double>& v) {
                            const Record& r = k.records[ri++];
                            for (size_t i = 0; i < r.items.size(); ++i) v.push_back(num(r, i, 0));
                        };
                        if (!inj) { t.wfr = upper(str(h, 3, "WCT")); t.gfr = upper(str(h, 4, "GOR")); t.alq = upper(str(h, 6, "")); }
                        axis(t.floV); axis(t.thpV);
                        if (!inj) { axis(t.wfrV); axis(t.gfrV); axis(t.alqV); }
                        const size_t nf = t.floV.size(), nt = t.thpV.size(), nw = inj ? 1 : t.wfrV.size(), ng = inj ? 1 : t.gfrV.size(), na = inj ? 1 : t.alqV.size();
                        if (!nf || !nt || !nw || !ng || !na) fail("empty axis");
                        t.bhp.assign(nf * nt * nw * ng * na, std::nan(""));
                        for (; ri < k.records.size(); ++ri) {
                            const Record& r = k.records[ri];
                            if (r.items.empty()) continue;
                            const size_t ni = inj ? 1 : 4;
                            const int it = inum(r, 0, 1) - 1, iw = inj ? 0 : inum(r, 1, 1) - 1, ig = inj ? 0 : inum(r, 2, 1) - 1, ia = inj ? 0 : inum(r, 3, 1) - 1;
                            if (it < 0 || it >= int(nt) || iw < 0 || iw >= int(nw) || ig < 0 || ig >= int(ng) || ia < 0 || ia >= int(na)) fail("table index out of range");
                            const size_t base = ((((size_t(it) * nw + size_t(iw)) * ng + size_t(ig)) * na + size_t(ia)) * nf);
                            for (size_t f = 0; f < nf && ni + f < r.items.size(); ++f) t.bhp[base + f] = num(r, ni + f, 0);
                        }
                        size_t missing = 0;
                        for (double v : t.bhp) if (std::isnan(v)) ++missing;
                        if (missing) warn(std::to_string(missing) + " BHP values missing in " + kw + " table " + std::to_string(t.id) + " - treated as 1e10 (unreachable)");
                        for (double& v : t.bhp) if (std::isnan(v)) v = 1e10;
                        bool replaced = false;
                        for (auto& o : c.vfp) if (o.id == t.id && o.inj == t.inj) { o = t; replaced = true; }
                        if (!replaced) c.vfp.push_back(std::move(t));
                    }
                    keep(k, "group / lift data (simulated, kept verbatim)");
                    continue;
                }
                if (kw == "WELPI") {   // C34: connection factors rescaled to the target PI (ResSim engine); kept verbatim
                    for (const Record& r : k.records) if (r.has(0) && r.has(1)) pendingWelpi.push_back({str(r, 0), num(r, 1, 0)});
                    keep(k, "well productivity index: modelled by the ResSim engine and kept verbatim");
                    continue;
                }
                if (kw == "FBHPDEF") {   // C44: default BHP limits of the well controls that follow (ECLIPSE / OPM)
                    if (!k.records.empty() && !k.records[0].items.empty()) {
                        const Record& r = k.records[0];
                        fbhpProd = r.has(0) ? num(r, 0, -1) : -1; fbhpInj = r.has(1) ? num(r, 1, -1) : -1;
                    }
                    keep(k, "default BHP limits");
                    grouped(0, "FBHPDEF: default producer / injector BHP limits applied by the ResSim engine to the well controls that follow (C44) and kept verbatim");
                    continue;
                }
                if (kw == "COMPLUMP")   // C46: completion numbers for WELOPEN items 6-7 (the keyword itself stays verbatim below)
                    for (const Record& r : k.records)
                        if (r.has(0) && r.has(5)) complump.push_back({str(r, 0), {inum(r, 1, 0), inum(r, 2, 0), inum(r, 3, 0), inum(r, 4, 0), inum(r, 5, 0)}});
                if (schKws.count(kw)) { schedule(k); continue; }
                if (kw == "WLIST") { wlistRecords(k); keep(k, "well lists"); grouped(0, "WLIST: well lists are expanded by the ResSim engine (C33) and kept verbatim"); continue; }
                if (kw == "WCYCLE") { keep(k, "well cycling"); grouped(0, "WCYCLE: automatic well opening / closing cycles applied by the ResSim engine (C37: on / off periods, linear start-up efficiency ramp, time steps landing on the switches) and kept verbatim"); continue; }
                if (kw == "WEFAC") { keep(k, "well efficiency factors"); grouped(0, "WEFAC: efficiency factors applied by the ResSim engine (C33: group / field rates, cumulatives and the reservoir source use factor x rate) and kept verbatim"); continue; }
                if (kw == "UDT") { keep(k, "user-defined tables"); grouped(0, "UDT: user-defined tables evaluated in UDQ expressions by the ResSim engine (C38: TU_NAME[x], 1-D, NV / LC / LL interpolation) and kept verbatim"); continue; }
                if (kw == "UDQ") { keep(k, "user-defined quantities"); grouped(0, "UDQ: evaluated by the ResSim engine (C31: DEFINE / ASSIGN / UPDATE / UNITS, UDA targets) and kept verbatim"); continue; }
                if (kw == "ACTIONX") { keep(k, "action block"); grouped(0, "ACTIONX: evaluated by the ResSim engine at the end of every report step (C31) and kept verbatim"); continue; }
                if (kw == "LIFTOPT" || kw == "WLIFTOPT" || kw == "GLIFTOPT") {
                    keep(k, "gas-lift optimisation");
                    grouped(0, "gas-lift optimisation (LIFTOPT / WLIFTOPT / GLIFTOPT): evaluated by the ResSim engine (C35: lift gas allocated in increments by "
                               "the weighted incremental gradient, group lift gas / total gas limits) and kept verbatim");
                    continue;
                }
                if (kw == "GSATPROD" || kw == "GCONSUMP" || kw == "GCONSALE") {
                    keep(k, "facilities");
                    grouped(0, "satellite production / gas consumption / sales (GSATPROD / GCONSUMP / GCONSALE): modelled by the ResSim engine (C35: satellite rates in the group / "
                               "field rates and group control, gas consumption / import -> FGCR / FGIMR / FGSR, taken out of the network flow and the "
                               "re-injected gas, sales target met by the group gas injection; sales limits / procedures are not acted on) and kept verbatim");
                    continue;
                }
                if (kw == "PYACTION") {
                    keep(k, "action block");
                    grouped(0, "PYACTION (Python actions): run by the ResSim engine (C37) in a Python 3 child process with an opm_embedded shim (summary state, "
                               "schedule open / shut / insert_keywords, legacy run()) at the end of every report step - OPM Flow calls them every time step; kept verbatim");
                    continue;
                }
                if (kw == "WELSEGS" || kw == "COMPSEGS" || kw == "WSEGVALV" || kw == "WSEGSICD" || kw == "WSEGAICD" || kw == "WSEGLABY" || kw == "WSEGDFPA" || kw == "WSEGDFMD" ||
                    kw == "WSEGFLIM" || kw == "WSEGPULL") {
                    keep(k, "multisegment well");
                    grouped(0, "multisegment well data: evaluated by the ResSim engine (segment network: hydrostatic + friction + acceleration 'HFA', drift flux 'DF' "
                               "(Shi et al. 2005 gas-liquid + oil-water, WSEGDFPA; WSEGDFMD: ORIGINAL / H-K only), tubing flash, segment storage, valves / SICD / AICD / "
                               "labyrinth, flow-limiting valves WSEGFLIM, pumps WSEGPULL) and kept verbatim");
                    continue;
                }
                if (kw == "WSEGDFIN" || kw == "WSEGTABL") {
                    keep(k, "multisegment well");
                    grouped(1, kw + ": this segment device / model is not evaluated by the ResSim engine (segment treated without it) - kept for OPM Flow");
                    continue;
                }
                {   // C35 chemical EOR data: kept verbatim, interpreted by the black-oil engine (chem.cpp)
                    static const std::set<std::string> chemKw = {"SURFVISC", "SURFST", "SURFCAPD", "SURFADS", "SURFROCK", "WSURFACT", "ALKADS", "ALKROCK",
                        "ALSURFST", "ALSURFAD", "ALPOLADS", "WALKALIN", "FOAMOPTS", "FOAMMOB", "FOAMMOBP", "FOAMMOBS", "FOAMADS", "FOAMROCK",
                        "FOAMDCYW", "FOAMDCYO", "WFOAM", "SURFOPTS",
                        "PVDS", "SDENSITY", "SSFN", "WSOLVENT", "MISC", "PMISC", "MSFN", "SORWMIS", "SGCWMIS", "TLMIXPAR", "TLPMIXPA"};   // + C35 SOLVENT
                    if (chemKw.count(kw)) { keep(k, "C35 chemical EOR data (ResSim black-oil engine)"); continue; }
                }
                if (compositionalKeywords().count(kw)) { keep(k, "compositional data (ResSim compositional core)"); continue; }
                if (thermalKeywords().count(kw)) { keep(k, "thermal data (ResSim thermal core)"); continue; }
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
                const bool located = w.find(" (") != std::string::npos && w.compare(0, kw.size(), kw) == 0;
                // only the grid dimensions are indispensable; anything else is kept verbatim so the original deck
                // still runs in OPM Flow / ECLIPSE (and the rest of the model is imported for viewing)
                if (kw == "DIMENS" || kw == "SPECGRID") { if (located) throw; fail(w); }
                std::string why = w;
                if (located) { const size_t q = w.find("): "); if (q != std::string::npos) why = w.substr(q + 3); }
                keep(k, "could not be interpreted by ResSim");
                unsup(kw + ": " + why);
            }
        }
        K = nullptr;
        for (auto& [name, cf] : unsupported) {
            K = cf.second;
            if (name == "UDQ") { msg(0, "UDQ: evaluated by the ResSim engine (C31: DEFINE / ASSIGN / UPDATE / UNITS, UDA targets) and kept verbatim"); continue; }
            if (name == "CECON") { msg(0, "CECON: connection economic limits (WCT / GOR / WGR; CON, +CON, WELL workovers) modelled by the ResSim engine (C41) and kept verbatim"); continue; }
            if (name == "COMPLUMP") { msg(0, "COMPLUMP: completion lumps used by the WECON connection workovers (C34) and kept verbatim"); continue; }
            std::string t = (name == "WECON" || name == "WTEST") ? "economic limits / well tests (C27; C34 CON / +CON / PLUG workovers, secondary water cut, POTN): modelled by the ResSim engine and kept verbatim for OPM Flow"
                          : cf.second->known ? "not modelled by the ResSim engine - kept verbatim for OPM Flow (Advanced keywords)"
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
        if (T < sgwfn.size() && R.swof.empty() && R.sgof.empty()) {   // SGWFN -> SWOF (Sw = 1 - Sg, krw) + SGOF (Sg, krg, Pcgw)
            const auto& v = sgwfn[T];
            for (size_t i = v.size(); i >= 4; i -= 4) R.swof.push_back({1.0 - v[i - 4], v[i - 2], 0.0, 0.0});
            for (size_t i = 0; i + 4 <= v.size(); i += 4) R.sgof.push_back({v[i], v[i + 1], 0.0, v[i + 3]});
        }
        if (R.sgof.empty() && T < slgof.size()) {
            const auto& v = slgof[T];
            for (size_t i = v.size(); i >= 4; i -= 4) R.sgof.push_back({1.0 - v[i - 4], v[i - 3], v[i - 2], v[i - 1]});
        }
        const bool fam2 = T < swfn.size() || T < sgfn.size();
        if (fam2) {
            auto kroOf = [&](double so, bool gasOil) -> double {
                if (T < sof3.size()) { const auto x = col(sof3[T], 3, 0); return lin(x, col(sof3[T], 3, gasOil ? 2 : 1), so); }
                if (T < sof2.size()) { const auto x = col(sof2[T], 2, 0); return lin(x, col(sof2[T], 2, 1), so); }
                if (!c.phases.oil) return 0.0;   // gas-water system (CO2STORE): no oil phase
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
        leaveLgr();
        for (const auto& [wn, v] : ressimPiMult)   // C44: the exported WPIMULT records carry piMult x mult of every open period
            for (auto& w : c.wells) if (w.name == wn) { w.piMult = v; for (auto& q : w.completions) if (q.lgr < 0) q.mult /= v; }
        if (n == 0) throw std::runtime_error("Deck: DIMENS is missing");
        for (size_t l = 0; l < c.grid.lgrs.size() && l < lgrArr.size(); ++l) {
            LgrInput& L = c.grid.lgrs[l];
            for (auto& kv : lgrArr[l]) {
                std::string nm = kv.first;
                if (nm == "TRAN-EDIT") continue;
                bool any = false; for (double v : kv.second) if (!std::isnan(v)) any = true;
                if (any) L.props[nm] = kv.second;
            }
            std::string why;
            if (!L.valid(&why)) { if (L.parent == "GLOBAL") unsup(why); }
        }
        if (!c.grid.lgrs.empty()) {
            int cells = 0; for (auto& L : c.grid.lgrs) cells += L.size();
            msg(0, std::to_string(c.grid.lgrs.size()) + " local grid refinement(s) (CARFIN), " + std::to_string(cells) + " local cells, solved in place");
            if (c.grid.dual) unsup("local grid refinements in a dual porosity model");
        }
        if (!c.phases.oil && !c.phases.water) { unsup("gas-only system (no OIL or WATER phase)"); c.phases.water = true; }
        if (c.phases.gas && c.phases.water && !c.phases.oil) unsup("gas-water system without oil (SGWFN / SWFN+SGFN)");
        {   // C49: OPM's SOLVENT + CO2SOL / H2SOL model (solvent properties from the CO2 / H2 tables, no PVDS) is not modelled
            bool solv = false, sol = false, pvds = false;
            for (const auto& k : D.keywords) { solv |= k.name == "SOLVENT"; sol |= k.name == "CO2SOL" || k.name == "H2SOL"; pvds |= k.name == "PVDS"; }
            if (solv && sol && !pvds) unsup("SOLVENT with CO2SOL / H2SOL (solvent from CO2 / H2 property tables, no PVDS)");
        }
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
                const double ph = nm[1] == 'Z' ? 10.0 : 100.0;
                if (!take(nm, dst)) { unsup(std::string("grid geometry: ") + nm + " is missing (binary GDFILE / unsupported grid input) - placeholder cells shown"); dst.assign(size_t(n), ph); }
                if (undefined(dst)) { unsup(std::string("grid geometry: ") + nm + " is not defined for every cell - placeholder values shown"); for (double& x : dst) if (std::isnan(x)) x = ph; }
            }
            if (arr.find("TOPS") == arr.end()) { unsup("grid geometry: TOPS is missing - placeholder depths shown"); arr["TOPS"].assign(size_t(n), NaN); for (int q = 0; q < nl; ++q) arr["TOPS"][size_t(q)] = 1000.0; }
            auto it = arr.find("TOPS");
            std::vector<double>& t = it->second;
            for (int q = 0; q < nl; ++q) if (std::isnan(t[size_t(q)])) { unsup("grid geometry: TOPS is missing in the top layer - placeholder depths shown"); t[size_t(q)] = 1000.0; }
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
        if (!take("PORO", g.poro)) { unsup("PORO is missing"); g.poro.assign(size_t(n), 0.2); }
        if (!take("PERMX", g.permx)) { unsup("PERMX is missing"); g.permx.assign(size_t(n), 100.0); }
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
        take("MULTTRX", g.tranMultX); take("MULTTRY", g.tranMultY); take("MULTTRZ", g.tranMultZ);   // C38 EDIT TRANX/Y/Z
        for (auto* v : {&g.ntg, &g.multx, &g.multy, &g.multz, &g.multpv, &g.tranMultX, &g.tranMultY, &g.tranMultZ}) for (double& x : *v) if (std::isnan(x)) x = 1.0;
        for (auto* v : {&g.tranMultX, &g.tranMultY, &g.tranMultZ})
            if (std::all_of(v->begin(), v->end(), [](double x) { return x == 1.0; })) v->clear();
        {   // C33 MULTX- / MULTY- / MULTZ-: the multiplier of the face to the minus neighbour (i-1, j-1, k-1) is folded into
            // MULTX / MULTY / MULTZ of that neighbour (ECLIPSE combines both multiplicatively on the same face; exact for
            // neighbour connections, fault NNCs use the multiplier of the low-index cell). Export writes the folded arrays.
            const int NX = c.grid.nx, NY = c.grid.ny, NZ = c.grid.nz;
            const char* nm[3] = {"MULTX-", "MULTY-", "MULTZ-"};
            std::vector<double>* dst[3] = {&g.multx, &g.multy, &g.multz};
            int folded = 0;
            for (int d = 0; d < 3; ++d) {
                std::vector<double> m;
                if (!take(nm[d], m) || int(m.size()) != NX * NY * NZ) continue;
                std::vector<double>& D = *dst[d];
                if (D.size() != m.size()) D.assign(m.size(), 1.0);
                for (int k = 0; k < NZ; ++k)
                    for (int j = 0; j < NY; ++j)
                        for (int i = 0; i < NX; ++i) {
                            const size_t u = size_t(i + NX * (j + NY * k));
                            const double f = std::isnan(m[u]) ? 1.0 : m[u];
                            if (f == 1.0) continue;
                            const int ii = i - (d == 0), jj = j - (d == 1), kk = k - (d == 2);
                            if (ii < 0 || jj < 0 || kk < 0) continue;
                            D[size_t(ii + NX * (jj + NY * kk))] *= f; ++folded;
                        }
            }
            if (folded) msg(0, "MULTX- / MULTY- / MULTZ-: " + std::to_string(folded) + " face multipliers folded into MULTX / MULTY / MULTZ of the minus neighbour");
        }

        // ---- dual porosity
        take("SIGMAV", g.sigmav); take("LX", g.lx); take("LY", g.ly); take("LZ", g.lz); take("MULTMF", g.multmf);
        take("DZMTRXV", g.dzmtrxv);
        if (take("DISPERC", g.disperc)) for (double& x : g.disperc) if (std::isnan(x) || x < 0) x = 0.0;
        if (g.dzmtrxv.empty() && dzmtrx > 0) g.dzmtrxv.assign(size_t(n), dzmtrx);
        g.dpnum.clear();
        if (!dpnumVals.empty()) {
            if (!g.dual) msg(1, "DPNUM given without DUALPORO / DUALPERM - ignored");
            else {
                const size_t half = size_t(n) / 2;
                if (dpnumVals.size() != half && dpnumVals.size() != size_t(n))
                    msg(1, "DPNUM: " + std::to_string(dpnumVals.size()) + " values, expected NX*NY*NZ/2 = " + std::to_string(half) + " - ignored");
                else {
                    g.dpnum.resize(half);
                    int nsp = 0;
                    for (size_t q = 0; q < half; ++q) {
                        const double v = dpnumVals[q];
                        g.dpnum[q] = std::isfinite(v) ? int(std::lround(v)) : 1;
                        nsp += g.dpnum[q] == 0;
                    }
                    msg(0, "DPNUM: " + std::to_string(nsp) + " single porosity blocks (their fracture cells carry the matrix rock)");
                    if (!nsp) g.dpnum.clear();
                }
            }
        }
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
            if (g.gravdr) msg(0, std::string("GRAVDR: gravity drainage between matrix and fracture, block height ") +
                                 (g.dzmtrxv.empty() ? "= cell thickness (no DZMTRX / DZMTRXV)" : "from DZMTRX / DZMTRXV"));
            if (g.minc > 1) {
                if (g.dual != 1) msg(2, "MINC needs DUALPORO (dual permeability with nested matrix shells is not supported) - the run will be rejected");
                else msg(0, "MINC: " + std::to_string(g.minc) + " nested matrix shells per block (" + std::to_string(g.mincDims) +
                            "-D fracture sets), transient matrix-fracture exchange; RESSIM own engine only (ECLIPSE / OPM Flow run the Warren-Root dual porosity)");
            }
        } else if (g.sigma > 0 || !g.sigmav.empty()) msg(1, "SIGMA/SIGMAV given without DUALPORO - ignored");
        if (!g.dual && g.minc > 1) msg(2, "MINC given without DUALPORO - add DUALPORO and the fracture layers");
        if (!g.hfracs.empty())
            msg(0, std::to_string(g.hfracs.size()) + " hydraulic fracture(s) (RESSIM HFRAC, embedded discrete fractures): RESSIM own engine only - "
                   "ECLIPSE / OPM Flow ignore the comment lines (export a refined / high-permeability equivalent for them)");
        {
            int ng = 0; for (const auto& F : g.hfracs) if (F.kic > 0) ++ng;
            if (ng) msg(c.comp.enabled ? 1 : 0, std::to_string(ng) + " fracture(s) with in-run LEFM growth (RESSIM HFGROW): black-oil engine only" +
                                                 (c.comp.enabled ? std::string(" - the compositional engine keeps the whole trace open") : std::string()));
        }
        if (g.langVL > 0 && g.langPL > 0) msg(0, "Langmuir adsorbed gas (RESSIM LANGMUIR): RESSIM own engine only");
        if (g.permGamma > 0) msg(0, "stress-dependent matrix permeability (RESSIM PERMSTRESS): RESSIM own engine only");

        // ---- MULTREGT: kept as records + the region set; the grid applies them to every connection (incl. NNCs)
        if (!multregt.empty()) {
            for (const MultRegT& m : multregt) {
                usedArrays.insert(m.set);
                GridInput::MultRegT x; x.r1 = m.r1; x.r2 = m.r2; x.m = m.m; x.set = m.set[0];
                x.dirs.clear(); for (char ch : m.dirs) if (ch == 'X' || ch == 'Y' || ch == 'Z') x.dirs += ch;
                if (x.dirs.empty()) x.dirs = "XYZ";
                x.nnc = m.nncMode == "NNC" ? 1 : m.nncMode == "NONNC" ? 2 : m.nncMode == "NOAQUNNC" ? 3 : 0;
                g.multregt.push_back(x);
                std::vector<int>& dst = m.set[0] == 'M' ? g.multnum : m.set[0] == 'F' ? g.fluxnum : g.opernum;
                if (dst.empty()) { const std::vector<double>& R = array(m.set); dst.resize(size_t(n)); for (int q = 0; q < n; ++q) dst[size_t(q)] = int(std::lround(R[size_t(q)])); }
            }
            msg(0, "MULTREGT: " + std::to_string(multregt.size()) + " region-pair multipliers (applied to neighbour and non-neighbour connections)");
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
        // C34 numerical aquifer cells: AQUNUM PVT / saturation table numbers
        for (const NumAquCell& a : g.numAquCells) {
            if (a.i < 0 || a.i >= g.nx || a.j < 0 || a.j >= g.ny || a.k < 0 || a.k >= g.nz) { msg(2, "AQUNUM cell outside the grid"); continue; }
            const size_t q = size_t(g.index(a.i, a.j, a.k));
            if (a.pvt > 0) { if (g.pvtnum.empty()) g.pvtnum.assign(size_t(n), 1); g.pvtnum[q] = a.pvt; }
            if (a.sat > 0) { if (g.satnum.empty()) g.satnum.assign(size_t(n), 1); g.satnum[q] = a.sat; }
        }
        if (!g.numAquCells.empty()) msg(0, "numerical aquifers (AQUNUM / AQUCON): " + std::to_string(g.numAquCells.size()) + " aquifer cells, " + std::to_string(g.numAquCons.size()) + " connection records");
        intArray("IMBNUM", g.imbnum, ntsfun, "saturation tables");
        if (!c.props.hyst.on) g.imbnum.clear();
        if (arr.count("ROCKNUM")) { intArray("ROCKNUM", g.rocknum, 0, ""); c.props.rockComp.regionArray = "ROCKNUM"; }   // C34 (an explicit ROCKNUM always wins, OPM)
        if (c.props.rockComp.active()) msg(0, std::string("rock compaction (ROCKCOMP): ") + (c.props.rockComp.tab2.empty() ? "ROCKTAB" : "ROCK2D / ROCK2DTR water-induced") +
                                              (c.props.rockComp.irreversible ? ", irreversible" : ", reversible") + ", regions from " + c.props.rockComp.regionArray +
                                              (c.props.rockComp.overburd.empty() ? "" : ", OVERBURD effective pressure"));
        // C34 directional relative permeability tables (used with SATOPTS DIRECT, as in ECLIPSE / OPM)
        intArray("KRNUMX", g.krnumx, ntsfun, "saturation tables"); intArray("KRNUMY", g.krnumy, ntsfun, "saturation tables"); intArray("KRNUMZ", g.krnumz, ntsfun, "saturation tables");
        intArray("SURFNUM", g.surfnum, ntsfun, "saturation tables");   // C35 miscible surfactant curves
        if (g.surfnum.empty() && arr.count("SURFNUM")) g.surfnum.assign(size_t(n), 1);   // SURFNUM 1 everywhere is still a table choice
        if (arr.count("SFOAM")) { usedArrays.insert("SFOAM"); g.sfoam = arr.at("SFOAM"); for (double& x : g.sfoam) if (std::isnan(x)) x = 0; }   // C35
        if (arr.count("SSOL")) { usedArrays.insert("SSOL"); g.ssol = arr.at("SSOL"); for (double& x : g.ssol) if (std::isnan(x)) x = 0; }      // C35 SOLVENT
        if (!satDirect && (!g.krnumx.empty() || !g.krnumy.empty() || !g.krnumz.empty())) { msg(1, "KRNUMX/Y/Z given without SATOPTS DIRECT - ignored (as in ECLIPSE)"); g.krnumx.clear(); g.krnumy.clear(); g.krnumz.clear(); }
        int npvtUsed = 1, nsatUsed = 1, neqlUsed = 1;
        for (int v : g.pvtnum) npvtUsed = std::max(npvtUsed, v);
        for (int v : g.satnum) nsatUsed = std::max(nsatUsed, v);
        for (int v : g.imbnum) nsatUsed = std::max(nsatUsed, std::min(v, ntsfun));   // imbibition tables (hysteresis)
        for (const auto* kr : {&g.krnumx, &g.krnumy, &g.krnumz, &g.surfnum}) for (int v : *kr) nsatUsed = std::max(nsatUsed, std::min(v, ntsfun));
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
            const bool needO = c.phases.oil && !c.phases.disgas && pvt[u].pvdo.empty(), needL = c.phases.disgas && pvt[u].pvto.empty();
            const bool needG = c.phases.gas && pvt[u].pvdg.empty();
            if (needO || needL || needG) {
                bool compDeck = false;   // E300 decks: the EOS fluid replaces the black-oil tables (placeholders for display)
                for (const auto& r : c.engineUnsupported) if (r.rfind("compositional model", 0) == 0 || r.rfind("CO2 storage", 0) == 0) compDeck = true;
                if (!compDeck) unsup(std::string(needO ? "PVDO (or PVCDO)" : needL ? "PVTO" : "PVDG (or PVTG)") + " is missing or not supported (e.g. PVZG, RSCONST)" + reg);
                Case tmp; placeholderFluids(tmp, c.units == UnitSystem::Metric);
                if (needO) { for (auto& r : tmp.props.pvto) pvt[u].pvdo.push_back({r.pb, r.bo, r.mu}); }
                if (needL) pvt[u].pvto = tmp.props.pvto;
                if (needG) pvt[u].pvdg = tmp.props.pvdg;
            }
            if (c.phases.vapoil && pvt[u].pvtg.size() < 2) { unsup("PVTG with at least 2 pressures is needed (VAPOIL)" + reg); c.phases.vapoil = false; }
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
            if ((c.phases.water && c.phases.oil && R.swof.size() < 2) || (c.phases.gas && c.phases.oil && R.sgof.size() < 2)) {
                unsup(std::string(R.swof.size() < 2 && c.phases.water ? "SWOF (or SWFN + SOF3/SOF2)" : "SGOF (or SGFN/SLGOF)") + " is missing or not supported" + reg);
                if (R.swof.size() < 2) R.swof = {{0.2, 0, 1, 0}, {0.5, 0.15, 0.25, 0}, {0.8, 0.5, 0, 0}, {1.0, 1.0, 0, 0}};
                if (R.sgof.size() < 2) R.sgof = {{0, 0, 1, 0}, {0.3, 0.2, 0.15, 0}, {0.8, 1.0, 0, 0}};
            }
            sats.push_back(R);
        }
        if (!swfn.empty() || !sgfn.empty() || !slgof.empty()) msg(0, "saturation functions converted to SWOF/SGOF form");
        P.swof = sats[0].swof; P.sgof = sats[0].sgof;
        P.satExtra.assign(sats.begin() + 1, sats.end());
        if (ntpvt > npvtUsed) msg(0, std::to_string(ntpvt - npvtUsed) + " PVT table(s) not referenced by PVTNUM - dropped");
        if (ntsfun > nsatUsed) msg(0, std::to_string(ntsfun - nsatUsed) + " saturation table(s) not referenced by SATNUM / IMBNUM / KRNUM / SURFNUM - dropped");

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
                if (m) unsup("PRESSURE is not defined for " + std::to_string(m) + " active cells");
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
            unsup("no initial conditions ResSim can use (EQUIL or PRESSURE/SWAT) in SOLUTION");
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
            // C34: OPM Flow applies end-point arrays even without ENDSCALE in RUNSPEC (opm-tests model2 6_ENDSCALE); ECLIPSE
            // needs ENDSCALE - follow OPM and say so (the export then writes ENDSCALE so both simulators agree)
            bool anyEps = false;
            for (int i = 0; i < 8; ++i) if (arr.count(epsNames[i])) anyEps = true;
            if (anyEps && !c.grid.endscale) {
                c.grid.endscale = true;
                msg(1, "end-point scaling arrays given without ENDSCALE in RUNSPEC - applied as OPM Flow does (ECLIPSE ignores them "
                       "without ENDSCALE; the exported deck adds ENDSCALE)");
            }
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
            if (c.grid.scalecrs && c.grid.endscale) msg(0, "SCALECRS: three-point relative permeability scaling");
            // hysteresis: imbibition regions and imbibition end points
            static const char* iNames[] = {"ISWL", "ISWCR", "ISWU", "ISGL", "ISGCR", "ISGU", "ISOWCR", "ISOGCR"};
            std::vector<double>* idst[] = {&c.grid.iswl, &c.grid.iswcr, &c.grid.iswu, &c.grid.isgl, &c.grid.isgcr, &c.grid.isgu, &c.grid.isowcr, &c.grid.isogcr};
            std::string ifound;
            for (int i = 0; i < 8; ++i) {
                if (!arr.count(iNames[i])) continue;
                usedArrays.insert(iNames[i]);
                ifound += std::string(ifound.empty() ? "" : ", ") + iNames[i];
                if (!c.grid.endscale || !c.props.hyst.on) continue;
                *idst[i] = arr.at(iNames[i]);
                for (double& x : *idst[i]) if (std::isnan(x)) x = -1;
            }
            if (c.props.hyst.on && !ehystr) { msg(1, "SATOPTS HYSTER without EHYSTR - default EHYSTR (Carlson, model 0) used"); }
            if (!c.props.hyst.on && ehystr && !noHyst) msg(1, "EHYSTR given without SATOPTS HYSTER - hysteresis is off (as in ECLIPSE)");
            if (c.props.hyst.on) {
                const char* mdl[] = {"Carlson non-wetting / drainage wetting", "Carlson non-wetting / imbibition wetting",
                                     "Killough non-wetting / drainage wetting", "Killough non-wetting / imbibition wetting"};
                msg(0, std::string("relative permeability hysteresis: ") + mdl[std::clamp(c.props.hyst.model, 0, 3)] +
                           (ifound.empty() ? std::string() : " - imbibition end points " + ifound));
            }
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
            if (aqTabOf.count(a.id)) {
                const int t = aqTabOf[a.id] - 2;
                if (t >= 0 && t < int(aqutab.size()) && aqutab[size_t(t)].size() >= 4) {
                    const auto& v = aqutab[size_t(t)];
                    a.tdTab.clear(); a.pdTab.clear();
                    for (size_t m = 0; m + 1 < v.size(); m += 2) { a.tdTab.push_back(v[m]); a.pdTab.push_back(v[m + 1]); }
                } else msg(1, "AQUCT aquifer " + std::to_string(a.id) + ": influence table " + std::to_string(aqTabOf[a.id]) + " not found in AQUTAB - built-in table used");
            }
            if (a.conns.empty()) msg(1, "aquifer " + std::to_string(a.id) + " has no AQUANCON connections - ignored");
            else keep.push_back(a);
        }
        c.aquifers = keep;

        // ---- E300 THERMAL rock heat arrays and initial temperature (global cells; used when THERMAL is present)
        for (const char* nm : {"HEATCR", "HEATCRT", "THCONR", "THCONSF", "TEMPI"}) {
            auto it = arr.find(nm);
            if (it == arr.end()) continue;
            std::vector<double>& dst = std::string(nm) == "HEATCR" ? c.thermal.heatcr : std::string(nm) == "HEATCRT" ? c.thermal.heatcrt
                                     : std::string(nm) == "THCONR" ? c.thermal.thconr : std::string(nm) == "THCONSF" ? c.thermal.thconsf : c.thermal.tempi;
            dst = it->second;
            usedArrays.insert(nm);
        }
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

        if (c.schedule.empty()) {
            unsup("SCHEDULE has no TSTEP or DATES ResSim can read");
            ReportStep st; st.length = 1.0; c.schedule.push_back(st);
        }
        // C33: the black-oil engine follows the type of the last control keyword at run time (Simulator::applyControls,
        // like OPM Flow / ECLIPSE): a WCONHIST with zero rates on an injector is a *stopped* producer (open, BHP reported),
        // an injector converted to a producer later keeps both periods. The thermal / compositional front ends still
        // have one type per well: there the periods as the other type are simulated as shut.
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
                    if (!thermalDeck) continue;
                    wc.status = WellStatus::Shut;
                    if (!prodWell) { wc.mode = Control::RATE; wc.rate = 0; } else { wc.mode = Control::BHP; }
                }
            for (auto& [name, cnt] : converted)
                msg(0, "well " + name + " changes between producer and injector (" + std::to_string(cnt) + " control records)" +
                       std::string(thermalDeck ? " - the periods as the other type are simulated as shut" : " - the type follows the control keyword at run time"));
        }
        std::string err = thermalDeck ? std::string() : validateCase(c);
        if (!err.empty() && co2storeDeck() && !c.engineUnsupported.empty() &&
            std::all_of(c.engineUnsupported.begin(), c.engineUnsupported.end(), [](const std::string& r) { return r.rfind("CO2 storage", 0) == 0 || r.rfind("gas-water system", 0) == 0; }))
            err.clear();   // C47: CO2STORE decks run on the CO2STORE core - not a validation problem
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

// E300 compositional decks run in ECLIPSE 300 / OPM Flow (pass-through); the EOS fluid is still read and checked with
// the own thermodynamics library so the import report shows the fluid and its saturation pressure.
void compositionalReport(const deck::RawDeck& D, bool field, std::vector<DeckMessage>& msgs) {
    try {
        eos::CompFluid fl; std::vector<std::string> notes;
        if (!eos::fluidFromDeck(D, field, fl, &notes)) return;
        auto add = [&](int level, const std::string& t) { DeckMessage m; m.level = level; m.keyword = "EOS"; m.text = t; msgs.push_back(m); };
        for (const auto& n : notes) add(1, n);
        const eos::Model& m = fl.eosRegions.front();
        std::string names; for (const auto& n : fl.names) names += (names.empty() ? "" : " ") + n;
        add(0, std::to_string(m.nc()) + "-component " + (m.type == eos::EosType::PR ? (fl.pr78 ? "Peng-Robinson (1978)" : "Peng-Robinson") : "Soave-Redlich-Kwong") +
                   " fluid: " + names + (fl.co2store ? " (CO2STORE)" : ""));
        std::vector<double> z = fl.zi;
        if (z.empty() && !fl.zmfvd.empty() && !fl.zmfvd[0].z.empty()) z = fl.zmfvd[0].z.front();
        bool okData = true;
        for (const auto& c : m.comp) if (!(c.tc > 0 && c.pc > 0 && c.mw > 0)) okData = false;
        if (fl.resT > 0 && int(z.size()) == m.nc() && okData) {
            const double ps = eos::saturationPressure(m, fl.resT, z, !eos::isVapourLike(m, z, fl.resT), 0);
            const bool vap = eos::isVapourLike(m, z, fl.resT);
            char b[200];
            if (ps > 0)
                std::snprintf(b, sizeof b, "initial fluid at reservoir temperature: %s, saturation pressure %.1f %s (own EOS flash)",
                              vap ? "gas (dew point)" : "oil (bubble point)", field ? eos::paToPsia(ps) : ps / 1e5, field ? "psia" : "bar");
            else std::snprintf(b, sizeof b, "initial fluid at reservoir temperature: single phase (%s-like) at all pressures checked", vap ? "gas" : "oil");
            add(0, b);
        }
    } catch (...) {}
}

}  // namespace

// Compositional / CO2STORE data of a deck -> Case::comp (fluid, WELLSTRE streams, FIELDSEP, WINJGAS per report step).
// The case then runs on the ResSim compositional core; the compositional keywords stay verbatim for export.
bool compositionalFromDeck(const deck::RawDeck& D, Case& c, std::vector<DeckMessage>* msgs) {
    const bool field = c.units == UnitSystem::Field;
    eos::CompFluid fl; std::vector<std::string> notes;
    try { if (!eos::fluidFromDeck(D, field, fl, &notes) || fl.eosRegions.empty()) return false; } catch (...) { return false; }
    CompositionalInput& C = c.comp;
    C = CompositionalInput();
    C.enabled = true; C.co2store = fl.co2store;
    if (C.co2store) {   // dedicated CO2-brine model: always CO2 + H2O
        const double T = fl.resT, m = fl.salinity;
        eos::CompFluid f = co2StoreFluid(T, m);
        f.tempvd = fl.tempvd; f.stdT = fl.stdT; f.stdP = fl.stdP;
        fl = f;
    }
    C.fluid = fl;
    const int nc = C.nc();
    auto tempToK = [&](double t) { return field ? (t + 459.67) * 5.0 / 9.0 : t + 273.15; };
    const double pf = field ? 6894.757293168 : 1e5;
    std::map<std::string, size_t> streamIdx;
    std::string section; int step = 0;
    const double fDiff = (field ? 0.09290304 : 1.0) / 86400.0;   // ft2/day | m2/day -> m2/s
    bool diffuse = false;
    static const char* kSections[] = {"RUNSPEC", "GRID", "EDIT", "PROPS", "REGIONS", "SOLUTION", "SUMMARY", "SCHEDULE"};
    auto upc = [](std::string x) { for (char& ch : x) ch = char(std::toupper(static_cast<unsigned char>(ch))); return x; };
    for (const deck::Keyword& k : D.keywords) {
        bool isSec = false;
        for (const char* q : kSections) if (k.name == q) { section = q; isSec = true; }
        if (isSec) continue;
        if (k.name == "FIELDSEP" && step == 0) {
            C.separator.clear();
            std::map<int, std::pair<double, double>> sep;
            for (const deck::Record& r : k.records) {
                if (r.terminator || r.items.empty()) continue;
                sep[r.integer(0, 1)] = {tempToK(r.num(1, field ? 60.0 : 15.56)), r.num(2, field ? 14.696 : 1.01325) * pf};
            }
            for (auto& kv : sep) C.separator.push_back(kv.second);
        } else if (k.name == "DIFFUSE") {
            C.trap.diffusion = true; diffuse = true;
        } else if ((k.name == "DIFFCWAT" || k.name == "DIFFCGAS") && !k.records.empty() && !k.records[0].terminator) {
            const deck::Record& r = k.records[0];
            double& a = k.name == "DIFFCWAT" ? C.trap.dCo2Water : C.trap.dCo2Gas;
            double& b = k.name == "DIFFCWAT" ? C.trap.dH2oWater : C.trap.dH2oGas;
            a = r.num(0, a / fDiff) * fDiff; b = r.num(1, b / fDiff) * fDiff;
        } else if (k.name == "DIFFC" && !k.records.empty() && !k.records[0].terminator) {
            const deck::Record& r = k.records[0];
            const double gg = r.num(2, 0.0), og = r.num(3, 0.0), go = r.num(4, 0.0), oo = r.num(5, 0.0);
            const double dG = (gg > 0 ? gg : og) * fDiff, dO = (go > 0 ? go : oo) * fDiff;
            if (dO > 0) C.diffOil.assign(size_t(nc), dO);
            if (dG > 0) C.diffGas.assign(size_t(nc), dG);
        } else if (k.name == "PRECSALT") {
            C.trap.saltPrecip = true;
        } else if (k.name == "VAPWAT") {
            C.waterVap = true;   // EOS mode: water vapour in the gas phase (CO2STORE: vaporised water for salt dry-out)
        } else if (k.name == "SALTSOL" && !k.records.empty() && !k.records[0].terminator) {
            const double cs = k.records[0].num(0, 6.15 * 58.44) * (field ? 2.853010 : 1.0);   // kg/sm3
            C.trap.saltSolubility = cs / 58.44;
        } else if (k.name == "PERMFACT") {
            C.trap.permFact.clear();
            for (const deck::Record& r : k.records) {
                if (r.terminator) break;
                for (size_t i = 0; i + 1 < r.items.size(); i += 2) C.trap.permFact.push_back({r.num(i, 1.0), r.num(i + 1, 1.0)});
                if (!C.trap.permFact.empty()) break;   // table 1 (all SATNUM regions use it)
            }
        } else if (k.name == "SEPCOND") {
            std::map<std::string, std::map<int, std::pair<double, double>>> tr;
            for (const deck::Record& r : k.records) {
                if (r.terminator || r.items.empty()) continue;
                tr[upc(r.str(0))][r.integer(2, 1)] = {tempToK(r.num(3, field ? 60.0 : 15.56)), r.num(4, field ? 14.696 : 1.01325) * pf};
            }
            for (auto& kv : tr) {
                CompositionalInput::SepTrain* T = nullptr;
                for (auto& x : C.sepTrains) if (x.name == kv.first) T = &x;
                if (!T) { C.sepTrains.push_back({kv.first, {}}); T = &C.sepTrains.back(); }
                T->stages.clear();
                for (auto& st : kv.second) T->stages.push_back(st.second);
            }
        } else if (k.name == "WSEPCOND") {
            for (const deck::Record& r : k.records) if (!r.terminator && !r.items.empty()) C.wellSep[r.str(0)] = upc(r.str(1));
            if (step > 0) notes.push_back("WSEPCOND after the first report step: applied from the start of the run");
        } else if (k.name == "__RESSIM") {
            for (const deck::Record& r : k.records) {
                if (r.items.empty()) continue;
                const std::string key = upc(r.str(0));
                if (key == "MINERAL") {
                    co2::MineralReaction m; m.name = r.str(1, m.name);
                    m.logk25 = r.num(2, m.logk25); m.Ea = r.num(3, m.Ea * 1e-3) * 1e3; m.area = r.num(4, m.area);
                    m.nuCO2 = r.num(5, m.nuCO2); m.nuH2O = r.num(6, m.nuH2O); m.dVsolid = r.num(7, m.dVsolid * 1e6) * 1e-6;
                    m.reactant0 = r.num(8, m.reactant0); m.vmReactant = r.num(9, m.vmReactant * 1e6) * 1e-6; m.mEq = r.num(10, m.mEq);
                    C.trap.minerals.push_back(m);
                } else if (key == "WATERVAP") C.waterVap = true;
                else if (key == "SGTMAX") C.trap.sgtMax = r.num(1, -1);
                else if (key == "STRUCTTOL") C.trap.structTol = r.num(1, C.trap.structTol);
                else if (key == "VERMAPRUESS") { C.trap.vpPhiR = r.num(1, C.trap.vpPhiR); C.trap.vpGamma = r.num(2, C.trap.vpGamma); }
                else if (key == "TORTUOSITY") C.trap.tortuosity = r.integer(1, 0);
                else if (key == "RHOSALT") C.trap.rhoSalt = r.num(1, C.trap.rhoSalt);
                else if (key == "DIFFCOMP") {
                    std::vector<double>& D = upc(r.str(1)) == "GAS" ? C.diffGas : C.diffOil;
                    D.clear(); for (size_t i = 2; i < r.items.size(); ++i) D.push_back(r.num(i, 0.0) * fDiff);
                    D.resize(size_t(nc), 0.0);
                }
            }
        } else if (k.name == "WELLSTRE") {
            for (const deck::Record& r : k.records) {
                if (r.terminator || r.items.empty()) continue;
                CompStream st; st.name = upc(r.str(0));
                for (size_t i = 1; i < r.items.size(); ++i) st.z.push_back(r.num(i, 0.0));
                st.z.resize(size_t(nc), 0.0);
                auto it = streamIdx.find(st.name);
                if (it == streamIdx.end()) { streamIdx[st.name] = C.streams.size(); C.streams.push_back(st); }
                else C.streams[it->second] = st;
            }
        } else if (k.name == "WINJGAS") {
            for (const deck::Record& r : k.records) {
                if (r.terminator || r.items.empty()) continue;
                const std::string how = upc(r.str(1, "STREAM"));
                if ((how == "STREAM" || how == "WELLSTRE") && size_t(step) < c.schedule.size())
                    c.schedule[size_t(step)].injStreams.push_back({r.str(0), upc(r.str(2))});
            }
        } else if (k.name == "TSTEP" && section == "SCHEDULE") {
            int nv = 0; for (const deck::Record& r : k.records) nv += int(r.items.size());
            step += nv;
        } else if (k.name == "DATES" && section == "SCHEDULE") {
            for (const deck::Record& r : k.records) if (!r.terminator && !r.items.empty()) ++step;
        }
    }
    if (!C.co2store) {
        if (!diffuse) { C.diffOil.clear(); C.diffGas.clear(); }
        else if (C.diffOil.empty() && C.diffGas.empty()) notes.push_back("DIFFUSE without DIFFC: no diffusion coefficients - diffusion is off");
        C.trap.diffusion = false;
        if (C.trap.saltPrecip) notes.push_back("PRECSALT: salt precipitation is modelled in CO2STORE runs only");
    }
    if (!C.co2store && !C.trap.minerals.empty()) notes.push_back("mineral trapping is modelled in CO2STORE runs only");
    if (!C.co2store && C.fluid.co2sol) notes.push_back("CO2SOL: CO2 dissolves in the water (ResSim compositional core / ECLIPSE 300); OPM Flow does not support CO2SOL");
    // the case is now fully represented: remove the "black-oil engine cannot model it" markers
    auto& U = c.engineUnsupported;
    U.erase(std::remove_if(U.begin(), U.end(), [&](const std::string& r) {
        return r.rfind("compositional model", 0) == 0 || (C.co2store && (r.rfind("CO2 storage", 0) == 0 || r.rfind("gas-water system", 0) == 0));
    }), U.end());
    if (msgs) {
        DeckMessage m; m.level = 0; m.keyword = C.co2store ? "CO2STORE" : "COMPS";
        char b[300];
        if (C.co2store)
            std::snprintf(b, sizeof b, "CO2 storage model (CO2 + brine, Spycher-Pruess / Duan-Sun): T %.1f C, salinity %.3g mol/kg - runs on the ResSim CO2STORE core",
                          fl.resT - 273.15, fl.salinity);
        else std::snprintf(b, sizeof b, "%d-component EOS model, %zu injection stream(s), %zu separator stage(s) - runs on the ResSim compositional core", nc, C.streams.size(), C.separator.size());
        m.text = b; msgs->push_back(m);
        for (const auto& n : notes) { DeckMessage w; w.level = 1; w.keyword = m.keyword; w.text = n; msgs->push_back(w); }
    }
    return true;
}

namespace {
Case importRaw(const deck::RawDeck& D, std::vector<std::string>* warnings, std::vector<DeckMessage>* report) {
    Importer im(D);
    try {
        bool field = false;   // ECLIPSE default unit system is METRIC
        for (const auto& k : D.keywords) if (k.name == "FIELD") field = true;
        bool thermal = false;
        for (const auto& k : D.keywords) if (k.name == "THERMAL") thermal = true;
        if (!thermal) compositionalReport(D, field, im.msgs);
        im.run();
        im.finish();
        if (thermal) thermalFromDeck(D, im.c, &im.msgs);
        else compositionalFromDeck(D, im.c, &im.msgs);
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

// ======================================================================== public: run-time schedule edits (C31 ACTIONX)
bool applyScheduleText(Case& cs, int step, const std::string& text, std::vector<std::string>* messages) {
    auto note = [&](const std::string& m) { if (messages) messages->push_back(m); };
    const deck::RawDeck D = deck::parseDeckText("SCHEDULE\n" + text, "", "ACTIONX");
    bool ok = true;
    for (const auto& m : D.messages) if (m.level == deck::Level::Error) { ok = false; note("parse: " + m.text); }
    const int S = std::max(0, std::min(step, int(cs.schedule.size())));
    if (S >= int(cs.schedule.size())) { note("no report step left - nothing applied"); return ok; }
    Importer im(D);
    im.c = std::move(cs);
    std::vector<ReportStep> tail(im.c.schedule.begin() + S, im.c.schedule.end());
    im.c.schedule.resize(size_t(S));
    im.anyStep = S > 0;
    im.n = im.c.grid.nx * im.c.grid.ny * im.c.grid.nz;
    im.section = "SCHEDULE";
    auto seed = [&](const WellControl& wc) { im.lastCtrl[wc.well] = wc; im.lastInj[wc.well] = wc.phase != 0 || wc.mode == Control::RATE; };
    for (int s = 0; s < S; ++s) for (const auto& wc : im.c.schedule[size_t(s)].controls) seed(wc);
    for (const auto& wc : tail[0].controls) seed(wc);
    const size_t nExtra = im.c.extraKeywords.size();
    // C33: well lists defined up to this step (kept as extra keywords) are rebuilt before the action's records
    for (const ExtraKeyword& e : im.c.extraKeywords) {
        if (e.name != "WLIST" || e.section != DeckSection::Schedule || e.step > S) continue;
        try {
            const deck::RawDeck W = deck::parseDeckText("SCHEDULE\nWLIST\n" + e.body, "", "WLIST");
            for (const auto& k : W.keywords) if (k.name == "WLIST") im.wlistRecords(k);
        } catch (...) {}
    }
    try { im.run(); }
    catch (const std::exception& e) { ok = false; note(e.what()); }
    for (const auto& m : im.msgs)
        if (m.level > 0 && m.keyword != "UDQ" && m.text.find("after the last TSTEP") == std::string::npos) note(m.keyword + ": " + m.text);   // (merged below)
    for (const std::string& key : im.groupOrder) {
        const auto& gr = im.groups[key];
        if (gr.level == 0 || gr.keyword == "UDQ") continue;   // UDQ records of an action are handled by the simulator
        std::string t = gr.keyword + ": " + key.substr(key.find('\x1f') + 1);
        for (size_t i = 0; i < gr.items.size() && i < 4; ++i) t += (i ? ", " : " - ") + gr.items[i];
        note(t);
    }
    std::vector<ExtraKeyword> runtime;   // C33: WEFAC from an action / script acts at run time (the simulator reads it)
    for (size_t i = nExtra; i < im.c.extraKeywords.size(); ++i) {
        ExtraKeyword& e = im.c.extraKeywords[i];
        if (e.name == "WEFAC" || e.name == "WCYCLE") { e.step = S; e.note = "runtime"; runtime.push_back(e); continue; }
        if (e.name != "UDQ") note(e.name + ": not applied by the ResSim engine at run time");
    }
    im.c.extraKeywords.resize(nExtra);
    for (auto& e : runtime) im.c.extraKeywords.push_back(std::move(e));
    // the action's records act from report step S on and win over the deck's records of that step
    ReportStep& t0 = tail[0];
    for (const auto& wc : im.pending) {
        bool done = false;
        for (auto& x : t0.controls) if (x.well == wc.well) { x = wc; done = true; }
        if (!done) t0.controls.push_back(wc);
    }
    for (const auto& g : im.pendingGprod) {
        bool done = false;
        for (auto& x : t0.gprod) if (x.group == g.group) { x = g; done = true; }
        if (!done) t0.gprod.push_back(g);
    }
    for (const auto& g : im.pendingGinj) {
        bool done = false;
        for (auto& x : t0.ginj) if (x.group == g.group && x.phase == g.phase) { x = g; done = true; }
        if (!done) t0.ginj.push_back(g);
    }
    for (const auto& p : im.pendingPoly) t0.wpolymer.push_back(p);
    if (im.pendingDrsdt > -1.5) { t0.drsdt = im.pendingDrsdt; t0.drsdtFree = im.pendingDrsdtFree; }
    for (auto& t : tail) im.c.schedule.push_back(std::move(t));
    cs = std::move(im.c);
    return ok;
}

// ======================================================================== public: keyword support
namespace deck {
Support ressimSupport(const std::string& keyword, std::string* note) {
    if (compositionalKeywords().count(keyword)) {
        if (note) *note = "E300 compositional data - used by the ResSim compositional core (EOS, FIM) and kept for OPM Flow / ECLIPSE";
        return Support::Full;
    }
    if (thermalKeywords().count(keyword) || keyword == "THERMAL") {
        if (note) *note = "E300 THERMAL data - used by the ResSim thermal core (K-values, IAPWS-IF97 steam, energy balance) and kept for ECLIPSE";
        return Support::Full;
    }
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
    Case c = importRaw(D, warnings, report);
    c.deckDir = baseDir;
    return c;
}

Case importEclipseDeck(const std::string& path, std::vector<std::string>* warnings, std::vector<DeckMessage>* report) {
    const deck::RawDeck D = deck::parseDeckFile(path);
    Case c = importRaw(D, warnings, report);
    {
        std::error_code ec;
        const std::filesystem::path pp = std::filesystem::absolute(std::filesystem::u8path(path), ec);
        c.deckDir = (ec ? std::filesystem::u8path(path) : pp).parent_path().u8string();
    }
    return c;
}

// ======================================================================== public: export
// '-- RESSIM:' lines of the trapping / interaction settings ECLIPSE has no keyword for (written in both export modes)
static std::string trapRessimLines(const Case& c) {
    const CompositionalInput& C = c.comp;
    if (!C.enabled) return "";
    const bool field = c.units == UnitSystem::Field;
    const double fDf = (field ? 0.09290304 : 1.0) / 86400.0;
    const co2::TrappingOptions& T = C.trap;
    std::ostringstream o;
    if (!C.co2store) {
        auto uniform = [](const std::vector<double>& v) { for (double x : v) if (std::fabs(x - v[0]) > 1e-15 * std::fabs(v[0])) return false; return true; };
        for (int ph = 0; ph < 2; ++ph) {
            const std::vector<double>& D = ph ? C.diffGas : C.diffOil;
            if (D.empty() || uniform(D)) continue;
            o << "-- RESSIM: DIFFCOMP " << (ph ? "GAS" : "OIL"); for (double v : D) o << ' ' << num(v / fDf); o << "\n";
        }
    }
    const co2::TrappingOptions D0;
    if (!C.co2store && C.waterVap) o << "-- RESSIM: WATERVAP\n";   // ECLIPSE E300 EOS mode has no VAPWAT
    if (T.sgtMax > 0) o << "-- RESSIM: SGTMAX " << num(T.sgtMax) << "\n";
    if (T.structTol != D0.structTol) o << "-- RESSIM: STRUCTTOL " << num(T.structTol) << "\n";
    if (T.tortuosity != 0) o << "-- RESSIM: TORTUOSITY " << T.tortuosity << "\n";
    if (C.co2store && (T.vpPhiR != D0.vpPhiR || T.vpGamma != D0.vpGamma)) o << "-- RESSIM: VERMAPRUESS " << num(T.vpPhiR) << ' ' << num(T.vpGamma) << "\n";
    if (C.co2store && T.rhoSalt != D0.rhoSalt) o << "-- RESSIM: RHOSALT " << num(T.rhoSalt) << "\n";
    if (C.co2store)
        for (const auto& m : T.minerals) {
        std::string nm = m.name; for (char& ch : nm) if (ch == '\'') ch = ' ';
        o << "-- RESSIM: MINERAL '" << nm << "' " << num(m.logk25) << ' ' << num(m.Ea * 1e-3) << ' ' << num(m.area) << ' ' << num(m.nuCO2) << ' '
          << num(m.nuH2O) << ' ' << num(m.dVsolid * 1e6) << ' ' << num(m.reactant0) << ' ' << num(m.vmReactant * 1e6) << ' ' << num(m.mEq) << "\n";
        }
    return o.str();
}

std::string compositionalKeywordsText(const Case& c, int section) {
    const CompositionalInput& C = c.comp;
    if (!C.enabled || C.fluid.eosRegions.empty()) return "";
    const bool field = c.units == UnitSystem::Field;
    const eos::CompFluid& F = C.fluid;
    const int nc = C.nc();
    std::ostringstream o;
    auto tOut = [&](double K) { return field ? K * 1.8 - 459.67 : K - 273.15; };
    const double fT = field ? 1.8 : 1.0, fP = field ? 6894.757293168 : 1e5, fV = field ? 0.028316846592 / 453.59237 : 1e-3, fD = field ? 0.3048 : 1.0;
    auto perComp = [&](const char* kw, auto get, bool always = true) {
        bool any = always;
        for (const auto& m : F.eosRegions) for (const auto& q : m.comp) if (get(q) != 0) any = true;
        if (!any) return;
        o << kw << "\n";
        for (const auto& m : F.eosRegions) { for (const auto& q : m.comp) o << ' ' << num(get(q)); o << " /\n"; }
        o << "\n";
    };
    const double fDf = (field ? 0.09290304 : 1.0) / 86400.0;
    const int npvt = c.props.numPvtRegions(), nsat = c.props.numSatRegions();
    const bool diffOn = C.co2store ? C.trap.diffusion : (!C.diffOil.empty() || !C.diffGas.empty());
    if (section == 0) {
        if (C.co2store) o << "CO2STORE\n\n";
        else o << "COMPS\n " << nc << " /\n\n";
        // CO2SOL (ECLIPSE 300): CO2 also dissolves in the water; without SOLUBILI / SOLUBILT ECLIPSE uses its default
        // (Chang-Coats-Nolen) data. OPM Flow does not support CO2SOL.
        if (!C.co2store && F.co2sol) o << "CO2SOL\n\n";
        if (diffOn) o << "DIFFUSE\n\n";
        if (C.co2store && C.trap.saltPrecip) o << "PRECSALT\n\nVAPWAT\n\n";
    } else if (section == 1) {
        if (!C.co2store && F.co2sol && F.salinity > 0) o << "SALINITY\n " << num(F.salinity) << " /    -- mol NaCl / kg H2O (CO2SOL)\n\n";
        if (!C.co2store && F.co2sol && F.solubili.size() >= 2) {
            const double frs = field ? 1000.0 / 5.614583333 : 1.0;
            o << "SOLUBILI\n-- P  RSW  BW  CW  VISC (CO2-saturated water)\n";
            for (const auto& r : F.solubili) o << ' ' << num(r.p / fP) << ' ' << num(r.rsw / frs) << ' ' << num(r.bw) << ' ' << num(r.cw * fP) << ' ' << num(r.muw * 1e3) << "\n";
            o << "/\n\n";
        }
        if (C.co2store) {
            o << "SALINITY\n " << num(F.salinity) << " /    -- mol NaCl / kg H2O\n\n";
        } else {
            o << "EOS\n"; for (const auto& m : F.eosRegions) o << ' ' << (m.type == eos::EosType::SRK ? "SRK" : "PR") << " /\n"; o << "\n";
            if (F.pr78) o << "PRCORR\n\n";
            o << "CNAMES\n"; for (const auto& n : F.names) o << ' ' << n; o << " /\n\n";
            perComp("TCRIT", [&](const eos::Component& q) { return q.tc * fT; });
            perComp("PCRIT", [&](const eos::Component& q) { return q.pc / fP; });
            perComp("VCRIT", [&](const eos::Component& q) { return q.vc / fV; });
            perComp("ZCRIT", [&](const eos::Component& q) { return q.zc; });
            perComp("ACF", [&](const eos::Component& q) { return q.omega; });
            perComp("MW", [&](const eos::Component& q) { return q.mw * 1e3; });
            perComp("SSHIFT", [&](const eos::Component& q) { return q.sshift; }, false);
            perComp("PARACHOR", [&](const eos::Component& q) { return q.parachor; }, false);
            perComp("VCRITVIS", [&](const eos::Component& q) { return q.vcVisc / fV; }, false);
            bool om = false; for (const auto& m : F.eosRegions) for (const auto& q : m.comp) om = om || q.omegaA >= 0 || q.omegaB >= 0;
            if (om) {
                const double dA = 0.457235529, dB = 0.077796074, sA = 0.42748023, sB = 0.08664035;
                o << "OMEGAA\n"; for (const auto& m : F.eosRegions) { for (const auto& q : m.comp) o << ' ' << num(q.omegaA >= 0 ? q.omegaA : m.type == eos::EosType::SRK ? sA : dA); o << " /\n"; } o << "\n";
                o << "OMEGAB\n"; for (const auto& m : F.eosRegions) { for (const auto& q : m.comp) o << ' ' << num(q.omegaB >= 0 ? q.omegaB : m.type == eos::EosType::SRK ? sB : dB); o << " /\n"; } o << "\n";
            }
            bool bic = false; for (const auto& m : F.eosRegions) bic = bic || !m.kij.empty();
            if (bic && nc > 1) {
                o << "BIC\n";
                for (const auto& m : F.eosRegions) {
                    for (int i = 1; i < nc; ++i) { o << ' '; for (int j = 0; j < i; ++j) o << ' ' << num(m.k(i, j)); o << "\n"; }
                    o << " /\n";
                }
                o << "\n";
            }
            const eos::Model def;
            const auto& m0 = F.eosRegions[0];
            bool lbcDef = true; for (int i = 0; i < 5; ++i) lbcDef = lbcDef && std::fabs(m0.lbc[i] - def.lbc[i]) < 1e-12;
            if (!lbcDef) { o << "LBCCOEF\n"; for (double v : m0.lbc) o << ' ' << num(v); o << " /\n\n"; }
            if (F.tempvd.size() > 1) { o << "TEMPVD\n"; for (auto& r : F.tempvd) o << ' ' << num(r.first / fD) << ' ' << num(tOut(r.second)) << "\n"; o << "/\n\n"; }
            else if (F.resT > 0) o << "RTEMP\n " << num(tOut(F.resT)) << " /\n\n";
            if (!F.zmfvd.empty()) {
                o << "ZMFVD\n";
                for (const auto& t : F.zmfvd) {
                    for (size_t r = 0; r < t.depth.size(); ++r) { o << ' ' << num(t.depth[r] / fD); for (double z : t.z[r]) o << ' ' << num(z); o << "\n"; }
                    o << "/\n";
                }
                o << "\n";
            }
        }
        o << "STCOND\n " << num(tOut(F.stdT)) << ' ' << num(F.stdP / fP) << " /\n\n";
        const co2::TrappingOptions& T = C.trap;
        if (C.co2store && diffOn) {
            o << "DIFFCWAT\n"; for (int r = 0; r < npvt; ++r) o << ' ' << num(T.dCo2Water / fDf) << ' ' << num(T.dH2oWater / fDf) << " /\n"; o << "\n";
            o << "DIFFCGAS\n"; for (int r = 0; r < npvt; ++r) o << ' ' << num(T.dCo2Gas / fDf) << ' ' << num(T.dH2oGas / fDf) << " /\n"; o << "\n";
        }
        if (!C.co2store && diffOn) {
            auto first = [](const std::vector<double>& v) { return v.empty() ? 0.0 : v[0]; };
            double mo = 0, mg = 0; const eos::Model& m0 = F.eosRegions[0];
            for (int i = 0; i < nc; ++i) { mo += m0.comp[size_t(i)].mw * 1e3 / nc; mg += m0.comp[size_t(i)].mw * 1e3 / nc; }
            o << "DIFFC\n";
            for (int r = 0; r < npvt; ++r)
                o << ' ' << num(mo) << ' ' << num(mg) << ' ' << num(first(C.diffGas) / fDf) << ' ' << num(first(C.diffGas) / fDf) << ' '
                  << num(first(C.diffOil) / fDf) << ' ' << num(first(C.diffOil) / fDf) << " /\n";
            o << "\n";
        }
        if (C.co2store && T.saltPrecip) {
            o << "SALTSOL\n"; for (int r = 0; r < npvt; ++r) o << ' ' << num(T.saltSolubility * 58.44 / (field ? 2.853010 : 1.0)) << " /\n"; o << "\n";
            if (!T.permFact.empty()) {
                o << "PERMFACT\n";
                for (int r = 0; r < nsat; ++r) { for (auto& pr : T.permFact) o << ' ' << num(pr.first) << ' ' << num(pr.second) << "\n"; o << "/\n"; }
                o << "\n";
            }
        }
        o << trapRessimLines(c);
        o << "\n";
    } else if (section == 2) {
        if (C.co2store) {
            o << "RTEMPVD\n";
            if (F.tempvd.size() > 1) for (auto& r : F.tempvd) o << ' ' << num(r.first / fD) << ' ' << num(tOut(r.second)) << "\n";
            else o << " 0 " << num(tOut(F.resT)) << "\n " << num(10000.0 / fD) << ' ' << num(tOut(F.resT)) << "\n";
            o << "/\n\n";
        } else {
            if (F.zmfvd.empty() && !F.zi.empty()) { o << "ZI\n"; for (double z : F.zi) o << ' ' << num(z); o << " /\n\n"; }
            if (!C.separator.empty()) {
                o << "FIELDSEP\n";
                for (size_t i = 0; i < C.separator.size(); ++i)
                    o << ' ' << i + 1 << ' ' << num(tOut(C.separator[i].first)) << ' ' << num(C.separator[i].second / fP) << " /\n";
                o << "/\n\n";
            }
        }
    } else if (section == 3) {
        if (!C.co2store && !C.streams.empty()) {
            o << "WELLSTRE\n";
            for (const auto& st : C.streams) { o << " '" << st.name << "'"; for (double z : st.z) o << ' ' << num(z); o << " /\n"; }
            o << "/\n\n";
        }
        if (!C.co2store && !C.sepTrains.empty()) {
            o << "SEPCOND\n";
            for (const auto& tr : C.sepTrains)
                for (size_t i = 0; i < tr.stages.size(); ++i)
                    o << " '" << tr.name << "' 'FIELD' " << i + 1 << ' ' << num(tOut(tr.stages[i].first)) << ' ' << num(tr.stages[i].second / fP) << " /\n";
            o << "/\n\n";
            if (!C.wellSep.empty()) {
                o << "WSEPCOND\n";
                for (const auto& kv : C.wellSep) o << " '" << kv.first << "' '" << kv.second << "' /\n";
                o << "/\n\n";
            }
        }
    }
    return o.str();
}

static std::string eclipseDeckToStringImpl(const Case& c);
std::string eclipseDeckToString(const Case& c) {
    if (c.thermal.enabled && c.thermal.structured) {   // structured thermal export: thermal keywords from c.thermal
        Case f = c;
        auto& X = f.extraKeywords;
        X.erase(std::remove_if(X.begin(), X.end(), [](const ExtraKeyword& e) {
            return thermalKeywords().count(e.name) > 0 || e.name == "WELLSTRE" || e.name == "WINJGAS";   // gas-injection streams: from c.thermal
        }), X.end());
        return eclipseDeckToStringImpl(f);
    }
    if (!(c.comp.enabled && c.comp.structured)) return eclipseDeckToStringImpl(c);
    // structured compositional export: the compositional keywords are written from c.comp, not from the verbatim copies
    Case f = c;
    auto& X = f.extraKeywords;
    X.erase(std::remove_if(X.begin(), X.end(), [](const ExtraKeyword& e) {
        return compositionalKeywords().count(e.name) || e.name == "CO2STORE" || e.name == "COMPS";
    }), X.end());
    return eclipseDeckToStringImpl(f);
}

static std::string eclipseDeckToStringImpl(const Case& c) {
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
            if (e.section == DeckSection::Schedule && e.note != "runtime" && (e.step == s || (last && e.step > s)) && (before < 0 || e.beforeControls == (before == 1)))
                out += extraText(e);
        return out;
    };
    const GridInput& g = c.grid;
    const PropsInput& p = c.props;
    const bool live = c.phases.gas && c.phases.oil && c.phases.disgas;
    const bool cp = g.cornerPoint();
    o << (c.thermal.enabled ? "-- ECLIPSE 300 THERMAL deck generated by ResSim Studio\n" : c.comp.enabled && !c.comp.co2store ? "-- ECLIPSE 300 (compositional) deck generated by ResSim Studio\n" : c.comp.co2store ? "-- ECLIPSE 300 / OPM Flow CO2STORE deck generated by ResSim Studio\n" : "-- ECLIPSE 100 deck generated by ResSim Studio\n");
    o << "-- Lines starting with '-- RESSIM:' carry settings ECLIPSE has no keyword for (ignored by ECLIPSE).\n\n";
    o << "RUNSPEC\n\nTITLE\n" << (c.title.empty() ? std::string("Untitled") : c.title) << "\n\n";
    o << "DIMENS\n " << g.nx << ' ' << g.ny << ' ' << g.nz << " /\n\n";
    if (c.phases.oil) o << "OIL\n";
    if (c.phases.water) o << "WATER\n";
    if (c.phases.gas) o << "GAS\n";
    if (live) o << "DISGAS\n";
    if (c.grid.nonnc) o << "NONNC\n";   // C38
    const bool vap = c.phases.gas && c.phases.oil && c.phases.vapoil;
    if (vap) o << "VAPOIL\n";
    if (c.grid.dual == 1) o << "DUALPORO\n";
    if (c.grid.dual == 2) o << "DUALPERM\n";
    if (c.grid.dual && c.grid.gravdr) o << "GRAVDR\n";
    if (!g.lgrs.empty()) {
        int maxc = 0; for (const auto& L : g.lgrs) maxc = std::max(maxc, L.size());
        o << "LGR\n " << g.lgrs.size() << ' ' << maxc << " /\n";
    }
    const bool polyOn = c.phases.polymer && c.phases.water;
    if (polyOn) o << "POLYMER\n";
    {
        bool keptEnd = false;
        for (const auto& x : c.extraKeywords) if (x.name == "ENDSCALE") keptEnd = true;
        if (g.hasEndPointScaling() && !keptEnd) o << "ENDSCALE\n NODIR REVERS /\n";
    }
    {
        const bool dirKr = !g.krnumx.empty() || !g.krnumy.empty() || !g.krnumz.empty();
        if (c.props.hyst.on || dirKr) o << "SATOPTS\n" << (c.props.hyst.on ? " 'HYSTER'" : "") << (dirKr ? " 'DIRECT'" : "") << " /\n";
    }
    const bool compS = c.comp.enabled && c.comp.structured, compOn = c.comp.enabled, co2 = compOn && c.comp.co2store;
    const bool thermOn = c.thermal.enabled;
    if (compS) o << compositionalKeywordsText(c, 0);
    if (thermOn) o << thermalKeywordsText(c, 0);
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
    {   // C38: UDTDIMS for the UDT tables (tables, rows, interpolation points, dimensions)
        size_t nt = 0, np = 0;
        for (const auto& e : c.extraKeywords) if (e.name == "UDT") {
            std::vector<udq::Table> T; std::string err; udq::parseTables(e.body, 0, T, err);
            nt += T.size(); for (const auto& t : T) np = std::max(np, t.x.size());
        }
        if (nt) o << "UDTDIMS\n " << nt << ' ' << std::max<size_t>(np, 2) << ' ' << std::max<size_t>(np, 2) << " 1 /\n\n";
    }
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
    arr("MULTNUM", g.multnum); arr("FLUXNUM", g.fluxnum); arr("OPERNUM", g.opernum);
    if (!g.disperc.empty()) arr("DISPERC", g.disperc);
    if (!g.multregt.empty()) {
        o << "MULTREGT\n";
        static const char* kNnc[4] = {"ALL", "NNC", "NONNC", "NOAQUNNC"};
        for (const auto& m : g.multregt) {
            auto reg = [](int r) { return r < 0 ? std::string("1*") : std::to_string(r); };
            o << ' ' << reg(m.r1) << ' ' << reg(m.r2) << ' ' << num(m.m) << ' ' << m.dirs << ' ' << kNnc[std::max(0, std::min(3, m.nnc))] << ' '
              << (m.set == 'F' ? "F" : m.set == 'O' ? "O" : "M") << " /\n";
        }
        o << "/\n\n";
    }
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
        if (g.minc > 1) {   // MINC has no ECLIPSE keyword: RESSIM extension line (dual porosity for ECLIPSE / OPM Flow)
            o << "-- RESSIM: MINC " << g.minc << ' ' << g.mincDims;
            for (double f : g.mincVol) o << ' ' << num(f);
            o << "\n\n";
        }
        if (!g.dzmtrxv.empty()) {
            bool uni = true; for (double v : g.dzmtrxv) if (v != g.dzmtrxv[0]) uni = false;
            if (uni) o << "DZMTRX\n " << num(g.dzmtrxv[0]) << " /\n\n"; else arr("DZMTRXV", g.dzmtrxv);
        }
        if (!g.dpnum.empty()) arr("DPNUM", g.dpnum);
        if (g.nodppm) o << "NODPPM\n\n";
    }
    // C15 hydraulic fractures / unconventional rock (ResSim extensions: ECLIPSE and OPM Flow ignore these comment lines)
    for (const GridInput::HFrac& F : g.hfracs)
        o << "-- RESSIM: HFRAC " << F.name << ' ' << num(F.x0) << ' ' << num(F.y0) << ' ' << num(F.x1) << ' ' << num(F.y1) << ' '
          << num(F.zTop) << ' ' << num(F.zBot) << ' ' << num(F.width) << ' ' << num(F.perm) << ' ' << num(F.poro) << ' '
          << num(F.gamma) << ' ' << F.satnum << "\n";
    for (const GridInput::HFrac& F : g.hfracs)
        if (F.kic > 0)
            o << "-- RESSIM: HFGROW " << F.name << ' ' << num(F.kic) << ' ' << num(F.sigma) << ' ' << num(F.sigmaGrad) << ' ' << num(F.initLen) << ' '
              << num(F.closedMult) << (F.subgrid ? (F.subLeak ? (F.subLeak == 2 ? " SUBGRID LEAKC" : " SUBGRID LEAK") : " SUBGRID") : "") << (F.feStress ? " FESTRESS" : "")
              << (F.ddmWidth ? (F.ddmE > 0 ? " DDMWIDTH " + num(F.ddmE) + " " + num(F.ddmNu) : std::string(" DDMWIDTH")) : std::string())
              << (F.phaseField ? " PHASEFIELD " + num(F.pfEta) + " " + num(F.pfH) + " " + num(F.pfE) + " " + num(F.pfNu) + (F.pfTwoWay ? " TWOWAY" : "") : std::string()) << "\n";
    if (!g.hfracs.empty()) o << "\n";
    if (g.langVL > 0 && g.langPL > 0) o << "-- RESSIM: LANGMUIR " << num(g.langVL) << ' ' << num(g.langPL) << ' ' << num(g.langRho) << "\n\n";
    if (g.permGamma > 0) o << "-- RESSIM: PERMSTRESS " << num(g.permGamma) << "\n\n";
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
    for (const auto& T : c.tracerProps)   // C27 tracer decay / adsorption (ResSim extension of the verbatim TRACER keyword)
    {
        o << "-- RESSIM: TRACERPROP " << T.name << ' ' << num(T.halfLife) << ' ' << num(T.kd) << ' ' << num(T.init);
        if (T.kPart > 0 || !T.parent.empty() || T.alphaL > 0 || T.dm > 0)
            o << ' ' << num(T.kPart) << ' ' << (T.parent.empty() ? std::string("1*") : T.parent) << ' ' << num(T.alphaL) << ' ' << num(T.dm);
        o << "\n";
    }
    if (!c.tracerProps.empty()) o << "\n";
    for (const auto& B : c.wellbores)   // C29 wellbore friction / ICD (ResSim extension)
        o << "-- RESSIM: WELLBORE " << B.well << ' ' << num(B.id) << ' ' << num(B.rough) << ' ' << num(B.icdK) << ' ' << num(B.rhoCal) << "\n";
    if (!c.mswOpt.flash || !c.mswOpt.storage || c.mswOpt.xflow || c.mswOpt.heads != -1) {   // C37 (defaults are not written)
        o << "-- RESSIM: MSWOPTS " << (c.mswOpt.flash ? 1 : 0) << ' ' << (c.mswOpt.storage ? 1 : 0);
        if (c.mswOpt.xflow || c.mswOpt.heads != -1) o << ' ' << (c.mswOpt.xflow ? 1 : 0);
        if (c.mswOpt.heads != -1) o << ' ' << c.mswOpt.heads;
        o << "\n";
    }
    if (c.hfField.on >= 0)
        o << "-- RESSIM: HFFIELD " << (c.hfField.on ? "ON" : "OFF") << ' ' << num(c.hfField.E) << ' ' << num(c.hfField.nu) << ' ' << num(c.hfField.radius) << ' ' << num(c.hfField.budget) << "\n";
    if (c.salt.on)                      // C29 EWASG-type salt (ResSim extension)
        o << "-- RESSIM: SALT " << num(c.salt.c0) << ' ' << num(c.salt.cSol) << ' ' << num(c.salt.rhoSalt) << ' ' << num(c.salt.phiR) << ' '
          << num(c.salt.gamma) << ' ' << num(c.salt.tempC) << ' ' << num(c.salt.cInj) << "\n";
    if (!c.wellbores.empty() || c.salt.on || c.hfField.on >= 0 || !c.mswOpt.flash || !c.mswOpt.storage) o << "\n";
    {   // C24 fault mechanics / seismicity and in-run fracture contact (ResSim extensions, comment lines for OPM / ECLIPSE)
        const FaultMechInput& F = c.faultMech;
        if (F.svGrad > 0 || !F.faults.empty()) {
            o << "-- RESSIM: STRESSFIELD " << num(F.svGrad) << ' ' << num(F.shRatio) << ' ' << num(F.sHRatio) << ' ' << num(F.azimuth) << ' '
              << num(F.stressPath) << ' ' << num(F.sv0) << (F.source == 1 ? " FE" : "");
            if (F.thE > 0 && F.thAlpha > 0) o << " TH " << num(F.thE) << ' ' << num(F.thAlpha) << ' ' << num(F.thNu);
            o << "\n";
            for (const auto& R : F.faults)
                o << "-- RESSIM: FAULTMECH " << R.name << ' ' << num(R.mu) << ' ' << num(R.coh) << ' ' << num(R.aSig) << ' ' << num(R.ta) << ' '
                  << num(R.r0) << ' ' << num(R.bValue) << ' ' << num(R.mMin) << ' ' << num(R.permFactor)
                  << (std::fabs(R.permRamp - 0.05) > 1e-12 ? " " + num(R.permRamp) : std::string()) << "\n";
            if (F.slipOn)
                o << "-- RESSIM: FAULTSLIP " << num(F.slipE) << ' ' << num(F.slipNu) << ' ' << num(F.slipDMu) << ' ' << num(F.slipRefMm) << ' '
                  << num(F.slipCutoff) << (F.slipHalfSpace ? "" : " FULLSPACE") << "\n";
            for (const auto& R : F.faults)
                if (R.rsfLayer > 0) o << "-- RESSIM: FAULTRSF " << R.name << ' ' << R.rsfLayer << ' ' << num(R.rsfA) << ' ' << num(R.rsfB) << ' ' << num(R.rsfDcMm) << "\n";
            o << "\n";
        }
        // C25 sequential FE geomechanics of black-oil / compositional cases (thermal cases write thermal.geo themselves)
        if (c.geoSeq.enabled && !c.thermal.enabled) o << geoKeywordsText(c.geoSeq) << "\n";
        // imported (verbatim) thermal decks: the thermal keyword writer only emits the RESSIM geomechanics lines in
        // structured mode, so the THM mechanics of an imported EGS deck is written here
        if (c.thermal.enabled && !c.thermal.structured && c.thermal.geo.enabled) o << geoKeywordsText(c.thermal.geo) << "\n";
        const GeoCoupleSpec& G = c.geoCouple;
        if (G.on) o << "-- RESSIM: GEOCOUPLE " << num(G.sxx) << ' ' << num(G.syy) << ' ' << num(G.eGPa) << ' ' << num(G.nu) << ' ' << num(G.friction) << ' ' << G.layer << "\n\n";
    }
    if (!g.nnc.empty()) {
        o << "NNC\n";
        for (auto& e : g.nnc)
            o << ' ' << e.i1 + 1 << ' ' << e.j1 + 1 << ' ' << e.k1 + 1 << "  " << e.i2 + 1 << ' ' << e.j2 + 1 << ' ' << e.k2 + 1 << "  " << num(e.trans) << " /\n";
        o << "/\n\n";
    }
    // local grid refinements: CARFIN + NXFIN.. + local grid arrays (resolved values) + ENDFIN; the other local arrays
    // go to their sections with REFINE (lgrBlock)
    Case lgrResolved;
    if (!g.lgrs.empty()) { lgrResolved = c; try { expandLgrArrays(lgrResolved); } catch (...) { lgrResolved.grid.lgrs.clear(); } }
    auto lgrArrays = [&](size_t l, const std::set<std::string>& names) {
        std::string out;
        if (lgrResolved.grid.lgrs.size() != g.lgrs.size()) return out;
        const LgrInput& L = g.lgrs[l];
        size_t off = size_t(g.size());
        for (size_t m = 0; m < l; ++m) off += size_t(g.lgrs[m].size());
        GridInput& R = lgrResolved.grid;
        const std::pair<const char*, const std::vector<double>*> dbl[] = {
            {"PORO", &R.poro}, {"PERMX", &R.permx}, {"PERMY", &R.permy}, {"PERMZ", &R.permz}, {"NTG", &R.ntg}, {"MULTX", &R.multx},
            {"MULTY", &R.multy}, {"MULTZ", &R.multz}, {"MULTPV", &R.multpv}, {"SWATINIT", &R.swatinit}, {"SWL", &R.swl}, {"SWCR", &R.swcr},
            {"SWU", &R.swu}, {"SGL", &R.sgl}, {"SGCR", &R.sgcr}, {"SGU", &R.sgu}, {"SOWCR", &R.sowcr}, {"SOGCR", &R.sogcr},
            {"ISWL", &R.iswl}, {"ISWCR", &R.iswcr}, {"ISWU", &R.iswu}, {"ISGL", &R.isgl}, {"ISGCR", &R.isgcr}, {"ISGU", &R.isgu},
            {"ISOWCR", &R.isowcr}, {"ISOGCR", &R.isogcr}, {"SPOLY", &R.spoly}, {"PRESSURE", &lgrResolved.init.pressure},
            {"SWAT", &lgrResolved.init.swat}, {"SGAS", &lgrResolved.init.sgas}, {"RS", &lgrResolved.init.rs}, {"RV", &lgrResolved.init.rv}};
        const std::pair<const char*, const std::vector<int>*> ints[] = {
            {"ACTNUM", &R.actnum}, {"SATNUM", &R.satnum}, {"PVTNUM", &R.pvtnum}, {"EQLNUM", &R.eqlnum}, {"FIPNUM", &R.fipnum},
            {"IMBNUM", &R.imbnum}, {"PLMIXNUM", &R.plmixnum}, {"MULTNUM", &R.multnum}, {"FLUXNUM", &R.fluxnum}, {"OPERNUM", &R.opernum}};
        std::ostringstream os;
        for (const auto& kv : L.props) {
            if (!names.count(kv.first)) continue;
            for (const auto& d : dbl)
                if (kv.first == d.first && d.second->size() >= off + size_t(L.size())) {
                    std::vector<double> v(d.second->begin() + std::ptrdiff_t(off), d.second->begin() + std::ptrdiff_t(off + size_t(L.size())));
                    os << d.first << "\n"; writeArray(os, v); os << "/\n";
                }
            for (const auto& d : ints)
                if (kv.first == d.first && d.second->size() >= off + size_t(L.size())) {
                    std::vector<int> v(d.second->begin() + std::ptrdiff_t(off), d.second->begin() + std::ptrdiff_t(off + size_t(L.size())));
                    os << d.first << "\n"; writeArray(os, v); os << "/\n";
                }
        }
        return os.str();
    };
    auto lgrBlock = [&](const std::set<std::string>& names) {
        for (size_t l = 0; l < g.lgrs.size(); ++l) {
            const std::string a = lgrArrays(l, names);
            if (!a.empty()) o << "REFINE\n '" << g.lgrs[l].name << "' /\n" << a << "ENDFIN\n\n";
        }
    };
    for (size_t l = 0; l < g.lgrs.size(); ++l) {
        const LgrInput& L = g.lgrs[l];
        int nwIn = 0;   // NWMAX must cover the wells declared in the refinement (WELSPECL)
        for (const auto& w : c.wells) if (w.lgr == int(l)) ++nwIn;
        const int nwmax = std::max(L.nwmax, nwIn);
        o << "CARFIN\n '" << L.name << "' " << L.i1 + 1 << ' ' << L.i2 + 1 << ' ' << L.j1 + 1 << ' ' << L.j2 + 1 << ' ' << L.k1 + 1 << ' ' << L.k2 + 1
          << ' ' << L.nx << ' ' << L.ny << ' ' << L.nz << ' ' << (nwmax > 0 ? std::to_string(nwmax) : std::string("1*"))
          << (L.parent != "GLOBAL" && !L.parent.empty() ? " '" + L.parent + "'" : std::string()) << " /\n";
        auto fin = [&](const char* nm, const auto& v) { if (!v.empty()) { o << nm << "\n"; writeArray(o, v); o << "/\n"; } };
        fin("NXFIN", L.nxfin); fin("NYFIN", L.nyfin); fin("NZFIN", L.nzfin);
        fin("HXFIN", L.hxfin); fin("HYFIN", L.hyfin); fin("HZFIN", L.hzfin);
        o << lgrArrays(l, {"PORO", "PERMX", "PERMY", "PERMZ", "NTG", "MULTX", "MULTY", "MULTZ", "MULTPV", "ACTNUM", "MULTNUM", "FLUXNUM", "OPERNUM"});
        o << "ENDFIN\n\n";
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
    if (thermOn) o << thermalKeywordsText(c, 4);
    writeExtra(DeckSection::Grid);
    if (!g.tranMultX.empty() || !g.tranMultY.empty() || !g.tranMultZ.empty()) {
        // C38 EDIT-section transmissibility edits: MULTIPLY records, runs of equal factors along I
        o << "EDIT\n\nMULTIPLY\n";
        const char* nm[3] = {"TRANX", "TRANY", "TRANZ"};
        const std::vector<double>* src[3] = {&g.tranMultX, &g.tranMultY, &g.tranMultZ};
        for (int d = 0; d < 3; ++d) {
            const auto& v = *src[d];
            if (int(v.size()) != g.nx * g.ny * g.nz) continue;
            for (int k = 0; k < g.nz; ++k)
                for (int j = 0; j < g.ny; ++j)
                    for (int i = 0; i < g.nx;) {
                        const double f = v[size_t(i + g.nx * (j + g.ny * k))];
                        int i2 = i;
                        while (i2 + 1 < g.nx && v[size_t(i2 + 1 + g.nx * (j + g.ny * k))] == f) ++i2;
                        if (f != 1.0)
                            o << " '" << nm[d] << "' " << num(f) << ' ' << i + 1 << ' ' << i2 + 1 << ' ' << j + 1 << ' ' << j + 1 << ' ' << k + 1 << ' ' << k + 1 << " /\n";
                        i = i2 + 1;
                    }
        }
        o << "/\n\n";
        writeExtra(DeckSection::Edit);
    } else writeExtra(DeckSection::Edit, "EDIT");
    o << "PROPS\n\n";
    std::vector<PvtRegion> PR; for (int r = 0; r < npvt; ++r) PR.push_back(p.pvtRegion(r));
    std::vector<SatRegion> SR; for (int r = 0; r < nsat; ++r) SR.push_back(p.satRegion(r));
    if (!co2) {   // CO2STORE: water / CO2 properties come from the CO2-brine model
        o << "PVTW\n";
        for (auto& R : PR) o << ' ' << num(R.pvtw.pref) << ' ' << num(R.pvtw.bw) << ' ' << num(R.pvtw.cw) << ' ' << num(R.pvtw.visc) << ' ' << num(R.pvtw.cv) << " /\n";
        o << "\n";
    }
    o << "ROCK\n";
    for (auto& R : PR) o << ' ' << num(R.rock.pref) << ' ' << num(R.rock.cr) << " /\n";
    if (!co2) {
        o << "\nDENSITY\n";
        for (auto& R : PR) o << ' ' << num(R.density.oil) << ' ' << num(R.density.water) << ' ' << num(R.density.gas) << " /\n";
    }
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
    if (compS) o << compositionalKeywordsText(c, 1);
    else if (compOn) o << trapRessimLines(c);
    if (thermOn) o << thermalKeywordsText(c, 1);
    if (compOn || thermOn) {   // the EOS / CO2-brine / thermal fluid replaces the black-oil tables
    } else if (c.phases.oil) {
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
    if (compOn || thermOn) {
    } else if (vap) {
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
    if (c.phases.gas && c.phases.water && !c.phases.oil) {   // gas-water (CO2STORE): SWFN + SGFN
        o << "SWFN\n";
        for (auto& S : SR) { for (auto& r : S.swof) o << " " << num(r.s) << ' ' << num(r.krw) << ' ' << num(r.pc) << "\n"; o << "/\n"; }
        o << "\nSGFN\n";
        for (auto& S : SR) { for (auto& r : S.sgof) o << " " << num(r.s) << ' ' << num(r.krw) << ' ' << num(r.pc) << "\n"; o << "/\n"; }
        o << "\n";
    }
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
        if (g.scalecrs) o << "SCALECRS\n 'YES' /\n\n";
        if (c.props.hyst.on) {
            // imbibition end points: cells without a value get the end point of their IMBNUM table
            auto iArr = [&](const char* name, const std::vector<double>& v, double SatFunctions::EndPoints::*m) {
                if (v.empty()) return;
                std::vector<double> w = v;
                for (size_t q = 0; q < w.size(); ++q) {
                    if (w[q] >= 0) continue;
                    const std::vector<int>& rn = g.imbnum.empty() ? g.satnum : g.imbnum;
                    const int r = q < rn.size() ? std::clamp(rn[q] - 1, 0, int(sf.size()) - 1) : 0;
                    w[q] = sfOk ? sf[size_t(r)].tableEndPoints().*m : 0.0;
                }
                arr(name, w);
            };
            iArr("ISWL", g.iswl, &EP::swl); iArr("ISWCR", g.iswcr, &EP::swcr); iArr("ISWU", g.iswu, &EP::swu);
            iArr("ISGL", g.isgl, &EP::sgl); iArr("ISGCR", g.isgcr, &EP::sgcr); iArr("ISGU", g.isgu, &EP::sgu);
            iArr("ISOWCR", g.isowcr, &EP::sowcr); iArr("ISOGCR", g.isogcr, &EP::sogcr);
        }
    }
    if (c.props.hyst.on) {
        const auto& h = c.props.hyst;
        o << "EHYSTR\n " << num(h.curv) << " " << h.model << " " << num(h.killough) << " 1* '" << h.what << "' /\n\n";
    }
    if (c.props.vap1 >= 0) o << "VAPPARS\n " << num(c.props.vap1) << " " << num(std::max(0.0, c.props.vap2)) << " /\n\n";
    lgrBlock({"SWATINIT", "SWL", "SWCR", "SWU", "SGL", "SGCR", "SGU", "SOWCR", "SOGCR", "ISWL", "ISWCR", "ISWU", "ISGL", "ISGCR", "ISGU", "ISOWCR", "ISOGCR"});

    writeExtra(DeckSection::Props);
    // ---------------- REGIONS
    const bool mixOut = polyOn && !g.plmixnum.empty();
    const bool imbOut = c.props.hyst.on && !g.imbnum.empty();
    const std::set<std::string> lgrRegionNames = {"SATNUM", "PVTNUM", "EQLNUM", "FIPNUM", "IMBNUM", "PLMIXNUM"};
    bool lgrRegions = false;
    for (const auto& L : g.lgrs) for (const auto& kv : L.props) if (lgrRegionNames.count(kv.first)) lgrRegions = true;
    const bool krOut = !g.krnumx.empty() || !g.krnumy.empty() || !g.krnumz.empty() || !g.rocknum.empty() || !g.surfnum.empty();
    if (!g.satnum.empty() || !g.pvtnum.empty() || !g.eqlnum.empty() || !g.fipnum.empty() || mixOut || imbOut || lgrRegions || krOut) {
        o << "REGIONS\n\n";
        arr("SATNUM", g.satnum); arr("PVTNUM", g.pvtnum); arr("EQLNUM", g.eqlnum); arr("FIPNUM", g.fipnum);
        if (imbOut) arr("IMBNUM", g.imbnum);
        if (!g.krnumx.empty()) arr("KRNUMX", g.krnumx);
        if (!g.krnumy.empty()) arr("KRNUMY", g.krnumy);
        if (!g.krnumz.empty()) arr("KRNUMZ", g.krnumz);
        if (!g.rocknum.empty()) arr("ROCKNUM", g.rocknum);
        if (mixOut) arr("PLMIXNUM", g.plmixnum);
        if (!g.surfnum.empty()) arr("SURFNUM", g.surfnum);   // C35
        lgrBlock(lgrRegionNames);
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
        lgrBlock({"PRESSURE", "SWAT", "SGAS", "RS", "RV"});
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
    if (!g.sfoam.empty()) arr("SFOAM", g.sfoam);   // C35 initial foam concentration
    if (!g.ssol.empty()) arr("SSOL", g.ssol);      // C35 initial solvent saturation
    if (polyOn) lgrBlock({"SPOLY"});
    if (compS) o << compositionalKeywordsText(c, 2);
    if (thermOn) o << thermalKeywordsText(c, 2);
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
    o << "-- RESSIM: TUNING " << num(t.dpMaxRel) << ' ' << t.threads << ' ' << t.linearSolver << ' ' << t.cprWeights << ' ' << (t.mixedPrecision ? 1 : 0) << ' ' << t.gpu << ' ' << t.amgReuse;
    if (t.dsTarget > 0 || t.dpTarget > 0 || t.leanMemory != -1) o << ' ' << num(t.dsTarget) << ' ' << num(t.dpTarget);
    if (t.leanMemory != -1) o << ' ' << t.leanMemory;
    o << "\n\n";
    // restart arrays: compositional runs also write the phase / overall compositions per component, densities,
    // viscosities and saturation pressure (E300 / OPM mnemonics) so the results viewers can show every detail
    if (c.comp.enabled && !c.comp.co2store) o << "RPTRST\n 'BASIC=2' 'XMF' 'YMF' 'ZMF' 'DEN' 'VISC' 'PSAT'" << (c.comp.fluid.co2sol ? " 'AMF' 'RSW'" : "") << " /\n\n";
    else if (c.comp.enabled && c.comp.co2store) o << "RPTRST\n 'BASIC=2' 'DEN' 'VISC' /\n\n";
    else o << "RPTRST\n 'BASIC=2' /\n\n";
    // COMPDAT for the connections that open (or shut) at report step s
    auto compRecord = [&](const WellInput& w, const Completion& q, bool open) {
        std::ostringstream r;
        r << " '" << w.name << "' ";
        if (q.lgr >= 0 && q.lgr < int(g.lgrs.size())) r << "'" << g.lgrs[size_t(q.lgr)].name << "' ";
        r << q.i + 1 << ' ' << q.j + 1 << ' ' << q.k1 + 1 << ' ' << q.k2 + 1 << ' ' << (open ? "OPEN" : "SHUT") << ' ' << (q.satTable > 0 ? std::to_string(q.satTable) : std::string("1*")) << ' ';
        if (q.cf > 0) r << num(q.cf); else r << "1*";
        r << ' ' << num(q.diameter) << ' ';
        if (q.kh > 0) r << num(q.kh); else r << "1*";
        r << ' ' << num(q.skin) << ' ' << (q.dfactor >= 0 ? num(q.dfactor) : std::string("1*")) << ' ' << char(std::toupper(q.dir));
        if (q.r0 > 0) r << ' ' << num(q.r0);
        r << " /\n";
        return r.str();
    };
    auto compEvents = [&](size_t s, bool local) {
        std::string out;
        for (auto& w : c.wells)
            for (auto& q : w.completions) {
                if ((q.lgr >= 0) != local) continue;
                if (q.defStep >= 0 && q.defStep < q.fromStep && (s == 0 ? q.defStep <= 0 : q.defStep == int(s))) out += compRecord(w, q, false);   // C46
                if (s == 0 ? q.fromStep <= 0 : q.fromStep == int(s)) out += compRecord(w, q, q.open);
                if (s > 0 && q.open && q.untilStep == int(s)) out += compRecord(w, q, false);
            }
        return out;
    };
    auto writeCompdat = [&](size_t s) {
        const std::string e = compEvents(s, false); if (!e.empty()) o << "COMPDAT\n" << e << "/\n\n";
        const std::string l = compEvents(s, true); if (!l.empty()) o << "COMPDATL\n" << l << "/\n\n";
        // C44: WPIMULT factors of the open periods starting at this step (the connections were just (re)defined above, so the
        // cumulative factor applies once); one record per layer
        std::string pm;
        for (auto& w : c.wells)
            for (auto& q : w.completions) {
                if ((q.mult == 1.0 && w.piMult == 1.0) || !q.open || q.lgr >= 0 || !(s == 0 ? q.fromStep <= 0 : q.fromStep == int(s))) continue;
                // the whole-run multiplier (UI / studies) is baked in: OPM applies it too. K = 0 (all layers of the column) unless
                // the well has another connection range in this column (then one record per layer)
                bool other = false;
                for (auto& x : w.completions) if (x.lgr < 0 && x.i == q.i && x.j == q.j && (x.k1 != q.k1 || x.k2 != q.k2)) other = true;
                const std::string head = " '" + w.name + "' " + num(w.piMult * q.mult) + " " + std::to_string(q.i + 1) + " " + std::to_string(q.j + 1) + " ";
                if (!other) pm += head + "0 /\n";
                else for (int kk = q.k1; kk <= q.k2; ++kk) pm += head + std::to_string(kk + 1) + " /\n";
            }
        if (!pm.empty()) o << "WPIMULT\n" << pm << "/\n\n";
    };
    if (!c.wells.empty()) {
        for (int local = 0; local < 2; ++local) {
            bool any = false;
            for (auto& w : c.wells) if ((w.lgr >= 0 && w.lgr < int(g.lgrs.size())) == (local == 1)) any = true;
            if (!any) continue;
            o << (local ? "WELSPECL\n" : "WELSPECS\n");
            for (auto& w : c.wells) {
                if ((w.lgr >= 0 && w.lgr < int(g.lgrs.size())) != (local == 1)) continue;
                const char* ph = w.type == WellType::Producer ? "OIL" : w.type == WellType::WaterInjector ? "WATER" : "GAS";
                o << " '" << w.name << "' '" << (w.group.empty() ? std::string("G") : w.group) << "' ";
                if (local) o << "'" << g.lgrs[size_t(w.lgr)].name << "' ";
                o << w.i + 1 << ' ' << w.j + 1 << ' ';
                if (w.refDepth > 0) o << num(w.refDepth); else o << "1*";
                o << " " << ph;
                if (!w.specsTail.empty()) o << ' ' << w.specsTail;
                o << " /\n";
            }
            o << "/\n\n";
        }
        writeCompdat(0);
        // C44: a whole-run PI multiplier is written into the per-connection WPIMULT records (also for connections opened later);
        // this comment lets ResSim restore it as the well multiplier on import
        for (auto& w : c.wells) if (w.piMult != 1.0) o << "-- RESSIM: PIMULT " << w.name << ' ' << num(w.piMult) << "\n";
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
    if (!c.propModels.empty()) {   // C40
        o << "-- Geostatistical property models re-simulated by Studies (ResSim extension)\n";
        for (const auto& m : c.propModels) {
            o << "-- RESSIM: PROPMODEL " << m.property << ' ' << m.method << ' ' << m.model << ' ' << num(m.nugget) << ' ' << num(m.sill) << ' ' << num(m.rangeMax) << ' '
              << num(m.rangeMin) << ' ' << num(m.rangeVert) << ' ' << num(m.azimuthDeg) << ' ' << m.seed << ' ' << m.maxNbr << ' ' << int(m.fromWells) << ' '
              << num(m.minV) << ' ' << num(m.maxV) << ' ' << int(m.permFromPoro) << ' ' << num(m.ppA) << ' ' << num(m.ppB) << ' ' << num(m.kvkh) << " '" << m.source << "'\n";
            if (m.property == 3) {   // C41
                o << "-- RESSIM: PROPFAC " << int(m.facSatnum);
                for (size_t k = 0; k < m.facProp.size(); ++k)
                    o << ' ' << num(m.facProp[k]) << ' ' << num(k < m.facPoro.size() ? m.facPoro[k] : -1) << ' ' << num(k < m.facPerm.size() ? m.facPerm[k] : -1);
                o << "\n";
            }
            for (size_t q = 0; q < m.data.size(); q += 4) {
                o << "-- RESSIM: PROPDATA";
                for (size_t k = q; k < std::min(m.data.size(), q + 4); ++k) o << ' ' << num(m.data[k][0]) << ' ' << num(m.data[k][1]) << ' ' << num(m.data[k][2]) << ' ' << num(m.data[k][3]);
                o << "\n";
            }
        }
        o << "\n";
    }
    auto findWell = [&](const std::string& nm) -> const WellInput* { for (auto& w : c.wells) if (w.name == nm) return &w; return nullptr; };
    auto dflt = [&](double v) { return v < 0 ? std::string("1*") : num(v); };
    auto udv = [&](const WellControl& wc, const char* t, double v) { const std::string* u = wc.udaOf(t); return u ? "'" + *u + "'" : dflt(v); };
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
        std::ostringstream prod, inj, hist, histLim;
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
                     << dflt(wc.grat) << ' ' << vfp << ' ' << dflt(wc.alq) << ' ' << dflt(wc.thp) << ' ' << dflt(wc.bhpHist > 0 ? wc.bhpHist : wc.mode == Control::BHP ? wc.bhp : -1.0) << " /\n";
            if (wc.history && isProd && wc.mode != Control::BHP && wc.bhp > 0 && !wc.bhpDefault)   // a BHP limit of a history well: WELTARG
                histLim << " '" << wc.well << "' BHP " << num(wc.bhp) << " /\n";
            else if (isProd)
                prod << " '" << wc.well << "' " << st << ' ' << mode << ' ' << udv(wc, "ORAT", wc.orat) << ' ' << udv(wc, "WRAT", wc.wrat) << ' '
                     << udv(wc, "GRAT", wc.grat) << ' ' << udv(wc, "LRAT", wc.lrat) << ' ' << udv(wc, "RESV", wc.resv) << ' ' << udv(wc, "BHP", wc.bhpDefault ? -1.0 : wc.bhp) << ' ' << udv(wc, "THP", wc.thp) << ' ' << vfp << ' '
                     << dflt(wc.alq) << " /\n";
            else {
                bool gas = wc.phase ? wc.phase == 2 : (w && w->type == WellType::GasInjector);
                if (!w && !wc.phase) for (auto& x : c.wells) if (matchName(wc.well, x.name)) { gas = x.type == WellType::GasInjector; break; }
                const bool resvStar = wc.deckMode == "RESV*";     // RESV injector without a surface-rate limit
                const std::string im = resvStar ? std::string("RESV") : mode;
                inj << " '" << wc.well << "' " << (gas ? "GAS" : "WATER") << ' ' << st << ' '
                    << im << ' ' << (resvStar ? std::string("1*") : udv(wc, "RATE", wc.rate)) << ' ' << udv(wc, "RESV", wc.resv) << ' ' << udv(wc, "BHP", wc.bhpDefault ? -1.0 : wc.bhp) << ' '
                    << udv(wc, "THP", wc.thp) << ' ' << vfp << " /\n";
            }
        }
        o << extraAtStep(int(s), false, 1);
        if (compS && s == 0) o << compositionalKeywordsText(c, 3);
        if (thermOn) o << thermalKeywordsText(c, 3, int(s));
        auto trimDefaults = [](std::string s) {          // " 1* 1* /" at the end of a record == " /"
            for (size_t p; (p = s.find(" 1* /\n")) != std::string::npos;) s.erase(p, 3);
            return s;
        };
        if (!prod.str().empty()) o << "WCONPROD\n" << trimDefaults(prod.str()) << "/\n\n";
        if (!hist.str().empty()) o << "WCONHIST\n" << hist.str() << "/\n\n";
        if (!histLim.str().empty()) o << "WELTARG\n" << histLim.str() << "/\n\n";
        if (!inj.str().empty()) o << "WCONINJE\n" << trimDefaults(inj.str()) << "/\n\n";
        {   // DRSDT and VAPPARS exclude each other: keep the deck order (the later keyword wins)
            const bool wd = c.schedule[s].drsdt > -1.5 && c.phases.disgas, wv = c.schedule[s].vap1 >= 0 && c.phases.oil && c.phases.gas;
            auto drs = [&] { o << "DRSDT\n " << (c.schedule[s].drsdt < 0 ? std::string("1e20") : num(c.schedule[s].drsdt)) << (c.schedule[s].drsdtFree ? " 'FREE'" : "") << " /\n\n"; };
            auto vap = [&] { o << "VAPPARS\n " << num(c.schedule[s].vap1) << " " << num(std::max(0.0, c.schedule[s].vap2)) << " /\n\n"; };
            if (wd && wv && c.schedule[s].vapAfterDrsdt) { drs(); vap(); }
            else { if (wv) vap(); if (wd) drs(); }
        }
        bool rawGuide = false, rawGprod = false, rawGinj = false;
        for (const ExtraKeyword& e : c.extraKeywords) {
            if (e.section != DeckSection::Schedule || e.step != int(s)) continue;
            if (e.name == "GUIDERAT") rawGuide = true;
            if (e.name == "GCONPROD" || e.name == "GCONPRI") rawGprod = true;
            if (e.name == "GCONINJE") rawGinj = true;
        }
        // group controls built in the case (not kept as deck keywords): written from the structures
        auto opt = [&](double v) { return v >= 0 ? num(v) : std::string("1*"); };
        if (!rawGprod) {
            std::ostringstream g;
            for (const GroupProdControl& gp : c.schedule[s].gprod) {
                if (gp.priority) continue;
                const char* m = gp.mode == Control::ORAT ? "ORAT" : gp.mode == Control::WRAT ? "WRAT" : gp.mode == Control::GRAT ? "GRAT" :
                                gp.mode == Control::LRAT ? "LRAT" : gp.mode == Control::RESV ? "RESV" : "NONE";
                g << " '" << gp.group << "' " << m << ' ' << opt(gp.orat) << ' ' << opt(gp.wrat) << ' ' << opt(gp.grat) << ' ' << opt(gp.lrat) << ' ' << gp.action;
                if (!gp.actW.empty() || !gp.actG.empty() || !gp.actL.empty() || gp.resv >= 0 || gp.guideRate >= 0 || !gp.guideDef.empty()) {
                    auto a = [](const std::string& x) { return x.empty() ? std::string("1*") : x; };
                    g << " 1* " << opt(gp.guideRate) << ' ' << a(gp.guideDef) << ' ' << a(gp.actW) << ' ' << a(gp.actG) << ' ' << a(gp.actL);
                    if (gp.resv >= 0) g << ' ' << num(gp.resv);
                }
                g << " /\n";
            }
            if (!g.str().empty()) o << "GCONPROD\n" << g.str() << "/\n\n";
        }
        if (!rawGinj && !c.schedule[s].ginj.empty()) {
            o << "GCONINJE\n";
            for (const GroupInjControl& gi : c.schedule[s].ginj)
                o << " '" << gi.group << "' " << (gi.phase == 2 ? "GAS" : "WATER") << ' ' << gi.mode << ' ' << opt(gi.rate) << ' ' << opt(gi.resv) << ' '
                  << opt(gi.rein) << ' ' << opt(gi.vrep) << " /\n";
            o << "/\n\n";
        }
        if (c.schedule[s].guiderat.set && !rawGuide) {
            const GuideRateSpec& g = c.schedule[s].guiderat;
            o << "GUIDERAT\n " << num(g.minInterval) << " '" << g.phase << "' " << num(g.A) << ' ' << num(g.B) << ' ' << num(g.C) << ' '
              << num(g.D) << ' ' << num(g.E) << ' ' << num(g.F) << ' ' << (g.allowIncrease ? "YES" : "NO") << ' ' << num(g.damping) << " /\n\n";
        }
        if (compS && !c.comp.co2store && !c.schedule[s].injStreams.empty()) {
            o << "WINJGAS\n";
            for (auto& ws : c.schedule[s].injStreams) o << " '" << ws.first << "' STREAM '" << ws.second << "' /\n";
            o << "/\n\n";
        }
        if (polyOn && !c.schedule[s].wpolymer.empty()) {
            o << "WPOLYMER\n";
            for (auto& wp : c.schedule[s].wpolymer) o << " '" << wp.well << "' " << num(wp.conc) << " 0 /\n";
            o << "/\n\n";
        }
        o << extraAtStep(int(s), false, 0);
        std::vector<double> lens{c.schedule[s].length};
        size_t e2 = s + 1;
        while (e2 < c.schedule.size() && c.schedule[e2].controls.empty() && c.schedule[e2].wpolymer.empty() && c.schedule[e2].drsdt < -1.5 && c.schedule[e2].injStreams.empty() && !c.schedule[e2].guiderat.set && c.schedule[e2].gprod.empty() && c.schedule[e2].ginj.empty() && compEvents(e2, false).empty() && compEvents(e2, true).empty() && extraAtStep(int(e2), false).empty())
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
