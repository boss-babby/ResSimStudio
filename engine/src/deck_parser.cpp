// Generic ECLIPSE / OPM Flow deck reader: lexer (INCLUDE, PATHS, SKIP, raw blocks) + keyword-structure parser.
#include "ressim/deck_parser.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <cstdio>
#include <stdexcept>

namespace ressim {
namespace deck {

// ============================================================================ keyword table
namespace {
const KeywordSpec kTable[] = {
#include "keyword_db.inc"
};
struct FamilySpec { const char* rx; uint16_t sections; KwSize kind; };
const FamilySpec kFamilies[] = {
#include "keyword_regex.inc"
};
const char* kSections[8] = {"RUNSPEC", "GRID", "EDIT", "PROPS", "REGIONS", "SOLUTION", "SUMMARY", "SCHEDULE"};

std::string upper(std::string s) { for (char& c : s) c = char(std::toupper((unsigned char)c)); return s; }
std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace((unsigned char)s[a])) ++a;
    while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
}
bool parseNumber(const std::string& s0, double& v) {
    if (s0.empty()) return false;
    const char c0 = s0[0];
    if (!(std::isdigit((unsigned char)c0) || c0 == '-' || c0 == '+' || c0 == '.')) return false;
    std::string s = s0;
    for (char& c : s) if (c == 'D' || c == 'd') c = 'E';
    char* end = nullptr;
    v = std::strtod(s.c_str(), &end);
    return end && *end == 0;
}
}  // namespace

const std::vector<KeywordSpec>& keywordTable() {
    static const std::vector<KeywordSpec> t(std::begin(kTable), std::end(kTable));
    return t;
}
const char* sectionName(int bit) { return (bit >= 0 && bit < 8) ? kSections[bit] : ""; }

const KeywordSpec* findKeyword(const std::string& name) {
    const auto& t = keywordTable();
    auto it = std::lower_bound(t.begin(), t.end(), name, [](const KeywordSpec& a, const std::string& b) { return std::strcmp(a.name, b.c_str()) < 0; });
    if (it != t.end() && name == it->name) return &*it;
    // summary-style families (WOPR, FOPT, ...) - matched against generated regular expressions
    static std::vector<std::pair<std::regex, KeywordSpec>> fam;
    static bool init = false;
    if (!init) {
        init = true;
        for (const auto& f : kFamilies) {
            try { fam.push_back({std::regex(f.rx), KeywordSpec{"(summary)", f.sections, f.kind, 1, "", -1, 1, 0, ""}}); }
            catch (...) {}
        }
    }
    if (name.size() < 2 || name.size() > 8) return nullptr;
    for (auto& f : fam) if (std::regex_match(name, f.first)) return &f.second;
    return nullptr;
}

bool Record::isNumber(size_t i) const { double v; return has(i) && parseNumber(items[i].s, v); }
double Record::num(size_t i, double def) const {
    if (!has(i)) return def;
    double v;
    if (parseNumber(items[i].s, v)) return v;
    throw std::runtime_error("expected a number in item " + std::to_string(i + 1) + " but found '" + items[i].s + "'");
}

std::string RawDeck::where(int file, int line) const {
    std::string f = (file >= 0 && file < int(files.size())) ? files[file] : std::string("deck");
    const size_t s = f.find_last_of("/\\");
    if (s != std::string::npos) f = f.substr(s + 1);
    return f + ":" + std::to_string(line);
}
void RawDeck::message(Level l, const Keyword& k, const std::string& text, int line) {
    messages.push_back(Message{l, k.name, k.file, line > 0 ? line : k.line, text});
}

// ============================================================================ lexer
namespace {

struct Tk {
    enum T : uint8_t { Kw, Data, Slash, Raw } t = Data;
    std::string s;
    bool quoted = false;
    bool known = true;          // keyword tokens: found in the table
    bool candidate = false;     // keyword tokens: unknown name on a keyword-like line
    int file = 0, line = 0;
    int array = -1;             // keyword tokens of Data kind: index into Lexer::arrays
    bool arrayClosed = false;
};

bool keywordLike(const std::string& w) {
    if (w.empty() || w.size() > 8 || !std::isupper((unsigned char)w[0])) return false;
    for (char ch : w) if (!(std::isupper((unsigned char)ch) || std::isdigit((unsigned char)ch) || ch == '_' || ch == '+' || ch == '-')) return false;
    return true;
}

std::string stripComment(const std::string& l) {
    char q = 0;
    for (size_t i = 0; i + 1 < l.size(); ++i) {
        if (q) { if (l[i] == q) q = 0; continue; }
        if (l[i] == '\'' || l[i] == '"') { q = l[i]; continue; }
        if (l[i] == '-' && l[i + 1] == '-') return l.substr(0, i);
    }
    return l;
}

class Lexer {
public:
    RawDeck& D;
    std::vector<Tk> toks;
    std::vector<std::vector<double>> arrays;
    std::map<std::string, std::string> paths;
    std::string rootDir;
    std::string section = "RUNSPEC";
    bool stop = false;

    explicit Lexer(RawDeck& d) : D(d) {}

    void err(Level l, int file, int line, const std::string& kw, const std::string& text) { D.messages.push_back(Message{l, kw, file, line, text}); }

    // keyword line?  -> name; known/candidate flags
    bool keywordLine(const std::string& code, std::string& name, const KeywordSpec*& spec) {
        size_t a = 0;
        while (a < code.size() && code[a] == ' ') ++a;
        if (a >= code.size()) return false;
        size_t e = a;
        while (e < code.size() && !std::isspace((unsigned char)code[e])) ++e;
        name = code.substr(a, e - a);
        size_t t = e; while (t < code.size() && std::isspace((unsigned char)code[t])) ++t;
        spec = nullptr;
        if (!name.empty() && std::isupper((unsigned char)name[0])) {
            spec = findKeyword(name);
            if (spec && std::strcmp(spec->name, "(summary)") == 0) {     // families only where the database allows them
                int bit = -1;
                for (int s = 0; s < 8; ++s) if (section == kSections[s]) bit = s;
                if (bit < 0 || !(spec->sections & (1u << bit))) spec = nullptr;
            }
        }
        if (t < code.size()) {
            // ECLIPSE reads a keyword from columns 1-8 and ignores the rest of the line. Accept trailing text for
            // section headers ("SOLUTION =====") and '=' decorations, and - for known keywords in column 1 - when
            // no record is open (the previous token closed a record or was a keyword without data). Old decks
            // carry generator stamps there, e.g. "PVDG       2 TABLES   20 P NODES ...   FIELD   12:13  2 JAN 90".
            // Inside an open record (e.g. RPTSOL mnemonics "PRESSURE SWAT ...") such a line stays data.
            bool isSection = false;
            for (auto* sn : kSections) if (name == sn) isSection = true;
            bool ok = a == 0 && (isSection || code[t] == '=');
            if (!ok && a == 0 && spec && name.size() <= 8 && (!arrayMode || toks.empty() || curKw < 0 || toks[curKw].arrayClosed) && recordsClosed()) ok = true;
            // ...but not inside a "records until an empty record" list that is still open when the line itself is a
            // complete record: EQUALS / COPY / ADD records start with an array name ("SWL 0.05 1 46 1 112 2 2 /")
            if (ok && !isSection && code[t] != '=' && slashListOpen() && code.find('/', t) != std::string::npos) ok = false;
            if (!ok) { spec = nullptr; return false; }
        }
        if (spec) return true;
        return a == 0 && keywordLike(name);          // unknown keywords must start in column 1
    }

    // true if a new keyword may start here: no tokens yet, the last token closed a record, or the last token is a
    // keyword that takes no data
    bool slashListOpen() const {
        if (curKw < 0 || toks.empty() || size_t(curKw) >= toks.size()) return false;
        const KeywordSpec* ks = findKeyword(toks[size_t(curKw)].s);
        if (!ks || ks->kind != KwSize::Slash) return false;
        const size_t n = toks.size();
        if (n - 1 == size_t(curKw)) return true;                         // keyword without records yet
        if (toks[n - 1].t != Tk::Slash) return true;                     // record still open
        return !(toks[n - 2].t == Tk::Slash || n - 2 == size_t(curKw));  // empty record terminates the list
    }
    bool recordsClosed() const {
        if (toks.empty()) return true;
        const Tk& b = toks.back();
        if (b.t == Tk::Slash) return true;
        if (b.t == Tk::Kw) {
            if (b.array >= 0) return b.arrayClosed;          // array keyword terminated by '/'
            const KeywordSpec* ks = findKeyword(b.s);
            return ks && ks->kind == KwSize::None;
        }
        return false;
    }

    std::string resolve(std::string f, const std::string& dir) {
        for (auto& [k, v] : paths) {
            const std::string key = "$" + k;
            size_t pos;
            while ((pos = f.find(key)) != std::string::npos) f.replace(pos, key.size(), v);
        }
        for (char& ch : f) if (ch == '\\') ch = '/';
        namespace fs = std::filesystem;
        auto exists = [](const std::string& p) { std::error_code ec; return fs::is_regular_file(fs::u8path(p), ec); };
        std::vector<std::string> cands;
        const bool abs = !f.empty() && (f[0] == '/' || (f.size() > 1 && f[1] == ':'));
        if (abs) cands.push_back(f);
        else { cands.push_back(dir + f); if (rootDir != dir) cands.push_back(rootDir + f); }
        for (auto& c : cands) if (exists(c)) return c;
        // case-insensitive fallback (decks written on Windows, read on case-sensitive file systems)
        for (auto& c : cands) {
            fs::path p = fs::u8path(c);
            fs::path acc;
            bool ok = true;
            for (auto& part : p) {
                const fs::path next = acc / part;
                std::error_code ec;
                if (acc.empty() || fs::exists(next, ec)) { acc = next; continue; }
                bool found = false;
                for (auto& de : fs::directory_iterator(acc, ec)) {
                    if (upper(de.path().filename().u8string()) == upper(part.u8string())) { acc = de.path(); found = true; break; }
                }
                if (!found) { ok = false; break; }
            }
            if (ok && exists(acc.u8string())) return acc.u8string();
        }
        return {};
    }

    void splitData(const std::string& l, int fi, int ln) {
        size_t i = 0;
        while (i < l.size()) {
            // ',' separates items like blanks (old fixed-format decks: 'LU1     ','ONE     ', 14,  8,7300.00,'OIL')
            if (std::isspace((unsigned char)l[i]) || l[i] == ',') { ++i; continue; }
            Tk t; t.file = fi; t.line = ln;
            if (l[i] == '\'' || l[i] == '"') {
                const char q = l[i];
                size_t e = l.find(q, i + 1);
                if (e == std::string::npos) { err(Level::Warning, fi, ln, "", "unterminated quote"); e = l.size(); }
                t.s = l.substr(i + 1, e - i - 1); t.quoted = true;
                // names padded to 8 characters ('I*      ', 'THP ') - trailing blanks are not significant
                while (!t.s.empty() && t.s.back() == ' ') t.s.pop_back();
                toks.push_back(t); i = e + 1; continue;
            }
            if (l[i] == '/') { t.t = Tk::Slash; toks.push_back(t); afterSlash(); return; }   // rest of the line is a comment
            size_t e = i;
            while (e < l.size() && !std::isspace((unsigned char)l[e]) && l[e] != '/' && l[e] != '\'' && l[e] != '"' && l[e] != ',') ++e;
            t.s = l.substr(i, e - i);
            toks.push_back(t); i = e;
        }
    }

    // numeric array of a Data keyword, lexed straight into doubles (ZCORN can hold millions of values)
    bool arrayLine(const std::string& l, int fi, int ln) {
        std::vector<double>& a = arrays[toks[curKw].array];
        size_t i = 0;
        while (i < l.size()) {
            if (std::isspace((unsigned char)l[i]) || l[i] == ',') { ++i; continue; }
            if (l[i] == '/') { toks[curKw].arrayClosed = true; return true; }
            size_t e = i;
            while (e < l.size() && !std::isspace((unsigned char)l[e]) && l[e] != '/' && l[e] != ',') ++e;
            std::string w = l.substr(i, e - i);
            i = e;
            if (!w.empty() && (w[0] == '\'' || w[0] == '"')) w = w.substr(1, w.size() > 1 ? w.size() - 2 : 0);
            const size_t star = w.find('*');
            long rep = 1;
            if (star != std::string::npos && star > 0 && std::all_of(w.begin(), w.begin() + long(star), ::isdigit)) {
                rep = std::atol(w.substr(0, star).c_str());
                w = w.substr(star + 1);
            }
            double v = std::numeric_limits<double>::quiet_NaN();
            if (!w.empty() && !parseNumber(w, v)) {
                if (!arrayWarned) { err(Level::Warning, fi, ln, toks[curKw].s, "non-numeric value '" + w + "' in array"); arrayWarned = true; }
                v = std::numeric_limits<double>::quiet_NaN();
            }
            if (rep > 0 && rep < 200000000) a.insert(a.end(), size_t(rep), v);
        }
        return false;
    }

    // INCLUDE / PATHS bookkeeping after each '/'
    void afterSlash() {
        if (curKw < 0) return;
        const std::string& k = toks[curKw].s;
        if (k == "INCLUDE" && includeReady) {
            // tokens: [kw] [name...] [/]
            std::string name;
            for (size_t i = size_t(curKw) + 1; i + 1 < toks.size(); ++i) name += toks[i].s;
            const int fi = toks[curKw].file, ln = toks[curKw].line;
            toks.resize(size_t(curKw));
            curKw = -1; includeReady = false;
            pendingInclude.push_back({name, fi, ln});
        } else if (k == "PATHS") {
            // record = tokens since the previous slash
            size_t s = toks.size() - 1;
            size_t b = s;
            while (b > size_t(curKw) + 1 && toks[b - 1].t != Tk::Slash) --b;
            if (s - b >= 2) {
                std::string a = toks[b].s, v = toks[b + 1].s;
                if (!v.empty() && v.back() != '/') v += '/';
                paths[a] = v;
            }
        }
    }

    struct Inc { std::string name; int file, line; };
    std::vector<Inc> pendingInclude;
    long curKw = -1;            // index of the current keyword token (-1 after INCLUDE was consumed)
    bool includeReady = false;
    bool arrayMode = false, arrayWarned = false;
    bool udqMode = false;

    void lexFile(const std::string& text, int fi, const std::string& dir, int depth) {
        std::istringstream in(text);
        std::string line;
        int ln = 0;
        int rawMode = 0;          // 1 = TITLE line, 2 = block until end marker
        std::string rawEnd, rawText;
        int rawLine = 0;
        bool skipping = false;
        while (!stop && std::getline(in, line)) {
            ++ln;
            if (!line.empty() && line.back() == '\r') line.pop_back();
            for (char& ch : line) if (ch == '\t') ch = ' ';
            if (skipping) {
                if (upper(trim(stripComment(line))) == "ENDSKIP") skipping = false;
                continue;
            }
            if (rawMode == 1) {
                Tk t; t.t = Tk::Raw; t.s = trim(line); t.file = fi; t.line = ln; toks.push_back(t);
                rawMode = 0; continue;
            }
            if (rawMode == 2) {
                const std::string tl = upper(trim(stripComment(line)));
                if (tl == rawEnd || tl.rfind(rawEnd + " ", 0) == 0) {
                    Tk t; t.t = Tk::Raw; t.s = rawText; t.file = fi; t.line = rawLine; toks.push_back(t);
                    rawMode = 0; rawText.clear(); continue;
                }
                rawText += line; rawText += '\n';
                continue;
            }
            if (line.rfind("-- RESSIM+", 0) == 0) {   // continuation of a long '-- RESSIM:' line (see wrapLongLines)
                if (!toks.empty() && toks.back().t == Tk::Slash && curKw >= 0 && toks[size_t(curKw)].s == "__RESSIM") {
                    const Tk sl = toks.back(); toks.pop_back();
                    splitData(line.substr(10), fi, ln);
                    toks.push_back(sl);
                }
                continue;
            }
            if (line.rfind("-- RESSIM:", 0) == 0) {
                Tk k; k.t = Tk::Kw; k.s = "__RESSIM"; k.file = fi; k.line = ln; toks.push_back(k); curKw = long(toks.size()) - 1;
                arrayMode = udqMode = false;
                splitData(line.substr(10), fi, ln);
                Tk sl; sl.t = Tk::Slash; sl.file = fi; sl.line = ln; toks.push_back(sl);
                continue;
            }
            const std::string code = stripComment(line);
            std::string name; const KeywordSpec* spec = nullptr;
            if (keywordLine(code, name, spec)) {
                if (arrayMode && !spec) { /* unknown word inside an array: treat as data */ }
                else {
                    if (arrayMode && !toks[curKw].arrayClosed)
                        err(Level::Warning, toks[curKw].file, toks[curKw].line, toks[curKw].s, "array not terminated by '/'");
                    arrayMode = false; udqMode = false;
                    if (name == "END") { stop = true; break; }
                    if (name == "ENDINC") break;
                    if (name == "SKIP" || name == "SKIP100") { skipping = true; continue; }
                    if (name == "SKIP300" || name == "ENDSKIP") continue;
                    for (int s = 0; s < 8; ++s) if (name == kSections[s]) section = name;
                    Tk k; k.t = Tk::Kw; k.s = name; k.file = fi; k.line = ln; k.known = spec != nullptr; k.candidate = !spec;
                    toks.push_back(k); curKw = long(toks.size()) - 1;
                    includeReady = name == "INCLUDE";
                    if (spec) {
                        if (spec->kind == KwSize::Raw) rawMode = 1;
                        else if (spec->kind == KwSize::Code || spec->kind == KwSize::Action) {
                            rawMode = 2; rawEnd = (spec->end && *spec->end) ? spec->end : "ENDACTIO"; rawLine = ln + 1;
                        } else if (spec->kind == KwSize::Data && section != "SUMMARY") {
                            toks[curKw].array = int(arrays.size()); arrays.emplace_back(); arrayMode = true; arrayWarned = false;
                        }
                        udqMode = name == "UDQ" || name == "UDT" || name == "UDADIMS" ? (name == "UDQ") : false;
                    }
                    flushIncludes(dir, depth);
                    continue;
                }
            }
            if (arrayMode) {
                if (arrayLine(code, fi, ln)) arrayMode = false;
                continue;
            }
            if (udqMode) {      // UDQ records contain '/' as division: the LAST '/' on a line terminates the record
                const std::string tc = trim(code);
                if (tc.empty()) continue;
                size_t last = std::string::npos; char q = 0;
                for (size_t i = 0; i < tc.size(); ++i) {
                    if (q) { if (tc[i] == q) q = 0; continue; }
                    if (tc[i] == '\'' || tc[i] == '"') q = tc[i]; else if (tc[i] == '/') last = i;
                }
                if (last == std::string::npos) { Tk t; t.s = tc; t.quoted = true; t.file = fi; t.line = ln; toks.push_back(t); continue; }
                const std::string body = trim(tc.substr(0, last));
                if (!body.empty()) { Tk t; t.s = body; t.quoted = true; t.file = fi; t.line = ln; toks.push_back(t); }
                Tk sl; sl.t = Tk::Slash; sl.file = fi; sl.line = ln; toks.push_back(sl);
                continue;
            }
            splitData(code, fi, ln);
            flushIncludes(dir, depth);
        }
        if (rawMode == 2) {
            err(Level::Warning, fi, rawLine, "", "block not closed by " + rawEnd);
            Tk t; t.t = Tk::Raw; t.s = rawText; t.file = fi; t.line = rawLine; toks.push_back(t);
        }
        if (arrayMode && depth == 0 && curKw >= 0 && !toks[curKw].arrayClosed)
            err(Level::Warning, toks[curKw].file, toks[curKw].line, toks[curKw].s, "array not terminated by '/'");
    }

    void flushIncludes(const std::string& dir, int depth) {
        while (!pendingInclude.empty()) {
            Inc inc = pendingInclude.front();
            pendingInclude.erase(pendingInclude.begin());
            const std::string path = resolve(inc.name, dir);
            if (path.empty()) { err(Level::Error, inc.file, inc.line, "INCLUDE", "cannot open include file '" + inc.name + "'"); continue; }
            if (depth > 20) { err(Level::Error, inc.file, inc.line, "INCLUDE", "INCLUDE nested too deep"); continue; }
            std::ifstream is(std::filesystem::u8path(path), std::ios::binary);
            if (!is) { err(Level::Error, inc.file, inc.line, "INCLUDE", "cannot read include file '" + path + "'"); continue; }
            std::stringstream ss; ss << is.rdbuf();
            D.files.push_back(path);
            const int nfi = int(D.files.size()) - 1;
            std::string ndir;
            const size_t s = path.find_last_of('/');
            if (s != std::string::npos) ndir = path.substr(0, s + 1);
            const bool am = arrayMode, um = udqMode; const long ck = curKw;
            arrayMode = udqMode = false; curKw = -1;
            lexFile(ss.str(), nfi, ndir, depth + 1);
            // an included file may continue the array / records of the keyword that preceded INCLUDE
            if (curKw < 0) { curKw = ck; arrayMode = am; udqMode = um; }
        }
    }
};

// ============================================================================ structural parser
class Structurer {
public:
    RawDeck& D;
    Lexer& L;
    size_t p = 0;
    std::map<std::string, Record> dims;   // first record of every keyword that sizes others
    std::string rockOpt;
    std::string section = "RUNSPEC";
    Structurer(RawDeck& d, Lexer& l) : D(d), L(l) {}

    const Tk& cur() const { return L.toks[p]; }
    bool done() const { return p >= L.toks.size(); }
    bool atKeyword() const { return !done() && cur().t == Tk::Kw && (cur().known || !inRecords); }
    bool inRecords = false;

    // reads one record; returns false if it ended at a keyword / end of deck without '/'
    bool record(Record& r) {
        r = Record();
        r.line = done() ? 0 : cur().line;
        while (!done()) {
            const Tk& t = cur();
            if (t.t == Tk::Kw && (t.known || t.s == "__RESSIM")) return false;
            if (t.t == Tk::Raw) return false;
            ++p;
            if (t.t == Tk::Slash) return true;
            if (t.t == Tk::Kw) { r.items.push_back({t.s, false, false}); continue; }   // unknown word on its own line -> data
            if (!t.quoted) {
                const size_t star = t.s.find('*');
                if (star != std::string::npos && (star == 0 || std::all_of(t.s.begin(), t.s.begin() + long(star), ::isdigit))) {
                    const int n = star == 0 ? 1 : std::atoi(t.s.substr(0, star).c_str());
                    const std::string v = t.s.substr(star + 1);
                    for (int k = 0; k < n && k < 10000000; ++k) r.items.push_back(v.empty() ? Item{"", true, false} : Item{v, false, false});
                    continue;
                }
            }
            r.items.push_back({t.s, false, t.quoted});
        }
        return false;
    }

    int sizeFrom(const KeywordSpec& s, const std::string& name) {
        std::string kw = s.sizeKw ? s.sizeKw : "";
        int item = s.sizeItem, def = s.sizeDefault;
        if (name == "ROCK" && !rockOpt.empty()) {
            if (rockOpt == "SATNUM") { kw = "TABDIMS"; item = 0; }
            else if (rockOpt == "ROCKNUM") {
                auto it = dims.find("TABDIMS");
                if (it != dims.end() && it->second.has(12)) { kw = "TABDIMS"; item = 12; }
            }
        }
        auto it = dims.find(kw);
        int n = def;
        if (it != dims.end() && item >= 0 && it->second.has(size_t(item))) {
            try { n = it->second.integer(size_t(item), def); } catch (...) {}
        }
        // keywords whose tables consist of several records each (OPM "records" definitions)
        static const char* twoRec[] = {"PLYSHLOG", "PLYMWINJ", "SKPRWAT", "SKPRPOLY"};
        int per = 1;
        for (const char* k2 : twoRec) if (name == k2) per = 2;
        return std::max(0, (n + s.shift) * per);
    }

    void run() {
        while (!done()) {
            const Tk& t = cur();
            if (t.t != Tk::Kw) {
                // data outside any keyword: report once per run of tokens, then skip
                const int fi = t.file, ln = t.line;
                std::string what = t.t == Tk::Slash ? "/" : t.s;
                while (!done() && cur().t != Tk::Kw) ++p;
                D.messages.push_back(Message{Level::Warning, "", fi, ln, "unexpected data '" + what + "' outside a keyword - skipped"});
                continue;
            }
            Keyword K;
            K.name = t.s; K.file = t.file; K.line = t.line; K.known = t.known || t.s == "__RESSIM";
            for (int sct = 0; sct < 8; ++sct) if (K.name == kSections[sct]) section = K.name;
            const size_t msgBefore = D.messages.size();
            const int arr = t.array;
            const bool arrClosed = t.arrayClosed;
            ++p;
            const KeywordSpec* spec = (K.name == "__RESSIM") ? nullptr : findKeyword(K.name);
            inRecords = true;
            Record r;
            auto readN = [&](int n, const char* what) {
                for (int i = 0; i < n; ++i) {
                    if (done() || (cur().t == Tk::Kw && (cur().known || cur().s == "__RESSIM"))) {
                        if (i < n) D.message(Level::Warning, K, "expected " + std::to_string(n) + " " + what + " but found " + std::to_string(i) +
                                                                   (spec && spec->kind == KwSize::FromKw ? " (record count set by " + std::string(spec->sizeKw) + ")" : ""));
                        return;
                    }
                    if (!record(r)) { D.message(Level::Warning, K, "record not terminated by '/'", r.line); if (!r.items.empty()) K.records.push_back(r); return; }
                    K.records.push_back(r);
                }
            };
            if (K.name == "__RESSIM") {
                if (record(r)) K.records.push_back(r);
            } else if (!spec) {
                while (!done() && !(cur().t == Tk::Kw)) { if (record(r) || !r.items.empty()) K.records.push_back(r); else break; }
            } else {
                switch (spec->kind) {
                case KwSize::None: break;
                case KwSize::Fixed: readN(spec->count, spec->count == 1 ? "record" : "records"); break;
                case KwSize::FromKw: readN(sizeFrom(*spec, K.name), "records"); break;
                case KwSize::Slash:
                    while (true) {
                        if (done() || (cur().t == Tk::Kw && cur().known)) { D.message(done() ? Level::Warning : Level::Info, K, "list not terminated by an empty '/' record (accepted: the next keyword follows)"); break; }
                        const bool ok = record(r);
                        if (ok && r.items.empty()) break;
                        K.records.push_back(r);
                        if (!ok) { D.message(Level::Warning, K, "record not terminated by '/'", r.line); break; }
                    }
                    break;
                case KwSize::Unknown:
                    while (!done() && !(cur().t == Tk::Kw && cur().known)) { const bool ok = record(r); if (!r.items.empty()) K.records.push_back(r); if (!ok) break; }
                    break;
                case KwSize::Tables: {
                    const int n = sizeFrom(*spec, K.name);
                    for (int tb = 0; tb < n; ++tb) {
                        bool any = false;
                        while (true) {
                            if (done() || (cur().t == Tk::Kw && cur().known)) {
                                if (tb < n) D.message(Level::Warning, K, "expected " + std::to_string(n) + " tables (" + std::string(spec->sizeKw) + ") but found " + std::to_string(tb + (any ? 1 : 0)));
                                tb = n; break;
                            }
                            const bool ok = record(r);
                            if (ok && r.items.empty()) { Record term; term.terminator = true; term.line = r.line; K.records.push_back(term); break; }
                            K.records.push_back(r); any = true;
                            if (!ok) break;
                        }
                    }
                    break;
                }
                case KwSize::Data:
                    if (arr >= 0) { K.values = std::move(L.arrays[size_t(arr)]); K.isArray = true; (void)arrClosed; }
                    else readN(1, "record");
                    break;
                case KwSize::Code: case KwSize::Action: case KwSize::Raw:
                    if (!done() && cur().t == Tk::Raw) { K.raw = cur().s; ++p; }
                    if (spec->kind == KwSize::Action) {   // ACTIONX first record lives on the raw lines; nothing else to read
                    }
                    break;
                }
            }
            inRecords = false;
            // trailing data after a complete keyword
            if (!done() && cur().t != Tk::Kw && cur().t != Tk::Raw && spec && spec->kind != KwSize::Unknown) {
                const int fi = cur().file, ln = cur().line;
                size_t q = p; int extra = 0;
                while (q < L.toks.size() && !(L.toks[q].t == Tk::Kw && (L.toks[q].known || L.toks[q].s == "__RESSIM"))) { if (L.toks[q].t == Tk::Slash) ++extra; ++q; }
                // stray empty records ("/" lines) after a complete keyword are harmless (OPM accepts them too)
                bool onlySlashes = true;
                for (size_t z = p; z < q; ++z) if (L.toks[z].t != Tk::Slash) onlySlashes = false;
                if (!onlySlashes)
                    D.messages.push_back(Message{Level::Warning, K.name, fi, ln,
                        "extra data after the keyword's " + std::to_string(K.records.size()) + " record(s) - ignored" +
                        (spec->kind == KwSize::FromKw || spec->kind == KwSize::Tables ? " (check " + std::string(spec->sizeKw) + ")" : "")});
                p = q;
            }
            if (spec && (K.name == "TABDIMS" || K.name == "EQLDIMS" || K.name == "REGDIMS" || K.name == "WELLDIMS" || K.name == "AQUDIMS" ||
                         K.name == "ROCKCOMP" || K.name == "ENDSCALE" || K.name == "VFPPDIMS" || K.name == "VFPIDIMS" || K.name == "NNEWTF" ||
                         K.name == "PIMTDIMS" || K.name == "NUPCOL" || K.name == "SMRYDIMS" || K.name == "ACTDIMS" || K.name == "UDQDIMS" ||
                         K.name == "UDADIMS" || K.name == "WSEGDIMS" || K.name == "TRACERS" || K.name == "MISCIBLE" || K.name == "NETWORK" ||
                         K.name == "PVTWSALT" || K.name == "FAULTDIM" || K.name == "GPTDIMS" || K.name == "NRSOUT" || K.name == "SATOPTS" ||
                         K.name == "ROCKOPTS" || K.name == "EHYSTR" || K.name == "PLMIXPAR" || K.name == "LGRCOPY" || K.name == "GRIDOPTS") &&
                !K.records.empty() && !dims.count(K.name))
                dims[K.name] = K.records.front();
            if (K.name == "ROCKOPTS" && !K.records.empty()) rockOpt = upper(K.records.front().str(2, ""));
            if (section == "SUMMARY")       // output requests only: structure problems there never matter for the model
                for (size_t z = msgBefore; z < D.messages.size(); ++z) D.messages[z].level = Level::Info;
            D.keywords.push_back(std::move(K));
        }
    }
};

}  // namespace

RawDeck parseDeckText(const std::string& text, const std::string& baseDir, const std::string& fileName) {
    RawDeck D;
    D.files.push_back(fileName);
    Lexer L(D);
    L.rootDir = baseDir;
    if (!L.rootDir.empty() && L.rootDir.back() != '/' && L.rootDir.back() != '\\') L.rootDir += '/';
    for (char& ch : L.rootDir) if (ch == '\\') ch = '/';
    L.lexFile(text, 0, L.rootDir, 0);
    Structurer S(D, L);
    S.run();
    return D;
}

RawDeck parseDeckFile(const std::string& path) {
    std::ifstream f(std::filesystem::u8path(path), std::ios::binary);
    if (!f) throw std::runtime_error("Cannot open " + path);
    std::stringstream ss; ss << f.rdbuf();
    std::string p = path;
    for (char& ch : p) if (ch == '\\') ch = '/';
    std::string dir;
    const size_t s = p.find_last_of('/');
    if (s != std::string::npos) dir = p.substr(0, s + 1);
    return parseDeckText(ss.str(), dir, path);
}

// ============================================================================ syntax classification
void classifyLine(const std::string& line, std::vector<TokSpan>& out, LineState& st) {
    out.clear();
    const int n = int(line.size());
    if (n == 0) return;
    if (st.raw == 1) { out.push_back({0, n, Tok::String}); st.raw = 0; return; }
    if (st.raw == 2) {
        std::string t = upper(trim(stripComment(line)));
        if (t == st.end || t.rfind(st.end + " ", 0) == 0) { st.raw = 0; out.push_back({0, n, Tok::Keyword}); return; }
        out.push_back({0, n, Tok::Text});
        return;
    }
    if (line.rfind("-- RESSIM:", 0) == 0 || line.rfind("-- RESSIM+", 0) == 0) { out.push_back({0, n, Tok::RessimExt}); return; }
    // keyword line?
    {
        int a = 0; while (a < n && line[a] == ' ') ++a;
        int e = a; while (e < n && !std::isspace((unsigned char)line[e])) ++e;
        int r = e; while (r < n && std::isspace((unsigned char)line[r])) ++r;
        const bool restComment = r >= n || (r + 1 < n && line[r] == '-' && line[r + 1] == '-');
        // trailing text after a column-1 keyword is ignored by ECLIPSE (generator stamps in old decks)
        const bool stamped = !restComment && a == 0 && e - a <= 8 && st.closed;
        if (e > a && (restComment || stamped) && std::isupper((unsigned char)line[a])) {
            const std::string w = line.substr(a, e - a);
            const KeywordSpec* s = findKeyword(w);
            bool section = false;
            for (auto* sn : kSections) if (w == sn) section = true;
            const bool accept = restComment ? (s || (a == 0 && keywordLike(w) && st.closed)) : (s || section);
            if (accept) {
                st.closed = section || (s && s->kind == KwSize::None);
                if (a > 0) out.push_back({0, a, Tok::Space});
                out.push_back({a, e, section ? Tok::Section : (s ? Tok::Keyword : Tok::KeywordUnknown)});
                if (r < n) out.push_back({r, n, Tok::Comment});
                if (s && s->kind == KwSize::Raw) st.raw = 1;
                if (s && (s->kind == KwSize::Code || s->kind == KwSize::Action)) { st.raw = 2; st.end = (s->end && *s->end) ? s->end : "ENDACTIO"; }
                return;
            }
        }
    }
    {   // record state for the next line
        std::string c = trim(stripComment(line));
        if (!c.empty()) {
            bool sl = false; char q = 0;
            for (char ch : c) { if (q) { if (ch == q) q = 0; continue; } if (ch == '\'' || ch == '"') q = ch; else if (ch == '/') { sl = true; break; } }
            st.closed = sl;
        }
    }
    int i = 0;
    while (i < n) {
        const char c = line[i];
        if (std::isspace((unsigned char)c)) { int e = i; while (e < n && std::isspace((unsigned char)line[e])) ++e; out.push_back({i, e, Tok::Space}); i = e; continue; }
        if (c == '-' && i + 1 < n && line[i + 1] == '-') { out.push_back({i, n, Tok::Comment}); return; }
        if (c == '\'' || c == '"') {
            int e = int(line.find(c, size_t(i) + 1)); e = e < 0 ? n : e + 1;
            out.push_back({i, e, Tok::String}); i = e; continue;
        }
        if (c == '/') { out.push_back({i, i + 1, Tok::Slash}); if (i + 1 < n) out.push_back({i + 1, n, Tok::Comment}); return; }
        int e = i;
        while (e < n && !std::isspace((unsigned char)line[e]) && line[e] != '/' && line[e] != '\'' && line[e] != '"') ++e;
        const std::string w = line.substr(size_t(i), size_t(e - i));
        double v;
        const size_t star = w.find('*');
        Tok k = Tok::Text;
        if (star != std::string::npos && std::all_of(w.begin(), w.begin() + long(star), ::isdigit)) k = (star + 1 == w.size()) ? Tok::Default : Tok::Number;
        else if (parseNumber(w, v)) k = Tok::Number;
        out.push_back({i, e, k});
        i = e;
    }
}

std::string keywordBodyText(const Keyword& k) {
    std::ostringstream o;
    const KeywordSpec* sp = findKeyword(k.name);
    auto item = [](const Item& it) -> std::string {
        if (it.def) return "1*";
        if (it.quoted || it.s.empty() || it.s.find_first_of(" \t/") != std::string::npos) return "'" + it.s + "'";
        return it.s;
    };
    if (sp && (sp->kind == KwSize::Code || sp->kind == KwSize::Action || sp->kind == KwSize::Raw)) {
        o << k.raw;
        if (!k.raw.empty() && k.raw.back() != '\n') o << '\n';
        if (sp->kind != KwSize::Raw) o << ((sp->end && *sp->end) ? sp->end : "ENDACTIO") << '\n';
        return o.str();
    }
    if (k.isArray) {
        int col = 0; size_t i = 0;
        o << ' ';
        while (i < k.values.size()) {
            size_t j = i + 1;
            const bool nan = std::isnan(k.values[i]);
            while (j < k.values.size() && (nan ? std::isnan(k.values[j]) : k.values[j] == k.values[i])) ++j;
            const size_t rep = j - i;
            char b[48];
            if (nan) std::snprintf(b, sizeof(b), "%zu*", rep);
            else if (rep > 1) std::snprintf(b, sizeof(b), "%zu*%.10g", rep, k.values[i]);
            else std::snprintf(b, sizeof(b), "%.10g", k.values[i]);
            o << b << ' ';
            if (++col % 8 == 0) o << "\n ";
            i = j;
        }
        o << "/\n";
        return o.str();
    }
    const bool udq = k.name == "UDQ";     // records are whole lines (the line text is one item)
    for (const Record& r : k.records) {
        if (r.terminator) { o << "/\n"; continue; }
        o << ' ';
        for (const Item& it : r.items) o << (udq ? it.s : item(it)) << ' ';
        o << "/\n";
    }
    if (sp && sp->kind == KwSize::Slash) o << "/\n";
    return o.str();
}

}  // namespace deck
}  // namespace ressim
