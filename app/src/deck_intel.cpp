// Deck intelligence: see deck_intel.h.
#include "deck_intel.h"
#include "ressim/deck_parser.h"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <set>
#include <sstream>

namespace ressim {
namespace deckintel {

namespace {

std::vector<std::string> splitTabs(const std::string& s) {
    std::vector<std::string> out; size_t a = 0;
    for (;;) {
        const size_t b = s.find('\t', a);
        out.push_back(s.substr(a, b == std::string::npos ? std::string::npos : b - a));
        if (b == std::string::npos) break;
        a = b + 1;
    }
    return out;
}
std::string unescapeNl(const std::string& s) {
    std::string o; o.reserve(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '\\' && i + 1 < s.size() && s[i + 1] == 'n') { o += '\n'; ++i; } else o += s[i];
    }
    return o;
}
std::vector<std::string> splitChar(const std::string& s, char c) {
    std::vector<std::string> out; std::string cur;
    for (char ch : s) { if (ch == c) { out.push_back(cur); cur.clear(); } else cur += ch; }
    out.push_back(cur);
    return out;
}
bool isWordStart(char c) { return std::isalpha((unsigned char)c) != 0; }
bool isWordChar(char c) { return std::isalnum((unsigned char)c) || c == '_' || c == '+' || c == '-'; }
bool startsWithNoCase(const std::string& s, const std::string& p) {
    if (p.size() > s.size()) return false;
    for (size_t i = 0; i < p.size(); ++i)
        if (std::toupper((unsigned char)s[i]) != std::toupper((unsigned char)p[i])) return false;
    return true;
}
bool isNumber(const std::string& t) {
    if (t.empty()) return false;
    char* e = nullptr;
    std::strtod(t.c_str(), &e);
    if (e && *e == 0) return true;
    // Fortran style exponents: 1.0D+03
    std::string u = t;
    for (char& c : u) if (c == 'd' || c == 'D') c = 'E';
    std::strtod(u.c_str(), &e);
    return e && *e == 0;
}
deck::KwSize kindOf(const std::string& kw) {
    const deck::KeywordSpec* s = deck::findKeyword(kw);
    if (s) return s->kind;
    const KwInfo* d = Docs::get().find(kw);
    if (d) {
        if (d->sizeKind == "none") return deck::KwSize::None;
        if (d->sizeKind == "array") return deck::KwSize::Data;
        if (d->sizeKind == "list") return deck::KwSize::Slash;
        if (d->sizeKind == "fixed") return deck::KwSize::Fixed;
    }
    return deck::KwSize::Unknown;
}
bool isListKind(deck::KwSize k) { return k == deck::KwSize::Slash || k == deck::KwSize::Tables; }
bool isRecordKind(deck::KwSize k) {
    return k == deck::KwSize::Fixed || k == deck::KwSize::FromKw || k == deck::KwSize::Slash || k == deck::KwSize::Unknown;
}
bool cleanDefault(const std::string& d) {
    if (d.empty() || d == "None" || d.size() > 14) return false;
    for (char c : d) if (c == ' ' || c == ',' || c == '(') return false;
    return true;
}

// Body items of a keyword with record boundaries: calls f(token, lineIndex, record, itemBefore) for every
// non-comment token in body order.
template <class F>
void walkBody(const std::vector<std::string>& lines, int first, int last, F&& f) {
    std::vector<Token> toks;
    int record = 0, item = 0;
    for (int l = first; l <= last && l < (int)lines.size(); ++l) {
        tokenize(lines[l], toks);
        for (const Token& t : toks) {
            if (t.kind == TK::Comment) continue;
            if (!f(t, l, record, item)) return;
            if (t.kind == TK::Slash) { ++record; item = 0; } else item += t.count;
        }
    }
}

} // namespace

// ================================================================== KwInfo / Docs
int KwInfo::maxRecord() const {
    int m = 0;
    for (const Param& p : params) m = std::max(m, p.record);
    return m;
}
int KwInfo::recordType(int r) const {
    const int m = maxRecord();
    return m == 0 ? 0 : std::min(r + 1, m);
}
std::vector<const Param*> KwInfo::recordParams(int r) const {
    const int t = recordType(r);
    std::vector<const Param*> out;
    for (const Param& p : params) if (t == 0 || p.record == 0 || p.record == t) out.push_back(&p);
    return out;
}
const Param* KwInfo::param(int r, int item) const {
    if (item < 0) return nullptr;
    if (variadic && !params.empty()) {
        const auto ps = recordParams(r);
        if (ps.empty()) return nullptr;
        return ps[item % ps.size()];
    }
    for (const Param* p : recordParams(r)) if (p->index == item + 1) return p;
    return nullptr;
}
bool KwInfo::validIn(const std::string& section) const {
    if (sections.empty() || section.empty()) return true;
    for (const std::string& s : splitChar(sections, ',')) if (s == section) return true;
    return false;
}

Docs& Docs::get() { static Docs d; return d; }

bool Docs::load(const std::string& path) {
    std::ifstream f(path);
    if (!f) return false;
    return loadStream(f);
}
bool Docs::loadText(std::string_view text) {
    if (text.empty()) return false;
    std::istringstream f{std::string(text)};
    return loadStream(f);
}
bool Docs::loadStream(std::istream& f) {
    map_.clear();
    std::string line; KwInfo* cur = nullptr;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.size() < 2 || line[0] == '#' || line[1] != '\t') continue;
        const auto t = splitTabs(line);
        if (line[0] == 'K' && t.size() >= 2) {
            KwInfo& d = map_[t[1]];
            d.name = t[1];
            if (t.size() > 2) d.sections = t[2];
            if (t.size() > 3) d.alias = t[3];
            if (t.size() > 4) d.summary = t[4];
            if (t.size() > 5) d.sizeKind = t[5];
            if (t.size() > 6 && !t[6].empty())
                for (const std::string& c : splitChar(t[6], ',')) d.recordColumns.push_back(std::atoi(c.c_str()));
            if (t.size() > 7) d.variadic = t[7] == "1";
            cur = &d;
        } else if (line[0] == 'P' && cur && t.size() >= 3) {
            Param p;
            p.index = std::atoi(t[1].c_str());
            p.name = t[2];
            if (t.size() > 3) p.manual = t[3];
            if (t.size() > 4) p.def = t[4] == "None" ? std::string() : t[4];
            if (t.size() > 5) {
                const size_t bar = t[5].find('|');
                p.unitsField = t[5].substr(0, bar);
                if (bar != std::string::npos) p.unitsMetric = t[5].substr(bar + 1);
            }
            if (t.size() > 6 && !t[6].empty()) p.options = splitChar(t[6], ',');
            if (t.size() > 7) p.desc = t[7];
            if (t.size() > 8) p.record = std::atoi(t[8].c_str());
            if (t.size() > 9) p.valueType = t[9];
            if (t.size() > 10) p.dimension = t[10];
            cur->params.push_back(p);
        } else if (line[0] == 'E' && cur && t.size() >= 2) {
            if (!cur->example.empty()) cur->example += "\n\n";
            cur->example += unescapeNl(t[1]);
        }
    }
    ok_ = true;
    return true;
}
const KwInfo* Docs::findRaw(const std::string& name) const {
    auto it = map_.find(name);
    return it == map_.end() ? nullptr : &it->second;
}
const KwInfo* Docs::find(const std::string& name) const {
    const KwInfo* k = findRaw(name);
    for (int guard = 0; k && !k->alias.empty() && guard < 4; ++guard) {
        const KwInfo* a = findRaw(k->alias);
        if (!a) break;
        k = a;
    }
    return k;
}

// ================================================================== lexical
std::string upper(std::string s) { for (char& c : s) c = (char)std::toupper((unsigned char)c); return s; }
std::string unquote(const std::string& s) {
    if (s.size() >= 2 && (s[0] == '\'' || s[0] == '"') && s.back() == s[0]) return s.substr(1, s.size() - 2);
    return s;
}
const std::vector<std::string>& sectionNames() {
    static const std::vector<std::string> s = {"RUNSPEC", "GRID", "EDIT", "PROPS", "REGIONS", "SOLUTION", "SUMMARY", "SCHEDULE"};
    return s;
}
bool isSectionName(const std::string& n) {
    for (const auto& s : sectionNames()) if (s == n) return true;
    return false;
}

void tokenize(const std::string& line, std::vector<Token>& out) {
    out.clear();
    const int n = (int)line.size();
    int i = 0;
    while (i < n) {
        const char c = line[i];
        if (c == ' ' || c == '\t' || c == '\r' || c == ',') { ++i; continue; }
        Token t; t.begin = i;
        if (c == '-' && i + 1 < n && line[i + 1] == '-') {
            t.kind = TK::Comment; t.end = n; t.text = line.substr(i);
            out.push_back(t); return;
        }
        if (c == '/') {
            t.kind = TK::Slash; t.end = i + 1; t.text = "/";
            out.push_back(t);
            if (i + 1 < n) { Token r; r.kind = TK::Comment; r.begin = i + 1; r.end = n; r.text = line.substr(i + 1); out.push_back(r); }
            return;
        }
        if (c == '\'' || c == '"') {
            int j = i + 1;
            while (j < n && line[j] != c) ++j;
            t.end = std::min(n, j + 1); t.kind = TK::String; t.text = line.substr(i, t.end - i);
            out.push_back(t); i = t.end; continue;
        }
        int j = i;
        while (j < n && line[j] != ' ' && line[j] != '\t' && line[j] != '\r' && line[j] != ',' && line[j] != '/' && line[j] != '\'' &&
               line[j] != '"' && !(line[j] == '-' && j + 1 < n && line[j + 1] == '-'))
            ++j;
        // n*'quoted string'
        if (j < n && (line[j] == '\'' || line[j] == '"') && j > i && line[j - 1] == '*') {
            const char q = line[j];
            int k = j + 1;
            while (k < n && line[k] != q) ++k;
            j = std::min(n, k + 1);
        }
        t.end = j; t.text = line.substr(i, j - i);
        const size_t star = t.text.find('*');
        if (star != std::string::npos && star > 0 && std::all_of(t.text.begin(), t.text.begin() + star, [](char ch) { return std::isdigit((unsigned char)ch); })) {
            t.count = std::max(1, std::atoi(t.text.substr(0, star).c_str()));
            t.kind = star + 1 == t.text.size() ? TK::Defaults : TK::Repeat;
        } else if (t.text == "*") {
            t.kind = TK::Defaults;
        } else if (isNumber(t.text)) {
            t.kind = TK::Number;
        } else {
            t.kind = TK::Word;
        }
        out.push_back(t);
        i = j;
    }
}

bool keywordLine(const std::string& line, std::string* name) {
    if (line.empty() || !isWordStart(line[0])) return false;
    size_t j = 0;
    while (j < line.size() && isWordChar(line[j]) && !(line[j] == '-' && j + 1 < line.size() && line[j + 1] == '-')) ++j;
    if (j > 8) return false;
    size_t k = j;
    while (k < line.size() && (line[k] == ' ' || line[k] == '\t' || line[k] == '\r')) ++k;
    if (k < line.size() && !(line[k] == '-' && k + 1 < line.size() && line[k + 1] == '-')) {
        // ECLIPSE ignores text after a keyword: same rule as the engine parser - accepted for section headers
        // ("PROPS =========") and for decorations that start with '='; otherwise the line is data
        if (j == 0 || k == j || !(line[k] == '=' || isSectionName(upper(line.substr(0, j))))) return false;
    }
    if (name) *name = line.substr(0, j);
    return true;
}

// ================================================================== outline
Outline buildOutline(const std::vector<std::string>& lines) {
    Outline o;
    std::string name;
    std::string rawEnd;   // inside ACTIONX / PYINPUT blocks: keyword lines are still keywords (ACTIONX) or raw code
    bool inCode = false;
    for (int l = 0; l < (int)lines.size(); ++l) {
        if (inCode) {
            if (lines[l].compare(0, rawEnd.size(), rawEnd) == 0) { inCode = false; o.items.push_back({l, rawEnd, false}); }
            continue;
        }
        if (!keywordLine(lines[l], &name)) continue;
        const std::string u = upper(name);
        o.items.push_back({l, name, isSectionName(u)});
        const deck::KeywordSpec* s = deck::findKeyword(u);
        if (s && s->kind == deck::KwSize::Code && s->end) { inCode = true; rawEnd = s->end; }
    }
    return o;
}
int Outline::itemAt(int line) const {
    auto it = std::upper_bound(items.begin(), items.end(), line, [](int l, const OutlineItem& a) { return l < a.line; });
    return int(it - items.begin()) - 1;
}
std::string Outline::sectionAt(int line) const {
    for (int i = itemAt(line); i >= 0; --i) if (items[i].section) return upper(items[i].name);
    return std::string();
}
int Outline::bodyEnd(int item, int lineCount) const {
    if (item < 0) return -1;
    return item + 1 < (int)items.size() ? items[item + 1].line - 1 : lineCount - 1;
}

// ================================================================== context
Context contextAt(const std::vector<std::string>& lines, const Outline& outline, int line, int col) {
    Context c;
    if (line < 0 || line >= (int)lines.size()) return c;
    c.section = outline.sectionAt(line);
    // comment at cursor?
    {
        std::vector<Token> toks; tokenize(lines[line], toks);
        for (const Token& t : toks) if (t.kind == TK::Comment && col >= t.begin + (t.text[0] == '-' ? 1 : 0)) c.inComment = true;
    }
    c.item = outline.itemAt(line);
    if (c.item < 0) return c;
    const OutlineItem& it = outline.items[c.item];
    c.kwLine = it.line;
    c.keyword = upper(it.name);
    c.info = Docs::get().find(c.keyword);
    if (line == it.line) { c.onKeywordLine = true; return c; }
    const deck::KwSize kind = kindOf(c.keyword);
    c.arrayKeyword = kind == deck::KwSize::Data;
    int record = 0, item = 0, found = -1;
    bool after = false;
    walkBody(lines, it.line + 1, line, [&](const Token& t, int l, int rec, int itemBefore) {
        record = rec; item = itemBefore;
        if (l == line && col <= t.end) {
            if (col >= t.begin) { found = itemBefore; if (t.kind == TK::Slash && col > t.begin) after = true; }
            return false;
        }
        if (t.kind == TK::Slash) { record = rec + 1; item = 0; } else item = itemBefore + t.count;
        return true;
    });
    if (after) { c.inComment = true; return c; }
    c.record = record;
    c.itemIndex = found >= 0 ? found : item;
    if (c.arrayKeyword) { c.arrayIndex = c.itemIndex; c.record = 0; }
    if (kind == deck::KwSize::Tables && c.info && c.info->params.size() >= 2) {
        // PVTO / PVTG style: item 0 is the outer value, then repeating rows of the remaining columns
        const int inner = (int)c.info->params.size() - 1;
        c.param = c.itemIndex == 0 ? &c.info->params[0] : &c.info->params[1 + (c.itemIndex - 1) % inner];
    } else if (c.info && !c.arrayKeyword) {
        c.param = c.info->param(c.record, c.itemIndex);
    } else if (c.info && c.arrayKeyword && !c.info->params.empty()) {
        c.param = &c.info->params[0];
    }
    return c;
}

// ================================================================== names
std::vector<std::string> wellNames(const std::vector<std::string>& lines, const Outline& o) {
    std::set<std::string> s;
    for (size_t i = 0; i < o.items.size(); ++i) {
        const std::string u = upper(o.items[i].name);
        if (u != "WELSPECS" && u != "WELSPECL") continue;
        walkBody(lines, o.items[i].line + 1, o.bodyEnd((int)i, (int)lines.size()), [&](const Token& t, int, int, int item) {
            if (item == 0 && t.kind != TK::Slash) s.insert(unquote(t.text));
            return true;
        });
    }
    return std::vector<std::string>(s.begin(), s.end());
}
std::vector<std::string> groupNames(const std::vector<std::string>& lines, const Outline& o) {
    std::set<std::string> s = {"FIELD"};
    for (size_t i = 0; i < o.items.size(); ++i) {
        const std::string u = upper(o.items[i].name);
        const bool spec = u == "WELSPECS" || u == "WELSPECL", tree = u == "GRUPTREE";
        if (!spec && !tree) continue;
        walkBody(lines, o.items[i].line + 1, o.bodyEnd((int)i, (int)lines.size()), [&](const Token& t, int, int, int item) {
            if (t.kind == TK::Slash || t.kind == TK::Defaults) return true;
            if ((spec && item == 1) || (tree && item <= 1)) s.insert(unquote(t.text));
            return true;
        });
    }
    return std::vector<std::string>(s.begin(), s.end());
}
int dimensCells(const std::vector<std::string>& lines, const Outline& o, int* nx, int* ny, int* nz) {
    for (size_t i = 0; i < o.items.size(); ++i) {
        const std::string u = upper(o.items[i].name);
        if (u != "DIMENS" && u != "SPECGRID") continue;
        int d[3] = {0, 0, 0};
        walkBody(lines, o.items[i].line + 1, o.bodyEnd((int)i, (int)lines.size()), [&](const Token& t, int, int rec, int item) {
            if (rec > 0 || t.kind == TK::Slash) return false;
            if (item < 3 && t.kind == TK::Number) d[item] = std::atoi(t.text.c_str());
            return item < 3;
        });
        if (nx) *nx = d[0];
        if (ny) *ny = d[1];
        if (nz) *nz = d[2];
        return d[0] * d[1] * d[2];
    }
    return 0;
}
std::map<std::string, std::string> pathAliases(const std::vector<std::string>& lines, const Outline& o) {
    std::map<std::string, std::string> m;
    for (size_t i = 0; i < o.items.size(); ++i) {
        if (upper(o.items[i].name) != "PATHS") continue;
        std::string alias;
        walkBody(lines, o.items[i].line + 1, o.bodyEnd((int)i, (int)lines.size()), [&](const Token& t, int, int, int item) {
            if (t.kind == TK::Slash) return true;
            if (item == 0) alias = unquote(t.text); else if (item == 1) m[alias] = unquote(t.text);
            return true;
        });
    }
    return m;
}
std::string includePath(const std::vector<std::string>& lines, const Outline& o, int line) {
    const int i = o.itemAt(line);
    if (i < 0) return std::string();
    const std::string u = upper(o.items[i].name);
    if (u != "INCLUDE" && u != "IMPORT" && u != "GDFILE" && u != "RESTART" && u != "LOAD") return std::string();
    std::string p;
    walkBody(lines, o.items[i].line + 1, o.bodyEnd(i, (int)lines.size()), [&](const Token& t, int, int, int) {
        if (t.kind != TK::Slash) p = unquote(t.text);
        return false;
    });
    if (!p.empty() && p[0] == '$') {
        const auto al = pathAliases(lines, o);
        const size_t sl = p.find_first_of("/\\");
        const std::string a = p.substr(1, sl == std::string::npos ? std::string::npos : sl - 1);
        auto it = al.find(a);
        if (it != al.end()) p = it->second + (sl == std::string::npos ? std::string() : p.substr(sl));
    }
    return p;
}

// ================================================================== diagnostics
std::vector<Diag> lint(const std::vector<std::string>& lines, const Outline& o) {
    std::vector<Diag> out;
    const Docs& docs = Docs::get();
    int nx = 0, ny = 0, nz = 0;
    const int cells = dimensCells(lines, o, &nx, &ny, &nz);
    bool box = false;
    std::vector<Token> toks;
    // indented keywords
    std::string nm;
    for (int l = 0; l < (int)lines.size(); ++l) {
        const std::string& s = lines[l];
        if (s.empty() || (s[0] != ' ' && s[0] != '\t')) continue;
        size_t a = s.find_first_not_of(" \t");
        if (a == std::string::npos) continue;
        if (!keywordLine(s.substr(a), &nm)) continue;
        const std::string u = upper(nm);
        if (u.size() >= 3 && (isSectionName(u) || deck::findKeyword(u)) && kindOf(u) != deck::KwSize::Unknown && docs.find(u))
            out.push_back({l, (int)a, (int)(a + nm.size()), 0, u + " looks like a keyword but does not start in column 1 - it is read as data"});
    }
    std::string section;
    for (size_t i = 0; i < o.items.size(); ++i) {
        const OutlineItem& it = o.items[i];
        const std::string u = upper(it.name);
        const int nameEnd = (int)it.name.size();
        if (it.section) { section = u; continue; }
        if (u != it.name) out.push_back({it.line, 0, nameEnd, 1, "keywords must be upper case (" + u + ")"});
        const KwInfo* info = docs.find(u);
        const deck::KeywordSpec* spec = deck::findKeyword(u);
        if (!info && !spec) {
            out.push_back({it.line, 0, nameEnd, 1, "unknown keyword " + u});
            continue;
        }
        if (info && !section.empty() && !info->validIn(section))
            out.push_back({it.line, 0, nameEnd, 1, u + " is not a " + section + " keyword (valid in: " + info->sections + ")"});
        if (u == "BOX") box = true;
        if (u == "ENDBOX") box = false;
        const deck::KwSize kind = kindOf(u);
        const int last = o.bodyEnd((int)i, (int)lines.size());
        if (kind == deck::KwSize::None || kind == deck::KwSize::Code || kind == deck::KwSize::Action || kind == deck::KwSize::Raw) continue;
        // walk the body
        int records = 0, values = 0, lastRecLine = -1, lastItems = 0, lastRecStart = -1;
        bool emptyRecord = false, anyData = false, open = false;
        int curItems = 0, curStartLine = -1, curStartCol = 0;
        std::vector<std::pair<int, int>> recSizes;   // (line, items) per record
        walkBody(lines, it.line + 1, last, [&](const Token& t, int l, int, int itemBefore) {
            anyData = true;
            if (t.kind == TK::Slash) {
                if (itemBefore == 0) emptyRecord = true; else emptyRecord = false;
                recSizes.push_back({curStartLine < 0 ? l : curStartLine, itemBefore});
                (void)curStartCol;
                ++records; lastRecLine = l; lastItems = itemBefore; lastRecStart = curStartLine;
                curStartLine = -1; open = false;
            } else {
                if (curStartLine < 0) { curStartLine = l; curStartCol = t.begin; }
                open = true;
                values += t.count;
            }
            curItems = itemBefore;
            return true;
        });
        (void)lastRecLine; (void)lastItems; (void)lastRecStart; (void)curItems;
        if (kind == deck::KwSize::Data) {
            if (open || records == 0) out.push_back({it.line, 0, nameEnd, 2, u + ": array is not terminated by '/'"});
            else if (cells > 0 && !box && (section == "GRID" || section == "EDIT" || section == "PROPS" || section == "REGIONS" || section == "SOLUTION")) {
                long expect = cells;
                if (u == "COORD") expect = 6L * (nx + 1) * (ny + 1);
                else if (u == "ZCORN") expect = 8L * cells;
                else if (u == "TOPS" || u == "DEPTHZ" || u == "MAPAXES" || u == "COORDSYS" || u == "TVDP") expect = -1;
                if (u == "TOPS" && values != nx * ny && values != cells && values > 0)
                    out.push_back({it.line, 0, nameEnd, 1, u + ": " + std::to_string(values) + " values; expected NX*NY = " + std::to_string(nx * ny) + " or " + std::to_string(cells)});
                if (expect > 0 && values > 0 && values != expect && records == 1)
                    out.push_back({it.line, 0, nameEnd, values > expect ? 2 : 1,
                                   u + ": " + std::to_string(values) + " values, the grid needs " + std::to_string(expect) +
                                       (values < expect ? " (fewer values are only valid inside BOX)" : "")});
            }
            continue;
        }
        if (open) {
            out.push_back({curStartLine, 0, (int)lines[curStartLine].size(), 2, u + ": record is not terminated by '/'"});
        } else if (isListKind(kind) && anyData && !emptyRecord) {
            out.push_back({it.line, 0, nameEnd, 2, u + ": list must be closed by an empty record '/'"});
        } else if (isListKind(kind) && !anyData) {
            out.push_back({it.line, 0, nameEnd, 2, u + ": missing records and the closing '/'"});
        }
        // too many items
        if (info && !info->variadic && isRecordKind(kind) && !info->params.empty()) {
            for (size_t r = 0; r < recSizes.size(); ++r) {
                if (recSizes[r].second == 0) continue;
                int maxIdx = 0;
                for (const Param* p : info->recordParams((int)r)) maxIdx = std::max(maxIdx, p->index);
                if (maxIdx > 0 && recSizes[r].second > maxIdx) {
                    const int l = recSizes[r].first;
                    out.push_back({l, 0, (int)lines[l].size(), 1,
                                   u + " record " + std::to_string(r + 1) + " has " + std::to_string(recSizes[r].second) + " items, at most " +
                                       std::to_string(maxIdx) + " are defined"});
                }
            }
        }
        if (kind == deck::KwSize::Fixed && spec && spec->count > 0 && records > spec->count && !(records == spec->count + 1 && emptyRecord))
            out.push_back({it.line, 0, nameEnd, 1, u + " takes " + std::to_string(spec->count) + " record(s), found " + std::to_string(records)});
        (void)toks;
    }
    std::sort(out.begin(), out.end(), [](const Diag& a, const Diag& b) { return a.line < b.line; });
    return out;
}

// ================================================================== autocomplete
std::vector<std::string> suggest(const std::vector<std::string>& lines, const Outline& o, int line, int col, const std::string& prefix, bool* keywordMode) {
    std::vector<std::string> out;
    if (keywordMode) *keywordMode = false;
    if (line < 0 || line >= (int)lines.size()) return out;
    const Docs& docs = Docs::get();
    const bool numeric = !prefix.empty() && (std::isdigit((unsigned char)prefix[0]) || prefix[0] == '.');
    if (col == 0 && !numeric) {
        if (keywordMode) *keywordMode = true;
        const std::string section = o.sectionAt(line);
        struct Cand { int rank; std::string name; };
        std::vector<Cand> c;
        const std::string up = upper(prefix);
        for (const auto& s : sectionNames()) if (startsWithNoCase(s, up)) c.push_back({0, s});
        for (const auto& kv : docs.all()) {
            if (!startsWithNoCase(kv.first, up) || isSectionName(kv.first)) continue;
            const bool valid = kv.second.validIn(section);
            if (!valid && up.empty()) continue;
            const deck::Support sup = deck::ressimSupport(kv.first);
            int rank = valid ? 1 : 3;
            if (sup == deck::Support::Full || sup == deck::Support::Partial) rank -= 0; else rank += 1;
            if (kv.first == up) rank = -1;
            c.push_back({rank, kv.first});
        }
        std::sort(c.begin(), c.end(), [](const Cand& a, const Cand& b) {
            if (a.rank != b.rank) return a.rank < b.rank;
            if (a.name.size() != b.name.size()) return a.name.size() < b.name.size();
            return a.name < b.name;
        });
        for (size_t i = 0; i < c.size() && i < 200; ++i) out.push_back(c[i].name);
        return out;
    }
    const Context ctx = contextAt(lines, o, line, col);
    if (!ctx.param || ctx.inComment || numeric) return out;
    const Param& p = *ctx.param;
    std::vector<std::string> pool;
    const std::string pn = upper(p.name);
    if (!p.options.empty()) pool = p.options;
    else if (pn.find("WELL") != std::string::npos && (p.isString() || p.valueType.empty())) pool = wellNames(lines, o);
    else if (pn.find("GROUP") != std::string::npos && (p.isString() || p.valueType.empty())) pool = groupNames(lines, o);
    else if (pn == "PHASE" || pn == "PREFERRED_PHASE" || pn == "TYPE") pool = {"OIL", "WATER", "GAS", "LIQ"};
    else if (pn == "STATUS" || pn == "STATE" || pn == "OPEN_SHUT_FLAG") pool = {"OPEN", "SHUT", "STOP", "AUTO"};
    std::string up = upper(prefix);
    while (!up.empty() && (up[0] == '\'' || up[0] == '"')) up.erase(0, 1);
    for (const auto& s : pool) {
        if (s.empty()) continue;
        if (startsWithNoCase(s, up)) out.push_back(s);
    }
    return out;
}

// ================================================================== text helpers
std::vector<std::string> splitLines(const std::string& text) {
    std::vector<std::string> v;
    size_t a = 0;
    while (a <= text.size()) {
        size_t b = text.find('\n', a);
        if (b == std::string::npos) b = text.size();
        std::string l = text.substr(a, b - a);
        if (!l.empty() && l.back() == '\r') l.pop_back();
        v.push_back(std::move(l));
        if (b == text.size()) break;
        a = b + 1;
    }
    return v;
}
std::string joinLines(const std::vector<std::string>& lines, size_t first, size_t last) {
    std::string s;
    last = std::min(last, lines.size() ? lines.size() - 1 : 0);
    for (size_t i = first; i < lines.size() && i <= last; ++i) { s += lines[i]; if (i < last) s += '\n'; }
    return s;
}
std::string externalizeLargeKeywords(const std::string& text, const std::string& relDir, const std::string& base, size_t threshold,
                                     std::vector<std::pair<std::string, std::string>>& files) {
    const auto lines = splitLines(text);
    const Outline o = buildOutline(lines);
    std::vector<std::string> out;
    out.reserve(lines.size());
    std::map<std::string, int> counter;
    size_t next = 0;   // next line to copy
    for (size_t i = 0; i < o.items.size(); ++i) {
        const int first = o.items[i].line, last = o.bodyEnd(int(i), int(lines.size()));
        size_t bytes = 0;
        for (int l = first; l <= last; ++l) bytes += lines[size_t(l)].size() + 1;
        if (bytes < threshold || o.items[i].section) continue;
        // keep trailing blank / comment lines in the main deck
        int end = last;
        while (end > first && lines[size_t(end)].find_first_not_of(" \t\r") == std::string::npos) --end;
        for (size_t l = next; l < size_t(first); ++l) out.push_back(lines[l]);
        const std::string kw = upper(o.items[i].name);
        const int n = ++counter[kw];
        std::string fname = base + "_" + kw + (n > 1 ? "_" + std::to_string(n) : std::string()) + ".INC";
        const std::string rel = relDir.empty() ? fname : relDir + "/" + fname;
        files.push_back({rel, joinLines(lines, size_t(first), size_t(end)) + "\n"});
        out.push_back("INCLUDE");
        out.push_back("  '" + rel + "' /");
        next = size_t(end) + 1;
    }
    for (size_t l = next; l < lines.size(); ++l) out.push_back(lines[l]);
    return joinLines(out);
}

// ================================================================== templates and alignment
namespace {

std::string headerName(const Param& p) {
    std::string n = p.name;
    if (n == "DATA" || n.find(' ') != std::string::npos) n = p.manual.empty() ? n : p.manual;
    for (char& ch : n) if (ch == ' ') ch = '_';
    return n;
}
std::string placeholder(const Param& p) {
    if (cleanDefault(p.def) && p.def != "1*") {
        const bool num = isNumber(p.def);
        if (p.isString() && !num) {
            bool word = true;
            for (char c : unquote(p.def)) if (!(std::isalnum((unsigned char)c) || c == '_')) word = false;
            if (word) return "'" + unquote(p.def) + "'";
        } else if (num) {
            return p.def;
        }
    }
    if (!p.def.empty()) return "1*";
    if (p.isString()) return "'" + headerName(p) + "'";
    return "1*";
}
// rows of tokens -> aligned lines; header (optional) occupies the same columns
std::vector<std::string> alignRows(const std::vector<std::vector<std::string>>& rows, const std::vector<std::string>& tails,
                                   const std::vector<std::string>* header) {
    std::vector<size_t> w;
    auto grow = [&](const std::vector<std::string>& r) {
        if (w.size() < r.size()) w.resize(r.size(), 0);
        for (size_t i = 0; i < r.size(); ++i) w[i] = std::max(w[i], r[i].size());
    };
    for (const auto& r : rows) grow(r);
    std::vector<std::string> out;
    if (header) {
        std::vector<std::string> h(*header);
        size_t ncols = 0;
        for (const auto& r : rows) ncols = std::max(ncols, r.size());
        if (ncols == 0) ncols = h.size();
        if (h.size() > ncols) h.resize(ncols);
        grow(h);
        std::string s = "--";
        for (size_t i = 0; i < h.size(); ++i) { s += ' '; s += h[i]; if (i + 1 < h.size()) s += std::string(w[i] - h[i].size(), ' '); }
        out.push_back(s);
    }
    for (size_t k = 0; k < rows.size(); ++k) {
        const auto& r = rows[k];
        std::string s = "  ";
        for (size_t i = 0; i < r.size(); ++i) { s += ' '; s += r[i]; if (i + 1 < r.size() || !tails[k].empty()) s += std::string(w[i] - r[i].size(), ' '); }
        if (!tails[k].empty()) { s += ' '; s += tails[k]; }
        while (!s.empty() && s.back() == ' ') s.pop_back();
        out.push_back(s);
    }
    return out;
}
std::vector<std::string> headerFor(const KwInfo* info, int recordIndex0) {
    std::vector<std::string> h;
    if (!info) return h;
    for (const Param* p : info->recordParams(recordIndex0)) h.push_back(headerName(*p));
    return h;
}
bool isHeaderComment(const std::string& s) {
    if (s.compare(0, 2, "--") != 0) return false;
    // generated headers are "-- NAME NAME ..." made of upper-case identifiers
    std::istringstream in(s.substr(2));
    std::string w; int n = 0;
    while (in >> w) {
        ++n;
        for (char c : w) if (!(std::isupper((unsigned char)c) || std::isdigit((unsigned char)c) || c == '_' || c == '/' || c == '-')) return false;
    }
    return n >= 2;
}

} // namespace

std::string recordTemplate(const std::string& kwIn, int cells) {
    const std::string kw = upper(kwIn);
    const KwInfo* info = Docs::get().find(kw);
    const deck::KwSize kind = kindOf(kw);
    std::string out;
    auto exampleBody = [&]() -> std::string {
        if (!info || info->example.empty()) return std::string();
        std::istringstream in(info->example);
        std::string l, body; bool seen = false;
        while (std::getline(in, l)) {
            if (!seen) { std::string n; if (keywordLine(l, &n) && upper(n) == kw) seen = true; continue; }
            std::string n;
            if (keywordLine(l, &n) && n != kw && Docs::get().find(upper(n)) && upper(n) != "ENDACTIO" && kind != deck::KwSize::Action) break;
            body += l + "\n";
        }
        while (body.size() > 1 && body[body.size() - 1] == '\n' && body[body.size() - 2] == '\n') body.pop_back();
        // examples are truncated in keyword_docs.txt: only use complete ones (last data line closes with '/')
        {
            std::istringstream bl(body); std::string x, lastData;
            while (std::getline(bl, x)) { std::vector<Token> tk; tokenize(x, tk); for (const Token& t : tk) if (t.kind != TK::Comment) lastData = x; }
            std::vector<Token> tk; tokenize(lastData, tk);
            bool closes = false;
            for (const Token& t : tk) if (t.kind == TK::Slash) closes = true;
            if (!closes && kind != deck::KwSize::Action && kind != deck::KwSize::Code && kind != deck::KwSize::Raw) return std::string();
        }
        return body;
    };
    switch (kind) {
    case deck::KwSize::None: return std::string();
    case deck::KwSize::Data: {
        std::string def = info && !info->params.empty() && cleanDefault(info->params[0].def) && isNumber(info->params[0].def) ? info->params[0].def : "0";
        if (cells > 0) out = "-- " + std::to_string(cells) + " values, one per cell (I fastest, then J, then K); n*value repeats\n   " + std::to_string(cells) + "*" + def + " /\n";
        else out = "-- one value per cell (I fastest, then J, then K); n*value repeats\n   " + def + " /\n";
        return out;
    }
    case deck::KwSize::Code: case deck::KwSize::Action: case deck::KwSize::Raw: case deck::KwSize::Tables: {
        std::string b = exampleBody();
        if (!b.empty()) return b;
        if (kind == deck::KwSize::Action) return "   'ACTION1' 1 0 /\n   FPR < 3000 /\n/\n-- keywords to run\nENDACTIO\n";
        if (kind == deck::KwSize::Raw) return "   text\n";
        return "/\n";
    }
    default: break;
    }
    if (!info || info->params.empty()) {
        const std::string b = exampleBody();
        return b.empty() ? std::string("   /\n") + (isListKind(kind) ? "/\n" : "") : b;
    }
    if (info->variadic && info->maxRecord() == 0) {
        const std::string b = exampleBody();
        if (!b.empty()) return b;
        std::vector<std::vector<std::string>> rows(1);
        for (size_t i = 0; i < info->params.size(); ++i) rows[0].push_back("0");
        const auto h = headerFor(info, 0);
        const auto lines = alignRows(rows, {std::string()}, &h);
        for (const auto& l : lines) out += l + "\n";
        return out + "/\n";
    }
    const int types = std::max(1, info->maxRecord());
    for (int r = 0; r < types; ++r) {
        const auto ps = info->recordParams(r);
        if (ps.empty()) continue;
        int lastReq = -1;
        for (size_t i = 0; i < ps.size(); ++i) if (ps[i]->required()) lastReq = (int)i;
        const size_t n = ps.size() <= 12 ? ps.size() : std::max<size_t>(lastReq + 1, 12);
        std::vector<std::string> vals, names;
        for (size_t i = 0; i < n && i < ps.size(); ++i) { vals.push_back(placeholder(*ps[i])); names.push_back(headerName(*ps[i])); }
        const auto lines = alignRows({vals}, {"/"}, &names);
        for (const auto& l : lines) out += l + "\n";
    }
    if (isListKind(kind)) out += "/\n";
    return out;
}

std::vector<std::string> alignBody(const std::string& kwIn, const std::vector<std::string>& body, bool addHeaders) {
    const std::string kw = upper(kwIn);
    const KwInfo* info = Docs::get().find(kw);
    const deck::KwSize kind = kindOf(kw);
    if (kind == deck::KwSize::Data || kind == deck::KwSize::Code || kind == deck::KwSize::Action || kind == deck::KwSize::Raw) return body;
    // group consecutive data lines of the same record type; comment lines break groups (header comments are regenerated)
    std::vector<std::string> out;
    std::vector<std::vector<std::string>> rows;
    std::vector<std::string> tails;
    int groupType = -2, groupRecord = 0;
    auto flush = [&]() {
        if (rows.empty()) return;
        std::vector<std::string> h;
        if (addHeaders && info && !(kind == deck::KwSize::Tables)) h = headerFor(info, groupRecord);
        const auto al = alignRows(rows, tails, h.empty() ? nullptr : &h);
        out.insert(out.end(), al.begin(), al.end());
        rows.clear(); tails.clear();
    };
    std::vector<Token> toks;
    int record = 0, item = 0;
    for (const std::string& line : body) {
        tokenize(line, toks);
        bool data = false;
        for (const Token& t : toks) if (t.kind != TK::Comment) data = true;
        if (!data) {
            const bool blank = line.find_first_not_of(" \t\r") == std::string::npos;
            if (addHeaders && isHeaderComment(line)) continue;
            flush(); groupType = -2;
            out.push_back(blank ? std::string() : line);
            continue;
        }
        // record type at line start
        const int type = info ? info->recordType(record) : 0;
        const bool lineStartsRecord = item == 0;
        if (type != groupType || !lineStartsRecord) { flush(); groupType = lineStartsRecord ? type : -3; groupRecord = record; }
        std::vector<std::string> row;
        std::string tail;
        bool emptyRecordLine = true;
        for (const Token& t : toks) {
            if (t.kind == TK::Comment) { if (!t.text.empty() && t.text.find_first_not_of(" \t\r") != std::string::npos) tail += (tail.empty() ? "" : " ") + std::string(t.text.substr(t.text.find_first_not_of(" \t\r"))); continue; }
            if (t.kind == TK::Slash) { tail = "/" + (tail.empty() ? std::string() : " " + tail); ++record; item = 0; continue; }
            emptyRecordLine = false;
            if (addHeaders && t.kind == TK::Defaults && t.count > 1 && !(info && info->variadic)) for (int k = 0; k < t.count; ++k) row.push_back("1*");
            else row.push_back(t.text);
            item += t.count;
        }
        // keep the slash right after the last token for a terminating comment "/ -- text"
        if (emptyRecordLine) { flush(); groupType = -2; out.push_back(tail.empty() ? "/" : tail); continue; }
        // trailing text after slash: tokenize put it in comment; rebuild as "/ text"
        rows.push_back(row);
        tails.push_back(tail);
        if (!lineStartsRecord) { flush(); groupType = -2; }
    }
    flush();
    return out;
}

// ================================================================== forms
bool Form::editable() const {
    if (!info) return false;
    const deck::KwSize k = kindOf(keyword);
    if (!isRecordKind(k)) return false;
    return !info->params.empty();
}

Form parseForm(const std::string& kwIn, const std::vector<std::string>& body) {
    Form f;
    f.keyword = upper(kwIn);
    f.info = Docs::get().find(f.keyword);
    if (f.info && f.info->variadic && !f.info->params.empty()) { f.table = true; f.columns = (int)f.info->params.size(); }
    std::vector<std::string> cur;
    std::vector<Token> toks;
    bool open = false;
    for (const std::string& line : body) {
        tokenize(line, toks);
        for (const Token& t : toks) {
            if (t.kind == TK::Comment) continue;
            if (t.kind == TK::Slash) {
                if (cur.empty() && !open) { f.terminated = true; continue; }
                f.records.push_back(cur); cur.clear(); open = false; continue;
            }
            if (f.terminated) continue;
            open = true;
            if (t.kind == TK::Defaults) for (int k = 0; k < t.count; ++k) cur.push_back("1*");
            else if (t.kind == TK::Repeat) {
                const std::string v = t.text.substr(t.text.find('*') + 1);
                for (int k = 0; k < t.count; ++k) cur.push_back(v);
            } else cur.push_back(t.text);
        }
    }
    if (!cur.empty()) f.records.push_back(cur);
    // drop trailing defaults so that rewriting does not grow records
    for (auto& r : f.records) while (!r.empty() && r.back() == "1*" && !f.table) r.pop_back();
    return f;
}

std::vector<std::string> writeForm(const Form& f) {
    std::vector<std::string> out;
    const deck::KwSize kind = kindOf(f.keyword);
    if (f.table) {
        const int nc = std::max(1, f.columns);
        for (size_t r = 0; r < f.records.size(); ++r) {
            std::vector<std::vector<std::string>> rows;
            std::vector<std::string> tails;
            const auto& v = f.records[r];
            for (size_t i = 0; i < v.size(); i += nc) {
                rows.emplace_back(v.begin() + i, v.begin() + std::min(v.size(), i + nc));
                tails.push_back(std::string());
            }
            if (!tails.empty()) tails.back() = "/";
            const auto h = headerFor(f.info, (int)r);
            auto al = alignRows(rows, tails, r == 0 ? &h : nullptr);
            if (rows.empty()) al.push_back("/");
            out.insert(out.end(), al.begin(), al.end());
        }
        return out;
    }
    // group consecutive records of the same type
    size_t r = 0;
    while (r < f.records.size()) {
        const int type = f.info ? f.info->recordType((int)r) : 0;
        size_t e = r;
        std::vector<std::vector<std::string>> rows;
        std::vector<std::string> tails;
        while (e < f.records.size() && (f.info ? f.info->recordType((int)e) : 0) == type) {
            std::vector<std::string> row = f.records[e];
            for (auto& s : row) if (s.empty()) s = "1*";
            rows.push_back(row); tails.push_back("/");
            ++e;
        }
        const auto h = headerFor(f.info, (int)r);
        const auto al = alignRows(rows, tails, h.empty() ? nullptr : &h);
        out.insert(out.end(), al.begin(), al.end());
        r = e;
    }
    if (isListKind(kind) || f.terminated) out.push_back("/");
    return out;
}

} // namespace deckintel
} // namespace ressim
