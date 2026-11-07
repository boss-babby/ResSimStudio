// ResSim - C31 user-defined quantities (UDQ / UDA) and ACTIONX: expression parser and evaluator (see udq.h)
#include "ressim/udq.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <sstream>

namespace ressim {
namespace udq {

namespace {
const double NaN = std::numeric_limits<double>::quiet_NaN();
std::string up(std::string s) { for (char& c : s) c = char(std::toupper(static_cast<unsigned char>(c))); return s; }
std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}
// strips a '--' comment outside quotes
std::string stripComment(const std::string& l) {
    bool q = false;
    for (size_t i = 0; i + 1 < l.size(); ++i) {
        if (l[i] == '\'') q = !q;
        if (!q && l[i] == '-' && l[i + 1] == '-') return l.substr(0, i);
    }
    return l;
}
// text of a record line up to its LAST '/' outside quotes (UDQ convention: '/' is also the division operator);
// returns false when the line has no terminating '/'
bool recordText(const std::string& line, std::string& rec) {
    bool q = false; long last = -1;
    for (size_t i = 0; i < line.size(); ++i) { if (line[i] == '\'') q = !q; else if (!q && line[i] == '/') last = long(i); }
    if (last < 0) { rec = line; return false; }
    rec = line.substr(0, size_t(last));
    return true;
}
bool isTimeWord(const std::string& u) { return u == "TIME" || u == "DAYS" || u == "DAY" || u == "MNTH" || u == "MONTH" || u == "YEAR" || u == "TIMESTEP"; }
int monthOf(const std::string& u) {
    static const char* m[] = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};
    for (int i = 0; i < 12; ++i) if (u == m[i]) return i + 1;
    if (u == "JLY") return 7;   // C47: ECLIPSE / OPM spelling of July (as in DATES; opm-tests ACTIONX_WCONINJE ACT-04)
    return 0;
}

// ---------------------------------------------------------------- tokens
enum class T { Num, Name, Str, Op, LP, RP, LB, RB, End };
// C38 union operators (OPM UADD ... UMAX): the result is defined where either operand is; one undefined operand -> the other
bool isUnionOp(const std::string& u) {
    static const char* ops[] = {"UADD", "USUB", "UMUL", "UDIV", "UPOW", "UMIN", "UMAX", "UEQ", "UNE", "UGT", "UGE", "ULT", "ULE"};
    for (const char* o : ops) if (u == o) return true;
    return false;
}
struct Tok { T t; std::string s; double v = 0; };
std::vector<Tok> lex(const std::string& s, std::string& err) {
    std::vector<Tok> out;
    size_t i = 0;
    while (i < s.size()) {
        const char c = s[i];
        if (std::isspace(static_cast<unsigned char>(c))) { ++i; continue; }
        if (c == '\'' || c == '"') {
            const size_t j = s.find(c, i + 1);
            if (j == std::string::npos) { err = "unterminated quote"; return out; }
            out.push_back({T::Str, s.substr(i + 1, j - i - 1)}); i = j + 1; continue;
        }
        if (std::isdigit(static_cast<unsigned char>(c)) || (c == '.' && i + 1 < s.size() && std::isdigit(static_cast<unsigned char>(s[i + 1])))) {
            size_t j = i;
            while (j < s.size() && (std::isdigit(static_cast<unsigned char>(s[j])) || s[j] == '.')) ++j;
            if (j < s.size() && (s[j] == 'e' || s[j] == 'E' || s[j] == 'd' || s[j] == 'D')) {
                size_t k = j + 1;
                if (k < s.size() && (s[k] == '+' || s[k] == '-')) ++k;
                if (k < s.size() && std::isdigit(static_cast<unsigned char>(s[k]))) { j = k; while (j < s.size() && std::isdigit(static_cast<unsigned char>(s[j]))) ++j; }
            }
            std::string num = s.substr(i, j - i);
            for (char& ch : num) if (ch == 'd' || ch == 'D') ch = 'e';
            out.push_back({T::Num, num, std::strtod(num.c_str(), nullptr)}); i = j; continue;
        }
        if (std::isalpha(static_cast<unsigned char>(c)) || c == '_' || c == '?') {
            size_t j = i;
            while (j < s.size()) {
                const char d = s[j];
                if (std::isalnum(static_cast<unsigned char>(d)) || d == '_' || d == '?' || d == ':' || d == '.') { ++j; continue; }
                // a '*' belongs to a well pattern when nothing operand-like follows it (OP*, 'P*' is quoted anyway)
                if (d == '*') {   // a well pattern star (OP*) unless an operand follows directly (WOPR*2, FOPR*(..))
                    const char e = j + 1 < s.size() ? s[j + 1] : ' ';
                    if (!(std::isalnum(static_cast<unsigned char>(e)) || e == '(' || e == '.' || e == '_')) { ++j; continue; }
                }
                break;
            }
            out.push_back({T::Name, s.substr(i, j - i)}); i = j; continue;
        }
        if (c == '(') { out.push_back({T::LP, "("}); ++i; continue; }
        if (c == '[') { out.push_back({T::LB, "["}); ++i; continue; }
        if (c == ']') { out.push_back({T::RB, "]"}); ++i; continue; }
        if (c == ')') { out.push_back({T::RP, ")"}); ++i; continue; }
        if (i + 1 < s.size()) {
            const std::string two = s.substr(i, 2);
            if (two == "==" || two == "!=" || two == "<=" || two == ">=" || two == "**") { out.push_back({T::Op, two == "**" ? "^" : two}); i += 2; continue; }
        }
        if (std::string("+-*/^<>=").find(c) != std::string::npos) { out.push_back({T::Op, std::string(1, c == '=' ? '=' : c)}); if (c == '=') out.back().s = "=="; ++i; continue; }
        err = std::string("unexpected character '") + c + "'";
        return out;
    }
    out.push_back({T::End, ""});
    return out;
}
}  // namespace

// ---------------------------------------------------------------- expression tree
struct Node {
    enum K { Num, Name, Bin, Un, Func, Tab } k = Num;
    double v = 0;
    std::string name, arg, op;     // Name: quantity / UDQ / time word, arg: well / group pattern; Bin / Un: op; Func: name
    bool hasArg = false;
    std::shared_ptr<Node> a, b;
};

namespace {
struct Parser {
    std::vector<Tok> t; size_t p = 0; std::string err;
    const Tok& cur() const { return t[std::min(p, t.size() - 1)]; }
    std::shared_ptr<Node> fail(const std::string& m) { if (err.empty()) err = m; return nullptr; }
    std::shared_ptr<Node> expr() { return uni(); }
    std::shared_ptr<Node> uni() {   // C38 union operators, lowest precedence, left associative
        auto l = cmp(); if (!l) return l;
        while (cur().t == T::Name && isUnionOp(up(cur().s))) {
            const std::string op = up(cur().s); ++p;
            auto r = cmp(); if (!r) return r;
            auto n = std::make_shared<Node>(); n->k = Node::Bin; n->op = op; n->a = l; n->b = r; l = n;
        }
        return l;
    }
    std::shared_ptr<Node> cmp() {
        auto l = add(); if (!l) return l;
        if (cur().t == T::Op && (cur().s == "==" || cur().s == "!=" || cur().s == "<" || cur().s == ">" || cur().s == "<=" || cur().s == ">=")) {
            const std::string op = cur().s; ++p;
            auto r = add(); if (!r) return r;
            auto n = std::make_shared<Node>(); n->k = Node::Bin; n->op = op; n->a = l; n->b = r; return n;
        }
        return l;
    }
    std::shared_ptr<Node> add() {
        auto l = mul(); if (!l) return l;
        while (cur().t == T::Op && (cur().s == "+" || cur().s == "-")) {
            const std::string op = cur().s; ++p;
            auto r = mul(); if (!r) return r;
            auto n = std::make_shared<Node>(); n->k = Node::Bin; n->op = op; n->a = l; n->b = r; l = n;
        }
        return l;
    }
    std::shared_ptr<Node> mul() {
        auto l = pw(); if (!l) return l;
        while (cur().t == T::Op && (cur().s == "*" || cur().s == "/")) {
            const std::string op = cur().s; ++p;
            auto r = pw(); if (!r) return r;
            auto n = std::make_shared<Node>(); n->k = Node::Bin; n->op = op; n->a = l; n->b = r; l = n;
        }
        return l;
    }
    std::shared_ptr<Node> pw() {
        auto l = unary(); if (!l) return l;
        if (cur().t == T::Op && cur().s == "^") {
            ++p; auto r = pw(); if (!r) return r;   // right associative
            auto n = std::make_shared<Node>(); n->k = Node::Bin; n->op = "^"; n->a = l; n->b = r; return n;
        }
        return l;
    }
    std::shared_ptr<Node> unary() {
        if (cur().t == T::Op && (cur().s == "-" || cur().s == "+")) {
            const std::string op = cur().s; ++p;
            auto x = unary(); if (!x) return x;
            if (op == "+") return x;
            auto n = std::make_shared<Node>(); n->k = Node::Un; n->op = "-"; n->a = x; return n;
        }
        return primary();
    }
    std::shared_ptr<Node> primary() {
        const Tok& k = cur();
        if (k.t == T::Num) { ++p; auto n = std::make_shared<Node>(); n->k = Node::Num; n->v = k.v; return n; }
        if (k.t == T::LP) {
            ++p; auto x = expr(); if (!x) return x;
            if (cur().t != T::RP) return fail("missing ')'");
            ++p; return x;
        }
        if (k.t == T::Name) {
            const std::string nm = up(k.s); ++p;
            if (cur().t == T::LB) {   // C38 UDT table lookup TU_NAME[expr]
                ++p; auto x = expr(); if (!x) return x;
                if (cur().t != T::RB) return fail("missing ']' after " + nm + "[");
                ++p;
                auto n = std::make_shared<Node>(); n->k = Node::Tab; n->name = nm; n->a = x; return n;
            }
            if (cur().t == T::LP) {   // function
                ++p; auto x = expr(); if (!x) return x;
                if (cur().t != T::RP) return fail("missing ')' after " + nm + "(");
                ++p;
                static const char* fn[] = {"ABS", "EXP", "LN", "LOG", "NINT", "SQRT", "DEF", "UNDEF", "IDV", "SORTA", "SORTD",
                                           "SUM", "AVEA", "AVEG", "AVEH", "MAX", "MIN", "PROD", "NORM1", "NORM2", "NORMI"};
                bool known = false; for (const char* f : fn) if (nm == f) known = true;
                if (!known) return fail("unknown UDQ function " + nm);
                auto n = std::make_shared<Node>(); n->k = Node::Func; n->name = nm; n->a = x; return n;
            }
            auto n = std::make_shared<Node>(); n->k = Node::Name; n->name = nm;
            // a well / group quantity may be followed by its well / group (quoted pattern or bare name)
            if ((nm[0] == 'W' || nm[0] == 'G') && (cur().t == T::Str || cur().t == T::Name) && !(cur().t == T::Name && t[std::min(p + 1, t.size() - 1)].t == T::LP) &&
                !(cur().t == T::Name && isUnionOp(up(cur().s)))) {
                n->arg = cur().s; n->hasArg = true; ++p;
            }
            return n;
        }
        if (k.t == T::Str) return fail("unexpected quoted text '" + k.s + "'");
        return fail(k.t == T::End ? "incomplete expression" : "unexpected '" + k.s + "'");
    }
};

bool hasWild(const std::string& s) { return s.find('*') != std::string::npos || s.find('?') != std::string::npos; }

Value binop(const std::string& op0, const Value& a, const Value& b) {
    const bool uni = op0.size() > 2 && op0[0] == 'U' && isUnionOp(op0);
    std::string op = op0;
    if (uni) {
        static const std::pair<const char*, const char*> m[] = {{"UADD", "+"}, {"USUB", "-"}, {"UMUL", "*"}, {"UDIV", "/"}, {"UPOW", "^"}, {"UMIN", "min"},
                                                                {"UMAX", "max"}, {"UEQ", "=="}, {"UNE", "!="}, {"UGT", ">"}, {"UGE", ">="}, {"ULT", "<"}, {"ULE", "<="}};
        for (const auto& [u, o] : m) if (op0 == u) op = o;
    }
    auto f = [&](double x, double y) -> double {
        if (uni && std::isnan(x)) return y;
        if (uni && std::isnan(y)) return x;
        if (std::isnan(x) || std::isnan(y)) return NaN;
        if (op == "min") return std::min(x, y);
        if (op == "max") return std::max(x, y);
        if (op == "+") return x + y;
        if (op == "-") return x - y;
        if (op == "*") return x * y;
        if (op == "/") return y == 0 ? NaN : x / y;
        if (op == "^") { const double r = std::pow(x, y); return std::isfinite(r) ? r : NaN; }
        if (op == "==") return x == y ? 1 : 0;
        if (op == "!=") return x != y ? 1 : 0;
        if (op == "<") return x < y ? 1 : 0;
        if (op == ">") return x > y ? 1 : 0;
        if (op == "<=") return x <= y ? 1 : 0;
        if (op == ">=") return x >= y ? 1 : 0;
        return NaN;
    };
    Value r;
    if (a.kind == Value::Scalar && b.kind == Value::Scalar) { r.s = f(a.s, b.s); return r; }
    if (a.kind != Value::Scalar && b.kind != Value::Scalar && (a.kind != b.kind || a.v.size() != b.v.size())) { r.s = NaN; return r; }
    r.kind = a.kind != Value::Scalar ? a.kind : b.kind;
    const size_t n = a.kind != Value::Scalar ? a.v.size() : b.v.size();
    r.v.resize(n);
    for (size_t i = 0; i < n; ++i) r.v[i] = f(a.kind == Value::Scalar ? a.s : a.v[i], b.kind == Value::Scalar ? b.s : b.v[i]);
    return r;
}

Value func(const std::string& fn, const Value& x) {
    static const char* setF[] = {"SUM", "AVEA", "AVEG", "AVEH", "MAX", "MIN", "PROD", "NORM1", "NORM2", "NORMI"};
    bool isSet = false; for (const char* f : setF) if (fn == f) isSet = true;
    if (isSet) {
        std::vector<double> d;
        if (x.kind == Value::Scalar) { if (!std::isnan(x.s)) d.push_back(x.s); }
        else for (double v : x.v) if (!std::isnan(v)) d.push_back(v);
        Value r;
        if (d.empty()) return r;   // undefined
        double acc = 0;
        if (fn == "SUM") { for (double v : d) acc += v; }
        else if (fn == "AVEA") { for (double v : d) acc += v; acc /= double(d.size()); }
        else if (fn == "AVEG") { for (double v : d) { if (v <= 0) return r; acc += std::log(v); } acc = std::exp(acc / double(d.size())); }
        else if (fn == "AVEH") { for (double v : d) { if (v == 0) return r; acc += 1 / v; } acc = double(d.size()) / acc; }
        else if (fn == "MAX") { acc = *std::max_element(d.begin(), d.end()); }
        else if (fn == "MIN") { acc = *std::min_element(d.begin(), d.end()); }
        else if (fn == "PROD") { acc = 1; for (double v : d) acc *= v; }
        else if (fn == "NORM1") { for (double v : d) acc += std::fabs(v); }
        else if (fn == "NORM2") { for (double v : d) acc += v * v; acc = std::sqrt(acc); }
        else if (fn == "NORMI") { for (double v : d) acc = std::max(acc, std::fabs(v)); }
        r.s = acc;
        return r;
    }
    auto f = [&](double v) -> double {
        if (fn == "DEF") return std::isnan(v) ? NaN : 1;
        if (fn == "UNDEF") return std::isnan(v) ? 1 : NaN;
        if (fn == "IDV") return std::isnan(v) ? 0 : 1;
        if (std::isnan(v)) return NaN;
        if (fn == "ABS") return std::fabs(v);
        if (fn == "EXP") { const double e = std::exp(v); return std::isfinite(e) ? e : NaN; }
        if (fn == "LN") return v > 0 ? std::log(v) : NaN;
        if (fn == "LOG") return v > 0 ? std::log10(v) : NaN;
        if (fn == "NINT") return std::round(v);          // half away from zero (Fortran NINT)
        if (fn == "SQRT") return v >= 0 ? std::sqrt(v) : NaN;
        return NaN;
    };
    Value r = x;
    if (fn == "SORTA" || fn == "SORTD") {   // rank (1 = smallest / largest), undefined stay undefined
        if (x.kind == Value::Scalar) { r.s = std::isnan(x.s) ? NaN : 1; return r; }
        std::vector<size_t> id;
        for (size_t i = 0; i < x.v.size(); ++i) if (!std::isnan(x.v[i])) id.push_back(i);
        std::stable_sort(id.begin(), id.end(), [&](size_t i, size_t j) { return fn == "SORTA" ? x.v[i] < x.v[j] : x.v[i] > x.v[j]; });
        for (double& v : r.v) v = NaN;
        for (size_t k = 0; k < id.size(); ++k) r.v[id[k]] = double(k + 1);
        return r;
    }
    if (r.kind == Value::Scalar) r.s = f(r.s); else for (double& v : r.v) v = f(v);
    return r;
}

Value evalNode(const Node& n, const Source& src, const Lookup& udqs) {
    Value r;
    switch (n.k) {
    case Node::Num: r.s = n.v; return r;
    case Node::Un: { Value x = evalNode(*n.a, src, udqs); if (x.kind == Value::Scalar) x.s = -x.s; else for (double& v : x.v) v = -v; return x; }
    case Node::Bin: return binop(n.op, evalNode(*n.a, src, udqs), evalNode(*n.b, src, udqs));
    case Node::Func: return func(n.name, evalNode(*n.a, src, udqs));
    case Node::Tab: {   // C38 UDT
        Value x = evalNode(*n.a, src, udqs);
        const Table* tb = nullptr;
        if (src.tables) for (const auto& t : *src.tables) if (t.name == n.name) tb = &t;   // the latest definition
        if (!tb) { Value u; if (x.kind != Value::Scalar) { u.kind = x.kind; u.v.assign(x.v.size(), NaN); } return u; }
        if (x.kind == Value::Scalar) x.s = tb->eval(x.s); else for (double& v : x.v) v = tb->eval(v);
        return x;
    }
    case Node::Name: {
        const std::string& nm = n.name;
        const Value* u = udqs ? udqs(nm) : nullptr;
        const bool isWell = nm[0] == 'W', isGroup = nm[0] == 'G';
        const std::vector<std::string>& names = isWell ? src.wells : src.groups;
        if (u) {
            if (!n.hasArg || u->kind == Value::Scalar) return *u;
            Value x = *u;
            size_t hit = 0, last = 0;
            for (size_t i = 0; i < x.v.size() && i < names.size(); ++i) { if (wildMatch(n.arg, names[i])) { ++hit; last = i; } else x.v[i] = NaN; }
            if (!hasWild(n.arg) && hit == 1) { Value s; s.s = x.v[last]; return s; }
            return x;
        }
        if (isTimeWord(nm)) {
            r.s = nm == "TIME" || nm == "DAYS" ? src.days : nm == "DAY" ? src.day : (nm == "MNTH" || nm == "MONTH") ? src.month : nm == "YEAR" ? src.year : NaN;
            return r;
        }
        if (isWell || isGroup) {
            r.kind = isWell ? Value::Well : Value::Group;
            r.v.assign(names.size(), NaN);
            size_t hit = 0, last = 0;
            for (size_t i = 0; i < names.size(); ++i) {
                if (n.hasArg && !wildMatch(n.arg, names[i])) continue;
                ++hit; last = i;
                r.v[i] = isWell ? (src.well ? src.well(nm, i) : NaN) : (src.group ? src.group(nm, i) : NaN);
            }
            if (n.hasArg && !hasWild(n.arg) && hit == 1) { Value s; s.s = r.v[last]; return s; }   // one named well / group: a scalar
            if (n.hasArg && hit == 0) { Value s; return s; }
            return r;
        }
        r.s = src.field ? src.field(nm) : NaN;
        return r;
    }
    }
    return r;
}
}  // namespace

bool wildMatch(const std::string& pat, const std::string& name) {
    const std::string P = up(pat), N = up(name);
    // iterative glob with backtracking
    size_t p = 0, n = 0, star = std::string::npos, mark = 0;
    while (n < N.size()) {
        if (p < P.size() && (P[p] == '?' || P[p] == N[n])) { ++p; ++n; }
        else if (p < P.size() && P[p] == '*') { star = p++; mark = n; }
        else if (star != std::string::npos) { p = star + 1; n = ++mark; }
        else return false;
    }
    while (p < P.size() && P[p] == '*') ++p;
    return p == P.size();
}

bool parseExpr(const std::string& text, Expr& e, std::string& err) {
    Parser P;
    P.t = lex(text, P.err);
    if (!P.err.empty()) { err = P.err; return false; }
    e.text = text;
    e.root = P.expr();
    if (!e.root) { err = P.err.empty() ? "syntax error" : P.err; return false; }
    if (P.cur().t != T::End) { err = "unexpected '" + P.cur().s + "' after the expression"; e.root.reset(); return false; }
    return true;
}

Value evaluate(const Expr& e, const Source& src, const Lookup& udqs) {
    if (!e.root) return Value{};
    return evalNode(*e.root, src, udqs);
}

// ---------------------------------------------------------------- UDQ set
double Table::eval(double v) const {
    if (std::isnan(v) || x.empty() || y.size() != x.size()) return NaN;
    const size_t n = x.size();
    if (n == 1) return y[0];
    if (interp == "NV") {   // nearest value
        size_t b = 0; for (size_t i = 1; i < n; ++i) if (std::fabs(x[i] - v) < std::fabs(x[b] - v)) b = i;
        return y[b];
    }
    size_t k = 0;
    if (v <= x[0]) { if (interp != "LL") return y[0]; k = 0; }
    else if (v >= x[n - 1]) { if (interp != "LL") return y[n - 1]; k = n - 2; }
    else { while (k + 2 < n && v > x[k + 1]) ++k; }
    const double dx = x[k + 1] - x[k];
    return dx != 0 ? y[k] + (y[k + 1] - y[k]) * (v - x[k]) / dx : y[k];
}

bool parseTables(const std::string& body, int step, std::vector<Table>& out, std::string& err) {
    // UDT records: 'NAME' dims / ; then per table 'INTERP' x1 x2 .. / y1 y2 .. / ; an empty record ends the table; an
    // empty record after that ends the keyword. Several tables may follow each other.
    std::vector<std::vector<std::string>> recs;
    {
        std::istringstream in(body); std::string line, acc;
        while (std::getline(in, line)) {
            line = stripComment(line);
            std::string rec;
            const bool term = recordText(line, rec);
            acc += " " + rec;
            if (!term) continue;
            std::istringstream rs(acc); std::vector<std::string> tk; std::string w;
            while (rs >> w) {
                if (w.size() >= 2 && (w.front() == '\'' || w.front() == '"') && w.back() == w.front()) w = w.substr(1, w.size() - 2);
                tk.push_back(w);
            }
            recs.push_back(tk); acc.clear();
        }
    }
    size_t r = 0; bool ok = true;
    auto num = [](const std::string& s, double& v) { char* e = nullptr; v = std::strtod(s.c_str(), &e); return e && *e == 0; };
    while (r < recs.size()) {
        if (recs[r].empty()) { ++r; continue; }
        Table t; t.step = step; t.name = up(recs[r][0]);
        const int dims = recs[r].size() > 1 ? std::atoi(recs[r][1].c_str()) : 1;
        ++r;
        if (dims != 1) { err += (err.empty() ? "" : "; ") + t.name + ": only 1-D UDT tables are supported"; ok = false; }
        if (r >= recs.size() || recs[r].empty()) { err += (err.empty() ? "" : "; ") + t.name + ": no interpolation / x record"; return false; }
        t.interp = up(recs[r][0]);
        if (t.interp != "NV" && t.interp != "LC" && t.interp != "LL") { err += (err.empty() ? "" : "; ") + t.name + ": interpolation " + t.interp + " (NV, LC or LL)"; ok = false; t.interp = "LC"; }
        for (size_t i = 1; i < recs[r].size(); ++i) { double v; if (num(recs[r][i], v)) t.x.push_back(v); }
        ++r;
        if (r < recs.size()) { for (const auto& s : recs[r]) { double v; if (num(s, v)) t.y.push_back(v); } ++r; }
        if (t.x.size() != t.y.size() || t.x.empty()) { err += (err.empty() ? "" : "; ") + t.name + ": x and y counts differ"; ok = false; }
        for (size_t i = 1; i < t.x.size(); ++i) if (t.x[i] <= t.x[i - 1]) { err += (err.empty() ? "" : "; ") + t.name + ": x not increasing"; ok = false; break; }
        if (!t.name.empty() && t.name[0] != 'T') { err += (err.empty() ? "" : "; ") + t.name + ": UDT names start with T"; }
        out.push_back(t);
        if (r < recs.size() && recs[r].empty()) ++r;   // end of this table
        if (r < recs.size() && recs[r].empty()) break; // end of the keyword
    }
    return ok;
}

bool Set::addTables(const std::string& body, int step, std::string& err) { return parseTables(body, step, tables_, err); }

bool Set::add(const std::string& body, int step, std::string& err) {
    std::istringstream in(body);
    std::string line, acc;
    bool ok = true;
    while (std::getline(in, line)) {
        line = stripComment(line);
        std::string rec;
        const bool term = recordText(line, rec);
        acc += " " + rec;
        if (!term) continue;
        std::string r = trim(acc);
        acc.clear();
        if (r.empty()) { if (trim(line) == "/") break; continue; }   // a lone '/' ends the keyword
        std::istringstream rs(r);
        std::string op, name;
        rs >> op >> name;
        op = up(op); name = up(name);
        if (!name.empty() && (name.front() == '\'' || name.front() == '"')) name = name.substr(1, name.size() - 2);
        std::string rest; std::getline(rs, rest); rest = trim(rest);
        auto bad = [&](const std::string& m) { ok = false; err += (err.empty() ? "" : "; ") + name + ": " + m; };
        if (name.size() < 3 || (name[0] != 'F' && name[0] != 'W' && name[0] != 'G') || name[1] != 'U') { bad("UDQ names start with FU, WU or GU"); continue; }
        if (op.rfind("DEFI", 0) == 0) {
            Definition d; d.name = name; d.var = name[0]; d.kind = 0; d.step = step;
            std::string e2;
            if (!parseExpr(rest, d.expr, e2)) { bad(e2); continue; }
            defs_.push_back(d);
        } else if (op.rfind("ASSI", 0) == 0) {
            Definition d; d.name = name; d.var = name[0]; d.kind = 1; d.step = step;
            std::vector<std::string> tk; { std::istringstream ts(rest); std::string x; while (ts >> x) tk.push_back(x); }
            if (tk.empty()) { bad("ASSIGN without a value"); continue; }
            char* endp = nullptr;
            std::string vs = tk.back(); for (char& ch : vs) if (ch == 'd' || ch == 'D') ch = 'e';
            d.value = std::strtod(vs.c_str(), &endp);
            if (!endp || *endp) { bad("ASSIGN value '" + tk.back() + "' is not a number"); continue; }
            if (tk.size() > 1) { d.selector = tk[0]; if (!d.selector.empty() && d.selector.front() == '\'') d.selector = d.selector.substr(1, d.selector.size() - 2); }
            defs_.push_back(d);
        } else if (op.rfind("UPDA", 0) == 0) {
            const std::string m = up(trim(rest));
            updRec_.push_back({name, m.rfind("OFF", 0) == 0 ? 1 : m.rfind("NEXT", 0) == 0 ? 2 : 0, step});
        } else if (op.rfind("UNIT", 0) == 0) {
            std::string u = trim(rest);
            if (u.size() >= 2 && (u.front() == '\'' || u.front() == '"')) u = u.substr(1, u.size() - 2);
            units_.push_back({name, u});
        } else bad("unknown UDQ operation " + op);
    }
    return ok;
}

Value& Set::slot(const std::string& name) {
    for (auto& v : vals_) if (v.first == name) return v.second;
    vals_.push_back({name, Value{}});
    return vals_.back().second;
}

void Set::evaluate(const Source& src0, int step, bool assignsOnly) {
    Lookup look = [this](const std::string& n) { return get(n); };
    // C38 UDT tables active at this step (latest definition last)
    std::vector<Table> act;
    for (const auto& t : tables_) if (t.step <= step) act.push_back(t);
    Source src = src0;
    src.tables = &act;
    // C38 order (OPM UDQConfig::eval): the due ASSIGNs in deck order, then the DEFINEs - well, group, field quantities, each
    // in the order of the name's first appearance in the UDQ input. So a UDQ that refers to a later-named one sees its
    // previous value (or undefined), e.g. ASSIGN X 0 / DEFINE X0 X / DEFINE X (..) UMIN X0: X is evaluated first.
    std::vector<size_t> order;
    {
        std::vector<std::string> first;
        for (const auto& d : defs_) if (std::find(first.begin(), first.end(), d.name) == first.end()) first.push_back(d.name);
        auto rank = [&](size_t i) {
            const Definition& d = defs_[i];
            const int vr = d.kind == 1 ? -1 : d.var == 'W' ? 0 : d.var == 'G' ? 1 : 2;
            const long fi = long(std::find(first.begin(), first.end(), d.name) - first.begin());
            return std::make_pair(vr, d.kind == 1 ? long(i) : fi);
        };
        for (size_t i = 0; i < defs_.size(); ++i) order.push_back(i);
        std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return rank(a) < rank(b); });
    }
    for (size_t oi = 0; oi < order.size(); ++oi) {
        const size_t i = order[oi];
        Definition& d = defs_[i];
        if (d.step > step) continue;
        const size_t nn = d.var == 'W' ? src.wells.size() : d.var == 'G' ? src.groups.size() : 0;
        const std::vector<std::string>& names = d.var == 'W' ? src.wells : src.groups;
        if (d.kind == 1) {   // ASSIGN at the report step it appears (once)
            if (d.assigned) continue;
            d.assigned = true;
            Value& v = slot(d.name);
            if (d.var == 'F') { v.kind = Value::Scalar; v.s = d.value; continue; }
            v.kind = d.var == 'W' ? Value::Well : Value::Group;
            if (v.v.size() != nn) v.v.assign(nn, NaN);
            for (size_t k = 0; k < nn; ++k) if (d.selector.empty() || wildMatch(d.selector, names[k])) v.v[k] = d.value;
            continue;
        }
        if (assignsOnly) continue;
        // the latest active DEFINE of this name only
        bool later = false;
        for (size_t j = i + 1; j < defs_.size(); ++j) if (defs_[j].name == d.name && defs_[j].kind == 0 && defs_[j].step <= step) later = true;
        if (later) continue;
        int mode = 0, recStep = -1;
        for (const Upd& u : updRec_) if (u.name == d.name && u.step <= step && u.step >= recStep) { mode = u.mode; recStep = u.step; }
        if (mode == 1) continue;
        if (mode == 2) {
            bool done = false;
            for (auto& nd : nextDone_) if (nd.first == d.name && nd.second == recStep) done = true;
            if (done) continue;
            nextDone_.push_back({d.name, recStep});
        }
        Value r = udq::evaluate(d.expr, src, look);
        Value& v = slot(d.name);
        if (d.var == 'F') {
            v.kind = Value::Scalar;
            v.s = r.kind == Value::Scalar ? r.s : NaN;     // a well / group set in a field UDQ needs a set function
        } else {
            v.kind = d.var == 'W' ? Value::Well : Value::Group;
            if (r.kind == Value::Scalar) v.v.assign(nn, r.s);
            else if ((r.kind == Value::Well) == (d.var == 'W') && r.v.size() == nn) v.v = r.v;
            else v.v.assign(nn, NaN);
        }
    }
}

const Value* Set::get(const std::string& name) const {
    for (const auto& v : vals_) if (v.first == name) return &v.second;
    return nullptr;
}
bool Set::has(const std::string& name) const {
    for (const auto& d : defs_) if (d.name == name) return true;
    return false;
}
double Set::scalar(const std::string& name) const {
    const Value* v = get(name);
    return v && v->kind == Value::Scalar ? v->s : NaN;
}
double Set::wellValue(const std::string& name, size_t w) const {
    const Value* v = get(name);
    if (!v) return NaN;
    if (v->kind == Value::Scalar) return v->s;
    return v->kind == Value::Well && w < v->v.size() ? v->v[w] : NaN;
}
double Set::groupValue(const std::string& name, size_t g) const {
    const Value* v = get(name);
    if (!v) return NaN;
    if (v->kind == Value::Scalar) return v->s;
    return v->kind == Value::Group && g < v->v.size() ? v->v[g] : NaN;
}
std::vector<std::string> Set::names() const {
    std::vector<std::string> n;
    for (const auto& d : defs_) if (std::find(n.begin(), n.end(), d.name) == n.end()) n.push_back(d.name);
    return n;
}
std::string Set::units(const std::string& name) const {
    std::string u;
    for (const auto& x : units_) if (x.first == name) u = x.second;
    return u;
}

std::string Set::save() const {
    std::ostringstream o;
    o.precision(17);
    o << defs_.size();
    for (const auto& d : defs_) o << ' ' << (d.assigned ? 1 : 0);
    o << ' ' << nextDone_.size();
    for (const auto& n : nextDone_) o << ' ' << n.first << ' ' << n.second;
    o << ' ' << vals_.size();
    for (const auto& v : vals_) {
        o << ' ' << v.first << ' ' << int(v.second.kind) << ' ' << v.second.s << ' ' << v.second.v.size();
        for (double x : v.second.v) o << ' ' << x;
    }
    return o.str();
}
void Set::load(const std::string& blob) {
    std::istringstream in(blob);
    auto rd = [&in]() { std::string s; in >> s; return s; };
    auto rdD = [&]() { const std::string s = rd(); return s == "nan" || s == "-nan" ? NaN : std::strtod(s.c_str(), nullptr); };
    size_t n = 0; in >> n;
    for (size_t i = 0; i < n; ++i) { int a = 0; in >> a; if (i < defs_.size()) defs_[i].assigned = a != 0; }
    size_t nd = 0; in >> nd; nextDone_.clear();
    for (size_t i = 0; i < nd; ++i) { std::string nm; int st = 0; in >> nm >> st; nextDone_.push_back({nm, st}); }
    size_t nv = 0; in >> nv; vals_.clear();
    for (size_t i = 0; i < nv; ++i) {
        Value v; std::string nm = rd(); int k = 0; in >> k; v.kind = Value::Kind(k); v.s = rdD();
        size_t m = 0; in >> m; v.v.resize(m);
        for (size_t j = 0; j < m; ++j) v.v[j] = rdD();
        vals_.push_back({nm, v});
    }
}

// ---------------------------------------------------------------- ACTIONX
struct CondNode {
    enum K { And, Or, Cmp } k = Cmp;
    std::shared_ptr<CondNode> a, b;
    std::shared_ptr<Node> l, r;
    std::string op;
};

namespace {
struct CondParser {
    std::vector<Tok> t; size_t p = 0; std::string err;
    const Tok& cur() const { return t[std::min(p, t.size() - 1)]; }
    bool word(const char* w) const { return cur().t == T::Name && up(cur().s) == w; }
    std::shared_ptr<CondNode> orE() {
        auto l = andE(); if (!l) return l;
        while (word("OR")) { ++p; auto r = andE(); if (!r) return r; auto n = std::make_shared<CondNode>(); n->k = CondNode::Or; n->a = l; n->b = r; l = n; }
        return l;
    }
    std::shared_ptr<CondNode> andE() {
        auto l = term(); if (!l) return l;
        while (word("AND")) { ++p; auto r = term(); if (!r) return r; auto n = std::make_shared<CondNode>(); n->k = CondNode::And; n->a = l; n->b = r; l = n; }
        return l;
    }
    std::shared_ptr<Node> operand() {
        const Tok& k = cur();
        if (k.t == T::Num) { ++p; auto n = std::make_shared<Node>(); n->k = Node::Num; n->v = k.v; return n; }
        if (k.t == T::Op && k.s == "-" && t[std::min(p + 1, t.size() - 1)].t == T::Num) { ++p; auto n = std::make_shared<Node>(); n->k = Node::Num; n->v = -cur().v; ++p; return n; }
        if (k.t == T::Name) {
            const std::string nm = up(k.s); ++p;
            if (const int m = monthOf(nm)) { auto n = std::make_shared<Node>(); n->k = Node::Num; n->v = m; return n; }
            auto n = std::make_shared<Node>(); n->k = Node::Name; n->name = nm;
            if ((nm[0] == 'W' || nm[0] == 'G') && (cur().t == T::Str || (cur().t == T::Name && up(cur().s) != "AND" && up(cur().s) != "OR"))) { n->arg = cur().s; n->hasArg = true; ++p; }
            return n;
        }
        if (err.empty()) err = k.t == T::End ? "incomplete condition" : "unexpected '" + k.s + "' in a condition";
        return nullptr;
    }
    std::shared_ptr<CondNode> term() {
        if (cur().t == T::LP) {
            ++p; auto x = orE(); if (!x) return x;
            if (cur().t != T::RP) { if (err.empty()) { err = "missing ')' in a condition"; } return nullptr; }
            ++p; return x;
        }
        auto l = operand(); if (!l) return nullptr;
        if (cur().t != T::Op || !(cur().s == "==" || cur().s == "!=" || cur().s == "<" || cur().s == ">" || cur().s == "<=" || cur().s == ">=")) {
            if (err.empty()) err = "comparison operator expected";
            return nullptr;
        }
        const std::string op = cur().s; ++p;
        auto r = operand(); if (!r) return nullptr;
        auto n = std::make_shared<CondNode>(); n->k = CondNode::Cmp; n->l = l; n->r = r; n->op = op;
        return n;
    }
};

struct CondVal { bool ok = false; bool wellSet = false; std::vector<char> w; };   // w: per well satisfied
CondVal evalCond(const CondNode& c, const Source& src, const Lookup& look) {
    CondVal r;
    if (c.k == CondNode::Cmp) {
        const Value a = evalNode(*c.l, src, look), b = evalNode(*c.r, src, look);
        const Value x = binop(c.op, a, b);
        if (x.kind == Value::Scalar) { r.ok = x.s == 1; return r; }
        if (x.kind == Value::Well) {
            r.wellSet = true; r.w.assign(x.v.size(), 0);
            for (size_t i = 0; i < x.v.size(); ++i) if (x.v[i] == 1) { r.w[i] = 1; r.ok = true; }
            return r;
        }
        for (double v : x.v) if (v == 1) r.ok = true;    // group set: any group
        return r;
    }
    const CondVal a = evalCond(*c.a, src, look), b = evalCond(*c.b, src, look);
    if (c.k == CondNode::And) {
        r.ok = a.ok && b.ok;
        if (a.wellSet && b.wellSet) { r.wellSet = true; r.w.resize(a.w.size()); for (size_t i = 0; i < r.w.size(); ++i) r.w[i] = char(a.w[i] && i < b.w.size() && b.w[i]); }
        else if (a.wellSet) { r.wellSet = true; r.w = a.w; } else if (b.wellSet) { r.wellSet = true; r.w = b.w; }
        if (!r.ok) for (char& v : r.w) v = 0;
    } else {
        r.ok = a.ok || b.ok;
        if (a.wellSet || b.wellSet) {
            r.wellSet = true;
            const size_t n = std::max(a.w.size(), b.w.size());
            r.w.assign(n, 0);
            for (size_t i = 0; i < n; ++i) r.w[i] = char((i < a.w.size() && a.w[i]) || (i < b.w.size() && b.w[i]));
        }
    }
    return r;
}
}  // namespace

bool parseAction(const std::string& raw, int step, Action& a, std::string& err) {
    std::istringstream in(raw);
    std::string line;
    int stage = 0;     // 0 header, 1 conditions, 2 body
    std::string cond, body;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::string clean = trim(stripComment(line));
        if (stage < 2 && clean.empty()) continue;
        if (stage == 0) {
            std::string rec; recordText(clean, rec);
            std::istringstream rs(rec);
            std::string nm; rs >> nm;
            if (!nm.empty() && nm.front() == '\'') nm = nm.substr(1, nm.size() - 2);
            if (nm.empty()) { err = "ACTIONX without a name"; return false; }
            a.name = nm;
            std::string x;
            if (rs >> x && x != "1*") a.maxRuns = std::max(1, std::atoi(x.c_str()));
            if (rs >> x && x != "1*") a.minWaitDays = std::max(0.0, std::strtod(x.c_str(), nullptr));
            stage = 1; continue;
        }
        if (stage == 1) {
            if (clean == "/") { stage = 2; continue; }
            std::string rec; recordText(clean, rec);
            cond += " " + rec;
            a.conditionText += (a.conditionText.empty() ? "" : "\n") + trim(rec);
            continue;
        }
        if (up(clean).rfind("ENDACTIO", 0) == 0) break;
        body += line + "\n";
    }
    a.fromStep = step;
    a.body = body;
    CondParser P;
    P.t = lex(cond, P.err);
    if (!P.err.empty()) { err = "ACTIONX " + a.name + ": " + P.err; return false; }
    a.cond = P.orE();
    if (!a.cond || P.cur().t != T::End) { err = "ACTIONX " + a.name + ": " + (P.err.empty() ? "unexpected '" + P.cur().s + "'" : P.err); a.cond.reset(); return false; }
    return true;
}

Trigger evaluateAction(const Action& a, const Source& src0, const Set& udqs) {
    Trigger t;
    if (!a.cond) return t;
    Source src = src0;
    if (!src.tables) src.tables = &udqs.tables();   // C38 UDT in conditions
    Lookup look = [&udqs](const std::string& n) { return udqs.get(n); };
    const CondVal v = evalCond(*a.cond, src, look);
    t.fired = v.ok;
    if (v.wellSet) for (size_t i = 0; i < v.w.size(); ++i) if (v.w[i]) t.wells.push_back(i);
    return t;
}

std::string expandBody(const std::string& body, const std::vector<std::string>& wells) {
    std::istringstream in(body);
    std::ostringstream o;
    std::string line;
    while (std::getline(in, line)) {
        const std::string t = trim(line);
        const bool q1 = t.rfind("'?'", 0) == 0, q2 = t.size() >= 2 && t[0] == '?' && std::isspace(static_cast<unsigned char>(t[1]));
        if (!q1 && !q2) { o << line << "\n"; continue; }
        const size_t at = line.find(q1 ? "'?'" : "?");
        for (const auto& w : wells) o << line.substr(0, at) << "'" << w << "'" << line.substr(at + (q1 ? 3 : 1)) << "\n";
    }
    return o.str();
}

}  // namespace udq
}  // namespace ressim
