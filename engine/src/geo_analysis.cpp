// C40 geomodel data analysis: log calculator, regression, variogram fitting, well tie (see geo_analysis.h)
#include "ressim/geo_analysis.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <memory>

namespace ressim {
namespace geo {

namespace {
const double kNaN = std::numeric_limits<double>::quiet_NaN();
std::string upper(std::string s) { for (char& c : s) c = char(std::toupper((unsigned char)c)); return s; }

// ------------------------------------------------------------------ expression parser -> tree, evaluated vector-wise
struct Node {
    enum Kind { Num, Var, Neg, Not, Bin, Fn } kind = Num;
    double num = 0;
    const std::vector<double>* var = nullptr;
    std::string op;                    // binary operator or function name
    std::vector<std::unique_ptr<Node>> a;
};
using NodeP = std::unique_ptr<Node>;

struct Parser {
    const std::string& s;
    const std::map<std::string, const std::vector<double>*>& vars;
    size_t p = 0;
    std::string err;
    Parser(const std::string& str, const std::map<std::string, const std::vector<double>*>& v) : s(str), vars(v) {}
    void ws() { while (p < s.size() && std::isspace((unsigned char)s[p])) ++p; }
    bool eat(const char* t) {
        ws();
        const size_t n = std::char_traits<char>::length(t);
        if (s.compare(p, n, t) == 0) { p += n; return true; }
        return false;
    }
    NodeP fail(const std::string& m) { if (err.empty()) err = m + " at position " + std::to_string(p + 1); return nullptr; }
    NodeP bin(const std::string& op, NodeP l, NodeP r) {
        if (!l || !r) return nullptr;
        NodeP n(new Node); n->kind = Node::Bin; n->op = op; n->a.push_back(std::move(l)); n->a.push_back(std::move(r)); return n;
    }
    NodeP parseOr() { NodeP l = parseAnd(); while (l && eat("||")) l = bin("||", std::move(l), parseAnd()); return l; }
    NodeP parseAnd() { NodeP l = parseCmp(); while (l && eat("&&")) l = bin("&&", std::move(l), parseCmp()); return l; }
    NodeP parseCmp() {
        NodeP l = parseAdd();
        while (l) {
            const char* ops[] = {"<=", ">=", "==", "!=", "<", ">"};
            const char* hit = nullptr;
            for (const char* o : ops) { const size_t save = p; if (eat(o)) { hit = o; break; } p = save; }
            if (!hit) break;
            l = bin(hit, std::move(l), parseAdd());
        }
        return l;
    }
    NodeP parseAdd() {
        NodeP l = parseMul();
        while (l) { if (eat("+")) l = bin("+", std::move(l), parseMul()); else if (eat("-")) l = bin("-", std::move(l), parseMul()); else break; }
        return l;
    }
    NodeP parseMul() {
        NodeP l = parseUnary();
        while (l) { if (eat("*")) l = bin("*", std::move(l), parseUnary()); else if (eat("/")) l = bin("/", std::move(l), parseUnary()); else break; }
        return l;
    }
    NodeP parseUnary() {
        if (eat("-")) { NodeP x = parseUnary(); if (!x) return nullptr; NodeP n(new Node); n->kind = Node::Neg; n->a.push_back(std::move(x)); return n; }
        if (eat("+")) return parseUnary();
        { const size_t save = p; if (eat("!") && !(p < s.size() && s[p] == '=')) { NodeP x = parseUnary(); if (!x) return nullptr; NodeP n(new Node); n->kind = Node::Not; n->a.push_back(std::move(x)); return n; } p = save; }
        return parsePow();
    }
    NodeP parsePow() {
        NodeP b = parseAtom();
        if (b && eat("^")) return bin("^", std::move(b), parseUnary());   // right associative, binds tighter than unary minus on the left
        return b;
    }
    NodeP parseAtom() {
        ws();
        if (p >= s.size()) return fail("unexpected end of the expression");
        if (eat("(")) { NodeP e = parseOr(); if (!e) return nullptr; if (!eat(")")) return fail("missing ')'"); return e; }
        const char c = s[p];
        if (std::isdigit((unsigned char)c) || c == '.') {
            char* end = nullptr;
            const double v = std::strtod(s.c_str() + p, &end);
            if (end == s.c_str() + p) return fail("bad number");
            p = size_t(end - s.c_str());
            NodeP n(new Node); n->kind = Node::Num; n->num = v; return n;
        }
        if (std::isalpha((unsigned char)c) || c == '_') {
            const size_t b = p;
            while (p < s.size() && (std::isalnum((unsigned char)s[p]) || s[p] == '_' || s[p] == '.')) ++p;
            const std::string name = upper(s.substr(b, p - b));
            ws();
            if (p < s.size() && s[p] == '(') {   // function call
                ++p;
                NodeP n(new Node); n->kind = Node::Fn; n->op = name;
                ws();
                if (!eat(")")) {
                    for (;;) {
                        NodeP x = parseOr(); if (!x) return nullptr;
                        n->a.push_back(std::move(x));
                        if (eat(",")) continue;
                        if (eat(")")) break;
                        return fail("expected ',' or ')'");
                    }
                }
                static const std::map<std::string, int> arity = {{"ABS", 1}, {"SQRT", 1}, {"LN", 1}, {"LOG", 1}, {"LOG10", 1}, {"EXP", 1}, {"ISNAN", 1},
                                                                 {"POW", 2}, {"MIN", 2}, {"MAX", 2}, {"CLAMP", 3}, {"IF", 3}};
                const auto it = arity.find(name);
                if (it == arity.end()) return fail("unknown function '" + name + "'");
                if (int(n->a.size()) != it->second) return fail(name + " needs " + std::to_string(it->second) + " argument(s)");
                return n;
            }
            if (name == "PI") { NodeP n(new Node); n->num = 3.14159265358979323846; return n; }
            if (name == "NAN" || name == "NULL") { NodeP n(new Node); n->num = kNaN; return n; }
            const auto it = vars.find(name);
            if (it == vars.end()) return fail("unknown curve '" + name + "'");
            NodeP n(new Node); n->kind = Node::Var; n->var = it->second; return n;
        }
        return fail(std::string("unexpected '") + c + "'");
    }
};

double evalAt(const Node& n, size_t i) {
    switch (n.kind) {
        case Node::Num: return n.num;
        case Node::Var: return i < n.var->size() ? (*n.var)[i] : kNaN;
        case Node::Neg: return -evalAt(*n.a[0], i);
        case Node::Not: { const double x = evalAt(*n.a[0], i); return std::isnan(x) ? kNaN : (x == 0 ? 1.0 : 0.0); }
        case Node::Bin: {
            const double x = evalAt(*n.a[0], i), y = evalAt(*n.a[1], i);
            const std::string& o = n.op;
            if (o == "+") return x + y;
            if (o == "-") return x - y;
            if (o == "*") return x * y;
            if (o == "/") return y == 0 ? kNaN : x / y;
            if (o == "^") return std::pow(x, y);
            if (std::isnan(x) || std::isnan(y)) {   // undefined (null) operands: null, except a definite short-circuit
                if (o == "||" && ((!std::isnan(x) && x != 0) || (!std::isnan(y) && y != 0))) return 1.0;
                if (o == "&&" && ((!std::isnan(x) && x == 0) || (!std::isnan(y) && y == 0))) return 0.0;
                return kNaN;
            }
            if (o == "<") return x < y;
            if (o == "<=") return x <= y;
            if (o == ">") return x > y;
            if (o == ">=") return x >= y;
            if (o == "==") return x == y;
            if (o == "!=") return x != y;
            if (o == "&&") return (x != 0 && y != 0) ? 1.0 : 0.0;
            if (o == "||") return (x != 0 || y != 0) ? 1.0 : 0.0;
            return kNaN;
        }
        case Node::Fn: {
            const std::string& f = n.op;
            if (f == "IF") { const double c = evalAt(*n.a[0], i); return std::isnan(c) ? kNaN : (c != 0 ? evalAt(*n.a[1], i) : evalAt(*n.a[2], i)); }
            const double x = evalAt(*n.a[0], i);
            if (f == "ABS") return std::fabs(x);
            if (f == "SQRT") return x < 0 ? kNaN : std::sqrt(x);
            if (f == "LN" || f == "LOG") return x <= 0 ? kNaN : std::log(x);
            if (f == "LOG10") return x <= 0 ? kNaN : std::log10(x);
            if (f == "EXP") return std::exp(x);
            if (f == "ISNAN") return std::isnan(x) ? 1.0 : 0.0;
            const double y = evalAt(*n.a[1], i);
            if (f == "POW") return std::pow(x, y);
            if (f == "MIN") return std::isnan(x) || std::isnan(y) ? kNaN : std::min(x, y);
            if (f == "MAX") return std::isnan(x) || std::isnan(y) ? kNaN : std::max(x, y);
            if (f == "CLAMP") { const double hi = evalAt(*n.a[2], i); return std::isnan(x) ? kNaN : std::min(std::max(x, y), hi); }
            return kNaN;
        }
    }
    return kNaN;
}
}  // namespace

bool evalExpression(const std::string& expr, const std::map<std::string, const std::vector<double>*>& vars, size_t n,
                    std::vector<double>& out, std::string& err) {
    std::map<std::string, const std::vector<double>*> uv;
    for (const auto& kv : vars) uv[upper(kv.first)] = kv.second;
    Parser P(expr, uv);
    NodeP root = P.parseOr();
    if (root) { P.ws(); if (P.p < expr.size()) { root.reset(); P.fail("unexpected text"); } }
    if (!root) { err = P.err.empty() ? std::string("syntax error") : P.err; return false; }
    out.assign(n, kNaN);
    for (size_t i = 0; i < n; ++i) out[i] = evalAt(*root, i);
    err.clear();
    return true;
}

std::vector<std::string> expressionNames(const WellLog& log) {
    std::vector<std::string> v{"MD", "TVD"};
    for (const LogCurve& c : log.curves) v.push_back(c.mnem);
    return v;
}

bool addCalculatedCurve(WellLog& log, const GeoWell* w, const std::string& mnem, const std::string& unit, const std::string& expr, std::string& err) {
    if (mnem.empty()) { err = "name of the new curve is empty"; return false; }
    const size_t n = log.md.size();
    std::vector<std::vector<double>> clean(log.curves.size());
    std::map<std::string, const std::vector<double>*> vars;
    for (size_t c = 0; c < log.curves.size(); ++c) {
        clean[c] = log.curves[c].v;
        for (double& x : clean[c]) if (std::fabs(x - log.nullValue) < 1e-6 * std::max(1.0, std::fabs(log.nullValue))) x = kNaN;
        vars[log.curves[c].mnem] = &clean[c];
    }
    std::vector<double> tvd(n, kNaN);
    if (w) for (size_t i = 0; i < n; ++i) tvd[i] = w->at(log.md[i])[2];
    vars["MD"] = &log.md;
    vars["TVD"] = &tvd;
    std::vector<double> out;
    if (!evalExpression(expr, vars, n, out, err)) return false;
    for (double& x : out) if (!std::isfinite(x)) x = log.nullValue;
    for (LogCurve& c : log.curves)
        if (upper(c.mnem) == upper(mnem)) { c.v = std::move(out); c.unit = unit; c.descr = "= " + expr; return true; }
    LogCurve c; c.mnem = upper(mnem); c.unit = unit; c.descr = "= " + expr; c.v = std::move(out);
    log.curves.push_back(std::move(c));
    return true;
}

std::vector<CalcPreset> calcPresets() {
    return {
        {"Shale volume, linear gamma ray", "VSH", "V/V", "clamp((GR - 30) / (120 - 30), 0, 1)", "IGR = (GR - GRclean) / (GRshale - GRclean); edit the 30 / 120 end points"},
        {"Shale volume, Larionov (tertiary rocks)", "VSH", "V/V", "clamp(0.083 * (pow(2, 3.7 * clamp((GR - 30) / 90, 0, 1)) - 1), 0, 1)", "Vsh = 0.083 (2^(3.7 IGR) - 1)"},
        {"Shale volume, Larionov (older rocks)", "VSH", "V/V", "clamp(0.33 * (pow(2, 2 * clamp((GR - 30) / 90, 0, 1)) - 1), 0, 1)", "Vsh = 0.33 (2^(2 IGR) - 1)"},
        {"Density porosity (sandstone matrix)", "PHID", "V/V", "clamp((2.65 - RHOB) / (2.65 - 1.0), 0, 0.5)", "PHID = (rho_ma - RHOB) / (rho_ma - rho_fl)"},
        {"Sonic porosity, Wyllie (us/ft)", "PHIS", "V/V", "clamp((DT - 55.5) / (189 - 55.5), 0, 0.5)", "PHIS = (DT - DTma) / (DTfl - DTma)"},
        {"Effective porosity from total porosity and Vsh", "PHIE", "V/V", "PHIT * (1 - VSH)", "needs PHIT and VSH"},
        {"Water saturation, Archie", "SW", "V/V", "clamp(pow(1 * 0.05 / (pow(PHIE, 2) * RT), 1 / 2), 0, 1)", "Sw = (a Rw / (phi^m Rt))^(1/n); a = 1, Rw = 0.05, m = n = 2"},
        {"Permeability, Timur", "PERM", "MD", "0.136 * pow(PHIE * 100, 4.4) / pow(max(SWIRR, 0.05) * 100, 2)", "k = 0.136 phi^4.4 / Swirr^2 (percent); needs SWIRR"},
        {"Net reservoir flag", "NET", "", "if(PHIE >= 0.08 && VSH <= 0.4, 1, 0)", "porosity and shale cut-offs"},
        {"Acoustic impedance (DT us/ft, RHOB g/cc)", "AI", "", "RHOB * 304800 / DT", "AI = rho * v (g/cc * m/s)"},
    };
}

// ------------------------------------------------------------------ regression
Regression regress(const std::vector<double>& x, const std::vector<double>& y, bool logY) {
    Regression r; r.logY = logY;
    double sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0; int n = 0;
    for (size_t i = 0; i < std::min(x.size(), y.size()); ++i) {
        double yv = y[i];
        if (!std::isfinite(x[i]) || !std::isfinite(yv)) continue;
        if (logY) { if (yv <= 0) continue; yv = std::log10(yv); }
        sx += x[i]; sy += yv; sxx += x[i] * x[i]; sxy += x[i] * yv; syy += yv * yv; ++n;
    }
    r.n = n;
    if (n < 2) return r;
    const double d = n * sxx - sx * sx;
    if (std::fabs(d) < 1e-300) return r;
    r.b = (n * sxy - sx * sy) / d;
    r.a = (sy - r.b * sx) / n;
    const double st = syy - sy * sy / n, se = st - r.b * (sxy - sx * sy / n);
    r.r2 = st > 0 ? std::max(0.0, 1.0 - se / st) : 0.0;
    return r;
}

// ------------------------------------------------------------------ variogram fitting
double structure(geostat::VarioModel m, double h) {
    if (h <= 0) return 0.0;
    switch (m) {
        case geostat::VarioModel::Spherical: return h >= 1 ? 1.0 : 1.5 * h - 0.5 * h * h * h;
        case geostat::VarioModel::Exponential: return 1.0 - std::exp(-3.0 * h);
        case geostat::VarioModel::Gaussian: return 1.0 - std::exp(-3.0 * h * h);
    }
    return 1.0;
}

namespace {
// nugget, sill >= 0 minimising sum w (g - n - s f)^2 for a fixed range; returns the weighted SSE
double solveNS(const geostat::ExpVario& ev, const std::vector<double>& w, geostat::VarioModel m, double range, double fixN, double& nug, double& sill) {
    double a11 = 0, a12 = 0, a22 = 0, b1 = 0, b2 = 0;
    for (size_t k = 0; k < ev.lag.size(); ++k) {
        if (w[k] <= 0) continue;
        const double f = structure(m, ev.lag[k] / range);
        a11 += w[k]; a12 += w[k] * f; a22 += w[k] * f * f; b1 += w[k] * ev.gamma[k]; b2 += w[k] * f * ev.gamma[k];
    }
    auto onlySill = [&](double n0) { return a22 > 0 ? std::max(0.0, (b2 - n0 * a12) / a22) : 0.0; };
    if (fixN >= 0) { nug = fixN; sill = onlySill(fixN); }
    else {
        const double d = a11 * a22 - a12 * a12;
        nug = d > 1e-300 ? (b1 * a22 - b2 * a12) / d : 0.0;
        sill = d > 1e-300 ? (a11 * b2 - a12 * b1) / d : 0.0;
        if (nug < 0) { nug = 0; sill = onlySill(0); }
        if (sill < 0) { sill = 0; nug = a11 > 0 ? std::max(0.0, b1 / a11) : 0.0; }
    }
    double sse = 0;
    for (size_t k = 0; k < ev.lag.size(); ++k) {
        if (w[k] <= 0) continue;
        const double r = ev.gamma[k] - nug - sill * structure(m, ev.lag[k] / range);
        sse += w[k] * r * r;
    }
    return sse;
}
}  // namespace

VarioFit fitVariogram(const geostat::ExpVario& ev, geostat::VarioModel m, double fixNugget, double maxRange) {
    VarioFit f;
    double hmax = 0; int used = 0;
    for (size_t k = 0; k < ev.lag.size(); ++k) if (k < ev.pairs.size() && ev.pairs[k] > 0 && std::isfinite(ev.gamma[k])) { hmax = std::max(hmax, ev.lag[k]); ++used; }
    if (used < 2 || hmax <= 0) return f;
    if (maxRange <= 0) maxRange = 2.0 * hmax;
    std::vector<double> w(ev.lag.size(), 0.0);
    for (size_t k = 0; k < ev.lag.size(); ++k)
        if (k < ev.pairs.size() && ev.pairs[k] > 0 && std::isfinite(ev.gamma[k]) && ev.lag[k] > 0) w[k] = double(ev.pairs[k]) / (ev.lag[k] * ev.lag[k]);
    // iterate: weights N / gamma_model^2 (Cressie 1985) starting from N / h^2
    double best = 0, bn = 0, bs = 0, br = hmax * 0.5;
    for (int it = 0; it < 3; ++it) {
        best = 1e300;
        const int ng = 80;
        for (int g = 0; g < ng; ++g) {
            const double r = hmax * 0.05 * std::pow(maxRange / (hmax * 0.05), double(g) / (ng - 1));
            double n0, s0;
            const double e = solveNS(ev, w, m, r, fixNugget, n0, s0);
            if (e < best) { best = e; bn = n0; bs = s0; br = r; }
        }
        // golden refinement around the best grid range
        double lo = br / 1.15, hi = br * 1.15;
        for (int g = 0; g < 30; ++g) {
            const double a = hi - 0.618 * (hi - lo), b = lo + 0.618 * (hi - lo);
            double n1, s1, n2, s2;
            if (solveNS(ev, w, m, a, fixNugget, n1, s1) < solveNS(ev, w, m, b, fixNugget, n2, s2)) hi = b; else lo = a;
        }
        br = 0.5 * (lo + hi);
        best = solveNS(ev, w, m, br, fixNugget, bn, bs);
        for (size_t k = 0; k < ev.lag.size(); ++k) {
            if (w[k] <= 0) continue;
            const double gm = std::max(1e-12 * (bn + bs + 1e-30), bn + bs * structure(m, ev.lag[k] / br));
            w[k] = double(ev.pairs[k]) / (gm * gm);
        }
    }
    double sse = 0; int n = 0;
    for (size_t k = 0; k < ev.lag.size(); ++k) {
        if (k >= ev.pairs.size() || ev.pairs[k] <= 0 || !std::isfinite(ev.gamma[k])) continue;
        const double r = ev.gamma[k] - bn - bs * structure(m, ev.lag[k] / br);
        sse += r * r; ++n;
    }
    f.nugget = bn; f.sill = bs; f.range = br; f.rmse = n ? std::sqrt(sse / n) : 0; f.ok = bs > 0;
    return f;
}

// ------------------------------------------------------------------ well tie
VelocityModel fitVelocityModel(const std::vector<double>& z, const std::vector<double>& t, double datum, double* rmsMs) {
    VelocityModel best; best.datum = datum;
    double bestE = 1e300;
    auto misfit = [&](const VelocityModel& v) {
        double e = 0; int n = 0;
        for (size_t i = 0; i < std::min(z.size(), t.size()); ++i) { const double r = depthToTime(v, z[i]) - t[i]; e += r * r; ++n; }
        return n ? e / n : 1e300;
    };
    // for each gradient k the best v0 by golden search (time decreases monotonically with v0)
    for (int kg = 0; kg <= 60; ++kg) {
        VelocityModel v; v.datum = datum; v.k = 2.0 * kg / 60.0;
        double lo = 100, hi = 20000;
        for (int g = 0; g < 60; ++g) {
            const double a = hi - 0.618 * (hi - lo), b = lo + 0.618 * (hi - lo);
            VelocityModel va = v, vb = v; va.v0 = a; vb.v0 = b;
            if (misfit(va) < misfit(vb)) hi = b; else lo = a;
        }
        v.v0 = 0.5 * (lo + hi);
        const double e = misfit(v);
        if (e < bestE) { bestE = e; best = v; }
    }
    if (rmsMs) *rmsMs = std::sqrt(bestE);
    return best;
}

TieResult wellTie(const GeoWell& w, const WellLog& log, const SeismicCube& cube, const VelocityModel& vm, const TieInput& in) {
    TieResult R;
    const LogCurve* dt = log.curve(in.dt);
    if (!dt) { R.msg = "no sonic curve '" + in.dt + "' in the log of " + w.name + " (add DT, or compute it with the log calculator)"; return R; }
    const LogCurve* rho = in.rhob.empty() ? nullptr : log.curve(in.rhob);
    if (cube.ns < 8 || cube.a.empty()) { R.msg = "no seismic volume loaded"; return R; }
    // samples with a valid slowness, in depth order
    std::vector<double> z, slow, den;
    double med = 0; { std::vector<double> v; for (double x : dt->v) if (std::isfinite(x) && std::fabs(x - log.nullValue) > 1e-6 && x > 0) v.push_back(x);
                      if (v.size() < 10) { R.msg = "fewer than 10 valid sonic samples"; return R; }
                      std::nth_element(v.begin(), v.begin() + long(v.size() / 2), v.end()); med = v[v.size() / 2]; }
    int unit = in.dtUnit;
    if (unit == 0) unit = med < 300 ? 1 : 2;   // us/ft typically 40-200, us/m 130-650
    const double lenPerM = in.field ? 3.280839895 : 1.0;   // project length units per metre
    for (size_t i = 0; i < log.md.size() && i < dt->v.size(); ++i) {
        const double s = dt->v[i];
        if (!std::isfinite(s) || std::fabs(s - log.nullValue) < 1e-6 || s <= 0) continue;
        const double sPerM = unit == 1 ? s * 1e-6 * 3.280839895 : s * 1e-6;   // seconds per metre
        double r = in.rhoConst;
        if (rho && i < rho->v.size() && std::isfinite(rho->v[i]) && std::fabs(rho->v[i] - log.nullValue) > 1e-6 && rho->v[i] > 0.5) r = rho->v[i];
        const double zz = w.at(log.md[i])[2];
        if (!z.empty() && zz <= z.back()) continue;
        z.push_back(zz); slow.push_back(sPerM / lenPerM); den.push_back(r);   // slowness per project length unit
    }
    if (z.size() < 10) { R.msg = "not enough sonic samples along increasing depth"; return R; }
    // time-depth: anchor at the first sample with the current velocity model, integrate the sonic below
    std::vector<double> t(z.size());
    t[0] = depthToTime(vm, z[0]);
    for (size_t i = 1; i < z.size(); ++i) t[i] = t[i - 1] + 2000.0 * (z[i] - z[i - 1]) * 0.5 * (slow[i] + slow[i - 1]);
    // impedance on a regular time axis
    const double ts = std::max(0.25, in.dtMs);
    const double t0 = t.front(), t1 = t.back();
    const int nt = int((t1 - t0) / ts) + 1;
    if (nt < 16) { R.msg = "logged interval shorter than 16 ms"; return R; }
    std::vector<double> tt(static_cast<size_t>(nt)), ai(static_cast<size_t>(nt));
    size_t q = 0;
    for (int k = 0; k < nt; ++k) {
        const double tk = t0 + k * ts;
        while (q + 1 < t.size() && t[q + 1] < tk) ++q;
        const size_t q1 = std::min(q + 1, t.size() - 1);
        const double f = t[q1] > t[q] ? std::clamp((tk - t[q]) / (t[q1] - t[q]), 0.0, 1.0) : 0.0;
        const double sl = slow[q] + f * (slow[q1] - slow[q]), de = den[q] + f * (den[q1] - den[q]);
        tt[size_t(k)] = tk;
        ai[size_t(k)] = de * (1.0 / sl) / lenPerM * 1.0;   // g/cc * m/s
    }
    std::vector<double> rc(size_t(nt), 0.0);
    for (int k = 0; k + 1 < nt; ++k) rc[size_t(k)] = (ai[size_t(k) + 1] - ai[size_t(k)]) / (ai[size_t(k) + 1] + ai[size_t(k)]);
    // Ricker wavelet, +-1.5 / f
    const double pi2 = 9.869604401089358;
    const int hw = int(std::ceil(1500.0 / std::max(1.0, in.freq) / ts));
    std::vector<double> wav(size_t(2 * hw + 1));
    for (int k = -hw; k <= hw; ++k) { const double a = pi2 * in.freq * in.freq * (k * ts / 1000.0) * (k * ts / 1000.0); wav[size_t(k + hw)] = (1 - 2 * a) * std::exp(-a); }
    std::vector<double> syn(size_t(nt), 0.0);
    for (int k = 0; k < nt; ++k) {
        double s = 0;
        for (int m = -hw; m <= hw; ++m) { const int j = k - m; if (j >= 0 && j < nt) s += rc[size_t(j)] * wav[size_t(m + hw)]; }
        syn[size_t(k)] = s;
    }
    // seismic trace at the well: plan position at the middle of the logged interval
    const double mdMid = log.md[log.md.size() / 2];
    const P3 pm = w.at(mdMid);
    double fi = 0, fj = 0;
    if (!cube.ij(pm[0], pm[1], fi, fj)) { R.msg = "the well is outside the seismic survey"; return R; }
    const int ci = std::clamp(int(std::lround(fi)), 0, cube.ni - 1), cj = std::clamp(int(std::lround(fj)), 0, cube.nj - 1);
    for (int s = 0; s < cube.ns; ++s) { R.traceT.push_back(cube.t0 + s * cube.dt); R.trace.push_back(cube.at(ci, cj, s)); }
    // correlation of the synthetic with the trace for a bulk shift (synthetic sampled at trace times)
    auto synAt = [&](double tq) {
        const double u = (tq - t0) / ts;
        if (u < 0 || u > nt - 1) return kNaN;
        const int k = std::min(int(u), nt - 2);
        const double f = u - k;
        return syn[size_t(k)] * (1 - f) + syn[size_t(k) + 1] * f;
    };
    auto corrFor = [&](double shift) {
        double sa = 0, sb = 0, saa = 0, sbb = 0, sab = 0; int n = 0;
        for (size_t s = 0; s < R.traceT.size(); ++s) {
            const double a = synAt(R.traceT[s] - shift);
            if (!std::isfinite(a)) continue;
            const double b = R.trace[s];
            sa += a; sb += b; saa += a * a; sbb += b * b; sab += a * b; ++n;
        }
        if (n < 8) return -2.0;
        const double va = saa - sa * sa / n, vb = sbb - sb * sb / n;
        return va > 0 && vb > 0 ? (sab - sa * sb / n) / std::sqrt(va * vb) : -2.0;
    };
    R.corr0 = corrFor(0);
    double bestC = -3, bestS = 0;
    for (double s = -in.maxShift; s <= in.maxShift + 1e-9; s += 0.5) { const double c = corrFor(s); if (c > bestC) { bestC = c; bestS = s; } }
    if (bestC < -1.5) { R.msg = "the logged interval does not overlap the seismic time window"; return R; }
    R.shift = bestS; R.corr = bestC;
    R.twt.resize(size_t(nt));
    for (int k = 0; k < nt; ++k) R.twt[size_t(k)] = tt[size_t(k)] + bestS;
    R.ai = std::move(ai); R.rc = std::move(rc); R.synth = std::move(syn);
    for (size_t i = 0; i < z.size(); i += std::max<size_t>(1, z.size() / 400)) { R.tdZ.push_back(z[i]); R.tdT.push_back(t[i] + bestS); }
    R.vmFit = fitVelocityModel(R.tdZ, R.tdT, vm.datum, &R.vmRmsMs);
    R.ok = true;
    char b[200];
    std::snprintf(b, sizeof b, "correlation %.2f -> %.2f with a bulk shift of %+.1f ms; calibrated v(z) = %.0f + %.3f (z - %.0f), misfit %.1f ms",
                  R.corr0, R.corr, R.shift, R.vmFit.v0, R.vmFit.k, R.vmFit.datum, R.vmRmsMs);
    R.msg = b;
    return R;
}

// ---------------------------------------------------------------- C42 calculators / statistics
std::string exprName(const std::string& name) {
    std::string r;
    for (char c : name) {
        const unsigned char u = (unsigned char)c;
        r += std::isalnum(u) ? char(std::toupper(u)) : '_';
    }
    if (r.empty() || std::isdigit((unsigned char)r[0])) r = "_" + r;
    return r;
}

bool surfaceCalculator(const std::vector<Surface>& hz, int onto, const std::string& expr, Surface& out, std::string& err) {
    if (onto < 0 || onto >= (int)hz.size() || !hz[size_t(onto)].valid()) { err = "the target lattice surface is not valid"; return false; }
    const Surface& L = hz[size_t(onto)];
    const size_t n = size_t(L.nx) * size_t(L.ny);
    std::vector<double> X(n), Y(n);
    for (int j = 0; j < L.ny; ++j)
        for (int i = 0; i < L.nx; ++i) L.nodeXY(i, j, X[size_t(i + L.nx * j)], Y[size_t(i + L.nx * j)]);
    std::vector<std::vector<double>> res;
    res.reserve(hz.size());
    std::map<std::string, const std::vector<double>*> vars;
    for (size_t h = 0; h < hz.size(); ++h) {
        if (!hz[h].valid()) { res.emplace_back(); continue; }
        if ((int)h == onto) res.push_back(L.z);
        else {
            std::vector<double> z(n);
            for (size_t k = 0; k < n; ++k) z[k] = hz[h].at(X[k], Y[k]);
            res.push_back(std::move(z));
        }
    }
    for (size_t h = 0; h < hz.size(); ++h)
        if (!res[h].empty()) vars[exprName(hz[h].name)] = &res[h];
    vars["X"] = &X;
    vars["Y"] = &Y;
    std::vector<double> v;
    if (!evalExpression(expr, vars, n, v, err)) return false;
    out = L;
    out.z = std::move(v);
    for (double& z : out.z)
        if (!std::isfinite(z)) z = std::numeric_limits<double>::quiet_NaN();
    return true;
}

Stats describe(const std::vector<double>& v) {
    Stats s;
    std::vector<double> f;
    f.reserve(v.size());
    for (double x : v) (std::isfinite(x) ? (void)f.push_back(x) : (void)++s.nNull);
    s.n = f.size();
    if (f.empty()) return s;
    std::sort(f.begin(), f.end());
    double sum = 0, sq = 0;
    for (double x : f) sum += x;
    s.sum = sum;
    s.mean = sum / double(f.size());
    for (double x : f) sq += (x - s.mean) * (x - s.mean);
    s.sd = f.size() > 1 ? std::sqrt(sq / double(f.size() - 1)) : 0.0;
    s.min = f.front();
    s.max = f.back();
    auto q = [&](double p) {
        const double r = p * double(f.size() - 1);
        const size_t a = size_t(r), b = std::min(a + 1, f.size() - 1);
        return f[a] + (r - double(a)) * (f[b] - f[a]);
    };
    s.p10 = q(0.10); s.p50 = q(0.50); s.p90 = q(0.90);
    return s;
}

}  // namespace geo
}  // namespace ressim
