// C32 multisegment wells: WELSEGS / COMPSEGS / WSEGVALV / WSEGSICD / WSEGAICD parsing and the segment pressure solve
// (see msw.h for the model).
#include "ressim/msw.h"
#include "ressim/case.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <map>
#include <sstream>

namespace ressim {
namespace msw {
namespace {

std::string up(std::string s) { for (char& c : s) c = char(std::toupper(static_cast<unsigned char>(c))); return s; }
std::string unq(const std::string& s) { return s.size() >= 2 && (s.front() == '\'' || s.front() == '"') ? s.substr(1, s.size() - 2) : s; }

// records of a keyword body: '/'-terminated, comments removed, N* defaults expanded to empty items, N*v repeated
std::vector<std::vector<std::string>> records(const std::string& body) {
    std::vector<std::vector<std::string>> out;
    std::vector<std::string> cur;
    std::istringstream is(body);
    std::string line;
    bool anyRec = false;
    while (std::getline(is, line)) {
        const size_t c = line.find("--");
        if (c != std::string::npos) line = line.substr(0, c);
        size_t i = 0;
        while (i < line.size()) {
            const char ch = line[i];
            if (std::isspace(static_cast<unsigned char>(ch))) { ++i; continue; }
            if (ch == '/') {
                if (cur.empty() && anyRec) { out.push_back({}); }   // empty record ends the keyword
                else out.push_back(cur);
                anyRec = true; cur.clear(); ++i;
                if (!out.empty() && out.back().empty()) return out;
                continue;
            }
            std::string tok;
            if (ch == '\'' || ch == '"') {
                const size_t e = line.find(ch, i + 1);
                tok = line.substr(i, (e == std::string::npos ? line.size() : e + 1) - i);
                i = e == std::string::npos ? line.size() : e + 1;
            } else {
                while (i < line.size() && !std::isspace(static_cast<unsigned char>(line[i])) && line[i] != '/') tok += line[i++];
            }
            const size_t st = tok.find('*');
            if (st != std::string::npos && st > 0 && tok[0] != '\'' && std::isdigit(static_cast<unsigned char>(tok[0]))) {
                const int n = std::atoi(tok.substr(0, st).c_str());
                const std::string v = tok.substr(st + 1);
                for (int r = 0; r < n; ++r) cur.push_back(v);
            } else cur.push_back(tok);
        }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}
bool has(const std::vector<std::string>& r, size_t i) { return i < r.size() && !r[i].empty() && r[i] != "1*"; }
double num(const std::vector<std::string>& r, size_t i, double def) { return has(r, i) ? std::atof(r[i].c_str()) : def; }
int inum(const std::vector<std::string>& r, size_t i, int def) { return has(r, i) ? std::atoi(r[i].c_str()) : def; }
bool wild(const std::string& pat, const std::string& name) {   // * and ? patterns
    const std::string p = up(pat), n = up(name);
    size_t pi = 0, ni = 0, star = std::string::npos, mark = 0;
    while (ni < n.size()) {
        if (pi < p.size() && (p[pi] == '?' || p[pi] == n[ni])) { ++pi; ++ni; }
        else if (pi < p.size() && p[pi] == '*') { star = pi++; mark = ni; }
        else if (star != std::string::npos) { pi = star + 1; ni = ++mark; }
        else return false;
    }
    while (pi < p.size() && p[pi] == '*') ++pi;
    return pi == p.size();
}
}  // namespace

int Well::indexOf(int id) const {
    for (size_t i = 0; i < segs.size(); ++i) if (segs[i].id == id) return int(i);
    return -1;
}

std::vector<std::vector<int>> children(const Well& W) {
    std::vector<std::vector<int>> ch(W.segs.size());
    for (size_t i = 0; i < W.segs.size(); ++i) {
        const int o = W.indexOf(W.segs[i].outlet);
        if (o >= 0 && size_t(o) != i) ch[size_t(o)].push_back(int(i));
    }
    return ch;
}
std::vector<int> topoOrder(const Well& W) {
    std::vector<int> ord;
    if (W.segs.empty()) return ord;
    const auto ch = children(W);
    std::vector<char> seen(W.segs.size(), 0);
    int top = W.indexOf(1); if (top < 0) top = 0;
    ord.push_back(top); seen[size_t(top)] = 1;
    for (size_t h = 0; h < ord.size(); ++h)
        for (int c : ch[size_t(ord[h])]) if (!seen[size_t(c)]) { seen[size_t(c)] = 1; ord.push_back(c); }
    return ord;   // segments not connected to the top are dropped (parse warns)
}

std::vector<Well> parse(const std::vector<ExtraKeyword>& kws, int step, bool field, std::string& err) {
    std::map<std::string, Well> wells;
    std::vector<std::string> order;
    // pass 1: latest WELSEGS / COMPSEGS per well
    for (const auto& e : kws) {
        if (e.section != DeckSection::Schedule || e.step > step) continue;
        if (e.name == "WELSEGS") {
            const auto R = records(e.body);
            if (R.empty() || R[0].empty()) { err += "WELSEGS without a header record; "; continue; }
            Well W; W.name = unq(R[0][0]);
            W.topDepth = num(R[0], 1, 0); W.topLength = num(R[0], 2, 0);
            const bool inc = !has(R[0], 4) || up(unq(R[0][4])) != "ABS";
            if (has(R[0], 5) && up(unq(R[0][5])) == "H--") W.frictionOn = false;
            W.accelOn = has(R[0], 5) && up(unq(R[0][5])) == "HFA";   // C36
            W.driftFlux = has(R[0], 6) && up(unq(R[0][6])) == "DF";
            Segment top; top.id = 1; top.branch = 1; top.outlet = 0; top.length = W.topLength; top.depth = W.topDepth;
            W.segs.push_back(top);
            for (size_t r = 1; r < R.size(); ++r) {
                const auto& x = R[r];
                if (x.empty()) break;
                const int s1 = inum(x, 0, 0), s2 = inum(x, 1, s1), br = inum(x, 2, 1), out = inum(x, 3, s1 - 1);
                const double L = num(x, 4, 0), D = num(x, 5, 0), dia = num(x, 6, 0.1), ro = num(x, 7, 1e-5);
                if (s1 < 2 || s2 < s1) { err += "WELSEGS " + W.name + ": bad segment range; "; continue; }
                const int oi = W.indexOf(out);
                if (oi < 0) { err += "WELSEGS " + W.name + ": outlet segment " + std::to_string(out) + " not defined before segment " + std::to_string(s1) + "; "; continue; }
                const double L0 = W.segs[size_t(oi)].length, D0 = W.segs[size_t(oi)].depth;
                const int n = s2 - s1 + 1;
                for (int s = s1; s <= s2; ++s) {
                    Segment g; g.id = s; g.branch = br; g.outlet = s == s1 ? out : s - 1; g.diam = dia; g.rough = ro;
                    g.volume = std::max(0.0, num(x, 9, 0.0));   // C37 item 10 (default: area x length)
                    const int m = s - s1 + 1;
                    if (inc) { g.length = L0 + m * L; g.depth = D0 + m * D; }
                    else { g.length = L0 + (L - L0) * double(m) / n; g.depth = D0 + (D - D0) * double(m) / n; }
                    const int ex = W.indexOf(s);
                    if (ex >= 0) W.segs[size_t(ex)] = g; else W.segs.push_back(g);
                }
            }
            if (!wells.count(W.name)) order.push_back(W.name);
            std::vector<Conn> keep = wells.count(W.name) ? wells[W.name].conns : std::vector<Conn>{};
            W.conns = keep;
            wells[W.name] = W;
        } else if (e.name == "COMPSEGS") {
            const auto R = records(e.body);
            if (R.empty() || R[0].empty()) continue;
            const std::string nm = unq(R[0][0]);
            if (!wells.count(nm)) { err += "COMPSEGS " + nm + " before its WELSEGS; "; continue; }
            Well& W = wells[nm];
            W.conns.clear();
            for (size_t r = 1; r < R.size(); ++r) {
                const auto& x = R[r];
                if (x.empty()) break;
                Conn c; c.i = inum(x, 0, 1) - 1; c.j = inum(x, 1, 1) - 1; c.k = inum(x, 2, 1) - 1; c.branch = inum(x, 3, 1);
                c.d1 = num(x, 4, 0); c.d2 = num(x, 5, c.d1);
                const int iseg = inum(x, 10, -1);
                if (iseg > 0) c.seg = iseg;
                // ECLIPSE item 7 / 8: a range of cells along a direction - expand
                const std::string dir = has(x, 6) ? up(unq(x[6])) : "";
                const int endR = inum(x, 7, -1);
                if (!dir.empty() && endR > 0) {
                    const int from = dir == "I" ? c.i : dir == "J" ? c.j : c.k, to = endR - 1;
                    const int n = std::abs(to - from) + 1, sg = to >= from ? 1 : -1;
                    const double len = (c.d2 - c.d1) / n;
                    for (int m = 0; m < n; ++m) {
                        Conn d = c;
                        (dir == "I" ? d.i : dir == "J" ? d.j : d.k) = from + sg * m;
                        d.d1 = c.d1 + m * len; d.d2 = d.d1 + len;
                        W.conns.push_back(d);
                    }
                } else W.conns.push_back(c);
            }
        }
    }
    // C36: WSEGDFPA drift-flux parameters (latest in force, all DF wells)
    DriftParams dfp;
    for (const auto& e : kws) {
        if (e.name != "WSEGDFPA" || e.step > step) continue;
        const auto R = records(e.body);
        if (R.empty() || R[0].empty()) continue;
        const auto& x = R[0];
        dfp.vdFactor = num(x, 0, 1.0); dfp.a1 = num(x, 1, 0.2); dfp.a2 = num(x, 2, 0.4);
        dfp.A = num(x, 3, 1.2); dfp.B = num(x, 4, 0.3); dfp.Fv = num(x, 5, 1.0);
        // C37 oil-water slip: items 7-11
        dfp.owFactor = num(x, 6, -1.0); dfp.owA = num(x, 7, 1.2); dfp.owB1 = num(x, 8, 0.4); dfp.owB2 = num(x, 9, 0.7); dfp.owN = num(x, 10, 2.0);
        if (!(dfp.a2 > dfp.a1) || dfp.A < 1 || dfp.B < 0 || dfp.B >= 1) { err += "WSEGDFPA: invalid parameters - defaults used; "; dfp = DriftParams(); }
        if (dfp.owFactor > 0 && (dfp.owA < 1 || !(dfp.owB2 > dfp.owB1) || dfp.owB1 < 0 || dfp.owB2 > 1 || dfp.owN < 0)) {
            err += "WSEGDFPA: invalid oil-water parameters (items 8-11) - defaults used; ";
            const double f = dfp.owFactor; const DriftParams d0; dfp.owA = d0.owA; dfp.owB1 = d0.owB1; dfp.owB2 = d0.owB2; dfp.owN = d0.owN; dfp.owFactor = f;
        }
    }
    for (auto& [nm, W] : wells) W.df = dfp;
    // pass 2: devices
    for (const auto& e : kws) {
        if (e.section != DeckSection::Schedule || e.step > step) continue;
        if (e.name != "WSEGVALV" && e.name != "WSEGSICD" && e.name != "WSEGAICD" && e.name != "WSEGLABY" && e.name != "WSEGFLIM" && e.name != "WSEGPULL") continue;
        for (const auto& x : records(e.body)) {
            if (x.empty()) break;
            const std::string pat = unq(x[0]);
            for (auto& [nm, W] : wells) {
                if (!wild(pat, nm)) continue;
                const int s1 = inum(x, 1, 0), s2 = e.name == "WSEGVALV" ? s1 : inum(x, 2, s1);
                for (int s = s1; s <= s2; ++s) {
                    const int i = W.indexOf(s);
                    if (i < 0) { err += e.name + " " + nm + ": segment " + std::to_string(s) + " not defined; "; continue; }
                    Segment& g = W.segs[size_t(i)];
                    if (e.name == "WSEGFLIM") {   // C37: WELL SEG1 SEG2 PHASE1 LIMIT1 A1 PHASE2 LIMIT2 A2 STATUS
                        auto ph = [](const std::string& t) { const std::string u = up(t); return u == "OIL" ? 1 : u == "WAT" || u == "WATER" ? 2 : u == "GAS" ? 3 : u == "LIQ" || u == "LIQUID" ? 4 : 0; };
                        for (int l = 0; l < 2; ++l) {
                            g.flimPh[l] = has(x, size_t(3 + 3 * l)) ? ph(unq(x[size_t(3 + 3 * l)])) : 0;
                            g.flimQ[l] = num(x, size_t(4 + 3 * l), 0.0); g.flimA[l] = std::max(1.0, num(x, size_t(5 + 3 * l), 8.0));
                            if (g.flimPh[l] && g.flimQ[l] <= 0) { err += "WSEGFLIM " + nm + ": limit " + std::to_string(l + 1) + " not positive - ignored; "; g.flimPh[l] = 0; }
                        }
                        if (has(x, 9) && up(unq(x[9])) == "SHUT") g.shut = true;
                        continue;
                    }
                    if (e.name == "WSEGPULL") {   // C37: WELL SEG1 SEG2 P DP OIL_RATE_LIMIT GRADIENT_FACTOR
                        g.pump = true; g.pumpP = num(x, 3, 1.0); g.pumpDp = num(x, 4, 0.0); g.pumpQo = std::max(0.0, num(x, 5, 0.0)); g.pumpG = std::max(0.0, num(x, 6, 0.0));
                        if (g.pumpDp <= 0 || g.pumpP <= 0) { g.pump = false; err += "WSEGPULL " + nm + ": zero boost - ignored; "; }
                        continue;
                    }
                    if (e.name == "WSEGVALV") {
                        g.device = 1; g.cv = num(x, 2, 0); g.area = num(x, 3, 0);
                        g.addLen = num(x, 4, 0); g.vDiam = num(x, 5, g.diam); g.vRough = num(x, 6, g.rough);
                        g.shut = has(x, 8) && up(unq(x[8])) == "SHUT";
                        if (g.area <= 0 && !g.shut) { err += "WSEGVALV " + nm + ": zero area - treated as SHUT; "; g.shut = true; }
                    } else if (e.name == "WSEGLABY") {   // C36: items 4-14 (defaults metric, converted for FIELD)
                        const double fl = field ? 1.0 / 0.3048 : 1.0;
                        g.device = 4; g.labN = std::max(1, inum(x, 3, 1)); g.labCh = std::max(1, inum(x, 4, 2));
                        g.labA = num(x, 5, 6e-5 * fl * fl); g.labL1 = num(x, 6, 0.1863 * fl); g.labL2 = num(x, 7, 0.2832 * fl);
                        g.labD = num(x, 8, 0.006316 * fl); g.labR = num(x, 9, 1e-5 * fl);
                        g.gIn = num(x, 10, 0.35); g.gOut = num(x, 11, 0.5); g.gLab = num(x, 12, 0.5);
                        g.shut = has(x, 13) && up(unq(x[13])) == "SHUT";
                        if ((g.labA <= 0 || g.labD <= 0) && !g.shut) { err += "WSEGLABY " + nm + ": zero channel area / diameter - treated as SHUT; "; g.shut = true; }
                    } else if (e.name == "WSEGSICD") {
                        g.device = 2; g.strength = num(x, 3, 0);
                        g.rhoCal = num(x, 5, field ? 62.4 : 1000.25); g.muCal = num(x, 6, 1.45); g.x = 2; g.y = 0.25;
                        g.shut = has(x, 11) && up(unq(x[11])) == "SHUT";
                    } else {
                        g.device = 3; g.strength = num(x, 3, 0);
                        g.rhoCal = num(x, 5, field ? 62.4 : 1000.25); g.muCal = num(x, 6, 1.45);
                        g.x = num(x, 12, 1.0); g.y = num(x, 13, 1.0);   // items 13 / 14: FLOW_RATE_EXPONENT, VISC_EXPONENT (OPM keyword JSON)
                        g.shut = has(x, 20) && up(unq(x[20])) == "SHUT";
                    }
                }
            }
        }
    }
    std::vector<Well> out;
    for (const auto& nm : order) {
        Well W = wells[nm];
        const auto ord = topoOrder(W);
        if (ord.size() != W.segs.size()) W.warn += std::to_string(W.segs.size() - ord.size()) + " segment(s) not connected to segment 1; ";
        // allocate connections: the node of the branch nearest to the perforation centre (ECLIPSE rule)
        for (auto& c : W.conns) {
            if (c.seg > 0 && W.indexOf(c.seg) >= 0) continue;
            const double mid = 0.5 * (c.d1 + c.d2);
            double best = 1e300; int bs = -1;
            for (const auto& g : W.segs) if (g.branch == c.branch && g.id != 1 && std::fabs(g.length - mid) < best) { best = std::fabs(g.length - mid); bs = g.id; }
            if (bs < 0) for (const auto& g : W.segs) if (std::fabs(g.length - mid) < best) { best = std::fabs(g.length - mid); bs = g.id; }
            c.seg = bs;
        }
        out.push_back(W);
    }
    return out;
}

double frictionFactor(double Re, double relRough) {
    if (Re <= 0) return 0;
    if (Re < 2300) return 64.0 / Re;
    const double t = -1.8 * std::log10(std::pow(relRough / 3.7, 1.11) + 6.9 / Re);
    return 1.0 / (t * t);
}

double labyrinthDp(const Segment& g, double qSI, double rhoSI, double muSI, double toM) {
    if (qSI <= 0 || g.labA <= 0 || g.labD <= 0) return 0.0;
    const double A = g.labA * toM * toM, D = g.labD * toM, L = (g.labL1 + (g.labCh - 1) * g.labL2) * toM;
    const double v = qSI / (g.labN * A), Re = rhoSI * v * D / std::max(muSI, 1e-6);
    const double K = g.gIn + g.gOut + (g.labCh - 1) * g.gLab + frictionFactor(Re, g.labR * toM / D) * L / D;
    return 0.5 * rhoSI * v * v * K;
}

double driftHoldup(double jg, double jl, double rhoG, double rhoL, double D, double theta, const DriftParams& P, double* C0o, double* vdo) {
    const double j = jg + jl;
    if (C0o) *C0o = 1.0;
    if (vdo) *vdo = 0.0;
    if (jg <= 0) return 0.0;
    if (jl <= 0 || !(rhoL > rhoG) || rhoG <= 0 || D <= 0) return jg / std::max(j, 1e-300);
    const double g = 9.80665, dr = rhoL - rhoG;
    const double Vc = std::pow(P.sigma * g * dr / (rhoL * rhoL), 0.25);
    const double Db = D * std::sqrt(g * dr / P.sigma);   // dimensionless diameter (sqrt of the bond number)
    const double CK = 142.0, Cw = 0.008;
    const double Ku = std::sqrt(CK / Db * (std::sqrt(1.0 + Db * Db / (CK * CK * Cw)) - 1.0));
    const double vSgf = Ku * std::sqrt(rhoL / rhoG) * Vc;   // flooding gas superficial velocity
    const double c = std::cos(theta), sn = std::sin(theta);
    const double m = std::fabs(c) < 1e-9 ? 0.0 : (c >= 0 ? 1.0 : -1.0) * std::sqrt(std::fabs(c)) * std::pow(1.0 + std::max(sn, 0.0), 1.2) * P.vdFactor;   // Hasan-Kabir
    auto eval = [&](double a, double& C0, double& vd) {
        const double beta = std::max(a, P.Fv * a * std::fabs(j) / vSgf);
        const double gam = std::min(1.0, std::max(0.0, (beta - P.B) / (1.0 - P.B)));
        C0 = P.A / (1.0 + (P.A - 1.0) * gam * gam);
        const double K = a <= P.a1 ? 1.53 / C0 : a >= P.a2 ? Ku : 1.53 / C0 + (Ku - 1.53 / C0) * (a - P.a1) / (P.a2 - P.a1);
        const double den = a * C0 * std::sqrt(rhoG / rhoL) + 1.0 - a * C0;
        vd = den > 1e-12 ? std::max(0.0, 1.0 - a * C0) * C0 * K * Vc * m / den : 0.0;
        return a * (C0 * j + vd) - jg;
    };
    double lo = 0, hi = 1, C0 = 1, vd = 0;
    for (int it = 0; it < 80; ++it) {
        const double a = 0.5 * (lo + hi);
        if (eval(a, C0, vd) > 0) hi = a; else lo = a;
    }
    const double a = 0.5 * (lo + hi);
    eval(a, C0, vd);
    if (C0o) *C0o = C0;
    if (vdo) *vdo = vd;
    return a;
}

// C37: Shi et al. (2005) oil-water drift flux inside the liquid (oil holdup of the liquid volume)
double owHoldup(double jo, double jw, double rhoO, double rhoW, double theta, const DriftParams& P, double* C0o, double* vdo) {
    if (C0o) *C0o = 1.0;
    if (vdo) *vdo = 0.0;
    const double j = jo + jw;
    if (j <= 0 || jo <= 0) return 0.0;
    if (jw <= 0) return 1.0;
    const double g = 9.80665, dr = rhoW - rhoO;
    const double Vc = dr > 0 && rhoW > 0 ? std::pow(P.owSigma * g * dr / (rhoW * rhoW), 0.25) : 0.0;
    const double c = std::cos(theta), sn = std::sin(theta);
    const double m = std::fabs(c) < 1e-9 ? 0.0 : (c >= 0 ? 1.0 : -1.0) * std::sqrt(std::fabs(c)) * std::pow(1.0 + std::max(sn, 0.0), 1.2) * std::max(0.0, P.owFactor);
    auto eval = [&](double a, double& C0, double& vd) {
        const double gam = std::min(1.0, std::max(0.0, (a - P.owB1) / std::max(1e-9, P.owB2 - P.owB1)));
        C0 = P.owA / (1.0 + (P.owA - 1.0) * gam * gam);
        vd = 1.53 * Vc * std::pow(std::max(0.0, 1.0 - a), P.owN) * m;
        return a * (C0 * j + vd) - jo;
    };
    double lo = 0, hi = 1, C0 = 1, vd = 0;
    for (int it = 0; it < 80; ++it) {
        const double a = 0.5 * (lo + hi);
        if (eval(a, C0, vd) > 0) hi = a; else lo = a;
    }
    const double a = 0.5 * (lo + hi);
    eval(a, C0, vd);
    if (C0o) *C0o = C0;
    if (vdo) *vdo = vd;
    return a;
}

double flowLimitFactor(double q0, double qlim, double a) {
    if (!(q0 > 0) || !(qlim > 0)) return 1.0;
    a = std::max(1.0, a);
    const double r = q0 / qlim, la = a * std::log(r);
    if (la < -60) return 1.0;
    if (la > 60) return 1.0 / r;
    return std::pow(1.0 + std::exp(la), -1.0 / a);
}

Solution solve(const Well& W, double pTop, const Inflow& in, const Units& u, const std::vector<double>& nonHydroPrev, const SolveOpts& opts) {
    const size_t n = W.segs.size();
    Solution S;
    S.p.assign(n, pTop); S.qSeg.assign(n, 0); S.rhoSeg.assign(n, 0); S.muSeg.assign(n, 1);
    S.dHydro.assign(n, 0); S.dFric.assign(n, 0); S.dDev.assign(n, 0); S.nonHydro.assign(n, 0); S.devFac.assign(n, 1.0);
    S.dAcc.assign(n, 0); S.holdup.assign(n, 0); S.noSlip.assign(n, 0); S.rhoHead.assign(n, 0);
    S.oilHold.assign(n, 0); S.dPump.assign(n, 0); S.flimFac.assign(n, 1.0); S.qoLoc.assign(n, 0); S.qwLoc.assign(n, 0); S.qgLoc.assign(n, 0); S.flimQ0.assign(n, 0.0);
    if (n == 0) return S;
    const auto ord = topoOrder(W);
    // subtree sums (children before parents = reverse topological order)
    std::vector<double> q(n, 0), rq(n, 0), mq(n, 0), A(n, 0), AP(n, 0), absq(n, 0), qg(n, 0), rgq(n, 0), rlq(n, 0);
    const bool phases = in.qg.size() == n && in.rhoG.size() == n && in.rhoL.size() == n;
    const bool fl = bool(in.flash) && in.qs.size() == 3 * n;   // C37 tubing flash data
    std::vector<double> qs(3 * n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        q[i] = in.q[i]; absq[i] = std::fabs(in.q[i]);
        rq[i] = absq[i] * in.rho[i]; mq[i] = absq[i] * in.mu[i];
        A[i] = in.a[i]; AP[i] = in.a[i] * in.pres[i];
        if (phases) { qg[i] = in.qg[i]; rgq[i] = std::fabs(in.qg[i]) * in.rhoG[i]; rlq[i] = std::fabs(in.q[i] - in.qg[i]) * in.rhoL[i]; }
        if (fl) for (int c = 0; c < 3; ++c) qs[3 * i + size_t(c)] = in.qs[3 * i + size_t(c)];
    }
    for (size_t h = ord.size(); h-- > 1;) {
        const int s = ord[h], o = W.indexOf(W.segs[size_t(s)].outlet);
        if (o < 0) continue;
        q[size_t(o)] += q[size_t(s)]; absq[size_t(o)] += absq[size_t(s)]; rq[size_t(o)] += rq[size_t(s)]; mq[size_t(o)] += mq[size_t(s)];
        A[size_t(o)] += A[size_t(s)]; AP[size_t(o)] += AP[size_t(s)];
        qg[size_t(o)] += qg[size_t(s)]; rgq[size_t(o)] += rgq[size_t(s)]; rlq[size_t(o)] += rlq[size_t(s)];
        for (int c = 0; c < 3; ++c) qs[3 * size_t(o) + size_t(c)] += qs[3 * size_t(s) + size_t(c)];
    }
    // densities with no flow: nearest flowing / own cell values (walk down from the top keeps the parent's value)
    for (size_t h = 0; h < ord.size(); ++h) {
        const size_t s = size_t(ord[h]);
        const int o = W.indexOf(W.segs[s].outlet);
        if (absq[s] > 0) { S.rhoSeg[s] = rq[s] / absq[s]; S.muSeg[s] = std::max(1e-3, mq[s] / absq[s]); }
        else if (in.rho[s] > 0) { S.rhoSeg[s] = in.rho[s]; S.muSeg[s] = std::max(1e-3, in.mu[s]); }
        else if (o >= 0) { S.rhoSeg[s] = S.rhoSeg[size_t(o)]; S.muSeg[s] = S.muSeg[size_t(o)]; }
        S.qSeg[s] = q[s];
    }
    // a top segment without its own density takes the first child's
    {
        const size_t t = size_t(ord[0]);
        if (S.rhoSeg[t] <= 0) for (size_t h = 1; h < ord.size(); ++h) if (S.rhoSeg[size_t(ord[h])] > 0) { S.rhoSeg[t] = S.rhoSeg[size_t(ord[h])]; break; }
    }
    const double PI = 3.14159265358979323846;
    // local (in-situ) volumetric rates and phase densities: inflow-weighted reservoir values, replaced by the flash
    std::vector<double> qv(q), qgv(qg), qov(n, 0.0), qwv(n, 0.0), rG(n, 0.0), rL(n, 0.0), rO(n, 0.0), rW(n, 0.0);
    for (size_t s = 0; s < n; ++s) {
        const double ql = std::fabs(q[s] - qg[s]);
        rG[s] = std::fabs(qg[s]) > 0 ? rgq[s] / std::fabs(qg[s]) : 0.0;
        rL[s] = ql > 0 ? rlq[s] / ql : 0.0;
    }
    const auto ch = children(W);
    const bool sto = opts.sto && opts.sto->size() == 2 * n && opts.dt > 0;
    const int passes = fl ? 3 : 1;
    std::vector<double> mom(n, 0.0);
    for (int pass = 0; pass < passes; ++pass) {
        if (pass > 0) {   // C37 tubing flash at the segment pressure of the previous pass (or the subtree reservoir pressure)
            for (size_t h = 0; h < ord.size(); ++h) {
                const size_t s = size_t(ord[h]);
                const double qq[3] = {std::fabs(qs[3 * s]), std::fabs(qs[3 * s + 1]), std::fabs(qs[3 * s + 2])};
                if (qq[0] + qq[1] + qq[2] <= 0) continue;
                const int o = W.indexOf(W.segs[s].outlet);
                double pe = S.p[s];
                if (!in.flashAtSegment) pe = A[s] > 0 ? AP[s] / A[s] : S.p[s];
                else if (h > 0 && o >= 0) pe = 0.5 * (S.p[s] + S.p[size_t(o)]);
                const Inflow::Flash F = in.flash(pe, qq);
                const double vol = F.qo + F.qw + F.qg;
                if (!(vol > 0)) continue;
                const double sg = q[s] >= 0 ? 1.0 : -1.0;
                qv[s] = sg * vol; qgv[s] = sg * F.qg; qov[s] = sg * F.qo; qwv[s] = sg * F.qw;
                rG[s] = F.rhoG; rO[s] = F.rhoO; rW[s] = F.rhoW;
                const double ql = F.qo + F.qw;
                rL[s] = ql > 0 ? (F.qo * F.rhoO + F.qw * F.rhoW) / ql : 0.0;
                S.rhoSeg[s] = (F.qo * F.rhoO + F.qw * F.rhoW + F.qg * F.rhoG) / vol;
                S.muSeg[s] = std::max(1e-3, (F.qo * F.muO + F.qw * F.muW + F.qg * F.muG) / vol);
            }
        }
        // C36 holdup (drift flux), head density and axial momentum flux per segment; C37 oil-water slip
        const bool havePh = phases || (fl && pass > 0);
        for (size_t h = 0; h < ord.size(); ++h) {
            const size_t s = size_t(ord[h]);
            const Segment& g = W.segs[s];
            S.rhoHead[s] = S.rhoSeg[s];
            S.qSeg[s] = qv[s];
            S.qoLoc[s] = qov[s]; S.qwLoc[s] = qwv[s]; S.qgLoc[s] = qgv[s];
            const double Dm = g.diam * u.toM, Ap = 0.25 * PI * Dm * Dm;
            const double qSI = qv[s] * u.toQ, rhoSI = S.rhoSeg[s] * u.toRho;
            mom[s] = Ap > 0 ? rhoSI * qSI * qSI / Ap : 0.0;
            S.holdup[s] = S.noSlip[s] = 0; S.oilHold[s] = 0;
            if (!havePh || qv[s] == 0) continue;
            const double fg = std::min(1.0, std::max(0.0, qgv[s] / qv[s]));
            S.noSlip[s] = S.holdup[s] = fg;
            const double ql = std::fabs(qv[s] - qgv[s]);
            S.oilHold[s] = ql > 0 ? std::min(1.0, std::fabs(qov[s]) / ql) : 0.0;
            if (!W.driftFlux || h == 0 || qv[s] <= 0 || Ap <= 0) continue;
            const int oi = W.indexOf(g.outlet);
            if (oi < 0) continue;
            const Segment& go = W.segs[size_t(oi)];
            const double dzUp = g.depth - go.depth, dL = std::fabs(g.length - go.length);   // flow s -> outlet rises if dzUp > 0
            const double ct = dL > std::fabs(dzUp) && dL > 0 ? dzUp / dL : (dzUp > 0 ? 1.0 : dzUp < 0 ? -1.0 : 0.0);
            const double theta = std::acos(std::max(-1.0, std::min(1.0, ct)));
            const double jg = std::fabs(qgv[s]) * u.toQ / Ap, jl = ql * u.toQ / Ap;
            double al = fg;
            if (fg > 0 && fg < 1 && rL[s] > rG[s] && rG[s] > 0) al = driftHoldup(jg, jl, rG[s] * u.toRho, rL[s] * u.toRho, Dm, theta, W.df);
            double rLin = rL[s];
            if (W.df.owFactor > 0 && std::fabs(qov[s]) > 0 && std::fabs(qwv[s]) > 0 && rW[s] > rO[s] && rO[s] > 0) {
                const double ao = owHoldup(std::fabs(qov[s]) * u.toQ / Ap, std::fabs(qwv[s]) * u.toQ / Ap, rO[s] * u.toRho, rW[s] * u.toRho, theta, W.df);
                S.oilHold[s] = ao;
                rLin = ao * rO[s] + (1 - ao) * rW[s];
            }
            const double rh = al * rG[s] + (1 - al) * rLin;
            if (!(rh > 0)) continue;
            S.holdup[s] = al;
            S.rhoHead[s] = rh;
            const double rGs = rG[s] * u.toRho, rLs = rLin * u.toRho;
            mom[s] = ((al > 1e-9 ? rGs * jg * jg / al : 0.0) + (al < 1 - 1e-9 ? rLs * jl * jl / (1 - al) : 0.0)) * Ap;
        }
        // C37 segment storage: the segment content moves toward the steady value with the volume swept in the step
        if (sto) {
            for (size_t h = 1; h < ord.size(); ++h) {
                const size_t s = size_t(ord[h]);
                const double pr = (*opts.sto)[s], ph = (*opts.sto)[n + s];
                if (!(pr > 0)) continue;
                const Segment& g = W.segs[s];
                const int oi = W.indexOf(g.outlet);
                const double Dm = g.diam * u.toM, Ap = 0.25 * PI * Dm * Dm;
                const double V = g.volume > 0 ? g.volume * u.toM * u.toM * u.toM : Ap * (oi >= 0 ? std::fabs(g.length - W.segs[size_t(oi)].length) : 0.0) * u.toM;
                if (!(V > 0)) continue;
                const double xs = opts.xfSweep && opts.xfSweep->size() == n ? (*opts.xfSweep)[s] : 0.0;   // C45 downward cross-flow
                const double sw = (std::fabs(qv[s]) + xs) * u.toQ * opts.dt * 86400.0, w = sw / (V + sw);
                S.rhoHead[s] = pr + w * (S.rhoHead[s] - pr);
                S.holdup[s] = std::min(1.0, std::max(0.0, ph + w * (S.holdup[s] - ph)));
            }
        }
        // pressures from the top node down
        for (size_t h = 1; h < ord.size(); ++h) {
            const size_t s = size_t(ord[h]);
            const Segment& g = W.segs[s];
            const int oi = W.indexOf(g.outlet);
            const size_t o = size_t(std::max(oi, 0));
            const Segment& go = W.segs[o];
            S.dHydro[s] = S.rhoHead[s] * u.g * (g.depth - go.depth);
            const double sg = q[s] >= 0 ? 1.0 : -1.0;   // flow toward the outlet: pressure rises away from the top
            const double rhoSI = S.rhoSeg[s] * u.toRho, muSI = S.muSeg[s] * 1e-3, qSI = std::fabs(qv[s]) * u.toQ;
            double fr = 0;
            auto pipe = [&](double L, double D, double eps) {
                if (L <= 0 || D <= 0 || qSI <= 0) return 0.0;
                const double Ar = 0.25 * PI * D * D, v = qSI / Ar, Re = rhoSI * v * D / muSI;
                return frictionFactor(Re, eps / D) * rhoSI * v * v * L / (2 * D) / u.toPa;
            };
            if (W.frictionOn) {
                fr = pipe(std::fabs(g.length - go.length) * u.toM, g.diam * u.toM, g.rough * u.toM);
                if (g.device == 1 && g.addLen > 0) fr += pipe(g.addLen * u.toM, (g.vDiam > 0 ? g.vDiam : g.diam) * u.toM, g.vRough * u.toM);
            }
            S.dFric[s] = sg * fr;
            S.dAcc[s] = 0;
            if (W.accelOn) {   // C36: momentum flux leaving through the outlet face minus that entering from the children
                const double Dm = g.diam * u.toM, Ap = 0.25 * PI * Dm * Dm;
                double mIn = 0;
                for (int c : ch[s]) mIn += mom[size_t(c)];
                if (Ap > 0) S.dAcc[s] = (mom[s] - mIn) / Ap / u.toPa;
            }
            // C37 WSEGPULL: pump boost for producing flow (lowers the upstream pressure)
            S.dPump[s] = 0;
            if (g.pump && q[s] > 0) {
                const double qo = fl ? std::fabs(qs[3 * s]) : 0.0;
                S.dPump[s] = -g.pumpP * std::max(0.0, g.pumpDp - g.pumpG * std::max(0.0, qo - g.pumpQo));
            }
            S.p[s] = S.p[o] + S.dHydro[s];   // provisional (for the device drawdown)
            S.nonHydro[s] = S.nonHydro[o] + S.dFric[s] + S.dAcc[s] + S.dPump[s];
            S.p[s] += S.dFric[s] + S.dAcc[s] + S.dPump[s];
            S.devFac[s] = S.devFac[o];
            S.flimFac[s] = 1.0;
            S.dDev[s] = 0;
            const double pRes = A[s] > 0 ? AP[s] / A[s] : 0.0;
            const double D0 = A[s] > 0 ? (sg > 0 ? pRes - S.p[s] : S.p[s] - pRes) : 0.0;   // drawdown of the subtree without the device
            double dp = 0;
            if (g.device > 0 && !g.shut && A[s] > 0) {
                // c in case pressure / (case rate)^x; C36 labyrinth: general dp(q)
                double c = 0, x = 2;
                if (g.device == 4) {
                    if (g.labA > 0 && g.labD > 0) c = 1.0;
                } else if (g.device == 1) {
                    if (g.cv > 0 && g.area > 0) {
                        const double Ac = g.area * u.toM * u.toM;
                        c = rhoSI / (2 * g.cv * g.cv * Ac * Ac) * u.toQ * u.toQ / u.toPa;
                    }
                } else if (g.device == 2) {
                    const double rr = g.rhoCal > 0 ? S.rhoSeg[s] / g.rhoCal : 1.0, mr = g.muCal > 0 ? S.muSeg[s] / g.muCal : 1.0;
                    c = g.strength * std::pow(std::max(rr, 0.0), 0.75) * std::pow(std::max(mr, 1e-9), 0.25);
                } else {
                    const double rr = g.rhoCal > 0 ? S.rhoSeg[s] / g.rhoCal : 1.0, mr = S.muSeg[s] > 0 ? g.muCal / S.muSeg[s] : 1.0;
                    c = g.strength * rr * std::pow(std::max(mr, 1e-9), g.y); x = std::max(0.1, g.x);
                }
                if (c > 0 && D0 > 0) {   // dp = c (A (D0 - dp))^x, monotone: bisection
                    double lo = 0, hi = D0;
                    for (int it = 0; it < 80; ++it) {
                        const double m = 0.5 * (lo + hi);
                        const double qd = A[s] * (D0 - m);
                        const double dq = g.device == 4 ? labyrinthDp(g, qd * u.toQ, rhoSI, muSI, u.toM) / u.toPa : c * std::pow(qd, x);
                        if (m - dq > 0) hi = m; else lo = m;
                    }
                    dp = 0.5 * (lo + hi);
                }
                S.dDev[s] = sg * dp;
                if (D0 > 0) S.devFac[s] *= std::max(1e-6, (D0 - dp) / D0);
            }
            // C37 WSEGFLIM: smooth cut of the subtree's potential phase rate to the limit
            const bool prevF = opts.flimPrev && opts.flimPrev->size() == n && (*opts.flimPrev)[s] > 0;
            if (!g.shut && fl && (g.flimPh[0] > 0 || g.flimPh[1] > 0) && std::fabs(q[s]) > 0 && (prevF || (A[s] > 0 && D0 > dp))) {
                auto phaseRate = [&](int ph) {
                    return ph == 1 ? std::fabs(qs[3 * s]) : ph == 2 ? std::fabs(qs[3 * s + 1]) : ph == 3 ? std::fabs(qs[3 * s + 2]) : std::fabs(qs[3 * s]) + std::fabs(qs[3 * s + 1]);
                };
                auto factor = [&](double qres0) {
                    double f = 1.0;
                    for (int l = 0; l < 2; ++l) if (g.flimPh[l] > 0) f *= flowLimitFactor(phaseRate(g.flimPh[l]) / std::fabs(q[s]) * qres0, g.flimQ[l], g.flimA[l]);
                    return std::max(1e-6, f);
                };
                // potential without the limiter (subtree reservoir rate): first solve - the device-free drawdown potential;
                // then the current rate over the factor in force (fixed point q = f(Q0) Q0), extrapolated one step with the
                // growth of the potential since the last solve (the factor acts with a one-step lag)
                double qres0 = A[s] * (D0 - dp) * S.devFac[o];
                if (prevF) {
                    const double used = (*opts.flimPrev)[s];
                    const double now = std::fabs(q[s]) / factor(used);
                    qres0 = now * std::min(1.25, std::max(0.8, now / used));
                }
                const double f = factor(qres0);
                S.flimQ0[s] = qres0;
                S.flimFac[s] = f;
                S.devFac[s] *= f;
            }
        }
    }
    // under-relaxation of the non-hydrostatic part with the previous step
    if (nonHydroPrev.size() == n) {
        for (size_t h = 0; h < ord.size(); ++h) {
            const size_t s = size_t(ord[h]);
            const double nh = 0.5 * (S.nonHydro[s] + nonHydroPrev[s]);
            S.p[s] += nh - S.nonHydro[s];
            S.nonHydro[s] = nh;
        }
    }
    return S;
}

}  // namespace msw
}  // namespace ressim
