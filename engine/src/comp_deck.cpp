// Compositional deck hook-up (milestone C3) - see comp_deck.h.
#include <unordered_map>
#include "ressim/binio.h"
#include "ressim/comp_deck.h"
#include "ressim/deck_parser.h"
#include "ressim/grid.h"
#include "ressim/io.h"
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <sstream>
#include <map>
#include <stdexcept>

namespace ressim {
namespace comp {

namespace {
constexpr double GRAV = 9.80665;
constexpr double MWW = 18.01528e-3;
struct UnitsSI {
    double len, p, liq, gas, dens, trans, comp, pv;
    explicit UnitsSI(bool field) {
        if (field) {
            len = 0.3048; p = 6894.757293168; liq = 0.158987294928; gas = 28.316846592; dens = 16.01846337;
            trans = 0.158987294928 * 1e-3 / (86400.0 * 6894.757293168); comp = 1.0 / 6894.757293168; pv = 0.158987294928;   // rb
        } else {
            len = 1; p = 1e5; liq = 1; gas = 1; dens = 1; trans = 1e-3 / (86400.0 * 1e5); comp = 1e-5; pv = 1;
        }
    }
};
std::string up(std::string s) { for (char& ch : s) ch = char(std::toupper(static_cast<unsigned char>(ch))); return s; }
double lin(const std::vector<double>& x, const std::vector<double>& y, double v) {
    if (x.empty()) return 0;
    if (v <= x.front()) return y.front();
    if (v >= x.back()) return y.back();
    const size_t i = size_t(std::upper_bound(x.begin(), x.end(), v) - x.begin());
    const double t = (v - x[i - 1]) / std::max(1e-300, x[i] - x[i - 1]);
    return y[i - 1] + t * (y[i] - y[i - 1]);
}
// composition at a depth from the ZMFVD table (linear, clamped) or ZI
std::vector<double> zAt(const eos::CompFluid& fl, double d) {
    if (!fl.zmfvd.empty() && !fl.zmfvd[0].depth.empty()) {
        const eos::ZmfTable& t = fl.zmfvd[0];
        if (d <= t.depth.front()) return t.z.front();
        if (d >= t.depth.back()) return t.z.back();
        const size_t i = size_t(std::upper_bound(t.depth.begin(), t.depth.end(), d) - t.depth.begin());
        const double w = (d - t.depth[i - 1]) / std::max(1e-300, t.depth[i] - t.depth[i - 1]);
        std::vector<double> z(t.z[i].size());
        for (size_t k = 0; k < z.size(); ++k) z[k] = (1 - w) * t.z[i - 1][k] + w * t.z[i][k];
        return z;
    }
    return fl.zi;
}
}  // namespace

// ------------------------------------------------------------------------------------------------- equilibration
bool equilibrate(const Case& c, const eos::CompFluid& fl, CompCase& cc, std::string* error) {
    const UnitsSI U(c.units == UnitSystem::Field);
    const int n = cc.size();
    if (fl.zmfvd.empty() && fl.zi.empty()) { if (error) *error = "initial composition: ZMFVD or ZI is required"; return false; }
    if (c.init.active()) {   // enumerated PRESSURE / SWAT with ZMFVD / ZI compositions
        cc.p0.resize(size_t(n)); cc.sw0.resize(size_t(n)); cc.z0.resize(size_t(n));
        return true;   // filled by the caller (needs the cell -> global map)
    }
    const EquilInput& E = c.equil;
    const double D0 = E.datumDepth * U.len, P0 = E.datumPressure * U.p, Dw = E.owc * U.len, pcw = E.pcowOWC * U.p;
    const eos::TParams tp = eos::atTemperature(cc.fluid, cc.T);
    double dmin = std::min(D0, Dw), dmax = std::max(D0, Dw);
    for (double d : cc.cellDepth) { dmin = std::min(dmin, d); dmax = std::max(dmax, d); }
    dmin -= 1.0; dmax += 1.0;
    const int NN = 400;
    const double h = (dmax - dmin) / NN;
    // mixture density of the local hydrocarbon at (p, depth)
    auto rhoHC = [&](double p, double d) {
        const std::vector<double> z = zAt(fl, d);
        const eos::FlashResult f = eos::flash(cc.fluid, tp, p, z);
        if (f.phases == 2) {
            const double vL = (1 - f.beta) * f.L.vm, vV = f.beta * f.V.vm;
            return ((1 - f.beta) * f.L.mw + f.beta * f.V.mw) / (vL + vV);
        }
        return f.vapourLabel ? f.V.rho : f.L.rho;
    };
    auto rhoW = [&](double p) {
        const double X = cc.cw * (p - cc.pwRef), bw = cc.bwRef / (1 + X + 0.5 * X * X);
        return cc.rhoWsc / bw;
    };
    // integrate dp/dD = rho g from (Dref, pref) over the node grid (RK2 / Heun)
    auto integrate = [&](double Dref, double pref, const std::function<double(double, double)>& rho) {
        std::vector<double> dep(size_t(NN) + 1), pr(size_t(NN) + 1);
        for (int i = 0; i <= NN; ++i) dep[size_t(i)] = dmin + h * i;
        int i0 = int(std::floor((Dref - dmin) / h)); i0 = std::clamp(i0, 0, NN);
        auto step = [&](double d, double p, double dd) {
            const double k1 = rho(p, d) * GRAV, k2 = rho(p + dd * k1, d + dd) * GRAV;
            return p + 0.5 * dd * (k1 + k2);
        };
        pr[size_t(i0)] = step(Dref, pref, dep[size_t(i0)] - Dref);
        for (int i = i0 + 1; i <= NN; ++i) pr[size_t(i)] = step(dep[size_t(i - 1)], pr[size_t(i - 1)], h);
        for (int i = i0 - 1; i >= 0; --i) pr[size_t(i)] = step(dep[size_t(i + 1)], pr[size_t(i + 1)], -h);
        return std::make_pair(dep, pr);
    };
    const auto po = integrate(D0, P0, [&](double p, double d) { return rhoHC(p, d); });
    const double poOWC = lin(po.first, po.second, Dw);
    const auto pw = integrate(Dw, poOWC - pcw, [&](double p, double) { return rhoW(p); });
    cc.p0.assign(size_t(n), 0); cc.sw0.assign(size_t(n), cc.swc); cc.z0.assign(size_t(n), {});
    const SatTables& S = cc.sat;
    bool pcTable = false; for (double v : S.pcow) pcTable = pcTable || std::fabs(v) > 0;
    const double swmax = S.sw.empty() ? 1.0 : S.sw.back();
    for (int q = 0; q < n; ++q) {
        const double d = cc.cellDepth[size_t(q)];
        const double p = lin(po.first, po.second, d), pwq = lin(pw.first, pw.second, d);
        cc.p0[size_t(q)] = p;
        cc.z0[size_t(q)] = zAt(fl, d);
        if (!cc.water) continue;
        double sw;
        if (!pcTable) sw = d > Dw ? swmax : cc.swc;
        else {   // invert Pcow(Sw) = p_o - p_w (Pcow non-increasing in Sw)
            const double pc = p - pwq;
            if (pc >= S.pcow.front()) sw = S.sw.front();
            else if (pc <= S.pcow.back()) sw = S.sw.back();   // below the free-water level of the table
            else {
                sw = S.sw.back();
                for (size_t i = 0; i + 1 < S.pcow.size(); ++i)
                    if (S.pcow[i] >= pc && S.pcow[i + 1] < pc) {
                        const double t = (S.pcow[i] - pc) / (S.pcow[i] - S.pcow[i + 1]);
                        sw = S.sw[i] + t * (S.sw[i + 1] - S.sw[i]);
                        break;
                    }
            }
        }
        cc.sw0[size_t(q)] = std::clamp(sw, 0.0, 1.0);
        // water zone (Sw at the table maximum): the oil phase is absent - the cell pressure follows the water gradient,
        // p_o = p_w + Pcow(Sw), so that the mobile water is in hydrostatic equilibrium
        const double pcS = pcTable ? lin(S.sw, S.pcow, cc.sw0[size_t(q)]) : 0.0;
        if (p - pwq < pcS || (!pcTable && d > Dw)) cc.p0[size_t(q)] = pwq + pcS;
    }
    return true;
}

// CO2STORE hydrostatic initialisation: brine below the contact (GOC for OIL + GAS decks, where the "oil" is the brine;
// GWC = EQUIL item 3 for GAS + WATER decks), CO2 above it at the connate brine saturation with CO2-saturated brine.
// Pressures: the brine gradient (Batzle-Wang / Garcia) and the CO2 gradient (Redlich-Kwong) integrated from the datum
// (RK2); the primary pressure is the brine pressure, p_w = p_CO2 - Pc(Sw).
bool equilibrateCo2(const Case& c, CompCase& cc, std::string* error) {
    const UnitsSI U(c.units == UnitSystem::Field);
    const int n = cc.size();
    const EquilInput& E = c.equil;
    const double D0 = E.datumDepth * U.len, P0 = E.datumPressure * U.p;
    const double Dc = (c.phases.oil ? E.goc : E.owc) * U.len;
    const co2::Fluid cf = co2::makeFluid(cc.T, cc.salinity);
    auto rhoB = [&](double p) { return co2::aqueousDensity(cf, std::max(p, 1e5), 0.0); };
    auto rhoG = [&](double p) { return co2::co2Density(cf, std::max(p, 1e5)); };
    auto pAt = [&](double Dref, double pref, double d, const std::function<double(double)>& rho) {
        const int ns = std::clamp(int(std::fabs(d - Dref) / 2.0) + 1, 1, 400);
        const double h = (d - Dref) / ns;
        double p = pref;
        for (int i = 0; i < ns; ++i) { const double k1 = rho(p) * GRAV, k2 = rho(p + h * k1) * GRAV; p += 0.5 * h * (k1 + k2); }
        return p;
    };
    if (c.init.active()) { if (error) *error = "CO2STORE: enumerated initial states are not supported yet (use EQUIL)"; return false; }
    double pBc, pGc;   // brine / CO2 pressure at the contact
    if (D0 >= Dc) { pBc = pAt(D0, P0, Dc, rhoB); pGc = pBc; }
    else { pGc = pAt(D0, P0, Dc, rhoG); pBc = pGc; }
    const SatTables& S = cc.sat;
    const double swc = S.sw.empty() ? cc.swc : S.sw.front();
    const double pcC = (S.pcow.empty() ? 0.0 : lin(S.sw, S.pcow, swc)) + (S.pcog.empty() ? 0.0 : lin(S.sg, S.pcog, 1 - swc));
    cc.p0.assign(size_t(n), 0.0); cc.z0.assign(size_t(n), {0.0, 1.0}); cc.sw0.assign(size_t(n), 1.0);
    for (int q = 0; q < n; ++q) {
        const double d = cc.cellDepth[size_t(q)];
        if (d >= Dc) { cc.p0[size_t(q)] = (D0 >= Dc) ? pAt(D0, P0, d, rhoB) : pAt(Dc, pBc, d, rhoB); continue; }
        const double pg = (D0 < Dc) ? pAt(D0, P0, d, rhoG) : pAt(Dc, pGc, d, rhoG);
        const double pw = pg - pcC;
        double x, y; co2::solubility(cf, pw, x, y);
        const double na = swc / co2::aqueousMolarVolume(cf, pw, x), ng = (1 - swc) / co2::gasMolarVolume(cf, pw);
        const double zc = (na * x + ng * y) / (na + ng);
        cc.p0[size_t(q)] = pw; cc.z0[size_t(q)] = {zc, 1 - zc}; cc.sw0[size_t(q)] = swc;
    }
    return true;
}

// ------------------------------------------------------------------------------------------------- deck -> model
// GUIDERAT record from deck text: MINTIME PHASE A B C D E F ALLOW_INCREASE DAMPING ... (n* defaults, quotes, '--' comments)
static GuideRateSpec parseGuideRat(const std::string& body) {
    std::vector<std::string> t;
    std::string line; std::istringstream is(body);
    bool end = false;
    while (!end && std::getline(is, line)) {
        const size_t cm = line.find("--"); if (cm != std::string::npos) line.resize(cm);
        std::istringstream ls(line); std::string w;
        while (ls >> w) {
            if (w == "/") { end = true; break; }
            if (!w.empty() && w.back() == '/') { w.pop_back(); end = true; }
            const size_t st = w.find('*');
            if (st != std::string::npos && st > 0 && std::all_of(w.begin(), w.begin() + long(st), ::isdigit)) {
                const int n = std::atoi(w.substr(0, st).c_str()); for (int k = 0; k < n; ++k) t.push_back("");
                if (st + 1 < w.size()) t.push_back(w.substr(st + 1));
            } else if (w == "1*" || w == "*") t.push_back("");
            else if (!w.empty()) { if (w.front() == '\'') w = w.substr(1); if (!w.empty() && w.back() == '\'') w.pop_back(); t.push_back(w); }
            if (end) break;
        }
    }
    auto num = [&](size_t i, double d) { return i < t.size() && !t[i].empty() ? std::atof(t[i].c_str()) : d; };
    auto str = [&](size_t i, const std::string& d) { std::string v = i < t.size() && !t[i].empty() ? t[i] : d; for (auto& ch : v) ch = char(std::toupper((unsigned char)ch)); return v; };
    GuideRateSpec g; g.set = true;
    g.minInterval = num(0, 0); g.phase = str(1, "NONE");
    g.A = num(2, 0); g.B = num(3, 0); g.C = num(4, 0); g.D = num(5, 0); g.E = num(6, 0); g.F = num(7, 0);
    g.allowIncrease = str(8, "YES") != "NO"; g.damping = num(9, 1.0);
    return g;
}

bool deckModel(const deck::RawDeck& D, const Case& c, DeckModel& out, std::string* error) {
    if (c.comp.enabled) return deckModel(c, out, error);
    Case cc = c;
    if (!compositionalFromDeck(D, cc)) { if (error) *error = "not a compositional deck (COMPS / NCOMPS and EOS data missing)"; return false; }
    return deckModel(cc, out, error);
}

bool deckModel(const Case& c, DeckModel& out, std::string* error) {
    auto fail = [&](const std::string& m) { if (error) *error = m; return false; };
    // C17 review fix: per-cell input arrays must cover every simulator cell (as Simulator::setup does): LGR cells (with their
    // local PORO / PERM / ... values), MINC shells, DPNUM = 0 folding and hydraulic-fracture cells. Before, LGR / MINC cells read
    // out of range or default values (well indices of LGR completions indexed PERMX past its end).
    {
        const GridInput& g0 = c.grid;
        const bool need = (!g0.lgrs.empty() && !g0.lgrExpanded) || (g0.minc > 1 && g0.dual == 1 && !g0.mincExpanded) ||
                          (g0.dual && !g0.dpnum.empty() && !g0.dpnumApplied) || (!g0.hfracs.empty() && !g0.hfExpanded);
        if (need) {
            Case cx = c;
            try { expandLgrArrays(cx); expandMincArrays(cx); applyDpnum(cx); prepareHydraulicFractures(cx); }
            catch (const std::exception& e) { return fail(std::string("grid: ") + e.what()); }
            return deckModel(cx, out, error);
        }
    }
    out = DeckModel();
    const bool field = c.units == UnitSystem::Field;
    out.field = field;
    const UnitsSI U(field);
    if (!c.comp.enabled || c.comp.fluid.eosRegions.empty()) return fail("the case has no compositional fluid (Fluid > Model: Compositional / CO2 storage)");
    const eos::CompFluid& fl = c.comp.fluid;
    const bool co2s = c.comp.co2store;
    CompCase& cc = out.base;
    cc.fluid = fl.eosRegions[0];
    cc.T = fl.resT;
    if (co2s) {
        if (cc.fluid.nc() != 2) return fail("CO2STORE: the fluid must be CO2 + H2O");
        cc.co2store = true; cc.salinity = fl.salinity;
        bool varT = false; for (auto& r : fl.tempvd) varT = varT || std::fabs(r.second - fl.tempvd.front().second) > 0.01;
        if (varT) out.notes.push_back("CO2STORE: temperature varies with depth (RTEMPVD) - the first value is used (isothermal)");
        if (!(cc.T > 273.15 && cc.T < 373.15 + 1e-9)) out.notes.push_back("CO2STORE: the Spycher-Pruess solubility model is calibrated for 12 - 100 C");
    }
    if (!(cc.T > 0)) return fail("reservoir temperature (RTEMP / RTEMPA / TEMPVD) is missing");
    if (!co2s && fl.co2sol) {   // CO2SOL: CO2 is an EOS component that also dissolves in the brine
        bool hasCo2 = false;
        for (const auto& q : cc.fluid.comp) { std::string u = q.name; for (char& ch : u) ch = char(std::toupper(static_cast<unsigned char>(ch))); hasCo2 = hasCo2 || u == "CO2"; }
        if (!c.phases.water) out.notes.push_back("CO2SOL: ignored - the model has no water phase");
        else if (!hasCo2) out.notes.push_back("CO2SOL: ignored - no component named CO2 (CNAMES)");
        else {
            cc.co2sol = true; cc.solSalinity = fl.salinity;
            if (fl.solubili.size() >= 2) {   // SOLUBILI: Rsw (sm3/sm3) -> saturated mole fraction x = r / (1 + r)
                const double molPerSm3Gas = cc.stdP / (eos::Rgas * cc.stdT), molPerSm3W = cc.rhoWsc / 0.018015;
                for (const auto& r : fl.solubili) {
                    const double rr = r.rsw * molPerSm3Gas / molPerSm3W;
                    cc.solTabP.push_back(r.p); cc.solTabX.push_back(rr / (1 + rr));
                }
                out.notes.push_back("CO2SOL: solubility from the SOLUBILI table (" + std::to_string(fl.solubili.size()) + " rows)");
            }
            out.notes.push_back("CO2SOL: CO2 dissolves in the aqueous phase (fugacity equality, Spycher-Pruess / Duan-Sun H(p) at " +
                                std::to_string(fl.salinity).substr(0, 5) + " mol/kg NaCl)");
        }
    }
    if (fl.eosRegions.size() > 1 && !co2s) out.notes.push_back("multiple EOS regions: the reservoir EOS (region 1) is used everywhere");
    cc.stdT = fl.stdT; cc.stdP = fl.stdP;
    // ---- grid (processed connection list, SI)
    Grid g;
    try { g.build(c); } catch (const std::exception& e) { return fail(std::string("grid: ") + e.what()); }
    const int na = g.numActive();
    // grid arrays (PORO, PERM*, NTG, *NUM, DISPERC) were expanded above and are indexed by the global cell directly;
    // arrays that are not expanded (enumerated initial state) map MINC shells, LGR and fracture cells to their host input cell
    std::unordered_map<int, int> hfHost;
    for (const Grid::HfCell& H : g.hfCells) hfHost[H.global] = H.host;
    auto hostQ = [&](int gq) {
        for (int it = 0; it < 3; ++it) {
            if (gq < g.numCells()) return gq;
            if (g.isMincShell(gq)) gq = g.mincHost[size_t(gq - g.numCells() - g.numLgrCells())];
            else if (g.isLocal(gq)) gq = g.lgrHost[size_t(gq - g.numCells())];
            else { auto h = hfHost.find(gq); if (h == hfHost.end()) return gq; gq = h->second; }
        }
        return gq;
    };
    if (g.numMincCells() > 0) out.notes.push_back("MINC: " + std::to_string(g.mincN) + " nested matrix shells per block (" + std::to_string(g.numMincCells()) + " shell cells)");
    cc.cellPV.resize(size_t(na)); cc.cellDepth.resize(size_t(na));
    for (int a = 0; a < na; ++a) { cc.cellPV[size_t(a)] = g.poreVolRef[size_t(a)] * U.pv; cc.cellDepth[size_t(a)] = g.depth[size_t(a)] * U.len; }
    cc.conns.clear();
    const GridInput& gr = c.grid;
    cc.cellPoro.assign(size_t(na), 0.2);
    for (int a = 0; a < na; ++a) { const int gq = g.globalOf[size_t(a)]; if (gq < int(gr.poro.size())) cc.cellPoro[size_t(a)] = gr.poro[size_t(gq)]; }
    // diffusive transmissibility phi A / d = T / k_harmonic * phi_harmonic (direction of the connection; NTG and
    // transmissibility multipliers act on diffusion too)
    const double mD = 9.869233e-16;
    auto permDir = [&](int gq, char d) {
        const std::vector<double>& k = (d == 'Y' && !gr.permy.empty()) ? gr.permy : (d == 'Z' && !gr.permz.empty()) ? gr.permz : gr.permx;
        return gq < int(k.size()) ? g.effPerm(k, gq) : 0.0;
    };
    for (const Connection& cn : g.conns) {
        const double T = cn.trans * U.trans;
        double Td = 0;
        if (cn.dir == 'M' && cn.geo > 0) {
            // C17: matrix-fracture / MINC shell diffusion over the shape-factor geometry with the matrix porosity
            const int ga = g.globalOf[size_t(cn.a)];
            Td = cn.geo * U.len * (g.isFracture(ga) ? cc.cellPoro[size_t(cn.b)] : cc.cellPoro[size_t(cn.a)]);
        } else if (cn.dir != 'M') {
            const int ga = g.globalOf[size_t(cn.a)], gb = g.globalOf[size_t(cn.b)];
            const double ka = permDir(ga, cn.dir), kb = permDir(gb, cn.dir);
            const double pa = cc.cellPoro[size_t(cn.a)], pb = cc.cellPoro[size_t(cn.b)];
            if (cn.geo > 0 && pa > 0 && pb > 0) {   // C17: exact geometry A / d (NTG of horizontal connections kept)
                double nt = 1;
                if (cn.dir == 'X' || cn.dir == 'Y') {
                    const int qa = ga, qb = gb;
                    const double na_ = qa < int(gr.ntg.size()) ? gr.ntg[size_t(qa)] : 1.0, nb_ = qb < int(gr.ntg.size()) ? gr.ntg[size_t(qb)] : 1.0;
                    nt = na_ > 0 && nb_ > 0 ? 2 * na_ * nb_ / (na_ + nb_) : 0.0;
                }
                Td = cn.geo * U.len * nt * (2 * pa * pb / (pa + pb));
            } else if (ka > 0 && kb > 0 && pa > 0 && pb > 0) Td = T / (2 * ka * kb / (ka + kb) * mD) * (2 * pa * pb / (pa + pb));
        }
        double Tdisp = 0;   // C17 DISPERC: harmonic face dispersivity / centre distance
        if (!gr.disperc.empty() && cn.dir != 'M') {
            const int qa = g.globalOf[size_t(cn.a)], qb = g.globalOf[size_t(cn.b)];
            const double aa = qa < int(gr.disperc.size()) ? gr.disperc[size_t(qa)] : 0.0, ab = qb < int(gr.disperc.size()) ? gr.disperc[size_t(qb)] : 0.0;
            const Vec3& ca = g.center[size_t(g.globalOf[size_t(cn.a)])]; const Vec3& cb = g.center[size_t(g.globalOf[size_t(cn.b)])];
            const double ex = double(ca.x) - cb.x, ey = double(ca.y) - cb.y, ez = double(ca.z) - cb.z;
            const double d = std::sqrt(ex * ex + ey * ey + ez * ez);
            if (aa > 0 && ab > 0 && d > 0) Tdisp = 2 * aa * ab / (aa + ab) / d;   // both in deck length units
        }
        cc.conns.push_back({cn.a, cn.b, T, Td, Tdisp});
    }
    out.cellGlobal = g.globalOf;
    // ---- rock, water, saturation functions
    cc.cr = c.props.rock.cr * U.comp; cc.pref = c.props.rock.pref * U.p;
    cc.water = c.phases.water && !co2s;
    if (!co2s && c.comp.waterVap) {
        if (!cc.water) out.notes.push_back("water vaporisation: ignored - the model has no mobile water phase");
        else if (cc.T >= 647.0) out.notes.push_back("water vaporisation: ignored - temperature above the water critical point");
        else { cc.waterVap = true; out.notes.push_back("water vaporisation: water vapour in the gas phase, y_w = Psat(T)/p (IAPWS-IF97)"); }
    }
    cc.rhoWsc = c.props.density.water * U.dens;
    cc.bwRef = c.props.pvtw.bw; cc.cw = c.props.pvtw.cw * U.comp; cc.pwRef = c.props.pvtw.pref * U.p; cc.muw = c.props.pvtw.visc * 1e-3;
    if (c.props.pvtw.cv != 0) out.notes.push_back("PVTW viscosibility is not used by the compositional core");
    if (c.grid.endscale) out.notes.push_back("end-point scaling is not applied by the compositional core");
    std::string satErr;
    // saturation tables of one SATNUM region: SWOF / SGOF (+ SOF3 on the union of oil saturations); CO2STORE: brine / CO2
    auto buildSat = [&](int region, SatTables& S) {
        const SatRegion R = c.props.satRegion(region);
        const auto& swof = R.swof; const auto& sgof = R.sgof;
        S = SatTables();
        if (co2s) {   // OIL + GAS decks (SGOF: oil = brine) or GAS + WATER decks (SWOF / SWFN + SGOF / SGFN / SGWFN)
            if (c.phases.oil) {
                for (size_t i = sgof.size(); i-- > 0;) { S.sw.push_back(1.0 - sgof[i].s); S.krw.push_back(sgof[i].kro); S.pcow.push_back(0.0); }
            } else for (const SatRow& r : swof) { S.sw.push_back(r.s); S.krw.push_back(r.krw); S.pcow.push_back(r.pc * U.p); }
            for (const SatRow& r : sgof) { S.sg.push_back(r.s); S.krg.push_back(r.krw); S.pcog.push_back(r.pc * U.p); }
            if (S.sw.empty() || S.sg.empty()) satErr = "CO2STORE: SGOF (OIL + GAS) or SWOF/SWFN + SGOF/SGFN/SGWFN (GAS + WATER) saturation tables are required";
            return;
        }
        for (const SatRow& r : swof) { S.sw.push_back(r.s); S.krw.push_back(r.krw); S.pcow.push_back(r.pc * U.p); }
        for (const SatRow& r : sgof) { S.sg.push_back(r.s); S.krg.push_back(r.krw); S.pcog.push_back(r.pc * U.p); }
        const double swc = S.sw.empty() ? 0.0 : S.sw.front();
        if (!swof.empty() || !sgof.empty()) {   // SOF3 on the union of the oil saturations of both tables
            std::vector<double> soW, krowV, soG, krogV;
            for (size_t i = swof.size(); i-- > 0;) { soW.push_back(1.0 - swof[i].s); krowV.push_back(swof[i].kro); }
            for (size_t i = sgof.size(); i-- > 0;) { soG.push_back(1.0 - swc - sgof[i].s); krogV.push_back(sgof[i].kro); }
            std::vector<double> so = soW; so.insert(so.end(), soG.begin(), soG.end());
            std::sort(so.begin(), so.end());
            std::vector<double> u;
            for (double v : so) if (u.empty() || v - u.back() > 1e-9) u.push_back(std::max(0.0, v));
            for (double v : u) {
                S.so.push_back(v);
                S.krow.push_back(soW.empty() ? (soG.empty() ? 0.0 : lin(soG, krogV, v)) : lin(soW, krowV, v));
                S.krog.push_back(soG.empty() ? lin(soW, krowV, v) : lin(soG, krogV, v));
            }
        }
        if (!cc.water) { S.sw.clear(); S.krw.clear(); S.pcow.clear(); }
    };
    SatTables& S = cc.sat;
    buildSat(0, S);
    if (!satErr.empty()) return fail(satErr);
    cc.swc = co2s ? S.sw.front() : (c.props.swof.empty() ? 0.0 : c.props.swof.front().s);
    const int nsr = c.props.numSatRegions();
    for (int r = 1; r < nsr; ++r) { SatTables T; buildSat(r, T); if (!satErr.empty()) return fail(satErr + " (SATNUM region " + std::to_string(r + 1) + ")"); cc.satRegions.push_back(T); }
    auto regionArray = [&](const std::vector<int>& g1) {   // 1-based global array -> 0-based per active cell
        std::vector<int> v;
        if (g1.empty()) return v;
        v.resize(size_t(na));
        for (int a = 0; a < na; ++a) { const int gq = g.globalOf[size_t(a)]; v[size_t(a)] = gq < int(g1.size()) ? std::clamp(g1[size_t(gq)] - 1, 0, nsr - 1) : 0; }
        return v;
    };
    if (nsr > 1) cc.satnum = regionArray(c.grid.satnum);
    // residual trapping: gas relative-permeability hysteresis (SATOPTS HYSTER + EHYSTR; IMBNUM imbibition tables)
    if (c.props.hyst.on) {
        cc.hyst.on = true;
        cc.hyst.model = c.props.hyst.model >= 2 ? 1 : 0;   // EHYSTR 0/1 Carlson, 2/3/4 Killough
        if (!c.grid.imbnum.empty()) cc.imbnum = regionArray(c.grid.imbnum);
        if (cc.imbnum.empty() && c.comp.trap.sgtMax <= 0) out.notes.push_back("hysteresis: no IMBNUM imbibition tables and no Land maximum trapped saturation - drainage curves only");
        if (c.props.hyst.model == 1 || c.props.hyst.model == 3) out.notes.push_back("hysteresis: wetting-phase hysteresis is not modelled (drainage curves for the wetting phase)");
    }
    cc.kroModel = c.props.kroModel == 2 ? 1 : 0;
    if (c.props.kroModel == 1) out.notes.push_back("STONE1: the compositional core uses the ECLIPSE default three-phase oil kr");
    if (!c.aquifers.empty()) out.notes.push_back("analytic aquifers are not modelled by the compositional core");
    // ---- separators (FIELDSEP + SEPCOND trains) and injection streams (WELLSTRE + WINJGAS per report step), from c.comp
    cc.separator = c.comp.separator;
    for (const auto& tr : c.comp.sepTrains) cc.sepTrains.push_back(tr.stages);
    auto trainOf = [&](const std::string& well) {
        auto it = c.comp.wellSep.find(well);
        auto match = [&](const std::string& pat) {   // exact (case-insensitive) or trailing-'*' pattern
            const std::string P = up(pat), W = up(well);
            if (!P.empty() && P.back() == '*') return W.compare(0, P.size() - 1, P, 0, P.size() - 1) == 0;
            return P == W;
        };
        if (it == c.comp.wellSep.end()) for (auto& kv : c.comp.wellSep) if (match(kv.first)) { it = c.comp.wellSep.find(kv.first); break; }
        if (it == c.comp.wellSep.end()) return 0;
        for (size_t k = 0; k < c.comp.sepTrains.size(); ++k) if (up(c.comp.sepTrains[k].name) == up(it->second)) return int(k) + 1;
        return 0;
    };
    if (co2s) cc.trap = c.comp.trap;
    else {
        cc.diffOil = c.comp.diffOil; cc.diffGas = c.comp.diffGas;
        cc.trap.sgtMax = c.comp.trap.sgtMax; cc.trap.tortuosity = c.comp.trap.tortuosity;
        if (c.comp.trap.saltPrecip || !c.comp.trap.minerals.empty()) out.notes.push_back("salt precipitation and mineral trapping are CO2STORE options (ignored in EOS mode)");
    }
    std::map<std::string, std::vector<double>> streams;
    for (const CompStream& st : c.comp.streams) { std::vector<double> z = st.z; z.resize(size_t(cc.fluid.nc()), 0.0); streams[up(st.name)] = z; }
    std::vector<std::map<std::string, std::vector<double>>> injComp(c.schedule.size());
    for (size_t s = 0; s < c.schedule.size(); ++s)
        for (const auto& ws : c.schedule[s].injStreams) {
            auto it = streams.find(up(ws.second));
            if (it == streams.end()) out.notes.push_back("WINJGAS " + ws.first + ": unknown stream " + ws.second);
            else injComp[s][up(ws.first)] = it->second;
        }
    // ---- initial state
    if (co2s) { if (!equilibrateCo2(c, cc, error)) return false; }
    else if (!equilibrate(c, fl, cc, error)) return false;
    if (c.init.active() && !co2s) {
        for (int a = 0; a < na; ++a) {
            const int gq = hostQ(g.globalOf[size_t(a)]);
            const size_t qq = size_t(std::min<int>(gq, int(c.init.pressure.size()) - 1));
            cc.p0[size_t(a)] = c.init.pressure[qq] * U.p;
            cc.sw0[size_t(a)] = c.init.swat.empty() ? cc.swc : c.init.swat[std::min(qq, c.init.swat.size() - 1)];
            cc.z0[size_t(a)] = zAt(fl, cc.cellDepth[size_t(a)]);
        }
        out.notes.push_back("enumerated initial state: PRESSURE / SWAT with ZMFVD / ZI compositions");
    }
    // ---- wells: connections (Peaceman, as the black-oil engine) and controls per report step
    const GridInput& gi = c.grid;
    const size_t nw = c.wells.size();
    for (const auto& w : c.wells) out.wellNames.push_back(w.name);
    std::vector<WellControl> ctrl(nw);
    std::vector<bool> haveCtrl(nw, false);
    std::vector<std::vector<double>> curInj(nw);
    double vgsc = eos::Rgas * cc.stdT / cc.stdP, vwsc = MWW / cc.rhoWsc;
    if (co2s) {   // CO2: real-gas surface volume; brine: surface volume per mol H2O (incl. salt)
        const co2::Fluid cf = co2::makeFluid(cc.T, cc.salinity);
        vgsc = co2::gasSurfaceVolume(cc.stdT, cc.stdP); vwsc = co2::brineSurfaceVolume(cf, cc.stdT, cc.stdP);
    }
    bool notedGroups = false, notedLimits = false;
    std::map<std::string, GroupProdControl> gprodState;
    GuideRateSpec guideState;
    std::map<std::string, GroupInjControl> ginjState[2];
    std::map<std::string, std::string> parentOf;
    for (const auto& pc : c.groupTree) parentOf[up(pc.first)] = up(pc.second);
    auto inGroup = [&](size_t w, const std::string& grp) {   // well w in the subtree of group grp (FIELD = all)
        if (grp == "FIELD") return true;
        std::string g = up(c.wells[w].group);
        for (int d = 0; d < 64 && !g.empty(); ++d) { if (g == grp) return true; auto it = parentOf.find(g); if (it == parentOf.end()) break; g = it->second; }
        return false;
    };
    for (size_t s = 0; s < c.schedule.size(); ++s) {
        const ReportStep& rs = c.schedule[s];
        out.stepDays.push_back(rs.length);
        for (const WellControl& wc : rs.controls)
            for (size_t w = 0; w < nw; ++w) if (c.wells[w].name == wc.well) { ctrl[w] = wc; haveCtrl[w] = true; }
        if (s < injComp.size()) for (auto& kv : injComp[s]) for (size_t w = 0; w < nw; ++w) if (up(c.wells[w].name) == kv.first) curInj[w] = kv.second;
        for (const GroupProdControl& gp : rs.gprod) gprodState[up(gp.group)] = gp;
        for (const GroupInjControl& gq : rs.ginj) ginjState[gq.phase == 2 ? 1 : 0][up(gq.group)] = gq;
        std::vector<CompWell> ws;
        for (size_t w = 0; w < nw; ++w) {
            const WellInput& in = c.wells[w];
            CompWell cw; cw.name = in.name;
            cw.refDepth = in.refDepth > 0 ? in.refDepth * U.len : -1;
            for (const Completion& cp : in.completions) {
                if (int(s) < std::max(0, cp.fromStep) || (cp.untilStep >= 0 && int(s) >= cp.untilStep) || !cp.open) continue;
                const int lg = cp.lgr;
                for (int k = std::min(cp.k1, cp.k2); k <= std::max(cp.k1, cp.k2); ++k) {
                    const int kk = (gi.dual == 1 && k < gi.nzMatrix()) ? k + gi.nzMatrix() : k;
                    std::vector<int> targets{lg >= 0 ? g.localCell(lg, cp.i, cp.j, k) : gi.index(cp.i, cp.j, kk)};
                    if (lg < 0 && g.isHost(targets[0])) targets = g.hostWellColumn(targets[0], cp.dir);
                    for (int q : targets) {
                        const int a = q < int(g.activeOf.size()) ? g.activeOf[size_t(q)] : -1;
                        if (a < 0) continue;
                        double wi = cp.cf;
                        if (wi <= 0) {
                            const double kx = q < int(gi.permx.size()) ? g.effPerm(gi.permx, q) : 0.0;
                            const double ky = q < int(gi.permy.size()) && !std::isnan(gi.permy[size_t(q)]) ? g.effPerm(gi.permy, q) : kx,
                                         kz = q < int(gi.permz.size()) && !std::isnan(gi.permz[size_t(q)]) ? g.effPerm(gi.permz, q) : kx;
                            const double DX = g.cellDx[size_t(q)], DY = g.cellDy[size_t(q)], DZ = g.cellDz[size_t(q)];
                            double k1, k2, d1, d2, hh;
                            if (cp.dir == 'X' || cp.dir == 'x') { k1 = ky; k2 = kz; d1 = DY; d2 = DZ; hh = DX; }
                            else if (cp.dir == 'Y' || cp.dir == 'y') { k1 = kx; k2 = kz; d1 = DX; d2 = DZ; hh = DY; }
                            else { k1 = kx; k2 = ky; d1 = DX; d2 = DY; hh = DZ * (q < int(gi.ntg.size()) ? gi.ntg[size_t(q)] : 1.0); }
                            const double kh = cp.kh > 0 ? cp.kh : std::sqrt(k1 * k2) * hh;
                            const double r12 = std::sqrt(k2 / k1), r21 = std::sqrt(k1 / k2);
                            const double ro = 0.28 * std::sqrt(r12 * d1 * d1 + r21 * d2 * d2) / (std::pow(k2 / k1, 0.25) + std::pow(k1 / k2, 0.25));
                            wi = Units::of(c.units).cDarcy * 2.0 * M_PI * kh / (std::log(ro / (0.5 * cp.diameter)) + cp.skin);
                        }
                        wi *= in.piMult * cp.mult;
                        cw.cells.push_back(a); cw.wi.push_back(wi * U.trans);
                        for (const Grid::HfCell& H : g.hfCells) {   // C15: the fractures of the host cell (EDFM well index)
                            if (H.host != q) continue;
                            const int fa = g.activeOf[size_t(H.global)];
                            if (fa < 0) continue;
                            const double rw = 0.5 * cp.diameter, re = 0.14 * std::sqrt(H.len * H.len + H.h * H.h);
                            const double wf = Units::of(c.units).cDarcy * 2 * M_PI * H.kf * H.w / std::max(std::log(std::max(re, 1.01 * rw) / rw), 1.0) * in.piMult * cp.mult;
                            cw.cells.push_back(fa); cw.wi.push_back(wf * U.trans);
                        }
                    }
                }
            }
            const WellControl& wc = ctrl[w];
            cw.sep = trainOf(in.name);
            cw.open = haveCtrl[w] && wc.status == WellStatus::Open && !cw.cells.empty();
            const bool inj = in.type != WellType::Producer;
            cw.injector = inj;
            if (inj) {
                const int ph = wc.phase ? wc.phase : (in.type == WellType::WaterInjector ? 1 : 2);
                cw.water = ph == 1;
                if (co2s) {   // CO2 (gas) or brine injection; both are components of the single flash
                    cw.zinj = ph == 1 ? std::vector<double>{0.0, 1.0} : std::vector<double>{1.0, 0.0};
                    cw.bhpMax = wc.bhp > 0 ? wc.bhp * U.p : 0.0;
                    if (wc.mode == Control::BHP) { cw.injBhpControl = true; cw.bhpMax = wc.bhp * U.p; }
                    else cw.rate = std::max(0.0, wc.rate) * (ph == 1 ? U.liq / vwsc : U.gas / vgsc) / 86400.0;
                    cw.water = false;
                    ws.push_back(cw);
                    continue;
                }
                cw.bhpMax = wc.bhp > 0 ? wc.bhp * U.p : 0.0;
                if (wc.mode == Control::BHP) { cw.injBhpControl = true; cw.bhpMax = wc.bhp * U.p; }
                else cw.rate = std::max(0.0, wc.rate) * (cw.water ? U.liq / vwsc : U.gas / vgsc) / 86400.0;
                if (!cw.water) {
                    cw.zinj = curInj[w];
                    if (cw.zinj.empty() && cw.open) {
                        cw.zinj.assign(size_t(cc.fluid.nc()), 0.0); cw.zinj[0] = 1.0;
                        out.notes.push_back("gas injector " + in.name + " without WINJGAS stream: pure " + (fl.names.empty() ? std::string("component 1") : fl.names[0]) + " injected");
                    }
                }
                if (cw.water && !cc.water) cw.open = false;
                if (wc.mode == Control::RESV && !notedLimits) { out.notes.push_back("injector RESV control: surface rate used"); notedLimits = true; }
            } else {
                cw.bhp = (wc.bhp > 0 ? wc.bhp : (field ? 14.696 : 1.01325)) * U.p;
                if (co2s) {   // the brine is the liquid: ORAT (OIL + GAS decks) / WRAT / LRAT all limit the brine rate
                    switch (wc.mode) {
                        case Control::ORAT: cw.control = 5; cw.target = wc.orat * U.liq / 86400.0; break;
                        case Control::WRAT: cw.control = 5; cw.target = wc.wrat * U.liq / 86400.0; break;
                        case Control::LRAT: cw.control = 5; cw.target = wc.lrat * U.liq / 86400.0; break;
                        case Control::GRAT: cw.control = 4; cw.target = wc.grat * U.gas / 86400.0; break;
                        case Control::RESV: cw.control = 2; cw.target = wc.resv * U.liq / 86400.0; break;
                        default: cw.control = 0; break;
                    }
                } else
                switch (wc.mode) {
                    case Control::ORAT: cw.control = 3; cw.target = wc.orat * U.liq / 86400.0; break;
                    case Control::GRAT: cw.control = 4; cw.target = wc.grat * U.gas / 86400.0; break;
                    case Control::WRAT: cw.control = 5; cw.target = wc.wrat * U.liq / 86400.0; break;
                    case Control::LRAT: cw.control = 6; cw.target = wc.lrat * U.liq / 86400.0; break;
                    case Control::RESV: cw.control = 2; cw.target = wc.resv * U.liq / 86400.0; break;
                    default: cw.control = 0; break;
                }
                if (cw.control != 0 && cw.target < 0) cw.control = 0;
                // secondary rate limits (WCONPROD items 4-8 other than the control mode; not for WCONHIST observed rates)
                if (!wc.history) {
                    auto lim = [&](int ctl, double v, double f) { if (v >= 0 && ctl != cw.control) cw.limits.push_back({ctl, v * f / 86400.0}); };
                    if (co2s) { lim(4, wc.grat, U.gas); lim(2, wc.resv, U.liq); lim(5, c.phases.oil ? wc.orat : wc.wrat, U.liq); lim(5, wc.lrat, U.liq); }
                    else {
                        lim(3, wc.orat, U.liq); lim(4, wc.grat, U.gas); lim(5, wc.wrat, U.liq); lim(6, wc.lrat, U.liq); lim(2, wc.resv, U.liq);
                    }
                }
            }
            ws.push_back(cw);
        }
        out.stepWells.push_back(ws);
        // group targets in force during this step (GCONPROD / GCONINJE RATE) for the wells under group control (GRUP)
        std::vector<CompGroup> gs;
        if (s < c.schedule.size() && c.schedule[s].guiderat.set) guideState = c.schedule[s].guiderat;
        for (const ExtraKeyword& e : c.extraKeywords)   // GUIDERAT kept as a deck keyword (edited in the group editor) wins
            if (e.name == "GUIDERAT" && e.section == DeckSection::Schedule && e.step == int(s)) guideState = parseGuideRat(e.body);
        for (const auto& kv : gprodState) {
            const GroupProdControl& gp = kv.second;
            int ctl = 0; double T = -1;
            switch (gp.mode) {
                case Control::ORAT: ctl = co2s ? 5 : 3; T = gp.orat * U.liq; break;
                case Control::GRAT: ctl = 4; T = gp.grat * U.gas; break;
                case Control::WRAT: ctl = 5; T = gp.wrat * U.liq; break;
                case Control::LRAT: ctl = co2s ? 5 : 6; T = gp.lrat * U.liq; break;
                case Control::RESV: ctl = 2; T = gp.resv * U.liq; break;
                default: break;
            }
            if (ctl == 0 || T < 0) continue;
            CompGroup G; G.name = kv.first; G.control = ctl; G.target = T / 86400.0;
            {   // GUIDERAT formula (deck units: rates per day, GOR Mscf/stb or sm3/sm3)
                const std::string& ph = guideState.phase;
                G.grPhase = ph == "OIL" ? 0 : ph == "LIQ" ? 1 : ph == "GAS" ? 2 : ph == "RES" ? 3 : -1;
                G.grA = guideState.A; G.grB = guideState.B; G.grC = guideState.C; G.grD = guideState.D; G.grE = guideState.E; G.grF = guideState.F;
                const double liq = 86400.0 / U.liq, gas = 86400.0 / U.gas, gor = U.liq / U.gas;   // SI -> deck
                G.potScale = G.grPhase == 2 ? gas : liq;
                if (G.grPhase == 2) { G.r1Scale = 1.0 / gor; G.r2Scale = 1.0 / gor; }   // WGR, OGR: liquid per gas
                else { G.r1Scale = 1.0; G.r2Scale = gor; }                              // WOR / WCT dimensionless; GOR / GLR gas per liquid
            }
            for (size_t w = 0; w < nw; ++w)
                if (ws[w].open && !ws[w].injector && haveCtrl[w] && ctrl[w].deckMode == "GRUP" && inGroup(w, kv.first)) G.wells.push_back(int(w));
            if (!G.wells.empty()) gs.push_back(G);
        }
        for (int ph = 0; ph < 2; ++ph)
            for (const auto& kv : ginjState[ph]) {
                const GroupInjControl& gq = kv.second;
                if (gq.rate < 0) { if (gq.mode != "NONE" && !notedGroups) { out.notes.push_back("group injection " + gq.mode + " targets are not applied by the compositional core (RATE only)"); notedGroups = true; } continue; }
                CompGroup G; G.name = kv.first; G.injector = true;
                G.target = ph == 1 ? gq.rate * U.gas / 86400.0 / vgsc : gq.rate * U.liq / 86400.0 / vwsc;
                for (size_t w = 0; w < nw; ++w) {
                    const bool gasInj = c.wells[w].type == WellType::GasInjector, watInj = c.wells[w].type == WellType::WaterInjector;
                    if (!ws[w].open || !ws[w].injector || !haveCtrl[w] || ctrl[w].deckMode != "GRUP" || !inGroup(w, kv.first)) continue;
                    if ((ph == 1 && gasInj) || (ph == 0 && watInj)) G.wells.push_back(int(w));
                }
                if (!G.wells.empty()) gs.push_back(G);
            }
        out.stepGroups.push_back(gs);
    }
    return true;
}

bool deckModelFromFile(const std::string& path, DeckModel& out, std::string* error) {
    try {
        const deck::RawDeck D = deck::parseDeckFile(path);
        const Case c = importEclipseDeck(path);
        return deckModel(D, c, out, error);
    } catch (const std::exception& e) {
        if (error) *error = e.what();
        return false;
    }
}

// ------------------------------------------------------------------------------------------------- run
bool runDeckModel(const DeckModel& m, std::vector<CompSummary>& out, const std::function<bool(const CompSummary&, CompSimulator&)>& cb,
                  std::string* error, int maxSteps, CompStats* stats, const DeckRestart* rst) {
    const UnitsSI U(m.field);
    CompSimulator sim;
    try { sim.setup(m.base); } catch (const std::exception& e) { if (error) *error = e.what(); return false; }
    out.clear();
    const size_t nw = m.wellNames.size();
    double t = 0;
    std::vector<CompWellReport> prev(nw);
    CompSummary last;
    auto field = [&](CompSummary& S) {
        double o, g, w; sim.inPlaceSurface(o, g, w);
        S.foip = o / U.liq; S.fgip = g / U.gas; S.fwip = w / U.liq; S.fpr = sim.hcpvPressure() / U.p;
        if (m.base.co2store) {
            const Co2Inventory I = sim.co2Trapping();
            S.co2Mob = I.mobile + I.structural; S.co2Immob = I.residual; S.co2Diss = I.dissolved;
            S.co2Struct = I.structural; S.co2Min = I.mineral; S.co2KrImmob = I.immobileKr; S.saltSolid = I.saltSolid;
        } else if (m.base.co2sol) {   // CO2SOL: CO2 in the hydrocarbon phases ("mobile") and dissolved in the brine (kg)
            double hc = 0, aq = 0; sim.co2SolSplit(hc, aq);
            S.co2Mob = hc; S.co2Diss = aq;
        }
    };
    // runner bookkeeping saved with the engine state: time and the cumulative field totals of the last row
    auto rowIO = [](auto& ar, CompSummary& S) {
        ar(S.time, S.fopr, S.fgpr, S.fwpr, S.fgir, S.fwir, S.fopt, S.fgpt, S.fwpt, S.fgit, S.fwit, S.fpr, S.foip, S.fgip, S.fwip);
        ar(S.co2Mob, S.co2Immob, S.co2Diss, S.co2Struct, S.co2Min, S.co2KrImmob, S.saltSolid);
    };
    size_t s0 = 0;
    if (rst && rst->startStep > 0) {
        try {
            StateIn in(rst->state);
            if (in.r.str() != "RESSIM-COMP-RUN-1") throw std::runtime_error("not a compositional run restart");
            CompSummary S; rowIO(in, S);
            const size_t sr = size_t(rst->startStep);
            if (sr >= 1 && sr - 1 < m.stepWells.size()) sim.wells() = m.stepWells[sr - 1];   // well set of the restart step
            if (sr >= 1 && sr - 1 < m.stepGroups.size()) sim.groups() = m.stepGroups[sr - 1];
            sim.loadState(in.r.str());
            S.wells.assign(nw, CompSummary::Well());
            t = S.time; s0 = size_t(rst->startStep);
            if (s0 > m.stepDays.size()) throw std::runtime_error("restart step beyond the schedule");
            const auto& R0 = sim.wellReports();
            for (size_t w = 0; w < nw && w < R0.size(); ++w) prev[w] = R0[w];
            out.push_back(S);
        } catch (const std::exception& e) { if (error) *error = std::string("restart: ") + e.what(); return false; }
    } else {   // initial row
        CompSummary S; S.wells.resize(nw); field(S); out.push_back(S);
        if (cb && !cb(S, sim)) return true;
    }
    const size_t ns = maxSteps >= 0 ? std::min(m.stepDays.size(), size_t(maxSteps)) : m.stepDays.size();
    for (size_t s = s0; s < ns; ++s) {
        sim.wells() = m.stepWells[s];
        sim.groups() = s < m.stepGroups.size() ? m.stepGroups[s] : std::vector<CompGroup>();
        const double len = m.stepDays[s];
        t += len;
        if (!sim.advanceTo(t * 86400.0, std::min(len, 30.0) * 86400.0)) {
            if (error) *error = "time step collapsed at day " + std::to_string(sim.time() / 86400.0);
            if (stats) *stats = sim.stats();
            return false;
        }
        CompSummary S = out.back(); S.time = t; S.wells.assign(nw, CompSummary::Well());
        S.fopr = S.fgpr = S.fwpr = S.fgir = S.fwir = 0;
        const auto& R = sim.wellReports();
        for (size_t w = 0; w < nw && w < R.size(); ++w) {
            const double dO = (R[w].cumO - prev[w].cumO) / U.liq, dG = (R[w].cumG - prev[w].cumG) / U.gas, dW = (R[w].cumW - prev[w].cumW) / U.liq;
            const double dGI = (R[w].cumGI - prev[w].cumGI) / U.gas, dWI = (R[w].cumWI - prev[w].cumWI) / U.liq;
            CompSummary::Well& W = S.wells[w];
            W.oil = dO / len; W.gas = dG / len; W.water = dW / len; W.ginj = dGI / len; W.winj = dWI / len;
            W.bhp = R[w].bhp / U.p; W.open = m.stepWells[s][w].open; W.onBhp = R[w].onBhp;
            W.molar.assign(R[w].cumMolP.size(), 0.0);
            for (size_t i = 0; i < R[w].cumMolP.size(); ++i) {
                const double p0 = i < prev[w].cumMolP.size() ? prev[w].cumMolP[i] : 0.0, i0 = i < prev[w].cumMolI.size() ? prev[w].cumMolI[i] : 0.0;
                W.molar[i] = ((R[w].cumMolP[i] - p0) - (R[w].cumMolI[i] - i0)) / (len * 86400.0);
            }
            S.fopr += W.oil; S.fgpr += W.gas; S.fwpr += W.water; S.fgir += W.ginj; S.fwir += W.winj;
            S.fopt += dO; S.fgpt += dG; S.fwpt += dW; S.fgit += dGI; S.fwit += dWI;
            prev[w] = R[w];
        }
        field(S);
        out.push_back(S);
        if (rst && rst->wants(int(s + 1), s + 1 == ns)) {
            StateOut o; o.w.str("RESSIM-COMP-RUN-1"); rowIO(o, S); o.w.str(sim.saveState());
            rst->onSave(int(s + 1), t, o.w.buf);
        }
        if (cb && !cb(S, sim)) break;
    }
    if (stats) *stats = sim.stats();
    return true;
}

}  // namespace comp
}  // namespace ressim
