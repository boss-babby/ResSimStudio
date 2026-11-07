// Built-in example library (see examples.h). All geology is generated procedurally with a portable RNG so that
// every platform builds bit-identical cases.
#include "ressim/gasdata.h"
#include "ressim/examples.h"
#include "ressim/geomech.h"
#include "ressim/tools.h"
#include "ressim/simulator.h"
#include "ressim/steam.h"
#include "ressim/io.h"
#include "ressim/hydrofrac.h"
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <algorithm>
#include <array>

namespace ressim {

namespace {

const double kPi = 3.14159265358979323846;

// ------------------------------------------------------------------ portable random numbers (splitmix64)
struct Rng {
    uint64_t s;
    explicit Rng(uint64_t seed) : s(seed * 0x9E3779B97F4A7C15ull + 0x1234567ull) {}
    uint64_t next() {
        uint64_t z = (s += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }
    double uni() { return double(next() >> 11) * (1.0 / 9007199254740992.0); }   // [0,1)
    double range(double a, double b) { return a + (b - a) * uni(); }
};

// Smooth, zero-mean, unit-variance random field on (x, y) built from M random plane waves.
struct SmoothField {
    struct Wave { double kx, ky, ph; };
    std::vector<Wave> w;
    SmoothField(uint64_t seed, double corrLen, int M = 24) {
        Rng r(seed);
        for (int m = 0; m < M; ++m) {
            const double ang = r.range(0, kPi), kmag = 2 * kPi / corrLen * r.range(0.4, 1.6);
            w.push_back({kmag * std::cos(ang), kmag * std::sin(ang), r.range(0, 2 * kPi)});
        }
    }
    double operator()(double x, double y) const {
        double s = 0;
        for (const Wave& v : w) s += std::cos(v.kx * x + v.ky * y + v.ph);
        return s * std::sqrt(2.0 / double(w.size()));
    }
};

// ------------------------------------------------------------------ fluids
// SPE1 (Odeh 1981) black-oil PVT, in FIELD units or converted to METRIC.
void spe1Fluids(Case& c, bool metric) {
    PropsInput& p = c.props;
    p.pvtw = {4017.55, 1.038, 3.22e-6, 0.318, 0.0};
    p.rock = {14.7, 3e-6};
    p.density = {53.66, 64.49, 0.0533};
    p.pvto = {
        {0.001, 14.7, 1.062, 1.04, {}},  {0.0905, 264.7, 1.15, 0.975, {}}, {0.18, 514.7, 1.207, 0.91, {}},
        {0.371, 1014.7, 1.295, 0.83, {}}, {0.636, 2014.7, 1.435, 0.695, {}}, {0.775, 2514.7, 1.5, 0.641, {}},
        {0.93, 3014.7, 1.565, 0.594, {}}, {1.27, 4014.7, 1.695, 0.51, {{9014.7, 1.579, 0.74}}},
        {1.618, 5014.7, 1.827, 0.449, {{9014.7, 1.737, 0.631}}}};
    p.pvdg = {{14.7, 166.666, 0.008}, {264.7, 12.093, 0.0096}, {514.7, 6.274, 0.0112}, {1014.7, 3.197, 0.014},
              {2014.7, 1.614, 0.0189}, {2514.7, 1.294, 0.0208}, {3014.7, 1.08, 0.0228}, {4014.7, 0.811, 0.0268},
              {5014.7, 0.649, 0.0309}, {9014.7, 0.386, 0.047}};
    if (metric) {
        const double bar = 0.0689475729, rs = 178.1076, bg = 0.158987295 / 28.3168466, dens = 16.0184634;
        p.pvtw.pref *= bar; p.pvtw.cw /= bar; p.rock.pref *= bar; p.rock.cr /= bar;
        p.density.oil *= dens; p.density.water *= dens; p.density.gas *= dens;
        for (auto& r : p.pvto) { r.rs *= rs; r.pb *= bar; for (auto& u : r.under) u.p *= bar; }
        for (auto& r : p.pvdg) { r.p *= bar; r.B *= bg; }
    }
}

// Dead oil + water (two-phase cases). Oil viscosity muo at the reference pressure.
void deadOilFluids(Case& c, double muo, bool metric) {
    PropsInput& p = c.props;
    c.phases = Phases{true, true, false, false};
    p.pvtw = {3000, 1.01, 3e-6, 0.5, 0.0};
    p.rock = {3000, 4e-6};
    p.density = {52.0, 63.0, 0.06};
    p.pvdo = {{14.7, 1.040, muo * 0.92}, {1014.7, 1.035, muo * 0.95}, {3014.7, 1.025, muo},
              {5014.7, 1.015, muo * 1.05}, {9014.7, 0.995, muo * 1.15}};
    if (metric) {
        const double bar = 0.0689475729, dens = 16.0184634;
        p.pvtw.pref *= bar; p.pvtw.cw /= bar; p.rock.pref *= bar; p.rock.cr /= bar;
        p.density.oil *= dens; p.density.water *= dens; p.density.gas *= dens;
        for (auto& r : p.pvdo) r.p *= bar;
    }
}

// SPE1 saturation tables (SWOF / SGOF) as published with the case.
void spe1Scal(Case& c) {
    c.props.swof = {{0.12, 0, 1, 0}, {0.18, 4.64876e-8, 1, 0}, {0.24, 1.86e-7, 0.997, 0}, {0.3, 4.18388e-7, 0.98, 0},
                    {0.36, 7.43802e-7, 0.7, 0}, {0.42, 1.16219e-6, 0.35, 0}, {0.48, 1.67355e-6, 0.2, 0},
                    {0.54, 2.27789e-6, 0.09, 0}, {0.6, 2.97521e-6, 0.021, 0}, {0.66, 3.7655e-6, 0.01, 0},
                    {0.72, 4.64876e-6, 0.001, 0}, {0.78, 5.625e-6, 0.0001, 0}, {0.84, 6.69421e-6, 0, 0},
                    {0.91, 8.05914e-6, 0, 0}, {1, 1e-5, 0, 0}};
    c.props.sgof = {{0, 0, 1, 0}, {0.001, 0, 1, 0}, {0.02, 0, 0.997, 0}, {0.05, 0.005, 0.98, 0}, {0.12, 0.025, 0.7, 0},
                    {0.2, 0.075, 0.35, 0}, {0.25, 0.125, 0.2, 0}, {0.3, 0.19, 0.09, 0}, {0.4, 0.41, 0.021, 0},
                    {0.45, 0.6, 0.01, 0}, {0.5, 0.72, 0.001, 0}, {0.6, 0.87, 0.0001, 0}, {0.7, 0.94, 0, 0},
                    {0.85, 0.98, 0, 0}, {0.88, 0.984, 0, 0}};
}

// Generic water-wet sandstone curves for three-phase cases with water movement (water flood / aquifer).
// (The SPE1 tables have krw <= 1e-5 - water is practically immobile there - so they suit only gas injection.)
void sandstoneScal(Case& c, double swc = 0.18) {
    coreySWOF(c.props.swof, swc, 0.22, 2.5, 2.0, 0.4, 1.0, 20);
    coreySGOF(c.props.sgof, swc, 0.04, 0.15, 2.0, 2.5, 0.8, 20);
}

// ------------------------------------------------------------------ grid helpers
void boxGrid(Case& c, int nx, int ny, int nz, double dx, double dy, const std::vector<double>& dzLayer, double top) {
    GridInput& g = c.grid;
    g.nx = nx; g.ny = ny; g.nz = nz;
    const int n = g.size();
    g.dx.assign(n, dx); g.dy.assign(n, dy); g.dz.resize(n);
    for (int k = 0; k < nz; ++k) for (int q = 0; q < nx * ny; ++q) g.dz[size_t(k) * nx * ny + q] = dzLayer[size_t(k)];
    g.tops.assign(size_t(nx) * ny, top);
    g.poro.assign(n, 0.2); g.permx.assign(n, 100); g.permy.assign(n, 100); g.permz.assign(n, 10);
}

void setLayerProps(Case& c, const std::vector<double>& kLayer, const std::vector<double>& phiLayer, double kvkh) {
    GridInput& g = c.grid;
    for (int k = 0; k < g.nz; ++k)
        for (int q = 0; q < g.nx * g.ny; ++q) {
            const size_t a = size_t(k) * g.nx * g.ny + q;
            g.permx[a] = g.permy[a] = kLayer[size_t(k)];
            g.permz[a] = kvkh * kLayer[size_t(k)];
            g.poro[a] = phiLayer[size_t(k)];
        }
}

// Log-normal permeability per layer from a smooth field; porosity correlated with log-permeability.
void heterogeneous(Case& c, uint64_t seed, double cellX, double cellY, double corrLen, const std::vector<double>& kMean,
                   double sigmaLn, double phiMean, double phiSpread, double kvkh) {
    GridInput& g = c.grid;
    for (int k = 0; k < g.nz; ++k) {
        SmoothField F(seed * 131 + uint64_t(k) * 7 + 1, corrLen);
        for (int j = 0; j < g.ny; ++j)
            for (int i = 0; i < g.nx; ++i) {
                const double f = F((i + 0.5) * cellX, (j + 0.5) * cellY);
                const size_t a = size_t(g.index(i, j, k));
                const double kk = kMean[size_t(k)] * std::exp(sigmaLn * f - 0.5 * sigmaLn * sigmaLn);
                g.permx[a] = g.permy[a] = kk;
                g.permz[a] = kvkh * kk;
                g.poro[a] = std::max(0.03, std::min(0.35, phiMean + phiSpread * f));
            }
    }
}

// ------------------------------------------------------------------ wells / schedule helpers
WellInput well(const std::string& name, WellType t, int i, int j, int k1, int k2, double diameter = 0.5) {
    WellInput w; w.name = name; w.type = t; w.i = i; w.j = j;
    w.group = t == WellType::Producer ? "PROD" : "INJ";
    Completion cp; cp.i = i; cp.j = j; cp.k1 = k1; cp.k2 = k2; cp.diameter = diameter;
    w.completions.push_back(cp);
    return w;
}
WellControl prodCtl(const std::string& w, Control m, double target, double bhp) {
    WellControl c; c.well = w; c.mode = m; c.bhp = bhp;
    switch (m) {
        case Control::ORAT: c.orat = target; break; case Control::WRAT: c.wrat = target; break;
        case Control::GRAT: c.grat = target; break; case Control::LRAT: c.lrat = target; break;
        case Control::RESV: c.resv = target; break; case Control::BHP: c.bhp = target; break;
        default: break;
    }
    return c;
}
WellControl injCtl(const std::string& w, double rate, double bhp) {
    WellControl c; c.well = w; c.mode = Control::RATE; c.rate = rate; c.bhp = bhp; return c;
}
void steps(Case& c, int count, double len) { for (int s = 0; s < count; ++s) { ReportStep r; r.length = len; c.schedule.push_back(r); } }
void at(Case& c, int step, const WellControl& wc) { c.schedule[size_t(step)].controls.push_back(wc); }

// ====================================================================== the examples
// ---- BASIC -------------------------------------------------------------------------------------------------------
Case exBuckleyLeverett() {
    Case c = makeBuckleyLeverettCase(100);
    c.title = "Basic 1: 1D Buckley-Leverett waterflood";
    return c;
}

Case exQuarterFiveSpot() {
    Case c; c.title = "Basic 2: Quarter five-spot waterflood (25x25x1)";
    c.units = UnitSystem::Field;
    boxGrid(c, 25, 25, 1, 40, 40, {30}, 5000);
    setLayerProps(c, {150}, {0.2}, 0.1);
    deadOilFluids(c, 1.2, false);
    coreySWOF(c.props.swof, 0.2, 0.25, 2.0, 2.0, 0.4, 1.0, 20);
    c.equil.datumDepth = 5015; c.equil.datumPressure = 3000; c.equil.owc = 6000; c.equil.goc = 0;
    c.wells = {well("INJ", WellType::WaterInjector, 0, 0, 0, 0), well("PROD", WellType::Producer, 24, 24, 0, 0)};
    steps(c, 60, 30.4375);
    at(c, 0, injCtl("INJ", 300, 5000));
    at(c, 0, prodCtl("PROD", Control::BHP, 1500, 1500));
    return c;
}

Case exDepletion() {
    Case c; c.title = "Basic 3: Solution-gas drive depletion (15x15x3)";
    c.units = UnitSystem::Field;
    boxGrid(c, 15, 15, 3, 400, 400, {20, 30, 40}, 8000);
    setLayerProps(c, {120, 40, 250}, {0.18, 0.15, 0.22}, 0.1);
    spe1Fluids(c, false); sandstoneScal(c);
    c.equil.datumDepth = 8050; c.equil.datumPressure = 4000;
    c.equil.owc = 9000; c.equil.goc = 7000; c.equil.initialRs = 0.93;      // bubble point ~3015 psia
    c.wells = {well("PROD", WellType::Producer, 7, 7, 0, 2)};
    steps(c, 40, 91.3125);
    at(c, 0, prodCtl("PROD", Control::ORAT, 3000, 500));
    return c;
}

// ---- INTERMEDIATE ------------------------------------------------------------------------------------------------
Case exConing() {
    Case c; c.title = "Intermediate 1: Bottom-water coning with an aquifer (15x15x12)";
    c.units = UnitSystem::Field;
    const int nx = 15, ny = 15, nz = 12;
    const double w[15] = {400, 250, 150, 90, 55, 35, 25, 20, 25, 35, 55, 90, 150, 250, 400};   // refined towards the well
    boxGrid(c, nx, ny, nz, 100, 100, std::vector<double>(nz, 10.0), 6000);
    GridInput& g = c.grid;
    for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        const size_t a = size_t(g.index(i, j, k)); g.dx[a] = w[i]; g.dy[a] = w[j];
    }
    setLayerProps(c, std::vector<double>(nz, 300.0), std::vector<double>(nz, 0.22), 0.1);
    deadOilFluids(c, 2.0, false);
    coreySWOF(c.props.swof, 0.2, 0.2, 2.5, 2.0, 0.5, 1.0, 20);
    c.equil.datumDepth = 6040; c.equil.datumPressure = 2800; c.equil.owc = 6085; c.equil.goc = 0;
    AquiferInput aq; aq.id = 1; aq.type = AquiferType::Fetkovich; aq.depth = 6120; aq.ct = 1e-5;
    aq.volume = 2e9; aq.pi = 30;
    AquiferConnection cn; cn.i1 = 0; cn.i2 = nx - 1; cn.j1 = 0; cn.j2 = ny - 1; cn.k1 = cn.k2 = nz - 1; cn.face = 'Z'; cn.minus = false;
    aq.conns.push_back(cn);
    c.aquifers.push_back(aq);
    c.wells = {well("PROD", WellType::Producer, 7, 7, 0, 2)};
    steps(c, 60, 30.4375);
    at(c, 0, prodCtl("PROD", Control::ORAT, 1200, 800));
    return c;
}

Case exSPE1() {
    Case c; c.title = "Intermediate 2: SPE1 gas injection (Odeh 1981, 10x10x3)";
    c.units = UnitSystem::Field;
    boxGrid(c, 10, 10, 3, 1000, 1000, {20, 30, 50}, 8325);
    setLayerProps(c, {500, 50, 200}, {0.3, 0.3, 0.3}, 1.0);
    spe1Fluids(c, false); spe1Scal(c);
    c.equil.datumDepth = 8400; c.equil.datumPressure = 4800; c.equil.owc = 8450; c.equil.goc = 8300;
    c.equil.initialRs = 1.27;
    WellInput p = well("PROD", WellType::Producer, 9, 9, 2, 2); p.refDepth = 8400;
    WellInput i = well("INJ", WellType::GasInjector, 0, 0, 0, 0); i.refDepth = 8335;
    c.wells = {p, i};
    const double months[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    for (int y = 0; y < 10; ++y) for (double m : months) { ReportStep r; r.length = m; c.schedule.push_back(r); }
    at(c, 0, prodCtl("PROD", Control::ORAT, 20000, 1000));
    at(c, 0, injCtl("INJ", 100000, 9014));
    c.startDay = 1; c.startMonth = 1; c.startYear = 2015;
    return c;
}

// Gas condensate with gas cycling: the black-oil version of SPE3 (Kenyon & Behie 1987) as used by OPM Flow
Case exCondensate() {
    Case c; c.title = "Advanced 5: Gas condensate cycling (SPE3 black-oil, 9x9x4)";
    c.units = UnitSystem::Field;
    c.phases.water = c.phases.oil = c.phases.gas = c.phases.disgas = c.phases.vapoil = true;
    boxGrid(c, 9, 9, 4, 293.3, 293.3, {30, 30, 50, 50}, 7315);
    setLayerProps(c, {130, 40, 20, 150}, {0.13, 0.13, 0.13, 0.13}, 0.1);
    PropsInput& P = c.props;
    P.rock = {3550, 4e-6};
    P.density = {43.33, 62.37, 0.05850};
    P.pvtw = {3427.6, 1.02629, 0.30698e-5, 0.31107, 0.59410e-5};
    // live oil (the condensate that drops out): PVTO of SPE3
    const std::vector<std::pair<double, std::vector<std::array<double, 3>>>> pvto = {
        {0.189473, {{500.0, 1.20936, 0.219}, {1000.0, 1.19352, 0.240}, {1500.0, 1.17997, 0.260}, {2000.0, 1.16819, 0.279}, {2500.0, 1.15780, 0.299}, {3000.0, 1.14853, 0.318}, {3427.6, 1.14135, 0.334}, {3500.0, 1.14019, 0.337}, {4000.0, 1.13263, 0.356}, {4500.0, 1.12574, 0.374}, {5000.0, 1.11941, 0.393}}},
        {0.479365, {{1000.0, 1.40215, 0.149}, {1500.0, 1.37754, 0.164}, {2000.0, 1.35701, 0.179}, {2500.0, 1.33949, 0.194}, {3000.0, 1.32426, 0.208}, {3427.6, 1.31270, 0.220}, {3500.0, 1.31085, 0.222}, {4000.0, 1.29891, 0.237}, {4500.0, 1.28817, 0.251}, {5000.0, 1.27845, 0.265}}},
        {0.826985, {{1500.0, 1.62448, 0.113}, {2000.0, 1.58913, 0.124}, {2500.0, 1.56020, 0.135}, {3000.0, 1.53584, 0.146}, {3427.6, 1.51777, 0.156}, {3500.0, 1.51492, 0.157}, {4000.0, 1.49666, 0.168}, {4500.0, 1.48052, 0.179}, {5000.0, 1.46610, 0.190}}},
        {1.250165, {{2000.0, 1.89110, 0.089}, {2500.0, 1.84197, 0.098}, {3000.0, 1.80234, 0.107}, {3427.6, 1.77380, 0.114}, {3500.0, 1.76937, 0.115}, {4000.0, 1.74128, 0.124}, {4500.0, 1.71695, 0.132}, {5000.0, 1.69557, 0.141}}},
        {1.794854, {{2500.0, 2.23421, 0.073}, {3000.0, 2.16656, 0.080}, {3427.6, 2.11974, 0.086}, {3500.0, 2.11259, 0.087}, {4000.0, 2.06805, 0.093}, {4500.0, 2.03039, 0.100}, {5000.0, 1.99794, 0.106}}},
        {2.549781, {{3000.0, 2.71411, 0.071}, {3427.6, 2.63221, 0.076}, {3500.0, 2.61998, 0.077}, {4000.0, 2.54537, 0.083}, {4500.0, 2.48418, 0.089}, {5000.0, 2.43271, 0.094}}},
        {2.951016, {{3427.6, 2.96669, 0.068}, {3500.0, 2.94777, 0.069}, {4000.0, 2.83495, 0.075}, {4500.0, 2.74550, 0.080}, {5000.0, 2.67220, 0.085}}},
        {3.033715, {{3500.0, 3.01897, 0.068}, {4000.0, 2.90343, 0.073}, {4500.0, 2.81181, 0.079}, {5000.0, 2.73675, 0.084}}},
        {3.605023, {{4000.0, 3.38017, 0.065}, {4500.0, 3.27351, 0.070}, {5000.0, 3.18612, 0.075}}},
    };
    P.pvto.clear();
    for (auto& [rs, rows] : pvto) {
        PVTORecord o; o.rs = rs; o.pb = rows[0][0]; o.bo = rows[0][1]; o.mu = rows[0][2];
        for (size_t j = 1; j < rows.size(); ++j) o.under.push_back({rows[j][0], rows[j][1], rows[j][2]});
        P.pvto.push_back(o);
    }
    // PVTG: saturated Rv vs pressure (dew point 3427.6 psia), undersaturated branch at 4000 psia
    const double sat[][4] = {{500, 0.0382886993, 6.419987, 0.012999}, {1000, 0.0314227763, 3.007969, 0.014122},
                             {1500, 0.0357629796, 1.922602, 0.015848}, {2000, 0.0473304978, 1.408289, 0.018550},
                             {2500, 0.0679314005, 1.123337, 0.022691}, {3000, 0.1043958616, 0.959755, 0.029087},
                             {3427.6, 0.1670518252, 0.893267, 0.038650}, {3500, 0.1670518252, 0.884653, 0.039158},
                             {4000, 0.1670518252, 0.834785, 0.042494}};
    for (auto& r : sat) { PVTGRecord g; g.pg = r[0]; g.rows.push_back({r[1], r[2], r[3]}); P.pvtg.push_back(g); }
    const double und[][3] = {{0.1043958616, 0.809080, 0.035345}, {0.0679314005, 0.798204, 0.031649}, {0.0473304978, 0.792752, 0.029766},
                             {0.0357629796, 0.789094, 0.028868}, {0.0314227763, 0.786258, 0.028754}, {0.0, 0.773146, 0.026778}};
    for (auto& u : und) P.pvtg.back().rows.push_back({u[0], u[1], u[2]});
    P.pvdg.clear(); for (auto& g : P.pvtg) P.pvdg.push_back({g.pg, g.rows.front().B, g.rows.front().mu});
    const double sg[][4] = {{0, 0, 0.8, 0}, {0.04, 0.005, 0.65, 0}, {0.08, 0.013, 0.513, 0}, {0.12, 0.026, 0.4, 0}, {0.16, 0.04, 0.315, 0},
                            {0.2, 0.058, 0.25, 0}, {0.24, 0.078, 0.196, 0}, {0.28, 0.1, 0.15, 0}, {0.32, 0.126, 0.112, 0}, {0.36, 0.156, 0.082, 0},
                            {0.4, 0.187, 0.06, 0}, {0.44, 0.222, 0.04, 0}, {0.48, 0.26, 0.024, 0}, {0.52, 0.3, 0.012, 0}, {0.56, 0.348, 0.005, 0},
                            {0.6, 0.4, 0, 0}, {0.64, 0.45, 0, 0}, {0.68, 0.505, 0, 0}, {0.72, 0.562, 0, 0}, {0.76, 0.62, 0, 0},
                            {0.8, 0.68, 0, 0}, {0.84, 0.74, 0, 0}};
    const double sw[][4] = {{0.16, 0, 0.8, 50}, {0.2, 0.002, 0.65, 32}, {0.24, 0.01, 0.513, 21}, {0.28, 0.02, 0.4, 15.5},
                            {0.32, 0.033, 0.315, 12}, {0.36, 0.049, 0.25, 9.2}, {0.4, 0.066, 0.196, 7}, {0.44, 0.09, 0.15, 5.3},
                            {0.48, 0.119, 0.112, 4.2}, {0.52, 0.15, 0.082, 3.4}, {0.56, 0.186, 0.06, 2.7}, {0.6, 0.227, 0.04, 2.1},
                            {0.64, 0.277, 0.024, 1.7}, {0.68, 0.33, 0.012, 1.3}, {0.72, 0.39, 0.005, 1}, {0.76, 0.462, 0, 0.7},
                            {0.8, 0.54, 0, 0.5}, {0.84, 0.62, 0, 0.4}, {0.88, 0.71, 0, 0.3}, {0.92, 0.8, 0, 0.2}, {0.96, 0.9, 0, 0.1}, {1, 1, 0, 0}};
    P.sgof.clear(); for (auto& r : sg) P.sgof.push_back({r[0], r[1], r[2], r[3]});
    P.swof.clear(); for (auto& r : sw) P.swof.push_back({r[0], r[1], r[2], r[3]});
    // gas directly on water: GOC = OWC = datum, Rv saturated at the contact
    c.equil.datumDepth = 7500; c.equil.datumPressure = 3550; c.equil.owc = 7500; c.equil.goc = 7500;
    c.equil.initialRs = -1; c.equil.initialRv = -1;
    WellInput p = well("PROD", WellType::Producer, 6, 6, 2, 3, 2.0); p.refDepth = 7375;
    WellInput i = well("INJ", WellType::GasInjector, 0, 0, 0, 1, 2.0); i.refDepth = 7315;
    c.wells = {p, i};
    const double months[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    for (int y = 0; y < 15; ++y) for (double m : months) { ReportStep r; r.length = m; c.schedule.push_back(r); }
    at(c, 0, prodCtl("PROD", Control::GRAT, 6200, 500));
    at(c, 0, injCtl("INJ", 4700, 4000));       // ten years of lean-gas recycling ...
    at(c, 120, injCtl("INJ", 0, 4000));        // ... then blow-down
    c.startDay = 1; c.startMonth = 1; c.startYear = 2015;
    return c;
}

Case exLayeredXsec() {
    Case c; c.title = "Intermediate 3: Layered cross-section, gravity and viscous fingering (30x1x10)";
    c.units = UnitSystem::Field;
    boxGrid(c, 30, 1, 10, 50, 500, std::vector<double>(10, 8.0), 5000);
    setLayerProps(c, {500, 50, 800, 20, 300, 100, 1000, 10, 400, 200},
                  {0.24, 0.17, 0.26, 0.14, 0.22, 0.19, 0.27, 0.12, 0.23, 0.2}, 0.05);
    deadOilFluids(c, 4.0, false);
    coreySWOF(c.props.swof, 0.18, 0.22, 2.0, 2.5, 0.5, 1.0, 20);
    c.equil.datumDepth = 5040; c.equil.datumPressure = 2500; c.equil.owc = 6000; c.equil.goc = 0;
    c.wells = {well("INJ", WellType::WaterInjector, 0, 0, 0, 9), well("PROD", WellType::Producer, 29, 0, 0, 9)};
    steps(c, 72, 30.4375);
    at(c, 0, injCtl("INJ", 800, 6000));
    at(c, 0, prodCtl("PROD", Control::BHP, 1500, 1500));
    return c;
}

Case exPolymerFlood() {
    Case c; c.title = "Advanced: Polymer flood of a viscous-oil layered section (30x1x10)";
    c.units = UnitSystem::Field;
    boxGrid(c, 30, 1, 10, 50, 500, std::vector<double>(10, 8.0), 5000);
    setLayerProps(c, {500, 50, 800, 20, 300, 100, 1000, 10, 400, 200},
                  {0.24, 0.17, 0.26, 0.14, 0.22, 0.19, 0.27, 0.12, 0.23, 0.2}, 0.05);
    deadOilFluids(c, 25.0, false);
    coreySWOF(c.props.swof, 0.18, 0.22, 2.0, 2.5, 0.5, 1.0, 20);
    c.equil.datumDepth = 5040; c.equil.datumPressure = 2500; c.equil.owc = 6000; c.equil.goc = 0;
    c.wells = {well("INJ", WellType::WaterInjector, 0, 0, 0, 9), well("PROD", WellType::Producer, 29, 0, 0, 9)};
    steps(c, 96, 30.4375);
    at(c, 0, injCtl("INJ", 600, 4500));
    at(c, 0, prodCtl("PROD", Control::BHP, 1500, 1500));
    // polymer (FIELD units: lb/stb, rock density lb/rb), data after opm-tests POLYMER-02
    c.phases.polymer = true;
    PolymerInput& P = c.polymer;
    const double cs = 0.350507;
    P.plyvisc = {{0, 1, 0.5 * cs, 3, 1 * cs, 6, 1.5 * cs, 12, 2 * cs, 24, 3 * cs, 48}};
    P.plyads = {{0, 0, 0.25 * cs, 1.2e-5, 0.5 * cs, 1.6e-5, 1 * cs, 2e-5, 1.5 * cs, 2.3e-5, 1.75 * cs, 2.5e-5, 3 * cs, 2.5e-5}};
    PolymerRock R; R.deadPoreVolume = 0.05; R.rrf = 1.5; R.rockDensity = 2500 * cs; R.adsIndex = 2; R.maxAdsorption = 2.5e-5;
    P.plyrock = {R};
    P.cmax = {3 * cs}; P.cmaxSalt = {0}; P.mixing = {0.7};
    PolymerInput::Shear sh; sh.refConc = 1.5 * cs;
    for (double v : {1e-6, 1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0}) { sh.rows.push_back(v * 3.28084); }
    const double mult[] = {1.0, 0.95, 0.85, 0.7, 0.55, 0.45, 0.4};   // shear thinning
    std::vector<double> rows;
    for (size_t k = 0; k < sh.rows.size(); ++k) { rows.push_back(sh.rows[k]); rows.push_back(mult[k]); }
    sh.rows = rows;
    P.plyshlog = {sh};
    c.schedule[6].wpolymer.push_back({"INJ", 1.0 * cs});    // 0.5 yr: start a 2-year polymer slug (1 kg/sm3)
    c.schedule[30].wpolymer.push_back({"INJ", 0.0});        // 2.5 yr: chase water
    return c;
}

Case exGasCapAquifer() {
    Case c; c.title = "Intermediate 4: Dipping oil rim - gas cap + Fetkovich aquifer (30x10x6)";
    c.units = UnitSystem::Field;
    const int nx = 30, ny = 10, nz = 6;
    boxGrid(c, nx, ny, nz, 200, 300, std::vector<double>(nz, 15.0), 7800);
    GridInput& g = c.grid;
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) g.tops[size_t(j) * nx + i] = 7800 + 20.0 * i;   // ~5.7 deg dip
    g.newtran = true;
    heterogeneous(c, 11, 200, 300, 2500, std::vector<double>(nz, 250.0), 0.6, 0.21, 0.025, 0.1);
    spe1Fluids(c, false); sandstoneScal(c);
    c.equil.datumDepth = 8000; c.equil.datumPressure = 3600; c.equil.owc = 8250; c.equil.goc = 7930;
    c.equil.initialRs = -1;                                       // saturated at the GOC
    AquiferInput aq; aq.id = 1; aq.type = AquiferType::Fetkovich; aq.depth = 8450; aq.ct = 1e-5; aq.volume = 3e9; aq.pi = 80;
    AquiferConnection cn; cn.i1 = cn.i2 = nx - 1; cn.j1 = 0; cn.j2 = ny - 1; cn.k1 = 0; cn.k2 = nz - 1; cn.face = 'X'; cn.minus = false;
    aq.conns.push_back(cn);
    c.aquifers.push_back(aq);
    c.wells = {well("P1", WellType::Producer, 12, 3, 2, 5), well("P2", WellType::Producer, 12, 7, 2, 5)};
    steps(c, 60, 91.3125);
    at(c, 0, prodCtl("P1", Control::ORAT, 2500, 1200));
    at(c, 0, prodCtl("P2", Control::ORAT, 2500, 1200));
    // uncertainty parameters for the tornado / Monte Carlo studies
    StudyParam a; a.kind = ParamKind::OWC; a.base = 8250; a.low = 8200; a.high = 8300; a.log = false;
    StudyParam b; b.kind = ParamKind::AquiferStrength; b.target = 0; b.base = 1; b.low = 0.2; b.high = 5;
    StudyParam d; d.kind = ParamKind::PermZ; d.base = 1; d.low = 0.1; d.high = 3;
    StudyParam e; e.kind = ParamKind::KrwMax; e.base = 1; e.low = 0.5; e.high = 1.5;
    StudyParam f; f.kind = ParamKind::PermXY; f.base = 1; f.low = 0.5; f.high = 2;
    c.studyParams = {a, b, d, e, f};
    return c;
}

// ---- ADVANCED ----------------------------------------------------------------------------------------------------
Case exSPE9Like() {
    Case c; c.title = "Advanced 1: SPE9-type heterogeneous dipping reservoir, 26 wells (24x25x15)";
    c.units = UnitSystem::Field;
    const int nx = 24, ny = 25, nz = 15;
    const std::vector<double> dz = {20, 15, 26, 15, 16, 14, 8, 8, 18, 12, 19, 18, 20, 50, 100};
    boxGrid(c, nx, ny, nz, 300, 300, dz, 9000);
    GridInput& g = c.grid;
    const double dip = std::tan(10.0 * kPi / 180.0) * 300.0;
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) g.tops[size_t(j) * nx + i] = 9000 + dip * i;
    g.newtran = true;   // displaced columns connect through their real overlaps (OPM Flow behaviour, so both engines agree)
    heterogeneous(c, 9, 300, 300, 2400, {100, 50, 300, 30, 150, 80, 500, 20, 200, 60, 250, 40, 120, 300, 80}, 1.1,
                  0.13, 0.03, 0.1);
    spe1Fluids(c, false); sandstoneScal(c);
    c.equil.datumDepth = 9035; c.equil.datumPressure = 3600; c.equil.owc = 9950; c.equil.goc = 8800;
    c.equil.initialRs = -1;
    WellInput inj = well("INJ", WellType::WaterInjector, 23, 24, 10, 14);
    c.wells.push_back(inj);
    int n = 0;
    for (int a = 0; a < 5; ++a)
        for (int b = 0; b < 5; ++b) {
            int i = 1 + 3 * a + (b % 2), j = 1 + 5 * b;      // columns 1..14: completions stay above the OWC on the dip
            if (i >= nx) i = nx - 1;
            if (i == 23 && j == 24) continue;
            char name[16]; std::snprintf(name, sizeof(name), "P%02d", ++n);
            c.wells.push_back(well(name, WellType::Producer, i, j, 1, 3));
        }
    steps(c, 30, 30);
    at(c, 0, injCtl("INJ", 5000, 5500));
    for (size_t w = 1; w < c.wells.size(); ++w) {
        at(c, 0, prodCtl(c.wells[w].name, Control::ORAT, 1500, 1000));
        at(c, 10, prodCtl(c.wells[w].name, Control::ORAT, 100, 1000));    // day 300: rate cut (as in SPE9)
        at(c, 12, prodCtl(c.wells[w].name, Control::ORAT, 1500, 1000));   // day 360: back to full rate
    }
    return c;
}

Case exFaultedAnticline() {
    Case c; c.title = "Advanced 2: Faulted corner-point anticline with gas cap (40x20x5, METRIC)";
    c.units = UnitSystem::Metric;
    const int nx = 40, ny = 20, nz = 5;
    const double DX = 100, DY = 100;
    const double dz[5] = {6, 10, 8, 12, 6};
    GridInput& g = c.grid;
    g.nx = nx; g.ny = ny; g.nz = nz;
    const int n = g.size();
    // vertical pillars
    g.coord.resize(size_t(6) * (nx + 1) * (ny + 1));
    for (int j = 0; j <= ny; ++j) for (int i = 0; i <= nx; ++i) {
        double* p = &g.coord[size_t(6) * (j * (nx + 1) + i)];
        p[0] = p[3] = i * DX; p[1] = p[4] = j * DY; p[2] = 1800; p[5] = 2300;
    }
    const double Lx = nx * DX, Ly = ny * DY;
    auto topAt = [&](double x, double y) {
        const double u = (x - 0.5 * Lx) / (0.5 * Lx), v = (y - 0.5 * Ly) / (0.5 * Ly);
        return 1980.0 + 70.0 * u * u + 45.0 * v * v;
    };
    // fault throws: main fault between i = 19 | 20 (east block down 20 m), minor fault j = 11 | 12 in the west block
    auto throwOf = [&](int i, int j) { return (i >= 20 ? 20.0 : 0.0) + (i < 20 && j >= 12 ? 8.0 : 0.0); };
    g.zcorn.resize(size_t(8) * n);
    for (int k = 0; k < nz; ++k) {
        double above = 0; for (int q = 0; q < k; ++q) above += dz[q];
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
            const double t = throwOf(i, j);
            for (int b = 0; b < 2; ++b) for (int jy = 0; jy < 2; ++jy) for (int ix = 0; ix < 2; ++ix) {
                const double z = topAt((i + ix) * DX, (j + jy) * DY) + t + above + (b ? dz[k] : 0.0);
                g.zcorn[(size_t(2 * k + b) * 2 * ny + size_t(2 * j + jy)) * 2 * nx + size_t(2 * i + ix)] = z;
            }
        }
    }
    g.dx.clear(); g.dy.clear(); g.dz.clear(); g.tops.clear();
    g.poro.assign(n, 0.2); g.permx.assign(n, 100); g.permy.assign(n, 100); g.permz.assign(n, 10);
    heterogeneous(c, 21, DX, DY, 1500, {400, 150, 600, 250, 80}, 0.7, 0.22, 0.03, 0.1);
    FaultInput f1; f1.name = "MAIN"; f1.mult = 0.05;
    { FaultFace ff; ff.i1 = ff.i2 = 19; ff.j1 = 0; ff.j2 = ny - 1; ff.k1 = 0; ff.k2 = nz - 1; ff.face = 'X'; f1.faces.push_back(ff); }
    FaultInput f2; f2.name = "MINOR"; f2.mult = 0.5;
    { FaultFace ff; ff.i1 = 0; ff.i2 = 19; ff.j1 = ff.j2 = 11; ff.k1 = 0; ff.k2 = nz - 1; ff.face = 'Y'; f2.faces.push_back(ff); }
    g.faults = {f1, f2};
    spe1Fluids(c, true); sandstoneScal(c);
    c.equil.datumDepth = 2030; c.equil.datumPressure = 205; c.equil.owc = 2060; c.equil.goc = 1995;
    c.equil.initialRs = -1;
    c.wells = {well("P1", WellType::Producer, 13, 9, 0, 4, 0.2), well("P2", WellType::Producer, 26, 9, 0, 4, 0.2),
               well("I1", WellType::WaterInjector, 2, 9, 2, 4, 0.2), well("I2", WellType::WaterInjector, 37, 9, 2, 4, 0.2)};
    steps(c, 32, 91.3125);
    at(c, 0, prodCtl("P1", Control::ORAT, 600, 100));
    at(c, 0, prodCtl("P2", Control::ORAT, 600, 100));
    at(c, 0, injCtl("I1", 800, 350));
    at(c, 0, injCtl("I2", 800, 350));
    return c;
}

Case exChannels() {
    Case c; c.title = "Advanced 3: Channelised fluvial reservoir, line-drive waterflood (60x60x3, METRIC)";
    c.units = UnitSystem::Metric;
    const int nx = 60, ny = 60, nz = 3;
    boxGrid(c, nx, ny, nz, 25, 25, {4, 4, 4}, 1500);
    GridInput& g = c.grid;
    Rng r(33);
    for (int k = 0; k < nz; ++k) {
        struct Ch { double y0, amp, lambda, phase, width; };
        std::vector<Ch> chans;
        for (int m = 0; m < 4; ++m)
            chans.push_back({r.range(150, 1350), r.range(60, 220), r.range(500, 1100), r.range(0, 2 * kPi), r.range(55, 110)});
        SmoothField F(77 + uint64_t(k), 400);
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
            const double x = (i + 0.5) * 25, y = (j + 0.5) * 25;
            double inside = 0;
            for (const Ch& ch : chans) {
                const double yc = ch.y0 + ch.amp * std::sin(2 * kPi * x / ch.lambda + ch.phase);
                const double d = std::fabs(y - yc) / (0.5 * ch.width);
                inside = std::max(inside, d < 1 ? 1.0 - d * d * 0.3 : 0.0);
            }
            const double f = F(x, y);
            const size_t a = size_t(g.index(i, j, k));
            const double kk = inside > 0 ? 1200.0 * inside * std::exp(0.3 * f) : 15.0 * std::exp(0.5 * f);
            g.permx[a] = g.permy[a] = kk; g.permz[a] = 0.1 * kk;
            g.poro[a] = inside > 0 ? 0.24 + 0.03 * inside : 0.13 + 0.01 * f;
        }
    }
    spe1Fluids(c, true); sandstoneScal(c);
    c.equil.datumDepth = 1506; c.equil.datumPressure = 180; c.equil.owc = 1600; c.equil.goc = 1400;
    c.equil.initialRs = 0.636 * 178.1076;                          // bubble point ~139 bar: undersaturated start
    const int js[5] = {6, 18, 30, 42, 54};
    for (int m = 0; m < 5; ++m) {
        char a[16], b[16]; std::snprintf(a, sizeof(a), "I%d", m + 1); std::snprintf(b, sizeof(b), "P%d", m + 1);
        c.wells.push_back(well(a, WellType::WaterInjector, 1, js[m], 0, 2, 0.2));
        c.wells.push_back(well(b, WellType::Producer, 58, js[m], 0, 2, 0.2));
    }
    steps(c, 24, 91.3125);
    for (auto& w : c.wells)
        at(c, 0, w.type == WellType::Producer ? prodCtl(w.name, Control::LRAT, 300, 120) : injCtl(w.name, 320, 300));
    return c;
}

Case exMultiRegion() {
    Case c; c.title = "Advanced 4: Two compartments - EQLNUM/SATNUM regions across a sealing fault (24x12x5, METRIC)";
    c.units = UnitSystem::Metric;
    const int nx = 24, ny = 12, nz = 5;
    boxGrid(c, nx, ny, nz, 100, 100, std::vector<double>(nz, 6.0), 2040);
    GridInput& g = c.grid;
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) g.tops[size_t(j) * nx + i] = 2040 + 1.5 * i;
    heterogeneous(c, 5, 100, 100, 1200, {300, 120, 450, 200, 90}, 0.5, 0.22, 0.025, 0.1);
    const int n = g.size();
    g.eqlnum.assign(n, 1); g.satnum.assign(n, 1); g.fipnum.assign(n, 1);
    for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 12; i < nx; ++i) {
        const size_t a = size_t(g.index(i, j, k));
        g.eqlnum[a] = g.satnum[a] = g.fipnum[a] = 2;
        g.permx[a] *= 0.4; g.permy[a] *= 0.4; g.permz[a] *= 0.4;   // tighter facies in the east block
    }
    FaultInput f; f.name = "SEAL"; f.mult = 0.001;
    { FaultFace ff; ff.i1 = ff.i2 = 11; ff.j1 = 0; ff.j2 = ny - 1; ff.k1 = 0; ff.k2 = nz - 1; ff.face = 'X'; f.faces.push_back(ff); }
    g.faults = {f};
    spe1Fluids(c, true);
    coreySWOF(c.props.swof, 0.18, 0.2, 2.0, 2.0, 0.45, 1.0, 20);
    coreySGOF(c.props.sgof, 0.18, 0.05, 0.15, 2.0, 2.0, 0.85, 20);
    SatRegion s2;
    coreySWOF(s2.swof, 0.30, 0.28, 3.0, 2.5, 0.3, 0.9, 20);
    coreySGOF(s2.sgof, 0.30, 0.08, 0.2, 2.5, 2.5, 0.7, 20);
    c.props.satExtra = {s2};
    c.equil.datumDepth = 2055; c.equil.datumPressure = 210; c.equil.owc = 2068; c.equil.goc = 2000;
    c.equil.initialRs = 0.636 * 178.1076;
    EquilInput e2 = c.equil;                                        // east block: deeper contact, over-pressured
    e2.datumPressure = 228; e2.owc = 2085;
    c.equilExtra = {e2};
    c.wells = {well("PW", WellType::Producer, 5, 6, 0, 2, 0.2), well("PE", WellType::Producer, 17, 6, 0, 2, 0.2),
               well("IE", WellType::WaterInjector, 22, 6, 3, 4, 0.2)};
    steps(c, 28, 91.3125);
    at(c, 0, prodCtl("PW", Control::ORAT, 350, 110));
    at(c, 0, prodCtl("PE", Control::ORAT, 250, 110));
    at(c, 4, injCtl("IE", 300, 320));                              // injection starts after one year
    return c;
}

// Naturally fractured carbonate (ECLIPSE dual porosity): 15 x 15 blocks, 3 matrix + 3 fracture layers. The water
// injected into the fractures imbibes into the water-wet matrix (capillary pressure in SATNUM 1); the fractures use
// straight-line relative permeability without capillary pressure (SATNUM 2). Sigma = 4 (1/lx^2 + 1/ly^2 + 1/lz^2)
// for 12 x 12 x 8 ft matrix blocks (Kazemi), fracture permeability 4000 mD x fracture porosity 0.006 = 24 mD effective.
Case exDualPorosity() {
    Case c; c.title = "Advanced 7: Naturally fractured carbonate - dual porosity waterflood (15x15x3 blocks)";
    c.units = UnitSystem::Field;
    const int nx = 15, ny = 15, nzm = 3;
    boxGrid(c, nx, ny, 2 * nzm, 100, 100, std::vector<double>(2 * nzm, 30.0), 7000);
    GridInput& g = c.grid;
    const int nl = nx * ny, n = g.size();
    g.dual = 1;
    g.satnum.assign(n, 1);
    for (int q = 0; q < n; ++q) {
        const bool frac = q >= nzm * nl;
        g.poro[size_t(q)] = frac ? 0.006 : 0.14 + 0.02 * ((q / nl) % nzm);
        g.permx[size_t(q)] = g.permy[size_t(q)] = frac ? 4000.0 : 0.8;
        g.permz[size_t(q)] = frac ? 1500.0 : 0.3;
        if (frac) g.satnum[size_t(q)] = 2;
    }
    g.sigma = 4.0 * (1.0 / 144 + 1.0 / 144 + 1.0 / 64);
    deadOilFluids(c, 1.8, false);
    coreySWOF(c.props.swof, 0.22, 0.30, 3.0, 2.0, 0.35, 0.9, 20);
    for (auto& r : c.props.swof) {                 // water-wet matrix: spontaneous-imbibition capillary pressure
        const double se = std::max(0.0, std::min(1.0, (r.s - 0.22) / (1.0 - 0.22 - 0.30)));
        r.pc = 12.0 * std::pow(1.0 - se, 2.5) - 0.5 * se;
    }
    SatRegion fr;                                   // fractures: straight lines, no capillary pressure
    fr.swof = {{0.0, 0.0, 1.0, 0.0}, {1.0, 1.0, 0.0, 0.0}};
    c.props.satExtra = {fr};
    c.equil.datumDepth = 7045; c.equil.datumPressure = 3200; c.equil.owc = 7300; c.equil.goc = 0;
    c.wells = {well("INJ", WellType::WaterInjector, 0, 0, 3, 5), well("PROD", WellType::Producer, 14, 14, 3, 5)};
    steps(c, 40, 91.3125);
    at(c, 0, injCtl("INJ", 1500, 4500));
    at(c, 0, prodCtl("PROD", Control::LRAT, 1500, 1000));
    return c;
}

// History-matching exercise: the hidden truth has a high-permeability east half and stronger water relative
// permeability; the base model does not. Observations (WOPR, WWCT, WBHP) are produced by running the truth.
Case hmBase(bool truth) {
    Case c; c.title = "Advanced 5: History-matching exercise (15x15x3) - find the hidden truth";
    c.units = UnitSystem::Field;
    boxGrid(c, 15, 15, 3, 300, 300, {20, 20, 20}, 7000);
    heterogeneous(c, 3, 300, 300, 2500, {150, 60, 200}, 0.5, 0.2, 0.02, 0.1);
    deadOilFluids(c, 1.5, false);
    coreySWOF(c.props.swof, 0.2, 0.2, 2.0, 2.0, 0.5, 1.0, 20);
    if (truth) {
        GridInput& g = c.grid;
        for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 8; i < g.nx; ++i) {
            const size_t a = size_t(g.index(i, j, k)); g.permx[a] *= 2.5; g.permy[a] *= 2.5; g.permz[a] *= 2.5;
        }
        for (auto& r : c.props.swof) r.krw = std::min(1.0, r.krw * 1.3);
    }
    c.equil.datumDepth = 7030; c.equil.datumPressure = 3200; c.equil.owc = 8000; c.equil.goc = 0;
    c.wells = {well("INJ", WellType::WaterInjector, 7, 7, 0, 2), well("PW", WellType::Producer, 1, 7, 0, 2),
               well("PE", WellType::Producer, 13, 7, 0, 2)};
    steps(c, 24, 91.3125);
    at(c, 0, injCtl("INJ", 4000, 6000));
    at(c, 0, prodCtl("PW", Control::LRAT, 2000, 1000));
    at(c, 0, prodCtl("PE", Control::LRAT, 2000, 1000));
    return c;
}

Case exHistoryMatch() {
    Case c = hmBase(false);
    Case t = hmBase(true);
    Simulator s;
    s.setup(t);
    s.initialize();
    while (s.runReportStep()) {
        const SummaryPoint& S = s.summary().back();
        for (size_t w = 0; w < t.wells.size() && w < S.wells.size(); ++w) {
            if (t.wells[w].type != WellType::Producer) continue;
            const WellRates& R = S.wells[w];
            const double liq = R.oil + R.water;
            c.observations.push_back({ObsQty::WOPR, t.wells[w].name, S.time, R.oil, -1});
            c.observations.push_back({ObsQty::WWCT, t.wells[w].name, S.time, liq > 0 ? R.water / liq : 0.0, -1});
            c.observations.push_back({ObsQty::WBHP, t.wells[w].name, S.time, R.bhp, -1});
        }
    }
    StudyParam a; a.kind = ParamKind::PermXY; a.label = "PERM west (i 1-8)"; a.i1 = 0; a.i2 = 7; a.base = 1; a.low = 0.3; a.high = 4;
    StudyParam b; b.kind = ParamKind::PermXY; b.label = "PERM east (i 9-15)"; b.i1 = 8; b.i2 = 14; b.base = 1; b.low = 0.3; b.high = 4;
    StudyParam k; k.kind = ParamKind::KrwMax; k.base = 1; k.low = 0.5; k.high = 2;
    StudyParam v; v.kind = ParamKind::WaterVisc; v.base = 1; v.low = 0.7; v.high = 1.4;
    c.studyParams = {a, b, k, v};
    return c;
}


// ---- local grid refinement: near-well LGRs (graded, layered) and a sub-grid high-permeability streak
Case exLocalRefinement() {
    Case c; c.title = "Advanced: near-well local grid refinement (21x21x4 + 3 LGRs)";
    c.units = UnitSystem::Field;
    boxGrid(c, 21, 21, 4, 60, 60, {15, 10, 20, 15}, 6000);
    setLayerProps(c, {250, 40, 400, 120}, {0.22, 0.16, 0.24, 0.2}, 0.1);
    deadOilFluids(c, 2.0, false);
    coreySWOF(c.props.swof, 0.2, 0.25, 2.0, 2.0, 0.4, 1.0, 20);
    c.equil.datumDepth = 6030; c.equil.datumPressure = 3200; c.equil.owc = 7000; c.equil.goc = 0;
    c.wells = {well("INJ", WellType::WaterInjector, 2, 2, 0, 3), well("PROD", WellType::Producer, 18, 18, 0, 3)};
    GridInput& g = c.grid;
    auto graded = [](int hosts, int per, int wellHost) {   // HXFIN: geometric widths, smallest next to the well
        std::vector<double> w;
        for (int h = 0; h < hosts; ++h)
            for (int m = 0; m < per; ++m) {
                const double d = std::abs((h * per + m) - (wellHost * per + 0.5 * (per - 1)));
                w.push_back(std::pow(1.5, d));
            }
        return w;
    };
    // producer: 3 x 3 host columns, 5 x 5 children each, graded towards the well, every layer split in two
    LgrInput P; P.name = "LPROD"; P.i1 = 17; P.i2 = 19; P.j1 = 17; P.j2 = 19; P.k1 = 0; P.k2 = 3;
    P.nx = 15; P.ny = 15; P.nz = 8; P.hxfin = graded(3, 5, 1); P.hyfin = P.hxfin;
    g.lgrs.push_back(P);
    // injector: 3 x 3 x 1 refinement of its column neighbourhood
    LgrInput I; I.name = "LINJ"; I.i1 = 1; I.i2 = 3; I.j1 = 1; I.j2 = 3; I.k1 = 0; I.k2 = 3; I.nx = 9; I.ny = 9; I.nz = 4;
    g.lgrs.push_back(I);
    // a 20 ft wide, 10 ft thick channel inside layer 3 between the wells - thinner than one host cell
    LgrInput S; S.name = "STREAK"; S.i1 = 6; S.i2 = 14; S.j1 = 10; S.j2 = 10; S.k1 = 2; S.k2 = 2; S.nx = 9; S.ny = 3; S.nz = 2;
    std::vector<double> kx(size_t(S.size()), NAN);
    for (int i = 0; i < S.nx; ++i) kx[size_t(S.index(i, 1, 1))] = 3000;
    S.props["PERMX"] = kx; S.props["PERMY"] = kx;
    g.lgrs.push_back(S);
    moveWellToLgr(c, 0, 1);
    moveWellToLgr(c, 1, 0);
    steps(c, 40, 45.65625);
    at(c, 0, injCtl("INJ", 1000, 6000));
    at(c, 0, prodCtl("PROD", Control::LRAT, 1000, 1500));
    return c;
}

struct Entry { ExampleInfo info; Case (*make)(); };

// ---- compositional / CO2 storage ------------------------------------------------------------------------------------
// Dedicated CO2STORE model: CO2 injection into a dipping saline aquifer (2-D section), 10 years injection + 40 years
// post-injection. Plume migration up-dip under the seal, dissolution and residual (krg = 0) trapping.
Case exCo2Storage() {
    Case c; c.title = "Advanced: CO2 storage in a dipping saline aquifer (CO2STORE, 40x1x20)";
    c.units = UnitSystem::Metric;
    const int nx = 40, nz = 20;
    boxGrid(c, nx, 1, nz, 50, 1000, std::vector<double>(nz, 5.0), 1200);
    GridInput& g = c.grid;
    const double dip = std::tan(2.0 * M_PI / 180.0);
    for (int i = 0; i < nx; ++i) g.tops[size_t(i)] = 1270.0 - dip * 50.0 * (i + 0.5);   // rises towards the closed up-dip end
    std::vector<double> kL(nz, 400.0), pL(nz, 0.24);
    for (int k = 0; k < nz; ++k) { kL[size_t(k)] = 250.0 + 350.0 * (0.5 + 0.5 * std::sin(1.7 * k)); pL[size_t(k)] = 0.20 + 0.06 * (0.5 + 0.5 * std::sin(1.7 * k)); }
    setLayerProps(c, kL, pL, 0.1);
    g.multpv.assign(size_t(g.size()), 1.0);
    for (int k = 0; k < nz; ++k) g.multpv[size_t(g.index(0, 0, k))] = 2000.0;   // open down-dip boundary (large aquifer volume)
    c.phases = Phases{true, false, true, false};
    PropsInput& p = c.props;
    p.rock = {150.0, 4.5e-5};
    coreySWOF(p.swof, 0.2, 0.0, 3.0, 2.0, 1.0, 0.0, 20);   // brine: Corey 3 from Swc = 0.2
    coreySGOF(p.sgof, 0.2, 0.05, 0.0, 2.0, 2.0, 0.7, 20);  // CO2: Corey 2, critical Sg 0.05 (residual trapping)
    for (auto& r : p.swof) r.kro = 0;
    for (auto& r : p.sgof) r.kro = 0;
    c.comp.enabled = true; c.comp.co2store = true; c.comp.structured = true;
    c.comp.fluid = co2StoreFluid(273.15 + 55.0, 1.0);
    c.equil.datumDepth = 1300; c.equil.datumPressure = 130; c.equil.owc = 900; c.equil.goc = 900;   // all brine (contact above the top)
    WellInput inj = well("CO2INJ", WellType::GasInjector, 4, 0, 14, 19, 0.2);
    c.wells = {inj};
    for (int y = 0; y < 10; ++y) { ReportStep r; r.length = 365.25; c.schedule.push_back(r); }
    for (int y = 0; y < 8; ++y) { ReportStep r; r.length = 5 * 365.25; c.schedule.push_back(r); }
    at(c, 0, injCtl("CO2INJ", 160000, 260));   // ~300 t/d of CO2
    WellControl stop = injCtl("CO2INJ", 0, 260); stop.status = WellStatus::Shut;
    at(c, 10, stop);
    c.startYear = 2030;
    return c;
}

// CO2 storage with every trapping mechanism and their interactions: caprock seal (structural), Land hysteresis
// (residual), dissolution + diffusion (solubility), salt precipitation (dry-out) and kinetic carbonation (mineral).
Case exCo2Trapping() {
    Case c = exCo2Storage();
    c.title = "Advanced: CO2 trapping mechanisms - structural, residual, solubility, mineral, salt (CO2STORE)";
    GridInput& g = c.grid;
    const int nx = g.nx, nz = g.nz;
    for (int i = 0; i < nx; ++i)
        for (int k = 0; k < 2; ++k) {   // two-layer caprock: tight, low porosity
            const size_t q = size_t(g.index(i, 0, k));
            g.permx[q] = 0.01; if (!g.permy.empty()) g.permy[q] = 0.01; if (!g.permz.empty()) g.permz[q] = 0.001; g.poro[q] = 0.08;
        }
    (void)nz;
    c.props.hyst.on = true; c.props.hyst.model = 0;                  // Carlson scanning curves
    c.comp.fluid.salinity = 2.0;                                      // ~105 000 ppm NaCl
    co2::TrappingOptions& T = c.comp.trap;
    T.sgtMax = 0.25;                                                  // Land: Sgt,max = 0.25
    T.diffusion = true;
    T.saltPrecip = true;
    co2::MineralReaction m = co2::mineralPresets().front(); m.area = 20.0;
    T.minerals = {m};
    c.wells[0].completions[0].k1 = 12;                                // perforate the lower part only (below the seal)
    return c;
}

// Compositional EOS model: SPE5 six-component fluid, gas injection (lean gas stream) into an undersaturated oil.
Case exCompositional() {
    Case c; c.title = "Advanced: Compositional lean-gas injection (SPE5 fluid, PR EOS, 7x7x3)";
    c.units = UnitSystem::Field;
    boxGrid(c, 7, 7, 3, 500, 500, {20, 30, 50}, 8325);
    setLayerProps(c, {500, 50, 200}, {0.3, 0.3, 0.3}, 0.1);
    spe1Fluids(c, false);
    c.phases = Phases{true, true, true, false};
    sandstoneScal(c, 0.2);
    c.props.rock = {4000.0, 5e-6};
    c.comp.enabled = true; c.comp.structured = true;
    c.comp.fluid = spe5CompFluid();
    c.comp.streams = {{"LEANGAS", {0.77, 0.20, 0.03, 0.0, 0.0, 0.0}}};
    // E300 FIELDSEP semantics: the last stage is the stock tank -> explicit 2-stage train (separator + stock tank)
    c.comp.separator = {{eos::fahrenheitToK(80), eos::psiaToPa(315)}, {eos::fahrenheitToK(60), eos::psiaToPa(14.696)}};
    c.equil.datumDepth = 8400; c.equil.datumPressure = 4000; c.equil.owc = 8500; c.equil.goc = 7000;
    WellInput p = well("PROD", WellType::Producer, 6, 6, 2, 2); p.refDepth = 8400;
    WellInput i = well("INJ", WellType::GasInjector, 0, 0, 0, 0); i.refDepth = 8335;
    c.wells = {p, i};
    steps(c, 40, 91.3125);
    at(c, 0, prodCtl("PROD", Control::ORAT, 3000, 1000));
    at(c, 4, injCtl("INJ", 8000, 4500));
    c.schedule[4].injStreams.push_back({"INJ", "LEANGAS"});
    return c;
}


// ------------------------------------------------------------------ E300 THERMAL examples (milestone C7)
// Heavy-oil thermal fluid (METRIC): HEAVY (non-volatile) + LIGHT (distillable fraction, Crookston K-values of a C10-like
// pseudo-component), OILVISCT, IAPWS-IF97 water / steam, rock heat capacity / conductivity, heat loss to the over- and
// underburden.
void thermalHeavyOil(Case& c, double tInit) {
    spe1Fluids(c, true);                                  // placeholders (the thermal fluid replaces the PVT tables)
    c.phases = Phases{true, true, true, false};
    ThermalInput& T = c.thermal;
    T = ThermalInput(); T.enabled = true; T.structured = true;
    ThermalComponent h; h.name = "HEAVY"; h.mw = 550; h.tcrit = 1e5; h.pcrit = 1; h.dead = true;
    h.dref = 995; h.cref = 7e-5; h.thermex1 = 6.5e-4; h.specha = 2.1; h.zfactor = 0.96;
    ThermalComponent l; l.name = "LIGHT"; l.mw = 180; l.tcrit = 650; l.pcrit = 18; l.dead = false;
    l.kvcr[0] = 0; l.kvcr[1] = 83740; l.kvcr[2] = 0; l.kvcr[3] = 5537; l.kvcr[4] = 0;
    l.dref = 820; l.cref = 1e-4; l.thermex1 = 9e-4; l.specha = 2.2; l.heatvap = 290; l.heatvape = 0.38; l.zfactor = 0.96;
    T.comps = {h, l};
    T.oilViscT = {{10, 8000, 8}, {40, 1500, 4}, {80, 200, 2}, {120, 45, 1.2}, {160, 15, 0.8}, {200, 6, 0.6}, {250, 3, 0.45}, {300, 1.6, 0.35}};
    T.gasViscT = {{10, 0.008, 0.008}, {300, 0.016, 0.016}};
    T.stdT = 15.56; T.stdP = 1.01325;
    T.heatcr = {2347}; T.thconr = {150}; T.thconsf = {0};
    T.heatLoss = true; T.obCv = 2347; T.obK = 150;
    T.zi = {0.75, 0.25};
    T.tempi = {tInit};
    c.props.rock = {20.0, 5e-5};
}
// enumerated initial state of a box model: hydrostatic pressure (water gradient) from (datum depth, pressure), uniform Sw
// (above connate, so that the water is mobile: a cold heavy-oil reservoir at Swc has almost no injectivity)
void thermalEnumInit(Case& c, double datum, double pDatum, double sw) {
    const GridInput& g = c.grid;
    const int nxy = g.nx * g.ny, n = g.size();
    InitialState& I = c.init;
    I.pressure.assign(size_t(n), pDatum); I.swat.assign(size_t(n), sw); I.sgas.assign(size_t(n), 0.0);
    for (int q = 0; q < nxy; ++q) {
        double d = g.tops[size_t(q)];
        for (int k = 0; k < g.nz; ++k) {
            const double dz = g.dz[size_t(k) * nxy + q], mid = d + 0.5 * dz;
            I.pressure[size_t(k) * nxy + q] = pDatum + 1000.0 * 9.80665 * (mid - datum) / 1e5;
            d += dz;
        }
    }
}
void winjtemp(Case& c, int step, const std::string& well, double quality, double tC) {
    ThermalWinj w; w.well = well; w.step = step; w.quality = quality; w.temp = tC;
    w.pres = steam::psat(tC + 273.15) / 1e5;             // saturated: pressure explicit (ECLIPSE gets the same state)
    c.thermal.winj.push_back(w);
}

// Steamflood: quarter of an inverted five-spot, 80 % quality steam at 240 C into the lower two layers
Case exSteamflood() {
    Case c; c.title = "Thermal: Steamflood of a heavy-oil reservoir (quarter five-spot, E300 THERMAL)";
    c.units = UnitSystem::Metric;
    boxGrid(c, 9, 9, 3, 8, 8, {4, 4, 4}, 300);
    setLayerProps(c, {3000, 3000, 3000}, {0.33, 0.33, 0.33}, 0.5);
    sandstoneScal(c, 0.15);
    thermalHeavyOil(c, 25);
    c.equil.datumDepth = 300; c.equil.datumPressure = 15; c.equil.owc = 400; c.equil.goc = 200;
    WellInput p = well("PROD", WellType::Producer, 8, 8, 0, 2, 0.18); p.refDepth = 300;
    WellInput i = well("INJ", WellType::WaterInjector, 0, 0, 1, 2, 0.18); i.refDepth = 304;
    c.wells = {p, i};
    steps(c, 20, 91.3125);                                // 5 years, quarterly reports
    at(c, 0, prodCtl("PROD", Control::BHP, 5, 5));
    at(c, 0, injCtl("INJ", 60, 40));
    winjtemp(c, 0, "INJ", 0.8, 240);
    return c;
}

// Cyclic steam stimulation: one vertical well (injector and producer completed in the same column) - 3 cycles of
// 30 days steam, 10 days soak and 325 days production; 174 x 174 m drainage area refined to 2 m at the well, mobile
// initial water (Sw 0.25)
Case exCss() {
    Case c; c.title = "Thermal: Cyclic steam stimulation (huff and puff), 3 cycles (E300 THERMAL)";
    c.units = UnitSystem::Metric;
    boxGrid(c, 11, 11, 4, 20, 20, {3, 3, 3, 3}, 350);
    {   // geometric refinement towards the well column (5, 5): 2 m near-well cells, 40 m at the boundary (174 m x 174 m)
        const double w[11] = {40, 25, 12, 6, 3, 2, 3, 6, 12, 25, 40};
        GridInput& g = c.grid;
        for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) {
            const size_t q = size_t(k) * g.nx * g.ny + size_t(j) * g.nx + size_t(i);
            g.dx[q] = w[i]; g.dy[q] = w[j];
        }
    }
    setLayerProps(c, {2500, 2500, 2500, 2500}, {0.32, 0.32, 0.32, 0.32}, 0.5);
    sandstoneScal(c, 0.15);
    thermalHeavyOil(c, 20);
    c.equil.datumDepth = 350; c.equil.datumPressure = 30; c.equil.owc = 450; c.equil.goc = 250;
    thermalEnumInit(c, 350, 30, 0.25);
    WellInput p = well("PROD", WellType::Producer, 5, 5, 0, 3, 0.18); p.refDepth = 350;
    WellInput i = well("STEAM", WellType::WaterInjector, 5, 5, 0, 3, 0.18); i.refDepth = 350;
    c.wells = {p, i};
    const double len[] = {30, 10, 100, 100, 125};
    for (int cyc = 0; cyc < 3; ++cyc) for (double d : len) { ReportStep r; r.length = d; c.schedule.push_back(r); }
    for (int cyc = 0; cyc < 3; ++cyc) {
        const int s0 = cyc * 5;
        WellControl inj = injCtl("STEAM", 150, 80);
        WellControl pshut = prodCtl("PROD", Control::BHP, 5, 5); pshut.status = WellStatus::Shut;
        WellControl ishut = injCtl("STEAM", 0, 80); ishut.status = WellStatus::Shut;
        at(c, s0, inj); at(c, s0, pshut);
        at(c, s0 + 1, ishut);                             // soak: both shut
        at(c, s0 + 2, prodCtl("PROD", Control::BHP, 5, 5));
    }
    winjtemp(c, 0, "STEAM", 0.8, 280);
    return c;
}

// SAGD: 2D cross-section normal to a 500 m horizontal well pair (injector 6 m above the producer)
Case exSagd() {
    Case c; c.title = "Thermal: SAGD well pair, 2D cross-section (E300 THERMAL)";
    c.units = UnitSystem::Metric;
    boxGrid(c, 21, 1, 12, 2, 500, std::vector<double>(12, 2.0), 400);
    setLayerProps(c, std::vector<double>(12, 4000), std::vector<double>(12, 0.33), 0.6);
    sandstoneScal(c, 0.15);
    thermalHeavyOil(c, 12);
    c.equil.datumDepth = 400; c.equil.datumPressure = 25; c.equil.owc = 500; c.equil.goc = 300;
    thermalEnumInit(c, 400, 25, 0.25);                  // mobile initial water: start-up communication between the wells
    WellInput p = well("PROD", WellType::Producer, 10, 0, 10, 10, 0.18); p.refDepth = 421; p.completions[0].dir = 'Y';
    WellInput i = well("INJ", WellType::WaterInjector, 10, 0, 7, 7, 0.18); i.refDepth = 415; i.completions[0].dir = 'Y';
    c.wells = {p, i};
    steps(c, 16, 91.3125);                                // 4 years
    at(c, 0, prodCtl("PROD", Control::LRAT, 500, 20));   // liquid-rate limit
    c.thermal.traps.push_back({"PROD", 0, 10.0});        // steam trap: produced fluid >= 10 C subcooled (BHP floor)
    at(c, 0, injCtl("INJ", 400, 30));                     // 400 sm3/d CWE, BHP <= 30 bar
    winjtemp(c, 0, "INJ", 0.9, 232);
    return c;
}

// ------------------------------------------------------------------ C8: chemistry, combustion, UCG, THM
namespace {
ThermalComponent gasComp(const char* name, double mw, double tc, double pc, double cpg) {
    ThermalComponent q; q.name = name; q.mw = mw; q.tcrit = tc; q.pcrit = pc; q.dead = false; q.type = 1;
    q.kvcr[0] = 1e4; q.dref = 800; q.spechg = cpg; q.specha = cpg; q.zfactor = 1.0;
    return q;
}
ThermalComponent solidComp(const char* name, double mw, double rho, double cp) {
    ThermalComponent q; q.name = name; q.mw = mw; q.tcrit = 1e5; q.pcrit = 1; q.dead = true; q.type = 2;
    q.dref = rho; q.specha = cp; q.zfactor = 0.96;
    return q;
}
ThermalReaction reaction(const char* name, size_t nk, std::initializer_list<std::pair<int, double>> sto,
                         std::initializer_list<std::pair<int, double>> prod, std::initializer_list<std::pair<int, double>> order,
                         double A, double Ea, double dH, int o2pp = -1) {
    ThermalReaction r; r.name = name; r.sto.assign(nk, 0.0); r.prod.assign(nk, 0.0); r.order.assign(nk, 0.0);
    for (auto& x : sto) r.sto[size_t(x.first)] = x.second;
    for (auto& x : prod) r.prod[size_t(x.first)] = x.second;
    for (auto& x : order) r.order[size_t(x.first)] = x.second;
    r.A = A; r.Ea = Ea; r.dH = dH; r.o2pp = o2pp;
    return r;
}
// cells (i, all j, all k) of a box: hot ignition zone (dry, gas-filled)
void ignite(Case& c, int i0, int i1, double tC, double sg) {
    const GridInput& g = c.grid;
    ThermalInput& T = c.thermal;
    const double t0 = T.tempi.empty() ? 20.0 : T.tempi[0];
    if (T.tempi.size() != size_t(g.size())) T.tempi.assign(size_t(g.size()), t0);
    for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = i0; i <= i1; ++i) {
        const size_t q = size_t(g.index(i, j, k));
        T.tempi[q] = tC; c.init.swat[q] = 0.0; c.init.sgas[q] = sg;
    }
}
}  // namespace

// In-situ combustion (dry forward combustion) of a heavy oil in a laboratory combustion tube (1.8 m, adiabatic - guard
// heaters): STARS-style 4-reaction kinetics (thermal cracking to coke, heavy / light oil and coke oxidation with O2
// partial-pressure kinetics), air injection 30 sm3/m2/h, electric igniter at the inlet
Case exCombustion() {
    Case c; c.title = "Thermal: In-situ combustion tube (dry forward combustion, 4-reaction kinetics)";
    c.units = UnitSystem::Metric;
    boxGrid(c, 36, 1, 1, 0.05, 0.1, {0.1}, 500);
    setLayerProps(c, {2000}, {0.32}, 0.5);
    sandstoneScal(c, 0.15);
    thermalHeavyOil(c, 50);
    ThermalInput& T = c.thermal;
    // components: HEAVY LIGHT N2 O2 CO2 COKE (+ water)
    T.comps.push_back(gasComp("N2", 28.01, 126.2, 33.9, 1.04));
    T.comps.push_back(gasComp("O2", 32.0, 154.6, 50.4, 0.92));
    T.comps.push_back(gasComp("CO2", 44.01, 304.1, 73.8, 0.84));
    T.comps.push_back(solidComp("COKE", 13.0, 1200, 1.0));
    T.oilViscT = {{10, 8000, 8, 1, 1, 1, 1}, {40, 1500, 4, 1, 1, 1, 1}, {80, 200, 2, 1, 1, 1, 1}, {120, 45, 1.2, 1, 1, 1, 1},
                  {160, 15, 0.8, 1, 1, 1, 1}, {200, 6, 0.6, 1, 1, 1, 1}, {250, 3, 0.45, 1, 1, 1, 1}, {300, 1.6, 0.35, 1, 1, 1, 1}};
    T.gasViscT = {{10, 0.008, 0.008, 0.018, 0.02, 0.015, 0.02}, {800, 0.03, 0.03, 0.04, 0.045, 0.035, 0.04}};
    T.zi = {0.8, 0.2, 0, 0, 0, 0};
    const size_t nk = 7;   // nc + 1
    T.reactions = {
        reaction("CRACK", nk, {{0, 1}}, {{1, 2}, {5, 14.6}}, {{0, 1}}, 1.2e5, 1.3e5, -5e4),
        reaction("HEAVYBURN", nk, {{0, 1}, {3, 56.52}}, {{4, 40.37}, {6, 32.30}}, {{0, 1}, {3, 1}}, 1.0, 8.0e4, 2.37e7, 3),
        reaction("LIGHTBURN", nk, {{1, 1}, {3, 19.25}}, {{4, 12.832}, {6, 12.832}}, {{1, 1}, {3, 1}}, 1.0, 8.0e4, 8.1e6, 3),
        reaction("COKEBURN", nk, {{5, 1}, {3, 1.245}}, {{4, 1}, {6, 0.49}}, {{5, 1}, {3, 1}}, 0.6, 8.0e4, 5.2e5, 3)};
    T.heatLoss = false;                                   // adiabatic tube
    CompStream air; air.name = "AIR"; air.z = {0, 0, 0.79, 0.21, 0, 0};
    T.streams = {air};
    thermalEnumInit(c, 500, 20, 0.2);
    const int n = c.grid.size();
    c.init.sgas.assign(size_t(n), 0.15);
    T.ymf.assign(size_t(n) * 6, 0.0);
    for (int q = 0; q < n; ++q) T.ymf[size_t(2) * size_t(n) + size_t(q)] = 1.0;   // initial gas: N2
    T.xmf.assign(size_t(n) * 6, 0.0);
    for (int q = 0; q < n; ++q) { T.xmf[size_t(q)] = 0.8; T.xmf[size_t(n) + size_t(q)] = 0.2; }
    ignite(c, 0, 1, 400, 0.5);                            // electric igniter: hot, dry start-up zone at the inlet
    WellInput p = well("PROD", WellType::Producer, 35, 0, 0, 0, 0.02); p.completions[0].dir = 'X';
    WellInput i = well("AIRINJ", WellType::GasInjector, 0, 0, 0, 0, 0.02); i.completions[0].dir = 'X';
    c.wells = {p, i};
    steps(c, 12, 1.0 / 24);                               // 12 hours, hourly reports
    at(c, 0, prodCtl("PROD", Control::BHP, 19.5, 19.5));
    at(c, 0, injCtl("AIRINJ", 7.2, 60));                  // 0.3 sm3/h = 30 sm3/m2/h
    c.schedule[0].injStreams.push_back({"AIRINJ", "AIR"});
    winjtemp(c, 0, "AIRINJ", 0, 50);
    T.winj.back().quality = 0; T.winj.back().pres = 0;
    return c;
}

// Underground coal gasification (UCG): linked-well channel in a coal seam, O2-enriched air; coal pyrolysis (volatiles),
// char combustion / partial oxidation, steam and CO2 (Boudouard) gasification, reversible water-gas shift, gas-phase
// CO / H2 / CH4 burning; coal is a solid component - the burnt cavity opens pore space (Kozeny-Carman permeability)
Case exUcg() {
    Case c; c.title = "Thermal: Underground coal gasification (linked vertical wells, THMC chemistry)";
    c.units = UnitSystem::Metric;
    boxGrid(c, 25, 1, 3, 2, 5, {1.5, 1.5, 1.5}, 600);
    setLayerProps(c, {5, 5, 2000}, {0.85, 0.85, 0.85}, 0.5);   // reactive pore space (coal + void); channel in the base layer
    sandstoneScal(c, 0.10);
    thermalHeavyOil(c, 30);
    ThermalInput& T = c.thermal;
    T.comps.clear();
    T.comps.push_back(solidComp("COAL", 12.5, 1350, 1.3));
    T.comps.push_back(solidComp("CHAR", 12.0, 1500, 1.2));
    T.comps.push_back(gasComp("O2", 32.0, 154.6, 50.4, 0.92));
    T.comps.push_back(gasComp("N2", 28.01, 126.2, 33.9, 1.04));
    T.comps.push_back(gasComp("CO2", 44.01, 304.1, 73.8, 0.84));
    T.comps.push_back(gasComp("CO", 28.01, 132.9, 35.0, 1.04));
    T.comps.push_back(gasComp("H2", 2.016, 33.2, 13.0, 14.3));
    T.comps.push_back(gasComp("CH4", 16.04, 190.6, 46.0, 2.2));
    T.oilViscT.clear();
    T.gasViscT = {{10, 0.02, 0.02, 0.02, 0.018, 0.015, 0.018, 0.009, 0.011}, {1200, 0.05, 0.05, 0.055, 0.05, 0.048, 0.05, 0.025, 0.035}};
    T.zi = {0, 0, 0, 1, 0, 0, 0, 0};
    T.heatcr = {2347}; T.thconr = {130};
    T.solidConc = {1.0165e5, 0, 0, 0, 0, 0, 0, 0};       // 80 % of the bulk volume is coal (mol / m3 of pore volume)
    T.permLaw = 1; T.permMultMax = 1e4;
    const size_t nk = 9;   // COAL CHAR O2 N2 CO2 CO H2 CH4 + water(8)
    std::vector<ThermalReaction> R = {
        reaction("PYROLYSIS", nk, {{0, 1}}, {{1, 0.7}, {7, 0.1}, {6, 0.1}, {5, 0.05}, {8, 0.05}}, {{0, 1}}, 1e2, 1.0e5, -1.5e3),
        reaction("CHARBURN", nk, {{1, 1}, {2, 1}}, {{4, 1}}, {{1, 1}, {2, 1}}, 0.036, 1.3e5, 3.935e5, 2),
        reaction("PARTOX", nk, {{1, 2}, {2, 1}}, {{5, 2}}, {{1, 1}, {2, 1}}, 1.0, 1.6e5, 2.21e5, 2),
        reaction("STEAMGASIF", nk, {{1, 1}, {8, 1}}, {{5, 1}, {6, 1}}, {{1, 1}, {8, 1}}, 1.0, 2.0e5, -1.313e5),
        reaction("BOUDOUARD", nk, {{1, 1}, {4, 1}}, {{5, 2}}, {{1, 1}, {4, 1}}, 1.0, 2.1e5, -1.725e5),
        reaction("SHIFT", nk, {{5, 1}, {8, 1}}, {{4, 1}, {6, 1}}, {{5, 1}, {8, 1}}, 0.25, 1.0e5, 4.11e4),
        reaction("COBURN", nk, {{5, 1}, {2, 0.5}}, {{4, 1}}, {{5, 1}, {2, 1}}, 90, 1.2e5, 2.83e5, 2),
        reaction("H2BURN", nk, {{6, 1}, {2, 0.5}}, {{8, 1}}, {{6, 1}, {2, 1}}, 90, 1.2e5, 2.42e5, 2),
        reaction("CH4BURN", nk, {{7, 1}, {2, 2}}, {{4, 1}, {8, 2}}, {{7, 1}, {2, 1}}, 90, 1.4e5, 8.02e5, 2)};
    R[3].phase.assign(nk, -1); R[3].phase[8] = 2;         // steam gasification: water vapour
    R[5].phase.assign(nk, -1); R[5].phase[8] = 2;
    R[5].reversible = true; R[5].keqA = -4.33; R[5].keqB = 4577.8;   // Moe (1962) shift equilibrium
    T.reactions = R;
    CompStream ox; ox.name = "OXYAIR"; ox.z = {0, 0, 0.4, 0.6, 0, 0, 0, 0};
    T.streams = {ox};
    thermalEnumInit(c, 600, 30, 0.4);
    const int n = c.grid.size();
    c.init.sgas.assign(size_t(n), 0.6);
    T.ymf.assign(size_t(n) * 8, 0.0);
    for (int q = 0; q < n; ++q) T.ymf[size_t(3) * size_t(n) + size_t(q)] = 1.0;   // initial gas: N2
    T.xmf.assign(size_t(n) * 8, 0.0);
    for (int q = 0; q < n; ++q) T.xmf[size_t(3) * size_t(n) + size_t(q)] = 1.0;
    ignite(c, 0, 1, 650, 1.0);                            // ignition zone at the injection well
    WellInput p = well("PROD", WellType::Producer, 24, 0, 2, 2, 0.2);
    WellInput i = well("INJ", WellType::GasInjector, 0, 0, 2, 2, 0.2);
    c.wells = {p, i};
    steps(c, 15, 2);                                      // 30 days
    at(c, 0, prodCtl("PROD", Control::BHP, 25, 25));
    at(c, 0, injCtl("INJ", 1500, 60));
    c.schedule[0].injStreams.push_back({"INJ", "OXYAIR"});
    winjtemp(c, 0, "INJ", 0, 30);
    T.winj.back().quality = 0; T.winj.back().pres = 0;
    return c;
}

// UCG + CO2 storage in one thermal run (C10, unified framework): the UCG channel of exUcg with dissolved gases (IAPWS
// Henry constants), a Peng-Robinson gas, Land gas hysteresis, ash CaO carbonation by aqueous CO2 (mineral trapping) and,
// after 30 days of gasification, CO2 injection into the cavity with the producer shut, then a shut-in period
Case exUcgCo2() {
    Case c = exUcg();
    c.title = "Thermal: UCG syngas production + CO2 storage in the cavity (unified trapping, THMC)";
    ThermalInput& T = c.thermal;
    const int n = c.grid.size();
    const size_t ncOld = T.comps.size();   // 8: COAL CHAR O2 N2 CO2 CO H2 CH4
    T.comps.push_back(solidComp("CAO", 56.08, 3340, 0.75));
    T.comps.push_back(solidComp("CALCITE", 100.09, 2710, 0.84));
    T.comps[0].formula = "C0.85H0.7O0.1";   // coal pseudo-formula consistent with the pyrolysis stoichiometry
    T.comps[1].formula = "C";
    T.comps[8].formula = "CaO";
    T.comps[9].formula = "CaCO3"; T.comps[9].trackEq = 1;
    for (auto& q : T.comps) {
        if (q.type == 2) continue;
        const gasdata::GasInfo* gi = gasdata::find(q.name);
        if (gi && gi->hA != 0) { q.aqModel = 1; q.hA = gi->hA; q.hB = gi->hB; q.hC = gi->hC; q.vInf = gi->vInf; q.ksalt = gi->ksalt; }
    }
    T.gasPR = true; T.salinity = 0.1; T.hystGas = true; T.sgtMax = 0.25; T.track = "CO2"; T.structTol = 0.01;
    for (auto& r : T.gasViscT) { r.push_back(0.02); r.push_back(0.02); }
    T.zi.resize(ncOld + 2, 0.0);
    T.solidConc.resize(ncOld + 2, 0.0); T.solidConc[8] = 500;   // ash CaO, mol / m3 of pore volume
    for (auto& st : T.streams) st.z.resize(ncOld + 2, 0.0);
    T.xmf.resize(size_t(n) * (ncOld + 2), 0.0); T.ymf.resize(size_t(n) * (ncOld + 2), 0.0);
    auto widen = [&](std::vector<double>& v, double fill) { if (v.size() == ncOld + 1) v.insert(v.begin() + long(ncOld), 2, fill); };
    for (auto& r : T.reactions) {
        widen(r.sto, 0.0); widen(r.prod, 0.0); widen(r.order, 0.0);
        if (r.phase.size() == ncOld + 1) r.phase.insert(r.phase.begin() + long(ncOld), 2, -1);
    }
    const size_t nk = ncOld + 3;
    ThermalReaction carb = reaction("CARBONATION", nk, {{4, 1}, {8, 1}}, {{9, 1}}, {{4, 1}, {8, 1}}, 1.0, 5.0e4, 1.6e5);
    carb.phase.assign(nk, -1); carb.phase[4] = 0;   // CO2 reacts from the aqueous phase (dissolved CO2 + ash: CarbFix-type)
    T.reactions.push_back(carb);
    CompStream co2; co2.name = "CO2INJ"; co2.z.assign(ncOld + 2, 0.0); co2.z[4] = 1.0;
    T.streams.push_back(co2);
    // schedule: 30 d gasification (exUcg); 30 d water quench of the hot cavity (oxidant off, steam / syngas produced);
    // 60 d CO2 injection into the cooled cavity (producer shut); 150 d shut-in (dissolution, residual, mineral trapping)
    WellInput q = well("QUENCH", WellType::WaterInjector, 1, 0, 2, 2, 0.2);
    c.wells.push_back(q);
    WellControl qshut = injCtl("QUENCH", 0, 60); qshut.status = WellStatus::Shut;
    at(c, 0, qshut);
    const int s0 = int(c.schedule.size());
    steps(c, 3, 10);
    steps(c, 6, 10);
    steps(c, 5, 30);
    WellControl ishut = injCtl("INJ", 0, 120); ishut.status = WellStatus::Shut;
    at(c, s0, ishut);
    at(c, s0, injCtl("QUENCH", 15, 60));
    winjtemp(c, s0, "QUENCH", 0, 20);
    T.winj.back().quality = -1; T.winj.back().pres = 0;   // cold liquid water
    WellControl qoff = injCtl("QUENCH", 0, 60); qoff.status = WellStatus::Shut;
    WellControl pshut = prodCtl("PROD", Control::BHP, 25, 25); pshut.status = WellStatus::Shut;
    at(c, s0 + 3, qoff); at(c, s0 + 3, pshut);
    at(c, s0 + 3, injCtl("INJ", 80, 120));
    c.schedule[size_t(s0 + 3)].injStreams.push_back({"INJ", "CO2INJ"});
    winjtemp(c, s0 + 3, "INJ", 0, 30);
    T.winj.back().quality = 0; T.winj.back().pres = 0;
    at(c, s0 + 9, ishut);
    return c;
}

// C12: CO2-enhanced coal-bed methane (ECBM) - Langmuir sorption on a shared site species (extended / competitive
// isotherm, van 't Hoff temperature dependence, sorption time constant). One year of primary depletion (dewatering +
// CH4 desorption), then four years of CO2 injection: CO2 (higher affinity) displaces sorbed CH4.
Case exEcbm() {
    Case c; c.title = "Thermal core: CO2-ECBM - competitive Langmuir desorption of coal-bed methane by CO2";
    c.units = UnitSystem::Metric;
    boxGrid(c, 11, 11, 1, 30, 30, {4}, 700);
    setLayerProps(c, {5}, {0.02}, 0.1);                   // cleat system: 2 % porosity, 5 mD
    sandstoneScal(c, 0.10);
    thermalHeavyOil(c, 35);
    ThermalInput& T = c.thermal;
    T.comps = {gasComp("CH4", 16.04, 190.6, 46.0, 2.2), gasComp("CO2", 44.01, 304.1, 73.8, 0.84)};
    T.oilViscT.clear();
    T.gasViscT = {{10, 0.0105, 0.0140}, {300, 0.0190, 0.0260}};
    T.zi = {1, 0};
    T.heatLoss = false;
    // Langmuir (Busch & Gensterblum 2011 ranges): V_L 8.77 m3/t (CH4) x 1.4 t/m3 coal / 0.02241 m3/mol = 548 mol/m3 bulk;
    // P_L 6.69 MPa (CH4), 2.9 MPa (CO2 - selectivity ~2.3); heats of adsorption 15 / 25 kJ/mol; sorption time 2 d
    ThermalSorb a; a.comp = "CH4"; a.capacity = 548; a.b = 1 / 6690.0; a.tRef = 308.15; a.dH = 15e3; a.tau = 2 * 86400.0;
    ThermalSorb b = a; b.comp = "CO2"; b.b = 1 / 2900.0; b.dH = 25e3;
    T.sorb = {a, b};
    thermalEnumInit(c, 700, 70, 0.4);
    const int n = c.grid.size();
    c.init.sgas.assign(size_t(n), 0.6);
    T.ymf.assign(size_t(n) * 2, 0.0); T.xmf.assign(size_t(n) * 2, 0.0);
    for (int q = 0; q < n; ++q) { T.ymf[size_t(q)] = 1.0; T.xmf[size_t(q)] = 1.0; }
    CompStream co2; co2.name = "CO2INJ"; co2.z = {0, 1};
    T.streams = {co2};
    WellInput p = well("PROD", WellType::Producer, 10, 10, 0, 0, 0.2);
    WellInput i = well("INJ", WellType::GasInjector, 0, 0, 0, 0, 0.2);
    c.wells = {p, i};
    steps(c, 4, 91.25);                                   // primary depletion
    steps(c, 16, 91.25);                                  // CO2 injection
    at(c, 0, prodCtl("PROD", Control::BHP, 10, 10));
    WellControl ishut = injCtl("INJ", 0, 150); ishut.status = WellStatus::Shut;
    at(c, 0, ishut);
    at(c, 4, injCtl("INJ", 15000, 150));
    c.schedule[4].injStreams.push_back({"INJ", "CO2INJ"});
    winjtemp(c, 4, "INJ", 0, 35);
    T.winj.back().quality = 0; T.winj.back().pres = 0;
    return c;
}

// EGS doublet with embedded discrete fractures (C12): cold water circulated through a stimulated fracture network in hot
// granite. The main fracture links the wells; two oblique fractures cross it (f-f intersection transmissibility). Barton-
// Bandis normal closure and Coulomb shear dilation respond to the pore pressure and the thermal unloading of the rock.
Case exEgsEdfm() {
    Case c; c.title = "Thermal core: EGS doublet in granite with embedded discrete fractures (EDFM + Barton-Bandis + hydro-shear)";
    c.units = UnitSystem::Metric;
    boxGrid(c, 21, 11, 3, 25, 25, {50, 50, 50}, 3000);
    setLayerProps(c, {0.01, 0.01, 0.01}, {0.01, 0.01, 0.01}, 1.0);   // tight granite: 0.01 mD, 1 % porosity
    sandstoneScal(c, 0.0);
    thermalHeavyOil(c, 200);
    ThermalInput& T = c.thermal;
    T.comps = {gasComp("N2", 28.01, 126.2, 33.9, 1.04)};            // inert trace gas (absent initially): water-only circulation
    T.oilViscT.clear();
    T.gasViscT = {{10, 0.0175}, {300, 0.0290}};
    T.zi = {1};
    T.heatLoss = false;
    T.heatcr = {2400}; T.thconr = {216}; T.thconsf = {0};             // granite: 2.4 MJ/m3/K, 2.5 W/m/K (216 kJ/m/d/K)
    c.props.rock = {300.0, 1e-6};
    thermalEnumInit(c, 3075, 300, 1.0);
    // fractures (metres from the box corner; full height when zTop / zBot < 0): aperture 0.4 mm, k = w^2 / 12
    ThermalFracture f1; f1.name = "MAIN"; f1.x0 = 40; f1.y0 = 137.5; f1.x1 = 485; f1.y1 = 137.5; f1.aperture = 0.4e-3;
    f1.mech = true; f1.kn0 = 20e9; f1.vm = 0.4e-3; f1.mu = 0.6; f1.ks = 2e9; f1.dilationDeg = 5;
    ThermalFracture f2 = f1; f2.name = "OBL1"; f2.x0 = 150; f2.y0 = 40; f2.x1 = 230; f2.y1 = 235; f2.aperture = 0.3e-3;
    ThermalFracture f3 = f1; f3.name = "OBL2"; f3.x0 = 300; f3.y0 = 235; f3.x1 = 380; f3.y1 = 40; f3.aperture = 0.3e-3;
    T.fractures = {f1, f2, f3};
    T.fracStress.sHmaxGrad = 26e3; T.fracStress.sHminGrad = 17e3; T.fracStress.azimuthDeg = 30;
    T.fracStress.E = 40e9; T.fracStress.nu = 0.25; T.fracStress.alphaT = 8e-6;
    WellInput p = well("PROD", WellType::Producer, 18, 5, 0, 2, 0.2);
    WellInput i = well("INJ", WellType::WaterInjector, 2, 5, 0, 2, 0.2);
    c.wells = {p, i};
    steps(c, 20, 91.25);                                  // 5 years, quarterly
    at(c, 0, prodCtl("PROD", Control::BHP, 280, 280));
    at(c, 0, injCtl("INJ", 2000, 450));                   // ~23 kg/s
    winjtemp(c, 0, "INJ", 0, 50);
    T.winj.back().quality = 0; T.winj.back().pres = 0;
    return c;
}

// C13: parametric FORGE-style multi-stage EGS doublet (see EgsDesign in examples.h) - public, outside the anonymous namespace
}  // namespace
bool egsStimulation(const EgsDesign& d, hf::Result& r) {
    hf::Input in;
    in.model = d.stimModel == 1 ? 1 : 0;
    in.E = d.rockE; in.nu = d.rockNu; in.height = std::max(1.0, d.fracHeight); in.mu = d.stimMu; in.rate = d.stimRate;
    in.cl = d.stimCl; in.kic = d.stimKic; in.sigmaMin = d.sHminGrad * d.depth; in.tEnd = d.stimTime;
    in.cells = in.model == 0 ? 300 : 120; in.outputs = 12;
    return hf::simulateAdaptive(in, r);
}
Case makeEgsCase(const EgsDesign& d0) {
    EgsDesign d = d0;
    double stimL = 0, stimAp = 0;
    if (d.stimulate) {
        hf::Result r;
        if (!egsStimulation(d, r)) throw std::runtime_error("EGS stimulation sizing failed: " + r.error);
        stimL = r.last().L;
        stimAp = std::max(1e-6, d.stimRetain * r.meanWidth());
    }
    d.stages = std::max(1, std::min(40, d.stages));
    d.cell = std::max(5.0, d.cell); d.nz = std::max(3, std::min(30, d.nz));
    d.spacing = std::max(2.0 * d.refineFirst, d.spacing);
    Case c;
    char t[200]; std::snprintf(t, sizeof t, "Thermal core: multi-stage EGS doublet (%d EDFM fractures, %.0f m apart, wells %.0f m apart, %.0f C)",
                               d.stages, d.spacing, d.wellSep, d.tRes);
    c.title = t;
    c.units = UnitSystem::Metric;
    const double margX = std::max(100.0, 2 * d.cell), margZ = 50.0;
    const double margY = d.stimulate ? std::max(125.0, stimL + 3 * d.cell) : std::max(125.0, d.fracExt + 3 * d.cell);
    const double lat = d.stages * d.spacing;
    // column widths along x: margins with uniform cells; around every plane a cell of refineFirst holding it, growing by
    // refineRatio up to the coarse size until half way to the next plane
    std::vector<double> wx, fracX;
    auto uniform = [&](double len) { const int m = std::max(1, int(std::ceil(len / d.cell - 1e-9))); for (int k = 0; k < m; ++k) wx.push_back(len / m); };
    uniform(margX);
    for (int s = 0; s < d.stages; ++s) {
        const double half = 0.5 * d.spacing;
        std::vector<double> side;   // widths from the plane outwards, filling half - refineFirst / 2
        if (d.refine) {
            double rem = half - 0.5 * d.refineFirst, w = d.refineFirst * d.refineRatio;
            while (rem > 1e-9) {
                double ww = std::min({w, d.cell, rem});
                if (rem - ww < 0.5 * std::min(w, d.cell)) ww = rem;   // no sliver at the midpoint
                side.push_back(ww); rem -= ww; w *= d.refineRatio;
            }
            for (size_t k = side.size(); k-- > 0;) wx.push_back(side[k]);
            fracX.push_back(margX + (s + 0.5) * d.spacing);
            wx.push_back(d.refineFirst);
            for (double v : side) wx.push_back(v);
        } else {
            const int m = std::max(1, int(std::ceil(d.spacing / d.cell - 1e-9)));
            for (int k = 0; k < m; ++k) wx.push_back(d.spacing / m);
            fracX.push_back(margX + (s + 0.5) * d.spacing + 1e-3 * d.cell);   // off the cell faces
        }
    }
    uniform(margX);
    const int nx = int(wx.size());
    const int ny = std::max(3, int(std::ceil((d.wellSep + 2 * margY) / d.cell)));
    std::vector<double> dzl = {margZ};
    for (int k = 0; k < d.nz - 2; ++k) dzl.push_back(d.fracHeight / (d.nz - 2));
    dzl.push_back(margZ);
    const double H = d.fracHeight + 2 * margZ;
    const double top = d.depth - 0.5 * H;
    boxGrid(c, nx, ny, d.nz, d.cell, d.cell, dzl, top);
    for (size_t q = 0; q < c.grid.dx.size(); ++q) c.grid.dx[q] = wx[q % size_t(nx)];
    setLayerProps(c, std::vector<double>(size_t(d.nz), d.kMatrix), std::vector<double>(size_t(d.nz), d.poro), 1.0);
    sandstoneScal(c, 0.0);
    thermalHeavyOil(c, d.tRes);
    ThermalInput& T = c.thermal;
    T.comps = {gasComp("N2", 28.01, 126.2, 33.9, 1.04)};
    T.oilViscT.clear();
    T.gasViscT = {{10, 0.0175}, {300, 0.0290}};
    T.zi = {1};
    T.heatLoss = false;
    T.heatcr = {2400}; T.thconr = {216}; T.thconsf = {0};
    c.props.rock = {300.0, 1e-6};
    const double pTop = 1.01325 + 1000.0 * 9.80665 * top / 1e5;
    thermalEnumInit(c, top, pTop, 1.0);
    // wells: horizontal along x through the stage region, in the layer at the wells' depth
    const double yI = margY, yP = margY + d.wellSep;
    const int jI = std::min(ny - 1, int(yI / d.cell)), jP = std::min(ny - 1, int(yP / d.cell));
    int kW = 0; { double z = top; for (int k = 0; k < d.nz; ++k) { if (d.depth >= z && d.depth < z + dzl[size_t(k)]) kW = k; z += dzl[size_t(k)]; } }
    int i0 = nx - 1, i1 = 0;
    { double x = 0; for (int i = 0; i < nx; ++i) { const double xc = x + 0.5 * wx[size_t(i)]; if (xc > margX && xc < margX + lat) { i0 = std::min(i0, i); i1 = std::max(i1, i); } x += wx[size_t(i)]; } }
    auto horiz = [&](const std::string& nm, WellType ty, int j) {
        WellInput w; w.name = nm; w.type = ty; w.i = i0; w.j = j; w.group = ty == WellType::Producer ? "PROD" : "INJ";
        for (int i = i0; i <= i1; ++i) { Completion cp; cp.i = i; cp.j = j; cp.k1 = cp.k2 = kW; cp.diameter = 0.2; cp.dir = 'X'; w.completions.push_back(cp); }
        w.refDepth = d.depth;
        return w;
    };
    c.wells = {horiz("PROD", WellType::Producer, jP), horiz("INJ", WellType::WaterInjector, jI)};
    // transverse fractures (planes x = const), centred on the wells' depth
    for (int s = 0; s < d.stages; ++s) {
        ThermalFracture f; f.name = "STAGE" + std::to_string(s + 1);
        f.x0 = f.x1 = fracX[size_t(s)]; f.y0 = yI - d.fracExt; f.y1 = yP + d.fracExt;
        if (d.stimulate) { f.y0 = yI - stimL; f.y1 = yI + stimL; }
        f.zTop = d.depth - 0.5 * d.fracHeight; f.zBot = d.depth + 0.5 * d.fracHeight;
        f.aperture = d.stimulate ? stimAp : d.aperture; f.mech = d.mech; f.kn0 = 20e9; f.vm = 0.4e-3; f.mu = 0.6; f.ks = 2e9; f.dilationDeg = 5;
        T.fractures.push_back(f);
    }
    T.fracStress.sHmaxGrad = d.sHmaxGrad; T.fracStress.sHminGrad = d.sHminGrad; T.fracStress.azimuthDeg = 90;
    T.fracStress.E = d.rockE; T.fracStress.nu = d.rockNu; T.fracStress.alphaT = 8e-6;
    const int nrep = std::max(1, int(std::lround(d.years * 4)));
    steps(c, nrep, 365.25 * d.years / nrep);
    const double pW = 1.01325 + 1000.0 * 9.80665 * d.depth / 1e5;
    at(c, 0, prodCtl("PROD", Control::BHP, pW - d.drawdown, pW - d.drawdown));
    at(c, 0, injCtl("INJ", d.rate, pW + 150));
    winjtemp(c, 0, "INJ", 0, d.tInj);
    T.winj.back().quality = 0; T.winj.back().pres = 0;
    return c;
}
namespace {

// C15: shale gas, multistage hydraulically fractured horizontal well (EDFM fractures in the black-oil grid, Langmuir
// adsorbed gas, stress-dependent fracture conductivity)
Case exShaleMultistage() {
    const char* deck = R"(RUNSPEC
TITLE
Shale gas - 10-stage fractured horizontal well
DIMENS
 61 31 3 /
OIL
WATER
GAS
METRIC
TABDIMS
 1 1 20 20 /
WELLDIMS
 1 60 1 1 /
START
 1 JAN 2025 /
GRID
DX
 5673*20 /
DY
 5673*10 /
DZ
 5673*10 /
TOPS
 1891*2500 1891*2510 1891*2520 /
PORO
 5673*0.06 /
PERMX
 5673*0.0005 /
PERMY
 5673*0.0005 /
PERMZ
 5673*0.0001 /
-- RESSIM: LANGMUIR 2.8 55 2.5
PROPS
PVDO
 20 1.05 1.0
 400 1.00 1.0 /
PVDG
 20 0.0600 0.0140
 50 0.0240 0.0150
 100 0.0120 0.0170
 150 0.0081 0.0190
 200 0.0062 0.0210
 250 0.0051 0.0230
 300 0.0044 0.0250
 350 0.0039 0.0270 /
PVTW
 300 1.0 4e-5 0.4 0 /
DENSITY
 800 1020 0.75 /
ROCK
 300 4e-5 /
SWOF
 0.30 0 1 0
 0.50 0.05 0.3 0
 1.00 1 0 0 /
SGOF
 0 0 1 0
 0.20 0.10 0.2 0
 0.70 0.80 0 0 /
SOLUTION
EQUIL
 2515 300 3000 0 2700 0 /
SCHEDULE
WELSPECS
 HZ1 FIELD 6 16 2515 GAS /
/
COMPDAT
 HZ1 6 16 2 2 OPEN 1* 1* 0.12 3* X /
 HZ1 7 16 2 2 OPEN 1* 1* 0.12 3* X /
/
WCONPROD
 HZ1 OPEN BHP 5* 50 /
/
TSTEP
 10*3 10*30 12*91.25 6*365 /
END
)";
    Case c = eclipseDeckFromString(deck, "");
    WellInput& w = c.wells[0];
    w.completions.clear();
    for (int i = 5; i <= 55; ++i) {   // 1000 m lateral along x in layer 2
        Completion cp; cp.i = i; cp.j = 15; cp.k1 = cp.k2 = 1; cp.dir = 'X'; cp.diameter = 0.12;
        w.completions.push_back(cp);
    }
    hf::StageLayout L;
    L.well = w.name; L.stages = 10; L.halfLength = 110; L.height = 30; L.width = 0.003;
    L.perm = 50000; L.poro = 0.35; L.gamma = 0.002;
    c.grid.hfracs = hf::layoutStages(c, L);
    c.title = "Shale gas - 10-stage fractured horizontal well";
    return c;
}
Case exEgsMultistage() { return makeEgsCase(EgsDesign{}); }

// THM: SAGD with coupled geomechanics - thermal expansion and pore-pressure increase heave the overburden (fixed-stress
// iterative coupling, FEM mesh padded with 400 m overburden, 100 m underburden and 200 m side-burden)
Case exSagdThm() {
    Case c = exSagd();
    c.title = "Thermal + geomechanics (THM): SAGD heave and stress change";
    c.schedule.resize(8);                                 // 2 years
    geomech::Input& G = c.thermal.geo;
    G.enabled = true; G.coupling = 1;
    G.E = {0.5e9}; G.nu = {0.3}; G.biot = {1.0}; G.alphaT = {2e-5};
    G.overLayers = 4; G.overThick = 400; G.underLayers = 2; G.underThick = 100; G.sidePad = 3; G.sideThick = 200;
    G.padE = 1.5e9; G.padNu = 0.3; G.padAlphaT = 1e-5;
    G.sigmaVGrad = 21e3; G.k0 = 0.7; G.cohesion = 0.5e6; G.frictionDeg = 30; G.tensile = 0.2e6;
    c.thermal.permLaw = 2; c.thermal.permC = 30; c.thermal.permMultMax = 50;   // dilation-enhanced permeability
    return c;
}

// CO2 EOR: SPE5 reservoir and oil with CO2 as a 7th EOS component, water-alternating-CO2 (WAG) and CO2SOL (CO2 dissolves
// in the brine), E300 style - no separate CO2 module
Case exCo2Eor() {
    Case c; c.title = "Advanced: CO2 EOR - WAG with CO2 as an EOS component (SPE5 + CO2, CO2SOL, 7x7x3)";
    c.units = UnitSystem::Field;
    boxGrid(c, 7, 7, 3, 500, 500, {20, 30, 50}, 8325);
    setLayerProps(c, {500, 50, 200}, {0.3, 0.3, 0.3}, 0.1);
    spe1Fluids(c, false);
    c.phases = Phases{true, true, true, false};
    sandstoneScal(c, 0.2);
    c.props.rock = {4000.0, 5e-6};
    c.comp.enabled = true; c.comp.structured = true;
    c.comp.fluid = spe5Co2CompFluid(1.0);
    c.comp.streams = {{"CO2", {0, 0, 0, 0, 0, 0, 1.0}}};
    c.comp.separator = {{eos::fahrenheitToK(80), eos::psiaToPa(315)}, {eos::fahrenheitToK(60), eos::psiaToPa(14.696)}};
    c.equil.datumDepth = 8400; c.equil.datumPressure = 4000; c.equil.owc = 8450; c.equil.goc = 7000;
    WellInput p = well("PROD", WellType::Producer, 6, 6, 2, 2); p.refDepth = 8400;
    WellInput i = well("INJ", WellType::GasInjector, 0, 0, 0, 0); i.refDepth = 8335;
    c.wells = {p, i};
    steps(c, 40, 91.3125);   // 10 years, quarterly reports
    at(c, 0, prodCtl("PROD", Control::ORAT, 3000, 1000));
    // year 1: primary depletion; then 1-year CO2 / water half-cycles (WAG ratio 1:1)
    for (int s = 4; s < 40; s += 4) {
        const bool gas = ((s / 4) % 2) == 1;
        WellControl w = injCtl("INJ", gas ? 6000 : 5000, 6000);
        w.phase = gas ? 2 : 1;
        at(c, s, w);
        if (gas) c.schedule[size_t(s)].injStreams.push_back({"INJ", "CO2"});
    }
    return c;
}

const std::vector<Entry>& registry() {
    static const std::vector<Entry> R = {
        {{"bl1d", "1D Buckley-Leverett waterflood", ExampleLevel::Basic,
          "100-cell water injection into a dead-oil column - compare with the analytical front.",
          "Two-phase immiscible displacement, Corey curves, shock front and rarefaction. Watch the water saturation front "
          "move in the 3D view and the water cut at the producer.", "100 x 1 x 1 Cartesian", 0.2}, exBuckleyLeverett},
        {{"qfs", "Quarter five-spot waterflood", ExampleLevel::Basic,
          "Classic pattern flood: injector and producer in opposite corners of a homogeneous layer.",
          "Areal sweep, breakthrough time and the water-cut curve. Try the streamline view to see the sweep pattern and "
          "change the oil viscosity (PVDO) to see the effect of mobility ratio.", "25 x 25 x 1 Cartesian", 0.5}, exQuarterFiveSpot},
        {{"depletion", "Solution-gas drive depletion", ExampleLevel::Basic,
          "One producer depletes an undersaturated live-oil reservoir below its bubble point.",
          "Pressure decline, bubble point crossing, free-gas saturation build-up and the rising producing GOR - the "
          "textbook solution-gas drive signature.", "15 x 15 x 3 Cartesian", 1}, exDepletion},
        {{"coning", "Bottom-water coning", ExampleLevel::Intermediate,
          "A partially completed producer above a strong Fetkovich aquifer; locally refined grid.",
          "Vertical-to-horizontal permeability, critical rate and water coning. Look at a vertical slice through the "
          "well, and at the aquifer influx curve. Try halving the oil rate.", "15 x 15 x 12, refined near the well", 7}, exConing},
        {{"spe1", "SPE1 gas injection (Odeh)", ExampleLevel::Intermediate,
          "The first SPE comparative solution project: dry-gas injection into an undersaturated oil reservoir.",
          "Three-phase black-oil, gas override in the high-permeability top layer, Rs re-dissolution. Results can be "
          "compared with OPM Flow / ECLIPSE (FPR ~3830 psia after 10 years).", "10 x 10 x 3 Cartesian", 1}, exSPE1},
        {{"xsec", "Layered cross-section", ExampleLevel::Intermediate,
          "Water flood of a strongly layered section with an adverse mobility ratio.",
          "Gravity under-running, thief zones and vertical crossflow (kv/kh = 0.05). The 2D cross-section view shows "
          "which layers flood first.", "30 x 1 x 10 Cartesian", 1}, exLayeredXsec},
        {{"gascap", "Dipping oil rim: gas cap + aquifer", ExampleLevel::Intermediate,
          "Combination drive in a dipping reservoir; includes uncertainty parameters for the Studies panel.",
          "Gas-cap expansion from up-dip, aquifer influx from down-dip, coning of both fluids to the producers. Run the "
          "tornado / Monte Carlo studies on OWC, aquifer strength, kv/kh, krw and permeability.", "30 x 10 x 6 Cartesian, 5.7 deg dip, NEWTRAN", 5}, exGasCapAquifer},
        {{"spe9like", "SPE9-type heterogeneous reservoir", ExampleLevel::Advanced,
          "Dipping, log-normal heterogeneous model with 1 injector and 25 producers and a mid-run rate cut.",
          "Field-scale three-phase behaviour, wells switching from rate to BHP control, gas liberation around the "
          "producers and CPR-AMG performance on 9000 cells.", "24 x 25 x 15 Cartesian, 10 deg dip, NEWTRAN", 12}, exSPE9Like},
        {{"faulted", "Faulted corner-point anticline", ExampleLevel::Advanced,
          "Dome-shaped corner-point grid with a 20 m throw main fault (MULTFLT 0.05) and a minor fault.",
          "Corner-point geometry (COORD/ZCORN), non-neighbour connections across fault throws, fault transmissibility "
          "multipliers, a small crestal gas cap and flank water injection in each fault block.", "40 x 20 x 5 corner-point", 4}, exFaultedAnticline},
        {{"channels", "Channelised fluvial reservoir", ExampleLevel::Advanced,
          "Sinuous high-permeability sand channels in a shale background, five-injector line drive.",
          "Strong heterogeneity, early breakthrough along channels, poor sweep of the background. Compare permeability "
          "and water-saturation maps and the streamlines.", "60 x 60 x 3 Cartesian (10 800 cells)", 25}, exChannels},
        {{"regions", "Two compartments with regions", ExampleLevel::Advanced,
          "EQLNUM/SATNUM/FIPNUM regions on both sides of a nearly sealing fault, different contacts and pressures.",
          "Multi-region equilibration (two OWCs and datum pressures), region-specific relative permeability, initial "
          "cross-fault disequilibrium and delayed pressure communication.", "24 x 12 x 5 Cartesian", 3}, exMultiRegion},
        {{"condensate", "Gas condensate cycling (SPE3)", ExampleLevel::Advanced,
          "Wet gas (VAPOIL / PVTG) reservoir below the dew point: condensate drop-out, lean-gas recycling, blow-down.",
          "Black-oil SPE3: all oil starts vaporised in the gas (Rv). Watch FOPR fall as liquid drops out near the producer "
          "(3D: SOIL, RV), the lean injection gas sweeping the rich gas, and the pressure decline after cycling stops at 10 years. "
          "Matches OPM Flow within 0.1 %.", "9 x 9 x 4 Cartesian", 1}, exCondensate},
        {{"polymer", "Polymer flood (EOR)", ExampleLevel::Advanced,
          "Viscous oil (25 cP) in a strongly layered section: water flood with a 2-year polymer slug and chase water.",
          "Mobility control with polymer: Todd-Longstaff mixing, adsorption (irreversible), residual resistance factor, "
          "inaccessible pore volume and shear thinning near the injector. Compare FOPR / FWCT with a pure water flood "
          "(EOR section: untick POLYMER) and look at the POLYMER and PADS maps.", "30 x 1 x 10 Cartesian", 1}, exPolymerFlood},
        {{"dualporo", "Fractured carbonate (dual porosity)", ExampleLevel::Advanced,
          "ECLIPSE-style DUALPORO model: tight water-wet matrix blocks drained by a high-permeability fracture network.",
          "Water races through the fractures (early breakthrough) while oil is recovered slowly by capillary imbibition "
          "from the matrix (Kazemi sigma, T = CDARCY V sigma K). Compare SWAT in the matrix and fracture layers, try a "
          "smaller SIGMA (bigger blocks), or switch to DUALPERM in Grid & rock.", "15 x 15 x (3 matrix + 3 fracture)", 2}, exDualPorosity},
        {{"lgr", "Local grid refinement (LGR)", ExampleLevel::Advanced,
          "Near-well CARFIN refinements around both wells (graded HXFIN/HYFIN, split layers) and a sub-grid channel.",
          "ECLIPSE / OPM-compatible Cartesian LGRs: WELSPECL / COMPDATL wells, NXFIN-style grading towards the producer, a "
          "high-permeability streak thinner than one host cell (local PERMX in CARFIN ... ENDFIN). Compare the water front in the "
          "refined region with the coarse cells, or delete an LGR in Grid & rock and rerun.", "21 x 21 x 4 + 3 LGRs (2178 local cells)", 3}, exLocalRefinement},
        {{"comp", "Compositional lean-gas injection (EOS)", ExampleLevel::Advanced,
          "SPE5 six-component Peng-Robinson fluid; lean gas (WELLSTRE / WINJGAS) injected after one year of depletion.",
          "Fully implicit compositional (natural variables, EOS flash). Open Fluid > Compositional to see the phase envelope, "
          "flash a composition and edit the BIC matrix; in Results look at the oil / gas saturation and the producing GOR "
          "as the injected gas vaporises the intermediate components.", "7 x 7 x 3 Cartesian, 6 components", 15}, exCompositional},
        {{"co2", "CO2 storage in a saline aquifer (CO2STORE)", ExampleLevel::Advanced,
          "300 t/d CO2 for 10 years into a dipping brine aquifer, then 40 years of post-injection migration.",
          "Dedicated CO2STORE model: Spycher-Pruess / Duan-Sun solubility, brine density increase on dissolution, "
          "residual trapping (Sgc 0.05). Results > CO2 storage shows the mobile / trapped / dissolved inventory; the XCO2 "
          "property shows the dissolved plume sinking below the gas tongue.", "40 x 1 x 20 cross-section, 2 deg dip", 20}, exCo2Storage},
        {{"co2trap", "CO2 trapping mechanisms and interactions (CO2STORE)", ExampleLevel::Advanced,
          "The CO2 aquifer with a caprock, gas hysteresis (Land), diffusion, salt precipitation and mineral carbonation.",
          "Every trapping mechanism at once: structural (caprock), residual (Carlson / Land hysteresis, Sgt,max 0.25), "
          "solubility (dissolution + molecular diffusion), mineral (anorthite carbonation) and dry-out salt near the injector. "
          "Fluid > CO2 storage > Trapping & interactions holds all settings; Results > CO2 storage shows the five-class "
          "inventory, and the SGTRAP / SGMAX / MINCO2 / SSALT / PORFAC maps the per-cell state.", "40 x 1 x 20 cross-section, caprock", 45}, exCo2Trapping},
        {{"co2eor", "CO2 EOR: WAG with CO2 as an EOS component (CO2SOL)", ExampleLevel::Advanced,
          "SPE5 oil + CO2 (7-component PR EOS): one year of depletion, then 1-year CO2 / water alternating injection.",
          "CO2 is simply the 7th EOS component - it develops miscibility with the oil by vaporising / condensing exchange and "
          "also dissolves in the brine (CO2SOL, Spycher-Pruess / Duan-Sun Henry function). Results: ZMF / XMF / YMF / AMF of "
          "CO2 per cell, producing GOR and the CO2 breakthrough; Fluid > Compositional > CO2 EOR tools: MMP, swelling test, "
          "WAG schedule generator.", "7 x 7 x 3 Cartesian, 7 components + water", 30}, exCo2Eor},
        {{"steamflood", "Thermal: steamflood of a heavy-oil reservoir (E300 THERMAL)", ExampleLevel::Advanced,
          "Quarter five-spot, 80 % quality steam at 240 C into a 25 C heavy oil (2-component K-value fluid).",
          "ResSim thermal core (STARS-like fully implicit p / N / U formulation): steam zone at Tsat(p), steam distillation of "
          "the LIGHT component, viscosity reduction (OILVISCT), heat loss to over / underburden. Results: TEMP, SGAS (steam), "
          "oil rate response and the cumulative steam-oil ratio.", "9 x 9 x 3 Cartesian, 2 components + water", 80}, exSteamflood},
        {{"css", "Thermal: cyclic steam stimulation (huff and puff)", ExampleLevel::Advanced,
          "Three cycles of 30 d steam injection, 10 d soak and 325 d production from one vertical well.",
          "Injector and producer are completed in the same column and alternate (WCONINJE / WCONPROD + SHUT). Watch the heated "
          "zone grow and cool between cycles and the oil-rate peaks after each soak.", "11 x 11 x 4 Cartesian (refined to 2 m at the well), 2 components + water", 60}, exCss},
        {{"sagd", "Thermal: SAGD well pair (2D cross-section)", ExampleLevel::Advanced,
          "Steam-assisted gravity drainage: 500 m horizontal injector 6 m above the producer, 4 years.",
          "The steam chamber rises and spreads above the injector; heated oil and condensate drain by gravity to the producer. "
          "Results: TEMP / SGAS cross-section and the steam-oil ratio.", "21 x 1 x 12 cross-section, 2 components + water", 120}, exSagd},
        {{"combustion", "Thermal: in-situ combustion (air injection, 4-reaction kinetics)", ExampleLevel::Advanced,
          "Laboratory combustion tube: dry forward combustion of a heavy oil - cracking to coke, oil and coke burning (O2 kinetics).",
          "STARS-style chemistry on the ResSim thermal core (C8): solid COKE component, non-condensable N2 / O2 / CO2, heat of "
          "reaction, ignition zone at the air injector. Results: the combustion front (TEMP > 450 C), O2 consumption, CO2 in the "
          "produced gas, coke deposition ahead of the front.", "36 x 1 x 1 tube (5 cm cells), 6 components + water", 30}, exCombustion},
        {{"ucg", "Thermal: underground coal gasification (UCG)", ExampleLevel::Advanced,
          "Linked-well UCG channel: coal pyrolysis, char combustion / gasification, water-gas shift, syngas production.",
          "9 reactions (pyrolysis, char burning and partial oxidation, steam and Boudouard gasification, reversible shift, CO / H2 / "
          "CH4 burning); coal is a solid component - the cavity grows as it is consumed and the permeability follows Kozeny-"
          "Carman. Results: temperature (> 900 C), cavity porosity, CO / H2 / CH4 rates.", "25 x 1 x 3 cross-section, 8 components + water", 60}, exUcg},
        {{"ucg_co2", "Thermal: UCG syngas + CO2 storage in the cavity (unified trapping)", ExampleLevel::Advanced,
          "UCG channel (pyrolysis, char burning / gasification, shift) for 30 d, then CO2 injected into the hot cavity, then shut-in.",
          "One thermal run combines syngas production and CO2 storage with the same trapping classes as the CO2STORE core: CO2 / "
          "H2 / CO / CH4 dissolve in water (IAPWS G7-04 Henry, Sechenov salting-out), Peng-Robinson gas, Land gas hysteresis "
          "(residual), structural, and mineral trapping by ash CaO carbonation with dissolved CO2 (CALCITE, FORMULA trackEq). "
          "Results: CO2 storage tab, XMF / YMF / AMF / ZMF maps per component over time, SGTRAP, SOLID_CALCITE.",
          "25 x 1 x 3 cross-section, 10 components + water", 120}, exUcgCo2},
        {{"ecbm", "Thermal core: CO2-ECBM (competitive Langmuir sorption, C12)", ExampleLevel::Advanced,
          "Coal seam: one year of dewatering / CH4 desorption, then four years of CO2 injection displacing sorbed methane.",
          "SORB (Langmuir on a shared site species) gives the extended Langmuir isotherm with van 't Hoff temperature dependence "
          "and a sorption time constant; adsorbed gas is ~8x the free cleat gas. Results: SOLID_CH4_ADS / SOLID_CO2_ADS / "
          "SOLID_SITES maps, CH4 rate response and the CO2 breakthrough at the producer, YMF of CO2.",
          "11 x 11 x 1 Cartesian, 2 gases + 3 sorption solids + water", 40}, exEcbm},
        {{"egs_edfm", "Thermal core: EGS doublet with embedded discrete fractures (C12)", ExampleLevel::Advanced,
          "Cold water circulated for 5 years between two wells through a stimulated fracture network in 200 C granite.",
          "EDFM: fracture cells appended to the matrix grid, f-m transmissibility 2 A k / <d>, f-f along traces and at "
          "intersections. Barton-Bandis normal closure and irreversible Coulomb shear dilation under the in-situ stress "
          "(FRACSTRESS) with thermo-elastic unloading. Results: TEMP cooling plume along the fractures, FRACAP aperture (mm), "
          "FRACP / FRACT on the host cells, producer temperature (thermal drawdown).",
          "21 x 11 x 3 Cartesian + 3 fractures, water + trace N2", 30}, exEgsEdfm},
        {{"egs_multistage", "Thermal: multi-stage EGS doublet (FORGE-style, 5 EDFM fractures)", ExampleLevel::Advanced,
          "Two parallel horizontal wells 150 m apart at 2.5 km in 210 C granite, connected by 5 transverse stimulated fractures.",
          "Parametric EGS designer case (Thermal > Fractures & EGS > EGS designer): EDFM fractures with Barton-Bandis closure and \n"
          "Coulomb hydro-shear under an anisotropic stress state; 50 C water circulated for 10 years. Results: Geothermal tab - \n"
          "production temperature / thermal drawdown, net thermal power (MW), cumulative heat, injectivity / impedance, fracture \n"
          "aperture and temperature evolution per fracture, cooled-volume growth; 3D: fracture planes coloured by aperture / T.",
          "44 x 16 x 3 Cartesian (refined at the planes) + 5 fractures", 60}, exEgsMultistage},
        {{"shale_hf", "Unconventional: shale gas, 10-stage fractured horizontal well (EDFM)", ExampleLevel::Advanced,
          "1000 m lateral with 10 transverse propped fractures in a 500 nD shale; Langmuir adsorbed gas; 10 years at 50 bar BHP.",
          "C15 hydraulic fractures in the black-oil grid: each HFRAC plane is cut into fracture cells (EDFM) connected to the \n"
          "matrix, to each other and to the well. Langmuir desorption (V_L 2.8 sm3/t, P_L 55 bar) adds gas as the pressure \n"
          "falls; fracture conductivity declines with drawdown (gamma 0.002 1/bar). Grid & rock > Hydraulic fractures shows \n"
          "the stage designer (PKN / KGD sizing from the pump schedule). Results: gas rate decline (linear then boundary-\n"
          "dominated flow), pressure depletion around the fractures, adsorbed vs free gas in the run log.",
          "61 x 31 x 3 Cartesian + 10 fractures", 20}, exShaleMultistage},
        {{"sagd_thm", "Thermal + geomechanics: SAGD heave (THM)", ExampleLevel::Advanced,
          "SAGD with coupled linear thermo-poroelastic FEM: surface heave, effective-stress change, shear-failure index.",
          "Fixed-stress iterative coupling (Kim-Tchelepi-Juanes): the pore volume follows the mechanics, dilation raises the "
          "permeability. Results: Geomechanics tab - uplift history, vertical displacement / mean effective stress / failure maps.",
          "21 x 1 x 12 + FEM padding, 2 components + water", 90}, exSagdThm},
        {{"hm", "History-matching exercise", ExampleLevel::Advanced,
          "Observed rates, water cuts and pressures from a hidden truth model - recover it with manual or automatic HM.",
          "Open the History Matching panel: the observations and four parameters are pre-loaded. Try the manual sliders "
          "first, then Levenberg-Marquardt or ES-MDA. (Hint: one half of the field is better than it looks.)",
          "15 x 15 x 3 Cartesian", 2}, exHistoryMatch},
    };
    return R;
}

}  // namespace

const char* toString(ExampleLevel l) {
    switch (l) { case ExampleLevel::Basic: return "Basic"; case ExampleLevel::Intermediate: return "Intermediate"; default: return "Advanced"; }
}

std::vector<ExampleInfo> listExamples() {
    std::vector<ExampleInfo> v;
    for (const Entry& e : registry()) v.push_back(e.info);
    return v;
}

Case makeExample(const std::string& id) {
    for (const Entry& e : registry()) if (e.info.id == id) return e.make();
    throw std::runtime_error("Unknown example '" + id + "'");
}

// placeholder fluids for imported decks the ResSim engine cannot model (display only; never run by the own engine)
void placeholderFluids(Case& c, bool metric) { spe1Fluids(c, metric); }


// public: heavy-oil thermal fluid preset (only c.thermal is replaced; values converted to the case unit system)
void thermalPreset(Case& c, double tInitC) {
    Case t; t.units = UnitSystem::Metric;
    thermalHeavyOil(t, tInitC);
    ThermalInput T = t.thermal;
    if (c.units == UnitSystem::Field) {
        auto tF = [](double C) { return C * 1.8 + 32; };
        for (auto& q : T.comps) {
            q.dref *= 0.0624279606; q.cref /= 14.5037738; q.thermex1 /= 1.8;
            q.specha /= 4.1868; q.spechb /= 4.1868 * 1.8; q.spechg /= 4.1868; q.spechh /= 4.1868 * 1.8;
            q.heatvap /= 2.326; q.heatvaps /= 2.326;
            q.kvcr[1] *= 14.5037738; q.kvcr[2] /= 14.5037738; q.kvcr[3] *= 1.8; q.kvcr[4] *= 1.8;
            q.pcrit *= 14.5037738; q.tcrit = q.tcrit >= 1e4 ? q.tcrit : q.tcrit * 1.8;
        }
        for (auto& r : T.oilViscT) r[0] = tF(r[0]);
        for (auto& r : T.gasViscT) r[0] = tF(r[0]);
        T.stdT = 60; T.stdP = 14.696;
        for (double& v : T.heatcr) v /= 67.0661;
        for (double& v : T.thconr) v /= 6.2306;
        T.obCv /= 67.0661; T.obK /= 6.2306;
        for (double& v : T.tempi) v = tF(v);
    }
    c.thermal = T;
}

} // namespace ressim
