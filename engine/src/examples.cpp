// Built-in example library (see examples.h). All geology is generated procedurally with a portable RNG so that
// every platform builds bit-identical cases.
#include "ressim/examples.h"
#include "ressim/simulator.h"
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

struct Entry { ExampleInfo info; Case (*make)(); };

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

} // namespace ressim
