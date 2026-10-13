#include "ressim/case.h"
#include <set>
#include <cmath>
#include <algorithm>

namespace ressim {

Units Units::of(UnitSystem u) {
    Units r;
    if (u == UnitSystem::Metric) {
        r.cDarcy = 0.00852702; r.gravity = 0.0000980665; r.gasDensity = 1.0; r.poreVolume = 1.0; r.atm = 1.01325;
        r.p = "bar"; r.len = "m"; r.liq = "sm3"; r.gas = "sm3"; r.rs = "sm3/sm3"; r.dens = "kg/m3"; r.rv = "sm3/sm3";
    } else {
        r.cDarcy = 0.00112712; r.gravity = 0.00694444; r.gasDensity = 1000.0 / 5.614583; r.poreVolume = 1.0 / 5.614583;
        r.atm = 14.696;
        r.p = "psia"; r.len = "ft"; r.liq = "stb"; r.gas = "Mscf"; r.rs = "Mscf/stb"; r.dens = "lb/ft3"; r.rv = "stb/Mscf";
    }
    return r;
}

const char* toString(DeckSection s) {
    static const char* n[] = {"RUNSPEC", "GRID", "EDIT", "PROPS", "REGIONS", "SOLUTION", "SUMMARY", "SCHEDULE"};
    const int i = int(s); return i >= 0 && i < 8 ? n[i] : "SCHEDULE";
}

const char* toString(Control c) {
    switch (c) {
        case Control::ORAT: return "ORAT"; case Control::WRAT: return "WRAT"; case Control::GRAT: return "GRAT";
        case Control::LRAT: return "LRAT"; case Control::RESV: return "RESV"; case Control::BHP: return "BHP";
        case Control::RATE: return "RATE";
    }
    return "?";
}
const char* toString(WellType t) {
    switch (t) {
        case WellType::Producer: return "Producer"; case WellType::WaterInjector: return "Water injector";
        case WellType::GasInjector: return "Gas injector";
    }
    return "?";
}

static bool increasing(const std::vector<double>& v) {
    for (size_t i = 1; i < v.size(); ++i) if (!(v[i] > v[i - 1])) return false;
    return true;
}

std::string validateCase(const Case& c) {
    const GridInput& g = c.grid;
    if (g.nx < 1 || g.ny < 1 || g.nz < 1) return "GRID: NX, NY, NZ must all be >= 1";
    if (double(g.nx) * g.ny * g.nz > 5e7) return "GRID: more than 50 million cells is not supported";
    const size_t n = size_t(g.size());
    auto chk = [&](const std::vector<double>& v, const char* name, bool req) -> std::string {
        if (v.empty()) return req ? std::string("GRID: ") + name + " is missing" : "";
        if (v.size() != n) return std::string("GRID: ") + name + " must have NX*NY*NZ values";
        for (double x : v) if (!std::isfinite(x) || x < 0) return std::string("GRID: ") + name + " has negative or invalid values";
        return "";
    };
    const bool cp = g.cornerPoint();
    if (cp) {
        if (g.coord.size() != size_t(6) * (g.nx + 1) * (g.ny + 1)) return "GRID: COORD must have 6*(NX+1)*(NY+1) values";
        if (g.zcorn.size() != 8 * n) return "GRID: ZCORN must have 8*NX*NY*NZ values";
        for (double x : g.zcorn) if (!std::isfinite(x)) return "GRID: ZCORN has invalid values";
        for (double x : g.coord) if (!std::isfinite(x)) return "GRID: COORD has invalid values";
    }
    for (auto r : {chk(g.dx, "DX", !cp), chk(g.dy, "DY", !cp), chk(g.dz, "DZ", !cp), chk(g.poro, "PORO", true),
                   chk(g.permx, "PERMX", true), chk(g.permy, "PERMY", false), chk(g.permz, "PERMZ", false),
                   chk(g.ntg, "NTG", false)})
        if (!r.empty()) return r;
    for (double v : g.poro) if (v > 1) return "GRID: PORO values must be <= 1";
    if (!cp && g.tops.size() != size_t(g.nx * g.ny) && g.tops.size() != n) return "GRID: TOPS must have NX*NY or NX*NY*NZ values";
    for (auto& f : g.faults) {
        if (!(f.mult >= 0)) return "GRID: fault " + f.name + " has a negative multiplier";
        for (auto& ff : f.faces)
            if (ff.i1 > ff.i2 || ff.j1 > ff.j2 || ff.k1 > ff.k2 || ff.i1 < 0 || ff.j1 < 0 || ff.k1 < 0 || ff.i2 >= g.nx || ff.j2 >= g.ny || ff.k2 >= g.nz)
                return "GRID: fault " + f.name + " has a face box outside the grid";
    }
    for (auto& e : g.nnc)
        if (e.i1 < 0 || e.i1 >= g.nx || e.j1 < 0 || e.j1 >= g.ny || e.k1 < 0 || e.k1 >= g.nz ||
            e.i2 < 0 || e.i2 >= g.nx || e.j2 < 0 || e.j2 >= g.ny || e.k2 < 0 || e.k2 >= g.nz || e.trans < 0)
            return "GRID: NNC outside the grid or negative transmissibility";
    for (auto& a : c.aquifers) {
        const std::string nm = "AQUIFER " + std::to_string(a.id);
        if (!c.phases.water) return nm + ": aquifers need the water phase";
        if (!(a.ct > 0)) return nm + ": total compressibility must be positive";
        if (a.type == AquiferType::Fetkovich && (!(a.volume > 0) || !(a.pi > 0))) return nm + ": Fetkovich volume and PI must be positive";
        if (a.type == AquiferType::CarterTracy && (!(a.perm > 0) || !(a.poro > 0) || !(a.rInner > 0) || !(a.thickness > 0) || !(a.angle > 0)))
            return nm + ": Carter-Tracy k, phi, r0, h and angle must be positive";
        if (a.conns.empty()) return nm + ": no connections (AQUANCON)";
        for (auto& q : a.conns)
            if (q.i1 > q.i2 || q.j1 > q.j2 || q.k1 > q.k2 || q.i1 < 0 || q.j1 < 0 || q.k1 < 0 || q.i2 >= g.nx || q.j2 >= g.ny || q.k2 >= g.nz)
                return nm + ": connection box outside the grid";
    }
    if (!g.actnum.empty() && g.actnum.size() != n) return "GRID: ACTNUM must have NX*NY*NZ values";

    const Phases& ph = c.phases;
    const PropsInput& p = c.props;
    if (!ph.water && !ph.oil) return "RUNSPEC: water or oil must be active";
    if (ph.water && ph.gas && !ph.oil) return "RUNSPEC: water-gas without oil is not supported";
    if (ph.water && ph.oil) {
        if (p.swof.size() < 2) return "PROPS: SWOF needs at least 2 rows";
        for (size_t i = 1; i < p.swof.size(); ++i) if (!(p.swof[i].s > p.swof[i - 1].s)) return "PROPS: SWOF saturations must increase";
    }
    if (ph.gas && ph.oil) {
        if (p.sgof.size() < 2) return "PROPS: SGOF needs at least 2 rows";
        for (size_t i = 1; i < p.sgof.size(); ++i) if (!(p.sgof[i].s > p.sgof[i - 1].s)) return "PROPS: SGOF saturations must increase";
        if (p.pvdg.size() < 2) return "PROPS: PVDG needs at least 2 rows";
        std::vector<double> pp; for (auto& r : p.pvdg) pp.push_back(r.p);
        if (!increasing(pp)) return "PROPS: PVDG pressures must increase";
    }
    if (ph.gas && ph.vapoil) {
        if (!ph.oil) return "RUNSPEC: VAPOIL needs both OIL and GAS";
        if (p.pvtg.size() < 2) return "PROPS: PVTG needs at least 2 pressure records (VAPOIL)";
        for (size_t i = 0; i < p.pvtg.size(); ++i) {
            const auto& r = p.pvtg[i];
            if (r.rows.empty()) return "PROPS: PVTG record " + std::to_string(i + 1) + " is empty";
            if (i > 0 && !(r.pg > p.pvtg[i - 1].pg)) return "PROPS: PVTG gas pressures must increase";
            for (size_t j = 0; j < r.rows.size(); ++j) {
                if (!(r.rows[j].B > 0) || !(r.rows[j].mu > 0)) return "PROPS: PVTG Bg and viscosity must be positive";
                if (r.rows[j].rv < 0) return "PROPS: PVTG Rv must be >= 0";
                if (j > 0 && !(r.rows[j].rv < r.rows[j - 1].rv)) return "PROPS: PVTG Rv must decrease within a record (saturated row first)";
            }
        }
        bool under = false; for (auto& r : p.pvtg) if (r.rows.size() > 1) under = true;
        if (!under) return "PROPS: at least one PVTG record needs undersaturated rows (Rv below saturation)";
    }
    if (ph.oil) {
        if (ph.gas && ph.disgas) {
            if (p.pvto.empty()) return "PROPS: PVTO needs at least one record";
            if (p.pvto.size() < 2 && p.pvto[0].under.empty()) return "PROPS: a single PVTO record needs undersaturated rows";
            for (size_t i = 1; i < p.pvto.size(); ++i)
                if (!(p.pvto[i].rs > p.pvto[i - 1].rs) || !(p.pvto[i].pb > p.pvto[i - 1].pb))
                    return "PROPS: PVTO Rs and bubble-point pressure must increase";
            if (p.pvto.back().under.empty()) return "PROPS: the last PVTO record needs undersaturated rows";
        } else {
            if (p.pvdo.size() < 2) return "PROPS: PVDO needs at least 2 rows";
            std::vector<double> pp; for (auto& r : p.pvdo) pp.push_back(r.p);
            if (!increasing(pp)) return "PROPS: PVDO pressures must increase";
        }
    }
    if (p.pvtw.bw <= 0 || p.pvtw.visc <= 0) return "PROPS: PVTW Bw and viscosity must be positive";
    if (p.density.oil <= 0 || p.density.water <= 0 || p.density.gas <= 0) return "PROPS: densities must be positive";

    if (c.schedule.empty()) return "SCHEDULE: at least one report step is required";
    for (auto& s : c.schedule) if (!(s.length > 0)) return "SCHEDULE: report step lengths must be positive";
    std::set<std::string> names;
    for (auto& w : c.wells) {
        if (w.name.empty()) return "SCHEDULE: a well has no name";
        if (!names.insert(w.name).second) return "SCHEDULE: duplicate well name " + w.name;
        if (w.type == WellType::WaterInjector && !ph.water) return "Well " + w.name + ": water injector but water is not active";
        if (w.type == WellType::GasInjector && !ph.gas) return "Well " + w.name + ": gas injector but gas is not active";
        for (auto& cp : w.completions) {
            if (cp.i < 0 || cp.i >= g.nx || cp.j < 0 || cp.j >= g.ny || cp.k1 < 0 || cp.k2 >= g.nz || cp.k1 > cp.k2)
                return "Well " + w.name + ": completion indices outside the grid";
            if (cp.cf <= 0 && cp.diameter <= 0) return "Well " + w.name + ": wellbore diameter must be positive";
        }
    }
    for (auto& s : c.schedule)
        for (auto& wc : s.controls) {
            bool ok = false;
            for (auto& w : c.wells) {
                const std::string& pat = wc.well;
                if (!pat.empty() && pat.back() == '*') { if (w.name.compare(0, pat.size() - 1, pat, 0, pat.size() - 1) == 0) ok = true; }
                else if (w.name == pat) {
                    ok = true;
                    const bool prod = w.type == WellType::Producer;
                    if (prod && wc.mode == Control::RATE) return "Well " + w.name + ": RATE control is for injectors";
                    if (!prod && wc.mode != Control::RATE && wc.mode != Control::BHP) return "Well " + w.name + ": injectors use RATE or BHP control";
                    if (wc.status == WellStatus::Open) {
                        double t = -1;
                        switch (wc.mode) {
                            case Control::ORAT: t = wc.orat; break; case Control::WRAT: t = wc.wrat; break;
                            case Control::GRAT: t = wc.grat; break; case Control::LRAT: t = wc.lrat; break;
                            case Control::RESV: t = wc.resv; break; case Control::RATE: t = wc.rate; break;
                            case Control::BHP: t = wc.bhp; break;
                        }
                        if (wc.history && wc.mode == Control::RESV) t = 0;   // computed from the observed rates
                        if (t < 0) return "Well " + w.name + ": no target given for control mode " + toString(wc.mode);
                    }
                }
            }
            if (!ok) return "SCHEDULE: control refers to unknown well " + wc.well;
        }
    if (c.phases.polymer) {
        const PolymerInput& P = c.polymer;
        if (!c.phases.water) return "POLYMER: the polymer option needs the water phase";
        if (P.plyvisc.empty() || P.plyads.empty() || P.plyrock.empty()) return "POLYMER: PLYVISC, PLYADS and PLYROCK are required";
        auto tab = [](const std::vector<double>& t, const char* kw, double first) -> std::string {
            if (t.size() < 4 || t.size() % 2) return std::string("POLYMER: ") + kw + " needs at least 2 rows";
            for (size_t i = 2; i < t.size(); i += 2) if (!(t[i] > t[i - 2])) return std::string("POLYMER: ") + kw + " concentrations must increase";
            if (std::fabs(t[0]) > 1e-12) return std::string("POLYMER: ") + kw + " must start at concentration 0";
            if (first >= 0 && std::fabs(t[1] - first) > 1e-9) return std::string("POLYMER: ") + kw + (first == 1 ? " multiplier at c = 0 must be 1" : " adsorption at c = 0 must be 0");
            for (size_t i = 1; i < t.size(); i += 2) if (!(t[i] >= 0) || !std::isfinite(t[i])) return std::string("POLYMER: ") + kw + " has negative or invalid values";
            return "";
        };
        for (auto& t : P.plyvisc) { auto r = tab(t, "PLYVISC", 1.0); if (!r.empty()) return r; }
        for (auto& t : P.plyads) { auto r = tab(t, "PLYADS", 0.0); if (!r.empty()) return r; }
        for (auto& R : P.plyrock) {
            if (!(R.rrf >= 1)) return "POLYMER: PLYROCK residual resistance factor must be >= 1";
            if (!(R.deadPoreVolume >= 0 && R.deadPoreVolume < 1)) return "POLYMER: PLYROCK inaccessible pore volume must be in [0, 1)";
            if (!(R.rockDensity >= 0) || !(R.maxAdsorption >= 0)) return "POLYMER: PLYROCK rock density / maximum adsorption must be >= 0";
        }
        const int nm = P.numMixRegions();
        if (nm < 1) return "POLYMER: PLYMAX / PLMIXPAR need at least one record";
        for (int k = 0; k < nm; ++k) {
            if (!(P.cmax[size_t(k)] > 0)) return "POLYMER: PLYMAX maximum concentration must be > 0";
            if (!(P.mixing[size_t(k)] >= 0 && P.mixing[size_t(k)] <= 1)) return "POLYMER: PLMIXPAR mixing parameter must be in [0, 1]";
        }
        if (!g.plmixnum.empty()) {
            if (g.plmixnum.size() != n) return "GRID: PLMIXNUM must have NX*NY*NZ values";
            for (int v : g.plmixnum) if (v < 1 || v > nm) return "POLYMER: PLMIXNUM refers to a region without PLYMAX / PLMIXPAR data";
        }
        if (!g.spoly.empty() && g.spoly.size() != n) return "GRID: SPOLY must have NX*NY*NZ values";
        for (auto& s : P.plyshlog) {
            if (s.rows.size() < 2 || s.rows.size() % 2) return "POLYMER: PLYSHLOG needs at least 1 row";
            for (size_t i = 2; i < s.rows.size(); i += 2) if (!(s.rows[i] > s.rows[i - 2])) return "POLYMER: PLYSHLOG velocities must increase";
            for (size_t i = 0; i < s.rows.size(); i += 2) if (!(s.rows[i] > 0) || !(s.rows[i + 1] > 0)) return "POLYMER: PLYSHLOG velocities and multipliers must be > 0";
        }
        for (auto& st : c.schedule)
            for (auto& wp : st.wpolymer) {
                if (!(wp.conc >= 0)) return "SCHEDULE: WPOLYMER concentration must be >= 0";
                bool ok = false;
                for (auto& w : c.wells) {
                    const std::string& pat = wp.well;
                    if (!pat.empty() && pat.back() == '*') { if (w.name.compare(0, pat.size() - 1, pat, 0, pat.size() - 1) == 0) ok = true; }
                    else if (w.name == pat) ok = true;
                }
                if (!ok) return "SCHEDULE: WPOLYMER refers to unknown well " + wp.well;
            }
    }
    const Tuning& t = c.tuning;
    if (!(t.initialDt > 0) || !(t.maxDt >= t.minDt) || !(t.minDt > 0)) return "TUNING: invalid time step limits";
    if (!(t.cutFactor > 0 && t.cutFactor < 1) || !(t.growth >= 1)) return "TUNING: invalid growth / cut factors";
    return "";
}

// ---------------------------------------------------------------- demo helpers
void coreySWOF(std::vector<SatRow>& out, double swc, double sor, double nw, double no, double krwMax, double kroMax, int n) {
    out.clear();
    for (int i = 0; i <= n; ++i) {
        const double sw = swc + (1.0 - swc) * i / n;
        const double se = std::max(0.0, std::min(1.0, (sw - swc) / (1.0 - swc - sor)));
        const double so = std::max(0.0, std::min(1.0, (1.0 - sw - sor) / (1.0 - swc - sor)));
        out.push_back({sw, krwMax * std::pow(se, nw), kroMax * std::pow(so, no), 0.0});
    }
}
void coreySGOF(std::vector<SatRow>& out, double swc, double sgc, double sorg, double ng, double nog, double krgMax, int n) {
    out.clear();
    const double sgMax = 1.0 - swc;
    for (int i = 0; i <= n; ++i) {
        const double sg = sgMax * i / n;
        const double sge = std::max(0.0, std::min(1.0, (sg - sgc) / (1.0 - swc - sgc - sorg)));
        const double so = std::max(0.0, std::min(1.0, (1.0 - swc - sg - sorg) / (1.0 - swc - sorg)));
        out.push_back({sg, krgMax * std::pow(sge, ng), std::pow(so, nog), 0.0});
    }
}

Case makeDemoCase() {
    Case c;
    c.title = "Demo: SPE1-type live oil waterflood (20x20x4)";
    c.units = UnitSystem::Field;
    c.phases = Phases{true, true, true, true};
    GridInput& g = c.grid;
    g.nx = 20; g.ny = 20; g.nz = 4;
    const int n = g.size();
    g.dx.assign(n, 300); g.dy.assign(n, 300); g.dz.assign(n, 25);
    g.tops.assign(g.nx * g.ny, 8325);
    const double kL[4] = {200, 50, 500, 100}, phiL[4] = {0.22, 0.18, 0.25, 0.20};
    g.poro.resize(n); g.permx.resize(n); g.permz.resize(n);
    for (int k = 0; k < g.nz; ++k)
        for (int j = 0; j < g.ny; ++j)
            for (int i = 0; i < g.nx; ++i) {
                const int q = g.index(i, j, k);
                // mild deterministic areal variation so maps are not uniform
                const double f = 1.0 + 0.35 * std::sin(0.45 * i + 0.3 * k) * std::cos(0.37 * j - 0.2 * k);
                g.permx[q] = kL[k] * f; g.permz[q] = 0.1 * g.permx[q];
                g.poro[q] = phiL[k] * (0.92 + 0.08 * f);
            }
    g.permy = g.permx;

    PropsInput& p = c.props;
    p.pvtw = {4014.7, 1.029, 3.13e-6, 0.31, 0.0};
    p.rock = {14.7, 3e-6};
    p.density = {49.1, 64.79, 0.06054};
    p.pvto = {
        {0.001, 14.7, 1.062, 1.04, {}},  {0.0905, 264.7, 1.15, 0.975, {}}, {0.18, 514.7, 1.207, 0.91, {}},
        {0.371, 1014.7, 1.295, 0.83, {}}, {0.636, 2014.7, 1.435, 0.695, {}}, {0.775, 2514.7, 1.5, 0.641, {}},
        {0.93, 3014.7, 1.565, 0.594, {}}, {1.27, 4014.7, 1.695, 0.51, {{9014.7, 1.579, 0.74}}},
        {1.618, 5014.7, 1.827, 0.449, {{9014.7, 1.737, 0.631}}}};
    p.pvdg = {{14.7, 166.666, 0.008}, {264.7, 12.093, 0.0096}, {514.7, 6.274, 0.0112}, {1014.7, 3.197, 0.014},
              {2014.7, 1.614, 0.0189}, {2514.7, 1.294, 0.0208}, {3014.7, 1.08, 0.0228}, {4014.7, 0.811, 0.0268},
              {5014.7, 0.649, 0.0309}, {9014.7, 0.386, 0.047}};
    coreySWOF(p.swof, 0.2, 0.2, 2.0, 2.0, 0.6, 1.0);
    coreySGOF(p.sgof, 0.2, 0.05, 0.2, 2.0, 2.0, 0.8);

    c.equil.datumDepth = 8400; c.equil.datumPressure = 4800;
    c.equil.owc = 9000; c.equil.goc = 8000; c.equil.initialRs = 1.27;

    WellInput prod; prod.name = "PROD"; prod.type = WellType::Producer; prod.i = 19; prod.j = 19;
    { Completion cp; cp.i = 19; cp.j = 19; cp.k1 = 0; cp.k2 = 1; prod.completions.push_back(cp); }
    WellInput inj; inj.name = "INJ"; inj.type = WellType::WaterInjector; inj.i = 0; inj.j = 0;
    { Completion cp; cp.i = 0; cp.j = 0; cp.k1 = 2; cp.k2 = 3; inj.completions.push_back(cp); }
    c.wells = {prod, inj};

    WellControl wp; wp.well = "PROD"; wp.mode = Control::ORAT; wp.orat = 10000; wp.bhp = 1000;
    WellControl wi; wi.well = "INJ"; wi.mode = Control::RATE; wi.rate = 12000; wi.bhp = 6000;
    for (int s = 0; s < 40; ++s) {
        ReportStep r; r.length = 91.3125;   // 10 years in 40 quarterly steps
        if (s == 0) r.controls = {wp, wi};
        c.schedule.push_back(r);
    }
    return c;
}

Case makeBuckleyLeverettCase(int n) {
    Case c;
    c.title = "Validation: 1D Buckley-Leverett waterflood";
    c.units = UnitSystem::Field;
    c.phases = Phases{true, true, false, false};
    GridInput& g = c.grid;
    g.nx = n; g.ny = 1; g.nz = 1;
    g.dx.assign(n, 1000.0 / n); g.dy.assign(n, 100); g.dz.assign(n, 10);
    g.tops.assign(n, 5000);
    g.poro.assign(n, 0.2); g.permx.assign(n, 100);
    PropsInput& p = c.props;
    p.pvtw = {1000, 1.0, 1e-7, 1.0, 0.0};
    p.rock = {1000, 1e-7};
    p.density = {50, 50, 0.05};
    p.pvdo = {{500, 1.0, 1.0}, {5000, 0.9999, 1.0}};
    coreySWOF(p.swof, 0.0, 0.0, 2.0, 2.0, 1.0, 1.0, 40);
    c.equil.datumDepth = 5005; c.equil.datumPressure = 1000; c.equil.owc = 10000; c.equil.goc = 0;
    WellInput inj; inj.name = "INJ"; inj.type = WellType::WaterInjector;
    { Completion cp; cp.i = 0; cp.cf = 50; inj.completions.push_back(cp); }
    WellInput prod; prod.name = "PROD"; prod.type = WellType::Producer; prod.i = n - 1;
    { Completion cp; cp.i = n - 1; cp.cf = 50; prod.completions.push_back(cp); }
    c.wells = {inj, prod};
    WellControl wi; wi.well = "INJ"; wi.mode = Control::RATE; wi.rate = 100; wi.bhp = 100000;
    WellControl wp; wp.well = "PROD"; wp.mode = Control::BHP; wp.bhp = 1000;
    for (int s = 0; s < 10; ++s) { ReportStep r; r.length = 36; if (s == 0) r.controls = {wi, wp}; c.schedule.push_back(r); }
    c.tuning.maxDt = 2.0;
    return c;
}


const char* toString(ObsQty q) {
    static const char* n[] = {"WOPR", "WWPR", "WGPR", "WWIR", "WBHP", "WWCT", "WGOR", "FOPR", "FWPR", "FGPR", "FPR", "FWCT", "FGOR"};
    return n[int(q)];
}
bool obsIsWell(ObsQty q) { return int(q) <= int(ObsQty::WGOR); }

} // namespace ressim
