#include <cstdio>
#include <cstdlib>
#include "ressim/io.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iomanip>
#include <cmath>

namespace ressim {

// ------------------------------------------------------------------ helpers
static std::string num(double v) {
    char b[40];
    for (int p = 10; p <= 17; ++p) { std::snprintf(b, sizeof(b), "%.*g", p, v); if (std::strtod(b, nullptr) == v) break; }
    return b;
}
// run-length compressed array: "n*v" ECLIPSE style
template <class T>
static void writeArray(std::ostream& o, const std::vector<T>& v, int perLine = 8) {
    size_t i = 0; int col = 0;
    while (i < v.size()) {
        size_t j = i + 1;
        while (j < v.size() && v[j] == v[i]) ++j;
        const size_t cnt = j - i;
        if (cnt > 1) o << cnt << '*' << num(double(v[i])); else o << num(double(v[i]));
        o << ((++col % perLine == 0) ? "\n" : " ");
        i = j;
    }
    if (col % perLine != 0) o << "\n";
}

static const char* ctrlName(Control c) { return toString(c); }
static Control ctrlFrom(const std::string& s) {
    for (Control c : {Control::ORAT, Control::WRAT, Control::GRAT, Control::LRAT, Control::RESV, Control::BHP, Control::RATE})
        if (s == toString(c)) return c;
    throw std::runtime_error("Unknown control mode " + s);
}
static const char* typeTag(WellType t) {
    return t == WellType::Producer ? "PROD" : t == WellType::WaterInjector ? "WINJ" : "GINJ";
}
static WellType typeFrom(const std::string& s) {
    if (s == "PROD") return WellType::Producer;
    if (s == "WINJ") return WellType::WaterInjector;
    if (s == "GINJ") return WellType::GasInjector;
    throw std::runtime_error("Unknown well type " + s);
}
static const char* statusTag(WellStatus s) { return s == WellStatus::Open ? "OPEN" : s == WellStatus::Shut ? "SHUT" : "STOP"; }
static WellStatus statusFrom(const std::string& s) {
    return s == "OPEN" ? WellStatus::Open : s == "STOP" ? WellStatus::Stop : WellStatus::Shut;
}

// ------------------------------------------------------------------ native format
std::string projectToString(const Case& c) {
    std::ostringstream o;
    o << "RESSIM_PROJECT 1\n";
    o << "TITLE " << c.title << "\n";
    o << "UNITS " << (c.units == UnitSystem::Field ? "FIELD" : "METRIC") << "\n";
    o << "PHASES " << c.phases.water << ' ' << c.phases.oil << ' ' << c.phases.gas << ' ' << c.phases.disgas << "\n";
    const GridInput& g = c.grid;
    o << "DIMENS " << g.nx << ' ' << g.ny << ' ' << g.nz << "\n";
    auto arr = [&](const char* name, const auto& v) {
        if (v.empty()) return;
        o << "ARRAY " << name << ' ' << v.size() << "\n";
        writeArray(o, v);
    };
    arr("DX", g.dx); arr("DY", g.dy); arr("DZ", g.dz); arr("TOPS", g.tops);
    arr("PORO", g.poro); arr("PERMX", g.permx); arr("PERMY", g.permy); arr("PERMZ", g.permz); arr("NTG", g.ntg);
    arr("ACTNUM", g.actnum); arr("MULTX", g.multx); arr("MULTY", g.multy); arr("MULTZ", g.multz);
    arr("COORD", g.coord); arr("ZCORN", g.zcorn); arr("SWATINIT", g.swatinit);
    for (auto& t : c.thpres) o << "THPRES " << t.r1 << ' ' << t.r2 << ' ' << num(t.value) << "\n";
    o << "START " << c.startDay << ' ' << c.startMonth << ' ' << c.startYear << "\n";
    for (auto& f : g.faults) {
        o << "FAULT " << f.name << ' ' << num(f.mult) << ' ' << f.faces.size() << "\n";
        for (auto& ff : f.faces)
            o << "  FACE " << ff.i1 << ' ' << ff.i2 << ' ' << ff.j1 << ' ' << ff.j2 << ' ' << ff.k1 << ' ' << ff.k2 << ' ' << ff.face << ' ' << ff.minus << "\n";
    }
    for (auto& e : g.nnc)
        o << "NNC " << e.i1 << ' ' << e.j1 << ' ' << e.k1 << ' ' << e.i2 << ' ' << e.j2 << ' ' << e.k2 << ' ' << num(e.trans) << "\n";
    for (auto& a : c.aquifers) {
        o << "AQUIFER " << a.id << ' ' << (a.type == AquiferType::Fetkovich ? "FETKOVICH" : "CARTERTRACY") << ' ' << num(a.depth) << ' '
          << num(a.p0) << ' ' << num(a.ct) << ' ' << num(a.volume) << ' ' << num(a.pi) << ' ' << num(a.perm) << ' ' << num(a.poro) << ' '
          << num(a.rInner) << ' ' << num(a.thickness) << ' ' << num(a.angle) << ' ' << num(a.visc) << ' ' << a.conns.size() << "\n";
        for (auto& q : a.conns)
            o << "  ACONN " << q.i1 << ' ' << q.i2 << ' ' << q.j1 << ' ' << q.j2 << ' ' << q.k1 << ' ' << q.k2 << ' ' << q.face << ' ' << q.minus
              << ' ' << num(q.mult) << "\n";
    }
    const PropsInput& p = c.props;
    o << "PVTW " << num(p.pvtw.pref) << ' ' << num(p.pvtw.bw) << ' ' << num(p.pvtw.cw) << ' ' << num(p.pvtw.visc) << ' ' << num(p.pvtw.cv) << "\n";
    o << "ROCK " << num(p.rock.pref) << ' ' << num(p.rock.cr) << "\n";
    o << "DENSITY " << num(p.density.oil) << ' ' << num(p.density.water) << ' ' << num(p.density.gas) << "\n";
    auto pvd = [&](const char* name, const std::vector<PVDRow>& t) {
        o << name << ' ' << t.size() << "\n";
        for (auto& r : t) o << "  " << num(r.p) << ' ' << num(r.B) << ' ' << num(r.mu) << "\n";
    };
    pvd("PVDO", p.pvdo); pvd("PVDG", p.pvdg);
    o << "PVTO " << p.pvto.size() << "\n";
    for (auto& r : p.pvto) {
        o << "  " << num(r.rs) << ' ' << num(r.pb) << ' ' << num(r.bo) << ' ' << num(r.mu) << ' ' << r.under.size();
        for (auto& u : r.under) o << ' ' << num(u.p) << ' ' << num(u.B) << ' ' << num(u.mu);
        o << "\n";
    }
    auto sat = [&](const char* name, const std::vector<SatRow>& t) {
        o << name << ' ' << t.size() << "\n";
        for (auto& r : t) o << "  " << num(r.s) << ' ' << num(r.krw) << ' ' << num(r.kro) << ' ' << num(r.pc) << "\n";
    };
    sat("SWOF", p.swof); sat("SGOF", p.sgof);
    const EquilInput& e = c.equil;
    o << "EQUIL " << num(e.datumDepth) << ' ' << num(e.datumPressure) << ' ' << num(e.owc) << ' ' << num(e.pcowOWC) << ' '
      << num(e.goc) << ' ' << num(e.pcogGOC) << ' ' << e.nacc << ' ' << num(e.initialRs) << "\n";
    for (auto& w : c.wells) {
        o << "WELL " << w.name << ' ' << w.group << ' ' << typeTag(w.type) << ' ' << w.i << ' ' << w.j << ' ' << num(w.refDepth)
          << ' ' << w.completions.size() << "\n";
        for (auto& cp : w.completions)
            o << "  COMP " << cp.i << ' ' << cp.j << ' ' << cp.k1 << ' ' << cp.k2 << ' ' << num(cp.cf) << ' ' << num(cp.kh) << ' '
              << num(cp.skin) << ' ' << num(cp.diameter) << ' ' << cp.dir << ' ' << cp.open << "\n";
    }
    for (auto& w : c.wells) if (w.piMult != 1.0) o << "WPIMULT " << w.name << ' ' << num(w.piMult) << "\n";
    for (auto& ob : c.observations)
        o << "OBS " << toString(ob.qty) << ' ' << (ob.well.empty() ? std::string("-") : ob.well) << ' ' << num(ob.time) << ' ' << num(ob.value)
          << ' ' << num(ob.sigma) << "\n";
    for (auto& p : c.studyParams)
        o << "PARAM " << int(p.kind) << ' ' << p.target << ' ' << p.i1 << ' ' << p.i2 << ' ' << p.j1 << ' ' << p.j2 << ' ' << p.k1 << ' ' << p.k2
          << ' ' << num(p.base) << ' ' << num(p.low) << ' ' << num(p.high) << ' ' << int(p.log) << ' ' << int(p.active) << ' ' << p.label << "\n";
    for (auto& s : c.schedule) {
        o << "STEP " << num(s.length) << ' ' << s.controls.size() << "\n";
        for (auto& wc : s.controls)
            o << "  CTRL " << wc.well << ' ' << statusTag(wc.status) << ' ' << ctrlName(wc.mode) << ' ' << num(wc.orat) << ' '
              << num(wc.wrat) << ' ' << num(wc.grat) << ' ' << num(wc.lrat) << ' ' << num(wc.resv) << ' ' << num(wc.rate) << ' '
              << num(wc.bhp) << ' ' << wc.history << "\n";
    }
    const Tuning& t = c.tuning;
    o << "TUNING " << num(t.initialDt) << ' ' << num(t.maxDt) << ' ' << num(t.minDt) << ' ' << num(t.growth) << ' '
      << num(t.cutFactor) << ' ' << t.maxNewton << ' ' << t.maxLinear << ' ' << num(t.cnvTol) << ' ' << num(t.mbTol) << ' '
      << num(t.linTol) << ' ' << num(t.dsMax) << ' ' << num(t.dpMaxRel) << ' ' << t.threads << "\n";
    o << "SOLVER " << t.linearSolver << ' ' << t.cprWeights << ' ' << (t.mixedPrecision ? 1 : 0) << ' ' << t.gpu << ' ' << t.leanMemory << "\n";
    o << "END\n";
    return o.str();
}

namespace {
struct Tok {
    std::istringstream in;
    explicit Tok(const std::string& s) : in(s) {}
    std::string word() { std::string w; if (!(in >> w)) throw std::runtime_error("Unexpected end of project file"); return w; }
    double d() {
        const std::string w = word();
        try { size_t pos; double v = std::stod(w, &pos); if (pos != w.size()) throw 0; return v; }
        catch (...) { throw std::runtime_error("Expected a number but found '" + w + "'"); }
    }
    int i() { return int(std::lround(d())); }
    std::string rest() { std::string l; std::getline(in, l); size_t a = l.find_first_not_of(" \t"); return a == std::string::npos ? "" : l.substr(a); }
};
template <class T>
void readArray(Tok& t, std::vector<T>& v, size_t n) {
    v.clear(); v.reserve(n);
    while (v.size() < n) {
        const std::string w = t.word();
        const size_t star = w.find('*');
        if (star == std::string::npos) v.push_back(T(std::stod(w)));
        else {
            const size_t cnt = std::stoul(w.substr(0, star));
            const T val = T(std::stod(w.substr(star + 1)));
            for (size_t k = 0; k < cnt; ++k) v.push_back(val);
        }
    }
    if (v.size() != n) throw std::runtime_error("Array length mismatch in project file");
}
}

Case projectFromString(const std::string& text) {
    Tok t(text);
    if (t.word() != "RESSIM_PROJECT") throw std::runtime_error("Not a ResSim project file");
    t.d();
    Case c;
    c.wells.clear(); c.schedule.clear();
    c.props.pvdo.clear(); c.props.pvdg.clear(); c.props.pvto.clear(); c.props.swof.clear(); c.props.sgof.clear();
    for (;;) {
        const std::string k = t.word();
        if (k == "END") break;
        if (k == "TITLE") c.title = t.rest();
        else if (k == "UNITS") c.units = t.word() == "FIELD" ? UnitSystem::Field : UnitSystem::Metric;
        else if (k == "PHASES") { c.phases.water = t.i(); c.phases.oil = t.i(); c.phases.gas = t.i(); c.phases.disgas = t.i(); }
        else if (k == "DIMENS") { c.grid.nx = t.i(); c.grid.ny = t.i(); c.grid.nz = t.i(); }
        else if (k == "ARRAY") {
            const std::string name = t.word(); const size_t n = size_t(t.d());
            GridInput& g = c.grid;
            if (name == "ACTNUM") readArray(t, g.actnum, n);
            else {
                std::vector<double>* dst = name == "DX" ? &g.dx : name == "DY" ? &g.dy : name == "DZ" ? &g.dz : name == "TOPS" ? &g.tops :
                                           name == "PORO" ? &g.poro : name == "PERMX" ? &g.permx : name == "PERMY" ? &g.permy :
                                           name == "PERMZ" ? &g.permz : name == "NTG" ? &g.ntg : name == "MULTX" ? &g.multx :
                                           name == "MULTY" ? &g.multy : name == "MULTZ" ? &g.multz : name == "COORD" ? &g.coord :
                                           name == "ZCORN" ? &g.zcorn : name == "SWATINIT" ? &g.swatinit : nullptr;
                if (!dst) throw std::runtime_error("Unknown array " + name);
                readArray(t, *dst, n);
            }
        }
        else if (k == "THPRES") { ThresholdPressure x; x.r1 = t.i(); x.r2 = t.i(); x.value = t.d(); c.thpres.push_back(x); }
        else if (k == "SOLVER") { c.tuning.linearSolver = t.i(); const std::string r = t.rest(); if (!r.empty()) { int a = 1, m = 0, g = 0, lm = -1; const int nr = std::sscanf(r.c_str(), "%d %d %d %d", &a, &m, &g, &lm);
                if (nr >= 1) { c.tuning.cprWeights = a ? 1 : 0; } if (nr >= 2) { c.tuning.mixedPrecision = m != 0; } if (nr >= 3) { c.tuning.gpu = g; }
                if (nr >= 4) { c.tuning.leanMemory = std::clamp(lm, -1, 1); } } }
        else if (k == "WPIMULT") { const std::string n = t.word(); const double m = t.d(); for (auto& w : c.wells) if (w.name == n) w.piMult = m; }
        else if (k == "PARAM") {
            StudyParam p; p.kind = ParamKind(t.i()); p.target = t.i();
            p.i1 = t.i(); p.i2 = t.i(); p.j1 = t.i(); p.j2 = t.i(); p.k1 = t.i(); p.k2 = t.i();
            p.base = t.d(); p.low = t.d(); p.high = t.d(); p.log = t.i() != 0; p.active = t.i() != 0; p.label = t.rest();
            if (int(p.kind) < 0 || int(p.kind) >= kNumParamKinds) throw std::runtime_error("Bad PARAM kind");
            c.studyParams.push_back(p);
            continue;
        }
        else if (k == "OBS") {
            Observation ob; const std::string q = t.word(); bool f = false;
            for (int i = 0; i <= int(ObsQty::FGOR); ++i) if (q == toString(ObsQty(i))) { ob.qty = ObsQty(i); f = true; }
            if (!f) throw std::runtime_error("Unknown observation quantity " + q);
            ob.well = t.word(); if (ob.well == "-") ob.well.clear();
            ob.time = t.d(); ob.value = t.d(); ob.sigma = t.d();
            c.observations.push_back(ob);
        }
        else if (k == "START") { c.startDay = t.i(); c.startMonth = t.i(); c.startYear = t.i(); }
        else if (k == "FAULT") {
            FaultInput f; f.name = t.word(); f.mult = t.d(); const int nf = t.i();
            for (int q = 0; q < nf; ++q) {
                if (t.word() != "FACE") throw std::runtime_error("Expected FACE in fault " + f.name);
                FaultFace ff; ff.i1 = t.i(); ff.i2 = t.i(); ff.j1 = t.i(); ff.j2 = t.i(); ff.k1 = t.i(); ff.k2 = t.i();
                ff.face = t.word()[0]; ff.minus = t.i() != 0; f.faces.push_back(ff);
            }
            c.grid.faults.push_back(f);
        }
        else if (k == "NNC") { NNCInput e; e.i1 = t.i(); e.j1 = t.i(); e.k1 = t.i(); e.i2 = t.i(); e.j2 = t.i(); e.k2 = t.i(); e.trans = t.d(); c.grid.nnc.push_back(e); }
        else if (k == "AQUIFER") {
            AquiferInput a; a.id = t.i(); a.type = t.word() == "FETKOVICH" ? AquiferType::Fetkovich : AquiferType::CarterTracy;
            a.depth = t.d(); a.p0 = t.d(); a.ct = t.d(); a.volume = t.d(); a.pi = t.d(); a.perm = t.d(); a.poro = t.d();
            a.rInner = t.d(); a.thickness = t.d(); a.angle = t.d(); a.visc = t.d();
            const int nc = t.i();
            for (int q = 0; q < nc; ++q) {
                if (t.word() != "ACONN") throw std::runtime_error("Expected ACONN in aquifer");
                AquiferConnection ac; ac.i1 = t.i(); ac.i2 = t.i(); ac.j1 = t.i(); ac.j2 = t.i(); ac.k1 = t.i(); ac.k2 = t.i();
                ac.face = t.word()[0]; ac.minus = t.i() != 0; ac.mult = t.d(); a.conns.push_back(ac);
            }
            c.aquifers.push_back(a);
        }
        else if (k == "PVTW") { auto& p = c.props.pvtw; p.pref = t.d(); p.bw = t.d(); p.cw = t.d(); p.visc = t.d(); p.cv = t.d(); }
        else if (k == "ROCK") { c.props.rock.pref = t.d(); c.props.rock.cr = t.d(); }
        else if (k == "DENSITY") { auto& d = c.props.density; d.oil = t.d(); d.water = t.d(); d.gas = t.d(); }
        else if (k == "PVDO" || k == "PVDG") {
            auto& tab = k == "PVDO" ? c.props.pvdo : c.props.pvdg;
            const int n = t.i(); tab.clear();
            for (int r = 0; r < n; ++r) { PVDRow row; row.p = t.d(); row.B = t.d(); row.mu = t.d(); tab.push_back(row); }
        }
        else if (k == "PVTO") {
            const int n = t.i();
            for (int r = 0; r < n; ++r) {
                PVTORecord rec; rec.rs = t.d(); rec.pb = t.d(); rec.bo = t.d(); rec.mu = t.d();
                const int nu = t.i();
                for (int u = 0; u < nu; ++u) { PVDRow row; row.p = t.d(); row.B = t.d(); row.mu = t.d(); rec.under.push_back(row); }
                c.props.pvto.push_back(rec);
            }
        }
        else if (k == "SWOF" || k == "SGOF") {
            auto& tab = k == "SWOF" ? c.props.swof : c.props.sgof;
            const int n = t.i();
            for (int r = 0; r < n; ++r) { SatRow row; row.s = t.d(); row.krw = t.d(); row.kro = t.d(); row.pc = t.d(); tab.push_back(row); }
        }
        else if (k == "EQUIL") {
            auto& e = c.equil; e.datumDepth = t.d(); e.datumPressure = t.d(); e.owc = t.d(); e.pcowOWC = t.d();
            e.goc = t.d(); e.pcogGOC = t.d(); e.nacc = t.i(); e.initialRs = t.d();
        }
        else if (k == "WELL") {
            WellInput w; w.name = t.word(); w.group = t.word(); w.type = typeFrom(t.word()); w.i = t.i(); w.j = t.i();
            w.refDepth = t.d();
            const int nc = t.i();
            for (int q = 0; q < nc; ++q) {
                if (t.word() != "COMP") throw std::runtime_error("Expected COMP in well " + w.name);
                Completion cp; cp.i = t.i(); cp.j = t.i(); cp.k1 = t.i(); cp.k2 = t.i(); cp.cf = t.d(); cp.kh = t.d();
                cp.skin = t.d(); cp.diameter = t.d(); cp.dir = t.word()[0]; cp.open = t.i() != 0;
                w.completions.push_back(cp);
            }
            c.wells.push_back(w);
        }
        else if (k == "STEP") {
            ReportStep s; s.length = t.d();
            const int nc = t.i();
            for (int q = 0; q < nc; ++q) {
                if (t.word() != "CTRL") throw std::runtime_error("Expected CTRL in schedule");
                WellControl wc; wc.well = t.word(); wc.status = statusFrom(t.word()); wc.mode = ctrlFrom(t.word());
                wc.orat = t.d(); wc.wrat = t.d(); wc.grat = t.d(); wc.lrat = t.d(); wc.resv = t.d(); wc.rate = t.d();
                wc.bhp = t.d(); wc.history = t.i() != 0;
                s.controls.push_back(wc);
            }
            c.schedule.push_back(s);
        }
        else if (k == "TUNING") {
            Tuning& u = c.tuning; u.initialDt = t.d(); u.maxDt = t.d(); u.minDt = t.d(); u.growth = t.d(); u.cutFactor = t.d();
            u.maxNewton = t.i(); u.maxLinear = t.i(); u.cnvTol = t.d(); u.mbTol = t.d(); u.linTol = t.d(); u.dsMax = t.d();
            u.dpMaxRel = t.d(); u.threads = t.i();
        }
        else throw std::runtime_error("Unknown keyword '" + k + "' in project file");
    }
    return c;
}

// Project format 2 (.rsim): a complete ECLIPSE deck with "-- RESSIM:" lines (study parameters, observations, solver
// settings). The deck round trip is exact, so regions, RSVD tables, enumerated initial state, MULTPV/MINPV, NEWTRAN and
// per-control injection phases are all preserved, and the file opens in OPM Flow / ECLIPSE as is.
// Format 1 (RESSIM_PROJECT 1 ...) is still read.
void saveProject(const Case& c, const std::string& path) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot write " + path);
    f << "-- RESSIM_PROJECT 2\n" << eclipseDeckToString(c);
    if (!f) throw std::runtime_error("Error while writing " + path);
}

std::string projectText(const Case& c) { return "-- RESSIM_PROJECT 2\n" + eclipseDeckToString(c); }
Case projectFromText(const std::string& text, const std::string& baseDir) {
    const size_t a = text.find_first_not_of(" \t\r\n");
    if (a != std::string::npos && text.compare(a, 14, "RESSIM_PROJECT") == 0) return projectFromString(text);
    return eclipseDeckFromString(text, baseDir);
}

Case loadProject(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot open " + path);
    std::stringstream ss; ss << f.rdbuf();
    const std::string text = ss.str();
    size_t a = text.find_first_not_of(" \t\r\n");
    if (a != std::string::npos && text.compare(a, 14, "RESSIM_PROJECT") == 0) return projectFromString(text);
    const size_t sl = path.find_last_of("/\\");
    return eclipseDeckFromString(text, sl == std::string::npos ? std::string() : path.substr(0, sl));
}

} // namespace ressim
