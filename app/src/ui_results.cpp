// 3D viewer panel, result plots and run monitor.
#include "app.h"
#include "icons.h"
#include "widgets.h"
#include "i18n.h"
#include "implot.h"
#include "plotkit.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <chrono>
#include <functional>

using namespace ressim;

// ----------------------------------------------------------------------------- property catalogue
namespace {
enum PropKind { Static, Initial, Dynamic };
struct PropDef { const char* name; const char* label; PropKind kind; };
const PropDef kProps[] = {
    {"PORO", "Porosity", Static},          {"PERMX", "Permeability X", Static}, {"PERMZ", "Permeability Z", Static},
    {"NTG", "Net-to-gross", Static},       {"DEPTH", "Cell depth", Static},     {"PORV", "Pore volume", Static},
    {"PRESSURE_0", "Initial pressure", Initial}, {"SWAT_0", "Initial water saturation", Initial},
    {"SGAS_0", "Initial gas saturation", Initial},
    {"PRESSURE", "Pressure", Dynamic},     {"SWAT", "Water saturation", Dynamic}, {"SOIL", "Oil saturation", Dynamic},
    {"SGAS", "Gas saturation", Dynamic},   {"RS", "Solution GOR", Dynamic},
    {"RV", "Vaporised oil-gas ratio", Dynamic},
    {"POLYMER", "Polymer concentration", Dynamic}, {"PADS", "Adsorbed polymer", Dynamic},
    {"TOF", "Time of flight (streamlines)", Dynamic}, {"DRAINAGE", "Drainage region (producer)", Dynamic},
    {"SWEEP", "Sweep region (injector)", Dynamic},
    {"XCO2", "CO2 mole fraction (brine / overall)", Dynamic},
    // appended (C5) - index 9 = PRESSURE is used as the default dynamic property elsewhere
    {"ZMF", "Overall mole fraction (component)", Dynamic}, {"XMF", "Liquid / brine mole fraction (component)", Dynamic},
    {"YMF", "Vapour / CO2-phase mole fraction (component)", Dynamic},
    {"SGTRAP", "Trapped gas saturation (hysteresis)", Dynamic}, {"SGMAX", "Historical maximum gas saturation", Dynamic},
    {"SSALT", "Precipitated salt (solid volume / PV0)", Dynamic}, {"MINCO2", "Mineral-trapped CO2", Dynamic},
    {"PORFAC", "Porosity factor phi / phi0", Dynamic}, {"PERMFAC", "Permeability factor k / k0", Dynamic},
    // appended (C6): aqueous-phase composition (E300 CO2SOL), any further engine restart array (generic)
    {"AMF", "Aqueous mole fraction (component)", Dynamic},
    {"ARRAY", "Engine output array", Dynamic},
    // appended (C13): adjoint per-cell history-matching maps (listed only after 'Show in the 3D viewer')
    {"HMMULTK", "HM permeability multiplier (adjoint)", Static}, {"HMMULTPHI", "HM porosity multiplier (adjoint)", Static},
    {"HMSENSK", "HM sensitivity d misfit / d ln k", Static},
};
bool isHmProp(const std::string& n) { return n == "HMMULTK" || n == "HMMULTPHI" || n == "HMSENSK"; }
std::string gArraySel;   // name of the generic engine array shown with "ARRAY"
int gCompSel = 0;   // component of the ZMF / XMF / YMF maps
bool isCompProp(const std::string& n) { return n == "ZMF" || n == "XMF" || n == "YMF" || n == "AMF"; }
std::string propKey(const char* name) {
    const std::string n = name;
    if (n == "ARRAY") return gArraySel;
    return isCompProp(n) ? n + std::to_string(gCompSel + 1) : n;
}
bool isExtraProp(const std::string& n) { return isCompProp(n) || n == "SGTRAP" || n == "SGMAX" || n == "SSALT" || n == "KSALT" || n == "MINCO2" || n == "PORFAC" || n == "PERMFAC" || n == "ARRAY"; }
std::string baseName(std::string b) { while (!b.empty() && std::isdigit(static_cast<unsigned char>(b.back()))) b.pop_back(); return b; }
// units of generic ECLIPSE / OPM restart mnemonics
const char* arrayUnit(const std::string& n0, UnitSystem us) {
    const std::string n = baseName(n0);
    const bool fld = us == UnitSystem::Field;
    if (n == "XMF" || n == "YMF" || n == "ZMF" || n == "AMF" || n == "VMF" || n == "XCO2") return "mol/mol";
    if (n == "PRESSURE" || n == "PSAT" || n == "POIL" || n == "PGAS" || n == "PWAT" || n == "PBUB" || n == "PDEW" || n == "PCOW" || n == "PCOG") return fld ? "psia" : "bar";
    if (n == "DENO" || n == "DENG" || n == "DENW" || n == "OIL_DEN" || n == "GAS_DEN" || n == "WAT_DEN") return fld ? "lb/ft3" : "kg/m3";
    if (n == "VOIL" || n == "VGAS" || n == "VWAT" || n == "OIL_VISC" || n == "GAS_VISC" || n == "WAT_VISC") return "cP";
    if (n == "RS" || n == "RSW" || n == "RSSAT") return fld ? "Mscf/stb" : "sm3/sm3";
    if (n == "RV" || n == "RVW" || n == "RVSAT") return fld ? "stb/Mscf" : "sm3/sm3";
    if (n == "TEMP") return fld ? "degF" : "degC";
    // C30 fracture / phase-field / geomechanics field arrays
    if (n == "HFPNET" || n.rfind("HFDS", 0) == 0 || n == "HFDTAU" || n.rfind("PFDS", 0) == 0 || n.rfind("SIGEFF", 0) == 0 || n == "SIGMIN" || n == "SQRTJ2") return fld ? "psi" : "bar";
    if (n == "HFWIDTH" || n == "PFWIDTH" || n.rfind("HFDISP", 0) == 0 || n == "PFDISP") return "mm";
    if (n == "UPLIFT" || n == "DISPX" || n == "DISPY" || n == "DISPMAG") return "m";
    if (n == "HFSTRV") return "microstrain";
    if (n == "SOIL" || n == "SWAT" || n == "SGAS" || n.rfind("KR", 0) == 0 || n == "SGTRH" || n == "SGTRAP" || n == "SGMAX") return "";
    return "";
}
bool isStreamlineProp(const char* n) { const std::string s = n; return s == "TOF" || s == "DRAINAGE" || s == "SWEEP"; }
const int kNumProps = int(sizeof(kProps) / sizeof(kProps[0]));

const char* unitOf(const char* name, UnitSystem us) {
    const Units U = Units::of(us);
    const std::string n = name;
    if (n.rfind("PERM", 0) == 0) return "mD";
    if (n == "DEPTH") return U.len;
    if (n == "PORV") return us == UnitSystem::Field ? "rb" : "rm3";
    if (n.rfind("PRESSURE", 0) == 0) return U.p;
    if (n == "RS") return U.rs;
    if (n == "RV") return U.rv;
    if (n == "POLYMER") return us == UnitSystem::Field ? "lb/stb" : "kg/sm3";
    if (n == "PADS") return us == UnitSystem::Field ? "lb/lb" : "kg/kg";
    if (n == "TOF") return "days";
    if (n == "XCO2" || isCompProp(n)) return "mol/mol";
    if (n == "ARRAY") return arrayUnit(gArraySel, us);
    if (n == "MINCO2") return "kg/m3 PV";
    return "";
}
bool isSaturation(const char* name) { const std::string n = name; return n.find("SWAT") == 0 || n.find("SGAS") == 0 || n == "SOIL"; }

// summary vectors
enum Vec { FOPR, FWPR, FGPR, FWIR, FGIR, FOPT, FWPT, FGPT, FWIT, FGIT, FPR, FWCT, FGOR, FOIP, FWIP, FGIP,
           FCIR, FCPR, FCIT, FCPT, FCIP, FCAD, FGMMO, FGMTR, FGMDS, FGMIP, FGMST, FGMMN, FGMIM, FSALT,
           FAVGT, FMAXT, FSTPV, FHLT, FHIN, FCSOR, FRHR, FRHT, FUPL, FUPR, FMXF, DT, NEWTON, LINEAR,
           WOPR, WWPR, WGPR, WWIR, WGIR, WBHP, WWCT, WGOR, WCIR, WCPR, WOPRH, WWPRH, WGPRH, NumVec };
[[maybe_unused]] const char* kVecNames[NumVec] = {"FOPR", "FWPR", "FGPR", "FWIR", "FGIR", "FOPT", "FWPT", "FGPT", "FWIT", "FGIT", "FPR", "FWCT",
                                 "FGOR", "FOIP", "FWIP", "FGIP", "FCIR", "FCPR", "FCIT", "FCPT", "FCIP", "FCAD",
                                 "FGMMO", "FGMTR", "FGMDS", "FGMIP", "FGMST", "FGMMN", "FGMIM", "FSALT",
                                 "FAVGT", "FMAXT", "FSTPV", "FHLT", "FHIN", "FCSOR", "FRHR", "FRHT", "FUPL", "FUPR", "FMXF",
                                 "DT", "NEWTON", "LINEAR", "WOPR", "WWPR", "WGPR", "WWIR",
                                 "WGIR", "WBHP", "WWCT", "WGOR", "WCIR", "WCPR", "WOPRH", "WWPRH", "WGPRH"};
double vecValue(const SummaryPoint& s, int v, int well) {
    const WellRates* w = (well >= 0 && well < int(s.wells.size())) ? &s.wells[well] : nullptr;
    switch (v) {
        case FOPR: return s.fopr; case FWPR: return s.fwpr; case FGPR: return s.fgpr; case FWIR: return s.fwir; case FGIR: return s.fgir;
        case FOPT: return s.fopt; case FWPT: return s.fwpt; case FGPT: return s.fgpt; case FWIT: return s.fwit; case FGIT: return s.fgit;
        case FPR: return s.fpr; case FWCT: return s.fwct; case FGOR: return s.fgor;
        case FOIP: return s.oip; case FWIP: return s.wip; case FGIP: return s.gip;
        case FCIR: return s.fcir; case FCPR: return s.fcpr; case FCIT: return s.fcit; case FCPT: return s.fcpt;
        case FCIP: return s.fcip; case FCAD: return s.fcad;
        case FGMMO: return s.co2Mob; case FGMTR: return s.co2Trap; case FGMDS: return s.co2Diss; case FGMIP: return s.co2Mob + s.co2Trap + s.co2Diss + s.co2Min + s.co2Oil;
        case FGMST: return s.co2Struct; case FGMMN: return s.co2Min; case FGMIM: return s.co2KrImmob; case FSALT: return s.saltSolid;
        case FAVGT: return s.favgT; case FMAXT: return s.fmaxT; case FSTPV: return s.fsteamPV; case FHLT: return s.fhlt;
        case FHIN: return s.fenergyIn; case FCSOR: return s.csor;
        case FRHR: return s.freactRate * 1e-6; case FRHT: return s.freactCum * 1e-9; case FUPL: return s.fuplift * 1e3; case FUPR: return s.fupliftRes * 1e3;
        case FMXF: return s.fmaxFail;
        case DT: return s.dt; case NEWTON: return s.newton; case LINEAR: return s.linear;
        default: break;
    }
    if (!w) return NAN;
    switch (v) {
        case WOPR: return w->oil; case WWPR: return w->water; case WGPR: return w->gas; case WWIR: return w->waterInj;
        case WGIR: return w->gasInj; case WBHP: return w->open ? w->bhp : NAN;
        case WWCT: { const double l = w->oil + w->water; return l > 0 ? w->water / l : 0; }
        case WGOR: return w->oil > 0 ? w->gas / w->oil : 0;
        case WCIR: return w->polyInj; case WCPR: return w->polyProd;
        case WOPRH: return w->oilH >= 0 ? w->oilH : NAN;
        case WWPRH: return w->waterH >= 0 ? w->waterH : NAN;
        case WGPRH: return w->gasH >= 0 ? w->gasH : NAN;
    }
    return NAN;
}
struct Series { const std::vector<SummaryPoint>* s; int vec; int well; double xs; int skipFirst; };
ImPlotPoint seriesGetter(int idx, void* data) {
    const Series* d = static_cast<const Series*>(data);
    const SummaryPoint& p = (*d->s)[idx + d->skipFirst];
    return ImPlotPoint(p.time * d->xs, vecValue(p, d->vec, d->well));
}
void plotVec(const char* label, const std::vector<SummaryPoint>& s, int vec, int well, double xs, ImU32 col, float weight = 2.0f,
             bool scatter = false) {
    // rates are step averages: skip the t=0 point for rate/ratio vectors
    const bool isRate = (vec <= FGIR) || vec == FCIR || vec == FCPR || vec == FWCT || vec == FGOR || (vec >= WOPR && vec != WBHP) || vec == DT || vec == NEWTON || vec == LINEAR;
    Series d{&s, vec, well, xs, (isRate && s.size() > 1) ? 1 : 0};
    const int n = int(s.size()) - d.skipFirst;
    if (n <= 0) return;
    ImPlotSpec spec;
    spec.LineColor = col4(col);
    spec.LineWeight = weight * g_plotLineScale;
    if (fig::recorder().recording()) {
        std::vector<double> X(static_cast<size_t>(n)), Y(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) { const ImPlotPoint q = seriesGetter(i, &d); X[size_t(i)] = q.x; Y[size_t(i)] = q.y; }
        recSeries(label, X.data(), Y.data(), n, col, weight, scatter);
    }
    if (scatter) {
        spec.Marker = ImPlotMarker_Circle; spec.MarkerSize = 3.0f;
        spec.MarkerFillColor = col4(col); spec.MarkerLineColor = col4(col);
        ImPlot::PlotScatterG(label, seriesGetter, &d, n, spec);
    } else ImPlot::PlotLineG(label, seriesGetter, &d, n, spec);
}
ImU32 wellColor(int i) { return seriesColor(i); }
}
int propIndexByName(const char* name) {
    for (int i = 0; i < int(sizeof(kProps) / sizeof(kProps[0])); ++i) if (std::string(kProps[i].name) == name) return i;
    return -1;
}

std::vector<LineVertex> streamlineVertices(const StreamlineResult& S, int colorBy, Colormap cm, const std::string& focus, bool focusOnly);

// ----------------------------------------------------------------------------- view data
// C13: EDFM fracture planes in the Model viewer - run segments coloured by the fracture value of the current dynamic
// property (TEMP / PRESSURE / FRACAP / ...), else the input planes
void App::updateFracturePanels() {
    if (!viewOk_) return;
    const PropDef& pd = kProps[std::min(prop_, kNumProps - 1)];
    const std::string n = pd.kind == Dynamic ? propKey(pd.name) : std::string();
    std::lock_guard<std::mutex> lock(runner_.mutex());
    const RunResults& R = runner_.results();
    const bool useRun = pd.kind == Dynamic && R.valid && ((!R.fracSegs.empty() && !R.fracSnap.empty()) || (!R.input.grid.hfracs.empty() && !R.snaps.empty()));
    const Case& c = useRun ? R.input : case_;
    // C24: fault faces / seismic events of the coupled fault mechanics (results of the run, else the model's fault objects)
    const bool fmRun = R.valid && R.faultMech && !R.snaps.empty() && R.grid.get() == static_cast<const Grid*>(viewGrid_);
    char key[384];
    std::snprintf(key, sizeof(key), "%d|%u|%u|%d|%s|%p|%zu|%g|%g|%d|%d|%d|%g|%zu", useRun ? 1 : 0, caseVersion_, (useRun || fmRun) ? runner_.version() : 0u, (useRun || fmRun) ? snap_ : -1,
                  n.c_str(), viewGrid_, c.thermal.fractures.size() + 1000 * c.grid.hfracs.size(), double(key_.a), double(key_.b),
                  fmRun ? 1 : 0, faultProp_, faultEvents_ ? 1 : 0, double(faultEventScale_), case_.faultMech.faults.size() + (case_.faultMech.on() ? 1000 : 0));
    if (fracKeyMv_ == key) return;
    fracKeyMv_ = key;
    std::vector<PanelQuad> P = fracturePanels(c, static_cast<const Grid*>(viewGrid_), useRun ? &R : nullptr, snap_, n);
    std::vector<PanelQuad> F = faultPanels(fmRun ? R.input : case_, static_cast<const Grid*>(viewGrid_), fmRun ? &R : nullptr, snap_);
    P.insert(P.end(), std::make_move_iterator(F.begin()), std::make_move_iterator(F.end()));
    view_.setPanels(std::move(P));
}

void App::updateViewData() {
    if (!viewOk_) return;
    const PropDef& pd = kProps[std::min(prop_, kNumProps - 1)];
    // choose grid source
    std::shared_ptr<const Grid> results;
    std::vector<float> vals;
    const Grid* grid = nullptr;
    const Case* src = &case_;
    unsigned runVer = runner_.version();
    std::string caption;
    {
        std::lock_guard<std::mutex> lock(runner_.mutex());
        const RunResults& R = runner_.results();
        if (pd.kind == Dynamic && R.valid && !R.snaps.empty()) {
            results = R.grid;
            grid = results.get();
            src = nullptr;
        }
        ViewKey k;
        k.prop = prop_ + 1000 * (isCompProp(pd.name) ? gCompSel : 0) + (std::string(pd.name) == "ARRAY" ? int(std::hash<std::string>{}(gArraySel) & 0x7fff) * 65536 : 0); k.snap = snap_; k.version = caseVersion_ + (isHmProp(pd.name) ? (hmViewVer_ << 20) : 0u); k.runVer = (pd.kind == Dynamic) ? runVer : 0;
        k.grid = grid ? (const void*)grid : (pd.kind == Initial && init_.valid ? (const void*)init_.grid.get() : (const void*)preview_.get());
        k.autoR = autoRange_; k.a = rmin_; k.b = rmax_; k.log = view_.settings().logScale; k.allSteps = settings_.rangeAllSteps;
        if (k.prop == key_.prop && k.snap == key_.snap && k.version == key_.version && k.runVer == key_.runVer && k.grid == key_.grid &&
            k.autoR == key_.autoR && k.a == key_.a && k.b == key_.b && k.log == key_.log && k.allSteps == key_.allSteps)
            return;
        key_ = k;
        if (pd.kind == Dynamic) {
            if (!grid) { vals.clear(); }
            else {
                snap_ = std::max(0, std::min(snap_, int(R.snaps.size()) - 1));
                const Snapshot& s = R.snaps[snap_];
                const std::string n = propKey(pd.name);
                if (isStreamlineProp(pd.name)) {
                    vals.assign(s.pressure.size(), NAN);
                    if (sl_.res && sl_.snap == snap_) {
                        const StreamlineResult& S = *sl_.res;
                        for (size_t a = 0; a < vals.size() && a < S.cellTof.size(); ++a) {
                            if (n == "TOF") vals[a] = S.cellTof[a] >= 0 ? S.cellTof[a] : NAN;
                            else if (n == "DRAINAGE") vals[a] = S.cellSink[a] >= 0 ? float(S.cellSink[a]) : NAN;
                            else vals[a] = S.cellSource[a] >= 0 ? float(S.cellSource[a]) : NAN;
                        }
                    }
                } else {
                    vals = s.field(n);
                    // C18.1 all-steps colour range (min / max over every report step), cached per run + property
                    grOk_ = false;
                    if (autoRange_ && settings_.rangeAllSteps) {
                        const std::string gk = std::to_string(runVer) + "|" + n + "|" + std::to_string(R.snaps.size());
                        if (gk != grKey_) {
                            float glo = FLT_MAX, ghi = -FLT_MAX;
                            for (const Snapshot& t : R.snaps)
                                for (float v : t.field(n)) if (std::isfinite(v)) { glo = std::min(glo, v); ghi = std::max(ghi, v); }
                            grKey_ = gk; grLo_ = glo; grHi_ = ghi;
                        }
                        grOk_ = grLo_ <= grHi_;
                    }
                }
                if (vals.size() != size_t(grid->numActive())) vals.assign(size_t(grid->numActive()), NAN);
                char b[128];
                std::snprintf(b, sizeof(b), "t = %.1f days (%.2f years)  |  report %d of %d", s.time, s.time / 365.25, snap_,
                              int(R.snaps.size()) - 1);
                caption = b;
                if (resultsCaseVersion_ != caseVersion_) caption += "  |  inputs changed since this run";
                if (grid != viewGrid_) {
                    view_.setGrid(grid);
                    viewHold_ = results;
                    viewGrid_ = grid;
                }
                view_.setWells(wellGlyphs(*grid, R.input, std::max(0, snap_ - 1)));
            }
        }
    }
    if (pd.kind == Initial) {
        if (init_.valid) {
            grid = init_.grid.get();
            const std::string n = pd.name;
            vals = n == "PRESSURE_0" ? init_.p : n == "SWAT_0" ? init_.sw : init_.sg;
            caption = "Initial state (equilibrium)";
            if (init_.caseVersion != caseVersion_) caption += "  |  inputs changed - recompute";
            if (grid != viewGrid_) {
                view_.setGrid(grid);
                viewHold_ = init_.grid;
                viewGrid_ = grid;
            }
            view_.setWells(wellGlyphs(*grid, case_));
        }
    }
    if (!grid) {    // static or fallback: model preview
        grid = preview_.get();
        if (grid && pd.kind == Static && grid->numCells() == (src ? src->grid : case_.grid).size()) {
            const GridInput& gi = src ? src->grid : case_.grid;
            const int na = grid->numActive();
            vals.resize(na);
            const std::string n = pd.name;
            for (int a = 0; a < na; ++a) {
                const int q = grid->globalOf[a];
                double v = 0;
                if (q >= int(gi.poro.size()) || q >= int(gi.permx.size())) {   // LGR / MINC / hydraulic fracture cells: not in the input arrays
                    v = NAN;
                    if (grid->isHf(q)) {
                        const Grid::HfCell& H = grid->hfCells[size_t(q - grid->hfBase())];
                        const bool okF = H.frac >= 0 && size_t(H.frac) < gi.hfracs.size();
                        if (n == "PERMX" || n == "PERMZ") v = H.kf;
                        else if (n == "PORO") v = okF ? gi.hfracs[size_t(H.frac)].poro : NAN;
                    }
                    if (n == "DEPTH") v = grid->depth[a];
                    else if (n == "PORV") v = grid->poreVolRef[a];
                    vals[a] = float(v);
                    continue;
                }
                if (n == "PORO") v = gi.poro[q];
                else if (n == "PERMX") v = gi.permx[q];
                else if (n == "PERMZ") v = gi.permz.empty() ? gi.permx[q] : gi.permz[q];
                else if (n == "NTG") v = gi.ntg.empty() ? 1.0 : gi.ntg[q];
                else if (n == "DEPTH") v = grid->depth[a];
                else if (n == "PORV") v = grid->poreVolRef[a];
                else if (isHmProp(n)) {
                    const std::vector<float>& H = n == "HMMULTK" ? hmMultK_ : n == "HMMULTPHI" ? hmMultPhi_ : hmSensK_;
                    v = size_t(q) < H.size() ? double(H[size_t(q)]) : NAN;
                }
                vals[a] = float(v);
            }
            caption = "Model input (" + std::to_string(grid->nx) + " x " + std::to_string(grid->ny) + " x " + std::to_string(grid->nz) + ")";
        } else if (pd.kind != Static && vals.empty()) {
            caption = pd.kind == Dynamic ? "No results yet - press Run (F5)" : "Press 'Compute initial state' (Initialization)";
        }
        if (grid && grid != viewGrid_) {
            view_.setGrid(grid);
            viewHold_ = preview_;
            view_.setWells(wellGlyphs(*grid, case_));
            viewGrid_ = grid;
        } else if (grid) view_.setWells(wellGlyphs(*grid, case_));
    }
    // range
    float lo = FLT_MAX, hi = -FLT_MAX;
    for (float v : vals) if (std::isfinite(v)) { lo = std::min(lo, v); hi = std::max(hi, v); }
    if (pd.kind == Dynamic && grOk_ && autoRange_ && settings_.rangeAllSteps && !isStreamlineProp(pd.name)) { lo = grLo_; hi = grHi_; }
    if (lo > hi) { lo = 0; hi = 1; }
    if (isSaturation(pd.name)) { lo = 0; hi = 1; }
    if (hi - lo < 1e-9f * std::max(1.0f, std::fabs(hi))) { hi = lo + std::max(1e-6f, 0.01f * std::fabs(lo)); }
    if (autoRange_) { rmin_ = lo; rmax_ = hi; key_.a = rmin_; key_.b = rmax_; }
    view_.setValues(vals, rmin_, rmax_, std::string(pd.name) + " - " + pd.label, unitOf(pd.name, case_.units));
    view_.setCaption(caption);
}

// ----------------------------------------------------------------------------- streamlines
void App::updateStreamlines() {
    ViewSettings& vs = view_.settings();
    const bool needProp = isStreamlineProp(kProps[prop_].name);
    vs.showCells = !(sl_.show && sl_.hideCells);
    vs.linesXray = sl_.xray;
    auto clearLines = [&]() { if (sl_.builtLines != -1) { view_.setLines({}); sl_.builtLines = -1; } };
    if (!sl_.show && !needProp) { clearLines(); return; }
    if (runner_.running()) return;
    std::shared_ptr<Simulator> sim; std::shared_ptr<const Grid> grid; int nsnap = 0;
    const unsigned ver = runner_.version();
    {
        std::lock_guard<std::mutex> l(runner_.mutex());
        const RunResults& R = runner_.results();
        sim = R.sim; grid = R.grid; nsnap = int(R.snaps.size());
    }
    if (!sim || nsnap == 0) { clearLines(); sl_.res.reset(); return; }
    const int s = std::clamp(snap_, 0, nsnap - 1);
    bool rebuilt = false;
    static int tracedFor = -1;
    if (!sl_.res || sl_.snap != s || sl_.runVer != ver || tracedFor != sl_.lines) {
        const auto t0 = std::chrono::steady_clock::now();
        try {
            sim->restoreReport(s);
            StreamlineOptions o; o.lines = sl_.lines;
            sl_.res = std::make_shared<StreamlineResult>(traceStreamlines(*sim, o));
        } catch (const std::exception& e) { runner_.addLog(std::string("Streamline tracing failed: ") + e.what()); sl_.res = std::make_shared<StreamlineResult>(); }
        sl_.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        sl_.snap = s; sl_.runVer = ver; sl_.grid = grid; tracedFor = sl_.lines;
        key_.prop = -1;   // refresh cell colours (TOF / regions)
        rebuilt = true;
    }
    if (!sl_.show) { clearLines(); return; }
    const auto& W = view_.wells();
    const std::string focus = vs.focusWell >= 0 && vs.focusWell < int(W.size()) ? W[size_t(vs.focusWell)].name : std::string();
    const std::string fkey = focus + (sl_.focusOnly ? "|1" : "|0");
    if (!rebuilt && sl_.builtLines == int(sl_.res->lines.size()) && sl_.builtColor == sl_.colorBy + 10 * int(vs.cmap) && sl_.builtFocus == fkey) return;
    const StreamlineResult& S = *sl_.res;
    std::vector<LineVertex> v = streamlineVertices(S, sl_.colorBy, vs.cmap, focus, sl_.focusOnly);
    sl_.builtFocus = fkey;
    view_.setLines(std::move(v));
    sl_.builtLines = int(S.lines.size());
    sl_.builtColor = sl_.colorBy + 10 * int(vs.cmap);
}

// ----------------------------------------------------------------------------- 3D viewer panel
void App::panelViewer() {
    ImGui::SetNextWindowSize(ImVec2(900, 600), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("3D Viewer")) { ImGui::End(); return; }
    if (!viewOk_) { ImGui::TextColored(col4(Pal::Err), "3D renderer failed to initialise (Direct3D 11 / shader compiler)."); ImGui::End(); return; }

    // ---- toolbar row
    std::set<std::string> extraHave;   // extra result arrays of this run (ZMF1 -> ZMF)
    std::vector<std::string> arrayNames;   // all extra arrays of the last report (generic group)
    {
        std::lock_guard<std::mutex> l(runner_.mutex());
        const RunResults& R = runner_.results();
        if (!R.snaps.empty()) for (const auto& e : R.snaps.back().extra) { extraHave.insert(baseName(e.first)); arrayNames.push_back(e.first); }
    }
    // generic engine arrays: everything not already offered as a named property
    {
        std::vector<std::string> keep;
        for (const auto& a : arrayNames) {
            const std::string b = baseName(a);
            bool named = false;
            for (int i = 0; i < kNumProps; ++i) if (kProps[i].kind == Dynamic && std::string(kProps[i].name) != "ARRAY" && (a == kProps[i].name || (isCompProp(kProps[i].name) && b == kProps[i].name))) named = true;
            if (!named) keep.push_back(a);
        }
        arrayNames.swap(keep);
        if (!arrayNames.empty()) extraHave.insert("ARRAY");
    }
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 13);
    const std::string propLabel = std::string(kProps[prop_].name) == "ARRAY" ? gArraySel + "  (engine array)" : std::string(kProps[prop_].label);
    if (ImGui::BeginCombo("##prop", propLabel.c_str(), ImGuiComboFlags_HeightLarge)) {
        const char* groups[3] = {"Static (input)", "Initial state", "Dynamic (results)"};
        for (int gk = 0; gk < 3; ++gk) {
            ImGui::SeparatorText(groups[gk]);
            for (int i = 0; i < kNumProps; ++i) {
                if (kProps[i].kind != gk) continue;
                if (std::string(kProps[i].name) == "RV" && !case_.phases.vapoil) continue;
                if ((std::string(kProps[i].name) == "POLYMER" || std::string(kProps[i].name) == "PADS") && !case_.phases.polymer) continue;
                if (std::string(kProps[i].name) == "XCO2" && !case_.comp.enabled) continue;
                if (isExtraProp(kProps[i].name) && !extraHave.count(kProps[i].name)) continue;
                if (std::string(kProps[i].name) == "HMMULTK" && hmMultK_.empty()) continue;
                if (std::string(kProps[i].name) == "HMMULTPHI" && hmMultPhi_.empty()) continue;
                if (std::string(kProps[i].name) == "HMSENSK" && hmSensK_.empty()) continue;
                if (std::string(kProps[i].name) == "ARRAY") continue;   // listed below by name
                char b[96]; std::snprintf(b, sizeof(b), "%-11s %s", kProps[i].name, kProps[i].label);
                if (ImGui::Selectable(b, prop_ == i)) { prop_ = i; autoRange_ = true; }
            }
        }
        if (!arrayNames.empty()) {
            ImGui::SeparatorText(TR("Engine output arrays (all restart arrays)"));
            int ai = -1; for (int i = 0; i < kNumProps; ++i) if (std::string(kProps[i].name) == "ARRAY") ai = i;
            for (const auto& a : arrayNames) {
                char b[96]; const char* u = arrayUnit(a, case_.units);
                std::snprintf(b, sizeof(b), "%-11s %s", a.c_str(), *u ? u : "");
                if (ImGui::Selectable(b, prop_ == ai && gArraySel == a)) { prop_ = ai; gArraySel = a; autoRange_ = true; }
            }
        }
        ImGui::EndCombo();
    }
    if (isCompProp(kProps[prop_].name)) {   // component of the mole-fraction map (run component list; else the case EOS names)
        std::vector<std::string> names;
        {
            std::lock_guard<std::mutex> l(runner_.mutex());
            for (const auto& ci : runner_.results().components) names.push_back(ci.name);
        }
        if (names.empty()) {
            if (case_.comp.co2store) names = {"CO2", "H2O"};
            else names = case_.comp.fluid.names;
            if (names.empty()) for (const auto& q : case_.thermal.comps) names.push_back(q.name);
        }
        const int nc = int(names.size());
        if (nc > 0) {
            gCompSel = std::clamp(gCompSel, 0, nc - 1);
            auto nm = [&](int i) { return names[size_t(i)]; };
            ImGui::SameLine(); ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6);
            if (ImGui::BeginCombo("##compsel", nm(gCompSel).c_str())) {
                for (int i = 0; i < nc; ++i) if (ImGui::Selectable(nm(i).c_str(), i == gCompSel)) { gCompSel = i; autoRange_ = true; }
                ImGui::EndCombo();
            }
        }
    }
    int nsnap = 0;
    {
        std::lock_guard<std::mutex> l(runner_.mutex());
        nsnap = int(runner_.results().snaps.size());
    }
    if (kProps[prop_].kind == Dynamic && nsnap > 0) {
        ImGui::SameLine();
        if (ImGui::Button(play_ ? ICON_PAUSE "##play" : ICON_PLAY "##play")) { play_ = !play_; if (play_ && snap_ >= nsnap - 1) snap_ = 0; }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip(play_ ? "Pause animation" : "Play animation through report steps");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 14);
        ImGui::SliderInt("##snap", &snap_, 0, std::max(0, nsnap - 1), "report %d");
        if (play_) {
            playAcc_ += ImGui::GetIO().DeltaTime;
            if (playAcc_ > 0.12) { playAcc_ = 0; if (++snap_ >= nsnap) { snap_ = nsnap - 1; play_ = false; } }
        }
    }
    ViewerBinding vb;
    vb.view = &view_; vb.autoRange = &autoRange_; vb.rmin = &rmin_; vb.rmax = &rmax_; vb.sl = &sl_;
    if (kProps[prop_].kind == Dynamic) vb.rangeAllSteps = &settings_.rangeAllSteps;
    vb.id = "model"; vb.propName = std::string(kProps[prop_].name) == "ARRAY" ? gArraySel.c_str() : kProps[prop_].name; vb.units = unitOf(kProps[prop_].name, case_.units);
    {
        std::lock_guard<std::mutex> l(runner_.mutex());
        const RunResults& R = runner_.results();
        vb.slAvailable = R.sim && !R.snaps.empty() && !runner_.running();
        vb.slHint = runner_.running() ? "Simulation running..." : R.valid && R.engine != "ResSim" ? "External engine results: open them in the Results workspace (streamlines are traced there from the restart file)."
                                                                : "Run the case with the ResSim engine first (F5).";
    }
    vb.grid = static_cast<const Grid*>(viewGrid_);
    if (!vb.grid) vb.grid = preview_.get();
    vb.pngName = case_.title.empty() ? std::string("model_view") : case_.title + "_" + kProps[prop_].name;
    if (kProps[prop_].kind == Dynamic && !isStreamlineProp(kProps[prop_].name)) {
        vb.cellHistory = [this](int cell, std::vector<double>& t, std::vector<double>& v) {
            const std::string n = propKey(kProps[prop_].name);
            std::lock_guard<std::mutex> l(runner_.mutex());
            for (const Snapshot& s : runner_.results().snaps) {
                const std::vector<float>& a = s.field(n);
                if (size_t(cell) >= a.size()) continue;
                t.push_back(s.time); v.push_back(a[size_t(cell)]);
            }
            return !t.empty();
        };
    }
    viewerToolbar(vb, mvUi_);

    updateStreamlines();
    updateViewData();
    updateFracturePanels();
    vb.grid = static_cast<const Grid*>(viewGrid_);
    viewerCanvas(vb, mvUi_, [this](int cell) {
        const Grid* g = static_cast<const Grid*>(viewGrid_);
        ImGui::Text("Cell %s", g->cellLabel(g->globalOf[size_t(cell)]).c_str());
        ImGui::Separator();
        ImGui::Text("%s = %.5g %s", std::string(kProps[prop_].name) == "ARRAY" ? gArraySel.c_str() : kProps[prop_].name, view_.valueOf(cell), unitOf(kProps[prop_].name, case_.units));
        ImGui::TextDisabled("depth %.1f %s, PV %.4g", g->depth[size_t(cell)], Units::of(case_.units).len, g->poreVolRef[size_t(cell)]);
    });
    ImGui::End();
}

// ----------------------------------------------------------------------------- plots
void App::panelPlots() {
    ImGui::SetNextWindowSize(ImVec2(900, 600), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Plots")) { ImGui::End(); return; }
    plotsBody();
    ImGui::End();
}

void App::plotsBody() {
    static bool years = false;
    std::lock_guard<std::mutex> lock(runner_.mutex());
    const RunResults& R = runner_.results();
    if (!R.valid || R.summary.empty()) {   // C40 empty state with the next actions
        const bool running = runner_.running();
        const int a = ui::emptyState(running ? ICON_ACTIVITY : ICON_CHART_LINE, running ? "Simulation running" : "No results yet",
                                     running ? "Profiles appear here as soon as the first report step is written. Progress: Run Monitor."
                                             : "Run the current model, or open results written earlier by ResSim, OPM Flow or ECLIPSE. "
                                               "Every finished run is kept in the Case Explorer (Objects > Runs) for comparison.",
                                     running ? std::vector<const char*>{} : std::vector<const char*>{ICON_PLAY "  Run simulation (F5)", ICON_FOLDER_OPEN "  Open ResSim results...", ICON_CHART_LINE "  Open ECLIPSE / OPM results..."});
        // actions need the runner lock released: run them at the end of the frame
        if (a == 0) rvDeferred_ = [this] { runSimulation(); };
        else if (a == 1) rvDeferred_ = [this] { std::string p; if (fileDialog(false, "ResSim results (*.rsres)\0*.rsres\0All files (*.*)\0*.*\0", "rsres", p)) openResultsFile(p); };
        else if (a == 2) rvDeferred_ = [this] { std::string p; if (fileDialog(false, "ECLIPSE / OPM Flow results\0*.SMSPEC;*.FSMSPEC;*.EGRID;*.UNRST;*.UNSMRY;*.DATA;*.data\0All files (*.*)\0*.*\0", "SMSPEC", p)) rvOpen(p); };
        return;
    }
    const auto& S = R.summary;
    const double xs = years ? 1.0 / 365.25 : 1.0;
    const Units U = Units::of(R.input.units);
    char liq[32], gas[32], pr[32];
    std::snprintf(liq, sizeof(liq), "%s/d", U.liq); std::snprintf(gas, sizeof(gas), "%s/d", U.gas); std::snprintf(pr, sizeof(pr), "%s", U.p);
    const bool hasGas = R.input.phases.gas, hasWater = R.input.phases.water || R.input.comp.co2store;   // CO2STORE: brine = water vectors
    double tSnap = -1;
    if (snap_ >= 0 && snap_ < int(R.snaps.size())) tSnap = R.snaps[snap_].time * xs;

    ImGui::Checkbox(TRL("Time in years"), &years);
    ImGui::SameLine();
    ImGui::TextDisabled("%s", TR("  Double-click a plot to move the 3D time cursor there."));
    static int tabNow = 0;
    static const char* tabFig[] = {"run.field", "run.well", "run.custom", "run.conv", "run.co2", "run.therm", "run.comp", "run.geo", "run.faults", "run.fields", "run.udq", "run.msw", "run.network"};
    ImGui::SameLine();
    exportButton(tabFig[tabNow]);
    auto cursor = [&]() {
        if (tSnap >= 0) {
            ImPlotSpec sp; sp.LineColor = ImVec4(1, 1, 1, 0.35f); sp.LineWeight = 1.0f; sp.Flags = ImPlotItemFlags_NoLegend;
            ImPlot::PlotInfLines("##cursor", &tSnap, 1, sp);
            recVLine(tSnap);
        }
        if (ImPlot::IsPlotHovered() && ImGui::IsMouseDoubleClicked(0) && !R.snaps.empty()) {
            const double t = ImPlot::GetPlotMousePos().x / xs;
            int best = 0;
            for (int i = 0; i < int(R.snaps.size()); ++i) if (std::fabs(R.snaps[i].time - t) < std::fabs(R.snaps[best].time - t)) best = i;
            snap_ = best;
            if (kProps[prop_].kind != Dynamic) prop_ = 9;   // PRESSURE
        }
    };
    const char* xl = years ? "Time (years)" : "Time (days)";
    const ImPlotAxisFlags af = ImPlotAxisFlags_AutoFit;

    if (ImGui::BeginTabBar("plots")) {
        if (ui::tabItem("Field")) {
            tabNow = 0;
            const ImVec2 av = ImGui::GetContentRegionAvail();
            const float fs = ImGui::GetFontSize();
            const PlotGrid G = fitGrid(4, av, fs * 22, fs * 11, 1.6f);
            if (G.scroll) ImGui::BeginChild("##fieldscroll", av);
            fig::recorder().beginFigure("run.field", G.rows, G.cols);
            if (ImPlot::BeginSubplots("##field", G.rows, G.cols, ImVec2(G.scroll ? av.x - ImGui::GetStyle().ScrollbarSize : av.x, G.size.y), ImPlotSubplotFlags_LinkAllX)) {
                if (ImPlot::BeginPlot(TRL("Production / injection rates"))) {
                    recPanel("Production / injection rates", xl, liq, hasGas ? gas : "");
                    ImPlot::SetupAxes(xl, liq, af, af); setupSI(ImAxis_Y1);
                    if (hasGas) { ImPlot::SetupAxis(ImAxis_Y2, gas, af | ImPlotAxisFlags_AuxDefault); setupSI(ImAxis_Y2); }
                    plotVec("Oil", S, FOPR, -1, xs, Pal::Oil);
                    if (hasWater) { plotVec("Water", S, FWPR, -1, xs, Pal::Water); plotVec("Water inj", S, FWIR, -1, xs, IM_COL32(120, 190, 255, 255), 1.5f); }
                    if (hasGas) { useY2(true); plotVec("Gas", S, FGPR, -1, xs, Pal::Gas); if (R.input.wells.size()) plotVec("Gas inj", S, FGIR, -1, xs, IM_COL32(255, 150, 140, 255), 1.5f); useY2(false); }
                    cursor();
                    recEnd();
                    ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot(TRL("Cumulative volumes"))) {
                    recPanel("Cumulative volumes", xl, U.liq, hasGas ? U.gas : "");
                    ImPlot::SetupAxes(xl, U.liq, af, af); setupSI(ImAxis_Y1);
                    if (hasGas) { ImPlot::SetupAxis(ImAxis_Y2, U.gas, af | ImPlotAxisFlags_AuxDefault); setupSI(ImAxis_Y2); }
                    plotVec("FOPT", S, FOPT, -1, xs, Pal::Oil);
                    if (hasWater) { plotVec("FWPT", S, FWPT, -1, xs, Pal::Water); plotVec("FWIT", S, FWIT, -1, xs, IM_COL32(120, 190, 255, 255), 1.5f); }
                    if (hasGas) { useY2(true); plotVec("FGPT", S, FGPT, -1, xs, Pal::Gas); useY2(false); }
                    cursor();
                    recEnd();
                    ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot(TRL("Pressure"))) {
                    recPanel("Pressure", xl, pr, "");
                    ImPlot::SetupAxes(xl, pr, af, af); setupSI(ImAxis_Y1);
                    plotVec("Average (FPR)", S, FPR, -1, xs, Pal::Pressure, 2.5f);
                    for (int w = 0; w < int(R.wellNames.size()); ++w) {
                        const std::string l = "BHP " + R.wellNames[w];
                        plotVec(l.c_str(), S, WBHP, w, xs, wellColor(w), 1.5f);
                    }
                    cursor();
                    recEnd();
                    ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot(TRL("Water cut / GOR"))) {
                    recPanel("Water cut / GOR", xl, "Water cut (-)", hasGas ? U.rs : "");
                    ImPlot::SetupAxes(xl, "Water cut (-)", af, 0);
                    ImPlot::SetupAxisLimits(ImAxis_Y1, 0, 1, ImPlotCond_Once);
                    if (hasGas) { ImPlot::SetupAxis(ImAxis_Y2, U.rs, af | ImPlotAxisFlags_AuxDefault); setupSI(ImAxis_Y2); }
                    if (hasWater) plotVec("FWCT", S, FWCT, -1, xs, Pal::Water);
                    if (hasGas) { useY2(true); plotVec("FGOR", S, FGOR, -1, xs, Pal::Gas); useY2(false); }
                    cursor();
                    recEnd();
                    ImPlot::EndPlot();
                }
                ImPlot::EndSubplots();
            }
            fig::recorder().endFigure();
            if (G.scroll) ImGui::EndChild();
            ImGui::EndTabItem();
        }
        if (ui::tabItem("Wells")) {
            tabNow = 1;
            const int nw = int(R.wellNames.size());
            if (nw == 0) ui::hint("The run has no wells.");
            else {
                resultsWell_ = std::min(resultsWell_, nw - 1);
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10);
                if (ImGui::BeginCombo("Well", R.wellNames[resultsWell_].c_str())) {
                    for (int w = 0; w < nw; ++w) if (ImGui::Selectable(R.wellNames[w].c_str(), w == resultsWell_)) resultsWell_ = w;
                    ImGui::EndCombo();
                }
                const int w = resultsWell_;
                const bool prod = R.wellTypes[w] == WellType::Producer;
                ImGui::SameLine();
                ImGui::TextDisabled("%s", toString(R.wellTypes[w]));
                if (!S.empty() && w < int(S.back().wells.size())) {
                    ImGui::SameLine();
                    ImGui::TextDisabled("  |  control now: %s", toString(S.back().wells[w].mode));
                }
                const PlotGrid G = fitGrid(3, ImGui::GetContentRegionAvail(), ImGui::GetFontSize() * 22, ImGui::GetFontSize() * 11, 1.5f);
                if (G.scroll) ImGui::BeginChild("##wellscroll", ImGui::GetContentRegionAvail());
                fig::recorder().beginFigure("run.well", G.rows, G.cols);
                if (ImPlot::BeginSubplots("##well", G.rows, G.cols, ImVec2(G.scroll ? -ImGui::GetStyle().ScrollbarSize : -1, G.size.y), ImPlotSubplotFlags_LinkAllX)) {
                    if (ImPlot::BeginPlot(prod ? "Rates (lines = simulated, dots = observed)" : "Injection rate")) {
                        recPanel(prod ? "Rates" : "Injection rate", xl, liq, prod && hasGas ? gas : "");
                        ImPlot::SetupAxes(xl, liq, af, af); setupSI(ImAxis_Y1);
                        if (prod) {
                            if (hasGas) { ImPlot::SetupAxis(ImAxis_Y2, gas, af | ImPlotAxisFlags_AuxDefault); setupSI(ImAxis_Y2); }
                            plotVec("WOPR", S, WOPR, w, xs, Pal::Oil);
                            plotVec("WOPRH", S, WOPRH, w, xs, Pal::Oil, 1, true);
                            if (hasWater) { plotVec("WWPR", S, WWPR, w, xs, Pal::Water); plotVec("WWPRH", S, WWPRH, w, xs, Pal::Water, 1, true); }
                            if (hasGas) { useY2(true); plotVec("WGPR", S, WGPR, w, xs, Pal::Gas); plotVec("WGPRH", S, WGPRH, w, xs, Pal::Gas, 1, true); useY2(false); }
                        } else {
                            if (R.wellTypes[w] == WellType::WaterInjector) plotVec("WWIR", S, WWIR, w, xs, Pal::Water);
                            else plotVec("WGIR", S, WGIR, w, xs, Pal::Gas);
                        }
                        cursor();
                        recEnd();
                        ImPlot::EndPlot();
                    }
                    if (ImPlot::BeginPlot(TRL("Bottom-hole pressure"))) {
                        recPanel("Bottom-hole pressure", xl, pr, "");
                        ImPlot::SetupAxes(xl, pr, af, af); setupSI(ImAxis_Y1);
                        plotVec("WBHP", S, WBHP, w, xs, wellColor(w));
                        plotVec("FPR", S, FPR, -1, xs, Pal::Pressure, 1.0f);
                        cursor();
                        recEnd();
                        ImPlot::EndPlot();
                    }
                    if (ImPlot::BeginPlot(TRL("Water cut / GOR"))) {
                        recPanel("Water cut / GOR", xl, "WWCT", hasGas ? U.rs : "");
                        ImPlot::SetupAxes(xl, "WWCT", af, 0);
                        ImPlot::SetupAxisLimits(ImAxis_Y1, 0, 1, ImPlotCond_Once);
                        if (hasGas) { ImPlot::SetupAxis(ImAxis_Y2, U.rs, af | ImPlotAxisFlags_AuxDefault); setupSI(ImAxis_Y2); }
                        if (prod && hasWater) plotVec("WWCT", S, WWCT, w, xs, Pal::Water);
                        if (prod && hasGas) { useY2(true); plotVec("WGOR", S, WGOR, w, xs, Pal::Gas); useY2(false); }
                        cursor();
                        recEnd();
                        ImPlot::EndPlot();
                    }
                    ImPlot::EndSubplots();
                }
                fig::recorder().endFigure();
                if (G.scroll) ImGui::EndChild();
            }
            ImGui::EndTabItem();
        }
        const bool co2solRun = !R.input.comp.co2store && R.input.comp.fluid.co2sol;
        const bool trackRun = R.input.thermal.enabled && !R.trackName.empty();   // thermal: inventory of the tracked component
        const std::string tn = trackRun ? R.trackName : std::string("CO2");
        const std::string co2Title = trackRun ? tn + " storage###co2tab" : co2solRun ? std::string("CO2 (EOR / CO2SOL)###co2tab") : std::string("CO2 storage###co2tab");
        if ((R.input.comp.co2store || co2solRun || trackRun) && ui::tabItem(co2Title.c_str())) {
            tabNow = 4;
            if (trackRun)
                ImGui::TextDisabled("Thermal run: %s inventory with the same classes as the compositional / CO2STORE cores (Henry dissolution in water, "
                                    "oil-phase solution, Land residual trapping, structural, mineral = solids with FORMULA trackEq).", tn.c_str());
            if (co2solRun)
                ImGui::TextDisabled("%s", TR("CO2SOL run: CO2 is an EOS component. 'Mobile' = CO2 held in the hydrocarbon phases (oil + gas), "
                                    "'dissolved' = CO2 dissolved in the brine (AMF / XCO2AQ maps in the 3D viewer)."));
            // trapping inventory (OPM FGMMO / FGMTR / FGMDS + ResSim structural / mineral): every CO2 molecule is counted once:
            // mineral + dissolved + residually trapped gas + structurally held gas + mobile gas = total
            const SummaryPoint* now = nullptr;
            for (const auto& p : S) if (tSnap < 0 || p.time * xs <= tSnap + 1e-9) now = &p;
            bool anyMin = false, anySalt = false, anyOil = false;
            for (const auto& p : S) { anyMin = anyMin || p.co2Min != 0; anySalt = anySalt || p.saltSolid > 0; anyOil = anyOil || p.co2Oil > 0; }
            if (now) {
                const double mob = now->co2Mob - now->co2Struct, tot = now->co2Mob + now->co2Trap + now->co2Diss + now->co2Min + now->co2Oil;
                auto pc = [&](double v) { return tot > 0 ? 100 * v / tot : 0.0; };
                ImGui::Text("Day %.0f: %s in place %.4g t  |  mobile %.4g t (%.1f %%)  |  structural %.4g t (%.1f %%)  |  residual %.4g t (%.1f %%)  |  "
                            "dissolved %.4g t (%.1f %%)  |  mineral %.4g t (%.1f %%)",
                            now->time, tn.c_str(), tot / 1e3, mob / 1e3, pc(mob), now->co2Struct / 1e3, pc(now->co2Struct), now->co2Trap / 1e3, pc(now->co2Trap),
                            now->co2Diss / 1e3, pc(now->co2Diss), now->co2Min / 1e3, pc(now->co2Min));
                if (anyOil) ImGui::TextDisabled("Dissolved in oil: %.4g t (%.1f %%)", now->co2Oil / 1e3, pc(now->co2Oil));
                if (anySalt || now->co2KrImmob > 0)
                    ImGui::TextDisabled("Gas in cells with krg = 0: %.4g t  |  precipitated NaCl: %.4g t", now->co2KrImmob / 1e3, now->saltSolid / 1e3);
            }
            const ImVec2 av = ImGui::GetContentRegionAvail();
            const float fs = ImGui::GetFontSize();
            const int npl = (anyMin || anySalt) ? 4 : 3;
            const PlotGrid G = fitGrid(npl, av, fs * 22, fs * 11, 1.6f);
            if (G.scroll) ImGui::BeginChild("##co2scroll", av);
            fig::recorder().beginFigure("run.co2", G.rows, G.cols);
            if (ImPlot::BeginSubplots("##co2", G.rows, G.cols, ImVec2(G.scroll ? av.x - ImGui::GetStyle().ScrollbarSize : av.x, G.size.y), ImPlotSubplotFlags_LinkAllX)) {
                struct Mech { const char* label; ImU32 col; };
                static const Mech kM[6] = {{"Mineral (FGMMN)", IM_COL32(160, 120, 220, 255)}, {"Dissolved (FGMDS)", IM_COL32(80, 170, 255, 255)},
                                           {"Residually trapped (FGMTR)", IM_COL32(250, 190, 60, 255)}, {"Structural (FGMST)", IM_COL32(120, 200, 120, 255)},
                                           {"Mobile", IM_COL32(240, 90, 80, 255)}, {"In oil", IM_COL32(150, 110, 60, 255)}};
                auto mechVal = [](const SummaryPoint& p, int m) {
                    switch (m) { case 0: return p.co2Min; case 1: return p.co2Diss; case 2: return p.co2Trap; case 3: return p.co2Struct; case 5: return p.co2Oil; default: return p.co2Mob - p.co2Struct; }
                };
                const std::string trapTitle = tn + " trapping (mass)", massAx = tn + " mass (kg)";
                if (ImPlot::BeginPlot(trapTitle.c_str())) {
                    recPanel(trapTitle.c_str(), xl, "kg", "");
                    ImPlot::SetupAxes(xl, massAx.c_str(), af, af); setupSI(ImAxis_Y1);
                    const int n = int(S.size());
                    std::vector<double> t(static_cast<size_t>(n)), lo(size_t(n), 0.0), hi(size_t(n), 0.0);
                    for (int i = 0; i < n; ++i) t[size_t(i)] = S[size_t(i)].time * xs;
                    for (int m = 0; m < 6 && n > 0; ++m) {
                        bool any = false;
                        for (int i = 0; i < n; ++i) { lo[size_t(i)] = hi[size_t(i)]; hi[size_t(i)] += std::max(0.0, mechVal(S[size_t(i)], m)); any = any || mechVal(S[size_t(i)], m) > 0; }
                        if (!any && m != 1 && m != 4) continue;
                        ImPlotSpec sp; sp.FillAlpha = 0.45f; sp.FillColor = col4(kM[m].col); sp.LineColor = sp.FillColor;
                        ImPlot::PlotShaded(kM[m].label, t.data(), lo.data(), hi.data(), n, sp);
                        if (fig::recorder().recording()) recSeries(kM[m].label, t.data(), hi.data(), n, kM[m].col, 2.0f, false);
                    }
                    plotVec("Total (FGMIP)", S, FGMIP, -1, xs, IM_COL32(230, 230, 230, 255), 1.5f);
                    cursor();
                    recEnd();
                    ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot(TRL("Trapping fractions"))) {
                    recPanel("Trapping fractions", xl, (std::string("fraction of ") + tn + " in place").c_str(), "");
                    ImPlot::SetupAxes(xl, "Fraction (-)", af, 0);
                    ImPlot::SetupAxisLimits(ImAxis_Y1, 0, 1, ImPlotCond_Once);
                    const int n = int(S.size());
                    for (int m = 0; m < 6; ++m) {
                        std::vector<double> t, f; bool any = false;
                        for (int i = 0; i < n; ++i) {
                            const SummaryPoint& p = S[size_t(i)];
                            const double tot = p.co2Mob + p.co2Trap + p.co2Diss + p.co2Min + p.co2Oil;
                            if (tot <= 0) continue;
                            t.push_back(p.time * xs); f.push_back(mechVal(p, m) / tot); any = any || mechVal(p, m) > 0;
                        }
                        if (t.empty() || (!any && m != 1 && m != 4)) continue;
                        ImPlotSpec sp; sp.LineColor = col4(kM[m].col); sp.LineWeight = 2.0f * g_plotLineScale;
                        ImPlot::PlotLine(kM[m].label, t.data(), f.data(), int(t.size()), sp);
                        if (fig::recorder().recording()) recSeries(kM[m].label, t.data(), f.data(), int(t.size()), kM[m].col, 2.0f, false);
                    }
                    cursor();
                    recEnd();
                    ImPlot::EndPlot();
                }
                if (npl == 4 && ImPlot::BeginPlot(TRL("Rock alteration"))) {
                    recPanel("Rock alteration", xl, "kg", "");
                    ImPlot::SetupAxes(xl, "Mass (kg)", af, af); setupSI(ImAxis_Y1);
                    if (anySalt) plotVec("Precipitated NaCl (FSALT)", S, FSALT, -1, xs, IM_COL32(220, 220, 220, 255), 2.0f);
                    if (anyMin) plotVec("CO2 in carbonate minerals (FGMMN)", S, FGMMN, -1, xs, IM_COL32(160, 120, 220, 255), 2.0f);
                    cursor();
                    recEnd();
                    ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot(TRL("Injection and pressure"))) {
                    recPanel("Injection and pressure", xl, gas, pr);
                    ImPlot::SetupAxes(xl, gas, af, af); setupSI(ImAxis_Y1);
                    ImPlot::SetupAxis(ImAxis_Y2, pr, af | ImPlotAxisFlags_AuxDefault); setupSI(ImAxis_Y2);
                    plotVec(trackRun ? "Gas injection (FGIR)" : "CO2 injection (FGIR)", S, FGIR, -1, xs, Pal::Gas);
                    plotVec(trackRun ? "Water production (FWPR)" : "Brine production (FWPR)", S, FWPR, -1, xs, Pal::Water, 1.5f);
                    if (trackRun) plotVec("Gas production (FGPR)", S, FGPR, -1, xs, IM_COL32(200, 160, 255, 255), 1.5f);
                    useY2(true);
                    plotVec("FPR", S, FPR, -1, xs, Pal::Pressure, 2.0f);
                    for (int w = 0; w < int(R.wellNames.size()); ++w) { const std::string l = "BHP " + R.wellNames[w]; plotVec(l.c_str(), S, WBHP, w, xs, wellColor(w), 1.2f); }
                    useY2(false);
                    cursor();
                    recEnd();
                    ImPlot::EndPlot();
                }
                ImPlot::EndSubplots();
            }
            fig::recorder().endFigure();
            if (G.scroll) ImGui::EndChild();
            ImGui::EndTabItem();
        }
        if (R.input.thermal.enabled && ui::tabItem("Thermal###thermtab")) {
            tabNow = 5;
            const bool fldU = R.input.units == UnitSystem::Field;
            const char* tU = fldU ? "Temperature (degF)" : "Temperature (degC)";
            const char* eU = fldU ? "Energy (Btu)" : "Energy (kJ)";
            const SummaryPoint* now = nullptr;
            for (const auto& p : S) if (tSnap < 0 || p.time * xs <= tSnap + 1e-9) now = &p;
            if (now)
                ImGui::TextDisabled("t = %.1f d | avg T %.1f, max T %.1f %s | steam PV %.2f %% | cSOR %.2f | heat lost %.3g of %.3g %s injected (%.1f %%)",
                                    now->time, now->favgT, now->fmaxT, fldU ? "degF" : "degC", 100 * now->fsteamPV, now->csor, now->fhlt,
                                    now->fenergyIn, fldU ? "Btu" : "kJ", now->fenergyIn > 0 ? 100 * now->fhlt / now->fenergyIn : 0.0);
            const bool hasRx = !R.input.thermal.reactions.empty(), hasGeo = R.input.thermal.geo.enabled;
            const PlotGrid G = fitGrid(4 + int(hasRx) + int(hasGeo), ImGui::GetContentRegionAvail(), ImGui::GetFontSize() * 22, ImGui::GetFontSize() * 11, 1.6f);
            if (G.scroll) ImGui::BeginChild("##thermscroll", ImGui::GetContentRegionAvail());
            fig::recorder().beginFigure("run.therm", G.rows, G.cols);
            if (ImPlot::BeginSubplots("##therm", G.rows, G.cols, ImVec2(G.scroll ? -ImGui::GetStyle().ScrollbarSize : -1, G.size.y), ImPlotSubplotFlags_LinkAllX)) {
                if (ImPlot::BeginPlot(TRL("Reservoir temperature"))) {
                    recPanel("Reservoir temperature", xl, tU, "");
                    ImPlot::SetupAxes(xl, tU, af, af);
                    plotVec("Average (FAVGT)", S, FAVGT, -1, xs, IM_COL32(255, 160, 60, 255), 2.5f);
                    plotVec("Maximum (FMAXT)", S, FMAXT, -1, xs, IM_COL32(230, 70, 60, 255), 1.5f);
                    cursor(); recEnd(); ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot(TRL("Steam zone"))) {
                    recPanel("Steam zone", xl, "Steam PV fraction (-)", "");
                    ImPlot::SetupAxes(xl, "Steam PV fraction (-)", af, af);
                    plotVec("FSTPV", S, FSTPV, -1, xs, IM_COL32(200, 200, 255, 255), 2.0f);
                    cursor(); recEnd(); ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot(TRL("Energy balance"))) {
                    recPanel("Energy balance", xl, eU, "");
                    ImPlot::SetupAxes(xl, eU, af, af); setupSI(ImAxis_Y1);
                    plotVec("Net enthalpy injected (FHIN)", S, FHIN, -1, xs, IM_COL32(255, 120, 80, 255), 2.0f);
                    plotVec("Heat lost to over/underburden (FHLT)", S, FHLT, -1, xs, IM_COL32(150, 150, 150, 255), 2.0f);
                    cursor(); recEnd(); ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot(TRL("Steam-oil ratio / oil rate"))) {
                    recPanel("Steam-oil ratio / oil rate", xl, "cSOR (CWE / oil)", liq);
                    ImPlot::SetupAxes(xl, "cSOR (-)", af, af);
                    ImPlot::SetupAxis(ImAxis_Y2, liq, af | ImPlotAxisFlags_AuxDefault); setupSI(ImAxis_Y2);
                    plotVec("cSOR", S, FCSOR, -1, xs, IM_COL32(120, 190, 255, 255), 2.0f);
                    useY2(true);
                    plotVec("Oil (FOPR)", S, FOPR, -1, xs, Pal::Oil, 1.5f);
                    plotVec("Steam CWE (FWIR)", S, FWIR, -1, xs, IM_COL32(255, 200, 120, 255), 1.5f);
                    useY2(false);
                    cursor(); recEnd(); ImPlot::EndPlot();
                }
                if (hasRx && ImPlot::BeginPlot(TRL("Reactions"))) {
                    recPanel("Reactions", xl, "Reaction heat rate (MW)", "Cumulative (GJ)");
                    ImPlot::SetupAxes(xl, "Heat rate (MW)", af, af);
                    ImPlot::SetupAxis(ImAxis_Y2, "Cumulative (GJ)", af | ImPlotAxisFlags_AuxDefault);
                    plotVec("Reaction heat rate (FRHR)", S, FRHR, -1, xs, IM_COL32(255, 90, 40, 255), 2.0f);
                    useY2(true);
                    plotVec("Cumulative reaction heat (FRHT)", S, FRHT, -1, xs, IM_COL32(255, 210, 90, 255), 1.5f);
                    useY2(false);
                    cursor(); recEnd(); ImPlot::EndPlot();
                }
                if (hasGeo && ImPlot::BeginPlot(TRL("Geomechanics"))) {
                    recPanel("Geomechanics", xl, "Heave (mm)", "Max failure index (-)");
                    ImPlot::SetupAxes(xl, "Heave (mm)", af, af);
                    ImPlot::SetupAxis(ImAxis_Y2, "Failure index (-)", af | ImPlotAxisFlags_AuxDefault);
                    plotVec("Surface heave (FUPL)", S, FUPL, -1, xs, IM_COL32(120, 220, 140, 255), 2.0f);
                    plotVec("Reservoir-top heave (FUPR)", S, FUPR, -1, xs, IM_COL32(60, 160, 255, 255), 1.5f);
                    useY2(true);
                    plotVec("Max failure index (FMXF)", S, FMXF, -1, xs, IM_COL32(230, 70, 60, 255), 1.5f);
                    useY2(false);
                    cursor(); recEnd(); ImPlot::EndPlot();
                }
                ImPlot::EndSubplots();
            }
            fig::recorder().endFigure();
            if (G.scroll) ImGui::EndChild();
            ImGui::EndTabItem();
        }
        if (R.input.thermal.enabled && R.wellT.size() >= 2 && ui::tabItem("Geothermal###geotab")) {
            tabNow = 7;
            geothermalTab(R, xs, tSnap, xl, cursor);
            ImGui::EndTabItem();
        }
        if ((R.faultMech || R.geoCouple || R.hfDdm || R.hfPf) && ui::tabItem("Faults & seismicity###faulttab")) {
            tabNow = 8;
            faultsTab(R, xs, tSnap, xl, cursor);
            ImGui::EndTabItem();
        }
        const bool trRun = !R.components.empty() && R.components[0].kind == 4;   // C27 black-oil tracer run
        if (!R.components.empty() && ui::tabItem(trRun ? "Tracers###comptab" : "Components###comptab")) {
            tabNow = 6;
            if (trRun) tracersTab(R, xs, tSnap, xl, cursor); else componentsTab(R, xs, tSnap, xl, cursor);
            ImGui::EndTabItem();
        }
        if (!R.wellNames.empty() && !R.input.comp.enabled && !R.input.thermal.enabled && ui::tabItem("Production analysis###prodtab")) {
            prodAnalysisTab(R);
            ImGui::EndTabItem();
        }
        if (R.grid && !R.snaps.empty() && ui::tabItem("Fields & sections###fieldtab")) {   // C30
            tabNow = 9;
            fieldsTab(R);
            ImGui::EndTabItem();
        }
        {   // C31 UDQ vectors and fired ACTIONX blocks
            const bool anyUdq = !R.actions.empty() || (!R.summary.empty() && !R.summary.back().udq.empty()) ||
                                std::any_of(R.input.extraKeywords.begin(), R.input.extraKeywords.end(), [](const ressim::ExtraKeyword& e) { return e.name == "UDQ" || e.name == "ACTIONX"; });
            if (anyUdq && ui::tabItem("UDQ & actions###udqtab")) {
                tabNow = 10;
                udqTab(R, xs, xl);
                ImGui::EndTabItem();
            }
        }
        if (!R.msw.empty() && R.mswRec && ui::tabItem("Multisegment wells###mswtab")) {   // C32
            tabNow = 11;
            mswTab(R, xs, xl);
            ImGui::EndTabItem();
        }
        if (!R.netNodes.empty() && R.netRec && ui::tabItem("Surface network###nettab")) {   // C35-4b
            tabNow = 12;
            networkTab(R, xs, xl);
            ImGui::EndTabItem();
        }
        if (ui::tabItem("Custom plots###customtab")) {
            tabNow = 2;
            customPlotsTab(R, xs, tSnap, xl, cursor);
            ImGui::EndTabItem();
        }
        if (ui::tabItem("Convergence")) {
            tabNow = 3;
            fig::recorder().beginFigure("run.conv", 1, 2);
            if (ImPlot::BeginSubplots("##conv", 1, 2, ImVec2(-1, std::max(ImGui::GetContentRegionAvail().y, ImGui::GetFontSize() * 13.0f)), ImPlotSubplotFlags_LinkAllX)) {
                if (ImPlot::BeginPlot(TRL("Time-step size"))) {
                    recPanel("Time-step size", xl, "dt (days)", "");
                    ImPlot::SetupAxes(xl, "dt (days)", af, af); setupSI(ImAxis_Y1);
                    plotVec("dt", S, DT, -1, xs, Pal::Accent, 1.5f, true);
                    recEnd();
                    ImPlot::EndPlot();
                }
                if (ImPlot::BeginPlot(TRL("Iterations per time step"))) {
                    recPanel("Iterations per time step", xl, "Newton", "Linear");
                    ImPlot::SetupAxes(xl, "Newton", af, af);
                    ImPlot::SetupAxis(ImAxis_Y2, "Linear", af | ImPlotAxisFlags_AuxDefault);
                    plotVec("Newton", S, NEWTON, -1, xs, Pal::Warn, 1.5f, true);
                    useY2(true);
                    plotVec("Linear", S, LINEAR, -1, xs, Pal::Muted, 1.0f);
                    useY2(false);
                    recEnd();
                    ImPlot::EndPlot();
                }
                ImPlot::EndSubplots();
            }
            fig::recorder().endFigure();
            ImGui::EndTabItem();
        }
        // C14: user plots marked "own tab" (star) - each its own page after any run
        for (int i = 0; i < int(cplots_.size()); ++i) {
            if (!cplots_[size_t(i)].ownTab) continue;
            char lab[300];
            std::snprintf(lab, sizeof(lab), ICON_STAR " %s###cptab%d", cplots_[size_t(i)].name.c_str(), i);
            if (ui::tabItem(lab)) {
                tabNow = 2;
                fig::recorder().beginFigure("run.custom", 1, 1);
                drawCustomPlot(cplots_[size_t(i)], R, xs, xl, cursor, ImVec2(-1, -1), 1000 + i);
                fig::recorder().endFigure();
                ImGui::EndTabItem();
            }
        }
        ImGui::EndTabBar();
    }
}

// ----------------------------------------------------------------------------- run monitor
// ============================================================================ multi-component product analysis
// Per-component production / injection of any multi-component run (compositional, CO2STORE, thermal / combustion /
// UCG): stream selector (field or a well), component rate / composition / cumulative plots, product identification
// (H2, syngas CO + H2, H2/CO, CO2 / O2 slip), heating values (LHV / HHV power and MJ/sm3) and a per-well breakdown.
namespace {
void plotXY(const char* label, const std::vector<double>& X, const std::vector<double>& Y, ImU32 col, float w) {
    const int n = int(std::min(X.size(), Y.size()));
    if (n <= 0) return;
    if (fig::recorder().recording()) recSeries(label, X.data(), Y.data(), n, col, w);
    ImPlotSpec sp; sp.LineColor = col4(col); sp.LineWeight = w * g_plotLineScale;
    ImPlot::PlotLine(label, X.data(), Y.data(), n, sp);
}
std::string fmtMol(double mol) {   // cumulative moles -> readable
    char b[48];
    const double a = std::fabs(mol);
    if (a >= 1e9) std::snprintf(b, sizeof(b), "%.4g Gmol", mol * 1e-9);
    else if (a >= 1e6) std::snprintf(b, sizeof(b), "%.4g Mmol", mol * 1e-6);
    else if (a >= 1e3) std::snprintf(b, sizeof(b), "%.4g kmol", mol * 1e-3);
    else std::snprintf(b, sizeof(b), "%.4g mol", mol);
    return b;
}
} // namespace

void App::componentsTab(const RunResults& R, double xs, double tSnap, const char* xl, const std::function<void()>& cursor) {
    static int stream = 0;   // 0 field, 1.. well
    static int basis = 0;    // 0 kmol/d, 1 t/d, 2 sm3/d (gas-like) / kmol/d
    const auto& S = R.summary;
    const auto& C = R.components;
    const int nc = int(C.size());
    const float fs = ImGui::GetFontSize();
    std::vector<std::string> streams{"Field"};
    for (const auto& n : R.wellNames) streams.push_back("Well " + n);
    if (stream >= int(streams.size())) stream = 0;
    ImGui::SetNextItemWidth(fs * 12);
    if (ImGui::BeginCombo("Stream", streams[size_t(stream)].c_str())) {
        for (int i = 0; i < int(streams.size()); ++i) if (ImGui::Selectable(streams[size_t(i)].c_str(), stream == i)) stream = i;
        ImGui::EndCombo();
    }
    ImGui::SameLine();
    static const char* kBasis[] = {"kmol/d", "tonne/d", "sm3/d"};
    ImGui::SetNextItemWidth(fs * 8);
    ImGui::Combo("Rate basis", &basis, kBasis, 3);
    ImGui::SameLine();
    ImGui::TextDisabled("  %d components: %s", nc, [&]() { static std::string s; s.clear(); for (const auto& c : C) { if (!s.empty()) s += ", "; s += c.name; } return s.c_str(); }());
    // rates of the stream at summary point i (production > 0 vector and injection > 0 vector)
    auto ratesAt = [&](size_t i, std::vector<double>& prod, std::vector<double>& inj) {
        prod.assign(size_t(nc), 0.0); inj.assign(size_t(nc), 0.0);
        const SummaryPoint& p = S[i];
        if (stream == 0) {
            for (int k = 0; k < nc; ++k) {
                if (size_t(k) < p.fcompP.size()) prod[size_t(k)] = p.fcompP[size_t(k)];
                if (size_t(k) < p.fcompI.size()) inj[size_t(k)] = p.fcompI[size_t(k)];
            }
        } else if (size_t(stream - 1) < p.wells.size()) {
            const auto& q = p.wells[size_t(stream - 1)].comp;
            for (int k = 0; k < nc && size_t(k) < q.size(); ++k) { if (q[size_t(k)] > 0) prod[size_t(k)] = q[size_t(k)]; else inj[size_t(k)] = -q[size_t(k)]; }
        }
    };
    auto cumAt = [&](size_t i, std::vector<double>& cp, std::vector<double>& ci) {
        cp.assign(size_t(nc), 0.0); ci.assign(size_t(nc), 0.0);
        const SummaryPoint& p = S[i];
        const std::vector<double>* a = stream == 0 ? &p.fcompPT : size_t(stream - 1) < p.wells.size() ? &p.wells[size_t(stream - 1)].compCumP : nullptr;
        const std::vector<double>* b = stream == 0 ? &p.fcompIT : size_t(stream - 1) < p.wells.size() ? &p.wells[size_t(stream - 1)].compCumI : nullptr;
        for (int k = 0; k < nc; ++k) {
            if (a && size_t(k) < a->size()) cp[size_t(k)] = (*a)[size_t(k)];
            if (b && size_t(k) < b->size()) ci[size_t(k)] = (*b)[size_t(k)];
        }
    };
    auto conv = [&](int k, double molps) {   // mol/s -> plot basis
        if (basis == 1) return molps * C[size_t(k)].mw * 86.4;                 // t/d
        if (basis == 2) return molps * 86400.0 / kMolPerSm3;                   // sm3/d (ideal gas)
        return molps * 86.4;                                                   // kmol/d
    };
    // time series of the selected stream (rates are report-step averages: skip t = 0)
    const size_t i0 = S.size() > 1 ? 1 : 0;
    std::vector<double> X;
    std::vector<std::vector<double>> Q(static_cast<size_t>(nc)), Y(static_cast<size_t>(nc)), CP(static_cast<size_t>(nc));
    std::vector<double> h2, co, syn, ratio, lhvP, hhvP, lhvG, co2s, o2s;
    size_t nowIdx = 0;
    std::vector<double> prod, inj, cp, ci;
    for (size_t i = i0; i < S.size(); ++i) {
        ratesAt(i, prod, inj); cumAt(i, cp, ci);
        const StreamAnalysis A = analyseStream(C, prod, &inj);
        X.push_back(S[i].time * xs);
        for (int k = 0; k < nc; ++k) { Q[size_t(k)].push_back(conv(k, prod[size_t(k)])); Y[size_t(k)].push_back(A.yDry[size_t(k)]); CP[size_t(k)].push_back(cp[size_t(k)] * 1e-3); }
        h2.push_back(A.h2 * 86.4); co.push_back(A.co * 86.4); syn.push_back(A.syngas * 86.4); ratio.push_back(A.h2co);
        lhvP.push_back(A.lhvRate * 1e-6); hhvP.push_back(A.hhvRate * 1e-6); lhvG.push_back(A.lhvGas);
        co2s.push_back(A.co2Slip); o2s.push_back(A.o2Slip);
        if (tSnap < 0 || S[i].time * xs <= tSnap + 1e-9) nowIdx = i;
    }
    if (X.empty()) { ui::hint("No report steps yet."); return; }
    ratesAt(nowIdx, prod, inj); cumAt(nowIdx, cp, ci);
    const StreamAnalysis N = analyseStream(C, prod, &inj);
    // ---- product identification at the time cursor
    {
        std::vector<double> cP, cI; cumAt(S.size() - 1, cP, cI);
        double cumH2 = 0, cumCO = 0, cumLHV = 0;
        for (int k = 0; k < nc; ++k) {
            if (C[size_t(k)].species == SP_H2) cumH2 += cp[size_t(k)];
            if (C[size_t(k)].species == SP_CO) cumCO += cp[size_t(k)];
            cumLHV += cp[size_t(k)] * C[size_t(k)].lhv;
        }
        ImGui::TextColored(col4(Pal::Accent), ICON_FLASK_CONICAL " t = %.1f d", S[nowIdx].time);
        ImGui::SameLine();
        if (N.h2 > 0 || N.co > 0)
            ImGui::Text("| H2 %.4g kmol/d (%.4g t/d, %.4g sm3/d) | syngas CO+H2 %.4g kmol/d | H2/CO %.3f | H2 in dry gas %.2f %%",
                        N.h2 * 86.4, N.h2Mass * 86.4, N.h2 * 86400 / kMolPerSm3, N.syngas * 86.4, N.h2co, N.gas > 0 ? 100 * N.h2 / N.gas : 0.0);
        else ImGui::Text("| dry production %.4g kmol/d", N.dry * 86.4);
        ImGui::TextDisabled("Heating value of the produced stream: LHV %.4g MW, HHV %.4g MW | produced gas LHV %.3f / HHV %.3f MJ/sm3 | CO2 slip %.1f %% | O2 slip %.1f %% | "
                            "cumulative H2 %s, CO %s, LHV energy %.4g GJ",
                            N.lhvRate * 1e-6, N.hhvRate * 1e-6, N.lhvGas, N.hhvGas, 100 * N.co2Slip, 100 * N.o2Slip,
                            fmtMol(cumH2).c_str(), fmtMol(cumCO).c_str(), cumLHV * 1e-9);
    }
    // ---- component table at the cursor
    if (ImGui::CollapsingHeader(TRL("Component table"), ImGuiTreeNodeFlags_DefaultOpen) &&
        ui::beginTable("##ctab", 10, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_ScrollY,
                          ImVec2(0, std::min(fs * 1.6f * float(nc + 2), fs * 14)))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        const char* hd[] = {"Component", "Species", "MW (g/mol)", "Production (kmol/d)", "Mass (t/d)", "Injection (kmol/d)",
                            "Dry mole frac.", "Gas mole frac.", "Cum. produced", "Cum. injected"};
        for (const char* h : hd) ui::tableColumn(h);
        ImGui::TableHeadersRow();
        for (int k = 0; k < nc; ++k) {
            const ComponentInfo& c = C[size_t(k)];
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::PushStyleColor(ImGuiCol_Text, col4(seriesColor(k))); ImGui::TextUnformatted(c.name.c_str()); ImGui::PopStyleColor();
            ImGui::TableNextColumn();
            if (c.kind == 2) ImGui::TextDisabled("%s", TR("solid"));
            else if (c.species >= 0) ImGui::TextUnformatted(speciesFormula(c.species));
            else ImGui::TextDisabled(c.kind == 1 ? "gas" : "pseudo");
            ImGui::TableNextColumn(); ImGui::Text("%.2f", c.mw * 1e3);
            ImGui::TableNextColumn(); ImGui::Text("%.5g", prod[size_t(k)] * 86.4);
            ImGui::TableNextColumn(); ImGui::Text("%.5g", prod[size_t(k)] * c.mw * 86.4);
            ImGui::TableNextColumn(); ImGui::Text("%.5g", inj[size_t(k)] * 86.4);
            ImGui::TableNextColumn(); ImGui::Text("%.4f", N.yDry[size_t(k)]);
            ImGui::TableNextColumn(); if (N.yGas[size_t(k)] > 0) ImGui::Text("%.4f", N.yGas[size_t(k)]); else ImGui::TextDisabled("%s", TR("-"));
            ImGui::TableNextColumn(); ImGui::TextUnformatted(fmtMol(cp[size_t(k)]).c_str());
            ImGui::TableNextColumn(); ImGui::TextUnformatted(fmtMol(ci[size_t(k)]).c_str());
        }
        ui::endTable();
    }
    // ---- per-well breakdown (cumulative production per component + products)
    if (R.wellNames.size() > 1 && ImGui::CollapsingHeader(TRL("Per-well breakdown (cumulative to the cursor)"))) {
        if (ui::beginTable("##cwell", nc + 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_ScrollX)) {
            ui::tableColumn("Well");
            for (const auto& c : C) ui::tableColumn((c.name + " (kmol)").c_str());
            ui::tableColumn("H2 (t)"); ui::tableColumn("Syngas (kmol)"); ui::tableColumn("LHV (GJ)");
            ImGui::TableHeadersRow();
            const SummaryPoint& p = S[nowIdx];
            for (size_t w = 0; w < R.wellNames.size() && w < p.wells.size(); ++w) {
                const auto& cpw = p.wells[w].compCumP;
                double tot = 0; for (double v : cpw) tot += v;
                if (tot <= 0) continue;
                ImGui::TableNextRow(); ImGui::TableNextColumn(); ImGui::TextUnformatted(R.wellNames[w].c_str());
                double h2m = 0, sy = 0, e = 0;
                for (int k = 0; k < nc; ++k) {
                    const double v = size_t(k) < cpw.size() ? cpw[size_t(k)] : 0.0;
                    ImGui::TableNextColumn(); ImGui::Text("%.5g", v * 1e-3);
                    if (C[size_t(k)].species == SP_H2) { h2m += v * C[size_t(k)].mw * 1e-3; sy += v; }
                    if (C[size_t(k)].species == SP_CO) sy += v;
                    e += v * C[size_t(k)].lhv;
                }
                ImGui::TableNextColumn(); ImGui::Text("%.5g", h2m);
                ImGui::TableNextColumn(); ImGui::Text("%.5g", sy * 1e-3);
                ImGui::TableNextColumn(); ImGui::Text("%.5g", e * 1e-9);
            }
            ui::endTable();
        }
    }
    // ---- plots
    const bool products = std::any_of(C.begin(), C.end(), [](const ComponentInfo& c) { return c.species == SP_H2 || c.species == SP_CO; });
    const bool heat = std::any_of(C.begin(), C.end(), [](const ComponentInfo& c) { return c.lhv > 0; });
    const int np = 3 + int(products) + int(heat);
    const PlotGrid G = fitGrid(np, ImGui::GetContentRegionAvail(), fs * 22, fs * 11, 1.6f);
    if (G.scroll) ImGui::BeginChild("##compscroll", ImGui::GetContentRegionAvail());
    fig::recorder().beginFigure("run.comp", G.rows, G.cols);
    const ImPlotAxisFlags af = ImPlotAxisFlags_AutoFit;
    char ql[48]; std::snprintf(ql, sizeof(ql), "Rate (%s)", kBasis[basis]);
    if (ImPlot::BeginSubplots("##comp", G.rows, G.cols, ImVec2(G.scroll ? -ImGui::GetStyle().ScrollbarSize : -1, G.size.y), ImPlotSubplotFlags_LinkAllX)) {
        if (ImPlot::BeginPlot(TRL("Component production rates"))) {
            recPanel("Component production rates", xl, ql, "");
            ImPlot::SetupAxes(xl, ql, af, af); setupSI(ImAxis_Y1);
            for (int k = 0; k < nc; ++k) if (C[size_t(k)].kind != 2) plotXY(C[size_t(k)].name.c_str(), X, Q[size_t(k)], seriesColor(k), 2.0f);
            cursor(); recEnd(); ImPlot::EndPlot();
        }
        if (ImPlot::BeginPlot(TRL("Produced composition (dry basis)"))) {
            recPanel("Produced composition (dry basis)", xl, "Mole fraction (-)", "");
            ImPlot::SetupAxes(xl, "Mole fraction (-)", af, af);
            for (int k = 0; k < nc; ++k) if (C[size_t(k)].kind != 2 && C[size_t(k)].kind != 3 && C[size_t(k)].species != SP_H2O)
                plotXY(C[size_t(k)].name.c_str(), X, Y[size_t(k)], seriesColor(k), 2.0f);
            cursor(); recEnd(); ImPlot::EndPlot();
        }
        if (ImPlot::BeginPlot(TRL("Cumulative component production"))) {
            recPanel("Cumulative component production", xl, "Cumulative (kmol)", "");
            ImPlot::SetupAxes(xl, "Cumulative (kmol)", af, af); setupSI(ImAxis_Y1);
            for (int k = 0; k < nc; ++k) if (C[size_t(k)].kind != 2) plotXY(C[size_t(k)].name.c_str(), X, CP[size_t(k)], seriesColor(k), 2.0f);
            cursor(); recEnd(); ImPlot::EndPlot();
        }
        if (products && ImPlot::BeginPlot(TRL("Products: hydrogen / syngas"))) {
            recPanel("Products: hydrogen / syngas", xl, "Rate (kmol/d)", "H2/CO (-)");
            ImPlot::SetupAxes(xl, "Rate (kmol/d)", af, af); setupSI(ImAxis_Y1);
            ImPlot::SetupAxis(ImAxis_Y2, "H2/CO (-)", af | ImPlotAxisFlags_AuxDefault);
            plotXY("H2", X, h2, IM_COL32(90, 200, 255, 255), 2.5f);
            plotXY("CO", X, co, IM_COL32(255, 150, 60, 255), 2.0f);
            plotXY("Syngas (CO + H2)", X, syn, IM_COL32(180, 120, 255, 255), 2.0f);
            useY2(true);
            plotXY("H2/CO ratio", X, ratio, IM_COL32(230, 230, 230, 255), 1.5f);
            useY2(false);
            cursor(); recEnd(); ImPlot::EndPlot();
        }
        if (heat && ImPlot::BeginPlot(TRL("Heating value of the produced stream"))) {
            recPanel("Heating value of the produced stream", xl, "Power (MW)", "Gas LHV (MJ/sm3)");
            ImPlot::SetupAxes(xl, "Power (MW)", af, af);
            ImPlot::SetupAxis(ImAxis_Y2, "Gas LHV (MJ/sm3)", af | ImPlotAxisFlags_AuxDefault);
            plotXY("LHV power", X, lhvP, IM_COL32(255, 200, 80, 255), 2.5f);
            plotXY("HHV power", X, hhvP, IM_COL32(255, 120, 60, 255), 1.5f);
            useY2(true);
            plotXY("Produced gas LHV", X, lhvG, IM_COL32(120, 220, 140, 255), 1.5f);
            useY2(false);
            cursor(); recEnd(); ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
    fig::recorder().endFigure();
    if (G.scroll) ImGui::EndChild();
}

// C27 passive tracers of the black-oil engine (TRACER / WTRACER / RESSIM TRACERPROP): breakthrough curves, produced
// concentration, cumulative recovery of the injected tracer - the classic interwell tracer test evaluation
void App::tracersTab(const RunResults& R, double xs, double tSnap, const char* xl, const std::function<void()>& cursor) {
    static int stream = 0;   // 0 field, 1.. well
    const auto& S = R.summary;
    const auto& C = R.components;
    const int nt = int(C.size());
    const float fs = ImGui::GetFontSize();
    std::vector<std::string> streams{"Field"};
    for (const auto& n : R.wellNames) streams.push_back("Well " + n);
    if (stream >= int(streams.size())) stream = 0;
    ImGui::SetNextItemWidth(fs * 12);
    if (ImGui::BeginCombo("Stream", streams[size_t(stream)].c_str())) {
        for (int i = 0; i < int(streams.size()); ++i) if (ImGui::Selectable(streams[size_t(i)].c_str(), stream == i)) stream = i;
        ImGui::EndCombo();
    }
    ImGui::SameLine();
    ImGui::TextDisabled("  %d tracer(s); amounts in the tracer unit of the deck, concentration = amount per surface volume of the carrier phase", nt);
    const Units U = Units::of(R.input.units);
    auto phaseRate = [&](const SummaryPoint& p, int ph) {
        if (stream == 0) return ph == 0 ? p.fwpr : ph == 1 ? p.fopr : p.fgpr;
        if (size_t(stream - 1) >= p.wells.size()) return 0.0;
        const WellRates& w = p.wells[size_t(stream - 1)];
        return ph == 0 ? w.water : ph == 1 ? w.oil : w.gas;
    };
    auto rateOf = [&](const SummaryPoint& p, int k, bool inj) {
        if (stream == 0) { const auto& v = inj ? p.fcompI : p.fcompP; return size_t(k) < v.size() ? v[size_t(k)] : 0.0; }
        if (size_t(stream - 1) >= p.wells.size()) return 0.0;
        const auto& q = p.wells[size_t(stream - 1)].comp;
        const double v = size_t(k) < q.size() ? q[size_t(k)] : 0.0;
        return inj ? std::max(0.0, -v) : std::max(0.0, v);
    };
    auto cumOf = [&](const SummaryPoint& p, int k, bool inj) {
        const std::vector<double>* v = stream == 0 ? (inj ? &p.fcompIT : &p.fcompPT)
                                     : size_t(stream - 1) < p.wells.size() ? (inj ? &p.wells[size_t(stream - 1)].compCumI : &p.wells[size_t(stream - 1)].compCumP) : nullptr;
        return v && size_t(k) < v->size() ? (*v)[size_t(k)] : 0.0;
    };
    // field injected total (recovery is measured against the field injection, also for a single producer)
    auto fieldInj = [&](const SummaryPoint& p, int k) { return size_t(k) < p.fcompIT.size() ? p.fcompIT[size_t(k)] : 0.0; };
    const size_t i0 = S.size() > 1 ? 1 : 0;
    std::vector<double> X;
    const size_t ntu = static_cast<size_t>(nt);
    std::vector<std::vector<double>> Q(ntu), CC(ntu), CP(ntu), REC(ntu);
    std::vector<double> tBreak(ntu, -1.0), peakC(ntu, 0.0), tPeak(ntu, 0.0);
    size_t nowIdx = 0;
    for (size_t i = i0; i < S.size(); ++i) {
        X.push_back(S[i].time * xs);
        for (int k = 0; k < nt; ++k) {
            const double q = rateOf(S[i], k, false) * 86400.0;                 // amount / d
            const double ph = phaseRate(S[i], C[size_t(k)].species);
            const double c = ph > 0 ? q / ph : 0.0;
            Q[size_t(k)].push_back(q); CC[size_t(k)].push_back(c);
            CP[size_t(k)].push_back(cumOf(S[i], k, false));
            const double fi = fieldInj(S[i], k);
            REC[size_t(k)].push_back(fi > 0 ? 100.0 * cumOf(S[i], k, false) / fi : 0.0);
            if (c > peakC[size_t(k)]) { peakC[size_t(k)] = c; tPeak[size_t(k)] = S[i].time; }
        }
        if (tSnap < 0 || S[i].time * xs <= tSnap + 1e-9) nowIdx = i;
    }
    if (X.empty()) { ui::hint("No report steps yet."); return; }
    for (int k = 0; k < nt; ++k)   // breakthrough: produced concentration first above 1 % of its peak
        for (size_t j = 0; j < CC[size_t(k)].size(); ++j)
            if (peakC[size_t(k)] > 0 && CC[size_t(k)][j] >= 0.01 * peakC[size_t(k)]) { tBreak[size_t(k)] = S[i0 + j].time; break; }
    static const char* kPh[] = {"water", "oil", "gas"};
    if (ImGui::CollapsingHeader(TRL("Tracer table (to the time cursor)"), ImGuiTreeNodeFlags_DefaultOpen) &&
        ui::beginTable("##trtab", 9, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_ScrollY,
                       ImVec2(0, std::min(fs * 1.6f * float(nt + 2), fs * 12)))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        const char* hd[] = {"Tracer", "Phase", "Production (amount/d)", "Produced conc.", "Cum. produced", "Cum. injected (stream)",
                            "Recovery of field injection", "Breakthrough (1 % of peak)", "Peak conc. at"};
        for (const char* h : hd) ui::tableColumn(h);
        ImGui::TableHeadersRow();
        const SummaryPoint& p = S[nowIdx];
        for (int k = 0; k < nt; ++k) {
            const ComponentInfo& c = C[size_t(k)];
            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::PushStyleColor(ImGuiCol_Text, col4(seriesColor(k))); ImGui::TextUnformatted(c.name.c_str()); ImGui::PopStyleColor();
            ImGui::TableNextColumn(); ImGui::TextUnformatted(kPh[std::clamp(c.species, 0, 2)]);
            const double q = rateOf(p, k, false) * 86400.0, ph = phaseRate(p, c.species);
            ImGui::TableNextColumn(); ImGui::Text("%.5g", q);
            ImGui::TableNextColumn(); ImGui::Text("%.5g /%s", ph > 0 ? q / ph : 0.0, c.species == 2 ? U.gas : U.liq);
            ImGui::TableNextColumn(); ImGui::Text("%.5g", cumOf(p, k, false));
            ImGui::TableNextColumn(); ImGui::Text("%.5g", cumOf(p, k, true));
            const double fi = fieldInj(p, k);
            ImGui::TableNextColumn(); if (fi > 0) ImGui::Text("%.2f %%", 100.0 * cumOf(p, k, false) / fi); else ImGui::TextDisabled("%s", TR("-"));
            ImGui::TableNextColumn(); if (tBreak[size_t(k)] >= 0) ImGui::Text("%.4g d", tBreak[size_t(k)]); else ImGui::TextDisabled("%s", TR("not yet"));
            ImGui::TableNextColumn(); if (peakC[size_t(k)] > 0) ImGui::Text("%.4g d (%.4g)", tPeak[size_t(k)], peakC[size_t(k)]); else ImGui::TextDisabled("%s", TR("-"));
        }
        ui::endTable();
    }
    ImGui::TextDisabled("%s", TR("3D: Results viewer -> property 'Engine output array' -> TR_<name> (cell concentration per report step)."));
    const PlotGrid G = fitGrid(3, ImGui::GetContentRegionAvail(), fs * 22, fs * 11, 1.6f);
    if (G.scroll) ImGui::BeginChild("##trscroll", ImGui::GetContentRegionAvail());
    fig::recorder().beginFigure("run.tracer", G.rows, G.cols);
    const ImPlotAxisFlags af = ImPlotAxisFlags_AutoFit;
    if (ImPlot::BeginSubplots("##tracer", G.rows, G.cols, ImVec2(G.scroll ? -ImGui::GetStyle().ScrollbarSize : -1, G.size.y), ImPlotSubplotFlags_LinkAllX)) {
        if (ImPlot::BeginPlot(TRL("Produced tracer concentration (breakthrough curve)"))) {
            recPanel("Produced tracer concentration (breakthrough curve)", xl, "Concentration (amount / surface volume)", "");
            ImPlot::SetupAxes(xl, "Concentration (amount / surface volume)", af, af); setupSI(ImAxis_Y1);
            for (int k = 0; k < nt; ++k) plotXY(C[size_t(k)].name.c_str(), X, CC[size_t(k)], seriesColor(k), 2.0f);
            cursor(); recEnd(); ImPlot::EndPlot();
        }
        if (ImPlot::BeginPlot(TRL("Tracer production rate"))) {
            recPanel("Tracer production rate", xl, "Rate (amount/d)", "");
            ImPlot::SetupAxes(xl, "Rate (amount/d)", af, af); setupSI(ImAxis_Y1);
            for (int k = 0; k < nt; ++k) plotXY(C[size_t(k)].name.c_str(), X, Q[size_t(k)], seriesColor(k), 2.0f);
            cursor(); recEnd(); ImPlot::EndPlot();
        }
        if (ImPlot::BeginPlot(TRL("Cumulative recovery of the injected tracer"))) {
            recPanel("Cumulative recovery of the injected tracer", xl, "Recovered (% of field injection)", "");
            ImPlot::SetupAxes(xl, "Recovered (% of field injection)", af, af);
            for (int k = 0; k < nt; ++k) plotXY(C[size_t(k)].name.c_str(), X, REC[size_t(k)], seriesColor(k), 2.0f);
            cursor(); recEnd(); ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
    fig::recorder().endFigure();
    if (G.scroll) ImGui::EndChild();
}

static std::string bigNum(double v) {
    char b[64];
    const double a = std::fabs(v);
    if (a >= 1e9) std::snprintf(b, sizeof(b), "%.3f B", v / 1e9);
    else if (a >= 1e6) std::snprintf(b, sizeof(b), "%.3f MM", v / 1e6);
    else if (a >= 1e3) std::snprintf(b, sizeof(b), "%.2f M", v / 1e3);
    else std::snprintf(b, sizeof(b), "%.4g", v);
    return b;
}

void App::panelRun() {
    if (!ImGui::Begin("Run Monitor")) { ImGui::End(); return; }
    const bool running = runner_.running();
    const double prog = runner_.progress();
    char ov[96];
    std::snprintf(ov, sizeof(ov), running ? "Running  %.1f%%" : (prog >= 0.999 ? "Finished" : "Idle  %.0f%%"), prog * 100);
    ImGui::PushStyleColor(ImGuiCol_PlotHistogram, running ? col4(Pal::Accent) : col4(Pal::Ok));
    ImGui::ProgressBar(float(prog), ImVec2(-1, 0), ov);
    ImGui::PopStyleColor();
    const std::string err = runner_.error();
    if (!err.empty()) { ImGui::PushTextWrapPos(0); ImGui::TextColored(col4(Pal::Err), "%s", err.c_str()); ImGui::PopTextWrapPos(); }

    std::lock_guard<std::mutex> lock(runner_.mutex());
    const RunResults& R = runner_.results();
    if (!R.valid) { ui::hint("Run statistics, current rates and recovery appear here during a simulation."); ImGui::End(); return; }
    const Units U = Units::of(R.input.units);
    const SummaryPoint& s = R.summary.back();
    if (ui::beginTable("kpi", 4, ImGuiTableFlags_SizingStretchSame)) {
        auto kpi = [&](const char* label, const std::string& value, ImU32 col) {
            ImGui::TableNextColumn();
            ImGui::TextDisabled("%s", label);
            ImGui::PushFont(g_fonts.title ? g_fonts.title : g_fonts.ui, ImGui::GetStyle().FontSizeBase * 1.25f);
            ImGui::TextColored(col4(col), "%s", value.c_str());
            ImGui::PopFont();
        };
        char b[64];
        std::snprintf(b, sizeof(b), "%.0f / %.0f d", s.time, R.totalTime); kpi("Simulated time", b, IM_COL32(230, 232, 240, 255));
        std::snprintf(b, sizeof(b), "%.1f", s.fpr); kpi((std::string("Avg pressure (") + U.p + ")").c_str(), b, Pal::Pressure);
        std::snprintf(b, sizeof(b), "%.1f %%", R.initial.oil > 0 ? 100.0 * s.fopt / R.initial.oil : 0.0);
        kpi("Oil recovery factor", b, Pal::Oil);
        std::snprintf(b, sizeof(b), "%.2f s", runner_.wallSeconds()); kpi("Wall-clock time", b, IM_COL32(230, 232, 240, 255));
        ui::endTable();
    }
    ImGui::Separator();
    if (ui::beginTable("rates", 2, ImGuiTableFlags_SizingStretchProp)) {
        auto row = [&](const char* k, const std::string& v) { ImGui::TableNextRow(); ImGui::TableNextColumn(); ImGui::TextDisabled("%s", k); ImGui::TableNextColumn(); ImGui::TextUnformatted(v.c_str()); };
        char b[96];
        std::snprintf(b, sizeof(b), "%.1f %s/d", s.fopr, U.liq); row("Oil rate", b);
        std::snprintf(b, sizeof(b), "%.1f %s/d   (WCT %.3f)", s.fwpr, U.liq, s.fwct); row("Water rate", b);
        std::snprintf(b, sizeof(b), "%.1f %s/d   (GOR %.4g)", s.fgpr, U.gas, s.fgor); row("Gas rate", b);
        std::snprintf(b, sizeof(b), "%.1f %s/d", s.fwir, U.liq); row("Water injection", b);
        row("Cumulative oil", bigNum(s.fopt) + " " + U.liq);
        row("STOIIP (initial)", bigNum(R.initial.oil + R.initial.oilVaporized) + " " + U.liq);
        if (R.initial.oilVaporized > 0) row("  of which vaporised in gas", bigNum(R.initial.oilVaporized) + " " + U.liq);
        row("Gas in place (initial)", bigNum(R.initial.gasFree + R.initial.gasDissolved) + " " + U.gas);
        std::snprintf(b, sizeof(b), "%d steps, %d cuts, %d Newton, %d linear", R.stats.timesteps, R.stats.cuts, R.stats.newton, R.stats.linear);
        row("Solver", b);
        ui::endTable();
    }
    ImGui::End();
}
