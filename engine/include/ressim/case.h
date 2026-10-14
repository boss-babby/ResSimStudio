// Input data model of a simulation case (what the user edits in the UI).
// Mirrors the ECLIPSE sections: RUNSPEC / GRID / PROPS / SOLUTION / SCHEDULE.
#pragma once
#include <string>
#include <vector>
#include <array>
#include <algorithm>

namespace ressim {

// ---------------------------------------------------------------- RUNSPEC
enum class UnitSystem { Metric, Field };

struct Units {
    double cDarcy;        // transmissibility constant
    double gravity;       // pressure gradient per (density * depth)
    double gasDensity;    // converts Rs*rho_gsc and b_g*rho_gsc into reservoir density units
    double poreVolume;    // converts length^3 into reservoir volume unit (rb or rm3)
    double atm;           // 1 atmosphere in pressure units
    const char* p; const char* len; const char* liq; const char* gas; const char* rs; const char* dens; const char* rv = "";
    static Units of(UnitSystem u);
};

struct Phases {
    bool water = true, oil = true, gas = true, disgas = true;
    bool vapoil = false;              // VAPOIL: oil vaporised in the gas phase (wet gas / gas condensate, Rv, PVTG)
    bool polymer = false;             // POLYMER: polymer component transported in the water phase (EOR)
};

// ---------------------------------------------------------------- GRID
// Fault face list (FAULTS): cells I1..I2, J1..J2, K1..K2 (zero-based), face X/Y/Z, minus = X-/Y-/Z-
struct FaultFace { int i1 = 0, i2 = 0, j1 = 0, j2 = 0, k1 = 0, k2 = 0; char face = 'X'; bool minus = false; };
struct FaultInput { std::string name = "F1"; double mult = 1.0; std::vector<FaultFace> faces; };   // FAULTS + MULTFLT
struct NNCInput { int i1 = 0, j1 = 0, k1 = 0, i2 = 0, j2 = 0, k2 = 0; double trans = 0; };       // explicit NNC (zero-based)

struct GridInput {
    int nx = 10, ny = 10, nz = 3;
    // per-cell arrays (size nx*ny*nz, natural ECLIPSE ordering i fastest, then j, then k)
    std::vector<double> dx, dy, dz, tops;     // tops: depth of the top face of each cell
    std::vector<double> poro, permx, permy, permz, ntg;
    std::vector<int> actnum;
    std::vector<double> multx, multy, multz;  // transmissibility multipliers (optional, size 0 or n)
    // corner-point geometry (ECLIPSE COORD/ZCORN). When zcorn is non-empty it replaces DX/DY/DZ/TOPS.
    std::vector<double> coord;                // 6*(nx+1)*(ny+1): pillar top (x,y,z) and bottom (x,y,z)
    std::vector<double> zcorn;                // 8*nx*ny*nz corner depths, ECLIPSE ordering
    std::vector<FaultInput> faults;
    std::vector<NNCInput> nnc;
    // region numbers (1-based, size 0 = every cell in region 1)
    std::vector<int> satnum, pvtnum, eqlnum, fipnum;
    std::vector<double> multpv;               // MULTPV pore-volume multipliers (empty = 1)
    std::vector<double> swatinit;             // SWATINIT initial water saturation (empty = none, < 0 = not set for the cell):
                                              // equilibration scales the water-oil capillary pressure per cell to honour it
    // saturation end-point scaling (ENDSCALE + SWL / SWCR / SWU / SGL / SGCR / SGU / SOWCR / SOGCR, ECLIPSE two-point
    // scaling of kr, Pc scaled between SWL and SWU / SGL and SGU). Empty = table end points; NaN or < 0 = table value.
    bool endscale = false;
    std::vector<double> swl, swcr, swu, sgl, sgcr, sgu, sowcr, sogcr;
    std::vector<double> spoly;                // SPOLY initial polymer concentration (empty = 0)
    std::vector<int> plmixnum;                // PLMIXNUM polymer mixing regions (empty = 1)
    bool hasEndPointScaling() const {
        return endscale && !(swl.empty() && swcr.empty() && swu.empty() && sgl.empty() && sgcr.empty() && sgu.empty() && sowcr.empty() && sogcr.empty());
    }
    double minpv = 0;                         // MINPV: cells with a smaller pore volume are made inactive
    // transmissibility of block-centred (DX/DY/DZ/TOPS) grids: false = OLDTRAN (logical i,j,k neighbours, ECLIPSE
    // default), true = NEWTRAN (geometric face overlap incl. displaced columns -> NNCs, what OPM Flow always does)
    bool newtran = false;
    // dual porosity (ECLIPSE DUALPORO / DUALPERM): NZ is doubled, layers 1..NZ/2 are the matrix and NZ/2+1..NZ the
    // fracture system of the same blocks (the fracture cells take the matrix geometry). Each matrix cell couples to its
    // fracture cell through T = CDARCY * V * sigma * K (Warren-Root / Kazemi shape factor). DUALPORO: matrix cells have
    // no mutual connections and wells connect to fracture cells only; DUALPERM: the matrix grid also flows.
    // PINCH: cells thinner than 'threshold' are pinched out (inactive) and vertical non-neighbour connections join the
    // active cells above and below a run of pinched / inactive cells (ECLIPSE item 2 GAP: inactive cells of any
    // thickness are bridged, NOGAP: only thin ones; item 3: maximum total gap; item 4 TOPBOT / ALL; item 5 MULTZ TOP / ALL)
    struct Pinch { bool on = false; double threshold = 0.001; bool gap = true; double maxGap = 1e20; bool topbot = true; bool multzAll = false; } pinch;
    int dual = 0;                             // 0 single porosity, 1 DUALPORO, 2 DUALPERM
    double sigma = 0;                         // SIGMA: shape factor (1/length^2) for every block
    std::vector<double> sigmav;               // SIGMAV per cell (matrix half used; empty = SIGMA; NaN = SIGMA)
    bool ltosigma = false;                    // LTOSIGMA: sigma from matrix block sizes LX/LY/LZ
    double ltoF[3] = {4.0, 4.0, 4.0};         // LTOSIGMA factors fx, fy, fz (Kazemi: 4)
    std::vector<double> lx, ly, lz;           // matrix block sizes (length)
    std::vector<double> multmf;               // MULTMF matrix-fracture transmissibility multipliers
    bool nodppm = false;                      // NODPPM: fracture permeability is NOT multiplied by fracture porosity
    bool gravdr = false;                      // GRAVDR requested (gravity drainage term; not modelled, kept for export)
    std::vector<double> dzmtrxv;              // DZMTRXV matrix block heights (kept for export / gravity drainage)
    int nzMatrix() const { return dual ? nz / 2 : nz; }
    bool cornerPoint() const { return !zcorn.empty(); }
    int size() const { return nx * ny * nz; }
    int index(int i, int j, int k) const { return i + nx * (j + ny * k); }
};

// ---------------------------------------------------------------- PROPS
struct PVTW { double pref = 4014.7, bw = 1.029, cw = 3.13e-6, visc = 0.31, cv = 0.0; };
struct RockC { double pref = 14.7, cr = 3.0e-6; };
struct DensitySC { double oil = 49.1, water = 64.79, gas = 0.06054; };

struct PVDRow { double p, B, mu; };          // PVDO / PVDG row
struct PVTORecord {                          // one PVTO record (saturated row + optional undersaturated rows)
    double rs, pb, bo, mu;
    std::vector<PVDRow> under;               // p > pb rows (p, Bo, mu)
};

struct PVTGRow { double rv, B, mu; };       // PVTG row (Rv, Bg, mu_g)
struct PVTGRecord {                          // one PVTG record: gas pressure, saturated row first, then undersaturated rows
    double pg = 0;                           // (decreasing Rv, usually down to Rv = 0)
    std::vector<PVTGRow> rows;
};

struct SatRow { double s, krw, kro, pc; };   // SWOF row (Sw, krw, krow, Pcow) or SGOF row (Sg, krg, krog, Pcog)

// PVT / rock data of one PVTNUM region and saturation functions of one SATNUM region
struct PvtRegion {
    PVTW pvtw; RockC rock; DensitySC density;
    std::vector<PVDRow> pvdo, pvdg;
    std::vector<PVTORecord> pvto;
    std::vector<PVTGRecord> pvtg;
};
struct SatRegion { std::vector<SatRow> swof, sgof; };

struct PropsInput {
    // ---- region 1 (what the input panels edit)
    PVTW pvtw;
    RockC rock;
    DensitySC density;
    std::vector<PVDRow> pvdo;        // dead oil (used when disgas == false)
    std::vector<PVTORecord> pvto;    // live oil (used when disgas == true)
    std::vector<PVDRow> pvdg;        // dry gas (with PVTG: the saturated rows, kept in sync)
    std::vector<PVTGRecord> pvtg;    // wet gas (used when vapoil == true)
    std::vector<SatRow> swof;
    std::vector<SatRow> sgof;
    // three-phase oil relative permeability: 0 = ECLIPSE default (saturation-weighted), 1 = STONE1, 2 = STONE2
    int kroModel = 0;
    double stone1Exp = 1.0;          // STONE1EX exponent (STONE1 only)
    // ---- regions 2..N (PVTNUM / SATNUM), usually from imported decks
    std::vector<PvtRegion> pvtExtra;
    std::vector<SatRegion> satExtra;

    int numPvtRegions() const { return 1 + int(pvtExtra.size()); }
    int numSatRegions() const { return 1 + int(satExtra.size()); }
    PvtRegion pvtRegion(int r) const {           // r: 0-based
        if (r > 0 && r <= int(pvtExtra.size())) return pvtExtra[size_t(r - 1)];
        return PvtRegion{pvtw, rock, density, pvdo, pvdg, pvto, pvtg};
    }
    SatRegion satRegion(int r) const {
        if (r > 0 && r <= int(satExtra.size())) return satExtra[size_t(r - 1)];
        return SatRegion{swof, sgof};
    }
};

// ---------------------------------------------------------------- SOLUTION
struct EquilInput {
    double datumDepth = 8400, datumPressure = 4800;
    double owc = 8500, pcowOWC = 0;
    double goc = 8200, pcogGOC = 0;
    int nacc = 0;                     // 0 = cell centre; |n| -> 2|n| slices per cell
    // live oil initial solution GOR: <0 -> saturated at GOC pressure (ECLIPSE default), else constant Rs
    double initialRs = -1;
    std::vector<double> rsvdDepth, rsvdRs;   // RSVD table (Rs vs depth); overrides initialRs when present
    // wet gas initial Rv: <0 -> saturated at the GOC pressure (ECLIPSE default, capped by the local dew point), else constant
    double initialRv = -1;
    std::vector<double> rvvdDepth, rvvdRv;   // RVVD table (Rv vs depth); overrides initialRv when present
};

// Enumerated initial conditions (PRESSURE / SWAT / SGAS / RS arrays) - used instead of EQUIL when present.
struct InitialState {
    std::vector<double> pressure, swat, sgas, rs, rv;   // global cell order, size n (sgas / rs / rv optional)
    bool active() const { return !pressure.empty(); }
};

// ---------------------------------------------------------------- AQUIFERS (AQUFETP / AQUCT + AQUANCON)
enum class AquiferType { Fetkovich, CarterTracy };
struct AquiferConnection { int i1 = 0, i2 = 0, j1 = 0, j2 = 0, k1 = 0, k2 = 0; char face = 'X'; bool minus = false; double mult = 1.0; };
struct AquiferInput {
    int id = 1;
    AquiferType type = AquiferType::Fetkovich;
    double depth = 0;              // datum depth of the aquifer
    double p0 = -1;                // initial pressure (<0: equilibrium with the reservoir at the aquifer datum)
    double ct = 1e-5;              // total (rock + water) compressibility, 1/pressure
    double volume = 1e9;           // Fetkovich: initial water volume (reservoir volume units)
    double pi = 10;                // Fetkovich: productivity index (res. volume / day / pressure)
    double perm = 100, poro = 0.2; // Carter-Tracy: permeability (mD) and porosity
    double rInner = 3000;          // Carter-Tracy: inner radius (reservoir equivalent radius)
    double thickness = 50;         // Carter-Tracy: thickness
    double angle = 360;            // Carter-Tracy: angle of influence (degrees)
    double visc = -1;              // water viscosity (<0: from PVTW)
    std::vector<AquiferConnection> conns;
};

// ---------------------------------------------------------------- SCHEDULE
enum class WellType { Producer, WaterInjector, GasInjector };
enum class Control  { ORAT, WRAT, GRAT, LRAT, RESV, BHP, RATE };
enum class WellStatus { Open, Shut, Stop };

struct Completion {
    int i = 0, j = 0, k1 = 0, k2 = 0;   // zero-based indices; k1..k2 inclusive
    double cf = -1;                     // connection factor (<=0 -> computed, Peaceman)
    double kh = -1;                     // <=0 -> computed
    double skin = 0;
    double diameter = 0.5;
    char dir = 'Z';
    bool open = true;
    int fromStep = 0;                   // report step at which the connection is added (COMPDAT later in SCHEDULE)
    int untilStep = -1;                 // report step at which it is shut (COMPDAT ... SHUT later), -1 = never
    int satTable = 0;                   // COMPDAT item 6: saturation table for the connection (0 = cell's SATNUM); OPM only
    double dfactor = -1, r0 = -1;       // COMPDAT items 11, 13: D-factor and pressure-equivalent radius (<0 default); OPM only
};

struct WellInput {
    std::string name = "P1";
    std::string group = "G";
    WellType type = WellType::Producer;
    int i = 0, j = 0;                  // wellhead location (zero-based)
    double refDepth = -1;              // BHP reference depth (<=0: top completion)
    double piMult = 1.0;               // WPIMULT: multiplier on all connection factors
    std::vector<Completion> completions;
    std::string specsTail;             // WELSPECS items 7.. (drainage radius, inflow eq., auto shut-in, crossflow, ...) as deck text
};

// A control record (WCONPROD / WCONINJE / WCONHIST), applied from a report step onwards.
struct WellControl {
    std::string well;
    WellStatus status = WellStatus::Open;
    Control mode = Control::ORAT;
    double orat = -1, wrat = -1, grat = -1, lrat = -1, resv = -1;  // producer targets/limits (<0: none)
    double rate = -1;                  // injector surface rate target
    double bhp = -1;                   // BHP target (BHP mode) or limit (rate modes)
    bool history = false;              // WCONHIST: rates are observed data
    int phase = 0;                     // injectors: injected phase for this period (0 = the well's type, 1 = water, 2 = gas)
    // deck items the ResSim engine does not use but OPM Flow / ECLIPSE do (kept for export)
    std::string deckMode;              // original control mode when not modelled (GRUP, THP, NONE, CRAT, RESV for injectors)
    double thp = -1;                   // THP target / limit
    int vfp = 0;                       // VFP table number (0 = none)
    double alq = -1;                   // artificial lift quantity
};

// WPOLYMER: polymer concentration of the injected water (well name or pattern; stays in force until changed)
struct WellPolymer { std::string well; double conc = 0; };

struct ReportStep {
    double length = 30;                // days
    std::vector<WellControl> controls; // applied at the start of this step
    std::vector<WellPolymer> wpolymer; // WPOLYMER records applied at the start of this step
    // DRSDT set at the start of this step: < -1.5 unchanged, < 0 off, >= 0 maximum Rs increase per day
    double drsdt = -2;
    bool drsdtFree = false;            // DRSDT 'FREE': only cells with free gas are limited (default 'ALL')
};

// ---------------------------------------------------------------- polymer flooding (ECLIPSE POLYMER option)
// Polymer is a water-phase component: c [mass / surface water volume]. Effects modelled as in ECLIPSE / OPM Flow:
// viscosity increase with Todd-Longstaff mixing (PLYVISC, PLMIXPAR, PLYMAX), adsorption on the rock (PLYADS,
// PLYROCK, optional irreversibility), permeability reduction RRF and the inaccessible (dead) pore volume.
struct PolymerRock {
    double deadPoreVolume = 0;         // PLYROCK 1: inaccessible pore volume fraction
    double rrf = 1;                    // PLYROCK 2: residual resistance factor (>= 1)
    double rockDensity = 0;            // PLYROCK 3: in-situ rock density (mass / reservoir volume)
    int adsIndex = 1;                  // PLYROCK 4: 1 = desorption allowed, 2 = irreversible adsorption
    double maxAdsorption = 0;          // PLYROCK 5: maximum adsorbed polymer (mass / rock mass)
};
struct PolymerInput {
    std::vector<std::vector<double>> plyvisc;   // per PVTNUM region: rows (c, viscosity multiplier), c ascending, mult(0) = 1
    std::vector<std::vector<double>> plyads;    // per SATNUM region: rows (c, adsorbed mass / rock mass), ads(0) = 0
    std::vector<PolymerRock> plyrock;           // per SATNUM region
    // per PLMIXNUM region (index 0 = region 1)
    std::vector<double> cmax{0.0};              // PLYMAX 1: maximum polymer concentration (Todd-Longstaff reference)
    std::vector<double> cmaxSalt{0.0};          // PLYMAX 2 (kept for export)
    std::vector<double> mixing{1.0};            // PLMIXPAR: Todd-Longstaff mixing parameter omega (1 = full mixing)
    // shear thinning / thickening (PLYSHLOG), per PVTNUM region; empty = no shear effect
    struct Shear { double refConc = 0, refSalt = -1, refTemp = -1; std::vector<double> rows; };   // rows: (water velocity, multiplier)
    std::vector<Shear> plyshlog;
    int numMixRegions() const { return int(std::min({cmax.size(), mixing.size()})); }
};

// ---------------------------------------------------------------- numerical controls (TUNING-like)
struct Tuning {
    double initialDt = 1.0, maxDt = 60.0, minDt = 1e-5, growth = 2.0, cutFactor = 0.3;
    int maxNewton = 14, maxLinear = 300;
    double cnvTol = 1e-3, mbTol = 1e-7, linTol = 1e-5;
    double dsMax = 0.2, dpMaxRel = 0.3;              // OPM Flow default dp-max-rel 0.3
    int threads = 0;                   // 0 = hardware concurrency
    int linearSolver = 2;              // 0 = ILU(0)-BiCGSTAB, 1 = CPR-AMG, 2 = auto (CPR for larger models)
    int cprWeights = 1;                // CPR decoupling: 0 = quasi-IMPES (Jacobian diagonal block), 1 = true-IMPES (accumulation)
};

// ---------------------------------------------------------------- observed data (history matching)
enum class ObsQty { WOPR, WWPR, WGPR, WWIR, WBHP, WWCT, WGOR, FOPR, FWPR, FGPR, FPR, FWCT, FGOR };
struct Observation {
    ObsQty qty = ObsQty::FOPR;
    std::string well;                  // empty for field quantities
    double time = 0, value = 0;
    double sigma = -1;                 // measurement error (<=0: automatic)
    double weight = 1;                 // relative weight in the objective function (0 = ignore)
};

// ---------------------------------------------------------------- uncertainty / history-matching parameters
enum class ParamKind {
    PermXY,          // multiplier on PERMX/PERMY in a box
    PermZ,           // multiplier on PERMZ in a box
    Poro,            // multiplier on PORO in a box
    FaultMult,       // MULTFLT value of fault 'target'
    AquiferStrength, // multiplier on aquifer 'target' size (Fetkovich V and PI, Carter-Tracy k and h)
    KrwMax,          // multiplier on krw (SWOF column 2)
    KrwExp,          // curvature: krw' = krwMax * (krw/krwMax)^value
    KroExp,          // curvature of krow
    WaterVisc,       // multiplier on PVTW viscosity
    RockComp,        // multiplier on rock compressibility
    OWC,             // absolute OWC depth
    DatumPressure,   // absolute EQUIL datum pressure
    WellSkin,        // absolute skin of all completions of well 'target'
    WellPI,          // multiplier on the connection factors of well 'target'
    OilVisc,         // multiplier on oil viscosity (PVDO / PVTO, all PVT regions)
    GOC,             // absolute GOC depth
    Ntg,             // multiplier on NTG in a box
    MultZ,           // vertical transmissibility multiplier (MULTZ) in a box
    KrgExp,          // curvature of krg (SGOF column 2)
    KrogExp,         // curvature of krog (SGOF column 3)
    KrowMax,         // multiplier on the krow end point (SWOF column 3)
    PcMult,          // multiplier on oil-water capillary pressure (SWOF column 4)
    InjRate,         // multiplier on the injection rate targets of well 'target' (-1 = all injectors)
    ProdRate         // multiplier on the production rate targets of well 'target' (-1 = all producers)
};
constexpr int kNumParamKinds = int(ParamKind::ProdRate) + 1;
struct StudyParam {
    ParamKind kind = ParamKind::PermXY;
    std::string label;                 // empty -> generated
    int target = 0;                    // fault / aquifer / well index
    int i1 = -1, i2 = -1, j1 = -1, j2 = -1, k1 = -1, k2 = -1;   // box (zero-based, -1 = full range)
    double base = 1, low = 0.5, high = 2;
    bool log = true;                   // sample / search in log space
    bool active = true;
};

// THPRES: threshold pressure between two EQLNUM regions (flow across the boundary only once the potential difference
// exceeds it; value < 0 = computed from the initial equilibrium, the largest potential difference across the boundary)
struct ThresholdPressure { int r1 = 1, r2 = 2; double value = -1; };

// One keyword kept as deck text (see Case::extraKeywords)
enum class DeckSection : int { Runspec = 0, Grid, Edit, Props, Regions, Solution, Summary, Schedule };
struct ExtraKeyword {
    std::string name;                         // upper case, e.g. "GRUPTREE"
    DeckSection section = DeckSection::Schedule;
    int step = 0;                             // SCHEDULE only: report step it precedes (0 = before the first TSTEP/DATES)
    std::string body;                         // the records as deck text ('/'-terminated), without the keyword line
    std::string note;                         // why ResSim keeps it (import report); informational
    bool beforeControls = true;               // SCHEDULE: written before the well controls (WCONPROD, ...) of its step
};
const char* toString(DeckSection s);

struct Case {
    std::string title = "Untitled";
    UnitSystem units = UnitSystem::Field;
    // non-empty: the deck uses features the ResSim engine cannot model (why). The case was imported leniently
    // (defaults filled in for display); it can only be run by OPM Flow / ECLIPSE on the original deck.
    std::string externalOnly;
    Phases phases;
    GridInput grid;
    PropsInput props;
    EquilInput equil;                         // EQLNUM region 1
    std::vector<EquilInput> equilExtra;       // EQLNUM regions 2..N
    InitialState init;                        // enumerated initial state (overrides equilibration)
    std::vector<ThresholdPressure> thpres;    // EQLOPTS THPRES / THPRES
    std::vector<AquiferInput> aquifers;
    std::vector<WellInput> wells;
    std::vector<ReportStep> schedule;
    Tuning tuning;
    int startDay = 1, startMonth = 1, startYear = 2025;   // START
    std::vector<Observation> observations;
    std::vector<StudyParam> studyParams;      // sensitivity / history-matching parameters
    // Keywords the ResSim engine does not model (groups, VFP, MSW, ACTIONX, tracers, ...), kept verbatim so that the
    // case survives import -> edit -> export and runs unchanged in OPM Flow / ECLIPSE. Written back in deck order.
    std::vector<ExtraKeyword> extraKeywords;
    PolymerInput polymer;                     // used when phases.polymer
};

// Validates the case; returns empty string if OK, otherwise a human-readable message.
std::string validateCase(const Case& c);

// Ready-to-run demo cases.
Case makeDemoCase();                 // SPE1-like 3-phase live oil, water injector + producer (FIELD)
Case makeBuckleyLeverettCase(int n); // 1D waterflood for validation (FIELD)

// Corey-type saturation tables (SWOF: Sw from swc to 1; SGOF: Sg from 0 to 1 - swc), n + 1 rows
void coreySWOF(std::vector<SatRow>& out, double swc, double sor, double nw, double no, double krwMax, double kroMax, int n = 16);
void coreySGOF(std::vector<SatRow>& out, double swc, double sgc, double sorg, double ng, double nog, double krgMax, int n = 16);

const char* toString(Control c);
const char* toString(ObsQty q);
bool obsIsWell(ObsQty q);
const char* toString(WellType t);

} // namespace ressim
