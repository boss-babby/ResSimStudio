// Input data model of a simulation case (what the user edits in the UI).
// Mirrors the ECLIPSE sections: RUNSPEC / GRID / PROPS / SOLUTION / SCHEDULE.
#pragma once
#include "ressim/geomech.h"
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <map>
#include "ressim/eos.h"
#include "ressim/co2brine.h"

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
// C34 numerical aquifers (AQUNUM / AQUCON, OPM NumericalAquifer): AQUNUM grid cells become a 1-D chain of aquifer
// cells (pore volume A L phi, depth, optional initial pressure) connected in record order; AQUCON faces of reservoir
// cells connect to the FIRST cell of the aquifer with T = Tcon Tcell / (Tcon + Tcell) x multiplier
struct NumAquCell { int id = 1, i = 0, j = 0, k = 0; double area = 0, length = 0, poro = 0, perm = 0, depth = 0, p0 = -1; int pvt = 0, sat = 0; };
struct NumAquCon { int id = 1, i1 = 0, i2 = 0, j1 = 0, j2 = 0, k1 = 0, k2 = 0; char face = 'X'; bool minus = false; double mult = 1; int transOpt = 0; bool allowInternal = false; };

// Local grid refinement (ECLIPSE CARFIN / NXFIN.. / HXFIN.. / ENDFIN, OPM Flow CARFIN). A box of host cells I1..I2,
// J1..J2, K1..K2 (zero-based) of the global grid is replaced by NX x NY x NZ local cells. NXFIN (one entry per host
// column, sum NX) says into how many local cells each host column is split (empty = NX / (I2-I1+1) each), HXFIN gives the
// relative sizes of the NX local cells (empty = equal within each host cell). Local cells inherit every property of
// their host cell unless 'props' overrides it (arrays given between CARFIN and ENDFIN or after REFINE; local natural
// order, NaN = inherit). The refinement is solved in place, fully implicitly with the global grid (no local time
// stepping), as ECLIPSE does by default and OPM Flow always does.
struct LgrInput {
    std::string name = "LGR1";
    std::string parent = "GLOBAL";            // nested refinements are not modelled (kept for OPM Flow / ECLIPSE)
    int i1 = 0, i2 = 0, j1 = 0, j2 = 0, k1 = 0, k2 = 0;
    int nx = 1, ny = 1, nz = 1;
    int nwmax = 0;                            // CARFIN item 10 (kept for export)
    std::vector<int> nxfin, nyfin, nzfin;
    std::vector<double> hxfin, hyfin, hzfin;
    std::map<std::string, std::vector<double>> props;
    int size() const { return nx * ny * nz; }
    int index(int i, int j, int k) const { return i + nx * (j + ny * k); }
    // host column / row / layer (0-based, global) and the fraction range inside it of local cell li / lj / lk
    void split(int dir, std::vector<int>& host, std::vector<double>& f0, std::vector<double>& f1) const;
    bool valid(std::string* why = nullptr) const;
};

struct GridInput {
    int nx = 10, ny = 10, nz = 3;
    // per-cell arrays (size nx*ny*nz, natural ECLIPSE ordering i fastest, then j, then k)
    std::vector<double> dx, dy, dz, tops;     // tops: depth of the top face of each cell
    std::vector<double> poro, permx, permy, permz, ntg;
    std::vector<int> actnum;
    std::vector<double> multx, multy, multz;  // transmissibility multipliers (optional, size 0 or n)
    // C38 EDIT-section MULTIPLY on TRANX / TRANY / TRANZ (size 0 or n): unlike MULTX they act on the Cartesian neighbour
    // connection (i,j,k)-(i+1,j,k) only, not on fault juxtaposition NNCs across the same face (OPM Transmissibility)
    std::vector<double> tranMultX, tranMultY, tranMultZ;
    // corner-point geometry (ECLIPSE COORD/ZCORN). When zcorn is non-empty it replaces DX/DY/DZ/TOPS.
    std::vector<double> coord;                // 6*(nx+1)*(ny+1): pillar top (x,y,z) and bottom (x,y,z)
    std::vector<double> zcorn;                // 8*nx*ny*nz corner depths, ECLIPSE ordering
    std::vector<FaultInput> faults;
    std::vector<NNCInput> nnc;
    // C34 EDIT-section EDITNNC (trans = multiplier) / EDITNNCR (trans = replacement value) of non-neighbour connections
    std::vector<NNCInput> editnnc, editnncr;
    // region numbers (1-based, size 0 = every cell in region 1)
    std::vector<int> satnum, pvtnum, eqlnum, fipnum;
    std::vector<double> multpv;               // MULTPV pore-volume multipliers (empty = 1)
    // MULTREGT: transmissibility multipliers between regions of MULTNUM / FLUXNUM / OPERNUM, applied by the grid to every
    // connection (neighbours, fault NNCs, pinch-outs) like ECLIPSE: directions (item 4) restrict neighbour connections,
    // item 5 selects NNC / NONNC / ALL / NOAQUNNC, and the last record for a region pair wins (not cumulative)
    struct MultRegT { int r1 = -1, r2 = -1; double m = 1.0; std::string dirs = "XYZ"; int nnc = 0; char set = 'M'; };   // nnc: 0 ALL 1 NNC 2 NONNC 3 NOAQUNNC
    std::vector<MultRegT> multregt;
    std::vector<NumAquCell> numAquCells;      // C34 AQUNUM
    std::vector<NumAquCon> numAquCons;        // C34 AQUCON
    std::vector<int> multnum, fluxnum, opernum;   // region sets referenced by MULTREGT (empty = not given)
    std::vector<double> swatinit;             // SWATINIT initial water saturation (empty = none, < 0 = not set for the cell):
                                              // equilibration scales the water-oil capillary pressure per cell to honour it
    // saturation end-point scaling (ENDSCALE + SWL / SWCR / SWU / SGL / SGCR / SGU / SOWCR / SOGCR, ECLIPSE two-point
    // scaling of kr, Pc scaled between SWL and SWU / SGL and SGU). Empty = table end points; NaN or < 0 = table value.
    bool endscale = false;
    std::vector<double> swl, swcr, swu, sgl, sgcr, sgu, sowcr, sogcr;
    bool nonnc = false;                       // C38 NONNC (RUNSPEC): no non-neighbour connections across fault throws
    bool scalecrs = false;                    // SCALECRS YES: three-point kr scaling (critical, residual of the other phase, max)
    // hysteresis (SATOPTS HYSTER + EHYSTR): imbibition table regions and imbibition end points (ISWL ... ISOGCR)
    std::vector<int> imbnum;                  // 1-based, empty = SATNUM
    std::vector<int> rocknum;                 // C34 ROCKNUM (1-based)
    std::vector<int> krnumx, krnumy, krnumz;  // C34 SATOPTS DIRECT: relative permeability table per flow direction (1-based, empty = SATNUM)
    std::vector<double> iswl, iswcr, iswu, isgl, isgcr, isgu, isowcr, isogcr;
    std::vector<double> spoly;                // SPOLY initial polymer concentration (empty = 0)
    std::vector<int> plmixnum;                // PLMIXNUM polymer mixing regions (empty = 1)
    std::vector<int> surfnum;                 // C35 SURFNUM: saturation table of the miscible (low IFT) surfactant curves (1-based, empty = none)
    std::vector<double> sfoam;                // C35 SFOAM: initial foam concentration (empty = 0)
    std::vector<double> ssol;                 // C35 SSOL: initial solvent saturation (empty = 0)
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
    std::vector<double> disperc;              // DISPERC (OPM): mechanical dispersivity per cell (length), compositional engine
    // MINC (multiple interacting continua, Pruess & Narasimhan 1985) on top of DUALPORO: every matrix block is split into
    // 'minc' nested shells - the matrix cell is the outermost shell (touching the fractures), the inner shells are
    // appended after the grid cells. mincVol: volume fractions outer -> inner (empty = geometric 1 : 2 : 4 ...);
    // fracture sets in mincDims directions (1: Z, 2: X+Y, 3: X+Y+Z) with spacings LX / LY / LZ, else a block of size
    // L = sqrt(4 dims / sigma) from SIGMA / SIGMAV (Kazemi). Shell boundaries follow the proximity function
    // V(x) = 1 - prod(1 - 2x / l_d); T(j-1, j) = CDARCY * V * sum_d k_d dV_d/dx(x_{j-1}) / (y_j - y_{j-1}) * MULTMF
    int minc = 0;
    int mincDims = 3;
    std::vector<double> mincVol;
    int mincExpanded = 0;                     // internal: per-cell arrays already carry the shell cells
    bool nodppm = false;                      // NODPPM: fracture permeability is NOT multiplied by fracture porosity
    bool gravdr = false;                      // GRAVDR: gravity drainage term in the matrix-fracture flow (oil-water and gas-oil)
    std::vector<double> dzmtrxv;              // DZMTRX / DZMTRXV matrix block heights for GRAVDR (empty / 0 = cell thickness)
    std::vector<int> dpnum;                   // DPNUM per matrix cell (size NX*NY*NZ/2, empty = all dual): 0 = single porosity
    int dpnumApplied = 0;                     // internal: single porosity cells already folded into their fracture cell
    // ---- C15: hydraulic fractures in the black-oil grid (embedded discrete fractures, EDFM). Each fracture is a vertical
    // planar rectangle: trace (x0,y0)-(x1,y1) in grid coordinates (case length units; the same x / y as the 3-D view and
    // COORD for corner-point grids), depth interval zTop..zBot (< 0: whole grid), width (case length units), permeability
    // (mD; <= 0: cubic law w^2 / 12) and porosity (propped pack). Each fracture piece inside a host cell becomes an extra
    // cell (global index >= cells + LGR + MINC cells) connected to its host (T = CDARCY 2 A k_n / <d>, Li & Lee 2008 /
    // Moinfar 2014), to the neighbouring pieces of the same fracture, to crossing fractures (intersection), and to every
    // well completed in the host (EDFM well index). gamma > 0: stress-dependent conductivity k = k0 exp(-gamma (p_i - p))
    // (1 / pressure units). '-- RESSIM: HFRAC name x0 y0 x1 y1 zTop zBot width perm poro [gamma] [satnum]'
    struct HFrac {
        std::string name = "HF1";
        double x0 = 0, y0 = 0, x1 = 0, y1 = 0, zTop = -1, zBot = -1;
        double width = 0.003, perm = 0, poro = 0.5, gamma = 0;
        int satnum = 0;                       // 0: the host's SATNUM
        // C18 in-run growth (linear elastic fracture mechanics): kic > 0 enables it. The trace is the MAXIMUM extent; at the
        // start only the pieces within initLen (case length units) of the injection point (the pieces in the hosts of the
        // completed wells, else the trace centre) are open; closed pieces conduct like the host rock (fracture-fracture
        // transmissibility x k_host / k_f x closedMult). After every
        // accepted time step each wing's tip piece is tested: K_I = (p_tip - sigma(z)) sqrt(pi c), c = min(L_wing, H / 2)
        // (KGD while short, PKN once L > H/2), and the next piece opens when K_I >= K_IC.
        // '-- RESSIM: HFGROW name KIC[MPa m^0.5] sigmaMin[pressure] [sigmaGrad pressure/length] [initLen] [closedMult]'
        double kic = 0, sigma = 0, sigmaGrad = 0, initLen = 0, closedMult = 1;
        // C21 sub-grid tip (H2): continuous wing length, the tip piece is open by the fraction it covers (transmissibility
        // interpolated), and the length is solved implicitly in the time step (K_I = K_IC at the tip). Deck: '... SUBGRID'
        int subgrid = 0;
        // C25: the closure stress of the growth criterion is the smallest principal total stress of the finite-element
        // geomechanics (GEOMECH) at the piece's host cell instead of sigmaMin + sigmaGrad (z - zMid). Deck: '... FESTRESS'
        int feStress = 0;
        // C26: mechanical aperture of the open pieces from a 3-D displacement-discontinuity solve (ddm.h, half space,
        // opening under the net pressure p - sigma_min with contact) instead of the fixed EDFM width; the piece
        // conductivity follows the cubic law k_f w = w^3 / 12 relative to the input width (perm <= 0) or k_f w / w_input.
        // Deck: '... DDMWIDTH [E_GPa nu]' (defaults: the GEOMECH moduli, else 20 GPa / 0.25)
        int ddmWidth = 0;
        double ddmE = 0, ddmNu = 0.25;
        // C28 in-run phase-field growth (phasefield.h, 2-D plane strain, volume control): at the end of every report step the
        // crack volume per unit height V = eta * (cumulative injection of the fracture's wells, reservoir volume) / H drives
        // one load step of a persistent phase-field model of the fracture plane (E', G_c = K_IC^2 / E', closure stress
        // sigmaMin); its half-length opens the pieces of both wings (fluid flow in the crack = the EDFM pieces of the flow
        // model). Replaces the LEFM piece test for this fracture. Deck: '... PHASEFIELD [eta [h_m [E_GPa nu]]]'
        int phaseField = 0;
        double pfEta = 0.5, pfH = 0, pfE = 0, pfNu = 0.25;
        // C28 H2b sub-grid leak-off of the SUBGRID tip piece (exposed part of the host cell as its own storage). Deck: '... SUBGRID LEAK'
        int subLeak = 0;
        // C29 two-way phase-field coupling ('... PHASEFIELD [eta [h [E nu]]] TWOWAY'): the crack volume of each load step follows
        // the flow model's net pressure (V = C p_net, C = V / p_net of the previous phase-field step; eta only for the first
        // step) and the piece conductance follows the Sneddon aperture of the phase-field crack, (w / w_ref)^3
        int pfTwoWay = 0;
    };
    std::vector<HFrac> hfracs;
    struct HfPiece {                          // internal: one fracture piece (computed from the host geometry)
        int frac = 0, host = 0;               // fracture index, host global cell
        double xa = 0, ya = 0, xb = 0, yb = 0, zt = 0, zb = 0, ta = 0, tb = 0, dAvg = 1, kn = 0;
    };
    std::vector<HfPiece> hfPieces;
    int hfExpanded = 0;                       // internal: per-cell arrays already carry the fracture cells
    // ---- C15 unconventional rock: Langmuir adsorbed gas (shale gas / CBM; black-oil gas equation, instantaneous
    // equilibrium): V_ads = V_bulk (1 - phi) rho_grain V_L p / (P_L + p)  [V_L sm3/t (METRIC) or scf/ton (FIELD), P_L in
    // pressure units, rho_grain in g/cm3]; '-- RESSIM: LANGMUIR VL PL rhoGrain'
    double langVL = 0, langPL = 0, langRho = 2.6;
    // matrix stress-dependent permeability k = k0 exp(-gamma (p_i - p)) (1 / pressure units); '-- RESSIM: PERMSTRESS gamma'
    double permGamma = 0;
    std::vector<LgrInput> lgrs;               // local grid refinements (CARFIN)
    int lgrExpanded = 0;                      // internal: per-cell arrays already carry the local cells (expandLgrArrays)
    int lgrIndex(const std::string& name) const {
        for (size_t l = 0; l < lgrs.size(); ++l) if (lgrs[l].name == name) return int(l);
        return -1;
    }
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

// C34 rock compaction (ROCKCOMP / ROCKTAB / ROCK2D / ROCK2DTR / ROCKWNOD / OVERBURD / ROCKOPTS, OPM FlowProblem):
// pore volume and transmissibility multipliers vs effective pressure p - p_overburden(depth) (minus the initial
// pressure with ROCKOPTS STORE); IRREVERS uses min(p, lowest p seen); ROCK2D(TR) add the water-induced compaction
// axis dSw = max(Sw, highest Sw seen) - Sw_initial. Tables per ROCKNUM region (else ROCKOPTS item 3, else PVTNUM).
struct RockCompInput {
    bool on = false, irreversible = false, water = false, store = false;
    std::string regionArray = "PVTNUM";
    struct Tab1 { std::vector<double> p, pv, tm; };
    struct Tab2 { std::vector<double> p; std::vector<std::vector<double>> pv, tm; };   // rows per pressure: values per ROCKWNOD node
    std::vector<Tab1> tab;                   // ROCKTAB (with ROCKTABH: the deflation curve = first row of every elastic curve)
    bool hyster = false;                     // C41 ROCKCOMP HYSTER: elastic reload curves (ROCKTABH) above the turning pressure
    std::vector<std::vector<Tab1>> tabH;     // C41 ROCKTABH elastic curves per region (each starts on the deflation curve)
    std::vector<Tab2> tab2;                  // ROCK2D (pv) + ROCK2DTR (tm)
    std::vector<std::vector<double>> swNodes;   // ROCKWNOD per region
    std::vector<std::vector<double>> overburd;  // OVERBURD per region: depth, pressure pairs
    bool active() const { return on && (!tab.empty() || !tab2.empty()); }
};

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
    // relative permeability hysteresis (SATOPTS HYSTER + EHYSTR). model: EHYSTR item 2 (0 Carlson non-wetting /
    // drainage wetting, 1 Carlson / imbibition wetting, 2 Killough / drainage wetting, 3 Killough / imbibition wetting)
    struct Hysteresis { bool on = false; int model = 0; double curv = 0.1, killough = 0.1; std::string what = "BOTH"; } hyst;
    RockCompInput rockComp;          // C34 rock compaction tables
    // VAPPARS (oil vaporisation / gas re-dissolution limits) given before SCHEDULE; < 0 = not set
    double vap1 = -1, vap2 = -1;
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
struct AquiferConnection { int i1 = 0, i2 = 0, j1 = 0, j2 = 0, k1 = 0, k2 = 0; char face = 'X'; bool minus = false; double mult = 1.0; double coeff = -1; };   // coeff: AQUANCON influx coefficient (< 0: face area)
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
    std::vector<double> tdTab, pdTab;   // C34 Carter-Tracy influence table (AQUTAB; empty: built-in infinite-acting table 1)
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
    int defStep = -1;                   // C46: step a connection defined SHUT was defined at (opened at fromStep by WELOPEN);
                                        // WELPI scales it from then on, the exporter writes the SHUT definition there (-1 = fromStep)
    int satTable = 0;                   // COMPDAT item 6: saturation table for the connection (0 = cell's SATNUM); OPM only
    double dfactor = -1, r0 = -1;       // COMPDAT items 11, 13: D-factor and pressure-equivalent radius (<0 default); OPM only
    int lgr = -1;                       // COMPDATL: index into GridInput::lgrs (i, j, k1, k2 are then local cell indices)
    double mult = 1.0;                  // C44: cumulative WPIMULT factor of this open period (ECLIPSE / OPM: applied to the
                                        // connections existing at the WPIMULT step, reset when COMPDAT redefines the connection)
};

struct WellInput {
    std::string name = "P1";
    std::string group = "G";
    WellType type = WellType::Producer;
    int i = 0, j = 0;                  // wellhead location (zero-based)
    double refDepth = -1;              // BHP reference depth (<=0: top completion)
    double piMult = 1.0;               // PI multiplier on all connection factors of the whole run (UI / studies; deck WPIMULT
                                       // is applied per connection and step: Completion::mult)
    std::string prefPhase;             // C34 WELSPECS preferred phase (OIL, WATER, GAS, LIQ; empty: from the type)
    std::vector<Completion> completions;
    std::string specsTail;             // WELSPECS items 7.. (drainage radius, inflow eq., auto shut-in, crossflow, ...) as deck text
    int lgr = -1;                      // WELSPECL: the well lives in this local grid (i, j are local); -1 = global
};

// A control record (WCONPROD / WCONINJE / WCONHIST), applied from a report step onwards.
struct WellControl {
    std::string well;
    WellStatus status = WellStatus::Open;
    Control mode = Control::ORAT;
    double orat = -1, wrat = -1, grat = -1, lrat = -1, resv = -1;  // producer targets/limits (<0: none)
    double rate = -1;                  // injector surface rate target
    double bhp = -1;                   // BHP target (BHP mode) or limit (rate modes)
    bool bhpDefault = false;           // C44: bhp is the FBHPDEF default limit (written defaulted, FBHPDEF is kept in the deck)
    bool history = false;              // WCONHIST: rates are observed data
    bool reopen = false;               // C37: set by a well-level WELOPEN OPEN at this step (restarts a WCYCLE on-period even
                                       // when the well is already open - OPM REQUEST_OPEN_WELL event); not written to the deck
    double bhpHist = -1;               // observed BHP (WCONHIST item 10 / WCONINJH item 5): history data, NOT a limit
                                       // (ECLIPSE / OPM: the limit stays 1 atm / 100000 psi unless WELTARG; BHP mode -> target)
    int phase = 0;                     // injectors: injected phase for this period (0 = the well's type, 1 = water, 2 = gas)
    // deck items the ResSim engine does not use but OPM Flow / ECLIPSE do (kept for export)
    std::string deckMode;              // original control mode when not modelled (GRUP, THP, NONE, CRAT, RESV for injectors)
    double thp = -1;                   // THP target / limit
    int vfp = 0;                       // VFP table number (0 = none)
    double alq = -1;                   // artificial lift quantity
    // C31 UDA: targets given by a UDQ name (ORAT WRAT GRAT LRAT RESV BHP THP RATE -> UDQ); the simulator sets the value
    // each time step from the UDQ (udq.h), the deck writer writes the name back
    std::vector<std::pair<std::string, std::string>> uda;
    const std::string* udaOf(const char* target) const { for (auto& u : uda) if (u.first == target) return &u.second; return nullptr; }
};

// WPOLYMER: polymer concentration of the injected water (well name or pattern; stays in force until changed)
struct WellPolymer { std::string well; double conc = 0; };

// ---------------------------------------------------------------- groups (GRUPTREE / GCONPROD / GCONINJE / GCONPRI)
// Group targets are allocated to the wells of the group's subtree at the start of each time step (explicit), in
// proportion to the well potentials: wells in GRUP mode always share the target, individually controlled wells are
// switched to group control when the group exceeds its target (ECLIPSE default WGRUPCON YES).
struct GroupProdControl {
    std::string group;
    Control mode = Control::BHP;             // BHP = NONE (limits only)
    double orat = -1, wrat = -1, grat = -1, lrat = -1, resv = -1;
    std::string action = "NONE";             // procedure on exceeding the limits (RATE: the violated phase controls; NONE: no action)
    std::string actW, actG, actL;            // C35 items 11-13: procedure for the water / gas / liquid limits (empty = item 7)
    const std::string& actionFor(Control q) const {
        const std::string& a = q == Control::WRAT ? actW : q == Control::GRAT ? actG : q == Control::LRAT ? actL : action;
        return a.empty() ? action : a;
    }
    bool priority = false;                   // from GCONPRI (priority ranking approximated by potentials)
    double guideRate = -1;                   // C34 GCONPROD item 9: group guide rate (< 0: sum of the well guide rates)
    std::string guideDef;                    // C34 item 10: OIL WAT GAS LIQ RES (with guideRate), FORM (GUIDERAT on the group potentials), POTN
    // C48 UDA: rates given by a UDQ name (ORAT WRAT GRAT LRAT RESV -> FU_/GU_ name), set every time step from the UDQ
    std::vector<std::pair<std::string, std::string>> uda;
};
struct GroupInjControl {
    std::string group;
    int phase = 1;                           // 1 water, 2 gas
    std::string mode = "NONE";               // NONE RATE RESV REIN VREP
    double rate = -1, resv = -1, rein = -1, vrep = -1;
    std::string reinGroup, vrepGroup;        // default: the group itself
    std::vector<std::pair<std::string, std::string>> uda;   // C48 UDA: RATE RESV REIN VREP -> FU_/GU_ name
};
// VFPPROD / VFPINJ lift tables (BHP at the datum depth vs flow, THP, water / gas fractions and ALQ)
struct VfpTable {
    int id = 0;
    bool inj = false;
    double datum = 0;
    std::string flo = "OIL", wfr = "WCT", gfr = "GOR", alq;
    std::vector<double> floV, thpV, wfrV, gfrV, alqV;
    std::vector<double> bhp;                 // prod: ((((it*nw+iw)*ng+ig)*na+ia)*nf+if), inj: it*nf+if
    double eval(double flo, double thp, double wfr, double gfr, double alq) const;
};

// GUIDERAT: general guide-rate formula GR = POT^A / (B + C R1^D + E R2^F) for the nominated phase
struct GuideRateSpec {
    bool set = false;                  // GUIDERAT given at the start of this step
    std::string phase = "NONE";        // OIL, LIQ, GAS, RES, COMB, NONE (NONE: well potentials)
    double minInterval = 0, A = 0, B = 0, C = 0, D = 0, E = 0, F = 0;
    bool allowIncrease = true;
    double damping = 1.0;
};
struct ReportStep {
    GuideRateSpec guiderat;            // GUIDERAT in force from this step (set = true)
    std::vector<std::pair<std::string, double>> welpi;   // C34 WELPI (well / pattern, target PI of the preferred phase)
    double length = 30;                // days
    std::vector<WellControl> controls; // applied at the start of this step
    std::vector<WellPolymer> wpolymer; // WPOLYMER records applied at the start of this step
    // DRSDT set at the start of this step: < -1.5 unchanged, < 0 off, >= 0 maximum Rs increase per day
    double drsdt = -2;
    bool drsdtFree = false;            // DRSDT 'FREE': only cells with free gas are limited (default 'ALL')
    double vap1 = -1, vap2 = -1;       // VAPPARS set at the start of this step (< 0 unchanged)
    // DRSDT / DRVDT and VAPPARS exclude each other (OPM OilVaporizationProperties: the last keyword sets the type, so
    // VAPPARS switches DRSDT off and DRSDT switches VAPPARS off). Order when both appear in the same step:
    bool vapAfterDrsdt = false;
    std::vector<GroupProdControl> gprod;  // GCONPROD / GCONPRI records applied at the start of this step
    std::vector<GroupInjControl> ginj;    // GCONINJE records
    std::vector<std::pair<std::string, std::string>> injStreams;   // compositional: WINJGAS (well, WELLSTRE stream) from this step
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
    bool relaxOscillations = true;     // damp Newton updates when residuals oscillate / stagnate (OPM-style)
    bool blockSolver = true;           // Schur-eliminate wells and solve the block (BSR) system with BILU0 / block CPR
    bool mixedPrecision = false;       // C16: BILU0 (standalone / CPR second stage) stored and applied in float32
    int leanMemory = -1;               // C45: low-memory linear algebra - BILU0 factorised straight into float32 storage
                                       // (no double copy of the matrix); -1 = auto (>= 500 000 active cells), 0 off, 1 on
    int gpu = 0;                       // C16: 1 = Krylov solve of the block system on the GPU backend (if available)
    int amgReuse = 10;                 // C26: keep the CPR-AMG hierarchy across the time steps of one report step for up to N
                                       // steps (0 = rebuild every time step); rebuilt early when the linear iterations rise 50 %
                                       // above the fresh-build count and always at a new report step (bitwise restarts). Egg -3 %
    // C27: target-change time-step control (CMG NORM PRESS / SATUR, Aziz-Settari): the next step is scaled by
    // (1 + w) eta / (d + w eta), w = 1, d = largest change of the step, eta = target; 0 = off (iteration-count rule only)
    double dsTarget = 0;               // target max saturation change per time step (e.g. 0.05)
    double dpTarget = 0;               // target max pressure change per time step (case pressure units, e.g. 50 bar / 700 psi)
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
    ProdRate,        // multiplier on the production rate targets of well 'target' (-1 = all producers)
    // C40 geomodel parameters: re-simulate the property model(s) of the case (Case::propModels; target = index, -1 = all)
    GeoRange,        // multiplier on the variogram ranges (major, minor, vertical)
    GeoSeed          // realisation number (rounded to an integer; the random seed of the simulation)
};
constexpr int kNumParamKinds = int(ParamKind::GeoSeed) + 1;
struct StudyParam {
    ParamKind kind = ParamKind::PermXY;
    std::string label;                 // empty -> generated
    int target = 0;                    // fault / aquifer / well index
    int i1 = -1, i2 = -1, j1 = -1, j2 = -1, k1 = -1, k2 = -1;   // box (zero-based, -1 = full range)
    double base = 1, low = 0.5, high = 2;
    bool log = true;                   // sample / search in log space
    bool active = true;
};

// C40: property model (geostatistics recipe) kept with the case, so that Studies can re-simulate it with other
// variogram ranges / seeds (ParamKind::GeoRange, GeoSeed). Written by the Model geostatistics panel and by the Geomodel
// workspace (Send to Model: blocked well data + Petrophysics variogram + poro-perm transform). Persisted as RESSIM lines.
struct PropModelInput {
    int property = 0;                  // 0 PORO, 1 PERMX (log10; PERMY = PERMX, PERMZ = kvkh PERMX), 2 NTG, 3 FACIES (C41)
    int method = 1;                    // 0 kriging, 1 SGS
    int model = 0;                     // variogram: 0 spherical, 1 exponential, 2 Gaussian
    double nugget = 0, sill = 1, rangeMax = 1000, rangeMin = 1000, rangeVert = 20, azimuthDeg = 0;
    unsigned seed = 1;
    int maxNbr = 16;
    bool fromWells = false;            // also condition on the current values at the completed cells
    std::vector<std::array<double, 4>> data;   // hard data x, y, z (depth), value (property units)
    double minV = -1e300, maxV = 1e300;
    bool permFromPoro = false;         // PORO models: PERMX = PERMY = 10^(ppA + ppB phi), PERMZ = kvkh PERMX
    double ppA = 0, ppB = 0, kvkh = 0.1;
    std::string source;                // where it came from (shown in Studies)
    // C41 FACIES models (property 3): truncated Gaussian simulation - one standard Gaussian SGS field (variogram
    // normalised to a unit sill) truncated at the thresholds of the cumulative target proportions, conditioned on
    // the facies codes (data value = code 1..K; well cells: the current SATNUM when facSatnum); each cell then gets
    // its facies' PORO / PERMX (PERMY = PERMX, PERMZ = kvkh PERMX) and, with facSatnum, SATNUM = facies code.
    std::vector<double> facProp, facPoro, facPerm;   // per facies (code k = index + 1); poro / perm < 0 = keep
    bool facSatnum = false;
};

// THPRES: threshold pressure between two EQLNUM regions (flow across the boundary only once the potential difference
// exceeds it; value < 0 = computed from the initial equilibrium, the largest potential difference across the boundary)
struct ThresholdPressure { int r1 = 1, r2 = 2; double value = -1; };

// Compositional model (E300 COMPS decks or built in the UI) and the dedicated CO2 storage model (CO2STORE). Run by the
// ResSim compositional core (EOS / CO2-brine FIM); exported as E300 / OPM keywords.
struct CompStream { std::string name; std::vector<double> z; };   // WELLSTRE: injection composition (mole fractions)
struct CompositionalInput {
    bool enabled = false;          // the case is compositional (black-oil PVT tables are placeholders)
    bool co2store = false;         // CO2STORE: CO2 + H2O with the Spycher-Pruess / Duan-Sun brine model
    // true: the exporter writes the compositional keywords from this data (UI edits, UI-built cases); false: the
    // imported keywords are written back verbatim (Case::extraKeywords) so the original deck is untouched
    bool structured = false;
    eos::CompFluid fluid;          // EOS components, BIC, RTEMP / TEMPVD, ZMFVD / ZI, STCOND, SALINITY (fluid.salinity)
    std::vector<CompStream> streams;
    std::vector<std::pair<double, double>> separator;   // FIELDSEP stages (T K, P Pa); last stage liquid -> stock tank
    // SEPCOND separator trains (T K, P Pa per stage) and the WSEPCOND assignment of wells to trains
    struct SepTrain { std::string name; std::vector<std::pair<double, double>> stages; };
    std::vector<SepTrain> sepTrains;
    std::map<std::string, std::string> wellSep;         // well -> train name (others use FIELDSEP)
    // molecular diffusion in EOS mode (DIFFUSE + DIFFC-style per-component coefficients, m2/s; empty = off)
    std::vector<double> diffOil, diffGas;
    // CO2STORE trapping mechanisms and interactions (SI units): diffusion (DIFFUSE / DIFFCWAT / DIFFCGAS), salt
    // precipitation (PRECSALT / VAPWAT / SALTSOL / PERMFACT), mineral trapping and Land sgtMax ('-- RESSIM:' lines).
    // Residual trapping uses SATOPTS HYSTER / EHYSTR (props.hyst) and IMBNUM (grid.imbnum).
    co2::TrappingOptions trap;
    bool waterVap = false;         // EOS mode: water vaporisation into the gas phase (VAPWAT / '-- RESSIM: WATERVAP')
    int nc() const { return fluid.eosRegions.empty() ? 0 : fluid.eosRegions[0].nc(); }
};
// E300 THERMAL model (milestone C7): steam / hot-water processes on the ResSim thermal core (K-value phase behaviour,
// water / steam from IAPWS-IF97). Values are kept in DECK units so the keywords round-trip unchanged:
//   METRIC: bar, deg C (TEMPI, TEMPVD, viscosity tables, WINJTEMP, STCOND), K (TCRIT, KVCR D / E), kg/m3, 1/bar, 1/K,
//           kJ/kg/K (SPECH*), kJ/kg (HEATVAP*), kJ/m3/K (HEATCR), kJ/m/day/K (THCONR), cP
//   FIELD:  psia, deg F, deg R, lb/ft3, 1/psi, 1/R, Btu/lb/R, Btu/lb, Btu/ft3/R, Btu/ft/day/R, cP
struct ThermalComponent {
    std::string name = "OIL";
    double mw = 400;                        // MW (g/mol)
    double tcrit = 1e5, pcrit = 1;          // TCRIT (TCRIT >= 1e4 or CVTYPE DEAD: non-volatile), PCRIT
    bool dead = true;                       // CVTYPE DEAD
    double kvcr[5] = {0, 0, 0, 0, 0};       // KVCR A B C D E: K = (A + B / p + C p) exp(-D / (T - E))
    double dref = 950, cref = 0, thermex1 = 0;          // DREF, CREF, THERMEX1
    double specha = 0, spechb = 0, spechg = 0, spechh = 0;   // SPECHA / B (oil), SPECHG / H (gas)
    double heatvaps = 0, heatvap = 0, heatvape = 0.38;       // HEATVAPS, HEATVAP, HEATVAPE
    double zfactor = 0.96, zfact1 = 0;      // ZFACTOR, ZFACT1 (ZFACT1 kept for export; not used by the core)
    double gvA = 0, gvB = 0;                // GASVISCF: mu_g = A T^B (cP; T in K | R); 0 = GASVISCT / default
    // component type (C8): 0 K-value, 1 non-condensable gas (CVTYPE GAS), 2 solid (coke / char / coal / ash; density
    // DREF / CREF / THERMEX1, heat capacity SPECHA / SPECHB; written as CVTYPE DEAD + '-- RESSIM: CTYPE name SOLID')
    int type = 0;
    // C10 (SI; '-- RESSIM:' lines AQSOL / CRIT / FORMULA): aqueous solubility (thermal::Component::aqModel ...),
    // Peng-Robinson critical constants (tcG <= 0: gas library value by name when GASEOS PR), element formula
    int aqModel = 0;
    double hA = 0, hB = 0, hC = 0, kh298 = 0, hdT = 0, vInf = 0, ksalt = 0;
    double tcG = 0, pcG = 0, omegaG = 0;
    std::string formula;
    double trackEq = 0;
};
// Chemical reaction of the thermal model (C8; SI units: A in the concentration units mol/m3 of bulk volume or kPa for
// the O2PP component, Ea / dH J/mol; written as '-- RESSIM: REACTION' lines - see thermal::Reaction)
struct ThermalReaction {
    std::string name = "R1";
    std::vector<double> sto, prod, order;   // nc + 1 (water last); order empty = sto
    std::vector<int> phase;                 // nc + 1; empty = auto
    double A = 0, Ea = 0, dH = 0, tMin = 0;
    int o2pp = -1;
    bool reversible = false;
    double keqA = 0, keqB = 0;
};
struct ThermalWinj {                        // WINJTEMP record in force from report step 'step'
    std::string well;
    int step = 0;
    double quality = 1, temp = -1e9, pres = 0, enth = 0;   // defaulted temperature: < -1e8
};
struct ThermalTrap {                        // '-- RESSIM: STEAMTRAP well dT [steamMax]' (SCHEDULE) in force from report step 'step'
    std::string well;
    int step = 0;
    double subcool = 0;                     // degC | degF temperature difference; <= 0: constraint off
    double steamMax = 0;                    // C17: max produced steam, STB/d | sm3/d cold-water equivalent; <= 0: off
};
struct ThermalFracture {                    // '-- RESSIM: FRACTURE name x0 y0 x1 y1 zTop zBot aperture perm [MECH kn0 vm mu ks dil]' (SI)
    std::string name = "F";
    double x0 = 0, y0 = 0, x1 = 0, y1 = 0, zTop = -1, zBot = -1, aperture = 1e-3, perm = 0;
    bool mech = false;
    double kn0 = 20e9, vm = 0.4e-3, mu = 0.6, ks = 0, dilationDeg = 5;
};
struct ThermalFracStress {                  // '-- RESSIM: FRACSTRESS sHmaxGrad sHminGrad azimuth E nu alphaT' (Pa/m, deg, Pa, -, 1/K)
    double sHmaxGrad = 0, sHminGrad = 0, azimuthDeg = 0, E = 30e9, nu = 0.25, alphaT = 0;
};
struct ThermalSorb {                        // '-- RESSIM: SORB comp capacity b tRef dH tau [AUTO|GAS|AQ] [adsRho] [site] [EQ|NOEQ]' (SI)
    std::string comp, site = "SITES";
    double capacity = 0, b = 0, tRef = 298.15, dH = 0, tau = 3600, adsRho = 0;
    int phase = -1;                         // -1 auto, 0 aqueous, 2 gas
    bool initEq = true;
};
struct ThermalInput {
    bool enabled = false;                   // THERMAL case (black-oil PVT tables are placeholders)
    bool structured = false;                // export the thermal keywords from this data (false: verbatim copies)
    bool deadOil = false;                   // DEADOIL (RUNSPEC)
    std::vector<ThermalComponent> comps;    // hydrocarbon components (water is implicit)
    std::vector<std::vector<double>> oilViscT;   // OILVISCT rows: T, mu_1 .. mu_nc
    std::vector<std::vector<double>> gasViscT;   // GASVISCT rows: T, mu_1 .. mu_nc
    std::vector<std::pair<double, double>> watViscT;   // WATVISCT (T, mu); empty: IAPWS-2008 (written as a table)
    double stdT = -1e9, stdP = -1;          // STCOND (defaults 15.56 C / 60 F, 1.01325 bar / 14.696 psia)
    // rock (GRID arrays: size 0 = default, 1 = uniform, else NX*NY*NZ)
    std::vector<double> heatcr, heatcrt, thconr, thconsf;
    // initial state: TEMPI (cells) or TEMPVD (depth, T); oil / gas compositions ZI | XMF / YMF (cells x nc)
    std::vector<double> tempi;
    std::vector<std::pair<double, double>> tempvd;
    std::vector<double> zi;                 // initial oil-phase composition (nc)
    std::vector<double> xmf, ymf;           // per cell, component-major (nc blocks of NX*NY*NZ)
    std::vector<ThermalWinj> winj;
    std::vector<ThermalTrap> traps;         // producer steam-trap (subcool) constraints (STARS OPERATE STEAMTRAP)
    // heat loss to over / underburden ('-- RESSIM: HEATLOSS on cv k' lines; ECLIPSE ROCKCON / ROCKPROP are not written):
    // cv kJ/m3/K | Btu/ft3/R, k kJ/m/day/K | Btu/ft/day/R
    bool heatLoss = false;
    double obCv = 0, obK = 0;
    // ---- C8: chemistry, solids, gas injection, permeability alteration, geomechanics (SI units)
    std::vector<ThermalReaction> reactions;
    std::vector<double> solidConc;          // initial solid concentrations, mol per m3 of pore volume (nc; solids only)
    std::vector<CompStream> streams;        // WELLSTRE (nc mole fractions) for gas injectors (WINJGAS per report step)
    int permLaw = 0;                        // 0 off, 1 Kozeny-Carman, 2 exponential (thermal::Case::permLaw)
    double permC = 20, permMultMax = 1e5;
    geomech::Input geo;                     // '-- RESSIM: GEOMECH' (OPM MECH / YMODULE / PRATIO / BIOTCOEF / THERMEXR imported)
    // ---- C10: dissolved gases, PR gas, gas hysteresis, tracked component ('-- RESSIM:' SALINITY / GASEOS / HYSTGAS / TRACK)
    double salinity = 0;                    // mol/kg NaCl (SALINITY keyword of CO2STORE decks is imported too)
    bool gasPR = false;
    bool hystGas = false;
    double sgtMax = 0.3;
    std::string track;                      // tracked component name ("" = CO2 when present)
    double structTol = 0.01;
    std::vector<ThermalSorb> sorb;         // C12 Langmuir sorption (thermal::addSorption)
    std::vector<ThermalFracture> fractures; // C12 embedded discrete fractures (thermal::addFractures; Cartesian boxes)
    ThermalFracStress fracStress;
    int nc() const { return int(comps.size()); }
};

// Builds the two-component (CO2, H2O) fluid of a CO2STORE case at temperature T (K) and salinity m (mol/kg)
eos::CompFluid co2StoreFluid(double T, double salinity);
// SPE5 (Killough & Kossack 1987) six-component Peng-Robinson fluid with its initial oil composition (ZI), 160 F
eos::CompFluid spe5CompFluid();
// SPE5 fluid + CO2 as the 7th EOS component (PR BICs with CO2 from the literature), CO2SOL on (CO2 dissolves in brine)
eos::CompFluid spe5Co2CompFluid(double salinity = 1.0);

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

// C24: fault mechanics + induced seismicity inside the flow run ('-- RESSIM: STRESSFIELD / FAULTMECH / FAULTRSF', see faultmech.h)
struct FaultMechInput {
    // in-situ total stress (compression +), case units: pressure / length; svGrad <= 0 = no stress field (feature off)
    double svGrad = 0, sv0 = 0;               // Sv = sv0 + svGrad * depth
    double shRatio = 0.7, sHRatio = 0.9;      // Shmin / Sv, SHmax / Sv (total stresses)
    double azimuth = 0;                       // SHmax direction, degrees from grid +x towards grid +y
    double stressPath = 0.5;                  // horizontal total-stress change per unit pore-pressure change (poroelastic)
    // C25 stress source: 0 analytic field above; 1 FE = the total stress tensor of the finite-element geomechanics (black oil /
    // compositional: GEOMECH sequential coupling, thermal: the THM mechanics) per face cell. Thermal runs with the analytic
    // field add the thermo-elastic horizontal stress dS_h = -E alpha / (1 - nu) dT (thE in GPa, thAlpha 1/K, 0 = off).
    // Deck: 'STRESSFIELD ... [FE] [thE thAlpha thNu]'
    int source = 0;
    double thE = 0, thAlpha = 0, thNu = 0.25;
    struct Row {
        std::string name;                     // FAULTS name
        double mu = 0.6, coh = 0;             // Coulomb friction, cohesion (pressure units)
        double aSig = 0;                      // Dieterich A sigma (pressure units); 0 = 0.003 x initial sigma'_n per face
        double ta = 100;                      // aftershock time (years) -> background stressing rate A sigma / ta
        double r0 = 0.1;                      // background rate of events >= mMin per km2 of fault per year
        double bValue = 1.0, mMin = 0.0;      // Gutenberg-Richter
        double permFactor = 1.0;              // two-way coupling: transmissibility factor of a reactivated face (1 = one-way)
        double permRamp = 0.05;               // C49: CFF over which the factor ramps in, as a fraction of sigma'_n (item 10)
        int rsfLayer = 0;                     // FAULTRSF: 1-based layer of the quasi-dynamic rate-and-state trace (0 = off)
        double rsfA = 0.010, rsfB = 0.015, rsfDcMm = 1.0;
    };
    std::vector<Row> faults;
    // C26 elastic stress transfer between the fault faces (3-D DDM, ddm.h) and quasi-static frictional slip: faces whose
    // CFF reaches 0 slip along their initial shear direction until CFF = 0 again (static friction mu, or the dynamic
    // friction mu - dMu once failed: weakening -> cascading slip events); the slip loads / unloads the other faces.
    // Deck: '-- RESSIM: FAULTSLIP E_GPa nu dMu [slipRef_mm] [cutoff] [FULLSPACE]'
    bool slipOn = false;
    double slipE = 20, slipNu = 0.25, slipDMu = 0;
    double slipRefMm = 10;                    // two-way coupling: slip at which a face reaches permFactor (slip-induced permeability)
    double slipCutoff = 25;                   // interaction cut-off distance in mean element sizes (<= 0: all pairs)
    int slipHalfSpace = 1;                    // 1: free surface at depth 0 (half space), 0: full space
    bool on() const { return svGrad > 0 && !faults.empty(); }
};
// C24: G1-c fracture contact coupling inside the main run ('-- RESSIM: GEOCOUPLE sxx syy E_GPa nu friction layer')
struct GeoCoupleSpec {
    bool on = false;
    double sxx = 0, syy = 0;                  // far-field total stresses, case pressure units, compression +
    double eGPa = 20, nu = 0.25, friction = 0.6;
    int layer = 1;                            // 1-based grid layer of the contact slice
};

struct Case {
    std::string title = "Untitled";
    UnitSystem units = UnitSystem::Field;
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
    // groups (simulation data; the deck keywords are kept verbatim for export)
    std::vector<std::pair<std::string, std::string>> groupTree;   // (child, parent) from GRUPTREE
    std::vector<VfpTable> vfp;                                    // VFPPROD / VFPINJ tables
    Tuning tuning;
    int startDay = 1, startMonth = 1, startYear = 2025;   // START
    std::vector<Observation> observations;
    std::vector<StudyParam> studyParams;      // sensitivity / history-matching parameters
    std::vector<PropModelInput> propModels;   // C40 geostatistical property models (re-simulated by GeoRange / GeoSeed)
    // Keywords the ResSim engine does not model (groups, VFP, MSW, ACTIONX, tracers, ...), kept verbatim so that the
    // case survives import -> edit -> export and runs unchanged in OPM Flow / ECLIPSE. Written back in deck order.
    std::vector<ExtraKeyword> extraKeywords;
    // C46: MULTFLT in the SCHEDULE section - cumulative fault transmissibility factors applied from the start of report
    // step 'step' (OPM / ECLIPSE semantics). Parsed from the verbatim ExtraKeyword, which is what the exporter writes.
    struct SchedMultflt { int step = 0; std::string fault; double mult = 1.0; };
    std::vector<SchedMultflt> schedMultflt;
    // imported decks: features the ResSim engine cannot model (gas-water systems, missing PVT/SCAL, binary grid
    // files, ...). The case still holds the grid / wells / schedule for display (placeholders elsewhere); it can only be
    // run by OPM Flow / ECLIPSE on the original deck. validateCase() rejects it for the own engine.
    std::vector<std::string> engineUnsupported;
    PolymerInput polymer;                     // used when phases.polymer
    CompositionalInput comp;                  // compositional / CO2STORE model (comp.enabled)
    ThermalInput thermal;                     // E300 THERMAL model (thermal.enabled)
    FaultMechInput faultMech;                 // C24 fault mechanics + seismicity in the run (black oil; C25 + thermal / comp)
    // C27 passive tracers of the ResSim black-oil engine: the ECLIPSE keywords TRACER (PROPS) / WTRACER (SCHEDULE) stay
    // verbatim extra keywords (exact round trip, OPM Flow runs them too); ResSim extensions per tracer (TOUGH EOS1 / EOS7R
    // style): '-- RESSIM: TRACERPROP name halfLife_days kd [init]' - first-order decay, linear equilibrium adsorption
    // (adsorbed amount per pore volume = kd * c), uniform initial concentration
    // C28: kPart = partition coefficient into the partner phase (oil for WAT / GAS tracers, water for OIL tracers; ratio of
    // the concentrations per reservoir volume, PITT), parent = decay-chain parent tracer (its decay feeds this tracer, exact
    // Bateman for decay alone), alphaL = longitudinal dispersivity [length], dm = molecular diffusion [length^2 / day]
    struct TracerProp { std::string name; double halfLife = 0, kd = 0, init = 0, kPart = 0, alphaL = 0, dm = 0; std::string parent; };
    std::vector<TracerProp> tracerProps;
    // C29 discretised-wellbore friction + inflow-control devices (FlexWell / T2Well-lite, black oil): per well, the wellbore
    // pressure profile adds the friction drop of the mixture flowing up (producers) / down (injectors) the tubing between the
    // connections (Haaland friction factor, along-hole lengths from the connection-cell centres in completion order) and a
    // per-connection ICD drop dp = K * rho_mix / rho_cal * q^2 (q = connection reservoir rate, case units per day). Explicit
    // (lagged with the hydrostatic heads, under-relaxed). Deck: '-- RESSIM: WELLBORE well ID roughness [icdK [rhoCal]]'
    // (ID, roughness in case length units; K in case pressure / (case reservoir rate)^2)
    struct WellboreSpec { std::string well; double id = 0.1, rough = 1.5e-5, icdK = 0, rhoCal = 1000; };
    std::vector<WellboreSpec> wellbores;
    // C37 multisegment-well options (ResSim extension, '-- RESSIM: MSWOPTS flash storage', 1 / 0): tubing flash of the segment
    // fluid at the segment pressure, segment storage (transient wellbore content); both on by default (ECLIPSE-like)
    // C45: third token = producer cross-flow in multisegment wells (experimental, off by default: with the explicit
    // start-of-step heads it can oscillate - see ENGINES_AND_MATH §9.v; RESSIM_MSW_XFLOW=1 switches it on as well)
    // C49: fourth token = connection heads of multisegment wells: 1 converged (the segment solution is re-evaluated from the
    // end-of-step state and the step re-solved until the connection heads stop moving - a Picard coupling, not a full
    // Jacobian), 0 explicit (start of step, as before C49), -1 default = converged (no case of the 24 OPM MSW references was
    // worse; TEST_WSEGSICD 93.4 -> 34.0 %)
    struct MswOptions { bool flash = true, storage = true, xflow = false; int heads = -1; };
    MswOptions mswOpt;
    // C37 PYACTION: directory of the imported deck (script files are relative to it; set by the deck import, not
    // written to decks) and the Python interpreter (empty: RESSIM_PYTHON, python3 / python on PATH)
    std::string deckDir;
    std::string pythonExe;
    // C29 salt precipitation / dry-out (TOUGH EWASG-type, black oil): dissolved NaCl is a water-phase species (internal tracer
    // NACL, kg per m3 of water); above the solubility it precipitates as solid S_s (fraction of the pore volume); gas entering
    // a water-bearing cell undersaturated in water vapour evaporates water (local equilibrium, vapour capacity omega =
    // p_sat(T) M_w / (R T rho_w), lagged sink) which concentrates the brine; porosity phi (1 - S_s), permeability by Verma &
    // Pruess (1988) tubes-in-series k/k0 = theta^2 (1 - G + G / w^2) / (1 - G + G (theta / (theta + w - 1))^2),
    // theta = (1 - S_s - phiR) / (1 - phiR), w = 1 + (1 / G) / (1 / phiR - 1).
    // Deck: '-- RESSIM: SALT c0 [cSol [rhoSalt [phiR [Gamma [T_C [cInj]]]]]]' (kg/m3 water, kg/m3 water, kg/m3, -, -, degC, kg/m3)
    struct SaltInput { bool on = false; double c0 = 0, cSol = 317, rhoSalt = 2165, phiR = 0.8, gamma = 0.8, tempC = 60, cInj = 0; };
    SaltInput salt;
    // C30 one-way field output of the in-run hydraulic fractures (hf_field.h): DDM opening under the net pressure, stress /
    // displacement / strain change of the matrix around the open pieces (display only - no feedback to the flow model).
    // on: -1 auto (fractures with HFGROW / DDMWIDTH), 0 off, 1 on. E (GPa) / nu: 0 = the DDMWIDTH / GEOMECH / 20 GPa
    // defaults; radius: region around the fractures as a multiple of the largest fracture extent; budget: millions of
    // triangle evaluations per report (cells x pieces cap). Deck: '-- RESSIM: HFFIELD ON|OFF [E_GPa nu] [radius] [budget]'
    struct HfFieldInput { int on = -1; double E = 0, nu = 0.25, radius = 1.0, budget = 6.0; };
    HfFieldInput hfField;
    // C25 sequential finite-element geomechanics of black-oil / compositional runs ('-- RESSIM: GEOMECH / MECHBC / GEOPLAST'
    // of a non-thermal deck; thermal decks use thermal.geo). Cartesian and corner-point grids (isoparametric hexahedra).
    // coupling 0: one-way (stress / strain / failure output, FE stress for HFGROW FESTRESS and STRESSFIELD FE);
    // 1: two-way explicit (after every report step: pore volume x phi_FE / phi_flow, matrix transmissibility x k-multiplier)
    geomech::Input geoSeq;
    GeoCoupleSpec geoCouple;                  // C24 fracture contact coupling in the run (black oil + HFRAC)
};

// Validates the case; returns empty string if OK, otherwise a human-readable message.
std::string validateCase(const Case& c);
std::string validateCompositional(const Case& c);
std::string validateThermal(const Case& c);         // c.thermal checks (called by validateCase when c.thermal.enabled)

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
