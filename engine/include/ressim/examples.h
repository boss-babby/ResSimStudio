// Built-in example library: ready-to-run cases from a 1D Buckley-Leverett check to faulted corner-point and
// channelised models. Every example is generated in code (deterministic), so it can also be exported as an
// ECLIPSE deck (see tools/export_examples or File > Examples > Export in the app).
#pragma once
#include "ressim/hydrofrac.h"
#include "ressim/case.h"
#include <string>
#include <vector>

namespace ressim {

enum class ExampleLevel { Basic, Intermediate, Advanced };

struct ExampleInfo {
    std::string id;            // stable identifier, also the exported deck name (upper-cased)
    std::string title;
    ExampleLevel level = ExampleLevel::Basic;
    std::string summary;       // one line
    std::string details;       // what it demonstrates / what to look at
    std::string grid;          // e.g. "25 x 25 x 1 Cartesian"
    double runSeconds = 1;     // rough run time on a 4-core desktop
};

const char* toString(ExampleLevel l);
void placeholderFluids(Case& c, bool metric);    // SPE1 PVT (internal: placeholders for unsupported decks)
std::vector<ExampleInfo> listExamples();
// Builds the example. Throws std::runtime_error for an unknown id. The history-matching exercise runs its hidden
// "truth" model once to create the observed data (a few seconds).
Case makeExample(const std::string& id);
// Heavy-oil thermal fluid (HEAVY dead + LIGHT distillable, OILVISCT, rock heat, heat loss) into c.thermal (enabled),
// converted to the case unit system; tInitC: initial temperature in deg C
void thermalPreset(Case& c, double tInitC = 25);

// C13: parametric EGS designer (FORGE-style multi-stage doublet): two parallel horizontal wells (injector / producer,
// along x, 'wellSep' apart in y) connected by 'stages' transverse vertical EDFM fractures (Barton-Bandis closure +
// Coulomb hydro-shear) in a tight hot crystalline rock; water circulation with a cold injection temperature. METRIC.
struct EgsDesign {
    int stages = 5;                 // number of transverse fractures
    double spacing = 80;            // m between stages
    double wellSep = 150;           // m between the injector and the producer
    double fracExt = 40;            // m the fractures extend beyond each well
    double fracHeight = 150;        // m
    double aperture = 0.5e-3;       // m (initial hydraulic aperture)
    bool mech = true;               // stress-dependent aperture / hydro-shear
    double depth = 2500;            // m, depth of the wells (centre of the fractured zone)
    double tRes = 210;              // deg C, initial rock temperature
    double kMatrix = 0.01;          // mD
    double poro = 0.01;             // -
    double rate = 2500;             // sm3/d water injection (~29 kg/s)
    double tInj = 50;               // deg C
    double drawdown = 20;           // bar below hydrostatic at the producer (BHP control)
    double years = 10;
    double cell = 25;               // m, horizontal cell size away from the fractures
    int nz = 3;                     // layers: 50 m cap + (nz - 2) layers over the fractured zone + 50 m base
    bool refine = true;             // geometric near-fracture refinement normal to the planes (early-time heat exchange:
                                    // see the Gringarten verification - uniform coarse cells over-predict early cooling)
    double refineFirst = 2.0, refineRatio = 1.8;   // m width of the cell holding the plane, growth ratio
    double sHmaxGrad = 26e3, sHminGrad = 17e3;   // Pa/m (S_hmin along the wells: fractures normal to x)
    // C14: stimulation sizing - when on, every stage is a hydraulic fracture grown from the injector by a PKN / KGD
    // propagation simulation (hydrofrac.h; height = fracHeight, E / nu of the rock, sigma_min = S_hmin at the wells):
    // the planes span the injector +- the simulated half-length (fracExt ignored; a fracture shorter than the well
    // separation does not reach the producer) and the initial hydraulic aperture is stimRetain x the mean stimulated
    // width (residual opening after shut-in from shear dilation / self-propping). Stages are sized independently
    // (no stress shadowing).
    bool stimulate = false;
    int stimModel = 0;              // 0 PKN (height-contained, L >> H), 1 KGD (plane strain, H >> L)
    double stimRate = 0.1;          // m3/s per stage (0.1 m3/s ~ 38 bbl/min)
    double stimTime = 3600;         // s pumped per stage
    double stimMu = 1e-3;           // Pa s (slickwater ~ 1 cP)
    double stimCl = 0;              // Carter leak-off coefficient, m / sqrt(s)
    double stimKic = 1.5e6;         // Pa sqrt(m), granite ~ 1-2 MPa sqrt(m)
    double stimRetain = 0.1;        // residual aperture / mean stimulated width
    double rockE = 40e9, rockNu = 0.25;
};
Case makeEgsCase(const EgsDesign& d);
// C14: the stimulation simulation of one stage for the design (false + r.error on failure)
bool egsStimulation(const EgsDesign& d, hf::Result& r);

} // namespace ressim
