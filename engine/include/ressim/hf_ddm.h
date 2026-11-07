// C26 (C25-4 "FE apertures"): mechanical apertures of the in-run hydraulic fractures (HFRAC + HFGROW ... DDMWIDTH) from a
// 3-D displacement-discontinuity solve on the fracture pieces (ddm.h, half space with free surface at depth 0).
//
// After every report step (black oil): the open pieces of all DDMWIDTH fractures are the elements (vertical rectangles of
// the EDFM pieces, metres), loaded by the net pressure p_piece - sigma_closure (the closure stress of the HFGROW criterion:
// sigmaMin + sigmaGrad (z - zMid) or the FE minimum stress with FESTRESS). The opening w >= 0 solves K_o w = -p_net with
// contact (closed elements carry a compressive net traction), including the interaction of all pieces of all fractures
// (stress shadow between parallel fractures). Effective aperture w_eff = max(w, 0.1 w_input) (propped / residual aperture);
// conductance multiplier per piece m = (w_eff / w_input)^3 (cubic law, HFRAC perm <= 0) or w_eff / w_input (propped pack of
// fixed permeability), applied through Simulator::setHfConductance (fracture-fracture connections: harmonic mean).
// Limits: elastic homogeneous rock (no layering of E), explicit one-report-step lag, the fracture pore volume keeps the input
// width, the growth criterion is still the HFGROW K_I tip test (the aperture does not move the tip).
#pragma once
#include <string>
#include <vector>

namespace ressim {
struct Case;
class Simulator;
namespace hfddm {

struct Output {
    std::vector<double> width;     // per HF cell (piece): DDM aperture (m), 0 = closed / not a DDMWIDTH fracture
    std::vector<double> mult;      // per HF cell: conductance multiplier applied
    std::vector<double> pnet;      // per HF cell: net pressure (case units)
    double maxWidth = 0, volume = 0;   // m, m^3 (DDM opening volume of all pieces)
    int open = 0, elements = 0, updates = 0;
    double seconds = 0;
    struct Step { double t = 0, maxWidth = 0, volume = 0, minMult = 1, maxMult = 1, maxPnet = 0; int open = 0; std::vector<float> width; };
    std::vector<Step> steps;       // one per update (report step): history for the results plots / along-fracture profile
};

class Aperture {
public:
    // false when no fracture asks for DDMWIDTH (err explains); uses sim.caseData() (expanded pieces) and sim.grid()
    bool setup(const Simulator& sim, std::string& err);
    bool active() const { return active_; }
    void update(Simulator& sim);
    const Output& out() const { return out_; }
private:
    bool active_ = false, halfSpace_ = true;
    std::vector<int> pieces_;      // HF cell indices handled
    std::vector<double> Ko_;       // dense opening matrix over pieces_
    std::vector<double> wIn_, zMid_;
    std::vector<char> cubic_;
    double Mlen_ = 1, Pa_ = 1e5;
    Output out_;
};

}  // namespace hfddm
}  // namespace ressim
