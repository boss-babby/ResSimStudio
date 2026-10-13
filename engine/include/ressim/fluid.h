// Black-oil PVT, rock compressibility and saturation functions.
#pragma once
#include "ressim/ad.h"
#include "ressim/case.h"
#include <vector>

namespace ressim {

// Piecewise-linear table y(x). Extrapolates linearly (PVT) or clamps (saturation tables).
struct Table1D {
    std::vector<double> x, y;
    bool clampEnds = false;
    double eval(double xv, double* dydx = nullptr) const;
    double inverse(double yv) const;   // for monotone tables
    bool empty() const { return x.empty(); }
    template <int N> AD<N> operator()(const AD<N>& xv) const { double d; double f = eval(xv.v, &d); return chain(xv, f, d); }
};

class FluidModel {
public:
    void setup(const Case& c, int pvtRegion = 0);   // 0-based PVTNUM region

    // ---- water (PVTW, ECLIPSE formulas)
    AD3 bw(const AD3& p) const;            // 1/Bw
    AD3 bwMu(const AD3& p) const;          // 1/(Bw*mu_w)
    // ---- oil
    AD3 bo(const AD3& p, const AD3& rs) const;     // 1/Bo
    AD3 boMu(const AD3& p, const AD3& rs) const;   // 1/(Bo*mu_o)
    double rsSat(double p, double* drsdp = nullptr) const;
    AD3 rsSat(const AD3& p) const { double d; double r = rsSat(p.v, &d); return chain(p, r, d); }
    double pbOfRs(double rs) const;
    // ---- gas (dry gas; with vaporised oil these return the saturated values)
    AD3 bg(const AD3& p) const;
    AD3 bgMu(const AD3& p) const;
    // ---- wet gas (PVTG): 1/Bg and 1/(Bg mu_g) as functions of pressure and Rv
    AD3 bg(const AD3& p, const AD3& rv) const;
    AD3 bgMu(const AD3& p, const AD3& rv) const;
    double rvSat(double p, double* drvdp = nullptr) const;
    AD3 rvSat(const AD3& p) const { double d; double r = rvSat(p.v, &d); return chain(p, r, d); }
    double wetGas(double p, double rv, int which, double* dp, double* drv) const;   // which: 0 -> 1/Bg, 1 -> mu_g
    // ---- rock
    AD3 pvMult(const AD3& p) const;
    double pvMult(double p) const;

    // scalar convenience versions
    double bwS(double p) const { return bw(AD3(p)).v; }
    double boS(double p, double rs) const { return bo(AD3(p), AD3(rs)).v; }
    double bgS(double p) const { return bg(AD3(p)).v; }

    bool live() const { return live_; }
    bool vap() const { return vap_; }
    const DensitySC& rhoSC() const { return rho_; }
    double gasFactor() const { return gasFactor_; }

private:
    double liveOil(double p, double rs, int which) const; // which: 0 -> b, 1 -> b/mu
    AD3 liveOilAD(const AD3& p, const AD3& rs, int which) const;

    PVTW pvtw_; RockC rock_; DensitySC rho_;
    double gasFactor_ = 1.0;
    bool live_ = false, vap_ = false;
    // wet gas: per pressure node, 1/Bg and mu_g vs Rv (ascending Rv), saturated Rv of each node, saturated Rv vs p
    std::vector<double> pg_, rvMax_;
    std::vector<Table1D> gB_, gMu_;
    Table1D rvOfP_;
    Table1D deadB_, deadBMu_;            // dead oil: 1/B and 1/(B mu) vs p
    Table1D gasB_, gasBMu_;              // dry gas
    // live oil
    std::vector<double> rs_, pb_, bSat_, bmuSat_;
    std::vector<Table1D> uB_, uBMu_;     // per saturated node: ratio vs (p - pb)
    Table1D rsOfP_, pbOfRs_;
};

class SatFunctions {
public:
    // saturation end points of a cell (end-point scaling) or of the tables
    struct EndPoints { double swl = 0, swcr = 0, swu = 1, sgl = 0, sgcr = 0, sgu = 1, sowcr = 0, sogcr = 0; };
    void setup(const Case& c, int satRegion = 0);   // 0-based SATNUM region
    const EndPoints& tableEndPoints() const { return ep_; }
    // scaled versions (e == nullptr: unscaled tables), ECLIPSE two-point scaling
    AD3 krw(const AD3& sw, const EndPoints* e) const;
    AD3 krg(const AD3& sg, const EndPoints* e) const;
    AD3 kro(const AD3& sw, const AD3& sg, const EndPoints* e) const;
    AD3 pcow(const AD3& sw, const EndPoints* e) const;
    AD3 pcog(const AD3& sg, const EndPoints* e) const;
    double swInverseFromPc(double pc, const EndPoints* e) const;
    double sgInverseFromPc(double pc, const EndPoints* e) const;
    double swInverseFromPcGW(double pcgw, const EndPoints* e) const;
    // relative permeabilities and capillary pressures
    AD3 krw(const AD3& sw) const { return krw_(sw); }
    AD3 krg(const AD3& sg) const { return krg_(sg); }
    AD3 kro(const AD3& sw, const AD3& sg) const;     // ECLIPSE default 3-phase model
    AD3 pcow(const AD3& sw) const { return pcow_(sw); }
    AD3 pcog(const AD3& sg) const { return pcog_(sg); }
    double swco() const { return swco_; }
    double swInverseFromPc(double pc) const;   // equilibrium
    double sgInverseFromPc(double pc) const;
    double swInverseFromPcGW(double pcgw) const;   // gas-water contact: Pcow(Sw) + Pcog(1 - Sw) = Pcgw
    double sgMax() const { return sgmax_; }
    double swMax() const { return swmax_; }
private:
    Table1D krw_, krow_, pcow_, krg_, krog_, pcog_;
    double swco_ = 0, swmax_ = 1, sgmax_ = 1;
    int kroModel_ = 0;              // 0 default, 1 STONE1, 2 STONE2
    double stoneEta_ = 1.0;
    AD3 stone(const AD3& sw, const AD3& sg, const AD3& krow, const AD3& krog, const AD3& krw, const AD3& krg,
              double krocw, double swco) const;
    bool water_ = true, gas_ = true;
    EndPoints ep_;
};

} // namespace ressim
