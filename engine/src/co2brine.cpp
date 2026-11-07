// CO2 - brine thermodynamics (CO2STORE). See co2brine.h for the models and references.
#include "ressim/co2brine.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace ressim {
namespace co2 {

namespace {
constexpr double R_BAR = 83.1447;            // bar cm3 / (mol K)
constexpr double B_CO2 = 27.80, B_H2O = 18.18, A_H2O_CO2 = 7.89e7;   // cm3/mol, bar cm6 K^0.5 / mol2
constexpr double V_H2O = 18.1, V_CO2_G = 32.6, V_CO2_L = 32.0;       // average partial molar volumes, cm3/mol
constexpr double TC_CO2 = 304.1282;          // K

// real roots of V^3 + c2 V^2 + c1 V + c0 = 0 (ascending)
int cubic(double c2, double c1, double c0, double r[3]) {
    const double q = (3 * c1 - c2 * c2) / 9.0, rr = (9 * c2 * c1 - 27 * c0 - 2 * c2 * c2 * c2) / 54.0;
    const double D = q * q * q + rr * rr;
    if (D > 0) {
        const double s = std::cbrt(rr + std::sqrt(D)), t = std::cbrt(rr - std::sqrt(D));
        r[0] = s + t - c2 / 3.0;
        return 1;
    }
    const double th = std::acos(std::max(-1.0, std::min(1.0, rr / std::sqrt(-q * q * q))));
    const double m = 2 * std::sqrt(-q);
    r[0] = m * std::cos(th / 3) - c2 / 3; r[1] = m * std::cos((th + 2 * M_PI) / 3) - c2 / 3; r[2] = m * std::cos((th + 4 * M_PI) / 3) - c2 / 3;
    std::sort(r, r + 3);
    return 3;
}

// Redlich-Kwong molar volume of CO2 (cm3/mol) at P (bar), T (K); liquid = true if the liquid root was chosen
double rkVolume(const Fluid& f, double Pbar, bool* liquid = nullptr) {
    const double T = f.T, a = f.aCO2, b = B_CO2, sT = std::sqrt(T);
    double r[3];
    const int n = cubic(-R_BAR * T / Pbar, -(R_BAR * T * b / Pbar - a / (Pbar * sT) + b * b), -a * b / (Pbar * sT), r);
    if (liquid) *liquid = false;
    if (n == 1) return r[0];
    const double vg = r[2], vl = r[0];
    if (vl <= b) return vg;
    // Spycher et al. (2003) root selection: w2 - w1 > 0 -> gas root
    const double w1 = Pbar * (vg - vl);
    const double w2 = R_BAR * T * std::log((vg - b) / (vl - b)) + a / (sT * b) * std::log((vg + b) * vl / ((vl + b) * vg));
    if (w2 - w1 > 0) return vg;
    if (liquid) *liquid = true;
    return vl;
}

// Duan & Sun (2003) parameter form; T K, P bar
double duanSun(const double c[11], double T, double P) {
    return c[0] + c[1] * T + c[2] / T + c[3] * T * T + c[4] / (630 - T) + c[5] * P + c[6] * P * std::log(T) + c[7] * P / T +
           c[8] * P / (630 - T) + c[9] * P * P / ((630 - T) * (630 - T)) + c[10] * T * std::log(std::max(P, 1e-6));
}
const double kLambda[11] = {-0.411370585, 6.07632013e-4, 97.5347708, 0, 0, 0, 0, -0.0237622469, 0.0170656236, 0, 1.41335834e-5};
const double kZeta[11] = {3.36389723e-4, -1.98298980e-5, 0, 0, 0, 0, 0, 2.12220830e-3, -5.24873303e-3, 0, 0};

// closed form (salt-free aqueous x, gas y_H2O) at P (Pa)
void solve(const Fluid& f, double P, double& xsf, double& yw) {
    const double Pb = std::max(P * 1e-5, 1e-3), T = f.T, sT15 = std::pow(T, 1.5);
    bool liq = false;
    const double V = rkVolume(f, Pb, &liq), b = B_CO2, a = f.aCO2;
    const double lv = std::log((V + b) / V), Z = Pb * V / (R_BAR * T);
    const double base = std::log(V / (V - b)) - std::log(Z);
    const double lnPhiC = base + b / (V - b) - 2 * a / (R_BAR * sT15 * b) * lv + a * b / (R_BAR * sT15 * b * b) * (lv - b / (V + b));
    const double lnPhiW = base + B_H2O / (V - b) - 2 * A_H2O_CO2 / (R_BAR * sT15 * b) * lv + a * B_H2O / (R_BAR * sT15 * b * b) * (lv - b / (V + b));
    const double lnK0c = liq ? std::log(10.0) * (1.169 + 1.368e-2 * f.tC - 5.380e-5 * f.tC * f.tC) : f.lnK0c;
    const double vco2 = liq ? V_CO2_L : V_CO2_G;
    const double A = std::exp(f.lnK0w - lnPhiW + (Pb - 1) * V_H2O / (R_BAR * T)) / Pb;
    const double B = std::exp(lnPhiC - lnK0c - (Pb - 1) * vco2 / (R_BAR * T)) * Pb / 55.508;
    double Bp = B, xs = 0;
    const double m = f.molality;
    if (m > 0) {
        const double lam = duanSun(kLambda, T, Pb), zet = duanSun(kZeta, T, Pb);
        const double gm = std::exp(2 * lam * m + zet * m * m);
        Bp = B / (gm * (1 + 2 * m / 55.508));
    }
    double y = 0, x = 0, mc = 0;
    for (int it = 0; it < (m > 0 ? 4 : 1); ++it) {
        xs = m > 0 ? 2 * m / (55.508 + 2 * m + mc) : 0.0;
        y = (1 - Bp - xs) / (1 / A - Bp);
        y = std::max(0.0, std::min(1.0, y));
        x = Bp * (1 - y);
        mc = 55.508 * x / std::max(1e-12, 1 - x - xs);
    }
    xsf = x / std::max(1e-12, 1 - xs);
    yw = y;
}
}  // namespace

Fluid makeFluid(double T, double molality) {
    Fluid f; f.T = T; f.molality = std::max(0.0, molality);
    const double t = T - 273.15; f.tC = t;
    f.aCO2 = 7.54e7 - 4.13e4 * T;
    f.lnK0w = std::log(10.0) * (-2.209 + 3.097e-2 * t - 1.098e-4 * t * t + 2.048e-7 * t * t * t);
    f.lnK0c = std::log(10.0) * (1.189 + 1.304e-2 * t - 5.446e-5 * t * t);
    f.salt = f.molality * MW_NACL / (1 + f.molality * MW_NACL);
    f.vphi = (37.51 - 9.585e-2 * t + 8.740e-4 * t * t - 5.044e-7 * t * t * t) * 1e-6;
    f.mwBrinePerWater = MW_H2O * (1 + f.molality * MW_NACL);
    const double S = f.salt;
    f.muBrine = 1e-3 * (0.1 + 0.333 * S + (1.65 + 91.9 * S * S * S) * std::exp(-(0.42 * std::pow(std::pow(S, 0.8) - 0.17, 2) + 0.045) * std::pow(std::max(t, 0.0), 0.8)));
    return f;
}

double molalityFromPpm(double ppm) { const double w = ppm * 1e-6; return w <= 0 ? 0.0 : w / ((1 - w) * MW_NACL); }
double ppmFromMolality(double m) { const double s = m * MW_NACL; return 1e6 * s / (1 + s); }

void solubility(const Fluid& f, double P, double& x, double& y) { solve(f, P, x, y); }
void solubility(const Fluid& f, double P, double& x, double& y, double& dx, double& dy) {
    solve(f, P, x, y);
    const double h = std::max(1.0, 1e-5 * P);
    double xp, yp, xm, ym;
    solve(f, P + h, xp, yp); solve(f, P - h, xm, ym);
    dx = (xp - xm) / (2 * h); dy = (yp - ym) / (2 * h);
}
double co2Molality(const Fluid& f, double P) {
    double x, y; solve(f, P, x, y);
    return 55.508 * x / (1 - x);
}

double gasMolarVolume(const Fluid& f, double P, double* dVdP) {
    const double Pb = std::max(P * 1e-5, 1e-3);
    const double V = rkVolume(f, Pb);
    if (dVdP) {   // implicit derivative of G = P V^3 - RT V^2 - V (RT b - a / sqrt T + b^2 P) - a b / sqrt T
        const double T = f.T, a = f.aCO2, b = B_CO2, sT = std::sqrt(T);
        const double dGdV = 3 * Pb * V * V - 2 * R_BAR * T * V - (R_BAR * T * b - a / sT + b * b * Pb);
        const double dGdP = V * V * V - b * b * V;
        *dVdP = -dGdP / dGdV * 1e-6 * 1e-5;   // cm3/(mol bar) -> m3/(mol Pa)
    }
    return V * 1e-6;
}
double co2Density(const Fluid& f, double P) { return MW_CO2 / gasMolarVolume(f, P); }

double co2Viscosity(double T, double rho, double* dmu) {
    static const double a[5] = {0.235156, -0.491266, 5.211155e-2, 5.347906e-2, -1.537102e-2};
    const double Ts = T / 251.196, lt = std::log(Ts);
    double lo = 0, p = 1; for (double ai : a) { lo += ai * p; p *= lt; }
    const double eta0 = 1.00697 * std::sqrt(T) / std::exp(lo);
    const double d11 = 0.4071119e-2, d21 = 0.7198037e-4, d64 = 0.2411697e-16, d81 = 0.2971072e-22, d82 = -0.1627888e-22;
    const double r2 = rho * rho, r6 = r2 * r2 * r2, r8 = r6 * r2;
    const double de = d11 * rho + d21 * r2 + d64 * r6 / (Ts * Ts * Ts) + d81 * r8 + d82 * r8 / Ts;
    if (dmu) *dmu = 1e-6 * (d11 + 2 * d21 * rho + 6 * d64 * r6 / rho / (Ts * Ts * Ts) + 8 * d81 * r8 / rho + 8 * d82 * r8 / rho / Ts);
    return 1e-6 * (eta0 + de);
}

double brineDensity(const Fluid& f, double P, double* drho) {
    const double t = f.tC, Pm = P * 1e-6, S = f.salt;
    const double rw = 1 + 1e-6 * (-80 * t - 3.3 * t * t + 0.00175 * t * t * t + 489 * Pm - 2 * t * Pm + 0.016 * t * t * Pm - 1.3e-5 * t * t * t * Pm -
                                  0.333 * Pm * Pm - 0.002 * t * Pm * Pm);
    const double rb = rw + S * (0.668 + 0.44 * S + 1e-6 * (300 * Pm - 2400 * Pm * S + t * (80 + 3 * t - 3300 * S - 13 * Pm + 47 * Pm * S)));
    if (drho) {
        const double dw = 1e-6 * (489 - 2 * t + 0.016 * t * t - 1.3e-5 * t * t * t - 0.666 * Pm - 0.004 * t * Pm);
        const double db = dw + S * 1e-6 * (300 - 2400 * S + t * (-13 + 47 * S));
        *drho = db * 1000.0 * 1e-6;   // g/cm3 per MPa -> kg/m3 per Pa
    }
    return rb * 1000.0;
}
double brineViscosity(const Fluid& f) { return f.muBrine; }

double aqueousMolarVolume(const Fluid& f, double P, double x, double* dP, double* dx) {
    double drb = 0;
    const double rb = brineDensity(f, P, &drb);
    const double vw = f.mwBrinePerWater / rb;   // per mol H2O
    if (dP) *dP = (1 - x) * (-f.mwBrinePerWater / (rb * rb)) * drb;
    if (dx) *dx = f.vphi - vw;
    return (1 - x) * vw + x * f.vphi;
}
double aqueousDensity(const Fluid& f, double P, double x) {
    return ((1 - x) * f.mwBrinePerWater + x * MW_CO2) / aqueousMolarVolume(f, P, x);
}

FlashResult flash(const Fluid& f, double P, double z) {
    FlashResult r;
    double x, yw, dx, dy;
    solubility(f, P, x, yw, dx, dy);
    r.x = x; r.y = 1 - yw; r.dx_dP = dx; r.dy_dP = -dy;
    if (z <= r.x) { r.phases = 1; r.gas = false; r.beta = 0; }
    else if (z >= r.y) { r.phases = 1; r.gas = true; r.beta = 1; }
    else { r.phases = 2; r.beta = (z - r.x) / (r.y - r.x); }
    return r;
}

double gasSurfaceVolume(double Tstd, double Pstd) { return gasMolarVolume(makeFluid(Tstd, 0), Pstd); }
double brineSurfaceVolume(const Fluid& f, double Tstd, double Pstd) {
    Fluid s = makeFluid(Tstd, f.molality);
    return s.mwBrinePerWater / brineDensity(s, Pstd);
}

std::string describe(const Fluid& f) {
    char b[256];
    std::snprintf(b, sizeof(b), "CO2-brine (Spycher-Pruess / Duan-Sun, Batzle-Wang, Garcia, Fenghour): T %.2f C, NaCl %.3f mol/kg (%.0f ppm)%s", f.tC,
                  f.molality, ppmFromMolality(f.molality), (f.tC < 12 || f.tC > 100) ? " - outside the 12-100 C calibration range" : "");
    return b;
}


std::vector<MineralReaction> mineralPresets() {
    std::vector<MineralReaction> v;
    MineralReaction a;   // CaAl2Si2O8 + CO2 + 2 H2O -> CaCO3 + Al2Si2O5(OH)4
    a.name = "Anorthite -> calcite + kaolinite"; a.logk25 = -9.12; a.Ea = 17.8e3; a.nuCO2 = 1; a.nuH2O = 2;
    a.dVsolid = (36.93 + 99.52 - 100.79) * 1e-6; a.vmReactant = 100.79e-6; v.push_back(a);
    MineralReaction b;   // NaAlSi3O8 + CO2 + H2O -> NaAlCO3(OH)2 + 3 SiO2
    b.name = "Albite -> dawsonite + quartz"; b.logk25 = -12.56; b.Ea = 69.8e3; b.nuCO2 = 1; b.nuH2O = 1;
    b.dVsolid = (59.30 + 3 * 22.69 - 100.07) * 1e-6; b.vmReactant = 100.07e-6; b.reactant0 = 0.10; v.push_back(b);
    MineralReaction c;   // fast reference: calcite precipitation from dissolved CO2 with a calcium source (fast kinetics)
    c.name = "Calcite (Ca-rich brine, fast)"; c.logk25 = -5.81; c.Ea = 23.5e3; c.nuCO2 = 1; c.nuH2O = 0;
    c.dVsolid = 36.93e-6; c.vmReactant = 36.93e-6; c.reactant0 = 0.01; c.area = 100; v.push_back(c);
    return v;
}

double permFactor(const TrappingOptions& t, double phiRatio) {
    phiRatio = std::clamp(phiRatio, 0.0, 1.0);
    if (!t.permFact.empty()) {
        const auto& T = t.permFact;
        if (phiRatio <= T.front().first) return T.front().second;
        if (phiRatio >= T.back().first) return T.back().second;
        for (size_t i = 1; i < T.size(); ++i)
            if (phiRatio <= T[i].first) {
                const double w = (phiRatio - T[i - 1].first) / std::max(1e-12, T[i].first - T[i - 1].first);
                return T[i - 1].second + w * (T[i].second - T[i - 1].second);
            }
        return T.back().second;
    }
    // Verma & Pruess (1988) tubes-in-series: theta = (phi/phi0 - phir) / (1 - phir), omega = 1 + (1/G) / (1/phir - 1)
    const double phir = std::clamp(t.vpPhiR, 0.0, 0.999), G = std::max(1e-6, t.vpGamma);
    const double th = (phiRatio - phir) / (1 - phir);
    if (th <= 0) return 0.0;
    const double om = 1 + (1 / G) / (1 / phir - 1);
    const double r = th / (th + om - 1);
    return th * th * (1 - G + G / (om * om)) / (1 - G + G * r * r);
}
}  // namespace co2
}  // namespace ressim
