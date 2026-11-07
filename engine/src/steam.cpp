#include "ressim/steam.h"
#include <algorithm>
#include <cmath>

namespace ressim {
namespace steam {

namespace {
static const int kR1I[] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32};
static const int kR1J[] = {-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41};
static const double kR1n[] = {0.14632971213167001, -0.84548187169113997, -3.756360367204, 3.3855169168385002, -0.95791963387872003, 0.15772038513228001, -0.016616417199501, 0.00081214629983567997, 0.00028319080123804, -0.00060706301565873996, -0.018990068218419, -0.032529748770504997, -0.021841717175413999, -5.2838357969930002e-05, -0.00047184321073266998, -0.00030001780793025999, 4.7661393906987001e-05, -4.4141845330845997e-06, -7.2694996297594001e-16, -3.1679644845054002e-05, -2.8270797985312e-06, -8.5205128120103004e-10, -2.2425281907999999e-06, -6.5171222895601002e-07, -1.4341729937923999e-13, -4.0516996860117e-07, -1.2734301741640999e-09, -1.7424871230634001e-10, -6.8762131295530996e-19, 1.4478307828521001e-20, 2.6335781662795e-23, -1.1947622640071e-23, 1.8228094581404e-24, -9.3537087292457998e-26};
static const int kR2I[] = {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24};
static const int kR2J[] = {0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58};
static const double kR2n[] = {-0.0017731742473212999, -0.017834862292357999, -0.045996013696365003, -0.057581259083432, -0.050325278727930002, -3.3032641670203e-05, -0.00018948987516315, -0.0039392777243355001, -0.043797295650572998, -2.6674547914087001e-05, 2.0481737692308999e-08, 4.3870667284435001e-07, -3.2277677238570002e-05, -0.0015033924542148, -0.040668253562648998, -7.8847309559367001e-10, 1.2790717852285001e-08, 4.8225372718507002e-07, 2.2922076337661001e-06, -1.6714766451061001e-11, -0.0021171472321354998, -23.895741934103999, -5.9059564324270004e-18, -1.2621808899101e-06, -0.038946842435739003, 1.1256211360459e-11, -8.2311340897998004, 1.9809712802088e-08, 1.0406965210174e-19, -1.0234747095929e-13, -1.0018179379511e-09, -8.0882908646984998e-11, 0.10693031879409, -0.33662250574170999, 8.9185845355420999e-25, 3.0629316876231997e-13, -4.2002467698208001e-06, -5.9056029685639003e-26, 3.7826947613457002e-06, -1.2768608934681e-15, 7.3087610595061e-29, 5.5414715350778001e-17, -9.4369707241209998e-07};
static const int kR2J0[] = {0, 1, -5, -4, -3, -2, -1, 2, 3};
static const double kR2n0[] = {-9.6927686500216996, 10.086655968018, -0.0056087911283020002, 0.071452738081454997, -0.40710498223927999, 1.4240819171443999, -4.3839511319450004, -0.28408632460771999, 0.021268463753307001};

const double kN4[] = {0, 0.11670521452767E+04, -0.72421316703206E+06, -0.17073846940092E+02, 0.12020824702470E+05,
                      -0.32325550322333E+07, 0.14915108613530E+02, -0.48232657361591E+04, 0.40511340542057E+06,
                      -0.23855557567849E+00, 0.65017534844798E+03};
inline double ipow(double x, int n) {
    if (n == 0) return 1.0;
    double r = 1.0, b = n > 0 ? x : 1.0 / x; int m = n > 0 ? n : -n;
    while (m) { if (m & 1) r *= b; b *= b; m >>= 1; }
    return r;
}
} // namespace

double psat(double T) {
    T = std::clamp(T, 273.15, kTc);
    const double* n = kN4;
    const double th = T + n[9] / (T - n[10]);
    const double A = th * th + n[1] * th + n[2], B = n[3] * th * th + n[4] * th + n[5], C = n[6] * th * th + n[7] * th + n[8];
    const double x = 2 * C / (-B + std::sqrt(B * B - 4 * A * C));
    return x * x * x * x * 1e6;
}
double tsat(double p) {
    p = std::clamp(p, 611.212677, kPc) * 1e-6;
    const double* n = kN4;
    const double b = std::pow(p, 0.25);
    const double E = b * b + n[3] * b + n[6], F = n[1] * b * b + n[4] * b + n[7], G = n[2] * b * b + n[5] * b + n[8];
    const double D = 2 * G / (-F - std::sqrt(F * F - 4 * E * G));
    return (n[10] + D - std::sqrt((n[10] + D) * (n[10] + D) - 4 * (n[9] + n[10] * D))) / 2;
}
double dpsatdT(double T) {
    const double h = 1e-3 * std::max(1.0, T * 1e-3);
    return (psat(T + h) - psat(T - h)) / (2 * h);
}

static Props region1(double T, double p) {
    const double pi = std::max(p, 1.0) / 16.53e6, tau = 1386.0 / T;
    double gp = 0, gt = 0, gtt = 0;
    const double a = 7.1 - pi, b = tau - 1.222;
    for (int k = 0; k < 34; ++k) {
        const int I = kR1I[k], J = kR1J[k];
        const double n = kR1n[k];
        const double aI = ipow(a, I), bJ = ipow(b, J);
        if (I) gp -= n * I * ipow(a, I - 1) * bJ;
        if (J) { const double bJ1 = ipow(b, J - 1); gt += n * aI * J * bJ1; gtt += n * aI * J * (J - 1) * ipow(b, J - 2); }
    }
    Props r;
    r.v = pi * gp * kR * T / std::max(p, 1.0);
    r.h = tau * gt * kR * T;
    r.u = r.h - std::max(p, 1.0) * r.v;
    r.cp = -kR * tau * tau * gtt;
    return r;
}

// ------------------------------------------------------------------------------------------------ IF97 region 3
// Helmholtz form f(rho, T) / RT = n1 ln(delta) + sum n_i delta^I_i tau^J_i, delta = rho / 322, tau = 647.096 / T (IF97
// eq. 28, table 30). Covers 623.15 K <= T <= T_B23(p): the compressed liquid between 623.15 K and Tc, the near-critical and
// the dense supercritical fluid, where the region-1 / region-2 Gibbs equations are only (badly) extrapolated.
namespace {
const int kR3I[] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11};
const int kR3J[] = {0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, 2, 4, 16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26, 2, 26, 2, 26, 0, 1, 26};
const double kR3n1 = 1.0658070028513;
const double kR3n[] = {-0.15732845290239e2, 0.20944396974307e2, -0.76867707878716e1, 0.26185947787954e1, -0.28080781148620e1,
                       0.12053369696517e1, -0.84566812812502e-2, -0.12654315477714e1, -0.11524407806681e1, 0.88521043984318,
                       -0.64207765181607, 0.38493460186671, -0.85214708824206, 0.48972281541877e1, -0.30502617256965e1,
                       0.39420536879154e-1, 0.12558408424308, -0.27999329698710, 0.13899799569460e1, -0.20189915023570e1,
                       -0.82147637173963e-2, -0.47596035734923, 0.43984074473500e-1, -0.44476435428739, 0.90572070719733,
                       0.70522450087967, 0.10770512626332, -0.32913623258954, -0.50871062041158, -0.22175400873096e-1,
                       0.94260751665092e-1, 0.16436278447961, -0.13503372241348e-1, -0.14834345352472e-1, 0.57922953628084e-3,
                       0.32308904703711e-2, 0.80964802996215e-4, -0.16557679795037e-3, -0.44923899061815e-4};
struct R3 { double p = 0, dpdr = 0, h = 0, cp = 0; };
R3 region3(double rho, double T) {
    const double d = rho / kRhoc, t = kTc / T;
    double fd = kR3n1 / d, fdd = -kR3n1 / (d * d), ft = 0, ftt = 0, fdt = 0;
    for (int k = 0; k < 39; ++k) {
        const int I = kR3I[k], J = kR3J[k]; const double n = kR3n[k];
        const double dI = ipow(d, I), tJ = ipow(t, J);
        if (I) { const double a = n * I * ipow(d, I - 1); fd += a * tJ; fdd += a * (I - 1) / d * tJ; if (J) fdt += a * J * ipow(t, J - 1); }
        if (J) { const double b = n * J * ipow(t, J - 1); ft += b * dI; ftt += b * (J - 1) / t * dI; }
    }
    R3 r;
    r.p = rho * kR * T * d * fd;
    r.dpdr = kR * T * (2 * d * fd + d * d * fdd);
    r.h = kR * T * (t * ft + d * fd);
    const double num = d * fd - d * t * fdt, den = 2 * d * fd + d * d * fdd;
    r.cp = kR * (-t * t * ftt + (den > 1e-12 ? num * num / den : 0.0));
    return r;
}
// saturated densities, Wagner & Pruss (1993) auxiliary equations (IAPWS-95 consistent; used to bracket the region-3
// roots below Tc, where a geometric walk can step over the narrow van der Waals loop near the critical point)
double rhoSatVap(double T) {
    const double t = std::max(1.0 - T / kTc, 0.0);
    const double c[] = {-2.03150240, -2.68302940, -5.38626492, -17.2991605, -44.7586581, -63.9201063};
    const double e[] = {2.0 / 6, 4.0 / 6, 8.0 / 6, 18.0 / 6, 37.0 / 6, 71.0 / 6};
    double s = 0; for (int k = 0; k < 6; ++k) s += c[k] * std::pow(t, e[k]);
    return kRhoc * std::exp(s);
}
double rhoSatLiq(double T) {
    const double t = std::max(1.0 - T / kTc, 0.0);
    const double b[] = {1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -6.74694450e5};
    const double e[] = {1.0 / 3, 2.0 / 3, 5.0 / 3, 16.0 / 3, 43.0 / 3, 110.0 / 3};
    double s = 1; for (int k = 0; k < 6; ++k) s += b[k] * std::pow(t, e[k]);
    return kRhoc * s;
}
// density root of p(rho, T) = p on the liquid-like (dense) or vapour-like branch; false when the branch has no root
// (beyond its spinodal below Tc). Below Tc the stable vapour branch is rho <= rho''(T) and the liquid branch rho >= rho'(T)
// (both monotone): bracket there, extend in small steps for metastable states. Above Tc the isotherm is monotone (one
// root). Safeguarded Newton inside the bracket, started from the branch side.
bool region3Rho(double T, double p, bool liquidLike, double& rho) {
    double lo = 0.7 * p / (kR * T);
    if (!(lo > 0)) return false;
    for (int k = 0; k < 20 && region3(lo, T).p > p; ++k) lo *= 0.5;
    double hi = 1.03 / region1(std::min(T, 623.15), p).v;
    for (int k = 0; k < 20 && region3(hi, T).p < p; ++k) hi *= 1.05;
    if (region3(lo, T).p > p || region3(hi, T).p < p) return false;
    if (T < kTc) {
        if (!liquidLike) {
            double b = std::max(lo * 1.0001, std::min(hi, rhoSatVap(T) * 1.0005));
            R3 r = region3(b, T);
            for (int k = 0; k < 400 && r.p < p; ++k) {   // metastable vapour: step on until p is crossed or the spinodal
                const double bn = b * 1.002; const R3 rn = region3(bn, T);
                if (!(rn.dpdr > 0)) return false;
                b = bn; r = rn;
            }
            if (r.p < p) return false;
            hi = b;
        } else {
            double b = std::min(hi / 1.0001, std::max(lo, rhoSatLiq(T) * 0.9995));
            R3 r = region3(b, T);
            for (int k = 0; k < 400 && r.p > p; ++k) {   // metastable (stretched) liquid
                const double bn = b / 1.002; const R3 rn = region3(bn, T);
                if (!(rn.dpdr > 0)) return false;
                b = bn; r = rn;
            }
            if (r.p > p) return false;
            lo = b;
        }
    }
    double x = liquidLike ? hi : lo;
    for (int it = 0; it < 200; ++it) {
        const R3 c = region3(x, T);
        const double F = c.p - p;
        if (F > 0) hi = x; else lo = x;
        double xn = c.dpdr > 0 ? x - F / c.dpdr : 0.5 * (lo + hi);
        if (!(xn > lo && xn < hi)) xn = 0.5 * (lo + hi);
        if (std::fabs(xn - x) <= 1e-12 * x || hi - lo <= 1e-12 * hi) { x = xn; break; }
        x = xn;
    }
    rho = x;
    return true;
}
double pB23(double T) { return (0.34805185628969e3 - 0.11671859879975e1 * T + 0.10192970039326e-2 * T * T) * 1e6; }
inline double smooth01(double s) { s = std::clamp(s, 0.0, 1.0); return s * s * (3 - 2 * s); }
} // namespace

// liquid: region 1 to 618.15 K; region 3 (liquid-like root) above 628.15 K, C1 blend in between (IF97 1/3 boundary at
// 623.15 K agrees to its consistency tolerance). Above Tc "liquid" is the single supercritical fluid (= vapour()).
Props liquid(double T, double p) {
    if (T <= 618.15) return region1(T, p);
    if (T > kTc) return vapour(T, p);
    const double pp = std::max(p, 1.0);
    double rho; R3 st;
    if (pp < pB23(std::max(T, 623.15)) * 0.98 || !region3Rho(T, pp, true, rho)) return region1(T, p);
    st = region3(rho, T);
    Props r; r.v = 1.0 / rho; r.h = st.h; r.u = st.h - pp * r.v; r.cp = st.cp;
    if (T >= 628.15) return r;
    const Props a = region1(T, p);
    const double w = smooth01((T - 618.15) / 10.0);
    Props m; m.v = (1 - w) * a.v + w * r.v; m.h = (1 - w) * a.h + w * r.h; m.u = (1 - w) * a.u + w * r.u; m.cp = (1 - w) * a.cp + w * r.cp;
    return m;
}

static void region2(double T, double p, double& zf, double& h, double& cp) {
    const double pi = std::max(p, 0.0) / 1e6, tau = 540.0 / T;
    double got = 0, gott = 0;
    for (int k = 0; k < 9; ++k) {
        const int J = kR2J0[k]; const double n = kR2n0[k];
        if (J) { got += n * J * ipow(tau, J - 1); gott += n * J * (J - 1) * ipow(tau, J - 2); }
    }
    double grp = 0, grt = 0, grtt = 0;
    const double b = tau - 0.5;
    for (int k = 0; k < 43; ++k) {
        const int I = kR2I[k], J = kR2J[k]; const double n = kR2n[k];
        const double pI = ipow(pi, I), bJ = ipow(b, J);
        grp += n * I * ipow(pi, I - 1) * bJ;
        if (J) { grt += n * pI * J * ipow(b, J - 1); grtt += n * pI * J * (J - 1) * ipow(b, J - 2); }
    }
    zf = 1 + pi * grp;
    h = tau * (got + grt) * kR * T;
    cp = -kR * tau * tau * (gott + grtt);
}
// IF97 region 5 (high-temperature steam, 1073.15 K - 2273.15 K, p <= 50 MPa; 2007 revision coefficients)
static void region5(double T, double p, double& zf, double& h, double& cp) {
    static const int J0[] = {0, 1, -3, -2, -1, 2};
    static const double n0[] = {-13.179983674201, 6.8540841634434, -0.024805148933466, 0.36901534980333, -3.1161318213925, -0.32961626538917};
    static const int I[] = {1, 1, 1, 2, 2, 3};
    static const int J[] = {1, 2, 3, 3, 9, 7};
    static const double n[] = {0.15736404855259e-2, 0.90153761673944e-3, -0.50270077677648e-2, 0.22440037409485e-5, -0.41163275453471e-5, 0.37919454822955e-7};
    const double pi = std::max(p, 0.0) / 1e6, tau = 1000.0 / T;
    double got = 0, gott = 0, grp = 0, grt = 0, grtt = 0;
    for (int k = 0; k < 6; ++k) {
        if (J0[k]) { got += n0[k] * J0[k] * std::pow(tau, J0[k] - 1); gott += n0[k] * J0[k] * (J0[k] - 1) * std::pow(tau, J0[k] - 2); }
        grp += n[k] * I[k] * ipow(pi, I[k] - 1) * ipow(tau, J[k]);
        grt += n[k] * ipow(pi, I[k]) * J[k] * ipow(tau, J[k] - 1);
        grtt += n[k] * ipow(pi, I[k]) * J[k] * (J[k] - 1) * ipow(tau, J[k] - 2);
    }
    zf = 1 + pi * grp;
    h = tau * (got + grt) * kR * T;
    cp = -kR * tau * tau * (gott + grtt);
}
// vapour: region 2 below 1063.15 K, region 5 above 1083.15 K, C1 smooth-step blend in between (the regions agree to
// IF97 consistency tolerances at 1073.15 K; the blend keeps E(T) smooth for the energy flash of combustion / UCG cells)
// region 3 above the B23 line (623.15 K <= T <= 863.15 K): vapour-like root, C1 blend with region 2 across B23 +- 0.3 MPa
static void vapourState(double T, double p, double& zf, double& h, double& cp) {
    if (T >= 623.15 && T < 863.15) {
        const double pb = pB23(T), band = 0.3e6;
        if (p > pb - band) {
            double rho; R3 st;
            if (region3Rho(T, p, false, rho) || (T > kTc && region3Rho(T, p, true, rho))) {
                st = region3(rho, T);
                const double z3 = p / (rho * kR * T);
                const double w = smooth01((p - (pb - band)) / (2 * band));
                if (w >= 1) { zf = z3; h = st.h; cp = st.cp; return; }
                double z2, h2, c2; region2(T, p, z2, h2, c2);
                zf = (1 - w) * z2 + w * z3; h = (1 - w) * h2 + w * st.h; cp = (1 - w) * c2 + w * st.cp;
                return;
            }
        }
    }
    if (T <= 1063.15) { region2(T, p, zf, h, cp); return; }
    if (T >= 1083.15) { region5(T, p, zf, h, cp); return; }
    double z2, h2, c2, z5, h5, c5;
    region2(T, p, z2, h2, c2); region5(T, p, z5, h5, c5);
    const double s = (T - 1063.15) / 20.0, w = s * s * (3 - 2 * s);
    zf = (1 - w) * z2 + w * z5; h = (1 - w) * h2 + w * h5; cp = (1 - w) * c2 + w * c5;
}
// Pressure on the isochore rho = 0.75 rho_c (IF97 region 3) above 623.15 K, +inf below: the thermal flash evaluates the
// enthalpy of water dissolved in a gas mixture (ideal mixing at its partial pressure) at min(p_w, this) so that it never
// approaches the critical point, where the saturated-vapour enthalpy h''(T) falls with an infinite slope and the mixture
// energy U(T) at fixed composition became non-monotone (two flash temperatures for one energy -> Newton oscillation).
double vapourMixCap(double T) {
    if (T < 623.15) return 1e300;
    return region3(0.75 * kRhoc, T).p;
}
double vapourZ(double T, double p) { double z, h, cp; vapourState(T, p, z, h, cp); return z; }
Props vapour(double T, double p) {
    double z, h, cp; vapourState(T, p, z, h, cp);
    Props r;
    r.v = z * kR * T / std::max(p, 1e-3);
    r.h = h; r.cp = cp;
    r.u = h - z * kR * T;
    return r;
}

double viscosity(double T, double rho) {
    const double Tr = T / kTc, Dr = rho / kRhoc;
    static const double H[] = {1.67752, 2.20462, 0.6366564, -0.241605};
    double s = 0; for (int i = 0; i < 4; ++i) s += H[i] / ipow(Tr, i);
    const double mu0 = 100 * std::sqrt(Tr) / s;
    static const int li[] = {0, 1, 2, 3, 0, 1, 2, 3, 5, 0, 1, 2, 3, 4, 0, 1, 0, 3, 4, 3, 5};
    static const int lj[] = {0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6};
    static const double Hij[] = {0.520094, 0.850895e-1, -0.108374e1, -0.289555, 0.222531, 0.999115, 0.188797e1, 0.126613e1,
                                 0.120573, -0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.161913, 0.257399,
                                 -0.325372e-1, 0.698452e-1, 0.872102e-2, -0.435673e-2, -0.593264e-3};
    double t = 0; for (int k = 0; k < 21; ++k) t += ipow(1 / Tr - 1, li[k]) * Hij[k] * ipow(Dr - 1, lj[k]);
    return mu0 * std::exp(Dr * t) * 1e-6;
}

} // namespace steam
} // namespace ressim
