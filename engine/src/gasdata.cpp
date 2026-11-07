// Gas / solute property library and element-formula utilities (see gasdata.h for the sources).
#include "ressim/gasdata.h"
#include <algorithm>
#include <cctype>
#include <cmath>

namespace ressim {
namespace gasdata {

namespace {
constexpr double kTc = 647.096, kPc = 22.064e6;
const std::vector<GasInfo>& table() {
    //                name    formula   mw         tc       pc        omega    hA          hB        hC         vInf      ksalt  cp300  cpSlope
    static const std::vector<GasInfo> t = {
        {"H2",   "H2",   0.0020159, 33.19,  1.313e6, -0.216, -4.73284,  6.08954,  6.06066,  25.2e-6, 0.23,  28.85, 0.00193},
        {"CO",   "CO",   0.0280101, 132.85, 3.494e6,  0.045, -10.52862, 5.13259,  12.01421, 35.0e-6, 0.30,  29.14, 0.00577},
        {"CO2",  "CO2",  0.0440095, 304.13, 7.3773e6, 0.2239, -8.55445, 4.01195,  9.52345,  -1.0,    0.237, 37.22, 0.02441},
        {"CH4",  "CH4",  0.0160425, 190.56, 4.599e6,  0.011, -10.44708, 4.66491,  12.12986, 37.3e-6, 0.288, 35.75, 0.05150},
        {"C2H6", "C2H6", 0.0300690, 305.32, 4.872e6,  0.099, -19.67563, 4.51222,  20.62567, 51.2e-6, 0.37,  52.85, 0.09980},
        {"C3H8", "C3H8", 0.0440956, 369.83, 4.248e6,  0.152, 0, 0, 0,                       67.0e-6, 0.40,  73.90, 0.14390},
        {"N2",   "N2",   0.0280134, 126.20, 3.398e6,  0.037, -9.67578,  4.72162,  11.70585, 35.7e-6, 0.30,  29.12, 0.00511},
        {"O2",   "O2",   0.0319988, 154.58, 5.043e6,  0.022, -9.44833,  4.43822,  11.42005, 31.0e-6, 0.31,  29.39, 0.00783},
        {"H2S",  "H2S",  0.0340809, 373.10, 8.963e6,  0.090, -4.51499,  5.23538,  4.42126,  35.0e-6, 0.14,  34.20, 0.01630},
        {"AR",   "Ar",   0.039948,  150.86, 4.898e6, -0.002, -8.40954,  4.29587,  10.52779, 32.0e-6, 0.32,  20.79, 0.0},
        {"HE",   "He",   0.0040026, 5.19,   0.227e6, -0.390, -3.52839,  7.12983,  4.47770,  31.0e-6, 0.25,  20.79, 0.0},
        {"NH3",  "NH3",  0.0170305, 405.40, 11.333e6, 0.256, 0, 0, 0,                       24.0e-6, 0.0,   35.60, 0.02990},
        {"SO2",  "SO2",  0.0640638, 430.80, 7.884e6,  0.245, 0, 0, 0,                       39.0e-6, 0.0,   39.90, 0.02060},
        {"H2O",  "H2O",  0.0180153, 647.096, 22.064e6, 0.3443, 0, 0, 0,                     18.0e-6, 0.0,   33.60, 0.01100},
    };
    return t;
}
std::string up(std::string s) { for (char& c : s) c = char(std::toupper(static_cast<unsigned char>(c))); return s; }
}  // namespace

const std::vector<GasInfo>& all() { return table(); }

const GasInfo* find(const std::string& name0) {
    std::string n = up(name0);
    n.erase(std::remove_if(n.begin(), n.end(), [](char c) { return c == ' ' || c == '_' || c == '-'; }), n.end());
    static const std::map<std::string, std::string> alias = {
        {"C1", "CH4"}, {"METHANE", "CH4"}, {"C2", "C2H6"}, {"ETHANE", "C2H6"}, {"C3", "C3H8"}, {"PROPANE", "C3H8"},
        {"HYDROGEN", "H2"}, {"CARBONMONOXIDE", "CO"}, {"CARBONDIOXIDE", "CO2"}, {"NITROGEN", "N2"}, {"OXYGEN", "O2"},
        {"ARGON", "AR"}, {"HELIUM", "HE"}, {"WATER", "H2O"}, {"STEAM", "H2O"}, {"AMMONIA", "NH3"}, {"HYDROGENSULFIDE", "H2S"}};
    auto a = alias.find(n);
    if (a != alias.end()) n = a->second;
    for (const GasInfo& g : table()) if (n == g.name) return &g;
    return nullptr;
}

double psatWagnerPruss(double T) {
    static const double a[6] = {-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502};
    static const double b[6] = {1.0, 1.5, 3.0, 3.5, 4.0, 7.5};
    T = std::clamp(T, 273.16, kTc);
    const double tr = T / kTc, tau = 1 - tr;
    double s = 0; for (int i = 0; i < 6; ++i) s += a[i] * std::pow(tau, b[i]);
    return kPc * std::exp(s / tr);
}

namespace {
double clampT(double T) { return std::clamp(T, 273.16, 0.999 * kTc); }
}

double henryIapws(double A, double B, double C, double T) {
    T = clampT(T);
    const double tr = T / kTc, tau = 1 - tr;
    const double f = A / tr + B * std::pow(tau, 0.355) / tr + C * std::pow(tr, -0.41) * std::exp(tau);
    return psatWagnerPruss(T) * std::exp(f);
}

double henryIapwsDlnDT(double A, double B, double C, double T) {
    if (T < 273.16 || T > 0.999 * kTc) return 0.0;   // clamped region: constant
    static const double a[6] = {-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502};
    static const double b[6] = {1.0, 1.5, 3.0, 3.5, 4.0, 7.5};
    const double tr = T / kTc, tau = 1 - tr;
    // d/dTr of ln p* = (1/Tr) sum a tau^b
    double s = 0, ds = 0;
    for (int i = 0; i < 6; ++i) { s += a[i] * std::pow(tau, b[i]); ds += a[i] * b[i] * std::pow(tau, b[i] - 1); }
    const double dlnp = -s / (tr * tr) - ds / tr;
    const double t355 = std::pow(tau, 0.355);
    const double df = -A / (tr * tr) + B * (-0.355 * std::pow(tau, -0.645) / tr - t355 / (tr * tr))
                      + C * std::exp(tau) * (-0.41 * std::pow(tr, -1.41) - std::pow(tr, -0.41));
    return (dlnp + df) / kTc;
}

double vInfCO2(double T) {
    const double t = T - 273.15;
    return (37.51 - 9.585e-2 * t + 8.740e-4 * t * t - 5.044e-7 * t * t * t) * 1e-6;
}

bool parseFormula(const std::string& f, std::map<std::string, double>& el) {
    el.clear();
    size_t i = 0;
    auto number = [&](double def) {
        size_t j = i;
        while (j < f.size() && (std::isdigit(static_cast<unsigned char>(f[j])) || f[j] == '.')) ++j;
        if (j == i) return def;
        const double v = std::atof(f.substr(i, j - i).c_str());
        i = j;
        return v;
    };
    std::vector<std::map<std::string, double>> stack(1);
    while (i < f.size()) {
        const char c = f[i];
        if (c == ' ') { ++i; continue; }
        if (c == '(' || c == '[') { stack.emplace_back(); ++i; continue; }
        if (c == ')' || c == ']') {
            ++i;
            if (stack.size() < 2) return false;
            const double m = number(1.0);
            std::map<std::string, double> top = stack.back(); stack.pop_back();
            for (auto& kv : top) stack.back()[kv.first] += kv.second * m;
            continue;
        }
        if (std::isupper(static_cast<unsigned char>(c))) {
            std::string e(1, c); ++i;
            while (i < f.size() && std::islower(static_cast<unsigned char>(f[i]))) e += f[i++];
            stack.back()[e] += number(1.0);
            continue;
        }
        return false;
    }
    if (stack.size() != 1) return false;
    el = stack[0];
    return !el.empty();
}

double carbonateCO2(const std::string& f) {
    double n = 0;
    for (size_t p = f.find("CO3"); p != std::string::npos; p = f.find("CO3", p + 3)) {
        // CO3 followed by a multiplier inside parentheses: (CO3)2
        double m = 1;
        if (p > 0 && f[p - 1] == '(' && p + 3 < f.size() && f[p + 3] == ')') {
            size_t j = p + 4; std::string d;
            while (j < f.size() && (std::isdigit(static_cast<unsigned char>(f[j])) || f[j] == '.')) d += f[j++];
            if (!d.empty()) m = std::atof(d.c_str());
        }
        n += m;
    }
    return n;
}

double elementImbalance(const std::vector<std::string>& formulas, const std::vector<double>& sto, const std::vector<double>& prod,
                        std::map<std::string, double>* perElement, std::vector<int>* unknown) {
    std::map<std::string, double> bal;
    for (size_t i = 0; i < formulas.size(); ++i) {
        const double nu = (i < prod.size() ? prod[i] : 0.0) - (i < sto.size() ? sto[i] : 0.0);
        if (nu == 0) continue;
        std::map<std::string, double> el;
        if (formulas[i].empty() || !parseFormula(formulas[i], el)) { if (unknown) unknown->push_back(int(i)); continue; }
        for (auto& kv : el) bal[kv.first] += nu * kv.second;
    }
    double worst = 0;
    for (auto& kv : bal) worst = std::max(worst, std::fabs(kv.second));
    if (perElement) *perElement = bal;
    return worst;
}

} // namespace gasdata
} // namespace ressim
