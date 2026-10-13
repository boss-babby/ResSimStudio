#include "ressim/fluid.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace ressim {

// ------------------------------------------------------------------ Table1D
double Table1D::eval(double xv, double* dydx) const {
    const size_t n = x.size();
    if (n == 0) { if (dydx) *dydx = 0; return 0; }
    if (n == 1) { if (dydx) *dydx = 0; return y[0]; }
    if (clampEnds) {
        if (xv <= x.front()) { if (dydx) *dydx = 0; return y.front(); }
        if (xv >= x.back())  { if (dydx) *dydx = 0; return y.back(); }
    }
    size_t i;
    if (xv <= x[1]) i = 0;
    else if (xv >= x[n - 2]) i = n - 2;
    else i = size_t(std::upper_bound(x.begin(), x.end(), xv) - x.begin()) - 1;
    const double dx = x[i + 1] - x[i];
    const double s = dx != 0 ? (y[i + 1] - y[i]) / dx : 0.0;
    if (dydx) *dydx = s;
    return y[i] + s * (xv - x[i]);
}

double Table1D::inverse(double yv) const {
    // assumes monotone y; returns x where y == yv (clamped to range)
    const size_t n = x.size();
    if (n == 0) return 0;
    const bool inc = y.back() >= y.front();
    for (size_t i = 0; i + 1 < n; ++i) {
        double a = y[i], b = y[i + 1];
        if ((inc && yv >= a && yv <= b) || (!inc && yv <= a && yv >= b)) {
            if (b == a) return x[i];
            return x[i] + (yv - a) / (b - a) * (x[i + 1] - x[i]);
        }
    }
    if (inc) return yv < y.front() ? x.front() : x.back();
    return yv > y.front() ? x.front() : x.back();
}

// ------------------------------------------------------------------ FluidModel
void FluidModel::setup(const Case& c, int region) {
    const PvtRegion PR = c.props.pvtRegion(region);
    pvtw_ = PR.pvtw; rock_ = PR.rock; rho_ = PR.density;
    gasFactor_ = Units::of(c.units).gasDensity;
    live_ = c.phases.oil && c.phases.gas && c.phases.disgas;

    if (c.phases.oil && !live_) {
        if (PR.pvdo.size() < 2) throw std::runtime_error("PVDO needs at least 2 rows");
        deadB_ = {}; deadBMu_ = {};
        for (auto& r : PR.pvdo) { deadB_.x.push_back(r.p); deadB_.y.push_back(1.0 / r.B);
                                       deadBMu_.x.push_back(r.p); deadBMu_.y.push_back(1.0 / (r.B * r.mu)); }
    }
    if (c.phases.gas && !(c.phases.oil && c.phases.vapoil)) {
        if (PR.pvdg.size() < 2) throw std::runtime_error("PVDG needs at least 2 rows");
        gasB_ = {}; gasBMu_ = {};
        for (auto& r : PR.pvdg) { gasB_.x.push_back(r.p); gasB_.y.push_back(1.0 / r.B);
                                       gasBMu_.x.push_back(r.p); gasBMu_.y.push_back(1.0 / (r.B * r.mu)); }
    }
    vap_ = c.phases.oil && c.phases.gas && c.phases.vapoil;
    if (vap_) {
        // PVTG (as OPM Flow): records sorted by pressure; a record with only the saturated row is extended with the
        // undersaturated steps of the next higher-pressure record that has them (same Rv steps, same relative change
        // of Bg and mu per step)
        auto t = PR.pvtg;
        if (t.size() < 2) throw std::runtime_error("PVTG needs at least 2 pressure records (VAPOIL)");
        std::sort(t.begin(), t.end(), [](const PVTGRecord& a, const PVTGRecord& b) { return a.pg < b.pg; });
        const int n = int(t.size());
        pg_.clear(); rvMax_.clear(); gB_.assign(size_t(n), {}); gMu_.assign(size_t(n), {}); rvOfP_ = {};
        for (int k = 0; k < n; ++k) {
            if (t[size_t(k)].rows.empty()) throw std::runtime_error("PVTG: empty record");
            pg_.push_back(t[size_t(k)].pg); rvMax_.push_back(t[size_t(k)].rows.front().rv);
            rvOfP_.x.push_back(t[size_t(k)].pg); rvOfP_.y.push_back(t[size_t(k)].rows.front().rv);
        }
        for (int k = 0; k < n; ++k) {
            std::vector<PVTGRow> rows = t[size_t(k)].rows;
            if (rows.size() < 2) {
                int src = -1;
                for (int m = k + 1; m < n && src < 0; ++m) if (t[size_t(m)].rows.size() >= 2) src = m;
                for (int m = k - 1; m >= 0 && src < 0; --m) if (t[size_t(m)].rows.size() >= 2) src = m;
                if (src >= 0) {
                    const auto& M = t[size_t(src)].rows;
                    for (size_t j = 1; j < M.size(); ++j) {
                        const PVTGRow last = rows.back();
                        const double xB = (M[j].B - M[j - 1].B) / (0.5 * (M[j].B + M[j - 1].B));
                        const double xM = (M[j].mu - M[j - 1].mu) / (0.5 * (M[j].mu + M[j - 1].mu));
                        rows.push_back({last.rv + (M[j].rv - M[j - 1].rv), last.B * (1 + xB / 2) / (1 - xB / 2), last.mu * (1 + xM / 2) / (1 - xM / 2)});
                    }
                }
            }
            std::sort(rows.begin(), rows.end(), [](const PVTGRow& a, const PVTGRow& b) { return a.rv < b.rv; });
            Table1D a, m;
            for (auto& r : rows) {
                if (!a.x.empty() && r.rv <= a.x.back() + 1e-14) continue;
                a.x.push_back(r.rv); a.y.push_back(1.0 / r.B);
                m.x.push_back(r.rv); m.y.push_back(r.mu);
            }
            gB_[size_t(k)] = a; gMu_[size_t(k)] = m;
        }
        // the dry-gas tables become the saturated curve (used for summaries / initial estimates)
        gasB_ = {}; gasBMu_ = {};
        for (auto& r : t) { gasB_.x.push_back(r.pg); gasB_.y.push_back(1.0 / r.rows.front().B);
                            gasBMu_.x.push_back(r.pg); gasBMu_.y.push_back(1.0 / (r.rows.front().B * r.rows.front().mu)); }
    }
    if (live_) {
        auto t = PR.pvto;
        if (t.size() == 1) {
            // a single saturated record (e.g. SPE3): Rs falls linearly to zero below its bubble point
            PVTORecord z = t[0]; z.rs = 0; z.pb = 0; z.under.clear();
            t.insert(t.begin(), z);
        }
        if (t.size() < 2) throw std::runtime_error("PVTO needs at least 1 saturated record");
        rs_.clear(); pb_.clear(); bSat_.clear(); bmuSat_.clear(); uB_.clear(); uBMu_.clear();
        rsOfP_ = {}; pbOfRs_ = {};
        for (auto& r : t) {
            rs_.push_back(r.rs); pb_.push_back(r.pb);
            bSat_.push_back(1.0 / r.bo); bmuSat_.push_back(1.0 / (r.bo * r.mu));
            rsOfP_.x.push_back(r.pb); rsOfP_.y.push_back(r.rs);
            pbOfRs_.x.push_back(r.rs); pbOfRs_.y.push_back(r.pb);
        }
        // undersaturated ratio curves; records without data borrow from the next record (or last with data)
        const int n = int(t.size());
        int lastWith = -1;
        for (int k = n - 1; k >= 0; --k) if (!t[k].under.empty()) { lastWith = k; break; }
        if (lastWith < 0) throw std::runtime_error("PVTO: undersaturated data must be given for the highest Rs");
        uB_.resize(n); uBMu_.resize(n);
        for (int k = 0; k < n; ++k) {
            int src = -1;
            for (int m = k; m < n; ++m) if (!t[m].under.empty()) { src = m; break; }
            if (src < 0) src = lastWith;
            const auto& R = t[src];
            Table1D a, b;
            a.x.push_back(0); a.y.push_back(1.0); b.x.push_back(0); b.y.push_back(1.0);
            const double b0 = 1.0 / R.bo, bm0 = 1.0 / (R.bo * R.mu);
            for (auto& u : R.under) {
                a.x.push_back(u.p - R.pb); a.y.push_back((1.0 / u.B) / b0);
                b.x.push_back(u.p - R.pb); b.y.push_back((1.0 / (u.B * u.mu)) / bm0);
            }
            uB_[k] = a; uBMu_[k] = b;
        }
    }
}

AD3 FluidModel::bw(const AD3& p) const {
    AD3 X = (p - pvtw_.pref) * pvtw_.cw;
    return (1.0 + X + X * X * 0.5) / pvtw_.bw;
}
AD3 FluidModel::bwMu(const AD3& p) const {
    AD3 Y = (p - pvtw_.pref) * (pvtw_.cw - pvtw_.cv);
    return (1.0 + Y + Y * Y * 0.5) / (pvtw_.bw * pvtw_.visc);
}
AD3 FluidModel::pvMult(const AD3& p) const {
    AD3 X = (p - rock_.pref) * rock_.cr;
    return 1.0 + X + X * X * 0.5;
}
double FluidModel::pvMult(double p) const { double X = rock_.cr * (p - rock_.pref); return 1 + X + 0.5 * X * X; }

AD3 FluidModel::bg(const AD3& p) const { return gasB_(p); }
AD3 FluidModel::bgMu(const AD3& p) const { return gasBMu_(p); }

double FluidModel::rvSat(double p, double* drvdp) const {
    if (!vap_) { if (drvdp) *drvdp = 0; return 0; }
    double d; double r = rvOfP_.eval(p, &d);
    if (r < 0) { r = 0; d = 0; }
    if (drvdp) *drvdp = d;
    return r;
}

// OPM Flow "right extreme" interpolation: the distance below the saturated Rv is kept constant while interpolating
// linearly in pressure between the two bracketing records (extrapolated beyond the table)
double FluidModel::wetGas(double p, double rv, int which, double* dp, double* drv) const {
    const int n = int(pg_.size());
    int k;
    if (p <= pg_[1]) k = 0;
    else if (p >= pg_[size_t(n - 2)]) k = n - 2;
    else k = int(std::upper_bound(pg_.begin(), pg_.end(), p) - pg_.begin()) - 1;
    const auto& T = which == 0 ? gB_ : gMu_;
    const double h = pg_[size_t(k + 1)] - pg_[size_t(k)];
    const double al = (p - pg_[size_t(k)]) / h;
    const double r0 = rvMax_[size_t(k)], r1 = rvMax_[size_t(k + 1)];
    const double d = rv - ((1 - al) * r0 + al * r1);
    double d0, d1;
    const double v0 = T[size_t(k)].eval(r0 + d, &d0), v1 = T[size_t(k + 1)].eval(r1 + d, &d1);
    const double dvdrv = (1 - al) * d0 + al * d1;
    if (drv) *drv = dvdrv;
    if (dp) *dp = ((v1 - v0) - (r1 - r0) * dvdrv) / h;
    return v0 + al * (v1 - v0);
}
AD3 FluidModel::bg(const AD3& p, const AD3& rv) const {
    if (!vap_) return gasB_(p);
    double fp, fr; const double f = wetGas(p.v, rv.v, 0, &fp, &fr);
    return chain2(p, rv, f, fp, fr);
}
AD3 FluidModel::bgMu(const AD3& p, const AD3& rv) const {
    if (!vap_) return gasBMu_(p);
    double bp, br, mp, mr;
    const double b = wetGas(p.v, rv.v, 0, &bp, &br), mu = wetGas(p.v, rv.v, 1, &mp, &mr);
    return chain2(p, rv, b / mu, (bp * mu - b * mp) / (mu * mu), (br * mu - b * mr) / (mu * mu));
}

double FluidModel::rsSat(double p, double* drsdp) const {
    if (!live_) { if (drsdp) *drsdp = 0; return 0; }
    double d; double r = rsOfP_.eval(p, &d);
    if (r < 0) { r = 0; d = 0; }
    if (drsdp) *drsdp = d;
    return r;
}
double FluidModel::pbOfRs(double rs) const { return live_ ? pbOfRs_.eval(rs) : 0.0; }

double FluidModel::liveOil(double p, double rs, int which) const {
    const int n = int(rs_.size());
    int k;
    if (rs <= rs_[1]) k = 0;
    else if (rs >= rs_[n - 2]) k = n - 2;
    else k = int(std::upper_bound(rs_.begin(), rs_.end(), rs) - rs_.begin()) - 1;
    const double pb = pbOfRs_.eval(rs);
    const double dp = p - pb;
    const auto& sat = which == 0 ? bSat_ : bmuSat_;
    const auto& ur  = which == 0 ? uB_ : uBMu_;
    const double v0 = sat[k] * ur[k].eval(dp);
    const double v1 = sat[k + 1] * ur[k + 1].eval(dp);
    const double w = (rs - rs_[k]) / (rs_[k + 1] - rs_[k]);
    return v0 + w * (v1 - v0);
}

AD3 FluidModel::liveOilAD(const AD3& p, const AD3& rs, int which) const {
    const double f = liveOil(p.v, rs.v, which);
    const double hp = 1e-4 * std::max(1.0, std::fabs(p.v));
    const double hr = 1e-5 * std::max(1e-2, std::fabs(rs.v));
    const double fp = (liveOil(p.v + hp, rs.v, which) - liveOil(p.v - hp, rs.v, which)) / (2 * hp);
    const double fr = (liveOil(p.v, rs.v + hr, which) - liveOil(p.v, rs.v - hr, which)) / (2 * hr);
    return chain2(p, rs, f, fp, fr);
}

AD3 FluidModel::bo(const AD3& p, const AD3& rs) const { return live_ ? liveOilAD(p, rs, 0) : deadB_(p); }
AD3 FluidModel::boMu(const AD3& p, const AD3& rs) const { return live_ ? liveOilAD(p, rs, 1) : deadBMu_(p); }

// ------------------------------------------------------------------ SatFunctions
void SatFunctions::setup(const Case& c, int region) {
    const SatRegion PR = c.props.satRegion(region);
    water_ = c.phases.water && c.phases.oil;
    gas_ = c.phases.gas && c.phases.oil;
    kroModel_ = c.props.kroModel; stoneEta_ = c.props.stone1Exp > 0 ? c.props.stone1Exp : 1.0;
    auto mk = [](const std::vector<SatRow>& t, int col) {
        Table1D tb; tb.clampEnds = true;
        for (auto& r : t) { tb.x.push_back(r.s); tb.y.push_back(col == 0 ? r.krw : col == 1 ? r.kro : r.pc); }
        return tb;
    };
    if (water_) {
        if (PR.swof.size() < 2) throw std::runtime_error("SWOF needs at least 2 rows");
        krw_ = mk(PR.swof, 0); krow_ = mk(PR.swof, 1); pcow_ = mk(PR.swof, 2);
        swco_ = PR.swof.front().s; swmax_ = PR.swof.back().s;
    } else { swco_ = 0; swmax_ = 0; }
    if (gas_) {
        if (PR.sgof.size() < 2) throw std::runtime_error("SGOF needs at least 2 rows");
        krg_ = mk(PR.sgof, 0); krog_ = mk(PR.sgof, 1); pcog_ = mk(PR.sgof, 2);
        sgmax_ = PR.sgof.back().s;
    }
    // table end points: critical = largest saturation with kr = 0 before kr rises (ECLIPSE definition)
    ep_ = EndPoints{};
    if (water_) {
        const auto& t = PR.swof;
        ep_.swl = t.front().s; ep_.swu = t.back().s; ep_.swcr = ep_.swl;
        for (size_t i = 0; i < t.size() && t[i].krw <= 0; ++i) ep_.swcr = t[i].s;
        double swAtKroZero = t.back().s;                       // lowest Sw with krow = 0
        for (size_t i = t.size(); i-- > 0;) { if (t[i].kro <= 0) swAtKroZero = t[i].s; else break; }
        ep_.sowcr = std::max(0.0, 1.0 - swAtKroZero);
    }
    if (gas_) {
        const auto& t = PR.sgof;
        ep_.sgl = t.front().s; ep_.sgu = t.back().s; ep_.sgcr = ep_.sgl;
        for (size_t i = 0; i < t.size() && t[i].krw <= 0; ++i) ep_.sgcr = t[i].s;
        double sgAtKroZero = t.back().s;
        for (size_t i = t.size(); i-- > 0;) { if (t[i].kro <= 0) sgAtKroZero = t[i].s; else break; }
        ep_.sogcr = std::max(0.0, 1.0 - ep_.swl - sgAtKroZero);
    }
}

// ---- end-point scaling (two-point: [critical, maximum] of the cell -> [critical, maximum] of the table)
namespace {
inline AD3 mapLin(const AD3& s, double a, double b, double at, double bt) {
    if (std::fabs(b - a) < 1e-9 || std::fabs(bt - at) < 1e-12) return s + (at - a);
    return at + (s - a) * ((bt - at) / (b - a));
}
inline double mapBack(double st, double a, double b, double at, double bt) {   // table -> cell (inverse of mapLin)
    if (std::fabs(bt - at) < 1e-12 || std::fabs(b - a) < 1e-9) return st - at + a;
    return a + (st - at) * ((b - a) / (bt - at));
}
}
AD3 SatFunctions::krw(const AD3& sw, const EndPoints* e) const {
    if (!e) return krw_(sw);
    return krw_(mapLin(sw, e->swcr, e->swu, ep_.swcr, ep_.swu));
}
AD3 SatFunctions::krg(const AD3& sg, const EndPoints* e) const {
    if (!e) return krg_(sg);
    return krg_(mapLin(sg, e->sgcr, e->sgu, ep_.sgcr, ep_.sgu));
}
AD3 SatFunctions::pcow(const AD3& sw, const EndPoints* e) const {
    if (!e) return pcow_(sw);
    return pcow_(mapLin(sw, e->swl, e->swu, ep_.swl, ep_.swu));
}
AD3 SatFunctions::pcog(const AD3& sg, const EndPoints* e) const {
    if (!e) return pcog_(sg);
    return pcog_(mapLin(sg, e->sgl, e->sgu, ep_.sgl, ep_.sgu));
}
AD3 SatFunctions::kro(const AD3& sw, const AD3& sg, const EndPoints* e) const {
    if (!e) return kro(sw, sg);
    const AD3 so = 1.0 - sw - sg;
    // krow(So): So in [SOWCR, 1 - SWL] -> table range; evaluated at Sw' = 1 - So'
    auto krowS = [&](const AD3& soIn) {
        const AD3 soT = mapLin(soIn, e->sowcr, 1.0 - e->swl, ep_.sowcr, 1.0 - ep_.swl);
        return krow_(1.0 - soT);
    };
    // krog(So): So in [SOGCR, 1 - SWL - SGL] -> table range; evaluated at Sg' = 1 - SWL_t - So'
    auto krogS = [&](const AD3& soIn) {
        const AD3 soT = mapLin(soIn, e->sogcr, 1.0 - e->swl - e->sgl, ep_.sogcr, 1.0 - ep_.swl - ep_.sgl);
        return krog_(1.0 - ep_.swl - soT);
    };
    if (water_ && !gas_) return krowS(so);
    if (!water_ && gas_) return krogS(so);
    if (!water_ && !gas_) return AD3(1.0);
    if (kroModel_ != 0) {
        // Stone models (OPM EclStone1/2Material): krow at Sw, krog at the gas-oil oil saturation 1 - Swl - Sg
        const AD3 krocw = krowS(AD3(1.0 - e->swl));
        return stone(sw, sg, krowS(1.0 - sw), krogS(1.0 - e->swl - sg), krw(sw, e), krg(sg, e), krocw.v, e->swl);
    }
    AD3 dw = sw - e->swl;
    if (dw.v < 0) dw = AD3(0.0);
    const AD3 den = sg + dw;
    if (den.v < 1e-9) return krowS(so);
    return (sg * krogS(so) + dw * krowS(so)) / den;
}
double SatFunctions::swInverseFromPc(double pc, const EndPoints* e) const {
    const double s = swInverseFromPc(pc);
    if (!e || !water_) return s;
    return std::clamp(mapBack(s, e->swl, e->swu, ep_.swl, ep_.swu), e->swl, e->swu);
}
double SatFunctions::sgInverseFromPc(double pc, const EndPoints* e) const {
    const double s = sgInverseFromPc(pc);
    if (!e || !gas_) return s;
    return std::clamp(mapBack(s, e->sgl, e->sgu, ep_.sgl, ep_.sgu), e->sgl, e->sgu);
}
double SatFunctions::swInverseFromPcGW(double pc, const EndPoints* e) const {
    if (!e) return swInverseFromPcGW(pc);
    if (!water_ || !gas_) return swInverseFromPc(pc, e);
    auto f = [&](double sw) { return pcow(AD3(sw), e).v + pcog(AD3(1.0 - sw), e).v - pc; };
    double lo = e->swl, hi = e->swu;
    if (f(lo) <= 0) return lo;
    if (f(hi) >= 0) return hi;
    for (int it = 0; it < 60; ++it) { const double m = 0.5 * (lo + hi); (f(m) > 0 ? lo : hi) = m; }
    return 0.5 * (lo + hi);
}

AD3 SatFunctions::kro(const AD3& sw, const AD3& sg) const {
    if (water_ && !gas_) return krow_(sw);
    if (!water_ && gas_) return krog_(sg);
    if (!water_ && !gas_) return AD3(1.0);
    // ECLIPSE default model: kro = (Sg*krog + (Sw-Swco)*krow) / (Sg + Sw - Swco), where both two-phase
    // curves are evaluated at the actual oil saturation So: krow(Sw' = 1-So), krog(Sg' = 1-So-Swco).
    // This guarantees kro -> 0 as So -> residual / zero.
    if (kroModel_ != 0) return stone(sw, sg, krow_(sw), krog_(sg), krw_(sw), krg_(sg), krow_(AD3(swco_)).v, swco_);
    AD3 dw = sw - swco_;
    if (dw.v < 0) dw = AD3(0.0);
    const AD3 swEff = sw + sg;                 // 1 - So
    const AD3 sgEff = sw + sg - swco_;         // 1 - So - Swco
    AD3 den = sg + dw;
    if (den.v < 1e-9) return krow_(swEff);
    return (sg * krog_(sgEff) + dw * krow_(swEff)) / den;
}

AD3 SatFunctions::stone(const AD3& sw, const AD3& sg, const AD3& krow, const AD3& krog, const AD3& krw, const AD3& krg,
                        double krocw, double swco) const {
    if (krocw <= 1e-12) return AD3(0.0);
    if (kroModel_ == 2) {   // STONE2: kro = krocw ((krow/krocw + krw)(krog/krocw + krg) - krw - krg)
        AD3 k = krocw * ((krow / krocw + krw) * (krog / krocw + krg) - krw - krg);
        return k.v < 0 ? AD3(0.0) : k;
    }
    // STONE1 (Aziz & Settari normalisation, Swco only): kro = beta krow krog / krocw
    AD3 beta(1.0);
    if (sw.v > swco) {
        const AD3 ssw = (sw - swco) / (1.0 - swco), ssg = sg / (1.0 - swco), sso = 1.0 - ssw - ssg;
        if (ssw.v < 1.0 && ssg.v < 1.0) {
            beta = sso / ((1.0 - ssw) * (1.0 - ssg));
            if (stoneEta_ != 1.0) beta = beta.v > 0 ? exp(stoneEta_ * log(beta)) : AD3(0.0);
        }
    }
    AD3 k = beta * krow * krog / krocw;
    if (k.v < 0) return AD3(0.0);
    if (k.v > 1) return AD3(1.0);
    return k;
}

double SatFunctions::swInverseFromPc(double pc) const {
    if (!water_) return 0;
    // Pcow decreasing with Sw
    const auto& t = pcow_;
    if (pc >= t.y.front()) return t.x.front();
    if (pc <= t.y.back()) return t.x.back();
    return t.inverse(pc);
}
// gas directly on water (overlapping transition zones): Sw with Pcow(Sw) + Pcog(1 - Sw) = Pcgw (OPM / ECLIPSE rule)
double SatFunctions::swInverseFromPcGW(double pc) const {
    if (!water_ || !gas_) return swInverseFromPc(pc);
    auto f = [&](double sw) { return pcow_.eval(sw) + pcog_.eval(1.0 - sw) - pc; };   // decreasing in Sw
    double lo = swco_, hi = swmax_;
    if (f(lo) <= 0) return lo;
    if (f(hi) >= 0) return hi;
    for (int it = 0; it < 60; ++it) { const double m = 0.5 * (lo + hi); (f(m) > 0 ? lo : hi) = m; }
    return 0.5 * (lo + hi);
}

double SatFunctions::sgInverseFromPc(double pc) const {
    if (!gas_) return 0;
    const auto& t = pcog_;          // increasing with Sg
    if (pc <= t.y.front()) return t.x.front();
    if (pc >= t.y.back()) return t.x.back();
    return t.inverse(pc);
}

} // namespace ressim
