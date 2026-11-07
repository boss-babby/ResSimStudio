// C24: fault mechanics and induced seismicity inside the flow run - see faultmech.h
#include "ressim/faultmech.h"
#include "ressim/case.h"
#include "ressim/ddm.h"
#include "ressim/grid.h"
#include "ressim/rsf.h"
#include "ressim/simulator.h"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <map>
#include <random>

namespace ressim {
namespace faultmech {
namespace {
constexpr double kYearDays = 365.25;
using V3 = std::array<double, 3>;
V3 sub(const std::array<float, 3>& a, const std::array<float, 3>& b) { return {double(a[0]) - b[0], double(a[1]) - b[1], double(a[2]) - b[2]}; }
V3 cross(const V3& a, const V3& b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }
double norm(const V3& a) { return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }
// corner indices (bit 0 = +x, bit 1 = +y, bit 2 = +z / deeper) of a cell face, in loop order
const int kFace[6][4] = {{0, 2, 6, 4}, {1, 3, 7, 5}, {0, 1, 5, 4}, {2, 3, 7, 6}, {0, 1, 3, 2}, {4, 5, 7, 6}};
// moment magnitude of a full-fault circular rupture with stress drop dSigma (Pa): M0 = 16/7 dSigma (A / pi)^1.5
double mwFullRupture(double areaM2, double dSigma = 3e6) {
    const double M0 = 16.0 / 7.0 * dSigma * std::pow(std::max(1.0, areaM2) / M_PI, 1.5);
    return 2.0 / 3.0 * (std::log10(M0) - 9.1);
}
}  // namespace

bool Monitor::setup(const Case& c, const Grid& g, std::string& err) {
    active_ = false; out_ = Output{};
    const FaultMechInput& F = c.faultMech;
    if (!F.on()) { err = "no STRESSFIELD / FAULTMECH"; return false; }
    const bool fld = c.units == UnitSystem::Field;
    out_.paPerUnit = fld ? 6894.757 : 1e5; out_.mPerUnit = fld ? 0.3048 : 1.0;
    out_.engine = c.thermal.enabled ? "thermal" : c.comp.enabled ? (c.comp.co2store ? "CO2 storage" : "compositional") : "black oil";
    source_ = F.source;
    thCoef_ = c.thermal.enabled && F.thE > 0 && F.thAlpha > 0 ? F.thE * 1e9 * F.thAlpha / (1 - F.thNu) / out_.paPerUnit : 0.0;
    out_.thermalStress = thCoef_ > 0 && source_ == 0;
    svGrad_ = F.svGrad; sv0_ = F.sv0; rh_ = F.shRatio; rH_ = F.sHRatio; az_ = F.azimuth * M_PI / 180; path_ = F.stressPath;
    const int nx = g.nx, ny = g.ny, nz = g.nz, nCells = g.numCells();
    std::string missing;
    for (const FaultMechInput::Row& R : F.faults) {
        const FaultInput* FI = nullptr;
        for (const auto& q : c.grid.faults) if (q.name == R.name) FI = &q;
        if (!FI) { missing += (missing.empty() ? "" : ", ") + R.name; continue; }
        const int fi = int(out_.names.size());
        out_.names.push_back(R.name);
        mu_.push_back(R.mu); coh_.push_back(R.coh); r0_.push_back(R.r0); ta_.push_back(R.ta); perm_.push_back(R.permFactor); ramp_.push_back(R.permRamp);
        bval_.push_back(R.bValue); mmin_.push_back(R.mMin);
        rsfLayer_.push_back(R.rsfLayer); rsfA_.push_back(R.rsfA); rsfB_.push_back(R.rsfB); rsfDc_.push_back(R.rsfDcMm * 1e-3);
        if (R.permFactor != 1.0) twoWay_ = true;
        std::map<std::pair<int, int>, bool> seen;   // (global cell, face id) - FAULTS records may overlap
        for (const FaultFace& ff : FI->faces) {
            const int fid = (ff.face == 'X' ? 0 : ff.face == 'Y' ? 2 : 4) + (ff.minus ? 0 : 1);
            const int di = ff.face == 'X' ? (ff.minus ? -1 : 1) : 0, dj = ff.face == 'Y' ? (ff.minus ? -1 : 1) : 0, dk = ff.face == 'Z' ? (ff.minus ? -1 : 1) : 0;
            for (int k = std::max(0, ff.k1); k <= std::min(nz - 1, ff.k2); ++k)
                for (int j = std::max(0, ff.j1); j <= std::min(ny - 1, ff.j2); ++j)
                    for (int i = std::max(0, ff.i1); i <= std::min(nx - 1, ff.i2); ++i) {
                        const int gl = i + nx * (j + ny * k);
                        if (gl >= nCells || size_t(gl) >= g.corners.size()) continue;
                        const int a = g.activeOf[size_t(gl)];
                        if (a < 0 || seen[{gl, fid}]) continue;
                        seen[{gl, fid}] = true;
                        Face f; f.fault = fi; f.g = gl; f.a = a; f.k = k;
                        const auto& C8 = g.corners[size_t(gl)];
                        for (int v = 0; v < 4; ++v) { const Vec3& q = C8[size_t(kFace[fid][v])]; f.quad[size_t(v)] = {q.x, q.y, q.z}; }
                        for (int d = 0; d < 3; ++d) f.c[size_t(d)] = 0.25f * (f.quad[0][size_t(d)] + f.quad[1][size_t(d)] + f.quad[2][size_t(d)] + f.quad[3][size_t(d)]);
                        const V3 nv = cross(sub(f.quad[2], f.quad[0]), sub(f.quad[3], f.quad[1]));
                        const double ln = norm(nv);
                        f.area = 0.5 * ln;
                        if (!(ln > 0)) continue;
                        for (int d = 0; d < 3; ++d) f.n[size_t(d)] = float(nv[size_t(d)] / ln);
                        // connections across: neighbours in the adjacent column (any layer: NNCs of displaced columns) / layer
                        double bestT = -1;
                        for (int q = g.adjStart[size_t(a)]; q < g.adjStart[size_t(a) + 1]; ++q) {
                            const int nb = g.adjNbr[size_t(q)], cn = g.adjConn[size_t(q)];
                            const int gb = g.globalOf[size_t(nb)];
                            if (gb >= nCells) continue;
                            const int ib = gb % nx, jb = (gb / nx) % ny, kb = gb / (nx * ny);
                            const bool across = dk != 0 ? (ib == i && jb == j && kb == k + dk) : (ib == i + di && jb == j + dj);
                            if (!across) continue;
                            f.conns.push_back(cn);
                            const double T = g.conns[size_t(cn)].trans;
                            if (T > bestT) { bestT = T; f.b = nb; }
                        }
                        out_.faces.push_back(std::move(f));
                    }
        }
    }
    if (out_.faces.empty()) {
        err = missing.empty() ? "the FAULTMECH faults have no faces on active cells" : "FAULTMECH names without FAULTS: " + missing;
        return false;
    }
    out_.series.assign(out_.names.size(), Series{});
    out_.mMax.assign(out_.names.size(), 0.0);
    {
        std::vector<double> A(out_.names.size(), 0.0);
        for (const Face& f : out_.faces) A[size_t(f.fault)] += f.area * out_.mPerUnit * out_.mPerUnit;
        for (size_t q = 0; q < A.size(); ++q) out_.mMax[q] = mwFullRupture(A[q]);
    }
    // Dieterich A sigma per face (set at t0 from sigma'_n when not given)
    for (Face& f : out_.faces) {
        const FaultMechInput::Row* R = nullptr;
        for (const auto& q : F.faults) if (q.name == out_.names[size_t(f.fault)]) R = &q;
        f.aSig = R ? R->aSig : 0;
    }
    gamma_.assign(out_.faces.size(), 0.0);
    mult_.assign(out_.faces.size(), 1.0);
    if (!missing.empty()) out_.error = "FAULTMECH names without FAULTS (ignored): " + missing;
    setupSlip(c);
    active_ = true; out_.ok = true;
    return true;
}

// ---------------------------------------------------------------------------------------------------------- C26 DDM slip
void Monitor::setupSlip(const Case& c) {
    const FaultMechInput& F = c.faultMech;
    slipOn_ = F.slipOn; out_.slip = slipOn_;
    ddm_.reset(); slip_.clear(); rake_.clear(); failed_.clear();
    if (!slipOn_) return;
    out_.shearModulus = F.slipE * 1e9 / (2 * (1 + F.slipNu));
    slipNu_ = std::clamp(F.slipNu, 0.0, 0.49); dMu_ = std::max(0.0, F.slipDMu);
    slipRef_ = std::max(1e-6, F.slipRefMm * 1e-3); slipCut_ = F.slipCutoff;
    out_.halfSpace = F.slipHalfSpace != 0;
    // half space only if the whole fault set lies below the datum depth 0
    for (const Face& f : out_.faces) for (const auto& q : f.quad) if (q[2] < 0) out_.halfSpace = false;
    slip_.assign(out_.faces.size(), 0.0);
    failed_.assign(out_.faces.size(), 0);
    if (out_.faces.size() > 12000) {
        out_.slipNote = "FAULTSLIP: " + std::to_string(out_.faces.size()) + " fault faces - interaction matrix limited by the cut-off (" +
                        std::to_string(int(slipCut_)) + " element sizes)";
        if (slipCut_ <= 0 || slipCut_ > 25) slipCut_ = 25;
    }
}

void Monitor::slipStep(const std::vector<double>& pf, const std::vector<double>& snR, const std::vector<std::array<double, 3>>& shR,
                       std::vector<double>& sn, std::vector<double>& tau, StepOut& S, double tDays, double tPrev) {
    const size_t nf = out_.faces.size();
    const double Pa = out_.paPerUnit, M = out_.mPerUnit, G = out_.shearModulus;
    S.slip.assign(nf, 0.0f); S.dtauSlip.assign(nf, 0.0f);
    const bool first = out_.steps.empty();
    // computation frame: metres, z up (x, y, -depth) - a reflection; vectors transform alike, scalars are invariant
    auto toC = [&](const std::array<float, 3>& q) { return ddm::V3{double(q[0]) * M, double(q[1]) * M, -double(q[2]) * M}; };
    if (first || !ddm_) {
        rake_.assign(nf, {1, 0, 0});
        std::vector<ddm::Element> el(nf);
        std::vector<ddm::V3> dir(nf);
        for (size_t i = 0; i < nf; ++i) {
            const Face& f = out_.faces[i];
            const ddm::V3 nref = {f.n[0], f.n[1], -double(f.n[2])};
            el[i] = ddm::quadElement(toC(f.quad[0]), toC(f.quad[1]), toC(f.quad[2]), toC(f.quad[3]));
            if (el[i].n[0] * nref[0] + el[i].n[1] * nref[1] + el[i].n[2] * nref[2] < 0)
                el[i] = ddm::quadElement(toC(f.quad[0]), toC(f.quad[3]), toC(f.quad[2]), toC(f.quad[1]));
            const double t = std::sqrt(shR[i][0] * shR[i][0] + shR[i][1] * shR[i][1] + shR[i][2] * shR[i][2]);
            if (t > 0) rake_[i] = {shR[i][0] / t, shR[i][1] / t, shR[i][2] / t};
            else {   // no resolved shear: dip direction
                const double nx = f.n[0], ny = f.n[1], nz = f.n[2];
                double d[3] = {-nx * nz, -ny * nz, nx * nx + ny * ny};
                const double l = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
                rake_[i] = l > 0 ? std::array<double, 3>{d[0] / l, d[1] / l, d[2] / l} : std::array<double, 3>{1, 0, 0};
            }
            dir[i] = {rake_[i][0], rake_[i][1], -rake_[i][2]};
        }
        auto K = std::make_shared<ddm::ShearMatrix>(ddm::buildShear(el, dir, G, slipNu_, out_.halfSpace, slipCut_));
        out_.ddmPairs = K->pairs; out_.ddmDropped = K->dropped; out_.ddmSeconds = K->seconds;
        ddm_ = K;
    }
    const ddm::ShearMatrix& K = *static_cast<const ddm::ShearMatrix*>(ddm_.get());
    // slip-induced stress (Pa) along the rake and normal (compression +) from the cumulative slip
    auto induced = [&](std::vector<double>& dt, std::vector<double>& dn) {
        dt.assign(nf, 0.0); dn.assign(nf, 0.0);
        for (size_t i = 0; i < nf; ++i)
            for (int q = K.start[i]; q < K.start[i + 1]; ++q) {
                const double s = slip_[size_t(K.col[size_t(q)])];
                if (s == 0) continue;
                dt[i] += K.kt[size_t(q)] * s; dn[i] += K.kn[size_t(q)] * s;
            }
    };
    std::vector<double> dt, dn;
    if (!first) {
        induced(dt, dn);
        // excess over the strength (Pa): shear along the rake - mu sigma'_n - c - max(0, CFF0) (faces critical at t0 are taken
        // as being at their strength, so slip is driven by the change since t0)
        std::vector<double> ex(nf), muS(nf), muUse(nf), ds(nf, 0.0);
        auto excess = [&](size_t i, double mu) {
            const Face& f = out_.faces[i];
            const size_t fi = size_t(f.fault);
            const double tauD = (shR[i][0] * rake_[i][0] + shR[i][1] * rake_[i][1] + shR[i][2] * rake_[i][2]) * Pa + dt[i];
            const double sne = (snR[i] - pf[i]) * Pa + dn[i];
            return tauD - mu * std::max(0.0, sne) - coh_[fi] * Pa - std::max(0.0, f.cff0) * Pa;
        };
        std::fill(failed_.begin(), failed_.end(), 0);   // friction heals between report steps
        int iterTot = 0;
        for (int pass = 0; pass < 60; ++pass) {
            for (size_t i = 0; i < nf; ++i) {
                muS[i] = mu_[size_t(out_.faces[i].fault)];
                muUse[i] = failed_[i] ? std::max(0.0, muS[i] - dMu_) : muS[i];
                ex[i] = excess(i, muUse[i]);
            }
            std::fill(ds.begin(), ds.end(), 0.0);
            const ddm::SlipResult r = ddm::solveSlip(K, ex, muUse, ds, 1e-3 * Pa, slipIterMax_);
            iterTot += r.iterations;
            if (dMu_ <= 0) break;
            bool grew = false;   // weakening cascade: faces that slipped now carry the dynamic friction
            for (size_t i = 0; i < nf; ++i) if (ds[i] > 0 && !failed_[i]) { failed_[i] = 1; grew = true; }
            if (!grew) break;
        }
        out_.maxSlipIter = std::max(out_.maxSlipIter, iterTot);
        double m0 = 0;
        std::vector<double> mF(out_.names.size(), 0.0), xs(out_.names.size(), 0.0), ys(out_.names.size(), 0.0), zs(out_.names.size(), 0.0);
        std::vector<int> nF(out_.names.size(), 0);
        for (size_t i = 0; i < nf; ++i) {
            if (!(ds[i] > 0)) continue;
            slip_[i] += ds[i];
            const Face& f = out_.faces[i];
            const double dm = G * f.area * M * M * ds[i];
            m0 += dm; ++S.slipping;
            const size_t q = size_t(f.fault);
            mF[q] += dm; xs[q] += dm * f.c[0]; ys[q] += dm * f.c[1]; zs[q] += dm * f.c[2]; ++nF[q];
        }
        S.momentStep = m0;
        if (dMu_ > 0)   // seismic slip events (weakening): one per fault and report step, at the moment centroid
            for (size_t q = 0; q < mF.size(); ++q) {
                if (!(mF[q] > 0)) continue;
                Event ev; ev.kind = 2; ev.fault = int(q); ev.faces = nF[q];
                ev.t = 0.5 * (tPrev + tDays);
                ev.x = float(xs[q] / mF[q]); ev.y = float(ys[q] / mF[q]); ev.z = float(zs[q] / mF[q]);
                ev.mw = float(2.0 / 3.0 * (std::log10(mF[q]) - 9.1));
                out_.events.push_back(ev);
            }
        induced(dt, dn);
    } else {
        dt.assign(nf, 0.0); dn.assign(nf, 0.0);
    }
    S.moment = (out_.steps.empty() ? 0.0 : out_.steps.back().moment) + S.momentStep;
    for (size_t i = 0; i < nf; ++i) {
        double v[3];
        for (int d = 0; d < 3; ++d) v[d] = shR[i][size_t(d)] + dt[i] / Pa * rake_[i][size_t(d)];
        tau[i] = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        sn[i] = snR[i] + dn[i] / Pa;
        S.slip[i] = float(slip_[i]); S.dtauSlip[i] = float(dt[i] / Pa);
    }
}

void Monitor::stress(const Face& f, double p, double dT, const double* dS6, double& sn, double& tau, double* shear) const {
    const double z = f.c[2];
    const double Sv = sv0_ + svGrad_ * z;
    const double dp = p - f.p0;
    // in-situ field + change: poroelastic path (+ thermo-elastic horizontal stress), or the FE tensor change (C25)
    const double dTh = dS6 ? 0.0 : thCoef_ * dT;
    const double dH = dS6 ? 0.0 : path_ * dp + dTh;
    const double SH = rH_ * Sv + dH, Sh = rh_ * Sv + dH;
    const double e1[3] = {std::cos(az_), std::sin(az_), 0}, e2[3] = {-std::sin(az_), std::cos(az_), 0};
    double S[3][3];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) S[i][j] = SH * e1[i] * e1[j] + Sh * e2[i] * e2[j] + (i == 2 && j == 2 ? Sv : 0.0);
    if (dS6) {
        const double D[3][3] = {{dS6[0], dS6[3], dS6[5]}, {dS6[3], dS6[1], dS6[4]}, {dS6[5], dS6[4], dS6[2]}};
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) S[i][j] += D[i][j];
    }
    double t[3] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) t[i] += S[i][j] * f.n[size_t(j)];
    sn = t[0] * f.n[0] + t[1] * f.n[1] + t[2] * f.n[2];
    const double sx = t[0] - sn * f.n[0], sy = t[1] - sn * f.n[1], sz = t[2] - sn * f.n[2];
    tau = std::sqrt(sx * sx + sy * sy + sz * sz);
    if (shear) { shear[0] = sx; shear[1] = sy; shear[2] = sz; }
}

void Monitor::step(double tDays, const std::vector<float>& p, const std::vector<float>* temperature) {
    if (!active_) return;
    const size_t nf = out_.faces.size();
    const bool first = out_.steps.empty();
    StepOut S; S.t = tDays;
    S.cff.resize(nf); S.dcfs.resize(nf); S.st.resize(nf); S.sne.resize(nf); S.rate.resize(nf); S.cum.resize(nf); S.mult.resize(nf);
    const StepOut* P = first ? nullptr : &out_.steps.back();
    const double dtYr = first ? 0.0 : std::max(0.0, (tDays - P->t) / kYearDays);
    std::mt19937 rng(seed_ + unsigned(out_.steps.size()) * 7919u);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    const size_t nflt = out_.names.size();
    std::vector<double> maxD(nflt, -1e300), maxC(nflt, -1e300), maxSt(nflt, 0), react(nflt, 0), area(nflt, 0), expct(nflt, 0);
    std::vector<double> nEv(nflt, 0), maxMw(nflt, -10);
    if (P) for (size_t q = 0; q < nflt; ++q) {   // carry cumulative counts
        const Series& R = out_.series[q];
        if (!R.events.empty()) nEv[q] = R.events.back();
        if (!R.maxMw.empty()) maxMw[q] = R.maxMw.back();
    }
    const double areaKm2 = out_.mPerUnit * out_.mPerUnit * 1e-6;
    std::vector<double> pfA(nf), snR(nf), snA(nf), tauA(nf);
    std::vector<std::array<double, 3>> shR(nf);
    for (size_t i = 0; i < nf; ++i) {
        Face& f = out_.faces[i];
        const double pa = size_t(f.a) < p.size() ? p[size_t(f.a)] : 0.0;
        const double pb = f.b >= 0 && size_t(f.b) < p.size() ? p[size_t(f.b)] : pa;
        const double pf = 0.5 * (pa + pb);
        if (first) f.p0 = pf;
        double Tf = 0;
        if (temperature) {
            const double ta = size_t(f.a) < temperature->size() ? (*temperature)[size_t(f.a)] : 0.0;
            const double tb = f.b >= 0 && size_t(f.b) < temperature->size() ? (*temperature)[size_t(f.b)] : ta;
            Tf = 0.5 * (ta + tb);
            if (first) f.T0 = Tf;
        }
        // C25 FE stress change of the face (mean of the two face cells, Pa -> case units)
        double dS6[6] = {0, 0, 0, 0, 0, 0}; bool useFe = false;
        if (source_ == 1 && fe_) {
            double Sa[6], Sb[6];
            if (fe_(f.a, Sa)) {
                if (!(f.b >= 0 && fe_(f.b, Sb))) for (int r = 0; r < 6; ++r) Sb[r] = Sa[r];
                double S6[6]; for (int r = 0; r < 6; ++r) S6[r] = 0.5 * (Sa[r] + Sb[r]) / out_.paPerUnit;
                if (first || !f.feOk) { for (int r = 0; r < 6; ++r) f.fe0[size_t(r)] = S6[r]; f.feOk = true; }
                for (int r = 0; r < 6; ++r) dS6[r] = S6[r] - f.fe0[size_t(r)];
                useFe = true;
            }
        }
        double sn = 0, tau = 0;
        stress(f, pf, temperature ? Tf - f.T0 : 0.0, useFe ? dS6 : nullptr, sn, tau, shR[i].data());
        pfA[i] = pf; snR[i] = sn; snA[i] = sn; tauA[i] = tau;
    }
    // C26: elastic stress transfer between the faces and quasi-static frictional slip (FAULTSLIP)
    if (slipOn_) slipStep(pfA, snR, shR, snA, tauA, S, tDays, P ? P->t : tDays);
    for (size_t i = 0; i < nf; ++i) {
        Face& f = out_.faces[i];
        const size_t fi = size_t(f.fault);
        const double pf = pfA[i], sn = snA[i], tau = tauA[i];
        const double sne = sn - pf;
        const double cff = tau - mu_[fi] * sne - coh_[fi];
        if (first) {
            f.sn0 = sn; f.tau0 = tau; f.cff0 = cff;
            if (!(f.aSig > 0)) f.aSig = 0.003 * std::max(sne, 0.0);
            if (!(f.aSig > 0)) f.aSig = 1e-3;
            gamma_[i] = ta_[fi] / f.aSig;   // steady background: gamma = 1 / Sr
            if (cff >= 0) ++out_.criticalAtStart;
        }
        const double Sr = f.aSig / ta_[fi];                // pressure units per year
        const double r = r0_[fi] * f.area * areaKm2;       // events per year (background)
        double nExp = 0;
        if (!first && dtYr > 0) {
            // constant stressing rate over the report step: tectonic background + simulated Coulomb change
            const double Sdot = Sr + (cff - (double(P->dcfs[i]) + f.cff0)) / dtYr;
            const int ns = 16;
            const double h = dtYr / ns;
            double g = gamma_[i];
            auto adv = [&](double gg, double hh) {
                if (std::fabs(Sdot) * hh / f.aSig < 1e-9) return gg + hh / f.aSig * (1 - gg * Sdot);
                const double e = std::exp(std::clamp(-Sdot * hh / f.aSig, -700.0, 700.0));
                return (gg - 1 / Sdot) * e + 1 / Sdot;
            };
            double acc = 0;
            for (int s = 0; s < ns; ++s) {   // Simpson on R = r / (gamma Sr) over each sub-step
                const double gm = std::max(1e-300, adv(g, 0.5 * h)), g1 = std::max(1e-300, adv(g, h));
                acc += h / 6 * (r / (std::max(1e-300, g) * Sr) + 4 * r / (gm * Sr) + r / (g1 * Sr));
                g = g1;
            }
            gamma_[i] = g;
            nExp = std::isfinite(acc) ? acc : 0;
        }
        const double cumPrev = P ? double(P->cum[i]) : 0.0;
        S.cff[i] = float(cff); S.dcfs[i] = float(cff - f.cff0); S.st[i] = float(sne > 0 ? tau / sne : 99.0); S.sne[i] = float(sne);
        S.rate[i] = float(1.0 / std::max(1e-300, gamma_[i] * Sr));
        S.cum[i] = float(cumPrev + nExp);
        // two-way coupling: irreversible transmissibility enhancement of reactivated faces
        if (perm_[fi] != 1.0) {
            double x = std::clamp(cff / std::max(1e-12, ramp_[fi] * std::max(sne, 1e-9)), 0.0, 1.0);   // C49: ramp width item 10
            if (slipOn_ && i < slip_.size()) x = std::max(x, std::clamp(slip_[i] / std::max(1e-9, slipRef_), 0.0, 1.0));   // C26
            mult_[i] = std::max(mult_[i], std::pow(perm_[fi], x));
        }
        S.mult[i] = float(mult_[i]);
        // stochastic realisation
        if (nExp > 0) {
            int k = 0;
            if (nExp < 50) { std::poisson_distribution<int> PD(nExp); k = PD(rng); }
            else k = std::max(0, int(std::lround(nExp + std::sqrt(nExp) * std::normal_distribution<double>(0, 1)(rng))));
            const double b = bval_[fi], mlo = mmin_[fi], mhi = std::max(mlo + 0.1, out_.mMax[fi]);
            const double cdfMax = 1 - std::pow(10.0, -b * (mhi - mlo));
            for (int e = 0; e < k && out_.events.size() < 50000; ++e) {
                Event ev; ev.fault = int(fi); ev.kind = 0;
                ev.t = P->t + U(rng) * (tDays - P->t);
                const double u = U(rng), v = U(rng);
                for (int d = 0; d < 3; ++d) {
                    const double x0 = f.quad[0][size_t(d)] + u * (f.quad[1][size_t(d)] - f.quad[0][size_t(d)]);
                    const double x1 = f.quad[3][size_t(d)] + u * (f.quad[2][size_t(d)] - f.quad[3][size_t(d)]);
                    (d == 0 ? ev.x : d == 1 ? ev.y : ev.z) = float(x0 + v * (x1 - x0));
                }
                ev.mw = float(mlo - std::log10(1 - U(rng) * cdfMax) / b);
                maxMw[fi] = std::max(maxMw[fi], double(ev.mw));
                out_.events.push_back(ev);
                nEv[fi] += 1;
            }
        }
        maxD[fi] = std::max(maxD[fi], double(S.dcfs[i])); maxC[fi] = std::max(maxC[fi], cff); maxSt[fi] = std::max(maxSt[fi], double(S.st[i]) < 90 ? double(S.st[i]) : 0.0);
        area[fi] += f.area; if (cff >= 0) react[fi] += f.area;
        expct[fi] += S.cum[i];
    }
    for (size_t q = 0; q < nflt; ++q) {
        if (!(area[q] > 0)) { maxD[q] = maxC[q] = 0; }
        Series& R = out_.series[q];
        R.t.push_back(tDays); R.maxDcfs.push_back(maxD[q]); R.maxCff.push_back(maxC[q]); R.maxSt.push_back(maxSt[q]);
        R.reactFrac.push_back(area[q] > 0 ? react[q] / area[q] : 0); R.expected.push_back(expct[q]); R.events.push_back(nEv[q]); R.maxMw.push_back(maxMw[q]);
        if (slipOn_) {   // C26
            double ms = 0, mo = 0;
            const double G = out_.shearModulus, M2 = out_.mPerUnit * out_.mPerUnit;
            for (size_t i = 0; i < nf; ++i) if (size_t(out_.faces[i].fault) == q) { ms = std::max(ms, slip_[i]); mo += G * out_.faces[i].area * M2 * slip_[i]; }
            R.maxSlip.push_back(ms); R.moment.push_back(mo);
        }
    }
    out_.steps.push_back(std::move(S));
}

std::map<int, double> Monitor::connMultipliers() const {
    std::map<int, double> m;
    if (!active_ || !twoWay_) return m;
    for (size_t i = 0; i < out_.faces.size(); ++i) {
        if (mult_[i] == 1.0) continue;
        for (int cn : out_.faces[i].conns) { auto it = m.find(cn); if (it == m.end()) m[cn] = mult_[i]; else it->second = std::max(it->second, mult_[i]); }
    }
    return m;
}

void Monitor::applyCoupling(Simulator& sim) {
    if (!active_ || !twoWay_) return;
    const std::map<int, double> m = connMultipliers();
    if (m.empty()) return;
    std::vector<int> cs; std::vector<double> ms;
    for (const auto& e : m) { cs.push_back(e.first); ms.push_back(e.second); }
    sim.setConnectionMultipliers(cs, ms);
    ++out_.twoWayUpdates;
}

// C47: per-face fault history as CSV (CLI --fault-csv / env RESSIM_FAULT_CSV): one row per report step and face
bool writeFaultCsv(const Output& o, const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) return false;
    std::fprintf(f, "t,face,fault,k,x,y,depth,area,aSig,cff,dcfs,rate,cum,sne,mult\n");   // C49: + sigma'_n, two-way T multiplier
    for (const auto& S : o.steps)
        for (size_t i = 0; i < o.faces.size(); ++i) {
            const Face& F = o.faces[i];
            auto at = [&](const std::vector<float>& v) { return i < v.size() ? double(v[i]) : 0.0; };
            std::fprintf(f, "%.10g,%zu,%s,%d,%.6g,%.6g,%.6g,%.8g,%.8g,%.10g,%.10g,%.10g,%.10g,%.8g,%.8g\n", S.t, i,
                         F.fault >= 0 && size_t(F.fault) < o.names.size() ? o.names[size_t(F.fault)].c_str() : "?", F.k + 1,
                         double(F.c[0]), double(F.c[1]), double(F.c[2]), F.area, F.aSig, at(S.cff), at(S.dcfs), at(S.rate), at(S.cum), at(S.sne), S.mult.empty() ? 1.0 : at(S.mult));
        }
    std::fclose(f);
    return true;
}

void Monitor::finish(std::atomic<bool>* cancel, const std::function<void(const std::string&)>& log) {
    if (active_ && !out_.steps.empty()) {
        if (const char* p = std::getenv("RESSIM_FAULT_CSV"); p && *p) {
            if (writeFaultCsv(out_, p)) { if (log) log(std::string("wrote fault history ") + p); }
            else if (log) log(std::string("NOTE: cannot write ") + p);
        }
    }
    if (!active_ || out_.steps.size() < 2) return;
    const double Pa = out_.paPerUnit, M = out_.mPerUnit, day = 86400.0;
    for (size_t fi = 0; fi < out_.names.size(); ++fi) {
        if (rsfLayer_[fi] <= 0) continue;
        if (cancel && cancel->load()) break;
        RsfInfo info; info.fault = int(fi); info.layer = rsfLayer_[fi];
        std::vector<int> idx;
        for (size_t i = 0; i < out_.faces.size(); ++i) if (out_.faces[i].fault == int(fi) && out_.faces[i].k == rsfLayer_[fi] - 1) idx.push_back(int(i));
        info.faces = int(idx.size());
        if (idx.size() < 3) { info.error = "fewer than 3 faces in layer " + std::to_string(rsfLayer_[fi]); out_.rsf.push_back(info); continue; }
        // trace direction = principal axis of the face centres (map view)
        double mx = 0, my = 0;
        for (int i : idx) { mx += out_.faces[size_t(i)].c[0]; my += out_.faces[size_t(i)].c[1]; }
        mx /= double(idx.size()); my /= double(idx.size());
        double sxx = 0, sxy = 0, syy = 0;
        for (int i : idx) { const double dx = out_.faces[size_t(i)].c[0] - mx, dy = out_.faces[size_t(i)].c[1] - my; sxx += dx * dx; sxy += dx * dy; syy += dy * dy; }
        const double th = 0.5 * std::atan2(2 * sxy, sxx - syy), ux = std::cos(th), uy = std::sin(th);
        std::vector<std::pair<double, int>> sp;
        for (int i : idx) sp.push_back({(out_.faces[size_t(i)].c[0] - mx) * ux + (out_.faces[size_t(i)].c[1] - my) * uy, i});
        std::sort(sp.begin(), sp.end());
        double halfW = 0;   // half the along-strike face width
        {
            const Face& f = out_.faces[size_t(sp.front().second)];
            for (int v = 0; v < 4; ++v) halfW = std::max(halfW, std::fabs((f.quad[size_t(v)][0] - f.c[0]) * ux + (f.quad[size_t(v)][1] - f.c[1]) * uy));
        }
        const double s0 = sp.front().first - halfW, s1 = sp.back().first + halfW, L = (s1 - s0) * M;
        rsf::Input in;
        in.x0 = 0; in.y0 = 0; in.x1 = L; in.y1 = 0;
        in.n = std::clamp(int(sp.size()) * 8, 64, 256);
        in.aVW = rsfA_[fi]; in.aVS = rsfA_[fi] + 0.010; in.b = rsfB_[fi]; in.Dc = rsfDc_[fi];
        in.vwFrom = 0.05; in.vwTo = 0.95;
        in.Vpl = 0; in.vInit = 1e-12;
        const int ne = in.n;
        // per element: interpolation of the face values along the trace
        auto interp = [&](double sm, const std::function<double(int)>& val) {
            const double sCase = s0 + sm / M;
            if (sCase <= sp.front().first) return val(sp.front().second);
            if (sCase >= sp.back().first) return val(sp.back().second);
            size_t q = 1; while (q < sp.size() && sp[q].first < sCase) ++q;
            const double w = (sCase - sp[q - 1].first) / std::max(1e-12, sp[q].first - sp[q - 1].first);
            return (1 - w) * val(sp[q - 1].second) + w * val(sp[q].second);
        };
        in.sigmaN = 0; in.p0 = 0;
        in.sigmaNProf.resize(size_t(ne)); in.tau0Prof.resize(size_t(ne));
        const double dxm = L / ne;
        const double fss = in.f0 + (in.aVW - in.b) * std::log(in.vInit / in.V0);
        for (int e = 0; e < ne; ++e) {
            const double sm = (e + 0.5) * dxm;
            const double sne0 = std::max(0.1e6, interp(sm, [&](int i) { return (out_.faces[size_t(i)].sn0 - out_.faces[size_t(i)].p0) * Pa; }));
            const double t0 = interp(sm, [&](int i) { return out_.faces[size_t(i)].tau0 * Pa; });
            in.sigmaNProf[size_t(e)] = sne0;
            in.tau0Prof[size_t(e)] = std::min(t0, 0.98 * fss * sne0);   // a pre-stress above the strength would fail at t = 0
        }
        // forcing table: effective-normal-stress decrease (as dp) and shear change from the simulated report steps
        for (size_t s = 0; s < out_.steps.size(); ++s) {
            const StepOut& S = out_.steps[s];
            in.forceT.push_back(S.t * day);
            std::vector<double> dp(static_cast<size_t>(ne)), dt(static_cast<size_t>(ne));
            for (int e = 0; e < ne; ++e) {
                const double sm = (e + 0.5) * dxm;
                const double dsne = interp(sm, [&](int i) { return (double(S.sne[size_t(i)]) - (out_.faces[size_t(i)].sn0 - out_.faces[size_t(i)].p0)) * Pa; });
                const double dtau = interp(sm, [&](int i) { return (double(S.dcfs[size_t(i)]) + mu_[fi] * (double(S.sne[size_t(i)]) - (out_.faces[size_t(i)].sn0 - out_.faces[size_t(i)].p0))) * Pa; });
                dp[size_t(e)] = -dsne; dt[size_t(e)] = dtau;
            }
            in.forceP.push_back(dp); in.forceTau.push_back(dt);
        }
        in.tEnd = out_.steps.back().t * day;
        in.dtMax = std::max(60.0, in.tEnd / 200);
        in.maxSteps = 600000;
        rsf::Output ro;
        const auto t0 = std::chrono::steady_clock::now();
        bool ok = false;
        try { ok = rsf::run(in, ro, cancel); } catch (const std::exception& e) { ro.error = e.what(); }
        info.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        info.hStar = ro.hStar; info.Lb = ro.Lb; info.dx = ro.dx; info.events = int(ro.events.size());
        if (!ok) info.error = ro.error.empty() ? "rate-and-state run failed" : ro.error;
        // events -> catalogue (location at the hypocentre on the trace, depth of the layer)
        double zMean = 0; for (int i : idx) zMean += out_.faces[size_t(i)].c[2]; zMean /= double(idx.size());
        // rsf.h assumes a square rupture (out-of-plane width = rupture length); the fault of the model is only H high
        double zLo = 1e300, zHi = -1e300;
        for (const Face& f : out_.faces) if (f.fault == int(fi)) for (const auto& q : f.quad) { zLo = std::min(zLo, double(q[2])); zHi = std::max(zHi, double(q[2])); }
        const double H = std::max(0.0, zHi - zLo) * M;
        for (const rsf::Event& E : ro.events) {
            Event ev; ev.kind = 1; ev.fault = int(fi); ev.t = E.tStart / day;
            const double M0 = E.moment * (E.length > 0 && H > 0 ? std::min(1.0, H / E.length) : 1.0);
            ev.mw = float(M0 > 0 ? 2.0 / 3.0 * (std::log10(M0) - 9.1) : E.Mw);
            const double sc = s0 + E.sHypo / M;
            ev.x = float(mx + sc * ux); ev.y = float(my + sc * uy); ev.z = float(zMean);
            out_.events.push_back(ev);
        }
        if (log) {
            char b[300];
            std::snprintf(b, sizeof b, "    fault %s, layer %d: rate-and-state trace %.0f m, %d elements, h* %.0f m, Lb/dx %.1f -> %d events (%.1f s)%s%s",
                          out_.names[fi].c_str(), info.layer, L, ne, ro.hStar, ro.dx > 0 ? ro.Lb / ro.dx : 0.0, info.events, info.seconds,
                          info.error.empty() ? "" : " - ", info.error.c_str());
            log(b);
        }
        out_.rsf.push_back(info);
    }
    std::sort(out_.events.begin(), out_.events.end(), [](const Event& a, const Event& b) { return a.t < b.t; });
}

}  // namespace faultmech
}  // namespace ressim
