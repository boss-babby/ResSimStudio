// C25: sequential FE geomechanics and report-step fracture growth for the black-oil / compositional engines - see geoseq.h
#include "ressim/geoseq.h"
#include "ressim/case.h"
#include "ressim/grid.h"
#include "ressim/simulator.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <chrono>
#include <set>

namespace ressim {
namespace geoseq {

std::vector<double> gridNodes(const Grid& g, double lenToM) {
    const int nx = g.nx, ny = g.ny, nz = g.nz;
    if (nx <= 0 || ny <= 0 || nz <= 0 || g.corners.size() < size_t(nx) * size_t(ny) * size_t(nz)) return {};
    const size_t nn = size_t(nx + 1) * size_t(ny + 1) * size_t(nz + 1);
    std::vector<double> X(3 * nn, 0.0), w(nn, 0.0);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const auto& C = g.corners[size_t(i + nx * (j + ny * k))];
                for (int v = 0; v < 8; ++v) {
                    const int I = i + (v & 1), J = j + ((v >> 1) & 1), K = k + ((v >> 2) & 1);
                    const size_t q = (size_t(K) * size_t(ny + 1) + size_t(J)) * size_t(nx + 1) + size_t(I);
                    X[3 * q] += double(C[size_t(v)].x); X[3 * q + 1] += double(C[size_t(v)].y); X[3 * q + 2] += double(C[size_t(v)].z);
                    w[q] += 1;
                }
            }
    for (size_t q = 0; q < nn; ++q)
        for (int d = 0; d < 3; ++d) X[3 * q + size_t(d)] = w[q] > 0 ? X[3 * q + size_t(d)] / w[q] * lenToM : 0.0;
    return X;
}

bool tensorBox(const Grid& g, double tol) {
    const int nx = g.nx, ny = g.ny, nz = g.nz;
    const size_t n = size_t(nx) * size_t(ny) * size_t(nz);
    if (g.cellDx.size() < n || g.cellDy.size() < n || g.cellDz.size() < n || g.corners.size() < n) return false;
    auto id = [&](int i, int j, int k) { return size_t(i + nx * (j + ny * k)); };
    auto near = [&](double a, double b) { return std::fabs(a - b) <= tol * (1 + std::fabs(a) + std::fabs(b)); };
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const size_t c = id(i, j, k);
                if (!near(g.cellDx[c], g.cellDx[id(i, 0, 0)]) || !near(g.cellDy[c], g.cellDy[id(0, j, 0)]) || !near(g.cellDz[c], g.cellDz[id(0, 0, k)])) return false;
            }
    // flat top and vertical, axis-aligned cells (corner z of the top layer equal, x of the -x face corners equal)
    const double z0 = g.corners[0][0].z;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const auto& C = g.corners[id(i, j, 0)];
            for (int v = 0; v < 4; ++v) if (!near(C[size_t(v)].z, z0)) return false;
            if (!near(C[0].x, C[4].x) || !near(C[0].y, C[4].y) || !near(C[0].y, C[1].y) || !near(C[0].x, C[2].x)) return false;
        }
    return true;
}

bool Coupler::setup(const geomech::Input& in, const Grid& g, const std::vector<int>& cellGlobal, double lenToM,
                    const std::vector<double>& pPa, const std::vector<double>& phi, std::string& err, bool generalHex) {
    ready_ = false; steps_ = 0;
    in_ = in; in_.enabled = true;
    const int nx = g.nx, ny = g.ny, nz = g.nz;
    nBox_ = nx * ny * nz;
    const size_t nf = cellGlobal.size();
    if (nBox_ <= 0 || nf == 0 || pPa.size() != nf) { err = "geomechanics: empty grid / pressure"; return false; }
    // uniform property vectors only (per-box arrays are accepted when they match the box)
    auto fix = [&](std::vector<double>& v) { if (v.size() > 1 && v.size() != size_t(nBox_)) v.resize(1); };
    fix(in_.E); fix(in_.nu); fix(in_.biot); fix(in_.alphaT);
    box_.assign(nf, -1); p0_ = pPa;
    double pm = 0; int np = 0;
    for (size_t q = 0; q < nf; ++q) {
        const int gq = cellGlobal[q];
        if (gq >= 0 && gq < nBox_ && !g.isFracture(gq)) { box_[q] = gq; pm += pPa[q]; ++np; }
    }
    if (np == 0) { err = "geomechanics: no flow cell lies in the grid box"; return false; }
    pm /= np;
    std::vector<double> pb(size_t(nBox_), pm), fb(size_t(nBox_), 0.2);
    for (size_t q = 0; q < nf; ++q) if (box_[q] >= 0) { pb[size_t(box_[q])] = pPa[q]; if (q < phi.size() && phi[q] > 1e-4) fb[size_t(box_[q])] = phi[q]; }
    phi0b_ = fb;
    std::string e;
    const bool hex = generalHex || !tensorBox(g);
    bool ok;
    if (hex) {
        const std::vector<double> X = gridNodes(g, lenToM);
        ok = S_.setupHex(in_, nx, ny, nz, X, pb, fb, &e);
    } else {
        std::vector<double> dx(static_cast<size_t>(nBox_)), dy(static_cast<size_t>(nBox_)), dz(static_cast<size_t>(nBox_));
        for (int c = 0; c < nBox_; ++c) { dx[size_t(c)] = g.cellDx[size_t(c)] * lenToM; dy[size_t(c)] = g.cellDy[size_t(c)] * lenToM; dz[size_t(c)] = g.cellDz[size_t(c)] * lenToM; }
        const double top = double(g.corners[0][0].z) * lenToM;
        ok = S_.setup(in_, nx, ny, nz, dx, dy, dz, top, pb, fb, &e);
    }
    if (!ok) { err = e.empty() ? "geomechanics setup failed" : e; return false; }
    ready_ = true;
    return true;
}

bool Coupler::step(const std::vector<double>& pPa, const std::vector<double>* dTK) {
    if (!ready_) return false;
    std::vector<double> dp(size_t(nBox_), 0.0), dT(size_t(nBox_), 0.0);
    for (size_t q = 0; q < box_.size() && q < pPa.size(); ++q) {
        const int b = box_[q];
        if (b < 0) continue;
        dp[size_t(b)] = pPa[q] - p0_[q];
        if (dTK && q < dTK->size()) dT[size_t(b)] = (*dTK)[q];
    }
    const bool ok = S_.solve(dp, dT);
    if (ok) { S_.commit(); ++steps_; }
    return ok;
}

std::vector<double> Coupler::pvFactor(const std::vector<double>& pPa, double crPerPa) const {
    std::vector<double> m(box_.size(), 1.0);
    if (!ready_ || steps_ == 0) return m;
    const auto& phi = S_.porosity();
    for (size_t q = 0; q < box_.size(); ++q) {
        const int b = box_[q];
        if (b < 0 || size_t(b) >= phi.size()) continue;
        const double dp = q < pPa.size() ? pPa[q] - p0_[q] : 0.0;
        const double r = phi[size_t(b)] / std::max(1e-6, phi0b_[size_t(b)]) * std::exp(-crPerPa * dp);
        m[q] = std::isfinite(r) ? std::clamp(r, 0.5, 2.0) : 1.0;
    }
    return m;
}

std::vector<double> Coupler::permFactor() const {
    std::vector<double> m(box_.size(), 1.0);
    if (!ready_ || steps_ == 0) return m;
    const auto& k = S_.permMultiplier();
    for (size_t q = 0; q < box_.size(); ++q) {
        const int b = box_[q];
        if (b >= 0 && size_t(b) < k.size() && std::isfinite(k[size_t(b)])) m[q] = k[size_t(b)];
    }
    return m;
}

std::vector<double> Coupler::connFactor(const std::vector<std::pair<int, int>>& conns) const {
    std::vector<double> out(conns.size(), 1.0);
    const std::vector<double> k = permFactor();
    for (size_t e = 0; e < conns.size(); ++e) {
        const int a = conns[e].first, b = conns[e].second;
        if (a < 0 || b < 0 || size_t(a) >= k.size() || size_t(b) >= k.size() || box_[size_t(a)] < 0 || box_[size_t(b)] < 0) continue;
        const double ka = k[size_t(a)], kb = k[size_t(b)];
        out[e] = ka + kb > 0 ? 2 * ka * kb / (ka + kb) : 1.0;
    }
    return out;
}

std::vector<double> Coupler::minStress() const {
    std::vector<double> s(box_.size(), 0.0);
    if (!ready_ || steps_ == 0) return s;
    for (size_t q = 0; q < box_.size(); ++q) if (box_[q] >= 0) s[q] = S_.minTotalStress(box_[q]);
    return s;
}

bool Coupler::totalStress(int cell, double S[6]) const {
    for (int r = 0; r < 6; ++r) S[r] = 0;
    const int b = boxOf(cell);
    if (!ready_ || b < 0 || steps_ == 0) return false;
    S_.totalStress(b, S);
    return true;
}

std::vector<float> Coupler::map(const std::vector<double>& v, double scale) const {
    std::vector<float> a(box_.size(), 0.0f);
    for (size_t q = 0; q < box_.size(); ++q) { const int b = box_[q]; if (b >= 0 && size_t(b) < v.size()) a[q] = float(v[size_t(b)] * scale); }
    return a;
}

double Coupler::maxFailure() const {
    double m = 0; for (double f : S_.failure()) m = std::max(m, f);
    return m;
}

// ------------------------------------------------------------------------------------------------ report-step HF growth
bool HfSequential::setup(const Case& cx, const Grid& gx, std::string& err) {
    active_ = false;
    const GridInput& g = cx.grid;
    const int nH = gx.numHfCells();
    bool any = false;
    for (const auto& F : g.hfracs) if (F.kic > 0) any = true;
    if (!any) { err = "no HFGROW fracture"; return false; }
    if (nH == 0 || int(g.hfPieces.size()) != nH) { err = "hydraulic fractures not prepared"; return false; }
    const Units U = Units::of(cx.units);
    toPa_ = 101325.0 / U.atm; toM_ = cx.units == UnitSystem::Field ? 0.3048 : 1.0;
    nx_ = g.nx; ny_ = g.ny; nz_ = g.nz;
    F_.clear();
    for (const auto& F : g.hfracs) { Frac f; f.kic = F.kic; f.sigma = F.sigma; f.sigmaGrad = F.sigmaGrad; f.zTop = F.zTop; f.zBot = F.zBot; f.closedMult = F.closedMult; f.fe = F.feStress; f.name = F.name; F_.push_back(f); }
    frac_.assign(size_t(nH), 0); act_.assign(size_t(nH), -1); hostAct_.assign(size_t(nH), -1);
    s_.assign(size_t(nH), 0.0); half_.assign(size_t(nH), 0.0); h_.assign(size_t(nH), 0.0); depth_.assign(size_t(nH), 0.0);
    kf_.assign(size_t(nH), 1.0); kn_.assign(size_t(nH), 0.0); open_.assign(size_t(nH), 1);
    std::vector<char> wellHost(size_t(gx.hfBase()), 0);
    for (const auto& w : cx.wells)
        for (const auto& cp : w.completions)
            for (int k = cp.k1; k <= cp.k2; ++k) {
                if (cp.i < 0 || cp.j < 0 || k < 0 || cp.i >= g.nx || cp.j >= g.ny || k >= g.nz) continue;
                const int q = g.index(cp.i, cp.j, k);
                if (q < int(wellHost.size())) wellHost[size_t(q)] = 1;
            }
    for (int k = 0; k < nH; ++k) {
        const auto& P = g.hfPieces[size_t(k)];
        const auto& HC = gx.hfCells[size_t(k)];
        frac_[size_t(k)] = P.frac;
        act_[size_t(k)] = gx.activeOf[size_t(HC.global)];
        hostAct_[size_t(k)] = P.host >= 0 && size_t(P.host) < gx.activeOf.size() ? gx.activeOf[size_t(P.host)] : -1;
        h_[size_t(k)] = HC.h; kf_[size_t(k)] = std::max(HC.kf, 1e-30); kn_[size_t(k)] = std::max(P.kn, 0.0);
        depth_[size_t(k)] = act_[size_t(k)] >= 0 && size_t(act_[size_t(k)]) < gx.depth.size() ? gx.depth[size_t(act_[size_t(k)])] : 0.5 * (P.zt + P.zb);
    }
    for (size_t f = 0; f < g.hfracs.size(); ++f) {
        const auto& F = g.hfracs[f];
        const double L = std::hypot(F.x1 - F.x0, F.y1 - F.y0), tx = (F.x1 - F.x0) / std::max(L, 1e-12), ty = (F.y1 - F.y0) / std::max(L, 1e-12);
        double sInj = 0; int nInj = 0;
        for (int k = 0; k < nH; ++k) {
            const auto& P = g.hfPieces[size_t(k)];
            if (P.frac != int(f)) continue;
            s_[size_t(k)] = (0.5 * (P.xa + P.xb) - F.x0) * tx + (0.5 * (P.ya + P.yb) - F.y0) * ty;
            half_[size_t(k)] = 0.5 * std::hypot(P.xb - P.xa, P.yb - P.ya);
            if (P.host < int(wellHost.size()) && wellHost[size_t(P.host)]) { sInj += s_[size_t(k)]; ++nInj; }
        }
        sInj = nInj ? sInj / nInj : 0.5 * L;
        for (int k = 0; k < nH; ++k) {
            if (g.hfPieces[size_t(k)].frac != int(f)) continue;
            s_[size_t(k)] -= sInj;
            if (F.kic > 0) open_[size_t(k)] = (std::fabs(s_[size_t(k)]) - half_[size_t(k)] <= F.initLen + 1e-9) ? 1 : 0;
        }
    }
    // fracture-fracture connections
    std::vector<int> pieceOf(size_t(gx.numActive()), -1);
    for (int k = 0; k < nH; ++k) if (act_[size_t(k)] >= 0) pieceOf[size_t(act_[size_t(k)])] = k;
    ff_.clear();
    for (const auto& cn : gx.conns) {
        const int ka = pieceOf[size_t(cn.a)], kb = pieceOf[size_t(cn.b)];
        if (ka >= 0 && kb >= 0) ff_.push_back({cn.a, cn.b, ka, kb});
    }
    active_ = true;
    return true;
}

double HfSequential::wingLength(int frac, int side) const {
    double L = 0;
    for (size_t k = 0; k < open_.size(); ++k)
        if (open_[k] && frac_[k] == frac && (side < 0 ? s_[k] <= 0 : s_[k] >= 0)) L = std::max(L, std::fabs(s_[k]) + half_[k]);
    return L;
}

int HfSequential::grow(double tDays, const std::vector<double>& pCase, const std::vector<double>& sMinCase) {
    if (!active_) return 0;
    const double pi = 3.14159265358979323846;
    int opened = 0;
    for (size_t f = 0; f < F_.size(); ++f) {
        const Frac& F = F_[f];
        if (!(F.kic > 0)) continue;
        for (int side : {-1, 1}) {
            double sTip = -1; int tip = -1;
            for (size_t k = 0; k < open_.size(); ++k) {
                if (frac_[k] != int(f) || !open_[k]) continue;
                const double s = side * s_[k];
                if (s >= -half_[k] && s > sTip) { sTip = s; tip = int(k); }
            }
            if (tip < 0) continue;
            double sNext = 1e300;
            for (size_t k = 0; k < open_.size(); ++k)
                if (frac_[k] == int(f) && !open_[k] && side * s_[k] > sTip && side * s_[k] < sNext) sNext = side * s_[k];
            if (sNext > 1e299) continue;
            double pTip = -1e300, zTip = 0, H = 0; size_t kTip = 0;
            for (size_t k = 0; k < open_.size(); ++k) {
                if (frac_[k] != int(f) || !open_[k] || std::fabs(side * s_[k] - sTip) > 1e-6 * (1 + std::fabs(sTip))) continue;
                const int a = act_[k];
                if (a < 0 || size_t(a) >= pCase.size()) continue;
                if (pCase[size_t(a)] > pTip) { pTip = pCase[size_t(a)]; zTip = depth_[k]; kTip = k; }
                H += h_[k];
            }
            if (pTip < -1e299) continue;
            double sigma;
            const int ha = hostAct_[kTip];
            if (F.fe && ha >= 0 && size_t(ha) < sMinCase.size() && sMinCase[size_t(ha)] > 0) sigma = sMinCase[size_t(ha)];
            else {
                const double zMid = 0.5 * (F.zTop + F.zBot) > 0 ? 0.5 * (F.zTop + F.zBot) : zTip;
                sigma = F.sigma + F.sigmaGrad * (zTip - zMid);
            }
            double pnet = pTip - sigma;
            if (F.fe)   // FE closure stress: the largest net pressure of the tip column drives it (see Simulator::growFractures)
                for (size_t k = 0; k < open_.size(); ++k) {
                    if (frac_[k] != int(f) || !open_[k] || std::fabs(side * s_[k] - sTip) > 1e-6 * (1 + std::fabs(sTip))) continue;
                    const int a = act_[k], hk = hostAct_[k];
                    if (a < 0 || size_t(a) >= pCase.size() || hk < 0 || size_t(hk) >= sMinCase.size() || !(sMinCase[size_t(hk)] > 0)) continue;
                    pnet = std::max(pnet, pCase[size_t(a)] - sMinCase[size_t(hk)]);
                }
            const double Lw = sTip + half_[size_t(tip)];
            const double c = std::min(Lw, 0.5 * H) * toM_;
            const double KI = pnet > 0 ? pnet * toPa_ * std::sqrt(pi * c) : 0.0;
            if (KI < F.kic * 1e6) continue;
            for (size_t k = 0; k < open_.size(); ++k)
                if (frac_[k] == int(f) && !open_[k] && std::fabs(side * s_[k] - sNext) <= 1e-6 * (1 + std::fabs(sNext))) open_[k] = 1;
            Event ev; ev.t = tDays; ev.frac = int(f); ev.side = side; ev.length = wingLength(int(f), side); ev.KI = KI; ev.pnet = pnet;
            ev_.push_back(ev);
            ++opened;
        }
    }
    return opened;
}

void HfSequential::connMultipliers(std::map<std::pair<int, int>, double>& m) const {
    if (!active_) return;
    for (const FF& e : ff_) {
        double v = 1.0;
        for (int k : {e.ka, e.kb}) {
            if (open_[size_t(k)]) continue;
            const double mk = std::min(1.0, kn_[size_t(k)] / kf_[size_t(k)]) * F_[size_t(frac_[size_t(k)])].closedMult;
            v = std::min(v, mk);
        }
        const std::pair<int, int> key{std::min(e.a, e.b), std::max(e.a, e.b)};
        auto it = m.find(key);
        if (it == m.end()) m[key] = v; else it->second *= v;
    }
}

void mechMaps(Snapshot& sn, const geomech::Solver& g, const std::function<int(size_t)>& box, size_t n, double pu, bool plastic) {
    auto mapv = [&](const std::vector<double>& v, double sc) {
        std::vector<float> a(n, 0.0f);
        for (size_t q = 0; q < n; ++q) { const int k = box(q); if (k >= 0 && size_t(k) < v.size()) a[q] = float(v[size_t(k)] * sc); }
        return a;
    };
    std::vector<std::pair<std::string, std::vector<float>>> add;
    add.push_back({"UPLIFT", mapv(g.dispZ(), -1.0)});            // m, + up
    add.push_back({"DISPX", mapv(g.dispX(), 1.0)});
    add.push_back({"DISPY", mapv(g.dispY(), 1.0)});
    {
        std::vector<float> m(n, 0.0f);
        const auto &X = g.dispX(), &Y = g.dispY(), &Z = g.dispZ();
        for (size_t q = 0; q < n; ++q) {
            const int k = box(q);
            if (k >= 0 && size_t(k) < X.size() && size_t(k) < Y.size() && size_t(k) < Z.size())
                m[q] = float(std::sqrt(X[size_t(k)] * X[size_t(k)] + Y[size_t(k)] * Y[size_t(k)] + Z[size_t(k)] * Z[size_t(k)]));
        }
        add.push_back({"DISPMAG", std::move(m)});
    }
    add.push_back({"STRAINV", mapv(g.volStrain(), 1.0)});
    add.push_back({"SIGEFFM", mapv(g.sigEffMean(), 1.0 / pu)});   // pressure units
    add.push_back({"SIGEFFV", mapv(g.sigEffV(), 1.0 / pu)});
    add.push_back({"SIGEFFH", mapv(g.sigEffH(), 1.0 / pu)});
    add.push_back({"FAILURE", mapv(g.failure(), 1.0)});
    add.push_back({"POROM", mapv(g.porosity(), 1.0)});
    {   // C25 smallest principal total stress and max shear of the full tensor
        std::vector<float> smin(n, 0.0f), tmax(n, 0.0f);
        const auto& F = g.sigEffFull();
        for (size_t q = 0; q < n; ++q) {
            const int k = box(q);
            if (k < 0 || 6 * size_t(k) + 5 >= F.size()) continue;
            smin[q] = float(g.minTotalStress(k) / pu);
            const double* s = &F[6 * size_t(k)];
            const double m = (s[0] + s[1] + s[2]) / 3, dxx = s[0] - m, dyy = s[1] - m, dzz = s[2] - m;
            tmax[q] = float(std::sqrt(0.5 * (dxx * dxx + dyy * dyy + dzz * dzz) + s[3] * s[3] + s[4] * s[4] + s[5] * s[5]) / pu);   // sqrt(J2)
        }
        add.push_back({"SIGMIN", std::move(smin)});
        add.push_back({"SQRTJ2", std::move(tmax)});
    }
    if (plastic) {
        add.push_back({"EPSPLAS", mapv(g.plasticStrain(), 1.0)});
        add.push_back({"EPSPV", mapv(g.plasticVolStrain(), 1.0)});
    }
    add.push_back({"KDAMAGE", mapv(g.permMultiplier(), 1.0)});   // geomechanics permeability multiplier
    for (auto& a : add) {
        bool rep = false;
        for (auto& e : sn.extra) if (e.first == a.first) { e.second = std::move(a.second); rep = true; break; }
        if (!rep) sn.extra.push_back(std::move(a));
    }
}

PostResult postProcess(const geomech::Input& inIn, const Grid& g, bool fieldUnits, std::vector<Snapshot>& snaps,
                       const std::function<bool(int, int)>& progress) {
    PostResult R;
    const auto t0 = std::chrono::steady_clock::now();
    if (snaps.empty()) { R.error = "no report steps"; return R; }
    const size_t na = size_t(g.numActive());
    if (snaps[0].pressure.size() != na) { R.error = "the snapshots do not match the grid (active cells)"; return R; }
    const double lenToM = fieldUnits ? 0.3048 : 1.0, paU = fieldUnits ? 6894.757 : 1e5;
    std::vector<double> p0(na), phi(na, 0.2);
    for (size_t a = 0; a < na; ++a) {
        p0[a] = double(snaps[0].pressure[a]) * paU;
        const int gq = g.globalOf[a];
        if (gq >= 0 && size_t(gq) < g.bulkVol.size() && g.bulkVol[size_t(gq)] > 0 && a < g.poreVolRef.size())
            phi[a] = std::clamp(g.poreVolRef[a] * (fieldUnits ? 5.614583 : 1.0) / g.bulkVol[size_t(gq)], 0.01, 0.6);   // rb vs ft3 (C35 fix)
    }
    geomech::Input in = inIn;
    in.enabled = true; in.coupling = 0;
    Coupler C;
    std::string err;
    if (!C.setup(in, g, g.globalOf, lenToM, p0, phi, err)) { R.error = err.empty() ? "mechanics setup failed" : err; return R; }
    const std::vector<float>& T0 = snaps[0].field("TEMP");
    const bool thermal = T0.size() == na;
    const double dTscale = fieldUnits ? 5.0 / 9.0 : 1.0;
    const Coupler* cp = &C;
    auto box = [cp](size_t q) { return cp->boxOf(int(q)); };
    std::vector<double> pPa(na), dT(thermal ? na : 0);
    for (size_t s = 0; s < snaps.size(); ++s) {
        if (progress && !progress(int(s), int(snaps.size()))) { R.error = "cancelled"; break; }
        Snapshot& sn = snaps[s];
        if (sn.pressure.size() != na) { R.error = "snapshot " + std::to_string(s) + " has a different cell count"; break; }
        for (size_t a = 0; a < na; ++a) pPa[a] = double(sn.pressure[a]) * paU;
        if (thermal) {
            const std::vector<float>& T = sn.field("TEMP");
            for (size_t a = 0; a < na; ++a) dT[a] = T.size() == na ? (double(T[a]) - double(T0[a])) * dTscale : 0.0;
        }
        if (!C.step(pPa, thermal ? &dT : nullptr)) { R.error = "the mechanics solve did not converge at snapshot " + std::to_string(s); break; }
        mechMaps(sn, C.solver(), box, na, paU, in.plasticity == 1);
        for (float u : sn.field("UPLIFT")) { R.maxUplift = std::max(R.maxUplift, double(u)); R.maxSubsidence = std::max(R.maxSubsidence, -double(u)); }   // reservoir cells
        R.maxFailure = std::max(R.maxFailure, C.maxFailure());
        ++R.steps;
    }
    R.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return R;
}

}  // namespace geoseq
}  // namespace ressim

namespace ressim {
namespace geoseq {
IterStep iterativeReportStep(Simulator& sim, Coupler& gs, std::vector<double>& pvApplied, double paPerUnit, double crPerPa,
                             int maxIt, double tol, const std::function<void()>& onRewind) {
    IterStep r;
    if (sim.finished()) return r;
    const std::string st0 = sim.saveState();
    const Coupler gs0 = gs;
    const std::vector<double> m0 = pvApplied;
    std::vector<double> mCur = m0;
    maxIt = std::max(1, maxIt);
    if (!(tol > 0)) tol = 1e-4;
    for (int it = 1;; ++it) {
        if (!sim.runReportStep()) return r;
        r.ran = true; r.passes = it;
        const auto& P = sim.snapshots().back().pressure;
        std::vector<double> pPa(P.size());
        for (size_t a = 0; a < P.size(); ++a) pPa[a] = double(P[a]) * paPerUnit;
        if (!gs.step(pPa)) { r.feOk = false; return r; }
        const std::vector<double> m = gs.pvFactor(pPa, crPerPa);
        double diff = 0;
        for (size_t a = 0; a < m.size(); ++a) diff = std::max(diff, std::fabs(m[a] - (a < mCur.size() ? mCur[a] : 1.0)));
        r.diff = diff;
        if (diff <= tol || it >= maxIt) {
            sim.setPoreVolumeMultipliers(m); pvApplied = m; r.converged = diff <= tol;
            return r;
        }
        sim.setPoreVolumeMultipliers(m0.empty() ? std::vector<double>(m.size(), 1.0) : m0);
        sim.loadState(st0);
        if (onRewind) onRewind();
        gs = gs0;
        sim.setPoreVolumeMultipliers(m);
        mCur = m;
    }
}
}  // namespace geoseq
}  // namespace ressim
