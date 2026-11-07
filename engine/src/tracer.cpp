// ResSim - C27 passive tracers in the black-oil engine (CMG IMEX "tracer module" / TOUGH3 EOS1 + EOS7R-style decay)
//
// Input: the ECLIPSE keywords TRACER (PROPS: name, phase WAT / OIL / GAS), WTRACER (SCHEDULE: well, concentration) and
// TBLKF<name> (SOLUTION: initial concentration per cell) stay verbatim extra keywords - the deck round trip is exact and
// OPM Flow runs the same deck. The ResSim extension '-- RESSIM: TRACERPROP name halfLife kd [init]' adds first-order
// decay (half-life, days / time unit of the case), linear equilibrium adsorption (adsorbed amount per unit pore volume
// = kd * c) and a uniform initial concentration.
//
// Model: a tracer is carried by the free phase it is declared in (concentration c = amount per surface volume of that
// phase). After every converged flow step one linear, fully implicit, upwind transport equation per tracer is solved
// on the converged fluxes (sequential transport, as the ECLIPSE / OPM tracer models):
//     [ (V c)^{n+1} + kd PV^{n+1} c^{n+1} - (V c)^n - kd PV^n c^n ] / dt + sum_out F c_a - sum_in F c_b
//         + sum_prod q c_a - sum_inj q c_inj + (e^{lambda dt} - 1) / dt (V + kd PV)^{n+1} c^{n+1} = 0
// The decay term is the integrating-factor form of lambda M: exact for decay alone (any half-life / step ratio, as needed
// for short-lived radionuclide tracers), lambda M + O(dt) otherwise; it only adds to the diagonal (M-matrix kept).
// V = PV b S of the phase (surface volume), F = upwind phase surface flux of the converged state, q = connection rates
// of the converged well state. Upwinding + implicit decay give an M-matrix: c >= 0, no time step limit, monotone.
#include "ressim/simulator.h"
#include "ressim/linsolve.h"
#include <cmath>
#include <algorithm>
#include <cstdio>

namespace ressim {
namespace {
bool trMatch(const std::string& pat, const std::string& name) {
    std::string a = pat, b = name;
    for (char& c : a) c = char(std::toupper(static_cast<unsigned char>(c)));
    for (char& c : b) c = char(std::toupper(static_cast<unsigned char>(c)));
    if (!a.empty() && a.back() == '*') return b.compare(0, a.size() - 1, a, 0, a.size() - 1) == 0;
    return a == b;
}
std::string trUpper(std::string s) { for (char& c : s) c = char(std::toupper(static_cast<unsigned char>(c))); return s; }
}  // namespace

void Simulator::tracerSetup() {
    trInfo_.clear(); trC_.clear(); trVolOld_.clear(); trInj_.clear(); trDecayCum_.clear();
    const int na = grid_.numActive();
    for (const ExtraKeyword& k : case_.extraKeywords) {
        if (k.name != "TRACER") continue;
        for (const auto& r : verbatimRecords(k.body)) {
            if (r.size() < 2 || r[0].empty()) continue;
            TracerInfo T; T.name = trUpper(r[0]);
            const std::string ph = trUpper(r[1]);
            T.phase = ph.rfind("WAT", 0) == 0 ? PW : ph.rfind("OIL", 0) == 0 ? PO : ph.rfind("GAS", 0) == 0 ? PG : -1;
            if (T.phase < 0 || eq_[T.phase] < 0) { log_("NOTE: tracer " + T.name + ": phase " + ph + " is not active / not supported - ignored"); continue; }
            bool dup = false; for (auto& q : trInfo_) if (q.name == T.name) dup = true;
            if (!dup) trInfo_.push_back(T);
        }
    }
    if (case_.salt.on && eq_[PW] >= 0) {   // C29 dissolved NaCl of SALT: internal water tracer (kg / m3 water)
        bool dup = false; for (auto& q : trInfo_) if (q.name == "NACL") dup = true;
        if (!dup) { TracerInfo T; T.name = "NACL"; T.phase = PW; trInfo_.push_back(T); }
    }
    if (trInfo_.empty()) {
        if (!case_.tracerProps.empty()) log_("NOTE: RESSIM TRACERPROP without a TRACER keyword - no tracer simulated");
        return;
    }
    const size_t nt = trInfo_.size();
    trC_.assign(nt * size_t(na), 0.0);
    trInj_.assign(nt, std::vector<double>(wells_.size(), 0.0));
    trDecayCum_.assign(nt, 0.0);
    for (size_t t = 0; t < nt; ++t) {
        TracerInfo& T = trInfo_[t];
        double init = 0;
        for (const auto& P : case_.tracerProps)
            if (P.name == T.name) { T.halfLife = P.halfLife; T.kd = P.kd; init = P.init; T.kPart = P.kPart; T.alphaL = P.alphaL; T.dm = P.dm; }
        if (T.kPart > 0) {   // partner phase of a partitioning tracer
            T.partner = T.phase == PO ? PW : PO;
            if (eq_[T.partner] < 0) { log_("NOTE: tracer " + T.name + ": partition partner phase not active - partitioning ignored"); T.kPart = 0; T.partner = -1; }
        }
        if (case_.salt.on && T.name == "NACL") {
            init = case_.salt.c0;
            for (size_t w = 0; w < wells_.size(); ++w) trInj_[t][w] = std::max(0.0, case_.salt.cInj);
        }
        for (int a = 0; a < na; ++a) trC_[t * size_t(na) + size_t(a)] = init;
        // TBLKF<name> / TBLKS<name>: per-cell initial concentration (global cell order, n*v expanded)
        for (const ExtraKeyword& k : case_.extraKeywords) {
            if (k.name != "TBLKF" + T.name && k.name != "TBLKS" + T.name && k.name != "TBLK" + T.name) continue;
            std::vector<double> v;
            for (const auto& r : verbatimRecords(k.body)) for (const auto& s : r) v.push_back(std::atof(s.c_str()));
            for (int a = 0; a < na; ++a) {
                const size_t g = size_t(grid_.globalOf[size_t(a)]);
                if (g < v.size()) trC_[t * size_t(na) + size_t(a)] = std::max(0.0, v[g]);
            }
        }
        char hl[64], b[256];
        if (T.halfLife > 0) std::snprintf(hl, sizeof hl, "%.4g d", T.halfLife); else std::snprintf(hl, sizeof hl, "none (stable)");
        std::snprintf(b, sizeof b, "Tracer %s in the %s phase: half-life %s, adsorption kd %.4g, initial %.4g", T.name.c_str(),
                      T.phase == PW ? "water" : T.phase == PO ? "oil" : "gas", hl, T.kd, init);
        std::string extra;
        if (T.kPart > 0) { char e[96]; std::snprintf(e, sizeof e, ", partition K %.4g into %s", T.kPart, T.partner == PO ? "oil" : "water"); extra += e; }
        if (T.alphaL > 0 || T.dm > 0) { char e[96]; std::snprintf(e, sizeof e, ", dispersivity %.4g, diffusion %.4g", T.alphaL, T.dm); extra += e; }
        log_(b + extra);
    }
    // C28 decay chains: parent index, solve order parents first (a cycle breaks the chain with a note)
    for (size_t t = 0; t < nt; ++t)
        for (const auto& P : case_.tracerProps)
            if (P.name == trInfo_[t].name && !P.parent.empty()) {
                for (size_t q = 0; q < nt; ++q) if (trInfo_[q].name == P.parent) trInfo_[t].parent = int(q);
                if (trInfo_[t].parent < 0) log_("NOTE: tracer " + trInfo_[t].name + ": decay-chain parent " + P.parent + " is not a declared tracer - ignored");
                else log_("Tracer " + trInfo_[t].name + " is the decay daughter of " + P.parent);
            }
    trOrder_.clear();
    {
        std::vector<char> done(nt, 0);
        for (size_t pass = 0; pass <= nt && trOrder_.size() < nt; ++pass)
            for (size_t t = 0; t < nt; ++t)
                if (!done[t] && (trInfo_[t].parent < 0 || done[size_t(trInfo_[t].parent)])) { done[t] = 1; trOrder_.push_back(int(t)); }
        for (size_t t = 0; t < nt; ++t)
            if (!done[t]) { log_("NOTE: tracer " + trInfo_[t].name + ": decay-chain cycle - parent link ignored"); trInfo_[t].parent = -1; trOrder_.push_back(int(t)); }
    }
    // scalar pattern: self + neighbours (sorted, unique) - the polymer transport layout
    trJ_ = CSR{}; trJ_.n = na; trJ_.rowPtr.assign(size_t(na) + 1, 0);
    trSelf_.assign(size_t(na), 0); trAdj_.assign(grid_.adjNbr.size(), 0);
    std::vector<int> cols;
    for (int a = 0; a < na; ++a) {
        cols.clear(); cols.push_back(a);
        for (int e = grid_.adjStart[a]; e < grid_.adjStart[a + 1]; ++e) cols.push_back(grid_.adjNbr[size_t(e)]);
        std::sort(cols.begin(), cols.end()); cols.erase(std::unique(cols.begin(), cols.end()), cols.end());
        const int base = int(trJ_.col.size());
        for (int q : cols) trJ_.col.push_back(q);
        trJ_.rowPtr[size_t(a) + 1] = int(trJ_.col.size());
        auto pos = [&](int q) { return base + int(std::lower_bound(cols.begin(), cols.end(), q) - cols.begin()); };
        trSelf_[size_t(a)] = pos(a);
        for (int e = grid_.adjStart[a]; e < grid_.adjStart[a + 1]; ++e) trAdj_[size_t(e)] = pos(grid_.adjNbr[size_t(e)]);
    }
    trJ_.val.assign(trJ_.col.size(), 0.0);
    trJ_.finalizeDiag();
    tracerVolumes();
    for (auto& w : wells_) w.rates.comp.assign(nt, 0.0);
}

void Simulator::tracerControls(int step) {
    if (trInfo_.empty()) return;
    for (const ExtraKeyword& k : case_.extraKeywords) {
        if (k.section != DeckSection::Schedule || k.step != step || k.name != "WTRACER") continue;
        for (const auto& r : verbatimRecords(k.body)) {
            if (r.size() < 3) continue;
            const std::string tn = trUpper(r[1]);
            const double c = std::max(0.0, std::atof(r[2].c_str()));
            for (size_t t = 0; t < trInfo_.size(); ++t) {
                if (trInfo_[t].name != tn) continue;
                for (size_t w = 0; w < wells_.size(); ++w)
                    if (trMatch(r[0], case_.wells[size_t(wells_[w].input)].name)) trInj_[t][w] = c;
            }
        }
    }
}

// capacity (V + kd PV) per tracer and cell at the current state
void Simulator::tracerVolumes() {
    const int na = grid_.numActive();
    const size_t nt = trInfo_.size();
    trVolOld_.assign(nt * size_t(na), 0.0);
    bool needFull = false;
    for (const auto& T : trInfo_) if (T.phase != PW || T.kPart > 0 || T.dm > 0) needFull = true;
    if (needFull) { trB_.assign(3 * size_t(na), 0.0); trS_.assign(3 * size_t(na), 0.0); }
    for (int a = 0; a < na; ++a) {
        const double pv = fl(a).pvMult(AD3(p_[size_t(a)])).v * grid_.poreVolRef[size_t(a)];
        double V[3] = {cf_[size_t(a)].acc[PW].v, 0.0, 0.0};
        if (needFull) {
            const CellFull C = evalCell(a);
            for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) { trB_[3 * size_t(a) + size_t(p)] = C.b[p].v; trS_[3 * size_t(a) + size_t(p)] = C.S[p].v; }
            if (eq_[PO] >= 0) V[PO] = pv * C.b[PO].v * C.S[PO].v;
            if (eq_[PG] >= 0) V[PG] = pv * C.b[PG].v * C.S[PG].v;
            if (eq_[PW] >= 0) V[PW] = pv * C.b[PW].v * C.S[PW].v;
        }
        for (size_t t = 0; t < nt; ++t) {
            const TracerInfo& T = trInfo_[t];
            double cap = V[T.phase] + T.kd * pv;
            // partitioning: amount in the partner phase = PV S_partner * K * (c b_carrier)   (K per reservoir volume)
            if (T.kPart > 0) cap += T.kPart * pv * trS_[3 * size_t(a) + size_t(T.partner)] * trB_[3 * size_t(a) + size_t(T.phase)];
            trVolOld_[t * size_t(na) + size_t(a)] = cap;
        }
    }
}

double Simulator::tracerInPlace(size_t t) const {
    const size_t na = size_t(grid_.numActive());
    double m = 0;
    if (t >= trInfo_.size() || trVolOld_.size() < (t + 1) * na) return 0.0;
    for (size_t a = 0; a < na; ++a) m += trVolOld_[t * na + a] * trC_[t * na + a];
    return m;
}

void Simulator::tracerTransport(double dt) {
    if (trInfo_.empty() || dt <= 0) return;
    const int na = grid_.numActive();
    const size_t nt = trInfo_.size();
    const double g = U_.gravity;
    const int refP = O_ ? PO : PW;
    const std::vector<double> capOld = trVolOld_;
    tracerVolumes();                                    // trVolOld_ := capacity at the new state
    const size_t nc = grid_.conns.size();
    const size_t nau = static_cast<size_t>(na);
    std::vector<double> F(nc), rhs(nau), x(nau), cap(nau);
    // converged well connection rates of the tracer phase: q > 0 production (out of the cell), q < 0 injection
    struct CQ { int a; double q[3]; };
    std::vector<std::vector<CQ>> wq(wells_.size());
    for (size_t w = 0; w < wells_.size(); ++w) {
        const WellRT& W = wells_[w];
        const WellRates& R = W.rates;
        if (!W.open || (R.oil == 0 && R.water == 0 && R.gas == 0 && R.waterInj == 0 && R.gasInj == 0)) continue;
        const bool prod = W.type == WellType::Producer;
        for (size_t k = 0; k < W.cell.size(); ++k) {
            if (W.wi[k] <= 0) continue;
            const int a = W.cell[k];
            const CellFlux& C = cf_[size_t(a)];
            CQ e{a, {0, 0, 0}};
            if (prod) {
                const double dd = C.P[refP].v - W.pConn(k, g);
                if (dd <= 0) continue;
                for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) e.q[p] = W.wi[k] * C.mob[p].v * dd;
                if (eq_[PW] >= 0 && k < W.wmulConn.size()) e.q[PW] *= W.wmulConn[k];
            } else {
                const int ip = W.type == WellType::GasInjector ? PG : PW;
                const double dd = W.pConn(k, g) - C.P[ip].v;
                if (dd <= 0 || eq_[ip] < 0) continue;
                const CellFull Cf = evalCell(a);
                double lt = 0;
                for (int p = 0; p < 3; ++p)
                    if (eq_[p] >= 0 && Cf.b[p].v > 0) lt += Cf.mob[p].v / Cf.b[p].v * (p == PW && k < W.wmulConn.size() ? W.wmulConn[k] : 1.0);
                e.q[ip] = -W.wi[k] * lt * Cf.b[ip].v * dd;
            }
            for (int p = 0; p < 3; ++p) e.q[p] *= W.eff;   // C33 WEFAC: the reservoir sees eff x the flowing rate
            wq[w].push_back(e);
        }
    }
    for (auto& w : wells_) { if (w.rates.comp.size() < nt) w.rates.comp.resize(nt, 0.0); std::fill_n(w.rates.comp.begin(), nt, 0.0); }   // C35: chemical entries follow
    // upwind surface flux a -> b of phase p (converged state)
    auto phaseFlux = [&](int p, std::vector<double>& Fp) {
        for (size_t ci = 0; ci < nc; ++ci) {
            const Connection& cn = grid_.conns[ci];
            const CellFlux& A = cf_[size_t(cn.a)]; const CellFlux& B = cf_[size_t(cn.b)];
            const double rho = 0.5 * (A.rho[p].v + B.rho[p].v);
            double dphi = A.P[p].v - B.P[p].v + rho * g * cn.dDepth;
            if (!thp_.empty() && thp_[ci] > 0) { const double th = thp_[ci]; dphi = dphi > th ? dphi - th : dphi < -th ? dphi + th : 0.0; }
            Fp[ci] = (dphi >= 0 ? A.mob[p].v : B.mob[p].v) * dphi * cn.trans;
        }
    };
    std::vector<double> Fpart(nc), dispC(nc), om(nau, 0.0);
    std::vector<std::vector<double>> decTrue(nt), srcCell(nt);   // C28 per-cell decayed amount of the step / ingrowth source (decay chains)
    for (int t : trOrder_) {
        const TracerInfo& T = trInfo_[size_t(t)];
        const int p = T.phase, pp = T.partner;
        const double lam = T.halfLife > 0 ? std::log(2.0) / T.halfLife : 0.0;
        const double decF = std::expm1(lam * dt);   // exponential (integrating-factor) decay: exact without flow
        double* c = &trC_[size_t(t) * size_t(na)];
        phaseFlux(p, F);
        // C28 partitioning: the partner phase carries K c b_carrier / b_partner per surface volume of the partner phase
        if (pp >= 0) {
            phaseFlux(pp, Fpart);
            for (int a = 0; a < na; ++a) {
                const double bp = trB_[3 * size_t(a) + size_t(pp)];
                om[size_t(a)] = bp > 0 ? T.kPart * trB_[3 * size_t(a) + size_t(p)] / bp : 0.0;
            }
        }
        // C28 dispersion + diffusion conductance per connection (amount per concentration difference and day):
        //   alphaL |F| / d  +  (phi S b)_avg Dm A / d      (d = centre distance, A / d from the conduction geometry)
        const bool disp = T.alphaL > 0 || T.dm > 0;
        if (disp) {
            for (size_t ci = 0; ci < nc; ++ci) {
                const Connection& cn = grid_.conns[ci];
                const Vec3& xa = grid_.center[size_t(grid_.globalOf[size_t(cn.a)])];
                const Vec3& xb = grid_.center[size_t(grid_.globalOf[size_t(cn.b)])];
                const double d = std::sqrt((xa.x - xb.x) * (xa.x - xb.x) + (xa.y - xb.y) * (xa.y - xb.y) + (xa.z - xb.z) * (xa.z - xb.z));
                double cdc = 0;
                if (T.alphaL > 0 && d > 0) cdc += T.alphaL * std::fabs(F[ci]) / d;
                if (T.dm > 0 && cn.geo > 0 && !trB_.empty()) {
                    auto phiSb = [&](int a) {
                        const double vol = grid_.bulkVol.empty() ? 0.0 : grid_.bulkVol[size_t(grid_.globalOf[size_t(a)])];
                        const double phi = vol > 0 ? grid_.poreVolRef[size_t(a)] / vol : 0.2;   // FIELD: rb / ft3 - this ratio also converts the ft3/day flux to rb/day
                        return phi * trS_[3 * size_t(a) + size_t(p)] * trB_[3 * size_t(a) + size_t(p)];
                    };
                    // geo = A / (d_a + d_b) in grid length units -> A / d
                    cdc += 0.5 * (phiSb(cn.a) + phiSb(cn.b)) * T.dm * cn.geo;
                }
                dispC[ci] = cdc;
            }
        }
        std::fill(trJ_.val.begin(), trJ_.val.end(), 0.0);
        bool any = false;
        for (int a = 0; a < na; ++a) {
            const size_t ai = size_t(a);
            cap[ai] = trVolOld_[size_t(t) * size_t(na) + ai];
            trJ_.val[size_t(trSelf_[ai])] += cap[ai] * (1.0 + decF) / dt;
            rhs[ai] = capOld[size_t(t) * size_t(na) + ai] * c[ai] / dt;
            if (c[ai] > 0) any = true;
            for (int e = grid_.adjStart[a]; e < grid_.adjStart[a + 1]; ++e) {
                const size_t ci = size_t(grid_.adjConn[size_t(e)]);
                const int b2 = grid_.adjNbr[size_t(e)];
                const double Fo = grid_.conns[ci].a == a ? F[ci] : -F[ci];   // out of a
                if (Fo > 0) trJ_.val[size_t(trSelf_[ai])] += Fo;
                else trJ_.val[size_t(trAdj_[size_t(e)])] += Fo;              // inflow from the neighbour: -|F| c_b
                if (pp >= 0) {
                    const double Fq = grid_.conns[ci].a == a ? Fpart[ci] : -Fpart[ci];
                    if (Fq > 0) trJ_.val[size_t(trSelf_[ai])] += Fq * om[ai];
                    else trJ_.val[size_t(trAdj_[size_t(e)])] += Fq * om[size_t(b2)];
                }
                if (disp && dispC[ci] > 0) { trJ_.val[size_t(trSelf_[ai])] += dispC[ci]; trJ_.val[size_t(trAdj_[size_t(e)])] -= dispC[ci]; }
            }
        }
        // C28 decay chain: the parent's decayed amount D_p of the step feeds the daughter. The ingrown atoms decay during the
        // step too: with the integrating factor the source is  D_p g e^{l2 dt},  g = [l1 (e^{-l1 dt} - e^{-l2 dt}) / (l2 - l1)]
        // / (1 - e^{-l1 dt}) = surviving fraction of the atoms that decayed from the parent within the step (exponential timing)
        // -> exact Bateman for decay alone of a root parent, any step; a stable daughter receives exactly D_p (conservation)
        srcCell[size_t(t)].clear();
        if (T.parent >= 0 && !decTrue[size_t(T.parent)].empty()) {
            const TracerInfo& Pp = trInfo_[size_t(T.parent)];
            const double l1 = Pp.halfLife > 0 ? std::log(2.0) / Pp.halfLife : 0.0;
            if (l1 > 0) {
                const double dl = lam - l1;
                const double num = l1 * std::exp(-l1 * dt) * (std::fabs(dl * dt) < 1e-12 ? dt : -std::expm1(-dl * dt) / dl);
                const double den = -std::expm1(-l1 * dt);
                const double fac = (den > 0 ? num / den : 1.0) * std::exp(lam * dt);
                srcCell[size_t(t)].assign(nau, 0.0);
                for (int a = 0; a < na; ++a) {
                    const double src = decTrue[size_t(T.parent)][size_t(a)] * fac;
                    if (src > 0) { rhs[size_t(a)] += src / dt; any = true; srcCell[size_t(t)][size_t(a)] = src; }
                }
            }
        }
        for (size_t w = 0; w < wells_.size(); ++w)
            for (const CQ& e : wq[w]) {
                if (e.q[p] > 0) trJ_.val[size_t(trSelf_[size_t(e.a)])] += e.q[p];
                else if (e.q[p] < 0 && trInj_[size_t(t)][w] > 0) { rhs[size_t(e.a)] += -e.q[p] * trInj_[size_t(t)][w]; any = true; }
                if (pp >= 0 && e.q[pp] > 0) trJ_.val[size_t(trSelf_[size_t(e.a)])] += e.q[pp] * om[size_t(e.a)];
            }
        if (!any) { std::fill(c, c + na, 0.0); decTrue[size_t(t)].assign(nau, 0.0); continue; }   // nothing in place, nothing injected
        // row scaling by dt / PV (concentration-like rows), then ILU0 + BiCGSTAB (GMRES fallback)
        for (int a = 0; a < na; ++a) {
            const double s = dt / std::max(grid_.poreVolRef[size_t(a)], 1e-300);
            rhs[size_t(a)] *= s;
            for (int q = trJ_.rowPtr[size_t(a)]; q < trJ_.rowPtr[size_t(a) + 1]; ++q) trJ_.val[size_t(q)] *= s;
        }
        for (int a = 0; a < na; ++a) x[size_t(a)] = c[a];
        trIlu_.factor(trJ_);
        LinResult lr = bicgstab(trJ_, trIlu_, rhs, x, 1e-12, 1000);
        if (!lr.converged) {
            LinResult lg = gmres(trJ_, trIlu_, rhs, x, 1e-12, 1000, 60);
            if (!lg.converged && lg.relResidual > 1e-6) log_("WARNING: tracer " + T.name + " transport solve did not converge (residual " + std::to_string(lg.relResidual) + ")");
        }
        double dec = 0, cRef = 0;
        for (size_t w = 0; w < wells_.size(); ++w) cRef = std::max(cRef, trInj_[size_t(t)][w]);
        for (int a = 0; a < na; ++a) cRef = std::max(cRef, std::max(c[a], x[size_t(a)]));
        // the upwind tails decay geometrically: values below 1e-30 of the concentration scale are set to zero, otherwise they
        // reach the denormal range, where every floating-point operation is ~100x slower (measured: 0.2 s per solve)
        const double cTiny = 1e-30 * cRef;
        decTrue[size_t(t)].assign(nau, 0.0);
        for (int a = 0; a < na; ++a) {
            c[a] = std::isfinite(x[size_t(a)]) && x[size_t(a)] > cTiny ? x[size_t(a)] : 0.0;
            // decayed amount of the step: integrating-factor term, minus the part of an ingrowth source that only stands for
            // the decay of atoms born within the step (source - parent decay): the chain conserves the total amount
            double d = decF * cap[size_t(a)] * c[a];
            if (!srcCell[size_t(t)].empty()) d -= srcCell[size_t(t)][size_t(a)] - decTrue[size_t(T.parent)][size_t(a)];
            decTrue[size_t(t)][size_t(a)] = d;
            dec += d;
        }
        trDecayCum_[size_t(t)] += dec;
        // well tracer rates (amount per second, the component-rate convention of the summary): + produced, - injected
        for (size_t w = 0; w < wells_.size(); ++w) {
            double r = 0;
            for (const CQ& e : wq[w]) {
                if (e.q[p] > 0) r += e.q[p] * c[e.a];
                else if (e.q[p] < 0) r += e.q[p] * trInj_[size_t(t)][w];
                if (pp >= 0 && e.q[pp] > 0) r += e.q[pp] * om[size_t(e.a)] * c[e.a];
            }
            wells_[w].rates.comp[size_t(t)] = r / 86400.0 / (wells_[w].eff > 0 ? wells_[w].eff : 1.0);   // flowing-rate convention (summary applies eff)
        }
    }
}

}  // namespace ressim
