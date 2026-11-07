// ResSim - C35 chemical EOR in the black-oil engine: surfactant (SURFACT), alkaline (ALKALINE) and foam (FOAM).
//
// Input: the ECLIPSE keywords stay verbatim extra keywords (exact deck round trip, OPM Flow / ECLIPSE run the same deck)
// and are interpreted here:
//   SURFACT (SURFACTW)  SURFVISC SURFST SURFCAPD SURFADS SURFROCK + REGIONS SURFNUM, SCHEDULE WSURFACT
//   ALKALINE            ALKADS ALKROCK ALSURFST ALSURFAD ALPOLADS, SCHEDULE WALKALIN
//   FOAM                FOAMOPTS FOAMMOB FOAMMOBP FOAMMOBS FOAMADS FOAMROCK FOAMDCYW FOAMDCYO + SOLUTION SFOAM, SCHEDULE WFOAM
//
// Transport (each component, ECLIPSE / OPM Flow storage form; c = mass per surface volume of the carrier phase):
//     M = V c + V_rock rho_r Ca(c) m_ads ,  V = PV b S of the carrier phase (water; foam: gas or water)
//     [M^{n+1} - M^n] / dt + sum_out F c_a - sum_in F c_b + sum_prod q c_a - sum_inj q c_inj + d (V c) = 0
//     d = (e^{lambda dt} - 1) / dt   (foam decay, lambda = ln2 / T_w(Sw) + ln2 / T_o(So), integrating-factor form)
// solved fully implicitly (Newton on the adsorption isotherm, ILU(0)-BiCGSTAB) with the converged upwind phase fluxes
// of the flow step; Ca = max(historical maximum, Ca(c)) for irreversible adsorption (ROCK index 2).
//
// Physical effects (lagged inside the time step and iterated with the flow solve until they stop changing - the same
// outer coupling as the polymer model, so the converged step is the fully coupled solution):
//   surfactant   water viscosity  mu_w(p, c) = mu_w(p) SURFVISC(c) / SURFVISC(0)
//                interfacial tension sigma = SURFST(c) * ALSURFST(c_alk)
//                capillary number  Nc = |K grad Phi_w| / sigma (SI), miscibility F = SURFCAPD(log10 Nc) in [0, 1]
//                relative permeability: water-oil end points (SWCR, SOWCR, SOGCR) interpolated between the immiscible
//                SATNUM table and the miscible SURFNUM table, both curves scaled to them and blended with F (ECLIPSE)
//                capillary pressure Pcow * sigma / SURFST(0)
//                adsorption SURFADS(c) * ALSURFAD(c_alk)
//   alkaline     ALSURFST (IFT), ALSURFAD (surfactant adsorption), ALPOLADS (polymer adsorption) multipliers; ALKADS
//   foam         gas mobility  M = 1 - (1 - M_c(c)) (1 - M_p(p_o)) (1 - M_s(u_g)),  M_c = FOAMMOB, M_p = FOAMMOBP,
//                M_s = FOAMMOBS (0 = no weakening of the foam; OPM Flow uses M_c only - identical when the optional
//                tables are absent); decay half-lives FOAMDCYW / FOAMDCYO
#include "ressim/simulator.h"
#include "ressim/linsolve.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <stdexcept>

namespace ressim {
namespace {
std::string chUp(std::string s) { for (char& c : s) c = char(std::toupper(static_cast<unsigned char>(c))); return s; }
bool chMatch(const std::string& pat, const std::string& name) {
    const std::string a = chUp(pat), b = chUp(name);
    if (!a.empty() && a.back() == '*') return b.compare(0, a.size() - 1, a, 0, a.size() - 1) == 0;
    return a == b;
}
// piecewise-linear (x, y) table, constant beyond the ends
double ctab(const std::vector<double>& t, double x, double* d = nullptr) {
    const size_t n = t.size() / 2;
    if (d) *d = 0;
    if (n == 0) return 0;
    if (n == 1 || x <= t[0]) return t[1];
    if (x >= t[2 * (n - 1)]) return t[2 * n - 1];
    size_t i = 0;
    while (i + 2 < n && x > t[2 * (i + 1)]) ++i;
    const double x0 = t[2 * i], x1 = t[2 * i + 2], y0 = t[2 * i + 1], y1 = t[2 * i + 3];
    const double s = x1 != x0 ? (y1 - y0) / (x1 - x0) : 0.0;
    if (d) *d = s;
    return y0 + s * (x - x0);
}
const std::vector<double>& pick(const std::vector<std::vector<double>>& T, int r) {
    static const std::vector<double> none;
    if (T.empty()) return none;
    return T[std::min(size_t(std::max(r, 0)), T.size() - 1)];
}
char chBuf[512];
template <class... A> std::string cfmt(const char* f, A... a) { std::snprintf(chBuf, sizeof chBuf, f, a...); return chBuf; }
}  // namespace

void Simulator::chemSetup() {
    chem_.clear(); chSurf_ = chAlk_ = chFoam_ = -1;
    chWmul_.clear(); chGmul_.clear(); chMisc_.clear(); chPcMul_.clear(); chLogNc_.clear(); chVg_.clear(); polyAdsMul_.clear();
    chMisSat_.clear(); chArea_.clear();
    for (auto* t : {&surfVisc_, &surfSt_, &surfCapd_, &alSurfSt_, &alSurfAd_, &alPolAds_, &foamMob_, &foamMobP_, &foamMobS_}) t->clear();
    bool surf = false, surfW = false, alk = false, foam = false;
    for (const ExtraKeyword& k : case_.extraKeywords) {
        if (k.section != DeckSection::Runspec) continue;
        if (k.name == "SURFACT") surf = true;
        if (k.name == "SURFACTW") surf = surfW = true;
        if (k.name == "ALKALINE") alk = true;
        if (k.name == "FOAM") foam = true;
    }
    if (!surf && !alk && !foam) return;
    if (!W_) throw std::runtime_error("chemical EOR (SURFACT / ALKALINE / FOAM) needs the WATER phase");
    // tables of a (non-SCHEDULE) keyword: one record per table, an empty record repeats the previous table
    auto tables = [&](const char* name, bool* found = nullptr) {
        std::vector<std::vector<double>> out;
        if (found) *found = false;
        for (const ExtraKeyword& k : case_.extraKeywords) {
            if (k.name != name || k.section == DeckSection::Schedule) continue;
            if (found) *found = true;
            out.clear();
            for (const auto& r : verbatimRecords(k.body)) {
                std::vector<double> v;
                for (const auto& s : r) if (!s.empty()) v.push_back(std::atof(s.c_str()));
                if (v.empty()) { if (!out.empty()) out.push_back(out.back()); continue; }
                if (v.size() % 2) throw std::runtime_error(std::string(name) + ": every table row needs 2 values");
                for (size_t i = 2; i < v.size(); i += 2)
                    if (!(v[i] > v[i - 2])) throw std::runtime_error(std::string(name) + ": the first column must increase");
                out.push_back(v);
            }
        }
        return out;
    };
    auto records = [&](const char* name) {
        std::vector<std::vector<std::string>> out;
        for (const ExtraKeyword& k : case_.extraKeywords)
            if (k.name == name && k.section != DeckSection::Schedule) out = verbatimRecords(k.body);
        return out;
    };
    auto rockRecs = [&](const char* name, ChemSp& S) {
        std::vector<std::string> last;
        for (const auto& r : records(name)) {
            const std::vector<std::string>& u = r.empty() ? last : r;
            if (u.empty()) continue;
            S.adsIdx.push_back(verbatimNum(u, 0, 1) >= 1.5 ? 2 : 1);
            // rock density per bulk volume unit (ft3 | m3): SURFROCK / ALKROCK are mass / reservoir volume (lb/rb), FOAMROCK
            // is a density (lb/ft3) - the OPM keyword dimensions
            const bool perRb = std::string(name) != "FOAMROCK";
            S.rho.push_back(verbatimNum(u, 1, 0) * (perRb ? U_.poreVolume : 1.0));
            last = u;
        }
    };
    const int na = grid_.numActive();
    std::string info;
    if (alk) {
        ChemSp S; S.kind = 1; S.name = "ALKALINE"; S.phase = PW;
        S.ads = tables("ALKADS"); rockRecs("ALKROCK", S);
        alSurfSt_ = tables("ALSURFST"); alSurfAd_ = tables("ALSURFAD"); alPolAds_ = tables("ALPOLADS");
        if (!S.ads.empty() && S.rho.empty()) log_("NOTE: ALKADS without ALKROCK - alkaline adsorption needs the rock density, ignored");
        chAlk_ = int(chem_.size()); chem_.push_back(S);
        info += "alkaline" + std::string(S.ads.empty() || S.rho.empty() ? "" : " (adsorption)") +
                (alSurfSt_.empty() ? "" : ", ALSURFST") + (alSurfAd_.empty() ? "" : ", ALSURFAD") + (alPolAds_.empty() ? "" : ", ALPOLADS");
        if (!alPolAds_.empty() && !poly_) log_("NOTE: ALPOLADS without POLYMER - ignored");
    }
    if (surf) {
        ChemSp S; S.kind = 0; S.name = "SURFACT"; S.phase = PW;
        surfVisc_ = tables("SURFVISC"); surfSt_ = tables("SURFST"); surfCapd_ = tables("SURFCAPD");
        S.ads = tables("SURFADS"); rockRecs("SURFROCK", S);
        if (surfVisc_.empty()) throw std::runtime_error("SURFACT: SURFVISC (water viscosity versus surfactant concentration) is required");
        if (surfSt_.empty()) throw std::runtime_error("SURFACT: SURFST (water-oil surface tension versus surfactant concentration) is required");
        for (const auto& t : surfSt_) for (size_t i = 1; i < t.size(); i += 2) if (!(t[i] > 0)) throw std::runtime_error("SURFACT: SURFST surface tensions must be positive");
        for (const auto& t : surfVisc_) for (size_t i = 1; i < t.size(); i += 2) if (!(t[i] > 0)) throw std::runtime_error("SURFACT: SURFVISC viscosities must be positive");
        if (!S.ads.empty() && S.rho.empty()) throw std::runtime_error("SURFACT: SURFADS needs SURFROCK (adsorption index, rock density)");
        chSurf_ = int(chem_.size()); chem_.push_back(S);
        // miscible (low-IFT) saturation tables
        const auto& sn = case_.grid.surfnum;
        int nMis = 0;
        if (!sn.empty() && !surfCapd_.empty()) {
            chMisSat_.assign(size_t(na), -1);
            for (int a = 0; a < na; ++a) {
                const size_t g = size_t(grid_.globalOf[size_t(a)]);
                const int t = g < sn.size() ? sn[g] - 1 : -1;
                if (t >= int(sats_.size())) throw std::runtime_error(cfmt("SURFACT: SURFNUM %d exceeds the number of saturation tables (%d)", t + 1, int(sats_.size())));
                chMisSat_[size_t(a)] = t;
                if (t >= 0) ++nMis;
            }
        } else if (surfCapd_.empty()) log_("NOTE: SURFACT without SURFCAPD - the relative permeabilities stay immiscible (viscosity, IFT -> Pc and adsorption modelled)");
        else log_("NOTE: SURFACT without SURFNUM - the miscible relative permeability curves are not defined, immiscible curves used");
        info += std::string(info.empty() ? "" : "; ") + "surfactant (SURFVISC, SURFST" + (surfCapd_.empty() ? "" : ", SURFCAPD") +
                (nMis ? cfmt(", SURFNUM miscible curves in %d cells", nMis) : std::string()) + (S.ads.empty() ? "" : ", adsorption") + ")";
        if (surfW) log_("NOTE: SURFACTW - wettability alteration (SURFADDW / SURFWNUM) is not modelled; the SURFACT model is used");
    }
    if (foam) {
        ChemSp S; S.kind = 2; S.name = "FOAM"; S.phase = PG;
        std::string model = "TAB";
        for (const auto& r : records("FOAMOPTS")) {
            const std::string ph = chUp(verbatimStr(r, 0, "GAS"));
            S.phase = ph.rfind("WAT", 0) == 0 ? PW : ph.rfind("SOL", 0) == 0 ? 3 : PG;   // 3: the solvent (SOLVENT option)
            model = chUp(verbatimStr(r, 1, "TAB"));
            break;
        }
        if (S.phase == PG && !G_) throw std::runtime_error("FOAM: the foam is transported in the gas phase (FOAMOPTS) but GAS is not active");
        if (S.phase == 3 && !solv_) throw std::runtime_error("FOAM: FOAMOPTS SOLVENT needs the SOLVENT option in RUNSPEC");
        foamMob_ = tables("FOAMMOB"); foamMobP_ = tables("FOAMMOBP"); foamMobS_ = tables("FOAMMOBS");
        S.ads = tables("FOAMADS"); rockRecs("FOAMROCK", S);
        S.dcyW = tables("FOAMDCYW"); S.dcyO = tables("FOAMDCYO");
        for (auto* T : {&S.dcyW, &S.dcyO}) for (const auto& t : *T) for (size_t i = 1; i < t.size(); i += 2)
            if (!(t[i] > 0)) throw std::runtime_error("FOAM: FOAMDCYW / FOAMDCYO half-lives must be positive");
        if (model == "FUNC") {
            if (foamMob_.empty()) throw std::runtime_error("FOAM: FOAMOPTS FUNC (functional STARS-type model) is not modelled by the ResSim engine - give FOAMMOB (TAB model)");
            log_("NOTE: FOAMOPTS FUNC is not modelled - the FOAMMOB table (TAB model) is used");
        }
        if (foamMob_.empty()) throw std::runtime_error("FOAM: FOAMMOB (gas mobility multiplier versus foam concentration) is required");
        if (!S.ads.empty() && S.rho.empty()) throw std::runtime_error("FOAM: FOAMADS needs FOAMROCK (adsorption index, rock density)");
        chFoam_ = int(chem_.size()); chem_.push_back(S);
        info += std::string(info.empty() ? "" : "; ") + "foam in the " + (S.phase == PG ? "gas" : S.phase == 3 ? "solvent" : "water") + " phase (FOAMMOB" +
                (foamMobP_.empty() ? "" : ", FOAMMOBP") + (foamMobS_.empty() ? "" : ", FOAMMOBS") + (S.ads.empty() ? "" : ", adsorption") +
                (S.dcyW.empty() && S.dcyO.empty() ? "" : ", decay") + ")";
    }
    for (ChemSp& S : chem_) {
        S.c.assign(size_t(na), 0.0); S.adsMax.assign(size_t(na), 0.0); S.massOld.assign(size_t(na), 0.0);
        S.adsMul.assign(size_t(na), 1.0); S.vol.assign(size_t(na), 0.0); S.lam.assign(size_t(na), 0.0);
        S.inj.assign(wells_.size(), 0.0);
        S.cumInj = S.cumProd = S.cumDecay = 0;
    }
    if (chFoam_ >= 0 && !case_.grid.sfoam.empty())
        for (int a = 0; a < na; ++a) {
            const size_t g = size_t(grid_.globalOf[size_t(a)]);
            if (g < case_.grid.sfoam.size()) chem_[size_t(chFoam_)].c[size_t(a)] = std::max(0.0, case_.grid.sfoam[g]);
        }
    chWmul_.assign(size_t(na), 1.0); chGmul_.assign(size_t(na), 1.0); chPcMul_.assign(size_t(na), 1.0);
    chMisc_.assign(size_t(na), 0.0); chLogNc_.assign(size_t(na), -20.0); chVg_.assign(size_t(na), 0.0);
    if (chAlk_ >= 0 && poly_ && !alPolAds_.empty()) polyAdsMul_.assign(size_t(na), 1.0);
    // face areas (capillary number, foam shear velocity)
    const size_t nc = grid_.conns.size();
    chArea_.assign(nc, 0.0);
    for (size_t ci = 0; ci < nc; ++ci) {
        const Connection& cn = grid_.conns[ci];
        const int ga = grid_.globalOf[size_t(cn.a)], gb = grid_.globalOf[size_t(cn.b)];
        const char d = (cn.dir == 'X' || cn.dir == 'Y' || cn.dir == 'Z') ? cn.dir : 'Z';
        const bool aLow = ga < gb;
        chArea_[ci] = std::max(1e-12, 0.5 * (grid_.faceArea(ga, d, !aLow) + grid_.faceArea(gb, d, aLow)));
    }
    // scalar Jacobian pattern: self + neighbours
    chJ_ = CSR{}; chJ_.n = na; chJ_.rowPtr.assign(size_t(na) + 1, 0);
    chSelf_.assign(size_t(na), 0); chAdj_.assign(grid_.adjNbr.size(), 0);
    std::vector<int> cols;
    for (int a = 0; a < na; ++a) {
        cols.clear(); cols.push_back(a);
        for (int e = grid_.adjStart[a]; e < grid_.adjStart[a + 1]; ++e) cols.push_back(grid_.adjNbr[size_t(e)]);
        std::sort(cols.begin(), cols.end()); cols.erase(std::unique(cols.begin(), cols.end()), cols.end());
        const int base = int(chJ_.col.size());
        for (int q : cols) chJ_.col.push_back(q);
        chJ_.rowPtr[size_t(a) + 1] = int(chJ_.col.size());
        auto pos = [&](int q) { return base + int(std::lower_bound(cols.begin(), cols.end(), q) - cols.begin()); };
        chSelf_[size_t(a)] = pos(a);
        for (int e = grid_.adjStart[a]; e < grid_.adjStart[a + 1]; ++e) chAdj_[size_t(e)] = pos(grid_.adjNbr[size_t(e)]);
    }
    chJ_.val.assign(chJ_.col.size(), 0.0);
    chJ_.finalizeDiag();
    log_("Chemical EOR: " + info);
}

void Simulator::chemControls(int step) {
    if (chem_.empty()) return;
    for (const ExtraKeyword& k : case_.extraKeywords) {
        if (k.section != DeckSection::Schedule || k.step != step) continue;
        const int s = k.name == "WSURFACT" ? chSurf_ : k.name == "WALKALIN" ? chAlk_ : k.name == "WFOAM" ? chFoam_ : -2;
        if (s == -2) continue;
        if (s < 0) { log_("NOTE: " + k.name + " without the matching RUNSPEC option - ignored"); continue; }
        for (const auto& r : verbatimRecords(k.body)) {
            if (r.size() < 2 || r[0].empty()) continue;
            const double c = std::max(0.0, verbatimNum(r, 1, 0));
            for (size_t w = 0; w < wells_.size(); ++w)
                if (chMatch(r[0], case_.wells[size_t(wells_[w].input)].name)) chem_[size_t(s)].inj[w] = c;
        }
    }
}

double Simulator::chemMass(const ChemSp& S, int a, double c, double* dM) const {
    const size_t u = size_t(a);
    const double pv = grid_.poreVolRef[u];
    const double V = std::max(S.vol[u], 1e-10 * pv);
    double m = V * c, d = V;
    if (!S.ads.empty() && !S.rho.empty()) {
        const int sr = satReg_[u];
        const double rho = S.rho[std::min(size_t(sr), S.rho.size() - 1)];
        const int idx = S.adsIdx[std::min(size_t(sr), S.adsIdx.size() - 1)];
        const double rockVol = std::max(0.0, grid_.bulkVol[size_t(grid_.globalOf[u])] - pv / U_.poreVolume);   // bulk units
        double da = 0;
        double ca = ctab(pick(S.ads, sr), std::max(c, 0.0), &da) * S.adsMul[u];
        da *= S.adsMul[u];
        if (idx == 2 && S.adsMax[u] > ca) { ca = S.adsMax[u]; da = 0; }
        m += rockVol * rho * ca; d += rockVol * rho * da;
    }
    if (dM) *dM = d;
    return m;
}

double Simulator::chemInPlace(size_t i, double* adsorbed) const {
    if (i >= chem_.size()) return 0;
    const ChemSp& S = chem_[i];
    double free = 0, tot = 0;
    for (size_t a = 0; a < S.c.size(); ++a) { free += S.vol[a] * S.c[a]; tot += chemMass(S, int(a), S.c[a], nullptr); }
    if (adsorbed) *adsorbed = std::max(0.0, tot - free);
    return tot;
}

void Simulator::chemBackup() {
    if (chem_.empty()) return;
    chBak_.c.clear();
    for (const ChemSp& S : chem_) chBak_.c.push_back(S.c);
    chBak_.w = chWmul_; chBak_.g = chGmul_; chBak_.m = chMisc_; chBak_.pc = chPcMul_; chBak_.pa = polyAdsMul_;
}

void Simulator::chemRestore() {
    if (chem_.empty() || chBak_.c.size() != chem_.size()) return;
    for (size_t s = 0; s < chem_.size(); ++s) chem_[s].c = chBak_.c[s];
    chWmul_ = chBak_.w; chGmul_ = chBak_.g; chMisc_ = chBak_.m; chPcMul_ = chBak_.pc; polyAdsMul_ = chBak_.pa;
}

// concentrations + current flow state -> coupling multipliers; returns their largest change
double Simulator::chemEffects() {
    if (chem_.empty()) return 0;
    const int na = grid_.numActive();
    const size_t nc = grid_.conns.size();
    double ch = 0;
    auto upd = [&](double& old, double v, bool rel) {
        const double d = rel ? std::fabs(v - old) / std::max(std::fabs(old), 1e-12) : std::fabs(v - old);
        ch = std::max(ch, d); old = v;
    };
    const double* ca = chAlk_ >= 0 ? chem_[size_t(chAlk_)].c.data() : nullptr;
    // alkaline -> adsorption multipliers
    if (ca && chSurf_ >= 0 && !alSurfAd_.empty()) {
        ChemSp& S = chem_[size_t(chSurf_)];
        for (int a = 0; a < na; ++a) S.adsMul[size_t(a)] = std::max(0.0, ctab(pick(alSurfAd_, satReg_[size_t(a)]), ca[a]));
    }
    if (ca && !polyAdsMul_.empty())
        for (int a = 0; a < na; ++a) upd(polyAdsMul_[size_t(a)], std::max(0.0, ctab(pick(alPolAds_, satReg_[size_t(a)]), ca[a])), true);
    // per-cell |K grad Phi_w| (directional averages of the face values) and gas Darcy velocity
    std::vector<double> gs, gn, vs;
    const bool needNc = chSurf_ >= 0 && !surfCapd_.empty() && !chMisSat_.empty();
    const bool needVg = chFoam_ >= 0 && !foamMobS_.empty() && G_;
    if (needNc || needVg) {
        const double g = U_.gravity;
        if (needNc) { gs.assign(size_t(na) * 3, 0.0); gn.assign(size_t(na) * 3, 0.0); }
        if (needVg) vs.assign(size_t(na) * 3, 0.0);
        std::vector<double> vn(needVg ? size_t(na) * 3 : 0, 0.0);
        for (size_t ci = 0; ci < nc; ++ci) {
            const Connection& cn = grid_.conns[ci];
            const int d = cn.dir == 'X' ? 0 : cn.dir == 'Y' ? 1 : cn.dir == 'Z' ? 2 : -1;
            if (d < 0) continue;
            const CellFlux& A = cf_[size_t(cn.a)]; const CellFlux& B = cf_[size_t(cn.b)];
            if (needNc) {
                const double rho = 0.5 * (A.rho[PW].v + B.rho[PW].v);
                const double dphi = A.P[PW].v - B.P[PW].v + rho * g * cn.dDepth;
                const double kg = cn.trans * std::fabs(dphi) / (U_.cDarcy * chArea_[ci]);   // mD * pressure / length
                for (int x : {cn.a, cn.b}) { gs[size_t(x) * 3 + size_t(d)] += kg; gn[size_t(x) * 3 + size_t(d)] += 1; }
            }
            if (needVg) {
                const double rho = 0.5 * (A.rho[PG].v + B.rho[PG].v);
                const double dphi = A.P[PG].v - B.P[PG].v + rho * g * cn.dDepth;
                const CellFlux& Up = dphi >= 0 ? A : B;
                const int up = dphi >= 0 ? cn.a : cn.b;
                const double Fs = Up.mob[PG].v * dphi * cn.trans;                         // surface gas / day
                const double bg = std::max(1e-12, Up.rho[PG].v / (fl(up).rhoSC().gas * fl(up).gasFactor()));
                const double u = std::fabs(Fs) / bg / U_.poreVolume / chArea_[ci];        // length / day (Darcy)
                for (int x : {cn.a, cn.b}) { vs[size_t(x) * 3 + size_t(d)] += u; vn[size_t(x) * 3 + size_t(d)] += 1; }
            }
        }
        if (needVg)
            for (int a = 0; a < na; ++a) {
                double s2 = 0;
                for (int d = 0; d < 3; ++d) { const size_t q = size_t(a) * 3 + size_t(d); if (vn[q] > 0) { const double m = vs[q] / vn[q]; s2 += m * m; } }
                chVg_[size_t(a)] = std::sqrt(s2);
            }
    }
    const bool field = case_.units == UnitSystem::Field;
    const double kgSI = 9.869233e-16 * (field ? 6894.757 / 0.3048 : 1e5);   // mD * (psi/ft | bar/m) -> N/m
    const double stSI = field ? 175.1268 : 1.0;                               // lbf/in | N/m -> N/m
    double maxRel = 0;
    for (int a = 0; a < na; ++a) {
        const size_t u = size_t(a);
        double wm = 1.0, gm = 1.0, pcm = 1.0, F = 0.0;
        if (chSurf_ >= 0) {
            const double cs = chem_[size_t(chSurf_)].c[u];
            const int pr = pvtReg_[u];
            const auto& vt = pick(surfVisc_, pr);
            wm = ctab(vt, 0.0) / std::max(1e-12, ctab(vt, cs));
            const auto& st = pick(surfSt_, pr);
            double sig = ctab(st, cs);
            if (ca && !alSurfSt_.empty()) sig *= std::max(1e-12, ctab(pick(alSurfSt_, pr), ca[a]));
            sig = std::max(sig, 1e-30);
            pcm = sig / std::max(1e-30, ctab(st, 0.0));
            if (needNc && chMisSat_[u] >= 0) {
                double s2 = 0;
                for (int d = 0; d < 3; ++d) { const size_t q = u * 3 + size_t(d); if (gn[q] > 0) { const double m = gs[q] / gn[q]; s2 += m * m; } }
                const double nc_ = std::sqrt(s2) * kgSI / (sig * stSI);
                chLogNc_[u] = std::log10(std::max(nc_, 1e-30));
                F = std::clamp(ctab(pick(surfCapd_, satReg_[u]), chLogNc_[u]), 0.0, 1.0);
            }
        }
        if (chFoam_ >= 0) {
            const double cf = chem_[size_t(chFoam_)].c[u];
            const int pr = pvtReg_[u];
            const double mc = std::clamp(ctab(pick(foamMob_, pr), cf), 0.0, 1.0);
            const double mp = foamMobP_.empty() ? 0.0 : std::clamp(ctab(pick(foamMobP_, pr), cf_[u].P[O_ ? PO : PW].v), 0.0, 1.0);
            const double ms = foamMobS_.empty() ? 0.0 : std::clamp(ctab(pick(foamMobS_, pr), chVg_[u]), 0.0, 1.0);
            gm = 1.0 - (1.0 - mc) * (1.0 - mp) * (1.0 - ms);
        }
        upd(chWmul_[u], wm, true); upd(chGmul_[u], gm, true); upd(chPcMul_[u], pcm, true); upd(chMisc_[u], F, false);
        maxRel = std::max(maxRel, std::fabs(1.0 - pcm));
    }
    chMaxRelPc_ = maxRel;
    return ch;
}

bool Simulator::chemSolve(ChemSp& S, double dt, const std::vector<double>& F, const std::vector<std::vector<std::array<double, 2>>>& wq) {
    const int na = grid_.numActive();
    const size_t nau = size_t(na);
    double scale = S.scale;
    for (double v : S.inj) scale = std::max(scale, v);
    for (double v : S.c) scale = std::max(scale, v);
    if (!(scale > 0)) return true;                       // nothing in place, nothing injected
    std::vector<double> c = S.c, R(nau), rhs(nau), dx(nau), dec(nau);
    for (int a = 0; a < na; ++a) dec[size_t(a)] = S.lam[size_t(a)] > 0 ? std::expm1(S.lam[size_t(a)] * dt) / dt * S.vol[size_t(a)] : 0.0;
    bool ok = false;
    // residual (and Jacobian) at concentrations cc; returns the scaled max residual
    auto resid = [&](const std::vector<double>& cc, std::vector<double>& RR, bool jac) {
        if (jac) std::fill(chJ_.val.begin(), chJ_.val.end(), 0.0);
        for (int a = 0; a < na; ++a) {
            const size_t u = size_t(a);
            double dM = 0; const double M = chemMass(S, a, cc[u], &dM);
            RR[u] = (M - S.massOld[u]) / dt + dec[u] * cc[u];
            if (jac) chJ_.val[size_t(chSelf_[u])] += dM / dt + dec[u];
            for (int e = grid_.adjStart[a]; e < grid_.adjStart[a + 1]; ++e) {
                const size_t ci = size_t(grid_.adjConn[size_t(e)]);
                const double Fo = grid_.conns[ci].a == a ? F[ci] : -F[ci];   // out of a
                const int b = grid_.adjNbr[size_t(e)];
                if (Fo >= 0) { RR[u] += Fo * cc[u]; if (jac) chJ_.val[size_t(chSelf_[u])] += Fo; }
                else { RR[u] += Fo * cc[size_t(b)]; if (jac) chJ_.val[size_t(chAdj_[size_t(e)])] += Fo; }
            }
        }
        for (size_t w = 0; w < wq.size(); ++w)
            for (const auto& e : wq[w]) {
                const size_t u = size_t(e[0]); const double q = e[1];
                if (q > 0) { RR[u] += q * cc[u]; if (jac) chJ_.val[size_t(chSelf_[u])] += q; }
                else RR[u] += q * S.inj[w];
            }
        double rm = 0;
        for (int a = 0; a < na; ++a) {
            const size_t u = size_t(a);
            rm = std::max(rm, std::fabs(RR[u]) * dt / std::max(S.vol[u], 1e-10 * grid_.poreVolRef[u]));
        }
        return rm;
    };
    std::vector<double> cTry(nau), Rtry(nau);
    double rmax = resid(c, R, true);
    for (int it = 0; it < 60; ++it) {
        if (rmax <= 1e-7 * scale) { ok = true; break; }
        for (int a = 0; a < na; ++a) {
            const size_t u = size_t(a);
            const double s = dt / grid_.poreVolRef[u];
            rhs[u] = -R[u] * s;
            for (int q = chJ_.rowPtr[u]; q < chJ_.rowPtr[u + 1]; ++q) chJ_.val[size_t(q)] *= s;
        }
        chIlu_.factor(chJ_);
        std::fill(dx.begin(), dx.end(), 0.0);
        LinResult lr = bicgstab(chJ_, chIlu_, rhs, dx, 1e-10, 500);
        if (!lr.converged) { LinResult lg = gmres(chJ_, chIlu_, rhs, dx, 1e-10, 500, 40); if (!lg.converged && lg.relResidual > 1e-4) break; }
        bool bad = false;
        for (int a = 0; a < na; ++a) {
            double& d = dx[size_t(a)];
            if (!std::isfinite(d)) { bad = true; break; }
            d = std::clamp(d, -0.5 * scale, 0.5 * scale);
        }
        if (bad) break;
        // backtracking line search on the residual (steep / saturating adsorption isotherms make plain Newton cycle)
        double lam = 1.0, rTry = 0, dmax = 0;
        for (int ls = 0; ls < 12; ++ls) {
            dmax = 0;
            for (int a = 0; a < na; ++a) {
                const size_t u = size_t(a);
                cTry[u] = std::max(0.0, c[u] + lam * dx[u]);
                dmax = std::max(dmax, std::fabs(cTry[u] - c[u]));
            }
            rTry = resid(cTry, Rtry, false);
            if (rTry < rmax || rTry <= 1e-7 * scale) break;
            lam *= 0.5;
        }
        c.swap(cTry);
        rmax = resid(c, R, true);
        if (dmax <= 1e-11 * scale) { ok = true; break; }
    }
    if (!ok) return false;
    S.c.swap(c);
    S.decayPending = 0;
    for (int a = 0; a < na; ++a) S.decayPending += dec[size_t(a)] * S.c[size_t(a)] * dt;   // committed when the step is accepted
    return true;
}

double Simulator::chemTransport(double dt) {
    if (chem_.empty() || dt <= 0) return 0;
    const int na = grid_.numActive();
    const size_t nc = grid_.conns.size(), nt = trInfo_.size();
    const double g = U_.gravity;
    const int refP = O_ ? PO : PW;
    bool needGas = false;
    for (const ChemSp& S : chem_) if (S.phase == PG || S.phase == 3 || !S.dcyO.empty() || !S.dcyW.empty()) needGas = true;
    // carrier volumes and decay rates at the current (converged) state
    for (int a = 0; a < na; ++a) {
        const size_t u = size_t(a);
        double Vg = 0, Vs = 0, sw = 0, so = 0;
        if (needGas) {
            const CellFull C = evalCell(a);
            const double bw = C.b[PW].v, swv = C.S[PW].v;
            const double pv = (bw * swv > 1e-14) ? C.acc[PW].v / (bw * swv)
                                                 : fl(a).pvMult(AD3(p_[u])).v * grid_.poreVolRef[u];
            if (G_) Vg = pv * C.b[PG].v * C.S[PG].v;
            if (solv_) Vs = zs_[u] * solventTotalNT(a);   // solvent amount consistent with the solvent transport fluxes
            sw = swv; so = O_ ? C.S[PO].v : 0.0;
        }
        for (ChemSp& S : chem_) {
            S.vol[u] = S.phase == PW ? cf_[u].acc[PW].v : S.phase == 3 ? Vs : Vg;
            double lam = 0;
            if (!S.dcyW.empty()) lam += std::log(2.0) / std::max(1e-12, ctab(pick(S.dcyW, satReg_[u]), sw));
            if (!S.dcyO.empty() && O_) lam += std::log(2.0) / std::max(1e-12, ctab(pick(S.dcyO, satReg_[u]), so));
            S.lam[u] = lam;
        }
    }
    auto phaseFlux = [&](int p, std::vector<double>& Fp) {
        Fp.assign(nc, 0.0);
        for (size_t ci = 0; ci < nc; ++ci) {
            const Connection& cn = grid_.conns[ci];
            const CellFlux& A = cf_[size_t(cn.a)]; const CellFlux& B = cf_[size_t(cn.b)];
            const double rho = 0.5 * (A.rho[p].v + B.rho[p].v);
            double dphi = A.P[p].v - B.P[p].v + rho * g * cn.dDepth;
            if (!thp_.empty() && thp_[ci] > 0) { const double th = thp_[ci]; dphi = dphi > th ? dphi - th : dphi < -th ? dphi + th : 0.0; }
            Fp[ci] = (dphi >= 0 ? A.mob[p].v : B.mob[p].v) * dphi * cn.trans;
        }
    };
    std::vector<double> Fw, Fg;
    phaseFlux(PW, Fw);
    if (needGas && G_) phaseFlux(PG, Fg);
    // converged well connection rates of water / gas (q > 0 production, q < 0 injection; x WEFAC)
    std::vector<std::vector<std::array<double, 2>>> qW(wells_.size()), qG(wells_.size());
    for (size_t w = 0; w < wells_.size(); ++w) {
        const WellRT& W = wells_[w];
        if (!W.open) continue;
        const bool prod = W.type == WellType::Producer;
        for (size_t k = 0; k < W.cell.size(); ++k) {
            if (W.wi[k] <= 0) continue;
            const int a = W.cell[k];
            const CellFlux& C = cf_[size_t(a)];
            double q[3] = {0, 0, 0};
            if (prod) {
                const double dd = C.P[refP].v - W.pConn(k, g);
                if (dd <= 0) continue;
                for (int p = 0; p < 3; ++p) if (eq_[p] >= 0) q[p] = W.wi[k] * C.mob[p].v * dd;
                if (k < W.wmulConn.size()) q[PW] *= W.wmulConn[k];
            } else {
                const int ip = W.type == WellType::GasInjector ? PG : PW;
                const double dd = W.pConn(k, g) - C.P[ip].v;
                if (dd <= 0 || eq_[ip] < 0) continue;
                const CellFull Cf = evalCell(a);
                double lt = 0;
                for (int p = 0; p < 3; ++p)
                    if (eq_[p] >= 0 && Cf.b[p].v > 0) lt += Cf.mob[p].v / Cf.b[p].v * (p == PW && k < W.wmulConn.size() ? W.wmulConn[k] : 1.0);
                q[ip] = -W.wi[k] * lt * Cf.b[ip].v * dd;
            }
            if (k < W.qwConn.size() && W.qwConn[k] != 0.0) q[PW] = W.qwConn[k];   // the well model's own water rate (FWIT / FWPT consistent)
            if (q[PW] != 0) qW[w].push_back({double(a), q[PW] * W.eff});
            if (q[PG] != 0) qG[w].push_back({double(a), q[PG] * W.eff});
        }
        // gas injectors: honour the well model's injected gas rate exactly (same formula, last-assembly state)
        if (W.type == WellType::GasInjector && W.rates.gasInj > 0 && !qG[w].empty()) {
            double sq = 0;
            for (const auto& e : qG[w]) sq -= e[1];
            if (sq > 0) for (auto& e : qG[w]) e[1] *= W.rates.gasInj * W.eff / sq;
        }
    }
    // alkaline first (its multipliers act on the surfactant adsorption), then surfactant, then foam
    std::vector<int> order;
    for (int s : {chAlk_, chSurf_, chFoam_}) if (s >= 0) order.push_back(s);
    for (int s : order) {
        ChemSp& S = chem_[size_t(s)];
        if (s == chSurf_ && chAlk_ >= 0 && !alSurfAd_.empty())
            for (int a = 0; a < na; ++a) S.adsMul[size_t(a)] = std::max(0.0, ctab(pick(alSurfAd_, satReg_[size_t(a)]), chem_[size_t(chAlk_)].c[size_t(a)]));
        const bool gas = S.phase == PG, sol = S.phase == 3;
        const std::vector<double>& Fc = sol ? svF_ : gas ? Fg : Fw;
        const auto& wq = sol ? svWq_ : gas ? qG : qW;
        if (!chemSolve(S, dt, Fc, wq)) { log_("  " + S.name + " transport did not converge"); if (std::getenv("RESSIM_DEBUG")) std::fprintf(stderr, "    [%s transport failed]\n", S.name.c_str()); return -1.0; }
        // well rates (amount / s, flowing-rate convention as the tracers: + produced, - injected)
        for (size_t w = 0; w < wells_.size(); ++w) {
            auto& comp = wells_[w].rates.comp;
            if (comp.size() < nt + chem_.size()) comp.resize(nt + chem_.size(), 0.0);
            double r = 0;
            for (const auto& e : wq[w]) r += e[1] > 0 ? e[1] * S.c[size_t(e[0])] : e[1] * S.inj[w];
            comp[nt + size_t(s)] = r / 86400.0 / (wells_[w].eff > 0 ? wells_[w].eff : 1.0);
        }
    }
    return chemEffects();
}

void Simulator::chemAccept() {
    if (chem_.empty()) return;
    const int na = grid_.numActive();
    for (ChemSp& S : chem_) { S.cumDecay += S.decayPending; S.decayPending = 0; }   // decay of the accepted step
    // carrier volumes at the accepted state (the old-time mass of the next step)
    bool needGas = false;
    for (const ChemSp& S : chem_) if (S.phase == PG || S.phase == 3) needGas = true;
    for (int a = 0; a < na; ++a) {
        const size_t u = size_t(a);
        double Vg = 0, Vs = 0;
        if (needGas && G_) {
            const CellFull C = evalCell(a);
            const double bw = C.b[PW].v, swv = C.S[PW].v;
            const double pv = (bw * swv > 1e-14) ? C.acc[PW].v / (bw * swv) : fl(a).pvMult(AD3(p_[u])).v * grid_.poreVolRef[u];
            Vg = pv * C.b[PG].v * C.S[PG].v;
            if (solv_) Vs = zs_[u] * solventTotalNT(a);   // solvent amount consistent with the solvent transport fluxes
        }
        for (ChemSp& S : chem_) {
            S.vol[u] = S.phase == PW ? cf_[u].acc[PW].v : S.phase == 3 ? Vs : Vg;
            if (!S.ads.empty() && !S.rho.empty()) {
                const int sr = satReg_[u];
                if (S.adsIdx[std::min(size_t(sr), S.adsIdx.size() - 1)] == 2)
                    S.adsMax[u] = std::max(S.adsMax[u], ctab(pick(S.ads, sr), S.c[u]) * S.adsMul[u]);
            }
            S.massOld[u] = chemMass(S, a, S.c[u], nullptr);
        }
    }
}

} // namespace ressim
