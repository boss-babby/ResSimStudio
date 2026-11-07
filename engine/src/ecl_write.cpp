// C14: ECLIPSE binary result writer (EGRID / INIT / SMSPEC / UNSMRY / UNRST), see ecl.h
#include "ressim/ecl.h"
#include "ressim/grid.h"
#include "ressim/simulator.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>

namespace ressim {
namespace ecl {
namespace {

// Fortran unformatted sequential records, big-endian; data split in blocks of 1000 numbers / 105 strings
class Out {
public:
    explicit Out(const std::string& path) : f_(std::filesystem::u8path(path), std::ios::binary | std::ios::trunc) {
        if (!f_) throw std::runtime_error("cannot write " + path);
    }
    ~Out() = default;
    bool ok() const { return bool(f_); }
    void inte(const char* name, const std::vector<int>& v) { head(name, v.size(), "INTE"); blocks(v, 1000, [](const int& x, unsigned char* p) { be32(uint32_t(x), p); }, 4); }
    void real(const char* name, const std::vector<float>& v) { head(name, v.size(), "REAL"); blocks(v, 1000, [](const float& x, unsigned char* p) { uint32_t u; std::memcpy(&u, &x, 4); be32(u, p); }, 4); }
    void doub(const char* name, const std::vector<double>& v) {
        head(name, v.size(), "DOUB");
        blocks(v, 1000, [](const double& x, unsigned char* p) { uint64_t u; std::memcpy(&u, &x, 8); be32(uint32_t(u >> 32), p); be32(uint32_t(u), p + 4); }, 8);
    }
    void logi(const char* name, const std::vector<int>& v) { head(name, v.size(), "LOGI"); blocks(v, 1000, [](const int& x, unsigned char* p) { be32(x ? 0xFFFFFFFFu : 0u, p); }, 4); }
    void chars(const char* name, const std::vector<std::string>& v) {
        head(name, v.size(), "CHAR");
        blocks(v, 105, [](const std::string& s, unsigned char* p) { std::memset(p, ' ', 8); std::memcpy(p, s.data(), std::min<size_t>(8, s.size())); }, 8);
    }
    void mess(const char* name) { head(name, 0, "MESS"); }
private:
    static void be32(uint32_t x, unsigned char* p) { p[0] = uint8_t(x >> 24); p[1] = uint8_t(x >> 16); p[2] = uint8_t(x >> 8); p[3] = uint8_t(x); }
    void rec(const unsigned char* d, uint32_t n) {
        unsigned char m[4]; be32(n, m);
        f_.write(reinterpret_cast<char*>(m), 4); f_.write(reinterpret_cast<const char*>(d), n); f_.write(reinterpret_cast<char*>(m), 4);
    }
    void head(const char* name, size_t n, const char* type) {
        unsigned char h[16]; std::memset(h, ' ', 8);
        std::memcpy(h, name, std::min<size_t>(8, std::strlen(name)));
        be32(uint32_t(n), h + 8); std::memcpy(h + 12, type, 4);
        rec(h, 16);
    }
    template <class T, class F> void blocks(const std::vector<T>& v, size_t blk, F put, size_t sz) {
        std::vector<unsigned char> b;
        for (size_t i0 = 0; i0 < v.size(); i0 += blk) {
            const size_t n = std::min(blk, v.size() - i0);
            b.assign(n * sz, 0);
            for (size_t i = 0; i < n; ++i) put(v[i0 + i], b.data() + i * sz);
            rec(b.data(), uint32_t(b.size()));
        }
    }
    std::ofstream f_;
};

// civil date from (start date + days) (proleptic Gregorian, H. Hinnant's algorithms)
long daysFromCivil(int y, int m, int d) {
    y -= m <= 2; const long era = (y >= 0 ? y : y - 399) / 400; const long yoe = y - era * 400;
    const long doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1; const long doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return era * 146097 + doe - 719468;
}
void civilFromDays(long z, int& y, int& m, int& d) {
    z += 719468; const long era = (z >= 0 ? z : z - 146096) / 146097; const long doe = z - era * 146097;
    const long yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365; const long doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    const long mp = (5 * doy + 2) / 153; d = int(doy - (153 * mp + 2) / 5 + 1); m = int(mp < 10 ? mp + 3 : mp - 9); y = int(yoe + era * 400 + (m <= 2));
}

std::string kwName(std::string s) {
    for (char& ch : s) ch = (ch >= 'a' && ch <= 'z') ? char(ch - 32) : (std::isalnum(static_cast<unsigned char>(ch)) || ch == '_') ? ch : '_';
    if (s.size() > 8) s = s.substr(0, 8);
    return s;
}

} // namespace

bool writeResultCase(const std::string& base, const Case& c, const Grid& G, const std::vector<SummaryPoint>& summary,
                     const std::vector<Snapshot>& snaps, const std::vector<std::string>& wellNames, std::string* err, WriteReport* rep) {
    WriteReport localRep; WriteReport& R = rep ? *rep : localRep;
    try {
        const int nx = G.nx, ny = G.ny, nz = G.nz, n = nx * ny * nz;
        if (n <= 0) throw std::runtime_error("empty grid");
        // active ordering of the main grid (ECLIPSE: increasing global index); ResSim active index of each
        std::vector<int> act;                 // ECLIPSE active -> ResSim active
        std::vector<int> actnum(size_t(n), 0);
        for (int g = 0; g < n; ++g) {
            const int a = g < int(G.activeOf.size()) ? G.activeOf[size_t(g)] : -1;
            if (a >= 0) { actnum[size_t(g)] = 1; act.push_back(a); }
        }
        const int na = int(act.size());
        if (G.numActive() > na) R.warnings.push_back("local grid refinement / fracture cells are not written (" + std::to_string(G.numActive() - na) + " cells)");
        const bool field = c.units == UnitSystem::Field;
        const int unit = field ? 2 : c.units == UnitSystem::Metric ? 1 : 3;

        // ---- EGRID
        std::vector<float> coord, zcorn;
        const GridInput& gi = c.grid;
        if (gi.cornerPoint() && gi.nx == nx && gi.ny == ny && gi.nz == nz && gi.coord.size() == size_t(6 * (nx + 1) * (ny + 1)) &&
            gi.zcorn.size() == size_t(8 * n)) {
            coord.assign(gi.coord.begin(), gi.coord.end());
            zcorn.assign(gi.zcorn.begin(), gi.zcorn.end());
        } else {
            if (G.corners.size() < size_t(n)) throw std::runtime_error("grid has no corner geometry");
            // pillars through the cell corners of the top and bottom layers (vertical when the cells are)
            coord.assign(size_t(6 * (nx + 1) * (ny + 1)), 0.f);
            for (int j = 0; j <= ny; ++j)
                for (int i = 0; i <= nx; ++i) {
                    const int ci = std::min(i, nx - 1), cj = std::min(j, ny - 1), vi = i > ci ? 1 : 0, vj = j > cj ? 1 : 0;
                    const Vec3& t = G.corners[size_t(ci + nx * cj)][size_t(vi + 2 * vj)];
                    const Vec3& b = G.corners[size_t(ci + nx * cj + nx * ny * (nz - 1))][size_t(4 + vi + 2 * vj)];
                    float* q = &coord[size_t(6 * (i + (nx + 1) * j))];
                    q[0] = t.x; q[1] = t.y; q[2] = t.z; q[3] = b.x; q[4] = b.y; q[5] = std::max(b.z, t.z + 1.f);
                }
            zcorn.assign(size_t(8 * n), 0.f);
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i) {
                        const auto& C = G.corners[size_t(i + nx * (j + ny * k))];
                        for (int t = 0; t < 2; ++t)
                            for (int vj = 0; vj < 2; ++vj)
                                for (int vi = 0; vi < 2; ++vi)
                                    zcorn[size_t(((2 * k + t) * 2 * ny + 2 * j + vj) * 2 * nx + 2 * i + vi)] = C[size_t(4 * t + 2 * vj + vi)].z;
                    }
        }
        {
            Out o(base + ".EGRID");
            std::vector<int> fh(100, 0); fh[0] = 3; fh[1] = 2007; fh[5] = G.nzMatrix > 0 ? (G.dualMode == 2 ? 2 : 1) : 0;
            fh[6] = gi.cornerPoint() ? 1 : 2;
            o.inte("FILEHEAD", fh);
            o.chars("MAPUNITS", {field ? "FEET" : "METRES"});
            o.chars("GRIDUNIT", {field ? "FEET" : "METRES", ""});
            std::vector<int> gh(100, 0); gh[0] = 1; gh[1] = nx; gh[2] = ny; gh[3] = nz; gh[24] = 1; gh[25] = 1;
            o.inte("GRIDHEAD", gh);
            o.real("COORD", coord);
            o.real("ZCORN", zcorn);
            o.inte("ACTNUM", actnum);
            o.inte("ENDGRID", {});
            if (!o.ok()) throw std::runtime_error("write error (EGRID)");
            R.files.push_back(base + ".EGRID");
        }

        // ---- headers shared by INIT / UNRST
        const int phases = 1 * (c.phases.oil ? 1 : 0) + 2 * (c.phases.water ? 1 : 0) + 4 * (c.phases.gas ? 1 : 0);
        const long day0 = daysFromCivil(c.startYear, c.startMonth, c.startDay);
        auto intehead = [&](double days) {
            std::vector<int> ih(411, 0);
            ih[1] = 201702; ih[2] = unit; ih[8] = nx; ih[9] = ny; ih[10] = nz; ih[11] = na; ih[14] = phases ? phases : 7;
            ih[16] = 0;   // wells are not written to the restart arrays (summary vectors carry them)
            int y, m, d; civilFromDays(day0 + long(std::floor(days + 1e-9)), y, m, d);
            ih[64] = d; ih[65] = m; ih[66] = y; ih[94] = 100;
            return ih;
        };
        auto doubhead = [](double days) { std::vector<double> dh(229, 0.0); dh[0] = days; return dh; };
        auto perActive = [&](const std::vector<double>& v, double def) {
            std::vector<float> o(static_cast<size_t>(na), static_cast<float>(def));
            if (v.size() == size_t(n)) for (int e = 0; e < na; ++e) o[size_t(e)] = float(v[size_t(G.globalOf[size_t(act[size_t(e)])])]);
            return o;
        };
        auto perActiveI = [&](const std::vector<int>& v) {
            std::vector<int> o(static_cast<size_t>(na), 1);
            if (v.size() == size_t(n)) for (int e = 0; e < na; ++e) o[size_t(e)] = v[size_t(G.globalOf[size_t(act[size_t(e)])])];
            return o;
        };

        // ---- INIT
        {
            Out o(base + ".INIT");
            o.inte("INTEHEAD", intehead(0));
            o.logi("LOGIHEAD", std::vector<int>(121, 0));
            o.doub("DOUBHEAD", doubhead(0));
            std::vector<float> porv(size_t(n), 0.f);
            for (int e = 0; e < na; ++e) {
                const int a = act[size_t(e)];
                if (size_t(a) < G.poreVolRef.size()) porv[size_t(G.globalOf[size_t(a)])] = float(G.poreVolRef[size_t(a)]);
            }
            o.real("PORV", porv);
            auto cellArr = [&](const std::vector<double>& v) {
                std::vector<float> r(size_t(na), 0.f);
                for (int e = 0; e < na; ++e) { const size_t g = size_t(G.globalOf[size_t(act[size_t(e)])]); if (g < v.size()) r[size_t(e)] = float(v[g]); }
                return r;
            };
            o.real("DX", cellArr(G.cellDx)); o.real("DY", cellArr(G.cellDy)); o.real("DZ", cellArr(G.cellDz));
            o.real("PERMX", perActive(gi.permx, 0)); o.real("PERMY", perActive(gi.permy.empty() ? gi.permx : gi.permy, 0));
            o.real("PERMZ", perActive(gi.permz.empty() ? gi.permx : gi.permz, 0));
            o.real("PORO", perActive(gi.poro, 0)); o.real("NTG", perActive(gi.ntg, 1));
            std::vector<float> depth(size_t(na), 0.f), tops(size_t(na), 0.f);
            for (int e = 0; e < na; ++e) {
                const int a = act[size_t(e)], g = G.globalOf[size_t(a)];
                depth[size_t(e)] = size_t(a) < G.depth.size() ? float(G.depth[size_t(a)]) : G.center[size_t(g)].z;
                float zt = 0; for (int v = 0; v < 4; ++v) zt += G.corners[size_t(g)][size_t(v)].z;
                tops[size_t(e)] = zt / 4;
            }
            o.real("DEPTH", depth); o.real("TOPS", tops);
            o.inte("SATNUM", perActiveI(gi.satnum)); o.inte("PVTNUM", perActiveI(gi.pvtnum));
            o.inte("EQLNUM", perActiveI(gi.eqlnum)); o.inte("FIPNUM", perActiveI(gi.fipnum));
            if (!o.ok()) throw std::runtime_error("write error (INIT)");
            R.files.push_back(base + ".INIT");
        }

        // ---- report steps: summary points that end at a snapshot time
        std::vector<double> repT;
        for (const Snapshot& s : snaps) repT.push_back(s.time);

        // ---- SMSPEC / UNSMRY
        if (!summary.empty()) {
            Summary S = summaryFromRun(summary, wellNames, c.units, componentNames(c));
            const size_t nv = S.vars.size();
            Out o(base + ".SMSPEC");
            o.inte("INTEHEAD", {unit, 100});
            o.chars("RESTART", std::vector<std::string>(9, ""));
            o.inte("DIMENS", {int(nv), nx, ny, nz, 0, -1});
            std::vector<std::string> kw(nv), wg(nv), un(nv); std::vector<int> nums(nv, 0);
            for (size_t i = 0; i < nv; ++i) {
                kw[i] = S.vars[i].keyword; wg[i] = S.vars[i].wgname.empty() ? ":+:+:+:+" : S.vars[i].wgname; un[i] = S.vars[i].unit; nums[i] = S.vars[i].num;
                if (S.vars[i].wgname.size() > 8) R.warnings.push_back("well name truncated to 8 characters: " + S.vars[i].wgname);
            }
            o.chars("KEYWORDS", kw); o.chars("WGNAMES", wg); o.inte("NUMS", nums); o.chars("UNITS", un);
            o.inte("STARTDAT", {c.startDay, c.startMonth, c.startYear, 0, 0, 0});
            if (!o.ok()) throw std::runtime_error("write error (SMSPEC)");
            R.files.push_back(base + ".SMSPEC");

            Out d(base + ".UNSMRY");
            const size_t nt = S.time.size();
            size_t r = 0; int seq = 0; bool open = false;
            while (r < repT.size() && repT[r] <= 1e-12) ++r;    // initial snapshot is not a report end
            for (size_t t = 0; t < nt; ++t) {
                if (!open) { d.inte("SEQHDR", {seq++}); open = true; }
                d.inte("MINISTEP", {int(t)});
                std::vector<float> par(nv);
                for (size_t i = 0; i < nv; ++i) par[i] = t < S.values[i].size() ? S.values[i][t] : 0.f;
                d.real("PARAMS", par);
                const double tt = S.time[t];
                if (r < repT.size() && tt >= repT[r] - 1e-6 * std::max(1.0, std::fabs(repT[r]))) {
                    open = false;
                    while (r < repT.size() && repT[r] <= tt + 1e-6 * std::max(1.0, std::fabs(tt))) ++r;
                }
            }
            if (!d.ok()) throw std::runtime_error("write error (UNSMRY)");
            R.files.push_back(base + ".UNSMRY");
        }

        // ---- UNRST
        if (!snaps.empty()) {
            Out o(base + ".UNRST");
            const int naR = G.numActive();
            auto cells = [&](const std::vector<float>& v, std::vector<float>& out) {
                if (v.size() == size_t(naR)) { out.resize(size_t(na)); for (int e = 0; e < na; ++e) out[size_t(e)] = v[size_t(act[size_t(e)])]; return true; }
                if (v.size() == size_t(na)) { out = v; return true; }
                return false;
            };
            std::set<std::string> skipped;
            int seq = 0;
            for (const Snapshot& s : snaps) {
                o.inte("SEQNUM", {seq++});
                o.inte("INTEHEAD", intehead(s.time));
                o.logi("LOGIHEAD", std::vector<int>(121, 0));
                o.doub("DOUBHEAD", doubhead(s.time));
                o.mess("STARTSOL");
                std::set<std::string> used;
                std::vector<float> buf;
                auto put = [&](const std::string& nm, const std::vector<float>& v) {
                    if (v.empty()) return;
                    std::string k = kwName(nm);
                    for (int dup = 2; used.count(k); ++dup) { const std::string sfx = std::to_string(dup); k = kwName(nm).substr(0, 8 - sfx.size()) + sfx; }
                    if (!cells(v, buf)) { skipped.insert(nm); return; }
                    used.insert(k);
                    o.real(k.c_str(), buf);
                };
                put("PRESSURE", s.pressure); put("SWAT", s.sw); put("SGAS", s.sg); put("SOIL", s.so); put("RS", s.rs); put("RV", s.rv);
                put("POLYMER", s.poly); put("PADS", s.polyAds); put("XCO2", s.xco2);
                for (const auto& e : s.extra) put(e.first, e.second);
                o.mess("ENDSOL");
            }
            for (const auto& k : skipped) R.warnings.push_back("restart array " + k + " skipped (not one value per active cell)");
            if (!o.ok()) throw std::runtime_error("write error (UNRST)");
            R.files.push_back(base + ".UNRST");
        }
        return true;
    } catch (const std::exception& e) {
        if (err) *err = e.what();
        return false;
    }
}

} // namespace ecl
} // namespace ressim
