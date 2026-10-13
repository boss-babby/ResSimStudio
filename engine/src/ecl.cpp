// ECLIPSE / OPM Flow result file reader. See ecl.h.
#include "ressim/ecl.h"
#include "ressim/grid.h"
#include "ressim/simulator.h"
#include <algorithm>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <functional>
#include <limits>
#include <map>
#include <sstream>

namespace fs = std::filesystem;

namespace ressim {
namespace ecl {

namespace {

std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace((unsigned char)s[a])) ++a;
    while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
}
std::string upper(std::string s) { for (char& c : s) c = char(std::toupper((unsigned char)c)); return s; }

uint32_t be32(const unsigned char* p) { return (uint32_t(p[0]) << 24) | (uint32_t(p[1]) << 16) | (uint32_t(p[2]) << 8) | uint32_t(p[3]); }
uint64_t be64(const unsigned char* p) { return (uint64_t(be32(p)) << 32) | be32(p + 4); }

Type typeOf(const std::string& t, int& strLen) {
    strLen = 8;
    if (t == "INTE") return Type::INTE;
    if (t == "REAL") return Type::REAL;
    if (t == "DOUB") return Type::DOUB;
    if (t == "LOGI") return Type::LOGI;
    if (t == "CHAR") return Type::CHAR;
    if (t == "MESS") return Type::MESS;
    if (t.size() == 4 && t[0] == 'C' && std::isdigit((unsigned char)t[1]) && std::isdigit((unsigned char)t[2]) && std::isdigit((unsigned char)t[3])) {
        strLen = std::atoi(t.c_str() + 1);
        return Type::CHAR;
    }
    return Type::Other;
}
int elemSize(const ArrayInfo& a) {
    switch (a.type) {
    case Type::INTE: case Type::REAL: case Type::LOGI: return 4;
    case Type::DOUB: return 8;
    case Type::CHAR: return a.strLen;
    default: return 0;
    }
}

bool fileExists(const std::string& p) { std::error_code ec; return fs::is_regular_file(fs::u8path(p), ec); }

// header line of a formatted file:   'KEYWORD '          123 'INTE'
bool parseFormattedHeader(const std::string& line, ArrayInfo& a) {
    size_t p = line.find('\'');
    if (p == std::string::npos) return false;
    for (size_t i = 0; i < p; ++i) if (!std::isspace((unsigned char)line[i])) return false;
    const size_t q = line.find('\'', p + 1);
    if (q == std::string::npos || q - p - 1 > 8) return false;
    const std::string name = line.substr(p + 1, q - p - 1);
    size_t r = q + 1;
    while (r < line.size() && std::isspace((unsigned char)line[r])) ++r;
    size_t s = r;
    if (s < line.size() && line[s] == '-') ++s;
    while (s < line.size() && std::isdigit((unsigned char)line[s])) ++s;
    if (s == r) return false;
    const long long cnt = std::atoll(line.substr(r, s - r).c_str());
    while (s < line.size() && std::isspace((unsigned char)line[s])) ++s;
    if (s >= line.size() || line[s] != '\'') return false;
    const size_t e = line.find('\'', s + 1);
    if (e == std::string::npos || e - s - 1 != 4) return false;
    for (size_t i = e + 1; i < line.size(); ++i) if (!std::isspace((unsigned char)line[i])) return false;
    a.name = trim(name);
    a.typeName = line.substr(s + 1, 4);
    a.type = typeOf(a.typeName, a.strLen);
    a.count = cnt;
    return true;
}

std::string extOf(const std::string& p) {
    const size_t sl = p.find_last_of("/\\"), d = p.find_last_of('.');
    if (d == std::string::npos || (sl != std::string::npos && d < sl)) return "";
    return upper(p.substr(d + 1));
}
std::string stripExt(const std::string& p) {
    const size_t sl = p.find_last_of("/\\"), d = p.find_last_of('.');
    if (d == std::string::npos || (sl != std::string::npos && d < sl)) return p;
    return p.substr(0, d);
}
bool isFormattedExt(const std::string& e) {
    if (e.empty()) return false;
    if (e == "FEGRID" || e == "FGRID" || e == "FINIT" || e == "FUNRST" || e == "FSMSPEC" || e == "FUNSMRY" || e == "FRFT") return true;
    if (e.size() == 5 && (e[0] == 'F' || e[0] == 'A') && std::isdigit((unsigned char)e[1])) return true;   // F0001 / A0001
    return false;
}

} // namespace

// ================================================================== File
bool File::open(const std::string& path, std::string* err) {
    path_.clear(); idx_.clear();
    if (in_.is_open()) in_.close();
    in_.clear();
    in_.open(fs::u8path(path), std::ios::binary);
    if (!in_) { if (err) *err = "cannot open " + path; return false; }
    // formatted? first non-blank char is a quote and no binary bytes in the first 64 bytes
    char head[64] = {0};
    in_.read(head, sizeof(head));
    const std::streamsize got = in_.gcount();
    in_.clear(); in_.seekg(0);
    bool fmt = isFormattedExt(extOf(path));
    if (!fmt && got > 4) {
        const unsigned char* u = reinterpret_cast<const unsigned char*>(head);
        if (be32(u) != 16) {   // binary files start with the 16-byte header record
            int i = 0; while (i < got && std::isspace((unsigned char)head[i])) ++i;
            fmt = i < got && head[i] == '\'';
        }
    }
    formatted_ = fmt;
    if (!fmt) {
        unsigned char m[4], h[16];
        for (;;) {
            if (!in_.read(reinterpret_cast<char*>(m), 4)) break;
            const uint32_t len = be32(m);
            if (len != 16) { if (err) *err = path + ": not an ECLIPSE binary file (unexpected record length)"; if (idx_.empty()) return false; break; }
            if (!in_.read(reinterpret_cast<char*>(h), 16)) break;
            in_.read(reinterpret_cast<char*>(m), 4);
            ArrayInfo a;
            a.name = trim(std::string(reinterpret_cast<char*>(h), 8));
            a.count = int32_t(be32(h + 8));
            a.typeName = std::string(reinterpret_cast<char*>(h + 12), 4);
            a.type = typeOf(a.typeName, a.strLen);
            a.offset = int64_t(in_.tellg());
            // skip the data records
            const int es = elemSize(a);
            int64_t remaining = a.type == Type::MESS || a.type == Type::Other ? 0 : a.count * int64_t(es);
            if (a.type == Type::Other && a.count > 0) {
                // unknown type: skip records until the next 16-byte header (best effort)
            }
            bool bad = false;
            while (remaining > 0) {
                if (!in_.read(reinterpret_cast<char*>(m), 4)) { bad = true; break; }
                const uint32_t rl = be32(m);
                in_.seekg(int64_t(rl) + 4, std::ios::cur);
                remaining -= rl;
                if (!in_) { bad = true; break; }
            }
            if (bad) { idx_.push_back(a); break; }   // truncated file (still being written) - keep what is complete
            idx_.push_back(a);
        }
    } else {
        std::string line;
        int64_t pos = 0;
        while (std::getline(in_, line)) {
            const int64_t next = int64_t(in_.tellg());
            ArrayInfo a;
            if (line.size() < 40 && parseFormattedHeader(line, a)) { a.offset = next; idx_.push_back(a); }
            pos = next; (void)pos;
        }
    }
    in_.clear();
    path_ = path;
    if (idx_.empty()) { if (err) *err = path + ": no arrays found"; path_.clear(); return false; }
    return true;
}

int File::find(const std::string& name, int occurrence, int from) const {
    for (int i = std::max(0, from); i < int(idx_.size()); ++i)
        if (idx_[size_t(i)].name == name) { if (occurrence-- == 0) return i; }
    return -1;
}

template <class F> void File::readValues(int i, F&& put) const {
    if (i < 0 || i >= int(idx_.size())) return;
    const ArrayInfo& a = idx_[size_t(i)];
    if (a.count <= 0 || a.type == Type::MESS || a.type == Type::Other) return;
    in_.clear();
    in_.seekg(a.offset);
    if (!formatted_) {
        const int es = elemSize(a);
        int64_t done = 0;
        std::vector<unsigned char> buf;
        unsigned char m[4];
        while (done < a.count) {
            if (!in_.read(reinterpret_cast<char*>(m), 4)) break;
            const uint32_t rl = be32(m);
            buf.resize(rl);
            if (!in_.read(reinterpret_cast<char*>(buf.data()), rl)) break;
            in_.read(reinterpret_cast<char*>(m), 4);
            const int64_t n = rl / es;
            for (int64_t k = 0; k < n && done < a.count; ++k, ++done) {
                const unsigned char* p = buf.data() + k * es;
                switch (a.type) {
                case Type::INTE: put(done, double(int32_t(be32(p))), nullptr); break;
                case Type::LOGI: put(done, be32(p) ? 1.0 : 0.0, nullptr); break;
                case Type::REAL: { uint32_t u = be32(p); float f; std::memcpy(&f, &u, 4); put(done, double(f), nullptr); break; }
                case Type::DOUB: { uint64_t u = be64(p); double d; std::memcpy(&d, &u, 8); put(done, d, nullptr); break; }
                case Type::CHAR: { const std::string s = trim(std::string(reinterpret_cast<const char*>(p), size_t(es))); put(done, 0.0, &s); break; }
                default: break;
                }
            }
        }
    } else {
        int64_t done = 0;
        std::string line;
        while (done < a.count && std::getline(in_, line)) {
            size_t p = 0;
            if (a.type == Type::CHAR) {
                while (done < a.count) {
                    const size_t q = line.find('\'', p);
                    if (q == std::string::npos) break;
                    const size_t r = line.find('\'', q + 1);
                    if (r == std::string::npos) break;
                    const std::string s = trim(line.substr(q + 1, r - q - 1));
                    put(done++, 0.0, &s);
                    p = r + 1;
                }
            } else {
                while (done < a.count) {
                    while (p < line.size() && std::isspace((unsigned char)line[p])) ++p;
                    if (p >= line.size()) break;
                    size_t e = p;
                    while (e < line.size() && !std::isspace((unsigned char)line[e])) ++e;
                    std::string tok = line.substr(p, e - p);
                    p = e;
                    if (a.type == Type::LOGI) { put(done++, (tok[0] == 'T' || tok[0] == 't') ? 1.0 : 0.0, nullptr); continue; }
                    for (char& c : tok) if (c == 'D' || c == 'd') c = 'E';
                    put(done++, std::strtod(tok.c_str(), nullptr), nullptr);
                }
            }
        }
    }
    in_.clear();
}

std::vector<double> File::numeric(int i) const {
    std::vector<double> v;
    if (i < 0 || i >= int(idx_.size())) return v;
    v.assign(size_t(std::max<int64_t>(0, idx_[size_t(i)].count)), 0.0);
    readValues(i, [&](int64_t k, double x, const std::string*) { v[size_t(k)] = x; });
    return v;
}
std::vector<float> File::floats(int i) const {
    std::vector<float> v;
    if (i < 0 || i >= int(idx_.size())) return v;
    v.assign(size_t(std::max<int64_t>(0, idx_[size_t(i)].count)), 0.0f);
    readValues(i, [&](int64_t k, double x, const std::string*) { v[size_t(k)] = float(x); });
    return v;
}
std::vector<int> File::ints(int i) const {
    std::vector<int> v;
    if (i < 0 || i >= int(idx_.size())) return v;
    v.assign(size_t(std::max<int64_t>(0, idx_[size_t(i)].count)), 0);
    readValues(i, [&](int64_t k, double x, const std::string*) { v[size_t(k)] = int(std::lround(x)); });
    return v;
}
std::vector<std::string> File::strings(int i) const {
    std::vector<std::string> v;
    if (i < 0 || i >= int(idx_.size())) return v;
    v.assign(size_t(std::max<int64_t>(0, idx_[size_t(i)].count)), std::string());
    readValues(i, [&](int64_t k, double, const std::string* s) { if (s) v[size_t(k)] = *s; });
    return v;
}

// ================================================================== grid
int EGrid::numActive() const {
    if (actnum.empty()) return nx * ny * nz;
    int n = 0; for (int a : actnum) if (a) ++n;
    return n;
}

bool readEGrid(const std::string& path, EGrid& g, std::string* err) {
    File f;
    if (!f.open(path, err)) return false;
    const auto gh = f.ints("GRIDHEAD");
    if (gh.size() >= 4) { g.nx = gh[1]; g.ny = gh[2]; g.nz = gh[3]; }
    else {
        const auto dim = f.ints("DIMENS");
        if (dim.size() < 3) { if (err) *err = path + ": GRIDHEAD missing (only EGRID corner-point files are supported)"; return false; }
        g.nx = dim[0]; g.ny = dim[1]; g.nz = dim[2];
    }
    g.coord = f.numeric("COORD");
    g.zcorn = f.numeric("ZCORN");
    g.actnum = f.ints("ACTNUM");
    const auto mu = f.strings(f.find("MAPUNITS")); if (!mu.empty()) g.mapUnits = mu[0];
    const auto gu = f.strings(f.find("GRIDUNIT")); if (!gu.empty()) g.gridUnits = gu[0];
    const size_t n = size_t(g.nx) * g.ny * g.nz;
    if (g.coord.size() != size_t(6) * (g.nx + 1) * (g.ny + 1) || g.zcorn.size() != 8 * n) {
        if (err) *err = path + ": COORD/ZCORN sizes do not match GRIDHEAD";
        return false;
    }
    if (!g.actnum.empty() && g.actnum.size() != n) g.actnum.clear();
    return true;
}

std::shared_ptr<Grid> displayGrid(const EGrid& e) {
    auto G = std::make_shared<Grid>();
    Grid& g = *G;
    g.nx = e.nx; g.ny = e.ny; g.nz = e.nz;
    g.cornerPoint = true;
    const int nx = e.nx, ny = e.ny, nz = e.nz;
    const size_t n = size_t(nx) * ny * nz;
    // shift the map origin away so float corner coordinates keep their precision
    double x0 = DBL_MAX, y0 = DBL_MAX;
    for (size_t p = 0; p + 5 < e.coord.size(); p += 6) { x0 = std::min({x0, e.coord[p], e.coord[p + 3]}); y0 = std::min({y0, e.coord[p + 1], e.coord[p + 4]}); }
    if (x0 == DBL_MAX) x0 = y0 = 0;
    auto pillar = [&](int pi, int pj, double z, double& x, double& y) {
        const double* p = &e.coord[6 * (size_t(pj) * (nx + 1) + pi)];
        const double dz = p[5] - p[2];
        if (std::fabs(dz) < 1e-12) { x = p[0]; y = p[1]; return; }
        const double t = (z - p[2]) / dz;
        x = p[0] + t * (p[3] - p[0]); y = p[1] + t * (p[4] - p[1]);
    };
    g.corners.assign(n, {});
    g.center.assign(n, {});
    g.activeOf.assign(n, -1); g.globalOf.clear();
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const size_t q = size_t(i) + size_t(nx) * (size_t(j) + size_t(ny) * k);
                double cx = 0, cy = 0, cz = 0;
                for (int kk = 0; kk < 2; ++kk)
                    for (int jj = 0; jj < 2; ++jj)
                        for (int ii = 0; ii < 2; ++ii) {
                            const size_t zi = (size_t(2 * k + kk) * 2 * ny + (2 * j + jj)) * 2 * nx + 2 * i + ii;
                            const double z = e.zcorn[zi];
                            double x, y; pillar(i + ii, j + jj, z, x, y);
                            g.corners[q][size_t(ii + 2 * jj + 4 * kk)] = {float(x - x0), float(y - y0), float(z)};
                            cx += x - x0; cy += y - y0; cz += z;
                        }
                g.center[q] = {float(cx / 8), float(cy / 8), float(cz / 8)};
                if (e.actnum.empty() || e.actnum[q]) { g.activeOf[q] = int(g.globalOf.size()); g.globalOf.push_back(int(q)); }
            }
    g.bbMin = {FLT_MAX, FLT_MAX, FLT_MAX}; g.bbMax = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
    const int na = g.numActive();
    g.depth.resize(size_t(na)); g.poreVolRef.assign(size_t(na), 0.0); g.dzCell.resize(size_t(na));
    for (int a = 0; a < na; ++a) {
        const int q = g.globalOf[size_t(a)];
        for (auto& v : g.corners[size_t(q)]) {
            g.bbMin.x = std::min(g.bbMin.x, v.x); g.bbMin.y = std::min(g.bbMin.y, v.y); g.bbMin.z = std::min(g.bbMin.z, v.z);
            g.bbMax.x = std::max(g.bbMax.x, v.x); g.bbMax.y = std::max(g.bbMax.y, v.y); g.bbMax.z = std::max(g.bbMax.z, v.z);
        }
        const auto& c = g.corners[size_t(q)];
        g.depth[size_t(a)] = g.center[size_t(q)].z;
        g.dzCell[size_t(a)] = 0.25 * ((c[4].z - c[0].z) + (c[5].z - c[1].z) + (c[6].z - c[2].z) + (c[7].z - c[3].z));
    }
    if (na == 0) { g.bbMin = {0, 0, 0}; g.bbMax = {1, 1, 1}; }
    return G;
}

// ================================================================== summary
namespace {
char categoryOf(const std::string& kw) {
    static const char* timeKw[] = {"TIME", "YEARS", "DAY", "MONTH", "YEAR", "DAYS", "ELAPSED", "TIMESTEP", "TCPU", "TCPUDAY", "TCPUTS", "TELAPLIN"};
    static const char* numKw[] = {"NEWTON", "NLINEARS", "NLINSMIN", "NLINSMAX", "MLINEARS", "MSUMLINS", "MSUMNEWT", "STEPTYPE", "NAIMFRAC", "MAXDPR", "MAXDSO", "MAXDSG", "MAXDSW"};
    for (auto* t : timeKw) if (kw == t) return 'T';
    for (auto* t : numKw) if (kw == t) return 'N';
    if (kw.empty()) return 'M';
    const char c = kw[0];
    if (c == 'F' || c == 'W' || c == 'G' || c == 'R' || c == 'B' || c == 'C' || c == 'A' || c == 'S') return c;
    return 'M';
}
bool validName(const std::string& s) { return !s.empty() && s != ":+:+:+:+"; }
std::string ijkString(int num, int nx, int ny) {
    if (num <= 0 || nx <= 0 || ny <= 0) return std::to_string(num);
    const int g = num - 1;
    return std::to_string(g % nx + 1) + "," + std::to_string((g / nx) % ny + 1) + "," + std::to_string(g / (nx * ny) + 1);
}
} // namespace

int Summary::find(const std::string& key) const {
    for (size_t i = 0; i < vars.size(); ++i) if (vars[i].key == key) return int(i);
    return -1;
}
std::vector<double> Summary::series(int v) const {
    std::vector<double> out;
    if (v < 0 || v >= int(values.size())) return out;
    out.assign(values[size_t(v)].begin(), values[size_t(v)].end());
    return out;
}
double Summary::at(int v, double t) const {
    if (v < 0 || v >= int(values.size()) || time.empty()) return std::numeric_limits<double>::quiet_NaN();
    const auto& y = values[size_t(v)];
    if (y.size() != time.size()) return std::numeric_limits<double>::quiet_NaN();
    if (t <= time.front()) return y.front();
    if (t >= time.back()) return y.back();
    const size_t k = size_t(std::upper_bound(time.begin(), time.end(), t) - time.begin());
    const double t0 = time[k - 1], t1 = time[k];
    const double w = t1 > t0 ? (t - t0) / (t1 - t0) : 1.0;
    return y[k - 1] + w * (y[k] - y[k - 1]);
}

bool readSummary(const std::string& smspec, Summary& s, std::string* err) {
    s = Summary{};
    File spec;
    if (!spec.open(smspec, err)) return false;
    const auto kws = spec.strings(spec.find("KEYWORDS"));
    std::vector<std::string> names = spec.strings(spec.find("WGNAMES"));
    if (names.empty()) names = spec.strings(spec.find("NAMES"));
    const auto nums = spec.ints("NUMS");
    const auto units = spec.strings(spec.find("UNITS"));
    const auto dim = spec.ints("DIMENS");
    const auto sd = spec.ints("STARTDAT");
    if (dim.size() >= 4) { s.nx = dim[1]; s.ny = dim[2]; s.nz = dim[3]; }
    if (sd.size() >= 3) { s.startDay = sd[0]; s.startMonth = sd[1]; s.startYear = sd[2]; }
    const size_t nv = kws.size();
    if (nv == 0) { if (err) *err = smspec + ": KEYWORDS missing"; return false; }
    std::map<std::string, int> seen;
    s.vars.resize(nv);
    for (size_t i = 0; i < nv; ++i) {
        SummaryVar& v = s.vars[i];
        v.keyword = kws[i];
        v.wgname = i < names.size() ? names[i] : std::string();
        v.num = i < nums.size() ? nums[i] : 0;
        v.unit = i < units.size() ? units[i] : std::string();
        v.category = categoryOf(v.keyword);
        std::string key = v.keyword;
        switch (v.category) {
        case 'W': case 'G': key += validName(v.wgname) ? ":" + v.wgname : std::string(); break;
        case 'B': key += ":" + ijkString(v.num, s.nx, s.ny); break;
        case 'R': case 'A': key += ":" + std::to_string(v.num); break;
        case 'C': key += ":" + v.wgname + ":" + ijkString(v.num, s.nx, s.ny); break;
        case 'S': key += ":" + v.wgname + ":" + std::to_string(v.num); break;
        default: break;
        }
        if ((v.category == 'W' || v.category == 'G' || v.category == 'C' || v.category == 'S') && !validName(v.wgname)) key.clear();
        if (!key.empty() && seen.count(key)) key.clear();   // duplicates: keep the first
        if (!key.empty()) seen[key] = int(i);
        v.key = key;
    }
    // data files
    const std::string base = stripExt(smspec);
    const bool fmt = spec.formatted();
    std::vector<std::string> dataFiles;
    for (const char* e : fmt ? std::initializer_list<const char*>{"FUNSMRY", "funsmry"} : std::initializer_list<const char*>{"UNSMRY", "unsmry"})
        if (fileExists(base + "." + e)) { dataFiles.push_back(base + "." + e); break; }
    if (dataFiles.empty()) {
        for (int k = 0; k <= 9999; ++k) {
            char e[8]; std::snprintf(e, sizeof(e), "%c%04d", fmt ? 'A' : 'S', k);
            std::string p = base + "." + e;
            if (!fileExists(p)) { std::string lo = e; lo[0] = char(std::tolower((unsigned char)lo[0])); p = base + "." + lo; }
            if (fileExists(p)) dataFiles.push_back(p);
            else if (k > 0 && !dataFiles.empty() && k > 2 + std::atoi(dataFiles.back().substr(dataFiles.back().size() - 4).c_str())) break;
            else if (k > 10 && dataFiles.empty()) break;
        }
    }
    if (dataFiles.empty()) { if (err) *err = "no summary data file (UNSMRY / S0001...) next to " + smspec; return false; }
    s.values.assign(nv, {});
    for (const std::string& df : dataFiles) {
        File d;
        if (!d.open(df, nullptr)) continue;
        const auto& A = d.arrays();
        bool sawStep = false;
        for (int i = 0; i < int(A.size()); ++i) {
            if (A[size_t(i)].name == "SEQHDR") {
                if (sawStep && !s.time.empty()) s.reportIndex.push_back(int(s.time.size()) - 1);
                sawStep = false;
            } else if (A[size_t(i)].name == "PARAMS") {
                const auto p = d.floats(i);
                if (p.size() < nv) continue;
                for (size_t v = 0; v < nv; ++v) s.values[v].push_back(p[v]);
                s.time.push_back(0);
                sawStep = true;
            }
        }
        if (sawStep) s.reportIndex.push_back(int(s.time.size()) - 1);
    }
    // time axis
    int tv = -1;
    for (size_t i = 0; i < nv; ++i) if (s.vars[i].keyword == "TIME") { tv = int(i); break; }
    if (tv >= 0) for (size_t k = 0; k < s.time.size(); ++k) s.time[k] = s.values[size_t(tv)][k];
    else {
        int yv = -1;
        for (size_t i = 0; i < nv; ++i) if (s.vars[i].keyword == "YEARS") { yv = int(i); break; }
        if (yv < 0) { if (err) *err = smspec + ": no TIME vector"; return false; }
        for (size_t k = 0; k < s.time.size(); ++k) s.time[k] = s.values[size_t(yv)][k] * 365.25;
    }
    // unit system from the TIME / FPR units
    const int pv = s.find("FPR");
    if (pv >= 0) {
        const std::string u = upper(s.vars[size_t(pv)].unit);
        s.unitSystem = u == "PSIA" ? "FIELD" : u == "BARSA" ? "METRIC" : u == "ATMA" ? "LAB" : "";
    }
    return true;
}

Summary summaryFromRun(const std::vector<SummaryPoint>& pts, const std::vector<std::string>& wellNames, UnitSystem units) {
    Summary s;
    const bool field = units == UnitSystem::Field;
    s.unitSystem = field ? "FIELD" : "METRIC";
    const char* liqR = field ? "STB/DAY" : "SM3/DAY";
    const char* gasR = field ? "MSCF/DAY" : "SM3/DAY";
    const char* liqT = field ? "STB" : "SM3";
    const char* gasT = field ? "MSCF" : "SM3";
    const char* pr = field ? "PSIA" : "BARSA";
    const char* gor = field ? "MSCF/STB" : "SM3/SM3";
    struct Def { std::string kw, wg, unit; std::function<double(const SummaryPoint&)> f; };
    std::vector<Def> defs = {
        {"TIME", "", "DAYS", [](const SummaryPoint& p) { return p.time; }},
        {"FOPR", "", liqR, [](const SummaryPoint& p) { return p.fopr; }},
        {"FWPR", "", liqR, [](const SummaryPoint& p) { return p.fwpr; }},
        {"FGPR", "", gasR, [](const SummaryPoint& p) { return p.fgpr; }},
        {"FLPR", "", liqR, [](const SummaryPoint& p) { return p.fopr + p.fwpr; }},
        {"FWIR", "", liqR, [](const SummaryPoint& p) { return p.fwir; }},
        {"FGIR", "", gasR, [](const SummaryPoint& p) { return p.fgir; }},
        {"FOPT", "", liqT, [](const SummaryPoint& p) { return p.fopt; }},
        {"FWPT", "", liqT, [](const SummaryPoint& p) { return p.fwpt; }},
        {"FGPT", "", gasT, [](const SummaryPoint& p) { return p.fgpt; }},
        {"FWIT", "", liqT, [](const SummaryPoint& p) { return p.fwit; }},
        {"FGIT", "", gasT, [](const SummaryPoint& p) { return p.fgit; }},
        {"FPR", "", pr, [](const SummaryPoint& p) { return p.fpr; }},
        {"FPRP", "", pr, [](const SummaryPoint& p) { return p.fprp; }},
        {"FWCT", "", "", [](const SummaryPoint& p) { return p.fwct; }},
        {"FGOR", "", gor, [](const SummaryPoint& p) { return p.fgor; }},
        {"FOIP", "", liqT, [](const SummaryPoint& p) { return p.oip; }},
        {"FWIP", "", liqT, [](const SummaryPoint& p) { return p.wip; }},
        {"FGIP", "", gasT, [](const SummaryPoint& p) { return p.gip; }},
        {"FCIR", "", field ? "LB/DAY" : "KG/DAY", [](const SummaryPoint& p) { return p.fcir; }},
        {"FCPR", "", field ? "LB/DAY" : "KG/DAY", [](const SummaryPoint& p) { return p.fcpr; }},
        {"FCIT", "", field ? "LB" : "KG", [](const SummaryPoint& p) { return p.fcit; }},
        {"FCPT", "", field ? "LB" : "KG", [](const SummaryPoint& p) { return p.fcpt; }},
        {"FCIP", "", field ? "LB" : "KG", [](const SummaryPoint& p) { return p.fcip; }},
        {"FCAD", "", field ? "LB" : "KG", [](const SummaryPoint& p) { return p.fcad; }},
        {"NEWTON", "", "", [](const SummaryPoint& p) { return double(p.newton); }},
        {"NLINEARS", "", "", [](const SummaryPoint& p) { return double(p.linear); }},
    };
    for (size_t w = 0; w < wellNames.size(); ++w) {
        auto W = [w](const SummaryPoint& p) -> const WellRates* { return w < p.wells.size() ? &p.wells[w] : nullptr; };
        const std::string& nm = wellNames[w];
        defs.push_back({"WOPR", nm, liqR, [W](const SummaryPoint& p) { auto r = W(p); return r ? r->oil : 0.0; }});
        defs.push_back({"WWPR", nm, liqR, [W](const SummaryPoint& p) { auto r = W(p); return r ? r->water : 0.0; }});
        defs.push_back({"WGPR", nm, gasR, [W](const SummaryPoint& p) { auto r = W(p); return r ? r->gas : 0.0; }});
        defs.push_back({"WLPR", nm, liqR, [W](const SummaryPoint& p) { auto r = W(p); return r ? r->oil + r->water : 0.0; }});
        defs.push_back({"WWIR", nm, liqR, [W](const SummaryPoint& p) { auto r = W(p); return r ? r->waterInj : 0.0; }});
        defs.push_back({"WGIR", nm, gasR, [W](const SummaryPoint& p) { auto r = W(p); return r ? r->gasInj : 0.0; }});
        defs.push_back({"WBHP", nm, pr, [W](const SummaryPoint& p) { auto r = W(p); return r && r->open ? r->bhp : 0.0; }});
        defs.push_back({"WWCT", nm, "", [W](const SummaryPoint& p) { auto r = W(p); return r && r->oil + r->water > 0 ? r->water / (r->oil + r->water) : 0.0; }});
        defs.push_back({"WCIR", nm, field ? "LB/DAY" : "KG/DAY", [W](const SummaryPoint& p) { auto r = W(p); return r ? r->polyInj : 0.0; }});
        defs.push_back({"WCPR", nm, field ? "LB/DAY" : "KG/DAY", [W](const SummaryPoint& p) { auto r = W(p); return r ? r->polyProd : 0.0; }});
        defs.push_back({"WGOR", nm, gor, [W](const SummaryPoint& p) { auto r = W(p); return r && r->oil > 0 ? r->gas / r->oil : 0.0; }});
    }
    s.vars.resize(defs.size());
    s.values.assign(defs.size(), {});
    for (size_t i = 0; i < defs.size(); ++i) {
        SummaryVar& v = s.vars[i];
        v.keyword = defs[i].kw; v.wgname = defs[i].wg; v.unit = defs[i].unit;
        v.category = categoryOf(v.keyword);
        v.key = v.wgname.empty() ? v.keyword : v.keyword + ":" + v.wgname;
        s.values[i].reserve(pts.size());
        for (const auto& p : pts) s.values[i].push_back(float(defs[i].f(p)));
    }
    for (const auto& p : pts) s.time.push_back(p.time);
    for (size_t k = 0; k < pts.size(); ++k) s.reportIndex.push_back(int(k));
    return s;
}

// ================================================================== restart
std::vector<SummaryPoint> runFromSummary(const Summary& s, const std::vector<std::string>& wellNames) {
    std::vector<SummaryPoint> pts(s.time.size());
    auto col = [&](const std::string& key) -> const std::vector<float>* { int v = s.find(key); return v < 0 ? nullptr : &s.values[size_t(v)]; };
    auto get = [](const std::vector<float>* c, size_t k) { return c && k < c->size() ? double((*c)[k]) : 0.0; };
    const auto *fopr = col("FOPR"), *fwpr = col("FWPR"), *fgpr = col("FGPR"), *fwir = col("FWIR"), *fgir = col("FGIR");
    const auto *fopt = col("FOPT"), *fwpt = col("FWPT"), *fgpt = col("FGPT"), *fwit = col("FWIT"), *fgit = col("FGIT");
    const auto *fpr = col("FPR"), *fprp = col("FPRP"), *fwct = col("FWCT"), *fgor = col("FGOR");
    const auto *foip = col("FOIP"), *fwip = col("FWIP"), *fgip = col("FGIP"), *newt = col("NEWTON");
    const auto* lin = col("NLINEARS"); if (!lin) lin = col("MLINEARS");
    std::vector<const std::vector<float>*> aqr, aqt;
    for (size_t v = 0; v < s.vars.size(); ++v) {
        if (s.vars[v].keyword == "AAQR") aqr.push_back(&s.values[v]);
        if (s.vars[v].keyword == "AAQT") aqt.push_back(&s.values[v]);
    }
    struct WC { const std::vector<float> *opr, *wpr, *gpr, *wir, *gir, *bhp, *vpr, *oprh, *wprh, *gprh; };
    std::vector<WC> wc;
    for (auto& nm : wellNames) {
        auto c = [&](const char* kw) { return col(std::string(kw) + ":" + nm); };
        wc.push_back({c("WOPR"), c("WWPR"), c("WGPR"), c("WWIR"), c("WGIR"), c("WBHP"), c("WVPR"), c("WOPRH"), c("WWPRH"), c("WGPRH")});
    }
    for (size_t k = 0; k < pts.size(); ++k) {
        SummaryPoint& p = pts[k];
        p.time = s.time[k]; p.dt = k ? s.time[k] - s.time[k - 1] : s.time[k];
        p.newton = int(get(newt, k)); p.linear = int(get(lin, k));
        p.fopr = get(fopr, k); p.fwpr = get(fwpr, k); p.fgpr = get(fgpr, k); p.fwir = get(fwir, k); p.fgir = get(fgir, k);
        p.fopt = get(fopt, k); p.fwpt = get(fwpt, k); p.fgpt = get(fgpt, k); p.fwit = get(fwit, k); p.fgit = get(fgit, k);
        p.fpr = get(fpr, k); p.fprp = fprp ? get(fprp, k) : p.fpr;
        p.fwct = fwct ? get(fwct, k) : (p.fopr + p.fwpr > 0 ? p.fwpr / (p.fopr + p.fwpr) : 0.0);
        p.fgor = fgor ? get(fgor, k) : (p.fopr > 0 ? p.fgpr / p.fopr : 0.0);
        p.oip = get(foip, k); p.wip = get(fwip, k); p.gip = get(fgip, k);
        for (auto* a : aqr) p.faqr += get(a, k);
        for (auto* a : aqt) p.faqt += get(a, k);
        p.faqrS = p.faqr; p.faqtS = p.faqt;
        p.wells.resize(wellNames.size());
        for (size_t w = 0; w < wc.size(); ++w) {
            WellRates& r = p.wells[w]; const WC& c = wc[w];
            r.oil = get(c.opr, k); r.water = get(c.wpr, k); r.gas = get(c.gpr, k);
            r.waterInj = get(c.wir, k); r.gasInj = get(c.gir, k); r.bhp = get(c.bhp, k); r.resv = get(c.vpr, k);
            if (c.oprh) r.oilH = get(c.oprh, k);
            if (c.wprh) r.waterH = get(c.wprh, k);
            if (c.gprh) r.gasH = get(c.gprh, k);
            r.open = r.bhp > 0 || r.oil + r.water + r.gas + r.waterInj + r.gasInj > 0;
            r.mode = Control::BHP;
        }
    }
    return pts;
}

bool Restart::open(const std::vector<std::string>& files, std::string* err) {
    files_.clear(); steps_.clear();
    for (size_t fi = 0; fi < files.size(); ++fi) {
        auto f = std::make_unique<File>();
        std::string e;
        if (!f->open(files[fi], &e)) { if (err) *err = e; continue; }
        const auto& A = f->arrays();
        const int fileIdx = int(files_.size());
        // step boundaries: SEQNUM (unified) or the whole file (non-unified)
        std::vector<int> starts;
        for (int i = 0; i < int(A.size()); ++i) if (A[size_t(i)].name == "SEQNUM") starts.push_back(i);
        int fileReport = 0;
        if (starts.empty()) {
            starts.push_back(0);
            const std::string ext = extOf(files[fi]);
            if (ext.size() == 5) fileReport = std::atoi(ext.c_str() + 1);
        }
        for (size_t s = 0; s < starts.size(); ++s) {
            RestartStep st;
            st.file = fileIdx;
            st.first = starts[s];
            st.last = s + 1 < starts.size() ? starts[s + 1] : int(A.size());
            if (A[size_t(st.first)].name == "SEQNUM") { auto q = f->ints(st.first); st.report = q.empty() ? int(s) : q[0]; }
            else st.report = fileReport;
            const int ih = f->find("INTEHEAD", 0, st.first);
            if (ih >= 0 && ih < st.last) {
                const auto h = f->ints(ih);
                if (h.size() > 66) { st.day = h[64]; st.month = h[65]; st.year = h[66]; }
            }
            const int dh = f->find("DOUBHEAD", 0, st.first);
            if (dh >= 0 && dh < st.last) { const auto d = f->numeric(dh); if (!d.empty()) st.days = d[0]; }
            steps_.push_back(st);
        }
        files_.push_back(std::move(f));
    }
    std::stable_sort(steps_.begin(), steps_.end(), [](const RestartStep& a, const RestartStep& b) { return a.report < b.report; });
    if (steps_.empty()) { if (err && err->empty()) *err = "no restart steps"; return false; }
    return true;
}

std::vector<std::string> Restart::arrayNames(int step) const {
    std::vector<std::string> out;
    if (step < 0 || step >= int(steps_.size())) return out;
    const RestartStep& st = steps_[size_t(step)];
    const auto& A = files_[size_t(st.file)]->arrays();
    bool sw = false, sg = false, so = false;
    for (int i = st.first; i < st.last; ++i) {
        const ArrayInfo& a = A[size_t(i)];
        if (a.count != na_ || na_ == 0) continue;
        if (a.type != Type::REAL && a.type != Type::DOUB && a.type != Type::INTE) continue;
        if (std::find(out.begin(), out.end(), a.name) != out.end()) continue;
        out.push_back(a.name);
        sw |= a.name == "SWAT"; sg |= a.name == "SGAS"; so |= a.name == "SOIL";
    }
    if (!so && (sw || sg)) out.push_back("SOIL");     // derived: 1 - Sw - Sg
    return out;
}

std::vector<float> Restart::load(int step, const std::string& name) const {
    std::vector<float> v;
    if (step < 0 || step >= int(steps_.size())) return v;
    const RestartStep& st = steps_[size_t(step)];
    const File& f = *files_[size_t(st.file)];
    const auto& A = f.arrays();
    auto get = [&](const std::string& nm) -> std::vector<float> {
        for (int i = st.first; i < st.last; ++i) if (A[size_t(i)].name == nm && A[size_t(i)].count == na_) return f.floats(i);
        return {};
    };
    v = get(name);
    if (v.empty() && name == "SOIL") {
        const auto sw = get("SWAT"), sg = get("SGAS");
        if (sw.empty() && sg.empty()) return v;
        v.assign(size_t(na_), 1.0f);
        for (size_t i = 0; i < v.size(); ++i) v[i] -= (sw.empty() ? 0.0f : sw[i]) + (sg.empty() ? 0.0f : sg[i]);
    }
    return v;
}

// ================================================================== case
CaseFiles locateCase(const std::string& anyFile) {
    CaseFiles c;
    c.base = stripExt(anyFile);
    auto pick = [&](std::initializer_list<const char*> exts) -> std::string {
        for (const char* e : exts) {
            if (fileExists(c.base + "." + e)) return c.base + "." + e;
            std::string lo(e); for (char& ch : lo) ch = char(std::tolower((unsigned char)ch));
            if (fileExists(c.base + "." + lo)) return c.base + "." + lo;
        }
        return std::string();
    };
    c.egrid = pick({"EGRID", "FEGRID"});
    c.init = pick({"INIT", "FINIT"});
    c.smspec = pick({"SMSPEC", "FSMSPEC"});
    const std::string unrst = pick({"UNRST", "FUNRST"});
    if (!unrst.empty()) c.restart.push_back(unrst);
    else {
        // non-unified restart files CASE.X0000 / CASE.F0000
        std::error_code ec;
        const fs::path b = fs::u8path(c.base);
        const fs::path dir = b.has_parent_path() ? b.parent_path() : fs::path(".");
        const std::string stem = upper(b.filename().u8string());
        for (auto& de : fs::directory_iterator(dir, ec)) {
            if (!de.is_regular_file(ec)) continue;
            const std::string fn = de.path().filename().u8string();
            const std::string up = upper(fn);
            if (up.size() != stem.size() + 6 || up.compare(0, stem.size(), stem) != 0 || up[stem.size()] != '.') continue;
            const std::string e = up.substr(stem.size() + 1);
            if ((e[0] == 'X' || e[0] == 'F') && std::all_of(e.begin() + 1, e.end(), [](char ch) { return std::isdigit((unsigned char)ch); }))
                c.restart.push_back(de.path().u8string());
        }
        std::sort(c.restart.begin(), c.restart.end());
    }
    return c;
}

std::vector<float> ResultCase::initArray(const std::string& name) const {
    std::vector<float> v;
    if (!hasInit) return v;
    const int i = init.find(name);
    if (i < 0) return v;
    v = init.floats(i);
    const int na = grid ? grid->numActive() : egrid.numActive();
    if (int(v.size()) == egrid.nx * egrid.ny * egrid.nz && int(v.size()) != na && grid) {   // global array (PORV)
        std::vector<float> a(static_cast<size_t>(na));
        for (int k = 0; k < na; ++k) a[size_t(k)] = v[size_t(grid->globalOf[size_t(k)])];
        v.swap(a);
    }
    return v;
}

bool openResultCase(const std::string& anyFile, ResultCase& rc, std::string* err) {
    rc = ResultCase{};
    rc.path = anyFile;
    rc.files = locateCase(anyFile);
    {
        const fs::path b = fs::u8path(rc.files.base);
        rc.name = b.filename().u8string();
    }
    std::string e;
    if (!rc.files.egrid.empty()) {
        if (readEGrid(rc.files.egrid, rc.egrid, &e)) { rc.hasGrid = true; rc.grid = displayGrid(rc.egrid); }
        else rc.warnings.push_back(e);
    }
    const int na = rc.hasGrid ? rc.grid->numActive() : 0;
    if (!rc.files.init.empty()) {
        if (rc.init.open(rc.files.init, &e)) {
            rc.hasInit = true;
            const auto ih = rc.init.ints("INTEHEAD");
            if (ih.size() > 2) rc.unitSystem = ih[2];
            for (const auto& a : rc.init.arrays()) {
                if (a.type != Type::REAL && a.type != Type::DOUB && a.type != Type::INTE) continue;
                const bool ok = (na > 0 && a.count == na) || (a.name == "PORV" && rc.hasGrid && a.count == int64_t(rc.egrid.nx) * rc.egrid.ny * rc.egrid.nz);
                if (ok && std::find(rc.initArrays.begin(), rc.initArrays.end(), a.name) == rc.initArrays.end()) rc.initArrays.push_back(a.name);
            }
        } else rc.warnings.push_back(e);
    }
    if (!rc.files.smspec.empty()) {
        if (readSummary(rc.files.smspec, rc.summary, &e)) rc.hasSummary = true;
        else rc.warnings.push_back(e);
    }
    if (!rc.files.restart.empty()) {
        rc.restart.setNumActive(na);
        if (rc.restart.open(rc.files.restart, &e)) rc.hasRestart = na > 0;
        else rc.warnings.push_back(e);
        if (rc.unitSystem == 0 && rc.hasRestart) {
            // unit system from the first INTEHEAD of the restart is not needed for display
        }
    }
    if (!rc.hasGrid && !rc.hasSummary) {
        if (err) *err = "no ECLIPSE result files found for " + rc.files.base + " (looked for .EGRID/.SMSPEC/.UNRST/.INIT, binary or formatted)";
        return false;
    }
    return true;
}

// ================================================================== comparison
std::vector<Snapshot> snapshotsFromRestart(const ResultCase& rc, int maxSteps) {
    std::vector<Snapshot> out;
    if (!rc.hasRestart) return out;
    const auto& st = rc.restart.steps();
    const int n = maxSteps < 0 ? int(st.size()) : std::min<int>(maxSteps, int(st.size()));
    const size_t na = rc.grid ? size_t(rc.grid->numActive()) : size_t(rc.egrid.numActive());
    for (int k = 0; k < n; ++k) {
        Snapshot sn; sn.time = st[size_t(k)].days;
        sn.pressure = rc.restart.load(k, "PRESSURE");
        if (sn.pressure.size() != na) continue;
        sn.sw = rc.restart.load(k, "SWAT"); sn.sg = rc.restart.load(k, "SGAS"); sn.rs = rc.restart.load(k, "RS");
        if (sn.sw.size() != na) sn.sw.assign(na, 0.f);
        if (sn.sg.size() != na) sn.sg.assign(na, 0.f);
        if (sn.rs.size() != na) sn.rs.assign(na, 0.f);
        sn.rv = rc.restart.load(k, "RV");
        if (sn.rv.size() != na) sn.rv.clear();
        sn.poly = rc.restart.load(k, "POLYMER");
        if (sn.poly.size() != na) sn.poly.clear();
        sn.polyAds = rc.restart.load(k, "PADS");
        if (sn.polyAds.size() != na) sn.polyAds.clear();
        sn.so.resize(na);
        for (size_t c = 0; c < na; ++c) sn.so[c] = std::max(0.f, 1.f - sn.sw[c] - sn.sg[c]);
        out.push_back(std::move(sn));
    }
    return out;
}

std::vector<VectorDiff> compareSummaries(const Summary& ref, const Summary& cmp, const std::vector<std::string>& keysIn) {
    std::vector<std::string> keys = keysIn;
    if (keys.empty()) {
        static const char* fieldKw[] = {"FOPR", "FWPR", "FGPR", "FWIR", "FGIR", "FOPT", "FWPT", "FGPT", "FWIT", "FGIT", "FPR", "FWCT", "FGOR",
                                         "FCIR", "FCPR", "FCIT", "FCPT", "FCIP", "FCAD"};
        for (auto* k : fieldKw) if (ref.find(k) >= 0 && cmp.find(k) >= 0) keys.push_back(k);
        static const char* wellKw[] = {"WBHP", "WOPR", "WWPR", "WGPR", "WWIR", "WGIR", "WWCT", "WGOR", "WCIR", "WCPR"};
        for (auto* wk : wellKw)
            for (const auto& v : ref.vars)
                if (v.keyword == wk && !v.key.empty() && cmp.find(v.key) >= 0) keys.push_back(v.key);
    }
    std::vector<VectorDiff> out;
    const double tEnd = cmp.time.empty() ? 0 : cmp.time.back();
    for (const auto& k : keys) {
        const int a = ref.find(k), b = cmp.find(k);
        if (a < 0 || b < 0 || ref.time.empty()) continue;
        VectorDiff d; d.key = k;
        double sum2 = 0, refMax = 0;
        std::vector<int> idx = ref.reportIndex;
        if (idx.empty()) for (int i = 0; i < int(ref.time.size()); ++i) idx.push_back(i);
        for (int i : idx) {
            const double t = ref.time[size_t(i)];
            if (t <= 0 || t > tEnd + 1e-6) continue;
            const double r = ref.values[size_t(a)][size_t(i)], c = cmp.at(b, t);
            if (!std::isfinite(c)) continue;
            const double e = c - r;
            d.maxAbs = std::max(d.maxAbs, std::fabs(e));
            sum2 += e * e; refMax = std::max(refMax, std::fabs(r));
            ++d.samples;
            d.refFinal = r; d.cmpFinal = c;
        }
        if (d.samples == 0) continue;
        d.rms = std::sqrt(sum2 / d.samples);
        double cmpMax = 0;
        for (int i : idx) { const double t = ref.time[size_t(i)]; if (t > 0 && t <= tEnd + 1e-6) { const double c = cmp.at(b, t); if (std::isfinite(c)) cmpMax = std::max(cmpMax, std::fabs(c)); } }
        const double scale = std::max(refMax, cmpMax);
        d.relRms = scale > 1e-4 ? d.rms / scale : 0.0;      // both ~zero (e.g. FWCT 0 vs 1e-8): no difference
        out.push_back(d);
    }
    return out;
}

} // namespace ecl
} // namespace ressim
