// C36 seismic for the 3-D geomodel scene: out-of-core SEG-Y, complex-trace / coherence attributes, flattening
#include "ressim/seismic3d.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <limits>
#include <mutex>
#include <numeric>
#include <random>

namespace ressim {
namespace geo {

namespace {
const double kPi = 3.14159265358979323846;
int fseek64(std::FILE* f, long long off) {
#ifdef _WIN32
    return _fseeki64(f, off, SEEK_SET);
#else
    return fseeko(f, off_t(off), SEEK_SET);
#endif
}
uint16_t rd16(const unsigned char* p, bool le) { return le ? uint16_t(p[0] | (p[1] << 8)) : uint16_t((p[0] << 8) | p[1]); }
int32_t rd32(const unsigned char* p, bool le) {
    return le ? int32_t(uint32_t(p[0]) | (uint32_t(p[1]) << 8) | (uint32_t(p[2]) << 16) | (uint32_t(p[3]) << 24))
              : int32_t((uint32_t(p[0]) << 24) | (uint32_t(p[1]) << 16) | (uint32_t(p[2]) << 8) | uint32_t(p[3]));
}
void wr16(unsigned char* p, int v, bool le) {
    if (le) { p[0] = (unsigned char)(v & 255); p[1] = (unsigned char)((v >> 8) & 255); }
    else { p[0] = (unsigned char)((v >> 8) & 255); p[1] = (unsigned char)(v & 255); }
}
void wr32(unsigned char* p, int32_t v, bool le) {
    const uint32_t u = uint32_t(v);
    if (le) { p[0] = (unsigned char)u; p[1] = (unsigned char)(u >> 8); p[2] = (unsigned char)(u >> 16); p[3] = (unsigned char)(u >> 24); }
    else { p[0] = (unsigned char)(u >> 24); p[1] = (unsigned char)(u >> 16); p[2] = (unsigned char)(u >> 8); p[3] = (unsigned char)u; }
}
float ibmToFloat(uint32_t u) {
    if ((u & 0x7fffffff) == 0) return 0.0f;
    const int sign = (u >> 31) ? -1 : 1, expo = int((u >> 24) & 0x7f) - 64;
    return float(sign * (double(u & 0x00ffffff) / 16777216.0) * std::pow(16.0, expo));
}
uint32_t floatToIbm(float f) {
    if (f == 0.0f || !std::isfinite(f)) return 0;
    uint32_t sign = f < 0 ? 0x80000000u : 0;
    double a = std::fabs(double(f));
    int e = 0;
    while (a >= 1.0) { a /= 16.0; ++e; }
    while (a < 0.0625) { a *= 16.0; --e; }
    uint32_t frac = uint32_t(std::lround(a * 16777216.0));
    if (frac >= 16777216u) { frac >>= 4; ++e; }
    const int ex = std::clamp(e + 64, 0, 127);
    return sign | (uint32_t(ex) << 24) | (frac & 0x00ffffffu);
}
bool validFmt(int f) { return f == 1 || f == 2 || f == 3 || f == 5 || f == 8; }
int bytesOf(int f) { return f == 3 ? 2 : f == 8 ? 1 : 4; }
char ebcdic(unsigned char c) {
    static const char* tbl =
        "                                                                "
        "           .<(+|&         !$*); -/        ,%_>?         `:#@'=\""
        " abcdefghi       jklmnopqr       ~stuvwxyz                      "
        "{ABCDEFGHI      }JKLMNOPQR      \\ STUVWXYZ      0123456789      ";
    return tbl[c];
}
int gcdPos(int a, int b) { a = std::abs(a); b = std::abs(b); while (b) { const int t = a % b; a = b; b = t; } return a; }

// streaming SEG-Y writer (regular lattice, any supported format / byte order)
struct SegyWriter {
    std::FILE* f = nullptr; int fmt = 5; bool le = false; int ns = 0; double dt = 4, t0 = 0; int seq = 0;
    bool open(const std::string& path, const std::string& text, int ns_, double dt_, double t0_, int fmt_, bool le_) {
        f = std::fopen(path.c_str(), "wb");
        if (!f) return false;
        ns = ns_; dt = dt_; t0 = t0_; fmt = fmt_; le = le_;
        std::vector<unsigned char> th(3200, ' ');
        for (size_t k = 0; k < 3200 && k < text.size(); ++k) th[k] = (unsigned char)text[k];
        std::fwrite(th.data(), 1, 3200, f);
        std::vector<unsigned char> bh(400, 0);
        wr16(&bh[16], int(std::lround(dt * 1000)), le); wr16(&bh[20], ns, le); wr16(&bh[24], fmt, le); wr16(&bh[28], 1, le);
        wr16(&bh[300], 0x0100, le); wr16(&bh[302], 1, le);
        std::fwrite(bh.data(), 1, 400, f);
        return true;
    }
    void trace(int il, int xl, double x, double y, const float* d) {
        std::vector<unsigned char> h(240, 0);
        ++seq;
        wr32(&h[0], seq, le); wr32(&h[4], seq, le); wr32(&h[8], il, le); wr32(&h[20], xl, le);
        wr16(&h[70], -100, le);
        wr32(&h[72], int32_t(std::lround(x * 100)), le); wr32(&h[76], int32_t(std::lround(y * 100)), le);
        wr32(&h[180], int32_t(std::lround(x * 100)), le); wr32(&h[184], int32_t(std::lround(y * 100)), le);
        wr32(&h[188], il, le); wr32(&h[192], xl, le);
        wr16(&h[108], int(std::lround(t0)), le); wr16(&h[114], ns, le); wr16(&h[116], int(std::lround(dt * 1000)), le);
        std::fwrite(h.data(), 1, 240, f);
        const int b = bytesOf(fmt);
        std::vector<unsigned char> o(size_t(b) * size_t(ns));
        for (int s = 0; s < ns; ++s) {
            unsigned char* p = &o[size_t(b) * size_t(s)];
            const float v = d[s];
            if (fmt == 5) { uint32_t u; std::memcpy(&u, &v, 4); wr32(p, int32_t(u), le); }
            else if (fmt == 1) wr32(p, int32_t(floatToIbm(v)), le);
            else if (fmt == 2) wr32(p, int32_t(std::lround(std::clamp(double(v), -2.1e9, 2.1e9))), le);
            else if (fmt == 3) wr16(p, int(std::lround(std::clamp(double(v), -32768.0, 32767.0))), le);
            else p[0] = (unsigned char)(int8_t(std::lround(std::clamp(double(v), -128.0, 127.0))));
        }
        std::fwrite(o.data(), 1, o.size(), f);
    }
    void close() { if (f) std::fclose(f); f = nullptr; }
};
}  // namespace

// ============================================================================ SegyVolume
SegyVolume::~SegyVolume() { close(); }

void SegyVolume::close() {
    if (f_) std::fclose(f_);
    f_ = nullptr;
    off_.clear(); lru_.clear(); cmap_.clear();
    delete static_cast<std::mutex*>(mtx_); mtx_ = nullptr;
}

bool SegyVolume::ij(double x, double y, double& i, double& j) const {
    const double det = ix * jy - iy * jx;
    if (std::fabs(det) < 1e-12) return false;
    const double u = x - x0, v = y - y0;
    i = (u * jy - v * jx) / det; j = (-u * iy + v * ix) / det;
    return i >= -1e-9 && j >= -1e-9 && i <= ni - 1 + 1e-9 && j <= nj - 1 + 1e-9;
}

bool SegyVolume::open(const std::string& p, std::string& err, const std::function<bool(double)>& progress) {
    close();
    path = p;
    f_ = std::fopen(p.c_str(), "rb");
    if (!f_) { err = "cannot open " + p; return false; }
    mtx_ = new std::mutex();
#ifdef _WIN32
    _fseeki64(f_, 0, SEEK_END); fileBytes = _ftelli64(f_);
#else
    fseeko(f_, 0, SEEK_END); fileBytes = (long long)ftello(f_);
#endif
    if (fileBytes < 3600 + 240) { err = p + ": too small for SEG-Y"; close(); return false; }
    std::vector<unsigned char> th(3200), bh(400);
    fseek64(f_, 0);
    if (std::fread(th.data(), 1, 3200, f_) != 3200 || std::fread(bh.data(), 1, 400, f_) != 400) { err = p + ": cannot read the headers"; close(); return false; }
    bool ascii = true;
    for (int k = 0; k < 80; ++k) if (th[size_t(k)] >= 0x80) ascii = false;
    textHeader.resize(3200);
    for (int k = 0; k < 3200; ++k) textHeader[size_t(k)] = ascii ? char(th[size_t(k)]) : ebcdic(th[size_t(k)]);
    int fmt = rd16(&bh[24], false);
    littleEndian = false;
    if (!validFmt(fmt) && validFmt(rd16(&bh[24], true))) { littleEndian = true; fmt = rd16(&bh[24], true); }
    if (!validFmt(fmt)) { err = p + ": sample format " + std::to_string(fmt) + " not supported (1, 2, 3, 5, 8)"; close(); return false; }
    const bool le = littleEndian;
    format = fmt; bps_ = bytesOf(fmt);
    int nsH = rd16(&bh[20], le);
    dt = rd16(&bh[16], le) / 1000.0;
    int ext = int(int16_t(rd16(&bh[304], le)));
    if (ext < 0 || ext > 100) ext = 0;
    start_ = 3600 + 3200LL * ext;
    std::vector<unsigned char> h(240);
    fseek64(f_, start_);
    if (std::fread(h.data(), 1, 240, f_) != 240) { err = p + ": no traces"; close(); return false; }
    if (nsH <= 0) nsH = rd16(&h[114], le);
    if (dt <= 0) dt = rd16(&h[116], le) / 1000.0;
    if (nsH <= 0) { err = p + ": zero samples per trace"; close(); return false; }
    if (dt <= 0) dt = 4;
    ns = nsH;
    t0 = int(int16_t(rd16(&h[108], le)));
    tlen_ = 240 + (long long)bps_ * ns;
    nTraces = (fileBytes - start_) / tlen_;
    struct Tr { int il, xl; double x, y; };
    std::vector<Tr> tr; tr.reserve(size_t(nTraces));
    // sequential header scan in blocks (one read per block of traces)
    const long long blk = std::max<long long>(1, (1 << 22) / tlen_);
    std::vector<unsigned char> buf(size_t(blk * tlen_));
    for (long long t = 0; t < nTraces; t += blk) {
        const long long n = std::min(blk, nTraces - t);
        fseek64(f_, start_ + t * tlen_);
        if (std::fread(buf.data(), 1, size_t(n * tlen_), f_) != size_t(n * tlen_)) { err = p + ": truncated file"; close(); return false; }
        for (long long q = 0; q < n; ++q) {
            const unsigned char* hh = &buf[size_t(q * tlen_)];
            auto B = [](int byte1) { return size_t(std::max(1, std::min(237, byte1)) - 1); };
            int il = rd32(&hh[B(byteIl)], le), xl = rd32(&hh[B(byteXl)], le);
            if (il == 0 && xl == 0) { il = rd32(&hh[8], le); xl = rd32(&hh[20], le); }
            const int sc = int(int16_t(rd16(&hh[B(byteScalar)], le)));
            const double scale = sc == 0 ? 1.0 : sc > 0 ? double(sc) : 1.0 / double(-sc);
            double x = rd32(&hh[B(byteX)], le) * scale, y = rd32(&hh[B(byteY)], le) * scale;
            if (x == 0 && y == 0) { x = rd32(&hh[72], le) * scale; y = rd32(&hh[76], le) * scale; }
            tr.push_back({il, xl, x, y});
        }
        if (progress && !progress(0.8 * double(t + n) / double(nTraces))) { err = "cancelled"; close(); return false; }
    }
    if (tr.empty()) { err = p + ": no traces"; close(); return false; }
    int ilMin = tr[0].il, ilMax = tr[0].il, xlMin = tr[0].xl, xlMax = tr[0].xl;
    for (const Tr& t : tr) { ilMin = std::min(ilMin, t.il); ilMax = std::max(ilMax, t.il); xlMin = std::min(xlMin, t.xl); xlMax = std::max(xlMax, t.xl); }
    int gi = 0, gj = 0;
    for (const Tr& t : tr) { gi = gcdPos(gi, t.il - ilMin); gj = gcdPos(gj, t.xl - xlMin); }
    dil = std::max(1, gi); dxl = std::max(1, gj);
    il0 = ilMin; xl0 = xlMin;
    ni = (ilMax - ilMin) / dil + 1; nj = (xlMax - xlMin) / dxl + 1;
    if (double(ni) * double(nj) > 5e7) { err = p + ": lattice too large (" + std::to_string(ni) + " x " + std::to_string(nj) + ")"; close(); return false; }
    off_.assign(size_t(ni) * size_t(nj), -1);
    // least-squares plan geometry: x = x0 + i ix + j jx, y = y0 + i iy + j jy
    double A[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, bx[3] = {0, 0, 0}, by[3] = {0, 0, 0};
    for (size_t t = 0; t < tr.size(); ++t) {
        const int i = (tr[t].il - il0) / dil, j = (tr[t].xl - xl0) / dxl;
        off_[size_t(i) * size_t(nj) + size_t(j)] = start_ + (long long)t * tlen_ + 240;
        const double r[3] = {1.0, double(i), double(j)};
        for (int a = 0; a < 3; ++a) { for (int b = 0; b < 3; ++b) A[a][b] += r[a] * r[b]; bx[a] += r[a] * tr[t].x; by[a] += r[a] * tr[t].y; }
    }
    // degenerate axes (a single inline / crossline): regularise that column so the solve stays defined
    if (ni == 1) A[1][1] += 1.0;
    if (nj == 1) A[2][2] += 1.0;
    auto solve3 = [](double M[3][3], const double rhs[3], double out[3]) {
        double a[3][4];
        for (int r = 0; r < 3; ++r) { for (int c = 0; c < 3; ++c) a[r][c] = M[r][c]; a[r][3] = rhs[r]; }
        for (int c = 0; c < 3; ++c) {
            int pv = c;
            for (int r = c + 1; r < 3; ++r) if (std::fabs(a[r][c]) > std::fabs(a[pv][c])) pv = r;
            for (int k = 0; k < 4; ++k) std::swap(a[c][k], a[pv][k]);
            if (std::fabs(a[c][c]) < 1e-300) { out[0] = out[1] = out[2] = 0; return; }
            for (int r = 0; r < 3; ++r) if (r != c) { const double fct = a[r][c] / a[c][c]; for (int k = c; k < 4; ++k) a[r][k] -= fct * a[c][k]; }
        }
        for (int r = 0; r < 3; ++r) out[r] = a[r][3] / a[r][r];
    };
    double sx[3], sy[3];
    solve3(A, bx, sx); solve3(A, by, sy);
    x0 = sx[0]; ix = sx[1]; jx = sx[2]; y0 = sy[0]; iy = sy[1]; jy = sy[2];
    if (ni == 1 && nj > 1) { ix = jy; iy = -jx; }     // perpendicular unit-length axis for a 2-D line
    if (nj == 1 && ni > 1) { jx = -iy; jy = ix; }
    if (ni == 1 && nj == 1) { ix = 1; iy = 0; jx = 0; jy = 1; }
    double e2 = 0;
    for (const Tr& t : tr) {
        const int i = (t.il - il0) / dil, j = (t.xl - xl0) / dxl;
        double x, y; xy(i, j, x, y);
        e2 += (x - t.x) * (x - t.x) + (y - t.y) * (y - t.y);
    }
    fitError = std::sqrt(e2 / double(tr.size()));
    missing = (long long)ni * nj - (long long)tr.size();
    if (missing < 0) missing = 0;
    // amplitude statistics from up to 2000 evenly spaced traces
    {
        const long long step = std::max<long long>(1, (long long)tr.size() / 2000);
        std::vector<float> v(static_cast<size_t>(ns)), all;
        double s2 = 0; long long n = 0;
        for (long long t = 0; t < (long long)tr.size(); t += step) {
            std::vector<unsigned char> raw(size_t(bps_) * size_t(ns));
            if (!readRaw(start_ + t * tlen_ + 240, (long long)raw.size(), raw.data())) continue;
            decode(raw.data(), v.data());
            for (float a : v) { s2 += double(a) * a; ++n; all.push_back(std::fabs(a)); }
        }
        ampRms = n ? float(std::sqrt(s2 / double(n))) : 0.0f;
        ampAbsMax = 0;
        for (float a : all) ampAbsMax = std::max(ampAbsMax, a);
        if (progress) progress(1.0);
    }
    const size_t sl = p.find_last_of("/\\");
    name = p.substr(sl == std::string::npos ? 0 : sl + 1);
    return true;
}

bool SegyVolume::readRaw(long long off, long long bytes, unsigned char* dst) const {
    if (!f_) return false;
    std::lock_guard<std::mutex> lk(*static_cast<std::mutex*>(mtx_));
    if (fseek64(f_, off) != 0) return false;
    return std::fread(dst, 1, size_t(bytes), f_) == size_t(bytes);
}

void SegyVolume::decode(const unsigned char* raw, float* out) const {
    const bool le = littleEndian;
    for (int s = 0; s < ns; ++s) {
        const unsigned char* p = raw + size_t(bps_) * size_t(s);
        switch (format) {
            case 1: out[s] = ibmToFloat(uint32_t(rd32(p, le))); break;
            case 2: out[s] = float(rd32(p, le)); break;
            case 3: out[s] = float(int16_t(rd16(p, le))); break;
            case 8: out[s] = float(int8_t(p[0])); break;
            default: { const uint32_t u = uint32_t(rd32(p, le)); float v; std::memcpy(&v, &u, 4); out[s] = std::isfinite(v) ? v : 0.0f; }
        }
    }
}

bool SegyVolume::trace(int i, int j, float* out) const {
    if (i < 0 || j < 0 || i >= ni || j >= nj || !f_) { std::fill(out, out + ns, 0.0f); return false; }
    const long long key = (long long)i * nj + j;
    const long long off = off_[size_t(key)];
    if (off < 0) { std::fill(out, out + ns, 0.0f); return false; }
    {
        std::lock_guard<std::mutex> lk(*static_cast<std::mutex*>(mtx_));
        auto it = cmap_.find(key);
        if (it != cmap_.end()) {
            lru_.splice(lru_.begin(), lru_, it->second);
            std::copy(it->second->second.begin(), it->second->second.end(), out);
            return true;
        }
    }
    std::vector<unsigned char> raw(size_t(bps_) * size_t(ns));
    if (!readRaw(off, (long long)raw.size(), raw.data())) { std::fill(out, out + ns, 0.0f); return false; }
    decode(raw.data(), out);
    std::lock_guard<std::mutex> lk(*static_cast<std::mutex*>(mtx_));
    if (cacheTraces > 0 && cmap_.find(key) == cmap_.end()) {
        lru_.emplace_front(key, std::vector<float>(out, out + ns));
        cmap_[key] = lru_.begin();
        while (lru_.size() > cacheTraces) { cmap_.erase(lru_.back().first); lru_.pop_back(); }
    }
    return true;
}

std::vector<float> SegyVolume::inlineSlice(int i) const {
    std::vector<float> out(size_t(nj) * size_t(ns), 0.0f);
    if (i < 0 || i >= ni || !f_) return out;
    // runs of traces that are contiguous in the file are read with one call (inline-sorted files: the whole line)
    std::vector<unsigned char> raw;
    int j = 0;
    while (j < nj) {
        const long long o = off_[size_t(i) * size_t(nj) + size_t(j)];
        if (o < 0) { ++j; continue; }
        int e = j + 1;
        while (e < nj && off_[size_t(i) * size_t(nj) + size_t(e)] == o + (long long)(e - j) * tlen_) ++e;
        const long long bytes = (long long)(e - j - 1) * tlen_ + (long long)bps_ * ns;
        raw.resize(size_t(bytes));
        if (readRaw(o, bytes, raw.data()))
            for (int q = j; q < e; ++q) decode(&raw[size_t((long long)(q - j) * tlen_)], &out[size_t(q) * size_t(ns)]);
        j = e;
    }
    return out;
}

std::vector<float> SegyVolume::xlineSlice(int j) const {
    std::vector<float> out(size_t(ni) * size_t(ns), 0.0f);
    if (j < 0 || j >= nj) return out;
    for (int i = 0; i < ni; ++i) trace(i, j, &out[size_t(i) * size_t(ns)]);
    return out;
}

std::vector<float> SegyVolume::timeSlice(int s) const {
    std::vector<float> out(size_t(ni) * size_t(nj), 0.0f);
    if (s < 0 || s >= ns || !f_) return out;
    unsigned char raw[4];
    float tmp[1];
    SegyVolume const* self = this;
    for (size_t k = 0; k < off_.size(); ++k) {
        if (off_[k] < 0) continue;
        if (!readRaw(off_[k] + (long long)bps_ * s, bps_, raw)) continue;
        // decode one sample (decode() handles ns samples; do it inline)
        const bool le = self->littleEndian;
        switch (format) {
            case 1: tmp[0] = ibmToFloat(uint32_t(rd32(raw, le))); break;
            case 2: tmp[0] = float(rd32(raw, le)); break;
            case 3: tmp[0] = float(int16_t(rd16(raw, le))); break;
            case 8: tmp[0] = float(int8_t(raw[0])); break;
            default: { const uint32_t u = uint32_t(rd32(raw, le)); std::memcpy(&tmp[0], &u, 4); if (!std::isfinite(tmp[0])) tmp[0] = 0; }
        }
        out[k] = tmp[0];
    }
    return out;
}

SeismicCube SegyVolume::subCube(int i0, int i1, int j0, int j1, int s0, int s1) const {
    SeismicCube C;
    i0 = std::clamp(i0, 0, ni - 1); i1 = std::clamp(i1, i0, ni - 1);
    j0 = std::clamp(j0, 0, nj - 1); j1 = std::clamp(j1, j0, nj - 1);
    s0 = std::clamp(s0, 0, ns - 1); s1 = std::clamp(s1, s0, ns - 1);
    C.name = name; C.textHeader = textHeader; C.format = format;
    C.ni = i1 - i0 + 1; C.nj = j1 - j0 + 1; C.ns = s1 - s0 + 1;
    C.il0 = il0 + i0 * dil; C.xl0 = xl0 + j0 * dxl; C.dil = dil; C.dxl = dxl;
    C.t0 = t0 + s0 * dt; C.dt = dt;
    xy(i0, j0, C.x0, C.y0);
    C.ix = ix; C.iy = iy; C.jx = jx; C.jy = jy;
    C.a.assign(size_t(C.ni) * size_t(C.nj) * size_t(C.ns), 0.0f);
    for (int i = i0; i <= i1; ++i) {
        const std::vector<float> L = inlineSlice(i);
        for (int j = j0; j <= j1; ++j)
            std::copy(L.begin() + long(size_t(j) * size_t(ns) + size_t(s0)), L.begin() + long(size_t(j) * size_t(ns) + size_t(s1) + 1),
                      C.a.begin() + long((size_t(i - i0) * size_t(C.nj) + size_t(j - j0)) * size_t(C.ns)));
    }
    return C;
}

SeismicCube SegyVolume::proxy(int maxDimPlan, int maxDimTime, const std::function<bool(double)>& progress) const {
    SeismicCube C;
    if (!f_ || ni < 1) return C;
    const int sp = std::max(1, int(std::ceil(double(std::max(ni, nj)) / std::max(2, maxDimPlan))));
    const int st = std::max(1, int(std::ceil(double(ns) / std::max(2, maxDimTime))));
    C.name = name + " (proxy)"; C.textHeader = textHeader; C.format = format;
    C.ni = (ni - 1) / sp + 1; C.nj = (nj - 1) / sp + 1; C.ns = ns / st;
    if (C.ns < 1) C.ns = 1;
    C.il0 = il0; C.xl0 = xl0; C.dil = dil * sp; C.dxl = dxl * sp;
    C.dt = dt * st; C.t0 = t0 + 0.5 * (st - 1) * dt;   // box average over st samples: centre of the box
    C.x0 = x0; C.y0 = y0; C.ix = ix * sp; C.iy = iy * sp; C.jx = jx * sp; C.jy = jy * sp;
    C.a.assign(size_t(C.ni) * size_t(C.nj) * size_t(C.ns), 0.0f);
    for (int a = 0; a < C.ni; ++a) {
        const std::vector<float> L = inlineSlice(a * sp);
        for (int b = 0; b < C.nj; ++b) {
            const float* tr = &L[size_t(b * sp) * size_t(ns)];
            float* o = &C.a[(size_t(a) * size_t(C.nj) + size_t(b)) * size_t(C.ns)];
            for (int s = 0; s < C.ns; ++s) {
                double m = 0;
                for (int q = 0; q < st; ++q) m += tr[s * st + q];
                o[s] = float(m / st);
            }
        }
        if (progress && !progress(double(a + 1) / C.ni)) return SeismicCube();
    }
    return C;
}

// ============================================================================ attributes
const char* seisAttrName(SeisAttr a) {
    static const char* n[] = {"Amplitude", "Envelope (reflection strength)", "Instantaneous phase", "Instantaneous frequency", "Cosine of phase",
                              "Sweetness", "RMS amplitude", "Relative acoustic impedance", "Coherence (semblance)", "Variance (edge)"};
    const int k = int(a);
    return k >= 0 && k < int(SeisAttr::Count) ? n[k] : "?";
}
const char* seisAttrUnit(SeisAttr a) {
    switch (a) { case SeisAttr::Phase: return "deg"; case SeisAttr::Frequency: return "Hz"; case SeisAttr::Coherence: case SeisAttr::Variance: case SeisAttr::CosPhase: return "-"; default: return "amp"; }
}
bool seisAttrMultiTrace(SeisAttr a) { return a == SeisAttr::Coherence || a == SeisAttr::Variance; }

namespace {
void fft(std::vector<std::complex<double>>& a, bool inverse) {
    const size_t n = a.size();
    for (size_t i = 1, j = 0; i < n; ++i) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (size_t len = 2; len <= n; len <<= 1) {
        const double ang = 2 * kPi / double(len) * (inverse ? 1 : -1);
        const std::complex<double> wl(std::cos(ang), std::sin(ang));
        for (size_t i = 0; i < n; i += len) {
            std::complex<double> w(1);
            for (size_t k = 0; k < len / 2; ++k) {
                const std::complex<double> u = a[i + k], v = a[i + k + len / 2] * w;
                a[i + k] = u + v; a[i + k + len / 2] = u - v;
                w *= wl;
            }
        }
    }
    if (inverse) for (auto& x : a) x /= double(n);
}
}  // namespace

void hilbert(const float* x, int n, std::vector<float>& q) {
    q.assign(size_t(std::max(n, 0)), 0.0f);
    if (n < 2) return;
    size_t N = 1;
    while (N < size_t(2 * n)) N <<= 1;
    std::vector<std::complex<double>> a(N, 0.0);
    for (int i = 0; i < n; ++i) a[size_t(i)] = x[i];
    fft(a, false);
    for (size_t k = 1; k < N / 2; ++k) a[k] *= 2.0;
    for (size_t k = N / 2 + 1; k < N; ++k) a[k] = 0.0;
    fft(a, true);
    for (int i = 0; i < n; ++i) q[size_t(i)] = float(a[size_t(i)].imag());
}

void traceAttribute(SeisAttr att, const float* x, int n, double dtMs, int win, float* out) {
    if (n <= 0) return;
    win = std::max(1, win);
    if (att == SeisAttr::Amplitude) { std::copy(x, x + n, out); return; }
    if (att == SeisAttr::Rms) {
        std::vector<double> c(size_t(n) + 1, 0.0);
        for (int i = 0; i < n; ++i) c[size_t(i) + 1] = c[size_t(i)] + double(x[i]) * x[i];
        for (int i = 0; i < n; ++i) { const int a = std::max(0, i - win), b = std::min(n - 1, i + win); out[i] = float(std::sqrt((c[size_t(b) + 1] - c[size_t(a)]) / (b - a + 1))); }
        return;
    }
    if (att == SeisAttr::RelImpedance) {
        // band-limited integration of the trace (reflectivity ~ d ln Z / 2): running sum minus its moving average
        std::vector<double> cs(static_cast<size_t>(n)), c2(size_t(n) + 1, 0.0);
        double s = 0; for (int i = 0; i < n; ++i) { s += x[i]; cs[size_t(i)] = s; c2[size_t(i) + 1] = c2[size_t(i)] + s; }
        const int W = std::max(10, 8 * win);
        for (int i = 0; i < n; ++i) { const int a = std::max(0, i - W), b = std::min(n - 1, i + W); out[i] = float(cs[size_t(i)] - (c2[size_t(b) + 1] - c2[size_t(a)]) / (b - a + 1)); }
        return;
    }
    std::vector<float> q;
    hilbert(x, n, q);
    std::vector<double> env(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) env[size_t(i)] = std::sqrt(double(x[i]) * x[i] + double(q[size_t(i)]) * q[size_t(i)]);
    double emax = 0; for (double e : env) emax = std::max(emax, e);
    const double tiny = 1e-9 * std::max(emax, 1e-30);
    switch (att) {
        case SeisAttr::Envelope: for (int i = 0; i < n; ++i) out[i] = float(env[size_t(i)]); return;
        case SeisAttr::Phase: for (int i = 0; i < n; ++i) out[i] = float(std::atan2(double(q[size_t(i)]), double(x[i])) * 180.0 / kPi); return;
        case SeisAttr::CosPhase: for (int i = 0; i < n; ++i) out[i] = env[size_t(i)] > tiny ? float(x[i] / env[size_t(i)]) : 0.0f; return;
        default: break;
    }
    // instantaneous frequency, envelope^2-weighted smoothing over +-win samples
    const double dts = dtMs / 1000.0, nyq = 0.5 / std::max(dts, 1e-9);
    std::vector<double> fw(static_cast<size_t>(n)), w(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        // phase difference of the analytic signal across the neighbours, arg(z[b] conj(z[a])) / (2 pi (b - a) dt):
        // exact for a sinusoid (a central difference of x and q is biased by sin(w dt) / (w dt))
        const int a = std::max(0, i - 1), b = std::min(n - 1, i + 1);
        const double h = (b - a) * dts;
        const double re = double(x[b]) * x[a] + double(q[size_t(b)]) * q[size_t(a)], im = double(q[size_t(b)]) * x[a] - double(x[b]) * q[size_t(a)];
        const double e2 = env[size_t(i)] * env[size_t(i)];
        const double f = (e2 > tiny * tiny && h > 0) ? std::atan2(im, re) / (2 * kPi * h) : 0.0;
        w[size_t(i)] = e2; fw[size_t(i)] = e2 * std::clamp(f, 0.0, nyq);
    }
    std::vector<double> cw(size_t(n) + 1, 0.0), cf(size_t(n) + 1, 0.0);
    for (int i = 0; i < n; ++i) { cw[size_t(i) + 1] = cw[size_t(i)] + w[size_t(i)]; cf[size_t(i) + 1] = cf[size_t(i)] + fw[size_t(i)]; }
    for (int i = 0; i < n; ++i) {
        const int a = std::max(0, i - win), b = std::min(n - 1, i + win);
        const double sw = cw[size_t(b) + 1] - cw[size_t(a)];
        const double f = sw > 0 ? (cf[size_t(b) + 1] - cf[size_t(a)]) / sw : 0.0;
        if (att == SeisAttr::Frequency) out[i] = float(f);
        else out[i] = float(env[size_t(i)] / std::sqrt(std::max(f, 1.0)));   // sweetness = envelope / sqrt(frequency)
    }
}

double semblance(const std::vector<const float*>& tr, int s0, int s1) {
    const size_t K = tr.size();
    if (K == 0 || s1 < s0) return 0;
    double num = 0, den = 0;
    for (int s = s0; s <= s1; ++s) {
        double sum = 0, sq = 0;
        for (const float* t : tr) { sum += t[s]; sq += double(t[s]) * t[s]; }
        num += sum * sum; den += sq;
    }
    return den > 0 ? num / (double(K) * den) : 1.0;
}

std::vector<float> sectionAttribute(SeisAttr a, const std::vector<float>& sec, int ntr, int ns, double dtMs, int win,
                                    const std::vector<float>* prev, const std::vector<float>* next) {
    std::vector<float> out(size_t(ntr) * size_t(ns), 0.0f);
    if (ntr <= 0 || ns <= 0 || sec.size() < size_t(ntr) * size_t(ns)) return out;
    win = std::max(1, win);
    if (!seisAttrMultiTrace(a)) {
#pragma omp parallel for schedule(dynamic, 8)
        for (int t = 0; t < ntr; ++t) traceAttribute(a, &sec[size_t(t) * size_t(ns)], ns, dtMs, win, &out[size_t(t) * size_t(ns)]);
        return out;
    }
#pragma omp parallel for schedule(dynamic, 8)
    for (int t = 0; t < ntr; ++t) {
        std::vector<const float*> tr;
        for (int d = -1; d <= 1; ++d) {
            const int u = t + d;
            if (u < 0 || u >= ntr) continue;
            tr.push_back(&sec[size_t(u) * size_t(ns)]);
            if (prev && prev->size() >= sec.size()) tr.push_back(&(*prev)[size_t(u) * size_t(ns)]);
            if (next && next->size() >= sec.size()) tr.push_back(&(*next)[size_t(u) * size_t(ns)]);
        }
        // running sums over the window: num(s) = (sum_k u)^2, den(s) = sum_k u^2
        std::vector<double> cn(size_t(ns) + 1, 0.0), cd(size_t(ns) + 1, 0.0);
        for (int s = 0; s < ns; ++s) {
            double sum = 0, sq = 0;
            for (const float* p : tr) { sum += p[s]; sq += double(p[s]) * p[s]; }
            cn[size_t(s) + 1] = cn[size_t(s)] + sum * sum; cd[size_t(s) + 1] = cd[size_t(s)] + sq;
        }
        const double K = double(tr.size());
        for (int s = 0; s < ns; ++s) {
            const int lo = std::max(0, s - win), hi = std::min(ns - 1, s + win);
            const double num = cn[size_t(hi) + 1] - cn[size_t(lo)], den = cd[size_t(hi) + 1] - cd[size_t(lo)];
            const double S = den > 1e-30 ? std::clamp(num / (K * den), 0.0, 1.0) : 1.0;
            out[size_t(t) * size_t(ns) + size_t(s)] = float(a == SeisAttr::Coherence ? S : 1.0 - S);
        }
    }
    return out;
}

SeismicCube cubeAttribute(const SeismicCube& c, SeisAttr a, int win) {
    SeismicCube o = c;
    o.name = c.name + " - " + seisAttrName(a);
    const size_t L = size_t(c.nj) * size_t(c.ns);
    for (int i = 0; i < c.ni; ++i) {
        const std::vector<float> sec(c.a.begin() + long(size_t(i) * L), c.a.begin() + long(size_t(i + 1) * L));
        std::vector<float> pv, nx;
        if (seisAttrMultiTrace(a)) {
            if (i > 0) pv.assign(c.a.begin() + long(size_t(i - 1) * L), c.a.begin() + long(size_t(i) * L));
            if (i + 1 < c.ni) nx.assign(c.a.begin() + long(size_t(i + 1) * L), c.a.begin() + long(size_t(i + 2) * L));
        }
        const std::vector<float> r = sectionAttribute(a, sec, c.nj, c.ns, c.dt, win, pv.empty() ? nullptr : &pv, nx.empty() ? nullptr : &nx);
        std::copy(r.begin(), r.end(), o.a.begin() + long(size_t(i) * L));
    }
    return o;
}

bool volumeAttribute(const SegyVolume& v, SeisAttr a, int win, const std::string& outPath, std::string& err, const std::function<bool(double)>& progress) {
    if (!v.isOpen()) { err = "no SEG-Y volume open"; return false; }
    SegyWriter W;
    std::string txt = std::string("C 1 ResSim / ResGeo attribute volume: ") + seisAttrName(a) + " of " + v.name;
    if (!W.open(outPath, txt, v.ns, v.dt, v.t0, 5, false)) { err = "cannot write " + outPath; return false; }
    const bool multi = seisAttrMultiTrace(a);
    std::vector<float> prev, cur = v.inlineSlice(0), next;
    for (int i = 0; i < v.ni; ++i) {
        if (multi) next = i + 1 < v.ni ? v.inlineSlice(i + 1) : std::vector<float>();
        const std::vector<float> r = sectionAttribute(a, cur, v.nj, v.ns, v.dt, win, prev.empty() ? nullptr : &prev, next.empty() ? nullptr : &next);
        for (int j = 0; j < v.nj; ++j) { double x, y; v.xy(i, j, x, y); W.trace(v.il0 + i * v.dil, v.xl0 + j * v.dxl, x, y, &r[size_t(j) * size_t(v.ns)]); }
        if (multi) { prev.swap(cur); cur.swap(next); }
        else if (i + 1 < v.ni) cur = v.inlineSlice(i + 1);
        if (progress && !progress(double(i + 1) / v.ni)) { W.close(); err = "cancelled"; return false; }
    }
    W.close();
    return true;
}

std::vector<float> flattenSection(const std::vector<float>& sec, int ntr, int ns, const std::vector<double>& hz, double ref) {
    std::vector<float> out(size_t(ntr) * size_t(ns), 0.0f);
    for (int t = 0; t < ntr; ++t) {
        const double sh = (size_t(t) < hz.size() && std::isfinite(hz[size_t(t)])) ? hz[size_t(t)] - ref : 0.0;
        const float* in = &sec[size_t(t) * size_t(ns)];
        float* o = &out[size_t(t) * size_t(ns)];
        for (int s = 0; s < ns; ++s) {
            const double u = s + sh;
            if (u < 0 || u > ns - 1) continue;
            const int k = std::min(int(u), ns - 2 < 0 ? 0 : ns - 2);
            const double f = u - k;
            o[s] = ns > 1 ? float((1 - f) * in[k] + f * in[k + 1]) : in[0];
        }
    }
    return out;
}

float amplitudeClip(const std::vector<float>& v, double pct) {
    if (v.empty()) return 1.0f;
    const size_t step = std::max<size_t>(1, v.size() / 200000);
    std::vector<float> a;
    for (size_t k = 0; k < v.size(); k += step) if (std::isfinite(v[k])) a.push_back(std::fabs(v[k]));
    if (a.empty()) return 1.0f;
    const size_t m = std::min(a.size() - 1, size_t(std::clamp(pct, 0.0, 1.0) * double(a.size() - 1)));
    std::nth_element(a.begin(), a.begin() + long(m), a.end());
    return std::max(a[m], 1e-30f);
}

bool writeSyntheticSegy(const std::string& path, int ni, int nj, int ns, double dtMs, int fmt, bool le, std::string& err) {
    if (!validFmt(fmt)) { err = "unsupported format"; return false; }
    SegyWriter W;
    if (!W.open(path, "C 1 ResSim synthetic SEG-Y: faulted anticline, Ricker 25 Hz", ns, dtMs, 0, fmt, le)) { err = "cannot write " + path; return false; }
    // reflectors: 6 horizons, anticline crest at the centre, normal fault along i = j + ni/4 with 20 ms throw
    const double f0 = 25.0;
    auto ricker = [&](double tms) { const double a = kPi * f0 * tms / 1000.0; return (1 - 2 * a * a) * std::exp(-a * a); };
    const double amp = fmt == 8 ? 40.0 : fmt == 3 ? 8000.0 : fmt == 2 ? 1e6 : 1.0;
    std::vector<float> tr(static_cast<size_t>(ns));
    for (int i = 0; i < ni; ++i)
        for (int j = 0; j < nj; ++j) {
            std::fill(tr.begin(), tr.end(), 0.0f);
            const double u = (i - 0.5 * ni) / std::max(1.0, 0.35 * ni), w = (j - 0.5 * nj) / std::max(1.0, 0.35 * nj);
            const double crest = 60.0 * std::exp(-(u * u + w * w));
            const double thr = (i - j - ni / 4.0) > 0 ? 20.0 : 0.0;
            for (int k = 0; k < 6; ++k) {
                const double t = 0.18 * ns * dtMs + k * 0.12 * ns * dtMs - crest + thr + 4.0 * std::sin(0.3 * i + 0.2 * j + k);
                const double r = (k % 2 ? -1.0 : 1.0) * (0.6 + 0.08 * k);
                const int sc = int(t / dtMs);
                for (int s = std::max(0, sc - 12); s < std::min(ns, sc + 13); ++s) tr[size_t(s)] += float(amp * r * ricker(s * dtMs - t));
            }
            const double x = 1000.0 + 25.0 * i * 0.8 - 25.0 * j * 0.6, y = 2000.0 + 25.0 * i * 0.6 + 25.0 * j * 0.8;   // rotated survey
            W.trace(100 + 2 * i, 300 + j, x, y, tr.data());
        }
    W.close();
    return true;
}

}  // namespace geo
}  // namespace ressim
