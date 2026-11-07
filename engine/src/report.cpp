// C33 report generator (report.h): DOCX / PPTX / HTML / Markdown writers, native Office charts, PNG maps, ZIP container
#include "ressim/report.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>

namespace ressim {
namespace report {

// ================================================================ utilities
namespace {
std::string xesc(const std::string& s) {
    std::string o;
    o.reserve(s.size() + 8);
    for (unsigned char c : s) {
        switch (c) {
        case '&': o += "&amp;"; break;
        case '<': o += "&lt;"; break;
        case '>': o += "&gt;"; break;
        case '"': o += "&quot;"; break;
        default:
            if (c < 32 && c != '\t' && c != '\n') o += ' ';
            else o += char(c);
        }
    }
    return o;
}
std::string num(double v, int prec = 6) {
    if (!std::isfinite(v)) return "";
    char b[48];
    std::snprintf(b, sizeof b, "%.*g", prec, v);
    return b;
}
std::string lower(std::string s) { for (auto& c : s) c = char(std::tolower(static_cast<unsigned char>(c))); return s; }
std::string hex6(uint32_t rgb) { char b[8]; std::snprintf(b, sizeof b, "%06X", rgb & 0xFFFFFFu); return b; }
const char* B64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
std::string base64(const std::vector<uint8_t>& d) {
    std::string o;
    o.reserve((d.size() + 2) / 3 * 4);
    for (size_t i = 0; i < d.size(); i += 3) {
        uint32_t v = uint32_t(d[i]) << 16;
        if (i + 1 < d.size()) v |= uint32_t(d[i + 1]) << 8;
        if (i + 2 < d.size()) v |= d[i + 2];
        o += B64[(v >> 18) & 63]; o += B64[(v >> 12) & 63];
        o += i + 1 < d.size() ? B64[(v >> 6) & 63] : '=';
        o += i + 2 < d.size() ? B64[v & 63] : '=';
    }
    return o;
}
// at most n points per series (keeps the first / last and the extremes of each bucket)
void decimate(const std::vector<double>& x, const std::vector<double>& y, size_t n, std::vector<double>& xo, std::vector<double>& yo) {
    xo.clear(); yo.clear();
    const size_t m = std::min(x.size(), y.size());
    if (m <= n) { xo.assign(x.begin(), x.begin() + long(m)); yo.assign(y.begin(), y.begin() + long(m)); return; }
    const size_t buckets = std::max<size_t>(1, n / 2);
    for (size_t b = 0; b < buckets; ++b) {
        const size_t a0 = b * m / buckets, a1 = std::max(a0 + 1, (b + 1) * m / buckets);
        size_t imin = a0, imax = a0;
        for (size_t i = a0; i < a1; ++i) {
            if (!std::isfinite(y[i])) continue;
            if (!std::isfinite(y[imin]) || y[i] < y[imin]) imin = i;
            if (!std::isfinite(y[imax]) || y[i] > y[imax]) imax = i;
        }
        const size_t lo = std::min(imin, imax), hi = std::max(imin, imax);
        xo.push_back(x[lo]); yo.push_back(y[lo]);
        if (hi != lo) { xo.push_back(x[hi]); yo.push_back(y[hi]); }
    }
    if (xo.back() != x[m - 1]) { xo.push_back(x[m - 1]); yo.push_back(y[m - 1]); }
}
}  // namespace

uint32_t paletteColor(int i) {
    static const uint32_t P[] = {0x1F77B4, 0xD62728, 0x2CA02C, 0xFF7F0E, 0x9467BD, 0x8C564B, 0xE377C2, 0x17BECF, 0x7F7F7F, 0xBCBD22};
    return P[(i % 10 + 10) % 10];
}

// ================================================================ CRC / ZIP / PNG
uint32_t crc32(const uint8_t* p, size_t n, uint32_t crc) {
    static uint32_t T[256];
    static bool init = false;
    if (!init) {
        for (uint32_t i = 0; i < 256; ++i) { uint32_t c = i; for (int k = 0; k < 8; ++k) c = (c & 1) ? 0xEDB88320u ^ (c >> 1) : c >> 1; T[i] = c; }
        init = true;
    }
    crc = ~crc;
    for (size_t i = 0; i < n; ++i) crc = T[(crc ^ p[i]) & 0xFF] ^ (crc >> 8);
    return ~crc;
}

namespace {
// ---- raw DEFLATE (RFC 1951): greedy LZ77 (32 KiB window, hash chains) + the fixed Huffman code; inflate for stored /
// fixed blocks (what this writer produces; used by readZip in the tests)
const int kLenBase[29] = {3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258};
const int kLenExtra[29] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0};
const int kDistBase[30] = {1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577};
const int kDistExtra[30] = {0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13};
struct BitOut {
    std::string out; uint32_t acc = 0; int n = 0;
    void bits(uint32_t v, int k) { acc |= v << n; n += k; while (n >= 8) { out += char(acc & 0xFF); acc >>= 8; n -= 8; } }
    void huff(uint32_t code, int len) { uint32_t r = 0; for (int i = 0; i < len; ++i) r |= ((code >> i) & 1u) << (len - 1 - i); bits(r, len); }
    void flush() { if (n > 0) { out += char(acc & 0xFF); acc = 0; n = 0; } }
};
void litLen(BitOut& b, int sym) {
    if (sym < 144) b.huff(0x30 + uint32_t(sym), 8);
    else if (sym < 256) b.huff(0x190 + uint32_t(sym - 144), 9);
    else if (sym < 280) b.huff(uint32_t(sym - 256), 7);
    else b.huff(0xC0 + uint32_t(sym - 280), 8);
}
}  // namespace

std::string deflateRaw(const uint8_t* d, size_t n) {
    BitOut b;
    b.bits(1, 1); b.bits(1, 2);   // final block, fixed Huffman
    const int WBITS = 15, W = 1 << WBITS, HBITS = 15, HS = 1 << HBITS;
    std::vector<int> head(size_t(HS), -1), prev(size_t(W), -1);
    auto hash = [&](size_t i) { return int(((uint32_t(d[i]) << 10) ^ (uint32_t(d[i + 1]) << 5) ^ d[i + 2]) & (HS - 1)); };
    size_t i = 0;
    auto insert = [&](size_t k) { if (k + 2 < n) { const int h = hash(k); prev[k & (W - 1)] = head[size_t(h)]; head[size_t(h)] = int(k); } };
    while (i < n) {
        int bestLen = 0, bestDist = 0;
        if (i + 2 < n) {
            int cand = head[size_t(hash(i))];
            int chain = 48;
            const size_t maxLen = std::min<size_t>(258, n - i);
            while (cand >= 0 && chain-- > 0 && i - size_t(cand) <= size_t(W - 1)) {
                if (d[size_t(cand) + size_t(bestLen)] == d[i + size_t(bestLen)] || bestLen == 0) {
                    size_t l = 0;
                    while (l < maxLen && d[size_t(cand) + l] == d[i + l]) ++l;
                    if (int(l) > bestLen) { bestLen = int(l); bestDist = int(i - size_t(cand)); if (l == maxLen) break; }
                }
                const int nx = prev[size_t(cand) & (W - 1)];
                if (nx >= cand) break;
                cand = nx;
            }
        }
        if (bestLen >= 3) {
            int lc = 0; while (lc < 28 && kLenBase[lc + 1] <= bestLen) ++lc;
            litLen(b, 257 + lc);
            if (kLenExtra[lc]) b.bits(uint32_t(bestLen - kLenBase[lc]), kLenExtra[lc]);
            int dc = 0; while (dc < 29 && kDistBase[dc + 1] <= bestDist) ++dc;
            b.huff(uint32_t(dc), 5);
            if (kDistExtra[dc]) b.bits(uint32_t(bestDist - kDistBase[dc]), kDistExtra[dc]);
            for (int k = 0; k < bestLen; ++k) insert(i + size_t(k));
            i += size_t(bestLen);
        } else {
            litLen(b, d[i]);
            insert(i);
            ++i;
        }
    }
    litLen(b, 256);
    b.flush();
    return b.out;
}

bool inflateRaw(const std::string& in, std::string& out) {
    size_t pos = 0; uint32_t acc = 0; int n = 0;
    auto need = [&](int k) { while (n < k) { if (pos >= in.size()) return false; acc |= uint32_t(uint8_t(in[pos++])) << n; n += 8; } return true; };
    auto bits = [&](int k, uint32_t& v) { if (!need(k)) return false; v = acc & ((1u << k) - 1u); acc >>= k; n -= k; return true; };
    auto huffBit = [&](uint32_t& code, int k) { for (int i = 0; i < k; ++i) { uint32_t bt; if (!bits(1, bt)) return false; code = (code << 1) | bt; } return true; };
    for (;;) {
        uint32_t fin, type;
        if (!bits(1, fin) || !bits(2, type)) return false;
        if (type == 0) {
            acc = 0; n = 0;   // byte align
            if (pos + 4 > in.size()) return false;
            const uint32_t len = uint32_t(uint8_t(in[pos])) | uint32_t(uint8_t(in[pos + 1])) << 8;
            pos += 4;
            if (pos + len > in.size()) return false;
            out.append(in, pos, len); pos += len;
        } else if (type == 1) {
            for (;;) {
                uint32_t code = 0; int sym;
                if (!huffBit(code, 7)) return false;
                if (code <= 0x17) sym = int(code) + 256;
                else {
                    if (!huffBit(code, 1)) return false;
                    if (code >= 0x30 && code <= 0xBF) sym = int(code - 0x30);
                    else if (code >= 0xC0 && code <= 0xC7) sym = int(code - 0xC0) + 280;
                    else { if (!huffBit(code, 1)) return false; if (code < 0x190 || code > 0x1FF) return false; sym = int(code - 0x190) + 144; }
                }
                if (sym < 256) { out += char(sym); continue; }
                if (sym == 256) break;
                const int lc = sym - 257;
                if (lc > 28) return false;
                uint32_t e = 0;
                if (kLenExtra[lc] && !bits(kLenExtra[lc], e)) return false;
                const int len = kLenBase[lc] + int(e);
                uint32_t dc = 0;
                if (!huffBit(dc, 5) || dc > 29) return false;
                uint32_t de = 0;
                if (kDistExtra[dc] && !bits(kDistExtra[dc], de)) return false;
                const size_t dist = size_t(kDistBase[dc]) + de;
                if (dist > out.size()) return false;
                const size_t from = out.size() - dist;
                for (int k = 0; k < len; ++k) out += out[from + size_t(k)];
            }
        } else return false;
        if (fin) return true;
    }
}

namespace {
void put16(std::string& s, uint32_t v) { s += char(v & 0xFF); s += char((v >> 8) & 0xFF); }
void put32(std::string& s, uint32_t v) { put16(s, v & 0xFFFF); put16(s, v >> 16); }
uint32_t get16(const std::string& s, size_t o) { return uint32_t(uint8_t(s[o])) | uint32_t(uint8_t(s[o + 1])) << 8; }
uint32_t get32(const std::string& s, size_t o) { return get16(s, o) | get16(s, o + 2) << 16; }
}  // namespace

bool writeZip(const std::string& path, const std::vector<ZipEntry>& entries, std::string& err) {
    std::string out, cdir;
    const uint32_t dosTime = 0, dosDate = (uint32_t(2026 - 1980) << 9) | (1u << 5) | 1u;
    for (const ZipEntry& e : entries) {
        const uint32_t crc = crc32(reinterpret_cast<const uint8_t*>(e.data.data()), e.data.size());
        const uint32_t off = uint32_t(out.size());
        // PNG images are already compressed: stored; everything else deflated (when that is smaller)
        const bool png = e.name.size() > 4 && e.name.compare(e.name.size() - 4, 4, ".png") == 0;
        std::string z;
        if (!png && !e.data.empty()) z = deflateRaw(reinterpret_cast<const uint8_t*>(e.data.data()), e.data.size());
        const bool defl = !z.empty() && z.size() < e.data.size();
        const std::string& body = defl ? z : e.data;
        const uint32_t method = defl ? 8 : 0;
        out += "PK\x03\x04"; put16(out, 20); put16(out, 0); put16(out, method); put16(out, dosTime); put16(out, dosDate);
        put32(out, crc); put32(out, uint32_t(body.size())); put32(out, uint32_t(e.data.size()));
        put16(out, uint32_t(e.name.size())); put16(out, 0);
        out += e.name; out += body;
        cdir += "PK\x01\x02"; put16(cdir, 20); put16(cdir, 20); put16(cdir, 0); put16(cdir, method); put16(cdir, dosTime); put16(cdir, dosDate);
        put32(cdir, crc); put32(cdir, uint32_t(body.size())); put32(cdir, uint32_t(e.data.size()));
        put16(cdir, uint32_t(e.name.size())); put16(cdir, 0); put16(cdir, 0); put16(cdir, 0); put16(cdir, 0); put32(cdir, 0); put32(cdir, off);
        cdir += e.name;
    }
    const uint32_t cdOff = uint32_t(out.size());
    out += cdir;
    out += "PK\x05\x06"; put16(out, 0); put16(out, 0); put16(out, uint32_t(entries.size())); put16(out, uint32_t(entries.size()));
    put32(out, uint32_t(cdir.size())); put32(out, cdOff); put16(out, 0);
    std::ofstream f(path, std::ios::binary);
    if (!f) { err = "cannot write " + path; return false; }
    f.write(out.data(), std::streamsize(out.size()));
    if (!f) { err = "write failed: " + path; return false; }
    return true;
}

std::vector<ZipEntry> readZip(const std::string& path, std::string& err) {
    std::vector<ZipEntry> r;
    std::ifstream f(path, std::ios::binary);
    if (!f) { err = "cannot open " + path; return r; }
    std::stringstream ss; ss << f.rdbuf();
    const std::string s = ss.str();
    size_t o = 0;
    while (o + 30 <= s.size() && s.compare(o, 4, "PK\x03\x04") == 0) {
        const uint32_t method = get16(s, o + 8), crc = get32(s, o + 14), csz = get32(s, o + 18);
        const uint32_t nl = get16(s, o + 26), xl = get16(s, o + 28);
        ZipEntry e; e.name = s.substr(o + 30, nl);
        const size_t d = o + 30 + nl + xl;
        if (d + csz > s.size()) { err = "truncated zip"; return r; }
        if (method == 0) e.data = s.substr(d, csz);
        else if (method == 8) { if (!inflateRaw(s.substr(d, csz), e.data)) { err = "cannot inflate " + e.name; return r; } }
        else { err = "unsupported compression method in " + e.name; return r; }
        if (crc32(reinterpret_cast<const uint8_t*>(e.data.data()), e.data.size()) != crc) { err = "CRC mismatch in " + e.name; return r; }
        r.push_back(std::move(e));
        o = d + csz;
    }
    return r;
}

std::vector<uint8_t> encodePng(const std::vector<uint8_t>& rgb, int w, int h) {
    std::vector<uint8_t> raw;
    raw.reserve(size_t(h) * (size_t(w) * 3 + 1));
    for (int y = 0; y < h; ++y) {
        raw.push_back(0);
        raw.insert(raw.end(), rgb.begin() + long(size_t(y) * size_t(w) * 3), rgb.begin() + long(size_t(y + 1) * size_t(w) * 3));
    }
    // zlib stream (RFC 1950) around a raw deflate stream
    uint32_t a = 1, b = 0;
    for (uint8_t c : raw) { a = (a + c) % 65521u; b = (b + a) % 65521u; }
    const std::string def = deflateRaw(raw.data(), raw.size());
    std::vector<uint8_t> z = {0x78, 0x01};
    z.insert(z.end(), def.begin(), def.end());
    const uint32_t ad = (b << 16) | a;
    z.push_back(uint8_t(ad >> 24)); z.push_back(uint8_t(ad >> 16)); z.push_back(uint8_t(ad >> 8)); z.push_back(uint8_t(ad));
    std::vector<uint8_t> png = {0x89, 'P', 'N', 'G', 0x0D, 0x0A, 0x1A, 0x0A};
    auto chunk = [&](const char* type, const std::vector<uint8_t>& data) {
        const uint32_t n = uint32_t(data.size());
        png.push_back(uint8_t(n >> 24)); png.push_back(uint8_t(n >> 16)); png.push_back(uint8_t(n >> 8)); png.push_back(uint8_t(n));
        std::vector<uint8_t> td(type, type + 4);
        td.insert(td.end(), data.begin(), data.end());
        png.insert(png.end(), td.begin(), td.end());
        const uint32_t c = crc32(td.data(), td.size());
        png.push_back(uint8_t(c >> 24)); png.push_back(uint8_t(c >> 16)); png.push_back(uint8_t(c >> 8)); png.push_back(uint8_t(c));
    };
    std::vector<uint8_t> ihdr = {uint8_t(w >> 24), uint8_t(w >> 16), uint8_t(w >> 8), uint8_t(w), uint8_t(h >> 24), uint8_t(h >> 16), uint8_t(h >> 8), uint8_t(h), 8, 2, 0, 0, 0};
    chunk("IHDR", ihdr);
    chunk("IDAT", z);
    chunk("IEND", {});
    return png;
}

namespace {
// Google "turbo" colour map, polynomial approximation (A. Mikhailov 2019)
void turbo(double t, uint8_t& r, uint8_t& g, uint8_t& b) {
    t = std::min(1.0, std::max(0.0, t));
    const double R = 0.13572138 + t * (4.61539260 + t * (-42.66032258 + t * (132.13108234 + t * (-152.94239396 + t * 59.28637943))));
    const double G = 0.09140261 + t * (2.19418839 + t * (4.84296658 + t * (-14.18503333 + t * (4.27729857 + t * 2.82956604))));
    const double B = 0.10667330 + t * (12.64194608 + t * (-60.58204836 + t * (110.36276771 + t * (-89.90310912 + t * 27.34824973))));
    auto c = [](double v) { return uint8_t(std::lround(255.0 * std::min(1.0, std::max(0.0, v)))); };
    r = c(R); g = c(G); b = c(B);
}
}  // namespace

Image layerMap(const std::vector<float>& val, const std::vector<int>& globalOf, int nx, int ny, int k, double lo, double hi,
               int scale, const std::string& caption) {
    Image im;
    scale = std::max(1, scale);
    const int W = nx * scale, Hm = ny * scale, bar = std::max(10, 2 * scale), gap = std::max(6, scale);
    const int H = Hm + gap + bar;
    std::vector<uint8_t> rgb(size_t(W) * size_t(H) * 3, 255);
    std::vector<int> act(size_t(nx) * size_t(ny), -1);
    const int nxy = nx * ny;
    for (size_t a = 0; a < globalOf.size() && a < val.size(); ++a) {
        const int g = globalOf[a];
        if (g >= nxy * (k + 1) || g < nxy * k) continue;
        act[size_t(g - nxy * k)] = int(a);
    }
    const double span = hi > lo ? hi - lo : 1.0;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const int a = act[size_t(j * nx + i)];
            uint8_t r = 200, g = 200, b = 200;
            if (a >= 0 && std::isfinite(val[size_t(a)])) turbo((val[size_t(a)] - lo) / span, r, g, b);
            for (int yy = 0; yy < scale; ++yy)
                for (int xx = 0; xx < scale; ++xx) {
                    const size_t o = (size_t(j * scale + yy) * size_t(W) + size_t(i * scale + xx)) * 3;   // row j = 0 at the top (ECLIPSE J down)
                    rgb[o] = r; rgb[o + 1] = g; rgb[o + 2] = b;
                }
        }
    for (int x = 0; x < W; ++x) {
        uint8_t r, g, b;
        turbo(double(x) / std::max(1, W - 1), r, g, b);
        for (int y = Hm + gap; y < H; ++y) { const size_t o = (size_t(y) * size_t(W) + size_t(x)) * 3; rgb[o] = r; rgb[o + 1] = g; rgb[o + 2] = b; }
    }
    im.png = encodePng(rgb, W, H);
    im.w = W; im.h = H;
    im.caption = caption + " (colour bar " + num(lo, 4) + " .. " + num(hi, 4) + ", grey = inactive)";
    return im;
}

// ================================================================ SVG chart (HTML)
namespace {
void niceTicks(double lo, double hi, int target, std::vector<double>& t) {
    t.clear();
    if (!(hi > lo)) { hi = lo + 1; }
    const double raw = (hi - lo) / std::max(1, target);
    const double mag = std::pow(10.0, std::floor(std::log10(raw)));
    const double r = raw / mag;
    const double step = (r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10) * mag;
    for (double v = std::ceil(lo / step - 1e-9) * step; v <= hi + step * 1e-6; v += step) t.push_back(std::fabs(v) < step * 1e-9 ? 0.0 : v);
}
struct Range { double lo = 1e300, hi = -1e300; void add(double v) { if (std::isfinite(v)) { lo = std::min(lo, v); hi = std::max(hi, v); } } bool ok() const { return lo <= hi; } };
}  // namespace

std::string svgChart(const Chart& c, int W, int H) {
    const int L = 70, R = 70, T = 30, B = 62;
    const double pw = W - L - R, ph = H - T - B;
    Range rx, ry, ry2;
    bool hasY2 = false;
    for (const auto& s : c.series) {
        for (double v : s.x) rx.add(v);
        for (double v : s.y) { if (c.logY && !(v > 0)) continue; (s.y2 ? ry2 : ry).add(v); }
        hasY2 = hasY2 || s.y2;
    }
    if (!rx.ok()) { rx.lo = 0; rx.hi = 1; }
    auto fixR = [&](Range& r) {
        if (!r.ok()) { r.lo = 0; r.hi = 1; }
        if (c.logY) { r.lo = std::pow(10.0, std::floor(std::log10(std::max(1e-300, r.lo)))); r.hi = std::pow(10.0, std::ceil(std::log10(std::max(r.hi, r.lo * 10)))); return; }
        if (r.lo > 0 && r.lo < 0.35 * r.hi) r.lo = 0;
        if (r.hi == r.lo) { r.hi += std::fabs(r.hi) * 0.1 + 1; r.lo -= std::fabs(r.lo) * 0.1; }
        const double pad = 0.04 * (r.hi - r.lo); r.hi += pad; if (r.lo != 0) r.lo -= pad;
    };
    fixR(ry); fixR(ry2);
    if (rx.hi == rx.lo) rx.hi = rx.lo + 1;
    auto X = [&](double v) { return L + (v - rx.lo) / (rx.hi - rx.lo) * pw; };
    auto Yr = [&](const Range& r, double v) {
        if (c.logY) return T + ph - (std::log10(v) - std::log10(r.lo)) / (std::log10(r.hi) - std::log10(r.lo)) * ph;
        return T + ph - (v - r.lo) / (r.hi - r.lo) * ph;
    };
    std::ostringstream o;
    o << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 " << W << ' ' << H << "\" width=\"" << W << "\" height=\"" << H
      << "\" font-family=\"Segoe UI,Arial,sans-serif\" font-size=\"11\">";
    o << "<rect x=\"0\" y=\"0\" width=\"" << W << "\" height=\"" << H << "\" fill=\"#fff\"/>";
    if (!c.title.empty()) o << "<text x=\"" << W / 2 << "\" y=\"18\" text-anchor=\"middle\" font-size=\"13\" font-weight=\"600\">" << xesc(c.title) << "</text>";
    std::vector<double> tk;
    niceTicks(rx.lo, rx.hi, 8, tk);
    for (double v : tk) {
        if (v < rx.lo - 1e-12 || v > rx.hi + 1e-12) continue;
        o << "<line x1=\"" << X(v) << "\" y1=\"" << T << "\" x2=\"" << X(v) << "\" y2=\"" << T + ph << "\" stroke=\"#e6e6e6\"/>";
        o << "<text x=\"" << X(v) << "\" y=\"" << T + ph + 14 << "\" text-anchor=\"middle\" fill=\"#444\">" << num(v, 6) << "</text>";
    }
    auto yTicks = [&](const Range& r, bool right) {
        std::vector<double> t;
        if (c.logY) { for (double v = r.lo; v <= r.hi * 1.0001; v *= 10) t.push_back(v); }
        else niceTicks(r.lo, r.hi, 6, t);
        for (double v : t) {
            if (v < r.lo - 1e-12 * std::fabs(r.hi) || v > r.hi * (1 + 1e-12) + 1e-300) continue;
            const double y = Yr(r, v);
            if (!right) o << "<line x1=\"" << L << "\" y1=\"" << y << "\" x2=\"" << L + pw << "\" y2=\"" << y << "\" stroke=\"#e6e6e6\"/>";
            o << "<text x=\"" << (right ? L + pw + 6 : L - 6) << "\" y=\"" << y + 4 << "\" text-anchor=\"" << (right ? "start" : "end") << "\" fill=\"#444\">" << num(v, 5) << "</text>";
        }
    };
    yTicks(ry, false);
    if (hasY2) yTicks(ry2, true);
    o << "<rect x=\"" << L << "\" y=\"" << T << "\" width=\"" << pw << "\" height=\"" << ph << "\" fill=\"none\" stroke=\"#888\"/>";
    if (!c.xLabel.empty()) o << "<text x=\"" << L + pw / 2 << "\" y=\"" << T + ph + 30 << "\" text-anchor=\"middle\">" << xesc(c.xLabel) << "</text>";
    if (!c.yLabel.empty()) o << "<text transform=\"translate(14," << T + ph / 2 << ") rotate(-90)\" text-anchor=\"middle\">" << xesc(c.yLabel) << "</text>";
    if (hasY2 && !c.y2Label.empty()) o << "<text transform=\"translate(" << W - 12 << "," << T + ph / 2 << ") rotate(90)\" text-anchor=\"middle\">" << xesc(c.y2Label) << "</text>";
    o << "<defs><clipPath id=\"pc\"><rect x=\"" << L << "\" y=\"" << T << "\" width=\"" << pw << "\" height=\"" << ph << "\"/></clipPath></defs>";
    int legendX = L, legendY = H - 12;
    for (size_t si = 0; si < c.series.size(); ++si) {
        const Series& s = c.series[si];
        const std::string col = "#" + hex6(s.rgb ? s.rgb : paletteColor(int(si)));
        std::vector<double> xs, ys;
        decimate(s.x, s.y, 1500, xs, ys);
        std::ostringstream pts;
        bool open = false;
        o << "<g clip-path=\"url(#pc)\">";
        for (size_t i = 0; i < xs.size(); ++i) {
            const bool ok = std::isfinite(xs[i]) && std::isfinite(ys[i]) && (!c.logY || ys[i] > 0);
            if (!ok) { if (open) { o << "<polyline fill=\"none\" stroke=\"" << col << "\" stroke-width=\"1.8\"" << (s.dashed ? " stroke-dasharray=\"6 4\"" : "") << " points=\"" << pts.str() << "\"/>"; pts.str(""); open = false; } continue; }
            pts << X(xs[i]) << ',' << Yr(s.y2 ? ry2 : ry, ys[i]) << ' ';
            open = true;
            if (s.markers) o << "<circle cx=\"" << X(xs[i]) << "\" cy=\"" << Yr(s.y2 ? ry2 : ry, ys[i]) << "\" r=\"2.6\" fill=\"" << col << "\"/>";
        }
        if (open) o << "<polyline fill=\"none\" stroke=\"" << col << "\" stroke-width=\"1.8\"" << (s.dashed ? " stroke-dasharray=\"6 4\"" : "") << " points=\"" << pts.str() << "\"/>";
        o << "</g>";
        const std::string lab = s.name + (s.y2 ? " (right)" : "");
        const int lw = 30 + int(lab.size()) * 6;
        if (legendX + lw > W - 10) { legendX = L; legendY -= 0; }
        o << "<line x1=\"" << legendX << "\" y1=\"" << legendY - 4 << "\" x2=\"" << legendX + 20 << "\" y2=\"" << legendY - 4 << "\" stroke=\"" << col
          << "\" stroke-width=\"2\"" << (s.dashed ? " stroke-dasharray=\"5 3\"" : "") << "/>";
        o << "<text x=\"" << legendX + 24 << "\" y=\"" << legendY << "\">" << xesc(lab) << "</text>";
        legendX += lw;
    }
    o << "</svg>";
    return o.str();
}

// ================================================================ Office chart part (DOCX + PPTX)
namespace {
std::string richText(const std::string& t, int sz, bool bold, bool rot) {
    std::ostringstream o;
    o << "<c:tx><c:rich><a:bodyPr" << (rot ? " rot=\"-5400000\" vert=\"horz\"" : "") << "/><a:lstStyle/><a:p><a:pPr><a:defRPr sz=\"" << sz << "\" b=\"" << (bold ? 1 : 0)
      << "\"/></a:pPr><a:r><a:rPr lang=\"en-US\" sz=\"" << sz << "\" b=\"" << (bold ? 1 : 0) << "\"/><a:t>" << xesc(t) << "</a:t></a:r></a:p></c:rich></c:tx><c:overlay val=\"0\"/>";
    return o.str();
}
std::string chartXml(const Chart& c) {
    std::ostringstream o;
    o << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
         "<c:chartSpace xmlns:c=\"http://schemas.openxmlformats.org/drawingml/2006/chart\" xmlns:a=\"http://schemas.openxmlformats.org/drawingml/2006/main\" "
         "xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\"><c:roundedCorners val=\"0\"/><c:chart>";
    if (!c.title.empty()) o << "<c:title>" << richText(c.title, 1200, true, false) << "</c:title><c:autoTitleDeleted val=\"0\"/>";
    else o << "<c:autoTitleDeleted val=\"1\"/>";
    o << "<c:plotArea><c:layout/>";
    bool hasY2 = false;
    for (const auto& s : c.series) hasY2 = hasY2 || s.y2;
    int idx = 0;
    for (int pass = 0; pass < (hasY2 ? 2 : 1); ++pass) {
        o << "<c:scatterChart><c:scatterStyle val=\"lineMarker\"/><c:varyColors val=\"0\"/>";
        for (size_t si = 0; si < c.series.size(); ++si) {
            const Series& s = c.series[si];
            if (s.y2 != (pass == 1)) continue;
            std::vector<double> xs, ys;
            decimate(s.x, s.y, 400, xs, ys);
            const std::string col = hex6(s.rgb ? s.rgb : paletteColor(int(si)));
            o << "<c:ser><c:idx val=\"" << idx << "\"/><c:order val=\"" << idx << "\"/><c:tx><c:v>" << xesc(s.name + (s.y2 ? " (right)" : "")) << "</c:v></c:tx>";
            o << "<c:spPr><a:ln w=\"22225\" cap=\"rnd\"><a:solidFill><a:srgbClr val=\"" << col << "\"/></a:solidFill>" << (s.dashed ? "<a:prstDash val=\"dash\"/>" : "")
              << "<a:round/></a:ln></c:spPr>";
            if (s.markers) o << "<c:marker><c:symbol val=\"circle\"/><c:size val=\"4\"/><c:spPr><a:solidFill><a:srgbClr val=\"" << col << "\"/></a:solidFill></c:spPr></c:marker>";
            else o << "<c:marker><c:symbol val=\"none\"/></c:marker>";
            std::ostringstream xv, yv;
            int n = 0;
            for (size_t i = 0; i < xs.size(); ++i) {
                if (!std::isfinite(xs[i]) || !std::isfinite(ys[i]) || (c.logY && !(ys[i] > 0))) continue;
                xv << "<c:pt idx=\"" << n << "\"><c:v>" << num(xs[i], 9) << "</c:v></c:pt>";
                yv << "<c:pt idx=\"" << n << "\"><c:v>" << num(ys[i], 9) << "</c:v></c:pt>";
                ++n;
            }
            o << "<c:xVal><c:numLit><c:formatCode>General</c:formatCode><c:ptCount val=\"" << n << "\"/>" << xv.str() << "</c:numLit></c:xVal>";
            o << "<c:yVal><c:numLit><c:formatCode>General</c:formatCode><c:ptCount val=\"" << n << "\"/>" << yv.str() << "</c:numLit></c:yVal>";
            o << "<c:smooth val=\"0\"/></c:ser>";
            ++idx;
        }
        o << "<c:axId val=\"" << (pass ? 50030 : 50010) << "\"/><c:axId val=\"" << (pass ? 50040 : 50020) << "\"/></c:scatterChart>";
    }
    auto ax = [&](int id, int cross, const char* pos, const std::string& title, bool grid, bool del, bool log, bool crossMax) {
        o << "<c:valAx><c:axId val=\"" << id << "\"/><c:scaling>" << (log ? "<c:logBase val=\"10\"/>" : "") << "<c:orientation val=\"minMax\"/></c:scaling><c:delete val=\""
          << (del ? 1 : 0) << "\"/><c:axPos val=\"" << pos << "\"/>";
        if (grid) o << "<c:majorGridlines><c:spPr><a:ln w=\"6350\"><a:solidFill><a:srgbClr val=\"E6E6E6\"/></a:solidFill></a:ln></c:spPr></c:majorGridlines>";
        if (!title.empty() && !del) o << "<c:title>" << richText(title, 1000, false, *pos == 'l' || *pos == 'r') << "</c:title>";
        o << "<c:numFmt formatCode=\"General\" sourceLinked=\"0\"/><c:majorTickMark val=\"out\"/><c:minorTickMark val=\"none\"/><c:tickLblPos val=\"nextTo\"/>"
             "<c:spPr><a:ln w=\"9525\"><a:solidFill><a:srgbClr val=\"888888\"/></a:solidFill></a:ln></c:spPr>"
             "<c:txPr><a:bodyPr/><a:lstStyle/><a:p><a:pPr><a:defRPr sz=\"900\"/></a:pPr><a:endParaRPr lang=\"en-US\"/></a:p></c:txPr>"
          << "<c:crossAx val=\"" << cross << "\"/><c:crosses val=\"" << (crossMax ? "max" : "autoZero") << "\"/><c:crossBetween val=\"midCat\"/></c:valAx>";
    };
    ax(50010, 50020, "b", c.xLabel, false, false, false, false);
    ax(50020, 50010, "l", c.yLabel, true, false, c.logY, false);
    if (hasY2) {
        ax(50030, 50040, "b", "", false, true, false, false);
        ax(50040, 50030, "r", c.y2Label, false, false, c.logY, true);
    }
    o << "<c:spPr><a:noFill/><a:ln><a:noFill/></a:ln></c:spPr></c:plotArea><c:legend><c:legendPos val=\"b\"/><c:overlay val=\"0\"/>"
         "<c:txPr><a:bodyPr/><a:lstStyle/><a:p><a:pPr><a:defRPr sz=\"900\"/></a:pPr><a:endParaRPr lang=\"en-US\"/></a:p></c:txPr></c:legend>"
         "<c:plotVisOnly val=\"1\"/><c:dispBlanksAs val=\"gap\"/></c:chart>"
         "<c:spPr><a:solidFill><a:srgbClr val=\"FFFFFF\"/></a:solidFill><a:ln><a:noFill/></a:ln></c:spPr>"
         "<c:txPr><a:bodyPr/><a:lstStyle/><a:p><a:pPr><a:defRPr sz=\"1000\"><a:latin typeface=\"Calibri\"/></a:defRPr></a:pPr><a:endParaRPr lang=\"en-US\"/></a:p></c:txPr></c:chartSpace>";
    return o.str();
}
std::string corePropsXml(const Doc& d) {
    std::ostringstream o;
    o << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n<cp:coreProperties xmlns:cp=\"http://schemas.openxmlformats.org/package/2006/metadata/core-properties\" "
         "xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:dcterms=\"http://purl.org/dc/terms/\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">"
      << "<dc:title>" << xesc(d.title) << "</dc:title><dc:creator>" << xesc(d.author.empty() ? "ResSim Studio" : d.author) << "</dc:creator>"
      << "<cp:keywords>ResSim report</cp:keywords></cp:coreProperties>";
    return o.str();
}
std::string appPropsXml() {
    return "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n<Properties xmlns=\"http://schemas.openxmlformats.org/officeDocument/2006/extended-properties\" "
           "xmlns:vt=\"http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes\"><Application>ResSim Studio</Application></Properties>";
}
std::string today() {
    std::time_t t = std::time(nullptr);
    char b[32];
    std::strftime(b, sizeof b, "%Y-%m-%d", std::localtime(&t));
    return b;
}
}  // namespace

// ================================================================ DOCX
namespace {
struct DocxBuilder {
    std::ostringstream body;
    std::vector<ZipEntry> media;              // word/media/..., word/charts/...
    std::vector<std::pair<std::string, std::string>> rels;   // id, type|target
    std::vector<std::string> overrides;
    int nextRel = 10, docPr = 1, nImg = 0, nChart = 0;
    static constexpr long long TEXT_W = 6120000;   // EMU (A4, 2 cm margins)

    std::string run(const std::string& t, bool bold = false, int sz = 0, const char* color = nullptr) {
        std::ostringstream o;
        o << "<w:r>";
        if (bold || sz || color) {
            o << "<w:rPr>" << (bold ? "<w:b/>" : "");
            if (color) o << "<w:color w:val=\"" << color << "\"/>";
            if (sz) o << "<w:sz w:val=\"" << sz << "\"/><w:szCs w:val=\"" << sz << "\"/>";
            o << "</w:rPr>";
        }
        // keep line breaks
        size_t a = 0;
        while (true) {
            const size_t e = t.find('\n', a);
            o << "<w:t xml:space=\"preserve\">" << xesc(t.substr(a, e == std::string::npos ? std::string::npos : e - a)) << "</w:t>";
            if (e == std::string::npos) break;
            o << "<w:br/>";
            a = e + 1;
        }
        o << "</w:r>";
        return o.str();
    }
    void para(const std::string& text, const char* style = nullptr, bool bold = false) {
        body << "<w:p>";
        if (style) body << "<w:pPr><w:pStyle w:val=\"" << style << "\"/></w:pPr>";
        body << run(text, bold) << "</w:p>";
    }
    std::string rel(const std::string& type, const std::string& target) {
        const std::string id = "rId" + std::to_string(nextRel++);
        rels.push_back({id, type + "|" + target});
        return id;
    }
    void table(const Table& t) {
        const size_t nc = std::max<size_t>(1, t.header.size());
        body << "<w:tbl><w:tblPr><w:tblStyle w:val=\"ResTable\"/><w:tblW w:w=\"5000\" w:type=\"pct\"/>"
                "<w:tblBorders><w:top w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"A6A6A6\"/><w:left w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"A6A6A6\"/>"
                "<w:bottom w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"A6A6A6\"/><w:right w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"A6A6A6\"/>"
                "<w:insideH w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"D9D9D9\"/><w:insideV w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"D9D9D9\"/></w:tblBorders>"
                "<w:tblCellMar><w:left w:w=\"80\" w:type=\"dxa\"/><w:right w:w=\"80\" w:type=\"dxa\"/></w:tblCellMar><w:tblLook w:val=\"04A0\" w:firstRow=\"1\"/></w:tblPr><w:tblGrid>";
        for (size_t i = 0; i < nc; ++i) body << "<w:gridCol w:w=\"" << 9638 / nc << "\"/>";
        body << "</w:tblGrid>";
        auto cell = [&](const std::string& s, const char* fill, bool bold, const char* color, bool right) {
            body << "<w:tc><w:tcPr><w:tcW w:w=\"" << 5000 / nc << "\" w:type=\"pct\"/>";
            if (fill) body << "<w:shd w:val=\"clear\" w:color=\"auto\" w:fill=\"" << fill << "\"/>";
            body << "</w:tcPr><w:p><w:pPr><w:pStyle w:val=\"TableText\"/>" << (right ? "<w:jc w:val=\"right\"/>" : "") << "</w:pPr>" << run(s, bold, 0, color) << "</w:p></w:tc>";
        };
        if (!t.header.empty()) {
            body << "<w:tr><w:trPr><w:tblHeader/></w:trPr>";
            for (const auto& h : t.header) cell(h, "1F4E79", true, "FFFFFF", false);
            body << "</w:tr>";
        }
        for (size_t r = 0; r < t.rows.size(); ++r) {
            const bool hl = std::find(t.highlightRows.begin(), t.highlightRows.end(), int(r)) != t.highlightRows.end();
            body << "<w:tr><w:trPr><w:cantSplit/></w:trPr>";
            for (size_t c = 0; c < nc; ++c) {
                const std::string s = c < t.rows[r].size() ? t.rows[r][c] : std::string();
                const bool numeric = !s.empty() && (std::isdigit(static_cast<unsigned char>(s[0])) || ((s[0] == '-' || s[0] == '+') && s.size() > 1 && std::isdigit(static_cast<unsigned char>(s[1]))));
                cell(s, hl ? "FCE4D6" : (r % 2 ? "F2F6FA" : nullptr), false, nullptr, numeric && c > 0);
            }
            body << "</w:tr>";
        }
        body << "</w:tbl>";
        if (!t.caption.empty()) para(t.caption, "Caption");
        else body << "<w:p/>";
    }
    void chart(const Chart& c) {
        ++nChart;
        const std::string part = "charts/chart" + std::to_string(nChart) + ".xml";
        media.push_back({"word/" + part, chartXml(c)});
        overrides.push_back("<Override PartName=\"/word/" + part + "\" ContentType=\"application/vnd.openxmlformats-officedocument.drawingml.chart+xml\"/>");
        const std::string id = rel("http://schemas.openxmlformats.org/officeDocument/2006/relationships/chart", part);
        const long long cy = 3300000;
        body << "<w:p><w:pPr><w:jc w:val=\"center\"/></w:pPr><w:r><w:drawing><wp:inline distT=\"0\" distB=\"0\" distL=\"0\" distR=\"0\"><wp:extent cx=\"" << TEXT_W << "\" cy=\"" << cy
             << "\"/><wp:effectExtent l=\"0\" t=\"0\" r=\"0\" b=\"0\"/><wp:docPr id=\"" << docPr++ << "\" name=\"Chart " << nChart
             << "\"/><wp:cNvGraphicFramePr/><a:graphic><a:graphicData uri=\"http://schemas.openxmlformats.org/drawingml/2006/chart\"><c:chart r:id=\"" << id
             << "\"/></a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>";
    }
    void image(const Image& im) {
        ++nImg;
        const std::string part = "media/image" + std::to_string(nImg) + ".png";
        media.push_back({"word/" + part, std::string(im.png.begin(), im.png.end())});
        const std::string id = rel("http://schemas.openxmlformats.org/officeDocument/2006/relationships/image", part);
        long long cx = std::min<long long>(TEXT_W * 6 / 10, (long long)im.w * 9525LL * 2);   // 2x pixel size, at most 60 % of the text width
        long long cy = im.w > 0 ? cx * im.h / im.w : 0;
        if (cy > 3600000) { cx = cx * 3600000 / cy; cy = 3600000; }
        body << "<w:p><w:pPr><w:jc w:val=\"center\"/></w:pPr><w:r><w:drawing><wp:inline distT=\"0\" distB=\"0\" distL=\"0\" distR=\"0\"><wp:extent cx=\"" << cx << "\" cy=\"" << cy
             << "\"/><wp:effectExtent l=\"0\" t=\"0\" r=\"0\" b=\"0\"/><wp:docPr id=\"" << docPr << "\" name=\"Picture " << nImg
             << "\"/><wp:cNvGraphicFramePr><a:graphicFrameLocks noChangeAspect=\"1\"/></wp:cNvGraphicFramePr><a:graphic><a:graphicData uri=\"http://schemas.openxmlformats.org/drawingml/2006/picture\">"
             << "<pic:pic><pic:nvPicPr><pic:cNvPr id=\"" << docPr << "\" name=\"image" << nImg << ".png\"/><pic:cNvPicPr/></pic:nvPicPr><pic:blipFill><a:blip r:embed=\"" << id
             << "\"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill><pic:spPr><a:xfrm><a:off x=\"0\" y=\"0\"/><a:ext cx=\"" << cx << "\" cy=\"" << cy
             << "\"/></a:xfrm><a:prstGeom prst=\"rect\"><a:avLst/></a:prstGeom></pic:spPr></pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>";
        ++docPr;
        if (!im.caption.empty()) para(im.caption, "Caption");
    }
};

const char* kDocxStyles =
    "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
    "<w:styles xmlns:w=\"http://schemas.openxmlformats.org/wordprocessingml/2006/main\">"
    "<w:docDefaults><w:rPrDefault><w:rPr><w:rFonts w:ascii=\"Calibri\" w:hAnsi=\"Calibri\" w:eastAsia=\"Calibri\" w:cs=\"Calibri\"/><w:sz w:val=\"21\"/><w:szCs w:val=\"21\"/>"
    "<w:lang w:val=\"en-US\"/></w:rPr></w:rPrDefault><w:pPrDefault><w:pPr><w:spacing w:after=\"120\" w:line=\"264\" w:lineRule=\"auto\"/></w:pPr></w:pPrDefault></w:docDefaults>"
    "<w:style w:type=\"paragraph\" w:default=\"1\" w:styleId=\"Normal\"><w:name w:val=\"Normal\"/><w:qFormat/></w:style>"
    "<w:style w:type=\"paragraph\" w:styleId=\"Title\"><w:name w:val=\"Title\"/><w:basedOn w:val=\"Normal\"/><w:next w:val=\"Normal\"/><w:qFormat/>"
    "<w:pPr><w:spacing w:before=\"1800\" w:after=\"240\"/></w:pPr><w:rPr><w:b/><w:color w:val=\"1F4E79\"/><w:sz w:val=\"52\"/><w:szCs w:val=\"52\"/></w:rPr></w:style>"
    "<w:style w:type=\"paragraph\" w:styleId=\"Subtitle\"><w:name w:val=\"Subtitle\"/><w:basedOn w:val=\"Normal\"/><w:next w:val=\"Normal\"/><w:qFormat/>"
    "<w:rPr><w:color w:val=\"595959\"/><w:sz w:val=\"28\"/><w:szCs w:val=\"28\"/></w:rPr></w:style>"
    "<w:style w:type=\"paragraph\" w:styleId=\"Heading1\"><w:name w:val=\"heading 1\"/><w:basedOn w:val=\"Normal\"/><w:next w:val=\"Normal\"/><w:qFormat/>"
    "<w:pPr><w:keepNext/><w:spacing w:before=\"360\" w:after=\"120\"/><w:pBdr><w:bottom w:val=\"single\" w:sz=\"6\" w:space=\"2\" w:color=\"1F4E79\"/></w:pBdr><w:outlineLvl w:val=\"0\"/></w:pPr>"
    "<w:rPr><w:b/><w:color w:val=\"1F4E79\"/><w:sz w:val=\"32\"/><w:szCs w:val=\"32\"/></w:rPr></w:style>"
    "<w:style w:type=\"paragraph\" w:styleId=\"Heading2\"><w:name w:val=\"heading 2\"/><w:basedOn w:val=\"Normal\"/><w:next w:val=\"Normal\"/><w:qFormat/>"
    "<w:pPr><w:keepNext/><w:spacing w:before=\"240\" w:after=\"80\"/><w:outlineLvl w:val=\"1\"/></w:pPr><w:rPr><w:b/><w:color w:val=\"2E75B6\"/><w:sz w:val=\"26\"/><w:szCs w:val=\"26\"/></w:rPr></w:style>"
    "<w:style w:type=\"paragraph\" w:styleId=\"Heading3\"><w:name w:val=\"heading 3\"/><w:basedOn w:val=\"Normal\"/><w:next w:val=\"Normal\"/><w:qFormat/>"
    "<w:pPr><w:keepNext/><w:spacing w:before=\"200\" w:after=\"60\"/><w:outlineLvl w:val=\"2\"/></w:pPr><w:rPr><w:b/><w:color w:val=\"404040\"/><w:sz w:val=\"22\"/><w:szCs w:val=\"22\"/></w:rPr></w:style>"
    "<w:style w:type=\"paragraph\" w:styleId=\"Caption\"><w:name w:val=\"caption\"/><w:basedOn w:val=\"Normal\"/><w:next w:val=\"Normal\"/><w:qFormat/>"
    "<w:pPr><w:jc w:val=\"center\"/><w:spacing w:after=\"200\"/></w:pPr><w:rPr><w:i/><w:color w:val=\"595959\"/><w:sz w:val=\"18\"/><w:szCs w:val=\"18\"/></w:rPr></w:style>"
    "<w:style w:type=\"paragraph\" w:styleId=\"TableText\"><w:name w:val=\"Table Text\"/><w:basedOn w:val=\"Normal\"/><w:pPr><w:spacing w:before=\"20\" w:after=\"20\" w:line=\"240\" w:lineRule=\"auto\"/></w:pPr>"
    "<w:rPr><w:sz w:val=\"18\"/><w:szCs w:val=\"18\"/></w:rPr></w:style>"
    "<w:style w:type=\"paragraph\" w:styleId=\"ListBullet\"><w:name w:val=\"List Bullet\"/><w:basedOn w:val=\"Normal\"/><w:pPr><w:numPr><w:ilvl w:val=\"0\"/><w:numId w:val=\"1\"/></w:numPr>"
    "<w:spacing w:after=\"40\"/></w:pPr></w:style>"
    "<w:style w:type=\"table\" w:styleId=\"ResTable\"><w:name w:val=\"ResSim Table\"/><w:tblPr><w:tblInd w:w=\"0\" w:type=\"dxa\"/></w:tblPr></w:style>"
    "</w:styles>";
const char* kDocxNumbering =
    "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
    "<w:numbering xmlns:w=\"http://schemas.openxmlformats.org/wordprocessingml/2006/main\"><w:abstractNum w:abstractNumId=\"0\"><w:multiLevelType w:val=\"singleLevel\"/>"
    "<w:lvl w:ilvl=\"0\"><w:start w:val=\"1\"/><w:numFmt w:val=\"bullet\"/><w:lvlText w:val=\"\xE2\x80\xA2\"/><w:lvlJc w:val=\"left\"/><w:pPr><w:ind w:left=\"360\" w:hanging=\"240\"/></w:pPr></w:lvl>"
    "</w:abstractNum><w:num w:numId=\"1\"><w:abstractNumId w:val=\"0\"/></w:num></w:numbering>";
}  // namespace

bool writeDocx(const std::string& path, const Doc& d, std::string& err) {
    DocxBuilder B;
    B.para(d.title.empty() ? "ResSim report" : d.title, "Title");
    if (!d.subtitle.empty()) B.para(d.subtitle, "Subtitle");
    {
        std::string meta = (d.date.empty() ? today() : d.date);
        if (!d.author.empty()) meta += "  |  " + d.author;
        meta += "  |  generated by ResSim Studio";
        B.para(meta, "Subtitle");
    }
    for (const Block& b : d.blocks) {
        switch (b.kind) {
        case Block::Heading: B.para(b.text, b.level <= 1 ? "Heading1" : b.level == 2 ? "Heading2" : "Heading3"); break;
        case Block::Para: B.para(b.text); break;
        case Block::Bullets: for (const auto& it : b.items) B.para(it, "ListBullet"); break;
        case Block::TableB: B.table(b.table); break;
        case Block::ChartB: B.chart(b.chart); break;
        case Block::ImageB: B.image(b.image); break;
        case Block::PageBreak: B.body << "<w:p><w:r><w:br w:type=\"page\"/></w:r></w:p>"; break;
        }
    }
    std::ostringstream doc;
    doc << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
           "<w:document xmlns:w=\"http://schemas.openxmlformats.org/wordprocessingml/2006/main\" xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\" "
           "xmlns:wp=\"http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing\" xmlns:a=\"http://schemas.openxmlformats.org/drawingml/2006/main\" "
           "xmlns:pic=\"http://schemas.openxmlformats.org/drawingml/2006/picture\" xmlns:c=\"http://schemas.openxmlformats.org/drawingml/2006/chart\"><w:body>"
        << B.body.str()
        << "<w:sectPr><w:footerReference w:type=\"default\" r:id=\"rId3\"/><w:pgSz w:w=\"11906\" w:h=\"16838\"/><w:pgMar w:top=\"1134\" w:right=\"1134\" w:bottom=\"1134\" w:left=\"1134\" "
           "w:header=\"567\" w:footer=\"567\" w:gutter=\"0\"/></w:sectPr></w:body></w:document>";
    std::string footer = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n<w:ftr xmlns:w=\"http://schemas.openxmlformats.org/wordprocessingml/2006/main\">"
                         "<w:p><w:pPr><w:jc w:val=\"center\"/></w:pPr><w:r><w:rPr><w:color w:val=\"808080\"/><w:sz w:val=\"16\"/></w:rPr><w:t xml:space=\"preserve\">" +
                         xesc(d.title) + "  -  page </w:t></w:r><w:r><w:rPr><w:color w:val=\"808080\"/><w:sz w:val=\"16\"/></w:rPr><w:fldChar w:fldCharType=\"begin\"/></w:r>"
                         "<w:r><w:rPr><w:color w:val=\"808080\"/><w:sz w:val=\"16\"/></w:rPr><w:instrText xml:space=\"preserve\"> PAGE </w:instrText></w:r>"
                         "<w:r><w:rPr><w:color w:val=\"808080\"/><w:sz w:val=\"16\"/></w:rPr><w:fldChar w:fldCharType=\"separate\"/></w:r><w:r><w:rPr><w:color w:val=\"808080\"/><w:sz w:val=\"16\"/></w:rPr><w:t>1</w:t></w:r>"
                         "<w:r><w:rPr><w:color w:val=\"808080\"/><w:sz w:val=\"16\"/></w:rPr><w:fldChar w:fldCharType=\"end\"/></w:r></w:p></w:ftr>";
    std::ostringstream rels;
    rels << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
            "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles\" Target=\"styles.xml\"/>"
            "<Relationship Id=\"rId2\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/numbering\" Target=\"numbering.xml\"/>"
            "<Relationship Id=\"rId3\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/footer\" Target=\"footer1.xml\"/>";
    for (const auto& r : B.rels) {
        const size_t bar = r.second.find('|');
        rels << "<Relationship Id=\"" << r.first << "\" Type=\"" << r.second.substr(0, bar) << "\" Target=\"" << r.second.substr(bar + 1) << "\"/>";
    }
    rels << "</Relationships>";
    std::ostringstream ct;
    ct << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">"
          "<Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/><Default Extension=\"xml\" ContentType=\"application/xml\"/>"
          "<Default Extension=\"png\" ContentType=\"image/png\"/>"
          "<Override PartName=\"/word/document.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml\"/>"
          "<Override PartName=\"/word/styles.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml\"/>"
          "<Override PartName=\"/word/numbering.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.wordprocessingml.numbering+xml\"/>"
          "<Override PartName=\"/word/footer1.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml\"/>"
          "<Override PartName=\"/docProps/core.xml\" ContentType=\"application/vnd.openxmlformats-package.core-properties+xml\"/>"
          "<Override PartName=\"/docProps/app.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.extended-properties+xml\"/>";
    for (const auto& o : B.overrides) ct << o;
    ct << "</Types>";
    std::vector<ZipEntry> z;
    z.push_back({"[Content_Types].xml", ct.str()});
    z.push_back({"_rels/.rels", "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
                                "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" Target=\"word/document.xml\"/>"
                                "<Relationship Id=\"rId2\" Type=\"http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties\" Target=\"docProps/core.xml\"/>"
                                "<Relationship Id=\"rId3\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties\" Target=\"docProps/app.xml\"/>"
                                "</Relationships>"});
    z.push_back({"docProps/core.xml", corePropsXml(d)});
    z.push_back({"docProps/app.xml", appPropsXml()});
    z.push_back({"word/document.xml", doc.str()});
    z.push_back({"word/styles.xml", kDocxStyles});
    z.push_back({"word/numbering.xml", kDocxNumbering});
    z.push_back({"word/footer1.xml", footer});
    z.push_back({"word/_rels/document.xml.rels", rels.str()});
    for (auto& m : B.media) z.push_back(std::move(m));
    return writeZip(path, z, err);
}

// ================================================================ PPTX
namespace {
const long long SW = 12192000, SH = 6858000;   // 16:9
const char* kTheme =
    "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
    "<a:theme xmlns:a=\"http://schemas.openxmlformats.org/drawingml/2006/main\" name=\"ResSim\"><a:themeElements>"
    "<a:clrScheme name=\"ResSim\"><a:dk1><a:srgbClr val=\"000000\"/></a:dk1><a:lt1><a:srgbClr val=\"FFFFFF\"/></a:lt1><a:dk2><a:srgbClr val=\"1F4E79\"/></a:dk2>"
    "<a:lt2><a:srgbClr val=\"E7E6E6\"/></a:lt2><a:accent1><a:srgbClr val=\"1F77B4\"/></a:accent1><a:accent2><a:srgbClr val=\"D62728\"/></a:accent2>"
    "<a:accent3><a:srgbClr val=\"2CA02C\"/></a:accent3><a:accent4><a:srgbClr val=\"FF7F0E\"/></a:accent4><a:accent5><a:srgbClr val=\"9467BD\"/></a:accent5>"
    "<a:accent6><a:srgbClr val=\"8C564B\"/></a:accent6><a:hlink><a:srgbClr val=\"0563C1\"/></a:hlink><a:folHlink><a:srgbClr val=\"954F72\"/></a:folHlink></a:clrScheme>"
    "<a:fontScheme name=\"ResSim\"><a:majorFont><a:latin typeface=\"Calibri Light\"/><a:ea typeface=\"\"/><a:cs typeface=\"\"/></a:majorFont>"
    "<a:minorFont><a:latin typeface=\"Calibri\"/><a:ea typeface=\"\"/><a:cs typeface=\"\"/></a:minorFont></a:fontScheme>"
    "<a:fmtScheme name=\"ResSim\"><a:fillStyleLst><a:solidFill><a:schemeClr val=\"phClr\"/></a:solidFill><a:solidFill><a:schemeClr val=\"phClr\"/></a:solidFill>"
    "<a:solidFill><a:schemeClr val=\"phClr\"/></a:solidFill></a:fillStyleLst><a:lnStyleLst><a:ln w=\"6350\"><a:solidFill><a:schemeClr val=\"phClr\"/></a:solidFill></a:ln>"
    "<a:ln w=\"12700\"><a:solidFill><a:schemeClr val=\"phClr\"/></a:solidFill></a:ln><a:ln w=\"19050\"><a:solidFill><a:schemeClr val=\"phClr\"/></a:solidFill></a:ln></a:lnStyleLst>"
    "<a:effectStyleLst><a:effectStyle><a:effectLst/></a:effectStyle><a:effectStyle><a:effectLst/></a:effectStyle><a:effectStyle><a:effectLst/></a:effectStyle></a:effectStyleLst>"
    "<a:bgFillStyleLst><a:solidFill><a:schemeClr val=\"phClr\"/></a:solidFill><a:solidFill><a:schemeClr val=\"phClr\"/></a:solidFill>"
    "<a:solidFill><a:schemeClr val=\"phClr\"/></a:solidFill></a:bgFillStyleLst></a:fmtScheme></a:themeElements></a:theme>";
const char* kNs = "xmlns:a=\"http://schemas.openxmlformats.org/drawingml/2006/main\" xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\" "
                  "xmlns:p=\"http://schemas.openxmlformats.org/presentationml/2006/main\"";
const char* kEmptyTree = "<p:nvGrpSpPr><p:cNvPr id=\"1\" name=\"\"/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr><a:xfrm><a:off x=\"0\" y=\"0\"/><a:ext cx=\"0\" cy=\"0\"/>"
                         "<a:chOff x=\"0\" y=\"0\"/><a:chExt cx=\"0\" cy=\"0\"/></a:xfrm></p:grpSpPr>";

struct Slide {
    std::ostringstream xml;
    std::vector<std::pair<std::string, std::string>> rels;   // id -> type|target
    int nextId = 2, nextRel = 2;
    std::string rel(const std::string& type, const std::string& target) {
        const std::string id = "rId" + std::to_string(nextRel++);
        rels.push_back({id, type + "|" + target});
        return id;
    }
    void text(long long x, long long y, long long w, long long h, const std::vector<std::string>& lines, int sz, bool bold, const char* color,
              bool bullets = false, const char* fill = nullptr, bool center = false) {
        xml << "<p:sp><p:nvSpPr><p:cNvPr id=\"" << nextId << "\" name=\"Text " << nextId << "\"/><p:cNvSpPr txBox=\"1\"/><p:nvPr/></p:nvSpPr><p:spPr><a:xfrm><a:off x=\"" << x
            << "\" y=\"" << y << "\"/><a:ext cx=\"" << w << "\" cy=\"" << h << "\"/></a:xfrm><a:prstGeom prst=\"rect\"><a:avLst/></a:prstGeom>";
        if (fill) xml << "<a:solidFill><a:srgbClr val=\"" << fill << "\"/></a:solidFill>";
        else xml << "<a:noFill/>";
        xml << "</p:spPr><p:txBody><a:bodyPr wrap=\"square\" lIns=\"91440\" tIns=\"45720\" rIns=\"91440\" bIns=\"45720\" anchor=\"" << (center ? "ctr" : "t")
            << "\"><a:normAutofit/></a:bodyPr><a:lstStyle/>";
        for (const auto& l : lines) {
            xml << "<a:p>";
            if (bullets) xml << "<a:pPr marL=\"285750\" indent=\"-285750\"><a:buFont typeface=\"Arial\"/><a:buChar char=\"&#8226;\"/></a:pPr>";
            else if (center) xml << "<a:pPr algn=\"ctr\"/>";
            xml << "<a:r><a:rPr lang=\"en-US\" sz=\"" << sz << "\" b=\"" << (bold ? 1 : 0) << "\" dirty=\"0\"><a:solidFill><a:srgbClr val=\"" << (color ? color : "262626")
                << "\"/></a:solidFill></a:rPr><a:t>" << xesc(l) << "</a:t></a:r></a:p>";
        }
        xml << "</p:txBody></p:sp>";
        ++nextId;
    }
    void title(const std::string& t) {
        // band at the top
        xml << "<p:sp><p:nvSpPr><p:cNvPr id=\"" << nextId << "\" name=\"Band\"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr><p:spPr><a:xfrm><a:off x=\"0\" y=\"0\"/><a:ext cx=\"" << SW
            << "\" cy=\"91440\"/></a:xfrm><a:prstGeom prst=\"rect\"><a:avLst/></a:prstGeom><a:solidFill><a:srgbClr val=\"1F4E79\"/></a:solidFill><a:ln><a:noFill/></a:ln></p:spPr></p:sp>";
        ++nextId;
        text(457200, 180000, SW - 914400, 640000, {t}, 2600, true, "1F4E79");
    }
    void table(const Table& t, long long x, long long y, long long w, size_t r0, size_t r1) {
        const size_t nc = std::max<size_t>(1, t.header.size());
        const long long rowH = 300000;
        xml << "<p:graphicFrame><p:nvGraphicFramePr><p:cNvPr id=\"" << nextId++ << "\" name=\"Table\"/><p:cNvGraphicFramePr><a:graphicFrameLocks noGrp=\"1\"/></p:cNvGraphicFramePr><p:nvPr/>"
               "</p:nvGraphicFramePr><p:xfrm><a:off x=\"" << x << "\" y=\"" << y << "\"/><a:ext cx=\"" << w << "\" cy=\"" << rowH * long(r1 - r0 + 1)
            << "\"/></p:xfrm><a:graphic><a:graphicData uri=\"http://schemas.openxmlformats.org/drawingml/2006/table\"><a:tbl><a:tblPr firstRow=\"1\" bandRow=\"1\"/><a:tblGrid>";
        for (size_t c = 0; c < nc; ++c) xml << "<a:gridCol w=\"" << w / long(nc) << "\"/>";
        xml << "</a:tblGrid>";
        const int sz = nc > 7 ? 1000 : 1200;
        auto cell = [&](const std::string& s, const char* fill, bool bold, const char* color) {
            xml << "<a:tc><a:txBody><a:bodyPr/><a:lstStyle/><a:p><a:r><a:rPr lang=\"en-US\" sz=\"" << sz << "\" b=\"" << (bold ? 1 : 0) << "\" dirty=\"0\"><a:solidFill><a:srgbClr val=\""
                << color << "\"/></a:solidFill></a:rPr><a:t>" << xesc(s) << "</a:t></a:r></a:p></a:txBody><a:tcPr marL=\"68580\" marR=\"68580\" marT=\"34290\" marB=\"34290\">"
                << "<a:lnB w=\"6350\"><a:solidFill><a:srgbClr val=\"BFBFBF\"/></a:solidFill></a:lnB><a:solidFill><a:srgbClr val=\"" << fill << "\"/></a:solidFill></a:tcPr></a:tc>";
        };
        xml << "<a:tr h=\"" << rowH << "\">";
        for (size_t c = 0; c < nc; ++c) cell(c < t.header.size() ? t.header[c] : "", "1F4E79", true, "FFFFFF");
        xml << "</a:tr>";
        for (size_t r = r0; r < r1; ++r) {
            const bool hl = std::find(t.highlightRows.begin(), t.highlightRows.end(), int(r)) != t.highlightRows.end();
            xml << "<a:tr h=\"" << rowH << "\">";
            for (size_t c = 0; c < nc; ++c) cell(c < t.rows[r].size() ? t.rows[r][c] : "", hl ? "FCE4D6" : (r % 2 ? "F2F6FA" : "FFFFFF"), false, "262626");
            xml << "</a:tr>";
        }
        xml << "</a:tbl></a:graphicData></a:graphic></p:graphicFrame>";
    }
    void chart(const std::string& rid, long long x, long long y, long long w, long long h) {
        xml << "<p:graphicFrame><p:nvGraphicFramePr><p:cNvPr id=\"" << nextId++ << "\" name=\"Chart\"/><p:cNvGraphicFramePr/><p:nvPr/></p:nvGraphicFramePr><p:xfrm><a:off x=\"" << x
            << "\" y=\"" << y << "\"/><a:ext cx=\"" << w << "\" cy=\"" << h << "\"/></p:xfrm><a:graphic><a:graphicData uri=\"http://schemas.openxmlformats.org/drawingml/2006/chart\">"
            << "<c:chart xmlns:c=\"http://schemas.openxmlformats.org/drawingml/2006/chart\" r:id=\"" << rid << "\"/></a:graphicData></a:graphic></p:graphicFrame>";
    }
    void picture(const std::string& rid, long long x, long long y, long long w, long long h) {
        xml << "<p:pic><p:nvPicPr><p:cNvPr id=\"" << nextId++ << "\" name=\"Picture\"/><p:cNvPicPr><a:picLocks noChangeAspect=\"1\"/></p:cNvPicPr><p:nvPr/></p:nvPicPr>"
            << "<p:blipFill><a:blip r:embed=\"" << rid << "\"/><a:stretch><a:fillRect/></a:stretch></p:blipFill><p:spPr><a:xfrm><a:off x=\"" << x << "\" y=\"" << y
            << "\"/><a:ext cx=\"" << w << "\" cy=\"" << h << "\"/></a:xfrm><a:prstGeom prst=\"rect\"><a:avLst/></a:prstGeom></p:spPr></p:pic>";
    }
};
}  // namespace

bool writePptx(const std::string& path, const Doc& d, std::string& err) {
    std::vector<std::unique_ptr<Slide>> slides;
    std::vector<ZipEntry> parts;
    std::vector<std::string> overrides;
    int nChart = 0, nImg = 0;
    auto newSlide = [&](const std::string& title) -> Slide& {
        slides.push_back(std::make_unique<Slide>());
        Slide& s = *slides.back();
        s.rel("http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout", "../slideLayouts/slideLayout1.xml");
        s.rels.front().first = "rId1";
        s.nextRel = 2;
        if (!title.empty()) s.title(title);
        return s;
    };
    {   // title slide
        Slide& s = newSlide("");
        s.xml << "<p:sp><p:nvSpPr><p:cNvPr id=\"" << s.nextId++ << "\" name=\"Bg\"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr><p:spPr><a:xfrm><a:off x=\"0\" y=\"0\"/><a:ext cx=\"" << SW
              << "\" cy=\"" << SH << "\"/></a:xfrm><a:prstGeom prst=\"rect\"><a:avLst/></a:prstGeom><a:solidFill><a:srgbClr val=\"1F4E79\"/></a:solidFill><a:ln><a:noFill/></a:ln></p:spPr></p:sp>";
        s.text(762000, 2100000, SW - 1524000, 1100000, {d.title.empty() ? "ResSim report" : d.title}, 4000, true, "FFFFFF");
        std::vector<std::string> sub;
        if (!d.subtitle.empty()) sub.push_back(d.subtitle);
        sub.push_back((d.date.empty() ? today() : d.date) + (d.author.empty() ? "" : "  |  " + d.author));
        sub.push_back("generated by ResSim Studio");
        s.text(762000, 3300000, SW - 1524000, 1400000, sub, 1800, false, "DDEBF7");
    }
    const long long X0 = 457200, Y0 = 900000, CW = SW - 914400, CH = SH - Y0 - 380000;
    std::string section = "Overview";
    std::vector<std::string> pending;   // text lines of the current section
    auto flushText = [&]() {
        if (pending.empty()) return;
        // split long text into slides of ~12 lines
        for (size_t a = 0; a < pending.size(); a += 12) {
            Slide& s = newSlide(section);
            std::vector<std::string> L(pending.begin() + long(a), pending.begin() + long(std::min(pending.size(), a + 12)));
            s.text(X0, Y0, CW, CH, L, 1600, false, nullptr, true);
        }
        pending.clear();
    };
    const Chart* heldChart = nullptr;   // two consecutive charts share a slide
    Slide* imgSlide = nullptr;          // two consecutive images share a slide
    auto chartPart = [&](const Chart& c) {
        ++nChart;
        const std::string part = "ppt/charts/chart" + std::to_string(nChart) + ".xml";
        parts.push_back({part, chartXml(c)});
        overrides.push_back("<Override PartName=\"/" + part + "\" ContentType=\"application/vnd.openxmlformats-officedocument.drawingml.chart+xml\"/>");
        return "../charts/chart" + std::to_string(nChart) + ".xml";
    };
    auto flushChart = [&](const Chart* second) {
        if (!heldChart) return;
        Slide& s = newSlide(section);
        const std::string t1 = chartPart(*heldChart);
        const std::string r1 = s.rel("http://schemas.openxmlformats.org/officeDocument/2006/relationships/chart", t1);
        if (second) {
            const std::string t2 = chartPart(*second);
            const std::string r2 = s.rel("http://schemas.openxmlformats.org/officeDocument/2006/relationships/chart", t2);
            s.chart(r1, X0, Y0, CW / 2 - 60000, CH);
            s.chart(r2, X0 + CW / 2 + 60000, Y0, CW / 2 - 60000, CH);
        } else s.chart(r1, X0 + CW / 8, Y0, CW * 3 / 4, CH);
        heldChart = nullptr;
    };
    for (const Block& b : d.blocks) {
        if (b.kind != Block::ChartB) flushChart(nullptr);
        if (b.kind != Block::ImageB) imgSlide = nullptr;
        switch (b.kind) {
        case Block::Heading:
            if (b.level <= 2) { flushText(); section = b.text; }
            else pending.push_back(b.text);
            break;
        case Block::Para: if (!b.text.empty()) pending.push_back(b.text); break;
        case Block::Bullets: for (const auto& it : b.items) pending.push_back(it); break;
        case Block::TableB: {
            std::vector<std::string> intro;
            if (pending.size() <= 2) { intro = pending; pending.clear(); }   // short text goes above the table
            flushText();
            const size_t per = intro.empty() ? 13 : 11;
            for (size_t r0 = 0; r0 < std::max<size_t>(1, b.table.rows.size()); r0 += per) {
                Slide& s = newSlide(section + (b.table.rows.size() > per ? " (" + std::to_string(r0 / per + 1) + ")" : ""));
                long long ty = Y0;
                if (!intro.empty() && r0 == 0) { s.text(X0, Y0, CW, 700000, intro, 1400, false, nullptr); ty = Y0 + 760000; }
                s.table(b.table, X0, ty, CW, r0, std::min(b.table.rows.size(), r0 + per));
                if (!b.table.caption.empty()) s.text(X0, SH - 520000, CW, 360000, {b.table.caption}, 1100, false, "595959");
            }
            break;
        }
        case Block::ChartB:
            flushText();
            if (heldChart) flushChart(&b.chart);
            else heldChart = &b.chart;
            break;
        case Block::ImageB: {
            flushText();
            const bool second = imgSlide != nullptr && imgSlide == slides.back().get();
            Slide& s = second ? *imgSlide : newSlide(section);
            ++nImg;
            const std::string part = "ppt/media/image" + std::to_string(nImg) + ".png";
            parts.push_back({part, std::string(b.image.png.begin(), b.image.png.end())});
            const std::string rid = s.rel("http://schemas.openxmlformats.org/officeDocument/2006/relationships/image", "../media/image" + std::to_string(nImg) + ".png");
            const long long half = CW / 2 - 80000;
            long long h = CH - 520000, w = b.image.h > 0 ? h * b.image.w / b.image.h : h;
            if (w > half) { h = h * half / w; w = half; }
            const long long x = X0 + (second ? CW / 2 + 80000 : 0) + (half - w) / 2;
            s.picture(rid, x, Y0, w, h);
            if (!b.image.caption.empty()) s.text(X0 + (second ? CW / 2 + 80000 : 0), Y0 + h + 40000, half, 480000, {b.image.caption}, 1000, false, "595959", false, nullptr, true);
            imgSlide = second ? nullptr : &s;
            break;
        }
        case Block::PageBreak: flushText(); break;
        }
    }
    flushChart(nullptr);
    flushText();
    // slide parts
    std::ostringstream presRels, sldIds, ct;
    presRels << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
                "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideMaster\" Target=\"slideMasters/slideMaster1.xml\"/>"
                "<Relationship Id=\"rId2\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/theme\" Target=\"theme/theme1.xml\"/>";
    for (size_t i = 0; i < slides.size(); ++i) {
        Slide& s = *slides[i];
        // footer: slide number
        if (i > 0) s.text(SW - 1500000, SH - 380000, 1200000, 300000, {std::to_string(i + 1) + " / " + std::to_string(slides.size())}, 1000, false, "808080");
        std::ostringstream x;
        x << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n<p:sld " << kNs << "><p:cSld><p:spTree>" << kEmptyTree << s.xml.str()
          << "</p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sld>";
        const std::string n = std::to_string(i + 1);
        parts.push_back({"ppt/slides/slide" + n + ".xml", x.str()});
        std::ostringstream r;
        r << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">";
        for (const auto& e : s.rels) {
            const size_t bar = e.second.find('|');
            r << "<Relationship Id=\"" << e.first << "\" Type=\"" << e.second.substr(0, bar) << "\" Target=\"" << e.second.substr(bar + 1) << "\"/>";
        }
        r << "</Relationships>";
        parts.push_back({"ppt/slides/_rels/slide" + n + ".xml.rels", r.str()});
        presRels << "<Relationship Id=\"rId" << 10 + i << "\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/slide\" Target=\"slides/slide" << n << ".xml\"/>";
        sldIds << "<p:sldId id=\"" << 256 + i << "\" r:id=\"rId" << 10 + i << "\"/>";
        overrides.push_back("<Override PartName=\"/ppt/slides/slide" + n + ".xml\" ContentType=\"application/vnd.openxmlformats-officedocument.presentationml.slide+xml\"/>");
    }
    presRels << "</Relationships>";
    std::ostringstream pres;
    pres << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n<p:presentation " << kNs << " saveSubsetFonts=\"1\">"
         << "<p:sldMasterIdLst><p:sldMasterId id=\"2147483648\" r:id=\"rId1\"/></p:sldMasterIdLst><p:sldIdLst>" << sldIds.str() << "</p:sldIdLst>"
         << "<p:sldSz cx=\"" << SW << "\" cy=\"" << SH << "\"/><p:notesSz cx=\"6858000\" cy=\"9144000\"/></p:presentation>";
    std::ostringstream master;
    master << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n<p:sldMaster " << kNs << "><p:cSld><p:bg><p:bgPr><a:solidFill><a:srgbClr val=\"FFFFFF\"/></a:solidFill>"
           << "<a:effectLst/></p:bgPr></p:bg><p:spTree>" << kEmptyTree << "</p:spTree></p:cSld>"
           << "<p:clrMap bg1=\"lt1\" tx1=\"dk1\" bg2=\"lt2\" tx2=\"dk2\" accent1=\"accent1\" accent2=\"accent2\" accent3=\"accent3\" accent4=\"accent4\" accent5=\"accent5\" "
              "accent6=\"accent6\" hlink=\"hlink\" folHlink=\"folHlink\"/><p:sldLayoutIdLst><p:sldLayoutId id=\"2147483649\" r:id=\"rId1\"/></p:sldLayoutIdLst></p:sldMaster>";
    std::ostringstream layout;
    layout << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n<p:sldLayout " << kNs << " type=\"blank\" preserve=\"1\"><p:cSld name=\"Blank\"><p:spTree>" << kEmptyTree
           << "</p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sldLayout>";
    const std::string relHead = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">";
    ct << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">"
          "<Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/><Default Extension=\"xml\" ContentType=\"application/xml\"/>"
          "<Default Extension=\"png\" ContentType=\"image/png\"/>"
          "<Override PartName=\"/ppt/presentation.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.presentationml.presentation.main+xml\"/>"
          "<Override PartName=\"/ppt/slideMasters/slideMaster1.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.presentationml.slideMaster+xml\"/>"
          "<Override PartName=\"/ppt/slideLayouts/slideLayout1.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.presentationml.slideLayout+xml\"/>"
          "<Override PartName=\"/ppt/theme/theme1.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.theme+xml\"/>"
          "<Override PartName=\"/docProps/core.xml\" ContentType=\"application/vnd.openxmlformats-package.core-properties+xml\"/>"
          "<Override PartName=\"/docProps/app.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.extended-properties+xml\"/>";
    for (const auto& o : overrides) ct << o;
    ct << "</Types>";
    std::vector<ZipEntry> z;
    z.push_back({"[Content_Types].xml", ct.str()});
    z.push_back({"_rels/.rels", relHead +
                                    "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" Target=\"ppt/presentation.xml\"/>"
                                    "<Relationship Id=\"rId2\" Type=\"http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties\" Target=\"docProps/core.xml\"/>"
                                    "<Relationship Id=\"rId3\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties\" Target=\"docProps/app.xml\"/>"
                                    "</Relationships>"});
    z.push_back({"docProps/core.xml", corePropsXml(d)});
    z.push_back({"docProps/app.xml", appPropsXml()});
    z.push_back({"ppt/presentation.xml", pres.str()});
    z.push_back({"ppt/_rels/presentation.xml.rels", presRels.str()});
    z.push_back({"ppt/slideMasters/slideMaster1.xml", master.str()});
    z.push_back({"ppt/slideMasters/_rels/slideMaster1.xml.rels", relHead +
                                    "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout\" Target=\"../slideLayouts/slideLayout1.xml\"/>"
                                    "<Relationship Id=\"rId2\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/theme\" Target=\"../theme/theme1.xml\"/></Relationships>"});
    z.push_back({"ppt/slideLayouts/slideLayout1.xml", layout.str()});
    z.push_back({"ppt/slideLayouts/_rels/slideLayout1.xml.rels", relHead +
                                    "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideMaster\" Target=\"../slideMasters/slideMaster1.xml\"/></Relationships>"});
    z.push_back({"ppt/theme/theme1.xml", kTheme});
    for (auto& p : parts) z.push_back(std::move(p));
    return writeZip(path, z, err);
}

// ================================================================ HTML / Markdown
bool writeHtml(const std::string& path, const Doc& d, std::string& err) {
    std::ostringstream o;
    o << "<!DOCTYPE html><html><head><meta charset=\"utf-8\"><title>" << xesc(d.title) << "</title><style>"
         "body{font-family:'Segoe UI',Calibri,Arial,sans-serif;max-width:960px;margin:24px auto;padding:0 18px;color:#222;line-height:1.45}"
         "h1.t{color:#1F4E79;font-size:30px;margin-bottom:4px}.sub{color:#595959}h1{color:#1F4E79;border-bottom:2px solid #1F4E79;font-size:22px;margin-top:34px}"
         "h2{color:#2E75B6;font-size:17px}h3{font-size:14px}table{border-collapse:collapse;width:100%;font-size:12.5px;margin:8px 0}"
         "th{background:#1F4E79;color:#fff;text-align:left;padding:4px 6px}td{border-bottom:1px solid #ddd;padding:3px 6px}tr:nth-child(even) td{background:#F2F6FA}"
         "tr.hl td{background:#FCE4D6}.cap{color:#595959;font-style:italic;text-align:center;font-size:12px;margin-bottom:16px}"
         "figure{margin:10px 0;text-align:center}svg{max-width:100%;height:auto}img{max-width:100%;image-rendering:pixelated}"
         "</style></head><body>";
    o << "<h1 class=\"t\">" << xesc(d.title.empty() ? "ResSim report" : d.title) << "</h1>";
    if (!d.subtitle.empty()) o << "<div class=\"sub\">" << xesc(d.subtitle) << "</div>";
    o << "<div class=\"sub\">" << xesc(d.date.empty() ? today() : d.date) << (d.author.empty() ? "" : " | " + xesc(d.author)) << " | generated by ResSim Studio</div>";
    for (const Block& b : d.blocks) {
        switch (b.kind) {
        case Block::Heading: o << "<h" << std::min(3, std::max(1, b.level)) << ">" << xesc(b.text) << "</h" << std::min(3, std::max(1, b.level)) << ">"; break;
        case Block::Para: o << "<p>" << xesc(b.text) << "</p>"; break;
        case Block::Bullets: o << "<ul>"; for (const auto& it : b.items) o << "<li>" << xesc(it) << "</li>"; o << "</ul>"; break;
        case Block::TableB:
            o << "<table><tr>";
            for (const auto& h : b.table.header) o << "<th>" << xesc(h) << "</th>";
            o << "</tr>";
            for (size_t r = 0; r < b.table.rows.size(); ++r) {
                const bool hl = std::find(b.table.highlightRows.begin(), b.table.highlightRows.end(), int(r)) != b.table.highlightRows.end();
                o << (hl ? "<tr class=\"hl\">" : "<tr>");
                for (const auto& c : b.table.rows[r]) o << "<td>" << xesc(c) << "</td>";
                o << "</tr>";
            }
            o << "</table>";
            if (!b.table.caption.empty()) o << "<div class=\"cap\">" << xesc(b.table.caption) << "</div>";
            break;
        case Block::ChartB: o << "<figure>" << svgChart(b.chart) << "</figure>"; break;
        case Block::ImageB:
            o << "<figure><img width=\"" << b.image.w * 2 << "\" src=\"data:image/png;base64," << base64(b.image.png) << "\" alt=\"" << xesc(b.image.caption) << "\">";
            if (!b.image.caption.empty()) o << "<div class=\"cap\">" << xesc(b.image.caption) << "</div>";
            o << "</figure>";
            break;
        case Block::PageBreak: o << "<hr>"; break;
        }
    }
    o << "</body></html>\n";
    std::ofstream f(path, std::ios::binary);
    if (!f) { err = "cannot write " + path; return false; }
    f << o.str();
    return bool(f);
}

bool writeMarkdown(const std::string& path, const Doc& d, std::string& err) {
    std::ostringstream o;
    o << "# " << (d.title.empty() ? "ResSim report" : d.title) << "\n\n";
    if (!d.subtitle.empty()) o << "*" << d.subtitle << "*\n\n";
    o << (d.date.empty() ? today() : d.date) << (d.author.empty() ? "" : " | " + d.author) << "\n\n";
    auto cellEsc = [](std::string s) { std::string r; for (char c : s) { if (c == '|') r += "\\|"; else if (c == '\n') r += ' '; else r += c; } return r; };
    for (const Block& b : d.blocks) {
        switch (b.kind) {
        case Block::Heading: o << std::string(size_t(std::min(3, b.level) + 1), '#') << ' ' << b.text << "\n\n"; break;
        case Block::Para: o << b.text << "\n\n"; break;
        case Block::Bullets: for (const auto& it : b.items) o << "* " << it << "\n"; o << "\n"; break;
        case Block::TableB:
            o << "|"; for (const auto& h : b.table.header) o << ' ' << cellEsc(h) << " |"; o << "\n|";
            for (size_t i = 0; i < b.table.header.size(); ++i) o << "---|";
            o << "\n";
            for (const auto& r : b.table.rows) { o << "|"; for (const auto& c : r) o << ' ' << cellEsc(c) << " |"; o << "\n"; }
            if (!b.table.caption.empty()) o << "\n*" << b.table.caption << "*\n";
            o << "\n";
            break;
        case Block::ChartB: {
            o << "*Chart: " << b.chart.title << " (";
            for (size_t i = 0; i < b.chart.series.size(); ++i) o << (i ? ", " : "") << b.chart.series[i].name;
            o << ")*\n\n";
            break;
        }
        case Block::ImageB: o << "*Image: " << b.image.caption << "*\n\n"; break;
        case Block::PageBreak: o << "---\n\n"; break;
        }
    }
    std::ofstream f(path, std::ios::binary);
    if (!f) { err = "cannot write " + path; return false; }
    f << o.str();
    return bool(f);
}

bool writeDoc(const std::string& path, const Doc& d, std::string& err) {
    const size_t dot = path.find_last_of('.');
    const std::string ext = dot == std::string::npos ? "" : lower(path.substr(dot + 1));
    if (ext == "docx") return writeDocx(path, d, err);
    if (ext == "pptx") return writePptx(path, d, err);
    if (ext == "html" || ext == "htm") return writeHtml(path, d, err);
    if (ext == "md") return writeMarkdown(path, d, err);
    err = "unknown report format '." + ext + "' (use .docx, .pptx, .html or .md)";
    return false;
}

// ================================================================ simulation report
Options parseOptions(const std::string& s) {
    Options o;
    size_t a = 0;
    while (a < s.size()) {
        size_t e = s.find(';', a);
        if (e == std::string::npos) e = s.size();
        const std::string kv = s.substr(a, e - a);
        const size_t q = kv.find('=');
        if (q != std::string::npos) {
            const std::string k = lower(kv.substr(0, q)), v = kv.substr(q + 1);
            const bool on = !(v == "0" || lower(v) == "no" || lower(v) == "false" || lower(v) == "off");
            if (k == "title") o.title = v;
            else if (k == "author") o.author = v;
            else if (k == "maps") o.maps = on;
            else if (k == "wells") o.wellCharts = on;
            else if (k == "inputs") o.inputs = on;
            else if (k == "log") o.log = on;
            else if (k == "layer") o.mapLayer = std::atoi(v.c_str()) - 1;
            else if (k == "maxwells") o.maxWells = std::max(0, std::atoi(v.c_str()));
            else if (k == "notes") o.notes = v;
        }
        a = e + 1;
    }
    return o;
}

namespace {
struct Units { std::string liq, gas, liqR, gasR, pr, gor, len; };
Units unitsOf(const Case& c) {
    if (c.units == UnitSystem::Field) return {"stb", "Mscf", "stb/d", "Mscf/d", "psia", "Mscf/stb", "ft"};
    return {"sm3", "sm3", "sm3/d", "sm3/d", "bar", "sm3/sm3", "m"};
}
std::vector<double> ser(const ecl::Summary& s, const std::string& key) {
    const int i = s.find(key);
    if (i < 0) return {};
    return s.series(i);
}
bool nonzero(const std::vector<double>& v) { for (double x : v) if (std::isfinite(x) && std::fabs(x) > 0) return true; return false; }
std::string fmt(double v, int prec = 4) {
    if (!std::isfinite(v)) return "-";
    char b[48];
    const double a = std::fabs(v);
    if (a != 0 && (a >= 1e7 || a < 1e-3)) std::snprintf(b, sizeof b, "%.*e", prec - 1, v);
    else if (a >= 1000) std::snprintf(b, sizeof b, "%.0f", v);
    else std::snprintf(b, sizeof b, "%.*g", prec, v);
    return b;
}
struct StatRow { double mn = 1e300, mx = -1e300, sum = 0; long n = 0; };
}  // namespace

Doc simulationReport(const Case& c, const api::RunOutput& run, const ecl::Summary& sum, const Options& opt,
                     const ecl::Summary* ref, const std::string& refLabel) {
    Doc d;
    d.title = opt.title.empty() ? (c.title.empty() ? std::string("Simulation report") : c.title) : opt.title;
    d.subtitle = "Simulation report - " + (run.engine.empty() ? api::engineName(c) : run.engine);
    d.author = opt.author;
    const Units U = unitsOf(c);
    const std::vector<double> T = ser(sum, "TIME");
    auto at = [&](const std::string& k, int i) { const auto v = ser(sum, k); if (v.empty()) return std::nan(""); return v[size_t(i < 0 ? int(v.size()) + i : i)]; };

    // ---- 1 summary
    d.h(1, "1. Summary");
    {
        std::ostringstream p;
        p << (run.ok ? "The run finished" : "The run did NOT finish (" + run.error + ")") << ": " << fmt(T.empty() ? 0 : T.back()) << " days simulated in "
          << fmt(run.seconds, 3) << " s wall time (" << run.stats.timesteps << " time steps, " << run.stats.newton << " Newton iterations, " << run.stats.linear
          << " linear iterations, " << run.stats.cuts << " time-step cuts).";
        d.p(p.str());
        if (!opt.notes.empty()) d.p(opt.notes);
        Table t; t.header = {"Indicator", "Value", "Unit"}; t.caption = "Table 1 - key results";
        const double oip0 = run.initial.oil > 0 ? run.initial.oil : (run.summary.empty() ? 0 : run.summary.front().oip);
        const double fopt = at("FOPT", -1);
        auto add = [&](const std::string& n, double v, const std::string& u) { if (std::isfinite(v)) t.rows.push_back({n, fmt(v, 5), u}); };
        add("Initial oil in place", oip0, U.liq);
        if (run.initial.gasFree + run.initial.gasDissolved > 0) add("Initial gas in place (free + dissolved)", run.initial.gasFree + run.initial.gasDissolved, U.gas);
        add("Cumulative oil (FOPT)", fopt, U.liq);
        if (oip0 > 0 && std::isfinite(fopt)) add("Oil recovery factor", 100.0 * fopt / oip0, "%");
        add("Cumulative water produced (FWPT)", at("FWPT", -1), U.liq);
        add("Cumulative gas produced (FGPT)", at("FGPT", -1), U.gas);
        if (nonzero(ser(sum, "FWIT"))) add("Cumulative water injected (FWIT)", at("FWIT", -1), U.liq);
        if (nonzero(ser(sum, "FGIT"))) add("Cumulative gas injected (FGIT)", at("FGIT", -1), U.gas);
        {
            const auto r = ser(sum, "FOPR");
            if (!r.empty()) {
                size_t k = size_t(std::max_element(r.begin(), r.end()) - r.begin());
                add("Peak oil rate (FOPR)", r[k], U.liqR);
                add("  reached at", T[k], "days");
            }
        }
        add("Final water cut (FWCT)", at("FWCT", -1), "-");
        add("Final gas-oil ratio (FGOR)", at("FGOR", -1), U.gor);
        add("Average pressure, start (FPR)", at("FPR", 0), U.pr);
        add("Average pressure, end (FPR)", at("FPR", -1), U.pr);
        for (const auto& k : run.kpi) add(k.first, k.second, "");
        d.table(t);
    }
    // ---- 2 model
    if (opt.inputs) {
        d.h(1, "2. Model");
        const GridInput& g = c.grid;
        const long ncell = long(g.nx) * g.ny * g.nz;
        long nact = 0;
        for (long i = 0; i < ncell; ++i) if (g.actnum.empty() || (size_t(i) < g.actnum.size() && g.actnum[size_t(i)] != 0)) ++nact;
        int prod = 0, inj = 0;
        for (const auto& w : c.wells) (w.type == WellType::Producer ? prod : inj)++;
        Table t; t.header = {"Item", "Value"}; t.caption = "Table 2 - model description";
        t.rows.push_back({"Title", c.title});
        t.rows.push_back({"Unit system", c.units == UnitSystem::Field ? "FIELD" : "METRIC"});
        t.rows.push_back({"Grid", std::to_string(g.nx) + " x " + std::to_string(g.ny) + " x " + std::to_string(g.nz) + (g.zcorn.empty() ? " (block-centred)" : " (corner point)")});
        t.rows.push_back({"Cells (total / active by ACTNUM)", std::to_string(ncell) + " / " + std::to_string(nact)});
        if (!run.globalOf.empty()) t.rows.push_back({"Active cells in the simulation", std::to_string(run.globalOf.size())});
        t.rows.push_back({"Faults", std::to_string(g.faults.size())});
        t.rows.push_back({"Wells (producers / injectors)", std::to_string(c.wells.size()) + " (" + std::to_string(prod) + " / " + std::to_string(inj) + ")"});
        double tEnd = 0;
        for (const auto& s : c.schedule) tEnd += s.length;
        t.rows.push_back({"Report steps / schedule length", std::to_string(c.schedule.size()) + " / " + fmt(tEnd) + " days"});
        char st[32]; std::snprintf(st, sizeof st, "%04d-%02d-%02d", c.startYear, c.startMonth, c.startDay);
        t.rows.push_back({"Start date", st});
        t.rows.push_back({"Engine", run.engine.empty() ? api::engineName(c) : run.engine});
        d.table(t);
        // static property statistics over active cells
        Table s; s.header = {"Property", "Min", "Mean", "Max", "Unit"}; s.caption = "Table 3 - static properties over the active cells";
        auto stat = [&](const char* name, const std::vector<double>& v, const char* unit) {
            if (v.size() != size_t(ncell)) return;
            StatRow r;
            for (long i = 0; i < ncell; ++i) {
                if (!g.actnum.empty() && g.actnum[size_t(i)] == 0) continue;
                const double x = v[size_t(i)];
                if (!std::isfinite(x)) continue;
                r.mn = std::min(r.mn, x); r.mx = std::max(r.mx, x); r.sum += x; ++r.n;
            }
            if (r.n) s.rows.push_back({name, fmt(r.mn), fmt(r.sum / double(r.n)), fmt(r.mx), unit});
        };
        stat("PORO", g.poro, "-"); stat("NTG", g.ntg, "-"); stat("PERMX", g.permx, "mD"); stat("PERMY", g.permy, "mD"); stat("PERMZ", g.permz, "mD");
        stat("DZ", g.dz, U.len.c_str());
        if (!s.rows.empty()) d.table(s);
    }
    // ---- 3 field performance
    d.h(1, "3. Field performance");
    auto addSeries = [&](Chart& ch, const std::string& key, const std::string& label, int color, bool y2) {
        const auto v = ser(sum, key);
        if (v.empty()) return;
        Series s; s.name = label; s.x = T; s.y = v; s.rgb = paletteColor(color); s.y2 = y2;
        ch.series.push_back(s);
        if (ref) {
            const int ri = ref->find(key);
            if (ri >= 0) {
                Series r; r.name = label + " (" + refLabel + ")"; r.x = ref->time; r.y = ref->series(ri); r.rgb = paletteColor(color); r.dashed = true; r.y2 = y2;
                ch.series.push_back(r);
            }
        }
    };
    {
        Chart ch; ch.title = "Production rates"; ch.xLabel = "Time (days)"; ch.yLabel = "Liquid rate (" + U.liqR + ")"; ch.y2Label = "Gas rate (" + U.gasR + ")";
        addSeries(ch, "FOPR", "FOPR", 2, false); addSeries(ch, "FWPR", "FWPR", 0, false);
        if (nonzero(ser(sum, "FGPR"))) addSeries(ch, "FGPR", "FGPR", 1, true);
        d.chart(ch);
        Chart cu; cu.title = "Cumulative production"; cu.xLabel = "Time (days)"; cu.yLabel = "Liquid (" + U.liq + ")"; cu.y2Label = "Gas (" + U.gas + ")";
        addSeries(cu, "FOPT", "FOPT", 2, false); addSeries(cu, "FWPT", "FWPT", 0, false);
        if (nonzero(ser(sum, "FGPT"))) addSeries(cu, "FGPT", "FGPT", 1, true);
        d.chart(cu);
        Chart pr; pr.title = "Average reservoir pressure"; pr.xLabel = "Time (days)"; pr.yLabel = "FPR (" + U.pr + ")";
        addSeries(pr, "FPR", "FPR", 4, false);
        d.chart(pr);
        Chart wc; wc.title = "Water cut and gas-oil ratio"; wc.xLabel = "Time (days)"; wc.yLabel = "FWCT (-)"; wc.y2Label = "FGOR (" + U.gor + ")";
        addSeries(wc, "FWCT", "FWCT", 0, false);
        if (nonzero(ser(sum, "FGOR"))) addSeries(wc, "FGOR", "FGOR", 1, true);
        d.chart(wc);
        if (nonzero(ser(sum, "FWIR")) || nonzero(ser(sum, "FGIR"))) {
            Chart ij; ij.title = "Injection"; ij.xLabel = "Time (days)"; ij.yLabel = "Water (" + U.liqR + ")"; ij.y2Label = "Gas (" + U.gasR + ")";
            if (nonzero(ser(sum, "FWIR"))) addSeries(ij, "FWIR", "FWIR", 0, false);
            if (nonzero(ser(sum, "FGIR"))) addSeries(ij, "FGIR", "FGIR", 1, true);
            d.chart(ij);
        }
    }
    // ---- 4 wells
    if (!run.wellNames.empty()) {
        d.h(1, "4. Wells");
        Table t; t.header = {"Well", "Type", "Oil (" + U.liq + ")", "Water (" + U.liq + ")", "Gas (" + U.gas + ")", "Water inj. (" + U.liq + ")", "Gas inj. (" + U.gas + ")", "Final BHP (" + U.pr + ")"};
        t.caption = "Table 4 - well cumulatives at the end of the run";
        std::vector<std::pair<double, std::string>> rank;
        for (const auto& wn : run.wellNames) {
            std::string type = "-";
            for (const auto& w : c.wells) if (w.name == wn) type = w.type == WellType::Producer ? "producer" : "injector";
            const double opt_ = at("WOPT:" + wn, -1), wpt = at("WWPT:" + wn, -1), gpt = at("WGPT:" + wn, -1), wit = at("WWIT:" + wn, -1), git = at("WGIT:" + wn, -1);
            double bhp = std::nan("");
            const auto b = ser(sum, "WBHP:" + wn);
            for (size_t i = b.size(); i-- > 0;) if (b[i] > 0) { bhp = b[i]; break; }
            t.rows.push_back({wn, type, fmt(opt_), fmt(wpt), fmt(gpt), fmt(wit), fmt(git), fmt(bhp)});
            rank.push_back({(std::isfinite(opt_) ? opt_ : 0) + (std::isfinite(wit) ? wit : 0) * 1e-3, wn});
        }
        d.table(t);
        if (opt.wellCharts && opt.maxWells > 0) {
            std::sort(rank.begin(), rank.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
            for (size_t k = 0; k < rank.size() && int(k) < opt.maxWells; ++k) {
                const std::string& wn = rank[k].second;
                Chart ch; ch.title = "Well " + wn; ch.xLabel = "Time (days)"; ch.yLabel = "Rate (" + U.liqR + ")"; ch.y2Label = "BHP (" + U.pr + ")";
                if (nonzero(ser(sum, "WOPR:" + wn))) addSeries(ch, "WOPR:" + wn, "WOPR", 2, false);
                if (nonzero(ser(sum, "WWPR:" + wn))) addSeries(ch, "WWPR:" + wn, "WWPR", 0, false);
                if (nonzero(ser(sum, "WWIR:" + wn))) addSeries(ch, "WWIR:" + wn, "WWIR", 7, false);
                if (nonzero(ser(sum, "WGIR:" + wn))) addSeries(ch, "WGIR:" + wn, "WGIR", 1, false);
                addSeries(ch, "WBHP:" + wn, "WBHP", 4, true);
                // BHP of shut periods is written as 0: hide
                for (auto& s : ch.series) if (s.y2) for (auto& y : s.y) if (y <= 0) y = std::nan("");
                d.chart(ch);
            }
        }
    }
    // ---- 5 maps
    if (opt.maps && run.snaps.size() >= 1 && !run.globalOf.empty()) {
        const int nx = c.grid.nx, ny = c.grid.ny, nz = c.grid.nz;
        int layer = opt.mapLayer;
        if (layer < 0 || layer >= nz) {   // the layer with the largest initial oil pore volume
            std::vector<double> lo(size_t(nz), 0.0);
            const auto& so = run.snaps.front().so;
            for (size_t a = 0; a < run.globalOf.size() && a < run.porv.size(); ++a) {
                const int g = run.globalOf[a];
                if (g >= nx * ny * nz) continue;
                lo[size_t(g / (nx * ny))] += run.porv[a] * (so.size() > a ? so[a] : 1.0);
            }
            layer = int(std::max_element(lo.begin(), lo.end()) - lo.begin());
        }
        d.h(1, "5. Maps (layer K = " + std::to_string(layer + 1) + ")");
        d.p("Top view of layer " + std::to_string(layer + 1) + " (I to the right, J downwards) at the start and the end of the run; both maps of a property share one colour range.");
        const int scale = std::max(2, std::min(12, 480 / std::max(1, std::max(nx, ny))));
        const Snapshot& s0 = run.snaps.front();
        const Snapshot& s1 = run.snaps.back();
        for (const char* prop : {"PRESSURE", "SOIL", "SWAT", "SGAS"}) {
            const auto& a = s0.field(prop);
            const auto& b = s1.field(prop);
            if (a.empty() || b.empty()) continue;
            double lo = 1e300, hi = -1e300;
            for (float v : a) if (std::isfinite(v)) { lo = std::min(lo, double(v)); hi = std::max(hi, double(v)); }
            for (float v : b) if (std::isfinite(v)) { lo = std::min(lo, double(v)); hi = std::max(hi, double(v)); }
            if (!(hi > lo)) { if (hi - lo == 0 && std::string(prop) != "PRESSURE" && hi == 0) continue; hi = lo + 1e-6; }
            d.h(2, std::string(prop));
            d.image(layerMap(a, run.globalOf, nx, ny, layer, lo, hi, scale, std::string(prop) + " at t = " + fmt(s0.time) + " d"));
            d.image(layerMap(b, run.globalOf, nx, ny, layer, lo, hi, scale, std::string(prop) + " at t = " + fmt(s1.time) + " d"));
        }
    }
    // ---- 6 comparison
    if (ref) {
        d.h(1, "6. Comparison with " + refLabel);
        const auto diffs = ecl::compareSummaries(*ref, sum);
        Table t; t.header = {"Vector", refLabel + " final", "ResSim final", "Difference (%)", "rel. RMS (%)"};
        t.caption = "Deviation over the " + refLabel + " report times (rel. RMS = RMS difference / max |" + refLabel + "|)";
        int row = 0;
        for (const auto& v : diffs) {
            if (v.key == "TIME" || v.key.rfind("NEWTON", 0) == 0 || v.key.rfind("NLINEARS", 0) == 0) continue;
            if (std::fabs(v.refFinal) < 1e-12 && std::fabs(v.cmpFinal) < 1e-12) continue;
            const double dp = std::fabs(v.refFinal) > 0 ? 100.0 * (v.cmpFinal - v.refFinal) / std::fabs(v.refFinal) : std::nan("");
            t.rows.push_back({v.key, fmt(v.refFinal, 6), fmt(v.cmpFinal, 6), fmt(dp, 3), fmt(100.0 * v.relRms, 3)});
            if (v.relRms > 0.05) t.highlightRows.push_back(row);
            ++row;
            if (row >= 60) break;
        }
        d.table(t);
    }
    // ---- 7 log
    if (opt.log && !run.log.empty()) {
        std::vector<std::string> L;
        for (const auto& l : run.log) {
            if (l.find("Warning") != std::string::npos || l.find("ACTIONX") != std::string::npos || l.find("script") != std::string::npos ||
                l.find("WARNING") != std::string::npos || l.find("economic") != std::string::npos)
                L.push_back(l);
            if (L.size() >= 40) break;
        }
        if (!L.empty()) { d.h(1, "Run messages"); d.bullets(L); }
    }
    return d;
}

}  // namespace report
}  // namespace ressim
