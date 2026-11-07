// C39 portable PNG export of fig::Figure (Linux / macOS) through plutovg, the counterpart of figure_png_win.cpp
// (Direct2D + WIC). The layout uses the Helvetica / Times AFM metrics of figure.cpp; glyphs come from a metric-compatible
// system font (Liberation Sans / Serif, Arial, Times New Roman) when present, otherwise the embedded Inter. The PNG
// carries the export DPI in a pHYs chunk (as the WIC encoder writes it on Windows).
#include "figure_png.h"
#ifndef _WIN32
#include "resources.h"
#include <plutovg.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace fig {
namespace {

bool fileExists(const char* p) { FILE* f = std::fopen(p, "rb"); if (!f) return false; std::fclose(f); return true; }

plutovg_font_face_t* loadFace(int font) {
    static const char* sans[] = {
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", "/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf",
        "/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf", "/usr/share/fonts/liberation/LiberationSans-Regular.ttf",
        "/usr/share/fonts/TTF/LiberationSans-Regular.ttf", "/System/Library/Fonts/Supplemental/Arial.ttf", "/Library/Fonts/Arial.ttf",
        "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", nullptr};
    static const char* serif[] = {
        "/usr/share/fonts/truetype/liberation/LiberationSerif-Regular.ttf", "/usr/share/fonts/truetype/liberation2/LiberationSerif-Regular.ttf",
        "/usr/share/fonts/liberation-serif/LiberationSerif-Regular.ttf", "/usr/share/fonts/liberation/LiberationSerif-Regular.ttf",
        "/usr/share/fonts/TTF/LiberationSerif-Regular.ttf", "/System/Library/Fonts/Supplemental/Times New Roman.ttf", "/Library/Fonts/Times New Roman.ttf",
        "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf", "/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf", nullptr};
    for (const char** p = font == 1 ? serif : sans; *p; ++p)
        if (fileExists(*p)) if (plutovg_font_face_t* f = plutovg_font_face_load_from_file(*p, 0)) return f;
    const std::string_view d = appResource(RES_FONT_UI);
    if (!d.empty()) return plutovg_font_face_load_from_data(d.data(), unsigned(d.size()), 0, nullptr, nullptr);
    return nullptr;
}

class VGCanvas : public Canvas {
public:
    plutovg_canvas_t* c = nullptr;
    plutovg_font_face_t* face = nullptr;
    float ptScale = 96.f / 72.f;

    void rgba(uint32_t col) { plutovg_canvas_set_rgba(c, (col & 255) / 255.f, ((col >> 8) & 255) / 255.f, ((col >> 16) & 255) / 255.f, ((col >> 24) & 255) / 255.f); }
    void strokeStyle(double w, int dash) {
        plutovg_canvas_set_line_width(c, float(w));
        plutovg_canvas_set_line_cap(c, dash ? PLUTOVG_LINE_CAP_BUTT : PLUTOVG_LINE_CAP_ROUND);
        plutovg_canvas_set_line_join(c, PLUTOVG_LINE_JOIN_ROUND);
        // dash lengths in units of the stroke width (the SVG / PDF / Direct2D patterns)
        float d[4]; int n = 0;
        if (dash == 1) { d[0] = 4.f; d[1] = 2.2f; n = 2; }
        else if (dash == 2) { d[0] = 0.8f; d[1] = 1.6f; n = 2; }
        else if (dash == 3) { d[0] = 5.f; d[1] = 1.8f; d[2] = 0.8f; d[3] = 1.8f; n = 4; }
        for (int k = 0; k < n; ++k) d[k] *= float(std::max(w, 0.1));
        plutovg_canvas_set_dash(c, 0.f, n ? d : nullptr, n);
    }
    void begin(double, double, bool transparent) override {
        plutovg_canvas_scale(c, ptScale, ptScale);
        if (!transparent) { plutovg_canvas_set_rgba(c, 1, 1, 1, 1); plutovg_canvas_paint(c); }
    }
    void end() override {}
    void polyline(const double* xy, int n, uint32_t col, double w, int dash) override {
        if (n < 2) return;
        plutovg_canvas_new_path(c);
        plutovg_canvas_move_to(c, float(xy[0]), float(xy[1]));
        for (int i = 1; i < n; ++i) plutovg_canvas_line_to(c, float(xy[2 * i]), float(xy[2 * i + 1]));
        rgba(col); strokeStyle(w, dash);
        plutovg_canvas_stroke(c);
    }
    void fillRect(double x, double y, double w, double h, uint32_t col) override {
        rgba(col);
        plutovg_canvas_fill_rect(c, float(x), float(y), float(w), float(h));
    }
    void strokeRect(double x, double y, double w, double h, uint32_t col, double sw) override {
        rgba(col); strokeStyle(sw, 0);
        plutovg_canvas_stroke_rect(c, float(x), float(y), float(w), float(h));
    }
    void circle(double x, double y, double r, uint32_t fill, uint32_t sc, double sw) override {
        plutovg_canvas_new_path(c);
        plutovg_canvas_circle(c, float(x), float(y), float(r));
        rgba(fill); plutovg_canvas_fill_preserve(c);
        rgba(sc); strokeStyle(sw, 0);
        plutovg_canvas_stroke(c);
    }
    void run(const std::string& s, double size, double x, double y) {
        if (!face || s.empty()) return;
        plutovg_canvas_set_font(c, face, float(size));
        plutovg_canvas_fill_text(c, s.data(), int(s.size()), PLUTOVG_TEXT_ENCODING_UTF8, float(x), float(y));
    }
    void text(double x, double y, const std::string& s, double size, int ha, bool vertical, uint32_t col) override {
        rgba(col);
        const double w = textWidth(s, size);
        const double d = ha == 1 ? -w * 0.5 : ha == 2 ? -w : 0;
        plutovg_canvas_save(c);
        plutovg_canvas_translate(c, float(x), float(y));
        if (vertical) plutovg_canvas_rotate(c, float(-M_PI / 2));
        double cx = d;
        size_t i = 0;
        std::string cur;
        auto flush = [&]() { if (!cur.empty()) { run(cur, size, cx, 0); cx += textWidth(cur, size); cur.clear(); } };
        while (i < s.size()) {
            if (s[i] == '^' && i + 1 < s.size() && s[i + 1] == '{') {
                const size_t e = s.find('}', i + 2);
                if (e != std::string::npos) {
                    flush();
                    const std::string sup = s.substr(i + 2, e - i - 2);
                    run(sup, size * 0.68, cx, -size * 0.38);
                    cx += textWidth(sup, size * 0.68);
                    i = e + 1;
                    continue;
                }
            }
            cur += s[i++];
        }
        flush();
        plutovg_canvas_restore(c);
    }
    void fillTriangles(const double* xy, int nTri, uint32_t col) override {
        if (nTri <= 0 || ((col >> 24) & 255) == 0) return;
        plutovg_canvas_new_path(c);
        for (int t = 0; t < nTri; ++t) {
            const double* q = xy + 6 * t;
            plutovg_canvas_move_to(c, float(q[0]), float(q[1]));
            plutovg_canvas_line_to(c, float(q[2]), float(q[3]));
            plutovg_canvas_line_to(c, float(q[4]), float(q[5]));
            plutovg_canvas_close_path(c);
        }
        plutovg_canvas_set_fill_rule(c, PLUTOVG_FILL_RULE_NON_ZERO);
        rgba(col);
        plutovg_canvas_fill(c);
    }
    void textFit(double x, double y, const std::string& s, double size, double widthPt, bool vertical, uint32_t col) override {
        if (s.empty()) return;
        rgba(col);
        const double w = textWidth(s, size);
        const float sx = (w > 0 && widthPt > 0) ? float(std::max(0.4, std::min(2.5, widthPt / w))) : 1.f;
        plutovg_canvas_save(c);
        plutovg_canvas_translate(c, float(x), float(y));
        if (vertical) plutovg_canvas_rotate(c, float(-M_PI / 2));
        plutovg_canvas_scale(c, sx, 1.f);
        run(s, size, 0, 0);
        plutovg_canvas_restore(c);
    }
    void pushClip(double x, double y, double w, double h) override {
        plutovg_canvas_save(c);
        plutovg_canvas_clip_rect(c, float(x), float(y), float(w), float(h));
    }
    void popClip() override { plutovg_canvas_restore(c); }
};

uint32_t crc32(const unsigned char* p, size_t n, uint32_t crc = 0xffffffffu) {
    for (size_t i = 0; i < n; ++i) {
        crc ^= p[i];
        for (int k = 0; k < 8; ++k) crc = (crc >> 1) ^ (0xedb88320u & (0u - (crc & 1u)));
    }
    return crc;
}
void put32(std::vector<unsigned char>& v, uint32_t x) { for (int s = 24; s >= 0; s -= 8) v.push_back((unsigned char)(x >> s)); }
void writeToVec(void* closure, void* data, int size) {
    auto* v = static_cast<std::vector<unsigned char>*>(closure);
    v->insert(v->end(), static_cast<unsigned char*>(data), static_cast<unsigned char*>(data) + size);
}

} // namespace

// PNG bytes of a surface with a pHYs (pixels per metre) chunk inserted after IHDR
bool writeSurfacePNG(plutovg_surface_t* s, const std::string& path, double dpi, std::string* err) {
    std::vector<unsigned char> png;
    if (!plutovg_surface_write_to_png_stream(s, writeToVec, &png) || png.size() < 33) { if (err) *err = "PNG encoding failed"; return false; }
    if (dpi > 0) {
        const uint32_t ppm = uint32_t(std::lround(dpi / 0.0254));
        std::vector<unsigned char> ch;
        put32(ch, 9);
        const size_t body = ch.size();
        ch.insert(ch.end(), {'p', 'H', 'Y', 's'});
        put32(ch, ppm); put32(ch, ppm); ch.push_back(1);
        put32(ch, crc32(ch.data() + body, ch.size() - body) ^ 0xffffffffu);
        png.insert(png.begin() + 33, ch.begin(), ch.end());   // 8 signature + 25 IHDR
    }
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) { if (err) *err = "cannot open " + path; return false; }
    const bool ok = std::fwrite(png.data(), 1, png.size(), f) == png.size();
    std::fclose(f);
    if (!ok && err) *err = "write failed: " + path;
    return ok;
}

bool writePNG(const Figure& fg, const Style& st, const std::string& path, std::string* err) {
    const double dpi = std::max(72, std::min(st.dpi, 1200));
    const int W = int(std::lround(st.widthMM / 25.4 * dpi)), H = int(std::lround(st.heightMM / 25.4 * dpi));
    if (W <= 0 || H <= 0 || double(W) * H > 400e6) { if (err) *err = "invalid image size"; return false; }
    plutovg_surface_t* surf = plutovg_surface_create(W, H);
    if (!surf) { if (err) *err = "out of memory"; return false; }
    VGCanvas cv;
    cv.c = plutovg_canvas_create(surf);
    cv.ptScale = float(dpi / 72.0);
    cv.face = loadFace(st.font);
    cv.font = st.font;
    render(fg, st, cv);
    plutovg_canvas_destroy(cv.c);
    if (cv.face) plutovg_font_face_destroy(cv.face);
    const bool ok = writeSurfacePNG(surf, path, dpi, err);
    plutovg_surface_destroy(surf);
    return ok;
}

} // namespace fig
#endif  // !_WIN32
