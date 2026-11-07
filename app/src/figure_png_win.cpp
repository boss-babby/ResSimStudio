// PNG export of fig::Figure through Direct2D (WIC bitmap render target) + DirectWrite, with the DPI written to
// the PNG metadata. Arial / Times New Roman are metric compatible with the Helvetica / Times layout metrics.
#include "figure_png.h"
#include <windows.h>
#include <d2d1.h>
#include <dwrite.h>
#include <wincodec.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

namespace fig {
namespace {

template <class T> void release(T*& p) { if (p) { p->Release(); p = nullptr; } }
std::wstring widen(const std::string& s) {
    if (s.empty()) return std::wstring();
    const int n = MultiByteToWideChar(CP_UTF8, 0, s.data(), int(s.size()), nullptr, 0);
    std::wstring w(size_t(n), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, s.data(), int(s.size()), &w[0], n);
    return w;
}
D2D1_COLOR_F col(uint32_t c) {
    return D2D1::ColorF((c & 255) / 255.f, ((c >> 8) & 255) / 255.f, ((c >> 16) & 255) / 255.f, ((c >> 24) & 255) / 255.f);
}

class D2DCanvas : public Canvas {
public:
    ID2D1Factory* f = nullptr;
    IDWriteFactory* dw = nullptr;
    ID2D1RenderTarget* rt = nullptr;
    ID2D1SolidColorBrush* br = nullptr;
    std::map<int, ID2D1StrokeStyle*> strokes;
    std::map<int, IDWriteTextFormat*> formats;
    float ptScale = 96.f / 72.f;

    ~D2DCanvas() override {
        for (auto& kv : strokes) release(kv.second);
        for (auto& kv : formats) release(kv.second);
        release(br);
    }
    ID2D1StrokeStyle* stroke(int dash) {
        auto it = strokes.find(dash);
        if (it != strokes.end()) return it->second;
        D2D1_STROKE_STYLE_PROPERTIES p = D2D1::StrokeStyleProperties(D2D1_CAP_STYLE_ROUND, D2D1_CAP_STYLE_ROUND, D2D1_CAP_STYLE_ROUND,
                                                                     D2D1_LINE_JOIN_ROUND, 10.f, D2D1_DASH_STYLE_SOLID, 0.f);
        // dash lengths in units of the stroke width (matching the SVG / PDF patterns of figure.cpp)
        std::vector<float> d;
        if (dash == 1) d = {4.f, 2.2f};
        else if (dash == 2) d = {0.8f, 1.6f};
        else if (dash == 3) d = {5.f, 1.8f, 0.8f, 1.8f};
        if (!d.empty()) { p.dashStyle = D2D1_DASH_STYLE_CUSTOM; p.dashCap = D2D1_CAP_STYLE_FLAT; }
        ID2D1StrokeStyle* s = nullptr;
        f->CreateStrokeStyle(p, d.empty() ? nullptr : d.data(), UINT32(d.size()), &s);
        strokes[dash] = s;
        return s;
    }
    IDWriteTextFormat* format(double sizePt) {
        const int key = int(std::lround(sizePt * 100));
        auto it = formats.find(key);
        if (it != formats.end()) return it->second;
        IDWriteTextFormat* tf = nullptr;
        dw->CreateTextFormat(font == 1 ? L"Times New Roman" : L"Arial", nullptr, DWRITE_FONT_WEIGHT_NORMAL, DWRITE_FONT_STYLE_NORMAL,
                             DWRITE_FONT_STRETCH_NORMAL, float(sizePt), L"en-us", &tf);
        if (tf) tf->SetWordWrapping(DWRITE_WORD_WRAPPING_NO_WRAP);
        formats[key] = tf;
        return tf;
    }
    void begin(double, double, bool transparent) override {
        rt->BeginDraw();
        rt->SetTransform(D2D1::Matrix3x2F::Scale(ptScale, ptScale));
        rt->SetAntialiasMode(D2D1_ANTIALIAS_MODE_PER_PRIMITIVE);
        rt->SetTextAntialiasMode(D2D1_TEXT_ANTIALIAS_MODE_GRAYSCALE);
        rt->Clear(transparent ? D2D1::ColorF(1, 1, 1, 0) : D2D1::ColorF(1, 1, 1, 1));
        rt->CreateSolidColorBrush(D2D1::ColorF(0, 0, 0, 1), &br);
    }
    void end() override {}
    void polyline(const double* xy, int n, uint32_t c, double w, int dash) override {
        if (n < 2) return;
        br->SetColor(col(c));
        ID2D1PathGeometry* g = nullptr;
        if (FAILED(f->CreatePathGeometry(&g))) return;
        ID2D1GeometrySink* s = nullptr;
        g->Open(&s);
        s->BeginFigure(D2D1::Point2F(float(xy[0]), float(xy[1])), D2D1_FIGURE_BEGIN_HOLLOW);
        for (int i = 1; i < n; ++i) s->AddLine(D2D1::Point2F(float(xy[2 * i]), float(xy[2 * i + 1])));
        s->EndFigure(D2D1_FIGURE_END_OPEN);
        s->Close();
        release(s);
        rt->DrawGeometry(g, br, float(w), stroke(dash));
        release(g);
    }
    void fillRect(double x, double y, double w, double h, uint32_t c) override {
        br->SetColor(col(c));
        rt->FillRectangle(D2D1::RectF(float(x), float(y), float(x + w), float(y + h)), br);
    }
    void strokeRect(double x, double y, double w, double h, uint32_t c, double sw) override {
        br->SetColor(col(c));
        rt->DrawRectangle(D2D1::RectF(float(x), float(y), float(x + w), float(y + h)), br, float(sw), stroke(0));
    }
    void circle(double x, double y, double r, uint32_t fill, uint32_t sc, double sw) override {
        const D2D1_ELLIPSE e = D2D1::Ellipse(D2D1::Point2F(float(x), float(y)), float(r), float(r));
        br->SetColor(col(fill));
        rt->FillEllipse(e, br);
        br->SetColor(col(sc));
        rt->DrawEllipse(e, br, float(sw));
    }
    // draws one run with its baseline at (x, y) in the current transform
    void run(const std::string& s, double size, double x, double y) {
        IDWriteTextFormat* tf = format(size);
        if (!tf) return;
        const std::wstring w = widen(s);
        IDWriteTextLayout* lay = nullptr;
        if (FAILED(dw->CreateTextLayout(w.c_str(), UINT32(w.size()), tf, 10000.f, 1000.f, &lay))) return;
        DWRITE_LINE_METRICS lm{};
        UINT32 cnt = 0;
        lay->GetLineMetrics(&lm, 1, &cnt);
        rt->DrawTextLayout(D2D1::Point2F(float(x), float(y - lm.baseline)), lay, br, D2D1_DRAW_TEXT_OPTIONS_NONE);
        release(lay);
    }
    void text(double x, double y, const std::string& s, double size, int ha, bool vertical, uint32_t c) override {
        br->SetColor(col(c));
        const double w = textWidth(s, size);
        const double d = ha == 1 ? -w * 0.5 : ha == 2 ? -w : 0;
        D2D1_MATRIX_3X2_F old;
        rt->GetTransform(&old);
        if (vertical) rt->SetTransform(D2D1::Matrix3x2F::Rotation(-90.f) * D2D1::Matrix3x2F::Translation(float(x), float(y)) * old);
        else rt->SetTransform(D2D1::Matrix3x2F::Translation(float(x), float(y)) * old);
        // split "^{...}" superscripts
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
        rt->SetTransform(old);
    }
    void fillTriangles(const double* xy, int nTri, uint32_t c) override {
        if (nTri <= 0 || ((c >> 24) & 255) == 0) return;
        ID2D1PathGeometry* g = nullptr;
        if (FAILED(f->CreatePathGeometry(&g)) || !g) return;
        ID2D1GeometrySink* sk = nullptr;
        if (SUCCEEDED(g->Open(&sk)) && sk) {
            sk->SetFillMode(D2D1_FILL_MODE_WINDING);
            for (int t = 0; t < nTri; ++t) {
                const double* q = xy + 6 * t;
                sk->BeginFigure(D2D1::Point2F(float(q[0]), float(q[1])), D2D1_FIGURE_BEGIN_FILLED);
                sk->AddLine(D2D1::Point2F(float(q[2]), float(q[3])));
                sk->AddLine(D2D1::Point2F(float(q[4]), float(q[5])));
                sk->EndFigure(D2D1_FIGURE_END_CLOSED);
            }
            sk->Close();
            release(sk);
            br->SetColor(col(c));
            rt->FillGeometry(g, br);
        }
        release(g);
    }
    void textFit(double x, double y, const std::string& s, double size, double widthPt, bool vertical, uint32_t c) override {
        if (s.empty()) return;
        br->SetColor(col(c));
        const double w = textWidth(s, size);
        const float sx = (w > 0 && widthPt > 0) ? float(std::max(0.4, std::min(2.5, widthPt / w))) : 1.f;
        D2D1_MATRIX_3X2_F old;
        rt->GetTransform(&old);
        D2D1_MATRIX_3X2_F m = D2D1::Matrix3x2F::Scale(sx, 1.f);
        if (vertical) m = m * D2D1::Matrix3x2F::Rotation(-90.f);
        rt->SetTransform(m * D2D1::Matrix3x2F::Translation(float(x), float(y)) * old);
        run(s, size, 0, 0);
        rt->SetTransform(old);
    }
    void pushClip(double x, double y, double w, double h) override {
        rt->PushAxisAlignedClip(D2D1::RectF(float(x), float(y), float(x + w), float(y + h)), D2D1_ANTIALIAS_MODE_ALIASED);
    }
    void popClip() override { rt->PopAxisAlignedClip(); }
};

} // namespace

bool writePNG(const Figure& fg, const Style& st, const std::string& path, std::string* err) {
    auto fail = [&](const char* what, HRESULT hr) {
        if (err) { char b[160]; std::snprintf(b, sizeof(b), "%s failed (0x%08lX)", what, (unsigned long)hr); *err = b; }
        return false;
    };
    const HRESULT co = CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED);
    const double dpi = std::max(72, std::min(st.dpi, 1200));
    const UINT W = UINT(std::lround(st.widthMM / 25.4 * dpi)), H = UINT(std::lround(st.heightMM / 25.4 * dpi));
    IWICImagingFactory* wic = nullptr;
    IWICBitmap* bmp = nullptr;
    D2DCanvas cv;
    bool ok = false;
    HRESULT hr = CoCreateInstance(CLSID_WICImagingFactory, nullptr, CLSCTX_INPROC_SERVER, IID_IWICImagingFactory, (void**)&wic);
    if (FAILED(hr)) { fail("WIC", hr); goto done; }
    hr = wic->CreateBitmap(W, H, GUID_WICPixelFormat32bppPBGRA, WICBitmapCacheOnLoad, &bmp);
    if (FAILED(hr)) { fail("CreateBitmap", hr); goto done; }
    hr = D2D1CreateFactory(D2D1_FACTORY_TYPE_SINGLE_THREADED, __uuidof(ID2D1Factory), nullptr, (void**)&cv.f);
    if (FAILED(hr)) { fail("Direct2D", hr); goto done; }
    hr = DWriteCreateFactory(DWRITE_FACTORY_TYPE_SHARED, __uuidof(IDWriteFactory), (IUnknown**)&cv.dw);
    if (FAILED(hr)) { fail("DirectWrite", hr); goto done; }
    {
        // render target in pixels (96 dpi DIPs); the canvas transform maps points -> pixels at the requested dpi
        const D2D1_RENDER_TARGET_PROPERTIES p = D2D1::RenderTargetProperties(
            D2D1_RENDER_TARGET_TYPE_SOFTWARE, D2D1::PixelFormat(DXGI_FORMAT_B8G8R8A8_UNORM, D2D1_ALPHA_MODE_PREMULTIPLIED), 96.f, 96.f);
        hr = cv.f->CreateWicBitmapRenderTarget(bmp, p, &cv.rt);
        if (FAILED(hr)) { fail("CreateWicBitmapRenderTarget", hr); goto done; }
        cv.ptScale = float(dpi / 72.0);
        render(fg, st, cv);
        hr = cv.rt->EndDraw();
        if (FAILED(hr)) { fail("Direct2D drawing", hr); goto done; }
    }
    {
        IWICStream* stream = nullptr;
        IWICBitmapEncoder* enc = nullptr;
        IWICBitmapFrameEncode* frame = nullptr;
        WICPixelFormatGUID fmt = GUID_WICPixelFormat32bppPBGRA;
        const std::wstring wpath = widen(path);
        hr = wic->CreateStream(&stream);
        if (SUCCEEDED(hr)) hr = stream->InitializeFromFilename(wpath.c_str(), GENERIC_WRITE);
        if (SUCCEEDED(hr)) hr = wic->CreateEncoder(GUID_ContainerFormatPng, nullptr, &enc);
        if (SUCCEEDED(hr)) hr = enc->Initialize(stream, WICBitmapEncoderNoCache);
        if (SUCCEEDED(hr)) hr = enc->CreateNewFrame(&frame, nullptr);
        if (SUCCEEDED(hr)) hr = frame->Initialize(nullptr);
        if (SUCCEEDED(hr)) hr = frame->SetSize(W, H);
        if (SUCCEEDED(hr)) hr = frame->SetResolution(dpi, dpi);
        if (SUCCEEDED(hr)) hr = frame->SetPixelFormat(&fmt);
        if (SUCCEEDED(hr)) hr = frame->WriteSource(bmp, nullptr);
        if (SUCCEEDED(hr)) hr = frame->Commit();
        if (SUCCEEDED(hr)) hr = enc->Commit();
        release(frame); release(enc); release(stream);
        if (FAILED(hr)) fail("PNG encoding", hr);
        else ok = true;
    }
done:
    release(cv.rt);
    release(cv.dw);
    release(cv.f);
    release(bmp);
    release(wic);
    if (SUCCEEDED(co)) CoUninitialize();
    return ok;
}

} // namespace fig
