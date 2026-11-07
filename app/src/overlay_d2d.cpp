// C39 Direct2D / DirectWrite overlay painter (Windows) - see overlay_paint.h
#include "gfx.h"
#include "overlay_paint.h"
#ifdef RESSIM_D3D11
#include <string>

static std::wstring ovWiden(const std::string& s) {
    if (s.empty()) return std::wstring();
    const int n = MultiByteToWideChar(CP_UTF8, 0, s.c_str(), int(s.size()), nullptr, 0);
    std::wstring w(size_t(n > 0 ? n : 0), L'\0');
    if (n > 0) MultiByteToWideChar(CP_UTF8, 0, s.c_str(), int(s.size()), &w[0], n);
    return w;
}

D2DPainter::D2DPainter(ID2D1DeviceContext* d, IDWriteFactory* dw, IDWriteTextFormat* small, IDWriteTextFormat* label, IDWriteTextFormat* title)
    : d_(d), dw_(dw), tf_{small, label, title} {
    if (d_) d_->CreateSolidColorBrush(D2D1::ColorF(1, 1, 1, 1), &br_);
}
D2DPainter::~D2DPainter() { SafeRelease(br_); }
ID2D1SolidColorBrush* D2DPainter::brush(OvColor c) { br_->SetColor(D2D1::ColorF(c.r, c.g, c.b, c.a)); return br_; }
void D2DPainter::fillRoundRect(float x0, float y0, float x1, float y1, float r, OvColor c) {
    d_->FillRoundedRectangle(D2D1::RoundedRect(D2D1::RectF(x0, y0, x1, y1), r, r), brush(c));
}
void D2DPainter::fillRect(float x0, float y0, float x1, float y1, OvColor c) { d_->FillRectangle(D2D1::RectF(x0, y0, x1, y1), brush(c)); }
void D2DPainter::strokeRect(float x0, float y0, float x1, float y1, OvColor c, float w) { d_->DrawRectangle(D2D1::RectF(x0, y0, x1, y1), brush(c), w); }
void D2DPainter::line(float x0, float y0, float x1, float y1, OvColor c, float w) { d_->DrawLine(D2D1::Point2F(x0, y0), D2D1::Point2F(x1, y1), brush(c), w); }
void D2DPainter::fillCircle(float cx, float cy, float r, OvColor c) { d_->FillEllipse(D2D1::Ellipse(D2D1::Point2F(cx, cy), r, r), brush(c)); }
float D2DPainter::textWidth(const std::string& s, OvFont f) {
    IDWriteTextFormat* tf = tf_[int(f)];
    if (!tf || !dw_) return 0;
    const std::wstring t = ovWiden(s);
    IDWriteTextLayout* lay = nullptr; float wd = 0;
    if (SUCCEEDED(dw_->CreateTextLayout(t.c_str(), UINT32(t.size()), tf, 4000, 200, &lay))) {
        DWRITE_TEXT_METRICS m; lay->GetMetrics(&m); wd = m.widthIncludingTrailingWhitespace; lay->Release();
    }
    return wd;
}
void D2DPainter::text(const std::string& s, OvFont f, float x, float y, OvColor c) {
    IDWriteTextFormat* tf = tf_[int(f)];
    if (!tf) return;
    const std::wstring t = ovWiden(s);
    d_->DrawText(t.c_str(), UINT32(t.size()), tf, D2D1::RectF(x, y, x + 4000, y + 400), brush(c));
}
#endif
