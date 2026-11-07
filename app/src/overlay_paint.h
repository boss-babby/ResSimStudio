// C39: backend-neutral 2-D painter for the 3-D view overlays (title, legend, well / object labels, axis triad).
// The overlay layout is written once (Viewport3D::paintOverlay, GeoScene3D::paintOverlay) and drawn by
//   Windows        : Direct2D + DirectWrite into the resolved D3D11 texture (overlay_d2d.cpp)
//   Linux / macOS  : an ImDrawList rendered by the OpenGL 3 ImGui renderer into the resolved FBO (gl_view.cpp)
// so the overlay is part of the view texture and of the PNG snapshots on every platform.
// Coordinates are logical overlay units; the painter applies the snapshot scale.
#pragma once
#include <cstdio>
#include <string>

struct OvColor { float r = 1, g = 1, b = 1, a = 1; };
enum class OvFont { Small = 0, Label = 1, Title = 2 };

class OverlayPainter {
public:
    virtual ~OverlayPainter() = default;
    virtual void fillRoundRect(float x0, float y0, float x1, float y1, float radius, OvColor c) = 0;
    virtual void fillRect(float x0, float y0, float x1, float y1, OvColor c) = 0;
    virtual void strokeRect(float x0, float y0, float x1, float y1, OvColor c, float width) = 0;
    virtual void line(float x0, float y0, float x1, float y1, OvColor c, float width) = 0;
    virtual void fillCircle(float cx, float cy, float r, OvColor c) = 0;
    virtual float textWidth(const std::string& utf8, OvFont f) = 0;
    // text with its line box top-left at (x, y)
    virtual void text(const std::string& utf8, OvFont f, float x, float y, OvColor c) = 0;
};

// shared number format of the colour legends
inline std::string ovLegendNumber(float v) {
    char buf[32];
    const float av = v < 0 ? -v : v;
    if (av != 0 && (av < 0.01f || av >= 1e5f)) std::snprintf(buf, sizeof(buf), "%.2e", v);
    else if (av < 10) std::snprintf(buf, sizeof(buf), "%.3f", v);
    else std::snprintf(buf, sizeof(buf), "%.1f", v);
    return buf;
}

#if defined(_WIN32) && !defined(RESSIM_GL)
#include <d2d1_1.h>
#include <dwrite.h>
// Direct2D painter: draws into the device context's current target (BeginDraw / EndDraw are done by the caller)
class D2DPainter final : public OverlayPainter {
public:
    D2DPainter(ID2D1DeviceContext* d, IDWriteFactory* dw, IDWriteTextFormat* small, IDWriteTextFormat* label, IDWriteTextFormat* title);
    ~D2DPainter() override;
    void fillRoundRect(float x0, float y0, float x1, float y1, float radius, OvColor c) override;
    void fillRect(float x0, float y0, float x1, float y1, OvColor c) override;
    void strokeRect(float x0, float y0, float x1, float y1, OvColor c, float width) override;
    void line(float x0, float y0, float x1, float y1, OvColor c, float width) override;
    void fillCircle(float cx, float cy, float r, OvColor c) override;
    float textWidth(const std::string& utf8, OvFont f) override;
    void text(const std::string& utf8, OvFont f, float x, float y, OvColor c) override;
    bool ok() const { return br_ != nullptr; }

private:
    ID2D1SolidColorBrush* brush(OvColor c);
    ID2D1DeviceContext* d_;
    IDWriteFactory* dw_;
    IDWriteTextFormat* tf_[3];
    ID2D1SolidColorBrush* br_ = nullptr;
};
#endif
