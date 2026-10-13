// Publication-quality figure export: a small plotting model (panels of line / scatter series), a layout engine
// (nice ticks, scaled axes, legends, multi-panel grids) and vector writers (SVG, PDF 1.4 with base-14 fonts).
// PNG export renders the same layout through a raster Canvas (Direct2D + WIC on Windows, figure_png_win.cpp).
// The in-app plots record what they draw into fig::recorder(), so an export reproduces the visible plots.
// Portable C++17.
#pragma once
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace fig {

// colours are 0xAABBGGRR (same byte order as ImU32 / IM_COL32)
inline uint32_t rgba(int r, int g, int b, int a = 255) { return uint32_t(r) | (uint32_t(g) << 8) | (uint32_t(b) << 16) | (uint32_t(a) << 24); }

struct Series {
    std::string label;
    std::vector<double> x, y;
    uint32_t color = 0xFF000000;
    float width = 1.0f;        // relative to Style::lineWidthPt
    bool scatter = false;      // markers only
    bool markers = false;      // line + markers
    int dash = 0;              // 0 solid, 1 dashed, 2 dotted, 3 dash-dot
    int axis = 0;              // 0 left Y, 1 right Y
    bool legend = true;
};

struct Panel {
    std::string title, xlabel, ylabel, y2label;
    std::vector<Series> series;
    bool hasLimits = false, hasY2 = false;
    double x0 = 0, x1 = 1, y0 = 0, y1 = 1, y20 = 0, y21 = 1;
    std::vector<double> vlines;   // vertical marker lines (e.g. time cursor)
};

struct Figure {
    std::string name;
    std::vector<Panel> panels;
    int rows = 1, cols = 1;
    bool empty() const { return panels.empty(); }
};

enum class Legend : int { Auto = 0, InsideTopRight, InsideTopLeft, InsideBottomRight, InsideBottomLeft, OutsideRight, Below, None };

struct Style {
    double widthMM = 160, heightMM = 100;
    int font = 0;                   // 0 Helvetica / Arial, 1 Times / Times New Roman
    double fontPt = 9;
    double lineWidthPt = 1.2;
    bool grid = true, minorTicks = true, titles = true, box = true;
    Legend legend = Legend::Auto;
    int palette = -1;               // -1 keep plot colours, otherwise a palette index (palettes())
    int dpi = 600;                  // PNG
    bool transparent = false;       // PNG / SVG background
    bool useLimits = true;          // use the visible axis ranges of the app plots (else fit the data)
};

// Qualitative colour palettes for distinct curves. palettes()[i].colors; first entry is the default.
struct Palette { const char* name; std::vector<uint32_t> colors; };
const std::vector<Palette>& palettes();
// Distinct colour for curve i: palette colours first, then well separated golden-angle hues.
uint32_t distinctColor(int i, int palette = 0);

// Drawing target. Coordinates in points (1/72 inch), origin top-left, y down. Text may contain "^{...}"
// superscripts. halign: 0 left, 1 centre, 2 right; y is the text baseline; vertical text is rotated 90 deg CCW.
class Canvas {
public:
    virtual ~Canvas() = default;
    virtual void begin(double wPt, double hPt, bool transparent) = 0;
    virtual void end() = 0;
    virtual void polyline(const double* xy, int n, uint32_t color, double widthPt, int dash) = 0;
    virtual void fillRect(double x, double y, double w, double h, uint32_t color) = 0;
    virtual void strokeRect(double x, double y, double w, double h, uint32_t color, double widthPt) = 0;
    virtual void circle(double x, double y, double r, uint32_t fill, uint32_t stroke, double widthPt) = 0;
    virtual void text(double x, double y, const std::string& s, double sizePt, int halign, bool vertical, uint32_t color) = 0;
    virtual void pushClip(double x, double y, double w, double h) = 0;
    virtual void popClip() = 0;
    virtual double textWidth(const std::string& s, double sizePt) const;   // default: base-14 AFM metrics
    int font = 0;
};

// Font metrics (AFM widths of Helvetica / Times-Roman), handles "^{...}" superscripts.
double textWidthPt(const std::string& s, double sizePt, int font);

// Copy with the series colours mapped onto a palette (equal colours stay equal).
Figure recolored(const Figure& f, int palette);

// Lays out and draws the figure on a canvas of the style's size.
void render(const Figure& f, const Style& s, Canvas& c);

bool writeSVG(const Figure& f, const Style& s, const std::string& path, std::string* err = nullptr);
bool writePDF(const Figure& f, const Style& s, const std::string& path, std::string* err = nullptr);
std::string svgString(const Figure& f, const Style& s);
// CSV of all series (x, y columns) - the numbers behind a figure.
bool writeCSV(const Figure& f, const std::string& path, std::string* err = nullptr);

// Nice tick positions covering [a, b] (about 'target' ticks).
std::vector<double> niceTicks(double a, double b, int target, double* step = nullptr);
// Tick label formatting with a common power-of-ten factor for very large / small ranges.
std::string formatTick(double v, double step, int exp10);
int axisExponent(double a, double b);   // 0 or multiple of 3 when the axis needs a x10^n factor

// --------------------------------------------------------------------- recorder used by the app plots
class Recorder {
public:
    void beginFigure(const std::string& id, int rows = 1, int cols = 1);
    void beginPanel(const std::string& title, const std::string& xlabel, const std::string& ylabel, const std::string& y2label = std::string());
    void series(const std::string& label, const double* x, const double* y, int n, uint32_t color, float width, bool scatter = false,
                int axis = 0, bool markers = false, bool legend = true);
    Series* openSeries();       // the last added series (to append points), null if none
    void vline(double x);
    void limits(double x0, double x1, double y0, double y1, bool hasY2, double y20, double y21);
    void endPanel();
    void endFigure();
    const Figure* get(const std::string& id) const;
    bool recording() const { return active_; }
private:
    std::map<std::string, Figure> figs_;
    Figure cur_;
    Panel panel_;
    bool active_ = false, panelOpen_ = false;
};
Recorder& recorder();

} // namespace fig
