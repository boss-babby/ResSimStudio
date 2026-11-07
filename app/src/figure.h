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
    std::vector<double> ylo;   // C23: non-empty -> filled band between ylo and y (colour = fill, alpha kept)
    double barWidth = 0;       // C23: > 0 -> vertical bars of this width (data units) from 0 to y
    uint32_t barLine = 0;      //      bar outline colour (0 = none)
};

// C23 native publication layers (re-drawn layout): heatmaps, shapes in data coordinates, colour bar.
struct Heatmap {
    std::string label;
    int rows = 0, cols = 0;
    bool colMajor = false;
    std::vector<double> v;                  // rows x cols, row 0 at y1 (ImPlot convention)
    double x0 = 0, y0 = 0, x1 = 1, y1 = 1;  // bounds as passed to PlotHeatmap
    double vmin = 0, vmax = 1;
    std::vector<uint32_t> lut;              // colormap samples over [vmin, vmax]
    uint32_t colorAt(double val) const;
};
struct Shape {
    enum Kind : int { Polyline = 0, Polygon = 1, Marker = 2, Text = 3 };
    int kind = Polyline;
    std::vector<double> x, y;     // data coordinates (marker / text: first point)
    uint32_t stroke = 0xFF000000, fill = 0;
    double widthPx = 1;           // on-screen width (px) -> ~0.75 pt x Style::lineWidthPt / 1.2
    int dash = 0;
    int marker = 0;               // 0 circle, 1 square, 2 diamond, 3 triangle, 4 cross, 5 plus
    double sizePx = 4;            // marker radius on screen
    std::string text;
    double offX = 0, offY = 0;    // text offset (screen px)
    bool vertical = false;
    int halign = 1;               // text: 0 left, 1 centre, 2 right (vertically centred on the point)
    bool keepColor = false;       // false: light neutral strokes / text (made for the dark theme) print dark
    std::string label;            // legend entry (empty = none)
};
struct ColorBar {
    bool on = false;
    std::string label;
    double vmin = 0, vmax = 1;
    std::vector<uint32_t> lut;
};

struct Panel {
    std::string title, xlabel, ylabel, y2label;
    std::vector<Series> series;
    bool hasLimits = false, hasY2 = false;
    double x0 = 0, x1 = 1, y0 = 0, y1 = 1, y20 = 0, y21 = 1;
    std::vector<double> vlines;   // vertical marker lines (e.g. time cursor)
    bool logX = false, logY = false, logY2 = false;   // C17: logarithmic axes (decade ticks; non-positive values are skipped)
    // C23 native layers
    std::vector<Heatmap> heatmaps;    // drawn first (cell backgrounds)
    std::vector<Shape> shapes;        // polygons below the series, lines / markers / text above
    ColorBar cbar;                    // colour scale right of the frame
    bool equal = false;               // equal axis scaling (maps, plan views)
    bool invertX = false, invertY = false;
    bool cellEdges = false;           // thin outlines around heatmap cells
    bool noAxes = false;              // C23: no frame / ticks / axis labels (stereonets, schematic maps)
};

// C22: exact ("as on screen") snapshot of one or more app plots - the draw geometry of the plots exactly as the
// screen shows it (grid, heatmap / cell backgrounds, fills, custom markers, legends, labels), as filled triangles in
// screen pixels plus the text strings. Captured from the ImGui draw lists (plotkit.cpp), drawn by renderSnapshot().
struct SnapText {
    double x = 0, y = 0;          // baseline start (vertical: baseline x, start y - text runs upwards)
    double size = 13, width = 0;  // font size and the on-screen width of the string, px
    uint32_t color = 0xFF000000;
    bool vertical = false;
    int clip = -1;                // index into Snapshot::clips
    int order = 0;                // drawn just before runs[order] (keeps the on-screen stacking order)
    std::string s;
};
struct Snapshot {
    struct Run { uint32_t color; int clip; int tri0, ntri; };
    double x0 = 0, y0 = 0, x1 = 0, y1 = 0;   // bounds (screen px)
    uint32_t bg = 0xFF000000;                 // window background behind the plots
    std::vector<double> xy;                   // 6 per triangle (counter-clockwise in screen coordinates)
    std::vector<Run> runs;                    // consecutive triangles of one colour and clip rectangle
    std::vector<SnapText> texts;
    std::vector<double> clips;                // 4 per clip rectangle: x0, y0, x1, y1
    int plots = 0;
    std::vector<uint32_t> printFrom, printTo;   // theme colours -> light "print" colours (exact RGBA match)
    bool empty() const { return runs.empty() && texts.empty(); }
    double aspect() const { return x1 > x0 ? (y1 - y0) / (x1 - x0) : 0.6; }
};

struct Figure {
    std::string name;
    std::vector<Panel> panels;
    int rows = 1, cols = 1;
    Snapshot snap;                // C22: exact on-screen geometry (may be empty)
    bool empty() const { return panels.empty() && snap.empty(); }
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
    bool exact = false;             // C22: draw Figure::snap (exactly as on screen) instead of the re-drawn layout (C23: publication first)
    bool printColors = false;       // C22 exact mode: light background - neutral greys inverted, colours kept
    bool layers = true;             // C23 re-drawn layout: heatmaps, shapes, colour bars, bands, bars (natively, journal style)
    bool cellEdges = false;         // C23: outline heatmap cells (grid of the model cells)
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
    // C22 (exact snapshots): filled triangles (6 doubles each, counter-clockwise -> drawn as one non-zero-winding path),
    // and a plain string (no "^{}" parsing) stretched horizontally to the given width.
    virtual void fillTriangles(const double* xy, int nTri, uint32_t color);
    virtual void textFit(double x, double y, const std::string& s, double sizePt, double widthPt, bool vertical, uint32_t color);
    int font = 0;
};

// Font metrics (AFM widths of Helvetica / Times-Roman), handles "^{...}" superscripts.
double textWidthPt(const std::string& s, double sizePt, int font);

// Copy with the series colours mapped onto a palette (equal colours stay equal).
Figure recolored(const Figure& f, int palette);

// Lays out and draws the figure on a canvas of the style's size (Style::exact + a snapshot -> renderSnapshot).
void render(const Figure& f, const Style& s, Canvas& c);
// Draws an exact snapshot scaled to the style width (height = width x snapshot aspect, centred if the style differs).
void renderSnapshot(const Snapshot& sn, const Style& s, Canvas& c);
// exact mode "print colours": the theme colours of the snapshot (backgrounds, text, grid, borders) -> light ones,
// other neutral greys inverted, real colours (series, colour maps) kept
uint32_t printColor(uint32_t c, const Snapshot* sn = nullptr);

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
    Panel* openPanel();         // C23: panel being recorded (native layers), null if none
    Panel* lastPanel();         // C23: last finished panel of the current figure (colour scale drawn after EndPlot)
    void vline(double x);
    void limits(double x0, double x1, double y0, double y1, bool hasY2, double y20, double y21);
    void axisScale(bool logX, bool logY, bool logY2);   // C17
    void endPanel();
    void endFigure();
    const Figure* get(const std::string& id) const;
    bool recording() const { return active_; }
    const std::string& currentId() const { return cur_.name; }   // id of the figure being recorded
private:
    std::map<std::string, Figure> figs_;
    Figure cur_;
    Panel panel_;
    bool active_ = false, panelOpen_ = false;
};
Recorder& recorder();

} // namespace fig
