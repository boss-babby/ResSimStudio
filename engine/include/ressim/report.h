// C33 report generator: a run (or several) -> Word (.docx), PowerPoint (.pptx) or HTML report, written natively (no
// Office / Python needed). Charts are native, editable Office charts (DrawingML scatter charts with the data embedded
// as literals) in DOCX / PPTX and inline SVG in HTML; 2-D property maps are PNG images.
//
// The content is a simple block model (Doc) so that the app, the CLI, the Python API and the regression harness can
// build their own reports with the same writers.
#pragma once
#include "ressim/api.h"
#include "ressim/case.h"
#include "ressim/ecl.h"
#include <cstdint>
#include <string>
#include <vector>

namespace ressim {
namespace report {

struct Series {
    std::string name;
    std::vector<double> x, y;
    uint32_t rgb = 0;            // 0xRRGGBB, 0 = automatic palette colour
    bool dashed = false;
    bool y2 = false;             // plotted against the secondary axis
    bool markers = false;
};
struct Chart {
    std::string title, xLabel, yLabel, y2Label;
    std::vector<Series> series;
    bool logY = false;
};
struct Image {
    std::vector<uint8_t> png;    // encoded PNG
    int w = 0, h = 0;            // pixels
    std::string caption;
};
struct Table {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
    std::string caption;
    std::vector<int> highlightRows;      // drawn with a warning fill (e.g. failed regression cases)
};
struct Block {
    enum Kind { Heading, Para, Bullets, TableB, ChartB, ImageB, PageBreak } kind = Para;
    int level = 1;               // Heading level 1..3
    std::string text;
    std::vector<std::string> items;
    Table table;
    Chart chart;
    Image image;
};
struct Doc {
    std::string title, subtitle, author, date;
    std::vector<Block> blocks;
    // helpers
    void h(int level, const std::string& t) { Block b; b.kind = Block::Heading; b.level = level; b.text = t; blocks.push_back(b); }
    void p(const std::string& t) { Block b; b.kind = Block::Para; b.text = t; blocks.push_back(b); }
    void bullets(const std::vector<std::string>& v) { Block b; b.kind = Block::Bullets; b.items = v; blocks.push_back(b); }
    void table(const Table& t) { Block b; b.kind = Block::TableB; b.table = t; blocks.push_back(b); }
    void chart(const Chart& c) { Block b; b.kind = Block::ChartB; b.chart = c; blocks.push_back(b); }
    void image(const Image& i) { Block b; b.kind = Block::ImageB; b.image = i; blocks.push_back(b); }
    void pageBreak() { Block b; b.kind = Block::PageBreak; blocks.push_back(b); }
};

// ---- writers (format by extension: .docx .pptx .html/.htm .md); false + err on failure
bool writeDoc(const std::string& path, const Doc& d, std::string& err);
bool writeDocx(const std::string& path, const Doc& d, std::string& err);
bool writePptx(const std::string& path, const Doc& d, std::string& err);   // one slide per heading-1/2 section
bool writeHtml(const std::string& path, const Doc& d, std::string& err);   // self-contained (SVG charts, data-URI PNG)
bool writeMarkdown(const std::string& path, const Doc& d, std::string& err);   // tables + text (charts as CSV notes)

// ---- building blocks (exposed for tests and other generators)
std::string svgChart(const Chart& c, int w = 760, int h = 360);
std::vector<uint8_t> encodePng(const std::vector<uint8_t>& rgb, int w, int h);          // 8-bit RGB
uint32_t crc32(const uint8_t* p, size_t n, uint32_t crc = 0);
std::string deflateRaw(const uint8_t* data, size_t n);             // raw DEFLATE (LZ77 + fixed Huffman)
bool inflateRaw(const std::string& in, std::string& out);           // stored / fixed-Huffman blocks
struct ZipEntry { std::string name; std::string data; };
bool writeZip(const std::string& path, const std::vector<ZipEntry>& entries, std::string& err);   // deflated (PNG stored)
std::vector<ZipEntry> readZip(const std::string& path, std::string& err);                       // stored / fixed-Huffman entries (tests)
// map of one layer (k, 0-based) of a per-active-cell array -> PNG (turbo colour map, fixed range lo..hi; inactive grey)
Image layerMap(const std::vector<float>& activeValues, const std::vector<int>& globalOf, int nx, int ny, int k,
               double lo, double hi, int scale, const std::string& caption);
uint32_t paletteColor(int i);

// ---- simulation report
struct Options {
    std::string title;                  // default: case title
    std::string author;
    bool wellCharts = true;             // per-well rate / BHP charts (up to maxWells)
    int maxWells = 8;
    bool maps = true;                   // pressure / saturation maps of the first and last snapshot
    int mapLayer = -1;                  // -1: the layer with the largest initial oil pore volume
    bool inputs = true;                 // model description tables
    bool log = true;                    // engine messages (warnings, ACTIONX events)
    std::string notes;                  // free text paragraph after the summary
};
// reference: optional second summary (e.g. OPM Flow / ECLIPSE run or history) overlaid on the charts + deviation table
Doc simulationReport(const Case& c, const api::RunOutput& run, const ecl::Summary& sum, const Options& opt,
                     const ecl::Summary* reference = nullptr, const std::string& refLabel = "reference");
// options string "title=...;author=...;maps=0;wells=0;layer=3;maxwells=4;notes=..." -> Options
Options parseOptions(const std::string& s);

}  // namespace report
}  // namespace ressim
