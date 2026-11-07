// Deck Editor workspace ("writing mode"): a full-window ECLIPSE / OPM Flow deck editor with syntax colouring,
// keyword / value / well / group autocomplete, hover documentation with the active parameter, live diagnostics
// (squiggles + line markers), outline, keyword browser, documentation panel, schema-driven keyword forms,
// record templates, column alignment, INCLUDE navigation, and Apply-to-model / Run.
// Text editing is provided by ImGuiColorTextEdit (goossens, MIT) in third_party/textedit.
#include "app.h"
#include "platform.h"
#include "widgets.h"
#include "theme.h"
#include "gfx.h"
#include "deck_intel.h"
#include "TextEditor.h"
#include "ressim/deck_parser.h"
#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <set>
#include <sstream>
#include "resources.h"
#include "icons.h"
#include <unordered_set>

using namespace ressim;
namespace di = ressim::deckintel;
namespace fs = std::filesystem;

namespace {

std::string exeDir() { return plat::exeDir(); }
void ensureDocs() {
    static bool tried = false;
    if (tried) return;
    tried = true;
    // the manual is embedded in the executable; an external keyword_docs.txt (next to the exe) overrides it
    if (di::Docs::get().load(exeDir() + "keyword_docs.txt")) return;
    di::Docs::get().loadText(appResource(RES_KEYWORD_DOCS));
}
std::string fileName(const std::string& p) { const size_t s = p.find_last_of("\\/"); return s == std::string::npos ? p : p.substr(s + 1); }
std::string dirOf(const std::string& p) { const size_t s = p.find_last_of("\\/"); return s == std::string::npos ? std::string() : p.substr(0, s + 1); }
bool readFile(const std::string& p, std::string& out) {
    std::ifstream f(fs::u8path(p), std::ios::binary);
    if (!f) return false;
    std::stringstream ss;
    ss << f.rdbuf();
    out = ss.str();
    // CRLF -> LF (the editor works with LF; files are written back with LF which every simulator accepts)
    out.erase(std::remove(out.begin(), out.end(), '\r'), out.end());
    return true;
}
bool writeFile(const std::string& p, const std::string& s) {
    std::error_code ec;
    fs::create_directories(fs::u8path(p).parent_path(), ec);
    std::ofstream f(fs::u8path(p), std::ios::binary);
    if (!f) return false;
    f.write(s.data(), std::streamsize(s.size()));
    return bool(f);
}
const std::unordered_set<std::string>& knownKeywords() {
    static std::unordered_set<std::string> s;
    if (s.empty()) {
        for (const auto& kv : di::Docs::get().all()) s.insert(kv.first);
        for (const auto& k : deck::keywordTable()) s.insert(k.name);
    }
    return s;
}

// ECLIPSE language for the colorizer: keywords only in column 1, sections, numbers, n* repeats, '/' and
// everything after it on the line as comment.
const TextEditor::Language* eclipseLanguage() {
    static TextEditor::Language L;
    static bool init = false;
    if (init) return &L;
    init = true;
    L.name = "ECLIPSE";
    L.caseSensitive = true;
    L.singleLineComment = "--";
    L.hasSingleQuotedStrings = true;
    L.hasDoubleQuotedStrings = false;
    L.customTokenizer = [](TextEditor::Iterator start, TextEditor::Iterator end, TextEditor::Color& color) -> TextEditor::Iterator {
        using C = TextEditor::Color;
        auto& info = textEditorTokenizerInfo;
        if (info.lineFlags == 1) { color = C::comment; return end; }
        const ImWchar c = *start;
        auto isWord = [](ImWchar ch) { return (ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z') || (ch >= '0' && ch <= '9') || ch == '_' || ch == '+' || ch == '-'; };
        if (c == '/') { info.lineFlags = 1; color = C::punctuation; auto i = start; ++i; return i; }
        if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) {
            std::string w;
            auto i = start;
            while (i < end && isWord(*i)) {
                if (*i == '-') { auto j = i; ++j; if (j < end && *j == '-') break; }
                w += char(*i);
                ++i;
            }
            if (info.column == 0 && w.size() <= 8) {
                if (di::isSectionName(w)) color = C::declaration;
                else if (knownKeywords().count(w)) color = C::keyword;
                else color = C::knownIdentifier;
                // text after a keyword is ignored by ECLIPSE (generator stamps in old decks): show it as a comment
                // unless the line looks like a data record (contains '/')
                if (color != C::knownIdentifier && i < end) {
                    bool slash = false;
                    for (auto j = i; j < end; ++j) if (*j == '/') { slash = true; break; }
                    if (!slash) info.lineFlags = 1;
                }
            } else color = C::identifier;
            return i;
        }
        const bool num = (c >= '0' && c <= '9') || c == '.' || ((c == '-' || c == '+') && [&]() { auto j = start; ++j; return j < end && ((*j >= '0' && *j <= '9') || *j == '.'); }());
        if (num || c == '*') {
            auto i = start;
            bool star = false;
            while (i < end) {
                const ImWchar ch = *i;
                if ((ch >= '0' && ch <= '9') || ch == '.' || ch == 'e' || ch == 'E' || ch == 'd' || ch == 'D' || ch == '+' || ch == '-') { ++i; continue; }
                if (ch == '*') { star = true; ++i; continue; }
                break;
            }
            color = star ? C::preprocessor : C::number;
            return i;
        }
        return start;
    };
    return &L;
}

TextEditor::Palette eclipsePalette() {
    TextEditor::Palette p = TextEditor::GetDarkPalette();
    auto set = [&](TextEditor::Color c, ImU32 v) { p[size_t(c)] = v; };
    set(TextEditor::Color::keyword, IM_COL32(86, 170, 255, 255));
    set(TextEditor::Color::declaration, IM_COL32(255, 184, 72, 255));
    set(TextEditor::Color::knownIdentifier, IM_COL32(220, 130, 255, 255));
    set(TextEditor::Color::identifier, IM_COL32(214, 218, 226, 255));
    set(TextEditor::Color::number, IM_COL32(152, 220, 140, 255));
    set(TextEditor::Color::string, IM_COL32(236, 170, 120, 255));
    set(TextEditor::Color::punctuation, IM_COL32(255, 110, 120, 255));
    set(TextEditor::Color::preprocessor, IM_COL32(190, 150, 255, 255));
    set(TextEditor::Color::comment, IM_COL32(112, 128, 120, 255));
    set(TextEditor::Color::background, IM_COL32(20, 22, 27, 255));
    set(TextEditor::Color::lineNumber, IM_COL32(90, 96, 110, 255));
    set(TextEditor::Color::currentLineNumber, IM_COL32(200, 205, 215, 255));
    set(TextEditor::Color::currentLineHighlight, IM_COL32(40, 44, 54, 120));
    set(TextEditor::Color::cursor, IM_COL32(205, 225, 255, 255));   // bright caret
    return p;
}

const char* kSkeleton =
    "-- New ECLIPSE / OPM Flow deck. Type a keyword in column 1 (Ctrl+Space for suggestions);\n"
    "-- accepting a suggestion inserts a record template with column headers.\n"
    "RUNSPEC\n\nTITLE\n   New model\n\nDIMENS\n-- NX NY NZ\n   10 10 3 /\n\nOIL\nWATER\nGAS\nDISGAS\n\nFIELD\n\nSTART\n   1 'JAN' 2025 /\n\n"
    "GRID\n\nPROPS\n\nSOLUTION\n\nSUMMARY\n\nSCHEDULE\n\nEND\n";

} // namespace

// ============================================================================ editor state
struct EdDoc {
    std::string path, title;
    TextEditor ed;
    size_t savedUndo = 0;
    size_t analyzedUndo = size_t(-1), seenUndo = size_t(-1);
    double changedAt = 0;
    std::vector<std::string> lines;
    di::Outline outline;
    std::vector<di::Diag> diags, parserDiags;
    bool big = false;
    std::string prevLineText;
    int prevLine = -1;
    // snippet session of a freshly inserted record template: Tab / Shift+Tab walk its values
    struct Snippet { bool active = false; int first = 0, last = 0; } snip;
    bool modified() const { return ed.GetUndoIndex() != savedUndo; }
};

struct App::DeckEditorState {
    std::vector<std::unique_ptr<EdDoc>> docs;
    int active = -1, wantActive = -1;
    TextEditor::AutoCompleteConfig ac;
    bool showLeft = true, showRight = true, showProblems = true, autoTemplate = true, minimap = false, whitespace = false;
    int leftTab = 0, rightTab = 0;
    int wantRightTab = -1;                 // 0 documentation, 1 form (programmatic switch, e.g. from the context menu)
    char outlineFilter[64] = {0}, kwFilter[64] = {0};
    float fontScale = 1.0f;
    // keyword form
    di::Form form;
    int formLine = -1;
    size_t formUndo = size_t(-1);
    bool formDirty = false;
    int closeAsk = -1;
    std::string lastApplied;
    unsigned appliedVersion = 0;
    EdDoc* doc() { return active >= 0 && active < int(docs.size()) ? docs[size_t(active)].get() : nullptr; }
};

namespace {

void analyze(EdDoc& d, bool force) {
    const size_t u = d.ed.GetUndoIndex();
    if (u == d.analyzedUndo && !force) return;
    const double now = ImGui::GetTime();
    if (!force && d.analyzedUndo != size_t(-1)) {
        // debounce: analyse when typing pauses
        if (u != d.seenUndo) { d.seenUndo = u; d.changedAt = now; return; }
        if (now - d.changedAt < (d.big ? 1.2 : 0.35)) return;
    }
    d.analyzedUndo = u;
    d.changedAt = 0;
    d.lines = di::splitLines(d.ed.GetText());
    d.outline = di::buildOutline(d.lines);
    d.diags = di::lint(d.lines, d.outline);
    d.ed.ClearSquiggles(1);
    d.ed.ClearMarkers();
    std::set<int> marked;
    auto add = [&](const di::Diag& g) {
        const ImU32 c = g.severity >= 2 ? IM_COL32(255, 90, 90, 255) : g.severity == 1 ? IM_COL32(255, 186, 70, 255) : IM_COL32(120, 170, 255, 255);
        const int len = g.line < int(d.lines.size()) ? int(d.lines[size_t(g.line)].size()) : 0;
        const int a = std::min(g.col0, len), b = std::max(a + 1, std::min(std::max(g.col1, g.col0 + 1), std::max(len, 1)));
        d.ed.AddSquiggle(TextEditor::DocPos(size_t(g.line), size_t(a)), TextEditor::DocPos(size_t(g.line), size_t(b)), 1, c, g.text);
        if (g.severity >= 1 && marked.insert(g.line).second) d.ed.AddMarker(size_t(g.line), c, 0, g.text, std::string());
    };
    for (const auto& g : d.diags) add(g);
    for (const auto& g : d.parserDiags) add(g);
}

std::string unitSystem(const EdDoc& d) {
    for (const auto& it : d.outline.items) {
        const std::string u = di::upper(it.name);
        if (u == "METRIC") return "METRIC";
        if (u == "FIELD") return "FIELD";
        if (u == "LAB") return "LAB";
        if (di::upper(it.name) == "GRID") break;
    }
    return "METRIC";
}

void paramTooltip(const di::Param& p, const std::string& units) {
    ImGui::TextColored(col4(Pal::Accent), "%d  %s", p.index, p.name.c_str());
    if (!p.manual.empty() && p.manual != p.name) { ImGui::SameLine(); ImGui::TextDisabled("(%s)", p.manual.c_str()); }
    const std::string u = units == "FIELD" ? p.unitsField : p.unitsMetric;
    if (!u.empty()) ImGui::Text("Units (%s): %s", units.c_str(), u.c_str());
    if (!p.def.empty()) ImGui::Text("Default: %s", p.def.c_str());
    if (!p.valueType.empty()) ImGui::TextDisabled("Type: %s%s%s", p.valueType.c_str(), p.dimension.empty() ? "" : "   dimension: ", p.dimension.c_str());
    if (!p.options.empty()) {
        std::string o;
        for (size_t i = 0; i < p.options.size(); ++i) o += (i ? ", " : "") + p.options[i];
        ImGui::TextWrapped("Options: %s", o.c_str());
    }
    if (!p.desc.empty()) ImGui::TextWrapped("%s", p.desc.c_str());
}

} // namespace

// ============================================================================ documents
namespace { void deckFolds(size_t n, const std::function<std::string(size_t)>& text, std::vector<std::pair<size_t, size_t>>& folds);
            bool snippetTab(EdDoc& d, bool back); }
static EdDoc* newDoc(App::DeckEditorState& S, const std::string& path, const std::string& text, const TextEditor::AutoCompleteConfig* ac) {
    auto d = std::make_unique<EdDoc>();
    d->path = path;
    d->title = path.empty() ? "untitled.DATA" : fileName(path);
    d->big = text.size() > 4u * 1024 * 1024;
    d->ed.SetLanguage(eclipseLanguage());
    d->ed.SetPalette(eclipsePalette());
    d->ed.SetTabSize(4);
    d->ed.SetShowWhitespacesEnabled(S.whitespace);
    d->ed.SetInsertSpacesOnTabs(true);
    d->ed.SetShowMatchingBrackets(!d->big);
    if (!d->big) { d->ed.SetLineFoldingEnabled(true); d->ed.SetFoldProvider(deckFolds); }
    d->ed.SetCompletePairedGlyphs(false);
    d->ed.SetAutoIndentEnabled(false);
    d->ed.SetText(text);
    d->ed.SetAutoCompleteConfig(ac);
    d->savedUndo = d->ed.GetUndoIndex();
    EdDoc* dp = d.get();
    d->ed.SetTabHandler([dp](bool shift) { return snippetTab(*dp, shift); });
    S.docs.push_back(std::move(d));
    S.wantActive = int(S.docs.size()) - 1;
    S.active = S.wantActive;
    return S.docs.back().get();
}

void App::edEnsure() {
    if (ed_) return;
    ensureDocs();
    ed_ = std::make_shared<DeckEditorState>();
    DeckEditorState* S = ed_.get();
    S->ac.triggerOnTyping = true;
    S->ac.triggerInComments = false;
    S->ac.triggerInStrings = true;
    S->ac.triggerDelay = std::chrono::milliseconds(120);
    S->ac.suggestionWidth = 26;
    S->ac.callback = [S](TextEditor::AutoCompleteState& st) {
        EdDoc* d = S->doc();
        st.suggestions.clear();
        if (!d) return;
        const int line = int(st.searchTermStart.line);
        if (line >= int(d->lines.size())) d->lines.resize(size_t(line) + 1);
        d->lines[size_t(line)] = d->ed.GetLineText(size_t(line));
        st.suggestions = di::suggest(d->lines, d->outline, line, int(st.searchTermStart.index), st.searchTerm);
    };
}

void App::edOpenFile(const std::string& path) {
    edEnsure();
    for (size_t i = 0; i < ed_->docs.size(); ++i)
        if (ed_->docs[i]->path == path) { ed_->wantActive = int(i); workspace_ = Workspace::DeckEditor; return; }
    std::string text;
    std::error_code ec;
    const auto sz = fs::file_size(fs::u8path(path), ec);
    if (!ec && sz > 64ull * 1024 * 1024) { notify("File too large for the editor (" + std::to_string(sz >> 20) + " MB)"); return; }
    if (!readFile(path, text)) { notify("Cannot open " + path); return; }
    newDoc(*ed_, path, text, &ed_->ac);
    workspace_ = Workspace::DeckEditor;
}

void App::edOpenModelDeck() {
    edEnsure();
    std::string text;
    try { text = eclipseDeckToString(case_); }
    catch (const std::exception& e) { notify(std::string("Cannot write the model deck: ") + e.what()); return; }
    // large arrays go to include files next to the deck (the editor stays fast and the deck stays readable)
    const std::string dir = runDirFor(path_, case_.title) + "\\deck\\";
    const std::string base = caseBaseName(case_);
    std::vector<std::pair<std::string, std::string>> files;
    const std::string main = di::externalizeLargeKeywords(text, "INCLUDE", base, 256 * 1024, files);
    bool ok = writeFile(dir + base + ".DATA", main + "\n");
    for (const auto& f : files) ok = ok && writeFile(dir + f.first, f.second);
    if (!ok) { notify("Cannot write to " + dir); return; }
    // replace an already open copy
    for (size_t i = 0; i < ed_->docs.size(); ++i) {
        if (ed_->docs[i]->path == dir + base + ".DATA") {
            EdDoc& d = *ed_->docs[i];
            d.ed.SetText(main + "\n");
            d.savedUndo = d.ed.GetUndoIndex();
            d.analyzedUndo = size_t(-1);
            ed_->wantActive = int(i);
            workspace_ = Workspace::DeckEditor;
            notify("Model deck regenerated");
            return;
        }
    }
    edOpenFile(dir + base + ".DATA");
    notify("Model deck written to " + dir + (files.empty() ? std::string() : " (" + std::to_string(files.size()) + " include files)"));
}

bool App::edSave(int idx, bool saveAs) {
    if (!ed_ || idx < 0 || idx >= int(ed_->docs.size())) return false;
    EdDoc& d = *ed_->docs[size_t(idx)];
    std::string p = d.path;
    if (p.empty() || saveAs) {
        if (p.empty()) p = d.title;
        if (!fileDialog(true, "ECLIPSE deck (*.DATA)\0*.DATA\0Include file (*.INC)\0*.INC;*.inc;*.GRDECL\0All files (*.*)\0*.*\0", "DATA", p)) return false;
    }
    if (!writeFile(p, d.ed.GetText())) { notify("Cannot write " + p); return false; }
    d.path = p;
    d.title = fileName(p);
    d.savedUndo = d.ed.GetUndoIndex();
    notify("Saved " + d.title);
    return true;
}

// ============================================================================ editing commands
namespace {

struct BlockRange { int kwLine = -1, first = -1, last = -1; std::string keyword; };

// value tokens of lines [first, last]: whitespace separated, quoted strings are one token, comments and the
// record terminator "/" are skipped
struct Tok { int line, c0, c1; };
std::vector<Tok> snippetTokens(const TextEditor& ed, int first, int last) {
    std::vector<Tok> t;
    for (int l = first; l <= last && l < int(ed.GetLineCount()); ++l) {
        const std::string s = ed.GetLineText(size_t(l));
        size_t i = 0;
        while (i < s.size()) {
            if (s[i] == ' ' || s[i] == '\t' || s[i] == '\r') { ++i; continue; }
            if (s.compare(i, 2, "--") == 0) break;
            size_t j = i;
            if (s[i] == '\'') { j = s.find('\'', i + 1); j = j == std::string::npos ? s.size() : j + 1; }
            else while (j < s.size() && s[j] != ' ' && s[j] != '\t' && s[j] != '\r' && s[j] != '\'' && s.compare(j, 2, "--") != 0) ++j;
            if (j == i) j = i + 1;
            const std::string tok = s.substr(i, j - i);
            if (tok != "/") t.push_back({l, int(i), int(j)});
            i = j;
        }
    }
    return t;
}
bool snippetTab(EdDoc& d, bool back) {
    if (!d.snip.active) return false;
    const auto cur = d.ed.GetMainCursorPosition();
    const int cl = int(cur.line), cc = int(cur.index);
    if (cl < d.snip.first || cl > d.snip.last) { d.snip.active = false; return false; }
    const auto T = snippetTokens(d.ed, d.snip.first, d.snip.last);
    if (T.empty()) { d.snip.active = false; return false; }
    int at = -1;
    for (int k = 0; k < int(T.size()); ++k) if (T[size_t(k)].line == cl && T[size_t(k)].c0 <= cc && cc <= T[size_t(k)].c1) { at = k; break; }
    int next = -1;
    if (back) {
        if (at > 0) next = at - 1;
        else if (at < 0) for (int k = int(T.size()) - 1; k >= 0; --k) if (T[size_t(k)].line < cl || (T[size_t(k)].line == cl && T[size_t(k)].c1 < cc)) { next = k; break; }
        if (next < 0) next = 0;
    } else {
        if (at >= 0) next = at + 1;
        else for (int k = 0; k < int(T.size()); ++k) if (T[size_t(k)].line > cl || (T[size_t(k)].line == cl && T[size_t(k)].c0 > cc)) { next = k; break; }
        if (next < 0 || next >= int(T.size())) {
            // past the last value: finish the session at the end of the block
            d.snip.active = false;
            d.ed.SetCursor(TextEditor::DocPos(size_t(d.snip.last), d.ed.GetLineText(size_t(d.snip.last)).size()));
            return true;
        }
    }
    const Tok& k = T[size_t(next)];
    d.ed.SelectRegion(TextEditor::DocPos(size_t(k.line), size_t(k.c0)), TextEditor::DocPos(size_t(k.line), size_t(k.c1)));
    return true;
}

// folds: deck sections (header to the next section) and keyword blocks (keyword to its last data line)
void deckFolds(size_t n, const std::function<std::string(size_t)>& text, std::vector<std::pair<size_t, size_t>>& folds) {
    static const std::set<std::string> sections = {"RUNSPEC", "GRID", "EDIT", "PROPS", "REGIONS", "SOLUTION", "SUMMARY", "SCHEDULE", "OPTIMIZE"};
    int sec = -1, kw = -1, lastData = -1;
    auto closeKw = [&]() { if (kw >= 0 && lastData > kw) folds.push_back({size_t(kw), size_t(lastData)}); kw = -1; };
    auto closeSec = [&](int endLine) {
        int e = endLine;
        while (e > sec) { const std::string t = text(size_t(e)); if (t.find_first_not_of(" \t\r") != std::string::npos) break; --e; }
        if (sec >= 0 && e > sec) folds.push_back({size_t(sec), size_t(e)});
        sec = -1;
    };
    for (size_t l = 0; l < n; ++l) {
        const std::string t = text(l);
        const size_t a = t.find_first_not_of(" \t\r");
        if (a == std::string::npos || t.compare(a, 2, "--") == 0) continue;
        std::string k;
        if (a == 0 && di::keywordLine(t, &k) && knownKeywords().count(di::upper(k))) {
            closeKw();
            if (sections.count(di::upper(k))) { closeSec(int(l) - 1); sec = int(l); }
            else kw = int(l);
            lastData = int(l);
            continue;
        }
        lastData = int(l);
    }
    closeKw();
    closeSec(int(n) - 1);
}
BlockRange blockAt(const EdDoc& d, int line) {
    BlockRange b;
    const int i = d.outline.itemAt(line);
    if (i < 0) return b;
    b.kwLine = d.outline.items[size_t(i)].line;
    b.keyword = di::upper(d.outline.items[size_t(i)].name);
    b.first = b.kwLine + 1;
    b.last = d.outline.bodyEnd(i, int(d.lines.size()));
    // trailing blank lines stay outside the block
    while (b.last >= b.first && d.lines[size_t(b.last)].find_first_not_of(" \t\r") == std::string::npos) --b.last;
    return b;
}
void replaceLines(EdDoc& d, int first, int last, const std::vector<std::string>& lines) {
    const std::string text = di::joinLines(lines);
    if (last < first) {
        // insert after line first-1
        const size_t at = size_t(first - 1);
        const size_t len = d.ed.GetLineText(at).size();
        d.ed.ReplaceSectionText(TextEditor::DocPos(at, len), TextEditor::DocPos(at, len), "\n" + text);
        return;
    }
    const size_t len = d.ed.GetLineText(size_t(last)).size();
    d.ed.ReplaceSectionText(TextEditor::DocPos(size_t(first), 0), TextEditor::DocPos(size_t(last), len), text);
}
void insertTemplate(EdDoc& d, int kwLine, const std::string& kw) {
    int nx = 0, ny = 0, nz = 0;
    const int cells = di::dimensCells(d.lines, d.outline, &nx, &ny, &nz);
    std::string t = di::recordTemplate(kw, cells);
    if (t.empty()) return;
    while (!t.empty() && t.back() == '\n') t.pop_back();
    const size_t len = d.ed.GetLineText(size_t(kwLine)).size();
    d.ed.ReplaceSectionText(TextEditor::DocPos(size_t(kwLine), len), TextEditor::DocPos(size_t(kwLine), len), "\n" + t);
    // snippet session: the first value is selected, Tab / Shift+Tab move between the values, Esc ends
    const auto tl = di::splitLines(t);
    d.snip.active = true;
    d.snip.first = kwLine + 1;
    d.snip.last = kwLine + int(tl.size());
    d.ed.SetCursor(TextEditor::DocPos(size_t(kwLine) + 1, 0));
    const auto T = snippetTokens(d.ed, d.snip.first, d.snip.last);
    if (T.empty()) d.snip.active = false;
    else d.ed.SelectRegion(TextEditor::DocPos(size_t(T[0].line), size_t(T[0].c0)), TextEditor::DocPos(size_t(T[0].line), size_t(T[0].c1)));
}

} // namespace

// ============================================================================ workspace
// Schema-driven record editor of a keyword (records as tables with parameter names, units, defaults and option
// lists). Shared by the Deck Editor's Form tab and the Advanced keywords page. Returns true when F changed.
bool deckFormRecords(di::Form& F, const std::string& units) {
    const float fs = ImGui::GetFontSize();
    bool changed = false;
    if (ImGui::Button("+ Record")) {
        std::vector<std::string> r;
        if (F.table) r.assign(size_t(std::max(1, F.columns)), "0");
        else if (F.info) for (const auto* p : F.info->recordParams(int(F.records.size()))) r.push_back(p->required() ? (p->isString() ? "''" : "0") : "1*");
        F.records.push_back(r);
        changed = true;
    }
    ImGui::BeginChild("##form", ImVec2(0, 0), 0, ImGuiWindowFlags_HorizontalScrollbar);
    int del = -1;
    for (size_t r = 0; r < F.records.size(); ++r) {
        auto& rec = F.records[r];
        ImGui::PushID(int(r));
        const auto ps = F.info ? F.info->recordParams(int(r)) : std::vector<const di::Param*>();
        const int cols = F.table ? std::max(1, F.columns) : std::max<int>(int(ps.size()), int(rec.size()));
        if (!F.table) while (int(rec.size()) < cols) rec.push_back("1*");
        ImGui::TextDisabled(F.table ? "Table %zu" : "Record %zu", r + 1);
        ImGui::SameLine();
        if (ImGui::SmallButton("x")) del = int(r);
        const int rows = F.table ? (int(rec.size()) + cols - 1) / cols : 1;
        const float tableH = std::min(fs * 22.0f, (rows + 1) * ImGui::GetFrameHeightWithSpacing() + ImGui::GetStyle().ScrollbarSize + 4.0f);
        if (ui::beginTable("##rec", cols + (F.table ? 1 : 0), ImGuiTableFlags_Borders | ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_ScrollX | ImGuiTableFlags_ScrollY, ImVec2(0, tableH))) {
            ImGui::TableSetupScrollFreeze(0, 1);
            for (int c = 0; c < cols; ++c) {
                const di::Param* p = F.table ? (F.info && !F.info->params.empty() ? &F.info->params[size_t(c) % F.info->params.size()] : nullptr)
                                             : (c < int(ps.size()) ? ps[size_t(c)] : nullptr);
                ui::tableColumn(p ? p->name.c_str() : "?", ImGuiTableColumnFlags_WidthFixed, fs * 6.5f);
            }
            if (F.table) ui::tableColumn("", ImGuiTableColumnFlags_WidthFixed, fs * 1.5f);
            ImGui::TableNextRow(ImGuiTableRowFlags_Headers);
            for (int c = 0; c < cols; ++c) {
                ImGui::TableSetColumnIndex(c);
                const di::Param* p = F.table ? (F.info && !F.info->params.empty() ? &F.info->params[size_t(c) % F.info->params.size()] : nullptr)
                                             : (c < int(ps.size()) ? ps[size_t(c)] : nullptr);
                ImGui::TableHeader(p ? p->name.c_str() : "?");
                if (p && ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) {
                    ImGui::BeginTooltip();
                    ImGui::PushTextWrapPos(fs * 30);
                    paramTooltip(*p, units);
                    ImGui::PopTextWrapPos();
                    ImGui::EndTooltip();
                }
            }
            int delRow = -1;
            for (int row = 0; row < rows; ++row) {
                ImGui::TableNextRow();
                for (int c = 0; c < cols; ++c) {
                    const size_t k = size_t(row * (F.table ? cols : 0) + c);
                    if (k >= rec.size()) break;
                    ImGui::TableSetColumnIndex(c);
                    ImGui::PushID(int(k));
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    const di::Param* p = F.table ? nullptr : (c < int(ps.size()) ? ps[size_t(c)] : nullptr);
                    if (p && !p->options.empty()) {
                        if (ImGui::BeginCombo("##v", rec[k].c_str())) {
                            if (ImGui::Selectable("1*  (default)", rec[k] == "1*")) { rec[k] = "1*"; changed = true; }
                            for (const auto& o : p->options) {
                                const std::string q = "'" + o + "'";
                                if (ImGui::Selectable(o.c_str(), di::unquote(rec[k]) == o)) { rec[k] = q; changed = true; }
                            }
                            ImGui::EndCombo();
                        }
                    } else if (ui::inputString("##v", &rec[k])) changed = true;
                    ImGui::PopID();
                }
                if (F.table) {
                    ImGui::TableSetColumnIndex(cols);
                    ImGui::PushID(row + 100000);
                    if (ImGui::SmallButton("-")) delRow = row;
                    ImGui::PopID();
                }
            }
            ui::endTable();
            if (delRow >= 0) {
                const size_t a = size_t(delRow * cols), e = std::min(rec.size(), a + size_t(cols));
                rec.erase(rec.begin() + long(a), rec.begin() + long(e));
                changed = true;
            }
        }
        if (F.table && ImGui::SmallButton("+ row")) { for (int c = 0; c < cols; ++c) rec.push_back("0"); changed = true; }
        ImGui::PopID();
    }
    if (del >= 0) { F.records.erase(F.records.begin() + del); changed = true; }
    ImGui::EndChild();
    return changed;
}

void App::deckEditorWorkspace() {
    edEnsure();
    DeckEditorState& S = *ed_;
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(vp->WorkPos);
    ImGui::SetNextWindowSize(vp->WorkSize);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0);
    ImGui::Begin("##deckeditor", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoSavedSettings |
                                              ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoMove);
    ImGui::PopStyleVar(2);
    const float fs = ImGui::GetFontSize();
    EdDoc* d = S.doc();
    auto tip = [](const char* t) { if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("%s", TR(t)); };

    // ---------------- command bar
    if (ImGui::Button(TRL(ICON_FILE_PLUS " New"))) { newDoc(S, std::string(), kSkeleton, &S.ac); d = S.doc(); }
    tip("New deck from a skeleton (Ctrl+N)");
    ImGui::SameLine();
    if (ImGui::Button(TRL(ICON_FOLDER_OPEN " Open..."))) {
        std::string p;
        if (fileDialog(false, "ECLIPSE deck / include\0*.DATA;*.data;*.INC;*.inc;*.GRDECL;*.grdecl;*.SCH;*.sch\0All files (*.*)\0*.*\0", "DATA", p)) edOpenFile(p);
    }
    tip("Open a deck or include file (Ctrl+O)");
    ImGui::SameLine();
    if (ImGui::Button(TRL(ICON_BOX " Model deck"))) edOpenModelDeck();
    tip("Write the current model as an ECLIPSE deck (large arrays as INCLUDE files) and open it");
    ImGui::SameLine();
    ImGui::BeginDisabled(!d);
    if (ImGui::Button(TRL(ICON_SAVE " Save"))) edSave(S.active, false);
    tip("Save (Ctrl+S)");
    ImGui::SameLine();
    if (ImGui::Button(TRL(ICON_SAVE " Save as..."))) edSave(S.active, true);
    ImGui::SameLine();
    ImGui::TextDisabled("|");
    ImGui::SameLine();
    bool doCheck = ImGui::Button(TRL(ICON_LIST_CHECKS " Check"));
    tip("Full check with the deck parser (F7): record structure, INCLUDE files, keyword support");
    ImGui::SameLine();
    bool doAlign = ImGui::Button(TRL(ICON_TABLE " Align"));
    tip("Align the columns of the keyword at the cursor and add column headers (Ctrl+Shift+A)");
    ImGui::SameLine();
    bool doTemplate = ImGui::Button(TRL(ICON_WAND_SPARKLES " Template"));
    tip("Insert a record template for the keyword on the cursor line (Ctrl+Shift+T)");
    ImGui::SameLine();
    ImGui::TextDisabled("|");
    ImGui::SameLine();
    bool doApply = ui::accentButton(ICON_UPLOAD " Apply to model");
    tip("Save and import this deck into the model (Ctrl+Enter). The Model and Results workspaces use it.");
    ImGui::SameLine();
    bool doRun = ImGui::Button(TRL(ICON_PLAY " Run"));
    tip("Apply to model and run it with the selected engine (F5)");
    ImGui::EndDisabled();
    ImGui::SameLine();
    if (ImGui::Button(TRL(ICON_EYE " View"))) ImGui::OpenPopup("edview");
    if (ImGui::BeginPopup("edview")) {
        ImGui::Checkbox(TRL("Outline / keywords panel"), &S.showLeft);
        ImGui::Checkbox(TRL("Documentation / form panel"), &S.showRight);
        ImGui::Checkbox(TRL("Problems panel"), &S.showProblems);
        ImGui::Separator();
        ImGui::Checkbox(TRL("Insert record template after keyword completion"), &S.autoTemplate);
        if (ImGui::Checkbox(TRL("Mini map"), &S.minimap)) for (auto& x : S.docs) x->ed.SetShowMiniMapEnabled(S.minimap);
        if (ImGui::Checkbox(TRL("Show whitespace"), &S.whitespace)) for (auto& x : S.docs) x->ed.SetShowWhitespacesEnabled(S.whitespace);
        ImGui::SetNextItemWidth(fs * 8);
        ImGui::SliderFloat("Text size", &S.fontScale, 0.7f, 2.0f, "%.2fx");
        ImGui::EndPopup();
    }
    if (d) {
        const auto pos = d->ed.GetMainCursorPosition();
        char b[64];
        std::snprintf(b, sizeof(b), "Ln %zu, Col %zu", pos.line + 1, pos.index + 1);
        ImGui::SameLine(ImGui::GetWindowWidth() - ImGui::CalcTextSize(b).x - fs * 1.2f);
        ImGui::TextDisabled("%s", b);
    }

    // shortcuts (editor workspace only)
    const ImGuiInputFlags G = ImGuiInputFlags_RouteGlobal;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_S, G) && d) edSave(S.active, false);
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiMod_Shift | ImGuiKey_S, G) && d) edSave(S.active, true);
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_O, G)) {
        std::string p;
        if (fileDialog(false, "ECLIPSE deck / include\0*.DATA;*.data;*.INC;*.inc;*.GRDECL;*.grdecl\0All files (*.*)\0*.*\0", "DATA", p)) edOpenFile(p);
    }
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_N, G)) { newDoc(S, std::string(), kSkeleton, &S.ac); d = S.doc(); }
    if (ImGui::Shortcut(ImGuiKey_F7, G)) doCheck = true;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiMod_Shift | ImGuiKey_A, G)) doAlign = true;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiMod_Shift | ImGuiKey_T, G)) doTemplate = true;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_Enter, G)) doApply = true;
    if (ImGui::Shortcut(ImGuiKey_F5, G)) doRun = true;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_W, G) && d) S.closeAsk = S.active;

    // ---------------- document tabs
    if (ImGui::BeginTabBar("##eddocs", ImGuiTabBarFlags_Reorderable | ImGuiTabBarFlags_AutoSelectNewTabs | ImGuiTabBarFlags_FittingPolicyScroll)) {
        for (size_t i = 0; i < S.docs.size(); ++i) {
            EdDoc& x = *S.docs[i];
            bool open = true;
            ImGuiTabItemFlags fl = (x.modified() ? ImGuiTabItemFlags_UnsavedDocument : 0) | (S.wantActive == int(i) ? ImGuiTabItemFlags_SetSelected : 0);
            const std::string label = x.title + "###doc" + std::to_string(reinterpret_cast<uintptr_t>(&x));
            if (ui::tabItem(label.c_str(), &open, fl)) { S.active = int(i); ImGui::EndTabItem(); }
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("%s", x.path.empty() ? "not saved" : x.path.c_str());
            if (!open) S.closeAsk = int(i);
        }
        S.wantActive = -1;
        ImGui::EndTabBar();
    }
    if (S.closeAsk >= 0 && S.closeAsk < int(S.docs.size())) {
        EdDoc& x = *S.docs[size_t(S.closeAsk)];
        int r = IDNO;
        if (x.modified()) r = plat::msgBox(("Save changes to " + x.title + "?").c_str(), "Deck Editor", MB_YESNOCANCEL | MB_ICONQUESTION);
        if (r == IDYES && !edSave(S.closeAsk, false)) r = IDCANCEL;
        if (r != IDCANCEL) {
            S.docs.erase(S.docs.begin() + S.closeAsk);
            S.active = std::min(S.active, int(S.docs.size()) - 1);
            S.formLine = -1;
        }
        S.closeAsk = -1;
        d = S.doc();
    }
    d = S.doc();
    if (!d) {
        ImGui::Dummy(ImVec2(0, fs * 3));
        ImGui::Indent(fs * 2);
        if (g_fonts.title) ImGui::PushFont(g_fonts.title, 0.0f);
        ImGui::TextUnformatted(TR("Deck Editor"));
        if (g_fonts.title) ImGui::PopFont();
        ImGui::TextDisabled("%s", TR("Intelligent ECLIPSE / OPM Flow deck editing: keyword and value completion, record templates, hover"));
        ImGui::TextDisabled("%s", TR("documentation, live checks, keyword forms, INCLUDE navigation, and Apply / Run."));
        ImGui::Dummy(ImVec2(0, fs));
        if (ui::accentButton("Open the current model as a deck")) edOpenModelDeck();
        ImGui::SameLine();
        if (ImGui::Button(TRL("Open a deck file..."))) {
            std::string p;
            if (fileDialog(false, "ECLIPSE deck / include\0*.DATA;*.data;*.INC;*.inc;*.GRDECL;*.grdecl\0All files (*.*)\0*.*\0", "DATA", p)) edOpenFile(p);
        }
        ImGui::SameLine();
        if (ImGui::Button(TRL("New deck"))) newDoc(S, std::string(), kSkeleton, &S.ac);
        ImGui::Unindent(fs * 2);
        ImGui::End();
        return;
    }
    analyze(*d, false);

    // ---------------- context at the cursor
    const auto cur = d->ed.GetMainCursorPosition();
    const int cline = int(cur.line);
    if (cline < int(d->lines.size())) d->lines[size_t(cline)] = d->ed.GetLineText(cur.line);
    const di::Context ctx = di::contextAt(d->lines, d->outline, cline, int(cur.index));
    const std::string units = unitSystem(*d);

    // ---------------- commands needing the analysis
    if (doTemplate && cline < int(d->lines.size())) {
        std::string kw;
        if (di::keywordLine(d->lines[size_t(cline)], &kw)) insertTemplate(*d, cline, di::upper(kw));
        else notify("Put the cursor on a keyword line (keyword in column 1)");
    }
    if (doAlign) {
        analyze(*d, true);
        const BlockRange b = blockAt(*d, cline);
        if (b.kwLine < 0 || b.last < b.first) notify("No keyword records at the cursor");
        else {
            std::vector<std::string> body(d->lines.begin() + b.first, d->lines.begin() + b.last + 1);
            const auto al = di::alignBody(b.keyword, body, true);
            if (al != body) replaceLines(*d, b.first, b.last, al);
        }
    }
    if (doCheck) {
        analyze(*d, true);
        d->parserDiags.clear();
        try {
            const deck::RawDeck raw = deck::parseDeckText(d->ed.GetText(), dirOf(d->path), d->title);
            int n = 0;
            for (const auto& m : raw.messages) {
                if (m.level == deck::Level::Info) continue;
                if (m.file != 0) { ++n; continue; }
                d->parserDiags.push_back({std::max(0, m.line - 1), 0, 200, m.level == deck::Level::Error ? 2 : 1, (m.keyword.empty() ? "" : m.keyword + ": ") + m.text});
            }
            notify("Check: " + std::to_string(raw.keywords.size()) + " keywords, " + std::to_string(d->parserDiags.size() + d->diags.size()) +
                   " problem(s) in this file" + (n ? ", " + std::to_string(n) + " in include files" : std::string()));
        } catch (const std::exception& e) { notify(std::string("Check failed: ") + e.what()); }
        d->analyzedUndo = size_t(-1);
        analyze(*d, true);
    }
    if ((doApply || doRun) && d) {
        if (d->path.empty() || d->modified()) { if (!edSave(S.active, false)) doApply = doRun = false; }
        // the model is replaced by the deck (as File > Open of a deck): ask about unsaved model edits, unless the
        // model is still exactly what the last Apply produced
        if ((doApply || doRun) && !(S.appliedVersion == caseVersion_ && S.appliedVersion != 0) && !maybeSaveChanges()) doApply = doRun = false;
        if (doApply || doRun) {
            dirty_ = false;
            importDeck(d->path);
            S.lastApplied = d->path;
            S.appliedVersion = caseVersion_;
            if (doRun) runSimulation();
        }
    }

    // ---------------- layout: left | editor | right, problems below
    const float probH = S.showProblems ? fs * 7.5f : 0.0f;
    const float leftW = S.showLeft ? fs * 15 : 0, rightW = S.showRight ? fs * 24 : 0;
    const ImVec2 avail = ImGui::GetContentRegionAvail();
    const float mainH = avail.y - probH - (S.showProblems ? ImGui::GetStyle().ItemSpacing.y : 0);
    int jumpLine = -1;

    if (S.showLeft) {
        ImGui::BeginChild("##edleft", ImVec2(leftW, mainH), ImGuiChildFlags_Borders | ImGuiChildFlags_ResizeX);
        if (ImGui::BeginTabBar("##edlefttabs")) {
            if (ui::tabItem("Outline")) {
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::InputTextWithHint("##of", "filter keywords", S.outlineFilter, sizeof(S.outlineFilter));
                ImGui::BeginChild("##outline");
                const std::string f = di::upper(S.outlineFilter);
                std::vector<int> diagLines;
                for (const auto& g : d->diags) if (g.severity >= 1) diagLines.push_back(g.line);
                const int curItem = d->outline.itemAt(cline);
                bool sectionOpen = true;
                for (size_t i = 0; i < d->outline.items.size(); ++i) {
                    const auto& it = d->outline.items[i];
                    if (it.section) {
                        if (!f.empty()) { sectionOpen = true; continue; }
                        ImGui::SetNextItemOpen(true, ImGuiCond_Once);
                        sectionOpen = ImGui::TreeNodeEx((it.name + "##s" + std::to_string(i)).c_str(), ImGuiTreeNodeFlags_SpanAvailWidth | ImGuiTreeNodeFlags_NoTreePushOnOpen);
                        if (ImGui::IsItemClicked()) jumpLine = it.line;
                        continue;
                    }
                    if (!sectionOpen) continue;
                    if (!f.empty() && di::upper(it.name).find(f) == std::string::npos) continue;
                    const int end = d->outline.bodyEnd(int(i), int(d->lines.size()));
                    const bool bad = std::any_of(diagLines.begin(), diagLines.end(), [&](int l) { return l >= it.line && l <= end; });
                    ImGui::Indent(fs * 0.8f);
                    if (bad) ImGui::PushStyleColor(ImGuiCol_Text, col4(Pal::Warn));
                    if (ImGui::Selectable((it.name + "##k" + std::to_string(i)).c_str(), int(i) == curItem)) jumpLine = it.line;
                    if (bad) ImGui::PopStyleColor();
                    if (int(i) == curItem && ImGui::IsWindowAppearing()) ImGui::SetScrollHereY();
                    ImGui::SameLine(ImGui::GetContentRegionMax().x - fs * 3);
                    ImGui::TextDisabled("%d", it.line + 1);
                    ImGui::Unindent(fs * 0.8f);
                }
                if (d->outline.items.empty()) ImGui::TextDisabled("No keywords yet");
                ImGui::EndChild();
                ImGui::EndTabItem();
            }
            if (ui::tabItem("Keywords")) {
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::InputTextWithHint("##kf", "search name or summary", S.kwFilter, sizeof(S.kwFilter));
                const std::string section = ctx.section;
                ImGui::TextDisabled("%s keywords - double-click inserts", section.empty() ? "All" : section.c_str());
                ImGui::BeginChild("##kwlist");
                const std::string f = di::upper(S.kwFilter);
                for (const auto& kv : di::Docs::get().all()) {
                    const di::KwInfo& k = kv.second;
                    if (!k.alias.empty()) continue;
                    if (!section.empty() && !k.validIn(section) && f.empty()) continue;
                    if (!f.empty() && kv.first.find(f) == std::string::npos && di::upper(k.summary).find(f) == std::string::npos) continue;
                    const deck::Support sup = deck::ressimSupport(kv.first);
                    const bool supported = sup == deck::Support::Full || sup == deck::Support::Partial;
                    ImGui::PushStyleColor(ImGuiCol_Text, supported ? ImGui::GetStyleColorVec4(ImGuiCol_Text) : col4(Pal::Muted));
                    ImGui::Selectable(kv.first.c_str(), false, ImGuiSelectableFlags_AllowDoubleClick);
                    ImGui::PopStyleColor();
                    if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) {
                        ImGui::BeginTooltip();
                        ImGui::PushTextWrapPos(fs * 30);
                        ImGui::TextUnformatted(k.summary.c_str());
                        ImGui::TextDisabled("Sections: %s   ResSim: %s", k.sections.c_str(), deck::toString(sup));
                        ImGui::PopTextWrapPos();
                        ImGui::EndTooltip();
                    }
                    if (ImGui::IsItemHovered() && ImGui::IsMouseDoubleClicked(0)) {
                        // insert after the keyword block at the cursor (or on the empty cursor line)
                        int at = cline;
                        const bool empty = cline < int(d->lines.size()) && d->lines[size_t(cline)].find_first_not_of(" \t") == std::string::npos;
                        if (!empty) { const BlockRange b = blockAt(*d, cline); at = b.kwLine < 0 ? cline : std::max(b.last, b.kwLine); }
                        if (empty) {
                            const size_t len = d->ed.GetLineText(size_t(at)).size();
                            d->ed.ReplaceSectionText(TextEditor::DocPos(size_t(at), 0), TextEditor::DocPos(size_t(at), len), kv.first);
                        } else {
                            const size_t len = d->ed.GetLineText(size_t(at)).size();
                            d->ed.ReplaceSectionText(TextEditor::DocPos(size_t(at), len), TextEditor::DocPos(size_t(at), len), "\n\n" + kv.first);
                            at += 2;
                        }
                        analyze(*d, true);
                        insertTemplate(*d, at, kv.first);
                        d->ed.SetFocus();
                    }
                }
                ImGui::EndChild();
                ImGui::EndTabItem();
            }
            ImGui::EndTabBar();
        }
        ImGui::EndChild();
        ImGui::SameLine();
    }

    // ---------------- editor
    {
        const float w = ImGui::GetContentRegionAvail().x - (S.showRight ? rightW + ImGui::GetStyle().ItemSpacing.x : 0);
        if (g_fonts.mono) ImGui::PushFont(g_fonts.mono, g_fonts.mono->LegacySize * S.fontScale);
        EdDoc* dp = d;
        d->ed.SetTextHoverCallback([this, dp, units](TextEditor::PopupData& pd) {
            // close when the mouse left the word (and is not over the popup)
            const ImVec2 mp = ImGui::GetMousePos();
            if (!ImGui::IsWindowHovered()) {
                const auto at = dp->ed.FindWordStart(dp->ed.GetDocPosAtMousePos(mp), true);
                if (at.line != pd.pos.line || at.index != pd.pos.index) { ImGui::CloseCurrentPopup(); return; }
            }
            const int line = int(pd.pos.line);
            if (line >= int(dp->lines.size())) { ImGui::CloseCurrentPopup(); return; }
            const di::Context hc = di::contextAt(dp->lines, dp->outline, line, int(pd.pos.index));
            ImGui::PushTextWrapPos(ImGui::GetFontSize() * 34);
            bool shown = false;
            for (const auto* list : {&dp->diags, &dp->parserDiags})
                for (const auto& g : *list)
                    if (g.line == line && int(pd.pos.index) >= g.col0 && int(pd.pos.index) <= std::max(g.col1, g.col0 + 1)) {
                        ImGui::TextColored(col4(g.severity >= 2 ? Pal::Err : g.severity == 1 ? Pal::Warn : Pal::Accent), "%s", g.text.c_str());
                        shown = true;
                    }
            if (hc.onKeywordLine && hc.info) {
                if (shown) ImGui::Separator();
                if (g_fonts.bold) ImGui::PushFont(g_fonts.bold, 0.0f);
                ImGui::TextUnformatted(hc.keyword.c_str());
                if (g_fonts.bold) ImGui::PopFont();
                ImGui::SameLine();
                ImGui::TextDisabled("%s  |  ResSim: %s", hc.info->sections.c_str(), deck::toString(deck::ressimSupport(hc.keyword)));
                ImGui::TextUnformatted(hc.info->summary.c_str());
                if (!hc.info->params.empty()) {
                    ImGui::Separator();
                    int shownP = 0;
                    for (const auto& p : hc.info->params) {
                        if (++shownP > 14) { ImGui::TextDisabled("... %zu items (Documentation panel)", hc.info->params.size()); break; }
                        const std::string u = units == "FIELD" ? p.unitsField : p.unitsMetric;
                        ImGui::Text("%s%2d  %-18s %s", p.record > 0 ? ("R" + std::to_string(p.record) + " ").c_str() : "", p.index, p.name.c_str(), u.c_str());
                    }
                }
                shown = true;
            } else if (hc.param && !hc.inComment) {
                if (shown) ImGui::Separator();
                ImGui::TextDisabled("%s  record %d  item %d", hc.keyword.c_str(), hc.record + 1, hc.itemIndex + 1);
                paramTooltip(*hc.param, units);
                shown = true;
            } else if (hc.arrayKeyword && hc.arrayIndex >= 0 && !hc.inComment) {
                int nx = 0, ny = 0, nz = 0;
                di::dimensCells(dp->lines, dp->outline, &nx, &ny, &nz);
                if (shown) ImGui::Separator();
                ImGui::TextDisabled("%s  value %d", hc.keyword.c_str(), hc.arrayIndex + 1);
                if (nx > 0 && ny > 0) {
                    const int q = hc.arrayIndex;
                    ImGui::Text("cell (I, J, K) = (%d, %d, %d)", q % nx + 1, (q / nx) % ny + 1, q / (nx * ny) + 1);
                }
                if (hc.info) ImGui::TextUnformatted(hc.info->summary.c_str());
                shown = true;
            }
            ImGui::PopTextWrapPos();
            if (!shown) ImGui::CloseCurrentPopup();
        });
        // right-click menu on the text
        d->ed.SetTextContextMenuCallback([this, d, &S](TextEditor::PopupData& pd) {
            const int line = int(pd.pos.line);
            analyze(*d, false);
            const BlockRange b = line < int(d->lines.size()) ? blockAt(*d, line) : BlockRange();
            if (b.kwLine >= 0) {
                ImGui::TextDisabled("%s", b.keyword.c_str());
                const di::KwInfo* kd = di::Docs::get().find(b.keyword);
                if (kd && !kd->summary.empty()) { ImGui::SameLine(); ImGui::TextDisabled("- %.60s", kd->summary.c_str()); }
                if (ImGui::MenuItem(ICON_BOOK_OPEN "  Keyword documentation")) { d->ed.SetCursor(TextEditor::DocPos(size_t(line), 0)); S.showRight = true; S.wantRightTab = 0; }
                if (ImGui::MenuItem(ICON_TABLE "  Edit in form")) { d->ed.SetCursor(TextEditor::DocPos(size_t(b.kwLine), 0)); S.showRight = true; S.wantRightTab = 1; }
                if (ImGui::MenuItem(ICON_WAND_SPARKLES "  Insert record template", nullptr, false, !di::recordTemplate(b.keyword, 1).empty()))
                    insertTemplate(*d, b.kwLine, b.keyword);
                if (ImGui::MenuItem(ICON_TEXT_SEARCH "  Align columns", nullptr, false, b.last >= b.first)) {
                    std::vector<std::string> body(d->lines.begin() + b.first, d->lines.begin() + b.last + 1);
                    const auto al = di::alignBody(b.keyword, body, true);
                    if (al != body) replaceLines(*d, b.first, b.last, al);
                }
                if (ImGui::MenuItem(ICON_LIST "  Select keyword block")) d->ed.SelectLines(size_t(b.kwLine), size_t(std::max(b.kwLine, b.last)));
                ImGui::Separator();
            }
            if (ImGui::MenuItem(ICON_SCISSORS "  Cut", "Ctrl+X", false, d->ed.AnyCursorHasSelection())) d->ed.Cut();
            if (ImGui::MenuItem(ICON_COPY "  Copy", "Ctrl+C", false, d->ed.AnyCursorHasSelection())) d->ed.Copy();
            if (ImGui::MenuItem(ICON_CLIPBOARD_PASTE "  Paste", "Ctrl+V")) d->ed.Paste();
            if (ImGui::MenuItem(ICON_CODE "  Toggle comment", "Ctrl+/")) d->ed.ToggleComments();
            ImGui::Separator();
            if (d->ed.IsLineFoldingEnabled()) {
                if (ImGui::MenuItem(ICON_FOLD_VERTICAL "  Fold all")) d->ed.FoldAll();
                if (ImGui::MenuItem(ICON_UNFOLD_VERTICAL "  Unfold all")) d->ed.UnfoldAll();
            }
            if (ImGui::MenuItem(ICON_SEARCH "  Find / replace", "Ctrl+F")) d->ed.OpenFindReplaceWindow();
        });
        d->ed.Render("##text", ImVec2(w, mainH), ImGuiChildFlags_Borders);
        if (d->snip.active) {
            const char* h = "Template: Tab / Shift+Tab = next / previous value, Esc = done";
            const ImVec2 r1 = ImGui::GetItemRectMax(), ts = ImGui::CalcTextSize(h);
            const ImVec2 q(r1.x - ts.x - 24, r1.y - ts.y - 14);
            ImGui::GetWindowDrawList()->AddRectFilled(ImVec2(q.x - 8, q.y - 4), ImVec2(q.x + ts.x + 8, q.y + ts.y + 4), IM_COL32(40, 70, 130, 230), 4.0f);
            ImGui::GetWindowDrawList()->AddText(q, IM_COL32(235, 240, 255, 255), h);
        }
        // snippet session ends with Esc or when the cursor leaves the inserted records
        if (d->snip.active) {
            const int l = int(d->ed.GetMainCursorPosition().line);
            if (ImGui::IsKeyPressed(ImGuiKey_Escape, false) || l < d->snip.first || l > d->snip.last) d->snip.active = false;
        }
        // Ctrl+click on an INCLUDE opens the file
        if (ImGui::IsItemHovered() && ImGui::GetIO().KeyCtrl && ImGui::IsMouseClicked(0)) {
            const auto mp = d->ed.GetDocPosAtMousePos(ImGui::GetMousePos());
            std::string inc = di::includePath(d->lines, d->outline, int(mp.line));
            if (!inc.empty()) {
                std::replace(inc.begin(), inc.end(), '/', '\\');
                const std::string full = fs::u8path(inc).is_absolute() ? inc : dirOf(d->path) + inc;
                rvDeferred_ = [this, full] { edOpenFile(full); };
            }
        }
        if (g_fonts.mono) ImGui::PopFont();
        // auto template: a keyword was just completed (several characters inserted at once) on an empty line
        const auto pos = d->ed.GetMainCursorPosition();
        const std::string lt = d->ed.GetLineText(pos.line);
        if (S.autoTemplate && int(pos.line) == d->prevLine && lt != d->prevLineText && lt.size() >= d->prevLineText.size() + 2 && pos.index == lt.size()) {
            std::string kw;
            const std::string prevU = di::upper(d->prevLineText);
            if (di::keywordLine(lt, &kw) && kw == di::upper(kw) && di::Docs::get().find(kw) && di::upper(lt).compare(0, prevU.size(), prevU) == 0) {
                const bool nextFree = pos.line + 1 >= d->ed.GetLineCount() || [&]() {
                    const std::string nl = d->ed.GetLineText(pos.line + 1);
                    std::string k2;
                    return nl.find_first_not_of(" \t") == std::string::npos || di::keywordLine(nl, &k2);
                }();
                if (nextFree) { analyze(*d, true); insertTemplate(*d, int(pos.line), kw); }
            }
        }
        const auto pos2 = d->ed.GetMainCursorPosition();
        d->prevLine = int(pos2.line);
        d->prevLineText = d->ed.GetLineText(pos2.line);
    }

    // ---------------- right: documentation / form
    if (S.showRight) {
        ImGui::SameLine();
        ImGui::BeginChild("##edright", ImVec2(0, mainH), ImGuiChildFlags_Borders);
        if (ImGui::BeginTabBar("##edrighttabs")) {
            const int wantTab = S.wantRightTab; S.wantRightTab = -1;
            if (ui::tabItem("Documentation", nullptr, wantTab == 0 ? ImGuiTabItemFlags_SetSelected : 0)) {
                ImGui::BeginChild("##docs");
                if (!ctx.info) {
                    if (ctx.keyword.empty()) ui::hint("Move the cursor into a keyword to see its documentation.");
                    else ImGui::TextColored(col4(Pal::Warn), "%s: no documentation (unknown keyword)", ctx.keyword.c_str());
                } else {
                    const di::KwInfo& k = *ctx.info;
                    if (g_fonts.title) ImGui::PushFont(g_fonts.title, 0.0f);
                    ImGui::TextUnformatted(k.name.c_str());
                    if (g_fonts.title) ImGui::PopFont();
                    const deck::Support sup = deck::ressimSupport(k.name);
                    ImGui::TextDisabled("Sections: %s", k.sections.c_str());
                    ImGui::TextColored(col4(sup == deck::Support::Full ? Pal::Ok : sup == deck::Support::Partial ? Pal::Warn : Pal::Muted), "ResSim: %s", deck::toString(sup));
                    if (sup != deck::Support::Full && sup != deck::Support::Partial) {
                        ImGui::SameLine();
                        ImGui::TextDisabled("(passed through to OPM Flow)");
                    }
                    ImGui::Separator();
                    ImGui::TextWrapped("%s", k.summary.c_str());
                    if (!ctx.onKeywordLine && ctx.itemIndex >= 0 && !ctx.arrayKeyword)
                        ImGui::TextColored(col4(Pal::Accent), "Cursor: record %d, item %d%s%s", ctx.record + 1, ctx.itemIndex + 1, ctx.param ? " - " : "",
                                           ctx.param ? ctx.param->name.c_str() : "");
                    if (!k.params.empty() && ui::beginTable("##params", 4, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerH | ImGuiTableFlags_SizingFixedFit)) {
                        ui::tableColumn("#");
                        ui::tableColumn("Item");
                        ui::tableColumn("Units");
                        ui::tableColumn("Default", ImGuiTableColumnFlags_WidthStretch);
                        ImGui::TableHeadersRow();
                        int lastRec = -1;
                        for (const auto& p : k.params) {
                            if (p.record != lastRec && k.maxRecord() > 0) {
                                lastRec = p.record;
                                ImGui::TableNextRow();
                                ImGui::TableSetColumnIndex(1);
                                ImGui::TextDisabled("Record %d", p.record);
                            }
                            ImGui::TableNextRow();
                            const bool act = ctx.param == &p;
                            if (act) ImGui::TableSetBgColor(ImGuiTableBgTarget_RowBg0, IM_COL32(64, 132, 246, 90));
                            ImGui::TableSetColumnIndex(0);
                            ImGui::Text("%d", p.index);
                            ImGui::TableSetColumnIndex(1);
                            ImGui::TextUnformatted(p.name.c_str());
                            if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) {
                                ImGui::BeginTooltip();
                                ImGui::PushTextWrapPos(fs * 32);
                                paramTooltip(p, units);
                                ImGui::PopTextWrapPos();
                                ImGui::EndTooltip();
                            }
                            if (act && ImGui::IsWindowAppearing()) ImGui::SetScrollHereY();
                            ImGui::TableSetColumnIndex(2);
                            ImGui::TextUnformatted((units == "FIELD" ? p.unitsField : p.unitsMetric).c_str());
                            ImGui::TableSetColumnIndex(3);
                            ImGui::TextUnformatted(p.def.c_str());
                        }
                        ui::endTable();
                    }
                    if (ctx.param) {
                        ImGui::Separator();
                        paramTooltip(*ctx.param, units);
                    }
                    if (!k.example.empty() && ImGui::CollapsingHeader("Example")) {
                        if (g_fonts.mono) ImGui::PushFont(g_fonts.mono, 0.0f);
                        ImGui::TextUnformatted(k.example.c_str());
                        if (g_fonts.mono) ImGui::PopFont();
                    }
                }
                ImGui::EndChild();
                ImGui::EndTabItem();
            }
            if (ui::tabItem("Form", nullptr, wantTab == 1 ? ImGuiTabItemFlags_SetSelected : 0)) {
                // schema-driven form of the keyword at the cursor
                const BlockRange b = blockAt(*d, cline);
                if (b.kwLine < 0) ui::hint("Move the cursor into a keyword to edit its records in a form.");
                else {
                    if (b.kwLine != S.formLine || (!S.formDirty && S.formUndo != d->ed.GetUndoIndex())) {
                        std::vector<std::string> body;
                        if (b.last >= b.first) body.assign(d->lines.begin() + b.first, d->lines.begin() + b.last + 1);
                        S.form = di::parseForm(b.keyword, body);
                        S.formLine = b.kwLine;
                        S.formUndo = d->ed.GetUndoIndex();
                        S.formDirty = false;
                    }
                    di::Form& F = S.form;
                    ImGui::Text("%s", F.keyword.c_str());
                    ImGui::SameLine();
                    ImGui::TextDisabled("line %d", b.kwLine + 1);
                    if (!F.table && !F.editable()) {
                        ui::hint("This keyword has no record schema for a form (arrays, tables of tables, code blocks): edit it as text.");
                    } else {
                        ImGui::BeginDisabled(!S.formDirty);
                        if (ui::accentButton("Apply to text")) {
                            const auto lines = di::writeForm(F);
                            replaceLines(*d, b.first, b.last, lines);
                            S.formDirty = false;
                            S.formUndo = size_t(-1);
                        }
                        ImGui::EndDisabled();
                        ImGui::SameLine();
                        ImGui::SameLine();
                        if (S.formDirty && ImGui::Button("Revert")) S.formLine = -1;
                        ImGui::SameLine();
                        if (deckFormRecords(F, units)) S.formDirty = true;
                    }
                }
                ImGui::EndTabItem();
            }
            ImGui::EndTabBar();
        }
        ImGui::EndChild();
    }

    // ---------------- problems
    if (S.showProblems) {
        ImGui::BeginChild("##edprob", ImVec2(0, 0), ImGuiChildFlags_Borders);
        const size_t nd = d->diags.size() + d->parserDiags.size();
        int ne = 0, nw = 0;
        for (const auto* l : {&d->diags, &d->parserDiags}) for (const auto& g : *l) { if (g.severity >= 2) ++ne; else if (g.severity == 1) ++nw; }
        ImGui::TextColored(col4(ne ? Pal::Err : nw ? Pal::Warn : Pal::Ok), "%s", nd == 0 ? "\xE2\x97\x8F No problems" : "\xE2\x97\x8F Problems");
        ImGui::SameLine();
        ImGui::TextDisabled("%d error(s), %d warning(s), %zu hint(s)   |   %d lines, %zu keywords   |   F7: full check with the parser", ne, nw,
                            nd - size_t(ne + nw), int(d->lines.size()), d->outline.items.size());
        if (ui::beginTable("##probs", 3, ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingFixedFit, ImVec2(0, -FLT_MIN))) {
            ui::tableColumn("", ImGuiTableColumnFlags_WidthFixed, fs * 1.2f);
            ui::tableColumn("Line", ImGuiTableColumnFlags_WidthFixed, fs * 3.5f);
            ui::tableColumn("Message", ImGuiTableColumnFlags_WidthStretch);
            int idx = 0;
            for (const auto* l : {&d->parserDiags, &d->diags}) {
                for (const auto& g : *l) {
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0);
                    ImGui::TextColored(col4(g.severity >= 2 ? Pal::Err : g.severity == 1 ? Pal::Warn : Pal::Accent), "%s", g.severity >= 2 ? "\xE2\x9C\x96" : g.severity == 1 ? "!" : "i");
                    ImGui::TableSetColumnIndex(1);
                    ImGui::Text("%d", g.line + 1);
                    ImGui::TableSetColumnIndex(2);
                    ImGui::PushID(idx++);
                    if (ImGui::Selectable(g.text.c_str(), false, ImGuiSelectableFlags_SpanAllColumns)) jumpLine = g.line;
                    ImGui::PopID();
                }
            }
            ui::endTable();
        }
        ImGui::EndChild();
    }

    if (jumpLine >= 0) {
        d->ed.SetCursor(TextEditor::DocPos(size_t(jumpLine), 0));
        d->ed.ScrollToLine(size_t(jumpLine), TextEditor::Scroll::alignMiddle);
        d->ed.SetFocus();
    }
    ImGui::End();
}

bool App::edConfirmClose() {
    if (!ed_) return true;
    for (size_t i = 0; i < ed_->docs.size(); ++i) {
        EdDoc& x = *ed_->docs[i];
        if (!x.modified()) continue;
        const int r = plat::msgBox(("Save changes to " + x.title + " (Deck Editor)?").c_str(), "ResSim Studio", MB_YESNOCANCEL | MB_ICONQUESTION);
        if (r == IDCANCEL) return false;
        if (r == IDYES && !edSave(int(i), false)) return false;
    }
    return true;
}
