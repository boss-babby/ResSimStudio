// Generic ECLIPSE / OPM Flow deck reader.
//
// Stage 1 of deck import: turns deck text (with INCLUDE / PATHS resolved) into a list of keywords with their
// records, using the record-structure rules of ~2900 ECLIPSE 100 / 300 / OPM keywords (generated from OPM's
// keyword metadata, see tools/gen_keywords.py). Every keyword - supported by ResSim or not - is consumed with
// the correct number of records, so unknown keywords can never derail the parse. Stage 2 (deck.cpp) maps the
// keywords onto a ressim::Case.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace ressim {
namespace deck {

enum class KwSize : uint8_t {
    None,     // no data (OIL, WATER, ...)
    Fixed,    // a fixed number of records
    FromKw,   // number of records given by an item of another keyword (e.g. SWOF: TABDIMS item 1)
    Slash,    // records until an empty record "/"
    Unknown,  // records until the next keyword (VFPPROD, ...)
    Tables,   // N tables, each a list of records closed by an empty record (PVTO, PVTG)
    Data,     // one record holding an array (PORO, PERMX, ...)
    Code,     // raw text lines until an end marker (PYINPUT)
    Action,   // ACTIONX block: raw lines until ENDACTIO
    Raw       // one raw text line (TITLE)
};

struct KeywordSpec {
    const char* name;
    uint16_t sections;       // bit mask: RUNSPEC GRID EDIT PROPS REGIONS SOLUTION SUMMARY SCHEDULE
    KwSize kind;
    int16_t count;           // Fixed
    const char* sizeKw;      // FromKw / Tables
    int8_t sizeItem;         // 0-based item of sizeKw
    int16_t sizeDefault;
    int8_t shift;
    const char* end;         // Code / Action end marker
};

const KeywordSpec* findKeyword(const std::string& name);     // exact names and summary-style families
const std::vector<KeywordSpec>& keywordTable();               // all known keywords (sorted by name)
const char* sectionName(int bit);                             // 0..7

// How ResSim treats a keyword (deck viewer colouring, keyword reference, import report).
enum class Support : uint8_t { Full, Partial, Ignored, Unsupported, Unknown };
Support ressimSupport(const std::string& keyword, std::string* note = nullptr);
const char* toString(Support s);

struct Item {
    std::string s;
    bool def = false;        // defaulted (1* / trailing omission)
    bool quoted = false;
};
struct Record {
    std::vector<Item> items;
    int line = 0;
    bool terminator = false; // empty "/" closing a table (Tables kind)
    bool has(size_t i) const { return i < items.size() && !items[i].def; }
    std::string str(size_t i, const std::string& def = "") const { return has(i) ? items[i].s : def; }
    bool isNumber(size_t i) const;
    double num(size_t i, double def) const;            // throws on a non-numeric item
    int integer(size_t i, int def) const { return int(num(i, def) + (num(i, def) < 0 ? -0.5 : 0.5)); }
};
struct Keyword {
    std::string name;
    int file = 0;            // index into RawDeck::files
    int line = 0;
    std::vector<Record> records;
    std::string raw;         // Raw / Code / Action text
    std::vector<double> values;   // Data keywords (PORO, ZCORN, ...): numbers with repeats expanded, NaN = defaulted
    bool isArray = false;
    bool known = true;       // found in the keyword table
};

enum class Level : uint8_t { Info, Warning, Error };
struct Message {
    Level level = Level::Warning;
    std::string keyword;
    int file = 0, line = 0;
    std::string text;
};

struct RawDeck {
    std::vector<std::string> files;      // [0] = main deck
    std::vector<Keyword> keywords;
    std::vector<Message> messages;
    std::string where(int file, int line) const;   // "file:line"
    void message(Level l, const Keyword& k, const std::string& text, int line = 0);
};

// Parse deck text. baseDir is used to resolve relative INCLUDE paths (with '/' at the end or empty).
// Never throws for syntax problems - they are reported in RawDeck::messages (Error level) instead.
// Throws std::runtime_error only when the main file cannot be read.
RawDeck parseDeckText(const std::string& text, const std::string& baseDir, const std::string& fileName = "deck");
RawDeck parseDeckFile(const std::string& path);

// Deck text of a parsed keyword's data (records, arrays, raw blocks with their end marker), without the keyword
// line. Parsing "NAME\n" + keywordBodyText(k) gives the same records again.
std::string keywordBodyText(const Keyword& k);

// Lightweight line classifier for syntax highlighting (no INCLUDE expansion).
enum class Tok : uint8_t { Space, Comment, Keyword, KeywordUnknown, Section, Number, String, Default, Slash, Text, RessimExt };
struct TokSpan { int begin, end; Tok kind; };
// 'state' carries multi-line context (TITLE line, ACTIONX/PYINPUT blocks) between consecutive calls.
struct LineState { int raw = 0; std::string end; };
void classifyLine(const std::string& line, std::vector<TokSpan>& out, LineState& state);

}  // namespace deck
}  // namespace ressim
