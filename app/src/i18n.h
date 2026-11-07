// C41 (G13) interface translation.
//
// TR("text") returns the translation of a UI string in the selected language, or the string itself. Catalogues are
// plain UTF-8 files "<exe folder>/lang/<code>.txt" (also <exe>/../../lang for development builds):
//     # LANGUAGE <code> <native name>
//     msgid<TAB>msgstr        (\t, \n escapes; '#' comment lines)
// Leading icon glyphs (Unicode private-use area) and spaces, and an ImGui "##id" / "###id" suffix, are kept and
// not looked up, so `TR(ICON_SAVE "  Save##x")` uses the catalogue entry "Save" and the ImGui ID stays stable.
// Results are cached per string pointer: TR must be given string literals (or other storage that outlives a frame).
// Limitation: Dear ImGui draws one glyph per code point without complex-script shaping, so scripts that need
// shaping (Bengali, Devanagari, Arabic) cannot be displayed correctly; only shaping-free scripts are offered.
#pragma once
#include <string>
#include <vector>

namespace i18n {
struct Language { std::string code, name, path; };
std::vector<Language> available();               // "en" (built in) + catalogues found on disk
bool setLanguage(const std::string& code);       // false: catalogue not found (falls back to English)
const std::string& language();
const char* tr(const char* s);
// C45: content-keyed variant for strings that may live in reused buffers (std::string::c_str() of temporaries, rows of
// tables): the cache is keyed by the text, so equal pointers with different contents never return a stale entry
const char* trs(const char* s);
// C45: a widget label translated with a stable ImGui ID - "Translated###English" when the label has no "##" id of its own
std::string trLabel(const char* s);
int entries();                                   // entries of the loaded catalogue
// parse a catalogue text into pairs (tests / tools); returns the number of entries, err set on malformed lines
int parseCatalogue(const std::string& text, std::vector<std::pair<std::string, std::string>>& out, std::string& err, std::string* name = nullptr);
}  // namespace i18n
#define TR(s) ::i18n::tr(s)
#define TRS(s) ::i18n::trs(s)
#define TRL(s) ::i18n::trLabel(s).c_str()   // widget label: translated text, English ImGui id
#define TRN(s) s   // marks a literal for tools/i18n_extract.py; translate it later with TR(variable)
