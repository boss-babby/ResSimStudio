// C41 (G13) interface translation - see i18n.h
#include "i18n.h"
#include "platform.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <cstring>

namespace fs = std::filesystem;

namespace i18n {
namespace {
std::string g_lang = "en";
std::unordered_map<std::string, std::string> g_cat;
std::unordered_map<const char*, std::string> g_cache;
std::unordered_map<std::string, std::string> g_cacheS;   // C45 content-keyed cache (trs)

std::string unescape(const std::string& s) {
    std::string o;
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '\\' && i + 1 < s.size()) {
            const char c = s[++i];
            o += c == 'n' ? '\n' : c == 't' ? '\t' : c;
        } else o += s[i];
    }
    return o;
}
std::vector<fs::path> dirs() {
    const fs::path exe = fs::u8path(plat::exePath()).parent_path();
    return {exe / "lang", exe.parent_path().parent_path() / "lang"};
}
// length of the leading icon glyphs (U+E000..U+F8FF, 3-byte UTF-8 EE / EF 80..A3) and spaces
size_t prefixLen(const char* s) {
    size_t i = 0;
    for (;;) {
        const unsigned char c = (unsigned char)s[i];
        if (c == ' ') { ++i; continue; }
        if ((c == 0xEE || (c == 0xEF && (unsigned char)s[i + 1] <= 0xA3)) && s[i + 1] && s[i + 2]) { i += 3; continue; }
        break;
    }
    return i;
}
}  // namespace

int parseCatalogue(const std::string& text, std::vector<std::pair<std::string, std::string>>& out, std::string& err, std::string* name) {
    out.clear(); err.clear();
    std::istringstream is(text);
    std::string line;
    int ln = 0;
    while (std::getline(is, line)) {
        ++ln;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (ln == 1 && line.size() >= 3 && (unsigned char)line[0] == 0xEF && (unsigned char)line[1] == 0xBB && (unsigned char)line[2] == 0xBF) line.erase(0, 3);
        if (line.rfind("# LANGUAGE ", 0) == 0) {
            if (name) { std::istringstream ls(line.substr(11)); std::string code; ls >> code; std::getline(ls, *name); while (!name->empty() && (*name)[0] == ' ') name->erase(0, 1); }
            continue;
        }
        if (line.empty() || line[0] == '#') continue;
        const size_t t = line.find('\t');
        if (t == std::string::npos || t == 0) { if (err.empty()) err = "line " + std::to_string(ln) + ": expected msgid<TAB>msgstr"; continue; }
        const std::string id = unescape(line.substr(0, t)), str = unescape(line.substr(t + 1));
        if (!str.empty()) out.emplace_back(id, str);
    }
    return int(out.size());
}

std::vector<Language> available() {
    std::vector<Language> v{{"en", "English", ""}};
    std::error_code ec;
    for (const fs::path& d : dirs()) {
        if (!fs::is_directory(d, ec)) continue;
        for (const auto& e : fs::directory_iterator(d, ec)) {
            if (e.path().extension() != ".txt") continue;
            const std::string code = e.path().stem().u8string();
            bool dup = false; for (const auto& l : v) dup = dup || l.code == code;
            if (dup) continue;
            std::ifstream f(e.path(), std::ios::binary);
            std::string head; std::getline(f, head);
            std::string name = code;
            std::vector<std::pair<std::string, std::string>> tmp; std::string err;
            parseCatalogue(head + "\n", tmp, err, &name);
            v.push_back({code, name, e.path().u8string()});
        }
    }
    return v;
}

bool setLanguage(const std::string& code) {
    g_cat.clear(); g_cache.clear(); g_cacheS.clear(); g_lang = "en";
    if (code.empty() || code == "en") return true;
    for (const Language& l : available()) {
        if (l.code != code || l.path.empty()) continue;
        std::ifstream f(fs::u8path(l.path), std::ios::binary);
        std::stringstream ss; ss << f.rdbuf();
        std::vector<std::pair<std::string, std::string>> p; std::string err;
        parseCatalogue(ss.str(), p, err);
        for (auto& kv : p) g_cat[kv.first] = kv.second;
        g_lang = code;
        return true;
    }
    return false;
}

const std::string& language() { return g_lang; }
int entries() { return int(g_cat.size()); }

namespace {
std::string translate(const std::string& all) {
    const size_t p = prefixLen(all.c_str());
    size_t e = all.find("##", p);
    if (e == std::string::npos) e = all.size();
    std::string mid = all.substr(p, e - p);
    size_t tail = mid.size(); while (tail > 0 && mid[tail - 1] == ' ') --tail;
    auto it = g_cat.find(mid.substr(0, tail));
    return it == g_cat.end() ? all : all.substr(0, p) + it->second + mid.substr(tail) + all.substr(e);
}
}  // namespace
const char* trs(const char* s) {
    if (!s || g_cat.empty()) return s;
    std::string k(s);
    auto c = g_cacheS.find(k);
    if (c != g_cacheS.end()) return c->second.c_str();
    if (g_cacheS.size() > 20000) g_cacheS.clear();   // unbounded dynamic labels (cell values, file names) must not grow it forever
    std::string out = translate(k);
    return g_cacheS.emplace(std::move(k), std::move(out)).first->second.c_str();
}
std::string trLabel(const char* s) {
    if (!s) return std::string();
    const char* t = trs(s);
    if (t == s || std::strcmp(t, s) == 0 || std::strstr(s, "##")) return std::string(t);
    return std::string(t) + "###" + s;
}
const char* tr(const char* s) {
    if (!s || g_cat.empty()) return s;
    auto c = g_cache.find(s);
    if (c != g_cache.end()) return c->second.c_str();
    const std::string all(s);
    const size_t p = prefixLen(s);
    size_t e = all.find("##", p);
    if (e == std::string::npos) e = all.size();
    std::string mid = all.substr(p, e - p);
    size_t tail = mid.size(); while (tail > 0 && mid[tail - 1] == ' ') --tail;   // trailing spaces kept
    const std::string key = mid.substr(0, tail);
    auto it = g_cat.find(key);
    std::string out = it == g_cat.end() ? all : all.substr(0, p) + it->second + mid.substr(tail) + all.substr(e);
    return g_cache.emplace(s, std::move(out)).first->second.c_str();
}
}  // namespace i18n
