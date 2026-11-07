// C39 platform layer (see platform.h)
#include "platform.h"
#include "gfx.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <vector>

namespace fs = std::filesystem;

std::string safeFileName(const std::string& name);   // app.cpp

bool fileDialog(bool save, const char* filter, const char* defExt, std::string& path) { return plat::fileDialog(save, filter, defExt, path); }

#ifdef _WIN32
// ============================================================================ Windows
#include <commdlg.h>
#include <shellapi.h>
#include <shlobj.h>

namespace plat {

static std::wstring widen(const std::string& s) {
    if (s.empty()) return std::wstring();
    const int n = MultiByteToWideChar(CP_UTF8, 0, s.data(), int(s.size()), nullptr, 0);
    std::wstring w(size_t(std::max(0, n)), L'\0');
    if (n > 0) MultiByteToWideChar(CP_UTF8, 0, s.data(), int(s.size()), &w[0], n);
    return w;
}

int msgBox(const std::string& text, const std::string& caption, unsigned flags) {
    return MessageBoxW(g_gfx.hwnd, widen(text).c_str(), widen(caption).c_str(), flags);
}

// Common open/save dialog. `path` may be a bare file name, a full path, or empty. Only the file-name part is
// sanitised (a drive letter's ':' stays). If the dialog refuses the initial name (GetSaveFileName fails silently
// with an invalid lpstrFile) it is reopened without one.
bool fileDialog(bool save, const char* filter, const char* defExt, std::string& path) {
    std::string dir, name = path;
    const size_t slash = path.find_last_of("\\/");
    if (slash != std::string::npos) { dir = path.substr(0, slash); name = path.substr(slash + 1); }
    if (!name.empty()) name = safeFileName(name);
    if (!dir.empty()) {
        const DWORD attr = GetFileAttributesA(dir.c_str());
        if (attr == INVALID_FILE_ATTRIBUTES || !(attr & FILE_ATTRIBUTE_DIRECTORY)) dir.clear();
    }
    for (int attempt = 0; attempt < 2; ++attempt) {
        std::vector<char> buf(32768, 0);
        if (attempt == 0 && !name.empty() && name.size() < buf.size()) std::memcpy(buf.data(), name.data(), name.size());
        OPENFILENAMEA ofn{};
        ofn.lStructSize = sizeof(ofn);
        ofn.hwndOwner = g_gfx.hwnd;
        ofn.lpstrFilter = filter;
        ofn.nFilterIndex = 1;
        ofn.lpstrFile = buf.data();
        ofn.nMaxFile = DWORD(buf.size());
        ofn.lpstrInitialDir = dir.empty() ? nullptr : dir.c_str();
        ofn.lpstrDefExt = defExt;
        ofn.Flags = OFN_EXPLORER | OFN_NOCHANGEDIR | OFN_HIDEREADONLY | OFN_PATHMUSTEXIST |
                    (save ? OFN_OVERWRITEPROMPT : OFN_FILEMUSTEXIST);
        const BOOL ok = save ? GetSaveFileNameA(&ofn) : GetOpenFileNameA(&ofn);
        if (ok) { path = buf.data(); return true; }
        const DWORD err = CommDlgExtendedError();
        if (err == 0) return false;                            // the user cancelled
        if (attempt == 0) { dir.clear(); continue; }           // e.g. FNERR_INVALIDFILENAME (0x3002): retry with no name
        char msg[160];
        std::snprintf(msg, sizeof(msg), "The file dialog could not be opened (error 0x%04lX).", (unsigned long)err);
        MessageBoxA(g_gfx.hwnd, msg, "ResSim Studio", MB_OK | MB_ICONWARNING);
    }
    return false;
}

// C42 multi-selection open dialog: "dir\0file1\0file2\0\0" (several) or "fullpath\0\0" (one)
bool fileDialogMulti(const char* filter, std::vector<std::string>& paths) {
    paths.clear();
    std::vector<char> buf(1 << 18, 0);
    OPENFILENAMEA ofn{};
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = g_gfx.hwnd;
    ofn.lpstrFilter = filter;
    ofn.nFilterIndex = 1;
    ofn.lpstrFile = buf.data();
    ofn.nMaxFile = DWORD(buf.size());
    ofn.Flags = OFN_EXPLORER | OFN_NOCHANGEDIR | OFN_HIDEREADONLY | OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST | OFN_ALLOWMULTISELECT;
    if (!GetOpenFileNameA(&ofn)) return false;
    const char* p = buf.data();
    const std::string first = p;
    p += first.size() + 1;
    if (!*p) { paths.push_back(first); return true; }
    while (*p) { const std::string f = p; paths.push_back(first + "\\" + f); p += f.size() + 1; }
    return !paths.empty();
}

std::string folderDialog(const std::string& title, const std::string& initial) {
    (void)initial;
    BROWSEINFOW bi{};
    const std::wstring t = widen(title);
    bi.hwndOwner = g_gfx.hwnd; bi.lpszTitle = t.c_str(); bi.ulFlags = BIF_RETURNONLYFSDIRS | BIF_NEWDIALOGSTYLE;
    PIDLIST_ABSOLUTE id = SHBrowseForFolderW(&bi);
    if (!id) return std::string();
    wchar_t buf[MAX_PATH] = {0};
    const bool ok = SHGetPathFromIDListW(id, buf);
    CoTaskMemFree(id);
    return ok ? fs::path(buf).u8string() : std::string();
}

void openPath(const std::string& p) { ShellExecuteW(nullptr, L"open", widen(p).c_str(), nullptr, nullptr, SW_SHOWNORMAL); }

std::string exePath() {
    wchar_t buf[MAX_PATH] = {0};
    GetModuleFileNameW(nullptr, buf, MAX_PATH);
    return fs::path(buf).u8string();
}

std::string configDir() {
    wchar_t buf[MAX_PATH];
    if (SUCCEEDED(SHGetFolderPathW(nullptr, CSIDL_APPDATA, nullptr, 0, buf))) {
        fs::path d = fs::path(buf) / L"ResSimStudio";
        std::error_code ec; fs::create_directories(d, ec);
        return d.u8string();
    }
    return std::string(".");
}

std::string getEnv(const char* name) {
    std::vector<wchar_t> buf(32768);
    const DWORD n = GetEnvironmentVariableW(widen(name).c_str(), buf.data(), DWORD(buf.size()));
    if (n == 0 || n >= buf.size()) return std::string();
    return fs::path(std::wstring(buf.data(), n)).u8string();
}
void setEnv(const char* name, const std::string& value) { SetEnvironmentVariableW(widen(name).c_str(), widen(value).c_str()); }
int screenWidth() { return GetSystemMetrics(SM_CXSCREEN); }
void setWindowTitle(const std::string& t) { SetWindowTextW(g_gfx.hwnd, widen(t).c_str()); }
const char* name() { return "Windows / Direct3D 11 + Direct2D"; }
const char* exeSuffix() { return ".exe"; }
long processId() { return long(GetCurrentProcessId()); }

}  // namespace plat

#else
// ============================================================================ Linux / macOS
#include "tinyfiledialogs.h"
#include <GLFW/glfw3.h>
#include <unistd.h>
#if defined(__APPLE__)
#include <mach-o/dyld.h>
#endif

namespace plat {

int msgBox(const std::string& text, const std::string& caption, unsigned flags) {
    const unsigned btn = flags & 0xFu, icon = flags & 0xF0u;
    const char* type = btn == MB_OKCANCEL ? "okcancel" : btn == MB_YESNO ? "yesno" : btn == MB_YESNOCANCEL ? "yesnocancel" : "ok";
    const char* ic = icon == MB_ICONERROR ? "error" : icon == MB_ICONWARNING ? "warning" : icon == MB_ICONQUESTION ? "question" : "info";
    // tinyfiledialogs rejects quotes in some back ends: replace them
    std::string t = text, c = caption;
    for (std::string* s : {&t, &c}) for (char& ch : *s) if (ch == '"' || ch == '\'' || ch == '`') ch = ' ';
    const int r = tinyfd_messageBox(c.c_str(), t.c_str(), type, ic, 1);
    if (btn == MB_OKCANCEL) return r == 1 ? IDOK : IDCANCEL;
    if (btn == MB_YESNO) return r == 1 ? IDYES : IDNO;
    if (btn == MB_YESNOCANCEL) return r == 1 ? IDYES : r == 2 ? IDNO : IDCANCEL;
    return IDOK;
}

// first filter group of a Win32 filter string -> description + patterns (both letter cases, case-sensitive dialogs)
static void parseFilter(const char* filter, std::string& desc, std::vector<std::string>& pats) {
    if (!filter || !*filter) return;
    desc = filter;
    const char* p = filter + desc.size() + 1;
    std::string spec = p;
    size_t a = 0;
    while (a <= spec.size()) {
        size_t b = spec.find(';', a);
        if (b == std::string::npos) b = spec.size();
        std::string s = spec.substr(a, b - a);
        while (!s.empty() && s.front() == ' ') s.erase(s.begin());
        if (!s.empty() && s != "*.*" && s != "*") {
            std::string lo = s, up = s;
            for (char& ch : lo) ch = char(std::tolower((unsigned char)ch));
            for (char& ch : up) ch = char(std::toupper((unsigned char)ch));
            for (const std::string& x : {s, lo, up}) if (std::find(pats.begin(), pats.end(), x) == pats.end()) pats.push_back(x);
        }
        a = b + 1;
    }
}

bool fileDialog(bool save, const char* filter, const char* defExt, std::string& path) {
    std::string dir, name = path;
    const size_t slash = path.find_last_of("\\/");
    if (slash != std::string::npos) { dir = path.substr(0, slash); name = path.substr(slash + 1); }
    if (!name.empty()) name = safeFileName(name);
    std::error_code ec;
    if (dir.empty() || !fs::is_directory(fs::u8path(dir), ec)) dir = fs::current_path(ec).u8string();
    std::string init = dir + "/" + name;
    std::string desc; std::vector<std::string> pats;
    parseFilter(filter, desc, pats);
    std::vector<const char*> pp; for (auto& s : pats) pp.push_back(s.c_str());
    const char* r = save ? tinyfd_saveFileDialog("Save as", init.c_str(), int(pp.size()), pp.empty() ? nullptr : pp.data(), desc.empty() ? nullptr : desc.c_str())
                         : tinyfd_openFileDialog("Open", init.c_str(), int(pp.size()), pp.empty() ? nullptr : pp.data(), desc.empty() ? nullptr : desc.c_str(), 0);
    if (!r || !*r) return false;
    path = r;
    if (save && defExt && *defExt && fs::u8path(path).extension().empty()) path += std::string(".") + defExt;
    return true;
}

bool fileDialogMulti(const char* filter, std::vector<std::string>& paths) {
    paths.clear();
    std::error_code ec;
    const std::string init = fs::current_path(ec).u8string() + "/";
    std::string desc; std::vector<std::string> pats;
    parseFilter(filter, desc, pats);
    std::vector<const char*> pp; for (auto& s : pats) pp.push_back(s.c_str());
    const char* r = tinyfd_openFileDialog("Open (several files allowed)", init.c_str(), int(pp.size()), pp.empty() ? nullptr : pp.data(), desc.empty() ? nullptr : desc.c_str(), 1);
    if (!r || !*r) return false;
    std::string s = r, cur;
    for (char c : s) { if (c == '|') { if (!cur.empty()) paths.push_back(cur); cur.clear(); } else cur += c; }
    if (!cur.empty()) paths.push_back(cur);
    return !paths.empty();
}

std::string folderDialog(const std::string& title, const std::string& initial) {
    const char* r = tinyfd_selectFolderDialog(title.c_str(), initial.empty() ? nullptr : initial.c_str());
    return r ? std::string(r) : std::string();
}

void openPath(const std::string& p) {
    std::string q = p;
    for (char& ch : q) if (ch == '\'') ch = '_';
#if defined(__APPLE__)
    const std::string cmd = "open '" + q + "' >/dev/null 2>&1 &";
#else
    const std::string cmd = "xdg-open '" + q + "' >/dev/null 2>&1 &";
#endif
    if (std::system(cmd.c_str()) != 0) std::fprintf(stderr, "could not open %s\n", p.c_str());
}

std::string exePath() {
#if defined(__APPLE__)
    char buf[4096]; uint32_t n = sizeof(buf);
    if (_NSGetExecutablePath(buf, &n) == 0) { std::error_code ec; return fs::weakly_canonical(fs::path(buf), ec).string(); }
    return std::string();
#else
    char buf[4096];
    const ssize_t n = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    return n > 0 ? std::string(buf, size_t(n)) : std::string();
#endif
}

std::string configDir() {
    fs::path d;
    const char* home = std::getenv("HOME");
#if defined(__APPLE__)
    d = fs::path(home ? home : ".") / "Library" / "Application Support" / "ResSimStudio";
#else
    const char* x = std::getenv("XDG_CONFIG_HOME");
    d = (x && *x) ? fs::path(x) / "ResSimStudio" : fs::path(home ? home : ".") / ".config" / "ResSimStudio";
#endif
    std::error_code ec; fs::create_directories(d, ec);
    return d.string();
}

std::string getEnv(const char* name) { const char* v = std::getenv(name); return v ? std::string(v) : std::string(); }
void setEnv(const char* name, const std::string& value) { setenv(name, value.c_str(), 1); }
int screenWidth() {
    if (GLFWmonitor* m = glfwGetPrimaryMonitor()) if (const GLFWvidmode* v = glfwGetVideoMode(m)) return v->width;
    return 0;
}
void setWindowTitle(const std::string& t) { if (g_gfx.window) glfwSetWindowTitle(g_gfx.window, t.c_str()); }
#if defined(__APPLE__)
const char* name() { return "macOS / OpenGL 3.3"; }
#else
const char* name() { return "Linux / OpenGL 3.3"; }
#endif
const char* exeSuffix() { return ""; }
long processId() { return long(getpid()); }

}  // namespace plat
#endif

std::string plat::exeDir() {
    const std::string p = exePath();
    const size_t s = p.find_last_of("\\/");
    return s == std::string::npos ? std::string() : p.substr(0, s + 1);
}
