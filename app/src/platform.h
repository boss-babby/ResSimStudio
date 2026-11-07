// C39 platform layer: everything the UI needs from the operating system, behind one portable interface.
//   Windows       : Win32 (common dialogs, MessageBox, ShellExecute, SHGetFolderPath) - the original behaviour
//   Linux / macOS : tinyfiledialogs (zenity / kdialog / osascript), xdg-open / open, $XDG_CONFIG_HOME / ~/Library
// Message-box flags and results use the Win32 values (MB_OK ..., IDOK ...) on every platform, so call sites stay the same.
#pragma once
#include <string>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
// Win32-compatible message box flags / results
#define MB_OK 0x0u
#define MB_OKCANCEL 0x1u
#define MB_YESNOCANCEL 0x3u
#define MB_YESNO 0x4u
#define MB_ICONERROR 0x10u
#define MB_ICONQUESTION 0x20u
#define MB_ICONWARNING 0x30u
#define MB_ICONINFORMATION 0x40u
#define IDOK 1
#define IDCANCEL 2
#define IDYES 6
#define IDNO 7
#endif

namespace plat {

// modal message box (Win32 semantics: flags MB_*, returns IDOK / IDCANCEL / IDYES / IDNO)
int msgBox(const std::string& text, const std::string& caption, unsigned flags);
// open / save dialog; filter in Win32 form ("Decks\0*.DATA;*.txt\0All\0*.*\0\0"); path: in = initial name / folder,
// out = chosen file (UTF-8). Returns false when cancelled.
bool fileDialog(bool save, const char* filter, const char* defExt, std::string& path);
// C42: open dialog with multiple selection (UTF-8 paths); false when cancelled
bool fileDialogMulti(const char* filter, std::vector<std::string>& paths);
// folder picker (UTF-8), empty = cancelled
std::string folderDialog(const std::string& title, const std::string& initial = std::string());
// open a file, folder or URL with the desktop's default handler
void openPath(const std::string& pathOrUrl);
// absolute path of the running executable and its folder (with trailing separator)
std::string exePath();
std::string exeDir();
// per-user settings folder (created): %APPDATA%\ResSimStudio, ~/.config/ResSimStudio, ~/Library/Application Support/ResSimStudio
std::string configDir();
// environment (UTF-8)
std::string getEnv(const char* name);
void setEnv(const char* name, const std::string& value);
// primary screen width in physical pixels (0 if unknown)
int screenWidth();
// main window title
void setWindowTitle(const std::string& title);
// name of the platform / graphics back end, e.g. "Windows / Direct3D 11", "Linux / OpenGL 3.3"
const char* name();
// default executable suffix (".exe" on Windows, "" elsewhere)
const char* exeSuffix();
long processId();

}  // namespace plat

// legacy entry point used throughout the UI (same as plat::fileDialog)
bool fileDialog(bool save, const char* filter, const char* defExt, std::string& path);
