// C32 ResGeo.exe: the standalone geomodelling front end. Starts ResSimStudio.exe (same folder) with --geo, which opens
// straight into the Geomodel workspace with the other workspaces hidden (window title "ResGeo"); arguments are forwarded.
#include <windows.h>
#include <string>

int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR cmdLine, int) {
    char exe[MAX_PATH] = {0};
    GetModuleFileNameA(nullptr, exe, MAX_PATH);
    std::string dir(exe);
    const size_t sl = dir.find_last_of("\\/");
    dir = sl == std::string::npos ? std::string() : dir.substr(0, sl + 1);
    const std::string studio = dir + "ResSimStudio.exe";
    std::string cmd = "\"" + studio + "\" --geo";
    if (cmdLine && *cmdLine) cmd += std::string(" ") + cmdLine;
    STARTUPINFOA si{}; si.cb = sizeof(si);
    PROCESS_INFORMATION pi{};
    if (!CreateProcessA(studio.c_str(), &cmd[0], nullptr, nullptr, FALSE, 0, nullptr, dir.empty() ? nullptr : dir.c_str(), &si, &pi)) {
        MessageBoxA(nullptr, ("ResGeo needs ResSimStudio.exe in the same folder:\n" + studio).c_str(), "ResGeo", MB_ICONERROR);
        return 1;
    }
    CloseHandle(pi.hThread); CloseHandle(pi.hProcess);
    return 0;
}
