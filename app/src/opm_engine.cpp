#include "opm_engine.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <thread>

#ifdef _WIN32
#include <windows.h>
#else
#include <csignal>
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace fs = std::filesystem;
using namespace ressim;

namespace {
bool isFile(const fs::path& p) { std::error_code ec; return fs::is_regular_file(p, ec); }
bool isDir(const fs::path& p) { std::error_code ec; return fs::is_directory(p, ec); }
std::string envStr(const char* name) {
#ifdef _WIN32
    wchar_t buf[32768];
    std::wstring wn(name, name + std::strlen(name));
    DWORD n = GetEnvironmentVariableW(wn.c_str(), buf, 32768);
    if (n == 0 || n >= 32768) return {};
    return fs::path(std::wstring(buf, n)).u8string();
#else
    const char* v = std::getenv(name); return v ? v : "";
#endif
}
fs::path u8(const std::string& s) { return fs::u8path(s); }

#ifdef _WIN32
const char* kExeNames[] = {"flow.exe", "flow_blackoil.exe"};
const char kPathSep = ';';
#else
const char* kExeNames[] = {"flow", "flow_blackoil"};
const char kPathSep = ':';
#endif

// flow executable inside a folder (folder itself, bin\, and one level of sub-folders, as in the portable zip)
std::string findInFolder(const fs::path& d, int depth) {
    if (!isDir(d)) return {};
    for (const char* n : kExeNames) {
        if (isFile(d / n)) return (d / n).u8string();
        if (isFile(d / "bin" / n)) return (d / "bin" / n).u8string();
    }
    if (depth <= 0) return {};
    std::error_code ec;
    for (auto& e : fs::directory_iterator(d, ec)) {
        if (!e.is_directory(ec)) continue;
        std::string nm = e.path().filename().u8string();
        std::transform(nm.begin(), nm.end(), nm.begin(), [](unsigned char c) { return char(std::tolower(c)); });
        if (nm.rfind("opm", 0) == 0 || nm.find("flow") != std::string::npos) {
            std::string r = findInFolder(e.path(), depth - 1);
            if (!r.empty()) return r;
        }
    }
    return {};
}
} // namespace

std::string detectFlowExecutable(const std::string& hint, std::vector<std::string>* searched) {
    auto note = [&](const std::string& s) { if (searched) searched->push_back(s); };
    if (!hint.empty()) {
        note(hint);
        if (isFile(u8(hint))) return hint;
        std::string r = findInFolder(u8(hint), 2);
        if (!r.empty()) return r;
    }
    for (const char* ev : {"FLOW_EXE", "OPM_FLOW"}) {
        std::string v = envStr(ev);
        if (v.empty()) continue;
        note(std::string("%") + ev + "%");
        if (isFile(u8(v))) return v;
        std::string r = findInFolder(u8(v), 1);
        if (!r.empty()) return r;
    }
    note("PATH");
    {
        const std::string path = envStr("PATH");
        size_t a = 0;
        while (a <= path.size()) {
            size_t b = path.find(kPathSep, a);
            if (b == std::string::npos) b = path.size();
            std::string dir = path.substr(a, b - a);
            if (!dir.empty() && dir.front() == '"' && dir.back() == '"' && dir.size() >= 2) dir = dir.substr(1, dir.size() - 2);
            if (!dir.empty())
                for (const char* n : kExeNames) if (isFile(u8(dir) / n)) return (u8(dir) / n).u8string();
            a = b + 1;
        }
    }
    std::vector<fs::path> roots;
#ifdef _WIN32
    for (const char* ev : {"ProgramW6432", "ProgramFiles", "ProgramFiles(x86)"}) { std::string v = envStr(ev); if (!v.empty()) roots.push_back(u8(v)); }
    { std::string v = envStr("LOCALAPPDATA"); if (!v.empty()) { roots.push_back(u8(v) / "Programs"); roots.push_back(u8(v)); } }
    { std::string v = envStr("USERPROFILE"); if (!v.empty()) { roots.push_back(u8(v)); roots.push_back(u8(v) / "Downloads"); roots.push_back(u8(v) / "Desktop"); } }
    roots.push_back("C:\\");
    roots.push_back("C:\\OPM");
#else
    roots = {"/usr/bin", "/usr/local/bin", "/opt"};
    { std::string v = envStr("HOME"); if (!v.empty()) roots.push_back(u8(v)); }
#endif
    for (auto& r : roots) {
        note(r.u8string());
        // direct "OPM Flow" folder (installer), otherwise OPM* / *flow* sub-folders up to two levels deep
        std::string f = findInFolder(r / "OPM Flow", 1);
        if (f.empty()) {
            std::error_code ec;
            if (!isDir(r)) continue;
            for (auto& e : fs::directory_iterator(r, ec)) {
                if (!e.is_directory(ec)) continue;
                std::string nm = e.path().filename().u8string();
                std::transform(nm.begin(), nm.end(), nm.begin(), [](unsigned char c) { return char(std::tolower(c)); });
                if (nm.rfind("opm", 0) != 0) continue;
                f = findInFolder(e.path(), 2);
                if (!f.empty()) break;
            }
        }
        if (!f.empty()) return f;
#ifndef _WIN32
        for (const char* n : kExeNames) if (isFile(r / n)) return (r / n).u8string();
#endif
    }
    return {};
}

std::string flowPrerequisiteProblem(const std::string& flowPath) {
    if (flowPath.empty()) return "OPM Flow executable not set.";
    if (!isFile(u8(flowPath))) return "OPM Flow executable not found: " + flowPath;
#ifdef _WIN32
    wchar_t sys[MAX_PATH]; UINT n = GetSystemDirectoryW(sys, MAX_PATH);
    const fs::path sysDir = n ? fs::path(std::wstring(sys, n)) : fs::path("C:\\Windows\\System32");
    const fs::path binDir = u8(flowPath).parent_path();
    auto onPath = [](const wchar_t* dll) {
        wchar_t out[MAX_PATH]; return SearchPathW(nullptr, dll, nullptr, MAX_PATH, out, nullptr) > 0;
    };
    auto anyDll = [&](const fs::path& d) {
        return isFile(d / "msmpi.dll") || isFile(d / "impi.dll") || isFile(d / "release" / "impi.dll");
    };
    bool mpi = isFile(sysDir / "msmpi.dll") || anyDll(binDir) || anyDll(binDir.parent_path() / "lib") ||
               anyDll(binDir.parent_path() / "redist") || onPath(L"msmpi.dll") || onPath(L"impi.dll");
    const std::string appdata = envStr("APPDATA"), iroot = envStr("I_MPI_ROOT");
    if (!appdata.empty()) mpi = mpi || anyDll(u8(appdata) / "Python" / "Library" / "bin");
    if (!iroot.empty()) mpi = mpi || anyDll(u8(iroot) / "bin") || anyDll(u8(iroot) / "opt" / "mpi" / "libfabric" / "bin");
    // Some builds link MPI statically: only a hint, never an error
    if (!mpi)
        return "Note: no MPI runtime DLL (msmpi.dll / impi.dll) was found next to flow.exe or on PATH. Builds that need one "
               "fail to start without it - then install Microsoft MPI ('winget install Microsoft.msmpi') or Intel MPI.";
#endif
    return {};
}

std::vector<std::string> opmCompatibilityNotes(const Case& c) {
    std::vector<std::string> out;
    const GridInput& g = c.grid;
    if (!g.cornerPoint() && !g.newtran && !g.tops.empty() && !g.dz.empty()) {
        // largest offset between laterally neighbouring columns relative to the thinnest layer
        const int nx = g.nx, ny = g.ny;
        auto top = [&](int i, int j) { return g.tops.size() >= size_t(nx) * ny ? g.tops[size_t(j) * nx + i] : g.tops[0]; };
        double off = 0;
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
            if (i + 1 < nx) off = std::max(off, std::fabs(top(i + 1, j) - top(i, j)));
            if (j + 1 < ny) off = std::max(off, std::fabs(top(i, j + 1) - top(i, j)));
        }
        double dzMin = 1e300;
        for (double v : g.dz) if (v > 0) dzMin = std::min(dzMin, v);
        if (off > 0.05 * dzMin) {
            char b[512];
            std::snprintf(b, sizeof(b),
                          "Block-centred grid with column offsets up to %.3g (thinnest layer %.3g): OPM Flow always connects "
                          "cells through their geometric overlap (NEWTRAN), ResSim currently uses OLDTRAN (i,j,k neighbours). "
                          "Tick 'Geometric transmissibility' on the Grid panel for a like-for-like comparison.", off, dzMin);
            out.push_back(b);
        }
    }
    if (c.tuning.maxDt > 0 && c.tuning.maxDt < 30)
        out.push_back("OPM Flow chooses its own time steps (TUNING is ignored unless 'Use TUNING' is enabled); set the maximum "
                      "time step in the OPM options for a closer comparison of fronts.");
    return out;
}

std::string runDirFor(const std::string& projectPath, const std::string& title, const std::string& tag) {
    if (!projectPath.empty()) {
        fs::path p = u8(projectPath);
        return (p.parent_path() / (p.stem().u8string() + "_" + tag)).u8string();
    }
    std::error_code ec;
    fs::path t = fs::temp_directory_path(ec);
    if (ec) t = ".";
    std::string s;
    for (char ch : title) s += std::isalnum((unsigned char)ch) ? ch : '_';
    if (s.size() > 40) s.resize(40);
    if (s.empty()) s = "case";
    return (t / "ResSimStudio" / tag / s).u8string();
}

std::string caseBaseName(const Case& c) {
    std::string s;
    for (char ch : c.title) {
        if (std::isalnum((unsigned char)ch)) s += char(std::toupper((unsigned char)ch));
        else if (!s.empty() && s.back() != '_') s += '_';
        if (s.size() >= 24) break;
    }
    while (!s.empty() && s.back() == '_') s.pop_back();
    return s.empty() ? std::string("RESSIM_CASE") : s;
}

// ------------------------------------------------------------------------------------------------ progress parsing
bool OpmProgress::parse(const std::string& line) {
    const char* p = line.c_str();
    if (const char* r = std::strstr(p, "Report step")) {
        int a = 0, b = 0; double d = 0, t = 0;
        if (std::sscanf(r, "Report step %d/%d at day %lf/%lf", &a, &b, &d, &t) == 4) {
            reportStep = a; reportSteps = b; day = d; totalDays = t; return true;
        }
    }
    if (std::strstr(p, "Starting time step")) {
        if (const char* r = std::strstr(p, " at day ")) {
            double d = 0;
            if (std::sscanf(r, " at day %lf", &d) == 1 && d >= day) { day = d; return true; }
        }
    }
    int v = 0;
    if (std::sscanf(p, "Number of timesteps: %d", &v) == 1) { timesteps = v; return true; }
    if (std::sscanf(p, "Overall Newton Iterations: %d", &v) == 1) { newton = v; return true; }
    if (std::sscanf(p, "Overall Linear Iterations: %d", &v) == 1) { linear = v; return true; }
    if (std::strstr(p, "End of simulation")) { finished = true; if (totalDays > 0) day = totalDays; return true; }
    // error lines (input problems, unknown options, aborts) are collected for the failure message
    if (std::strncmp(p, "Error", 5) == 0 || std::strstr(p, "Error:") || std::strstr(p, "Aborting") || std::strstr(p, "not supported") ||
        std::strstr(p, "Problem with keyword") || std::strstr(p, "Unrecoverable")) {
        std::string l = line;
        while (!l.empty() && (l[0] == ' ' || l[0] == '\t')) l.erase(0, 1);
        if (lastError.find(l) == std::string::npos && lastError.size() < 600) lastError += (lastError.empty() ? "" : " | ") + l;
        return true;
    }
    return false;
}

// ------------------------------------------------------------------------------------------------ child process
ChildProcess::~ChildProcess() {
    kill();
#ifdef _WIN32
    if (out_) CloseHandle(out_);
    if (thread_) CloseHandle(thread_);
    if (proc_) CloseHandle(proc_);
    if (job_) CloseHandle(job_);
#else
    if (fd_ >= 0) ::close(fd_);
    if (pid_ > 0) { int st; ::waitpid(pid_, &st, 0); }
#endif
}

#ifdef _WIN32
static std::wstring widen(const std::string& s) {
    if (s.empty()) return {};
    int n = MultiByteToWideChar(CP_UTF8, 0, s.data(), int(s.size()), nullptr, 0);
    std::wstring w(size_t(n), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, s.data(), int(s.size()), &w[0], n);
    return w;
}
static std::wstring quoteArg(const std::wstring& a) {
    if (!a.empty() && a.find_first_of(L" \t\"") == std::wstring::npos) return a;
    std::wstring q = L"\"";
    int bs = 0;
    for (wchar_t c : a) {
        if (c == L'\\') { ++bs; continue; }
        if (c == L'"') { q.append(size_t(bs) * 2 + 1, L'\\'); q += L'"'; }
        else { q.append(size_t(bs), L'\\'); q += c; }
        bs = 0;
    }
    q.append(size_t(bs) * 2, L'\\');
    q += L'"';
    return q;
}
#endif

bool ChildProcess::start(const std::string& exe, const std::vector<std::string>& args, const std::string& cwd, std::string* err,
                         const std::string& stdinText) {
    killed_ = false;
#ifdef _WIN32
    SECURITY_ATTRIBUTES sa{sizeof(sa), nullptr, TRUE};
    HANDLE rd = nullptr, wr = nullptr;
    if (!CreatePipe(&rd, &wr, &sa, 0)) { if (err) *err = "CreatePipe failed"; return false; }
    SetHandleInformation(rd, HANDLE_FLAG_INHERIT, 0);
    HANDLE nul = INVALID_HANDLE_VALUE, inWr = nullptr;
    if (!stdinText.empty()) {
        if (!CreatePipe(&nul, &inWr, &sa, 0)) { nul = INVALID_HANDLE_VALUE; inWr = nullptr; }
        else SetHandleInformation(inWr, HANDLE_FLAG_INHERIT, 0);
    }
    if (nul == INVALID_HANDLE_VALUE) nul = CreateFileW(L"NUL", GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE, &sa, OPEN_EXISTING, 0, nullptr);
    STARTUPINFOW si{}; si.cb = sizeof(si);
    si.dwFlags = STARTF_USESTDHANDLES; si.hStdOutput = wr; si.hStdError = wr; si.hStdInput = nul;
    std::wstring cmd = quoteArg(widen(exe));
    for (auto& a : args) { cmd += L' '; cmd += quoteArg(widen(a)); }
    // the flow folder (runtime DLLs) and a pip-installed Intel MPI runtime first on PATH for the child
    {
        std::wstring path;
        wchar_t buf[32767]; DWORD n = GetEnvironmentVariableW(L"PATH", buf, 32767);
        if (n && n < 32767) path.assign(buf, n);
        std::wstring add = fs::u8path(exe).parent_path().wstring();
        const std::string appdata = envStr("APPDATA");
        if (!appdata.empty()) {
            fs::path impi = fs::u8path(appdata) / "Python" / "Library" / "bin";
            if (isDir(impi)) add += L";" + impi.wstring();
        }
        if (path.find(add) == std::wstring::npos) SetEnvironmentVariableW(L"PATH", (add + L";" + path).c_str());
    }
    PROCESS_INFORMATION pi{};
    std::wstring wcwd = widen(cwd);
    const BOOL ok = CreateProcessW(nullptr, &cmd[0], nullptr, nullptr, TRUE, CREATE_NO_WINDOW | BELOW_NORMAL_PRIORITY_CLASS | CREATE_SUSPENDED,
                                   nullptr, wcwd.empty() ? nullptr : wcwd.c_str(), &si, &pi);
    CloseHandle(wr);
    if (nul != INVALID_HANDLE_VALUE) CloseHandle(nul);
    if (!ok) {
        const DWORD e = GetLastError();
        CloseHandle(rd);
        if (inWr) CloseHandle(inWr);
        if (err) { char b[160]; std::snprintf(b, sizeof(b), "could not start the process (Windows error %lu)", (unsigned long)e); *err = b; }
        return false;
    }
    // job object: Stop must also end grandchildren (eclrun starts eclipse.exe, mpiexec starts flow ranks)
    job_ = CreateJobObjectW(nullptr, nullptr);
    if (job_) {
        JOBOBJECT_EXTENDED_LIMIT_INFORMATION li{};
        li.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
        SetInformationJobObject(job_, JobObjectExtendedLimitInformation, &li, sizeof(li));
        if (!AssignProcessToJobObject(job_, pi.hProcess)) { CloseHandle(job_); job_ = nullptr; }
    }
    ResumeThread(pi.hThread);
    if (inWr) {
        DWORD w = 0;
        WriteFile(inWr, stdinText.data(), DWORD(stdinText.size()), &w, nullptr);
        CloseHandle(inWr);
    }
    proc_ = pi.hProcess; thread_ = pi.hThread; out_ = rd;
    return true;
#else
    int p[2], in[2] = {-1, -1};
    if (::pipe(p) != 0) { if (err) *err = "pipe failed"; return false; }
    if (!stdinText.empty() && ::pipe(in) != 0) in[0] = in[1] = -1;
    pid_t pid = ::fork();
    if (pid < 0) { ::close(p[0]); ::close(p[1]); if (in[0] >= 0) { ::close(in[0]); ::close(in[1]); } if (err) *err = "fork failed"; return false; }
    if (pid == 0) {
        ::setpgid(0, 0);
        if (!cwd.empty() && ::chdir(cwd.c_str()) != 0) _exit(126);
        ::dup2(p[1], 1); ::dup2(p[1], 2); ::close(p[0]); ::close(p[1]);
        if (in[0] >= 0) { ::dup2(in[0], 0); ::close(in[0]); ::close(in[1]); }
        else { int nul = ::open("/dev/null", O_RDONLY); if (nul >= 0) { ::dup2(nul, 0); ::close(nul); } }
        std::vector<char*> av; av.push_back(const_cast<char*>(exe.c_str()));
        for (auto& a : args) av.push_back(const_cast<char*>(a.c_str()));
        av.push_back(nullptr);
        ::execvp(exe.c_str(), av.data());
        _exit(127);
    }
    ::close(p[1]);
    if (in[0] >= 0) {
        ::close(in[0]);
        ssize_t w = ::write(in[1], stdinText.data(), stdinText.size()); (void)w;
        ::close(in[1]);
    }
    pid_ = pid; fd_ = p[0];
    return true;
#endif
}

int ChildProcess::pump(const std::function<void(const std::string&)>& onLine) {
    std::string buf;
    auto flush = [&](bool all) {
        size_t a = 0;
        for (;;) {
            size_t b = buf.find('\n', a);
            if (b == std::string::npos) break;
            std::string l = buf.substr(a, b - a);
            if (!l.empty() && l.back() == '\r') l.pop_back();
            if (onLine) onLine(l);
            a = b + 1;
        }
        buf.erase(0, a);
        if (all && !buf.empty()) { if (onLine) onLine(buf); buf.clear(); }
    };
    char tmp[4096];
#ifdef _WIN32
    if (!proc_) return -1;
    for (;;) {
        DWORD n = 0;
        if (!ReadFile(out_, tmp, sizeof(tmp), &n, nullptr) || n == 0) break;
        buf.append(tmp, n);
        flush(false);
    }
    flush(true);
    WaitForSingleObject(proc_, INFINITE);
    DWORD code = 0;
    GetExitCodeProcess(proc_, &code);
    return int(code);
#else
    if (pid_ <= 0) return -1;
    for (;;) {
        ssize_t n = ::read(fd_, tmp, sizeof(tmp));
        if (n <= 0) { if (n < 0 && errno == EINTR) continue; break; }
        buf.append(tmp, size_t(n));
        flush(false);
    }
    flush(true);
    int st = 0;
    ::waitpid(pid_, &st, 0);
    pid_ = -1;
    if (WIFEXITED(st)) return WEXITSTATUS(st);
    return -1;
#endif
}

void ChildProcess::kill() {
    if (killed_.exchange(true)) return;
#ifdef _WIN32
    if (job_) TerminateJobObject(job_, 1);
    if (proc_) TerminateProcess(proc_, 1);
#else
    if (pid_ > 0) ::kill(-pid_, SIGKILL);
#endif
}

std::string flowVersion(const std::string& flowPath) {
    if (flowPath.empty()) return {};
    ChildProcess p;
    std::string err;
    std::error_code ec;
    if (!p.start(flowPath, {"--version"}, fs::temp_directory_path(ec).u8string(), &err)) return {};
    std::string first;
    const int rc = p.pump([&](const std::string& l) { if (first.empty() && !l.empty()) first = l; });
    return rc == 0 ? first : std::string();
}
