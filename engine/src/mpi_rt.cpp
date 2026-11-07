// C16: run-time loaded MPI (MPICH ABI: MS-MPI, MPICH, Intel MPI) and the distributed batch job used by ResSimCLI.
#include "ressim/mpi_rt.h"
#include "ressim/io.h"
#include "ressim/study.h"
#include "ressim/results_file.h"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <thread>
#include <vector>
#include <algorithm>
#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <dlfcn.h>
#include <unistd.h>
#endif

namespace ressim {
namespace mpi {
namespace {

// MPICH-ABI constants (identical in MS-MPI, which is derived from MPICH)
constexpr int kCommWorld = 0x44000000;
constexpr int kChar = 0x4c000101, kInt = 0x4c000405, kDouble = 0x4c00080b;
constexpr int kSum = 0x58000003, kMax = 0x58000001;
constexpr int kAnySource = -2, kAnyTag = -1;
struct Status { int countLo, countHiAndCancelled, source, tag, error; };

using FnInit = int (*)(int*, char***);
using FnInitialized = int (*)(int*);
using FnVoid = int (*)();
using FnCommInt = int (*)(int, int*);
using FnSend = int (*)(const void*, int, int, int, int, int);
using FnRecv = int (*)(void*, int, int, int, int, int, Status*);
using FnProbe = int (*)(int, int, int, Status*);
using FnIprobe = int (*)(int, int, int, int*, Status*);
using FnGetCount = int (*)(const Status*, int, int*);
using FnAllreduce = int (*)(const void*, void*, int, int, int, int);
using FnBarrier = int (*)(int);
using FnAbort = int (*)(int, int);
using FnProcName = int (*)(char*, int*);
using FnBcast = int (*)(void*, int, int, int, int);
using FnSendrecv = int (*)(const void*, int, int, int, int, void*, int, int, int, int, int, Status*);

struct Api {
    void* lib = nullptr;
    std::string name;
    FnInit init = nullptr; FnInitialized initialized = nullptr; FnVoid finalize = nullptr;
    FnCommInt rank = nullptr, size = nullptr;
    FnSend send = nullptr; FnRecv recv = nullptr; FnProbe probe = nullptr; FnIprobe iprobe = nullptr; FnGetCount getCount = nullptr;
    FnAllreduce allreduce = nullptr; FnBarrier barrier = nullptr; FnAbort abort = nullptr; FnProcName procName = nullptr;
    FnBcast bcast = nullptr; FnSendrecv sendrecv = nullptr;
    bool inited = false;
} g;

void* openLib(const std::string& n) {
#ifdef _WIN32
    return (void*)LoadLibraryA(n.c_str());
#else
    return dlopen(n.c_str(), RTLD_NOW | RTLD_GLOBAL);
#endif
}
void* sym(void* lib, const char* s) {
#ifdef _WIN32
    return (void*)GetProcAddress((HMODULE)lib, s);
#else
    return dlsym(lib, s);
#endif
}
template <class F> bool bind(F& f, const char* s) { f = reinterpret_cast<F>(sym(g.lib, s)); return f != nullptr; }

} // namespace

bool load(std::string* error) {
    if (g.lib) return true;
    std::vector<std::string> cand;
    if (const char* e = std::getenv("RESSIM_MPI_LIB"); e && *e) cand.push_back(e);
#ifdef _WIN32
    cand.push_back("msmpi.dll");
    if (const char* e = std::getenv("MSMPI_BIN"); e && *e) cand.push_back(std::string(e) + "\\..\\..\\..\\Windows\\System32\\msmpi.dll");
    cand.push_back("C:\\Windows\\System32\\msmpi.dll");
    cand.push_back("impi.dll");
#else
    cand.push_back("libmpi.so.12"); cand.push_back("libmpich.so.12"); cand.push_back("libmpi.so"); cand.push_back("libmpich.so");
#endif
    std::string tried;
    for (const auto& c : cand) {
        void* h = openLib(c);
        if (!h) { tried += (tried.empty() ? "" : ", ") + c; continue; }
        g.lib = h; g.name = c;
        if (sym(h, "ompi_mpi_comm_world")) {
            if (error) *error = c + " is Open MPI (pointer-handle ABI) - install MPICH / Intel MPI / MS-MPI or set RESSIM_MPI_LIB";
            g.lib = nullptr; return false;
        }
        bool ok = bind(g.init, "MPI_Init") && bind(g.initialized, "MPI_Initialized") && bind(g.finalize, "MPI_Finalize") &&
                  bind(g.rank, "MPI_Comm_rank") && bind(g.size, "MPI_Comm_size") && bind(g.send, "MPI_Send") &&
                  bind(g.recv, "MPI_Recv") && bind(g.probe, "MPI_Probe") && bind(g.getCount, "MPI_Get_count") &&
                  bind(g.allreduce, "MPI_Allreduce") && bind(g.barrier, "MPI_Barrier") && bind(g.abort, "MPI_Abort") &&
                  bind(g.procName, "MPI_Get_processor_name") && bind(g.bcast, "MPI_Bcast") && bind(g.sendrecv, "MPI_Sendrecv");
        if (!ok) { if (error) *error = c + ": missing MPI entry points"; g.lib = nullptr; return false; }
        bind(g.iprobe, "MPI_Iprobe");                // optional: polite (sleeping) waits, see politeWait()
        return true;
    }
    if (error) *error = "no MPI library found (tried " + tried + "). Windows: install Microsoft MPI (msmpisetup.exe)";
    return false;
}

bool loaded() { return g.lib != nullptr; }
std::string libraryName() { return g.name; }

bool init() {
    if (!g.lib) return false;
    if (g.inited) return true;
    int flag = 0; g.initialized(&flag);
    if (!flag && g.init(nullptr, nullptr) != 0) return false;
    g.inited = true;
    return true;
}
void finalize() { if (g.inited) { g.finalize(); g.inited = false; } }
int rank() { int r = 0; if (g.inited) g.rank(kCommWorld, &r); return r; }
int size() { int s = 1; if (g.inited) g.size(kCommWorld, &s); return s; }
void barrier() { if (g.inited) g.barrier(kCommWorld); }
void abort(int code) { if (g.inited) g.abort(kCommWorld, code); std::exit(code); }
std::string processorName() {
    if (!g.inited) return "";
    char buf[512] = {}; int len = 0; g.procName(buf, &len);
    return std::string(buf, size_t(std::max(0, std::min(len, 511))));
}
bool send(const std::string& d, int dest, int tag) {
    if (!g.inited || d.size() > size_t(0x7fffffff)) return false;
    return g.send(d.data(), int(d.size()), kChar, dest, tag, kCommWorld) == 0;
}
// C21: MPI implementations busy-poll inside blocking calls (MPI_Probe / MPI_Recv / MPI_Bcast), so an idle rank burns a full
// core - on a laptop that steals the turbo / thermal budget of the rank that actually works and shows as 100 % CPU in the
// Task Manager. Idle phases (worker waiting for the next solve, batch dispatcher waiting for results) first wait here:
// a short spin (back-to-back messages stay fast), then MPI_Iprobe + sleep with a 1 ms timer. RESSIM_MPI_SPIN=1 disables it.
bool politeWait(int source, int tag) {
    if (!g.inited || !g.iprobe) return true;
    static const bool spin = [] { const char* e = std::getenv("RESSIM_MPI_SPIN"); return e && *e == '1'; }();
    if (spin) return true;
#ifdef _WIN32
    static const bool fineTimer = [] {                  // Sleep(1) is ~15.6 ms without timeBeginPeriod(1)
        if (HMODULE w = LoadLibraryA("winmm.dll"))
            if (auto f = reinterpret_cast<unsigned (WINAPI*)(unsigned)>(reinterpret_cast<void*>(GetProcAddress(w, "timeBeginPeriod")))) f(1);
        return true;
    }();
    (void)fineTimer;
#endif
    const auto t0 = std::chrono::steady_clock::now();
    for (;;) {
        int flag = 0; Status st{};
        if (g.iprobe(source < 0 ? kAnySource : source, tag < 0 ? kAnyTag : tag, kCommWorld, &flag, &st) != 0) return false;
        if (flag) return true;
        const double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (el < 2e-4) continue;
        std::this_thread::sleep_for(std::chrono::microseconds(el < 0.02 ? 50 : 500));
    }
}
bool recv(std::string& d, int source, int tag, int* srcOut, int* tagOut) {
    if (!g.inited) return false;
    if (!politeWait(source, tag)) return false;
    Status st{};
    if (g.probe(source < 0 ? kAnySource : source, tag < 0 ? kAnyTag : tag, kCommWorld, &st) != 0) return false;
    int n = 0; g.getCount(&st, kChar, &n);
    d.resize(size_t(std::max(0, n)));
    Status st2{};
    if (g.recv(n > 0 ? &d[0] : nullptr, n, kChar, st.source, st.tag, kCommWorld, &st2) != 0) return false;
    if (srcOut) *srcOut = st.source;
    if (tagOut) *tagOut = st.tag;
    return true;
}
bool allreduceSumN(double* v, int n) {
    if (!g.inited || n <= 0) return n <= 0 || !g.inited;
    std::vector<double> o(static_cast<size_t>(n));
    if (g.allreduce(v, o.data(), n, kDouble, kSum, kCommWorld) != 0) return false;
    std::copy(o.begin(), o.end(), v);
    return true;
}
bool bcast(void* data, int bytes, int root) { return !g.inited || g.bcast(data, bytes, kChar, root, kCommWorld) == 0; }
bool sendDoubles(const double* v, int n, int dest, int tag) { return g.inited && g.send(v, n, kDouble, dest, tag, kCommWorld) == 0; }
bool recvDoubles(double* v, int n, int source, int tag) { Status st{}; return g.inited && g.recv(v, n, kDouble, source, tag, kCommWorld, &st) == 0; }
bool sendInts(const int* v, int n, int dest, int tag) { return g.inited && g.send(v, n, kInt, dest, tag, kCommWorld) == 0; }
bool recvInts(int* v, int n, int source, int tag) { Status st{}; return g.inited && g.recv(v, n, kInt, source, tag, kCommWorld, &st) == 0; }
bool sendrecvDoubles(const double* s, int ns, int dest, double* r, int nr, int source, int tag) {
    Status st{};
    return g.inited && g.sendrecv(s, ns, kDouble, dest, tag, r, nr, kDouble, source, tag, kCommWorld, &st) == 0;
}
bool sendrecvInts(const int* s, int ns, int dest, int* r, int nr, int source, int tag) {
    Status st{};
    return g.inited && g.sendrecv(s, ns, kInt, dest, tag, r, nr, kInt, source, tag, kCommWorld, &st) == 0;
}
double allreduceSum(double v) { double o = v; if (g.inited) g.allreduce(&v, &o, 1, kDouble, kSum, kCommWorld); return o; }
double allreduceMax(double v) { double o = v; if (g.inited) g.allreduce(&v, &o, 1, kDouble, kMax, kCommWorld); return o; }

// ============================================================================ distributed batch
namespace {
namespace fs = std::filesystem;
enum Tag { kReady = 1, kJob = 2, kStop = 3, kResult = 4 };

std::string readAll(const fs::path& p) { std::ifstream f(p, std::ios::binary); std::stringstream ss; ss << f.rdbuf(); return ss.str(); }
void writeAll(const fs::path& p, const std::string& d) {
    const fs::path tmp = p.string() + ".part";
    { std::ofstream f(tmp, std::ios::binary); f.write(d.data(), std::streamsize(d.size())); }
    std::error_code ec; fs::rename(tmp, p, ec);
    if (ec) { fs::remove(p, ec); fs::rename(tmp, p, ec); }
}

// run one case text (ResSim project text, or an ECLIPSE deck if it starts like one) -> .rsres bytes
std::string runCaseBytes(const std::string& name, const std::string& text, int threads) {
    ResultsFile rf;
    rf.engine = "ResSim";
    const auto t0 = std::chrono::steady_clock::now();
    try {
        std::string ext = fs::path(name).extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        Case c = (ext == ".data") ? eclipseDeckFromString(text) : projectFromText(text);
        if (threads > 0) c.tuning.threads = threads;
        RunResult r = runSingle(c, nullptr);
        rf.summary = r.summary; rf.stats = r.stats;
        rf.title = r.ok ? "OK" : "ERROR: " + r.error;
        // C28: the coupled-mechanics KPIs (seismicity / uplift study responses) travel back as the blob "KPI"
        if (!r.kpi.empty()) {
            std::string k; char ln[64];
            for (const auto& e : r.kpi) { std::snprintf(ln, sizeof ln, " %.17g\n", e.second); k += e.first + ln; }
            rf.blobs.push_back({"KPI", k});
        }
    } catch (const std::exception& e) { rf.title = std::string("ERROR: ") + e.what(); }
    rf.totalTime = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    rf.created = processorName();
    std::error_code ec;
    fs::path tmp = fs::temp_directory_path(ec);
    if (ec) tmp = ".";
#ifdef _WIN32
    const long pid = long(GetCurrentProcessId());
#else
    const long pid = long(getpid());
#endif
    tmp /= "ressim_mpi_" + std::to_string(pid) + "_" + std::to_string(rank()) + ".rsres";
    try { saveResultsFile(rf, tmp.string()); } catch (...) { return std::string(); }
    std::string bytes = readAll(tmp);
    fs::remove(tmp, ec);
    return bytes;
}
} // namespace

int runBatchJob(const std::string& jobdir, bool verbose) {
    const int me = rank(), np = size();
    int threads = 1;                                   // one rank per core is the usual layout; job.txt may override
    if (me == 0) {
        std::vector<fs::path> cases;
        std::error_code ec;
        for (const auto& e : fs::directory_iterator(jobdir, ec)) {
            const std::string fn = e.path().filename().string();
            if (fn.rfind("case_", 0) == 0) {
                std::string ext = e.path().extension().string();
                std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
                if (ext == ".data" || ext == ".rsim") cases.push_back(e.path());
            }
        }
        std::sort(cases.begin(), cases.end());
        if (fs::exists(fs::path(jobdir) / "job.txt")) {
            std::istringstream is(readAll(fs::path(jobdir) / "job.txt")); std::string k;
            while (is >> k) if (k == "threads") is >> threads;
        }
        const int n = int(cases.size());
        const fs::path prog = fs::path(jobdir) / "progress.txt";
        int done = 0;
        auto writeProgress = [&]() { writeAll(prog, std::to_string(done) + " " + std::to_string(n) + " " + std::to_string(np) + "\n"); };
        writeProgress();
        if (verbose) std::printf("[mpi] %d cases on %d ranks (%s)\n", n, np, libraryName().c_str());
        auto stem = [&](int i) { return (fs::path(jobdir) / cases[size_t(i)].stem()).string() + ".rsres"; };
        auto cancelled = [&]() { return fs::exists(fs::path(jobdir) / "cancel.txt"); };
        if (np == 1) {
            for (int i = 0; i < n && !cancelled(); ++i) {
                writeAll(stem(i), runCaseBytes(cases[size_t(i)].filename().string(), readAll(cases[size_t(i)]), 0));
                ++done; writeProgress();
            }
        } else {
            int next = 0, active = np - 1;
            while (active > 0) {
                std::string msg; int src = 0, tag = 0;
                if (!recv(msg, -1, -1, &src, &tag)) { std::fprintf(stderr, "[mpi] receive failed\n"); abort(3); }
                if (tag == kResult) {
                    const size_t nl = msg.find('\n');
                    const int idx = std::atoi(msg.substr(0, nl).c_str());
                    if (idx >= 0 && idx < n) writeAll(stem(idx), msg.substr(nl + 1));
                    ++done; writeProgress();
                    if (verbose) { std::printf("[mpi] case %d done on rank %d (%d/%d)\n", idx, src, done, n); std::fflush(stdout); }
                }
                if (next < n && !cancelled()) {
                    const std::string job = std::to_string(next) + "\n" + cases[size_t(next)].filename().string() + "\n" + std::to_string(threads) + "\n" +
                                            readAll(cases[size_t(next)]);
                    send(job, src, kJob); ++next;
                } else { send(std::string(), src, kStop); --active; }
            }
        }
        writeAll(fs::path(jobdir) / "done.txt", std::to_string(done) + " " + std::to_string(n) + "\n");
        return done == n ? 0 : 1;
    }
    // worker ranks
    send(std::string(), 0, kReady);
    for (;;) {
        std::string msg; int tag = 0;
        if (!recv(msg, 0, -1, nullptr, &tag)) return 3;
        if (tag != kJob) break;
        const size_t a = msg.find('\n'), b = msg.find('\n', a + 1), c = msg.find('\n', b + 1);
        if (a == std::string::npos || b == std::string::npos || c == std::string::npos) break;
        const std::string idx = msg.substr(0, a), name = msg.substr(a + 1, b - a - 1);
        const int thr = std::atoi(msg.substr(b + 1, c - b - 1).c_str());
        std::string res = idx + "\n" + runCaseBytes(name, msg.substr(c + 1), thr);
        send(res, 0, kResult);
    }
    return 0;
}

} // namespace mpi
} // namespace ressim
