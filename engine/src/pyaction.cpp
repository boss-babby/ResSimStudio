// C37 PYACTION: persistent Python child process, line-based JSON protocol, opm_embedded shim (see pyaction.h).
#include "ressim/pyaction.h"
#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <thread>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
#include <csignal>
#include <fcntl.h>
#include <poll.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace ressim {
namespace pyact {

// ============================================================================ JSON
const Json* Json::get(const std::string& k) const {
    for (const auto& kv : o) if (kv.first == k) return &kv.second;
    return nullptr;
}
static void dumpStr(const std::string& s, std::string& out) {
    out += '"';
    for (unsigned char c : s) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (c < 0x20) { char b[8]; std::snprintf(b, sizeof b, "\\u%04x", c); out += b; }
                else out += char(c);
        }
    }
    out += '"';
}
static void dumpTo(const Json& j, std::string& out) {
    switch (j.t) {
        case Json::Null: out += "null"; break;
        case Json::Bool: out += j.b ? "true" : "false"; break;
        case Json::Num: {
            if (!std::isfinite(j.n)) { out += "null"; break; }
            char b[40]; std::snprintf(b, sizeof b, "%.17g", j.n); out += b; break;
        }
        case Json::Str: dumpStr(j.s, out); break;
        case Json::Arr:
            out += '[';
            for (size_t i = 0; i < j.a.size(); ++i) { if (i) out += ','; dumpTo(j.a[i], out); }
            out += ']';
            break;
        case Json::Obj:
            out += '{';
            for (size_t i = 0; i < j.o.size(); ++i) { if (i) out += ','; dumpStr(j.o[i].first, out); out += ':'; dumpTo(j.o[i].second, out); }
            out += '}';
            break;
    }
}
std::string dump(const Json& j) { std::string s; dumpTo(j, s); return s; }

namespace {
struct Parser {
    const std::string& t;
    size_t i = 0;
    std::string err;
    explicit Parser(const std::string& s) : t(s) {}
    void ws() { while (i < t.size() && std::isspace(static_cast<unsigned char>(t[i]))) ++i; }
    bool fail(const std::string& m) { if (err.empty()) err = m + " at " + std::to_string(i); return false; }
    static void utf8(unsigned cp, std::string& o) {
        if (cp < 0x80) o += char(cp);
        else if (cp < 0x800) { o += char(0xC0 | (cp >> 6)); o += char(0x80 | (cp & 0x3F)); }
        else if (cp < 0x10000) { o += char(0xE0 | (cp >> 12)); o += char(0x80 | ((cp >> 6) & 0x3F)); o += char(0x80 | (cp & 0x3F)); }
        else { o += char(0xF0 | (cp >> 18)); o += char(0x80 | ((cp >> 12) & 0x3F)); o += char(0x80 | ((cp >> 6) & 0x3F)); o += char(0x80 | (cp & 0x3F)); }
    }
    bool str(std::string& o) {
        if (t[i] != '"') return fail("string expected");
        ++i;
        while (i < t.size() && t[i] != '"') {
            char c = t[i++];
            if (c != '\\') { o += c; continue; }
            if (i >= t.size()) return fail("bad escape");
            const char e = t[i++];
            switch (e) {
                case 'n': o += '\n'; break; case 't': o += '\t'; break; case 'r': o += '\r'; break;
                case 'b': o += '\b'; break; case 'f': o += '\f'; break;
                case 'u': {
                    if (i + 4 > t.size()) return fail("bad \\u");
                    unsigned cp = unsigned(std::strtoul(t.substr(i, 4).c_str(), nullptr, 16)); i += 4;
                    if (cp >= 0xD800 && cp < 0xDC00 && i + 6 <= t.size() && t[i] == '\\' && t[i + 1] == 'u') {
                        const unsigned lo = unsigned(std::strtoul(t.substr(i + 2, 4).c_str(), nullptr, 16));
                        if (lo >= 0xDC00 && lo < 0xE000) { cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00); i += 6; }
                    }
                    utf8(cp, o);
                    break;
                }
                default: o += e;
            }
        }
        if (i >= t.size()) return fail("unterminated string");
        ++i;
        return true;
    }
    bool value(Json& v) {
        ws();
        if (i >= t.size()) return fail("value expected");
        const char c = t[i];
        if (c == '{') {
            v.t = Json::Obj; ++i; ws();
            if (i < t.size() && t[i] == '}') { ++i; return true; }
            for (;;) {
                ws();
                std::string k;
                if (i >= t.size() || !str(k)) return fail("key expected");
                ws();
                if (i >= t.size() || t[i] != ':') return fail("':' expected");
                ++i;
                Json x;
                if (!value(x)) return false;
                v.o.emplace_back(k, std::move(x));
                ws();
                if (i < t.size() && t[i] == ',') { ++i; continue; }
                if (i < t.size() && t[i] == '}') { ++i; return true; }
                return fail("',' or '}' expected");
            }
        }
        if (c == '[') {
            v.t = Json::Arr; ++i; ws();
            if (i < t.size() && t[i] == ']') { ++i; return true; }
            for (;;) {
                Json x;
                if (!value(x)) return false;
                v.a.push_back(std::move(x));
                ws();
                if (i < t.size() && t[i] == ',') { ++i; continue; }
                if (i < t.size() && t[i] == ']') { ++i; return true; }
                return fail("',' or ']' expected");
            }
        }
        if (c == '"') { v.t = Json::Str; return str(v.s); }
        if (t.compare(i, 4, "true") == 0) { v.t = Json::Bool; v.b = true; i += 4; return true; }
        if (t.compare(i, 5, "false") == 0) { v.t = Json::Bool; v.b = false; i += 5; return true; }
        if (t.compare(i, 4, "null") == 0) { v.t = Json::Null; i += 4; return true; }
        if (t.compare(i, 3, "NaN") == 0) { v.t = Json::Null; i += 3; return true; }
        if (t.compare(i, 8, "Infinity") == 0 || t.compare(i, 9, "-Infinity") == 0) { v.t = Json::Null; i += t[i] == '-' ? 9 : 8; return true; }
        char* e = nullptr;
        const double d = std::strtod(t.c_str() + i, &e);
        if (e == t.c_str() + i) return fail("bad value");
        v.t = Json::Num; v.n = d; i = size_t(e - t.c_str());
        return true;
    }
};
}  // namespace

bool parse(const std::string& text, Json& out, std::string& err) {
    Parser p(text);
    out = Json();
    if (!p.value(out)) { err = p.err; return false; }
    p.ws();
    if (p.i != text.size()) { err = "trailing characters at " + std::to_string(p.i); return false; }
    return true;
}

// ============================================================================ PYACTION records
static std::vector<std::string> tokens(const std::string& rec) {
    std::vector<std::string> out;
    size_t i = 0;
    while (i < rec.size()) {
        if (std::isspace(static_cast<unsigned char>(rec[i]))) { ++i; continue; }
        if (rec[i] == '\'' || rec[i] == '"') {
            const char q = rec[i];
            const size_t e = rec.find(q, i + 1);
            out.push_back(rec.substr(i + 1, (e == std::string::npos ? rec.size() : e) - i - 1));
            i = e == std::string::npos ? rec.size() : e + 1;
        } else {
            size_t e = i;
            while (e < rec.size() && !std::isspace(static_cast<unsigned char>(rec[e]))) ++e;
            out.push_back(rec.substr(i, e - i));
            i = e;
        }
    }
    return out;
}
std::vector<Action> parseActions(const std::vector<std::pair<int, std::string>>& bodies, std::string& err) {
    std::vector<Action> out;
    for (const auto& [step, body] : bodies) {
        // strip comments, split at '/' outside quotes
        std::string clean;
        std::istringstream is(body);
        std::string line;
        while (std::getline(is, line)) {
            bool q = false; char qc = 0;
            for (size_t i = 0; i < line.size(); ++i) {
                const char c = line[i];
                if (q) { if (c == qc) q = false; continue; }
                if (c == '\'' || c == '"') { q = true; qc = c; continue; }
                if (c == '-' && i + 1 < line.size() && line[i + 1] == '-') { line = line.substr(0, i); break; }
            }
            clean += line + "\n";
        }
        std::vector<std::string> recs;
        std::string cur;
        bool q = false; char qc = 0;
        for (char c : clean) {
            if (q) { cur += c; if (c == qc) q = false; continue; }
            if (c == '\'' || c == '"') { q = true; qc = c; cur += c; continue; }
            if (c == '/') { recs.push_back(cur); cur.clear(); continue; }
            cur += c;
        }
        if (recs.size() < 2) { err += "PYACTION: two records (name / run count, script file) expected; "; continue; }
        const auto r1 = tokens(recs[0]), r2 = tokens(recs[1]);
        if (r1.empty() || r2.empty()) { err += "PYACTION: empty record; "; continue; }
        Action a;
        a.name = r1[0];
        std::string rc = r1.size() > 1 ? r1[1] : "SINGLE";
        for (char& c : rc) c = char(std::toupper(static_cast<unsigned char>(c)));
        if (rc == "1*") rc = "SINGLE";
        if (rc != "SINGLE" && rc != "UNLIMITED" && rc != "FIRST_TRUE") { err += "PYACTION " + a.name + ": unknown run count " + rc + " - SINGLE used; "; rc = "SINGLE"; }
        a.runCount = rc;
        a.file = r2[0];
        a.fromStep = step;
        // a later PYACTION with the same name replaces the earlier one
        bool rep = false;
        for (auto& b : out) if (b.name == a.name) { b = a; rep = true; }
        if (!rep) out.push_back(a);
    }
    return out;
}

// ============================================================================ child process
struct Host::Impl {
#ifdef _WIN32
    HANDLE proc = nullptr, toChild = nullptr, fromChild = nullptr;
#else
    pid_t pid = -1;
    int toChild = -1, fromChild = -1;
#endif
    std::string buf;   // bytes read past the last line
    std::string script;
};

Host::Host() : d_(new Impl) {}
Host::~Host() { stop(); }

bool Host::running() const {
#ifdef _WIN32
    if (!d_->proc) return false;
    DWORD code = 0;
    return GetExitCodeProcess(d_->proc, &code) && code == STILL_ACTIVE;
#else
    if (d_->pid <= 0) return false;
    int st = 0;
    return waitpid(d_->pid, &st, WNOHANG) == 0;
#endif
}

bool Host::writeLine(const std::string& s0) {
    const std::string s = s0 + "\n";
#ifdef _WIN32
    size_t off = 0;
    while (off < s.size()) {
        DWORD w = 0;
        if (!WriteFile(d_->toChild, s.data() + off, DWORD(s.size() - off), &w, nullptr)) return false;
        off += w;
    }
    return true;
#else
    size_t off = 0;
    while (off < s.size()) {
        const ssize_t w = ::write(d_->toChild, s.data() + off, s.size() - off);
        if (w < 0) { if (errno == EINTR) continue; return false; }
        off += size_t(w);
    }
    return true;
#endif
}

bool Host::readLine(std::string& s, double timeoutSec) {
    const auto t0 = std::chrono::steady_clock::now();
    for (;;) {
        const size_t nl = d_->buf.find('\n');
        if (nl != std::string::npos) { s = d_->buf.substr(0, nl); d_->buf.erase(0, nl + 1); if (!s.empty() && s.back() == '\r') s.pop_back(); return true; }
        const double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (el > timeoutSec) return false;
        char tmp[65536];
#ifdef _WIN32
        DWORD avail = 0;
        if (!PeekNamedPipe(d_->fromChild, nullptr, 0, nullptr, &avail, nullptr)) return false;
        if (avail == 0) {
            if (!running()) return false;
            std::this_thread::sleep_for(std::chrono::milliseconds(2));
            continue;
        }
        DWORD r = 0;
        if (!ReadFile(d_->fromChild, tmp, DWORD(std::min<size_t>(sizeof tmp, avail)), &r, nullptr) || r == 0) return false;
        d_->buf.append(tmp, r);
#else
        pollfd pf{d_->fromChild, POLLIN, 0};
        const int ms = int(std::min(1000.0, std::max(1.0, (timeoutSec - el) * 1000.0)));
        const int pr = ::poll(&pf, 1, ms);
        if (pr < 0) { if (errno == EINTR) continue; return false; }
        if (pr == 0) continue;
        const ssize_t r = ::read(d_->fromChild, tmp, sizeof tmp);
        if (r < 0) { if (errno == EINTR) continue; return false; }
        if (r == 0) return false;   // EOF: child gone
        d_->buf.append(tmp, size_t(r));
#endif
    }
}

bool Host::start(const std::string& python, std::string& err) {
    if (running()) return true;
    stop();
    // bootstrap script in the temp directory (one per process id: parallel runs do not clash)
    std::error_code ec;
    std::filesystem::path dir = std::filesystem::temp_directory_path(ec);
    if (ec) dir = ".";
#ifdef _WIN32
    const unsigned long me = GetCurrentProcessId();
#else
    const unsigned long me = (unsigned long)getpid();
#endif
    const std::filesystem::path sp = dir / ("ressim_pyaction_host_" + std::to_string(me) + ".py");
    {
        std::ofstream f(sp, std::ios::binary);
        if (!f) { err = "cannot write the PYACTION host script to " + sp.string(); return false; }
        f << hostScript();
    }
    d_->script = sp.string();
    std::vector<std::string> cands;
    if (!python.empty()) cands.push_back(python);
    if (const char* e = std::getenv("RESSIM_PYTHON")) if (*e) cands.push_back(e);
#ifdef _WIN32
    cands.push_back("python"); cands.push_back("python3"); cands.push_back("py");
#else
    cands.push_back("python3"); cands.push_back("python");
#endif
    std::string tried;
    for (const auto& exe : cands) {
        tried += (tried.empty() ? "" : ", ") + exe;
#ifdef _WIN32
        SECURITY_ATTRIBUTES sa{sizeof(SECURITY_ATTRIBUTES), nullptr, TRUE};
        HANDLE inR = nullptr, inW = nullptr, outR = nullptr, outW = nullptr;
        if (!CreatePipe(&inR, &inW, &sa, 0)) continue;
        if (!CreatePipe(&outR, &outW, &sa, 0)) { CloseHandle(inR); CloseHandle(inW); continue; }
        SetHandleInformation(inW, HANDLE_FLAG_INHERIT, 0);
        SetHandleInformation(outR, HANDLE_FLAG_INHERIT, 0);
        STARTUPINFOA si; ZeroMemory(&si, sizeof si); si.cb = sizeof si;
        si.dwFlags = STARTF_USESTDHANDLES; si.hStdInput = inR; si.hStdOutput = outW; si.hStdError = GetStdHandle(STD_ERROR_HANDLE);
        PROCESS_INFORMATION pi; ZeroMemory(&pi, sizeof pi);
        std::string cmd = "\"" + exe + "\" " + (exe == "py" ? "-3 " : "") + "-u \"" + d_->script + "\"";
        std::vector<char> cl(cmd.begin(), cmd.end()); cl.push_back(0);
        const BOOL ok = CreateProcessA(nullptr, cl.data(), nullptr, nullptr, TRUE, CREATE_NO_WINDOW, nullptr, nullptr, &si, &pi);
        CloseHandle(inR); CloseHandle(outW);
        if (!ok) { CloseHandle(inW); CloseHandle(outR); continue; }
        CloseHandle(pi.hThread);
        d_->proc = pi.hProcess; d_->toChild = inW; d_->fromChild = outR;
#else
        int in[2], out[2];
        if (::pipe(in) != 0) continue;
        if (::pipe(out) != 0) { ::close(in[0]); ::close(in[1]); continue; }
        static bool sigInit = false;
        if (!sigInit) { sigInit = true; struct sigaction old; if (sigaction(SIGPIPE, nullptr, &old) == 0 && old.sa_handler == SIG_DFL) std::signal(SIGPIPE, SIG_IGN); }
        const pid_t pid = ::fork();
        if (pid < 0) { ::close(in[0]); ::close(in[1]); ::close(out[0]); ::close(out[1]); continue; }
        if (pid == 0) {
            ::dup2(in[0], 0); ::dup2(out[1], 1);
            ::close(in[0]); ::close(in[1]); ::close(out[0]); ::close(out[1]);
            ::execlp(exe.c_str(), exe.c_str(), "-u", d_->script.c_str(), (char*)nullptr);
            ::_exit(127);
        }
        ::close(in[0]); ::close(out[1]);
        ::fcntl(in[1], F_SETFD, FD_CLOEXEC); ::fcntl(out[0], F_SETFD, FD_CLOEXEC);
        d_->pid = pid; d_->toChild = in[1]; d_->fromChild = out[0];
#endif
        d_->buf.clear();
        // handshake
        std::string line;
        if (writeLine("{\"cmd\":\"ping\"}") && readLine(line, 30)) {
            Json r; std::string pe;
            if (parse(line, r, pe) && r.get("pong")) {
                if (const Json* v = r.get("version")) version_ = v->s;
                return true;
            }
        }
        stop();
    }
    err = "no working Python 3 interpreter found (tried " + tried + ") - install Python 3 or set RESSIM_PYTHON / the Python path in the run options";
    return false;
}

bool Host::call(const Json& request, const Query& query, Json& response, std::string& err, double timeoutSec) {
    if (!running()) { err = "the Python process is not running"; return false; }
    if (!writeLine(dump(request))) { err = "cannot write to the Python process"; return false; }
    for (int guard = 0; guard < 10000000; ++guard) {
        std::string line;
        if (!readLine(line, timeoutSec)) { err = running() ? "the Python script did not answer within " + std::to_string(int(timeoutSec)) + " s" : "the Python process exited"; stop(); return false; }
        Json m; std::string pe;
        if (!parse(line, m, pe)) { err = "bad message from Python: " + pe + " - " + line.substr(0, 200); return false; }
        if (const Json* q = m.get("q")) {
            Json ans = Json::obj();
            ans.set("v", query ? query(*q) : Json());
            if (!writeLine(dump(ans))) { err = "cannot write to the Python process"; return false; }
            continue;
        }
        response = std::move(m);
        return true;
    }
    err = "too many queries";
    return false;
}

void Host::stop() {
#ifdef _WIN32
    if (d_->toChild) { writeLine("{\"cmd\":\"quit\"}"); CloseHandle(d_->toChild); d_->toChild = nullptr; }
    if (d_->fromChild) { CloseHandle(d_->fromChild); d_->fromChild = nullptr; }
    if (d_->proc) {
        if (WaitForSingleObject(d_->proc, 2000) != WAIT_OBJECT_0) TerminateProcess(d_->proc, 1);
        CloseHandle(d_->proc); d_->proc = nullptr;
    }
#else
    if (d_->toChild >= 0) { writeLine("{\"cmd\":\"quit\"}"); ::close(d_->toChild); d_->toChild = -1; }
    if (d_->fromChild >= 0) { ::close(d_->fromChild); d_->fromChild = -1; }
    if (d_->pid > 0) {
        int st = 0;
        bool gone = false;
        for (int i = 0; i < 200 && !gone; ++i) { if (waitpid(d_->pid, &st, WNOHANG) == d_->pid) gone = true; else std::this_thread::sleep_for(std::chrono::milliseconds(10)); }
        if (!gone) { ::kill(d_->pid, SIGKILL); waitpid(d_->pid, &st, 0); }
        d_->pid = -1;
    }
#endif
    if (!d_->script.empty()) { std::error_code ec; std::filesystem::remove(d_->script, ec); d_->script.clear(); }
    d_->buf.clear();
}

// ============================================================================ the child's bootstrap
const char* hostScript() {
    return R"PY(# ResSim PYACTION host (generated) - line-based JSON over stdin / stdout
import sys, json, types, datetime, traceback, io, fnmatch, math
_proto_in = sys.stdin
_proto_out = sys.stdout
def _send(obj):
    _proto_out.write(json.dumps(obj, allow_nan=False, default=str) + "\n"); _proto_out.flush()
def _recv():
    line = _proto_in.readline()
    if not line:
        sys.exit(0)
    return json.loads(line)
_cache = {}
def _query(kind, name, ent=""):
    key = (kind, name, ent)
    if key in _cache:
        return _cache[key]
    _send({"q": [kind, name, ent]})
    v = _recv().get("v")
    _cache[key] = v
    return v

mod = types.ModuleType("opm_embedded")
class OpmLog:
    _buf = []
    @staticmethod
    def info(m): OpmLog._buf.append(["info", str(m)])
    @staticmethod
    def note(m): OpmLog._buf.append(["note", str(m)])
    @staticmethod
    def warning(m): OpmLog._buf.append(["warning", str(m)])
    @staticmethod
    def error(m): OpmLog._buf.append(["error", str(m)])
    @staticmethod
    def problem(m): OpmLog._buf.append(["problem", str(m)])
    @staticmethod
    def bug(m): OpmLog._buf.append(["bug", str(m)])
    @staticmethod
    def debug(m): OpmLog._buf.append(["debug", str(m)])

_updates = {}
class SummaryState:
    def __init__(self, req):
        self._req = req
    def elapsed(self):
        return float(self._req.get("elapsed", 0.0))
    def __getitem__(self, key):
        key = str(key)
        if key in _updates:
            return _updates[key]
        v = _query("F", key)
        if v is None:
            raise KeyError(key)
        return v
    def __contains__(self, key):
        key = str(key)
        return key in _updates or _query("F", key) is not None
    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default
    def has(self, key):
        return key in self
    def update(self, key, value):
        _updates[str(key)] = float(value)
    def well_var(self, well, var):
        v = _query("W", str(var), str(well))
        if v is None:
            raise KeyError(str(well) + ":" + str(var))
        return v
    def group_var(self, group, var):
        v = _query("G", str(var), str(group))
        if v is None:
            raise KeyError(str(group) + ":" + str(var))
        return v
    def has_well_var(self, well, var=None):
        if var is None:
            return str(well) in self._req.get("wells", [])
        return _query("W", str(var), str(well)) is not None
    def has_group_var(self, group, var=None):
        if var is None:
            return str(group) in self._req.get("groups", [])
        return _query("G", str(var), str(group)) is not None
    def wells(self, pattern="*"):
        return [w for w in self._req.get("wells", []) if fnmatch.fnmatchcase(w.upper(), str(pattern).upper())]
    def groups(self, pattern="*"):
        return [g for g in self._req.get("groups", []) if fnmatch.fnmatchcase(g.upper(), str(pattern).upper())]
    def update_well_var(self, well, var, value):
        _updates[str(var) + ":" + str(well)] = float(value)
    def update_group_var(self, group, var, value):
        _updates[str(var) + ":" + str(group)] = float(value)

class Well:
    def __init__(self, name, status):
        self.name = name
        self._status = status
    def status(self):
        return self._status

class Schedule:
    def __init__(self, req):
        self._req = req
        self._kw = []
        y, m, d = req.get("start", [2000, 1, 1])
        self.start = datetime.datetime(int(y), int(m), int(d))
        self.end = None
    def _welopen(self, name, what):
        self._kw.append("WELOPEN\n '%s' %s /\n/\n" % (name, what))
    def open_well(self, name, step=None): self._welopen(name, "OPEN")
    def shut_well(self, name, step=None): self._welopen(name, "SHUT")
    def stop_well(self, name, step=None): self._welopen(name, "STOP")
    def insert_keywords(self, keywords, step=None, unit_system=None):
        if isinstance(keywords, (list, tuple)):
            keywords = "\n".join(str(k) for k in keywords)
        self._kw.append(str(keywords))
    def well_names(self, pattern="*"):
        return [w for w in self._req.get("wells", []) if fnmatch.fnmatchcase(w.upper(), str(pattern).upper())]
    def group_names(self, pattern="*"):
        return [g for g in self._req.get("groups", []) if fnmatch.fnmatchcase(g.upper(), str(pattern).upper())]
    def get_well(self, name, step=None):
        st = self._req.get("well_status", {}).get(str(name))
        if st is None:
            raise KeyError(str(name))
        return Well(str(name), st)
    def get_wells(self, step=None):
        return [Well(w, self._req.get("well_status", {}).get(w, "OPEN")) for w in self._req.get("wells", [])]
    def __contains__(self, name):
        return str(name) in self._req.get("wells", [])

class EclState:
    def __init__(self, req):
        self.title = req.get("title", "")
        self.units = req.get("units", "")

mod.OpmLog = OpmLog
mod.SummaryState = SummaryState
mod.Schedule = Schedule
mod.current_summary_state = None
mod.current_schedule = None
mod.current_report_step = 0
mod.current_ecl_state = None
sys.modules["opm_embedded"] = mod

_ns = {}
_code = {}
_legacy = {}
while True:
    req = _recv()
    cmd = req.get("cmd")
    if cmd == "quit":
        break
    if cmd == "ping":
        _send({"pong": True, "version": sys.version.split()[0]})
        continue
    if cmd == "check":
        try:
            compile(req.get("source", ""), req.get("file", "<script>"), "exec")
            _send({"ok": True})
        except SyntaxError as e:
            _send({"ok": False, "error": "%s (line %s, column %s)" % (e.msg, e.lineno, e.offset), "line": e.lineno or 0})
        except BaseException as e:
            _send({"ok": False, "error": str(e), "line": 0})
        continue
    if cmd == "reset":
        _ns.clear(); _code.clear(); _updates.clear()
        _send({"ok": True})
        continue
    resp = {"ok": True, "result": None, "keywords": "", "actionx": [], "log": []}
    _cache.clear()
    cap = io.StringIO()
    old_out = sys.stdout
    sys.stdout = cap
    try:
        ss = SummaryState(req); sch = Schedule(req); es = EclState(req)
        mod.current_summary_state = ss; mod.current_schedule = sch
        mod.current_report_step = int(req.get("step", 0)); mod.current_ecl_state = es
        name = req["name"]; path = req["file"]
        first = name not in _ns
        if first:
            with open(path, "r", encoding="utf-8") as f:
                src = f.read()
            _code[name] = compile(src, path, "exec")
            _ns[name] = {"__name__": "pyaction_" + name, "__file__": path, "__builtins__": __builtins__}
            exec(_code[name], _ns[name])
        ns = _ns[name]
        if callable(ns.get("run")):
            def _cb(action_name, wells=None):
                resp["actionx"].append([str(action_name), [str(w) for w in (wells or [])]])
            r = ns["run"](es, sch, mod.current_report_step, ss, _cb)
            resp["result"] = None if r is None else bool(r)
        elif not first:
            exec(_code[name], ns)
        resp["keywords"] = "\n".join(sch._kw)
    except BaseException:
        resp["ok"] = False
        resp["error"] = traceback.format_exc()
    finally:
        sys.stdout = old_out
    out = cap.getvalue()
    for ln in out.splitlines():
        OpmLog._buf.append(["print", ln])
    resp["log"] = OpmLog._buf
    OpmLog._buf = []
    _send(resp)
)PY";
}

}  // namespace pyact
}  // namespace ressim
