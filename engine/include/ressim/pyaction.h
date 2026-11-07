// C37: PYACTION (OPM Flow keyword) in the ResSim black-oil engine.
//
//   PYACTION
//    'NAME' 'SINGLE' | 'UNLIMITED' | 'FIRST_TRUE' /
//    'script.py' /                        -- relative to the deck directory
//
// The scripts run in ONE persistent Python child process (python3 / python, or RESSIM_PYTHON, or Case::pythonExe) that
// talks to the engine over pipes with line-based JSON (fork + exec on POSIX, CreateProcess on Windows). The child provides
// an 'opm_embedded' module like OPM Flow 2024.04+:
//   opm_embedded.current_summary_state   SummaryState: ss['FOPR'], ss.well_var('P1', 'WOPR'), ss.group_var('G', 'GOPR'),
//                                        ss.elapsed() (seconds), ss.wells(), ss.groups(), has_well_var / has_group_var,
//                                        ss.update(key, value) (script-defined field values, kept between calls)
//   opm_embedded.current_schedule        Schedule: open_well / shut_well / stop_well(name), insert_keywords(text),
//                                        well_names(pattern), get_well(name).status(), .start (datetime)
//   opm_embedded.current_report_step     int
//   opm_embedded.current_ecl_state       minimal (title, units)
//   opm_embedded.OpmLog                  info / warning / error / debug / note -> the run log
// Any summary vector the engine knows (the UDQ / ACTIONX source: F*, W*, G* vectors and UDQ values) is fetched on demand
// from the engine. Module globals persist between calls ('if "setup_done" not in locals():' works); the legacy form
// 'def run(ecl_state, schedule, report_step, summary_state, actionx_callback)' is supported (module body executed once,
// run() called every time, actionx_callback(name, wells) applies the body of ACTIONX 'name').
// print() output goes to the run log. Keywords returned by a call are applied at the report-step boundary through the
// ACTIONX machinery (same as an ACTIONX body; replayed on restart). RUN_COUNT: SINGLE = once, FIRST_TRUE = until the
// script (legacy run) returns True, UNLIMITED = every call.
// Differences to OPM Flow: the scripts are called at the end of every REPORT step (OPM: every time step), the Python
// state is not part of the restart file, ecl_state is minimal.
#pragma once
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace ressim {
namespace pyact {

// minimal JSON value (objects keep insertion order)
struct Json {
    enum Type { Null, Bool, Num, Str, Arr, Obj } t = Null;
    bool b = false;
    double n = 0;
    std::string s;
    std::vector<Json> a;
    std::vector<std::pair<std::string, Json>> o;
    const Json* get(const std::string& k) const;
    static Json str(const std::string& v) { Json j; j.t = Str; j.s = v; return j; }
    static Json num(double v) { Json j; j.t = Num; j.n = v; return j; }
    static Json boolean(bool v) { Json j; j.t = Bool; j.b = v; return j; }
    static Json arr() { Json j; j.t = Arr; return j; }
    static Json obj() { Json j; j.t = Obj; return j; }
    Json& set(const std::string& k, Json v) { o.emplace_back(k, std::move(v)); return *this; }
};
std::string dump(const Json& j);
bool parse(const std::string& text, Json& out, std::string& err);

struct Action {
    std::string name, file, runCount = "SINGLE";   // file: absolute or relative to the deck directory
    int fromStep = 0;
    int calls = 0;
    bool done = false;
};
// PYACTION records in the extra keywords (body text of each keyword at its step)
std::vector<Action> parseActions(const std::vector<std::pair<int, std::string>>& bodies, std::string& err);

// the Python child process
class Host {
public:
    Host();
    ~Host();
    Host(const Host&) = delete;
    Host& operator=(const Host&) = delete;
    // starts the child (python executable: 'python' argument, else RESSIM_PYTHON, else python3 / python); false + message
    bool start(const std::string& python, std::string& err);
    bool running() const;
    // one script call: request -> (queries answered by 'query') -> final response
    using Query = std::function<Json(const Json& q)>;
    bool call(const Json& request, const Query& query, Json& response, std::string& err, double timeoutSec = 600);
    void stop();
    std::string version() const { return version_; }
private:
    struct Impl;
    std::unique_ptr<Impl> d_;
    std::string version_;
    bool writeLine(const std::string& s);
    bool readLine(std::string& s, double timeoutSec);
};

// the bootstrap script run by the child (written to a temp file)
const char* hostScript();

}  // namespace pyact
}  // namespace ressim
