// C41 (Petrel parity 6): recorded workflows.
//
// A workflow is an ordered list of steps that replays commands of the ResSim Studio command registry (the same
// registry the command palette uses), so everything reachable from Ctrl+Shift+P can be recorded and replayed.
// Step kinds:
//   cmd  <label>   run the registry command with exactly this label
//   wait           wait until the simulator is idle (use after "Run simulation")
//   note <text>    comment, does nothing
//   for_variants   ... end   repeat the enclosed steps once per saved model variant (case manager); the loop body
//                  starts with the variant loaded. Loops do not nest.
// File format (.rswf, UTF-8 text): first line "RESSIM_WORKFLOW 1", then one step per line "<kind>\t<argument>".
// Lines starting with '#' and blank lines are ignored.
//
// Runner is a small state machine that executes at most one step per tick() so the UI stays responsive and
// commands that open windows or start background runs take effect between steps. Header-only, no UI dependency.
#pragma once
#include <functional>
#include <sstream>
#include <string>
#include <vector>

namespace ressim::wf {

enum class Kind { Command, Wait, Note, ForVariants, End };

struct Step {
    Kind kind = Kind::Command;
    std::string arg;
};

inline const char* kindName(Kind k) {
    switch (k) {
    case Kind::Command: return "cmd";
    case Kind::Wait: return "wait";
    case Kind::Note: return "note";
    case Kind::ForVariants: return "for_variants";
    case Kind::End: return "end";
    }
    return "?";
}

inline std::string serialize(const std::vector<Step>& steps) {
    std::string s = "RESSIM_WORKFLOW 1\n";
    for (const Step& st : steps) {
        s += kindName(st.kind);
        if (!st.arg.empty()) { s += '\t'; for (char c : st.arg) s += (c == '\n' || c == '\r') ? ' ' : c; }
        s += '\n';
    }
    return s;
}

// Structural check: balanced, non-nested loops. Returns "" when valid.
inline std::string validate(const std::vector<Step>& steps) {
    int depth = 0;
    for (size_t i = 0; i < steps.size(); ++i) {
        if (steps[i].kind == Kind::ForVariants) { if (depth) return "step " + std::to_string(i + 1) + ": loops cannot be nested"; depth = 1; }
        if (steps[i].kind == Kind::End) { if (!depth) return "step " + std::to_string(i + 1) + ": 'end' without a loop"; depth = 0; }
        if (steps[i].kind == Kind::Command && steps[i].arg.empty()) return "step " + std::to_string(i + 1) + ": empty command";
    }
    return depth ? std::string("loop is not closed with 'end'") : std::string();
}

inline bool parse(const std::string& text, std::vector<Step>& out, std::string& err) {
    out.clear(); err.clear();
    std::istringstream is(text);
    std::string line;
    int ln = 0;
    bool header = false;
    while (std::getline(is, line)) {
        ++ln;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (ln == 1 && line.size() >= 3 && (unsigned char)line[0] == 0xEF && (unsigned char)line[1] == 0xBB && (unsigned char)line[2] == 0xBF) line.erase(0, 3);
        if (line.empty() || line[0] == '#') continue;
        if (!header) {
            if (line.rfind("RESSIM_WORKFLOW", 0) != 0) { err = "line " + std::to_string(ln) + ": not a ResSim workflow (missing RESSIM_WORKFLOW header)"; return false; }
            header = true;
            continue;
        }
        const size_t tab = line.find('\t');
        const std::string k = line.substr(0, tab), a = tab == std::string::npos ? std::string() : line.substr(tab + 1);
        Step st;
        st.arg = a;
        if (k == "cmd") st.kind = Kind::Command;
        else if (k == "wait") st.kind = Kind::Wait;
        else if (k == "note") st.kind = Kind::Note;
        else if (k == "for_variants") st.kind = Kind::ForVariants;
        else if (k == "end") st.kind = Kind::End;
        else { err = "line " + std::to_string(ln) + ": unknown step '" + k + "'"; return false; }
        out.push_back(st);
    }
    if (!header) { err = "empty file (missing RESSIM_WORKFLOW header)"; return false; }
    err = validate(out);
    if (!err.empty()) { out.clear(); return false; }
    return true;
}

struct Runner {
    std::vector<Step> steps;
    int pc = -1;              // next step, -1 = idle
    int loopStart = -1, variant = 0;
    int ticksSinceCommand = 0;
    int executed = 0;         // number of executed command steps
    std::string error, log;

    bool active() const { return pc >= 0; }
    bool start(const std::vector<Step>& s) {
        error = validate(s);
        if (!error.empty()) { pc = -1; return false; }
        steps = s; pc = steps.empty() ? -1 : 0; loopStart = -1; variant = 0; ticksSinceCommand = 0; executed = 0; log.clear();
        return true;
    }
    void stop(const std::string& why = "") { if (!why.empty()) error = why; pc = -1; }

    // exec(label) -> false when the label is not in the registry (or the command is disabled)
    // busy()      -> true while a simulation (or another background job) runs
    // nVariants() -> number of saved model variants; loadVariant(i) loads variant i as the current model
    void tick(const std::function<bool(const std::string&)>& exec, const std::function<bool()>& busy,
              const std::function<int()>& nVariants, const std::function<bool(int)>& loadVariant) {
        if (pc < 0) return;
        if (pc >= int(steps.size())) { pc = -1; return; }
        ++ticksSinceCommand;
        const Step& s = steps[size_t(pc)];
        switch (s.kind) {
        case Kind::Note: ++pc; break;
        case Kind::Wait:
            // give a command at least two frames to start its background job before we look at busy()
            if (ticksSinceCommand < 3 || busy()) return;
            ++pc;
            break;
        case Kind::Command:
            if (busy()) return;   // never fire a command into a running simulation
            if (!exec(s.arg)) { stop("step " + std::to_string(pc + 1) + ": command not available: " + s.arg); return; }
            log += "ran: " + s.arg + "\n";
            ++executed; ticksSinceCommand = 0; ++pc;
            break;
        case Kind::ForVariants: {
            const int n = nVariants ? nVariants() : 0;
            if (n <= 0) {   // skip the body
                int j = pc + 1;
                while (j < int(steps.size()) && steps[size_t(j)].kind != Kind::End) ++j;
                pc = j + 1;
                log += "loop skipped: no variants\n";
                break;
            }
            loopStart = pc; variant = 0;
            if (!loadVariant(0)) { stop("could not load variant 1"); return; }
            log += "variant 1 of " + std::to_string(n) + "\n";
            ticksSinceCommand = 0; ++pc;
            break;
        }
        case Kind::End: {
            if (loopStart < 0) { ++pc; break; }
            const int n = nVariants ? nVariants() : 0;
            if (busy()) return;
            if (++variant < n) {
                if (!loadVariant(variant)) { stop("could not load variant " + std::to_string(variant + 1)); return; }
                log += "variant " + std::to_string(variant + 1) + " of " + std::to_string(n) + "\n";
                ticksSinceCommand = 0; pc = loopStart + 1;
            } else { loopStart = -1; ++pc; }
            break;
        }
        }
        if (pc >= int(steps.size())) pc = -1;
    }
};

}  // namespace ressim::wf
