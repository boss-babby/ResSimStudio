// C31: user-defined quantities (ECLIPSE / OPM UDQ, UDA) and actions (ACTIONX) evaluated by the ResSim black-oil engine.
//
// UDQ  (SCHEDULE keyword, kept verbatim - exact round trip, OPM Flow runs the same deck)
//   DEFINE name expression /    FU* field scalar, WU* one value per well, GU* one value per group
//   ASSIGN name [selector] value /   UPDATE name ON|OFF|NEXT /   UNITS name 'text' /
//   expression: numbers, summary quantities (F*, W* ['pattern' | well], G* ['group']), other UDQs, TIME / DAYS / DAY /
//   MNTH / YEAR, + - * / ^, comparisons == != < > <= >= (1 / 0), parentheses, elementwise functions ABS EXP LN LOG NINT
//   SQRT DEF UNDEF IDV SORTA SORTD, set functions SUM AVEA AVEG AVEH MAX MIN PROD NORM1 NORM2 NORMI (well / group set ->
//   scalar). Undefined values (wells outside a pattern, division by zero) propagate as NaN; set functions skip them.
//   Evaluated after every accepted time step in deck order (DEFINE with UPDATE ON every step, NEXT once, ASSIGN when its
//   report step starts).
// UDA  a UDQ name in a numeric item of WCONPROD (ORAT WRAT GRAT LRAT RESV BHP THP), WCONINJE (RATE RESV BHP THP) or
//   WELTARG: the target follows the UDQ (well UDQ: the well's value, field UDQ: the scalar) from the next time step on;
//   undefined values leave the previous target in force.
// ACTIONX  ACTIONX / name [max runs [min interval days]] / condition lines 'lhs op rhs [AND|OR] /' / '/' / body keywords /
//   ENDACTIO. lhs / rhs: numbers, summary quantities with optional well / group pattern, UDQs, DAY MNTH (JAN..DEC) YEAR
//   TIME; '(' ')' grouping; AND binds tighter than OR. A well-set comparison is true when at least one well satisfies it
//   and those wells replace '?' in the body (AND: intersection, OR: union of the well sets). Evaluated at the end of
//   every report step (OPM Flow practice); the body (WELOPEN, WCONPROD, WCONINJE, WELTARG, GCONPROD, COMPDAT on existing
//   connections, UDQ ...) acts from the next report step on.
#pragma once
#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace ressim {
namespace udq {

// what an expression can refer to (rates / pressures in case units, cumulatives in case surface units, time in days)
// C38 UDT user-defined table (1-D): TU_NAME[x] in UDQ expressions. Interpolation 'NV' nearest value, 'LC' linear with
// constant extrapolation, 'LL' linear with linear extrapolation (OPM Flow UDT semantics)
struct Table {
    std::string name, interp = "LC";
    std::vector<double> x, y;
    int step = 0;
    double eval(double v) const;
};
bool parseTables(const std::string& body, int step, std::vector<Table>& out, std::string& err);

struct Source {
    std::vector<std::string> wells, groups;
    std::function<double(const std::string& kw)> field;                  // F* vectors (NaN = unknown)
    std::function<double(const std::string& kw, size_t well)> well;       // W* vectors
    std::function<double(const std::string& kw, size_t group)> group;     // G* vectors
    double days = 0;                                                      // simulated time
    int day = 1, month = 1, year = 2000;                                  // calendar date
    const std::vector<Table>* tables = nullptr;                           // C38 UDT tables (set by Set::evaluate)
};

struct Value {                 // a scalar or one value per well / group (NaN = undefined)
    enum Kind { Scalar, Well, Group } kind = Scalar;
    double s = NAN;
    std::vector<double> v;
};

struct Node;
struct Expr { std::shared_ptr<Node> root; std::string text; };
// parses a UDQ expression (tokens of a DEFINE record after the name); false + message on a syntax error
bool parseExpr(const std::string& text, Expr& e, std::string& err);
using Lookup = std::function<const Value*(const std::string&)>;   // UDQ values by name
Value evaluate(const Expr& e, const Source& src, const Lookup& udqs);

bool wildMatch(const std::string& pattern, const std::string& name);   // ECLIPSE well patterns: * ? (case-insensitive)

struct Definition {
    std::string name;          // FU.. / WU.. / GU..
    char var = 'F';            // F W G
    int kind = 0;              // 0 DEFINE, 1 ASSIGN
    Expr expr;                 // DEFINE
    double value = 0;          // ASSIGN
    std::string selector;      // ASSIGN: well / group pattern ("" = all)
    int step = 0;              // report step the record appears at
    bool assigned = false;     // ASSIGN done
};

class Set {
public:
    // body of one UDQ keyword (verbatim deck text), active from report step 'step'; false + message on errors (valid
    // records are kept)
    bool add(const std::string& body, int step, std::string& err);
    bool addTables(const std::string& udtBody, int step, std::string& err);   // C38 UDT keyword
    const std::vector<Table>& tables() const { return tables_; }
    // evaluates the records active at report step 'step' (in deck order) with the current state
    // assignsOnly: only the ASSIGN records due at this step (start of a report step - DEFINEs keep their last values)
    void evaluate(const Source& src, int step, bool assignsOnly = false);
    const Value* get(const std::string& name) const;
    double scalar(const std::string& name) const;            // FU value (NaN undefined / unknown)
    double wellValue(const std::string& name, size_t w) const;   // WU value of well w; an FU name returns its scalar
    double groupValue(const std::string& name, size_t g) const;
    bool has(const std::string& name) const;
    std::vector<std::string> names() const;
    std::string units(const std::string& name) const;
    bool empty() const { return defs_.empty(); }
    const std::vector<Definition>& definitions() const { return defs_; }
    // restart: assignment flags, update states and current values
    std::string save() const;
    void load(const std::string& blob);
private:
    std::vector<Definition> defs_;
    std::vector<Table> tables_;
    std::vector<std::pair<std::string, int>> update_;          // name -> 0 ON, 1 OFF, 2 NEXT (latest active record)
    struct Upd { std::string name; int mode; int step; };
    std::vector<Upd> updRec_;
    std::vector<std::pair<std::string, std::string>> units_;
    std::vector<std::pair<std::string, Value>> vals_;
    std::vector<std::pair<std::string, int>> nextDone_;        // UPDATE NEXT consumed (name, step of the record)
    Value& slot(const std::string& name);
};

// ---------------------------------------------------------------- ACTIONX
struct CondNode;
struct Action {
    std::string name;
    int maxRuns = 1;
    double minWaitDays = 0;
    int fromStep = 0;              // report step of the ACTIONX keyword
    std::string conditionText;     // the condition lines (for the log / UI)
    std::string body;              // keywords between the conditions and ENDACTIO (deck text)
    std::shared_ptr<CondNode> cond;
    int runs = 0;
    double lastRun = -1e300;       // simulated time of the last run (days)
};
struct Trigger { bool fired = false; std::vector<size_t> wells; std::string why; };
// raw: the verbatim ACTIONX body (header record, condition lines, '/', keywords, optional ENDACTIO)
bool parseAction(const std::string& raw, int step, Action& a, std::string& err);
Trigger evaluateAction(const Action& a, const Source& src, const Set& udqs);
// body with the well-list marker '?' replaced by the wells that satisfied the condition (one record per well)
std::string expandBody(const std::string& body, const std::vector<std::string>& wells);

}  // namespace udq
}  // namespace ressim
