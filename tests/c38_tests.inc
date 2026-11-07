// C38 tests: production-data import / QC / WCONHIST generation and deck merge; float32 level-parallel BILU0; NONNC
#include "ressim/proddata.h"
#include "ressim/bsr.h"
#include <chrono>
namespace c38 {
inline void setEnv(const char* k, const char* v) {
#ifdef _WIN32
    _putenv_s(k, v ? v : "");
#else
    if (v) setenv(k, v, 1); else unsetenv(k);
#endif
}
inline std::string slurp(const std::string& p) { std::ifstream f(p, std::ios::binary); std::stringstream ss; ss << f.rdbuf(); return ss.str(); }
inline bool hasIssue(const std::vector<ressim::pdata::Issue>& v, const char* text, int sev) {
    for (const auto& i : v) if (i.severity == sev && i.msg.find(text) != std::string::npos) return true;
    return false;
}
}  // namespace c38

static void testC38ProdEnds();
static void testC38ProdData() {
    using namespace ressim;
    std::printf("== C38 production-data import: dates, QC, rates / WEFAC, WCONHIST merge into SPE1\n");
    // ---- dates
    {
        long d = 0; bool ok = true; int y, m, dd;
        auto is = [&](const char* s, bool mf, int Y, int M, int D) {
            long x = 0; if (!pdata::parseDate(s, mf, x)) return false;
            pdata::civilFromDays(x, y, m, dd); return y == Y && m == M && dd == D;
        };
        ok = ok && is("2020-02-29", false, 2020, 2, 29) && is("01-FEB-2020", false, 2020, 2, 1) && is("1 FEB 2020", false, 2020, 2, 1);
        ok = ok && is("02/01/2020", false, 2020, 1, 2) && is("02/01/2020", true, 2020, 2, 1) && is("31.12.2019", false, 2019, 12, 31);
        ok = ok && is("2020-03", false, 2020, 3, 1) && is("JAN-2021", false, 2021, 1, 1) && is("2020/07/15T00:00:00", false, 2020, 7, 15);
        ok = ok && !pdata::parseDate("31/04/2020", false, d) && !pdata::parseDate("2019-02-29", false, d) && !pdata::parseDate("hello", false, d);
        check(ok, "production table dates: ISO, DD-MON-YYYY, D MON YYYY, DMY / MDY, dotted, YYYY-MM, MON-YYYY, ISO with time; invalid dates rejected");
        check(pdata::columnRole("Oil (sm3)") == pdata::C_OIL && pdata::columnRole("WGPRH") == pdata::C_GAS && pdata::columnRole("Well Name") == pdata::C_WELL &&
              pdata::columnRole("hours_on") == pdata::C_HOURS && pdata::columnRole("BHP_bar") == pdata::C_BHP && pdata::columnRole("comment") == pdata::C_NONE,
              "production table column aliases (with unit suffixes) recognised, unknown columns ignored");
    }
    // ---- QC
    {
        const std::string tab =
            "WELL;DATE;OIL;WATER;GAS;DAYS\n"
            "A;01.01.2020;100;10;1000;31\n"
            "A;01.02.2020;-5;10;1000;29\n"
            "A;01.02.2020;50;0;0;\n"
            "A;xx;1;1;1;1\n"
            "A;01.03.2020;100;10;100000;40\n"
            "VERYLONGWELLNAME;01.03.2020;1;0;0;31\n";
        pdata::Options o; pdata::Table t;
        const bool ok = pdata::parseTable(tab, o, t);
        for (const auto& i : t.issues) std::printf("    [%d] line %d %s %s %s\n", i.severity, i.row, i.well.c_str(), i.date.c_str(), i.msg.c_str());
        const auto& I = t.issues;
        check(ok && t.nErrors() == 1 && c38::hasIssue(I, "negative", 1) && c38::hasIssue(I, "duplicate", 1) && c38::hasIssue(I, "clamped", 1) &&
              c38::hasIssue(I, "GOR jumps", 1) && c38::hasIssue(I, "longer than 8", 1) && c38::hasIssue(I, "not understood", 2),
              "QC: bad date dropped (error), negative value, duplicate (summed), uptime > period, GOR jump, long well name flagged");
        // duplicate volumes summed: Feb oil = 0 + 50 over 29 days on (WEFAC 1, Feb 2020 has 29 days); March: 31 d period
        double febO = -1, marEf = -1, marO = -1;
        for (const auto& x : t.wp) if (x.well == "A") { if (x.period == 1) febO = x.orat; if (x.period == 2) { marEf = x.wefac; marO = x.orat; } }
        std::printf("    Feb oil rate %.6g (50/29 = %.6g), March WEFAC %.4g rate %.6g\n", febO, 50.0 / 29, marEf, marO);
        check(std::fabs(febO - 50.0 / 29) < 1e-9 && std::fabs(marEf - 1) < 1e-12 && std::fabs(marO - 100.0 / 31) < 1e-9,
              "volumes -> rates: duplicates summed, uptime clamped to the period, last monthly period closed at the next month");
    }
    // ---- merge into SPE1 at the existing TSTEP dates (monthly) with uptime -> WEFAC; ORAT control -> FOPT = table volumes
    {
        const std::string deck = c38::slurp(findExample("SPE1.DATA"));
        const std::string tab =
            "WELL,DATE,OIL,GAS,GINJ,DAYS\n"
            "PROD,2015-01-01,62000,80000,0,31\n"
            "PROD,2015-02-01,54000,70000,0,14\n"
            "PROD,2015-03-01,55000,72000,0,31\n"
            "PROD,2015-04-01,48000,64000,0,30\n"
            "INJ,2015-02-01,0,0,1400000,14\n"
            "INJ,2015-03-01,0,0,3100000,31\n";
        pdata::Options o; o.prodControl = "ORAT"; o.startYear = 2015;
        pdata::Table t;
        bool ok = !deck.empty() && pdata::parseTable(tab, o, t);
        std::vector<pdata::Issue> mi;
        const std::string merged = ok ? pdata::mergeIntoDeck(deck, t, o, mi, {"PROD", "INJ"}) : std::string();
        for (const auto& i : mi) std::printf("    merge [%d] %s %s\n", i.severity, i.well.c_str(), i.msg.c_str());
        ok = ok && merged.find("WCONHIST") != std::string::npos && merged.find("WCONINJH") != std::string::npos && merged.find("WEFAC") != std::string::npos &&
             merged.find("WCONPROD") == std::string::npos && merged.find("WCONINJE") == std::string::npos && c38::hasIssue(mi, "4 period(s) at existing", 0);
        check(ok, "merge: history appended at the deck's own TSTEP dates, WCONPROD / WCONINJE of the imported wells removed, WEFAC for downtime");
        double fopt120 = -1, fgit = -1; size_t nsteps = 0;
        try {
            std::vector<DeckMessage> rep;
            const Case c = eclipseDeckFromString(merged, "", nullptr, &rep);
            nsteps = c.schedule.size();
            Simulator s; s.setup(c); s.initialize();
            while (s.reportStep() < 4 && s.runReportStep()) {}
            for (const auto& p : s.summary()) if (std::fabs(p.time - 120) < 1e-6) { fopt120 = p.fopt; fgit = p.fgit; }
        } catch (const std::exception& e) { std::printf("    %s\n", e.what()); }
        const double want = 62000 + 54000 + 55000 + 48000;
        std::printf("    merged SPE1: %zu report steps, FOPT(120 d) %.8g (table %.8g), FGIT %.8g (table %.8g)\n", nsteps, fopt120, want, fgit, 4.5e6);
        check(nsteps == 120 && std::fabs(fopt120 - want) < 1e-3 * want, "merged SPE1 runs; ORAT history with WEFAC 0.5 in February reproduces the table oil volume");
        check(std::fabs(fgit - 4.5e6) < 2e-3 * 4.5e6, "WCONINJH gas injection with WEFAC: injected volume = table volume");
    }
    // ---- dates between the deck's report dates -> new DATES; the schedule is extended to the end of the last period
    {
        const std::string deck = c38::slurp(findExample("SPE1.DATA"));
        std::string shortDeck = deck;
        const size_t p = shortDeck.find("TSTEP");
        if (p != std::string::npos) { const size_t e = shortDeck.find('/', p); shortDeck.replace(p, e - p + 1, "TSTEP\n 31 28 /"); }
        const std::string tab = "WELL\tDATE\tOIL\nPROD\t15-JAN-2015\t31000\nPROD\t15-FEB-2015\t28000\nPROD\t15-MAR-2015\t31000\n";
        pdata::Options o; o.prodControl = "ORAT";
        pdata::Table t; std::vector<pdata::Issue> mi;
        bool ok = pdata::parseTable(tab, o, t);
        const std::string merged = ok ? pdata::mergeIntoDeck(shortDeck, t, o, mi) : std::string();
        size_t nsteps = 0; double tEnd = 0, fopt = 0;
        try {
            const Case c = eclipseDeckFromString(merged);
            nsteps = c.schedule.size();
            for (const auto& r : c.schedule) tEnd += r.length;
            Simulator s; s.setup(c); s.initialize(); while (s.runReportStep()) {}
            fopt = s.summary().back().fopt;
        } catch (const std::exception& e) { std::printf("    %s\n", e.what()); }
        // boundaries: 14 (15 Jan), 31, 45 (15 Feb), 59, 73 (15 Mar), 104 (15 Apr = end of the last monthly period)
        std::printf("    %zu report steps, schedule length %.6g d, FOPT %.8g (table 90000)\n", nsteps, tEnd, fopt);
        check(ok && nsteps == 6 && std::fabs(tEnd - 104) < 1e-9 && std::fabs(fopt - 90000) < 1e-3 * 90000,
              "history dates between report dates: new DATES inserted, schedule extended to the end of the last period, oil volume honoured");
    }
    // ---- stand-alone include + observed CSV
    {
        pdata::Options o; o.startYear = 2015; pdata::Table t;
        pdata::parseTable("WELL,DATE,OIL,WATER,BHP\nP1,2015-01-01,3100,0,250\nP1,2015-02-01,2800,280,240\n", o, t);
        const std::string inc = pdata::historyInclude(t, o), obs = pdata::observedCsv(t, o);
        const bool ok = inc.find("DATES\n  1 'FEB' 2015 /") != std::string::npos && inc.find("DATES\n  1 'MAR' 2015 /") != std::string::npos &&
                        inc.find("DATES\n  1 'JAN' 2015") == std::string::npos && inc.find("'P1' OPEN RESV 100 0 0 1* 1* 1* 250 /") != std::string::npos &&
                        obs.find("P1,2015-03-01,59,100,10,0,0,0,240") != std::string::npos;
        if (!ok) std::printf("%s\n%s\n", inc.c_str(), obs.c_str());
        check(ok, "history include (no DATES at START, RESV + observed BHP) and observed CSV (period-end dates, calendar rates)");
    }
    testC38ProdEnds();
}

// C38 end-of-period dates are inclusive (31 JAN = all of January) and an implausible GOR (scf in a FIELD deck) is flagged
static void testC38ProdEnds() {
    using namespace ressim;
    const std::string csv = "WELL,DATE,OIL,GAS\nPROD,2015-01-31,620000,800000000\nPROD,2015-02-28,560000,700000000\n";
    pdata::Options o; o.datesAreEnds = true;
    pdata::Table t; const bool ok = pdata::parseTable(csv, o, t);
    bool per = ok && t.periods.size() == 2;
    if (per) per = pdata::dateString(t.periods[0].d0) == pdata::dateString(pdata::daysFromCivil(2015, 1, 1)) &&
                   t.periods[0].d1 == pdata::daysFromCivil(2015, 2, 1) && t.periods[1].d1 == pdata::daysFromCivil(2015, 3, 1);
    std::vector<pdata::Issue> iss;
    const std::string deck = c38::slurp(findExample("SPE1.DATA"));
    pdata::mergeIntoDeck(deck, t, o, iss);
    bool gorWarn = false, early = false;
    for (auto& i : iss) { gorWarn = gorWarn || i.msg.find("implausible") != std::string::npos; early = early || i.msg.find("before START") != std::string::npos; }
    pdata::Options o2 = o; o2.gasFactor = 0.001; pdata::Table t2; pdata::parseTable(csv, o2, t2);
    std::vector<pdata::Issue> iss2; pdata::mergeIntoDeck(deck, t2, o2, iss2);
    bool warn2 = false; for (auto& i : iss2) warn2 = warn2 || i.msg.find("implausible") != std::string::npos;
    std::printf("    end dates: first period %s .. %s; GOR warning scf %d / Mscf %d\n", per ? pdata::dateString(t.periods[0].d0).c_str() : "?",
                per ? pdata::dateString(t.periods[0].d1).c_str() : "?", gorWarn, warn2);
    check(per && !early, "production dates as period ends: 31 JAN covers 1 JAN .. 1 FEB 00:00 (no shift before START)");
    check(gorWarn && !warn2, "production import: gas in scf for a FIELD deck flagged (GOR implausible), Mscf accepted");
}

static void testC38FloatIlu() {
    using namespace ressim;
    std::printf("== C38 float32 factors in the level-parallel (contiguous) BILU0\n");
    const int nx = 30, ny = 30, nz = 8, nb = 3, n = nx * ny * nz;
    BSR A; A.n = n; A.nb = nb; A.rowPtr.assign(size_t(n) + 1, 0);
    unsigned long long seed = 12345;
    auto rnd = [&]() { seed = seed * 6364136223846793005ULL + 1442695040888963407ULL; return double((seed >> 11) & 0xFFFFF) / double(0xFFFFF); };
    auto id = [&](int i, int j, int k) { return (k * ny + j) * nx + i; };
    for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        std::vector<int> cols;
        if (k > 0) cols.push_back(id(i, j, k - 1));
        if (j > 0) cols.push_back(id(i, j - 1, k));
        if (i > 0) cols.push_back(id(i - 1, j, k));
        cols.push_back(id(i, j, k));
        if (i < nx - 1) cols.push_back(id(i + 1, j, k));
        if (j < ny - 1) cols.push_back(id(i, j + 1, k));
        if (k < nz - 1) cols.push_back(id(i, j, k + 1));
        const int r = id(i, j, k);
        for (int c : cols) {
            A.col.push_back(c);
            double blk[9];
            for (int e = 0; e < 9; ++e) blk[e] = c == r ? 0.3 * (rnd() - 0.5) : -(e % 4 == 0 ? 1.0 + rnd() : 0.2 * (rnd() - 0.5)) * (c / (nx * ny) != r / (nx * ny) ? 0.1 : 1.0);
            if (c == r) for (int e = 0; e < 3; ++e) blk[e * 4] += 7.5 + rnd();
            A.val.insert(A.val.end(), blk, blk + 9);
        }
        A.rowPtr[size_t(r) + 1] = int(A.col.size());
    }
    A.finalizeDiag();
    std::vector<double> rhs(size_t(n) * nb), z1(rhs.size()), z2(rhs.size());
    for (auto& x : rhs) x = rnd() - 0.5;
    c38::setEnv("RESSIM_LEVEL_ILU", "2");
    BILU0 d, f; f.setFloatFactors(true);
    d.factor(A); f.factor(A);
    const bool contig = d.contiguousLevels() && f.contiguousLevels() && f.floatFactors() && !d.floatFactors();
    auto time = [&](const BILU0& M, std::vector<double>& z) {
        const auto t0 = std::chrono::steady_clock::now();
        for (int r = 0; r < 20; ++r) M.apply(rhs.data(), z.data());
        return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count() / 20;
    };
    const double td = time(d, z1), tf = time(f, z2);
    double num = 0, den = 0; for (size_t i = 0; i < z1.size(); ++i) { num += (z1[i] - z2[i]) * (z1[i] - z2[i]); den += z1[i] * z1[i]; }
    const double rel = std::sqrt(num / std::max(den, 1e-300));
    std::printf("    %d cells x 3: contiguous %d, factors %zu B (double) vs %zu B (float), apply %.3f ms vs %.3f ms, |dz|/|z| %.3g\n", n, int(contig),
                d.factorBytes(), f.factorBytes(), td, tf, rel);
    check(contig && f.factorBytes() * 2 == d.factorBytes(), "float BILU0: level-ordered factors kept in float32 (half the bytes)");
    check(rel < 1e-5, "float BILU0: preconditioner output equals the double one to float rounding");
    // refactorisation with new values (same pattern) keeps working (the double work array is re-created)
    for (auto& v : A.val) v *= 2.0;
    f.factor(A); d.factor(A);
    f.apply(rhs.data(), z2.data()); d.apply(rhs.data(), z1.data());
    num = den = 0; for (size_t i = 0; i < z1.size(); ++i) { num += (z1[i] - z2[i]) * (z1[i] - z2[i]); den += z1[i] * z1[i]; }
    check(std::sqrt(num / den) < 1e-5 && f.floatFactors(), "float BILU0: re-factorisation of new values with the same pattern");
    c38::setEnv("RESSIM_LEVEL_ILU", nullptr);
}

static void testC38Nonnc() {
    using namespace ressim;
    std::printf("== C38 NONNC: no fault-throw non-neighbour connections\n");
    const std::string base = c38::slurp(findExample("FAULTED.DATA"));
    std::string d2 = base;
    const size_t g = d2.find("\nGRID");
    bool ok = g != std::string::npos;
    if (ok) d2.insert(d2.find('\n', g + 1) + 1, "NONNC\n");
    size_t c1 = 0, c2 = 0; int nn = 0; bool exp = false, run = false;
    try {
        const Case a = eclipseDeckFromString(base), b = eclipseDeckFromString(d2);
        Grid G1, G2; G1.build(a); G2.build(b);
        c1 = G1.conns.size(); c2 = G2.conns.size(); nn = G2.numNonnc;
        exp = eclipseDeckToString(b).find("\nNONNC\n") != std::string::npos;
        Simulator s; s.setup(b); s.initialize(); run = s.runReportStep();
    } catch (const std::exception& e) { std::printf("    %s\n", e.what()); ok = false; }
    std::printf("    connections %zu -> %zu, NONNC dropped %d\n", c1, c2, nn);
    check(ok && nn > 0 && c2 < c1 && exp && run, "NONNC: fault-throw connections dropped, keyword exported, the model runs");
}

static void testC38Udt() {
    using namespace ressim;
    std::printf("== C38 UDT tables, UDQ union operators, OPM evaluation order\n");
    std::vector<udq::Table> T; std::string err;
    const bool ok = udq::parseTables(" 'TU_FBHP' 1 /\n 'LC'  100.0  500.0 / -- x\n 100.0  180.0 / -- y\n/\n 'TU_B' 1 /\n 'LL' 0 10 20 /\n 0 5 20 /\n/\n/\n", 0, T, err);
    bool vals = ok && T.size() == 2;
    if (vals) {
        const udq::Table& a = T[0]; const udq::Table& b = T[1];
        vals = std::fabs(a.eval(368.337) - 153.6674) < 1e-3 && a.eval(50) == 100 && a.eval(900) == 180 &&
               std::fabs(b.eval(15) - 12.5) < 1e-12 && std::fabs(b.eval(30) - 35) < 1e-12 && std::fabs(b.eval(-10) + 5) < 1e-12;
        udq::Table c = a; c.interp = "NV";
        vals = vals && c.eval(250) == 100 && c.eval(350) == 180;
    }
    std::printf("    %zu tables %s\n", T.size(), err.c_str());
    check(vals, "UDT: two tables parsed; LC clamps, LL extrapolates linearly, NV nearest value");
    // union operators: an undefined operand gives the other one
    udq::Set s;
    std::string e2;
    s.addTables(" 'TU_FBHP' 1 /\n 'LC' 100 500 /\n 100 180 /\n/\n/\n", 0, e2);
    const bool added = s.add(" ASSIGN FU_X 0 /\n DEFINE FU_X0 FU_X /\n DEFINE FU_X (TU_FBHP[FOPR] UMIN WBHP 'P1') UMIN FU_X0 /\n"
                             " DEFINE WU_Y TU_FBHP[WOPR] UMAX 120 /\n DEFINE FU_U FU_NOPE UADD 7 /\n/\n", 0, e2);
    udq::Source src; src.wells = {"P1", "P2"};
    double fopr = 500, bhp = 190;
    src.field = [&](const std::string& k) { return k == "FOPR" ? fopr : NAN; };
    src.well = [&](const std::string& k, size_t w) { return k == "WBHP" ? bhp : k == "WOPR" ? (w == 0 ? 100.0 : 400.0) : NAN; };
    s.evaluate(src, 0, true); s.evaluate(src, 0, false);
    const double x1 = s.scalar("FU_X"), x01 = s.scalar("FU_X0");
    fopr = 368.337; bhp = 180;
    s.evaluate(src, 0, false);
    const double x2 = s.scalar("FU_X");
    std::printf("    FU_X %.6g (FU_X0 %.6g) -> %.6g; WU_Y %.4g %.4g; FU_U %.4g %s\n", x1, x01, x2, s.wellValue("WU_Y", 0), s.wellValue("WU_Y", 1), s.scalar("FU_U"), e2.c_str());
    check(added && std::fabs(x1 - 180) < 1e-9 && std::fabs(x01 - 180) < 1e-9 && std::fabs(x2 - 153.6674) < 1e-3,
          "OPM order: ASSIGN first, then DEFINEs by first appearance - X sees the undefined X0 (UMIN ignores it): 180, then the table value 153.67");
    check(std::fabs(s.wellValue("WU_Y", 0) - 120) < 1e-9 && std::fabs(s.wellValue("WU_Y", 1) - 160) < 1e-9 && std::fabs(s.scalar("FU_U") - 7) < 1e-12,
          "UMAX on well sets, UADD with an undefined operand returns the defined one");
    // the deck keeps UDT (engine-evaluated) and writes UDTDIMS
    {
        std::string d = c38::slurp(findExample("SPE1.DATA"));
        const size_t p = d.find("TSTEP");
        d.insert(p, "UDT\n 'TU_LIM' 1 /\n 'LC' 0 40000 /\n 4000 1500 /\n/\n/\nUDQ\n DEFINE WU_BHP TU_LIM[WOPR] /\n/\n");
        bool ok2 = false; size_t nudt = 0;
        try {
            const Case c = eclipseDeckFromString(d);
            for (const auto& e : c.extraKeywords) nudt += e.name == "UDT";
            const std::string out = eclipseDeckToString(c);
            ok2 = nudt == 1 && out.find("UDTDIMS\n 1 2 2 1 /") != std::string::npos && out.find("'TU_LIM'") != std::string::npos;
            Simulator sim; sim.setup(c); sim.initialize(); sim.runReportStep();
            double v = NAN; for (const auto& [n, x] : sim.summary().back().udq) if (n == "WU_BHP:PROD") v = x;
            const double q = sim.summary().back().wells[0].oil;
            std::printf("    SPE1 + UDT: WU_BHP:PROD %.6g for WOPR %.6g (expected %.6g)\n", v, q, 4000 - 2500 * q / 40000);
            ok2 = ok2 && std::fabs(v - (4000 - 2500 * q / 40000)) < 1e-6 * 4000;
        } catch (const std::exception& e) { std::printf("    %s\n", e.what()); }
        check(ok2, "UDT kept by the importer, UDTDIMS written by the exporter, TU_LIM[WOPR] evaluated in the run");
    }
}

// C38 EDIT-section MULTIPLY on TRANX (neighbour connections only, fault NNCs untouched, exported) and SWATINIT below SWL
static void testC38TranEdit() {
    using namespace ressim;
    std::printf("== C38 EDIT TRANX edits and SWATINIT < SWL\n");
    const std::string base = c38::slurp(findExample("FAULTED.DATA"));
    std::string d2 = base;
    const size_t pr = d2.find("\nPROPS");
    bool ok = pr != std::string::npos;
    if (ok) d2.insert(pr + 1, "EDIT\nMULTIPLY\n 'TRANX' 0.5 /\n/\n");
    int nNb = 0, nNbOk = 0, nNnc = 0, nNncOk = 0; bool exp = false, rt = false;
    try {
        const Case a = eclipseDeckFromString(base), b = eclipseDeckFromString(d2);
        Grid G1, G2; G1.build(a); G2.build(b);
        std::map<std::pair<int, int>, double> t1;
        for (auto& c : G1.conns) t1[{std::min(c.a, c.b), std::max(c.a, c.b)}] = c.trans;
        const int nx = a.grid.nx, nl = a.grid.nx * a.grid.ny;
        for (auto& c : G2.conns) {
            if (c.dir != 'X') continue;
            auto it = t1.find({std::min(c.a, c.b), std::max(c.a, c.b)});
            if (it == t1.end()) continue;
            int ga = G2.globalOf[size_t(c.a)], gb = G2.globalOf[size_t(c.b)];
            if (ga > gb) std::swap(ga, gb);
            const bool nb = gb - ga == 1 && ga % nx + 1 < nx;
            const double r = c.trans / it->second;
            if (nb) { ++nNb; nNbOk += std::fabs(r - 0.5) < 1e-9; }
            else { ++nNnc; nNncOk += std::fabs(r - 1.0) < 1e-9; }
        }
        (void)nl;
        const std::string ex = eclipseDeckToString(b);
        exp = ex.find("\nEDIT\n") != std::string::npos && ex.find("'TRANX' 0.5") != std::string::npos;
        const Case b2 = eclipseDeckFromString(ex);
        Grid G3; G3.build(b2);
        rt = G3.conns.size() == G2.conns.size();
        for (size_t q = 0; rt && q < G3.conns.size(); ++q) rt = std::fabs(G3.conns[q].trans - G2.conns[q].trans) <= 1e-9 * std::max(1e-12, G2.conns[q].trans);
    } catch (const std::exception& e) { std::printf("    %s\n", e.what()); ok = false; }
    std::printf("    X neighbours halved %d/%d, fault NNCs across X faces unchanged %d/%d, export %d, round trip %d\n", nNbOk, nNb, nNncOk, nNnc, exp, rt);
    check(ok && nNb > 0 && nNbOk == nNb && nNnc > 0 && nNncOk == nNnc, "EDIT MULTIPLY TRANX: Cartesian neighbours scaled, fault juxtaposition NNCs untouched");
    check(exp && rt, "EDIT MULTIPLY TRANX: written as an EDIT section and re-imported to identical transmissibilities");

    // SWATINIT below the connate saturation: Sw = SWL (OPM applySwatinit), Pc scaled with Pc(SWL)
    bool ok2 = false;
    try {
        Case c = eclipseDeckFromString(c38::slurp(findExample("SPE1.DATA")));
        double swl = 1;
        for (auto& r : c.props.swof) swl = std::min(swl, r.s);
        for (auto& r : c.props.swof) { const double se = (r.s - swl) / std::max(1e-9, c.props.swof.back().s - swl); r.pc = 3.0 * std::pow(1 - se, 3); }
        c.grid.swatinit.assign(size_t(c.grid.size()), 0.5 * swl);
        Simulator s; s.log = [](const std::string&) {}; s.setup(c); s.initialize();
        const auto& sn = s.snapshots();
        double smin = 1, smax = 0;
        for (float v : sn[0].sw) { smin = std::min(smin, double(v)); smax = std::max(smax, double(v)); }
        std::printf("    SWATINIT %.4f (SWL %.4f): initial Sw %.4f..%.4f\n", 0.5 * swl, swl, smin, smax);
        ok2 = std::fabs(smin - swl) < 1e-6;
    } catch (const std::exception& e) { std::printf("    %s\n", e.what()); }
    check(ok2, "SWATINIT < SWL: the initial water saturation is raised to SWL");
}
