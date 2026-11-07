// C38 production-data workflow: import of measured well data (monthly volumes / daily rates / pressures from a CSV,
// TSV or whitespace table, e.g. a production-database export), quality control, and generation of the history
// SCHEDULE (WCONHIST / WCONINJH / WEFAC at DATES) merged into an ECLIPSE / OPM deck.
//
//   table      header row with column names (case-insensitive aliases, see columnRole()), one row per well and date;
//              dates as 2020-01-31, 2020/01/31, 31.01.2020, 31/01/2020 (day first; option monthFirst), 01-JAN-2020,
//              1 JAN 2020, 2020-01 (month = day 1), or separate YEAR / MONTH / DAY columns
//   periods    the unique dates define the history periods; a row is either the START of its period (default: rates
//              from this date to the next date) or the END (datesAreEnds: the period from the previous date)
//   volumes    volumes per period (volumes = true) are turned into rates: calendar-average rate = V / period length,
//              or with uptime (DAYS / HOURS column) and useWefac: flowing rate = V / days on and WEFAC = days on / period
//   QC         unparseable rows, negative values (set to 0), duplicates (summed for volumes, averaged for rates), uptime
//              longer than the period (clamped), water cut / GOR outliers (jump x10 against the previous period of the
//              well), producing and injecting in the same period, wells missing in some periods (shut there), wells
//              not defined in the deck (merge)
//   merge      the deck's SCHEDULE is split into time blocks at its DATES / TSTEP; the history keywords are appended to
//              the block starting at the same time or to a new block (DATES for whole days, TSTEP otherwise); the
//              schedule is extended to the end of the last period. Optionally the existing WCONPROD / WCONHIST /
//              WCONINJE / WCONINJH records of the imported wells are removed so that the history controls them alone.
#pragma once
#include <string>
#include <vector>

namespace ressim {
namespace pdata {

enum Col { C_NONE = -1, C_WELL, C_DATE, C_YEAR, C_MONTH, C_DAY, C_OIL, C_WATER, C_GAS, C_WINJ, C_GINJ, C_BHP, C_THP,
           C_DAYS, C_HOURS, C_PRES, C_NCOL };
const char* colName(int role);
int columnRole(const std::string& header);          // alias -> role (C_NONE if unknown)

struct Record {
    std::string well;
    long day = 0;                                     // days since 1970-01-01 (civil)
    double v[C_NCOL];                                 // values by role (NaN = missing)
    int row = 0;                                      // source line (1-based)
    Record();
};

struct Issue {
    int severity = 0;                                 // 0 info, 1 warning, 2 error (row dropped)
    int row = 0;                                      // source line, 0 = general
    std::string well, date, msg;
};

struct Options {
    bool volumes = true;                              // rate columns hold volumes per period (false: rates)
    bool datesAreEnds = false;                        // the date closes the period, inclusive: 31 JAN = up to 1 FEB 00:00 (false: opens it)
    bool monthFirst = false;                          // 01/02/2020 = 2 Jan (US) instead of 1 Feb
    bool useWefac = true;                             // with an uptime column: flowing rate + WEFAC
    bool replaceControls = true;                      // drop the deck's WCONPROD / WCONHIST / WCONINJE / WCONINJH of imported wells
    bool shutMissing = true;                          // a well absent in a period (after its first record) is shut there
    bool writeBhp = true, writeThp = true;            // measured BHP / THP into WCONHIST / WCONINJH (observed values)
    std::string prodControl = "RESV";                 // WCONHIST control: ORAT WRAT GRAT LRAT RESV BHP
    double liqFactor = 1, gasFactor = 1;              // unit factors applied to liquid / gas columns
    double outlierFactor = 10;                        // WCT / GOR jump factor flagged as outlier
    int startYear = 0, startMonth = 1, startDay = 1;  // deck START when it is not found in the deck text
    std::vector<int> roleOverride;                    // per column: -2 = automatic (columnRole), C_NONE = ignore, else the role
};

struct Period { long d0 = 0, d1 = 0; };               // [d0, d1) days since epoch

struct WellPeriod {                                   // one well in one period, after QC
    std::string well;
    int period = 0;
    bool inj = false, shut = false;
    double orat = 0, wrat = 0, grat = 0, winj = 0, ginj = 0;   // rates (flowing rates with WEFAC)
    double bhp = 0, thp = 0, pres = 0;                // 0 = not measured
    double wefac = 1;
};

struct Table {
    std::vector<std::string> columns;                 // header as read
    std::vector<int> roles;                           // role per column
    std::vector<Record> recs;                         // parsed rows (sorted by well, date)
    std::vector<std::string> wells;                   // in first-appearance order
    std::vector<long> dates;                          // unique sorted dates
    std::vector<Period> periods;
    std::vector<WellPeriod> wp;                       // processed (per well, per period)
    std::vector<Issue> issues;
    int nErrors() const, nWarnings() const;
};

// parse + QC + rates; returns false only if nothing usable was found (Table::issues tells why)
bool parseTable(const std::string& text, const Options& o, Table& t);
// the history keywords of one period (DATES not included)
std::string periodKeywords(const Table& t, const Options& o, int period);
// merge into a deck text; deckWells = wells defined in the deck (empty = no check); returns the new deck text
std::string mergeIntoDeck(const std::string& deck, const Table& t, const Options& o, std::vector<Issue>& issues,
                          const std::vector<std::string>& deckWells = {});
// stand-alone history include file (DATES + keywords for every period; the deck START must precede the first date)
std::string historyInclude(const Table& t, const Options& o);
// observed series for history matching: CSV "WELL,DATE,DAYS,WOPRH,WWPRH,WGPRH,WWIRH,WGIRH,WBHPH" (calendar rates)
std::string observedCsv(const Table& t, const Options& o);

long daysFromCivil(int y, int m, int d);
void civilFromDays(long z, int& y, int& m, int& d);
std::string dateString(long day);                    // 1 'JAN' 2020
bool parseDate(const std::string& s, bool monthFirst, long& day);

}  // namespace pdata
}  // namespace ressim
