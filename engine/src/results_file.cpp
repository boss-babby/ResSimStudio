// C14 results + restart file (.rsres), see results_file.h
#include "ressim/results_file.h"
#include "ressim/binio.h"
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace ressim {
namespace {

const char kMagic[8] = {'R', 'S', 'R', 'E', 'S', 0x1a, 0, 1};
const uint32_t kVersion = 2;      // 2 (C18): RunStats stored size-prefixed (forward / backward compatible)
uint32_t g_readVer = kVersion;
// C17 files (version 1) stored the RunStats struct raw, which then ended at gpuFallbacks
constexpr size_t kStatsV1 = offsetof(RunStats, gpuWallSec);

template <class Ar> void wellRatesIO(Ar& ar, WellRates& r) {
    ar(r.oil, r.water, r.gas, r.waterInj, r.gasInj, r.resv, r.polyInj, r.polyProd, r.bhp, r.oilH, r.waterH, r.gasH, r.mode, r.open);
    ar(r.comp, r.compCumP, r.compCumI);
}
template <class Ar> void summaryIO(Ar& ar, SummaryPoint& s) {
    ar(s.time, s.dt, s.newton, s.linear);
    ar(s.fopr, s.fwpr, s.fgpr, s.fwir, s.fgir, s.fopt, s.fwpt, s.fgpt, s.fwit, s.fgit, s.fpr, s.fprp, s.fwct, s.fgor);
    ar(s.oip, s.wip, s.gip, s.faqr, s.faqt, s.faqrS, s.faqtS, s.fcir, s.fcpr, s.fcit, s.fcpt, s.fcip, s.fcad);
    ar(s.co2Mob, s.co2Trap, s.co2Diss, s.co2Struct, s.co2Min, s.co2KrImmob, s.saltSolid, s.co2Oil);
    ar(s.favgT, s.fmaxT, s.fsteamPV, s.fhlt, s.fenergyIn, s.csor, s.freactRate, s.freactCum, s.fuplift, s.fupliftRes, s.fmaxFail);
    ar(s.fcompP, s.fcompI, s.fcompPT, s.fcompIT);
    uint64_t nw = s.wells.size();
    ar(nw);
    if (Ar::reading) s.wells.resize(size_t(nw));
    for (WellRates& w : s.wells) wellRatesIO(ar, w);
}
template <class Ar> void snapIO(Ar& ar, Snapshot& s) {
    ar(s.time, s.pressure, s.sw, s.sg, s.so, s.rs, s.rv, s.poly, s.polyAds, s.xco2);
    uint64_t ne = s.extra.size();
    ar(ne);
    if (Ar::reading) s.extra.resize(size_t(ne));
    for (auto& e : s.extra) ar(e.first, e.second);
}
template <class Ar> void metaIO(Ar& ar, ResultsFile& r) {
    ar(r.engine, r.title, r.created, r.wellNames, r.wellTypes, r.trackName, r.totalTime);
    ar(r.initial.oil, r.initial.water, r.initial.gasFree, r.initial.gasDissolved, r.initial.oilVaporized);
    if constexpr (Ar::reading) {
        uint64_t ns = kStatsV1;
        if (g_readVer >= 2) ar(ns);
        std::string raw(size_t(ns), '\0');
        if (ns) ar.r.raw(&raw[0], size_t(ns));
        r.stats = RunStats{};
        std::memcpy(&r.stats, raw.data(), std::min<size_t>(size_t(ns), sizeof(RunStats)));
    } else {
        uint64_t ns = sizeof(RunStats);
        ar(ns);
        ar.w.raw(&r.stats, sizeof(RunStats));
    }
    uint64_t nc = r.components.size();
    ar(nc);
    if (Ar::reading) r.components.resize(size_t(nc));
    for (ComponentInfo& c : r.components) ar(c.name, c.mw, c.kind, c.species, c.hhv, c.lhv);
}

void chunk(std::ofstream& f, const char* tag, const std::string& payload) {
    char t[8]; std::memset(t, ' ', 8); std::memcpy(t, tag, std::min<size_t>(8, std::strlen(tag)));
    const uint64_t n = payload.size();
    f.write(t, 8); f.write(reinterpret_cast<const char*>(&n), 8); f.write(payload.data(), std::streamsize(n));
}

ResultsFile readFile(const std::string& path, bool snaps, bool states) {
    std::ifstream f(std::filesystem::u8path(path), std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    char m[8]; f.read(m, 8);
    if (!f || std::memcmp(m, kMagic, 8) != 0) throw std::runtime_error(path + " is not a ResSim results file (.rsres)");
    uint32_t ver = 0; f.read(reinterpret_cast<char*>(&ver), 4);
    if (ver > kVersion) throw std::runtime_error(path + " was written by a newer ResSim (format version " + std::to_string(ver) + ")");
    g_readVer = ver;
    ResultsFile r;
    std::string buf;
    for (;;) {
        char t[8]; uint64_t n = 0;
        f.read(t, 8); if (!f) break;
        f.read(reinterpret_cast<char*>(&n), 8); if (!f) throw std::runtime_error(path + ": truncated chunk header");
        const std::string tag(t, 8);
        const bool want = tag == "CASE    " || tag == "META    " || tag == "SUMMARY " || tag == "BLOB    " || tag == "RESTART " ||
                          (snaps && tag == "SNAP    ");
        if (!want) { f.seekg(std::streamoff(n), std::ios::cur); continue; }
        buf.resize(size_t(n));
        f.read(&buf[0], std::streamsize(n));
        if (!f) throw std::runtime_error(path + ": truncated chunk " + tag);
        StateIn in(buf);
        if (tag == "CASE    ") r.caseText = in.r.str();
        else if (tag == "META    ") metaIO(in, r);
        else if (tag == "SUMMARY ") {
            const uint64_t ns = in.r.u64();
            r.summary.resize(size_t(ns));
            for (SummaryPoint& s : r.summary) summaryIO(in, s);
        } else if (tag == "SNAP    ") { r.snaps.emplace_back(); snapIO(in, r.snaps.back()); }
        else if (tag == "BLOB    ") { std::string k = in.r.str(); std::string v = in.r.str(); r.blobs.emplace_back(std::move(k), std::move(v)); }
        else if (tag == "RESTART ") {
            ResultsFile::Restart R;
            uint64_t sc = 0;
            in(R.report, R.time, R.snap, sc, R.kind);
            R.summaryCount = size_t(sc);
            if (states) R.state = in.r.str();
            r.restarts.push_back(std::move(R));
        }
        if (tag == "END     ") break;
    }
    return r;
}

} // namespace

void saveResultsFile(const ResultsFile& r0, const std::string& path) {
    ResultsFile& r = const_cast<ResultsFile&>(r0);   // the IO visitors are symmetric (non-const); nothing is modified
    const std::string tmp = path + ".part";
    {
        std::ofstream f(std::filesystem::u8path(tmp), std::ios::binary | std::ios::trunc);
        if (!f) throw std::runtime_error("cannot write " + path);
        f.write(kMagic, 8);
        f.write(reinterpret_cast<const char*>(&kVersion), 4);
        { StateOut o; o.w.str(r.caseText); chunk(f, "CASE", o.w.buf); }
        { StateOut o; metaIO(o, r); chunk(f, "META", o.w.buf); }
        { StateOut o; o.w.u64(r.summary.size()); for (SummaryPoint& s : r.summary) summaryIO(o, s); chunk(f, "SUMMARY", o.w.buf); }
        for (Snapshot& s : r.snaps) { StateOut o; snapIO(o, s); chunk(f, "SNAP", o.w.buf); }
        for (const auto& b : r.blobs) { StateOut o; o.w.str(b.first); o.w.str(b.second); chunk(f, "BLOB", o.w.buf); }
        for (ResultsFile::Restart& R : r.restarts) {
            StateOut o; uint64_t sc = R.summaryCount;
            o(R.report, R.time, R.snap, sc, R.kind);
            o.w.str(R.state);
            chunk(f, "RESTART", o.w.buf);
        }
        chunk(f, "END", std::string());
        if (!f) throw std::runtime_error("write error on " + path);
    }
    std::error_code ec;
    std::filesystem::rename(std::filesystem::u8path(tmp), std::filesystem::u8path(path), ec);
    if (ec) {   // Windows: rename over an existing file can fail -> remove and retry
        std::filesystem::remove(std::filesystem::u8path(path), ec);
        std::filesystem::rename(std::filesystem::u8path(tmp), std::filesystem::u8path(path), ec);
        if (ec) throw std::runtime_error("cannot replace " + path + ": " + ec.message());
    }
}

ResultsFile loadResultsFile(const std::string& path, bool loadSnapshots) { return readFile(path, loadSnapshots, true); }
ResultsFile peekResultsFile(const std::string& path) { return readFile(path, false, false); }

} // namespace ressim
