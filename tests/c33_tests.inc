// C33 tests: scripting API (api.h Session), report building blocks (DEFLATE / ZIP / PNG / DOCX / PPTX), the regression
// harness (regress.h) and the engine items found by the OPM opm-tests regression: WEFAC efficiency factors, WLIST well
// lists with shell-style patterns, wells created at run time (WELSPECS + COMPDAT in a schedule script / ACTIONX).
#include "ressim/api.h"
#include "ressim/report.h"
#include "ressim/regress.h"

namespace c33 {
static std::string tmpDir() {
    const std::filesystem::path p = std::filesystem::temp_directory_path() / "ressim_c33";
    std::error_code ec; std::filesystem::create_directories(p, ec);
    return p.string();
}
static double last(const ressim::ecl::Summary& s, const std::string& key) {
    const int v = s.find(key);
    if (v < 0) return std::nan("");
    const auto x = s.series(v);
    return x.empty() ? std::nan("") : x.back();
}
// exported example deck with extra schedule text inserted before the first TSTEP / DATES
static std::string deckWith(const std::string& exampleId, const std::string& insert) {
    ressim::api::Session s;
    if (!s.open("example:" + exampleId)) return "";
    const std::string path = tmpDir() + "/" + exampleId + "_c33.DATA";
    if (!s.exportDeck(path)) return "";
    std::ifstream f(path); std::stringstream ss; ss << f.rdbuf();
    std::string t = ss.str();
    size_t sch = t.find("\nSCHEDULE");
    if (sch == std::string::npos) return "";
    size_t at = std::min(t.find("\nTSTEP", sch), t.find("\nDATES", sch));
    if (at == std::string::npos) return "";
    t.insert(at + 1, insert);
    return t;
}
}  // namespace c33

static void testC33Report() {
    std::printf("== C33 report building blocks: CRC-32, DEFLATE, PNG, DOCX / PPTX containers\n");
    using namespace ressim::report;
    const std::string nine = "123456789";
    check(crc32(reinterpret_cast<const uint8_t*>(nine.data()), nine.size()) == 0xCBF43926u, "CRC-32 check value of '123456789' = CBF43926");
    // DEFLATE round trips: repetitive text, random bytes, empty input, long runs
    std::string text;
    for (int i = 0; i < 4000; ++i) text += "WOPR:PROD" + std::to_string(i % 37) + " 1500.0 stb/d; ";
    std::mt19937 rng(3); std::string rnd(50000, '\0');
    for (char& ch : rnd) ch = char(rng() & 0xFF);
    const std::string runs(70000, 'A');
    bool rt = true; double ratioText = 1;
    for (const std::string* in : std::vector<const std::string*>{&text, &rnd, &runs}) {
        const std::string z = deflateRaw(reinterpret_cast<const uint8_t*>(in->data()), in->size());
        std::string back;
        if (!inflateRaw(z, back) || back != *in) rt = false;
        if (in == &text) ratioText = double(z.size()) / double(in->size());
    }
    { std::string back; const std::string z = deflateRaw(nullptr, 0); if (!inflateRaw(z, back) || !back.empty()) rt = false; }
    std::printf("    text compressed to %.1f %% of %zu bytes\n", 100 * ratioText, text.size());
    check(rt, "DEFLATE: compress -> inflate restores repetitive text, random bytes, long runs and empty input exactly");
    check(ratioText < 0.25, "DEFLATE: repetitive report text compresses below 25 %");
    // PNG: zlib stream (CMF 0x78), filtered scanlines, Adler-32
    {
        const int w = 16, h = 8;
        std::vector<uint8_t> rgb(size_t(w) * h * 3);
        for (int y = 0; y < h; ++y) for (int x = 0; x < w; ++x) { uint8_t* p = &rgb[(size_t(y) * w + x) * 3]; p[0] = uint8_t(16 * x); p[1] = uint8_t(32 * y); p[2] = 77; }
        const std::vector<uint8_t> png = encodePng(rgb, w, h);
        bool ok = png.size() > 57 && png[0] == 0x89 && png[1] == 'P' && png[2] == 'N' && png[3] == 'G';
        // locate IDAT
        size_t pos = 8; std::string idat;
        while (ok && pos + 8 <= png.size()) {
            const uint32_t len = (uint32_t(png[pos]) << 24) | (uint32_t(png[pos + 1]) << 16) | (uint32_t(png[pos + 2]) << 8) | png[pos + 3];
            const std::string type(reinterpret_cast<const char*>(&png[pos + 4]), 4);
            if (type == "IDAT") idat.append(reinterpret_cast<const char*>(&png[pos + 8]), len);
            pos += 12 + len;
        }
        std::string raw;
        ok = ok && idat.size() > 6 && uint8_t(idat[0]) == 0x78 && inflateRaw(idat.substr(2, idat.size() - 6), raw) && raw.size() == size_t(h) * (1 + 3 * w);
        if (ok) for (int y = 0; y < h && ok; ++y) {
            ok = raw[size_t(y) * (1 + 3 * w)] == 0;
            for (int x = 0; x < w * 3 && ok; ++x) ok = uint8_t(raw[size_t(y) * (1 + 3 * w) + 1 + x]) == rgb[size_t(y) * w * 3 + x];
        }
        if (ok) {   // Adler-32 trailer
            uint32_t a = 1, b = 0;
            for (unsigned char ch : raw) { a = (a + ch) % 65521; b = (b + a) % 65521; }
            const size_t n = idat.size();
            const uint32_t ad = (uint32_t(uint8_t(idat[n - 4])) << 24) | (uint32_t(uint8_t(idat[n - 3])) << 16) | (uint32_t(uint8_t(idat[n - 2])) << 8) | uint8_t(idat[n - 1]);
            ok = ad == ((b << 16) | a);
        }
        check(ok, "PNG: zlib-wrapped DEFLATE IDAT decodes to the filtered RGB scanlines, Adler-32 matches");
    }
    // simulation report of the quarter five-spot as DOCX and PPTX: valid OOXML containers, compact
    ressim::api::Session s;
    const bool ran = s.open("example:qfs") && s.run();
    bool okD = false, okP = false; size_t szD = 0, szP = 0; int slides = 0;
    if (ran) {
        const Doc d = simulationReport(s.caseData(), s.output(), s.summary(), Options());
        std::string err;
        const std::string pd = c33::tmpDir() + "/qfs.docx", pp = c33::tmpDir() + "/qfs.pptx";
        if (writeDocx(pd, d, err)) {
            auto z = readZip(pd, err);
            bool doc = false, ct = false, chart = false;
            for (auto& e : z) { doc |= e.name == "word/document.xml" && e.data.find("<w:body>") != std::string::npos; ct |= e.name == "[Content_Types].xml"; chart |= e.name.find("charts/chart") != std::string::npos; }
            szD = size_t(std::filesystem::file_size(pd));
            okD = doc && ct && chart && szD < 250000;
        }
        if (writePptx(pp, d, err)) {
            auto z = readZip(pp, err);
            bool pres = false;
            for (auto& e : z) { pres |= e.name == "ppt/presentation.xml"; if (e.name.rfind("ppt/slides/slide", 0) == 0 && e.name.find(".xml") != std::string::npos && e.name.find("_rels") == std::string::npos) ++slides; }
            szP = size_t(std::filesystem::file_size(pp));
            okP = pres && slides >= 6 && szP < 300000;
        }
    }
    std::printf("    docx %zu bytes, pptx %zu bytes, %d slides\n", szD, szP, slides);
    check(ran && okD, "report: DOCX with document body, content types and native charts, deflated (< 250 KB)");
    check(ran && okP, "report: PPTX with >= 6 slides, deflated (< 300 KB)");
}

static void testC33Api() {
    std::printf("== C33 scripting API: whole vs step-wise run, schedule steering, well cumulatives\n");
    using ressim::api::Session;
    Session a; const bool okA = a.open("example:qfs") && a.run();
    Session b; bool okB = b.open("example:qfs") && b.begin();
    int n = 0; while (okB && b.step() > 0) ++n;
    const double fA = okA ? c33::last(a.summary(), "FOPT") : 0, fB = okB ? c33::last(b.summary(), "FOPT") : -1;
    std::printf("    whole run FOPT %.6g, step-wise %.6g (%d report steps)\n", fA, fB, n + 1);
    check(okA && okB && std::fabs(fA - fB) <= 1e-9 * std::max(1.0, std::fabs(fA)), "Session: step-wise run (begin / step) reproduces the whole run exactly");
    const double wopt = okA ? c33::last(a.summary(), "WOPT:PROD") : 0;
    check(okA && std::fabs(wopt - fA) <= 1e-6 * fA, "summary: producer WOPT equals FOPT (single producer)");
    // steering: shut the producer half-way
    Session c; bool okC = c.open("example:qfs") && c.begin();
    const int half = okC ? c.numReportSteps() / 2 : 0;
    for (int i = 0; okC && i < half; ++i) okC = c.step() > 0;
    std::vector<std::string> msgs;
    okC = okC && c.applySchedule("WELOPEN\n 'PROD' SHUT /\n/\n", &msgs);
    while (okC && c.step() > 0) {}
    const double fC = okC ? c33::last(c.summary(), "FOPT") : 0, rC = okC ? c33::last(c.summary(), "FOPR") : -1;
    std::printf("    producer shut at step %d: FOPT %.6g (base %.6g), final FOPR %.3g\n", half, fC, fA, rC);
    check(okC && fC < 0.8 * fA && rC == 0.0, "Session::applySchedule: WELOPEN SHUT mid-run stops production from that step on");
}

static void testC33Regress() {
    std::printf("== C33 regression harness: manifest, self-reference run, deviation detection\n");
    using namespace ressim;
    const std::string dir = c33::tmpDir() + "/reg";
    std::error_code ec; std::filesystem::create_directories(dir, ec);
    api::Session s;
    bool ok = s.open("example:qfs") && s.exportDeck(dir + "/QFS.DATA");
    api::Session r; ok = ok && r.open(dir + "/QFS.DATA") && r.run() && r.writeEclipse(dir + "/REF");
    { std::ofstream m(dir + "/manifest.txt");
      m << "# C33 test manifest\ntol 0.02\nkeys FOPR,FOPT,FWCT,FPR,WBHP:*,WOPR:*\n"
        << "case QFS deck=QFS.DATA ref=REF.SMSPEC group=basic\n"
        << "skip OTHER not part of the test\n"; }
    regress::Manifest M; std::string err;
    const bool okM = ok && regress::readManifest(dir + "/manifest.txt", M, err);
    const bool parsed = okM && std::fabs(M.tol - 0.02) < 1e-12 && M.cases.size() == 2 && M.keys.size() == 6 && M.cases[0].group == "basic" && !M.cases[1].skipReason.empty();
    check(parsed, "manifest: tol / keys / case / skip directives");
    regress::Options o; o.outDir = dir + "/out";
    std::vector<regress::CaseResult> res;
    if (parsed) res = regress::runAll(M, o);
    const bool pass = res.size() == 2 && res[0].status == "PASS" && res[0].worst < 1e-4 && res[1].status == "SKIP";
    std::printf("    %s worst %.2e (%s), %zu vectors\n", res.empty() ? "-" : res[0].status.c_str(), res.empty() ? 0.0 : res[0].worst,
                res.empty() ? "" : res[0].worstKey.c_str(), res.empty() ? size_t(0) : res[0].vectors.size());
    check(pass, "runAll: a deck against its own ECLIPSE-format output passes (relRMS < 1e-4), skip entries reported");
    if (pass) {
        regress::CaseResult c = res[0];
        const int v = c.ref.find("FOPR");
        if (v >= 0) for (float& x : c.ref.values[size_t(v)]) x *= 1.1f;
        regress::compare(c, {"FOPR", "FPR"}, 0.02);
        std::printf("    FOPR reference +10 %%: %s worst %.4f (%s)\n", c.status.c_str(), c.worst, c.worstKey.c_str());
        check(c.status == "FAIL" && c.worstKey == "FOPR" && c.worst > 0.03 && c.worst < 0.1, "compare: a 10 % deviation of one vector is flagged with relRMS 3-10 %");
        const int bad = regress::writeOutputs(res, o, err);
        check(bad == 0 && std::filesystem::exists(o.outDir + "/report.html") && std::filesystem::exists(o.outDir + "/junit.xml") &&
              std::filesystem::exists(o.outDir + "/report.docx"), "writeOutputs: HTML + DOCX report, JUnit XML, history");
        // C36: cancel - before the start nothing runs; set after the first case, the rest is not run
        regress::Manifest M3 = M; M3.cases = {M.cases[0], M.cases[0], M.cases[0]};
        std::atomic<bool> stop{true};
        regress::Options oc = o; oc.cancel = &stop;
        const auto r0 = regress::runAll(M3, oc);
        stop = false;
        oc.log = [&](const std::string& l) { if (l.rfind("PASS", 0) == 0) stop = true; };
        const auto r1 = regress::runAll(M3, oc);
        std::printf("    cancel: before start %zu cases, after the first %zu case(s)\n", r0.size(), r1.size());
        check(r0.empty() && r1.size() == 1 && r1[0].status == "PASS", "C36: deck regression can be cancelled (no further cases run, finished results kept)");
    }
}

static void testC33ScheduleItems() {
    std::printf("== C33 opm-tests items: WEFAC, WLIST + patterns, wells created at run time\n");
    using ressim::api::Session;
    // WEFAC 0.5 on both wells: the reservoir sees half the flowing rates -> in-place water change = FWIT - FWPT
    const std::string deck = c33::deckWith("qfs", "WEFAC\n 'INJ' 0.5 /\n 'P*' 0.5 /\n/\n");
    Session e; const bool okE = !deck.empty() && e.openText(deck) && e.run();
    bool massOk = false, rateOk = false;
    if (okE) {
        const auto& pts = e.output().summary;
        const auto& p0 = pts.front(); const auto& p1 = pts.back();
        const double dWip = p1.wip - pts.front().wip + (p0.fwit - p0.fwpt);   // from the first recorded point
        const double net = p1.fwit - p1.fwpt;
        massOk = std::fabs(dWip - net) < 0.01 * std::max(1.0, std::fabs(net));
        const auto& w = p1.wells;
        double wwir = 0; for (auto& r : w) wwir += r.waterInj;
        rateOk = std::fabs(p1.fwir - 0.5 * wwir) < 1e-9 * std::max(1.0, wwir) && wwir > 0;
        std::printf("    in-place water change %.6g vs FWIT - FWPT %.6g; FWIR %.4g = 0.5 x WWIR %.4g\n", dWip, net, p1.fwir, wwir);
    }
    check(okE && rateOk, "WEFAC: field rates are efficiency x flowing well rates");
    check(okE && massOk, "WEFAC: the reservoir source uses efficiency x rate (water balance against FWIT - FWPT within 1 %)");
    // WLIST + pattern '*PR?D' style: a list name with '?' wildcard, applied by a schedule script
    Session l; bool okL = l.open("example:qfs") && l.begin();
    for (int i = 0; okL && i < 5; ++i) okL = l.step() > 0;
    std::vector<std::string> msgs;
    okL = okL && l.applySchedule("WLIST\n '*PLIST' NEW PROD /\n/\nWELOPEN\n '*PL?ST' SHUT /\n/\n", &msgs);
    while (okL && l.step() > 0) {}
    const double rL = okL ? c33::last(l.summary(), "FOPR") : -1;
    check(okL && rL == 0.0, "WLIST: '*PL?ST' (wildcards in a list name) expands to the list's wells; WELOPEN SHUT applied");
    // a new producer created mid-run by WELSPECS + COMPDAT + WCONPROD
    Session base; const bool okBase = base.open("example:qfs") && base.run();
    Session nw; bool okN = nw.open("example:qfs") && nw.begin();
    for (int i = 0; okN && i < 10; ++i) okN = nw.step() > 0;
    okN = okN && nw.applySchedule("WELSPECS\n 'P2' 'G' 25 1 1* OIL /\n/\nCOMPDAT\n 'P2' 25 1 1 1 OPEN 1* 1* 0.5 /\n/\n"
                                  "WCONPROD\n 'P2' OPEN ORAT 100 4* 1500 /\n/\n", &msgs);
    while (okN && nw.step() > 0) {}
    const double q2 = okN ? c33::last(nw.summary(), "WOPR:P2") : 0, fN = okN ? c33::last(nw.summary(), "FOPT") : 0,
                 fB = okBase ? c33::last(base.summary(), "FOPT") : 0;
    std::printf("    new well P2: final WOPR %.4g, FOPT %.6g vs base %.6g\n", q2, fN, fB);
    const double wopt2 = okN ? c33::last(nw.summary(), "WOPT:P2") : 0;
    bool created = false; for (auto& m : msgs) created |= m.find("new well P2 created") != std::string::npos;
    std::printf("    WOPT:P2 %.6g, message '%s'\n", wopt2, msgs.empty() ? "" : msgs.back().c_str());
    check(okN && okBase && created && q2 > 1.0 && wopt2 > 1000 && std::fabs(fN - fB) > 1e-3 * fB,
          "run-time well: WELSPECS + COMPDAT + WCONPROD from a script create and open a new producer (own summary vectors)");
}
