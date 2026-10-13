// ResSim Studio application: case editing UI, run control, 3D viewer and result plots.
#pragma once
#include "ressim/case.h"
#include "ressim/grid.h"
#include "ressim/tools.h"
#include "viewport3d.h"
#include "sim_runner.h"
#include "ressim/study.h"
#include "ressim/streamline.h"
#include "ressim/io.h"
#include "ressim/ecl.h"
#include "opm_engine.h"
#include "figure.h"
#include "deck_intel.h"
#include <atomic>
#include <functional>
#include <thread>
#include "imgui.h"
#include <memory>
#include <string>
#include <vector>

bool fileDialog(bool save, const char* filter, const char* defExt, std::string& path);
// schema-driven record table of a keyword (ui_editor.cpp); returns true when the form changed
bool deckFormRecords(ressim::deckintel::Form& F, const std::string& units);
std::string safeFileName(const std::string& name);

enum class Section { Runspec, Grid, PVT, SCAL, Init, Aquifers, EOR, Wells, Schedule, Numerics, Advanced, Count };

class App {
public:
    void init();
    void shutdown();
    void frame();
    void openFile(const std::string& path);
    bool confirmClose();
    bool quitRequested() const { return quit_; }
    bool wantsContinuousFrames() const;

    // workspaces: Model (docked case editing + run), Deck Editor ("writing mode"), Results
    enum class Workspace { Model, DeckEditor, Results, Studies };
    struct DeckEditorState;   // ui_editor.cpp

private:
    // ---- case & builders
    ressim::Case case_;
    std::string path_;
    bool dirty_ = false;
    ressim::GeomInput geom_;
    ressim::RockGenInput rock_;
    ressim::CoreyInput corey_;
    ressim::CorrelationInput corr_;
    bool autoBuild_ = true;

    // ---- model preview (grid built from the current inputs)
    std::shared_ptr<ressim::Grid> preview_;
    std::string gridError_;
    bool gridDirty_ = true;
    unsigned caseVersion_ = 1;          // bumps on every edit (for view cache keys)

    struct InitPreview {
        bool valid = false;
        unsigned caseVersion = 0;
        std::shared_ptr<ressim::Grid> grid;
        std::vector<float> p, sw, sg, so, rs;
        ressim::InPlace ip;
        std::vector<double> depth;      // per active cell (for P-vs-depth plot)
        std::string error;
    } init_;

    // ---- run
    SimRunner runner_;
    unsigned resultsCaseVersion_ = 0;   // case version that produced the current results

    // ---- 3D view
    Viewport3D view_;
    bool viewOk_ = false;
    int prop_ = 0;
    int snap_ = 0;
    bool play_ = false;
    double playAcc_ = 0;
    bool autoRange_ = true;
    float rmin_ = 0, rmax_ = 1;
    int hoverCell_ = -1;
    struct ViewKey { int prop = -1, snap = -1; unsigned version = 0, runVer = 0; const void* grid = nullptr; bool autoR = true; float a = 0, b = 0; bool log = false; };
    ViewKey key_;
    const void* viewGrid_ = nullptr;
    std::shared_ptr<const ressim::Grid> viewHold_;   // keeps the grid shown in the viewport alive

    // ---- UI state
    Section section_ = Section::Grid;
    int selWell_ = 0, selStep_ = 0;
    int resultsWell_ = 0;
    bool quit_ = false, showAbout_ = false, showImGuiDemo_ = false, resetLayout_ = false;
    std::string toast_;
    double toastUntil_ = 0;
    std::vector<char> titleBuf_;
    bool logAutoScroll_ = true;

    // ---- studies (sensitivity / history matching) - background worker
    struct StudyState {
        std::thread th;
        std::atomic<bool> running{false}, cancel{false};
        std::atomic<int> done{0}, total{0};
        std::mutex m;
        std::vector<std::string> log;
        std::string error, task;
        double t0 = 0, seconds = 0;
        // results
        bool hasTornado = false, hasMC = false, hasHM = false, hasManual = false;
        ressim::TornadoResult tornado; ressim::Response tornadoResp = ressim::Response::FOPT;
        ressim::MonteCarloResult mc; ressim::Response mcResp = ressim::Response::RF;
        ressim::HMResult hm; int hmMethod = 0;  // 0 LM, 1 ES-MDA
        std::vector<ressim::StudyParam> hmParams; // parameters used by the HM result
        ressim::RunResult manual; double manualMisfit = 0; std::vector<double> manualVals;
        // options
        int response = 0, mcSamples = 64, workers = 0, seed = 1;
        ressim::Economics eco;
        ressim::LMOptions lm;
        ressim::ESMDAOptions es;
        int synthMask = (1 << 4) | (1 << 5) | (1 << 10);   // WBHP, WWCT, FPR
        float synthNoise = 0.03f;
        bool manualAuto = true; bool manualPending = false; double manualEditTime = 0;
        int plotObs = 0;
        unsigned manualVersion = 0;
        std::thread manualTh;
        std::atomic<bool> manualRunning{false};
        std::function<void()> deferred;        // actions requested while the results lock is held
        // ---- run archive (every simulation of every study; Studies workspace)
        struct ArchRun {
            int id = 0; std::string tag; bool ok = false;
            std::vector<double> values;        // one per archParams entry
            ressim::RunResult run;
            double misfit = NAN, response = NAN;
            std::vector<double> terms;         // objective terms (termNames)
        };
        std::vector<ArchRun> archive;
        std::vector<ressim::StudyParam> archParams;
        std::vector<std::string> termNames;
        int nextId = 1, studyFirstId = 1, selRun = -1;
        unsigned archiveVersion = 0;
        int grouping = 0;                      // ObjectiveGrouping of the objective terms
        // ---- experimental design + proxy
        bool hasProxy = false, proxyQuad = true, proxyMax = true;
        int designKind = 0, designN = 24, proxyResp = 0;
        ressim::ProxyModel proxy;
        std::vector<std::vector<double>> designX;
        std::vector<double> designY, proxyOpt;
        std::vector<ressim::StudyParam> proxyParams;
        double proxyPred = NAN, proxyActual = NAN;
        // ---- workspace UI
        int wsMethod = 0, ovVec = 4, ovWell = 0, ovFilter = 0, cpX = 0, cpY = 0, parA = 0, parB = 1, centerTab = -1;
        bool ovBand = true, ovArchive = true, ovObs = true;
        int wsMode = 0;                        // 0 sensitivity & uncertainty, 1 history matching
        int sumHistory = 1; float sumSpacing = 30; unsigned sumMask = 0x1FFFu; bool sumAppend = false;   // summary import
        int mqView = 0, mqRun = 0;             // match quality: 0 series x runs, 1 wells x quantities (of run mqRun)
    } study_;
    struct MatchCache {
        unsigned key = 0; double chk = -1;
        std::vector<std::string> cols; std::vector<const ressim::RunResult*> runs;
        std::vector<std::vector<double>> rms;  // [run][series] normalised RMS
        std::vector<double> total;             // [run] misfit
    } mq_;
    void obsImportSummaryPopup();
    void hmColumns(std::vector<std::string>& names, std::vector<const ressim::RunResult*>& runs, ressim::RunResult& lastHolder);
    void hmMatchQuality();
    void hmParameters();
    struct ArchiveCtx {
        ressim::Case base;
        std::vector<double> sigma;
        std::vector<ressim::ObjectiveTerm> terms;
        ressim::Response resp = ressim::Response::FOPT;
        ressim::Economics eco;
    };
    void archiveRun(const ArchiveCtx& ctx, const char* tag, const std::vector<double>& values, const ressim::RunResult& run);
    void studiesWorkspace();
    void studiesSetup();
    void studiesCenter();
    void studiesRuns();
    void applyStudyValues(const std::vector<ressim::StudyParam>& params, const std::vector<double>& values, const std::string& what);
    void runManual();
    void startStudy(const std::string& task, std::function<void(App::StudyState&, const ressim::StudyCallbacks&)> fn);
    void panelStudies();
    void studyParams();
    void studyObservations();
    void studySensitivity();
    void studyHistoryMatch();
    void plotObservationFit(const char* id, const ressim::Case& c, const std::vector<const ressim::RunResult*>& runs,
                            const std::vector<const char*>& labels, int obsGroup, float height);
    bool showStudies_ = true;

    // ---- deck tools (import report + keyword reference), ui_deck.cpp
    void panelDeckTools();
    void keywordDocView(const std::string& kw);
    void openKeywordReference(const std::string& kw);
    bool showDeckTools_ = false, deckShowInfo_ = true, deckShowWarn_ = true, kwSupportedOnly_ = false, deckTabSet_ = false;
    int deckTab_ = 0, deckSel_ = -1;
    std::vector<ressim::DeckMessage> deckReport_;
    std::string deckReportPath_, kwSel_;
    char kwSearch_[64] = {0};

    // ---- simulation engine + persistent settings (ui_viewer.cpp)
    struct Settings {
        float uiScale = -1.0f;             // interface size (fonts + spacing), 0.6-1.6; <= 0 = automatic from the screen size
        int engine = 0;                    // 0 ResSim (built-in), 1 OPM Flow (external), 2 ECLIPSE (external, licensed)
        std::string flowPath;              // empty = auto-detect
        int opmThreads = 0;                // 0 = auto
        float opmMaxDt = 0;                // days, 0 = OPM default
        bool opmTuning = false;
        bool opmKeepWells = true;          // never shut wells that fail to converge (OPM default shuts them)
        bool opmAutoView = true;           // open OPM results in the Result Viewer
        std::string opmArgs;
        std::string eclPath;               // eclrun.exe / eclipse.exe or its folder, empty = auto-detect
        std::string eclVersion;            // empty = latest installed
        int eclProgram = 0;                // 0 E100 (eclipse), 1 E300, 2 FrontSim
        std::string eclArgs;
        int deckRunMode = 0;               // external engines, imported decks: 0 run the original deck while unedited, 1 always the generated deck
    } settings_;
    std::string srcDeck_;                  // imported ECLIPSE deck (original file) and the case version right after the import
    unsigned srcDeckVersion_ = 0;
    void applySourceDeck(OpmOptions& o);
    std::string eclDetected_;
    std::vector<std::string> eclVersions_;
    bool eclDetectDone_ = false;
    void detectEclipse(bool force);
    static const char* engineName(int e) { return e == 1 ? "OPM Flow" : e == 2 ? "ECLIPSE" : "ResSim"; }
    std::string settingsPath_, flowDetected_, flowVersion_;
    bool flowDetectDone_ = false;
    void loadSettings();
    void saveSettings();
    void detectFlow(bool force);
    void opmSettingsWindow();
    bool showOpmSettings_ = false;

    // ---- Result Viewer: external ECLIPSE / OPM Flow results and ResSim runs, overlay + comparison (ui_viewer.cpp)
    struct SLState {
        bool show = false, hideCells = false, xray = true;
        int lines = 300, colorBy = 0;          // 0 TOF, 1 injector, 2 producer
        bool focusOnly = true;                 // with a focused well: only the lines that start or end in it
        std::string builtFocus;
        int snap = -1; unsigned runVer = 0; int builtLines = -1, builtColor = -1; float builtZ = -1;
        std::shared_ptr<ressim::StreamlineResult> res;
        std::shared_ptr<const ressim::Grid> grid;
        double ms = 0;
        std::string source;                    // where the fluxes came from (shown in the streamline popup)
    } sl_;
    // Binding of one 3D viewer (Model or Results) to the shared toolbar / canvas (viewer_ui.cpp): both viewers
    // offer the same display, filter, slicing, transparency, well focus, streamline and snapshot options.
    struct ViewerBinding {
        Viewport3D* view = nullptr;
        bool* autoRange = nullptr; float* rmin = nullptr; float* rmax = nullptr;
        SLState* sl = nullptr;
        bool slAvailable = false;              // streamlines can be traced for the current data
        std::string slHint;                    // why not (shown in the popup)
        std::function<void()> slPrepare;       // optional: build what tracing needs (e.g. load the deck of external results)
        const ressim::Grid* grid = nullptr;
        std::string id, pngName, units;
        std::string propName;
        // value history of one active cell over the report steps (time in days); false if unavailable
        std::function<bool(int cell, std::vector<double>& t, std::vector<double>& v)> cellHistory;
    };
    struct ViewerUi {                          // per-viewer UI state of the shared controls
        int sliceMode = 0, slice = 0;          // 0 none, 1 I, 2 J, 3 K single slice
        int histCell = -1; std::vector<double> histT, histV; std::string histProp; bool histOpen = false;
        int ctxCell = -1, ctxWell = -1;
        const void* dualGrid = nullptr;        // dual porosity grid for which the medium filter was initialised
    };
    void viewerToolbar(ViewerBinding& b, ViewerUi& u);
    void viewerCanvas(ViewerBinding& b, ViewerUi& u, const std::function<void(int cell)>& tooltip);
    ViewerUi mvUi_, rvUi_;

    struct RVCase {
        std::string name, path;
        bool show = true, internal = false;
        ImU32 color = 0;
        std::shared_ptr<ressim::ecl::ResultCase> rc;    // external files (null for a ResSim run)
        ressim::ecl::Summary summary;                   // ResSim run (external cases use rc->summary)
        std::shared_ptr<const ressim::Grid> grid;       // ResSim run grid
        std::vector<ressim::Snapshot> snaps;            // ResSim run report steps
        std::shared_ptr<const ressim::Case> input;      // model input (ResSim / in-app runs, or the deck found next to the results)
        std::shared_ptr<ressim::Simulator> sim;         // for streamlines (ResSim run: the run's simulator; external: built from the deck)
        bool deckTried = false;
        std::string simNote;                            // why streamlines are unavailable
        const ressim::ecl::Summary& sum() const { return rc ? rc->summary : summary; }
        const ressim::Grid* grid3d() const { return rc ? rc->grid.get() : grid.get(); }
    };
    struct RVState {
        std::vector<RVCase> cases;
        std::vector<std::string> selected;     // vector keys
        char filter[64] = {0};
        int category = 0;
        bool years = false, perVector = true, onlyCommon = false;
        int ref = 0, cmp = 1, cmpRef = -1, cmpCmp = -1;
        bool cmpFieldOnly = false, cmpWells = true;
        std::vector<ressim::ecl::VectorDiff> diffs;
        int tab = 0, wantTab = -1;
        // 3D
        Viewport3D view; bool viewInit = false, viewOk = false;
        int c3d = 0, step = 0; std::string prop = "PRESSURE";
        std::string key3d; const void* grid3d = nullptr;
        bool autoRange = true; float rmin = 0, rmax = 1;
        std::vector<float> vals;
        std::string caption;
        int colorIdx = 0;
        int colorMode = 0, gridCols = 0;
        SLState sl;
        bool wellsSet = false;
    } rv_;
    bool rvEnsureInput(RVCase& c, bool wantSim);
    void rvUpdateStreamlines(RVCase& c, int step);
    bool showResultViewer_ = false;
    std::function<void()> rvDeferred_;       // Result Viewer actions that must run outside of the panel loop / results lock
    void panelResultViewer();
    void rvOpen(const std::string& path, bool quiet = false);
    void rvAddRun();
    void rvPlots();
    void rvCompare();
    void rv3D();
    void rvInfo();

    // ---- streamlines
    void updateStreamlines();
    void editAquifers();
    void editEOR();
    int eorPvtSel_ = 0, eorSatSel_ = 0, eorStepSel_ = 0;
    void importDeck(const std::string& path);

    // ---- workspaces (ui_workspace.cpp, ui_editor.cpp)
    Workspace workspace_ = Workspace::Model;
    std::shared_ptr<DeckEditorState> ed_;
    void edEnsure();
    void edOpenFile(const std::string& path);
    void edOpenModelDeck();
    bool edSave(int idx, bool saveAs);
    bool edConfirmClose();
    void deckEditorWorkspace();
    void resultsWorkspace();
    void workspaceSwitch();
    int resultsTab_ = 0, resultsWantTab_ = -1;
    void plotsBody();          // run dashboard (Field / Wells / Custom / Convergence)
    void rvBody();             // Result Viewer contents (cases, vectors, plots, compare, 3D, info)
    // publication figure export (ui_workspace.cpp)
    struct FigExport {
        bool open = false;
        std::string id;
        fig::Figure fig;
        fig::Style style;
        int preset = 1, cols = 0;
        std::string lastDir;
    } figx_;
    void exportButton(const char* figId);
    void figureExportWindow();

    // ---- frame parts
    void dockSpace();
    void menuBar();
    void toolbar();
    void statusBar();
    void panelCase();
    void panelEditor();
    void panelViewer();
    void panelPlots();
    void panelLog();
    void panelRun();
    void aboutWindow();

    // ---- section editors
    void editRunspec();
    void editGrid();
    void editPVT();
    void editSCAL();
    void editEndPointScaling();
    void editGroups();
    std::string selGroup_;
    void editInit();
    void editWells();
    void editSchedule();
    void editNumerics();
    void editAdvanced();          // extra (verbatim) keywords for OPM Flow: groups, VFP, MSW, ACTIONX, summary requests, ...
    int advSel_ = -1;
    std::string advFilter_, advNewName_;
    int advNewSection_ = 7, advCheckSel_ = -2;
    std::string advCheckBody_, advCheckName_;
    std::vector<std::string> advCheckMsgs_;
    int advBigEdit_ = -1;
    bool advFormView_ = true;                  // Advanced keywords: schema form (true) or text
    int advFormKey_ = -1;
    std::string advFormSrc_;
    ressim::deckintel::Form advForm_;                      // index of a very large keyword the user chose to edit anyway

    // ---- actions
    void newCase(int which, const std::string& example = std::string());   // 3 = built-in example
    void examplesMenu(bool popup);
    void save(bool saveAs);
    void exportDeck();
    void runSimulation();
    void computeInit();
    void rebuildPreview();
    void markDirty(bool gridChanged = false);
    void markStudyDirty() { dirty_ = true; }   // study params / observed data: project modified, results still valid
    void syncBuildersFromCase();
    void notify(const std::string& msg);
    bool maybeSaveChanges();                      // false = user cancelled
    void openDialog();
    bool wasRunning_ = false;
    std::string sectionStatus(Section s) const;   // empty = OK, otherwise problem
    const std::string& validation() const;          // cached validateCase() message
    mutable unsigned validatedVersion_ = 0;
    mutable std::string validation_;

    // ---- viewer helpers
    void updateViewData();
    std::vector<WellGlyph> wellGlyphs(const ressim::Grid& g, const ressim::Case& c, int step = -1) const;
};
