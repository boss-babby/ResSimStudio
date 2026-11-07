// ResSim Studio application: case editing UI, run control, 3D viewer and result plots.
#pragma once
#include "ressim/workflow.h"
#include "ressim/case.h"
#include "ressim/grid.h"
#include "ressim/tools.h"
#include "viewport3d.h"
#include "sim_runner.h"
#include "custom_plots.h"
#include "ressim/study.h"
#include "ressim/streamline.h"
#include "ressim/io.h"
#include "ressim/ecl.h"
#include "opm_engine.h"
#include "figure.h"
#include "deck_intel.h"
#include <atomic>
#include <mutex>
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

int propIndexByName(const char* name);   // index of a Model-viewer property (ui_results.cpp), -1 = unknown

class App {
public:
    struct GeoState;   // C32 Geomodel workspace state (ui_geomodel.cpp)
    void init();
    void shutdown();
    void frame();
    void openFile(const std::string& path);
    bool confirmClose();
    bool quitRequested() const { return quit_; }
    bool wantsContinuousFrames() const;

    // workspaces: Model (docked case editing + run), Deck Editor ("writing mode"), Results
    enum class Workspace { Model, DeckEditor, Results, Studies, Geomodel };
    void setGeoMode(bool on);   // C32 ResGeo: standalone geomodelling front end (--geo / ResGeo.exe)
    bool smokeGeoDemo_ = false;
    std::string smokeAction_;   // C40 smoke tests: wizard | palette | run | objects (scripted UI state before the screenshot)
    int smokeFrame_ = 0;
    int geoStepReq_ = -1, geoViewReq_ = -1;   // C40 palette deep links into the Geomodel workspace   // C39 smoke test: load the Geomodel demo and open its 3-D scene
    void setWorkspaceIndex(int w) { if (w >= 0 && w <= 4) { workspace_ = Workspace(w); showStart_ = false; } }   // C39 smoke tests
    struct DeckEditorState;   // ui_editor.cpp

private:
    // ---- case & builders
    ressim::Case case_;
    std::string path_;
    bool dirty_ = false;
    ressim::GeomInput geom_;
    ressim::RockGenInput rock_;
    ressim::CoreyInput corey_;
    ressim::VanGenuchtenInput vg_;          // C27 van Genuchten - Mualem generator (TOUGH style)
    int vgTarget_ = 0;                      // 0 SWOF, 1 SGOF
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
    struct ViewKey { int prop = -1, snap = -1; unsigned version = 0, runVer = 0; const void* grid = nullptr; bool autoR = true; float a = 0, b = 0; bool log = false; bool allSteps = true; };
    // C18.1: colour range over ALL report steps of the run (cached per run version + property) - the colours then mean the
    // same value at every time step; 'This step' restores the old per-step auto range
    std::string grKey_; float grLo_ = 0, grHi_ = 1; bool grOk_ = false;
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
        ressim::HMResult hm; int hmMethod = 0;  // 0 LM / gradient family, 1 ES-MDA
        int hmAlgo = 0;                          // gradient family: 0 LM (FD Jacobian), 1 LM (analytic), 2 adjoint L-BFGS
        int hmOpt = 0;                           // selected gradient-family optimiser (same coding as hmAlgo)
        // C13: dedicated adjoint HM method (wsMethod 6): 0 scalar study parameters (L-BFGS), 1 per-cell field;
        // result views (centre tab 11): which result is shown, selected iteration / sub-view / field map settings
        int adjMode = 0, adjLast = 0, adjIter = -1, adjFieldQty = 0, adjLayer = 0, adjResp = 0;
        bool adjAllResp = false;
        ressim::AdjointHMOptions adj;
        ressim::FieldHMOptions fieldOpt; ressim::FieldHMResult field; bool hasField = false;   // C11 per-cell adjoint HM
        ressim::GradientCheck gradCheck; bool hasGradCheck = false; int gradMask = 7; double gradFdStep = 1e-3;
        std::vector<ressim::StudyParam> gradParams;
        std::vector<ressim::StudyParam> hmParams; // parameters used by the HM result
        ressim::HMStatistics hmStats; bool hasStats = false; float statsFd = 0.02f;   // C27 iTOUGH2-style inversion statistics
        unsigned statsPredMask = 1u | (1u << 7);   // C28 FOSM predictions (bit = ressim::Response; default FOPT + NPV)
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
        int proxyKind = 0;   // C28: 0 polynomial response surface, 1 Gaussian process (kriging / ML proxy)
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
    // ---- C32 Geomodel workspace (ui_geomodel.cpp): data -> seismic -> structure -> grid -> facies -> petrophysics ->
    // volumetrics -> upscaling -> simulation case
    std::shared_ptr<GeoState> geo_;
    bool geoMode_ = false;
    void geomodelWorkspace();
    void geomodel3D(GeoState& S, ImVec2 size);          // C36 3-D scene (ui_geo3d.cpp)
    void geoAnalysis(GeoState& S, ImVec2 size);         // C40 data analysis: log calculator, cross-plot, histogram, variogram fit, well tie (ui_geoanalysis.cpp)
    void geoSeismicSection(GeoState& S, ImVec2 size);   // C36 inline / crossline section, attributes, flattening, picking
    void studiesSetup();
    void gradientCheckPanel(int activeParams);
    void componentsTab(const RunResults& R, double xs, double tSnap, const char* xl, const std::function<void()>& cursor);
    void tracersTab(const RunResults& R, double xs, double tSnap, const char* xl, const std::function<void()>& cursor);   // C27
    void prodAnalysisTab(const RunResults& R);   // C28 decline curves + nodal analysis (ui_prod.cpp)
    // C13 (ui_egs.cpp): geothermal / EDFM
    void geothermalTab(const RunResults& R, double xs, double tSnap, const char* xl, const std::function<void()>& cursor);
    // C14 user-defined plot library (ui_customplots.cpp)
    void customPlotsTab(const RunResults& R, double xs, double tSnap, const char* xl, const std::function<void()>& cursor);
    void drawCustomPlot(const CustomPlot& P, const RunResults& R, double xs, const char* xl, const std::function<void()>& cursor, ImVec2 size, int id);
    void loadCustomPlots();
    void saveCustomPlots();
    std::vector<CustomPlot> cplots_;
    int cplotSel_ = 0;
    std::string cplotsPath_;
    bool egsEditor();                                // Fractures & EGS tab of the thermal editor; false = the case was replaced
    std::vector<PanelQuad> fracturePanels(const ressim::Case& c, const ressim::Grid* g, const RunResults* R, int snap, const std::string& prop);
    void updateFracturePanels();                     // Model viewer (ui_results.cpp)
    // C24 (ui_faultmech.cpp): fault mechanics / induced seismicity of the main model
    void editFaultMechanics();
    void editReservoirGeomech();                          // C25: case.geoSeq (black oil / compositional), FE stress options
    bool editGeomechInput(ressim::geomech::Input& G, bool thermal); // shared FE geomechanics property sheet
    std::vector<PanelQuad> faultPanels(const ressim::Case& c, const ressim::Grid* g, const RunResults* R, int snap);
    void faultsTab(const RunResults& R, double xs, double tSnap, const char* xl, const std::function<void()>& cursor);
    // C30 fields / sections / fracture faces / phase field / post-processed geomechanics (ui_fieldviz.cpp)
    void fieldsTab(const RunResults& R);
    void udqTab(const RunResults& R, double xs, const char* xl);
    void networkTab(const RunResults& R, double xs, const char* xl);   // C35-4b surface network node pressures (ui_msw.cpp)
    void mswTab(const RunResults& R, double xs, const char* xl);   // C32 multisegment wells (ui_msw.cpp)
    void mswGenerator(int wellIndex);                             // C32 Model -> Wells: WELSEGS / COMPSEGS generator   // C31 UDQ vectors + ACTIONX log (ui_udq.cpp)
    void postGeoControls(int target, std::shared_ptr<const ressim::Grid> g, bool fieldUnits, const std::vector<ressim::Snapshot>* snaps);
    void pollPostGeo();
    struct PostGeoJob;
    std::shared_ptr<PostGeoJob> postGeo_;
    ressim::geomech::Input postGeoIn_;
    bool postGeoInit_ = false;
    std::string postGeoMsg_;
    int faultProp_ = 1;                              // 3D fault colour quantity (kFaultProps)
    bool faultEvents_ = true;
    float faultEventScale_ = 1.0f;
    double faultRange_[2] = {0, 1};
    std::string fracKeyMv_, fracKeyRv_;
    // C13: adjoint HM maps for the Model 3D viewer (per global cell; empty = not listed)
    std::vector<float> hmMultK_, hmMultPhi_, hmSensK_;
    unsigned hmViewVer_ = 0;
    void adjointHMView();                            // Studies centre tab 11 (ui_studies.cpp)
    const char* hmName(bool shortName) const;
    void studiesCenter();
    void studiesRuns();
    void applyStudyValues(const std::vector<ressim::StudyParam>& params, const std::vector<double>& values, const std::string& what);
    void runManual();
    void startStudy(const std::string& task, std::function<void(App::StudyState&, const ressim::StudyCallbacks&)> fn);
    void hmStatsPanel();                        // C27: linearised inversion statistics of the HM result (hmParameters view)
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
        std::string pythonExe;             // C37 PYACTION interpreter (empty: RESSIM_PYTHON, python / python3 on PATH)
        std::string eclPath;               // eclrun.exe / eclipse.exe or its folder, empty = auto-detect
        std::string eclVersion;            // empty = latest installed
        int eclProgram = 0;                // 0 E100 (eclipse), 1 E300, 2 FrontSim
        std::string eclArgs;
        int deckRunMode = 0;               // external engines, imported decks: 0 run the original deck while unedited, 1 always the generated deck
        // C14 results / restart output of the built-in engine
        int rstEvery = 5;                  // save the engine restart state every N report steps (0 = never)
        bool rstLast = true;               // ... and at the end of every run
        bool autoSaveResults = false;      // write <project>.rsres after every finished run
        bool eclOutput = false;            // also write ECLIPSE result files (EGRID / INIT / UNRST / SMSPEC / UNSMRY)
        // C16 high-performance computing
        bool mpiEnabled = false;           // studies / ensembles: distribute the runs with mpiexec + ResSimCLI.exe
        std::string mpiExec;               // mpiexec.exe, empty = Microsoft MPI default location / PATH
        bool mpiDD = false;                // C18: single black-oil runs domain-decomposed (mpiexec + ResSimCLI --mpi-dd)
        int mpiDDProcs = 2;
        bool mpiDistAsm = true;            // C28: --mpi-da (every rank assembles its own rows) instead of --mpi-dd (rank 0 assembles)
        bool rangeAllSteps = true;   // C18.1: 3D colour range over all report steps (false: per step)
                        // ranks of a domain-decomposed run
        int mpiProcs = 0;                  // MPI ranks (0 = cores + 1: rank 0 only dispatches)
        std::string mpiArgs;               // extra mpiexec arguments, e.g. -hosts 2 node1 8 node2 8
        std::string mpiJobDir;             // shared job folder (UNC path for several machines), empty = %TEMP%
        int gpuPrec = 1;                   // C18 device preconditioner: 1 multicolour BILU0, 0 block Jacobi
        bool gpuEnabled = false;           // register the sokol_gfx compute linear solver (cases with Numerics -> GPU)
        // C26 usability
        std::vector<std::string> recent;   // recently opened / saved files (projects, decks, results), newest first
        bool showStartPage = true;         // start page at start-up
        bool runHistory = true;            // keep every finished ResSim run in the Result Viewer compare list
        int runHistoryMax = 6;             // ... at most this many automatic entries
        int theme = 0;                     // C39: 0 dark, 1 light, 2 high contrast
        int accent = 0;                    // C39: 0 blue, 1 teal, 2 violet, 3 orange, 4 green
        int helpStyle = 0;                 // C42: 0 compact (i) icons with hover / click details, 1 inline help text
        bool autosave = true;              // C39: recovery copy of the unsaved project
        int autosaveSec = 120;             // ... every N seconds while there are unsaved changes
        int bg3d = 0;                      // C40: reservoir 3-D background (theme.h set3dBackground)
        bool timeCursor = true;            // C40: report-time cursor on every time plot (drag = change report step)
        std::string language = "en";       // C41 (G13): interface language (lang/<code>.txt catalogue)
    } settings_;
    // ---- C26 usability (ui_usability.cpp)
    bool paletteOpen_ = false, paletteFocus_ = false; char paletteBuf_[128] = {0}; int paletteSel_ = 0;
    bool showStart_ = false; char startFilter_[64] = {0};
    struct Note { long long when = 0; std::string text; bool error = false; };
    std::vector<Note> notes_; int notesUnread_ = 0; bool notesPopupReq_ = false;
    struct DeckFile { std::string path, name, folder, desc; };
    std::vector<DeckFile> deckFiles_; bool deckFilesScanned_ = false;
    int runCounter_ = 0; std::string lastRunText_;
    void commandPalette();
    struct PaletteCmd { std::string label, hint; std::function<void()> fn; bool enabled = true; };
    std::vector<PaletteCmd> buildCommands(bool withKeywords);
    void commandPaletteWindow(std::vector<PaletteCmd>& C);
    // ---- C41 workflow recorder (ui_workflow.cpp, engine ressim/workflow.h)
    bool showWorkflow_ = false, wfRecording_ = false;
    std::vector<ressim::wf::Step> wfSteps_;
    ressim::wf::Runner wfRun_;
    std::string wfPath_;
    void wfRecord(const std::string& label);
    void wfTick();
    void panelWorkflow();
    void startPage();
    void workflowChecklist();
    void notificationButton();
    void pushNote(const std::string& msg);
    void addRecent(const std::string& path);
    void recentMenu();
    void scanDeckFiles();
    void deckFilesMenu();
    void rvAddRunHistory();
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
        bool* rangeAllSteps = nullptr;          // C18.1: auto range over all time steps (nullptr: static data, no choice)
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
        bool autoHist = false;                          // C26: automatic run-history entry
        std::string changes;                            // C26: input changes against the previous run
        ImU32 color = 0;
        std::shared_ptr<ressim::ecl::ResultCase> rc;    // external files (null for a ResSim run)
        ressim::ecl::Summary summary;                   // ResSim run (external cases use rc->summary)
        std::vector<ressim::Snapshot> post;             // C30: post-processed geomechanics arrays per restart step (external cases)
        std::shared_ptr<const ressim::Grid> grid;       // ResSim run grid
        std::vector<ressim::Snapshot> snaps;            // ResSim run report steps
        std::shared_ptr<const ressim::Case> input;      // model input (ResSim / in-app runs, or the deck found next to the results)
        std::shared_ptr<ressim::Simulator> sim;         // for streamlines (ResSim run: the run's simulator; external: built from the deck)
        bool deckTried = false;
        std::string simNote;                            // why streamlines are unavailable
        RunResults egs;                                 // C13: EDFM segments / states of a ResSim thermal run (fracture fields only)
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
        std::string grKey; float grLo = 0, grHi = 1; bool grOk = false;   // C18.1 all-steps range cache
        std::vector<float> vals;
        std::string caption;
        int colorIdx = 0;
        int colorMode = 0, gridCols = 0;
        SLState sl;
        bool wellsSet = false;
    } rv_;
    std::vector<float> rvLoad(RVCase& c, int step, const std::string& name);   // C30 restart / post-processed array
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
    void editChemEOR();   // C35 surfactant / alkaline / foam (ui_chem.cpp)
    void editGasLift();   // C35 LIFTOPT / WLIFTOPT / GLIFTOPT (ui_gaslift.cpp)
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
        // C22 exact ("as on screen") export
        ImGuiID plotId = 0;          // captured single plot (re-capture)
        bool snapWait = false;       // snapshot of a recorded figure requested, waiting for the next frame
        bool snapTried = false;      // a request came back (empty -> the plots were not visible)
        int pxScale = 2;             // PNG pixel scale relative to the screen (0 = use dpi)
    } figx_;
    void exportButton(const char* figId);
    bool figDeferredFocus_ = false;   // focus the export window after a single-plot capture
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
    void editLgrs();
    // ---- C15 (ui_hfrac.cpp): hydraulic fractures / unconventional rock, geomechanical stability window, user guide
    void editHydraulicFractures();
    void stabilityWindow();
    void openUserGuide();
    bool showStability_ = false;
    // ---- C20 (ui_seis.cpp): rate-and-state seismicity + convex minimum-seismicity injection designer
    void seismicityWindow();
    bool showSeis_ = false, seisDock_ = false;
    void efemWindow();                          // C20 G1: 2-D embedded fractures with Coulomb contact (ui_efem.cpp)
    bool showEfem_ = false, efemDock_ = false;
    bool stabDock_ = false;       // dock the stability window as a tab next to the 3D viewer when it is opened
    void editPVT();
    bool editFluidModel();        // model selector; true = black oil (draw the PVT tables)
    void editCompositional();     // EOS components, BIC, ZI / ZMFVD, streams, separator, phase envelope, flash
    void editThermal();           // E300 THERMAL: components / K-values, viscosity tables, rock heat, heat loss, TEMPI / TEMPVD, WINJTEMP
    void editCo2Store();          // CO2STORE: temperature, salinity, brine phase, contact, property plots
    void editTrapping(bool co2);  // trapping mechanisms + interactions (hysteresis, solubility / diffusion, salt, minerals, structural)
    void editHysteresis(bool compCore);   // SATOPTS HYSTER / EHYSTR / IMBNUM (+ Land Sgt,max for the compositional core), curve preview
    void editSCAL();
    void editEndPointScaling();
    void editGroups();
    std::string selGroup_;
    void editInit();
    void editWells();
    void editSchedule();
    void editRunLogic();     // C31 UDQ / ACTIONX / UDA editor (ui_logic.cpp)
    void editProductionData();   // C38 production-data import, QC, WCONHIST / WCONINJH / WEFAC generation (ui_proddata.cpp)
    void editPyActions();    // C37 PYACTION editor: records, script editor, templates, syntax check (ui_pyaction.cpp)
    void editNumerics();
    void editAdvanced();          // extra (verbatim) keywords for OPM Flow: groups, VFP, MSW, ACTIONX, summary requests, ...
    int advSel_ = -1;
    std::string advFilter_, advNewName_;
    std::string advTracerName_; int advTracerPhase_ = 0;   // C27 tracer quick-add
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
    std::unique_ptr<ressim::Case> pendingCase_;   // newCase(4): load this case (built by a designer)
    void examplesMenu(bool popup);
    void save(bool saveAs);
    void exportDeck();
    // ---- C14: results files (.rsres), ECLIPSE result export, restart from a report step (ui_restart.cpp)
    void openResultsFile(const std::string& path);
    void saveResultsDialog();
    void exportResultsEclipseDialog();
    void afterRunOutputs();               // auto-save .rsres / ECLIPSE files after a finished ResSim run
    void restartWindow();
    // ---- C33 report generator + deck regression runner (ui_report.cpp)
    void reportWindow();
    void regressionWindow();
    bool showReport_ = false, showRegress_ = false;
    struct ReportUi {
        int format = 0; char title[160] = ""; char author[96] = ""; char notes[2048] = ""; char ref[512] = ""; char refLabel[64] = "OPM Flow";
        bool inputs = true, wellCharts = true, maps = true, log = true, openAfter = true; int maxWells = 8, layer = 0;
        std::string last;
    } rep_;
    struct RegressUi { char manifest[512] = ""; char outDir[512] = ""; char filter[96] = ""; double tolPct = 0, maxSec = 600; } reg_;
    struct RegressJob;
    std::shared_ptr<RegressJob> regJob_;
    void restartSettingsUi(bool& changed);
    // ---- C16: MPI cluster batches and GPU solver (ui_hpc.cpp)
    void hpcSettingsUi(bool& changed);
    void applyHpcSettings();
    bool startMpiDD();                     // C18: launch the current black-oil case domain-decomposed (false = not applicable)
    std::string hpcMpiReport_, hpcGpuReport_;
    std::shared_ptr<std::atomic<int>> hpcBusy_ = std::make_shared<std::atomic<int>>(0);
    std::shared_ptr<std::mutex> hpcMx_ = std::make_shared<std::mutex>();
    bool showRestart_ = false;
    int restartSel_ = -1;
    std::shared_ptr<RunResults> restartBase_;   // base results loaded from a file for a restart (else the current run)
    std::string restartBaseFile_;
    void runSimulation();
    void computeInit();
    void rebuildPreview();
    void markDirty(bool gridChanged = false);
    void markStudyDirty() { dirty_ = true; }   // study params / observed data: project modified, results still valid
    void syncBuildersFromCase();
    void objectTree();                 // C40 Case Explorer object tree (ui_objects.cpp)
    void newModelWizard();             // C40 new-model wizard (ui_wizard.cpp)
    bool showWizard_ = false;
    // ---- C40 (G10) case manager (ui_cases.cpp)
    struct CaseVariant {
        std::string name, note, run;               // run: name of the linked entry of the run history (rv_.cases)
        std::shared_ptr<const ressim::Case> c;
        std::string text;                          // projectText(*c), for diffs and run linking
        long long created = 0;
        bool base = false;
        const void* diffBase = nullptr; const void* diffSelf = nullptr; int nDiff = 0;   // cached line count of the diff vs the base
    };
    std::vector<CaseVariant> variants_;
    int casesSel_ = -1;
    bool showCases_ = false;
    std::string inspFilter_;                  // C40 (G11) Inspector filter text
    void panelCases();
    void casesAdd(const std::string& name, const ressim::Case& c);
    void casesLoad(int i, bool run);
    void casesLinkRun(const std::string& inputText, std::string& runName);
    void casesWrite(const std::string& project);
    void casesRead(const std::string& project);
    int casesFindRun(const std::string& runName) const;
    char objFilter_[64] = {0};

    // ---- C39 undo / redo of model edits (snapshots of the whole case, one per edit gesture) and autosave / recovery
    struct UndoState {
        std::vector<std::shared_ptr<const ressim::Case>> stack;
        int cur = -1;
        unsigned version = 0;              // caseVersion_ of stack[cur]
        size_t maxLevels = 60;
    };
    UndoState undo_;
    void undoReset();                       // new / opened / imported case: history restarts here
    void undoTrack();                       // end of frame: snapshot after a finished edit gesture
    void undoApply(int dir);                // -1 undo, +1 redo
    bool canUndo() const { return undo_.cur > 0; }
    bool canRedo() const { return undo_.cur >= 0 && undo_.cur + 1 < int(undo_.stack.size()); }
    double autosaveAt_ = 0;
    unsigned autosavedVersion_ = 0;
    std::atomic<bool> autosaveBusy_{false};
    std::thread autosaveThr_;
    bool recoveryPrompt_ = false;
    std::string recoveryInfo_, recoveryOrig_;
    void autosaveTick();
    void recoveryCheck();                   // start-up: a recovery copy left by a crash?
    void recoveryClear();                   // saved / closed normally
    void recoveryWindow();
    std::string recoveryPath() const;
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
