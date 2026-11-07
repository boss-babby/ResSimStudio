// C20 geomechanics G1: Tools > Geomechanics: fracture contact & slip (2-D embedded fractures, Coulomb friction).
// Plane-strain section with far-field stresses and any number of straight fractures (pressure, friction, cohesion);
// results: opening / slip / contact tractions along each fracture, stress maps with the fractures drawn on top.
#include "app.h"
#include "widgets.h"
#include "icons.h"
#include "theme.h"
#include "implot.h"
#include "plotkit.h"
#include "imgui_internal.h"
#include "ressim/efem.h"
#include "ressim/geocouple.h"
#include "ressim/rsf.h"
#include "ressim/phasefield.h"
#include "ressim/simulator.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <memory>
#include <mutex>
#include <thread>

using namespace ressim;

namespace {
struct EfemState {
    efem::Input in;
    struct Row { std::string name = "F1"; double x0 = 4, y0 = 4.5, x1 = 6, y1 = 5.5, pMPa = 0, mu = 0.6, cMPa = 0; };
    std::vector<Row> rows;
    double sxxMPa = 0, syyMPa = -10, eGPa = 10;
    int field = 1;
    std::shared_ptr<std::mutex> mx = std::make_shared<std::mutex>();
    std::shared_ptr<std::atomic<bool>> busy = std::make_shared<std::atomic<bool>>(false);
    std::shared_ptr<efem::Result> res;
    std::shared_ptr<efem::Input> resIn;
    bool init = false;
    // C21 G1-c: coupled flow run of the open case
    struct Coupled {
        geocouple::Input in;
        double fricC = 0.6, cohMPa = 0, aResMm = 0.02, vmMm = 0.3, kn0GPaPerM = 20, dilDeg = 0, aRefMm = 0;
        bool compare = true;
        std::shared_ptr<std::atomic<bool>> busy = std::make_shared<std::atomic<bool>>(false), cancel = std::make_shared<std::atomic<bool>>(false);
        std::shared_ptr<std::string> status = std::make_shared<std::string>();
        struct Res { geocouple::Output out; std::vector<double> t, bhpC, tU, bhpU; std::string well; double seconds = 0; bool field = false; };
        std::shared_ptr<Res> res;
        bool select = false;
    } cp;
    // C23 G3: rate-and-state friction on one embedded fault
    struct Rsf {
        rsf::Input in;
        double sigMPa = 50, p0MPa = 0, dcMm = 1, vplMmYr = 31.6, years = 60, lengthM = 2000, gGPa = 30;
        int stress = 0;          // 0 uniform, 1 from the contact solve (fracture 'frac')
        int frac = 0;
        int pressure = 0;        // 0 none, 1 injection (line source), 2 last coupled flow run (G1-c)
        double injAt = 0.5, injOff = 50, injPMPa = 1.0, injDiff = 0.05, injStopD = 200;
        std::shared_ptr<std::atomic<bool>> busy = std::make_shared<std::atomic<bool>>(false), cancel = std::make_shared<std::atomic<bool>>(false);
        std::shared_ptr<rsf::Output> res;
        bool select = false;
    } rs;
    // C23 G4: phase-field pressurised crack (near-wellbore / lab scale)
    struct Pf {
        phasefield::Input in;
        double eGPa = 10, aCm = 10, hMm = 10;
        double sxxM = 0, syyM = 0, sxyM = 0, vMaxCm2 = 0;   // C26: in-situ stress (MPa, compression +), injected volume (cm2 per m, 0 = auto)
        int show = 1;            // phase field shown: 0 final, 1 at the critical pressure
        std::shared_ptr<std::atomic<bool>> busy = std::make_shared<std::atomic<bool>>(false), cancel = std::make_shared<std::atomic<bool>>(false);
        std::shared_ptr<std::string> status = std::make_shared<std::string>();
        std::shared_ptr<phasefield::Output> res;
        bool select = false;
    } pf;
};
EfemState g_ef;
}  // namespace

void App::efemWindow() {
    if (!showEfem_) return;
    EfemState& S = g_ef;
    if (!S.init) {
        S.init = true;
        S.in.nx = S.in.ny = 100;
        S.rows.push_back({});
        S.rows[0].x0 = 3.8; S.rows[0].x1 = 6.2; S.rows[0].y0 = 5.037 - 1.2 * std::tan(30 * M_PI / 180); S.rows[0].y1 = 5.037 + 1.2 * std::tan(30 * M_PI / 180); S.rows[0].mu = 0.3;
    }
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetFontSize() * 70, ImGui::GetFontSize() * 44), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_FirstUseEver, ImVec2(0.5f, 0.5f));
    if (efemDock_) {
        efemDock_ = false;
        if (workspace_ == Workspace::Model)
            if (ImGuiDockNode* cn = ImGui::DockBuilderGetCentralNode(ImGui::GetID("MainDock"))) ImGui::SetNextWindowDockID(cn->ID, ImGuiCond_Always);
        ImGui::SetNextWindowFocus();
    }
    if (!ImGui::Begin(ICON_ACTIVITY "  Geomechanics: fracture contact & slip (2-D)###efemwin", &showEfem_)) { ImGui::End(); return; }
    const bool busy = S.busy->load();
    if (ImGui::BeginTable("eflay", 2, ImGuiTableFlags_Resizable | ImGuiTableFlags_BordersInnerV)) {
        ImGui::TableSetupColumn("in", ImGuiTableColumnFlags_WidthStretch, 0.36f);
        ImGui::TableSetupColumn("out", ImGuiTableColumnFlags_WidthStretch, 0.64f);
        ImGui::TableNextColumn();
        ImGui::BeginChild("efin");
        ImGui::BeginDisabled(busy);
        ui::header("Section and rock");
        ui::hint("Plane-strain section, rollers on the left and bottom, far-field normal stresses on the right and top (tension +, "
                 "compression -). Fractures are embedded in the mesh (no meshing along them): continuous displacement jumps by "
                 "nodal Heaviside enrichment, penalty contact and Coulomb friction with return mapping.");
        if (ui::beginProps("efp", 0.5f)) {
            ui::propD("Width (m)", &S.in.Lx, "%.4g"); ui::propD("Height (m)", &S.in.Ly, "%.4g");
            ui::propI("Elements x", &S.in.nx, 4, 400); ui::propI("Elements y", &S.in.ny, 4, 400);
            ui::propD("Young's modulus (GPa)", &S.eGPa, "%.4g"); ui::propD("Poisson's ratio", &S.in.nu, "%.3f");
            ui::propD("Far-field sigma_xx (MPa)", &S.sxxMPa, "%.4g"); ui::propD("Far-field sigma_yy (MPa)", &S.syyMPa, "%.4g");
            static const char* jm[] = {"Continuous jumps (accurate)", "Element-constant jumps (classic EFEM)"};
            ui::propCombo("Jump interpolation", &S.in.jumpModel, jm, 2, "Element-constant jumps over-predict slip by 10-25 % on meshes oblique to the fracture (research finding F9)");
            ui::endProps();
        }
        ui::header("Fractures");
        ui::hint("End points on element edges give the most accurate tips. Fractures must be at least one element apart.");
        {
            ui::RecordTable<EfemState::Row> t;
            t.text("Name", [](EfemState::Row& r) -> std::string& { return r.name; })
             .num("x0 (m)", [](EfemState::Row& r) -> double& { return r.x0; }).num("y0 (m)", [](EfemState::Row& r) -> double& { return r.y0; })
             .num("x1 (m)", [](EfemState::Row& r) -> double& { return r.x1; }).num("y1 (m)", [](EfemState::Row& r) -> double& { return r.y1; })
             .num("Pressure (MPa)", [](EfemState::Row& r) -> double& { return r.pMPa; })
             .num("Friction", [](EfemState::Row& r) -> double& { return r.mu; })
             .num("Cohesion (MPa)", [](EfemState::Row& r) -> double& { return r.cMPa; });
            t.draw("efrac", S.rows, ImGui::GetFontSize() * 8);
        }
        ImGui::EndDisabled();
        if (!busy) {
            ImGui::BeginDisabled(S.rows.empty());
            if (ui::accentButton(ICON_PLAY " Solve")) {
                auto in = std::make_shared<efem::Input>(S.in);
                in->E = S.eGPa * 1e9; in->sxx = S.sxxMPa * 1e6; in->syy = S.syyMPa * 1e6; in->fracs.clear();
                for (const auto& r : S.rows) { efem::Fracture F; F.x0 = r.x0; F.y0 = r.y0; F.x1 = r.x1; F.y1 = r.y1; F.pressure = r.pMPa * 1e6; F.friction = r.mu; F.cohesion = r.cMPa * 1e6; in->fracs.push_back(F); }
                S.busy->store(true);
                auto mx = S.mx; auto bz = S.busy;
                std::thread([mx, bz, in]() {
                    auto r = std::make_shared<efem::Result>();
                    efem::solve(*in, *r);
                    { std::lock_guard<std::mutex> l(*mx); g_ef.res = r; g_ef.resIn = in; }
                    bz->store(false);
                }).detach();
            }
            ImGui::EndDisabled();
        } else ImGui::TextColored(col4(Pal::Warn), "solving ...");
        // ---- C21 G1-c: staggered flow <-> contact-mechanics run of the open case
        ui::header("Coupled flow run with the open case (G1-c)");
        ui::hint("Runs the open black-oil case and, after every report step, loads the fractures (HFRAC) of one layer with their "
                 "simulated pressure profile, solves the contact problem above (same rock, stresses and mesh size over the grid's x-y "
                 "extent) and sets each fracture piece's conductance from the Barton-Bandis aperture: open a = a_res + V_m + w_n, "
                 "closed a = a_res + V_m - V_m s'n / (K_n0 V_m + s'n), plus shear dilation; multiplier (a / a_ref)^3. Lagged by one "
                 "report step; total-stress elasticity (no poroelastic matrix); the map-view opening is scaled by H / L (PKN).");
        {
            auto& C = S.cp;
            const bool cb = C.busy->load();
            ImGui::BeginDisabled(cb);
            if (C.in.layer < 1) C.in.layer = 1;
            if (ui::beginProps("efcp", 0.5f)) {
                ui::propI("Layer (k, 1-based)", &C.in.layer, 1, std::max(1, case_.grid.nz));
                ui::propD("Fracture friction", &C.fricC, "%.3f"); ui::propD("Cohesion (MPa)", &C.cohMPa, "%.4g");
                ui::propD("Residual aperture a_res (mm)", &C.aResMm, "%.4g");
                ui::propD("Max closure V_m (mm)", &C.vmMm, "%.4g");
                ui::propD("Initial normal stiffness K_n0 (GPa/m)", &C.kn0GPaPerM, "%.4g");
                ui::propD("Dilation angle (deg)", &C.dilDeg, "%.3g");
                ui::propD("Reference aperture (mm, 0 = HFRAC width)", &C.aRefMm, "%.4g");
                ui::propD("Fracture height for H/L (m, -1 auto, 0 off)", &C.in.height, "%.4g");
                ui::propBool("Also run uncoupled (BHP comparison)", &C.compare);
                ui::endProps();
            }
            ImGui::EndDisabled();
            const bool can = !case_.grid.hfracs.empty() && !case_.comp.enabled && !case_.thermal.enabled;
            if (!cb) {
                ImGui::BeginDisabled(!can);
                if (ui::accentButton(ICON_PLAY " Run coupled")) {
                    geocouple::Input gi = C.in;
                    gi.layer = std::max(0, C.in.layer - 1);
                    gi.E = S.eGPa * 1e9; gi.nu = S.in.nu; gi.sxx = S.sxxMPa * 1e6; gi.syy = S.syyMPa * 1e6; gi.nx = S.in.nx; gi.ny = S.in.ny;
                    gi.friction = C.fricC; gi.cohesion = C.cohMPa * 1e6;
                    gi.bb.aRes = C.aResMm * 1e-3; gi.bb.Vm = C.vmMm * 1e-3; gi.bb.Kn0 = C.kn0GPaPerM * 1e9; gi.bb.dilationDeg = C.dilDeg;
                    gi.aRef = C.aRefMm * 1e-3;
                    auto base = std::make_shared<Case>(case_);
                    const bool cmp = C.compare;
                    C.busy->store(true); C.cancel->store(false); *C.status = "starting ...";
                    auto bz = C.busy; auto cn = C.cancel; auto st = C.status; auto mx = S.mx;
                    std::thread([base, gi, cmp, bz, cn, st, mx]() {
                        auto R = std::make_shared<EfemState::Coupled::Res>();
                        R->field = base->units == UnitSystem::Field;
                        int wi = 0;
                        for (size_t w = 0; w < base->wells.size(); ++w) if (base->wells[w].type != WellType::Producer) { wi = int(w); break; }
                        if (!base->wells.empty()) R->well = base->wells[size_t(wi)].name;
                        const auto t0 = std::chrono::steady_clock::now();
                        try {
                            Simulator sim; sim.log = [](const std::string&) {};
                            sim.setup(*base); sim.initialize();
                            geocouple::runCoupled(sim, gi, R->out, cn.get(), [st, mx](const std::string& m) { std::lock_guard<std::mutex> l(*mx); *st = m; });
                            for (const auto& p : sim.summary()) if (size_t(wi) < p.wells.size()) { R->t.push_back(p.time); R->bhpC.push_back(p.wells[size_t(wi)].bhp); }
                            if (cmp && R->out.ok && !cn->load()) {
                                { std::lock_guard<std::mutex> l(*mx); *st = "uncoupled comparison run ..."; }
                                Simulator u; u.log = [](const std::string&) {};
                                u.setup(*base); u.initialize();
                                while (!cn->load() && u.runReportStep()) {}
                                for (const auto& p : u.summary()) if (size_t(wi) < p.wells.size()) { R->tU.push_back(p.time); R->bhpU.push_back(p.wells[size_t(wi)].bhp); }
                            }
                        } catch (const std::exception& e) { R->out.ok = false; R->out.error = e.what(); }
                        R->seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                        { std::lock_guard<std::mutex> l(*mx); g_ef.cp.res = R; g_ef.cp.select = true; }
                        bz->store(false);
                    }).detach();
                }
                ImGui::EndDisabled();
                if (!can) ui::hint(case_.grid.hfracs.empty() ? "The open case has no HFRAC fractures (Model > Fractures)." : "Black-oil cases only.");
            } else {
                std::string m; { std::lock_guard<std::mutex> l(*S.mx); m = *C.status; }
                ImGui::TextColored(col4(Pal::Warn), "running: %s", m.c_str());
                if (ImGui::Button("Cancel")) C.cancel->store(true);
            }
        }
        // ---- C23 G3: rate-and-state friction on an embedded fault
        ui::header("Rate-and-state fault: earthquake cycles & induced events (G3)");
        ui::hint("Quasi-dynamic rate-and-state friction (aging law, radiation damping) on one straight fault: a velocity-weakening patch "
                 "inside creeping (velocity-strengthening) ends, loaded at a long-term slip rate, plus pore pressure from an injection "
                 "well or from the last coupled flow run. Initial stresses either uniform or taken from the contact solve above along "
                 "one fracture (the embedded fault). Output: slip-rate history, space-time slip-rate map and an event catalogue (Mw).");
        {
            auto& Q = S.rs;
            const bool rb = Q.busy->load();
            ImGui::BeginDisabled(rb);
            if (ui::beginProps("efrs", 0.55f)) {
                static const char* srcs[] = {"Uniform", "From the contact solve"};
                ui::propCombo("Initial stresses", &Q.stress, srcs, 2, "From the contact solve: sigma_n and |tau| along the chosen fracture (Solve first); tau is capped just below the steady-state strength");
                if (Q.stress == 1) ui::propI("Fracture (1-based)", &Q.frac, 1, std::max(1, int(S.rows.size())));
                else {
                    ui::propD("Fault length (m)", &Q.lengthM, "%.4g");
                    ui::propD("Normal stress (MPa)", &Q.sigMPa, "%.4g");
                    ui::propD("Initial pore pressure (MPa)", &Q.p0MPa, "%.4g");
                }
                ui::propI("Fault elements", &Q.in.n, 8, 512);
                ui::propD("a (weakening patch)", &Q.in.aVW, "%.4g"); ui::propD("a (creeping ends)", &Q.in.aVS, "%.4g");
                ui::propD("b", &Q.in.b, "%.4g"); ui::propD("Dc (mm)", &Q.dcMm, "%.4g");
                ui::propD("Patch from (fraction)", &Q.in.vwFrom, "%.3f"); ui::propD("Patch to (fraction)", &Q.in.vwTo, "%.3f");
                ui::propD("Loading slip rate (mm/yr)", &Q.vplMmYr, "%.4g");
                ui::propD("Shear modulus (GPa)", &Q.gGPa, "%.4g");
                static const char* ps[] = {"None", "Injection well (line source)", "Last coupled flow run (G1-c)"};
                ui::propCombo("Pore pressure", &Q.pressure, ps, 3);
                if (Q.pressure == 1) {
                    ui::propD("Well at (fraction of length)", &Q.injAt, "%.3f"); ui::propD("Well distance from fault (m)", &Q.injOff, "%.4g");
                    ui::propD("Pressure scale q mu / 4 pi k h (MPa)", &Q.injPMPa, "%.4g"); ui::propD("Diffusivity (m2/s)", &Q.injDiff, "%.4g");
                    ui::propD("Shut-in after (days)", &Q.injStopD, "%.4g");
                }
                ui::propD("Simulated time (years)", &Q.years, "%.4g");
                ui::endProps();
            }
            ImGui::EndDisabled();
            if (!rb) {
                std::shared_ptr<efem::Result> R0; std::shared_ptr<efem::Input> I0; std::shared_ptr<EfemState::Coupled::Res> C0;
                { std::lock_guard<std::mutex> l(*S.mx); R0 = S.res; I0 = S.resIn; C0 = S.cp.res; }
                const bool canStress = Q.stress == 0 || (R0 && R0->ok && I0);
                const bool canP = Q.pressure != 2 || (C0 && C0->out.ok);
                ImGui::BeginDisabled(!canStress || !canP);
                if (ui::accentButton(ICON_PLAY " Run fault")) {
                    rsf::Input in = Q.in;
                    const double yr = 3.15576e7;
                    in.G = Q.gGPa * 1e9; in.Dc = Q.dcMm * 1e-3; in.Vpl = Q.vplMmYr * 1e-3 / yr; in.vInit = in.Vpl > 0 ? in.Vpl : 1e-12; in.tEnd = Q.years * yr;
                    in.x0 = 0; in.y0 = 0; in.x1 = Q.lengthM; in.y1 = 0; in.sigmaN = Q.sigMPa * 1e6; in.p0 = Q.p0MPa * 1e6;
                    std::string err;
                    bool ok = true;
                    if (Q.stress == 1) ok = rsf::fromEfem(*I0, *R0, Q.frac - 1, in, err);
                    if (ok && Q.pressure == 1) { in.inj.on = true; in.inj.at = Q.injAt; in.inj.offset = Q.injOff; in.inj.pStar = Q.injPMPa * 1e6; in.inj.diffusivity = Q.injDiff; in.inj.tStop = Q.injStopD * 86400; }
                    if (ok && Q.pressure == 2) ok = rsf::fromCoupled(C0->out, std::max(0, Q.frac - 1), in, err);
                    if (!ok) notify("Rate-and-state: " + err);
                    else {
                        Q.busy->store(true); Q.cancel->store(false);
                        auto bz = Q.busy; auto cn = Q.cancel; auto mx = S.mx;
                        std::thread([in, bz, cn, mx]() {
                            auto r = std::make_shared<rsf::Output>();
                            try { rsf::run(in, *r, cn.get()); } catch (const std::exception& e) { r->ok = false; r->error = e.what(); }
                            { std::lock_guard<std::mutex> l(*mx); g_ef.rs.res = r; g_ef.rs.select = true; }
                            bz->store(false);
                        }).detach();
                    }
                }
                ImGui::EndDisabled();
                if (!canStress) ui::hint("Solve the contact problem first (initial stresses from the contact solve).");
                else if (!canP) ui::hint("Run the coupled flow run first (pore pressure from G1-c).");
            } else {
                ImGui::TextColored(col4(Pal::Warn), "running the fault ...");
                ImGui::SameLine();
                if (ImGui::Button("Cancel##rsf")) Q.cancel->store(true);
            }
        }
        // ---- C23 G4: phase-field fracture
        ui::header("Phase-field fracture: pressurised crack (G4, experimental)");
        ui::hint("Variational AT2 phase-field model of a fluid-pressurised crack in 2-D plane strain (Bourdin-Chukwudozie-Yoshioka): "
                 "the crack is a diffuse band of width ~ l, no crack path is prescribed. The pressure is ramped to 1.3 x the Griffith "
                 "pressure; the run reports the crack-opening volume (vs Sneddon) and the pressure at which the crack jumps (vs Griffith). "
                 "Errors are O(l / a): with l = 2h and l/a = 0.1 expect ~ +40 % volume and ~ -7 % critical pressure; halve h to check. "
                 "C26: any crack orientation, anisotropic in-situ stress (net pressure p - sigma_n and shear on the faces) and VOLUME control - "
                 "the injected volume is ramped and the pressure follows, so growth after the peak is stable and is followed (toughness regime).");
        {
            auto& P = S.pf;
            const bool pb = P.busy->load();
            ImGui::BeginDisabled(pb);
            if (ui::beginProps("efpf", 0.55f)) {
                ui::propD("Initial half-length a (cm)", &P.aCm, "%.4g");
                ui::propD("Young's modulus (GPa)", &P.eGPa, "%.4g");
                ui::propD("Poisson's ratio", &P.in.nu, "%.3f");
                ui::propD("Fracture toughness Gc (J/m2)", &P.in.Gc, "%.4g");
                ui::propD("Element size h (mm)", &P.hMm, "%.4g", "Core element size; l = 2h. Nodes grow like (3a/h)(a/h): h = a/20 ~ 16k nodes, ~4 min on 2 cores");
                ui::propD("Domain half-width / a", &P.in.domain, "%.4g");
                static const char* const ctl[] = {"Pressure ramp (unstable jump at p_c)", "Injected volume (stable growth)"};
                ui::propCombo("Loading", &P.in.control, ctl, 2, "Volume control: p is solved from V(u) = V_injected every staggered iteration");
                ui::propI(P.in.control == 1 ? "Volume steps" : "Pressure steps", &P.in.steps, 8, 200);
                if (P.in.control == 1) ui::propD("Injected volume (cm2 per m, 0 = auto)", &P.vMaxCm2, "%.4g", "0: 6 x the Sneddon volume at the Griffith pressure");
                ui::propD("Crack angle from x (deg)", &P.in.angleDeg, "%.1f", "Initial crack direction; the uniform core box is rotated to contain it");
                ui::propD("In-situ Sxx (MPa, compression +)", &P.sxxM, "%.4g");
                ui::propD("In-situ Syy (MPa, compression +)", &P.syyM, "%.4g");
                ui::propD("In-situ Sxy (MPa)", &P.sxyM, "%.4g");
                ui::endProps();
            }
            ImGui::EndDisabled();
            if (!pb) {
                if (ui::accentButton(ICON_PLAY " Run phase field")) {
                    phasefield::Input in = P.in;
                    in.crackA = P.aCm * 1e-2; in.E = P.eGPa * 1e9; in.h = P.hMm * 1e-3; in.ell = 0;
                    in.sxx = P.sxxM * 1e6; in.syy = P.syyM * 1e6; in.sxy = P.sxyM * 1e6; in.vMax = P.vMaxCm2 * 1e-4;
                    if (!(in.h > 0) || in.h > 0.25 * in.crackA) notify("Phase field: element size must be > 0 and <= a/4");
                    else {
                        P.busy->store(true); P.cancel->store(false); *P.status = "starting";
                        auto bz = P.busy; auto cn = P.cancel; auto mx = S.mx; auto st = P.status;
                        std::thread([in, bz, cn, mx, st]() {
                            auto r = std::make_shared<phasefield::Output>();
                            try {
                                phasefield::run(in, *r, cn.get(), [mx, st](const std::string& m) { std::lock_guard<std::mutex> l(*mx); *st = m; });
                            } catch (const std::exception& e) { r->ok = false; r->error = e.what(); }
                            { std::lock_guard<std::mutex> l(*mx); g_ef.pf.res = r; g_ef.pf.select = true; }
                            bz->store(false);
                        }).detach();
                    }
                }
            } else {
                std::string m; { std::lock_guard<std::mutex> l(*S.mx); m = *P.status; }
                ImGui::TextColored(col4(Pal::Warn), "running: %s", m.c_str());
                ImGui::SameLine();
                if (ImGui::Button("Cancel##pf")) P.cancel->store(true);
            }
        }
        ImGui::EndChild();
        ImGui::TableNextColumn();
        ImGui::BeginChild("efout");
        std::shared_ptr<efem::Result> R; std::shared_ptr<efem::Input> RI;
        { std::lock_guard<std::mutex> l(*S.mx); R = S.res; RI = S.resIn; }
        bool selC = false; { std::lock_guard<std::mutex> l(*S.mx); selC = S.cp.select; S.cp.select = false; }
        if (ImGui::BeginTabBar("efmode")) {
        if (ui::tabItem("Contact solve")) {
        if (!R) ui::hint("Set up the section and the fractures, then Solve. Example: a 30 deg fracture with friction 0.3 under 10 MPa vertical "
                         "compression slips (tan 30 > 0.3); with friction 0.8 it stays locked. Give a fracture a pressure to open it (Sneddon).");
        else if (!R->ok) ImGui::TextColored(col4(Pal::Err), "Failed: %s", R->error.c_str());
        else {
            int nOpen = 0, nStick = 0, nSlip = 0;
            for (const auto& c : R->cut) (c.state == 0 ? nOpen : c.state == 1 ? nStick : nSlip)++;
            ImGui::TextDisabled("%d dofs, %zu fracture elements (open %d, stick %d, slip %d), %d Newton iterations, %.2f s", R->dofs, R->cut.size(), nOpen, nStick, nSlip, R->newton, R->seconds);
            if (ImGui::BeginTabBar("eftabs")) {
                const float ph = std::max(ImGui::GetFontSize() * 15, ImGui::GetContentRegionAvail().y - ImGui::GetFrameHeightWithSpacing() * 1.4f);
                if (ui::tabItem("Along the fractures")) {
                    if (ImPlot::BeginSubplots("##efsub", 2, 1, ImVec2(-1, ph), ImPlotSubplotFlags_LinkAllX)) {
                        if (ImPlot::BeginPlot("Opening and slip")) {
                            ImPlot::SetupAxes("distance from the fracture centre (m)", "mm", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                            for (size_t f = 0; RI && f < RI->fracs.size(); ++f) {
                                std::vector<std::pair<double, size_t>> o;
                                for (size_t k = 0; k < R->cut.size(); ++k) if (R->cut[k].frac == int(f)) o.push_back({R->cut[k].s, k});
                                std::sort(o.begin(), o.end());
                                std::vector<double> s, wn, wt;
                                for (auto& p : o) { s.push_back(p.first); wn.push_back(R->cut[p.second].wn * 1e3); wt.push_back(R->cut[p.second].wt * 1e3); }
                                const std::string nm = f < S.rows.size() ? S.rows[f].name : std::to_string(f + 1);
                                ImPlotSpec a; a.LineColor = col4(seriesColor(int(2 * f))); a.LineWeight = 2 * g_plotLineScale;
                                ImPlot::PlotLine((nm + " opening").c_str(), s.data(), wn.data(), int(s.size()), a);
                                ImPlotSpec b; b.LineColor = col4(seriesColor(int(2 * f + 1))); b.LineWeight = 2 * g_plotLineScale;
                                ImPlot::PlotLine((nm + " slip").c_str(), s.data(), wt.data(), int(s.size()), b);
                            }
                            ImPlot::EndPlot();
                        }
                        if (ImPlot::BeginPlot("Contact tractions")) {
                            ImPlot::SetupAxes("distance from the fracture centre (m)", "MPa", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                            for (size_t f = 0; RI && f < RI->fracs.size(); ++f) {
                                std::vector<std::pair<double, size_t>> o;
                                for (size_t k = 0; k < R->cut.size(); ++k) if (R->cut[k].frac == int(f)) o.push_back({R->cut[k].s, k});
                                std::sort(o.begin(), o.end());
                                std::vector<double> s, tn, tt;
                                for (auto& p : o) { s.push_back(p.first); tn.push_back(R->cut[p.second].tn * 1e-6); tt.push_back(R->cut[p.second].tt * 1e-6); }
                                const std::string nm = f < S.rows.size() ? S.rows[f].name : std::to_string(f + 1);
                                ImPlotSpec a; a.LineColor = col4(seriesColor(int(2 * f))); a.LineWeight = 2 * g_plotLineScale;
                                ImPlot::PlotLine((nm + " normal").c_str(), s.data(), tn.data(), int(s.size()), a);
                                ImPlotSpec b; b.LineColor = col4(seriesColor(int(2 * f + 1))); b.LineWeight = 2 * g_plotLineScale;
                                ImPlot::PlotLine((nm + " shear").c_str(), s.data(), tt.data(), int(s.size()), b);
                            }
                            ImPlot::EndPlot();
                        }
                        ImPlot::EndSubplots();
                    }
                    ImGui::EndTabItem();
                }
                if (ui::tabItem("Stress map")) {
                    static const char* comps[] = {"sigma_xx", "sigma_yy", "sigma_xy", "max shear"};
                    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10);
                    ImGui::Combo("##efcomp", &S.field, comps, 4);
                    const int nx = RI ? RI->nx : 1, ny = RI ? RI->ny : 1;
                    std::vector<double> v(static_cast<size_t>(nx) * ny);
                    double lo = 1e300, hi = -1e300;
                    for (int j = 0; j < ny; ++j)
                        for (int i = 0; i < nx; ++i) {
                            const size_t e = size_t(i + nx * j);
                            const double a = R->sxx[e], b = R->syy[e], c = R->sxy[e];
                            const double val = (S.field == 0 ? a : S.field == 1 ? b : S.field == 2 ? c : std::sqrt(0.25 * (a - b) * (a - b) + c * c)) * 1e-6;
                            v[size_t(i + nx * (ny - 1 - j))] = val;   // heat map rows top-down
                            lo = std::min(lo, val); hi = std::max(hi, val);
                        }
                    if (!(hi > lo)) hi = lo + 1;
                    const float cw = ImGui::GetFontSize() * 5;
                    ImPlot::PushColormap(S.field == 3 ? ImPlotColormap_Viridis : ImPlotColormap_Jet);   // continuous (the default map is categorical)
                    if (ImPlot::BeginPlot("##efmap", ImVec2(ImGui::GetContentRegionAvail().x - cw - ImGui::GetStyle().ItemSpacing.x, ph - ImGui::GetFrameHeightWithSpacing()), ImPlotFlags_Equal)) {
                        ImPlot::SetupAxes("x (m)", "y (m)");
                        ImPlot::PlotHeatmap("##stress", v.data(), ny, nx, lo, hi, nullptr, ImPlotPoint(0, 0), ImPlotPoint(RI->Lx, RI->Ly));
                        for (const auto& F : RI->fracs) {
                            const double xs[2] = {F.x0, F.x1}, ys[2] = {F.y0, F.y1};
                            ImPlotSpec sp; sp.LineColor = ImVec4(1, 1, 1, 1); sp.LineWeight = 2.5f;
                            ImPlot::PlotLine("##fr", xs, ys, 2, sp);
                        }
                        ImPlot::EndPlot();
                    }
                    ImGui::SameLine();
                    ImPlot::ColormapScale("MPa", lo, hi, ImVec2(cw, ph - ImGui::GetFrameHeightWithSpacing()), "%.3g");
                    ImPlot::PopColormap();
                    ImGui::EndTabItem();
                }
                ImGui::EndTabBar();
            }
        }
        ImGui::EndTabItem();
        }
        if (ui::tabItem("Coupled flow run (G1-c)", nullptr, selC ? ImGuiTabItemFlags_SetSelected : 0)) {
            ui::hint("C24: this coupling now also runs inside the normal Run - Grid & rock > Faults > Fracture contact coupling in the run (GEOCOUPLE). This tab stays as a stand-alone study.");
            std::shared_ptr<EfemState::Coupled::Res> C; { std::lock_guard<std::mutex> l(*S.mx); C = S.cp.res; }
            if (!C) ui::hint("Open a case with HFRAC fractures (e.g. examples/HFGROW_INJECTION.DATA), set the stresses on the left "
                             "(sigma_yy acts across a fracture along x) and press Run coupled.");
            else if (!C->out.ok || C->out.steps.empty()) ImGui::TextColored(col4(Pal::Err), "Failed: %s", C->out.error.c_str());
            else {
                const auto& L = C->out.steps.back();
                ImGui::TextDisabled("%zu coupling steps, %.1f s (contact solves %.1f s); last: %d open / %d slipping pieces, max opening %.3g mm, conductance x %.3g .. %.3g",
                                    C->out.steps.size(), C->seconds, C->out.efemSeconds, L.openPieces, L.slipPieces, L.maxOpening * 1e3, L.minMult, L.maxMult);
                const float ph = std::max(ImGui::GetFontSize() * 15, ImGui::GetContentRegionAvail().y - ImGui::GetFrameHeightWithSpacing() * 0.5f);
                if (ImPlot::BeginSubplots("##cpsub", 3, 1, ImVec2(-1, ph))) {
                    if (ImPlot::BeginPlot("Hydraulic aperture along the fractures (last step)")) {
                        ImPlot::SetupAxes("distance from the fracture centre (m)", "mm", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                        int nf = 0; for (int f : C->out.frac) nf = std::max(nf, f + 1);
                        for (int f = 0; f < nf; ++f) {
                            std::vector<std::pair<double, double>> o;
                            for (size_t k = 0; k < C->out.s.size() && k < L.a.size(); ++k) if (C->out.frac[k] == f) o.push_back({C->out.s[k], L.a[k] * 1e3});
                            std::sort(o.begin(), o.end());
                            std::vector<double> x, y; for (auto& q : o) { x.push_back(q.first); y.push_back(q.second); }
                            ImPlotSpec a; a.LineColor = col4(seriesColor(f)); a.LineWeight = 2 * g_plotLineScale; a.Marker = ImPlotMarker_Circle; a.MarkerSize = 3;
                            const std::string nm = f < int(case_.grid.hfracs.size()) ? case_.grid.hfracs[size_t(f)].name : std::to_string(f + 1);
                            ImPlot::PlotLine(nm.c_str(), x.data(), y.data(), int(x.size()), a);
                        }
                        ImPlot::EndPlot();
                    }
                    if (ImPlot::BeginPlot("Conductance multiplier range")) {
                        ImPlot::SetupAxes("time (days)", "(a / a_ref)^3", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                        std::vector<double> t, lo, hi;
                        for (const auto& q : C->out.steps) { t.push_back(q.t); lo.push_back(q.minMult); hi.push_back(q.maxMult); }
                        ImPlotSpec a; a.LineColor = col4(seriesColor(0)); a.LineWeight = 2 * g_plotLineScale;
                        ImPlot::PlotLine("min", t.data(), lo.data(), int(t.size()), a);
                        ImPlotSpec b; b.LineColor = col4(seriesColor(1)); b.LineWeight = 2 * g_plotLineScale;
                        ImPlot::PlotLine("max", t.data(), hi.data(), int(t.size()), b);
                        ImPlot::EndPlot();
                    }
                    if (ImPlot::BeginPlot(("Bottom-hole pressure " + C->well).c_str())) {
                        ImPlot::SetupAxes("time (days)", C->field ? "psia" : "bar", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                        ImPlotSpec a; a.LineColor = col4(seriesColor(0)); a.LineWeight = 2 * g_plotLineScale;
                        ImPlot::PlotLine("coupled", C->t.data(), C->bhpC.data(), int(C->t.size()), a);
                        if (!C->tU.empty()) { ImPlotSpec b; b.LineColor = col4(seriesColor(1)); b.LineWeight = 2 * g_plotLineScale; ImPlot::PlotLine("uncoupled", C->tU.data(), C->bhpU.data(), int(C->tU.size()), b); }
                        ImPlot::EndPlot();
                    }
                    ImPlot::EndSubplots();
                }
            }
            ImGui::EndTabItem();
        }
        bool selR = false; { std::lock_guard<std::mutex> l(*S.mx); selR = S.rs.select; S.rs.select = false; }
        if (ui::tabItem("Rate-and-state fault (G3)", nullptr, selR ? ImGuiTabItemFlags_SetSelected : 0)) {
            ui::hint("C24: faults of the model run rate-and-state on their simulated pressure inside the normal Run - Grid & rock > Faults > Fault mechanics (RSF layer). This tab stays as a stand-alone study.");
            std::shared_ptr<rsf::Output> Q; { std::lock_guard<std::mutex> l(*S.mx); Q = S.rs.res; }
            if (!Q) ui::hint("Set up the fault on the left (defaults: 2 km fault, 1.2 km weakening patch, 31.6 mm/yr loading) and press Run fault. "
                             "Add an injection well to see the clock advance of the next event and the induced catalogue.");
            else if (!Q->ok && Q->snapT.empty()) ImGui::TextColored(col4(Pal::Err), "Failed: %s", Q->error.c_str());
            else {
                const double yr = 3.15576e7;
                ImGui::TextDisabled("%zu events, %ld steps (%ld rejected), %.1f s | nucleation length h* %.0f m, Lb %.1f m, element %.1f m (Lb/dx %.1f)%s%s",
                                    Q->events.size(), Q->steps, Q->rejected, Q->seconds, Q->hStar, Q->Lb, Q->dx, Q->Lb / std::max(1e-12, Q->dx),
                                    Q->error.empty() ? "" : " - stopped: ", Q->error.c_str());
                if (Q->Lb / std::max(1e-12, Q->dx) < 3) ui::hint("Under-resolved: Lb / dx < 3 - use more fault elements or a larger Dc.");
                const float ph = std::max(ImGui::GetFontSize() * 15, ImGui::GetContentRegionAvail().y - ImGui::GetFrameHeightWithSpacing() * 0.5f);
                if (ImGui::BeginTabBar("rstabs")) {
                    if (ui::tabItem("Slip rate & events")) {
                        if (ImPlot::BeginSubplots("##rssub", 2, 1, ImVec2(-1, ph), ImPlotSubplotFlags_LinkAllX)) {
                            if (ImPlot::BeginPlot("Maximum slip rate")) {
                                ImPlot::SetupAxes("time (years)", "max V (m/s)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                                ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                                std::vector<double> t(Q->t.size()); for (size_t i = 0; i < t.size(); ++i) t[i] = Q->t[i] / yr;
                                ImPlotSpec a; a.LineColor = col4(seriesColor(0)); a.LineWeight = 2 * g_plotLineScale;
                                ImPlot::PlotLine("max slip rate", t.data(), Q->vmax.data(), int(t.size()), a);
                                ImPlot::EndPlot();
                            }
                            if (ImPlot::BeginPlot("Event magnitudes")) {
                                ImPlot::SetupAxes("time (years)", "Mw", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                                std::vector<double> t, m; for (const auto& e : Q->events) { t.push_back(e.tStart / yr); m.push_back(e.Mw); }
                                double lo = 1e9; for (double v : m) lo = std::min(lo, v);
                                ImPlotSpec a; a.LineColor = col4(seriesColor(1)); a.Marker = ImPlotMarker_Circle; a.MarkerFillColor = a.LineColor; a.MarkerSize = 4;
                                if (!t.empty()) ImPlot::PlotStems("events", t.data(), m.data(), int(t.size()), lo - 0.5, a);
                                ImPlot::EndPlot();
                            }
                            ImPlot::EndSubplots();
                        }
                        ImGui::EndTabItem();
                    }
                    if (ui::tabItem("Space-time slip rate")) {
                        const int ns = int(Q->snapLogV.size()), nel = int(Q->s.size());
                        if (ns > 1 && nel > 0) {
                            std::vector<double> v(static_cast<size_t>(ns) * nel);
                            for (int r = 0; r < ns; ++r) for (int c = 0; c < nel; ++c) v[size_t(r) * nel + c] = Q->snapLogV[size_t(ns - 1 - r)][size_t(c)];   // rows top-down
                            const double Lf = Q->s.back() + 0.5 * Q->dx;
                            const float cw = ImGui::GetFontSize() * 5;
                            ImPlot::PushColormap(ImPlotColormap_Viridis);
                            if (ImPlot::BeginPlot("Slip rate along the fault##rsmap", ImVec2(ImGui::GetContentRegionAvail().x - cw - ImGui::GetStyle().ItemSpacing.x, ph))) {
                                ImPlot::SetupAxes("along the fault (m)", "snapshot (events sampled densely)");
                                ImPlot::SetupAxesLimits(0, Lf, 0, ns, ImGuiCond_Always);
                                ImPlot::PlotHeatmap("##lv", v.data(), ns, nel, -12, 0, nullptr, ImPlotPoint(0, 0), ImPlotPoint(Lf, ns));
                                ImPlot::EndPlot();
                            }
                            ImGui::SameLine();
                            ImPlot::ColormapScale("log10 V (m/s)", -12, 0, ImVec2(cw, ph), "%.0f");
                            ImPlot::PopColormap();
                        }
                        ImGui::EndTabItem();
                    }
                    if (ui::tabItem("Event catalogue")) {
                        if (ui::beginTable("rscat", 8, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY, ImVec2(0, ph))) {
                            for (const char* c : {"#", "Start (yr)", "Duration (s)", "Mw", "Max slip (m)", "Rupture (m)", "Stress drop (MPa)", "Hypocentre (m)"}) ui::tableColumn(c);
                            ImGui::TableSetupScrollFreeze(0, 1);
                            ImGui::TableHeadersRow();
                            for (size_t k = 0; k < Q->events.size(); ++k) {
                                const auto& e = Q->events[k];
                                ImGui::TableNextRow();
                                ImGui::TableNextColumn(); ImGui::Text("%zu", k + 1);
                                ImGui::TableNextColumn(); ImGui::Text("%.4f", e.tStart / yr);
                                ImGui::TableNextColumn(); ImGui::Text("%.3g", e.tEnd - e.tStart);
                                ImGui::TableNextColumn(); ImGui::Text("%.2f", e.Mw);
                                ImGui::TableNextColumn(); ImGui::Text("%.3g", e.slipMax);
                                ImGui::TableNextColumn(); ImGui::Text("%.0f", e.length);
                                ImGui::TableNextColumn(); ImGui::Text("%.2f", e.stressDrop * 1e-6);
                                ImGui::TableNextColumn(); ImGui::Text("%.0f", e.sHypo);
                            }
                            ui::endTable();
                        }
                        ImGui::EndTabItem();
                    }
                    ImGui::EndTabBar();
                }
            }
            ImGui::EndTabItem();
        }
        bool selP = false; { std::lock_guard<std::mutex> l(*S.mx); selP = S.pf.select; S.pf.select = false; }
        if (ui::tabItem("Phase-field crack (G4)", nullptr, selP ? ImGuiTabItemFlags_SetSelected : 0)) {
            std::shared_ptr<phasefield::Output> Q; { std::lock_guard<std::mutex> l(*S.mx); Q = S.pf.res; }
            if (!Q) ui::hint("Press Run phase field on the left (defaults: Sneddon-Griffith benchmark, a = 10 cm, E = 10 GPa, Gc = 100 J/m2).");
            else if (!Q->ok && Q->steps.empty()) ImGui::TextColored(col4(Pal::Err), "Failed: %s", Q->error.c_str());
            else {
                ImGui::TextDisabled("%d nodes, h %.3g mm, l %.3g mm, %.1f s | Sneddon ratio (first step) %.3f | p_c %s | Griffith %.3f MPa%s%s",
                                    Q->nodes, Q->h * 1e3, Q->ell * 1e3, Q->seconds, Q->sneddonRatio,
                                    Q->pcNumeric > 0 ? (std::to_string(Q->pcBracketLo * 1e-6).substr(0, 5) + "-" + std::to_string(Q->pcBracketHi * 1e-6).substr(0, 5) + " MPa").c_str() : "not reached",
                                    Q->pcTheory * 1e-6, Q->error.empty() ? "" : " - stopped: ", Q->error.c_str());
                if (Q->controlUsed == 1 || Q->sigmaN != 0 || Q->tauS != 0)
                    ImGui::TextDisabled("%s | on the initial crack: sigma_n %.3g MPa, tau %.3g MPa | peak pressure %.4g MPa | + tip direction %.1f deg",
                                        Q->controlUsed == 1 ? "Volume control" : "Pressure control", Q->sigmaN * 1e-6, Q->tauS * 1e-6, Q->pcNumeric * 1e-6, Q->tipAngleDeg);
                const float ph = std::max(ImGui::GetFontSize() * 15, ImGui::GetContentRegionAvail().y - ImGui::GetFrameHeightWithSpacing() * 0.5f);
                if (ImGui::BeginTabBar("pftabs")) {
                    if (ui::tabItem("Pressure - volume")) {
                        if (ImPlot::BeginSubplots("##pfsub", 1, 2, ImVec2(-1, ph))) {
                            std::vector<double> p, V, L, Vs;
                            const double a = Q->x.empty() ? 0.1 : S.pf.aCm * 1e-2;
                            for (const auto& st : Q->steps) { p.push_back(st.p * 1e-6); V.push_back(st.volume * 1e6); L.push_back(st.lengthHalf * 1e2); Vs.push_back(2 * M_PI * (st.p - Q->sigmaN) * a * a / Q->Eprime * 1e6); }
                            if (Q->controlUsed == 1) {   // C26 volume control: pressure and half-length vs injected volume, toughness-regime reference
                                std::vector<double> pT;
                                const double Kc = std::sqrt(S.pf.in.Gc * Q->Eprime);
                                for (const auto& st : Q->steps) pT.push_back((Q->sigmaN + Kc / std::sqrt(M_PI * std::max(a, st.lengthHalf))) * 1e-6);
                                if (ImPlot::BeginPlot("Pressure vs injected volume")) {
                                    ImPlot::SetupAxes("crack volume (cm2 per m)", "pressure (MPa)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                                    ImPlotSpec a1; a1.LineColor = col4(seriesColor(0)); a1.LineWeight = 2 * g_plotLineScale; a1.Marker = ImPlotMarker_Circle; a1.MarkerSize = 3;
                                    ImPlot::PlotLine("phase field", V.data(), p.data(), int(p.size()), a1);
                                    ImPlotSpec a2; a2.LineColor = col4(seriesColor(1)); a2.LineWeight = 1.5f * g_plotLineScale;
                                    ImPlot::PlotLine("sigma_n + K_IC / sqrt(pi a)", V.data(), pT.data(), int(p.size()), a2);
                                    ImPlot::EndPlot();
                                }
                                if (ImPlot::BeginPlot("Crack half-length vs volume")) {
                                    ImPlot::SetupAxes("crack volume (cm2 per m)", "half-length (cm)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                                    ImPlotSpec a1; a1.LineColor = col4(seriesColor(0)); a1.LineWeight = 2 * g_plotLineScale; a1.Marker = ImPlotMarker_Circle; a1.MarkerSize = 3;
                                    ImPlot::PlotLine("half-length (v < 0.5)", V.data(), L.data(), int(p.size()), a1);
                                    ImPlot::EndPlot();
                                }
                                ImPlot::EndSubplots();
                                goto pfPlotsDone;
                            }
                            if (ImPlot::BeginPlot("Crack volume")) {
                                ImPlot::SetupAxes("pressure (MPa)", "volume (cm2 per m)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                                ImPlotSpec a1; a1.LineColor = col4(seriesColor(0)); a1.LineWeight = 2 * g_plotLineScale; a1.Marker = ImPlotMarker_Circle; a1.MarkerSize = 3;
                                ImPlot::PlotLine("phase field", p.data(), V.data(), int(p.size()), a1);
                                ImPlotSpec a2; a2.LineColor = col4(seriesColor(1)); a2.LineWeight = 1.5f * g_plotLineScale;
                                ImPlot::PlotLine("Sneddon 2 pi p a^2 / E'", p.data(), Vs.data(), int(p.size()), a2);
                                if (Q->pcTheory > 0) { double xg = Q->pcTheory * 1e-6; ImPlotSpec g; g.LineColor = col4(seriesColor(2)); ImPlot::PlotInfLines("Griffith p_c", &xg, 1, g); }
                                ImPlot::EndPlot();
                            }
                            if (ImPlot::BeginPlot("Crack half-length")) {
                                ImPlot::SetupAxes("pressure (MPa)", "half-length (cm)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                                ImPlotSpec a1; a1.LineColor = col4(seriesColor(0)); a1.LineWeight = 2 * g_plotLineScale; a1.Marker = ImPlotMarker_Circle; a1.MarkerSize = 3;
                                ImPlot::PlotLine("half-length (v < 0.5)", p.data(), L.data(), int(p.size()), a1);
                                if (Q->pcTheory > 0) { double xg = Q->pcTheory * 1e-6; ImPlotSpec g; g.LineColor = col4(seriesColor(2)); ImPlot::PlotInfLines("Griffith p_c", &xg, 1, g); }
                                ImPlot::EndPlot();
                            }
                            ImPlot::EndSubplots();
                        }
                        pfPlotsDone:
                        ImGui::EndTabItem();
                    }
                    if (ui::tabItem("Phase field")) {
                        static const char* sh[] = {"Final", "At the critical pressure"};
                        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 14);
                        ImGui::Combo("##pfshow", &S.pf.show, sh, 2);
                        const std::vector<double>& vv = (S.pf.show == 1 && !Q->vAtPc.empty()) ? Q->vAtPc : Q->v;
                        const int nx1 = int(Q->x.size()), ny1 = int(Q->y.size());
                        if (nx1 > 1 && ny1 > 1 && vv.size() == size_t(nx1) * ny1) {
                            // resample the graded rectilinear mesh onto a uniform raster over the core region (bilinear)
                            const double a = S.pf.aCm * 1e-2, thp = S.pf.in.angleDeg * M_PI / 180;
                            const double ca = std::fabs(std::cos(thp)), sa = std::fabs(std::sin(thp));
                            const bool incl = sa > 1e-9;
                            const double X = incl ? std::min(Q->x.back(), S.pf.in.coreX * a * ca + S.pf.in.coreY * a * sa) : S.pf.in.coreX * a;
                            const double Y = incl ? std::min(Q->y.back(), std::max(S.pf.in.coreX * a * sa + S.pf.in.coreY * a * ca, 0.3 * X)) : std::min(S.pf.in.coreX * a, Q->y.back());
                            const int rx = 240, ry = std::max(20, int(rx * Y / X));
                            std::vector<double> img(size_t(rx) * ry);
                            auto locate = [](const std::vector<double>& c, double q, int& i, double& t) {
                                i = int(std::upper_bound(c.begin(), c.end(), q) - c.begin()) - 1;
                                i = std::clamp(i, 0, int(c.size()) - 2);
                                t = std::clamp((q - c[size_t(i)]) / (c[size_t(i) + 1] - c[size_t(i)]), 0.0, 1.0);
                            };
                            for (int r = 0; r < ry; ++r) {
                                const double yq = Y - (r + 0.5) * 2 * Y / ry; int j; double ty; locate(Q->y, yq, j, ty);
                                for (int c = 0; c < rx; ++c) {
                                    const double xq = -X + (c + 0.5) * 2 * X / rx; int i; double tx; locate(Q->x, xq, i, tx);
                                    auto V = [&](int ii, int jj) { return vv[size_t(ii) + size_t(nx1) * jj]; };
                                    img[size_t(r) * rx + c] = (1 - ty) * ((1 - tx) * V(i, j) + tx * V(i + 1, j)) + ty * ((1 - tx) * V(i, j + 1) + tx * V(i + 1, j + 1));
                                }
                            }
                            const float cw = ImGui::GetFontSize() * 5;
                            const float h2 = ph - ImGui::GetFrameHeightWithSpacing();
                            ImPlot::PushColormap(ImPlotColormap_Viridis);
                            if (ImPlot::BeginPlot("Phase field v (0 = broken, 1 = intact)##pfmap", ImVec2(ImGui::GetContentRegionAvail().x - cw - ImGui::GetStyle().ItemSpacing.x, h2), ImPlotFlags_Equal)) {
                                ImPlot::SetupAxes("x (m)", "y (m)");
                                ImPlot::PlotHeatmap("##pfv", img.data(), ry, rx, 0, 1, nullptr, ImPlotPoint(-X, -Y), ImPlotPoint(X, Y));
                                double xs[2] = {-a * std::cos(thp), a * std::cos(thp)}, ys[2] = {-a * std::sin(thp), a * std::sin(thp)};
                                ImPlotSpec m; m.LineColor = ImVec4(1, 1, 1, 0.8f); m.LineWeight = 1.5f; m.Marker = ImPlotMarker_Diamond; m.MarkerSize = 4; m.MarkerFillColor = m.LineColor;
                                ImPlot::PlotScatter("initial crack tips", xs, ys, 2, m);
                                ImPlot::EndPlot();
                            }
                            ImGui::SameLine();
                            ImPlot::ColormapScale("v", 0, 1, ImVec2(cw, h2), "%.1f");
                            ImPlot::PopColormap();
                        }
                        ImGui::EndTabItem();
                    }
                    ImGui::EndTabBar();
                }
            }
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
        }
        ImGui::EndChild();
        ImGui::EndTable();
    }
    ImGui::End();
}
