// C40 new-model wizard (G7): four short pages - reservoir, fluid & rock-fluid, wells, schedule - with a live summary,
// then a complete, runnable case (grid, rock, correlation PVT, Corey curves, equilibration, wells, controls, report steps).
// Everything it creates stays editable in the Model workspace; the wizard only removes the blank-page problem.
#include "app.h"
#include "icons.h"
#include "theme.h"
#include "widgets.h"
#include "ressim/tools.h"
#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace ressim;

namespace {
struct Wiz {
    int page = 0;
    char title[96] = "New waterflood model";
    // reservoir
    int nx = 25, ny = 25, nz = 5;
    double dx = 200, dy = 200, thick = 100, top = 8000, dip = 0;
    double perm = 150, kvkh = 0.1, poro = 0.22, sigma = 0.5;
    // fluid
    int fluid = 0;                       // 0 three-phase live oil, 1 dead oil + water
    double api = 35, gasGrav = 0.75, tempF = 200, pinit = 4500, pb = 3500;
    bool gasCap = false; double goc = 7950, owc = 8150;
    double swc = 0.2, sorw = 0.2, nw = 2, no = 2;
    // wells
    int pattern = 2;                     // 0 none, 1 single producer, 2 pair, 3 five-spot, 4 inverted five-spot, 5 line drive
    int injFluid = 0;                    // 0 water, 1 gas
    double qo = 2000, qi = 3000, bhpP = 1500, bhpI = 6000;
    // schedule
    double years = 10; int every = 1;    // 0 monthly, 1 quarterly, 2 half-yearly, 3 yearly
};
Wiz g_wiz;
const char* kPatterns[] = {"No wells yet", "Single producer (depletion)", "Injector - producer pair (corners)", "Five-spot (4 injectors, 1 producer)",
                           "Inverted five-spot (1 injector, 4 producers)", "Line drive (injector row west, producer row east)"};

struct WellSpot { std::string name; WellType type; int i, j; };
std::vector<WellSpot> wellLayout(const Wiz& w) {
    const int nx = w.nx, ny = w.ny;
    const WellType inj = w.injFluid == 1 && w.fluid == 0 ? WellType::GasInjector : WellType::WaterInjector;
    const int ci0 = nx / 2, cj0 = ny / 2, e = std::max(0, std::min(nx, ny) / 10);
    std::vector<WellSpot> s;
    auto add = [&](const std::string& n, WellType t, int i, int j) { s.push_back({n, t, std::clamp(i, 0, nx - 1), std::clamp(j, 0, ny - 1)}); };
    switch (w.pattern) {
        case 1: add("P1", WellType::Producer, ci0, cj0); break;
        case 2: add("I1", inj, e, e); add("P1", WellType::Producer, nx - 1 - e, ny - 1 - e); break;
        case 3:
            add("P1", WellType::Producer, ci0, cj0);
            add("I1", inj, e, e); add("I2", inj, nx - 1 - e, e); add("I3", inj, e, ny - 1 - e); add("I4", inj, nx - 1 - e, ny - 1 - e);
            break;
        case 4:
            add("I1", inj, ci0, cj0);
            add("P1", WellType::Producer, e, e); add("P2", WellType::Producer, nx - 1 - e, e);
            add("P3", WellType::Producer, e, ny - 1 - e); add("P4", WellType::Producer, nx - 1 - e, ny - 1 - e);
            break;
        case 5:
            for (int q = 0; q < 3; ++q) {
                const int j = (2 * q + 1) * ny / 6;
                add("I" + std::to_string(q + 1), inj, e, j);
                add("P" + std::to_string(q + 1), WellType::Producer, nx - 1 - e, j);
            }
            break;
        default: break;
    }
    return s;
}

Case buildWizardCase(const Wiz& w) {
    Case c;
    c.title = w.title;
    c.units = UnitSystem::Field;
    // ---- grid and rock
    GeomInput g; g.nx = w.nx; g.ny = w.ny; g.nz = w.nz; g.dx = w.dx; g.dy = w.dy; g.top = w.top; g.dipX = w.dip; g.dipY = 0;
    g.dz.assign(size_t(w.nz), w.thick / std::max(1, w.nz));
    generateGeometry(g, c.grid);
    RockGenInput r;
    r.layerPerm = {w.perm}; r.layerPoro = {w.poro}; r.layerKvKh = {w.kvkh}; r.layerNtg = {1.0}; r.sigmaLn = w.sigma;
    generateRock(r, c.grid);
    // ---- fluid
    c.phases = Phases{};
    if (w.fluid == 1) { c.phases.gas = false; c.phases.disgas = false; }
    CorrelationInput ci; ci.api = w.api; ci.gasGravity = w.gasGrav; ci.temperatureF = w.tempF;
    ci.pbMax = std::max(w.pb * 1.3, 500.0); ci.pMax = std::max({w.pinit * 1.6, w.bhpI * 1.3, 2000.0});
    makeCorrelationPVT(ci, c.units, c.props, w.fluid == 1);
    CoreyInput k; k.swc = w.swc; k.sorw = w.sorw; k.sorg = w.sorw; k.nw = w.nw; k.now = w.no;
    makeCoreySWOF(k, c.props.swof);
    if (c.phases.gas) makeCoreySGOF(k, c.props.sgof);
    // ---- equilibration: pressure given at the top of the structure
    c.equil.datumDepth = w.top; c.equil.datumPressure = w.pinit;
    c.equil.owc = w.owc;
    c.equil.goc = (c.phases.gas && w.gasCap) ? w.goc : w.top - 200;   // no gas cap: contact above the reservoir
    c.equil.initialRs = -1;
    if (c.phases.gas && !w.gasCap) {
        // undersaturated oil: constant Rs at the requested bubble point (Standing, the correlation family of the table)
        const double yg = 0.00091 * w.tempF - 0.0125 * w.api;
        c.equil.initialRs = w.gasGrav * std::pow((w.pb / 18.2 + 1.4) / std::pow(10.0, yg), 1.2048);   // scf/stb
    }
    // ---- wells
    for (const WellSpot& s : wellLayout(w)) {
        WellInput wi; wi.name = s.name; wi.group = s.type == WellType::Producer ? "PROD" : "INJ"; wi.type = s.type;
        wi.i = s.i; wi.j = s.j;
        Completion cp; cp.i = wi.i; cp.j = wi.j; cp.k1 = 0; cp.k2 = w.nz - 1;
        wi.completions.push_back(cp);
        c.wells.push_back(wi);
    }
    // ---- schedule: controls at the first step, constant report interval
    static const double kLen[] = {30.4375, 91.3125, 182.625, 365.25};
    const double len = kLen[std::clamp(w.every, 0, 3)];
    const int nsteps = std::max(1, int(std::lround(w.years * 365.25 / len)));
    ReportStep rs; rs.length = len;
    c.schedule.assign(size_t(nsteps), rs);
    int np = 0, ni = 0;
    for (const WellInput& wi : c.wells) (wi.type == WellType::Producer ? np : ni)++;
    for (const WellInput& wi : c.wells) {
        WellControl wc; wc.well = wi.name; wc.status = WellStatus::Open;
        if (wi.type == WellType::Producer) { wc.mode = Control::ORAT; wc.orat = w.qo / std::max(1, np); wc.bhp = w.bhpP; }
        else { wc.mode = Control::RATE; wc.rate = w.qi / std::max(1, ni) * (wi.type == WellType::GasInjector ? 1000.0 : 1.0); wc.bhp = w.bhpI; }
        c.schedule[0].controls.push_back(wc);
    }
    return c;
}
}  // namespace

void App::newModelWizard() {
    if (!showWizard_) return;
    Wiz& w = g_wiz;
    const float fs = ImGui::GetFontSize();
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + vp->WorkSize.x * 0.5f, vp->WorkPos.y + vp->WorkSize.y * 0.5f), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSize(ImVec2(std::min(vp->WorkSize.x * 0.9f, fs * 52), std::min(vp->WorkSize.y * 0.9f, fs * 34)), ImGuiCond_Appearing);
    if (!ImGui::IsPopupOpen("New model wizard")) ImGui::OpenPopup("New model wizard");
    bool open = true;
    if (!ImGui::BeginPopupModal("New model wizard", &open, ImGuiWindowFlags_NoSavedSettings)) { if (!open) showWizard_ = false; return; }
    static const char* kPages[] = {ICON_LAYERS "  Reservoir", ICON_DROPLET "  Fluid & rock-fluid", ICON_TARGET "  Wells", ICON_HISTORY "  Schedule"};
    // ---- stepper
    ImGui::BeginChild("##wzsteps", ImVec2(fs * 12, -ImGui::GetFrameHeightWithSpacing() * 1.4f));
    for (int p = 0; p < 4; ++p) {
        ImGui::PushID(p);
        if (p < w.page) ImGui::PushStyleColor(ImGuiCol_Text, col4(Pal::Ok));
        if (ImGui::Selectable(kPages[p], w.page == p, 0, ImVec2(0, fs * 2))) w.page = p;
        if (p < w.page) ImGui::PopStyleColor();
        ImGui::PopID();
    }
    ImGui::Spacing();
    ImGui::SeparatorText("Summary");
    ImGui::PushTextWrapPos(0);
    const int nwells = int(wellLayout(w).size());
    ImGui::TextDisabled("%d x %d x %d = %d cells", w.nx, w.ny, w.nz, w.nx * w.ny * w.nz);
    const double pv = double(w.nx) * w.dx * w.ny * w.dy * w.thick * w.poro / 5.6146;   // rb
    ImGui::TextDisabled("pore volume %.1f MMrb", pv / 1e6);
    ImGui::TextDisabled("%s", w.fluid == 0 ? "3-phase live oil" : "dead oil + water");
    ImGui::TextDisabled("%d wells, %.0f years", nwells, w.years);
    ImGui::PopTextWrapPos();
    ImGui::EndChild();
    ImGui::SameLine();
    // ---- page
    ImGui::BeginChild("##wzpage", ImVec2(0, -ImGui::GetFrameHeightWithSpacing() * 1.4f), ImGuiChildFlags_Borders);
    if (w.page == 0) {
        ui::header("Reservoir geometry and rock");
        ImGui::SetNextItemWidth(-FLT_MIN);
        ImGui::InputText("##title", w.title, sizeof w.title);
        if (ui::beginProps("wz0", 0.45f)) {
            ui::label("Cells in X / Y / Z");
            ImGui::InputInt3("##n", &w.nx);
            w.nx = std::clamp(w.nx, 1, 500); w.ny = std::clamp(w.ny, 1, 500); w.nz = std::clamp(w.nz, 1, 200);
            ui::propD("Cell size X (ft)", &w.dx, "%.0f");
            ui::propD("Cell size Y (ft)", &w.dy, "%.0f");
            ui::propD("Gross thickness (ft)", &w.thick, "%.0f");
            ui::propD("Depth of the top (ft)", &w.top, "%.0f");
            ui::propD("Dip along X (deg)", &w.dip, "%.1f");
            ui::propD("Permeability (mD)", &w.perm, "%.0f");
            ui::propD("Kv / Kh", &w.kvkh, "%.2f");
            ui::propD("Porosity", &w.poro, "%.3f");
            ui::propD("Heterogeneity (std. dev. of ln k)", &w.sigma, "%.2f", "0 = homogeneous; 0.5 moderate; 1 strong. Correlated random field, porosity follows permeability.");
            ui::endProps();
        }
        if (double(w.nx) * w.ny * w.nz > 2e6) ImGui::TextColored(col4(Pal::Warn), ICON_TRIANGLE_ALERT "  more than 2 million cells - runs will be slow on a desktop");
    } else if (w.page == 1) {
        ui::header("Fluid and saturation functions");
        ImGui::RadioButton("Three-phase live oil (water, oil, gas)", &w.fluid, 0);
        ImGui::RadioButton("Dead oil + water", &w.fluid, 1);
        if (ui::beginProps("wz1", 0.45f)) {
            ui::propD("Oil gravity (API)", &w.api, "%.1f");
            if (w.fluid == 0) ui::propD("Gas specific gravity (air = 1)", &w.gasGrav, "%.3f");
            ui::propD("Reservoir temperature (F)", &w.tempF, "%.0f");
            ui::propD("Initial pressure at the top (psia)", &w.pinit, "%.0f");
            if (w.fluid == 0) {
                ui::propD("Bubble-point pressure (psia)", &w.pb, "%.0f", "Undersaturated oil: initial Rs at this bubble point (Standing). With a gas cap the oil is saturated at the contact.");
                ui::label("Gas cap");
                ImGui::Checkbox("##gc", &w.gasCap);
                if (w.gasCap) ui::propD("Gas-oil contact (ft)", &w.goc, "%.0f");
            }
            ui::propD("Oil-water contact (ft)", &w.owc, "%.0f", "Below the base of the reservoir: no aquifer zone in the grid.");
            ui::propD("Connate water Swc", &w.swc, "%.2f");
            ui::propD("Residual oil Sor", &w.sorw, "%.2f");
            ui::propD("Corey exponent water", &w.nw, "%.1f");
            ui::propD("Corey exponent oil", &w.no, "%.1f");
            ui::endProps();
        }
        if (w.fluid == 0 && !w.gasCap && w.pb > w.pinit) ImGui::TextColored(col4(Pal::Warn), ICON_TRIANGLE_ALERT "  bubble point above the initial pressure: the oil starts saturated (free gas forms at once)");
        if (w.owc < w.top) ImGui::TextColored(col4(Pal::Warn), ICON_TRIANGLE_ALERT "  oil-water contact above the top: the reservoir is full of water");
    } else if (w.page == 2) {
        ui::header("Well pattern and controls");
        for (int p = 0; p < 6; ++p) ImGui::RadioButton(kPatterns[p], &w.pattern, p);
        // pattern preview
        {
            const ImVec2 p0 = ImGui::GetCursorScreenPos();
            const float S = fs * 7;
            ImDrawList* dl = ImGui::GetWindowDrawList();
            ImGui::SetCursorScreenPos(ImVec2(p0.x + fs * 22, p0.y - fs * 9.5f));
            const ImVec2 q0 = ImGui::GetCursorScreenPos();
            dl->AddRect(q0, ImVec2(q0.x + S, q0.y + S * float(w.ny * w.dy) / float(std::max(1.0, w.nx * w.dx))), ImGui::GetColorU32(ImGuiCol_Border), 3.0f);
            const float sy = S * float(w.ny * w.dy) / float(std::max(1.0, w.nx * w.dx));
            for (const WellSpot& wi : wellLayout(w)) {
                const ImVec2 c(q0.x + (wi.i + 0.5f) / w.nx * S, q0.y + (wi.j + 0.5f) / w.ny * sy);
                dl->AddCircleFilled(c, fs * 0.3f, wi.type == WellType::Producer ? Pal::Oil : wi.type == WellType::WaterInjector ? Pal::Water : Pal::Gas);
            }
            ImGui::SetCursorScreenPos(p0);
        }
        if (w.pattern >= 2 && w.fluid == 0) { ImGui::TextUnformatted("Injected fluid"); ImGui::SameLine(); ImGui::RadioButton("water", &w.injFluid, 0); ImGui::SameLine(); ImGui::RadioButton("gas", &w.injFluid, 1); }
        if (ui::beginProps("wz2", 0.45f)) {
            if (w.pattern >= 1) {
                ui::propD("Total oil rate target (stb/d)", &w.qo, "%.0f", "Shared equally by the producers");
                ui::propD("Producer minimum BHP (psia)", &w.bhpP, "%.0f");
            }
            if (w.pattern >= 2) {
                ui::propD(w.injFluid == 1 && w.fluid == 0 ? "Total gas injection (Mscf/d)" : "Total water injection (stb/d)", &w.qi, "%.0f", "Shared equally by the injectors");
                ui::propD("Injector maximum BHP (psia)", &w.bhpI, "%.0f");
            }
            ui::endProps();
        }
    } else {
        ui::header("Schedule");
        if (ui::beginProps("wz3", 0.45f)) {
            ui::propD("Simulated period (years)", &w.years, "%.1f");
            w.years = std::clamp(w.years, 0.1, 200.0);
            ui::label("Report interval");
            ImGui::Combo("##ev", &w.every, "monthly\0quarterly\0half-yearly\0yearly\0");
            ui::endProps();
        }
        ImGui::Spacing();
        ImGui::TextWrapped("The model is created with every input in place: grid and rock, correlation PVT (Standing / Vasquez-Beggs type), "
                           "Corey relative permeability, equilibration, wells with full-thickness completions and their controls. "
                           "Press Run (F5) straight away, or refine any page in the Model workspace. Ctrl+Z undoes later edits.");
    }
    ImGui::EndChild();
    // ---- buttons
    if (ImGui::Button(ICON_X " Cancel")) { showWizard_ = false; ImGui::CloseCurrentPopup(); }
    ImGui::SameLine(ImGui::GetContentRegionMax().x - fs * 17);
    ImGui::BeginDisabled(w.page == 0);
    if (ImGui::Button("< Back", ImVec2(fs * 5, 0))) --w.page;
    ImGui::EndDisabled();
    ImGui::SameLine();
    if (w.page < 3) { if (ui::accentButton("Next >", ImVec2(fs * 5, 0))) ++w.page; }
    else if (ui::accentButton(ICON_CHECK " Create", ImVec2(fs * 5, 0))) {
        showWizard_ = false;
        ImGui::CloseCurrentPopup();
        try {
            pendingCase_ = std::make_unique<Case>(buildWizardCase(w));
            if (maybeSaveChanges()) { newCase(4); w.page = 0; notify("Model created - press Run (F5) or refine it in the Model workspace"); }
            else pendingCase_.reset();
        } catch (const std::exception& e) { notify(std::string("Could not create the model: ") + e.what()); }
    }
    ImGui::EndPopup();
    if (!open) showWizard_ = false;
}
