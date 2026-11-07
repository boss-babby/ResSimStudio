// C39 Linux / macOS entry point: GLFW window + OpenGL 3.3 core context, Dear ImGui (GLFW + OpenGL3 back ends).
// The Windows entry point (main.cpp: Win32 + Direct3D 11 + Direct2D + DirectWrite) is unchanged; both drive the same App.
#include "i18n.h"
#include "gfx.h"
#ifdef RESSIM_GL
#include "app.h"
#include "gl_api.h"
#include "gl_view.h"
#include "gpu_compute.h"
#include "platform.h"
#include "theme.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "implot.h"
#include <GLFW/glfw3.h>
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

Gfx g_gfx;
static App* g_app = nullptr;

static void onGlfwError(int code, const char* msg) { std::fprintf(stderr, "GLFW error %d: %s\n", code, msg); }
static void onDrop(GLFWwindow*, int n, const char** paths) { if (n > 0 && g_app) g_app->openFile(paths[0]); }
static void onClose(GLFWwindow* w) { if (g_app && !g_app->confirmClose()) glfwSetWindowShouldClose(w, GLFW_FALSE); }

int main(int argc, char** argv) {
    std::string fileArg, smokePng, smokeAction, langArg;
    bool geo = false, geoDemo = false;
    int smokeFrames = 0, smokeWorkspace = -1;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        // C18 headless linear-solver benchmark:  ResSimStudio --gpu-bench [maxCells] [report.txt]
        if (a == "--gpu-bench") {
            int maxCells = 1000000; std::string file = "gpu_benchmark.txt";
            if (i + 1 < argc && std::isdigit((unsigned char)argv[i + 1][0])) maxCells = std::atoi(argv[++i]);
            if (i + 1 < argc) file = argv[++i];
            const std::string r = gpu::benchmark(maxCells, 0);
            bool ok = false;
            const std::string st = gpu::selfTest(ok);
            std::ofstream(file) << r << "self test: " << (ok ? "PASS " : "FAIL ") << st << "\n";
            std::printf("%s\nself test: %s %s\nreport written to %s\n", r.c_str(), ok ? "PASS" : "FAIL", st.c_str(), file.c_str());
            return ok ? 0 : 2;
        }
        if (a == "--geo") { geo = true; continue; }
        if (a == "--geo-demo") { geoDemo = true; continue; }
        if (a == "--lang" && i + 1 < argc) { langArg = argv[++i]; continue; }   // C43: interface language for this session
        if (a == "--smoke-action" && i + 1 < argc) { smokeAction = argv[++i]; continue; }   // C40
        // C39 CI / headless check:  --smoke-test [frames] [screenshot.png] [workspace 0-4]  renders and exits
        if (a == "--smoke-test") {
            smokeFrames = 90;
            if (i + 1 < argc && std::isdigit((unsigned char)argv[i + 1][0])) smokeFrames = std::max(2, std::atoi(argv[++i]));
            if (i + 1 < argc && argv[i + 1][0] != '-' && std::strstr(argv[i + 1], ".png")) smokePng = argv[++i];
            if (i + 1 < argc && std::isdigit((unsigned char)argv[i + 1][0])) smokeWorkspace = std::atoi(argv[++i]);
            continue;
        }
        if (!a.empty() && a[0] != '-') fileArg = a;
    }
    {   // C32 ResGeo: an executable named ResGeo* starts the standalone geomodelling front end
        std::string en = std::filesystem::path(argc > 0 ? argv[0] : "").filename().string();
        for (char& ch : en) ch = char(std::tolower((unsigned char)ch));
        if (en.rfind("resgeo", 0) == 0) geo = true;
    }

    glfwSetErrorCallback(onGlfwError);
    if (!glfwInit()) { std::fprintf(stderr, "ResSim Studio: GLFW initialisation failed (no display?)\n"); return 1; }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);   // required on macOS
    glfwWindowHint(GLFW_MAXIMIZED, GLFW_TRUE);
#if defined(__APPLE__)
    glfwWindowHint(GLFW_COCOA_RETINA_FRAMEBUFFER, GLFW_TRUE);
#endif
    // UI scale: monitor content scale (X11 Xft.dpi / Wayland / Windows-style scaling); macOS works in points (scale 1)
    float scale = 1.0f;
#if !defined(__APPLE__)
    if (GLFWmonitor* mon = glfwGetPrimaryMonitor()) scale = ImGui_ImplGlfw_GetContentScaleForMonitor(mon);
    if (const char* s = std::getenv("RESSIM_UI_SCALE")) { const float v = float(std::atof(s)); if (v >= 0.5f && v <= 4.0f) scale = v; }
#endif
    scale = std::max(0.75f, std::min(4.0f, scale));
    GLFWwindow* win = glfwCreateWindow(int(1600 * scale), int(960 * scale), geo ? "ResGeo - geomodelling (ResSim Studio)" : "ResSim Studio",
                                       nullptr, nullptr);
    if (!win) {
        std::fprintf(stderr, "ResSim Studio: could not create an OpenGL 3.3 core window\n");
        glfwTerminate();
        return 1;
    }
    glfwMakeContextCurrent(win);
    glfwSwapInterval(1);
    glfwSetWindowSizeLimits(win, int(960 * scale), int(600 * scale), GLFW_DONT_CARE, GLFW_DONT_CARE);
    g_gfx.window = win;
    g_gfx.dpiScale = scale;
    glGetIntegerv(GL_MAJOR_VERSION, &g_gfx.glMajor);
    glGetIntegerv(GL_MINOR_VERSION, &g_gfx.glMinor);
    if (const GLubyte* r = glGetString(GL_RENDERER)) g_gfx.renderer = reinterpret_cast<const char*>(r);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    GImPlotTranslate = &i18n::trs;   // C48: plot axis labels, legends and titles follow the interface language
    GImGuiTitleTranslate = &i18n::trs;   // C49: window titles and dock tabs (display only; IDs and the saved layout unchanged)
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard | ImGuiConfigFlags_DockingEnable;
    io.ConfigWindowsMoveFromTitleBarOnly = true;
    // dock layout in the per-user settings folder (the install folder may be read-only); first start: shipped default
    static std::string iniPath;
    {
        namespace fs = std::filesystem;
        std::error_code ec;
        iniPath = (fs::u8path(plat::configDir()) / "ressim_layout.ini").u8string();
        const fs::path shipped = fs::u8path(plat::exeDir()) / "ressim_layout.ini";
        if (!fs::exists(fs::u8path(iniPath), ec) && fs::exists(shipped, ec)) fs::copy_file(shipped, fs::u8path(iniPath), ec);
    }
    io.IniFilename = iniPath.c_str();
    applyTheme(scale);
    ImGui_ImplGlfw_InitForOpenGL(win, true);
#if defined(__APPLE__)
    ImGui_ImplOpenGL3_Init("#version 150");
#else
    ImGui_ImplOpenGL3_Init("#version 330 core");
#endif
    glfwSetDropCallback(win, onDrop);
    glfwSetWindowCloseCallback(win, onClose);

    App app;
    g_app = &app;
    app.init();
    if (geo) app.setGeoMode(true);
    if (!langArg.empty()) i18n::setLanguage(langArg);
    if (!fileArg.empty()) app.openFile(fileArg);
    if (smokeWorkspace >= 0) app.setWorkspaceIndex(smokeWorkspace);
    app.smokeGeoDemo_ = geoDemo;
    app.smokeAction_ = smokeAction;
    int frame = 0;

    while (!glfwWindowShouldClose(win)) {
        // idle throttling: when nothing is running and there is no input, wait instead of spinning the GPU
        if (app.wantsContinuousFrames() || smokeFrames > 0) glfwPollEvents();
        else glfwWaitEventsTimeout(0.05);
        if (glfwGetWindowAttrib(win, GLFW_ICONIFIED)) { ImGui_ImplGlfw_Sleep(10); continue; }
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        app.frame();
        if (app.quitRequested() && app.confirmClose()) glfwSetWindowShouldClose(win, GLFW_TRUE);
        ImGui::Render();
        int fw = 0, fh = 0;
        glfwGetFramebufferSize(win, &fw, &fh);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glViewport(0, 0, fw, fh);
        glDisable(GL_SCISSOR_TEST);
        glClearColor(0.07f, 0.075f, 0.085f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        if (smokeFrames > 0 && ++frame >= smokeFrames) {
            if (!smokePng.empty()) {
                std::vector<unsigned char> px(size_t(fw) * size_t(fh) * 4), flip(px.size());
                glPixelStorei(GL_PACK_ALIGNMENT, 1);
                glReadPixels(0, 0, fw, fh, GL_RGBA, GL_UNSIGNED_BYTE, px.data());
                for (int y = 0; y < fh; ++y) std::memcpy(&flip[size_t(y) * fw * 4], &px[size_t(fh - 1 - y) * fw * 4], size_t(fw) * 4);
                for (size_t k = 3; k < flip.size(); k += 4) flip[k] = 255;
                const bool ok = glv::writePngRGBA(smokePng, fw, fh, flip.data());
                std::printf("smoke test: %d frames, %dx%d, GL %d.%d (%s), screenshot %s %s\n", frame, fw, fh, g_gfx.glMajor, g_gfx.glMinor,
                            g_gfx.renderer, smokePng.c_str(), ok ? "written" : "FAILED");
            } else std::printf("smoke test: %d frames OK, GL %d.%d (%s)\n", frame, g_gfx.glMajor, g_gfx.glMinor, g_gfx.renderer);
            break;
        }
        glfwSwapBuffers(win);
    }

    app.shutdown();
    g_app = nullptr;
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    glfwDestroyWindow(win);
    glfwTerminate();
    return 0;
}
#endif  // RESSIM_GL
