// ResSim Studio - entry point: Win32 window, D3D11 swap chain, Direct2D/DirectWrite interop, Dear ImGui.
#include "i18n.h"
#include "gfx.h"
#include "app.h"
#include "theme.h"
#include "imgui.h"
#include "implot.h"
#include "backends/imgui_impl_win32.h"
#include "backends/imgui_impl_dx11.h"
#include "gpu_compute.h"
#include <shellapi.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <objbase.h>
#include <string>
#include <cctype>

Gfx g_gfx;
static UINT g_resizeW = 0, g_resizeH = 0;
static bool g_occluded = false;
static App* g_app = nullptr;

extern IMGUI_IMPL_API LRESULT ImGui_ImplWin32_WndProcHandler(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);

static void createBackBufferRTV() {
    ID3D11Texture2D* back = nullptr;
    g_gfx.swap->GetBuffer(0, IID_PPV_ARGS(&back));
    if (back) { g_gfx.dev->CreateRenderTargetView(back, nullptr, &g_gfx.backRTV); back->Release(); }
}

static bool createDevice(HWND hwnd) {
    DXGI_SWAP_CHAIN_DESC sd = {};
    sd.BufferCount = 2;
    sd.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    sd.BufferDesc.RefreshRate.Numerator = 60;
    sd.BufferDesc.RefreshRate.Denominator = 1;
    sd.Flags = DXGI_SWAP_CHAIN_FLAG_ALLOW_MODE_SWITCH;
    sd.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
    sd.OutputWindow = hwnd;
    sd.SampleDesc.Count = 1;
    sd.Windowed = TRUE;
    sd.SwapEffect = DXGI_SWAP_EFFECT_DISCARD;
    UINT flags = D3D11_CREATE_DEVICE_BGRA_SUPPORT;      // required for Direct2D interop
#ifdef _DEBUG
    flags |= D3D11_CREATE_DEVICE_DEBUG;
#endif
    const D3D_FEATURE_LEVEL levels[] = {D3D_FEATURE_LEVEL_11_1, D3D_FEATURE_LEVEL_11_0, D3D_FEATURE_LEVEL_10_0};
    D3D_FEATURE_LEVEL got;
    HRESULT hr = D3D11CreateDeviceAndSwapChain(nullptr, D3D_DRIVER_TYPE_HARDWARE, nullptr, flags, levels, 3, D3D11_SDK_VERSION,
                                               &sd, &g_gfx.swap, &g_gfx.dev, &got, &g_gfx.ctx);
#ifdef _DEBUG
    if (FAILED(hr)) {   // debug layer not installed
        flags &= ~D3D11_CREATE_DEVICE_DEBUG;
        hr = D3D11CreateDeviceAndSwapChain(nullptr, D3D_DRIVER_TYPE_HARDWARE, nullptr, flags, levels, 3, D3D11_SDK_VERSION,
                                           &sd, &g_gfx.swap, &g_gfx.dev, &got, &g_gfx.ctx);
    }
#endif
    if (hr == DXGI_ERROR_UNSUPPORTED)   // fall back to WARP software rasterizer
        hr = D3D11CreateDeviceAndSwapChain(nullptr, D3D_DRIVER_TYPE_WARP, nullptr, flags, levels, 3, D3D11_SDK_VERSION, &sd,
                                           &g_gfx.swap, &g_gfx.dev, &got, &g_gfx.ctx);
    if (FAILED(hr)) return false;
    createBackBufferRTV();

    // Direct2D on the same device
    D2D1_FACTORY_OPTIONS fo = {};
    if (FAILED(D2D1CreateFactory(D2D1_FACTORY_TYPE_SINGLE_THREADED, __uuidof(ID2D1Factory1), &fo, (void**)&g_gfx.d2dFactory))) return false;
    IDXGIDevice* dxgiDev = nullptr;
    g_gfx.dev->QueryInterface(IID_PPV_ARGS(&dxgiDev));
    if (!dxgiDev || FAILED(g_gfx.d2dFactory->CreateDevice(dxgiDev, &g_gfx.d2dDevice))) { SafeRelease(dxgiDev); return false; }
    dxgiDev->Release();
    if (FAILED(g_gfx.d2dDevice->CreateDeviceContext(D2D1_DEVICE_CONTEXT_OPTIONS_NONE, &g_gfx.d2d))) return false;
    g_gfx.d2d->SetTextAntialiasMode(D2D1_TEXT_ANTIALIAS_MODE_GRAYSCALE);
    if (FAILED(DWriteCreateFactory(DWRITE_FACTORY_TYPE_SHARED, __uuidof(IDWriteFactory), (IUnknown**)&g_gfx.dwrite))) return false;
    return true;
}

static void destroyDevice() {
    SafeRelease(g_gfx.dwrite); SafeRelease(g_gfx.d2d); SafeRelease(g_gfx.d2dDevice); SafeRelease(g_gfx.d2dFactory);
    SafeRelease(g_gfx.backRTV); SafeRelease(g_gfx.swap); SafeRelease(g_gfx.ctx); SafeRelease(g_gfx.dev);
}

static LRESULT WINAPI WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    if (ImGui_ImplWin32_WndProcHandler(hWnd, msg, wParam, lParam)) return true;
    switch (msg) {
        case WM_SIZE:
            if (wParam == SIZE_MINIMIZED) return 0;
            g_resizeW = (UINT)LOWORD(lParam); g_resizeH = (UINT)HIWORD(lParam);
            return 0;
        case WM_GETMINMAXINFO: {
            auto* mm = reinterpret_cast<MINMAXINFO*>(lParam);
            mm->ptMinTrackSize.x = 960; mm->ptMinTrackSize.y = 600;
            return 0;
        }
        case WM_DROPFILES: {
            HDROP drop = (HDROP)wParam;
            char path[MAX_PATH];
            if (DragQueryFileA(drop, 0, path, MAX_PATH) && g_app) g_app->openFile(path);
            DragFinish(drop);
            return 0;
        }
        case WM_CLOSE:
            if (g_app && !g_app->confirmClose()) return 0;
            break;
        case WM_SYSCOMMAND:
            if ((wParam & 0xfff0) == SC_KEYMENU) return 0;   // disable ALT application menu
            break;
        case WM_DESTROY:
            PostQuitMessage(0);
            return 0;
    }
    return DefWindowProcW(hWnd, msg, wParam, lParam);
}

int WINAPI WinMain(HINSTANCE hInst, HINSTANCE, LPSTR cmdLine, int) {
    // C18 headless linear-solver benchmark:  ResSimStudio.exe --gpu-bench [maxCells] [report.txt]
    if (cmdLine && std::strstr(cmdLine, "--gpu-bench")) {
        int maxCells = 1000000; std::string file = "gpu_benchmark.txt";
        const char* a = std::strstr(cmdLine, "--gpu-bench") + 11;
        while (*a == ' ') ++a;
        if (*a >= '0' && *a <= '9') { maxCells = std::atoi(a); while (*a && *a != ' ') ++a; while (*a == ' ') ++a; }
        if (*a) { file = a; while (!file.empty() && (file.back() == ' ' || file.back() == '"')) file.pop_back(); if (!file.empty() && file[0] == '"') file.erase(0, 1); }
        const bool con = AttachConsole(ATTACH_PARENT_PROCESS) != 0;
        if (con) { FILE* fo = nullptr; freopen_s(&fo, "CONOUT$", "w", stdout); }
        const std::string r = gpu::benchmark(maxCells, 0);
        bool ok = false;
        const std::string st = gpu::selfTest(ok);
        std::ofstream(file) << r << "self test: " << (ok ? "PASS " : "FAIL ") << st << "\n";
        if (con) { std::printf("%s\nself test: %s %s\nreport written to %s\n", r.c_str(), ok ? "PASS" : "FAIL", st.c_str(), file.c_str()); std::fflush(stdout); }
        return ok ? 0 : 2;
    }
    // the common file dialogs (shell namespace) expect a single-threaded COM apartment on the UI thread
    const HRESULT comInit = CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED | COINIT_DISABLE_OLE1DDE);
    ImGui_ImplWin32_EnableDpiAwareness();
    WNDCLASSEXW wc = {sizeof(wc), CS_CLASSDC, WndProc, 0, 0, hInst, nullptr, LoadCursor(nullptr, IDC_ARROW), nullptr, nullptr,
                      L"ResSimStudio", nullptr};
    wc.hIcon = LoadIcon(hInst, MAKEINTRESOURCE(1));
    wc.hIconSm = (HICON)LoadImageW(hInst, MAKEINTRESOURCEW(1), IMAGE_ICON, GetSystemMetrics(SM_CXSMICON), GetSystemMetrics(SM_CYSMICON), 0);
    if (!wc.hIcon) wc.hIcon = LoadIcon(nullptr, IDI_APPLICATION);
    RegisterClassExW(&wc);
    const float scale = ImGui_ImplWin32_GetDpiScaleForMonitor(MonitorFromPoint(POINT{0, 0}, MONITOR_DEFAULTTOPRIMARY));
    HWND hwnd = CreateWindowW(wc.lpszClassName, L"ResSim Studio", WS_OVERLAPPEDWINDOW, 100, 60, int(1600 * scale),
                              int(960 * scale), nullptr, nullptr, wc.hInstance, nullptr);
    g_gfx.hwnd = hwnd;
    g_gfx.dpiScale = scale;
    if (!createDevice(hwnd)) {
        MessageBoxW(hwnd, L"Could not create a Direct3D 11 device with Direct2D support.", L"ResSim Studio", MB_ICONERROR);
        destroyDevice();
        UnregisterClassW(wc.lpszClassName, wc.hInstance);
        return 1;
    }
    DragAcceptFiles(hwnd, TRUE);
    // dark title bar on Windows 10/11 (DWMWA_USE_IMMERSIVE_DARK_MODE = 20), loaded dynamically
    if (HMODULE dwm = LoadLibraryA("dwmapi.dll")) {
        typedef HRESULT(WINAPI * SetAttr)(HWND, DWORD, LPCVOID, DWORD);
        if (auto f = (SetAttr)(void*)GetProcAddress(dwm, "DwmSetWindowAttribute")) { BOOL on = TRUE; f(hwnd, 20, &on, sizeof(on)); }
    }
    ShowWindow(hwnd, SW_SHOWMAXIMIZED);
    UpdateWindow(hwnd);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    GImPlotTranslate = &i18n::trs;   // C48: plot axis labels, legends and titles follow the interface language
    GImGuiTitleTranslate = &i18n::trs;   // C49: window titles and dock tabs (display only; IDs and the saved layout unchanged)
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard | ImGuiConfigFlags_DockingEnable;
    io.ConfigWindowsMoveFromTitleBarOnly = true;
    io.IniFilename = "ressim_layout.ini";
    applyTheme(scale);
    ImGui_ImplWin32_Init(hwnd);
    ImGui_ImplDX11_Init(g_gfx.dev, g_gfx.ctx);

    App app;
    g_app = &app;
    app.init();
    // C32 ResGeo: 'ResSimStudio.exe --geo' (ResGeo.exe launcher) or an executable named ResGeo*.exe
    std::string cmd = cmdLine ? cmdLine : "";
    {
        char exe[MAX_PATH] = {0};
        GetModuleFileNameA(nullptr, exe, MAX_PATH);
        std::string en(exe);
        const size_t sl = en.find_last_of("\\/");
        if (sl != std::string::npos) en = en.substr(sl + 1);
        for (char& ch : en) ch = char(std::tolower((unsigned char)ch));
        const size_t lg = cmd.find("--lang ");   // C43: --lang <code> = interface language for this session
        if (lg != std::string::npos) {
            size_t a = lg + 7; while (a < cmd.size() && cmd[a] == ' ') ++a;
            size_t b = a; while (b < cmd.size() && cmd[b] != ' ') ++b;
            i18n::setLanguage(cmd.substr(a, b - a));
            cmd.erase(lg, b - lg);
        }
        const size_t g = cmd.find("--geo");
        if (g != std::string::npos || en.rfind("resgeo", 0) == 0) {
            if (g != std::string::npos) cmd.erase(g, 5);
            while (!cmd.empty() && cmd.front() == ' ') cmd.erase(0, 1);
            while (!cmd.empty() && cmd.back() == ' ') cmd.pop_back();
            app.setGeoMode(true);
            SetWindowTextW(hwnd, L"ResGeo - geomodelling (ResSim Studio)");
        }
    }
    if (!cmd.empty()) {
        std::string p(cmd);
        if (p.size() > 1 && p.front() == '"') p = p.substr(1, p.find('"', 1) - 1);
        app.openFile(p);
    }

    bool done = false;
    while (!done) {
        MSG msg;
        while (PeekMessage(&msg, nullptr, 0U, 0U, PM_REMOVE)) {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
            if (msg.message == WM_QUIT) done = true;
        }
        if (done) break;
        if (g_occluded && g_gfx.swap->Present(0, DXGI_PRESENT_TEST) == DXGI_STATUS_OCCLUDED) { Sleep(10); continue; }
        g_occluded = false;
        if (g_resizeW != 0 && g_resizeH != 0) {
            SafeRelease(g_gfx.backRTV);
            g_gfx.swap->ResizeBuffers(0, g_resizeW, g_resizeH, DXGI_FORMAT_UNKNOWN, 0);
            g_resizeW = g_resizeH = 0;
            createBackBufferRTV();
        }
        // idle throttling: when nothing is running and no input, don't spin the GPU
        if (!app.wantsContinuousFrames()) MsgWaitForMultipleObjects(0, nullptr, FALSE, 50, QS_ALLINPUT);

        ImGui_ImplDX11_NewFrame();
        ImGui_ImplWin32_NewFrame();
        ImGui::NewFrame();
        app.frame();
        if (app.quitRequested()) PostMessage(hwnd, WM_CLOSE, 0, 0);
        ImGui::Render();
        const float clear[4] = {0.07f, 0.075f, 0.085f, 1.0f};
        g_gfx.ctx->OMSetRenderTargets(1, &g_gfx.backRTV, nullptr);
        g_gfx.ctx->ClearRenderTargetView(g_gfx.backRTV, clear);
        ImGui_ImplDX11_RenderDrawData(ImGui::GetDrawData());
        const HRESULT hr = g_gfx.swap->Present(1, 0);
        g_occluded = (hr == DXGI_STATUS_OCCLUDED);
    }

    app.shutdown();
    g_app = nullptr;
    ImGui_ImplDX11_Shutdown();
    ImGui_ImplWin32_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    destroyDevice();
    DestroyWindow(hwnd);
    UnregisterClassW(wc.lpszClassName, wc.hInstance);
    if (SUCCEEDED(comInit)) CoUninitialize();
    return 0;
}
