// Shared graphics context owned by the entry point.
//   Windows (default)          : D3D11 + Direct2D + DirectWrite (main.cpp)          -> RESSIM_D3D11
//   Linux / macOS (or -DRESSIM_GL): GLFW window + OpenGL 3.3 core (main_glfw.cpp)   -> RESSIM_GL
#pragma once
#if defined(_WIN32) && !defined(RESSIM_GL)
#define RESSIM_D3D11 1
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <d3d11.h>
#include <dxgi1_2.h>
#include <d2d1_1.h>
#include <dwrite.h>

struct Gfx {
    HWND hwnd = nullptr;
    ID3D11Device* dev = nullptr;
    ID3D11DeviceContext* ctx = nullptr;
    IDXGISwapChain* swap = nullptr;
    ID3D11RenderTargetView* backRTV = nullptr;
    ID2D1Factory1* d2dFactory = nullptr;
    ID2D1Device* d2dDevice = nullptr;
    ID2D1DeviceContext* d2d = nullptr;
    IDWriteFactory* dwrite = nullptr;
    float dpiScale = 1.0f;
};
template <class T> inline void SafeRelease(T*& p) { if (p) { p->Release(); p = nullptr; } }
#else
#ifndef RESSIM_GL
#define RESSIM_GL 1
#endif
struct GLFWwindow;
struct Gfx {
    GLFWwindow* window = nullptr;
    float dpiScale = 1.0f;
    int glMajor = 3, glMinor = 3;
    const char* renderer = "";
};
#endif
extern Gfx g_gfx;
