// Shared graphics context (D3D11 + Direct2D + DirectWrite) owned by main.cpp.
#pragma once
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
extern Gfx g_gfx;

template <class T> inline void SafeRelease(T*& p) { if (p) { p->Release(); p = nullptr; } }
