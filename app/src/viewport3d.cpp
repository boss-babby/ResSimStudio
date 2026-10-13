#include "viewport3d.h"
#include <d3dcompiler.h>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <algorithm>
#include <string>
#include <numeric>
#include <wincodec.h>

// ----------------------------------------------------------------------------- shaders
static const char* kShader = R"(
cbuffer CB : register(b0) {
    float4x4 VP;
    float4 center;      // xyz model centre, w = vertical exaggeration
    float4 range;       // x = vmin, y = 1/(vmax-vmin), z = edge strength, w = log scale flag
    float4 lightDir;    // xyz world-space light direction, w = discrete colour bands (0 = continuous)
    float4 eye;         // xyz camera position (world), w = cell opacity
};
Buffer<float> vals : register(t0);
Texture1D cmap : register(t1);
SamplerState samp : register(s0);

float3 toWorld(float3 p) { return float3(p.x - center.x, -(p.z - center.z) * center.w, p.y - center.y); }
float3 nrmWorld(float3 n) { return normalize(float3(n.x, -n.z / center.w, n.y)); }

struct VSIn  { float3 pos : POSITION; float3 nrm : NORMAL; float2 uv : TEXCOORD0; uint cell : TEXCOORD1; };
struct PSIn  { float4 pos : SV_POSITION; float3 nrm : NORMAL; float2 uv : TEXCOORD0; float val : TEXCOORD1; float3 w : TEXCOORD2; };

PSIn VSMain(VSIn i) {
    PSIn o;
    float3 w = toWorld(i.pos);
    o.pos = mul(float4(w, 1), VP);
    o.nrm = nrmWorld(i.nrm);
    o.uv = i.uv;
    o.val = vals.Load(i.cell);
    o.w = w;
    return o;
}

float3 shade(float3 base, float3 n, float3 w) {
    float3 v = normalize(eye.xyz - w);
    if (dot(n, v) < 0) n = -n;
    float d1 = saturate(dot(n, normalize(lightDir.xyz)));
    float d2 = saturate(dot(n, v));
    float3 h = normalize(normalize(lightDir.xyz) + v);
    float sp = pow(saturate(dot(n, h)), 48) * 0.12;
    return base * (0.30 + 0.45 * d1 + 0.30 * d2) + sp;
}

float4 PSMain(PSIn i) : SV_TARGET {
    float3 base;
    if (i.val < -1e29) base = float3(0.45, 0.47, 0.50);
    else {
        float v = i.val;
        if (range.w > 0.5) v = log10(max(v, 1e-12));
        float t = saturate((v - range.x) * range.y);
        if (lightDir.w > 1.5) t = (min(floor(t * lightDir.w), lightDir.w - 1) + 0.5) / lightDir.w;
        base = cmap.SampleLevel(samp, t, 0).rgb;
    }
    float3 c = shade(base, normalize(i.nrm), i.w);
    float2 fw = max(fwidth(i.uv), 1e-6);
    float2 a = smoothstep(float2(0, 0), fw * 1.3, i.uv) * smoothstep(float2(0, 0), fw * 1.3, 1 - i.uv);
    float edge = (1 - min(a.x, a.y)) * range.z;
    c = lerp(c, c * 0.25, edge);
    return float4(c, saturate(eye.w + edge * 0.35 * (1 - eye.w)));
}

struct VSInC { float3 pos : POSITION; float3 nrm : NORMAL; float4 col : COLOR; };
struct PSInC { float4 pos : SV_POSITION; float3 nrm : NORMAL; float4 col : COLOR; float3 w : TEXCOORD0; };
PSInC VSColor(VSInC i) {
    PSInC o; float3 w = toWorld(i.pos);
    o.pos = mul(float4(w, 1), VP); o.nrm = (dot(i.nrm, i.nrm) > 0.01) ? nrmWorld(i.nrm) : float3(0, 0, 0);
    o.col = i.col; o.w = w; return o;
}
float4 PSColor(PSInC i) : SV_TARGET {
    if (dot(i.nrm, i.nrm) < 0.01) return i.col;
    return float4(shade(i.col.rgb, normalize(i.nrm), i.w), i.col.a);
}
)";

struct CBData { float vp[16]; float center[4]; float range[4]; float light[4]; float eye[4]; };

static ID3DBlob* compile(const char* entry, const char* target) {
    ID3DBlob* code = nullptr; ID3DBlob* err = nullptr;
    HRESULT hr = D3DCompile(kShader, std::strlen(kShader), "viewport.hlsl", nullptr, nullptr, entry, target,
                            D3DCOMPILE_OPTIMIZATION_LEVEL3, 0, &code, &err);
    if (FAILED(hr)) {
        if (err) { OutputDebugStringA((const char*)err->GetBufferPointer()); MessageBoxA(nullptr, (const char*)err->GetBufferPointer(), "Shader error", MB_ICONERROR); err->Release(); }
        return nullptr;
    }
    SafeRelease(err);
    return code;
}

// ----------------------------------------------------------------------------- colormaps
void colormapRGB(Colormap m, float t, float rgb[3]) {
    t = std::min(1.0f, std::max(0.0f, t));
    switch (m) {
        case Colormap::Turbo: {   // polynomial approximation (Google, Apache-2.0)
            const float t2 = t * t, t3 = t2 * t, t4 = t3 * t, t5 = t4 * t;
            rgb[0] = 0.13572138f + 4.61539260f * t - 42.66032258f * t2 + 132.13108234f * t3 - 152.94239396f * t4 + 59.28637943f * t5;
            rgb[1] = 0.09140261f + 2.19418839f * t + 4.84296658f * t2 - 14.18503333f * t3 + 4.27729857f * t4 + 2.82956604f * t5;
            rgb[2] = 0.10667330f + 12.64194608f * t - 60.58204836f * t2 + 110.36276771f * t3 - 89.90310912f * t4 + 27.34824973f * t5;
        } break;
        case Colormap::Viridis: {
            const float c0[3] = {0.2777273272234177f, 0.005407344544966578f, 0.3340998053353061f};
            const float c1[3] = {0.1050930431085774f, 1.404613529898575f, 1.384590162594685f};
            const float c2[3] = {-0.3308618287255563f, 0.214847559468213f, 0.09509516302823659f};
            const float c3[3] = {-4.634230498983486f, -5.799100973351585f, -19.33244095627987f};
            const float c4[3] = {6.228269936347081f, 14.17993336680509f, 56.69055260068105f};
            const float c5[3] = {4.776384997670288f, -13.74514537774601f, -65.35303263337234f};
            const float c6[3] = {-5.435455855934631f, 4.645852612178535f, 26.3124352495832f};
            for (int k = 0; k < 3; ++k) rgb[k] = c0[k] + t * (c1[k] + t * (c2[k] + t * (c3[k] + t * (c4[k] + t * (c5[k] + t * c6[k])))));
        } break;
        case Colormap::Jet: {
            rgb[0] = std::min(1.0f, std::max(0.0f, 1.5f - std::fabs(4.0f * t - 3.0f)));
            rgb[1] = std::min(1.0f, std::max(0.0f, 1.5f - std::fabs(4.0f * t - 2.0f)));
            rgb[2] = std::min(1.0f, std::max(0.0f, 1.5f - std::fabs(4.0f * t - 1.0f)));
        } break;
        case Colormap::CoolWarm: {
            const float a[3] = {0.230f, 0.299f, 0.754f}, b[3] = {0.865f, 0.865f, 0.865f}, c[3] = {0.706f, 0.016f, 0.150f};
            for (int k = 0; k < 3; ++k) rgb[k] = t < 0.5f ? a[k] + (b[k] - a[k]) * (t * 2) : b[k] + (c[k] - b[k]) * (t * 2 - 1);
        } break;
        case Colormap::Grey: rgb[0] = rgb[1] = rgb[2] = 0.15f + 0.8f * t; break;
    }
    for (int k = 0; k < 3; ++k) rgb[k] = std::min(1.0f, std::max(0.0f, rgb[k]));
}

static void normalize3(float v[3]) { const float l = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); if (l > 0) { v[0] /= l; v[1] /= l; v[2] /= l; } }
static void cross3(const float a[3], const float b[3], float o[3]) { o[0] = a[1] * b[2] - a[2] * b[1]; o[1] = a[2] * b[0] - a[0] * b[2]; o[2] = a[0] * b[1] - a[1] * b[0]; }

// ----------------------------------------------------------------------------- init / shutdown
bool Viewport3D::init() {
    ID3D11Device* d = g_gfx.dev;
    ID3DBlob* vsb = compile("VSMain", "vs_5_0"); ID3DBlob* psb = compile("PSMain", "ps_5_0");
    ID3DBlob* vcb = compile("VSColor", "vs_5_0"); ID3DBlob* pcb = compile("PSColor", "ps_5_0");
    if (!vsb || !psb || !vcb || !pcb) return false;
    d->CreateVertexShader(vsb->GetBufferPointer(), vsb->GetBufferSize(), nullptr, &vs_);
    d->CreatePixelShader(psb->GetBufferPointer(), psb->GetBufferSize(), nullptr, &ps_);
    d->CreateVertexShader(vcb->GetBufferPointer(), vcb->GetBufferSize(), nullptr, &vsC_);
    d->CreatePixelShader(pcb->GetBufferPointer(), pcb->GetBufferSize(), nullptr, &psC_);
    const D3D11_INPUT_ELEMENT_DESC il[] = {
        {"POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0},
        {"NORMAL", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D11_INPUT_PER_VERTEX_DATA, 0},
        {"TEXCOORD", 0, DXGI_FORMAT_R32G32_FLOAT, 0, 24, D3D11_INPUT_PER_VERTEX_DATA, 0},
        {"TEXCOORD", 1, DXGI_FORMAT_R32_UINT, 0, 32, D3D11_INPUT_PER_VERTEX_DATA, 0}};
    d->CreateInputLayout(il, 4, vsb->GetBufferPointer(), vsb->GetBufferSize(), &il_);
    const D3D11_INPUT_ELEMENT_DESC ilc[] = {
        {"POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0},
        {"NORMAL", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D11_INPUT_PER_VERTEX_DATA, 0},
        {"COLOR", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 0, 24, D3D11_INPUT_PER_VERTEX_DATA, 0}};
    d->CreateInputLayout(ilc, 3, vcb->GetBufferPointer(), vcb->GetBufferSize(), &ilC_);
    vsb->Release(); psb->Release(); vcb->Release(); pcb->Release();

    D3D11_BUFFER_DESC bd = {};
    bd.ByteWidth = sizeof(CBData); bd.Usage = D3D11_USAGE_DYNAMIC; bd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
    bd.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
    d->CreateBuffer(&bd, nullptr, &cb_);

    D3D11_SAMPLER_DESC sd = {};
    sd.Filter = D3D11_FILTER_MIN_MAG_MIP_LINEAR;
    sd.AddressU = sd.AddressV = sd.AddressW = D3D11_TEXTURE_ADDRESS_CLAMP;
    sd.MaxLOD = D3D11_FLOAT32_MAX;
    d->CreateSamplerState(&sd, &samp_);
    D3D11_RASTERIZER_DESC rd = {};
    rd.FillMode = D3D11_FILL_SOLID; rd.CullMode = D3D11_CULL_NONE; rd.DepthClipEnable = TRUE; rd.MultisampleEnable = TRUE;
    rd.AntialiasedLineEnable = TRUE;
    d->CreateRasterizerState(&rd, &rs_);
    D3D11_DEPTH_STENCIL_DESC dd = {};
    dd.DepthEnable = TRUE; dd.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ALL; dd.DepthFunc = D3D11_COMPARISON_LESS_EQUAL;
    d->CreateDepthStencilState(&dd, &dss_);
    D3D11_BLEND_DESC bld = {};
    bld.RenderTarget[0].BlendEnable = TRUE;
    bld.RenderTarget[0].SrcBlend = D3D11_BLEND_SRC_ALPHA; bld.RenderTarget[0].DestBlend = D3D11_BLEND_INV_SRC_ALPHA;
    bld.RenderTarget[0].BlendOp = D3D11_BLEND_OP_ADD;
    bld.RenderTarget[0].SrcBlendAlpha = D3D11_BLEND_ONE; bld.RenderTarget[0].DestBlendAlpha = D3D11_BLEND_INV_SRC_ALPHA;
    bld.RenderTarget[0].BlendOpAlpha = D3D11_BLEND_OP_ADD;
    bld.RenderTarget[0].RenderTargetWriteMask = D3D11_COLOR_WRITE_ENABLE_ALL;
    d->CreateBlendState(&bld, &bs_);

    UINT q = 0;
    if (FAILED(d->CheckMultisampleQualityLevels(DXGI_FORMAT_B8G8R8A8_UNORM, 4, &q)) || q == 0) msaa_ = 1;

    // DirectWrite text formats
    const float s = g_gfx.dpiScale;
    g_gfx.dwrite->CreateTextFormat(L"Segoe UI", nullptr, DWRITE_FONT_WEIGHT_NORMAL, DWRITE_FONT_STYLE_NORMAL,
                                   DWRITE_FONT_STRETCH_NORMAL, 12.0f * s, L"en-us", &tfSmall_);
    g_gfx.dwrite->CreateTextFormat(L"Segoe UI", nullptr, DWRITE_FONT_WEIGHT_SEMI_BOLD, DWRITE_FONT_STYLE_NORMAL,
                                   DWRITE_FONT_STRETCH_NORMAL, 13.0f * s, L"en-us", &tfLabel_);
    g_gfx.dwrite->CreateTextFormat(L"Segoe UI", nullptr, DWRITE_FONT_WEIGHT_SEMI_BOLD, DWRITE_FONT_STYLE_NORMAL,
                                   DWRITE_FONT_STRETCH_NORMAL, 17.0f * s, L"en-us", &tfTitle_);
    return vs_ && ps_ && il_ && ilC_ && cb_;
}

void Viewport3D::shutdown() {
    releaseTargets();
    SafeRelease(vs_); SafeRelease(ps_); SafeRelease(vsC_); SafeRelease(psC_); SafeRelease(il_); SafeRelease(ilC_);
    SafeRelease(cb_); SafeRelease(vb_); SafeRelease(ib_); SafeRelease(wvb_); SafeRelease(wib_); SafeRelease(lvb_); SafeRelease(slvb_); SafeRelease(dssOff_); SafeRelease(dssNoWrite_);
    SafeRelease(valSRV_); SafeRelease(valBuf_); SafeRelease(cmapSRV_); SafeRelease(cmapTex_); SafeRelease(samp_);
    SafeRelease(rs_); SafeRelease(dss_); SafeRelease(bs_);
    SafeRelease(tfSmall_); SafeRelease(tfLabel_); SafeRelease(tfTitle_);
}

void Viewport3D::releaseTargets() {
    if (g_gfx.d2d) g_gfx.d2d->SetTarget(nullptr);
    SafeRelease(d2dTarget_); SafeRelease(colorSRV_); SafeRelease(color_);
    SafeRelease(msDSV_); SafeRelease(msDepth_); SafeRelease(msRTV_); SafeRelease(msColor_);
    w_ = h_ = 0;
}

void Viewport3D::createTargets(int w, int h) {
    releaseTargets();
    ID3D11Device* d = g_gfx.dev;
    D3D11_TEXTURE2D_DESC td = {};
    td.Width = w; td.Height = h; td.MipLevels = 1; td.ArraySize = 1; td.Format = DXGI_FORMAT_B8G8R8A8_UNORM;
    td.SampleDesc.Count = msaa_; td.Usage = D3D11_USAGE_DEFAULT; td.BindFlags = D3D11_BIND_RENDER_TARGET;
    d->CreateTexture2D(&td, nullptr, &msColor_);
    d->CreateRenderTargetView(msColor_, nullptr, &msRTV_);
    td.Format = DXGI_FORMAT_D32_FLOAT; td.BindFlags = D3D11_BIND_DEPTH_STENCIL;
    d->CreateTexture2D(&td, nullptr, &msDepth_);
    d->CreateDepthStencilView(msDepth_, nullptr, &msDSV_);
    td.Format = DXGI_FORMAT_B8G8R8A8_UNORM; td.SampleDesc.Count = 1;
    td.BindFlags = D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE;
    d->CreateTexture2D(&td, nullptr, &color_);
    d->CreateShaderResourceView(color_, nullptr, &colorSRV_);
    // Direct2D bitmap on the resolved texture
    IDXGISurface* surf = nullptr;
    color_->QueryInterface(IID_PPV_ARGS(&surf));
    if (surf) {
        D2D1_BITMAP_PROPERTIES1 bp = {};
        bp.pixelFormat = {DXGI_FORMAT_B8G8R8A8_UNORM, D2D1_ALPHA_MODE_PREMULTIPLIED};
        bp.dpiX = bp.dpiY = 96.0f;
        bp.bitmapOptions = D2D1_BITMAP_OPTIONS_TARGET | D2D1_BITMAP_OPTIONS_CANNOT_DRAW;
        g_gfx.d2d->CreateBitmapFromDxgiSurface(surf, &bp, &d2dTarget_);
        surf->Release();
    }
    w_ = w; h_ = h;
}

// ----------------------------------------------------------------------------- data
void Viewport3D::setGrid(const ressim::Grid* g) {
    grid_ = g;
    geomDirty_ = true; wellsDirty_ = true; valuesDirty_ = true;
    if (g) {
        center_[0] = 0.5f * (g->bbMin.x + g->bbMax.x);
        center_[1] = 0.5f * (g->bbMin.y + g->bbMax.y);
        center_[2] = 0.5f * (g->bbMin.z + g->bbMax.z);
    }
}

void Viewport3D::setValues(const std::vector<float>& v, float vmin, float vmax, const std::string& title, const std::string& units) {
    values_ = v; vmin_ = vmin; vmax_ = vmax; title_ = title; units_ = units;
    valuesDirty_ = true;
    if (set_.filterValue) geomDirty_ = true;
}

void Viewport3D::worldOf(float x, float y, float z, float o[3]) const {
    o[0] = x - center_[0]; o[1] = -(z - center_[2]) * set_.zScale; o[2] = y - center_[1];
}

void Viewport3D::buildGeometry() {
    SafeRelease(vb_); SafeRelease(ib_); SafeRelease(lvb_); indexCount_ = 0; lineCount_ = 0;
    if (!grid_) return;
    const ressim::Grid& g = *grid_;
    const int nx = g.nx, ny = g.ny, nz = g.nz;
    const ViewSettings& s = set_;
    const int i0 = std::max(0, s.iMin), i1 = std::min(nx - 1, s.iMax);
    const int j0 = std::max(0, s.jMin), j1 = std::min(ny - 1, s.jMax);
    const int k0 = std::max(0, s.kMin), k1 = std::min(nz - 1, s.kMax);
    auto visible = [&](int i, int j, int k) -> bool {
        if (i < i0 || i > i1 || j < j0 || j > j1 || k < k0 || k > k1) return false;
        const int a = g.activeOf[i + nx * (j + ny * k)];
        if (a < 0) return false;
        if (s.filterValue && a < int(values_.size())) {
            const float v = values_[a];
            if (!(v >= s.filterMin && v <= s.filterMax)) return false;
        }
        return true;
    };
    std::vector<Vtx> vtx; std::vector<unsigned> idx;
    faceCenters_.clear(); sortedCam_[0] = 1e30f;
    vtx.reserve(size_t(nx) * ny * 8); idx.reserve(size_t(nx) * ny * 12);
    // face corner indices (into the 8 corners) and model-space normals
    static const int F[6][4] = {{0, 2, 6, 4}, {1, 5, 7, 3}, {0, 4, 5, 1}, {2, 3, 7, 6}, {0, 1, 3, 2}, {4, 6, 7, 5}};
    static const float N[6][3] = {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};
    static const int D[6][3] = {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};
    static const float UV[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    for (int k = k0; k <= k1; ++k)
        for (int j = j0; j <= j1; ++j)
            for (int i = i0; i <= i1; ++i) {
                if (!visible(i, j, k)) continue;
                const int q = i + nx * (j + ny * k);
                const unsigned cell = unsigned(g.activeOf[q]);
                const auto& cr = g.corners[q];
                for (int f = 0; f < 6; ++f) {
                    if (visible(i + D[f][0], j + D[f][1], k + D[f][2])) continue;
                    const unsigned base = unsigned(vtx.size());
                    float fc[3] = {0, 0, 0};
                    for (int c = 0; c < 4; ++c) {
                        const auto& p = cr[F[f][c]];
                        vtx.push_back({{p.x, p.y, p.z}, {N[f][0], N[f][1], N[f][2]}, {UV[c][0], UV[c][1]}, cell});
                        fc[0] += 0.25f * p.x; fc[1] += 0.25f * p.y; fc[2] += 0.25f * p.z;
                    }
                    faceCenters_.insert(faceCenters_.end(), {fc[0], fc[1], fc[2]});
                    idx.insert(idx.end(), {base, base + 1, base + 2, base, base + 2, base + 3});
                }
            }
    // extent for camera
    extent_ = std::max({g.bbMax.x - g.bbMin.x, g.bbMax.y - g.bbMin.y, (g.bbMax.z - g.bbMin.z) * set_.zScale, 1.0f});
    if (!vtx.empty()) {
        D3D11_BUFFER_DESC bd = {};
        bd.ByteWidth = UINT(vtx.size() * sizeof(Vtx)); bd.Usage = D3D11_USAGE_IMMUTABLE; bd.BindFlags = D3D11_BIND_VERTEX_BUFFER;
        D3D11_SUBRESOURCE_DATA sd = {vtx.data(), 0, 0};
        g_gfx.dev->CreateBuffer(&bd, &sd, &vb_);
        bd.ByteWidth = UINT(idx.size() * sizeof(unsigned)); bd.BindFlags = D3D11_BIND_INDEX_BUFFER;
        bd.Usage = D3D11_USAGE_DEFAULT;          // re-ordered back-to-front for transparency
        sd.pSysMem = idx.data();
        g_gfx.dev->CreateBuffer(&bd, &sd, &ib_);
        indexCount_ = unsigned(idx.size());
        sortedIdx_ = std::move(idx);
    }
    // bounding box lines
    {
        const float x0 = g.bbMin.x, x1 = g.bbMax.x, y0 = g.bbMin.y, y1 = g.bbMax.y, z0 = g.bbMin.z, z1 = g.bbMax.z;
        const float P[8][3] = {{x0, y0, z0}, {x1, y0, z0}, {x0, y1, z0}, {x1, y1, z0}, {x0, y0, z1}, {x1, y0, z1}, {x0, y1, z1}, {x1, y1, z1}};
        const int E[12][2] = {{0, 1}, {2, 3}, {4, 5}, {6, 7}, {0, 2}, {1, 3}, {4, 6}, {5, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
        std::vector<VtxC> lv;
        for (auto& e : E)
            for (int t = 0; t < 2; ++t) {
                VtxC v = {{P[e[t]][0], P[e[t]][1], P[e[t]][2]}, {0, 0, 0}, {0.55f, 0.6f, 0.7f, 0.55f}};
                lv.push_back(v);
            }
        D3D11_BUFFER_DESC bd = {};
        bd.ByteWidth = UINT(lv.size() * sizeof(VtxC)); bd.Usage = D3D11_USAGE_IMMUTABLE; bd.BindFlags = D3D11_BIND_VERTEX_BUFFER;
        D3D11_SUBRESOURCE_DATA sd = {lv.data(), 0, 0};
        g_gfx.dev->CreateBuffer(&bd, &sd, &lvb_);
        lineCount_ = unsigned(lv.size());
    }
    geomDirty_ = false;
}

void Viewport3D::buildWellGeometry() {
    SafeRelease(wvb_); SafeRelease(wib_); wIndexCount_ = 0;
    wellsDirty_ = false;
    builtWellZ_ = set_.zScale;
    if (!grid_ || wells_.empty()) return;
    const ressim::Grid& g = *grid_;
    const float cellSize = std::min((g.bbMax.x - g.bbMin.x) / std::max(1, g.nx), (g.bbMax.y - g.bbMin.y) / std::max(1, g.ny));
    const float r = 0.22f * cellSize * std::max(0.1f, set_.wellRadius);
    const float zs = std::max(1e-6f, set_.zScale);
    std::vector<VtxC> v; std::vector<unsigned> id;
    const int seg = 14;
    // tube between two model points; the cross-section is circular in WORLD space (unaffected by the vertical
    // exaggeration), vertices are mapped back to model coordinates for the shared shader
    auto toModel = [&](const float w[3], float m[3]) { m[0] = w[0] + center_[0]; m[1] = w[2] + center_[1]; m[2] = center_[2] - w[1] / zs; };
    auto nrmModel = [&](const float nw[3], float n[3]) { n[0] = nw[0]; n[1] = nw[2]; n[2] = -nw[1] * zs; };
    auto tube = [&](const std::array<float, 3>& A, const std::array<float, 3>& B, float rad, const float col[4], bool caps,
                    bool extA = false, bool extB = false) {
        float a[3], b[3]; worldOf(A[0], A[1], A[2], a); worldOf(B[0], B[1], B[2], b);
        float ax[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        const float len = std::sqrt(ax[0] * ax[0] + ax[1] * ax[1] + ax[2] * ax[2]);
        if (len < 1e-6f) return;
        normalize3(ax);
        // overlap at bends so consecutive segments join without gaps
        if (extA) for (int k = 0; k < 3; ++k) a[k] -= ax[k] * rad * 0.95f;
        if (extB) for (int k = 0; k < 3; ++k) b[k] += ax[k] * rad * 0.95f;
        const float ref[3] = {std::fabs(ax[1]) < 0.9f ? 0.f : 1.f, std::fabs(ax[1]) < 0.9f ? 1.f : 0.f, 0.f};
        float u[3]; cross3(ax, ref, u); normalize3(u);
        float w[3]; cross3(ax, u, w); normalize3(w);
        const unsigned base = unsigned(v.size());
        for (int k = 0; k <= seg; ++k) {
            const float t = 2.0f * 3.14159265f * k / seg, c = std::cos(t), sn = std::sin(t);
            const float nw[3] = {u[0] * c + w[0] * sn, u[1] * c + w[1] * sn, u[2] * c + w[2] * sn};
            float n[3]; nrmModel(nw, n);
            for (const float* e : {a, b}) {
                const float pw[3] = {e[0] + rad * nw[0], e[1] + rad * nw[1], e[2] + rad * nw[2]};
                float pm[3]; toModel(pw, pm);
                v.push_back({{pm[0], pm[1], pm[2]}, {n[0], n[1], n[2]}, {col[0], col[1], col[2], col[3]}});
            }
        }
        for (int k = 0; k < seg; ++k) { const unsigned q = base + 2 * k; id.insert(id.end(), {q, q + 1, q + 2, q + 2, q + 1, q + 3}); }
        if (!caps) return;
        for (int end = 0; end < 2; ++end) {
            const float* e = end ? b : a;
            const float nw[3] = {end ? ax[0] : -ax[0], end ? ax[1] : -ax[1], end ? ax[2] : -ax[2]};
            float n[3]; nrmModel(nw, n);
            const unsigned c0 = unsigned(v.size());
            float pm[3]; toModel(e, pm);
            v.push_back({{pm[0], pm[1], pm[2]}, {n[0], n[1], n[2]}, {col[0], col[1], col[2], col[3]}});
            for (int k = 0; k <= seg; ++k) {
                const float t = 2.0f * 3.14159265f * k / seg, c = std::cos(t), sn = std::sin(t);
                const float pw[3] = {e[0] + rad * (u[0] * c + w[0] * sn), e[1] + rad * (u[1] * c + w[1] * sn), e[2] + rad * (u[2] * c + w[2] * sn)};
                toModel(pw, pm);
                v.push_back({{pm[0], pm[1], pm[2]}, {n[0], n[1], n[2]}, {col[0], col[1], col[2], col[3]}});
            }
            for (int k = 0; k < seg; ++k) id.insert(id.end(), {c0, c0 + 1 + unsigned(k), c0 + 2 + unsigned(k)});
        }
    };
    for (int wi = 0; wi < int(wells_.size()); ++wi) {
        const auto& wl = wells_[size_t(wi)];
        const bool focus = set_.focusWell >= 0, me = wi == set_.focusWell;
        const float dim = focus && !me ? 0.35f : 1.0f;
        const float rs = me ? 1.6f : 1.0f;
        float pipe[4] = {0.78f * dim, 0.80f * dim, 0.84f * dim, 1.0f};
        float colP[4] = {0.18f, 0.80f, 0.44f, 1.0f}, colI[4] = {0.20f, 0.60f, 0.96f, 1.0f}, colX[4] = {0.45f, 0.45f, 0.48f, 1.0f};
        if (me) { pipe[0] = 1.0f; pipe[1] = 0.85f; pipe[2] = 0.30f; }
        for (float* cc : {colP, colI, colX}) for (int k = 0; k < 3; ++k) cc[k] *= dim;
        const float* perfCol = !wl.open ? colX : wl.injector ? colI : colP;
        for (size_t k = 1; k < wl.path.size(); ++k)
            tube(wl.path[k - 1], wl.path[k], r * 0.5f * rs, pipe, true, k > 1, k + 1 < wl.path.size());
        for (const auto& pf : wl.perfs) tube(pf.a, pf.b, r * rs, pf.open ? perfCol : colX, true);
    }
    if (v.empty()) return;
    D3D11_BUFFER_DESC bd = {};
    bd.ByteWidth = UINT(v.size() * sizeof(VtxC)); bd.Usage = D3D11_USAGE_IMMUTABLE; bd.BindFlags = D3D11_BIND_VERTEX_BUFFER;
    D3D11_SUBRESOURCE_DATA sd = {v.data(), 0, 0};
    g_gfx.dev->CreateBuffer(&bd, &sd, &wvb_);
    bd.ByteWidth = UINT(id.size() * sizeof(unsigned)); bd.BindFlags = D3D11_BIND_INDEX_BUFFER; sd.pSysMem = id.data();
    g_gfx.dev->CreateBuffer(&bd, &sd, &wib_);
    wIndexCount_ = unsigned(id.size());
}

void Viewport3D::uploadValues() {
    valuesDirty_ = false;
    if (!grid_) return;
    const int na = grid_->numActive();
    std::vector<float> v(na, -1e30f);                 // -1e30 = undefined (drawn grey)
    for (int a = 0; a < na && a < int(values_.size()); ++a) if (std::isfinite(values_[a])) v[a] = values_[a];
    if (na > valCap_ || !valBuf_) {
        SafeRelease(valSRV_); SafeRelease(valBuf_);
        D3D11_BUFFER_DESC bd = {};
        bd.ByteWidth = UINT(std::max(1, na) * sizeof(float)); bd.Usage = D3D11_USAGE_DEFAULT; bd.BindFlags = D3D11_BIND_SHADER_RESOURCE;
        g_gfx.dev->CreateBuffer(&bd, nullptr, &valBuf_);
        D3D11_SHADER_RESOURCE_VIEW_DESC sd = {};
        sd.Format = DXGI_FORMAT_R32_FLOAT; sd.ViewDimension = D3D11_SRV_DIMENSION_BUFFER;
        sd.Buffer.FirstElement = 0; sd.Buffer.NumElements = UINT(std::max(1, na));
        g_gfx.dev->CreateShaderResourceView(valBuf_, &sd, &valSRV_);
        valCap_ = na;
    }
    if (na > 0) {
        D3D11_BOX box = {0, 0, 0, UINT(na * sizeof(float)), 1, 1};
        g_gfx.ctx->UpdateSubresource(valBuf_, 0, &box, v.data(), 0, 0);
    }
}

void Viewport3D::updateColormap() {
    if (builtCmap_ == set_.cmap && cmapSRV_) return;
    SafeRelease(cmapSRV_); SafeRelease(cmapTex_);
    unsigned char px[256 * 4];
    for (int i = 0; i < 256; ++i) {
        float rgb[3]; colormapRGB(set_.cmap, i / 255.0f, rgb);
        px[i * 4 + 0] = (unsigned char)(rgb[0] * 255); px[i * 4 + 1] = (unsigned char)(rgb[1] * 255);
        px[i * 4 + 2] = (unsigned char)(rgb[2] * 255); px[i * 4 + 3] = 255;
    }
    D3D11_TEXTURE1D_DESC td = {};
    td.Width = 256; td.MipLevels = 1; td.ArraySize = 1; td.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    td.Usage = D3D11_USAGE_IMMUTABLE; td.BindFlags = D3D11_BIND_SHADER_RESOURCE;
    D3D11_SUBRESOURCE_DATA sd = {px, 0, 0};
    g_gfx.dev->CreateTexture1D(&td, &sd, &cmapTex_);
    g_gfx.dev->CreateShaderResourceView(cmapTex_, nullptr, &cmapSRV_);
    builtCmap_ = set_.cmap;
}

// ----------------------------------------------------------------------------- camera

std::array<float, 16> Viewport3D::viewProj() const {
    const float d = dist_ * extent_;
    const float tgt[3] = {panX_, panY_, panZ_};
    const float eye[3] = {tgt[0] + d * std::cos(pitch_) * std::sin(yaw_), tgt[1] + d * std::sin(pitch_), tgt[2] + d * std::cos(pitch_) * std::cos(yaw_)};
    float z[3] = {eye[0] - tgt[0], eye[1] - tgt[1], eye[2] - tgt[2]}; normalize3(z);
    const float up[3] = {0, 1, 0};
    float x[3]; cross3(up, z, x); normalize3(x);
    float y[3]; cross3(z, x, y);
    const float V[16] = {x[0], x[1], x[2], -(x[0] * eye[0] + x[1] * eye[1] + x[2] * eye[2]),
                         y[0], y[1], y[2], -(y[0] * eye[0] + y[1] * eye[1] + y[2] * eye[2]),
                         z[0], z[1], z[2], -(z[0] * eye[0] + z[1] * eye[1] + z[2] * eye[2]),
                         0, 0, 0, 1};
    const float aspect = h_ > 0 ? float(w_) / float(h_) : 1.0f;
    const float f = 1.0f / std::tan(0.5f * 0.75f);   // ~43 deg vertical FOV
    const float zn = 0.01f * d, zf = d + 4.0f * extent_;
    const float P[16] = {f / aspect, 0, 0, 0, 0, f, 0, 0, 0, 0, zf / (zn - zf), zn * zf / (zn - zf), 0, 0, -1, 0};
    std::array<float, 16> M{};
    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 4; ++c) {
            float s = 0;
            for (int k = 0; k < 4; ++k) s += P[r * 4 + k] * V[k * 4 + c];
            M[r * 4 + c] = s;
        }
    return M;
}

bool Viewport3D::project(float x, float y, float z, float& sx, float& sy) const {
    float w[3]; worldOf(x, y, z, w);
    const auto M = viewProj();
    float c[4];
    for (int r = 0; r < 4; ++r) c[r] = M[r * 4 + 0] * w[0] + M[r * 4 + 1] * w[1] + M[r * 4 + 2] * w[2] + M[r * 4 + 3];
    if (c[3] <= 1e-6f) return false;
    sx = (c[0] / c[3] * 0.5f + 0.5f) * w_;
    sy = (1.0f - (c[1] / c[3] * 0.5f + 0.5f)) * h_;
    return true;
}

void Viewport3D::orbit(float dx, float dy) {
    yaw_ -= dx * 0.008f;
    pitch_ = std::min(1.55f, std::max(-1.55f, pitch_ + dy * 0.008f));
}
void Viewport3D::pan(float dx, float dy) {
    const float d = dist_ * extent_ * 0.0015f;
    // camera right/up vectors
    const float rx = std::cos(yaw_), rz = -std::sin(yaw_);
    panX_ -= dx * d * rx; panZ_ -= dx * d * rz;
    panY_ += dy * d * std::cos(pitch_);
}
void Viewport3D::zoom(float wheel) { dist_ = std::min(20.0f, std::max(0.05f, dist_ * std::pow(0.88f, wheel))); }
void Viewport3D::resetCamera() { yaw_ = 0.7f; pitch_ = 0.55f; dist_ = 2.0f; panX_ = panY_ = panZ_ = 0; }
void Viewport3D::viewPreset(int which) {
    panX_ = panY_ = panZ_ = 0; dist_ = 1.9f;
    switch (which) {
        case 0: yaw_ = 0; pitch_ = 1.55f; break;
        case 1: yaw_ = 0; pitch_ = 0.0f; break;
        case 2: yaw_ = 1.5708f; pitch_ = 0.0f; break;
        default: yaw_ = 0.7f; pitch_ = 0.55f; break;
    }
}

int Viewport3D::pick(float px, float py) const {
    if (!grid_ || w_ <= 0) return -1;
    // inverse of view-projection via ray from the camera through the pixel
    const float d = dist_ * extent_;
    const float tgt[3] = {panX_, panY_, panZ_};
    const float eye[3] = {tgt[0] + d * std::cos(pitch_) * std::sin(yaw_), tgt[1] + d * std::sin(pitch_), tgt[2] + d * std::cos(pitch_) * std::cos(yaw_)};
    float z[3] = {eye[0] - tgt[0], eye[1] - tgt[1], eye[2] - tgt[2]}; normalize3(z);
    const float up[3] = {0, 1, 0};
    float x[3]; cross3(up, z, x); normalize3(x);
    float y[3]; cross3(z, x, y);
    const float aspect = float(w_) / float(h_), th = std::tan(0.5f * 0.75f);
    const float nx = (2.0f * px / w_ - 1.0f) * th * aspect, ny = (1.0f - 2.0f * py / h_) * th;
    float dir[3] = {nx * x[0] + ny * y[0] - z[0], nx * x[1] + ny * y[1] - z[1], nx * x[2] + ny * y[2] - z[2]};
    normalize3(dir);
    const ressim::Grid& g = *grid_;
    const ViewSettings& s = set_;
    float best = FLT_MAX; int bestA = -1;
    for (int a = 0; a < g.numActive(); ++a) {
        const int q = g.globalOf[a];
        int i, j, k; g.ijk(q, i, j, k);
        if (i < s.iMin || i > s.iMax || j < s.jMin || j > s.jMax || k < s.kMin || k > s.kMax) continue;
        if (s.filterValue && a < int(values_.size()) && !(values_[a] >= s.filterMin && values_[a] <= s.filterMax)) continue;
        const auto& cr = g.corners[q];
        float lo[3], hi[3];
        worldOf(cr[0].x, cr[0].y, cr[7].z, lo); worldOf(cr[7].x, cr[7].y, cr[0].z, hi);
        float t0 = 0, t1 = best;
        bool hit = true;
        for (int ax = 0; ax < 3 && hit; ++ax) {
            const float mn = std::min(lo[ax], hi[ax]), mx = std::max(lo[ax], hi[ax]);
            if (std::fabs(dir[ax]) < 1e-12f) { if (eye[ax] < mn || eye[ax] > mx) hit = false; continue; }
            float ta = (mn - eye[ax]) / dir[ax], tb = (mx - eye[ax]) / dir[ax];
            if (ta > tb) std::swap(ta, tb);
            t0 = std::max(t0, ta); t1 = std::min(t1, tb);
            if (t0 > t1) hit = false;
        }
        if (hit && t0 < best) { best = t0; bestA = a; }
    }
    return bestA;
}

// ----------------------------------------------------------------------------- transparency: back-to-front faces
void Viewport3D::sortFaces(bool force) {
    const size_t nf = faceCenters_.size() / 3;
    if (!ib_ || nf == 0 || sortedIdx_.size() != nf * 6) return;
    const float cam[7] = {yaw_, pitch_, dist_, panX_, panY_, panZ_, set_.zScale};
    bool same = !force;
    for (int k = 0; k < 7 && same; ++k) same = std::fabs(cam[k] - sortedCam_[k]) < 1e-5f;
    if (same) return;
    std::memcpy(sortedCam_, cam, sizeof(cam));
    const float d = dist_ * extent_;
    const float eye[3] = {panX_ + d * std::cos(pitch_) * std::sin(yaw_), panY_ + d * std::sin(pitch_), panZ_ + d * std::cos(pitch_) * std::cos(yaw_)};
    std::vector<float> key(nf);
    for (size_t f = 0; f < nf; ++f) {
        float w[3]; worldOf(faceCenters_[3 * f], faceCenters_[3 * f + 1], faceCenters_[3 * f + 2], w);
        const float dx = w[0] - eye[0], dy = w[1] - eye[1], dz = w[2] - eye[2];
        key[f] = dx * dx + dy * dy + dz * dz;
    }
    std::vector<unsigned> order(nf);
    std::iota(order.begin(), order.end(), 0u);
    std::sort(order.begin(), order.end(), [&](unsigned a, unsigned b) { return key[a] > key[b]; });
    for (size_t i = 0; i < nf; ++i) {
        const unsigned b = order[i] * 4;
        unsigned* o = &sortedIdx_[i * 6];
        o[0] = b; o[1] = b + 1; o[2] = b + 2; o[3] = b; o[4] = b + 2; o[5] = b + 3;
    }
    g_gfx.ctx->UpdateSubresource(ib_, 0, nullptr, sortedIdx_.data(), 0, 0);
}

int Viewport3D::pickWell(float px, float py, float tol) const {
    if (!set_.wells) return -1;
    int best = -1; float bd = tol * tol;
    auto segDist2 = [](float ax, float ay, float bx, float by, float x, float y) {
        const float vx = bx - ax, vy = by - ay, l2 = vx * vx + vy * vy;
        float t = l2 > 0 ? ((x - ax) * vx + (y - ay) * vy) / l2 : 0.f;
        t = std::max(0.f, std::min(1.f, t));
        const float dx = ax + t * vx - x, dy = ay + t * vy - y;
        return dx * dx + dy * dy;
    };
    for (int w = 0; w < int(wells_.size()); ++w) {
        const auto& P = wells_[size_t(w)].path;
        if (set_.labels && !P.empty()) {   // the name tag above the wellhead is clickable too
            float hx, hy;
            if (project(P[0][0], P[0][1], P[0][2], hx, hy)) {
                const float half = 10.0f + 4.2f * float(wells_[size_t(w)].name.size());
                if (std::fabs(px - hx) <= half && py >= hy - 32.0f && py <= hy + 2.0f) return w;
            }
        }
        for (size_t k = 1; k < P.size(); ++k) {
            float ax, ay, bx, by;
            if (!project(P[k - 1][0], P[k - 1][1], P[k - 1][2], ax, ay) || !project(P[k][0], P[k][1], P[k][2], bx, by)) continue;
            const float d2 = segDist2(ax, ay, bx, by, px, py);
            if (d2 < bd) { bd = d2; best = w; }
        }
    }
    return best;
}

void Viewport3D::focusWellCamera(int idx) {
    if (idx < 0 || idx >= int(wells_.size()) || !grid_) return;
    float lo[3] = {FLT_MAX, FLT_MAX, FLT_MAX}, hi[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
    auto add = [&](const std::array<float, 3>& p) { float w[3]; worldOf(p[0], p[1], p[2], w); for (int k = 0; k < 3; ++k) { lo[k] = std::min(lo[k], w[k]); hi[k] = std::max(hi[k], w[k]); } };
    const auto& wl = wells_[size_t(idx)];
    for (size_t k = 1; k < wl.path.size(); ++k) add(wl.path[k]);    // skip the wellhead above the reservoir
    for (const auto& pf : wl.perfs) { add(pf.a); add(pf.b); }
    if (lo[0] > hi[0]) return;
    panX_ = 0.5f * (lo[0] + hi[0]); panY_ = 0.5f * (lo[1] + hi[1]); panZ_ = 0.5f * (lo[2] + hi[2]);
    const float span = std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2], 0.30f * extent_});
    dist_ = std::min(20.0f, std::max(0.05f, 2.2f * span / extent_));
}

bool Viewport3D::saveImage(const std::string& path, int scale, std::string* err) {
    auto fail = [&](const char* what) { if (err) *err = what; return false; };
    const int w0 = w_, h0 = h_;
    if (w0 <= 0 || h0 <= 0) return fail("view not rendered yet");
    scale = std::max(1, std::min(4, scale));
    overlayScale_ = float(scale);
    render(w0 * scale, h0 * scale);
    overlayScale_ = 1.0f;
    D3D11_TEXTURE2D_DESC td = {};
    color_->GetDesc(&td);
    td.Usage = D3D11_USAGE_STAGING; td.BindFlags = 0; td.CPUAccessFlags = D3D11_CPU_ACCESS_READ; td.MiscFlags = 0;
    ID3D11Texture2D* st = nullptr;
    if (FAILED(g_gfx.dev->CreateTexture2D(&td, nullptr, &st))) { render(w0, h0); return fail("staging texture"); }
    g_gfx.ctx->CopyResource(st, color_);
    D3D11_MAPPED_SUBRESOURCE ms = {};
    bool ok = false;
    if (SUCCEEDED(g_gfx.ctx->Map(st, 0, D3D11_MAP_READ, 0, &ms))) {
        const UINT W = td.Width, H = td.Height;
        std::vector<unsigned char> px(size_t(W) * H * 4);
        for (UINT y = 0; y < H; ++y) {
            std::memcpy(&px[size_t(y) * W * 4], static_cast<const unsigned char*>(ms.pData) + size_t(y) * ms.RowPitch, size_t(W) * 4);
            for (UINT x = 0; x < W; ++x) px[(size_t(y) * W + x) * 4 + 3] = 255;
        }
        g_gfx.ctx->Unmap(st, 0);
        const HRESULT co = CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED);
        IWICImagingFactory* wic = nullptr; IWICBitmap* bmp = nullptr; IWICStream* stream = nullptr;
        IWICBitmapEncoder* enc = nullptr; IWICBitmapFrameEncode* frame = nullptr;
        HRESULT hr = CoCreateInstance(CLSID_WICImagingFactory, nullptr, CLSCTX_INPROC_SERVER, IID_IWICImagingFactory, (void**)&wic);
        if (SUCCEEDED(hr)) hr = wic->CreateBitmapFromMemory(W, H, GUID_WICPixelFormat32bppBGRA, W * 4, UINT(px.size()), px.data(), &bmp);
        const int n = MultiByteToWideChar(CP_UTF8, 0, path.c_str(), -1, nullptr, 0);
        std::wstring wp(size_t(std::max(1, n)), L'\0');
        MultiByteToWideChar(CP_UTF8, 0, path.c_str(), -1, &wp[0], n);
        WICPixelFormatGUID fmt = GUID_WICPixelFormat32bppBGRA;
        if (SUCCEEDED(hr)) hr = wic->CreateStream(&stream);
        if (SUCCEEDED(hr)) hr = stream->InitializeFromFilename(wp.c_str(), GENERIC_WRITE);
        if (SUCCEEDED(hr)) hr = wic->CreateEncoder(GUID_ContainerFormatPng, nullptr, &enc);
        if (SUCCEEDED(hr)) hr = enc->Initialize(stream, WICBitmapEncoderNoCache);
        if (SUCCEEDED(hr)) hr = enc->CreateNewFrame(&frame, nullptr);
        if (SUCCEEDED(hr)) hr = frame->Initialize(nullptr);
        if (SUCCEEDED(hr)) hr = frame->SetSize(W, H);
        if (SUCCEEDED(hr)) hr = frame->SetPixelFormat(&fmt);
        if (SUCCEEDED(hr)) hr = frame->WriteSource(bmp, nullptr);
        if (SUCCEEDED(hr)) hr = frame->Commit();
        if (SUCCEEDED(hr)) hr = enc->Commit();
        SafeRelease(frame); SafeRelease(enc); SafeRelease(stream); SafeRelease(bmp); SafeRelease(wic);
        if (SUCCEEDED(co)) CoUninitialize();
        ok = SUCCEEDED(hr);
        if (!ok) fail("PNG encoding failed");
    } else fail("map failed");
    st->Release();
    render(w0, h0);
    return ok;
}

// ----------------------------------------------------------------------------- rendering
ID3D11ShaderResourceView* Viewport3D::render(int w, int h) {
    w = std::max(16, w); h = std::max(16, h);
    if (w != w_ || h != h_ || !colorSRV_) createTargets(w, h);
    {
        const ViewSettings& a = set_; const ViewSettings& b = built_;
        if (a.iMin != b.iMin || a.iMax != b.iMax || a.jMin != b.jMin || a.jMax != b.jMax || a.kMin != b.kMin || a.kMax != b.kMax ||
            a.filterValue != b.filterValue || (a.filterValue && (a.filterMin != b.filterMin || a.filterMax != b.filterMax)))
            geomDirty_ = true;
        if (grid_ && a.zScale != builtZ_) {
            extent_ = std::max({grid_->bbMax.x - grid_->bbMin.x, grid_->bbMax.y - grid_->bbMin.y, (grid_->bbMax.z - grid_->bbMin.z) * a.zScale, 1.0f});
            builtZ_ = a.zScale;
        }
    }
    if (geomDirty_) { buildGeometry(); built_ = set_; }
    if (set_.zScale != builtWellZ_ || set_.focusWell != builtFocus_ || set_.wellRadius != builtRadius_) wellsDirty_ = true;
    if (wellsDirty_) { buildWellGeometry(); builtFocus_ = set_.focusWell; builtRadius_ = set_.wellRadius; }
    if (valuesDirty_) uploadValues();
    updateColormap();

    ID3D11DeviceContext* c = g_gfx.ctx;
    const float bg[4] = {0.085f, 0.09f, 0.105f, 1.0f};
    c->ClearRenderTargetView(msRTV_, bg);
    c->ClearDepthStencilView(msDSV_, D3D11_CLEAR_DEPTH, 1.0f, 0);
    c->OMSetRenderTargets(1, &msRTV_, msDSV_);
    D3D11_VIEWPORT vp = {0, 0, float(w), float(h), 0, 1};
    c->RSSetViewports(1, &vp);
    c->RSSetState(rs_);
    c->OMSetDepthStencilState(dss_, 0);
    const float bf[4] = {0, 0, 0, 0};
    c->OMSetBlendState(bs_, bf, 0xffffffff);

    // constants
    CBData cbd = {};
    const auto M = viewProj();
    std::memcpy(cbd.vp, M.data(), sizeof(cbd.vp));
    cbd.center[0] = center_[0]; cbd.center[1] = center_[1]; cbd.center[2] = center_[2]; cbd.center[3] = set_.zScale;
    float lo = vmin_, hi = vmax_;
    if (set_.logScale) { lo = std::log10(std::max(vmin_, 1e-12f)); hi = std::log10(std::max(vmax_, 1e-12f)); }
    cbd.range[0] = lo; cbd.range[1] = (hi > lo) ? 1.0f / (hi - lo) : 0.0f; cbd.range[2] = set_.edges ? 0.85f : 0.0f;
    cbd.range[3] = set_.logScale ? 1.0f : 0.0f;
    cbd.light[0] = 0.4f; cbd.light[1] = 0.85f; cbd.light[2] = 0.35f; cbd.light[3] = float(std::max(0, set_.bands));
    const float d = dist_ * extent_;
    cbd.eye[0] = panX_ + d * std::cos(pitch_) * std::sin(yaw_); cbd.eye[1] = panY_ + d * std::sin(pitch_);
    cbd.eye[2] = panZ_ + d * std::cos(pitch_) * std::cos(yaw_);
    const float opacity = std::min(1.0f, std::max(0.05f, set_.opacity));
    cbd.eye[3] = opacity;
    D3D11_MAPPED_SUBRESOURCE ms;
    if (SUCCEEDED(c->Map(cb_, 0, D3D11_MAP_WRITE_DISCARD, 0, &ms))) { std::memcpy(ms.pData, &cbd, sizeof(cbd)); c->Unmap(cb_, 0); }
    c->VSSetConstantBuffers(0, 1, &cb_);
    c->PSSetConstantBuffers(0, 1, &cb_);

    if (linesDirty_) {
        linesDirty_ = false;
        SafeRelease(slvb_); slCount_ = 0;
        if (!lines_.empty()) {
            std::vector<VtxC> lv(lines_.size());
            for (size_t i = 0; i < lines_.size(); ++i) {
                const LineVertex& L = lines_[i];
                lv[i] = VtxC{{L.x, L.y, L.z}, {0, 0, 0}, {L.r, L.g, L.b, L.a}};
            }
            D3D11_BUFFER_DESC bd = {};
            bd.ByteWidth = UINT(lv.size() * sizeof(VtxC)); bd.Usage = D3D11_USAGE_IMMUTABLE; bd.BindFlags = D3D11_BIND_VERTEX_BUFFER;
            D3D11_SUBRESOURCE_DATA sd = {lv.data(), 0, 0};
            if (SUCCEEDED(g_gfx.dev->CreateBuffer(&bd, &sd, &slvb_))) slCount_ = unsigned(lv.size());
        }
    }
    if (!dssOff_) {
        D3D11_DEPTH_STENCIL_DESC dd = {};
        dd.DepthEnable = FALSE; dd.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ZERO; dd.DepthFunc = D3D11_COMPARISON_ALWAYS;
        g_gfx.dev->CreateDepthStencilState(&dd, &dssOff_);
    }
    if (!dssNoWrite_) {
        D3D11_DEPTH_STENCIL_DESC dd = {};
        dd.DepthEnable = TRUE; dd.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ZERO; dd.DepthFunc = D3D11_COMPARISON_LESS_EQUAL;
        g_gfx.dev->CreateDepthStencilState(&dd, &dssNoWrite_);
    }
    const bool transparent = opacity < 0.999f;
    auto drawCells = [&]() {
        if (!(set_.showCells && vb_ && indexCount_)) return;
        if (transparent) { sortFaces(false); if (dssNoWrite_) c->OMSetDepthStencilState(dssNoWrite_, 0); }
        UINT stride = sizeof(Vtx), off = 0;
        c->IASetInputLayout(il_);
        c->IASetVertexBuffers(0, 1, &vb_, &stride, &off);
        c->IASetIndexBuffer(ib_, DXGI_FORMAT_R32_UINT, 0);
        c->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
        c->VSSetShader(vs_, nullptr, 0);
        c->VSSetShaderResources(0, 1, &valSRV_);
        c->PSSetShader(ps_, nullptr, 0);
        ID3D11ShaderResourceView* ps[2] = {valSRV_, cmapSRV_};
        c->PSSetShaderResources(0, 2, ps);
        c->PSSetSamplers(0, 1, &samp_);
        c->DrawIndexed(indexCount_, 0, 0);
        c->OMSetDepthStencilState(dss_, 0);
    };
    auto drawWells = [&]() {
        if (!(set_.wells && wvb_ && wIndexCount_)) return;
        c->IASetInputLayout(ilC_);
        c->VSSetShader(vsC_, nullptr, 0);
        c->PSSetShader(psC_, nullptr, 0);
        if (set_.wellsXray && !transparent && dssOff_) c->OMSetDepthStencilState(dssOff_, 0);
        UINT stride = sizeof(VtxC), off = 0;
        c->IASetVertexBuffers(0, 1, &wvb_, &stride, &off);
        c->IASetIndexBuffer(wib_, DXGI_FORMAT_R32_UINT, 0);
        c->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
        c->DrawIndexed(wIndexCount_, 0, 0);
        c->OMSetDepthStencilState(dss_, 0);
    };
    auto drawLines = [&]() {
        c->IASetInputLayout(ilC_);
        c->VSSetShader(vsC_, nullptr, 0);
        c->PSSetShader(psC_, nullptr, 0);
        if (set_.bbox && lvb_ && lineCount_) {
            UINT stride = sizeof(VtxC), off = 0;
            c->IASetVertexBuffers(0, 1, &lvb_, &stride, &off);
            c->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_LINELIST);
            c->Draw(lineCount_, 0);
        }
        if (slvb_ && slCount_) {
            UINT stride = sizeof(VtxC), off = 0;
            if (set_.linesXray && dssOff_) c->OMSetDepthStencilState(dssOff_, 0);
            c->IASetVertexBuffers(0, 1, &slvb_, &stride, &off);
            c->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_LINELIST);
            c->Draw(slCount_, 0);
            c->OMSetDepthStencilState(dss_, 0);
        }
    };
    if (transparent) {
        // opaque geometry first (depth written), then the see-through cells sorted back to front
        drawWells();
        drawLines();
        drawCells();
    } else {
        drawCells();
        if (!set_.wellsXray) drawWells();
        drawLines();
        if (set_.wellsXray) drawWells();
    }
    ID3D11ShaderResourceView* nulls[2] = {nullptr, nullptr};
    c->PSSetShaderResources(0, 2, nulls);
    c->VSSetShaderResources(0, 1, nulls);
    // resolve MSAA -> single-sample texture shared with Direct2D
    if (msaa_ > 1) c->ResolveSubresource(color_, 0, msColor_, 0, DXGI_FORMAT_B8G8R8A8_UNORM);
    else c->CopyResource(color_, msColor_);
    ID3D11RenderTargetView* nullRTV = nullptr;
    c->OMSetRenderTargets(1, &nullRTV, nullptr);
    drawOverlay();
    return colorSRV_;
}

// ----------------------------------------------------------------------------- Direct2D overlay
static std::wstring widen(const std::string& s) {
    if (s.empty()) return std::wstring();
    const int n = MultiByteToWideChar(CP_UTF8, 0, s.c_str(), int(s.size()), nullptr, 0);
    std::wstring w(n, L'\0');
    MultiByteToWideChar(CP_UTF8, 0, s.c_str(), int(s.size()), &w[0], n);
    return w;
}

void Viewport3D::drawOverlay() {
    if (!d2dTarget_) return;
    ID2D1DeviceContext* d = g_gfx.d2d;
    d->SetTarget(d2dTarget_);
    d->BeginDraw();
    // snapshots render at a multiple of the view size: scale the whole overlay so it keeps its proportions
    const float S = overlayScale_;
    const float W = float(w_) / S, H = float(h_) / S;
    d->SetTransform(D2D1::Matrix3x2F::Scale(S, S));
    ID2D1SolidColorBrush* text = nullptr; ID2D1SolidColorBrush* panel = nullptr; ID2D1SolidColorBrush* br = nullptr;
    d->CreateSolidColorBrush(D2D1::ColorF(0.93f, 0.94f, 0.96f, 1.0f), &text);
    d->CreateSolidColorBrush(D2D1::ColorF(0.06f, 0.065f, 0.075f, 0.72f), &panel);
    d->CreateSolidColorBrush(D2D1::ColorF(1, 1, 1, 1), &br);
    const float s = g_gfx.dpiScale;
    auto drawText = [&](const std::wstring& t, IDWriteTextFormat* f, float x, float y, float wmax, ID2D1Brush* b) {
        d->DrawText(t.c_str(), UINT32(t.size()), f, D2D1::RectF(x, y, x + wmax, y + 200 * s), b);
    };
    auto textWidth = [&](const std::wstring& t, IDWriteTextFormat* f) {
        IDWriteTextLayout* lay = nullptr; float wd = 0;
        if (SUCCEEDED(g_gfx.dwrite->CreateTextLayout(t.c_str(), UINT32(t.size()), f, 4000, 200, &lay))) {
            DWRITE_TEXT_METRICS m; lay->GetMetrics(&m); wd = m.widthIncludingTrailingWhitespace; lay->Release();
        }
        return wd;
    };

    if (grid_ && text && panel && br) {
        // title + caption (top-left)
        if (!title_.empty()) {
            std::wstring t = widen(title_ + (units_.empty() ? "" : "  [" + units_ + "]"));
            const float tw = textWidth(t, tfTitle_);
            std::wstring cap = widen(caption_);
            const float cw = textWidth(cap, tfSmall_);
            d->FillRoundedRectangle(D2D1::RoundedRect(D2D1::RectF(10 * s, 10 * s, 10 * s + std::max(tw, cw) + 24 * s, 10 * s + (cap.empty() ? 34 : 54) * s), 8 * s, 8 * s), panel);
            drawText(t, tfTitle_, 22 * s, 15 * s, 2000, text);
            if (!cap.empty()) { br->SetColor(D2D1::ColorF(0.65f, 0.7f, 0.78f, 1)); drawText(cap, tfSmall_, 22 * s, 40 * s, 2000, br); }
        }
        // legend (right side)
        if (set_.legend && !title_.empty()) {
            const float lh = std::min(H * 0.55f, 320 * s), lw = 16 * s;
            const float x0 = W - 90 * s, y0 = 0.5f * (H - lh);
            d->FillRoundedRectangle(D2D1::RoundedRect(D2D1::RectF(x0 - 12 * s, y0 - 14 * s, W - 8 * s, y0 + lh + 14 * s), 8 * s, 8 * s), panel);
            const int nseg = set_.bands > 1 ? set_.bands : 96;
            for (int i = 0; i < nseg; ++i) {
                float rgb[3]; colormapRGB(set_.cmap, 1.0f - (i + 0.5f) / nseg, rgb);
                br->SetColor(D2D1::ColorF(rgb[0], rgb[1], rgb[2], 1));
                d->FillRectangle(D2D1::RectF(x0, y0 + lh * i / nseg, x0 + lw, y0 + lh * (i + 1) / nseg + 0.8f), br);
            }
            br->SetColor(D2D1::ColorF(1, 1, 1, 0.35f));
            d->DrawRectangle(D2D1::RectF(x0, y0, x0 + lw, y0 + lh), br, 1.0f);
            const int nt = 6;
            for (int t = 0; t < nt; ++t) {
                const float f = float(t) / (nt - 1);
                float v;
                if (set_.logScale) v = std::pow(10.0f, std::log10(std::max(vmin_, 1e-12f)) + f * (std::log10(std::max(vmax_, 1e-12f)) - std::log10(std::max(vmin_, 1e-12f))));
                else v = vmin_ + f * (vmax_ - vmin_);
                char buf[32];
                const float av = std::fabs(v);
                if (av != 0 && (av < 0.01f || av >= 1e5f)) std::snprintf(buf, sizeof(buf), "%.2e", v);
                else if (av < 10) std::snprintf(buf, sizeof(buf), "%.3f", v);
                else std::snprintf(buf, sizeof(buf), "%.1f", v);
                const float y = y0 + lh * (1.0f - f);
                br->SetColor(D2D1::ColorF(1, 1, 1, 0.6f));
                d->DrawLine(D2D1::Point2F(x0 + lw, y), D2D1::Point2F(x0 + lw + 4 * s, y), br, 1.0f);
                drawText(widen(buf), tfSmall_, x0 + lw + 7 * s, y - 8 * s, 80 * s, text);
            }
        }
        // well labels
        if (set_.wells && set_.labels) {
            for (int wi = 0; wi < int(wells_.size()); ++wi) {
                const auto& wl = wells_[size_t(wi)];
                if (wl.path.empty()) continue;
                if (set_.focusWell >= 0 && wi != set_.focusWell && wells_.size() > 12) continue;   // declutter while focused
                float sx, sy;
                if (!project(wl.path[0][0], wl.path[0][1], wl.path[0][2], sx, sy)) continue;
                sx /= S; sy /= S;
                if (sx < 0 || sy < 0 || sx > W || sy > H) continue;
                const std::wstring t = widen(wl.name);
                const float tw = textWidth(t, tfLabel_);
                const D2D1_RECT_F r = D2D1::RectF(sx - tw * 0.5f - 7 * s, sy - 30 * s, sx + tw * 0.5f + 7 * s, sy - 9 * s);
                br->SetColor(wl.injector ? D2D1::ColorF(0.12f, 0.42f, 0.85f, 0.92f) : D2D1::ColorF(0.10f, 0.55f, 0.30f, 0.92f));
                if (!wl.open) br->SetColor(D2D1::ColorF(0.35f, 0.35f, 0.38f, 0.92f));
                if (wi == set_.focusWell) br->SetColor(D2D1::ColorF(0.92f, 0.62f, 0.10f, 0.97f));
                else if (set_.focusWell >= 0) br->SetColor(D2D1::ColorF(0.25f, 0.27f, 0.30f, 0.80f));
                d->FillRoundedRectangle(D2D1::RoundedRect(r, 5 * s, 5 * s), br);
                d->DrawLine(D2D1::Point2F(sx, sy - 9 * s), D2D1::Point2F(sx, sy), br, 2.0f * s);
                drawText(t, tfLabel_, r.left + 7 * s, r.top + 2 * s, tw + 10 * s, text);
            }
        }
        // axis triad (bottom-left): project unit directions from a virtual origin
        if (set_.axes) {
            const float ox = 55 * s, oy = H - 55 * s, L = 38 * s;
            // camera basis (same construction as viewProj)
            float z[3] = {std::cos(pitch_) * std::sin(yaw_), std::sin(pitch_), std::cos(pitch_) * std::cos(yaw_)};
            const float up[3] = {0, 1, 0};
            float x[3]; cross3(up, z, x); normalize3(x);
            float y[3]; cross3(z, x, y);
            struct Ax { float w[3]; const wchar_t* name; D2D1::ColorF col; };
            const Ax axes[3] = {{{1, 0, 0}, L"X (I)", D2D1::ColorF(0.95f, 0.36f, 0.36f, 1)},
                                {{0, 0, 1}, L"Y (J)", D2D1::ColorF(0.36f, 0.86f, 0.46f, 1)},
                                {{0, -1, 0}, L"Depth", D2D1::ColorF(0.40f, 0.62f, 1.0f, 1)}};
            d->FillEllipse(D2D1::Ellipse(D2D1::Point2F(ox, oy), 48 * s, 48 * s), panel);
            for (const auto& a : axes) {
                const float vx = x[0] * a.w[0] + x[1] * a.w[1] + x[2] * a.w[2];
                const float vy = y[0] * a.w[0] + y[1] * a.w[1] + y[2] * a.w[2];
                const float ex = ox + vx * L, ey = oy - vy * L;
                br->SetColor(a.col);
                d->DrawLine(D2D1::Point2F(ox, oy), D2D1::Point2F(ex, ey), br, 2.2f * s);
                d->FillEllipse(D2D1::Ellipse(D2D1::Point2F(ex, ey), 2.5f * s, 2.5f * s), br);
                drawText(a.name, tfSmall_, ex + 3 * s, ey - 8 * s, 60 * s, br);
            }
        }
    }
    SafeRelease(text); SafeRelease(panel); SafeRelease(br);
    d->EndDraw();
    d->SetTarget(nullptr);
}
