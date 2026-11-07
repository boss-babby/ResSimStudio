// C36 3-D geomodel / seismic scene - see geoscene3d.h
#include "geoscene3d.h"
#include "geoscene3d_cpu.h"
#include "overlay_paint.h"
#ifdef RESSIM_D3D11
#include <d3dcompiler.h>
#include <wincodec.h>
#endif
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <string>

using namespace geo3d;

#ifdef RESSIM_D3D11
namespace {
const char* kGeoShader = R"(
cbuffer Frame : register(b0) {
    float4x4 VP;
    float4x4 invVP;
    float4 center;      // xyz model centre, w = vertical exaggeration
    float4 eye;         // xyz camera (world)
    float4 light;       // xyz light direction (world)
    float4 vpx;         // width, height, 1/width, 1/height
};
cbuffer Obj : register(b1) {
    float4 range;       // lo, 1/(hi-lo), opacity, contour interval
    float4 flags;       // alpha mode, edges, -, -
    float4 volN;        // ni, nj, ns, steps
    float4 volP;        // x0, y0, threshold, zMode
    float4 volA;        // inverse plan matrix a11 a12 a21 a22  (i = a11 dx + a12 dy, j = a21 dx + a22 dy)
    float4 volZ;        // t0, dt, zScale (mode 0), -
    float4 volV;        // v0, k, datum, -
    float4 volB0;       // iMin, jMin, zMin, -
    float4 volB1;       // iMax, jMax, zMax, -
};
Texture1D cmap : register(t1);
Texture2D<float> stex : register(t0);
Texture3D<float> vol : register(t2);
#if MSDEPTH
Texture2DMS<float> sceneDepth : register(t3);
#else
Texture2D<float> sceneDepth : register(t3);
#endif
SamplerState samp : register(s0);

float3 toWorld(float3 p) { return float3(p.x - center.x, -(p.z - center.z) * center.w, p.y - center.y); }
float3 nrmWorld(float3 n) { return normalize(float3(n.x, -n.z / center.w, n.y)); }
bool isNaNf(float v) { return (asuint(v) & 0x7fffffff) > 0x7f800000; }
float3 shade(float3 base, float3 n, float3 w) {
    float3 v = normalize(eye.xyz - w);
    if (dot(n, v) < 0) n = -n;
    float d1 = saturate(dot(n, normalize(light.xyz)));
    float d2 = saturate(dot(n, v));
    float3 h = normalize(normalize(light.xyz) + v);
    float sp = pow(saturate(dot(n, h)), 40) * 0.15;
    return base * (0.32 + 0.45 * d1 + 0.28 * d2) + sp;
}

// ---- value meshes (horizons, grid faces)
struct VSInM { float3 pos : POSITION; float3 nrm : NORMAL; float val : TEXCOORD0; float2 uv : TEXCOORD1; };
struct PSInM { float4 pos : SV_POSITION; float3 nrm : NORMAL; float val : TEXCOORD0; float2 uv : TEXCOORD1; float3 w : TEXCOORD2; };
PSInM VSMesh(VSInM i) {
    PSInM o; float3 w = toWorld(i.pos);
    o.pos = mul(float4(w, 1), VP); o.nrm = nrmWorld(i.nrm); o.val = i.val; o.uv = i.uv; o.w = w; return o;
}
float4 PSMesh(PSInM i) : SV_TARGET {
    float3 base = float3(0.5, 0.52, 0.55);
    if (!isNaNf(i.val)) base = cmap.SampleLevel(samp, saturate((i.val - range.x) * range.y), 0).rgb;
    float3 c = shade(base, normalize(i.nrm), i.w);
    // derivatives evaluated in uniform control flow (contour lines, cell edges)
    float f = i.val / max(range.w, 1e-30);
    float fwf = max(fwidth(f), 1e-6);
    float2 fw = max(fwidth(i.uv), 1e-6);
    if (range.w > 0 && !isNaNf(i.val)) {
        float d = abs(frac(f - 0.5) - 0.5) / fwf;
        c = lerp(c, c * 0.30, (1 - saturate(d - 0.4)) * 0.85);
    }
    if (flags.y > 0.5) {
        float2 a = smoothstep(float2(0, 0), fw * 1.2, i.uv) * smoothstep(float2(0, 0), fw * 1.2, 1 - i.uv);
        c = lerp(c, c * 0.3, (1 - min(a.x, a.y)) * 0.8);
    }
    return float4(c, range.z);
}

// ---- textured slices (seismic sections / time slices)
struct VSInT { float3 pos : POSITION; float2 uv : TEXCOORD0; };
struct PSInT { float4 pos : SV_POSITION; float2 uv : TEXCOORD0; };
PSInT VSTex(VSInT i) { PSInT o; o.pos = mul(float4(toWorld(i.pos), 1), VP); o.uv = i.uv; return o; }
float4 PSTex(PSInT i) : SV_TARGET {
    float v = stex.Sample(samp, i.uv);
    if (isNaNf(v)) discard;
    float t = saturate((v - range.x) * range.y);
    float3 c = cmap.SampleLevel(samp, t, 0).rgb;
    float a = range.z;
    if (flags.x > 1.5) a *= saturate((1 - t) * 1.6);
    else if (flags.x > 0.5) a *= saturate(abs(t - 0.5) * 3.0);
    if (a < 0.004) discard;
    return float4(c, a);
}

// ---- lit coloured solids and lines
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

// ---- ray-marched volume (full-screen triangle, premultiplied output)
float4 VSQuad(uint id : SV_VertexID) : SV_POSITION {
    float2 p = float2((id << 1) & 2, id & 2);
    return float4(p * float2(2, -2) + float2(-1, 1), 0, 1);
}
float zToS(float z) {
    float t;
    if (volP.w < 0.5) t = z / max(volZ.z, 1e-9);
    else {
        float dz = z - volV.z;
        t = abs(volV.y) > 1e-9 ? 2000.0 * log(max(1e-6, 1 + volV.y * dz / volV.x)) / volV.y : 2000.0 * dz / volV.x;
    }
    return (t - volZ.x) / volZ.y;
}
float4 PSVol(float4 pos : SV_POSITION) : SV_TARGET {
    float2 ndc = float2(pos.x * vpx.z * 2 - 1, 1 - pos.y * vpx.w * 2);
    float4 pn = mul(float4(ndc, 0, 1), invVP); pn /= pn.w;
    float4 pf = mul(float4(ndc, 1, 1), invVP); pf /= pf.w;
#if MSDEPTH
    float dz = sceneDepth.Load(int2(pos.xy), 0);
#else
    float dz = sceneDepth.Load(int3(pos.xy, 0));
#endif
    float3 ro = pn.xyz, rd = pf.xyz - pn.xyz;
    float tEnd = 1;
    if (dz < 1) { float4 ps = mul(float4(ndc, dz, 1), invVP); ps /= ps.w; tEnd = dot(ps.xyz - ro, rd) / dot(rd, rd); }
    float3 mo = float3(ro.x + center.x, ro.z + center.y, center.z - ro.y / center.w);
    float3 md = float3(rd.x, rd.z, -rd.y / center.w);
    float2 q = mo.xy - volP.xy;
    float i0 = volA.x * q.x + volA.y * q.y, di = volA.x * md.x + volA.y * md.y;
    float j0 = volA.z * q.x + volA.w * q.y, dj = volA.z * md.x + volA.w * md.y;
    float t0 = 0, t1 = tEnd;
    // slabs: i, j, z
    float3 o3 = float3(i0, j0, mo.z), d3 = float3(di, dj, md.z);
    [unroll] for (int a = 0; a < 3; ++a) {
        float lo = volB0[a], hi = volB1[a];
        if (abs(d3[a]) < 1e-12) { if (o3[a] < lo || o3[a] > hi) return float4(0, 0, 0, 0); }
        else {
            float ta = (lo - o3[a]) / d3[a], tb = (hi - o3[a]) / d3[a];
            t0 = max(t0, min(ta, tb)); t1 = min(t1, max(ta, tb));
        }
    }
    if (t1 <= t0) return float4(0, 0, 0, 0);
    int N = (int)volN.w;
    float dt = (t1 - t0) / N;
    float3 C = 0; float A = 0;
    float corr = 128.0 / max(volN.w, 1.0);
    // jitter the start (removes slicing artefacts)
    float jit = frac(sin(dot(pos.xy, float2(12.9898, 78.233))) * 43758.5453);
    [loop] for (int k = 0; k < N; ++k) {
        float t = t0 + (k + jit) * dt;
        float3 p = o3 + t * d3;
        float s = zToS(p.z);
        float3 uvw = float3((s + 0.5) / volN.z, (p.y + 0.5) / volN.y, (p.x + 0.5) / volN.x);
        float v = vol.SampleLevel(samp, uvw, 0);
        float tc = saturate((v - range.x) * range.y);
        float strength = abs(tc - 0.5) * 2;
        float al = saturate((strength - volP.z) / max(1 - volP.z, 1e-3));
        al = 1 - pow(saturate(1 - al * al * range.z), corr);
        if (al > 0.001) {
            float3 c = cmap.SampleLevel(samp, tc, 0).rgb;
            C += (1 - A) * al * c; A += (1 - A) * al;
            if (A > 0.985) break;
        }
    }
    return float4(C, A);
}
)";

struct FrameCB { float vp[16]; float inv[16]; float center[4]; float eye[4]; float light[4]; float vpx[4]; };

ID3DBlob* compileGeo(const char* entry, const char* target, bool msDepth) {
    ID3DBlob* code = nullptr; ID3DBlob* err = nullptr;
    const D3D_SHADER_MACRO defs[] = {{"MSDEPTH", msDepth ? "1" : "0"}, {nullptr, nullptr}};
    const HRESULT hr = D3DCompile(kGeoShader, std::strlen(kGeoShader), "geoscene.hlsl", defs, nullptr, entry, target, D3DCOMPILE_OPTIMIZATION_LEVEL3, 0, &code, &err);
    if (FAILED(hr)) {
        if (err) { OutputDebugStringA((const char*)err->GetBufferPointer()); MessageBoxA(nullptr, (const char*)err->GetBufferPointer(), "Geo scene shader error", MB_ICONERROR); err->Release(); }
        return nullptr;
    }
    SafeRelease(err);
    return code;
}
std::wstring widenG(const std::string& s) {
    if (s.empty()) return std::wstring();
    const int n = MultiByteToWideChar(CP_UTF8, 0, s.c_str(), int(s.size()), nullptr, 0);
    std::wstring w(size_t(n), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, s.c_str(), int(s.size()), &w[0], n);
    return w;
}
template <class V> ID3D11Buffer* makeBuffer(const std::vector<V>& v, UINT bind) {
    if (v.empty()) return nullptr;
    D3D11_BUFFER_DESC bd = {};
    bd.ByteWidth = UINT(v.size() * sizeof(V)); bd.Usage = D3D11_USAGE_IMMUTABLE; bd.BindFlags = bind;
    D3D11_SUBRESOURCE_DATA sd = {v.data(), 0, 0};
    ID3D11Buffer* b = nullptr;
    if (FAILED(g_gfx.dev->CreateBuffer(&bd, &sd, &b))) return nullptr;
    return b;
}
}  // namespace
#endif  // RESSIM_D3D11

namespace geo3d {
bool invert4(const float m[16], float out[16]) {
    double inv[16], det;
    const double* a = nullptr; double M[16];
    for (int k = 0; k < 16; ++k) M[k] = m[k];
    a = M;
    inv[0] = a[5] * a[10] * a[15] - a[5] * a[11] * a[14] - a[9] * a[6] * a[15] + a[9] * a[7] * a[14] + a[13] * a[6] * a[11] - a[13] * a[7] * a[10];
    inv[4] = -a[4] * a[10] * a[15] + a[4] * a[11] * a[14] + a[8] * a[6] * a[15] - a[8] * a[7] * a[14] - a[12] * a[6] * a[11] + a[12] * a[7] * a[10];
    inv[8] = a[4] * a[9] * a[15] - a[4] * a[11] * a[13] - a[8] * a[5] * a[15] + a[8] * a[7] * a[13] + a[12] * a[5] * a[11] - a[12] * a[7] * a[9];
    inv[12] = -a[4] * a[9] * a[14] + a[4] * a[10] * a[13] + a[8] * a[5] * a[14] - a[8] * a[6] * a[13] - a[12] * a[5] * a[10] + a[12] * a[6] * a[9];
    inv[1] = -a[1] * a[10] * a[15] + a[1] * a[11] * a[14] + a[9] * a[2] * a[15] - a[9] * a[3] * a[14] - a[13] * a[2] * a[11] + a[13] * a[3] * a[10];
    inv[5] = a[0] * a[10] * a[15] - a[0] * a[11] * a[14] - a[8] * a[2] * a[15] + a[8] * a[3] * a[14] + a[12] * a[2] * a[11] - a[12] * a[3] * a[10];
    inv[9] = -a[0] * a[9] * a[15] + a[0] * a[11] * a[13] + a[8] * a[1] * a[15] - a[8] * a[3] * a[13] - a[12] * a[1] * a[11] + a[12] * a[3] * a[9];
    inv[13] = a[0] * a[9] * a[14] - a[0] * a[10] * a[13] - a[8] * a[1] * a[14] + a[8] * a[2] * a[13] + a[12] * a[1] * a[10] - a[12] * a[2] * a[9];
    inv[2] = a[1] * a[6] * a[15] - a[1] * a[7] * a[14] - a[5] * a[2] * a[15] + a[5] * a[3] * a[14] + a[13] * a[2] * a[7] - a[13] * a[3] * a[6];
    inv[6] = -a[0] * a[6] * a[15] + a[0] * a[7] * a[14] + a[4] * a[2] * a[15] - a[4] * a[3] * a[14] - a[12] * a[2] * a[7] + a[12] * a[3] * a[6];
    inv[10] = a[0] * a[5] * a[15] - a[0] * a[7] * a[13] - a[4] * a[1] * a[15] + a[4] * a[3] * a[13] + a[12] * a[1] * a[7] - a[12] * a[3] * a[5];
    inv[14] = -a[0] * a[5] * a[14] + a[0] * a[6] * a[13] + a[4] * a[1] * a[14] - a[4] * a[2] * a[13] - a[12] * a[1] * a[6] + a[12] * a[2] * a[5];
    inv[3] = -a[1] * a[6] * a[11] + a[1] * a[7] * a[10] + a[5] * a[2] * a[11] - a[5] * a[3] * a[10] - a[9] * a[2] * a[7] + a[9] * a[3] * a[6];
    inv[7] = a[0] * a[6] * a[11] - a[0] * a[7] * a[10] - a[4] * a[2] * a[11] + a[4] * a[3] * a[10] + a[8] * a[2] * a[7] - a[8] * a[3] * a[6];
    inv[11] = -a[0] * a[5] * a[11] + a[0] * a[7] * a[9] + a[4] * a[1] * a[11] - a[4] * a[3] * a[9] - a[8] * a[1] * a[7] + a[8] * a[3] * a[5];
    inv[15] = a[0] * a[5] * a[10] - a[0] * a[6] * a[9] - a[4] * a[1] * a[10] + a[4] * a[2] * a[9] + a[8] * a[1] * a[6] - a[8] * a[2] * a[5];
    det = a[0] * inv[0] + a[1] * inv[4] + a[2] * inv[8] + a[3] * inv[12];
    if (std::fabs(det) < 1e-300) return false;
    for (int k = 0; k < 16; ++k) out[k] = float(inv[k] / det);
    return true;
}
}  // namespace geo3d

#ifdef RESSIM_D3D11

// ============================================================================ init / shutdown
bool GeoScene3D::init() {
    ID3D11Device* d = g_gfx.dev;
    UINT q = 0;
    if (FAILED(d->CheckMultisampleQualityLevels(DXGI_FORMAT_B8G8R8A8_UNORM, 4, &q)) || q == 0) msaa_ = 1;
    const bool ms = msaa_ > 1;
    ID3DBlob* b[9] = {compileGeo("VSMesh", "vs_5_0", ms), compileGeo("PSMesh", "ps_5_0", ms), compileGeo("VSTex", "vs_5_0", ms), compileGeo("PSTex", "ps_5_0", ms),
                      compileGeo("VSColor", "vs_5_0", ms), compileGeo("PSColor", "ps_5_0", ms), compileGeo("VSQuad", "vs_5_0", ms), compileGeo("PSVol", "ps_5_0", ms), nullptr};
    for (int k = 0; k < 8; ++k) if (!b[k]) { for (int m = 0; m < 8; ++m) SafeRelease(b[m]); return false; }
    d->CreateVertexShader(b[0]->GetBufferPointer(), b[0]->GetBufferSize(), nullptr, &vs_);
    d->CreatePixelShader(b[1]->GetBufferPointer(), b[1]->GetBufferSize(), nullptr, &psMesh_);
    d->CreateVertexShader(b[2]->GetBufferPointer(), b[2]->GetBufferSize(), nullptr, &vsT_);
    d->CreatePixelShader(b[3]->GetBufferPointer(), b[3]->GetBufferSize(), nullptr, &psT_);
    d->CreateVertexShader(b[4]->GetBufferPointer(), b[4]->GetBufferSize(), nullptr, &vsC_);
    d->CreatePixelShader(b[5]->GetBufferPointer(), b[5]->GetBufferSize(), nullptr, &psC_);
    d->CreateVertexShader(b[6]->GetBufferPointer(), b[6]->GetBufferSize(), nullptr, &vsQ_);
    d->CreatePixelShader(b[7]->GetBufferPointer(), b[7]->GetBufferSize(), nullptr, &psV_);
    const D3D11_INPUT_ELEMENT_DESC ilm[] = {
        {"POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0},
        {"NORMAL", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D11_INPUT_PER_VERTEX_DATA, 0},
        {"TEXCOORD", 0, DXGI_FORMAT_R32_FLOAT, 0, 24, D3D11_INPUT_PER_VERTEX_DATA, 0},
        {"TEXCOORD", 1, DXGI_FORMAT_R32G32_FLOAT, 0, 28, D3D11_INPUT_PER_VERTEX_DATA, 0}};
    d->CreateInputLayout(ilm, 4, b[0]->GetBufferPointer(), b[0]->GetBufferSize(), &il_);
    const D3D11_INPUT_ELEMENT_DESC ilt[] = {
        {"POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0},
        {"TEXCOORD", 0, DXGI_FORMAT_R32G32_FLOAT, 0, 12, D3D11_INPUT_PER_VERTEX_DATA, 0}};
    d->CreateInputLayout(ilt, 2, b[2]->GetBufferPointer(), b[2]->GetBufferSize(), &ilT_);
    const D3D11_INPUT_ELEMENT_DESC ilc[] = {
        {"POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0},
        {"NORMAL", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D11_INPUT_PER_VERTEX_DATA, 0},
        {"COLOR", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 0, 24, D3D11_INPUT_PER_VERTEX_DATA, 0}};
    d->CreateInputLayout(ilc, 3, b[4]->GetBufferPointer(), b[4]->GetBufferSize(), &ilC_);
    for (int k = 0; k < 8; ++k) SafeRelease(b[k]);

    D3D11_BUFFER_DESC bd = {};
    bd.Usage = D3D11_USAGE_DYNAMIC; bd.BindFlags = D3D11_BIND_CONSTANT_BUFFER; bd.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
    bd.ByteWidth = sizeof(FrameCB); d->CreateBuffer(&bd, nullptr, &cbFrame_);
    bd.ByteWidth = sizeof(ObjCB); d->CreateBuffer(&bd, nullptr, &cbObj_);
    D3D11_SAMPLER_DESC sd = {};
    sd.Filter = D3D11_FILTER_MIN_MAG_MIP_LINEAR;
    sd.AddressU = sd.AddressV = sd.AddressW = D3D11_TEXTURE_ADDRESS_CLAMP;
    sd.MaxLOD = D3D11_FLOAT32_MAX;
    d->CreateSamplerState(&sd, &samp_);
    sd.Filter = D3D11_FILTER_MIN_MAG_MIP_POINT;
    d->CreateSamplerState(&sd, &sampPt_);
    D3D11_RASTERIZER_DESC rd = {};
    rd.FillMode = D3D11_FILL_SOLID; rd.CullMode = D3D11_CULL_NONE; rd.DepthClipEnable = TRUE; rd.MultisampleEnable = TRUE; rd.AntialiasedLineEnable = TRUE;
    d->CreateRasterizerState(&rd, &rs_);
    D3D11_DEPTH_STENCIL_DESC dd = {};
    dd.DepthEnable = TRUE; dd.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ALL; dd.DepthFunc = D3D11_COMPARISON_LESS_EQUAL;
    d->CreateDepthStencilState(&dd, &dss_);
    dd.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ZERO;
    d->CreateDepthStencilState(&dd, &dssNoWrite_);
    D3D11_BLEND_DESC bl = {};
    bl.RenderTarget[0].BlendEnable = TRUE;
    bl.RenderTarget[0].SrcBlend = D3D11_BLEND_SRC_ALPHA; bl.RenderTarget[0].DestBlend = D3D11_BLEND_INV_SRC_ALPHA; bl.RenderTarget[0].BlendOp = D3D11_BLEND_OP_ADD;
    bl.RenderTarget[0].SrcBlendAlpha = D3D11_BLEND_ONE; bl.RenderTarget[0].DestBlendAlpha = D3D11_BLEND_INV_SRC_ALPHA; bl.RenderTarget[0].BlendOpAlpha = D3D11_BLEND_OP_ADD;
    bl.RenderTarget[0].RenderTargetWriteMask = D3D11_COLOR_WRITE_ENABLE_ALL;
    d->CreateBlendState(&bl, &bs_);
    bl.RenderTarget[0].SrcBlend = D3D11_BLEND_ONE;
    d->CreateBlendState(&bl, &bsPre_);
    const float s = g_gfx.dpiScale;
    g_gfx.dwrite->CreateTextFormat(L"Segoe UI", nullptr, DWRITE_FONT_WEIGHT_NORMAL, DWRITE_FONT_STYLE_NORMAL, DWRITE_FONT_STRETCH_NORMAL, 11.5f * s, L"en-us", &tfSmall_);
    g_gfx.dwrite->CreateTextFormat(L"Segoe UI", nullptr, DWRITE_FONT_WEIGHT_SEMI_BOLD, DWRITE_FONT_STYLE_NORMAL, DWRITE_FONT_STRETCH_NORMAL, 12.5f * s, L"en-us", &tfLabel_);
    g_gfx.dwrite->CreateTextFormat(L"Segoe UI", nullptr, DWRITE_FONT_WEIGHT_SEMI_BOLD, DWRITE_FONT_STYLE_NORMAL, DWRITE_FONT_STRETCH_NORMAL, 16.0f * s, L"en-us", &tfTitle_);
    ready_ = vs_ && psMesh_ && vsT_ && psT_ && vsC_ && psC_ && vsQ_ && psV_ && il_ && ilT_ && ilC_ && cbFrame_ && cbObj_;
    return ready_;
}

void GeoScene3D::releaseObjects() {
    for (Slice& s : slices_) { SafeRelease(s.vb); SafeRelease(s.ib); SafeRelease(s.srv); SafeRelease(s.tex); SafeRelease(s.cmSrv); SafeRelease(s.cm); }
    for (Mesh& m : meshes_) { SafeRelease(m.vb); SafeRelease(m.ib); SafeRelease(m.cmSrv); SafeRelease(m.cm); }
    for (Solid& s : solids_) { SafeRelease(s.vb); SafeRelease(s.ib); }
    slices_.clear(); meshes_.clear(); solids_.clear();
    SafeRelease(lvb_); lineCount_ = 0;
}
void GeoScene3D::releaseVolume() { SafeRelease(volSrv_); SafeRelease(volTex_); SafeRelease(volCmSrv_); SafeRelease(volCm_); hasVol_ = false; volCmapBuilt_ = Colormap(-1); }

void GeoScene3D::shutdown() {
    releaseObjects(); releaseVolume(); releaseTargets();
    SafeRelease(boxVb_);
    SafeRelease(vs_); SafeRelease(psMesh_); SafeRelease(il_); SafeRelease(vsT_); SafeRelease(psT_); SafeRelease(ilT_);
    SafeRelease(vsC_); SafeRelease(psC_); SafeRelease(ilC_); SafeRelease(vsQ_); SafeRelease(psV_);
    SafeRelease(cbFrame_); SafeRelease(cbObj_); SafeRelease(samp_); SafeRelease(sampPt_); SafeRelease(rs_);
    SafeRelease(dss_); SafeRelease(dssNoWrite_); SafeRelease(bs_); SafeRelease(bsPre_);
    SafeRelease(tfSmall_); SafeRelease(tfLabel_); SafeRelease(tfTitle_);
    ready_ = false;
}

void GeoScene3D::clear() { releaseObjects(); releaseVolume(); labels_.clear(); }

void GeoScene3D::releaseTargets() {
    if (g_gfx.d2d) g_gfx.d2d->SetTarget(nullptr);
    SafeRelease(d2dTarget_); SafeRelease(colorSRV_); SafeRelease(color_);
    SafeRelease(depthSRV_); SafeRelease(msDSV_); SafeRelease(msDepth_); SafeRelease(msRTV_); SafeRelease(msColor_);
    w_ = h_ = 0;
}

void GeoScene3D::createTargets(int w, int h) {
    releaseTargets();
    ID3D11Device* d = g_gfx.dev;
    D3D11_TEXTURE2D_DESC td = {};
    td.Width = UINT(w); td.Height = UINT(h); td.MipLevels = 1; td.ArraySize = 1; td.Format = DXGI_FORMAT_B8G8R8A8_UNORM;
    td.SampleDesc.Count = msaa_; td.Usage = D3D11_USAGE_DEFAULT; td.BindFlags = D3D11_BIND_RENDER_TARGET;
    d->CreateTexture2D(&td, nullptr, &msColor_);
    d->CreateRenderTargetView(msColor_, nullptr, &msRTV_);
    // depth readable by the volume pass (typeless storage: D32 for the depth test, R32 for sampling)
    td.Format = DXGI_FORMAT_R32_TYPELESS; td.BindFlags = D3D11_BIND_DEPTH_STENCIL | D3D11_BIND_SHADER_RESOURCE;
    d->CreateTexture2D(&td, nullptr, &msDepth_);
    D3D11_DEPTH_STENCIL_VIEW_DESC dv = {};
    dv.Format = DXGI_FORMAT_D32_FLOAT; dv.ViewDimension = msaa_ > 1 ? D3D11_DSV_DIMENSION_TEXTURE2DMS : D3D11_DSV_DIMENSION_TEXTURE2D;
    d->CreateDepthStencilView(msDepth_, &dv, &msDSV_);
    D3D11_SHADER_RESOURCE_VIEW_DESC sv = {};
    sv.Format = DXGI_FORMAT_R32_FLOAT;
    if (msaa_ > 1) sv.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE2DMS;
    else { sv.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE2D; sv.Texture2D.MipLevels = 1; }
    d->CreateShaderResourceView(msDepth_, &sv, &depthSRV_);
    td.Format = DXGI_FORMAT_B8G8R8A8_UNORM; td.SampleDesc.Count = 1; td.BindFlags = D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE;
    d->CreateTexture2D(&td, nullptr, &color_);
    d->CreateShaderResourceView(color_, nullptr, &colorSRV_);
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

ID3D11ShaderResourceView* GeoScene3D::makeColormap(Colormap cm, bool rev, ID3D11Texture1D** tex) {
    float px[256 * 4];
    for (int i = 0; i < 256; ++i) { float rgb[3]; colormapRGB(cm, i / 255.0f, rgb, rev); px[4 * i] = rgb[0]; px[4 * i + 1] = rgb[1]; px[4 * i + 2] = rgb[2]; px[4 * i + 3] = 1; }
    D3D11_TEXTURE1D_DESC td = {};
    td.Width = 256; td.MipLevels = 1; td.ArraySize = 1; td.Format = DXGI_FORMAT_R32G32B32A32_FLOAT; td.Usage = D3D11_USAGE_IMMUTABLE; td.BindFlags = D3D11_BIND_SHADER_RESOURCE;
    D3D11_SUBRESOURCE_DATA sd = {px, 0, 0};
    ID3D11ShaderResourceView* srv = nullptr;
    if (SUCCEEDED(g_gfx.dev->CreateTexture1D(&td, &sd, tex))) g_gfx.dev->CreateShaderResourceView(*tex, nullptr, &srv);
    return srv;
}

// ============================================================================ content
void GeoScene3D::setSlices(const std::vector<GeoSliceDesc>& in) {
    for (Slice& s : slices_) { SafeRelease(s.vb); SafeRelease(s.ib); SafeRelease(s.srv); SafeRelease(s.tex); SafeRelease(s.cmSrv); SafeRelease(s.cm); }
    slices_.clear();
    for (const GeoSliceDesc& D : in) {
        Slice S;
        std::vector<VtxT> v; std::vector<unsigned> idx;
        if (!sliceGeometry(D, v, idx)) continue;
        S.vb = makeBuffer(v, D3D11_BIND_VERTEX_BUFFER); S.ib = makeBuffer(idx, D3D11_BIND_INDEX_BUFFER); S.n = unsigned(idx.size());
        D3D11_TEXTURE2D_DESC td = {};
        td.Width = UINT(D.tw); td.Height = UINT(D.th); td.MipLevels = 1; td.ArraySize = 1; td.Format = DXGI_FORMAT_R32_FLOAT;
        td.SampleDesc.Count = 1; td.Usage = D3D11_USAGE_IMMUTABLE; td.BindFlags = D3D11_BIND_SHADER_RESOURCE;
        D3D11_SUBRESOURCE_DATA sd = {D.tex.data(), UINT(D.tw * 4), 0};
        if (SUCCEEDED(g_gfx.dev->CreateTexture2D(&td, &sd, &S.tex))) g_gfx.dev->CreateShaderResourceView(S.tex, nullptr, &S.srv);
        S.cmSrv = makeColormap(D.cmap, D.reverse, &S.cm);
        S.lo = D.lo; S.hi = D.hi; S.opacity = D.opacity; S.alphaMode = D.alphaMode;
        slices_.push_back(S);
    }
}

void GeoScene3D::setMeshes(const std::vector<GeoMeshDesc>& in) {
    for (Mesh& m : meshes_) { SafeRelease(m.vb); SafeRelease(m.ib); SafeRelease(m.cmSrv); SafeRelease(m.cm); }
    meshes_.clear();
    for (const GeoMeshDesc& D : in) {
        std::vector<VtxM> v; float cz[3];
        if (!meshGeometry(D, v, cz)) continue;
        const size_t n = v.size();
        Mesh M;
        M.vb = makeBuffer(v, D3D11_BIND_VERTEX_BUFFER); M.ib = makeBuffer(D.idx, D3D11_BIND_INDEX_BUFFER); M.n = unsigned(D.idx.size());
        M.cmSrv = makeColormap(D.cmap, D.reverse, &M.cm);
        M.lo = D.lo; M.hi = D.hi; M.opacity = D.opacity; M.contour = D.contour; M.edges = D.edges && D.uv.size() >= 2 * n;
        for (int m = 0; m < 3; ++m) M.cz[m] = cz[m];
        meshes_.push_back(M);
    }
}

void GeoScene3D::setSolids(const std::vector<GeoSolidDesc>& in) {
    for (Solid& s : solids_) { SafeRelease(s.vb); SafeRelease(s.ib); }
    solids_.clear();
    for (const GeoSolidDesc& D : in) {
        std::vector<VtxC> v;
        if (!solidGeometry(D, v)) continue;
        Solid S;
        S.vb = makeBuffer(v, D3D11_BIND_VERTEX_BUFFER); S.ib = makeBuffer(D.idx, D3D11_BIND_INDEX_BUFFER); S.n = unsigned(D.idx.size()); S.transparent = D.transparent;
        solids_.push_back(S);
    }
}

void GeoScene3D::setLines(const std::vector<GeoLineVtx>& l) {
    SafeRelease(lvb_); lineCount_ = 0;
    if (l.size() < 2) return;
    const std::vector<VtxC> v = lineVertices(l);
    lvb_ = makeBuffer(v, D3D11_BIND_VERTEX_BUFFER);
    if (lvb_) lineCount_ = unsigned(v.size() & ~size_t(1));
}

void GeoScene3D::setVolume(const GeoVolumeDesc* v) {
    releaseVolume();
    if (!v || v->ni < 2 || v->nj < 2 || v->ns < 2 || v->data.size() < size_t(v->ni) * size_t(v->nj) * size_t(v->ns)) return;
    D3D11_TEXTURE3D_DESC td = {};
    td.Width = UINT(v->ns); td.Height = UINT(v->nj); td.Depth = UINT(v->ni); td.MipLevels = 1; td.Format = DXGI_FORMAT_R32_FLOAT;
    td.Usage = D3D11_USAGE_IMMUTABLE; td.BindFlags = D3D11_BIND_SHADER_RESOURCE;
    D3D11_SUBRESOURCE_DATA sd = {v->data.data(), UINT(v->ns * 4), UINT(v->ns * v->nj * 4)};
    if (FAILED(g_gfx.dev->CreateTexture3D(&td, &sd, &volTex_))) return;
    g_gfx.dev->CreateShaderResourceView(volTex_, nullptr, &volSrv_);
    vol_ = *v;
    vol_.data.clear(); vol_.data.shrink_to_fit();   // the GPU keeps the samples
    hasVol_ = volSrv_ != nullptr;
}

#endif  // RESSIM_D3D11

void GeoScene3D::setVolumeStyle(const GeoVolumeDesc& v) {
    if (!hasVol_) return;
    vol_.lo = v.lo; vol_.hi = v.hi; vol_.cmap = v.cmap; vol_.reverse = v.reverse; vol_.opacity = v.opacity; vol_.threshold = v.threshold; vol_.steps = v.steps;
    vol_.iMin = v.iMin; vol_.iMax = v.iMax; vol_.jMin = v.jMin; vol_.jMax = v.jMax; vol_.sMin = v.sMin; vol_.sMax = v.sMax;
    vol_.zMode = v.zMode; vol_.zScale = v.zScale; vol_.v0 = v.v0; vol_.k = v.k; vol_.datum = v.datum;
}

void GeoScene3D::setBounds(const float lo[3], const float hi[3]) {
    for (int k = 0; k < 3; ++k) { lo_[k] = std::min(lo[k], hi[k]); hi_[k] = std::max(lo[k], hi[k]); center_[k] = 0.5f * (lo_[k] + hi_[k]); }
    builtBoxZ_ = -1;
}

// ============================================================================ camera
void GeoScene3D::worldOf(float x, float y, float z, float o[3]) const { o[0] = x - center_[0]; o[1] = -(z - center_[2]) * zScale; o[2] = y - center_[1]; }
void GeoScene3D::eyePos(float e[3]) const {
    const float d = dist_ * extent_;
    e[0] = panX_ + d * std::cos(pitch_) * std::sin(yaw_); e[1] = panY_ + d * std::sin(pitch_); e[2] = panZ_ + d * std::cos(pitch_) * std::cos(yaw_);
}
std::array<float, 16> GeoScene3D::viewProj() const {
    const float d = dist_ * extent_;
    const float tgt[3] = {panX_, panY_, panZ_};
    float eye[3]; eyePos(eye);
    float z[3] = {eye[0] - tgt[0], eye[1] - tgt[1], eye[2] - tgt[2]}; norm3(z);
    const float up[3] = {0, 1, 0};
    float x[3]; cross(up, z, x); norm3(x);
    if (x[0] == 0 && x[1] == 0 && x[2] == 0) { x[0] = 1; }
    float y[3]; cross(z, x, y);
    const float V[16] = {x[0], x[1], x[2], -(x[0] * eye[0] + x[1] * eye[1] + x[2] * eye[2]),
                         y[0], y[1], y[2], -(y[0] * eye[0] + y[1] * eye[1] + y[2] * eye[2]),
                         z[0], z[1], z[2], -(z[0] * eye[0] + z[1] * eye[1] + z[2] * eye[2]), 0, 0, 0, 1};
    const float aspect = h_ > 0 ? float(w_) / float(h_) : 1.0f;
    const float f = 1.0f / std::tan(0.5f * 0.75f);
    const float zn = 0.01f * d, zf = d + 4.0f * extent_;
    const float P[16] = {f / aspect, 0, 0, 0, 0, f, 0, 0, 0, 0, zf / (zn - zf), zn * zf / (zn - zf), 0, 0, -1, 0};
    std::array<float, 16> M{};
    for (int r = 0; r < 4; ++r) for (int c = 0; c < 4; ++c) { float s = 0; for (int k = 0; k < 4; ++k) s += P[r * 4 + k] * V[k * 4 + c]; M[size_t(r * 4 + c)] = s; }
    return M;
}
bool GeoScene3D::project(float x, float y, float z, float& sx, float& sy) const {
    float w[3]; worldOf(x, y, z, w);
    const auto M = viewProj();
    float c[4];
    for (int r = 0; r < 4; ++r) c[r] = M[size_t(r * 4)] * w[0] + M[size_t(r * 4 + 1)] * w[1] + M[size_t(r * 4 + 2)] * w[2] + M[size_t(r * 4 + 3)];
    if (c[3] <= 1e-6f) return false;
    sx = (c[0] / c[3] * 0.5f + 0.5f) * float(w_);
    sy = (1.0f - (c[1] / c[3] * 0.5f + 0.5f)) * float(h_);
    return true;
}
bool GeoScene3D::pixelRay(float px, float py, double o[3], double d[3]) const {
    if (w_ <= 0 || h_ <= 0) return false;
    const auto M = viewProj();
    float inv[16];
    if (!invert4(M.data(), inv)) return false;
    const float nx = px / float(w_) * 2 - 1, ny = 1 - py / float(h_) * 2;
    auto un = [&](float z, double out[3]) {
        double c[4];
        for (int r = 0; r < 4; ++r) c[r] = double(inv[r * 4]) * nx + double(inv[r * 4 + 1]) * ny + double(inv[r * 4 + 2]) * z + double(inv[r * 4 + 3]);
        for (int k = 0; k < 3; ++k) out[k] = c[k] / c[3];
    };
    double a[3], b[3];
    un(0.0f, a); un(1.0f, b);
    // world -> model: x = wx + cx, y = wz + cy, z = cz - wy / zScale
    o[0] = a[0] + center_[0]; o[1] = a[2] + center_[1]; o[2] = center_[2] - a[1] / zScale;
    double e[3] = {b[0] + center_[0], b[2] + center_[1], center_[2] - b[1] / zScale};
    double l = 0;
    for (int k = 0; k < 3; ++k) { d[k] = e[k] - o[k]; l += d[k] * d[k]; }
    l = std::sqrt(l);
    if (l <= 0) return false;
    for (int k = 0; k < 3; ++k) d[k] /= l;
    return true;
}
void GeoScene3D::orbit(float dx, float dy) { yaw_ -= dx * 0.008f; pitch_ = std::min(1.55f, std::max(-1.55f, pitch_ + dy * 0.008f)); }
void GeoScene3D::pan(float dx, float dy) {
    const float d = dist_ * extent_ * 0.0015f;
    panX_ -= dx * d * std::cos(yaw_); panZ_ += dx * d * std::sin(yaw_);
    panY_ += dy * d * std::cos(pitch_);
}
void GeoScene3D::zoom(float wheel) { dist_ = std::min(20.0f, std::max(0.03f, dist_ * std::pow(0.88f, wheel))); }
void GeoScene3D::resetCamera() { yaw_ = 0.7f; pitch_ = 0.5f; dist_ = 2.0f; panX_ = panY_ = panZ_ = 0; }
void GeoScene3D::viewPreset(int which) {
    panX_ = panY_ = panZ_ = 0; dist_ = 1.9f;
    switch (which) {
        case 0: yaw_ = 0; pitch_ = 1.55f; break;
        case 1: yaw_ = 0; pitch_ = 0.0f; break;
        case 2: yaw_ = 1.5708f; pitch_ = 0.0f; break;
        default: yaw_ = 0.7f; pitch_ = 0.5f; break;
    }
}

// ============================================================================ render
#ifdef RESSIM_D3D11
ImTextureID GeoScene3D::render(int w, int h) {
    w = std::max(16, w); h = std::max(16, h);
    if (w != w_ || h != h_ || !colorSRV_) createTargets(w, h);
    if (!msRTV_ || !msDSV_) return ImTextureID(0);
    extent_ = std::max({hi_[0] - lo_[0], hi_[1] - lo_[1], (hi_[2] - lo_[2]) * zScale, 1.0f});
    if (builtBoxZ_ != zScale || !boxVb_) {   // bounding box lines with tick marks
        SafeRelease(boxVb_); boxCount_ = 0;
        const std::vector<VtxC> b = boxLines(lo_, hi_);
        boxVb_ = makeBuffer(b, D3D11_BIND_VERTEX_BUFFER); boxCount_ = unsigned(b.size());
        builtBoxZ_ = zScale;
    }
    ID3D11DeviceContext* c = g_gfx.ctx;
    const float bg[4] = {background[0], background[1], background[2], 1.0f};
    c->ClearRenderTargetView(msRTV_, bg);
    c->ClearDepthStencilView(msDSV_, D3D11_CLEAR_DEPTH, 1.0f, 0);
    c->OMSetRenderTargets(1, &msRTV_, msDSV_);
    D3D11_VIEWPORT vp = {0, 0, float(w), float(h), 0, 1};
    c->RSSetViewports(1, &vp);
    c->RSSetState(rs_);
    const float bf[4] = {0, 0, 0, 0};
    c->OMSetBlendState(bs_, bf, 0xffffffff);
    c->OMSetDepthStencilState(dss_, 0);
    FrameCB F = {};
    const auto M = viewProj();
    std::memcpy(F.vp, M.data(), sizeof(F.vp));
    if (!invert4(M.data(), F.inv)) std::memset(F.inv, 0, sizeof(F.inv));
    F.center[0] = center_[0]; F.center[1] = center_[1]; F.center[2] = center_[2]; F.center[3] = zScale;
    eyePos(F.eye);
    F.light[0] = 0.4f; F.light[1] = 0.85f; F.light[2] = 0.35f;
    F.vpx[0] = float(w); F.vpx[1] = float(h); F.vpx[2] = 1.0f / float(w); F.vpx[3] = 1.0f / float(h);
    D3D11_MAPPED_SUBRESOURCE ms;
    if (SUCCEEDED(c->Map(cbFrame_, 0, D3D11_MAP_WRITE_DISCARD, 0, &ms))) { std::memcpy(ms.pData, &F, sizeof(F)); c->Unmap(cbFrame_, 0); }
    auto setObj = [&](const ObjCB& o) { if (SUCCEEDED(c->Map(cbObj_, 0, D3D11_MAP_WRITE_DISCARD, 0, &ms))) { std::memcpy(ms.pData, &o, sizeof(o)); c->Unmap(cbObj_, 0); } };
    ID3D11Buffer* cbs[2] = {cbFrame_, cbObj_};
    c->VSSetConstantBuffers(0, 2, cbs);
    c->PSSetConstantBuffers(0, 2, cbs);
    c->PSSetSamplers(0, 1, &samp_);

    auto drawMesh = [&](const Mesh& m) {
        ObjCB o = {};
        o.range[0] = m.lo; o.range[1] = m.hi > m.lo ? 1.0f / (m.hi - m.lo) : 0.0f; o.range[2] = m.opacity; o.range[3] = m.contour; o.flags[1] = m.edges ? 1.0f : 0.0f;
        setObj(o);
        UINT st = sizeof(VtxM), of = 0;
        c->IASetInputLayout(il_); c->IASetVertexBuffers(0, 1, &m.vb, &st, &of); c->IASetIndexBuffer(m.ib, DXGI_FORMAT_R32_UINT, 0);
        c->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
        c->VSSetShader(vs_, nullptr, 0); c->PSSetShader(psMesh_, nullptr, 0);
        c->PSSetShaderResources(1, 1, &m.cmSrv);
        c->DrawIndexed(m.n, 0, 0);
    };
    auto drawSlice = [&](const Slice& s) {
        ObjCB o = {};
        o.range[0] = s.lo; o.range[1] = s.hi > s.lo ? 1.0f / (s.hi - s.lo) : 0.0f; o.range[2] = s.opacity; o.flags[0] = float(s.alphaMode);
        setObj(o);
        UINT st = sizeof(VtxT), of = 0;
        c->IASetInputLayout(ilT_); c->IASetVertexBuffers(0, 1, &s.vb, &st, &of); c->IASetIndexBuffer(s.ib, DXGI_FORMAT_R32_UINT, 0);
        c->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
        c->VSSetShader(vsT_, nullptr, 0); c->PSSetShader(psT_, nullptr, 0);
        ID3D11ShaderResourceView* r[2] = {s.srv, s.cmSrv};
        c->PSSetShaderResources(0, 2, r);
        c->DrawIndexed(s.n, 0, 0);
    };
    auto drawSolid = [&](const Solid& s) {
        UINT st = sizeof(VtxC), of = 0;
        c->IASetInputLayout(ilC_); c->IASetVertexBuffers(0, 1, &s.vb, &st, &of); c->IASetIndexBuffer(s.ib, DXGI_FORMAT_R32_UINT, 0);
        c->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
        c->VSSetShader(vsC_, nullptr, 0); c->PSSetShader(psC_, nullptr, 0);
        c->DrawIndexed(s.n, 0, 0);
    };
    auto drawLines = [&](ID3D11Buffer* vb, unsigned n) {
        if (!vb || n < 2) return;
        UINT st = sizeof(VtxC), of = 0;
        c->IASetInputLayout(ilC_); c->IASetVertexBuffers(0, 1, &vb, &st, &of);
        c->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_LINELIST);
        c->VSSetShader(vsC_, nullptr, 0); c->PSSetShader(psC_, nullptr, 0);
        c->Draw(n, 0);
    };
    // 1) opaque geometry
    for (const Mesh& m : meshes_) if (m.opacity >= 0.999f && m.vb) drawMesh(m);
    for (const Solid& s : solids_) if (!s.transparent && s.vb) drawSolid(s);
    for (const Slice& s : slices_) if (s.opacity >= 0.999f && s.alphaMode == 0 && s.vb && s.srv) drawSlice(s);
    drawLines(lvb_, lineCount_);
    if (showBox) drawLines(boxVb_, boxCount_);
    // 2) volume: ray-marched against the opaque depth
    if (hasVol_ && volSrv_ && depthSRV_) {
        if (volCmapBuilt_ != vol_.cmap || volRevBuilt_ != vol_.reverse || !volCmSrv_) {
            SafeRelease(volCmSrv_); SafeRelease(volCm_);
            volCmSrv_ = makeColormap(vol_.cmap, vol_.reverse, &volCm_);
            volCmapBuilt_ = vol_.cmap; volRevBuilt_ = vol_.reverse;
        }
        const ObjCB o = volumeConstants(vol_);
        setObj(o);
        c->OMSetRenderTargets(1, &msRTV_, nullptr);   // the depth buffer is read, not tested
        c->OMSetBlendState(bsPre_, bf, 0xffffffff);
        c->IASetInputLayout(nullptr);
        c->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
        c->VSSetShader(vsQ_, nullptr, 0); c->PSSetShader(psV_, nullptr, 0);
        ID3D11ShaderResourceView* r[4] = {nullptr, volCmSrv_, volSrv_, depthSRV_};
        c->PSSetShaderResources(0, 4, r);
        c->Draw(3, 0);
        ID3D11ShaderResourceView* nulls[4] = {nullptr, nullptr, nullptr, nullptr};
        c->PSSetShaderResources(0, 4, nulls);
        c->OMSetBlendState(bs_, bf, 0xffffffff);
        c->OMSetRenderTargets(1, &msRTV_, msDSV_);
    }
    // 3) transparent geometry (depth tested, not written), far objects first
    c->OMSetDepthStencilState(dssNoWrite_, 0);
    float eye[3]; eyePos(eye);
    std::vector<std::pair<float, int>> order;
    for (int k = 0; k < int(meshes_.size()); ++k) if (meshes_[size_t(k)].opacity < 0.999f && meshes_[size_t(k)].vb) {
        float wc[3]; worldOf(meshes_[size_t(k)].cz[0], meshes_[size_t(k)].cz[1], meshes_[size_t(k)].cz[2], wc);
        order.push_back({-((wc[0] - eye[0]) * (wc[0] - eye[0]) + (wc[1] - eye[1]) * (wc[1] - eye[1]) + (wc[2] - eye[2]) * (wc[2] - eye[2])), k});
    }
    std::sort(order.begin(), order.end());
    for (const auto& pr : order) drawMesh(meshes_[size_t(pr.second)]);
    for (const Slice& s : slices_) if ((s.opacity < 0.999f || s.alphaMode != 0) && s.vb && s.srv) drawSlice(s);
    for (const Solid& s : solids_) if (s.transparent && s.vb) drawSolid(s);
    c->OMSetDepthStencilState(dss_, 0);
    ID3D11ShaderResourceView* nulls[4] = {nullptr, nullptr, nullptr, nullptr};
    c->PSSetShaderResources(0, 4, nulls);
    if (msaa_ > 1) c->ResolveSubresource(color_, 0, msColor_, 0, DXGI_FORMAT_B8G8R8A8_UNORM);
    else c->CopyResource(color_, msColor_);
    ID3D11RenderTargetView* nullRTV = nullptr;
    c->OMSetRenderTargets(1, &nullRTV, nullptr);
    drawOverlay();
    return (ImTextureID)(intptr_t)colorSRV_;
}

// ============================================================================ overlay
void GeoScene3D::drawOverlay() {
    if (!d2dTarget_) return;
    ID2D1DeviceContext* d = g_gfx.d2d;
    d->SetTarget(d2dTarget_);
    d->BeginDraw();
    d->SetTransform(D2D1::Matrix3x2F::Scale(overlayScale_, overlayScale_));
    {
        D2DPainter p(d, g_gfx.dwrite, tfSmall_, tfLabel_, tfTitle_);
        if (p.ok()) paintOverlay(p);
    }
    d->EndDraw();
    d->SetTarget(nullptr);
}
#endif  // RESSIM_D3D11

// overlay layout (all back ends)
void GeoScene3D::paintOverlay(OverlayPainter& P) {
    const float S = overlayScale_;
    const float W = float(w_) / S, H = float(h_) / S;
    const OvColor text{0.93f, 0.94f, 0.96f, 1.0f}, panel{0.06f, 0.065f, 0.075f, 0.74f};
    const float s = g_gfx.dpiScale;
    if (showLabels) {
        for (const GeoLabel& L : labels_) {
            float sx, sy;
            if (!project(L.p[0], L.p[1], L.p[2], sx, sy)) continue;
            sx /= S; sy /= S;
            if (sx < 0 || sy < 0 || sx > W || sy > H) continue;
            const OvFont f = L.small ? OvFont::Small : OvFont::Label;
            const float tw = P.textWidth(L.text, f);
            const float hgt = (L.small ? 17 : 21) * s, lift = (L.small ? 6 : 9) * s;
            const float rl = sx - tw * 0.5f - 6 * s, rt = sy - lift - hgt, rr = sx + tw * 0.5f + 6 * s, rb = sy - lift;
            const OvColor c{L.rgba[0], L.rgba[1], L.rgba[2], L.rgba[3]};
            P.fillRoundRect(rl, rt, rr, rb, 4 * s, c);
            P.line(sx, sy - lift, sx, sy, c, 1.5f * s);
            P.text(L.text, f, rl + 6 * s, rt + 1.5f * s, text);
        }
    }
    if (!title_.empty()) {
        const float tw = std::max(P.textWidth(title_, OvFont::Title), P.textWidth(caption_, OvFont::Small));
        P.fillRoundRect(10 * s, 10 * s, 10 * s + tw + 24 * s, 10 * s + (caption_.empty() ? 32 : 52) * s, 8 * s, panel);
        P.text(title_, OvFont::Title, 22 * s, 14 * s, text);
        if (!caption_.empty()) P.text(caption_, OvFont::Small, 22 * s, 37 * s, OvColor{0.65f, 0.7f, 0.78f, 1});
    }
    if (showLegend && !legTitle_.empty()) {
        const float lh = std::min(H * 0.5f, 300 * s), lw = 15 * s;
        const float x0 = W - 92 * s, y0 = 0.5f * (H - lh);
        P.fillRoundRect(x0 - 12 * s, y0 - 34 * s, W - 8 * s, y0 + lh + 14 * s, 8 * s, panel);
        P.text(legTitle_, OvFont::Small, x0 - 4 * s, y0 - 28 * s, text);
        for (int i = 0; i < 96; ++i) {
            float rgb[3]; colormapRGB(legCmap_, 1.0f - (i + 0.5f) / 96.0f, rgb, legRev_);
            P.fillRect(x0, y0 + lh * i / 96.0f, x0 + lw, y0 + lh * (i + 1) / 96.0f + 0.8f, OvColor{rgb[0], rgb[1], rgb[2], 1});
        }
        P.strokeRect(x0, y0, x0 + lw, y0 + lh, OvColor{1, 1, 1, 0.35f}, 1.0f);
        for (int t = 0; t < 6; ++t) {
            const float f = float(t) / 5.0f, v = legLo_ + f * (legHi_ - legLo_);
            const float y = y0 + lh * (1.0f - f);
            P.line(x0 + lw, y, x0 + lw + 4 * s, y, OvColor{1, 1, 1, 0.6f}, 1.0f);
            P.text(ovLegendNumber(v), OvFont::Small, x0 + lw + 7 * s, y - 8 * s, text);
        }
    }
    if (showAxes) {
        const float ox = 55 * s, oy = H - 55 * s, L = 38 * s;
        float z[3] = {std::cos(pitch_) * std::sin(yaw_), std::sin(pitch_), std::cos(pitch_) * std::cos(yaw_)};
        const float up[3] = {0, 1, 0};
        float x[3]; cross(up, z, x); norm3(x);
        float y[3]; cross(z, x, y);
        struct Ax { float w[3]; const char* name; OvColor col; };
        const Ax axes[3] = {{{1, 0, 0}, "E", OvColor{0.95f, 0.36f, 0.36f, 1}}, {{0, 0, 1}, "N", OvColor{0.36f, 0.86f, 0.46f, 1}},
                            {{0, -1, 0}, "Z", OvColor{0.40f, 0.62f, 1.0f, 1}}};
        P.fillCircle(ox, oy, 46 * s, panel);
        for (const auto& a : axes) {
            const float vx = x[0] * a.w[0] + x[1] * a.w[1] + x[2] * a.w[2], vy = y[0] * a.w[0] + y[1] * a.w[1] + y[2] * a.w[2];
            const float ex = ox + vx * L, ey = oy - vy * L;
            P.line(ox, oy, ex, ey, a.col, 2.2f * s);
            P.text(a.name, OvFont::Small, ex + 3 * s, ey - 8 * s, a.col);
        }
        char zs[32]; std::snprintf(zs, sizeof(zs), "VE x%.1f", zScale);
        P.text(zs, OvFont::Small, ox - 22 * s, oy + 48 * s, OvColor{0.7f, 0.74f, 0.8f, 1});
    }
}

#ifdef RESSIM_D3D11
// ============================================================================ snapshot
bool GeoScene3D::saveImage(const std::string& path, int scale, std::string* err) {
    auto fail = [&](const char* what) { if (err) *err = what; return false; };
    const int w0 = w_, h0 = h_;
    if (w0 <= 0 || h0 <= 0) return fail("view not rendered yet");
    scale = std::max(1, std::min(4, scale));
    overlayScale_ = float(scale);
    render(w0 * scale, h0 * scale);
    overlayScale_ = 1.0f;
    if (!color_) { render(w0, h0); return fail("render target"); }
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
        IWICImagingFactory* wic = nullptr; IWICBitmap* bmp = nullptr; IWICStream* stream = nullptr; IWICBitmapEncoder* enc = nullptr; IWICBitmapFrameEncode* frame = nullptr;
        HRESULT hr = CoCreateInstance(CLSID_WICImagingFactory, nullptr, CLSCTX_INPROC_SERVER, IID_IWICImagingFactory, (void**)&wic);
        if (SUCCEEDED(hr)) hr = wic->CreateBitmapFromMemory(W, H, GUID_WICPixelFormat32bppBGRA, W * 4, UINT(px.size()), px.data(), &bmp);
        const std::wstring wp = widenG(path);
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
#endif  // RESSIM_D3D11
