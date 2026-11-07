// C16: sokol_gfx compute backend for the engine's GPU block-solver interface (see gpu_compute.h)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <d3d11.h>
#include <dxgi.h>
#include "gpu_compute.h"
#include "ressim/bsr.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <random>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#define SOKOL_IMPL
#define SOKOL_D3D11
#if !defined(NDEBUG) && !defined(RESSIM_SOKOL_VALIDATE)
#define NDEBUG_SOKOL_LOCAL
#define NDEBUG
#endif
#include "../../third_party/sokol/sokol_gfx.h"
#ifdef NDEBUG_SOKOL_LOCAL
#undef NDEBUG
#endif

namespace gpu {
namespace {

std::mutex g_mx;
bool g_tried = false, g_ok = false;
std::string g_info, g_err, g_lastLog;
ID3D11Device* g_dev = nullptr;
ID3D11DeviceContext* g_ctx = nullptr;

void sgLog(const char*, uint32_t level, uint32_t item, const char* msg, uint32_t line, const char*, void*) {
    // keep errors and the shader-compiler output (logged as a message after the failure) for the diagnosis line
    if (level <= 1 || (msg && std::strlen(msg) > 0)) {
        char b[512]; std::snprintf(b, sizeof(b), "sokol_gfx: %s (item %u, line %u)", msg ? msg : "error", item, line);
        if (g_lastLog.size() < 1500) g_lastLog += (g_lastLog.empty() ? "" : " | ") + std::string(b);
    }
}

// ---------------------------------------------------------------- HLSL compute kernels (cs_5_0, raw buffers)
const char* kCommon = R"(
cbuffer P : register(b0) { uint n; uint nbr; uint nb; uint ng; float ca; float cb; float cc; float cd; };
uint gidx(uint3 g, uint3 l, uint w) { return (g.y * 65535u + g.x) * w + l.x; }
)";
const char* kSpmv = R"(
ByteAddressBuffer rowPtr : register(t0);
ByteAddressBuffer colIdx : register(t1);
ByteAddressBuffer vals   : register(t2);
ByteAddressBuffer xv     : register(t3);
RWByteAddressBuffer yv   : register(u0);
[numthreads(64, 1, 1)]
void main(uint3 g : SV_GroupID, uint3 l : SV_GroupThreadID) {
    uint i = gidx(g, l, 64u);
    if (i >= nbr) return;
    float acc[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    uint k0 = rowPtr.Load(i * 4), k1 = rowPtr.Load(i * 4 + 4);
    for (uint k = k0; k < k1; ++k) {
        uint j = colIdx.Load(k * 4);
        uint base = k * nb * nb;
        for (uint e = 0; e < nb; ++e) {
            float s = 0;
            for (uint f = 0; f < nb; ++f) s += asfloat(vals.Load((base + e * nb + f) * 4)) * asfloat(xv.Load((j * nb + f) * 4));
            acc[e] += s;
        }
    }
    for (uint e2 = 0; e2 < nb; ++e2) yv.Store((i * nb + e2) * 4, asuint(acc[e2]));
}
)";
const char* kPrec = R"(
ByteAddressBuffer dinv : register(t0);
ByteAddressBuffer xv   : register(t1);
RWByteAddressBuffer zv : register(u0);
[numthreads(64, 1, 1)]
void main(uint3 g : SV_GroupID, uint3 l : SV_GroupThreadID) {
    uint i = gidx(g, l, 64u);
    if (i >= nbr) return;
    uint base = i * nb * nb;
    for (uint e = 0; e < nb; ++e) {
        float s = 0;
        for (uint f = 0; f < nb; ++f) s += asfloat(dinv.Load((base + e * nb + f) * 4)) * asfloat(xv.Load((i * nb + f) * 4));
        zv.Store((i * nb + e) * 4, asuint(s));
    }
}
)";
const char* kAxpbypcz = R"(
ByteAddressBuffer X : register(t0);
ByteAddressBuffer Y : register(t1);
ByteAddressBuffer Z : register(t2);
RWByteAddressBuffer O : register(u0);
[numthreads(256, 1, 1)]
void main(uint3 g : SV_GroupID, uint3 l : SV_GroupThreadID) {
    uint i = gidx(g, l, 256u);
    if (i >= n) return;
    float v = ca * asfloat(X.Load(i * 4)) + cb * asfloat(Y.Load(i * 4)) + cc * asfloat(Z.Load(i * 4));
    O.Store(i * 4, asuint(v));
}
)";
// two dot products sharing X: partial sums per group of 256 threads x 4 elements
const char* kDot2 = R"(
ByteAddressBuffer X : register(t0);
ByteAddressBuffer Y : register(t1);
ByteAddressBuffer Z : register(t2);
RWByteAddressBuffer O : register(u0);
groupshared float sa[256];
groupshared float sb[256];
[numthreads(256, 1, 1)]
void main(uint3 g : SV_GroupID, uint3 l : SV_GroupThreadID) {
    uint grp = g.y * 65535u + g.x;
    float a = 0, b = 0;
    for (uint q = 0; q < 4; ++q) {
        uint i = grp * 1024u + q * 256u + l.x;
        if (i < n) { float x = asfloat(X.Load(i * 4)); a += x * asfloat(Y.Load(i * 4)); b += x * asfloat(Z.Load(i * 4)); }
    }
    sa[l.x] = a; sb[l.x] = b;
    GroupMemoryBarrierWithGroupSync();
    for (uint s = 128; s > 0; s >>= 1) {
        if (l.x < s) { sa[l.x] += sa[l.x + s]; sb[l.x] += sb[l.x + s]; }
        GroupMemoryBarrierWithGroupSync();
    }
    if (l.x == 0) { O.Store(grp * 4, asuint(sa[0])); O.Store((ng + grp) * 4, asuint(sb[0])); }
}
)";

// C18 multicolour block ILU(0) sweeps in the colour-permuted space. One thread per block row of the current colour:
// rows [ng, ng + nbr) (cbuffer fields reused: ng = first row of the colour, nbr = rows in the colour). Rows of one
// colour never reference each other, and earlier (forward) / later (backward) colours were finished by the previous
// dispatches (D3D11 orders UAV accesses between dispatches), so each colour is one fully parallel dispatch.
const char* kIluFwd = R"(
ByteAddressBuffer rowPtr : register(t0);
ByteAddressBuffer colIdx : register(t1);
ByteAddressBuffer lu     : register(t2);
ByteAddressBuffer diagIx : register(t3);
ByteAddressBuffer rv     : register(t4);
RWByteAddressBuffer yv   : register(u0);
@LDY@
[numthreads(64, 1, 1)]
void main(uint3 g : SV_GroupID, uint3 l : SV_GroupThreadID) {
    uint q = gidx(g, l, 64u);
    if (q >= nbr) return;
    uint i = ng + q;
    float acc[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    for (uint e0 = 0; e0 < nb; ++e0) acc[e0] = asfloat(rv.Load((i * nb + e0) * 4));
    uint k0 = rowPtr.Load(i * 4), kd = diagIx.Load(i * 4);
    for (uint k = k0; k < kd; ++k) {
        uint j = colIdx.Load(k * 4);
        uint base = k * nb * nb;
        for (uint e = 0; e < nb; ++e) {
            float s = 0;
            for (uint f = 0; f < nb; ++f) s += asfloat(lu.Load((base + e * nb + f) * 4)) * LDY((j * nb + f) * 4);
            acc[e] -= s;
        }
    }
    for (uint e2 = 0; e2 < nb; ++e2) yv.Store((i * nb + e2) * 4, asuint(acc[e2]));
}
)";
const char* kIluBwd = R"(
ByteAddressBuffer rowPtr : register(t0);
ByteAddressBuffer colIdx : register(t1);
ByteAddressBuffer lu     : register(t2);
ByteAddressBuffer diagIx : register(t3);
RWByteAddressBuffer yv   : register(u0);
@LDY@
[numthreads(64, 1, 1)]
void main(uint3 g : SV_GroupID, uint3 l : SV_GroupThreadID) {
    uint q = gidx(g, l, 64u);
    if (q >= nbr) return;
    uint i = ng + q;
    float acc[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    for (uint e0 = 0; e0 < nb; ++e0) acc[e0] = LDY((i * nb + e0) * 4);
    uint kd = diagIx.Load(i * 4), k1 = rowPtr.Load(i * 4 + 4);
    for (uint k = kd + 1; k < k1; ++k) {
        uint j = colIdx.Load(k * 4);
        uint base = k * nb * nb;
        for (uint e = 0; e < nb; ++e) {
            float s = 0;
            for (uint f = 0; f < nb; ++f) s += asfloat(lu.Load((base + e * nb + f) * 4)) * LDY((j * nb + f) * 4);
            acc[e] -= s;
        }
    }
    uint db = kd * nb * nb;
    for (uint e2 = 0; e2 < nb; ++e2) {
        float s2 = 0;
        for (uint f2 = 0; f2 < nb; ++f2) s2 += asfloat(lu.Load((db + e2 * nb + f2) * 4)) * acc[f2];
        yv.Store((i * nb + e2) * 4, asuint(s2));
    }
}
)";

// copy rows [ng, ng + nbr) (all nb entries) from T to O - used by the copy-based ILU fallback
const char* kCopyRows = R"(
ByteAddressBuffer T : register(t0);
RWByteAddressBuffer O : register(u0);
[numthreads(256, 1, 1)]
void main(uint3 g : SV_GroupID, uint3 l : SV_GroupThreadID) {
    uint q = gidx(g, l, 256u);
    if (q >= nbr * nb) return;
    uint i = ng * nb + q;
    O.Store(i * 4, T.Load(i * 4));
}
)";
// fallback reduction without groupshared memory / barriers (for HLSL compilers lacking them, e.g. vkd3d in Wine):
// every thread sums 32 consecutive elements, the CPU adds the partials
const char* kDot2s = R"(
ByteAddressBuffer X : register(t0);
ByteAddressBuffer Y : register(t1);
ByteAddressBuffer Z : register(t2);
RWByteAddressBuffer O : register(u0);
[numthreads(64, 1, 1)]
void main(uint3 g : SV_GroupID, uint3 l : SV_GroupThreadID) {
    uint t = gidx(g, l, 64u);
    if (t >= ng) return;
    float a = 0, b = 0;
    for (uint q = 0; q < 32; ++q) {
        uint i = t * 32u + q;
        if (i < n) { float x = asfloat(X.Load(i * 4)); a += x * asfloat(Y.Load(i * 4)); b += x * asfloat(Z.Load(i * 4)); }
    }
    O.Store(t * 4, asuint(a)); O.Store((ng + t) * 4, asuint(b));
}
)";

struct Params { uint32_t n, nbr, nb, ng; float a, b, c, d; };

struct Kernel { sg_shader sh{}; sg_pipeline pip{}; int nIn = 0; };
Kernel g_spmv, g_prec, g_axpy, g_dot, g_fwd, g_bwd, g_copy;
bool g_iluCopy = false;         // ILU via SRV reads + scratch copy (compiler fallback)
bool g_dotSimple = false;       // groupshared reduction unavailable -> per-thread partials
bool g_mcOk = false;            // multicolour ILU kernels compiled
int g_precMode = 1;             // 1 = MC-BILU0 (default), 0 = block Jacobi, 2 = C29 CPR-AMG (device V-cycle)

bool makeKernel(Kernel& k, const char* body, int nIn) {
    g_lastLog.clear();
    const std::string src = std::string(kCommon) + body;
    sg_shader_desc d{};
    d.compute_func.source = src.c_str();
    d.compute_func.entry = "main";
    d.compute_func.d3d11_target = "cs_5_0";
    d.uniform_blocks[0].stage = SG_SHADERSTAGE_COMPUTE;
    d.uniform_blocks[0].size = sizeof(Params);
    d.uniform_blocks[0].hlsl_register_b_n = 0;
    for (int i = 0; i < nIn; ++i) {
        d.views[i].storage_buffer.stage = SG_SHADERSTAGE_COMPUTE;
        d.views[i].storage_buffer.readonly = true;
        d.views[i].storage_buffer.hlsl_register_t_n = uint8_t(i);
    }
    d.views[nIn].storage_buffer.stage = SG_SHADERSTAGE_COMPUTE;
    d.views[nIn].storage_buffer.readonly = false;
    d.views[nIn].storage_buffer.hlsl_register_u_n = 0;
    d.label = "ressim-cs";
    k.sh = sg_make_shader(&d);
    if (sg_query_shader_state(k.sh) != SG_RESOURCESTATE_VALID) return false;
    sg_pipeline_desc p{};
    p.compute = true;
    p.shader = k.sh;
    k.pip = sg_make_pipeline(&p);
    k.nIn = nIn;
    return sg_query_pipeline_state(k.pip) == SG_RESOURCESTATE_VALID;
}

bool initLocked() {
    if (g_tried) return g_ok;
    g_tried = true;
    // feature level 11_0 on purpose: the solver needs <= 8 UAV slots, and sokol binds all 64 D3D11.1 UAV slots on an
    // 11_1 device, which some drivers / translation layers (Wine) silently reject -> kernels would write nothing
    const D3D_FEATURE_LEVEL want[] = {D3D_FEATURE_LEVEL_11_0};
    D3D_FEATURE_LEVEL got{};
    HRESULT hr = D3D11CreateDevice(nullptr, D3D_DRIVER_TYPE_HARDWARE, nullptr, 0, want, 1, D3D11_SDK_VERSION, &g_dev, &got, &g_ctx);
    if (FAILED(hr) || !g_dev) { g_err = "no Direct3D 11 hardware device"; return false; }
    std::string adapter = "unknown adapter";
    IDXGIDevice* dx = nullptr;
    if (SUCCEEDED(g_dev->QueryInterface(__uuidof(IDXGIDevice), (void**)&dx)) && dx) {
        IDXGIAdapter* ad = nullptr;
        if (SUCCEEDED(dx->GetAdapter(&ad)) && ad) {
            DXGI_ADAPTER_DESC ds{};
            if (SUCCEEDED(ad->GetDesc(&ds))) {
                char nm[256] = {};
                WideCharToMultiByte(CP_UTF8, 0, ds.Description, -1, nm, sizeof(nm) - 1, nullptr, nullptr);
                adapter = nm;
                char mem[64]; std::snprintf(mem, sizeof(mem), ", %.0f MB", double(ds.DedicatedVideoMemory) / 1048576.0);
                adapter += mem;
            }
            ad->Release();
        }
        dx->Release();
    }
    sg_desc sd{};
    sd.environment.d3d11.device = g_dev;
    sd.environment.d3d11.device_context = g_ctx;
    sd.environment.defaults.color_format = SG_PIXELFORMAT_RGBA8;
    sd.environment.defaults.depth_format = SG_PIXELFORMAT_NONE;
    sd.environment.defaults.sample_count = 1;
    sd.logger.func = sgLog;
    sd.uniform_buffer_size = 4 * 1024 * 1024;
    sg_setup(&sd);
    if (!sg_isvalid()) { g_err = "sokol_gfx setup failed " + g_lastLog; return false; }
    if (!sg_query_features().compute) { g_err = "device has no compute-shader support"; return false; }
    const char* failed = !makeKernel(g_spmv, kSpmv, 4) ? "spmv" : !makeKernel(g_prec, kPrec, 2) ? "jacobi" : !makeKernel(g_axpy, kAxpbypcz, 3) ? "axpy"
                       : nullptr;
    if (!failed && !makeKernel(g_dot, kDot2, 3)) { g_dotSimple = true; if (!makeKernel(g_dot, kDot2s, 3)) failed = "dot"; }
    if (failed) {
        g_err = std::string("compute shader '") + failed + "' creation failed (d3dcompiler_47.dll?) " + g_lastLog;
        return false;
    }
    // reads of the in-place UAV: plain Load (fxc / real drivers); fallback = atomic read (InterlockedAdd 0) for HLSL
    // compilers without RWByteAddressBuffer::Load (vkd3d in Wine) - same algorithm, only the read primitive differs
    // fallback (vkd3d in Wine has neither RW Load nor RW atomics on raw buffers): read the vector through an extra SRV
    // and write the colour's rows into a scratch buffer, then copy them back (2 dispatches per colour, no RW reads)
    auto ilu = [](const char* body, int mode, int srvSlot) {
        std::string s2 = body;
        const std::string def = mode == 0 ? "float LDY(uint a) { return asfloat(yv.Load(a)); }"
                                          : "ByteAddressBuffer ys : register(t" + std::to_string(srvSlot) + ");\nfloat LDY(uint a) { return asfloat(ys.Load(a)); }";
        s2.replace(s2.find("@LDY@"), 5, def);
        return s2;
    };
    static std::string fw0, bw0, fw1, bw1;
    fw0 = ilu(kIluFwd, 0, 0); bw0 = ilu(kIluBwd, 0, 0); fw1 = ilu(kIluFwd, 1, 5); bw1 = ilu(kIluBwd, 1, 4);
    g_mcOk = makeKernel(g_fwd, fw0.c_str(), 5) && makeKernel(g_bwd, bw0.c_str(), 4);
    if (!g_mcOk) {
        g_iluCopy = true;
        g_mcOk = makeKernel(g_fwd, fw1.c_str(), 6) && makeKernel(g_bwd, bw1.c_str(), 5) && makeKernel(g_copy, kCopyRows, 1);
    }
    char fl[32]; std::snprintf(fl, sizeof(fl), "FL %d.%d", (int(got) >> 12) & 0xF, (int(got) >> 8) & 0xF);
    g_info = adapter + " (" + fl + ", sokol_gfx D3D11 compute, float32" + (g_mcOk ? ", MC-BILU0 kernels" : ", block-Jacobi only: ILU kernels failed " + g_lastLog) +
             (g_dotSimple ? ", simple reductions" : "") + (g_iluCopy ? ", copy-based ILU sweeps" : "") + ")";
    g_ok = true;
    return true;
}

// ---------------------------------------------------------------- device buffers
struct Buf {
    sg_buffer b{}; sg_view v{}; size_t bytes = 0;
    void make(const void* data, size_t nbytes) {
        bytes = std::max<size_t>(nbytes, 16);
        std::vector<uint8_t> zero;
        if (!data || nbytes < bytes) { zero.assign(bytes, 0); if (data) std::memcpy(zero.data(), data, nbytes); data = zero.data(); }
        sg_buffer_desc d{};
        d.size = bytes; d.usage.storage_buffer = true; d.usage.vertex_buffer = false;
        d.data.ptr = data; d.data.size = bytes;
        b = sg_make_buffer(&d);
        sg_view_desc vd{}; vd.storage_buffer.buffer = b;
        v = sg_make_view(&vd);
    }
    void destroy() { if (v.id) sg_destroy_view(v); if (b.id) sg_destroy_buffer(b); v = {}; b = {}; }
    bool ok() const { return sg_query_buffer_state(b) == SG_RESOURCESTATE_VALID && sg_query_view_state(v) == SG_RESOURCESTATE_VALID; }
};

bool readBack(const Buf& src, void* dst, size_t nbytes) {
    ID3D11Buffer* sb = (ID3D11Buffer*)sg_d3d11_query_buffer_info(src.b).buf;
    if (!sb) return false;
    D3D11_BUFFER_DESC d{};
    d.ByteWidth = UINT(src.bytes); d.Usage = D3D11_USAGE_STAGING; d.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
    ID3D11Buffer* st = nullptr;
    if (FAILED(g_dev->CreateBuffer(&d, nullptr, &st)) || !st) return false;
    D3D11_BOX box{0, 0, 0, UINT(nbytes), 1, 1};
    g_ctx->CopySubresourceRegion(st, 0, 0, 0, 0, sb, 0, &box);
    D3D11_MAPPED_SUBRESOURCE m{};
    const bool ok = SUCCEEDED(g_ctx->Map(st, 0, D3D11_MAP_READ, 0, &m));
    if (ok) { std::memcpy(dst, m.pData, nbytes); g_ctx->Unmap(st, 0); }
    st->Release();
    return ok;
}

void dispatch(const Kernel& k, const Params& p, std::initializer_list<const Buf*> in, const Buf& out, uint32_t threads, uint32_t perGroup) {
    sg_bindings bnd{};
    int i = 0;
    for (const Buf* b : in) bnd.views[i++] = b->v;
    bnd.views[k.nIn] = out.v;
    sg_apply_pipeline(k.pip);
    sg_apply_bindings(&bnd);
    sg_range r{&p, sizeof(p)};
    sg_apply_uniforms(0, &r);
    const uint32_t groups = (threads + perGroup - 1) / perGroup;
    sg_dispatch(int(std::min<uint32_t>(groups, 65535u)), int((groups + 65534u) / 65535u), 1);
}

// C18 real device timing: D3D11 timestamp queries around every compute pass (disjoint query gives the tick rate);
// busy = sum of pass durations measured on the GPU clock, independent of CPU-side submission and readback latency
struct GpuTimer {
    ID3D11Query* disj = nullptr;
    std::vector<ID3D11Query*> qs;
    int used = 0; bool on = false;
    void begin() {
        D3D11_QUERY_DESC d{D3D11_QUERY_TIMESTAMP_DISJOINT, 0};
        on = SUCCEEDED(g_dev->CreateQuery(&d, &disj)) && disj;
        if (on) g_ctx->Begin(disj);
        used = 0;
    }
    void stamp() {
        if (!on) return;
        if (used >= 8192) return;
        if (used >= int(qs.size())) {
            D3D11_QUERY_DESC d{D3D11_QUERY_TIMESTAMP, 0}; ID3D11Query* q = nullptr;
            if (FAILED(g_dev->CreateQuery(&d, &q)) || !q) return;
            qs.push_back(q);
        }
        g_ctx->End(qs[size_t(used++)]);
    }
    // pairs (stamp before / after each pass); returns busy seconds, -1 when unavailable
    double end() {
        if (!on) return -1;
        g_ctx->End(disj);
        D3D11_QUERY_DATA_TIMESTAMP_DISJOINT dj{};
        auto t0 = std::chrono::steady_clock::now();
        while (g_ctx->GetData(disj, &dj, sizeof(dj), 0) == S_FALSE) {
            if (std::chrono::steady_clock::now() - t0 > std::chrono::seconds(5)) break;
            Sleep(0);
        }
        double busy = -1;
        if (!dj.Disjoint && dj.Frequency > 0) {
            busy = 0;
            for (int k = 0; k + 1 < used; k += 2) {
                UINT64 a = 0, b = 0;
                if (g_ctx->GetData(qs[size_t(k)], &a, sizeof(a), 0) != S_OK || g_ctx->GetData(qs[size_t(k) + 1], &b, sizeof(b), 0) != S_OK) { busy = -1; break; }
                if (b > a) busy += double(b - a) / double(dj.Frequency);
            }
        }
        disj->Release(); disj = nullptr; on = false;
        return busy;
    }
    ~GpuTimer() { for (auto* q : qs) q->Release(); if (disj) disj->Release(); }
};

bool invertSmall(double* a, int n) {   // Gauss-Jordan with partial pivoting, row-major in place
    std::vector<double> inv(size_t(n) * n, 0.0);
    for (int i = 0; i < n; ++i) inv[size_t(i) * n + i] = 1;
    for (int c = 0; c < n; ++c) {
        int piv = c;
        for (int r = c + 1; r < n; ++r) if (std::fabs(a[r * n + c]) > std::fabs(a[piv * n + c])) piv = r;
        if (!(std::fabs(a[piv * n + c]) > 0)) return false;
        if (piv != c) for (int k = 0; k < n; ++k) { std::swap(a[c * n + k], a[piv * n + k]); std::swap(inv[size_t(c) * n + k], inv[size_t(piv) * n + k]); }
        const double d = a[c * n + c];
        for (int k = 0; k < n; ++k) { a[c * n + k] /= d; inv[size_t(c) * n + k] /= d; }
        for (int r = 0; r < n; ++r) if (r != c) {
            const double f = a[r * n + c];
            if (f == 0) continue;
            for (int k = 0; k < n; ++k) { a[r * n + k] -= f * a[c * n + k]; inv[size_t(r) * n + k] -= f * inv[size_t(c) * n + k]; }
        }
    }
    std::memcpy(a, inv.data(), sizeof(double) * size_t(n) * n);
    return true;
}

class SokolSolver : public ressim::GpuBlockSolver {
public:
    SolveTiming last;
    const char* name() const override {
        return g_precMode == 2 ? "sokol_gfx D3D11 compute (float32 BiCGSTAB + CPR-AMG device V-cycle)"
             : useMc() ? "sokol_gfx D3D11 compute (float32 BiCGSTAB + multicolour BILU0)" : "sokol_gfx D3D11 compute (float32 BiCGSTAB + block Jacobi)";
    }
    bool useMc() const { return g_precMode == 1 && g_mcOk; }
    bool solve(const ressim::BSR& A0, const std::vector<double>& bd0, std::vector<double>& xd, double tol, int maxIt, int& iterations) override {
        std::lock_guard<std::mutex> lk(g_mx);
        iterations = 0;
        last = SolveTiming{};
        if (!initLocked() || A0.nb > 8 || A0.n <= 0) return false;
        using clk = std::chrono::steady_clock;
        auto secs = [](clk::time_point a, clk::time_point b) { return std::chrono::duration<double>(b - a).count(); };
        const auto tStart = clk::now();
        const bool mc = useMc();
        last.multicolour = mc;
        // multicolour: factor P A P^T on the CPU in double (OpenMP over the rows of each colour); the device iterates
        // entirely in the permuted space (A, b permuted on upload, x un-permuted after read-back)
        std::vector<double> bperm;
        if (mc) {
            mcPrec_.keepPermuted = true;
            mcPrec_.factor(A0);
            const ressim::Colouring& C = mcPrec_.colouring();
            bperm.resize(bd0.size());
            for (int k = 0; k < A0.n; ++k) for (int e = 0; e < A0.nb; ++e) bperm[size_t(k) * A0.nb + e] = bd0[size_t(C.perm[size_t(k)]) * A0.nb + e];
            last.colours = C.colours;
        }
        const ressim::BSR& A = mc ? mcPrec_.permuted() : A0;
        const std::vector<double>& bd = mc ? bperm : bd0;
        const uint32_t nbr = uint32_t(A.n), nb = uint32_t(A.nb), N = nbr * nb;
        const uint32_t ngDot = g_dotSimple ? (N + 31) / 32 : (N + 1023) / 1024;
        // block-diagonal inverse (double on the CPU, stored float); MC mode: float copy of the ILU factors
        std::vector<float> dinv, luF;
        std::vector<int> diagIx;
        std::vector<double> blk(size_t(nb) * nb);
        if (mc) {
            const ressim::BSR& F = mcPrec_.factors();
            luF.resize(F.val.size());
            for (size_t e = 0; e < luF.size(); ++e) luF[e] = float(F.val[e]);
            diagIx = F.diag;
        } else dinv.resize(size_t(nbr) * nb * nb);
        for (uint32_t i = 0; i < nbr && !mc; ++i) {
            const double* D = A.blk(A.diag[i]);
            std::copy(D, D + nb * nb, blk.begin());
            if (!invertSmall(blk.data(), int(nb))) { std::fill(blk.begin(), blk.end(), 0.0); for (uint32_t e = 0; e < nb; ++e) blk[e * nb + e] = 1.0 / std::max(1e-30, std::fabs(D[e * nb + e])); }
            for (uint32_t e = 0; e < nb * nb; ++e) dinv[size_t(i) * nb * nb + e] = float(blk[e]);
        }
        std::vector<float> vals(A.val.size()), bf(N);
        for (size_t e = 0; e < vals.size(); ++e) vals[e] = float(A.val[e]);
        for (uint32_t e = 0; e < N; ++e) bf[e] = float(bd[e]);
        const auto tPrep = clk::now();
        last.setup = secs(tStart, tPrep);
        Buf rp, ci, va, di, r, r0, p, p2, v, ph, s, sh, t, x, x2, part;
        rp.make(A.rowPtr.data(), A.rowPtr.size() * 4); ci.make(A.col.data(), A.col.size() * 4); va.make(vals.data(), vals.size() * 4);
        Buf lub, dgb, scr;
        scr.make(nullptr, (mc && g_iluCopy) ? size_t(N) * 4 : 16);
        if (mc) { lub.make(luF.data(), luF.size() * 4); dgb.make(diagIx.data(), diagIx.size() * 4); di.make(nullptr, 16); }
        else { di.make(dinv.data(), dinv.size() * 4); lub.make(nullptr, 16); dgb.make(nullptr, 16); }
        const size_t vb = size_t(N) * 4;
        r.make(bf.data(), vb); r0.make(bf.data(), vb);
        for (Buf* q : {&p, &p2, &v, &ph, &s, &sh, &t, &x, &x2}) q->make(nullptr, vb);
        part.make(nullptr, size_t(2) * ngDot * 4);
        bool good = true;
        for (Buf* q : {&rp, &ci, &va, &di, &r, &r0, &p, &p2, &v, &ph, &s, &sh, &t, &x, &x2, &part, &lub, &dgb, &scr}) good = good && q->ok();
        // C29 CPR-AMG on the device: host hierarchy (ressim::DeviceCprAmg, the algorithm the float32 emulation tests),
        // every level's operator / prolongation / restriction uploaded as nb = 1 CSR; the V-cycle is a fixed sequence of
        // the existing SpMV / Jacobi / axpy kernels (no host round trip inside the preconditioner)
        struct DM { Buf rp, ci, va; uint32_t rows = 0; void destroy() { rp.destroy(); ci.destroy(); va.destroy(); } };
        struct DL { DM A, P, R; Buf dinv, x, x2, b, r, z, t; uint32_t n = 0; };
        ressim::DeviceCprAmg dc;
        const bool cpr = g_precMode == 2 && !mc && dc.build(A0);
        std::vector<DL> dl;
        DM dW, dE, dCi;
        Buf dz1, dt2, dt3;
        auto up = [&](DM& M, const ressim::DevCsr& C) {
            M.rows = uint32_t(C.rows);
            M.rp.make(C.rp.data(), C.rp.size() * 4); M.ci.make(C.ci.data(), C.ci.size() * 4); M.va.make(C.v.data(), C.v.size() * 4);
            good = good && M.rp.ok() && M.ci.ok() && M.va.ok();
        };
        if (cpr) {
            up(dW, dc.W); up(dE, dc.E);
            if (!dc.coarseInv.rp.empty()) up(dCi, dc.coarseInv);
            dl.resize(dc.lv.size());
            for (size_t l = 0; l < dc.lv.size(); ++l) {
                DL& V = dl[l];
                V.n = uint32_t(dc.lv[l].A.rows);
                up(V.A, dc.lv[l].A);
                if (!dc.lv[l].P.rp.empty()) { up(V.P, dc.lv[l].P); up(V.R, dc.lv[l].R); }
                std::vector<float> wd(dc.lv[l].dinv);
                for (float& q : wd) q *= dc.omega;                      // omega D^-1 (damped Jacobi)
                V.dinv.make(wd.data(), wd.size() * 4);
                for (Buf* q : {&V.x, &V.x2, &V.b, &V.r, &V.z, &V.t}) { q->make(nullptr, size_t(V.n) * 4); good = good && q->ok(); }
                good = good && V.dinv.ok();
            }
            for (Buf* q : {&dz1, &dt2, &dt3}) { q->make(nullptr, vb); good = good && q->ok(); }
        }
        auto sp1 = [&](const DM& M, const Buf& X, const Buf& Y) { Params pr{M.rows, M.rows, 1, 0, 0, 0, 0, 0}; dispatch(g_spmv, pr, {&M.rp, &M.ci, &M.va, &X}, Y, M.rows, 64); };
        auto ax1 = [&](uint32_t n1, const Buf& X, const Buf& Y, const Buf& Z, const Buf& O, float a, float b, float c) {
            Params pr{n1, n1, 1, 0, a, b, c, 0}; dispatch(g_axpy, pr, {&X, &Y, &Z}, O, n1, 256);
        };
        auto jac = [&](uint32_t n1, const Buf& D, const Buf& X, const Buf& Z) { Params pr{n1, n1, 1, 0, 0, 0, 0, 0}; dispatch(g_prec, pr, {&D, &X}, Z, n1, 64); };
        auto smooth = [&](DL& V) {
            sp1(V.A, V.x, V.t); ax1(V.n, V.b, V.t, V.t, V.r, 1.f, -1.f, 0.f); jac(V.n, V.dinv, V.r, V.z);
            ax1(V.n, V.x, V.z, V.z, V.x2, 1.f, 1.f, 0.f); std::swap(V.x, V.x2);
        };
        GpuTimer tm;
        tm.begin();
        std::vector<float> pv(size_t(2) * ngDot);
        auto dot2 = [&](const Buf& X, const Buf& Y, const Buf& Z, double& a, double& b) {
            Params pr{N, nbr, nb, ngDot, 0, 0, 0, 0};
            sg_pass ps{}; ps.compute = true;
            sg_begin_pass(&ps);
            tm.stamp();
            if (g_dotSimple) dispatch(g_dot, pr, {&X, &Y, &Z}, part, ngDot, 64);
            else dispatch(g_dot, pr, {&X, &Y, &Z}, part, ngDot * 256, 256);
            tm.stamp();
            sg_end_pass();
            sg_commit();
            ++last.readbacks;
            if (!readBack(part, pv.data(), pv.size() * 4)) return false;
            a = 0; b = 0;
            for (uint32_t gI = 0; gI < ngDot; ++gI) { a += pv[gI]; b += pv[ngDot + gI]; }
            return std::isfinite(a) && std::isfinite(b);
        };
        auto axpy = [&](const Buf& X, const Buf& Y, const Buf& Z, const Buf& O, double a, double b, double c) {
            Params pr{N, nbr, nb, ngDot, float(a), float(b), float(c), 0};
            dispatch(g_axpy, pr, {&X, &Y, &Z}, O, N, 256);
        };
        auto spmv = [&](const Buf& X, const Buf& Y) { Params pr{N, nbr, nb, ngDot, 0, 0, 0, 0}; dispatch(g_spmv, pr, {&rp, &ci, &va, &X}, Y, nbr, 64); };
        const ressim::Colouring* Cc = mc ? &mcPrec_.colouring() : nullptr;
        auto cprPrec = [&](const Buf& R, const Buf& Z) {
            const size_t L = dl.size();
            sp1(dW, R, dl[0].b);
            for (size_t l = 0; l + 1 < L; ++l) {
                DL& V = dl[l];
                jac(V.n, V.dinv, V.b, V.x);
                for (int q = 1; q < dc.pre; ++q) smooth(V);
                sp1(V.A, V.x, V.t); ax1(V.n, V.b, V.t, V.t, V.r, 1.f, -1.f, 0.f);
                sp1(V.R, V.r, dl[l + 1].b);
            }
            DL& C = dl[L - 1];
            if (dCi.rows) sp1(dCi, C.b, C.x);
            else { jac(C.n, C.dinv, C.b, C.x); for (int q = 0; q < 30; ++q) smooth(C); }
            for (size_t l = L - 1; l-- > 0;) {
                DL& V = dl[l];
                sp1(V.P, dl[l + 1].x, V.t); ax1(V.n, V.x, V.t, V.t, V.x2, 1.f, 1.f, 0.f); std::swap(V.x, V.x2);
                for (int q = 0; q < dc.post; ++q) smooth(V);
            }
            sp1(dE, dl[0].x, dz1);
            spmv(dz1, dt2);
            axpy(R, dt2, dt2, dt3, 1.0, -1.0, 0.0);
            { Params pr{N, nbr, nb, ngDot, 0, 0, 0, 0}; dispatch(g_prec, pr, {&di, &dt3}, dt2, nbr, 64); }
            axpy(dz1, dt2, dt2, Z, 1.0, 1.0, 0.0);
        };
        auto prec = [&](const Buf& X, const Buf& Z) {
            if (cpr) { cprPrec(X, Z); return; }
            if (!mc) { Params pr{N, nbr, nb, ngDot, 0, 0, 0, 0}; dispatch(g_prec, pr, {&di, &X}, Z, nbr, 64); return; }
            for (int c = 0; c < Cc->colours; ++c) {
                const uint32_t s0 = uint32_t(Cc->start[size_t(c)]), cnt = uint32_t(Cc->start[size_t(c) + 1]) - s0;
                if (!cnt) continue;
                Params pr{N, cnt, nb, s0, 0, 0, 0, 0};
                if (g_iluCopy) { dispatch(g_fwd, pr, {&rp, &ci, &lub, &dgb, &X, &Z}, scr, cnt, 64); dispatch(g_copy, pr, {&scr}, Z, cnt * nb, 256); }
                else dispatch(g_fwd, pr, {&rp, &ci, &lub, &dgb, &X}, Z, cnt, 64);
            }
            for (int c = Cc->colours - 1; c >= 0; --c) {
                const uint32_t s0 = uint32_t(Cc->start[size_t(c)]), cnt = uint32_t(Cc->start[size_t(c) + 1]) - s0;
                if (!cnt) continue;
                Params pr{N, cnt, nb, s0, 0, 0, 0, 0};
                if (g_iluCopy) { dispatch(g_bwd, pr, {&rp, &ci, &lub, &dgb, &Z}, scr, cnt, 64); dispatch(g_copy, pr, {&scr}, Z, cnt * nb, 256); }
                else dispatch(g_bwd, pr, {&rp, &ci, &lub, &dgb}, Z, cnt, 64);
            }
        };
        double bn = 0; for (uint32_t e = 0; e < N; ++e) bn += bd[e] * bd[e];
        bn = std::sqrt(bn);
        xd.assign(N, 0.0);
        bool conv = false;
        if (good && bn > 0) {
            double rho = 1, alpha = 1, omega = 1;
            const int itMax = std::max(maxIt, 400);
            for (int it = 0; it < itMax; ++it) {
                double rhoN = 0, rr = 0;
                if (!dot2(r, r0, r, rhoN, rr)) { good = false; break; }
                // a device that silently drops the UAV writes (vkd3d / FL 11_1 drivers) returns r = 0 at the first check
                // although b != 0 (x0 = 0 -> r0 = b): treat that as a device failure, never as convergence
                if (!std::isfinite(rr) || (it == 0 && rr <= 1e-300)) { good = false; g_err = "GPU returned an invalid residual (device writes lost)"; break; }
                if (std::sqrt(rr) <= tol * bn) { conv = true; break; }
                if (rhoN == 0) break;
                iterations = it + 1;
                const double beta = (rhoN / rho) * (alpha / omega);
                sg_pass ps{}; ps.compute = true;
                sg_begin_pass(&ps);
                tm.stamp();
                axpy(r, p, v, p2, 1.0, beta, -beta * omega);
                std::swap(p, p2);
                prec(p, ph);
                spmv(ph, v);
                tm.stamp();
                sg_end_pass();
                double r0v = 0, dummy = 0;
                if (!dot2(r0, v, v, r0v, dummy) || r0v == 0) { good = r0v != 0; break; }
                alpha = rhoN / r0v;
                sg_begin_pass(&ps);
                tm.stamp();
                axpy(r, v, v, s, 1.0, -alpha, 0.0);
                prec(s, sh);
                spmv(sh, t);
                tm.stamp();
                sg_end_pass();
                double tt = 0, ts = 0;
                if (!dot2(t, t, s, tt, ts)) { good = false; break; }
                omega = tt > 0 ? ts / tt : 0.0;
                sg_begin_pass(&ps);
                tm.stamp();
                axpy(x, ph, sh, x2, 1.0, alpha, omega);
                std::swap(x, x2);
                axpy(s, t, t, r, 1.0, -omega, 0.0);
                tm.stamp();
                sg_end_pass();
                rho = rhoN;
                if (omega == 0) break;
            }
            std::vector<float> xf(N);
            if (good && readBack(x, xf.data(), vb)) {
                if (mc) { const auto& pm = mcPrec_.colouring().perm; for (uint32_t k = 0; k < nbr; ++k) for (uint32_t e = 0; e < nb; ++e) xd[size_t(pm[k]) * nb + e] = xf[size_t(k) * nb + e]; }
                else for (uint32_t e = 0; e < N; ++e) xd[e] = xf[e];
            } else good = false;
        }
        last.deviceBusy = tm.end();
        if (last.deviceBusy >= 0) busyTotal_ += last.deviceBusy;
        for (Buf* q : {&rp, &ci, &va, &di, &r, &r0, &p, &p2, &v, &ph, &s, &sh, &t, &x, &x2, &part, &lub, &dgb, &scr}) q->destroy();
        for (DL& V : dl) { V.A.destroy(); V.P.destroy(); V.R.destroy(); for (Buf* q : {&V.dinv, &V.x, &V.x2, &V.b, &V.r, &V.z, &V.t}) q->destroy(); }
        dW.destroy(); dE.destroy(); dCi.destroy(); dz1.destroy(); dt2.destroy(); dt3.destroy();
        sg_commit();
        last.total = secs(tStart, clk::now());
        last.iterate = last.total - last.setup;
        last.iterations = iterations;
        last.converged = conv;
        for (double vv : xd) if (!std::isfinite(vv)) return false;
        return good;
    }
    double deviceSeconds() const override { return busyTotal_; }
private:
    ressim::MCBILU0 mcPrec_;
    double busyTotal_ = 0;
};
SokolSolver g_solver;

} // namespace

bool available(std::string* info) {
    std::lock_guard<std::mutex> lk(g_mx);
    const bool ok = initLocked();
    if (info) *info = ok ? g_info : g_err;
    return ok;
}
std::string deviceInfo() { std::lock_guard<std::mutex> lk(g_mx); return g_ok ? g_info : g_err; }
ressim::GpuBlockSolver* solver() { return &g_solver; }
SolveTiming lastTiming() { std::lock_guard<std::mutex> lk(g_mx); return g_solver.last; }
void setPreconditioner(int mode) { std::lock_guard<std::mutex> lk(g_mx); g_precMode = std::clamp(mode, 0, 2); }
int preconditioner() { std::lock_guard<std::mutex> lk(g_mx); return g_precMode; }

namespace {
// pressure-dominated block system on a heterogeneous 2D (nz = 1) or 3D grid, 2 decades of permeability (the same kind
// as the engine's Schur-reduced Newton system)
ressim::BSR testMatrix(int nx, int ny, int nz, int nb, unsigned seed) {
    std::mt19937 rng(seed); std::uniform_real_distribution<double> U(-1, 1);
    ressim::BSR A; A.n = nx * ny * nz; A.nb = nb; A.rowPtr.assign(1, 0);
    A.col.reserve(size_t(A.n) * 7); A.val.reserve(size_t(A.n) * 7 * nb * nb);
    std::vector<double> kx(size_t(A.n)), kz(size_t(A.n));
    for (auto& v : kx) v = std::pow(10.0, U(rng));
    for (auto& v : kz) v = 0.1 * std::pow(10.0, U(rng));
    std::vector<std::pair<int, double>> row;
    for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        const int c = (k * ny + j) * nx + i;
        row.clear();
        double sum = 0;
        auto add = [&](int o, double t) { row.push_back({o, t}); sum += t; };
        if (k > 0) add(c - nx * ny, std::sqrt(kz[size_t(c)] * kz[size_t(c - nx * ny)]));
        if (j > 0) add(c - nx, std::sqrt(kx[size_t(c)] * kx[size_t(c - nx)]));
        if (i > 0) add(c - 1, std::sqrt(kx[size_t(c)] * kx[size_t(c - 1)]));
        if (i < nx - 1) add(c + 1, std::sqrt(kx[size_t(c)] * kx[size_t(c + 1)]));
        if (j < ny - 1) add(c + nx, std::sqrt(kx[size_t(c)] * kx[size_t(c + nx)]));
        if (k < nz - 1) add(c + nx * ny, std::sqrt(kz[size_t(c)] * kz[size_t(c + nx * ny)]));
        row.push_back({c, sum * 1.05});
        std::sort(row.begin(), row.end());
        for (auto& q : row) {
            A.col.push_back(q.first);
            for (int e = 0; e < nb; ++e) for (int f = 0; f < nb; ++f)
                A.val.push_back(q.first == c ? (e == f ? q.second * (1 + 0.2 * e) : 0.03 * q.second) : (e == f ? -q.second * (e ? 0.7 : 1.0) : 0.0));
        }
        A.rowPtr.push_back(int(A.col.size()));
    }
    A.finalizeDiag();
    return A;
}
double relErr(const std::vector<double>& x, const std::vector<double>& xt) {
    if (x.size() != xt.size()) return 1e30;
    double e = 0, n = 0;
    for (size_t k = 0; k < x.size(); ++k) { e = std::max(e, std::fabs(x[k] - xt[k])); n = std::max(n, std::fabs(xt[k])); }
    return n > 0 ? e / n : e;
}
} // namespace

std::string selfTest(bool& ok) {
    ok = false;
    std::string info;
    if (!available(&info)) return "GPU not available: " + info;
    const ressim::BSR A = testMatrix(120, 100, 1, 3, 5);
    const int N = A.size();
    std::vector<double> xt(static_cast<size_t>(N)), b(static_cast<size_t>(N));
    for (int e = 0; e < N; ++e) xt[size_t(e)] = std::cos(0.003 * e) + 0.2 * (e % 3);
    A.spmv(xt.data(), b.data());
    using clk = std::chrono::steady_clock;
    std::string rep = info;
    bool allOk = true;
    const int keep = preconditioner();
    for (int mode : {1, 0, 2}) {
        if (mode == 1 && !g_mcOk) continue;
        setPreconditioner(mode);
        auto t0 = clk::now();
        std::vector<double> xg;
        const ressim::LinResult rg = ressim::refineWithDevice(g_solver, A, b, xg, 1e-10, 8, 2000);
        const double tg = std::chrono::duration<double>(clk::now() - t0).count();
        const SolveTiming lt = lastTiming();
        const double eg = relErr(xg, xt);
        allOk = allOk && rg.converged && eg < 1e-8;
        char buf[320];
        std::snprintf(buf, sizeof(buf), " | GPU %s: %d its %.3f s (last inner solve: device busy %s, wall %.1f ms) err %.1e%s", mode == 2 ? "CPR-AMG" : mode ? "MC-BILU0" : "Jacobi",
                      rg.iterations, tg, lt.deviceBusy >= 0 ? (std::to_string(int(lt.deviceBusy * 1e3 + 0.5)) + " ms").c_str() : "n/a", lt.total * 1e3, eg,
                      rg.converged ? "" : " NOT converged");
        rep += buf;
    }
    setPreconditioner(keep);
    auto t0 = clk::now();
    ressim::BILU0 M; M.factor(A);
    ressim::BsrOperator op(A);
    std::vector<double> xc;
    const ressim::LinResult rc = ressim::bicgstab(op, M, b, xc, 1e-10, 2000);
    const double tc = std::chrono::duration<double>(clk::now() - t0).count();
    char buf[200];
    std::snprintf(buf, sizeof(buf), " | CPU BILU0 %d its %.3f s err %.1e (%d unknowns)", rc.iterations, tc, relErr(xc, xt), N);
    rep += buf;
    ok = allOk;
    return rep;
}

std::string benchmark(int maxCells, int threads) {
    std::string info;
    const bool have = available(&info);
    std::string out = "ResSim linear-solver benchmark (BiCGSTAB to 1e-8 relative residual, double-accurate answers)\n";
    out += std::string("device: ") + info + "\n";
#ifdef _OPENMP
    if (threads > 0) omp_set_num_threads(threads);
    char th[64]; std::snprintf(th, sizeof(th), "CPU threads: %d\n", omp_get_max_threads()); out += th;
#else
    (void)threads;
#endif
    out += "case                     unknowns | solver            colours  its   wall [s]  device busy [s]  readbacks  rel.err\n";
    struct Sz { int nx, ny, nz; };
    const Sz sizes[] = {{100, 100, 1}, {60, 60, 10}, {300, 300, 1}, {100, 100, 20}, {600, 600, 1}, {150, 150, 30}};
    using clk = std::chrono::steady_clock;
    const int keep = have ? preconditioner() : 1;
    for (const Sz& z : sizes) {
        if (z.nx * z.ny * z.nz > maxCells) continue;
        const ressim::BSR A = testMatrix(z.nx, z.ny, z.nz, 3, 7);
        const int N = A.size();
        std::vector<double> xt(static_cast<size_t>(N)), b(static_cast<size_t>(N));
        for (int e = 0; e < N; ++e) xt[size_t(e)] = std::cos(0.001 * e) + 0.2 * (e % 3);
        A.spmv(xt.data(), b.data());
        char nm[64]; std::snprintf(nm, sizeof(nm), "%s %dx%dx%d nb3", z.nz > 1 ? "3D 7-pt" : "2D 5-pt", z.nx, z.ny, z.nz);
        auto line = [&](const char* solver, int colours, int its, double wall, double busy, int rb, double err) {
            char l[300];
            std::snprintf(l, sizeof(l), "%-24s %9d | %-17s %7s %5d %10.3f %16s %10s  %.1e\n", nm, N, solver, colours > 0 ? std::to_string(colours).c_str() : "-", its, wall,
                          busy >= 0 ? std::to_string(busy).substr(0, 8).c_str() : "-", rb >= 0 ? std::to_string(rb).c_str() : "-", err);
            out += l;
        };
        ressim::BsrOperator op(A);
        {
            auto t0 = clk::now();
            ressim::BILU0 M; M.factor(A);
            std::vector<double> x; const ressim::LinResult r = ressim::bicgstab(op, M, b, x, 1e-8, 5000);
            line("CPU BILU0", 0, r.iterations, std::chrono::duration<double>(clk::now() - t0).count(), -1, -1, relErr(x, xt));
        }
        {
            auto t0 = clk::now();
            ressim::MCBILU0 M; M.factor(A);
            std::vector<double> x; const ressim::LinResult r = ressim::bicgstab(op, M, b, x, 1e-8, 5000);
            line("CPU MC-BILU0", M.colours(), r.iterations, std::chrono::duration<double>(clk::now() - t0).count(), -1, -1, relErr(x, xt));
        }
        if (!have) continue;
        for (int mode : {0, 1}) {
            if (mode == 1 && !g_mcOk) continue;
            setPreconditioner(mode);
            double busy = 0; int rb = 0; bool busyOk = true; int colours = 0;
            // accumulate device timing over all refinement inner solves
            struct Acc : ressim::GpuBlockSolver {
                ressim::GpuBlockSolver* inner; double* busy; int* rb; bool* ok; int* col;
                const char* name() const override { return inner->name(); }
                bool solve(const ressim::BSR& A2, const std::vector<double>& b2, std::vector<double>& x2, double tol, int maxIt, int& its) override {
                    const bool r = inner->solve(A2, b2, x2, tol, maxIt, its);
                    const SolveTiming t = lastTiming();
                    if (t.deviceBusy >= 0) *busy += t.deviceBusy; else *ok = false;
                    *rb += t.readbacks; *col = t.colours;
                    return r;
                }
            } acc;
            acc.inner = &g_solver; acc.busy = &busy; acc.rb = &rb; acc.ok = &busyOk; acc.col = &colours;
            auto t0 = clk::now();
            std::vector<double> x;
            const ressim::LinResult r = ressim::refineWithDevice(acc, A, b, x, 1e-8, 8, 5000);
            line(mode ? "GPU MC-BILU0 f32" : "GPU Jacobi f32", mode ? colours : 0, r.iterations, std::chrono::duration<double>(clk::now() - t0).count(), busyOk ? busy : -1,
                 rb, r.converged ? relErr(x, xt) : 1.0);
        }
    }
    if (have) setPreconditioner(keep);
    out += "wall = factor + solve incl. upload / read-back; device busy = sum of GPU timestamp intervals of the compute passes\n";
    return out;
}

} // namespace gpu
