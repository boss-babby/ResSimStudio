// GPU linear-solver accelerator: Direct3D 11 compute shaders.
//
// The engine solves J dx = -R with a double-precision flexible GMRES on the CPU. When the user selects the GPU
// solver, the preconditioner of that FGMRES is this device solve: the cell system B' (3x3 blocks, wells folded in by
// a Schur complement) is block-Jacobi scaled in double on the CPU (A^ = D^-1 B'), uploaded once per linear solve in
// single precision, and every preconditioner application runs a few BiCGStab iterations on the GPU.
// Single precision on the device + double outer iteration = mixed-precision iterative refinement: the final
// accuracy is that of the CPU FGMRES, the bandwidth-heavy work (SpMV, vector updates, reductions) runs on the GPU.
#include "gpu_solver.h"
#include "ressim/blocksolve.h"

#include <windows.h>
#include <d3d11.h>
#include <d3dcompiler.h>
#include <dxgi.h>

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstring>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace {

const char* kCS = R"(
cbuffer P : register(b0) { uint n; uint nb; uint na; uint mode; uint oX; uint oY; uint oZ; uint oW; float a; float b; float pad0; float pad1; };
StructuredBuffer<uint>  RP : register(t0);
StructuredBuffer<uint>  CI : register(t1);
StructuredBuffer<float> BV : register(t2);
RWStructuredBuffer<float> V    : register(u0);
RWStructuredBuffer<float> PART : register(u1);

[numthreads(128, 1, 1)]
void CSSpmv(uint3 id : SV_DispatchThreadID) {
    uint i = id.x;
    if (i >= na) return;
    float acc[4] = { 0, 0, 0, 0 };
    uint bs = nb * nb;
    for (uint j = RP[i]; j < RP[i + 1]; ++j) {
        uint c = CI[j] * nb, bo = j * bs;
        [unroll] for (uint k = 0; k < 4; ++k) {
            if (k < nb) {
                float s = 0;
                [unroll] for (uint l = 0; l < 4; ++l) if (l < nb) s += BV[bo + k * nb + l] * V[oX + c + l];
                acc[k] += s;
            }
        }
    }
    [unroll] for (uint k = 0; k < 4; ++k) if (k < nb) V[oY + i * nb + k] = acc[k];
}

[numthreads(256, 1, 1)]
void CSOp(uint3 id : SV_DispatchThreadID) {
    uint e = id.x;
    if (e >= n) return;
    if (mode == 0)      V[oY + e] = V[oX + e] + a * (V[oY + e] - b * V[oZ + e]);
    else if (mode == 1) V[oY + e] = V[oX + e] - a * V[oZ + e];
    else if (mode == 2) V[oY + e] = V[oY + e] + a * V[oX + e] + b * V[oZ + e];
    else if (mode == 4) V[oY + e] = V[oX + e];
    else                V[oY + e] = 0;
}

groupshared float s1[256];
groupshared float s2[256];
[numthreads(256, 1, 1)]
void CSDot(uint3 id : SV_DispatchThreadID, uint3 gt : SV_GroupThreadID, uint3 g : SV_GroupID) {
    uint e = id.x, t = gt.x;
    float x1 = 0, x2 = 0;
    if (e < n) { x1 = V[oX + e] * V[oY + e]; x2 = V[oZ + e] * V[oW + e]; }
    s1[t] = x1; s2[t] = x2;
    GroupMemoryBarrierWithGroupSync();
    for (uint st = 128; st > 0; st >>= 1) {
        if (t < st) { s1[t] += s1[t + st]; s2[t] += s2[t + st]; }
        GroupMemoryBarrierWithGroupSync();
    }
    if (t == 0) { PART[g.x * 2] = s1[0]; PART[g.x * 2 + 1] = s2[0]; }
}
)";

struct Params { UINT n, nb, na, mode, oX, oY, oZ, oW; float a, b, pad0, pad1; };

template <class T> void rel(T*& p) { if (p) { p->Release(); p = nullptr; } }

// vector slots inside the single device vector buffer
enum { X = 0, R, RH, PV, VV, SV, TV, NVEC };

class D3D11Accelerator final : public ressim::BlockAccelerator {
public:
    ~D3D11Accelerator() override { freeMatrix(); freeDevice(); }
    const char* name() const override { return adapter_.c_str(); }

    bool setup(int nCells, int nb, const std::vector<int>& rp, const std::vector<int>& ci, const std::vector<double>& blocks) override {
        if (nb < 1 || nb > 4 || nCells <= 0) return false;
        if (!claim()) return false;
        if (!ensureDevice()) return false;
        na_ = nCells; nb_ = nb; n_ = nCells * nb;
        const int bs = nb * nb;
        // block-Jacobi scaling in double: Dinv_i, A^_ij = Dinv_i B'_ij
        dinv_.assign(size_t(na_) * size_t(bs), 0.0);
        std::vector<float> scaled(blocks.size());
        for (int i = 0; i < na_; ++i) {
            const double* D = nullptr;
            for (int j = rp[size_t(i)]; j < rp[size_t(i) + 1]; ++j) if (ci[size_t(j)] == i) { D = &blocks[size_t(j) * size_t(bs)]; break; }
            double* Di = &dinv_[size_t(i) * size_t(bs)];
            if (!D || !invert(D, Di, nb)) for (int k = 0; k < nb; ++k) Di[k * nb + k] = 1.0;
            for (int j = rp[size_t(i)]; j < rp[size_t(i) + 1]; ++j) {
                const double* Bj = &blocks[size_t(j) * size_t(bs)];
                float* Sj = &scaled[size_t(j) * size_t(bs)];
                for (int r = 0; r < nb; ++r)
                    for (int c = 0; c < nb; ++c) {
                        double s = 0; for (int k = 0; k < nb; ++k) s += Di[r * nb + k] * Bj[k * nb + c];
                        Sj[r * nb + c] = float(s);
                    }
            }
        }
        std::vector<UINT> urp(rp.begin(), rp.end()), uci(ci.begin(), ci.end());
        const size_t nnzb = ci.size();
        if (nnzb == 0) return false;
        bool ok = true;
        if (na_ != capN_ || nnzb > capNnz_ || nb_ != capNb_) {
            freeMatrix();
            ok = makeSrv(urp.data(), sizeof(UINT), urp.size(), bRP_, sRP_) && makeSrv(uci.data(), sizeof(UINT), uci.size(), bCI_, sCI_) &&
                 makeSrv(scaled.data(), sizeof(float), scaled.size(), bBV_, sBV_) &&
                 makeUav(sizeof(float), size_t(n_) * NVEC, bV_, uV_) &&
                 makeUav(sizeof(float), size_t(2 * groups(n_, 256)), bPart_, uPart_) && makeStaging();
            capN_ = na_; capNnz_ = nnzb; capNb_ = nb_;
        } else {
            ctx_->UpdateSubresource(bRP_, 0, nullptr, urp.data(), 0, 0);
            ctx_->UpdateSubresource(bCI_, 0, nullptr, uci.data(), 0, 0);
            ctx_->UpdateSubresource(bBV_, 0, nullptr, scaled.data(), 0, 0);
        }
        if (!ok) { freeMatrix(); return false; }
        return true;
    }

    bool apply(const double* r, double* z) override {
        std::vector<float> rhs(static_cast<size_t>(n_));
        const int bs = nb_ * nb_;
        double bn = 0;
        for (int i = 0; i < na_; ++i) {
            const double* Di = &dinv_[size_t(i) * size_t(bs)];
            for (int k = 0; k < nb_; ++k) {
                double s = 0; for (int l = 0; l < nb_; ++l) s += Di[k * nb_ + l] * r[size_t(i) * size_t(nb_) + size_t(l)];
                rhs[size_t(i) * size_t(nb_) + size_t(k)] = float(s); bn += s * s;
            }
        }
        if (bn == 0) { std::memset(z, 0, sizeof(double) * size_t(n_)); return true; }
        // r = b, rhat = b, x = p = v = 0
        upload(R, rhs.data());
        op(5, X, X, X); op(5, PV, PV, PV); op(5, VV, VV, VV); op(4, R, RH, R);
        double rho = 1, alpha = 1, omega = 1;
        const double tol2 = 1e-4 * bn;                     // 1e-2 relative: an inexact (flexible) preconditioner
        bool ok = true;
        for (int it = 0; it < maxIt_; ++it) {
            double rho1, rr; dot(RH, R, R, R, rho1, rr);
            if (!std::isfinite(rho1) || rho1 == 0) { ok = it > 0; break; }
            const double beta = (rho1 / rho) * (alpha / omega);
            op(0, R, PV, VV, float(beta), float(omega));             // p = r + beta (p - omega v)
            spmv(PV, VV);                                            // v = A p
            double rv, dummy; dot(RH, VV, RH, RH, rv, dummy);
            if (!std::isfinite(rv) || rv == 0) { ok = it > 0; break; }
            alpha = rho1 / rv;
            op(1, R, SV, VV, float(alpha));                          // s = r - alpha v
            spmv(SV, TV);                                            // t = A s
            double ts, tt; dot(TV, SV, TV, TV, ts, tt);
            omega = tt > 0 ? ts / tt : 0.0;
            op(2, PV, X, SV, float(alpha), float(omega));            // x += alpha p + omega s
            op(1, SV, R, TV, float(omega));                          // r = s - omega t
            rho = rho1;
            if (omega == 0 || !std::isfinite(omega)) break;
            if (it % 4 == 3) { double r2, d2; dot(R, R, R, R, r2, d2); if (r2 <= tol2) break; }
        }
        if (!ok) return false;
        std::vector<float> xs(static_cast<size_t>(n_));
        if (!download(X, xs.data())) return false;
        for (int e = 0; e < n_; ++e) { z[e] = double(xs[size_t(e)]); if (!std::isfinite(z[e])) return false; }
        return true;
    }

private:
    ID3D11Device* dev_ = nullptr;
    ID3D11DeviceContext* ctx_ = nullptr;
    ID3D11ComputeShader *csSpmv_ = nullptr, *csOp_ = nullptr, *csDot_ = nullptr;
    ID3D11Buffer* cb_ = nullptr;
    ID3D11Buffer *bRP_ = nullptr, *bCI_ = nullptr, *bBV_ = nullptr, *bV_ = nullptr, *bPart_ = nullptr, *stPart_ = nullptr, *stVec_ = nullptr;
    ID3D11ShaderResourceView *sRP_ = nullptr, *sCI_ = nullptr, *sBV_ = nullptr;
    ID3D11UnorderedAccessView *uV_ = nullptr, *uPart_ = nullptr;
    bool tried_ = false, okDev_ = false;
    int na_ = 0, nb_ = 0, n_ = 0, capN_ = -1, capNb_ = 0, maxIt_ = 30;
    size_t capNnz_ = 0;
    std::vector<double> dinv_;
    std::string adapter_ = "GPU (D3D11)";
    std::mutex claimM_;
    std::thread::id owner_;
    std::chrono::steady_clock::time_point lastUse_;

    static int groups(int n, int g) { return (n + g - 1) / g; }

    // one simulation thread at a time owns the device (studies run several simulators in parallel)
    bool claim() {
        std::lock_guard<std::mutex> lk(claimM_);
        const auto now = std::chrono::steady_clock::now(); const auto me = std::this_thread::get_id();
        if (owner_ != std::thread::id() && owner_ != me && now - lastUse_ < std::chrono::seconds(3)) return false;
        owner_ = me; lastUse_ = now;
        return true;
    }

    static bool invert(const double* A, double* Ai, int n) {
        double M[4][8] = {};
        for (int r = 0; r < n; ++r) { for (int c = 0; c < n; ++c) M[r][c] = A[r * n + c]; M[r][n + r] = 1; }
        for (int c = 0; c < n; ++c) {
            int piv = c; for (int r = c + 1; r < n; ++r) if (std::fabs(M[r][c]) > std::fabs(M[piv][c])) piv = r;
            if (std::fabs(M[piv][c]) < 1e-300) return false;
            if (piv != c) for (int k = 0; k < 2 * n; ++k) std::swap(M[c][k], M[piv][k]);
            const double inv = 1.0 / M[c][c];
            for (int k = 0; k < 2 * n; ++k) M[c][k] *= inv;
            for (int r = 0; r < n; ++r) if (r != c && M[r][c] != 0) { const double f = M[r][c]; for (int k = 0; k < 2 * n; ++k) M[r][k] -= f * M[c][k]; }
        }
        for (int r = 0; r < n; ++r) for (int c = 0; c < n; ++c) Ai[r * n + c] = M[r][n + c];
        return true;
    }

    bool ensureDevice() {
        if (tried_) return okDev_;
        tried_ = true;
        const D3D_FEATURE_LEVEL fl[] = { D3D_FEATURE_LEVEL_11_1, D3D_FEATURE_LEVEL_11_0 };
        D3D_FEATURE_LEVEL got;
        HRESULT hr = D3D11CreateDevice(nullptr, D3D_DRIVER_TYPE_HARDWARE, nullptr, D3D11_CREATE_DEVICE_SINGLETHREADED, fl, 2,
                                       D3D11_SDK_VERSION, &dev_, &got, &ctx_);
        if (FAILED(hr)) hr = D3D11CreateDevice(nullptr, D3D_DRIVER_TYPE_HARDWARE, nullptr, D3D11_CREATE_DEVICE_SINGLETHREADED, fl + 1, 1,
                                               D3D11_SDK_VERSION, &dev_, &got, &ctx_);
        if (FAILED(hr) || !dev_) return false;
        IDXGIDevice* dx = nullptr;
        if (SUCCEEDED(dev_->QueryInterface(__uuidof(IDXGIDevice), reinterpret_cast<void**>(&dx)))) {
            IDXGIAdapter* ad = nullptr;
            if (SUCCEEDED(dx->GetAdapter(&ad))) {
                DXGI_ADAPTER_DESC d; if (SUCCEEDED(ad->GetDesc(&d))) {
                    char buf[256]; WideCharToMultiByte(CP_UTF8, 0, d.Description, -1, buf, sizeof buf, nullptr, nullptr);
                    adapter_ = std::string("GPU: ") + buf;
                }
                ad->Release();
            }
            dx->Release();
        }
        auto compile = [&](const char* entry, ID3D11ComputeShader*& cs) {
            ID3DBlob *code = nullptr, *err = nullptr;
            HRESULT h = D3DCompile(kCS, std::strlen(kCS), "gpu_solver.hlsl", nullptr, nullptr, entry, "cs_5_0", D3DCOMPILE_OPTIMIZATION_LEVEL3, 0, &code, &err);
            if (err) { OutputDebugStringA(static_cast<const char*>(err->GetBufferPointer())); err->Release(); }
            if (FAILED(h) || !code) return false;
            h = dev_->CreateComputeShader(code->GetBufferPointer(), code->GetBufferSize(), nullptr, &cs);
            code->Release();
            return SUCCEEDED(h);
        };
        if (!compile("CSSpmv", csSpmv_) || !compile("CSOp", csOp_) || !compile("CSDot", csDot_)) { freeDevice(); return false; }
        D3D11_BUFFER_DESC bd = {}; bd.ByteWidth = sizeof(Params); bd.Usage = D3D11_USAGE_DEFAULT; bd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
        if (FAILED(dev_->CreateBuffer(&bd, nullptr, &cb_))) { freeDevice(); return false; }
        okDev_ = true;
        return true;
    }

    bool makeSrv(const void* data, UINT stride, size_t count, ID3D11Buffer*& b, ID3D11ShaderResourceView*& v) {
        D3D11_BUFFER_DESC bd = {}; bd.ByteWidth = UINT(stride * count); bd.Usage = D3D11_USAGE_DEFAULT;
        bd.BindFlags = D3D11_BIND_SHADER_RESOURCE; bd.MiscFlags = D3D11_RESOURCE_MISC_BUFFER_STRUCTURED; bd.StructureByteStride = stride;
        D3D11_SUBRESOURCE_DATA sd = {}; sd.pSysMem = data;
        if (FAILED(dev_->CreateBuffer(&bd, &sd, &b))) return false;
        D3D11_SHADER_RESOURCE_VIEW_DESC vd = {}; vd.Format = DXGI_FORMAT_UNKNOWN; vd.ViewDimension = D3D11_SRV_DIMENSION_BUFFER;
        vd.Buffer.FirstElement = 0; vd.Buffer.NumElements = UINT(count);
        return SUCCEEDED(dev_->CreateShaderResourceView(b, &vd, &v));
    }
    bool makeUav(UINT stride, size_t count, ID3D11Buffer*& b, ID3D11UnorderedAccessView*& v) {
        D3D11_BUFFER_DESC bd = {}; bd.ByteWidth = UINT(stride * count); bd.Usage = D3D11_USAGE_DEFAULT;
        bd.BindFlags = D3D11_BIND_UNORDERED_ACCESS; bd.MiscFlags = D3D11_RESOURCE_MISC_BUFFER_STRUCTURED; bd.StructureByteStride = stride;
        if (FAILED(dev_->CreateBuffer(&bd, nullptr, &b))) return false;
        D3D11_UNORDERED_ACCESS_VIEW_DESC vd = {}; vd.Format = DXGI_FORMAT_UNKNOWN; vd.ViewDimension = D3D11_UAV_DIMENSION_BUFFER;
        vd.Buffer.FirstElement = 0; vd.Buffer.NumElements = UINT(count);
        return SUCCEEDED(dev_->CreateUnorderedAccessView(b, &vd, &v));
    }
    bool makeStaging() {
        D3D11_BUFFER_DESC bd = {}; bd.Usage = D3D11_USAGE_STAGING; bd.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
        bd.ByteWidth = UINT(sizeof(float) * 2 * size_t(groups(n_, 256)));
        if (FAILED(dev_->CreateBuffer(&bd, nullptr, &stPart_))) return false;
        bd.ByteWidth = UINT(sizeof(float) * size_t(n_));
        return SUCCEEDED(dev_->CreateBuffer(&bd, nullptr, &stVec_));
    }

    void bind(ID3D11ComputeShader* cs, const Params& p) {
        ctx_->UpdateSubresource(cb_, 0, nullptr, &p, 0, 0);
        ctx_->CSSetShader(cs, nullptr, 0);
        ctx_->CSSetConstantBuffers(0, 1, &cb_);
        ID3D11ShaderResourceView* srv[3] = { sRP_, sCI_, sBV_ };
        ctx_->CSSetShaderResources(0, 3, srv);
        ID3D11UnorderedAccessView* uav[2] = { uV_, uPart_ };
        ctx_->CSSetUnorderedAccessViews(0, 2, uav, nullptr);
    }
    Params base() const { Params p = {}; p.n = UINT(n_); p.nb = UINT(nb_); p.na = UINT(na_); return p; }
    UINT off(int slot) const { return UINT(slot) * UINT(n_); }

    void op(int mode, int x, int y, int z, float a = 0, float b = 0) {
        Params p = base(); p.mode = UINT(mode); p.oX = off(x); p.oY = off(y); p.oZ = off(z); p.a = a; p.b = b;
        bind(csOp_, p); ctx_->Dispatch(UINT(groups(n_, 256)), 1, 1);
    }
    void spmv(int x, int y) {
        Params p = base(); p.oX = off(x); p.oY = off(y);
        bind(csSpmv_, p); ctx_->Dispatch(UINT(groups(na_, 128)), 1, 1);
    }
    void dot(int a, int b, int c, int d, double& r1, double& r2) {
        Params p = base(); p.oX = off(a); p.oY = off(b); p.oZ = off(c); p.oW = off(d);
        bind(csDot_, p);
        const int g = groups(n_, 256);
        ctx_->Dispatch(UINT(g), 1, 1);
        ctx_->CopyResource(stPart_, bPart_);
        D3D11_MAPPED_SUBRESOURCE m;
        r1 = r2 = std::nan("");
        if (FAILED(ctx_->Map(stPart_, 0, D3D11_MAP_READ, 0, &m))) return;
        const float* f = static_cast<const float*>(m.pData);
        double s1 = 0, s2 = 0;
        for (int i = 0; i < g; ++i) { s1 += double(f[2 * i]); s2 += double(f[2 * i + 1]); }
        ctx_->Unmap(stPart_, 0);
        r1 = s1; r2 = s2;
    }
    void upload(int slot, const float* data) {
        D3D11_BOX box = { UINT(off(slot) * sizeof(float)), 0, 0, UINT((off(slot) + UINT(n_)) * sizeof(float)), 1, 1 };
        ctx_->UpdateSubresource(bV_, 0, &box, data, 0, 0);
    }
    bool download(int slot, float* data) {
        D3D11_BOX box = { UINT(off(slot) * sizeof(float)), 0, 0, UINT((off(slot) + UINT(n_)) * sizeof(float)), 1, 1 };
        ctx_->CopySubresourceRegion(stVec_, 0, 0, 0, 0, bV_, 0, &box);
        D3D11_MAPPED_SUBRESOURCE m;
        if (FAILED(ctx_->Map(stVec_, 0, D3D11_MAP_READ, 0, &m))) return false;
        std::memcpy(data, m.pData, sizeof(float) * size_t(n_));
        ctx_->Unmap(stVec_, 0);
        return true;
    }

    void freeMatrix() {
        rel(sRP_); rel(sCI_); rel(sBV_); rel(uV_); rel(uPart_);
        rel(bRP_); rel(bCI_); rel(bBV_); rel(bV_); rel(bPart_); rel(stPart_); rel(stVec_);
        capN_ = -1; capNnz_ = 0;
    }
    void freeDevice() { rel(csSpmv_); rel(csOp_); rel(csDot_); rel(cb_); rel(ctx_); rel(dev_); okDev_ = false; }
};

D3D11Accelerator* g_gpu = nullptr;

} // namespace

void registerGpuSolver() {
    if (!g_gpu) g_gpu = new D3D11Accelerator();       // lives for the process; device created lazily on first use
    ressim::setBlockAccelerator(g_gpu);
}
