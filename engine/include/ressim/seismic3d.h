// C36 seismic for the 3-D geomodel scene:
//   * SegyVolume - out-of-core SEG-Y: only the trace headers are scanned (inline / crossline / CDP X / CDP Y, plan
//     geometry by a least-squares fit over all traces); traces, inline / crossline / time slices and sub-cubes are read
//     on demand from the file through a small trace cache, so volumes far larger than the memory can be browsed.
//     Sample formats 1 (IBM float), 2 (int32), 3 (int16), 5 (IEEE float), 8 (int8); big-endian (standard) and
//     byte-swapped little-endian files (detected from the format code); missing traces read as zeros.
//   * proxy(): a decimated in-memory SeismicCube (stride per axis) for volume rendering / quick looks, built in one
//     sequential pass with a progress callback.
//   * Trace attributes (complex-trace analysis, Taner et al. 1979): envelope, instantaneous phase / frequency, cosine
//     of phase, sweetness, RMS amplitude, relative acoustic impedance; multi-trace: semblance coherence (Marfurt et al.
//     1998) and variance - per trace window, on slices or whole volumes (streamed inline by inline to a new SEG-Y file).
//   * Flattening of a section on a picked horizon (each trace shifted so the horizon is flat at a reference time).
#pragma once
#include "ressim/geomodel.h"
#include <cstdint>
#include <cstdio>
#include <functional>
#include <list>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace ressim {
namespace geo {

class SegyVolume {
public:
    SegyVolume() = default;
    ~SegyVolume();
    SegyVolume(const SegyVolume&) = delete;
    SegyVolume& operator=(const SegyVolume&) = delete;

    bool open(const std::string& path, std::string& err, const std::function<bool(double)>& progress = {});
    void close();
    bool isOpen() const { return f_ != nullptr; }

    // geometry (same conventions as SeismicCube)
    std::string path, name, textHeader;
    int ni = 0, nj = 0, ns = 0;
    int il0 = 1, xl0 = 1, dil = 1, dxl = 1;
    double t0 = 0, dt = 4;
    double x0 = 0, y0 = 0, ix = 1, iy = 0, jx = 0, jy = 1;
    int format = 5;                    // sample format code of the file
    bool littleEndian = false;
    long long nTraces = 0, missing = 0;   // traces in the file / lattice positions without a trace
    double fitError = 0;               // RMS of the plan-geometry fit over all traces (coordinate units)
    float ampAbsMax = 0, ampRms = 0;   // from a sparse scan (every k-th trace)
    long long fileBytes = 0;
    // C42 trace-header byte locations (1-based; set before open())
    int byteIl = 189, byteXl = 193, byteX = 181, byteY = 185, byteScalar = 71;

    // reading (thread-safe; false / zeros where the lattice has no trace)
    bool trace(int i, int j, float* out) const;
    std::vector<float> inlineSlice(int i) const;           // [j][s]  (nj x ns)
    std::vector<float> xlineSlice(int j) const;            // [i][s]  (ni x ns)
    std::vector<float> timeSlice(int s) const;             // [i][j]  (ni x nj) - one sample per trace
    SeismicCube subCube(int i0, int i1, int j0, int j1, int s0, int s1) const;   // inclusive ranges, in memory
    // decimated cube: the largest axis gets at most maxDim samples (same stride on the plan axes, separate stride in time)
    SeismicCube proxy(int maxDimPlan, int maxDimTime, const std::function<bool(double)>& progress = {}) const;
    SeismicCube toCube() const { return subCube(0, ni - 1, 0, nj - 1, 0, ns - 1); }   // whole volume (small files only)
    void xy(double i, double j, double& x, double& y) const { x = x0 + i * ix + j * jx; y = y0 + i * iy + j * jy; }
    bool ij(double x, double y, double& i, double& j) const;
    size_t cacheTraces = 4096;         // trace cache capacity

private:
    std::FILE* f_ = nullptr;
    long long start_ = 3600, tlen_ = 0;
    int bps_ = 4;                      // bytes per sample
    std::vector<long long> off_;       // data offset per lattice position (i*nj+j), -1 = missing
    mutable std::unordered_map<long long, std::list<std::pair<long long, std::vector<float>>>::iterator> cmap_;
    mutable std::list<std::pair<long long, std::vector<float>>> lru_;
    mutable void* mtx_ = nullptr;      // std::mutex (kept opaque in the header)
    void decode(const unsigned char* raw, float* out) const;
    bool readRaw(long long off, long long bytes, unsigned char* dst) const;
};

// ---------------------------------------------------------------- attributes
enum class SeisAttr { Amplitude = 0, Envelope, Phase, Frequency, CosPhase, Sweetness, Rms, RelImpedance, Coherence, Variance, Count };
const char* seisAttrName(SeisAttr a);
const char* seisAttrUnit(SeisAttr a);   // "", "deg", "Hz", ...
bool seisAttrMultiTrace(SeisAttr a);    // needs neighbouring traces (coherence, variance)

// analytic signal of one trace (FFT Hilbert transform, zero-padded to a power of two): quadrature 'q'
void hilbert(const float* x, int n, std::vector<float>& q);
// single-trace attribute (dtMs = sample interval, winSamples = half window for RMS / sweetness smoothing)
void traceAttribute(SeisAttr a, const float* x, int n, double dtMs, int winSamples, float* out);
// attribute of a vertical section: 'sec' [tr][s] (ntr x ns); for multi-trace attributes the traces must be ordered
// along the section; 'nbrPrev' / 'nbrNext' optional parallel sections (3-D 3x3 coherence), same layout
std::vector<float> sectionAttribute(SeisAttr a, const std::vector<float>& sec, int ntr, int ns, double dtMs, int winSamples,
                                    const std::vector<float>* nbrPrev = nullptr, const std::vector<float>* nbrNext = nullptr);
// attribute cube of an in-memory volume (OpenMP over inlines)
SeismicCube cubeAttribute(const SeismicCube& c, SeisAttr a, int winSamples);
// streamed attribute volume of an out-of-core SEG-Y written to 'outPath' (IEEE float, same headers / geometry)
bool volumeAttribute(const SegyVolume& v, SeisAttr a, int winSamples, const std::string& outPath, std::string& err,
                     const std::function<bool(double)>& progress = {});
// semblance of K traces over samples [s0, s1]: (sum_t (sum_k u_k)^2) / (K sum_t sum_k u_k^2), in [0, 1]
double semblance(const std::vector<const float*>& tr, int s0, int s1);

// ---------------------------------------------------------------- flattening
// shift every trace of a section [tr][s] so that sample position hzSample[tr] moves to refSample (linear
// interpolation, NaN horizon = trace left unshifted); returns the flattened section (same size, outside -> 0)
std::vector<float> flattenSection(const std::vector<float>& sec, int ntr, int ns, const std::vector<double>& hzSample, double refSample);

// ---------------------------------------------------------------- helpers for the 3-D scene
// robust amplitude clip: the |a| percentile (e.g. 0.99) of a sample of the values
float amplitudeClip(const std::vector<float>& v, double pct);
// writes a synthetic SEG-Y (faulted anticline reflectivity convolved with a Ricker wavelet) - tests / demo
bool writeSyntheticSegy(const std::string& path, int ni, int nj, int ns, double dtMs, int format, bool littleEndian, std::string& err);

}  // namespace geo
}  // namespace ressim
