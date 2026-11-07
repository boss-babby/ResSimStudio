// C14: little binary serialisation kit used by the results / restart files (.rsres) and the engine restart states.
//   BinWriter / BinReader   append / consume PODs, strings and vectors (little-endian host layout, size-prefixed)
//   StateOut / StateIn      symmetric "archive" visitors: one templated stateIO(ar) per engine lists every dynamic
//                           member once, and the same list is used for writing and reading (no drift between them)
#pragma once
#include <cstdint>
#include <cstring>
#include <functional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace ressim {

struct BinWriter {
    std::string buf;
    void raw(const void* p, size_t n) { buf.append(static_cast<const char*>(p), n); }
    template <class T> void pod(const T& v) { static_assert(std::is_trivially_copyable<T>::value, "POD only"); raw(&v, sizeof(T)); }
    void u64(uint64_t v) { pod(v); }
    void str(const std::string& s) { u64(s.size()); raw(s.data(), s.size()); }
    template <class T> void vec(const std::vector<T>& v) {
        static_assert(std::is_trivially_copyable<T>::value, "POD vector only");
        u64(v.size()); if (!v.empty()) raw(v.data(), v.size() * sizeof(T));
    }
    void vecBool(const std::vector<bool>& v) { u64(v.size()); for (bool b : v) { const uint8_t c = b ? 1 : 0; pod(c); } }
    void strs(const std::vector<std::string>& v) { u64(v.size()); for (const auto& s : v) str(s); }
};

struct BinReader {
    const char* p = nullptr;
    const char* end = nullptr;
    BinReader() = default;
    BinReader(const char* b, size_t n) : p(b), end(b + n) {}
    explicit BinReader(const std::string& s) : p(s.data()), end(s.data() + s.size()) {}
    size_t left() const { return size_t(end - p); }
    void need(size_t n) const { if (size_t(end - p) < n) throw std::runtime_error("truncated or corrupt binary data"); }
    void raw(void* d, size_t n) { need(n); std::memcpy(d, p, n); p += n; }
    template <class T> void pod(T& v) { raw(&v, sizeof(T)); }
    template <class T> T get() { T v{}; pod(v); return v; }
    uint64_t u64() { return get<uint64_t>(); }
    std::string str() { const uint64_t n = u64(); need(size_t(n)); std::string s(p, size_t(n)); p += n; return s; }
    template <class T> void vec(std::vector<T>& v) {
        const uint64_t n = u64();
        if (n > uint64_t(left() / (sizeof(T) ? sizeof(T) : 1))) throw std::runtime_error("corrupt vector length in binary data");
        v.resize(size_t(n)); if (n) raw(v.data(), size_t(n) * sizeof(T));
    }
    void vecBool(std::vector<bool>& v) { const uint64_t n = u64(); need(size_t(n)); v.resize(size_t(n)); for (size_t i = 0; i < n; ++i) { uint8_t c = 0; pod(c); v[i] = c != 0; } }
    void strs(std::vector<std::string>& v) { const uint64_t n = u64(); if (n > left()) throw std::runtime_error("corrupt string list"); v.resize(size_t(n)); for (auto& s : v) s = str(); }
};

// ---- archive visitors: ar(a, b, c, ...) writes or reads every argument in order
struct StateOut {
    BinWriter w;
    static constexpr bool reading = false;
    template <class T> void one(const T& v) {
        if constexpr (std::is_trivially_copyable<T>::value) w.pod(v);
        else static_assert(std::is_trivially_copyable<T>::value, "add an overload");
    }
    template <class T> void one(const std::vector<T>& v) {
        if constexpr (std::is_trivially_copyable<T>::value) w.vec(v);
        else { w.u64(v.size()); for (const auto& e : v) one(e); }
    }
    void one(const std::vector<bool>& v) { w.vecBool(v); }
    void one(const std::string& s) { w.str(s); }
    template <class... A> void operator()(const A&... a) { (one(a), ...); }
};

struct StateIn {
    BinReader r;
    static constexpr bool reading = true;
    explicit StateIn(const std::string& s) : r(s) {}
    template <class T> void one(T& v) {
        if constexpr (std::is_trivially_copyable<T>::value) r.pod(v);
        else static_assert(std::is_trivially_copyable<T>::value, "add an overload");
    }
    template <class T> void one(std::vector<T>& v) {
        if constexpr (std::is_trivially_copyable<T>::value) r.vec(v);
        else { const uint64_t n = r.u64(); if (n > r.left()) throw std::runtime_error("corrupt state"); v.resize(size_t(n)); for (auto& e : v) one(e); }
    }
    void one(std::vector<bool>& v) { r.vecBool(v); }
    void one(std::string& s) { s = r.str(); }
    template <class... A> void operator()(A&... a) { (one(a), ...); }
};


// Restart control of the deck runners (comp::runDeckModel / thermal::runDeckModel): save the complete run state at the
// end of chosen report steps and / or continue a run from such a state. The blob is self-contained (engine state +
// runner bookkeeping) and is stored in the RESTART chunks of the .rsres file.
struct DeckRestart {
    int startStep = 0;                  // > 0: continue after report step startStep (1-based count of completed steps)
    std::string state;                  // blob saved by onSave at that step
    int saveEvery = 0;                  // save every N report steps (0 = only saveSteps)
    std::vector<int> saveSteps;         // explicit report steps to save
    bool saveLast = false;              // always save the final report step
    std::function<void(int step, double days, const std::string& blob)> onSave;
    bool wants(int step, bool last) const {
        if (!onSave) return false;
        if (saveLast && last) return true;
        if (saveEvery > 0 && step % saveEvery == 0) return true;
        for (int s : saveSteps) if (s == step) return true;
        return false;
    }
};

} // namespace ressim
