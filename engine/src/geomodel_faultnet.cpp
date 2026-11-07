// C41 (Petrel parity 5): fault network - fault-fault junctions (X / Y / lambda), truncation of minor faults at their
// major fault on every depth level (new terminal stick from the trace intersections), and horizon-fault throw profiles.
// See geomodel.h for the conventions. Used by the Geomodel Structure step before buildGrid: a minor fault that ends
// short of its major leaves the fault block open (both sides stay one connected block); after truncation the tip
// lies on the major surface and the block closes.
#include "ressim/geomodel.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace ressim {
namespace geo {
namespace {

using Pt = std::array<double, 2>;
using Line = std::vector<Pt>;

double cross2(double ax, double ay, double bx, double by) { return ax * by - ay * bx; }

// segment a-b vs c-d: true with the parameters t (on a-b) and u (on c-d) in [0, 1]
bool segHit(const Pt& a, const Pt& b, const Pt& c, const Pt& d, double& t, double& u) {
    const double rx = b[0] - a[0], ry = b[1] - a[1], sx = d[0] - c[0], sy = d[1] - c[1];
    const double den = cross2(rx, ry, sx, sy);
    if (std::fabs(den) < 1e-12 * (std::hypot(rx, ry) * std::hypot(sx, sy) + 1e-300)) return false;
    const double qx = c[0] - a[0], qy = c[1] - a[1];
    t = cross2(qx, qy, sx, sy) / den;
    u = cross2(qx, qy, rx, ry) / den;
    return t >= -1e-12 && t <= 1 + 1e-12 && u >= -1e-12 && u <= 1 + 1e-12;
}

double length(const Line& L) { double s = 0; for (size_t k = 1; k < L.size(); ++k) s += std::hypot(L[k][0] - L[k - 1][0], L[k][1] - L[k - 1][1]); return s; }

// nearest point on the polyline: distance, signed side (+1 left of the direction of travel, -1 right), arc position
struct Near { double d = 1e300, s = 0; int side = 0; Pt p{0, 0}; size_t seg = 0; };
Near nearest(const Line& L, const Pt& q) {
    Near r;
    double acc = 0;
    for (size_t k = 0; k + 1 < L.size(); ++k) {
        const double ax = L[k][0], ay = L[k][1], bx = L[k + 1][0] - ax, by = L[k + 1][1] - ay;
        const double L2 = bx * bx + by * by, Ls = std::sqrt(L2);
        const double t = L2 > 0 ? std::clamp(((q[0] - ax) * bx + (q[1] - ay) * by) / L2, 0.0, 1.0) : 0.0;
        const double px = ax + t * bx, py = ay + t * by, d = std::hypot(q[0] - px, q[1] - py);
        if (d < r.d) {
            r.d = d; r.s = acc + t * Ls; r.p = {px, py}; r.seg = k;
            const double c = cross2(bx, by, q[0] - ax, q[1] - ay);
            r.side = c > 0 ? 1 : c < 0 ? -1 : 0;
        }
        acc += Ls;
    }
    return r;
}

// first crossing of segment a-b with the polyline M (smallest t)
bool hitLine(const Pt& a, const Pt& b, const Line& M, double& tBest, Pt& p) {
    bool any = false; tBest = 2;
    for (size_t k = 0; k + 1 < M.size(); ++k) {
        double t, u;
        if (segHit(a, b, M[k], M[k + 1], t, u) && t < tBest) { tBest = t; any = true; p = {a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])}; }
    }
    return any;
}

void depthSpan(const Fault& f, double& zt, double& zb) {
    zt = 1e300; zb = -1e300;
    for (const auto& s : f.sticks) for (const auto& p : s) { zt = std::min(zt, p[2]); zb = std::max(zb, p[2]); }
}

}  // namespace

std::vector<FaultJunction> detectJunctions(const std::vector<Fault>& faults, double zRef, double tol) {
    std::vector<FaultJunction> out;
    std::vector<Line> T(faults.size());
    for (size_t f = 0; f < faults.size(); ++f) T[f] = faults[f].traceAt(zRef);
    for (size_t i = 0; i < faults.size(); ++i)
        for (size_t j = i + 1; j < faults.size(); ++j) {
            const Line &A = T[i], &B = T[j];
            if (A.size() < 2 || B.size() < 2) continue;
            const double LA = length(A), LB = length(B);
            // 1) crossings
            bool crossed = false;
            double accA = 0;
            for (size_t a = 0; a + 1 < A.size() && !crossed; ++a) {
                const double la = std::hypot(A[a + 1][0] - A[a][0], A[a + 1][1] - A[a][1]);
                double accB = 0;
                for (size_t b = 0; b + 1 < B.size(); ++b) {
                    const double lb = std::hypot(B[b + 1][0] - B[b][0], B[b + 1][1] - B[b][1]);
                    double t, u;
                    if (segHit(A[a], A[a + 1], B[b], B[b + 1], t, u)) {
                        const double sA = accA + t * la, sB = accB + u * lb;
                        const double restA = std::min(sA, LA - sA), restB = std::min(sB, LB - sB);
                        FaultJunction J;
                        J.x = A[a][0] + t * (A[a + 1][0] - A[a][0]); J.y = A[a][1] + t * (A[a + 1][1] - A[a][1]);
                        if (restA > tol && restB > tol) { J.type = 'X'; J.major = int(LA >= LB ? i : j); J.minor = int(LA >= LB ? j : i); J.truncate = false; }
                        else {
                            J.type = 'Y';
                            const bool aEnds = restA <= tol, bEnds = restB <= tol;
                            const bool aMinor = aEnds && (!bEnds || LA <= LB);
                            J.minor = int(aMinor ? i : j); J.major = int(aMinor ? j : i);
                        }
                        out.push_back(J);
                        crossed = true;
                        break;
                    }
                    accB += lb;
                }
                accA += la;
            }
            if (crossed) continue;
            // 2) tips short of the other trace: extend the end segment by tol
            FaultJunction best; best.gap = 1e300;
            auto tryTips = [&](const Line& M, const Line& O, int minor, int major) {
                for (int end = 0; end < 2; ++end) {
                    const Pt tip = end ? M.back() : M.front(), prev = end ? M[M.size() - 2] : M[1];
                    const double dx = tip[0] - prev[0], dy = tip[1] - prev[1], L = std::hypot(dx, dy);
                    if (L <= 0) continue;
                    const Pt far{tip[0] + dx / L * tol, tip[1] + dy / L * tol};
                    double t; Pt p;
                    if (hitLine(tip, far, O, t, p) && t * tol < best.gap) {
                        best.gap = t * tol; best.minor = minor; best.major = major; best.x = p[0]; best.y = p[1]; best.type = 'Y';
                    }
                }
            };
            tryTips(A, B, int(i), int(j));
            tryTips(B, A, int(j), int(i));
            if (best.minor >= 0) out.push_back(best);
        }
    return out;
}

std::vector<Fault> truncateFaults(const std::vector<Fault>& faults, const std::vector<FaultJunction>& rules, double tol, FaultNetReport* rep) {
    std::vector<Fault> out = faults;
    FaultNetReport R;
    for (const FaultJunction& J : rules) {
        if (!J.truncate || J.minor < 0 || J.major < 0 || J.minor == J.major || size_t(J.minor) >= out.size() || size_t(J.major) >= out.size()) continue;
        const Fault minorF = out[size_t(J.minor)];   // copy: geometry before this rule
        const Fault& majorF = out[size_t(J.major)];
        std::vector<std::vector<P3>> sticks;
        for (const auto& s : minorF.sticks) if (!s.empty()) sticks.push_back(s);
        if (sticks.size() < 2) { R.notes.push_back(minorF.name + ": fewer than two sticks - not truncated"); continue; }
        double zt, zb; depthSpan(minorF, zt, zb);
        if (!(zb > zt)) { R.notes.push_back(minorF.name + ": no depth extent - not truncated"); continue; }
        const double zm = 0.5 * (zt + zb);
        Fault mf = minorF; mf.sticks = sticks;
        const Line trm = mf.traceAt(zm), M = majorF.traceAt(zm);
        if (M.size() < 2) { R.notes.push_back(majorF.name + ": no trace at " + std::to_string(zm) + " - " + minorF.name + " not truncated"); continue; }
        // side of every stick (at mid depth) relative to the major
        std::vector<int> sd(trm.size());
        int np = 0, nm = 0;
        for (size_t k = 0; k < trm.size(); ++k) {
            const Near n = nearest(M, trm[k]);
            sd[k] = n.d < 1e-9 * (1 + tol) ? 0 : n.side;
            np += sd[k] > 0; nm += sd[k] < 0;
        }
        const int keep = J.keepSide != 0 ? (J.keepSide > 0 ? 1 : -1) : (np >= nm ? 1 : -1);
        // contiguous run of kept sticks
        int a = -1, b = -1;
        for (int k = 0; k < int(sd.size()); ++k) if (sd[size_t(k)] == keep || sd[size_t(k)] == 0) { if (a < 0) a = k; b = k; }
        if (a < 0) { R.notes.push_back(minorF.name + ": no stick on the kept side of " + majorF.name); continue; }
        for (int k = a; k <= b; ++k)
            if (sd[size_t(k)] == -keep) { R.notes.push_back(minorF.name + " crosses " + majorF.name + " more than once - not truncated"); a = -1; break; }
        if (a < 0) continue;
        bool clipFront;
        if (a > 0 && b < int(sd.size()) - 1) { R.notes.push_back(minorF.name + ": both ends lie beyond " + majorF.name + " - not truncated"); continue; }
        if (a > 0) clipFront = true;
        else if (b < int(sd.size()) - 1) clipFront = false;
        else clipFront = nearest(M, trm.front()).d < nearest(M, trm.back()).d;   // gap case: the end nearer the major
        // terminal stick: intersection of the minor's trace (from the last kept stick outward, then extended by tol)
        // with the major's trace at a series of depths
        const int nLev = 9;
        std::vector<P3> term;
        for (int l = 0; l < nLev; ++l) {
            const double z = zt + (zb - zt) * double(l) / double(nLev - 1);
            const Line tr = mf.traceAt(z), Mz = majorF.traceAt(z);
            if (tr.size() != sticks.size() || Mz.size() < 2) continue;
            Line path;   // from the boundary kept stick toward the clipped tip, then the extension
            if (clipFront) for (int k = a; k >= 0; --k) path.push_back(tr[size_t(k)]);
            else for (int k = b; k < int(tr.size()); ++k) path.push_back(tr[size_t(k)]);
            if (path.size() == 1) path.insert(path.begin(), clipFront ? tr[1] : tr[tr.size() - 2]);   // direction from the neighbour stick
            const Pt tip = path.back(), prev = path[path.size() - 2];
            const double dx = tip[0] - prev[0], dy = tip[1] - prev[1], L = std::hypot(dx, dy);
            if (L > 0) path.push_back({tip[0] + dx / L * tol, tip[1] + dy / L * tol});
            for (size_t k = 0; k + 1 < path.size(); ++k) {
                double t; Pt p;
                if (hitLine(path[k], path[k + 1], Mz, t, p)) { term.push_back({p[0], p[1], z}); break; }
            }
        }
        if (term.size() < 2) { R.notes.push_back(minorF.name + ": does not reach " + majorF.name + " within the tolerance - not truncated"); continue; }
        std::vector<std::vector<P3>> ns;
        if (clipFront) { ns.push_back(term); for (int k = a; k <= b; ++k) ns.push_back(sticks[size_t(k)]); }
        else { for (int k = a; k <= b; ++k) ns.push_back(sticks[size_t(k)]); ns.push_back(term); }
        R.sticksRemoved += int(sticks.size()) - (b - a + 1);
        R.sticksAdded += 1;
        ++R.truncated;
        out[size_t(J.minor)].sticks = ns;
        R.notes.push_back(minorF.name + " truncated at " + majorF.name + " (" + std::to_string(int(sticks.size()) - (b - a + 1)) + " stick(s) removed, terminal stick with " +
                          std::to_string(term.size()) + " points)");
    }
    if (rep) *rep = R;
    return out;
}

ThrowProfile faultThrowProfile(const Fault& f, const Surface& h, double offset, int n) {
    ThrowProfile P;
    if (!h.valid() || n < 2) return P;
    double z0 = 0; int c = 0;
    for (double v : h.z) if (std::isfinite(v)) { z0 += v; ++c; }
    if (!c) return P;
    z0 /= c;
    const Line tr = f.traceAt(z0), deep = f.traceAt(z0 + 50), shal = f.traceAt(z0 - 50);
    if (tr.size() < 2) return P;
    const double L = length(tr);
    if (!(L > 0)) return P;
    double prevSign = 0;
    for (int q = 0; q < n; ++q) {
        const double s = L * double(q) / double(n - 1);
        // point and tangent at arc position s
        double acc = 0; size_t k = 0;
        while (k + 2 < tr.size() && acc + std::hypot(tr[k + 1][0] - tr[k][0], tr[k + 1][1] - tr[k][1]) < s) { acc += std::hypot(tr[k + 1][0] - tr[k][0], tr[k + 1][1] - tr[k][1]); ++k; }
        const double sl = std::hypot(tr[k + 1][0] - tr[k][0], tr[k + 1][1] - tr[k][1]);
        const double t = sl > 0 ? std::clamp((s - acc) / sl, 0.0, 1.0) : 0.0;
        const double px = tr[k][0] + t * (tr[k + 1][0] - tr[k][0]), py = tr[k][1] + t * (tr[k + 1][1] - tr[k][1]);
        const double tx = sl > 0 ? (tr[k + 1][0] - tr[k][0]) / sl : 1, ty = sl > 0 ? (tr[k + 1][1] - tr[k][1]) / sl : 0;
        const double nx = -ty, ny = tx;   // left normal
        double hw = 1;
        if (deep.size() == tr.size() && shal.size() == tr.size()) {
            const double dxv = (1 - t) * (deep[k][0] - shal[k][0]) + t * (deep[k + 1][0] - shal[k + 1][0]);
            const double dyv = (1 - t) * (deep[k][1] - shal[k][1]) + t * (deep[k + 1][1] - shal[k + 1][1]);
            const double dd = dxv * nx + dyv * ny;
            if (std::fabs(dd) > 1e-9) hw = dd > 0 ? 1 : -1;   // the fault dips toward the hanging wall
        }
        const double zh = h.at(px + hw * nx * offset, py + hw * ny * offset), zf = h.at(px - hw * nx * offset, py - hw * ny * offset);
        if (!std::isfinite(zh) || !std::isfinite(zf)) continue;
        const double th = zh - zf;
        P.s.push_back(s); P.throwV.push_back(th); P.x.push_back(px); P.y.push_back(py);
        P.maxThrow = std::max(P.maxThrow, std::fabs(th));
    }
    for (double th : P.throwV) {
        if (std::fabs(th) < 0.05 * P.maxThrow) continue;
        const double sg = th > 0 ? 1 : -1;
        if (prevSign != 0 && sg != prevSign) P.reverses = true;
        prevSign = sg;
    }
    if (P.throwV.size() >= 3) P.tipsTaper = std::fabs(P.throwV.front()) < 0.5 * P.maxThrow && std::fabs(P.throwV.back()) < 0.5 * P.maxThrow;
    return P;
}

}  // namespace geo
}  // namespace ressim
