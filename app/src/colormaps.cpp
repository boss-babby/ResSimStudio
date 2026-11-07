// Colour maps - see colormaps.h. Turbo / Viridis: polynomial fits (Google Turbo Apache-2.0; viridis fit by M. Zucker);
// C30 maps: piecewise-linear through the matplotlib / ColorBrewer control colours (perceptually uniform maps sampled at
// 9 points, ColorBrewer diverging maps at their 11 class colours).
#include "colormaps.h"
#include <algorithm>
#include <cmath>

namespace {
struct Stops { int n; const unsigned char (*c)[3]; };
const unsigned char kPlasma[9][3] = {{13, 8, 135}, {75, 3, 161}, {125, 3, 168}, {168, 34, 150}, {203, 70, 121}, {229, 107, 93}, {248, 148, 65}, {253, 195, 40}, {240, 249, 33}};
const unsigned char kInferno[9][3] = {{0, 0, 4}, {31, 12, 72}, {85, 15, 109}, {136, 34, 106}, {186, 54, 85}, {227, 89, 51}, {249, 140, 10}, {249, 201, 50}, {252, 255, 164}};
const unsigned char kMagma[9][3] = {{0, 0, 4}, {28, 16, 68}, {79, 18, 123}, {129, 37, 129}, {181, 54, 122}, {229, 80, 100}, {251, 135, 97}, {254, 194, 135}, {252, 253, 191}};
const unsigned char kCividis[9][3] = {{0, 34, 78}, {35, 61, 110}, {67, 78, 108}, {93, 96, 110}, {120, 118, 116}, {147, 141, 120}, {177, 165, 115}, {210, 192, 99}, {254, 232, 56}};
const unsigned char kRainbow[9][3] = {{110, 0, 220}, {40, 40, 255}, {0, 150, 255}, {0, 220, 200}, {60, 230, 60}, {200, 240, 20}, {255, 200, 0}, {255, 110, 0}, {220, 0, 0}};
const unsigned char kRdBu[11][3] = {{5, 48, 97}, {33, 102, 172}, {67, 147, 195}, {146, 197, 222}, {209, 229, 240}, {247, 247, 247}, {253, 219, 199}, {244, 165, 130}, {214, 96, 77}, {178, 24, 43}, {103, 0, 31}};
const unsigned char kSpectral[11][3] = {{94, 79, 162}, {50, 136, 189}, {102, 194, 165}, {171, 221, 164}, {230, 245, 152}, {255, 255, 191}, {254, 224, 139}, {253, 174, 97}, {244, 109, 67}, {213, 62, 79}, {158, 1, 66}};
const unsigned char kSeismic[5][3] = {{0, 0, 76}, {0, 0, 255}, {255, 255, 255}, {255, 0, 0}, {128, 0, 0}};
const unsigned char kTerrain[6][3] = {{51, 51, 153}, {0, 153, 255}, {0, 204, 102}, {255, 255, 153}, {128, 92, 84}, {255, 255, 255}};
const unsigned char kHot[5][3] = {{11, 0, 0}, {160, 0, 0}, {255, 60, 0}, {255, 220, 0}, {255, 255, 255}};

void lerpStops(const unsigned char (*c)[3], int n, float t, float rgb[3]) {
    const float x = t * float(n - 1);
    const int i = std::min(n - 2, std::max(0, int(std::floor(x))));
    const float f = x - float(i);
    for (int k = 0; k < 3; ++k) rgb[k] = (float(c[i][k]) + f * (float(c[i + 1][k]) - float(c[i][k]))) / 255.0f;
}
}  // namespace

void colormapRGB(Colormap m, float t, float rgb[3], bool reverse) {
    if (!(t == t)) t = 0;
    t = std::min(1.0f, std::max(0.0f, t));
    if (reverse) t = 1.0f - t;
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
        case Colormap::Plasma: lerpStops(kPlasma, 9, t, rgb); break;
        case Colormap::Inferno: lerpStops(kInferno, 9, t, rgb); break;
        case Colormap::Magma: lerpStops(kMagma, 9, t, rgb); break;
        case Colormap::Cividis: lerpStops(kCividis, 9, t, rgb); break;
        case Colormap::Rainbow: lerpStops(kRainbow, 9, t, rgb); break;
        case Colormap::RdBu: lerpStops(kRdBu, 11, t, rgb); break;
        case Colormap::Spectral: lerpStops(kSpectral, 11, t, rgb); break;
        case Colormap::Seismic: lerpStops(kSeismic, 5, t, rgb); break;
        case Colormap::Terrain: lerpStops(kTerrain, 6, t, rgb); break;
        case Colormap::Hot: lerpStops(kHot, 5, t, rgb); break;
        default: rgb[0] = rgb[1] = rgb[2] = t; break;
    }
    for (int k = 0; k < 3; ++k) rgb[k] = std::min(1.0f, std::max(0.0f, rgb[k]));
}

const char* colormapName(Colormap m) {
    static const char* n[] = {"Turbo", "Viridis", "Jet", "Cool-warm", "Grey", "Plasma", "Inferno", "Magma", "Cividis",
                              "Rainbow", "Blue-white-red", "Spectral", "Seismic", "Terrain", "Hot"};
    const int i = int(m);
    return i >= 0 && i < int(Colormap::Count) ? n[i] : "?";
}

bool colormapDiverging(Colormap m) {
    return m == Colormap::CoolWarm || m == Colormap::RdBu || m == Colormap::Spectral || m == Colormap::Seismic;
}
