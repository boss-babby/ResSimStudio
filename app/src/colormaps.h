// Colour maps of the 3-D viewers, the panels and the plots (pure C++, no graphics API - also used by the tests).
#pragma once

enum class Colormap { Turbo, Viridis, Jet, CoolWarm, Grey,
                      // C30
                      Plasma, Inferno, Magma, Cividis, Rainbow, RdBu, Spectral, Seismic, Terrain, Hot, Count };

// t in [0, 1] -> rgb in [0, 1]; reverse: 1 - t
void colormapRGB(Colormap m, float t, float rgb[3], bool reverse = false);
const char* colormapName(Colormap m);
// sequential / perceptual (good for magnitudes) vs diverging (good for signed changes around 0)
bool colormapDiverging(Colormap m);
