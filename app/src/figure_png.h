// PNG export of a fig::Figure (Windows: Direct2D + WIC, figure_png_win.cpp).
#pragma once
#include "figure.h"

namespace fig {
bool writePNG(const Figure& f, const Style& s, const std::string& path, std::string* err = nullptr);
}
