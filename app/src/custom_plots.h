// C14: user-defined result plots. A plot is a named list of curves; each curve is a key into the run results:
//   F:<VEC>                    field vector (FOPR, FPR, FGMIP, FAVGT, DT ...)
//   W:<VEC>:<well>             well vector; <well> = name, "*" (all wells), "*P" (producers), "*I" (injectors)
//   C:<VEC>:<component>        field component vector (FCPR / FCIR mol/s, FCPT / FCIT mol, FCMPR kg/d, FCMPT kg)
//   WC:<VEC>:<component>:<well> well component vector (WCPR mol/s, WCPT / WCIT mol)
//   G:<VEC>[:<well>]           geothermal: WTEMP / WMW per well, COOLED (bulk-volume fraction cooled > 10 K)
//   B:<ARRAY>:<i>,<j>,<k>      cell value of a report-step array (1-based i, j, k)
//   O:<QTY>[:<well>]           observed (history-matching) data
// The plot library is saved in the user's profile (custom_plots.txt) and can be exported / imported (.rsplots);
// each plot can appear as its own tab of the Plots window after any run.
#pragma once
#include <string>
#include <vector>

struct CustomCurve {
    std::string key;
    std::string label;          // empty: automatic
    int axis = 0;               // 0 = left, 1 = right axis
    int style = 0;              // 0 line, 1 thin line, 2 markers, 3 line + markers
    double scale = 1;           // multiplier applied to the values
};

struct CustomPlot {
    std::string name = "My plot";
    std::vector<CustomCurve> curves;
    bool ownTab = true;         // show as its own tab in the Plots window
    bool logY = false;
    std::string ylabel, y2label;   // empty: unit of the first curve on that axis
};

bool writeCustomPlots(const std::string& path, const std::vector<CustomPlot>& plots);
bool readCustomPlots(const std::string& path, std::vector<CustomPlot>& plots);
std::string customPlotsToText(const std::vector<CustomPlot>& plots);
std::vector<CustomPlot> customPlotsFromText(const std::string& text);
