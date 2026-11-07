// Data embedded in ResSimStudio.exe as Win32 RCDATA resources (app/res/ressim.rc) on Windows and
// assembler .incbin blobs on Linux / macOS (resources.cpp), so the application has no
// external data-file dependencies.
#pragma once
#include <string_view>

enum AppResourceId : int {
    RES_KEYWORD_DOCS = 101,  // app/res/keyword_docs.txt - OPM/ECLIPSE keyword manual (3000+ keywords)
    RES_ICON_FONT = 102,     // app/res/icons.ttf - Lucide icon subset (ISC licence)
    RES_FONT_UI = 103,       // app/res/inter_regular.ttf  - Inter 400, latin (SIL OFL) - UI font where Segoe UI is absent
    RES_FONT_UI_BOLD = 104,  // app/res/inter_semibold.ttf - Inter 600, latin (SIL OFL)
    RES_FONT_MONO = 105,     // app/res/jbmono_regular.ttf - JetBrains Mono 400, latin (SIL OFL)
};

// Returns the bytes of an embedded resource (empty view if absent). The memory stays valid for the process lifetime.
std::string_view appResource(int id);
