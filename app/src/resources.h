// Data embedded in ResSimStudio.exe as Win32 RCDATA resources (app/res/ressim.rc), so the application has no
// external data-file dependencies.
#pragma once
#include <string_view>

enum AppResourceId : int {
    RES_KEYWORD_DOCS = 101,  // app/res/keyword_docs.txt - OPM/ECLIPSE keyword manual (3000+ keywords)
    RES_ICON_FONT = 102,     // app/res/icons.ttf - Lucide icon subset (ISC licence)
};

// Returns the bytes of an embedded resource (empty view if absent). The memory stays valid for the process lifetime.
std::string_view appResource(int id);
