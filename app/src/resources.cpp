#include "resources.h"
#include <windows.h>

std::string_view appResource(int id) {
    HMODULE m = GetModuleHandleW(nullptr);
    HRSRC r = FindResourceW(m, MAKEINTRESOURCEW(id), MAKEINTRESOURCEW(10) /* RT_RCDATA */);
    if (!r) return {};
    HGLOBAL g = LoadResource(m, r);
    if (!g) return {};
    const void* p = LockResource(g);
    const DWORD n = SizeofResource(m, r);
    return p && n ? std::string_view(static_cast<const char*>(p), size_t(n)) : std::string_view();
}
