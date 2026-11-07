#include "resources.h"

#ifdef _WIN32
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

#else
// C39 Linux / macOS: the same files are linked into the executable with the assembler's .incbin (RESSIM_RES_DIR = absolute
// path of app/res, set by CMake), so the portable build has no external data-file dependency either.
#ifndef RESSIM_RES_DIR
#error "RESSIM_RES_DIR must point to app/res"
#endif
#define RS_STR2(x) #x
#define RS_STR(x) RS_STR2(x)
#if defined(__APPLE__)
#define RS_SECTION ".const_data\n"
#define RS_SYM(n) "_" #n
#else
#define RS_SECTION ".section .rodata\n"
#define RS_SYM(n) #n
#endif
#define RS_INCBIN(name, file)                                                                  \
    __asm__(RS_SECTION ".global " RS_SYM(name##_begin) "\n.balign 16\n" RS_SYM(name##_begin) ":\n" \
            ".incbin \"" RS_STR(RESSIM_RES_DIR) "/" file "\"\n"                                   \
            ".global " RS_SYM(name##_end) "\n" RS_SYM(name##_end) ":\n.byte 0\n.balign 16\n.text\n");  \
    extern "C" const char name##_begin[];                                                      \
    extern "C" const char name##_end[];

RS_INCBIN(ressim_res_kwdocs, "keyword_docs.txt")
RS_INCBIN(ressim_res_icons, "icons.ttf")
RS_INCBIN(ressim_res_ui, "inter_regular.ttf")
RS_INCBIN(ressim_res_uib, "inter_semibold.ttf")
RS_INCBIN(ressim_res_mono, "jbmono_regular.ttf")

std::string_view appResource(int id) {
    switch (id) {
        case RES_KEYWORD_DOCS: return std::string_view(ressim_res_kwdocs_begin, size_t(ressim_res_kwdocs_end - ressim_res_kwdocs_begin));
        case RES_ICON_FONT: return std::string_view(ressim_res_icons_begin, size_t(ressim_res_icons_end - ressim_res_icons_begin));
        case RES_FONT_UI: return std::string_view(ressim_res_ui_begin, size_t(ressim_res_ui_end - ressim_res_ui_begin));
        case RES_FONT_UI_BOLD: return std::string_view(ressim_res_uib_begin, size_t(ressim_res_uib_end - ressim_res_uib_begin));
        case RES_FONT_MONO: return std::string_view(ressim_res_mono_begin, size_t(ressim_res_mono_end - ressim_res_mono_begin));
        default: return {};
    }
}
#endif
