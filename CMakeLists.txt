cmake_minimum_required(VERSION 3.20)
project(ResSimStudio LANGUAGES C CXX)

# ---------------------------------------------------------------------------
# ResSim Studio - black-oil reservoir simulator
#   ressim_engine : portable C++17 simulation engine (no UI dependencies)
#   ResSimStudio  : desktop app - Windows: Direct3D 11 + Direct2D/DirectWrite + Dear ImGui
#                                 Linux / macOS (C39): GLFW + OpenGL 3.3 core + Dear ImGui + plutovg + tinyfiledialogs
#   ResSimCLI     : console front end (single runs, batches, regression)
#   ressim_tests  : engine regression tests (console)
# ---------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(RESSIM_AVX2 "Compile with AVX2 instructions" ON)
option(RESSIM_BUILD_TESTS "Build the engine regression tests" ON)

find_package(OpenMP)

# ---- engine ----------------------------------------------------------------
file(GLOB ENGINE_SRC CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/engine/src/*.cpp)
add_library(ressim_engine STATIC ${ENGINE_SRC})
target_include_directories(ressim_engine PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/engine/include)
if(OpenMP_CXX_FOUND)
  target_link_libraries(ressim_engine PUBLIC OpenMP::OpenMP_CXX)
endif()

if(MSVC)
  # /fp:precise keeps the solver deterministic; /MP parallel compile; /utf-8 for the UI strings
  set(RESSIM_FLAGS /W3 /MP /utf-8 /permissive- /Zc:__cplusplus $<$<CONFIG:Release,RelWithDebInfo>:/O2 /Ob2 /GL>)
  if(RESSIM_AVX2)
    list(APPEND RESSIM_FLAGS /arch:AVX2)
  endif()
  set(RESSIM_LINK $<$<CONFIG:Release,RelWithDebInfo>:/LTCG>)
else()
  set(RESSIM_FLAGS -Wall -Wextra $<$<CONFIG:Release>:-O3>)
  if(RESSIM_AVX2)
    list(APPEND RESSIM_FLAGS -mavx2 -mfma)
  endif()
  set(RESSIM_LINK "")
endif()
target_compile_options(ressim_engine PRIVATE ${RESSIM_FLAGS})
target_compile_definitions(ressim_engine PRIVATE _USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS)

# ---- tests -----------------------------------------------------------------
if(RESSIM_BUILD_TESTS)
  add_executable(ressim_tests tests/run_tests.cpp)
  target_link_libraries(ressim_tests PRIVATE ressim_engine)
  target_compile_options(ressim_tests PRIVATE ${RESSIM_FLAGS})
  target_link_options(ressim_tests PRIVATE ${RESSIM_LINK})
endif()

# ---- developer tools: WPWE event timelines vs an OPM reference (C44) --------------
option(RESSIM_BUILD_DEVTOOLS "Build developer comparison tools (wpwe_compare, vcmp)" OFF)
if(RESSIM_BUILD_DEVTOOLS)
  foreach(t wpwe_compare dev/vcmp)
    get_filename_component(tn ${t} NAME)
    add_executable(${tn} tools/${t}.cpp)
    target_link_libraries(${tn} PRIVATE ressim_engine)
    target_compile_options(${tn} PRIVATE ${RESSIM_FLAGS})
  endforeach()
endif()

# ---- console front end ---------------------------------------------------------
add_executable(ResSimCLI tools/ressim_cli.cpp)
target_link_libraries(ResSimCLI PRIVATE ressim_engine)
target_compile_options(ResSimCLI PRIVATE ${RESSIM_FLAGS})
target_compile_definitions(ResSimCLI PRIVATE _USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS)
if(NOT WIN32)
  target_link_libraries(ResSimCLI PRIVATE ${CMAKE_DL_LIBS})
endif()

# ---- desktop application ------------------------------------------------------
set(IMGUI_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/imgui)
set(IMPLOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/implot)
if(WIN32)
  add_library(imgui_lib STATIC
    ${IMGUI_DIR}/imgui.cpp ${IMGUI_DIR}/imgui_draw.cpp ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp ${IMGUI_DIR}/imgui_demo.cpp
    ${IMGUI_DIR}/backends/imgui_impl_win32.cpp ${IMGUI_DIR}/backends/imgui_impl_dx11.cpp
    ${IMPLOT_DIR}/implot.cpp ${IMPLOT_DIR}/implot_items.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/textedit/TextEditor.cpp)
  target_include_directories(imgui_lib PUBLIC ${IMGUI_DIR} ${IMGUI_DIR}/backends ${IMPLOT_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/third_party/textedit)
  if(MSVC)
    target_compile_options(imgui_lib PRIVATE /MP /utf-8 /bigobj $<$<CONFIG:Release,RelWithDebInfo>:/O2>)
  endif()

  file(GLOB APP_SRC CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/app/src/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/app/src/*.h)
  list(FILTER APP_SRC EXCLUDE REGEX "/(main_glfw|gl_view|viewport3d_gl|geoscene3d_gl|figure_png_vg|gpu_compute_stub)\\.cpp$")
  add_executable(ResSimStudio WIN32 ${APP_SRC} ${CMAKE_CURRENT_SOURCE_DIR}/app/res/ressim.rc)
  target_include_directories(ResSimStudio PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/app/src)
  target_link_libraries(ResSimStudio PRIVATE ressim_engine imgui_lib
    d3d11 dxgi d3dcompiler d2d1 dwrite dwmapi shell32 comdlg32 ole32 user32 gdi32 imm32 windowscodecs uuid)
  target_compile_definitions(ResSimStudio PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES)
  target_compile_options(ResSimStudio PRIVATE ${RESSIM_FLAGS})
  target_link_options(ResSimStudio PRIVATE ${RESSIM_LINK})
  set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ResSimStudio)
  set_target_properties(ResSimStudio PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:ResSimStudio>")
else()
  # C39 Linux / macOS: same App, portable renderer (OpenGL 3.3 core, no loader: system libGL / OpenGL.framework)
  find_package(OpenGL REQUIRED)
  find_package(glfw3 3.3 QUIET)
  if(NOT glfw3_FOUND)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GLFW REQUIRED IMPORTED_TARGET glfw3)
    add_library(glfw ALIAS PkgConfig::GLFW)
  endif()
  set(PLUTOVG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/plutovg)
  file(GLOB PLUTOVG_SRC CONFIGURE_DEPENDS ${PLUTOVG_DIR}/source/*.c)
  add_library(plutovg STATIC ${PLUTOVG_SRC})
  target_include_directories(plutovg PUBLIC ${PLUTOVG_DIR}/include)
  target_compile_definitions(plutovg PUBLIC PLUTOVG_BUILD_STATIC)
  set_target_properties(plutovg PROPERTIES C_STANDARD 11 POSITION_INDEPENDENT_CODE ON)
  if(NOT WIN32)
    target_link_libraries(plutovg PRIVATE m)
  endif()
  add_library(tinyfd STATIC ${CMAKE_CURRENT_SOURCE_DIR}/third_party/tinyfd/tinyfiledialogs.c)
  target_include_directories(tinyfd PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/third_party/tinyfd)

  add_library(imgui_lib STATIC
    ${IMGUI_DIR}/imgui.cpp ${IMGUI_DIR}/imgui_draw.cpp ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp ${IMGUI_DIR}/imgui_demo.cpp
    ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp
    ${IMPLOT_DIR}/implot.cpp ${IMPLOT_DIR}/implot_items.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/textedit/TextEditor.cpp)
  target_include_directories(imgui_lib PUBLIC ${IMGUI_DIR} ${IMGUI_DIR}/backends ${IMPLOT_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/third_party/textedit)
  target_link_libraries(imgui_lib PUBLIC glfw OpenGL::GL)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # implot_items.cpp is a template-heavy unit: keep GCC's peak memory low (machines with 2 GB)
    set_source_files_properties(${IMPLOT_DIR}/implot_items.cpp PROPERTIES COMPILE_OPTIONS "-O1;--param;ggc-min-expand=20;--param;ggc-min-heapsize=65536")
  endif()

  file(GLOB APP_SRC CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/app/src/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/app/src/*.h)
  list(FILTER APP_SRC EXCLUDE REGEX "/(main|figure_png_win|gpu_compute|overlay_d2d)\\.cpp$")
  add_executable(ResSimStudio ${APP_SRC})
  target_include_directories(ResSimStudio PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/app/src)
  target_link_libraries(ResSimStudio PRIVATE ressim_engine imgui_lib plutovg tinyfd glfw OpenGL::GL ${CMAKE_DL_LIBS})
  target_compile_definitions(ResSimStudio PRIVATE _USE_MATH_DEFINES RESSIM_GL RESSIM_RES_DIR=${CMAKE_CURRENT_SOURCE_DIR}/app/res)
  target_compile_options(ResSimStudio PRIVATE ${RESSIM_FLAGS} -Wno-unknown-pragmas)
  # resources.cpp embeds app/res with .incbin: rebuild it when a resource changes
  file(GLOB RES_FILES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/app/res/*.ttf ${CMAKE_CURRENT_SOURCE_DIR}/app/res/*.txt)
  set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/app/src/resources.cpp PROPERTIES OBJECT_DEPENDS "${RES_FILES}")
  if(APPLE)
    target_link_libraries(ResSimStudio PRIVATE "-framework Cocoa" "-framework IOKit")
  endif()
  # C32 ResGeo: the same executable started under the name ResGeo opens the standalone geomodelling front end
  add_custom_command(TARGET ResSimStudio POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink ResSimStudio $<TARGET_FILE_DIR:ResSimStudio>/ResGeo)
  include(GNUInstallDirs)
  install(TARGETS ResSimStudio ResSimCLI RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/packaging/ressim-studio.desktop DESTINATION ${CMAKE_INSTALL_DATADIR}/applications OPTIONAL)
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/app/res/ressim.svg DESTINATION ${CMAKE_INSTALL_DATADIR}/icons/hicolor/scalable/apps RENAME ressim-studio.svg)
endif()
