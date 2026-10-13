cmake_minimum_required(VERSION 3.20)
project(ResSimStudio LANGUAGES CXX)

# ---------------------------------------------------------------------------
# ResSim Studio - black-oil reservoir simulator
#   ressim_engine : portable C++17 simulation engine (no UI dependencies)
#   ResSimStudio  : Windows desktop app (Direct3D 11 + Direct2D/DirectWrite + Dear ImGui)
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

# ---- desktop application (Windows only) -------------------------------------
if(WIN32)
  set(IMGUI_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/imgui)
  set(IMPLOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/implot)
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
  add_executable(ResSimStudio WIN32 ${APP_SRC} ${CMAKE_CURRENT_SOURCE_DIR}/app/res/ressim.rc)
  target_include_directories(ResSimStudio PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/app/src)
  target_link_libraries(ResSimStudio PRIVATE ressim_engine imgui_lib
    d3d11 dxgi d3dcompiler d2d1 dwrite dwmapi shell32 comdlg32 ole32 user32 gdi32 imm32 windowscodecs uuid)
  target_compile_definitions(ResSimStudio PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES)
  target_compile_options(ResSimStudio PRIVATE ${RESSIM_FLAGS})
  target_link_options(ResSimStudio PRIVATE ${RESSIM_LINK})
  set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ResSimStudio)
  set_target_properties(ResSimStudio PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:ResSimStudio>")
endif()
