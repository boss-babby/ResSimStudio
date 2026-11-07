// OpenGL 3.3 core declarations for the portable renderer (Linux: system libGL prototypes, macOS: OpenGL.framework).
#pragma once
#if defined(__APPLE__)
#define GL_SILENCE_DEPRECATION
#include <OpenGL/gl3.h>
#else
#ifndef GL_GLEXT_PROTOTYPES
#define GL_GLEXT_PROTOTYPES 1
#endif
#include <GL/gl.h>
#include <GL/glext.h>
#endif
