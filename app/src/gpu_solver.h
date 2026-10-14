#pragma once
// Registers the Direct3D 11 compute accelerator with the engine (ressim::setBlockAccelerator). The device is
// created lazily on the first GPU linear solve; if no D3D11 compute device exists the engine falls back to CPR.
void registerGpuSolver();
