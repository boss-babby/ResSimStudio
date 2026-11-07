#!/usr/bin/env python3
"""Compile-check the HLSL shaders embedded in the app sources with DirectXShaderCompiler (Linux dxc, SM 6.0 targets,
HLSL 2018 semantics close to fxc / SM 5). The app compiles them at run time with D3DCompile (SM 5.0), which cannot run
on the Linux build machine.
  dxc: https://github.com/microsoft/DirectXShaderCompiler/releases (linux_dxc_*.tar.gz) -> DXC=/path/to/bin/dxc
Usage: DXC=/tmp/sd/dxc/bin/dxc python3 tools/check_hlsl.py"""
import os, subprocess, sys
DXC = os.environ.get('DXC', '/tmp/sd/dxc/bin/dxc')
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
JOBS = [('app/src/geoscene3d.cpp', 'kGeoShader', ['VSMesh:vs', 'PSMesh:ps', 'VSTex:vs', 'PSTex:ps', 'VSColor:vs', 'PSColor:ps', 'VSQuad:vs', 'PSVol:ps'], ['MSDEPTH=0', 'MSDEPTH=1']),
        ('app/src/viewport3d.cpp', 'kShader', ['VSMain:vs', 'PSMain:ps', 'VSColor:vs', 'PSColor:ps'], [''])]
def extract(path, name):
    s = open(os.path.join(ROOT, path)).read(); i = s.index(name); j = s.index('R"(', i) + 3; return s[j:s.index(')"', j)]
bad = 0
tmp = os.path.join(os.environ.get('TMPDIR', '/tmp'), 'ressim_hlsl'); os.makedirs(tmp, exist_ok=True)
for path, name, eps, defs in JOBS:
    f = os.path.join(tmp, name + '.hlsl'); open(f, 'w').write(extract(path, name))
    for d in defs:
        for ep in eps:
            e, t = ep.split(':')
            cmd = [DXC, '-HV', '2018', '-T', t + '_6_0', '-E', e, f, '-Fo', os.devnull] + (['-D', d] if d else [])
            r = subprocess.run(cmd, capture_output=True, text=True); out = (r.stdout + r.stderr).strip()
            if r.returncode != 0 or 'warning' in out or 'error' in out:
                bad += 1; print(f'{path} {e} {d}:\n{out[:2000]}')
print(f'HLSL check: {bad} problem(s)'); sys.exit(1 if bad else 0)
