#!/usr/bin/env python3
"""C41 (G13): extract the TR("...") strings of ResSim Studio into a translation catalogue.

    python3 tools/i18n_extract.py                    # print the template (msgid<TAB>) to stdout
    python3 tools/i18n_extract.py --update lang/es.txt   # add missing ids to a catalogue, keep translations,
                                                         # report ids no longer used
Keys follow i18n::tr(): leading ICON_ macros and spaces, a '##' suffix and trailing spaces are not part of the key.
"""
import re, sys, glob, os

LIT = re.compile(r'TR[NLS]?\(\s*((?:ICON_[A-Z0-9_]+\s*)?(?:"(?:[^"\\]|\\.)*"\s*)+)\)')
STR = re.compile(r'"((?:[^"\\]|\\.)*)"')

# C45: helpers that translate their text arguments internally (widgets.h, nav.cpp tabItem, ui_table.cpp tableColumn,
# viewer_ui.cpp tip): every string-literal argument except printf formats and pure ImGui ids is a key
HELPER = re.compile(r'(?<![A-Za-z0-9_])(?:ui::)?(propD|propI|propF|propSliderD|propBool|propCombo|propText|propString|header|hint|help|label|tabItem|tableColumn|emptyState|tip|accentButton|dangerButton|secondaryButton)\(')
FMT = re.compile(r'^%[-+ #0-9.]*[a-zA-Z]+$')

def call_span(src, i):
    depth, j, q = 1, i, False
    while j < len(src) and depth:
        ch = src[j]
        if q:
            if ch == '\\': j += 1
            elif ch == '"': q = False
        elif ch == '"': q = True
        elif ch == '(': depth += 1
        elif ch == ')': depth -= 1
        j += 1
    return src[i:j - 1]

def split_args(body):
    args, depth, q, cur, k = [], 0, False, '', 0
    while k < len(body):
        ch = body[k]
        if q:
            cur += ch
            if ch == '\\': cur += body[k + 1]; k += 1
            elif ch == '"': q = False
        elif ch == '"': q = True; cur += ch
        elif ch in '([{': depth += 1; cur += ch
        elif ch in ')]}': depth -= 1; cur += ch
        elif ch == ',' and depth == 0: args.append(cur); cur = ''
        else: cur += ch
        k += 1
    args.append(cur)
    return args

def norm(s):
    s = bytes(s, 'utf-8').decode('unicode_escape').encode('latin-1').decode('utf-8')
    s = s.lstrip(' ')
    if '##' in s: s = s[:s.index('##')]
    return s.rstrip(' ')

def keys(root):
    out = []
    for f in sorted(glob.glob(os.path.join(root, 'app', 'src', '*.cpp')) + glob.glob(os.path.join(root, 'app', 'src', '*.h'))):
        src = open(f, encoding='utf-8').read()
        for m in HELPER.finditer(src):
            if src[max(0, m.start() - 5):m.start()].rstrip().endswith(('void', 'bool', 'int')): continue   # the definition
            for a in split_args(call_span(src, m.end())):
                a = a.strip()
                mm = re.fullmatch(r'(?:ICON_[A-Z0-9_]+\s*)?((?:"(?:[^"\\]|\\.)*"\s*)+)', a)
                if not mm: continue
                s = norm(''.join(STR.findall(mm.group(1))))
                if not s or FMT.match(s) or not re.search('[A-Za-z]{2}', s): continue
                if s not in out: out.append(s)
        for m in LIT.finditer(src):
            body = m.group(1)
            s = ''.join(STR.findall(body))
            s = bytes(s, 'utf-8').decode('unicode_escape').encode('latin-1').decode('utf-8')
            s = s.lstrip(' ')
            if '##' in s: s = s[:s.index('##')]
            s = s.rstrip(' ')
            if s and s not in out: out.append(s)
    return out

def esc(s): return s.replace('\\', '\\\\').replace('\t', '\\t').replace('\n', '\\n')

def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ks = keys(root)
    if len(sys.argv) >= 3 and sys.argv[1] == '--update':
        path = sys.argv[2]
        lines = open(path, encoding='utf-8').read().splitlines() if os.path.exists(path) else ['# LANGUAGE xx Name']
        have = set()
        for l in lines:
            if l and not l.startswith('#') and '\t' in l: have.add(l.split('\t', 1)[0])
        add = [k for k in ks if esc(k) not in have]
        unused = [h for h in have if h not in set(esc(k) for k in ks)]
        with open(path, 'w', encoding='utf-8', newline='\n') as f:
            f.write('\n'.join(lines).rstrip('\n') + '\n')
            for k in add: f.write(esc(k) + '\t\n')
        print(f'{path}: {len(ks)} ids in the source, {len(add)} added (untranslated), {len(unused)} unused')
        for u in unused: print('  unused:', u)
    else:
        for k in ks: print(esc(k) + '\t')

if __name__ == '__main__':
    main()
