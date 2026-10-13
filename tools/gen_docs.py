#!/usr/bin/env python3
"""Generate app/res/keyword_docs.txt (keyword documentation loaded at run time by the deck editor,
keyword reference, hover tooltips, autocomplete and keyword forms) from the OPM Flow reference manual
index of opm-flow-editor-support (vscode-extension/data/keyword_index_compact.json).

usage: python3 tools/gen_docs.py EDITOR_INDEX_JSON [OUT]
Format (TAB separated, one record per line; trailing columns were added in v2, readers ignore extras):
  K name sections alias_of summary size_kind expected_columns variadic
  P index name manual_name default units(field|metric) options description record value_type dimension
  E example (\\n separated)
The documentation comes from the Open Porous Media project (https://opm-project.org), GPL-3.0.
"""
import json, os, re, sys

def clean(t, lim=None):
    t = re.sub(r"\s+", " ", (t or "")).strip()
    if lim and len(t) > lim:
        cut = t.rfind(". ", 0, lim)
        t = t[: cut + 1] if cut > lim * 0.5 else t[:lim].rstrip() + "..."
    return t.replace("\t", " ")

def main(idxfile, out):
    idx = json.load(open(idxfile, encoding="utf-8"))
    lines = ["# ResSim keyword documentation v2. Generated from the OPM Flow reference manual index\n",
             "# (opm-flow-editor-support / opm-reference-manual, https://opm-project.org, GPL-3.0).\n",
             "# K<TAB>name<TAB>sections<TAB>alias_of<TAB>summary<TAB>size_kind<TAB>expected_columns(per record, comma separated)<TAB>variadic\n",
             "# P<TAB>index<TAB>name<TAB>manual_name<TAB>default<TAB>units(field|metric)<TAB>options<TAB>description<TAB>record<TAB>value_type<TAB>dimension\n",
             "# E<TAB>example line (\\n separated)\n"]
    for name in sorted(idx):
        v = idx[name]
        cols = str(v.get("expected_columns") or "")
        if v.get("records_meta"):
            cols = ",".join(str(m.get("expected_columns") or "") for m in v["records_meta"])
        lines.append("K\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (name, ",".join(v.get("sections") or []), v.get("alias_of") or "",
                     clean(v.get("summary"), 1200), v.get("size_kind") or "", cols, "1" if v.get("variadic_record") else ""))
        if v.get("alias_of"):
            continue
        for p in v.get("parameters") or []:
            u = p.get("units") or {}
            units = "%s|%s" % (u.get("field", ""), u.get("metric", ""))
            lines.append("P\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                p.get("index", ""), p.get("name", ""), p.get("manual_name", "") or "",
                clean(str(p.get("default", "") or ""), 60), clean(units, 80), ",".join(p.get("options") or [])[:200],
                clean(p.get("description"), 420), p.get("record", "") or "", p.get("value_type", "") or "",
                clean(str(p.get("dimension", "") or ""), 60)))
        ex = (v.get("example") or "").strip()
        if ex:
            ex = "\n".join(l for l in ex.splitlines()[:14])
            lines.append("E\t%s\n" % ex.replace("\t", " ").replace("\n", "\\n"))
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    open(out, "w", encoding="utf-8").write("".join(lines))
    print("docs:", len(idx), "entries,", sum(len(l) for l in lines) // 1024, "KB")

if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(here), "app/res/keyword_docs.txt"))
