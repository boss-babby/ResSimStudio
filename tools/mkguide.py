#!/usr/bin/env python3
"""docs/USER_GUIDE.md -> release/USER_GUIDE.html (self-contained, no external resources).
Run after every change of the guide:  python3 tools/mkguide.py"""
import os, re, sys
import markdown

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src = os.path.join(root, "docs", "USER_GUIDE.md")
md = open(src, encoding="utf-8").read()
body = markdown.markdown(md, extensions=["tables", "fenced_code", "toc", "sane_lists"])

# images next to the guide are embedded as data URIs (the HTML stays self-contained)
import base64
def _img(m):
    f = os.path.join(root, "docs", m.group(1))
    if not os.path.exists(f): return m.group(0)
    return 'src="data:image/png;base64,' + base64.b64encode(open(f, "rb").read()).decode() + '" style="max-width:100%"'
body = re.sub(r'src="([^":/]+\.png)"', _img, body)

# table of contents from the h2 headings
toc = []
for m in re.finditer(r'<h2 id="([^"]+)">(.*?)</h2>', body):
    toc.append(f'<a href="#{m.group(1)}">{re.sub("<[^>]+>", "", m.group(2))}</a>')

css = """
:root{--bg:#15171c;--panel:#1d2027;--text:#dfe3ea;--muted:#9aa3b2;--acc:#4084f6;--line:#2d313a}
@media (prefers-color-scheme: light){:root{--bg:#f7f8fa;--panel:#fff;--text:#1d2129;--muted:#5b6472;--acc:#1f66e0;--line:#dfe3ea}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--text);font:15px/1.6 "Segoe UI",system-ui,sans-serif}
nav{position:fixed;top:0;left:0;bottom:0;width:250px;overflow:auto;background:var(--panel);border-right:1px solid var(--line);padding:18px 14px}
nav b{display:block;margin-bottom:10px;font-size:16px}nav a{display:block;color:var(--muted);text-decoration:none;padding:4px 6px;border-radius:5px;font-size:14px}
nav a:hover{background:var(--line);color:var(--text)}
nav input{width:100%;padding:6px 8px;margin-bottom:12px;border-radius:6px;border:1px solid var(--line);background:var(--bg);color:var(--text)}
main{margin-left:250px;padding:28px 44px;max-width:1150px}
h1{font-size:28px;margin-top:0}h2{border-bottom:1px solid var(--line);padding-bottom:6px;margin-top:40px}h3{margin-top:28px;color:var(--acc)}
a{color:var(--acc)}table{border-collapse:collapse;width:100%;margin:12px 0;font-size:14px}
th,td{border:1px solid var(--line);padding:6px 9px;text-align:left;vertical-align:top}th{background:var(--panel)}
tr:nth-child(even) td{background:rgba(127,127,127,.05)}
code{background:var(--panel);padding:1px 5px;border-radius:4px;font-family:Consolas,monospace;font-size:13px}
pre{background:var(--panel);padding:12px 14px;border-radius:8px;overflow:auto;border:1px solid var(--line)}pre code{padding:0;background:none}
.hide{display:none}
@media (max-width:800px){nav{position:static;width:auto}main{margin-left:0;padding:16px}}
"""
js = """
function filt(q){q=q.toLowerCase();document.querySelectorAll('main tbody tr, main li, main p').forEach(function(e){
 e.classList.toggle('hide', q.length>0 && e.textContent.toLowerCase().indexOf(q)<0);});}
"""
html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>ResSim Studio - User Guide</title><style>{css}</style></head><body>
<nav><b>ResSim Studio guide</b><input placeholder="Filter (e.g. fracture, CO2)" oninput="filt(this.value)">{''.join(toc)}</nav>
<main>{body}</main><script>{js}</script></body></html>"""
out = os.path.join(root, "release", "USER_GUIDE.html")
os.makedirs(os.path.dirname(out), exist_ok=True)
open(out, "w", encoding="utf-8").write(html)
print("wrote", out, len(html), "bytes,", len(toc), "sections")
