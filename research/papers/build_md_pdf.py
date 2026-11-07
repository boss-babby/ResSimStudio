#!/usr/bin/env python3
"""Shared draft-PDF builder for all manuscripts: build_md_pdf.py SRC.md OUT.pdf "footer text".
Figures: markdown images pointing at .pdf (vector, journal files) are rasterised at 250 dpi for the layout; .png used as is.
Paths are relative to the directory of SRC.md."""
import os, re, sys, shutil, markdown, pymupdf
src, out = sys.argv[1], sys.argv[2]; foot = sys.argv[3] if len(sys.argv) > 3 else "draft"
base = os.path.dirname(os.path.abspath(src)); md = open(src, encoding="utf-8").read()
tmp = os.path.join("/tmp", "_mdpdf_%d" % os.getpid()); os.makedirs(tmp, exist_ok=True)
def fig(m):
    p = os.path.join(base, m.group(2)); name = os.path.splitext(os.path.basename(p))[0]
    if p.endswith(".pdf"):
        pymupdf.open(p)[0].get_pixmap(dpi=250).save(os.path.join(tmp, name + ".jpg"), jpg_quality=92); fn = name + ".jpg"
    else:
        fn = os.path.basename(p); shutil.copy(p, os.path.join(tmp, fn))
    return f'<p class="fig"><img src="{fn}" width="100%"/></p>'
md = re.sub(r"!\[([^\]]*)\]\(([^)]+\.(?:pdf|png))\)", fig, md)
# Python-Markdown needs a blank line before a list; authors (and GitHub) do not. Insert one where a list follows text.
_L = re.compile(r"^\s*(?:[*+-]|\d+\.)\s+")
_out = []
for _ln in md.split("\n"):
    if _L.match(_ln) and _out and _out[-1].strip() and not _L.match(_out[-1]) and not _out[-1].startswith(("  ", "|", ">")):
        _out.append("")
    _out.append(_ln)
md = "\n".join(_out).replace("𝒢", "<i>G</i>")
html = markdown.markdown(md, extensions=["tables"])
css = """* { font-family: sans-serif; } body { font-size: 9.5pt; line-height: 1.35; }
h1 { font-size: 15pt; margin-bottom: 4pt; } h2 { font-size: 12pt; margin-top: 10pt; color: #1a3d7c; } h3 { font-size: 10.5pt; margin-top: 8pt; }
p { margin: 3pt 0 3pt 0; text-align: justify; } blockquote { background-color: #fff4dc; padding: 4pt; margin: 4pt 0; font-size: 9pt; }
table { border-collapse: collapse; margin: 4pt 0; } th, td { border: 0.5px solid #999; padding: 2pt 4pt; font-size: 8pt; }
code { font-family: monospace; font-size: 8.5pt; } li { margin: 1pt 0; }"""
story = pymupdf.Story(html=html, user_css=css, archive=pymupdf.Archive(tmp))
writer = pymupdf.DocumentWriter(out); mb = pymupdf.paper_rect("a4"); where = mb + (50, 50, -50, -55); more = 1
while more:
    dev = writer.begin_page(mb); more, _ = story.place(where); story.draw(dev); writer.end_page()
writer.close(); doc = pymupdf.open(out)
for i, pg in enumerate(doc):
    pg.insert_text((50, mb.height - 30), f"{foot}    page {i + 1}/{len(doc)}", fontsize=7, color=(0.4, 0.4, 0.4))
doc.saveIncr(); print("pages", len(doc)); shutil.rmtree(tmp)
