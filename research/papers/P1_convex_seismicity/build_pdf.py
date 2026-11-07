"""manuscript.md + figs/*.pdf -> P1_manuscript_draft.pdf (A4) with PyMuPDF Story. Figures: the vector PDFs written by the
ResSim figure engine are rasterised at 300 dpi for the layout (the vector files stay in figs/ for the journal)."""
import os, re, markdown, pymupdf
import sys
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "manuscript.md")   # e.g. jgr/P1_JGR_long_v0.5.md
md = open(SRC, encoding="utf-8").read()
tmp = os.path.join(HERE, "figs", "_png"); os.makedirs(tmp, exist_ok=True)
def fig(m):
    pdf = m.group(2); name = os.path.splitext(os.path.basename(pdf))[0]
    png = os.path.join(tmp, name + ".jpg")
    pymupdf.open(os.path.join(HERE, pdf))[0].get_pixmap(dpi=250).save(png, jpg_quality=92)
    return f'<p class="fig"><img src="{name}.jpg" width="100%"/></p>'
md = re.sub(r"!\[([^\]]*)\]\(([^)]+\.pdf)\)", fig, md)
html = markdown.markdown(md, extensions=["tables"])
css = """
* { font-family: sans-serif; }
body { font-size: 9.5pt; line-height: 1.35; }
h1 { font-size: 15pt; margin-bottom: 4pt; }
h2 { font-size: 12pt; margin-top: 10pt; color: #1a3d7c; }
h3 { font-size: 10.5pt; margin-top: 8pt; }
p { margin: 3pt 0 3pt 0; text-align: justify; }
blockquote { background-color: #fff4dc; padding: 4pt; margin: 4pt 0; font-size: 9pt; }
table { border-collapse: collapse; margin: 4pt 0; }
th, td { border: 0.5px solid #999; padding: 2pt 4pt; font-size: 8.5pt; }
code { font-family: monospace; font-size: 8.5pt; }
li { margin: 1pt 0; }
"""
story = pymupdf.Story(html=html, user_css=css, archive=pymupdf.Archive(tmp))
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "P1_manuscript_draft.pdf")
writer = pymupdf.DocumentWriter(out)
mb = pymupdf.paper_rect("a4"); where = mb + (50, 50, -50, -55)
more = 1
while more:
    dev = writer.begin_page(mb)
    more, _ = story.place(where)
    story.draw(dev)
    writer.end_page()
writer.close()
doc = pymupdf.open(out)
for i, pg in enumerate(doc):   # page numbers + draft footer
    pg.insert_text((50, mb.height - 30), f"P1 draft v0.5 (JGR long version; pre-submission) - figures rendered by ResSim Studio's figure engine    page {i + 1}/{len(doc)}", fontsize=7, color=(0.4, 0.4, 0.4))
doc.saveIncr()
print("pages", len(doc))
for f in os.listdir(tmp): os.remove(os.path.join(tmp, f))
os.rmdir(tmp)
