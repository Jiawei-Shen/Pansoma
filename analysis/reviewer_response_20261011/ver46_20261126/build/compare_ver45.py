#!/wanglab/jshen/anaconda3/bin/python
"""Every difference between ver4.5 and ver4.6: body paragraphs (visible text), table cells, embedded images.
usage: /wanglab/jshen/anaconda3/bin/python -I compare_ver45.py OLD.docx NEW.docx"""
import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "ver4_20261123/build"))
from docx_helpers import Docx, cell_text, para_text, w  # noqa: E402

a, b = Docx(sys.argv[1]), Docx(sys.argv[2])
pa, pb = [para_text(p) for p in a.paragraphs()], [para_text(p) for p in b.paragraphs()]
print(f"body paragraphs: {len(pa)} -> {len(pb)}")
assert len(pa) == len(pb)
n = 0
for i, (x, y) in enumerate(zip(pa, pb)):
    if x != y:
        n += 1
        # common prefix / suffix trimmed
        k = 0
        while k < min(len(x), len(y)) and x[k] == y[k]:
            k += 1
        j = 0
        while j < min(len(x), len(y)) - k and x[-1 - j] == y[-1 - j]:
            j += 1
        print(f"\n[paragraph {i}] {x[:60]!r}\n   - {x[max(0, k - 30):len(x) - j]!r}\n   + {y[max(0, k - 30):len(y) - j]!r}")
print(f"\nchanged body paragraphs: {n}")
ta, tb = a.tables(), b.tables()
assert len(ta) == len(tb)
kids = list(b.body)
for i, (x, y) in enumerate(zip(ta, tb)):
    ra = [[cell_text(c) for c in r.findall(w("tc"))] for r in a.rows(x)]
    rb = [[cell_text(c) for c in r.findall(w("tc"))] for r in b.rows(y)]
    if ra != rb:
        title = ""
        k = kids.index(y)
        for prev in kids[max(0, k - 3):k][::-1]:
            if prev.tag == w("p") and para_text(prev).strip():
                title = para_text(prev).strip()[:60]
                break
        same_shape = len(ra) == len(rb)
        cells = sum(u != v for p, q in zip(ra, rb) for u, v in zip(p, q)) if same_shape else "-"
        print(f"table {i} [{title}]: rows {len(ra)} -> {len(rb)}, changed cells {cells}")
for name in sorted(set(a.parts) | set(b.parts)):
    if name.startswith("word/media/"):
        ha = hashlib.md5(a.parts[name]).hexdigest() if name in a.parts else None
        hb = hashlib.md5(b.parts[name]).hexdigest() if name in b.parts else None
        if ha != hb:
            print("image", name, ha, "->", hb)
for name in ("word/comments.xml",):
    print(name, "identical" if a.parts.get(name) == b.parts.get(name) else "DIFFERENT")
