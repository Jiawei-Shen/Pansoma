#!/wanglab/jshen/anaconda3/bin/python
"""Reviewer_Letter_Response_ver4.6_changes_vs_ver4.5.md: what ver4.6 changes in ver4.5, item by item (the user's list of
2026-11-26). Tables are read from the letters themselves (old = ver4.5, new = ver4.6) so the file shows what is in the
Word files; the text around them is written here.
usage: /wanglab/jshen/anaconda3/bin/python -I write_changes.py > Reviewer_Letter_Response_ver4.6_changes_vs_ver4.5.md"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "ver4_20261123/build"))
from docx_helpers import Docx, cell_text, w  # noqa: E402

LET = Path("/wanglab/jshen/scratch_jshen/Pansoma_papers/ver1/Pansoma_Resubmit")
OLD, NEW = Docx(str(LET / "Reviewer_Letter_Response_ver4.5.docx")), Docx(str(HERE / "Reviewer_Letter_Response_ver4.6.docx"))


def table(d, title, cols=None, rows=None):
    t = [[cell_text(c).replace("\n", " ") for c in r.findall(w("tc"))] for r in d.rows(d.table_by_title(title))]
    if cols:
        t = [[r[i] for i in cols] for r in t]
    if rows:
        t = [t[0]] + [t[i] for i in rows]
    return "\n".join(["| " + " | ".join(t[0]) + " |", "|" + "---|" * len(t[0])] + ["| " + " | ".join(r) + " |" for r in t[1:]])


def both(title, **kw):
    return f"**ver4.5：**\n\n{table(OLD, title, **kw)}\n\n**ver4.6：**\n\n{table(NEW, title, **kw)}\n"


print(Path(HERE / "changes_template.md").read_text()
      .replace("{{T2}}", both("Response Table 1.2"))
      .replace("{{T3}}", both("Response Table 1.5"))
      .replace("{{T7}}", both("Response Table 2.8"))
      .replace("{{T9}}", "**ver4.5（最后 4 行，已删）：**\n\n" + table(OLD, "Response Table 3.2", rows=[19, 20, 21, 22])
               + "\n\n**ver4.6（新增 9 行）：**\n\n" + table(NEW, "Response Table 3.2", rows=list(range(19, 28))) + "\n")
      .replace("{{T10}}", both("Response Table 3.9", cols=[0, 1, 5])), end="")
