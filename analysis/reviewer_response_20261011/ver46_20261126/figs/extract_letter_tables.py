#!/wanglab/jshen/anaconda3/bin/python
"""letter_tables.json: the cell texts of Response Tables 1.2 and 3.2 of the ver4.6 letter
(../build/Reviewer_Letter_Response_ver4.6.docx), read with ver4_20261123/build/docx_helpers.py. The figure scripts
run in the sambcfenv environment (matplotlib), which has no lxml, so the tables are extracted here first.
usage: /wanglab/jshen/anaconda3/bin/python -I extract_letter_tables.py"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "ver4_20261123/build"))
from docx_helpers import Docx, cell_text, w  # noqa: E402

src = HERE.parent / "build/Reviewer_Letter_Response_ver4.6.docx"
d = Docx(str(src))
out = {"source": str(src)}
for t in ("Response Table 1.2", "Response Table 3.2"):
    out[t] = [[cell_text(c) for c in tr.findall(w("tc"))] for tr in d.rows(d.table_by_title(t))]
assert len(out["Response Table 1.2"]) == 10 and len(out["Response Table 3.2"]) == 28
(HERE / "letter_tables.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print("wrote letter_tables.json")
