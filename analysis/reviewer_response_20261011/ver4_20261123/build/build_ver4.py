#!/wanglab/jshen/anaconda3/bin/python
"""Build Reviewer_Letter_Response_ver4.docx from ver3 (ver3_input.docx, md5 08c42978...): replace the data of the
tables and figures that were not on the current pipeline, and the captions / sentences that restate them.
Edits are plain (untracked) edits of word/document.xml; existing tracked changes and comments of ver3 are untouched.
Every replacement value comes from a generated file (see SOURCES); nothing numeric is typed here except inside the
caption / sentence strings of edits_text.py, which the verification step checks against the same files.

usage: /wanglab/jshen/anaconda3/bin/python -I build_ver4.py [--out ver4.docx]
"""
import argparse
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from docx_helpers import Docx, cell_text, check, copy_row, delete_row, para_text, replace_image, replace_in_paragraph, set_cell_text, w  # noqa: E402

V4 = HERE.parent
LOG = []


def log(msg):
    LOG.append(msg)


def nrows(d, tbl):
    return len(d.rows(tbl))


def fit_rows(d, tbl, n_data, header_rows=1, footer_rows=0):
    """Make the table hold exactly n_data data rows between header and footer rows, by copying or deleting the last
    data row."""
    have = nrows(d, tbl) - header_rows - footer_rows
    while have < n_data:
        copy_row(d, tbl, header_rows + have - 1)
        have += 1
    while have > n_data:
        delete_row(tbl, header_rows + have - 1)
        have -= 1


def fill(d, title, matrix, col0=0, header_rows=1, footer_rows=0, resize=True):
    """Write matrix (list of rows; None = keep the cell) into the data rows of the table titled `title`,
    starting at column col0. Returns the number of cells whose text changed."""
    tbl = d.table_by_title(title)
    if resize:
        fit_rows(d, tbl, len(matrix), header_rows, footer_rows)
    rows = d.rows(tbl)[header_rows:len(d.rows(tbl)) - footer_rows]
    assert len(rows) == len(matrix), (title, len(rows), len(matrix))
    changed = 0
    for tr, vals in zip(rows, matrix):
        cells = tr.findall(w("tc"))
        assert col0 + len(vals) <= len(cells), (title, len(cells), len(vals))
        for tc, v in zip(cells[col0:], vals):
            if v is None:
                continue
            if cell_text(tc) != v:
                set_cell_text(d, tc, v)
                changed += 1
    log(f"{title}: {len(matrix)} data rows, {changed} cells changed")
    return changed


def header(d, title, col, text, row=0):
    tbl = d.table_by_title(title)
    tc = d.cell(tbl, row, col)
    if cell_text(tc) != text:
        set_cell_text(d, tc, text)
        log(f"{title}: header cell {col} -> {text!r}")


def para(d, regex):
    ps = d.paragraphs(regex)
    assert len(ps) == 1, (regex, len(ps))
    return ps[0]


def sub(d, regex, old, new, occurrence=None):
    """Replace old by new in the single body paragraph matching regex."""
    p = para(d, regex)
    replace_in_paragraph(d, p, old, new, occurrence=occurrence)
    log(f"text [{regex[:40]}]: {old[:60]!r} -> {new[:60]!r}")


def sub_all(d, regex, old, new, expect):
    """Replace old by new in every body paragraph matching regex (table copies have identical captions)."""
    ps = d.paragraphs(regex)
    assert len(ps) == expect, (regex, len(ps), expect)
    for p in ps:
        replace_in_paragraph(d, p, old, new)
    log(f"text x{expect} [{regex[:40]}]: {old[:60]!r} -> {new[:60]!r}")


def load(rel):
    return json.loads((V4 / rel).read_text())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE / "Reviewer_Letter_Response_ver4.docx"))
    ap.add_argument("--only", default="", help="comma list of step names (debug)")
    a = ap.parse_args()
    d = Docx(str(HERE / "ver3_input.docx"))
    import steps  # the edit steps live in steps.py so they can be reviewed apart from the machinery
    only = set(x for x in a.only.split(",") if x)
    for name, fn in steps.STEPS:
        if only and name not in only:
            continue
        before = len(LOG)
        fn(d, sys.modules[__name__])
        print(f"[{name}] {len(LOG) - before} operations")
    d.save(a.out)
    check(a.out)
    (Path(a.out).with_suffix(".buildlog.txt")).write_text("\n".join(LOG) + "\n")
    print("wrote", a.out, "operations:", len(LOG))


if __name__ == "__main__":
    main()
