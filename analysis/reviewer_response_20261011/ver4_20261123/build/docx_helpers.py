#!/wanglab/jshen/anaconda3/bin/python
"""Helpers to edit Reviewer_Letter_Response_ver*.docx with lxml + zipfile only
(python-docx is not installed anywhere; lxml exists only in /wanglab/jshen/anaconda3/bin/python).

Every part except word/document.xml (and, for images, word/media/* and
word/_rels/document.xml.rels) is copied byte for byte, in the original zip order.

    from docx_helpers import Docx, set_cell_text, copy_row, replace_image, replace_in_paragraph
    d = Docx(src)                                   # never the user's file: work on a copy
    t = d.table_by_title("Response Table 1.2")      # or d.tables()[1]
    set_cell_text(d.cell(t, 1, 3), "0.100 / 0.200 / 0.133")
    d.save(dst); check(dst)

"Visible" text = the text with all tracked changes accepted: runs inside w:del / w:moveFrom
are ignored, runs inside w:ins / w:moveTo are read and edited in place.

tracked=None  : plain edit (text that sits inside an existing w:ins stays inside it).
tracked=(author, "2026-11-23T16:00:00Z") : the old text becomes a tracked deletion and the new
text a tracked insertion of that author, so Word's Review pane lists the edit.
"""
import copy
import io
import random
import re
import struct
import zipfile

from lxml import etree

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
W14 = "http://schemas.microsoft.com/office/word/2010/wordml"
R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
WP = "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"
A = "http://schemas.openxmlformats.org/drawingml/2006/main"
PIC = "http://schemas.openxmlformats.org/drawingml/2006/picture"
PKG_REL = "http://schemas.openxmlformats.org/package/2006/relationships"
XML_SPACE = "{http://www.w3.org/XML/1998/namespace}space"
IMAGE_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships/image"
NS = {"w": W, "w14": W14, "r": R, "wp": WP, "a": A, "pic": PIC}

EMU_PER_INCH = 914400
TEXT_WIDTH_EMU = 5943600      # 6.5 in: Letter page, 1 in margins (sectPr of ver3)
TEXT_HEIGHT_EMU = 8229600     # 9.0 in


def w(tag):
    return "{%s}%s" % (W, tag)


_HIDDEN = (w("del"), w("moveFrom"))
_ANNOT = ("ins", "del", "moveFrom", "moveTo", "moveFromRangeStart", "moveFromRangeEnd",
          "moveToRangeStart", "moveToRangeEnd", "bookmarkStart", "bookmarkEnd",
          "commentRangeStart", "commentRangeEnd", "commentReference")


class Docx:
    def __init__(self, path):
        self.path = path
        with zipfile.ZipFile(path) as z:
            self.infos = z.infolist()
            self.parts = {i.filename: z.read(i.filename) for i in self.infos}
        self.order = [i.filename for i in self.infos]
        self.doc = etree.fromstring(self.parts["word/document.xml"])
        self.rels = etree.fromstring(self.parts["word/_rels/document.xml.rels"])
        self.body = self.doc.find(w("body"))
        self._rels_dirty = False

    # ---- lookup -------------------------------------------------------------
    def tables(self):
        """Top-level tables in document order (ver3 has no nested tables)."""
        return self.body.findall(w("tbl"))

    def table_by_title(self, title):
        """The table that follows the title paragraph, e.g. 'Response Table 1.2'.
        The match is on the paragraph start and refuses a longer number (1.1 vs 1.10)."""
        pat = re.compile(r"^\s*" + re.escape(title) + r"(?![0-9])")
        kids = list(self.body)
        hits = []
        for i, k in enumerate(kids):
            if k.tag == w("p") and pat.match(para_text(k)):
                for nxt in kids[i + 1:i + 4]:
                    if nxt.tag == w("tbl"):
                        hits.append(nxt)
                        break
                    if nxt.tag == w("p") and para_text(nxt).strip():
                        break
        if len(hits) != 1:
            raise KeyError("%r: %d tables found" % (title, len(hits)))
        return hits[0]

    @staticmethod
    def rows(tbl):
        return tbl.findall(w("tr"))

    @staticmethod
    def cell(tbl, r, c):
        return tbl.findall(w("tr"))[r].findall(w("tc"))[c]

    def paragraphs(self, regex=None):
        """Body-level paragraphs (not table cells), optionally those whose visible text matches."""
        ps = self.body.findall(w("p"))
        if regex is None:
            return ps
        rx = re.compile(regex)
        return [p for p in ps if rx.search(para_text(p))]

    def rel_target(self, rid):
        for r in self.rels:
            if r.get("Id") == rid:
                return r.get("Target")
        raise KeyError(rid)

    def drawings(self):
        """[(w:drawing element, rId, 'media/imageN.png')] in document order."""
        out = []
        for d in self.doc.iter(w("drawing")):
            blip = d.find(".//{%s}blip" % A)
            rid = blip.get("{%s}embed" % R)
            out.append((d, rid, self.rel_target(rid)))
        return out

    # ---- ids ----------------------------------------------------------------
    def next_annotation_id(self):
        m = 0
        for tag in _ANNOT:
            for e in self.doc.iter(w(tag)):
                v = e.get(w("id"))
                if v is not None and v.lstrip("-").isdigit():
                    m = max(m, int(v))
        self._next_id = max(getattr(self, "_next_id", 0), m) + 1
        return str(self._next_id)

    def new_para_id(self):
        used = getattr(self, "_para_ids", None)
        if used is None:
            used = self._para_ids = {e.get("{%s}paraId" % W14) for e in self.doc.iter() if e.get("{%s}paraId" % W14)}
        while True:
            v = "%08X" % random.randrange(0x10000000, 0x7FFFFFFF)
            if v not in used:
                used.add(v)
                return v

    # ---- output -------------------------------------------------------------
    def save(self, out_path):
        self.parts["word/document.xml"] = etree.tostring(
            self.doc, xml_declaration=True, encoding="UTF-8", standalone=True)
        if self._rels_dirty:
            self.parts["word/_rels/document.xml.rels"] = etree.tostring(
                self.rels, xml_declaration=True, encoding="UTF-8", standalone=True)
        with zipfile.ZipFile(out_path, "w", zipfile.ZIP_DEFLATED) as z:
            for name in self.order:
                z.writestr(zipfile.ZipInfo(name, date_time=(2026, 1, 1, 0, 0, 0)), self.parts[name],
                           compress_type=zipfile.ZIP_DEFLATED)
        return out_path


# ---- text -------------------------------------------------------------------
def _hidden(el):
    return any(a.tag in _HIDDEN for a in el.iterancestors())


def visible_runs(p):
    """Runs of one paragraph with tracked changes accepted (deleted / moved-away runs left out)."""
    return [r for r in p.iter(w("r")) if not _hidden(r)]


def para_text(p):
    """Visible text; w:br -> '\\n', w:tab -> '\\t'. The string replace_in_paragraph() matches on
    is this one without the '\\n' / '\\t' characters of w:br / w:tab (they are not w:t text)."""
    out = []
    for r in visible_runs(p):
        for c in r:
            if c.tag == w("t"):
                out.append(c.text or "")
            elif c.tag in (w("br"), w("cr")):
                out.append("\n")
            elif c.tag == w("tab"):
                out.append("\t")
    return "".join(out)


def cell_text(tc):
    return "\n".join(para_text(p) for p in tc.findall(w("p")))


def _set_t(t, text):
    t.text = text
    if text != text.strip() or "  " in text:
        t.set(XML_SPACE, "preserve")


def _fill_run(run, text, deleted=False):
    """Replace the content of a run (rPr kept) by text; '\\n' -> w:br."""
    for c in list(run):
        if c.tag != w("rPr"):
            run.remove(c)
    for i, line in enumerate(text.split("\n")):
        if i:
            etree.SubElement(run, w("br"))
        if line or i == 0:
            _set_t(etree.SubElement(run, w("delText" if deleted else "t")), line)


def _clean_rpr(rpr):
    """Copy of an rPr without revision markers (for a new run)."""
    if rpr is None:
        return None
    rpr = copy.deepcopy(rpr)
    for c in list(rpr):
        if c.tag in (w("ins"), w("del"), w("moveFrom"), w("moveTo"), w("rPrChange")):
            rpr.remove(c)
    return rpr if len(rpr) else None


def _new_run(rpr, text):
    r = etree.Element(w("r"))
    if rpr is not None:
        r.append(rpr)
    _fill_run(r, text)
    return r


def _track_el(docx, kind, tracked):
    e = etree.Element(w(kind))
    e.set(w("id"), docx.next_annotation_id())
    e.set(w("author"), tracked[0])
    e.set(w("date"), tracked[1])
    return e


def _mark_deleted(docx, run, tracked):
    """Wrap one run in w:del (inside its w:ins if it has one) and turn w:t into w:delText."""
    for t in run.findall(w("t")):
        t.tag = w("delText")
    d = _track_el(docx, "del", tracked)
    parent = run.getparent()
    parent.replace(run, d)
    d.append(run)
    return d


def _top_in_paragraph(el, p):
    while el.getparent() is not p:
        el = el.getparent()
    return el


def _split_wrappers_after(docx, node, p):
    """Close every wrapper (w:ins, w:moveTo, w:hyperlink) of `node` right after it: what
    followed inside the wrapper moves to a copy of the wrapper placed after it. Returns the
    child of `p` that now ends with `node`."""
    while node.getparent() is not p:
        wrap = node.getparent()
        rest = list(node.itersiblings())
        if rest:
            w2 = etree.Element(wrap.tag, attrib=dict(wrap.attrib))
            if w2.get(w("id")) is not None and wrap.tag != w("hyperlink"):
                w2.set(w("id"), docx.next_annotation_id())
            for c in rest:
                w2.append(c)
            wrap.addnext(w2)
        node = wrap
    return node


def set_cell_text(docx, tc, text, tracked=None):
    """Replace the text of one table cell. Kept: tcPr (width, shading, vMerge / gridSpan), the
    paragraph and its pPr (style, alignment, keepNext), and the rPr of the cell's first visible
    run (bold, highlight, colour, fonts), which the whole new text takes. '\\n' in text -> w:br.
    A cell whose runs differ in format (Response Table 3.17: symbol-font runs; cells with a
    hyperlink field) should be edited with replace_in_paragraph() instead.
    Returns the old text."""
    ps = tc.findall(w("p"))
    if len(ps) != 1:
        raise ValueError("cell has %d paragraphs" % len(ps))
    p = ps[0]
    if p.findall(".//" + w("fldChar")) or p.findall(".//" + w("hyperlink")) or p.findall(".//" + w("drawing")):
        raise ValueError("cell holds a field, hyperlink or drawing: use replace_in_paragraph")
    old = para_text(p)
    runs = visible_runs(p)
    if runs:
        rpr = runs[0].find(w("rPr"))
    else:                                   # empty cell: the paragraph-mark format
        rpr = p.find(w("pPr") + "/" + w("rPr"))
    if tracked is None:
        if runs:
            keep = runs[0]
            for r in runs[1:]:
                par = r.getparent()
                par.remove(r)
                if par.tag in (w("ins"), w("moveTo")) and len(par) == 0:
                    par.getparent().remove(par)
            _fill_run(keep, text)
        else:
            p.append(_new_run(_clean_rpr(rpr), text))
    else:
        last = None
        for r in runs:
            last = _top_in_paragraph(_mark_deleted(docx, r, tracked), p)
        ins = _track_el(docx, "ins", tracked)
        ins.append(_new_run(_clean_rpr(rpr), text))
        if last is not None:
            last.addnext(ins)
        else:
            p.append(ins)
    return old


def copy_row(docx, tbl, src, texts=None, after=True, tracked=None):
    """Insert a copy of row `src` (index) of `tbl` right after (or before) it and return the new
    w:tr. Formats come with the copy. Paragraph ids and revision ids are renewed; a copied
    vMerge 'restart' cell becomes a continuation cell (empty) when the copy goes below its
    source, so the merged block grows by one row; w:tblHeader is dropped from the copy.
    texts: optional list, one string per cell (None = leave the copied text)."""
    rows = tbl.findall(w("tr"))
    tr = rows[src]
    new = copy.deepcopy(tr)
    if new.get("{%s}paraId" % W14) is not None:          # the row element carries its own paraId
        new.set("{%s}paraId" % W14, docx.new_para_id())
    for p in new.iter(w("p")):
        if p.get("{%s}paraId" % W14) is not None:
            p.set("{%s}paraId" % W14, docx.new_para_id())
    for tag in ("ins", "del", "moveFrom", "moveTo"):
        for e in new.iter(w(tag)):
            if e.get(w("id")) is not None:
                e.set(w("id"), docx.next_annotation_id())
    trpr = new.find(w("trPr"))
    if trpr is not None:
        for e in trpr.findall(w("tblHeader")):
            trpr.remove(e)
    if after:
        tr.addnext(new)
        for tc in new.findall(w("tc")):
            vm = tc.find(w("tcPr") + "/" + w("vMerge"))
            if vm is not None and vm.get(w("val")) == "restart":
                del vm.attrib[w("val")]
                set_cell_text(docx, tc, "")
    else:
        tr.addprevious(new)
    if tracked is not None:
        if trpr is None:
            trpr = etree.Element(w("trPr"))
            new.insert(0, trpr)
        for e in trpr.findall(w("ins")):
            trpr.remove(e)
        trpr.append(_track_el(docx, "ins", tracked))
    if texts is not None:
        cells = new.findall(w("tc"))
        if len(texts) != len(cells):
            raise ValueError("%d texts for %d cells" % (len(texts), len(cells)))
        for tc, s in zip(cells, texts):
            if s is not None:
                set_cell_text(docx, tc, s)
    return new


def delete_row(tbl, idx):
    tr = tbl.findall(w("tr"))[idx]
    tbl.remove(tr)


# ---- sentences split over runs ------------------------------------------------
def _t_map(p):
    """[(w:t element, start offset)] over the visible runs, and the joined text."""
    items, pos, parts = [], 0, []
    for r in visible_runs(p):
        for t in r.findall(w("t")):
            s = t.text or ""
            items.append((t, pos))
            parts.append(s)
            pos += len(s)
    return items, "".join(parts)


def _split_run_at(t, offset):
    """Split the run of w:t `t` so that a new run starts at character `offset` of t (the new
    run takes a copy of the rPr, the text from `offset` on and every later child of the run).
    Returns the first w:t of the second run (t itself when no split is needed)."""
    run = t.getparent()
    s = t.text or ""
    kids = [c for c in run if c.tag != w("rPr")]
    if offset == 0 and kids[0] is t:
        return t
    rest = kids[kids.index(t) + 1:]
    if offset >= len(s) and not rest:
        return None
    r2 = etree.Element(w("r"), attrib=dict(run.attrib))
    rpr = run.find(w("rPr"))
    if rpr is not None:
        r2.append(copy.deepcopy(rpr))
    t2 = None
    if offset < len(s):
        t2 = etree.SubElement(r2, w("t"))
        _set_t(t2, s[offset:])
    for c in rest:
        r2.append(c)
    if offset == 0:
        run.remove(t)
    else:
        _set_t(t, s[:offset])
    run.addnext(r2)
    return t2


def replace_in_paragraph(docx, p, old, new, tracked=None, occurrence=None):
    """Replace `old` by `new` in the visible text of paragraph `p` even when `old` runs across
    several runs (and across w:ins wrappers). `old` must occur exactly once unless
    `occurrence` (0-based) is given; otherwise ValueError. w:br / w:tab are not part of the
    matched string (see para_text). Text outside the match keeps its runs and formats; the new
    text takes the format of the run in which the match starts. Returns the match offset."""
    if not old:
        raise ValueError("empty search string")
    items, full = _t_map(p)
    hits = [m.start() for m in re.finditer(re.escape(old), full)]
    if occurrence is None:
        if len(hits) != 1:
            raise ValueError("%d occurrences of %r" % (len(hits), old))
        start = hits[0]
    else:
        start = hits[occurrence]
    end = start + len(old)
    if tracked is None:
        placed = False
        for t, lo in items:
            s = t.text or ""
            hi = lo + len(s)
            a, b = max(start, lo), min(end, hi)
            if not placed:
                if lo <= start < hi:
                    _set_t(t, s[:start - lo] + new + s[b - lo:])
                    placed = True
            elif a < b:
                _set_t(t, s[:a - lo] + s[b - lo:])
        return start
    # tracked: put run boundaries at both ends of the match, then w:del the runs in between
    for t, lo in items:                                   # end first: offsets before it stay valid
        if lo < end <= lo + len(t.text or ""):
            _split_run_at(t, end - lo)
            break
    first_t = None
    for t, lo in _t_map(p)[0]:
        if lo <= start < lo + len(t.text or ""):
            first_t = _split_run_at(t, start - lo)
            break
    doomed = []
    for t, lo in _t_map(p)[0]:
        n = len(t.text or "")
        if n and start <= lo and lo + n <= end and t.getparent() not in doomed:
            doomed.append(t.getparent())
    if not doomed or first_t is None or doomed[0] is not first_t.getparent():
        raise RuntimeError("run split failed")
    rpr = _clean_rpr(doomed[0].find(w("rPr")))
    for r in doomed:
        last_del = _mark_deleted(docx, r, tracked)
    anchor = _split_wrappers_after(docx, last_del, p)
    if new:
        ins = _track_el(docx, "ins", tracked)
        ins.append(_new_run(rpr, new))
        anchor.addnext(ins)
    return start


# ---- images -----------------------------------------------------------------
def png_size(data):
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError("not a PNG")
    return struct.unpack(">II", data[16:24])


def _resize_drawing(d, wpx, hpx, max_height):
    ext = d.find(".//{%s}extent" % WP)
    cx = int(ext.get("cx"))
    cy = int(round(cx * hpx / float(wpx)))
    scale = 1.0
    if max_height and cy > max_height:
        scale = max_height / float(cy)
        cx, cy = int(round(cx * scale)), int(max_height)
    old = (int(ext.get("cx")), int(ext.get("cy")))
    ext.set("cx", str(cx))
    ext.set("cy", str(cy))
    aext = d.find(".//{%s}spPr/{%s}xfrm/{%s}ext" % (PIC, A, A))
    if aext is not None:                      # may differ slightly from wp:extent: keep its own width
        acx = int(round(int(aext.get("cx")) * scale))
        aext.set("cx", str(acx))
        aext.set("cy", str(int(round(acx * hpx / float(wpx)))))
    return old, (cx, cy)


def replace_image(docx, target, png_path, only=None, max_height=TEXT_HEIGHT_EMU - EMU_PER_INCH):
    """Put the PNG at png_path in place of the embedded image `target` ('media/image5.png').
    Width on the page is kept; height follows the new pixel aspect ratio in wp:extent and in
    pic:spPr/a:xfrm/a:ext (if that would exceed max_height, default 8 in, both are scaled down).
    only=None : the media part is overwritten, so EVERY drawing that uses it changes (ver3 shows
                image4 three times and image5, image6, image7 twice each).
    only=[i,..]: 0-based indexes into docx.drawings(); the PNG is added as a new media part with
                a new relationship and only these drawings are pointed at it.
    Returns [(drawing index, old (cx, cy), new (cx, cy))]."""
    with open(png_path, "rb") as f:
        data = f.read()
    wpx, hpx = png_size(data)
    draws = docx.drawings()
    if "word/" + target not in docx.parts:
        raise KeyError(target)
    if only is None:
        docx.parts["word/" + target] = data
        sel = [i for i, (_, _, tg) in enumerate(draws) if tg == target]
        new_rid = None
    else:
        for i in only:
            if draws[i][2] != target:
                raise ValueError("drawing %d uses %s" % (i, draws[i][2]))
        n = 1
        while "word/media/image%d.png" % n in docx.parts:
            n += 1
        name = "media/image%d.png" % n
        docx.parts["word/" + name] = data
        docx.order.append("word/" + name)
        ids = [int(r.get("Id")[3:]) for r in docx.rels if r.get("Id", "").startswith("rId") and r.get("Id")[3:].isdigit()]
        new_rid = "rId%d" % (max(ids) + 1)
        rel = etree.SubElement(docx.rels, "{%s}Relationship" % PKG_REL)
        rel.set("Id", new_rid)
        rel.set("Type", IMAGE_REL)
        rel.set("Target", name)
        docx._rels_dirty = True
        sel = list(only)
    out = []
    for i in sel:
        d = draws[i][0]
        if new_rid:
            d.find(".//{%s}blip" % A).set("{%s}embed" % R, new_rid)
        old, newsz = _resize_drawing(d, wpx, hpx, max_height)
        out.append((i, old, newsz))
    return out


# ---- checks -----------------------------------------------------------------
def check(path):
    """Zip CRC test, XML well-formedness of every .xml / .rels part, every image relationship
    resolves to a part, every drawing's rId exists, revision / bookmark ids unique, PNG default
    content type present. Returns a dict of counts; raises on a failure."""
    with zipfile.ZipFile(path) as z:
        bad = z.testzip()
        if bad:
            raise ValueError("bad CRC: " + bad)
        names = z.namelist()
        if len(names) != len(set(names)):
            raise ValueError("duplicate zip entries")
        trees = {}
        for n in names:
            if n.endswith(".xml") or n.endswith(".rels"):
                trees[n] = etree.fromstring(z.read(n))
        rels = {r.get("Id"): r for r in trees["word/_rels/document.xml.rels"]}
        if len(rels) != len(trees["word/_rels/document.xml.rels"]):
            raise ValueError("duplicate relationship ids")
        for r in rels.values():
            if r.get("TargetMode") != "External" and "word/" + r.get("Target") not in names:
                raise ValueError("relationship target missing: " + r.get("Target"))
        doc = trees["word/document.xml"]
        nd = 0
        for d in doc.iter(w("drawing")):
            nd += 1
            rid = d.find(".//{%s}blip" % A).get("{%s}embed" % R)
            if rid not in rels:
                raise ValueError("drawing rId missing: " + rid)
            png_size(z.read("word/" + rels[rid].get("Target")))
            e = d.find(".//{%s}extent" % WP)
            if int(e.get("cx")) > TEXT_WIDTH_EMU or int(e.get("cy")) > TEXT_HEIGHT_EMU:
                raise ValueError("drawing larger than the text area")
        ids = [e.get(w("id")) for tag in ("ins", "del", "moveFrom", "moveTo") for e in doc.iter(w(tag))]
        if len(ids) != len(set(ids)):
            raise ValueError("duplicate revision ids")
        for tbl in doc.iter(w("tbl")):
            ncol = len(tbl.find(w("tblGrid")).findall(w("gridCol")))
            for tr in tbl.findall(w("tr")):
                span = 0
                for tc in tr.findall(w("tc")):
                    gs = tc.find(w("tcPr") + "/" + w("gridSpan"))
                    span += int(gs.get(w("val"))) if gs is not None else 1
                    if not tc.findall(w("p")):
                        raise ValueError("table cell without a paragraph")
                if span != ncol:
                    raise ValueError("row spans %d of %d grid columns" % (span, ncol))
        if b'Extension="png"' not in z.read("[Content_Types].xml"):
            raise ValueError("no PNG content type")
        return {"parts": len(names), "xml_parts": len(trees), "drawings": nd,
                "tables": len(list(doc.iter(w("tbl")))), "revisions": len(ids),
                "paragraphs": len(list(doc.iter(w("p"))))}
