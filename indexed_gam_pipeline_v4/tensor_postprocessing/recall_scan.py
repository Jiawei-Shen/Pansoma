"""The truth recall's scan of the filtered-candidate streams in parallel; label_run with that scan.

truth_labels.scan_filtered reads <kind>/filtered_candidates.ndjson in one stream (0.1-6.4 h per genome set; the
ONT streams are hundreds of GB). scan_filtered here returns the same dict:
  * the file is cut into byte ranges that start at line starts, so every range holds whole lines, exactly the
    lines of `for line in stream`;
  * WORKERS processes scan the ranges with truth_labels.CANDIDATE and REASONS (the same per-line searches) and keep
    the lines whose candidate_id has its 64-bit BLAKE2b hash in the sorted hashes of the truth keys;
  * the parent keeps the hits whose candidate_id is a truth key and applies them in file order, so the last line of
    an id gives its reasons and the dict's order is that of first appearance.
The workers come from a forkserver (Slurm counts the RSS of every process: forked copies of a parent that holds the
truth sets would each count them again) and get the hash array once, as their initializer's argument. As with
spawn, a worker imports the caller's main module (`python -m <package>.orchestrate` and `python -c` are
import-safe), so a script that calls this keeps its work under `if __name__ == "__main__":`. When the pool cannot
work (a main module read from stdin, a TMPDIR too long for the forkserver's socket, a worker killed), the scan
falls back to truth_labels.scan_filtered: the same dict, read in one stream.

label_run is truth_labels.label_run verbatim; the scan_filtered it calls is the one here, every other name it uses
is truth_labels' own (tests/test_recall_scan.py pins both). The labels, labels.manifest.json (rules_sha256 is still
the SHA-256 of truth_labels.py) and the recall reports are those of truth_labels.label_run.
"""
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
import hashlib
import multiprocessing
import os
from pathlib import Path
import sys

import numpy as np

from ..common import sha256_file, write_json
from .reference_path import ReferencePath
from . import truth_labels
from .truth_labels import CANDIDATE, REASONS, Locator, TruthSet, label_directory, recall

WORKERS = 8  # processes of one scan: 8 streams reach a node's read rate on /scratch, 16 are no faster
PARTS = 4  # byte ranges per worker (ranges are handed out as workers free up)
BLOCK = 16 << 20  # bytes per read
_TABLE = None  # in a worker: the sorted truth-key hashes (uint64)


def key_hash(key):
    """64-bit BLAKE2b of a candidate_id's bytes."""
    return int.from_bytes(hashlib.blake2b(key, digest_size=8).digest(), "little")


def byte_ranges(path, parts):
    """[(start, stop)] cutting `path` into at most `parts` ranges of about equal size; every start is a line start
    (the first line start at or after each even cut)."""
    size = os.path.getsize(path)
    cuts = [0]
    with open(path, "rb") as stream:
        for k in range(1, parts):
            at = size * k // parts
            if at > cuts[-1]:
                stream.seek(at - 1)
                stream.readline()  # through the end of the line holding byte at - 1
                cuts.append(stream.tell())
    cuts.append(size)
    return [(a, b) for a, b in zip(cuts, cuts[1:]) if b > a]


def load_table(table):
    """Pool initializer: the worker's truth-key hashes."""
    global _TABLE
    _TABLE = table


def match(lines, hits):
    """Append (candidate_id, reasons bytes or None) of each line whose candidate_id hash is in _TABLE to `hits`."""
    keys, found = [], []
    for line in lines:
        m = CANDIDATE.search(line)
        if m:
            key = m.group(1)
            if not key.isascii():
                key.decode()  # the serial scan decodes every id: invalid UTF-8 raises the same error here
            keys.append(key)
            found.append(line)
    if not keys:
        return
    hashes = np.fromiter((key_hash(k) for k in keys), dtype=np.uint64, count=len(keys))
    at = np.searchsorted(_TABLE, hashes)
    known = at < len(_TABLE)
    known[known] = _TABLE[at[known]] == hashes[known]
    for i in np.flatnonzero(known).tolist():
        r = REASONS.search(found[i])
        hits.append((keys[i].decode(), r.group(1) if r else None))


def scan_range(path, start, stop, block=BLOCK):
    """Worker job: match() over the lines of bytes [start, stop) of `path` (start a line start), in file order."""
    hits, carry = [], b""
    with open(path, "rb") as stream:
        stream.seek(start)
        left = stop - start
        while left:
            data = stream.read(min(block, left))
            if not data:
                raise ValueError(f"{path} ends before byte {stop}")
            left -= len(data)
            lines = (carry + data).split(b"\n")
            carry = lines.pop()  # the line not ended yet (b"" after a newline)
            match(lines, hits)
    if carry:  # the file's last line, without a newline
        match([carry], hits)
    return hits


def scan_filtered(directory, keys, workers=WORKERS, parts=None, block=BLOCK):
    """{candidate_id: reasons} for filtered candidates whose id is a truth key: truth_labels.scan_filtered's dict
    (the same items in the same order), the file's byte ranges (`parts`, default PARTS per worker) scanned by
    `workers` processes."""
    found = {}
    path = Path(directory) / "filtered_candidates.ndjson"
    if not path.exists():
        return found
    jobs = byte_ranges(path, parts or PARTS * workers)
    if not jobs:
        return found
    hashes = np.fromiter((key_hash(k.encode()) for k in keys), dtype=np.uint64, count=len(keys))
    hashes.sort()  # in place (np.unique holds several copies: +1.3 GiB at 21 M keys); a repeated hash still matches
    try:
        with ProcessPoolExecutor(min(workers, len(jobs)), mp_context=multiprocessing.get_context("forkserver"),
                                 initializer=load_table, initargs=(hashes,)) as pool:
            results = list(pool.map(scan_range, *zip(*[(str(path), a, b, block) for a, b in jobs])))
    except (BrokenProcessPool, OSError) as error:  # the pool, not the data: the serial scan reads it (or raises)
        print(f"recall_scan: parallel scan of {path} failed ({type(error).__name__}: {error}); scanning it in one "
              f"stream", file=sys.stderr, flush=True)
        return truth_labels.scan_filtered(directory, keys)
    for hits in results:
        for cid, r in hits:  # ranges in file order: the last line of an id wins
            if cid in keys:
                found[cid] = [x.strip().strip('"') for x in r.decode().split(",")] if r is not None else []
    return found


# truth_labels.label_run, verbatim (scan_filtered above).
def label_run(tensors, kinds, reference_path, fasta, somatic_vcf, somatic_bed, germline_vcf, germline_bed, truth_dir,
              recall_dir=None, snv_min_af=None, indel_min_af=None):
    """Build both truth sets, label every merged autosome directory of `kinds`, write recall reports.

    Truth tables (<set>.graph.tsv) depend only on the graph and the VCFs and go to `truth_dir`;
    recall reports depend on the run and go to `recall_dir` (default: `tensors`).
    """
    truth_dir = Path(truth_dir)
    truth_dir.mkdir(parents=True, exist_ok=True)
    recall_dir = Path(recall_dir or tensors)
    recall_dir.mkdir(parents=True, exist_ok=True)
    path = ReferencePath(reference_path)
    locator = Locator(path)
    somatic = TruthSet("somatic", somatic_vcf, somatic_bed, fasta, locator)
    germline = TruthSet("germline", germline_vcf, germline_bed, fasta, locator)
    confident = somatic.bed.intersect(germline.bed)
    for truth in (somatic, germline):
        truth.write_table(truth_dir / f"{truth.name}.graph.tsv")
    provenance = dict(truth={t.name: dict(vcf=str(t.vcf), vcf_sha256=sha256_file(t.vcf), bed=str(t.bed_path),
                                          bed_sha256=sha256_file(t.bed_path), stats=dict(t.stats),
                                          table=str(truth_dir / f"{t.name}.graph.tsv"))
                             for t in (somatic, germline)},
                      confident_region=dict(definition="somatic BED ∩ germline BED", bp=confident.size()),
                      fasta=str(fasta), reference_path=str(path.directory))
    reports, matched = {}, {somatic.name: defaultdict(list), germline.name: defaultdict(list)}
    for kind in kinds:
        reports[kind], found = label_directory(Path(tensors) / kind, somatic, germline, path, confident, provenance, fasta,
                                               snv_min_af, indel_min_af)
        for name, items in found.items():
            for tid, hits in items.items():
                matched[name][tid] += hits
    keys = {k for truth in (somatic, germline) for a in truth.alleles for k in a["keys"]}
    filtered = {}
    for kind in kinds:
        filtered.update(scan_filtered(Path(tensors) / kind, keys))
    recall_report = {t.name: recall(t, matched[t.name], filtered, recall_dir) for t in (somatic, germline)}
    write_json(recall_dir / "truth_recall.json", recall_report)
    return dict(labels={k: dict(totals=r["totals"], reasons=r["reasons"]) for k, r in reports.items()},
                recall=recall_report, truth_stats={t.name: dict(t.stats) for t in (somatic, germline)})
