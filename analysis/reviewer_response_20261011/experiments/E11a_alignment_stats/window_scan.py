"""E11a window design: the reads of fixed GRCh38 windows (windows.tsv) in one arm, plus the arm's unmapped reads.

GAM arm: per window, the GRCh38-path nodes whose path start lies in [start - FLANK, end + FLANK), plus every node ID
between two consecutive (sorted) such IDs at most GAP apart (the off-reference nodes of the bubbles in between);
the GAI runs of the bins holding those IDs are walked (as IndexedGam._walk_hits does) and every record visiting one
of the nodes is kept. Unmapped reads: the sorted GAM's head (the runs of GAI bin [0, 1]: records without a path).
BAM arm: per window fetch(chrom, start - FLANK, end + FLANK); primary records kept, supplementary records counted
per read; unmapped reads: flag 0x4 records there plus fetch('*').
Name rule (--mod M --keep K): only names with zlib.crc32(name) % M < K are kept (both mates), in every arm
(short reads: 10 / 1 within windows; long reads: 100 / 10, the rule of scan_gam.py's genome-wide sample).
Flank: --short 1 kb, else 25 kb. Windows: windows.tsv (core, index 0-47) and, with --extra, windows_extra.tsv
(index 48 + i).
Per-read columns: scan_gam.ROW (GAM) / scan_bam.BROW (BAM) plus `window` (index into windows.tsv).
Output: <out>.npz with rows, unmapped (keys name\\tmate of unmapped reads that pass the name rule), contigs; and
<out>.summary.json (per window: records walked / kept; head: records, unmapped records counted).

Run: python window_scan.py --kind gam|bam --path P [--refpath R] --out PREFIX [--short] [--procs 8]
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import multiprocessing
import os
from pathlib import Path
import sys
import time
import zlib

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, '/scratch/jshen/Github/Pansoma')
import scan_gam  # noqa: E402
from scan_gam import ROW, head, read_stats, load_ref  # noqa: E402
from indexed_gam_pipeline_v4.gam_reader import BgzfReader, IndexedGam, decode, read_gai  # noqa: E402

GAP = 64
FLANK_LONG, FLANK_SHORT = 25_000, 1_000
WROW = np.dtype(ROW.descr + [('window', 'u2')])
READERS = {}


def windows(extra=False):
    out = []
    for name in ['windows.tsv'] + (['windows_extra.tsv'] if extra else []):
        with open(HERE / name) as f:
            next(f)
            for line in f:
                w, c, s, e, k = line.rstrip('\n').split('\t')[:5]
                out.append((w, c, int(s), int(e), k))
    return out


def window_nodes(refdir, chrom, start, end, flank):
    d = Path(refdir)
    meta = json.load(open(d / 'meta.json'))
    ci = [c for c in meta['contigs'] if c['name'] == chrom][0]
    sl = slice(ci['offset'], ci['offset'] + ci['nodes'])
    ps = np.load(d / 'path_starts.npy', mmap_mode='r')[sl]
    pn = np.load(d / 'path_nodes.npy', mmap_mode='r')[sl]
    a, b = np.searchsorted(ps, start - flank), np.searchsorted(ps, end + flank)
    ids = np.unique(np.asarray(pn[a:b]))
    if not len(ids):
        return ids
    gaps = np.diff(ids)
    fill = [np.arange(x + 1, y) for x, y, g in zip(ids[:-1], ids[1:], gaps) if 1 < g <= GAP]
    return np.unique(np.concatenate([ids] + fill)) if fill else ids


RULE = dict(mod=1, keep=1)


def keep_name(name, short=None):
    return zlib.crc32(name) % RULE['mod'] < RULE['keep']


def gam_window(job):
    path, refdir, wi, chrom, start, end, short = job
    ref = load_ref(refdir)
    flank = FLANK_SHORT if short else FLANK_LONG
    nodes = window_nodes(refdir, chrom, start, end, flank)
    t0 = time.time()
    reader = READERS.get(path) or IndexedGam(path, path + '.gai', cache_bytes=1)
    ranges = reader.ranges(nodes.tolist())
    wanted = set(nodes.tolist())
    rows, walked = [], 0
    with BgzfReader(path) as stream:
        for s, e in ranges:
            stream.limit(e)
            stream.seek(s)
            while stream.tell() < e:
                msgs = stream.read_group()
                if msgs is None:
                    break
                for raw in msgs:
                    walked += 1
                    name, has_path, _ = head(raw)
                    if not has_path or not keep_name(name):
                        continue
                    a = decode(raw)
                    if not any(m.position.node_id in wanted for m in a.path.mapping):
                        continue
                    r = read_stats(a, ref)
                    r['window'] = wi
                    rows.append(r)
    out = np.zeros(len(rows), dtype=WROW)
    for k in WROW.names:
        if rows:
            out[k] = [r.get(k, 0) for r in rows]
    # one read can be in two GAI runs only once (runs are merged), but keep the first copy to be safe
    return dict(kind='window', window=wi, rows=out, walked=walked, nodes=int(len(nodes)), runs=len(ranges),
                seconds=round(time.time() - t0, 1))


def gam_head(job):
    path, short = job
    t0 = time.time()
    size = Path(path).stat().st_size
    _, bins, _ = read_gai(path + '.gai', size)
    runs = sorted({(s, e) for lo, hi, rr in bins if lo == 0 and hi <= 1 for s, e in rr})
    keys, records, unmapped = [], 0, 0
    with BgzfReader(path) as stream:
        for s, e in runs:
            stream.limit(e)
            stream.seek(s)
            while stream.tell() < e:
                msgs = stream.read_group()
                if msgs is None:
                    break
                for raw in msgs:
                    records += 1
                    name, has_path, _ = head(raw)
                    if has_path:
                        continue
                    unmapped += 1
                    if not keep_name(name, short):
                        continue
                    a = decode(raw)
                    mate = 1 if a.fragment_next.name else 2 if a.fragment_prev.name else 0
                    keys.append(name + b'\t' + str(mate).encode())
    return dict(kind='head', keys=keys, records=records, unmapped=unmapped, runs=len(runs),
                seconds=round(time.time() - t0, 1))


def bam_window(job):
    import pysam
    from scan_bam import BROW, primary_row
    path, wi, chrom, start, end, short = job
    flank = FLANK_SHORT if short else FLANK_LONG
    t0 = time.time()
    rows, supp, unmapped, walked = [], {}, [], 0
    with pysam.AlignmentFile(path, 'rb', threads=2) as f:
        for r in f.fetch(chrom, max(0, start - flank), end + flank):
            walked += 1
            flag = r.flag
            if flag & 0x100:
                continue
            name = r.query_name.encode()
            if not keep_name(name, short):
                continue
            mate = 1 if flag & 0x40 else 2 if flag & 0x80 else 0
            key = name + b'\t' + str(mate).encode()
            if flag & 0x800:
                s = supp.setdefault(key, [0, 0, set()])
                s[0] += 1
                s[1] += r.query_alignment_length
                s[2].add(r.reference_id)
                continue
            if flag & 0x4:
                unmapped.append(key)
                continue
            row = primary_row(r)
            row['window'] = wi
            rows.append(row)
    dt = np.dtype(BROW.descr + [('window', 'u2')])
    out = np.zeros(len(rows), dtype=dt)
    for k in dt.names:
        if rows:
            out[k] = [r.get(k, 0) for r in rows]
    for i in range(len(out)):
        s = supp.get(out['name'][i] + b'\t' + str(int(out['mate'][i])).encode())
        if s:
            out['supp'][i], out['supp_bp'][i] = s[0], s[1]
            other = len(s[2] - {int(out['chrom'][i])})
            out['supp_other_chrom'][i] = other
            out['nchrom'][i] = 1 + other
    return dict(kind='window', window=wi, rows=out, walked=walked, unmapped_keys=unmapped,
                supp_only=len(set(supp) - {out['name'][i] + b'\t' + str(int(out['mate'][i])).encode() for i in range(len(out))}),
                seconds=round(time.time() - t0, 1))


def bam_head(job):
    import pysam
    path, short = job
    t0 = time.time()
    keys, records = [], 0
    with pysam.AlignmentFile(path, 'rb', threads=2) as f:
        for r in f.fetch(contig='*', until_eof=True) if f.has_index() else []:
            records += 1
            if r.flag & 0x900:
                continue
            name = r.query_name.encode()
            if not keep_name(name, short):
                continue
            mate = 1 if r.flag & 0x40 else 2 if r.flag & 0x80 else 0
            keys.append(name + b'\t' + str(mate).encode())
    return dict(kind='head', keys=keys, records=records, unmapped=records, runs=1, seconds=round(time.time() - t0, 1))


def run(job):
    kind = job[0]
    return {'gw': gam_window, 'gh': gam_head, 'bw': bam_window, 'bh': bam_head}[kind](job[1:])


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--kind', choices=['gam', 'bam'], required=True)
    p.add_argument('--path', required=True)
    p.add_argument('--refpath')
    p.add_argument('--out', required=True)
    p.add_argument('--short', action='store_true')
    p.add_argument('--procs', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', 4)))
    p.add_argument('--windows', help='test: comma-separated window indices')
    p.add_argument('--no-head', action='store_true')
    p.add_argument('--extra', action='store_true', help='also windows_extra.tsv')
    p.add_argument('--mod', type=int, required=True)
    p.add_argument('--keep', type=int, required=True)
    args = p.parse_args()
    RULE.update(mod=args.mod, keep=args.keep)
    t0 = time.time()
    wins = windows(args.extra)
    pick = range(len(wins)) if not args.windows else [int(x) for x in args.windows.split(',')]
    jobs = []
    if not args.no_head:
        jobs.append(('gh', args.path, args.short) if args.kind == 'gam' else ('bh', args.path, args.short))
    for wi in pick:
        w, c, s, e, k = wins[wi]
        if args.kind == 'gam':
            jobs.append(('gw', args.path, args.refpath, wi, c, s, e, args.short))
        else:
            jobs.append(('bw', args.path, wi, c, s, e, args.short))
    if args.kind == 'gam':
        load_ref(args.refpath)
        READERS[args.path] = IndexedGam(args.path, args.path + '.gai', cache_bytes=1)
    with ProcessPoolExecutor(args.procs, mp_context=multiprocessing.get_context('fork')) as pool:
        results = list(pool.map(run, jobs))
    wres = [r for r in results if r['kind'] == 'window']
    hres = [r for r in results if r['kind'] == 'head']
    rows = np.concatenate([r['rows'] for r in wres])
    # a read found in two windows (only possible for very long reads near two close windows): keep the first
    keys = np.char.add(np.char.add(rows['name'], b'\t'), rows['mate'].astype('S1'))
    _, first = np.unique(keys, return_index=True)
    dup = len(rows) - len(first)
    rows = rows[np.sort(first)]
    unm = [k for r in hres for k in r['keys']] + [k for r in wres for k in r.get('unmapped_keys', [])]
    if args.kind == 'gam':
        contigs = scan_gam.REF['contigs']
    else:
        import pysam
        with pysam.AlignmentFile(args.path, 'rb') as f:
            contigs = list(f.references)
    np.savez_compressed(args.out + '.npz', rows=rows, unmapped=np.array(sorted(set(unm)), dtype='S90'),
                        contigs=np.array(json.dumps(contigs)))
    summ = dict(kind=args.kind, path=args.path, refpath=args.refpath, short=args.short, mod=args.mod, keep=args.keep,
                extra=args.extra, windows=len(wres),
                complete=not args.windows and not args.no_head, rows=int(len(rows)), duplicate_rows_dropped=int(dup),
                unmapped_keys=len(set(unm)),
                head=[{k: v for k, v in r.items() if k != 'keys'} for r in hres],
                per_window=[{k: v for k, v in r.items() if k not in ('rows', 'unmapped_keys')} for r in wres],
                wall_seconds=round(time.time() - t0, 1))
    json.dump(summ, open(args.out + '.summary.json', 'w'), indent=1)
    print(json.dumps({k: v for k, v in summ.items() if k != 'per_window'}))


if __name__ == '__main__':
    main()
