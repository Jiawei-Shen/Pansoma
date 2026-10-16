"""E11a: one full parallel walk over a sorted GAM; whole-file counts for every record and per-read statistics for a
fixed name-hash sample of reads (the same reads in every arm, because every arm aligned the same FASTQ).

Selection: zlib.crc32(read name) % MOD < KEEP (mates share the name, so both mates of a selected pair are kept).
Whole file (every record): records, records with a path, MAPQ histogram (all records / records with a path).
Selected reads: decoded with the pipeline's vg_pb2 (indexed_gam_pipeline_v4.gam_reader.decode, imported, not changed);
per read:
  mate        0 unpaired, 1 = fragment_next set (first mate), 2 = fragment_prev set (second mate)
  read_len    len(sequence); mapped = path has >= 1 mapping; mapq
  match / mm  bases of edits with from_length == to_length without / with sequence; mm_ev = runs of mismatch edits
  ins / ins_ev, dl / del_ev   inserted / deleted bases and events (consecutive indel edits of one kind, also across
              node boundaries, are one event); an insertion that is the first edit of the first mapping or the last
              edit of the last mapping is a soft clip (sc_l / sc_r), as vg writes soft clips
  cx          bases of complex edits (from != to, both > 0), counted as mismatch min(from,to) + indel |from - to|
  off_visits / off_bp / off_edited   mappings on nodes off the GRCh38 path (reference_path chrom == -1), the graph
              bases they cover (sum of from_length) and how many of them carry any edit
  ref_visits  mappings on GRCh38-path nodes
  chrom, gstart, gend   GRCh38 placement: the contig holding most anchored (GRCh38-path) graph bases; span = min
              start0 .. max start0 + node length over that contig's anchored nodes (node granularity)
  nchrom      distinct contigs among the anchored nodes
  identity    giraffe's Alignment.identity (cross-check only)
Output: <out>.npz (rows: structured array; contigs: JSON list for chrom) and <out>.summary.json.

Run: python scan_gam.py --gam G --refpath DIR --out PREFIX [--mod 1000 --keep 1] [--procs 16] [--max-segments N]
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

sys.path.insert(0, '/scratch/jshen/Github/Pansoma')
from indexed_gam_pipeline_v4 import gam_reader  # noqa: E402
from indexed_gam_pipeline_v4.gam_reader import BgzfReader, IndexedGam, decode  # noqa: E402

ROW = np.dtype([('name', 'S80'), ('mate', 'u1'), ('read_len', 'u4'), ('mapped', 'u1'), ('mapq', 'u1'),
                ('nmap', 'u4'), ('match', 'u4'), ('mm', 'u4'), ('mm_ev', 'u4'), ('ins', 'u4'), ('ins_ev', 'u4'),
                ('dl', 'u4'), ('del_ev', 'u4'), ('cx', 'u4'), ('sc_l', 'u4'), ('sc_r', 'u4'),
                ('off_visits', 'u4'), ('off_bp', 'u4'), ('off_edited', 'u4'), ('ref_visits', 'u4'),
                ('chrom', 'i2'), ('gstart', 'i8'), ('gend', 'i8'), ('nchrom', 'u1'), ('identity', 'f4')])

REF = {}


def load_ref(refdir):
    if not REF:
        d = Path(refdir)
        REF['chrom'] = np.load(d / 'chrom.npy', mmap_mode='r')
        REF['start0'] = np.load(d / 'start0.npy', mmap_mode='r')
        REF['lengths'] = np.load(d / 'lengths.npy', mmap_mode='r')
        REF['contigs'] = [c['name'] for c in json.load(open(d / 'meta.json'))['contigs']]
    return REF


def head(raw):
    """(name, has_path, mapq) from the wire format: fields 1-5 come first in vg's serialization."""
    n, i, name, has_path, mapq = len(raw), 0, b'', False, 0
    while i < n:
        key = raw[i]
        i += 1
        if key & 0x80:  # field number >= 16: past the fields we need
            break
        field, wire = key >> 3, key & 7
        if field > 5:
            break
        if wire == 2:
            size = shift = 0
            while True:
                b = raw[i]
                i += 1
                size |= (b & 0x7F) << shift
                if not b & 0x80:
                    break
                shift += 7
            if field == 3:
                name = raw[i:i + size]
            elif field == 2:
                has_path = size > 0
            i += size
        elif wire == 0:
            v = shift = 0
            while True:
                b = raw[i]
                i += 1
                v |= (b & 0x7F) << shift
                if not b & 0x80:
                    break
                shift += 7
            if field == 5:
                mapq = v
        elif wire == 1:
            i += 8
        elif wire == 5:
            i += 4
        else:
            break
    return name, has_path, mapq


def read_stats(a, ref):
    maps = a.path.mapping
    mate = 1 if a.fragment_next.name else 2 if a.fragment_prev.name else 0
    row = dict(name=a.name.encode()[:80], mate=mate, read_len=len(a.sequence), mapped=0, mapq=a.mapping_quality,
               nmap=len(maps), identity=a.identity, chrom=-1, gstart=-1, gend=-1)
    if not maps:
        return row
    row['mapped'] = 1
    match = mm = mm_ev = ins = ins_ev = dl = del_ev = cx = 0
    prev = 0
    nodes, fl_node, ed_node = [], [], []
    for m in maps:
        f_sum = 0
        edited = False
        for e in m.edit:
            f, t = e.from_length, e.to_length
            if f == t:
                if e.sequence:
                    mm += f
                    if prev != 1:
                        mm_ev += 1
                    prev = 1
                    edited = True
                else:
                    match += f
                    prev = 0
            elif t == 0:
                dl += f
                if prev != 3:
                    del_ev += 1
                prev = 3
                edited = True
            elif f == 0:
                ins += t
                if prev != 2:
                    ins_ev += 1
                prev = 2
                edited = True
            else:
                cx += max(f, t)
                mm += min(f, t)
                if f > t:
                    dl += f - t
                    del_ev += 1
                else:
                    ins += t - f
                    ins_ev += 1
                prev = 4
                edited = True
            f_sum += f
        nodes.append(m.position.node_id)
        fl_node.append(f_sum)
        ed_node.append(edited)
    sc_l = sc_r = 0
    first = maps[0].edit[0] if len(maps[0].edit) else None
    last = maps[-1].edit[-1] if len(maps[-1].edit) else None
    if first is not None and first.from_length == 0 and first.to_length > 0:
        sc_l = first.to_length
    if last is not None and last.from_length == 0 and last.to_length > 0 and not (len(maps) == 1 and len(maps[0].edit) == 1):
        sc_r = last.to_length
    # soft clips were counted as insertion events / bases above
    if sc_l:
        ins -= sc_l
        ins_ev -= 1
    if sc_r:
        ins -= sc_r
        ins_ev -= 1
    nodes = np.asarray(nodes, dtype=np.int64)
    fl_node = np.asarray(fl_node, dtype=np.int64)
    ed_node = np.asarray(ed_node, dtype=bool)
    ch = ref['chrom'][nodes]
    off = ch < 0
    row.update(match=match, mm=mm, mm_ev=mm_ev, ins=max(ins, 0), ins_ev=max(ins_ev, 0), dl=dl, del_ev=del_ev, cx=cx,
               sc_l=sc_l, sc_r=sc_r, off_visits=int(off.sum()), off_bp=int(fl_node[off].sum()),
               off_edited=int((off & ed_node).sum()), ref_visits=int((~off).sum()))
    if (~off).any():
        cn, nb = ch[~off], fl_node[~off]
        uniq, inv = np.unique(cn, return_inverse=True)
        bases = np.bincount(inv, weights=nb + 1)
        c = uniq[int(np.argmax(bases))]
        sel = nodes[ch == c]
        st = ref['start0'][sel]
        row.update(chrom=int(c), gstart=int(st.min()), gend=int((st + ref['lengths'][sel]).max()), nchrom=len(uniq))
    return row


def scan_segment(job):
    gam, start, end, mod, keep, refdir = job
    ref = load_ref(refdir)
    hist_all = np.zeros(256, dtype=np.int64)
    hist_mapped = np.zeros(256, dtype=np.int64)
    records = with_path = groups = 0
    rows = []
    t0 = time.time()
    crc = zlib.crc32
    with BgzfReader(gam) as reader:
        target = None if end is None else reader.position(end)
        reader.seek(start)
        while True:
            here = reader.position()
            if here is None or here == target:
                break
            if target is not None and here > target:
                raise ValueError(f'walk from {start} passed {end}')
            msgs = reader.read_group()
            if msgs is None:
                break
            groups += 1
            for raw in msgs:
                name, has_path, mapq = head(raw)
                records += 1
                hist_all[mapq] += 1
                if has_path:
                    with_path += 1
                    hist_mapped[mapq] += 1
                if crc(name) % mod < keep:
                    rows.append(read_stats(decode(raw), ref))
    out = np.zeros(len(rows), dtype=ROW)
    for k in ROW.names:
        if rows:
            out[k] = [r.get(k, 0) for r in rows]
    return dict(start=start, records=records, with_path=with_path, groups=groups, hist_all=hist_all,
                hist_mapped=hist_mapped, rows=out, seconds=time.time() - t0)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--gam', required=True)
    p.add_argument('--gai')
    p.add_argument('--refpath', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--mod', type=int, default=1000)
    p.add_argument('--keep', type=int, default=1)
    p.add_argument('--procs', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', 4)))
    p.add_argument('--segments-per-proc', type=int, default=8)
    p.add_argument('--max-segments', type=int, default=0, help='test: scan only the first N segments')
    p.add_argument('--skip-segments', type=int, default=0, help='test: skip the first N segments')
    args = p.parse_args()
    t0 = time.time()
    gai = args.gai or args.gam + '.gai'
    bins = IndexedGam(args.gam, gai, cache_bytes=1).bins
    size = Path(args.gam).stat().st_size
    cuts = gam_reader.segments(size, bins, max(1, args.procs * args.segments_per_proc))
    if args.max_segments:
        cuts = cuts[args.skip_segments:args.skip_segments + args.max_segments]
    jobs = [(args.gam, s, e, args.mod, args.keep, args.refpath) for s, e in cuts]
    load_ref(args.refpath)
    if args.procs > 1:
        with ProcessPoolExecutor(args.procs, mp_context=multiprocessing.get_context('fork')) as pool:
            results = list(pool.map(scan_segment, jobs))
    else:
        results = [scan_segment(j) for j in jobs]
    rows = np.concatenate([r['rows'] for r in results]) if results else np.zeros(0, dtype=ROW)
    contigs = REF['contigs']
    np.savez_compressed(args.out + '.npz', rows=rows, contigs=np.array(json.dumps(contigs)))
    summ = dict(gam=args.gam, gam_bytes=size, refpath=args.refpath, mod=args.mod, keep=args.keep,
                segments=len(jobs), complete=not args.max_segments,
                records=sum(r['records'] for r in results), with_path=sum(r['with_path'] for r in results),
                groups=sum(r['groups'] for r in results), selected=int(len(rows)),
                mapq_hist_all=np.sum([r['hist_all'] for r in results], axis=0).tolist(),
                mapq_hist_with_path=np.sum([r['hist_mapped'] for r in results], axis=0).tolist(),
                segment_seconds=[round(r['seconds'], 1) for r in results], wall_seconds=round(time.time() - t0, 1))
    json.dump(summ, open(args.out + '.summary.json', 'w'), indent=1)
    print(json.dumps({k: v for k, v in summ.items() if not k.startswith('mapq_hist') and k != 'segment_seconds'}))


if __name__ == '__main__':
    main()
