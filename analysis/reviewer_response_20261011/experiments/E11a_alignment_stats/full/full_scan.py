"""E11a full data: per-read alignment counts of one whole GAM or BAM, by GRCh38 placement contig.

GAM: the sorted GAM is cut into segments at GAI group starts (gam_reader.segments; ends normalized with
BgzfReader.position), each segment is counted by gam_fullstats (C, same per-read rules as scan_gam.read_stats; checked
against it on Illumina, ONT and unmapped-head segments, 2026-10-19), and the counts are summed.
BAM: bam_fullstats (C) per contig chr1-chr22 (primary records only).
Output: <out>.json: per contig name {records, perfect, sub, indel, clip, read_bases, aligned_bases, sub_bases,
ins_bases, del_bases, indel_events, clip_bases}, the keys "unmapped" / "unplaced" (GAM), "chr1_22" (the sum over
chr1-chr22) and "other" (GAM: total records, groups, secondary; BAM: unmapped placed, secondary, supplementary,
duplicates, records without NM).
Finished segments / contigs are kept in <out>.parts/ and reused by a rerun with the same --procs.
Run: python full_scan.py --kind gam|bam --path P [--refpath R] --out PREFIX [--procs 16]
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, '/scratch/jshen/Github/Pansoma')
from indexed_gam_pipeline_v4 import gam_reader  # noqa: E402
from indexed_gam_pipeline_v4.gam_reader import BgzfReader, IndexedGam  # noqa: E402

BIN = Path('/scratch/jshen/data/reviewer_experiments_20261013/E11a_alignment_stats/full/bin')
FIELDS = ['records', 'perfect', 'sub', 'indel', 'clip', 'read_bases', 'aligned_bases', 'sub_bases', 'ins_bases',
          'del_bases', 'indel_events', 'clip_bases']
AUTOSOMES = [f'chr{i}' for i in range(1, 23)]


def run(job):
    """(stdout, seconds) of one C call; its stdout is kept in `part` so that a resubmitted scan reuses it."""
    cmd, part = job
    if part.exists():
        return part.read_text(), 0.0
    t0 = time.time()
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode:
        raise RuntimeError(f'{cmd}: exit {out.returncode}: {out.stderr.strip()}')
    tmp = part.with_suffix('.tmp')
    tmp.write_text(out.stdout)
    tmp.rename(part)
    return out.stdout, time.time() - t0


def add(acc, key, vals):
    row = acc.setdefault(key, dict.fromkeys(FIELDS, 0))
    for k, v in zip(FIELDS, vals):
        row[k] += v


def parts_dir(args):
    """<out>.parts/: one file per finished segment (GAM) or contig (BAM); the same --procs gives the same segments."""
    d = Path(args.out + '.parts')
    d.mkdir(exist_ok=True)
    return d


def scan_gam(args):
    contigs = [c['name'] for c in json.load(open(Path(args.refpath) / 'meta.json'))['contigs']]
    bins = IndexedGam(args.path, args.path + '.gai', cache_bytes=1).bins
    cuts = gam_reader.segments(Path(args.path).stat().st_size, bins, args.procs * args.segments_per_proc)
    jobs, parts = [], parts_dir(args)
    with BgzfReader(args.path) as r:
        for s, e in cuts:
            p = None if e is None else r.position(e)
            end = -1 if e is None else (-1 if p is None else (p[0] << 16) | p[1])
            jobs.append(([str(BIN / 'gam_fullstats'), args.path, str(Path(args.refpath) / 'chrom.npy'), str(len(contigs)),
                          str(s), str(end)], parts / f'{s}_{end}.tsv'))
    del bins, cuts  # the GAI bins of an Illumina GAM take several GB
    acc, other, secs = {}, dict(records=0, groups=0, secondary=0), []
    with ThreadPoolExecutor(args.procs) as pool:  # threads only wait on the C processes; no forked copies of the bins
        for out, sec in pool.map(run, jobs):
            secs.append(round(sec, 1))
            for line in out.splitlines():
                f = line.split('\t')
                if f[0] == '#total':
                    for k, v in zip(('records', 'groups', 'secondary'), f[1:]):
                        other[k] += int(v)
                    continue
                i = int(f[0])
                add(acc, 'unmapped' if i == -2 else 'unplaced' if i == -1 else contigs[i], map(int, f[1:]))
    return acc, other, dict(segments=len(jobs), segment_seconds=secs)


def scan_bam(args):
    parts = parts_dir(args)
    jobs = [([str(BIN / 'bam_fullstats'), args.path, c], parts / f'{c}.tsv') for c in AUTOSOMES]
    acc, other, secs = {}, dict(unmapped_placed=0, secondary=0, supplementary=0, duplicates=0, no_nm=0), []
    with ThreadPoolExecutor(args.procs) as pool:
        for out, sec in pool.map(run, jobs):
            secs.append(round(sec, 1))
            for line in out.splitlines():
                f = line.split('\t')
                if f[0] == '#other':
                    for k, v in zip(other, f[1:]):
                        other[k] += int(v)
                else:
                    add(acc, f[0], map(int, f[1:]))
    return acc, other, dict(contig_seconds=secs)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--kind', choices=['gam', 'bam'], required=True)
    p.add_argument('--path', required=True)
    p.add_argument('--refpath')
    p.add_argument('--out', required=True)
    p.add_argument('--procs', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', 4)))
    p.add_argument('--segments-per-proc', type=int, default=16)
    args = p.parse_args()
    t0 = time.time()
    acc, other, run_info = (scan_gam if args.kind == 'gam' else scan_bam)(args)
    tot = dict.fromkeys(FIELDS, 0)
    for c in AUTOSOMES:
        for k in FIELDS:
            tot[k] += acc.get(c, {}).get(k, 0)
    res = dict(kind=args.kind, path=args.path, refpath=args.refpath, per_contig=acc, chr1_22=tot, other=other,
               run=run_info, wall_seconds=round(time.time() - t0, 1))
    Path(args.out + '.json').write_text(json.dumps(res, indent=1) + '\n')
    print(json.dumps(dict(chr1_22=tot, other=other, wall_seconds=res['wall_seconds'])))


if __name__ == '__main__':
    main()
