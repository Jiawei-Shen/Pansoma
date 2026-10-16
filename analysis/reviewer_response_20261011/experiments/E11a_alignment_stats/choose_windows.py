"""E11a window design: fixed GRCh38 windows (the same for every set and arm), and each window's node-ID range in each
graph (from the pipeline's reference-path arrays).

Windows (seed 20261013, chr1-22 and chrX; chrY left out):
  random  24 x 500 kb, uniform over the genome, >= 95 % covered by the SMaHT v2 union (no assembly gaps), >= 1 Mb apart
  segdup  24 x 100 kb, centred on a random base of GRCh38_segdups.bed, >= 50 % segdup, >= 1 Mb from every other window
Node range of a window in a graph: the GRCh38-path nodes whose path start lies in [start - flank, end + flank); lo / hi =
their 0.5 / 99.5 percentile node ID (a few far-away IDs would otherwise widen the range); reported: the ID span and
how many GRCh38 nodes fall outside [lo, hi].
Extra windows (long-read BAM arm only; seed 20261014): 100 x 500 kb random windows by the same rules, >= 1 Mb from
the core windows and from each other -> windows_extra.tsv.
Writes windows.tsv, windows_extra.tsv and window_nodes.json in this folder.
"""
import json
from pathlib import Path
import random

import numpy as np

HERE = Path(__file__).resolve().parent
SEED = 20261013
REFS = {'d9': '/scratch/jshen/data/pansoma_v2_tensors/graph_index/hprc-v1.1-mc-grch38.d9.grch38_path',
        'd46': '/scratch/jshen/data/pansoma_v2_tensors/graph_index_hprc-v2.1-d46/hprc-v2.1-mc-grch38.d46.grch38_path',
        'linear': '/scratch/jshen/data/Linear_Reference_VG_Indexes/graph_index/hg38_linear.giraffe.grch38_path'}
FLANK = 25_000  # the largest flank used by window_scan (long reads); short reads use a subrange


def main():
    import sys
    sys.path.insert(0, str(HERE))
    from compare import Bed, BEDS
    union = [Bed(BEDS[k]) for k in ('smaht_easy', 'smaht_difficult', 'smaht_extreme')]
    seg = Bed(BEDS['segdup'])
    contigs = [f'chr{i}' for i in range(1, 23)] + ['chrX']
    lengths = {}
    for c in contigs:
        lengths[c] = max(int(b.iv[c][1][-1]) for b in union if c in b.iv)
    rng = random.Random(SEED)
    total = sum(lengths.values())
    wins = []

    def far(c, s, e):
        return all(not (w[0] == c and s < w[2] + 1_000_000 and w[1] < e + 1_000_000) for w in wins)

    def cov(bed, c, s, e):
        return float(bed.fraction(np.array([c]), np.array([s]), np.array([e]))[0]) if c in bed.iv else 0.0

    while sum(1 for w in wins if w[3] == 'random') < 24:
        x = rng.randrange(total)
        for c in contigs:
            if x < lengths[c]:
                break
            x -= lengths[c]
        s = x - x % 1000
        e = s + 500_000
        if e > lengths[c] or not far(c, s, e):
            continue
        if sum(cov(b, c, s, e) for b in union) < 0.95:
            continue
        wins.append((c, s, e, 'random'))
    sc, ss, se = [], [], []
    for c in contigs:
        if c in seg.iv:
            ms, me, _ = seg.iv[c]
            sc += [c] * len(ms)
            ss += ms.tolist()
            se += me.tolist()
    seglen = np.array(se) - np.array(ss)
    cum = np.cumsum(seglen)
    tries = 0
    while sum(1 for w in wins if w[3] == 'segdup') < 24 and tries < 100000:
        tries += 1
        x = rng.randrange(int(cum[-1]))
        k = int(np.searchsorted(cum, x, side='right'))
        pos = ss[k] + (x - (cum[k - 1] if k else 0))
        c = sc[k]
        s = max(0, pos - 50_000)
        s -= s % 1000
        e = s + 100_000
        if not far(c, s, e) or cov(seg, c, s, e) < 0.5 or sum(cov(b, c, s, e) for b in union) < 0.95:
            continue
        wins.append((c, s, e, 'segdup'))
    wins.sort(key=lambda w: (contigs.index(w[0]), w[1]))
    core = list(wins)
    rng2 = random.Random(SEED + 1)
    extra = []
    while len(extra) < 100:
        x = rng2.randrange(total)
        for c in contigs:
            if x < lengths[c]:
                break
            x -= lengths[c]
        s = x - x % 1000
        e = s + 500_000
        if e > lengths[c] or not all(not (w[0] == c and s < w[2] + 1_000_000 and w[1] < e + 1_000_000)
                                     for w in core + extra):
            continue
        if sum(cov(b, c, s, e) for b in union) < 0.95:
            continue
        extra.append((c, s, e, 'random_extra'))
    extra.sort(key=lambda w: (contigs.index(w[0]), w[1]))
    with open(HERE / 'windows_extra.tsv', 'w') as f:
        f.write('window\tchrom\tstart\tend\tkind\tsegdup_frac\tsmaht_easy\tsmaht_difficult\tsmaht_extreme\n')
        for i, (c, s, e, k) in enumerate(extra):
            f.write(f'x{i:03d}\t{c}\t{s}\t{e}\t{k}\t{cov(seg, c, s, e):.3f}\t' +
                    '\t'.join(f'{cov(b, c, s, e):.3f}' for b in union) + '\n')
    with open(HERE / 'windows.tsv', 'w') as f:
        f.write('window\tchrom\tstart\tend\tkind\tsegdup_frac\tsmaht_easy\tsmaht_difficult\tsmaht_extreme\n')
        for i, (c, s, e, k) in enumerate(wins):
            f.write(f'w{i:02d}\t{c}\t{s}\t{e}\t{k}\t{cov(seg, c, s, e):.3f}\t' +
                    '\t'.join(f'{cov(b, c, s, e):.3f}' for b in union) + '\n')
    nodes = {}
    for arm, d in REFS.items():
        d = Path(d)
        meta = json.load(open(d / 'meta.json'))
        cidx = {c['name']: c for c in meta['contigs']}
        pn = np.load(d / 'path_nodes.npy', mmap_mode='r')
        ps = np.load(d / 'path_starts.npy', mmap_mode='r')
        nodes[arm] = {}
        for i, (c, s, e, k) in enumerate(wins):
            ci = cidx[c]
            sl = slice(ci['offset'], ci['offset'] + ci['nodes'])
            st = np.asarray(ps[sl])
            out = {}
            for flank in (1000, FLANK):
                a, b = np.searchsorted(st, s - flank), np.searchsorted(st, e + flank)
                ids = np.asarray(pn[sl][a:b])
                lo, hi = (int(np.percentile(ids, 0.5)), int(np.percentile(ids, 99.5))) if len(ids) else (0, -1)
                out[str(flank)] = dict(lo=lo, hi=hi, ref_nodes=int(len(ids)), id_span=hi - lo + 1,
                                       ref_outside=int(((ids < lo) | (ids > hi)).sum()))
            nodes[arm][f'w{i:02d}'] = out
    json.dump(nodes, open(HERE / 'window_nodes.json', 'w'), indent=1)
    for arm in nodes:
        spans = [v[str(FLANK)]['id_span'] / max(1, v[str(FLANK)]['ref_nodes']) for v in nodes[arm].values()]
        print(arm, 'id span / ref nodes: median', np.median(spans), 'max', max(spans))


if __name__ == '__main__':
    main()
