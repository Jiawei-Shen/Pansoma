"""How many read alignments in a GAM lie on graph branch nodes, and how many of them are perfect (no edit)?

Input: the set's discovery/node_stats.json (indexed_gam_pipeline_v4 discovery: one full scan of the sorted GAM; per node,
`perfect` / `not_perfect` = MAPQ > 5 read visits of that node without / with an edit, after the builder's decoding and
indel left-normalization; an edit = any alignment column that is not a match: mismatch, insertion, deletion) and the
pipeline's GRCh38 reference-path index of the graph (graph_index/<graph>.grch38_path: chrom / visits per node).
GRCh38 node = on the GRCh38 path, visited once (the classify.py ref() rule); branch node = every other node (an allele
of the other haplotypes of the graph: HPRC samples, CHM13). Counts are per node visit ("mapping": one read on one node;
a read crossing 10 nodes counts 10 times).
Output: <set>.json in this folder (totals; GRCh38 vs branch: visits, perfect visits, nodes observed, nodes whose visits
are all perfect / <= 5 % edited / discovery targets; the same per node-length bin, because a visit's chance of an edit
grows with the node's length and branch nodes are mostly 1 bp) and, with `summary`, branch_node_alignments.md.
Run: python s1_branch_node_stats.py SET (Slurm, ~10 GB for a 5 GB node_stats.json); python s1_branch_node_stats.py summary
"""
import json, sys
from pathlib import Path
import numpy as np
T = Path('/scratch/jshen/data/pansoma_v2_tensors')
GP = T / 'graph_index/hprc-v1.1-mc-grch38.d9.grch38_path'
OUT = Path(__file__).resolve().parent
SETS = ['HG008T_PacBio', 'HG008T_ONT', 'HG008T_Illumina', 'COLO829T_fiberseq', 'COLO829T_ONT', 'COLO829T_Illumina']
BINS = [('1 bp', 1, 1), ('2-5 bp', 2, 5), ('6-32 bp', 6, 32), ('33-256 bp', 33, 256), ('>256 bp', 257, 10 ** 9)]
DIGITS = bytes(b if 48 <= b <= 57 else 32 for b in range(256))      # every non-digit byte -> space


def load(path):
    """node, perfect, not_perfect arrays of node_stats.json (records of 4 numbers: node, perfect, not_perfect, max_len)."""
    a = np.fromstring(path.read_bytes().translate(DIGITS), dtype=np.int64, sep=' ')
    assert a.size % 4 == 0
    a = a.reshape(-1, 4)
    return a[:, 0], a[:, 1], a[:, 2]


def stats(s):
    rep = json.load(open(T / s / 'discovery/discovery_report.json'))
    node, p, q = load(T / s / 'discovery/node_stats.json')
    assert node.size == rep['nodes_observed'], (node.size, rep['nodes_observed'])
    chrom, visits = np.load(GP / 'chrom.npy', mmap_mode='r'), np.load(GP / 'visits.npy', mmap_mode='r')
    length = np.asarray(np.load(GP / 'lengths.npy', mmap_mode='r')[node])
    ref = (np.asarray(chrom[node]) >= 0) & (np.asarray(visits[node]) == 1)
    n = p + q
    out = dict(set=s, gam=rep['gam'], alignments_scanned=rep['alignments_scanned'], alignments_mapq_gt5=rep['alignments_passing_mapq'],
               graph_nodes=int(len(visits) - 1), graph_grch38_nodes=int(((np.asarray(chrom) >= 0) & (np.asarray(visits) == 1)).sum()))
    for name, m in (('all', np.ones_like(ref)), ('grch38', ref), ('branch', ~ref)):
        nn, pp, qq = n[m], p[m], q[m]
        out[name] = dict(visits=int(nn.sum()), perfect_visits=int(pp.sum()), edited_visits=int(qq.sum()),
                         perfect_fraction=float(pp.sum() / max(1, nn.sum())), nodes_observed=int(m.sum()),
                         nodes_all_perfect=int((qq == 0).sum()), nodes_edited_le_5pct=int((qq <= 0.05 * nn).sum()),
                         nodes_target=int(((qq >= 1) & (qq > 0.05 * nn)).sum()),
                         visits_on_all_perfect_nodes=int(nn[qq == 0].sum()),
                         median_visits_per_node=float(np.median(nn)) if nn.size else 0.0)
    out['branch_share_of_visits'] = out['branch']['visits'] / out['all']['visits']
    out['by_length'] = {}
    for lab, lo, hi in BINS:
        for name, m in (('grch38', ref), ('branch', ~ref)):
            k = m & (length >= lo) & (length <= hi)
            out['by_length'][f'{name} {lab}'] = dict(nodes=int(k.sum()), visits=int(n[k].sum()), perfect_fraction=float(p[k].sum() / max(1, n[k].sum())))
    json.dump(out, open(OUT / f'{s}.json', 'w'), indent=1)
    print(json.dumps(out, indent=1))


def summary():
    R = [json.load(open(OUT / f'{s}.json')) for s in SETS if (OUT / f'{s}.json').exists()]
    f = lambda x: f'{x:,}'
    pc = lambda a, b: f'{100 * a / b:.1f}%'
    L = ['# Read alignments on HPRC v1.1 d9 branch nodes (from discovery node_stats.json)', '',
         'Visit = one MAPQ > 5 read on one node; perfect = no edit on that node (after decoding + indel left-normalization). '
         'Branch node = not on the GRCh38 path.', '',
         '| set | MAPQ>5 alignments | node visits | **on branch nodes** | share | **perfect on branch** | % perfect (branch) | % perfect (GRCh38 nodes) |',
         '|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in R:
        b, g = r['branch'], r['grch38']
        L.append(f"| {r['set']} | {f(r['alignments_mapq_gt5'])} | {f(r['all']['visits'])} | {f(b['visits'])} | {pc(b['visits'], r['all']['visits'])} | "
                 f"{f(b['perfect_visits'])} | {pc(b['perfect_visits'], b['visits'])} | {pc(g['perfect_visits'], g['visits'])} |")
    L += ['', '| set | branch nodes in the graph | branch nodes with reads | all visits perfect | <= 5 % edited (not a target) | discovery targets | '
          'GRCh38 nodes with reads | GRCh38 all perfect | GRCh38 <= 5 % edited |', '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in R:
        b, g = r['branch'], r['grch38']
        nb = r['graph_nodes'] - r['graph_grch38_nodes']
        L.append(f"| {r['set']} | {f(nb)} | {f(b['nodes_observed'])} ({pc(b['nodes_observed'], nb)}) | {f(b['nodes_all_perfect'])} ({pc(b['nodes_all_perfect'], b['nodes_observed'])}) | "
                 f"{f(b['nodes_edited_le_5pct'])} ({pc(b['nodes_edited_le_5pct'], b['nodes_observed'])}) | {f(b['nodes_target'])} | "
                 f"{f(g['nodes_observed'])} | {pc(g['nodes_all_perfect'], g['nodes_observed'])} | {pc(g['nodes_edited_le_5pct'], g['nodes_observed'])} |")
    L += ['', '### Perfect fraction of visits by node length (branch vs GRCh38 nodes; visits in parentheses)', '',
          '| set | ' + ' | '.join(f'{b[0]} branch | {b[0]} GRCh38' for b in BINS) + ' |', '|---|' + '---:|' * (2 * len(BINS))]
    for r in R:
        bl = r['by_length']
        L.append(f"| {r['set']} | " + ' | '.join(f"{100 * bl[f'{nm} {b[0]}']['perfect_fraction']:.1f}% ({bl[f'{nm} {b[0]}']['visits'] / 1e6:.0f} M)" for b in BINS for nm in ('branch', 'grch38')) + ' |')
    open(OUT / 'branch_node_alignments.md', 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    summary() if sys.argv[1] == 'summary' else stats(sys.argv[1])
