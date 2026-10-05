"""COLO829T step 1: one d9 graph window per truth allele of variant_set.tsv (the shared contract of every later step).

Same rule as analysis/graph_absorbed_somatic_20261001/s1_loci.py (HG008T). Input: $D/variant_set.tsv (c0, 497 truths).
Window = the GRCh38 nodes overlapping [pos0 - 60, pos0 + len(REF) + 60], widened to the event window +-20 bp when that
is larger (long tandem repeats; classify.py window_nodes, d9 pipeline reference-path index), first to last of them:
  first_node / last_node   - anchor GRCh38 nodes (forward); every local path is read between them
  win_start0 / win_end0    - GRCh38 span of the window, 0-based, end inclusive (first node start .. last node end)
  ref_hap / alt_hap        - GRCh38 sequence of the window, and the same with the truth allele applied
  ev_lo0 / ev_hi0          - event window, 0-based half-open GRCh38: the union of every periodic stretch (period 1-6,
                             length >= max(2p, p + 3)) touching the event, +-5 bp (the a3_windows.py rule), so a length
                             change anywhere in the surrounding tandem repeat counts as the same event
Differences from HG008: no asm_event column (the COLO829T truth VCF has no patient-frame INFO such as
HG008Nv62SOMATICVARIANT; the COLO829BL frame comes from the donor assembly in a later step); VAF_Ill / VAF_PB / RGN /
n_platforms are carried along.
Output: $D/loci.tsv (same rows and order as variant_set.tsv; status no_window if no GRCh38 node / a reverse GRCh38 node,
event_window_exceeds_graph_window if the event window is not inside the graph window).
Assumptions: the d9 graph (hprc-v1.1-mc-grch38.d9) is the graph the COLO829T tensors were built on; REF is checked
against GRCh38 (assert).
"""
import csv
from collections import Counter
import pysam
D = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005'
C = '/scratch/jshen/data/pansoma_net_v2_runs/indel_graph_paths_20260930/classify.py'
ns = {}
exec(open(C).read().split("fa = pysam.FastaFile(FASTA)\nchroms")[0], ns)
window_nodes = ns['window_nodes']
fa = pysam.FastaFile(ns['FASTA'])
FLANK = 60
rows = list(csv.DictReader(open(f'{D}/variant_set.tsv'), delimiter='\t'))


def event_window(chrom, pos0, vref, valt):
    lref = len(vref)
    if lref == 1 and len(valt) == 1:
        e_lo, e_hi = pos0, pos0 + 1                                   # SNV: the base itself
    else:
        e_lo, e_hi = pos0 + 1, max(pos0 + lref, pos0 + 1)              # INDEL after the padding base (INS: a point)
    s0 = max(0, e_lo - 300)
    seq = fa.fetch(chrom, s0, e_hi + 300).upper()
    lo, hi = e_lo - s0, e_hi - s0
    L, R = lo, hi
    for p in range(1, 7):
        i = 0
        while i < len(seq) - p:
            if seq[i] != seq[i + p]:
                i += 1; continue
            j = i
            while j < len(seq) - p and seq[j] == seq[j + p]:
                j += 1
            a, b = i, j + p
            if b - a >= max(2 * p, p + 3) and a <= hi + 1 and b >= lo - 1:
                L, R = min(L, a), max(R, b)
            i = j + 1
    return s0 + L - 5, s0 + R + 5


cols = ['truth_id', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind2', 'status', 'first_node', 'last_node', 'n_ref_nodes',
        'win_start0', 'win_end0', 'ev_lo0', 'ev_hi0', 'VAF_Ill', 'VAF_PB', 'RGN', 'n_platforms', 'ref_hap', 'alt_hap']
st = Counter()
with open(f'{D}/loci.tsv', 'w') as w:
    w.write('\t'.join(cols) + '\n')
    for r in rows:
        chrom, pos1, vref, valt = r['chrom'], int(r['vcf_pos']), r['vcf_ref'].upper(), r['vcf_alt'].upper()
        pos0 = pos1 - 1
        ev = event_window(chrom, pos0, vref, valt)
        nodes = window_nodes(chrom, min(pos0 - FLANK, ev[0] - 20), max(pos0 + len(vref) + FLANK, ev[1] + 20))
        out = dict(r, status='ok', ev_lo0=ev[0], ev_hi0=ev[1])
        if not nodes:
            out['status'] = 'no_window'
        else:
            first, last = nodes[0], nodes[-1]
            region = fa.fetch(chrom, first[1], last[2] + 1).upper()
            off = pos0 - first[1]
            assert region[off:off + len(vref)] == vref, (chrom, pos1)
            out.update(first_node=first[0], last_node=last[0], n_ref_nodes=len(nodes), win_start0=first[1], win_end0=last[2],
                       ref_hap=region, alt_hap=region[:off] + valt + region[off + len(vref):])
            if not (first[1] <= ev[0] and ev[1] <= last[2] + 1):
                out['status'] = 'event_window_exceeds_graph_window'
        st[(r['kind2'], out['status'])] += 1
        w.write('\t'.join(str(out.get(c, '')) for c in cols) + '\n')
print(len(rows), 'loci', sorted(st.items()))
L = list(csv.DictReader(open(f'{D}/loci.tsv'), delimiter='\t'))
ok = [x for x in L if x['status'] == 'ok']
span = sorted(int(x['win_end0']) - int(x['win_start0']) + 1 for x in ok)
ev = sorted(int(x['ev_hi0']) - int(x['ev_lo0']) for x in ok)
print('window bp min/median/max', span[0], span[len(span) // 2], span[-1], '| event window bp', ev[0], ev[len(ev) // 2], ev[-1],
      '| GRCh38 nodes median', sorted(int(x['n_ref_nodes']) for x in ok)[len(ok) // 2])

# Result (2026-10-05, login node, 3.4 min, 0.18 GB):
#   497 loci, all status ok (INDEL 430, SNV 67); no no_window / event_window_exceeds_graph_window
#   graph window 124 / 317 / 1,361 bp (min / median / max), event window 10 / 28 / 146 bp, median 7 GRCh38 nodes
