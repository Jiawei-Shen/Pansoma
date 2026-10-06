"""Step 1: one d9 graph window per truth allele of variant_set.tsv (the shared contract of every later step).

Window = the GRCh38 nodes overlapping [pos0 - 60, pos0 + len(REF) + 60], widened to the event window +-20 bp when that
is larger (long tandem repeats; classify.py window_nodes, d9 pipeline reference-path index), first to last of them:
  first_node / last_node   - anchor GRCh38 nodes (forward); every local path is read between them
  win_start0 / win_end0    - GRCh38 span of the window, 0-based, end inclusive (first node start .. last node end)
  ref_hap / alt_hap        - GRCh38 sequence of the window, and the same with the truth allele applied
  ev_lo0 / ev_hi0          - event window, 0-based half-open GRCh38: the union of every periodic stretch (period 1-6,
                             length >= max(2p, p + 3)) touching the event, +-5 bp (the a3_windows.py rule), so a length
                             change anywhere in the surrounding tandem repeat counts as the same event
  asm_event                - truth INFO HG008Nv62SOMATICVARIANT (normal-assembly haplotype:pos(0-based)-REF-ALT)
Output: loci.tsv (same rows and order as variant_set.tsv; status no_window if no GRCh38 node / reverse node).
"""
import csv, json
import pysam
D = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
C = '/scratch/jshen/data/pansoma_net_v2_runs/indel_graph_paths_20260930/classify.py'
ns = {}
exec(open(C).read().split("fa = pysam.FastaFile(FASTA)\nchroms")[0], ns)
window_nodes = ns['window_nodes']
fa = pysam.FastaFile(ns['FASTA'])
FLANK = 60
# the tensors' labels manifest (names the somatic truth VCF); since the 2026-10-01 tensor rebuild the 2026-09-30 sets
# (whose no-candidate statuses s0 uses) are in backup_ch6_linear100_20261001
man = json.load(open('/scratch/jshen/data/pansoma_v2_tensors/backup_ch6_linear100_20261001/HG008T_Illumina/tensors/INDEL/labels.manifest.json'))
rows = list(csv.DictReader(open(f'{D}/variant_set.tsv'), delimiter='\t'))
keys = {(r['chrom'], int(r['vcf_pos']), r['vcf_ref'], r['vcf_alt']) for r in rows}
asm = {}
for rec in pysam.VariantFile(man['truth']['somatic']['vcf']):
    for a in rec.alts or ():
        if (rec.chrom, rec.pos, rec.ref, a) in keys:
            v = rec.info.get('HG008Nv62SOMATICVARIANT'); asm[(rec.chrom, rec.pos, rec.ref, a)] = (v[0] if isinstance(v, tuple) else v) or ''


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
        'win_start0', 'win_end0', 'ev_lo0', 'ev_hi0', 'asm_event', 'ref_hap', 'alt_hap']
n_bad = 0
with open(f'{D}/loci.tsv', 'w') as w:
    w.write('\t'.join(cols) + '\n')
    for r in rows:
        chrom, pos1, vref, valt = r['chrom'], int(r['vcf_pos']), r['vcf_ref'].upper(), r['vcf_alt'].upper()
        pos0 = pos1 - 1
        ev = event_window(chrom, pos0, vref, valt)
        nodes = window_nodes(chrom, min(pos0 - FLANK, ev[0] - 20), max(pos0 + len(vref) + FLANK, ev[1] + 20))
        out = dict(r, status='ok', asm_event=asm.get((chrom, pos1, r['vcf_ref'], r['vcf_alt']), ''), ev_lo0=ev[0], ev_hi0=ev[1])
        if not nodes:
            out['status'] = 'no_window'; n_bad += 1
        else:
            first, last = nodes[0], nodes[-1]
            region = fa.fetch(chrom, first[1], last[2] + 1).upper()
            off = pos0 - first[1]
            assert region[off:off + len(vref)] == vref, (chrom, pos1)
            out.update(first_node=first[0], last_node=last[0], n_ref_nodes=len(nodes), win_start0=first[1], win_end0=last[2],
                       ref_hap=region, alt_hap=region[:off] + valt + region[off + len(vref):])
            if not (first[1] <= ev[0] and ev[1] <= last[2] + 1):
                out['status'] = 'event_window_exceeds_graph_window'
        w.write('\t'.join(str(out.get(c, '')) for c in cols) + '\n')
print(len(rows), 'loci', n_bad, 'without window')
