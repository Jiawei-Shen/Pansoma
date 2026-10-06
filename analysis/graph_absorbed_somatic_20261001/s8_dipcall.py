"""Step 8: HG008-N dipcall (GRCh38 germline truth, PASS) around every locus of loci.tsv, SNVs and INDELs
(generalizes indel_graph_paths_20260930/i1_germline.py, which covered the Illumina I1 INDELs only).

Core span [lo, hi) of the truth: SNV = the base; INDEL = the VCF span (padding base included) united with the tandem-repeat
stretches touching the event (loci.tsv ev_lo0 + 5 .. ev_hi0 - 5; i1 used the graph bubble span instead). PASS germline
records of the dipcall VCF (labels manifest truth.germline) inside the graph window, relation in this order:
  identical germline allele      - one germline allele applied alone to the window gives alt_hap (any representation)
  germline variant at the same site (other allele)
                                 - INDEL truth: a length-changing germline allele overlapping [lo - 1, hi + 1];
                                   SNV truth: a germline allele that changes the truth base (other base, MNP, deletion)
  germline INDEL within 10 bp / germline SNV within 10 bp
                                 - germline records overlapping [lo - 10, hi + 10]
  none within 10 bp
gt = GT of the deciding record(s) (dipcall GT is phased hap1|hap2 of the same v6.2 assembly, so ident_haps / site_haps say
which HG008-N haplotype carries the identical / same-site allele); in_dip_bed = [lo - 10, hi + 10] inside one dip.bed
interval (as i1), in_dip_bed_window = the whole graph window inside one interval.
Output: dipcall_relation.tsv (same rows and order as loci.tsv). The cross-check with i1_germline.tsv is printed.
"""
import bisect, csv, json, os
from collections import Counter, defaultdict
import pysam
D = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
MAN = '/scratch/jshen/data/pansoma_v2_tensors/HG008T_Illumina/tensors/INDEL/labels.manifest.json'
if not os.path.exists(MAN):                       # tensors were being rebuilt on 2026-10-01; same truth files
    MAN = '/scratch/jshen/data/pansoma_v2_tensors/backup_ch6_linear100_20261001/HG008T_Illumina/tensors/INDEL/labels.manifest.json'
man = json.load(open(MAN))['truth']['germline']
L = list(csv.DictReader(open(f'{D}/loci.tsv'), delimiter='\t'))
IV = defaultdict(list)
for r in L:
    IV[r['chrom']].append((int(r['win_start0']) - 200, int(r['win_end0']) + 200))
for c, v in IV.items():                           # merge
    v.sort(); m = [list(v[0])]
    for a, b in v[1:]:
        if a <= m[-1][1]: m[-1][1] = max(m[-1][1], b)
        else: m.append([a, b])
    IV[c] = m
IVS = {c: [a for a, b in v] for c, v in IV.items()}
GERM = defaultdict(list)                          # chrom -> [(pos0, ref, alts, gt tuple)], one pass (VCF not indexed)
for rec in pysam.VariantFile(man['vcf']):
    if rec.chrom not in IV or (rec.filter.keys() and 'PASS' not in rec.filter.keys()):
        continue
    i = bisect.bisect_right(IVS[rec.chrom], rec.pos - 1) - 1
    if i >= 0 and IV[rec.chrom][i][1] >= rec.pos - 1:
        GERM[rec.chrom].append((rec.pos - 1, rec.ref.upper(), tuple(a.upper() for a in rec.alts or ()), rec.samples[0]['GT']))
GPOS = {c: [g[0] for g in v] for c, v in GERM.items()}
DIP = defaultdict(list)
for line in open(man['bed']):
    c, a, b = line.split()[:3]
    DIP[c].append((int(a), int(b)))
DIP = {c: sorted(v) for c, v in DIP.items()}
DST = {c: [a for a, b in v] for c, v in DIP.items()}


def in_dip(chrom, a, b):
    i = bisect.bisect_right(DST.get(chrom, []), a) - 1
    return i >= 0 and DIP[chrom][i][1] >= b


def apply(region, base, p0, a, b):
    off = p0 - base
    return None if off < 0 or off + len(a) > len(region) or region[off:off + len(a)] != a else region[:off] + b + region[off + len(a):]


def changed(p0, a, b):
    """GRCh38 bases a germline allele changes (same length: mismatching bases; INDEL: everything after a shared first base)."""
    if len(a) == len(b):
        return {p0 + i for i in range(len(a)) if a[i] != b[i]}
    k = 1 if a[:1] == b[:1] else 0
    return set(range(p0 + k, p0 + len(a)))


def haps(gt, idx):
    return '+'.join(h for h, g in zip(('hap1', 'hap2'), gt) if g in idx)


cols = ['truth_id', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind2', 'relation', 'records', 'gt', 'ident_haps', 'site_haps',
        'in_dip_bed', 'in_dip_bed_window', 'n_near']
out = []
for r in L:
    chrom, pos0, vref, valt = r['chrom'], int(r['vcf_pos']) - 1, r['vcf_ref'].upper(), r['vcf_alt'].upper()
    ws, we = int(r['win_start0']), int(r['win_end0'])
    region, alt_hap = r['ref_hap'], r['alt_hap']
    snv = r['kind2'] == 'SNV'
    lo, hi = (pos0, pos0 + 1) if snv else (min(pos0, int(r['ev_lo0']) + 5), max(pos0 + len(vref), int(r['ev_hi0']) - 5))
    G, P = GERM.get(chrom, []), GPOS.get(chrom, [])
    win = [g for g in G[bisect.bisect_left(P, ws - 300):bisect.bisect_right(P, we + 1)] if g[0] >= ws and g[0] + len(g[1]) <= we + 1]
    near = [g for g in win if g[0] + len(g[1]) > lo - 10 and g[0] < hi + 10]
    ident = [(g, i + 1) for g in win for i, a in enumerate(g[2]) if apply(region, ws, g[0], g[1], a) == alt_hap]
    if snv:
        same = [(g, i + 1) for g in near for i, a in enumerate(g[2]) if pos0 in changed(g[0], g[1], a)]
    else:
        same = [(g, i + 1) for g in near for i, a in enumerate(g[2]) if len(a) != len(g[1]) and g[0] <= hi + 1 and g[0] + len(g[1]) >= lo - 1]
    fmt = lambda gs: ';'.join(f'{g[0] + 1}:{g[1][:12]}>{",".join(a[:12] for a in g[2])}' for g in gs)
    gts = lambda gs: ';'.join('|'.join('.' if x is None else str(x) for x in g[3]) for g in gs)
    if ident:
        rel, dec = 'identical germline allele', [g for g, _ in ident]
    elif same:
        rel, dec = 'germline variant at the same site (other allele)', list(dict.fromkeys(g for g, _ in same))
    elif any(len(a) != len(g[1]) for g in near for a in g[2]):
        rel, dec = 'germline INDEL within 10 bp', near
    elif near:
        rel, dec = 'germline SNV within 10 bp', near
    else:
        rel, dec = 'none within 10 bp', []
    out.append(dict(r, relation=rel, records=fmt(list(dict.fromkeys(near + [g for g, _ in ident]))), gt=gts(dec),
                    ident_haps=';'.join(haps(g[3], {i}) for g, i in ident), site_haps=';'.join(haps(g[3], {i}) for g, i in same),
                    in_dip_bed=in_dip(chrom, lo - 10, hi + 10), in_dip_bed_window=in_dip(chrom, ws, we + 1), n_near=len(near)))
with open(f'{D}/dipcall_relation.tsv', 'w') as w:
    w.write('\t'.join(cols) + '\n')
    for o in out:
        w.write('\t'.join(str(o[c]) for c in cols) + '\n')
for kd in ('SNV', 'INDEL'):
    s = [o for o in out if o['kind2'] == kd]
    print(kd, len(s), dict(Counter(o['relation'] for o in s).most_common()), '| in_dip_bed', sum(o['in_dip_bed'] for o in s))
I1 = {(x['chrom'], x['vcf_pos'], x['vcf_ref'], x['vcf_alt']): x['germline_relation'].replace('INDEL at the same site', 'variant at the same site')
      for x in csv.DictReader(open('/scratch/jshen/data/pansoma_net_v2_runs/indel_graph_paths_20260930/i1_germline.tsv'), delimiter='\t')}
cmp = Counter((I1[k], o['relation']) for o in out if (k := (o['chrom'], o['vcf_pos'], o['vcf_ref'], o['vcf_alt'])) in I1)
print('vs i1_germline.tsv (overlapping truths', sum(cmp.values()), ', agree', sum(v for (a, b), v in cmp.items() if a == b), '):')
for (a, b), v in cmp.most_common():
    if a != b: print(f'  {v:4d}  i1: {a}  ->  now: {b}')
