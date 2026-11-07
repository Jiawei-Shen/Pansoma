"""E5d step 5: does a normal haplotype carry the truth ALT together with other germline differences in the same array?

Why (2026-11-05): s8b (HG008T) and c3 (COLO829T) call a normal haplotype ALT only when its whole tandem array between the
common anchors equals GRCh38 + truth (s8b also allows the hap's dipcall records outside the array). A haplotype that has
the truth allele and another germline difference inside the array (an SNV next to the event, a second length change in
the STR) was called 'other' / 'germline_*', and the allele counted as absent from the normal. E20b found it: HG008-N
reads carry the ALT at 28.1% of the reads at the absorbed SNVs (COLO829BL 2.2%), concentrated at loci where the HG008-N
dipcall has the same ALT at the site.
Rule: per normal hap and sequence source (assembly S between the anchors; dipcall D = GRCh38 with the hap's PASS records
applied), 'ALT_plus_germline' when the sequence equals R with the truth allele plus the hap's own dipcall records applied,
leaving out a set Y of at most two of those records (the records that spell the truth allele in dipcall's representation;
Y may be empty for S). Records must lie inside [a, b) and not overlap each other or the truth (edit() refuses overlaps).
Only haps that are not already ALT are tested; the old call is kept otherwise.
Inputs: HG008T: loci.tsv + array_alleles.tsv of 2026-10-01 (OLD) and of E5d e2 ($DATA/hg008t); HG008-N dipcall (no index,
  one pass). COLO829T: c3 normal_frame.tsv (R_seq, hapX/hapY_seq, src, asm/dip calls); COLO829BL dipcall (tabix, GT
  hapY|hapX).
Output: $DATA/alt_with_germline.tsv (sample, truth_id, hap, src, asm_old, asm_new, dip_old, dip_new, n_recs, Y).
Run: python e5_alt_with_germline.py   (login node, ~3 min: the HG008-N dipcall pass)
"""
import collections
import csv
import gzip
import itertools
import sys

import pysam

DATA = '/scratch/jshen/data/reviewer_experiments_20261013/E5d_absorbed_normal_d46'
OLD = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
COLO = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005'
GR = '/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta'
HDIP = '/scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.vcf.gz'
CDIP = '/scratch/qfu/COLO829BL_DSA/dipcall_hg38/dipcall_hg38.dip.vcf.gz'
csv.field_size_limit(sys.maxsize)
rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))  # noqa: E731


def edit(R, a, eds):
    """s8b / c3 edit(): R (starting at a) with non-overlapping (p0, ref, alt) edits; None on overlap or REF mismatch."""
    out, cur = [], a
    for p0, r, x in sorted(eds):
        if p0 < cur or R[p0 - a:p0 - a + len(r)] != r:
            return None
        out += [R[cur - a:p0 - a], x]
        cur = p0 + len(r)
    return ''.join(out) + R[cur - a:]


def hap_recs(recs, h, a, b):
    """the hap's applied PASS records inside [a, b) (GT index h), or None when one is non-PASS, crosses a / b or lacks a GT."""
    out = []
    for p0, r, alts, gt, ok in recs:
        g = gt[h] if gt and len(gt) > h else None
        if not ok or g is None:
            return None
        if g == 0 or alts[g - 1] == '*':
            continue
        if p0 < a or p0 + len(r) > b:
            return None
        out.append((p0, r, alts[g - 1]))
    return out


def with_germline(seq, R, a, truth, hr):
    """Y (records left out) when seq == R + truth + (hr - Y), |Y| <= 2, else None."""
    if seq is None or hr is None:
        return None
    for k in range(3):
        for Y in itertools.combinations(hr, k):
            X = [e for e in hr if e not in Y]
            if edit(R, a, [truth] + X) == seq:
                return list(Y)
    return None


def new_call(old, seq, R, a, truth, hr):
    if old in ('ALT', 'ALT_plus_flank', 'NA') or seq is None:
        return old, None
    Y = with_germline(seq, R, a, truth, hr)
    return ('ALT_plus_germline', Y) if Y is not None else (old, None)


def read_vcf_stream(path, want):
    """{truth_id: [(p0, ref, alts, gt, pass)]} for records overlapping each [a, b) of want {chrom: [(a, b, id)]} (one pass)."""
    out = collections.defaultdict(list)
    iv = {c: sorted(v) for c, v in want.items()}
    span = max((b - a for v in iv.values() for a, b, _ in v), default=0) + 1
    import bisect
    st = {c: [x[0] for x in v] for c, v in iv.items()}
    with gzip.open(path, 'rt') as f:
        for ln in f:
            if ln[0] == '#':
                continue
            x = ln.split('\t', 10)
            v = iv.get(x[0])
            if not v:
                continue
            q0 = int(x[1]) - 1
            ref = x[3].upper()
            q1 = q0 + len(ref)
            for a, b, t in v[bisect.bisect_left(st[x[0]], q0 - span):bisect.bisect_left(st[x[0]], q1)]:
                if q0 < b and q1 > a:
                    gt = tuple(None if g == '.' else int(g) for g in x[9].split(':')[0].replace('/', '|').split('|'))
                    alts = tuple(z.upper() for z in x[4].split(','))
                    out[t].append((q0, ref, alts, gt, x[6] in ('PASS', '.')))
    return out


def hg008t(fa):
    rows, want, info = [], collections.defaultdict(list), {}
    for d in (OLD, f'{DATA}/hg008t'):
        loci = {r['truth_id']: r for r in rd(f'{d}/loci.tsv')}
        for ar in rd(f'{d}/array_alleles.tsv'):
            t = ar['truth_id']
            if ar['anchor_a0'] == '' or ar['anchor_b0'] == '' or t in info:
                continue
            L = loci[t]
            a, b = int(ar['anchor_a0']), int(ar['anchor_b0'])
            info[t] = (L, ar, a, b)
            want[L['chrom']].append((a, b, t))
    recs = read_vcf_stream(HDIP, want)
    chk = collections.Counter()
    for t, (L, ar, a, b) in info.items():
        R = fa.fetch(L['chrom'], a, b).upper()
        truth = (int(L['vcf_pos']) - 1, L['vcf_ref'].upper(), L['vcf_alt'].upper())
        for h, hi in (('hap1', 0), ('hap2', 1)):
            hr = hap_recs(recs[t], hi, a, b)
            D = edit(R, a, hr) if hr is not None and ar[f'{h}_dip_na'] == '' else None
            S = ar[f'{h}_seq'] if ar[f'{h}_asm_call'] != 'NA' else None
            # check: the recomputed dipcall sequence reproduces s8b's dipcall call
            if ar[f'{h}_dip_call'] in ('REF', 'ALT'):
                A = edit(R, a, [truth])
                chk[(ar[f'{h}_dip_call'], D == (R if ar[f'{h}_dip_call'] == 'REF' else A))] += 1
            an, ay = new_call(ar[f'{h}_asm_call'], S, R, a, truth, hr)
            dn, dy = new_call(ar[f'{h}_dip_call'], D, R, a, truth, hr)
            rows.append(dict(sample='HG008T', truth_id=t, hap=h, kind='SNV' if len(truth[1]) == len(truth[2]) == 1 else 'INDEL',
                             asm_old=ar[f'{h}_asm_call'], asm_new=an, dip_old=ar[f'{h}_dip_call'], dip_new=dn,
                             n_recs='' if hr is None else len(hr), Y=';'.join(f'{p + 1}:{r}>{x}' for p, r, x in (ay or dy or []))))
    print('HG008T dipcall recompute check (s8b call, reproduced):', dict(chk))
    return rows


def colo():
    rows, chk = [], collections.Counter()
    vcf = pysam.VariantFile(CDIP)
    for L in rd(f'{COLO}/normal_frame.tsv'):
        if L['anchor_a0'] == '' or L['anchor_b0'] == '' or not L['R_seq']:
            continue
        a, b, R = int(L['anchor_a0']), int(L['anchor_b0']), L['R_seq']
        truth = (int(L['vcf_pos']) - 1, L['vcf_ref'].upper(), L['vcf_alt'].upper())
        recs = [(r.start, r.ref.upper(), tuple((x or '').upper() for x in r.alts or ()), r.samples[0]['GT'],
                 not r.filter.keys() or 'PASS' in r.filter.keys()) for r in vcf.fetch(L['chrom'], a, b)]
        for h, hi in (('hapY', 0), ('hapX', 1)):
            hr = hap_recs(recs, hi, a, b)
            D = edit(R, a, hr) if hr is not None and L[f'{h}_dip_na'] == '' else None
            S = L[f'{h}_seq'] if L[f'{h}_src'] == 'asm' else None
            if L[f'{h}_dip_call'] in ('REF', 'ALT'):
                chk[(L[f'{h}_dip_call'], D == (R if L[f'{h}_dip_call'] == 'REF' else L['A_seq']))] += 1
            an, ay = new_call(L[f'{h}_asm_call'], S, R, a, truth, hr)
            dn, dy = new_call(L[f'{h}_dip_call'], D, R, a, truth, hr)
            rows.append(dict(sample='COLO829T', truth_id=L['truth_id'], hap=h,
                             kind='SNV' if len(truth[1]) == len(truth[2]) == 1 else 'INDEL',
                             asm_old=L[f'{h}_asm_call'], asm_new=an, dip_old=L[f'{h}_dip_call'], dip_new=dn,
                             n_recs='' if hr is None else len(hr), Y=';'.join(f'{p + 1}:{r}>{x}' for p, r, x in (ay or dy or []))))
    print('COLO829T dipcall recompute check (c3 call, reproduced):', dict(chk))
    return rows


def main():
    rows = hg008t(pysam.FastaFile(GR)) + colo()
    with open(f'{DATA}/alt_with_germline.tsv', 'w') as w:
        w.write('\t'.join(rows[0]) + '\n')
        for r in rows:
            w.write('\t'.join(str(r[k]) for k in rows[0]) + '\n')
    for s in ('HG008T', 'COLO829T'):
        x = [r for r in rows if r['sample'] == s]
        c = collections.Counter((r['kind'], r['asm_old'], r['asm_new'], r['dip_old'], r['dip_new']) for r in x
                                if 'ALT_plus_germline' in (r['asm_new'], r['dip_new']))
        print(s, 'haps', len(x), 'upgraded', sum(c.values()))
        for k, v in c.most_common():
            print('  ', k, v)


if __name__ == '__main__':
    main()
