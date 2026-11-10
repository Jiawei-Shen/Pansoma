"""s4a: loci and k-mer pairs for the read check of one tumor.

usage: s4_loci.py TUMOR

Rows of s3 rows.tsv.gz:
  - every control row (group ctrl_*),
  - off-reference rows with asm_class tumor_only,
  - off-reference rows with asm_class neither and normal_side ref_exact / site_only that pass the gate for the s5
    reads-only tier: n_a1 >= 4, n_split 0, not a homopolymer / STR INDEL (only Illumina evidence counts there, so it
    cannot reach two platforms), and the key tumor_only / neither in >= 2 platform sets or AF >= 0.2; the others get the
    locus 'none:gate',
  - off-reference rows whose key is a known off-reference somatic key (KNOWN below), whatever their class,
  - off-reference rows with asm_class normal and at most two strict-hit loci in the normal (n_alt_loci <= 2: one per
    haplotype), for the germline read check of the s6 training labels (purpose 'germ'; the rows above: 'cand').
Rows without a usable query (asm_class no_query / short_flank) get the locus 'none:<asm_class>'.
Anchor (GRCh38 interval of the site): the s3 anchor_mapped (site of ref_q / alt_q on its primary GRCh38 alignment,
MAPQ >= 20) when there is no label anchor or it lies within the label anchor +- 1 kb ('mapped'); else the label anchor
('label', or 'label_conflict' when the mapped anchor is elsewhere). No anchor: 'none:no_anchor'; wider than 20 kb:
'none:wide'.
k-mer pair. ALT = alt_q[wl - fl : wh + fr] around the event window [wl, wh) of alt_q. REF = the same around the REF
window of the patient's normal sequence (s3 norm_q / norm_win; ref_src 'normal') for normal_side site_only /
site_plus_other, else of ref_q (graph REF; ref_src 'graph'). fl and fr grow separately (0..40, shortest pair first,
both k-mers >= 21 bp) until (1) the ALT k-mer (either strand) is not in the REF sequence and the REF k-mer not in alt_q,
and (2) neither k-mer occurs in GRCh38 within anchor +- 150 bp except over the anchor itself (s4_count only searches
the part of a read aligned near the anchor). specific = 1 when both hold; else 0 with note grch38_copy (only (1)
holds) or not_unique (neither).
Writes DATA/s4/<TUMOR>/loci.tsv (locus_id, chrom, start, end = anchor +- 200 bp fetch region, a_start, a_end, alt_k,
ref_k, ref_src, fl, fr, specific, note, cand, germ: 1 when a row of that purpose uses the locus), row_loci.tsv.gz (set,
kind, candidate_id, group, known, anchor_src, purpose, locus) and regions.bed (merged fetch regions of every locus).
Locus ids are content hashes (chrom, anchor, k-mers), so counts of an earlier run stay valid for the same loci.
"""
import csv
import gzip
import hashlib
import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, GRCH38, context, platform, rc, windows  # noqa: E402

R = '/scratch/jshen/data/reviewer_experiments_20261013'
KNOWN = {
    'HG008T': [('asm_truth', f'{R}/E1.1_ch6_counterfactual/asm_d46/offref_somatic_d46.tsv')],
    'COLO829T': [('smaht_truth', f'{R}/E9j_offref_colo829t/offref_somatic_d46.tsv'),
                 ('dsa_somatic', f'{R}/E9k_dsa_somatic_offref/all/offref_somatic_d46.tsv')],
}
PAD = 200
MAX_ANCHOR = 20_000
AGREE = 1000
REACH = 150
MIN_K = 21
MAX_F = 40
PAIRS = sorted(((fl, fr) for fl in range(MAX_F + 1) for fr in range(MAX_F + 1)), key=lambda p: (p[0] + p[1], abs(p[0] - p[1])))


def known_keys(tumor):
    out = {}
    for tag, path in KNOWN[tumor]:
        for r in csv.DictReader(open(path), delimiter='\t'):
            out.setdefault(r['key'], tag)
    return out


def iv(s):
    chrom, x = s.rsplit(':', 1)
    a, b = map(int, x.split('-'))
    return chrom, a, b


def pick_anchor(r):
    lab = iv(r['anchor']) if r['anchor'] else None
    mp = iv(r['anchor_mapped']) if r.get('anchor_mapped') else None
    if mp and (lab is None or (mp[0] == lab[0] and mp[1] < lab[2] + AGREE and mp[2] > lab[1] - AGREE)):
        return mp, 'mapped'
    if lab:
        return lab, 'label_conflict' if mp else 'label'
    return None, ''


def occurrences(hay, k):
    out, i = [], hay.find(k)
    while i >= 0:
        out.append(i)
        i = hay.find(k, i + 1)
    return out


def kmers(r, fasta, chrom, s, e):
    alt = r['alt_q']
    (wl, wh), (rl, rh) = windows(r)
    if r['normal_side'] in ('site_only', 'site_plus_other') and r['norm_q']:
        ref, src = r['norm_q'], 'normal'
        rl, rh = map(int, r['norm_win'].split('-'))
    else:
        ref, src = r['ref_q'], 'graph'
    g0 = max(0, s - REACH)
    g = fasta.fetch(chrom, g0, e + REACH).upper()
    first_local = last = None
    for fl, fr in PAIRS:
        if wl - fl < 0 or rl - fl < 0 or wh + fr > len(alt) or rh + fr > len(ref):
            continue
        ak, rk = alt[wl - fl:wh + fr], ref[rl - fl:rh + fr]
        if min(len(ak), len(rk)) < MIN_K:
            continue
        last = (ak, rk, fl, fr)
        if ak in ref or rc(ak) in ref or rk in alt or rc(rk) in alt:
            continue
        if first_local is None:
            first_local = last
        conflict = False
        for k in (ak, rc(ak), rk, rc(rk)):
            for x in occurrences(g, k):
                # an occurrence over the anchor itself is the GRCh38 allele at the site: a valid measurement
                if not (g0 + x < e + len(k) and g0 + x + len(k) > s - len(k)):
                    conflict = True
                    break
            if conflict:
                break
        if not conflict:
            return (*last, src, 1, '')
    if first_local:
        return (*first_local, src, 0, 'grch38_copy')
    if last:
        return (*last, src, 0, 'not_unique')
    return None


def main():
    tumor = sys.argv[1]
    known = known_keys(tumor)
    fasta = pysam.FastaFile(GRCH38)
    out = f'{DATA}/s4/{tumor}'
    os.makedirs(out, exist_ok=True)
    rows_gz = f'{DATA}/s3/{tumor}/rows.tsv.gz'
    keyplat = {}
    for r in csv.DictReader(gzip.open(rows_gz, 'rt'), delimiter='\t'):
        if r['group'] == 'offref' and r['asm_class'] in ('tumor_only', 'neither'):
            keyplat.setdefault(r['candidate_id'], set()).add(platform(r['set']))
    loci, flags = {}, {}
    n = {}
    with gzip.open(f'{out}/row_loci.tsv.gz.tmp', 'wt') as w:
        w.write('set\tkind\tcandidate_id\tgroup\tknown\tanchor_src\tpurpose\tlocus\n')
        for r in csv.DictReader(gzip.open(rows_gz, 'rt'), delimiter='\t'):
            kn = known.get(r['candidate_id'], '') if r['group'] == 'offref' else ''
            sel = (r['group'] != 'offref' or kn or r['asm_class'] == 'tumor_only' or
                   (r['asm_class'] == 'neither' and r['normal_side'] in ('ref_exact', 'site_only')))
            germ = (r['group'] == 'offref' and r['asm_class'] == 'normal' and int(r['n_alt_loci'] or 0) <= 2)
            if not (sel or germ):
                continue
            purpose = 'cand' if sel else 'germ'
            src = ''
            gate = (r['group'] != 'offref' or kn or r['asm_class'] != 'neither' or
                    (int(r['n_a1']) >= 4 and r['n_split'] == '0' and
                     not (r['kind'] == 'INDEL' and context(r) in ('hp', 'str')) and
                     (len(keyplat.get(r['candidate_id'], ())) >= 2 or float(r['af']) >= 0.2)))
            if r['asm_class'] in ('no_query', 'short_flank'):
                lid = f'none:{r["asm_class"]}'
            elif not gate:
                lid = 'none:gate'
            else:
                anc, src = pick_anchor(r)
                if anc is None:
                    lid = 'none:no_anchor'
                elif anc[2] - anc[1] > MAX_ANCHOR:
                    lid = 'none:wide'
                elif anc[0] not in fasta.references:
                    lid = 'none:contig'
                else:
                    chrom, s, e = anc
                    km = kmers(r, fasta, chrom, s, e)
                    if km is None:
                        lid = 'none:no_kmer'
                    else:
                        ak, rk, fl, fr, rsrc, spec, note = km
                        rec = (chrom, max(0, s - PAD), e + PAD, s, e, ak, rk, rsrc, fl, fr, spec, note)
                        lid = 'L' + hashlib.sha1(f'{chrom}|{s}|{e}|{ak}|{rk}'.encode()).hexdigest()[:12]
                        loci[lid] = rec
                        flags.setdefault(lid, set()).add(purpose)
            key = (purpose, lid if lid.startswith('none') else 'locus')
            n[key] = n.get(key, 0) + 1
            w.write(f'{r["set"]}\t{r["kind"]}\t{r["candidate_id"]}\t{r["group"]}\t{kn}\t{src}\t{purpose}\t{lid}\n')
    order = sorted(loci.items(), key=lambda kv: (kv[1][0], kv[1][1], kv[1][2], kv[0]))
    with open(f'{out}/loci.tsv.tmp', 'w') as w:
        w.write('locus_id\tchrom\tstart\tend\ta_start\ta_end\talt_k\tref_k\tref_src\tfl\tfr\tspecific\tnote\tcand\tgerm\n')
        for lid, rec in order:
            fl_ = flags[lid]
            w.write(lid + '\t' + '\t'.join(map(str, rec)) + f'\t{int("cand" in fl_)}\t{int("germ" in fl_)}\n')
    with open(f'{out}/regions.bed.tmp', 'w') as w:
        cur = None
        for _, (c, s, e, *_r) in order:
            if cur and cur[0] == c and s <= cur[2] + 100:
                cur[2] = max(cur[2], e)
            else:
                if cur:
                    w.write(f'{cur[0]}\t{cur[1]}\t{cur[2]}\n')
                cur = [c, s, e]
        if cur:
            w.write(f'{cur[0]}\t{cur[1]}\t{cur[2]}\n')
    os.replace(f'{out}/row_loci.tsv.gz.tmp', f'{out}/row_loci.tsv.gz')
    os.replace(f'{out}/loci.tsv.tmp', f'{out}/loci.tsv')
    os.replace(f'{out}/regions.bed.tmp', f'{out}/regions.bed')
    spec = sum(1 for v in loci.values() if v[10] == 1)
    n_c = sum(1 for f in flags.values() if 'cand' in f)
    print(tumor, 'rows', dict(sorted(n.items())), 'unique loci', len(loci), 'cand', n_c, 'germ only', len(loci) - n_c,
          'specific', spec)


if __name__ == '__main__':
    main()
