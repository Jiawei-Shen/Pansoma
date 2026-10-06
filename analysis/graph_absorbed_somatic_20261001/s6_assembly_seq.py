"""Step 6: which allele does each HG008-N (normal) and HG008-T (tumor) assembly haplotype carry at each locus? (minimap2)

All loci of loci.tsv (SNVs and INDELs; the earlier Illumina-INDEL-only version is indel_graph_paths_20260930/assembly/
a1_queries.py + a2_parse.py + a3_windows.py, whose rules are kept).
  python s6_assembly_seq.py queries   -> asm_queries.fa / asm_queries.tsv
  (job s6_minimap2.sh, run from tmp/graph_absorbed_somatic_20261001/: minimap2 2.28 -ax sr --secondary=yes -N 20 --cs
   vs each assembly -> asm_normal.sam / asm_tumor.sam)
  python s6_assembly_seq.py parse     -> asm_seq_status.tsv / asm_seq_summary.tsv
Queries: '<truth_id>|ref' = GRCh38 [pos0 - F, pos0 + len(REF) + F), '<truth_id>|alt' = the same with the truth allele
applied; F = max(100, room for the event window (loci.tsv ev_lo0..ev_hi0 = tandem-repeat union touching the event +-5 bp,
the a3 rule) +-40 bp).
Per (truth, assembly, contig of the truth's chromosome; tumor contigs can be fused, chr3_chr13_hap1; unplaced
haplotype*-contigs are not used), per query side the best alignment ranked (covers, AS, clean) as in a2_parse.py:
covers = aligned query part spans the event window with >= 25 bases on each side; clean = no cs difference (mismatch /
insertion / deletion) inside the event window. Status: ALT (ALT clean, REF not) / REF / ambiguous (both clean) /
other (covered, neither clean = another allele at the site) / none (nothing covers the window or no alignment).
Extra columns: AS and target span (start0-end0 strand) of the best REF / ALT alignments, spans_agree (the two overlap; if
not, they come from different copies, e.g. a segmental duplication), and graph_win_t0/t1 = the graph window
(loci.tsv win_start0..win_end0) projected onto the contig through the alignment used (REF if it covers, else ALT; CIGAR
walk, linear extrapolation beyond the aligned part; ext_bases = how much was extrapolated) for s7_assembly_paths.py;
hap_diffs = the haplotype's differences from GRCh38 inside the event window read off the best REF alignment's cs when it
covers the window ('pos1:XG>A' base change, 'pos1:DSEQ' GRCh38 bases missing from the haplotype, 'pos1:ISEQ' haplotype bases
inserted after pos1), hap_len_change = their net length change, site_base = the haplotype base at an SNV truth ('-' deleted),
site_call = a coarser call that splits 'other' (see site_call()).
asm_seq_summary.tsv: per truth the two normal haplotype statuses, the event haplotype (asm_event) and the other one,
the tumor contig statuses, spans_disagree = (assembly:contig) whose REF / ALT best alignments are different copies.
Assumptions: the tumor assembly is read from the uncompressed copy of an earlier session (the .fasta.gz has no index);
best contig per normal haplotype = the single chrN_hapK contig; per truth the event haplotype comes from asm_event
(truth INFO HG008Nv62SOMATICVARIANT, 26 empty).
"""
import csv, re, sys
from collections import Counter, defaultdict
import pysam
D = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
FASTA = '/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta'
ASM = {'normal': '/scratch/jshen/data/HG008_GIAB/HG008N_curatedv6_250714_bothhaps_polished6.2.fasta.gz',
       'tumor': '/scratch/jshen/data/pansoma_net_v2_runs/indel_graph_paths_20260930/assembly/HG008T_v3.2.fasta'}
CS = re.compile(r'(:\d+|\*[a-z][a-z]|[+-][a-z]+)')
CIG = re.compile(r'(\d+)([MIDNSHP=X])')
STATS = ['ALT', 'REF', 'ambiguous', 'other', 'none']


def queries():
    fa = pysam.FastaFile(FASTA)
    with open(f'{D}/asm_queries.fa', 'w') as fq, open(f'{D}/asm_queries.tsv', 'w') as ft:
        ft.write('truth_id\tchrom\tvcf_pos\tvcf_ref\tvcf_alt\tkind2\tasm_event\tq_start0\tref_len\talt_len\tev_ref_lo\tev_ref_hi\t'
                 'ev_alt_lo\tev_alt_hi\twin_start0\twin_end0\n')
        for r in csv.DictReader(open(f'{D}/loci.tsv'), delimiter='\t'):
            pos0, vref, valt = int(r['vcf_pos']) - 1, r['vcf_ref'].upper(), r['vcf_alt'].upper()
            lo, hi = int(r['ev_lo0']), int(r['ev_hi0'])
            F = max(100, pos0 - lo + 40, hi - (pos0 + len(vref)) + 40)
            s0 = pos0 - F
            refq = fa.fetch(r['chrom'], s0, pos0 + len(vref) + F).upper()
            altq = refq[:F] + valt + refq[F + len(vref):]
            lt = len(valt) - len(vref)
            fq.write(f">{r['truth_id']}|ref\n{refq}\n>{r['truth_id']}|alt\n{altq}\n")
            ft.write('\t'.join(map(str, [r['truth_id'], r['chrom'], r['vcf_pos'], vref, valt, r['kind2'], r['asm_event'], s0, len(refq),
                                         len(altq), lo - s0, hi - s0, lo - s0, hi - s0 + lt, r['win_start0'], r['win_end0']])) + '\n')


def evaluate(flag, cigar, cs, qlen, window):
    """(covers, clean) of one alignment; window [lo, hi) in original query coordinates (a2_parse.py)."""
    ops = CIG.findall(cigar)
    clip_l = int(ops[0][0]) if ops and ops[0][1] in 'SH' else 0
    clip_r = int(ops[-1][0]) if ops and ops[-1][1] in 'SH' else 0
    lo, hi = window
    if flag & 16:                                     # reverse strand: cs runs along the reverse complement
        lo, hi = qlen - hi, qlen - lo
    a0, a1 = clip_l, qlen - clip_r
    covers = a0 <= lo - 25 and a1 >= hi + 25
    q, clean = a0, True
    for op in CS.findall(cs):
        if op[0] == ':':
            q += int(op[1:])
        elif op[0] == '*':
            if lo <= q < hi: clean = False
            q += 1
        elif op[0] == '+':
            if q < hi and q + len(op) - 1 > lo: clean = False
            q += len(op) - 1
        elif lo <= q <= hi:                           # deletion from the query, between q-1 and q
            clean = False
    return covers, clean


def diffs(flag, cigar, cs, qlen, window):
    """Haplotype differences from the GRCh38 (REF) context inside window [lo, hi) (original query coordinates):
    [(query offset, 'X' GRCh38>hap base | 'D' GRCh38 bases the haplotype lacks | 'I' hap bases inserted after the offset)]."""
    ops = CIG.findall(cigar)
    q = int(ops[0][0]) if ops and ops[0][1] in 'SH' else 0
    lo, hi = window
    rv = bool(flag & 16)
    rc = lambda x: x.upper().translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
    out = []
    for op in CS.findall(cs):
        if op[0] == ':':
            q += int(op[1:]); continue
        if op[0] == '*':
            o, v = (qlen - 1 - q, f'{rc(op[2])}>{rc(op[1])}') if rv else (q, f'{op[2].upper()}>{op[1].upper()}'); q += 1
            if lo <= o < hi: out.append((o, 'X', v))
        elif op[0] == '+':
            n = len(op) - 1
            o, v = (qlen - q - n, rc(op[1:])) if rv else (q, op[1:].upper()); q += n
            if o < hi and o + n > lo: out.append((o, 'D', v))
        else:
            o, v = (qlen - q - 1, rc(op[1:])) if rv else (q - 1, op[1:].upper())
            if lo - 1 <= o < hi: out.append((o, 'I', v))
    return sorted(out)


def project(pos0, flag, cigar, qlen, x):
    """Target coordinate of query offset x (original query orientation; may lie outside the aligned part) -> (t, ext)."""
    ops = [(int(n), op) for n, op in CIG.findall(cigar)]
    if flag & 16:
        x = qlen - 1 - x
    clip_l = ops[0][0] if ops[0][1] in 'SH' else 0
    if x < clip_l:
        return pos0 - (clip_l - x), clip_l - x
    t, q = pos0, clip_l
    for n, op in ops:
        if op in 'SH':
            continue
        if op in 'M=X':
            if q <= x < q + n:
                return t + x - q, 0
            q += n; t += n
        elif op == 'I':
            if q <= x < q + n:
                return t, 0
            q += n
        elif op in 'DN':
            t += n
    return t - 1 + (x - q + 1), x - q + 1             # beyond the aligned end


def best_alignments(assembly, Q):
    """(truth_id, contig) -> side -> (covers, AS, clean, start0, end0, strand, flag, cigar)."""
    best = defaultdict(dict)
    for line in open(f'{D}/asm_{assembly}.sam'):
        if line[0] == '@':
            continue
        f = line.rstrip('\n').split('\t')
        if int(f[1]) & 4:
            continue
        tid, side = f[0].split('|'); q = Q[tid]
        if q['chrom'] not in f[2].split('_'):
            continue
        tags = dict(t.split(':', 2)[0::2] for t in f[11:])
        if 'cs' not in tags:
            continue
        win = (int(q['ev_ref_lo']), int(q['ev_ref_hi'])) if side == 'ref' else (int(q['ev_alt_lo']), int(q['ev_alt_hi']))
        covers, clean = evaluate(int(f[1]), f[5], tags['cs'], int(q[f'{side}_len']), win)
        span = sum(int(n) for n, op in CIG.findall(f[5]) if op in 'MDN=X')
        cand = (covers, int(tags.get('AS', 0)), clean, int(f[3]) - 1, int(f[3]) - 1 + span, '-' if int(f[1]) & 16 else '+', int(f[1]), f[5], tags['cs'])
        cur = best[(tid, f[2])].get(side)
        if cur is None or cand[:3] > cur[:3]:
            best[(tid, f[2])][side] = cand
    return best


def status(r, a):
    if a[0] and a[2] and not (r[0] and r[2]): return 'ALT'
    if r[0] and r[2] and not (a[0] and a[2]): return 'REF'
    if a[0] and r[0] and a[2] and r[2]: return 'ambiguous'
    return 'other' if a[0] or r[0] else 'none'


def site_call(q, st, dv, sb, covered):
    """Coarser allele call at the truth: SNV = the haplotype base (ALT / REF / other base / deleted);
    INDEL = ALT / REF, else the net length change in the event window vs the truth's ('same length change' / 'GRCh38 length' /
    'other length change')."""
    if not covered and st not in ('ALT', 'REF'):
        return 'not covered' if st == 'none' else st
    if q['kind2'] == 'SNV' and sb:
        return 'ALT' if sb == q['vcf_alt'] else 'REF' if sb == q['vcf_ref'] else 'deleted' if sb == '-' else 'other base'
    if st in ('ALT', 'REF', 'ambiguous'):
        return st
    n = sum(len(v) * (t == 'I') - len(v) * (t == 'D') for o, t, v in dv)
    return 'same length change' if n == len(q['vcf_alt']) - len(q['vcf_ref']) else 'GRCh38 length' if n == 0 else 'other length change'


def parse():
    Q = {r['truth_id']: r for r in csv.DictReader(open(f'{D}/asm_queries.tsv'), delimiter='\t')}
    NONE = (False, '', False, '', '', '', 0, '', '')
    rows = defaultdict(list)                           # truth_id -> [(assembly, contig, status)]
    with open(f'{D}/asm_seq_status.tsv', 'w') as w:
        w.write('truth_id\tassembly\tcontig\tstatus\tref_cover\tref_clean\talt_cover\talt_clean\tref_AS\talt_AS\tref_span\talt_span\t'
                'spans_agree\tgraph_win_t0\tgraph_win_t1\text_bases\tproj_side\thap_diffs\thap_len_change\tsite_base\tsite_call\n')
        for assembly in ASM:
            best = best_alignments(assembly, Q)
            for (tid, contig), b in sorted(best.items(), key=lambda x: (int(x[0][0]), x[0][1])):
                q = Q[tid]; r = b.get('ref', NONE); a = b.get('alt', NONE)
                st = status(r, a)
                agree = '' if not (r[3] != '' and a[3] != '') else r[3] < a[4] and a[3] < r[4]
                side = 'ref' if r[0] or not a[0] and r[3] != '' else 'alt'
                al = b[side]
                lt = int(q['alt_len']) - int(q['ref_len']) if side == 'alt' else 0
                s0 = int(q['q_start0'])
                x0 = int(q['win_start0']) - s0; x1 = int(q['win_end0']) - s0 + lt
                (t0, e0), (t1, e1) = (project(al[3], al[6], al[7], int(q[f'{side}_len']), x) for x in (x0, x1))
                dv, sb = [], ''
                if r[0]:
                    dv = diffs(r[6], r[7], r[8], int(q['ref_len']), (int(q['ev_ref_lo']), int(q['ev_ref_hi'])))
                    if q['kind2'] == 'SNV':
                        F = int(q['vcf_pos']) - 1 - s0
                        sb = next((v[-1] for o, t, v in dv if t == 'X' and o == F), None) or \
                            ('-' if any(t == 'D' and o <= F < o + len(v) for o, t, v in dv) else q['vcf_ref'])
                w.write('\t'.join(map(str, [tid, assembly, contig, st, r[0], r[2], a[0], a[2], r[1], a[1],
                                            f'{r[3]}-{r[4]}{r[5]}' if r[3] != '' else '', f'{a[3]}-{a[4]}{a[5]}' if a[3] != '' else '',
                                            agree, min(t0, t1), max(t0, t1), e0 + e1, side,
                                            ';'.join(f'{s0 + o + 1}:{t}{v}' for o, t, v in dv), sum(len(v) * (t == 'I') - len(v) * (t == 'D') for o, t, v in dv) if r[0] else '',
                                            sb, site_call(q, st, dv, sb, r[0])])) + '\n')
                rows[tid].append((assembly, contig, st if agree is not False else st + '*', site_call(q, st, dv, sb, r[0])))
    with open(f'{D}/asm_seq_summary.tsv', 'w') as w:
        w.write('truth_id\tchrom\tvcf_pos\tvcf_ref\tvcf_alt\tkind2\tnormal_hap1\tnormal_hap2\tevent_hap\tevent_hap_status\t'
                'other_hap_status\tnormal_genotype\ttumor_statuses\ttumor_any_ALT\tspans_disagree\tnormal_hap1_site\tnormal_hap2_site\n')
        for tid, q in Q.items():
            n = {c: s.rstrip('*') for a, c, s, x in rows[tid] if a == 'normal'}
            ns_ = {c: x for a, c, s, x in rows[tid] if a == 'normal'}
            h1, h2 = n.get(f"{q['chrom']}_hap1", 'none'), n.get(f"{q['chrom']}_hap2", 'none')
            ev = q['asm_event'].split(':')[0]
            evs = {'hap1': h1, 'hap2': h2}.get(ev[-4:], '') if ev else ''
            oth = {'hap1': h2, 'hap2': h1}.get(ev[-4:], '') if ev else ''
            tum = [f'{c}:{s.rstrip("*")}' for a, c, s, x in rows[tid] if a == 'tumor']
            dis = ','.join(f'{a}:{c}' for a, c, s, x in rows[tid] if s.endswith('*'))
            w.write('\t'.join(map(str, [tid, q['chrom'], q['vcf_pos'], q['vcf_ref'], q['vcf_alt'], q['kind2'], h1, h2, ev, evs, oth,
                                        f'{h1}/{h2}', ';'.join(tum), any(t.endswith(':ALT') for t in tum), dis,
                                        ns_.get(f"{q['chrom']}_hap1", 'not covered'), ns_.get(f"{q['chrom']}_hap2", 'not covered')])) + '\n')
    S = list(csv.DictReader(open(f'{D}/asm_seq_summary.tsv'), delimiter='\t'))
    for kd in ('SNV', 'INDEL'):
        s = [x for x in S if x['kind2'] == kd]
        print(f'== {kd} {len(s)}')
        for col in ('normal_hap1', 'normal_hap2', 'event_hap_status', 'other_hap_status'):
            print(f'  {col:17s}', {k: sum(x[col] == k for x in s) for k in STATS + ['']})
        for col in ('normal_hap1_site', 'normal_hap2_site'):
            print(f'  {col:17s}', dict(Counter(x[col] for x in s).most_common()))
        print('  a normal hap ALT:', sum('ALT' in (x['normal_hap1'], x['normal_hap2']) for x in s), '| both:',
              sum(x['normal_hap1'] == x['normal_hap2'] == 'ALT' for x in s), '| tumor any ALT:', sum(x['tumor_any_ALT'] == 'True' for x in s))
    R = list(csv.DictReader(open(f'{D}/asm_seq_status.tsv'), delimiter='\t'))
    print('spans_agree:', dict(Counter((x['assembly'], x['spans_agree']) for x in R)))


if __name__ == '__main__':
    {'queries': queries, 'parse': parse}[sys.argv[1]]()
