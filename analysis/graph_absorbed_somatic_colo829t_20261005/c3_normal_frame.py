"""COLO829T step 3: whole tandem-array alleles of the COLO829BL (matched normal) haplotypes and the patient-frame event.

Port of analysis/graph_absorbed_somatic_20261001/s8b_array_alleles.py (array, common anchors, assembly vs dipcall) and
audit/b5_normal_frame.py (repeat class in the patient's own sequence). COLO829T has no truth INFO event in the normal's
frame (HG008's HG008Nv62SOMATICVARIANT), so the event is derived here from the normal haplotypes.
Loci: the 497 absorbed truths of $D/variant_set.tsv (c0) and every COLO829T chr1-22 truth INDEL of per_truth_COLO829T.tsv
(1,912, the background for absorption rates; the 430 absorbed INDELs are among them) = 1,979 loci (SNV background: c1,
GRCh38 frame only).
Inputs: GRCh38 no-alt FASTA; COLO829BL verkko 2.1 ragtag scaffolds hapX (chr1-22, X) / hapY (chr1-22, Y), read only from
/scratch/qfu (never written); dipcall of COLO829BL vs GRCh38 (dip.vcf.gz + .tbi, dip.bed; run as 'run-dipcall GRCh38
hapY.fa hapX.fa' on the raw verkko contigs, so GT = hapY|hapX: checked below on phased het SNVs).
Method
 1 array [arr_lo, arr_hi) (GRCh38, 0-based half-open), s8b rules: exact periodic stretches (seq[i] == seq[i + p] runs)
   touching it with 1-bp tolerance, round 0 from the event core (SNV base / INDEL bases after the padding base) with
   period 1-6 length >= max(2p, p + 3) and period 7-60 >= 2 copies, then period 1-60 >= max(2p, 10) until stable, cap 3 kb;
   united with the s1 event window (period 1-6 stretches touching the event, +-5 bp).
 1b imperfect tandem repeat (audit fix; the exact-period array misses minisatellites / VNTRs with indels between copies):
   GRCh38 +-600 bp around the truth; offset i is 'tandem' when its 12-mer recurs 7-150 bp away (consecutive occurrences),
   flags smoothed over +-15 bp (>= 60 %), the run holding the truth (within 10 bp) >= 40 bp; ivntr_spacing = the commonest
   recurrence distance in it; in_ivntr = spacing >= 7 and >= 20 bp of the run outside the array (the independent audit's a6
   detector and thresholds).
 2 placement: query = GRCh38 [arr_lo - 5 kb, arr_hi + 5 kb) ('queries' mode -> $D/c3_queries.fa; job $T/c3_minimap2.sbatch:
   minimap2 2.28 -x asm5 -c --secondary=yes -N 10 vs each ragtag primary FASTA -> $D/c3_{hapX,hapY}.paf). Per hap, the hit
   on the scaffold of the same chromosome whose query span covers the array +-200 bp, most matching bases first (the HG008
   audit a3 rule); its cg CIGAR projects the array centre; the scaffold stretch centre +-(2.5 kb + array length) is cut and
   oriented to GRCh38 (reverse complement on '-' hits). n_cover_other = covering hits on other scaffolds (duplication flag).
   When no hit covers the array +-200 bp (the asm5 chain often breaks inside the repeat), place: cover0 = a hit covering
   the array itself; bridge = two colinear same-strand hits on both sides of the array (target gap -500..20,000 bp), the
   stretch spans both inner ends +-(2.5 kb + array length).
   Fallback contigs: the same queries vs the ragtag unlocalized FASTAs (login node, $D/c3_{hapX,hapY}_unloc.paf) rescue no
   locus (hapX 0 hits, hapY only secondary copies), and the raw verkko contigs hold the same sequence as primary +
   unlocalized (total length differs by the ragtag N gaps only), so neither is used for placement; unloc_better = an
   unlocalized contig covers the array +-200 bp with more matching bases than the chosen primary hit (wrong-copy flag).
 2b paralog fix (audit): a hap copy is suspect when its ragtag hit has MAPQ < 20, another scaffold also covers the array
   +-200 bp, an isolated phased het SNV of dipcall (step 3 phase check on the ragtag copies) contradicts it, or it has no
   ragtag placement. A suspect copy is replaced by the contig copy dipcall itself aligned there: dipcall's primary
   alignment (hap1.paf = hapY, hap2.paf = hapX, tp:A:P) covering the array +-200 bp (highest MAPQ, then longest), the array
   centre projected through its cg CIGAR, the raw verkko contig centre +-(2.5 kb + array length) oriented to GRCh38
   (copy 'dipcall'). A suspect copy with no such alignment is kept, unless dipcall's phase contradicts it: then the hap is
   NA (copy 'dropped', place_na dropped_paralog). The phase check is redone on the copies used (phase_mismatch). GRCh38 24-mers ending at a <= arr_lo / starting at b >= arr_hi, nearest first up to 400 bp out,
   not low-complexity (< 3 bases, or a period <= 6 matching >= 80 % of the k-mer), once in GRCh38 +-5.5 kb, exactly once in
   each placed hap stretch (right anchor after the left one); nearest one both haps have, else the nearest one hap has
   (the other hap is then NA 'anchor'; no hap: the nearest GRCh38 ones, dipcall only). Audit fix: when a placed hap lacks
   that pair, the search is retried out to 2 kb requiring every placed hap (anchor_wide; e.g. a hap deletion or VNTR length
   difference over the near anchors). S_h = the hap bases between the anchors, R = GRCh38[a, b), A = R with the truth.
   asm call ALT (S == A) / REF (S == R) / other; dipcall: R with the hap's PASS records applied by GT (hapY = GT[0],
   hapX = GT[1]); NA if [a, b) is not inside one dip.bed interval, a non-PASS record overlaps, a record crosses a or b, a GT
   is missing or applied records overlap; asm_eq_dip. Hap sequence = assembly S, else the dipcall one (src).
 4 flank: germline differences of S_h from R outside the array (an anchor pushed past a hap SNP) are not part of the array
   allele: S_h vs R by affine-gap global alignment (mismatch 4, gap 5 + 1/base; after the common prefix / suffix; HG008
   audit a3 aligner), maximal difference runs touching [arr_lo, arr_hi] are the array allele, the rest flank
   (n_flank_diffs); core_h = R with the array runs only. allele_h: REF / ALT / germline_len (other length = a germline
   repeat-length allele) / germline_seq (same length, other bases) / NA; len_change_h = len(core_h) - len(R);
   total_len_change_h = len(S_h) - len(R), flank_len_h = their difference (flag flank_len>=4: in a long / imperfect repeat
   the aligner can put repeat-copy indels at the array edge or anchor, and len_change then understates the hap's change).
   alt_len_haps: INDEL truths, haps whose core has exactly len(A) but is not A (the normal may carry the ALT length).
 5 patient-frame event (the somatic change relative to the normal): event hap = the hap whose core is closest to A
   (Levenshtein); tie -> the one equal to R, else hapX (event_tie). core == A: normal_carries_ALT (germline-identical, no
   somatic change visible). Else the edit core -> A after trimming the common suffix, then prefix (left-normalized): one
   side empty -> INDEL (len change = len(A) - len(core); unit u = minimal period of the inserted / deleted bases X);
   both one base -> SNV; anything else -> complex. Context = the hap's own 300 bp beyond the anchors (GRCh38 when the hap
   sequence is dipcall-only) + core.
   INDEL tract L = the maximal exact period-u stretch of the normal that contains the deleted bases (insertion: the
   stretch of core + X containing X, minus len(X)) = the normal's repeat length; tandem = a unit of X sits next to the
   event. Class (b5 / s14 thresholds): HP>=7 (u 1, L >= 7) / HP4-6 (u 1, L 4-6) / STR2-6>=3copies (u 2-6, L >= 3u) /
   VNTR>6 (u > 6, X repeated next to the event, b5 rule) / in_STR (otherwise, inside a period 1-6 stretch >= 10 bp
   touching the event) / none. pf_class_b5 = the b5 rule verbatim (L = longest period-u stretch merely touching the event,
   in_STR only from that period-u stretch for u <= 6); it differs where a repeat of other bases / another motif touches
   the event (b5 then calls e.g. a non-tandem G deletion next to A x 10 HP>=7) and where only another period reaches 10 bp.
   SNV event: b5 SNV rule (longest period 1-6 stretch touching the base), label 'SNV in <class>' or 'SNV'.
   Units changed = |len change| / u (repeat INDEL classes only). pf_tract_b5 = the b5 rule's tract (HG008 comparison).
 5b imperfect VNTR (audit fix): a 'none' event with in_ivntr -> pf_class imperfect_VNTR, label 'imperfect VNTR' / 'SNV in
   imperfect VNTR' (pf_class_b5 unchanged: the b5 rule has no such check).
 5c single hap (audit fix): only one hap has a sequence and its event is not normal_carries_ALT -> pf_label 'single hap'
   (the missing hap may carry the ALT; the one-hap label in pf_label_one_hap).
 6 other_hap_* = the same derivation from the other hap: Levenshtein charges a k-bp indel k, so the closest hap can give
   'complex' where the other one gives a single INDEL / SNV (flag other_hap_single_event; the rule is kept, the
   alternative recorded). flags (per locus): <hap>:NA / dip_only / asm_ne_dip / mapq<20 / other_scaffold_hit / cover0 /
   bridge / dipcall_copy / dropped_paralog / flank_diffs>=3 / flank_len>=4, wide_anchor, tie_REF_chosen / tie_other,
   other_hap_single_event.
Outputs: $D/normal_frame.tsv (one row per locus; sequences included, '' when NA), $D/c3_dipcall_phase.tsv (GT-order check:
isolated phased het SNVs within the hap stretches, 21-mer REF / ALT context found in hapX / hapY), counts printed.
Usage: python c3_normal_frame.py queries | python c3_normal_frame.py [id,...] (ids: test printout only) |
       python c3_normal_frame.py show id,... (prints S_hapX, S_hapY, R, A and the event from normal_frame.tsv)
Assumptions: the same-chromosome scaffold hit holds the right copy unless suspect (step 2b); for suspect copies dipcall's
long colinear alignment holds the patient's own copy (not checked against a segmental-duplication track); exact periodic
runs miss imperfect copies (both haps are still compared over the same stretch; step 1b only flags / relabels 'none'); difference runs are placed by one alignment (gap placement inside a short repeat at the array
edge can move a flank difference into the array); event_tie loci have two equally close haps and the event is the one of
the chosen hap.
Results (2026-10-05, after the audit fixes; job 380926 minimap2 5 min 40 s, peak RSS 12.7 GB by /usr/bin/time, sacct
MaxRSS 7.4 GB; this script on the login node 3 min 54 s, 0.26 GB, $T/c3_normal_frame.fix.login.log; the first version's
output is frozen in $D/prefix_20261005/, its log $T/c3_normal_frame.login.log)
 - 1,979 loci (INDEL 1,912 of which 430 absorbed; absorbed SNV 67), all with GRCh38 anchors, no array capped.
 - placement (ragtag) hapX cover 1,971 / bridge 6 / cover0 1 / no hit 1; hapY 1,974 / 2 / 2 / 1; unloc_better 0.
 - 2b paralog fix: suspect hapX 35 / hapY 26 copies; dipcall copy used hapX 29 / hapY 22 (35 loci, 12 absorbed; hap
   sequence differs from the ragtag copy at 10 loci); suspect without a dipcall alignment kept hapX 3 (35157, 35158, 43396)
   / hapY 3 (33233, 35151, 35160), dropped (dipcall phase contradicted) hapX 2 (110, 35171); 41860 has no copy at all. Isolated phased het SNVs (6,232 distinct): consistent 6,068,
   swapped 0, other 0, undetermined 164 (before the fix: rows 6,002 / 0 / 122 / 169; loci with a contradicted copy 10 -> 0).
   Labels changed by it: 2561 none -> SNV, 2562 STR -> SNV, 2563 HP4-6 -> normal carries ALT, 35166 none -> normal
   carries ALT (absorbed); 35160 none -> normal carries ALT, 35167 normal carries ALT -> in_STR (background). The 8 haps
   whose copy changed reproduce the independent audit's dipcall-copy alleles (audit a7) exactly.
 - wide anchors (2 kb retry) 5 loci: 29401, 30786, 34621, 40078, 40329 (40329 none -> normal carries ALT, hapX ALT as in
   dipcall; 34621 stays HP>=7: hapX carries a 1,955-bp deletion over the homopolymer, dipcall record at 90210597).
 - Hap sequence NA: both haps 1142 (chr1:121746315), 15472 (chr5:49786031) (no common anchor), 41860 (chr20:25162371, no
   hit in either hap or in dipcall's alignments; absorbed); hapX only 110, 35171 (dropped_paralog), 29403 (anchor) ->
   'single hap' (3, all background).
 - 1b imperfect tandem repeat: in_ivntr 36 loci (INDEL 33, 9 absorbed; absorbed SNV 3), identical to audit a6; relabelled
   'none' -> imperfect VNTR: absorbed INDEL 29401, 39443, background 24415, 40334 (and 29403's one-hap label), absorbed
   SNV 2570, 6942, 6943 -> 'SNV in imperfect VNTR'.
 - GT order: 559 het loci (hapX != hapY, both asm + dip): asm == dip for both haps 540 as hapY|hapX, 0 swapped;
   asm_eq_dip False hapX 21 / hapY 13.
 - array alleles: REF hapX 1,471 / hapY 1,483, ALT 47 / 9, germline_seq 61 / 80, NA 6 / 3; flank_len>=4 hapX 19 / hapY 20.
 - event basis both 1,973 (tie_same_seq 1,430, tie_REF_chosen 59, tie_other 14), only_hapY 3, none 3; flags on 168 loci.
 - checks: every derived event re-applied to the event hap gives A (INDEL 1,625, SNV 126, complex 59; 110 event haps with
   flank differences not re-checked); about 40 loci read by hand (the audits' loci plus 2570, 6942, 32518, 1166, 2569).
 Patient-frame label of the absorbed truths (per platform perfect set; union):
   INDEL                    fiberseq 218   ONT 170      Illumina 396  union 430
   HP>=7                    79 (36.2%)    51 (30.0%)   235 (59.3%)   241 (56.0%)
   STR2-6>=3copies          57 (26.1%)    44 (25.9%)    75 (18.9%)    83 (19.3%)
   VNTR>6                   1             1              1             1
   in_STR / imperfect VNTR  12 / 2        13 / 1        15 / 1        17 / 2
   SNV in STR / HP4-6       10 / 1        11 / 1        11 / 1        15 / 1
   SNV                      6             5             5             6
   none                     2             3             3             3
   complex                  21 (9.6%)     18 (10.6%)    21 (5.3%)     27 (6.3%)
   normal carries ALT       27 (12.4%)    22 (12.9%)    27 (6.8%)     33 (7.7%)
   NA                       0             0             1             1
   SNV (union 67): SNV 38, SNV in STR 10, SNV in HP4-6 7, SNV in HP>=7 3, SNV in imperfect VNTR 3, complex 3, in_STR 2,
   normal carries ALT 1.
   absorbed repeat INDELs: 1 unit 316, 2-3 units 7, > 3 units 2; HP>=7 normal tract 7-9 6, 10-14 86, 15-19 94, 20-29 55.
   alt_len_haps (a hap with the ALT length, not the ALT), absorbed INDELs not 'normal carries ALT': 35 (SNV in STR 15,
   complex 9, SNV 6, HP>=7 2, none 1, STR 1, SNV in HP4-6 1).
 All 1,912 chr1-22 truth INDELs, absorbed / all by patient-frame label: HP>=7 241 / 1,103 (21.8%), STR 83 / 189 (43.9%),
   HP4-6 0 / 67, VNTR>6 1 / 13, in_STR 17 / 92, imperfect VNTR 2 / 4, SNV in STR 15 / 46, SNV in HP>=7 0 / 16, SNV 6 / 11,
   complex 27 / 72, none 3 / 235 (1.3%), normal carries ALT 33 / 55, single hap 0 / 3, NA 1 / 3.
"""
import bisect, collections, csv, gzip, re, sys
import numpy as np
import pysam
csv.field_size_limit(10 ** 9)
D = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005'
TRUTH = '/scratch/jshen/Github/Pansoma/analysis/tensor_recall_20260930/per_truth_COLO829T.tsv'
GR = '/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta'
Q = '/scratch/qfu/COLO829BL_DSA'
HAP = {'hapX': f'{Q}/ragtag_hg38/ragtag_SMHTCOLO829BL-X-X-M45-XX-uwsc-SMASFTH1B2J6-verkko_2.1_hapX_primary_normalized.fa',
       'hapY': f'{Q}/ragtag_hg38/ragtag_SMHTCOLO829BL-X-X-M45-XX-uwsc-SMASFK4V33V8-verkko_2.1_hapY_primary_normalized.fa'}
GTI = {'hapY': 0, 'hapX': 1}                                    # dipcall GT hap1|hap2 = hapY|hapX
DIP = f'{Q}/dipcall_hg38/dipcall_hg38.dip'
RAW = {'hapX': f'{Q}/SMHTCOLO829BL-X-X-M45-XX-uwsc-SMASFTH1B2J6-verkko_2.1_hapX.fa',          # raw verkko contigs haplotype1-*
       'hapY': f'{Q}/SMHTCOLO829BL-X-X-M45-XX-uwsc-SMASFK4V33V8-verkko_2.1_hapY.fa'}          # haplotype2-*
DPAF = {'hapY': f'{Q}/dipcall_hg38/dipcall_hg38.hap1.paf.gz', 'hapX': f'{Q}/dipcall_hg38/dipcall_hg38.hap2.paf.gz'}
K, CTX, SHIFT, SHIFT_WIDE, MAXARR, QF, FL, CAP = 24, 2500, 400, 2000, 3000, 5000, 300, 2_500_000
IVF, IVK, IVMIN, IVOUT = 600, 12, 40, 20                        # imperfect tandem repeat: window, k, run, bp beyond the array
AUTO = {f'chr{i}' for i in range(1, 23)}
RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
rc = lambda s: s.translate(RC)[::-1]
CG = re.compile(r'(\d+)([MID])')


# ---- 1 array (s1 / s8b) ----
def stretches(g, g0):
    """[(start, end, p)] exact periodic stretches of g (absolute coordinates), period 1-60, >= 2 copies and >= p + 3."""
    a, out = np.frombuffer(g.encode(), dtype=np.uint8), []
    for p in range(1, 61):
        d = np.diff(np.concatenate(([0], (a[:-p] == a[p:]).astype(np.int8), [0])))
        s, e = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
        m = (e - s) >= max(p, 3)
        out += [(g0 + i, g0 + j + p, p) for i, j in zip(s[m], e[m])]
    return out


def grow(st, lo, hi):
    ok1 = lambda x, y, p: (y - x >= max(2 * p, p + 3)) if p <= 6 else y - x >= 2 * p
    ok2 = lambda x, y, p: y - x >= max(2 * p, 10)
    for rnd in range(100):
        f = ok1 if rnd == 0 else ok2
        n = [(x, y) for x, y, p in st if x <= hi + 1 and y >= lo - 1 and f(x, y, p)]
        nlo, nhi = min([lo] + [x for x, _ in n]), max([hi] + [y for _, y in n])
        if nhi - nlo > MAXARR:
            return lo, hi, True
        if (nlo, nhi) == (lo, hi) and rnd > 0:
            return lo, hi, False
        lo, hi = nlo, nhi
    return lo, hi, True


def event_window(st, lo, hi):
    """s1: union of period 1-6 stretches of length >= max(2p, p + 3) touching the core, +-5 bp."""
    L, R = lo, hi
    for x, y, p in st:
        if p <= 6 and y - x >= max(2 * p, p + 3) and x <= hi + 1 and y >= lo - 1:
            L, R = min(L, x), max(R, y)
    return L - 5, R + 5


def loci():
    vs = {r['truth_id']: r for r in csv.DictReader(open(f'{D}/variant_set.tsv'), delimiter='\t')}
    out = []
    for r in csv.DictReader(open(TRUTH), delimiter='\t'):
        kd = 'SNV' if r['kind'] == 'SNP' else 'INDEL'
        if r['truth_id'] in vs or (kd == 'INDEL' and r['chrom'] in AUTO):
            out.append(dict(truth_id=r['truth_id'], chrom=r['chrom'], vcf_pos=r['vcf_pos'], vcf_ref=r['vcf_ref'].upper(),
                            vcf_alt=r['vcf_alt'].upper(), kind2=kd, absorbed=r['truth_id'] in vs))
    assert sum(o['absorbed'] for o in out) == len(vs), 'variant_set truth missing from per_truth'
    return out


def array(gfa, L):
    ch, p0, vr, va = L['chrom'], int(L['vcf_pos']) - 1, L['vcf_ref'], L['vcf_alt']
    g0 = max(0, p0 - CTX - MAXARR); g = gfa.fetch(ch, g0, p0 + len(vr) + CTX + MAXARR).upper()
    assert g[p0 - g0:p0 - g0 + len(vr)] == vr, L['truth_id']
    core = (p0, p0 + 1) if len(vr) == 1 and len(va) == 1 else (p0 + 1, max(p0 + len(vr), p0 + 1))
    st = stretches(g, g0)
    lo, hi, capped = grow(st, *core)
    ev = event_window(st, *core)
    return g, g0, min(lo, ev[0]), max(hi, ev[1]), capped


# ---- 2/3 placement, anchors, dipcall (s8b / a3) ----
def project(h, x):
    """query offset x -> target position through the hit's cg CIGAR (HG008 audit a3)."""
    q, t = (h['qs'], h['ts']) if h['st'] == '+' else (h['ql'] - h['qe'], h['ts'])
    xx = x if h['st'] == '+' else h['ql'] - 1 - x
    for n, op in CG.findall(h['cg']):
        n = int(n)
        if op == 'M':
            if q <= xx < q + n:
                return t + xx - q
            q += n; t += n
        elif op == 'I':
            if q <= xx < q + n:
                return t
            q += n
        else:
            t += n
    return None


def lowc(k):
    return len(set(k)) < 3 or any(sum(k[i] == k[i + p] for i in range(K - p)) >= 0.8 * (K - p) for p in range(1, 7))


def occ(s, k, start=0):
    i = s.find(k, start)
    return [] if i < 0 else [i] if s.find(k, i + 1) < 0 else [i, s.find(k, i + 1)]


def anchors(g, g0, lo, hi, segs, shift=SHIFT, need_all=False):
    """(a, b, {hap: (left offset, right offset)}) - docstring 3 (s8b anchors with the two haps as the 'normals');
    need_all: only anchors every hap has (the wide retry), else (None, None, {})."""
    def pick(cands):
        best = None
        for pos, ok in cands:
            if len(ok) == len(segs):
                return pos, ok
            if not need_all and (best is None or len(ok) > len(best[1])):
                best = (pos, ok)
        return best or (None, {})                               # no hap has it: GRCh38 anchors only (dipcall-only S)
    def left():
        for e in range(lo, lo - shift, -1):
            k = g[e - K - g0:e - g0]
            if len(k) < K or 'N' in k or lowc(k) or len(occ(g, k)) != 1:
                continue
            yield e, {c: o[0] for c, s in segs.items() for o in [occ(s, k)] if len(o) == 1}
    a, lo_ok = pick(left())
    if a is None:
        return None, None, {}
    def right():
        for b in range(hi, hi + shift):
            k = g[b - g0:b - g0 + K]
            if len(k) < K or 'N' in k or lowc(k) or len(occ(g, k)) != 1:
                continue
            ok = {}
            for c, i in lo_ok.items():
                o = occ(segs[c], k, i + K)
                if len(o) == 1 and segs[c].find(k) == o[0]:
                    ok[c] = (i + K, o[0])
            yield b, ok
    b, ok = pick(right())
    return (a, b, ok) if b is not None else (a, None, {})


def edit(R, a, eds):
    out, cur = [], a
    for p0, r, x in sorted(eds):
        if p0 < cur or R[p0 - a:p0 - a + len(r)] != r:
            return None
        out += [R[cur - a:p0 - a], x]; cur = p0 + len(r)
    return ''.join(out) + R[cur - a:]


def read_dpaf(h):
    """dipcall's own primary contig -> GRCh38 alignments of hap h: {chrom: [(ts, te, contig, qs, qe, strand, mapq, cg)]}."""
    out = collections.defaultdict(list)
    with gzip.open(DPAF[h], 'rt') as f:
        for ln in f:
            x = ln.rstrip('\n').split('\t')
            if len(x) < 12 or x[5] not in AUTO or 'tp:A:P' not in x[12:]:
                continue
            cg = next((t[5:] for t in x[12:] if t.startswith('cg:Z:')), '')
            out[x[5]].append((int(x[7]), int(x[8]), x[0], int(x[2]), int(x[3]), x[4], int(x[11]), cg))
    return out


def tproject(x, t):
    """GRCh38 position t -> contig position through a dipcall PAF record (target = GRCh38, query = contig)."""
    ts, te, qn, qs, qe, st, mq, cg = x
    tp, qo = ts, 0
    for n, op in CG.findall(cg):
        n = int(n)
        if op == 'M':
            if tp <= t < tp + n:
                qo += t - tp; break
            tp += n; qo += n
        elif op == 'I':
            qo += n
        else:
            if tp <= t < tp + n:
                break
            tp += n
    return qs + qo if st == '+' else qe - 1 - qo


def imperfect_repeat(g, c, k=IVK):
    """(spacing, lo, hi) of an imperfect tandem repeat around offset c of g (HG008-independent audit a6 detector): offset i is
    'tandem' when its k-mer recurs 7-150 bp away (consecutive occurrences); flags smoothed over +-15 bp (>= 60 %); the run
    holding c (or within 10 bp), >= IVMIN bp; spacing = the commonest recurrence distance inside it (unit or a multiple)."""
    pos = collections.defaultdict(list)
    for i in range(len(g) - k + 1):
        pos[g[i:i + k]].append(i)
    flag, dist = np.zeros(len(g), dtype=np.int32), {}
    for km, ps in pos.items():
        if len(ps) < 2 or 'N' in km:
            continue
        for x, y in zip(ps, ps[1:]):
            if 7 <= y - x <= 150:
                flag[x:x + k] = 1; flag[y:y + k] = 1; dist[x] = y - x
    cs = np.concatenate(([0], np.cumsum(flag))); i = np.arange(len(g))
    l, r = np.clip(i - 15, 0, len(g)), np.clip(i + 16, 0, len(g))
    ok = (cs[r] - cs[l]) / (r - l) >= 0.6
    near = [j for j in range(max(0, c - 10), min(len(g), c + 11)) if ok[j]]
    if not near:
        return 0, c, c
    lo = hi = min(near, key=lambda j: abs(j - c))
    while lo > 0 and ok[lo - 1]: lo -= 1
    while hi < len(g) - 1 and ok[hi + 1]: hi += 1
    if hi + 1 - lo < IVMIN:
        return 0, c, c
    ds = collections.Counter(d for x, d in dist.items() if lo <= x <= hi)
    return (ds.most_common(1)[0][0] if ds else 0), lo, hi + 1


def dip_hap(R, a, b, recs, h, inbed):
    """(sequence, NA reason) for GT index h (s8b)."""
    if not inbed:
        return None, 'outside_dip_bed'
    eds = []
    for p0, r, alts, gt, ok in recs:
        if not ok:
            return None, 'dip_nonpass_record'
        g = gt[h] if gt and len(gt) > h else None
        if g is None:
            return None, 'dip_gt_missing'
        if g == 0 or alts[g - 1] == '*':
            continue
        if p0 < a or p0 + len(r) > b:
            return None, 'dip_record_crosses_anchor'
        eds.append((p0, r, alts[g - 1]))
    s = edit(R, a, eds)
    return (s, '') if s is not None else (None, 'dip_overlapping_records')


# ---- 4 flank vs array (a3 aligner) ----
def align(a, b, mm=4, go=5, ge=1):
    """affine-gap global alignment of a (hap) to b (GRCh38); columns (i_b or None, j_a or None) (HG008 audit a3)."""
    n, m, INF = len(a), len(b), 10 ** 9
    M = [[INF] * (m + 1) for _ in range(n + 1)]; X = [[INF] * (m + 1) for _ in range(n + 1)]; Y = [[INF] * (m + 1) for _ in range(n + 1)]
    M[0][0] = 0
    for i in range(1, n + 1):
        X[i][0] = go + ge * i
    for j in range(1, m + 1):
        Y[0][j] = go + ge * j
    for i in range(1, n + 1):
        ai, Mi, Mp, Xi, Xp, Yi, Yp = a[i - 1], M[i], M[i - 1], X[i], X[i - 1], Y[i], Y[i - 1]
        for j in range(1, m + 1):
            Mi[j] = min(Mp[j - 1], Xp[j - 1], Yp[j - 1]) + (0 if ai == b[j - 1] else mm)
            Xi[j] = min(Mp[j] + go + ge, Xp[j] + ge, Yp[j] + go + ge)
            Yi[j] = min(Mi[j - 1] + go + ge, Yi[j - 1] + ge, Xi[j - 1] + go + ge)
    i, j, cols = n, m, []
    st = min((M[n][m], 0), (X[n][m], 1), (Y[n][m], 2))[1]
    while i > 0 or j > 0:
        if st == 0:
            cols.append((j - 1, i - 1)); c = 0 if a[i - 1] == b[j - 1] else mm
            v = M[i][j] - c; i -= 1; j -= 1
            st = 0 if M[i][j] == v else 1 if X[i][j] == v else 2
        elif st == 1:
            cols.append((None, i - 1)); v = X[i][j]; i -= 1
            st = 1 if X[i][j] + ge == v else 0 if M[i][j] + go + ge == v else 2
        else:
            cols.append((j - 1, None)); v = Y[i][j]; j -= 1
            st = 2 if Y[i][j] + ge == v else 0 if M[i][j] + go + ge == v else 1
    return cols[::-1]


def trim(S, R, suffix_first=False):
    """(p, q): common prefix / suffix lengths (bounded so that they do not overlap)."""
    m = min(len(S), len(R))
    if suffix_first:
        q = 0
        while q < m and S[-1 - q] == R[-1 - q]: q += 1
        p = 0
        while p < m - q and S[p] == R[p]: p += 1
        return p, q
    p = 0
    while p < m and S[p] == R[p]: p += 1
    q = 0
    while q < m - p and S[-1 - q] == R[-1 - q]: q += 1
    return p, q


def diff_runs(S, R):
    """maximal difference runs of S vs R as (r0, r1, s_seq) in R coordinates; None if the middle is too long to align."""
    p, q = trim(S, R)
    s, r = S[p:len(S) - q], R[p:len(R) - q]
    if not s and not r:
        return []
    if not s or not r:
        return [(p, len(R) - q, s)]
    if len(s) * len(r) > CAP:
        return None
    out, r0, buf, rpos = [], None, [], 0
    for bi, aj in align(s, r):
        if bi is not None and aj is not None and s[aj] == r[bi]:
            if r0 is not None:
                out.append((p + r0, p + rpos, ''.join(buf))); r0, buf = None, []
        else:
            if r0 is None:
                r0 = rpos
            if aj is not None:
                buf.append(s[aj])
        if bi is not None:
            rpos = bi + 1
    if r0 is not None:
        out.append((p + r0, p + rpos, ''.join(buf)))
    return out


def lev(x, y):
    p, q = trim(x, y)
    x, y = x[p:len(x) - q], y[p:len(y) - q]
    if not x or not y:
        return len(x) + len(y)
    prev = list(range(len(y) + 1))
    for i, cx in enumerate(x, 1):
        cur = [i]
        for j, cy in enumerate(y, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (cx != cy)))
        prev = cur
    return prev[-1]


# ---- 5 patient-frame event (b5) ----
def tract(seq, lo, hi, periods):
    """b5: longest exact periodic stretch (period p) touching [lo, hi): (p, start, end)."""
    best = (0, lo, lo)
    for p in periods:
        i = max(0, lo - 200)
        while i < min(len(seq) - p, hi + 200):
            if seq[i] != seq[i + p]: i += 1; continue
            j = i
            while j < len(seq) - p and seq[j] == seq[j + p]: j += 1
            a, b = i, j + p
            if a <= hi and b >= lo and b - a > best[2] - best[1] and b - a >= 2 * p: best = (p, a, b)
            i = j + 1
    return best


def minper(s):
    for p in range(1, len(s) + 1):
        if len(s) % p == 0 and s[:p] * (len(s) // p) == s: return p
    return len(s)


def run_in(seq, i0, i1, u):
    """length of the maximal exact period-u stretch containing seq[i0:i1] (itself u-periodic)."""
    a, b = i0, i1
    while a > 0 and a - 1 + u < len(seq) and seq[a - 1] == seq[a - 1 + u]: a -= 1
    while b < len(seq) and b - u >= 0 and seq[b] == seq[b - u]: b += 1
    return b - a


def hp_class(u, L):
    return ('HP>=7' if u == 1 and L >= 7 else 'HP4-6' if u == 1 and L >= 4 else
            'STR2-6>=3copies' if 2 <= u <= 6 and L >= 3 * u else '')


def event(core, A, ctx_l, ctx_r):
    """patient-frame event core -> A (docstring 5)."""
    if core == A:
        return dict(pf_kind='normal_carries_ALT', pf_class='normal_carries_ALT', pf_class_b5='normal_carries_ALT',
                    pf_label='normal carries ALT')
    p, q = trim(core, A, suffix_first=True)
    x, y = core[p:len(core) - q], A[p:len(A) - q]
    seq = ctx_l + core + ctx_r; lo = len(ctx_l) + p; hi = lo + len(x)
    pad = seq[lo - 1]
    o = dict(pf_event=f'{p}:{pad + x}>{pad + y}' if not x or not y else f'{p}:{x}>{y}', pf_len_change=len(y) - len(x))
    sp, sa, sb = tract(seq, lo - (not x), hi + (not x), range(1, 7))
    o['pf_stretch'] = f'p{sp}:{sb - sa}' if sp else ''
    if len(x) == 1 and len(y) == 1:
        c = hp_class(sp, sb - sa) or ('STR2-6>=3copies' if sp >= 2 and sb - sa >= 3 * sp else 'none')
        return dict(o, pf_kind='SNV', pf_class=c, pf_class_b5=c, pf_label='SNV' if c == 'none' else f'SNV in {c}',
                    pf_cpg_ti=(x == 'C' and y == 'T' and seq[hi] == 'G') or (x == 'G' and y == 'A' and seq[lo - 1] == 'C'),
                    pf_titv='Ti' if {x, y} in ({'A', 'G'}, {'C', 'T'}) else 'Tv')
    if x and y:
        return dict(o, pf_kind='complex', pf_class='complex', pf_class_b5='complex', pf_label='complex')
    X = x or y; u = minper(X)
    if x:                                                       # deletion: the stretch of the normal holding X
        L = run_in(seq, lo, hi, u)
        tandem = seq[hi:hi + u] == X[:u] or seq[lo - u:lo] == X[-u:]
    else:                                                       # insertion: the stretch of normal + X holding X, minus X
        L = run_in(seq[:lo] + X + seq[lo:], lo, lo + len(X), u) - len(X)
        tandem = seq[lo:lo + u] == X[:u] or seq[lo - u:lo] == X[-u:]
    vntr = u > 6 and (seq[hi:hi + len(X)] == X or seq[lo - len(X):lo] == X)
    other = f'in_STR(p{sp})' if sb - sa >= 10 else 'none'
    c = hp_class(u, L) or ('VNTR>6' if vntr else other)
    bp, ba, bb = tract(seq, lo - (not x), hi + (not x), [u] if u <= 6 else range(1, 7))   # b5 verbatim
    cb = hp_class(u, bb - ba) or ('VNTR>6' if vntr else f'in_STR(p{bp})' if bb - ba >= 10 else 'none')
    rep = c in ('HP>=7', 'HP4-6', 'STR2-6>=3copies', 'VNTR>6')
    return dict(o, pf_kind='INDEL', pf_unit=u, pf_tract=L, pf_tandem=tandem, pf_class=c, pf_class_b5=cb, pf_tract_b5=bb - ba,
                pf_label='in_STR' if c.startswith('in_STR') else c, pf_units_changed=len(X) // u if rep else '')


# ---- main ----
def queries():
    gfa = pysam.FastaFile(GR); n = 0
    with open(f'{D}/c3_queries.fa', 'w') as w:
        for L in loci():
            g, g0, lo, hi, _ = array(gfa, L)
            q0 = max(0, lo - QF)
            w.write(f">{L['truth_id']}|{q0}\n{gfa.fetch(L['chrom'], q0, hi + QF).upper()}\n"); n += 1
    print('queries', n)


def read_paf(name):
    hits = collections.defaultdict(list)
    for ln in open(f'{D}/c3_{name}.paf'):
        f = ln.rstrip('\n').split('\t')
        tg = dict(x.split(':', 2)[0::2] for x in f[12:])
        t, q0 = f[0].split('|')
        hits[t].append(dict(q0=int(q0), ql=int(f[1]), qs=int(f[2]), qe=int(f[3]), st=f[4], tn=f[5], tl=int(f[6]), ts=int(f[7]),
                            te=int(f[8]), nm=int(f[9]), mapq=int(f[11]), tp=tg.get('tp', ''), cg=tg.get('cg', '')))
    return hits


def place(hs, ch, elo, ehi):
    """(tier, best hit, target points) on the scaffold ch: cover (array +-200 bp in one hit) / cover0 (array only) /
    bridge (two colinear same-strand hits on both sides of the array, target gap <= 20 kb); else (NA reason, None, None)."""
    own = [x for x in hs if x['tn'] == ch]
    for tier, m in (('cover', 200), ('cover0', 0)):
        c = sorted([x for x in own if x['qs'] <= elo - m and x['qe'] >= ehi + m], key=lambda x: -x['nm'])
        if c:
            tc = project(c[0], (elo + ehi) // 2)
            return (tier, c[0], (tc, tc)) if tc is not None else ('no_projection', None, None)
    br = []
    for l in own:
        for r in own:
            if l is r or l['st'] != r['st'] or not (l['qs'] <= elo - 200 and l['qe'] < ehi + 200 and r['qs'] > elo - 200 and r['qe'] >= ehi + 200):
                continue
            tl, tr = project(l, l['qe'] - 1), project(r, r['qs'])
            if tl is not None and tr is not None and -500 <= (tr - tl if l['st'] == '+' else tl - tr) <= 20000:
                br.append((l['nm'] + r['nm'], l, (min(tl, tr), max(tl, tr))))
    if br:
        _, l, pts = max(br, key=lambda z: z[0])
        return 'bridge', l, pts
    return ('no_hits' if not hs else 'no_hit_on_scaffold' if not own else 'no_covering_hit'), None, None


def phase_check(g, g0, segs, near):
    """isolated PASS phased het SNVs (no other record within 10 bp) in `near`: per SNV the 21-mer REF / ALT context found
    once in each placed hap stretch -> obs[h] '0' / '1' / '?' and the GT expectation exp[h] (hapY = GT[0], hapX = GT[1])."""
    out = []
    for r in near:
        gt = r.samples[0]['GT']
        if r.filter.keys() or None in gt or sorted(gt) != [0, 1] or len(r.ref) != 1 or len(r.alts[0]) != 1:
            continue
        q0 = r.pos - 1; cx = g[q0 - 10 - g0:q0 + 11 - g0]
        if any(z.pos - 1 - 10 <= q0 < z.pos - 1 + len(z.ref) + 10 for z in near if z.pos != r.pos):
            continue
        obs = {}
        for h in segs:
            sn = [cx[:10] + al + cx[11:] in segs[h] for al in (r.ref.upper(), r.alts[0].upper())]
            obs[h] = '?' if sum(sn) != 1 else str(sn.index(True))
        out.append(dict(r=r, gt=gt, obs=obs, exp={h: str(gt[GTI[h]]) for h in HAP}))
    return out


def run(ids=None):
    gfa = pysam.FastaFile(GR); fas = {h: pysam.FastaFile(p) for h, p in HAP.items()}
    raw = {h: pysam.FastaFile(p) for h, p in RAW.items()}
    hits = {h: read_paf(h) for h in HAP}
    unloc = {h: read_paf(f'{h}_unloc') for h in HAP}
    dpaf = {h: read_dpaf(h) for h in HAP}
    dvcf = pysam.VariantFile(f'{DIP}.vcf.gz')
    bed = collections.defaultdict(list)
    for l in open(f'{DIP}.bed'):
        c, s, e = l.split()[:3]; bed[c].append((int(s), int(e)))
    for c in bed:
        bed[c].sort()
    inbed = lambda c, s, e: (lambda i: i >= 0 and bed[c][i][0] <= s and e <= bed[c][i][1])(bisect.bisect_right(bed[c], (s, 10 ** 12)) - 1)
    LS = loci()
    if ids:
        LS = [L for L in LS if L['truth_id'] in ids]
    rows, phase = [], []
    for n, L in enumerate(LS):
        t, ch, p0, vr, va = L['truth_id'], L['chrom'], int(L['vcf_pos']) - 1, L['vcf_ref'], L['vcf_alt']
        g, g0, lo, hi, capped = array(gfa, L)
        o = dict(L, arr_lo0=lo, arr_hi0=hi, arr_len=hi - lo, arr_capped=capped)
        # 1b imperfect tandem repeat (minisatellite / VNTR) around the truth, GRCh38 +-IVF bp
        w0 = max(g0, p0 - IVF); ivp, ivl, ivh = imperfect_repeat(g[w0 - g0:p0 + IVF - g0], p0 - w0)
        ivl, ivh = w0 + ivl, w0 + ivh
        beyond = max(0, lo - ivl) + max(0, ivh - hi) if ivp else 0
        o.update(ivntr_spacing=ivp, ivntr_lo0=ivl if ivp else '', ivntr_hi0=ivh if ivp else '', ivntr_beyond=beyond,
                 in_ivntr=ivp >= 7 and beyond >= IVOUT)
        segs, meta = {}, {}
        for h in HAP:
            hs = hits[h].get(t, [])
            q0 = hs[0]['q0'] if hs else 0; elo, ehi = lo - q0, hi - q0
            cov = [x for x in hs if x['qs'] <= elo - 200 and x['qe'] >= ehi + 200]
            o[f'{h}_n_cover'] = sum(x['tn'] == ch for x in cov); o[f'{h}_n_cover_other'] = sum(x['tn'] != ch for x in cov)
            tier, x, pts = place(hs, ch, elo, ehi)
            un = [y['nm'] for y in unloc[h].get(t, []) if y['qs'] <= elo - 200 and y['qe'] >= ehi + 200]
            o[f'{h}_unloc_better'] = bool(un) and (x is None or max(un) > x['nm'])
            if x is None:
                o[f'{h}_place_na'] = tier; continue
            o.update({f'{h}_place': tier, f'{h}_strand': x['st'], f'{h}_mapq': x['mapq'], f'{h}_tp': x['tp']})
            s0, s1 = max(0, pts[0] - CTX - (hi - lo)), min(x['tl'], pts[1] + CTX + (hi - lo))
            s = fas[h].fetch(ch, s0, s1).upper()
            segs[h] = s if x['st'] == '+' else rc(s); meta[h] = (ch, s0, s1, x['st'])
        # 2b paralog fix: a suspect hap copy (ragtag hit MAPQ < 20, a covering hit on another scaffold, a phased het SNV of
        # dipcall contradicted, or no ragtag placement) is replaced by the copy dipcall itself aligned there (its primary
        # contig alignment covering the array +-200 bp; highest MAPQ, then longest)
        near = list(dvcf.fetch(ch, max(0, lo - 2000), hi + 2000))
        ph0 = phase_check(g, g0, segs, near)
        for h in HAP:
            bad = sum(x['obs'][h] not in ('?', x['exp'][h]) for x in ph0 if h in x['obs'])
            sus = [k for k, v in (('unplaced', h not in segs), ('mapq<20', int(o.get(f'{h}_mapq', 99)) < 20),
                                  ('other_scaffold_hit', o[f'{h}_n_cover_other'] > 0), ('phase_mismatch', bad > 0)) if v]
            o.update({f'{h}_phase_mismatch_ragtag': bad if h in segs else '', f'{h}_suspect': ','.join(sus),
                      f'{h}_copy': 'ragtag' if h in segs else ''})
            if not sus:
                continue
            cov = [x for x in dpaf[h].get(ch, []) if x[0] <= lo - 200 and x[1] >= hi + 200]
            if not cov:                                         # no dipcall copy: a copy dipcall's phase contradicts is dropped
                o[f'{h}_dpaf'] = 'none'
                if bad > 0:
                    segs.pop(h, None); meta.pop(h, None); o.update({f'{h}_copy': 'dropped', f'{h}_place_na': 'dropped_paralog'})
                continue
            x = max(cov, key=lambda x: (x[6], x[4] - x[3]))
            cp = tproject(x, (lo + hi) // 2)
            s0, s1 = max(0, cp - CTX - (hi - lo)), min(raw[h].get_reference_length(x[2]), cp + CTX + (hi - lo))
            s = raw[h].fetch(x[2], s0, s1).upper()
            segs[h] = s if x[5] == '+' else rc(s); meta[h] = (x[2], s0, s1, x[5])
            o.update({f'{h}_copy': 'dipcall', f'{h}_dpaf': f'{x[2]}:{x[3]}-{x[4]}{x[5]} mapq {x[6]}'})
        ph = phase_check(g, g0, segs, near) if any(o[f'{h}_copy'] == 'dipcall' for h in HAP) else ph0
        for h in HAP:
            o[f'{h}_phase_mismatch'] = sum(x['obs'][h] not in ('?', x['exp'][h]) for x in ph if h in x['obs']) if h in segs else ''
        a, b, ok = anchors(g, g0, lo, hi, segs)
        o['anchor_wide'] = False
        if segs and len(ok) < len(segs):                       # a placed hap lacks the common anchors: wide retry, all haps
            w = anchors(g, g0, lo, hi, segs, SHIFT_WIDE, need_all=True)
            if w[1] is not None:
                (a, b, ok), o['anchor_wide'] = w, True
        o.update(anchor_a0=a if a is not None else '', anchor_b0=b if b is not None else '')
        R = A = None
        if a is not None and b is not None:
            R = g[a - g0:b - g0]; A = edit(R, a, [(p0, vr, va)])
            recs = [(r.pos - 1, r.ref.upper(), tuple((z or '').upper() for z in r.alts or ()), r.samples[0]['GT'],
                     not r.filter.keys() or 'PASS' in r.filter.keys()) for r in dvcf.fetch(ch, a, b)]
            ib = inbed(ch, a, b)
            o['len_R'] = len(R)
        seqs, ctx = {}, {}
        for h in HAP:
            S = None
            if h in ok:
                i, j = ok[h]; S = segs[h][i:j]; nm, s0, s1, st = meta[h]
                c0, c1 = (s0 + i, s0 + j) if st == '+' else (s1 - j, s1 - i)
                o[f'{h}_contig_span'] = f'{nm}:{c0}-{c1}{st}'
                ctx[h] = (segs[h][max(0, i - K - FL):i], segs[h][j:j + K + FL])
            ds, why = dip_hap(R, a, b, recs, GTI[h], ib) if R is not None else (None, 'no_anchors')
            call = lambda s: 'NA' if s is None else 'ALT' if s == A else 'REF' if s == R else 'other'
            o.update({f'{h}_asm_call': call(S), f'{h}_asm_len_change': len(S) - len(R) if S is not None else '',
                      f'{h}_asm_na': '' if S is not None else 'no_anchors' if R is None else o.get(f'{h}_place_na') or 'anchor',
                      f'{h}_dip_call': call(ds), f'{h}_dip_len_change': len(ds) - len(R) if ds is not None else '',
                      f'{h}_dip_na': why, f'{h}_asm_eq_dip': (S == ds) if S is not None and ds is not None else ''})
            if S is not None and ds is not None:                # GT-order check: the other GT index
                dx, _ = dip_hap(R, a, b, recs, 1 - GTI[h], ib)
                o[f'{h}_asm_eq_dip_swapped'] = S == dx if dx is not None else ''
            seqs[h] = S if S is not None else ds
            o[f'{h}_src'] = 'asm' if S is not None else 'dip' if ds is not None else 'NA'
            o[f'{h}_seq'] = seqs[h] or ''
            if h not in ctx and R is not None:
                ctx[h] = (g[a - g0 - K - FL:a - g0], g[b - g0:b - g0 + K + FL])
        o['R_seq'], o['A_seq'] = R or '', A or ''
        # 4 flank vs array, allele classes, distances
        cores = {}
        for h in HAP:
            S = seqs[h]
            if S is None:
                o[f'{h}_allele'] = 'NA'; continue
            dr = diff_runs(S, R)
            if dr is None:
                core, nf = S, 'unresolved'
            else:
                arr = [d for d in dr if d[0] <= hi - a and d[1] >= lo - a]
                core = edit(R, 0, [(r0, R[r0:r1], x) for r0, r1, x in arr]); nf = len(dr) - len(arr)
            cores[h] = core
            o.update({f'{h}_n_flank_diffs': nf, f'{h}_len_change': len(core) - len(R), f'{h}_dist_A': lev(core, A),
                      f'{h}_allele': 'ALT' if core == A else 'REF' if core == R else
                      f'germline_len{len(core) - len(R):+d}' if len(core) != len(R) else 'germline_seq',
                      f'{h}_total_len_change': len(S) - len(R), f'{h}_flank_len': len(S) - len(core)})
        # a hap whose array allele has exactly the ALT length but is not the ALT (INDEL truths): germline-like reading
        o['alt_len_haps'] = ','.join(h for h in HAP if h in cores and cores[h] != A and len(cores[h]) == len(A) != len(R))
        # 5 event
        if cores:
            hs = sorted(cores, key=lambda h: (o[f'{h}_dist_A'], cores[h] != R, h))
            e = hs[0]
            tie = len(hs) == 2 and o[f'{hs[0]}_dist_A'] == o[f'{hs[1]}_dist_A']
            o.update(event_hap=e, event_basis='both' if len(cores) == 2 else f'only_{e}',
                     event_tie='' if not tie else 'tie_same_seq' if cores[hs[0]] == cores[hs[1]] else
                     'tie_REF_chosen' if cores[e] == R else 'tie_other')
            o.update(event(cores[e], A, *ctx[e]))
            if len(hs) == 2:                                    # the other hap's event (when the closest is not parsimonious)
                ev2 = event(cores[hs[1]], A, *ctx[hs[1]])
                o.update(other_hap_kind=ev2['pf_kind'], other_hap_event=ev2.get('pf_event', ''), other_hap_label=ev2['pf_label'])
        else:
            o.update(event_basis='none', pf_kind='NA', pf_class='NA', pf_class_b5='NA', pf_label='NA')
        # 5b imperfect tandem repeat: a 'none' event inside one (reaching >= IVOUT bp beyond the exact-period array) is not a
        # non-repeat event; pf_class_b5 stays the b5 rule verbatim
        if o.get('pf_class') == 'none' and o['in_ivntr']:
            o.update(pf_class='imperfect_VNTR', pf_label='SNV in imperfect VNTR' if o['pf_kind'] == 'SNV' else 'imperfect VNTR')
        # 5c one hap only (the other NA): the missing hap may carry the ALT, so no confident class (kept in pf_label_one_hap)
        if o.get('event_basis', '').startswith('only_') and o['pf_kind'] != 'normal_carries_ALT':
            o.update(pf_label_one_hap=o['pf_label'], pf_label='single hap')
        # dipcall GT order: isolated PASS phased het SNVs within array +-2 kb, 21-mer REF / ALT context in each final hap stretch
        if all(h in segs for h in HAP):
            for x in ph:
                r, obs, exp = x['r'], x['obs'], x['exp']
                phase.append(dict(truth_id=t, pos=r.pos, ref=r.ref, alt=r.alts[0], gt='|'.join(map(str, x['gt'])),
                                  obs_hapY_hapX=f"{obs['hapY']}|{obs['hapX']}",
                                  verdict='undetermined' if '?' in obs.values() else 'consistent' if obs == exp else
                                  'swapped' if obs == {h: exp[o2] for h, o2 in (('hapX', 'hapY'), ('hapY', 'hapX'))} else 'other'))
        fl = [f'{h}:{k}' for h in HAP for k, bad in
              (('NA', o[f'{h}_src'] == 'NA'), ('dip_only', o[f'{h}_src'] == 'dip'), ('asm_ne_dip', o[f'{h}_asm_eq_dip'] is False),
               ('mapq<20', int(o.get(f'{h}_mapq', 99)) < 20), ('other_scaffold_hit', o.get(f'{h}_n_cover_other', 0) > 0),
               (o.get(f'{h}_place', ''), o.get(f'{h}_place', '') in ('cover0', 'bridge')),
               ('dipcall_copy', o[f'{h}_copy'] == 'dipcall'), ('dropped_paralog', o[f'{h}_copy'] == 'dropped'),
               ('flank_diffs>=3', isinstance(o.get(f'{h}_n_flank_diffs'), int) and o[f'{h}_n_flank_diffs'] >= 3),
               ('flank_len>=4', abs(o.get(f'{h}_flank_len', 0)) >= 4)) if bad]
        fl += ['wide_anchor'] if o['anchor_wide'] else []
        fl += [o['event_tie']] if o.get('event_tie') in ('tie_REF_chosen', 'tie_other') else []
        fl += ['other_hap_single_event'] if o.get('pf_kind') == 'complex' and o.get('other_hap_kind') in ('INDEL', 'SNV') else []
        o['flags'] = ','.join(fl)
        rows.append(o)
        if n % 250 == 0:
            print('loci', n, len(LS), flush=True)
    return rows, phase


COLS = ['truth_id', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind2', 'absorbed', 'arr_lo0', 'arr_hi0', 'arr_len', 'arr_capped',
        'ivntr_spacing', 'ivntr_lo0', 'ivntr_hi0', 'ivntr_beyond', 'in_ivntr', 'anchor_a0', 'anchor_b0', 'anchor_wide', 'len_R']
for _h in HAP:
    COLS += [f'{_h}_{c}' for c in ('place', 'place_na', 'unloc_better', 'strand', 'mapq', 'tp', 'n_cover', 'n_cover_other',
                                   'phase_mismatch_ragtag', 'suspect', 'copy', 'dpaf', 'phase_mismatch', 'contig_span', 'asm_call',
                                   'asm_len_change', 'asm_na', 'dip_call', 'dip_len_change', 'dip_na', 'asm_eq_dip',
                                   'asm_eq_dip_swapped', 'src', 'n_flank_diffs', 'allele', 'len_change', 'total_len_change',
                                   'flank_len', 'dist_A')]
COLS += ['alt_len_haps', 'event_hap', 'event_basis', 'event_tie', 'pf_kind', 'pf_event', 'pf_len_change', 'pf_unit', 'pf_tract',
         'pf_tract_b5', 'pf_tandem', 'pf_stretch', 'pf_class', 'pf_class_b5', 'pf_label', 'pf_label_one_hap', 'pf_units_changed',
         'pf_cpg_ti', 'pf_titv', 'other_hap_kind', 'other_hap_event', 'other_hap_label', 'flags', 'R_seq', 'A_seq',
         'hapX_seq', 'hapY_seq']


def show(ids):
    for r in csv.DictReader(open(f'{D}/normal_frame.tsv'), delimiter='\t'):
        if r['truth_id'] in ids:
            print(f"== {r['truth_id']} {r['chrom']}:{r['vcf_pos']} {r['vcf_ref']}>{r['vcf_alt']} array {r['arr_lo0']}-{r['arr_hi0']} "
                  f"anchors {r['anchor_a0']}-{r['anchor_b0']} absorbed {r['absorbed']}")
            for k in ('R_seq', 'A_seq', 'hapX_seq', 'hapY_seq'):
                print(f'  {k:9s}{r[k][:300]}')
            for h in HAP:
                print(f"  {h}: src {r[f'{h}_src']} asm {r[f'{h}_asm_call']} dip {r[f'{h}_dip_call']} eq {r[f'{h}_asm_eq_dip']} "
                      f"allele {r[f'{h}_allele']} flank {r[f'{h}_n_flank_diffs']} dist_A {r[f'{h}_dist_A']} na {r[f'{h}_asm_na']}/{r[f'{h}_dip_na']}")
            print(f"  event {r['event_hap']} {r['event_tie']} {r['pf_kind']} {r['pf_event']} len {r['pf_len_change']} u {r['pf_unit']} "
                  f"L {r['pf_tract']} tandem {r['pf_tandem']} stretch {r['pf_stretch']} class {r['pf_class']} (b5 {r['pf_class_b5']}) "
                  f"units {r['pf_units_changed']}")


def main():
    if sys.argv[1:2] == ['queries']:
        return queries()
    if sys.argv[1:2] == ['show']:
        return show(set(sys.argv[2].split(',')))
    ids = set(sys.argv[1].split(',')) if sys.argv[1:] else None
    rows, phase = run(ids)
    if ids:
        for r in rows:
            print('\n'.join(f'  {c}: {str(r.get(c, ""))[:300]}' for c in COLS)); print()
        print(phase)
        return
    with open(f'{D}/normal_frame.tsv', 'w') as w:
        w.write('\t'.join(COLS) + '\n')
        for r in rows:
            w.write('\t'.join(str(r.get(c, '')) for c in COLS) + '\n')
    with open(f'{D}/c3_dipcall_phase.tsv', 'w') as w:
        pc = ['truth_id', 'pos', 'ref', 'alt', 'gt', 'obs_hapY_hapX', 'verdict']
        w.write('\t'.join(pc) + '\n')
        for x in phase:
            w.write('\t'.join(str(x[c]) for c in pc) + '\n')
    C = collections.Counter
    print('loci', len(rows), C((r['kind2'], r['absorbed']) for r in rows), 'anchored', sum(r['anchor_b0'] != '' for r in rows),
          'capped', sum(r['arr_capped'] for r in rows))
    for h in HAP:
        print(h, 'place', C(r.get(f'{h}_place', '') for r in rows), 'place_na', C(r.get(f'{h}_place_na', '') for r in rows),
              'unloc_better', C(r[f'{h}_unloc_better'] for r in rows), 'mapq<20', sum(int(r.get(f'{h}_mapq') or 99) < 20 for r in rows))
        print(h, 'asm', C(r[f'{h}_asm_call'] for r in rows),
              'dip', C(r[f'{h}_dip_call'] for r in rows))
        print(h, 'asm_eq_dip', C(r[f'{h}_asm_eq_dip'] for r in rows), 'swapped', C(r.get(f'{h}_asm_eq_dip_swapped', '') for r in rows),
              'asm_na', C(r[f'{h}_asm_na'] for r in rows), 'dip_na', C(r[f'{h}_dip_na'] for r in rows))
        print(h, 'src', C(r[f'{h}_src'] for r in rows), 'allele', C(r[f'{h}_allele'] for r in rows).most_common(12),
              'flank', C(r.get(f'{h}_n_flank_diffs', '') for r in rows))
    for h in HAP:
        print(h, 'suspect', C(r[f'{h}_suspect'] for r in rows if r[f'{h}_suspect']).most_common(), 'copy', C(r[f'{h}_copy'] for r in rows),
              'dpaf none', sum(r.get(f'{h}_dpaf') == 'none' for r in rows),
              'phase mismatch ragtag -> final', sum(bool(r[f'{h}_phase_mismatch_ragtag']) for r in rows), sum(bool(r[f'{h}_phase_mismatch']) for r in rows),
              'flank_len>=4', sum(abs(r.get(f'{h}_flank_len', 0)) >= 4 for r in rows))
    print('anchor_wide', sum(r['anchor_wide'] for r in rows), 'in_ivntr', C((r['absorbed'], r['kind2']) for r in rows if r['in_ivntr']),
          'single hap', C((r['absorbed'], r['pf_label_one_hap']) for r in rows if r['pf_label'] == 'single hap'),
          'alt_len_haps (absorbed INDEL)', C(r['pf_label'] for r in rows if r['alt_len_haps'] and r['absorbed'] and r['kind2'] == 'INDEL'))
    vd = collections.defaultdict(set)
    for x in phase:
        vd[(x['pos'], x['ref'], x['alt'])].add(x['verdict'])
    print('phase rows (locus x SNV)', C(x['verdict'] for x in phase), 'distinct SNVs', len(vd),
          C(next(iter(s)) if len(s) == 1 else 'conflicting' for s in vd.values()))
    print('event', C((r.get('event_basis'), r.get('event_tie', '')) for r in rows))
    for ab in (True, False):
        for kd in ('INDEL', 'SNV'):
            sel = [r for r in rows if r['absorbed'] == ab and r['kind2'] == kd]
            if sel:
                print('absorbed' if ab else 'not_absorbed', kd, len(sel), C(r['pf_label'] for r in sel).most_common())
    print('flags', C(f for r in rows for f in r['flags'].split(',') if f).most_common(), 'loci with flags', sum(r['flags'] != '' for r in rows))
    print('complex with a single-event other hap', C((r['absorbed'], r['other_hap_label']) for r in rows if 'other_hap_single_event' in r['flags']))
    print('class differs from b5', sum(r.get('pf_class') != r.get('pf_class_b5') for r in rows),
          C((r.get('pf_class_b5'), r.get('pf_class')) for r in rows if r.get('pf_class') != r.get('pf_class_b5')))


if __name__ == '__main__':
    main()
