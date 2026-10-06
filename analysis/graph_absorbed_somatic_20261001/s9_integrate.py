"""Step 9: integration. Per truth allele of loci.tsv: the somatic graph elements the ALT reads walk, who else walks
them (HPRC haplotypes, the HG008-N assembly haplotypes), and a category: the patient's own germline allele already in
the graph (germline filtering) vs an allele / nodes of other HPRC haplotypes that hide the somatic mutation.

Inputs ($D): variant_set / loci (s0, s1); graph_paths / graph_elements (s2); read_paths_P.tsv.gz + read_path_summary_P.tsv
(s3; every platform whose summary file exists is used, the others count as without read data, so a re-run picks up PacBio /
ONT once merged); hprc_local_paths.tsv.gz / hprc_haplotypes / hprc_window_coverage (s4); hprc_vcf_alleles /
hprc_vcf_haplotypes.tsv.gz / hprc_sample_metadata (s5, s5b); asm_seq_status / asm_seq_summary / asm_local_paths /
dipcall_relation (s6-s8); array_alleles.tsv (s8b: whole-array alleles of HG008-N hap1 / hap2 and HG008-T contigs, required);
the HG008-N dipcall VCF of s2 (PASS alleles inside each window, nearest 6 to the truth = the s2 rule, with GT);
analysis/tensor_recall_20260930/pon/truth_pon_HG008T.tsv (rule, gnomAD / CoLoRSdb AF; all truths, also the non-absorbed);
GIAB v0.2 nogermlineoverlap BED (scope). Graph helpers: s2_graph_paths.py imported (classify.py ref / oseq / seq; s2
walk(), elements(), apply(), lev()). Fixes after the four audits of 2026-10-01: $D/audit/decisions.md.
Method
 0 Read walks (s3, no-edit spanning ALT reads, cut between GRCh38 anchor nodes): spelled over the walk's GRCh38 anchor span
   and classed alt (= GRCh38 + truth) / alt+germline (+ a subset of the 6 nearest PASS dipcall alleles, or + all PASS
   records of one normal hap inside the span) / ref / germline_only (a subset / one whole normal hap) / other (also when the
   truth overlaps a germline record of the same hap: no edit-wise composition). Truth-spelling reads = alt or alt+germline (s3 'alt_like' also admits germline other-length
   STR reads). truth_read_support = yes / no_truth_spelling_read (read data on a perfect platform, none spells the truth) /
   no_perfect_platform_read_data.
 1 Somatic elements (s2 element rules on any walk: 'B:x:y:>n..' = run of non-GRCh38 nodes between GRCh38 nodes x, y;
   'S:x:y' = consecutive GRCh38 nodes not adjacent in GRCh38; ref() of classify.py):
     exact          every element of the chosen exact path. Several paths (s2 walk() re-run with the same caps, count and
                    primary checked against graph_paths): the one whose set is contained in the most pooled truth-spelling
                    reads (all platforms with read data); ties: the first in s2 order (the primary when tied); no truth
                    read: the s2 primary
     with_germline  the same choice, phase-consistent paths first (every germline allele of the path on the truth's event
                    hap by dipcall GT; none -> other_ambiguous), minus the elements whose sequence alone (the element
                    applied alone to GRCh38 over the window) = GRCh38 + a non-empty subset of the path's germline alleles
     closest        the commonest walk of the pooled truth-spelling reads (element_source reads_truth_spelling); without
                    such reads the commonest pooled walk (element_source reads) or the s2 closest primary
                    (graph_closest_no_reads) is reported but the locus is other_ambiguous (the elements do not encode the
                    truth; closest_kind = read_allele_nearer_truth / read_allele_not_truth by Levenshtein of the kept set
                    applied to GRCh38 vs truth(+germline) and GRCh38(+germline) haplotypes); minus germline-alone elements
   contains(walk, element): branch = x, the node run, y consecutive in order; skip = x then y consecutive (GRCh38-forward
   walks). A read tests only the elements between its anchor nodes; it contains the set = every testable element
   contained (>= 1 testable). Per platform on all good reads and on truth-spelling reads: n, n_testable, n_contain; the
   majority walk (truth-spelling reads if any) vs the graph ALT path elements inside that walk's span.
 2 HPRC (88 haplotypes = 44 x 2; CHM13 separately; GRCh38 never): s4 complete traversals only (a haplotype with several
   complete traversals counts if any of them does). Element level: the traversal contains the element (locus: min over
   the elements, elem_freq_min); set level: all elements of the chosen path (set_carriers) / of any candidate ALT path
   (set_carriers_any_path, representation-independent); allele level: bracket = the last GRCh38 node of the ALT path
   starting before lo and the first after it ending at / after hi that the haplotype walk also visits (lo, hi = the s8b
   whole array united with the event window and the elements' spans, clipped to the graph window); carrier = same
   spelled sequence from L to R as the ALT path (allele) or as the truth haplotype (truth_seq); exact allele = truth_seq,
   or the path allele for with_germline (truth + the patient's germline alleles). Frequencies over the complete
   haplotypes (an upper-biased conditional frequency: partial haplotypes carry less, audit) and over all 88 (lower bound).
   hprc_support = exact_allele (>= 1 exact-allele carrier) / element_set_other_allele (>= 1 haplotype walks a full ALT
   element set but spells another allele over the array) / recombinant_pieces (every element has carriers, no haplotype
   has a full set) / none. RARE_AF = 0.20 on exact_allele_freq (broad / rare = near the d9 floor; freq_class_agree = the
   element / set-any-path / exact / full-VCF frequencies fall on the same side of 0.20).
   VCF cross-check: carriers_any of the d9 and full VCFs vs allele-level carriers, on the haplotypes complete in the GFA.
   d9 floor: haplotypes (HPRC + CHM13 + GRCh38) whose s4 rows (complete or partial paths: a lower bound) visit each branch
   node of a somatic element.
 3 HG008-N hap1 / hap2: array call (s8b): the hap's sequence between common unique anchors around the whole tandem array
   (assembly; dipcall reconstruction where the assembly is NA) == GRCh38 + truth (or + the hap's own flank records) ->
   ALT, else no; assembly and dipcall disagreeing on ALT -> unresolved (conflict); both NA -> the narrow s6/s7 call when
   it is 'no', else unresolved (array_NA). Narrow call (kept as columns): minimap2 ALT (accepted only when the hap's own
   hap_diffs give alt_hap), giraffe ALT on the elements, ALT_plus_other (one or two hap_diffs give alt_hap); spans_agree
   False -> unresolved. Event hap = truth INFO HG008Nv62SOMATICVARIANT contig; T_info_class (s8b) = the INFO event applied
   to the event hap's array sequence vs ALT / the other hap. Pattern: other_hap_only / event_hap_only / both /
   no_event_hap (ALT on a hap), other_hap_only:patient_frame (no hap == GRCh38 + truth, but the tumor event-hap allele by
   the INFO event == the other hap's allele: they differ only by germline differences shared by both haps), unresolved,
   neither. Allele class per hap from the array length change: ALT / REF / STR_other_length / same_length_other_seq /
   substitution / unresolved; locus: carries_ALT / both_REF / SNV_site_deleted_on_event_hap (s6 site_call) /
   germline_STR_other_length / other / unresolved. Tumor (s8b): tumor_verdict, tumor_any_ALT_array.
 4 Category, first rule that applies -> other_ambiguous with sub_reason:
     closest without truth-spelling reads (closest_no_read_element / closest:read_allele_*); no kept element
     (reads_walk_only_germline_branches / no_somatic_element); with_germline:germline_allele_off_event_hap;
     HPRC_no_complete_haplotype; evidence_conflict:ALT_on_both_normal_haps; evidence_conflict:ALT_on_normal_event_hap
     (the GRCh38 ALT is the normal event hap's own allele: truth representation conflict);
   HG008N_present_broad / _rare  ALT on the other hap / no event hap / patient frame, exact_allele_freq >= / < RARE_AF
                                 (unless T_info_class novel -> other_ambiguous evidence_conflict:germline_ALT_vs_truth_INFO_event_novel)
   then HG008N_hap_unresolved[:conflict / :array_NA]; no_HPRC_carrier[:CHM13_carries] (hprc_support none);
   HG008N_absent_HPRC_other      sub_class = hprc_support
Outputs ($D): per_variant.tsv, per_node.tsv, hprc_membership.tsv.gz, summary.json, tables.md, spot_checks.txt.
Run: python s9_integrate.py [n_workers] [truth_id,...] (default SLURM_CPUS_PER_TASK or 4; with ids: printout only, no files;
sbatch tmp/graph_absorbed_somatic_20261001/s9.sbatch; the 25 GB node sqlite on beegfs is IO-bound when cold, hence the pool).
Assumptions: PacBio / ONT read data may be absent (recorded per platform); HPRC partial traversals are outside the
denominators; allele level for closest = the read-walked graph allele; dipcall GT order = hap1|hap2.
"""
import collections, csv, gzip, importlib.util, itertools, json, math, os, random, re, statistics, sys
from bisect import bisect_left, bisect_right
from multiprocessing import Pool
import pysam
D = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
A = '/scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_20261001'
PON = '/scratch/jshen/Github/Pansoma/analysis/tensor_recall_20260930/pon/truth_pon_HG008T.tsv'
NOG = '/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_nogermlineoverlap.bed'
PLATS, RARE_AF, SUPER = ('PacBio', 'ONT', 'Illumina'), 0.20, ('AFR', 'AMR', 'EAS', 'SAS')
CATS = ('HG008N_present_broad', 'HG008N_present_rare', 'HG008N_absent_HPRC_other', 'other_ambiguous')
TRUTHC = ('alt', 'alt+germline')
csv.field_size_limit(sys.maxsize)
_s = importlib.util.spec_from_file_location('s2', f'{A}/s2_graph_paths.py'); S2 = importlib.util.module_from_spec(_s); _s.loader.exec_module(S2)
ref, oseq, seq, tag = S2.ref, S2.oseq, S2.ns['seq'], S2.tag
rd = lambda f: list(csv.DictReader((gzip.open(f, 'rt') if f.endswith('.gz') else open(f)), delimiter='\t'))
toks = lambda w: tuple((int(n), '+' if o == '>' else '-') for o, n in re.findall(r'([<>])(\d+)', w or ''))
fcls = lambda f: 'NA' if f is None or f == '' else '<0.2' if f < RARE_AF else '0.2-0.5' if f < 0.5 else '>=0.5'
_sp = {}


def spell(w):
    if w not in _sp:
        _sp[w] = ''.join(oseq(n, o) for n, o in w).upper()
    return _sp[w]


def prep(e):
    """s2 element + token pattern x, run, y and the GRCh38 start of x / end of y."""
    e = dict(e, x=int(e['x']), y=int(e['y']))
    e['pat'] = ((e['x'], '+'),) + toks(e['nodes']) + ((e['y'], '+'),)
    e['xs'], e['ye'] = ref(e['x'])[1], ref(e['y'])[2]
    return e


def has(w, e):
    p = e['pat']; k = len(p)
    return any(w[i:i + k] == p for i in range(len(w) - k + 1) if w[i] == p[0])


def testable(w, e):
    a, z = (ref(w[0][0]) if w[0][1] == '+' else None), (ref(w[-1][0]) if w[-1][1] == '+' else None)
    return bool(a and z) and a[1] <= e['xs'] and e['ye'] <= z[2]


def contains_set(w, els):
    t = [e for e in els if testable(w, e)]
    return None if not t else all(has(w, e) for e in t)


def alone(L, e):
    b = int(L['win_start0']); s, t = e['start0'] - b, e['end0'] - b
    return L['ref_hap'][:s] + e['seq'] + L['ref_hap'][t:] if s <= t else None


def som(L):
    return int(L['vcf_pos']) - 1, L['vcf_ref'].upper(), L['vcf_alt'].upper()


def gseqs(L, alleles):
    """window sequences with a non-empty subset of the germline alleles applied (without / with the truth allele)."""
    b, g, tg = int(L['win_start0']), set(), set()
    for k in range(1, len(alleles) + 1):
        for c in itertools.combinations(alleles, k):
            g.add(S2.apply(L['ref_hap'], b, list(c))); tg.add(S2.apply(L['ref_hap'], b, [som(L)] + list(c)))
    return g - {None}, tg - {None}


def rclass(L, w, cache):
    """what a read walk spells over its GRCh38 anchor span: alt / alt+germline / ref / germline_only / other."""
    if w not in cache:
        a = ref(w[0][0]) if w and w[0][1] == '+' else None; z = ref(w[-1][0]) if w and w[-1][1] == '+' else None
        b, e = int(L['win_start0']), int(L['win_end0'])
        if not (a and z) or a[1] < b or z[2] > e or z[2] < a[1]:
            cache[w] = 'other'; return cache[w]
        lo, hi = a[1], z[2] + 1
        r, s, v = L['ref_hap'][lo - b:hi - b], spell(w), som(L)
        gs = [x for x in L['germ'] if lo <= x[0] and x[0] + len(x[1]) <= hi]
        hk = [[x for x in L['germ_all'] if L['germ_gt'][x][k] and lo <= x[0] and x[0] + len(x[1]) <= hi] for k in (0, 1)]
        inside = lo <= v[0] and v[0] + len(v[1]) <= hi
        c = 'alt' if inside and s == S2.apply(r, lo, [v]) else 'ref' if s == r else 'other'
        if c == 'other':                             # each whole normal hap over the span (all its PASS records), +- truth
            c = 'alt+germline' if inside and s in {S2.apply(r, lo, [v] + h) for h in hk if h} else \
                'germline_only' if s in {S2.apply(r, lo, h) for h in hk if h} else 'other'
        for k in range(1, len(gs) + 1) if c in ('other', 'germline_only') else ():
            for cb in itertools.combinations(gs, k):
                if inside and s == S2.apply(r, lo, [v] + list(cb)):
                    c = 'alt+germline'; break
                if c == 'other' and s == S2.apply(r, lo, list(cb)):
                    c = 'germline_only'
            if c == 'alt+germline':
                break
        cache[w] = c
    return cache[w]


def bed(path):
    b = collections.defaultdict(list)
    for l in open(path):
        c, s, e = l.split()[:3]; b[c].append((int(s), int(e)))
    for c in b:
        b[c].sort()
    return b


def inbed(b, c, s, e):
    i = bisect_right(b.get(c, []), (s, 10 ** 12)) - 1
    return i >= 0 and b[c][i][0] <= s and e <= b[c][i][1]


def germline(loci):
    """PASS HG008-N dipcall alleles inside each window, the 6 nearest to the truth (s2 run() rule, same order), with
    germ_gt[allele] = (on hap1, on hap2) from the phased GT; germ_all = every such allele with a diploid GT."""
    iv = collections.defaultdict(list)
    for L in loci:
        L['germ'], L['germ_gt'] = [], {}; iv[L['chrom']].append((int(L['win_start0']), int(L['win_end0']), L))
    for c in iv:
        iv[c].sort(key=lambda z: z[0])
    st = {c: [z[0] for z in v] for c, v in iv.items()}
    spn = max(b - a for v in iv.values() for a, b, _ in v)
    for rec in pysam.VariantFile(S2.GERM_VCF):
        v = iv.get(rec.chrom)
        if not v or (rec.filter.keys() and 'PASS' not in rec.filter.keys()):
            continue
        p0, a, gt = rec.pos - 1, rec.ref.upper(), rec.samples[0]['GT'] or ()
        for s, e, L in v[bisect_left(st[rec.chrom], p0 - spn):bisect_left(st[rec.chrom], p0 + 1)]:
            if s <= p0 and p0 + len(a) <= e + 1:
                for i, b in enumerate(rec.alts or (), 1):
                    if b and b[0] not in '<*':
                        al = (p0, a, b.upper()); L['germ'].append(al); L['germ_gt'][al] = tuple(x == i for x in gt)
    for L in loci:
        p0 = int(L['vcf_pos']) - 1
        L['germ_all'] = [x for x in L['germ'] if len(L['germ_gt'][x]) == 2]
        L['germ'] = sorted(L['germ'], key=lambda v: abs(v[0] - p0))[:6]


def phase_ok(L, used, evh):
    """True / False: every germline allele of the path lies on the event hap (dipcall GT); None: no event hap / no GT."""
    if not used or evh not in ('hap1', 'hap2'):
        return None
    g = [L['germ_gt'].get(tuple(al)) for al in used]
    if any(x is None or len(x) < 2 for x in g):
        return None
    return all(x[int(evh[-1]) - 1] for x in g)


def candidates(L, g):
    """[(path tokens, germline alleles of the path)] for exact / with_germline, the s2 primary first."""
    prim = toks(g['primary_path'])
    used = [(int(p) - 1, *ab.split('>')) for p, ab in (x.split(':') for x in g['germline_used'].split(';'))] if g['germline_used'] else []
    if int(g['n_alt_paths']) <= 1:
        return [(prim, used)]
    LL = dict(L, first=int(L['first_node']), last=int(L['last_node']), lo=int(L['win_start0']), hi=int(L['win_end0']))
    lab = {L['alt_hap']: []} if g['match'] == 'exact' else {}
    for c in ([] if g['match'] == 'exact' else itertools.combinations(L['germ'], len(used))):
        lab.setdefault(S2.apply(L['ref_hap'], LL['lo'], [som(L)] + list(c)), list(c))
    lab.pop(None, None)
    hits = S2.walk(LL, sorted(lab))[0]
    paths = [(tuple(p), lab[t]) for p, t in hits]
    assert len(paths) == int(g['n_alt_paths']) and prim in [p for p, _ in paths], L['truth_id']
    return [x for x in paths if x[0] == prim] + [x for x in paths if x[0] != prim]


def load():
    X = dict(V={r['truth_id']: r for r in rd(f'{D}/variant_set.tsv')}, L=rd(f'{D}/loci.tsv'),
             G={r['truth_id']: r for r in rd(f'{D}/graph_paths.tsv')}, PONALL=rd(PON),
             COV={r['truth_id']: r for r in rd(f'{D}/hprc_window_coverage.tsv')}, ARR={r['truth_id']: r for r in rd(f'{D}/array_alleles.tsv')},
             SUM={r['truth_id']: r for r in rd(f'{D}/asm_seq_summary.tsv')}, DIP={r['truth_id']: r for r in rd(f'{D}/dipcall_relation.tsv')},
             AS={(r['truth_id'], r['assembly'], r['contig']): r for r in rd(f'{D}/asm_seq_status.tsv')},
             AP={(r['truth_id'], r['assembly'], r['contig']): r for r in rd(f'{D}/asm_local_paths.tsv')},
             VA={(r['truth_id'], r['source']): r for r in rd(f'{D}/hprc_vcf_alleles.tsv')}, VH=collections.defaultdict(dict),
             META={r['sample']: r for r in rd(f'{D}/hprc_sample_metadata.tsv')}, HP=collections.defaultdict(list), R={}, SM={}, NOG=bed(NOG))
    X['PON'] = {r['truth_id']: r for r in X['PONALL']}
    assert [L['truth_id'] for L in X['L']] == list(X['ARR']), 'array_alleles.tsv (s8b) must match loci.tsv'
    ev = collections.Counter(L['asm_event'] for L in X['L'] if L['asm_event'])
    for L in X['L']:
        L['n_truths_same_asm_event'] = ev[L['asm_event']] if L['asm_event'] else ''
    for r in rd(f'{D}/hprc_vcf_haplotypes.tsv.gz'):
        X['VH'][(r['truth_id'], r['source'])][f"{r['sample']}#{r['hap']}"] = r['allele_index']
    for r in rd(f'{D}/hprc_local_paths.tsv.gz'):
        X['HP'][r['truth_id']].append((f"{r['sample']}#{r['hap']}", r['status'], r['local_path'], r['partial_path']))
    X['HAPS'] = [f"{r['sample']}#{r['hap']}" for r in rd(f'{D}/hprc_haplotypes.tsv') if r['reference'] == 'False']
    assert len(X['HAPS']) == 88
    for p in PLATS:
        if not os.path.exists(f'{D}/read_path_summary_{p}.tsv'):
            continue
        X['SM'][p] = {r['truth_id']: r for r in rd(f'{D}/read_path_summary_{p}.tsv')}
        X['R'][p] = collections.defaultdict(list)
        for r in rd(f'{D}/read_paths_{p}.tsv.gz'):
            if r['call'].startswith('alt') and r['has_edit_in_window'] == 'False' and r['spans_event'] == 'True':
                X['R'][p][r['truth_id']].append(r['local_path'])
        for t, s in X['SM'][p].items():
            assert len(X['R'][p].get(t, [])) == int(s['n_alt_no_edit_spanning']), (p, t)
    germline(X['L'])
    return X


def bracket(w, P, gi, lo, hi):
    """(walk seq, ALT path seq, L, R), L..R inclusive: L = the last GRCh38 node of P starting before lo, R = the first after
    it ending at / after hi, both visited by the walk (the s3 cut rule; shared nodes spell the same, so they may reach
    into [lo, hi)); None if there is no such pair."""
    first = {}
    for i, x in enumerate(w):
        first.setdefault(x, i)
    ls = [i for i, r in gi if r[1] < lo and P[i] in first]
    if not ls:
        return None
    iL = max(ls); jL = first[P[iL]]
    for i, r in gi:
        if i > iL and r[2] >= hi and P[i] in first:
            jR = next((j for j in range(jL + 1, len(w)) if w[j] == P[i]), None)
            if jR is not None:
                return spell(w[jL:jR + 1]), spell(P[iL:i + 1]), P[iL][0], P[i][0]
    return None


def diffs(L, s):
    """the hap's event-window differences from GRCh38 (s6 hap_diffs: 'pos1:XG>A' base, 'pos1:DSEQ' deleted from pos1,
    'pos1:ISEQ' inserted after pos1) as (p0, ref, alt) edits checked against GRCh38."""
    b, e = int(L['win_start0']), []
    for d in filter(None, (s or {}).get('hap_diffs', '').split(';')):
        p, v = d.split(':', 1); p0 = int(p) - 1
        x = tuple(v[1:].split('>')) if v[0] == 'X' else (v[1:], '') if v[0] == 'D' else ('', v[1:])
        p0 += v[0] == 'I'
        if L['ref_hap'][p0 - b:p0 - b + len(x[0])] == x[0] and b <= p0 <= b + len(L['ref_hap']):
            e.append((p0, *x))
    return e


def truth_in_diffs(L, s):
    """one or two of the hap's event-window differences applied alone give alt_hap."""
    b, e = int(L['win_start0']), diffs(L, s)
    return any(S2.apply(L['ref_hap'], b, list(c)) == L['alt_hap'] for k in (1, 2) for c in itertools.combinations(e, k))


def carries(s, p, wall, L):
    """narrow (s6 event window) call of a HG008-N hap: (carries_alt, flag)."""
    mm, gi = (s or {}).get('status', 'none'), (p or {}).get('asm_event_allele', '')
    if s and s.get('spans_agree') == 'False':
        return 'unresolved', 'spans_disagree'
    cons = lambda: S2.apply(L['ref_hap'], int(L['win_start0']), diffs(L, s)) == L['alt_hap'] or truth_in_diffs(L, s)
    if mm == 'ALT':
        if not cons():
            return 'unresolved', 'conflict:mm2_ALT_vs_hap_diffs'
        return 'ALT', 'mm2_ALT_giraffe_REF' if gi == 'REF' else '' if gi == 'ALT' else f'mm2_ALT_giraffe_{gi or "none"}'
    if gi == 'ALT' and wall:
        if mm == 'REF':
            return 'unresolved', 'conflict:mm2_REF_giraffe_ALT'
        return ('ALT', f'giraffe_only:mm2_{mm}') if mm in ('ambiguous', 'none') or cons() else ('unresolved', 'conflict:giraffe_ALT_vs_hap_diffs')
    if gi == 'ALT':
        return 'unresolved', f'conflict:giraffe_ALT_off_elements:mm2_{mm}'
    if mm == 'other' and truth_in_diffs(L, s):
        return 'ALT_plus_other', 'truth_allele_in_hap_diffs'
    if mm in ('REF', 'other'):
        return 'no', ''
    return ('no', f'mm2_{mm}_giraffe_REF') if gi == 'REF' else ('unresolved', f'mm2_{mm}_giraffe_{gi or "none"}')


def hap_call(ar, h, narrow):
    """final HG008-N hap call from the s8b array alleles (assembly, dipcall), narrow call only to confirm 'no'."""
    x, y = ar[f'{h}_asm_call'], ar[f'{h}_dip_call']
    isalt = lambda c: c in ('ALT', 'ALT_plus_flank')
    if x != 'NA' and y != 'NA' and isalt(x) != isalt(y):
        return 'unresolved', f'conflict:array_asm_{x}_vs_dipcall_{y}'
    c, src = (x, 'asm') if x != 'NA' else (y, 'dipcall')
    if c != 'NA':
        return ('ALT' if isalt(c) else 'no'), f'array_{src}:{c}'
    if narrow == 'no':
        return 'no', 'array_NA:narrow_no'
    return 'unresolved', f'array_NA:narrow_{narrow}'


def allele_class(L, s, c):
    """narrow allele class (s6 event window), used when the array has no sequence for the hap."""
    if c in ('ALT', 'ALT_plus_other'):
        return 'ALT'
    if c == 'unresolved' or not s:
        return 'unresolved'
    if s['status'] in ('REF', 'ambiguous', 'none'):
        return 'REF'
    n, tn = int(s['hap_len_change'] or 0), len(L['vcf_alt']) - len(L['vcf_ref'])
    return 'STR_other_length' if n not in (0, tn) else 'same_length_other_seq' if n == tn and n != 0 else 'substitution'


def event_kind(e):
    m = re.match(r'.*:\d+-([ACGTN]*)-([ACGTN]*)$', (e or '').upper())
    if not m:
        return 'none'
    a, b = m.groups()
    return 'SNV' if len(a) == len(b) == 1 else 'MNV' if len(a) == len(b) else 'INS' if b.startswith(a) else 'DEL' if a.startswith(b) else 'complex'


def locus(X, L):
    t, ev_lo, ev_hi = L['truth_id'], int(L['ev_lo0']), int(L['ev_hi0'])
    g, v, ar = X['G'][t], X['V'][t], X['ARR'][t]
    b0 = int(L['win_start0'])
    out = dict(v, match=g['match'], n_alt_paths=g['n_alt_paths'], germline_used=g['germline_used'], event_net=g['event_net'], truth_net=g['truth_net'])
    RW = {p: [toks(w) for w in X['R'][p].get(t, [])] for p in X['R']}
    rc = {}
    pooled = [w for p in RW for w in RW[p]]
    truth_pool = [w for w in pooled if rclass(L, w, rc) in TRUTHC]
    gs6 = gseqs(L, L['germ'])
    evh = ar['event_hap']
    # ---- 1 somatic elements
    closest_kind, phase = '', ''
    if g['match'] in ('exact', 'with_germline'):
        sc = []
        for i, (P, used) in enumerate(candidates(L, g)):
            els = [prep(e) for e in S2.elements(P, ev_lo, ev_hi)]
            gu = gseqs(L, used)[0] if used else set()
            keep = [e for e in els if alone(L, e) not in gu]
            ph = phase_ok(L, used, evh)
            sc.append((ph is not False, sum(contains_set(w, keep) is True for w in truth_pool), -i, P, keep, els, used, ph))
        best = max(sc, key=lambda z: z[:3])
        ph_ok, n_best, mi, P, keep, els, used, ph = best
        ties = sum(z[:2] == best[:2] for z in sc)
        choice = 'single_path' if len(sc) == 1 else ('truth_read_supported' if n_best else 'no_truth_read_support') + \
            (':primary' if mi == 0 else ':alternative') + (':tie' if ties > 1 and n_best else '') + (':phase' if not sc[0][0] and ph_ok else '')
        phase = '' if g['match'] == 'exact' else 'no_event_hap_or_GT' if ph is None else 'consistent' if ph else 'off_event_hap'
        src, cmp_els = 'graph', els
        alts = [z[4] for z in sc if z[4]]
        out.update(n_candidate_paths=len(sc), read_support_of_paths=';'.join(str(z[1]) for z in sc))
        dropped = [e for e in els if e not in keep]
    else:
        base = truth_pool or pooled
        maj = collections.Counter(base).most_common(1)
        P = maj[0][0] if maj else toks(g['primary_path'])
        src = 'reads_truth_spelling' if truth_pool else 'reads' if pooled else 'graph_closest_no_reads'
        choice = f"{'truth_read' if truth_pool else 'pooled'}_majority:{maj[0][1]}/{len(base)}" if maj else 'closest_primary'
        els = [prep(e) for e in S2.elements(P, ev_lo, ev_hi)]
        keep = [e for e in els if alone(L, e) not in gs6[0]]
        dropped = [e for e in els if e not in keep]
        cmp_els = [prep(e) for e in S2.elements(toks(g['primary_path']), ev_lo, ev_hi)]
        alts = [keep] if keep else []
        sp = S2.apply(L['ref_hap'], b0, [(e['start0'], L['ref_hap'][e['start0'] - b0:e['end0'] - b0], e['seq']) for e in keep])
        if truth_pool:
            closest_kind = 'truth_spelling_reads'
        elif sp is not None:
            dt = min(S2.lev(sp, x) for x in {L['alt_hap']} | gs6[1]); dn = min(S2.lev(sp, x) for x in {L['ref_hap']} | gs6[0])
            closest_kind = ('' if pooled else 's2_primary:') + ('read_allele_nearer_truth' if dt < dn else 'read_allele_not_truth')
            out.update(closest_lev_truth=dt, closest_lev_nontruth=dn)
        out.update(n_candidate_paths=g['n_alt_paths'], read_support_of_paths='')
    eclass = lambda e: (lambda s: 'truth' if s == L['alt_hap'] else 'germline' if s in gs6[0] else 'truth+germline' if s in gs6[1] else 'partial')(alone(L, e))
    out.update(element_source=src, path_choice=choice, closest_kind=closest_kind, germline_phase=phase, alt_path=tag(P),
               somatic_elements=';'.join(e['id'] for e in keep), n_somatic_elements=len(keep),
               element_subtypes=';'.join(e['subtype'] for e in keep), element_alone=';'.join(eclass(e) for e in keep),
               element_in_event=';'.join(str(e['in_event']) for e in keep), dropped_germline_elements=';'.join(e['id'] for e in dropped))
    # ---- read support per platform (all good reads and truth-spelling reads)
    node_reads, perfect_data, truth_any = {e['id']: {} for e in keep}, [], []
    for p in PLATS:
        if p not in X['R'] or t not in X['SM'][p]:
            out[f'{p}_read_data'] = False; continue
        ws = RW[p]; cl = [rclass(L, w, rc) for w in ws]; tr = [w for w, c in zip(ws, cl) if c in TRUTHC]; cc = collections.Counter(cl)
        res, rt = [contains_set(w, keep) for w in ws], [contains_set(w, keep) for w in tr]
        frac = lambda z: round(sum(r is True for r in z) / sum(r is not None for r in z), 3) if any(r is not None for r in z) else ''
        mj, ma = collections.Counter(tr or ws).most_common(1), collections.Counter(ws).most_common(1)
        if mj:
            me = {e['id'] for e in S2.elements(mj[0][0], ev_lo, ev_hi)}
            ge = {e['id'] for e in cmp_els if testable(mj[0][0], prep(e))}
        out.update({f'{p}_read_data': True, f'{p}_n_good': len(ws), f'{p}_n_truth_reads': len(tr), f'{p}_n_alt_reads': cc['alt'],
                    f'{p}_n_alt_germline_reads': cc['alt+germline'], f'{p}_n_germline_only_reads': cc['germline_only'],
                    f'{p}_n_ref_reads': cc['ref'], f'{p}_n_other_reads': cc['other'],
                    f'{p}_n_testable': sum(r is not None for r in res), f'{p}_n_contain_set': sum(r is True for r in res), f'{p}_frac_contain_set': frac(res),
                    f'{p}_truth_n_testable': sum(r is not None for r in rt), f'{p}_truth_n_contain_set': sum(r is True for r in rt),
                    f'{p}_truth_frac_contain_set': frac(rt), f'{p}_majority_from': ('truth_reads' if tr else 'all_reads') if mj else '',
                    f'{p}_majority_walk': tag(mj[0][0]) if mj else '', f'{p}_majority_fraction': round(mj[0][1] / len(tr or ws), 3) if mj else '',
                    f'{p}_majority_class': rclass(L, mj[0][0], rc) if mj else '', f'{p}_majority_vs_graph': ('same' if me == ge else 'different') if mj else '',
                    f'{p}_all_reads_majority_class': rclass(L, ma[0][0], rc) if ma else '',
                    f'{p}_majority_seq': X['SM'][p][t]['majority_seq']})
        if v[f'{p}_perfect'] == 'True':
            perfect_data.append(p)
            if tr:
                truth_any.append(p)
        for e in keep:
            te, tt = [w for w in ws if testable(w, e)], [w for w in tr if testable(w, e)]
            node_reads[e['id']][p] = (len(ws), len(te), sum(has(w, e) for w in te), len(tr), len(tt), sum(has(w, e) for w in tt))
    out.update(truth_read_support='yes' if truth_any else 'no_truth_spelling_read' if perfect_data else 'no_perfect_platform_read_data',
               truth_read_platforms='+'.join(truth_any), perfect_platforms_with_reads='+'.join(perfect_data))
    # ---- 2 HPRC membership
    sp_ = [z for e in keep for z in (e['start0'], e['end0'])]
    w_lo, w_hi = b0, int(L['win_end0'])
    a_lo, a_hi = int(ar['arr_lo0']), int(ar['arr_hi0'])
    lo, hi = min([ev_lo, max(a_lo, w_lo + 1)] + sp_), max([ev_hi, min(a_hi, w_hi)] + sp_)
    out.update(allele_lo0=lo, allele_hi0=hi, allele_window_clipped=a_lo <= w_lo or a_hi > w_hi)
    gi = [(i, ref(x[0])) for i, x in enumerate(P) if x[1] == '+' and ref(x[0])]
    comp, present, anyw = collections.defaultdict(list), set(), collections.defaultdict(list)
    for k, stt, lp, pp in X['HP'].get(t, []):
        present.add(k); lp, pp = toks(lp), toks(pp)
        if stt == 'complete':
            comp[k].append(lp)
        if lp or pp:
            anyw[k].append(lp or pp)
    anyn = {k: {n for w in anyw[k] for n, _ in w} for k in anyw}
    dl = len(L['vcf_alt']) - len(L['vcf_ref'])
    H = {}
    for k in X['HAPS'] + ['CHM13#0']:
        ws = comp.get(k)
        if not ws:
            H[k] = dict(status='partial' if k in present else 'absent'); continue
        br = [bracket(w, P, gi, lo, hi) for w in ws]
        h = H[k] = dict(status='complete', el={e['id']: any(has(w, e) for w in ws) for e in keep},
                        joint=bool(keep) and any(all(has(w, e) for e in keep) for w in ws),
                        anyp=any(all(has(w, e) for e in ks) for ks in alts for w in ws),
                        allele=any(b is not None and b[0] == b[1] for b in br), no_bracket=all(b is None for b in br),
                        truth_seq=any(b is not None and b[0] == L['alt_hap'][ref(b[2])[1] - b0:ref(b[3])[2] + 1 - b0 + dl] for b in br))
        h['exact'] = bool(keep) and (h['truth_seq'] or (h['allele'] and g['match'] == 'with_germline'))
    C = [k for k in X['HAPS'] if H[k]['status'] == 'complete']
    cov = X['COV'][t]
    assert len(C) == int(cov['n_complete']), t
    sp = lambda k: X['META'][k.split('#')[0]]['superpopulation']
    nC = len(C)
    car = {e['id']: [k for k in C if H[k]['el'][e['id']]] for e in keep}
    jc, ac = [k for k in C if H[k]['joint']], [k for k in C if H[k]['allele']]
    pc, xc = [k for k in C if H[k]['anyp']], [k for k in C if H[k]['exact']]
    emin = min((len(c) for c in car.values()), default=None)
    fr = lambda n: round(n / nC, 4) if nC else ''
    efreq, xfreq = (fr(emin) if emin is not None else None), (fr(len(xc)) if keep and nC else None)
    support = 'none' if not keep or not emin else 'exact_allele' if xc else 'element_set_other_allele' if pc else 'recombinant_pieces'
    vaf = {so: X['VA'][(t, so)]['AF'] for so in ('d9', 'full')}
    sides = {f >= RARE_AF for f in (efreq, fr(len(pc)) if keep and nC else None, xfreq,
                                    float(vaf['full']) if vaf['full'] else None) if f is not None and f != ''}

    def indiv(cs):
        s, hs = collections.Counter(k.split('#')[0] for k in cs), set(cs)
        hom = sum(n == 2 for n in s.values())
        het = sum(1 for x, n in s.items() if n == 1 and all((f'{x}#{h}' in hs) or H[f'{x}#{h}']['status'] == 'complete' for h in '12'))
        return len(s), hom, het, len(s) - hom - het
    ch = H['CHM13#0']
    out.update(hprc_n_complete=nC, hprc_n_partial=sum(H[k]['status'] == 'partial' for k in X['HAPS']),
               hprc_n_absent=sum(H[k]['status'] == 'absent' for k in X['HAPS']), hprc_n_multi_complete=cov['n_multi_complete'],
               hprc_support=support, elem_min_carriers=emin if emin is not None else '', elem_freq_min=efreq if efreq is not None else '',
               elem_freq_class=fcls(efreq), elem_carriers_each=';'.join(str(len(c)) for c in car.values()),
               set_carriers=len(jc), set_freq=fr(len(jc)), set_carriers_any_path=len(pc), set_freq_any_path=fr(len(pc)),
               set_freq_any_path_all88=round(len(pc) / 88, 4),
               allele_carriers=len(ac), allele_freq=fr(len(ac)), allele_eq_truth_seq=sum(H[k]['truth_seq'] for k in C),
               exact_allele_carriers=len(xc), exact_allele_freq=xfreq if xfreq is not None else '', exact_allele_freq_all88=round(len(xc) / 88, 4),
               exact_freq_class=fcls(xfreq), freq_class_agree=len(sides) <= 1, d9_vcf_AF=vaf['d9'], full_vcf_AF=vaf['full'],
               allele_no_bracket=sum(H[k]['no_bracket'] for k in C), n_set_carriers_other_allele=sum(H[k]['joint'] and not H[k]['allele'] for k in C),
               allele_carrier_haps=','.join(ac), exact_allele_carrier_haps=','.join(xc))
    for nm, cs in (('allele', ac), ('set', jc), ('exact', xc)):
        out.update(dict(zip([f'{nm}_n_individuals', f'{nm}_n_hom', f'{nm}_n_het', f'{nm}_n_het_or_unknown'], indiv(cs))))
    if not keep:                                     # no somatic element: no ALT allele to compare (the path is GRCh38 there)
        for c in [c for c in out if c.startswith(('set_', 'allele_', 'exact_', 'n_set_carriers'))]:
            out[c] = ''
        ac = jc = pc = xc = []
    for s_ in SUPER:
        out.update({f'{s_}_complete': sum(sp(k) == s_ for k in C), f'{s_}_allele_carriers': sum(sp(k) == s_ for k in ac),
                    f'{s_}_set_carriers': sum(sp(k) == s_ for k in jc), f'{s_}_anypath_carriers': sum(sp(k) == s_ for k in pc),
                    f'{s_}_exact_carriers': sum(sp(k) == s_ for k in xc)})
    out.update(CHM13_status=ch['status'], CHM13_set=ch.get('joint', '') if keep else '', CHM13_set_any_path=ch.get('anyp', '') if keep else '',
               CHM13_allele=ch.get('allele', '') if keep else '', CHM13_exact=ch.get('exact', '') if keep else '')
    # VCF cross-check (allele-level GFA carriers vs carriers_any, on the GFA-complete haplotypes)
    for so in ('d9', 'full') if keep else ():
        va = X['VA'][(t, so)]
        Vs = set(va['carriers_any'].split(',')) - {'', 'CHM13#0'}
        Cs, Gs = set(C), set(ac)
        why = collections.Counter()
        for k in Gs - Vs:
            ai = X['VH'].get((t, so), {}).get(k)
            why['gfa_only:' + ('no_matching_record' if va['n_records_matching'] == '0' else 'vcf_gt_missing' if ai == '' else 'vcf_other_allele')] += 1
        for k in (Vs & Cs) - Gs:
            why['vcf_only:' + ('gfa_set_other_allele' if H[k]['joint'] else 'gfa_other_walk')] += 1
        out.update({f'{so}_matching_record': va['n_records_matching'] != '0', f'{so}_n_carriers_any': len(Vs),
                    f'{so}_agree': Gs == Vs & Cs, f'{so}_n_vcf_carriers_not_complete': len(Vs - Cs),
                    f'{so}_disagreements': ';'.join(f'{a}:{b}' for a, b in sorted(why.items())),
                    f'{so}_CHM13_carrier': 'CHM13#0' in va['carriers_any'].split(',')})
    # ---- 3 HG008-N (array call, narrow call kept)
    hn = {}
    for h in ('hap1', 'hap2'):
        key = (t, 'normal', f"{L['chrom']}_{h}")
        s, p = X['AS'].get(key), X['AP'].get(key)
        w, wsrc = (toks(p['local_path']), 'local_path') if p and p['status'] == 'ok' else \
            (toks(p['local_path_ext']), 'local_path_ext') if p and p['local_path_ext'] else (None, '')
        nin = sum(has(w, e) for e in keep) if w else None
        wel = 'not_covered' if w is None else 'no_elements' if not keep else 'all' if nin == len(keep) else 'some' if nin else 'none'
        cn, fn = carries(s, p, wel == 'all', L)
        c, fl = hap_call(ar, h, cn)
        x = ar[f'{h}_asm_call'] if ar[f'{h}_asm_call'] != 'NA' else ar[f'{h}_dip_call']
        lc = ar[f'{h}_asm_len_change'] if ar[f'{h}_asm_call'] != 'NA' else ar[f'{h}_dip_len_change']
        tn = len(L['vcf_alt']) - len(L['vcf_ref'])
        acls = 'ALT' if c == 'ALT' else 'unresolved' if c == 'unresolved' else allele_class(L, s, cn) if x == 'NA' else 'REF' if x == 'REF' else \
            'STR_other_length' if int(lc) not in (0, tn) else 'same_length_other_seq' if int(lc) == tn != 0 else 'substitution'
        hn[h] = dict(c=c, w=w, fl=fl)
        out.update({f'{h}_walk_elements': wel, f'{h}_walk_source': wsrc, f'{h}_mm2_status': s['status'] if s else 'none',
                    f'{h}_site_call': s['site_call'] if s else '', f'{h}_hap_len_change': s['hap_len_change'] if s else '',
                    f'{h}_hap_diffs': s['hap_diffs'] if s else '', f'{h}_spans_agree': s['spans_agree'] if s else '',
                    f'{h}_giraffe_event_allele': p['asm_event_allele'] if p else '', f'{h}_giraffe_ext': p['asm_ext'] if p else '',
                    f'{h}_narrow_carries_alt': cn, f'{h}_narrow_flag': fn, f'{h}_array_asm_call': ar[f'{h}_asm_call'],
                    f'{h}_array_dip_call': ar[f'{h}_dip_call'], f'{h}_array_len_change': lc, f'{h}_asm_eq_dip': ar[f'{h}_asm_eq_dip'],
                    f'{h}_carries_alt': c, f'{h}_flag': fl, f'{h}_allele_class': acls})
    cars = [h for h in ('hap1', 'hap2') if hn[h]['c'] == 'ALT']
    unres = [h for h in ('hap1', 'hap2') if hn[h]['c'] == 'unresolved']
    tic = ar['T_info_class']
    if cars:
        pat = ('no_event_hap' if not evh else 'both' if len(cars) == 2 else 'event_hap_only' if cars == [evh] else 'other_hap_only') + \
              ('+other_unresolved' if unres else '')
    elif tic == 'eq_other_hap':
        pat, cars = 'other_hap_only:patient_frame', ['hap2' if evh == 'hap1' else 'hap1']
    else:
        pat = 'unresolved' if unres else 'neither'
    cls = [out[f'{h}_allele_class'] for h in ('hap1', 'hap2')]
    evs = out.get(f'{evh}_site_call', '') if evh else ''
    ncls = 'carries_ALT' if cars else 'unresolved' if 'unresolved' in cls else 'both_REF' if cls == ['REF', 'REF'] else \
        'SNV_site_deleted_on_event_hap' if L['kind2'] == 'SNV' and evs == 'deleted' else \
        'germline_STR_other_length' if 'STR_other_length' in cls else 'other'
    sm, dp, pn = X['SUM'][t], X['DIP'][t], X['PON'][t]
    out.update(normal_carrying_haps='+'.join(cars), normal_evidence='/'.join(hn[h]['c'] for h in ('hap1', 'hap2')),
               normal_narrow_evidence='/'.join(out[f'{h}_narrow_carries_alt'] for h in ('hap1', 'hap2')), normal_allele_class=ncls,
               event_hap=evh, event_hap_pattern=pat, asm_event=L['asm_event'], asm_event_kind=event_kind(L['asm_event']),
               n_truths_same_asm_event=L['n_truths_same_asm_event'], T_info_class=tic, tumor_verdict=ar['tumor_verdict'],
               tumor_copies=ar['tumor_copies'], tumor_any_ALT_array=ar['tumor_any_ALT_array'],
               array_lo0=a_lo, array_hi0=a_hi, array_len=ar['arr_len'], HG008T_any_ALT_narrow=sm['tumor_any_ALT'],
               HG008T_statuses=sm['tumor_statuses'], spans_disagree=sm['spans_disagree'],
               dipcall_relation=dp['relation'], dipcall_records=dp['records'], dipcall_gt=dp['gt'], dipcall_ident_haps=dp['ident_haps'],
               dipcall_site_haps=dp['site_haps'], in_dip_bed=dp['in_dip_bed'],
               in_nogermline_bed=inbed(X['NOG'], L['chrom'], int(L['vcf_pos']) - 1, int(L['vcf_pos']) - 1 + len(L['vcf_ref'])),
               PoN_rule=pn['rule'], PoN_tagged=pn['rule'] not in ('', 'none'), gnomAD_AF=pn['gnomAD_AF'], CoLoRSdb_AF=pn['CoLoRSdb_AF'],
               dbSNP_AF=pn['dbSNP_AF'])
    # ---- 4 category
    sub_class, why, cat = '', '', ''
    if g['match'] == 'closest' and src != 'reads_truth_spelling':
        why = 'closest_no_read_element' if src == 'graph_closest_no_reads' else f'closest:{closest_kind or "read_allele_not_truth"}'
    elif not keep:
        why = 'reads_walk_only_germline_branches' if dropped else 'no_somatic_element'
    elif phase == 'off_event_hap':
        why = 'with_germline:germline_allele_off_event_hap'
    elif nC == 0:
        why = 'HPRC_no_complete_haplotype'
    elif pat.startswith('both'):
        why = 'evidence_conflict:ALT_on_both_normal_haps'
    elif pat.startswith('event_hap_only'):
        why = 'evidence_conflict:ALT_on_normal_event_hap'
    elif cars:
        if tic == 'novel':
            why = 'evidence_conflict:germline_ALT_vs_truth_INFO_event_novel'
        else:
            cat = 'HG008N_present_broad' if xfreq >= RARE_AF else 'HG008N_present_rare'; sub_class = pat
    elif pat == 'unresolved':
        why = 'HG008N_hap_unresolved' + (':conflict' if any(hn[h]['fl'].startswith('conflict') for h in unres) else ':array_NA')
    elif support == 'none':
        why = 'no_HPRC_carrier' + (':CHM13_carries' if ch.get('anyp') else '')
    else:
        cat, sub_class = 'HG008N_absent_HPRC_other', support
    if why:
        cat = 'other_ambiguous'
    out.update(category=cat, sub_class=sub_class, sub_reason=why, interpretation=interpret(cat, sub_class, why),
               tumor_shows_no_change=ar['tumor_verdict'].startswith(('no_change', 'all_eq')) if cat != 'HG008N_present_broad' and cat != 'HG008N_present_rare' else '')
    nodes = []
    for e in keep:
        cs = car[e['id']]
        nodes.append(dict(truth_id=t, element_id=e['id'], type=e['type'], subtype=e['subtype'], x=e['x'], y=e['y'], nodes=e['nodes'],
                          node_seqs=','.join(f'{n}{o}:{oseq(n, o)}' for n, o in toks(e['nodes'])), seq=e['seq'], replaced=e['replaced'],
                          grch38_start0=e['start0'], grch38_end0=e['end0'], in_event=e['in_event'], role=e['role'], alone=eclass(e),
                          element_source=src, category=cat, **{f'{p}_{a}': node_reads[e['id']].get(p, ('',) * 6)[i] for p in PLATS
                                                               for i, a in enumerate(('n_good', 'n_testable', 'n_contain', 'n_truth',
                                                                                      'n_truth_testable', 'n_truth_contain'))},
                          hprc_n_complete=nC, hprc_carriers=len(cs), hprc_freq=round(len(cs) / nC, 4) if nC else '',
                          **{f'{s_}_carriers': sum(sp(k) == s_ for k in cs) for s_ in SUPER},
                          hprc_any_row_carriers=sum(any(has(w, e) for w in anyw[k]) for k in X['HAPS']),
                          min_node_coverage=min((sum(n in anyn.get(k, ()) for k in X['HAPS'] + ['CHM13#0', 'GRCh38#0'])
                                                 for n, _ in toks(e['nodes'])), default=''),
                          CHM13=ch.get('el', {}).get(e['id'], ch['status']),
                          HG008N_hap1_walked='' if hn['hap1']['w'] is None else has(hn['hap1']['w'], e),
                          HG008N_hap2_walked='' if hn['hap2']['w'] is None else has(hn['hap2']['w'], e), carrier_haps=','.join(cs)))
    mem = []
    for k in X['HAPS'] + ['CHM13#0']:
        s_, h = k.split('#'); m = X['META'][s_]; hk = H[k]
        pop = (m['population_code'], m['superpopulation'] if s_ != 'CHM13' else 'CHM13', hk['status'])
        for e in keep:
            mem.append((t, e['id'], s_, h) + pop + (hk['el'][e['id']] if 'el' in hk else '', '', '', ''))
        mem.append((t, 'allele', s_, h) + pop + tuple(hk.get(x, '') if keep else '' for x in ('joint', 'allele', 'anyp', 'exact')))
    return out, nodes, mem


INTERP = {
    'other_hap_only': 'the ALT allele over the whole tandem array is the patient\'s germline allele on the non-event hap; the tumor '
                      'event hap now carries it (local conversion or recurrent slippage to the other hap\'s allele; a tumor-only '
                      'caller sees a germline allele): germline filtering by the graph',
    'other_hap_only:patient_frame': 'the tumor event-hap allele (truth INFO event applied to the normal event hap) equals the '
                                    'patient\'s other-hap germline allele over the array; the GRCh38 ALT differs from it only by '
                                    'germline differences shared by both haps: germline filtering by the graph',
    'no_event_hap': 'the ALT allele over the array is a germline allele of one normal hap; no event hap in the truth INFO',
    'exact_allele': 'the patient does not carry it; >= 1 complete HPRC haplotype spells the exact ALT allele over the array '
                    '(somatic slippage / substitution recreating a population allele kept by the d9 graph): pangenome-induced false negative',
    'element_set_other_allele': 'the patient does not carry it; HPRC haplotypes walk all somatic elements but spell another '
                                'allele over the array (e.g. another repeat length): the graph spells the somatic allele from '
                                'population nodes, not a population allele: pangenome-induced false negative',
    'recombinant_pieces': 'the patient does not carry it; every somatic element is walked by some HPRC haplotype but none walks a '
                          'whole ALT path: the ALT path is a recombination of population nodes: pangenome-induced false negative',
}
WHY = {
    'closest_no_read_element': 'no graph path spells the truth and no read walk to take the elements from',
    'closest:read_allele_not_truth': 'no graph path spells the truth; the reads walk a graph allele nearer GRCh38 / the patient\'s germline than the truth',
    'closest:read_allele_nearer_truth': 'no graph path spells the truth; the reads walk a graph allele nearer the truth, but not the truth',
    'reads_walk_only_germline_branches': 'every element of the read walk alone spells a germline allele',
    'no_somatic_element': 'the ALT path has no non-GRCh38 element',
    'with_germline:germline_allele_off_event_hap': 'the ALT path needs a germline allele that dipcall phases to the non-event hap',
    'HPRC_no_complete_haplotype': 'no HPRC haplotype traverses the window completely',
    'evidence_conflict:ALT_on_both_normal_haps': 'both normal haps carry the GRCh38 ALT over the array: homozygous germline, truth conflict',
    'evidence_conflict:ALT_on_normal_event_hap': 'the GRCh38 ALT is the normal event hap\'s own allele; the somatic change relative to the normal is the INFO event: truth representation conflict',
    'evidence_conflict:germline_ALT_vs_truth_INFO_event_novel': 'the GRCh38 ALT equals a normal hap\'s allele, but the truth INFO event gives a tumor allele that neither normal hap has',
    'HG008N_hap_unresolved:conflict': 'a normal hap: assembly and dipcall disagree over the array',
    'HG008N_hap_unresolved:array_NA': 'a normal hap has no array sequence and the narrow call is not a clear no',
    'no_HPRC_carrier': 'no complete HPRC haplotype walks the somatic elements',
    'no_HPRC_carrier:CHM13_carries': 'no complete HPRC haplotype walks the somatic elements; CHM13 does',
}


def interpret(cat, sub, why):
    if cat.startswith('HG008N_present'):
        return 'HG008-N germline allele already a graph path (' + ('rare' if cat.endswith('rare') else 'common') + \
            ' exact allele in HPRC): ' + INTERP.get(sub.split('+')[0], sub)
    if cat == 'HG008N_absent_HPRC_other':
        return INTERP[sub]
    return 'unresolved: ' + WHY.get(why, why)


def md(title, header, rows, note=''):
    s = f'\n### {title}\n\n' + (note + '\n\n' if note else '') + '| ' + ' | '.join(header) + ' |\n|' + '---|' * len(header) + '\n'
    return s + ''.join('| ' + ' | '.join(str(x) for x in r) + ' |\n' for r in rows)


def tables(rows, nodes, X):
    SC = [('chr1-22', lambda r: True), ('chr2-22', lambda r: r['chrom'] != 'chr1'), ('chr1 BED', lambda r: r['chrom'] == 'chr1' and r['in_bed'] == 'True'),
          ('chr1-22 nogermline BED', lambda r: r['in_nogermline_bed'] is True)]
    SETS = [(p, lambda r, p=p: r[f'{p}_perfect'] == 'True') for p in PLATS] + [('union', lambda r: True)]
    K2 = ('SNV', 'INDEL')
    out, summ = ['# s9 integration: graph-absorbed HG008T somatic truth alleles (d9 HPRC v1.1)\n',
                 f'RARE_AF = {RARE_AF} on the exact-allele HPRC frequency (complete haplotypes spelling the ALT allele over the '
                 'tandem array). HG008-N calls: whole-array alleles (s8b). Read data: ' +
                 ', '.join(f'{p} {len(X["SM"].get(p, {}))} loci' for p in PLATS) + '. Fixes: audit/decisions.md.\n'], {}
    sel = lambda f, kd=None: [r for r in rows if f(r) and (kd is None or r['kind2'] == kd)]
    key = lambda r: r['sub_class'] or r['sub_reason']
    for scn, scf in SC:
        hdr, body, tot = ['category'], {c: [c] for c in CATS}, ['total']
        for sn, sf in SETS:
            for kd in K2:
                x = [r for r in rows if scf(r) and sf(r) and r['kind2'] == kd]
                hdr.append(f'{sn} {kd} {scn}'); tot.append(len(x))
                cc = collections.Counter(r['category'] for r in x)
                for c in CATS:
                    body[c].append(cc[c])
                summ[f'{sn} {kd} {scn}'] = dict(cc, total=len(x))
        out.append(md(f'Categories per perfect set, {scn}', hdr, list(body.values()) + [tot]))
    for scn, scf in (SC[0], SC[2], SC[3]):
        keys = sorted({(r['category'], key(r)) for r in rows}, key=lambda k: (CATS.index(k[0]), k[1]))
        hdr = ['category', 'sub_class / sub_reason'] + [f'{sn} {kd} {scn}' for sn, _ in SETS for kd in K2]
        body = [[c, s] + [sum(1 for r in rows if scf(r) and sf(r) and r['kind2'] == kd and r['category'] == c and key(r) == s)
                          for _, sf in SETS for kd in K2] for c, s in keys]
        out.append(md(f'Categories with sub-class / sub-reason, {scn}', hdr, body))
    keys = sorted({(r['category'], key(r), r['interpretation']) for r in rows}, key=lambda k: (CATS.index(k[0]), k[1]))
    out.append(md('Interpretation per category and sub-class (union, chr1-22)', ['category', 'sub_class / sub_reason', 'SNV', 'INDEL', 'interpretation'],
                  [[c, s, len(sel(lambda r: (r['category'], key(r)) == (c, s), 'SNV')), len(sel(lambda r: (r['category'], key(r)) == (c, s), 'INDEL')), it]
                   for c, s, it in keys]))
    fb = [('0', lambda f: f == 0), ('(0,0.05)', lambda f: 0 < f < 0.05), ('[0.05,0.1)', lambda f: 0.05 <= f < 0.1), ('[0.1,0.2)', lambda f: 0.1 <= f < 0.2),
          ('[0.2,0.5)', lambda f: 0.2 <= f < 0.5), ('[0.5,0.8)', lambda f: 0.5 <= f < 0.8), ('[0.8,1]', lambda f: f >= 0.8)]
    for col, nm in (('exact_allele_freq', 'exact allele (spells the ALT allele over the array; RARE_AF level)'),
                    ('set_freq_any_path', 'set level, any ALT path (walks all somatic elements)'), ('elem_freq_min', 'element level (min over elements)'),
                    ('full_vcf_AF', 'full-graph VCF AF (all called haplotypes, matching record)')):
        body = []
        for c in CATS:
            for kd in K2:
                f = [float(r[col]) for r in sel(lambda r: r['category'] == c, kd) if r[col] != '']
                body.append([c, kd, len(f)] + [sum(g(x) for x in f) for _, g in fb] + [round(statistics.median(f), 3) if f else ''])
        out.append(md(f'HPRC frequency per category, {nm} (union, chr1-22)', ['category', 'kind', 'n'] + [b for b, _ in fb] + ['median'], body))
    body = []
    for c in CATS + ('all',):
        x = [r for r in rows if (c == 'all' or r['category'] == c) and r['n_somatic_elements'] and r['hprc_n_complete']]
        for nm, lab in (('exact', 'exact allele'), ('anypath', 'set, any ALT path')):
            O, E, V = collections.Counter(), collections.Counter(), collections.Counter()
            for r in x:
                n, tot = int(r['hprc_n_complete']), sum(int(r[f'{s}_{nm}_carriers']) for s in SUPER)
                for s in SUPER:
                    q = int(r[f'{s}_complete']) / n
                    O[s] += int(r[f'{s}_{nm}_carriers']); E[s] += tot * q
                    V[s] += tot * q * (1 - q) * (n - tot) / (n - 1) if n > 1 else 0
            body.append([c, lab, len(x), sum(O.values())] + [f"{O[s] / E[s]:.3f} (z {(O[s] - E[s]) / math.sqrt(V[s]):+.1f})" if E[s] and V[s] else '' for s in SUPER] +
                        [sum(r['CHM13_' + ('exact' if nm == 'exact' else 'set_any_path')] is True for r in x)])
    out.append(md('Superpopulation of HPRC carriers per category: observed / expected (union, chr1-22)',
                  ['category', 'level', 'loci', 'carrier haps'] + [f'{s} O/E (z)' for s in SUPER] + ['CHM13 carries'], body,
                  'expected per locus = carriers x the superpopulation\'s share of the complete haplotypes at that locus (hypergeometric variance '
                  'for z; ignores the pairing of haplotypes within individuals and linkage between loci, so z overstates significance). '
                  'HG008 donor: European; HPRC v1.1 has no EUR sample (AFR 23, AMR 16, EAS 4, SAS 1).'))
    for col, nm in (('hprc_support', 'HPRC support'), ('exact_freq_class', 'exact-allele frequency class'), ('elem_freq_class', 'element-level frequency class'),
                    ('freq_class_agree', 'frequency class agrees across element / set / exact / full-VCF'),
                    ('event_hap_pattern', 'Event-hap pattern'), ('normal_evidence', 'HG008-N carries_alt hap1/hap2 (array)'),
                    ('normal_narrow_evidence', 'HG008-N narrow (s6/s7 event window) call hap1/hap2'), ('normal_allele_class', 'HG008-N allele class (locus)'),
                    ('T_info_class', 'truth INFO event applied to the event hap (s8b)'), ('tumor_verdict', 'HG008-T contigs vs normal haps over the array (s8b)'),
                    ('tumor_any_ALT_array', 'HG008-T any contig == ALT over the array'), ('tumor_shows_no_change', 'HG008-T shows no change at the array'),
                    ('asm_event_kind', 'kind of the truth INFO event (normal-assembly frame)'), ('truth_read_support', 'truth-spelling reads on a perfect platform'),
                    ('dipcall_relation', 'Dipcall relation'), ('in_dip_bed', 'in_dip_bed'), ('in_nogermline_bed', 'in GIAB v0.2 nogermlineoverlap BED'),
                    ('PoN_tagged', 'PoN tagged (repo rule)'), ('PoN_rule', 'PoN rule tag'), ('element_source', 'element source'), ('closest_kind', 'closest: what the reads walk'),
                    ('germline_phase', 'with_germline: germline alleles on the event hap'), ('path_choice', 'path choice (read counts dropped)'), ('match', 's2 match')):
        val = (lambda r: re.sub(r':\d+/\d+$', '', r[col])) if col == 'path_choice' else (lambda r: r[col])
        vals = sorted({val(r) for r in rows}, key=str)
        out.append(md(f'{nm} x category (union, chr1-22)', [col] + [f'{c} {kd}' for c in CATS for kd in K2],
                      [[vv] + [sum(1 for r in rows if val(r) == vv and r['category'] == c and r['kind2'] == kd) for c in CATS for kd in K2] for vv in vals]))
    ids = {r['truth_id'] for r in rows}
    paf = lambda r: max([float(x) for x in (r['gnomAD_AF'], r['CoLoRSdb_AF']) if x] or [0.0])
    groups = [(f'{c}', lambda r, c=c: r['category'] == c) for c in CATS] + [('all absorbed', lambda r: True)]
    body = []
    for kd in K2:
        for nm, f in groups:
            x = sel(f, kd)
            body.append([nm, kd, len(x)] + [f'{100 * sum(g(r) for r in x) / len(x):.1f}' if x else '' for g in
                                            (lambda r: r['PoN_rule'] not in ('', 'none'), lambda r: paf(r) >= 1e-3, lambda r: paf(r) >= 0.01, lambda r: paf(r) >= 0.05)])
        y = [r for r in X['PONALL'] if r['truth_id'] not in ids and r['kind'] == kd and r['chrom'] in {f'chr{i}' for i in range(1, 23)}]   # chrX/Y rows (no truth_id) are not chr1-22 controls
        body.append(['not absorbed (other HG008T truths)', kd, len(y)] + [f'{100 * sum(g(r) for r in y) / len(y):.1f}' for g in
                                                                           (lambda r: r['rule'] not in ('', 'none'), lambda r: paf(r) >= 1e-3, lambda r: paf(r) >= 0.01, lambda r: paf(r) >= 0.05)])
    out.append(md('PoN: % tagged per category vs the non-absorbed truths (chr1-22)', ['group', 'kind', 'n', 'repo rule %', 'pop AF >= 1e-3 %',
                                                                                   'pop AF >= 0.01 %', 'pop AF >= 0.05 %'], body,
                  'repo rule = allele match, gnomAD / CoLoRSdb AF >= 1e-4, dbSNP non-somatic, 1000G (SNV; the repo INDEL path has no PoN). '
                  'pop AF = max(gnomAD, CoLoRSdb). The graph-specific loss = absorbed truths a linear tumor-only caller with that filter would still report.'))
    body = []
    for c in CATS:
        for kd in K2:
            x = sel(lambda r: r['category'] == c, kd)
            gn = [float(r['gnomAD_AF']) for r in x if r['gnomAD_AF']]; co = [float(r['CoLoRSdb_AF']) for r in x if r['CoLoRSdb_AF']]
            body.append([c, kd, len(x), sum(r['PoN_tagged'] for r in x), len(gn), round(statistics.median(gn), 4) if gn else '', sum(a >= 1e-4 for a in gn),
                         len(co), round(statistics.median(co), 4) if co else '', sum(a >= 1e-4 for a in co)])
    out.append(md('PoN rule and population AF per category (union, chr1-22)', ['category', 'kind', 'n', 'PoN tagged', 'gnomAD AF n', 'gnomAD median', 'gnomAD >=1e-4',
                                                                               'CoLoRSdb AF n', 'CoLoRSdb median', 'CoLoRSdb >=1e-4'], body))
    body = []
    for p in PLATS:
        for kd in K2:
            x = [r for r in rows if r[f'{p}_perfect'] == 'True' and r['kind2'] == kd]
            y = [r for r in x if r.get(f'{p}_read_data') is True]
            z = [r for r in y if r.get(f'{p}_n_good')]
            zt = [r for r in z if r.get(f'{p}_n_truth_reads')]
            fr = [float(r[f'{p}_truth_frac_contain_set']) for r in zt if r[f'{p}_truth_frac_contain_set'] != '']
            mv = collections.Counter(r[f'{p}_majority_vs_graph'] for r in zt)
            mc = collections.Counter(r[f'{p}_all_reads_majority_class'] for r in z)
            body.append([p, kd, len(x), len(y), len(z), len(zt), mc['germline_only'], sum(not r[f'{p}_truth_n_testable'] for r in zt), mv['same'], mv['different'],
                         round(statistics.median(fr), 3) if fr else '', sum(f >= 0.9 for f in fr), sum(f < 0.5 for f in fr)])
    out.append(md('Read-vs-graph agreement per platform (perfect loci)', ['platform', 'kind', 'perfect loci', 'with read data', 'with no-edit spanning ALT reads',
                                                                          'with truth-spelling reads', 'majority of all reads spells germline only',
                                                                          'no truth read can test an element', 'truth majority walk = graph elements', 'different',
                                                                          'median frac truth reads containing the set', 'frac >= 0.9', 'frac < 0.5'], body))
    body = []
    for so in ('d9', 'full'):
        for kd in K2:
            x = sel(lambda r: True, kd); why = collections.Counter()
            for r in x:
                for d in filter(None, r.get(f'{so}_disagreements', '').split(';')):
                    a, n = d.rsplit(':', 1); why[a] += int(n)
            x = [r for r in x if r['n_somatic_elements']]
            body.append([so, kd, len(x), sum(r[f'{so}_matching_record'] is True for r in x), sum(r[f'{so}_agree'] is True for r in x),
                         sum(r[f'{so}_agree'] is True for r in x if r[f'{so}_matching_record'] is True), '; '.join(f'{a} {n}' for a, n in why.most_common())])
    out.append(md('HPRC carriers: GFA allele level vs VCF carriers_any (on GFA-complete haplotypes; loci with somatic elements)', ['VCF', 'kind', 'loci', 'with matching record', 'identical carrier sets',
                                                                                                         'identical, with record', 'hap-level disagreements'], body))
    cov = [int(n['min_node_coverage']) for n in nodes if n['min_node_coverage'] != '']
    ce = [int(n['hprc_carriers']) for n in nodes]; ae = [int(n['hprc_any_row_carriers']) for n in nodes]
    full = {r['truth_id'] for r in rows if r['hprc_n_partial'] == 0 and r['hprc_n_absent'] == 0}
    fn = [n for n in nodes if n['truth_id'] in full]
    fcov = [int(n['min_node_coverage']) for n in fn if n['min_node_coverage'] != '']
    floor = dict(full_windows=len(full), full_window_elements=len(fn), full_window_branch_node_coverage_min=min(fcov) if fcov else None,
                 full_window_node_coverage_hist=dict(sorted(collections.Counter(min(c, 15) for c in fcov).items())),
                 full_window_hprc_carriers_min=min((int(n['hprc_carriers']) for n in fn), default=None),
                 full_window_hprc_carriers_hist=dict(sorted(collections.Counter(min(int(n['hprc_carriers']), 15) for n in fn).items())),
                 full_window_not_on_CHM13_hprc_carriers_min=min((int(n['hprc_carriers']) for n in fn if n['CHM13'] is not True), default=None),
                 full_window_on_CHM13_below9=sum(int(n['hprc_carriers']) < 9 for n in fn if n['CHM13'] is True),
                 all_not_on_CHM13_node_coverage_lt9=sum(n['min_node_coverage'] != '' and int(n['min_node_coverage']) < 9 for n in nodes if n['CHM13'] is not True),
                 n_elements=len(nodes), branch_elements_node_coverage_min=min(cov) if cov else None,
                 node_coverage_lt9=sum(c < 9 for c in cov), node_coverage_hist=dict(sorted(collections.Counter(min(c, 15) for c in cov).items())),
                 complete_carriers_min=min(ce), complete_carriers_lt9=sum(c < 9 for c in ce), any_row_carriers_min=min(ae), any_row_carriers_lt9=sum(c < 9 for c in ae))
    out.append(md('d9 frequency-filter floor (somatic elements)', ['measure', 'value'], [[k, v] for k, v in floor.items()],
                  'node coverage = HPRC + CHM13 + GRCh38 haplotypes whose s4 rows (complete or partial path) visit the branch node (min over the nodes of the element); a lower bound '
                  '(s4 partial paths stop 500 nodes from an anchor). full_window = loci where all 88 HPRC haplotypes have a complete traversal, where the counts are exact. '
                  'A whole-GFA scan (audit hprc_membership a5) found >= 7 HPRC haplotypes for every branch node not on a reference path.'))
    body = [[p, len(X['SM'].get(p, {})), sum(r.get(f'{p}_read_data') is True and bool(r.get(f'{p}_n_good')) for r in rows),
             sum(r[f'{p}_perfect'] == 'True' for r in rows), sum(r[f'{p}_perfect'] == 'True' and r.get(f'{p}_read_data') is True for r in rows)] for p in PLATS]
    out.append(md('Read data per platform', ['platform', 'loci with read data', 'loci with no-edit spanning ALT reads', 'perfect loci', 'perfect loci with read data'], body))
    return ''.join(out), summ, floor


def work(i):
    return locus(X, X['L'][i])


def main(nproc, test=None):
    global X
    X = load()
    if test:
        X['L'] = [L for L in X['L'] if L['truth_id'] in test]
    print('loaded', len(X['L']), 'loci; read data', {p: len(v) for p, v in X['SM'].items()}, flush=True)
    rows, nodes, mem = [], [], []
    with Pool(nproc, initializer=S2.init) as pool:                 # one sqlite connection per worker (s2 init)
        for i, (r, n, m) in enumerate(pool.imap(work, range(len(X['L'])), chunksize=2)):
            rows.append(r); nodes += n; mem += m
            if i % 250 == 0:
                print(i, len(X['L']), flush=True)
    if test:
        for r in rows:
            print('\n'.join(f'  {k}: {str(v)[:300]}' for k, v in r.items() if k not in X['V'][r['truth_id']] or k == 'truth_id')); print()
        return
    cols = list(X['V'][X['L'][0]['truth_id']].keys())
    for r in rows:
        cols += [c for c in r if c not in cols]
    with open(f'{D}/per_variant.tsv', 'w') as w:
        w.write('\t'.join(cols) + '\n')
        for r in rows:
            w.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
    with open(f'{D}/per_node.tsv', 'w') as w:
        nc = list(nodes[0].keys()); w.write('\t'.join(nc) + '\n')
        for n in nodes:
            w.write('\t'.join(str(n[c]) for c in nc) + '\n')
    with gzip.open(f'{D}/hprc_membership.tsv.gz', 'wt') as w:
        w.write('truth_id\telement_id\tsample\thap\tpopulation_code\tsuperpopulation\ttraversal_status\tcarries_element\tcarries_alt_allele\t'
                'carries_set_any_path\tcarries_exact_allele\n')
        for m in mem:
            w.write('\t'.join(map(str, m)) + '\n')
    txt, summ, floor = tables(rows, nodes, X)
    open(f'{D}/tables.md', 'w').write(txt)
    perfect = {p: {kd: sum(r[f'{p}_perfect'] == 'True' and r['kind2'] == kd for r in rows) for kd in ('SNV', 'INDEL')} for p in PLATS}
    for p in PLATS:
        for kd in ('SNV', 'INDEL'):
            assert summ[f'{p} {kd} chr1-22']['total'] == perfect[p][kd] == sum(summ[f'{p} {kd} chr1-22'].get(c, 0) for c in CATS)
    cnt = lambda col, f=lambda r: True: {f'{a}|{b}|{k}': n for (a, b, k), n in sorted(collections.Counter((r['category'], r[col], r['kind2']) for r in rows if f(r)).items())}
    S = dict(n_loci=len(rows), RARE_AF=RARE_AF, perfect=perfect, read_data={p: len(X['SM'].get(p, {})) for p in PLATS},
             categories=summ, floor=floor,
             vcf_agree={so: {kd: sum(r.get(f'{so}_agree') is True for r in rows if r['kind2'] == kd) for kd in ('SNV', 'INDEL')} for so in ('d9', 'full')},
             element_source=dict(collections.Counter(r['element_source'] for r in rows)), path_choice=dict(collections.Counter(r['path_choice'].split(':')[0] for r in rows)),
             sub_class_or_reason={f'{a}|{b}|{k}': n for (a, b, k), n in sorted(collections.Counter((r['category'], r['sub_class'] or r['sub_reason'], r['kind2'])
                                                                                                 for r in rows).items())},
             event_hap_pattern=cnt('event_hap_pattern', lambda r: r['category'].startswith('HG008N_present')),
             tumor_verdict=cnt('tumor_verdict'), T_info_class=cnt('T_info_class'), truth_read_support=cnt('truth_read_support'),
             asm_event_kind=cnt('asm_event_kind'), normal_allele_class=cnt('normal_allele_class'),
             narrow_vs_array={f'{a}->{b}': n for (a, b), n in sorted(collections.Counter((r[f'{h}_narrow_carries_alt'], r[f'{h}_carries_alt'])
                                                                                        for r in rows for h in ('hap1', 'hap2')).items())},
             freq_class_agree_present=sum(r['freq_class_agree'] is True for r in rows if r['category'].startswith('HG008N_present')))
    json.dump(S, open(f'{D}/summary.json', 'w'), indent=1, default=str)
    rng, lines = random.Random(9), []
    for c in CATS:
        x = [r for r in rows if r['category'] == c]
        pick = (rng.sample([r for r in x if r['kind2'] == 'SNV'], min(1, sum(r['kind2'] == 'SNV' for r in x))) +
                rng.sample([r for r in x if r['kind2'] == 'INDEL'], min(3, sum(r['kind2'] == 'INDEL' for r in x))))[:3]
        for r in pick:
            lines.append(f"== {c} truth {r['truth_id']} {r['chrom']}:{r['vcf_pos']} {r['vcf_ref']}>{r['vcf_alt']} {r['kind2']} perfect P/O/I "
                         f"{r['PacBio_perfect']}/{r['ONT_perfect']}/{r['Illumina_perfect']}")
            for k in ('match', 'element_source', 'path_choice', 'closest_kind', 'somatic_elements', 'element_alone', 'dropped_germline_elements', 'alt_path',
                      'Illumina_n_good', 'Illumina_n_truth_reads', 'Illumina_n_germline_only_reads', 'Illumina_truth_n_contain_set', 'Illumina_majority_class',
                      'truth_read_support', 'hprc_n_complete', 'elem_carriers_each', 'elem_freq_min', 'set_carriers_any_path', 'exact_allele_carriers',
                      'exact_allele_freq', 'hprc_support', 'exact_allele_carrier_haps', 'd9_n_carriers_any', 'd9_agree', 'CHM13_exact',
                      'array_lo0', 'array_hi0', 'hap1_array_asm_call', 'hap1_array_dip_call', 'hap1_array_len_change', 'hap1_narrow_carries_alt', 'hap1_carries_alt',
                      'hap2_array_asm_call', 'hap2_array_dip_call', 'hap2_array_len_change', 'hap2_narrow_carries_alt', 'hap2_carries_alt',
                      'event_hap', 'event_hap_pattern', 'T_info_class', 'tumor_verdict', 'tumor_copies', 'normal_allele_class', 'asm_event', 'asm_event_kind',
                      'dipcall_relation', 'dipcall_records', 'dipcall_gt', 'PoN_rule', 'gnomAD_AF', 'CoLoRSdb_AF', 'sub_class', 'sub_reason'):
                lines.append(f'  {k}: {str(r.get(k, ""))[:300]}')
    open(f'{D}/spot_checks.txt', 'w').write('\n'.join(lines) + '\n')
    print(json.dumps({k: S[k] for k in ('perfect', 'read_data', 'element_source', 'path_choice', 'vcf_agree')}, default=str))
    for k in ('PacBio SNV chr1-22', 'PacBio INDEL chr1-22', 'ONT SNV chr1-22', 'ONT INDEL chr1-22', 'Illumina SNV chr1-22', 'Illumina INDEL chr1-22',
              'union SNV chr1-22', 'union INDEL chr1-22', 'union SNV chr1 BED', 'union INDEL chr1 BED'):
        print(k, summ[k])


if __name__ == '__main__':
    main(int(sys.argv[1]) if sys.argv[1:] else int(os.environ.get('SLURM_CPUS_PER_TASK', 4)), set(sys.argv[2].split(',')) if sys.argv[2:] else None)
