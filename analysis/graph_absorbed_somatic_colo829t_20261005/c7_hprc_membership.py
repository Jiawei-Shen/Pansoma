"""COLO829T step 7: which HPRC v1.1 haplotypes walk the somatic graph elements and carry the ALT allele of each absorbed
truth, and the cross-evaluation with the patient's normal (COLO829BL): germline filtering vs pangenome-induced false
negative vs ambiguous.

Port of the HPRC part of analysis/graph_absorbed_somatic_20261001/s9_integrate.py (HG008T: method steps 1, 2 and 4, same
rules, same RARE_AF). The HG008-N part (s6-s8b assemblies, giraffe walks of the normal, truth INFO event, tumor contigs)
is replaced by c3's COLO829BL whole-array alleles. There is no read part: the COLO829T reads were not re-decoded.
Inputs ($D): variant_set / loci (c0, c1b); graph_paths (c2); normal_frame.tsv (c3: arr_lo0 / arr_hi0, hapX / hapY array
allele, event hap, pf_label, alt_len_haps); hprc_local_paths.tsv.gz / hprc_haplotypes / hprc_window_coverage (c5);
hprc_vcf_alleles / hprc_vcf_haplotypes.tsv.gz / hprc_sample_metadata (c6); the COLO829BL dipcall VCF of c2 (PASS alleles
inside each window, the 6 nearest to the truth = the c2 rule; GT hap1|hap2 = hapY|hapX, checked in c3). Graph helpers:
c2_graph_paths.py imported (classify.py ref / oseq; c2 walk(), elements(), apply(), lev()).
Method
 1 Somatic elements (c2 element rules on the whole path: 'B:x:y:>n..' = run of non-GRCh38 nodes between GRCh38 nodes x, y;
   'S:x:y' = consecutive GRCh38 nodes not adjacent in GRCh38):
     exact          every element of the chosen ALT path. Several paths (c2 walk() re-run with the same caps; count and
                    primary checked against graph_paths): no read data decides, so the c2 primary (path_choice
                    no_read_data:primary; n_candidate_paths recorded; the set level below uses every path)
     with_germline  the same, phase-consistent paths first (every germline allele of the path on the c3 event hap by the
                    dipcall GT; none -> other_ambiguous), minus the elements whose sequence alone (the element applied
                    alone to GRCh38 over the window) = GRCh38 + a non-empty subset of the path's germline alleles
     closest        the c2 closest primary minus germline-alone elements (of the 6 nearest alleles); reported, but the
                    locus is other_ambiguous (no d9 path spells the truth); closest_kind = graph_allele_nearer_truth /
                    graph_allele_not_truth by Levenshtein of the kept set applied to GRCh38 vs truth(+germline) and
                    GRCh38(+germline) haplotypes
   alone (per element): truth / germline / truth+germline / partial (one piece of a multi-element path).
 2 HPRC (88 haplotypes = 44 x 2; CHM13 separately; GRCh38 never): c5 complete traversals only (a haplotype with several
   complete traversals counts if any of them does). Element level: the traversal contains the element (x, the node run,
   y consecutive; locus: min over the elements, elem_freq_min); set level: all elements of the chosen path
   (set_carriers) / of any candidate ALT path (set_carriers_any_path); allele level: bracket = the last GRCh38 node of the
   ALT path starting before lo and the first after it ending at / after hi that the haplotype walk also visits (lo, hi =
   the c3 tandem array united with the event window and the elements' spans, clipped to the graph window); carrier = same
   spelled sequence from L to R as the ALT path (allele) or as GRCh38 + truth (truth_seq); COLO829T addition (a strict
   superset of the HG008 rule): a walk spelling the whole ALT path / GRCh38 + truth over the whole window also counts
   (exact_window_only = haplotypes found only this way: compensating graph edits outside the bracket, e.g. 29562); and
   (core rule, added after the HPRC audit) a walk whose core() over [lo, hi) equals that of the ALT path / of GRCh38 +
   truth: core() maps the walk onto GRCh38 (forward GRCh38 nodes clipped to [lo, hi), plus each non-GRCh38 run whose
   replaced span overlaps it), so a germline SNV / indel just outside [lo, hi) that moves the node bracket outward no
   longer hides a carrier (exact_core_only = haplotypes found only this way). Exact allele = truth_seq, or the path
   allele for with_germline (truth + the patient's germline alleles), per traversal; a haplotype counts if any of its
   complete traversals does. PSV flag: exact_carriers_alt_ref_copy = carriers with another complete traversal that
   spells GRCh38 over [lo, hi) (an ALT copy and a REF copy); multicopy_psv_like = every exact carrier is such a
   haplotype (a paralogous sequence variant rather than an allele of the locus). Frequencies over the complete haplotypes
   (conditional, upper-biased) and over all 88 (lower bound). hprc_support, first that applies: exact_allele (>= 1
   exact carrier) / element_set_other_allele (>= 1 haplotype walks a full element set of any candidate ALT path but
   spells another allele over the array) / recombinant_pieces (every element of the chosen path has carriers, no
   haplotype has a full set) / none. RARE_AF = 0.20 on exact_allele_freq; freq_class_agree =
   the element / set-any-path / exact / full-VCF frequencies fall on the same side of 0.20.
   VCF cross-check: carriers_any of the d9 and full VCFs (c6) vs allele-level carriers, on the haplotypes complete in the
   GFA. d9 floor: haplotypes (HPRC + CHM13 + GRCh38) whose c5 rows (complete or partial paths: a lower bound) visit each
   branch node of a somatic element.
   Window-level check (per traversal at loci with somatic elements; violations counted and printed, 0 expected): a
   complete walk spelling ref_hap over the whole window contains no element and is no allele / truth_seq carrier
   (window_alt_outside_bracket = walks spelling alt_hap whose bracket does not: the cases the addition above rescues).
 3 COLO829BL (c3; no node-level walk of the normal exists): carrying haps = array allele ALT (GRCh38 + truth over the
   whole tandem array, germline differences outside it reverted); unresolved = allele NA. Pattern hapX_only / hapY_only /
   both (+other_unresolved) / unresolved / neither. with_germline loci (added after the audit, the analogue of HG008T's
   other_hap_only:patient_frame): a hap whose array sequence (c3 hap*_seq) equals the ALT path's allele over the array
   (c3 R_seq + truth + the path's germline alleles) carries it too (normal_carries_via = path_allele, sub_class
   <hap>_only:path_allele). normal_allele_class = carries_ALT / germline_other_length /
   germline_same_length_other_seq / both_REF / one_hap_NA / both_NA (c4 germline status, same order).
   germline_like_alt_len = alt_len_haps (a hap has exactly the ALT length, other bases) and no hap carries the ALT.
 4 Category, first rule that applies -> other_ambiguous with sub_reason:
     closest:<closest_kind> (no d9 path spells the truth); no kept element (only_germline_elements / no_somatic_element);
     with_germline:germline_allele_off_event_hap; HPRC_no_complete_haplotype; evidence_conflict:ALT_on_both_normal_haps;
   normal_present_broad / _rare   a COLO829BL hap carries the ALT ('normal carries ALT'), exact_allele_freq >= / < RARE_AF
   then normal_hap_unresolved:<one_hap_NA / both_NA>; no_HPRC_carrier[:CHM13_carries] (hprc_support none);
   normal_absent_HPRC_other       sub_class = hprc_support (flag germline_like_alt_len)
Outputs ($D): hprc_per_variant.tsv, hprc_per_node.tsv, hprc_membership.tsv.gz, hprc_tables.md, hprc_summary.json,
hprc_spot_checks.txt; check mode: hprc_check_kmer.tsv (kcheck(): exact-allele carriers re-derived by sequence between unique
16-mer anchors, independent of the node bracket, compared per locus with hprc_per_variant).
Run: python c7_hprc_membership.py [n_workers] [truth_id,...] (default 4; with ids: printout only, no files; login node);
python c7_hprc_membership.py check [n_workers] (after the main run).
O/E tables: + 95 % interval of O/E from 2,000 locus bootstrap resamples (oe_ci, seed 20261005).
Complete traversals include c5's bypass rows (a germline variant over the window edge skips an anchor node).
Assumptions: HPRC partial traversals are outside the denominators; the event hap is c3's derived one (the hap closest to
GRCh38 + truth), not a truth INFO event, so 'ALT on the event hap' is no conflict here (the carrying hap is the closest
by construction); dipcall GT order = hapY|hapX; multi-path loci take the c2 primary (no reads to choose).
"""
import collections, csv, gzip, importlib.util, itertools, json, math, os, random, re, statistics, sys
from bisect import bisect_left
from multiprocessing import Pool
import pysam
sys.dont_write_bytecode = True
D = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005'
A = '/scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_colo829t_20261005'
PLATS, RARE_AF, SUPER, NH = ('fiberseq', 'ONT', 'Illumina'), 0.20, ('AFR', 'AMR', 'EAS', 'SAS'), ('hapX', 'hapY')
CATS = ('normal_present_broad', 'normal_present_rare', 'normal_absent_HPRC_other', 'other_ambiguous')
GTI = {'hapY': 0, 'hapX': 1}                                   # dipcall GT = hapY|hapX
csv.field_size_limit(sys.maxsize)
_s = importlib.util.spec_from_file_location('c2', f'{A}/c2_graph_paths.py'); S2 = importlib.util.module_from_spec(_s); _s.loader.exec_module(S2)
ref, oseq, tag = S2.ref, S2.oseq, S2.tag
rd = lambda f: list(csv.DictReader((gzip.open(f, 'rt') if f.endswith('.gz') else open(f)), delimiter='\t'))
toks = lambda w: tuple((int(n), '+' if o == '>' else '-') for o, n in re.findall(r'([<>])(\d+)', w or ''))
fcls = lambda f: 'NA' if f is None or f == '' else '<0.2' if f < RARE_AF else '0.2-0.5' if f < 0.5 else '>=0.5'
_sp = {}


def spell(w):
    if w not in _sp:
        _sp[w] = ''.join(oseq(n, o) for n, o in w).upper()
    return _sp[w]


def prep(e):
    e = dict(e, x=int(e['x']), y=int(e['y']))
    e['pat'] = ((e['x'], '+'),) + toks(e['nodes']) + ((e['y'], '+'),)
    return e


def has(w, e):
    p = e['pat']; k = len(p)
    return any(w[i:i + k] == p for i in range(len(w) - k + 1) if w[i] == p[0])


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


def germline(loci):
    """PASS COLO829BL dipcall alleles inside each window, the 6 nearest to the truth (c2 rule, same order), with
    germ_gt[allele] = (on hapY, on hapX) from the phased GT (hap1|hap2 = hapY|hapX)."""
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
        L['germ'] = sorted(L['germ'], key=lambda v: abs(v[0] - p0))[:6]


def phase_ok(L, used, evh):
    """True / False: every germline allele of the path lies on the event hap (dipcall GT); None: no event hap / no GT."""
    if not used or evh not in GTI:
        return None
    g = [L['germ_gt'].get(tuple(al)) for al in used]
    if any(x is None or len(x) < 2 for x in g):
        return None
    return all(x[GTI[evh]] for x in g)


def candidates(L, g):
    """[(path tokens, germline alleles of the path)] for exact / with_germline, the c2 primary first."""
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
             G={r['truth_id']: r for r in rd(f'{D}/graph_paths.tsv')}, NF={r['truth_id']: r for r in rd(f'{D}/normal_frame.tsv')},
             COV={r['truth_id']: r for r in rd(f'{D}/hprc_window_coverage.tsv')},
             VA={(r['truth_id'], r['source']): r for r in rd(f'{D}/hprc_vcf_alleles.tsv')}, VH=collections.defaultdict(dict),
             META={r['sample']: r for r in rd(f'{D}/hprc_sample_metadata.tsv')}, HP=collections.defaultdict(list))
    assert [L['truth_id'] for L in X['L']] == list(X['COV']) and set(X['V']) == set(X['G']) == {L['truth_id'] for L in X['L']}
    for r in rd(f'{D}/hprc_vcf_haplotypes.tsv.gz'):
        X['VH'][(r['truth_id'], r['source'])][f"{r['sample']}#{r['hap']}"] = r['allele_index']
    for r in rd(f'{D}/hprc_local_paths.tsv.gz'):
        X['HP'][r['truth_id']].append((f"{r['sample']}#{r['hap']}", r['status'], r['local_path'], r['partial_path']))
    X['HAPS'] = [f"{r['sample']}#{r['hap']}" for r in rd(f'{D}/hprc_haplotypes.tsv') if r['reference'] == 'False']
    assert len(X['HAPS']) == 88
    germline(X['L'])
    return X


def bracket(w, P, gi, lo, hi):
    """(walk seq, ALT path seq, L, R), L..R inclusive: L = the last GRCh38 node of P starting before lo, R = the first after
    it ending at / after hi, both visited by the walk; None if there is no such pair."""
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


def core(w, chrom, lo, hi):
    """the walk's sequence over GRCh38 [lo, hi): its forward GRCh38 nodes clipped to [lo, hi), plus each run of other nodes
    whose replaced GRCh38 span (end of the GRCh38 node before it .. start of the one after it) overlaps [lo, hi), or for an
    insertion lies in [lo, hi] (a run crossing lo / hi is taken whole); None unless the walk has a GRCh38 node starting
    before lo and one ending after hi. Germline variants outside [lo, hi) do not change it."""
    out, run, prev, left, right = [], [], None, False, False
    for n, o in w:
        r = ref(n) if o == '+' else None
        if not r or r[0] != chrom:
            run.append(oseq(n, o)); continue
        st, en = r[1], r[2] + 1
        if prev is not None:
            a, b = min(prev, st), max(prev, st)
            if (a < hi and b > lo) or (a == b and lo <= a <= hi):
                out.append(''.join(run))
        run, prev = [], en                                  # a run before the first GRCh38 node lies before lo (left)
        left |= st < lo; right |= en > hi
        if en > lo and st < hi:
            out.append(oseq(n, o)[max(st, lo) - st:min(en, hi) - st])
    return ''.join(out).upper() if left and right else None


def normal_class(a):
    return ('carries_ALT' if 'ALT' in a else 'germline_other_length' if any(x.startswith('germline_len') for x in a)
            else 'germline_same_length_other_seq' if 'germline_seq' in a else 'both_REF' if a == ['REF', 'REF']
            else 'both_NA' if a == ['NA', 'NA'] else 'one_hap_NA')


def locus(X, L):
    t, ev_lo, ev_hi = L['truth_id'], int(L['ev_lo0']), int(L['ev_hi0'])
    g, v, nf = X['G'][t], X['V'][t], X['NF'][t]
    b0 = int(L['win_start0'])
    out = dict(v, match=g['match'], n_alt_paths=g['n_alt_paths'], germline_used=g['germline_used'], germline_gt=g['germline_gt'],
               event_net=g['event_net'], truth_net=g['truth_net'])
    gs6 = gseqs(L, L['germ'])
    evh = nf['event_hap']
    # ---- 1 somatic elements
    closest_kind, phase = '', ''
    if g['match'] in ('exact', 'with_germline'):
        sc = []
        for i, (P, used) in enumerate(candidates(L, g)):
            els = [prep(e) for e in S2.elements(P, ev_lo, ev_hi)]
            gu = gseqs(L, used)[0] if used else set()
            ph = phase_ok(L, used, evh)
            sc.append((ph is not False, -i, P, [e for e in els if alone(L, e) not in gu], els, used, ph))
        ph_ok, mi, P, keep, els, used, ph = max(sc, key=lambda z: z[:2])
        choice = 'single_path' if len(sc) == 1 else 'no_read_data:' + ('primary' if mi == 0 else 'alternative:phase')
        phase = '' if g['match'] == 'exact' else 'no_event_hap_or_GT' if ph is None else 'consistent' if ph else 'off_event_hap'
        src, alts = 'graph', [z[3] for z in sc if z[3]]
        out.update(n_candidate_paths=len(sc), n_phase_consistent_paths=sum(z[0] for z in sc) if g['match'] == 'with_germline' else '')
    else:
        P, src, choice = toks(g['primary_path']), 'graph_closest_primary', 'closest_primary'
        els = [prep(e) for e in S2.elements(P, ev_lo, ev_hi)]
        keep = [e for e in els if alone(L, e) not in gs6[0]]
        alts = [keep] if keep else []
        sp = S2.apply(L['ref_hap'], b0, [(e['start0'], L['ref_hap'][e['start0'] - b0:e['end0'] - b0], e['seq']) for e in keep])
        if sp is not None:
            dt = min(S2.lev(sp, x) for x in {L['alt_hap']} | gs6[1]); dn = min(S2.lev(sp, x) for x in {L['ref_hap']} | gs6[0])
            closest_kind = 'graph_allele_nearer_truth' if dt < dn else 'graph_allele_not_truth'
            out.update(closest_lev_truth=dt, closest_lev_nontruth=dn)
        out.update(n_candidate_paths=g['n_alt_paths'], n_phase_consistent_paths='')
    dropped = [e for e in els if e not in keep]
    eclass = lambda e: (lambda s: 'truth' if s == L['alt_hap'] else 'germline' if s in gs6[0] else 'truth+germline' if s in gs6[1] else 'partial')(alone(L, e))
    out.update(element_source=src, path_choice=choice, closest_kind=closest_kind, germline_phase=phase, alt_path=tag(P),
               somatic_elements=';'.join(e['id'] for e in keep), n_somatic_elements=len(keep),
               element_subtypes=';'.join(e['subtype'] for e in keep), element_alone=';'.join(eclass(e) for e in keep),
               element_in_event=';'.join(str(e['in_event']) for e in keep), dropped_germline_elements=';'.join(e['id'] for e in dropped))
    # ---- 2 HPRC membership
    sp_ = [z for e in keep for z in (e['start0'], e['end0'])]
    w_lo, w_hi = b0, int(L['win_end0'])
    a_lo, a_hi = int(nf['arr_lo0']), int(nf['arr_hi0'])
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
    H, chk = {}, collections.Counter()
    tseq = lambda b: b is not None and b[0] == L['alt_hap'][ref(b[2])[1] - b0:ref(b[3])[2] + 1 - b0 + dl]
    chrom, wg = L['chrom'], g['match'] == 'with_germline'
    rcore = L['ref_hap'][lo - b0:hi - b0]
    tcore, pcore = S2.apply(rcore, lo, [som(L)]), core(P, chrom, lo, hi)
    for k in X['HAPS'] + ['CHM13#0']:
        ws = comp.get(k)
        if not ws:
            H[k] = dict(status='partial' if k in present else 'absent'); continue
        br, wa, co = [bracket(w, P, gi, lo, hi) for w in ws], [spell(w) for w in ws], [core(w, chrom, lo, hi) for w in ws]
        a0 = [(b is not None and b[0] == b[1]) or s == spell(P) for b, s in zip(br, wa)]          # HG008 rule + whole window
        t0 = [tseq(b) or s == L['alt_hap'] for b, s in zip(br, wa)]
        a1 = [a or (c is not None and c == pcore) for a, c in zip(a0, co)]                       # + the core rule
        t1 = [t or (c is not None and c == tcore) for t, c in zip(t0, co)]
        ex = [t or (a and wg) for t, a in zip(t1, a1)]                                            # per traversal
        h = H[k] = dict(status='complete', el={e['id']: any(has(w, e) for w in ws) for e in keep},
                        joint=bool(keep) and any(all(has(w, e) for e in keep) for w in ws),
                        anyp=any(all(has(w, e) for e in ks) for ks in alts for w in ws),
                        allele=any(a1), no_bracket=all(b is None for b in br), truth_seq=any(t1),
                        window_only=not any(tseq(b) for b in br) and L['alt_hap'] in wa,
                        core_only=any(ex) and not any(t or (a and wg) for t, a in zip(t0, a0)), n_trav=len(ws),
                        ref_copy=any(not e and (c == rcore or s == L['ref_hap']) for e, c, s in zip(ex, co, wa)))
        h['exact'] = bool(keep) and any(ex)
        for w, b, s, c in zip(ws, br, wa, co) if keep else ():   # window-level checks (per traversal; violations counted, 0 expected)
            if s == L['alt_hap']:
                chk['window_alt'] += 1; chk['window_alt_outside_bracket'] += not tseq(b)
            if s == L['ref_hap']:
                chk['window_ref'] += 1
                chk['VIOLATION_window_ref_carries'] += (any(has(w, e) for e in keep) or (b is not None and b[0] == b[1]) or tseq(b)
                                                        or (c is not None and c in (tcore, pcore)))
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
    support = ('none' if not keep else 'exact_allele' if xc else 'element_set_other_allele' if pc else 'recombinant_pieces' if emin
               else 'none')
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
               exact_window_only=sum(H[k]['window_only'] for k in C), exact_core_only=sum(H[k]['core_only'] for k in C if H[k]['exact']),
               exact_carriers_multicopy=sum(H[k]['n_trav'] > 1 for k in xc), exact_carriers_alt_ref_copy=sum(H[k]['ref_copy'] for k in xc),
               multicopy_psv_like=bool(xc) and all(H[k]['ref_copy'] for k in xc),
               allele_carrier_haps=','.join(ac), exact_allele_carrier_haps=','.join(xc),
               window_checks=';'.join(f'{a}:{b}' for a, b in sorted(chk.items()) if b))
    for nm, cs in (('allele', ac), ('set', jc), ('exact', xc)):
        out.update(dict(zip([f'{nm}_n_individuals', f'{nm}_n_hom', f'{nm}_n_het', f'{nm}_n_het_or_unknown'], indiv(cs))))
    if not keep:                                     # no somatic element: no ALT allele to compare
        for c in [c for c in out if c.startswith(('set_', 'allele_', 'exact_', 'n_set_carriers'))]:
            out[c] = ''
        ac = jc = pc = xc = []
    for s_ in SUPER:
        out.update({f'{s_}_complete': sum(sp(k) == s_ for k in C), f'{s_}_allele_carriers': sum(sp(k) == s_ for k in ac),
                    f'{s_}_set_carriers': sum(sp(k) == s_ for k in jc), f'{s_}_anypath_carriers': sum(sp(k) == s_ for k in pc),
                    f'{s_}_exact_carriers': sum(sp(k) == s_ for k in xc)})
    out.update(CHM13_status=ch['status'], CHM13_set=ch.get('joint', '') if keep else '', CHM13_set_any_path=ch.get('anyp', '') if keep else '',
               CHM13_allele=ch.get('allele', '') if keep else '', CHM13_exact=ch.get('exact', '') if keep else '')
    for so in ('d9', 'full') if keep else ():      # VCF cross-check (allele-level GFA carriers vs carriers_any)
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
    # ---- 3 COLO829BL (c3 array alleles)
    al = [nf[f'{h}_allele'] for h in NH]
    cars = [h for h, a in zip(NH, al) if a == 'ALT']
    via = 'ALT' if cars else ''
    if not cars and g['match'] == 'with_germline' and used:  # the ALT path's own allele over the array (truth + its germline alleles)
        pa = S2.apply(nf['R_seq'], a_lo, [som(L)] + list(used))
        cars = [h for h in NH if pa is not None and nf[f'{h}_seq'] == pa]
        via = 'path_allele' if cars else ''
    unres = [h for h in NH if nf[f'{h}_allele'] == 'NA' and h not in cars]
    pat = ('both' if len(cars) == 2 else f'{cars[0]}_only') + ('+other_unresolved' if unres else '') if cars else 'unresolved' if unres else 'neither'
    out.update(hapX_allele=al[0], hapY_allele=al[1], hapX_src=nf['hapX_src'], hapY_src=nf['hapY_src'], event_hap=evh,
               pf_label=nf['pf_label'], pf_event=nf['pf_event'], alt_len_haps=nf['alt_len_haps'], normal_carrying_haps='+'.join(cars), normal_carries_via=via,
               normal_pattern=pat, normal_allele_class=normal_class(al), germline_like_alt_len=bool(nf['alt_len_haps']) and not cars,
               array_lo0=a_lo, array_hi0=a_hi, array_len=nf['arr_len'], nf_flags=nf['flags'])
    # ---- 4 category
    sub_class, why, cat = '', '', ''
    if g['match'] == 'closest':
        why = f'closest:{closest_kind or "no_element"}'
    elif not keep:
        why = 'only_germline_elements' if dropped else 'no_somatic_element'
    elif phase == 'off_event_hap':
        why = 'with_germline:germline_allele_off_event_hap'
    elif nC == 0:
        why = 'HPRC_no_complete_haplotype'
    elif len(cars) == 2:
        why = 'evidence_conflict:ALT_on_both_normal_haps'
    elif cars:
        cat = 'normal_present_broad' if xfreq >= RARE_AF else 'normal_present_rare'; sub_class = pat + (':path_allele' if via == 'path_allele' else '')
    elif unres:
        why = 'normal_hap_unresolved:' + ('both_NA' if len(unres) == 2 else 'one_hap_NA')
    elif support == 'none':
        why = 'no_HPRC_carrier' + (':CHM13_carries' if ch.get('anyp') else '')
    else:
        cat, sub_class = 'normal_absent_HPRC_other', support
    if why:
        cat = 'other_ambiguous'
    out.update(category=cat, sub_class=sub_class, sub_reason=why, interpretation=interpret(cat, sub_class, why))
    nodes = []
    for e in keep:
        cs = car[e['id']]
        nodes.append(dict(truth_id=t, element_id=e['id'], type=e['type'], subtype=e['subtype'], x=e['x'], y=e['y'], nodes=e['nodes'],
                          node_seqs=','.join(f'{n}{o}:{oseq(n, o)}' for n, o in toks(e['nodes'])), seq=e['seq'], replaced=e['replaced'],
                          grch38_start0=e['start0'], grch38_end0=e['end0'], in_event=e['in_event'], role=e['role'], alone=eclass(e),
                          element_source=src, category=cat, hprc_n_complete=nC, hprc_carriers=len(cs), hprc_freq=round(len(cs) / nC, 4) if nC else '',
                          **{f'{s_}_carriers': sum(sp(k) == s_ for k in cs) for s_ in SUPER},
                          hprc_any_row_carriers=sum(any(has(w, e) for w in anyw[k]) for k in X['HAPS']),
                          min_node_coverage=min((sum(n in anyn.get(k, ()) for k in X['HAPS'] + ['CHM13#0', 'GRCh38#0'])
                                                 for n, _ in toks(e['nodes'])), default=''),
                          CHM13=ch.get('el', {}).get(e['id'], ch['status']), carrier_haps=','.join(cs)))
    mem = []
    for k in X['HAPS'] + ['CHM13#0']:
        s_, h = k.split('#'); m = X['META'][s_]; hk = H[k]
        pop = (m['population_code'], m['superpopulation'] if s_ != 'CHM13' else 'CHM13', hk['status'])
        for e in keep:
            mem.append((t, e['id'], s_, h) + pop + (hk['el'][e['id']] if 'el' in hk else '', '', '', ''))
        mem.append((t, 'allele', s_, h) + pop + tuple(hk.get(x, '') if keep else '' for x in ('joint', 'allele', 'anyp', 'exact')))
    return out, nodes, mem


INTERP = {
    'hapX_only': 'the ALT allele over the whole tandem array is a COLO829BL germline allele (hapX); a tumor-only caller sees a '
                 'germline allele: germline filtering by the graph',
    'hapY_only': 'the ALT allele over the whole tandem array is a COLO829BL germline allele (hapY); a tumor-only caller sees a '
                 'germline allele: germline filtering by the graph',
    'hapX_only:path_allele': 'the d9 ALT path (GRCh38 + truth + the germline alleles it needs) spells over the whole tandem array '
                             'the COLO829BL hapX allele; the perfectly aligned reads carry the patient\'s own germline allele: germline '
                             'filtering by the graph (patient frame)',
    'hapY_only:path_allele': 'the d9 ALT path (GRCh38 + truth + the germline alleles it needs) spells over the whole tandem array '
                             'the COLO829BL hapY allele; the perfectly aligned reads carry the patient\'s own germline allele: germline '
                             'filtering by the graph (patient frame)',
    'exact_allele': 'the patient does not carry it; >= 1 complete HPRC haplotype spells the exact ALT allele over the array '
                    '(somatic slippage / substitution recreating a population allele kept by the d9 graph): pangenome-induced false negative',
    'element_set_other_allele': 'the patient does not carry it; HPRC haplotypes walk all somatic elements but spell another '
                                'allele over the array (e.g. another repeat length): the graph spells the somatic allele from '
                                'population nodes, not a population allele: pangenome-induced false negative',
    'recombinant_pieces': 'the patient does not carry it; every somatic element is walked by some HPRC haplotype but none walks a '
                          'whole ALT path: the ALT path is a recombination of population nodes: pangenome-induced false negative',
}
WHY = {
    'closest:graph_allele_not_truth': 'no d9 path spells the truth; the closest graph allele is nearer GRCh38 / the patient\'s germline than the truth (no read data to say what the reads walk)',
    'closest:graph_allele_nearer_truth': 'no d9 path spells the truth; the closest graph allele is nearer the truth, but not the truth (no read data)',
    'closest:no_element': 'no d9 path spells the truth; the closest graph path has no somatic element',
    'only_germline_elements': 'every element of the path alone spells a COLO829BL germline allele',
    'no_somatic_element': 'the ALT path has no non-GRCh38 element',
    'with_germline:germline_allele_off_event_hap': 'the ALT path needs a germline allele that dipcall phases to the other hap than the c3 event hap',
    'HPRC_no_complete_haplotype': 'no HPRC haplotype traverses the window completely',
    'evidence_conflict:ALT_on_both_normal_haps': 'both COLO829BL haps carry the GRCh38 ALT over the array: homozygous germline, truth conflict',
    'normal_hap_unresolved:one_hap_NA': 'one COLO829BL hap has no array sequence, the other does not carry the ALT',
    'normal_hap_unresolved:both_NA': 'no COLO829BL hap has an array sequence',
    'no_HPRC_carrier': 'no complete HPRC haplotype walks the somatic elements',
    'no_HPRC_carrier:CHM13_carries': 'no complete HPRC haplotype walks the somatic elements; CHM13 does',
}


def interpret(cat, sub, why):
    if cat.startswith('normal_present'):
        return 'COLO829BL germline allele already a graph path (' + ('rare' if cat.endswith('rare') else 'common') + \
            ' exact allele in HPRC): ' + INTERP.get(sub.split('+')[0].split(':')[0] + (':path_allele' if sub.endswith(':path_allele') else ''), sub)
    if cat == 'normal_absent_HPRC_other':
        return INTERP[sub]
    return 'unresolved: ' + WHY.get(why, why)


def md(title, header, rows, note=''):
    s = f'\n### {title}\n\n' + (note + '\n\n' if note else '') + '| ' + ' | '.join(header) + ' |\n|' + '---|' * len(header) + '\n'
    return s + ''.join('| ' + ' | '.join(str(x) for x in r) + ' |\n' for r in rows)


def oe(rows, nm):
    """observed / expected carrier haplotypes per superpopulation (expected per locus = carriers x the superpopulation's
    share of the complete haplotypes; hypergeometric variance) -> {SP: (O, E, z)}."""
    O, E, V = collections.Counter(), collections.Counter(), collections.Counter()
    for r in rows:
        n, tot = int(r['hprc_n_complete']), sum(int(r[f'{s}_{nm}_carriers']) for s in SUPER)
        for s in SUPER:
            q = int(r[f'{s}_complete']) / n
            O[s] += int(r[f'{s}_{nm}_carriers']); E[s] += tot * q
            V[s] += tot * q * (1 - q) * (n - tot) / (n - 1) if n > 1 else 0
    return {s: (O[s], E[s], (O[s] - E[s]) / math.sqrt(V[s]) if V[s] else None) for s in SUPER}


def oe_ci(rows, nm, B=2000, seed=20261005):
    """locus bootstrap of oe(): {SP: (2.5 %, 97.5 %) of O/E} over B resamples of the loci (with replacement)."""
    per = [oe([r], nm) for r in rows]
    rng, res = random.Random(seed), {s: [] for s in SUPER}
    for _ in range(B if rows else 0):
        pick = [per[rng.randrange(len(per))] for _ in per]
        for s in SUPER:
            e = sum(p[s][1] for p in pick)
            if e:
                res[s].append(sum(p[s][0] for p in pick) / e)
    q = lambda v, f: sorted(v)[min(len(v) - 1, int(f * len(v)))]
    return {s: (q(v, 0.025), q(v, 0.975)) if v else None for s, v in res.items()}


def tables(rows, nodes, X):
    SC = [('chr1-22', lambda r: True), ('chr1', lambda r: r['chrom'] == 'chr1'), ('chr2-22', lambda r: r['chrom'] != 'chr1')]
    SETS = [(p, lambda r, p=p: r[f'{p}_perfect'] == 'True') for p in PLATS] + [('union', lambda r: True)]
    K2 = ('SNV', 'INDEL')
    out, summ = ['# c7: HPRC membership and cross-evaluation of the graph-absorbed COLO829T somatic truth alleles (d9 HPRC v1.1)\n',
                 f'RARE_AF = {RARE_AF} on the exact-allele HPRC frequency (complete haplotypes spelling the ALT allele over the '
                 'tandem array). COLO829BL calls: c3 whole-array alleles. No read data.\n'], {}
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
    for scn, scf in SC[:2]:
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
                    ('exact_allele_freq_all88', 'exact allele over all 88 HPRC haplotypes (lower bound)'),
                    ('set_freq_any_path', 'set level, any ALT path (walks all somatic elements)'), ('elem_freq_min', 'element level (min over elements)'),
                    ('d9_vcf_AF', 'd9 VCF AF (called haplotypes, matching record)'), ('full_vcf_AF', 'full-graph VCF AF (all called haplotypes, matching record)')):
        body = []
        num = dict(exact_allele_freq='exact_allele_carriers', exact_allele_freq_all88='exact_allele_carriers', set_freq_any_path='set_carriers_any_path',
                   elem_freq_min='elem_min_carriers').get(col)      # medians from the counts (the stored frequencies are rounded)
        val = lambda r: (float(r[col]) if num is None else int(r[num]) / (88 if col.endswith('all88') else int(r['hprc_n_complete'])))
        for c in CATS:
            for kd in K2:
                f = [val(r) for r in sel(lambda r: r['category'] == c, kd) if r[col] != '']
                body.append([c, kd, len(f)] + [sum(g(x) for x in f) for _, g in fb] + [round(statistics.median(f), 3) if f else ''])
        out.append(md(f'HPRC frequency per category, {nm} (union, chr1-22)', ['category', 'kind', 'n'] + [b for b, _ in fb] + ['median'], body))
    body = []
    for c in CATS + ('all',):
        x = [r for r in rows if (c == 'all' or r['category'] == c) and r['n_somatic_elements'] and r['hprc_n_complete']]
        for nm, lab in (('exact', 'exact allele'), ('anypath', 'set, any ALT path')):
            o, ci = oe(x, nm), oe_ci(x, nm)
            body.append([c, lab, len(x), sum(v[0] for v in o.values())] +
                        [f"{o[s][0] / o[s][1]:.3f} (z {o[s][2]:+.1f}; {ci[s][0]:.2f}-{ci[s][1]:.2f})" if o[s][1] and o[s][2] is not None and ci[s] else ''
                         for s in SUPER] +
                        [sum(r['CHM13_' + ('exact' if nm == 'exact' else 'set_any_path')] is True for r in x)])
    out.append(md('Superpopulation of HPRC carriers per category: observed / expected (union, chr1-22)',
                  ['category', 'level', 'loci', 'carrier haps'] + [f'{s} O/E (z; 95% CI)' for s in SUPER] + ['CHM13 carries'], body,
                  'expected per locus = carriers x the superpopulation\'s share of the complete haplotypes at that locus (hypergeometric variance '
                  'for z; ignores the pairing of haplotypes within individuals and linkage between loci, so z overstates significance); '
                  'after z: 95 % interval of O/E from 2,000 locus bootstrap resamples (seed 20261005), which keeps the loci\'s carriers together. '
                  'COLO829 donor: European (white male); HPRC v1.1 has no EUR sample (AFR 23, AMR 16, EAS 4, SAS 1).'))
    for col, nm in (('hprc_support', 'HPRC support'), ('exact_freq_class', 'exact-allele frequency class'), ('elem_freq_class', 'element-level frequency class'),
                    ('freq_class_agree', 'frequency class agrees across element / set / exact / full-VCF'),
                    ('normal_pattern', 'COLO829BL hap carrying the ALT (c3 array allele, or the with_germline path allele)'),
                    ('normal_carries_via', 'COLO829BL carries: the ALT (GRCh38 + truth) / the with_germline ALT path allele'),
                    ('multicopy_psv_like', 'every exact carrier has a GRCh38 copy too (paralog-like)'), ('normal_allele_class', 'COLO829BL allele class (locus)'),
                    ('event_hap', 'c3 event hap'), ('germline_like_alt_len', 'a COLO829BL hap has the ALT length, other bases (no hap ALT)'),
                    ('germline_phase', 'with_germline: germline alleles on the event hap'), ('path_choice', 'path choice'),
                    ('element_source', 'element source'), ('closest_kind', 'closest: the graph allele'), ('match', 'c2 match'),
                    ('allele_window_clipped', 'tandem array reaches beyond the graph window')):
        vals = sorted({r[col] for r in rows}, key=str)
        out.append(md(f'{nm} x category (union, chr1-22)', [col] + [f'{c} {kd}' for c in CATS for kd in K2],
                      [[vv] + [sum(1 for r in rows if r[col] == vv and r['category'] == c and r['kind2'] == kd) for c in CATS for kd in K2] for vv in vals]))
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
    out.append(md('HPRC carriers: GFA allele level vs VCF carriers_any (on GFA-complete haplotypes; loci with somatic elements)',
                  ['VCF', 'kind', 'loci', 'with matching record', 'identical carrier sets', 'identical, with record', 'hap-level disagreements'], body))
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
                 n_elements=len(nodes), branch_elements=len(cov), branch_elements_node_coverage_min=min(cov) if cov else None,
                 node_coverage_lt7=sum(c < 7 for c in cov), node_coverage_lt9=sum(c < 9 for c in cov),
                 node_coverage_hist=dict(sorted(collections.Counter(min(c, 15) for c in cov).items())),
                 complete_carriers_min=min(ce), complete_carriers_lt9=sum(c < 9 for c in ce), any_row_carriers_min=min(ae), any_row_carriers_lt9=sum(c < 9 for c in ae))
    out.append(md('d9 frequency-filter floor (somatic elements)', ['measure', 'value'], [[k, v] for k, v in floor.items()],
                  'node coverage = HPRC + CHM13 + GRCh38 haplotypes whose c5 rows (complete or partial path) visit the branch node (min over the nodes of the element); a lower bound '
                  '(c5 partial paths stop 500 nodes from an anchor). full_window = loci where all 88 HPRC haplotypes have a complete traversal, where the counts are exact.'))
    return ''.join(out), summ, floor


def work(i):
    return locus(X, X['L'][i])


K = 16


def kcheck(i):
    """sequence-level exact-allele carriers of one locus (independent of the node bracket): anchors = the nearest 16-mers
    of the window left of allele_lo0 / right of allele_hi0 that occur once in GRCh38 (ref_hap), once in GRCh38 + truth
    (alt_hap) and, with_germline, once in the ALT path's sequence; a complete walk is anchored when both occur once in its
    spelled window, in order; carrier = its sequence between them equals alt_hap's (or the ALT path's for with_germline).
    A haplotype is compared only when all its complete traversals are anchored (a SNP inside an anchor hides the allele)."""
    L, r = CK['L'][i], CK['R'][CK['L'][i]['truth_id']]
    t, R, Al, b0 = L['truth_id'], L['ref_hap'], L['alt_hap'], int(L['win_start0'])
    if not r['somatic_elements'] or not int(r['hprc_n_complete']):
        return None
    lo, hi, dl = int(r['allele_lo0']) - b0, int(r['allele_hi0']) - b0, len(Al) - len(R)
    Ps = spell(toks(r['alt_path'])) if r['match'] == 'with_germline' else None
    once = lambda s, m: s.count(m) == 1 and s.find(m) == s.rfind(m)
    ok = lambda m: once(R, m) and once(Al, m) and (Ps is None or once(Ps, m))
    a = next((x for x in range(lo - K, -1, -1) if ok(R[x:x + K])), None)
    b = next((x for x in range(hi, len(R) - K + 1) if ok(R[x:x + K])), None)
    if a is None or b is None:
        return dict(truth_id=t, status='no_anchor')
    la, rb = R[a:a + K], R[b:b + K]
    want = {Al[Al.find(la) + K:Al.find(rb)]} | ({Ps[Ps.find(la) + K:Ps.find(rb)]} if Ps else set())
    seq_c, unanch = set(), set()
    for k, ws in CK['W'][t].items():
        hit = False; anch = True                     # anchored = every complete traversal of the hap is anchored
        for w in ws:
            s = spell(w)
            if once(s, la) and once(s, rb) and s.find(la) + K <= s.find(rb):
                hit |= s[s.find(la) + K:s.find(rb)] in want
            else:
                anch = False
        if hit:
            seq_c.add(k)
        if not anch:
            unanch.add(k)
    c7 = set(filter(None, r['exact_allele_carrier_haps'].split(',')))
    return dict(truth_id=t, status='ok', anchor_left0=b0 + a, anchor_right0=b0 + b, n_complete=len(CK['W'][t]), n_unanchored=len(unanch),
                c7_exact=len(c7), kmer_exact=len(seq_c), agree=(c7 - unanch) == (seq_c - unanch), c7_only=','.join(sorted(c7 - seq_c)),
                kmer_only=','.join(sorted(seq_c - c7)), c7_only_unanchored=len((c7 - seq_c) & unanch))


def check(nproc):
    """independent re-derivation of the exact-allele carriers by sequence (kcheck) -> $D/hprc_check_kmer.tsv + summary."""
    global CK
    CK = dict(L=rd(f'{D}/loci.tsv'), R={r['truth_id']: r for r in rd(f'{D}/hprc_per_variant.tsv')}, W=collections.defaultdict(lambda: collections.defaultdict(list)))
    for r in rd(f'{D}/hprc_local_paths.tsv.gz'):
        if r['status'] == 'complete' and r['sample'] not in ('GRCh38', 'CHM13'):
            CK['W'][r['truth_id']][f"{r['sample']}#{r['hap']}"].append(toks(r['local_path']))
    with Pool(nproc, initializer=S2.init) as pool:
        res = [x for x in pool.map(kcheck, range(len(CK['L'])), chunksize=4) if x]
    cols = ['truth_id', 'status', 'anchor_left0', 'anchor_right0', 'n_complete', 'n_unanchored', 'c7_exact', 'kmer_exact', 'agree',
            'c7_only', 'kmer_only', 'c7_only_unanchored']
    with open(f'{D}/hprc_check_kmer.tsv', 'w') as w:
        w.write('\t'.join(cols) + '\n')
        for x in res:
            w.write('\t'.join(str(x.get(c, '')) for c in cols) + '\n')
    ok = [x for x in res if x['status'] == 'ok']
    print('loci checked', len(res), '| no anchor', sum(x['status'] == 'no_anchor' for x in res), '| anchored loci', len(ok),
          '| identical carrier sets (on anchored haps)', sum(x['agree'] for x in ok),
          '| locus-hap pairs', sum(x['n_complete'] for x in ok), 'unanchored', sum(x['n_unanchored'] for x in ok),
          '| carriers c7', sum(x['c7_exact'] for x in ok), 'kmer', sum(x['kmer_exact'] for x in ok),
          '| c7-only haps', sum(len(x['c7_only'].split(',')) for x in ok if x['c7_only']), '(unanchored', sum(x['c7_only_unanchored'] for x in ok), ')',
          '| kmer-only haps', sum(len(x['kmer_only'].split(',')) for x in ok if x['kmer_only']))
    print('disagreeing loci:', [(x['truth_id'], x['c7_exact'], x['kmer_exact'], x['c7_only'][:80], x['kmer_only'][:80]) for x in ok if not x['agree']])


def main(nproc, test=None):
    global X
    X = load()
    if test:
        X['L'] = [L for L in X['L'] if L['truth_id'] in test]
    print('loaded', len(X['L']), 'loci', flush=True)
    rows, nodes, mem = [], [], []
    with Pool(nproc, initializer=S2.init) as pool:                 # one sqlite connection per worker (c2 init)
        for i, (r, n, m) in enumerate(pool.imap(work, range(len(X['L'])), chunksize=2)):
            rows.append(r); nodes += n; mem += m
            if i % 100 == 0:
                print(i, len(X['L']), flush=True)
    if test:
        for r in rows:
            print('\n'.join(f'  {k}: {str(v)[:300]}' for k, v in r.items() if k not in X['V'][r['truth_id']] or k == 'truth_id')); print()
        return
    cols = list(X['V'][X['L'][0]['truth_id']].keys())
    for r in rows:
        cols += [c for c in r if c not in cols]
    with open(f'{D}/hprc_per_variant.tsv', 'w') as w:
        w.write('\t'.join(cols) + '\n')
        for r in rows:
            w.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
    with open(f'{D}/hprc_per_node.tsv', 'w') as w:
        nc = list(nodes[0].keys()); w.write('\t'.join(nc) + '\n')
        for n in nodes:
            w.write('\t'.join(str(n[c]) for c in nc) + '\n')
    with gzip.open(f'{D}/hprc_membership.tsv.gz', 'wt') as w:
        w.write('truth_id\telement_id\tsample\thap\tpopulation_code\tsuperpopulation\ttraversal_status\tcarries_element\tcarries_alt_allele\t'
                'carries_set_any_path\tcarries_exact_allele\n')
        for m in mem:
            w.write('\t'.join(map(str, m)) + '\n')
    txt, summ, floor = tables(rows, nodes, X)
    open(f'{D}/hprc_tables.md', 'w').write(txt)
    perfect = {p: {kd: sum(r[f'{p}_perfect'] == 'True' and r['kind2'] == kd for r in rows) for kd in ('SNV', 'INDEL')} for p in PLATS}
    for p in PLATS:
        for kd in ('SNV', 'INDEL'):
            assert summ[f'{p} {kd} chr1-22']['total'] == perfect[p][kd] == sum(summ[f'{p} {kd} chr1-22'].get(c, 0) for c in CATS)
    wc = collections.Counter()
    for r in rows:
        for d in filter(None, r['window_checks'].split(';')):
            a, n = d.split(':'); wc[a] += int(n)
    S = dict(n_loci=len(rows), RARE_AF=RARE_AF, perfect=perfect, categories=summ, floor=floor, window_checks=dict(wc),
             vcf_agree={so: {kd: sum(r.get(f'{so}_agree') is True for r in rows if r['kind2'] == kd) for kd in ('SNV', 'INDEL')} for so in ('d9', 'full')},
             path_choice=dict(collections.Counter(r['path_choice'] for r in rows)),
             exact_core_only={'haplotypes': sum(int(r['exact_core_only'] or 0) for r in rows), 'loci': sorted((r['truth_id'], r['exact_core_only']) for r in rows if r['exact_core_only'])},
             exact_window_only={'haplotypes': sum(int(r['exact_window_only'] or 0) for r in rows), 'loci': sorted((r['truth_id'], r['exact_window_only']) for r in rows if r['exact_window_only'])},
             multicopy_psv_like=sorted((r['truth_id'], r['category'], r['exact_allele_carriers'], r['exact_carriers_alt_ref_copy']) for r in rows if r['multicopy_psv_like']),
             normal_path_allele=sorted((r['truth_id'], r['category'], r['sub_class']) for r in rows if r['normal_carries_via'] == 'path_allele'),
             bypass={'loci': sorted((r['truth_id'], int(X['COV'][r['truth_id']]['n_bypass'])) for r in rows if int(X['COV'][r['truth_id']]['n_bypass']))},
             sub_class_or_reason={f'{a}|{b}|{k}': n for (a, b, k), n in sorted(collections.Counter((r['category'], r['sub_class'] or r['sub_reason'], r['kind2'])
                                                                                                 for r in rows).items())},
             germline_like_alt_len={f'{a}|{k}': n for (a, k), n in sorted(collections.Counter((r['category'], r['kind2']) for r in rows if r['germline_like_alt_len']).items())},
             freq_class_agree_present=sum(r['freq_class_agree'] is True for r in rows if r['category'].startswith('normal_present')))
    json.dump(S, open(f'{D}/hprc_summary.json', 'w'), indent=1, default=str)
    rng, lines = random.Random(9), []
    for c in CATS:
        x = [r for r in rows if r['category'] == c]
        pick = (rng.sample([r for r in x if r['kind2'] == 'SNV'], min(1, sum(r['kind2'] == 'SNV' for r in x))) +
                rng.sample([r for r in x if r['kind2'] == 'INDEL'], min(3, sum(r['kind2'] == 'INDEL' for r in x))))[:3]
        for r in pick:
            lines.append(f"== {c} truth {r['truth_id']} {r['chrom']}:{r['vcf_pos']} {r['vcf_ref']}>{r['vcf_alt']} {r['kind2']} perfect F/O/I "
                         f"{r['fiberseq_perfect']}/{r['ONT_perfect']}/{r['Illumina_perfect']}")
            for k in ('match', 'element_source', 'path_choice', 'closest_kind', 'somatic_elements', 'element_alone', 'dropped_germline_elements', 'alt_path',
                      'hprc_n_complete', 'elem_carriers_each', 'elem_freq_min', 'set_carriers_any_path', 'exact_allele_carriers', 'exact_allele_freq',
                      'hprc_support', 'exact_allele_carrier_haps', 'd9_n_carriers_any', 'd9_agree', 'CHM13_exact', 'array_lo0', 'array_hi0',
                      'hapX_allele', 'hapY_allele', 'event_hap', 'pf_label', 'pf_event', 'alt_len_haps', 'normal_pattern', 'sub_class', 'sub_reason'):
                lines.append(f'  {k}: {str(r.get(k, ""))[:300]}')
    open(f'{D}/hprc_spot_checks.txt', 'w').write('\n'.join(lines) + '\n')
    print(json.dumps({k: S[k] for k in ('perfect', 'path_choice', 'vcf_agree', 'window_checks', 'germline_like_alt_len', 'exact_core_only', 'exact_window_only',
                                         'multicopy_psv_like', 'normal_path_allele', 'bypass')}, default=str))
    print('window-check violations:', [(r['truth_id'], r['window_checks']) for r in rows if 'VIOLATION' in r['window_checks']])
    for k in [f'{p} {kd} {sc}' for sc in ('chr1-22', 'chr1') for p in PLATS + ('union',) for kd in ('SNV', 'INDEL')]:
        print(k, summ[k])
    print(txt)


if __name__ == '__main__':
    if sys.argv[1:2] == ['check']:
        check(int(sys.argv[2]) if sys.argv[2:] else 4); sys.exit()
    main(int(sys.argv[1]) if sys.argv[1:] else int(os.environ.get('SLURM_CPUS_PER_TASK', 4)), set(sys.argv[2].split(',')) if sys.argv[2:] else None)

# Results (2026-10-05, after the HPRC audit fixes: c5 bypass rule, core rule, support order, path-allele normal rule;
# login node, $T/c7_hprc_membership.login.log: 149 s with 6 processes, 0.15 GB per process; check 40 s,
# $T/c7_hprc_check_kmer.login.log). 497 truths: path choice single 438, c2 primary of several 34, phase-consistent
# alternative 1, closest 24. Categories (union, chr1-22, SNV / INDEL): normal_present_broad 0 / 25 (hapX 23, hapY 1,
# hapX:path_allele 1 = 2562, frequency 14 / 69 = 0.203), normal_present_rare 1 / 8, normal_absent_HPRC_other 63 / 367
# (exact_allele 58 / 332, element_set_other_allele 5 / 33, recombinant_pieces 0 / 2), other_ambiguous 3 / 30 (closest 3 /
# 21, no_HPRC_carrier 0 / 6 of which CHM13 only 0 / 2, off-phase germline 0 / 2, normal hap unresolved 0 / 1). chr1:
# present 0 / 14, absent 4 / 22, ambiguous 1 / 2. Changed vs the first version: 2562 ambiguous -> present_broad
# (path allele), 2570 / 38059 ambiguous -> absent element_set_other_allele, 29389 ambiguous -> absent exact_allele (31 /
# 79), 21200 element_set_other_allele -> exact_allele (10 / 74); carriers / complete changed at 2133 (6 -> 19 / 77),
# 21808, 22598, 22742, 28014, 29562, 29579, 37389, 39443, 41860. Core rule: +54 carriers (21200 10, 2133 13, 29389 31);
# whole-window rule: 1 (29562). PSV-like: 35165, 35166 (every exact carrier has an ALT and a REF copy). Median
# exact-allele frequency, absent SNV / INDEL 0.188 / 0.159 (0.205 / 0.171 at loci with a carrier), present 0.316.
# O/E (absent, exact allele) AFR 1.17 (z +17.4; locus bootstrap 1.14-1.21), AMR 0.83, EAS 0.80, SAS 0.86; present_broad
# AFR 1.08 (0.96-1.19), AMR 0.90 (0.77-1.04). VCF carriers_any = GFA allele carriers: d9 SNV 52 / 65, INDEL 308 / 427;
# full 52 / 65, 312 / 427. Window checks: 17,311 GRCh38 walks carry nothing; all 6,297 GRCh38 + truth walks are exact
# carriers. k-mer check: identical carrier sets at 474 / 475 anchored loci (34,215 locus-haplotype pairs; 17 loci without
# unique anchors); 29562: c7 11, sequence 14 (the 4 carry a germline +4 of the next (TTTC)n repeat, which their walk
# places at the allele-window end). d9 floor: 12 of 413 branch elements have < 9 visiting haplotypes in the window rows, 10 on CHM13; 3346 (8)
# and 37542 (4) are window lower bounds. The first version's outputs are in $D/prefix_hprc_20261005/.
