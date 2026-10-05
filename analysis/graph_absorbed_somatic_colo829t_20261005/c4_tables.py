"""COLO829T step 4: repeat-context tables, the per-variant file and the README of the absorbed (perfect-bypass) truths,
with the HG008T side-by-side, and the README sections 7-8 (HPRC haplotypes / individuals / populations, cross-evaluation
with the normal) from the c5-c8 outputs. Run it last (after c8).

Inputs (all read only; $D = /scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005):
  $D/variant_set.tsv      c0: the 497 absorbed truths (per platform: perfect flag, status, miss class, read-level reason)
  $D/grch38_context.tsv   c1: GRCh38-frame context of every chr1-22 truth allele (43,947; b1 rules; b5 rule in GRCh38)
  $D/loci.tsv             c1b: d9 graph window and event window per absorbed truth
  $D/graph_paths.tsv, $D/graph_elements.tsv   c2: does a d9 path spell the truth (exact / with_germline / closest), elements
  $D/normal_frame.tsv     c3: COLO829BL hapX / hapY array alleles and the patient-frame event (497 absorbed + all 1,912
                          chr1-22 truth INDELs)
  $D/prefix_20261005/normal_frame.tsv   c3 before the audit fixes (frozen copy; only for the 'label changed' row of 2h)
  analysis/tensor_recall_20260930/per_truth_COLO829T.tsv (no-candidate counts per platform)
  HG008T (side-by-side only): analysis/graph_absorbed_somatic_20261001/{per_variant.tsv, indel_repeat_context.tsv},
  analysis/tensor_recall_20260930/per_truth_HG008T.tsv, $H/audit/biology/{b1_truth_context.tsv, b5_normal_frame.tsv}
  ($H = /scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001); HG008 classes recomputed with the s14 /
  s11 functions and checked against the published indel_repeat_context.md numbers (asserts).
  HPRC (sections 7-8): hprc_per_variant.tsv, hprc_per_node.tsv, hprc_individuals.tsv, hprc_populations.tsv,
  per_locus_populations.tsv (c8, this folder); $D/hprc_summary.json, $D/hprc_check_kmer.tsv (c7), $D/hprc_window_coverage.tsv
  (c5), $D/hprc_columns.tsv (c8, appended to columns.tsv); c8_populations.py imported for oe(); HG008T per_locus_populations /
  hprc_individuals.tsv; analysis/tensor_recall_20260930/pon/truth_pon_COLO829T.tsv (PoN tags, gnomAD / CoLoRSdb AF).
  HG008T comparison numbers recomputed with the same functions and checked against its README / tables.md (asserts).
Classes
  patient-frame label (c3 pf_label) of a truth: the somatic change from the COLO829BL hap closest to GRCh38 + truth over the
  tandem array: HP>=7 / HP4-6 / STR2-6>=3copies / VNTR>6 / in_STR / imperfect VNTR (INDEL event); 'SNV in <class>' / SNV
  (substitution); complex; normal carries ALT (a COLO829BL hap already has GRCh38 + truth over the array); single hap (one
  hap only, not the ALT); NA (no hap sequence).
  group: HP>=7 | STR | other repeat (HP4-6, VNTR>6, in_STR, imperfect VNTR) | SNV in repeat (SNV in HP>=7 / HP4-6 / STR /
  imperfect VNTR) | non-repeat (none, SNV) | normal carries ALT | complex | unresolved (NA, single hap; HG008: no INFO event).
  repeat-unit changes = HP>=7, STR, HP4-6, VNTR>6 events; in or next to a repeat = in_STR, imperfect VNTR.
  label_b5 = the label with pf_class_b5 (audit b5 rule verbatim, the HG008 rule; no imperfect-VNTR relabel).
  ALT-length hap (c3 alt_len_haps): a hap with exactly the ALT length that is not the ALT (germline-like reading).
  units changed = |pf_len_change| / pf_unit for HP>=7 / HP4-6 / STR / VNTR>6 events (1 unit / 2-3 units / >3 units).
  HP tract bin = the normal's homopolymer length (pf_tract) of HP>=7 events: 7-9 / 10-14 / 15-19 / 20-29 / >=30.
  germline status (normal frame, array allele of hapX / hapY): normal carries ALT > germline other length (a hap's array
  has another length) > germline same length, other bases > both REF > one hap NA, other REF > both NA.
  mechanism class (background rates, INDEL records): 'HP>=7 1-unit, tract <bin>', 'STR 1-unit', 'repeat multi-unit'
  (HP>=7 / STR with >= 2 units), HP4-6, VNTR>6, in_STR, imperfect VNTR, SNV in repeat, non-repeat INDEL / SNV, complex,
  normal carries ALT, unresolved. COLO829T: c3 label + c3 tract (mechanism_class) and b5 label + b5 tract
  (mechanism_class_b5); HG008T: the same function on its b5 table. Every COLO829T vs HG008T rate / share uses b5 on both.
  VAF bin (truth VAF_Ill): <0.1, 0.1-0.25, 0.25-0.4, >=0.4. RGN = SMaHT region (Easy / Difficult / Extreme).
Outputs (this folder): per_variant.tsv (one row per absorbed truth, 497), columns.tsv (meaning and source of every column,
the HPRC files included), repeat_context.md (the repeat-context tables), README.md (report; every number filled in from this
script; tables h1-h9 of sections 7-8 are printed, not written to repeat_context.md); everything printed.
Assumptions: per-set tables count truth alleles (one GRCh38 record = one row, as in HG008); the union set = perfect on >= 1
platform; INDEL / SNV = the GRCh38 record kind (kind2), so 'SNV in repeat' rows of an INDEL table are GRCh38 INDEL records
whose patient-frame event is a substitution.
"""
import collections, csv, importlib.util, json, re, statistics, sys
from pathlib import Path
C = collections.Counter
csv.field_size_limit(10 ** 9)
sys.dont_write_bytecode = True
A = Path(__file__).resolve().parent
D = Path('/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005')
TR = Path('/scratch/jshen/Github/Pansoma/analysis/tensor_recall_20260930')
HA = Path('/scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_20261001')
HB = Path('/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001/audit/biology')
PLAT = ('fiberseq', 'ONT', 'Illumina')
rd = lambda p: list(csv.DictReader(open(p), delimiter='\t'))
byid = lambda rows: {r['truth_id']: r for r in rows}

# ---------------------------------------------------------------- inputs + consistency checks
vs = rd(D / 'variant_set.tsv'); VS = byid(vs); ABS = set(VS)
g38 = byid(rd(D / 'grch38_context.tsv')); loc = byid(rd(D / 'loci.tsv')); gp = byid(rd(D / 'graph_paths.tsv'))
nf = byid(rd(D / 'normal_frame.tsv')); pt = rd(TR / 'per_truth_COLO829T.tsv')
PRE = {t: r['pf_label'] for t, r in byid(rd(D / 'prefix_20261005' / 'normal_frame.tsv')).items()}   # c3 before the audit fixes
gel = collections.defaultdict(list)
for r in rd(D / 'graph_elements.tsv'):
    gel[r['truth_id']].append(r)
assert len(vs) == 497 and set(loc) == ABS and set(gp) == ABS and set(gel) <= ABS
assert {t for t, r in nf.items() if r['absorbed'] == 'True'} == ABS and {t for t, r in g38.items() if r['absorbed'] == 'True'} == ABS
assert len(g38) == len(pt) == 43947 and set(g38) == {r['truth_id'] for r in pt}
assert all(r['in_bed'] == 'True' and r['passed'] == 'True' for r in pt)          # README: all in the BED, all PASS
_ph = collections.defaultdict(set)
for r in rd(D / 'c3_dipcall_phase.tsv'):
    _ph[(r['pos'], r['ref'], r['alt'])].add(r['verdict'])
PHASE = C(next(iter(v)) if len(v) == 1 else 'conflicting' for v in _ph.values())       # per distinct SNV (a SNV near several loci once)
INDELS = [t for t, r in g38.items() if r['kind2'] == 'INDEL']
assert len(INDELS) == 1912 and set(INDELS) <= set(nf) and len(nf) == 1912 + sum(VS[t]['kind2'] == 'SNV' for t in ABS)
for t in ABS:
    assert VS[t]['kind2'] == g38[t]['kind2'] == nf[t]['kind2'] == loc[t]['kind2'] == gp[t]['kind2'], t
    assert (VS[t]['chrom'], VS[t]['vcf_pos'], VS[t]['vcf_ref'], VS[t]['vcf_alt']) == \
           (nf[t]['chrom'], nf[t]['vcf_pos'], nf[t]['vcf_ref'], nf[t]['vcf_alt']), t
    assert loc[t]['status'] == 'ok'
    assert {e['id'] for e in json.loads(gp[t]['primary_elements'])} == {e['element_id'] for e in gel[t] if e['in_primary'] == 'True'}, t
perf = lambda t, p: VS[t][f'{p}_perfect'] == 'True'

# ---------------------------------------------------------------- derived per-locus fields
REP = ('HP>=7', 'HP4-6', 'STR2-6>=3copies', 'VNTR>6')
GROUPS = ['HP>=7', 'STR', 'other repeat', 'SNV in repeat', 'non-repeat', 'normal carries ALT', 'complex', 'unresolved']
LABELS = ['HP>=7', 'STR2-6>=3copies', 'HP4-6', 'VNTR>6', 'in_STR', 'imperfect VNTR', 'SNV in HP>=7', 'SNV in STR2-6>=3copies',
          'SNV in HP4-6', 'SNV in imperfect VNTR', 'SNV', 'none', 'complex', 'MNV', 'normal carries ALT', 'single hap', 'NA', 'no event']
HPB = ['7-9', '10-14', '15-19', '20-29', '>=30']
UNB = ['1 unit', '2-3 units', '>3 units']
GERM = ['normal carries ALT', 'germline other length', 'germline same length, other bases', 'both REF', 'one hap NA, other REF',
        'both NA']
VAFB = ['<0.1', '0.1-0.25', '0.25-0.4', '>=0.4']
RGN = ['Easy', 'Difficult', 'Extreme']
G38 = ['HP>=7', 'STR2-6>=3copies', 'HP4-6', 'VNTR>6', 'complex', 'none']


def group(lab):
    return ('HP>=7' if lab == 'HP>=7' else 'STR' if lab == 'STR2-6>=3copies' else 'other repeat' if lab in ('HP4-6', 'VNTR>6', 'in_STR',
            'imperfect VNTR') else 'SNV in repeat' if lab.startswith('SNV in') else 'non-repeat' if lab in ('none', 'SNV')
            else 'unresolved' if lab in ('NA', 'no event', 'single hap', '') else 'complex' if lab in ('complex', 'MNV') else lab)


hbin = lambda L: '7-9' if L < 10 else '10-14' if L < 15 else '15-19' if L < 20 else '20-29' if L < 30 else '>=30'
ubin = lambda u: '1 unit' if u == 1 else '2-3 units' if u <= 3 else '>3 units'
vbin = lambda v: '<0.1' if v < 0.1 else '0.1-0.25' if v < 0.25 else '0.25-0.4' if v < 0.4 else '>=0.4'


def units_bin(n):
    return ubin(abs(int(n['pf_len_change'])) // int(n['pf_unit'])) if n['pf_label'] in REP else ''


def hp_bin(n):
    return hbin(int(n['pf_tract'])) if n['pf_label'] == 'HP>=7' else ''


def label_b5(n):
    """the b5 rule verbatim (HG008): pf_class_b5, without c3's imperfect-VNTR relabel; single-hap / NA / normal carries ALT kept."""
    if n['pf_label'] in ('single hap', 'NA', 'normal carries ALT') or n['pf_kind'] not in ('INDEL', 'SNV'):
        return n['pf_label']
    c = n['pf_class_b5']
    if n['pf_kind'] == 'SNV':
        return 'SNV' if c == 'none' else f'SNV in {c}'
    return 'in_STR' if c.startswith('in_STR') else c


def germ(n):
    a = [n['hapX_allele'], n['hapY_allele']]
    s = ('normal carries ALT' if 'ALT' in a else 'germline other length' if any(x.startswith('germline_len') for x in a)
         else 'germline same length, other bases' if 'germline_seq' in a else 'both REF' if a == ['REF', 'REF']
         else 'both NA' if a == ['NA', 'NA'] else 'one hap NA, other REF')
    assert (s == 'normal carries ALT') == (n['pf_label'] == 'normal carries ALT'), n['truth_id']
    return s


def mech(lab, units, tract):
    """background class of an INDEL record (patient frame); same rule for COLO829T (c3) and HG008T (b5)."""
    if lab in ('HP>=7', 'STR2-6>=3copies'):
        return ('repeat multi-unit' if units != 1 else f'HP>=7 1-unit, tract {hbin(tract)}' if lab == 'HP>=7' else 'STR 1-unit')
    if lab in ('HP4-6', 'VNTR>6', 'in_STR', 'imperfect VNTR', 'complex', 'normal carries ALT'):
        return lab
    return ('SNV in repeat' if lab.startswith('SNV in') else 'non-repeat INDEL' if lab == 'none' else 'non-repeat SNV' if lab == 'SNV'
            else 'complex' if lab == 'MNV' else 'unresolved')


MECH = [f'HP>=7 1-unit, tract {b}' for b in HPB] + ['STR 1-unit', 'repeat multi-unit', 'HP4-6', 'VNTR>6', 'in_STR', 'imperfect VNTR',
                                                     'SNV in repeat', 'non-repeat INDEL', 'non-repeat SNV', 'complex', 'normal carries ALT',
                                                     'unresolved']
cmech = lambda n: mech(n['pf_label'], abs(int(n['pf_len_change'])) // int(n['pf_unit']) if n['pf_kind'] == 'INDEL' else 0,
                       int(n['pf_tract'] or 0))
cmech5 = lambda n: mech(label_b5(n), abs(int(n['pf_len_change'])) // int(n['pf_unit']) if n['pf_kind'] == 'INDEL' else 0,
                        int(n['pf_tract_b5'] or 0))                # the b5 rule (HG008 comparison)
RC = str.maketrans('ACGT', 'TGCA')
spec = lambda r, a: f'{r}>{a}' if r in 'CT' else f'{r.translate(RC)}>{a.translate(RC)}'      # pyrimidine-normalised SNV


def compact(e):
    s, e0 = (e['start0'], e['end0'])
    return f"{e.get('id') or e['element_id']}|{e['subtype']}|{s}-{e0}|{e['seq'] or '-'}"


# ---------------------------------------------------------------- per_variant.tsv
COLS = []          # (column, source, meaning)
def col(name, src, meaning): COLS.append((name, src, meaning))


col('truth_id', 'variant_set', 'row id of the truth allele in analysis/tensor_recall_20260930/per_truth_COLO829T.tsv')
for c, m in (('chrom', ''), ('vcf_pos', 'truth VCF POS (1-based)'), ('vcf_ref', ''), ('vcf_alt', ''),
             ('kind2', 'SNV or INDEL: the GRCh38 record kind of the truth'), ('VAF_Ill', 'truth VAF from Illumina (truth VCF)'),
             ('VAF_PB', 'truth VAF from PacBio (truth VCF)'), ('RGN', 'SMaHT region of the truth: Easy / Difficult / Extreme'),
             ('n_platforms', 'number of platforms on which the truth is perfect bypass')):
    col(c, 'variant_set', m)
col('vaf_bin', 'derived', 'VAF_Ill bin: <0.1 / 0.1-0.25 / 0.25-0.4 / >=0.4')
col('perfect_platforms', 'derived', 'platforms on which the truth is perfect bypass, comma-separated')
for p in PLAT:
    col(f'{p}_perfect', 'variant_set', f'{p}: status no_candidate and read-level reason site_bypassed:{{branch,skip_edge}}:no_edit')
    col(f'{p}_status', 'variant_set', f'{p}: truth status in the 2026-09-30 tensor set (per_truth_COLO829T.tsv)')
    col(f'{p}_class', 'variant_set', f'{p}: miss class of the COLO829T read-level miss analysis (I1-I7 INDEL, S1-S4c SNV)')
    col(f'{p}_reason', 'variant_set', f'{p}: read-level reason (majority of the ALT reads)')
    col(f'{p}_alt_reads', 'variant_set', f'{p}: ALT-exact + ALT-like reads of the miss analysis')
    col(f'{p}_spanning', 'variant_set', f'{p}: reads spanning the site in the miss analysis')
for c, m in (('first_node', 'd9 GRCh38 anchor node at the window start'), ('last_node', 'd9 GRCh38 anchor node at the window end'),
             ('n_ref_nodes', 'GRCh38 nodes in the window'), ('win_start0', 'window start, 0-based GRCh38'),
             ('win_end0', 'window end, 0-based inclusive'),
             ('ev_lo0', 'event window start (period 1-6 stretches touching the event +-5 bp), 0-based'),
             ('ev_hi0', 'event window end, 0-based half-open')):
    col(c, 'loci', m)
for c, m in (('match', 'a d9 path between the anchors spells GRCh38 + truth: exact / with_germline (only with nearby COLO829BL '
                       'dipcall PASS alleles applied too) / closest (no such path; the nearest graph haplotypes)'),
             ('germline_used', 'with_germline: dipcall alleles applied (pos:REF>ALT)'),
             ('germline_gt', 'their dipcall GT (hap1|hap2 = hapY|hapX)'),
             ('closest_edit_distance', 'closest: edit distance of the nearest graph haplotype to the ALT haplotype'),
             ('grch38_edit_distance', 'closest: edit distance of GRCh38 itself to the ALT haplotype (equal = nothing closer)'),
             ('n_alt_paths', 'distinct matching paths kept (<= 50)'), ('cap_hit', 'an enumeration cap was hit'),
             ('cap_info', 'which cap'), ('n_completed', 'closest: paths completed'),
             ('n_event_sets', 'distinct sets of in-event elements over the kept paths'),
             ('primary_n_nodes', 'nodes of the primary path (fewest in-event elements, then fewest nodes)'),
             ('primary_event_ids', 'in-event element ids of the primary path'),
             ('primary_d_grch38', 'edit distance of the primary path to GRCh38'),
             ('primary_local_alt', 'the primary path spells the ALT sequence of the event window +-10 bp'),
             ('primary_local_ref', 'the primary path spells the GRCh38 sequence of the event window +-10 bp'),
             ('truth_net', 'len(ALT) - len(REF)'), ('event_net', 'net length of the primary path\'s in-event elements'),
             ('primary_path', 'primary path, oriented d9 node ids (>n forward, <n reverse)')):
    col(c, 'graph_paths', m)
col('graph_primary_elements', 'graph_paths primary_elements (compact)',
    "every element of the primary path: id|subtype|start0-end0|seq ('B:x:y:>n..' = non-GRCh38 node run between GRCh38 "
    "nodes x and y; 'S:x:y' = edge skipping GRCh38 bases), ';'-separated")
col('graph_n_primary_elements', 'derived', 'elements of the primary path')
col('graph_in_event_subtypes', 'derived', 'subtypes of the primary path\'s in-event elements (ins_branch / del_skip / snv_branch / '
    'mnv_branch / replacing_branch / back_branch / back_edge), sorted, \'+\'-joined')
col('graph_event_element_sets', 'graph_paths event_element_sets (compact)', 'in-event element sets over the kept paths: '
    'elements \'+\'-joined, (n_paths), sets \';\'-separated')
for c, s, m in (('grch38_ins_del', 'ins_del', 'INS / DEL (INDEL records)'),
                ('grch38_indel_seq', 'indel_seq', 'inserted / deleted bases after the shared prefix'),
                ('grch38_unit', 'indel_unit', 'minimal period of the inserted / deleted bases'),
                ('grch38_len_change', 'indel_len', 'len(ALT) - len(REF)'),
                ('grch38_period', 'ctx_period', 'period of the GRCh38 stretch used for the class (0: none found)'),
                ('grch38_tract', 'ctx_tract', 'length of that stretch (bp)'),
                ('grch38_class', 'ctx_class', 'GRCh38-frame class (b1 rule): HP>=7 / HP4-6 / STR2-6>=3copies / VNTR>6 / complex / none'),
                ('grch38_b5rule_class', 'b5rule_class', 'the patient-frame (b5) rule applied in GRCh38 (adds in_STR(pP))'),
                ('grch38_b5rule_period', 'b5rule_period', 'its period'), ('grch38_b5rule_tract', 'b5rule_tract', 'its tract (bp)'),
                ('cpg', 'cpg', 'SNV: the base is the C or G of a GRCh38 CpG'), ('cpg_ti', 'cpg_ti', 'SNV: C>T / G>A at a CpG'),
                ('titv', 'titv', 'SNV: Ti / Tv'), ('alt_hp_len', 'alt_hp_len', 'SNV: homopolymer of the ALT base it creates'),
                ('tri', 'tri', 'SNV: GRCh38 trinucleotide'), ('snv_spectrum', None, 'SNV: pyrimidine-normalised substitution')):
    col(c, f'grch38_context {s}' if s else 'derived', m)
NFC = [c for c in rd(D / 'normal_frame.tsv')[0] if c not in ('truth_id', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind2', 'absorbed',
                                                             'R_seq', 'A_seq', 'hapX_seq', 'hapY_seq', 'flags')]
NF_MEAN = {'arr_lo0': 'tandem array start (0-based GRCh38; the event window grown by every periodic stretch of period 1-60 '
                      'touching it, cap 3 kb)', 'arr_hi0': 'tandem array end (half-open)', 'arr_len': 'array length (bp)',
           'arr_capped': 'array hit the 3 kb cap', 'anchor_a0': 'left common unique 24-mer anchor (its end, 0-based)',
           'anchor_b0': 'right common anchor start (0-based)', 'len_R': 'GRCh38 length between the anchors',
           'event_hap': 'the COLO829BL hap whose array allele is closest to GRCh38 + truth (Levenshtein; tie -> the REF one, '
                        'else hapX)', 'event_basis': 'both / only_<hap> / none: haps with a sequence',
           'event_tie': 'tie_same_seq (both haps identical) / tie_REF_chosen / tie_other / empty',
           'pf_kind': 'patient-frame event kind: INDEL / SNV / complex / normal_carries_ALT / NA',
           'pf_event': 'event in the array frame: offset:REF>ALT (left-normalised)', 'pf_len_change': 'event length change',
           'pf_unit': 'minimal period u of the inserted / deleted bases', 'pf_tract': "the normal's repeat length L (bp): the "
           "exact period-u stretch holding the event's bases", 'pf_tandem': 'a unit of the inserted / deleted bases sits next to the event',
           'pf_stretch': 'longest period 1-6 stretch touching the event (pP:L)',
           'pf_class': 'patient-frame class: HP>=7 / HP4-6 / STR2-6>=3copies / VNTR>6 / in_STR(pP) / imperfect_VNTR / none; SNV '
                       'events: the class of the stretch around the base; complex; normal_carries_ALT; NA',
           'pf_class_b5': 'the same with the audit b5 rule verbatim (the HG008 rule; any stretch touching the event)',
           'pf_label': 'the class label of the tables (INDEL: pf_class with in_STR(pP) -> in_STR, imperfect_VNTR -> imperfect VNTR; '
                       "SNV: SNV in <class> / SNV; 'single hap' when only one hap has a sequence and it is not the ALT)",
           'pf_units_changed': '|pf_len_change| / pf_unit (repeat INDEL events)', 'pf_cpg_ti': 'SNV event: C>T / G>A at a CpG',
           'ivntr_spacing': 'imperfect tandem repeat around the truth (GRCh38 +-600 bp; 12-mers recurring 7-150 bp apart, >= 60 % '
                            'over +-15 bp, run >= 40 bp holding the truth): commonest recurrence distance (unit or a multiple); 0 = none',
           'ivntr_lo0': 'its start (0-based GRCh38)', 'ivntr_hi0': 'its end (half-open)',
           'ivntr_beyond': 'bp of it outside the exact-period array [arr_lo0, arr_hi0)',
           'in_ivntr': 'ivntr_spacing >= 7 and ivntr_beyond >= 20: the truth sits in an imperfect VNTR / minisatellite that the '
                       'exact-period array misses (a none event there is labelled imperfect VNTR)',
           'anchor_wide': 'no common anchor pair within 400 bp for all placed haps, found within 2 kb (wide retry)',
           'pf_tract_b5': 'the tract of the b5 rule (longest period-u stretch merely touching the event; u > 6: periods 1-6)',
           'pf_label_one_hap': "single hap: the label the one available hap gives (pf_label is then 'single hap')",
           'alt_len_haps': 'INDEL truths: haps whose array allele has exactly the ALT length but is not the ALT (the normal may '
                           'already carry the ALT length, differing by substitutions)',
           'pf_titv': 'SNV event: Ti / Tv', 'other_hap_kind': 'the event derived from the other hap: kind',
           'other_hap_event': 'the event from the other hap', 'other_hap_label': 'its label'}
for c in NFC:
    h = c[:4] if c.startswith('hap') else ''
    m = NF_MEAN.get(c) or (f'{h}: ' + {'place': 'placement: cover / cover0 / bridge / empty (none)', 'place_na': 'why not placed',
                                       'unloc_better': 'an unlocalized contig covers the array better (wrong-copy flag)',
                                       'strand': 'scaffold strand of the hit', 'mapq': 'minimap2 MAPQ of the hit', 'tp': 'PAF tp tag',
                                       'n_cover': 'covering hits on the same-chromosome scaffold',
                                       'n_cover_other': 'covering hits on other scaffolds (duplication flag)',
                                       'contig_span': 'scaffold span between the anchors', 'asm_call': 'assembly array allele: '
                                       'ALT / REF / other / NA', 'asm_len_change': 'assembly allele length - len_R',
                                       'asm_na': 'why no assembly allele', 'dip_call': 'dipcall-reconstructed allele: ALT / REF / other / NA',
                                       'dip_len_change': 'dipcall allele length - len_R', 'dip_na': 'why no dipcall allele',
                                       'asm_eq_dip': 'assembly allele == dipcall allele',
                                       'asm_eq_dip_swapped': 'assembly allele == the other GT index (GT-order check)',
                                       'phase_mismatch_ragtag': 'isolated phased het SNVs of dipcall contradicted by the ragtag copy',
                                       'suspect': 'why the ragtag copy was suspect: unplaced / mapq<20 / other_scaffold_hit / '
                                       'phase_mismatch (empty: kept)', 'copy': 'hap copy used: ragtag (same-chromosome scaffold hit) / '
                                       "dipcall (the contig dipcall aligned there; for suspect copies)",
                                       'dpaf': "the dipcall primary alignment used (contig:qs-qe strand, MAPQ); 'none': no such alignment",
                                       'phase_mismatch': 'the same for the copy used',
                                       'total_len_change': 'hap length between the anchors - len_R (before reverting flank differences)',
                                       'flank_len': 'total_len_change - len_change: length of the reverted flank differences',
                                       'src': 'sequence used: asm / dip / NA', 'n_flank_diffs': 'germline differences outside the '
                                       'array (reverted to GRCh38)', 'allele': 'array allele: REF / ALT / germline_len<+-n> (other '
                                       'length) / germline_seq (same length, other bases) / NA', 'len_change': 'array allele '
                                       'length - len_R', 'dist_A': 'Levenshtein distance of the array allele to GRCh38 + truth'}[c[5:]])
    col(c, f'normal_frame {c}', m)
col('nf_flags', 'normal_frame flags', 'per-locus flags: <hap>:NA / dip_only / asm_ne_dip / mapq<20 / other_scaffold_hit / cover0 / '
    'bridge / flank_diffs>=3, tie_REF_chosen / tie_other, other_hap_single_event')
for c, m in (('pf_group', 'group of pf_label: HP>=7 / STR / other repeat / SNV in repeat / non-repeat / normal carries ALT / complex / '
              'unresolved (NA, single hap)'), ('label_b5', 'pf_label with pf_class_b5 for INDEL events (the HG008 rule)'),
             ('units_bin', '1 unit / 2-3 units / >3 units (HP>=7, HP4-6, STR, VNTR>6 events)'),
             ('hp_tract_bin', "the normal's homopolymer length of HP>=7 events: 7-9 / 10-14 / 15-19 / 20-29 / >=30"),
             ('germline_status', 'COLO829BL array alleles: normal carries ALT / germline other length / germline same length, other '
              'bases / both REF / one hap NA, other REF / both NA'),
             ('mechanism_class', "background class (HP>=7 1-unit by tract, STR 1-unit, ...; README 'Definitions')"),
             ('mechanism_class_b5', 'the same with label_b5 and the b5 tract (the HG008 rule; HG008T comparison)')):
    col(c, 'derived', m)
for c in ('R_seq', 'A_seq', 'hapX_seq', 'hapY_seq'):
    col(c, f'normal_frame {c}', {'R_seq': 'GRCh38 between the anchors', 'A_seq': 'R_seq with the truth applied',
                                 'hapX_seq': 'hapX between the anchors (assembly, else dipcall)', 'hapY_seq': 'hapY, same'}[c])

pv = []
for t in sorted(ABS, key=int):
    v, g, l, q, n = VS[t], g38[t], loc[t], gp[t], nf[t]
    o = {c: v[c] for c in ('truth_id', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind2', 'VAF_Ill', 'VAF_PB', 'RGN', 'n_platforms')}
    o['vaf_bin'] = vbin(float(v['VAF_Ill'])); o['perfect_platforms'] = ','.join(p for p in PLAT if perf(t, p))
    for p in PLAT:
        o.update({f'{p}_{c}': v[f'{p}_{c}'] for c in ('perfect', 'status', 'class', 'reason', 'alt_reads', 'spanning')})
    o.update({c: l[c] for c in ('first_node', 'last_node', 'n_ref_nodes', 'win_start0', 'win_end0', 'ev_lo0', 'ev_hi0')})
    o.update({c: q[c] for c in ('match', 'germline_used', 'germline_gt', 'closest_edit_distance', 'grch38_edit_distance', 'n_alt_paths',
                                'cap_hit', 'cap_info', 'n_completed', 'n_event_sets', 'primary_n_nodes', 'primary_event_ids',
                                'primary_d_grch38', 'primary_local_alt', 'primary_local_ref', 'truth_net', 'event_net', 'primary_path')})
    els = json.loads(q['primary_elements'])
    o['graph_primary_elements'] = ';'.join(compact(e) for e in els); o['graph_n_primary_elements'] = len(els)
    o['graph_in_event_subtypes'] = '+'.join(sorted({e['subtype'] for e in els if e['in_event']}))
    o['graph_event_element_sets'] = ';'.join('+'.join(s['elements']) + f"({s['n_paths']})" for s in json.loads(q['event_element_sets'] or '[]'))
    o.update({c: g[s] for c, s in (('grch38_ins_del', 'ins_del'), ('grch38_indel_seq', 'indel_seq'), ('grch38_unit', 'indel_unit'),
                                   ('grch38_len_change', 'indel_len'), ('grch38_period', 'ctx_period'), ('grch38_tract', 'ctx_tract'),
                                   ('grch38_class', 'ctx_class'), ('grch38_b5rule_class', 'b5rule_class'),
                                   ('grch38_b5rule_period', 'b5rule_period'), ('grch38_b5rule_tract', 'b5rule_tract'),
                                   ('cpg', 'cpg'), ('cpg_ti', 'cpg_ti'), ('titv', 'titv'), ('alt_hp_len', 'alt_hp_len'), ('tri', 'tri'))})
    o['snv_spectrum'] = spec(v['vcf_ref'], v['vcf_alt']) if v['kind2'] == 'SNV' else ''
    o.update({c: n[c] for c in NFC}); o['nf_flags'] = n['flags']
    o.update(pf_group=group(n['pf_label']), label_b5=label_b5(n), units_bin=units_bin(n), hp_tract_bin=hp_bin(n),
             germline_status=germ(n), mechanism_class=cmech(n) if v['kind2'] == 'INDEL' else '',
             mechanism_class_b5=cmech5(n) if v['kind2'] == 'INDEL' else '')
    o.update({c: n[c] for c in ('R_seq', 'A_seq', 'hapX_seq', 'hapY_seq')})
    pv.append(o)
names = [c for c, _, _ in COLS]
assert len(names) == len(set(names)) and all(set(o) == set(names) for o in pv), set(pv[0]) ^ set(names)
with open(A / 'per_variant.tsv', 'w') as w:
    w.write('\t'.join(names) + '\n')
    for o in pv:
        assert not any('\t' in str(o[c]) or '\n' in str(o[c]) for c in names)
        w.write('\t'.join(str(o[c]) for c in names) + '\n')
HCOLS = rd(D / 'hprc_columns.tsv')                         # c8: columns of the HPRC files of this folder
with open(A / 'columns.tsv', 'w') as w:
    w.write('file\tcolumn\tsource\tmeaning\n')
    for c, s, m in COLS:
        w.write(f'per_variant.tsv\t{c}\t{s}\t{m}\n')
    for r in HCOLS:
        w.write(f"{r['file']}\t{r['column']}\t{r['source']}\t{r['meaning']}\n")
PV = byid(pv)
print(f'per_variant.tsv {len(pv)} rows x {len(names)} columns; columns.tsv {len(COLS) + len(HCOLS)} rows ({len(HCOLS)} of the HPRC files, c8)')

# ---------------------------------------------------------------- numbers + tables
N = {}                                    # every number the README quotes
out = ['# COLO829T absorbed truths: repeat context (c4)', '',
       'Tables of `c4_tables.py` (definitions in its docstring and in README.md). Set = the truths that are perfect bypass on that '
       'platform; union = on >= 1 platform. Patient frame = COLO829BL (matched normal) verkko 2.1 hapX / hapY; GRCh38 frame = the '
       'truth VCF record. chr1-22 only (per_truth_COLO829T.tsv holds the chr1-22 alleles).', '']
pct = lambda a, b: f'{a:,} ({100 * a / b:.1f}%)' if b else f'{a:,}'
rate = lambda a, b: f'{a:,} / {b:,} ({100 * a / b:.1f}%)' if b else f'{a:,} / 0'


TBL = {}                                  # table id ('1', '2a', ...) -> markdown lines without the title (README)


def md(title, header, rows, note=None):
    body = ['| ' + ' | '.join(header) + ' |', '|---|' + '---:|' * (len(header) - 1)] + ['| ' + ' | '.join(str(x) for x in r) + ' |' for r in rows]
    TBL[title.split('.')[0]] = body + ([''] + [note] if note else [])
    out.extend([f'### {title}', ''] + body + [''] + ([note, ''] if note else []))


def per_set(title, kind, key, order, store=None):
    """count table of key(row) per perfect set (fiberseq / ONT / Illumina / union) for the absorbed truths of a kind."""
    sets = [(p, [o for o in pv if o['kind2'] == kind and o[f'{p}_perfect'] == 'True']) for p in PLAT]
    sets.append(('union', [o for o in pv if o['kind2'] == kind]))
    cnt = {s: C(key(o) for o in xs) for s, xs in sets}
    keys = [k for k in order if any(cnt[s][k] for s, _ in sets)] + sorted({k for s, _ in sets for k in cnt[s]} - set(order))
    md(title, ['class'] + [f'{s} ({kind})' for s, _ in sets],
       [[k] + [pct(cnt[s][k], len(xs)) for s, xs in sets] for k in keys] + [['total'] + [f'{len(xs):,}' for _, xs in sets]])
    if store:
        N[store] = {s: C({**cnt[s], 'total': len(xs)}) for s, xs in sets}
    return cnt


# 1. Q1: how many, per platform
KIND = {'SNP': 'SNV', 'INS': 'INDEL', 'DEL': 'INDEL'}
allk = C(KIND[r['kind']] for r in pt)
q1 = []
for s in PLAT + ('union',):
    row = {}
    for k in ('SNV', 'INDEL'):
        sel = [o for o in pv if o['kind2'] == k and (s == 'union' or o[f'{s}_perfect'] == 'True')]
        row[k] = dict(nc=sum(KIND[r['kind']] == k and r[f'{s}_status'] == 'no_candidate' for r in pt) if s != 'union' else None,
                      perfect=len(sel), spells=sum(o['match'] in ('exact', 'with_germline') for o in sel),
                      chr1=sum(o['chrom'] == 'chr1' for o in sel), match=C(o['match'] for o in sel),
                      cls=C(o[f'{s}_class'] for o in sel) if s != 'union' else None)
    q1.append((s, row))
N['q1'] = {s: r for s, r in q1}; N['all_truth'] = dict(allk)
N['all3'] = C(o['kind2'] for o in pv if o['n_platforms'] == '3')
md('1. Perfect bypass per platform (chr1-22; the chr1 column is a subset)',
   ['platform', 'SNV no candidate', 'SNV perfect bypass', 'of which a d9 path spells the truth', 'INDEL no candidate',
    'INDEL perfect bypass', 'of which a d9 path spells the truth', 'chr1 SNV / INDEL perfect bypass'],
   [[s, f"{r['SNV']['nc']:,}" if r['SNV']['nc'] is not None else '', pct(r['SNV']['perfect'], r['SNV']['nc']) if r['SNV']['nc'] else f"{r['SNV']['perfect']:,}",
     r['SNV']['spells'], f"{r['INDEL']['nc']:,}" if r['INDEL']['nc'] is not None else '',
     pct(r['INDEL']['perfect'], r['INDEL']['nc']) if r['INDEL']['nc'] else f"{r['INDEL']['perfect']:,}", r['INDEL']['spells'],
     f"{r['SNV']['chr1']} / {r['INDEL']['chr1']}"] for s, r in q1],
   f"All chr1-22 truth alleles: SNV {allk['SNV']:,}, INDEL {allk['INDEL']:,}. Perfect bypass on all three platforms: SNV "
   f"{N['all3']['SNV']}, INDEL {N['all3']['INDEL']}. Share of all truth: " +
   '; '.join(f"{s} SNV {100 * r['SNV']['perfect'] / allk['SNV']:.2f}% / INDEL {100 * r['INDEL']['perfect'] / allk['INDEL']:.1f}%" for s, r in q1))
MATCH = ['exact', 'with_germline', 'closest']
md('1b. Miss classes and d9 match of the perfect-bypass truths', ['platform'] + [f'SNV {c}' for c in ('S3', 'S4b', 'S4a', 'S2')] +
   ['INDEL I1'] + [f'SNV {m}' for m in MATCH] + [f'INDEL {m}' for m in MATCH],
   [[s] + ([sum(v for k, v in r['SNV']['cls'].items() if k.startswith(c + '_')) for c in ('S3', 'S4b', 'S4a', 'S2')] +
           [sum(v for k, v in r['INDEL']['cls'].items() if k.startswith('I1_'))] if r['SNV']['cls'] is not None else [''] * 5) +
    [r['SNV']['match'][m] for m in MATCH] + [r['INDEL']['match'][m] for m in MATCH] for s, r in q1],
   'S3 = the ALT is an existing graph SNP allele; S4b = repeat, absorbed by a graph path; S4a = repeat, the SNV becomes an edit on '
   'another branch; S2 = ALT not seen in same-length reads; I1 = the INDEL allele is fully in the graph. exact / with_germline / '
   'closest: c2 (a d9 path spells GRCh38 + truth; only with COLO829BL dipcall alleles applied too; no such path).')
for s, r in q1:
    if r['SNV']['cls'] is not None:
        assert sum(r['SNV']['cls'].values()) == r['SNV']['perfect'] and set(r['INDEL']['cls']) == {'I1_allele_fully_in_graph'}

# 2. INDELs, patient frame
cl = per_set('2a. INDELs, patient frame (COLO829BL): class label per perfect set', 'INDEL', lambda o: o['pf_label'], LABELS, 'pf_label')
per_set('2b. INDELs, patient frame: grouped', 'INDEL', lambda o: o['pf_group'], GROUPS, 'pf_group')
per_set('2c. INDELs, patient frame: units changed (repeat INDEL events only)', 'INDEL',
        lambda o: o['units_bin'] or '(not a repeat INDEL event)', UNB + ['(not a repeat INDEL event)'], 'units')
per_set("2d. INDELs, patient frame: the normal's homopolymer length (HP>=7 events only)", 'INDEL',
        lambda o: o['hp_tract_bin'] or '(not HP>=7)', HPB + ['(not HP>=7)'], 'hp_bin')
per_set('2e. INDELs, GRCh38 frame (b1 rule): class per perfect set', 'INDEL', lambda o: o['grch38_class'], G38, 'g38_indel')
per_set('2f. INDELs, patient frame with the b5 rule (HG008 rule): class label per perfect set', 'INDEL', lambda o: o['label_b5'],
        LABELS, 'label_b5')
ai = [o for o in pv if o['kind2'] == 'INDEL']
N['g38none_pf'] = C(o['pf_label'] for o in ai if o['grch38_class'] == 'none')
N['g38none_b5rule'] = C(o['grch38_b5rule_class'].split('(')[0] for o in ai if o['grch38_class'] == 'none')
x = C((o['grch38_class'], o['pf_group']) for o in ai)
md('2g. Union INDELs: GRCh38-frame class x patient-frame group', ['GRCh38 class \\ patient frame'] + GROUPS + ['total'],
   [[g] + [x[(g, k)] for k in GROUPS] + [sum(x[(g, k)] for k in GROUPS)] for g in G38 if any(x[(g, k)] for k in GROUPS)],
   'GRCh38-frame "none" INDELs (union) in the patient frame: ' + ', '.join(f'{k} {v}' for k, v in N['g38none_pf'].most_common()) +
   '. The same INDELs under the b5 rule in GRCh38 (c1 b5rule_class): ' + ', '.join(f'{k} {v}' for k, v in N['g38none_b5rule'].most_common()) + '.')
N['ncalt_other'] = C(o['other_hap_label'] for o in ai if o['pf_label'] == 'normal carries ALT')
N['ncalt_hap'] = C(('hapX' if o['hapX_allele'] == 'ALT' else '') + ('hapY' if o['hapY_allele'] == 'ALT' else '') for o in ai
                   if o['pf_label'] == 'normal carries ALT')
N['complex_other'] = C(o['other_hap_label'] for o in ai if o['pf_label'] == 'complex')
N['sign'] = C('same' if int(o['pf_len_change']) == int(o['grch38_len_change']) else 'opposite sign'
              if int(o['pf_len_change']) * int(o['grch38_len_change']) < 0 else 'other size' for o in ai if o['pf_kind'] == 'INDEL')
N['eventhap_allele'] = C(re.sub(r'[+-]\d+$', '', o[f"{o['event_hap']}_allele"]) for o in ai if o['event_hap'])
def ev_bases(o):                          # inserted / deleted bases of a patient-frame INDEL event ('off:pad+X>pad')
    x, y = o['pf_event'].split(':', 1)[1].split('>')
    return (x if len(x) > len(y) else y)[1:]


def motif(X):                             # canonical STR unit: smallest rotation over both strands
    u = X[:int(minp(X))]
    return min(min(r[i:] + r[:i] for i in range(len(r))) for r in (u, u[::-1].translate(RC)))


minp = lambda X: next(p for p in range(1, len(X) + 1) if len(X) % p == 0 and X[:p] * (len(X) // p) == X)
N['hp_base'] = C('A/T' if ev_bases(o)[0] in 'AT' else 'C/G' for o in ai if o['pf_label'] == 'HP>=7')
N['str_motif'] = C(motif(ev_bases(o)) for o in ai if o['pf_label'] == 'STR2-6>=3copies')
nca = [o for o in ai if o['pf_label'] == 'normal carries ALT']
N['nca'] = dict(rgn=C(o['RGN'] for o in nca), vaf=C(o['vaf_bin'] for o in nca), flagged=sum(bool(o['nf_flags']) for o in nca),
                asm_dip=sum(any(o[f'{h}_asm_call'] == o[f'{h}_dip_call'] == 'ALT' for h in ('hapX', 'hapY')) for o in nca),
                seg=C(f"{o['chrom']}:{int(o['vcf_pos']) // 10 ** 6}Mb" for o in nca))
N['complex_rgn'] = C(o['RGN'] for o in ai if o['pf_label'] == 'complex')
N['nca_other_ref'] = sum('REF' in (o['hapX_allele'], o['hapY_allele']) for o in nca)
N['hp_insdel'] = C('INS' if int(o['pf_len_change']) > 0 else 'DEL' for o in ai if o['pf_label'] == 'HP>=7')
N['str_unit'] = C(int(o['pf_unit']) for o in ai if o['pf_label'] == 'STR2-6>=3copies')
N['altlen'] = C(o['pf_label'] for o in ai if o['alt_len_haps'] and o['pf_label'] != 'normal carries ALT')
N['ivntr'] = C(o['pf_label'] for o in ai if o['in_ivntr'] == 'True')
N['dcopy'] = C(o['pf_label'] for o in ai if 'dipcall' in (o['hapX_copy'], o['hapY_copy']))
N['changed'] = sorted(((o['truth_id'], o['kind2'], PRE[o['truth_id']], o['pf_label']) for o in pv if PRE[o['truth_id']] != o['pf_label']), key=lambda z: int(z[0]))
N['single'] = C(o['pf_label_one_hap'] for o in ai if o['pf_label'] == 'single hap')
N['wide'] = C(o['pf_label'] for o in ai if o['anchor_wide'] == 'True')
N['flanklen'] = sum(any(abs(int(o[f'{h}_flank_len'] or 0)) >= 4 for h in ('hapX', 'hapY')) for o in ai)
N['flanklen_ev'] = sum(abs(int(o[f"{o['event_hap']}_flank_len"] or 0)) >= 4 for o in ai if o['event_hap'])
md('2h. Union INDELs: miscellany', ['item', 'counts'],
   [['normal carries ALT: other hap\'s event', ', '.join(f'{k} {v}' for k, v in N['ncalt_other'].most_common())],
    ['normal carries ALT: hap carrying it', ', '.join(f'{k} {v}' for k, v in N['ncalt_hap'].most_common())],
    ['complex: other hap\'s event', ', '.join(f'{k} {v}' for k, v in N['complex_other'].most_common())],
    ['patient-frame INDEL events vs the GRCh38 length change', ', '.join(f'{k} {v}' for k, v in N['sign'].most_common())],
    ['array allele of the event hap', ', '.join(f'{k} {v}' for k, v in N['eventhap_allele'].most_common())],
    ['HP>=7 events: ins / del', ', '.join(f'{k} {v}' for k, v in N['hp_insdel'].most_common())],
    ['HP>=7 events: base', ', '.join(f'{k} {v}' for k, v in N['hp_base'].most_common())],
    ['STR events: unit length', ', '.join(f'{k} {v}' for k, v in sorted(N['str_unit'].items()))],
    ['STR events: motif (canonical, both strands)', ', '.join(f'{k} {v}' for k, v in N['str_motif'].most_common())],
    ['normal carries ALT: SMaHT region', ', '.join(f'{k} {v}' for k, v in N['nca']['rgn'].most_common())],
    ['normal carries ALT: VAF_Ill', ', '.join(f'{k} {v}' for k, v in sorted(N['nca']['vaf'].items()))],
    ['normal carries ALT: assembly and dipcall both ALT on a hap / any c3 flag',
     f"{N['nca']['asm_dip']} / {N['nca']['flagged']} of {len(nca)}"],
    ['normal carries ALT: loci per chromosome Mb (>= 2)', ', '.join(f'{k} {v}' for k, v in N['nca']['seg'].most_common() if v >= 2)],
    ['complex: SMaHT region', ', '.join(f'{k} {v}' for k, v in N['complex_rgn'].most_common())],
    ['a normal hap has the ALT length but not the ALT (alt_len_haps), by label', ', '.join(f'{k} {v}' for k, v in N['altlen'].most_common())],
    ['in an imperfect VNTR / minisatellite (in_ivntr), by label', ', '.join(f'{k} {v}' for k, v in N['ivntr'].most_common())],
    ["a hap taken from dipcall's own copy (paralog fix), by label", ', '.join(f'{k} {v}' for k, v in N['dcopy'].most_common())],
    ['wide anchors (2 kb retry), by label', ', '.join(f'{k} {v}' for k, v in N['wide'].most_common())],
    ['single hap: the one hap\'s label', ', '.join(f'{k} {v}' for k, v in N['single'].most_common()) or '(none absorbed)'],
    ['label changed by the audit fixes, all absorbed truths (truth kind: before -> after)', '; '.join(f'{t} {k}: {a} -> {b}' for t, k, a, b in N['changed'])],
    ['a hap with >= 4 bp of reverted flank length (flank_len) / the event hap', f"{N['flanklen']} / {N['flanklen_ev']}"]])

# 3. background absorption rates
nfI = [nf[t] for t in INDELS]
bgl = C(n['pf_label'] for n in nfI); abl = C(n['pf_label'] for n in nfI if n['truth_id'] in ABS)
bgg = C(group(n['pf_label']) for n in nfI); abg = C(group(n['pf_label']) for n in nfI if n['truth_id'] in ABS)
N['bg_pf_label'] = (bgl, abl); N['bg_pf_group'] = (bgg, abg)
md('3a. All 1,912 chr1-22 truth INDELs, patient frame: share absorbed (union)', ['class', 'truth INDELs', 'absorbed', '% absorbed'],
   [[k, f'{bgl[k]:,}', abl[k], f'{100 * abl[k] / bgl[k]:.1f}%'] for k in LABELS if bgl[k]] +
   [['total', f'{len(nfI):,}', sum(abl.values()), f'{100 * sum(abl.values()) / len(nfI):.1f}%']])
bgm = C(cmech(n) for n in nfI); abm = C(cmech(n) for n in nfI if n['truth_id'] in ABS)
N['bg_mech'] = (bgm, abm)
md('3b. All chr1-22 truth INDELs, patient-frame mechanism class: share absorbed (union)', ['class', 'truth INDELs', 'absorbed', '%'],
   [[k, f'{bgm[k]:,}', abm[k], f'{100 * abm[k] / bgm[k]:.1f}%'] for k in MECH if bgm[k]])
for kd in ('INDEL', 'SNV'):
    allx = [g for g in g38.values() if g['kind2'] == kd]
    b = C(g['ctx_class'] for g in allx); a = C(g['ctx_class'] for g in allx if g['absorbed'] == 'True')
    N[f'bg_g38_{kd}'] = (b, a)
    md(f'3{"c" if kd == "INDEL" else "d"}. All chr1-22 truth {kd}s, GRCh38 frame (b1 rule): share absorbed (union)',
       ['class', f'truth {kd}s', 'absorbed', '%'],
       [[k, f'{b[k]:,}', a[k], f'{100 * a[k] / b[k]:.2f}%'] for k in G38 if b[k]] +
       [['total', f'{len(allx):,}', sum(a.values()), f'{100 * sum(a.values()) / len(allx):.2f}%']])

# 3e. SNVs
asnv = [o for o in pv if o['kind2'] == 'SNV']
per_set('3e. Absorbed SNVs, GRCh38-frame class per perfect set', 'SNV', lambda o: o['grch38_class'], G38, 'g38_snv')
per_set('3f. Absorbed SNVs, patient frame (COLO829BL) per perfect set', 'SNV', lambda o: o['pf_label'], LABELS, 'pf_snv')
per_set('3g. Absorbed SNVs, CpG transition (GRCh38) per perfect set', 'SNV',
        lambda o: 'CpG Ti' if o['cpg_ti'] == 'True' else 'other', ['CpG Ti', 'other'], 'cpg_snv')
allS = [g for g in g38.values() if g['kind2'] == 'SNV']
cb = C(g['cpg_ti'] for g in allS); ca = C(g['cpg_ti'] for g in allS if g['absorbed'] == 'True')
spb = C(spec(g['vcf_ref'], g['vcf_alt']) for g in allS); spa = C(spec(g['vcf_ref'], g['vcf_alt']) for g in allS if g['absorbed'] == 'True')
N['snv_cpg'] = (cb, ca); N['snv_spec'] = (spb, spa)
xc = C((o['grch38_class'], 'CpG Ti' if o['cpg_ti'] == 'True' else 'other', o['pf_label']) for o in asnv)
md('3h. All chr1-22 truth SNVs: share absorbed by CpG transition and substitution', ['class', 'truth SNVs', 'absorbed', '%'],
   [['CpG Ti', f"{cb['True']:,}", ca['True'], f"{100 * ca['True'] / cb['True']:.2f}%"],
    ['not CpG Ti', f"{cb['False']:,}", ca['False'], f"{100 * ca['False'] / cb['False']:.2f}%"]] +
   [[k, f'{spb[k]:,}', spa[k], f'{100 * spa[k] / spb[k]:.2f}%'] for k in ('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G')])
md('3i. Absorbed SNVs (union): GRCh38 class x CpG Ti x patient frame', ['GRCh38 class', 'CpG Ti', 'patient frame', 'SNVs'],
   [[a, b, c, v] for (a, b, c), v in sorted(xc.items(), key=lambda kv: -kv[1])])
N['snv_cls_union'] = C(re.sub(r'_.*', '', o[f'{p}_class']) for o in asnv for p in PLAT if o[f'{p}_perfect'] == 'True')
N['snv_xc'] = xc
N['snv_match'] = C(o['match'] for o in asnv)

# 4. COLO-specific strata (union INDELs; rates over all chr1-22 truth INDELs where the class is the patient-frame group)
def strata(title, key_all, key_abs, order, store):
    b = C((key_all(t), group(nf[t]['pf_label'])) for t in INDELS)
    a = C((key_all(t), group(nf[t]['pf_label'])) for t in INDELS if t in ABS)
    bt = C(key_all(t) for t in INDELS); at = C(key_all(t) for t in INDELS if t in ABS)
    rows = [[g] + [rate(a[(s, g)], b[(s, g)]) for s in order] for g in GROUPS if any(b[(s, g)] for s in order)]
    rows.append(['all INDELs'] + [rate(at[s], bt[s]) for s in order])
    md(title, ['patient-frame group \\ ' + store] + order, rows,
       'Cells: absorbed (union) / all chr1-22 truth INDELs of the class and stratum (% absorbed).')
    N[f'strata_{store}'] = (b, a, bt, at)


strata('4a. SMaHT region x patient-frame group, INDELs', lambda t: g38[t]['RGN'], None, RGN, 'RGN')
strata('4b. Truth VAF (Illumina) x patient-frame group, INDELs', lambda t: vbin(float(g38[t]['VAF_Ill'])), None, VAFB, 'VAF')
strata('4c. Normal-frame germline status x patient-frame group, INDELs', lambda t: germ(nf[t]), None, GERM, 'germline status')
hpr = C((g38[t]['RGN'], cmech(nf[t])) for t in INDELS); hpra = C((g38[t]['RGN'], cmech(nf[t])) for t in INDELS if t in ABS)
HPM = [f'HP>=7 1-unit, tract {b}' for b in HPB]
md('4d. HP>=7 1-unit events by tract length x SMaHT region, INDELs', ['tract \\ RGN'] + RGN,
   [[m.split(', ')[1]] + [rate(hpra[(r, m)], hpr[(r, m)]) for r in RGN] for m in HPM if any(hpr[(r, m)] for r in RGN)],
   'Cells: absorbed (union) / all chr1-22 truth INDELs whose patient-frame event is a 1-base change of a homopolymer of that length.')
N['hp_rgn'] = (hpr, hpra)
for kd in ('SNV',):
    allx = [g for g in g38.values() if g['kind2'] == kd]
    for nm, f, order in (('RGN', lambda g: g['RGN'], RGN), ('VAF', lambda g: vbin(float(g['VAF_Ill'])), VAFB)):
        b = C(f(g) for g in allx); a = C(f(g) for g in allx if g['absorbed'] == 'True')
        N[f'snv_{nm}'] = (b, a)
        out.append(f'SNVs by {nm}: ' + '; '.join(f'{s} {rate(a[s], b[s])}' for s in order))
out.append('')
x = C((o['germline_status'], o['units_bin'] or o['pf_group']) for o in ai)
N['germ_abs'] = C(o['germline_status'] for o in ai)
N['germ_hp'] = C((o['germline_status'], o['pf_group']) for o in ai)
N['pf_vs_g38_hp'] = C((o['pf_group'], o['hp_tract_bin'], o['germline_status']) for o in ai if o['pf_group'] == 'HP>=7')

# 5. HG008T side-by-side
hpv = rd(HA / 'per_variant.tsv'); hirc = rd(HA / 'indel_repeat_context.tsv'); hpt = rd(TR / 'per_truth_HG008T.tsv')
hb1 = byid(rd(HB / 'b1_truth_context.tsv')); hb5 = byid(rd(HB / 'b5_normal_frame.tsv'))
HABS = {r['truth_id'] for r in hpv}; HABSI = {r['truth_id'] for r in hpv if r['kind'] == 'INDEL'}


def hpclass(t):                                            # s14 pclass, unchanged
    n = hb5.get(t, {})
    if not n.get('nf_kind'):
        return 'no event'
    c = n['nf_class']
    if n['nf_kind'] != 'INDEL':
        return f"{n['nf_kind']} in {c}" if c not in ('none', 'complex', '') else n['nf_kind']
    return 'in_STR' if c.startswith('in_STR') else c


def hmech(t):
    n = hb5.get(t, {}); lab = hpclass(t)
    return mech(lab, abs(int(n['nf_len'])) // max(1, int(n['nf_unit'])) if n.get('nf_kind') == 'INDEL' else 0, int(n.get('nf_tract') or 0))


AUTO = {f'chr{i}' for i in range(1, 23)}
hI = [t for t, g in hb1.items() if g['kind2'] == 'INDEL' and g['chrom'] in AUTO]
hS = [t for t, g in hb1.items() if g['kind2'] == 'SNV' and g['chrom'] in AUTO]
assert len(hirc) == len(HABSI) == 2277 and len(hI) == 8496 and len(hS) == 8690
assert all(hpclass(r['truth_id']) == r['patient_class'] for r in hirc)
hcl = C(r['patient_class'] for r in hirc); hun = C(r['units_changed'] for r in hirc); hhp = C(r['hp_tract_bin'] for r in hirc)
assert (hcl['HP>=7'], hcl['STR2-6>=3copies'], hun['1 unit'], hhp['20-29']) == (1853, 292, 2103, 806)     # indel_repeat_context.md
hpres = C(r['category'] for r in hpv if r['kind'] == 'INDEL')
hplat = {p: (sum(r[f'{p}_perfect'] == 'True' and r['patient_class'] == 'HP>=7' for r in hirc), sum(r[f'{p}_perfect'] == 'True' for r in hirc))
         for p in ('PacBio', 'ONT', 'Illumina')}
hq1 = {}
for p in ('PacBio', 'ONT', 'Illumina'):
    for k in ('SNV', 'INDEL'):
        hq1[(p, k)] = (sum(KIND[r['kind']] == k and r[f'{p}_status'] == 'no_candidate' for r in hpt),
                       sum(r['kind'] == k and r[f'{p}_perfect'] == 'True' for r in hpv))
assert hq1[('PacBio', 'SNV')] == (398, 170) and hq1[('Illumina', 'INDEL')] == (3954, 1916)
cn = len(ai); hn = len(hirc)
cg = C(o['pf_group'] for o in ai); cgb = C(group(o['label_b5']) for o in ai); hg = C(group(r['patient_class']) for r in hirc)
md('5a. Perfect bypass, COLO829T vs HG008T (chr1-22)', ['', 'COLO829T', 'HG008T'],
   [['truth alleles SNV / INDEL', f"{allk['SNV']:,} / {allk['INDEL']:,}", f'{len(hS):,} / {len(hI):,}'],
    ['fiberseq / PacBio HiFi: SNV / INDEL perfect (% of no candidate)',
     ' / '.join(pct(N['q1']['fiberseq'][k]['perfect'], N['q1']['fiberseq'][k]['nc']) for k in ('SNV', 'INDEL')),
     ' / '.join(pct(hq1[('PacBio', k)][1], hq1[('PacBio', k)][0]) for k in ('SNV', 'INDEL'))],
    ['ONT: SNV / INDEL perfect', ' / '.join(pct(N['q1']['ONT'][k]['perfect'], N['q1']['ONT'][k]['nc']) for k in ('SNV', 'INDEL')),
     ' / '.join(pct(hq1[('ONT', k)][1], hq1[('ONT', k)][0]) for k in ('SNV', 'INDEL'))],
    ['Illumina: SNV / INDEL perfect', ' / '.join(pct(N['q1']['Illumina'][k]['perfect'], N['q1']['Illumina'][k]['nc']) for k in ('SNV', 'INDEL')),
     ' / '.join(pct(hq1[('Illumina', k)][1], hq1[('Illumina', k)][0]) for k in ('SNV', 'INDEL'))],
    ['union SNV (% of all truth SNVs)', pct(len(asnv), allk['SNV']), pct(len(HABS - HABSI), len(hS))],
    ['union INDEL (% of all truth INDELs)', pct(cn, allk['INDEL']), pct(hn, len(hI))]])
md('5b. Union absorbed INDELs, patient-frame group: COLO829T vs HG008T',
   ['group', 'COLO829T (c3 rule)', 'COLO829T (b5 rule)', 'HG008T (b5 rule)'],
   [[k, pct(cg[k], cn), pct(cgb[k], cn), pct(hg[k], hn)] for k in GROUPS] + [['total', cn, cn, f'{hn:,}']],
   'HG008T patient frame = the GIAB truth INFO event (HG008Nv62SOMATICVARIANT) in the HG008-N v6.2 event haplotype, always a '
   "simple event, so HG008T has no 'complex' or 'normal carries ALT' rows; 'unresolved' = no INFO event there. HG008T truths whose "
   f"tumor allele is the normal's other-haplotype allele (category HG008N_present_*): {hpres['HG008N_present_broad'] + hpres['HG008N_present_rare']} "
   f"INDELs ({100 * (hpres['HG008N_present_broad'] + hpres['HG008N_present_rare']) / hn:.1f}%), classed by their INFO event inside the rows above.")
cl5 = C(o['pf_label'] for o in ai); clb = C(o['label_b5'] for o in ai)
md('5c. Union absorbed INDELs, patient-frame class label: COLO829T vs HG008T',
   ['class', 'COLO829T (c3 rule)', 'COLO829T (b5 rule)', 'HG008T (b5 rule)'],
   [[k, pct(cl5[k], cn), pct(clb[k], cn), pct(hcl[k], hn)] for k in LABELS if cl5[k] or clb[k] or hcl[k]])
cu = C(o['units_bin'] for o in ai); ch = C(o['hp_tract_bin'] for o in ai)
cu5 = C(ubin(abs(int(o['pf_len_change'])) // int(o['pf_unit'])) for o in ai if o['label_b5'] in REP)
ch5 = C(hbin(int(o['pf_tract_b5'])) for o in ai if o['label_b5'] == 'HP>=7')
nrep_c = sum(cu[u] for u in UNB); nrep_h = sum(hun[u] for u in UNB); nhp_c = sum(ch[b] for b in HPB); nhp_h = sum(hhp[b] for b in HPB)
nrep_c5 = sum(cu5.values()); nhp_c5 = sum(ch5.values())
md('5d. Union absorbed INDELs: units changed and homopolymer length, COLO829T vs HG008T',
   ['', 'COLO829T (c3 rule): of all absorbed INDELs', 'COLO829T (c3): of repeat / HP>=7 events', 'COLO829T (b5 rule): of all',
    'COLO829T (b5): of repeat / HP>=7 events', 'HG008T (b5): of all absorbed INDELs', 'HG008T (b5): of repeat / HP>=7 events'],
   [[u, pct(cu[u], cn), pct(cu[u], nrep_c), pct(cu5[u], cn), pct(cu5[u], nrep_c5), pct(hun[u], hn), pct(hun[u], nrep_h)] for u in UNB] +
   [[f'HP {b} bp', pct(ch[b], cn), pct(ch[b], nhp_c), pct(ch5[b], cn), pct(ch5[b], nhp_c5), pct(hhp[b], hn), pct(hhp[b], nhp_h)] for b in HPB],
   'b5 rule = the HG008 rule on both sides (the same-rule comparison); the c3 rule is the COLO829T default of sections 2-5.')
N['side'] = dict(hplat=hplat, cn=cn, hn=hn, cg=cg, cgb=cgb, hg=hg, cu=cu, hun=hun, ch=ch, hhp=hhp, nrep_c=nrep_c, nrep_h=nrep_h, nhp_c=nhp_c,
                 nhp_h=nhp_h, cu5=cu5, ch5=ch5, nrep_c5=nrep_c5, nhp_c5=nhp_c5, hpres=hpres, hcl=hcl, cl5=cl5, clb=clb, hS=len(hS), hI=len(hI), hq1=hq1, habs_snv=len(HABS - HABSI))
hbg = C(hb1[t]['ctx_class'] for t in hI); hab = C(hb1[t]['ctx_class'] for t in hI if t in HABSI)
cbg, cab = N['bg_g38_INDEL']
md('5e. Background, GRCh38 frame (b1 rule): share of all chr1-22 truth INDELs absorbed', ['class', 'COLO829T', 'HG008T'],
   [[k, rate(cab[k], cbg[k]), rate(hab[k], hbg[k])] for k in G38 if cbg[k] or hbg[k]] +
   [['total', rate(sum(cab.values()), sum(cbg.values())), rate(sum(hab.values()), len(hI))]])
hbs = C(hb1[t]['ctx_class'] for t in hS); has_ = C(hb1[t]['ctx_class'] for t in hS if t in HABS)
csb, csa = N['bg_g38_SNV']
hcb = C(hb1[t]['cpg_ti'] for t in hS); hca = C(hb1[t]['cpg_ti'] for t in hS if t in HABS)
md('5f. Background, GRCh38 frame: share of all chr1-22 truth SNVs absorbed', ['class', 'COLO829T', 'HG008T'],
   [[k, rate(csa[k], csb[k]), rate(has_[k], hbs[k])] for k in G38 if csb[k] or hbs[k]] +
   [['CpG Ti', rate(N['snv_cpg'][1]['True'], N['snv_cpg'][0]['True']), rate(hca['True'], hcb['True'])],
    ['total', rate(sum(csa.values()), sum(csb.values())), rate(sum(has_.values()), len(hS))]])
hpl = C(hpclass(t) for t in hI); hpa = C(hpclass(t) for t in hI if t in HABSI)
bgl5 = C(label_b5(nf[t]) for t in INDELS); abl5 = C(label_b5(nf[t]) for t in INDELS if t in ABS)
md('5g. Background, patient frame: share of all chr1-22 truth INDELs absorbed',
   ['class', 'COLO829T (c3 rule)', 'COLO829T (b5 rule)', 'HG008T (b5 rule)'],
   [[k, rate(abl[k], bgl[k]), rate(abl5[k], bgl5[k]), rate(hpa[k], hpl[k])] for k in LABELS if bgl[k] or bgl5[k] or hpl[k]],
   'Compare the two b5 columns (the same rule). HG008T has no imperfect-VNTR check, no single-hap label and no '
   "'complex' / 'normal carries ALT' events (its INFO events are simple); 'no event' = no INFO event.")
hm = C(hmech(t) for t in hI); hma = C(hmech(t) for t in hI if t in HABSI)
bgm5 = C(cmech5(nf[t]) for t in INDELS); abm5 = C(cmech5(nf[t]) for t in INDELS if t in ABS)
md('5h. Background, patient-frame mechanism class (INDEL records): share absorbed',
   ['class', 'COLO829T (c3 rule)', 'COLO829T (b5 rule)', 'HG008T (b5 rule)'],
   [[k, rate(abm[k], bgm[k]), rate(abm5[k], bgm5[k]), rate(hma[k], hm[k])] for k in MECH if bgm[k] or bgm5[k] or hm[k]],
   'Compare the two b5 columns (the same rule, tract = the b5 tract). HG008T row counts differ from mechanism.md, which counts SNV '
   'and INDEL records together; here INDEL records only.')
N['side'].update(hbg=hbg, hab=hab, hbs=hbs, has_=has_, hcb=hcb, hca=hca, hpl=hpl, hpa=hpa, hm=hm, hma=hma, bgl5=bgl5, abl5=abl5,
                 bgm5=bgm5, abm5=abm5)

open(A / 'repeat_context.md', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))

# ---------------------------------------------------------------- README.md (every number from the counters above)
q = N['q1']; S = N['side']
r_ = lambda a, b: f'{a:,} / {b:,} ({100 * a / b:.1f}%)' if b else f'{a:,} / 0'
p_ = lambda a, b: f'{100 * a / b:.1f}%'
pfu = N['pf_group']['union']; plu = N['pf_label']['union']; un = N['units']['union']; hb = N['hp_bin']['union']
g38u = N['g38_indel']['union']; bgm, abm = N['bg_mech']; bgl, abl = N['bg_pf_label']
cbg, cab = N['bg_g38_INDEL']; csb, csa = N['bg_g38_SNV']; cpb, cpa = N['snv_cpg']; spb, spa = N['snv_spec']
gb, ga, gbt, gat = N['strata_germline status']; rb, ra, rbt, rat = N['strata_RGN']; vb, va, vbt, vat = N['strata_VAF']
hpr, hpra = N['hp_rgn']; snr = N['snv_RGN']
nI, nS = len(ai), len(asnv)
unit_ch = sum(plu[k] for k in REP)                      # repeat-unit changes (HP>=7, STR, HP4-6, VNTR>6 events)
near_rep = plu['in_STR'] + plu['imperfect VNTR']         # in / next to a repeat, not a whole-unit change
altlen_n = sum(N['altlen'].values())
nrep = sum(un[u] for u in UNB)
g38none = N['g38none_pf']; g38n = sum(g38none.values())
g38n_unit = sum(v for k, v in g38none.items() if k in REP); g38n_near = g38none['in_STR'] + g38none['imperfect VNTR']
g38n_snvrep = sum(v for k, v in g38none.items() if k.startswith('SNV in'))
s3 = {p: sum(v for k, v in q[p]['SNV']['cls'].items() if k.startswith('S3_')) for p in PLAT}
xc = N['snv_xc']
snv_nonrep = sum(v for (a, b, c), v in xc.items() if a == 'none' and c == 'SNV')
snv_nonrep_cpg = sum(v for (a, b, c), v in xc.items() if a == 'none' and c == 'SNV' and b == 'CpG Ti')
snv_pf = N['pf_snv']['union']
snv_rep_pf = sum(v for k, v in snv_pf.items() if k.startswith('SNV in'))
hmb = lambda b: (S['hma'][f'HP>=7 1-unit, tract {b}'], S['hm'][f'HP>=7 1-unit, tract {b}'])
cmb = lambda b: (abm[f'HP>=7 1-unit, tract {b}'], bgm[f'HP>=7 1-unit, tract {b}'])
cmb5 = lambda b: (S['abm5'][f'HP>=7 1-unit, tract {b}'], S['bgm5'][f'HP>=7 1-unit, tract {b}'])     # b5 rule (HG008 comparison)
h_present = S['hpres']['HG008N_present_broad'] + S['hpres']['HG008N_present_rare']
h_hp_truth = S['hpl']['HP>=7']; c_hp_truth = bgl['HP>=7']; c_hp_truth5 = S['bgl5']['HP>=7']
h_snv_rep = S['has_']['HP>=7'] + S['has_']['STR2-6>=3copies']
diff_b5 = sum(o['pf_label'] != o['label_b5'] for o in ai)
ties = C(o['event_tie'] for o in pv)
odd = [o for o in ai if o['pf_group'] in ('other repeat', 'SNV in repeat', 'complex')]
odd_germ = sum(o['germline_status'] not in ('both REF', 'one hap NA, other REF', 'both NA') for o in odd)
gref_snv, gref_cx = (sum(gb[(g, k)] for g in ('both REF', 'one hap NA, other REF', 'both NA')) for k in ('SNV in repeat', 'complex'))
for t, lab in (('949', 'complex'), ('2719', 'in_STR'), ('1996', 'SNV in STR2-6>=3copies')):     # README examples
    assert PV[t]['pf_label'] == lab, (t, PV[t]['pf_label'])
assert 'hapX' in PV['949']['alt_len_haps'] and 'hapX' in PV['1996']['alt_len_haps']
cx_match = C(o['match'] for o in ai if o['pf_label'] == 'complex')
in2rep = sum(o['pf_label'] == 'in_STR' and o['label_b5'] in ('HP>=7', 'STR2-6>=3copies') for o in ai)
NAWHY = {'anchor': 'no common anchor', 'no_hits': 'no assembly hit', 'no_anchors': 'no GRCh38 anchors'}
dhp = [100 * cmb5(b)[0] / cmb5(b)[1] - 100 * hmb(b)[0] / hmb(b)[1] for b in HPB[:4]]      # COLO829T - HG008T, both b5
dstr = 100 * S['abm5']['STR 1-unit'] / S['bgm5']['STR 1-unit'] - 100 * S['hma']['STR 1-unit'] / S['hm']['STR 1-unit']
na_loci = [(t, nf[t]['absorbed'] == 'True', nf[t]['hapX_place_na'] or nf[t]['hapX_asm_na'], nf[t]['hapY_place_na'] or nf[t]['hapY_asm_na'])
           for t in sorted(nf, key=int) if nf[t]['pf_label'] == 'NA']
single_loci = [(t, nf[t]['absorbed'] == 'True', nf[t]['pf_label_one_hap']) for t in sorted(nf, key=int) if nf[t]['pf_label'] == 'single hap']
dcl = [n for n in nf.values() if 'dipcall' in (n['hapX_copy'], n['hapY_copy'])]
n_dcopy_loci, n_dcopy_abs = len(dcl), sum(n['absorbed'] == 'True' for n in dcl)
n_dcopy_changed = sum(PRE[n['truth_id']] != n['pf_label'] for n in dcl if n['absorbed'] == 'True')
n_phase_bad = sum(any(int(n[f'{h}_phase_mismatch'] or 0) for h in ('hapX', 'hapY')) for n in nf.values())
n_phase_bad0 = sum(any(int(n[f'{h}_phase_mismatch_ragtag'] or 0) for h in ('hapX', 'hapY')) for n in nf.values())
dupf = lambda o: any(f.split(':')[-1] in ('mapq<20', 'other_scaffold_hit') for f in o['nf_flags'].split(','))
nca_dup = sum(dupf(o) for o in nca); nca_dup_dc = sum(dupf(o) and 'dipcall' in (o['hapX_copy'], o['hapY_copy']) for o in nca)
nca_146 = [o for o in nca if o['chrom'] == 'chr1' and 146_000_000 <= int(o['vcf_pos']) < 150_000_000]; nca_146_dup = sum(dupf(o) for o in nca_146)
C3_RUN = '3 min 54 s, 0.26 GB'
C4_RUN = '11 s, 0.44 GB'
_pos = collections.defaultdict(list)
for r in pt:
    _pos[r['chrom']].append((int(r['vcf_pos']) - 1, r['truth_id']))
shared = [(t, [u for x, u in _pos[nf[t]['chrom']] if u != t and int(nf[t]['arr_lo0']) <= x < int(nf[t]['arr_hi0'])]) for t in sorted(ABS, key=int)]
shared = [(t, us) for t, us in shared if us]
t_ = lambda k: '\n'.join(TBL[k])

# ---------------------------------------------------------------- 7-8. HPRC haplotypes / populations and the cross-evaluation (c5-c8)
_s8 = importlib.util.spec_from_file_location('c8', A / 'c8_populations.py'); C8 = importlib.util.module_from_spec(_s8); _s8.loader.exec_module(C8)
SPP, HC = C8.SP, C8.CATS                                   # superpopulations, categories (c7)
hv = byid(rd(A / 'hprc_per_variant.tsv')); hlp = byid(rd(A / 'per_locus_populations.tsv')); HV = [hv[t] for t in sorted(hv, key=int)]
hsum = json.load(open(D / 'hprc_summary.json')); kck = rd(D / 'hprc_check_kmer.tsv'); hcov = byid(rd(D / 'hprc_window_coverage.tsv'))
hind = rd(A / 'hprc_individuals.tsv'); hpop = rd(A / 'hprc_populations.tsv'); hnode = rd(A / 'hprc_per_node.tsv')
assert set(hv) == set(hlp) == set(hcov) == ABS
for t in ABS:
    assert hv[t]['kind2'] == VS[t]['kind2'] and all(hv[t][f'{p}_perfect'] == VS[t][f'{p}_perfect'] for p in PLAT), t
    assert (hv[t]['COLO829BL_hapX_allele'], hv[t]['COLO829BL_hapY_allele'], hv[t]['pf_label']) == (nf[t]['hapX_allele'], nf[t]['hapY_allele'], nf[t]['pf_label']), t
    assert not hv[t]['category'].startswith('normal_present') or nf[t]['pf_label'] == 'normal carries ALT', t
    assert (hv[t]['germline_like_alt_len'] == 'True') == (bool(nf[t]['alt_len_haps']) and nf[t]['pf_label'] != 'normal carries ALT'), t
fl = lambda x: float(x) if x not in ('', None) else None
med = lambda xs: statistics.median(xs) if xs else float('nan')
hsel = lambda f, k=None, rows=None: [r for r in (HV if rows is None else rows) if f(r) and (k is None or r['kind2'] == k)]
iscat = lambda *cs: (lambda r: r['category'] in cs)
PB, PR, AB, AM = HC
PRES = (PB, PR)
HN = {}                                                    # numbers of sections 7-8


def tb(key, title, header, rows, note=None):
    """a README-only table (sections 7-8; repeat_context.md stays the repeat-context tables), printed."""
    body = ['| ' + ' | '.join(header) + ' |', '|---|' + '---:|' * (len(header) - 1)] + ['| ' + ' | '.join(str(x) for x in r) + ' |' for r in rows]
    TBL[key] = body + ([''] + [note] if note else [])
    print('\n'.join([f'### {key}. {title}', ''] + TBL[key] + ['']))


def ccell(rs):
    c = C(r['category'] for r in rs)
    return ' / '.join(f'{c[x]:,}' for x in HC) + f' ({len(rs):,})'


rows = []
for sn in PLAT + ('union',):
    f = (lambda r: True) if sn == 'union' else (lambda r, sn=sn: r[f'{sn}_perfect'] == 'True')
    row = [f'**{sn}**' if sn == 'union' else sn]
    for sc in ('chr1-22', 'chr1'):
        for k in ('SNV', 'INDEL'):
            rs = [r for r in HV if f(r) and r['kind2'] == k and (sc == 'chr1-22' or r['chrom'] == 'chr1')]
            assert C(r['category'] for r in rs) == C({c: n for c, n in hsum['categories'][f'{sn} {k} {sc}'].items() if c != 'total'}), (sn, k, sc)
            row.append(f'**{ccell(rs)}**' if sn == 'union' and sc == 'chr1-22' else ccell(rs))
    rows.append(row)
tb('h1', 'Categories per perfect set (chr1-22 and chr1)', ['set', 'SNV chr1-22', 'INDEL chr1-22', 'SNV chr1', 'INDEL chr1'], rows,
   'Cells: normal_present_broad / normal_present_rare / normal_absent_HPRC_other / other_ambiguous (total). chr1 = the chr1 truths '
   '(every COLO829T truth is inside the BED).')
RULES = [
    ('1', AM, 'no d9 path spells the truth (c2 `closest`; no read data to take the elements from)', lambda r: r['sub_reason'].startswith('closest')),
    ('1', AM, 'no somatic element (every path element alone spells a germline allele, or none)', lambda r: r['sub_reason'] in ('only_germline_elements', 'no_somatic_element')),
    ('1', AM, 'the `with_germline` path needs a germline allele that dipcall phases to the other haplotype than the c3 event hap',
     lambda r: r['sub_reason'] == 'with_germline:germline_allele_off_event_hap'),
    ('1', AM, 'no HPRC haplotype traverses the window completely', lambda r: r['sub_reason'] == 'HPRC_no_complete_haplotype'),
    ('1', AM, 'both COLO829BL haplotypes carry the ALT (homozygous germline; truth conflict)', lambda r: r['sub_reason'] == 'evidence_conflict:ALT_on_both_normal_haps'),
    ('2', PB, "a COLO829BL haplotype carries the ALT over the whole tandem array ('normal carries ALT'); exact-allele HPRC frequency >= 0.20", iscat(PB)),
    ('2', PR, 'the same, frequency < 0.20', iscat(PR)),
    ('3', AM, 'a COLO829BL haplotype has no array sequence and the other one does not carry the ALT', lambda r: r['sub_reason'].startswith('normal_hap_unresolved')),
    ('3', AM, 'no complete HPRC haplotype walks every somatic element (CHM13 may)', lambda r: r['sub_reason'].startswith('no_HPRC_carrier')),
    ('4', AB, 'both haplotypes resolved, neither carries the ALT; sub_class `exact_allele`: >= 1 complete HPRC haplotype spells the exact allele',
     lambda r: r['category'] == AB and r['sub_class'] == 'exact_allele'),
    ('4', AB, '`element_set_other_allele`: HPRC haplotypes walk all elements of an ALT path but spell another allele over the array',
     lambda r: r['category'] == AB and r['sub_class'] == 'element_set_other_allele'),
    ('4', AB, '`recombinant_pieces`: every element is walked by some HPRC haplotype, a whole ALT path by none',
     lambda r: r['category'] == AB and r['sub_class'] == 'recombinant_pieces')]
rows = [[o, f'`{c}`', txt, len(hsel(f, 'SNV')), len(hsel(f, 'INDEL'))] for o, c, txt, f in RULES]
assert sum(r[3] + r[4] for r in rows) == len(HV) and all(sum(f(r) for _, _, _, f in RULES) == 1 for r in HV)
tb('h2', 'Category rules, in order (union, chr1-22)', ['order', 'category', 'rule', 'SNV', 'INDEL'], rows)
fq = lambda col, f, k=None, pos=False: [fl(r[col]) for r in hsel(f, k) if fl(r[col]) is not None and (not pos or fl(r[col]) > 0)]
lv = [('**exact allele** (the categories use this)', 'hprc_exact_allele_freq', 'the walk spells GRCh38 + truth over the tandem array (with_germline: + the germline alleles)'),
      ('element set', 'hprc_element_set_freq', 'the walk contains all elements of any ALT path'),
      ('element (min)', 'hprc_element_min_freq', 'the walk contains each element (minimum over the elements)')]
rows = []
for nm, col, what in lv:
    a = [med(fq(col, iscat(AB), k)) for k in ('SNV', 'INDEL')]
    ap = [med(fq(col, iscat(AB), k, True)) for k in ('SNV', 'INDEL')]
    rows.append([nm, what, f'{a[0]:.3f} / {a[1]:.3f}' + (f' ({ap[0]:.3f} / {ap[1]:.3f} at loci with >= 1 carrier)' if col == 'hprc_exact_allele_freq' else ''),
                 f"{med(fq(col, iscat(*PRES))):.3f}"])
tb('h3', 'HPRC membership levels: median frequency over the complete traversals (union, chr1-22)',
   ['level', 'carrier =', 'normal_absent_HPRC_other SNV / INDEL', 'normal_present (SNV + INDEL)'], rows)
rows = []
for sc in ('exact_allele', 'element_set_other_allele', 'recombinant_pieces'):
    for k in ('SNV', 'INDEL'):
        x = hsel(lambda r: r['category'] == AB and r['sub_class'] == sc, k)
        if not x:
            continue
        pr_ = C(hlp[r['truth_id']]['private'] for r in x)
        rows.append([f'`{sc}`', k, len(x), f"{med(fq('hprc_exact_allele_freq', lambda r: r in x)):.3f}", f"{med(fq('hprc_element_set_freq', lambda r: r in x)):.3f}",
                     f"{med([int(hlp[r['truth_id']]['n_individuals']) for r in x]):g}", pr_['shared'], pr_['one superpopulation'],
                     pr_['one population'], pr_['one individual'], pr_['none'], sum(r['CHM13_carries_exact'] == 'True' for r in x)])
tb('h4', 'normal_absent_HPRC_other by HPRC support (union, chr1-22)',
   ['sub_class', 'kind', 'loci', 'median exact-allele freq', 'median element-set freq', 'median carrier individuals', '>= 2 superpops',
    'one superpop', 'one population', 'one individual', 'none', 'CHM13 exact'], rows,
   "Carrier individuals / superpopulations (per_locus_populations.tsv): exact-allele carriers; for `element_set_other_allele` / "
   "`recombinant_pieces` the haplotypes walking all elements of an ALT path. First match of one individual > one population > one superpopulation.")
oerows = lambda f, rows_=None: [r for r in (HV if rows_ is None else rows_) if f(r) and r['somatic_elements'] and r['hprc_n_complete'] != '0']
fmt_oe = lambda o, s: f'{o[s][0] / o[s][1]:.2f} (z {o[s][2]:+.1f})' if o[s][1] and o[s][2] is not None else ''
rows, HN['oe'] = [], {}
for c in HC + ['all']:
    for nm, lab in (('exact', 'exact allele'), ('element_set', 'element set')):
        x = oerows(lambda r: c == 'all' or r['category'] == c)
        o = C8.oe(x, nm); HN['oe'][(c, nm)] = o
        rows.append([f'`{c}`' if c != 'all' else 'all', lab, len(x), sum(v[0] for v in o.values())] + [fmt_oe(o, s) for s in SPP] +
                    [sum(r['CHM13_carries_exact' if nm == 'exact' else 'CHM13_carries_element_set'] == 'True' for r in x)])
tb('h5', 'Superpopulation of the HPRC carriers: observed / expected carrier haplotypes (union, chr1-22)',
   ['category', 'level', 'loci', 'carrier haplotypes'] + [f'{s} O/E (z)' for s in SPP] + ['CHM13 carries'], rows,
   "Expected per locus = carriers x the superpopulation's share of that locus' complete haplotypes; z hypergeometric, ignoring the "
   "pairing of haplotypes within individuals and linkage between loci (it overstates significance). HPRC v1.1: AFR 23, AMR 16, "
   "EAS 4, SAS 1 samples, no EUR; the COLO829 donor is a European (white) male.")
assert abs(HN['oe'][(AB, 'exact')]['AFR'][0] / HN['oe'][(AB, 'exact')]['AFR'][1] - 1.178) < 5e-4      # c7 hprc_tables.md (oe() reused)
rows = []
for p_r in sorted(hpop, key=lambda r: (r['superpopulation'] == 'reference', r['superpopulation'], r['population_code'])):
    pres_e = sum(int(p_r[f'{c}_INDEL_exact_carrier_haps']) for c in PRES); pres_n = sum(int(p_r[f'{c}_INDEL_complete_haps']) for c in PRES)
    rows.append([p_r['superpopulation'] if p_r['superpopulation'] != 'reference' else 'CHM13', p_r['population_code'] or 'CHM13',
                 p_r['population_name'][:40], p_r['haplotypes'], f"{float(p_r[f'{AB}_INDEL_carrier_fraction']):.3f}",
                 f"{float(p_r[f'{AB}_SNV_carrier_fraction']):.3f}", f'{pres_e / pres_n:.3f}' if pres_n else ''])
tb('h6', '1000 Genomes populations: exact-allele carrier haplotypes / complete traversals, summed over the loci (hprc_populations.tsv)',
   ['superpop', 'population', 'name', 'HPRC haplotypes', 'absent INDEL', 'absent SNV', 'present INDEL (broad + rare)'], rows)
hapr = [r for r in hind if r['sample'] != 'CHM13']; chm = next(r for r in hind if r['sample'] == 'CHM13')
ex_ai = lambda r: int(r[f'{AB}_INDEL_exact']); sh_ai = lambda r: int(r[f'{AB}_INDEL_exact']) / int(r[f'{AB}_INDEL_complete'])
HN['ind'] = dict(med=med([ex_ai(r) for r in hapr]), share=med([sh_ai(r) for r in hapr]), lo=min(hapr, key=ex_ai), hi=max(hapr, key=ex_ai),
                 snv=med([int(r[f'{AB}_SNV_exact']) for r in hapr]), n_ai=len(hsel(iscat(AB), 'INDEL')),
                 sp={s: (med([ex_ai(r) for r in hapr if r['superpopulation'] == s]), med([sh_ai(r) for r in hapr if r['superpopulation'] == s])) for s in SPP},
                 chm=(ex_ai(chm), sh_ai(chm)))
grp = sorted({PV[t]['pf_group'] for t in ABS}, key=GROUPS.index)
rows = [[g] + [f'{len(hsel(lambda r: r["category"] == c and PV[r["truth_id"]]["pf_group"] == g, k))}' for k in ('SNV', 'INDEL') for c in HC] for g in grp]
tb('h7', 'Patient-frame group (section 2) x category (union, chr1-22)', ['patient-frame group'] + [f'{k} {c}' for k in ('SNV', 'INDEL') for c in HC], rows)
gl = hsel(lambda r: r['germline_like_alt_len'] == 'True')
HN['gl'] = C((r['category'], r['sub_class'] or r['sub_reason']) for r in gl)
# HG008T (analysis/graph_absorbed_somatic_20261001: per_variant.tsv, per_locus_populations.tsv, hprc_individuals.tsv)
HM = {'HG008N_present_broad': PB, 'HG008N_present_rare': PR, 'HG008N_absent_HPRC_other': AB, 'other_ambiguous': AM}
hh = [dict(r, category=HM[r['category']], kind2=r['kind']) for r in hpv]
hhl = byid(rd(HA / 'per_locus_populations.tsv')); hhind = [r for r in rd(HA / 'hprc_individuals.tsv') if r['sample'] != 'CHM13']
hcat = C((r['category'], r['kind2']) for r in hh)
assert (hcat[(PB, 'SNV')] + hcat[(PR, 'SNV')], hcat[(PB, 'INDEL')] + hcat[(PR, 'INDEL')], hcat[(AB, 'SNV')], hcat[(AB, 'INDEL')]) == (5, 105, 150, 1935)   # HG008 README
ho = C8.oe(oerows(iscat(AB), hh), 'exact'); hob = C8.oe(oerows(iscat(PB), hh), 'exact')
assert round(ho['AFR'][0] / ho['AFR'][1], 3) == 1.107 and round(hob['AMR'][0] / hob['AMR'][1], 3) == 1.098       # HG008 tables.md


def side(rows_, lp, ind, nm):
    """the comparison numbers of one sample (rows with c7 / s9-style columns, its per_locus_populations, its haplotypes)."""
    n = C((r['category'], r['kind2']) for r in rows_); tot = C(r['kind2'] for r in rows_)
    pc = lambda cs, k: f'{sum(n[(c, k)] for c in cs):,} ({100 * sum(n[(c, k)] for c in cs) / tot[k]:.1f}%)'
    ex = [r for r in rows_ if r['category'] == AB and r['sub_class'] == 'exact_allele']
    pv_ = C(lp[r['truth_id']]['private'] for r in ex)
    fx = lambda f, k=None, pos=False: med([fl(r['hprc_exact_allele_freq']) for r in rows_ if f(r) and (k is None or r['kind2'] == k)
                                           and fl(r['hprc_exact_allele_freq']) is not None and (not pos or fl(r['hprc_exact_allele_freq']) > 0)])
    o, ob = C8.oe(oerows(iscat(AB), rows_), 'exact'), C8.oe(oerows(iscat(PB), rows_), 'exact')
    return [f"{tot['SNV']:,} / {tot['INDEL']:,}", f'{pc(PRES, "SNV")} / {pc(PRES, "INDEL")}', f'{pc((AB,), "SNV")} / {pc((AB,), "INDEL")}',
            f"{sum(r['kind2'] == 'SNV' for r in ex):,} / {sum(r['kind2'] == 'INDEL' for r in ex):,}", f'{pc((AM,), "SNV")} / {pc((AM,), "INDEL")}',
            f"{med([int(r['hprc_n_complete']) for r in rows_]):g}",
            f"{fx(iscat(AB), 'SNV'):.3f} / {fx(iscat(AB), 'INDEL'):.3f}", f"{fx(iscat(AB), 'SNV', True):.3f} / {fx(iscat(AB), 'INDEL', True):.3f}",
            f'{fx(iscat(*PRES)):.3f}', f"{pv_['shared']:,} ({100 * pv_['shared'] / len(ex):.1f}%) / {pv_['one individual']} ({100 * pv_['one individual'] / len(ex):.1f}%)",
            ' / '.join(f'{o[s][0] / o[s][1]:.2f} ({o[s][2]:+.1f})' for s in SPP), ' / '.join(f'{ob[s][0] / ob[s][1]:.2f} ({ob[s][2]:+.1f})' for s in SPP[:2]),
            f"{sum(r['CHM13_carries_exact'] == 'True' for r in rows_ if r['category'] == AB):,}",
            f"{med([int(r[ind[1]]) / int(r[ind[2]]) for r in ind[0]]):.3f} ({med([int(r[ind[1]]) for r in ind[0]]):g} of {n[(AB, 'INDEL')]:,})"]


SIDE_ROWS = ['absorbed truth alleles, SNV / INDEL', 'germline filtering (the normal carries the ALT: present broad + rare), SNV / INDEL',
             'pangenome-induced false negative (absent from the normal, HPRC carries), SNV / INDEL', '  of which >= 1 HPRC haplotype has the exact allele, SNV / INDEL',
             'ambiguous, SNV / INDEL', 'median complete HPRC traversals per locus (of 88)', 'absent: median exact-allele frequency, SNV / INDEL',
             'absent: the same at loci with >= 1 carrier, SNV / INDEL', 'present: median exact-allele frequency',
             'absent exact-allele loci carried in >= 2 superpopulations / by one HPRC individual only',
             'absent, exact allele: O/E (z) AFR / AMR / EAS / SAS', 'present broad, exact allele: O/E (z) AFR / AMR',
             'CHM13 has the exact allele (absent loci)', 'per HPRC haplotype: median share of the absent INDELs it traverses whose exact allele it carries (median count)']
cside = side(HV, hlp, (hapr, f'{AB}_INDEL_exact', f'{AB}_INDEL_complete'), 'COLO829T')
hside = side(hh, hhl, (hhind, 'HG008N_absent_HPRC_other_INDEL_exact', 'HG008N_absent_HPRC_other_INDEL_complete'), 'HG008T')
tb('h8', 'HG008T vs COLO829T: HPRC membership and cross-evaluation (union, chr1-22)', ['', 'COLO829T', 'HG008T'],
   [[a, b, c] for a, b, c in zip(SIDE_ROWS, cside, hside)],
   'HG008T from analysis/graph_absorbed_somatic_20261001 (s9 / s13; recomputed here with the same functions and checked against its '
   'README / tables.md). The same category rules: HG008T `HG008N_*` = COLO829T `normal_*`; closest loci are ambiguous on both sides; '
   'HG008T chose among several ALT paths with re-decoded reads, COLO829T takes the c2 primary. % = of the absorbed alleles of that kind.')
kok = [r for r in kck if r['status'] == 'ok']
HN['kmer'] = dict(loci=len(kck), no_anchor=len(kck) - len(kok), ok=len(kok), agree=sum(r['agree'] == 'True' for r in kok),
                  pairs=sum(int(r['n_complete']) for r in kok), unanch=sum(int(r['n_unanchored']) for r in kok),
                  c7=sum(int(r['c7_exact']) for r in kok), km=sum(int(r['kmer_exact']) for r in kok),
                  dis=[(r['truth_id'], r['c7_exact'], r['kmer_exact'], r['n_complete']) for r in kok if r['agree'] != 'True'])
fl9 = hsum['floor']
low9 = [n for n in hnode if n['min_node_coverage'] and int(n['min_node_coverage']) < 9]
HN['floor'] = dict(branch=fl9['branch_elements'], low=len(low9), low_chm=sum(n['CHM13'] == 'True' for n in low9),
                   low_other=sorted((n['truth_id'], int(n['min_node_coverage'])) for n in low9 if n['CHM13'] != 'True'))
HN['vcf'] = {(so, k): (sum(r[f'{so}_vcf_agree'] == 'True' for r in hsel(lambda r: r['n_somatic_elements'] != '0', k)),
                       len(hsel(lambda r: r['n_somatic_elements'] != '0', k))) for so in ('d9', 'full') for k in ('SNV', 'INDEL')}
nc_ = [int(hcov[t]['n_complete']) for t in ABS]
HN['cov'] = dict(all=med(nc_), SNV=med([int(hcov[t]['n_complete']) for t in ABS if VS[t]['kind2'] == 'SNV']),
                 INDEL=med([int(hcov[t]['n_complete']) for t in ABS if VS[t]['kind2'] == 'INDEL']), lt44=sum(n < 44 for n in nc_),
                 zero=[t for t in ABS if hcov[t]['n_complete'] == '0'])
hc = C((r['category'], r['kind2']) for r in HV); hct = C(r['category'] for r in HV); nH = len(HV)
sr = lambda pre: sum(r['sub_reason'].startswith(pre) for r in HV)
HN['am'] = dict(closest=sr('closest'), near=sr('closest:graph_allele_nearer_truth'), nohprc=sr('no_HPRC_carrier'), nohprc_chm=sr('no_HPRC_carrier:CHM13_carries'),
                off=sr('with_germline:germline_allele_off_event_hap'), unres=sr('normal_hap_unresolved'))
exl = hsel(lambda r: r['category'] == AB and r['sub_class'] == 'exact_allele')
HN['ex'] = dict(n=len(exl), snv=sum(r['kind2'] == 'SNV' for r in exl), indel=sum(r['kind2'] == 'INDEL' for r in exl),
                shared=sum(hlp[r['truth_id']]['private'] == 'shared' for r in exl), one_ind=sum(hlp[r['truth_id']]['private'] == 'one individual' for r in exl),
                ind={k: med([int(hlp[r['truth_id']]['n_individuals']) for r in exl if r['kind2'] == k]) for k in ('SNV', 'INDEL')},
                fr={k: med(fq('hprc_exact_allele_freq', lambda r: r in exl, k)) for k in ('SNV', 'INDEL')})
HN['nca_am'] = [r['truth_id'] for r in HV if r['pf_label'] == 'normal carries ALT' and r['category'] == AM]
HN['ab_grp'] = {k: C(PV[r['truth_id']]['pf_group'] for r in hsel(iscat(AB), k)) for k in ('SNV', 'INDEL')}
HN['pop_rng'] = {s: (min(float(r[f'{AB}_INDEL_carrier_fraction']) for r in hpop if r['superpopulation'] == s),
                     max(float(r[f'{AB}_INDEL_carrier_fraction']) for r in hpop if r['superpopulation'] == s)) for s in SPP}
HN['present_rgn'] = C(VS[r['truth_id']]['RGN'] for r in hsel(iscat(*PRES)))
HN['pr_zero'] = [r['truth_id'] for r in hsel(iscat(PR)) if r['hprc_exact_allele_carriers'] == '0']
oeab, oepb = HN['oe'][(AB, 'exact')], HN['oe'][(PB, 'exact')]
oef = lambda o, s: f'{o[s][0] / o[s][1]:.2f} (z {o[s][2]:+.1f})'
pp_ = lambda a, b: f'{a:,} ({100 * a / b:.1f}%)'
hh_ab = sum(1 for r in hh if r['category'] == AB); hh_pr = sum(1 for r in hh if r['category'] in PRES)
nca_n = len(nca)
R = []
R += [f"""# COLO829T somatic truths absorbed by the HPRC graph: repeat context, and whose alleles they are (COLO829BL or other HPRC individuals)

2026-10-05. COLO829T version of `analysis/graph_absorbed_somatic_20261001` (HG008T): the sequence context of the lost truths
(sections 1-6) and the HPRC haplotypes, individuals and populations that carry them (sections 7-8). Pansoma makes a
candidate (and a tensor) only from the edits a read has relative to the graph nodes it aligns to. When a somatic allele is
already a path of the HPRC v1.1 d9 graph, its reads follow that branch or skip edge with no edit, and the truth gets no
tensor. This analysis takes every COLO829T somatic truth allele lost this way and asks:

1. how many SNVs and INDELs, per platform (fiberseq, ONT, Illumina);
2. in what sequence context they sit: homopolymer, STR, other repeats or none, in the patient's own germline sequence
   (the matched normal COLO829BL, verkko 2.1 hapX / hapY) and in GRCh38;
3. how the absorption rate depends on that context, the SMaHT region, the truth VAF and the normal's own genotype;
4. how this compares with HG008T;
5. which HPRC v1.1 haplotypes, individuals and populations carry each absorbed allele and its graph nodes, and whether the
   allele is the patient's own germline allele (COLO829BL: germline filtering) or other people's (a pangenome-induced false
   negative).

## Short answer

- **Count.** The ALT reads align perfectly (no edit) to a graph branch or skip edge for:
  - fiberseq: {q['fiberseq']['SNV']['perfect']} SNVs / {q['fiberseq']['INDEL']['perfect']} INDELs; ONT: {q['ONT']['SNV']['perfect']} / {q['ONT']['INDEL']['perfect']}; Illumina: {q['Illumina']['SNV']['perfect']} / {q['Illumina']['INDEL']['perfect']}.
  - That is {min(100 * q[p]['SNV']['perfect'] / q[p]['SNV']['nc'] for p in PLAT):.1f}-{max(100 * q[p]['SNV']['perfect'] / q[p]['SNV']['nc'] for p in PLAT):.1f}% of each platform's no-candidate SNVs and {min(100 * q[p]['INDEL']['perfect'] / q[p]['INDEL']['nc'] for p in PLAT):.1f}-{max(100 * q[p]['INDEL']['perfect'] / q[p]['INDEL']['nc'] for p in PLAT):.1f}% of its no-candidate INDELs. Every
    such INDEL is class I1 of the miss analysis.
  - Union of the three platforms: **{nS} SNVs + {nI} INDELs = {nS + nI} truth alleles** ({p_(nI, allk['INDEL'])} of all {allk['INDEL']:,} chr1-22 truth
    INDELs, {100 * nS / allk['SNV']:.2f}% of the {allk['SNV']:,} SNVs); chr1: {q['union']['SNV']['chr1']} + {q['union']['INDEL']['chr1']}.
  - A d9 path spells the truth (alone or with the normal's nearby germline alleles) at {q['union']['SNV']['spells']} SNVs / {q['union']['INDEL']['spells']} INDELs.
- **The absorbed INDELs are repeat changes.** In the patient frame (union, {nI}):
  - **{unit_ch} ({p_(unit_ch, nI)}) change the length of a tandem repeat by whole units**: homopolymer `HP>=7` {pct(pfu['HP>=7'], nI)}, `STR` {pct(pfu['STR'], nI)},
    `HP4-6` {plu['HP4-6']}, `VNTR>6` {plu['VNTR>6']}. {un['1 unit']} of these {nrep} ({p_(un['1 unit'], nrep)}) change exactly one unit (one base of a homopolymer, one
    copy of an STR unit).
  - {near_rep} more sit in or next to a repeat without being a whole-unit change (`in_STR` {plu['in_STR']}, inside an imperfect VNTR /
    minisatellite {plu['imperfect VNTR']}), and {pfu['SNV in repeat']} are substitutions inside a repeat (the GRCh38 INDEL record absorbs a germline length
    difference). Only {pfu['non-repeat']} ({p_(pfu['non-repeat'], nI)}) are non-repeat.
  - The homopolymers are mostly 10-29 bp: {', '.join(f'{b} bp {hb[b]}' for b in HPB)}; {N['hp_base']['A/T']} of {pfu['HP>=7']} are A/T. STR motifs: {', '.join(f'({k})n {v}' for k, v in N['str_motif'].most_common(3))}.
  - The rest: the normal already carries the ALT ({pfu['normal carries ALT']}, {p_(pfu['normal carries ALT'], nI)}), a complex change ({pfu['complex']}, {p_(pfu['complex'], nI)}), unresolved {pfu['unresolved']}
    (no sequence, or only one haplotype).
  - **Germline-like or somatic?** At {altlen_n} more absorbed INDELs ({', '.join(f'`{k}` {v}' for k, v in N['altlen'].most_common())}) a COLO829BL
    haplotype already has exactly the ALT length but other bases than GRCh38 + truth (mostly a germline SNP in the
    array, or a tie where the other haplotype is GRCh38). Either the tumor
    changed the patient's repeat by one unit, or the truth is a germline length allele written against GRCh38; the tumor
    allele was not assembled, so this is not decided. So up to {pfu['normal carries ALT'] + altlen_n} ({p_(pfu['normal carries ALT'] + altlen_n, nI)}) of the absorbed INDELs may be
    germline-like, not only the {pfu['normal carries ALT']} ({p_(pfu['normal carries ALT'], nI)}) where the normal carries the ALT exactly (section 2).
- **The absorption rate rises with homopolymer length.** Over all {allk['INDEL']:,} truth INDELs, a one-base change of a homopolymer is
  absorbed at {', '.join(f'{100 * cmb(b)[0] / cmb(b)[1]:.1f}% ({b} bp)' for b in HPB[:4])}; a one-unit STR change at
  {p_(abm['STR 1-unit'], bgm['STR 1-unit'])}; a non-repeat INDEL at {p_(abm['non-repeat INDEL'], bgm['non-repeat INDEL'])} ({abm['non-repeat INDEL']} / {bgm['non-repeat INDEL']}; INDELs inside an imperfect VNTR are not
  counted as non-repeat: {abm['imperfect VNTR']} / {bgm['imperfect VNTR']} absorbed).
- **The patient's germline matters.**
  - At {pfu['normal carries ALT']} INDELs a COLO829BL haplotype already has the truth allele over the whole tandem array; {N['nca']['rgn']['Extreme']} of them are
    in SMaHT Extreme regions. There the graph acts as germline filtering, as for HG008T's {h_present} `HG008N_present` INDELs.
    {pct(abl['normal carries ALT'], bgl['normal carries ALT'])} of the {bgl['normal carries ALT']} truth INDELs whose allele the normal carries are absorbed.
  - A germline repeat-length allele in the normal (another array length on >= 1 haplotype) raises the rate for STRs
    ({p_(ga[('germline other length', 'STR')], gb[('germline other length', 'STR')])} vs {p_(ga[('both REF', 'STR')], gb[('both REF', 'STR')])} when both haplotypes equal GRCh38), hardly for homopolymers ({p_(ga[('germline other length', 'HP>=7')], gb[('germline other length', 'HP>=7')])} vs {p_(ga[('both REF', 'HP>=7')], gb[('both REF', 'HP>=7')])}).
- **SNVs are graph SNP alleles, not repeat events.** S3 (the ALT is an existing graph SNP allele) is {', '.join(f'{s3[p]} of {q[p]["SNV"]["perfect"]} ({p})' for p in PLAT)}.
  {snv_nonrep} of the {nS} are non-repeat in both frames; {sum(v for (a, b, c), v in xc.items() if b == 'CpG Ti')} are CpG transitions, which are absorbed at
  {100 * cpa['True'] / cpb['True']:.2f}% against {100 * cpa['False'] / cpb['False']:.2f}% for the other truth SNVs.
- **Region, not VAF.** {r_(rat['Easy'], rbt['Easy'])} of the INDELs in SMaHT Easy regions are absorbed, against {r_(rat['Difficult'], rbt['Difficult'])} (Difficult) and
  {r_(rat['Extreme'], rbt['Extreme'])} (Extreme). Easy regions hold the short homopolymers, but the gap remains at the same length
  (10-14 bp: {p_(hpra[('Easy', 'HP>=7 1-unit, tract 10-14')], hpr[('Easy', 'HP>=7 1-unit, tract 10-14')])} Easy vs {p_(hpra[('Difficult', 'HP>=7 1-unit, tract 10-14')], hpr[('Difficult', 'HP>=7 1-unit, tract 10-14')])} Difficult). The truth VAF hardly matters: {', '.join(f'{p_(vat[b], vbt[b])} ({b})' for b in VAFB if vbt[b])}.
- **Whose allele is it? Mostly other people's** (sections 7-8; HPRC v1.1 walks of the d9 graph, complete traversals only).
  - **Pangenome-induced false negative: {pp_(hct[AB], nH)}** ({hc[(AB, 'SNV')]} SNVs, {hc[(AB, 'INDEL')]} INDELs). Neither COLO829BL haplotype carries the
    ALT, but HPRC haplotypes walk the same nodes. At {HN['ex']['n']} of them ({HN['ex']['snv']} SNVs / {HN['ex']['indel']} INDELs) at least one HPRC haplotype
    carries exactly the allele: median exact-allele frequency {HN['ex']['fr']['SNV']:.3f} / {HN['ex']['fr']['INDEL']:.3f} (SNV / INDEL), carried by a median of
    {HN['ex']['ind']['SNV']:g} / {HN['ex']['ind']['INDEL']:g} HPRC individuals, in >= 2 superpopulations at {HN['ex']['shared']} loci and by a single individual at {HN['ex']['one_ind']}.
    The INDELs are the repeat slippages of section 2 ({HN['ab_grp']['INDEL']['HP>=7']} `HP>=7`, {HN['ab_grp']['INDEL']['STR']} `STR`); the SNVs are common graph SNP alleles.
  - **Germline filtering: {pp_(hct[PB] + hct[PR], nH)}**, the 'normal carries ALT' loci of section 2 ({hct[PB]} common in HPRC, {hct[PR]} rare).
  - **Ambiguous: {pp_(hct[AM], nH)}**, mainly the {HN['am']['closest']} 'closest' loci (no d9 path spells the truth).
  - **Populations.** The carriers of the absent alleles lean AFR, the panel's largest and most diverse group: observed / expected
    carrier haplotypes AFR {oef(oeab, 'AFR')}, AMR {oef(oeab, 'AMR')}, EAS {oef(oeab, 'EAS')}, SAS {oef(oeab, 'SAS')}. HPRC v1.1 has no
    EUR haplotype; the COLO829 donor is a European (white) male.
  - **Same as HG008T**: {p_(hc[(AB, 'INDEL')], nI)} vs {p_(sum(1 for r in hh if r['category'] == AB and r['kind2'] == 'INDEL'), S['hn'])} of the absorbed INDELs are pangenome-induced false negatives,
    {p_(hc[(PB, 'INDEL')] + hc[(PR, 'INDEL')], nI)} vs {p_(sum(1 for r in hh if r['category'] in PRES and r['kind2'] == 'INDEL'), S['hn'])} germline filtering; each HPRC haplotype carries the exact allele of a median
    {100 * HN['ind']['share']:.1f}% of the absent INDELs it traverses (HG008T 20.0%).
- **Compared with HG008T: the same mechanism, a different truth set** (both sides classed with the HG008 b5 rule here).
  - Per-length rates are close up to 19 bp: one-base homopolymer change {', '.join(f'{100 * cmb5(b)[0] / cmb5(b)[1]:.1f} vs {100 * hmb(b)[0] / hmb(b)[1]:.1f}% ({b} bp)' for b in HPB[:4])}
    (COLO829T vs HG008T); COLO829T is lower at 20-29 bp. One-unit STR change {p_(S['abm5']['STR 1-unit'], S['bgm5']['STR 1-unit'])} vs {p_(S['hma']['STR 1-unit'], S['hm']['STR 1-unit'])}.
  - Homopolymers are a smaller share of the absorbed INDELs ({p_(S['cgb']['HP>=7'], nI)} vs {p_(S['hg']['HP>=7'], S['hn'])}) mainly because the COLO829T truth has fewer
    and shorter homopolymer INDELs: {p_(c_hp_truth5, allk['INDEL'])} of its truth INDELs are `HP>=7` events (HG008T {p_(h_hp_truth, S['hI'])}), and only {S['bgm5']['HP>=7 1-unit, tract >=30']} is a
    one-base change of a homopolymer >= 30 bp (HG008T {S['hm']['HP>=7 1-unit, tract >=30']:,}). Overall INDEL absorption: {p_(nI, allk['INDEL'])} vs {p_(S['hn'], S['hI'])}.
  - Absorbed SNVs: {100 * nS / allk['SNV']:.2f}% of truth SNVs vs {100 * S['habs_snv'] / S['hS']:.1f}%. HG008T's were mostly homopolymer / STR events written as GRCh38 SNVs
    ({h_snv_rep} of {S['habs_snv']} in a GRCh38 `HP>=7` or `STR`); COLO829T's are mostly non-repeat graph SNP alleles ({csa['HP>=7'] + csa['STR2-6>=3copies']} of {nS} in such a repeat).
"""]
R += [f"""## Definitions

| Term | Meaning |
|---|---|
| d9 graph | HPRC v1.1 Minigraph-Cactus GRCh38 graph, frequency-filtered (nodes on fewer than ~9 of the 90 haplotypes removed). The COLO829T tensors were built on it |
| perfect bypass | status `no_candidate` on the platform (tensor sets of 2026-09-30) and the majority of the ALT-like reads leave GRCh38 through a non-GRCh38 node run or an edge that skips GRCh38 bases, with no edit within +-5 bp: read-level reason `site_bypassed:branch:no_edit` / `site_bypassed:skip_edge:no_edit` of the COLO829T miss analysis. **Absorbed** = perfect bypass on >= 1 platform (the union) |
| perfect set | the truths that are perfect bypass on one platform; per-set tables count truth alleles |
| d9 match | (c2) a d9 path between the window anchors spells GRCh38 + truth: `exact`; or only with nearby COLO829BL dipcall PASS alleles applied too: `with_germline`; no such path, the nearest graph haplotypes: `closest` |
| element | one departure of the ALT path from GRCh38: `B:x:y:>n1>n2` = a run of non-GRCh38 nodes between GRCh38 nodes x and y; `S:x:y` = an edge x->y that skips GRCh38 bases. Subtypes `ins_branch`, `del_skip`, `snv_branch`, `mnv_branch`, `replacing_branch` |
| GRCh38 frame | the truth VCF record in GRCh38, classed with the HG008 audit b1 rule (c1): unit = minimal period of the inserted / deleted bases, tract = the longest exact period-unit stretch touching the event |
| tandem array | (c3, the HG008 s8b rule) the event window grown by every exact periodic stretch of period 1-60 touching it (cap 3 kb), bracketed by GRCh38 24-mers that are unique and found once in each placed COLO829BL haplotype |
| hap copy | (c3) where a COLO829BL haplotype's sequence is cut from: the best hit of a 10-kb GRCh38 window on the same-chromosome ragtag scaffold (minimap2 asm5); where that copy is suspect (MAPQ < 20, a covering hit on another scaffold, a dipcall phased het SNV contradicted, or no hit) the contig copy dipcall itself aligned there (`dipcall copy`, the paralog fix). Anchors are searched within 400 bp of the array, else within 2 kb (`wide anchor`) so that both haplotypes share them |
| array allele | a COLO829BL haplotype's sequence between the anchors (assembly; the dipcall reconstruction where the assembly cannot be anchored), with germline differences outside the array reverted: `REF` (= GRCh38), `ALT` (= GRCh38 + truth), `germline_len` (another length), `germline_seq` (same length, other bases), `NA`. `flank_len` = the reverted length (>= 4 bp flagged; the hap's total change is `total_len_change`) |
| event hap / patient frame | the haplotype whose array allele is closest to GRCh38 + truth (Levenshtein; tie -> the one equal to GRCh38, else hapX). The **patient-frame event** is the edit from that allele to GRCh38 + truth, left-normalised. COLO829T has no truth INFO event (HG008T used the GIAB `HG008Nv62SOMATICVARIANT`), so the event is derived from the normal |
| `HP>=7` / `HP4-6` | INDEL event with unit 1 in a homopolymer of >= 7 / 4-6 bp (the normal's length, the stretch holding the event's own bases) |
| `STR2-6>=3copies` (`STR`) | unit 2-6 with >= 3 copies |
| `VNTR>6` | unit > 6, the inserted / deleted sequence repeated next to the event |
| `in_STR` | none of these, but the event sits inside or touches a period 1-6 stretch >= 10 bp without being a whole-unit change of it: e.g. 2719, a 2-bp TC insertion inside a (GT)n; most absorbed ones are 1-bp changes of another base next to a homopolymer (the b5 rule calls these `HP>=7`) |
| `imperfect VNTR` | a `none` event inside an imperfect tandem repeat (VNTR / minisatellite with unit 7-150 bp and indels between copies) that reaches >= 20 bp beyond the exact-period array: 12-mers recurring 7-150 bp apart cover >= 60 % of every 31-bp window of a >= 40-bp run holding the truth (GRCh38 +-600 bp; `SNV in imperfect VNTR` for substitutions) |
| `none` | INDEL event not in a repeat |
| `SNV in <class>` / `SNV` | the patient-frame event is a substitution: inside an `HP>=7` / `HP4-6` / `STR` stretch, or not. On a GRCh38 INDEL record this means the record absorbs a germline length difference |
| `complex` | the event changes REF and ALT bases at once (no single INDEL or SNV turns the closest haplotype into the ALT) |
| `normal carries ALT` | a COLO829BL haplotype already has GRCh38 + truth over the whole array |
| `NA` / `single hap` (unresolved) | no haplotype sequence (no common anchor, no assembly hit) / only one haplotype has a sequence and it is not the ALT: the missing one may carry the ALT, so no class is given (the one-hap class is in `pf_label_one_hap`) |
| ALT-length hap | (`alt_len_haps`) an INDEL truth where a COLO829BL haplotype's array allele has exactly the ALT length but is not the ALT (other bases, mostly a germline SNP): the normal may already carry the ALT length |
| group | `HP>=7` / `STR` / other repeat (`HP4-6`, `VNTR>6`, `in_STR`, `imperfect VNTR`) / SNV in repeat / non-repeat (`none`, `SNV`) / normal carries ALT / complex / unresolved (`NA`, `single hap`). Repeat-unit changes = `HP>=7`, `STR`, `HP4-6`, `VNTR>6` events; in or next to a repeat = `in_STR`, `imperfect VNTR` |
| b5 rule | the HG008 audit b5 rule verbatim (`pf_class_b5`, `label_b5`, `pf_tract_b5`): the tract is any period-unit stretch merely touching the event; no imperfect-VNTR check. The default (c3) rule requires the stretch to hold the event's own bases. HG008T numbers use the b5 rule, so every COLO829T vs HG008T rate comparison uses the b5 rule on both sides |
| units changed | abs(length change) / unit, for `HP>=7`, `HP4-6`, `STR`, `VNTR>6` events |
| mechanism class | (background rates) `HP>=7 1-unit, tract <bin>` = one-base change of a homopolymer of that length; `STR 1-unit`; `repeat multi-unit` = `HP>=7` / `STR` with >= 2 units; the other labels as above. COLO829T: c3 rule (`mechanism_class`) and b5 rule (`mechanism_class_b5`, label and tract of the b5 rule); HG008T: the same function on its b5 table |
| germline status | from the two array alleles, first that applies: normal carries ALT; germline other length (`germline_len` on >= 1 hap); germline same length, other bases; both REF; one hap NA, other REF; both NA |
| RGN | SMaHT region of the truth (truth VCF): Easy / Difficult / Extreme |
| VAF bin | the truth's Illumina VAF (`VAF_Ill`): < 0.1, 0.1-0.25, 0.25-0.4, >= 0.4 |
| miss classes | S3 = the ALT is an existing graph SNP allele; S4b = repeat, absorbed by a graph path; S4a = repeat, the SNV becomes an edit on another branch; S2 = ALT not seen in same-length reads; I1 = the INDEL allele is fully in the graph |
| somatic elements | (c7) `exact`: every element of the chosen ALT path (the c2 primary; no read data to choose); `with_germline`: minus the elements that alone spell the path's germline alleles; `closest`: the c2 closest primary's elements minus germline-alone ones (the locus is ambiguous) |
| complete traversal | an HPRC haplotype walk (c5, from the d9 GFA) that visits both window anchors in order. Only these count in HPRC frequencies: 88 HPRC haplotypes (44 samples x 2); CHM13 apart; GRCh38 never |
| exact allele | (c7) a haplotype's sequence over the tandem array, between the nearest GRCh38 nodes it shares with the ALT path (or over the whole window), equals GRCh38 + truth; for `with_germline` also the ALT path's allele (GRCh38 + truth + the patient's germline alleles) |
| element set | a walk contains every element of an ALT path (node level; any candidate path) |
| category | (c7, the HG008 s9 rules, section 8) `normal_present_broad` / `_rare` (germline filtering), `normal_absent_HPRC_other` (pangenome-induced false negative; sub_class `exact_allele` / `element_set_other_allele` / `recombinant_pieces`), `other_ambiguous` |
| RARE_AF | 0.20 on the exact-allele frequency among complete traversals (broad >= 0.20 > rare) |
| O/E | observed / expected carrier haplotypes of a superpopulation: expected = carriers x its share of the locus' complete haplotypes, summed over the loci |

## Data

| Item | Source |
|---|---|
| graph | `/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.*` (d9 sqlite / reference-path index via `indel_graph_paths_20260930/classify.py`) |
| truth | `/scratch/jshen/data/pansoma_v2_tensors/COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz`; chr1-22 alleles in `analysis/tensor_recall_20260930/per_truth_COLO829T.tsv` (all in the BED, all PASS; VAF_Ill, VAF_PB, RGN) |
| no-candidate statuses | `per_truth_COLO829T.tsv` (tensor sets of 2026-09-30, now in `pansoma_v2_tensors/backup_ch6_linear100_20261001`) |
| read-level reasons, miss classes | `analysis/tensor_recall_20260930/colo829t_miss/{{fiberseq,ONT,Illumina}}/{{snv,indel}}_no_candidate.tsv` |
| COLO829BL | SMaHT verkko 2.1 donor-specific assembly, ragtag GRCh38 scaffolds hapX / hapY (`/scratch/qfu/COLO829BL_DSA/ragtag_hg38/`, read only); dipcall vs GRCh38 `dipcall_hg38/dipcall_hg38.dip.vcf.gz` + `.dip.bed`, GT = hapY\\|hapX (c3 check on {sum(PHASE.values()):,} distinct isolated phased het SNVs near the loci: {PHASE['consistent']:,} consistent, {PHASE['swapped']} swapped, {PHASE['other']} other, {PHASE['undetermined']} undetermined, {PHASE['conflicting']} with different verdicts at two loci); dipcall's own contig alignments `dipcall_hg38.hap{{1,2}}.paf.gz` (hap1 = hapY) and the raw verkko contigs for the dipcall copy |
| GRCh38 | `GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta` |
| HG008T (comparison) | `analysis/graph_absorbed_somatic_20261001/{{per_variant,indel_repeat_context,per_locus_populations,hprc_individuals}}.tsv`, `per_truth_HG008T.tsv`, audit tables `b1_truth_context.tsv` / `b5_normal_frame.tsv` |
| HPRC walks | the d9 GFA `hprc-v1.1-mc-grch38.d9.gfa` (P and W lines; c5) |
| HPRC genotypes | the d9 VCF and the full-graph `vg deconstruct` VCF (`hprc-v1.1-mc-grch38.raw.vcf.gz`, read through the HG008 `full_vcf_index`), 44 sample columns + CHM13 (c6) |
| populations | HPRC v1.1 sample -> 1000 Genomes population / superpopulation: a byte-identical copy of the HG008 analysis' `hprc_sample_metadata.tsv` (HPRC Year 1 metadata, checked with IGSR / Coriell). AFR 23 / AMR 16 / EAS 4 / SAS 1 / **EUR 0** samples. The COLO829 donor is a 45-year-old white male (ATCC) |

- The 2026-10-01 tensor rebuild (channel 6, build AF 0.08) cannot change this set: a read with no edit gives no candidate
  in any build.
"""]
R += ['## 1. How many truths are lost because the reads align perfectly to a graph branch', '', t_('1'), '',
      f"- **Platforms.** {N['all3']['INDEL']} INDELs and {N['all3']['SNV']} SNVs are perfect bypass on all three platforms.",
      f"- **Miss classes and d9 match** (`c2_graph_paths.py`):", '', t_('1b'), '',
      f"- At the {q['union']['INDEL']['match']['closest']} INDEL and {q['union']['SNV']['match']['closest']} SNV 'closest' loci no d9 path spells the truth, so the reads that align perfectly there "
      "follow another allele (most likely a germline repeat length; the reads were not re-decoded for COLO829T).", '',
      '## 2. Repeat context of the absorbed INDELs', '',
      'Each INDEL truth is classed in the **patient frame**: the change from the closest COLO829BL haplotype to GRCh38 + truth '
      'over the tandem array (`c3_normal_frame.py`). This is the frame that matters: the GRCh38 record of a repeat INDEL '
      "depends on the GRCh38 repeat length, while the tumor changed the patient's own repeat.", '',
      '### Patient frame, grouped (per perfect set)', '', t_('2b'), '',
      f"- **Per platform.** Illumina's perfect set is {p_(N['pf_group']['Illumina']['HP>=7'], N['pf_group']['Illumina']['total'])} homopolymer; fiberseq and ONT are "
      f"{p_(N['pf_group']['fiberseq']['HP>=7'], N['pf_group']['fiberseq']['total'])} / {p_(N['pf_group']['ONT']['HP>=7'], N['pf_group']['ONT']['total'])} homopolymer and have relatively more STR, "
      f"normal-carries-ALT and complex loci (the same direction as HG008T: {', '.join(f'{k} {p_(a, b)}' for k, (a, b) in S['hplat'].items())} homopolymer).",
      f"- **Detailed labels** (`repeat_context.md` 2a): `in_STR` {plu['in_STR']}, `HP4-6` {plu['HP4-6']}, `VNTR>6` {plu['VNTR>6']}; substitution events "
      f"{', '.join(f'`{k}` {plu[k]}' for k in ('SNV in STR2-6>=3copies', 'SNV in HP4-6', 'SNV in HP>=7') if plu[k])}; non-repeat `none` {plu['none']}, `SNV` {plu['SNV']}.",
      f"- **Size of the change.** {un['1 unit']} of the {nrep} repeat INDEL events ({p_(un['1 unit'], nrep)}) change one unit, {un['2-3 units']} change 2-3 units, {un['>3 units']} more.",
      f"- **Homopolymer length** (the normal's): {', '.join(f'{b} bp {hb[b]}' for b in HPB)}. Bases A/T {N['hp_base']['A/T']}, C/G {N['hp_base']['C/G']}; "
      f"insertions {N['hp_insdel']['INS']}, deletions {N['hp_insdel']['DEL']}.",
      f"- **STR motifs** (canonical unit, both strands): {', '.join(f'{k} {v}' for k, v in N['str_motif'].most_common())}.",
      f"- **Patient frame vs GRCh38 record.** Of the {sum(N['sign'].values())} patient-frame INDEL events, {N['sign']['same']} have the GRCh38 record's length change, {N['sign']['other size']} "
      f"another size and {N['sign']['opposite sign']} the opposite sign (e.g. GRCh38 says a deletion, but the patient's haplotypes are shorter, so the "
      f"tumor inserted). The event haplotype's array allele is GRCh38 at {N['eventhap_allele']['REF']} loci, another germline length at "
      f"{N['eventhap_allele']['germline_len']}, same length other bases at {N['eventhap_allele']['germline_seq']} and ALT at {N['eventhap_allele']['ALT']}.",
      f"- **Normal carries ALT** ({pfu['normal carries ALT']}): the carrying haplotype is {', '.join(f'{k} {v}' for k, v in N['ncalt_hap'].most_common())}; "
      f"{N['nca']['rgn']['Extreme']} are in SMaHT Extreme regions; VAF_Ill >= 0.25 at {N['nca']['vaf']['0.25-0.4'] + N['nca']['vaf']['>=0.4']}; assembly and dipcall both give "
      f"ALT on a haplotype at {N['nca']['asm_dip']}; {N['nca']['flagged']} carry a c3 flag (`nf_flags`), {nca_dup} a duplicated-region sign (MAPQ < 20 or an other-scaffold hit). "
      f"The other haplotype's change to the ALT: {', '.join(f'{k} {v}' for k, v in N['ncalt_other'].most_common())}; at the {N['ncalt_other']['none']} 'none' loci the normal is heterozygous for the truth allele, a non-repeat "
      f"INDEL (one haplotype ALT, the other one non-repeat INDEL away; the other haplotype is GRCh38 at {N['nca_other_ref']} of all {pfu['normal carries ALT']}).",
      f"- **Complex** ({pfu['complex']}; {N['complex_rgn']['Extreme']} in Extreme regions): the other haplotype gives {', '.join(f'{k} {v}' for k, v in N['complex_other'].most_common())}; "
      f"d9 match {', '.join(f'{k} {v}' for k, v in cx_match.most_common())}.",
      f"- **The odd labels come with a germline difference in the array.** {odd_germ} of the {len(odd)} other-repeat, SNV-in-repeat and complex INDELs "
      "have a germline allele other than GRCh38 in the array. A = GRCh38 + truth keeps GRCh38's bases there, so the change from "
      "the patient's haplotype to A also reverts germline differences. Examples: 949 (`complex`; a homozygous germline C>G two "
      "bases from a (TG)n insertion, the c2 `with_germline` allele; hapX already has the +14 length), 2719 (`in_STR`; GRCh38 + GTGT "
      "keeps GRCh38's C where the patient has G, so the event becomes a 2-bp TC insertion inside (GT)n), 1996 (`SNV in STR`; "
      "hapX is (CT)13(AT)8, the length of GRCh38 + truth).",
      f"- **Two readings.** At {altlen_n} absorbed INDELs a COLO829BL haplotype has exactly the ALT length but other bases than A "
      f"({', '.join(f'`{k}` {v}' for k, v in N['altlen'].most_common())}; column `alt_len_haps`). (1) The tumor changed the patient's "
      "repeat by one unit and the truth record was written against GRCh38; or (2) the truth is a germline length allele the "
      "normal already carries, like the 'normal carries ALT' group. No tumor assembly was used, so the two are not told apart: "
      f"the germline-like share of the absorbed INDELs is between {p_(pfu['normal carries ALT'], nI)} ({pfu['normal carries ALT']}) and {p_(pfu['normal carries ALT'] + altlen_n, nI)} ({pfu['normal carries ALT'] + altlen_n}).",
      '', '### GRCh38 frame (b1 rule, per perfect set)', '', t_('2e'), '',
      f"- {g38u['none']} of the union ({p_(g38u['none'], nI)}) look non-repeat in GRCh38. In the patient frame {g38n_unit} of them are whole-unit repeat changes, "
      f"{g38n_near} sit in or next to a repeat (`in_STR`, imperfect VNTR), {g38n_snvrep} are substitutions in a repeat, "
      f"{g38none['normal carries ALT']} normal carries ALT, {g38none['complex']} complex, {g38none['single hap'] + g38none['NA']} unresolved and {g38none['none'] + g38none['SNV']} non-repeat "
      f"({', '.join(f'{k} {v}' for k, v in g38none.most_common())}). The patient's germline repeat length or an interruption hides the "
      "repeat in GRCh38. Union crosstab (`repeat_context.md` 2g):", '', t_('2g'), '',
      '## 3. Background: absorption rate by context, all chr1-22 truth INDELs', '',
      f"Patient frame (c3 covers all {allk['INDEL']:,} truth INDELs):", '', t_('3b'), '',
      f"- Same by label (`repeat_context.md` 3a): `HP>=7` {r_(abl['HP>=7'], bgl['HP>=7'])}, `STR` {r_(abl['STR2-6>=3copies'], bgl['STR2-6>=3copies'])}, "
      f"`none` {r_(abl['none'], bgl['none'])}, `normal carries ALT` {r_(abl['normal carries ALT'], bgl['normal carries ALT'])}.",
      '', 'GRCh38 frame:', '', t_('3c'), '',
      '## 4. SNVs', '',
      f"- {nS} absorbed SNVs: S3 {', '.join(f'{s3[p]} ({p})' for p in PLAT)}; S4b {', '.join(str(sum(v for k, v in q[p]['SNV']['cls'].items() if k.startswith('S4b_'))) for p in PLAT)}; "
      f"d9 match exact {N['snv_match']['exact']}, closest {N['snv_match']['closest']}.",
      f"- GRCh38 context: `none` {N['g38_snv']['union']['none']}, `STR` {N['g38_snv']['union']['STR2-6>=3copies']}, `HP4-6` {N['g38_snv']['union']['HP4-6']}, `HP>=7` {N['g38_snv']['union']['HP>=7']}. "
      f"Patient frame: `SNV` (non-repeat) {snv_pf['SNV']}, substitution in a repeat {snv_rep_pf}, `complex` {snv_pf['complex']}, `in_STR` INDEL event {snv_pf['in_STR']}, normal carries ALT {snv_pf['normal carries ALT']}.",
      f"- CpG transitions: {sum(v for (a, b, c), v in xc.items() if b == 'CpG Ti')} of {nS}; {snv_nonrep_cpg} of the {snv_nonrep} non-repeat ones.", '',
      'All truth SNVs, GRCh38 frame:', '', t_('3d'), '', t_('3h'), '',
      f"- The truth SNVs are UV-dominated (C>T {spb['C>T']:,} of {allk['SNV']:,}) and almost never absorbed ({100 * spa['C>T'] / spb['C>T']:.2f}% of C>T). The absorbed "
      "ones are enriched for CpG transitions, the commonest population SNPs, which the d9 graph keeps as SNP bubbles (S3).", '',
      'Absorbed SNVs, GRCh38 class x CpG transition x patient frame:', '', t_('3i'), '',
      '## 5. Region, VAF and the normal\'s genotype', '',
      'Cells are absorbed (union) / all chr1-22 truth INDELs of the patient-frame group and stratum.', '',
      '### SMaHT region', '', t_('4a'), '', 'Homopolymer one-base changes by length and region:', '', t_('4d'), '',
      f"- SNVs: {'; '.join(f'{s} {r_(snr[1][s], snr[0][s])}' for s in RGN)}.",
      '', '### Truth VAF (Illumina)', '', t_('4b'), '',
      '### Normal-frame germline status', '', t_('4c'), '',
      f"- Overall, `germline other length` loci are absorbed more often ({p_(gat['germline other length'], gbt['germline other length'])}) than `both REF` loci "
      f"({p_(gat['both REF'], gbt['both REF'])}), but much of that is class mix: `both REF` holds {gb[('both REF', 'non-repeat')]} of the {sum(gb[(g, 'non-repeat')] for g in GERM)} non-repeat INDELs. Within a class the "
      f"germline length allele matters for STRs ({p_(ga[('germline other length', 'STR')], gb[('germline other length', 'STR')])} vs {p_(ga[('both REF', 'STR')], gb[('both REF', 'STR')])}) and "
      f"other repeats ({p_(ga[('germline other length', 'other repeat')], gb[('germline other length', 'other repeat')])} vs {p_(ga[('both REF', 'other repeat')], gb[('both REF', 'other repeat')])}), hardly for homopolymers "
      f"({p_(ga[('germline other length', 'HP>=7')], gb[('germline other length', 'HP>=7')])} vs {p_(ga[('both REF', 'HP>=7')], gb[('both REF', 'HP>=7')])}). "
      f"None of the truth INDELs labelled 'SNV in repeat' or 'complex' has both haplotypes equal to GRCh38 in the array ({gref_snv} / {gref_cx} such loci).", '',
      '## 6. Side by side with HG008T', '',
      'HG008T numbers are recomputed by `c4_tables.py` from the HG008 files (checked against `indel_repeat_context.md`). '
      'HG008T patient frame = the GIAB INFO event in the HG008-N v6.2 event haplotype, classed with the b5 rule. For a same-rule '
      'comparison every rate and share below uses the b5 rule for COLO829T too (the tables show the c3 rule beside it). The '
      'events themselves still differ in origin: COLO829T derives them from the closest normal haplotype, HG008T takes the INFO event.', '',
      t_('5a'), '', t_('5b'), '', t_('5d'), '', '### Background rates', '', t_('5e'), '', t_('5h'), '', t_('5f'), '',
      f"- **Shares (b5 rule on both sides).** Homopolymers are {p_(S['cgb']['HP>=7'], nI)} of COLO829T's absorbed INDELs vs {p_(S['hg']['HP>=7'], S['hn'])} of HG008T's; STRs {p_(S['cgb']['STR'], nI)} vs "
      f"{p_(S['hg']['STR'], S['hn'])}. One-unit changes are {p_(S['cu5']['1 unit'], S['nrep_c5'])} vs {p_(S['hun']['1 unit'], S['nrep_h'])} of the repeat events. COLO829T has no "
      f"absorbed homopolymer >= 30 bp ({p_(S['hhp']['>=30'], S['nhp_h'])} of HG008T's `HP>=7` events) and fewer at 20-29 bp ({p_(S['ch5']['20-29'], S['nhp_c5'])} vs {p_(S['hhp']['20-29'], S['nhp_h'])}).",
      f"- **Rates (b5 rule on both sides).** COLO829T minus HG008T, percentage points: one-base homopolymer changes "
      f"{', '.join(f'{d:+.1f} ({b} bp)' for d, b in zip(dhp, HPB))}, one-unit STR changes {dstr:+.1f}. Up to 19 bp the rates are close; at 20-29 bp "
      f"COLO829T is lower ({S['abm5']['HP>=7 1-unit, tract 20-29']} / {S['bgm5']['HP>=7 1-unit, tract 20-29']}), and STRs are absorbed more often. The lower homopolymer share follows the truth set: "
      f"`HP>=7` events are {p_(c_hp_truth5, allk['INDEL'])} of COLO829T truth INDELs and {p_(h_hp_truth, S['hI'])} of HG008T's (c3 rule for COLO829T: {p_(c_hp_truth, allk['INDEL'])}).",
      f"- **SNVs.** The no-candidate SNV counts are similar ({min(q[p]['SNV']['nc'] for p in PLAT)}-{max(q[p]['SNV']['nc'] for p in PLAT)} vs "
      f"{min(S['hq1'][(p, 'SNV')][0] for p in ('PacBio', 'ONT', 'Illumina'))}-{max(S['hq1'][(p, 'SNV')][0] for p in ('PacBio', 'ONT', 'Illumina'))} per platform), but a much smaller share is perfect bypass "
      f"({min(100 * q[p]['SNV']['perfect'] / q[p]['SNV']['nc'] for p in PLAT):.1f}-{max(100 * q[p]['SNV']['perfect'] / q[p]['SNV']['nc'] for p in PLAT):.1f}% vs "
      f"{min(100 * S['hq1'][(p, 'SNV')][1] / S['hq1'][(p, 'SNV')][0] for p in ('PacBio', 'ONT', 'Illumina')):.1f}-{max(100 * S['hq1'][(p, 'SNV')][1] / S['hq1'][(p, 'SNV')][0] for p in ('PacBio', 'ONT', 'Illumina')):.1f}%): "
      f"HG008T's GRCh38 SNVs in homopolymers / STRs are absorbed at {p_(S['has_']['HP>=7'], S['hbs']['HP>=7'])} / {p_(S['has_']['STR2-6>=3copies'], S['hbs']['STR2-6>=3copies'])}, COLO829T's at "
      f"{p_(csa['HP>=7'], csb['HP>=7'])} / {p_(csa['STR2-6>=3copies'], csb['STR2-6>=3copies'])}.",
      f"- **Rules.** With the b5 rule COLO829T has `HP>=7` {S['cgb']['HP>=7']} and `STR` {S['cgb']['STR']} ({diff_b5} absorbed INDELs change label, {in2rep} of them c3 `in_STR` -> b5 `HP>=7` / `STR`).",
      f"- **Germline.** COLO829T 'normal carries ALT' {pfu['normal carries ALT']} ({p_(pfu['normal carries ALT'], nI)}) corresponds to HG008T's {h_present} `HG008N_present_*` INDELs "
      f"({p_(h_present, S['hn'])}), which HG008T's table classes by their INFO event (inside the `HP>=7` / `STR` rows). HG008T has no 'complex' row because its INFO events are simple.", '']
# sections 7-8 (numbers from c5-c8; tables h1-h8 above)
c5w, kk, flr, ind_ = hsum['window_checks'], HN['kmer'], HN['floor'], HN['ind']
pres = hsel(iscat(*PRES)); psub = C(r['sub_class'] for r in pres)
glc = HN['gl']; gl_ab = sum(v for (c, _), v in glc.items() if c == AB)
hh_snv_noex = sum(1 for r in hh if r['category'] == AB and r['kind2'] == 'SNV' and r['sub_class'] != 'exact_allele')
hh_ab_snv = sum(1 for r in hh if r['category'] == AB and r['kind2'] == 'SNV')
hh_closest_snv = sum(1 for r in hpv if r['graph_match'] == 'closest' and r['kind'] == 'SNV')
chm_ex = C(r['category'] for r in HV if r['CHM13_carries_exact'] == 'True')
pon = byid(rd(TR / 'pon' / 'truth_pon_COLO829T.tsv'))
for t in ABS:
    assert (pon[t]['chrom'], pon[t]['vcf_pos'], pon[t]['vcf_ref'], pon[t]['vcf_alt']) == (VS[t]['chrom'], VS[t]['vcf_pos'], VS[t]['vcf_ref'], VS[t]['vcf_alt']), t
paf = lambda r: max([float(x) for x in (r['gnomAD_AF'], r['CoLoRSdb_AF']) if x] or [0.0])
PONF = [('repo rule %', lambda r: r['rule'] not in ('', 'none')), ('pop AF >= 1e-3 %', lambda r: paf(r) >= 1e-3),
        ('pop AF >= 0.01 %', lambda r: paf(r) >= 0.01), ('pop AF >= 0.05 %', lambda r: paf(r) >= 0.05)]
rows = []
for k in ('SNV', 'INDEL'):
    for nm, f in [(f'`{c}`', iscat(c)) for c in HC] + [('all absorbed', lambda r: True)]:
        x = [pon[r['truth_id']] for r in hsel(f, k)]
        rows.append([nm, k, len(x)] + [f'{100 * sum(g(r) for r in x) / len(x):.1f}' if x else '' for _, g in PONF])
    y = [r for r in pon.values() if r['truth_id'] not in ABS and r['kind'] == k and r['chrom'] in AUTO]
    rows.append(['not absorbed (other COLO829T truths)', k, f'{len(y):,}'] + [f'{100 * sum(g(r) for r in y) / len(y):.1f}' for _, g in PONF])
    HN[f'pon_{k}'] = {nm: (sum(g(pon[r['truth_id']]) for r in hsel(iscat(AB), k)), hc[(AB, k)], sum(g(r) for r in y), len(y)) for nm, g in PONF}
tb('h9', 'PoN: % tagged per category vs the non-absorbed truths (chr1-22)', ['group', 'kind', 'n'] + [nm for nm, _ in PONF], rows,
   'repo rule = `scripts/filter_panel_of_normals.py` (allele match; gnomAD / CoLoRSdb AF >= 1e-4; dbSNP non-somatic; 1000G), from '
   '`analysis/tensor_recall_20260930/pon/truth_pon_COLO829T.tsv`; the repo INDEL path has no PoN. pop AF = max(gnomAD, CoLoRSdb).')
pS, pI = HN['pon_SNV'], HN['pon_INDEL']
pq = lambda d, nm: f'{100 * d[nm][0] / d[nm][1]:.1f}% ({100 * d[nm][2] / d[nm][3]:.1f}%)'
low_o = ', '.join(f'{t}: {c}' for t, c in flr['low_other'])
R += ['## 7. HPRC haplotypes, individuals, populations', '',
      f"`c5_hprc_walks.py` (the HG008 s4 code) scans the d9 GFA once (P and W lines, 24 processes) and extracts every haplotype's walk "
      f"between the window anchors. A median of {HN['cov']['all']:g} of the 88 HPRC haplotypes per locus have a complete traversal (SNV "
      f"{HN['cov']['SNV']:g}, INDEL {HN['cov']['INDEL']:g}); {HN['cov']['lt44']} loci have fewer than 44, and truth {', '.join(HN['cov']['zero'])} has none (a 'closest' "
      "locus). `c6_hprc_vcf.py` reads the HPRC genotypes of the d9 and the full-graph VCF. `c7_hprc_membership.py` (the HPRC part of "
      "HG008 s9) counts the carriers at three levels:", '', t_('h3'), '',
      '- **Absent loci by HPRC support** (`sub_class`; individuals and superpopulations from `per_locus_populations.tsv`):', '', t_('h4'), '',
      f"  - Unlike HG008T ({hh_snv_noex} of {hh_ab_snv} absent SNVs without an exact HPRC carrier), {HN['ex']['snv']} of the {hc[(AB, 'SNV')]} absent COLO829T SNVs have one: "
      "they are mostly non-repeat graph SNP alleles (section 4), so the population allele is the same allele.",
      f"  - `element_set_other_allele`: HPRC haplotypes walk all the somatic nodes but spell another allele over the array, mostly another repeat length.",
      f"- **'Rare'** means an exact-allele frequency below 0.20 among complete traversals. {hct[PR]} present loci are rare; {len(HN['pr_zero'])} of them "
      f"({', '.join(HN['pr_zero'])}) has no exact HPRC carrier although HPRC haplotypes walk its nodes. Rarity in large populations is "
      "gnomAD / CoLoRSdb (section 8, PoN table).",
      f"- **d9 floor.** The d9 graph keeps nodes walked by >= ~9 of the 90 haplotypes. Of the {flr['branch']} somatic branch elements, {flr['low']} have a branch "
      f"node that fewer than 9 haplotypes visit in the c5 window rows: {flr['low_chm']} of them lie on CHM13 (a reference path, which the filter keeps); the other "
      f"{len(flr['low_other'])} ({low_o} haplotypes) are window counts, a lower bound (no whole-GFA scan was run for them).",
      "- **Frequencies are conditional** on complete traversals and on the d9 graph: a haplotype cut inside the window does not count. "
      "`hprc_exact_allele_freq_over88` is the lower bound.",
      '- **Cross-checks.**',
      f"  - d9 / full-graph VCF (c6 `carriers_any` vs the GFA allele-level carriers, on GFA-complete haplotypes): identical carrier sets at "
      f"{HN['vcf'][('d9', 'SNV')][0]} / {HN['vcf'][('d9', 'SNV')][1]} SNV and {HN['vcf'][('d9', 'INDEL')][0]} / {HN['vcf'][('d9', 'INDEL')][1]} INDEL loci (d9), "
      f"{HN['vcf'][('full', 'SNV')][0]} / {HN['vcf'][('full', 'SNV')][1]} and {HN['vcf'][('full', 'INDEL')][0]} / {HN['vcf'][('full', 'INDEL')][1]} (full). The VCF counts snarl "
      "alleles: most differences are haplotypes with the record's allele but another allele elsewhere in the array, or loci without a "
      "matching record (`with_germline`, compound repeats; `hprc_tables.md`).",
      f"  - Sequence-level re-derivation (`c7_hprc_membership.py check`: the nearest 16-mers left and right of the allele window that are "
      f"unique in GRCh38, in GRCh38 + truth and in each walk; carrier = the walk's sequence between them equals GRCh38 + truth): identical "
      f"carrier sets at {kk['agree']} of {kk['ok']} anchored loci ({kk['pairs'] - kk['unanch']:,} locus-haplotype pairs compared; {kk['unanch']:,} skipped "
      f"because a germline variant sits in an anchor; {kk['no_anchor']} loci have no unique anchor in the window). The exception is "
      + '; '.join(f'{t} ({c} by c7, {m} by sequence, of {n})' for t, c, m, n in kk['dis']) +
      ": these haplotypes spell GRCh38 + truth between the anchors but compensate beyond the array bracket, so c7 undercounts there "
      "(the category does not change).",
      f"  - Window check: none of the {c5w['window_ref']:,} complete walks that spell GRCh38 over the whole window contains a somatic element "
      f"or counts as a carrier; all {c5w['window_alt']:,} walks that spell GRCh38 + truth over the whole window are exact carriers "
      f"({c5w.get('window_alt_outside_bracket', 0)} only through the whole-window rule, at 29562).",
      f"- **Populations.** The alleles that hide somatic events (absent, exact allele) lean AFR: observed / expected carrier haplotypes AFR "
      f"{oef(oeab, 'AFR')}, AMR {oef(oeab, 'AMR')}, EAS {oef(oeab, 'EAS')}, SAS {oef(oeab, 'SAS')}. AFR is the panel's largest and most diverse "
      f"group. The patient's own germline alleles (present_broad, {hct[PB]} loci) lean slightly AFR too (AFR {oef(oepb, 'AFR')}, AMR "
      f"{oef(oepb, 'AMR')}), unlike HG008T's (AFR 0.92, AMR 1.10). HPRC v1.1 has no EUR haplotype, so the donor's ancestry is not in the panel. "
      f"CHM13 has the exact allele at {chm_ex[AB]} absent and {chm_ex[PB] + chm_ex[PR]} present loci.", '', t_('h5'), '',
      '  Per 1000 Genomes population (`hprc_populations.tsv`; summed over the loci):', '', t_('h6'), '',
      f"  The absent-INDEL exact-allele fraction is {HN['pop_rng']['AFR'][0]:.3f}-{HN['pop_rng']['AFR'][1]:.3f} in the AFR populations, "
      f"{HN['pop_rng']['AMR'][0]:.3f}-{HN['pop_rng']['AMR'][1]:.3f} AMR, {HN['pop_rng']['EAS'][0]:.3f}-{HN['pop_rng']['EAS'][1]:.3f} EAS, {HN['pop_rng']['SAS'][0]:.3f} SAS "
      "(`populations.md`: loci per superpopulation, patterns, shared / private, one-superpopulation loci by population).",
      f"- **Per individual** (`hprc_individuals.tsv`). Each HPRC haplotype carries the exact allele of a median {ind_['med']:g} absent INDELs "
      f"({100 * ind_['share']:.1f}% of those it traverses completely; {100 * ind_['med'] / ind_['n_ai']:.1f}% of all {ind_['n_ai']}) and {ind_['snv']:g} absent SNVs. "
      f"By superpopulation: " + ', '.join(f"{s} {ind_['sp'][s][0]:g} ({100 * ind_['sp'][s][1]:.1f}%)" for s in SPP) +
      f". The range is {ex_ai(ind_['lo'])} ({ind_['lo']['sample']}#{ind_['lo']['hap']}, {ind_['lo']['population_code']}) to {ex_ai(ind_['hi'])} "
      f"({ind_['hi']['sample']}#{ind_['hi']['hap']}, {ind_['hi']['population_code']}); CHM13 {ind_['chm'][0]}. `hprc_carriers.tsv.gz` lists the carrier "
      "haplotypes of every element and allele.", '',
      '## 8. Cross-evaluation with the normal: germline filtering or pangenome-induced false negative', '',
      "Categories (`c7_hprc_membership.py`: the HG008 s9 rules, with the COLO829BL whole-array alleles of c3 in place of HG008-N), "
      "RARE_AF = 0.20 on the exact-allele HPRC frequency. Rules apply in this order:", '', t_('h2'), '',
      '| category | interpretation |', '|---|---|',
      "| `normal_present_*` | the patient's germline allele (a COLO829BL haplotype carries GRCh38 + truth over the whole tandem array), "
      "already a graph path. A tumor-only caller sees a germline allele there and would filter it anyway: the graph acts as **germline filtering** |",
      "| `normal_absent_HPRC_other` | **pangenome-induced false negative**: the somatic change recreates an allele of other individuals "
      "(`exact_allele`), or is spelled by their nodes (`element_set_other_allele`, `recombinant_pieces`) |",
      "| `other_ambiguous` | not resolved; mostly no d9 path spells the truth ('closest'; the COLO829T reads were not re-decoded) |", '',
      'Counts per category (present_broad / present_rare / absent_HPRC_other / other_ambiguous):', '', t_('h1'), '',
      "- **Sub-classes per platform** (chr1-22 and chr1): `hprc_tables.md`.",
      f"- **Present loci** ({len(pres)}): the carrying haplotype is " + ', '.join(f'{k.replace("_only", "")} {v}' for k, v in psub.most_common()) +
      f"; {HN['present_rgn']['Extreme']} are in SMaHT Extreme regions. They inherit the caveats of section 2 ('Normal carries ALT': dipcall "
      f"copy, duplicated-region signs). The one other 'normal carries ALT' locus, {', '.join(HN['nca_am'])}, is a 'closest' locus (rule 1).",
      f"- **Germline-like reading** (section 2, 'Two readings'). At {gl_ab} absent loci (" +
      ', '.join(f'`{s}` {v}' for (c, s), v in sorted(glc.items(), key=lambda kv: -kv[1]) if c == AB) +
      f") a COLO829BL haplotype already has the ALT length with other bases (`germline_like_alt_len`); {sum(glc.values()) - gl_ab} more such loci are ambiguous. "
      f"If these are germline length alleles, the pangenome-induced INDEL losses are between {hc[(AB, 'INDEL')] - gl_ab} and {hc[(AB, 'INDEL')]}.",
      '- **Patient-frame group per category:**', '', t_('h7'), '',
      '### Is the graph helping or hurting?', '',
      f"- **Germline filtering: {pp_(len(pres), nH)}.** The tumor allele is a COLO829BL germline allele. COLO829T has no truth INFO event, so "
      "whether the tumor's other haplotype changed to it (HG008T's case: local conversion or recurrent slippage) or the truth is a "
      "germline allele is not decided here; either way a tumor-only caller sees a germline allele, so the graph costs nothing a "
      "tumor-only design could keep.",
      f"- **Interference: {pp_(hct[AB], nH)}.**",
      f"  - INDELs ({hc[(AB, 'INDEL')]}), patient frame: " + ', '.join(f'`{g}` {v}' for g, v in HN['ab_grp']['INDEL'].most_common()) +
      ". Mostly a one-unit slippage that recreates a length allele found in HPRC, which the d9 graph keeps as a branch or skip edge "
      "(section 3: the absorption rate rises with homopolymer length).",
      f"  - SNVs ({hc[(AB, 'SNV')]}): " + ', '.join(f'`{g}` {v}' for g, v in HN['ab_grp']['SNV'].most_common()) +
      f". {HN['ex']['snv']} have an exact HPRC carrier (median frequency {HN['ex']['fr']['SNV']:.3f}): common population SNP alleles (S3, CpG "
      "transitions enriched, section 4).",
      f"- **Ambiguous: {pp_(hct[AM], nH)}**: {HN['am']['closest']} closest ({HN['am']['near']} with a graph allele nearer the truth than GRCh38), "
      f"{HN['am']['nohprc']} with no complete HPRC haplotype walking the elements ({HN['am']['nohprc_chm']} CHM13 only), {HN['am']['off']} off-phase "
      f"germline, {HN['am']['unres']} normal unresolved.",
      f"- **What a linear-reference tumor-only caller would lose anyway (population PoN).** SNV: the repo PoN rule tags "
      f"{pq(pS, 'repo rule %')} of the absent SNVs (in brackets: the other chr1-22 truth SNVs), pop-AF >= 0.05 {pq(pS, 'pop AF >= 0.05 %')}. A linear "
      f"caller would lose most of these SNVs too. INDEL: the repo INDEL path has no PoN; a population-AF filter would tag the absent INDELs "
      f"{100 * pI['pop AF >= 1e-3 %'][0] / pI['pop AF >= 1e-3 %'][1]:.1f}% / {100 * pI['pop AF >= 0.01 %'][0] / pI['pop AF >= 0.01 %'][1]:.1f}% / "
      f"{100 * pI['pop AF >= 0.05 %'][0] / pI['pop AF >= 0.05 %'][1]:.1f}% at AF 1e-3 / 0.01 / 0.05, against "
      f"{100 * pI['pop AF >= 1e-3 %'][2] / pI['pop AF >= 1e-3 %'][3]:.1f}% / {100 * pI['pop AF >= 0.01 %'][2] / pI['pop AF >= 0.01 %'][3]:.1f}% / "
      f"{100 * pI['pop AF >= 0.05 %'][2] / pI['pop AF >= 0.05 %'][3]:.1f}% of the other truth INDELs.", '', t_('h9'), '',
      '### HG008T vs COLO829T', '', t_('h8'), '',
      f"- **The same picture for INDELs.** {p_(hc[(AB, 'INDEL')], nI)} vs {p_(sum(1 for r in hh if r['category'] == AB and r['kind2'] == 'INDEL'), S['hn'])} of the "
      "absorbed INDELs are pangenome-induced false negatives; the exact allele is about as common in HPRC and as widely shared, and "
      "each HPRC haplotype carries it at about 20% of the absent INDELs it traverses on both sides.",
      f"- **More germline filtering in COLO829T** ({p_(hc[(PB, 'INDEL')] + hc[(PR, 'INDEL')], nI)} vs {p_(sum(1 for r in hh if r['category'] in PRES and r['kind2'] == 'INDEL'), S['hn'])} "
      f"of the INDELs): the 'normal carries ALT' loci, {HN['present_rgn']['Extreme']} of {len(pres)} in SMaHT Extreme regions.",
      f"- **SNVs differ.** COLO829T's absent SNVs are graph SNP alleles that HPRC haplotypes carry exactly ({HN['ex']['snv']} of {hc[(AB, 'SNV')]}); "
      f"HG008T's were repeat events, {hh_snv_noex} of {hh_ab_snv} without an exact carrier. COLO829T has fewer "
      f"ambiguous SNVs ({p_(hc[(AM, 'SNV')], nS)} vs {p_(sum(1 for r in hh if r['category'] == AM and r['kind2'] == 'SNV'), S['habs_snv'])}; HG008T had {hh_closest_snv} closest SNVs).",
      f"- **Populations.** The absent alleles lean AFR on both sides (AFR O/E {oeab['AFR'][0] / oeab['AFR'][1]:.2f} vs {ho['AFR'][0] / ho['AFR'][1]:.2f}). "
      f"The patient's germline alleles (present_broad) lean AMR for HG008T ({hob['AMR'][0] / hob['AMR'][1]:.2f}, z {hob['AMR'][2]:+.1f}) but not for COLO829T "
      f"({oepb['AMR'][0] / oepb['AMR'][1]:.2f}, z {oepb['AMR'][2]:+.1f}; {hct[PB]} loci). Neither donor's (European) ancestry is in HPRC v1.1.", '']
R += [f"""## Files

This folder (`columns.tsv` explains every column of `per_variant.tsv` and of the HPRC files, and gives its source file and column):

| File | Content |
|---|---|
| `README.md` | this report (written by `c4_tables.py`) |
| `per_variant.tsv` | one row per absorbed truth allele ({len(pv)} rows, {len(names)} columns): truth and VAF / RGN; per platform perfect flag, status, miss class, read-level reason, read counts; d9 window, match, primary path and its elements (compact), in-event element subtypes; GRCh38-frame class / unit / tract (and the b5 rule), SNV CpG / Ti-Tv / trinucleotide; COLO829BL array, per-hap placement, assembly / dipcall calls, array allele and its distance to ALT; event hap, patient-frame event, unit, tract, class (c3 and b5 rule), group, units changed, homopolymer bin, germline status, mechanism class, the other hap's event, flags; R / A / hapX / hapY array sequences |
| `columns.tsv` | meaning and source of every column of `per_variant.tsv`, `hprc_per_variant.tsv`, `hprc_per_node.tsv`, `hprc_carriers.tsv.gz`, `hprc_individuals.tsv`, `hprc_populations.tsv`, `per_locus_populations.tsv` |
| `repeat_context.md` | every repeat-context table of `c4_tables.py` (2a-2h, 3a-3i, 4a-4d, 5a-5h), printed by the script too (the tables of sections 7-8 are printed only) |
| `hprc_per_variant.tsv` | one row per absorbed truth allele ({len(hv)} rows, {len(next(iter(hv.values())))} columns): perfect flags, d9 match, chosen ALT path and its somatic elements, allele window; HPRC complete / partial / absent haplotypes, exact-allele / element-set / element frequencies (and over 88), carrier individuals, d9 / full VCF AF and agreement, per-superpopulation complete / carrier haplotypes, CHM13, the carrier haplotypes; COLO829BL hapX / hapY array alleles, event hap, pattern, germline-like flag; category, sub_class, sub_reason, interpretation |
| `hprc_per_node.tsv` | one row per (truth allele, somatic element), {len(hnode)} rows: node ids and sequences, GRCh38 interval, subtype, `alone`; HPRC carriers / frequency, per superpopulation, node coverage (d9 floor), CHM13, carrier haplotypes |
| `hprc_carriers.tsv.gz` | one row per (truth_id, level, carrying haplotype with a complete traversal; HPRC and CHM13): `element` rows (= `hprc_per_node` carriers) and `allele` rows (`carries_set_any_path`, `carries_exact_allele` = `hprc_per_variant` `hprc_exact_carrier_haps`); denominators = `hprc_n_complete` / `{{SP}}_complete` |
| `hprc_individuals.tsv` | per HPRC haplotype (and CHM13): loci per category and kind it traverses completely / carries (exact allele, element set) |
| `hprc_populations.tsv` | per 1000 Genomes population: haplotypes, exact-allele carrier haplotypes and complete traversals per category and kind |
| `per_locus_populations.tsv`, `populations.md` | per truth allele: carrier / complete haplotypes per superpopulation and per 1000G population, superpopulation pattern, shared vs private (one superpopulation / population / individual); `populations.md` = the tables per group (+ O/E, per-haplotype summary) |
| `hprc_tables.md` | every table of `c7_hprc_membership.py` (categories per platform and scope, sub-classes, frequency distributions per level, O/E, crosstabs, VCF cross-check, d9 floor; `$D` column names) |
| `audit_decisions.md` | the two audits of the first version (numbers; independent normal-frame check): every finding, the decision and what changed |
| `c0_*.py` ... `c8_populations.py` | the scripts (below) |

Intermediate data, in `$D = /scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005/`:
`variant_set.tsv` (c0), `grch38_context.tsv` (c1, all {len(g38):,} truth alleles), `loci.tsv` (c1b), `graph_paths.tsv` and
`graph_elements.tsv` (c2), `normal_frame.tsv` (c3, {len(nf):,} loci: the absorbed truths + all truth INDELs), `c3_dipcall_phase.tsv`
(GT-order check), `c3_queries.fa`, `c3_hapX.paf` / `c3_hapY.paf` / `c3_*_unloc.paf` and their minimap2 logs;
`prefix_20261005/` (c3 / c4 outputs before the audit fixes; c4 reads its `normal_frame.tsv` for the 'label changed' row of 2h);
`audit/` (the audits' own scripts and tables); `hprc_local_paths.tsv.gz`, `hprc_haplotypes.tsv`, `hprc_window_coverage.tsv` (c5);
`hprc_vcf_alleles.tsv`, `hprc_vcf_haplotypes.tsv.gz`, `hprc_sample_metadata.tsv` (c6); `hprc_per_variant.tsv` (all c7 columns),
`hprc_per_node.tsv`, `hprc_membership.tsv.gz` (every haplotype x locus x element / allele, carriers or not), `hprc_tables.md`,
`hprc_summary.json`, `hprc_spot_checks.txt`, `hprc_check_kmer.tsv` (c7); `hprc_columns.tsv` (c8). Job scripts and logs:
`tmp/graph_absorbed_somatic_colo829t_20261005/`.

Scripts (each has a docstring with inputs, method, outputs and assumptions, and its results at the end):

| Step | Script / job | What it does |
|---|---|---|
| c0 | `c0_variant_set.py` | the {nS + nI} perfect-bypass truth alleles (`variant_set.tsv`) |
| c1 | `c1_grch38_context.py` | GRCh38-frame repeat context of every chr1-22 truth allele, b1 rule (`grch38_context.tsv`); `--check-hg008` compares the port with the HG008 audit tables |
| c1b | `c1b_loci.py` | d9 graph window, anchors and event window per locus (`loci.tsv`) |
| c2 | `c2_graph_paths.py`, `c2_graph_paths.sbatch` | d9 paths spelling the truth ALT; their elements (`graph_paths.tsv`, `graph_elements.tsv`) |
| c3 | `c3_normal_frame.py`, `c3_minimap2.sbatch` | COLO829BL array alleles (ragtag copy, or dipcall's copy where the ragtag one is suspect), imperfect-VNTR check, the patient-frame event and class (`normal_frame.tsv`, `c3_dipcall_phase.tsv`) |
| c4 | `c4_tables.py` | `per_variant.tsv`, `columns.tsv`, `repeat_context.md`, this README (run last: it reads the c7 / c8 outputs) |
| c5 | `c5_hprc_walks.py`, `c5_hprc_walks.sbatch` | HPRC haplotype walks through every window from the d9 GFA (`hprc_local_paths.tsv.gz`, `hprc_haplotypes.tsv`, `hprc_window_coverage.tsv`; HG008 s4) |
| c6 | `c6_hprc_vcf.py` | d9 / full-graph VCF carriers by sequence (`hprc_vcf_alleles.tsv`, `hprc_vcf_haplotypes.tsv.gz`); population metadata copy (HG008 s5) |
| c7 | `c7_hprc_membership.py` | somatic elements, HPRC membership (element / set / exact allele), COLO829BL calls, categories (`hprc_per_variant.tsv`, `hprc_per_node.tsv`, `hprc_membership.tsv.gz`, `hprc_tables.md`, `hprc_summary.json`); `check` = the sequence-level re-derivation (`hprc_check_kmer.tsv`) (HG008 s9, HPRC part) |
| c8 | `c8_populations.py` | the HPRC files of this folder, `per_locus_populations.tsv`, `populations.md`, `$D/hprc_columns.tsv` (HG008 s10 HPRC part + s13) |

## Reproduce

All inputs are read only (the COLO829BL files belong to another user: never write or index there). Use
`/wanglab/jshen/anaconda3/bin/python` with `PYTHONDONTWRITEBYTECODE=1`; minimap2 2.28 at `/opt/apps/minimap2/2.28/minimap2`.

1. `python c0_variant_set.py` (login, seconds); `python c1_grch38_context.py` (login, 469 s, 0.19 GB);
   `python c1b_loci.py` (login, 3.4 min, 0.18 GB).
2. `sbatch tmp/graph_absorbed_somatic_colo829t_20261005/c2_graph_paths.sbatch` (32 CPUs, 3G; 6 min, MaxRSS 1.6 GB).
3. `python c3_normal_frame.py queries`; `sbatch tmp/graph_absorbed_somatic_colo829t_20261005/c3_minimap2.sbatch` (16 CPUs,
   24G; 5 min 40 s, peak RSS 12.7 GB by /usr/bin/time, sacct MaxRSS 7.4 GB); the unlocalized check on the login node:
   `minimap2 -x asm5 -c --secondary=yes -N 10 -t 4 <ragtag ..._hap{{X,Y}}_unlocalized_normalized.fa> $D/c3_queries.fa >
   $D/c3_hap{{X,Y}}_unloc.paf` (16 s, 1.3 GB); then `python c3_normal_frame.py` (login, {C3_RUN}; it also reads dipcall's
   `hap{{1,2}}.paf.gz` and the raw verkko contigs, read only).
4. `sbatch tmp/graph_absorbed_somatic_colo829t_20261005/c5_hprc_walks.sbatch` (24 CPUs, 16G; 8 min 6 s, MaxRSS 10.1 GB).
5. `python c6_hprc_vcf.py` (login, 2.8 min, 2.9 GB).
6. `python c7_hprc_membership.py 6` (login, 1.5 min, 0.16 GB per process; it imports `c2_graph_paths.py` for the d9 node index);
   `python c7_hprc_membership.py check 6` (login, 43 s).
7. `python c8_populations.py` (login, 2 s, 26 MB).
8. `python c4_tables.py` (login, {C4_RUN}).

## Assumptions and caveats

- **The patient frame is derived, not given.** HG008T used the GIAB truth INFO event in the HG008-N assembly; COLO829T's
  truth has none, so c3 takes the COLO829BL haplotype closest to GRCh38 + truth and derives the event. Levenshtein
  charges a k-bp indel k, so at some loci the closest haplotype gives `complex` while the other gives a single event
  (flag `other_hap_single_event`, {sum('other_hap_single_event' in o['nf_flags'] for o in ai)} absorbed INDELs; the alternative is in `other_hap_*`). Ties: {ties['tie_REF_chosen']} absorbed truths
  chose the GRCh38 haplotype, {ties['tie_other']} fell back to hapX.
- **A = GRCh38 + truth.** The patient's germline differences inside the array are not applied to A (HG008T's INFO event
  was already in the patient's frame), so loci with such differences become `complex`, `in_STR` or `SNV in <class>`
  (section 2, {odd_germ} of {len(odd)}). The `with_germline` d9 match ({q['union']['INDEL']['match']['with_germline']} INDELs) is the graph-side view of the same effect.
  At {altlen_n} of them a normal haplotype already has the ALT length: these may be germline length alleles rather than
  somatic unit changes (section 2, 'Two readings').
- **The tumor is not checked.** Only the normal's assembly is used. Whether the tumor haplotype carries exactly GRCh38 +
  truth over the array (HG008T: `novel_ALT` vs `novel_other`) is not tested here.
- **Paralog fix.** Where the ragtag copy of a haplotype is suspect, c3 uses the contig copy dipcall aligned there
  ({n_dcopy_loci} of the {len(nf):,} loci, {n_dcopy_abs} absorbed; the label changed at {n_dcopy_changed} absorbed truths, `repeat_context.md` 2h).
  That dipcall's long colinear alignment holds the patient's own copy is an assumption. After the fix a copy used
  contradicts an isolated phased het SNV of dipcall at {n_phase_bad} loci (before: {n_phase_bad0}). No segmental-duplication track was used.
- **Normal carries ALT.** {pfu['normal carries ALT']} absorbed INDELs, {N['nca']['rgn']['Extreme']} in Extreme regions, mostly on hapX
  ({N['ncalt_hap']['hapX']}). {nca_dup} of them have a duplicated-region sign (a ragtag hit with MAPQ < 20 or a covering hit on
  another scaffold); of the {len(nca_146)} at chr1:146-149 Mb, {nca_146_dup} do. At {nca_dup_dc} of the {nca_dup} the haplotypes come from dipcall's copy and agree with
  dipcall by construction; a paralogous copy cannot be excluded. dipcall agreement is not independent evidence (same assembly).
- **Class rules.** `pf_class` requires the repeat to hold the event's own bases; `pf_class_b5` (the HG008 rule) does not. They
  differ at {diff_b5} absorbed INDELs; section 6 shows both. SNV events use the b5 rule.
- **Germline differences outside the array** are reverted to GRCh38 by one affine-gap alignment before the event is derived;
  with >= 3 such differences (flag `flank_diffs>=3`) the array / flank split is less reliable. In long or imperfect repeats
  the aligner can put repeat-copy indels at the array edge or the anchor, where they are reverted: the haplotype's array
  allele then understates its repeat-length change (e.g. 40275 hapX: array +96, whole region +475). `total_len_change`
  keeps the whole change and flag `flank_len>=4` marks {N['flanklen']} absorbed INDELs ({N['flanklen_ev']} on the event haplotype).
  Event labels at the loci read by hand were not affected.
- **Imperfect repeats.** The arrays are exact-period stretches. A k-mer recurrence check (12-mers 7-150 bp apart; thresholds
  chosen by the audit, not tuned) finds the truth inside a longer imperfect VNTR / minisatellite at {sum(N['ivntr'].values())} absorbed INDELs
  ({', '.join(f'{k} {v}' for k, v in N['ivntr'].most_common())}); only `none` events are relabelled (`imperfect VNTR`). The b5 rule and
  the HG008T side have no such check.
- **Wide anchors.** {sum(N['wide'].values())} absorbed INDELs needed the 2-kb anchor retry ({', '.join(f'{k} {v}' for k, v in N['wide'].most_common())}); their haplotype
  sequences span longer, often repeat-rich stretches.
- **Unresolved.** {len(na_loci)} loci have no haplotype sequence: {'; '.join(f"{t} ({'absorbed' if a else 'background'}; hapX {NAWHY.get(x, x)}, hapY {NAWHY.get(y, y)})" for t, a, x, y in na_loci)}.
  {len(single_loci)} have one haplotype only (`single hap`): {'; '.join(f"{t} ({'absorbed' if a else 'background'}; one-hap label {l})" for t, a, l in single_loci)}.
- **Truths sharing an array.** {len(shared)} absorbed truths have another truth allele inside their tandem array ({'; '.join(f"{t} {nf[t]['kind2']} `{nf[t]['pf_label']}` with {', '.join(us)}" for t, us in shared)}).
  A applies only the truth itself, so the derived event ignores the other allele (40275 + 40276 applied together still
  match neither haplotype).
- **'closest' loci** ({q['union']['INDEL']['match']['closest']} INDEL, {q['union']['SNV']['match']['closest']} SNV): no d9 path spells the truth; the perfectly aligned reads were not re-decoded.
- **HPRC frequencies** are over complete traversals, so they are conditional on the d9 graph and biased upward (a haplotype cut
  inside the window does not count); `hprc_exact_allele_freq_over88` is the lower bound. 'Rare' (< 0.20) is relative to HPRC
  v1.1, which has no EUR haplotype; the COLO829 donor is European. The O/E z-scores ignore haplotype pairing and linkage, so
  they overstate significance.
- **No read data for the HPRC part.** HG008T used re-decoded reads to choose among several ALT paths; COLO829T keeps the c2
  primary at the {hsum['path_choice'].get('no_read_data:primary', 0)} multi-path loci (the element-set level uses every path; the exact-allele level compares
  sequence over the array, so it hardly depends on the path). All {HN['am']['closest']} 'closest' loci stay ambiguous.
- **Event hap in the categories.** The `with_germline` phase check uses c3's derived event hap (the closest haplotype to GRCh38 +
  truth), not a truth INFO event, and 'ALT on the event hap' is no conflict here (the carrying haplotype is the closest one by
  construction). {HN['am']['off']} loci fail the phase check.
- **Exact-allele rule.** c7 adds one rule to the HG008 one (a walk that spells GRCh38 + truth over the whole window counts); it
  adds 1 haplotype at 29562. The sequence-level check still finds 4 more carriers there (c7 counts 10 of 43, the sequence 14).
- **Present = c3 'normal carries ALT'**, so the present categories inherit c3's caveats (dipcall copy, duplicated-region signs;
  {HN['present_rgn']['Extreme']} of {len(pres)} in SMaHT Extreme regions). The {gl_ab} absent loci flagged `germline_like_alt_len` may be germline length
  alleles (section 8).
- **d9 floor.** {len(flr['low_other'])} somatic branch elements not on CHM13 have fewer than 9 visiting haplotypes in the c5 window rows ({low_o});
  these are lower bounds, not checked with a whole-GFA scan.
- **GRCh38 frame quirks** (b1 kept unchanged so the classes match HG008T): an INDEL with no stretch found has tract 0; for a
  unit > 6 the tract counts copies only to the right; a deletion such as ATATAT>A gets unit 5.
- **Counting.** One GRCh38 record = one truth allele, as in HG008T. The SNV background is GRCh38 frame only (c3 covers truth
  INDELs and the absorbed SNVs). No COLO829T truth INDEL has VAF_Ill < 0.1.
"""]
open(A / 'README.md', 'w').write('\n'.join(R))
print(f'README.md written ({sum(x.count(chr(10)) + 1 for x in R)} lines)')

# Results (2026-10-05, after the audit fixes; login node 6.8 s, 0.35 GB; two runs byte-identical): per_variant.tsv 497 x 180,
# columns.tsv 180; 166 sourced columns re-checked against their source files (82,502 cells, 0 differences); label, group,
# unit, HP-bin, background (c3 and b5) and RGN tables recounted with independent code (same; the b5 background equals the
# numbers audit's same-rule table X1); HG008T numbers reproduce indel_repeat_context.md (asserts).
# Union INDELs (430), patient frame (c3 rule): HP>=7 241 (56.0%), STR 83 (19.3%), other repeat 20 (in_STR 17, imperfect VNTR
# 2, VNTR>6 1), SNV in repeat 16, non-repeat 9 (none 3, SNV 6), normal carries ALT 33, complex 27, NA 1; whole-unit repeat
# changes 325 (75.6%), 1 unit 316 / 325; HP tract 7-9 6, 10-14 86, 15-19 94, 20-29 55, >=30 0; a hap with the ALT length
# (not ALT) at 35 more -> germline-like share 33-68 (7.7-15.8%). Background (1,912 truth INDELs), c3: HP 1-unit 5.2 / 20.3 /
# 25.1 / 32.5% by tract bin, STR 1-unit 46.1%, non-repeat 1.3% (3 / 235); b5 (HG008 rule): 4.7 / 20.0 / 25.6 / 31.0%, STR
# 45.8%; HG008T (b5) 4.6 / 21.5 / 28.8 / 38.4%, 37.8% -> deltas +0.0 / -1.5 / -3.2 / -7.3, STR +8.0. HP>=7 share of truth
# INDELs (b5) 60.7% vs 79.1%. RGN Easy 5.9%, Difficult 29.0%, Extreme 29.8%. SNVs 67 (S3 57 / 51 / 50), CpG Ti absorbed
# 1.04% vs 0.10%. Hand checks: 138, 16258, 415, 521, 949, 1996, 2710, 2719 (first version); 2561-2563, 35166, 40329, 29401,
# 39443, 34621, 2570, 6942, 32518, 1166, 2569 (fixes).
# Sections 7-8 (2026-10-05, after c5-c8; login node 11 s, 0.44 GB): per_variant.tsv and repeat_context.md byte-identical to
# the run before; columns.tsv + 189 HPRC rows (370 lines). Union (chr1-22) SNV / INDEL: present 1 / 32, absent_HPRC_other
# 62 / 365 (exact allele 58 / 330), ambiguous 4 / 33; HG008T numbers recomputed and asserted against its README / tables.md
# (present 5 / 105, absent 150 / 1,935; O/E absent exact AFR 1.107, present broad AMR 1.098).
