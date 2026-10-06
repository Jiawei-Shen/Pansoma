"""COLO829T step 8: the HPRC result tables of this folder and which HPRC individuals / populations carry the graph allele of
each absorbed truth (port of the HPRC part of analysis/graph_absorbed_somatic_20261001/s10_deliverables.py and of
s13_populations.py, same rules).

Inputs ($D): hprc_per_variant.tsv, hprc_per_node.tsv, hprc_membership.tsv.gz, hprc_tables.md, hprc_summary.json (c7);
hprc_sample_metadata.tsv (c6: HPRC v1.1 sample -> 1000 Genomes population / superpopulation; AFR 23 / AMR 16 / EAS 4 /
SAS 1 samples, no EUR).
Method
  carrier of a locus = a haplotype with a complete traversal that spells the exact ALT allele over the tandem array (c7
  level 'exact'); for the absent loci without any exact carrier (sub_class element_set_other_allele / recombinant_pieces)
  the haplotypes that walk all somatic elements of an ALT path (level 'nodes'). Population / superpopulation frequency =
  carriers / complete haplotypes of that group at those loci (pooled). CHM13 is reported apart (reference, not a
  population). Private: one individual > one population > one superpopulation (first match); shared = >= 2
  superpopulations; none = no carrier. O/E per superpopulation = c7 oe() (expected per locus = carriers x the
  superpopulation's share of the complete haplotypes; hypergeometric z, ignores haplotype pairing and linkage) with a
  95 % interval from 2,000 locus bootstrap resamples (oe_ci).
Outputs (this folder):
  hprc_per_variant.tsv   one row per absorbed truth (497): the columns of PV below (renamed from the $D names)
  hprc_per_node.tsv      one row per (truth, somatic graph element): $D hprc_per_node.tsv + the locus' per-superpopulation
                         complete-traversal counts
  hprc_carriers.tsv.gz   one row per (truth, level, haplotype with a complete traversal) that carries it (HPRC + CHM13):
                         level 'element' (one element) or 'allele' (element set of any ALT path; carries_exact_allele =
                         also the exact allele); non-carriers left out (denominators: hprc_n_complete, {SP}_complete)
  hprc_individuals.tsv   per haplotype: loci of each category / kind it traverses completely / carries (exact, element set)
  hprc_populations.tsv   per 1000G population: haplotypes, exact-allele carrier haplotypes / complete traversals per category
  per_locus_populations.tsv, populations.md   s13: per truth carriers / complete per superpopulation and population,
                         pattern, private flags; the tables per group (+ O/E per superpopulation, per-haplotype summary)
  hprc_tables.md         copy of the c7 tables ($D column names; columns.tsv maps them)
  $D/hprc_columns.tsv    meaning of every column of the files above (c4 appends it to columns.tsv)
Run: python c8_populations.py (login node, seconds).
"""
import collections, csv, gzip, math, random, shutil, statistics, sys
from pathlib import Path
sys.dont_write_bytecode = True
A = Path(__file__).resolve().parent
D = Path('/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005')
csv.field_size_limit(10 ** 9)
SP = ('AFR', 'AMR', 'EAS', 'SAS')
P3 = ('fiberseq', 'ONT', 'Illumina')
CATS = ['normal_present_broad', 'normal_present_rare', 'normal_absent_HPRC_other', 'other_ambiguous']
rd = lambda p: list(csv.DictReader(open(p), delimiter='\t'))


def oe(rows, nm):
    """c7 oe(): {SP: (observed, expected, z)} over the rows (needs hprc_n_complete, {SP}_complete, {SP}_{nm}_carriers)."""
    O, E, V = collections.Counter(), collections.Counter(), collections.Counter()
    for r in rows:
        n, tot = int(r['hprc_n_complete']), sum(int(r[f'{s}_{nm}_carriers']) for s in SP)
        for s in SP:
            q = int(r[f'{s}_complete']) / n
            O[s] += int(r[f'{s}_{nm}_carriers']); E[s] += tot * q
            V[s] += tot * q * (1 - q) * (n - tot) / (n - 1) if n > 1 else 0
    return {s: (O[s], E[s], (O[s] - E[s]) / math.sqrt(V[s]) if V[s] else None) for s in SP}


def oe_ci(rows, nm, B=2000, seed=20261005):
    """c7 oe_ci(): locus bootstrap of oe(), {SP: (2.5 %, 97.5 %) of O/E} over B resamples of the loci."""
    per = [oe([r], nm) for r in rows]
    rng, res = random.Random(seed), {s: [] for s in SP}
    for _ in range(B if rows else 0):
        pick = [per[rng.randrange(len(per))] for _ in per]
        for s in SP:
            e = sum(p[s][1] for p in pick)
            if e:
                res[s].append(sum(p[s][0] for p in pick) / e)
    q = lambda v, f: sorted(v)[min(len(v) - 1, int(f * len(v)))]
    return {s: (q(v, 0.025), q(v, 0.975)) if v else None for s, v in res.items()}


PV = [  # ($D column, output column, meaning)
    ('truth_id', 'truth_id', 'row id of the truth allele in analysis/tensor_recall_20260930/per_truth_COLO829T.tsv (as per_variant.tsv)'),
    ('chrom', 'chrom', ''), ('vcf_pos', 'vcf_pos', 'truth VCF POS (1-based)'), ('vcf_ref', 'vcf_ref', ''), ('vcf_alt', 'vcf_alt', ''),
    ('kind2', 'kind2', 'SNV or INDEL (GRCh38 record kind of the truth)'), ('RGN', 'RGN', 'SMaHT region: Easy / Difficult / Extreme'),
    ('VAF_Ill', 'VAF_Ill', 'truth VAF from Illumina')]
PV += [(f'{p}_perfect', f'{p}_perfect', f'{p}: the truth is perfect bypass on this platform (in its set)') for p in P3]
PV += [
    ('match', 'match', 'c2: a d9 path spells GRCh38 + truth: exact / with_germline (only with nearby COLO829BL dipcall alleles applied too) / closest (none)'),
    ('germline_used', 'germline_used', 'with_germline: the dipcall alleles the path needs (pos:REF>ALT)'),
    ('germline_gt', 'germline_gt', 'their dipcall GT (hapY|hapX)'),
    ('germline_phase', 'germline_phase', 'with_germline: consistent (every used germline allele on the c3 event hap by the GT) / off_event_hap / no_event_hap_or_GT'),
    ('n_candidate_paths', 'n_candidate_paths', 'ALT paths of the d9 graph (c2 n_alt_paths; closest: tied closest paths)'),
    ('path_choice', 'path_choice', 'single_path / no_read_data:primary (several paths, c2 primary kept) / no_read_data:alternative:phase (the first phase-consistent path) / closest_primary'),
    ('element_source', 'element_source', 'graph (elements of the chosen exact / with_germline path) / graph_closest_primary (closest: the c2 closest primary)'),
    ('somatic_elements', 'somatic_elements', "kept somatic elements: 'B:x:y:>n1>n2' = branch of non-GRCh38 nodes between GRCh38 nodes x and y; 'S:x:y' = edge x->y skipping GRCh38 bases"),
    ('n_somatic_elements', 'n_somatic_elements', 'number of kept somatic elements'),
    ('element_subtypes', 'element_subtypes', 'subtype per kept element (ins_branch / del_skip / snv_branch / mnv_branch / replacing_branch / back_branch / back_edge)'),
    ('element_alone', 'element_alone', 'per kept element, GRCh38 with only that element applied equals: truth / germline (a subset of the 6 nearest COLO829BL dipcall alleles) / truth+germline / partial'),
    ('dropped_germline_elements', 'dropped_germline_elements', 'path elements dropped because alone they spell the path\'s germline alleles (with_germline) / nearby germline (closest)'),
    ('closest_kind', 'closest_kind', 'closest: graph_allele_nearer_truth / graph_allele_not_truth (Levenshtein of the kept elements applied to GRCh38)'),
    ('alt_path', 'alt_path', 'the chosen ALT path, first to last anchor node, oriented node ids'),
    ('allele_lo0', 'allele_start0', 'allele-level window start (c3 tandem array united with the event window and the element spans, clipped to the graph window), 0-based'),
    ('allele_hi0', 'allele_end0', 'its end'),
    ('hprc_n_complete', 'hprc_n_complete', 'HPRC haplotypes (of 88) whose walk visits both window anchor nodes (complete traversal) = frequency denominator'),
    ('hprc_n_partial', 'hprc_n_partial', 'HPRC haplotypes with only part of the window (walk cut by the d9 filter or leaving the window)'),
    ('hprc_n_absent', 'hprc_n_absent', 'HPRC haplotypes with no walk through the window'),
    ('hprc_support', 'hprc_support', 'first that applies: exact_allele (>= 1 complete haplotype spells the exact ALT allele) / element_set_other_allele (walks all elements of any candidate ALT path, other allele) / recombinant_pieces (each element of the chosen path walked by someone, a whole set by no one) / none'),
    ('exact_allele_carriers', 'hprc_exact_allele_carriers', 'complete HPRC haplotypes whose sequence over the tandem array (bracketed by the nearest GRCh38 nodes shared with the ALT path; or mapped onto GRCh38 [allele_start0, allele_end0), germline variants outside it ignored), or over the whole window, equals GRCh38 + truth (with_germline: or the path allele = + the germline alleles)'),
    ('exact_allele_freq', 'hprc_exact_allele_freq', 'hprc_exact_allele_carriers / hprc_n_complete (the category RARE_AF = 0.20 uses this)'),
    ('exact_allele_freq_all88', 'hprc_exact_allele_freq_over88', 'hprc_exact_allele_carriers / 88 (lower bound)'),
    ('exact_window_only', 'hprc_exact_window_only', 'carriers found only by the whole-window rule (compensating graph edits outside the array bracket)'),
    ('exact_core_only', 'hprc_exact_core_only', 'carriers found only by the core rule (sequence mapped onto GRCh38 [allele_start0, allele_end0): a germline variant just outside it moved the node bracket)'),
    ('exact_carriers_multicopy', 'hprc_exact_carriers_multicopy', 'exact carriers with >= 2 complete traversals of the window'),
    ('exact_carriers_alt_ref_copy', 'hprc_exact_carriers_alt_ref_copy', 'exact carriers with another complete traversal spelling GRCh38 over the allele window (an ALT copy and a REF copy)'),
    ('multicopy_psv_like', 'multicopy_psv_like', 'every exact carrier has an ALT copy and a REF copy: a paralogous sequence variant rather than an allele of the locus'),
    ('exact_n_individuals', 'hprc_exact_n_individuals', 'HPRC samples with >= 1 exact-allele haplotype'),
    ('exact_n_hom', 'hprc_exact_n_hom', 'HPRC samples with both haplotypes exact-allele carriers'),
    ('set_carriers_any_path', 'hprc_element_set_carriers', 'complete HPRC haplotypes that walk all somatic elements of any ALT path (node level)'),
    ('set_freq_any_path', 'hprc_element_set_freq', 'hprc_element_set_carriers / hprc_n_complete'),
    ('n_set_carriers_other_allele', 'hprc_set_carriers_other_allele', 'complete haplotypes walking all elements of the chosen path but spelling another allele over the array'),
    ('elem_min_carriers', 'hprc_element_min_carriers', 'minimum over the somatic elements of the complete haplotypes walking that element'),
    ('elem_freq_min', 'hprc_element_min_freq', 'hprc_element_min_carriers / hprc_n_complete'),
    ('elem_carriers_each', 'hprc_element_carriers_each', 'complete haplotypes walking each element, in somatic_elements order'),
    ('exact_freq_class', 'hprc_exact_freq_class', '<0.2 / 0.2-0.5 / >=0.5 / NA of hprc_exact_allele_freq'),
    ('freq_class_agree', 'freq_class_agree', 'element / set / exact / full-VCF frequencies on the same side of 0.20'),
    ('d9_vcf_AF', 'd9_vcf_AF', 'AF of the d9 VCF record whose allele alone spells the ALT haplotype (c6; snarl level, over called haplotypes)'),
    ('full_vcf_AF', 'full_graph_vcf_AF', 'same in the full (unfiltered) HPRC v1.1 vg deconstruct VCF'),
    ('d9_agree', 'd9_vcf_agree', 'd9 VCF carriers_any (on GFA-complete haplotypes) == the GFA allele-level carriers'),
    ('d9_disagreements', 'd9_vcf_disagreements', 'haplotype-level reasons when not: gfa_only:<no_matching_record / vcf_gt_missing / vcf_other_allele>, vcf_only:<gfa_set_other_allele / gfa_other_walk>'),
    ('full_agree', 'full_vcf_agree', 'the same for the full-graph VCF')]
for s in SP:
    PV += [(f'{s}_complete', f'{s}_complete', f'{s} HPRC haplotypes with a complete traversal'),
           (f'{s}_exact_carriers', f'{s}_exact_carriers', f'{s} exact-allele carrier haplotypes'),
           (f'{s}_anypath_carriers', f'{s}_element_set_carriers', f'{s} haplotypes walking all elements of an ALT path')]
PV += [
    ('CHM13_status', 'CHM13_status', 'CHM13 traversal: complete / partial / absent'),
    ('CHM13_exact', 'CHM13_carries_exact', 'CHM13 spells the exact ALT allele (not counted in the HPRC numbers)'),
    ('CHM13_set_any_path', 'CHM13_carries_element_set', 'CHM13 walks all elements of an ALT path'),
    ('exact_allele_carrier_haps', 'hprc_exact_carrier_haps', 'the exact-allele carrier haplotypes, sample#hap'),
    ('hapX_allele', 'COLO829BL_hapX_allele', 'c3 array allele of hapX: REF / ALT / germline_len<+-n> / germline_seq / NA'),
    ('hapY_allele', 'COLO829BL_hapY_allele', 'same for hapY'),
    ('event_hap', 'event_hap', 'c3: the COLO829BL hap whose array allele is closest to GRCh38 + truth'),
    ('pf_label', 'pf_label', 'c3 patient-frame label (README section 2)'),
    ('alt_len_haps', 'alt_len_haps', 'c3: haps whose array allele has the ALT length but other bases'),
    ('normal_pattern', 'COLO829BL_pattern', 'which COLO829BL hap carries the ALT: hapX_only / hapY_only / both (+other_unresolved) / neither / unresolved'),
    ('normal_carries_via', 'COLO829BL_carries_via', "ALT (the hap's array allele = GRCh38 + truth, c3 'normal carries ALT') / path_allele (with_germline: = the ALT path's allele, truth + the path's germline alleles) / empty"),
    ('normal_allele_class', 'COLO829BL_allele_class', 'carries_ALT / germline_other_length / germline_same_length_other_seq / both_REF / one_hap_NA / both_NA'),
    ('germline_like_alt_len', 'germline_like_alt_len', 'alt_len_haps non-empty and no hap carries the ALT: the normal may already carry the ALT length (README section 2, Two readings)'),
    ('array_lo0', 'array_start0', 'c3 tandem array start, 0-based'), ('array_hi0', 'array_end0', 'c3 tandem array end, 0-based exclusive'),
    ('category', 'category', 'normal_present_broad / normal_present_rare / normal_absent_HPRC_other / other_ambiguous (README section 8)'),
    ('sub_class', 'sub_class', 'present: the carrying hap (:path_allele = via the with_germline path allele); absent: hprc_support'),
    ('sub_reason', 'sub_reason', 'other_ambiguous: why'),
    ('interpretation', 'interpretation', 'one sentence per category / sub-class')]
PN_DESC = {
    'truth_id': 'truth allele', 'element_id': "'B:x:y:>n1>n2' branch / 'S:x:y' skip edge", 'type': 'branch / skip',
    'subtype': 'ins_branch / del_skip / snv_branch / mnv_branch / replacing_branch / back_branch / back_edge', 'x': 'GRCh38 node before the element',
    'y': 'GRCh38 node after it', 'nodes': 'the branch nodes, oriented', 'node_seqs': 'their sequences', 'seq': 'branch sequence',
    'replaced': 'GRCh38 bases between x and y that the element replaces (skip: deleted bases; < 0: a back edge re-walking a repeat copy)',
    'grch38_start0': 'GRCh38 interval replaced, 0-based start', 'grch38_end0': 'end (exclusive; = start for a pure insertion)',
    'in_event': 'overlaps the event window (c1b)', 'role': 'snv (same length) / indel',
    'alone': 'GRCh38 with only this element applied equals: truth / germline (a subset of the 6 nearest COLO829BL dipcall alleles) / truth+germline / partial',
    'element_source': 'graph / graph_closest_primary', 'category': 'the locus category',
    'hprc_n_complete': 'complete HPRC traversals of the locus window', 'hprc_carriers': 'complete HPRC haplotypes whose walk contains the element (node level)',
    'hprc_freq': 'hprc_carriers / hprc_n_complete', 'hprc_any_row_carriers': 'HPRC haplotypes containing the element in any traversal row, partial ones included',
    'min_node_coverage': 'lower bound of haplotypes (HPRC + CHM13 + GRCh38) visiting the least-visited branch node (empty for skip edges)',
    'CHM13': 'CHM13 walk contains it: True / False / partial / absent', 'carrier_haps': 'carrier haplotypes sample#hap'}
for s in SP:
    PN_DESC[f'{s}_carriers'] = f'{s} haplotypes carrying the element'
    PN_DESC[f'{s}_complete'] = f'{s} haplotypes with a complete traversal of the locus (denominator)'


def main():
    COL = []                                                       # (file, column, source, meaning) -> $D/hprc_columns.tsv
    rows = rd(D / 'hprc_per_variant.tsv')
    missing = [c for c, _, _ in PV if c not in rows[0]]
    assert not missing, missing
    with open(A / 'hprc_per_variant.tsv', 'w') as w:
        w.write('\t'.join(o for _, o, _ in PV) + '\n')
        for r in rows:
            w.write('\t'.join(r[c] for c, _, _ in PV) + '\n')
    COL += [('hprc_per_variant.tsv', o, f'c7 {c}', m) for c, o, m in PV]
    by = {r['truth_id']: r for r in rows}
    pn = rd(D / 'hprc_per_node.tsv')
    pcols = list(pn[0])
    i = pcols.index('SAS_carriers') + 1
    pcols = pcols[:i] + [f'{s}_complete' for s in SP] + pcols[i:]
    with open(A / 'hprc_per_node.tsv', 'w') as w:
        w.write('\t'.join(pcols) + '\n')
        for r in pn:
            r.update({f'{s}_complete': by[r['truth_id']][f'{s}_complete'] for s in SP})
            w.write('\t'.join(r[c] for c in pcols) + '\n')
    assert not [c for c in pcols if c not in PN_DESC], [c for c in pcols if c not in PN_DESC]
    COL += [('hprc_per_node.tsv', c, f'c7 hprc_per_node {c}', PN_DESC[c]) for c in pcols]
    shutil.copyfile(D / 'hprc_tables.md', A / 'hprc_tables.md')

    cat = {r['truth_id']: (r['category'], r['kind2']) for r in rows}
    meta = {r['sample']: r for r in rd(D / 'hprc_sample_metadata.tsv')}
    ind = collections.defaultdict(collections.Counter)            # (sample, hap) -> Counter
    L = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))   # s13: truth -> level -> Counter
    with gzip.open(D / 'hprc_membership.tsv.gz', 'rt') as f, gzip.open(A / 'hprc_carriers.tsv.gz', 'wt') as w:
        cols = ['truth_id', 'level', 'element_id', 'sample', 'hap', 'population_code', 'superpopulation', 'carries_element',
                'carries_set_any_path', 'carries_exact_allele']
        w.write('\t'.join(cols) + '\n')
        for r in csv.DictReader(f, delimiter='\t'):
            if r['traversal_status'] != 'complete':
                continue
            allele = r['element_id'] == 'allele'
            t, s, k = r['truth_id'], r['sample'], (r['sample'], r['hap'])
            c, kind = cat[t]
            if allele:
                ind[k][('complete', c, kind)] += 1
                ind[k][('exact', c, kind)] += r['carries_exact_allele'] == 'True'
                ind[k][('set', c, kind)] += r['carries_set_any_path'] == 'True'
                if s == 'CHM13':
                    L[t]['chm13']['exact'] += r['carries_exact_allele'] == 'True'
                    L[t]['chm13']['nodes'] += r['carries_set_any_path'] == 'True'
                else:
                    pop, sp = meta[s]['population_code'], meta[s]['superpopulation']
                    for key in (pop, sp):
                        L[t]['complete'][key] += 1
                        L[t]['exact'][key] += r['carries_exact_allele'] == 'True'
                        L[t]['nodes'][key] += r['carries_set_any_path'] == 'True'
                    L[t]['exact_ind'][s] += r['carries_exact_allele'] == 'True'
                    L[t]['nodes_ind'][s] += r['carries_set_any_path'] == 'True'
            if not (r['carries_element'] == 'True' or r['carries_set_any_path'] == 'True' or r['carries_exact_allele'] == 'True'):
                continue
            w.write('\t'.join([t, 'allele' if allele else 'element', '' if allele else r['element_id'], s, r['hap'],
                               r['population_code'], r['superpopulation'], r['carries_element'], r['carries_set_any_path'],
                               r['carries_exact_allele']]) + '\n')
    COL += [('hprc_carriers.tsv.gz', c, 'c7 hprc_membership', m) for c, m in (
        ('truth_id', 'truth allele'), ('level', 'element (one somatic element) / allele (the locus: element set of any ALT path, exact allele)'),
        ('element_id', 'the element (level element)'), ('sample', 'HPRC sample or CHM13'), ('hap', 'haplotype 1 / 2 (CHM13 0)'),
        ('population_code', '1000 Genomes population'), ('superpopulation', 'AFR / AMR / EAS / SAS / CHM13'),
        ('carries_element', 'level element: the complete walk contains the element'),
        ('carries_set_any_path', 'level allele: walks all elements of an ALT path'),
        ('carries_exact_allele', 'level allele: spells the exact ALT allele (= hprc_per_variant hprc_exact_carrier_haps)'))]
    with open(A / 'hprc_individuals.tsv', 'w') as w:
        cols = ['sample', 'hap', 'population_code', 'population_name', 'superpopulation']
        for c in CATS:
            for kind in ('SNV', 'INDEL'):
                cols += [f'{c}_{kind}_complete', f'{c}_{kind}_exact', f'{c}_{kind}_element_set']
        w.write('\t'.join(cols) + '\n')
        for (s, h), cnt in sorted(ind.items()):
            m = meta.get(s, {})
            out = [s, h, m.get('population_code', ''), m.get('population_name', ''), m.get('superpopulation', '')]
            for c in CATS:
                for kind in ('SNV', 'INDEL'):
                    out += [cnt[('complete', c, kind)], cnt[('exact', c, kind)], cnt[('set', c, kind)]]
            w.write('\t'.join(map(str, out)) + '\n')
    COL += [('hprc_individuals.tsv', c, 'c7 hprc_membership', m) for c, m in (
        ('sample', 'HPRC sample (or CHM13)'), ('hap', 'haplotype'), ('population_code', '1000G population'), ('population_name', ''),
        ('superpopulation', ''), ('<category>_<kind>_complete', 'loci of that category and kind the haplotype traverses completely'),
        ('<category>_<kind>_exact', 'of those, loci where it spells the exact ALT allele'),
        ('<category>_<kind>_element_set', 'of those, loci where it walks all elements of an ALT path'))]
    pop = collections.defaultdict(collections.Counter)
    for (s, h), cnt in ind.items():
        m = meta.get(s, {})
        key = (m.get('population_code') or s, m.get('population_name', ''), m.get('superpopulation', ''))
        pop[key]['haplotypes'] += 1
        for (lvl, c, kind), v in cnt.items():
            if lvl in ('complete', 'exact'):
                pop[key][(lvl, c, kind)] += v
    with open(A / 'hprc_populations.tsv', 'w') as w:
        cols = ['population_code', 'population_name', 'superpopulation', 'haplotypes']
        for c in CATS[:3]:
            for kind in ('SNV', 'INDEL'):
                cols += [f'{c}_{kind}_exact_carrier_haps', f'{c}_{kind}_complete_haps', f'{c}_{kind}_carrier_fraction']
        w.write('\t'.join(cols) + '\n')
        for key, cnt in sorted(pop.items(), key=lambda x: (x[0][2], x[0][0])):
            out = list(key) + [cnt['haplotypes']]
            for c in CATS[:3]:
                for kind in ('SNV', 'INDEL'):
                    e, n = cnt[('exact', c, kind)], cnt[('complete', c, kind)]
                    out += [e, n, f'{e / n:.3f}' if n else '']
            w.write('\t'.join(map(str, out)) + '\n')
    COL += [('hprc_populations.tsv', c, 'c7 hprc_membership', m) for c, m in (
        ('population_code', '1000G population (CHM13: the sample)'), ('population_name', ''), ('superpopulation', ''),
        ('haplotypes', 'HPRC haplotypes of the population'),
        ('<category>_<kind>_exact_carrier_haps', 'exact-allele carrier haplotypes summed over the loci of the category and kind'),
        ('<category>_<kind>_complete_haps', 'complete traversals summed over those loci'), ('<category>_<kind>_carrier_fraction', 'their ratio'))]
    print('hprc_per_variant', len(rows), 'x', len(PV), '| hprc_per_node', len(pn), 'x', len(pcols), '| haplotypes', len(ind), '| populations', len(pop))

    # ---------------------------------------------------------------- s13: populations per locus and per group
    POPS = sorted({(r['superpopulation'], r['population_code']) for r in meta.values() if r['superpopulation'] in SP})
    PNAME = {r['population_code']: r['population_name'] for r in meta.values()}
    NHAP = collections.Counter(r['population_code'] for r in meta.values() if r['superpopulation'] in SP for _ in (1, 2))
    assert sum(NHAP.values()) == 88

    def group(r):
        c = r['category']
        if c == 'normal_absent_HPRC_other':
            return 'absent: exact allele in HPRC' if r['sub_class'] == 'exact_allele' else 'absent: only the nodes in HPRC'
        return {'normal_present_broad': 'COLO829BL present, broad', 'normal_present_rare': 'COLO829BL present, rare'}.get(c, 'ambiguous')
    GROUPS = ['COLO829BL present, broad', 'COLO829BL present, rare', 'absent: exact allele in HPRC', 'absent: only the nodes in HPRC', 'ambiguous']
    level = lambda t: 'nodes' if group(by[t]) == 'absent: only the nodes in HPRC' else 'exact'
    pattern = lambda t: '+'.join(s for s in SP if L[t][level(t)][s] > 0) or 'none'
    out_rows = []
    for t, r in by.items():
        lv = level(t)
        pops_c = [p for _, p in POPS if L[t][lv][p] > 0]
        ii = [s for s, n in L[t][lv + '_ind'].items() if n]
        assert sum(L[t]['exact'][s] for s in SP) == (int(r['exact_allele_carriers']) if r['exact_allele_carriers'] else 0), t
        assert sum(L[t]['complete'][s] for s in SP) == int(r['hprc_n_complete']), t
        out_rows.append(dict(truth_id=t, chrom=r['chrom'], vcf_pos=r['vcf_pos'], kind2=r['kind2'], category=r['category'], sub_class=r['sub_class'],
                             group=group(r), carrier_level=lv, sp_pattern=pattern(t), n_superpops=sum(L[t][lv][s] > 0 for s in SP),
                             n_populations=len(pops_c), carrier_populations=','.join(pops_c), n_individuals=len(ii),
                             private=('one individual' if len(ii) == 1 else 'one population' if len(pops_c) == 1 else 'one superpopulation'
                                      if pattern(t) in SP else 'none' if not ii else 'shared'),
                             CHM13=L[t]['chm13'][lv] > 0, **{f'{p}_perfect': r[f'{p}_perfect'] for p in P3},
                             **{f'{s}_carriers': L[t][lv][s] for s in SP}, **{f'{s}_complete': L[t]['complete'][s] for s in SP},
                             **{f'{p}_carriers': L[t][lv][p] for _, p in POPS}, **{f'{p}_complete': L[t]['complete'][p] for _, p in POPS}))
    cols = list(out_rows[0])
    with open(A / 'per_locus_populations.tsv', 'w') as w:
        w.write('\t'.join(cols) + '\n')
        for x in out_rows:
            w.write('\t'.join(str(x[c]) for c in cols) + '\n')
    PLM = dict(truth_id='truth allele', chrom='', vcf_pos='', kind2='SNV / INDEL', category='c7 category', sub_class='c7 sub_class',
               group='COLO829BL present, broad / rare; absent: exact allele in HPRC / only the nodes in HPRC; ambiguous',
               carrier_level="exact (haplotypes spelling the exact ALT allele) / nodes (group 'absent: only the nodes': haplotypes walking all elements of an ALT path)",
               sp_pattern='superpopulations with >= 1 carrier, + joined (none: no carrier)', n_superpops='', n_populations='1000G populations with >= 1 carrier',
               carrier_populations='those populations', n_individuals='HPRC samples with >= 1 carrier haplotype',
               private='one individual / one population / one superpopulation / shared (>= 2 superpopulations) / none (first match)',
               CHM13='CHM13 carries at the same level')
    for c in cols:
        m = PLM.get(c) or (f'{c[:-8]} perfect bypass' if c.endswith('_perfect') else f'{c.split("_")[0]} carrier haplotypes' if c.endswith('_carriers')
                           else f'{c.split("_")[0]} haplotypes with a complete traversal')
        COL.append(('per_locus_populations.tsv', c, 'c8', m))

    out = ['# HPRC populations carrying the graph allele of the absorbed COLO829T truths (c8)', '',
           'Carrier = complete HPRC haplotype spelling the exact ALT allele over the tandem array; for "absent: only the nodes in '
           'HPRC" = haplotype walking all somatic elements of an ALT path (no haplotype has the exact allele there). Panel: '
           + ', '.join(f'{p} {NHAP[p]}' for _, p in POPS) + ' haplotypes (' + ', '.join(f'{s} {sum(NHAP[p] for q, p in POPS if q == s)}' for s in SP) +
           '; no EUR: the COLO829 donor is a 45-year-old white male). Frequencies pooled over the loci (carriers / complete haplotypes). '
           'Loci counts with >= 1 carrier grow with the number of haplotypes of the population, so compare the frequencies.', '']
    K = ('SNV', 'INDEL', 'all')
    sel = lambda g, k: [x for x in out_rows if x['group'] == g and (k == 'all' or x['kind2'] == k)]
    LV = lambda x: L[x['truth_id']][level(x['truth_id'])]
    out += ['## 1. Superpopulations: loci with >= 1 carrier (share) and pooled carrier frequency', '',
            '| group | kind | loci | ' + ' | '.join(f'{s} loci' for s in SP) + ' | ' + ' | '.join(f'{s} freq' for s in SP) + ' | CHM13 loci |',
            '|' + '---|' * (4 + 2 * len(SP))]
    for g in GROUPS:
        for k in K:
            x = sel(g, k)
            if not x:
                continue
            cells = [f"{sum(LV(r)[s] > 0 for r in x)} ({100 * sum(LV(r)[s] > 0 for r in x) / len(x):.0f}%)" for s in SP]
            fr = []
            for s in SP:
                c = sum(LV(r)[s] for r in x); n = sum(L[r['truth_id']]['complete'][s] for r in x)
                fr.append(f'{c / n:.3f}' if n else '')
            out.append(f"| {g} | {k} | {len(x)} | " + ' | '.join(cells) + ' | ' + ' | '.join(fr) + f" | {sum(r['CHM13'] for r in x)} |")
    out.append('')
    out += ['## 2. 1000 Genomes populations: loci with >= 1 carrier / pooled carrier frequency', '',
            '| superpop | population | haplotypes | ' + ' | '.join(f'{g} ({k})' for g in GROUPS[:4] for k in ('SNV', 'INDEL')) + ' |',
            '|' + '---|' * (3 + 8)]
    for s, p in POPS:
        cells = []
        for g in GROUPS[:4]:
            for k in ('SNV', 'INDEL'):
                x = sel(g, k)
                c = sum(LV(r)[p] for r in x); n = sum(L[r['truth_id']]['complete'][p] for r in x)
                cells.append(f"{sum(LV(r)[p] > 0 for r in x)} / {c / n:.3f}" if n else '')
        out.append(f'| {s} | {p} ({PNAME[p]}) | {NHAP[p]} | ' + ' | '.join(cells) + ' |')
    out.append('| | n loci | | ' + ' | '.join(str(len(sel(g, k))) for g in GROUPS[:4] for k in ('SNV', 'INDEL')) + ' |')
    out.append('')
    pats = collections.Counter((x['group'], x['sp_pattern']) for x in out_rows)
    allp = sorted({p for _, p in pats}, key=lambda p: (-sum(v for (g, q), v in pats.items() if q == p), p))
    out += ['## 3. Which superpopulations carry the allele (pattern per locus)', '',
            '| pattern | ' + ' | '.join(GROUPS) + ' |', '|' + '---|' * (1 + len(GROUPS))]
    for p in allp:
        out.append(f'| {p} | ' + ' | '.join(str(pats[(g, p)]) for g in GROUPS) + ' |')
    out.append('')
    pr = collections.Counter((x['group'], x['private']) for x in out_rows)
    out += ['## 4. Shared or private', '', '| | ' + ' | '.join(GROUPS) + ' |', '|' + '---|' * (1 + len(GROUPS))]
    for p in ('shared', 'one superpopulation', 'one population', 'one individual', 'none'):
        out.append(f'| {p} | ' + ' | '.join(str(pr[(g, p)]) for g in GROUPS) + ' |')
    out += ['', 'shared = carriers in >= 2 superpopulations; one superpopulation = >= 2 populations of one superpopulation '
            '(first match: one individual, one population, one superpopulation).', '']
    one = [x for x in out_rows if x['n_superpops'] == 1]
    c5 = collections.Counter((x['group'], x['sp_pattern'], x['carrier_populations']) for x in one)
    out += ['## 5. Loci carried by one superpopulation only: by superpopulation and populations', '',
            '| group | superpop | carrier populations | loci |', '|---|---|---|---|']
    for (g, s, pp), v in sorted(c5.items(), key=lambda kv: (GROUPS.index(kv[0][0]), kv[0][1], -kv[1])):
        out.append(f'| {g} | {s} | {pp} | {v} |')
    out += ['', '## 6. Observed / expected carrier haplotypes per superpopulation (c7 oe(); exact allele / element set)', '',
            'expected per locus = carriers x the superpopulation\'s share of the complete haplotypes at that locus; z hypergeometric, '
            'ignoring haplotype pairing and linkage (overstates significance); after z the 95 % interval of O/E from 2,000 locus '
            'bootstrap resamples (oe_ci, seed 20261005: the loci, not the haplotypes, are resampled).', '',
            '| group | kind | level | loci | carrier haps | ' + ' | '.join(f'{s} O/E (z; 95% CI)' for s in SP) + ' |', '|' + '---|' * (5 + len(SP))]
    for g in GROUPS:
        for k in K:
            x = [by[r['truth_id']] for r in sel(g, k) if by[r['truth_id']]['n_somatic_elements'] != '0' and by[r['truth_id']]['hprc_n_complete'] != '0']
            if not x:
                continue
            for nm, lab in (('exact', 'exact allele'), ('anypath', 'element set')):
                o = oe(x, nm)
                if not sum(v[0] for v in o.values()):
                    continue
                ci = oe_ci(x, nm)
                out.append(f'| {g} | {k} | {lab} | {len(x)} | {sum(v[0] for v in o.values())} | ' +
                           ' | '.join(f'{o[s][0] / o[s][1]:.3f} (z {o[s][2]:+.1f}; {ci[s][0]:.2f}-{ci[s][1]:.2f})' if o[s][1] and o[s][2] is not None and ci[s] else ''
                                      for s in SP) + ' |')
    out += ['', '## 7. Per HPRC haplotype: exact-allele carriers summed over the loci (hprc_individuals.tsv)', '',
            '| superpop | haplotypes | absent INDEL: median carried (range) | median share of its complete traversals | absent SNV: median carried | '
            'present INDEL: median carried |', '|---|---|---|---|---|---|']
    hap = {k: v for k, v in ind.items() if k[0] != 'CHM13'}
    for s in SP + ('all',):
        hs = [v for k, v in hap.items() if s == 'all' or meta[k[0]]['superpopulation'] == s]
        ai = [v[('exact', 'normal_absent_HPRC_other', 'INDEL')] for v in hs]
        sh = [v[('exact', 'normal_absent_HPRC_other', 'INDEL')] / v[('complete', 'normal_absent_HPRC_other', 'INDEL')] for v in hs]
        asn = [v[('exact', 'normal_absent_HPRC_other', 'SNV')] for v in hs]
        prs = [v[('exact', 'normal_present_broad', 'INDEL')] + v[('exact', 'normal_present_rare', 'INDEL')] for v in hs]
        out.append(f'| {s} | {len(hs)} | {statistics.median(ai):g} ({min(ai)}-{max(ai)}) | {statistics.median(sh):.3f} | {statistics.median(asn):g} | {statistics.median(prs):g} |')
    ch = ind[('CHM13', '0')]
    out.append(f"| CHM13 | 1 | {ch[('exact', 'normal_absent_HPRC_other', 'INDEL')]} | "
               f"{ch[('exact', 'normal_absent_HPRC_other', 'INDEL')] / ch[('complete', 'normal_absent_HPRC_other', 'INDEL')]:.3f} | "
               f"{ch[('exact', 'normal_absent_HPRC_other', 'SNV')]} | {ch[('exact', 'normal_present_broad', 'INDEL')] + ch[('exact', 'normal_present_rare', 'INDEL')]} |")
    mx = max(hap.items(), key=lambda kv: kv[1][('exact', 'normal_absent_HPRC_other', 'INDEL')])
    mn = min(hap.items(), key=lambda kv: kv[1][('exact', 'normal_absent_HPRC_other', 'INDEL')])
    out += ['', f"Most / fewest absent INDEL exact alleles: {mx[0][0]}#{mx[0][1]} ({meta[mx[0][0]]['population_code']}) "
            f"{mx[1][('exact', 'normal_absent_HPRC_other', 'INDEL')]} / {mn[0][0]}#{mn[0][1]} ({meta[mn[0][0]]['population_code']}) "
            f"{mn[1][('exact', 'normal_absent_HPRC_other', 'INDEL')]} of {sum(r['category'] == 'normal_absent_HPRC_other' and r['kind2'] == 'INDEL' for r in rows)} absent INDELs."]
    open(A / 'populations.md', 'w').write('\n'.join(out) + '\n')
    with open(D / 'hprc_columns.tsv', 'w') as w:
        w.write('file\tcolumn\tsource\tmeaning\n')
        for x in COL:
            assert not any('\t' in y for y in x)
            w.write('\t'.join(x) + '\n')
    print('\n'.join(out))
    print('per_locus_populations', len(out_rows), 'x', len(cols), '| hprc_columns', len(COL))


if __name__ == '__main__':
    main()

# Results (2026-10-05, after the HPRC audit fixes; login node 5.1 s, 26 MB; $T/c8_populations.login.log):
# hprc_per_variant.tsv 497 x 85, hprc_per_node.tsv 623 x 32, 89 haplotypes (88 HPRC + CHM13), 14 populations (13 + CHM13),
# per_locus_populations.tsv 497 x 52, $D/hprc_columns.tsv 194 rows. Absent loci with an exact HPRC carrier (390): shared by
# >= 2 superpopulations 336, one superpopulation 42, one population 1, one individual 11 (carried in one superpopulation
# only: AFR 42, AMR 11, EAS 1). Pooled exact-allele frequency, absent INDEL: AFR 0.253, AMR 0.179, EAS 0.171, SAS 0.191;
# present broad INDEL: AFR 0.445, AMR 0.367. O/E present broad INDEL AFR 1.076 (z +3.0, locus bootstrap 0.96-1.19), AMR
# 0.902 (0.77-1.04): no clear lean. Per HPRC haplotype: median 62 absent INDEL exact alleles (20.0% of its complete
# traversals), AFR 69, AMR 54, EAS 52, SAS 57.5; range 38 (HG00673#1, CHS) - 90 (HG03453#2, MSL); CHM13 68.
