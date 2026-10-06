"""Step 10: the result tables of this folder, cut from the s9 outputs in $D (run after s9_integrate.py).

Inputs ($D): per_variant.tsv (266 columns), per_node.tsv, hprc_membership.tsv.gz, hprc_sample_metadata.tsv, tables.md,
summary.json.
Outputs (this folder):
  per_variant.tsv       one row per truth allele (2,478), the columns of PV below (renamed from the $D names)
  per_node.tsv          one row per (truth allele, somatic graph element) = the nodes / skip edge the ALT reads walk;
                        $D per_node.tsv plus the locus' per-superpopulation complete-traversal counts, with the two
                        HG008-N columns renamed to what they are (giraffe-walk containment, not allele identity)
  hprc_carriers.tsv.gz  one row per (truth allele, level, haplotype) with a complete traversal of the window that
                        carries it (HPRC haplotypes and CHM13): level 'element' (one somatic element) or 'allele' (the
                        whole element set of any ALT path; carries_exact_allele says whether it also spells the exact
                        ALT allele over the tandem array); non-carriers are left out (denominators: per_variant
                        hprc_n_complete and {SP}_complete)
  hprc_individuals.tsv  per haplotype: loci of each category it traverses completely / carries (exact allele, element set)
  hprc_populations.tsv  per 1000G population: haplotypes, and exact-allele carrier haplotypes summed over the loci of
                        each category / summed complete traversals
  columns.tsv           what every column of per_variant.tsv and per_node.tsv means (and its $D name)
  tables.md, summary.json  copies of the s9 tables / numbers (tables.md uses the $D column names; columns.tsv maps them)
"""
import collections, csv, gzip, json, shutil
from pathlib import Path
D = Path('/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001')
A = Path(__file__).resolve().parent
csv.field_size_limit(10**9)
SP = ('AFR', 'AMR', 'EAS', 'SAS')
P3 = ('PacBio', 'ONT', 'Illumina')
PV = [  # ($D column, output column, meaning)
    ('truth_id', 'truth_id', 'row id of the truth allele in analysis/tensor_recall_20260930/per_truth_HG008T.tsv'),
    ('chrom', 'chrom', ''), ('vcf_pos', 'vcf_pos', 'truth VCF POS (1-based)'), ('vcf_ref', 'vcf_ref', ''), ('vcf_alt', 'vcf_alt', ''),
    ('kind2', 'kind', 'SNV or INDEL (GRCh38 representation of the truth)'),
    ('in_bed', 'in_bed', 'truth inside GIAB v0.2 _all.bed'),
    ('in_nogermline_bed', 'in_nogermline_bed', 'truth span inside GIAB v0.2 nogermlineoverlap.bed (stricter)')]
for p in P3:
    PV += [(f'{p}_perfect', f'{p}_perfect', f'{p}: status no_candidate and read-level reason site_bypassed:{{branch,skip_edge}}:no_edit (= in the set for this platform)'),
           (f'{p}_class', f'{p}_miss_class', f'{p}: miss class of the earlier read-level analysis (I1-I7 INDEL, S1-S4c SNV; empty if the platform built the truth)')]
PV += [
    ('match', 'graph_match', 'd9 path spelling the truth ALT haplotype of the window: exact / with_germline (only with nearby HG008-N germline alleles applied too) / closest (none; nearest graph haplotypes kept)'),
    ('germline_used', 'germline_used', 'with_germline: the dipcall alleles that had to be applied (pos:REF>ALT)'),
    ('element_source', 'element_source', 'graph (elements of the chosen exact / with_germline ALT path) / reads (closest: majority walk of the reads) / graph_closest_no_reads'),
    ('somatic_elements', 'somatic_elements', "kept somatic elements: 'B:x:y:>n1>n2' = branch of non-GRCh38 nodes n1,n2 between GRCh38 nodes x and y; 'S:x:y' = edge x->y skipping GRCh38 bases"),
    ('element_subtypes', 'element_subtypes', 'ins_branch / del_skip / snv_branch / mnv_branch / replacing_branch per element (see README section 2)'),
    ('alt_path', 'alt_path', 'the chosen ALT path, first to last anchor node, oriented node ids'),
]
for p in P3:
    PV += [(f'{p}_read_data', f'{p}_read_data', f'{p}: the re-decoded reads (s3) cover this locus'),
           (f'{p}_n_truth_reads', f'{p}_n_truth_reads', f'{p}: no-edit spanning ALT reads whose walk spells GRCh38+ALT (or + the patient germline) over its anchor span'),
           (f'{p}_n_germline_only_reads', f'{p}_n_germline_only_reads', f'{p}: no-edit spanning ALT-like reads whose walk spells only HG008-N germline alleles'),
           (f'{p}_truth_frac_contain_set', f'{p}_truth_frac_contain_set', f'{p}: fraction of truth reads whose walk contains all somatic elements')]
PV += [
    ('truth_read_support', 'truth_read_support', 'yes = >= 1 truth read on a perfect platform with read data / no_truth_spelling_read / no_perfect_platform_read_data'),
    ('hprc_n_complete', 'hprc_n_complete', 'HPRC haplotypes (of 88) whose walk visits both window anchor nodes (complete traversal) = frequency denominator'),
    ('hprc_n_partial', 'hprc_n_partial', 'HPRC haplotypes with only one anchor in the window (walk cut by the d9 filter or a deletion of an anchor)'),
    ('hprc_n_absent', 'hprc_n_absent', 'HPRC haplotypes with no walk through the window'),
    ('hprc_support', 'hprc_support', 'exact_allele (>= 1 complete haplotype spells the exact ALT allele) / element_set_other_allele (walks all elements, other allele) / recombinant_pieces (each element walked by someone, the whole set by no one) / none'),
    ('exact_allele_carriers', 'hprc_exact_allele_carriers', 'complete HPRC haplotypes whose sequence over the tandem array (bracketed by the nearest GRCh38 nodes shared with the ALT path) equals the ALT allele (with_germline: ALT + the germline alleles)'),
    ('exact_allele_freq', 'hprc_exact_allele_freq', 'hprc_exact_allele_carriers / hprc_n_complete (category RARE_AF = 0.20 uses this)'),
    ('exact_allele_freq_all88', 'hprc_exact_allele_freq_over88', 'hprc_exact_allele_carriers / 88 (lower bound)'),
    ('exact_n_individuals', 'hprc_exact_n_individuals', 'HPRC samples with >= 1 exact-allele haplotype'),
    ('exact_n_hom', 'hprc_exact_n_hom', 'HPRC samples with both haplotypes exact-allele carriers'),
    ('set_carriers_any_path', 'hprc_element_set_carriers', 'complete HPRC haplotypes that walk all somatic elements of any ALT path (node level)'),
    ('set_freq_any_path', 'hprc_element_set_freq', 'hprc_element_set_carriers / hprc_n_complete'),
    ('elem_min_carriers', 'hprc_element_min_carriers', 'minimum over the somatic elements of the complete haplotypes walking that element'),
    ('elem_freq_min', 'hprc_element_min_freq', 'hprc_element_min_carriers / hprc_n_complete'),
    ('exact_freq_class', 'hprc_exact_freq_class', '<0.2 / 0.2-0.5 / >=0.5 of hprc_exact_allele_freq'),
    ('d9_vcf_AF', 'd9_vcf_AF', 'AF of the d9 VCF record whose allele alone spells the ALT haplotype (snarl level)'),
    ('full_vcf_AF', 'full_graph_vcf_AF', 'same in the full (unfiltered) HPRC v1.1 vg deconstruct VCF, over all called haplotypes'),
]
for s in SP:
    PV += [(f'{s}_complete', f'{s}_complete', f'{s} HPRC haplotypes with a complete traversal'),
           (f'{s}_exact_carriers', f'{s}_exact_carriers', f'{s} exact-allele carrier haplotypes')]
PV += [
    ('CHM13_exact', 'CHM13_carries_exact', 'CHM13 spells the exact ALT allele (not counted in the HPRC numbers)'),
    ('exact_allele_carrier_haps', 'hprc_exact_carrier_haps', 'the exact-allele carrier haplotypes, sample#hap'),
    ('event_hap', 'HG008N_event_hap', 'the HG008-N v6.2 haplotype the truth INFO HG008Nv62SOMATICVARIANT names (the somatic event is on it)'),
    ('asm_event', 'HG008N_event_INFO', 'truth INFO HG008Nv62SOMATICVARIANT: hap contig:pos(0-based)-REF-ALT = the event in the patient frame'),
    ('asm_event_kind', 'HG008N_event_kind', 'SNV / INS / DEL of that patient-frame event'),
    ('n_truths_same_asm_event', 'n_truths_same_event', 'truth alleles that share this patient-frame event (one event written as several GRCh38 records)'),
    ('hap1_carries_alt', 'HG008N_hap1_carries_ALT', 'final allele-level call for HG008-N hap1: ALT / no / unresolved (whole tandem array, s8b)'),
    ('hap2_carries_alt', 'HG008N_hap2_carries_ALT', 'same for hap2'),
    ('hap1_array_asm_call', 'HG008N_hap1_array_call', 'hap1 sequence between the common unique anchors vs GRCh38 (REF) / GRCh38+truth (ALT) / ALT + the hap\'s own dipcall variants outside the array (ALT_plus_flank) / other / NA (not anchored)'),
    ('hap2_array_asm_call', 'HG008N_hap2_array_call', 'same for hap2'),
    ('hap1_array_len_change', 'HG008N_hap1_array_len_change', 'hap1 array length minus GRCh38 array length (bp)'),
    ('hap2_array_len_change', 'HG008N_hap2_array_len_change', 'same for hap2'),
    ('hap1_walk_elements', 'HG008N_hap1_walks_elements', 'node level: the giraffe walk (vg 1.65, d9) of the hap1 assembly window contains all / some / none of the somatic elements (not_covered: an anchor missing; no_elements). Walking the nodes is not carrying the allele'),
    ('hap2_walk_elements', 'HG008N_hap2_walks_elements', 'same for hap2'),
    ('normal_allele_class', 'HG008N_allele_class', 'the normal\'s own allele: carries_ALT / both_REF / germline_STR_other_length (>= 1 hap has another length in the array) / SNV_site_deleted_on_event_hap / other / unresolved'),
    ('event_hap_pattern', 'event_hap_pattern', 'which HG008-N hap carries the ALT: other_hap_only / other_hap_only:patient_frame (truth INFO event on the event hap = the other hap\'s array allele) / event_hap_only / both / neither / no_event_hap / unresolved'),
    ('T_info_class', 'T_info_class', 'truth INFO event applied to the event hap\'s array: eq_ALT (= GRCh38+truth) / eq_ALT=other_hap / eq_other_hap / novel (a third allele) / event_outside_stretch / NA'),
    ('tumor_verdict', 'HG008T_verdict', 'HG008-T v3.2 contigs over the array: all_eq_hap1/2 (every tumor copy = that normal hap) / novel_ALT (a copy = GRCh38+truth, not in the normal) / novel_other / no_change_* / NA'),
    ('tumor_any_ALT_array', 'HG008T_any_ALT_array', 'a tumor contig equals GRCh38+truth over the array'),
    ('array_lo0', 'array_start0', 'tandem array (s8b) GRCh38 start, 0-based'), ('array_hi0', 'array_end0', 'array end, 0-based exclusive'),
    ('dipcall_relation', 'dipcall_relation', 'HG008-N dipcall PASS records near the truth: identical germline allele / germline variant at the same site (other allele) / germline INDEL or SNV within 10 bp / none within 10 bp'),
    ('dipcall_records', 'dipcall_records', 'those records, pos:REF>ALT(GT hap1|hap2)'),
    ('in_dip_bed', 'in_dip_bed', 'core span +-10 bp inside the dipcall BED'),
    ('PoN_rule', 'PoN_rule', 'which of the four PoNs tag the truth allele under the repo rule (none = not tagged)'),
    ('gnomAD_AF', 'gnomAD_AF', ''), ('CoLoRSdb_AF', 'CoLoRSdb_AF', ''),
    ('category', 'category', 'HG008N_present_broad / HG008N_present_rare / HG008N_absent_HPRC_other / other_ambiguous (README section 5)'),
    ('sub_class', 'sub_class', 'present: event_hap_pattern; absent: hprc_support'),
    ('sub_reason', 'sub_reason', 'other_ambiguous: why'),
    ('interpretation', 'interpretation', 'one sentence per category / sub-class')]

PN_RENAME = {'HG008N_hap1_walked': 'HG008N_hap1_giraffe_walk_contains', 'HG008N_hap2_walked': 'HG008N_hap2_giraffe_walk_contains'}
PN_DESC = {
    'truth_id': 'truth allele', 'element_id': "'B:x:y:>n1>n2' branch / 'S:x:y' skip edge", 'type': 'branch / skip',
    'subtype': 'ins_branch / del_skip / snv_branch / mnv_branch / replacing_branch', 'x': 'GRCh38 node before the element',
    'y': 'GRCh38 node after it', 'nodes': 'the branch nodes, oriented', 'node_seqs': 'their sequences', 'seq': 'branch sequence',
    'replaced': 'GRCh38 bases between x and y that the element replaces (skip: deleted bases)',
    'grch38_start0': 'GRCh38 interval replaced, 0-based start', 'grch38_end0': 'end (exclusive; = start for a pure insertion)',
    'in_event': 'overlaps the event window (s1)', 'role': 'snv (same length) / indel',
    'alone': 'GRCh38 with only this element applied equals: truth (GRCh38+ALT) / germline (GRCh38 + a subset of the 6 nearest HG008-N dipcall alleles) / truth+germline / partial (one piece of a multi-element ALT path)',
    'element_source': 'graph / reads / graph_closest_no_reads', 'category': 'the locus category',
    'hprc_n_complete': 'complete HPRC traversals of the locus window', 'hprc_carriers': 'complete HPRC haplotypes whose walk contains the element (node level)',
    'hprc_freq': 'hprc_carriers / hprc_n_complete', 'hprc_any_row_carriers': 'haplotypes containing the element in any traversal row, partial ones included',
    'min_node_coverage': 'lower bound of haplotypes (HPRC + CHM13 + GRCh38) visiting the least-visited branch node (empty for skip edges)',
    'CHM13': 'CHM13 walk contains it: True / False / partial / absent', 'carrier_haps': 'carrier haplotypes sample#hap',
    'HG008N_hap1_giraffe_walk_contains': 'node-level HG008-N membership: the giraffe walk of the hap1 assembly window contains this element (graph containment; the allele call is per_variant HG008N_hap1_carries_ALT)',
    'HG008N_hap2_giraffe_walk_contains': 'same for hap2'}
for p in P3:
    PN_DESC.update({f'{p}_n_good': f'{p}: no-edit spanning ALT-like reads', f'{p}_n_testable': f'{p}: of those, reads whose walk reaches the element',
                    f'{p}_n_contain': f'{p}: of those, reads whose walk contains it', f'{p}_n_truth': f'{p}: truth-spelling reads',
                    f'{p}_n_truth_testable': f'{p}: truth reads reaching it', f'{p}_n_truth_contain': f'{p}: truth reads containing it'})
for s in SP:
    PN_DESC[f'{s}_carriers'] = f'{s} haplotypes carrying the element'
    PN_DESC[f'{s}_complete'] = f'{s} haplotypes with a complete traversal of the locus (denominator)'

rows = list(csv.DictReader(open(D / 'per_variant.tsv'), delimiter='\t'))
missing = [c for c, _, _ in PV if c not in rows[0]]
assert not missing, missing
with open(A / 'per_variant.tsv', 'w') as w:
    w.write('\t'.join(o for _, o, _ in PV) + '\n')
    for r in rows:
        w.write('\t'.join(r[c] for c, _, _ in PV) + '\n')
by = {r['truth_id']: r for r in rows}
pn = list(csv.DictReader(open(D / 'per_node.tsv'), delimiter='\t'))
pcols = [PN_RENAME.get(c, c) for c in pn[0]]
i = pcols.index('SAS_carriers') + 1
pcols = pcols[:i] + [f'{s}_complete' for s in SP] + pcols[i:]
with open(A / 'per_node.tsv', 'w') as w:
    w.write('\t'.join(pcols) + '\n')
    for r in pn:
        r = {PN_RENAME.get(k, k): v for k, v in r.items()}
        r.update({f'{s}_complete': by[r['truth_id']][f'{s}_complete'] for s in SP})
        w.write('\t'.join(r[c] for c in pcols) + '\n')
assert not [c for c in pcols if c not in PN_DESC], [c for c in pcols if c not in PN_DESC]
with open(A / 'columns.tsv', 'w') as w:
    w.write('file\tcolumn\tD_column\tmeaning\n')
    for c, o, m in PV:
        w.write(f'per_variant.tsv\t{o}\t{c}\t{m}\n')
    inv = {v: k for k, v in PN_RENAME.items()}
    for c in pcols:
        w.write(f'per_node.tsv\t{c}\t{inv.get(c, c)}\t{PN_DESC[c]}\n')
for f in ('tables.md', 'summary.json'):
    shutil.copyfile(D / f, A / f)

cat = {r['truth_id']: (r['category'], r['kind2']) for r in rows}
meta = {r['sample']: r for r in csv.DictReader(open(D / 'hprc_sample_metadata.tsv'), delimiter='\t')}
ind = collections.defaultdict(collections.Counter)            # (sample, hap) -> Counter
with gzip.open(D / 'hprc_membership.tsv.gz', 'rt') as f, gzip.open(A / 'hprc_carriers.tsv.gz', 'wt') as w:
    cols = ['truth_id', 'level', 'element_id', 'sample', 'hap', 'population_code', 'superpopulation', 'carries_element',
            'carries_set_any_path', 'carries_exact_allele']
    w.write('\t'.join(cols) + '\n')
    for r in csv.DictReader(f, delimiter='\t'):
        if r['traversal_status'] != 'complete':
            continue
        allele = r['element_id'] == 'allele'
        k = (r['sample'], r['hap'])
        c, kind = cat[r['truth_id']]
        if allele:
            ind[k][('complete', c, kind)] += 1
            ind[k][('exact', c, kind)] += r['carries_exact_allele'] == 'True'
            ind[k][('set', c, kind)] += r['carries_set_any_path'] == 'True'
        if not (r['carries_element'] == 'True' or r['carries_set_any_path'] == 'True' or r['carries_exact_allele'] == 'True'):
            continue
        w.write('\t'.join([r['truth_id'], 'allele' if allele else 'element', '' if allele else r['element_id'], r['sample'], r['hap'],
                           r['population_code'], r['superpopulation'], r['carries_element'], r['carries_set_any_path'],
                           r['carries_exact_allele']]) + '\n')

CATS = ['HG008N_present_broad', 'HG008N_present_rare', 'HG008N_absent_HPRC_other', 'other_ambiguous']
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
print('per_variant', len(rows), 'columns', len(PV), '| per_node', len(pn), 'columns', len(pcols), '| haplotypes', len(ind), '| populations', len(pop))
