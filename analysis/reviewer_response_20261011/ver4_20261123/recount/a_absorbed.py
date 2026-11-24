#!/usr/bin/env python3
"""ver4 recount (a): the absorbed-allele passage of R2.3 point 2 / R2.4 / R3.7 on the FINAL E5d classification.

    PYTHONDONTWRITEBYTECODE=1 /wanglab/jshen/anaconda3/bin/python -B a_absorbed.py

Read only. Sources (all existing; nothing is recomputed from reads or assemblies):
  PA   $E5D/per_allele.tsv                      E5d e4_summary.py output: per absorbed allele `category`
                                                (present / absent / ambiguous), `kind2`, HG008T `tumor_verdict`
  NF   graph_absorbed_somatic_colo829t_20261005/normal_frame.tsv   c3: `alt_len_haps`, `hapX_allele`, `hapY_allele`
  AWG  $E5D/alt_with_germline.tsv               E5d step 5 (2026-11-05): haps that carry the truth plus another
                                                germline difference (`asm_new` / `dip_new` = ALT_plus_germline)
  PL   $E5E/per_locus.tsv                       E5e: tumor hifiasm assembly at the COLO829T normal-absent alleles
                                                (`tumor_verdict`, `r37_outcome`)
  A22  E27 v3 colo829t_altlen_tumor_asm/altlen_tumor_asm_22_summary.tsv   the 22 old "same length, other bases" loci,
                                                column `outcome` (tumor assembly, R3.7 rule)
The E5d / E5e classification code is imported unchanged (colo_alleles, hg008t_alleles, asm_verdicts). The pre-step-5
state (the numbers in letter ver3) is reproduced by running the same functions with the step-5 table emptied; the final
state must equal per_allele.tsv allele by allele (asserted).
Writes a_absorbed.md and a_absorbed.json next to this script.
"""
import collections
import csv
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = Path('/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments')
E5D = '/scratch/jshen/data/reviewer_experiments_20261013/E5d_absorbed_normal_d46'
E5E = '/scratch/jshen/data/reviewer_experiments_20261013/E5e_colo829t_tumor_evidence'
NF = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_colo829t_20261005/normal_frame.tsv'
A22 = EXP / 'E27_assessment_20261030/v3/colo829t_altlen_tumor_asm/altlen_tumor_asm_22_summary.tsv'
csv.field_size_limit(10 ** 9)
rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))  # noqa: E731


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def pct(a, b, d=1):
    return f'{100 * a / b:.{d}f}%'


TRUTH = 'tumor carries the truth allele (not in normal)'
NORMAL = 'tumor carries only normal haplotype sequences'


def main():
    d4 = load('e5d_e4', EXP / 'E5d_absorbed_normal_d46/e4_summary.py')
    e4 = load('e5e_e4', EXP / 'E5e_colo829t_tumor_evidence/e4_summary.py')
    sets = rd(f'{E5D}/sets.tsv')
    pa = rd(f'{E5D}/per_allele.tsv')
    S, L = {}, []

    # ------------------------------------------------------------------ final and pre-step-5 classification
    awg = d4.AWG
    Cf, Hf = d4.colo_alleles(sets), d4.hg008t_alleles()            # final (with step 5)
    d4.AWG = {}
    Co, Ho = d4.colo_alleles(sets), d4.hg008t_alleles()            # before step 5 = letter ver3
    d4.AWG = awg
    for sample, A in (('COLO829T', Cf), ('HG008T', Hf)):           # the imported code reproduces per_allele.tsv
        rows = [r for r in pa if r['sample'] == sample]
        assert rows and all(A[r['truth_id']]['category'] == r['category'] for r in rows), sample
    S['check_code_equals_per_allele_tsv'] = True

    def cats(sample, A, kind=None):
        x = [r for r in sets if r['sample'] == sample and r['arm'] == 'd46' and (kind is None or r['kind2'] == kind)]
        c = collections.Counter(A[r['truth_id']]['category'] for r in x)
        return dict(n=len(x), present=c['present'], absent=c['absent'], ambiguous=c['ambiguous'])

    for sample, Af, Ao in (('HG008T', Hf, Ho), ('COLO829T', Cf, Co)):
        for kind in (None, 'SNV', 'INDEL'):
            S[f'{sample}_d46_{kind or "all"}_final'] = cats(sample, Af, kind)
            S[f'{sample}_d46_{kind or "all"}_before_step5'] = cats(sample, Ao, kind)
    # per_allele.tsv itself, counted directly (the file and column behind every final number)
    for sample in ('HG008T', 'COLO829T'):
        for kind in ('SNV', 'INDEL'):
            c = collections.Counter(r['category'] for r in pa if r['sample'] == sample and r['arm'] == 'd46' and r['kind2'] == kind)
            f = S[f'{sample}_d46_{kind}_final']
            assert (c['present'], c['absent'], c['ambiguous']) == (f['present'], f['absent'], f['ambiguous'])

    # ------------------------------------------------------------------ HG008T: tumor assembly at normal-absent INDELs
    def hg_tumor(A):
        ab = [r for r in sets if r['sample'] == 'HG008T' and r['arm'] == 'd46' and r['kind2'] == 'INDEL' and A[r['truth_id']]['category'] == 'absent']
        tv = collections.Counter(A[r['truth_id']]['tumor_verdict'] for r in ab)
        return dict(absent_INDEL=len(ab), novel_ALT=tv['novel_ALT'], novel_any=tv['novel_ALT'] + tv['novel_other'])
    S['HG008T_tumor_final'], S['HG008T_tumor_before_step5'] = hg_tumor(Hf), hg_tumor(Ho)
    x = [r for r in pa if r['sample'] == 'HG008T' and r['arm'] == 'd46' and r['kind2'] == 'INDEL' and r['category'] == 'absent']
    assert len(x) == S['HG008T_tumor_final']['absent_INDEL']
    assert sum(r['tumor_verdict'] == 'novel_ALT' for r in x) == S['HG008T_tumor_final']['novel_ALT']

    # ------------------------------------------------------------------ COLO829T INDELs
    cind = [r for r in sets if r['sample'] == 'COLO829T' and r['arm'] == 'd46' and r['kind2'] == 'INDEL']
    ids = [r['truth_id'] for r in cind]
    assert len(ids) == len(set(ids))
    gl_old = [t for t in ids if Co[t]['germline_like']]
    gl_new = [t for t in ids if Cf[t]['germline_like']]
    a22 = {r['truth_id']: r for r in rd(A22)}
    assert set(gl_old) == set(a22), 'the 22 loci of E27 v3 are the pre-step-5 same-length loci'
    assert set(gl_new) <= set(gl_old)
    moved = sorted((t for t in ids if Co[t]['category'] != Cf[t]['category']), key=int)
    S['COLO829T_INDEL_changed_by_step5'] = {t: f"{Co[t]['category']}->{Cf[t]['category']}" for t in moved}
    S['COLO829T_same_length_before_step5'] = len(gl_old)
    S['COLO829T_same_length_final'] = len(gl_new)
    by = collections.Counter((Cf[t]['category'], a22[t]['outcome']) for t in gl_old)
    S['old22_by_final_category_and_tumor_asm'] = {f'{a}: {b}': n for (a, b), n in sorted(by.items())}
    S['old22_final_category'] = dict(collections.Counter(Cf[t]['category'] for t in gl_old))
    S['same_length_final_category'] = dict(collections.Counter(Cf[t]['category'] for t in gl_new))

    # E5e per_locus.tsv (final absent loci) against the E27 outcome at the still-absent same-length loci
    pl = {r['truth_id']: r for r in rd(f'{E5E}/per_locus.tsv') if r['set'] == 'absent_INDEL'}
    assert set(pl) == {t for t in ids if Cf[t]['category'] == 'absent'}
    for t in gl_new:
        if Cf[t]['category'] == 'absent':
            assert pl[t]['r37_outcome'] == a22[t]['outcome'], t
    oc_new = collections.Counter(a22[t]['outcome'] for t in gl_new)
    n_truth, n_norm = oc_new[TRUTH], oc_new[NORMAL]
    n_unres = len(gl_new) - n_truth - n_norm
    S['same_length_final_tumor_asm'] = dict(oc_new)
    oc_old = collections.Counter(a22[t]['outcome'] for t in gl_old)
    S['same_length_before_tumor_asm'] = dict(oc_old)
    o_truth, o_norm = oc_old[TRUTH], oc_old[NORMAL]
    o_unres = len(gl_old) - o_truth - o_norm
    pres_new, pres_old = S['COLO829T_d46_INDEL_final']['present'], S['COLO829T_d46_INDEL_before_step5']['present']
    n_ind = len(ids)
    S['germline_like_final'] = dict(present=pres_new, normal_only=n_norm, unresolved=n_unres, lo=pres_new + n_norm,
                                    hi=pres_new + n_norm + n_unres, of=n_ind)
    S['germline_like_before_step5'] = dict(present=pres_old, normal_only=o_norm, unresolved=o_unres, lo=pres_old + o_norm,
                                           hi=pres_old + o_norm + o_unres, of=n_ind)
    S['upper_bound_E5d_table'] = dict(final=pres_new + len(gl_new), before=pres_old + len(gl_old), of=n_ind)

    # tumor assembly at the loci that step 5 moved to 'present' (E5e hits cover the 319 old absent loci)
    Vm, _ = e4.asm_verdicts([dict(truth_id=t, in_nf='True', set='absent_INDEL') for t in moved])
    S['moved_to_present_tumor_asm'] = {t: dict(verdict=Vm[t]['tumor_verdict'], r37=Vm[t]['r37_outcome'],
                                               was_same_length=t in gl_old) for t in moved}
    # R3.7 / R2.4: tumor assembly = GRCh38 + truth at the normal-absent INDELs, final and before step 5
    tv = collections.Counter(r['tumor_verdict'] for r in pl.values())
    r37 = collections.Counter(r['r37_outcome'] for r in pl.values())
    S['COLO829T_tumor_final'] = dict(absent_INDEL=len(pl), novel_ALT=tv['novel_ALT'], novel_any=tv['novel_ALT'] + tv['novel_other'],
                                     verdicts=dict(tv), r37=dict(r37))
    old_abs = [t for t in ids if Co[t]['category'] == 'absent']
    Vo, _ = e4.asm_verdicts([dict(truth_id=t, in_nf='True', set='absent_INDEL') for t in old_abs])
    tvo = collections.Counter(v['tumor_verdict'] for v in Vo.values())
    S['COLO829T_tumor_before_step5'] = dict(absent_INDEL=len(old_abs), novel_ALT=tvo['novel_ALT'], novel_any=tvo['novel_ALT'] + tvo['novel_other'])
    for t in pl:                                                   # the imported verdict code reproduces per_locus.tsv
        assert Vo[t]['tumor_verdict'] == pl[t]['tumor_verdict'] and Vo[t]['r37_outcome'] == pl[t]['r37_outcome'], t

    # ------------------------------------------------------------------ reproduce the letter-ver3 numbers
    b = S['germline_like_before_step5']
    old_ok = dict(
        colo_present_26=S['COLO829T_d46_INDEL_before_step5']['present'] == 26,
        colo_snv_present_1=S['COLO829T_d46_SNV_before_step5']['present'] == 1,
        same_length_22=len(gl_old) == 22, truth_7=o_truth == 7, normal_only_11=o_norm == 11, unresolved_4=o_unres == 4,
        range_37_41=(b['lo'], b['hi']) == (37, 41),
        hg008t_1153_of_1351=(S['HG008T_tumor_before_step5']['novel_ALT'], S['HG008T_tumor_before_step5']['absent_INDEL']) == (1153, 1351),
        colo_226_of_319=(S['COLO829T_tumor_before_step5']['novel_ALT'], S['COLO829T_tumor_before_step5']['absent_INDEL']) == (226, 319),
        hg008t_present_110=S['HG008T_d46_all_before_step5']['present'] == 110,
        colo_present_27=S['COLO829T_d46_all_before_step5']['present'] == 27)
    S['reproduces_letter_ver3'] = old_ok
    assert all(old_ok.values()), old_ok

    # ------------------------------------------------------------------ tables
    f, o = S['germline_like_final'], S['germline_like_before_step5']
    hf, ho = S['HG008T_tumor_final'], S['HG008T_tumor_before_step5']
    cf, co = S['COLO829T_tumor_final'], S['COLO829T_tumor_before_step5']
    Hall, Call = S['HG008T_d46_all_final'], S['COLO829T_d46_all_final']
    Hold, Cold = S['HG008T_d46_all_before_step5'], S['COLO829T_d46_all_before_step5']
    ci, cs = S['COLO829T_d46_INDEL_final'], S['COLO829T_d46_SNV_final']
    L += ['# ver4 recount (a): absorbed alleles, HPRC v2.1 strict subset, final E5d classification', '',
          'Generated by `a_absorbed.py` from the per-allele tables of E5d / E5e (no reads or assemblies re-read). '
          '"letter ver3" = the same code with E5d step 5 switched off (every ver3 number is reproduced, asserted); '
          '"final" = with step 5 (2026-11-05: a normal haplotype with the truth allele plus another germline difference '
          'in the same array counts as carrying it), equal to `per_allele.tsv` allele by allele (asserted).', '',
          '## 1. Sentence numbers', '',
          '| quantity | letter ver3 | final | file and column |', '|---|---|---|---|',
          f"| HG008T absorbed alleles carried by the matched normal | {Hold['present']:,} of {Hold['n']:,} ({pct(Hold['present'], Hold['n'])}) | "
          f"{Hall['present']:,} of {Hall['n']:,} ({pct(Hall['present'], Hall['n'])}) | E5d `per_allele.tsv`, sample HG008T, arm d46, `category` = present |",
          f"| COLO829T absorbed alleles carried by the matched normal | {Cold['present']:,} of {Cold['n']:,} ({pct(Cold['present'], Cold['n'])}) | "
          f"{Call['present']:,} of {Call['n']:,} ({pct(Call['present'], Call['n'])}) | same, sample COLO829T |",
          f"| absent from both normal haplotypes, HG008T / COLO829T | {Hold['absent']:,} ({pct(Hold['absent'], Hold['n'])}) / {Cold['absent']:,} ({pct(Cold['absent'], Cold['n'])}) | "
          f"{Hall['absent']:,} ({pct(Hall['absent'], Hall['n'])}) / {Call['absent']:,} ({pct(Call['absent'], Call['n'])}) | same, `category` = absent |",
          f"| HG008T normal-absent INDELs where the HG008-T assembly carries exactly the benchmark allele | {ho['novel_ALT']:,} of {ho['absent_INDEL']:,} ({pct(ho['novel_ALT'], ho['absent_INDEL'])}) | "
          f"{hf['novel_ALT']:,} of {hf['absent_INDEL']:,} ({pct(hf['novel_ALT'], hf['absent_INDEL'])}) | E5d `per_allele.tsv`, HG008T d46 INDEL, `category` = absent, `tumor_verdict` = novel_ALT |",
          f"| COLO829T absorbed alleles / INDELs | {Cold['n']} / {n_ind} | {Call['n']} / {n_ind} | E5d `per_allele.tsv`, `kind2` |",
          f"| COLO829T INDELs carried by the normal | {pres_old} | {pres_new} | `category` = present, `kind2` = INDEL |",
          f"| COLO829T SNVs carried by the normal | {S['COLO829T_d46_SNV_before_step5']['present']} | {cs['present']} | `category` = present, `kind2` = SNV |",
          f"| COLO829T INDELs absent from the normal / ambiguous | {S['COLO829T_d46_INDEL_before_step5']['absent']} / {S['COLO829T_d46_INDEL_before_step5']['ambiguous']} | {ci['absent']} / {ci['ambiguous']} | `category` |",
          f"| further INDELs where a normal haplotype has an allele of the truth's length with other bases (no haplotype carries the truth) | {len(gl_old)} | {len(gl_new)} | c3 `normal_frame.tsv` `alt_len_haps` non-empty and no hap ALT after `alt_with_germline.tsv` (E5d `colo_alleles`, `germline_like`) |",
          f"| at these sites, tumor hifiasm assembly: truth allele / only normal-haplotype sequences / unresolved | {o_truth} / {o_norm} / {o_unres} | {n_truth} / {n_norm} / {n_unres} | E27 v3 `altlen_tumor_asm_22_summary.tsv` `outcome` (= E5e `per_locus.tsv` `r37_outcome` at the still-absent loci, asserted) |",
          f"| germline-like INDELs: carried by the normal + normal-only (+ unresolved) | {o['lo']}–{o['hi']} of {n_ind} ({pct(o['lo'], n_ind)}–{pct(o['hi'], n_ind)}) | "
          f"**{f['lo']}–{f['hi']} of {n_ind} ({pct(f['lo'], n_ind)}–{pct(f['hi'], n_ind)})** | {pres_new} + {n_norm} (+ {n_unres}) |",
          f"| upper bound of E5d `tables.md` §3 (carried + every same-length locus) | {S['upper_bound_E5d_table']['before']} of {n_ind} ({pct(S['upper_bound_E5d_table']['before'], n_ind)}) | "
          f"{S['upper_bound_E5d_table']['final']} of {n_ind} ({pct(S['upper_bound_E5d_table']['final'], n_ind)}) | {pres_new} + {len(gl_new)} |",
          f"| COLO829T normal-absent INDELs where the tumor assembly carries exactly the benchmark allele (R2.4, R3.7) | {co['novel_ALT']} of {co['absent_INDEL']} ({pct(co['novel_ALT'], co['absent_INDEL'])}) | "
          f"**{cf['novel_ALT']} of {cf['absent_INDEL']} ({pct(cf['novel_ALT'], cf['absent_INDEL'])})** | E5e `per_locus.tsv`, set absent_INDEL, `tumor_verdict` = novel_ALT |",
          '', '## 2. The 22 old "same length, different bases" loci on the final classification', '',
          '| final E5d category | tumor hifiasm assembly (E27 v3 `outcome`) | loci | truth_id |', '|---|---|---:|---|']
    for (cat, oc), n in sorted(by.items()):
        tids = ', '.join(sorted((t for t in gl_old if (Cf[t]['category'], a22[t]['outcome']) == (cat, oc)), key=int))
        L.append(f'| {cat} | {oc} | {n} | {tids} |')
    n_moved22 = sum(Cf[t]['category'] == 'present' for t in gl_old)
    L += ['', f"{n_moved22} of the 22 are now among the {pres_new} INDELs carried by the normal (the normal haplotype has the truth allele plus "
          f"another germline difference); {len(gl_new)} stay normal-absent same-length loci (final categories {S['same_length_final_category']}).", '',
          '## 3. Every COLO829T INDEL whose category step 5 changed', '',
          '| truth_id | chrom:pos REF>ALT | before -> final | one of the 22 | tumor assembly verdict (E5e rule) | R3.7 outcome |', '|---|---|---|---|---|---|']
    pos = {r['truth_id']: r for r in cind}
    for t in moved:
        r, m = pos[t], S['moved_to_present_tumor_asm'][t]
        L.append(f"| {t} | {r['chrom']}:{r['vcf_pos']} {r['vcf_ref']}>{r['vcf_alt']} | {S['COLO829T_INDEL_changed_by_step5'][t]} | "
                 f"{'yes' if m['was_same_length'] else 'no'} | {m['verdict']} | {m['r37']} |")
    n_m = len(moved)
    n_m_norm = sum(S['moved_to_present_tumor_asm'][t]['r37'] == NORMAL for t in moved)
    n_m_alt = sum(S['moved_to_present_tumor_asm'][t]['verdict'] == 'novel_ALT' for t in moved)
    L += ['', f"{n_m} INDELs moved (all absent -> present): {n_moved22} of the 22 and {n_m - n_moved22} other; the tumor assembly carries only "
          f"normal-haplotype sequences at {n_m_norm} of them and a non-normal copy equal to GRCh38 + truth at {n_m_alt}, so the numerator "
          f"{cf['novel_ALT']} is unchanged and only the denominator moves ({co['absent_INDEL']} -> {cf['absent_INDEL']}).", '',
          '## 4. Tumor hifiasm assembly at the final groups (COLO829T INDELs)', '',
          '| group | n | truth allele (not in normal) | only normal-haplotype sequences | other / unresolved | source |', '|---|---:|---:|---:|---:|---|']
    absn = [t for t in ids if Cf[t]['category'] == 'absent']
    rest = [t for t in absn if t not in gl_new]

    def row(name, tids, get, src):
        c = collections.Counter(get(t) for t in tids)
        L.append(f'| {name} | {len(tids)} | {c[TRUTH]} | {c[NORMAL]} | {len(tids) - c[TRUTH] - c[NORMAL]} | {src} |')
    row('normal-absent, all', absn, lambda t: pl[t]['r37_outcome'], 'E5e `per_locus.tsv` `r37_outcome`')
    row('normal-absent, same length as a normal haplotype allele', [t for t in gl_new if t in pl], lambda t: pl[t]['r37_outcome'], 'same')
    row('normal-absent, other', rest, lambda t: pl[t]['r37_outcome'], 'same')
    row('carried by the normal: the loci step 5 moved', moved, lambda t: S['moved_to_present_tumor_asm'][t]['r37'], 'E5e `asm_verdicts` on the E5e hits')
    L += ['', f"The R3.7 per-haplotype outcome 'truth allele' ({r37[TRUTH]}) equals the E5d-style verdict novel_ALT ({tv['novel_ALT']}) at the {len(pl)} loci. "
          f"The {pres_old} INDELs that were already 'carried by the normal' before step 5 and the {ci['ambiguous']} ambiguous INDELs have no tumor-assembly "
          'run (E5e searched the normal-absent loci only).', '',
          '## 5. Checks', '',
          '- E5d `colo_alleles` / `hg008t_alleles` (imported unchanged) give the `category` of `per_allele.tsv` for every allele of both samples and both graphs.',
          f"- With step 5 off, the same code gives every letter-ver3 number: {', '.join(k for k in old_ok)} (all True).",
          '- The pre-step-5 same-length set equals the 22 truth ids of E27 v3 `altlen_tumor_asm_22_summary.tsv`; the final set is a subset of it.',
          '- E5e `asm_verdicts` (imported unchanged) reproduces `tumor_verdict` and `r37_outcome` of `per_locus.tsv` at all '
          f"{len(pl)} final normal-absent INDELs; at the still-absent same-length loci `r37_outcome` equals the E27 v3 `outcome`.",
          f"- `per_locus.tsv` set absent_INDEL has exactly the {len(pl)} truth ids with `category` = absent in `per_allele.tsv`."]
    (HERE / 'a_absorbed.md').write_text('\n'.join(L) + '\n')
    (HERE / 'a_absorbed.json').write_text(json.dumps(S, indent=1) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
