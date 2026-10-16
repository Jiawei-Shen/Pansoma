"""E5b step 5: classes of somatic truth alleles lost on the pangenome graphs, graph spelling, PoN share; tables.

Inputs (data dir): status_matrix.tsv.gz (s1), spell_d9.tsv.gz / spell_d46.tsv.gz (s4), and the PoN tags of the truth VCFs
(scripts/filter_panel_of_normals.py; since 2026-10-16 the current INDEL rule, AF >= 0.01 in homopolymers >= 8 bp, else
AF >= 1e-4, on the same inputs as before: tmp/reviewer_cpu_20261013/pon_hp8_20261016/truth_pon.sh):
  HG008T   tmp/reviewer_cpu_20261013/pon_hp8_20261016/truth/HG008T.truth.pon.vcf.gz (multi-allelic records split per ALT)
  COLO829T tmp/reviewer_cpu_20261013/pon_hp8_20261016/truth/COLO829T.truth.pon.vcf.gz (sites only)
  HCC1395T tmp/reviewer_cpu_20261013/pon_hp8_20261016/truth/HCC1395T.truth.pon.vcf.gz (SEQC2 v1.2.1 truth)
The tags of the earlier rule (HP >= 7: 0.01, else 1e-3) stay in tmp/graph_absorbed_pon_20261012 and DATA/hcc1395.truth.pon.vcf.gz.
Classes, per platform, for an arm X in (d9, d46) against the GRCh38-only arm (linear):
  lost_X        = a tensor (tensor_representative) on linear, none on X
    absorbed    = X status no_candidate and the ALT is spelled by an X graph path (s4)
    nocand_other= X status no_candidate, the ALT not spelled by X
    filtered    = X status filtered (a candidate below the AF / min-variant cut)
    other_allele= X status tensor_non_representative_allele
  gained_X      = a tensor on X, none on linear
d46 vs d9: d46_only = tensor on d46 not on d9; d9_only = the reverse (sub-classed by the losing arm's status as above).
A truth allele is PoN-tagged when PANSOMA_PON is set on its record (allele-exact). 'after PoN' = not tagged.
Writes tables.md and summary.json (experiment dir), per_truth.tsv.gz (data dir; all alleles) and
lost_alleles.tsv.gz (experiment dir; alleles lost on d9 or d46 on at least one platform).
Run: python s5_tables.py
"""
import csv, gzip, json, os, sys
from collections import Counter, defaultdict
import pysam
sys.path.insert(0, '/scratch/jshen/Github/Pansoma/scripts')
from filter_panel_of_normals import homopolymer_run      # the PoN's own homopolymer rule (library use)
from e5b_common import DATA, FASTA

EXP = os.path.dirname(os.path.abspath(__file__))
PON_DIR = '/scratch/jshen/Github/Pansoma/tmp/reviewer_cpu_20261013/pon_hp8_20261016/truth'
PON = {s: f'{PON_DIR}/{s}.truth.pon.vcf.gz' for s in ('HG008T', 'COLO829T', 'HCC1395T')}
PLATS = {'HG008T': ['Illumina', 'PacBio', 'ONT'], 'COLO829T': ['Illumina', 'fiberseq', 'ONT'], 'HCC1395T': ['Illumina', 'PacBio', 'ONT']}
ABSORBED_D9 = {'HG008T': '/scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_20261001/per_variant.tsv',
               'COLO829T': '/scratch/jshen/Github/Pansoma/analysis/graph_absorbed_somatic_colo829t_20261005/per_variant.tsv'}
csv.field_size_limit(1 << 30)
MISSING = []


def pon_tags(path):
    tags = {}
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            x = line.rstrip('\n').split('\t', 8)
            pon = next((v.split('=', 1)[1] for v in x[7].split(';') if v.startswith('PANSOMA_PON=')), '')
            for a in x[4].split(','):
                tags[(x[0], int(x[1]), x[3].upper(), a.upper())] = pon
    return tags


def load_spell(g):
    out = {}
    if not os.path.exists(f'{DATA}/spell_{g}.tsv.gz'):
        print(f'WARNING: no spell_{g}.tsv.gz yet: every ALT counted as not spelled on {g} (preliminary tables)', file=sys.stderr)
        MISSING.append(g)
        return out
    with gzip.open(f'{DATA}/spell_{g}.tsv.gz', 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            out[(r['sample'], r['truth_id'])] = r
    return out


def cls(status, spelled):
    return {'N': 'absorbed' if spelled else 'nocand_other', 'F': 'filtered', 'A': 'other_allele'}.get(status, 'other_' + status)


def main():
    fa = pysam.FastaFile(FASTA)
    spell = {g: load_spell(g) for g in ('d9', 'd46')}
    rows = list(csv.DictReader(gzip.open(f'{DATA}/status_matrix.tsv.gz', 'rt'), delimiter='\t'))
    tags = {s: pon_tags(p) for s, p in PON.items()}
    absorbed_d9 = {}
    for s, p in ABSORBED_D9.items():
        absorbed_d9[s] = {(r['chrom'], int(r['vcf_pos']), r['vcf_ref'].upper(), r['vcf_alt'].upper()) for r in csv.DictReader(open(p), delimiter='\t')}
    missing_pon = Counter()
    ref_fail = Counter()
    for r in rows:
        s = r['sample']
        k = (r['chrom'], int(r['vcf_pos']), r['vcf_ref'].upper(), r['vcf_alt'].upper())
        if k not in tags[s]:
            missing_pon[s] += 1
        r['pon'] = tags[s].get(k, '')
        r['tagged'] = bool(r['pon'])
        r['kind2'] = 'SNV' if r['kind'] == 'SNP' else 'INDEL'
        r['hp7'] = r['kind2'] == 'INDEL' and homopolymer_run(fa, r['chrom'], int(r['vcf_pos']), r['vcf_ref'], r['vcf_alt']) >= 7
        r['absorbed_d9_reads'] = k in absorbed_d9.get(s, ())
        for g in ('d9', 'd46'):
            sp = spell[g].get((s, r['truth_id']), {})
            r[f'{g}_spelled'] = sp.get('alt_spelled') == '1'
            r[f'{g}_ref_ok'] = sp.get('ref_ok', '')
            r[f'{g}_elements'] = sp.get('elements', '')
            r[f'{g}_alt_min_paths'] = sp.get('alt_min_paths', '')
            if sp.get('status') == 'ok' and sp.get('ref_ok') != '1':
                ref_fail[(s, g, sp.get('info', ''))] += 1
    print('PoN keys missing per sample:', dict(missing_pon), file=sys.stderr)
    print('GRCh38 path not spelled (subgraph check):', dict(ref_fail), file=sys.stderr)

    # per-truth outputs
    base = ['sample', 'truth_id', 'chrom', 'vcf_pos', 'vcf_ref', 'vcf_alt', 'kind', 'passed', 'in_bed']
    stat_cols = [c for c in rows[0] if '_' in c and c.split('_')[-1] in ('d9', 'd46', 'linear') and c not in base]
    extra = ['pon', 'hp7', 'absorbed_d9_reads', 'd9_spelled', 'd46_spelled', 'd9_ref_ok', 'd46_ref_ok', 'd9_elements',
             'd46_elements', 'd9_alt_min_paths', 'd46_alt_min_paths']
    with gzip.open(f'{DATA}/per_truth.tsv.gz', 'wt') as w, gzip.open(f'{EXP}/lost_alleles.tsv.gz', 'wt') as w2:
        w.write('\t'.join(base + stat_cols + extra) + '\n')
        w2.write('\t'.join(base + stat_cols + extra) + '\n')
        for r in rows:
            line = '\t'.join(str(int(r[c]) if isinstance(r[c], bool) else r[c]) for c in base + stat_cols + extra) + '\n'
            w.write(line)
            if any(r.get(f'{p}_{g}', 'T') not in ('T', '') for p in PLATS[r['sample']] for g in ('d9', 'd46')):
                w2.write(line)

    S = {}
    L = ['# E5b tables (generated by s5_tables.py)', '']
    if MISSING:
        L += [f'**PRELIMINARY: graph spelling missing for {MISSING}; "absorbed" columns are 0 and "not spelled" holds every no-candidate allele.**', '']
    pct = lambda a, b: f'{100 * a / b:.1f} %' if b else '-'

    def tagged_str(sub):
        n = len(sub)
        t = sum(r['tagged'] for r in sub)
        return n, t, f'{n:,} ({t:,} tagged, {pct(t, n)}; {n - t:,} after PoN)'

    # ---- Table 1: candidate status per arm (somatic truth alleles with a tensor) ----
    L += ['## 1. Truth alleles with a tensor (tensor_representative), chr1-22, all truth records; and those not PoN-tagged', '',
          'Not PoN-tagged with a tensor = the most a caller on that arm can still report after the PoN (tensor-level ceiling).', '',
          '| sample | platform | kind | n (not tagged) | tensor: linear / d9 / d46 | tensor and not tagged: linear / d9 / d46 | d9 - linear, d46 - linear (not tagged) |',
          '|---|---|---|---|---|---|---|']
    for s in ('HG008T', 'COLO829T', 'HCC1395T'):
        for p in PLATS[s]:
            for kind in ('SNV', 'INDEL'):
                sub = [r for r in rows if r['sample'] == s and r['kind2'] == kind]
                has = {g: bool(sub[0].get(f'{p}_{g}')) for g in ('linear', 'd9', 'd46')}
                v = {g: sum(r.get(f'{p}_{g}') == 'T' for r in sub) for g in has}
                u = {g: sum(r.get(f'{p}_{g}') == 'T' and not r['tagged'] for r in sub) for g in has}
                nt = sum(not r['tagged'] for r in sub)
                f1 = ' / '.join(f'{v[g]:,} ({v[g] / len(sub):.3f})' if has[g] else '-' for g in ('linear', 'd9', 'd46'))
                f2 = ' / '.join(f'{u[g]:,} ({u[g] / nt:.3f})' if has[g] else '-' for g in ('linear', 'd9', 'd46'))
                f3 = (f'{u["d9"] - u["linear"]:+,} ({(u["d9"] - u["linear"]) / len(sub):+.1%} of truth), '
                      f'{u["d46"] - u["linear"]:+,} ({(u["d46"] - u["linear"]) / len(sub):+.1%})') if has['linear'] else '-'
                L.append(f'| {s} | {p} | {kind} | {len(sub):,} ({nt:,}) | {f1} | {f2} | {f3} |')
                S[f'{s}/{p}/{kind}/tensor'] = dict(n=len(sub), not_tagged=nt, tensor=v, tensor_not_tagged=u)
    L.append('')

    # ---- Table 2: lost vs linear, per platform, per class, PoN ----
    L += ['## 2. Lost on a pangenome arm against linear (tensor on linear, none on the arm), per platform', '',
          'Each cell: alleles (PoN-tagged, share tagged; not tagged = graph-specific loss after the PoN).', '',
          '| sample | platform | kind | arm | lost (all) | absorbed (no candidate, ALT spelled by the graph) | no candidate, ALT not spelled | filtered | other allele | gained (tensor on arm, not on linear) |',
          '|---|---|---|---|---|---|---|---|---|---|']
    for s in ('HG008T', 'COLO829T'):
        for kind in ('SNV', 'INDEL'):
            for p in PLATS[s]:
                for g in ('d9', 'd46'):
                    sub = [r for r in rows if r['sample'] == s and r['kind2'] == kind]
                    lost = [r for r in sub if r[f'{p}_linear'] == 'T' and r[f'{p}_{g}'] != 'T']
                    by = defaultdict(list)
                    for r in lost:
                        by[cls(r[f'{p}_{g}'], r[f'{g}_spelled'])].append(r)
                    gained = [r for r in sub if r[f'{p}_linear'] != 'T' and r[f'{p}_{g}'] == 'T']
                    cells = [tagged_str(lost)[2]] + [tagged_str(by[c])[2] for c in ('absorbed', 'nocand_other', 'filtered', 'other_allele')] + [tagged_str(gained)[2]]
                    L.append(f'| {s} | {p} | {kind} | {g} | ' + ' | '.join(cells) + ' |')
                    key = f'{s}/{p}/{kind}/{g}'
                    S[key] = dict(lost=len(lost), lost_after_pon=sum(not r['tagged'] for r in lost),
                                  **{f'{c}': len(by[c]) for c in ('absorbed', 'nocand_other', 'filtered', 'other_allele')},
                                  **{f'{c}_after_pon': sum(not r['tagged'] for r in by[c]) for c in ('absorbed', 'nocand_other', 'filtered', 'other_allele')},
                                  gained=len(gained), gained_after_pon=sum(not r['tagged'] for r in gained))
    L.append('')

    # ---- Table 2b: d9 losses split by the read-level absorbed set (no s4 needed) ----
    L += ['## 2b. Lost on d9 against linear, split by the read-level absorbed set of the earlier analyses (ALT reads align with '
          'no edit to a d9 branch / skip edge on >= 1 platform; graph_absorbed_somatic_20261001 / _colo829t_20261005)', '',
          '| sample | platform | kind | lost on d9 | in the read-level absorbed set | other | of the read-level absorbed set: tensor on linear | of the read-level set: tensor on d9 now |',
          '|---|---|---|---|---|---|---|---|']
    for s in ('HG008T', 'COLO829T'):
        for kind in ('SNV', 'INDEL'):
            for p in PLATS[s]:
                sub = [r for r in rows if r['sample'] == s and r['kind2'] == kind]
                lost = [r for r in sub if r[f'{p}_linear'] == 'T' and r[f'{p}_d9'] != 'T']
                a = [r for r in lost if r['absorbed_d9_reads']]
                o = [r for r in lost if not r['absorbed_d9_reads']]
                rl = [r for r in sub if r['absorbed_d9_reads']]
                L.append(f'| {s} | {p} | {kind} | {tagged_str(lost)[2]} | {tagged_str(a)[2]} | {tagged_str(o)[2]} | '
                         f'{sum(r[f"{p}_linear"] == "T" for r in rl):,} of {len(rl):,} | {sum(r[f"{p}_d9"] == "T" for r in rl):,} |')
                S[f'{s}/{p}/{kind}/d9_readlevel'] = dict(lost=len(lost), lost_readlevel=len(a), lost_readlevel_after_pon=sum(not r['tagged'] for r in a),
                                                       lost_other=len(o), lost_other_after_pon=sum(not r['tagged'] for r in o))
    L.append('')

    # ---- Table 3: union over the three platforms ----
    L += ['## 3. Union over the three platforms (lost on >= 1 platform; absorbed on >= 1 platform), and lost on all three', '',
          '| sample | kind | arm | truth n (tagged) | lost on >= 1 platform | absorbed on >= 1 platform | lost on all 3 | absorbed on all 3 | lost on all 3, not absorbed anywhere | gained on all 3 (tensor on the arm, not on linear) |',
          '|---|---|---|---|---|---|---|---|---|---|']
    for s in ('HG008T', 'COLO829T'):
        for kind in ('SNV', 'INDEL'):
            sub = [r for r in rows if r['sample'] == s and r['kind2'] == kind]
            for g in ('d9', 'd46'):
                lost_any = [r for r in sub if any(r[f'{p}_linear'] == 'T' and r[f'{p}_{g}'] != 'T' for p in PLATS[s])]
                abs_any = [r for r in sub if any(r[f'{p}_linear'] == 'T' and r[f'{p}_{g}'] == 'N' and r[f'{g}_spelled'] for p in PLATS[s])]
                lost_all = [r for r in sub if all(r[f'{p}_linear'] == 'T' and r[f'{p}_{g}'] != 'T' for p in PLATS[s])]
                abs_all = [r for r in sub if all(r[f'{p}_linear'] == 'T' and r[f'{p}_{g}'] == 'N' for p in PLATS[s]) and r[f'{g}_spelled']]
                not_abs = [r for r in lost_all if not (r[f'{g}_spelled'] and any(r[f'{p}_{g}'] == 'N' for p in PLATS[s]))]
                gained_all = [r for r in sub if all(r[f'{p}_linear'] != 'T' and r[f'{p}_{g}'] == 'T' for p in PLATS[s])]
                L.append(f'| {s} | {kind} | {g} | {tagged_str(sub)[2]} | {tagged_str(lost_any)[2]} | {tagged_str(abs_any)[2]} | '
                         f'{tagged_str(lost_all)[2]} | {tagged_str(abs_all)[2]} | {tagged_str(not_abs)[2]} | {tagged_str(gained_all)[2]} |')
                S[f'{s}/union/{kind}/{g}'] = dict(n=len(sub), n_tagged=sum(r['tagged'] for r in sub),
                                                  lost_any=len(lost_any), lost_any_after_pon=sum(not r['tagged'] for r in lost_any),
                                                  absorbed_any=len(abs_any), absorbed_any_after_pon=sum(not r['tagged'] for r in abs_any),
                                                  lost_all=len(lost_all), lost_all_after_pon=sum(not r['tagged'] for r in lost_all),
                                                  absorbed_all=len(abs_all), absorbed_all_after_pon=sum(not r['tagged'] for r in abs_all))
    L.append('')

    # ---- Table 4: d46 vs d9 ----
    L += ['## 4. d46 against d9, per platform', '',
          'd46_only = tensor on d46, none on d9 (the d9 status class); d9_only = tensor on d9, none on d46 (the d46 status class).', '',
          '| sample | platform | kind | d46_only (all) | of which d9 absorbed | d9 no cand., not spelled | d9 filtered | d9_only (all) | of which d46 absorbed | d46 no cand., not spelled | d46 filtered | net d46 - d9 after PoN |',
          '|---|---|---|---|---|---|---|---|---|---|---|---:|']
    for s in ('HG008T', 'COLO829T', ):
        for kind in ('SNV', 'INDEL'):
            for p in PLATS[s]:
                sub = [r for r in rows if r['sample'] == s and r['kind2'] == kind]
                a = [r for r in sub if r[f'{p}_d46'] == 'T' and r[f'{p}_d9'] != 'T']
                b = [r for r in sub if r[f'{p}_d9'] == 'T' and r[f'{p}_d46'] != 'T']
                ac = defaultdict(list); bc = defaultdict(list)
                for r in a:
                    ac[cls(r[f'{p}_d9'], r['d9_spelled'])].append(r)
                for r in b:
                    bc[cls(r[f'{p}_d46'], r['d46_spelled'])].append(r)
                net = sum(not r['tagged'] for r in a) - sum(not r['tagged'] for r in b)
                L.append(f'| {s} | {p} | {kind} | {tagged_str(a)[2]} | {tagged_str(ac["absorbed"])[2]} | {tagged_str(ac["nocand_other"])[2]} | {tagged_str(ac["filtered"])[2]} | '
                         f'{tagged_str(b)[2]} | {tagged_str(bc["absorbed"])[2]} | {tagged_str(bc["nocand_other"])[2]} | {tagged_str(bc["filtered"])[2]} | {net:+,} |')
                S[f'{s}/{p}/{kind}/d46_vs_d9'] = dict(d46_only=len(a), d46_only_after_pon=sum(not r['tagged'] for r in a),
                                                     d46_only_d9_absorbed=len(ac['absorbed']), d9_only=len(b),
                                                     d9_only_after_pon=sum(not r['tagged'] for r in b),
                                                     d9_only_d46_absorbed=len(bc['absorbed']), net_after_pon=net)
    L.append('')

    # ---- Table 5: graph spelling of all truth alleles (independent of reads), and validation against the d9 read analysis ----
    L += ['## 5. Truth ALT spelled by a graph path at its locus (all truth alleles, independent of the reads)', '',
          '| sample | kind | n | spelled by d9 | spelled by d46 | d9 only | d46 only | both | spelled by d46, PoN-tagged |',
          '|---|---|---:|---|---|---:|---:|---:|---|']
    for s in ('HG008T', 'COLO829T', 'HCC1395T'):
        for kind in ('SNV', 'INDEL'):
            sub = [r for r in rows if r['sample'] == s and r['kind2'] == kind]
            n = len(sub)
            a = sum(r['d9_spelled'] for r in sub); b = sum(r['d46_spelled'] for r in sub)
            o9 = sum(r['d9_spelled'] and not r['d46_spelled'] for r in sub); o46 = sum(r['d46_spelled'] and not r['d9_spelled'] for r in sub)
            both = sum(r['d9_spelled'] and r['d46_spelled'] for r in sub)
            sp46 = [r for r in sub if r['d46_spelled']]
            L.append(f'| {s} | {kind} | {n:,} | {a:,} ({pct(a, n)}) | {b:,} ({pct(b, n)}) | {o9:,} | {o46:,} | {both:,} | {tagged_str(sp46)[2]} |')
            S[f'{s}/spelled/{kind}'] = dict(n=n, d9=a, d46=b, d9_only=o9, d46_only=o46, both=both,
                                            d46_tagged=sum(r['tagged'] for r in sp46))
    L.append('')
    L += ['Validation against the read-level d9 analysis (graph_absorbed_somatic_20261001 / _colo829t_20261005: ALT reads align '
          'with no edit to a d9 branch / skip edge, union of 3 platforms, 2026-09-30 tensor sets):', '',
          '| sample | kind | read-level absorbed (d9) | of which ALT spelled by d9 (this check) | of which no candidate on d9 on >= 1 platform now | this check: no candidate on d9 on >= 1 platform and spelled | of which in the read-level set |',
          '|---|---|---:|---:|---:|---:|---:|']
    for s in ('HG008T', 'COLO829T'):
        for kind in ('SNV', 'INDEL'):
            sub = [r for r in rows if r['sample'] == s and r['kind2'] == kind]
            ab = [r for r in sub if r['absorbed_d9_reads']]
            sp = sum(r['d9_spelled'] for r in ab)
            nc = sum(any(r[f'{p}_d9'] == 'N' for p in PLATS[s]) for r in ab)
            mine = [r for r in sub if r['d9_spelled'] and any(r[f'{p}_d9'] == 'N' for p in PLATS[s])]
            L.append(f'| {s} | {kind} | {len(ab):,} | {sp:,} | {nc:,} | {len(mine):,} | {sum(r["absorbed_d9_reads"] for r in mine):,} |')
            S[f'{s}/validation/{kind}'] = dict(read_level=len(ab), spelled=sp, no_cand_now=nc, mine=len(mine),
                                               mine_in_read_level=sum(r['absorbed_d9_reads'] for r in mine))
    L.append('')

    # ---- Table 6: INDEL absorbed by homopolymer context ----
    L += ['## 6. INDELs lost on the arm against linear, by GRCh38 homopolymer context (the PoN rule\'s own test: the ALT changes a homopolymer >= 7 bp)', '',
          '| sample | arm | class | HP>=7: n (tagged) | other INDEL: n (tagged) |', '|---|---|---|---|---|']
    for s in ('HG008T', 'COLO829T'):
        for g in ('d9', 'd46'):
            sub = [r for r in rows if r['sample'] == s and r['kind2'] == 'INDEL']
            lost_any = [r for r in sub if any(r[f'{p}_linear'] == 'T' and r[f'{p}_{g}'] != 'T' for p in PLATS[s])]
            abs_any = [r for r in lost_any if any(r[f'{p}_linear'] == 'T' and r[f'{p}_{g}'] == 'N' for p in PLATS[s]) and r[f'{g}_spelled']]
            for name, x in (('all truth', sub), ('lost on >= 1 platform', lost_any), ('absorbed on >= 1 platform', abs_any)):
                L.append(f'| {s} | {g} | {name} | {tagged_str([r for r in x if r["hp7"]])[2]} | {tagged_str([r for r in x if not r["hp7"]])[2]} |')
    L.append('')

    # ---- Table 7: HCC1395T on d46 ----
    L += ['## 7. HCC1395T (SEQC2 v1.2.1 truth), d46 only (no d9 / linear tensor sets)', '',
          '| platform | kind | n (tagged) | tensor | no candidate | no candidate, ALT spelled by d46 (absorbed) | no candidate, not spelled | filtered |',
          '|---|---|---|---|---|---|---|---|']
    for p in PLATS['HCC1395T']:
        for kind in ('SNV', 'INDEL'):
            sub = [r for r in rows if r['sample'] == 'HCC1395T' and r['kind2'] == kind]
            n = len(sub)
            T = [r for r in sub if r[f'{p}_d46'] == 'T']
            N = [r for r in sub if r[f'{p}_d46'] == 'N']
            Na = [r for r in N if r['d46_spelled']]
            No = [r for r in N if not r['d46_spelled']]
            F = [r for r in sub if r[f'{p}_d46'] == 'F']
            L.append(f'| {p} | {kind} | {tagged_str(sub)[2]} | {len(T):,} ({len(T) / n:.3f}) | {tagged_str(N)[2]} | {tagged_str(Na)[2]} ({pct(len(Na), n)} of truth) | {tagged_str(No)[2]} | {tagged_str(F)[2]} |')
            S[f'HCC1395T/{p}/{kind}/d46'] = dict(n=n, tensor=len(T), no_candidate=len(N), absorbed=len(Na),
                                                absorbed_after_pon=sum(not r['tagged'] for r in Na), filtered=len(F))
    # the same absorbed-share definition on HG008T / COLO829T d46, for comparison
    L += ['', 'Same definition (no candidate on d46 and ALT spelled by d46, share of all truth alleles of the kind) on every sample:', '',
          '| sample | platform | SNV absorbed (share) | INDEL absorbed (share) | INDEL absorbed after PoN (share of truth INDELs) |', '|---|---|---|---|---|']
    for s in ('HG008T', 'COLO829T', 'HCC1395T'):
        for p in PLATS[s]:
            cells = []
            for kind in ('SNV', 'INDEL'):
                sub = [r for r in rows if r['sample'] == s and r['kind2'] == kind]
                a = [r for r in sub if r[f'{p}_d46'] == 'N' and r['d46_spelled']]
                cells.append(f'{len(a):,} ({pct(len(a), len(sub))})')
                if kind == 'INDEL':
                    na = sum(not r['tagged'] for r in a)
                    cells.append(f'{na:,} ({pct(na, len(sub))})')
            L.append(f'| {s} | {p} | ' + ' | '.join(cells) + ' |')
    L.append('')

    # ---- Table 8: chr1 in the BED, INDEL: lost on d9 / d46 against linear, after PoN (the vcfeval scope) ----
    L += ['## 8. chr1 in the benchmark BED (the vcfeval test scope), INDEL: truth with a tensor and lost alleles after the PoN', '',
          '| sample | platform | truth INDEL (not tagged) | tensor & not tagged: linear / d9 / d46 | lost on d9 vs linear, not tagged (absorbed) | lost on d46 vs linear, not tagged (absorbed) |',
          '|---|---|---|---|---|---|']
    for s in ('HG008T', 'COLO829T'):
        for p in PLATS[s]:
            sub = [r for r in rows if r['sample'] == s and r['kind2'] == 'INDEL' and r['chrom'] == 'chr1' and r['in_bed'] == 'True']
            nt = [r for r in sub if not r['tagged']]
            tv = {g: sum(r[f'{p}_{g}'] == 'T' for r in nt) for g in ('linear', 'd9', 'd46')}
            cells = []
            for g in ('d9', 'd46'):
                lost = [r for r in nt if r[f'{p}_linear'] == 'T' and r[f'{p}_{g}'] != 'T']
                ab = [r for r in lost if r[f'{p}_{g}'] == 'N' and r[f'{g}_spelled']]
                cells.append(f'{len(lost):,} ({len(ab):,})')
                S[f'{s}/{p}/INDEL/chr1bed/{g}'] = dict(n=len(sub), not_tagged=len(nt), lost_after_pon=len(lost), absorbed_after_pon=len(ab),
                                                       tensor_not_tagged=tv)
            L.append(f'| {s} | {p} | {len(sub):,} ({len(nt):,}) | {tv["linear"]:,} / {tv["d9"]:,} / {tv["d46"]:,} | ' + ' | '.join(cells) + ' |')
    L.append('')
    open(f'{EXP}/tables.md', 'w').write('\n'.join(L) + '\n')
    json.dump(S, open(f'{EXP}/summary.json', 'w'), indent=1)
    print('\n'.join(L))


if __name__ == '__main__':
    main()
