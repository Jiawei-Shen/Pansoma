"""s6: training labels for the off-reference tensors of one tumor, written as new tensor sets <SET>_offref.

usage: s6_dataset.py TUMOR [--dry-run]

Evidence per off-reference key (s3 rows of every platform set, s4 read counts, s5 tiers):
  somatic   s5 tier T1 -> label 1, reason offref_somatic_t1; T2 -> 1, offref_somatic_t2 (tumor_only in the tumor
            assembly with the REF or only a site difference in the normal assembly, no ALT in the matched-normal reads,
            tumor read support at the locus; T1 also on >= 2 platforms and clean)
  germline  the A1 context is in a normal haplotype at <= 2 loci (asm_class normal, n_alt_loci <= 2: one per haplotype)
            on the tensor's chromosome (HG008-N contigs chrN_hapK; COLO829BL DSA contigs by the SMaHT GRCh38ToDSA
            chain; a hit on another chromosome or an unplaced contig is a paralog / collapsed copy and does not count),
            and the Illumina normal reads carry the ALT at the locus (ALT >= 3 and fraction >= 0.15; the germline rule
            of s5's repeat contexts) -> 2, offref_germline
  non       every usable row is asm_class neither (the ALT context is in neither assembly), the key has a tensor in one
            platform set only and tumor read support (where counted) on <= 1 platform -> 0, offref_non; the same with
            ALT in the normal reads (s5 rejected_normal: a systematic error in tumor and normal) -> 0, offref_non_normal_alt
  else      no label (uncertain: reads-only, tumor_only without confirmation, normal-assembly hit without read support,
            multi-copy, ALT context in neither assembly but tensors on several platforms, ...)
Only tensors labelled -1 with reason off_reference_no_truth_match (anchored inside the confident region, no truth
allele) and not on the test chromosome chr1 (unless --relabel-test: training never reads chr1, and keeping it as in the
source keeps the per-tensor chr1 test comparable) are relabelled, the region where the existing labels call untruthed tensors non at evaluation; outside the
region and without an anchor the labels stay -1, as for every other tensor there. Uncertain keys keep their label (-1)
and reason, so evaluation still counts them as non (data.EVAL_AS_NON), as before. Every off-reference record gets the
field "offref" (proposed label and class, s5 tier, platform sets) and relabelled records also "reason_original".

Output per target set (S6_SETS) and kind (the kinds the source set has; COLO829T *_d46_snvaf05 keep SNV only): TENSORS/<SET>_offref/tensors/<KIND>/ with the shards' data, the summaries
and manifest.json linked to the source set, new <chrom>_shard_NNNNN_labels.npy, <chrom>_labels.ndjson and
labels.manifest.json (format, region BEDs and the rest of the source's; rules_sha256 of this file, created now,
counts per chromosome, totals and reasons recounted, an offref block); the set's recall tables, non_autosomal/ and
shared/ linked. Site features (train --site-features) do not depend on labels: <SET>_offref.<KIND>.npz is linked to
the source set's file in SF_DIRS (the snvaf05 sets' INDEL to the d46 INDEL file: same tensors). The new sets are then
indexed with pansoma_net_v2.data.KindIndex and their site features loaded as a check. DATA/s6/<TUMOR>/key_labels.tsv.gz
holds the key evidence and controls.md the same rules applied to the on-reference controls (against their truth
labels), with each control key's platforms taken from the tensors of the target sets, as for off-reference keys.
"""
import csv
import datetime
import glob
import gzip
import hashlib
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, platform  # noqa: E402

REPO = '/scratch/jshen/Github/Pansoma'
TENSORS = '/scratch/jshen/data/pansoma_v2_tensors'
S6_SETS = {'HG008T': ['HG008T_Illumina_d46', 'HG008T_PacBio_d46', 'HG008T_ONT_d46'],
           'COLO829T': ['COLO829T_Illumina_d46', 'COLO829T_fiberseq_d46', 'COLO829T_ONT_d46',
                        'COLO829T_Illumina_d46_snvaf05', 'COLO829T_fiberseq_d46_snvaf05']}
RELABEL_REASON = 'off_reference_no_truth_match'
LABEL_NAME = {-1: 'ignore', 0: 'non', 1: 'somatic', 2: 'germline'}
GERM_ALT, GERM_FRAC, ABS_DEPTH, ABS_FRAC = 3, 0.15, 10, 0.10
TEST_CHROMS = ('chr1',)
SF_DIRS = {'SNV': '/scratch/jshen/data/snv_af05_20261103/site_features',
           'INDEL': '/scratch/jshen/data/indel_gap_20261020/site_features'}
DSA_CHAIN = '/scratch/qfu/COLO829BL_DSA/SMHTCOLO829BL-X-X-M45-XX-uwsc-SMASF9YNBDYM-GRCh38ToDSA.chain'
_CID = __import__('re').compile(r'"candidate_id": "([^"]+)"')
_RSN = __import__('re').compile(r'"reason": "([^"]+)"')


def contig_chroms(tumor):
    """normal-assembly contig -> GRCh38 chromosome ('' when unplaced or unknown)."""
    if tumor == 'HG008T':
        return None            # by name: chrN_hapK
    span = defaultdict(Counter)
    for line in open(DSA_CHAIN):
        if line.startswith('chain'):
            f = line.split()
            span[f[2]][f[7]] += int(f[6]) - int(f[5])
    return {c: max(v, key=v.get) for c, v in span.items()}


def normal_chrom(normal_aln, table):
    if not normal_aln:
        return ''
    contig = normal_aln.split(':')[1]
    if table is None:
        return contig.rsplit('_hap', 1)[0] if contig.startswith('chr') and '_hap' in contig else ''
    return table.get(contig, '')


def n_illumina_counts(tumor):
    d = {}
    for p in glob.glob(f'{DATA}/s4/{tumor}/counts/N_Illumina.*.tsv'):
        for r in csv.DictReader(open(p), delimiter='\t'):
            d[r['locus_id']] = (int(r['alt']), int(r['ref']))
    return d


def germ_call(c):
    if c is None:
        return 'no_count'
    alt, ref = c
    frac = alt / (alt + ref) if alt + ref else 0.0
    if alt >= GERM_ALT and frac >= GERM_FRAC:
        return 'present'
    if alt + ref >= ABS_DEPTH and frac < ABS_FRAC:
        return 'absent'
    return 'low_depth'


def key_label(rows, tier, s5, plats=None):
    """(label, reason, note) of one key from its rows (s3 + locus + germline call), its s5 tier and s5 record; plats =
    the platforms with a tensor of the key (default: those of its rows)."""
    usable = [r for r in rows if r['asm_class'] not in ('no_query', 'short_flank')]
    classes = {r['asm_class'] for r in usable}
    if tier == 'T1':
        return 1, 'offref_somatic_t1', ''
    if tier == 'T2':
        return 1, 'offref_somatic_t2', ''
    one = [r for r in usable if r['asm_class'] == 'normal' and int(r['n_alt_loci'] or 0) <= 2]
    germ_rows = [r for r in one if r['nchrom'] and r['nchrom'] == r['chrom']]
    calls = {r['germ'] for r in germ_rows}
    if 'present' in calls:
        return 2, 'offref_germline', ''
    if not usable:
        return -1, '', 'no_query'
    if classes == {'neither'}:
        sets = plats if plats is not None else {platform(r['set']) for r in rows}
        support = {x for x in (s5.get('support') or '').split(',') if x} if s5 else set()
        if tier == 'T3_reads_only':
            return -1, '', 'reads_only'
        if len(sets) > 1:
            return -1, '', 'neither_multiplatform'
        if len(support) > 1:
            return -1, '', 'neither_supported'
        if tier == 'rejected_normal':
            return 0, 'offref_non_normal_alt', ''
        return 0, 'offref_non', ''
    if 'tumor_only' in classes:
        return -1, '', 'tumor_only_normal_alt' if tier == 'rejected_normal' else 'tumor_only_unconfirmed'
    if 'normal' in classes:
        if not one:
            return -1, '', 'asm_normal_multicopy'
        if not germ_rows:
            return -1, '', 'asm_normal_other_chrom'
        if 'absent' in calls:
            return -1, '', 'asm_normal_reads_absent'
        return -1, '', 'asm_normal_no_read_check'
    return -1, '', 'mixed'


def evidence(tumor, ctrl_plats=None):
    """key -> dict(label, reason, note, tier, sets) for the off-reference keys, plus the control keys' (group) labels
    (ctrl_plats: control key -> platforms with a tensor of it)."""
    nill = n_illumina_counts(tumor)
    table = contig_chroms(tumor)
    rl = {}
    for r in csv.DictReader(gzip.open(f'{DATA}/s4/{tumor}/row_loci.tsv.gz', 'rt'), delimiter='\t'):
        rl[(r['set'], r['kind'], r['candidate_id'])] = r['locus']
    s5 = {}
    for r in csv.DictReader(gzip.open(f'{DATA}/s5/{tumor}/keys.tsv.gz', 'rt'), delimiter='\t'):
        s5[(r['group'], r['key'])] = r
    rows = defaultdict(list)
    keep = ('set', 'kind', 'chrom', 'candidate_id', 'group', 'label', 'asm_class', 'n_alt_loci')
    for r in csv.DictReader(gzip.open(f'{DATA}/s3/{tumor}/rows.tsv.gz', 'rt'), delimiter='\t'):
        x = {k: r[k] for k in keep}
        lid = rl.get((r['set'], r['kind'], r['candidate_id']), '')
        x['germ'] = germ_call(nill.get(lid)) if lid.startswith('L') else 'no_locus'
        x['nchrom'] = normal_chrom(r['normal_aln'], table) if r['asm_class'] == 'normal' else ''
        rows[(r['group'], r['candidate_id'])].append(x)
    out = {}
    for (group, key), rs in rows.items():
        rec = s5.get((group, key))
        tier = rec['tier'] if rec else ''
        plats = ctrl_plats.get(key) if (ctrl_plats is not None and group != 'offref') else None
        label, reason, note = key_label(rs, tier, rec, plats)
        out[(group, key)] = dict(label=label, reason=reason, note=note, tier=tier,
                                 sets=','.join(sorted({platform(r['set']) for r in rs})),
                                 truth=rs[0]['label'], kind=rs[0]['kind'])
    return out


def controls_report(tumor, ev, path):
    c = Counter()
    for (group, key), e in ev.items():
        if group == 'offref':
            continue
        c[(group, e['kind'], e['reason'] or 'uncertain:' + e['note'])] += 1
    with open(path, 'w') as w:
        w.write(f'# {tumor}: s6 rules applied to the on-reference controls (keys)\n\n')
        w.write('| group (truth) | kind | keys | offref_somatic_t1 | offref_somatic_t2 | offref_germline | offref_non | '
                'offref_non_normal_alt | uncertain |\n|---|---|---|---|---|---|---|---|---|\n')
        for g in ('ctrl_somatic', 'ctrl_germline', 'ctrl_non'):
            for k in ('SNV', 'INDEL'):
                n = sum(v for kk, v in c.items() if kk[:2] == (g, k))
                if not n:
                    continue
                row = [c[(g, k, r)] for r in ('offref_somatic_t1', 'offref_somatic_t2', 'offref_germline', 'offref_non',
                                             'offref_non_normal_alt')]
                unc = n - sum(row)
                w.write(f'| {g} | {k} | {n} | ' + ' | '.join(f'{x} ({100 * x / n:.1f} %)' for x in row + [unc]) + ' |\n')
        w.write('\nControl keys get the platforms of every target set with a tensor of them, as off-reference keys do. '
                'The snvaf05-only SNV keys (AF 0.05-0.08, COLO829T) have no controls of their own.\n')
        w.write('\nUncertain, by note:\n\n')
        for g in ('ctrl_somatic', 'ctrl_germline', 'ctrl_non'):
            notes = Counter({kk[2]: v for kk, v in c.items() if kk[0] == g and kk[2].startswith('uncertain')})
            w.write(f'- {g}: ' + ', '.join(f'{k.split(":", 1)[1]} {v}' for k, v in notes.most_common()) + '\n')


def sf_source(set_name, kind):
    """the source set's site-feature file (resolved) for <set>_offref, or None."""
    base = set_name.replace('_snvaf05', '') if kind == 'INDEL' else set_name   # snvaf05 INDEL tensors = d46 INDEL
    p = Path(SF_DIRS[kind]) / f'{base}.{kind}.npz'
    return os.path.realpath(p) if p.exists() else None


def kinds_of(set_name):
    """the kinds the source set has (the snvaf05 sets keep SNV only; their INDEL was the d46 INDEL and was deleted)."""
    return [k for k in ('SNV', 'INDEL') if (Path(TENSORS) / set_name / 'tensors' / k / 'manifest.json').exists()]


def link_site_features(set_name):
    out = {}
    for kind in kinds_of(set_name):
        src = sf_source(set_name, kind)
        dst = Path(SF_DIRS[kind]) / f'{set_name}_offref.{kind}.npz'
        if src is None:
            out[kind] = 'missing'
            continue
        link(src, str(dst))
        out[kind] = str(dst)
    return out


def link(src, dst):
    if os.path.lexists(dst):
        if os.path.islink(dst) and os.readlink(dst) == src:
            return
        raise SystemExit(f'{dst} exists and is not a link to {src}')
    os.symlink(src, dst)


def write_set(set_name, ev, rules_sha, dry, relabel_test=False, ctrl_keys=None, ctrl_plats=None):
    src_root = Path(TENSORS) / set_name
    dst_root = Path(TENSORS) / f'{set_name}_offref'
    summary = Counter()
    for kind in kinds_of(set_name):
        src = src_root / 'tensors' / kind
        dst = dst_root / 'tensors' / kind
        if not dry:
            dst.mkdir(parents=True, exist_ok=True)
        manifest = json.loads((src / 'manifest.json').read_text())
        lman = json.loads((src / 'labels.manifest.json').read_text())
        counts = {}
        reasons = Counter()
        for chrom, info in manifest['chromosomes'].items():
            shards = {s['file']: np.load(src / s['file'].replace('_data.npy', '_labels.npy')).copy()
                      for s in info['shards']}
            cc = Counter()
            tmp = dst / f'{chrom}_labels.ndjson.tmp'
            with open(src / f'{chrom}_labels.ndjson') as fin, (open(tmp, 'w') if not dry else open(os.devnull, 'w')) as fo:
                for line in fin:
                    if '"grch38": null' not in line:
                        fo.write(line)
                        reasons[_RSN.search(line).group(1)] += 1
                        if ctrl_keys:
                            cid = _CID.search(line).group(1)
                            if cid in ctrl_keys:
                                ctrl_plats.setdefault(cid, set()).add(platform(set_name))
                        continue
                    r = json.loads(line)
                    e = ev.get(('offref', r['candidate_id']))
                    if e is None:
                        summary[(kind, r['reason'], 'not_analysed')] += 1
                        r['offref'] = dict(label=None, **{'class': 'not_analysed'})
                    else:
                        r['offref'] = {'label': e['label'], 'class': e['reason'] or 'uncertain:' + e['note'],
                                       'tier': e['tier'], 'sets': e['sets']}
                        if (r['label'] == -1 and r['reason'] == RELABEL_REASON and e['label'] != -1 and
                                (relabel_test or chrom not in TEST_CHROMS)):
                            if shards[r['shard_file']][r['index_within_shard']] != -1:
                                raise SystemExit(f'{set_name} {kind} {r["candidate_id"]}: labels.npy and ndjson disagree')
                            shards[r['shard_file']][r['index_within_shard']] = e['label']
                            r['reason_original'] = r['reason']
                            r['label'], r['label_name'], r['reason'] = e['label'], LABEL_NAME[e['label']], e['reason']
                            summary[(kind, RELABEL_REASON, e['reason'])] += 1
                        else:
                            summary[(kind, r['reason'], r['offref']['class'])] += 1
                    reasons[r['reason']] += 1
                    fo.write(json.dumps(r) + '\n')
            # per-chromosome label counts from the arrays
            for f, a in shards.items():
                for v, n in zip(*np.unique(a, return_counts=True)):
                    cc[LABEL_NAME[int(v)]] += int(n)
            counts[chrom] = {k: v for k, v in cc.items() if k}
            if not dry:
                os.replace(tmp, dst / f'{chrom}_labels.ndjson')
                for f, a in shards.items():
                    out = dst / f.replace('_data.npy', '_labels.npy')
                    np.save(str(out) + '.tmp.npy', a)
                    os.replace(str(out) + '.tmp.npy', out)
                    link(str(src / f), str(dst / f))
                link(str(src / info['summary']), str(dst / info['summary']))
        if not dry:
            for p in src.iterdir():   # other kind files (manifest, filtered / unsupported candidates, reports)
                if p.name.endswith(('_data.npy', '_labels.npy', '_labels.ndjson', '_variant_summary.ndjson')) or \
                        p.name == 'labels.manifest.json':
                    continue
                link(str(p), str(dst / p.name))
            totals = Counter()
            for v in counts.values():
                totals.update(v)
            sf_src = sf_source(set_name, kind)
            lman.update(rules_sha256=rules_sha, created=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                        counts=counts, totals=dict(totals), reasons=dict(reasons),
                        source_labels=str(src / 'labels.manifest.json'),
                        offref=dict(script=f'{REPO}/analysis/offref_somatic_candidates_20261109/s6_dataset.py',
                                    relabelled_reason=RELABEL_REASON, test_chroms_relabelled=relabel_test,
                                    test_chroms=list(TEST_CHROMS), site_features=sf_src,
                                    source_totals=lman.get('totals'), source_reasons=lman.get('reasons'),
                                    reasons=['offref_somatic_t1', 'offref_somatic_t2', 'offref_germline', 'offref_non',
                                             'offref_non_normal_alt']))
            tmpm = dst / 'labels.manifest.json.tmp'
            tmpm.write_text(json.dumps(lman, indent=2) + '\n')
            os.replace(tmpm, dst / 'labels.manifest.json')
    if not dry:
        for p in (src_root / 'tensors').iterdir():
            if p.name not in ('SNV', 'INDEL'):
                link(str(p), str(dst_root / 'tensors' / p.name))
        for p in src_root.iterdir():
            if p.name not in ('tensors', 'README.txt'):
                link(str(p), str(dst_root / p.name))
        (dst_root / 'README.txt').write_text(
            f'{set_name}_offref: {set_name} with training labels for its off-reference tensors (anchored in the confident '
            f'region, no truth allele: reason {RELABEL_REASON}) from tumor / normal assemblies and matched-normal reads.\n'
            f'New reasons: offref_somatic_t1, offref_somatic_t2 (1), offref_germline (2), offref_non, '
            f'offref_non_normal_alt (0); every other tensor as in {set_name}. Shards, summaries and manifests are links '
            f'to {src_root}.\nRules and evidence: {REPO}/analysis/offref_somatic_candidates_20261109 (s6_dataset.py, '
            f'README.md).\n')
    return summary


def check_index(set_name, cache_dir):
    sys.path.insert(0, f'{REPO}/machine_learning')
    from pansoma_net_v2.data import KindIndex
    from pansoma_net_v2 import site_features
    out = {}
    for kind in kinds_of(set_name):
        index = KindIndex(Path(TENSORS) / f'{set_name}_offref' / 'tensors' / kind, cache_dir=cache_dir, rebuild=True)
        arrays, meta = index.arrays, index.meta
        try:
            sf = site_features.load(SF_DIRS[kind], index)
            sf_ok = f"ok {sf['path']}"
        except ValueError as e:
            sf_ok = f'FAILED {e}'
        reasons = meta['reasons']
        lab = arrays['label']
        out[kind] = dict(tensors=int(len(lab)), labels={LABEL_NAME[int(v)]: int((lab == v).sum()) for v in (-1, 0, 1, 2)},
                         offref_labelled={r: int(((arrays['reason'] == reasons.index(r)) & (lab != -1)).sum())
                                          for r in reasons if r.startswith('offref_')},
                         eval_labelled_offref=int(((arrays['eval_label'] != -1) & arrays['off_reference']).sum()),
                         site_features=sf_ok)
    return out


def main():
    tumor = sys.argv[1]
    dry = '--dry-run' in sys.argv
    relabel_test = '--relabel-test' in sys.argv
    rules_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    ev = evidence(tumor)
    od = Path(f'{DATA}/s6/{tumor}')
    od.mkdir(parents=True, exist_ok=True)
    ctrl_keys = {k for (g, k) in ev if g != 'offref'}
    ctrl_plats = {}
    report = {}
    for s in S6_SETS[tumor]:
        summ = write_set(s, ev, rules_sha, dry, relabel_test, ctrl_keys, ctrl_plats)
        report[s] = {f'{k} | {a} | {b}': v for (k, a, b), v in sorted(summ.items())}
        if not dry:
            report[s + ' site_features'] = link_site_features(s)
        print(s, 'written' if not dry else 'dry run', flush=True)
    # controls with the platforms of their tensors, as the off-reference keys get them
    ev_c = evidence(tumor, ctrl_plats)
    for k, e in ev_c.items():
        if k[0] != 'offref':
            ev[k] = e
    with gzip.open(od / 'key_labels.tsv.gz.tmp', 'wt') as w:
        w.write('group\tkey\tkind\ttruth_label\tlabel\treason\tnote\ttier\tsets\n')
        for (g, k), e in sorted(ev.items()):
            w.write(f"{g}\t{k}\t{e['kind']}\t{e['truth']}\t{e['label']}\t{e['reason']}\t{e['note']}\t{e['tier']}\t{e['sets']}\n")
    os.replace(od / 'key_labels.tsv.gz.tmp', od / 'key_labels.tsv.gz')
    controls_report(tumor, ev, od / 'controls.md')
    kc = Counter((e['kind'], e['reason'] or 'uncertain:' + e['note']) for (g, _), e in ev.items() if g == 'offref')
    print(tumor, 'off-reference keys by proposed label:', dict(sorted(kc.items())), flush=True)
    if not dry:
        for s in S6_SETS[tumor]:
            report[s + ' index'] = check_index(s, od / 'index_cache')
            print(s, 'index', json.dumps(report[s + ' index']), flush=True)
    (od / 'datasets.json').write_text(json.dumps(report, indent=1) + '\n')
    print('done ->', od)


if __name__ == '__main__':
    main()
