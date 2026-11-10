"""s1: read-consensus contexts of every off-reference tensor of one (set, kind, chromosome), plus on-reference controls.

usage: s1_extract.py SET KIND CHROM OUT_DIR

Off-reference tensor = a record of tensors/<KIND>/<CHROM>_labels.ndjson with "grch38": null (the candidate node has no
GRCh38 position). Controls (on-reference, evenly spaced over the chromosome's records):
  ctrl_somatic   label 1, representative_allele_is_somatic_truth   (up to 40 per chromosome)
  ctrl_germline  label 2, representative_allele_is_germline_truth  (up to 40)
  ctrl_non       label 0, confident_no_truth_allele                (up to 20)

Per tensor, from the A1 rows (variant_summary row_groups): per column the plurality read base (channel 0, codes
A1 C2 G3 T4 N5 gap6) over the rows that cover it (channel 4 > 0), when >= 2 rows cover it (E9h build_queries.consensus;
ties -> lowest code). The site columns are [anchor, anchor + insertion_slots + span); every one needs a consensus and the
consensus bases there (gaps dropped) must spell the key's ALT site string (SNP: alt; INS: alt + reference; DEL:
reference[len(ref):]). The flanks extend from the site while columns have a consensus, up to 50 columns per side;
gap columns add no base. alt_q = left + ALT site + right, site interval [a, b) in alt_q (a == b for a pure deletion);
ref_q = left + reference + right (the same read haplotype with only the site set back to the graph REF).
QC: min_share (lowest plurality share over columns with >= 4 covering rows) and n_split (such columns with share < 0.7).
"""
import json
import os
import sys

import numpy as np

TENSORS = '/scratch/jshen/data/pansoma_v2_tensors'
BASE = np.array(list('?ACGTN'))
HALF = 50
MIN_ROWS = 2
CTRL = {'ctrl_somatic': (1, 'representative_allele_is_somatic_truth', 40),
        'ctrl_germline': (2, 'representative_allele_is_germline_truth', 40),
        'ctrl_non': (0, 'confident_no_truth_allele', 20)}
COLS = ['set', 'kind', 'chrom', 'candidate_id', 'group', 'label', 'reason', 'anchor', 'event_type', 'af', 'alt_count',
        'coverage', 'allele_count', 'n_a1', 'status', 'alt_q', 'site', 'ref_q', 'min_share', 'n_split']


def site_strings(r):
    lay = r['site_layout']
    ref_span = lay['reference']
    et = r['event_type']
    if et == 'SNP':
        return r['alt'], ref_span
    if et == 'INS':
        return r['alt'] + ref_span, ref_span
    if et == 'DEL':
        return ref_span[len(r['ref']):], ref_span
    raise ValueError(et)


def consensus(x, rows):
    """plurality code per column (1..6) or 0, plus share and covering-row count."""
    v = x[0, rows]
    cov = x[4, rows] > 0
    n = np.stack([((v == k) & cov).sum(0) for k in range(1, 7)])        # (6, width)
    tot = n.sum(0)
    best = n.argmax(0) + 1
    share = np.where(tot > 0, n.max(0) / np.maximum(tot, 1), 0.0)
    cons = np.where(tot >= MIN_ROWS, best, 0)
    return cons, share, tot


def contexts(x, r):
    a1 = next((g for g in r.get('row_groups', []) if g['allele'] == 'A1'), None)
    if a1 is None or a1['end_row'] - a1['start_row'] < 1:
        return 'no_a1', None
    cons, share, tot = consensus(x, slice(a1['start_row'], a1['end_row']))
    lay = r['site_layout']
    s0 = r.get('anchor_column', 50)
    s1 = s0 + lay['insertion_slots'] + lay['span']
    if (cons[s0:s1] == 0).any():
        return 'no_site_consensus', None
    alt_site, ref_site = site_strings(r)
    site_read = ''.join(BASE[c] for c in cons[s0:s1] if 1 <= c <= 5)
    if site_read != alt_site:
        return 'site_consensus_not_alt', None
    width = x.shape[2]
    lo = s0
    while lo > max(0, s0 - HALF) and cons[lo - 1] != 0:
        lo -= 1
    hi = s1                                                               # exclusive
    while hi < min(width, s1 + HALF) and cons[hi] != 0:
        hi += 1
    left = ''.join(BASE[c] for c in cons[lo:s0] if 1 <= c <= 5)
    right = ''.join(BASE[c] for c in cons[s1:hi] if 1 <= c <= 5)
    span = slice(lo, hi)
    deep = tot[span] >= 4
    min_share = float(share[span][deep].min()) if deep.any() else 1.0
    n_split = int((share[span][deep] < 0.7).sum())
    a = len(left)
    return 'ok', dict(alt_q=left + alt_site + right, site=f'{a}-{a + len(alt_site)}', ref_q=left + ref_site + right,
                      min_share=f'{min_share:.3f}', n_split=n_split)


def runs(idx, gap=4, maxlen=64):
    """coalesce sorted tensor indices into (lo, hi) read runs."""
    out = []
    for i in idx:
        if out and i - out[-1][1] <= gap and i + 1 - out[-1][0] <= maxlen:
            out[-1][1] = i + 1
        else:
            out.append([i, i + 1])
    return out


def main():
    set_name, kind, chrom, out_dir = sys.argv[1:5]
    tdir = f'{TENSORS}/{set_name}/tensors/{kind}'
    lab_path, sum_path = f'{tdir}/{chrom}_labels.ndjson', f'{tdir}/{chrom}_variant_summary.ndjson'
    # pass 1: control candidates by line number
    ctrl_lines = {g: [] for g in CTRL}
    n_lines = 0
    with open(lab_path) as f:
        for i, line in enumerate(f):
            n_lines += 1
            if '"grch38": null' in line:
                continue
            for g, (lab, reason, _) in CTRL.items():
                if f'"label": {lab},' in line and f'"reason": "{reason}"' in line:
                    ctrl_lines[g].append(i)
    want = {}
    for g, lines in ctrl_lines.items():
        k = CTRL[g][2]
        if len(lines) > k:
            lines = [lines[int(j * len(lines) / k)] for j in range(k)]
        for i in lines:
            want[i] = g
    # pass 2: records
    recs = []
    with open(lab_path) as fl, open(sum_path) as fs:
        for i, (ll, sl) in enumerate(zip(fl, fs)):
            off = '"grch38": null' in ll
            if not off and i not in want:
                continue
            lab = json.loads(ll)
            s = json.loads(sl)
            if lab['candidate_id'] != s['candidate_id']:
                raise SystemExit(f'line {i}: labels / summary out of step ({lab["candidate_id"]} / {s["candidate_id"]})')
            recs.append((('offref' if off else want[i]), lab, s))
    # tensors, shard by shard, in coalesced runs
    by_shard = {}
    for j, (_, lab, s) in enumerate(recs):
        by_shard.setdefault(s['shard_file'], []).append((s['index_within_shard'], j))
    rows_out = [None] * len(recs)
    for shard, items in sorted(by_shard.items()):
        items.sort()
        mm = np.load(f'{tdir}/{shard}', mmap_mode='r')
        row_bytes = int(np.prod(mm.shape[1:])) * mm.dtype.itemsize
        pos = {i: j for i, j in items}
        with open(f'{tdir}/{shard}', 'rb') as fh:
            for lo, hi in runs([i for i, _ in items]):
                fh.seek(mm.offset + lo * row_bytes)
                buf = np.frombuffer(fh.read((hi - lo) * row_bytes), dtype=mm.dtype).reshape((hi - lo,) + mm.shape[1:])
                for i in range(lo, hi):
                    j = pos.get(i)
                    if j is None:
                        continue
                    group, lab, s = recs[j]
                    status, ctx = contexts(buf[i - lo], s)
                    an = lab.get('anchor')
                    g38 = lab.get('grch38')
                    if g38:   # on-reference controls: the GRCh38 interval of the candidate
                        an = dict(chrom=g38['chrom'], start=g38['pos0'], end=g38['pos0'] + max(1, len(g38['ref'])))
                    a1 = next((g for g in s.get('row_groups', []) if g['allele'] == 'A1'), None)
                    rows_out[j] = dict(
                        set=set_name, kind=kind, chrom=chrom, candidate_id=lab['candidate_id'], group=group,
                        label=lab['label'], reason=lab['reason'],
                        anchor=(f'{an["chrom"]}:{an["start"]}-{an["end"]}' if an else ''),
                        event_type=s['event_type'], af=f'{s["af"]:.4f}', alt_count=s['alt_count'],
                        coverage=s['coverage'], allele_count=s.get('allele_count', ''),
                        n_a1=(a1['end_row'] - a1['start_row'] if a1 else 0), status=status,
                        alt_q='', site='', ref_q='', min_share='', n_split='')
                    if ctx:
                        rows_out[j].update(ctx)
        del mm
    os.makedirs(out_dir, exist_ok=True)
    out = f'{out_dir}/{set_name}.{kind}.{chrom}.tsv'
    with open(out + '.tmp', 'w') as w:
        w.write('\t'.join(COLS) + '\n')
        for r in rows_out:
            w.write('\t'.join(str(r[c]) for c in COLS) + '\n')
    os.replace(out + '.tmp', out)
    n = {}
    for r in rows_out:
        n[(r['group'], r['status'])] = n.get((r['group'], r['status']), 0) + 1
    print(set_name, kind, chrom, 'labels lines', n_lines, 'records', len(rows_out),
          ' '.join(f'{g}/{s}={v}' for (g, s), v in sorted(n.items())), flush=True)


if __name__ == '__main__':
    main()
