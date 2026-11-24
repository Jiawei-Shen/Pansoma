"""E9c per run: contexts, bubble locus and graph-native population tags of the unplaced PASS calls of one vcf dir,
plus controls (placed PASS calls: rtg TP / FP after the PoN, and PASS calls the GRCh38 PoN tags) and the
assembly-verified off-reference somatic keys (gold) present in the run's tensor set.

usage: e9c_run.py VCFDIR NAME [--controls N] [--force]
Output: DATA/runs/NAME/calls.tsv (one row per (group, candidate)) and DATA/runs/NAME/run.json.

Graph-native population tags (each for the graph-path context `g` and the read-consensus context `r`):
  ingraph   a walk of the run's own graph (d9 / d46 edges, node sequences) within the call's neighbourhood spells the
            ALT context (flank f each side of the event): the allele is a graph allele (>= 10 % of HPRC haplotypes
            for d9 / d46) that the reads did not take
  hprc      a haplotype of the full HPRC v1.1 graph (raw deconstruct VCF: 44 samples + CHM13, every allele, no
            frequency filter) or GRCh38 itself spells it inside the call's GRCh38 bubble interval +- 100 bp
Flank f: fixed 5 / 10 / 15 (11 / 21 / 31-mers for an SNV) and adaptive (10, widened by 5 up to 40 while the ALT query
is a substring of the REF path or REF read context, i.e. inside a repeat).
"""
import argparse
import json
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e9c_lib as L  # noqa: E402

FLANKS = (5, 10, 15)


def first_id(field):
    return field.split(';')[0]


def calls_from_vcfdir(vcfdir):
    unplaced = [p for p in os.listdir(vcfdir) if p.endswith('.linear.unplaced.vcf.gz')]
    if len(unplaced) != 1:
        raise SystemExit(f'{vcfdir}: expected one *.linear.unplaced.vcf.gz, found {unplaced}')
    return f'{vcfdir}/{unplaced[0]}'


def process_row(row, rec, ts, g, hprc):
    """Contexts, locus, in-graph and HPRC tags of one row (in place)."""
    ctx = L.tensor_contexts(ts.tensor(rec), rec)
    for k in ('alt_read', 'alt_graph', 'ref_graph', 'ref_read'):
        v = ctx.get(k) if ctx else None
        row[k] = v[0] if v else ''
        row[f'{k}_site'] = f'{v[1]}-{v[2]}' if v else ''
    if not ctx or not ctx.get('alt_graph'):
        return
    s, a0, b0 = ctx['ref_graph'] if ctx.get('ref_graph') else ctx['alt_graph']
    row['homopolymer'] = L.homopolymer_at(s, a0, b0)
    unit, copies = L.str_at(s, a0, b0)
    row['str_unit'], row['str_copies'] = unit, copies
    node = row['node']
    loc = g.locus(node, row['chrom'])
    row.update({f'locus_{k}': v for k, v in loc.items()})
    row['node_path_count'] = g.path_count(node)
    refs = [ctx['ref_graph'][0] if ctx.get('ref_graph') else '', ctx['ref_read'][0] if ctx.get('ref_read') else '']
    queries = {}
    for src in ('alt_graph', 'alt_read'):
        if not ctx.get(src):
            continue
        for f in FLANKS:
            q, lf, rf = L.query_at(ctx[src], f)
            if lf >= f and rf >= f:
                queries[f'{src[4]}{f}'] = q
        q, f, amb = L.adaptive_query(ctx[src], refs)
        queries[f'{src[4]}A'] = q
        row[f'{src[4]}A_flank'] = f
        row[f'{src[4]}A_repeat_ambiguous'] = amb
    # in-graph walk
    reach = max(loc['left_bp'] or 0, loc['right_bp'] or 0)
    longest = max((len(q) for q in queries.values()), default=0)
    hs = g.neighborhood(node, min(5000, reach + longest + 50))
    row['nbhd_nodes'] = len(hs) // 2
    for qk, q in queries.items():
        walk = g.spell(q, hs)
        row[f'ingraph_{qk}'] = 'capped' if walk == 'capped' else int(bool(walk))
        if walk and walk != 'capped' and qk in ('g10', 'gA'):
            row[f'ingraph_{qk}_walk_offref'] = int(any(not g.on_ref(w) for w in walk))
            row[f'ingraph_{qk}_walk_has_node'] = int(node in walk)
    # HPRC v1.1 haplotypes at the GRCh38 locus
    if loc['status'] in ('bubble', 'one_side', 'on_ref') and loc['lo'] is not None:
        haps, lo, hi, nsites, skipped = hprc.window(loc['chrom'], loc['lo'] - 100, loc['hi'] + 100)
        row['hprc_window'] = f'{loc["chrom"]}:{lo}-{hi}'
        row['hprc_sites'] = nsites
        row['hprc_haps'] = len(haps) - 1
        for qk, q in queries.items():
            c = L.HprcHaplotypes.carriers(haps, q)
            row[f'hprc_{qk}'] = len([h for h in c if h != 'GRCh38'])
            row[f'hprc_{qk}_grch38'] = int('GRCh38' in c)


def write_rows(rows, path):
    cols = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(f'{path}.tmp', 'w') as w:
        w.write('\t'.join(cols) + '\n')
        for r in rows:
            w.write('\t'.join(str(r.get(k, '')) for k in cols) + '\n')
    os.replace(f'{path}.tmp', path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('vcfdir')
    ap.add_argument('name')
    ap.add_argument('--controls', type=int, default=150)
    ap.add_argument('--max-unplaced', type=int, default=3000, help='random sample (seed 0) of the unplaced PASS calls above this')
    ap.add_argument('--force', action='store_true')
    a = ap.parse_args()
    out = f'{L.DATA}/runs/{a.name}'
    if os.path.exists(f'{out}/calls.tsv') and not a.force:
        print(f'{out}/calls.tsv exists; skip')
        return
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    upath = calls_from_vcfdir(a.vcfdir)
    hdr = L.vcf_header(upath)
    graph_name = L.graph_of(hdr['pansoma_graph_index'])
    if graph_name not in L.GRAPHS:
        raise SystemExit(f'{a.vcfdir}: graph {graph_name} has no off-reference nodes to test')
    kind = hdr['pansoma_kind']
    tensor_dir = hdr['pansoma_tensors']
    sample = L.sample_of(tensor_dir)
    chroms = hdr.get('pansoma_chromosomes', 'chr1').split(',')
    rows = {}   # (group, cid) -> dict

    def add(group, cid, **kw):
        rows[(group, cid)] = dict(group=group, candidate_id=cid, **kw)

    # unplaced PASS calls
    unplaced_all = {}
    for f in L.vcf_records(upath):
        inf = L.info(f[7])
        cid = first_id(f[2])
        unplaced_all[cid] = f[6]
        if f[6] == 'PASS':
            add('unplaced_pass', cid, node=int(f[0]), start=int(inf['START']), kind_event=inf['KIND'],
                evref=inf.get('EVREF', ''), evalt=inf.get('EVALT', ''), chrom=inf['BLOCK'],
                p_somatic=inf.get('P_SOMATIC'), offref_p_somatic=inf.get('OFFREF_P_SOMATIC', ''),
                af=f[9].split(':')[-1], filter=f[6])
    unplaced_pass_total = sum(g == 'unplaced_pass' for g, _ in rows)
    if unplaced_pass_total > a.max_unplaced:
        keep = set(random.Random(0).sample([k for k in rows if k[0] == 'unplaced_pass'], a.max_unplaced))
        rows = {k: v for k, v in rows.items() if k[0] != 'unplaced_pass' or k in keep}
    # placed controls
    rnd = random.Random(0)
    pon_calls = f'{a.vcfdir}/eval_nobed/vcfeval/pon_calls'
    ctrl = {}
    for g, fn in (('ctrl_tp', 'tp.vcf.gz'), ('ctrl_fp', 'fp.vcf.gz')):
        p = f'{pon_calls}/{fn}'
        if os.path.exists(p):
            ctrl[g] = [first_id(f[2]) for f in L.vcf_records(p)]
    ponv = [p for p in os.listdir(f'{a.vcfdir}/eval_bed') if p.endswith('.linear.pon.vcf.gz')] if os.path.isdir(f'{a.vcfdir}/eval_bed') else []
    if ponv:
        ctrl['ctrl_pon'] = [first_id(f[2]) for f in L.vcf_records(f'{a.vcfdir}/eval_bed/{ponv[0]}') if f[6] == 'PanelOfNormals']
    for g, ids in ctrl.items():
        for cid in rnd.sample(ids, min(a.controls, len(ids))):
            node, start, ev, allele = cid.split('@')[0].split(':')
            ref, alt = allele.split('>')
            add(g, cid, node=int(node), start=int(start), kind_event=ev, evref=ref, evalt=alt, chrom=chroms[0],
                p_somatic='', offref_p_somatic='', af='', filter='PASS')
    # gold keys on the tested chromosomes (exact candidate key, any label)
    import csv
    gold = {}
    gold_tsv = {'d9': L.GOLD, 'd46': f'{L.DATA}/gold_d46_mapping/offref_somatic_d46.tsv'}[graph_name]
    for r in (csv.DictReader(open(gold_tsv), delimiter='\t') if os.path.exists(gold_tsv) else []):
        if r['tensor_exact'] == 'True' and r['contig'].split('_')[0] in chroms:
            gold[r['tensor']] = r
    print(f'{a.name}: graph {graph_name}, kind {kind}, sample {sample}; unplaced PASS '
          f'{sum(g == "unplaced_pass" for g, _ in rows)}, controls {dict((g, len(v)) for g, v in ctrl.items())}', flush=True)

    # tensors
    ts = L.TensorSet(tensor_dir, kind)
    want = {cid for _, cid in rows} | set(gold)
    summ = {}
    for chrom in chroms:
        summ.update(ts.locations(chrom, want))   # labels ndjson; layout inferred from the tensor
    print(f'  tensor locations: {len(summ):,} of {len(want):,} wanted, {time.time() - t0:.0f}s', flush=True)
    for cid, r in gold.items():
        if cid in summ:   # gold keys in this graph's node coordinates
            node, start, ev, allele = cid.split(':')
            ref, alt = allele.split('>')
            add('gold', cid, node=int(node), start=int(start), kind_event=ev, evref=ref, evalt=alt,
                chrom=r['contig'].split('_')[0], p_somatic='', offref_p_somatic='', af='',
                filter=unplaced_all.get(cid, 'not_unplaced_record'), gold_event=r['event'])
    g = L.Graph(graph_name)
    hprc = L.HprcHaplotypes()
    for chrom in chroms:
        print(f'  focus {chrom}: {g.focus(chrom):,} node IDs in RAM, {time.time() - t0:.0f}s', flush=True)
    if len(rows) > 100:
        for chrom in chroms:
            print(f'  preload {chrom}: {g.preload(chrom):,} node IDs, {hprc.preload(chrom):,} HPRC top-level records, '
                  f'{time.time() - t0:.0f}s', flush=True)
    n = 0
    # node order = roughly GRCh38 order inside a chromosome component: sequential reads of every index
    for key, row in sorted(rows.items(), key=lambda kv: (kv[1]['chrom'], kv[1]['node'])):
        cid = row['candidate_id']
        rec = summ.get(cid)
        row.update(tensor=bool(rec))
        n += 1
        if not rec:
            continue
        process_row(row, rec, ts, g, hprc)
        if n % 50 == 0:
            print(f'  {n:,}/{len(rows):,} {time.time() - t0:.0f}s', flush=True)
    write_rows(list(rows.values()), f'{out}/calls.tsv')
    report = {}
    for ev in ('eval_nobed', 'eval_bed'):
        p = f'{a.vcfdir}/{ev}/report.json'
        if os.path.exists(p):
            rj = json.load(open(p))
            report[ev] = dict(pon=rj['pon']['at_threshold'], raw=rj['raw']['at_threshold'],
                              unplaced_pass=rj.get('unplaced_pass'), bed=rj.get('bed'), truth=rj.get('truth'),
                              baseline=rj['pon'].get('baseline'))
    json.dump(dict(name=a.name, vcfdir=a.vcfdir, unplaced_vcf=upath, graph=graph_name, kind=kind, sample=sample,
                   tensors=tensor_dir, chromosomes=chroms, controls={k: len(v) for k, v in ctrl.items()},
                   unplaced_records=len(unplaced_all), unplaced_pass=sum(v == 'PASS' for v in unplaced_all.values()),
                   unplaced_pass_processed=sum(g == 'unplaced_pass' for g, _ in rows),
                   report=report, seconds=round(time.time() - t0, 1)), open(f'{out}/run.json', 'w'), indent=1)
    print(f'{a.name}: {len(rows):,} rows in {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
