"""E5d step 3 (HG008T): tumor read alleles at the set alleles not covered by the 2026-10-01 biology audit (b2), with the
audit's own code: analysis data dir graph_absorbed_somatic_20261001/audit/biology/b2_tumor_reads.py work() (imported, not
copied: linear GRCh38 BAMs of HG008-T, Illumina 195x and PacBio HiFi 116x; repeat window widened to the overlapping
HG008-N dipcall records; <= 300 reads with MAPQ >= 20 per locus). The audit's 400 + 400 random non-absorbed controls are
reused from its b2_tumor_reads.tsv.
Then the b4 rule (b4_tumor_alleles.py, re-implemented line for line below) on old + new rows: EXP = normal event-hap
length + the normal-frame event length (b5_normal_frame.tsv nf_len, every truth); EXP 'seen' = >= max(3, 15 %) of the
reads; loci with >= 10 reads and a length-changing event only.
Outputs ($DATA/hg008t/): b2_new.tsv (new loci), b4_all.tsv (old + new, one row per locus x platform).
Run: python e3_hg008t_tumor_reads.py [workers]   (Slurm: e3_tumor_reads.sbatch)
"""
import collections, csv, importlib.util, multiprocessing as mp, sys

OLD = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
BIO = f'{OLD}/audit/biology'
OUT = '/scratch/jshen/data/reviewer_experiments_20261013/E5d_absorbed_normal_d46/hg008t'
spec = importlib.util.spec_from_file_location('b2', f'{BIO}/b2_tumor_reads.py')
b2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b2)


def work(t):
    return b2.work(t)   # a picklable top-level wrapper (the audit module is loaded from a path)


def b4_row(r, b5):
    """b4_tumor_alleles.py, per (locus, platform) row of b2."""
    t, n = r['truth_id'], int(r['n'])
    o = dict(truth_id=t, platform=r['platform'], set=r['set'], kind2=r['kind2'], n=n, status='')
    if n < 10:
        o['status'] = 'lt10_reads'
        return o
    if r['NE_len'] in ('None', '') or r['NO_len'] in ('None', ''):
        o['status'] = 'normal_hap_unresolved'
        return o
    prof = collections.Counter()
    for x in r['top_len'].split(';'):
        if x:
            d, k = x.rsplit(':', 1)
            prof[int(d)] += int(k)
    seen = lambda L: prof[L] >= max(3, 0.15 * n)  # noqa: E731
    NE = int(r['NE_len'])
    nf = b5[t]['nf_len']
    if nf in ('', None) or int(nf) == 0:
        o['status'] = 'event_not_length_changing'
        return o
    EXP = NE + int(nf)
    o.update(EXP=EXP, f_EXP=round(prof[EXP] / n, 3), EXP_seen=seen(EXP), status='ok')
    return o


def main(workers):
    loci = list(csv.DictReader(open(f'{OUT}/loci.tsv'), delimiter='\t'))
    todo = []
    for L in loci:
        ev = L['asm_event']
        eh = ev.split(':')[0].split('_')[-1] if ev else ''
        todo.append(dict(L, set='absorbed', category='e5d', event_hap=eh))
    with mp.Pool(workers) as pool:
        res = pool.map(work, todo, chunksize=4)
    new = [o for rr in res for o in rr]
    cols = list(new[0])
    with open(f'{OUT}/b2_new.tsv', 'w') as w:
        w.write('\t'.join(cols) + '\n')
        for o in new:
            w.write('\t'.join(str(o[k]) for k in cols) + '\n')
    old = list(csv.DictReader(open(f'{BIO}/b2_tumor_reads.tsv'), delimiter='\t'))
    b5 = {r['truth_id']: r for r in csv.DictReader(open(f'{BIO}/b5_normal_frame.tsv'), delimiter='\t')}
    newids = {L['truth_id'] for L in loci}
    rows = [r for r in old if not (r['set'] == 'absorbed' and r['truth_id'] in newids)] + \
           [{k: str(v) for k, v in o.items()} for o in new]
    out = [b4_row(r, b5) for r in rows]
    cols = ['truth_id', 'platform', 'set', 'kind2', 'n', 'status', 'EXP', 'f_EXP', 'EXP_seen']
    with open(f'{OUT}/b4_all.tsv', 'w') as w:
        w.write('\t'.join(cols) + '\n')
        for o in out:
            w.write('\t'.join(str(o.get(k, '')) for k in cols) + '\n')
    print('new loci', len(loci), 'rows', len(out))


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 8)
