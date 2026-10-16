"""E5b step 4: is the truth ALT haplotype a path of the graph at its locus? (d9 or d46; one graph per run)

For every locus of loci.tsv.gz with a window on the graph: ref_hap = GRCh38 from the first to the last anchor node,
alt_hap = ref_hap with the truth allele applied. A path from (first,+) to (last,+), whose GRCh38 nodes all lie inside
the window, that spells alt_hap exactly = the graph encodes the allele ('spelled'; graph_absorbed_somatic_20261001 s2's
'exact' class, without its germline-background and closest-allele searches). The GRCh38 path must spell ref_hap
(checked: ref_ok; a failure means the subgraph is incomplete).
For a spelled ALT: its departures from GRCh38 (elements: branch = run of non-GRCh38 nodes, skip = edge over GRCh38 bases),
and the smallest distinct_path_count (haplotypes, pipeline graph_index sqlite) over its non-GRCh38 nodes
(alt_min_paths; empty when the ALT uses only skip edges).
Inputs: loci.tsv.gz (s2), subgraph_<graph>.gfa (s3_extract.sh). Output: spell_<graph>.tsv.gz.
Run: python s4_spell.py d9|d46 [workers]
"""
import csv, gzip, sqlite3, sys, time
from multiprocessing import Pool
import pysam
from e5b_common import DATA, FASTA, GRAPHS, RefPath, SubGraph, spells, elements

g = sys.argv[1]
W = int(sys.argv[2]) if len(sys.argv) > 2 else 1
R = G = fa = db = None


def init():
    global fa, db
    fa = pysam.FastaFile(FASTA)
    db = sqlite3.connect(f"file:{GRAPHS[g]['sqlite']}?mode=ro", uri=True)


def one(r):
    out = dict(sample=r['sample'], truth_id=r['truth_id'])
    if r[f'{g}_status'] != 'ok':
        out['status'] = r[f'{g}_status']
        return out
    chrom, pos0, vref, valt = r['chrom'], int(r['vcf_pos']) - 1, r['vcf_ref'].upper(), r['vcf_alt'].upper()
    first, last = int(r[f'{g}_first']), int(r[f'{g}_last'])
    lo, hi = int(r[f'{g}_win_start0']), int(r[f'{g}_win_end0'])
    region = fa.fetch(chrom, lo, hi + 1).upper()
    off = pos0 - lo
    assert region[off:off + len(vref)] == vref
    alt_hap = region[:off] + valt + region[off + len(vref):]
    ok_ref, _, info_ref = spells(G, R, chrom, first, last, lo, hi, region)
    ok_alt, path, info = spells(G, R, chrom, first, last, lo, hi, alt_hap)
    out.update(status='ok', ref_ok=int(ok_ref), alt_spelled=int(ok_alt), info=info or info_ref)
    if ok_alt:
        el = elements(R, path)
        out['elements'] = ';'.join(':'.join(map(str, e)) for e in el)
        br = [n for n, o in path if R.ref(n) is None]
        if br:
            q = ','.join('?' * len(br))
            c = [x[0] for x in db.execute(f'select distinct_path_count from nodes where node_id in ({q})', br)]
            out['alt_min_paths'] = min(c) if c else ''
            out['alt_branch_nodes'] = len(br)
    return out


if __name__ == '__main__':
    t0 = time.time()
    R = RefPath(g)
    G = SubGraph(f'{DATA}/subgraph_{g}.gfa')
    print(g, 'subgraph', len(G.seq), 'nodes', round(time.time() - t0), 's', file=sys.stderr, flush=True)
    rows = list(csv.DictReader(gzip.open(f'{DATA}/loci.tsv.gz', 'rt'), delimiter='\t'))
    cols = ['sample', 'truth_id', 'status', 'ref_ok', 'alt_spelled', 'info', 'elements', 'alt_min_paths', 'alt_branch_nodes']
    with gzip.open(f'{DATA}/spell_{g}.tsv.gz', 'wt') as w, Pool(W, initializer=init) as pool:   # fork: R and G shared
        w.write('\t'.join(cols) + '\n')
        for i, o in enumerate(pool.imap(one, rows, chunksize=200)):
            w.write('\t'.join(str(o.get(c, '')) for c in cols) + '\n')
            if (i + 1) % 10000 == 0:
                print(i + 1, round(time.time() - t0), 's', file=sys.stderr, flush=True)
    print('done', len(rows), round(time.time() - t0), 's', file=sys.stderr)
