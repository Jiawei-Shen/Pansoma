"""Gold set of off-reference somatic variants: the HG008-T somatic events that the HG008-N v6.2 assembly places on a
non-reference d9 node (assembly_offref_20260929/offref_somatic.tsv, 233 with an exact HG008T_Illumina tensor key),
on every chromosome. For each HG008T d9 tensor set that holds the key: contexts, bubble locus and the graph-native
population tags of e9c_run.process_row, i.e. the false-tag rate of the graph PoN on known somatic variants.

usage: gold.py [--graph d46 --gold GOLD_TSV] SET [SET ...]      e.g. HG008T_Illumina HG008T_PacBio HG008T_ONT
(d46: the gold TSV of gold_graph.py, i.e. the same events aligned to d46 and keyed in d46 node coordinates)
Output: DATA/gold/SET/calls.tsv, run.json (sample HG008T, so e9d_asm.py also checks the contexts in the assemblies).
"""
import csv
import json
import os
import sys
import time
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e9c_lib as L  # noqa: E402
from e9c_run import process_row, write_rows  # noqa: E402


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--graph', default='d9')
    ap.add_argument('--gold', default=L.GOLD)
    ap.add_argument('sets', nargs='+')
    a = ap.parse_args()
    gold = [r for r in csv.DictReader(open(a.gold), delimiter='\t') if r['tensor_exact'] == 'True']
    g = L.Graph(a.graph)
    hprc = L.HprcHaplotypes()
    for name in a.sets:
        t0 = time.time()
        out = f'{L.DATA}/gold/{name}'
        os.makedirs(out, exist_ok=True)
        tdir = f'{L.TENSORS}/{name}/tensors'
        by = defaultdict(lambda: defaultdict(dict))   # kind dir -> chrom -> cid -> gold row
        for r in gold:
            kind = 'SNV' if ':SNP:' in r['tensor'] else 'INDEL'
            by[kind][r['contig'].split('_')[0]][r['tensor']] = r
        rows = []
        for kind, chroms in by.items():
            ts = L.TensorSet(tdir, kind)
            for chrom, cids in sorted(chroms.items()):
                summ = ts.locations(chrom, set(cids))   # labels ndjson; layout inferred from the tensor
                for cid, r in sorted(cids.items(), key=lambda kv: int(kv[0].split(':')[0])):
                    node, start, ev, allele = cid.split(':')
                    ref, alt = allele.split('>')
                    row = dict(group='gold', candidate_id=cid, node=int(node), start=int(start), kind_event=ev,
                               evref=ref, evalt=alt, chrom=chrom, gold_event=r['event'],
                               gold_assembly=f"{r['contig']}:{r['pos']}:{r['ref']}>{r['alt']}",
                               normal_offref_at_site=r['normal_offref_at_site'], label_reason=r['tensor_reason'],
                               tensor=cid in summ)
                    if cid in summ:
                        process_row(row, summ[cid], ts, g, hprc)
                    rows.append(row)
            print(f'{name} {kind}: {sum(len(c) for c in chroms.values())} gold keys, {time.time() - t0:.0f}s', flush=True)
        write_rows(rows, f'{out}/calls.tsv')
        json.dump(dict(name=f'gold__{name}', set=name, graph=a.graph, gold_tsv=a.gold, sample='HG008T', gold_keys=len(gold),
                       with_tensor=sum(r['tensor'] for r in rows), seconds=round(time.time() - t0, 1)),
                  open(f'{out}/run.json', 'w'), indent=1)
        print(f'{name}: {sum(r["tensor"] for r in rows)} of {len(rows)} gold keys have a tensor; {time.time() - t0:.0f}s',
              flush=True)


if __name__ == '__main__':
    main()
