"""The off-reference gold set on another graph (d46): the assembly somatic events of assembly_offref_20260929 (normal
haplotype window N and the same window with the somatic allele T, 100 bp flanks, windows.fq) aligned with vg giraffe
to that graph, classified as in assembly_offref_20260929/s3_classify.py and keyed as in s4_tensors.py:
  A  a somatic edit (a difference of T's alignment that N's lacks) on a GRCh38 node
  B  the somatic edit on a node off the GRCh38 path (off-reference somatic: the gold set)
  C  no edit, T and N walk different nodes (the allele is a graph branch); D unmapped / MAPQ < 20 / same alignment
Class B events get their node-forward candidate key and the tensor of that key in the given tensor sets.

usage: gold_graph.py GRAPH GAM_JSON OUT_TSV SET [SET ...]
"""
import csv
import json
import re
import sys
from collections import Counter

sys.path.insert(0, '/scratch/jshen/Github/Pansoma/analysis/reviewer_response_20261011/experiments/E9c_offref_filter_validation')
import e9c_lib as L  # noqa: E402

A = '/scratch/jshen/data/pansoma_net_v2_runs/assembly_offref_20260929'


def walk(a):
    out, cur = [], 0
    for m in a.get('path', {}).get('mapping', []):
        p = m['position']
        node, rev, off = int(p['node_id']), bool(p.get('is_reverse', False)), int(p.get('offset', 0))
        for e in m.get('edit', []):
            fl, tl, seq = int(e.get('from_length', 0)), int(e.get('to_length', 0)), e.get('sequence', '')
            out.append((cur, cur + tl, node, rev, off, fl, tl, seq, fl != tl or bool(seq)))
            cur += tl
            off += fl
    return out


def nodes_at(w, x):
    return {s[2] for s in w if s[0] <= x < s[1] or (s[6] == 0 and s[0] == x)}


def main():
    graph, gam_json, out_tsv, sets = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
    g = L.Graph(graph, in_memory=False)
    on_ref = lambda n: g.on_ref(n) and int(g.rp.visits[n]) == 1   # noqa: E731
    flank = {}
    with open(f'{A}/windows.fq') as f:
        for k, line in enumerate(f):
            if k % 4 == 0:
                name, _, fl = line[1:].split()
                flank[name] = int(fl.split('=')[1])
    events = {r['event']: r for r in csv.DictReader(open(f'{A}/events.tsv'), delimiter='\t')}
    aln = {}
    for line in open(gam_json):
        a = json.loads(line)
        aln[a['name']] = a
    rows, cls = [], Counter()
    for eid, ev in events.items():
        t, n = aln.get(f'{eid}_T'), aln.get(f'{eid}_N')
        F = flank[f'{eid}_T']
        Ln = max(1, len(ev['alt']))
        row = dict(event=eid, contig=ev['contig'], pos=ev['pos'], ref=ev['ref'], alt=ev['alt'],
                   grch38_records=ev['grch38_records'], cls='', edit_node='', edit_rev='', edit_offset='', edit='',
                   background_on_ref='', mapq_T=(t or {}).get('mapping_quality', 0), mapq_N=(n or {}).get('mapping_quality', 0))
        if not t or not n or not t.get('path') or not n.get('path'):
            row['cls'] = 'D_unmapped'
        elif int(row['mapq_T']) < 20 or int(row['mapq_N']) < 20:
            row['cls'] = 'D_mapq<20'
        else:
            wt, wn = walk(t), walk(n)
            bg = nodes_at(wn, F)
            row['background_on_ref'] = all(on_ref(x) for x in bg) if bg else ''
            ndiff = {(s[2], s[4], s[5], s[7]) for s in wn if s[8]}
            tonly = [s for s in wt if s[8] and (s[2], s[4], s[5], s[7]) not in ndiff]
            som = sorted(tonly, key=lambda s: abs((s[0] + s[1]) / 2 - (F + Ln / 2)))
            if som:
                s = som[0]
                row.update(cls='A' if on_ref(s[2]) else 'B', edit_node=s[2], edit_rev=s[3], edit_offset=s[4],
                           edit=f'{s[5]}>{s[6]}:{s[7]}')
            elif nodes_at(wt, F) != bg or [s[2] for s in wt] != [s[2] for s in wn]:
                row['cls'] = 'C'
            else:
                row['cls'] = 'D_same_alignment'
        rows.append(row)
        cls[row['cls']] += 1
    print(graph, 'events', len(rows), dict(sorted(cls.items())), flush=True)
    # class B -> node-forward key
    B = [r for r in rows if r['cls'] == 'B']
    seqs = g.seqs([int(r['edit_node']) for r in B])
    for r in B:
        node, rev, off = int(r['edit_node']), r['edit_rev'] is True, int(r['edit_offset'])
        fl_tl, seq = r['edit'].split(':', 1)
        fl, tl = map(int, fl_tl.split('>'))
        s = seqs[node]
        Ln = len(s)
        if fl == tl == 1:
            start = Ln - off - 1 if rev else off
            kind, ref, alt = 'SNP', s[start], (seq.translate(L.COMP) if rev else seq)
        elif fl == 0:
            start = Ln - off if rev else off
            kind, ref, alt = 'INS', '', (L.rc(seq) if rev else seq)
        elif tl == 0:
            start = Ln - off - fl if rev else off
            kind, ref, alt = 'DEL', s[start:start + fl], ''
        else:
            start, kind, ref, alt = (Ln - off - fl if rev else off), 'COMPLEX', '', ''
        r.update(node=node, start=start, kind=kind, key=f'{node}:{start}:{kind}:{ref}>{alt}', chrom=r['contig'].split('_')[0])
    cid_re = re.compile(r'"candidate_id": ?"([^"]+)"')
    keys = {r['key']: r for r in B}
    found = {}
    for st in sets:
        for kind_dir in ('SNV', 'INDEL'):
            for chrom in sorted({r['chrom'] for r in B}):
                try:
                    f = open(f'{L.TENSORS}/{st}/tensors/{kind_dir}/{chrom}_labels.ndjson')
                except FileNotFoundError:
                    continue
                with f:
                    for line in f:
                        cid = cid_re.search(line).group(1)
                        if cid.split('@')[0] in keys and cid not in found:
                            lab = json.loads(line)
                            found[cid] = (st, lab['label'], lab['reason'])
        print(st, 'exact keys found so far', len(found), flush=True)
    with open(out_tsv, 'w') as w:
        cols = ['event', 'contig', 'pos', 'ref', 'alt', 'grch38_records', 'normal_offref_at_site', 'node', 'start', 'kind',
                'key', 'tensor', 'tensor_exact', 'tensor_label', 'tensor_reason', 'tensor_sets']
        w.write('\t'.join(cols) + '\n')
        for r in B:
            hit = [(c, v) for c, v in found.items() if c.split('@')[0] == r['key']]
            c, v = hit[0] if hit else ('', ('', '', ''))
            w.write('\t'.join(map(str, [r['event'], r['contig'], r['pos'], r['ref'], r['alt'], r['grch38_records'],
                                        r['background_on_ref'] is False, r['node'], r['start'], r['kind'], r['key'], c,
                                        bool(hit), v[1], v[2], ','.join(sorted({vv[0] for cc, vv in hit}))])) + '\n')
    with open(out_tsv.replace('.tsv', '.classes.tsv'), 'w') as w:
        cols = list(rows[0])
        w.write('\t'.join(cols) + '\n')
        for r in rows:
            w.write('\t'.join(str(r.get(k, '')) for k in cols) + '\n')
    print(f'{graph}: class B {len(B)}, with an exact tensor key {sum(1 for r in B if any(c.split("@")[0] == r["key"] for c in found))}')


if __name__ == '__main__':
    main()
