#!/usr/bin/env python3
"""Audit debug example tensors against indexed records, graph bases and row groups."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src')]
from indexed_gam_pipeline.candidates import BASES, rc
from indexed_gam_pipeline.gam_reader import IndexedGam
from indexed_gam_pipeline.run import node_records


def check(condition, message):
    if not condition:
        raise ValueError(message)


def validate(folder, gam, index, node_sqlite, walk_counts=None):
    folder=Path(folder)
    manifest=json.loads((folder/'manifest.json').read_text())
    metadata=[json.loads(s) for s in (folder/'variant_summary.ndjson').read_text().splitlines()]
    targets={m['node_id'] for m in metadata}
    graph_nodes=targets | {c['node_id'] for m in metadata for row in m['rows'] for c in row['columns'] if c}
    records=node_records(argparse.Namespace(node_json=None,gfa=None,node_sqlite=node_sqlite),graph_nodes)
    expected_walks = None
    if manifest['tensor_format_version'] == 'indexed-gam-candidate-v3':
        from indexed_gam_pipeline.walk_counts import WalkCounts
        check(walk_counts is not None, 'V3 validation requires --walk-counts')
        with WalkCounts(walk_counts) as lookup:
            expected_walks = lookup.get_counts(graph_nodes)
    covered={m['candidate_id']:Counter() for m in metadata}
    query={}
    # Independently count record overlap from original mapping intervals. This
    # deliberately does not call candidates.overlap or use tensor row counts.
    for alignment in IndexedGam(gam,index).fetch(targets,query):
        digest=hashlib.sha256(alignment.SerializeToString(deterministic=True)).hexdigest()
        for meta in metadata:
            if alignment.mapping_quality <= meta['parameters']['min_mapq']:
                continue
            hit=False
            for mapping in alignment.path.mapping:
                if mapping.position.node_id != meta['node_id']:
                    continue
                if not any(e.from_length or e.to_length for e in mapping.edit):
                    continue
                length=len(records[meta['node_id']]['sequence'])
                span=sum(e.from_length for e in mapping.edit)
                start=mapping.position.offset
                if mapping.position.is_reverse:
                    start=length-start-span
                end=start+span
                if meta['event_type']=='INS':
                    hit=start <= meta['start'] <= end
                else:
                    hit=start < meta['end'] and end > meta['start']
                if hit:
                    break
            if hit:
                covered[meta['candidate_id']][digest]+=1
    results=[]
    for meta in metadata:
        x=np.load(folder/f"shard_{meta['shard_index']:05d}_data.npy")[meta['index_within_shard']]
        label=meta['candidate_id']
        n=meta['selected_alignments']
        check(list(x.shape)==manifest['shape'] and str(x.dtype)==manifest['dtype'],f'{label}: tensor shape/dtype')
        check(len(meta['rows'])==n,f'{label}: debug row count')
        check(not x[:,n:].any(),f'{label}: unused row padding')
        check(meta['coverage']==sum(meta[f'{c}_count'] for c in ('alt','ref','other')),f'{label}: counts')
        check(meta['af']==meta['alt_count']/meta['coverage'],f'{label}: AF')
        check(sum(covered[label].values())==meta['coverage'],f'{label}: independent source coverage')
        selected=Counter(row['record_sha256'] for row in meta['rows'])
        check(not selected-covered[label],f'{label}: selected source records')
        if n==meta['coverage']:
            check(selected==covered[label],f'{label}: full source record multiset')
        check(Counter(row['support'] for row in meta['rows'])==Counter(meta['selected_counts']),f'{label}: selected support counts')
        lo,hi=meta['candidate_columns']
        check(bool(np.all(x[2,:n,lo:hi]&2)),f'{label}: full candidate flags')
        if meta['event_type']=='INS':
            check(bool(np.all(x[5,:n,lo:hi]==6)),f'{label}: insertion reference gaps')
        check(bool(np.all(x[1,:n][x[0,:n]==6]==-1)),f'{label}: gap quality')
        actual_paths=[]
        checked=0
        for ri,row in enumerate(meta['rows']):
            path=[]; previous=None
            for ci,col in enumerate(row['columns']):
                if col is None:
                    if expected_walks is not None:
                        check(int(x[6,ri,ci]) == 0, f'{label}: padded walk count')
                    continue
                if expected_walks is not None:
                    check(int(x[6,ri,ci]) == expected_walks[col['node_id']], f'{label}: walk count at row {ri}, column {ci}')
                if col['mapping_index'] != previous:
                    path.append((col['node_id'],col['reverse']))
                    previous=col['mapping_index']
                if 'offset' in col:
                    base=records[col['node_id']]['sequence'][col['offset']]
                    if col['reverse']:
                        base=rc(base)
                    check(int(x[5,ri,ci])==BASES.get(base,5),f'{label}: graph base at row {ri}, column {ci}')
                    checked+=1
            actual_paths.append(tuple(path))
        check(actual_paths==sorted(actual_paths),f'{label}: node-path sorting')
        next_row=0
        for group in meta['row_groups']:
            check(group['start_row']==next_row and group['end_row']>next_row,f'{label}: group boundaries')
            key=tuple((v['node_id'],v['reverse']) for v in group['path'])
            check(all(p==key for p in actual_paths[next_row:group['end_row']]),f'{label}: group membership')
            next_row=group['end_row']
        check(next_row==n,f'{label}: all rows grouped')
        results.append(dict(candidate_id=label,event_type=meta['event_type'],coverage=meta['coverage'],
            selected_alignments=n,node_path_groups=len(meta['row_groups']),
            reference_columns_checked=checked,passed=True))
    check(len(results)==manifest['tensors'],'Manifest tensor count')
    return dict(passed=True,examples=len(results),event_types=dict(Counter(r['event_type'] for r in results)),
        reference_columns_checked=sum(r['reference_columns_checked'] for r in results),
        selected_rows_checked=sum(r['selected_alignments'] for r in results),
        indexed_query=query,source_gam=gam,source_index=index,graph_node_sqlite=node_sqlite,
        checks=['shard/manifest shape and dtype','counts and AF','independent raw-mapping coverage',
                'selected records versus original GAM multiset','padding and gap qualities',
                'candidate flags and insertion reference gaps','oriented graph-reference bases',
                'node-path sorting and group memberships'],
        walk_counts_checked=expected_walks is not None, walk_count_cache=walk_counts,
        candidates=results,full_genome_run=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder')
    parser.add_argument('--gam',required=True)
    parser.add_argument('--index')
    parser.add_argument('--node-sqlite',required=True)
    parser.add_argument('--walk-counts')
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    report=validate(args.folder,args.gam,args.index,args.node_sqlite,args.walk_counts)
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
