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


def validate(folder, gam, index, node_sqlite=None, walk_counts=None, gbz=None, gbz_query=None, occurrence_cache=None, graph_index=None):
    folder=Path(folder)
    manifest=json.loads((folder/'manifest.json').read_text())
    metadata=[json.loads(s) for s in (folder/'variant_summary.ndjson').read_text().splitlines()]
    targets={m['node_id'] for m in metadata}
    graph_nodes=targets | {c['node_id'] for m in metadata for row in m['rows'] for c in row['columns'] if c}
    expected_walks = None
    if graph_index:
        from indexed_gam_pipeline.graph_index import GraphIndex
        check(not node_sqlite, 'Choose graph-index or node-sqlite')
        with GraphIndex(graph_index) as lookup:
            records = lookup.get_nodes(graph_nodes)
        expected_walks = {n:r['distinct_path_count'] for n,r in records.items()}
    else:
        check(node_sqlite is not None, 'Validation requires graph-index or node-sqlite')
        records=node_records(argparse.Namespace(node_json=None,gfa=None,node_sqlite=node_sqlite),graph_nodes)
    if manifest['tensor_format_version'] == 'indexed-gam-candidate-v3':
        from indexed_gam_pipeline.walk_counts import WalkCounts
        check(walk_counts is not None, 'V3 validation requires --walk-counts')
        with WalkCounts(walk_counts) as lookup:
            expected_walks = lookup.get_counts(graph_nodes)
    if manifest['tensor_format_version'] == 'indexed-gam-candidate-v4' and (not graph_index or gbz):
        from indexed_gam_pipeline.gbz_counts import GBZCounts
        check(all((gbz, gbz_query, occurrence_cache)), 'V4 validation requires GBZ, helper and cache')
        with GBZCounts(gbz, gbz_query, occurrence_cache) as lookup:
            expected_walks = lookup.get_counts(graph_nodes, {n:r['sequence'] for n,r in records.items()})
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
        x=np.load(folder/f"shard_{meta['shard_index']:05d}_data.npy", mmap_mode="r")[meta['index_within_shard']]
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
        anchor_window = meta.get('window_encoding_version') == 'anchor-centered-columns-v1'
        if not anchor_window:
            check(bool(np.all(x[2,:n,lo:hi]&2)),f'{label}: full candidate flags')
        if meta['event_type']=='INS' and not anchor_window:
            check(bool(np.all(x[5,:n,lo:hi]==6)),f'{label}: insertion reference gaps')
        check(bool(np.all(x[1,:n][x[0,:n]==6]==-1)),f'{label}: gap quality')
        actual_paths=[]
        checked=0
        for ri,row in enumerate(meta['rows']):
            path=[]; previous=None
            for ci,col in enumerate(row['columns']):
                if col is None:
                    if anchor_window:
                        check(not x[:,ri,ci].any(),f'{label}: all-zero missing evidence')
                    if expected_walks is not None:
                        check(int(x[6,ri,ci]) == 0, f'{label}: padded walk count')
                    continue
                if anchor_window:
                    anchored = col['mapping_index'] == row['anchor_mapping_index']
                    boundary = 'boundary' in col
                    pos = col.get('boundary', col.get('offset'))
                    if meta['event_type'] == 'INS':
                        expected_flag = anchored and boundary and pos == meta['start']
                    else:
                        expected_flag = anchored and (meta['start'] < pos < meta['end'] if boundary
                                                      else meta['start'] <= pos < meta['end'])
                    check(bool(x[2,ri,ci]&2)==expected_flag,f'{label}: per-row candidate flags')
                    if ci == meta['anchor_column'] and not boundary:
                        check(anchored and pos == meta['start'],f'{label}: anchor coordinate')
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
        if meta.get('row_selection_version') == 'window-edit-bp-group-uniform-v1':
            audit = meta['selection_audit']
            check(Counter(a['record_sha256'] for a in audit)==covered[label],
                  f'{label}: all preselection records')
            ordered_groups = {}
            for a in sorted(audit, key=lambda a: (-a['window_mismatch_bp'],
                    -a['mapping_quality'], a['record_sha256'], a['anchor_mapping_index'])):
                key = tuple(tuple(p) for p in a['path'])
                ordered_groups.setdefault(key, []).append(a)
            ranked = [a for group in ordered_groups.values() for a in group]
            k = min(manifest['shape'][1], len(ranked))
            indices = ([len(ranked)//2] if k == 1 else
                       [i*(len(ranked)-1)//(k-1) for i in range(k)] if k else [])
            check(meta['selected_grouped_ranks']==indices, f'{label}: uniform sampling ranks')
            check(n==k, f'{label}: sampled depth')
            for ri, idx in enumerate(indices):
                a = ranked[idx]
                row = meta['rows'][ri]
                check(row['record_sha256']==a['record_sha256'] and
                      row['anchor_mapping_index']==a['anchor_mapping_index'] and
                      row['support']==a['support'], f'{label}: sampled record ordering')
                check(actual_paths[ri]==tuple(tuple(p) for p in a['path']),
                      f'{label}: sampled path')
                mismatch = int(np.count_nonzero((x[0,ri] != 0) & (x[5,ri] != 0) &
                                               (x[0,ri] != x[5,ri]) & (x[4,ri] != 6)))
                check(mismatch==a['window_mismatch_bp']==row['window_mismatch_bp']==
                      meta['window_mismatch_bp'][ri], f'{label}: visible edit bp count')
        else:
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
        indexed_query=query,source_gam=gam,source_index=index,graph_node_sqlite=node_sqlite,graph_index=graph_index,
        checks=['shard/manifest shape and dtype','counts and AF','independent raw-mapping coverage',
                'selected records versus original GAM multiset','padding and gap qualities',
                'candidate flags and insertion reference gaps','oriented graph-reference bases',
                'node-path ordering and group memberships',
                *(['window edit bp and uniform sampling from recorded preselection windows']
                  if manifest.get('row_selection_version') else [])],
        walk_counts_checked=expected_walks is not None, walk_count_cache=walk_counts, gbz=gbz, occurrence_cache=occurrence_cache,
        candidates=results,full_genome_run=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder')
    parser.add_argument('--gam',required=True)
    parser.add_argument('--index')
    parser.add_argument('--node-sqlite')
    parser.add_argument('--graph-index')
    parser.add_argument('--walk-counts')
    parser.add_argument('--gbz')
    parser.add_argument('--gbz-query')
    parser.add_argument('--occurrence-cache')
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    report=validate(args.folder,args.gam,args.index,args.node_sqlite,args.walk_counts,args.gbz,args.gbz_query,args.occurrence_cache,args.graph_index)
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
