"""Prepare a frozen full run using an existing complete discovery, without rescanning GAM."""
import argparse
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time
from indexed_gam_pipeline.full_run import stamp, write_json
from indexed_gam_pipeline.unified_full_run import digest, partition_nodes

ROOT=Path(__file__).resolve().parents[1]


def prepare(prior_config, output, processes, builder, query):
    start=time.perf_counter()
    old=json.loads(Path(prior_config).read_text())
    c=old['config'].copy()
    for key in ('gam','index','gbz'):
        if stamp(c[key]) != old['inputs'][key]:
            raise ValueError('Discovery run input has changed: '+key)
    discovery=Path(old['prior_run'])/'discovery'
    report=json.loads((discovery/'discovery_report.json').read_text())
    if report.get('exploratory') or report.get('max_alignments') or report.get('max_nodes') or report['gam']!=c['gam']:
        raise ValueError('Cannot reuse incomplete or mismatched discovery')
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
    source=output/'source';source.mkdir()
    for folder in ('indexed_gam_pipeline','src'):
        shutil.copytree(ROOT/folder,source/folder,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    for name,path in [('gbz_graph_index',builder),('gbz_node_counts',query)]:
        shutil.copy2(path,output/name)
    c.update(unified_graph_index=str(output/'all_graph_nodes.gbz.sqlite'),
             gbz_index_builder=str(output/'gbz_graph_index'),gbz_query=str(output/'gbz_node_counts'),
             discovery_nodes=str(discovery/'target_nodes.txt'),discovery_report=str(discovery/'discovery_report.json'))
    parts=partition_nodes(c['discovery_nodes'],output,processes,report['nodes_selected'])
    config=dict(config=c,prior_run=old['prior_run'],parts=parts,total_nodes=report['nodes_selected'],
                inputs={k:stamp(c[k]) for k in ('gam','index','gbz','discovery_nodes','discovery_report')},
                params=dict(gam_cache_mb_per_process=8192,processes=processes,early_alt_filter=True),
                source_sha256={str(p):digest(p) for p in source.rglob('*') if p.is_file()},
                source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                discovery_reused=True,preparation_seconds=time.perf_counter()-start)
    config['source_sha256'].update({str(output/name):digest(output/name) for name in ('gbz_graph_index','gbz_node_counts')})
    write_json(output/'config.json',config)
    script='\n'.join(['#!/bin/bash','set -euo pipefail',
        'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1',
        'export PYTHONUNBUFFERED=1',f'cd {shlex.quote(str(source))}',
        'exec /usr/bin/time -v -o '+shlex.quote(str(output/'job.resources.txt'))+' '+shlex.quote(sys.executable)+
        ' -m indexed_gam_pipeline.unified_full_run --root '+shlex.quote(str(output)),''])
    (output/'run.sh').write_text(script)
    print(json.dumps(dict(root=str(output),nodes=config['total_nodes'],processes=processes,
                         min_nodes=min(p['nodes'] for p in parts),max_nodes=max(p['nodes'] for p in parts),
                         preparation_seconds=config['preparation_seconds']),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prior-config',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--processes',type=int,default=20)
    p.add_argument('--builder',required=True)
    p.add_argument('--query',required=True)
    a=p.parse_args();prepare(a.prior_config,a.output,a.processes,a.builder,a.query)
