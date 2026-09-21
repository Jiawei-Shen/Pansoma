#!/usr/bin/env python3
"""Controlled serial/parallel candidate benchmark, with exact tensor equivalence."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from indexed_gam_pipeline.full_run import write_json,validate_shards


def process_memory(pid):
    """Sample descendants externally, never traverse Python caches in the hot loop."""
    todo=[pid];seen=set();rss=pss=0
    while todo:
        current=todo.pop()
        if current in seen:continue
        seen.add(current)
        root=Path('/proc')/str(current)
        try:
            for task in (root/'task').iterdir():
                try:todo.extend(map(int,(task/'children').read_text().split()))
                except (OSError,ValueError):pass
            for line in (root/'status').read_text().splitlines():
                if line.startswith('VmRSS:'):rss+=int(line.split()[1])
            for line in (root/'smaps_rollup').read_text().splitlines():
                if line.startswith('Pss:'):pss+=int(line.split()[1])
        except OSError:pass
    return rss,pss


def digest_file(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for data in iter(lambda:stream.read(4*1024*1024),b''):h.update(data)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prior-run',type=Path,required=True)
    p.add_argument('--verified-test',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();prior=args.prior_run.resolve();test=args.verified_test.resolve();out=args.output.resolve()
    cfg=json.loads((prior/'config.json').read_text());out.mkdir(parents=True,exist_ok=True)
    with (prior/'discovery/target_nodes.txt').open() as stream:
        nodes=[next(stream).strip() for _ in range(1024)]
    (out/'nodes.txt').write_text('\n'.join(nodes)+'\n')
    report=dict(status='running',nodes=1024,width=101,shard_size=2048,steps=[],
                original_full_job='362252 remains SIGSTOP; benchmark does not resume it')
    def save():write_json(out/'status.json',report)
    def run(name,command):
        item=dict(name=name,status='running',started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),command=list(map(str,command)))
        report['steps'].append(item);save();start=time.perf_counter();peak_rss=peak_pss=0;sample_seconds=0
        with (out/(name+'.log')).open('w') as log:
            proc=subprocess.Popen(['/usr/bin/time','-v','-o',str(out/(name+'.resources.txt')),*map(str,command)],stdout=log,stderr=subprocess.STDOUT,cwd=ROOT)
            while proc.poll() is None:
                t=time.perf_counter();rss,pss=process_memory(proc.pid);sample_seconds+=time.perf_counter()-t
                peak_rss=max(peak_rss,rss);peak_pss=max(peak_pss,pss)
                try:proc.wait(timeout=5)
                except subprocess.TimeoutExpired:pass
        item.update(wall_seconds=time.perf_counter()-start,returncode=proc.returncode,
                    sampled_peak_tree_rss_kib=peak_rss,sampled_peak_tree_pss_kib=peak_pss,
                    memory_sampling_seconds=sample_seconds,status='complete' if proc.returncode==0 else 'failed')
        save()
        if proc.returncode:raise RuntimeError(name+' failed')
        return item
    cases=[('baseline_before',False,0,1),('early_alt_only',True,0,1),
           ('indexed_one_worker',True,1000,1),('indexed_two_workers',True,1000,2),
           ('baseline_after',False,0,1)]
    baseline_source=test/'source'
    try:
        save()
        for name,early,index_nodes,workers in cases:
            source=baseline_source if name.startswith('baseline') else ROOT
            cache=out/(name+'_gbwt.sqlite');shutil.copyfile(test/'gbwt_counts.sqlite',cache)
            cmd=[sys.executable,source/'indexed_gam_pipeline/run.py','build','--format','candidate-v4',
                '--gam',cfg['gam'],'--index',cfg['index'],'--nodes',out/'nodes.txt',
                '--node-sqlite',prior/'all_graph_nodes.sqlite','--gbz',cfg['gbz'],'--gbz-query',cfg['gbz_query'],
                '--occurrence-cache',cache,'--output',out/name,'--batch-nodes','512','--max-node-span','10000',
                '--gam-cache-mb','8192','--max-batch-segments','20000','--shard-size','2048',
                '--rows','200','--width','101','--min-mapq','10','--min-af','0.05','--min-variants','3',
                '--min-allele-bq','10','--max-indel-len','50','--variant-type','all','--debug-rows']
            if not name.startswith('baseline'):
                cmd+=['--workers',str(workers),'--node-index-cache-nodes',str(index_nodes),'--node-index-cache-mb','64']
                if early:cmd.append('--early-alt-filter')
            item=run(name,cmd)
            item['validation']=validate_shards(out/name,2048,width=101)
            m=json.loads((out/name/'manifest.json').read_text())
            item['pipeline_timing']=m['timing'];item['optimization']=m.get('candidate_optimization')
            item['worker_summed_timing']=m.get('candidate_worker_summed_timing')
            files=sorted((out/name).glob('shard_*_data.npy'))+[out/name/'variant_summary.ndjson']
            item['output_sha256']={f.name:digest_file(f) for f in files}
            if name!='baseline_before':
                if item['output_sha256']!=report['steps'][0]['output_sha256']:
                    raise RuntimeError(name+' tensors or complete candidate metadata differ from baseline')
                item['bitwise_equal_to_baseline']=True
            save()
        run('source_audit',[sys.executable,ROOT/'indexed_gam_pipeline/validate_examples.py',out/'indexed_two_workers',
            '--gam',cfg['gam'],'--index',cfg['index'],'--node-sqlite',prior/'all_graph_nodes.sqlite',
            '--gbz',cfg['gbz'],'--gbz-query',cfg['gbz_query'],'--occurrence-cache',out/'indexed_two_workers_gbwt.sqlite',
            '--output',out/'source_audit.json'])
        report['status']='complete'
    except BaseException as e:report.update(status='failed',error=str(e));raise
    finally:save()


if __name__=='__main__':main()
