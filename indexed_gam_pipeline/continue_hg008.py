#!/usr/bin/env python3
"""Resume full HG008 building from verified discovery with a tested source snapshot."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from indexed_gam_pipeline.full_run import stamp,write_json,validate_shards


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior-run',type=Path,required=True)
    parser.add_argument('--verified-test',type=Path,required=True)
    parser.add_argument('--work',type=Path,required=True)
    args=parser.parse_args()
    prior,test,work=(p.resolve() for p in (args.prior_run,args.verified_test,args.work))
    cfg=json.loads((prior/'config.json').read_text())
    verified=json.loads((test/'status.json').read_text())
    if verified['status']!='complete' or not verified['hg008_100_validation']['passed']:
        raise ValueError('HG008 window regression has not passed')
    for filename in ('candidates.py','build_v2.py','run.py','gam_reader.py','gbz_counts.py','validate_examples.py'):
        rel=Path('indexed_gam_pipeline')/filename
        if hashlib.sha256((ROOT/rel).read_bytes()).digest()!=hashlib.sha256((test/'source'/rel).read_bytes()).digest():
            raise ValueError(f'Source differs from tested version: {filename}')
    old=json.loads((prior/'status.json').read_text())
    if not any(s['name']=='discovery' and s['status']=='complete' for s in old['steps']):
        raise ValueError('Prior discovery is incomplete')
    for k,v in old['inputs'].items():
        if stamp(cfg[k])!=v:raise ValueError(f'Input changed since discovery: {k}')
    graph=prior/'all_graph_nodes.sqlite'
    graph_meta=json.loads(Path(str(graph)+'.json').read_text())
    if graph_meta['status']!='complete' or graph_meta['source']!=stamp(cfg['gfa']):
        raise ValueError('Graph sequence index source differs')
    output=Path(cfg['output'])
    building=output/'.building_tensors_window101'
    if building.exists():raise ValueError('Window101 build directory already exists; manual recovery required')
    cache=work/'gbwt_counts.sqlite'
    shutil.copyfile(test/'gbwt_counts.sqlite',cache)
    ledger=dict(status='running',slurm_job_id=os.environ.get('SLURM_JOB_ID'),steps=[],
        started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        inputs=old['inputs'],reused_discovery=str(prior/'discovery'),
        verified_test=str(test),width=101,shard_size=2048,batch_nodes=512,gam_cache_mb=8192,
        reused_stage_times={s['name']:s.get('wall_seconds') for s in old['steps'] if s['status']=='complete'},
        staging_directory=str(building),output_directory=str(output))
    started=time.perf_counter()
    def save():
        ledger['elapsed_seconds']=time.perf_counter()-started
        write_json(work/'status.json',ledger)
    def stage(name,command=None,fn=None):
        entry=dict(name=name,status='running',started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        if command:entry['command']=list(map(str,command))
        ledger['steps'].append(entry);save();t=time.perf_counter()
        try:
            if command:
                with (work/(name+'.log')).open('w') as log:
                    subprocess.run(['/usr/bin/time','-v','-o',str(work/(name+'.resources.txt')),*map(str,command)],
                        stdout=log,stderr=subprocess.STDOUT,check=True,cwd=ROOT)
            else:entry['result']=fn()
            entry['status']='complete'
        except BaseException as e:entry.update(status='failed',error=str(e));raise
        finally:entry['wall_seconds']=time.perf_counter()-t;save()
    try:
        save()
        stage('vg_version',[cfg['vg'],'version'])
        command=[sys.executable,ROOT/'indexed_gam_pipeline/run.py','build','--format','candidate-v4',
            '--gam',cfg['gam'],'--index',cfg['index'],'--nodes',prior/'discovery/target_nodes.txt',
            '--node-sqlite',graph,'--gbz',cfg['gbz'],'--gbz-query',cfg['gbz_query'],
            '--occurrence-cache',cache,'--output',building,'--batch-nodes','512','--max-node-span','10000',
            '--gam-cache-mb','8192','--max-batch-segments','20000','--shard-size','2048',
            '--rows','200','--width','101','--min-mapq','10','--min-af','0.05','--min-variants','3',
            '--min-allele-bq','10','--max-indel-len','50','--variant-type','all']
        stage('tensor_build',command)
        stage('shard_validation',fn=lambda:validate_shards(building,2048,width=101))
        def publish():
            for k,v in ledger['inputs'].items():
                if stamp(cfg[k])!=v:raise ValueError(f'Input changed during build: {k}')
            files=list(building.iterdir())
            for p in files:
                if (output/p.name).exists():raise ValueError(f'Output collision: {p.name}')
            for p in files:p.rename(output/p.name)
            building.rmdir()
            return dict(files=len(files),directory=str(output))
        stage('publish',fn=publish)
        ledger['status']='complete'
    except BaseException as e:ledger.update(status='failed',error=str(e));raise
    finally:save()


if __name__=='__main__':main()
