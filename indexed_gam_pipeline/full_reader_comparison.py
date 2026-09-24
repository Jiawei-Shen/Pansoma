"""Guarded whole-dataset reader comparison with two independent contiguous builders."""
import argparse
import contextlib
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from indexed_gam_pipeline.full_run import stamp,write_json,validate_shards
from indexed_gam_pipeline.benchmark_partition_reader import merge_parts
from indexed_gam_pipeline.benchmark_candidates import digest_file


def verify_inputs(config):
    expected = config['graph_index']
    selected = config['config'].get('unified_graph_index', expected['path'])
    if stamp(selected) != expected:
        raise ValueError('Graph index changed or provenance does not match selected index')
    for key,value in config['inputs'].items():
        if stamp(config['config'][key])!=value:raise ValueError('Input changed: '+key)
    for part in config['parts']:
        if digest_file(Path(part['nodes_file']))!=part['sha256']:raise ValueError('Partition changed')
    for path,expected in config['source_sha256'].items():
        if digest_file(ROOT/path)!=expected:raise ValueError('Source differs from prepared snapshot: '+path)


def build_command(config,nodes,dest,cache,backend,cache_mb=4096):
    c=config['config'];prior=Path(config['prior_run'])
    graph_args = (['--graph-index', c['unified_graph_index']] if c.get('unified_graph_index') else
        ['--node-sqlite',prior/'all_graph_nodes.sqlite','--gbz',c['gbz'],'--gbz-query',c['gbz_query'],'--occurrence-cache',cache])
    command = list(map(str,[sys.executable,ROOT/'indexed_gam_pipeline/run.py','build','--format','candidate-v4',
        '--gam',c['gam'],'--index',c['index'],'--nodes',nodes,*graph_args,'--output',dest,
        '--gam-reader',backend,'--vg',c['vg'],'--batch-nodes','512','--max-node-span','10000',
        '--gam-cache-mb',cache_mb,'--max-batch-segments','20000','--shard-size','2048','--rows','200','--width','101',
        '--min-mapq','10','--min-af','0.05','--min-variants','3','--min-allele-bq','10','--max-indel-len','50',
        '--variant-type','all','--early-alt-filter','--node-index-cache-nodes','0','--workers','1']))
    params = config.get('params', {})
    # Written explicitly so frozen snapshots do not depend on CLI defaults.
    command += ['--early-af-filter' if params.get('early_af_filter', True) else '--no-early-af-filter',
                '--candidate-unit', params.get('candidate_unit', 'site'),
                '--max-node-reads', str(params.get('max_node_reads', 800))]
    if params.get('max_tensors'):
        command += ['--max-tensors', str(params['max_tensors'])]
    return command


def run_builders(config,folder,nodefiles,backend):
    folder.mkdir(exist_ok=False);parts=[];commands=[]
    # Separate SQLite caches and files: no writes shared by independent builders.
    for i,nodes in enumerate(nodefiles):
        cache=folder/f'gbwt_{i}.sqlite'
        if not config.get('config', {}).get('unified_graph_index'):
            shutil.copyfile(config['seed_cache'],cache)
        part=folder/f'part_{i}';parts.append(part)
        commands.append(build_command(config,nodes,part,cache,backend,config.get('params',{}).get('gam_cache_mb_per_process',8192)))
    start=time.perf_counter();ledger=dict(status='running',backend=backend,job_id=os.environ.get('SLURM_JOB_ID'),
        commands=commands,processes=[],started_unix=time.time(),memory_sampler=False)
    def save():
        ledger['elapsed_seconds']=time.perf_counter()-start;write_json(folder/'status.json',ledger)
    procs=[]
    def interrupted(signum,frame):raise RuntimeError('Controller signal '+str(signum))
    previous={sig:signal.signal(sig,interrupted) for sig in (signal.SIGTERM,signal.SIGINT)}
    try:
        with contextlib.ExitStack() as stack:
            for i,cmd in enumerate(commands):
                log=stack.enter_context((folder/f'part_{i}.log').open('w'))
                proc=subprocess.Popen(['/usr/bin/time','-v','-o',str(folder/f'part_{i}.resources.txt'),*cmd],
                    stdout=log,stderr=subprocess.STDOUT,cwd=ROOT,start_new_session=True)
                procs.append(proc);ledger['processes'].append(dict(pid=proc.pid,part=i,status='running'))
            save();last_save=0
            while True:
                running=False
                for proc,item in zip(procs,ledger['processes']):
                    code=proc.poll()
                    if code is None:running=True;continue
                    if item['status']=='running':item.update(status='complete' if code==0 else 'failed',returncode=code,wall_seconds=time.perf_counter()-start)
                    if code:raise RuntimeError(f"{backend} partition {item['part']} failed: {code}")
                if not running:break
                if time.perf_counter()-last_save>=30:save();last_save=time.perf_counter()
                time.sleep(.5)
            ledger['build_wall_seconds']=time.perf_counter()-start
            ledger['status']='validating';save();t=time.perf_counter()
            ledger['validation']=[validate_shards(part,2048,width=101) for part in parts]
            ledger['validation_seconds']=time.perf_counter()-t
            ledger['pipeline_timing']=[json.loads((part/'manifest.json').read_text())['timing'] for part in parts]
            ledger['status']='complete';save()
    except BaseException as e:ledger.update(status='failed',error=str(e));save();raise
    finally:
        for proc in procs:
            if proc.poll() is None:
                try:os.killpg(proc.pid,signal.SIGTERM)
                except ProcessLookupError:pass
        for proc in procs:
            try:proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid,signal.SIGKILL);proc.wait()
        for sig,handler in previous.items():signal.signal(sig,handler)
    return parts


def logical_hash(parts):
    tensors=hashlib.sha256();metadata=hashlib.sha256();count=0
    for part in parts:
        for f in sorted(part.glob('shard_*_data.npy')):
            arr=np.load(f,mmap_mode='r',allow_pickle=False)
            for tensor in arr:tensors.update(tensor.tobytes());count+=1
            del arr
        with (part/'variant_summary.ndjson').open() as stream:
            for line in stream:
                m=json.loads(line);m.pop('shard_index');m.pop('index_within_shard')
                metadata.update((json.dumps(m,sort_keys=True,separators=(',',':'))+'\n').encode())
    return dict(tensors=count,tensor_sha256=tensors.hexdigest(),metadata_sha256=metadata.hexdigest())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--mode',choices=['preflight','python','vg','compare'],required=True)
    a=p.parse_args();root=a.root.resolve();config=json.loads((root/'comparison_config.json').read_text())
    verify_inputs(config)
    status=root/(a.mode+'_status.json');start=time.perf_counter()
    ledger=dict(status='running',mode=a.mode,job_id=os.environ.get('SLURM_JOB_ID'),started_unix=time.time(),stages=[])
    def save():ledger['elapsed_seconds']=time.perf_counter()-start;write_json(status,ledger)
    def stage(name,fn):
        x=dict(name=name,status='running');ledger['stages'].append(x);save();t=time.perf_counter()
        try:r=fn();x.update(status='complete',wall_seconds=time.perf_counter()-t);save();return r
        except BaseException as e:x.update(status='failed',error=str(e),wall_seconds=time.perf_counter()-t);save();raise
    try:
        save()
        if a.mode=='preflight':
            files=[Path(x['preflight_file']) for x in config['parts']]
            union=root/'preflight_all_nodes.txt';union.write_text(''.join(f.read_text() for f in files))
            if config.get('preflight_reference'):
                reference=[Path(config['preflight_reference'])]
                stage('validate_reused_serial_reference',lambda:validate_shards(reference[0],2048,width=101))
                if (reference[0]/'target_nodes.txt').read_text()!=union.read_text():raise ValueError('Reference nodes differ')
                if logical_hash(reference)!=config['preflight_reference_hash']:raise ValueError('Reference changed')
            else:
                reference=stage('serial_python_reference',lambda:run_builders(config,root/'preflight_reference',[union],'python'))
            expected=logical_hash(reference)
            if not expected['tensors']:raise ValueError('Preflight produced no tensors')
            ledger['reference']=expected
            for backend in ('python','vg'):
                parts=stage(backend+'_two_partitions',lambda:run_builders(config,root/('preflight_'+backend),files,backend))
                actual=logical_hash(parts);ledger[backend]=actual
                if actual!=expected:raise ValueError('Preflight output mismatch: '+backend)
                dest=root/('preflight_'+backend)/'tensors'
                stage(backend+'_merge',lambda:merge_parts(parts,dest))
                if logical_hash([dest])!=expected:raise ValueError('Merge changed tensors/metadata')
        elif a.mode in ('python','vg'):
            gate=json.loads((root/'preflight_status.json').read_text())
            if gate['status']!='complete':raise ValueError('Preflight is not complete')
            parts=stage('full_build',lambda:run_builders(config,root/a.mode,[Path(x['nodes_file']) for x in config['parts']],a.mode))
            dest=root/a.mode/'tensors'
            stage('merge_and_validate',lambda:merge_parts(parts,dest))
            verify_inputs(config)
            # Retain parts for audit/recovery; aligned full shards are hard-linked.
            ledger['tensors']=json.loads((dest/'manifest.json').read_text())['tensors'];ledger['output']=str(dest)
        else:
            for mode in ('python','vg'):
                if json.loads((root/(mode+'_status.json')).read_text())['status']!='complete':raise ValueError('Build not complete: '+mode)
            hashes={mode:stage('hash_'+mode,lambda mode=mode:logical_hash([root/mode/'tensors'])) for mode in ('python','vg')}
            ledger['hashes']=hashes
            if hashes['python']!=hashes['vg']:raise ValueError('Full outputs differ')
            ledger['equal']=True
        ledger['status']='complete'
    except BaseException as e:ledger.update(status='failed',error=str(e));raise
    finally:save()

if __name__=='__main__':main()
