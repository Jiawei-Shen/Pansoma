#!/usr/bin/env python3
"""Timed HG008 regression: old batch 111 failure, new windows, then 100 examples."""
import argparse
import datetime
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from indexed_gam_pipeline.full_run import write_json, validate_shards, stamp
from indexed_gam_pipeline.run import batches


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior-run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    prior, out = args.prior_run.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((prior/'config.json').read_text())
    with (prior/'discovery/target_nodes.txt').open() as stream:
        prefix = [int(line) for line in itertools.islice(stream, 8192)]
    failed_nodes = next(b for i,b in enumerate(batches(prefix,32,10000),1) if i == 111)
    (out/'batch111_nodes.txt').write_text(''.join(f'{n}\n' for n in failed_nodes))
    (out/'early_nodes.txt').write_text(''.join(f'{n}\n' for n in prefix[:1024]))
    shutil.copyfile(prior/'gbwt_counts.sqlite',out/'gbwt_counts.sqlite')
    py = sys.executable
    ledger = dict(status='running',slurm_job_id=os.environ.get('SLURM_JOB_ID'),
        window_encoding_version='anchor-centered-columns-v1',shape=[7,200,101],
        batch111_nodes=failed_nodes,steps=[],inputs={k:stamp(cfg[k]) for k in ('gam','index','gbz')})
    started = time.perf_counter()
    def save():
        ledger['elapsed_seconds']=time.perf_counter()-started
        write_json(out/'status.json',ledger)
    def stage(name,command,expected=0):
        entry=dict(name=name,status='running',command=list(map(str,command)),
                   started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        ledger['steps'].append(entry);save();t=time.perf_counter()
        with (out/(name+'.log')).open('w') as log:
            result=subprocess.run(['/usr/bin/time','-v','-o',str(out/(name+'.resources.txt')),
                *map(str,command)],stdout=log,stderr=subprocess.STDOUT,cwd=ROOT)
        entry.update(seconds=time.perf_counter()-t,returncode=result.returncode,
                     status='complete' if result.returncode==expected else 'failed')
        save()
        if result.returncode!=expected:raise RuntimeError(f'{name} exited {result.returncode}')
    def build(dest,nodes,old=False):
        source=prior/'source' if old else ROOT
        cmd=[py,source/'indexed_gam_pipeline/run.py','build','--format','candidate-v4',
            '--gam',cfg['gam'],'--index',cfg['index'],'--nodes',nodes,
            '--node-sqlite',prior/'all_graph_nodes.sqlite','--gbz',cfg['gbz'],
            '--gbz-query',cfg['gbz_query'],'--occurrence-cache',out/'gbwt_counts.sqlite',
            '--output',dest,'--batch-nodes','32' if old else '512','--max-node-span','10000',
            '--gam-cache-mb','8192','--max-batch-segments','20000','--shard-size','2048',
            '--rows','200','--width','100' if old else '101','--min-mapq','10','--min-af','0.05',
            '--min-variants','3','--min-allele-bq','10','--max-indel-len','50','--variant-type','all']
        if not old:cmd.append('--debug-rows')
        return cmd
    def audit(dest,label):
        stage(label,[py,ROOT/'indexed_gam_pipeline/validate_examples.py',dest,
            '--gam',cfg['gam'],'--index',cfg['index'],'--node-sqlite',prior/'all_graph_nodes.sqlite',
            '--gbz',cfg['gbz'],'--gbz-query',cfg['gbz_query'],
            '--occurrence-cache',out/'gbwt_counts.sqlite','--output',out/(label+'.json')])
        return validate_shards(dest,2048,width=101)
    try:
        save()
        probe=out/'old_probe.py'
        probe.write_text('''import json,sys
from pathlib import Path
source=Path(sys.argv.pop(1)); report=Path(sys.argv.pop(1))
sys.path.insert(0,str(source))
from indexed_gam_pipeline import build_v2,run
from indexed_gam_pipeline.candidates import oriented_columns,split_columns
original=build_v2.make_tensor
def traced(candidate,eligible,*args,**kwargs):
 try:return original(candidate,eligible,*args,**kwargs)
 except ValueError as e:
  if 'Tensor width cannot preserve' not in str(e):raise
  reads=[]
  for r,s,v in eligible:
   center=split_columns(oriented_columns(r,v,1),v,candidate)[1]
   reads.append(dict(read_name=r.name,record_sha256=r.digest,support=s,
      central_insertion_bp=sum(c.boundary for c in center)))
  report.write_text(json.dumps(dict(candidate=candidate.metadata(),error=str(e),
    reads=sorted(reads,key=lambda r:-r['central_insertion_bp'])),indent=2)+'\\n')
  raise
build_v2.make_tensor=traced
run.main()
''')
        cmd=build(out/'old_failed_batch',out/'batch111_nodes.txt',old=True)
        stage('reproduce_old_failure',[py,probe,prior/'source',out/'original_failure.json',*cmd[2:]],expected=1)
        failure=json.loads((out/'original_failure.json').read_text())
        if 'Tensor width cannot preserve' not in failure['error']:raise RuntimeError('Wrong reproduced error')
        ledger['original_failure_candidate']=failure['candidate'];save()
        stage('fixed_batch111_build',build(out/'fixed_batch111',out/'batch111_nodes.txt'))
        ledger['fixed_batch111_validation']=audit(out/'fixed_batch111','fixed_batch111_audit')
        metas=[json.loads(line) for line in (out/'fixed_batch111/variant_summary.ndjson').read_text().splitlines()]
        if not any(m['candidate_id']==failure['candidate']['candidate_id'] for m in metas):
            raise RuntimeError('Original failure candidate absent from fixed output')
        stage('hg008_100_build',build(out/'hg008_100',out/'early_nodes.txt')+['--max-tensors','100'])
        ledger['hg008_100_validation']=audit(out/'hg008_100','hg008_100_audit')
        if ledger['hg008_100_validation']['tensors']!=100:raise RuntimeError('Expected 100 tensors')
        if any(stamp(cfg[k])!=v for k,v in ledger['inputs'].items()):raise RuntimeError('Inputs changed')
        ledger['status']='complete'
    except BaseException as e:
        ledger.update(status='failed',error=str(e));raise
    finally:save()


if __name__=='__main__':main()
