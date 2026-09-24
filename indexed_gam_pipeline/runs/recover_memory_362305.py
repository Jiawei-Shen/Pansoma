"""Resume only unfinished recovery intervals and publish a nonoverlapping catalog."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import sys

REPO = Path(__file__).resolve().parents[2]
frozen_source = Path(__file__).resolve().parent/'source'
sys.path.insert(0, str(frozen_source if frozen_source.is_dir() else REPO))
from indexed_gam_pipeline.unified_full_run import verify, digest

BASE = Path('/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/run')
PRIOR = BASE/'recovery17_p32_512_nfilter_eoffix_20260922'
ROOT = BASE/'recovery_remaining121_p32_full_cache8g_release_20260922'

def read(p):
    return json.loads(p.read_text())

def check_entries(entries):
    assert len({e['path'] for e in entries}) == len(entries)
    def check(e):
        p = Path(e['path'])/'manifest.json'
        data=p.read_bytes()
        assert hashlib.sha256(data).hexdigest() == e['manifest_sha256'], str(p)
        assert json.loads(data)['status'] == 'complete'
    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(check,entries))

def prepare():
    c = read(PRIOR/'config.json')
    ledger = read(PRIOR/'python/status.json')
    inventory = read(PRIOR/'recovery_inventory.json')
    verify(c)
    done = [t['task'] for t in ledger['tasks'] if t['status']=='complete']
    remaining = [t['task'] for t in ledger['tasks'] if t['status']!='complete']
    assert len(done)==391 and len(remaining)==121
    entries = inventory['preserved_outputs'][:]
    def completed_entry(item):
        i,kind=item
        folder = PRIOR/kind/f'task_{i:04d}'
        data=(folder/'manifest.json').read_bytes(); m=json.loads(data)
        assert m['status']=='complete'
        return dict(kind=kind,path=str(folder),tensors=m['tensors'],
                    manifest_sha256=hashlib.sha256(data).hexdigest(),prior_recovery_task=i)
    with ThreadPoolExecutor(max_workers=8) as pool:
        entries.extend(pool.map(completed_entry,((i,kind) for i in done for kind in ('SNV','INDEL'))))
    check_entries(entries)
    print('Verified 495 original and 391 completed recovery tasks', flush=True)
    ROOT.mkdir(exist_ok=False)
    # Preserve the exact N/EOF behavior of the prior frozen source; replace only
    # memory-management/command code tested in the current workspace.
    shutil.copytree(PRIOR/'source', ROOT/'source', ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    for name in ('build_v2.py','full_reader_comparison.py'):
        shutil.copy2(REPO/'indexed_gam_pipeline'/name, ROOT/'source/indexed_gam_pipeline'/name)
    parts=[]; previous=-1; count=0
    for i in remaining:
        part=dict(c['parts'][i])
        for key in ('nodes_file','preflight_file'):
            src=Path(part[key]); dst=ROOT/src.name
            shutil.copy2(src,dst); part[key]=str(dst)
        for line in Path(part['nodes_file']).open():
            node=int(line); assert node>previous; previous=node; count+=1
        parts.append(part)
    assert count==129916
    c.update(parts=parts,total_nodes=count,replaced_recovery=str(PRIOR),
             source_sha256={str(p):digest(p) for p in (ROOT/'source').rglob('*') if p.is_file()})
    c['params'].update(processes=32,gam_cache_mb_per_process=8192,decode_policy='full',batch_nodes=512)
    (ROOT/'config.json').write_text(json.dumps(c,indent=2))
    inventory.update(status='prepared', recovery_run=str(ROOT), recovery_nodes=count,
                     recovery_tasks=121, preserved_outputs=entries,
                     recovery_task_mapping=remaining, replaced_recovery=str(PRIOR))
    (ROOT/'recovery_inventory.json').write_text(json.dumps(inventory,indent=2))
    shutil.copy2(Path(__file__),ROOT/'finalize.py')
    script=['#!/bin/bash','set -euo pipefail',
            'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1',
            'export PYTHONUNBUFFERED=1','cd '+shlex.quote(str(ROOT/'source')),
            '/usr/bin/time -v -o '+shlex.quote(str(ROOT/'job.resources.txt'))+' '+shlex.quote(sys.executable)+
            ' -m indexed_gam_pipeline.dynamic_tensor_run run --root '+shlex.quote(str(ROOT)),
            shlex.quote(sys.executable)+' '+shlex.quote(str(ROOT/'finalize.py'))+' finalize','']
    (ROOT/'run.sh').write_text('\n'.join(script))
    verify(c)
    print(json.dumps(dict(root=str(ROOT),nodes=count,tasks=121,processes=32)),flush=True)

def finalize():
    c=read(ROOT/'config.json');verify(c)
    inventory=read(ROOT/'recovery_inventory.json'); ledger=read(ROOT/'python/status.json')
    assert read(ROOT/'status.json')['status']=='complete'
    assert ledger['status']=='complete' and ledger['completed_tasks']==len(c['parts'])
    entries=inventory['preserved_outputs'][:];check_entries(entries)
    for i in range(len(c['parts'])):
        for kind in ('SNV','INDEL'):
            p=ROOT/kind/f'task_{i:04d}';m=read(p/'manifest.json')
            assert m['status']=='complete'
            entries.append(dict(kind=kind,path=str(p),tensors=m['tensors'],
                                manifest_sha256=digest(p/'manifest.json'), recovery_task=i))
    assert len({e['path'] for e in entries})==len(entries)
    totals={kind:sum(e['tensors'] for e in entries if e['kind']==kind) for kind in ('SNV','INDEL')}
    result=dict(status='complete',total_nodes=inventory['original_total_nodes'],
                tensors_by_type=totals,tensors=sum(totals.values()),outputs=entries,
                filtering_note=inventory['filtering_note'])
    target=ROOT/'combined_outputs.json';tmp=target.with_suffix('.tmp')
    tmp.write_text(json.dumps(result,indent=2));tmp.replace(target)
    print(json.dumps(totals),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','finalize'])
    args=parser.parse_args()
    (prepare if args.action=='prepare' else finalize)()
