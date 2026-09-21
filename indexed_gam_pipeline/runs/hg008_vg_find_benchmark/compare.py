import collections,hashlib,json,subprocess,sys,time
from pathlib import Path
sys.path.insert(0,'/scratch/jshen/Github/Pansoma')
from indexed_gam_pipeline.gam_reader import IndexedGam,scan_gam
out=Path(sys.argv[1]);out.mkdir(exist_ok=True)
root=Path('/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/run');cfg=json.loads((root/'config.json').read_text())
with (root/'discovery/target_nodes.txt').open() as f:nodes=[int(next(f)) for _ in range(1024)]
report={'status':'running','vg':'/scratch/jshen/bin/vg_v1.77.0','gam':cfg['gam'],'steps':[],'note':'One Slurm allocation. VG range query followed by exact-node membership filtering. OS cache uncontrolled. VG before and after Python. Python reader reused across two batches.'}
reference={}
def save(): (out/'status.json').write_text(json.dumps(report,indent=2)+'\n')
def fingerprint(records,wanted):
 h=hashlib.sha256();c=collections.Counter();raw=kept=mapq=0
 for a in records:
  raw+=1
  if not any(m.position.node_id in wanted for m in a.path.mapping):continue
  b=a.SerializeToString(deterministic=True);digest=hashlib.sha256(b).hexdigest();c[digest]+=1;h.update(len(b).to_bytes(8,'little'));h.update(b);kept+=1;mapq+=a.mapping_quality>10
 return dict(raw_returned=raw,exact_node_records=kept,passed_mapq=mapq,ordered_sha256=h.hexdigest(),record_multiset=dict(sorted(c.items())))
try:
 report['version']=subprocess.check_output([report['vg'],'version'],text=True,stderr=subprocess.STDOUT);save()
 for mode in ['vg_before','python','vg_after']:
  reader=IndexedGam(cfg['gam'],cfg['index'],8192*1024*1024) if mode=='python' else None
  for bi in range(2):
   wanted=set(nodes[bi*512:(bi+1)*512]);item=dict(mode=mode,batch=bi+1,first_node=min(wanted),last_node=max(wanted),status='running');report['steps'].append(item);save()
   started=time.perf_counter()
   if reader:
    metrics={};result=fingerprint(reader.fetch(wanted,metrics),wanted);item.update(total_seconds=time.perf_counter()-started,metrics=metrics,cache=dict(reader.cache_stats))
   else:
    dest=out/f'{mode}_{bi+1}.gam';cmd=[report['vg'],'find','-l',cfg['gam'],'-o',f'{min(wanted)}:{max(wanted)}'];item['command']=cmd
    with dest.open('wb') as f,(out/f'{mode}_{bi+1}.stderr').open('w') as err:
     subprocess.run(['/usr/bin/time','-v','-o',str(out/f'{mode}_{bi+1}.resources.txt'),*cmd],stdout=f,stderr=err,check=True)
    query=time.perf_counter()-started;t=time.perf_counter();result=fingerprint(scan_gam(dest),wanted)
    item.update(query_seconds=query,output_decode_filter_seconds=time.perf_counter()-t,total_seconds=time.perf_counter()-started,output_bytes=dest.stat().st_size)
   item.update(result);item['status']='complete'
   if bi not in reference:reference[bi]=result
   else:
    item['multiset_equal']=result['record_multiset']==reference[bi]['record_multiset'];item['order_equal']=result['ordered_sha256']==reference[bi]['ordered_sha256']
    if not item['multiset_equal']:raise RuntimeError(f'Record mismatch {mode} batch {bi+1}')
   save();print(json.dumps({k:v for k,v in item.items() if k!='record_multiset'}),flush=True)
  del reader
 report['status']='complete'
except BaseException as e:report.update(status='failed',error=str(e));raise
finally:save()
