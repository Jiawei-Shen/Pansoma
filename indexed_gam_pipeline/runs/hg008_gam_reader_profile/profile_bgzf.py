import json,sys,time
from pathlib import Path
sys.path.insert(0,'/scratch/jshen/Github/Pansoma')
from indexed_gam_pipeline import gam_reader as g
root=Path('/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/run');cfg=json.loads((root/'config.json').read_text())
with (root/'discovery/target_nodes.txt').open() as f:nodes=[int(next(f)) for _ in range(512)]
r=g.IndexedGam(cfg['gam'],cfg['index'],0);runs=r.ranges(nodes)
report={'runs':len(runs),'compressed_span_bytes':sum((b>>16)-(a>>16) for a,b in runs),'note':'Same ranges; OS cache uncontrolled. CPU versus wall helps distinguish execution from wait/descheduling, not disk latency alone.'}
t=time.perf_counter();cpu=time.process_time();records=groups=payload=0
with g.pysam.BGZFile(cfg['gam'],'rb') as stream:
 for start,end in runs:
  stream.seek(start)
  while stream.tell()<end:
   messages=g.group(stream)
   if messages is None:raise RuntimeError('EOF')
   groups+=1;records+=len(messages);payload+=sum(map(len,messages))
  if stream.tell()!=end:raise RuntimeError('boundary mismatch')
report['bgzf_only']=dict(wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-cpu,groups=groups,records=records,payload_bytes=payload)
print(json.dumps(report),flush=True)
t=time.perf_counter();cpu=time.process_time();n=0
with open(cfg['gam'],'rb',buffering=0) as f:
 for start,end in runs:
  f.seek(start>>16);left=(end>>16)-(start>>16)
  while left:
   b=f.read(min(left,4*1024*1024))
   if not b:raise RuntimeError('EOF')
   left-=len(b);n+=len(b)
report['compressed_read_after_bgzf']=dict(wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-cpu,bytes=n)
Path(sys.argv[1]).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
