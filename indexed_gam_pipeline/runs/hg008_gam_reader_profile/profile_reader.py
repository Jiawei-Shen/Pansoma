import collections, hashlib, inspect, json, sys, textwrap, time
from pathlib import Path
sys.path.insert(0,'/scratch/jshen/Github/Pansoma')
from indexed_gam_pipeline import gam_reader as g

# Instrument a private in-memory copy; production reader is unchanged.
source=textwrap.dedent(inspect.getsource(g.IndexedGam._indexed_group))
source=source.replace('messages = group(stream)', "t=time.perf_counter()\n    messages = group(stream)\n    self.profile['bgzf_group_read_seconds']+=time.perf_counter()-t")
source=source.replace('alignment = decode(raw)', "t=time.perf_counter()\n        alignment = decode(raw)\n        self.profile['protobuf_decode_seconds']+=time.perf_counter()-t\n        t=time.perf_counter()")
source=source.replace('postings[node].append(i)', "postings[node].append(i)\n        self.profile['node_membership_and_postings_seconds']+=time.perf_counter()-t")
source=source.replace('postings = {node: tuple(indices) for node, indices in postings.items()}', "t=time.perf_counter()\n    postings = {node: tuple(indices) for node, indices in postings.items()}\n    self.profile['freeze_postings_seconds']+=time.perf_counter()-t")
source=source.replace('size = (512', 't=time.perf_counter()\n    size = (512')
source=source.replace('if size <= self.cache_limit:', "self.profile['cache_accounting_seconds']+=time.perf_counter()-t\n    self.profile['raw_bytes']+=sum(map(len,messages))\n    self.profile['postings_entries']+=sum(map(len,postings.values()))\n    self.profile['node_keys']+=len(postings)\n    if size <= self.cache_limit:")
ns=dict(vars(g),time=time);exec(compile(source,'<profiled _indexed_group>','exec'),ns)
class Profiled(g.IndexedGam):
    _indexed_group=ns['_indexed_group']

root=Path('/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/run')
cfg=json.loads((root/'config.json').read_text())
with (root/'discovery/target_nodes.txt').open() as f:nodes=[int(next(f)) for _ in range(1024)]
r=Profiled(cfg['gam'],cfg['index'],8192*1024*1024);r.profile=collections.Counter()
report={'source_sha256':hashlib.sha256(Path(g.__file__).read_bytes()).hexdigest(),'batches':[],'note':'Cold application cache; OS/storage cache not controlled. Production source unchanged. Per-record profiling adds overhead.'}
for bi in range(2):
 metrics={};h=hashlib.sha256();count=0;passed=0;t=time.perf_counter();consume=0
 for a in r.fetch(nodes[bi*512:(bi+1)*512],metrics):
  tc=time.perf_counter();b=a.SerializeToString(deterministic=True);h.update(len(b).to_bytes(8,'little'));h.update(b);count+=1;passed+=a.mapping_quality>10;consume+=time.perf_counter()-tc
 item=dict(batch=bi+1,wall_seconds=time.perf_counter()-t,consumer_seconds=consume,metrics=metrics,returned=count,passed_mapq=passed,ordered_record_sha256=h.hexdigest(),cumulative_profile=dict(r.profile),cache=dict(r.cache_stats))
 report['batches'].append(item);print(json.dumps(item),flush=True)
Path(sys.argv[1]).write_text(json.dumps(report,indent=2)+'\n')
