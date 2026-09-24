"""Two independent end-to-end experiments against the saved early-ALT-only baseline."""
import argparse
import contextlib
import json
import os
import signal
from pathlib import Path
import shutil
import subprocess
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from indexed_gam_pipeline.benchmark_candidates import process_memory,digest_file
from indexed_gam_pipeline.full_run import write_json,validate_shards
from indexed_gam_pipeline.tensor_storage import manifest_dtype


def merge_parts(parts,dest,shard_size=2048):
    """Join contiguous partitions, repack shards, stream metadata in candidate order."""
    manifests=[json.loads((p/'manifest.json').read_text()) for p in parts]
    if not manifests:
        raise ValueError('No partitions to merge')
    dtype=manifest_dtype(manifests[0])
    if any(manifest_dtype(m) != dtype or m.get('tensor_storage_version') != manifests[0].get('tensor_storage_version') or m.get('encodings') != manifests[0].get('encodings') for m in manifests):
        raise ValueError('Cannot merge incompatible tensor storage encodings')
    dest.mkdir()
    total=sum(m['tensors'] for m in manifests);size=shard_size
    # Open one source shard at a time: full runs can contain thousands of shards.
    position=0;output=None
    for part in parts:
        for source in sorted(part.glob('shard_*_data.npy')):
            arr=np.load(source,mmap_mode='r',allow_pickle=False)
            if arr.dtype != dtype:
                raise ValueError('Shard dtype differs from manifest')
            if position%size==0 and len(arr)==min(size,total-position):
                os.link(source,dest/f'shard_{position//size:05d}_data.npy')
                position+=len(arr)
            else:
                offset=0
                while offset<len(arr):
                    if position%size==0:
                        output=np.lib.format.open_memmap(dest/f'shard_{position//size:05d}_data.npy',mode='w+',dtype=dtype,
                            shape=(min(size,total-position),7,200,101))
                    count=min(len(arr)-offset,len(output)-position%size)
                    output[position%size:position%size+count]=arr[offset:offset+count]
                    position+=count;offset+=count
                    if position%size==0 or position==total:output.flush();output=None
            del arr
    i=0
    with (dest/'variant_summary.ndjson').open('w') as target:
        for part in parts:
            with (part/'variant_summary.ndjson').open() as stream:
                for line in stream:
                    m=json.loads(line);m['shard_index']=i//size;m['index_within_shard']=i%size
                    target.write(json.dumps(m)+'\n');i+=1
    if i!=total:raise ValueError('Partition count mismatch')
    manifest=dict(manifests[0],tensors=total,shards=(total+size-1)//size,nodes=sum(m['nodes'] for m in manifests),
                  partition_sources=list(map(str,parts)),timing={},merged_contiguous_partitions=True)
    for key in ('filtered_candidates','unsupported_events'):
        if all(key in m for m in manifests):manifest[key]=sum(m[key] for m in manifests)
    for key in ('gam_group_cache','occurrence_performance','candidate_worker_summed_timing','candidate_optimization'):
        manifest.pop(key,None)
    manifest['partition_manifests']=[str(p/'manifest.json') for p in parts]
    manifest['partition_timings']=[m.get('timing',{}) for m in manifests]
    if 'arguments' in manifest:
        manifest['arguments']=dict(manifest['arguments'],output=str(dest),nodes=None)
    write_json(dest/'manifest.json',manifest)
    return validate_shards(dest,size,width=101)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prior-run',type=Path,required=True)
    p.add_argument('--baseline-run',type=Path,required=True)
    p.add_argument('--verified-test',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    prior=args.prior_run.resolve();base=args.baseline_run.resolve();test=args.verified_test.resolve()
    cfg=json.loads((prior/'config.json').read_text())
    baseline=base/'early_alt_only'
    nodes=[int(n) for n in (baseline/'target_nodes.txt').read_text().split()]
    if len(nodes)!=1024 or nodes!=sorted(set(nodes)):raise ValueError('Expected fixed 1024-node baseline')
    report=dict(status='running',baseline=str(baseline),nodes=len(nodes),early_alt_filter=True,node_index_cache_nodes=0,
        candidate_workers=1,debug_rows=True,steps=[],full_job='362252 remains stopped',
        source_comparison='Two treatments only; reuse previous exact-workload early-ALT-only result; no new baseline run')
    def save():write_json(out/'status.json',report)
    reference={f.name:digest_file(f) for f in sorted(baseline.glob('shard_*_data.npy'))+[baseline/'variant_summary.ndjson']}
    report['reference_sha256']=reference;save()
    def command(nodes_path,dest,cache,backend,cache_mb):
        return list(map(str,[sys.executable,ROOT/'indexed_gam_pipeline/run.py','build','--format','candidate-v4',
            '--gam',cfg['gam'],'--index',cfg['index'],'--nodes',nodes_path,'--node-sqlite',prior/'all_graph_nodes.sqlite',
            '--gbz',cfg['gbz'],'--gbz-query',cfg['gbz_query'],'--occurrence-cache',cache,'--output',dest,
            '--gam-reader',backend,'--vg',cfg['vg'],'--batch-nodes','512','--max-node-span','10000',
            '--gam-cache-mb',cache_mb,'--max-batch-segments','20000','--shard-size','2048','--rows','200','--width','101',
            '--min-mapq','10','--min-af','0.05','--min-variants','3','--min-allele-bq','10','--max-indel-len','50',
            '--variant-type','all','--early-alt-filter','--node-index-cache-nodes','0','--workers','1','--debug-rows']))
    try:
        for name,backend,partitions in [('python_two_partitions','python',2),('vg_one_process','vg',1)]:
            case=out/name;case.mkdir();items=[];parts=[]
            for i in range(partitions):
                subset=nodes[i*len(nodes)//partitions:(i+1)*len(nodes)//partitions]
                nf=case/f'nodes_{i}.txt';nf.write_text(''.join(f'{n}\n' for n in subset))
                cache=case/f'gbwt_{i}.sqlite';shutil.copyfile(test/'gbwt_counts.sqlite',cache)
                dest=case/f'part_{i}';parts.append(dest)
                cmd=command(nf,dest,cache,backend,8192//partitions)
                items.append(dict(command=cmd,output=str(dest),nodes=len(subset),first_node=subset[0],last_node=subset[-1]))
            step=dict(name=name,status='running',parts=items);report['steps'].append(step);save()
            start=time.perf_counter();peak_pss=peak_rss=0;sample_seconds=0;procs=[];next_sample=0
            with contextlib.ExitStack() as stack:
                try:
                    for i,item in enumerate(items):
                        log=stack.enter_context((case/f'part_{i}.log').open('w'))
                        proc=subprocess.Popen(['/usr/bin/time','-v','-o',str(case/f'part_{i}.resources.txt'),*item['command']],stdout=log,stderr=subprocess.STDOUT,cwd=ROOT,start_new_session=True)
                        procs.append(proc)
                    while any(x.poll() is None for x in procs):
                        if time.perf_counter()>=next_sample:
                            t=time.perf_counter();samples=[process_memory(x.pid) for x in procs if x.poll() is None]
                            sample_seconds+=time.perf_counter()-t;next_sample=time.perf_counter()+5
                            peak_rss=max(peak_rss,sum(x[0] for x in samples));peak_pss=max(peak_pss,sum(x[1] for x in samples))
                        for proc,item in zip(procs,items):
                            if proc.poll() is not None and 'completion_wall_seconds' not in item:
                                item.update(completion_wall_seconds=time.perf_counter()-start,returncode=proc.returncode)
                        if any(x.poll() not in (None,0) for x in procs):raise RuntimeError(name+' subprocess failed')
                        time.sleep(.2)
                    for proc,item in zip(procs,items):item.update(returncode=proc.returncode)
                    if any(x.returncode for x in procs):raise RuntimeError(name+' subprocess failed')
                finally:
                    for proc in procs:
                        if proc.poll() is None:os.killpg(proc.pid,signal.SIGTERM)
                    for proc in procs:proc.wait()
            step.update(build_wall_seconds=time.perf_counter()-start,sampled_peak_tree_pss_kib=peak_pss,
                        sampled_peak_tree_rss_kib=peak_rss,memory_sampling_seconds=sample_seconds)
            for item,part in zip(items,parts):
                item['validation']=validate_shards(part,2048,width=101)
                manifest=json.loads((part/'manifest.json').read_text())
                item['pipeline_timing']=manifest['timing'];item['gam_cache']=manifest['gam_group_cache']
            if partitions>1:
                t=time.perf_counter();dest=case/'merged';step['validation']=merge_parts(parts,dest)
                step['merge_and_validation_seconds']=time.perf_counter()-t
            else:dest=parts[0];step['validation']=items[0]['validation'];step['merge_and_validation_seconds']=0
            hashes={f.name:digest_file(f) for f in sorted(dest.glob('shard_*_data.npy'))+[dest/'variant_summary.ndjson']}
            step.update(output_sha256=hashes,bitwise_equal_to_baseline=hashes==reference)
            if hashes!=reference:raise RuntimeError(name+' differs from baseline tensors/metadata')
            step['status']='complete';save();print(json.dumps(step),flush=True)
        report['status']='complete'
    except BaseException as e:report.update(status='failed',error=str(e));raise
    finally:save()

if __name__=='__main__':main()
