"""Full Python-reader run: reuse discovery, build unified graph, contiguous processes."""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import signal
import time
import threading

from indexed_gam_pipeline.build_graph_index import build_index
from indexed_gam_pipeline.full_reader_comparison import run_builders
from indexed_gam_pipeline.full_run import stamp, write_json
from indexed_gam_pipeline.graph_index import GraphIndex


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def partition_nodes(source, folder, count, total, preflight_nodes=32):
    """Stream sorted discovery nodes; every node belongs to exactly one interval."""
    if count < 1 or total < count:
        raise ValueError('Need at least one node per process')
    parts = []
    previous = 0
    with Path(source).open() as stream:
        for i in range(count):
            size = total // count + (i < total % count)
            path = Path(folder)/f'nodes_{i}.txt'
            preflight = Path(folder)/f'preflight_nodes_{i}.txt'
            first = None
            with path.open('x') as out, preflight.open('x') as small:
                for j in range(size):
                    line = stream.readline()
                    if not line:
                        raise ValueError('Discovery node count is smaller than report')
                    node = int(line)
                    if node <= previous:
                        raise ValueError('Discovery nodes must be positive, sorted and unique')
                    if first is None:
                        first = node
                    out.write(f'{node}\n')
                    if j < preflight_nodes:
                        small.write(f'{node}\n')
                    previous = node
            parts.append(dict(nodes_file=str(path), preflight_file=str(preflight), nodes=size,
                              first_node=first, last_node=previous, sha256=digest(path),
                              preflight_sha256=digest(preflight)))
        if stream.read().strip():
            raise ValueError('Discovery node count is larger than report')
    return parts


def verify(config):
    for key, expected in config['inputs'].items():
        if stamp(config['config'][key]) != expected:
            raise ValueError('Input changed: '+key)
    for path, expected in config['source_sha256'].items():
        if digest(path) != expected:
            raise ValueError('Source changed: '+path)
    for part in config['parts']:
        for key, hashkey in [('nodes_file','sha256'), ('preflight_file','preflight_sha256')]:
            if digest(part[key]) != part[hashkey]:
                raise ValueError('Partition changed: '+part[key])
    if sum(p['nodes'] for p in config['parts']) != config['total_nodes']:
        raise ValueError('Partition total mismatch')


class MemoryRecorder:
    """Low-overhead process-tree RSS samples; not a PSS/physical-memory estimate."""
    def __init__(self, root, stage, interval=30):
        self.path=Path(root)/'memory.ndjson'
        self.stage=stage
        self.interval=interval
        self.stop=threading.Event()
        self.thread=threading.Thread(target=self.loop,daemon=True)
        self.peak_rss_kib=0
        self.samples=0

    def sample(self):
        pending=[os.getpid()];rows=[];seen=set()
        while pending:
            pid=pending.pop()
            if pid in seen: continue
            seen.add(pid)
            try:
                status={line.split(':',1)[0]:line.split(':',1)[1].strip()
                        for line in Path(f'/proc/{pid}/status').read_text().splitlines() if ':' in line}
                rows.append(dict(pid=pid,name=status.get('Name'),
                                 rss_kib=int(status.get('VmRSS','0 kB').split()[0]),
                                 hwm_kib=int(status.get('VmHWM','0 kB').split()[0])))
                pending.extend(map(int,Path(f'/proc/{pid}/task/{pid}/children').read_text().split()))
            except (FileNotFoundError,ProcessLookupError,PermissionError): pass
        total=sum(r['rss_kib'] for r in rows)
        self.peak_rss_kib=max(self.peak_rss_kib,total);self.samples+=1
        with self.path.open('a') as out:
            out.write(json.dumps(dict(unix=time.time(),stage=self.stage(),sum_rss_kib=total,processes=rows))+'\n')

    def loop(self):
        while not self.stop.wait(self.interval):
            self.sample()

    def __enter__(self):
        self.sample();self.thread.start();return self

    def __exit__(self,*args):
        self.stop.set();self.thread.join();self.sample()


def run(root):
    root = Path(root).resolve()
    config = json.loads((root/'config.json').read_text())
    c = config['config']
    if int(os.environ.get('SLURM_CPUS_PER_TASK', len(config['parts']))) < len(config['parts']):
        raise ValueError('Allocation has fewer CPUs than builder processes')
    start = time.perf_counter()
    ledger = dict(status='running', job_id=os.environ.get('SLURM_JOB_ID'), started_unix=time.time(),
                  processes=len(config['parts']), stages=[], discovery='reused; no GAM rescan')
    def save():
        ledger['wall_seconds'] = time.perf_counter()-start
        write_json(root/'status.json', ledger)
    def stage(name, fn):
        item = dict(name=name, status='running', started_unix=time.time())
        ledger['stages'].append(item); save(); t=time.perf_counter()
        print('START '+name, flush=True)
        try:
            result = fn()
            item['status']='complete'
            return result
        except BaseException as e:
            item.update(status='failed',error=str(e)); raise
        finally:
            item['wall_seconds']=time.perf_counter()-t; save()
            print(f"{item['status'].upper()} {name}: {item['wall_seconds']:.3f} seconds", flush=True)
    def graph():
        result = build_index(c['gbz'], c['unified_graph_index'], c['gbz_index_builder'])
        write_json(root/'graph_index_report.json', result)
    def audit():
        from indexed_gam_pipeline.gbz_counts import GBZCounts
        nodes = sorted({2753,51182} | {p['first_node'] for p in config['parts']})
        with GraphIndex(c['unified_graph_index']) as idx, GBZCounts(c['gbz'], c['gbz_query'], root/'audit_counts.sqlite') as oracle:
            records = idx.get_nodes(nodes)
            actual = {n:r['distinct_path_count'] for n,r in records.items()}
            expected = oracle.get_counts(nodes, {n:r['sequence'] for n,r in records.items()})
            if actual != expected:
                raise ValueError('Unified index differs from GBWT search/locate')
            write_json(root/'graph_audit.json', dict(passed=True,nodes=actual,scope='independent search/locate on known nodes and partition starts'))
    def preflight():
        cfg=copy.deepcopy(config);cfg['params']['max_tensors']=8
        parts=run_builders(cfg,root/'preflight',[Path(p['preflight_file']) for p in cfg['parts']],'python')
        total=sum(json.loads((p/'manifest.json').read_text())['tensors'] for p in parts)
        if total == 0:
            raise ValueError('Preflight generated no tensors')
        ledger['preflight_tensors']=total
    def interrupted(signum, frame):
        raise RuntimeError('Controller signal '+str(signum))
    old={sig:signal.signal(sig,interrupted) for sig in (signal.SIGTERM,signal.SIGINT)}
    memory=MemoryRecorder(root,lambda:ledger['stages'][-1]['name'] if ledger['stages'] else 'startup')
    memory.__enter__()
    try:
        stage('verify_inputs_and_partitions',lambda:verify(config))
        stage('build_unified_graph_index',graph)
        stage('independent_graph_count_audit',audit)
        stage('preflight_20_processes' if len(config['parts'])==20 else 'preflight_processes',preflight)
        stage('full_tensor_build_and_shard_validation',lambda:run_builders(config,root/'python',
              [Path(p['nodes_file']) for p in config['parts']],'python'))
        stage('verify_inputs_after_build',lambda:verify(config))
        ledger['tensors']=sum(json.loads((root/'python'/f'part_{i}'/'manifest.json').read_text())['tensors'] for i in range(len(config['parts'])))
        ledger['status']='complete'
    except BaseException as e:
        ledger.update(status='failed',error=str(e));raise
    finally:
        memory.__exit__()
        ledger['memory']=dict(samples=memory.samples,interval_seconds=30,peak_sampled_sum_rss_kib=memory.peak_rss_kib,
                              metric='sum of process RSS; shared pages may be counted more than once; sampling may miss short peaks')
        save()
        for sig,handler in old.items(): signal.signal(sig,handler)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',required=True)
    run(p.parse_args().root)
