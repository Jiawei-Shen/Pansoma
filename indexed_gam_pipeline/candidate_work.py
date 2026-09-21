"""Bounded node indices and ordered per-node candidate tasks; window code unchanged."""
from collections import Counter, OrderedDict
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
import multiprocessing
import time

from indexed_gam_pipeline.candidates import overlap, make_tensor


class NodeIndexCache:
    """FIFO; references only, batch-scoped so obsolete decoded reads never accumulate.

    Byte accounting is an O(1)-per-entry conservative estimate of index containers,
    not a measurement of the independently owned reads or total process RSS.
    """
    def __init__(self, max_nodes=1000, max_bytes=64*1024*1024):
        if not 0 <= max_nodes <= 1000 or max_bytes < 0:
            raise ValueError('Node index limit must be in [0,1000]; bytes must be nonnegative')
        self.max_nodes,self.max_bytes=max_nodes,max_bytes
        self.entries=OrderedDict();self.bytes=0
        self.stats=dict(hits=0,misses=0,evictions=0,peak_nodes=0,peak_estimated_bytes=0)

    def get(self,node,reads):
        if not self.max_nodes or not self.max_bytes:return None
        if node in self.entries:
            self.stats['hits']+=1
            return self.entries[node][0]  # FIFO: do not move on access.
        self.stats['misses']+=1
        result={};size=256
        for read in reads:
            visits=tuple(v for v in read.visits if v.node==node)
            evidence={}
            for obs in read.observations:
                if obs.candidate.node==node:
                    key=(obs.visit,obs.candidate)
                    evidence[key]=max(evidence.get(key,float('-inf')),obs.quality)
            result[id(read)]=(visits,evidence)
            size+=256+16*len(visits)+160*len(evidence)
            if size>self.max_bytes:
                # Do not retain an oversized node, or grow a giant temporary index.
                return None
        while self.entries and (len(self.entries)>=self.max_nodes or self.bytes+size>self.max_bytes):
            _,(_,removed)=self.entries.popitem(last=False)
            self.bytes-=removed;self.stats['evictions']+=1
        self.entries[node]=(result,size);self.bytes+=size
        self.stats['peak_nodes']=max(self.stats['peak_nodes'],len(self.entries))
        self.stats['peak_estimated_bytes']=max(self.stats['peak_estimated_bytes'],self.bytes)
        return result


def alt_support_bounds(reads,candidates,min_bq):
    """One vote per record/candidate, not per visit; identical GAM records retained."""
    counts=Counter()
    for read in reads:
        counts.update({o.candidate for o in read.observations
                       if o.quality>=min_bq and o.candidate in candidates})
    return counts


_STATE=None
_CACHE=None


def process_node(task):
    node,candidates=task
    by_node,args,walk_counts=_STATE
    reads=by_node[node]
    results=[];timing=Counter()
    t=time.perf_counter();index=_CACHE.get(node,reads)
    timing['node_index_seconds']+=time.perf_counter()-t
    for candidate in candidates:
        t=time.perf_counter();eligible=[]
        for read in reads:
            hit=overlap(read,candidate,args.min_allele_bq,
                        node_index=index.get(id(read)) if index is not None else None)
            if hit is not None:eligible.append((read,*hit))
        counts=Counter(s for _,s,_ in eligible)
        af=counts['alt']/len(eligible) if eligible else 0
        reasons=[]
        if counts['alt']<args.min_variants:reasons.append('min_variants')
        if af<args.min_af:reasons.append('min_af')
        if args.variant_type!='all' and (candidate.kind=='SNP')!=(args.variant_type=='snp'):
            reasons.append('variant_type')
        timing['support_count_seconds']+=time.perf_counter()-t
        if reasons:
            results.append((None,dict(candidate.metadata(),reasons=reasons,coverage=len(eligible),
                alt_count=counts['alt'],ref_count=counts['ref'],other_count=counts['other'],af=af)))
        else:
            t=time.perf_counter()
            tensor,meta=make_tensor(candidate,eligible,args.rows,args.width,args.debug_rows,node_walk_counts=walk_counts)
            timing['window_sort_encode_seconds']+=time.perf_counter()-t
            results.append((tensor,meta))
    return results,dict(timing),dict(_CACHE.stats)


@contextmanager
def candidate_results(tasks,by_node,args,walk_counts):
    """Share decoded inputs via Linux fork; bounded submissions, deterministic output.

    Workers never access GAM, SQLite or GBZ handles. Only the parent reads inputs
    and writes shards. At most `workers` node-task results are in flight. Process
    pools are batch-scoped to release obsolete decoded data and index references.
    """
    global _STATE,_CACHE
    workers=getattr(args,'workers',1)
    limit=getattr(args,'node_index_cache_nodes',0)
    budget=getattr(args,'node_index_cache_mb',64)*1024*1024
    if workers<1 or workers>max(1,limit) and limit:
        raise ValueError('workers must be positive and not exceed the total node cache limit')
    _STATE=(by_node,args,walk_counts)
    _CACHE=NodeIndexCache(limit//workers,budget//workers)
    pool=None
    try:
        if workers==1:
            yield map(process_node,tasks)
        else:
            if 'fork' not in multiprocessing.get_all_start_methods():
                raise ValueError('Candidate workers currently require Linux fork')
            pool=ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('fork'))
            iterator=iter(tasks);pending=[]
            for _ in range(workers):
                task=next(iterator,None)
                if task is not None:pending.append(pool.submit(process_node,task))
            def ordered():
                while pending:
                    future=pending.pop(0)
                    result=future.result()
                    task=next(iterator,None)
                    if task is not None:pending.append(pool.submit(process_node,task))
                    yield result
            yield ordered()
    finally:
        if pool is not None:pool.shutdown(wait=True,cancel_futures=True)
        _STATE=None;_CACHE=None
