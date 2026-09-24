"""Bounded node indices and ordered per-node candidate tasks; window code unchanged."""
from collections import Counter, OrderedDict, defaultdict
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
import multiprocessing
import time

from indexed_gam_pipeline.candidates import Candidate, overlap, make_tensor


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


def exact_coverage(candidates,by_node):
    """Exact overlap() coverage without materializing columns.

    overlap() returns a hit for a record iff one of its non-empty visits to the
    candidate node covers the candidate (closed interval for INS, half-open
    overlap otherwise); each record counts once. That test needs only visit
    intervals, so coverage here equals len(eligible) in process_node when both
    see the same records.
    """
    import numpy as np
    by_cand_node=defaultdict(list)
    for c in candidates:by_cand_node[c.node].append(c)
    result={}
    for node,items in by_cand_node.items():
        rows=[(ri,v.start,v.end) for ri,read in enumerate(by_node.get(node,()))
              for v in read.visits if v.node==node and v.first!=v.last]
        if not rows:
            result.update((c,0) for c in items);continue
        rid,vs,ve=(np.fromiter((r[k] for r in rows),dtype=np.int64,count=len(rows)) for k in range(3))
        for c in items:
            mask=(vs<=c.start)&(ve>=c.start) if c.kind=='INS' else (vs<c.end)&(ve>c.start)
            result[c]=int(np.unique(rid[mask]).size)
    return result


def af_threshold(candidate,args):
    threshold=getattr(args,'snv_min_af' if candidate.kind=='SNP' else 'indel_min_af',None)
    return args.min_af if threshold is None else threshold


_STATE=None
_CACHE=None


DEFAULT_MAX_NODE_READS=800
SITE_UNIT='site-v1'


def site_key(candidate):
    """SNVs and indels are separate sites (separate outputs and models); INS and DEL
    anchored at the same forward node position share one indel site."""
    return candidate.node,candidate.start,'SNV' if candidate.kind=='SNP' else 'INDEL'


def site_id(candidate):
    return '%d:%d:%s'%site_key(candidate)


def group_sites(candidates):
    """Sorted candidates -> sorted tuples of alleles, one tuple per site."""
    sites=defaultdict(list)
    for c in sorted(candidates):sites[site_key(c)].append(c)
    return [tuple(sites[k]) for k in sorted(sites)]


def node_reads(by_node,node,args):
    reads=by_node[node]
    cap=getattr(args,'max_node_reads',DEFAULT_MAX_NODE_READS)
    if cap and len(reads)>cap:
        # Deterministic subsample: order by record digest so the same records are
        # kept for every candidate on this node and on every rerun.
        reads=sorted(reads,key=lambda r:r.digest)[:cap]
    return reads


def evaluate(candidate,reads,index,args):
    eligible=[]
    for read in reads:
        hit=overlap(read,candidate,args.min_allele_bq,
                    node_index=index.get(id(read)) if index is not None else None)
        if hit is not None:eligible.append((read,*hit))
    counts=Counter(s for _,s,_ in eligible)
    af=counts['alt']/len(eligible) if eligible else 0
    reasons=[]
    if counts['alt']<args.min_variants:reasons.append('min_variants')
    if af<af_threshold(candidate,args):reasons.append('min_af')
    if args.variant_type!='all' and (candidate.kind=='SNP')!=(args.variant_type=='snp'):
        reasons.append('variant_type')
    summary=dict(candidate.metadata(),coverage=len(eligible),alt_count=counts['alt'],
                 ref_count=counts['ref'],other_count=counts['other'],af=af)
    return eligible,summary,reasons


ALLELE_FIELDS=('candidate_id','start','end','ref','alt','event_type','event_length',
               'coverage','alt_count','ref_count','other_count','af')


def process_node(task):
    """Candidate units on one node. Allele mode: each unit is one candidate (legacy
    output). Site mode: each unit holds every prefiltered allele of one site; every
    allele is evaluated exactly as in allele mode, and one tensor is emitted for the
    site's representative (highest ALT count among passing alleles; ties by allele
    order), so the tensor equals that allele's legacy tensor."""
    node,units=task
    by_node,args,walk_counts=_STATE
    reads=node_reads(by_node,node,args)
    sites=getattr(args,'candidate_unit','site')=='site'
    results=[];timing=Counter()
    t=time.perf_counter();index=_CACHE.get(node,reads)
    timing['node_index_seconds']+=time.perf_counter()-t
    for unit in units:
        unit=(unit,) if isinstance(unit,Candidate) else tuple(unit)
        t=time.perf_counter()
        evaluated=[(c,*evaluate(c,reads,index,args)) for c in unit]
        timing['support_count_seconds']+=time.perf_counter()-t
        passing=sorted((e for e in evaluated if not e[3]),key=lambda e:(-e[2]['alt_count'],e[0]))
        extra=dict(site_id=site_id(unit[0])) if sites else {}
        for candidate,eligible,summary,reasons in evaluated:
            if reasons:  # legacy key order: metadata, reasons, counts
                results.append((None,dict(candidate.metadata(),reasons=reasons,
                    **{k:summary[k] for k in ('coverage','alt_count','ref_count','other_count','af')},**extra)))
        if not passing:continue
        if not sites and len(passing)>1:raise AssertionError('allele-mode unit with several passing alleles')
        candidate,eligible,summary,_=passing[0]
        t=time.perf_counter()
        tensor,meta=make_tensor(candidate,eligible,args.rows,args.width,args.debug_rows,node_walk_counts=walk_counts)
        timing['window_sort_encode_seconds']+=time.perf_counter()-t
        if sites:
            alleles=[{k:e[2][k] for k in ALLELE_FIELDS} for e in passing]
            meta.update(sample_unit=SITE_UNIT,site_id=extra['site_id'],alleles=alleles,
                        allele_count=len(alleles),second_allele_af=alleles[1]['af'] if len(alleles)>1 else 0.0)
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
