"""VG-backed complete-alignment retrieval with exact target-node filtering."""
import subprocess
import tempfile
import time
from pathlib import Path

from indexed_gam_pipeline.gam_reader import scan_gam


class VgGam:
    def __init__(self, gam, index=None, vg='/scratch/jshen/bin/vg_v1.77.0'):
        self.gam = Path(gam).resolve()
        self.index = Path(index or str(self.gam)+'.gai').resolve()
        # vg find resolves its index beside the GAM; never silently ignore --index.
        if self.index != Path(str(self.gam)+'.gai').resolve():
            raise ValueError('vg reader requires the adjacent GAM.gai index')
        if not self.index.is_file():
            raise ValueError('vg reader requires an existing GAM.gai index')
        self.vg = str(Path(vg).resolve())
        self.version = 'vg-find-adjacent-gai'
        self._stamps = self._fingerprints()
        self.cache_stats = dict(backend='vg', queries=0, query_seconds=0., output_bytes=0,
                                group_hits=0, group_misses=0, limit_bytes=0)

    def _fingerprints(self):
        return [(p.stat().st_size,p.stat().st_mtime_ns) for p in (self.gam,self.index)]

    def fetch(self,nodes,metrics=None):
        wanted=set(nodes)
        if metrics is None:metrics={}
        metrics.update(backend='vg',decoded_alignments=0,returned_alignments=0)
        if not wanted:return
        if self._fingerprints()!=self._stamps:
            raise ValueError('GAM/index changed after opening vg reader')
        # File-backed transport bounds RAM, captures errors before yielding records,
        # and cleans up even when a consumer exits early. TMPDIR chooses storage.
        with tempfile.TemporaryDirectory(prefix='pansoma-vg-') as directory:
            output=Path(directory)/'query.gam'
            command=[self.vg,'find','-l',str(self.gam),'-o',f'{min(wanted)}:{max(wanted)}']
            t=time.perf_counter()
            with output.open('wb') as stream, tempfile.TemporaryFile() as error:
                proc=subprocess.run(command,stdout=stream,stderr=error)
                if proc.returncode:
                    error.seek(0,2);size=error.tell();error.seek(max(0,size-8192))
                    raise RuntimeError(f'vg find failed ({proc.returncode}): '+error.read().decode(errors='replace'))
            elapsed=time.perf_counter()-t;size=output.stat().st_size
            self.cache_stats['queries']+=1
            self.cache_stats['query_seconds']+=elapsed
            self.cache_stats['output_bytes']+=size
            metrics.update(vg_query_seconds=elapsed,vg_output_bytes=size)
            for alignment in scan_gam(output):
                metrics['decoded_alignments']+=1
                if any(m.position.node_id in wanted for m in alignment.path.mapping):
                    metrics['returned_alignments']+=1
                    yield alignment
            if self._fingerprints()!=self._stamps:
                raise ValueError('GAM/index changed during vg query')
