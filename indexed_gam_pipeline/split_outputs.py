"""Two bounded shard writers fed by one shared candidate decoding pass."""
from contextlib import ExitStack
from copy import deepcopy
from pathlib import Path
import json
import numpy as np


def split_enabled(args):
    values = [getattr(args, key, None) for key in
              ('snv_output', 'indel_output', 'snv_min_af', 'indel_min_af')]
    if all(v is None for v in values):
        return False
    if any(v is None for v in values):
        raise ValueError('Split generation requires both output directories and both AF thresholds')
    if not all(0 <= v <= 1 for v in values[2:]):
        raise ValueError('Split AF thresholds must be in [0,1]')
    if args.variant_type != 'all' or args.format != 'candidate-v4':
        raise ValueError('Split generation requires candidate-v4 and variant-type all')
    paths = [Path(p).resolve() for p in (args.output, *values[:2])]
    if len(set(paths)) != 3 or any(a in b.parents for a in paths for b in paths if a != b):
        raise ValueError('Shared, SNV and INDEL output directories must be distinct and non-nested')
    if any(p.exists() for p in paths):
        raise ValueError('Split output directories must not already exist')
    return True


class SplitOutputs:
    def __init__(self, args, manifest, timings):
        self.args, self.shared, self.timings = args, manifest, timings
        self.stack = ExitStack()
        self.outputs = {}

    def __enter__(self):
        from indexed_gam_pipeline.run import new_output, write_json
        try:
            for name, path, kind, af in [('SNV', self.args.snv_output, 'snp', self.args.snv_min_af),
                                         ('INDEL', self.args.indel_output, 'indel', self.args.indel_min_af)]:
                out = new_output(path)
                m = deepcopy(self.shared)
                m['parameters'].update(variant_type=kind, min_af=af)
                m['arguments'] = dict(m['arguments'], output=str(out), variant_type=kind, min_af=af)
                m['shared_build_directory'] = str(Path(self.args.output).resolve())
                m['timing_scope'] = 'Shared read/decode/build pass; do not sum across output types'
                m['output_layout'] = 'separate-type-shards'
                self.outputs[name] = dict(path=out, manifest=m, tensors=[], metadata=[],
                    streams={key: self.stack.enter_context((out/file).open('w')) for key,file in
                             [('summary','variant_summary.ndjson'),('filtered','filtered_candidates.ndjson'),
                              ('unsupported','unsupported_events.ndjson')]})
                write_json(out/'manifest.json',m)
                write_json(out/'run_report.json',m)
            self.shared['output_layout'] = 'split-coordinator-no-tensor-shards'
            self.shared['variant_outputs'] = {k:str(v['path'].resolve()) for k,v in self.outputs.items()}
            return self
        except BaseException:
            self.stack.close()
            raise

    def __exit__(self, *exc):
        return self.stack.__exit__(*exc)

    @staticmethod
    def kind(meta):
        return 'SNV' if meta.get('event_type') == 'SNP' else 'INDEL' if meta.get('event_type') in ('INS','DEL') else None

    @property
    def buffered(self):
        return sum(len(v['tensors']) for v in self.outputs.values())

    def record(self, stream, meta):
        name = self.kind(meta)
        if name is None:
            return  # Complex unsupported edits remain in the shared audit only.
        out = self.outputs[name]
        out['streams'][stream].write(json.dumps(meta)+'\n')
        key = 'filtered_candidates' if stream == 'filtered' else 'unsupported_events'
        out['manifest'][key] += 1
        if stream == 'filtered' and meta.get('coverage_not_evaluated'):
            out['manifest']['candidate_optimization']['early_rejected'] += 1
        if stream == 'filtered' and meta.get('support_not_evaluated'):
            out['manifest']['candidate_optimization']['early_af_rejected'] += 1

    def append(self, tensor, meta):
        name = self.kind(meta)
        out = self.outputs[name]
        out['tensors'].append(tensor); out['metadata'].append(meta)
        if len(out['tensors']) >= self.args.shard_size:
            self.flush(name)

    def flush(self, name=None):
        from indexed_gam_pipeline.run import write_json
        for key in ([name] if name else self.outputs):
            out = self.outputs[key]; m = out['manifest']
            if not out['tensors']:
                continue
            shard = m['shards']
            np.save(out['path']/f'shard_{shard:05d}_data.npy', np.stack(out['tensors']))
            for i,meta in enumerate(out['metadata']):
                out['streams']['summary'].write(json.dumps(dict(meta, schema_version=m['schema_version'],
                    channels=m['channels'], parameters=m['parameters'], shard_index=shard,index_within_shard=i))+'\n')
            m['shards'] += 1; m['tensors'] += len(out['tensors'])
            self.shared['shards'] += 1; self.shared['tensors'] += len(out['tensors'])
            out['tensors'].clear(); out['metadata'].clear(); out['streams']['summary'].flush()
            m['timing'] = dict(self.timings)
            write_json(out['path']/'manifest.json',m); write_json(out['path']/'run_report.json',m)
        self.shared['tensors_by_type'] = {k:v['manifest']['tensors'] for k,v in self.outputs.items()}
        self.shared['timing'] = dict(self.timings)
        write_json(Path(self.args.output)/'manifest.json', self.shared)
        write_json(Path(self.args.output)/'run_report.json', self.shared)

    def finish(self, shared):
        from indexed_gam_pipeline.run import write_json
        for out in self.outputs.values():
            m = out['manifest']
            for key in ('timing','gam_group_cache','occurrence_performance','occurrence_lookup','candidate_worker_summed_timing'):
                if key in shared:
                    m[key] = deepcopy(shared[key])
            m['status'] = 'complete'
            write_json(out['path']/'manifest.json',m); write_json(out['path']/'run_report.json',m)
