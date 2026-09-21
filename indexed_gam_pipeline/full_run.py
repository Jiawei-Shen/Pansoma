#!/usr/bin/env python3
"""Timed, fail-fast full-data orchestration; publish 2048-example NPY shards."""
import argparse
from collections import Counter
import datetime
import json
import os
import resource
from pathlib import Path
import sqlite3
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def stamp(path):
    path = Path(path).resolve()
    s = path.stat()
    return dict(path=str(path), size=s.st_size, mtime_ns=s.st_mtime_ns)


def graph_index(gfa, output):
    """All S records, including alternate nodes; bounded transactions and RAM."""
    initial = stamp(gfa)
    output = Path(output)
    if output.exists():
        raise ValueError(f'Graph index already exists: {output}')
    partial = output.with_suffix('.building.sqlite')
    count = 0
    with sqlite3.connect(partial) as db:
        db.execute('PRAGMA journal_mode=OFF')  # unpublished, rebuilt on any failure
        db.execute('PRAGMA synchronous=OFF')
        db.execute('PRAGMA cache_size=-131072')
        db.execute('CREATE TABLE nodes (node_id INTEGER PRIMARY KEY, seq TEXT NOT NULL)')
        batch = []
        with open(gfa, 'rb') as stream:
            for line in stream:
                if not line.startswith(b'S\t'):
                    continue
                fields = line.split(b'\t', 3)
                node, sequence = int(fields[1]), fields[2].strip().decode('ascii').upper()
                if node <= 0 or sequence == '*' or not sequence:
                    raise ValueError(f'Missing/invalid graph sequence for node {node}')
                batch.append((node, sequence))
                if len(batch) == 100000:
                    db.executemany('INSERT INTO nodes VALUES (?,?)', batch)
                    db.commit(); count += len(batch); batch.clear()
                    if count % 1000000 == 0:
                        print(f'Graph index: {count:,} nodes', flush=True)
        db.executemany('INSERT INTO nodes VALUES (?,?)', batch)
        db.commit(); count += len(batch)
        if not count or stamp(gfa) != initial:
            raise ValueError('Empty graph or GFA changed during indexing')
    partial.replace(output)
    write_json(str(output) + '.json', dict(source=initial, nodes=count, status='complete'))
    return count


def validate_shards(folder, shard_size, width=101):
    folder = Path(folder)
    manifest = json.loads((folder / 'manifest.json').read_text())
    if manifest['status'] != 'complete' or manifest['tensor_format_version'] != 'indexed-gam-candidate-v4':
        raise ValueError('Incomplete or unexpected tensor format')
    expected = manifest['shards']
    if len(list(folder.glob('shard_*_data.npy'))) != expected:
        raise ValueError('Shard file count mismatch')
    sizes = []
    for i in range(expected):
        path = folder / f'shard_{i:05d}_data.npy'
        x = np.load(path, mmap_mode='r', allow_pickle=False)
        if x.ndim != 4 or x.shape[1:] != (7, 200, width) or x.dtype != np.int32:
            raise ValueError(f'Unexpected shard shape/dtype: {path}')
        if not 0 < len(x) <= shard_size or (i < expected-1 and len(x) != shard_size):
            raise ValueError(f'Invalid shard length: {path}')
        if path.stat().st_size != x.offset + x.nbytes:
            raise ValueError(f'Truncated or overlong shard: {path}')
        sizes.append(len(x))
        del x
    seen = [0] * expected
    events = Counter()
    with (folder / 'variant_summary.ndjson').open() as stream:
        for line in stream:
            m = json.loads(line)
            shard = m['shard_index']
            if not 0 <= shard < expected or m['index_within_shard'] != seen[shard]:
                raise ValueError('Summary shard/index ordering mismatch')
            seen[shard] += 1
            if m['coverage'] != sum(m[k] for k in ('alt_count','ref_count','other_count')):
                raise ValueError('Candidate coverage mismatch')
            if m['af'] != m['alt_count']/m['coverage'] or m['selected_alignments'] != min(200,m['coverage']):
                raise ValueError('Candidate AF/row count mismatch')
            if m['tensor_format_version'] != 'indexed-gam-candidate-v4':
                raise ValueError('Summary format mismatch')
            events[m['event_type']] += 1
    if seen != sizes or sum(sizes) != manifest['tensors']:
        raise ValueError('Summary, shard, and manifest tensor counts differ')
    report = dict(passed=True, tensors=sum(sizes), shards=expected, shard_size=shard_size,
        last_shard_tensors=sizes[-1] if sizes else 0, shape_per_tensor=[7,200,width],
        dtype='int32', event_types=dict(events),
        scope='Every NPY header/file size and every summary record; independent source audit performed on the preflight tensors')
    write_json(folder / 'validation_report.json', report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    args = parser.parse_args()
    cfg = json.loads(Path(args.config).read_text())
    output = Path(cfg['output']).resolve()
    work = output / 'run'
    work.mkdir(parents=True, exist_ok=True)
    ledger = dict(status='running', started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        slurm_job_id=os.environ.get('SLURM_JOB_ID'), config=cfg, steps=[],
        inputs={key:stamp(cfg[key]) for key in ('gam','index','gfa','gbz')})
    started = time.perf_counter()
    def save():
        ledger['elapsed_seconds'] = time.perf_counter()-started
        write_json(work/'status.json', ledger)
    def stage(name, command=None, function=None):
        record = dict(name=name, status='running', started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        if command: record['command'] = list(map(str,command))
        ledger['steps'].append(record);save()
        t = time.perf_counter()
        cpu = resource.getrusage(resource.RUSAGE_SELF)
        try:
            if command:
                log = work/(name+'.log')
                record['log'] = str(log)
                with log.open('w') as stream:
                    subprocess.run(['/usr/bin/time','-v','-o',str(work/(name+'.resources.txt')),
                        *map(str,command)],stdout=stream,stderr=subprocess.STDOUT,check=True,cwd=ROOT)
            else:
                record['result'] = function()
            record['status'] = 'complete'
        except BaseException as e:
            record.update(status='failed',error=str(e));raise
        finally:
            record['wall_seconds'] = time.perf_counter()-t
            record['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
            current_cpu = resource.getrusage(resource.RUSAGE_SELF)
            record['controller_user_cpu_seconds'] = current_cpu.ru_utime-cpu.ru_utime
            record['controller_system_cpu_seconds'] = current_cpu.ru_stime-cpu.ru_stime
            save()
    py = sys.executable
    cli = ROOT/'indexed_gam_pipeline/run.py'
    db = work/'all_graph_nodes.sqlite'
    cache = work/'gbwt_counts.sqlite'
    targets = work/'discovery/target_nodes.txt'
    def build_command(dest, smoke=False):
        command = [py,cli,'build','--format','candidate-v4','--gam',cfg['gam'],'--index',cfg['index'],
            '--nodes',cfg['preflight_nodes'] if smoke else targets,'--node-sqlite',db,'--gbz',cfg['gbz'],'--gbz-query',cfg['gbz_query'],
            '--occurrence-cache',cache,'--output',dest,'--batch-nodes',str(cfg['batch_nodes']),
            '--max-node-span','10000','--gam-cache-mb',str(cfg['gam_cache_mb']),
            '--max-batch-segments','20000','--shard-size',str(8 if smoke else cfg['shard_size']),
            '--rows','200','--width','101','--min-mapq','10','--min-af','0.05',
            '--min-variants','3','--min-allele-bq','10','--max-indel-len','50','--variant-type','all']
        if smoke:command.extend(['--max-tensors','8','--debug-rows'])
        return command
    try:
        save()
        stage('vg_version',[cfg['vg'],'version'])
        stage('graph_sequence_index',function=lambda:graph_index(cfg['gfa'],db))
        stage('preflight_build',build_command(work/'preflight',True))
        preflight = json.loads((work/'preflight/manifest.json').read_text())
        if preflight.get('status') != 'complete' or preflight.get('tensors', 0) == 0:
            raise ValueError('Preflight must produce at least one completed tensor')
        stage('preflight_source_audit',[py,ROOT/'indexed_gam_pipeline/validate_examples.py',work/'preflight',
            '--gam',cfg['gam'],'--index',cfg['index'],'--node-sqlite',db,'--gbz',cfg['gbz'],
            '--gbz-query',cfg['gbz_query'],'--occurrence-cache',cache,'--output',work/'preflight_audit.json'])
        stage('discovery',[py,cli,'discover','--gam',cfg['gam'],'--output',work/'discovery',
            '--min-mapq','5','--node-alt','0.05'])
        if stamp(cfg['gam']) != ledger['inputs']['gam'] or stamp(cfg['index']) != ledger['inputs']['index']:
            raise ValueError('Input GAM/index changed during discovery')
        building = output/'.building_tensors'
        stage('tensor_build',build_command(building))
        stage('shard_validation',function=lambda:validate_shards(building,cfg['shard_size']))
        if any(stamp(cfg[key]) != ledger['inputs'][key] for key in ledger['inputs']):
            raise ValueError('An input changed during the run')
        def publish():
            files = list(building.iterdir())
            for path in files:
                if (output/path.name).exists():raise ValueError(f'Output already exists: {path.name}')
            for path in files:path.rename(output/path.name)
            building.rmdir()
            return dict(files=len(files),directory=str(output))
        stage('publish',function=publish)
        ledger['status']='complete'
    except BaseException as e:
        ledger.update(status='failed',error=str(e));raise
    finally:
        save()


if __name__ == '__main__':
    main()
