"""Rebuild tensors from a verified prior run using bounded, disposable task processes."""
import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import sys
import time

from indexed_gam_pipeline.full_run import stamp, write_json, validate_shards
from indexed_gam_pipeline.tensor_storage import STORAGE_VERSION
from indexed_gam_pipeline.full_reader_comparison import build_command
from indexed_gam_pipeline.unified_full_run import digest, partition_nodes, verify, MemoryRecorder


def prepare(prior, root, tasks=512, processes=32, cache_mb=8192, snv_min_af=None, indel_min_af=None,
            max_node_reads=800, candidate_unit='site'):
    prior, root = Path(prior).resolve(), Path(root).resolve()
    if tasks < 1 or not 1 <= processes <= tasks or cache_mb < 1:
        raise ValueError('Require tasks >= processes >= 1 and positive cache size')
    if max_node_reads < 0 or candidate_unit not in ('site', 'allele'):
        raise ValueError('Require read cap >= 0 and candidate unit site|allele')
    split = snv_min_af is not None or indel_min_af is not None
    if split and not all(v is not None and 0 <= v <= 1 for v in (snv_min_af, indel_min_af)):
        raise ValueError('Specify both SNV and INDEL AF thresholds in [0,1]')
    old = json.loads((prior/'config.json').read_text())
    verify(old)
    status = json.loads((prior/'status.json').read_text())
    stages = {s['name']: s['status'] for s in status['stages']}
    for name in ('verify_inputs_and_partitions', 'build_unified_graph_index',
                 'independent_graph_count_audit', 'preflight_20_processes'):
        if stages.get(name) != 'complete':
            raise ValueError('Cannot reuse incomplete stage: '+name)
    if not json.loads((prior/'graph_audit.json').read_text())['passed']:
        raise ValueError('Prior graph audit did not pass')
    root.mkdir(parents=True, exist_ok=False)
    # Freeze the current pipeline, including its storage encoding changes.
    source_root = Path(__file__).resolve().parents[1]
    (root/'source').mkdir()
    for folder in ('indexed_gam_pipeline', 'src'):
        shutil.copytree(source_root/folder, root/'source'/folder,
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    parts = partition_nodes(old['config']['discovery_nodes'], root, tasks, old['total_nodes'])
    config = dict(old, parts=parts, reused_from=str(prior),
                  tensor_storage_version=STORAGE_VERSION, source_root=str(source_root),
                  preflight_reused_from=str(prior),
                  params=dict(old['params'], processes=processes, gam_cache_mb_per_process=cache_mb,
                              max_node_reads=max_node_reads, candidate_unit=candidate_unit, early_af_filter=True),
                  graph_index=stamp(old['config']['unified_graph_index']),
                  source_sha256={str(p):digest(p) for p in (root/'source').rglob('*') if p.is_file()})
    config['params']['decode_policy'] = 'full'
    if split:
        config['generation_mode'] = 'single-decode-split-output-v1'
        config['variant_outputs'] = [dict(folder='SNV', variant_type='snp', min_af=snv_min_af),
                                     dict(folder='INDEL', variant_type='indel', min_af=indel_min_af)]
    config['prior_source_commit'] = config.pop('source_commit', None)
    # The copied file hashes, not the prior run's commit, identify this snapshot.
    write_json(root/'config.json', config)
    script = '\n'.join(['#!/bin/bash', 'set -euo pipefail',
        'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1',
        'export PYTHONUNBUFFERED=1', 'cd '+shlex.quote(str(root/'source')),
        'exec /usr/bin/time -v -o '+shlex.quote(str(root/'job.resources.txt'))+' '+shlex.quote(sys.executable)+
        ' -m indexed_gam_pipeline.dynamic_tensor_run run --root '+shlex.quote(str(root)), ''])
    (root/'run.sh').write_text(script)
    return config


def execute_queue(commands, folder, processes, poll_seconds=.5):
    """Launch another independent process whenever a slot becomes free; fail fast."""
    folder = Path(folder)
    folder.mkdir(exist_ok=False)
    ledger = dict(status='running', job_id=os.environ.get('SLURM_JOB_ID'),
                  processes=processes, total_tasks=len(commands), completed_tasks=0,
                  started_unix=time.time(), tasks=[dict(task=i, status='pending') for i in range(len(commands))])
    active = {}; next_task = 0
    def save():
        ledger['elapsed_seconds'] = time.time()-ledger['started_unix']
        write_json(folder/'status.json', ledger)
    def interrupted(signum, frame):
        raise RuntimeError('Controller signal '+str(signum))
    previous = {s:signal.signal(s, interrupted) for s in (signal.SIGTERM,signal.SIGINT)}
    save()
    try:
        while next_task < len(commands) or active:
            for i, (proc, log) in list(active.items()):
                code = proc.poll()
                if code is None:
                    continue
                log.close(); del active[i]
                item = ledger['tasks'][i]
                item.update(status='complete' if code == 0 else 'failed', returncode=code,
                            wall_seconds=time.time()-item['started_unix'])
                if code:
                    raise RuntimeError(f'Task {i} exited with code {code}')
                ledger['completed_tasks'] += 1
                print(f"COMPLETE task {i}: {ledger['completed_tasks']}/{len(commands)}", flush=True)
            while next_task < len(commands) and len(active) < processes:
                i = next_task
                log = (folder/f'task_{i:04d}.log').open('w')
                try:
                    proc = subprocess.Popen(commands[i], stdout=log, stderr=subprocess.STDOUT,
                                            start_new_session=True)
                except BaseException:
                    log.close(); raise
                active[i] = (proc, log)
                ledger['tasks'][i].update(status='running', pid=proc.pid, started_unix=time.time())
                next_task += 1
            save()
            if active:
                time.sleep(poll_seconds)
        ledger['status'] = 'complete'
    except BaseException as error:
        ledger.update(status='failed', error=str(error)); raise
    finally:
        for proc, log in active.values():
            if proc.poll() is None:
                try: os.killpg(proc.pid, signal.SIGTERM)
                except ProcessLookupError: pass
        deadline = time.monotonic()+10
        for i, (proc, log) in active.items():
            try: proc.wait(timeout=max(.01, deadline-time.monotonic()))
            except subprocess.TimeoutExpired:
                try: os.killpg(proc.pid, signal.SIGKILL)
                except ProcessLookupError: pass
                proc.wait()
            log.close()
            ledger['tasks'][i].update(status='cancelled', returncode=proc.returncode)
        save()
        for s, handler in previous.items(): signal.signal(s, handler)
    return ledger


def task(root, index):
    config = json.loads((root/'config.json').read_text())
    outputs = config.get('variant_outputs', [dict(folder='python', variant_type='all', min_af=.05)])
    dest = root/'python'/f'task_{index:04d}'
    dest.parent.mkdir(parents=True, exist_ok=True)
    command = build_command(config, config['parts'][index]['nodes_file'], dest, None,
                            'python', config['params']['gam_cache_mb_per_process'])
    if config['params'].get('decode_policy'):
        command += ['--decode-mode', config['params']['decode_policy']]
    if 'variant_outputs' in config:
        for output in outputs:
            flag = 'snv' if output['variant_type'] == 'snp' else 'indel'
            path = root/output['folder']/f'task_{index:04d}'
            command += [f'--{flag}-output', str(path), f'--{flag}-min-af', str(output['min_af'])]
    # One builder decodes both types, exits, then each output is validated.
    subprocess.run(['/usr/bin/time', '-v', '-o', str(dest)+'.resources.txt', *command], check=True)
    for output in outputs:
        path = root/output['folder']/f'task_{index:04d}'
        report = validate_shards(path, 2048, width=101)
        if config.get('tensor_storage_version') == STORAGE_VERSION and report['dtype'] != 'int8':
            raise ValueError('Task did not produce configured int8 storage')
        allowed = {'SNP'} if output['variant_type'] == 'snp' else {'INS', 'DEL'} if output['variant_type'] == 'indel' else {'SNP', 'INS', 'DEL'}
        if set(report['event_types']) - allowed:
            raise ValueError('Unexpected event type in separated output')


def run(root):
    config = json.loads((root/'config.json').read_text())
    processes = config['params']['processes']
    if int(os.environ.get('SLURM_CPUS_PER_TASK', processes)) < processes:
        raise ValueError('Not enough allocated CPUs')
    verify(config)
    if stamp(config['config']['unified_graph_index']) != config['graph_index']:
        raise ValueError('Reused graph index changed')
    status = dict(status='running', job_id=os.environ.get('SLURM_JOB_ID'),
                  started_unix=time.time(), stage='tensor_build_and_validation',
                  tasks=len(config['parts']), processes=processes, reused_from=config['reused_from'])
    write_json(root/'status.json', status)
    commands = [[sys.executable, '-m', 'indexed_gam_pipeline.dynamic_tensor_run',
                 'task', '--root', str(root), '--index', str(i)] for i in range(len(config['parts']))]
    try:
        with MemoryRecorder(root, lambda:status['stage']):
            execute_queue(commands, root/'python', processes)
            verify(config)
            if stamp(config['config']['unified_graph_index']) != config['graph_index']:
                raise ValueError('Reused graph index changed during build')
            totals = {output['folder']: sum(json.loads(
                (root/output['folder']/f'task_{i:04d}'/'manifest.json').read_text())['tensors']
                for i in range(len(commands)))
                for output in config.get('variant_outputs', [dict(folder='python')])}
            status.update(status='complete', tensors=sum(totals.values()), tensors_by_type=totals)
    except BaseException as error:
        status.update(status='failed', error=str(error)); raise
    finally:
        status['elapsed_seconds'] = time.time()-status['started_unix']
        write_json(root/'status.json', status)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'run', 'task'])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--prior', type=Path)
    parser.add_argument('--tasks', type=int, default=512)
    parser.add_argument('--processes', type=int, default=32)
    parser.add_argument('--cache-mb', type=int, default=8192)
    parser.add_argument('--index', type=int)
    parser.add_argument('--snv-min-af', type=float)
    parser.add_argument('--indel-min-af', type=float)
    parser.add_argument('--max-node-reads', type=int, default=800)
    parser.add_argument('--candidate-unit', choices=['site', 'allele'], default='site')
    args = parser.parse_args()
    if args.action == 'prepare':
        if args.prior is None: parser.error('--prior is required')
        prepare(args.prior, args.root, args.tasks, args.processes, args.cache_mb, args.snv_min_af, args.indel_min_af,
                args.max_node_reads, args.candidate_unit)
    elif args.action == 'task':
        if args.index is None: parser.error('--index is required')
        task(args.root, args.index)
    else:
        run(args.root)


if __name__ == '__main__':
    main()
