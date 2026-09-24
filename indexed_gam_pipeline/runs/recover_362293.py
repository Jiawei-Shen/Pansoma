"""Prepare isolated recovery partitions and publish an exclusive result catalog."""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import sys


def check(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(prior, root, repo):
    sys.path.insert(0, str(prior / 'source'))
    from indexed_gam_pipeline.unified_full_run import verify, partition_nodes
    from indexed_gam_pipeline.full_run import write_json
    config = read(prior / 'config.json')
    ledger = read(prior / 'python/status.json')
    check(ledger['status'] == 'failed', 'Prior run is not stopped/failed')
    verify(config)
    tasks = ledger['tasks']
    check([t['task'] for t in tasks] == list(range(512)), 'Unexpected task inventory')
    completed = [t['task'] for t in tasks if t['status'] == 'complete']
    remaining = [t['task'] for t in tasks if t['status'] != 'complete']
    check(len(completed) == 495 and len(remaining) == 17, 'Prior status changed')
    check(remaining == [38,175,258,263,270,271,276,285,393,440,441,445,486,489,490,498,511],
          'Unexpected recovery targets')
    # Original partitions must form one strictly increasing sequence. Recovery
    # membership comes exclusively from the unfinished entries in that sequence.
    previous = count = 0
    recovery_nodes = []
    for i, part in enumerate(config['parts']):
        n = 0
        with Path(part['nodes_file']).open() as stream:
            for line in stream:
                node = int(line)
                check(node > previous, 'Duplicate/unsorted original node')
                previous = node
                n += 1
                if i in remaining:
                    recovery_nodes.append(node)
        check(n == part['nodes'], 'Original partition count mismatch')
        count += n
    check(count == config['total_nodes'], 'Original node total mismatch')
    check(len(recovery_nodes) == 549850, 'Unexpected recovery node count')
    preserved = []
    for i in completed:
        for kind in ('SNV', 'INDEL'):
            folder = prior / kind / f'task_{i:04d}'
            manifest = read(folder / 'manifest.json')
            check(manifest['status'] == 'complete', 'Incomplete preserved output')
            preserved.append(dict(original_task=i, kind=kind, path=str(folder),
                                  manifest_sha256=sha(folder/'manifest.json'),
                                  tensors=manifest['tensors']))
    archive_targets = []
    for i in remaining:
        for kind in ('python', 'SNV', 'INDEL'):
            source = prior / kind / f'task_{i:04d}'
            check(source.is_dir() and not source.is_symlink(), f'Invalid archive target: {source}')
            archive_targets.append((kind, i, source))
    root.mkdir(parents=True, exist_ok=False)
    (root/'source').mkdir()
    for folder in ('indexed_gam_pipeline', 'src'):
        shutil.copytree(prior/'source'/folder, root/'source'/folder,
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    shutil.copy2(repo/'indexed_gam_pipeline/gam_reader.py',
                 root/'source/indexed_gam_pipeline/gam_reader.py')
    shutil.copy2(Path(__file__), root/'recovery.py')
    union = root/'recovery_nodes.txt'
    with union.open('x') as stream:
        stream.writelines(f'{node}\n' for node in recovery_nodes)
    parts = partition_nodes(union, root, 512, len(recovery_nodes))
    rebuilt = []
    for part in parts:
        rebuilt.extend(map(int, Path(part['nodes_file']).read_text().split()))
    check(rebuilt == recovery_nodes, 'Recovery coverage differs or duplicates nodes')
    config.update(parts=parts, total_nodes=len(recovery_nodes), reused_from=str(prior),
                  source_root=str(repo), recovery_of=str(prior), original_tasks=remaining,
                  source_sha256={str(p):sha(p) for p in (root/'source').rglob('*') if p.is_file()})
    config['params']['processes'] = 48
    write_json(root/'config.json', config)
    verify(config)
    inventory = dict(status='prepared', original_run=str(prior), recovery_run=str(root),
                     original_total_nodes=count, recovery_nodes=len(recovery_nodes),
                     recovery_tasks=512, original_tasks=remaining,
                     preserved_outputs=preserved, archived_outputs=[])
    write_json(root/'recovery_inventory.json', inventory)
    # Move only the 17 obsolete partial task directories. No completed directory
    # or underlying shard is modified. Archives are excluded from the catalog.
    for kind, i, source in archive_targets:
        destination = root/'excluded_partial_backup'/kind/f'task_{i:04d}'
        destination.parent.mkdir(parents=True, exist_ok=True)
        source.rename(destination)
        inventory['archived_outputs'].append(dict(source=str(source), backup=str(destination)))
        write_json(root/'recovery_inventory.json', inventory)
    script = '\n'.join([
        '#!/bin/bash', 'set -euo pipefail',
        'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1',
        'export PYTHONUNBUFFERED=1', 'cd '+shlex.quote(str(root/'source')),
        '/usr/bin/time -v -o '+shlex.quote(str(root/'job.resources.txt'))+' '+shlex.quote(sys.executable)+
        ' -m indexed_gam_pipeline.dynamic_tensor_run run --root '+shlex.quote(str(root)),
        shlex.quote(sys.executable)+' '+shlex.quote(str(root/'recovery.py'))+
        ' finalize --root '+shlex.quote(str(root)), ''])
    (root/'run.sh').write_text(script)
    print(json.dumps(dict(root=str(root), original_tasks=remaining, tasks=512, processes=48,
                          recovery_nodes=len(recovery_nodes), preserved_tasks=495,
                          archived_directories=len(archive_targets)), indent=2), flush=True)


def finalize(root):
    sys.path.insert(0, str(root/'source'))
    from indexed_gam_pipeline.unified_full_run import verify
    from indexed_gam_pipeline.full_run import write_json
    inventory = read(root/'recovery_inventory.json')
    config = read(root/'config.json')
    verify(config)
    ledger = read(root/'python/status.json')
    check(read(root/'status.json')['status'] == 'complete', 'Recovery did not complete')
    check(ledger['status'] == 'complete' and ledger['completed_tasks'] == 512,
          'Recovery tasks not all validated')
    entries = inventory['preserved_outputs'][:]
    for entry in entries:
        check(sha(Path(entry['path'])/'manifest.json') == entry['manifest_sha256'],
              'Preserved output manifest changed')
    for i in range(512):
        for kind in ('SNV', 'INDEL'):
            folder = root/kind/f'task_{i:04d}'
            manifest = read(folder/'manifest.json')
            check(manifest['status'] == 'complete', 'Incomplete recovery output')
            entries.append(dict(recovery_task=i, kind=kind, path=str(folder),
                                tensors=manifest['tensors'], manifest_sha256=sha(folder/'manifest.json')))
    check(len({e['path'] for e in entries}) == len(entries), 'Duplicate output directory')
    totals = {kind:sum(e['tensors'] for e in entries if e['kind'] == kind) for kind in ('SNV','INDEL')}
    write_json(root/'combined_outputs.json', dict(status='complete',
               total_nodes=inventory['original_total_nodes'], tensors_by_type=totals,
               tensors=sum(totals.values()), outputs=entries,
               excluded_original_tasks=inventory['original_tasks']))
    inventory['status'] = 'complete'
    write_json(root/'recovery_inventory.json', inventory)
    print(json.dumps(dict(combined_tensors=totals, catalog=str(root/'combined_outputs.json'))), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'finalize'))
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--prior', type=Path)
    parser.add_argument('--repo', type=Path)
    args = parser.parse_args()
    if args.action == 'prepare':
        prepare(args.prior.resolve(), args.root.resolve(), args.repo.resolve())
    else:
        finalize(args.root.resolve())
