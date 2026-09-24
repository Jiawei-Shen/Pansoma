"""Repeat the 17-task recovery with the same frozen code and 32 builders."""
import json
from pathlib import Path
import shutil
import sys

base = Path('/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/run')
prior = base/'recovery17_p48_512_nfilter_eoffix_20260922'
root = base/'recovery17_p32_512_nfilter_eoffix_20260922'
sys.path.insert(0, str(prior/'source'))
from indexed_gam_pipeline.unified_full_run import verify, digest

config = json.loads((prior/'config.json').read_text())
inventory = json.loads((prior/'recovery_inventory.json').read_text())
verify(config)
assert config['total_nodes'] == 549850 and len(config['parts']) == 512
assert len(inventory['preserved_outputs']) == 990
for entry in inventory['preserved_outputs']:
    assert digest(Path(entry['path'])/'manifest.json') == entry['manifest_sha256']
root.mkdir(exist_ok=False)
shutil.copytree(prior/'source', root/'source',
                ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
for part in config['parts']:
    for key in ('nodes_file', 'preflight_file'):
        source = Path(part[key])
        target = root/source.name
        shutil.copy2(source, target)
        part[key] = str(target)
config['params']['processes'] = 32
config['replaced_recovery'] = str(prior)
config['source_sha256'] = {str(p):digest(p) for p in (root/'source').rglob('*') if p.is_file()}
assert config['params']['gam_cache_mb_per_process'] == 8192
(root/'config.json').write_text(json.dumps(config, indent=2))
inventory.update(status='prepared', recovery_run=str(root), replaced_recovery=str(prior))
(root/'recovery_inventory.json').write_text(json.dumps(inventory, indent=2))
shutil.copy2(prior/'recovery.py', root/'recovery.py')
(root/'run.sh').write_text((prior/'run.sh').read_text().replace(str(prior), str(root)))
verify(config)
print(json.dumps(dict(root=str(root), nodes=549850, tasks=512, processes=32)))
