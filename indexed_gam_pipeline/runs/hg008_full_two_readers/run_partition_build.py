import json,sys
from pathlib import Path
root=Path(__file__).resolve().parent
sys.path.insert(0,str(root/'source'))
from indexed_gam_pipeline.full_reader_comparison import verify_inputs,run_builders
mode=sys.argv[1]
if mode not in ('python','vg'):raise ValueError('Invalid backend')
config=json.loads((root/'comparison_config.json').read_text())
if json.loads((root/'preflight_status.json').read_text())['status']!='complete':raise ValueError('Preflight incomplete')
verify_inputs(config)
run_builders(config,root/mode,[Path(p['nodes_file']) for p in config['parts']],mode)
