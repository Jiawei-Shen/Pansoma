# usage: python make_prepare.py SET  -> <work>/<SET>/prepare.sh
# The prepare command of <SET>_snvaf05: the run of pansoma_v2_tensors/<SET> (its run/config.json: GAM, discovery node
# list and node stats, graph index, haplotypes, task count, processes, every builder parameter, chromosome index,
# merge, reference path) and the truth of its current SNV labels (tensors/SNV/labels.manifest.json), with one change:
# --snv-min-af 0.05 (was 0.08). The discovery is the old one (node selection at imperfect fraction > 0.05 does not
# depend on the AF floor). check_config.py compares the two configs after prepare.
import json, shlex, sys
from pathlib import Path
S = sys.argv[1]
OLD = Path(f"/scratch/jshen/data/pansoma_v2_tensors/{S}"); NEW = Path(f"/scratch/jshen/data/pansoma_v2_tensors/{S}_snvaf05")
W = Path("/scratch/jshen/Github/Pansoma/tmp/snv_af05_20261103")
c = json.loads((OLD / "run/config.json").read_text()); m = json.loads((OLD / "tensors/SNV/labels.manifest.json").read_text())
b = c["builder"]; assert c["variant_outputs"] == {"SNV": 0.08, "INDEL": 0.08}, c["variant_outputs"]
gams = c["inputs"]["gam"] if isinstance(c["inputs"]["gam"], list) else [c["inputs"]["gam"]]
assert len(gams) == 1, gams
args = ["--root", NEW / "run", "--tensors", NEW / "tensors", "--gam", gams[0]["path"], "--index", c["inputs"]["index"]["path"],
        "--nodes", OLD / "discovery/target_nodes.txt", "--node-stats", OLD / "discovery/node_stats.json",
        "--graph-index", c["inputs"]["graph_index"]["path"], "--haplotypes", b["haplotypes"],
        "--tasks", c["tasks"], "--processes", c["processes"],
        "--snv-min-af", 0.05, "--indel-min-af", 0.08,
        "--rows", b["rows"], "--width", b["width"], "--gam-cache-mb", b["gam_cache_mb"], "--batch-nodes", b["batch_nodes"],
        "--max-node-span", b["max_node_span"], "--max-batch-alignments", b["max_batch_alignments"],
        "--shard-size", b["shard_size"], "--min-mapq", b["min_mapq"], "--min-af", b["min_af"],
        "--min-variants", b["min_variants"], "--min-allele-bq", b["min_allele_bq"], "--max-indel-len", b["max_indel_len"],
        "--max-node-reads", b["max_node_reads"], "--chromosomes", b["chromosomes"],
        "--chr-index", c["chromosome_selection"]["chr_index"]["path"],
        "--early-af-filter" if b["early_af_filter"] else "--no-early-af-filter", "--decoder", b["decoder"],
        "--merge-shard-size", c["postprocess"]["merge_shard_size"], "--reference-path", m["reference_path"],
        "--somatic-vcf", m["truth"]["somatic"]["vcf"], "--somatic-bed", m["truth"]["somatic"]["bed"],
        "--germline-vcf", m["truth"]["germline"]["vcf"], "--germline-bed", m["truth"]["germline"]["bed"],
        "--reference-fasta", m["fasta"], "--truth-dir", NEW / "truth"]
if m.get("ignore_repeat_indels"): args.append("--label-ignore-repeat-indels")
for k in ("snv_min_af", "indel_min_af"):
    assert m.get(k) is None, (k, m.get(k))  # no label floor in the current labels; none given here
(W / S).mkdir(parents=True, exist_ok=True)
(W / S / "prepare.sh").write_text(f"""#!/bin/bash
# prepare of {S}_snvaf05 (make_prepare.py): {OLD}/run with --snv-min-af 0.05; code {W}/code
set -euo pipefail
[ ! -e {NEW} ] || {{ echo "{NEW} exists"; exit 1; }}
mkdir -p {NEW}
cd {W}/code
/usr/bin/time -v /wanglab/jshen/anaconda3/bin/python -m indexed_gam_pipeline_v4.orchestrate prepare \\
  {' '.join(shlex.quote(str(a)) for a in args)}
/wanglab/jshen/anaconda3/bin/python {W}/check_config.py {S}
""")
print(W / S / "prepare.sh")
