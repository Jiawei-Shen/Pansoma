"""chr1 HG008T Illumina INDEL tensors that got label 1 by partial matching: figures + index.tsv + README tables.

Selects residual_partial_somatic_truth (GRCh38 nodes and branch nodes) and allele_partial_somatic_truth from
<tensors>/INDEL/chr1_labels.ndjson, renders each with scripts/visualize_tensor.py (title: label, partial truth,
context, model scores), saves the PNG with a 256-colour palette (about a third of the size, same look) under
png/<group>/, and writes index.tsv and tables.md (the per-group tables of README.md).

Run from the repository root with a Python that has matplotlib (the base env's matplotlib is broken under NumPy 2):
  /wanglab/jshen/anaconda3/envs/polymarket-btc-5m-bot/bin/python analysis/indel_partial_tensors_20260930/render.py
"""
import gzip
import json
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

from PIL import Image

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
TENSORS = Path("/scratch/jshen/data/pansoma_v2_tensors/HG008T_Illumina/tensors/INDEL")
FASTA = Path("/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta")
RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs")
MODELS = {"base": "HG008_Illumina_INDEL_base_b1024",           # partial somatic trained as 1
          "nopartial": "HG008_Illumina_INDEL_nopartial_b1024"}  # partial somatic ignored (--ignore-reasons)
GROUPS = [("residual_on_ref", "residual_partial_somatic_truth, GRCh38 node"),
          ("residual_off_ref", "residual_partial_somatic_truth, branch node (off-reference)"),
          ("allele_partial", "allele_partial_somatic_truth")]


class Fasta:
    """Minimal faidx reader (the plotting env has no pysam)."""

    def __init__(self, path):
        self.handle = open(path, "rb")
        self.index = {}
        for line in open(f"{path}.fai"):
            name, length, offset, bases, width = line.split("\t")[:5]
            self.index[name] = (int(length), int(offset), int(bases), int(width))

    def fetch(self, chrom, start, end):
        length, offset, bases, width = self.index[chrom]
        start, end = max(0, start), min(end, length)
        out, pos = [], start
        while pos < end:
            line, col = divmod(pos, bases)
            n = min(bases - col, end - pos)
            self.handle.seek(offset + line * width + col)
            out.append(self.handle.read(n).decode())
            pos += n
        return "".join(out).upper()


def context(fasta, truth):
    """Repeat context of the truth INDEL on GRCh38: HP<n>, STR <unit>x<copies> or other."""
    ref, alt, pos = truth["vcf_ref"], truth["vcf_alt"], truth["vcf_pos"]   # VCF: first base is the shared anchor
    event = ref[1:] if len(ref) > len(alt) else alt[1:]
    if not event:
        return "other"
    unit = next(event[:k] for k in range(1, len(event) + 1) if len(event) % k == 0 and event[:k] * (len(event) // k) == event)
    seq = fasta.fetch(truth["chrom"], pos - 200, pos + 200)   # seq[200] = the base after the anchor (0-based pos)
    left, right = 200, 200
    while right + len(unit) <= len(seq) and seq[right:right + len(unit)] == unit:
        right += len(unit)
    while left - len(unit) >= 0 and seq[left - len(unit):left] == unit:
        left -= len(unit)
    copies = (right - left) // len(unit)
    if len(unit) == 1:
        return f"HP{copies}"
    return f"STR {unit}x{copies}" if copies >= 2 else "other"


def load_scores():
    scores, thresholds = {}, {}
    for key, run in MODELS.items():
        test = RUNS / run / "test_chr1"
        thresholds[key] = json.load(open(next(test.glob("*.metrics.json"))))["threshold"]
        for line in gzip.open(next(test.glob("*.predictions.ndjson.gz")), "rt"):
            x = json.loads(line)
            scores.setdefault(x["candidate_id"], {})[key] = x["p_somatic"]
    return scores, thresholds


def main():
    fasta = Fasta(FASTA)
    scores, thresholds = load_scores()
    labels = [json.loads(line) for line in open(TENSORS / "chr1_labels.ndjson")]
    summaries = [json.loads(line) for line in open(TENSORS / "chr1_variant_summary.ndjson")]
    rows = []
    for lab, summ in zip(labels, summaries):
        assert lab["candidate_id"] == summ["candidate_id"]
        if lab["reason"] == "residual_partial_somatic_truth":
            group = "residual_on_ref" if lab["grch38"] else "residual_off_ref"
        elif lab["reason"] == "allele_partial_somatic_truth":
            group = "allele_partial"
        else:
            continue
        truth = lab["partial_truth"]
        locus = (f"{lab['grch38']['chrom']}:{lab['grch38']['pos0'] + 1}" if lab["grch38"]
                 else f"anchor {lab['anchor']['chrom']}:{lab['anchor']['start']}-{lab['anchor']['end']}")
        sort_pos = lab["grch38"]["pos0"] if lab["grch38"] else lab["anchor"]["start"]
        s = scores.get(lab["candidate_id"], {})
        rows.append(dict(group=group, sort=sort_pos, candidate_id=lab["candidate_id"], reason=lab["reason"],
                         overlap=lab["overlap"], locus=locus,
                         truth=f"{truth['chrom']}:{truth['vcf_pos']} {truth['vcf_ref']}>{truth['vcf_alt']} ({truth['gt']})",
                         context=context(fasta, truth), af=summ["af"], alt_count=summ["alt_count"],
                         coverage=summ["coverage"], shard=lab["shard_file"], index=lab["index_within_shard"],
                         p_base=s.get("base"), p_nopartial=s.get("nopartial")))
    rows.sort(key=lambda r: ([g for g, _ in GROUPS].index(r["group"]), r["sort"]))
    number = {}
    with tempfile.TemporaryDirectory() as tmp:
        for r in rows:
            k = number[r["group"]] = number.get(r["group"], 0) + 1
            name = r["candidate_id"].translate(str.maketrans(":>@+", "____"))
            r["file"] = f"png/{r['group']}/{k:03d}_{name}.png"
            called = {m: r[f"p_{m}"] is not None and r[f"p_{m}"] >= thresholds[m] for m in MODELS}
            title = (f"{r['candidate_id']} | label 1 {r['reason']} | overlap {r['overlap']:.3f}\n"
                     f"partial truth {r['truth']} | {r['context']} | {r['locus']}\n"
                     f"p_somatic: base_b1024 (partials trained as 1) {r['p_base']:.3f}"
                     f"{' called' if called['base'] else ''} (t {thresholds['base']:.3f}) | "
                     f"nopartial_b1024 (partials ignored) {r['p_nopartial']:.3f}"
                     f"{' called' if called['nopartial'] else ''} (t {thresholds['nopartial']:.3f})")
            r["called_base"], r["called_nopartial"] = called["base"], called["nopartial"]
            raw = Path(tmp) / "raw.png"
            subprocess.run([sys.executable, str(REPO / "scripts/visualize_tensor.py"), str(TENSORS / r["shard"]),
                            "-i", str(r["index"]), "-o", str(raw), "--title", title],
                           check=True, stdout=subprocess.DEVNULL)
            out = HERE / r["file"]
            out.parent.mkdir(parents=True, exist_ok=True)
            Image.open(raw).convert("RGB").quantize(256, method=Image.Quantize.MEDIANCUT,
                                                    dither=Image.Dither.NONE).save(out, optimize=True)
            print(r["file"], flush=True)
    cols = ["group", "file", "candidate_id", "reason", "overlap", "truth", "context", "locus", "af", "alt_count",
            "coverage", "p_base", "called_base", "p_nopartial", "called_nopartial", "shard", "index"]
    with open(HERE / "index.tsv", "w") as w:
        w.write("\t".join(cols) + "\n")
        for r in rows:
            w.write("\t".join(f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    with open(HERE / "tables.md", "w") as w:
        w.write("| 组 | n | overlap = 0.5 | 0.5–0.99 | 1.0 | HP ≥ 6 | AF 中位数 | A1 reads 中位数 | base 过 t | nopartial 过 t |\n"
                "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for g, desc in GROUPS:
            rs = [r for r in rows if r["group"] == g]
            hp6 = sum(r["context"].startswith("HP") and int(r["context"][2:]) >= 6 for r in rs)
            w.write(f"| {desc} | {len(rs)} | {sum(r['overlap'] <= 0.5 for r in rs)} | "
                    f"{sum(0.5 < r['overlap'] < 0.999 for r in rs)} | {sum(r['overlap'] >= 0.999 for r in rs)} | {hp6} | "
                    f"{statistics.median(r['af'] for r in rs):.2f} | {statistics.median(r['alt_count'] for r in rs):.0f} | "
                    f"{sum(r['called_base'] for r in rs)} | {sum(r['called_nopartial'] for r in rs)} |\n")
        for g, desc in GROUPS:
            rs = [r for r in rows if r["group"] == g]
            w.write(f"\n### {desc}（{len(rs)}）\n\n| # | 图 | partial truth | context | overlap | A1 AF（reads/coverage） | "
                    f"p base | p nopartial |\n|---:|---|---|---|---:|---|---:|---:|\n")
            for k, r in enumerate(rs, 1):
                w.write(f"| {k} | [{r['candidate_id']}]({r['file']}) | {r['truth']} | {r['context']} | {r['overlap']:.3f} | "
                        f"{r['af']:.2f}（{r['alt_count']}/{r['coverage']}） | {r['p_base']:.3f} | {r['p_nopartial']:.3f} |\n")


if __name__ == "__main__":
    main()
