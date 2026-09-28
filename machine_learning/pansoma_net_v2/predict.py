"""Class probabilities of a trained PansomaNetV2 for every tensor of merged tensor sets.

    cd machine_learning
    python -m pansoma_net_v2.predict --checkpoint <run>/best.pth --tensors <set> [...] --output <dir> [--chroms chr20]

Writes <output>/<sample>.<set>.<KIND>.predictions.tsv.gz (chrom, candidate_id, label, p_non, p_somatic,
p_germline, pred; every tensor, labels -1 included, in index order; pred uses the checkpoint's somatic
threshold, --threshold overrides it) and .metrics.json (metrics.report over the tensors labelled 0/1/2:
argmax and thresholded precision / recall / F1, somatic and germline average precision). The encoder uses the
checkpoint's statistics; they are not refitted on these tensors.
"""
import argparse
import gzip
import json
from pathlib import Path

import torch

from . import metrics
from .data import CLASSES, KINDS, EpochSampler, TensorDataset, load_parts
from .model import PansomaNetV2
from .train import describe, make_loader, predict_probs


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--checkpoint", required=True)
    p.add_argument("--tensors", nargs="+", required=True)
    p.add_argument("--kinds", nargs="+", default=list(KINDS), choices=KINDS)
    p.add_argument("--chroms", nargs="+", help="default: all chromosomes of the sets")
    p.add_argument("--output", required=True)
    p.add_argument("--cache-dir", help="index cache (default: <output>/index_cache)")
    p.add_argument("--batch-size", type=int, default=256)
    p.add_argument("--num-workers", type=int, default=8)
    p.add_argument("--amp", choices=["bf16", "off"], default="bf16")
    p.add_argument("--threshold", type=float, help="somatic threshold (default: the checkpoint's)")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    checkpoint = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    model = PansomaNetV2.from_checkpoint(checkpoint).to(device).eval()
    amp = args.amp == "bf16" and device.type == "cuda"
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    print(f"checkpoint {args.checkpoint} (epoch {checkpoint['epoch']}), statistics {json.dumps(model.encoder.stats())}",
          flush=True)
    for index, positions in load_parts(args.tensors, args.kinds, args.cache_dir or str(out / "index_cache"),
                                       chroms=args.chroms, labelled=False):
        dataset = TensorDataset([(index, positions)])
        loader = make_loader(dataset, EpochSampler(len(dataset), shuffle=False), args.num_workers, args.batch_size,
                             persistent=False)
        labels, probs = predict_probs(model, loader, device, amp)
        threshold = args.threshold if args.threshold is not None else checkpoint.get("somatic_threshold")
        pred = metrics.call(probs, threshold)
        name = f"{index.dir.parent.parent.name}.{index.dir.parent.name}.{index.kind}"
        chroms, candidates = index.meta["chroms"], index.candidates()
        with gzip.open(out / f"{name}.predictions.tsv.gz", "wt") as f:
            f.write("chrom\tcandidate_id\tlabel\t" + "\t".join(f"p_{c}" for c in CLASSES) + "\tpred\n")
            for k, pos in enumerate(positions):
                f.write(f"{chroms[index.arrays['chrom'][pos]]}\t{candidates[pos]}\t{labels[k]}\t"
                        + "\t".join(f"{p:.5f}" for p in probs[k]) + f"\t{CLASSES[pred[k]]}\n")
        report = metrics.report(labels, probs, threshold)
        report.update(all_tensors=int(len(labels)), ignored=int((labels < 0).sum()), directory=str(index.dir),
                      labels=index.meta["labels"], checkpoint=str(Path(args.checkpoint).resolve()), chroms=args.chroms)
        (out / f"{name}.metrics.json").write_text(json.dumps(report, indent=2) + "\n")
        print(f"{name}: {len(labels):,} tensors ({report['tensors']:,} labelled) | {describe(report)}", flush=True)


if __name__ == "__main__":
    main()
