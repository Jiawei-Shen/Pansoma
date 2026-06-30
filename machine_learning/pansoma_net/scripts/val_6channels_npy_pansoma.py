#!/usr/bin/env python3
import argparse
import torch
import sys
import os
import json
from tqdm import tqdm
from collections import defaultdict

# ------------------------------------------------------------
# Project imports
# ------------------------------------------------------------
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from mynet import ConvNeXtCBAMClassifier
from dataset_pansoma_npy_6ch import get_data_loader

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def print_and_log(msg, log_path):
    print(msg)
    with open(log_path, 'a', encoding='utf-8') as f:
        f.write(msg + '\n')


def _read_paths_file(file_path):
    paths = []
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith('#'):
                    continue
                paths.append(os.path.abspath(os.path.expanduser(s)))
    except Exception:
        pass
    return paths


@torch.no_grad()
def evaluate_model(model, data_loader, genotype_map, log_file):
    model.eval()

    correct = 0
    total = 0

    class_correct = defaultdict(int)
    class_total = defaultdict(int)

    # for P/R/F1
    tp = defaultdict(int)
    fp = defaultdict(int)
    fn = defaultdict(int)

    inference_results = defaultdict(list)
    idx_to_class = {v: k for k, v in genotype_map.items()}

    if not data_loader or len(data_loader) == 0:
        metrics = {
            'precision_macro': 0.0, 'recall_macro': 0.0, 'f1_macro': 0.0,
            'precision_weighted': 0.0, 'recall_weighted': 0.0, 'f1_weighted': 0.0
        }
        return 0.0, {}, {}, metrics

    for images, labels, paths in tqdm(data_loader, desc="Validating", leave=True):
        images, labels = images.to(device), labels.to(device)
        outputs = model(images)
        if isinstance(outputs, tuple):
            outputs = outputs[0]

        _, pred = torch.max(outputs, 1)
        correct += (pred == labels).sum().item()
        total += labels.size(0)

        for i in range(labels.size(0)):
            pred_idx = int(pred[i].item())
            true_idx = int(labels[i].item())
            path = paths[i]

            class_total[true_idx] += 1
            if pred_idx == true_idx:
                class_correct[true_idx] += 1
                tp[true_idx] += 1
            else:
                fp[pred_idx] += 1
                fn[true_idx] += 1

            predicted_class_name = idx_to_class.get(pred_idx, str(pred_idx))
            inference_results[predicted_class_name].append(os.path.basename(path))

    acc = (correct / total) * 100 if total > 0 else 0.0

    # per-class stats
    class_stats = {}
    for cname, cidx in genotype_map.items():
        ctot = class_total[cidx]
        ccorr = class_correct[cidx]
        cacc = (ccorr / ctot) * 100 if ctot > 0 else 0.0
        class_stats[cname] = {'acc': cacc, 'correct': ccorr, 'total': ctot, 'idx': cidx}

    # macro / weighted P/R/F1
    class_indices = list(genotype_map.values())
    precisions, recalls, f1s, supports = [], [], [], []
    for c in class_indices:
        tpc, fpc, fnc = tp[c], fp[c], fn[c]
        p = (tpc / (tpc + fpc)) if (tpc + fpc) > 0 else 0.0
        r = (tpc / (tpc + fnc)) if (tpc + fnc) > 0 else 0.0
        f1 = (2 * p * r / (p + r)) if (p + r) > 0 else 0.0
        precisions.append(p); recalls.append(r); f1s.append(f1); supports.append(tpc + fnc)

    if len(class_indices) > 0:
        precision_macro = sum(precisions) / len(precisions)
        recall_macro = sum(recalls) / len(recalls)
        f1_macro = sum(f1s) / len(f1s)
    else:
        precision_macro = recall_macro = f1_macro = 0.0

    total_support = sum(supports)
    if total_support > 0:
        precision_weighted = sum(p * s for p, s in zip(precisions, supports)) / total_support
        recall_weighted = sum(r * s for r, s in zip(recalls, supports)) / total_support
        f1_weighted = sum(f * s for f, s in zip(f1s, supports)) / total_support
    else:
        precision_weighted = recall_weighted = f1_weighted = 0.0

    metrics = {
        'precision_macro': precision_macro,
        'recall_macro': recall_macro,
        'f1_macro': f1_macro,
        'precision_weighted': precision_weighted,
        'recall_weighted': recall_weighted,
        'f1_weighted': f1_weighted,
    }

    return acc, class_stats, inference_results, metrics


def main():
    parser = argparse.ArgumentParser(description="Validate a saved classifier on 6-channel .npy dataset (VAL ONLY, no loss)")
    # Data (Mode A or Mode B)
    parser.add_argument("data_path", nargs="?", type=str,
                        help="Dataset root that contains 'val/' (Mode A).")
    parser.add_argument("--val_data_paths_file", type=str, default=None,
                        help="Text file listing VAL dataset roots (one per line). (Mode B)")

    # Model checkpoint
    parser.add_argument("--model_path", type=str, required=True,
                        help="Path to the saved checkpoint (.pth) with model_state_dict etc.")
    parser.add_argument("-o", "--output_path", default="./val_only_output", type=str,
                        help="Path to write logs / optional JSON.")

    # Dataloader
    parser.add_argument("--batch_size", type=int, default=32, help="Batch size for validation")
    parser.add_argument("--num_workers", type=int, default=8, help="#workers for dataloader")

    # Arch (used if checkpoint lacks dims/depths)
    parser.add_argument("--depths", type=int, nargs='+', default=[3, 3, 27, 3],
                        help="ConvNeXt stage depths, used if ckpt missing this info")
    parser.add_argument("--dims", type=int, nargs='+', default=[192, 384, 768, 1536],
        help="ConvNeXt dims, used if ckpt missing this info")

    # Output toggles
    parser.add_argument("--save_val_results_json", action="store_true",
                        help="If set, save inference_results + metrics JSON to output_path.")

    args = parser.parse_args()
    os.makedirs(args.output_path, exist_ok=True)
    log_file = os.path.join(args.output_path, "val_log_6ch.txt")
    if os.path.exists(log_file):
        os.remove(log_file)

    # ---------- Build VAL loader ----------
    if (args.data_path is None) == (args.val_data_paths_file is None):
        parser.error("Provide exactly one of: positional data_path (Mode A) OR --val_data_paths_file (Mode B).")

    if args.data_path is not None:
        val_source = os.path.abspath(os.path.expanduser(args.data_path))
    else:
        val_roots = _read_paths_file(args.val_data_paths_file)
        if not val_roots:
            parser.error(f"--val_data_paths_file is empty or unreadable: {args.val_data_paths_file}")
        val_source = val_roots

    try:
        val_loader, genotype_map = get_data_loader(
            data_dir=val_source,
            dataset_type="val",
            batch_size=args.batch_size,
            num_workers=args.num_workers,
            shuffle=False,
            return_paths=True
        )
    except Exception as e:
        print_and_log(f"FATAL: could not create VAL dataloader with return_paths=True\nError: {e}", log_file)
        sys.exit(1)

    if not genotype_map:
        print_and_log("Error: genotype_map is empty (from dataloader).", log_file)
        sys.exit(1)
    num_classes = len(genotype_map)
    print_and_log(f"Using device: {device}", log_file)
    print_and_log(f"Detected {num_classes} classes from dataloader.", log_file)

    # ---------- Load checkpoint ----------
    ckpt_path = os.path.abspath(os.path.expanduser(args.model_path))
    if not os.path.isfile(ckpt_path):
        print_and_log(f"Checkpoint not found: {ckpt_path}", log_file)
        sys.exit(1)

    ckpt = torch.load(ckpt_path, map_location="cpu")
    ckpt_in_channels = ckpt.get('in_channels', 6)
    ckpt_genotype_map = ckpt.get('genotype_map', None)
    ckpt_depths = ckpt.get('depths', None)
    ckpt_dims = ckpt.get('dims', None)

    # Prefer checkpoint's genotype_map if provided
    if ckpt_genotype_map:
        if ckpt_genotype_map != genotype_map:
            print_and_log("Warning: checkpoint genotype_map differs from dataloader map. Using checkpoint's map.", log_file)
        genotype_map = ckpt_genotype_map
        num_classes = len(genotype_map)

    depths = ckpt_depths if ckpt_depths is not None else args.depths
    dims = ckpt_dims if ckpt_dims is not None else args.dims

    # ---------- Build model and load weights ----------
    model = ConvNeXtCBAMClassifier(in_channels=ckpt_in_channels, class_num=num_classes,
                                   depths=depths, dims=dims).to(device)

    state = ckpt.get('model_state_dict', ckpt)  # allow raw state_dict file
    missing, unexpected = model.load_state_dict(state, strict=False)
    if missing:
        print_and_log(f"Note: missing keys when loading state_dict: {missing}", log_file)
    if unexpected:
        print_and_log(f"Note: unexpected keys when loading state_dict: {unexpected}", log_file)

    # ---------- Evaluate (no loss) ----------
    val_acc, class_stats, inference_results, metrics = evaluate_model(
        model, val_loader, genotype_map, log_file
    )

    # ---------- Report ----------
    print_and_log("\n=== Validation Summary (no loss) ===", log_file)
    print_and_log(f"Val Acc : {val_acc:.2f}%", log_file)
    print_and_log(
        f"Precision (macro): {metrics['precision_macro']*100:.2f}% | "
        f"Recall (macro): {metrics['recall_macro']*100:.2f}% | "
        f"F1 (macro): {metrics['f1_macro']*100:.2f}%", log_file
    )
    print_and_log(
        f"Precision (weighted): {metrics['precision_weighted']*100:.2f}% | "
        f"Recall (weighted): {metrics['recall_weighted']*100:.2f}% | "
        f"F1 (weighted): {metrics['f1_weighted']*100:.2f}%", log_file
    )

    if class_stats:
        print_and_log("\nPer-class Accuracy:", log_file)
        for cname, stats in sorted(class_stats.items(), key=lambda kv: kv[1]['idx']):
            print_and_log(
                f"  {cname} (idx {stats['idx']}): {stats['acc']:.2f}% "
                f"({stats['correct']}/{stats['total']})", log_file
            )

    if args.save_val_results_json:
        out_json = os.path.join(args.output_path, "validation_results.json")
        payload = {
            'val_acc': val_acc,
            'metrics': metrics,
            'class_stats': class_stats,
            'inference_results': inference_results,
            'genotype_map': genotype_map,
        }
        try:
            with open(out_json, 'w') as f:
                json.dump(payload, f, indent=2)
            print_and_log(f"\nSaved validation JSON to: {out_json}", log_file)
        except Exception as e:
            print_and_log(f"Error saving validation_results.json: {e}", log_file)


if __name__ == "__main__":
    main()
