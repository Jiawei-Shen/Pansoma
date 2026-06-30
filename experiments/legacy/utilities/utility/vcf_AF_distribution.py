#!/usr/bin/env python3
"""
AF distribution from VCF/BCF using pysam, plotted over 8 predefined AF bins.

Predefined AF bins (super-bins):
  0–1e-6, 1e-6–1e-5, 1e-5–1e-4, 1e-4–1e-3,
  1e-3–1e-2, 1e-2–0.1, 0.1–0.5, 0.5–1.0

Features:
- Streams AF from INFO (default: AF)
- Handles multi-allelic sites (one AF per ALT)
- Optional FILTER=PASS
- Optional region fetch (requires .tbi/.csi)
- Variant-type filter: snp / indel / all
- Streaming counting (low-mem) or accumulate mode
- Optional progress bar
- Optional threads (htslib) if supported

Outputs (based on the 8 super-bins):
- PNG bar chart
- CSV with AF_Class, bin_left, bin_right, count, proportion
- NPY arrays (edges + counts) for reproducibility

Author: revised to use predefined 8 AF bins.
"""

import argparse
import math
from pathlib import Path

import numpy as np
import pysam
import matplotlib.pyplot as plt


# ──────────────────────────────────────────────────────────────────────────────
# Predefined AF super-bins and labels
# (edges are used with right-closed intervals: (edge[i], edge[i+1]]; 0 maps into first bin)
AF_BINS = np.array([0.0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5, 1.0], dtype=float)
AF_LABELS = [
    "0–1e-6",
    "1e-6–1e-5",
    "1e-5–1e-4",
    "1e-4–1e-3",
    "1e-3–1e-2",
    "1e-2–0.1",
    "0.1–0.5",
    "0.5–1.0",
]
# ──────────────────────────────────────────────────────────────────────────────


def parse_args():
    p = argparse.ArgumentParser(
        description="Compute AF distribution from VCF/BCF using 8 predefined AF bins"
    )
    p.add_argument("vcf", help="Path to .vcf.gz/.vcf.bgz or .bcf")
    p.add_argument("--field", default="AF", help="INFO field to use (default: AF)")
    p.add_argument("--pass-only", action="store_true", help="Keep only FILTER=PASS")
    p.add_argument("--region", default=None,
                   help="Region like 'chr1:1-100000' (requires index)")
    p.add_argument("--variant-type", choices=["all", "snp", "indel"], default="all",
                   help="Variant type to include (default: all)")
    p.add_argument("--max-af", type=float, default=None,
                   help="Keep AF <= value (e.g., 0.05 to focus on rare variants)")
    # The following two are ignored for super-bins; kept for CLI compatibility.
    p.add_argument("--bins", type=int, default=100, help="(Ignored) histogram bins")
    p.add_argument("--logx", action="store_true",
                   help="(Ignored) log10 histogram; super-bins override")
    p.add_argument("--stream-hist", action="store_true",
                   help="Stream counting into the 8 super-bins (low memory)")
    p.add_argument("--out-prefix", default="af_distribution",
                   help="Output file prefix (default: af_distribution)")
    p.add_argument("--show-progress", action="store_true",
                   help="Display a tqdm progress bar while reading")
    p.add_argument("--threads", type=int, default=0,
                   help="htslib bgzf threads if supported (0=disabled)")
    p.add_argument("--save-values", action="store_true",
                   help="Write raw AF values to CSV (accumulate mode only; large!)")
    return p.parse_args()


def is_pass(rec) -> bool:
    """Accept PASS if filter set is empty or contains 'PASS'."""
    keys = rec.filter.keys()
    return (len(keys) == 0) or ("PASS" in keys)


def is_snp_record(rec) -> bool:
    """True if REF and all ALT alleles are length 1."""
    if rec.alts is None:
        return False
    if len(rec.ref) != 1:
        return False
    return all(len(alt) == 1 for alt in rec.alts)


def is_indel_record(rec) -> bool:
    """True if any ALT has length different from REF (insertions/deletions)."""
    if rec.alts is None:
        return False
    rlen = len(rec.ref)
    return any(len(alt) != rlen for alt in rec.alts)


def af_iterator(vf: pysam.VariantFile, field: str, pass_only: bool,
                variant_type: str, region: str | None, show_progress: bool):
    """
    Yield AF values (float) for each ALT allele of records passing filters.
    """
    it = vf.fetch(region) if region else vf.fetch()
    if show_progress:
        try:
            from tqdm import tqdm
            it = tqdm(it, unit="variants", mininterval=5, desc="Processing VCF")
        except Exception:
            pass

    for rec in it:
        if pass_only and not is_pass(rec):
            continue
        if variant_type == "snp" and not is_snp_record(rec):
            continue
        if variant_type == "indel" and not is_indel_record(rec):
            continue

        info_val = rec.info.get(field, None)
        if info_val is None:
            continue

        vals = info_val if isinstance(info_val, (tuple, list)) else (info_val,)
        for v in vals:
            if v is None:
                continue
            try:
                fv = float(v)
            except Exception:
                continue
            if math.isfinite(fv):
                yield fv


def bin_counts_superbins(values_iter, max_af=None):
    """
    Count AF values into the 8 super-bins.
    Uses right-closed intervals: (edge[i], edge[i+1]], and maps 0 into first bin.
    """
    counts = np.zeros(len(AF_LABELS), dtype=np.int64)
    edges = AF_BINS

    for fv in values_iter:
        if fv < 0:
            continue
        if max_af is not None and fv > max_af:
            continue
        # clip tiny zeros (e.g., exactly 0) into first bin
        if fv == 0.0:
            idx = 0
        else:
            # np.digitize with right=True gives (edge[i], edge[i+1]]
            idx = np.digitize(fv, edges, right=True) - 1
        if 0 <= idx < counts.size:
            counts[idx] += 1
        # values > edges[-1] or < edges[0] are ignored

    return counts


def save_superbin_artifacts(out_prefix: Path, counts: np.ndarray, title_note: str):
    """
    Save CSV, NPY, and PNG using the 8 super-bins.
    """
    total = int(counts.sum())
    proportions = counts / total if total > 0 else np.zeros_like(counts, dtype=float)

    # CSV: AF_Class, bin_left, bin_right, count, proportion
    csvp = out_prefix.with_suffix(".histogram.csv")
    with open(csvp, "w") as f:
        f.write("AF_Class,bin_left,bin_right,count,proportion\n")
        for i, lab in enumerate(AF_LABELS):
            left = AF_BINS[i]
            right = AF_BINS[i + 1]
            f.write(f"{lab},{left},{right},{int(counts[i])},{proportions[i]:.6f}\n")

    # Save NPY for reproducibility
    np.save(out_prefix.with_suffix(".hist_edges.npy"), AF_BINS)
    np.save(out_prefix.with_suffix(".hist_counts.npy"), counts)

    # Plot PNG: bar per super-bin
    # x positions and widths
    mids = (AF_BINS[:-1] + AF_BINS[1:]) / 2.0
    widths = AF_BINS[1:] - AF_BINS[:-1]

    plt.figure()
    plt.bar(mids, counts, width=widths, align="center")
    plt.xticks(mids, AF_LABELS, rotation=35, ha="right")
    plt.xlabel("Allele Frequency (AF) — 8 predefined bins")
    plt.ylabel("Count (ALT alleles)")
    ttl = "AF Distribution (8 bins)"
    if title_note:
        ttl += f" — {title_note}"
    plt.title(ttl)
    plt.tight_layout()
    fig_path = out_prefix.with_suffix(".png")
    plt.savefig(fig_path, dpi=180)
    plt.close()
    return csvp, fig_path


def main():
    args = parse_args()

    if args.logx or (args.bins and args.bins != 100):
        print("[info] Note: --logx/--bins are ignored; plotting uses the 8 predefined AF bins.")

    vcf_path = args.vcf
    out_prefix = Path(args.out_prefix)

    vf = pysam.VariantFile(vcf_path)
    if args.threads and args.threads > 0:
        try:
            vf.set_threads(args.threads)
            print(f"[info] Using htslib threads: {args.threads}")
        except Exception:
            print("[warn] VariantFile.set_threads() not available; proceeding single-threaded.")

    # ── STREAMING MODE: count directly into 8 bins ────────────────────────────
    if args.stream_hist:
        counts = bin_counts_superbins(
            af_iterator(vf, args.field, args.pass_only, args.variant_type, args.region, args.show_progress),
            max_af=args.max_af
        )
        note = f"type={args.variant_type}"
        csvp, figp = save_superbin_artifacts(out_prefix, counts, title_note=note)
        print(f"[ok] Saved (streaming) -> {csvp} and {figp}")
        return

    # ── ACCUMULATE MODE: (optionally save raw values), then bin ───────────────
    af_vals = []
    push = af_vals.append
    for fv in af_iterator(vf, args.field, args.pass_only, args.variant_type, args.region, args.show_progress):
        if args.max_af is not None and fv > args.max_af:
            continue
        push(fv)

    if not af_vals:
        print("[warn] No AF values collected. Check field/filters/region/variant-type.")
        return

    af = np.fromiter(af_vals, dtype=np.float64)
    if args.save_values:
        raw_csv = out_prefix.with_suffix(".af_values.csv")
        np.savetxt(raw_csv, af, delimiter=",", fmt="%.8g")
        print(f"[ok] Saved AF values -> {raw_csv}")

    counts = bin_counts_superbins(af, max_af=None)  # already applied max-af above
    note = f"type={args.variant_type}, n={af.size}"
    csvp, figp = save_superbin_artifacts(out_prefix, counts, title_note=note)
    print(f"[ok] Saved -> {csvp} and {figp}")


if __name__ == "__main__":
    main()
