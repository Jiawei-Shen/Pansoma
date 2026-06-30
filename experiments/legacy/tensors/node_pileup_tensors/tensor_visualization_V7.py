#!/usr/bin/env python3
import argparse
import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")  # safe for cluster / no GUI
import matplotlib.pyplot as plt


def log(msg: str) -> None:
    print(msg, flush=True)


def visualize_sample(sample, out_prefix):
    """
    sample:
      - (C, H, W)  → multiple channels
      - (H, W)     → single channel

    out_prefix: path without extension
    """
    if sample.ndim == 3:
        # (C, H, W)
        C, H, W = sample.shape
        # layout: 1xC or 2x3 etc.
        # keep it simple: up to 3 columns per row
        ncols = min(3, C)
        nrows = (C + ncols - 1) // ncols

        fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows))
        if not isinstance(axes, np.ndarray):
            axes = np.array([[axes]])
        axes = axes.reshape(nrows, ncols)

        for c in range(C):
            r = c // ncols
            k = c % ncols
            ax = axes[r, k]
            im = ax.imshow(sample[c], aspect="auto", interpolation="nearest")
            ax.set_title(f"Channel {c}")
            ax.set_xlabel("W")
            ax.set_ylabel("H")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        # hide any unused axes
        for c in range(C, nrows * ncols):
            r = c // ncols
            k = c % ncols
            axes[r, k].axis("off")

        fig.suptitle(f"Sample preview: shape={sample.shape}", fontsize=12)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    elif sample.ndim == 2:
        # (H, W)
        H, W = sample.shape
        fig, ax = plt.subplots(figsize=(6, 4))
        im = ax.imshow(sample, aspect="auto", interpolation="nearest")
        ax.set_title(f"Sample preview: shape={sample.shape}")
        ax.set_xlabel("W")
        ax.set_ylabel("H")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        fig.tight_layout()
    else:
        raise ValueError(f"Expected sample ndim 2 or 3, got {sample.ndim} with shape {sample.shape}")

    out_png = out_prefix + ".png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    log(f"  Saved {out_png}")


def main():
    parser = argparse.ArgumentParser(
        description="Visualize head K tensors from a single .npy file."
    )
    parser.add_argument(
        "npy_path",
        help="Path to the .npy file (expected shape (N, C, H, W) or (N, H, W)).",
    )
    parser.add_argument(
        "--head",
        type=int,
        default=5,
        help="Number of first tensors (samples) to visualize.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./npy_head_vis",
        help="Where to save PNGs.",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.npy_path):
        sys.exit(f"Error: file not found: {args.npy_path}")

    os.makedirs(args.output_dir, exist_ok=True)

    log(f"Loading {args.npy_path} ...")
    arr = np.load(args.npy_path)

    if arr.ndim not in (3, 4):
        sys.exit(f"Error: expected array ndim 3 or 4 (N,C,H,W or N,H,W), but got {arr.ndim} with shape {arr.shape}")

    N = arr.shape[0]
    k = min(args.head, N)

    log(f"Array shape: {arr.shape}")
    log(f"Total samples (N): {N}")
    log(f"Visualizing first {k} samples...")

    base_name = os.path.splitext(os.path.basename(args.npy_path))[0]

    for i in range(k):
        sample = arr[i]  # (C,H,W) or (H,W)
        out_prefix = os.path.join(args.output_dir, f"{base_name}_head{i:04d}")
        visualize_sample(sample, out_prefix)

    log("Done.")


if __name__ == "__main__":
    main()
