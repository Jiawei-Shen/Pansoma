#!/usr/bin/env python3
"""Render five-channel pileup tensors as one-row 21 bp locus-centered panels."""

import argparse
import json
import os
import re
import sys
from typing import Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
DEFAULT_OUTPUT_ROOT = os.path.join(PROJECT_DIR, "tensor_figures", "21bp")
FLANK_BP = 10
WINDOW_BP = FLANK_BP * 2 + 1

BASE_LABELS = {
    0: "Padding",
    1: "Deletion",
    5: "N",
    20: "A",
    30: "C",
    50: "G",
    70: "T",
    90: "Insertion",
}
BASE_COLORS = {
    0: "#ffffff",
    1: "#777777",
    5: "#a6a6a6",
    20: "#4daf4a",
    30: "#377eb8",
    50: "#ffb000",
    70: "#e41a1c",
    90: "#984ea3",
}
CIGAR_LABELS = {
    0: "Padding",
    10: "M",
    20: "N",
    30: "S",
    40: "I",
    50: "D",
    60: "H",
    70: "P",
    80: "=",
    90: "X",
}
CIGAR_COLORS = {
    0: "#ffffff",
    10: "#5ba66c",
    20: "#9c9c9c",
    30: "#efb04c",
    40: "#9655c7",
    50: "#d9534f",
    60: "#5e6472",
    70: "#a8a8a8",
    80: "#3788c8",
    90: "#d62828",
}
MISMATCH_LABELS = {-1: "Padding", 0: "Match", 1: "Mismatch"}
MISMATCH_COLORS = {-1: "#ffffff", 0: "#d9ead3", 1: "#cc0000"}
QUALITY_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "quality_blue",
    ["#e8f1f8", "#c3d9ea", "#82b0d0", "#3c78a8", "#083d67"],
)
QUALITY_COLOR_SCALE_MAX = 60
BASE_QUALITY_DISPLAY_MAX = 40
MAPPING_QUALITY_DISPLAY_MAX = 60
PANEL_LABELS = [
    "Bases",
    "Base quality",
    "Mapping quality",
    "Mismatch flag",
    "CIGAR",
]


def configure_style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": "#ffffff",
            "axes.facecolor": "#fbfcfd",
            "axes.edgecolor": "#c8d0d9",
            "axes.linewidth": 0.8,
            "axes.labelcolor": "#101828",
            "font.family": "DejaVu Sans",
            "font.size": 17,
            "xtick.color": "#475467",
            "ytick.color": "#475467",
        }
    )


def load_npy(npy_path: str) -> np.ndarray:
    array = np.load(npy_path, mmap_mode="r")
    if array.ndim not in (3, 4):
        raise ValueError(
            f"expected (N, C, H, W) or (C, H, W) NPY data, got shape {array.shape}"
        )
    return array


def load_single_tensor(array: np.ndarray, sample_index: int) -> Tuple[np.ndarray, Optional[int]]:
    if array.ndim == 4:
        if not 0 <= sample_index < array.shape[0]:
            raise IndexError(
                f"sample index {sample_index} is outside shard range 0..{array.shape[0] - 1}"
            )
        return np.asarray(array[sample_index]), sample_index
    if sample_index != 0:
        raise IndexError("--sample-index can only be 0 for a single tensor input")
    return np.asarray(array), None


def infer_shard_index(npy_path: str) -> int:
    match = re.search(r"_(\d+)_data\.npy$", os.path.basename(npy_path))
    if not match:
        raise ValueError(
            "could not infer shard index; expected an input name such as *_00000_data.npy"
        )
    return int(match.group(1))


def classified_sample_indices(
    summary_path: str, classification: str, shard_index: int, sample_count: int
) -> list[int]:
    indices = []
    with open(summary_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"invalid JSON in summary line {line_number}: {exc}") from exc
            if (
                record.get("classification") == classification
                and record.get("shard_index") == shard_index
            ):
                index = record.get("index_within_shard")
                if not isinstance(index, int) or not 0 <= index < sample_count:
                    raise ValueError(
                        f"invalid index_within_shard in summary line {line_number}: {index}"
                    )
                indices.append(index)
    return indices


def infer_marker_column(tensor: np.ndarray) -> int:
    marked_columns = np.flatnonzero(np.any(tensor[2] == 5, axis=0))
    return int(marked_columns[0]) if marked_columns.size else tensor.shape[2] // 2


def trim_padding_rows(tensor: np.ndarray, show_all_rows: bool) -> np.ndarray:
    if show_all_rows:
        return tensor
    non_padding_reads = np.flatnonzero(np.any(tensor[0, 1:, :] != 0, axis=1))
    final_row = int(non_padding_reads[-1] + 2) if non_padding_reads.size else 1
    return tensor[:, :final_row, :]


def locus_window(tensor: np.ndarray, marker_column: int) -> np.ndarray:
    fills = (0, 0, -1, -1, 0)
    window = np.empty((tensor.shape[0], tensor.shape[1], WINDOW_BP), dtype=tensor.dtype)
    for channel, fill in enumerate(fills):
        window[channel].fill(fill)
    source_start = max(0, marker_column - FLANK_BP)
    source_end = min(tensor.shape[2], marker_column + FLANK_BP + 1)
    destination_start = source_start - (marker_column - FLANK_BP)
    destination_end = destination_start + (source_end - source_start)
    window[:, :, destination_start:destination_end] = tensor[:, :, source_start:source_end]
    return window


def discrete_cmap(
    labels: dict[int, str], colors: dict[int, str]
) -> Tuple[mcolors.ListedColormap, mcolors.BoundaryNorm]:
    values = sorted(labels)
    boundaries = [values[0] - 0.5]
    boundaries.extend((left + right) / 2.0 for left, right in zip(values, values[1:]))
    boundaries.append(values[-1] + 0.5)
    cmap = mcolors.ListedColormap([colors[value] for value in values])
    return cmap, mcolors.BoundaryNorm(boundaries, cmap.N)


def legend_handles(
    values: np.ndarray, labels: dict[int, str], colors: dict[int, str]
) -> list[Line2D]:
    return [
        Line2D(
            [0],
            [0],
            marker="s",
            linestyle="",
            label=labels[value],
            markerfacecolor=colors[value],
            markeredgecolor="#555555",
            markersize=9,
        )
        for value in sorted(set(int(x) for x in np.unique(values)) & set(labels))
    ]


def add_discrete_legend(
    ax: plt.Axes, values: np.ndarray, labels: dict[int, str], colors: dict[int, str]
) -> None:
    handles = legend_handles(values, labels, colors)
    if handles:
        compact = labels is BASE_LABELS
        ax.legend(
            handles=handles,
            loc="lower center",
            bbox_to_anchor=(0.5, 1.075),
            ncol=len(handles),
            frameon=False,
            fontsize=12 if compact else 13,
            handletextpad=0.12 if compact else 0.2,
            columnspacing=0.30 if compact else 0.50,
            borderaxespad=0,
        )


def draw_discrete_panel(
    ax: plt.Axes, values: np.ndarray, labels: dict[int, str], colors: dict[int, str]
) -> None:
    cmap, norm = discrete_cmap(labels, colors)
    ax.imshow(values, cmap=cmap, norm=norm, interpolation="nearest", aspect="auto")
    add_discrete_legend(ax, values, labels, colors)


def draw_continuous_panel(
    fig: plt.Figure,
    ax: plt.Axes,
    values: np.ndarray,
    displayed_max: int,
) -> None:
    cmap = QUALITY_CMAP.copy()
    cmap.set_bad("#ffffff")
    image = ax.imshow(
        values,
        cmap=cmap,
        vmin=0,
        vmax=QUALITY_COLOR_SCALE_MAX,
        interpolation="nearest",
        aspect="auto",
    )
    cax = ax.inset_axes((0.05, 1.09, 0.90, 0.045))
    colorbar = fig.colorbar(image, cax=cax, orientation="horizontal")
    colorbar.ax.set_xlim(0, displayed_max)
    tick_step = 10 if displayed_max <= 40 else 20
    colorbar.set_ticks(list(range(0, displayed_max + 1, tick_step)))
    colorbar.ax.tick_params(labelsize=13, length=2, pad=2)
    colorbar.outline.set_edgecolor("#c8d0d9")


def visualize_tensor(
    tensor: np.ndarray,
    output: str,
    title: str,
    show_all_rows: bool,
    marker_column: int,
) -> int:
    if tensor.ndim != 3 or tensor.shape[0] != 5:
        raise ValueError(f"expected a 5-channel (5, H, W) tensor, got shape {tensor.shape}")

    configure_style()
    view = trim_padding_rows(locus_window(tensor, marker_column), show_all_rows)
    rows = view.shape[1]
    fig, axes = plt.subplots(1, 5, figsize=(18, 18 / 3.5), sharey=True)
    fig.subplots_adjust(left=0.06, right=0.985, top=0.73, bottom=0.19, wspace=0.08)

    draw_discrete_panel(axes[0], view[0], BASE_LABELS, BASE_COLORS)
    populated_reads = view[0] != 0
    populated_reads[0, :] = False
    base_quality = np.ma.masked_where(~populated_reads, view[1])
    mapping_quality = np.ma.masked_where((~populated_reads) | (view[3] < 0), view[3])
    draw_continuous_panel(fig, axes[1], base_quality, BASE_QUALITY_DISPLAY_MAX)
    draw_continuous_panel(fig, axes[2], mapping_quality, MAPPING_QUALITY_DISPLAY_MAX)
    mismatch_values = np.where(view[2] == 5, 1, view[2])
    draw_discrete_panel(axes[3], mismatch_values, MISMATCH_LABELS, MISMATCH_COLORS)
    draw_discrete_panel(axes[4], view[4], CIGAR_LABELS, CIGAR_COLORS)

    for panel_index, ax in enumerate(axes):
        ax.axvline(FLANK_BP, color="#101828", linewidth=0.9, linestyle=(0, (3, 3)))
        ax.set_xticks([0, FLANK_BP, WINDOW_BP - 1], ["-10", "0", "+10"])
        ax.set_xlabel(PANEL_LABELS[panel_index], fontsize=21, fontweight="semibold", labelpad=12)
        ax.tick_params(axis="x", labelsize=14)
        ax.grid(axis="x", color="#e9edf2", linewidth=0.5)
        if panel_index != 0:
            ax.tick_params(axis="y", left=False, labelleft=False)
    axes[0].set_ylabel("Read row", fontsize=19)
    axes[0].tick_params(axis="y", labelsize=14)
    fig.suptitle(title, fontsize=24, fontweight="semibold", color="#101828", y=0.97)

    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    fig.savefig(output, dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Render five pileup channels side by side over a 21 bp variant-centered window."
    )
    parser.add_argument("npy_path", help="Input .npy path: (N,5,H,W) shard or (5,H,W) tensor.")
    parser.add_argument("-i", "--sample-index", type=int, default=0)
    parser.add_argument("-o", "--output", help="Output PNG for a single sample.")
    parser.add_argument("--output-dir", help="Output directory override.")
    parser.add_argument("--all-samples", action="store_true")
    parser.add_argument("--classification", choices=("true", "false"))
    parser.add_argument(
        "--summary-path",
        help="Classified NDJSON path. Default: variant_summary_classified.ndjson beside input.",
    )
    parser.add_argument("--max-samples", type=int)
    parser.add_argument("--show-all-rows", action="store_true")
    parser.add_argument("--marker-column", type=int)
    parser.add_argument("--title")
    args = parser.parse_args()

    path = os.path.expanduser(args.npy_path)
    if not os.path.isfile(path):
        sys.exit(f"Error: input NPY does not exist: {path}")
    try:
        array = load_npy(path)
    except (OSError, ValueError) as exc:
        sys.exit(f"Error: {exc}")

    basename = os.path.splitext(os.path.basename(path))[0]
    default_output_dir = os.path.join(DEFAULT_OUTPUT_ROOT, basename)
    if args.classification:
        default_output_dir = os.path.join(default_output_dir, f"{args.classification}_variants")
    output_dir = os.path.expanduser(args.output_dir or default_output_dir)
    if args.classification and args.all_samples:
        parser.error("--classification cannot be combined with --all-samples")
    if (args.classification or args.all_samples) and args.output:
        parser.error("--output cannot be combined with batch rendering")
    if args.summary_path and not args.classification:
        parser.error("--summary-path requires --classification")
    if args.max_samples is not None and not (args.classification or args.all_samples):
        parser.error("--max-samples requires --classification or --all-samples")
    if args.max_samples is not None and args.max_samples <= 0:
        parser.error("--max-samples must be greater than zero")

    if args.classification:
        if array.ndim != 4:
            parser.error("--classification requires an (N,C,H,W) shard")
        summary_path = os.path.expanduser(
            args.summary_path
            or os.path.join(os.path.dirname(path), "variant_summary_classified.ndjson")
        )
        if not os.path.isfile(summary_path):
            parser.error(f"classified summary does not exist: {summary_path}")
        try:
            sample_indices = classified_sample_indices(
                summary_path, args.classification, infer_shard_index(path), array.shape[0]
            )
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
        if args.max_samples is not None:
            sample_indices = sample_indices[:args.max_samples]
        if not sample_indices:
            parser.error(f"no {args.classification} samples found for this shard")
    elif args.all_samples:
        if array.ndim != 4:
            parser.error("--all-samples requires an (N,C,H,W) shard")
        stop = array.shape[0]
        if args.max_samples is not None:
            stop = min(stop, args.sample_index + args.max_samples)
        sample_indices = list(range(args.sample_index, stop))
    else:
        sample_indices = [args.sample_index]

    batch = args.classification is not None or args.all_samples
    if batch:
        print(f"Rendering {len(sample_indices)} figures to: {os.path.abspath(output_dir)}")
    rendered = 0
    last_output = ""
    for sample_index in sample_indices:
        try:
            tensor, tensor_index = load_single_tensor(array, sample_index)
        except IndexError as exc:
            sys.exit(f"Error: {exc}")
        suffix = f"_sample_{sample_index:05d}_21bp" if tensor_index is not None else "_21bp"
        output = os.path.expanduser(
            args.output or os.path.join(output_dir, f"{basename}{suffix}.png")
        )
        classification_text = f" / {args.classification}" if args.classification else ""
        title = args.title or f"{basename}{classification_text} / sample {sample_index:05d} / 21 bp"
        marker_column = (
            infer_marker_column(tensor) if args.marker_column is None else args.marker_column
        )
        if not 0 <= marker_column < tensor.shape[2]:
            sys.exit(f"Error: marker column {marker_column} outside tensor width {tensor.shape[2]}")
        try:
            visualize_tensor(tensor, output, title, args.show_all_rows, marker_column)
        except ValueError as exc:
            sys.exit(f"Error: {exc}")
        rendered += 1
        last_output = os.path.abspath(output)
        if batch and (rendered == 1 or rendered % 100 == 0 or rendered == len(sample_indices)):
            print(f"Rendered {rendered}/{len(sample_indices)}")

    print(f"Input shape: {array.shape}; dtype: {array.dtype}; window: {WINDOW_BP} bp")
    print(f"Rendered figures: {rendered}")
    print(f"Saved in: {os.path.abspath(output_dir)}" if batch else f"Saved: {last_output}")


if __name__ == "__main__":
    main()
