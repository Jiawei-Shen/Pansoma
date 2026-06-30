#!/usr/bin/env python3
"""Visualize one 5-channel pileup tensor stored in an NPY file or shard."""

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
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import numpy as np


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
DEFAULT_OUTPUT_ROOT = os.path.join(PROJECT_DIR, "tensor_figures")
BASE_LABELS = {
    0: "Padding",
    1: "Deletion (*)",
    5: "N",
    20: "A",
    30: "C",
    50: "G",
    70: "T",
    90: "Insertion (I)",
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
CHANNEL_TITLES = [
    "Bases",
    "Base quality",
    "Mapping quality",
    "Mismatch flag",
    "CIGAR operation",
]


def configure_style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": "#ffffff",
            "axes.facecolor": "#fbfcfd",
            "axes.edgecolor": "#c8d0d9",
            "axes.linewidth": 0.8,
            "axes.labelcolor": "#344054",
            "axes.titlecolor": "#101828",
            "axes.titleweight": "semibold",
            "font.family": "DejaVu Sans",
            "font.size": 16,
            "xtick.color": "#475467",
            "ytick.color": "#475467",
            "xtick.labelsize": 15,
            "ytick.labelsize": 15,
            "xtick.major.size": 4,
            "ytick.major.size": 4,
        }
    )


def discrete_cmap(
    labels: dict[int, str], colors: dict[int, str]
) -> Tuple[mcolors.ListedColormap, mcolors.BoundaryNorm, Sequence[int]]:
    values = sorted(labels)
    boundaries = [values[0] - 0.5]
    for left, right in zip(values, values[1:]):
        boundaries.append((left + right) / 2.0)
    boundaries.append(values[-1] + 0.5)
    cmap = mcolors.ListedColormap([colors[value] for value in values])
    norm = mcolors.BoundaryNorm(boundaries, cmap.N)
    return cmap, norm, values


def legend_handles(
    present_values: np.ndarray, labels: dict[int, str], colors: dict[int, str]
) -> list[Line2D]:
    handles = []
    for value in sorted(set(int(x) for x in present_values) & set(labels)):
        handles.append(
            Line2D(
                [0],
                [0],
                marker="s",
                linestyle="",
                label=labels[value],
                markerfacecolor=colors[value],
                markeredgecolor="#555555",
                markersize=7,
            )
        )
    return handles


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
    else:
        if sample_index != 0:
            raise IndexError("--sample-index can only be 0 when the input is already one tensor")
        return np.asarray(array), None


def infer_shard_index(npy_path: str) -> int:
    match = re.search(r"_(\d+)_data\.npy$", os.path.basename(npy_path))
    if not match:
        raise ValueError(
            "could not infer shard index from input name; expected a name such as *_00000_data.npy"
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


def trim_padding_rows(tensor: np.ndarray, show_all_rows: bool) -> np.ndarray:
    if show_all_rows:
        return tensor
    non_padding_reads = np.flatnonzero(np.any(tensor[0, 1:, :] != 0, axis=1))
    last_row = int(non_padding_reads[-1] + 2) if non_padding_reads.size else 1
    return tensor[:, :last_row, :]


def infer_marker_column(tensor: np.ndarray) -> int:
    marked_columns = np.flatnonzero(np.any(tensor[2] == 5, axis=0))
    if marked_columns.size:
        return int(marked_columns[0])
    return tensor.shape[2] // 2


def plot_legend(
    ax: plt.Axes, values: np.ndarray, labels: dict[int, str], colors: dict[int, str]
) -> None:
    handles = legend_handles(np.unique(values), labels, colors)
    ax.axis("off")
    if handles:
        ax.legend(
            handles=handles,
            loc="center left",
            borderaxespad=0,
            fontsize=15,
            frameon=False,
            handletextpad=0.6,
            labelspacing=0.45,
        )


def plot_discrete_track(
    ax: plt.Axes, side_ax: plt.Axes, values: np.ndarray, labels: dict[int, str],
    colors: dict[int, str]
) -> None:
    cmap, norm, _ = discrete_cmap(labels, colors)
    ax.imshow(values, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
    plot_legend(side_ax, values, labels, colors)


def plot_continuous_track(
    fig: plt.Figure,
    ax: plt.Axes,
    side_ax: plt.Axes,
    values: np.ndarray,
    cmap: str | mcolors.Colormap,
    vmin: int,
    color_scale_max: int,
    displayed_max: int,
) -> None:
    colormap = plt.get_cmap(cmap).copy()
    colormap.set_bad("#ffffff")
    image = ax.imshow(
        values, cmap=colormap, vmin=vmin, vmax=color_scale_max,
        aspect="auto", interpolation="nearest"
    )
    side_ax.axis("off")
    cax = side_ax.inset_axes((0.03, 0.12, 0.10, 0.76))
    colorbar = fig.colorbar(image, cax=cax)
    colorbar.outline.set_edgecolor("#c8d0d9")
    colorbar.ax.set_ylim(vmin, displayed_max)
    tick_step = 10 if displayed_max <= 40 else 20
    ticks = list(range(vmin, displayed_max + 1, tick_step))
    if ticks[-1] != displayed_max:
        ticks.append(displayed_max)
    colorbar.set_ticks(ticks)
    colorbar.ax.tick_params(labelsize=14, length=3, colors="#475467")


def visualize_tensor(
    tensor: np.ndarray,
    out_path: str,
    title: str,
    show_all_rows: bool,
    marker_column: Optional[int],
) -> int:
    if tensor.ndim != 3 or tensor.shape[0] != 5:
        raise ValueError(f"expected a 5-channel (5, H, W) tensor, got shape {tensor.shape}")

    configure_style()
    view = trim_padding_rows(tensor, show_all_rows)
    _, rows, _ = view.shape
    fig = plt.figure(figsize=(15.0, 10.0), layout="constrained")
    grid = GridSpec(
        5, 2, figure=fig, width_ratios=(1.0, 0.14), hspace=0.10, wspace=0.03
    )
    axes = []
    side_axes = []
    for row in range(5):
        axes.append(fig.add_subplot(grid[row, 0], sharex=axes[0] if axes else None))
        side_axes.append(fig.add_subplot(grid[row, 1]))

    plot_discrete_track(axes[0], side_axes[0], view[0], BASE_LABELS, BASE_COLORS)
    populated_read = view[0] != 0
    populated_read[0, :] = False
    quality_values = np.ma.masked_where(~populated_read, view[1])
    plot_continuous_track(
        fig, axes[1], side_axes[1], quality_values, QUALITY_CMAP, 0,
        QUALITY_COLOR_SCALE_MAX, BASE_QUALITY_DISPLAY_MAX,
    )

    mapping_values = np.ma.masked_where((~populated_read) | (view[3] < 0), view[3])
    plot_continuous_track(
        fig, axes[2], side_axes[2], mapping_values, QUALITY_CMAP, 0,
        QUALITY_COLOR_SCALE_MAX, MAPPING_QUALITY_DISPLAY_MAX,
    )
    axes[2].set_facecolor("#ffffff")
    mismatch_values = np.where(view[2] == 5, 1, view[2])
    plot_discrete_track(
        axes[3], side_axes[3], mismatch_values, MISMATCH_LABELS, MISMATCH_COLORS
    )
    plot_discrete_track(axes[4], side_axes[4], view[4], CIGAR_LABELS, CIGAR_COLORS)

    for channel, ax in enumerate(axes):
        ax.set_title(CHANNEL_TITLES[channel], loc="left", fontsize=18, pad=6)
        ax.set_ylabel("Read row", fontsize=16)
        ax.grid(axis="x", color="#e9edf2", linewidth=0.5)
        ax.set_xlim(-0.5, tensor.shape[2] - 0.5)
        if marker_column is not None:
            ax.axvline(
                marker_column, color="#101828", linewidth=0.8, linestyle=(0, (3, 3)),
            )
        if channel != len(axes) - 1:
            ax.tick_params(labelbottom=False)
    axes[-1].set_xlabel("Window position (bp)", fontsize=17)
    fig.suptitle(title, fontsize=24, fontweight="semibold", color="#101828")
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    fig.savefig(out_path, dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Visualize one five-channel pileup tensor from a .npy tensor or shard."
    )
    parser.add_argument("npy_path", help="Input .npy path: (N,5,H,W) shard or (5,H,W) tensor.")
    parser.add_argument(
        "-i", "--sample-index", type=int, default=0, help="Sample index within a shard (default: 0)."
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output PNG for single-sample mode. Default: output-dir/<basename>_sample_<index>.png.",
    )
    parser.add_argument(
        "--output-dir",
        help="Figure directory. Default: <project>/tensor_figures/<input_basename>/.",
    )
    parser.add_argument(
        "--all-samples",
        action="store_true",
        help="Render every sample in a shard to output-dir. For large shards this can take substantial time and disk.",
    )
    parser.add_argument(
        "--classification",
        choices=("true", "false"),
        help="Render only classified samples in this shard, using variant_summary_classified.ndjson.",
    )
    parser.add_argument(
        "--summary-path",
        help="Classified NDJSON path for --classification. Default: beside the input NPY.",
    )
    parser.add_argument(
        "--max-samples",
        type=int,
        help="Limit rendered figures for --all-samples or --classification.",
    )
    parser.add_argument(
        "--show-all-rows",
        action="store_true",
        help="Show all padded read rows instead of cropping below the final populated row.",
    )
    parser.add_argument(
        "--marker-column",
        type=int,
        help="Column to mark; use -1 to hide. Default: infer from mismatch flag 5, else midpoint.",
    )
    parser.add_argument("--title", help="Optional figure title.")
    args = parser.parse_args()

    path = os.path.expanduser(args.npy_path)
    if not os.path.isfile(path):
        sys.exit(f"Error: input NPY does not exist: {path}")

    try:
        array = load_npy(path)
    except (ValueError, OSError) as exc:
        sys.exit(f"Error: {exc}")

    basename = os.path.splitext(os.path.basename(path))[0]
    default_output_dir = os.path.join(DEFAULT_OUTPUT_ROOT, basename)
    if args.classification:
        default_output_dir = os.path.join(
            default_output_dir, f"{args.classification}_variants"
        )
    output_dir = os.path.expanduser(
        args.output_dir or default_output_dir
    )
    if args.classification and args.all_samples:
        parser.error("--classification cannot be combined with --all-samples")
    if (args.all_samples or args.classification) and args.output:
        parser.error("--output cannot be combined with batch rendering; use --output-dir")
    if args.summary_path and not args.classification:
        parser.error("--summary-path requires --classification")
    if args.max_samples is not None and not (args.all_samples or args.classification):
        parser.error("--max-samples requires --all-samples or --classification")
    if args.max_samples is not None and args.max_samples <= 0:
        parser.error("--max-samples must be greater than zero")

    if args.classification:
        if array.ndim != 4:
            parser.error("--classification requires a shard with shape (N, C, H, W)")
        summary_path = os.path.expanduser(
            args.summary_path
            or os.path.join(os.path.dirname(path), "variant_summary_classified.ndjson")
        )
        if not os.path.isfile(summary_path):
            parser.error(f"classified summary does not exist: {summary_path}")
        try:
            shard_index = infer_shard_index(path)
            sample_indices = classified_sample_indices(
                summary_path, args.classification, shard_index, array.shape[0]
            )
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
        if args.max_samples is not None:
            sample_indices = sample_indices[:args.max_samples]
        if not sample_indices:
            parser.error(
                f"no {args.classification} variants found for shard index {shard_index}"
            )
    elif args.all_samples:
        if array.ndim != 4:
            parser.error("--all-samples requires a shard with shape (N, C, H, W)")
        if not 0 <= args.sample_index < array.shape[0]:
            parser.error(
                f"--sample-index must be within 0..{array.shape[0] - 1}"
            )
        stop_index = array.shape[0]
        if args.max_samples is not None:
            stop_index = min(stop_index, args.sample_index + args.max_samples)
        sample_indices = list(range(args.sample_index, stop_index))
    else:
        sample_indices = [args.sample_index]

    total_to_render = len(sample_indices)
    if args.all_samples or args.classification:
        print(f"Rendering {total_to_render} figures to: {os.path.abspath(output_dir)}")
    rendered = 0
    last_output = ""
    for sample_index in sample_indices:
        try:
            tensor, index = load_single_tensor(array, sample_index)
        except IndexError as exc:
            sys.exit(f"Error: {exc}")
        sample_suffix = f"_sample_{sample_index:05d}" if index is not None else ""
        output = os.path.expanduser(
            args.output or os.path.join(output_dir, f"{basename}{sample_suffix}.png")
        )
        classification_prefix = f"{args.classification} / " if args.classification else ""
        figure_title = args.title or f"{basename} / {classification_prefix}sample {sample_index:05d}"
        marker_column = (
            infer_marker_column(tensor) if args.marker_column is None else args.marker_column
        )
        marker_column = None if marker_column < 0 else marker_column
        if marker_column is not None and marker_column >= tensor.shape[2]:
            sys.exit(
                f"Error: marker column {marker_column} is outside width 0..{tensor.shape[2] - 1}"
            )

        try:
            displayed_rows = visualize_tensor(
                tensor, output, figure_title, args.show_all_rows, marker_column,
            )
        except ValueError as exc:
            sys.exit(f"Error: {exc}")
        rendered += 1
        last_output = os.path.abspath(output)
        if not (args.all_samples or args.classification):
            print(f"Displayed rows: {displayed_rows}/{tensor.shape[1]}")
        elif rendered == 1 or rendered % 100 == 0 or rendered == total_to_render:
            print(f"Rendered {rendered}/{total_to_render}")

    print(f"Input shape: {array.shape}; dtype: {array.dtype}")
    print(
        f"Quality color scale: 0-{QUALITY_COLOR_SCALE_MAX}; "
        f"displayed colorbars: base 0-{BASE_QUALITY_DISPLAY_MAX}, "
        f"mapping 0-{MAPPING_QUALITY_DISPLAY_MAX}"
    )
    print(f"Rendered figures: {rendered}")
    if args.all_samples or args.classification:
        print(f"Saved in: {os.path.abspath(output_dir)}")
    else:
        print(f"Saved: {last_output}")


if __name__ == "__main__":
    main()
