#!/usr/bin/env python3
import argparse
import os
import torch
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

BASE_TO_INDEX = {
    'A': 2, 'C': 3, 'G': 5, 'T': 7,   # Standard bases
    'N': 1,                           # Unknown or ambiguous base
    '*': 9,                           # Deletion character
    'PADDING': 0                   # Padding
}
INDEX_TO_BASE_VIS = {v: k for k, v in BASE_TO_INDEX.items()}
MAX_BASE_INDEX = max(BASE_TO_INDEX.values())


def load_tensor(filepath):
    ext = os.path.splitext(filepath)[1].lower()
    if ext == '.pth':
        tensor = torch.load(filepath, map_location='cpu')
        if isinstance(tensor, torch.Tensor):
            tensor = tensor.numpy()
    elif ext == '.npy':
        tensor = np.load(filepath)
    else:
        raise ValueError(f"Unsupported file type: {filepath}")

    if tensor.ndim != 3:
        raise ValueError(f"Expected 3D tensor, got shape {tensor.shape}")

    # Convert (C, W, H) → (C, H, W)
    # tensor = tensor[:, :, :].transpose(2, 0, 1)
    return tensor


def plot_on_ax(ax, channel_data, channel_index, tensor_basename, num_rows, window_size,
               user_cmap_name=None):
    cmap = user_cmap_name or 'gray'
    norm = None

    if channel_index == 0:
        # Custom color map: padding (0) = white, rest follow base indices
        base_colors = ['white', 'plum', 'lightcoral', 'lightskyblue', 'gray',
                       'lightgreen', 'orange', 'gold', 'tan', 'lightgray']
        cmap = mcolors.ListedColormap(base_colors[:MAX_BASE_INDEX + 1])
        norm = mcolors.BoundaryNorm(np.arange(-0.5, MAX_BASE_INDEX + 1.5, 1), cmap.N)

        legend_elements = []
        unique_indices = np.unique(channel_data)
        for i in unique_indices:
            label = INDEX_TO_BASE_VIS.get(i, f"Idx {i}")
            legend_elements.append(
                Line2D([0], [0], marker='s', color='w', label=f'{label} ({i})',
                       markerfacecolor=cmap(norm(i)), markersize=6)
            )
        ax.imshow(channel_data, cmap=cmap, norm=norm, interpolation='nearest', aspect='auto')
        if legend_elements:
            ax.legend(handles=legend_elements, title="Base Legend", fontsize='small',
                      title_fontsize='small', bbox_to_anchor=(1.02, 1), loc='upper left')
    else:
        im = ax.imshow(channel_data, cmap=cmap, interpolation='nearest', aspect='auto')
        ax.figure.colorbar(im, ax=ax, orientation='vertical')

    ax.set_title(f"Channel {channel_index} from {tensor_basename}", fontsize=10)
    ax.set_xlabel("Position", fontsize=9)
    ax.set_ylabel(f"Row (Total: {num_rows})", fontsize=9)
    ax.tick_params(axis='both', labelsize=8)

    if num_rows > 10:
        step = max(1, num_rows // (5 if num_rows <= 40 else 10))
        ax.set_yticks(np.arange(0, num_rows, step))

    # Draw downward arrow to center of the width
    h, w = channel_data.shape  # rows, cols
    center_x = w / 2.0 - 0.5   # align with imshow pixel center
    ax.annotate('▼', xy=(center_x, -1), xytext=(center_x, -3),
                ha='center', va='center',
                fontsize=14, color='black',
                arrowprops=dict(arrowstyle='-|>', color='black'))


def process_file(filepath, output_dir, args):
    try:
        tensor = load_tensor(filepath)
    except Exception as e:
        print(f"Skipping {filepath}: {e}")
        return

    basename = os.path.splitext(os.path.basename(filepath))[0]
    out_path = os.path.join(output_dir, args.output_image or f"{basename}.png")

    num_channels = tensor.shape[0]

    if args.channel is None:
        rows, cols = tensor[0].shape
        fig, axes = plt.subplots(num_channels, 1,
                                 figsize=(max(10, cols / 8), num_channels * max(3, rows / 10) + 2),
                                 sharex=True)
        if num_channels == 1:
            axes = [axes]
        for i in range(num_channels):
            data = tensor[i].astype(float)
            plot_on_ax(axes[i], data, i, basename, data.shape[0], data.shape[1],
                       user_cmap_name=args.cmap)
        fig.suptitle(args.title or f"{basename}", fontsize=16)
        fig.tight_layout(rect=[0, 0.03, 0.90, 0.95])
    else:
        if args.channel >= num_channels:
            print(f"Skipping {filepath}: channel {args.channel} out of bounds (max {num_channels - 1}).")
            return
        data = tensor[args.channel].astype(float)
        rows, cols = data.shape
        fig, ax = plt.subplots(figsize=(max(10, cols / 8), max(6, rows / 10) + 1))
        plot_on_ax(ax, data, args.channel, basename, rows, cols,
                   user_cmap_name=args.cmap)
        fig.tight_layout(rect=[0, 0, 0.90, 1])

    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Visualize multi-channel .npy or .pth tensors.")
    parser.add_argument("input_dir", help="Directory containing .npy or .pth files.")
    parser.add_argument("output_dir", help="Directory to save output images.")
    parser.add_argument("-ch", "--channel", type=int, default=None, help="Channel index to plot (default: all)")
    parser.add_argument("-o", "--output_image", help="Custom output image filename.")
    parser.add_argument("-t", "--title", help="Figure title.")
    parser.add_argument("-c", "--cmap", help="Custom matplotlib colormap.")
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        print(f"Error: input directory not found: {args.input_dir}")
        return
    os.makedirs(args.output_dir, exist_ok=True)

    for fname in os.listdir(args.input_dir):
        if fname.endswith((".npy", ".pth")):
            filepath = os.path.join(args.input_dir, fname)
            process_file(filepath, args.output_dir, args)


if __name__ == "__main__":
    main()
