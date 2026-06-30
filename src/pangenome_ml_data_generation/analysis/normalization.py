import numpy as np
import os
import argparse
from tqdm import tqdm


def calculate_channel_stats(directory):
    """
    Recursively finds all .npy files in a directory and calculates the mean and
    standard deviation for each channel across all files.

    This function uses Welford's algorithm for a numerically stable, one-pass
    calculation of the mean and variance, which is memory-efficient.

    Args:
        directory (str): The path to the directory to scan for .npy files.

    Returns:
        tuple: A tuple containing two numpy arrays (mean, std).
    """
    print(f"Scanning for .npy files in '{directory}'...")

    # Find all .npy files recursively
    npy_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".npy"):
                npy_files.append(os.path.join(root, file))

    if not npy_files:
        print("Error: No .npy files found in the specified directory.")
        return None, None

    print(f"Found {len(npy_files)} files. Starting calculation...")

    # Initialize variables for Welford's algorithm
    count = 0
    mean = None
    M2 = None  # Sum of squares of differences from the current mean

    # Use tqdm for a progress bar
    for file_path in tqdm(npy_files, desc="Processing files"):
        try:
            # Load the numpy array from the file
            data = np.load(file_path)

            # Initialize mean and M2 on the first file
            if mean is None:
                # Assuming data shape is (channels, height, width) or just (channels,)
                # We calculate stats along the channel axis (axis 0)
                num_channels = data.shape[0]
                mean = np.zeros(num_channels)
                M2 = np.zeros(num_channels)

            # Flatten the spatial dimensions to get a (channels, num_pixels) array
            # This handles both 1D and multi-dimensional channel data
            reshaped_data = data.reshape(data.shape[0], -1)

            # Iterate over each pixel/data point in the file
            for i in range(reshaped_data.shape[1]):
                pixel = reshaped_data[:, i]
                count += 1
                delta = pixel - mean
                mean += delta / count
                delta2 = pixel - mean
                M2 += delta * delta2

        except Exception as e:
            print(f"\nWarning: Could not process file {file_path}. Error: {e}")
            continue

    if count == 0:
        print("Error: Could not process any files successfully.")
        return None, None

    # Calculate variance and standard deviation
    variance = M2 / (count - 1)  # Use (count - 1) for sample variance
    std_dev = np.sqrt(variance)

    return mean, std_dev


def main():
    """Main function to parse arguments and run the calculation."""
    parser = argparse.ArgumentParser(
        description="Calculate mean and standard deviation from .npy files for dataset normalization.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "directory",
        type=str,
        help="The root directory containing the .npy files."
    )
    args = parser.parse_args()

    if not os.path.isdir(args.directory):
        print(f"Error: Directory not found at '{args.directory}'")
        return

    mean, std = calculate_channel_stats(args.directory)

    if mean is not None and std is not None:
        print("\n--- Calculation Complete ---")
        print(f"Calculated Mean:\n{list(mean)}")
        print(f"\nCalculated Std Dev:\n{list(std)}")
        print("\n--- PyTorch transforms.Normalize format ---")
        print(f"transforms.Normalize(\n    mean={list(mean)},\n    std={list(std)}\n)")


if __name__ == "__main__":
    main()
