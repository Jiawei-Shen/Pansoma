import numpy as np
import argparse
import os


def inspect_npy_file(file_path):
    """
    Loads an .npy file and prints its properties.
    """
    if not os.path.exists(file_path):
        print(f"Error: File not found at {file_path}")
        return

    if not file_path.lower().endswith(".npy"):
        print(f"Error: File {file_path} is not an .npy file.")
        return

    try:
        data = np.load(file_path)
        print(f"--- Inspecting: {file_path} ---")
        print(f"Shape: {data.shape}")
        print(f"Data type (dtype): {data.dtype}")

        # If the data is numerical, print some basic statistics
        if np.issubdtype(data.dtype, np.number):
            print(f"Min value: {np.min(data)}")
            print(f"Max value: {np.max(data)}")
            print(f"Mean value: {np.mean(data)}")

            # Print a small slice of the data
            if data.ndim == 1:
                print(f"Sample slice (first 10 elements): {data[:10]}")
            elif data.ndim == 2:
                print(f"Sample slice (up to first 5x5 elements):\n{data[:5, :5]}")
            elif data.ndim == 3:
                print(f"Sample slice (up to first 2x2x2 elements from corners):\n{data[:2, :2, :2]}")
                # If you suspect it's HxWxC or CxHxW, you can print more specific slices
                # For HxWxC (e.g., height, width, channels):
                # print(f"Slice assuming HxWxC (first pixel, all channels): {data[0, 0, :]}")
                # print(f"Slice assuming HxWxC (first channel, 2x2 patch): {data[:2, :2, 0]}")
            elif data.ndim > 3:
                print(f"Data has {data.ndim} dimensions. Showing slice of the first element along each dim.")
                slicing = tuple([slice(0, 2) for _ in range(data.ndim)])
                print(f"{data[slicing]}")
        else:
            print("Data is not a numerical type.")

        print("-" * 30)

    except Exception as e:
        print(f"Error loading or inspecting file {file_path}: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Inspect an .npy file.")
    parser.add_argument("file_path", type=str, help="Path to the .npy file to inspect.")

    args = parser.parse_args()
    inspect_npy_file(args.file_path)

    # --- How to get a list of your .npy files ---
    # If you want to inspect one file from your dataset,
    # you might need to find its full path first.
    # For example, if your data_path is '/home/jiawei/Documents/Dockers/PansomaNet/data/COLO829T_SNV_chr1_tensors'
    # and it has subdirectories 'train/false/*.npy'
    #
    # import glob
    # data_root = '/home/jiawei/Documents/Dockers/PansomaNet/data/COLO829T_SNV_chr1_tensors/train/false/'
    # npy_files = glob.glob(os.path.join(data_root, '*.npy'))
    # if npy_files:
    #     print("\nFound .npy files. Example to inspect:")
    #     print(f"python inspect_npy.py \"{npy_files[0]}\"")
    # else:
    #     print(f"\nNo .npy files found in {data_root} to provide an example path.")