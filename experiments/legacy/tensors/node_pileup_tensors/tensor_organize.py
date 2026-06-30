#!/usr/bin/env python3

import argparse
import os
import random
import shutil
from multiprocessing import Pool, cpu_count
from functools import partial


def prepare_split_jobs(src_dir, dest_dir, ratio, label, use_symlink):
    files = sorted(os.listdir(src_dir))
    random.shuffle(files)
    total = len(files)
    n_train = int(total * ratio[0])
    n_val = int(total * ratio[1])
    n_test = total - n_train - n_val

    splits = {
        "train": files[:n_train],
        "val": files[n_train:n_train + n_val],
        "test": files[n_train + n_val:]
    }

    jobs = []
    for split, split_files in splits.items():
        split_path = os.path.join(dest_dir, split, label)
        os.makedirs(split_path, exist_ok=True)
        for fname in split_files:
            src_path = os.path.join(src_dir, fname)
            dest_path = os.path.join(split_path, fname)
            jobs.append((src_path, dest_path, use_symlink, label, split))
    return jobs


def copy_or_link_file(job):
    src, dest, use_symlink, label, split = job
    try:
        if use_symlink:
            os.symlink(os.path.abspath(src), dest)
        else:
            shutil.copyfile(src, dest)
    except FileExistsError:
        pass
    return (label, split)


def copy_files_parallel(jobs, workers):
    copied = {}
    with Pool(processes=workers) as pool:
        for i, result in enumerate(pool.imap_unordered(copy_or_link_file, jobs), 1):
            label, split = result
            key = f"{label}/{split}"
            copied[key] = copied.get(key, 0) + 1
            if copied[key] % 10000 == 0:
                print(f"[{key}] Processed {copied[key]} files...")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, help="Source folder with true/ and false/")
    parser.add_argument("--output", required=True, help="Output folder")
    parser.add_argument("--organize", nargs=3, type=float, metavar=('TRAIN', 'VAL', 'TEST'),
                        help="Split ratio, e.g., --organize 0.6 0.2 0.2")
    parser.add_argument("--symlink", action="store_true", help="Use symbolic links instead of copying")
    parser.add_argument("-j", "--workers", type=int, default=cpu_count(), help="Number of parallel workers")

    args = parser.parse_args()

    assert abs(sum(args.organize) - 1.0) < 1e-6, "Split ratio must sum to 1.0"
    all_jobs = []

    for label in ['true', 'false']:
        src_subdir = os.path.join(args.source, label)
        assert os.path.isdir(src_subdir), f"Missing directory: {src_subdir}"
        jobs = prepare_split_jobs(src_subdir, args.output, args.organize, label, args.symlink)
        all_jobs.extend(jobs)

    copy_files_parallel(all_jobs, args.workers)

    print(f"\n✅ Done. Data organized into {args.output}/train|val|test with true/false splits.")

if __name__ == "__main__":
    main()
