#!/usr/bin/env python3
import argparse, json, os, shutil, sys, time, random
from multiprocessing import Pool, cpu_count
from typing import Dict, Iterable, List, Set, Tuple
import pysam
from concurrent.futures import ThreadPoolExecutor, as_completed

# ----------------------- worker globals (set via initializer) -----------------
G_VCF = None              # pysam.VariantFile handle
G_CHR = None              # chromosome string
G_NODE_POS: Dict[str,int] = {}
G_TRUE_DIR = None         # may be None when --organize is used (no direct output)
G_FALSE_DIR = None        # may be None when --organize is used (no direct output)
G_USE_SYMLINKS = False

def _init_worker(vcf_path: str, chrom: str, node_pos: Dict[str,int],
                 true_dir: str, false_dir: str, use_symlinks: bool):
    """Called once per worker; avoids pickling huge objects per task."""
    global G_VCF, G_CHR, G_NODE_POS, G_TRUE_DIR, G_FALSE_DIR, G_USE_SYMLINKS
    G_VCF = pysam.VariantFile(vcf_path)
    G_CHR = chrom
    G_NODE_POS = node_pos
    G_TRUE_DIR = true_dir
    G_FALSE_DIR = false_dir
    G_USE_SYMLINKS = use_symlinks
# ------------------------------------------------------------------------------

def format_time(seconds: float) -> str:
    seconds = int(seconds)
    h, r = divmod(seconds, 3600); m, s = divmod(r, 60)
    return f"{h:02}:{m:02}:{s:02}"

def list_node_dirs(base: str) -> Iterable[str]:
    for e in os.scandir(base):
        if e.is_dir():
            yield e.path

def load_needed_node_positions(json_path: str, needed_ids: Set[str]) -> Dict[str,int]:
    """Only keep node_ids that exist in the tensor folder."""
    # Try streaming if ijson is available
    try:
        import ijson  # type: ignore
        pos: Dict[str,int] = {}
        with open(json_path, "rb") as f:
            for node in ijson.items(f, "nodes.item"):
                nid = str(node.get("node_id"))
                if nid in needed_ids:
                    p = node.get("grch38_position_start")
                    if isinstance(p, int):
                        pos[nid] = p
        return pos
    except Exception:
        # Fallback: load once, then filter (uses more RAM)
        with open(json_path) as f:
            data = json.load(f)
        pos = {}
        for node in data.get("nodes", []):
            nid = str(node.get("node_id"))
            if nid in needed_ids:
                p = node.get("grch38_position_start")
                if isinstance(p, int):
                    pos[nid] = p
        return pos

# ---------- NEW: relaxed matching helpers (per-item VCF fetch) ----------------

def _parse_key(variant_key: str) -> Tuple[int, str, str, str]:
    """
    Expected canonical forms used earlier in your pipeline:
      SNP/other:   "{offset}_{type}_{ref}_{alt}"
      Deletion D:  "{anchor}_{D}_{anchor+deleted}_{anchor}"
      Insertion I: "{anchor}_{I}_{anchorBase}_{anchorBase+inserted}"

    Returns: (offset, vtype, REF, ALT) all uppercased.
    """
    parts = variant_key.split("_")
    if len(parts) < 4:
        raise ValueError(f"Unexpected variant_key: {variant_key}")
    offset = int(parts[0])
    vtype  = parts[1]
    ref    = parts[2].upper()
    alt    = parts[3].upper()
    return offset, vtype, ref, alt

def _vcf_has_partial_match(pos: int, v_ref: str, v_alt: str, v_type: str) -> bool:
    """
    Fetch records at POS and test:
      - Deletion (D):     VCF.REF contains v_ref  AND  any VCF.ALT contains v_alt
      - Insertion (I):    any VCF.ALT contains v_alt  AND  VCF.REF contains v_ref
      - Other (e.g. X):   exact REF==v_ref and ALT==v_alt
    """
    # pysam fetch uses 0-based start, half-open end; rec.pos is 1-based
    for rec in G_VCF.fetch(G_CHR, max(0, pos-1), pos+1):
        ref_truth = (rec.ref or "").upper()
        alts_truth = [(a or "").upper() for a in (rec.alts or [])]
        if v_type == "D":
            if ref_truth and (v_ref[len(v_alt):] in ref_truth):
                return True
        elif v_type == "I":
            if ref_truth and (v_ref in ref_truth):
                for a in alts_truth:
                    if v_alt in a:
                        return True
        else:
            for a in alts_truth:
                if ref_truth == v_ref and a == v_alt:
                    return True
    return False

# ------------------------------------------------------------------------------

def _copy_or_link(src: str, dst: str, use_symlinks: bool) -> None:
    """Create dst as a symlink to src, or copy if symlinks not requested/possible."""
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if use_symlinks:
        if os.path.lexists(dst):
            return
        try:
            os.symlink(os.path.abspath(src), dst)
            return
        except OSError:
            pass  # fall back to copy if symlink fails
    if not os.path.exists(dst):
        shutil.copy2(src, dst)

def _classify_node(node_dir: str):
    """
    Worker task: classify one node directory.
    Returns (records, t_cnt, f_cnt).
    Each record includes absolute tensor_path and classification; if G_TRUE_DIR/G_FALSE_DIR
    are set (no --organize), we also copy/link immediately into those dirs.
    """
    node_id = os.path.basename(node_dir)
    start_pos = G_NODE_POS.get(node_id)
    if start_pos is None:
        return [], 0, 0

    summary_path = os.path.join(node_dir, "variant_summary.json")
    try:
        with open(summary_path) as f:
            summary = json.load(f)
    except Exception:
        return [], 0, 0

    variants = summary.get("variants_passing_af_filter", [])
    if not variants:
        return [], 0, 0

    recs = []
    t = f = 0

    for v in variants:
        tf = v.get("tensor_file"); vk = v.get("variant_key")
        if not tf or not vk or not tf.endswith(".npy"):
            continue
        tpath = os.path.join(node_dir, tf)
        if not os.path.isfile(tpath):
            continue
        try:
            off, vtype, ref, alt = _parse_key(vk)
        except Exception:
            continue

        pos = start_pos + off  # absolute GRCh38 position (1-based to match pysam VariantRecord.pos)

        # relaxed matching via per-item fetch
        is_match = _vcf_has_partial_match(pos, ref, alt, vtype)

        label = "true" if is_match else "false"
        if is_match: t += 1
        else:        f += 1

        # If NOT organizing, place now.
        if G_TRUE_DIR and G_FALSE_DIR:
            dest_dir = G_TRUE_DIR if is_match else G_FALSE_DIR
            dst = os.path.join(dest_dir, f"{node_id}_{tf}")
            try:
                _copy_or_link(tpath, dst, G_USE_SYMLINKS)
            except OSError:
                pass

        recs.append({
            "node_id": node_id,
            "tensor_file": tf,
            "tensor_path": os.path.abspath(tpath),
            "genomic_position": pos,
            "ref": ref,
            "alt": alt,
            "type": vtype,
            "classification": label
        })

    return recs, t, f

def _parallel_place(tasks: List[Tuple[str, str, bool]], max_workers: int) -> None:
    """
    Run _copy_or_link over (src, dst, use_symlinks) tasks in parallel (threads).
    """
    total = len(tasks)
    done = 0
    start = time.monotonic()
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = [ex.submit(_copy_or_link, src, dst, use_symlinks) for (src, dst, use_symlinks) in tasks]
        for fut in as_completed(futures):
            _ = fut.result()
            done += 1
            if done % 500 == 0 or done == total:
                elapsed = time.monotonic() - start
                rate = done / elapsed if elapsed > 0 else 0.0
                print(f"\rOrganizing: {done}/{total} ({done/total:.1%})  {rate:.1f} files/s", end="")
    print()

def organize_classified_data_from_summary(summary_ndjson: str, base_dir: str,
                                          ratios, use_symlinks: bool,
                                          seed=None, max_workers: int = 32):
    """
    Read records from NDJSON and populate base_dir/{train,val,test}/{true,false}
    by SYMLINKING (or copying) from tensor_path. No true/false roots are created.
    Uses parallel I/O for placement.
    """
    if seed is not None:
        random.seed(seed)
    assert len(ratios) == 3 and abs(sum(ratios) - 1.0) < 1e-6

    grouped = {"true": [], "false": []}

    with open(summary_ndjson, "r") as f:
        for line in f:
            if not line.strip(): continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            label = rec.get("classification")
            tpath = rec.get("tensor_path")
            tf = rec.get("tensor_file")
            nid = rec.get("node_id")
            if label in grouped and tpath and tf and nid:
                grouped[label].append((tpath, f"{nid}_{tf}"))

    tasks: List[Tuple[str, str, bool]] = []
    for label, items in grouped.items():
        random.shuffle(items)
        n = len(items)
        n_tr = int(n * ratios[0])
        n_va = int(n * ratios[1])
        splits = {
            "train": items[:n_tr],
            "val":   items[n_tr:n_tr + n_va],
            "test":  items[n_tr + n_va:],
        }
        for sp, lst in splits.items():
            outd = os.path.join(base_dir, sp, label)
            os.makedirs(outd, exist_ok=True)
            for src_path, out_name in lst:
                dst_path = os.path.join(outd, out_name)
                tasks.append((src_path, dst_path, use_symlinks))

    if tasks:
        _parallel_place(tasks, max_workers=max_workers)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tensor_folder_path")
    ap.add_argument("vcf_file")
    ap.add_argument("node_pos_json")
    ap.add_argument("--output_folder", default="./classification_results")
    ap.add_argument("--chr", default="chr1")
    ap.add_argument("--use-symlinks", action="store_true",
                    help="Use symlinks instead of copying for output placement.")
    ap.add_argument("-j","--workers", type=int, default=max(1, min(8, cpu_count())))
    ap.add_argument("--chunksize", type=int, default=32)
    ap.add_argument("--maxtasksperchild", type=int, default=500)
    ap.add_argument("--organize", nargs=3, type=float, metavar=("TRAIN","VAL","TEST"),
                    help="If set, skip output_folder/true,false and directly fill train/val/test splits.")
    ap.add_argument("--organize-workers", type=int, default=32,
                    help="Parallel workers for the organize (I/O) step; threads are used.")
    ap.add_argument("--seed", type=int)
    args = ap.parse_args()

    os.makedirs(args.output_folder, exist_ok=True)

    organize_mode = args.organize is not None
    if not organize_mode:
        true_dir  = os.path.join(args.output_folder, "true")
        false_dir = os.path.join(args.output_folder, "false")
        os.makedirs(true_dir, exist_ok=True)
        os.makedirs(false_dir, exist_ok=True)
    else:
        true_dir = None
        false_dir = None

    node_dirs = list(list_node_dirs(args.tensor_folder_path))
    need_ids = {os.path.basename(p) for p in node_dirs}
    print(f"Found {len(node_dirs):,} node dirs.")

    node_pos = load_needed_node_positions(args.node_pos_json, need_ids)
    print(f"Loaded positions for {len(node_pos):,} nodes from JSON.")

    summary_path = os.path.join(args.output_folder, "classification_summary.ndjson")
    summary_f = open(summary_path, "w")

    total_nodes = len(node_dirs); total_true = total_false = 0
    start = time.monotonic()

    with Pool(processes=args.workers,
              initializer=_init_worker,
              initargs=(args.vcf_file, args.chr, node_pos,
                        true_dir, false_dir, args.use_symlinks),
              maxtasksperchild=args.maxtasksperchild) as pool:

        for i, (records, t_cnt, f_cnt) in enumerate(
            pool.imap_unordered(_classify_node, node_dirs, chunksize=args.chunksize)
        ):
            if records:
                for rec in records:
                    summary_f.write(json.dumps(rec) + "\n")

            total_true  += t_cnt
            total_false += f_cnt

            done = i + 1
            elapsed = time.monotonic() - start
            rate = done / elapsed if elapsed > 0 else 0.0
            eta = format_time((total_nodes - done)/rate) if rate > 0 else "..."
            print(f"\rProgress {done}/{total_nodes} "
                  f"({done/total_nodes:.1%})  true={total_true} false={total_false}  "
                  f"{rate:.2f} nodes/s  ETA {eta}    ", end="")

    summary_f.close()
    print("\nDone.")
    print(f"Summary written to: {summary_path}")
    print(f"True: {total_true}  False: {total_false}")

    if organize_mode:
        print("Organizing directly into train/val/test (parallel, no true/false roots)…")
        organize_classified_data_from_summary(
            summary_ndjson=summary_path,
            base_dir=args.output_folder,
            ratios=tuple(map(float, args.organize)),
            use_symlinks=args.use_symlinks,
            seed=args.seed,
            max_workers=max(1, args.organize_workers),
        )
        print("Organization complete.")

if __name__ == "__main__":
    main()
