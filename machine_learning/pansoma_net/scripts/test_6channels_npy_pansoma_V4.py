#!/usr/bin/env python3
"""
ONE-STOP pipeline (optimized):

Main speedups vs your original:
1) Vectorized keep-mask on GPU (loop only over kept indices, not every sample).
2) Faster variant_summary lookup: store as per-shard list (O(1) index), avoids (sidx,iws) tuple hashing.
3) Avoid per-sample shard_index parsing in the hot path: precompute per-file shard_index once in Dataset.
4) Move normalization to GPU (optional, enabled by default) to reduce DataLoader worker CPU load.
5) Use torch.inference_mode() + optional AMP autocast for faster inference on CUDA.

Behavior/outputs remain the same (no PoN filtering here; only inference->join->VCF outputs).
"""

import argparse
import csv
import json
import os
import re
import sys
import time
from typing import Dict, List, Tuple, Optional, Any, Iterator, Iterable, Set
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm

# Import your model the same way as in training
THIS_DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(THIS_DIR, "..")))
from mynet import ConvNeXtCBAMClassifier  # noqa: E402

# pysam for bgzip + tabix
try:
    import pysam
except Exception:
    pysam = None


# --- EXACT SAME NORMALIZATION AS VAL (COLO829T Testing) ---
VAL_MEAN = torch.tensor([
    18.417816162109375, 12.649129867553711, -0.5452527403831482,
    24.723854064941406, 4.690611362457275, 0.2813551473402196
], dtype=torch.float32)
VAL_STD = torch.tensor([
    25.028322219848633, 14.809632301330566, 0.6181337833404541,
    29.972835540771484, 7.9231791496276855, 0.7659083659074717
], dtype=torch.float32)


# ----------------------------- Utilities ------------------------------------- #

def ensure_pysam_or_die():
    if pysam is None:
        msg = (
            "ERROR: pysam is required to bgzip and index the VCF.\n"
            "Install with:  pip install pysam\n"
            "Or in conda:   conda install -c bioconda pysam"
        )
        print(msg, file=sys.stderr)
        sys.exit(2)


def escape_info(v) -> str:
    return (
        str(v)
        .replace(" ", "_")
        .replace(",", "%2C")
        .replace(";", "%3B")
        .replace("=", "%3D")
    )


def _scan_tree(start_dir: str) -> List[str]:
    """Iteratively scan one subtree with os.scandir (fast)."""
    stack = [start_dir]
    out: List[str] = []
    while stack:
        d = stack.pop()
        try:
            with os.scandir(d) as it:
                for e in it:
                    try:
                        if e.is_dir(follow_symlinks=False):
                            stack.append(e.path)
                        else:
                            if e.name.lower().endswith(".npy"):
                                out.append(os.path.abspath(e.path))
                    except OSError:
                        continue
        except (PermissionError, FileNotFoundError):
            continue
    return out


def _format_eta(seconds: float) -> str:
    seconds = int(max(0, seconds))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    if h > 0:
        return f"{h:d}h{m:02d}m{s:02d}s"
    elif m > 0:
        return f"{m:d}m{s:02d}s"
    else:
        return f"{s:d}s"


def list_npy_files_parallel(root: str, workers: int) -> List[str]:
    """Parallel file discovery with tqdm bar."""
    root = os.path.abspath(os.path.expanduser(root))
    top_dirs, files = [], []

    try:
        with os.scandir(root) as it:
            for e in it:
                try:
                    if e.is_dir(follow_symlinks=False):
                        top_dirs.append(e.path)
                    elif e.name.lower().endswith(".npy"):
                        files.append(os.path.abspath(e.path))
                except OSError:
                    continue
    except FileNotFoundError:
        return []

    if workers <= 1 or not top_dirs:
        for d in top_dirs:
            files.extend(_scan_tree(d))
        return sorted(set(files))

    start = time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex, \
         tqdm(total=len(top_dirs), desc="Scan", unit="dir", leave=False, dynamic_ncols=True) as bar:
        futures = [ex.submit(_scan_tree, d) for d in top_dirs]
        for fut in as_completed(futures):
            files.extend(fut.result())
            bar.update(1)
            done = bar.n
            elapsed = time.time() - start
            speed = done / max(1e-9, elapsed)
            eta = (len(top_dirs) - done) / max(1e-9, speed)
            bar.set_postfix_str(f"speed={speed:.1f} dir/s ETA={_format_eta(eta)}")

    return sorted(set(files))


def read_file_list(txt_path: str) -> List[str]:
    """Read one path per line."""
    with open(txt_path, "rb") as f:
        lines = f.read().splitlines()
    try:
        return [ln.decode("utf-8").strip() for ln in lines if ln.strip()]
    except UnicodeDecodeError:
        return [ln.strip() for ln in lines if ln.strip()]


def save_file_list(paths: List[str], txt_path: str) -> None:
    """Write one absolute path per line."""
    txt_path = os.path.abspath(os.path.expanduser(txt_path))
    os.makedirs(os.path.dirname(txt_path) or ".", exist_ok=True)
    with open(txt_path, "w") as f:
        for p in paths:
            f.write(p + "\n")
    print(f"Saved {len(paths)} paths to {txt_path}")


def print_data_info(files: List[str], mmap: bool, k: int) -> None:
    print("\n=== Data Information ===")
    print(f"Total .npy files: {len(files)} | mmap={mmap}")
    for i, p in enumerate(files[:k]):
        try:
            arr = np.load(p, mmap_mode="r" if mmap else None)
            print(f"  [{i+1}] {p} | shape={arr.shape}, dtype={arr.dtype}")
        except Exception as e:
            print(f"  [{i+1}] {p} | <error: {e}>")
    print("=" * 26)


def detect_input_mode(files: List[str], mmap: bool) -> str:
    """
    Auto-detect by inspecting the first readable file:
      - if ndim==4 -> shard
      - if ndim==3 -> single
    """
    for p in files:
        try:
            arr = np.load(p, mmap_mode="r" if mmap else None)
            if arr.ndim == 4:
                return "shard"
            if arr.ndim == 3:
                return "single"
        except Exception:
            continue
    return "single"


def to_chw_float32(x: np.ndarray, channels: int) -> np.ndarray:
    """
    Convert numpy array to float32 CHW:
      - CHW if already (C,H,W)
      - HWC if (H,W,C)
    """
    if x.ndim != 3:
        raise ValueError(f"Expected 3D tensor, got shape={x.shape}")
    if x.shape[0] == channels:
        chw = x
    elif x.shape[-1] == channels:
        chw = np.transpose(x, (2, 0, 1))
    else:
        raise ValueError(f"Unexpected shape {x.shape}; expected channels={channels}")
    return chw.astype(np.float32, copy=False)


# ---------------------------- Variant summary (latest shards) ---------------- #

def iter_variant_summary_ndjson(path: str) -> Iterator[dict]:
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or not s.startswith("{"):
                continue
            try:
                obj = json.loads(s)
                if isinstance(obj, dict):
                    yield obj
            except Exception:
                continue


def load_variant_summary_by_shard_list(ndjson_path: str) -> Dict[int, List[Optional[dict]]]:
    """
    Faster lookup structure:
      vs_by_shard[sidx][iws] -> meta dict (or None)

    This avoids (sidx,iws) tuple hashing overhead in the hot loop.
    """
    vs_by_shard: Dict[int, List[Optional[dict]]] = {}
    missing = 0
    total = 0

    with tqdm(desc="Load variant_summary", unit="rec", dynamic_ncols=True, leave=True) as bar:
        for row in iter_variant_summary_ndjson(ndjson_path):
            total += 1
            bar.update(1)
            try:
                sidx = int(row.get("shard_index"))
                iws = int(row.get("index_within_shard"))
            except Exception:
                missing += 1
                continue

            lst = vs_by_shard.get(sidx)
            if lst is None:
                lst = []
                vs_by_shard[sidx] = lst
            if iws >= len(lst):
                lst.extend([None] * (iws + 1 - len(lst)))
            lst[iws] = row

    if missing:
        tqdm.write(f"WARNING: {missing} variant_summary rows missing shard_index/index_within_shard; skipped.")
    kept = sum(1 for _sidx, lst in vs_by_shard.items() for x in lst if x is not None)
    tqdm.write(f"Loaded variant_summary: shards={len(vs_by_shard)} populated={kept} / total_rows={total}")
    return vs_by_shard


# ------------------------ TSV (ALT->REF) mapping ------------------------ #

_SHARD_BASENAME_RE = re.compile(r"(?:^|/)(?:shard_)?(?P<idx>\d+)(?:_data)?\.npy$", re.IGNORECASE)

def shard_index_from_path(shard_path: str) -> Optional[int]:
    if not shard_path:
        return None
    bn = os.path.basename(str(shard_path))
    m = _SHARD_BASENAME_RE.search(bn)
    if not m:
        m2 = re.search(r"(\d+)", bn)
        if not m2:
            return None
        try:
            return int(m2.group(1))
        except Exception:
            return None
    try:
        return int(m.group("idx"))
    except Exception:
        return None


def load_alt_to_ref_tsv(tsv_path: Optional[str]) -> Tuple[Dict[int, int], Set[int], Set[int]]:
    if not tsv_path:
        return {}, set(), set()

    alt_to_ref: Dict[int, int] = {}
    tsv_ref_nodes: Set[int] = set()
    tsv_alt_nodes: Set[int] = set()

    with open(tsv_path, "r", encoding="utf-8") as f:
        try:
            csv.field_size_limit(sys.maxsize)
        except OverflowError:
            csv.field_size_limit(2**31 - 1)
        reader = csv.DictReader(f, delimiter="\t")
        hdr = [c.strip() for c in (reader.fieldnames or [])]
        ref_node_key = None
        alt_nodes_key = None
        for k in hdr:
            lk = k.lower()
            if ref_node_key is None and lk in ("ref_node", "refnode", "ref"):
                ref_node_key = k
            if alt_nodes_key is None and lk in ("alt_node(s)", "alt_nodes", "altnodes", "alt"):
                alt_nodes_key = k
        if ref_node_key is None or alt_nodes_key is None:
            raise ValueError("TSV must have REF_NODE and ALT_NODE(S) columns (case-insensitive).")

        for row in reader:
            try:
                ref_node = int(str(row[ref_node_key]).strip())
                tsv_ref_nodes.add(ref_node)
            except Exception:
                continue

            alts_raw = str(row[alt_nodes_key]).strip()
            for token in re.findall(r"\d+", alts_raw):
                try:
                    alt_node = int(token)
                except Exception:
                    continue
                alt_to_ref[alt_node] = ref_node
                tsv_alt_nodes.add(alt_node)

    return alt_to_ref, tsv_ref_nodes, tsv_alt_nodes


# ------------------------ Node map (JSON/JSONL) ------------------------ #

def _first_present(d: dict, keys: Iterable[str], default=None):
    for k in keys:
        if k in d:
            return d[k]
    return default

def _to_int(x) -> Optional[int]:
    try:
        return int(x)
    except Exception:
        return None

def iter_node_map_records(json_path: str, jsonl: bool) -> Iterable[dict]:
    if jsonl:
        with open(json_path, "r", encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if not s or not s.startswith("{"):
                    continue
                try:
                    obj = json.loads(s)
                    if isinstance(obj, dict):
                        yield obj
                except Exception:
                    continue
    else:
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, list):
            for obj in data:
                if isinstance(obj, dict):
                    yield obj
        elif isinstance(data, dict):
            for _, v in data.items():
                if isinstance(v, dict):
                    yield v
        else:
            raise ValueError("Unsupported JSON structure for node map.")

def _has_key_case_insensitive(obj: dict, key_name: str) -> bool:
    kn = key_name.lower()
    for k, v in obj.items():
        if str(k).lower() == kn:
            return True
        if isinstance(v, dict) and _has_key_case_insensitive(v, key_name):
            return True
        if isinstance(v, (list, tuple)):
            for elt in v:
                if isinstance(elt, dict) and _has_key_case_insensitive(elt, key_name):
                    return True
    return False

def load_node_map(
    json_path: str,
    jsonl: bool,
    node_start_is_1_based: bool = True
) -> Tuple[
    Dict[int, Tuple[str, int, str, int]],
    Set[int],
    Set[int],
]:
    anchors: Dict[int, Tuple[str, int, str, int]] = {}
    nodes_seen: Set[int] = set()
    nodes_with_genomead_af: Set[int] = set()

    for rec in iter_node_map_records(json_path, jsonl):
        nid = _first_present(rec, ["node_id", "node", "id", "nid"])
        node_id = _to_int(nid)
        if node_id is None:
            continue
        nodes_seen.add(node_id)

        if _has_key_case_insensitive(rec, "genomead_af"):
            nodes_with_genomead_af.add(node_id)

        chrom = _first_present(rec, ["chrom", "chr", "chromosome", "contig"])
        start1 = _first_present(rec, ["grch38_position_start", "start_1based", "pos1"])
        start0 = _first_present(rec, ["start_0based", "pos0"])
        strand = str(_first_present(rec, ["strand_in_path", "strand"], "+") or "+")
        length = _to_int(_first_present(rec, ["length", "len", "node_length"]))

        if start0 is not None:
            s0 = _to_int(start0)
        elif start1 is not None:
            s1 = _to_int(start1)
            s0 = (s1 - 1) if s1 is not None else None
        else:
            generic = _first_present(rec, ["start", "pos", "start_pos"])
            if generic is not None:
                sg = _to_int(generic)
                if sg is not None:
                    s0 = (sg - 1) if node_start_is_1_based else sg
                else:
                    s0 = None
            else:
                s0 = None

        if chrom is not None and s0 is not None and length is not None:
            anchors[node_id] = (str(chrom), int(s0), strand if strand in ("+", "-") else "+", int(length))

    return anchors, nodes_seen, nodes_with_genomead_af


# ---------------------------- Dataset ---------------------------------------- #

class NpyShardsDataset(Dataset):
    """
    Optimized:
      - precompute per-file shard_index once
      - return (tensor, index_in_shard, shard_index_int)
      - shard mode always mmap
    """
    def __init__(
        self,
        files: List[str],
        channels: int,
        mmap: bool = False,
        input_mode: str = "shard",   # "single" or "shard"
        cache_size: int = 2,
    ):
        self.files = files
        self.channels = channels
        self.mmap = mmap
        self.input_mode = input_mode
        self.cache_size = max(1, int(cache_size))

        if self.channels != 6:
            raise ValueError(f"Expected 6 channels, got {channels}")

        # Precompute shard_index per file once (no regex in hot path)
        self._file_shard_index: List[Optional[int]] = [shard_index_from_path(p) for p in self.files]

        # Per-process cache (each DataLoader worker has its own dataset copy)
        self._arr_cache: Dict[str, np.ndarray] = {}
        self._arr_cache_order: List[str] = []

        # Build an index: global_sample_idx -> (file_idx, local_idx)
        self._offsets: List[int] = [0]
        self._counts: List[int] = []

        for fp in tqdm(self.files, desc="Index", unit="file", dynamic_ncols=True, leave=False):
            try:
                arr = np.load(fp, mmap_mode="r")
            except Exception:
                self._counts.append(0)
                self._offsets.append(self._offsets[-1])
                continue

            if self.input_mode == "single":
                self._counts.append(1 if arr.ndim == 3 else 0)
            else:
                self._counts.append(int(arr.shape[0]) if arr.ndim == 4 else 0)

            self._offsets.append(self._offsets[-1] + self._counts[-1])

        self.total_samples = self._offsets[-1]
        if self.total_samples == 0:
            raise RuntimeError("No valid samples found in provided .npy files.")

    def __len__(self):
        return self.total_samples

    def _locate(self, global_idx: int) -> Tuple[int, int]:
        lo, hi = 0, len(self._offsets) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if self._offsets[mid + 1] <= global_idx:
                lo = mid + 1
            else:
                hi = mid
        file_idx = lo
        local_idx = global_idx - self._offsets[file_idx]
        return file_idx, local_idx

    def _get_array(self, path: str) -> np.ndarray:
        if path in self._arr_cache:
            return self._arr_cache[path]

        if self.input_mode == "shard":
            arr = np.load(path, mmap_mode="r")
        else:
            arr = np.load(path, mmap_mode="r" if self.mmap else None)

        self._arr_cache[path] = arr
        self._arr_cache_order.append(path)

        if len(self._arr_cache_order) > self.cache_size:
            evict = self._arr_cache_order.pop(0)
            self._arr_cache.pop(evict, None)

        return arr

    def __getitem__(self, idx: int):
        file_idx, local_idx = self._locate(idx)
        shard_path = self.files[file_idx]
        shard_index = self._file_shard_index[file_idx]

        arr = self._get_array(shard_path)

        if self.input_mode == "single":
            x = arr
            index_in_shard = 0
        else:
            x = arr[local_idx]
            index_in_shard = int(local_idx)

        chw = to_chw_float32(x, self.channels)
        t = torch.from_numpy(chw)  # float32 tensor (CPU)

        # return fewer objects (no shard_path)
        return t, index_in_shard, (-1 if shard_index is None else int(shard_index))


# --------------------------- Model Loading ---------------------------------- #

def build_and_load_model(ckpt_path: str, device: torch.device,
                         depths: List[int], dims: List[int]):
    if not os.path.isfile(ckpt_path):
        raise FileNotFoundError(f"Checkpoint not found: {ckpt_path}")
    checkpoint = torch.load(ckpt_path, map_location=device)

    genotype_map = checkpoint.get("genotype_map", {})
    in_channels = int(checkpoint.get("in_channels", 6))

    model = ConvNeXtCBAMClassifier(
        in_channels=in_channels,
        class_num=len(genotype_map) or 2,
        depths=depths,
        dims=dims
    ).to(device)

    state = checkpoint.get("model_state_dict", checkpoint)
    model.load_state_dict(state, strict=True)
    model.eval()
    return model, genotype_map, in_channels


# --------------------------- Prediction helpers ------------------------------ #

def decide_keep_and_filter(pred_is_true: bool,
                           prob_true: Optional[float],
                           min_true_prob: Optional[float],
                           refcall_prob: Optional[float]) -> Tuple[bool, str]:
    if pred_is_true:
        if (min_true_prob is None) or (prob_true is None) or (prob_true >= min_true_prob):
            return True, "PASS"
        return False, "PASS"
    else:
        if (refcall_prob is not None) and (prob_true is not None) and (prob_true >= refcall_prob):
            return True, "RefCall"
        return False, "PASS"


# --------------------------- VCF builders/writers ---------------------------- #

def write_vcf_bgzip(
    records: List[Tuple[str, int, str, str, str, str, str, str]],
    out_vcfgz: str,
    meta_id: str,
    extra_header_lines: Optional[List[str]] = None,
    do_index: bool = True,
    sort_records: bool = True,
    chrom_sort_mode: str = "node",  # "node" or "linear"
):
    ensure_pysam_or_die()

    if not out_vcfgz.endswith(".vcf.gz"):
        if out_vcfgz.endswith(".vcf"):
            out_vcfgz = out_vcfgz + ".gz"
        else:
            out_vcfgz = out_vcfgz + ".vcf.gz"

    out_vcfgz = os.path.abspath(out_vcfgz)
    os.makedirs(os.path.dirname(out_vcfgz) or ".", exist_ok=True)
    tmp_vcf = out_vcfgz[:-3]  # strip ".gz" -> ".vcf"

    if sort_records or do_index:
        if chrom_sort_mode == "linear":
            def chrom_key_linear(c: str):
                m = re.fullmatch(r"(?:chr)?(\d+)", c, flags=re.IGNORECASE)
                if m:
                    return (0, int(m.group(1)))
                cl = c.lower()
                if cl in ("chrx", "x"): return (1, 23)
                if cl in ("chry", "y"): return (1, 24)
                if cl in ("chrm", "chrmt", "m", "mt"): return (1, 25)
                return (2, cl)
            records.sort(key=lambda r: (chrom_key_linear(r[0]), int(r[1])))
        else:
            def chrom_key_node(c: str):
                try:
                    return (0, int(c))
                except Exception:
                    return (1, c)
            records.sort(key=lambda r: (chrom_key_node(r[0]), int(r[1])))

    contigs: List[str] = []
    seen: Set[str] = set()
    for chrom, *_ in records:
        if chrom not in seen:
            seen.add(chrom)
            contigs.append(chrom)

    with open(tmp_vcf, "w", encoding="utf-8") as out:
        out.write("##fileformat=VCFv4.2\n")
        if extra_header_lines:
            for h in extra_header_lines:
                if h and h.startswith("##") and not h.lower().startswith("##fileformat"):
                    out.write(h.rstrip("\n") + "\n")

        out.write(f'##META=<ID={meta_id},Description="Pipeline: inference -> shard-join (variant_summary) -> VCF output.">\n')
        out.write("##INFO=<ID=PROB,Number=1,Type=Float,Description=\"Model probability for positive/'true' class\">\n")
        out.write("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Alt allele frequency from variant_summary\">\n")
        out.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Coverage at locus from variant_summary\">\n")
        out.write("##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths: ref,alt\">\n")
        out.write("##INFO=<ID=OTHER,Number=1,Type=Integer,Description=\"Other allele count at locus\">\n")
        out.write("##INFO=<ID=BQ,Number=1,Type=Float,Description=\"Mean alt allele base quality\">\n")
        out.write("##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant type from variant_summary (X/I/D)\">\n")
        out.write("##INFO=<ID=SHARD,Number=1,Type=Integer,Description=\"Shard index\">\n")
        out.write("##INFO=<ID=IDX,Number=1,Type=Integer,Description=\"Index within shard\">\n")
        out.write("##INFO=<ID=NID,Number=1,Type=String,Description=\"Original node_id\">\n")
        out.write("##INFO=<ID=NSO,Number=1,Type=String,Description=\"Original node offset (same convention as v_pos)\">\n")
        out.write("##FILTER=<ID=RefCall,Description=\"Non-true prediction retained because prob >= --refcall_prob\">\n")

        for c in contigs:
            out.write(f"##contig=<ID={c}>\n")

        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        with tqdm(total=len(records), desc=f"Write {meta_id}", unit="rec", dynamic_ncols=True, leave=True) as bar:
            for chrom, pos, vid, ref, alt, qual, filt, info in records:
                out.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\n")
                bar.update(1)

    if do_index:
        produced_gz = pysam.tabix_index(tmp_vcf, preset="vcf", force=True, keep_original=False)
        produced_tbi = produced_gz + ".tbi"
        if os.path.abspath(produced_gz) != os.path.abspath(out_vcfgz):
            try:
                if os.path.exists(out_vcfgz):
                    os.remove(out_vcfgz)
            except Exception:
                pass
            os.replace(produced_gz, out_vcfgz)
            if os.path.exists(produced_tbi):
                os.replace(produced_tbi, out_vcfgz + ".tbi")
        return out_vcfgz, out_vcfgz + ".tbi"
    else:
        try:
            if os.path.exists(out_vcfgz):
                os.remove(out_vcfgz)
        except Exception:
            pass
        try:
            pysam.tabix_compress(tmp_vcf, out_vcfgz, force=True)
        except TypeError:
            pysam.tabix_compress(tmp_vcf, out_vcfgz)
        try:
            os.remove(tmp_vcf)
        except Exception:
            pass
        return out_vcfgz, None


def build_node_vcf_record_from_meta(
    meta: dict,
    prob_true: Optional[float],
    filt: str,
    pos_is_1based: bool,
) -> Optional[Tuple[str, int, str, str, str, str, str, str]]:
    try:
        node_id = str(meta["node_id"])
    except Exception:
        return None

    try:
        vpos = int(meta.get("v_pos"))
    except Exception:
        return None

    pos_vcf = vpos if pos_is_1based else (vpos + 1)

    ref = str(meta.get("v_ref", "N"))
    alt = str(meta.get("v_alt", "N"))
    vtype = meta.get("v_type")

    try:
        sidx = int(meta.get("shard_index"))
    except Exception:
        sidx = None
    try:
        iws = int(meta.get("index_within_shard"))
    except Exception:
        iws = None

    info: Dict[str, Optional[str]] = {
        "PROB": f"{prob_true:.6g}" if prob_true is not None else None,
        "TYPE": str(vtype) if vtype is not None else None,
        "SHARD": str(sidx) if sidx is not None else None,
        "IDX": str(iws) if iws is not None else None,
        "NID": str(node_id),
        "NSO": str(vpos),
    }

    af = meta.get("alt_allele_frequency")
    cov = meta.get("coverage_at_locus")
    rc  = meta.get("ref_allele_count_at_locus")
    ac  = meta.get("alt_allele_count")
    oc  = meta.get("other_allele_count_at_locus")
    bq  = meta.get("mean_alt_allele_base_quality")

    if af is not None:
        try:
            info["AF"] = f"{float(af):.6g}"
        except Exception:
            info["AF"] = escape_info(af)
    if cov is not None:
        if isinstance(cov, (int, float)):
            try:
                info["DP"] = str(int(cov))
            except Exception:
                info["DP"] = escape_info(cov)
        else:
            info["DP"] = escape_info(cov)
    if rc is not None and ac is not None:
        try:
            info["AD"] = f"{int(rc)},{int(ac)}"
        except Exception:
            info["AD"] = f"{rc},{ac}"
    if oc is not None:
        try:
            info["OTHER"] = str(int(oc))
        except Exception:
            info["OTHER"] = escape_info(oc)
    if bq is not None:
        try:
            info["BQ"] = f"{float(bq):.3f}"
        except Exception:
            info["BQ"] = escape_info(bq)

    info_str = ";".join(f"{k}={escape_info(v)}" for k, v in info.items() if v is not None) or "."
    return (str(node_id), int(pos_vcf), ".", ref, alt, ".", filt, info_str)


def convert_node_record_to_linear(
    node_rec: Tuple[str, int, str, str, str, str, str, str],
    anchors: Dict[int, Tuple[str, int, str, int]],
    alt_to_ref: Dict[int, int],
    offset_is_0_based: bool,
) -> Optional[Tuple[str, int, str, str, str, str, str, str]]:
    chrom_node = node_rec[0]
    info = node_rec[7]

    try:
        node_id = int(chrom_node)
    except Exception:
        return None

    nso = None
    for field in info.split(";"):
        if field.startswith("NSO="):
            nso = field.split("=", 1)[1]
            break
    if nso is None:
        return None
    try:
        offset = int(nso)
    except Exception:
        return None

    base_node = node_id
    ref_node_candidate = alt_to_ref.get(base_node)
    anchor = anchors.get(ref_node_candidate) if ref_node_candidate is not None else None
    if anchor is None:
        anchor = anchors.get(base_node)
    if anchor is None:
        return None

    chrom_lin, start0, strand, length = anchor
    off0 = offset if offset_is_0_based else (offset - 1)
    if strand == "+":
        pos_lin = start0 + off0 + 1
    else:
        pos_lin = start0 + (length - 1 - off0) + 1

    return (chrom_lin, int(pos_lin), node_rec[2], node_rec[3], node_rec[4], node_rec[5], node_rec[6], node_rec[7])


# --------------------------- Inference Core (optimized) ---------------------- #

def run_inference_and_build_vcfs(
    model,
    dl,
    device,
    true_idx: Optional[int],
    total_samples: int,
    vs_by_shard: Dict[int, List[Optional[dict]]],
    min_true_prob: Optional[float],
    refcall_prob: Optional[float],
    pos_is_1based: bool,
    gpu_normalize: bool,
    use_amp: bool,
) -> Tuple[List[Tuple[str,int,str,str,str,str,str,str]], int, int, int]:
    """
    Returns:
      node_vcf_records_all, total_inferred, kept, missing_meta
    """
    node_records: List[Tuple[str,int,str,str,str,str,str,str]] = []

    total_inferred = 0
    kept = 0
    missing_meta = 0

    mean_dev = std_dev = None
    if gpu_normalize and device.type == "cuda":
        mean_dev = VAL_MEAN.to(device).view(1, -1, 1, 1)
        std_dev = VAL_STD.to(device).view(1, -1, 1, 1)

    t0 = time.time()
    with tqdm(total=total_samples, desc="Infer", unit="sample", dynamic_ncols=True, leave=True) as bar:
        with torch.inference_mode():
            for images, index_in_shard, shard_index in dl:
                bs = images.shape[0]
                total_inferred += bs

                images = images.to(device, non_blocking=True)
                if mean_dev is not None and std_dev is not None:
                    images = (images - mean_dev) / std_dev

                # forward + probs (AMP optional)
                if device.type == "cuda" and use_amp:
                    with torch.cuda.amp.autocast(dtype=torch.float16):
                        outputs = model(images)
                        if isinstance(outputs, tuple):
                            outputs = outputs[0]
                        probs = torch.softmax(outputs, dim=1)
                else:
                    outputs = model(images)
                    if isinstance(outputs, tuple):
                        outputs = outputs[0]
                    probs = torch.softmax(outputs, dim=1)

                # prob_true (prefer explicit "true" class if known)
                if true_idx is not None and true_idx >= 0 and true_idx < probs.shape[1]:
                    prob_true_vec = probs[:, true_idx]
                    pred_is_true_vec = (probs.argmax(dim=1) == true_idx)
                else:
                    # fallback: treat top prob as "prob_true" and all as "true-like"
                    prob_true_vec, _ = probs.max(dim=1)
                    pred_is_true_vec = torch.ones((bs,), device=probs.device, dtype=torch.bool)

                # vectorized keep mask
                keep_true = pred_is_true_vec
                if (min_true_prob is not None) and (true_idx is not None):
                    keep_true = keep_true & (prob_true_vec >= min_true_prob)

                keep_refcall = torch.zeros_like(keep_true)
                if (refcall_prob is not None) and (true_idx is not None):
                    keep_refcall = (~pred_is_true_vec) & (prob_true_vec >= refcall_prob)

                keep_mask = keep_true | keep_refcall
                keep_idx = torch.nonzero(keep_mask, as_tuple=False).squeeze(1)

                # move just what we need to CPU
                keep_idx_cpu = keep_idx.detach().cpu().numpy()
                prob_true_cpu = prob_true_vec.detach().cpu().numpy()
                pred_is_true_cpu = pred_is_true_vec.detach().cpu().numpy()

                iws_cpu = index_in_shard.detach().cpu().numpy() if torch.is_tensor(index_in_shard) else np.asarray(index_in_shard)
                sidx_cpu = shard_index.detach().cpu().numpy() if torch.is_tensor(shard_index) else np.asarray(shard_index)

                for j in keep_idx_cpu:
                    sidx = int(sidx_cpu[j])
                    if sidx < 0:
                        missing_meta += 1
                        continue
                    iws = int(iws_cpu[j])
                    prob_true = float(prob_true_cpu[j])

                    pred_is_true = bool(pred_is_true_cpu[j])
                    keep, filt = decide_keep_and_filter(pred_is_true, prob_true, min_true_prob, refcall_prob)
                    if not keep:
                        continue

                    lst = vs_by_shard.get(sidx)
                    if (lst is None) or (iws < 0) or (iws >= len(lst)):
                        missing_meta += 1
                        continue
                    meta = lst[iws]
                    if meta is None:
                        missing_meta += 1
                        continue

                    node_rec = build_node_vcf_record_from_meta(
                        meta=meta,
                        prob_true=prob_true,
                        filt=filt,
                        pos_is_1based=pos_is_1based,
                    )
                    if node_rec is None:
                        missing_meta += 1
                        continue

                    node_records.append(node_rec)
                    kept += 1

                bar.update(bs)
                elapsed = time.time() - t0
                speed = total_inferred / max(1e-9, elapsed)
                eta = (total_samples - bar.n) / max(1e-9, speed)
                bar.set_postfix_str(
                    f"kept={kept} miss_meta={missing_meta} speed={speed:.1f} samp/s ETA={_format_eta(eta)}"
                )

    return node_records, total_inferred, kept, missing_meta


# ----------------------------- Main ----------------------------------------- #

def main():
    p = argparse.ArgumentParser(description="Infer over .npy shards and directly emit VCFs (latest shard format).")

    # Input discovery
    p.add_argument("--input_dir", required=False, help="Directory to scan for .npy files (recursive).")
    p.add_argument("--file_list", required=False, help="Text file with one .npy path per line (skips scanning).")
    p.add_argument("--save_file_list", required=False, help="Write discovered .npy paths to this text file.")

    # Model
    p.add_argument("--ckpt", required=True)
    p.add_argument("--batch_size", type=int, default=64)
    p.add_argument("--num_workers", type=int, default=8, help="DataLoader workers.")
    p.add_argument("--scan_workers", type=int, default=None, help="Directory scanning workers (default: num_workers).")
    p.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    p.add_argument("--compile", action="store_true")

    # Data
    p.add_argument("--mmap", action="store_true", help="Use mmap for single mode; shard mode always uses mmap.")
    p.add_argument("--cache_size", type=int, default=2, help="Per-worker LRU cache size for opened shard memmaps.")
    p.add_argument("--input_mode", choices=["auto", "single", "shard"], default="auto",
                   help="auto: detect. single: each .npy is one sample. shard: each .npy contains (N,6,H,W).")
    p.add_argument("--sample_info_k", type=int, default=3)
    p.add_argument("--depths", type=int, nargs="+", default=[3, 3, 27, 3])
    p.add_argument("--dims", type=int, nargs="+", default=[192, 384, 768, 1536])

    # Latest shard meta join
    p.add_argument("--variant_summary", required=True,
                   help="variant_summary.ndjson (must contain shard_index and index_within_shard).")

    # Keep rules
    p.add_argument("--min_true_prob", type=float, default=None,
                   help="Keep only pred_class==true with prob>=this (if prob present).")
    p.add_argument("--refcall_prob", type=float, default=None,
                   help="Also keep pred_class!=true when prob>=this, FILTER=RefCall.")

    # v_pos convention
    p.add_argument("--pos_is_1based", action="store_true",
                   help="If set, treat meta v_pos as already 1-based; otherwise VCF POS=v_pos+1 (default).")

    # VCF outputs
    p.add_argument("--out_prefix", required=True,
                   help="Output prefix. Files: <prefix>.linear.vcf.gz, <prefix>.all.vcf.gz, <prefix>.ref.vcf.gz")
    p.add_argument("--emit", nargs="+", choices=["linear", "all", "ref"], default=["linear"],
                   help="Which VCF(s) to write. Default: linear.")
    p.add_argument("--sort", action="store_true", help="Sort VCF outputs (also enabled automatically for indexing).")
    p.add_argument("--no-index", action="store_true", help="Do not create Tabix indexes (.tbi).")

    # Linear conversion inputs
    p.add_argument("--map_json", default=None, help="Node map JSON (array or dict) for anchors and reference/ref-alt classification.")
    p.add_argument("--map_jsonl", action="store_true", help="Interpret --map_json as JSONL.")
    p.add_argument("--tsv", default=None, help="Optional TSV ALT_NODE(S)->REF_NODE mapping (ALT can borrow REF anchor).")
    p.add_argument("--offset_is_0_based", type=lambda s: str(s).lower() in ("1", "true", "t", "yes", "y"), default=True,
                   help="Offsets used for linear conversion (NSO) are 0-based (default: true). Set false if 1-based.")
    p.add_argument("--node_start_is_1_based", type=lambda s: str(s).lower() in ("1", "true", "t", "yes", "y"), default=True,
                   help="Node starts in --map_json are 1-based (default: true).")

    # New perf toggles
    p.add_argument("--gpu-normalize", action="store_true",
                   help="Normalize on GPU (recommended). If not set, inputs must already be normalized or you accept raw.")
    p.add_argument("--amp", action="store_true",
                   help="Use CUDA autocast (fp16) during inference (recommended on H100/A100).")

    args = p.parse_args()

    want_index = (not args.no_index)
    want_sort = args.sort or want_index

    needs_map = ("linear" in args.emit) or ("ref" in args.emit)
    if needs_map and not args.map_json:
        raise SystemExit("ERROR: --emit includes linear/ref, but --map_json was not provided.")

    # Device
    if args.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)

    if device.type == "cuda":
        torch.backends.cudnn.benchmark = True
        torch.backends.cuda.matmul.allow_tf32 = True

    # Build model
    print(f"Loading Model from {args.ckpt}")
    print(f"depths={args.depths}, dims={args.dims}")
    model, genotype_map, in_channels = build_and_load_model(args.ckpt, device, args.depths, args.dims)
    if in_channels != 6:
        raise ValueError(f"Expected 6-channel input, got {in_channels}")

    if args.compile:
        try:
            model = torch.compile(model)
        except Exception as e:
            print(f"torch.compile not enabled: {e}")

    # Determine true_idx
    true_idx = None
    if len(genotype_map) > 0:
        for cname, cidx in genotype_map.items():
            if str(cname).strip().lower() == "true":
                true_idx = int(cidx)
                break

    # Resolve file list vs directory scan
    if args.file_list:
        files = read_file_list(args.file_list)
        print(f"Loaded {len(files)} paths from --file_list")
    else:
        if not args.input_dir:
            raise SystemExit("Either --file_list or --input_dir must be provided.")
        scan_workers = args.scan_workers if args.scan_workers is not None else args.num_workers
        files = list_npy_files_parallel(args.input_dir, workers=scan_workers)
        if not files:
            raise SystemExit(f"No .npy files found in {args.input_dir}")
        if args.save_file_list:
            save_file_list(files, args.save_file_list)

    print_data_info(files, args.mmap, args.sample_info_k)

    # Decide input mode
    if args.input_mode == "auto":
        mode = detect_input_mode(files, mmap=args.mmap)
        print(f"Auto-detected input_mode={mode}")
    else:
        mode = args.input_mode
        print(f"Using input_mode={mode}")

    # Load variant_summary (fast structure)
    if not os.path.isfile(args.variant_summary):
        raise SystemExit(f"variant_summary not found: {args.variant_summary}")
    vs_by_shard = load_variant_summary_by_shard_list(args.variant_summary)

    # Load map + TSV (if needed)
    anchors: Dict[int, Tuple[str, int, str, int]] = {}
    nodes_seen: Set[int] = set()
    ref_nodes_with_af: Set[int] = set()
    alt_to_ref: Dict[int, int] = {}

    if needs_map:
        anchors, nodes_seen, ref_nodes_with_af = load_node_map(
            args.map_json,
            jsonl=args.map_jsonl,
            node_start_is_1_based=args.node_start_is_1_based
        )
        print(f"Node map loaded: anchors={len(anchors)}, nodes_seen={len(nodes_seen)}, ref_nodes={len(ref_nodes_with_af)}")

        if args.tsv:
            alt_to_ref, tsv_ref_nodes, tsv_alt_nodes = load_alt_to_ref_tsv(args.tsv)
            print(f"TSV loaded: ALT->REF mappings={len(alt_to_ref)} (TSV REF nodes={len(tsv_ref_nodes)}, TSV ALT nodes={len(tsv_alt_nodes)})")

    # Dataset & Loader
    print("Building Data Loader...")
    ds = NpyShardsDataset(
        files,
        channels=in_channels,
        mmap=args.mmap,
        input_mode=mode,
        cache_size=args.cache_size,
    )

    pin_memory = (device.type == "cuda")
    persistent_workers = (args.num_workers > 0)
    prefetch_factor = 4 if args.num_workers > 0 else None

    dl = DataLoader(
        ds,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        pin_memory=pin_memory,
        shuffle=False,
        persistent_workers=persistent_workers,
        prefetch_factor=prefetch_factor,
    )

    # Infer + build node records
    t0 = time.time()
    node_records_all, total_inferred, kept, miss_meta = run_inference_and_build_vcfs(
        model=model,
        dl=dl,
        device=device,
        true_idx=true_idx,
        total_samples=len(ds),
        vs_by_shard=vs_by_shard,
        min_true_prob=args.min_true_prob,
        refcall_prob=args.refcall_prob,
        pos_is_1based=args.pos_is_1based,
        gpu_normalize=args.gpu_normalize,
        use_amp=args.amp,
    )

    elapsed = time.time() - t0
    print(f"\nFinished inference: total={total_inferred} kept={kept} missing_meta={miss_meta} time={elapsed:.2f}s")

    if kept == 0:
        print("No qualifying predictions to write. Exiting.", file=sys.stderr)
        sys.exit(0)

    out_prefix = os.path.abspath(args.out_prefix)
    do_index = (not args.no_index)

    # all (node-form)
    if "all" in args.emit:
        out_all = out_prefix + ".all.vcf.gz"
        gz, tbi = write_vcf_bgzip(
            records=node_records_all,
            out_vcfgz=out_all,
            meta_id="ALL_NODE",
            extra_header_lines=None,
            do_index=do_index,
            sort_records=want_sort,
            chrom_sort_mode="node",
        )
        print(f"[all]   {gz}" + (f"  (index: {tbi})" if tbi else ""))

    # ref (node-form subset)
    if "ref" in args.emit:
        if not ref_nodes_with_af:
            print("WARNING: --emit ref requested but ref_nodes_with_af is empty (check --map_json).", file=sys.stderr)
        ref_set = set(ref_nodes_with_af)
        node_records_ref: List[Tuple[str,int,str,str,str,str,str,str]] = []
        for r in node_records_all:
            try:
                nid = int(r[0])
            except Exception:
                continue
            if nid in ref_set:
                node_records_ref.append(r)

        out_ref = out_prefix + ".ref.vcf.gz"
        gz, tbi = write_vcf_bgzip(
            records=node_records_ref,
            out_vcfgz=out_ref,
            meta_id="REF_NODE_ONLY",
            extra_header_lines=None,
            do_index=do_index,
            sort_records=want_sort,
            chrom_sort_mode="node",
        )
        print(f"[ref]   {gz}" + (f"  (index: {tbi})" if tbi else "") + f"  (records={len(node_records_ref)})")

    # linear (converted subset)
    if "linear" in args.emit:
        converted: List[Tuple[str,int,str,str,str,str,str,str]] = []
        unconverted_n = 0

        for nr in node_records_all:
            lin = convert_node_record_to_linear(
                node_rec=nr,
                anchors=anchors,
                alt_to_ref=alt_to_ref,
                offset_is_0_based=args.offset_is_0_based,
            )
            if lin is None:
                unconverted_n += 1
                continue
            converted.append(lin)

        if not converted:
            print("WARNING: No records could be converted to linear coordinates (anchors missing?).", file=sys.stderr)

        out_lin = out_prefix + ".linear.vcf.gz"
        gz, tbi = write_vcf_bgzip(
            records=converted,
            out_vcfgz=out_lin,
            meta_id="LINEAR_CONVERTED",
            extra_header_lines=None,
            do_index=do_index,
            sort_records=want_sort,
            chrom_sort_mode="linear",
        )
        print(f"[linear] {gz}" + (f"  (index: {tbi})" if tbi else "") +
              f"  (records={len(converted)}, dropped_unconverted={unconverted_n})")

    print("\nDone.")


if __name__ == "__main__":
    main()
