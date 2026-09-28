"""Merged Pansoma tensor sets -> samples (x int8 (8, 200, 101), blocks int16 (4,), label).

A tensor set is a directory with SNV/ and INDEL/ in the merged layout of tensor_postprocessing:
<chrom>_shard_NNNNN_data.npy (N, 8, 200, 101) int8, <chrom>_shard_NNNNN_labels.npy (N,) int8
(-1 ignore, 0 non, 1 somatic, 2 germline), <chrom>_variant_summary.ndjson (one record per tensor),
manifest.json and labels.manifest.json.

The index of one kind directory holds per tensor: chromosome, shard, row in the shard, label, A1 AF and the
ends of the site's row blocks A1, ALT (A2..Ak), REF, OTHER (rows are ordered in these blocks; the summary's
row_groups give them exactly, and rows from the OTHER end on are padding). Reading the summaries takes about
a minute per 2.5 M tensors, so the index is cached (cache_dir) and rebuilt when the labels change.
"""
import hashlib
import json
import math
import os
import re
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import Dataset, Sampler

CLASSES = ("non", "somatic", "germline")  # label values 0, 1, 2; -1 = ignore
KINDS = ("SNV", "INDEL")
TENSOR_SHAPE = (8, 200, 101)
MAX_OPEN_SHARDS = 256  # open shard files per DataLoader worker
INDEX_VERSION = 1
_ROW_GROUPS = re.compile(rb'"row_groups": (\[[^\]]*\])')
_CANDIDATE = re.compile(rb'"candidate_id": "([^"]+)"')
_SHARD = re.compile(rb'"shard_file": "([^"]+)"')
_ROW = re.compile(rb'"index_within_shard": (\d+)')
_AF = re.compile(rb'"af": ([-+0-9.eE]+)')
_BLOCK = {"A1": 0, "REF": 2, "OTHER": 3}  # A2, A3, ... -> 1 (ALT)


def block_ends(groups):
    """Row ends (A1, ALT, REF, OTHER) of a summary's row_groups; blocks must be contiguous from row 0."""
    ends, last = [None] * 4, 0
    for g in groups:
        if g["start_row"] != last:
            raise ValueError(f"row_groups not contiguous: {groups}")
        b = _BLOCK.get(g["allele"], 1)
        if any(e is not None for e in ends[b + 1:]):
            raise ValueError(f"row_groups out of order: {groups}")
        ends[b], last = g["end_row"], g["end_row"]
    out, prev = [], 0
    for e in ends:
        prev = prev if e is None else e
        out.append(prev)
    return out


def _labels_meta(kind_dir):
    m = json.loads((kind_dir / "labels.manifest.json").read_text())
    return {k: m.get(k) for k in ("version", "created", "tensors", "snv_min_af")}


def build_index(kind_dir):
    """Index arrays of one merged kind directory (every tensor, labels -1 included)."""
    kind_dir = Path(kind_dir)
    manifest = json.loads((kind_dir / "manifest.json").read_text())
    if manifest.get("layout", "").split("-")[0] != "chromosome" or list(manifest.get("shape", [])) != [8, 200, 101]:
        raise ValueError(f"{kind_dir}: not a merged (8, 200, 101) tensor directory")
    chroms, shards = list(manifest["chromosomes"]), []
    cols = {k: [] for k in ("chrom", "shard", "row", "label", "af", "blocks")}
    candidates = []
    for ci, chrom in enumerate(chroms):
        info = manifest["chromosomes"][chrom]
        local = {}
        for s in info["shards"]:
            local[s["file"]] = len(shards)
            shards.append(s["file"])
        labels = {f: np.load(kind_dir / f.replace("_data.npy", "_labels.npy")) for f in local}
        before = len(cols["chrom"])
        with open(kind_dir / info["summary"], "rb") as src:
            for line in src:
                f = _SHARD.search(line).group(1).decode()
                row = int(_ROW.search(line).group(1))
                cols["chrom"].append(ci)
                cols["shard"].append(local[f])
                cols["row"].append(row)
                cols["label"].append(int(labels[f][row]))
                cols["af"].append(float(_AF.search(line).group(1)))
                cols["blocks"].append(block_ends(json.loads(_ROW_GROUPS.search(line).group(1))))
                candidates.append(_CANDIDATE.search(line).group(1).decode())
        if len(cols["chrom"]) - before != info["tensors"]:
            raise ValueError(f"{kind_dir}/{info['summary']}: summary does not list {info['tensors']} tensors")
    arrays = dict(chrom=np.asarray(cols["chrom"], np.int16), shard=np.asarray(cols["shard"], np.int32),
                  row=np.asarray(cols["row"], np.int32), label=np.asarray(cols["label"], np.int8),
                  af=np.asarray(cols["af"], np.float32), blocks=np.asarray(cols["blocks"], np.int16).reshape(-1, 4))
    meta = dict(index_version=INDEX_VERSION, directory=str(kind_dir.resolve()), chroms=chroms, shards=shards,
                labels=_labels_meta(kind_dir), tensors=int(manifest["tensors"]))
    return arrays, candidates, meta


class KindIndex:
    """The cached index of one merged kind directory (e.g. <set>/SNV)."""

    def __init__(self, kind_dir, cache_dir=None, rebuild=False):
        self.dir = Path(kind_dir).resolve()
        cache = None
        if cache_dir is not None:
            key = hashlib.sha1(str(self.dir).encode()).hexdigest()[:12]
            cache = Path(cache_dir) / f"{self.dir.parent.name}.{self.dir.name}.{key}"
        current = _labels_meta(self.dir)
        if cache is not None and not rebuild and (cache.with_suffix(".npz")).exists():
            with np.load(cache.with_suffix(".npz")) as z:
                meta = json.loads(str(z["meta"]))
                if meta.get("index_version") == INDEX_VERSION and meta["labels"] == current:
                    self.arrays = {k: z[k] for k in z.files if k != "meta"}
                    self.meta, self._candidates_path = meta, cache.with_suffix(".candidates.txt")
                    return
        arrays, candidates, meta = build_index(self.dir)
        self.arrays, self.meta, self._candidates, self._candidates_path = arrays, meta, candidates, None
        if cache is not None:
            cache.parent.mkdir(parents=True, exist_ok=True)
            tmp = cache.with_suffix(".tmp.npz")
            np.savez(tmp, meta=np.array(json.dumps(meta)), **arrays)
            tmp.replace(cache.with_suffix(".npz"))
            text = cache.with_suffix(".candidates.txt")
            text.with_suffix(".tmp").write_text("\n".join(candidates) + "\n")
            text.with_suffix(".tmp").replace(text)
            self._candidates_path = text

    def __len__(self):
        return len(self.arrays["label"])

    @property
    def kind(self):
        return self.dir.name

    def candidates(self):
        if getattr(self, "_candidates", None) is None:
            self._candidates = self._candidates_path.read_text().split("\n")[:len(self)]
        return self._candidates

    def select(self, chroms=None, labelled=True):
        """Positions of the tensors on `chroms` (None = all) with a label >= 0 (labelled) or any label."""
        keep = np.ones(len(self), bool)
        if chroms is not None:
            wanted = [i for i, c in enumerate(self.meta["chroms"]) if c in set(chroms)]
            keep &= np.isin(self.arrays["chrom"], wanted)
        if labelled:
            keep &= self.arrays["label"] >= 0
        return np.flatnonzero(keep)


class TensorDataset(Dataset):
    """Samples from several KindIndex selections; returns (x int8 (8, 200, 101), blocks int16 (4,), label int)."""

    def __init__(self, parts):
        """parts: [(KindIndex, positions)]."""
        # 17 bytes per sample, pickled to every DataLoader worker
        self.files, shard, row, label, blocks = [], [], [], [], []
        for index, positions in parts:
            offset = len(self.files)
            self.files += [str(index.dir / f) for f in index.meta["shards"]]
            a = index.arrays
            shard.append(a["shard"][positions] + offset)
            row.append(a["row"][positions])
            label.append(a["label"][positions])
            blocks.append(a["blocks"][positions])
        cat = lambda xs, shape, dtype: np.concatenate(xs).astype(dtype) if xs else np.zeros(shape, dtype)  # noqa: E731
        self.shard, self.row = cat(shard, 0, np.int32), cat(row, 0, np.int32)
        self.label, self.blocks = cat(label, 0, np.int8), cat(blocks, (0, 4), np.int16)
        self._open = {}

    def __len__(self):
        return len(self.label)

    def _file(self, f):
        """(fd, offset of row 0) of shard f. Tensors are read with pread, not a memory map: mapped pages stay
        in every worker's RSS, and Slurm's summed-RSS limit counts them once per worker."""
        if f not in self._open:
            with open(self.files[f], "rb") as fh:
                version = np.lib.format.read_magic(fh)
                read = np.lib.format.read_array_header_1_0 if version == (1, 0) else np.lib.format.read_array_header_2_0
                shape, fortran, dtype = read(fh)
                offset = fh.tell()
            if tuple(shape[1:]) != TENSOR_SHAPE or dtype != np.int8 or fortran:
                raise ValueError(f"{self.files[f]}: {shape} {dtype} is not a C-order int8 (N, 8, 200, 101) array")
            if len(self._open) >= MAX_OPEN_SHARDS:  # stay far below the per-process file limit
                os.close(self._open.pop(next(iter(self._open)))[0])
            self._open[f] = (os.open(self.files[f], os.O_RDONLY), offset)
        return self._open[f]

    def __getitem__(self, i):
        fd, offset = self._file(self.shard[i])
        x = np.empty(TENSOR_SHAPE, np.int8)
        if os.preadv(fd, [x], offset + int(self.row[i]) * x.nbytes) != x.nbytes:
            raise IOError(f"{self.files[self.shard[i]]}: short read at row {self.row[i]}")
        return torch.from_numpy(x), torch.from_numpy(self.blocks[i].copy()), int(self.label[i])

    def __getstate__(self):  # workers open their own files
        state = dict(self.__dict__)
        state["_open"] = {}
        return state

    def class_counts(self):
        return np.bincount(self.label[self.label >= 0].astype(np.int64), minlength=len(CLASSES))


class EpochSampler(Sampler):
    """A fresh random (or fixed) order each epoch, optionally only num_samples of it, split over DDP ranks:
    shuffled training pads to equal per-rank lengths; evaluation (shuffle=False) gives disjoint slices."""

    def __init__(self, n, num_samples=None, shuffle=True, rank=0, world=1, seed=0):
        self.n, self.num_samples, self.shuffle = n, num_samples, shuffle
        self.rank, self.world, self.seed, self.epoch = rank, world, seed, 0

    def set_epoch(self, epoch):
        self.epoch = epoch

    def _order(self):
        if self.shuffle:
            g = torch.Generator()
            g.manual_seed(self.seed + self.epoch)
            order = torch.randperm(self.n, generator=g)
        else:
            order = torch.arange(self.n)
        if self.num_samples is not None:
            order = order[:self.num_samples]
        if self.shuffle and self.world > 1 and len(order):
            per = math.ceil(len(order) / self.world)
            order = torch.cat([order, order[:per * self.world - len(order)]])
        return order[self.rank::self.world]

    def __iter__(self):
        return iter(self._order().tolist())

    def __len__(self):
        m = self.n if self.num_samples is None else min(self.n, self.num_samples)
        if self.shuffle and self.world > 1:
            return math.ceil(m / self.world)
        return max(0, math.ceil((m - self.rank) / self.world))


def load_parts(tensor_sets, kinds, cache_dir, chroms=None, labelled=True, rebuild=False):
    """[(KindIndex, positions)] for every set x kind present."""
    parts = []
    for s in tensor_sets:
        for kind in kinds:
            d = Path(s) / kind
            if (d / "manifest.json").exists():
                index = KindIndex(d, cache_dir, rebuild=rebuild)
                parts.append((index, index.select(chroms, labelled)))
    if not parts:
        raise FileNotFoundError(f"no merged {'/'.join(kinds)} directory under {list(map(str, tensor_sets))}")
    return parts
