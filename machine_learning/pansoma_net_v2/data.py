"""Merged Pansoma tensor sets -> samples (x int8 (8, 200, 101), blocks int16 (4,), scalars float32 (9,), label).

A tensor set is a directory with SNV/ and INDEL/ in the merged layout of tensor_postprocessing:
<chrom>_shard_NNNNN_data.npy (N, 8, 200, 101) int8, <chrom>_shard_NNNNN_labels.npy (N,) int8
(-1 ignore, 0 non, 1 somatic, 2 germline), <chrom>_variant_summary.ndjson (one record per tensor),
manifest.json and labels.manifest.json.

The index of one kind directory holds per tensor: chromosome, node, shard, row in the shard, label, the ends of
the site's row blocks A1, ALT (A2..Ak), REF, OTHER (rows are ordered in these blocks; the summary's row_groups
give them exactly, and rows from the OTHER end on are padding) and the site's scalars (SCALARS: counts before
the 200-row cap, AFs, alleles, event length). Reading the summaries takes about a minute per 2.5 M tensors, so
the index is cached (cache_dir) and rebuilt when the labels change.

From <chrom>_labels.ndjson it also keeps each tensor's label reason and whether its node is off the GRCh38
path (grch38 null). Evaluation (validation and the test chromosome) uses `eval_label`: the label, except that
off-reference tensors without a truth match (reason off_reference_no_truth_match, -1 for training) count as
non, because test-time calling meets them inside the BED. The other -1 reasons are left out: outside the BED,
below the AF floors and without a GRCh38 position can all be dropped without truth.
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
INDEX_VERSION = 4  # 4: label reasons, off-reference flags, eval labels; 3: node ids, scalars, cache names
EVAL_AS_NON = ("off_reference_no_truth_match",)  # -1 for training, 0 (non) when evaluating
SCALARS = ("log_coverage", "log_site_coverage", "log_alt_count", "log_ref_count", "log_other_count", "af",
           "second_allele_af", "allele_count", "log_event_length")
_ROW_GROUPS = re.compile(rb'"row_groups": (\[[^\]]*\])')
_CANDIDATE = re.compile(rb'"candidate_id": "([^"]+)"')
_SHARD = re.compile(rb'"shard_file": "([^"]+)"')
_ROW = re.compile(rb'"index_within_shard": (\d+)')
_NUMBER = {k: re.compile(rb'"' + k.encode() + rb'": ([-+0-9.eE]+)')
           for k in ("node_id", "coverage", "alt_count", "ref_count", "other_count", "af", "site_coverage",
                     "allele_count", "second_allele_af", "event_length")}  # first match = the record's own (A1) field
_REASON = re.compile(rb'"reason": "([^"]+)"')
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


def site_scalars(v):
    """SCALARS from a record's numbers (dict)."""
    log = math.log1p
    return [log(v["coverage"]), log(v["site_coverage"]), log(v["alt_count"]), log(v["ref_count"]),
            log(v["other_count"]), v["af"], v["second_allele_af"], v["allele_count"], log(v["event_length"])]


def _labels_meta(kind_dir):
    m = json.loads((kind_dir / "labels.manifest.json").read_text())
    return {k: m.get(k) for k in ("format", "version", "rules_sha256", "created", "tensors", "snv_min_af", "indel_min_af")}


def build_index(kind_dir):
    """Index arrays of one merged kind directory (every tensor, labels -1 included)."""
    kind_dir = Path(kind_dir)
    manifest = json.loads((kind_dir / "manifest.json").read_text())
    if manifest.get("layout", "").split("-")[0] != "chromosome" or list(manifest.get("shape", [])) != [8, 200, 101]:
        raise ValueError(f"{kind_dir}: not a merged (8, 200, 101) tensor directory")
    chroms, shards = list(manifest["chromosomes"]), []
    cols = {k: [] for k in ("chrom", "node", "shard", "row", "label", "blocks", "scalars", "reason", "off_reference")}
    candidates, reasons = [], {}
    for ci, chrom in enumerate(chroms):
        info = manifest["chromosomes"][chrom]
        local = {}
        for s in info["shards"]:
            local[s["file"]] = len(shards)
            shards.append(s["file"])
        labels = {f: np.load(kind_dir / f.replace("_data.npy", "_labels.npy")) for f in local}
        before = len(cols["chrom"])
        label_lines = open(kind_dir / f"{chrom}_labels.ndjson", "rb")  # same order as the summary
        with open(kind_dir / info["summary"], "rb") as src, label_lines:
            for line in src:
                lab = label_lines.readline()
                if _CANDIDATE.search(lab).group(1) != _CANDIDATE.search(line).group(1):
                    raise ValueError(f"{kind_dir}/{chrom}_labels.ndjson is not in the summary's order")
                reason = _REASON.search(lab).group(1).decode()
                cols["reason"].append(reasons.setdefault(reason, len(reasons)))
                cols["off_reference"].append(b'"grch38": null' in lab)
                f = _SHARD.search(line).group(1).decode()
                row = int(_ROW.search(line).group(1))
                v = {k: float(p.search(line).group(1)) for k, p in _NUMBER.items()}
                cols["chrom"].append(ci)
                cols["node"].append(int(v["node_id"]))
                cols["shard"].append(local[f])
                cols["row"].append(row)
                cols["label"].append(int(labels[f][row]))
                cols["blocks"].append(block_ends(json.loads(_ROW_GROUPS.search(line).group(1))))
                cols["scalars"].append(site_scalars(v))
                candidates.append(_CANDIDATE.search(line).group(1).decode())
        if len(cols["chrom"]) - before != info["tensors"]:
            raise ValueError(f"{kind_dir}/{info['summary']}: summary does not list {info['tensors']} tensors")
    arrays = dict(chrom=np.asarray(cols["chrom"], np.int16), node=np.asarray(cols["node"], np.int64),
                  shard=np.asarray(cols["shard"], np.int32), row=np.asarray(cols["row"], np.int32),
                  label=np.asarray(cols["label"], np.int8), blocks=np.asarray(cols["blocks"], np.int16).reshape(-1, 4),
                  scalars=np.asarray(cols["scalars"], np.float32).reshape(-1, len(SCALARS)),
                  reason=np.asarray(cols["reason"], np.int16), off_reference=np.asarray(cols["off_reference"], bool))
    vocabulary = sorted(reasons, key=reasons.get)
    eval_label = arrays["label"].copy()
    eval_label[np.isin(arrays["reason"], [reasons[r] for r in EVAL_AS_NON if r in reasons])] = 0
    arrays["eval_label"] = eval_label
    meta = dict(index_version=INDEX_VERSION, directory=str(kind_dir.resolve()), chroms=chroms, shards=shards,
                labels=_labels_meta(kind_dir), tensors=int(manifest["tensors"]), scalars=list(SCALARS),
                reasons=vocabulary, eval_as_non=list(EVAL_AS_NON))
    return arrays, candidates, meta


class KindIndex:
    """The cached index of one merged kind directory (e.g. <sample>/<set>/SNV)."""

    def __init__(self, kind_dir, cache_dir=None, rebuild=False):
        self.dir = Path(kind_dir).resolve()
        cache = None
        if cache_dir is not None:  # e.g. COLO829T_Illumina.v3_tensors.SNV.<hash of the path>
            key = hashlib.sha1(str(self.dir).encode()).hexdigest()[:12]
            cache = Path(cache_dir) / f"{self.dir.parent.parent.name}.{self.dir.parent.name}.{self.dir.name}.{key}"
        npz, text = (Path(f"{cache}.npz"), Path(f"{cache}.candidates.txt")) if cache is not None else (None, None)
        current = _labels_meta(self.dir)
        if npz is not None and not rebuild and npz.exists() and text.exists():
            with np.load(npz) as z:
                meta = json.loads(str(z["meta"]))
                if (meta.get("index_version") == INDEX_VERSION and meta["labels"] == current
                        and meta["directory"] == str(self.dir)):
                    self.arrays = {k: z[k] for k in z.files if k != "meta"}
                    self.meta, self._candidates, self._candidates_path = meta, None, text
                    return
        arrays, candidates, meta = build_index(self.dir)
        self.arrays, self.meta, self._candidates, self._candidates_path = arrays, meta, candidates, None
        if npz is not None:
            npz.parent.mkdir(parents=True, exist_ok=True)
            tmp = Path(f"{cache}.{os.getpid()}.tmp.npz")  # per process: jobs sharing a cache may build at once
            np.savez(tmp, meta=np.array(json.dumps(meta)), **arrays)
            Path(f"{cache}.{os.getpid()}.candidates.tmp").write_text("\n".join(candidates) + "\n")
            Path(f"{cache}.{os.getpid()}.candidates.tmp").replace(text)
            tmp.replace(npz)  # the .npz last: a cache is used only when both files exist
            self._candidates_path = text

    def __len__(self):
        return len(self.arrays["label"])

    @property
    def kind(self):
        return self.dir.name

    def candidates(self):
        if self._candidates is None:
            self._candidates = self._candidates_path.read_text().split("\n")[:len(self)]
        return self._candidates

    def select(self, chroms=None, labelled=True, evaluation=False):
        """Positions of the tensors on `chroms` (None = all): with a training label >= 0 (labelled), with an
        evaluation label >= 0 (evaluation), or all (labelled=False)."""
        keep = np.ones(len(self), bool)
        if chroms is not None:
            wanted = [i for i, c in enumerate(self.meta["chroms"]) if c in set(chroms)]
            keep &= np.isin(self.arrays["chrom"], wanted)
        if evaluation:
            keep &= self.arrays["eval_label"] >= 0
        elif labelled:
            keep &= self.arrays["label"] >= 0
        return np.flatnonzero(keep)

    def reason(self, position):
        return self.meta["reasons"][int(self.arrays["reason"][position])]


def block_split(index, positions, fraction, block_nodes=20000, seed=0):
    """(train, validation) positions: whole blocks of `block_nodes` consecutive node IDs of one chromosome (about
    1 Mb on HPRC v1.1 d9) go to validation with probability `fraction`, so no validation tensor shares reads
    with a training tensor across a block edge. The choice depends only on (chromosome, block, seed)."""
    chrom = index.arrays["chrom"][positions].astype(np.int64)
    block = index.arrays["node"][positions] // block_nodes
    keys, inverse = np.unique(np.stack([chrom, block], 1), axis=0, return_inverse=True)
    names = index.meta["chroms"]
    u = np.array([int.from_bytes(hashlib.sha1(f"{seed}:{names[c]}:{b}".encode()).digest()[:8], "big") / 2 ** 64
                  for c, b in keys.tolist()])
    is_val = (u < fraction)[inverse.reshape(-1)]
    return positions[~is_val], positions[is_val]


class TensorDataset(Dataset):
    """Samples from several KindIndex selections: (x int8 (8, 200, 101), blocks int16 (4,), scalars float32 (9,),
    label int)."""

    def __init__(self, parts, labels="label"):
        """parts: [(KindIndex, positions)]; labels: "label" (training) or "eval_label" (evaluation)."""
        # 53 bytes per sample, pickled to every DataLoader worker
        self.files, shard, row, label, blocks, scalars = [], [], [], [], [], []
        for index, positions in parts:
            offset = len(self.files)
            self.files += [str(index.dir / f) for f in index.meta["shards"]]
            a = index.arrays
            shard.append(a["shard"][positions] + offset)
            row.append(a["row"][positions])
            label.append(a[labels][positions])
            blocks.append(a["blocks"][positions])
            scalars.append(a["scalars"][positions])
        cat = lambda xs, shape, dtype: np.concatenate(xs).astype(dtype) if xs else np.zeros(shape, dtype)  # noqa: E731
        self.shard, self.row = cat(shard, 0, np.int32), cat(row, 0, np.int32)
        self.label, self.blocks = cat(label, 0, np.int8), cat(blocks, (0, 4), np.int16)
        self.scalars = cat(scalars, (0, len(SCALARS)), np.float32)
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
        return (torch.from_numpy(x), torch.from_numpy(self.blocks[i].copy()), torch.from_numpy(self.scalars[i].copy()),
                int(self.label[i]))

    def __getstate__(self):  # workers open their own files
        state = dict(self.__dict__)
        state["_open"] = {}
        return state

    def class_counts(self):
        return np.bincount(self.label[self.label >= 0].astype(np.int64), minlength=len(CLASSES))


class EpochSampler(Sampler):
    """A fresh random (or fixed) order each epoch, optionally only num_samples of it, split over DDP ranks:
    shuffled training pads to equal per-rank lengths; evaluation (shuffle=False) gives disjoint slices.
    With labels and non_fraction < 1 each epoch takes every tensor of the other classes and a fresh random
    non_fraction of the non tensors (label 0)."""

    def __init__(self, n, num_samples=None, shuffle=True, rank=0, world=1, seed=0, labels=None, non_fraction=1.0):
        self.n, self.num_samples, self.shuffle = n, num_samples, shuffle
        self.rank, self.world, self.seed, self.epoch = rank, world, seed, 0
        self.non, self.other, self.keep_non = None, None, None
        if labels is not None and non_fraction < 1.0:
            labels = np.asarray(labels)
            self.non, self.other = torch.from_numpy(np.flatnonzero(labels == 0)), torch.from_numpy(np.flatnonzero(labels != 0))
            self.keep_non = int(round(non_fraction * len(self.non)))

    def set_epoch(self, epoch):
        self.epoch = epoch

    def epoch_size(self):
        return self.n if self.non is None else len(self.other) + self.keep_non

    def _order(self):
        if self.shuffle:
            g = torch.Generator()
            g.manual_seed(self.seed + self.epoch)
            if self.non is None:
                order = torch.randperm(self.n, generator=g)
            else:
                pool = torch.cat([self.other, self.non[torch.randperm(len(self.non), generator=g)[:self.keep_non]]])
                order = pool[torch.randperm(len(pool), generator=g)]
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
        m = self.epoch_size() if self.num_samples is None else min(self.epoch_size(), self.num_samples)
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
