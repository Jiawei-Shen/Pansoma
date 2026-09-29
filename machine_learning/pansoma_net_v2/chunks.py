"""Training batches read chunk by chunk, from a zstd-compressed copy of the shards.

On BeeGFS (/scratch) one random 161.6 KB pread takes ~26 ms under load (38 tensors/s per thread); a shard
streams at ~210 MB/s for one reader, but all readers of a node together got ~300 MB/s (tequila, 2026-09-28), which
is ~2,000 raw tensors/s. The tensors are 45 % zero bytes and zstd level 1 packs a 512-row chunk 28.8x
(decompression ~1.1 GB/s per thread), so the shards are copied once into a chunk cache (`build_cache`,
`python -m pansoma_net_v2.chunks`): per shard one file of compressed chunks of `chunk_rows` rows and an index of
their offsets. HG008 Illumina SNV (chr2-22): ~440 GB of shards, ~15 GB of chunks.

Each epoch ChunkLoader
- takes the epoch's tensors from an EpochSampler (every somatic / germline tensor, the sampled non);
- groups them by chunk and reads each chunk with one pread (the compressed chunk, or without a cache the raw rows
  from the first to the last chosen one) by a pool of `threads`, in a random chunk order;
- draws the batches from a shuffled pool of `window` chunks at a time (plus the rest of the previous pool), so
  a batch mixes tensors of ~window random places of the genome; the next window is read while the GPU trains.
With DDP the chunks are dealt to the ranks by load, and every rank yields the same number of tensors: a rank
short of its share repeats random ones of its own tensors, one with more drops its last ones.
"""
import argparse
import hashlib
import math
import os
import queue
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import torch
import zstandard

from .data import KINDS, TENSOR_SHAPE, TensorDataset, load_parts, read_shard_header

ROW_BYTES = int(np.prod(TENSOR_SHAPE))
ZSTD_LEVEL = 1


def cache_name(shard, chunk_rows):
    """The cache file stem of a shard: its name, and a hash of its real path, size, mtime and chunk_rows."""
    st = os.stat(shard)
    key = f"{os.path.realpath(shard)}|{st.st_size}|{st.st_mtime_ns}|{chunk_rows}|zstd{ZSTD_LEVEL}"
    return f"{Path(shard).name[:-len('_data.npy')]}.{hashlib.sha1(key.encode()).hexdigest()[:12]}"


def build_shard(shard, cache, chunk_rows):
    """Writes cache/<name>.zst (the compressed chunks) and .idx.npy (per chunk: offset, bytes, rows), each via a
    temporary file of this process and a rename; returns the stem."""
    stem = Path(cache) / cache_name(shard, chunk_rows)
    if Path(f"{stem}.idx.npy").exists():
        return stem
    offset, rows = read_shard_header(shard)
    tmp = f"{stem}.zst.tmp{os.getpid()}.{threading.get_ident()}"
    compressor = zstandard.ZstdCompressor(level=ZSTD_LEVEL)
    index, pos = [], 0
    fd = os.open(shard, os.O_RDONLY)
    try:
        with open(tmp, "wb") as out:
            for lo in range(0, rows, chunk_rows):
                buf = np.empty((min(chunk_rows, rows - lo), *TENSOR_SHAPE), np.int8)
                if os.preadv(fd, [buf], offset + lo * ROW_BYTES) != buf.nbytes:
                    raise IOError(f"{shard}: short read of rows {lo}-{lo + len(buf)}")
                z = compressor.compress(buf.data)
                out.write(z)
                index.append((pos, len(z), len(buf)))
                pos += len(z)
    finally:
        os.close(fd)
    os.replace(tmp, f"{stem}.zst")
    tmp_idx = f"{stem}.idx.tmp{os.getpid()}.{threading.get_ident()}.npy"
    np.save(tmp_idx, np.array(index, np.int64).reshape(-1, 3))
    os.replace(tmp_idx, f"{stem}.idx.npy")
    return stem


def build_cache(files, cache, chunk_rows=512, threads=16, log=print):
    """The chunk cache of every shard in `files` (built where missing, `threads` shards at a time)."""
    Path(cache).mkdir(parents=True, exist_ok=True)
    missing = [f for f in files if not Path(f"{Path(cache) / cache_name(f, chunk_rows)}.idx.npy").exists()]
    if missing:
        log(f"chunk cache {cache}: compressing {len(missing)} of {len(files)} shards with {threads} threads")
        t0, done = time.time(), 0
        with ThreadPoolExecutor(threads) as ex:
            for _ in ex.map(lambda f: build_shard(f, cache, chunk_rows), missing):
                done += 1
                if done % 10 == 0 or done == len(missing):
                    log(f"  {done}/{len(missing)} shards ({time.time() - t0:.0f} s)")
    return {f: Path(cache) / cache_name(f, chunk_rows) for f in files}


class ChunkLoader:
    def __init__(self, dataset, sampler, batch_size, threads=16, chunk_rows=512, window=32, rank=0, world=1,
                 seed=0, prefetch=4, pin_memory=None, cache=None):
        """dataset: a TensorDataset; sampler: an EpochSampler over it with rank 0 and world 1 (this loader splits
        the epoch over the ranks); cache: {shard file: cache stem} from build_cache (None: read the raw shards)."""
        self.ds, self.sampler, self.batch = dataset, sampler, batch_size
        self.threads, self.chunk_rows, self.window = threads, chunk_rows, window
        self.rank, self.world, self.seed, self.prefetch = rank, world, seed, prefetch
        self.pin = torch.cuda.is_available() if pin_memory is None else pin_memory
        self.headers = [None] * len(dataset.files)  # (row-0 offset, rows) per shard, read on first use
        self.cache = [cache.get(f) for f in dataset.files] if cache is not None else None  # None: not cached
        self.chunk_index = [None] * len(dataset.files)

    def set_epoch(self, epoch):
        self.sampler.set_epoch(epoch)

    def per_rank(self):
        return math.ceil(self.sampler.selection_size() / self.world)

    def __len__(self):
        return math.ceil(self.per_rank() / self.batch)

    def _header(self, f):
        if self.headers[f] is None:
            self.headers[f] = read_shard_header(self.ds.files[f])
        return self.headers[f]

    def plan(self):
        """(dataset indices sorted by shard and row, chunk start / end in them, this rank's chunks in reading
        order) for the current epoch."""
        sel = self.sampler.epoch_indices()
        order = np.lexsort((self.ds.row[sel], self.ds.shard[sel]))
        sel = sel[order]
        key = self.ds.shard[sel].astype(np.int64) * (1 << 32) + self.ds.row[sel] // self.chunk_rows
        starts = np.flatnonzero(np.r_[True, key[1:] != key[:-1]]) if len(sel) else np.zeros(0, np.int64)
        ends = np.r_[starts[1:], len(sel)].astype(np.int64)
        chunks = np.random.default_rng([self.seed, self.sampler.epoch, 7]).permutation(len(starts))
        if self.world > 1:  # deal the chunks (in their random order) to the least loaded rank
            load, mine = np.zeros(self.world), []
            for c in chunks.tolist():
                r = int(np.argmin(load))
                load[r] += ends[c] - starts[c]
                if r == self.rank:
                    mine.append(c)
            chunks = np.asarray(mine, np.int64)
        return sel, starts, ends, chunks

    def _read(self, idx):
        """x (len(idx), 8, 200, 101) of the dataset indices idx (one shard, one chunk, ascending rows)."""
        f = int(self.ds.shard[idx[0]])
        rows = self.ds.row[idx].astype(np.int64)
        if self.cache is not None and self.cache[f] is not None:
            return self._read_cached(f, rows)
        offset, _ = self._header(f)
        lo, hi = int(rows[0]), int(rows[-1]) + 1
        buf = np.empty((hi - lo, *TENSOR_SHAPE), np.int8)
        fd = os.open(self.ds.files[f], os.O_RDONLY)
        try:
            got = os.preadv(fd, [buf], offset + lo * ROW_BYTES)
        finally:
            os.close(fd)
        if got != buf.nbytes:
            raise IOError(f"{self.ds.files[f]}: short read of rows {lo}-{hi}")
        return buf if len(rows) == hi - lo else buf[rows - lo]

    def _read_cached(self, f, rows):
        stem = self.cache[f]
        if self.chunk_index[f] is None:
            index = np.load(f"{stem}.idx.npy")
            if (index[:-1, 2] != self.chunk_rows).any() or (len(index) and index[-1, 2] > self.chunk_rows):
                raise ValueError(f"{stem}: chunks of {index[0, 2]} rows, not {self.chunk_rows}")
            self.chunk_index[f] = index
        c = int(rows[0]) // self.chunk_rows
        pos, size, n = (int(v) for v in self.chunk_index[f][c])
        z = bytearray(size)
        fd = os.open(f"{stem}.zst", os.O_RDONLY)
        try:
            got = os.preadv(fd, [z], pos)
        finally:
            os.close(fd)
        if got != size:
            raise IOError(f"{stem}.zst: short read of chunk {c}")
        raw = zstandard.ZstdDecompressor().decompress(z, max_output_size=n * ROW_BYTES)
        if len(raw) != n * ROW_BYTES:
            raise IOError(f"{stem}.zst: chunk {c} has {len(raw)} bytes, not {n} rows")
        chunk = np.frombuffer(raw, np.int8).reshape(n, *TENSOR_SHAPE)
        return chunk[rows - c * self.chunk_rows].copy() if len(rows) < n else chunk.copy()

    def _batches(self, out, stop):
        """Producer: puts (x, blocks, scalars, y) batches on `out`, then None (or the exception)."""
        pool_exec = ThreadPoolExecutor(self.threads, thread_name_prefix="chunk-read")
        try:
            sel, starts, ends, chunks = self.plan()
            rng = np.random.default_rng([self.seed, self.sampler.epoch, 11 + self.rank])
            windows = [chunks[i:i + self.window] for i in range(0, len(chunks), self.window)]

            def submit(w):
                return [(sel[starts[c]:ends[c]], pool_exec.submit(self._read, sel[starts[c]:ends[c]])) for c in w]
            pending = [submit(w) for w in windows[:2]]
            quota, sent = self.per_rank(), 0
            left_idx, left_x = np.zeros(0, np.int64), []
            for k in range(len(windows)):
                reads = pending.pop(0)
                if k + 2 < len(windows):
                    pending.append(submit(windows[k + 2]))
                idx = np.concatenate([left_idx] + [i for i, _ in reads])
                parts = left_x + [fut.result() for _, fut in reads]
                # (part, row in part) of every pooled tensor, shuffled
                where = np.concatenate([np.stack([np.full(len(p), j), np.arange(len(p))], 1) for j, p in enumerate(parts)])
                perm = rng.permutation(len(idx))
                idx, where = idx[perm], where[perm]
                n_full = min(len(idx) // self.batch, max(0, quota - sent) // self.batch)
                for b in range(n_full):
                    if stop.is_set():
                        return
                    s = slice(b * self.batch, (b + 1) * self.batch)
                    out.put(self._assemble(idx[s], where[s], parts))
                    sent += self.batch
                rest = slice(n_full * self.batch, len(idx))
                left_idx = idx[rest]
                left_x = [self._gather(where[rest], parts)] if len(left_idx) else []
            # the end: the rest of the pool, filled up (or cut) to this rank's quota
            need = quota - sent
            if need > 0 and len(left_idx) < need:
                mine = sel[np.concatenate([np.arange(starts[c], ends[c]) for c in chunks])] if len(chunks) else sel
                extra = rng.choice(mine, need - len(left_idx))
                left_idx = np.concatenate([left_idx, extra])
                left_x = left_x + [np.stack([self._read(np.array([i]))[0] for i in extra])] if len(extra) else left_x
            left_idx = left_idx[:max(0, need)]
            if len(left_idx):
                x = np.concatenate(left_x)[:len(left_idx)]
                for b in range(0, len(left_idx), self.batch):
                    out.put(self._pack(left_idx[b:b + self.batch], x[b:b + self.batch]))
            out.put(None)
        except BaseException as e:  # noqa: BLE001 -- handed to the consumer
            out.put(e)
        finally:
            pool_exec.shutdown(wait=False, cancel_futures=True)

    @staticmethod
    def _gather(where, parts):
        x = np.empty((len(where), *TENSOR_SHAPE), np.int8)
        for j in np.unique(where[:, 0]).tolist():
            m = where[:, 0] == j
            x[m] = parts[j][where[m, 1]]
        return x

    def _assemble(self, idx, where, parts):
        return self._pack(idx, self._gather(where, parts))

    def _pack(self, idx, x):
        ds = self.ds
        t = (torch.from_numpy(x), torch.from_numpy(ds.blocks[idx]), torch.from_numpy(ds.scalars[idx]),
             torch.from_numpy(ds.label[idx].astype(np.int64)))
        return tuple(v.pin_memory() for v in t) if self.pin else t

    def __iter__(self):
        out, stop = queue.Queue(maxsize=self.prefetch), threading.Event()
        producer = threading.Thread(target=self._batches, args=(out, stop), daemon=True, name="chunk-batches")
        producer.start()
        try:
            while True:
                item = out.get()
                if item is None:
                    return
                if isinstance(item, BaseException):
                    raise item
                yield item
        finally:
            stop.set()
            while producer.is_alive():  # unblock a producer waiting on a full queue
                try:
                    out.get(timeout=0.1)
                except queue.Empty:
                    pass


def main(argv=None):
    p = argparse.ArgumentParser(description="Build the zstd chunk cache of merged tensor sets (training chromosomes)")
    p.add_argument("--tensors", nargs="+", required=True)
    p.add_argument("--kinds", nargs="+", default=list(KINDS), choices=KINDS)
    p.add_argument("--chroms", nargs="+", help="default: the autosomes but chr1 (train's default training chromosomes)")
    p.add_argument("--index-cache", required=True, help="the index cache (as train --cache-dir)")
    p.add_argument("--cache", required=True, help="chunk cache directory (as train --chunk-cache)")
    p.add_argument("--chunk-rows", type=int, default=512)
    p.add_argument("--threads", type=int, default=16)
    args = p.parse_args(argv)
    for index, pos in load_parts(args.tensors, args.kinds, args.index_cache, labelled=False):
        chroms = args.chroms or [c for c in index.meta["chroms"] if c != "chr1" and c[3:].isdigit()]
        ds = TensorDataset([(index, index.select(chroms, labelled=False))])
        files = [ds.files[f] for f in np.unique(ds.shard).tolist()]
        build_cache(files, args.cache, args.chunk_rows, args.threads, lambda m: print(m, flush=True))
        print(f"{index.dir}: {len(files)} shards cached", flush=True)


if __name__ == "__main__":
    main()
