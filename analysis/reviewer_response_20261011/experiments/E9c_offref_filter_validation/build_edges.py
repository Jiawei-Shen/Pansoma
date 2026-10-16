"""Edge table of a GFA (its L lines only) as a CSR over node handles, for bounded local graph searches.

usage: build_edges.py GFA OUTDIR

The L lines are located by binary search on the byte offset (both GFAs are S lines, then L and P/W blocks:
d46 has S, L, W; d9 has S, P/W, L), then parsed in chunks with numpy. A handle is 2 * node + rev (rev = 1 for '-').
An L line a(oa) -> b(ob) gives the directed handle edges (a,oa) -> (b,ob) and (b,!ob) -> (a,!oa).

Outputs in OUTDIR: offsets.npy (int32, or int64 past 2^31 edges; 2 * (max_node + 1) + 1), targets.npy (int32 handles; int64 when a
handle reaches 2^31), meta.json (source, byte range, counts).
"""
import json
import os
import sys
import time

import numpy as np

def line_kind(f, off):
    f.seek(off)
    if off:
        f.readline()
    line = f.readline()
    return line[:1].decode() if line else ''


def find_first(f, lo, hi, pred):
    """Smallest line start in [lo, hi) whose line satisfies pred (pred monotone False..True over the range)."""
    while hi - lo > 1 << 16:
        mid = (lo + hi) // 2
        if pred(line_kind(f, mid)):
            hi = mid
        else:
            lo = mid
    f.seek(lo)
    if lo:
        f.readline()
    while True:
        pos = f.tell()
        line = f.readline()
        if not line or pred(line[:1].decode()):
            return pos


def l_range(path):
    size = os.path.getsize(path)
    with open(path, 'rb') as f:
        last = line_kind(f, max(0, size - 4096))
        if last == 'L':   # d9 layout: S, P/W, L (L block at the end)
            start = find_first(f, 0, size, lambda k: k == 'L')
            return start, size
        # d46 layout: S, L, W. L block between the last S line and the first P/W line
        start = find_first(f, 0, size, lambda k: k not in ('H', 'S'))
        end = find_first(f, start, size, lambda k: k in ('P', 'W'))
        return start, end


def main():
    gfa, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    start, end = l_range(gfa)
    print(f'L block bytes {start:,}-{end:,} ({(end - start) / 1e9:.2f} GB), located in {time.time() - t0:.1f}s', flush=True)

    # "L<TAB>a<TAB>+<TAB>b<TAB>-<TAB>0M" -> translate to whitespace-separated integers "a 0 b 1 0" and parse in C
    table = bytes.maketrans(b'L+-M\t\n', b' 01   ')
    src, srev, dst, drev = [], [], [], []
    n = 0
    with open(gfa, 'rb') as f:
        f.seek(start)
        left = end - start
        while left > 0:
            data = f.read(min(left, 1 << 28))
            left -= len(data)
            if left > 0:
                tail = f.readline()
                left -= len(tail)
                data += tail
            assert data.count(b'L') == data.count(b'\n'), 'non-L line inside the L block'
            v = np.fromstring(data.translate(table), dtype=np.int64, sep=' ')
            assert len(v) % 5 == 0 and not v[4::5].any(), 'unexpected overlap field (not 0M)'
            v = v.reshape(-1, 5)
            src.append(v[:, 0].copy()); srev.append(v[:, 1].astype(bool)); dst.append(v[:, 2].copy())
            drev.append(v[:, 3].astype(bool))
            n += len(v)
            print(f'  {n:,} L lines, {time.time() - t0:.0f}s', flush=True)
    a = np.concatenate(src); b = np.concatenate(dst); ar = np.concatenate(srev); br = np.concatenate(drev)
    del src, dst, srev, drev
    max_node = int(max(a.max(), b.max()))
    assert 2 * (max_node + 1) < 2 ** 30
    ha = 2 * a + ar; hb = 2 * b + br
    del a, b, ar, br
    frm = np.concatenate([ha, hb ^ 1]); to = np.concatenate([hb, ha ^ 1])
    del ha, hb
    key = np.unique((frm << 30) | to)          # handles < 2^30 (max node 2.1e8); sorts and drops duplicates
    del frm, to
    frm = key >> 30; to = key & ((1 << 30) - 1)
    del key
    counts = np.bincount(frm, minlength=2 * (max_node + 1))
    offsets = np.zeros(2 * (max_node + 1) + 1, np.int64 if len(to) >= 2 ** 31 else np.int32)
    offsets[1:] = np.cumsum(counts)
    del counts
    dtype = np.int32 if 2 * (max_node + 1) < 2 ** 31 else np.int64
    np.save(f'{out}/offsets.npy', offsets)
    np.save(f'{out}/targets.npy', to.astype(dtype))
    meta = dict(source=gfa, source_size=os.path.getsize(gfa), l_bytes=[start, end], l_lines=n, max_node=max_node,
                directed_handle_edges=int(len(to)), seconds=round(time.time() - t0, 1),
                handle='2 * node + (1 if reverse)', note='successors of handle h: targets[offsets[h]:offsets[h + 1]]')
    json.dump(meta, open(f'{out}/meta.json', 'w'), indent=1)
    print(json.dumps(meta), flush=True)


if __name__ == '__main__':
    main()
