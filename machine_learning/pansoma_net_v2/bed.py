"""BED regions with the labeller's semantics (indexed_gam_pipeline_v4 tensor_postprocessing/truth_labels.py: Bed,
linear_interval and the confident-region test), so evaluation drops exactly the tensors a BED filter would."""
from collections import defaultdict

import numpy as np


class Bed:
    """Merged half-open intervals per chromosome."""

    def __init__(self, intervals):
        self.data = {}
        for chrom, items in intervals.items():
            merged = []
            for s, e in sorted(items):
                if merged and s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            self.data[chrom] = (np.array([m[0] for m in merged], dtype=np.int64),
                                np.array([m[1] for m in merged], dtype=np.int64))

    @classmethod
    def read(cls, path):
        intervals = defaultdict(list)
        with open(path) as stream:
            for line in stream:
                if line.strip() and not line.startswith(("#", "track", "browser")):
                    f = line.split("\t")
                    intervals[f[0]].append((int(f[1]), int(f[2])))
        return cls(intervals)

    def intersect(self, other):
        result = defaultdict(list)
        for chrom in self.data.keys() & other.data.keys():
            (a0, a1), (b0, b1) = self.data[chrom], other.data[chrom]
            i = j = 0
            while i < len(a0) and j < len(b0):
                s, e = max(a0[i], b0[j]), min(a1[i], b1[j])
                if s < e:
                    result[chrom].append((int(s), int(e)))
                if a1[i] < b1[j]:
                    i += 1
                else:
                    j += 1
        return Bed(result)

    def contains(self, chrom, s, e):
        if chrom not in self.data:
            return False
        starts, ends = self.data[chrom]
        k = int(np.searchsorted(starts, s, side="right")) - 1
        return bool(k >= 0 and ends[k] >= e)


def in_region(confident, event_type, grch38, anchor):
    """The labeller's test: a tensor on a GRCh38 node by its linear interval (an INS spans the bases around its
    insertion point), an off-reference tensor by the middle base of its anchor interval, else not in the region."""
    if grch38 is not None:
        pos0 = grch38["pos0"]
        s, e = (pos0 - 1, pos0 + 1) if event_type == "INS" else (pos0, pos0 + max(1, len(grch38["ref"])))
        return confident.contains(grch38["chrom"], s, e)
    if anchor is not None:
        m = (anchor["start"] + anchor["end"]) // 2
        return confident.contains(anchor["chrom"], m, m + 1)
    return False
