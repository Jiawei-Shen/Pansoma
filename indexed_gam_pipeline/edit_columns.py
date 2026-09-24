"""Edit-sized storage with bounded materialization of alignment-column windows."""
from bisect import bisect_right
from dataclasses import dataclass


@dataclass(slots=True)
class EditRun:
    start: int
    ref: str
    alt: str
    qualities: object
    op: str
    node: int
    cursor: int
    reverse: bool
    visit: int
    node_length: int

    def __len__(self):
        return max(len(self.ref), len(self.alt))

    def column(self, k):
        from indexed_gam_pipeline.candidates import Column
        f, t = len(self.ref), len(self.alt)
        boundary = k >= f
        pos = self.cursor + min(k, f)
        if self.reverse:
            pos = self.node_length - pos - (0 if boundary else 1)
        return Column(self.alt[k] if k < t else '-', self.ref[k] if k < f else '-',
                      self.qualities[k] if k < t else -1, self.op, self.node,
                      pos, self.reverse, self.visit, boundary)


class EditColumns:
    def __init__(self):
        self.runs = []
        self.starts = []
        self.length = 0
        self.materialized = 0

    def __len__(self):
        return self.length

    def add(self, ref, alt, qualities, op, node, cursor, reverse, visit, node_length):
        run = EditRun(self.length, ref, alt, qualities, op, node, cursor, reverse, visit, node_length)
        if len(run):
            self.starts.append(self.length)
            self.runs.append(run)
            self.length += len(run)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return [self[i] for i in range(*index.indices(len(self)))]
        if index < 0:
            index += len(self)
        if not 0 <= index < len(self):
            raise IndexError(index)
        run = self.runs[bisect_right(self.starts, index)-1]
        self.materialized += 1
        return run.column(index-run.start)

    def view(self, visit, context):
        lo = max(0, visit.first-context) if context is not None else 0
        hi = min(len(self), visit.last+context) if context is not None else len(self)
        return ColumnView(self, lo, hi, visit.reverse)


class ColumnView:
    def __init__(self, columns, lo, hi, reverse=False):
        self.columns, self.lo, self.hi, self.reverse = columns, lo, hi, reverse

    def __len__(self):
        return self.hi-self.lo

    def __getitem__(self, index):
        if isinstance(index, slice):
            lo, hi, step = index.indices(len(self))
            if step != 1:
                return [self[i] for i in range(lo, hi, step)]
            hi = max(lo, hi)
            return ColumnView(self.columns, self.hi-hi if self.reverse else self.lo+lo,
                              self.hi-lo if self.reverse else self.lo+hi, self.reverse)
        if index < 0:
            index += len(self)
        if not 0 <= index < len(self):
            raise IndexError(index)
        col = self.columns[self.hi-1-index if self.reverse else self.lo+index]
        return col.flipped() if self.reverse else col

    def _interval(self, run, k0, k1):
        lo, hi = max(self.lo, run.start+k0), min(self.hi, run.start+k1)
        if hi <= lo:
            return None
        return (self.hi-hi, self.hi-lo) if self.reverse else (lo-self.lo, hi-self.lo)

    def split(self, visit, candidate):
        """Find first/last matching columns using edit intervals, never base scans."""
        cols = self.columns
        ri = max(0, bisect_right(cols.starts, visit.first)-1)
        matches, cuts, local = [], [], []
        for index in range(ri, len(cols.runs)):
            run = cols.runs[index]
            if run.start >= visit.last:
                break
            if run.visit != visit.index:
                continue
            size, f = len(run), len(run.ref)
            interval = self._interval(run, 0, size)
            if interval is None:
                continue
            local.append(interval)
            # Reference positions decrease in read order for reverse mappings.
            p0 = run.node_length-run.cursor-1 if run.reverse else run.cursor
            boundary = run.node_length-run.cursor-f if run.reverse else run.cursor+f
            if candidate.kind == 'INS':
                match = self._interval(run, f, size) if boundary == candidate.start else None
            else:
                k0, k1 = ((p0-candidate.end+1, p0-candidate.start+1) if run.reverse
                          else (candidate.start-p0, candidate.end-p0))
                match = self._interval(run, max(0, k0), min(f, k1))
            if match:
                matches.append(match)
            # Fallback insertion point: first local column with pos >= start.
            k0, k1 = ((0, min(f, p0-candidate.start+1)) if run.reverse
                      else (max(0, candidate.start-p0), f))
            cut = self._interval(run, k0, k1)
            if cut:
                cuts.append(cut[0])
            if boundary >= candidate.start:
                cut = self._interval(run, f, size)
                if cut:
                    cuts.append(cut[0])
        if matches:
            lo, hi = min(x[0] for x in matches), max(x[1] for x in matches)
        else:
            lo = hi = min(cuts) if cuts else max((x[1] for x in local), default=0)
        return self[:lo], self[lo:hi], self[hi:]

    def reference_support(self, visit, candidate):
        """Check reference support from edit intervals without expanding insertions."""
        covered = 0
        cols = self.columns
        ri = max(0, bisect_right(cols.starts, visit.first)-1)
        for index in range(ri, len(cols.runs)):
            run = cols.runs[index]
            if run.start >= visit.last:
                break
            f = len(run.ref)
            boundary = run.node_length-run.cursor-f if run.reverse else run.cursor+f
            if len(run) > f and candidate.start < boundary < candidate.end:
                return False
            p0 = run.node_length-run.cursor-1 if run.reverse else run.cursor
            k0, k1 = ((p0-candidate.end+1, p0-candidate.start+1) if run.reverse
                      else (candidate.start-p0, candidate.end-p0))
            lo, hi = max(0, k0), min(f, k1)
            if hi > lo:
                if run.op not in ('M', 'X') or run.alt[lo:hi] != run.ref[lo:hi]:
                    return False
                covered += hi-lo
        return covered == len(candidate.ref)


def select_decode_mode(gam, requested='auto'):
    """Classify once per input from the first ten unfiltered GAM records."""
    if requested not in ('auto', 'full', 'window'):
        raise ValueError('Decode mode must be auto, full or window')
    lengths = []
    if requested == 'auto':
        from indexed_gam_pipeline.gam_reader import scan_gam
        from contextlib import closing
        from itertools import islice
        with closing(scan_gam(gam)) as records:
            lengths = [len(a.sequence) for a in islice(records, 10)]
    mean = sum(lengths)/len(lengths) if lengths else None
    selected = ('window' if mean is not None and mean >= 1000 else 'full') if requested == 'auto' else requested
    return dict(requested=requested, selected=selected, sample_lengths=lengths,
                mean_read_length=mean, threshold_bp=1000, statistic='arithmetic_mean',
                sample_scope='first 10 unfiltered GAM records; all records if fewer than 10')
