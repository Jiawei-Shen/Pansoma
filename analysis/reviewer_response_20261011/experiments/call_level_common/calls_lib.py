"""Shared code of the call-level analyses E4b, E19, E10 and E2.4 (reviewer response 2026-10-11).

Inputs are read only:
- call sets: <RUNS>/<run>/vcf_chr1*/ written by machine_learning/pansoma_net_v2/jobs/vcf.sh (graph VCF -> linear VCF ->
  PoN -> rtg vcfeval in the benchmark BED (eval_bed) and on the whole chromosome (eval_nobed));
- the tensor set's labels.manifest.json (somatic truth VCF / BED, germline truth VCF / BED, FASTA), chr1 labels and
  variant_summary;
- HPRC v1.1 haplotype alleles (vcfbub a100k + vcfwave VCF), the four PoN VCFs, local stratification BEDs.

Caches and per-call outputs go to DATA. Nothing in the pipeline is modified; RepeatContext / tract / minimal_period
(indexed_gam_pipeline_v4) and indel_min_af (scripts/filter_panel_of_normals.py) are imported.
"""
from __future__ import annotations

import bisect
import gzip
import json
import math
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

REPO = Path("/scratch/jshen/Github/Pansoma")
RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013")
TENSORS = Path("/scratch/jshen/data/pansoma_v2_tensors")
DATA = Path("/scratch/jshen/data/reviewer_experiments_20261013/call_level")
PON_DIR = Path("/scratch/jshen/data/Pansoma/panel_of_normal_VCFs")
PON_FILES = {"gnomAD": PON_DIR / "af-only-gnomad.hg38.vcf.gz",
             "dbSNP": PON_DIR / "Homo_sapiens_assembly38.dbsnp138.vcf.gz",
             "1000G": PON_DIR / "1000g_pon.hg38.vcf.gz",
             "CoLoRSdb": PON_DIR / "CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz"}
HPRC_VCF = Path("/scratch/jshen/data/hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz")
HPRC_TABLE = DATA / "hprc_v1.1_chr1_alleles.norm.tsv.gz"
SMAHT_DIR = Path("/scratch/qfu/reference/hg38_DAC_stratification")
E9A_DIR = Path("/scratch/jshen/data/reviewer_experiments_20261013/E9a_unplaced_projection")
CHROM = "chr1"
BUILD_AF_FLOOR = 0.08  # tensor build AF (pansoma_v2_tensors/README.txt: build AF 0.08 for SNV and INDEL)

sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ------------------------------------------------------------------ call sets
SAMPLES = {"HG008": "HG008T", "COLO829": "COLO829T"}


def set_parts(set_name):
    """(sample, platform, arm, purity) of a tensor set name, e.g. COLO829T_fiberseq_P25 -> COLO829T fiberseq d9 P25."""
    toks = set_name.split("_")
    sample, platform = toks[0], toks[1]
    arm = "d46" if set_name.endswith("_d46") else "linear" if set_name.endswith("_linear") else "d9"
    m = re.search(r"_P(\d+)$", set_name)
    purity = f"P{m.group(1)}" if m else "P100"
    return sample, platform, arm, purity


@dataclass
class CallSet:
    run: str
    vcf_dir: Path
    group: str          # insample | xs | joint
    variant: str        # main | norescue | valbest
    model: str          # source model run (RIS name)
    model_sample: str
    model_platform: str
    kind: str
    set: str
    sample: str         # sample of the tensor set (target)
    platform: str
    arm: str
    purity: str
    prefix: str         # <set>.<KIND>

    @property
    def id(self):
        return f"{self.run}__{self.vcf_dir.name}"

    @property
    def label(self):
        s = f"{self.sample} {self.platform} {self.kind} {self.arm}"
        if self.purity != "P100":
            s += f" {self.purity}"
        if self.group == "xs":
            s += f" (model {self.model_sample} {self.model_platform})"
        if self.group == "joint":
            s += " (joint model)"
        if self.variant != "main":
            s += f" [{self.variant}]"
        return s

    def path(self, ev, part, which="pon_calls"):
        return self.vcf_dir / ev / "vcfeval" / which / f"{part}.vcf.gz"

    def report(self, ev):
        return json.loads((self.vcf_dir / ev / "report.json").read_text())

    @property
    def linear_vcf(self):
        return self.vcf_dir / f"{self.prefix}.linear.vcf.gz"

    @property
    def unplaced_vcf(self):
        return self.vcf_dir / f"{self.prefix}.linear.unplaced.vcf.gz"

    @property
    def stats(self):
        p = self.vcf_dir / f"{self.prefix}.linear.stats.json"
        return json.loads(p.read_text()) if p.exists() else {}


def complete(vcf_dir):
    """vcf.sh finished: both reports and every rtg directory done."""
    for ev in ("eval_bed", "eval_nobed"):
        if not (vcf_dir / ev / "report.json").is_file():
            return False
        for which in ("raw_calls", "pon_calls"):
            if not (vcf_dir / ev / "vcfeval" / which / "done").exists():
                return False
    return True


def discover(groups=("insample", "xs", "joint"), variants=("main", "norescue", "valbest")):
    """Every finished call set under RUNS (see the README of each analysis for the naming rules)."""
    out = []
    for d in sorted(RUNS.iterdir()):
        if not d.is_dir():
            continue
        run = d.name
        variant = "main"
        base = run
        if run.endswith("__norescue"):
            variant, base = "norescue", run[: -len("__norescue")]
        if base.startswith("xs__"):
            group, model = "xs", base.split("__")[1]
        elif "_joint_" in base:
            group, model = "joint", base
        else:
            group, model = "insample", base
        toks = model.split("_")
        model_sample, model_platform, kind = SAMPLES.get(toks[0], toks[0]), toks[1], toks[2]
        if kind not in ("SNV", "INDEL"):
            continue
        for v in sorted(d.glob("vcf_chr1*")):
            if not v.is_dir() or not complete(v):
                continue
            suffix = v.name[len("vcf_chr1"):].lstrip("_")
            var = variant
            if suffix == "valbest":
                var = "valbest"
            lin = sorted(v.glob(f"*.{kind}.linear.vcf.gz"))
            if len(lin) != 1:
                continue
            prefix = lin[0].name[: -len(".linear.vcf.gz")]
            set_name = prefix[: -len(f".{kind}")]
            if suffix and suffix != "valbest" and suffix != set_name:
                continue
            if group not in groups or var not in variants:
                continue
            sample, platform, arm, purity = set_parts(set_name)
            out.append(CallSet(run=run, vcf_dir=v, group=group, variant=var, model=model, model_sample=model_sample,
                               model_platform=model_platform, kind=kind, set=set_name, sample=sample,
                               platform=platform, arm=arm, purity=purity, prefix=prefix))
    return out


def manifest(set_name, kind):
    m = json.loads((TENSORS / set_name / "tensors" / kind / "labels.manifest.json").read_text())
    t = m["truth"]
    return dict(truth_vcf=t["somatic"]["vcf"], truth_bed=t["somatic"]["bed"],
                germline_vcf=t.get("germline", {}).get("vcf"), germline_bed=t.get("germline", {}).get("bed"),
                fasta=m["fasta"], rules=m.get("rules_sha256", "")[:12])


ARM_ORDER = {"linear": 0, "d9": 1, "d46": 2}
PLATFORM_ORDER = {"Illumina": 0, "PacBio": 1, "fiberseq": 1, "ONT": 2}
KIND_ORDER = {"SNV": 0, "INDEL": 1}
GROUP_ORDER = {"insample": 0, "joint": 1, "xs": 2}
VARIANT_ORDER = {"main": 0, "norescue": 1, "valbest": 2}
PURITY_ORDER = {"P100": 0, "P75": 1, "P50": 2, "P25": 3}


def sort_key(cs):
    return (GROUP_ORDER[cs.group], KIND_ORDER[cs.kind], cs.sample, PLATFORM_ORDER.get(cs.platform, 9),
            cs.model_sample, PLATFORM_ORDER.get(cs.model_platform, 9), ARM_ORDER[cs.arm],
            PURITY_ORDER.get(cs.purity, 9), VARIANT_ORDER[cs.variant])


# ------------------------------------------------------------------ reference and alleles
_SEQ = {}


def chrom_seq(fasta, chrom=CHROM):
    key = (str(fasta), chrom)
    if key not in _SEQ:
        import pysam
        with pysam.FastaFile(str(fasta)) as fa:
            _SEQ[key] = fa.fetch(chrom).upper()
    return _SEQ[key]


def normalize(seq, pos, ref, alt):
    """Left-aligned, parsimonious (pos 1-based, ref, alt); None for symbolic / identical alleles. `seq` = the
    chromosome, upper case."""
    ref, alt = ref.upper(), alt.upper()
    if not alt or alt in ("*", ".") or alt.startswith("<") or "[" in alt or "]" in alt or ref == alt:
        return None
    while True:
        changed = False
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
            changed = True
        if not ref or not alt:
            if pos <= 1:
                break
            b = seq[pos - 2]
            ref, alt, pos = b + ref, b + alt, pos - 1
            changed = True
        if not changed:
            break
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def keys_of(seq, pos, ref, alt):
    """Normalized key plus, for an MNP, its single-base substitutions."""
    n = normalize(seq, pos, ref, alt)
    if n is None:
        return []
    out = [n]
    p, r, a = n
    if len(r) == len(a) > 1:
        out += [(p + i, r[i], a[i]) for i in range(len(r)) if r[i] != a[i]]
    return out


def interval(pos, ref, alt):
    """0-based half-open reference interval of a normalized allele: substitutions their bases, deletions the deleted
    bases, insertions the anchor and the next base, complex the REF bases."""
    if len(ref) == len(alt):
        return pos - 1, pos - 1 + len(ref)
    if len(alt) == 1 and ref[0] == alt[0]:
        return pos, pos + len(ref) - 1
    if len(ref) == 1 and alt[0] == ref[0]:
        return pos - 1, pos + 1
    return pos - 1, pos - 1 + len(ref)


def allele_kind(ref, alt):
    return "SNV" if len(ref) == len(alt) else "INDEL"


# ------------------------------------------------------------------ VCF reading
@dataclass
class Rec:
    pos: int
    ref: str
    alts: list
    id: str
    qual: str
    filt: str
    info: str
    fmt: dict = field(default_factory=dict)

    def info_get(self, key):
        for item in self.info.split(";"):
            k, _, v = item.partition("=")
            if k == key:
                return v
        return None


def read_records(path, chrom=CHROM, with_format=True):
    """Records of one chromosome (plain gzip scan; stops after the chromosome when the file is sorted)."""
    out = []
    seen = False
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t")
            if c[0] != chrom:
                if seen:
                    break
                continue
            seen = True
            fmt = {}
            if with_format and len(c) > 9:
                fmt = dict(zip(c[8].split(":"), c[9].split(":")))
            out.append(Rec(int(c[1]), c[3], c[4].split(","), c[2], c[5], c[6], c[7], fmt))
    return out


def rtg_parts(cs, ev="eval_bed", which="pon_calls"):
    """{'tp', 'fp', 'fn', 'tp-baseline'}: records of one rtg run (PASS calls after the PoN at the threshold for
    pon_calls; raw_calls without the PoN)."""
    return {p: read_records(cs.path(ev, p, which)) for p in ("tp", "fp", "fn", "tp-baseline")}


def call_af(rec):
    v = rec.fmt.get("AF")
    try:
        return float(v.split(",")[0])
    except (AttributeError, ValueError):
        return None


# ------------------------------------------------------------------ intervals
class Intervals:
    """Sorted, merged chromosome intervals (BED)."""

    def __init__(self, pairs):
        pairs = sorted(pairs)
        merged = []
        for a, b in pairs:
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        self.starts = [a for a, _ in merged]
        self.ends = [b for _, b in merged]
        self.bp = sum(b - a for a, b in merged)

    @classmethod
    def from_bed(cls, path, chrom=CHROM):
        pairs = []
        op = gzip.open if str(path).endswith(".gz") else open
        with op(path, "rt") as f:
            for line in f:
                if line.startswith(chrom + "\t"):
                    c = line.split("\t", 3)
                    pairs.append((int(c[1]), int(c[2])))
        return cls(pairs)

    def overlaps(self, a, b):
        i = bisect.bisect_right(self.starts, b - 1) - 1
        return i >= 0 and self.ends[i] > a

    def contains_point(self, x):
        i = bisect.bisect_right(self.starts, x) - 1
        return i >= 0 and self.ends[i] > x


_BEDS = {}


def bed(path):
    if path not in _BEDS:
        import hashlib
        tag = hashlib.sha1(str(Path(path).resolve()).encode()).hexdigest()[:8]
        cache = DATA / "beds" / (Path(path).name.replace(".bed", "") + f".{tag}.{CHROM}.json")
        if cache.exists():
            pairs = json.loads(cache.read_text())
            _BEDS[path] = Intervals([tuple(p) for p in pairs])
        else:
            iv = Intervals.from_bed(path)
            cache.parent.mkdir(parents=True, exist_ok=True)
            cache.write_text(json.dumps(list(zip(iv.starts, iv.ends))))
            _BEDS[path] = iv
    return _BEDS[path]


# ------------------------------------------------------------------ germline truth of the normal
class Germline:
    """chr1 alleles of the normal's germline truth (dipcall / DeepVariant): every allele a genotype carries, any FILTER.
    keys: normalized allele (+ MNP single bases) -> FILTER; records: one interval per record (union of its alleles)."""
    LONG = 1000

    def __init__(self, vcf, fasta):
        self.vcf = vcf
        seq = chrom_seq(fasta)
        self.keys = {}
        ivs = []
        for r in read_records(vcf):
            gt = re.split(r"[|/]", r.fmt.get("GT", "."))
            idx = {int(g) for g in gt if g not in (".", "") and g.isdigit() and int(g) > 0}
            if not idx:
                continue
            lo = hi = None
            for i in sorted(idx):
                if i - 1 >= len(r.alts):
                    continue
                ks = keys_of(seq, r.pos, r.ref, r.alts[i - 1])
                if not ks:
                    continue
                for k in ks:
                    self.keys.setdefault(k, r.filt)
                a, b = interval(*ks[0])
                lo = a if lo is None else min(lo, a)
                hi = b if hi is None else max(hi, b)
            if lo is not None:
                ivs.append((lo, hi, r.filt))
        ivs.sort()
        short = [(a, b) for a, b, _ in ivs if b - a <= self.LONG]
        self.s_start = np.array([a for a, _ in short], dtype=np.int64)
        self.s_end = np.array([b for _, b in short], dtype=np.int64)
        self.long = [(a, b) for a, b, _ in ivs if b - a > self.LONG]
        self.n_records = len(ivs)

    def _hits(self, a, b):
        """Records intersecting [a, b): (starts, ends) arrays."""
        lo = np.searchsorted(self.s_start, a - self.LONG, "left")
        hi = np.searchsorted(self.s_start, b, "left")
        s, e = self.s_start[lo:hi], self.s_end[lo:hi]
        m = e > a
        s, e = s[m], e[m]
        if self.long:
            ls = [(x, y) for x, y in self.long if x < b and y > a]
            if ls:
                s = np.concatenate([s, np.array([x for x, _ in ls])])
                e = np.concatenate([e, np.array([y for _, y in ls])])
        return s, e

    def overlapping(self, a, b):
        s, e = self._hits(a, b)
        return len(s)

    def count_near(self, a, b, w):
        """Records intersecting [a - w, b + w) that do not intersect [a, b) itself."""
        s, e = self._hits(a - w, b + w)
        own = (s < b) & (e > a)
        return int((~own).sum())

    def nearest(self, a, b, cap=100000):
        """Bases between [a, b) and the nearest record not intersecting it (None beyond cap)."""
        w = 16
        while w <= cap:
            s, e = self._hits(a - w, b + w)
            own = (s < b) & (e > a)
            s, e = s[~own], e[~own]
            if len(s):
                gaps = np.maximum(s - b, a - e)
                return int(gaps.min())
            w *= 4
        return None


_GERM = {}


def germline(vcf, fasta):
    if vcf not in _GERM:
        log(f"germline: loading {vcf}")
        _GERM[vcf] = Germline(vcf, fasta)
        log(f"germline: {_GERM[vcf].n_records} records, {len(_GERM[vcf].keys)} allele keys")
    return _GERM[vcf]


# ------------------------------------------------------------------ repeat context (pipeline classifier)
class Repeats:
    """GRCh38 repeat class of an allele. INDEL: indexed_gam_pipeline_v4 RepeatContext.classify on the event bases
    (shared prefix / suffix removed): HP>=7, HP4-6, STR, VNTR, none, complex. SNV / MNP (RepeatContext classifies
    INS / DEL only): the module's tract() of period 1-6 touching the substituted bases, with the same thresholds
    (HP>=7: period-1 run >= 7 bp; HP4-6; STR: period 2-6, >= 3 units; else none)."""

    def __init__(self, fasta):
        from indexed_gam_pipeline_v4.tensor_postprocessing.truth_labels import RepeatContext, tract
        self.rc = RepeatContext(fasta)
        self.tract = tract
        self.seq = chrom_seq(fasta)
        self.cache = {}

    def classify(self, pos, ref, alt):
        key = (pos, ref, alt)
        if key in self.cache:
            return self.cache[key]
        ref, alt = ref.upper(), alt.upper()
        i = 0
        while i < min(len(ref), len(alt)) and ref[i] == alt[i]:
            i += 1
        r, a = ref[i:], alt[i:]
        while r and a and r[-1] == a[-1]:
            r, a = r[:-1], a[:-1]
        pos0 = pos - 1 + i
        if len(r) == len(a):  # substitution
            lo, hi = pos0, pos0 + max(1, len(r))
            start = max(0, lo - 300)
            win = self.seq[start:hi + 300]
            l0, h0 = lo - start, hi - start
            hp = self.tract(win, l0, h0, 1)
            if hp[1] - hp[0] >= 7:
                c = "HP>=7"
            else:
                c = None
                for p in range(2, 7):
                    t = self.tract(win, l0, h0, p)
                    if t[1] - t[0] >= 3 * p:
                        c = "STR"
                        break
                if c is None:
                    c = "HP4-6" if hp[1] - hp[0] >= 4 else "none"
        elif r and a:
            c = "complex"
        else:
            c = self.rc.classify(CHROM, pos0, r, a)
        self.cache[key] = c
        return c


def repeat_group(c):
    return "HP>=7" if c == "HP>=7" else "STR" if c == "STR" else "other"


# ------------------------------------------------------------------ HPRC v1.1 haplotype alleles
class HPRC:
    """Normalized chr1 alleles with AC > 0 among the HPRC v1.1 haplotypes (44 diploid samples + CHM13; vcfbub a100k +
    vcfwave decomposition), table built by prep_hprc.py."""

    def __init__(self, table=HPRC_TABLE):
        import pysam
        if not Path(table).exists():
            raise SystemExit(f"{table} missing: run call_level_common/prep_hprc.py")
        self.t = pysam.TabixFile(str(table))

    def lookup(self, key, window=0):
        """AC of the allele (0 when absent)."""
        p, r, a = key
        best = 0
        for line in self.t.fetch(CHROM, max(0, p - 1 - window), p + window):
            c = line.split("\t")
            if int(c[1]) == p and c[2] == r and c[3] == a:
                best = max(best, int(c[4]))
        return best


# ------------------------------------------------------------------ population resources (PoN VCFs)
class Population:
    """Allele lookups in the four PoN VCFs, matched after normalization (the PoN script itself matches POS / REF /
    ALT exactly, so a normalized match can find an allele in another representation)."""

    def __init__(self, fasta):
        import pysam
        self.v = {k: pysam.VariantFile(str(p)) for k, p in PON_FILES.items()}
        self.contig = {k: (CHROM if CHROM in v.header.contigs else CHROM[3:]) for k, v in self.v.items()}
        self.seq = chrom_seq(fasta)
        from filter_panel_of_normals import indel_min_af, SNV_MIN_AF
        import pysam as _p
        self.fa = _p.FastaFile(str(fasta))
        self._indel_min_af = indel_min_af
        self.snv_floor = SNV_MIN_AF

    def floor(self, pos, ref, alt):
        """The PoN's gnomAD / CoLoRSdb AF floor for this call allele (raw VCF representation)."""
        if len(ref) == len(alt):
            return self.snv_floor
        return self._indel_min_af(self.fa, CHROM, pos, ref, alt)

    def lookup(self, pos, ref, alt):
        """{resource: dict(af=max AF of matching alleles or None, exact=bool, somatic_only=bool)} for resources that
        carry the allele."""
        key = normalize(self.seq, pos, ref, alt)
        if key is None:
            return {}
        a0, b0 = interval(*key)
        w = 10 if len(ref) == len(alt) else 200
        lo, hi = min(a0, pos - 1) - w, max(b0, pos - 1 + len(ref)) + w
        out = {}
        for name, v in self.v.items():
            hit = None
            for rec in v.fetch(self.contig[name], max(0, lo), hi):
                alts = rec.alts or ()
                for k, alt2 in enumerate(alts):
                    if alt2 is None:
                        continue
                    ks = keys_of(self.seq, rec.pos, rec.ref, alt2)
                    if key not in ks:
                        continue
                    exact = rec.pos == pos and rec.ref.upper() == ref.upper() and alt2.upper() == alt.upper()
                    af = None
                    if "AF" in rec.info:
                        x = rec.info["AF"]
                        x = x if isinstance(x, tuple) else (x,) * len(alts)
                        af = x[k] if k < len(x) else None
                    sao = rec.info.get("SAO") if "SAO" in rec.header.info else None
                    h = hit or dict(af=None, exact=False, somatic_only=True)
                    if af is not None:
                        h["af"] = af if h["af"] is None else max(h["af"], af)
                    h["exact"] = h["exact"] or exact
                    h["somatic_only"] = h["somatic_only"] and sao == 2
                    hit = h
            if hit:
                out[name] = hit
        return out


# ------------------------------------------------------------------ metrics
def prf(tpb, fn, tpc, fp):
    r = tpb / (tpb + fn) if tpb + fn else None
    p = tpc / (tpc + fp) if tpc + fp else None
    f = 2 * p * r / (p + r) if p is not None and r is not None and p + r > 0 else (0.0 if p is not None and r is not None else None)
    return p, r, f


def fmt(x, d=3):
    if x is None:
        return "-"
    if isinstance(x, float):
        if math.isnan(x):
            return "-"
        return f"{x:.{d}f}"
    return str(x)


def pct(n, d):
    return "-" if not d else f"{n} ({100 * n / d:.1f} %)"


def binom_two_sided(b, c):
    """Exact two-sided McNemar p (binomial, p = 0.5) for discordant counts b, c."""
    n = b + c
    if n == 0:
        return None
    k = min(b, c)
    s = sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n
    return min(1.0, 2 * s)


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(fmt(x) if not isinstance(x, str) else x for x in r) + " |" for r in rows]
    return out


def write_tsv_gz(path, header, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp{os.getpid()}")
    with gzip.open(tmp, "wt") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join("" if x is None else str(x) for x in r) + "\n")
    tmp.replace(path)


# ------------------------------------------------------------------ truth alleles with a tensor
def truth_tensor_keys(set_name, kind):
    """(vcf_pos, vcf_ref, vcf_alt) of the chr1 truth alleles that some tensor of the set matches exactly (the labels'
    `somatic` entries; partial matches not counted). Cached in DATA/tensor_keys."""
    lab = TENSORS / set_name / "tensors" / kind / f"{CHROM}_labels.ndjson"
    cache = DATA / "tensor_keys" / f"{set_name}.{kind}.{CHROM}.tsv"
    if cache.exists() and cache.stat().st_mtime > lab.stat().st_mtime:
        with open(cache) as f:
            return {(int(p), r, a) for p, r, a in (line.rstrip("\n").split("\t") for line in f)}
    keys = set()
    with open(lab) as f:
        for line in f:
            if '"somatic": [{' not in line:
                continue
            for e in json.loads(line)["somatic"]:
                keys.add((int(e["vcf_pos"]), e["vcf_ref"], e["vcf_alt"]))
    cache.parent.mkdir(parents=True, exist_ok=True)
    tmp = cache.with_suffix(f".tmp{os.getpid()}")
    tmp.write_text("".join(f"{p}\t{r}\t{a}\n" for p, r, a in sorted(keys)))
    tmp.replace(cache)
    return keys


def has_tensor(keys, rec):
    return any((rec.pos, rec.ref, alt) in keys for alt in rec.alts)
