"""Call-set discovery shared by E16 (block bootstrap), E17a (INDEL strata) and E18 (SBS96).

A "unit" is one chr1 call set scored by machine_learning/pansoma_net_v2/jobs/vcf.sh: a vcf_chr1*/ directory of a run
under /scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013/ that holds <set>.<KIND>.linear.vcf.gz and finished
rtg outputs. Each unit's PASS-after-PoN result is the rtg directory eval_<mode>/vcfeval/pon_calls/ (tp.vcf.gz =
call-side TPs, tp-baseline.vcf.gz = truth-side TPs, fp.vcf.gz, fn.vcf.gz); mode bed = inside the set's benchmark BED,
nobed = whole chr1. Nothing here writes to the run directories.

Run names (reviewer_cpu_20261013):
  <SAMPLE>_<PLATFORM>_<KIND>_{7fc7396 | d46_7fc7396 | linear_7fc7396}[...]   in-sample model (arm from the tensor set:
                                                                           <set>_d46 = HPRC v2.1 d46, <set>_linear =
                                                                           GRCh38-only graph, no suffix = HPRC v1.1 d9)
  *_joint_*                                                               COLO829T fiberseq joint model (4 purities),
                                                                           one vcf_chr1_<SET> per purity set
  xs__<model>__<target set>[__normds]                                      cross-sample prediction
  <run>__norescue                                                         SNV without predict's off-reference rescue
  <run>/vcf_chr1_valbest                                                  INDEL at the validation-chosen threshold
Comparators: rows of comparators.tsv (next to this file), rtg directories of score_comparator.py (E8) or any rtg
vcfeval output scored with the same truth / BED.
Environment (tests only): REVIEWER_RUNS overrides the run root, REVIEWER_COMPARATORS the comparator table.
"""
import csv
import gzip
import json
import os
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = Path(os.environ.get("REVIEWER_RUNS", "/scratch/jshen/data/pansoma_net_v2_runs/reviewer_cpu_20261013"))
TENSORS = Path("/scratch/jshen/data/pansoma_v2_tensors")
REPO = Path("/scratch/jshen/Github/Pansoma")
COMPARATORS = Path(os.environ.get("REVIEWER_COMPARATORS", HERE / "comparators.tsv"))
ARM_ORDER = {"linear": 0, "d9": 1, "d46": 2}
ARM_LABEL = {"linear": "GRCh38-only", "d9": "HPRC v1.1 d9", "d46": "HPRC v2.1 d46"}
ROLE_ORDER = {"insample": 0, "joint": 1, "xs": 2, "comparator": 3}
VARIANT_ORDER = {"standard": 0, "valbest": 1, "norescue": 2, "norescue+valbest": 3}
RTG_FILES = ("tp", "tp-baseline", "fp", "fn")


def arm_of(set_name):
    return "d46" if set_name.endswith("_d46") else "linear" if set_name.endswith("_linear") else "d9"


def base_set(set_name):
    for suffix in ("_d46", "_linear"):
        if set_name.endswith(suffix):
            return set_name[: -len(suffix)]
    return set_name


def sample_of(set_name):
    return set_name.split("_", 1)[0]


def manifest(set_name, kind):
    p = TENSORS / set_name / "tensors" / kind / "labels.manifest.json"
    return json.loads(p.read_text()) if p.exists() else None


@dataclass
class Unit:
    uid: str                      # <run>/<vcf dir> or comparator:<label>
    run: str
    set: str
    kind: str
    role: str                     # insample | joint | xs | comparator
    variant: str                  # standard | valbest | norescue | norescue+valbest
    model: str                    # the checkpoint's run (xs: the source model); comparator: its label
    normds: bool = False
    rtg: dict = field(default_factory=dict)       # mode -> rtg pon_calls directory
    report: dict = field(default_factory=dict)    # mode -> report.json (Pansoma) or {}
    truth: str = None
    bed: str = None
    truth_pon: Path = None        # PoN-tagged truth of the chromosome (eval_bed/truth.<KIND>.pon.vcf.gz)
    standard_uid: str = None      # the standard call set of the same predictions (variants only)
    note: str = ""

    @property
    def arm(self):
        return arm_of(self.set)

    @property
    def base_set(self):
        return base_set(self.set)

    @property
    def sample(self):
        return sample_of(self.set)

    @property
    def unplaced_pass(self):
        r = self.report.get("bed") or {}
        return r.get("unplaced_pass")

    @property
    def label(self):
        v = "" if self.variant == "standard" else f" [{self.variant}]"
        if self.role == "insample":
            return f"{self.base_set} {self.kind} {self.arm}{v}"
        if self.role == "joint":
            return f"{self.set} {self.kind} d9 joint-purity model{v}"
        if self.role == "xs":
            return f"{self.model.replace('_7fc7396', '')} -> {self.set}{' normds' if self.normds else ''}{v}"
        return f"{self.model} on {self.set} {self.kind}"

    def sort_key(self):
        return (ROLE_ORDER.get(self.role, 9), self.kind != "SNV", self.sample, self.base_set,
                ARM_ORDER.get(self.arm, 9), self.model, self.normds, VARIANT_ORDER.get(self.variant, 9))


def resolved(path):
    return str(Path(path).resolve()) if path else None


def _rtg_complete(w):
    """rtg finished: its 'done' marker and the four output VCFs."""
    return (w / "done").exists() and all((w / f"{k}.vcf.gz").exists() for k in RTG_FILES)


def _parse_run(name):
    """(role, model, variant, normds) of a run directory name."""
    parts = name.split("__")
    variant, normds = "standard", False
    if parts[-1] == "norescue":
        variant, parts = "norescue", parts[:-1]
    if parts[-1] == "normds":
        normds, parts = True, parts[:-1]
    if parts[0] == "xs" and len(parts) >= 3:
        return "xs", parts[1], variant, normds, "__".join(parts)
    if "_joint_" in name:
        return "joint", "__".join(parts), variant, normds, "__".join(parts)
    return "insample", "__".join(parts), variant, normds, "__".join(parts)


def discover(kinds=("SNV", "INDEL"), comparators=True):
    """Every finished unit (eval_bed rtg complete; eval_nobed when present), sorted for display."""
    units = []
    if RUNS.exists():
        for run_dir in sorted(p for p in RUNS.iterdir() if p.is_dir()):
            role, model, variant0, normds, base_run = _parse_run(run_dir.name)
            for vdir in sorted(run_dir.glob("vcf_chr1*")):
                linear = sorted(vdir.glob("*.linear.vcf.gz"))
                if len(linear) != 1 or not vdir.is_dir():
                    continue
                set_name, kind = linear[0].name.split(".")[:2]
                if kind not in kinds:
                    continue
                u = Unit(uid=f"{run_dir.name}/{vdir.name}", run=run_dir.name, set=set_name, kind=kind, role=role,
                         variant=variant0, model=model, normds=normds)
                if vdir.name.startswith("vcf_chr1_valbest"):
                    u.variant = "valbest" if variant0 == "standard" else f"{variant0}+valbest"
                for mode in ("bed", "nobed"):
                    d = vdir / f"eval_{mode}"
                    w = d / "vcfeval" / "pon_calls"
                    if (d / "report.json").exists() and _rtg_complete(w):
                        u.rtg[mode] = w
                        u.report[mode] = json.loads((d / "report.json").read_text())
                if "bed" not in u.rtg:
                    continue
                rb = u.report["bed"]
                u.truth, u.bed = resolved(rb.get("truth")), resolved(rb.get("bed"))
                tp = vdir / "eval_bed" / f"truth.{kind}.pon.vcf.gz"
                u.truth_pon = tp if tp.exists() else None
                if u.variant != "standard":
                    std_dir = "vcf_chr1" if vdir.name == "vcf_chr1_valbest" or u.variant == "norescue" else vdir.name
                    u.standard_uid = f"{base_run}/{std_dir}"
                units.append(u)
    if comparators:
        units += read_comparators(kinds)
    units.sort(key=Unit.sort_key)
    return units


def read_comparators(kinds=("SNV", "INDEL")):
    """Rows of comparators.tsv (label, set, kind, eval_bed, eval_nobed, truth_pon, note); '#' lines are comments.
    eval_bed / eval_nobed: rtg vcfeval output directories of the PASS-after-our-PoN calls (score_comparator.py:
    <out>/eval_bed/pass_pon, <out>/eval_nobed/pass_pon). Truth and BED: the set's labels manifest."""
    out = []
    if not COMPARATORS.exists():
        return out
    with open(COMPARATORS) as f:
        rows = [line for line in f if line.strip() and not line.startswith("#")]
    for row in csv.DictReader(rows, fieldnames=["label", "set", "kind", "eval_bed", "eval_nobed", "truth_pon", "note"],
                              delimiter="\t"):
        if row["label"] == "label" or row["kind"] not in kinds:
            continue
        u = Unit(uid=f"comparator:{row['label']}:{row['set']}:{row['kind']}", run=row["label"], set=row["set"],
                 kind=row["kind"], role="comparator", variant="standard", model=row["label"], note=row.get("note") or "")
        for mode in ("bed", "nobed"):
            p = row.get(f"eval_{mode}")
            if p and p != "-" and _rtg_complete(Path(p)):
                u.rtg[mode] = Path(p)
                u.report[mode] = {}
        if "bed" not in u.rtg:
            continue
        m = manifest(row["set"], row["kind"])
        if m:
            u.truth, u.bed = resolved(m["truth"]["somatic"]["vcf"]), resolved(m["truth"]["somatic"]["bed"])
        tp = row.get("truth_pon")
        u.truth_pon = Path(tp) if tp and tp != "-" and Path(tp).exists() else None
        out.append(u)
    return out


def read_records(vcf):
    """(chrom, pos, ref, alt, filter, info) of each record of a (gzipped) VCF; alt = the first ALT."""
    out = []
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t", 8)
            out.append((c[0], int(c[1]), c[3].upper(), c[4].split(",")[0].upper(), c[6], c[7]))
    return out


def fasta_path(units=None):
    """The labels' GRCh38 FASTA (the same file for every set)."""
    for u in units or ():
        m = manifest(u.set, u.kind)
        if m:
            return m["fasta"]
    return "/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"


def fmt(x, digits=3):
    return "-" if x is None or x != x else f"{x:.{digits}f}"
