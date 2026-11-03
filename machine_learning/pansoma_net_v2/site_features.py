"""Per-tensor site features from population databases: <DIR>/<SET>.<KIND>.npz, one file per tensor set and kind,
arrays aligned with the KindIndex order of that kind directory (manifest chromosome order, each chromosome in its
<chrom>_variant_summary.ndjson line order).

    cd machine_learning
    python -m pansoma_net_v2.site_features build --tensors <set>/tensors --fasta GRCh38.fa \\
        --gnomad af-only-gnomad.hg38.vcf.gz --colors CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz --out DIR \\
        [--kind INDEL] [--cache-dir R/index_cache] [--workers 8] [--bcftools bcftools]
    python -m pansoma_net_v2.site_features build --kind SNV --tensors <set>/tensors --fasta GRCh38.fa \\
        --gnomad ... --colors ... --dbsnp Homo_sapiens_assembly38.dbsnp138.vcf.gz --g1000 1000g_pon.hg38.vcf.gz --out DIR

File contents (`build`):
- 'chroms' (str, the manifest's chromosomes in order), 'chrom_counts' (int64, tensors per chromosome) and
  'candidates_sha1' (sha1 of the index's candidate ids joined by newlines, KindIndex.candidates() order): a file that
  disagrees with the index is refused (`load`), so a rebuilt set or another tensor order cannot silently pair features
  with the wrong tensors. The features do not depend on the labels: a relabel keeps the file valid.
- 'popaf_exact' float32: the largest population AF of the exact normalised (left-aligned) ALT allele over gnomAD
  af-only and CoLoRSdb; 0 when absent; -1 when the tensor has no GRCh38 allele (off-reference).
- 'popaf_pos' float32: the largest AF of any INDEL allele of either database at the tensor's left-aligned POS.
- 'pon_tag' int8: 1 when scripts/filter_panel_of_normals.py's rule tags the allele, 0 when not, -1 off-reference.
- 'hp_len' int16: filter_panel_of_normals.homopolymer_run of the allele (0 off-reference).
- 'popaf_gnomad', 'popaf_colors' (the two components of popaf_exact), 'rawaf_gnomad', 'rawaf_colors' and
  'pon_tag_gnomad', 'pon_tag_colors' (the AF and tag of each database's exact, unnormalised record, as the PoN
  script matches it), 'vcf_pos' / 'vcf_ref' / 'vcf_alt' (the left-aligned VCF allele linear_vcf writes; vcf_pos is
  checked against the index's GRCh38 position).
- 'meta' (JSON): set, kind, databases, FASTA, the PoN rule constants ('pon_rule') and the sha256 of
  scripts/filter_panel_of_normals.py ('pon_script_sha256'). `load` refuses a file whose PoN rule differs from the
  current script's and `pon_check` reports a file built with another version of the script; a file without 'pon_rule'
  (written before it was recorded) is checked by the script's sha256 only.
  Older files may also carry 'labels_sha' (the labels manifest they were built next to); it is not checked.

Per tensor with a GRCh38 allele (the summary's grch38, the projection linear_vcf uses):
- the VCF allele: linear_vcf.Genome.left_align + Genome.vcf_allele;
- hp_len = filter_panel_of_normals.homopolymer_run of that allele;
- raw AF per database: the largest AF of a record with that exact POS, REF and an ALT equal to the allele's ALT, which
  is what filter_panel_of_normals.record_matches tests (no normalisation of the database records);
  pon_tag_<db> = raw AF >= the INDEL floor * (1 - 1e-6) (0.01 when hp_len >= INDEL_HOMOPOLYMER, else 1e-4, the
  script's constants), pon_tag = either database;
- normalised AF per database (popaf_gnomad, popaf_colors, max = popaf_exact): records split per ALT, trimmed (suffix,
  then prefix down to an empty allele), left-aligned with Genome.left_align and padded with Genome.vcf_allele;
- popaf_pos: the largest normalised AF of any pure INDEL allele of either database at the allele's left-aligned POS.
Database records come from `bcftools query -r <chrom> -e 'TYPE="snp"'` (records with only SNV ALTs cannot hold an
INDEL ALT). This is the INDEL build (--kind INDEL, the default).

SNV sets (--kind SNV, with --dbsnp and --g1000): the same arrays and the same model input, with the SNV rule of
scripts/filter_panel_of_normals.py. The VCF allele is linear_vcf's (POS = pos0 + 1, no left-alignment); hp_len is 0. popaf_<db> = rawaf_<db> = the largest AF of a gnomAD / CoLoRSdb record with that POS,
REF and ALT (no normalisation, as record_matches); popaf_pos = the largest AF of any single-base ALT of either
database at the POS. pon_tag_gnomad / pon_tag_colors: that AF >= SNV_MIN_AF (1e-4); pon_tag_dbsnp: a dbSNP record of
the allele whose INFO SAO is not 2 (somatic); pon_tag_g1000: any 1000G PoN record of the allele; pon_tag = any of
the four, which is the PoN script's PanelOfNormals tag of the allele. Records are read with bcftools query -r <chrom>
-T <the chromosome's tensor positions> (only records at a tensor's POS can match).

The model input (`transform`, train --site-features) is a fixed map to [0, 1], not z-scored (a z-score of a rare
0/1 flag fitted on training data can blow up on another set):
    log_popaf_exact = 1 + log10(popaf_exact + 1e-6) / 6     (0: absent or off-reference, 0.33: AF 1e-4, 1: AF 1)
    log_popaf_pos   = 1 + log10(popaf_pos + 1e-6) / 6
    pon_tag         = 1 if tagged else 0
    hp_len          = min(hp_len, 20) / 20
    has_grch38      = 1 if the tensor has a GRCh38 allele else 0
`VERSION` names this map; a checkpoint records it with the directory, and predict refuses another version.

Off-reference check, both ways: a tensor the index marks off-reference must have pon_tag -1 and popaf_exact < 0, and a
tensor with a GRCh38 allele must have pon_tag 0 / 1 and popaf_exact >= 0.
`load` also returns the sha256 of the file's bytes as loaded ('sha256'): train records it per set in the checkpoint
(feature_files) and predict warns when a training set's file has changed since.
"""
import argparse
import hashlib
import importlib.util
import io
import json
import re
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from .vcfio import REPO

VERSION = "site-features-v1"
SITE_FEATURES = ("log_popaf_exact", "log_popaf_pos", "pon_tag", "hp_len", "has_grch38")
TRANSFORM = ("log_popaf_exact = 1 + log10(popaf_exact + 1e-6) / 6; log_popaf_pos = 1 + log10(popaf_pos + 1e-6) / 6 "
             "(AF < 0, i.e. off-reference, as 0); pon_tag = pon_tag == 1; hp_len = min(hp_len, 20) / 20; "
             "has_grch38 = pon_tag >= 0 and popaf_exact >= 0")
REQUIRED = ("chroms", "chrom_counts", "candidates_sha1", "popaf_exact", "popaf_pos", "pon_tag", "hp_len")
PER_TENSOR = ("popaf_exact", "popaf_pos", "pon_tag", "hp_len", "popaf_gnomad", "popaf_colors", "vcf_pos", "vcf_ref",
              "vcf_alt")
AF_EPS, HP_CAP = 1e-6, 20
VCF_POS_SLACK = (-2, 5000)  # pos0 - vcf_pos of a left-aligned INDEL: >= 0 up to the convention, < a long repeat
VCF_POS_AGREE = 0.99        # share of on-reference tensors whose vcf_pos must fall in VCF_POS_SLACK
PON_SCRIPT = REPO / "scripts" / "filter_panel_of_normals.py"
DBS = ("gnomad", "colors")


def set_name(index):
    """The tensor set's name (e.g. HCC1395T_ONT_d46): the directory above <set>/tensors/<KIND>, or the parent of the
    kind directory when it is not called tensors."""
    parent = Path(index.dir).parent
    return parent.parent.name if parent.name == "tensors" else parent.name


def feature_path(directory, index):
    return Path(directory) / f"{set_name(index)}.{index.kind}.npz"


def _text(value):
    v = np.asarray(value).reshape(-1)
    if len(v) != 1:
        raise ValueError(f"expected one string, got {len(v)}")
    v = v[0]
    return v.decode() if isinstance(v, bytes) else str(v)


def candidates_sha1(index):
    """sha1 of the index's candidate ids joined by newlines (the builder's candidates_sha1). Ids the index had not
    read yet are dropped again afterwards: millions of strings need not stay in memory for a whole run."""
    loaded = index._candidates is not None
    digest = hashlib.sha1("\n".join(index.candidates()).encode()).hexdigest()
    if not loaded:
        index._candidates = None
    return digest


def pon_module():
    """scripts/filter_panel_of_normals.py of this repository, imported from its path."""
    spec = importlib.util.spec_from_file_location("filter_panel_of_normals", PON_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def pon_rule(module=None, kind="INDEL"):
    """The constants of the PoN script's rule for kind (INDEL or SNV), as the builder records them."""
    m = module or pon_module()
    if kind == "SNV":
        return dict(snv_min_af=m.SNV_MIN_AF, snv_rules=[dict(r) for r in m.PON_RULES["SNV"]])
    return dict(indel_min_af=m.INDEL_MIN_AF, indel_min_af_homopolymer=m.INDEL_MIN_AF_HOMOPOLYMER,
                indel_homopolymer=m.INDEL_HOMOPOLYMER, indel_rules=list(m.PON_RULES["INDEL"]))


def file_meta(data):
    """The file's meta JSON ({} when it has none)."""
    return json.loads(_text(data["meta"])) if "meta" in data else {}


def pon_check(data):
    """Warnings for a feature file built with another version of scripts/filter_panel_of_normals.py than the current
    one (its rule constants are checked by `load`)."""
    meta, out = file_meta(data), []
    built = meta.get("pon_script_sha256")
    current = hashlib.sha256(PON_SCRIPT.read_bytes()).hexdigest()
    if built and built != current:
        out.append(f"{data.get('path')}: built with filter_panel_of_normals.py sha256 {built[:12]}, the current script "
                   f"is {current[:12]} ({'same INDEL rule constants' if meta.get('pon_rule') else 'rule not recorded'})")
    return out


def load(directory, index):
    """The raw per-tensor arrays of <directory>/<SET>.<KIND>.npz for `index` (dict of numpy arrays, length
    len(index)), after checking them against the index; ValueError when the file is missing or disagrees.
    Also 'path' (the file) and 'sha256' (of the bytes loaded)."""
    path = feature_path(directory, index)
    if not path.exists():
        raise ValueError(f"{index.dir}: no site features {path}")
    raw = path.read_bytes()  # hashed and loaded from the same bytes: the hash is that of the arrays used
    try:
        with np.load(io.BytesIO(raw), allow_pickle=False) as z:
            data = {k: z[k] for k in z.files}
    except ValueError as e:  # object arrays: save strings as numpy unicode arrays
        raise ValueError(f"{path}: {e}") from None
    missing = [k for k in REQUIRED if k not in data]
    if missing:
        raise ValueError(f"{path}: missing {', '.join(missing)}")
    chroms = [c.decode() if isinstance(c, bytes) else str(c) for c in np.asarray(data["chroms"]).reshape(-1)]
    if chroms != list(index.meta["chroms"]):
        raise ValueError(f"{path}: chromosomes {chroms} are not the index's {index.meta['chroms']}")
    counts = np.bincount(index.arrays["chrom"].astype(np.int64), minlength=len(chroms))
    if not np.array_equal(np.asarray(data["chrom_counts"], np.int64).reshape(-1), counts):
        raise ValueError(f"{path}: chrom_counts {np.asarray(data['chrom_counts']).tolist()} are not the index's "
                         f"{counts.tolist()}")
    n = len(index)
    for k in PER_TENSOR:
        if k in data and len(data[k]) != n:
            raise ValueError(f"{path}: {k} has {len(data[k])} values, the index {n} tensors")
    for k in ("popaf_exact", "popaf_pos", "popaf_gnomad", "popaf_colors"):
        if k in data and (not np.isfinite(data[k]).all() or (data[k] > 1.0 + 1e-6).any()):
            raise ValueError(f"{path}: {k} has values that are not finite or above 1")
    want = candidates_sha1(index)
    if _text(data["candidates_sha1"]) != want:
        raise ValueError(f"{path}: candidates_sha1 {_text(data['candidates_sha1'])[:12]} is not that of the "
                         f"index's candidate ids ({want[:12]}): stale features or another tensor order")
    rule = file_meta(data).get("pon_rule")
    if rule is not None and rule != pon_rule(kind=index.kind):
        raise ValueError(f"{path}: PoN tags of the rule {rule}, scripts/filter_panel_of_normals.py now has "
                         f"{pon_rule(kind=index.kind)}: stale features")
    off = index.arrays["off_reference"]
    no_allele = (data["pon_tag"] == -1) | (data["popaf_exact"] < 0)
    has_allele = (data["pon_tag"] != -1) | (data["popaf_exact"] >= 0)
    bad = int((off & has_allele).sum())
    if bad:  # the index says no GRCh38 allele, the file gives one: not this index's order
        raise ValueError(f"{path}: {bad:,} off-reference tensors have a pon_tag or popaf_exact other than -1")
    bad = int((~off & no_allele).sum())
    if bad:  # the index has a GRCh38 allele, the file none
        raise ValueError(f"{path}: {bad:,} tensors with a GRCh38 allele (not off-reference) have pon_tag or "
                         f"popaf_exact -1")
    if "vcf_pos" in data:
        on = ~off & (data["pon_tag"] >= 0) & (np.asarray(data["vcf_pos"]) > 0)
        if on.any():
            d = index.arrays["pos0"][on] - np.asarray(data["vcf_pos"], np.int64)[on]
            agree = float(((d >= VCF_POS_SLACK[0]) & (d <= VCF_POS_SLACK[1])).mean())
            if agree < VCF_POS_AGREE:
                raise ValueError(f"{path}: vcf_pos agrees with the index's GRCh38 position on {agree:.1%} of the "
                                 f"on-reference tensors only (pos0 - vcf_pos in {list(VCF_POS_SLACK)})")
    data["path"] = str(path)
    data["sha256"] = hashlib.sha256(raw).hexdigest()
    return data


def transform(data):
    """(n, len(SITE_FEATURES)) float32 model inputs of the raw arrays (module docstring)."""
    def log_af(v):
        v = np.clip(np.asarray(v, np.float64), 0.0, 1.0)
        return 1.0 + np.log10(v + AF_EPS) / 6.0
    pon = np.asarray(data["pon_tag"])
    has = (pon >= 0) & (np.asarray(data["popaf_exact"]) >= 0)
    out = np.stack([log_af(data["popaf_exact"]), log_af(data["popaf_pos"]), (pon == 1).astype(np.float64),
                    np.minimum(np.maximum(np.asarray(data["hp_len"], np.float64), 0), HP_CAP) / HP_CAP,
                    has.astype(np.float64)], 1).astype(np.float32)
    return np.clip(out, 0.0, 1.0)


def describe(index, data):
    """One log line of a set's features."""
    pon, off = np.asarray(data["pon_tag"]), index.arrays["off_reference"]
    return (f"{index.dir}: site features {data['path']}: {len(pon):,} tensors, {int((pon == 1).sum()):,} PoN-tagged, "
            f"{int((pon == -1).sum()):,} without a GRCh38 allele ({int(off.sum()):,} off-reference), "
            f"{int((np.asarray(data['popaf_exact']) > 0).sum()):,} with a population AF")


def attach(index, directory, model_input=True):
    """Load the set's features for `index`: index.arrays['pon_tag'] (int8, the raw tag) always and, with model_input,
    index.arrays['site_features'] (the transform), which TensorDataset appends to the scalars. Returns the raw arrays."""
    data = load(directory, index)
    index.arrays["pon_tag"] = np.asarray(data["pon_tag"], np.int8)
    if model_input:
        index.arrays["site_features"] = transform(data)
    return data


def key(index):
    """<SET>.<KIND> of an index: the name of its feature file without .npz, the key of the checkpoint's records."""
    return f"{set_name(index)}.{index.kind}"


def record(directory, loaded):
    """The checkpoint's record of --site-features: the directory, version, names, transform and, per set
    (loaded: {<SET>.<KIND>: load's arrays}), the file's candidates_sha1, PoN rule and PoN script sha256."""
    sets = {}
    for k, data in sorted(loaded.items()):
        meta = file_meta(data)
        sets[k] = dict(candidates_sha1=_text(data["candidates_sha1"]), pon_rule=meta.get("pon_rule"),
                       pon_script_sha256=meta.get("pon_script_sha256"))
    return dict(directory=str(Path(directory).resolve()), version=VERSION, names=list(SITE_FEATURES),
                transform=TRANSFORM, sets=sets)


def compare_files(loaded, recorded, what):
    """Warnings for the sets whose feature file differs from the one a checkpoint was trained with.
    loaded: {<SET>.<KIND>: (path, sha256)} of the files read now; recorded: the checkpoint's {<SET>.<KIND>: sha256}
    (feature_files[what]) or None. Returns {key: {path, sha256, checkpoint_sha256, changed}} and the warning lines;
    a set the checkpoint did not train on (or a checkpoint without the record) has changed None."""
    out, warnings = {}, []
    for k, (path, sha) in sorted(loaded.items()):
        trained = (recorded or {}).get(k)
        changed = None if trained is None else trained != sha
        out[k] = dict(path=str(path), sha256=sha, checkpoint_sha256=trained, changed=changed)
        if changed:
            warnings.append(f"warning: {what} file {path} of training set {k} has changed since the checkpoint was "
                            f"trained (sha256 {sha[:12]}, trained with {trained[:12]})")
    return out, warnings


# --- builder -----------------------------------------------------------------------------------------------------

_CANDIDATE = re.compile(rb'"candidate_id": "([^"]+)"')
_GRCH38 = re.compile(rb'"grch38": (\{[^}]*\}|null)')  # data.py's pattern


class _MemGenome:
    """linear_vcf.Genome over one upper-case chromosome held in memory (same get() semantics)."""

    def __init__(self, chrom, seq):
        from .linear_vcf import Genome
        self.chrom, self.seq, self.lengths = chrom, seq, {chrom: len(seq)}
        self.left_align = Genome.left_align.__get__(self)
        self.vcf_allele = Genome.vcf_allele.__get__(self)

    def get(self, chrom, start, end):
        return self.seq[max(0, start):min(end, len(self.seq))]


class _MemFasta:
    """pysam.FastaFile.fetch stand-in for filter_panel_of_normals.homopolymer_run."""

    def __init__(self, chrom, seq):
        self.chrom, self.seq = chrom, seq

    def fetch(self, contig, start, end):
        if contig != self.chrom:
            raise ValueError(f"{contig} is not {self.chrom}")
        return self.seq[max(0, start):end]


def _read_tensors(summary, chrom, genome, fasta, pon):
    """Per summary line: (candidate id, allele key or None, vcf_pos, vcf_ref, vcf_alt, hp_len); off-reference
    tensors (grch38 null) get None, -1, '', '', 0."""
    rows = []
    with open(summary, "rb") as f:
        for line in f:
            cid = _CANDIDATE.search(line).group(1).decode()
            g = json.loads(_GRCH38.search(line).group(1))
            if g is None:
                rows.append((cid, None, -1, "", "", 0))
                continue
            if g["chrom"] != chrom:
                raise ValueError(f"{summary} {cid}: grch38 chromosome {g['chrom']}")
            ref, alt = (g["ref"] or "").upper(), (g["alt"] or "").upper()
            if ref and not alt:
                kind = "DEL"
            elif alt and not ref:
                kind = "INS"
            else:
                raise ValueError(f"{summary} {cid}: not a pure INS / DEL: {g}")
            pos0 = int(g["pos0"])
            if genome.get(chrom, pos0, pos0 + len(ref)) != ref:  # linear_vcf stops on this
                raise ValueError(f"{summary} {cid}: REF {ref} differs from the FASTA")
            pos0, ref, alt, _ = genome.left_align(chrom, pos0, ref, alt, kind)
            vpos, vref, valt = genome.vcf_allele(chrom, pos0, ref, alt, kind)
            rows.append((cid, (vpos, vref, valt), vpos, vref, valt, pon.homopolymer_run(fasta, chrom, vpos, vref, valt)))
    return rows


def _normalise(genome, chrom, pos, ref, alt):
    """(POS, REF, ALT) of a database INDEL allele after trimming, left-alignment and padding; None if not a pure
    INDEL or its REF differs from the FASTA."""
    seq = genome.seq
    if seq[pos - 1:pos - 1 + len(ref)] != ref:
        return None
    # shared bases off both ends, down to an empty allele: also catches right-anchored INDELs of multi-allelic records
    r, a, p0 = ref, alt, pos - 1  # p0: 0-based start of r
    while r and a and r[-1] == a[-1]:
        r, a = r[:-1], a[:-1]
    while r and a and r[0] == a[0]:
        r, a, p0 = r[1:], a[1:], p0 + 1
    if r and a:
        return None
    kind, evref, evalt = ("DEL", r, "") if r else ("INS", "", a)
    pos0, evref, evalt, _ = genome.left_align(chrom, p0, evref, evalt, kind)
    return genome.vcf_allele(chrom, pos0, evref, evalt, kind)


def _stream_db(path, chrom, genome, alleles, positions, raw, norm, posmax, bcftools):
    """One pass over a database's non-SNV-only records of chrom: raw and normalised exact-allele AFs, AF per POS."""
    cmd = [bcftools, "query", "-r", chrom, "-e", 'TYPE="snp"', "-f", r"%POS\t%REF\t%ALT\t%AF\n", str(path)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in proc.stdout:
        pos, ref, alts, afs = line.rstrip("\n").split("\t")
        pos, ref, alts, afs = int(pos), ref.upper(), alts.split(","), afs.split(",")
        if len(afs) != len(alts):
            if afs != ["."]:
                raise ValueError(f"{path} {chrom}:{pos}: {len(alts)} ALTs, AF {afs}")
            afs = ["."] * len(alts)
        for alt, af in zip(alts, afs):
            if not alt or alt.startswith("<") or alt in ("*", ".") or len(alt) == len(ref):
                continue
            alt = alt.upper()
            af = None if af == "." else float(np.float32(af))  # pysam reads Type=Float as float32
            idx = alleles.get((pos, ref, alt))
            if idx is not None and af is not None and af > raw[idx]:
                raw[idx] = af
            n = _normalise(genome, chrom, pos, ref, alt)
            if n is None or af is None:
                continue
            idx = alleles.get(n)
            if idx is not None and af > norm[idx]:
                norm[idx] = af
            if n[0] in positions and af > posmax.get(n[0], 0.0):
                posmax[n[0]] = af
    if proc.wait() != 0:
        raise RuntimeError(f"{' '.join(cmd)} failed")


def build_chrom(kind_dir, chrom, summary, fasta_path, dbs, bcftools):
    """The per-tensor arrays of one chromosome's summary lines (dict of numpy arrays, plus the candidate ids)."""
    import pysam
    pon = pon_module()
    with pysam.FastaFile(str(fasta_path)) as f:
        seq = f.fetch(chrom).upper()
    genome, fasta = _MemGenome(chrom, seq), _MemFasta(chrom, seq)
    rows = _read_tensors(Path(kind_dir) / summary, chrom, genome, fasta, pon)
    alleles = {}
    for r in rows:
        if r[1] is not None:
            alleles.setdefault(r[1], len(alleles))
    positions = {k[0] for k in alleles}
    raw = {db: np.zeros(len(alleles)) for db in dbs}
    norm = {db: np.zeros(len(alleles)) for db in dbs}
    posmax = {}
    for db, path in dbs.items():
        _stream_db(path, chrom, genome, alleles, positions, raw[db], norm[db], posmax, bcftools)
    idx = np.array([alleles[r[1]] if r[1] is not None else -1 for r in rows], np.int64)
    on = idx >= 0
    hp = np.array([r[5] for r in rows], np.int16)
    out = dict(candidates=[r[0] for r in rows], vcf_pos=np.array([r[2] for r in rows], np.int64),
               vcf_ref=np.array([r[3] for r in rows], dtype="S"), vcf_alt=np.array([r[4] for r in rows], dtype="S"),
               hp_len=hp)
    floor = np.where(hp >= pon.INDEL_HOMOPOLYMER, pon.INDEL_MIN_AF_HOMOPOLYMER, pon.INDEL_MIN_AF)
    for db in dbs:
        a = np.full(len(rows), -1.0)
        a[on] = norm[db][idx[on]]
        out[f"popaf_{db}"] = a.astype(np.float32)
        r = np.full(len(rows), -1.0)
        r[on] = raw[db][idx[on]]
        out[f"rawaf_{db}"] = r.astype(np.float32)
        tag = np.full(len(rows), -1, np.int8)
        tag[on] = (r[on] >= floor[on] * (1 - 1e-6)).astype(np.int8)  # record_matches' comparison
        out[f"pon_tag_{db}"] = tag
    pp = np.full(len(rows), -1.0)
    pp[on] = [posmax.get(int(p), 0.0) for p in out["vcf_pos"][on]]
    out["popaf_pos"] = pp.astype(np.float32)
    return out


SNV_DBS = ("gnomad", "colors", "dbsnp", "g1000")


def _read_snv_tensors(summary, chrom, genome):
    """Per summary line of an SNV kind directory: (candidate id, (POS, REF, ALT) or None); off-reference tensors
    (grch38 null) get None."""
    rows = []
    with open(summary, "rb") as f:
        for line in f:
            cid = _CANDIDATE.search(line).group(1).decode()
            g = json.loads(_GRCH38.search(line).group(1))
            if g is None:
                rows.append((cid, None))
                continue
            if g["chrom"] != chrom:
                raise ValueError(f"{summary} {cid}: grch38 chromosome {g['chrom']}")
            ref, alt = (g["ref"] or "").upper(), (g["alt"] or "").upper()
            if len(ref) != 1 or len(alt) != 1 or ref == alt:
                raise ValueError(f"{summary} {cid}: not a single-base substitution: {g}")
            pos0 = int(g["pos0"])
            if genome.get(chrom, pos0, pos0 + 1) != ref:  # linear_vcf stops on this
                raise ValueError(f"{summary} {cid}: REF {ref} differs from the FASTA")
            rows.append((cid, genome.vcf_allele(chrom, pos0, ref, alt, "SNP")))
    return rows


def _stream_snv_db(db, path, chrom, targets, alleles, raw, tag, posmax, bcftools, min_af):
    """One pass over a database's records of chrom at the tensors' positions (bcftools -T targets): per allele the
    largest AF (gnomad / colors) and the PoN tag of filter_panel_of_normals.record_matches under the SNV rule;
    per POS the largest AF of any single-base ALT (gnomad / colors)."""
    field = {"gnomad": r"\t%AF", "colors": r"\t%AF", "dbsnp": r"\t%INFO/SAO", "g1000": ""}[db]
    cmd = [bcftools, "query", "-r", chrom, "-T", str(targets), "-f", rf"%POS\t%REF\t%ALT{field}\n", str(path)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in proc.stdout:
        x = line.rstrip("\n").split("\t")
        pos, ref, alts = int(x[0]), x[1].upper(), x[2].upper().split(",")
        if len(ref) != 1:  # record_matches needs REF equal to the allele's one base
            continue
        if db in ("gnomad", "colors"):
            afs = x[3].split(",")
            if len(afs) != len(alts):
                if afs != ["."]:
                    raise ValueError(f"{path} {chrom}:{pos}: {len(alts)} ALTs, AF {afs}")
                afs = ["."] * len(alts)
        for k, alt in enumerate(alts):
            if len(alt) != 1 or alt == ref or alt in ("*", "."):
                continue
            idx = alleles.get((pos, ref, alt))
            if db in ("gnomad", "colors"):
                af = None if afs[k] == "." else float(np.float32(afs[k]))  # pysam reads Type=Float as float32
                if af is None:
                    continue
                if af > posmax.get(pos, 0.0):
                    posmax[pos] = af
                if idx is not None:
                    raw[idx] = max(raw[idx], af)
                    if af >= min_af * (1 - 1e-6):  # record_matches' comparison
                        tag[idx] = 1
            elif idx is not None:
                if db == "dbsnp" and x[3] == "2":  # non_somatic: SAO == 2 does not match
                    continue
                tag[idx] = 1
    if proc.wait() != 0:
        raise RuntimeError(f"{' '.join(cmd)} failed")


def build_chrom_snv(kind_dir, chrom, summary, fasta_path, dbs, bcftools):
    """The per-tensor arrays of one chromosome of an SNV kind directory (dict of numpy arrays, plus the ids)."""
    import pysam
    import tempfile
    pon = pon_module()
    with pysam.FastaFile(str(fasta_path)) as f:
        seq = f.fetch(chrom).upper()
    genome = _MemGenome(chrom, seq)
    rows = _read_snv_tensors(Path(kind_dir) / summary, chrom, genome)
    alleles = {}
    for r in rows:
        if r[1] is not None:
            alleles.setdefault(r[1], len(alleles))
    raw = {db: np.zeros(len(alleles)) for db in ("gnomad", "colors")}
    tag = {db: np.zeros(len(alleles), np.int8) for db in SNV_DBS}
    posmax = {}
    with tempfile.NamedTemporaryFile("w", suffix=f".{chrom}.targets.tsv", delete=False) as t:
        for p in sorted({k[0] for k in alleles}):
            t.write(f"{chrom}\t{p}\n")
        targets = t.name
    try:
        for db in SNV_DBS:
            _stream_snv_db(db, dbs[db], chrom, targets, alleles, raw.get(db), tag[db], posmax if db in raw else {},
                           bcftools, pon.SNV_MIN_AF)
    finally:
        Path(targets).unlink()
    idx = np.array([alleles[r[1]] if r[1] is not None else -1 for r in rows], np.int64)
    on = idx >= 0
    out = dict(candidates=[r[0] for r in rows],
               vcf_pos=np.array([r[1][0] if r[1] else -1 for r in rows], np.int64),
               vcf_ref=np.array([r[1][1] if r[1] else "" for r in rows], dtype="S"),
               vcf_alt=np.array([r[1][2] if r[1] else "" for r in rows], dtype="S"),
               hp_len=np.zeros(len(rows), np.int16))
    for db in SNV_DBS:
        if db in raw:
            a = np.full(len(rows), -1.0)
            a[on] = raw[db][idx[on]]
            out[f"popaf_{db}"] = a.astype(np.float32)
            out[f"rawaf_{db}"] = a.astype(np.float32)
        t = np.full(len(rows), -1, np.int8)
        t[on] = tag[db][idx[on]]
        out[f"pon_tag_{db}"] = t
    pp = np.full(len(rows), -1.0)
    pp[on] = [posmax.get(int(p), 0.0) for p in out["vcf_pos"][on]]
    out["popaf_pos"] = pp.astype(np.float32)
    return out


def build(tensors, kind, fasta, dbs, out_dir, cache_dir=None, workers=1, bcftools="bcftools", log=print):
    """Build <out_dir>/<SET>.<KIND>.npz for the kind directory <tensors>/<KIND> (module docstring); returns its
    path. The KindIndex (from cache_dir, built there if missing) fixes the order and is checked against the summaries."""
    from .data import KindIndex
    kind_dir = Path(tensors) / kind
    manifest = json.loads((kind_dir / "manifest.json").read_text())
    chroms = list(manifest["chromosomes"])
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=max(1, workers)) as pool:
        fn = build_chrom_snv if kind == "SNV" else build_chrom
        futures = {c: pool.submit(fn, str(kind_dir), c, manifest["chromosomes"][c]["summary"], str(fasta),
                                  {k: str(v) for k, v in dbs.items()}, bcftools) for c in chroms}
        parts = []
        for c in chroms:
            parts.append(futures[c].result())
            log(f"{kind_dir} {c}: {len(parts[-1]['candidates']):,} tensors ({time.time() - t0:.0f} s)")
    counts = np.array([len(p["candidates"]) for p in parts], np.int64)
    for c, n in zip(chroms, counts):
        if n != manifest["chromosomes"][c]["tensors"]:
            raise ValueError(f"{kind_dir} {c}: {n} summary lines, the manifest lists "
                             f"{manifest['chromosomes'][c]['tensors']}")
    cands = [x for p in parts for x in p["candidates"]]
    cat = {k: np.concatenate([p[k] for p in parts]) for k in parts[0] if k != "candidates"}
    index = KindIndex(kind_dir, cache_dir=cache_dir)
    if len(index) != len(cands) or list(index.meta["chroms"]) != chroms or list(index.candidates()) != cands:
        raise ValueError(f"{kind_dir}: the summaries' candidate ids are not the KindIndex order")
    if not np.array_equal(cat["vcf_pos"] < 0, index.arrays["off_reference"]):
        raise ValueError(f"{kind_dir}: off-reference tensors differ from the KindIndex")
    pon = pon_module()
    af_dbs = ("gnomad", "colors") if kind == "SNV" else tuple(dbs)
    popaf_exact = np.max([cat[f"popaf_{db}"] for db in af_dbs], 0)
    pon_tag = np.max([cat[f"pon_tag_{db}"] for db in dbs], 0)
    name = Path(tensors).parent.name if Path(tensors).name == "tensors" else Path(tensors).name
    arrays = dict(chroms=np.array(chroms, dtype="U"), chrom_counts=counts,
                  candidates_sha1=np.array(hashlib.sha1("\n".join(cands).encode()).hexdigest()),
                  popaf_exact=popaf_exact.astype(np.float32), popaf_pos=cat["popaf_pos"],
                  pon_tag=pon_tag.astype(np.int8), hp_len=cat["hp_len"], vcf_pos=cat["vcf_pos"],
                  vcf_ref=cat["vcf_ref"], vcf_alt=cat["vcf_alt"],
                  meta=np.array(json.dumps(dict(
                      set=name, kind=kind, kind_dir=str(kind_dir.resolve()), tensors=int(len(index)),
                      pon_rule=pon_rule(pon, kind), pon_script_sha256=hashlib.sha256(PON_SCRIPT.read_bytes()).hexdigest(),
                      dbs={k: str(v) for k, v in dbs.items()}, fasta=str(fasta),
                      created=time.strftime("%Y-%m-%dT%H:%M:%S")))),
                  **{f"{a}_{db}": cat[f"{a}_{db}"] for db in dbs for a in ("popaf", "rawaf", "pon_tag")
                     if f"{a}_{db}" in cat})
    out = Path(out_dir) / f"{name}.{kind}.npz"
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name[:-len(".npz")] + ".tmp.npz")
    np.savez_compressed(tmp, **arrays)
    tmp.replace(out)
    on = cat["vcf_pos"] >= 0
    log(f"{out}: {len(index):,} tensors, {int(on.sum()):,} with a GRCh38 allele, {int((popaf_exact > 0).sum()):,} "
        f"with a population AF, {int((pon_tag == 1).sum()):,} PoN-tagged ({time.time() - t0:.0f} s)")
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="command", required=True)
    b = sub.add_parser("build", help="build <out>/<SET>.<KIND>.npz of one tensor set")
    b.add_argument("--tensors", required=True, help="the set's tensors directory (holding <KIND>/)")
    b.add_argument("--kind", default="INDEL", choices=["INDEL", "SNV"])
    b.add_argument("--fasta", required=True, help="the GRCh38 FASTA of the set's labels / linear_vcf")
    b.add_argument("--gnomad", required=True, help="gnomAD af-only VCF (the PoN script's PoN1)")
    b.add_argument("--colors", required=True, help="CoLoRSdb VCF (the PoN script's PoN4)")
    b.add_argument("--dbsnp", help="dbSNP VCF (the PoN script's PoN2; --kind SNV)")
    b.add_argument("--g1000", help="1000G PoN VCF (the PoN script's PoN3; --kind SNV)")
    b.add_argument("--out", required=True, help="output directory")
    b.add_argument("--cache-dir", help="KindIndex cache directory (as train / predict --cache-dir)")
    b.add_argument("--workers", type=int, default=1, help="chromosomes built in parallel")
    b.add_argument("--bcftools", default="bcftools")
    args = p.parse_args(argv)
    dbs = {"gnomad": args.gnomad, "colors": args.colors}
    if args.kind == "SNV":
        if not (args.dbsnp and args.g1000):
            p.error("--kind SNV needs --dbsnp and --g1000 (the SNV PoN rule uses all four databases)")
        dbs.update(dbsnp=args.dbsnp, g1000=args.g1000)
    build(args.tensors, args.kind, args.fasta, dbs, args.out,
          cache_dir=args.cache_dir, workers=args.workers, bcftools=args.bcftools,
          log=lambda m: print(time.strftime("%H:%M:%S"), m, flush=True))


if __name__ == "__main__":
    main()
