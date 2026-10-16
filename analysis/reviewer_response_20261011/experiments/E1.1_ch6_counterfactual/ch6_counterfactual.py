"""E1.1: does a trained PansomaNetV2 SNV model use channel 6 (ch6, the haplotype path count)?

Inference-time counterfactuals on the chr1 test tensors of one set, no retraining. The checkpoint, the tensors and
the label files are only read. ch6 is edited in memory, on the cells the encoder treats as valid (ch6 > 0); cells
without a read stay 0, every other channel is untouched.

Modes (codes are data.PATH_COUNT_ALL's int8-count-haplotypes100: H haplotypes -> 100):
  base           the tensors as they are
  all_H          every valid ch6 cell set to 100 ("every haplotype of the graph runs through this node")
  train_mean     every valid ch6 cell set to round(mean of the training tensors' valid ch6 cells): z ~ 0, the plane
                 carries no information (mean imputation)
  train_median   the same with the training median (skipped when it is 100: then it equals all_H up to the repeat
                 codes > 100)
  shuffle        ch6 of another test tensor of the set (a fixed random permutation): the donor's cell where the
                 donor has a read there, else the donor's column median (else the donor's median); the target's
                 coverage pattern is kept
  site0_low      only the site column (50): valid cells set to LOW
  site0_high     only the site column: cells below 100 raised to 100
  site10_low     columns 40-60 (site +- 10): valid cells set to LOW
  site10_high    columns 40-60: cells below 100 raised to 100 (predict's off-reference rescue transform, applied to
                 every tensor)
LOW = the median site code (median ch6 of the site column over the rows that have one) of the set's off-reference
test tensors (off_reference_no_truth_match: a node GRCh38 does not traverse); with fewer than 50 of them (the
GRCh38-only graph) the 5th percentile of the site codes of all test tensors.

Test tensors are predict's: eval_label >= 0 (labels 0/1/2 in the confident region, off_reference_no_truth_match as
0), with the checkpoint's --ignore-reasons left out (SNV: the partial somatic matches). Thresholds: the checkpoint's
validation recall-0.9 threshold (val.somatic_at_recall["0.9"], the PASS threshold of graph_vcf) and argmax.

Per mode: somatic AP (per tensor), somatic / germline / non tensors with p_somatic >= threshold, truth-level TP / FP /
P / R / F1 (metrics.truth_report: each chr1 somatic truth allele once), and the same after a per-tensor PoN (the
on-reference tensors' GRCh38 alleles tagged with scripts/filter_panel_of_normals.py's SNV rule; off-reference tensors
have no GRCh38 position and are left out there, as rtg leaves out unplaced calls). With --assembly, the scores of the
assembly off-reference somatic SNV keys present in the set (all autosomes; none of them is a training label: they
are -1 or an ignored partial match).

usage: python ch6_counterfactual.py --run RUN --set SET --out DIR [--threads 16] [--limit N]
"""
import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import torch

# A frozen copy of the repository code at 41f7830 (git archive; the code the RIS checkpoints were tested with, the
# same pansoma_net_v2 behaviour as 7fc7396): the shared checkout changes under running jobs (cbc94a7, 2026-10-13
# 11:20, changed TensorDataset and the encoder), and DataLoader workers import the package anew.
REPO = Path("/scratch/jshen/data/reviewer_experiments_20261013/E1.1_ch6_counterfactual/code_41f7830")
sys.path.insert(0, str(REPO / "machine_learning"))
sys.path.insert(0, str(REPO / "scripts"))
from pansoma_net_v2 import metrics  # noqa: E402
from pansoma_net_v2.data import (INDEX_VERSION, KindIndex, TensorDataset, EpochSampler, _labels_meta,  # noqa: E402
                                 somatic_truth)
from pansoma_net_v2.model import PansomaNetV2  # noqa: E402
from pansoma_net_v2.train import make_loader  # noqa: E402
import filter_panel_of_normals as pon_rules  # noqa: E402

RUNS = Path("/scratch/jshen/data/pansoma_net_v2_runs/ris_7fc7396_20261012")
TENSORS = Path("/scratch/jshen/data/pansoma_v2_tensors")
SHARED_CACHE = Path("/scratch/jshen/data/pansoma_net_v2_runs/index_cache")
ASSEMBLY = Path("/scratch/jshen/data/pansoma_net_v2_runs/assembly_offref_20260929/offref_somatic.tsv")
PON_DIR = Path("/scratch/jshen/data/Pansoma/panel_of_normal_VCFs")
PONS = [PON_DIR / f for f in ("af-only-gnomad.hg38.vcf.gz", "Homo_sapiens_assembly38.dbsnp138.vcf.gz",
                              "1000g_pon.hg38.vcf.gz", "CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz")]
PC, SITE, ALL = 6, 50, 100
MODES = ("base", "all_H", "train_mean", "train_median", "shuffle", "site0_low", "site0_high", "site10_low",
         "site10_high")
GRCH38_RE = re.compile(r'"grch38": (\{[^}]*\}|null)')


DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def open_index(kind_dir, own_cache):
    """KindIndex from the shared cache when its entry is current (read only), else built into own_cache."""
    d = Path(kind_dir).resolve()
    name = f"{d.parent.parent.name}.{d.parent.name}.{d.name}.{hashlib.sha1(str(d).encode()).hexdigest()[:12]}"
    npz, txt = SHARED_CACHE / f"{name}.npz", SHARED_CACHE / f"{name}.candidates.txt"
    if npz.exists() and txt.exists():
        with np.load(npz) as z:
            meta = json.loads(str(z["meta"]))
        if meta.get("index_version") == INDEX_VERSION and meta["labels"] == _labels_meta(d) and meta["directory"] == str(d):
            log(f"index: shared cache {npz.name}")
            return KindIndex(d, SHARED_CACHE)
    log(f"index: building into {own_cache}")
    return KindIndex(d, own_cache)


class Counterfactual(TensorDataset):
    """TensorDataset whose ch6 is edited by `mode` (see the module docstring)."""

    def __init__(self, parts, mode, params, donors=None):
        super().__init__(parts, labels="eval_label")
        self.mode, self.params, self.donors = mode, params, donors

    def __getitem__(self, i):
        x, blocks, scalars, label = super().__getitem__(i)
        pc = x[PC]
        m, p = self.mode, self.params
        valid = pc > 0
        if m == "all_H":
            pc[valid] = ALL
        elif m in ("train_mean", "train_median"):
            pc[valid] = p[m]
        elif m == "shuffle":
            d = super().__getitem__(int(self.donors[i]))[0][PC].to(torch.float32)
            has = d > 0
            dn = d.masked_fill(~has, float("nan"))
            col = torch.nanmedian(dn, dim=0).values
            overall = torch.nanmedian(dn)
            if torch.isnan(overall):
                overall = torch.tensor(float(p["train_median"]))
            col = torch.where(torch.isnan(col), overall, col)
            fill = torch.where(has, d, col[None, :].expand_as(d))
            pc[valid] = fill[valid].round().to(pc.dtype)
        elif m.startswith("site"):
            w = int(m[4:].split("_")[0])
            win = pc[:, SITE - w:SITE + w + 1]
            if m.endswith("_low"):
                win[win > 0] = p["low"]
            else:
                win[(win > 0) & (win < ALL)] = ALL
        elif m != "base":
            raise ValueError(m)
        return x, blocks, scalars, label


def site_profile(x):
    """Per tensor (numpy): site code (median ch6 of column 50 over rows with a value; 0 if none), the fraction of
    valid ch6 cells below 100, and the number of distinct ch6 codes in columns 40-60."""
    pc = x[:, PC].to(torch.float32)
    col = pc[:, :, SITE]
    site = torch.nanmedian(col.masked_fill(col <= 0, float("nan")), dim=1).values
    site = torch.nan_to_num(site, nan=0.0)
    valid = pc > 0
    below = ((pc < ALL) & valid).sum((1, 2)).float() / valid.sum((1, 2)).clamp(min=1).float()
    win = pc[:, :, SITE - 10:SITE + 11]
    distinct = torch.tensor([len(torch.unique(w[w > 0])) for w in win])
    return site.numpy(), below.numpy(), distinct.numpy()


def run_mode(model, ds, args, profile=False):
    loader = make_loader(ds, EpochSampler(len(ds), shuffle=False), args.workers, args.batch_size, persistent=False)
    probs, prof = [], []
    t0, done, nxt = time.time(), 0, 0
    amp = DEVICE.type == "cuda"  # bf16 autocast on a GPU, as predict does; fp32 on CPU
    with torch.no_grad():
        for batch in loader:  # (x, blocks, scalars, label[, part]); part: the dataset row of the BQ / MAPQ z-score
            x, blocks, scalars = batch[0], batch[1], batch[2]
            part = batch[4].to(DEVICE) if len(batch) > 4 else None
            if profile:
                prof.append(site_profile(x))
            with torch.autocast(DEVICE.type, dtype=torch.bfloat16, enabled=amp):
                inputs = [x.to(DEVICE, non_blocking=True), blocks.to(DEVICE, non_blocking=True),
                          scalars.to(DEVICE, non_blocking=True)] + ([part] if part is not None else [])
                logits = model(*inputs)
            probs.append(torch.softmax(logits.float(), 1).cpu().numpy())
            done += len(x)
            if done >= nxt:
                rate = done / max(time.time() - t0, 1e-9)
                log(f"  {ds.mode}: {done:,}/{len(ds):,} tensors, {rate:.1f}/s")
                nxt += 20 * args.batch_size
    log(f"  {ds.mode}: {done:,} tensors in {time.time() - t0:.0f}s ({done / max(time.time() - t0, 1e-9):.1f}/s)")
    probs = np.concatenate(probs) if probs else np.zeros((0, 3), np.float32)
    if profile:
        prof = tuple(np.concatenate([p[k] for p in prof]) for k in range(3))
        return probs, prof
    return probs


def train_ch6_stats(index, ck, n, seed=0):
    """Histogram, mean and median of the valid ch6 cells of n random training tensors (the checkpoint's training
    chromosomes, label >= 0, its ignore reasons left out)."""
    ignore = set(ck.get("args", {}).get("ignore_reasons") or [])
    pos = index.select(ck["chroms"]["train"], labelled=True)
    bad = np.isin(index.arrays["reason"][pos], [k for k, r in enumerate(index.meta["reasons"]) if r in ignore])
    pos = pos[~bad]
    rng = np.random.default_rng(seed)
    pick = np.sort(rng.choice(pos, min(n, len(pos)), replace=False))
    ds = TensorDataset([(index, pick)])
    hist = np.zeros(128, np.int64)
    for i in range(len(ds)):
        x = ds[i][0]
        v = x[PC][(x[4] > 0) & (x[PC] > 0)].numpy().astype(np.int64)
        hist += np.bincount(v, minlength=128)[:128]
    codes = np.arange(128)
    total = hist.sum()
    mean = float((hist * codes).sum() / total)
    median = int(codes[np.searchsorted(np.cumsum(hist), (total + 1) / 2)])
    return dict(tensors=int(len(pick)), cells=int(total), mean=mean, median=median,
                frac_100=float(hist[100] / total), frac_below_100=float(hist[1:100].sum() / total),
                frac_above_100=float(hist[101:].sum() / total), hist={int(k): int(v) for k, v in enumerate(hist) if v})


def grch38_alleles(index, positions):
    """{position: (chrom, pos1, ref, alt)} of the on-reference tensors among positions (chr1 only), from the summary."""
    a, chroms = index.arrays, index.meta["chroms"]
    out = {}
    want = set(int(p) for p in positions)
    manifest = json.loads((index.dir / "manifest.json").read_text())
    for ci, chrom in enumerate(chroms):
        cpos = np.flatnonzero(a["chrom"] == ci)
        mine = [p for p in want if a["chrom"][p] == ci]
        if not mine:
            continue
        first = int(cpos[0])
        assert np.all(np.diff(cpos) == 1), "positions of a chromosome are not contiguous"
        lines = {p - first: p for p in mine}
        cands = index.candidates()
        with open(index.dir / manifest["chromosomes"][chrom]["summary"]) as f:
            for k, line in enumerate(f):
                p = lines.get(k)
                if p is None:
                    continue
                cid = line[line.index('"candidate_id": "') + 17:].split('"', 1)[0]
                assert cid == cands[p], (cid, cands[p])
                g = GRCH38_RE.findall(line)[-1]
                if g != "null":
                    g = json.loads(g)
                    out[p] = (g["chrom"], g["pos0"] + 1, g["ref"], g["alt"])
    return out


def pon_tag(alleles):
    """{position: [PoN names matched]} with the SNV rule of scripts/filter_panel_of_normals.py."""
    import pysam
    vcfs = [pysam.VariantFile(str(p)) for p in PONS]
    contig = {}
    cache, out = {}, {}
    for p, (chrom, pos1, ref, alt) in alleles.items():
        key = (chrom, pos1, ref, alt)
        if key not in cache:
            hits = []
            kind = pon_rules.record_kind(ref, (alt,))
            for k, (vcf, rule) in enumerate(zip(vcfs, pon_rules.PON_RULES[kind])):
                if rule is None or rule == "indel_af":
                    continue
                if (k, chrom) not in contig:
                    contig[(k, chrom)] = pon_rules.resolve_contig(vcf, chrom)
                if contig[(k, chrom)] is None:
                    continue
                if pon_rules.record_matches(vcf, contig[(k, chrom)], pos1, ref, (alt,), **rule):
                    hits.append(pon_rules.PON_NAMES[k])
            cache[key] = hits
        out[p] = cache[key]
    return out


def summarize(mode, probs, ctx):
    """Metrics of one mode on the test tensors."""
    lab, thr, test = ctx["labels"], ctx["thr"], ctx["test"]
    pr = probs[test]
    lt = lab[test]
    ps = pr[:, 1]
    off = ctx["off"][test]
    am = pr.argmax(1) == 1
    call = ps >= thr
    base = ctx.get("base_probs")
    r = dict(mode=mode, tensors=int(len(lt)))
    r["somatic_ap"] = metrics.average_precision(lt == 1, ps)
    r["germline_ap"] = metrics.average_precision(lt == 2, pr[:, 2])
    for name, m in (("somatic", lt == 1), ("germline", lt == 2), ("non_onref", (lt == 0) & ~off),
                    ("non_offref", (lt == 0) & off)):
        r[name] = dict(n=int(m.sum()), ge_thr=int((call & m).sum()), argmax_somatic=int((am & m).sum()),
                       p_median=float(np.median(ps[m])) if m.any() else None,
                       p_mean=float(np.mean(ps[m])) if m.any() else None)
    r["calls_ge_thr"] = int(call.sum())
    tp_t = int((call & (lt == 1)).sum())
    r["per_tensor"] = dict(tp=tp_t, fp=int(call.sum()) - tp_t,
                           precision=tp_t / max(int(call.sum()), 1), recall=tp_t / max(int((lt == 1).sum()), 1))
    if base is not None:
        bt = base[test][:, 1]
        r["flips_vs_base"] = dict(new_calls=int((call & ~(bt >= thr)).sum()), lost_calls=int((~call & (bt >= thr)).sum()),
                                  somatic_delta_p_median=float(np.median(ps[lt == 1] - bt[lt == 1])) if (lt == 1).any() else None,
                                  germline_delta_p_median=float(np.median(ps[lt == 2] - bt[lt == 2])) if (lt == 2).any() else None)
    tr = metrics.truth_report(ctx["truth"], ctx["matches"], lt, pr, thr, ctx["other"])
    r["truth"] = dict(truth_alleles=tr["truth_alleles"], with_tensor=tr["with_tensor"], ap=tr["ap"],
                      **{k: v for k, v in tr["at_threshold"].items() if k != "threshold"})
    # after a per-tensor PoN: on-reference tensors only (off-reference calls are unplaced, never scored by rtg)
    keep = ~off & ~ctx["pon"][test]
    lab_pon = np.where(keep, lt, -1)
    trp = metrics.truth_report(ctx["truth"], ctx["matches"], lab_pon, pr, thr, ctx["other"])
    r["truth_after_pon"] = dict(ap=trp["ap"], with_tensor=trp["with_tensor"],
                                **{k: v for k, v in trp["at_threshold"].items() if k != "threshold"})
    r["germline_after_pon_ge_thr"] = int((call & keep & (lt == 2)).sum())
    r["non_after_pon_ge_thr"] = int((call & keep & (lt == 0)).sum())
    if ctx["asm"] is not None and len(ctx["asm"]):
        k = ctx["asm"]
        pa = probs[k]
        r["assembly"] = dict(n=int(len(k)), ge_thr=int((pa[:, 1] >= thr).sum()), argmax_somatic=int((pa.argmax(1) == 1).sum()),
                             p_median=float(np.median(pa[:, 1])), chr1_n=int(ctx["asm_chr1"].sum()),
                             chr1_ge_thr=int((pa[ctx["asm_chr1"], 1] >= thr).sum()),
                             chr1_argmax_somatic=int((pa[ctx["asm_chr1"]].argmax(1) == 1).sum()))
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--set", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--modes", nargs="+", default=list(MODES))
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--train-sample", type=int, default=3000)
    ap.add_argument("--limit", type=int, default=0, help="smoke test: only the first N evaluation tensors")
    ap.add_argument("--assembly", action="store_true", help="score the assembly off-reference somatic SNV keys")
    ap.add_argument("--no-pon", action="store_true")
    args = ap.parse_args()
    torch.set_num_threads(args.threads)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    own_cache = out.parent / "index_cache"
    log(f"run {args.run} set {args.set} threads {args.threads} host {os.uname().nodename} job {os.environ.get('SLURM_JOB_ID')}")

    ck_path = RUNS / args.run / "best.pth"
    ck = torch.load(ck_path, map_location="cpu", weights_only=False)
    model = PansomaNetV2.from_checkpoint(ck).to(DEVICE).eval()
    thr = float(ck["val"]["somatic_at_recall"]["0.9"]["threshold"])
    assert "path_count" not in model.config["drop_planes"]
    log(f"checkpoint epoch {ck['epoch']}, val recall-0.9 threshold {thr:.4f}, stats {json.dumps(model.encoder.stats())}")

    index = open_index(TENSORS / args.set / "tensors" / "SNV", own_cache)
    if index.path_count != ck.get("path_count_storage"):
        raise SystemExit(f"ch6 storage {index.path_count} != checkpoint {ck.get('path_count_storage')}")
    a = index.arrays
    pos = index.select(["chr1"], labelled=False)
    ignore = set(ck.get("args", {}).get("ignore_reasons") or [])
    reasons = index.meta["reasons"]
    lab_all = a["eval_label"].astype(np.int64).copy()
    lab_all[np.isin(a["reason"], [k for k, r in enumerate(reasons) if r in ignore])] = -1
    test_pos = pos[lab_all[pos] >= 0]
    cands = index.candidates()

    asm_pos, asm_rows = np.zeros(0, np.int64), []
    if args.assembly:
        rows = [r for r in csv.DictReader(open(ASSEMBLY), delimiter="\t") if r["kind"] == "SNP"]
        by_key = {r["key"]: r for r in rows}
        test_set = set(test_pos.tolist())
        found = {c: p for p, c in enumerate(cands) if c in by_key}
        asm_pos = np.array(sorted(found.values()), np.int64)
        asm_rows = [dict(key=cands[p], event=by_key[cands[p]]["event"], grch38=by_key[cands[p]]["grch38_records"],
                         chrom=index.meta["chroms"][a["chrom"][p]], label=int(a["label"][p]), reason=reasons[a["reason"][p]],
                         off_reference=bool(a["off_reference"][p]), in_test=bool(p in test_set))
                    for p in asm_pos]
        log(f"assembly off-reference somatic SNV keys: {len(rows)} in the table, {len(asm_pos)} with a tensor in this set "
            f"({sum(r['chrom'] == 'chr1' for r in asm_rows)} on chr1)")
    ev_pos = np.union1d(test_pos, asm_pos)
    if args.limit:
        asm_pos, asm_rows = asm_pos[:20], asm_rows[:20]
        ev_pos = np.union1d(test_pos[:args.limit], asm_pos)
    test = np.isin(ev_pos, test_pos)
    labels = lab_all[ev_pos]
    off = a["off_reference"][ev_pos].astype(bool)
    log(f"chr1 tensors {len(pos):,}; test {int(test.sum()):,} ({Counter(labels[test].tolist())}); evaluated {len(ev_pos):,}")

    # training ch6 statistics (mean, median) from the training chromosomes
    stats_file = out / "train_ch6_stats.json"
    if stats_file.exists():
        tstats = json.loads(stats_file.read_text())
    else:
        t0 = time.time()
        tstats = train_ch6_stats(index, ck, args.train_sample)
        tstats["checkpoint_mean"] = float(ck["stats"]["path_count"]["mean"])
        tstats["checkpoint_std"] = float(ck["stats"]["path_count"]["std"])
        stats_file.write_text(json.dumps(tstats, indent=1) + "\n")
        log(f"training ch6: mean {tstats['mean']:.2f} (checkpoint {tstats['checkpoint_mean']:.2f}) median {tstats['median']} "
            f"frac 100 {tstats['frac_100']:.3f} below {tstats['frac_below_100']:.3f} ({time.time() - t0:.0f}s)")
    params = dict(train_mean=int(round(tstats["mean"])), train_median=int(tstats["median"]))

    # PoN tags of the on-reference test tensors' GRCh38 alleles
    pon_file = out / "pon_tags.json"
    if args.no_pon:
        pon = np.zeros(len(ev_pos), bool)
    else:
        if pon_file.exists():
            tags = {int(k): v for k, v in json.loads(pon_file.read_text()).items()}
        else:
            t0 = time.time()
            alleles = grch38_alleles(index, ev_pos[test & ~off])
            tags = pon_tag(alleles)
            pon_file.write_text(json.dumps({str(k): v for k, v in tags.items()}) + "\n")
            log(f"PoN: {len(alleles):,} on-reference test alleles, {sum(bool(v) for v in tags.values()):,} tagged "
                f"({time.time() - t0:.0f}s)")
        pon = np.array([bool(tags.get(int(p))) for p in ev_pos])

    on = {"chr1"}
    truth = somatic_truth(index, chroms=on)
    other = somatic_truth(index, on, True)
    ctx = dict(labels=labels, thr=thr, test=test, off=off, pon=pon, truth=set(truth), other=set(other or {}),
               matches=[set(index.truth_of(p).tolist()) for p in ev_pos[test]],
               asm=np.flatnonzero(np.isin(ev_pos, asm_pos)) if len(asm_pos) else None)
    if ctx["asm"] is not None:
        ctx["asm_chr1"] = np.array([index.meta["chroms"][a["chrom"][p]] == "chr1" for p in ev_pos[ctx["asm"]]])

    pos_file = out / "positions.npy"
    if pos_file.exists() and not np.array_equal(np.load(pos_file), ev_pos):
        raise SystemExit(f"{pos_file} holds other tensors: remove the old outputs of {out} first")
    np.save(pos_file, ev_pos)
    parts = [(index, ev_pos)]
    rng = np.random.default_rng(1)
    test_idx = np.flatnonzero(test)
    donors = np.empty(len(ev_pos), np.int64)
    donors[test_idx] = test_idx[rng.permutation(len(test_idx))]
    donors[~test] = rng.choice(test_idx, int((~test).sum()))

    results, probs_all = {}, {}
    for mode in args.modes:
        f = out / f"probs_{mode}.npy"
        if mode == "train_median" and params["train_median"] == ALL:
            log("train_median is 100: the same as all_H on every cell below the repeat codes; skipped")
            continue
        if mode.endswith("_low") and "low" not in params:
            raise SystemExit("base must run before the _low modes")
        if f.exists() and (out / "profile.npz").exists():
            probs = np.load(f)
            if mode == "base":
                z = np.load(out / "profile.npz")
                prof = (z["site"], z["below"], z["distinct"])
            log(f"{mode}: reused {f.name}")
        else:
            log(f"{mode}: start")
            ds = Counterfactual(parts, mode, params, donors)
            if mode == "base":
                probs, prof = run_mode(model, ds, args, profile=True)
                np.savez(out / "profile.npz", site=prof[0], below=prof[1], distinct=prof[2])
            else:
                probs = run_mode(model, ds, args)
            np.save(f, probs.astype(np.float32))
        if mode == "base":
            ctx["base_probs"] = probs
            site = prof[0]
            offneg = test & off & (labels == 0) & (site > 0)
            if offneg.sum() >= 50:
                params["low"], params["low_rule"] = int(np.median(site[offneg])), f"median site code of {int(offneg.sum())} off-reference test tensors"
            else:
                params["low"], params["low_rule"] = int(np.percentile(site[test & (site > 0)], 5)), "5th percentile of the test tensors' site codes (fewer than 50 off-reference)"
            groups = {"somatic": test & (labels == 1), "germline": test & (labels == 2), "non_onref": test & (labels == 0) & ~off,
                      "non_offref": test & (labels == 0) & off}
            params["site_code"] = {g: dict(n=int(m.sum()), median=float(np.median(site[m])) if m.any() else None,
                                           p10=float(np.percentile(site[m], 10)) if m.any() else None,
                                           frac_site_100=float((site[m] == ALL).mean()) if m.any() else None,
                                           frac_cells_below_100=float(np.mean(prof[1][m])) if m.any() else None,
                                           distinct_codes_pm10_median=float(np.median(prof[2][m])) if m.any() else None)
                                   for g, m in groups.items()}
            log(f"LOW = {params['low']} ({params['low_rule']}); site codes {json.dumps(params['site_code'])}")
        probs_all[mode] = probs
        results[mode] = summarize(mode, probs, ctx)
        s = results[mode]
        log(f"{mode}: AP {s['somatic_ap']:.4f} | somatic >=thr {s['somatic']['ge_thr']}/{s['somatic']['n']} "
            f"(median p {s['somatic']['p_median']:.3f}) | germline >=thr {s['germline']['ge_thr']}/{s['germline']['n']} | "
            f"non on/off >=thr {s['non_onref']['ge_thr']}/{s['non_offref']['ge_thr']} | truth TP {s['truth']['tp']} FP {s['truth']['fp']} "
            f"F1 {s['truth']['f1']:.3f} | after PoN TP {s['truth_after_pon']['tp']} FP {s['truth_after_pon']['fp']} F1 {s['truth_after_pon']['f1']:.3f}"
            + (f" | assembly >=thr {s['assembly']['ge_thr']}/{s['assembly']['n']}" if "assembly" in s else ""))
    report = dict(run=args.run, set=args.set, checkpoint=str(ck_path), epoch=ck["epoch"], threshold_val_recall_0_9=thr,
                  checkpoint_somatic_threshold=float(ck["somatic_threshold"]), haplotypes=index.haplotypes,
                  train_ch6=tstats, params=params, chr1_tensors=int(len(pos)), test_tensors=int(test.sum()),
                  evaluated=int(len(ev_pos)), limit=args.limit, truth_alleles=len(truth),
                  pon_tagged_test=int((pon & test).sum()), results=results)
    if ctx["asm"] is not None:
        for r, k in zip(asm_rows, range(len(asm_rows))):
            j = int(np.flatnonzero(ev_pos == asm_pos[k])[0])
            r.update({f"p_{m}": round(float(probs_all[m][j, 1]), 4) for m in probs_all},
                     **{f"argmax_{m}": int(probs_all[m][j].argmax()) for m in probs_all})
        report["assembly_keys"] = asm_rows
    (out / "report.json").write_text(json.dumps(report, indent=1) + "\n")
    log("done")


if __name__ == "__main__":
    main()
