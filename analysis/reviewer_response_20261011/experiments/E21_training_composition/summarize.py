"""E21 step 3: replay each run's training / validation split from the index caches (build_index.py) and count
what the models were trained on: tensors per class, the -1 tensors by reason, the distinct somatic truth alleles behind
the training tensors, and the AF (VAF) of the somatic / germline / non training tensors.

    python summarize.py            # login node: reads the caches, writes the tables next to this script

The split is train.py's (machine_learning/pansoma_net_v2/train.py main(), imported functions only): chr1 left out;
the training chromosomes chr2-22 split by data.block_split (seed 0, --val-fraction 0.05, --val-block-nodes 20000;
every run of ris_7fc7396_20261012 used these, args.json); training = label >= 0, validation = eval_label >= 0
(off-reference tensors without a truth match count as non there); both without the --ignore-reasons of the recipes
(residual_partial_somatic_truth, allele_partial_somatic_truth). The replayed counts are checked against every
train.log (runs_replay_check.tsv).

Outputs (next to this script):
  set_split_counts.tsv      per set x kind: chr2-22 tensors, the training-block and validation-block tensors by
                            label, the -1 by reason, the training / validation tensors by class
  runs_replay_check.tsv     per run: replayed train / validation counts against the train.log counts
  somatic_truth_alleles.tsv per set x kind: PASS somatic truth alleles of the kind on chr2-22 (all / in the BED),
                            in the validation region, behind training tensors, behind validation tensors
  somatic_truth_unions.tsv  distinct somatic truth alleles behind training tensors, unions over platforms / arms
  vaf_quantiles.tsv         per set x kind x class (training split): n, mean, quantiles of the tensor AF (A1 AF)
  vaf_hist.tsv              the same as a histogram, bins of 0.02 on [0, 1]
  somatic_spectrum.tsv      SBS6 fractions of the SNV truth alleles behind the training tensors, with UV-type and
                            APOBEC-type context fractions (R1.7: the mutational processes the models are trained on)
"""
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = Path("/scratch/jshen/Github/Pansoma")
sys.path.insert(0, str(REPO / "machine_learning"))
from pansoma_net_v2.data import KindIndex, block_split, somatic_truth, validation_truth  # noqa: E402

TENSORS = Path("/scratch/jshen/data/pansoma_v2_tensors")
CACHE = Path("/scratch/jshen/data/reviewer_experiments_20261013/E21_training_composition/index_cache")
RIS = Path("/scratch/jshen/data/pansoma_net_v2_runs/ris_7fc7396_20261012")
SEED, VAL_FRACTION, BLOCK_NODES = 0, 0.05, 20000
IGNORE_REASONS = ("residual_partial_somatic_truth", "allele_partial_somatic_truth")
TRAIN_CHROMS = [f"chr{i}" for i in range(2, 23)]
CLASSES = ("non", "somatic", "germline")
AF_COLUMN = 5  # data.SCALARS: log_coverage, log_site_coverage, log_alt_count, log_ref_count, log_other_count, af, ...
QUANTILES = (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
BINS = np.linspace(0.0, 1.0, 51)
ONLY = set(sys.argv[1:])  # optional: restrict to these sets (testing)

SAMPLES = {"HG008T": ("Illumina", "PacBio", "ONT"), "COLO829T": ("Illumina", "fiberseq", "ONT")}


def sets():
    """(set, sample, platform, arm) of every set in the analysis."""
    out = []
    for sample, plats in SAMPLES.items():
        for p in plats:
            out.append((f"{sample}_{p}", sample, p, "d9"))
            out.append((f"{sample}_{p}_d46", sample, p, "d46"))
            out.append((f"{sample}_{p}_linear", sample, p, "linear"))
    for p in SAMPLES["COLO829T"]:
        for q in ("P75", "P50", "P25"):
            out.append((f"COLO829T_{p}_{q}", "COLO829T", p, f"d9 {q}"))
    for p in ("Illumina", "PacBio", "ONT"):
        out.append((f"HCC1395T_{p}_d46", "HCC1395T", p, "d46"))
    for s in ("HCC1954T", "H1437T", "H2009T"):  # CASTLE PacBio, labels of 2026-10-12 (snapshot)
        out.append((f"{s}_PacBio_d46", s, "PacBio", "d46"))
    return out


def truth_table(set_name):
    """{truth_id: (chrom, vcf_pos, vcf_ref, vcf_alt, kind)} of the set's somatic truth (truth/somatic.graph.tsv)."""
    out = {}
    with open(TENSORS / set_name / "truth" / "somatic.graph.tsv") as f:
        head = f.readline().rstrip("\n").split("\t")
        i = [head.index(c) for c in ("truth_id", "chrom", "vcf_pos", "vcf_ref", "vcf_alt", "kind")]
        for line in f:
            x = line.split("\t", max(i) + 1)  # the keys column can be very long; it is not read
            out[int(x[i[0]])] = (x[i[1]], int(x[i[2]]), x[i[3]], x[i[4]], x[i[5]])
    return out


def recall_rows(set_name):
    with open(TENSORS / set_name / "tensors" / "somatic.recall.tsv") as f:
        return list(csv.DictReader(f, delimiter="\t"))


KIND_TRUTH = {"SNV": ("SNP",), "INDEL": ("INS", "DEL")}


COMP = str.maketrans("ACGT", "TGCA")
FASTA = "/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"  # the labels' reference


def contexts(positions):
    """{(chrom, vcf_pos): 3-mer around it} for 1-based positions. Each chromosome is read from the FASTA in one large
    read at its .fai offset (htslib's small reads crawl on BeeGFS: about 1 MB/s against 16 MB/s)."""
    fai = {}
    with open(FASTA + ".fai") as f:
        for line in f:
            x = line.split("\t")
            fai[x[0]] = tuple(int(v) for v in x[1:5])  # length, offset, bases per line, bytes per line
    out, by_chrom = {}, defaultdict(set)
    for c, p in positions:
        by_chrom[c].add(p)
    with open(FASTA, "rb", buffering=0) as f:
        for c in sorted(by_chrom):
            length, offset, lb, lw = fai[c]
            f.seek(offset)
            raw = f.read((length // lb) * lw + length % lb + lw)
            seq = raw.replace(b"\n", b"").replace(b"\r", b"")[:length].decode().upper()
            for p in by_chrom[c]:
                out[(c, p)] = seq[p - 2:p + 1]
    return out


def spectrum(alleles, ctx_of):
    """SBS6 class counts of SNV truth alleles {(chrom, vcf_pos, ref, alt)} (pyrimidine reference) and two context
    markers: UV-type = C>T with a pyrimidine 5' of the C (CC / TC dinucleotide, the cyclobutane-dimer site), APOBEC-type
    = C>T or C>G in TCW (W = A/T; overlaps UV-type at TCA / TCT), TCW C>G alone (APOBEC, not UV), C>T at CpG."""
    out = defaultdict(int)
    for chrom, pos, ref, alt in alleles:
        if len(ref) != 1 or len(alt) != 1:
            continue
        ctx = ctx_of.get((chrom, pos), "")
        if len(ctx) != 3 or ctx[1] != ref:
            out["context_mismatch"] += 1
            continue
        if ref in "GA":
            ref, alt, ctx = ref.translate(COMP), alt.translate(COMP), ctx.translate(COMP)[::-1]
        out[f"{ref}>{alt}"] += 1
        out["snv"] += 1
        if ref == "C" and alt == "T" and ctx[0] in "CT":
            out["uv_CT_at_dipyrimidine"] += 1
        if ref == "C" and alt in "TG" and ctx[0] == "T" and ctx[2] in "AT":
            out["apobec_TCW"] += 1
        if ref == "C" and alt == "G" and ctx[0] == "T" and ctx[2] in "AT":
            out["apobec_TCW_CG"] += 1  # SBS13-like; UV does not make C>G
        if ref == "C" and alt == "T" and ctx[2] == "G":
            out["CpG_CT"] += 1  # SBS1-like (clock-like deamination of methylated CpG)
    return out


def split(index):
    """Replay of train.py's split: (train positions, validation positions, training-block positions, validation-block
    positions, partial left out train, partial left out val)."""
    a = index.arrays
    pos = index.select(TRAIN_CHROMS, labelled=False)
    tr_all, va_all = block_split(index, pos, VAL_FRACTION, BLOCK_NODES, SEED)
    tr = tr_all[a["label"][tr_all] >= 0]
    va = va_all[a["eval_label"][va_all] >= 0]
    codes = [k for k, r in enumerate(index.meta["reasons"]) if r in IGNORE_REASONS]
    keep_tr, keep_va = ~np.isin(a["reason"][tr], codes), ~np.isin(a["reason"][va], codes)
    return tr[keep_tr], va[keep_va], tr_all, va_all, int((~keep_tr).sum()), int((~keep_va).sum())


def class_counts(labels):
    return {c: int((labels == k).sum()) for k, c in enumerate(CLASSES)}


def main():
    rows_split, rows_truth, rows_q, rows_h, rows_spec = [], [], [], [], []
    per_set = {}  # (set, kind) -> dict for the run check and the unions
    for set_name, sample, platform, arm in sets():
        for kind in ("SNV", "INDEL"):
            d = TENSORS / set_name / "tensors" / kind
            if not (d / "labels.manifest.json").exists():
                print(f"{set_name} {kind}: no labels", flush=True)
                continue
            if ONLY and set_name not in ONLY:
                continue
            if not list(CACHE.glob(f"{set_name}.tensors.{kind}.*.npz")):
                print(f"{set_name} {kind}: no index cache yet (build_index.py), skipped", flush=True)
                continue
            try:
                index = KindIndex(d, cache_dir=CACHE)  # loads the cache of build_index.py (builds it if missing)
            except Exception as e:  # noqa: BLE001
                print(f"{set_name} {kind}: {e}", flush=True)
                continue
            lm = json.loads((d / "labels.manifest.json").read_text())
            a = index.arrays
            tr, va, tr_all, va_all, left_tr, left_va = split(index)
            reasons = index.meta["reasons"]
            neg_tr = defaultdict(int)
            for code, n in zip(*np.unique(a["reason"][tr_all[a["label"][tr_all] < 0]], return_counts=True)):
                neg_tr[reasons[code]] += int(n)
            neg_va = defaultdict(int)
            for code, n in zip(*np.unique(a["reason"][va_all[a["eval_label"][va_all] < 0]], return_counts=True)):
                neg_va[reasons[code]] += int(n)
            ctr = class_counts(a["label"][tr])
            cva = class_counts(a["eval_label"][va])
            n_off = int((a["eval_label"][va] != a["label"][va]).sum())
            rows_split.append(dict(
                set=set_name, sample=sample, platform=platform, arm=arm, kind=kind,
                rules_sha256=lm["rules_sha256"][:12], labels_created=lm["created"],
                ignore_repeat_indels=lm.get("ignore_repeat_indels"),
                chr2_22_tensors=len(tr_all) + len(va_all), training_block_tensors=len(tr_all),
                validation_block_tensors=len(va_all),
                train=len(tr), **{f"train_{c}": ctr[c] for c in CLASSES},
                val=len(va), **{f"val_{c}": cva[c] for c in CLASSES},
                val_offref_as_non=n_off, partial_left_out_train=left_tr, partial_left_out_val=left_va,
                ignored_train_blocks=int((a["label"][tr_all] < 0).sum()),
                ignored_train_blocks_by_reason=json.dumps(dict(sorted(neg_tr.items(), key=lambda x: -x[1]))),
                ignored_val_blocks_by_reason=json.dumps(dict(sorted(neg_va.items(), key=lambda x: -x[1])))))
            # somatic truth alleles behind the training / validation tensors
            som_tr = a["label"][tr] == 1
            ids_tr = set()
            for p in tr[som_tr]:
                ids_tr.update(index.truth_of(p).tolist())
            ids_va = set()
            for p in va[a["eval_label"][va] == 1]:
                ids_va.update(index.truth_of(p).tolist())
            rec = recall_rows(set_name)
            kinds = KIND_TRUTH[kind]
            on = [r for r in rec if r["chrom"] in TRAIN_CHROMS and r["kind"] in kinds and r["passed"] == "True"]
            t_all = somatic_truth(index, chroms=TRAIN_CHROMS)  # PASS, in the BED, of the kind
            vt = validation_truth(index, index.select(TRAIN_CHROMS, labelled=False), t_all, VAL_FRACTION, BLOCK_NODES, SEED)
            outside_val = set(t_all) - vt
            # status of the training-region truth alleles (PASS, in the BED, outside the validation blocks) in the
            # set's recall table: representative allele of a tensor / another allele of a site / filtered (reason) /
            # no candidate; this is what the positive set is drawn from
            status = defaultdict(int)
            for r in rec:
                t = int(r["truth_id"])
                if t in outside_val:
                    st = r["status"] + (f":{r['detail']}" if r["status"] == "filtered" and r["detail"] else "")
                    status[st] += 1
            table = truth_table(set_name)
            keys_tr = {table[t][:4] for t in ids_tr if t in table}
            per_set[(set_name, kind)] = dict(sample=sample, platform=platform, arm=arm, keys_tr=keys_tr,
                                            counts_tr=ctr, counts_va=cva, left=(left_tr, left_va), n_off=n_off)
            rows_truth.append(dict(
                set=set_name, sample=sample, platform=platform, arm=arm, kind=kind,
                truth_pass_chr2_22=len(on), truth_pass_in_bed_chr2_22=len(t_all),
                truth_in_validation_region=len(vt), truth_training_region=len(outside_val),
                train_somatic_tensors=int(som_tr.sum()), truth_behind_train_tensors=len(ids_tr),
                truth_behind_train_tensors_in_bed=len(ids_tr & set(t_all)),
                frac_training_region_truth_with_train_tensor=round(len(ids_tr & outside_val) / len(outside_val), 4)
                if outside_val else None,
                val_somatic_tensors=int((a["eval_label"][va] == 1).sum()), truth_behind_val_tensors=len(ids_va),
                training_region_truth_status=json.dumps(dict(sorted(status.items(), key=lambda x: -x[1])))))
            # AF of the training tensors per class
            af = a["scalars"][:, AF_COLUMN]
            for k, c in enumerate(CLASSES):
                x = af[tr[a["label"][tr] == k]].astype(np.float64)
                q = np.quantile(x, QUANTILES) if len(x) else [None] * len(QUANTILES)
                rows_q.append(dict(set=set_name, sample=sample, platform=platform, arm=arm, kind=kind, cls=c, n=len(x),
                                   mean=round(float(x.mean()), 4) if len(x) else None,
                                   min=round(float(x.min()), 4) if len(x) else None,
                                   **{f"q{int(round(qq * 100)):02d}": (round(float(v), 4) if v is not None else None)
                                      for qq, v in zip(QUANTILES, q)},
                                   frac_af_lt_0_10=round(float((x < 0.10).mean()), 4) if len(x) else None,
                                   frac_af_lt_0_20=round(float((x < 0.20).mean()), 4) if len(x) else None,
                                   frac_af_ge_0_40=round(float((x >= 0.40).mean()), 4) if len(x) else None))
                h, _ = np.histogram(np.clip(x, 0, 1), BINS)
                rows_h.append(dict(set=set_name, kind=kind, cls=c, n=len(x),
                                   **{f"{BINS[i]:.2f}-{BINS[i + 1]:.2f}": int(v) for i, v in enumerate(h)}))
            print(f"{set_name} {kind}: train {ctr} val {cva} truth behind train {len(ids_tr):,}", flush=True)

    def write(name, rows):
        if not rows:
            return
        cols = list(rows[0])
        for r in rows:
            cols += [k for k in r if k not in cols]
        with open(HERE / name, "w", newline="") as f:
            w = csv.DictWriter(f, cols, delimiter="\t")
            w.writeheader()
            for r in rows:
                w.writerow({k: ("" if r.get(k) is None else r[k]) for k in cols})

    write("set_split_counts.tsv", rows_split)
    write("somatic_truth_alleles.tsv", rows_truth)
    write("vaf_quantiles.tsv", rows_q)
    write("vaf_hist.tsv", rows_h)
    snv_sets = [(k, v) for k, v in per_set.items() if k[1] == "SNV"]
    ctx_of = contexts({(c, p) for _, v in snv_sets for c, p, _r, _a in v["keys_tr"]})
    for (set_name, _), v in snv_sets:
        sp = spectrum(v["keys_tr"], ctx_of)
        tot = sp.get("snv", 0) or 1
        rows_spec.append(dict(set=set_name, sample=v["sample"], platform=v["platform"], arm=v["arm"],
                              snv_alleles=sp.get("snv", 0),
                              **{c: round(sp.get(c, 0) / tot, 4) for c in ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")},
                              uv_CT_at_dipyrimidine=round(sp.get("uv_CT_at_dipyrimidine", 0) / tot, 4),
                              apobec_TCW=round(sp.get("apobec_TCW", 0) / tot, 4),
                              apobec_TCW_CG=round(sp.get("apobec_TCW_CG", 0) / tot, 4),
                              CpG_CT=round(sp.get("CpG_CT", 0) / tot, 4),
                              context_mismatch=sp.get("context_mismatch", 0)))
    write("somatic_spectrum.tsv", rows_spec)

    # unions of the distinct somatic truth alleles behind training tensors (allele = chrom, pos, ref, alt of the VCF)
    unions = []
    groups = defaultdict(set)
    for (s, kind), v in per_set.items():
        base_arm = v["arm"].split()[0]
        groups[(v["sample"], kind, f"{v['platform']} {v['arm']}")] |= v["keys_tr"]
        if v["arm"] in ("d9", "d46", "linear"):
            groups[(v["sample"], kind, f"all platforms, {v['arm']}")] |= v["keys_tr"]
            groups[(v["sample"], kind, f"{v['platform']}, all three arms")] |= v["keys_tr"]
            groups[(v["sample"], kind, "all platforms, all three arms")] |= v["keys_tr"]
        if base_arm == "d9" and v["sample"] == "COLO829T":
            groups[(v["sample"], kind, f"{v['platform']} joint P100+P75+P50+P25")] |= v["keys_tr"]
    for (sample, kind, what), keys in sorted(groups.items()):
        unions.append(dict(sample=sample, kind=kind, group=what, distinct_truth_alleles_behind_train_tensors=len(keys)))
    write("somatic_truth_unions.tsv", unions)

    # replay check against every train.log
    check = []
    sys.path.insert(0, str(HERE))
    from collect_runs import find_log, parse_log, runs  # noqa: E402
    for run, arm, sample, platform, kind in runs():
        log, _ = find_log(run)
        if log is None:
            continue
        o = parse_log(log.read_text(errors="replace"))
        sets_ = o.get("sets", [])
        rep_tr = {c: sum(per_set.get((s, kind), {}).get("counts_tr", {}).get(c, 0) for s in sets_) for c in CLASSES}
        rep_va = {c: sum(per_set.get((s, kind), {}).get("counts_va", {}).get(c, 0) for s in sets_) for c in CLASSES}
        left = [sum(per_set.get((s, kind), {}).get("left", (0, 0))[i] for s in sets_) for i in (0, 1)]
        n_off = sum(per_set.get((s, kind), {}).get("n_off", 0) for s in sets_)
        keys = set()
        for s in sets_:
            keys |= per_set.get((s, kind), {}).get("keys_tr", set())
        ok = (rep_tr == o.get("train_counts") and rep_va == o.get("val_counts")
              and left == [o.get("ignored_train"), o.get("ignored_val")] and n_off == o.get("val_offref_as_non"))
        if any((s, kind) not in per_set for s in sets_):
            ok = None  # a set of the run not summarized (no index cache yet)
        check.append(dict(run=run, arm=arm, sample=sample, platform=platform, kind=kind, sets="+".join(sets_),
                          log_train=json.dumps(o.get("train_counts")), replay_train=json.dumps(rep_tr),
                          log_val=json.dumps(o.get("val_counts")), replay_val=json.dumps(rep_va),
                          log_left_out=f"{o.get('ignored_train')}/{o.get('ignored_val')}", replay_left_out=f"{left[0]}/{left[1]}",
                          log_val_offref=o.get("val_offref_as_non"), replay_val_offref=n_off,
                          distinct_truth_alleles_behind_train_tensors=len(keys), match=ok))
    write("runs_replay_check.tsv", check)
    print(f"replay check: {sum(r['match'] is True for r in check)} of {len(check)} runs match their train.log, "
          f"{sum(r['match'] is False for r in check)} differ, {sum(r['match'] is None for r in check)} not checked")


if __name__ == "__main__":
    main()
