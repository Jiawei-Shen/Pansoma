"""PoN features (pon_features.py) of every INDEL truth allele and every INDEL tensor of one chromosome.

    python indel_scan.py chr1   -> scan/chr1.json

Same inputs and allele construction as ../pon_scan.py (truth: the VCF records as they are; tensors: the label
record's GRCh38 allele, left-aligned and padded as pansoma_net_v2.linear_vcf does; unplaced tensors kept apart).
Each PoN is streamed once over the chromosome and every record starting at a position of interest is kept.
"""
import json, sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
sys.path.insert(0, str(REPO / "machine_learning"))
sys.path.insert(0, str(HERE))
from pansoma_net_v2.linear_vcf import Genome  # noqa: E402
from pon_features import PON_DIR, PON_FILES, features, pon_record  # noqa: E402

pysam.set_verbosity(0)
DATA = Path("/scratch/jshen/data/pansoma_v2_tensors")
FASTA = "/scratch/jshen/data/HapMap/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
TRUTH_VCF = {"HG008T": "/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz",
             "COLO829T": str(DATA / "COLO829T_truth/COLO829T_somatic_snv_indel.vcf.gz")}
SETS = ["HG008T_PacBio", "HG008T_ONT", "HG008T_Illumina", "COLO829T_fiberseq", "COLO829T_ONT", "COLO829T_Illumina"]


def main(chrom):
    genome = Genome(FASTA)
    keys = {}

    def key_index(pos, ref, alt):
        return keys.setdefault((pos, ref.upper(), alt.upper()), len(keys))

    truth = []
    for sample, path in TRUTH_VCF.items():
        with pysam.VariantFile(path) as vcf:
            for v in vcf.fetch(chrom):
                for alt in v.alts:
                    if len(v.ref) == 1 and len(alt) == 1:
                        continue
                    truth.append(dict(sample=sample, vcf_pos=v.pos, vcf_ref=v.ref, vcf_alt=alt, key=key_index(v.pos, v.ref, alt)))
    tensors = Counter()
    for name in SETS:
        with open(DATA / name / "tensors" / "INDEL" / f"{chrom}_labels.ndjson") as stream:
            for line in stream:
                r = json.loads(line)
                lin = r.get("grch38")
                group = (name, r["label"], r["reason"])
                if not lin:
                    tensors[group + (-1,)] += 1
                    continue
                ev = r["candidate_id"].split(":")[2]
                c, pos0, ref, alt = lin["chrom"], lin["pos0"], lin["ref"].upper(), lin["alt"].upper()
                pos0, ref, alt, _ = genome.left_align(c, pos0, ref, alt, ev)
                pos, vref, valt = genome.vcf_allele(c, pos0, ref, alt, ev)
                tensors[group + (key_index(pos, vref, valt),)] += 1

    positions = {k[0] for k in keys}
    at = defaultdict(lambda: defaultdict(list))  # pos -> PoN -> records
    for p, fname in enumerate(PON_FILES):
        with pysam.VariantFile(f"{PON_DIR}/{fname}") as pon:
            contig = chrom if chrom in pon.header.contigs else chrom[3:]
            if contig not in pon.header.contigs:
                continue
            for rec in pon.fetch(contig):
                if rec.pos in positions:
                    at[rec.pos][p].append(pon_record(rec, p))
    feat = [None] * len(keys)
    for (pos, ref, alt), i in keys.items():
        feat[i] = features(at.get(pos, {}), ref, alt)

    agg = Counter()
    for (name, label, reason, i), c in tensors.items():
        agg[(name, label, reason, "unplaced" if i < 0 else ",".join(map(str, feat[i])))] += c
    out = dict(chrom=chrom, truth=[dict(t, feat=feat[t.pop("key")]) for t in truth],
               tensors=[[*k, v] for k, v in sorted(agg.items())])
    (HERE / "scan").mkdir(exist_ok=True)
    (HERE / "scan" / f"{chrom}.json").write_text(json.dumps(out))
    print(chrom, f"{len(truth):,} truth INDEL alleles, {sum(tensors.values()):,} INDEL tensors, {len(keys):,} alleles", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
