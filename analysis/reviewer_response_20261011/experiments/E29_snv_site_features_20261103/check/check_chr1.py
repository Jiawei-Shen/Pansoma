"""Check the SNV site-feature PoN tag against scripts/filter_panel_of_normals.py: build chr1 of one set with
build_chrom_snv and compare, allele by allele, with the PanelOfNormals FILTER of an existing PoN-tagged chr1 linear VCF
of that set (vcf.sh output; records with P_SOMATIC >= 0.01). usage: check_chr1.py SET PON_VCF"""
import gzip, sys, time
from collections import Counter
from pathlib import Path
sys.path.insert(0, "/scratch/jshen/Github/Pansoma/tmp/snv_sitefeat_20261103/code/machine_learning")
import json
from pansoma_net_v2 import site_features as sf
SET, PONVCF = sys.argv[1], sys.argv[2]
P = "/scratch/jshen/data/Pansoma/panel_of_normal_VCFs"
dbs = dict(gnomad=f"{P}/af-only-gnomad.hg38.vcf.gz", colors=f"{P}/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz",
           dbsnp=f"{P}/Homo_sapiens_assembly38.dbsnp138.vcf.gz", g1000=f"{P}/1000g_pon.hg38.vcf.gz")
kd = Path(f"/scratch/jshen/data/pansoma_v2_tensors/{SET}/tensors/SNV")
m = json.loads((kd / "labels.manifest.json").read_text())
t0 = time.time()
out = sf.build_chrom_snv(str(kd), "chr1", "chr1_variant_summary.ndjson", m["fasta"], dbs, sys.argv[3] if len(sys.argv) > 3 else "bcftools")
print(f"built chr1: {len(out['candidates']):,} tensors in {time.time() - t0:.0f} s", flush=True)
tag = {}
on = out["vcf_pos"] > 0
import numpy as np
anytag = np.max([out[f"pon_tag_{d}"] for d in sf.SNV_DBS], 0)
for p, r, a, t in zip(out["vcf_pos"][on], out["vcf_ref"][on], out["vcf_alt"][on], anytag[on]):
    tag.setdefault((int(p), r.decode(), a.decode()), set()).add(int(t))
c = Counter()
for l in gzip.open(PONVCF, "rt"):
    if l[0] == "#": continue
    x = l.split("\t")
    k = (int(x[1]), x[3], x[4])
    if len(x[3]) != 1 or len(x[4]) != 1: continue
    vcf = "PanelOfNormals" in x[6]
    s = tag.get(k)
    if s is None: c["vcf allele without a tensor allele"] += 1; continue
    if len(s) > 1: c["tensors of one allele disagree"] += 1
    c[f"vcf {'tagged' if vcf else 'untagged'} / features {'tagged' if 1 in s else 'untagged'}"] += 1
print(dict(c))
print("tagged per db (tensors with an allele):", {d: int((out[f'pon_tag_{d}'][on] == 1).sum()) for d in sf.SNV_DBS},
      "popaf>0:", int((np.maximum(out['popaf_gnomad'], out['popaf_colors'])[on] > 0).sum()), "of", int(on.sum()))
