# Pansoma Pipeline

This document records the end-to-end data generation and modeling workflow. Script names in the old project were not always literal, so the `scripts/` folder exposes clearer entry points while legacy copies remain under `experiments/legacy/`.

## 0. GBZ To GFA

The tensor pipeline needs a GFA with node sequences that matches the GBZ used
for alignment. If only the GBZ is available, export the GFA directly with
`vg`:

```bash
vg convert -f --no-translation \
  /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gbz \
  > /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gfa
```

Run this once per graph. `--no-translation` keeps the internal `vg` node IDs,
which ensures that the exported GFA is compatible with the GAM alignments.

## 1. FASTQ To GAM

Use `vg giraffe` directly. HiFi/long-read data requires `vg >= 1.63`.

```bash
vg giraffe \
  -Z /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gbz \
  -m /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.min \
  -d /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.dist \
  -f /scratch/qfu/COLO829BL_WGS/COLO829BL_1.fq.gz \
  -f /scratch/qfu/COLO829BL_WGS/COLO829BL_2.fq.gz \
  -b default \
  -t 12 \
  -p \
  > /scratch/jshen/data/COLO829T/illumina/GAM/COLO829BL.gam
```

For single-end long reads, use a long-read minimizer/zipcode index pair and
select `-b hifi` for PacBio HiFi or `-b r10` for ONT R10.

## 2. GAM To `.dat/.idx`

Find nodes containing imperfect reads:

```bash
python -u scripts/find_unperfect_nodes.py \
  /scratch/jshen/data/COLO829T/illumina/GAM/COLO829BL.gam \
  --output /scratch/jshen/data/COLO829T/illumina/GAM/unperfect_nodes_COLO829BL.pkl \
  --output_format pickle \
  --milestone 10000000 \
  --threads 12
```

Build the `fast_writer` extension once:

```bash
bash scripts/build_fast_writer.sh
```

Build the packed alignment store:

```bash
python -u scripts/build_dat_idx.py \
  /scratch/jshen/data/COLO829T/COLO829T_ONT/COLO829T_ONT_merged.gam \
  /scratch/jshen/data/COLO829T/COLO829T_ONT/unperfect_nodes_COLO829T_ONT_merged.pkl \
  /scratch/jshen/data/COLO829T/COLO829T_ONT/unperfect_nodes_COLO829T_ONT_merged \
  --milestone 1000000 \
  --threads 12
```

## 3. GRCh38 Node Mapping And Optional AF Annotation

Extract GRCh38 path nodes from GFA and optionally annotate gnomAD AF:

```bash
python -u scripts/build_grch38_path_json.py \
  /scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gfa \
  '^W\tGRCh38\t0\tchr1' \
  -o ./tmp/hprc-v1.1-mc-grch38.d9.genomeAD.chr1.json \
  --vcf /scratch/jshen/data/genomeAD/genomes/gnomad.genomes.v4.1.sites.chr1.vcf.bgz \
  --vcf-contig chr1
```

Filter that JSON by `.idx`, chromosome, and optional truth VCF:

```bash
python -u scripts/filter_node_json.py \
  /scratch/jshen/data/genomeAD/graph_genomeAD_JSON/hprc-v1.1-mc-grch38.d9.genomeAD.chr1.json \
  /scratch/jshen/data/SEQC2/SRR7890824_Tumor/SRR7890824_Ver2.idx \
  /scratch/jshen/data/Pansoma/SRR7890824_Tumor/tmp_INDEL/hprc-v1.1-mc-grch38.d9.genomeAD.chr1.SRR7890824_Ver2.json \
  --txt /scratch/jshen/data/Pansoma/SRR7890824_Tumor/tmp_INDEL/hprc-v1.1-mc-grch38.d9.genomeAD.chr1.SRR7890824_Ver2_nodeIDs.txt \
  --vcf /scratch/jshen/data/SEQC2/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz
```

For testing-set generation across whole chromosomes, build chromosome component and GRCh38 path filters:

```bash
cd /scratch/jshen/Github/Pansoma

GBZ=/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gbz \
GFA=/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.gfa \
OUTDIR=/scratch/jshen/Github/Pansoma/tmp/chr_component_vs_GRCh38_summary \
CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22" \
bash ./scripts/build_chr_node_filters.sh
```

## 4. Build, Label, And Organize Tensors

Generate tensors:

```bash
python -u scripts/generate_testing_tensors.py \
  sample.dat \
  sample.idx \
  output_tensor_dir \
  candidate_nodes.json \
  --chr_nodes chr1.component.nodes.raw.txt chr1.GRCh38_path.nodes.raw.txt \
  --num_workers 8 \
  --variant_type snp \
  --view 0 \
  --min_af 0.08 \
  --shard_size 32768
```

Label tensors against a truth VCF:

```bash
python -u scripts/label_tensors.py \
  output_tensor_dir/variant_summary.ndjson \
  candidate_nodes.json \
  /path/to/truth.vcf.gz \
  --chr chr1 \
  --data-dir output_tensor_dir
```

## 5. Train Model

Training code is bundled under `machine_learning/pansoma_net/`. Use the wrapper template:

```bash
sbatch scripts/slurm/train_pansoma_net.sh
```

The wrapper defaults to `machine_learning/pansoma_net/`. Override paths through environment variables such as `PANSOMA_NET_DIR`, `TRAIN_DATA_PATHS_FILE`, `VAL_DATA_PATHS_FILE`, and `OUT_DIR`.

## 6. Inference

Build testing tensor sets first, then infer using the bundled PansomaNet model code:

```bash
INPUT_DIR=/scratch/jshen/data/Pansoma/HG008_GIAB/AF_HPRC/5ch_testing_data_SNV_chr1 \
CKPT=/scratch/jshen/Github/PansomaNet/HG008T_GIAB_AF-HPRC_CE_Large_Model_V2_weight200/model_e053_f1_0.1733.pth \
OUT_PREFIX=/scratch/jshen/Pansoma_testing_results_V2/HG008_GIAB/AF_HPRC/pansoma_HG008T_WGS_Chr1/pansoma-to_SNV_pansoma_HG008T_WGS_Chr1 \
MAP_JSON=/scratch/jshen/data/Pansoma/HG008_GIAB/AF_HPRC/tmp/hprc-v1.1-mc-grch38.d9.ALL_node_pos.HapMap_HG008T_AF-HPRC.json \
VARIANT_SUMMARY=/scratch/jshen/data/Pansoma/HG008_GIAB/AF_HPRC/5ch_testing_data_SNV_chr1/variant_summary.ndjson \
sbatch scripts/slurm/infer_pansoma_net.sh
```

## 7. Panel Of Normals Filter

After inference, use the four indexed GRCh38 PoNs documented in the main
README to tag recurrent germline or technical-background sites:

```bash
python -u scripts/filter_panel_of_normals.py \
  pansoma_sample.vcf.gz \
  pansoma_sample.pon-tagged.vcf.gz \
  --pon \
    /path/to/pons/gnomad.r2.1.af-ge-0.001.sites.vcf.gz \
    /path/to/pons/dbsnp.b138.non-somatic.sites.vcf.gz \
    /path/to/pons/1000g-pon.sites.vcf.gz \
    /path/to/pons/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.af-ge-0.001.vcf.gz
```

PoN 1/2 use exact position/REF/ALT matching. PoN 3/4 use position matching.
Matched records are tagged by default; add `--drop-matched` to remove them.
