# E33: HCC1395T in-sample models on GRCh38-only, HPRC v1.1 and HPRC v2.1 (2026-11-09)

For R1.2a: one model per platform x graph, HCC1395T only (the one CASTLE line with GRCh38-only and v1.1 sets), trained
on chr2-22 without chr20, tested on chr1 and chr20, seed 0.

- SNV: PacBio HiFi and Illumina with the AF 0.05 recipe (`_snvaf05` sets, SNV site features, non tensors with AF <= 0.1
  sampled 25% per epoch at loss weight 1; scored at rule 2 and rule 1); ONT with the normfix recipe (AF 0.08 sets, no
  site features, every negative; scored at its PASS).
- INDEL: the final INDEL recipe (PASS ppr8), on the AF 0.08 sets, all three platforms.

Job scripts, design notes and results (`results.md`, `results_indel.md`, from `compare_graphs.py`) are in
`tmp/hcc1395_linear_d9_20261104/train_af05/` on the old cluster.

## Five SNV models trained on RIS

Illumina v1.1 / v2.1 and ONT GRCh38-only / v1.1 / v2.1 run on RIS's epigenome H100s instead of tequila's two H100s.
The code is the same: origin/develop's `machine_learning/pansoma_net_v2` Python code last changed in 1497adb, the code tequila
runs. RIS trains and predicts; the VCFs and vcfeval run on the old cluster, as for the other models.

1. `to_ris.sh` (wangcluster, in tmux): rsync the five sets' SNV tensors (4.4 TB apparent, 613 GB on disk), the HCC1395T
   truth and the two Illumina site-feature files to storage2, then compare file counts and bytes.
2. `ris_submit.sh` (RIS login node, after `git pull`): checks the code and the data, then submits the five trainings
   (`run.sh`, which predicts chr1) and the five chr20 predictions. `ris/submit.sh` bind-mounts
   `Dataset/CASTLE/<line>/truth` at the label manifests' `/scratch/jshen/data/CASTLE/<line>/truth`.
3. `from_ris.sh` (wangcluster, when the `pnv2_pred20_*` jobs are done): copies the finished runs back without
   `last.pth`, cancels the held tequila copies (474803, 474807, 474811, 474815, 474819 and their VCF jobs) and submits
   the chr1 / chr20 VCF jobs.
