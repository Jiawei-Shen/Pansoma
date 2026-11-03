# ClairS-TO and DeepSomatic tumor-only on COLO829T P100 / P75 / P50 / P25 (chr1), against Pansoma

The user (2026-10-29): compare the tumor-only callers at every purity; run what is missing (DeepSomatic on the
mixtures), on the epigenome H100s. Results: `tables.md` (`collect.py`), every score in `summary.json`.

## Call sets

| caller | P100 | P75 / P50 / P25 |
|---|---|---|
| ClairS-TO | the user's runs, 0.4.1 (`outputs/clairs-to_COLO829T_{WGS,Pacbio,ONT}`; Illumina model `ilmn`) | the user's runs, 0.4.2 (`outputs/clairs-to_INDEL_COLO829_{WGS,Pacbio,ONT}_chr1_P*`; Illumina `ilmn_ssrs`, ONT `--indel_min_af 0.1`) |
| DeepSomatic | the user's runs, 1.9.0 (`outputs/deepsomatic_COLO829T_*`) | new: 1.9.0-gpu, `run_deepsomatic.sh` (`submit_deepsomatic.sh`; jobs 592324-592352, exec-419 H100, 16 shards, 25-37 min each) |

All on chr1 of the same BAMs: the SMaHT COLO829T BAMs (P100) and the in-silico mixture BAMs
`Dataset/COLO829T_Tumor_Normal_inSilico_Mix/*/MIX_p{75,50,25}_*.mixed.sorted.bam`. DeepSomatic as the user's P100
runs: the platform's `*_TUMOR_ONLY` model, `GRCh38.primary.fasta`, `--regions chr1`, the default PoN filtering
(postprocess with 4 workers).

## Scoring

E8's `score_comparator.py` unchanged (`score_comparator_ris.py` swaps its wangcluster paths for the RIS copies;
`score.sh` / `submit.sh`): truth and BED of the tensor set's labels, our PoN (71cb51a rule, the files of `ris/vcf.sh`),
rtg vcfeval. Reported: the caller's PASS records after our PoN, and the best point of the curve of all its records
after the PoN (by QUAL). Check: the P100 scores reproduce E8's wangcluster scores exactly (both callers, SNV and INDEL;
`tables.md`, last section).

## Pansoma columns

- joint: one model per platform on P100-P25 (INDEL `COLO829_<P>_INDEL_joint_ppr8`; SNV `COLO829_<P>_SNV_joint_7fc7396_gpu4`).
- P100 model: the P100 single-set d9 models, seeds 0-2, on chr1 of each mixture (`../p100_on_dilutions/submit.sh`):
  INDEL `COLO829_<P>_INDEL_d9_ppr8{,_s1,_s2}`; SNV `COLO829_<P>_SNV_d9_normbq{,_s1,_s2}` (Illumina at epochs 6 / 5 / 6,
  `*_sweep_e06`, `_s1_sweep_e05`, `_s2_sweep_e06`). These models never saw a mixture, so their dilution scores carry
  no read-source BQ cue (the joint models' Illumina / ONT dilution scores likely do).
