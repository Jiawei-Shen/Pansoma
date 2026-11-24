# Response Figure 1.3 (ver4): with / without the haplotype-count channel, plotted values

COLO829T single-tumor SNV models of the current pipeline (`snv_newpipe_20261109`), chr1 (held out), benchmark regions, 2,443 truth SNVs in every cell. PASS calls. P / R / F1 = after the population filter (four resource VCFs); F1 before = the same PASS calls before the filter. Values are computed from the TP / FP / FN counts of `eval_bed/report.json`; the deltas are differences of the unrounded values, so a delta can differ by 0.001 from the difference of two rounded cells.

## Table A. Every plotted value

| COLO829T chr1 SNV | Graph | Model | TP / FP / FN (after) | Precision | Recall | F1 | F1 before the filter |
|---|---|---|---|---|---|---|---|
| Illumina | HPRC v2.1 | With channel, seed 0 | 1,933 / 476 / 510 | 0.802 | 0.791 | 0.797 | 0.797 |
| Illumina | HPRC v2.1 | With channel, seed 1 | 1,951 / 985 / 492 | 0.665 | 0.799 | 0.725 | 0.728 |
| Illumina | HPRC v2.1 | With channel, seed 2 | 1,966 / 408 / 477 | 0.828 | 0.805 | 0.816 | 0.817 |
| Illumina | HPRC v2.1 | Without channel, seed 0 | 1,936 / 886 / 507 | 0.686 | 0.792 | 0.735 | 0.735 |
| Illumina | HPRC v1.1 | With channel, seed 0 | 1,955 / 656 / 488 | 0.749 | 0.800 | 0.774 | 0.774 |
| Illumina | HPRC v1.1 | With channel, seed 1 | 1,939 / 1,051 / 504 | 0.648 | 0.794 | 0.714 | 0.716 |
| Illumina | HPRC v1.1 | With channel, seed 2 | 1,888 / 900 / 555 | 0.677 | 0.773 | 0.722 | 0.723 |
| Illumina | HPRC v1.1 | Without channel, seed 0 | 1,958 / 1,240 / 485 | 0.612 | 0.801 | 0.694 | 0.693 |
| Fiber-seq | HPRC v2.1 | With channel, seed 0 | 1,920 / 334 / 523 | 0.852 | 0.786 | 0.818 | 0.819 |
| Fiber-seq | HPRC v2.1 | With channel, seed 1 | 1,920 / 322 / 523 | 0.856 | 0.786 | 0.820 | 0.820 |
| Fiber-seq | HPRC v2.1 | With channel, seed 2 | 1,940 / 334 / 503 | 0.853 | 0.794 | 0.823 | 0.823 |
| Fiber-seq | HPRC v2.1 | Without channel, seed 0 | 1,899 / 311 / 544 | 0.859 | 0.777 | 0.816 | 0.817 |
| Fiber-seq | HPRC v1.1 | With channel, seed 0 | 1,926 / 336 / 517 | 0.851 | 0.788 | 0.819 | 0.820 |
| Fiber-seq | HPRC v1.1 | With channel, seed 1 | 1,915 / 319 / 528 | 0.857 | 0.784 | 0.819 | 0.820 |
| Fiber-seq | HPRC v1.1 | With channel, seed 2 | 1,910 / 280 / 533 | 0.872 | 0.782 | 0.825 | 0.826 |
| Fiber-seq | HPRC v1.1 | Without channel, seed 0 | 1,940 / 428 / 503 | 0.819 | 0.794 | 0.806 | 0.808 |

## Table B. Deltas of panel b (with, seed 0 − without) and the seed spread

| COLO829T chr1 SNV | Graph | ΔP | ΔR | ΔF1 | ΔF1 before the filter | F1 with: range of 3 seeds (mean) | F1 without | Without vs the seed range (P / R / F1 / F1 before) | ΔF1, mean of 3 seeds − without |
|---|---|---|---|---|---|---|---|---|---|
| Illumina | HPRC v2.1 | +0.116 | −0.001 | +0.061 | +0.062 | 0.725–0.816 (0.779) | 0.735 | inside / inside / inside / inside | +0.044 |
| Illumina | HPRC v1.1 | +0.136 | −0.001 | +0.079 | +0.082 | 0.714–0.774 (0.736) | 0.694 | below / above / below / below | +0.042 |
| Fiber-seq | HPRC v2.1 | −0.007 | +0.009 | +0.001 | +0.002 | 0.818–0.823 (0.820) | 0.816 | above / below / below / below | +0.004 |
| Fiber-seq | HPRC v1.1 | +0.032 | −0.006 | +0.012 | +0.012 | 0.819–0.825 (0.821) | 0.806 | below / above / below / below | +0.014 |

## Table C. Range of the three with-channel seeds per metric

| COLO829T chr1 SNV | Graph | Precision, min–max | Recall, min–max | F1, min–max | F1 before, min–max |
|---|---|---|---|---|---|
| Illumina | HPRC v2.1 | 0.665–0.828 | 0.791–0.805 | 0.725–0.816 | 0.728–0.817 |
| Illumina | HPRC v1.1 | 0.648–0.749 | 0.773–0.800 | 0.714–0.774 | 0.716–0.774 |
| Fiber-seq | HPRC v2.1 | 0.852–0.856 | 0.786–0.794 | 0.818–0.823 | 0.819–0.823 |
| Fiber-seq | HPRC v1.1 | 0.851–0.872 | 0.782–0.788 | 0.819–0.825 | 0.820–0.826 |

Largest |Δ| in panel b: 0.136.

Checks done by `pathcount_data.py` (asserted):
- truth SNVs (TP + FN = baseline = truth_records) = 2443 in all 16 reports, both stages
- P / R / F1 from the counts equal the 4-decimal values of every report (32 blocks)
- truth VCF, BED, kind SNV, chr1, four population resources identical in all 16 reports
- args.json: with-seeds differ from seed 0 only in output / seed; the without model only in output / drop_planes = [path_count]
- PASS rule identical: p_somatic >= <thr> (validation post-PoN recall 0.95 of the callable positives (val.somatic_at_recall_postpon_callable))

Runs: `/scratch/jshen/data/pansoma_net_v2_runs/snv_newpipe_20261109/COLO829_{Illumina,fiberseq}_SNV_{d46,d9}_newpipe{,_s1,_s2}` (with) and `COLO829_{Illumina,fiberseq}_SNV_nopc_{d46,d9}_newpipe` (without); VCF dirs `vcf_chr1_COLO829T_<P>_d46` (v2.1) and `vcf_chr1_COLO829T_<P>` (v1.1).
