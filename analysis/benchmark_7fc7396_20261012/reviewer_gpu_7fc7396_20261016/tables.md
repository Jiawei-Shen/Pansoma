# RIS GPU jobs of the 2026-10-15 reviewer response: results so far

(`collect.py`; chr1 after the PoN, in the benchmark BED; SNV: PASS F1 (best F1), INDEL: best F1 (PASS F1) with the 71cb51a INDEL PoN rule; one seed each unless stated; '-' = not finished.)


## Matched d9 (4 GPUs for SNV, as d46 / linear) against the earlier d9 runs and d46

| set | kind | d9 earlier | d9 matched | d46 | d46 - d9 matched |
|---|---|---|---|---|---|
| HG008T Illumina | SNV | - | 0.585 | 0.570 | -0.015 |
| HG008T Illumina | INDEL | - | 0.269 | 0.256 | -0.013 |
| HG008T PacBio | SNV | 0.621 | 0.608 | 0.605 | -0.003 |
| HG008T PacBio | INDEL | - | 0.263 | 0.275 | +0.012 |
| HG008T ONT | SNV | 0.596 | 0.600 | 0.575 | -0.025 |
| HG008T ONT | INDEL | - | 0.143 | 0.159 | +0.015 |
| COLO829T Illumina | SNV | 0.778 | 0.572 | 0.783 | +0.212 |
| COLO829T Illumina | INDEL | 0.202 | 0.161 | 0.168 | +0.008 |
| COLO829T fiberseq | SNV | 0.781 | 0.752 | 0.779 | +0.026 |
| COLO829T fiberseq | INDEL | 0.124 | 0.142 | 0.137 | -0.005 |
| COLO829T ONT | SNV | 0.706 | 0.744 | 0.737 | -0.007 |
| COLO829T ONT | INDEL | 0.074 | 0.095 | 0.089 | -0.006 |

## No ch6 (`--drop-planes path_count`) against with ch6

| set | kind | d9 with ch6 | d9 no ch6 | d9 change | d46 with ch6 | d46 no ch6 | d46 change |
|---|---|---|---|---|---|---|---|
| HG008T Illumina | SNV | 0.585 | 0.542 | -0.043 | 0.570 | 0.528 | -0.042 |
| HG008T PacBio | SNV | 0.621 | 0.580 | -0.041 | 0.605 | 0.594 | -0.011 |
| HG008T ONT | SNV | 0.596 | 0.562 | -0.033 | 0.575 | 0.584 | +0.009 |
| COLO829T Illumina | SNV | 0.778 | 0.760 | -0.019 | 0.783 | 0.472 | -0.311 |
| COLO829T fiberseq | SNV | 0.781 | 0.775 | -0.006 | 0.779 | 0.771 | -0.007 |
| COLO829T ONT | SNV | 0.706 | 0.693 | -0.012 | 0.737 | 0.675 | -0.062 |
| HG008T Illumina | INDEL | 0.269 | 0.259 | -0.010 | 0.256 | 0.252 | -0.004 |
| HG008T PacBio | INDEL | 0.263 | 0.229 | -0.034 | 0.275 | 0.245 | -0.030 |
| HG008T ONT | INDEL | 0.143 | 0.133 | -0.011 | 0.159 | 0.158 | -0.001 |
| COLO829T Illumina | INDEL | 0.161 | 0.172 | +0.012 | 0.168 | 0.178 | +0.010 |
| COLO829T fiberseq | INDEL | 0.142 | 0.113 | -0.029 | 0.137 | 0.086 | -0.051 |
| COLO829T ONT | INDEL | 0.095 | 0.082 | -0.013 | 0.089 | 0.102 | +0.013 |

d9 controls: the runs on the GPU count of the no-ch6 run (SNV: HG008T Illumina `HG008_Illumina_SNV_7fc7396_ddp4_b1024`, 4 GPUs; the others `HG008_*_SNV_7fc7396` / `COLO829_*_SNV_h90`, 1 GPU; INDEL `*_INDEL_7fc7396`, 1 GPU).

## `--norm dataset` against `--norm train` (d46, in-sample)

normds: BQ and MAPQ per dataset, predict fitted on chr1 only (2026-10-16); normds re-predicted: the same checkpoints with the checkpoint's statistics (predict of 2026-10-19); normbq: BQ only per dataset, statistics of all chromosomes (2026-10-19). Seed range: the three `--norm train` seeds.

| set | kind | --norm train (seed range) | normds | normds re-predicted | normbq | normbq - train |
|---|---|---|---|---|---|---|
| HG008T Illumina | SNV | 0.570 (0.565-0.589) | 0.564 | - | 0.568 | -0.002 |
| HG008T PacBio | SNV | 0.605 (0.604-0.611) | 0.604 | - | 0.603 | -0.002 |
| HG008T ONT | SNV | 0.575 (0.575-0.598) | 0.548 | 0.594 | 0.576 | +0.001 |
| COLO829T Illumina | SNV | 0.783 (0.704-0.783) | 0.782 | - | 0.782 | -0.001 |
| COLO829T fiberseq | SNV | 0.779 (0.774-0.789) | 0.779 | - | 0.777 | -0.001 |
| COLO829T ONT | SNV | 0.737 (0.718-0.737) | 0.520 | 0.746 | 0.753 | +0.016 |
| HG008T Illumina | INDEL | 0.256 | 0.256 | - | 0.256 | +0.000 |
| COLO829T Illumina | INDEL | 0.168 | 0.172 | - | 0.169 | +0.001 |

## Seeds (SNV PASS F1 (best F1); seed 0 = the existing run, seeds 1 and 2 new)

| set | arm | seed 0 | seed 1 | seed 2 | PASS mean (range) | best mean (range) |
|---|---|---|---|---|---|---|
| HG008T Illumina | linear | 0.469 (0.544) | 0.524 (0.547) | 0.520 (0.550) | 0.504 (0.469-0.524) | 0.547 (0.544-0.550) |
| HG008T Illumina | d9 | 0.585 (0.607) | 0.587 (0.613) | 0.585 (0.605) | 0.586 (0.585-0.587) | 0.608 (0.605-0.613) |
| HG008T Illumina | d46 | 0.570 (0.613) | 0.589 (0.612) | 0.565 (0.611) | 0.575 (0.565-0.589) | 0.612 (0.611-0.613) |
| HG008T PacBio | linear | 0.513 (0.596) | 0.592 (0.604) | 0.588 (0.602) | 0.564 (0.513-0.592) | 0.601 (0.596-0.604) |
| HG008T PacBio | d9 | 0.608 (0.631) | 0.572 (0.612) | 0.595 (0.622) | 0.592 (0.572-0.608) | 0.621 (0.612-0.631) |
| HG008T PacBio | d46 | 0.605 (0.622) | 0.604 (0.607) | 0.611 (0.617) | 0.607 (0.604-0.611) | 0.616 (0.607-0.622) |
| HG008T ONT | linear | 0.537 (0.570) | 0.509 (0.568) | 0.524 (0.549) | 0.523 (0.509-0.537) | 0.562 (0.549-0.570) |
| HG008T ONT | d9 | 0.600 (0.616) | 0.581 (0.606) | 0.579 (0.597) | 0.587 (0.579-0.600) | 0.606 (0.597-0.616) |
| HG008T ONT | d46 | 0.575 (0.604) | 0.596 (0.616) | 0.598 (0.611) | 0.589 (0.575-0.598) | 0.610 (0.604-0.616) |
| COLO829T Illumina | linear | 0.468 (0.554) | 0.386 (0.523) | 0.530 (0.579) | 0.461 (0.386-0.530) | 0.552 (0.523-0.579) |
| COLO829T Illumina | d9 | 0.572 (0.659) | 0.626 (0.709) | 0.566 (0.644) | 0.588 (0.566-0.626) | 0.671 (0.644-0.709) |
| COLO829T Illumina | d46 | 0.783 (0.786) | 0.704 (0.730) | 0.765 (0.769) | 0.751 (0.704-0.783) | 0.762 (0.730-0.786) |
| COLO829T fiberseq | linear | 0.680 (0.682) | 0.700 (0.706) | 0.692 (0.718) | 0.691 (0.680-0.700) | 0.702 (0.682-0.718) |
| COLO829T fiberseq | d9 | 0.752 (0.773) | 0.782 (0.788) | 0.771 (0.773) | 0.768 (0.752-0.782) | 0.778 (0.773-0.788) |
| COLO829T fiberseq | d46 | 0.779 (0.783) | 0.774 (0.782) | 0.789 (0.791) | 0.780 (0.774-0.789) | 0.785 (0.782-0.791) |
| COLO829T ONT | linear | 0.418 (0.511) | 0.427 (0.535) | 0.518 (0.572) | 0.454 (0.418-0.518) | 0.539 (0.511-0.572) |
| COLO829T ONT | d9 | 0.744 (0.745) | 0.749 (0.755) | 0.745 (0.749) | 0.746 (0.744-0.749) | 0.750 (0.745-0.755) |
| COLO829T ONT | d46 | 0.737 (0.739) | 0.718 (0.732) | 0.730 (0.736) | 0.729 (0.718-0.737) | 0.736 (0.732-0.739) |

## E13 on HG008T PacBio d46 (SNV PASS F1 (best F1))

| run | training blocks / model | PASS F1 | best F1 |
|---|---|---|---|
| HG008_PacBio_SNV_d46_7fc7396 | 100 %, 50.4 M | 0.605 | 0.622 |
| HG008_PacBio_SNV_d46_e13_keep10_rep1 | 10 %, rep 1 | 0.572 | 0.579 |
| HG008_PacBio_SNV_d46_e13_keep10_rep2 | 10 %, rep 2 | 0.575 | 0.579 |
| HG008_PacBio_SNV_d46_e13_keep25_rep1 | 25 %, rep 1 | 0.580 | 0.590 |
| HG008_PacBio_SNV_d46_e13_keep25_rep2 | 25 %, rep 2 | 0.583 | 0.585 |
| HG008_PacBio_SNV_d46_e13_keep50_rep1 | 50 %, rep 1 | 0.596 | 0.599 |
| HG008_PacBio_SNV_d46_e13_keep50_rep2 | 50 %, rep 2 | 0.597 | 0.601 |
| HG008_PacBio_SNV_d46_e13_keep75_rep1 | 75 %, rep 1 | 0.595 | 0.622 |
| HG008_PacBio_SNV_d46_e13_keep75_rep2 | 75 %, rep 2 | 0.604 | 0.610 |
| HG008_PacBio_SNV_d46_e13_8p7M | 100 %, 8.7 M | 0.616 | 0.632 |
| HG008_PacBio_SNV_d46_e13_199p6M | 100 %, 199.6 M | 0.603 | 0.608 |

## `--norm dataset` (BQ only per dataset, f344bca): SNV over three seeds

PASS F1 (best F1) per seed; mean (range). d46 seed 0 = `<S>_<P>_SNV_d46_normbq` (submit_normbq.sh), the rest submit_normbq_all.sh.

| set | arm | seed 0 | seed 1 | seed 2 | PASS mean (range) | best mean (range) |
|---|---|---|---|---|---|---|
| HG008T Illumina | linear | - | - | - | - | - |
| HG008T Illumina | d9 | 0.588 (0.608) | 0.588 (0.611) | 0.586 (0.605) | 0.588 (0.586-0.588) | 0.608 (0.605-0.611) |
| HG008T Illumina | d46 | 0.568 (0.613) | 0.585 (0.612) | 0.566 (0.611) | 0.573 (0.566-0.585) | 0.612 (0.611-0.613) |
| HG008T PacBio | linear | - | - | - | - | - |
| HG008T PacBio | d9 | 0.611 (0.630) | 0.555 (0.609) | 0.592 (0.624) | 0.586 (0.555-0.611) | 0.621 (0.609-0.630) |
| HG008T PacBio | d46 | 0.603 (0.623) | 0.609 (0.617) | 0.607 (0.618) | 0.606 (0.603-0.609) | 0.619 (0.617-0.623) |
| HG008T ONT | linear | - | - | - | - | - |
| HG008T ONT | d9 | 0.598 (0.618) | 0.575 (0.604) | 0.581 (0.606) | 0.585 (0.575-0.598) | 0.609 (0.604-0.618) |
| HG008T ONT | d46 | 0.576 (0.601) | 0.605 (0.613) | 0.607 (0.611) | 0.596 (0.576-0.607) | 0.608 (0.601-0.613) |
| COLO829T Illumina | linear | - | - | - | - | - |
| COLO829T Illumina | d9 | 0.566 (0.655) | 0.628 (0.708) | 0.568 (0.646) | 0.587 (0.566-0.628) | 0.670 (0.646-0.708) |
| COLO829T Illumina | d46 | 0.782 (0.787) | 0.705 (0.730) | 0.765 (0.769) | 0.751 (0.705-0.782) | 0.762 (0.730-0.787) |
| COLO829T fiberseq | linear | - | - | - | - | - |
| COLO829T fiberseq | d9 | 0.758 (0.774) | 0.783 (0.790) | 0.771 (0.773) | 0.771 (0.758-0.783) | 0.779 (0.773-0.790) |
| COLO829T fiberseq | d46 | 0.777 (0.783) | 0.773 (0.783) | 0.788 (0.792) | 0.779 (0.773-0.788) | 0.786 (0.783-0.792) |
| COLO829T ONT | linear | - | - | - | - | - |
| COLO829T ONT | d9 | 0.744 (0.744) | 0.752 (0.756) | 0.745 (0.748) | 0.747 (0.744-0.752) | 0.749 (0.744-0.756) |
| COLO829T ONT | d46 | 0.753 (0.754) | 0.726 (0.735) | 0.729 (0.737) | 0.736 (0.726-0.753) | 0.742 (0.735-0.754) |

Arm differences of the three-seed means (PASS F1 / best F1):

| set | d46 - d9 | d9 - linear | d46 - linear |
|---|---|---|---|
| HG008T Illumina | -0.015 / +0.004 | - | - |
| HG008T PacBio | +0.021 / -0.001 | - | - |
| HG008T ONT | +0.011 / -0.001 | - | - |
| COLO829T Illumina | +0.164 / +0.092 | - | - |
| COLO829T fiberseq | +0.009 / +0.007 | - | - |
| COLO829T ONT | -0.011 / -0.008 | - | - |

## `--norm dataset`: INDEL per arm (best F1 (PASS F1), one seed)

| set | linear | d9 | d46 |
|---|---|---|---|
| HG008T Illumina | - | - | 0.256 (0.253) |
| HG008T PacBio | - | - | - |
| HG008T ONT | - | - | - |
| COLO829T Illumina | - | - | 0.169 (0.117) |
| COLO829T fiberseq | - | - | - |
| COLO829T ONT | - | - | - |

## `--norm dataset`: no ch6 on d46 (SNV three seeds each side; INDEL one seed)

| set | kind | with ch6: mean (range) | no ch6: mean (range) | change of the means | no ch6 below every with-ch6 seed |
|---|---|---|---|---|---|
| HG008T Illumina | SNV (PASS) | 0.573 (0.566-0.585) | - | - | - |
| HG008T PacBio | SNV (PASS) | 0.606 (0.603-0.609) | - | - | - |
| HG008T ONT | SNV (PASS) | 0.596 (0.576-0.607) | - | - | - |
| COLO829T Illumina | SNV (PASS) | 0.751 (0.705-0.782) | - | - | - |
| COLO829T fiberseq | SNV (PASS) | 0.779 (0.773-0.788) | - | - | - |
| COLO829T ONT | SNV (PASS) | 0.736 (0.726-0.753) | - | - | - |
| HG008T Illumina | INDEL (best) | 0.256 | - | - | - |
| HG008T PacBio | INDEL (best) | - | - | - | - |
| HG008T ONT | INDEL (best) | - | - | - | - |
| COLO829T Illumina | INDEL (best) | 0.169 | - | - | - |
| COLO829T fiberseq | INDEL (best) | - | - | - | - |
| COLO829T ONT | INDEL (best) | - | - | - | - |
