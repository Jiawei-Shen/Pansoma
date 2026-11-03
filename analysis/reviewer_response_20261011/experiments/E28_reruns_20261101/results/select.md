| run | epochs (rerun / orig) | reproduces original (max Δloss; max Δselection; first differing epoch) | removed off-path negatives | original epoch (score) | rerun: selected epoch with off-path (score; score without) | rerun: selected epoch without off-path (score) | same epoch |
|---|---|---|---|---|---|---|---|
| release_s0_Illumina_INDEL | 12 / 12 | yes, bit-identical | 8813 | 2 (0.5449) | 2 (0.5449; 0.5523) | 3 (0.5842) | False |
| release_s0_PacBio_INDEL | 12 / 12 | no (4.2e-04; 0.0122; epoch 1) | 31560 | 7 (0.5734) | 6 (0.5808; 0.5978) | 6 (0.5978) | True |
| release_s0_ONT_INDEL | 12 / 12 | yes, bit-identical | 53317 | 5 (0.5335) | 5 (0.5335; 0.5399) | 5 (0.5399) | True |
| release_s1_Illumina_INDEL | 12 / 12 | no (1.0e-04; 0.0034; epoch 1) | 8927 | 1 (0.5194) | 1 (0.5194; 0.5233) | 6 (0.5887) | False |
| release_s1_PacBio_INDEL | 12 / 12 | no (1.6e-04; 0.0085; epoch 1) | 32686 | 6 (0.5645) | 6 (0.5612; 0.5802) | 8 (0.5829) | False |
| release_s1_ONT_INDEL | 12 / 12 | yes, bit-identical | 53774 | 6 (0.5433) | 6 (0.5433; 0.5456) | 7 (0.5463) | False |
| release_s2_Illumina_INDEL | 12 / 12 | no (3.1e-05; 0.0017; epoch 1) | 9314 | 4 (0.5139) | 4 (0.5142; 0.5677) | 6 (0.5779) | False |
| release_s2_PacBio_INDEL | 1 / 12 | no (6.0e-05; 0.0059; epoch 1) | 32576 | 8 (0.5667) | 1 (0.4949; 0.5068) | 1 (0.5068) | True |
| release_s2_ONT_INDEL | 0 / 12 | running or pending | | | | | |
