# E14b: COLO829T ONT purity joint SNV model, control vs per-read BQ centring (chr1, BED, PASS)

## after population filtering (pon)

| model | purity | P | R | F1 | best F1 | ceiling |
|---|---|---|---|---|---|---|
| control | P100 | 0.858 | 0.639 | 0.733 | 0.757 | 0.779 |
| control | P75 | 0.798 | 0.658 | 0.721 | 0.727 | 0.754 |
| control | P50 | 0.740 | 0.648 | 0.691 | 0.695 | 0.722 |
| control | P25 | 0.779 | 0.659 | 0.714 | 0.724 | 0.684 |
| readnorm | P100 | 0.804 | 0.659 | 0.724 | 0.731 | 0.789 |
| readnorm | P75 | 0.731 | 0.659 | 0.693 | 0.696 | 0.760 |
| readnorm | P50 | 0.690 | 0.643 | 0.666 | 0.669 | 0.727 |
| readnorm | P25 | 0.742 | 0.662 | 0.700 | 0.715 | 0.689 |

## before population filtering (raw)

| model | purity | P | R | F1 | best F1 | ceiling |
|---|---|---|---|---|---|---|
| control | P100 | 0.515 | 0.657 | 0.578 | 0.592 | 0.801 |
| control | P75 | 0.354 | 0.674 | 0.464 | 0.495 | 0.775 |
| control | P50 | 0.288 | 0.666 | 0.402 | 0.459 | 0.742 |
| control | P25 | 0.356 | 0.677 | 0.467 | 0.634 | 0.703 |
| readnorm | P100 | 0.445 | 0.677 | 0.537 | 0.587 | 0.812 |
| readnorm | P75 | 0.287 | 0.676 | 0.403 | 0.463 | 0.783 |
| readnorm | P50 | 0.251 | 0.661 | 0.364 | 0.423 | 0.750 |
| readnorm | P25 | 0.327 | 0.680 | 0.441 | 0.610 | 0.709 |

## Training

- control: done; best epoch 8 (f1 0.7251, somatic threshold 0.8417853116989136)
- readnorm: done; best epoch 6 (f1 0.7095, somatic threshold 0.879496693611145)
