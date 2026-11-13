# COLO829T dilutions, SNV, chr1: new pipeline vs ClairS-TO and DeepSomatic-TO

(`collect.py`; benchmark BED, after our PoN, rtg vcfeval; cells P / R / **F1** of the PASS calls. One threshold per Pansoma model for its four purities. Final PASS thresholds (the user, 2026-11-12): Illumina and fiberseq the validation call threshold - 0.10, ONT the recipe's threshold, the post-PoN validation recall 0.95 of the callable positives.)

| platform | purity | ClairS-TO | DeepSomatic-TO | Pansoma new pipeline, final threshold | Pansoma new pipeline, the recipe's threshold | Pansoma earlier joint (7fc7396) |
|---|---|---|---|---|---|---|
| Illumina | P100 | 0.431 / 0.764 / **0.551** | 0.427 / 0.835 / **0.565** | 0.588 / 0.733 / **0.653** | 0.751 / 0.643 / **0.693** | 0.931 / 0.530 / **0.675** |
| Illumina | P75 | 0.589 / 0.752 / **0.661** | 0.360 / 0.815 / **0.500** | 0.644 / 0.763 / **0.698** | 0.805 / 0.757 / **0.781** | 0.968 / 0.742 / **0.840** |
| Illumina | P50 | 0.542 / 0.686 / **0.606** | 0.346 / 0.779 / **0.479** | 0.624 / 0.737 / **0.675** | 0.791 / 0.736 / **0.762** | 0.966 / 0.733 / **0.833** |
| Illumina | P25 | 0.493 / 0.622 / **0.550** | 0.326 / 0.726 / **0.450** | 0.591 / 0.689 / **0.636** | 0.773 / 0.689 / **0.728** | 0.967 / 0.689 / **0.804** |
| fiberseq | P100 | 0.779 / 0.745 / **0.762** | 0.732 / 0.555 / **0.632** | 0.814 / 0.768 / **0.790** | 0.839 / 0.756 / **0.795** | 0.881 / 0.698 / **0.779** |
| fiberseq | P75 | 0.651 / 0.728 / **0.688** | 0.604 / 0.835 / **0.701** | 0.682 / 0.735 / **0.708** | 0.711 / 0.712 / **0.712** | 0.741 / 0.639 / **0.687** |
| fiberseq | P50 | 0.610 / 0.716 / **0.659** | 0.568 / 0.835 / **0.676** | 0.662 / 0.699 / **0.680** | 0.686 / 0.673 / **0.679** | 0.715 / 0.621 / **0.665** |
| fiberseq | P25 | 0.581 / 0.650 / **0.614** | 0.518 / 0.783 / **0.623** | 0.697 / 0.674 / **0.685** | 0.735 / 0.669 / **0.701** | 0.786 / 0.663 / **0.720** |
| ONT | P100 | 0.749 / 0.684 / **0.715** | 0.582 / 0.695 / **0.633** | 0.795 / 0.720 / **0.756** | the final threshold | 0.723 / 0.637 / **0.677** |
| ONT | P75 | 0.650 / 0.680 / **0.664** | 0.526 / 0.756 / **0.621** | 0.728 / 0.701 / **0.714** | the final threshold | 0.659 / 0.648 / **0.654** |
| ONT | P50 | 0.601 / 0.663 / **0.630** | 0.513 / 0.748 / **0.609** | 0.688 / 0.676 / **0.682** | the final threshold | 0.635 / 0.643 / **0.639** |
| ONT | P25 | 0.557 / 0.573 / **0.565** | 0.493 / 0.680 / **0.572** | 0.714 / 0.664 / **0.688** | the final threshold | 0.670 / 0.661 / **0.666** |

## Every caller and threshold (PASS after the PoN; truth 2,443 SNVs on chr1; `final` marks the Pansoma threshold in use)

| platform | purity | caller | PASS threshold | TP | FP | FN | P | R | F1 | TP, FP, F1 against the recipe's threshold | best F1 of the chr1 curve | recall ceiling |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Illumina | P100 | ClairS-TO | | 1,867 | 2,468 | 576 | 0.431 | 0.764 | 0.551 | | 0.719 | 0.789 |
| Illumina | P100 | DeepSomatic-TO | | 2,039 | 2,740 | 404 | 0.427 | 0.835 | 0.565 | | 0.771 | 0.849 |
| Illumina | P100 | earlier joint | | 1,294 | 96 | 1,149 | 0.931 | 0.530 | 0.675 | | 0.784 | 0.792 |
| Illumina | P100 | new pipeline, the recipe's threshold | 0.9165 | 1,571 | 521 | 872 | 0.751 | 0.643 | 0.693 |  | 0.694 | 0.797 |
| Illumina | P100 | new pipeline, - 0.10, final | 0.8165 | 1,790 | 1,253 | 653 | 0.588 | 0.733 | 0.653 | +219, +732, -0.040 | 0.694 | 0.797 |
| Illumina | P100 | new pipeline, - 0.15 | 0.7665 | 1,853 | 1,609 | 590 | 0.535 | 0.758 | 0.628 | +282, +1,088, -0.065 | 0.694 | 0.797 |
| Illumina | P75 | ClairS-TO | | 1,837 | 1,282 | 606 | 0.589 | 0.752 | 0.661 | | 0.690 | 0.775 |
| Illumina | P75 | DeepSomatic-TO | | 1,990 | 3,533 | 453 | 0.360 | 0.815 | 0.500 | | 0.685 | 0.817 |
| Illumina | P75 | earlier joint | | 1,812 | 60 | 631 | 0.968 | 0.742 | 0.840 | | 0.841 | 0.766 |
| Illumina | P75 | new pipeline, the recipe's threshold | 0.9165 | 1,850 | 447 | 593 | 0.805 | 0.757 | 0.781 |  | 0.830 | 0.769 |
| Illumina | P75 | new pipeline, - 0.10, final | 0.8165 | 1,864 | 1,031 | 579 | 0.644 | 0.763 | 0.698 | +14, +584, -0.082 | 0.830 | 0.769 |
| Illumina | P75 | new pipeline, - 0.15 | 0.7665 | 1,870 | 1,318 | 573 | 0.587 | 0.765 | 0.664 | +20, +871, -0.116 | 0.830 | 0.769 |
| Illumina | P50 | ClairS-TO | | 1,677 | 1,416 | 766 | 0.542 | 0.686 | 0.606 | | 0.628 | 0.706 |
| Illumina | P50 | DeepSomatic-TO | | 1,902 | 3,589 | 541 | 0.346 | 0.779 | 0.479 | | 0.656 | 0.779 |
| Illumina | P50 | earlier joint | | 1,790 | 64 | 653 | 0.966 | 0.733 | 0.833 | | 0.834 | 0.738 |
| Illumina | P50 | new pipeline, the recipe's threshold | 0.9165 | 1,798 | 475 | 645 | 0.791 | 0.736 | 0.762 |  | 0.832 | 0.739 |
| Illumina | P50 | new pipeline, - 0.10, final | 0.8165 | 1,800 | 1,086 | 643 | 0.624 | 0.737 | 0.675 | +2, +611, -0.087 | 0.832 | 0.739 |
| Illumina | P50 | new pipeline, - 0.15 | 0.7665 | 1,802 | 1,375 | 641 | 0.567 | 0.738 | 0.641 | +4, +900, -0.121 | 0.832 | 0.739 |
| Illumina | P25 | ClairS-TO | | 1,520 | 1,563 | 923 | 0.493 | 0.622 | 0.550 | | 0.560 | 0.638 |
| Illumina | P25 | DeepSomatic-TO | | 1,773 | 3,670 | 670 | 0.326 | 0.726 | 0.450 | | 0.642 | 0.726 |
| Illumina | P25 | earlier joint | | 1,682 | 57 | 761 | 0.967 | 0.689 | 0.804 | | 0.805 | 0.692 |
| Illumina | P25 | new pipeline, the recipe's threshold | 0.9165 | 1,682 | 495 | 761 | 0.773 | 0.689 | 0.728 |  | 0.804 | 0.692 |
| Illumina | P25 | new pipeline, - 0.10, final | 0.8165 | 1,684 | 1,165 | 759 | 0.591 | 0.689 | 0.636 | +2, +670, -0.092 | 0.804 | 0.692 |
| Illumina | P25 | new pipeline, - 0.15 | 0.7665 | 1,687 | 1,484 | 756 | 0.532 | 0.691 | 0.601 | +5, +989, -0.127 | 0.804 | 0.692 |
| fiberseq | P100 | ClairS-TO | | 1,819 | 515 | 624 | 0.779 | 0.745 | 0.762 | | 0.760 | 0.803 |
| fiberseq | P100 | DeepSomatic-TO | | 1,357 | 496 | 1,086 | 0.732 | 0.555 | 0.632 | | 0.792 | 0.894 |
| fiberseq | P100 | earlier joint | | 1,705 | 230 | 738 | 0.881 | 0.698 | 0.779 | | 0.790 | 0.795 |
| fiberseq | P100 | new pipeline, the recipe's threshold | 0.6039 | 1,847 | 355 | 596 | 0.839 | 0.756 | 0.795 |  | 0.797 | 0.797 |
| fiberseq | P100 | new pipeline, - 0.10, final | 0.5039 | 1,877 | 430 | 566 | 0.814 | 0.768 | 0.790 | +30, +75, -0.005 | 0.797 | 0.797 |
| fiberseq | P100 | new pipeline, - 0.15 | 0.4539 | 1,889 | 472 | 554 | 0.800 | 0.773 | 0.786 | +42, +117, -0.009 | 0.797 | 0.797 |
| fiberseq | P75 | ClairS-TO | | 1,779 | 952 | 664 | 0.651 | 0.728 | 0.688 | | 0.691 | 0.775 |
| fiberseq | P75 | DeepSomatic-TO | | 2,039 | 1,338 | 404 | 0.604 | 0.835 | 0.701 | | 0.727 | 0.874 |
| fiberseq | P75 | earlier joint | | 1,562 | 545 | 881 | 0.741 | 0.639 | 0.687 | | 0.703 | 0.767 |
| fiberseq | P75 | new pipeline, the recipe's threshold | 0.6039 | 1,740 | 708 | 703 | 0.711 | 0.712 | 0.712 |  | 0.712 | 0.770 |
| fiberseq | P75 | new pipeline, - 0.10, final | 0.5039 | 1,795 | 836 | 648 | 0.682 | 0.735 | 0.708 | +55, +128, -0.004 | 0.712 | 0.770 |
| fiberseq | P50 | ClairS-TO | | 1,749 | 1,118 | 694 | 0.610 | 0.716 | 0.659 | | 0.659 | 0.743 |
| fiberseq | P50 | DeepSomatic-TO | | 2,041 | 1,554 | 402 | 0.568 | 0.835 | 0.676 | | 0.684 | 0.841 |
| fiberseq | P50 | earlier joint | | 1,518 | 605 | 925 | 0.715 | 0.621 | 0.665 | | 0.671 | 0.740 |
| fiberseq | P50 | new pipeline, the recipe's threshold | 0.6039 | 1,644 | 752 | 799 | 0.686 | 0.673 | 0.679 |  | 0.681 | 0.742 |
| fiberseq | P50 | new pipeline, - 0.10, final | 0.5039 | 1,708 | 872 | 735 | 0.662 | 0.699 | 0.680 | +64, +120, +0.001 | 0.681 | 0.742 |
| fiberseq | P25 | ClairS-TO | | 1,588 | 1,144 | 855 | 0.581 | 0.650 | 0.614 | | 0.615 | 0.672 |
| fiberseq | P25 | DeepSomatic-TO | | 1,914 | 1,784 | 529 | 0.518 | 0.783 | 0.623 | | 0.638 | 0.785 |
| fiberseq | P25 | earlier joint | | 1,620 | 440 | 823 | 0.786 | 0.663 | 0.720 | | 0.729 | 0.688 |
| fiberseq | P25 | new pipeline, the recipe's threshold | 0.6039 | 1,635 | 589 | 808 | 0.735 | 0.669 | 0.701 |  | 0.714 | 0.689 |
| fiberseq | P25 | new pipeline, - 0.10, final | 0.5039 | 1,647 | 716 | 796 | 0.697 | 0.674 | 0.685 | +12, +127, -0.015 | 0.714 | 0.689 |
| fiberseq | P25 | new pipeline, - 0.15 | 0.4539 | 1,655 | 789 | 788 | 0.677 | 0.677 | 0.677 | +20, +200, -0.023 | 0.714 | 0.689 |
| ONT | P100 | ClairS-TO | | 1,672 | 561 | 771 | 0.749 | 0.684 | 0.715 | | 0.717 | 0.771 |
| ONT | P100 | DeepSomatic-TO | | 1,697 | 1,218 | 746 | 0.582 | 0.695 | 0.633 | | 0.697 | 0.810 |
| ONT | P100 | earlier joint | | 1,556 | 595 | 887 | 0.723 | 0.637 | 0.677 | | 0.681 | 0.782 |
| ONT | P100 | new pipeline, the recipe's threshold, final | 0.6324 | 1,759 | 454 | 684 | 0.795 | 0.720 | 0.756 |  | 0.756 | 0.792 |
| ONT | P75 | ClairS-TO | | 1,661 | 896 | 782 | 0.650 | 0.680 | 0.664 | | 0.666 | 0.747 |
| ONT | P75 | DeepSomatic-TO | | 1,848 | 1,665 | 595 | 0.526 | 0.756 | 0.621 | | 0.669 | 0.784 |
| ONT | P75 | earlier joint | | 1,583 | 818 | 860 | 0.659 | 0.648 | 0.654 | | 0.665 | 0.759 |
| ONT | P75 | new pipeline, the recipe's threshold, final | 0.6324 | 1,712 | 639 | 731 | 0.728 | 0.701 | 0.714 |  | 0.716 | 0.761 |
| ONT | P50 | ClairS-TO | | 1,620 | 1,077 | 823 | 0.601 | 0.663 | 0.630 | | 0.633 | 0.708 |
| ONT | P50 | DeepSomatic-TO | | 1,828 | 1,736 | 615 | 0.513 | 0.748 | 0.609 | | 0.641 | 0.750 |
| ONT | P50 | earlier joint | | 1,571 | 903 | 872 | 0.635 | 0.643 | 0.639 | | 0.654 | 0.726 |
| ONT | P50 | new pipeline, the recipe's threshold, final | 0.6324 | 1,652 | 750 | 791 | 0.688 | 0.676 | 0.682 |  | 0.694 | 0.728 |
| ONT | P25 | ClairS-TO | | 1,399 | 1,112 | 1,044 | 0.557 | 0.573 | 0.565 | | 0.577 | 0.616 |
| ONT | P25 | DeepSomatic-TO | | 1,661 | 1,708 | 782 | 0.493 | 0.680 | 0.572 | | 0.586 | 0.680 |
| ONT | P25 | earlier joint | | 1,616 | 797 | 827 | 0.670 | 0.661 | 0.666 | | 0.711 | 0.689 |
| ONT | P25 | new pipeline, the recipe's threshold, final | 0.6324 | 1,623 | 649 | 820 | 0.714 | 0.664 | 0.688 |  | 0.714 | 0.688 |

## Models

| run | best epoch | final PASS threshold (p_somatic) | the recipe's threshold | PASS calls without a GRCh38 position at the final threshold (P100 / P75 / P50 / P25) |
|---|---|---|---|---|
| `COLO829_Illumina_SNV_joint_newpipe_e5` | 4 | 0.8165 | 0.9165 | 43 / 43 / 38 / 41 |
| `COLO829_fiberseq_SNV_joint_newpipe_e5` | 3 | 0.5039 | 0.6039 | 30 / 44 / 50 / 46 |
| `COLO829_ONT_SNV_joint_newpipe` | 3 | 0.6324 | 0.6324 | 91 / 130 / 153 / 134 |
