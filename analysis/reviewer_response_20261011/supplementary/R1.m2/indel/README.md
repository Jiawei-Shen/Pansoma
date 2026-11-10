# Supplementary Fig. R1.m2, INDEL examples 1–10

Ten chr1 somatic truth INDEL tensors whose 101-column window spans several graph nodes, as INDEL companions of the SNV
example in `../fig_R1.m2_tensor_across_nodes.png`. Each is an exact match to a truth INDEL (label reason
`representative_allele_is_somatic_truth`, FILTER PASS, inside the benchmark BED). In each, the ALT and REF reads pass
through different nodes somewhere in the window.

Panels are drawn in the order of channels 0, 6, 1, 2, 3, 4, 5, 7, so the path count (channel 6) is panel 2. Panel
titles keep the 0-based channel number. Rows are reads, grouped by the allele they carry at the site (A1 = ALT, then
REF, then OTHER). Columns are the 101 positions of the window; the site starts at column 50.

## Files

| File | What |
|---|---|
| `make_figures.sh` | runs the four steps below (`SKIP_SELECT=1` reuses `candidates.tsv`) |
| `select_examples.py` → `candidates.tsv` | every exact, PASS, in-BED somatic truth INDEL tensor on chr1 of the twelve HPRC sets (HG008T / COLO829T × Illumina, PacBio HiFi or Fiber-seq, ONT × HPRC v1.1 d9 / v2.1 d46) with >= 5 ALT and REF rows: 1,661 tensors |
| `choose.py` → `picks.tsv` | the ten examples (rules below) |
| `annotate.py` → `annotations.md` | caption facts: GRCh38-path nodes under the window and, per branch column, the channel-6 codes of the ALT and REF rows decoded to haplotype counts and matched to nodes |
| `fig_R1.m2_INDEL_NN.png` | the figures (`scripts/visualize_tensor.py --panel-order 0,6,1,2,3,4,5,7`) |

Selection (`choose.py`): >= 10 ALT and >= 10 REF rows; the window over >= 3 GRCh38-path nodes; at least one column
outside the site where >= 80% of ALT rows and >= 80% of REF rows each share one path-count code, the two codes differ by
>= 10 and neither is the other allele's code within 3 columns (different nodes, not the same nodes shifted by the
INDEL). Of the tensors meeting the row and node thresholds, 61 had such a code difference and 27 passed the shift
test; ranked by the number of such columns, then
GRCh38-path nodes under the window, then ALT rows; one example per truth INDEL; first the best of each set, then a
second from the sets with more examples. Seven sets contributed: HG008T × all six, and COLO829T Fiber-seq v1.1 (the
other COLO829T sets had no tensor passing the rules with a new truth INDEL).

Haplotype counts are exact for HPRC v1.1 (H = 90). For HPRC v2.1 (H = 464), codes 1–50 are square-root binned, so the
counts below are the graph's own counts of the matched node (graph index sqlite), not decoded codes.

## Captions

*Supplementary Fig. R1.m2, INDEL 1.* HG008T Illumina, HPRC v1.1 (d9). Somatic insertion chr1:213,363,821 C>CT
(candidate 2543944:0:INS:>T); 32 ALT, 78 REF and 3 other reads (AF 0.28). The window covers 14 GRCh38-path nodes, most of
them 4-bp nodes of a TTCT repeat. At columns 41–43, all 32 ALT rows traverse node 2543941 (3 bp, on the GRCh38 path,
carried by 44 of the 90 haplotypes) and 85% of REF rows the competing 3-bp node 2543942 (36 haplotypes).

*Supplementary Fig. R1.m2, INDEL 2.* HG008T ONT, HPRC v1.1 (d9). Somatic deletion chr1:188,508,320 GT>G (candidate
2264449:0:DEL:T>); 26 ALT, 13 REF and 7 other reads (AF 0.57). The window covers 6 GRCh38-path nodes. At column 90,
40 columns right of the site, all 26 ALT rows traverse the single-base node 2264452 (on the GRCh38 path, 31 of the 90
haplotypes) and 92% of REF rows the competing single-base node 2264453 (59 haplotypes); 31 + 59 = 90.

*Supplementary Fig. R1.m2, INDEL 3.* HG008T PacBio HiFi, HPRC v1.1 (d9). Somatic deletion chr1:218,103,815 CA>C
(candidate 3693578:0:DEL:A>); 45 ALT, 37 REF and 16 other reads (AF 0.37). The window covers 10 GRCh38-path nodes.
At columns 47–49, just left of the site, 98% of ALT rows are on the GRCh38-path node 3693574 (82 bp, all 90
haplotypes), while 84% of REF rows traverse a three-base branch off the GRCh38 path (nodes 3693576 and 3693577, 59
haplotypes).

*Supplementary Fig. R1.m2, INDEL 4.* HG008T Illumina, HPRC v2.1 (d46). Somatic insertion chr1:159,111,805 A>AC
(candidate 19704905:0:INS:>C); 19 ALT, 123 REF and 1 other reads (AF 0.13). The window covers 13 GRCh38-path nodes.
At column 49, the base before the insertion, all 19 ALT rows traverse the single-base node 19704903 (on the GRCh38 path,
103 of the 464 haplotypes) and 88% of REF rows the competing single-base node 19704904 (345 haplotypes).

*Supplementary Fig. R1.m2, INDEL 5.* HG008T Illumina, HPRC v2.1 (d46). Somatic insertion chr1:233,387,578 C>CA
(candidate 10089773:17:INS:>A); 100 ALT and 103 REF reads (AF 0.49), of which the 200-row tensor holds 99 and 101. The window covers 9 GRCh38-path nodes. At column 53,
all 99 ALT rows traverse the single-base node 10089775 (off the GRCh38 path, 153 of the 464 haplotypes) and 85% of REF
rows the competing GRCh38-path node 10089774 (311 haplotypes); 153 + 311 = 464.

*Supplementary Fig. R1.m2, INDEL 6.* HG008T ONT, HPRC v2.1 (d46). Somatic insertion chr1:189,870,594 C>CA (candidate
8181742:24:INS:>A); 14 ALT, 29 REF and 3 other reads (AF 0.30). The window covers 8 GRCh38-path nodes. At column 53,
all 14 ALT rows traverse the single-base node 8181744 (off the GRCh38 path, 112 of the 464 haplotypes) and 97% of REF
rows the competing GRCh38-path node 8181743 (344 haplotypes).

*Supplementary Fig. R1.m2, INDEL 7.* HG008T ONT, HPRC v2.1 (d46). Somatic insertion chr1:191,234,698 T>TTA (candidate
8237854:0:INS:>TA) in a TA repeat; 14 ALT, 33 REF and 1 other reads (AF 0.29). The window covers 46 GRCh38-path nodes,
most of them single bases. At column 65, all 14 ALT rows traverse node 8237873 (146 of the 464 haplotypes) and 85% of
REF rows the competing node 8237871 (284 haplotypes). At columns 81–90, the ALT rows are on nodes carried by 463–464
haplotypes, while at least 76% of REF rows are still on TA-repeat branch nodes carried by about a quarter of the
haplotypes (path-count code 27, 113–121 paths).

*Supplementary Fig. R1.m2, INDEL 8.* HG008T PacBio HiFi, HPRC v2.1 (d46). Somatic insertion chr1:194,115,661 G>GCA
(candidate 8349226:45:INS:>CA), one more unit of a CA repeat; 27 ALT, 65 REF and 1 other reads (AF 0.29). The window
covers 14 GRCh38-path nodes. Each CA unit after the site is its own 2-bp node, carried by 210, 243, 299, 405 and 463 of
the 464 haplotypes (nodes 8349232–8349236). At columns 54–58, all ALT rows are on the units carried by 210–299
haplotypes and at least 86% of REF rows on nodes carried by 405–463.

*Supplementary Fig. R1.m2, INDEL 9.* HG008T PacBio HiFi, HPRC v2.1 (d46). Somatic deletion chr1:240,187,269 CA>C
(candidate 11316111:0:DEL:A>); 11 ALT and 42 REF reads (AF 0.21). The window covers 9 GRCh38-path nodes. At columns
68–72, all 11 ALT rows traverse the 5-bp node 11316117 (AAAAG, on the GRCh38 path, 230 of the 464 haplotypes), while
95% of REF rows are on the node that follows it, 11316118 (all 464 haplotypes, also starting AAAAG): the REF reads
skip node 11316117.

*Supplementary Fig. R1.m2, INDEL 10.* COLO829T Fiber-seq, HPRC v1.1 (d9). Somatic 19-bp insertion chr1:9,364,848
G>GTTGTTGTTGTTGGTTTTTC (candidate 266437:1:INS:>TTGTTGTTGTTGGTTTTTC); 30 ALT, 52 REF and 2 other reads (AF 0.36). The
window covers 7 GRCh38-path nodes. Right of the site, all ALT rows stay on nodes carried by 89–90 of the 90
haplotypes. At columns 70–71, at least 92% of REF rows leave the GRCh38 path through the single-base node 266438 (34
haplotypes) and node 266439 (2 bp, 33 haplotypes), and then carry 16 bases that are not in the graph (alignment
operation I, columns 72–87). In sequence these are a few bases after the site (columns 50–68 are the insertion's
slots): the REF reads carry a different insertion nearby.
