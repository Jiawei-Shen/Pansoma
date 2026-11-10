# Supplementary Fig. R1.m2, INDEL examples

Two sets: `fig_R1.m2_INDEL_01`–`10`, where the window spans several nodes and the ALT and REF reads take different
branches; and `fig_R1.m2_INDEL_span_01`–`08`, where the deleted bases themselves lie on two to six nodes (below).

Ten chr1 somatic truth INDEL tensors whose 101-column window spans several graph nodes, as INDEL companions of the SNV
example in `../fig_R1.m2_tensor_across_nodes.png`. Each is an exact match to a truth INDEL (label reason
`representative_allele_is_somatic_truth`, FILTER PASS, inside the benchmark BED). In each, the ALT and REF reads pass
through different nodes somewhere in the window.

Panels are drawn in the order of channels 0, 5, 6, 1, 2, 3, 4, 7: read bases, graph-reference bases (channel 5), path
count (channel 6), then base quality, site allele, mapping quality, alignment operation and strand. Rows are reads,
grouped by the allele they carry at the site (A1 = ALT, then REF, then OTHER). Columns are the 101 positions of the
window; the site starts at column 50.

## Files

| File | What |
|---|---|
| `make_figures.sh` | runs the four steps below (`SKIP_SELECT=1` reuses `candidates.tsv`) |
| `select_examples.py` → `candidates.tsv` | every exact, PASS, in-BED somatic truth INDEL tensor on chr1 of the twelve HPRC sets (HG008T / COLO829T × Illumina, PacBio HiFi or Fiber-seq, ONT × HPRC v1.1 d9 / v2.1 d46) with >= 5 ALT and REF rows: 1,661 tensors |
| `choose.py` → `picks.tsv` | the ten examples (rules below) |
| `annotate.py` → `annotations.md` | caption facts: GRCh38-path nodes under the window and, per branch column, the channel-6 codes of the ALT and REF rows decoded to haplotype counts and matched to nodes |
| `span_examples.py` → `span_picks.tsv`, `span_nodes.md` | the truth deletions whose allele spans >= 2 nodes, with each allele node's length and haplotype count |
| `fig_R1.m2_INDEL_NN.png`, `fig_R1.m2_INDEL_span_NN.png` | the figures (`scripts/visualize_tensor.py --panel-order 0,5,6,1,2,3,4,7`) |

Selection (`choose.py`): >= 10 ALT and >= 10 REF rows; the window over >= 3 GRCh38-path nodes; at least one column
outside the site where >= 80% of ALT rows and >= 80% of REF rows each share one path-count code, the two codes differ by
>= 10 and neither is the other allele's code within 3 columns (different nodes, not the same nodes shifted by the
INDEL); fewer than 10% of the REF rows carry the ALT reads' sequence (below). Of the tensors meeting the row and node
thresholds, 61 had such a code difference, 27 passed the shift test and 25 also the REF-row test. Ranked by the number
of such columns, then GRCh38-path nodes under the window, then ALT rows; one example per truth INDEL; first the best of
each set, then a second from the sets with more examples. All ten are HG008T (two from each of the six sets except
HG008T ONT v1.1 and HG008T PacBio HiFi v1.1, one each); no COLO829T tensor passed.

REF-row test. The site allele of a read is decided at the site, so a read that carries the same INDEL but is aligned
another way is counted as REF. The test compares each REF row's read bases from column 0 (gaps dropped) with the ALT
rows' most common sequence; a REF row matches when one is a prefix of the other over at least the site, the INDEL and
10 more bases. It removed the only COLO829T candidates, the 19-bp insertion chr1:9,364,848 G>GTTGTTGTTGTTGGTTTTTC in
Fiber-seq and ONT v1.1: 43 of the 52 Fiber-seq REF rows (83%; ONT 66%) carry the same duplicated sequence as the ALT
rows, but were aligned through another branch (nodes 266438 and 266439, 34 and 33 of the 90 haplotypes) with a C/T
mismatch and a 16-bp insertion three bases to the right, so the tensor gives AF 0.36 where the reads give
(30 + 43) / 84 = 0.87. Below the 10% cut: of the REF rows that start at column 0, INDEL 5 has 4 of 84 written as a mismatch
plus an insertion one base later (the same sequence as the ALT rows; Illumina reads ending in the window), and INDEL 6
has 6 of 67 whose reads end inside the A run, where the two alleles cannot be told apart.

Haplotype counts are exact for HPRC v1.1 (H = 90). For HPRC v2.1 (H = 464), codes 1–50 are square-root binned, so the
counts below are the graph's own counts of the matched node (graph index sqlite), not decoded codes.

## Captions

*Supplementary Fig. R1.m2, INDEL 1.* HG008T Illumina, HPRC v1.1 (d9). Somatic insertion chr1:213,363,821 C>CT
(candidate 2543944:0:INS:>T); 32 ALT, 78 REF and 3 other reads (AF 0.28). The window covers 14 GRCh38-path nodes, most of
them 4-bp nodes of a TTCT repeat. At columns 41–43, all 32 ALT rows traverse node 2543941 (3 bp, on the GRCh38 path,
carried by 44 of the 90 haplotypes) and 85% of REF rows the competing 3-bp node 2543942 (36 haplotypes).

*Supplementary Fig. R1.m2, INDEL 2.* HG008T Illumina, HPRC v1.1 (d9). Somatic insertion chr1:231,875,647 G>GT
(candidate 3042010:13:INS:>T); 87 ALT, 78 REF and 4 other reads (AF 0.51). The window covers 5 GRCh38-path nodes. At
column 29, 21 columns left of the site, all 87 ALT rows traverse the single-base node 3042014 (on the GRCh38 path, 40
of the 90 haplotypes) and 93% of REF rows the competing single-base node 3042015 (50 haplotypes); 40 + 50 = 90.

*Supplementary Fig. R1.m2, INDEL 3.* HG008T ONT, HPRC v1.1 (d9). Somatic deletion chr1:188,508,320 GT>G (candidate
2264449:0:DEL:T>); 26 ALT, 13 REF and 7 other reads (AF 0.57). The window covers 6 GRCh38-path nodes. At column 90,
40 columns right of the site, all 26 ALT rows traverse the single-base node 2264452 (on the GRCh38 path, 31 of the 90
haplotypes) and 92% of REF rows the competing single-base node 2264453 (59 haplotypes); 31 + 59 = 90.

*Supplementary Fig. R1.m2, INDEL 4.* HG008T PacBio HiFi, HPRC v1.1 (d9). Somatic deletion chr1:218,103,815 CA>C
(candidate 3693578:0:DEL:A>); 45 ALT, 37 REF and 16 other reads (AF 0.37). The window covers 10 GRCh38-path nodes.
At columns 47–49, just left of the site, 98% of ALT rows are on the GRCh38-path node 3693574 (82 bp, all 90
haplotypes), while 84% of REF rows traverse a three-base branch off the GRCh38 path (nodes 3693576 and 3693577, 59
haplotypes).

*Supplementary Fig. R1.m2, INDEL 5.* HG008T Illumina, HPRC v2.1 (d46). Somatic insertion chr1:159,111,805 A>AC
(candidate 19704905:0:INS:>C); 19 ALT, 123 REF and 1 other reads (AF 0.13). The window covers 13 GRCh38-path nodes.
At column 49, the base before the insertion, all 19 ALT rows traverse the single-base node 19704903 (on the GRCh38 path,
103 of the 464 haplotypes) and 88% of REF rows the competing single-base node 19704904 (345 haplotypes).

*Supplementary Fig. R1.m2, INDEL 6.* HG008T Illumina, HPRC v2.1 (d46). Somatic insertion chr1:233,387,578 C>CA
(candidate 10089773:17:INS:>A); 100 ALT and 103 REF reads (AF 0.49), of which the 200-row tensor holds 99 and 101. The
window covers 9 GRCh38-path nodes. At column 53, all 99 ALT rows traverse the single-base node 10089775 (off the GRCh38
path, 153 of the 464 haplotypes) and 85% of REF rows the competing GRCh38-path node 10089774 (311 haplotypes);
153 + 311 = 464.

*Supplementary Fig. R1.m2, INDEL 7.* HG008T ONT, HPRC v2.1 (d46). Somatic insertion chr1:189,870,594 C>CA (candidate
8181742:24:INS:>A); 14 ALT, 29 REF and 3 other reads (AF 0.30). The window covers 8 GRCh38-path nodes. At column 53,
all 14 ALT rows traverse the single-base node 8181744 (off the GRCh38 path, 112 of the 464 haplotypes) and 97% of REF
rows the competing GRCh38-path node 8181743 (344 haplotypes).

*Supplementary Fig. R1.m2, INDEL 8.* HG008T ONT, HPRC v2.1 (d46). Somatic insertion chr1:191,234,698 T>TTA (candidate
8237854:0:INS:>TA) in a TA repeat; 14 ALT, 33 REF and 1 other reads (AF 0.29). The window covers 46 GRCh38-path nodes,
most of them single bases. At column 65, all 14 ALT rows traverse node 8237873 (146 of the 464 haplotypes) and 85% of
REF rows the competing node 8237871 (284 haplotypes). At columns 81–90, the ALT rows are on nodes carried by 463–464
haplotypes, while at least 76% of REF rows are still on TA-repeat branch nodes carried by about a quarter of the
haplotypes (path-count code 27, 113–121 paths).

*Supplementary Fig. R1.m2, INDEL 9.* HG008T PacBio HiFi, HPRC v2.1 (d46). Somatic insertion chr1:194,115,661 G>GCA
(candidate 8349226:45:INS:>CA), one more unit of a CA repeat; 27 ALT, 65 REF and 1 other reads (AF 0.29). The window
covers 14 GRCh38-path nodes. Each CA unit after the site is its own 2-bp node, carried by 210, 243, 299, 405 and 463 of
the 464 haplotypes (nodes 8349232–8349236). At columns 54–58, all ALT rows are on the units carried by 210–299
haplotypes and at least 86% of REF rows on nodes carried by 405–463.

*Supplementary Fig. R1.m2, INDEL 10.* HG008T PacBio HiFi, HPRC v2.1 (d46). Somatic deletion chr1:240,187,269 CA>C
(candidate 11316111:0:DEL:A>); 11 ALT and 42 REF reads (AF 0.21). The window covers 9 GRCh38-path nodes. At columns
68–72, all 11 ALT rows traverse the 5-bp node 11316117 (AAAAG, on the GRCh38 path, 230 of the 464 haplotypes), while
95% of REF rows are on the node that follows it, 11316118 (all 464 haplotypes, also starting AAAAG): the REF reads
skip node 11316117.

## Deletions through several nodes (`fig_R1.m2_INDEL_span_01`–`08`)

`span_examples.py`: of the 1,661 truth INDEL tensors, 22 have an allele on two or more nodes (the candidate's node plus
the summary's allele path); all are deletions, since inserted bases are not graph nodes. Kept: >= 5 ALT and >= 5 REF
rows and the REF-row test above, which removed three: the 8-bp TA deletion chr1:23,017,177 in PacBio HiFi and ONT v2.1
(52% and 83% of REF rows carry the ALT sequence) and the 4-bp deletion chr1:75,902,673 GTATA>G in PacBio HiFi v1.1
(44%). One tensor per truth deletion, the one with the most allele nodes, then the most rows of the smaller allele:
eight deletions, all HG008T. In each, the ALT rows show the deletion (alignment operation D) over the columns of all
these nodes, and panels 2 and 3 show the graph bases and path counts of the deleted nodes. Several of the short
deletions lie in homopolymer runs where the tensor also holds rows of other deletion lengths (A2, A3, ...; mostly ONT).

*Supplementary Fig. R1.m2, INDEL span 1.* HG008T Illumina, HPRC v2.1 (d46). Somatic 8-bp deletion chr1:23,017,177
TTATATATA>T (candidate 1811743:0) in a TA repeat; 16 ALT, 27 REF and 5 other reads (AF 0.33). The deleted bases lie on
six nodes, carried by 20, 19, 27, 51, 78 and 77 of the 464 haplotypes (1811743, 1811744, 1811746, 1811747, 1811748,
1811749; 1, 1, 2, 2, 1 and 1 bp); the path counts rise along the repeat, so few haplotypes carry its full length. The
window covers 36 GRCh38-path nodes. At this site in PacBio HiFi and ONT, most REF rows carry the ALT sequence (see
above); in Illumina, fewer than 10% of the REF rows that start at column 0 do.

*Supplementary Fig. R1.m2, INDEL span 2.* HG008T PacBio HiFi, HPRC v2.1 (d46). Somatic 47-bp deletion chr1:5,023,298
(candidate 1697501:0); 54 ALT, 54 REF and 1 other reads (AF 0.50). The deleted bases lie on four nodes: 1697501 (1 bp,
455 of the 464 haplotypes), 1697503 (13 bp, 464), 1697504 (1 bp, 462) and the first 32 bp of 1697506 (63 bp, 464).
All ALT rows skip the four, all REF rows traverse them. The window covers 8 GRCh38-path nodes. On HPRC v1.1 the same
deletion spans two nodes.

*Supplementary Fig. R1.m2, INDEL span 3.* HG008T ONT, HPRC v1.1 (d9). Somatic 3-bp deletion chr1:224,674,272 GTTT>G
(candidate 3816171:0) in a T run; 9 ALT reads, 17 with other deletion lengths of the run (A2–A4: 4, 1 and 5 bp; 7, 5
and 5 reads), 15 REF and 7 other (AF 0.19). The three deleted T's are three single-base nodes (3816171,
3816172, 3816173), carried by 25, 29 and 33 of the 90 haplotypes. The window covers 19 GRCh38-path nodes.

*Supplementary Fig. R1.m2, INDEL span 4.* HG008T ONT, HPRC v1.1 (d9). Somatic 3-bp deletion chr1:242,964,213 CTTT>C
(candidate 3140807:0) in a T run; 8 ALT reads, 12 with other deletion lengths (A2, A3: 1 and 2 bp; 7 and 5 reads), 19
REF and 4 other (AF 0.19). The deleted bases are node 3140807 (1 bp, 48 of the 90
haplotypes) and the first 2 bp of node 3140808 (TTTT, 52). The window covers 10 GRCh38-path nodes.

*Supplementary Fig. R1.m2, INDEL span 5.* HG008T ONT, HPRC v1.1 (d9). Somatic 2-bp deletion chr1:11,619,716 CTT>C
(candidate 312358:0); 5 ALT, 43 REF and 5 other reads (AF 0.09). The deleted bases are node 312358 (1 bp, 82 of the 90
haplotypes) and the first base of node 312359 (71 bp, starting with a T run, 86). The window covers 8 GRCh38-path nodes.

*Supplementary Fig. R1.m2, INDEL span 6.* HG008T ONT, HPRC v1.1 (d9). Somatic 2-bp deletion chr1:166,733,167 TAA>T
(candidate 1884972:0); 9 ALT reads, 4 with a 1-bp deletion (A2), 29 REF and 7 other (AF 0.18). The deleted bases are node 1884972 (1 bp, 84 of the 90
haplotypes) and the first base of node 1884973 (25 bp, starting with an A run, 86). The window covers 10 GRCh38-path
nodes.

*Supplementary Fig. R1.m2, INDEL span 7.* HG008T Illumina, HPRC v1.1 (d9). Somatic 2-bp deletion chr1:217,084,133
TAA>T (candidate 3673609:0); 64 ALT, 96 REF and 8 other reads (AF 0.38). The deleted bases are node 3673609 (1 bp, 83 of the 90
haplotypes) and the first base of node 3673610 (AA, 85). The window covers 4 GRCh38-path nodes.

*Supplementary Fig. R1.m2, INDEL span 8.* HG008T ONT, HPRC v1.1 (d9). Somatic 2-bp deletion chr1:74,016,690 CAA>C
(candidate 1070945:0); 5 ALT reads, 7 with other deletion lengths (A2, A3: 4 and 1 bp; 4 and 3 reads), 17 REF and 3
other (AF 0.16). The deleted bases are node 1070945 (1 bp, 85 of the 90
haplotypes) and the first base of node 1070946 (146 bp, starting with an A run, 89). The window covers 6 GRCh38-path
nodes.
