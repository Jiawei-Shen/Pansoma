# COLO829T: 100 inspectable candidate-v4 examples

These 100 examples were selected from a 1,000-tensor run on `COLO829T_3M.sorted.gam`.
The set contains 70 SNPs, 20 insertions, and 10 deletions. Within each type,
selection favors higher coverage and distinct graph node IDs. This is an
inspection set; the source run used `--min-variants 1`, so low-support
candidates are present.

Open an image in the table or load the matching tensor from `tensors_100.npy.gz`:

```python
import gzip, numpy as np
with gzip.open("tensors_100.npy.gz", "rb") as f:
    tensors = np.load(f)  # (100, 7, 200, 100), int32
```

Rows in `summary.ndjson` and the tensor array use the table index below.
The `source_index` field points to the original 1,000-tensor shard.
The corrected 1,000-tensor generation took **245.76 seconds**, including initial
GBZ hashing/loading and cache creation, with debug row metadata enabled.
Rendering these **100 PNGs** with four workers took **38.21 seconds**.
Node discovery is excluded. Exact timings and build arguments are in `manifest.json`.

Channel 7 now counts distinct GBWT paths in **hprc-v1.1-mc-grch38.d9.gbz**.
Repeated visits and reverse-complement copies count once; this includes reference
paths. The old GFA W-record counts undercounted the paths available in this GBZ.
Among 368 candidate nodes in the 1,000-tensor run, the median is **90** and only
**one** node has count 1 (range 1–91). The first six channels are unchanged.

The [occurrence audit](occurrence_audit.json) checks all **109,800 occupied cells**
across 1,000 tensors and compares six nodes independently against `gbz-tool`:

| Node ID | Correct GBWT path count |
|---:|---:|
| 233 | 91 |
| 2753 | 86 |
| 2758 | 42 |
| 15395 | 1 |
| 38702 | 90 |
| 51182 | 24 |

The saved 100 examples also passed the independent GAM/graph auditor (coverage,
source records, graph bases, row grouping, padding, and occurrence values).
These node IDs refer specifically to the d9 graph: node 2753 is 114 bases and
node 51182 is `C`. Counts from a graph with different node sequences do not apply.

| # | Candidate | Coverage | ALT/REF/other | Image |
|---:|---|---:|---:|---|
| 0 | `2753:4:SNP:G>C` | 4 | 2/2/0 | [PNG](images/sample_000.png) |
| 1 | `2758:1:INS:>GAATGGAATGGAATGGAATGGAATAGATTATAATGGAAAAGAATAG` | 5 | 2/3/0 | [PNG](images/sample_001.png) |
| 2 | `15395:529:SNP:A>G` | 2 | 1/1/0 | [PNG](images/sample_002.png) |
| 3 | `22084:17:SNP:A>C` | 2 | 2/0/0 | [PNG](images/sample_003.png) |
| 4 | `28147:2:INS:>CAAATAAAATTATT` | 2 | 1/0/1 | [PNG](images/sample_004.png) |
| 5 | `38702:46:SNP:C>T` | 2 | 1/1/0 | [PNG](images/sample_005.png) |
| 6 | `44135:14:SNP:A>C` | 2 | 1/1/0 | [PNG](images/sample_006.png) |
| 7 | `46470:19:SNP:C>T` | 3 | 1/2/0 | [PNG](images/sample_007.png) |
| 8 | `46864:4:SNP:G>A` | 2 | 2/0/0 | [PNG](images/sample_008.png) |
| 9 | `46916:27:INS:>ATTG` | 2 | 1/1/0 | [PNG](images/sample_009.png) |
| 10 | `49313:5:SNP:C>T` | 2 | 1/1/0 | [PNG](images/sample_010.png) |
| 11 | `51180:0:DEL:G>` | 2 | 2/0/0 | [PNG](images/sample_011.png) |
| 12 | `51180:10:SNP:G>C` | 2 | 2/0/0 | [PNG](images/sample_012.png) |
| 13 | `51180:28:INS:>TGGAGGGAGATGGAGA` | 2 | 1/0/1 | [PNG](images/sample_013.png) |
| 14 | `51182:0:DEL:C>` | 2 | 2/0/0 | [PNG](images/sample_014.png) |
| 15 | `54263:179:SNP:G>A` | 3 | 1/2/0 | [PNG](images/sample_015.png) |
| 16 | `61228:108:SNP:C>T` | 2 | 2/0/0 | [PNG](images/sample_016.png) |
| 17 | `62736:149:SNP:G>T` | 2 | 1/1/0 | [PNG](images/sample_017.png) |
| 18 | `82150:104:SNP:C>T` | 2 | 2/0/0 | [PNG](images/sample_018.png) |
| 19 | `84774:255:SNP:G>A` | 2 | 1/1/0 | [PNG](images/sample_019.png) |
| 20 | `85641:60:SNP:C>T` | 2 | 1/1/0 | [PNG](images/sample_020.png) |
| 21 | `90138:28:SNP:C>T` | 2 | 2/0/0 | [PNG](images/sample_021.png) |
| 22 | `91193:8:SNP:C>A` | 3 | 1/2/0 | [PNG](images/sample_022.png) |
| 23 | `91193:10:INS:>CCCCCCCA` | 3 | 1/2/0 | [PNG](images/sample_023.png) |
| 24 | `93741:72:SNP:C>A` | 2 | 1/1/0 | [PNG](images/sample_024.png) |
| 25 | `93954:20:SNP:G>T` | 2 | 1/1/0 | [PNG](images/sample_025.png) |
| 26 | `94410:23:SNP:C>T` | 2 | 1/1/0 | [PNG](images/sample_026.png) |
| 27 | `112625:201:SNP:G>T` | 3 | 1/2/0 | [PNG](images/sample_027.png) |
| 28 | `112780:42:INS:>CACTCTTTCCCTACACGACGCTCTTCCGATCT` | 2 | 1/0/1 | [PNG](images/sample_028.png) |
| 29 | `113010:0:SNP:C>T` | 2 | 2/0/0 | [PNG](images/sample_029.png) |
| 30 | `114481:0:SNP:A>T` | 2 | 2/0/0 | [PNG](images/sample_030.png) |
| 31 | `143947:391:SNP:C>A` | 2 | 1/1/0 | [PNG](images/sample_031.png) |
| 32 | `146362:130:INS:>TTGATGGACTA` | 2 | 1/1/0 | [PNG](images/sample_032.png) |
| 33 | `149216:15:SNP:C>A` | 2 | 1/1/0 | [PNG](images/sample_033.png) |
| 34 | `155717:0:SNP:G>A` | 2 | 2/0/0 | [PNG](images/sample_034.png) |
| 35 | `157411:15:INS:>TCTACAGGCATGTGACTGGAGTTCAGACGTGTGCTCTTCNGATCT` | 2 | 1/0/1 | [PNG](images/sample_035.png) |
| 36 | `165925:89:SNP:C>A` | 2 | 1/1/0 | [PNG](images/sample_036.png) |
| 37 | `173894:26:SNP:G>T` | 2 | 1/1/0 | [PNG](images/sample_037.png) |
| 38 | `176921:8:SNP:C>T` | 2 | 1/1/0 | [PNG](images/sample_038.png) |
| 39 | `179139:0:SNP:C>T` | 2 | 2/0/0 | [PNG](images/sample_039.png) |
| 40 | `182413:44:SNP:G>A` | 2 | 2/0/0 | [PNG](images/sample_040.png) |
| 41 | `182416:107:SNP:C>G` | 5 | 5/0/0 | [PNG](images/sample_041.png) |
| 42 | `182416:124:INS:>AAGATCGGAAGAGCACACGTCTG` | 5 | 1/4/0 | [PNG](images/sample_042.png) |
| 43 | `204974:190:SNP:T>G` | 2 | 1/1/0 | [PNG](images/sample_043.png) |
| 44 | `207599:73:INS:>CGACGCTCTTCCGATCT` | 2 | 1/0/1 | [PNG](images/sample_044.png) |
| 45 | `212423:340:INS:>GGCATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT` | 2 | 1/0/1 | [PNG](images/sample_045.png) |
| 46 | `214080:45:SNP:C>A` | 2 | 1/1/0 | [PNG](images/sample_046.png) |
| 47 | `214776:28:SNP:C>A` | 2 | 1/1/0 | [PNG](images/sample_047.png) |
| 48 | `214776:31:INS:>AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATGCCTGTCAAATTAT` | 2 | 1/0/1 | [PNG](images/sample_048.png) |
| 49 | `223650:46:SNP:C>G` | 2 | 1/1/0 | [PNG](images/sample_049.png) |
| 50 | `231735:209:INS:>TCCGATCT` | 2 | 1/0/1 | [PNG](images/sample_050.png) |
| 51 | `231978:196:SNP:G>T` | 2 | 1/1/0 | [PNG](images/sample_051.png) |
| 52 | `244994:15:INS:>AGATCGGAAGAGCACACGT` | 4 | 1/2/1 | [PNG](images/sample_052.png) |
| 53 | `244994:100:SNP:A>G` | 2 | 1/1/0 | [PNG](images/sample_053.png) |
| 54 | `245969:391:SNP:C>A` | 2 | 1/1/0 | [PNG](images/sample_054.png) |
| 55 | `248124:0:SNP:T>G` | 2 | 2/0/0 | [PNG](images/sample_055.png) |
| 56 | `256956:98:INS:>CCCGTGCTTTGGTTGCAAGCGTACGTCACTTTCTCCCGCCCCTATA` | 2 | 1/1/0 | [PNG](images/sample_056.png) |
| 57 | `256958:65:INS:>CGTTGTTTCCCTCTTTTACCCCAGGTATTTGGGGG` | 2 | 1/1/0 | [PNG](images/sample_057.png) |
| 58 | `276604:72:SNP:C>A` | 3 | 1/2/0 | [PNG](images/sample_058.png) |
| 59 | `288184:44:SNP:G>T` | 2 | 2/0/0 | [PNG](images/sample_059.png) |
| 60 | `293382:239:SNP:G>T` | 2 | 2/0/0 | [PNG](images/sample_060.png) |
| 61 | `296192:33:SNP:T>G` | 2 | 2/0/0 | [PNG](images/sample_061.png) |
| 62 | `300558:3:SNP:C>A` | 3 | 1/2/0 | [PNG](images/sample_062.png) |
| 63 | `313660:119:DEL:C>` | 1 | 1/0/0 | [PNG](images/sample_063.png) |
| 64 | `314241:0:SNP:C>T` | 2 | 2/0/0 | [PNG](images/sample_064.png) |
| 65 | `328146:311:SNP:G>A` | 3 | 2/1/0 | [PNG](images/sample_065.png) |
| 66 | `330151:0:SNP:G>T` | 3 | 1/2/0 | [PNG](images/sample_066.png) |
| 67 | `334535:104:SNP:A>C` | 4 | 1/3/0 | [PNG](images/sample_067.png) |
| 68 | `344223:161:SNP:C>A` | 3 | 1/2/0 | [PNG](images/sample_068.png) |
| 69 | `392203:150:INS:>GACGTGTGCTCTTCCGATCT` | 3 | 1/1/1 | [PNG](images/sample_069.png) |
| 70 | `392203:209:SNP:A>C` | 3 | 1/2/0 | [PNG](images/sample_070.png) |
| 71 | `419835:374:DEL:T>` | 1 | 1/0/0 | [PNG](images/sample_071.png) |
| 72 | `467966:491:SNP:C>T` | 2 | 2/0/0 | [PNG](images/sample_072.png) |
| 73 | `471149:0:SNP:C>T` | 2 | 2/0/0 | [PNG](images/sample_073.png) |
| 74 | `479657:0:SNP:C>T` | 2 | 2/0/0 | [PNG](images/sample_074.png) |
| 75 | `483686:0:DEL:TTTCT>` | 2 | 2/0/0 | [PNG](images/sample_075.png) |
| 76 | `486268:179:SNP:A>G` | 2 | 2/0/0 | [PNG](images/sample_076.png) |
| 77 | `493284:395:SNP:G>T` | 3 | 1/2/0 | [PNG](images/sample_077.png) |
| 78 | `524006:1:DEL:CAA>` | 2 | 2/0/0 | [PNG](images/sample_078.png) |
| 79 | `524007:0:DEL:A>` | 2 | 2/0/0 | [PNG](images/sample_079.png) |
| 80 | `532106:202:SNP:C>T` | 2 | 2/0/0 | [PNG](images/sample_080.png) |
| 81 | `539049:0:DEL:AGA>` | 2 | 1/0/1 | [PNG](images/sample_081.png) |
| 82 | `539105:71:DEL:TTT>` | 2 | 2/0/0 | [PNG](images/sample_082.png) |
| 83 | `549033:70:INS:>ATCGGAAGAGC` | 3 | 1/2/0 | [PNG](images/sample_083.png) |
| 84 | `551344:367:SNP:C>A` | 3 | 1/2/0 | [PNG](images/sample_084.png) |
| 85 | `559168:0:SNP:G>A` | 3 | 3/0/0 | [PNG](images/sample_085.png) |
| 86 | `579428:0:SNP:C>T` | 2 | 2/0/0 | [PNG](images/sample_086.png) |
| 87 | `580942:0:SNP:A>T` | 2 | 2/0/0 | [PNG](images/sample_087.png) |
| 88 | `589725:119:SNP:A>G` | 2 | 2/0/0 | [PNG](images/sample_088.png) |
| 89 | `591464:122:SNP:C>A` | 2 | 2/0/0 | [PNG](images/sample_089.png) |
| 90 | `595006:90:SNP:G>C` | 3 | 1/1/1 | [PNG](images/sample_090.png) |
| 91 | `597319:153:SNP:C>T` | 3 | 1/2/0 | [PNG](images/sample_091.png) |
| 92 | `597657:107:SNP:T>G` | 2 | 2/0/0 | [PNG](images/sample_092.png) |
| 93 | `611744:18:INS:>ACAGGCATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT` | 3 | 1/1/1 | [PNG](images/sample_093.png) |
| 94 | `620919:649:INS:>AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAAATTTC` | 3 | 1/1/1 | [PNG](images/sample_094.png) |
| 95 | `652369:2:DEL:AC>` | 1 | 1/0/0 | [PNG](images/sample_095.png) |
| 96 | `654773:181:SNP:C>T` | 2 | 2/0/0 | [PNG](images/sample_096.png) |
| 97 | `656792:350:SNP:G>T` | 2 | 2/0/0 | [PNG](images/sample_097.png) |
| 98 | `663778:522:SNP:G>A` | 2 | 2/0/0 | [PNG](images/sample_098.png) |
| 99 | `685149:443:SNP:C>T` | 3 | 1/2/0 | [PNG](images/sample_099.png) |
