# Pansoma: response to the reviewers (draft 20261015)

The English responses of `response_point_by_point_20261015.md`, without the internal evidence. `[待补: …]` marks numbers still pending.

## Reviewer 1

### R1.1 "Pangenome-native" is not supported

**Comment（要点）。**

- The pangenome is used mainly at vg Giraffe alignment; the downstream representation is a local pileup, so the method resembles DeepVariant, DeepSomatic or ClairS.
- Neither the input nor the network models graph edges, bubbles / superbubbles, cross-node context or long-range haplotype structure.
- Replacing linear coordinates with node IDs is not enough; haplotype and path information is discarded.
- Request: substantial additional evidence, or a re-design, showing measurable benefit from pangenome-specific information beyond alignment.

**Response.**

We agree that "fully exploit the richness of pangenome data" overstated the original design and that node identifiers alone do not make a method pangenome-native; the phrase is removed and Pansoma is now described as pangenome-aware. We do not agree that the pangenome is used mainly at alignment, so that the method resembles DeepSomatic or ClairS: those callers work on reads aligned to GRCh38, whereas Pansoma generates candidates, builds tensors and calls variants on the graph, and its classifier reads haplotype information that only the graph provides. The reviewer's criterion, a measurable benefit from pangenome information beyond alignment, can be tested directly, and the revision does so.

1. *Path-following tensors.* In the first version a tensor was confined to the node carrying the candidate, so short nodes truncated the read context (R1.m2). Each tensor row now follows its read's walk through the graph: the 101-column window around a candidate extends across every node the read traverses, on reference and alternative nodes alike, and each column records the base of the node the read passes there (channel 5). The context is therefore no longer bounded by node length: on chr1, the window of 57–76% of somatic SNV candidates and 92–97% of somatic INDEL candidates reaches beyond the candidate's own node on the HPRC graphs.

2. *A graph-native channel, with ablation tests.* Channel 6 gives, for every column, the number of haplotype paths through the node the read traverses, relative to the graph's haplotype count, which only a pangenome provides. Alone, the site's path count separates somatic from germline chr1 SNV tensors with an AUC of 0.831–0.865 on the HPRC v1.1 graph (d9 filter) and 0.885–0.915 on the HPRC v2.1 graph (d46 filter), against 0.500 on a GRCh38-only graph. Ablating channel 6 at inference (training-mean or all-haplotype value) in the 12 HPRC SNV models lowers somatic average precision from 0.199–0.352 to 0.044–0.087 (HG008T) and from 0.639–0.829 to 0.329–0.520 (COLO829T); setting it to the all-haplotype count lets 3.4–7.6 times as many germline tensors pass the threshold. Retraining without channel 6 on both HPRC graphs, with the alignments unchanged, confirms this: chr1 somatic average precision falls from 0.199–0.353 to 0.065–0.089 on HG008T (COLO829T, v2.1: from 0.755–0.829 to 0.484–0.592), and SNV F1 after population filtering falls in 11 of 12 models (median 0.03, at most 0.31); on v2.1, where the models with channel 6 have three training seeds, all six models without it fall below that seed range by test-tuned best F1. Population filtering removes most of these extra germline calls, so channel 6 and the filter partly overlap. Channel 6 also suppresses somatic SNVs on off-reference nodes; we report every SNV result with and without an off-reference rescue (R3.3, R3.4).

3. *Three-reference ablation.* To measure what the pangenome itself contributes, the same reads were aligned to a GRCh38-only graph and to the frequency-filtered HPRC v1.1 and v2.1 graphs, and models were trained and tested on each. After population filtering both HPRC graphs lead on all six datasets (SNV F1 +0.028 to +0.292, means over three training seeds), mainly because most germline variants never become candidates; R1.2 gives the full comparison and its caveats.

Pansoma does not model bubbles, competing branches or long-range haplotypes explicitly; the revised text states this as a limitation.

### R1.2 Composition and size of the pangenome

**Comment（要点）。**

- No systematic analysis of how the number and diversity of haplotypes, haplotype walks and alternative paths affect node pileups and germline identification.
- The graph decides which alleles are represented, and so which variants are suppressed before the classifier.
- With HPRC2 (many more haplotypes): better population representation, but possibly more mapping ambiguity and more absorbed somatic mutations that overlap population alleles.

**Response.**

We compared three graph references with the same reads, tensor code and labels: a GRCh38-only graph, the HPRC v1.1 graph (90 haplotypes; d9 filter: non-reference nodes on at least 9 haplotypes) and the HPRC v2.1 graph (464 haplotypes; d46 filter) (R3.1, R3.m3).

*Pangenome versus GRCh38-only.* The HPRC graphs produce 84–91% fewer germline candidates (germline variants of the matched normal that become candidate tensors), because reads carrying a germline allele already in the graph follow its path without a mismatch. After population filtering both lead on all six datasets: SNV F1 +0.028 to +0.292 (benchmark BED, means over three training seeds; single-run paired block-bootstrap intervals exclude zero in all twelve pairs), and indel F1 is also higher.

*v1.1 versus v2.1.* The larger v2.1 graph produces a further 6–11% fewer germline candidates and keeps candidates for as many somatic alleles. Its gain is mainly before population filtering (SNV F1 +0.079 to +0.154, through precision); afterwards the two graphs are equal on all six datasets (three-seed means of PASS F1 differ by −0.015 to +0.021 and of test-tuned best F1 by at most 0.012; COLO829T Illumina v1.1 at epochs 5–6, R3.1), and indel F1 differs by at most 0.006.

*Mapping ambiguity.* Graph alignment with Giraffe is as accurate as linear mapping and reduces reference bias (Sirén et al., 2021). On the same 10%-filtered HPRC v1.1 graph, Giraffe-based calling made 34% fewer small-variant errors than GRCh38, with the largest gains in difficult regions (challenging medically relevant genes, regions where GRCh38 is erroneous, L1HS repeats), and mapping to GRCh38 left 6.6–8.5 times more reads unmapped in VNTRs (Liao et al., 2023).

*Absorbed somatic alleles.* The concern holds for indels. Relative to GRCh38-only, v2.1 leaves 24–42% of HG008T and COLO829T truth indels without a candidate; at 85–96% of these loci the alternative reads follow graph nodes, entirely (43–59%) or with one remaining mismatch or indel (30–52%), mostly in repeats. Population filtering removes 76–86% of these indels anyway, so the net loss is 3.5–7.0% of truth indels. Population filtering itself removes 53–63% of truth indels, and comparable panels are the standard germline filter of tumor-only callers (Mutect2, ClairS-TO, DeepSomatic), so the graph adds little to a loss that every tumor-only caller shares (R1.4).

### R1.3 Difficult regions and ONT ambiguity

**Comment（要点）。**

- 在重复、paralog、结构复杂和低 mappability 区域，pangenome 比对仍可能有歧义，ONT 这类高噪声 reads 尤其如此；haplotype 越多，竞争比对位置可能越多。
- 汇总统计（更多 perfect reads、更少 mismatch / indel / soft clip）不能说明改善发生在 pangenome 应该起作用的地方。
- ONT 的 graph-coordinate predictions 更多：是分辨率提高，还是额外的歧义？
- 需要按区域的 benchmark 和错误分析。

**Response.**

The region- and call-level analyses the reviewer asked for show that the HPRC graphs improve calling where the linear reference is hardest. Against the GRCh38-only (linear) graph with the same reads and pipeline, SNV F1 after population filtering rises in all 12 dataset–graph pairs in the GIAB union of difficult regions and in SMaHT's extreme regions (significantly in 11 of 12 in both, mostly through precision), and is never significantly lower in any difficult stratum. The graphs do add alignment ambiguity for ONT reads in extreme regions, and long-read INDEL calls in homopolymers lose precision, which we quantify below.

1. *Alignment by region* (same reads, six datasets, SMaHT region classes). In easy regions the graphs give cleaner alignments at the same locus: on the HPRC v1.1 graph (d9 filter), HG008T PacBio mismatched bases per kb fall from 1.02 to 0.16, and 99.93–99.99% of long reads keep their placement. In extreme regions they add ambiguity: on v1.1, 6.26–7.41% of long reads lose MAPQ ≥ 30 and 3.05–4.93% gain it; on the HPRC v2.1 graph (d46 filter), 13.97–17.30% lose it. For short reads, lower MAPQ is not more misplacement: of simulated HG008-N Illumina reads of known origin, v1.1 and v2.1 misplace 0.98% and 0.93% of placed reads against 1.24% on a GRCh38-only graph (R3.6). Simulated ONT reads differ: v1.1 misplaces slightly fewer (0.61% against 0.77%), but v2.1 about twice as many (1.48%), almost all in extreme regions, in line with the more frequent confident disagreements of real ONT reads on v2.1 (HG008T: 1,953 against 1,019 reads on v1.1).

2. *Calls by region* (chr1 SNVs, PASS after population filtering, benchmark regions; each HPRC graph against the GRCh38-only (linear) graph with the same reads and pipeline; paired 1-Mb block-bootstrap intervals in the Supplementary Table). On HG008T the arms are close in easy regions (F1 0.627–0.642 on GRCh38-only, 0.636–0.659 on the HPRC graphs; all six paired intervals include zero), whereas on COLO829T the HPRC graphs gain +0.053 to +0.154 there, through recall (all six exclude zero). In difficult regions (25 Mb of chr1) F1 changes by −0.005 to +0.259, significantly for COLO829T ONT on both graphs and HG008T PacBio on v1.1. With the GIAB v3.5 stratifications, SNV F1 rises in all 12 pairs in the union of difficult regions (46.8 Mb of chr1; +0.091 to +0.483) and significantly in 10–12 of 12 pairs in low-mappability, segmental-duplication and tandem-repeat regions, never significantly below. INDELs gain outside difficult regions (9 of 12) and with Illumina in homopolymers and tandem repeats, but lose precision there with long reads (HG008T PacBio F1 −0.065 / −0.056); every reference scores at most 0.24 INDEL F1 in these regions. In extreme regions the GRCh38-only graph makes 351–2,210 PASS calls per dataset at precision 0.056–0.253, against 102–382 at 0.161–0.764 on the HPRC graphs; the F1 gain excludes zero in 11 of 12 pairs, the exception being COLO829T Fiber-seq on v1.1 (wide interval; COLO829T Illumina v1.1 at epoch 6, R3.1). Placed ONT SNVs lose no more recall in extreme regions (HG008T, v2.1: 0.451 against 0.341–0.374), but on COLO829T their precision there is lower (v2.1: 0.452 against 0.693–0.764).

3. *ONT graph-coordinate calls.* Some calls lie on graph nodes that are not part of GRCh38, so they have only graph coordinates and cannot be scored against the GRCh38 truth sets. ONT produces the most of these (chr1, HG008T: 1,409 / 1,261 on v1.1 / v2.1, against 116–754 for Illumina and 211–244 for PacBio). Checked against the HPRC alleles and the patients' assemblies, most are population or germline alleles (49.8–91.4% match an HPRC allele, 4.2–45.4% are germline alleles of the patient), and only 1–4 per HG008T run are supported as somatic by the tumor assembly. The extra ONT calls therefore reflect alignment ambiguity, not higher resolution.

### R1.4 Candidate generation removes somatic variants already encoded in the graph

**Comment（要点）。**

- Candidates come only from mismatches / indels against the traversed node. A somatic allele that is already an alternative path, with tumor reads aligning perfectly to it, gives no candidate, tensor or call.
- Called a "fatal methodological issue".
- Supplementary Fig. 3 wrongly describes these false negatives as GRCh38 reference-bias artefacts.
- An allele on an alternative path is not necessarily germline: the same allele can be germline in one person and somatic in another.

**Response.**

We agree that the mechanism exists: an allele on a graph path is not necessarily germline, and tumor reads that align
perfectly to it give no candidate. We have withdrawn the reference-bias reading of Supplementary Fig. 3 and quantified
the loss. We respectfully disagree, however, that it is a fatal methodological issue. It is the sensitivity cost that
every tumor-only caller pays when it uses population information to remove germline variants. What differs in Pansoma
is where part of this suppression happens: at candidate generation, rather than at population filtering.

*The absorbed alleles are population alleles.* On the HPRC v2.1 graph (d46 filter), 1,467 HG008T and 347 COLO829T truth
INDELs (and 140 and 54 truth SNVs) get no candidate on at least one platform although a graph path spells their ALT (a
strict test, so a lower bound). Population filtering of the kind used by tumor-only callers tags 96.7–98.3% of these
INDELs and all of these SNVs (R2.6). The part specific to the graph, absorbed but not tagged, is 0.6% of HG008T and 0.3%
of COLO829T truth INDELs, and 1–2 INDELs per dataset in the scored region. The HPRC v1.1 graph (d9 filter) of the
original submission is similar (1.2% and 0.5% by a read-level test).

*Linear tumor-only callers lose the same alleles at filtering.* Our filter uses the four resources ClairS-TO uses by
default; ClairS-TO's own default rule tags even more truth INDELs (HG008T: 87.1% against our 62.7%). At the alleles the
HPRC v2.1 graph absorbs on chr1 (2–108 INDELs and 3–10 SNVs per dataset, HCC1395 included), DeepSomatic-TO and ClairS-TO
emit a record for many of the INDELs (3–100%), but their own filters pass 9 of 780 allele–caller pairs, and after
population filtering one (on v1.1: 7 of 882, and none). DeepSomatic's tumor-only models also give the network
population allele frequencies as an input channel, the role our haplotype-count channel plays.

*End to end, the trade-off favours the graph.* Against a GRCh38-only graph with the same reads and code (chr1–22), the
HPRC graphs give a tensor to 0.2–0.6% fewer truth SNVs and 1.7–7.7% fewer truth INDELs that filtering would keep, and
absorption explains only a small part of this (HG008T, v2.1: 38 / 31 / 19 of the 446 / 542 / 268 untagged INDELs lost). In exchange they
produce 84–91% fewer germline candidates (R1.2). After filtering, 24–27% of truth INDELs are called on every reference,
SNV F1 is higher on both HPRC graphs for all six datasets, and INDEL F1 on all six at the validation best-F1 threshold
(at the recall-0.9 threshold the GRCh38-only graph is higher on the four long-read datasets, where every reference
scores 0.133 or less; R3.1).

Pansoma is a variant caller on standard Giraffe alignments, not an aligner: absorption is a property of aligning reads
to a pangenome, and any caller working on such alignments inherits it. No existing somatic small-variant caller works
on pangenome alignments: Mutect2, DeepSomatic and ClairS-TO take reads aligned to a linear sequence, and pangenome
studies of somatic calling project the graph alignments onto GRCh38 or a personal assembly before calling (Fu et al.,
2026). We state the mechanism and the second use of the population prior by the haplotype-count
channel (off-reference somatic SNVs: 57–79% have read support, 0–4 are called; R1.1, R3.3) as limitations of
population-based tumor-only calling, and describe candidate generation on alternative paths as future work.

### R1.5 HapMap mixture as a somatic benchmark, and data leakage

**Comment（要点）。**

- HapMap 的 "somatic-like" positives 是少数个体（minor donors）的 germline allele；Pansoma 的前提是图会压掉 population germline allele，两者冲突，这个 benchmark 可能系统性地对 pangenome caller 不利，也不是获得性低 VAF 突变的合适检验。
- 数据泄漏：HPRC v1.1 含 6 个 donor 中的 4 个（HG00438、HG02257、HG02486、HG02622）；companion paper（Fu, Xin, Miao et al., bioRxiv 2026.01.04.697580）提过这个问题，但没有排除这些 donor。

**Response.**

We agree that the HapMap mixture does not establish accuracy for detecting acquired somatic mutations: its positives are
germline alleles of the five minor donors (evaluated at VAF 1–16.5%), which Pansoma's population prior, used by the
trained model (R2.1), is designed to suppress. The donors' rare alleles are also found in CoLoRSdb, one of our
population filtering resources, far more often than those of our cancer benchmarks' normals (chr20 SNVs: 14.5–25.9%
against 4.0–7.9%).

On data leakage, the submitted HapMap results used a graph built without any of the six mixture donors, so the donors'
presence in the standard HPRC graphs did not apply to this analysis. That graph had the same Minigraph-Cactus
construction and d9 filter but was built from 40 HPRC year-1 individuals plus GRCh38 and CHM13 (82 haplotypes; none of
the six mixture individuals HG005, HG00438, HG02257, HG02486, HG02622 and HG002); the Methods now state this. Four of
the donors are in the HPRC v1.1 graph and all six in the HPRC v2.1 graph (d46 filter), as our companion study noted,
which is why the donor-free graph was used. Alleles that a donor shares with other individuals of the graph can still be
graph paths; only 1.4% of the graph-carried SNV positives on the autosomes have a tensor, and the evaluation leaves all
graph-carried truth out of both recall and precision (R2.1). This removes the graph overlap from the evaluation but does not make the
remaining germline positives equivalent to somatic mutations.

We therefore removed the HapMap mixture from the main benchmark; the HapMap-trained low-VAF configuration (up to 500
reads per tensor) is kept in the Supplement only as an exploratory option, and its behaviour on somatic variants is
assessed on the COLO829T dilutions (R2.4). We keep the mixture in the Supplement because, sequenced to about 500x, it
has a large truth set at allele fractions of 1–16.5%, below those our tumor/normal dilutions reach; it is rerun with the
current pipeline on the donor-free graph and scored without population filtering, which would remove many of its
germline positives (R2.4). No claim about somatic accuracy rests on it.

Low-VAF performance is now measured on in-silico dilutions of COLO829T with COLO829BL on three platforms (nominal purity
75, 50 and 25%). One model per platform and variant type is trained on chr2–22 of all four purities and tested on chr1
of each, so the tumor itself is not held out. Two effects bear on these results. Admixed normal reads can reintroduce
germline alleles that the tumor lost through loss of heterozygosity, placing them in the somatic allele-fraction range;
the false positives of the 75 and 50% mixtures have a germline-like spectrum. In the Illumina and ONT mixtures, tumor-
and normal-derived reads also differ in base quality, a cue that a real tumor lacks; we evaluated per-read base-quality
centring, which in the matched ONT comparison reduced the dilutions' SNV F1 by 0.014–0.028, and no dilution scores above
the undiluted sample (Supplementary block). Results by purity and allele fraction are in R2.4. These low-purity results
come from one tumor.

**Supplementary (R1.5).**

- *Dilutions.* Measured purity 0.66 / 0.40 / 0.20 for Illumina and Fiber-seq, 0.70 / 0.44 / 0.22 for ONT (nominal 75 /
  50 / 25%). Fiber-seq SNV F1 after population filtering 0.665–0.779 across purities (R2.4). Cosine of the
  false-positive spectrum to the germline spectrum 0.788 (75%) and 0.768 (50%), against 0.536 at full purity.
- *Base-quality centring, Illumina.* The dilutions' INDEL best F1 (upper bound) falls to 0.133–0.146, below the 100%
  sample (0.182), and their SNV PASS F1 from 0.804–0.840 to 0.678–0.709. These runs also changed the quality
  normalization and, for SNVs, the batch size, so not all of the decrease can be attributed to the cue.
- *Base-quality centring, ONT.* Two models that differ only in this centring (same code, hardware, seed, recipe and
  normalization): SNV PASS F1 0.733 / 0.721 / 0.691 / 0.714 without centring and 0.724 / 0.693 / 0.666 / 0.700 with it
  (100 / 75 / 50 / 25% purity; chr1, after population filtering); centring costs 0.009 at 100% purity and 0.014–0.028 in
  the dilutions.

### R1.6 Low sensitivity

**Comment（要点）。**

- F1 的提升来自 precision；Supplementary Fig. 1 和 Tables 2–3 里 COLO829T、HG008T 的 SNV recall 低于 ClairS-TO 和 DeepSomatic-TO。
- INDEL recall 只有约 4–22 %（COLO829T）和 13–19 %（HG008T）。
- tumor-only 没有 matched normal 或第二阶段来补回漏掉的 variant，所以 sensitivity 很重要，F1 不能替代它。

**Response.**

We agree that sensitivity is essential in tumor-only calling and that F1 cannot stand in for it. The revision reports
recall separately, with full precision–recall curves for every caller, and addresses it in three ways (details in the
Supplementary block below).

1. *A new pipeline with higher recall.* Tensors now follow each read's path across graph nodes and carry the number of
   haplotypes through every node, and the training sets label each candidate as somatic, germline (from the matched
   normal's truth) or other, with exact allele matches only. On chr1 after population filtering (HPRC v2.1 graph, three
   training seeds), recall rose from 0.575–0.596 to 0.700–0.734 for COLO829T SNVs and from 0.759–0.828 to 0.814–0.859
   for HG008T SNVs (three platforms each; the three-seed mean is higher on every platform).
2. *Recall at equal precision.* Under the identical protocol (R3.m2), at each comparator's post-filter precision,
   Pansoma's SNV recall on HG008T is 0.856–0.898 (v2.1, three-seed means), against 0.806–0.889 for DeepSomatic-TO and
   ClairS-TO. It is lower only than DeepSomatic-TO's on Illumina (HG008T 0.871 against 0.889; COLO829T 0.781 against
   0.835).
3. *Where recall is lost.* A per-allele funnel places each missed truth variant at its first loss. For HG008T SNVs the
   largest loss is the threshold; on COLO829T, 15–16% of truth SNVs are subclonal mutations below the 0.08 candidate
   allele-fraction floor. Most truth INDELs are population alleles: under the same population-filtering policy, the
   maximum attainable recall is 0.342 (HG008T) and 0.382 (COLO829T). After that filtering, Pansoma's INDEL recall is
   0.11–0.18 at its operating threshold and up to 0.22–0.28 along its precision–recall curve, and the comparators reach
   0.015–0.183.

**Supplementary (R1.6).**

*Recall before and after the revision* (chr1, after population filtering; original submission: the chr1 results of
the submitted models on the HPRC v1.1 graph (d9 filter); revision: PASS at the validation threshold (SNVs: recall 0.9
before population filtering; INDELs: recall 0.8 of the validation INDELs that population filtering keeps), ranges over
three platforms and three training seeds; revised INDEL models on HPRC v2.1 only):

| | HG008T SNV | COLO829T SNV | HG008T INDEL | COLO829T INDEL |
|---|---|---|---|---|
| original submission (v1.1) | 0.759–0.828 | 0.575–0.596 | 0.13–0.19 | 0.04–0.22 |
| revision, HPRC v1.1 | 0.803–0.831 | 0.675–0.723 | – | – |
| revision, HPRC v2.1 (d46 filter) | 0.814–0.859 | 0.700–0.734 | 0.111–0.177 | 0.107–0.160 |

*Where recall is lost* (chr1, benchmark regions, HPRC v2.1 graph (d46 filter); SNV models trained on chr2–22 of the
same tumor, INDEL models on chr2–22 of six tumors including it; ranges over three platforms and three training seeds;
INDEL candidates include partial allele matches):

| Recall | SNV HG008T | SNV COLO829T | INDEL HG008T | INDEL COLO829T |
|---|---|---|---|---|
| truth allele reaches a candidate | 0.948–0.954 | 0.812–0.817 | 0.611–0.757 | 0.565–0.695 |
| PASS | 0.841–0.900 | 0.716–0.751 | 0.111–0.191 | 0.107–0.176 |
| PASS, after population filtering | 0.814–0.859 | 0.700–0.734 | 0.111–0.177 | 0.107–0.160 |
| any threshold, after population filtering | 0.889–0.900 | 0.788–0.795 | 0.230–0.259 | 0.221–0.275 |
| maximum attainable under the same population-filtering policy (whole chr1 truth) | 0.913 | 0.966 | 0.342 | 0.382 |

- *Operating point.* For SNVs, PASS is the recall-0.9 validation threshold, set on held-out validation blocks before
  population filtering; for INDELs, it is the threshold at which recall of the validation INDELs that population
  filtering keeps is 0.8. Records down to p = 0.01 are kept as LowQual and population matches are tagged, not deleted, so users can
  move along the precision–recall curve.
- *Funnel.* On all three references (GRCh38-only, HPRC v1.1, HPRC v2.1 with the d46 filter), HG008T SNVs are called at
  81.8%, 82.0% and 83.2% (HPRC graphs: mean over three seeds), the threshold being the largest loss on every reference
  and seed; candidate-stage losses are mostly graph alleles and repeat absorption that filtering would tag anyway
  (R1.4). On COLO829T, 17–18% of truth SNVs get no candidate, mostly subclonal mutations below the 0.08 allele-fraction
  floor (15.1–16.1% of chr1 truth SNVs), whereas recall at allele fraction 0.05–0.2 is 0.850–0.966 once a candidate
  exists (R2.3, R2.4).
- *INDELs.* 418 of the 635 HG008T chr1 truth INDELs (81 of 131 on COLO829T) are exact population alleles at the
  filtering frequencies, which sets the maximum attainable recall above. At the operating threshold, post-filter
  recall is 0.111–0.177 (HG008T) and 0.107–0.160 (COLO829T); along the post-filter curve it reaches 0.230–0.259 and
  0.221–0.275. A GRCh38-only graph gives more HG008T truth INDELs a candidate (0.75–0.87 against 0.61–0.76 on v2.1),
  but these are almost all population alleles, which the filter removes. On HCC1395 (SEQC2), only 7.9% of truth INDELs are population-tagged and 0.76–0.85
  reach a candidate (v2.1; HG008T: 0.43–0.50). The Abstract and Discussion state the INDEL limitation (R3.10).
- *Comparators* (R3.m2). At their own PASS, DeepSomatic-TO and ClairS-TO reach HG008T SNV recall 0.816–0.911 and
  0.845–0.859 (precision 0.167–0.359), and after the same population filtering 0.806–0.889 and 0.831–0.845, against
  Pansoma's 0.814–0.859 (v2.1); on COLO829T they reach 0.555–0.835 and Pansoma 0.700–0.734, DeepSomatic-TO on Illumina
  being the highest (R2.3). For INDELs their post-filter recall is 0.019–0.141 on HG008T and 0.015–0.183 on COLO829T,
  against Pansoma's 0.111–0.177 and 0.107–0.160 (R2.3).
- *Other tumors* (HPRC v2.1; R1.7, R1.m1). A model trained on four other lines keeps post-filter SNV recall at
  0.684–0.726 on the unseen HG008T, at higher precision than HG008T's own models; models trained on a single tumor lose
  recall on another mainly through the transferred threshold (COLO829T models on HCC1395 Illumina: recall 0.157 at
  precision 0.934).

### R1.7 Training scale and tumor diversity

**Comment（要点）。**

- 约 200 M 参数的 ConvNeXt 主要用 COLO829（40–50 k truth variants）训练，数据量够不够？
- 报告确切参数量，以及每个平台、每种 variant type 的正 / 负训练样本数，并说明数据量撑得住模型容量。
- COLO829 是单一的 UV 主导的 melanoma，模型可能学到 tumor-specific 特征；DeepSomatic 用了多个细胞系并做了大量 augmentation。
- 把训练扩展到不同 mutational process 的肿瘤，并直接证明模型不是 tumor-specific。

**Response.**

We agree that the original model was unnecessarily large and that training on a single tumor was too narrow a basis
for a general tumor-only caller. We reduced the model, report every training set exactly, measured overfitting and the
effect of training-set size on one controlled dataset, added four tumor lines with different mutational processes, and
tested a multi-tumor model on a tumor excluded from its training (details in the Supplementary block).

*Model size and training data.* Network size now follows the size of the training sets. On the pangenome graphs,
single-tumor SNV models see 7,098–34,110 somatic examples and INDEL models 971–3,761; multi-tumor models see
237,577–281,853 and 10,390–15,486. SNV models therefore use a 50.4 M-parameter network and INDEL models an 8.7 M one,
against 199.6 M in the original configuration. Supplementary Table
[待补: number] lists the exact training and validation examples per class (somatic, germline, other) of every model, by
tumor, platform and variant type.

*Model selection and overfitting.* We chose the smaller networks by testing on held-out data (Supplementary). Each
model is the checkpoint with the best validation somatic F1 on held-out blocks of chr2–22 (early stopping), and these
checkpoints are not overfitted: for the twelve HPRC v2.1 models of HG008T and COLO829T, the training and validation
losses at the selected epoch are 0.141–0.304 and 0.178–0.445 for SNVs and 0.250–0.715 and 0.260–0.864 for INDELs, and
the validation loss is 1.00–1.27 times its minimum over training (Supplementary Table [待补: number]). F1 rises with
training data, with diminishing gains (HG008T PacBio SNVs): 0.574 with 10% of the training blocks, 0.597 with half and
0.607 with all.

*Tumor diversity and cross-tumor generalization.* We added two breast (HCC1395, HCC1954) and two lung adenocarcinoma
lines (H1437, H2009), three platforms each (R1.m1), whose substitution spectra differ from COLO829T's UV-dominated one
(C>G-enriched HCC1395, C>A-enriched H1437 and H2009). Cross-tumor tests confirm the reviewer's concern for single-tumor
models: COLO829T-trained models transfer poorly to HG008T (SNV F1 0.359 / 0.557 / 0.487 against 0.570 / 0.605 / 0.575
for HG008T's own models, Illumina / PacBio / ONT) and show strongly context-dependent recall (on HG008T, 0.493–0.841
for UV-context C>T against 0.016–0.129 for T>C), consistent with a tumor-specific bias. We no longer claim that a model
trained on one tumor transfers in general.

*Multi-tumor training improves generalization.* The direct test excludes the target tumor: a model trained on HCC1395
and the three CASTLE lines, without HG008T or COLO829T, reaches SNV F1 0.605 / 0.614 / 0.587 on HG008T, above HG008T's
own models. A second test is within a tumor: with chr2–22 of HG008T and COLO829T added to training, the model scores
0.755 / 0.774 / 0.737 on COLO829T chr1, close to COLO829T's own models (0.783 / 0.779 / 0.736). It shows that training
on six tumors keeps the performance on a training tumor's held-out chromosome, not generalization across tumors.

**Supplementary (R1.7).**

- *Spectra of the training tumors.* COLO829T positives are 71% C>T, 65% at dipyrimidines; HG008T positives are 41% C>T,
  19% C>A and 8.4% C>G at TCW; HCC1395 24–25% C>G; H1437 and H2009 42% and 38% C>A.
- *Training sets* (Supplementary Table [待补: number]; E21). An example is one candidate site (the candidate alleles of one
  variant type at one graph position) in one dataset. Training uses chr2–22 (HCC1395: chr2–19, 21, 22); 5% of its
  20,000-node blocks are the validation set, disjoint from training; chr1 is used in neither. Somatic examples match a
  somatic truth allele exactly (partial matches are left out and counted in the table); germline examples match a
  germline truth allele or overlap one by more than 45%. Other examples are the remaining candidates on GRCh38 nodes
  inside the benchmark region, whose allele matches neither truth set: sequencing and alignment errors, including those
  next to true variants. Not used: candidates outside the region, candidates on non-reference nodes without a somatic
  match, and, for the four added lines, INDEL candidates without a truth allele in homopolymers of 7 bp or more or in
  short tandem repeats, where their truth sets miss many INDELs. Each training example is used once per epoch, without
  resampling or on-the-fly augmentation; class imbalance is handled by loss weights. In single-tumor SNV models each
  somatic example is a different truth allele; the joint COLO829T models contain the same alleles at four purities
  (below).
- *Training and model size.* Size tests: on HG008T PacBio SNVs (HPRC v2.1, one recipe), networks of 8.7 M, 50.4 M
  (three seeds) and 199.6 M parameters reach PASS F1 0.616, 0.604–0.611 and 0.603, and the larger the network, the
  earlier its validation loss rises after its minimum; on HG008T Illumina INDELs (HPRC v1.1, one pair, same recipe and
  population filter), the 8.7 M network scores 0.243 against 0.187 for the 199.6 M configuration. Early stopping
  selects the epoch with the best validation somatic F1: single-tumor models epochs 4–12 of 20 (SNV) and 2–7 of 12
  (INDEL); multi-tumor models, with more training data, 10–19 and 6–9. The validation loss also scores the
  off-reference candidates without a truth match as other, which training does not use, so it overstates the gap to the
  training loss. Learning curve: two replicates per fraction of the training blocks; the 8.7 M and 199.6 M SNV networks
  are single runs with the 50.4 M recipe.
- *Single-tumor models across tumors* (HPRC v2.1, chr1, F1 after population filtering at the recall-0.9 validation
  threshold of the training tumor, benchmark BED; Supplementary Table [待补]). Models trained on HG008T (about 7,200
  somatic examples) reach 0.708, 0.754 and 0.630 on COLO829T (Illumina, PacBio HiFi / Fiber-seq, ONT), 0.024–0.106 below
  COLO829T's own models (paired block-bootstrap intervals exclude zero), 0.689–0.790 on HCC1395 and 0.787–0.937 on H1437
  and H2009; their true positives keep each target's spectrum (R3.5). COLO829T-trained models lose partly through the
  transferred threshold (HCC1395 Illumina: recall 0.157 at precision 0.934; test-tuned best F1 0.810) and partly through
  the UV preference above (HG008T's own models: 0.905–0.925 for UV-context C>T, 0.613–0.645 for T>C). On HCC1954 every
  model has low precision (0.22–0.37): 93–99% of the false positives are germline alleles of its matched normal, while its
  truth set shows no germline content (R1.m1). These models standardise base and mapping qualities by the training
  tumor's statistics. [待补: the same tests with base quality standardised per dataset (R3.m1)]
- *Multi-tumor models* (one run each; Supplementary Table [待补]). Trained on HCC1395 (without chr1 and chr20), HCC1954,
  H1437 and H2009, the model reaches 0.845 / 0.820 / 0.658 on HCC1395 chr1 (Illumina / PacBio / ONT; a held-out
  chromosome of a training tumor), against at most 0.790 / 0.741 / 0.689 for a single-tumor model of another line. Across
  tumors, on HG008T it has higher precision and lower recall than HG008T's own models; on COLO829T, the only UV-driven
  tumor, it reaches 0.671 / 0.595 / 0.380, through lower recall.
- *Augmentation.* Joint COLO829T models add in-silico purity dilutions (75, 50, 25%; R1.5). In the Illumina and ONT
  dilutions, tumor and normal reads differ in base quality; with this cue removed, no dilution scores above the
  undiluted sample (Illumina SNVs and INDELs, ONT SNVs; R1.5).

### R1.8 Construction of the training set

**Comment（要点）。**

- 正样本：GRCh38 上的 truth 投影到 graph 坐标并按 exact allele 匹配；truth allele 本身就是 graph path、投影有歧义、graph 和 linear 表示不同（尤其 INDEL）时，正样本怎么定？
- graph-only 位置通过相邻的 1-bp SNV node 间接得到 GRCh38 坐标。
- 被吸收的位点永远成不了正样本（Supp. Fig. 3）。
- 要明确写出标签和投影过程，包括 graph 一侧是否做 normalisation，并说明正样本集合没有被系统性扭曲。
- 负样本：什么叫 "determinable"？germline 是否显式包含？artefact 怎么定义？VAF 和 purity 怎么取样？

**Response.**

We thank the reviewer and agree that the original description of the training set was insufficient. The revised Methods
and a Supplementary Note describe it in full. Because no somatic small-variant caller works on pangenome alignments yet
(R1.4), there is no established convention for labelling candidates on a graph. We therefore made every choice explicit
and conservative, and anchored every label to the GRCh38 path, on which current somatic benchmarks are defined: somatic
positives require exact allele identity, representations that differ between graph and linear reference are excluded
rather than forced into a label, and germline variants and other candidates are explicit classes.

*Positions are defined on the GRCh38 path of the graph.* The HPRC graphs embed GRCh38 as a path. Every node that this
path traverses once, in either orientation, gives each of its bases one GRCh38 coordinate (a node it traverses more than
once is treated as off the path), so positions on the path and on the linear GRCh38 reference correspond one to one. For
an SNV, the truth position gives the node on the GRCh38 path and the offset within it; a candidate exists there when the
reads aligned to that node carry the ALT base as a mismatch, and its tensor is built from these aligned reads. An INDEL
is located in the same way, by its node and start offset on the GRCh38 path, and extended along the reads' alignment,
now also across several nodes (a deletion over several nodes is one candidate). Candidate INDELs are left-aligned along
each read's alignment in forward node coordinates (the VCF convention), and truth INDELs in repeats are converted to
node keys at every equivalent placement (up to 5,000), so a match does not depend on where in a repeat the candidate's
left-alignment stops. Reads that pass through nodes off the GRCh38 path at a site do not change this: each read row of
the tensor records the graph base and the haplotype path count of the node that read is aligned to. A tensor is somatic
only when its representative allele equals a truth allele.

*Coverage of the positive set.* Positives are lost mostly before labelling, at candidate generation: a truth allele that
is itself a graph path gives no mismatch and therefore no candidate (the absorption discussed in R1.4; on the HPRC v2.1
graph, population filtering tags 96.7–98.3% of these INDELs). The projection does not force an ambiguous representation
into a positive label: non-exact matches (allele or read-haplotype overlap above 0.45) are recorded as partial and
excluded from training. In the training region, 94.8–95.9% (HG008T) and 89.5–90.5% (COLO829T) of truth SNVs stand behind
a training example on the HPRC graphs (GRCh38-only graph: 98.5–98.7% and 90.1–91.4%); most of the COLO829T remainder, on
every graph, is low-allele-fraction SNVs below the candidate floor (allele fraction 0.08, three ALT reads). Among truth
INDELs that population filtering keeps, the HPRC graphs give a tensor to 1.7–7.7% fewer than the GRCh38-only graph
(R1.4). Supplementary Table [待补] gives the status of every truth allele.

*Graph-only positions.* The original manuscript gave such positions GRCh38 coordinates indirectly, through neighbouring
1-bp nodes. The new pipeline assigns them no GRCh38 coordinate: a candidate is identified by its node and offset in the
graph, and its tensor is built from the reads aligned there. A candidate on a node off the GRCh38 path has no GRCh38
truth allele to be compared with and is not used as a training example; its calls are reported in graph coordinates
(R3.3, R3.4).

*Three labels.* The model learns germline variants and sequencing or alignment errors as two explicit classes:

1. *Somatic*: an exact match to a somatic truth allele.
2. *Germline*: a match to the matched normal's truth set (HG008-N and COLO829BL assembly calls; DeepVariant calls with
   QUAL ≥ 10 for the SEQC2 and CASTLE normals), for INDELs also an allele overlap above 0.45, such as a one-base length
   difference in a homopolymer; 256,583–381,159 germline SNV examples per pangenome model. Such an INDEL is typically
   the germline variant with a different repeat length, so it is labelled germline rather than other, a negative either
   way, whereas a somatic positive requires exact allele identity (excluding partial somatic matches raised HG008T
   Illumina's best F1 from 0.623 to 0.665 in an earlier model version).
3. *Other*: no matching truth allele inside the confident region (mostly sequencing and alignment errors, including
   errors next to true variants).

Not used for training are tensors without a truth match outside the confident region, tensors off the GRCh38 path and
partial somatic matches; these criteria replace the term "determinable".

*VAF and purity.* Each tumor keeps its observed allele-fraction distribution (median somatic SNV allele fraction
0.53–0.54 for HG008T, 0.46–0.51 for COLO829T). To cover lower purities, we built in-silico COLO829T mixtures of 75%, 50%
and 25% tumor (three platforms; median somatic SNV allele fraction 0.40, 0.28 and 0.15 on Illumina) and trained joint
models on all four purities, evaluated per purity in R1.5.

### R1.m1 Additional benchmarks

**Comment（要点）。**

- 在 SEQC2 的 HCC1395 tumor / normal 细胞系上做 benchmark。
- 在 DeepSomatic 研究里另外五个细胞系（CASTLE）上做 benchmark。

（审稿人原文不放进 repo：GitHub 仓库是 public 的；原文在 Pansoma_Resubmit/Reviewer_Letter_Response.docx。）

**Response.**

We added four of these lines: HCC1395 (SEQC2 v1.2.1 truth) and three CASTLE lines of the DeepSomatic study, HCC1954 (breast) and H1437 and H2009 (lung adenocarcinoma), each with Illumina, PacBio HiFi and ONT data processed on the HPRC v2.1 graph (R2.2). Because every line adds three whole-genome datasets, we took three of the five CASTLE lines, covering both cancer types of the panel.

**The CASTLE lines serve as training data, not as benchmarks.** With one method for all tumors, we measured how many somatic truth alleles are already alleles of the pangenome, that is, carried by at least 46 of the 463 HPRC v2.1 haplotypes other than GRCh38 and therefore kept in the graph we use. For truth SNVs, this share is 0.15–1.9% for COLO829T, HCC1395, HG008T and HCC1954, but 4.6% for H2009 and 24.5% for H1437 (truth INDELs: 39.5% for H1437, against 3.0–19.1% for the other five). In H1437 and H2009, 99% and 91% of these alleles lie in segments where the B-lymphoblastoid normal has lost heterozygosity. At 21 such sites (Supplementary Table [X]), the tumor carries a common population allele (carried by 47–386 HPRC haplotypes) on 11–71 PacBio reads, whereas the normal shows 0–5 ALT reads among 22–95. A new mutation rarely recreates a common population allele, so these are most likely inherited alleles lost from the normal cell line, which the tumor–normal comparison labeled somatic. In training on the graph, these alleles do little harm: reads carrying them follow the graph path, so 99.1–99.7% of H1437's 18,028 such SNVs produce no candidate, and they make up at most 0.11% of the somatic SNV training examples of the three CASTLE lines (at most 2.3% for INDELs). In an evaluation, however, they would count as missed somatic variants: population filtering tags all but two of them, which alone would cap the recall of every caller on H1437 at about 75%. In addition, the CASTLE truth sets were derived from DeepSomatic calls, and these lines are in both comparators' training data. We therefore evaluate on HCC1395, whose independent SEQC2 truth shows no such excess, and on HG008T and COLO829T.

**Results** (chr1, SNV F1 after population filtering, Illumina / PacBio / ONT; R1.7). A model trained on HCC1395 (without chr1 and chr20) and the three CASTLE lines reaches 0.845 / 0.820 / 0.658 on HCC1395, against 0.848 / 0.839 / 0.760 for DeepSomatic-TO, which was also trained on HCC1395 with chr1 held out, and 0.832 / 0.814 / 0.760 for ClairS-TO. On HG008T, which this model never saw, it reaches 0.605 / 0.614 / 0.587, above DeepSomatic-TO (0.529 / 0.557 / 0.457) and ClairS-TO (0.552 / 0.591 / 0.552); on the UV-dominated COLO829T, also unseen, it reaches 0.671 / 0.595 / 0.380. Adding chr2–22 of HG008T and COLO829T (six tumors) raises HCC1395 to 0.850 / 0.825 / 0.730 and COLO829T to 0.755 / 0.774 / 0.737, close to the COLO829T-specific models (0.783 / 0.779 / 0.736). [待补: INDEL results of the multi-tumor models.]

### R1.m2 Tensor windows truncated at node boundaries

**Comment（要点）。**

- 窗口限制在当前节点内、在节点边界截断，相邻节点和跨路径的证据丢掉了。
- 短节点、单碱基节点很难处理；很多 alternative 节点比窗口还短。
- 可信的 graph-native 方法需要跨相邻节点的 path-aware tensor，并显式建模 edges、bubbles、alternative paths 和 competing path support。

（审稿人原文不放进 repo：GitHub 仓库是 public 的；原文在 Pansoma_Resubmit/Reviewer_Letter_Response.docx。）

**Response.**

The node-boundary truncation of the original Methods (lines 463–464) has been removed.

1. **Rows follow read walks.** Each row of a tensor (8 channels × up to 200 reads × 101 columns) follows its read's walk across every node it traverses, so a single-base node keeps its neighbors' context. Such nodes are common: 83–86% of the non-GRCh38 nodes visited by reads on HPRC v1.1 are single bases. On chr1, the window of 57–76% of somatic SNV and 92–97% of somatic INDEL candidates reaches beyond the candidate's own node on the HPRC graphs.

2. **Edges and competing branches.** Each column records the graph base (channel 5) and the haplotype path count (channel 6) of the node the read traverses there. Rows are grouped by allele and clustered by their events, including every node transition, so reads on the same branch are adjacent. In Supplementary Fig. R1.m2 (a chr1 somatic SNV, HG008T PacBio HiFi, HPRC v1.1), all 38 ALT reads traverse one branch of a single-base bubble 21 bases from the site, carried by 18 of the 90 haplotypes, and all 40 REF reads the competing branch, carried by 72. R1.1 gives the channel-6 ablation.

Bubbles and competing branches are not modeled as explicit objects but through the reads' walks and per-node haplotype counts; the revised text states this as a limitation, and "fully exploit the richness of pangenome data" is removed (R1.1).

**Supplementary (R1.m2).**

![A chr1 SNV tensor whose rows follow their reads across node boundaries](supplementary/R1.m2/fig_R1.m2_tensor_across_nodes.png)

*Supplementary Fig. R1.m2,* the eight channels of one chr1 SNV tensor (HG008T PacBio HiFi, HPRC v1.1 graph with the d9
filter; somatic truth SNV chr1:242,529,453 G>A, exact match, allele fraction 0.49; 38 ALT and 40 REF reads). Rows are
reads, grouped by the allele they carry at the site (A1 = ALT, then REF); columns are the 101 positions centered on the
site (column 50). Panels 1–8 show channels 0–7. Each row follows its read's alignment walk, so the window spans several
graph nodes. Panel 6 (channel 5) gives the graph base of the node the read traverses at each column, and panel 7
(channel 6) the number of haplotype paths through that node. At column 71, 21 bases from the site, the reads cross a
single-base bubble: all ALT reads traverse the branch carried by 18 of the 90 haplotypes (graph base G), and all REF
reads the GRCh38 branch carried by 72 (graph base T); every other column lies on nodes carried by all haplotypes.
Rendered with `scripts/visualize_tensor.py` (`supplementary/R1.m2/make_figure.sh`).

## Reviewer 2

### R2.1 Why include the HapMap benchmark

**Comment（要点）。**

- HapMap 的 positives 是 germline variant，它们的 donor 在 HPRC 图里，所以理论上 pangenome 比对会把它们去掉。
- 请讨论为什么纳入这个 benchmark，以及它在模拟 somatic mutation 上的局限。

**Response.**

We thank the reviewer. Four of the six mixture donors are indeed in the HPRC v1.1 graph and all six in HPRC v2.1, but our HapMap analysis did not use the public HPRC graphs. For it we built a separate graph with the same Minigraph-Cactus construction and d9 filter from 40 HPRC year-1 individuals plus GRCh38 and CHM13 (82 haplotypes), leaving out all six mixture individuals (the HG005 background and the minor donors HG00438, HG02257, HG02486, HG02622 and HG002). No mixture donor genome was included in the graph, so the direct donor overlap underlying the reviewer's concern does not apply to this analysis. Alleles that a donor shares with other individuals of the graph can still be graph paths. Reads carrying such an allele mostly align along that path without an edit: only 1.4% of the graph-carried SNV positives on the autosomes have a tensor (2% of the training SNV positives), and the evaluation leaves all graph-carried truth out of both recall and precision (R2.4).

We included the mixture to examine performance at low allele fractions. Our tumor/normal dilutions reach truth allele fractions of only about 5%, whereas this cell-line mixture, sequenced to about 500x, has a large truth set at 1–16.5%.

We agree that it has clear limits. Its positives are germline alleles of the minor donors, not acquired mutations, and mixing cell lines does not reproduce tumor mutational processes, clonal evolution or copy-number changes. Population filtering would remove many of these benchmark positives because they are population germline alleles, so we do not apply it to the mixture; the analysis measures low-VAF detection but does not assess the ability to distinguish somatic mutations from germline variation. Considering these limits, we moved it to the Supplement (R2.4); somatic performance at lower purity is assessed on COLO829T tumor/normal dilutions (R1.5, R2.4).

### R2.2 How the truth sets were generated

**Comment（要点）。**

- 简要说明每个 benchmark 的 truth 是怎么生成的。
- 特别是要支撑文中 "truth set 可能含 germline variants 和测序 artefacts" 的说法。

**Response.**

We substantially revised the pipeline, including how training sets are built (Methods, Supplementary Note).

*Three labels.* Replacing the original binary labels, training examples now carry one of three: somatic (1), germline
(2) or other (0: sequencing and alignment errors); examples that cannot be labelled reliably are excluded (−1).

*Building a training set.* Tumor reads are aligned to the graph with vg Giraffe. On nodes where more than 5% of
alignments differ from the node, alleles with ≥ 3 supporting reads and allele fraction ≥ 0.08 become candidates, and
each candidate site gives one tensor of up to 200 reads along their graph walks (R1.1). Truth alleles are projected
into node coordinates along the graph's GRCh38 path. A tensor is somatic or germline when its representative allele
(most ALT reads) equals such a truth allele, and other when it lies in regions confident in both truths and matches
none. Tensors outside these regions, on off-reference nodes without a somatic match, or matching a somatic allele only
partially are excluded (R1.8). Models are trained on chr2–22 (5% of node blocks for validation) and tested on chr1.

*Datasets.* Six tumor cell lines from four cancer types, each sequenced with short and long reads; only tumor reads are
used, the matched normal only for the truth.

- *HG008T* (pancreatic ductal adenocarcinoma): GIAB PacBio HiFi Revio (116×), ONT ultra-long R10.4.1 (54×) and Illumina
  2 × 150 bp. Somatic truth: the GIAB draft benchmark v0.2 (9,609 SNVs, 9,343 INDELs in 2.81 Gb), a pre-release draft
  derived from comparing the tumor and normal (HG008-N v6.2) assemblies. Germline truth: dipcall of the HG008-N assembly.
- *COLO829T* (melanoma): Illumina, PacBio Fiber-seq and ONT, plus in-silico mixtures with COLO829BL reads at 25–75%
  purity. Somatic truth: the SMaHT COLO829 reference set (Ha et al., 2025; 44,005 SNVs, 2,059 INDELs), calls of four
  Illumina somatic callers kept when PacBio reads show the ALT in the tumor but not the normal; scored genome-wide.
  Germline truth: dipcall of the COLO829BL assembly.
- *HCC1395* (triple-negative breast cancer): CASTLE PacBio HiFi and ONT, SEQC2 Illumina (about 108×). Somatic truth:
  SEQC2 v1.2.1 (41,072 high- and medium-confidence variants in 2.49 Gb), from 378 Illumina call sets and therefore
  limited to regions analysable with short reads.
- *HCC1954* (breast), *H1437*, *H2009* (lung adenocarcinoma): CASTLE PacBio HiFi (about 60×), ONT (about 90×) and
  Illumina (100×). Somatic truth: CASTLE (26,578, 87,410 and 168,397 variants in 2.51 Gb), made with DeepSomatic models
  and so favouring DeepSomatic, as its authors note.

For these four lines, germline truth is hybrid DeepVariant calls of the matched normal (QUAL ≥ 10), and unmatched
INDELs in homopolymers ≥ 7 bp or STRs are excluded, as the SEQC2 truth has few (15.5% of its INDELs, against 74.9% for
HG008T). HG008T and COLO829T are processed on all three references (R1.2); the other four lines, on the HPRC v2.1 graph,
serve as cross-sample tests (R1.7). The SMaHT HAPMAP6 mixture (Illumina, about 500×) is reported in the Supplement only, as the
empirical sensitivity by allele fraction (R2.1, R2.4); its positives are the minor donors' germline alleles absent from
the HG005 background (evaluated at expected VAF 1–16.5%; truth v1.6 of Kong et al.).

*Germline alleles in the truth sets.* We removed the artefact claim, for which we have no direct evidence, and state the
germline claim only where measured. 209 HG008T truth alleles equal an HG008-N germline allele; of the 110 that the graph
absorbs and the normal carries, 109 are the allele of the normal's other haplotype, and flanking phased SNPs indicate
local changes rather than regional LOH. At 55 COLO829T truth INDELs a COLO829BL haplotype carries the truth allele over
the whole repeat. These are ambiguities at a few sites, not errors of the truth sets. In contrast, about 30% of H1437
truth SNVs are probably germline: 27.8% are common population alleles, 97.5% of them where the B-lymphoblastoid normal
lost heterozygosity (chr9q, chr1:9–12 Mb), so the tumor–normal comparison called them somatic (H2009 about 5%, HCC1954
none). On chr1, our test chromosome, this is negligible (19 of 4,266 truth SNVs).

### R2.3 Explaining the lower recall; normal sample at missed sites

**Comment（要点）。**

- recall 比 DeepSomatic 低约 20 %（COLO829T）和 10 %（HG008T），用 truth 错误来解释需要展开。
- 解释 Supp. Fig. 2 的 false-negative 类别 "variant node not ranked #1" 和 "first aligned node mismatch"。
- 对落在已有 graph node 上的 false negative，检查 normal 样本的 reads 是否优先比对到 alternative path。

**Response.**

We thank the reviewer for this important point. We agree that the recall deficit needed a clearer explanation and that
graph alignment alone cannot establish errors in the benchmark truth. We have withdrawn that interpretation, reassessed
recall under a common protocol, and replaced the original false-negative categories with an allele-level analysis of
where variants are lost.

*Recall under a common protocol.* We evaluated the revised Pansoma models and the comparator call sets with a common
scoring protocol (chr1, benchmark regions, the same additional population filter for all callers; R3.m2). The ranges
below span sequencing platforms and, for Pansoma's SNV models, three training seeds; they are not confidence intervals.
On HG008T, SNV recall is comparable: 0.814–0.859 for Pansoma on the HPRC v2.1 graph, against 0.806–0.889 for
DeepSomatic-TO and 0.831–0.845 for ClairS-TO, with Pansoma's models trained on chr2–22 of the same tumor. On COLO829T,
Pansoma reaches 0.700–0.734 and the comparators 0.555–0.835, and a substantial deficit remains against DeepSomatic-TO on
Illumina (0.705–0.717 against 0.835). Among truth SNVs whose allele fraction in the truth set (its Illumina VAF) is
below 0.08, DeepSomatic-TO recovers 30%, compared with 6–7% for Pansoma, whose 0.08 candidate floor applies to the
allele fraction observed in its own reads; at 0.08 or above, DeepSomatic-TO recovers 94%, compared with 83–85% for
Pansoma, so the floor does not explain the deficit fully (low-allele-fraction detection: R2.4). INDEL recall remains low
in absolute terms (0.11–0.18 for the revised INDEL models; R3.10). Under the same protocol it is higher than the
comparators' on the four long-read datasets and lower on Illumina (HG008T: 0.13 against 0.14 for DeepSomatic-TO;
COLO829T: 0.14 against 0.18 for DeepSomatic-TO and 0.15 for ClairS-TO).

*Where recall is lost (chr1).* A per-allele funnel assigns each missed truth allele to the first stage at which it is
lost (R1.6). On HG008T the largest loss on every reference and seed is the classification threshold (on the HPRC
v2.1 graph 5.6–10.9% of truth SNVs, against 4.4–5.0% without a candidate and 2.6–4.0% first lost at filtering). On COLO829T, candidate generation is the main
limitation: 17–18% of truth SNVs get no candidate, mainly because 15.1–16.1% have a truth-set allele fraction below 0.08.

*Revised false-negative classification.* We have withdrawn the original category assignments, including "variant node
not ranked #1" and "first aligned node mismatch", together with their statistics, and replaced them with an allele-level
classification based on read evidence and the stage of candidate loss (revised Supplementary Fig. 2; class definitions
in Supplementary Table [X]). We additionally examined candidate-generation failures across chr1–22 on the HPRC v2.1
graph (six datasets). For INDELs, 81–93% of candidate-generation failures involve complete or partial absorption into
graph paths, including cases that leave a residual edit with a different representation. SNVs without a candidate lie
mostly in repeats on HG008T (52–81%), whereas on COLO829T they have too few ALT reads (19–70%) or no read with MAPQ > 10
(34–42% on Illumina and Fiber-seq).

*The normal at missed sites on graph nodes.* We examined the matched normal directly at a strict subset of these
alleles: those without a candidate on the HPRC v2.1 graph on at least one platform, although the GRCh38-only analysis
produces one, and whose ALT an HPRC v2.1 graph path spells exactly (1,607 HG008T and 401 COLO829T truth alleles,
chr1–22). First, as the reviewer suggested, we realigned the COLO829BL normal's
ONT reads with the tumor's aligner settings at the 401 COLO829T loci. At these loci, only 1.1% of evaluable REF-like
normal reads traversed an ALT-like path, compared with 95% of ALT-like normal reads. At the absorbed SNVs, 2.2% of
normal and 45% of tumor reads carry the ALT. For the repeat-associated INDELs, the high ALT-like fractions in both normal
and tumor reads made read-level classification inconclusive; we therefore relied on assembly evidence to assess normal
allele status. Second, we used the haplotype-resolved matched-normal assemblies (HG008-N v6.2; COLO829BL). Only 110
(6.8%) HG008T and 27 (6.7%) COLO829T alleles are carried by the matched normal, and 1,481 (92.2%) and 371 (92.5%) are
absent from both normal haplotypes. At the 1,351 HG008T INDEL loci absent from the normal, the tumor assembly
supports the exact benchmark allele at 1,153. At 1,347 loci, it supports a tumor sequence absent from both
normal haplotypes, although this broader category does not necessarily confirm the exact benchmark allele. The
manuscript now describes graph absorption as a cause of missed somatic alleles, not as evidence of benchmark error. Most
of these alleles would also be tagged by the same allele-level population filter on any reference (R1.4), but this does
not establish that they are germline in the patient. COLO829T has no tumor assembly independent of the sequencing data
analysed here (the one used in R3.3 is built from the same long reads); its normal carries the exact ALT at 26 of the 347 absorbed INDELs and a different allele of the same
length at another 22, which we treat as ambiguous rather than germline.

**Supplementary (R2.3).**

- *Definitions.* Absorbed on v2.1: no candidate on at least one platform, a tensor on the GRCh38-only graph, and an HPRC
  v2.1 path that spells the ALT exactly. Carried by the normal: a normal haplotype has the ALT at the locus. Ambiguous
  (COLO829T): no normal haplotype carries the ALT, but one has an allele of exactly the ALT length with other bases. The
  v1.1 and v2.1 counts use different definitions (below) and are not a comparison of the two graphs.
- *Further v2.1 evidence.* The ALT paths of the absorbed alleles run through nodes walked by a median of 138 (HG008T)
  and 110 (COLO829T) HPRC haplotypes. At the 1,351 HG008T INDEL loci absent from the normal, the tumor's Illumina reads
  show the expected allele at 97.4% of the evaluable loci (PacBio 97.5%), against 95.7% at random non-absorbed truth
  INDELs. The population
  filter tags 97.0% of the HG008T and 98.5% of the COLO829T absorbed alleles; in the chr1 benchmark regions the
  GRCh38-only graph calls at most two of the 39–93 HG008T truth INDELs per platform that v2.1 absorbs.
- *Filtering overlap (chr1).* The population filter tags 61 of the 702 HG008T truth SNVs (8.7%); these would also be
  tagged if the same allele-level rules were applied to another caller's output, and those not counted at the
  filtering stage of the funnel are already lost before it.
- *Normal reads (COLO829BL ONT).* Reads overlapping the absorbed loci (+-500 bp) in the SMaHT GRCh38 BAMs, realigned
  with Giraffe to each graph; tumor (COLO829T ONT) reads processed identically; a read is labelled by its own bases
  (REF-like / ALT-like) and by the graph path it walks at the locus. Evaluable reads span the allele with MAPQ > 10;
  reads are selected by their GRCh38 placement, so a read placed elsewhere on GRCh38 is not examined, for normal and
  tumor alike. Median 145 normal and 84 tumor reads span a locus. At the absorbed INDELs (mostly homopolymer and STR
  length changes), 42% of normal and 59% of tumor reads are ALT-like.
  On the v1.1 graph (497 loci): 1.1% of REF-like normal reads take the ALT path; at SNVs 4.0% of normal against 45.2% of
  tumor reads carry the ALT. The realigned tumor reads agree with Pansoma's own tumor alignments (per-locus ALT fraction
  within 0.1 at 256 of 258 v2.1 loci).
- *The HPRC v1.1 graph (original submission).* By the read-level definition used originally (tumor reads align perfectly
  to an existing graph path), 2,478 HG008T and 497 COLO829T truth alleles are absorbed; 110 (4.4%) and 34 (6.8%) are
  carried by the matched normal and 2,085 (84.1%) and 430 (86.5%) are absent from both normal haplotypes. At the 1,935
  HG008T INDEL loci absent from the normal, the tumor assembly carries a sequence in neither normal haplotype at 1,931
  and exactly the truth allele at 1,652, and the tumor's Illumina reads show the expected allele at 95.4%. With the v2.1
  definition above, v1.1 absorbs 2,002 HG008T and 446 COLO829T truth alleles, of which 104 (5.2%) and 26 (5.8%) are
  carried by the normal; the GRCh38-only graph calls at most one of the 38–114 HG008T truth INDELs per platform that
  v1.1 absorbs (chr1, benchmark regions).

### R2.4 Low-VAF detection and the low-VAF model

**Comment（要点）。**

- 10 % 的 read-support 阈值会限制低 purity、高异质性真实样本里的低 VAF 检测：真正的低 VAF SNV 应该怎么测？
- 专门的 low-VAF model 描述不清：它和主模型有什么不同？
- 它能不能找回低 VAF 变异？主模型在多低的 VAF 开始变差？要求做 VAF 递减的 sensitivity 分析。
- 两个模型各自在什么情况下用？

**Response.**

1. We agree that low-VAF detection is a limitation of Pansoma. A candidate needs an observed allele fraction of at least 0.08 (lowered from the submitted 0.10) and three ALT-supporting reads, and performance falls as the VAF decreases: on chr1 of COLO829T, SNV recall is 0.67–0.90 at VAF 0.10–0.20 across three platforms, falls below 0.35 at 0.05–0.10 and is zero below 0.05, where no truth SNV had a candidate; in the 25% purity Fiber-seq dilution, 29% of truth SNVs have no candidate (Supplementary Table [待补: number]; R1.5).
2. The low-VAF model uses the same network but deeper tensors, up to 500 reads instead of 200, so that deep data keep more of their read evidence, and a lower candidate floor (0.015). It is trained on the SMaHT HAPMAP6 mixture of six cell lines (about 500x), whose positives are the minor donors' alleles absent from the HG005 background, at expected VAFs of 1–16.5%. These positives are germline alleles, not somatic mutations, and population filtering would remove many of them, so the mixture cannot assess tumor-only somatic calling; we therefore moved it to the Supplement (R2.1).
3. On the mixture itself (rerun with the current pipeline's 200-read tensors), apart from a clear drop below 2%, recall does not decrease with VAF (SNV recall 0.30 at 1–2% and 0.81 / 0.91 / 0.80 at 2–5 / 5–10 / 10–16.5%, precision 0.98–0.99). A companion study on the same mixture likewise found little difference across 1–16.5% and a clear drop below 1% (Fu et al., 2026) [待核实: these two observations in that paper]; because of computational limits, we did not explore VAFs below 1%.
4. For tumor data, we recommend the multi-tumor model trained on six cancer lines (R1.7). The HapMap model is an experimental model trained on germline alleles in a cell-line mixture and is not suited to real tumor data; we provide it only in the Supplement.

### R2.5 What drives the higher precision

**Comment（要点）。**

- Is the precision gain mostly reduced germline contamination, or are there other factors?
- Beyond modelling germline variation at the site, does graph-aware calling help when several germline variants are nearby?

**Response.**

We agree that the precision gain needs a mechanism. Comparing the same reads, tensors and labels on a GRCh38-only graph
and the HPRC v1.1 and v2.1 graphs (R1.2, R3.1), fewer germline candidates are a major contributor; after population
filtering, fewer false calls in the hardest regions and near clustered germline variants also contribute, in
proportions that differ by dataset and overlap.

*Germline contamination.* At candidate generation the HPRC v1.1 and v2.1 graphs produce 84–90% and 85–91% fewer
candidates at the normal's germline variants than the GRCh38-only graph, consistent with reads carrying a
graph-represented germline allele following its path without a mismatch. Before population filtering, the HG008T models
call 117,205–144,251 germline tensors somatic on the GRCh38-only graph, against 2,528–5,147 on the HPRC graphs, so most
germline contamination is removed before the population filter rather than by it. After filtering, the false calls that
the GRCh38-only graph adds over an HPRC graph are 60–73% germline (the normal's allele, an overlapping allele, or within
10 bp of one) on HG008T PacBio and ONT, 22–30% on HG008T Illumina, and none on COLO829T, where v2.1 makes more germline
false calls and 71–94% of the GRCh38-only false calls lie in SMaHT extreme regions at low allele fraction (median
0.15–0.29).

*Difficult regions.* In SMaHT extreme regions (R1.3) the HPRC graphs raise SNV precision from 0.056–0.253 to
0.161–0.764 in all twelve comparisons (six datasets, two graphs; COLO829T Illumina on v1.1 at epochs 5–6, as in R3.1);
COLO829T Illumina makes 39 false calls there on v2.1 against 2,087 on the GRCh38-only graph. Elsewhere the gain is less
uniform: in difficult regions F1 changes by −0.005 to +0.259 (higher in eleven of twelve), and in easy regions by −0.003
to +0.020 on HG008T but +0.054 to +0.154 on COLO829T, through recall (on ONT also precision).

*Clustered germline variation.* We grouped calls by the number of the normal's germline variants within 50 bp,
excluding variants that overlap the call itself (0, 1, 2+). For HG008T SNVs the HPRC graphs raise F1 most with two or
more nearby (+0.107 to +0.204, against +0.022 to +0.087 with none), through precision; recall in that stratum is
0.256–0.308 on every reference. The reads are consistent with better local sequence agreement: at truth SNVs with two
or more germline variants nearby, reads differ from the traversed reference at 1.8–3.0 positions each on GRCh38-only
against 0.4–1.8 on the HPRC graphs (at least two differences: 53–74% against 5–30% of reads). Mean MAPQ is 58–60 on
every reference, which alone does not show equally accurate placement. Of the 125–236 GRCh38-only false calls per
HG008T platform in that stratum, 82–97% are germline alleles of the normal or lie within 10 bp of one; 86–94% are not
called on the HPRC graphs, and 41–87% are not even candidates there. COLO829T agrees (94–97% of 23–49 such false calls
removed). The benefit therefore extends beyond the germline allele at the site.

*Limits.* For indels the graphs lose where germline variants cluster: only 5.9–21.6% of the truth indels in that
stratum get an exact candidate, against 62.7–94.1% on GRCh38-only (R1.4). The false calls left after filtering are
mostly the patient's rare germline alleles, which neither the graph nor the filter contains (R2.6).

### R2.6 Panels of normals

**Comment（要点）。**

- Some PoNs were built from common benchmark datasets and include variants of the tested samples (HapMap in particular), which may inflate the reduction in germline contamination.
- Fig. 2e shows that > 20% of germline variants are not in the AF-HPRC graph: were they removed by the PoNs?
- Only COLO829 was evaluated without the PoN; add HG008T without the PoN.

**Response.**

*1. Did the filter see the test samples?* We use the four resources ClairS-TO uses by default: the GATK allele-frequency file of the gnomAD v2 genome callset, about 20,000 genomes aggregated from population and disease studies (Karczewski et al., 2020); dbSNP build 138, an archive of variants submitted by many studies (Sherry et al., 2001); the GATK 1000 Genomes panel of normals, Mutect2 calls on 1000 Genomes Project exomes (1000 Genomes Project Consortium, 2015; Benjamin et al., 2019); and CoLoRSdb v1.1.0, allele frequencies of 1,381 PacBio HiFi genomes, mostly a rare-disease cohort plus HPRC, HGSVC, Solve-RD and HudsonAlpha samples (CoLoRSdb, 2024). None of our benchmark lines (HG008, COLO829, HCC1395, HCC1954, H1437, H2009) is among the listed samples or cohorts of the 1000 Genomes panel or CoLoRSdb. gnomAD and dbSNP publish no sample lists, so we tested membership genetically: a resource built from an individual holds that individual's private alleles. Of the benchmark normals' rare germline SNVs on chr20, gnomAD holds 0.6–2.5% and CoLoRSdb 4.0–7.9%, the lowest of the 12 individuals tested. The filter therefore meets these patients' private germline variants as it would a new patient's, and the benchmark does not overstate its effect; the private germline variants it misses are the main remaining false calls (point 2). The HapMap mixtures differ (point 4).

*2. Germline variants the graph does not carry.* The graphs used here are frequency filtered. On HPRC v1.1 (d9 filter), 13% of the matched normal's germline SNVs and 25% of its germline indels are not graph alleles (benchmark regions); on HPRC v2.1 (d46 filter) we estimate 12–13% and 22–23% from the allele carriers among its haplotypes (R1.2). These variants become candidates, labelled germline in training. On HG008T chr1 (SNVs; per-tensor estimates of the final models), the model rejects 79% of these germline candidates on v1.1 and 87–88% on v2.1. Population filtering, which tags 95–97% of all germline SNV candidates on either graph, removes 87.5–88.5% (v1.1) and 74.3–76.1% (v2.1) of those called somatic, so 2.4–3.2% of the germline candidates remain as false calls. These are rare germline variants that neither the graph nor the population resources contain, and they make up 86.8–97.9% (v1.1) and 82.7–95.2% (v2.1) of the false calls left after filtering (R2.5).

*3. HG008T without filtering*, now reported for every dataset and reference. HG008T SNV, v1.1, chr1, benchmark regions (graph-coordinate calls not scored; R3.4): F1 0.182–0.191 without filtering and 0.585–0.621 with it. INDEL F1 at the recall-0.9 validation threshold is 0.011–0.101 against 0.015–0.251, and test-tuned best F1 (upper bound) 0.090–0.150 against 0.144–0.270.

*4. HapMap mixtures.* [待补: the authors will revise this paragraph.] Here the concern applies: the six HapMap donors are HPRC samples, four of them also 1000 Genomes samples, and CoLoRSdb includes 127 HPRC genomes; their rare germline SNVs are 2–5 times more often in CoLoRSdb (14.5–25.9%) than the cancer lines'. The HapMap positives were already filtered in that benchmark (R1.5).

### R2.7 Alignment statistics overstated

**Comment（要点）。**

- 第 292–294 行说 graph 上 perfect reads 更多、soft clipping 更少，对照 Supplementary Fig. 4 看有些夸大。
- 要给出 read 数或比例。

**Response.**

We agree that lines 292–294 overstated the effect. They read: "Across both datasets, graph alignment increased perfect
alignments by up to 11% in WGS data and by more than four-fold in long-read datasets, while substantially reducing
soft-clipped bases, substitutions and indel events (Supplementary Fig. 4)." Supplementary Fig. 4 is replaced by
absolute read counts and proportions for both tumors on the three platforms, aligned to four references: GRCh38 (linear,
the data providers' alignments), a GRCh38-only graph, and the HPRC v1.1 and v2.1 graphs (Supplementary Table and Figure
below). Each proportion is among the reads assigned to chr1–22 on that reference; the denominators differ slightly and
are given for every dataset and reference, so the comparison is not read by read. In summary:

- *Illumina.* Perfectly aligned reads rise by about 10 percentage points on the HPRC graphs: HG008T from 1.87 billion of
  2.74 billion reads (68.3%) on linear GRCh38 to 2.12 billion of 2.70 billion (78.5–78.7%), COLO829T from 69.9% to
  79.8–80.5%. Reads with a substitution fall from 27–30% to 18–21%, and reads with an indel from 2.9–3.2% to 0.7–1.0%.
  The GRCh38-only graph gives nearly the same rates as linear GRCh38 (68.7% against 68.3%, 70.0% against 69.9%
  perfectly aligned), which suggests that the increase mainly reflects the additional allelic sequence in the HPRC
  graphs.
- *PacBio.* Perfectly aligned reads remain rare: 0.15–0.18% on GRCh38 and 1.4–1.5% on the HPRC graphs, so the
  "four-fold" increase was on a very small base. Reads with a substitution fall from 95–96% to 66–75%, while nearly
  every read still carries an indel (99% on GRCh38, 97% on the HPRC graphs).
- *ONT.* Perfectly aligned reads are 1.4% on GRCh38 and 2.7–2.8% on the HPRC graphs for HG008T, and 0.01% on every
  reference for COLO829T. Reads with a substitution fall from 93–94% on linear GRCh38 to 87–90%; indels change
  little.
- *Soft clipping* is not substantially reduced: 2.4–4.1% of Illumina reads on every reference, and for long reads it
  depends more on the aligner than on the reference.

These measures describe agreement with the traversed reference sequence, not the accuracy of read placement, which we
test with simulated reads of known origin (R3.6).

The revised text reads: "Alignment to the HPRC graphs increased the proportion of perfectly aligned Illumina reads by
about 10 percentage points in both tumors (HG008T, 68.3% on linear GRCh38 against 78.5–78.7%; COLO829T, 69.9% against
79.8–80.5%), consistent with reads matching alternative alleles represented in the graphs. Perfectly aligned long reads
remained uncommon (at most 2.8%), although the proportion with a substitution decreased (PacBio, 95–96% against 66–75%;
ONT, 93–94% against 87–90%). The proportion of long reads with an indel changed little, and soft clipping showed no
consistent reduction across references. These statistics measure agreement with the traversed reference sequence
rather than the accuracy of read placement (Supplementary Table X)."

**Supplementary (R2.7).**

![Reads on chr1–22 on four references, by platform](supplementary/R2.7/fig_R2.7.png)

Reads aligned to chr1–22 for each dataset and reference: the share perfectly aligned (no substitution, indel or soft
clipping) and the shares with at least one substitution, indel or soft-clipped end. Rows: Illumina; PacBio (HG008T HiFi, COLO829T
Fiber-seq); ONT. Note the different y-axis ranges.

Supplementary Table (reads aligned to chr1–22; N (% of these reads)):

| dataset | reference | reads on chr1–22 | perfectly aligned | with substitution | with indel | with soft clipping |
|---|---|---:|---:|---:|---:|---:|
| HG008T Illumina | GRCh38 (linear) | 2,739,599,264 | 1,869,903,341 (68.25%) | 810,918,802 (29.60%) | 78,905,626 (2.88%) | 96,334,574 (3.52%) |
| HG008T Illumina | GRCh38 graph | 2,710,092,615 | 1,861,828,942 (68.70%) | 800,943,101 (29.55%) | 79,136,287 (2.92%) | 74,218,393 (2.74%) |
| HG008T Illumina | HPRC v1.1 graph | 2,704,033,209 | 2,121,622,695 (78.46%) | 560,922,772 (20.74%) | 21,774,798 (0.81%) | 66,387,175 (2.46%) |
| HG008T Illumina | HPRC v2.1 graph | 2,695,745,607 | 2,121,713,347 (78.71%) | 552,985,555 (20.51%) | 19,851,148 (0.74%) | 63,829,495 (2.37%) |
| COLO829T Illumina | GRCh38 (linear) | 2,296,226,777 | 1,604,645,813 (69.88%) | 627,386,769 (27.32%) | 72,684,106 (3.17%) | 95,094,133 (4.14%) |
| COLO829T Illumina | GRCh38 graph | 2,283,443,818 | 1,599,146,760 (70.03%) | 629,366,832 (27.56%) | 71,551,146 (3.13%) | 78,056,743 (3.42%) |
| COLO829T Illumina | HPRC v1.1 graph | 2,285,859,727 | 1,823,482,724 (79.77%) | 423,169,962 (18.51%) | 23,805,599 (1.04%) | 78,789,217 (3.45%) |
| COLO829T Illumina | HPRC v2.1 graph | 2,265,562,422 | 1,822,701,086 (80.45%) | 410,679,337 (18.13%) | 20,250,428 (0.89%) | 66,561,782 (2.94%) |
| HG008T PacBio HiFi | GRCh38 (linear) | 12,931,641 | 18,818 (0.15%) | 12,264,874 (94.84%) | 12,862,468 (99.47%) | 961,155 (7.43%) |
| HG008T PacBio HiFi | GRCh38 graph | 12,675,374 | 18,819 (0.15%) | 12,020,862 (94.84%) | 12,606,854 (99.46%) | 314,848 (2.48%) |
| HG008T PacBio HiFi | HPRC v1.1 graph | 12,689,479 | 179,824 (1.42%) | 8,558,774 (67.45%) | 12,377,219 (97.54%) | 228,857 (1.80%) |
| HG008T PacBio HiFi | HPRC v2.1 graph | 12,977,109 | 184,297 (1.42%) | 8,605,339 (66.31%) | 12,666,519 (97.61%) | 422,918 (3.26%) |
| COLO829T PacBio Fiber-seq | GRCh38 (linear) | 25,308,248 | 44,846 (0.18%) | 24,325,270 (96.12%) | 25,138,324 (99.33%) | 2,529,764 (10.00%) |
| COLO829T PacBio Fiber-seq | GRCh38 graph | 25,026,903 | 44,873 (0.18%) | 24,011,605 (95.94%) | 24,858,107 (99.33%) | 936,884 (3.74%) |
| COLO829T PacBio Fiber-seq | HPRC v1.1 graph | 25,070,732 | 386,300 (1.54%) | 18,732,761 (74.72%) | 24,353,478 (97.14%) | 717,920 (2.86%) |
| COLO829T PacBio Fiber-seq | HPRC v2.1 graph | 26,066,512 | 397,244 (1.52%) | 19,343,766 (74.21%) | 25,350,225 (97.25%) | 1,538,959 (5.90%) |
| HG008T ONT | GRCh38 (linear) | 3,208,600 | 43,815 (1.37%) | 3,009,593 (93.80%) | 3,084,419 (96.13%) | 1,215,595 (37.89%) |
| HG008T ONT | GRCh38 graph | 3,100,754 | 43,748 (1.41%) | 2,869,956 (92.56%) | 2,950,998 (95.17%) | 930,380 (30.00%) |
| HG008T ONT | HPRC v1.1 graph | 3,099,600 | 85,847 (2.77%) | 2,691,454 (86.83%) | 2,929,825 (94.52%) | 913,057 (29.46%) |
| HG008T ONT | HPRC v2.1 graph | 3,087,083 | 84,623 (2.74%) | 2,704,935 (87.62%) | 2,944,975 (95.40%) | 899,235 (29.13%) |
| COLO829T ONT | GRCh38 (linear) | 12,636,129 | 723 (0.01%) | 11,794,165 (93.34%) | 11,928,948 (94.40%) | 12,590,223 (99.64%) |
| COLO829T ONT | GRCh38 graph | 11,839,534 | 778 (0.01%) | 10,956,152 (92.54%) | 11,047,335 (93.31%) | 11,723,760 (99.02%) |
| COLO829T ONT | HPRC v1.1 graph | 11,877,368 | 1,121 (0.01%) | 10,662,112 (89.77%) | 11,030,919 (92.87%) | 11,760,246 (99.01%) |
| COLO829T ONT | HPRC v2.1 graph | 11,789,372 | 1,098 (0.01%) | 10,645,364 (90.30%) | 11,021,767 (93.49%) | 11,674,637 (99.03%) |

Data: `supplementary/R2.7/alignment_counts.tsv` (made by `supplementary/R2.7/make_figures.py` from
`experiments/E11a_alignment_stats/full/summary.json`).

*Methods.*
1. *Alignments.* Each dataset's tumor reads were aligned with vg Giraffe to the GRCh38-only graph and to the HPRC v1.1
   (d9) and v2.1 (d46) graphs (versions as in R3.m1). Linear GRCh38 is the data providers' alignment: for HG008T (GIAB)
   BWA-MEM 0.7.17, pbmm2 1.10.0 and minimap2 2.26; for COLO829T (SMaHT) BWA-MEM 0.7.17 in Sentieon 202308.01 (with
   duplicate marking, indel realignment and base recalibration), pbmm2 1.13.0 and minimap2 2.26.
2. *Reads counted.* Every read of the dataset, each mate of a pair separately, whose alignment is on chr1–22: on linear
   GRCh38 the primary alignment (secondary and supplementary alignments are not counted); on a graph, the GRCh38
   chromosome holding most of the read's reference-path bases. Unmapped reads, reads on other contigs and reads aligned
   only to non-reference nodes are not counted.
3. *Measures.* Perfectly aligned: aligned end to end with no substitution, indel or soft clipping. Substitution: at least
   one mismatched base (linear: NM minus inserted and deleted bases; graph: mismatches against the sequence of the
   read's graph path, so a match to a non-reference allele is not a substitution). Indel: at least one inserted or
   deleted base. Soft clipping: read bases at either end not aligned to the reference (linear: CIGAR S in the primary
   alignment; hard clipping occurs only in supplementary alignments, which are not counted; graph: read bases at an
   alignment end not aligned to the path).
4. On linear GRCh38 a long read split by the aligner keeps its remainder as a supplementary alignment, so its primary
   alignment is soft-clipped; Giraffe writes one alignment per read.

### R2.m Minor comments 1–8

**Comment（要点）。**

1. Figure 2d、2f 是 Pansoma 的输出还是 truth / 旧的 linear 分析？如果是 truth，可以移到 Supplement。
2. 正文先引用 Fig. 2f、后引用 Fig. 2e。
3. Supplementary Fig. 2 里 COLO829T 的图没有标注。
4. Supplementary Fig. 1 的 purity / VAF 标签压在点上。
5. 第 13 页 "truth overlapping" 是什么意思。
6. Figure 4a 中间的 signature 柱状图也加上 truth set（同 R3.5）。
7. Figure 4 的 panel 顺序和正文引用顺序不一致。
8. Supplementary Fig. 5 里 "candidate sites" 和 "candidates" 有什么区别。

（审稿人原文不放进 repo：GitHub 仓库是 public 的；原文在 Pansoma_Resubmit/Reviewer_Letter_Response.docx。）

**Response.**

1. Neither panel uses Pansoma output; both summarize the published benchmark truth sets, whose construction R2.2 describes. Fig. 2d shows the allele-fraction distributions of the truth SNVs and INDELs of COLO829T, HG008T and the HapMap mixture, and Fig. 2f the SBS96 spectra and COSMIC signature decomposition of the COLO829T and HG008T truth SNVs. As the reviewer suggests, both panels are now in the Supplementary Information (Supplementary Fig. [X]), with a legend stating that they describe the truth sets.

2. After panels d and f moved to the Supplement, the remaining panels of Figure 2 were relettered in the order of their first citation, and the text and legend were updated accordingly.

3. In the revised Supplementary Fig. 2 (now the classification of missed truth alleles described in R2.3), every chart is labeled with its tumor, sequencing platform and variant type.

4. In the revised Supplementary Fig. 1, the purity and allele-fraction labels are given in the legend rather than next to the points, so no label overlaps the data.

5. We apologize for the unclear term. "Truth-overlapping calls" (page 13, line 349) meant PASS calls that match a benchmark truth SNV, that is, true positives, not calls that merely share a position with one. The "shared core" in the same sentence is the subset of these true positives found by all three callers. The revised text uses "true positives" and states the shared core separately.

6. We thank the reviewer; Reviewer 3 made the same request (R3.5). Rather than adding a truth bar to the signature decomposition, we replaced the decomposition with the SBS96 spectra themselves: for each tumor and sequencing platform, the revised Figure 4a shows the 96-channel spectrum of the chr1 truth SNVs beside those of Pansoma, DeepSomatic-TO and ClairS-TO, together with the matched-normal germline spectrum and each caller's cosine similarity to the truth (Supplementary Fig. R3.5a). The truth and each caller can thus be compared channel by channel, which shows directly which substitution classes and contexts drive agreement or disagreement. This matters most for HG008T, whose heterogeneous truth spectrum no single COSMIC signature describes well (best match SBS40a, cosine 0.78).

7. The panels of the revised Figure 4 are lettered in the order in which the text first cites them, and the text and legend were updated accordingly.

8. We apologize for leaving these terms undefined. In the original Supplementary Fig. 5, a candidate site was a position with at least one passing variant, and a candidate was one passing variant at that position, so a site carrying an SNV and a deletion counted once as a site and twice as candidates. The revised figure (R3.2) uses the units of the revised pipeline: a candidate is one allele (graph node, position, REF and ALT), and a candidate site groups the candidates at one node position and variant class (SNV or INDEL), each site giving one tensor. An SNV and an INDEL at the same position are therefore two sites. The legend and Methods define both units.

## Reviewer 3

### R3.1 Reference-level ablation

**Comment（要点）。**

- End-to-end benchmarking does not isolate the pangenome: Pansoma differs from the comparators in reference, alignment, candidate generation, tensorisation, architecture and post-processing.
- Run the Pansoma pipeline with a GRCh38 linear reference (or a single-path graph) and with the full HPRC graph. Keep candidate extraction, tensors, architecture, training and filtering thresholds identical, and show how much of the gain the pangenome drives.

**Response.**

For HG008T and COLO829T, three platforms each, the same reads were aligned with vg Giraffe to a GRCh38-only graph built from the same FASTA, the HPRC v1.1 graph (90 haplotypes; d9 filter) and the HPRC v2.1 graph (464 haplotypes; d46 filter); everything downstream is identical, with training on chr2–22 and testing on chr1 of the same tumor. SNV results are three-seed means after population filtering, in the benchmark regions.

1. *Both HPRC graphs outperform GRCh38-only.* SNV F1 is 0.586–0.768 (v1.1) and 0.575–0.780 (v2.1), against 0.454–0.691, with higher precision and recall on all six datasets (F1 +0.028 to +0.292; all paired block-bootstrap intervals exclude zero, R3.9). INDEL F1 at the validation best-F1 threshold is higher in all twelve comparisons; at the recall-0.9 threshold, where all arms are poor on long reads, GRCh38-only is higher on four datasets (R3.10).

2. *v2.1 against v1.1: a gain before population filtering, little after it.* Before filtering, v2.1 is higher on all six datasets (SNV F1 +0.079 to +0.154), mostly through precision (+0.089 to +0.179). After filtering, which removes the same germline variants, the two graphs are equal within the seed range on all six datasets (means differ by −0.015 to +0.021). COLO829T Illumina needs v1.1 at epochs 5–6, chosen on chr1; at the validation-selected epoch it calls unlabelled tensors in two chr1 segments and is 0.164 lower (R3.9).

3. *The gain comes from alignment and candidate generation.* Reads carrying a germline allele already in the graph follow its path without a mismatch (HG008T PacBio on v1.1: 0.16 instead of 1.02 mismatched bases per kb; R1.3). Only 7–13% of the matched normal's germline variants then become candidates, against 76–82% on GRCh38-only, while somatic SNVs keep theirs (89–96%, against 90–99%; R3.2). Retrained without the haplotype-count channel, both HPRC graphs still beat GRCh38-only (+0.011 to +0.299; R1.1).

### R3.2 Composition of the removed candidates

**Comment（要点）。**

- Fewer candidates help only if germline variants and artefacts are removed preferentially while true somatic variants are kept.
- Compare linear and graph candidate generation for (1) somatic truth, (2) matched-normal germline variants, (3) background / non-truth candidates.
- Stratify the reduction of false positives by whether the germline allele is encoded in the HPRC graph, to connect the 79 % / 73 % overlap of normal-sample SNVs / INDELs with the graph to the precision gain.

**Response.**

We thank the reviewer for this suggestion. We performed a reference-controlled candidate-composition analysis on all six datasets, with identical reads and the same candidate-generation code on the GRCh38-only, HPRC v1.1 (d9 filter) and HPRC v2.1 (d46 filter) graphs (R3.1).

*Selective reduction.* The HPRC graphs strongly reduce germline and background candidates while largely preserving somatic SNV candidates. With the HPRC graphs, 89.2–95.6% of somatic SNVs remain candidates, against 89.6–98.6% with GRCh38-only, whereas only 7.2–13.1% of matched-normal germline variants become candidates, against 75.7–82.2%. Non-truth SNV candidates fall by 5–61% (non-truth INDEL candidates by 75–84% on Illumina, but not on long reads, where most are sequencing errors). Somatic INDELs lose more candidates (insertions 0.36–0.58 against 0.53–0.96), but 76–86% of the lost INDELs are population alleles that filtering removes on any reference (R1.2, R1.4).

*Link to graph content.* In the benchmark regions, 87% of matched-normal germline SNVs and 75% of germline INDELs are alleles of the HPRC v1.1 graph, against 0.2–2.5% of somatic SNVs and 21–25% of somatic INDELs (exact allele matches; the 79% / 73% of the original text counted overlaps across all autosomes). Reads carrying germline alleles therefore traverse an existing alternative path without a mismatch, whereas most somatic SNVs remain candidates. Stratified by graph content, the false-positive reduction is concentrated on graph-encoded germline alleles: after population filtering, HG008T SNV false positives carrying a germline allele of the normal that is in HPRC v1.1 fall from 129 on GRCh38-only to 4 on either HPRC graph (three platforms; 97% fewer), whereas those absent from the graph fall only from 2,036 to 1,794–1,944 (5–12%).

*Remaining false positives.* After population filtering, GRCh38-only retains 1.4–7.8 times more SNV false positives than HPRC v1.1, and the HPRC graphs gain +0.028 to +0.292 SNV F1 (three-seed means; R3.1). On HG008T, 82.7–97.9% of the residual HPRC false positives carry exactly the normal's private germline allele. The excess GRCh38-only false positives lie at or near the normal's germline variants (HG008T PacBio and ONT: 60.3–73.3%) or are low-allele-fraction calls in SMaHT "extreme" regions (HG008T Illumina, COLO829T: 351–2,210 PASS calls at precision 0.056–0.253, against 102–281 at 0.161–0.764 on the HPRC graphs), suggesting reads misplaced on a reference that lacks the alternative sequence (R3.6).

Candidate reduction on the pangenome graphs is therefore selective rather than indiscriminate: germline and background candidates are preferentially suppressed while somatic SNV sensitivity is largely retained. This links graph allele content directly to the precision gain, together with fewer misplaced reads, in proportions that differ between datasets.

### R3.3 Validation of graph-coordinate variants

**Comment（要点）。**

- A subset of the final predictions exists only in graph coordinates. The evidence for them (signatures, read examples) is indirect.
- Validate a representative sample orthogonally: matched-normal alignment evidence, phased haplotype assemblies, or cross-platform data.

**Response.**

We agree that signatures and read examples are insufficient. We therefore checked the graph-coordinate SNV PASS calls against tumor and phased normal assemblies, matched-normal reads and the other two platforms, on chr1 of twelve runs (HG008T and COLO829T, three platforms, HPRC v1.1 and v2.1 graphs); for runs with more than 3,000 such calls (COLO829T ONT), we evaluated a random sample of 3,000. Pansoma calls variants that cannot be represented in GRCh38 coordinates, one of its design aims, and these analyses support a small subset of them as potentially somatic, including calls on sequence absent from GRCh38. However, population representation and matched-normal support are common, limiting the evidence that these calls broadly extend somatic detection. We have qualified our claims accordingly.

1. *Assemblies.* Each call's ALT context was aligned to a tumor assembly and to both normal haplotypes (HG008T: HG008-T v3.2 and HG008-N v6.2; COLO829T: a hifiasm assembly of the tumor's own HiFi and ONT reads, which shows consistency rather than independent support, and the COLO829BL donor-specific assembly). Calls are assigned in a fixed order to population-represented sequence (the ALT context is spelled by GRCh38 or an HPRC haplotype), matched-normal support (ALT in a normal haplotype), tumor-only assembly support, and unverifiable or unsupported. Under this ordered classification, 49.8–91.4% per run fall into the population-represented category, 4.2–45.4% into the matched-normal category and 0.2–3.3% into the tumor-only category; by the assemblies alone, without the population step, 0.6–5.2% are tumor-only, since representation in the population does not by itself establish the patient's germline status. On controls, the assembly classes alone place 130–149 of 150 placed true positives per run as tumor-only and 94–146 of 150 placed false positives in the normal.
2. *Matched-normal reads* (COLO829T, ONT). Per run, 66–93% of the calls have the ALT in at least two normal reads, and 0.7–3.1% have no ALT read in a normal with at least five REF reads while the tumor reads carry it; under the same condition, 89–94% of placed true positives lack the ALT in the normal. Absence of normal ALT support at this coverage is suggestive, rather than conclusive, evidence of tumor specificity.
3. *Other platforms.* Per tumor and graph, 37–81 calls are PASS on all three platforms; 70–85% of them are population-represented and 1–4 have tumor-only support.

Most graph-coordinate calls come from the optional off-reference rescue (115–4,569 per run with it, 6–162 without). The rescue recovers additional off-reference SNVs with tumor-only assembly and read support, but most of its additional calls lack evidence of somatic status, so it stays off by default and we report every SNV result with and without it (details in a Supplementary Table). Counting all graph-coordinate calls as false positives (R3.4) is a conservative sensitivity analysis, not a determination of each call's biological status. Reliable somatic classification of off-reference variants remains a limitation; whether assembly-derived training labels improve it requires further evaluation.

### R3.4 Metrics exclude graph-coordinate calls; "about 2%"; CHM13

**Comment（要点）。**

- TP / FP / FN are computed only for calls that project to GRCh38, so precision, recall and F1 describe a subset, not the complete output.
- The graph-coordinate fraction depends on the platform and is large: above 10 %, and sometimes above 20 % for ONT (Supplementary Table 6). This contradicts the Discussion's "about 2%".
- Are graph-coordinate calls truly graph-dependent, or just not projectable to GRCh38? What fraction can be represented on T2T-CHM13, and what remains graph-specific?

**Response.**

We agree, and the evaluation now covers the complete output.

1. *"About 2%" was wrong and has been removed.* Every benchmark table now reports the share of graph-coordinate calls and
F1 both for the GRCh38-evaluable calls and for the complete output, counting each graph-coordinate PASS call as a false
positive (a worst case; this replaces Supplementary Table 6). Chr1 SNVs after population filtering, no benchmark regions,
twelve runs (two tumors, three platforms, HPRC v1.1 and v2.1 graphs):

   | | without the off-reference rescue (default) | with the optional rescue |
   |---|---|---|
   | graph-coordinate share of PASS calls | 0.5–10.1% | 8.0–68.0% (highest for ONT) |
   | complete-output F1 (GRCh38-only graph: 0.197–0.680) | 0.524–0.780 | 0.342–0.736 |
   | graph-coordinate calls with a CHM13 allele | 76.6% (663 calls) | 36.7–65.9% per run |

   Without the rescue, the default, the HPRC graphs are above the GRCh38-only graph in all twelve comparisons except
   HG008T ONT on v2.1 inside the benchmark regions (0.533 against 0.537). With the rescue, they fall below it in four of
   the twelve, all long-read (HG008T ONT on both graphs, COLO829T ONT on v2.1, COLO829T Fiber-seq on v1.1).
2. *Why the rescue exists.* Our truth sets are on GRCh38, so the labelled training examples lie almost entirely on nodes
   of the GRCh38 path, which nearly every haplotype traverses (93.6–99.7% of the somatic examples; all germline and
   non-variant ones). Through the path-count channel, the classifier therefore scores variants on nodes carried by few
   haplotypes low (R1.1). The optional rescue, off by default, re-scores tensors on nodes off GRCh38 with the path count of a node that
   all haplotypes traverse, which yields more calls on alternative paths. It recovers read-supported off-reference
   somatic SNVs that rtg cannot credit (22–66 per HG008T dataset genome-wide, 0–4 without it; R3.3), but on chr1 these
   are 2–8 per HG008T run against 210–1,410 graph-coordinate calls, and it causes the long-read losses above. We report
   every SNV result with and without it.
3. *Graph-specific or merely unprojected?* The GRCh38 nodes flanking a call's bubble give a GRCh38 allele for 2.2–19.1%
   of graph-coordinate calls; after population filtering this adds seven true positives across the twelve runs and moves
   complete-output F1 by at most 0.009. Pansoma can also write its calls on the graphs' T2T-CHM13 reference path (an
   optional CHM13 VCF; Supplementary Methods): with either linear reference, 38.1–67.8% of the graph-coordinate calls
   per run have an allele (77.4% of all such calls without the rescue), and those on neither are at most 2.9% of a run's
   output without the rescue. For indels, 36.2–48.6% per run have a CHM13 allele. A linear position does not single out
   somatic calls: with or without one, about 95% of the SNV calls are HPRC or germline alleles, and under 1% are tumor-only
   in the assemblies (R3.3). A Supplementary Table lists each call's node, class, GRCh38
   interval, CHM13 position and assembly class. Our truth sets are on GRCh38, so CHM13 positions do not change any score.

**Supplementary (R3.4).**

*Methods: CHM13 coordinates.* Each graph's CHM13 reference path is extracted like its GRCh38 path (`graph_prep
ref-path-scan --reference-sample CHM13`) and checked base for base against the FASTA. In HPRC v1.1 the CHM13 sample's
chrY is GRCh38's chrY, so chrY is excluded and the path is checked against T2T-CHM13 v2.0; in HPRC v2.1 it matches
T2T-CHM13 v2.0 with the masked chrY and the rCRS mitochondrion. `linear_vcf` writes the graph VCF on this path as a
second VCF, without population filtering or evaluation, because the truth sets and the filter resources are on GRCh38.
A call has a CHM13 allele when its node lies on the CHM13 path, and a GRCh38 allele when its node lies on the GRCh38 path
or the projection of point 3 gives one.

### R3.5 Mutational signatures with the truth set

**Comment（要点）。**

- Figure 4a 里把 truth 的 spectrum 和 Pansoma、DeepSomatic-TO、ClairS-TO 并排画出来，而不只是当作 cosine similarity 的参照；HG008T 最需要，因为它的 profile 不单一，各 caller 的 cosine 都不高。
- 在 Supplement 里给出 SBS96 profiles，用来支持 "graph 坐标的 variants 是 biological" 的说法。

**Response.**

We agree and made both changes.

1. *Truth spectrum beside each caller (Figure 4a).* For each tumor and platform, Figure 4a now shows the SBS96 spectra of the chr1 truth set and of the matched normal's germline variants beside those of Pansoma, DeepSomatic-TO and ClairS-TO (PASS SNVs after the same population filtering, benchmark regions; Supplementary Fig. R3.5a). On five of six datasets Pansoma's spectrum is the closest to the truth (cosine 0.996–0.997 on COLO829T and 0.940–0.941 on HG008T, against 0.913–0.995 and 0.878–0.943 for the comparators) and the furthest from the germline (0.39–0.41 and 0.78–0.79, against 0.46–0.62 and 0.76–0.86); on HG008T PacBio, DeepSomatic-TO is as close to the truth and as far from the germline (0.943 and 0.763, against 0.941 and 0.776). The germline row matters most for HG008T, whose truth spectrum is heterogeneous (best single COSMIC v3.4 match SBS40a, cosine 0.780) and where precision after filtering is below 0.5 for every caller: each caller adds T>C changes, the hallmark of the germline spectrum (15–21% of calls, against 9% in the truth and 33% in the germline; R3.2).

2. *SBS96 profiles of the graph-coordinate calls (Supplementary).* The Supplement gives the SBS96 profiles of Pansoma's graph-coordinate SNV calls for each tumor and platform (chr1; trinucleotide context from the graph path; Supplementary Fig. R3.5b). They do not support the original claim on their own: on five of six datasets they resemble the matched normal's germline spectrum (cosine 0.64–0.91) more than the tumor truth (0.50–0.71). Only COLO829T Fiber-seq leans towards the tumor (0.67 against 0.47; UV-type C>T 42%, against 16% in the germline and 59% in the truth). We therefore no longer use spectral concordance as evidence that graph-coordinate calls are somatic. Direct validation against the tumor and normal assemblies and the population haplotypes (R3.3) shows that most are population or germline alleles, while the off-reference rescue recovers 22–66 assembly-validated somatic SNVs per HG008T dataset that GRCh38 cannot represent.

**Supplementary (R3.5).**

![Figure 4a revision, HG008T](supplementary/R3.5/fig_R3.5a_HG008T.png)

![Figure 4a revision, COLO829T](supplementary/R3.5/fig_R3.5a_COLO829T.png)

*R3.5a,* revised Figure 4a: SBS96 spectra (percent of SNVs per channel) of the chr1 somatic truth, of the matched normal's
germline SNVs on chr1 in the same regions (HG008-N / COLO829BL, Dipcall), and of the PASS SNV calls
after population filtering of Pansoma (HPRC v2.1 graph, d46 filter), DeepSomatic-TO (1.9.0) and ClairS-TO (0.4.1 / 0.4.2),
in the benchmark regions, three platforms separately.

![SBS96 of the graph-coordinate SNV calls](supplementary/R3.5/fig_R3.5b_graph_coordinate.png)

*R3.5b,* SBS96 spectra of Pansoma's graph-coordinate SNV PASS calls on chr1 (same runs; no GRCh38 position, so neither
the population filter nor the benchmark regions apply). The trinucleotide context is read from the graph: the bases
flanking the site in the REF graph-path context spelled from each call's reads (7 of 6,025 calls without a usable
context are left out). COLO829T ONT: a random 3,000 of its 4,569 calls.

Supplementary Table (UV-type C>T: C>T at a 5' pyrimidine):

| use | tumor | platform | call set | SNVs | cosine to truth | cosine to normal germline | UV-type C>T share |
|---|---|---|---|---:|---:|---:|---:|
| Fig. 4a | HG008T | Illumina | Pansoma | 1,386 | 0.940 | 0.786 | 0.240 |
| Fig. 4a | HG008T | Illumina | DeepSomatic-TO | 1,649 | 0.921 | 0.815 | 0.216 |
| Fig. 4a | HG008T | Illumina | ClairS-TO | 1,435 | 0.912 | 0.821 | 0.248 |
| Fig. 4a | HG008T | PacBio HiFi | Pansoma | 1,232 | 0.941 | 0.776 | 0.260 |
| Fig. 4a | HG008T | PacBio HiFi | DeepSomatic-TO | 1,321 | 0.943 | 0.763 | 0.245 |
| Fig. 4a | HG008T | PacBio HiFi | ClairS-TO | 1,292 | 0.918 | 0.827 | 0.236 |
| Fig. 4a | HG008T | ONT | Pansoma | 1,387 | 0.941 | 0.786 | 0.254 |
| Fig. 4a | HG008T | ONT | DeepSomatic-TO | 1,954 | 0.878 | 0.858 | 0.236 |
| Fig. 4a | HG008T | ONT | ClairS-TO | 1,399 | 0.899 | 0.847 | 0.238 |
| Fig. 4a | COLO829T | Illumina | Pansoma | 2,013 | 0.996 | 0.394 | 0.643 |
| Fig. 4a | COLO829T | Illumina | DeepSomatic-TO | 4,779 | 0.913 | 0.596 | 0.336 |
| Fig. 4a | COLO829T | Illumina | ClairS-TO | 4,335 | 0.913 | 0.620 | 0.414 |
| Fig. 4a | COLO829T | Fiber-seq | Pansoma | 1,996 | 0.997 | 0.387 | 0.652 |
| Fig. 4a | COLO829T | Fiber-seq | DeepSomatic-TO | 1,853 | 0.982 | 0.457 | 0.477 |
| Fig. 4a | COLO829T | Fiber-seq | ClairS-TO | 2,334 | 0.995 | 0.490 | 0.539 |
| Fig. 4a | COLO829T | ONT | Pansoma | 2,147 | 0.996 | 0.405 | 0.631 |
| Fig. 4a | COLO829T | ONT | DeepSomatic-TO | 2,915 | 0.984 | 0.549 | 0.463 |
| Fig. 4a | COLO829T | ONT | ClairS-TO | 2,233 | 0.993 | 0.494 | 0.546 |
| graph-coordinate | HG008T | Illumina | Pansoma graph-coordinate | 754 | 0.604 | 0.793 | 0.142 |
| graph-coordinate | HG008T | PacBio HiFi | Pansoma graph-coordinate | 240 | 0.496 | 0.635 | 0.171 |
| graph-coordinate | HG008T | ONT | Pansoma graph-coordinate | 1,254 | 0.618 | 0.753 | 0.181 |
| graph-coordinate | COLO829T | Illumina | Pansoma graph-coordinate | 516 | 0.711 | 0.736 | 0.260 |
| graph-coordinate | COLO829T | Fiber-seq | Pansoma graph-coordinate | 257 | 0.674 | 0.468 | 0.416 |
| graph-coordinate | COLO829T | ONT | Pansoma graph-coordinate | 2,997 | 0.646 | 0.906 | 0.227 |
| reference | HG008T | - | matched-normal germline | 296,786 | 0.585 | 1.000 | 0.158 |
| reference | HG008T | - | chr1 somatic truth | 697 | 1.000 | 0.585 | 0.288 |
| reference | COLO829T | - | matched-normal germline | 298,057 | 0.428 | 1.000 | 0.158 |
| reference | COLO829T | - | chr1 somatic truth | 2,443 | 1.000 | 0.428 | 0.585 |

Data and code: `supplementary/R3.5/sbs96.tsv` (96 counts of every spectrum), `cosine.tsv`, `make_sbs96.py` (caller and
truth counts from E8, graph-coordinate contexts from E9c, germline spectra from E18).

### R3.6 More perfect alignment is not more correct alignment

**Comment（要点）。**

- 更多 perfect alignment、更少 mismatch / soft clip 不等于比对更正确：图里可匹配的序列本来就更多。
- 比较 GRCh38、T2T-CHM13 和完整 HPRC graph 上的比对质量，把 "参考更完整" 和 "多 haplotype 建模" 两个因素分开。
- 如果可行：MAPQ 分布、phased-haplotype 一致性、基于组装的 truth、segmental duplications 和其他困难区域。

**Response.**

We agree and removed the inference from "more perfect" to "more correct".

1. *Real reads.* Of all HG008T Illumina reads aligned to chr1–22, 78.5% are perfectly aligned (end to end, with no
substitution, indel or soft clipping) on the HPRC v1.1 graph (d9 filter) against 68.7% on the GRCh38-only graph (R2.7),
and in reads from 48 GRCh38 windows, 99.74% of those that become perfectly aligned keep their GRCh38 position. This is
consistent with reads matching alleles that the sample shares with HPRC haplotypes, but it does not show more correct
placement. Mapping confidence also changes: for HG008T PacBio, 95.88% of mapped reads have MAPQ 60 on the GRCh38-only
and v1.1 graphs and 83.32% on the HPRC v2.1 graph (d46 filter), a drop that is smaller for ONT and absent for Illumina;
both HPRC graphs add some ambiguity in SMaHT extreme regions and segmental duplications (R1.3).

2. *Placement with known origin.* We simulated reads from the HG008-N v6.2 diploid assembly (in no graph) on chr1 and
chr20, where it aligns one-to-one to GRCh38, and aligned them to four references: the GRCh38-only graph, a CHM13-only
graph built the same way (both HPRC graphs contain T2T-CHM13 as an unfiltered reference path), and the two HPRC graphs.
On the reads with a one-to-one truth on both linear references (99.8% of Illumina reads, 99.4% of long reads):

| | GRCh38-only | CHM13-only | HPRC v1.1 | HPRC v2.1 |
|---|---|---|---|---|
| Illumina (4,988,641 reads): perfectly aligned | 63.7% | 64.8% | 73.0% | 73.1% |
| Illumina: perfectly aligned but misplaced | 32,724 | 31,494 | 32,180 | 31,012 |
| Illumina: misplaced, % of placed | 1.19% | 1.12% | 0.95% | 0.90% |
| Illumina: misplaced at MAPQ ≥ 30 | 4,240 | 4,574 | 823 | 567 |
| HiFi-like (39,743): misplaced, % of placed | 0.47% | 0.41% | 0.31% | 0.36% |
| ONT-like (39,744): misplaced, % of placed | 0.75% | 0.68% | 0.60% | 1.45% |

For Illumina and HiFi-like reads the HPRC graphs misplace fewer reads, and "perfect but wrong" alignments do not
increase. Noisy ONT-like reads misplace slightly less on v1.1 but about twice as often on v2.1, almost all in SMaHT
extreme regions, as the confident disagreements of real ONT reads on v2.1 suggested (R1.3). MAPQ overstates confidence
on every reference: at MAPQ 60 the observed Illumina error rate is 5.0 × 10⁻⁴ on GRCh38-only and 6–7 × 10⁻⁵ on the HPRC
graphs, against the nominal 10⁻⁶. Of the 18,851 Illumina reads spanning a structural variant of their haplotype, the
GRCh38-only graph places 25.1% correctly and 63.6% wrongly and leaves 11.3% unmapped; on v1.1 and v2.1, 21.2% and 21.7%
are correct and 22.1% and 18.5% wrong, while 55.8% and 58.4% align in the graph without a GRCh38 position and cannot be
scored against the truth. Far from any HG008-N variant, the HPRC graphs misplace slightly more (0.81% against 0.70%).

3. *Conclusions and limits.* Within these shared, evaluable regions, the HPRC graphs improve placement beyond what CHM13
alone provides, consistent with a benefit from multi-haplotype representation; for noisy long reads, v2.1 places worse.
The simulated regions leave out much of the sequence where CHM13 corrects GRCh38, such as the more than 9 Mbp of falsely
collapsed or duplicated sequence behind many of its gains in variant calling (Aganezov et al., 2022), so this experiment
cannot fully assess the completeness of CHM13; that sequence is in both HPRC graphs. The hardest regions are
under-represented, and we did not compare real tumor reads with their alignment to the patient's phased assembly.
Pansoma can report its calls in CHM13 coordinates (R3.4), but our truth sets are on GRCh38, so the calling ablation
(R3.1) has no CHM13 arm.

### R3.7 Caution on reclassifying benchmark false negatives

**Comment（要点）。**

- 单靠 graph alignment 不能推翻已建立的 truth 注释。
- 把 Supp. Fig. 3 的 false negatives 重新归类，需要正交证据：matched-normal 的 alignments 或 haplotype-resolved assemblies。
- 否则就要把说法收弱。

**Response.**

We agree, and we withdrew the interpretation that these false negatives are reference-bias artefacts or benchmark
errors. We tested them with orthogonal evidence: haplotype-resolved assemblies of the matched normals (HG008-N v6.2;
COLO829BL) and, for HG008T, of the tumor (HG008-T v3.2), together with the HPRC haplotype walks (R1.4).

- *HG008T.* On the HPRC v1.1 graph (d9 filter), 2,478 truth alleles (201 SNVs, 2,277 INDELs) are lost because their reads align
  perfectly to an existing graph path. Only 110 (4.4%) are carried by the matched normal, and in 109 of them the
  tumor allele is the germline allele of the other normal haplotype. For 2,085 (84.1%), neither normal haplotype carries
  the allele, but other HPRC haplotypes walk the same nodes (1,854 with exactly this allele). At 1,652 of the 1,935 such
  INDELs the tumor assembly carries exactly the truth allele over the repeat, and at 1,931 it carries a sequence found
  in neither normal haplotype. Illumina tumor reads show the expected allele at 96.2% of these INDEL loci, as at 95.7%
  of randomly chosen truth INDELs that were not absorbed.
- *COLO829T.* Of 497 such alleles, 34 (6.8%) are the normal's germline allele, and 430 (86.5%) are absent from both
  normal haplotypes but lie on nodes walked by other HPRC haplotypes (390 with exactly this allele). Without a tumor
  assembly we cannot resolve 35 further INDELs at which a normal haplotype has the ALT length with different bases, so
  the germline-like share of the absorbed INDELs lies between 7.7% and 15.8%.

The evidence therefore points the other way from our original interpretation: most of these alleles are genuine
somatic events hidden by population alleles in the graph, not truth-set errors. Supplementary Fig. 3 is reframed as
"somatic alleles absorbed by the graph", and the cost of this absorption is quantified in R1.4. We did not
realign matched-normal reads at these sites; the matched-normal evidence comes from the assemblies, whose errors in
long homopolymers we can bound only indirectly (HG008-N: about 6–8% of loci, equally at absorbed and other truth
alleles).

### R3.8 Scope: tumor-only versus somatic mosaicism

**Comment（要点）。**

- 验证以肿瘤为中心（两个肿瘤加一个合成 admixture），没有正常组织 mosaicism 数据集（如 SMaHT tissues）。
- 二选一：要么全文一致地把范围限定在 tumor-only calling，要么在独立表征过的 mosaicism 数据集上验证。

**Response.**

We agree and have narrowed the scope of the manuscript to tumor-only somatic variant calling in cancer samples. The title, Abstract, Results and Discussion have been revised accordingly, and broader claims about somatic variation in general or normal-tissue mosaicism have been removed.

Pansoma has not been trained or validated for somatic mosaicism in normal tissues. The models are trained on tumor-derived somatic variants, predominantly moderate- to high-allele-fraction events (Methods), and the standard candidate configuration requires at least three ALT-supporting reads and an allele fraction of at least 0.08, so the much lower allele fractions typical of mosaicism in normal tissue are not evaluated. A low allele fraction is not the same as normal-tissue mosaicism: within the cancer scope, low tumor fractions are tested with in-silico purity dilutions of COLO829T (75, 50 and 25%; R2.4), and these results are not extrapolated to normal tissue.

Validation on independently characterised normal-tissue datasets such as the SMaHT tissues would require dedicated low-allele-fraction candidate settings, training labels and models; the Discussion identifies this as future work.

Within the cancer scope, we have instead expanded validation beyond HG008T and COLO829T to four further cancer cell lines, two breast (HCC1395, HCC1954) and two lung adenocarcinoma lines (H1437, H2009), three platforms each, tested with models trained on other tumors and with multi-tumor models (R1.7, R1.m1).

### R3.9 Model stability and statistical uncertainty

**Comment（要点）。**

- 结果都是点估计；用多 seed 重训、基因组区间的 block bootstrap 置信区间或染色体间方差评估稳定性。
- 特别是：多 seed 重训，报告 held-out 染色体和独立 HG008T benchmark 上 precision / recall / F1 的分布。
- 用染色体或 block 级重抽样给 Pansoma 与比较器之间差值的置信区间。
- 报告对固定 0.5 阈值的敏感性。

**Response.**

We agree that single point estimates are insufficient. The revised manuscript therefore reports variation across independent training seeds, paired genomic-block bootstrap confidence intervals for performance differences, and sensitivity to the choice of operating threshold.

1. *Variation across training seeds.* Each SNV model of the three-reference comparison was trained three times with independent random seeds, which also determine the held-out validation blocks; models were trained on chr2–22 and evaluated on chr1. Supplementary Tables R3.9a and R3.9c report chr1 precision, recall and F1 after population filtering for every seed and reference, with mean ± SD. Across the HPRC models (v1.1, d9 filter; v2.1, d46 filter; six datasets), recall SD is 0.006–0.019, precision SD 0.003–0.083 and F1 SD 0.001–0.040 (median 0.012), so single-run differences of several percentage points can fall within training variation for some datasets. Nevertheless, the HPRC graphs outperform the seed-matched GRCh38-only model in 35 of 36 comparisons (R3.1). The three-seed analysis includes the HG008T models trained on chr2–22 and evaluated on held-out chr1; the independent cross-tumor HG008T experiments currently use one training seed and are reported separately (R1.7).

2. *Paired genomic-block bootstrap.* Both members of each comparison are scored on the same 2,000 resamples of the 1-Mb blocks of chr1, giving paired 95% confidence intervals for Δprecision, Δrecall and ΔF1 (Supplementary Table R3.9b). For HPRC versus GRCh38-only, ΔF1 is +0.038 to +0.326 and excludes zero in all 12 SNV comparisons; Δprecision excludes zero in 11 and Δrecall in eight. For Pansoma (v2.1) versus DeepSomatic-TO, SNV ΔF1 is +0.042 to +0.219 and excludes zero on all six datasets. Versus ClairS-TO, ΔF1 is positive in all six comparisons (+0.012 to +0.232), although the interval excludes zero in only two. For INDELs at the recall-0.9 operating point, Pansoma trades higher recall for lower precision relative to both callers (R3.10). These intervals quantify the genomic sampling uncertainty of chr1 and are reported separately from seed-to-seed training variation.

3. *Operating point.* The revised models are trained with new labels (somatic, germline and non-variant; R1.8) instead of the original single label, so the fixed 0.5 threshold no longer applies. PASS is set on the held-out validation blocks, never on chr1, and we show complete precision–recall curves for every model.

**Supplementary (R3.9).**

*Supplementary Table R3.9a.* chr1 SNV precision, recall and F1 after population filtering (benchmark regions, PASS at the recall-0.9 validation threshold) of the in-sample models, three training seeds each; mean ± SD (sample SD) in brackets. COLO829T Illumina on v1.1 is shown at epochs 6 / 5 / 6 of seeds 0 / 1 / 2 (from epoch 7, or 6 for seed 1, these models call unlabelled tensors in chr1:121–126 and 143–144 Mb; the validation F1 does not change) and at the validation-selected epoch 8.

| tumor | platform | graph | precision: seeds 0 / 1 / 2 (mean ± SD) | recall | F1 |
|---|---|---|---|---|---|
| HG008T | Illumina | HPRC v1.1 (d9) | 0.456 / 0.459 / 0.461 (0.459 ± 0.003) | 0.831 / 0.819 / 0.803 (0.818 ± 0.014) | 0.588 / 0.588 / 0.586 (0.588 ± 0.001) |
| HG008T | Illumina | HPRC v2.1 (d46) | 0.426 / 0.451 / 0.433 (0.437 ± 0.013) | 0.852 / 0.831 / 0.814 (0.832 ± 0.019) | 0.568 / 0.585 / 0.566 (0.573 ± 0.010) |
| HG008T | PacBio | HPRC v1.1 (d9) | 0.484 / 0.422 / 0.460 (0.455 ± 0.031) | 0.828 / 0.811 / 0.829 (0.823 ± 0.010) | 0.611 / 0.555 / 0.592 (0.586 ± 0.028) |
| HG008T | PacBio | HPRC v2.1 (d46) | 0.471 / 0.487 / 0.485 (0.481 ± 0.009) | 0.838 / 0.814 / 0.814 (0.822 ± 0.014) | 0.603 / 0.609 / 0.607 (0.606 ± 0.003) |
| HG008T | ONT | HPRC v1.1 (d9) | 0.472 / 0.443 / 0.447 (0.454 ± 0.016) | 0.814 / 0.819 / 0.831 (0.821 ± 0.009) | 0.598 / 0.575 / 0.581 (0.585 ± 0.012) |
| HG008T | ONT | HPRC v2.1 (d46) | 0.433 / 0.470 / 0.481 (0.461 ± 0.025) | 0.859 / 0.848 / 0.824 (0.844 ± 0.018) | 0.576 / 0.605 / 0.607 (0.596 ± 0.017) |
| COLO829T | Illumina | HPRC v1.1 (d9) at epochs 6 / 5 / 6 | 0.832 / 0.705 / 0.782 (0.773 ± 0.064) | 0.723 / 0.705 / 0.698 (0.708 ± 0.013) | 0.773 / 0.705 / 0.738 (0.739 ± 0.034) |
| COLO829T | Illumina | HPRC v2.1 (d46) | 0.867 / 0.706 / 0.821 (0.798 ± 0.083) | 0.712 / 0.705 / 0.717 (0.712 ± 0.006) | 0.782 / 0.705 / 0.765 (0.751 ± 0.040) |
| COLO829T | Illumina | HPRC v1.1 (d9), validation-selected epoch 8 | 0.476 / 0.565 / 0.477 (0.506 ± 0.051) | 0.698 / 0.706 / 0.703 (0.702 ± 0.004) | 0.566 / 0.628 / 0.568 (0.587 ± 0.035) |
| COLO829T | fiberseq | HPRC v1.1 (d9) | 0.798 / 0.876 / 0.843 (0.839 ± 0.039) | 0.722 / 0.708 / 0.710 (0.713 ± 0.008) | 0.758 / 0.783 / 0.771 (0.771 ± 0.012) |
| COLO829T | fiberseq | HPRC v2.1 (d46) | 0.863 / 0.848 / 0.852 (0.854 ± 0.008) | 0.707 / 0.710 / 0.734 (0.717 ± 0.015) | 0.777 / 0.773 / 0.788 (0.779 ± 0.008) |
| COLO829T | ONT | HPRC v1.1 (d9) | 0.793 / 0.822 / 0.830 (0.815 ± 0.020) | 0.700 / 0.693 / 0.675 (0.689 ± 0.013) | 0.744 / 0.752 / 0.745 (0.747 ± 0.004) |
| COLO829T | ONT | HPRC v2.1 (d46) | 0.814 / 0.753 / 0.740 (0.769 ± 0.039) | 0.700 / 0.701 / 0.718 (0.706 ± 0.010) | 0.753 / 0.726 / 0.729 (0.736 ± 0.015) |

*Supplementary Results.* On COLO829T Illumina, the v1.1 models begin after epoch 5–6 to call unlabelled tensors in chr1:121–126 and 143–144 Mb (1,075–1,639 false calls at the validation-selected epoch 8, 15–502 at epochs 6 / 5 / 6), with no change in validation F1; Table R3.9a lists both. Outside these segments the v1.1 and v2.1 calls are alike.

*Supplementary Methods (reproducibility).* Training is deterministic for a given seed and dataset (three models retrained with the same seed reproduced their chr1 calls exactly), and re-running prediction and evaluation of 30 in-sample checkpoints on a second cluster reproduced their chr1 F1 after population filtering within 0.002.

*Supplementary Table R3.9b.* Paired 1-Mb block-bootstrap 95% intervals (2,000 resamples) of the differences A − B in chr1 precision, recall and F1 after population filtering (benchmark regions). Pansoma call sets are the first-seed models (COLO829T Illumina on v1.1 at epoch 6, as in R3.1); comparators are DeepSomatic-TO 1.9.0 and ClairS-TO (0.4.2 on HG008T, 0.4.1 on COLO829T), tumor-only, scored with the same population filter.

| comparison | dataset | kind | A | B | Δprecision [95% CI] | Δrecall [95% CI] | ΔF1 [95% CI] |
|---|---|---|---|---|---|---|---|
| HPRC - GRCh38-only | COLO829T Illumina | SNV | HPRC v1.1 (d9) | GRCh38-only graph | +0.442 [+0.056, +0.602] | +0.140 [+0.118, +0.161] | +0.306 [+0.107, +0.441] |
| HPRC - GRCh38-only | COLO829T Illumina | SNV | HPRC v2.1 (d46) | GRCh38-only graph | +0.476 [+0.086, +0.634] | +0.132 [+0.109, +0.153] | +0.316 [+0.116, +0.453] |
| HPRC - GRCh38-only | COLO829T ONT | SNV | HPRC v1.1 (d9) | GRCh38-only graph | +0.454 [+0.283, +0.566] | +0.158 [+0.133, +0.181] | +0.326 [+0.219, +0.422] |
| HPRC - GRCh38-only | COLO829T ONT | SNV | HPRC v2.1 (d46) | GRCh38-only graph | +0.447 [+0.276, +0.569] | +0.149 [+0.122, +0.175] | +0.319 [+0.208, +0.417] |
| HPRC - GRCh38-only | COLO829T fiberseq | SNV | HPRC v1.1 (d9) | GRCh38-only graph | +0.039 [-0.109, +0.179] | +0.098 [+0.077, +0.118] | +0.073 [+0.002, +0.132] |
| HPRC - GRCh38-only | COLO829T fiberseq | SNV | HPRC v2.1 (d46) | GRCh38-only graph | +0.120 [+0.028, +0.211] | +0.083 [+0.064, +0.100] | +0.099 [+0.062, +0.141] |
| HPRC - GRCh38-only | HG008T Illumina | SNV | HPRC v1.1 (d9) | GRCh38-only graph | +0.117 [+0.061, +0.172] | +0.044 [+0.016, +0.075] | +0.116 [+0.062, +0.171] |
| HPRC - GRCh38-only | HG008T Illumina | SNV | HPRC v2.1 (d46) | GRCh38-only graph | +0.094 [+0.040, +0.144] | +0.066 [+0.039, +0.095] | +0.101 [+0.050, +0.152] |
| HPRC - GRCh38-only | HG008T ONT | SNV | HPRC v1.1 (d9) | GRCh38-only graph | +0.072 [+0.044, +0.107] | -0.016 [-0.039, +0.007] | +0.059 [+0.034, +0.091] |
| HPRC - GRCh38-only | HG008T ONT | SNV | HPRC v2.1 (d46) | GRCh38-only graph | +0.038 [+0.009, +0.070] | +0.017 [-0.004, +0.040] | +0.038 [+0.010, +0.069] |
| HPRC - GRCh38-only | HG008T PacBio | SNV | HPRC v1.1 (d9) | GRCh38-only graph | +0.125 [+0.057, +0.205] | +0.001 [-0.024, +0.028] | +0.108 [+0.046, +0.187] |
| HPRC - GRCh38-only | HG008T PacBio | SNV | HPRC v2.1 (d46) | GRCh38-only graph | +0.100 [+0.032, +0.182] | +0.009 [-0.019, +0.037] | +0.090 [+0.028, +0.172] |
| Pansoma - comparator | COLO829T Illumina | SNV | HPRC v2.1 (d46) | ClairS-TO 0.4.1 | +0.436 [+0.239, +0.566] | -0.050 [-0.065, -0.035] | +0.232 [+0.095, +0.353] |
| Pansoma - comparator | COLO829T Illumina | SNV | HPRC v2.1 (d46) | DeepSomatic-TO 1.9.0 | +0.440 [+0.249, +0.566] | -0.120 [-0.135, -0.105] | +0.219 [+0.073, +0.344] |
| Pansoma - comparator | COLO829T ONT | SNV | HPRC v2.1 (d46) | ClairS-TO 0.4.1 | +0.038 [-0.020, +0.082] | +0.007 [-0.015, +0.029] | +0.021 [-0.010, +0.045] |
| Pansoma - comparator | COLO829T ONT | SNV | HPRC v2.1 (d46) | DeepSomatic-TO 1.9.0 | +0.205 [+0.130, +0.268] | -0.003 [-0.025, +0.020] | +0.103 [+0.062, +0.141] |
| Pansoma - comparator | COLO829T fiberseq | SNV | HPRC v2.1 (d46) | ClairS-TO 0.4.1 | +0.086 [+0.053, +0.112] | -0.037 [-0.053, -0.020] | +0.017 [+0.000, +0.032] |
| Pansoma - comparator | COLO829T fiberseq | SNV | HPRC v2.1 (d46) | DeepSomatic-TO 1.9.0 | +0.133 [+0.099, +0.168] | +0.152 [+0.124, +0.182] | +0.147 [+0.123, +0.173] |
| Pansoma - comparator | HG008T Illumina | SNV | HPRC v2.1 (d46) | ClairS-TO 0.4.2 | +0.018 [-0.009, +0.048] | +0.007 [-0.016, +0.032] | +0.018 [-0.008, +0.046] |
| Pansoma - comparator | HG008T Illumina | SNV | HPRC v2.1 (d46) | DeepSomatic-TO 1.9.0 | +0.053 [+0.011, +0.098] | -0.037 [-0.054, -0.023] | +0.042 [+0.002, +0.088] |
| Pansoma - comparator | HG008T ONT | SNV | HPRC v2.1 (d46) | ClairS-TO 0.4.2 | +0.018 [-0.029, +0.075] | +0.029 [+0.004, +0.056] | +0.022 [-0.020, +0.076] |
| Pansoma - comparator | HG008T ONT | SNV | HPRC v2.1 (d46) | DeepSomatic-TO 1.9.0 | +0.122 [+0.039, +0.206] | -0.010 [-0.030, +0.010] | +0.118 [+0.033, +0.213] |
| Pansoma - comparator | HG008T PacBio | SNV | HPRC v2.1 (d46) | ClairS-TO 0.4.2 | +0.017 [+0.003, +0.031] | -0.009 [-0.032, +0.018] | +0.012 [-0.002, +0.027] |
| Pansoma - comparator | HG008T PacBio | SNV | HPRC v2.1 (d46) | DeepSomatic-TO 1.9.0 | +0.047 [+0.028, +0.066] | +0.029 [-0.001, +0.060] | +0.046 [+0.027, +0.067] |
| Pansoma - comparator | COLO829T Illumina | INDEL | HPRC v2.1 (d46) | ClairS-TO 0.4.1 | -0.029 [-0.090, +0.012] | +0.099 [+0.044, +0.160] | -0.005 [-0.051, +0.036] |
| Pansoma - comparator | COLO829T Illumina | INDEL | HPRC v2.1 (d46) | DeepSomatic-TO 1.9.0 | -0.121 [-0.223, -0.043] | +0.061 [-0.007, +0.134] | -0.073 [-0.142, -0.011] |
| Pansoma - comparator | COLO829T ONT | INDEL | HPRC v2.1 (d46) | ClairS-TO 0.4.1 | -0.118 [-0.185, -0.056] | +0.137 [+0.078, +0.205] | -0.101 [-0.158, -0.049] |
| Pansoma - comparator | COLO829T ONT | INDEL | HPRC v2.1 (d46) | DeepSomatic-TO 1.9.0 | -0.031 [-0.082, +0.001] | +0.214 [+0.144, +0.289] | -0.018 [-0.053, -0.007] |
| Pansoma - comparator | COLO829T fiberseq | INDEL | HPRC v2.1 (d46) | ClairS-TO 0.4.1 | -0.293 [-0.424, -0.161] | +0.153 [+0.085, +0.228] | -0.137 [-0.206, -0.067] |
| Pansoma - comparator | COLO829T fiberseq | INDEL | HPRC v2.1 (d46) | DeepSomatic-TO 1.9.0 | -0.303 [-0.467, -0.159] | +0.176 [+0.095, +0.257] | -0.112 [-0.180, -0.048] |
| Pansoma - comparator | HG008T Illumina | INDEL | HPRC v2.1 (d46) | ClairS-TO 0.4.2 | -0.111 [-0.194, -0.028] | +0.160 [+0.130, +0.191] | +0.121 [+0.085, +0.158] |
| Pansoma - comparator | HG008T Illumina | INDEL | HPRC v2.1 (d46) | DeepSomatic-TO 1.9.0 | -0.125 [-0.179, -0.077] | +0.099 [+0.074, +0.124] | +0.045 [+0.017, +0.073] |
| Pansoma - comparator | HG008T ONT | INDEL | HPRC v2.1 (d46) | ClairS-TO 0.4.2 | -0.404 [-0.509, -0.300] | +0.162 [+0.134, +0.192] | -0.112 [-0.149, -0.076] |
| Pansoma - comparator | HG008T ONT | INDEL | HPRC v2.1 (d46) | DeepSomatic-TO 1.9.0 | -0.164 [-0.266, -0.079] | +0.218 [+0.183, +0.254] | -0.019 [-0.041, -0.001] |
| Pansoma - comparator | HG008T PacBio | INDEL | HPRC v2.1 (d46) | ClairS-TO 0.4.2 | -0.388 [-0.493, -0.281] | +0.181 [+0.151, +0.213] | -0.030 [-0.065, +0.006] |
| Pansoma - comparator | HG008T PacBio | INDEL | HPRC v2.1 (d46) | DeepSomatic-TO 1.9.0 | -0.137 [-0.225, -0.064] | +0.225 [+0.192, +0.257] | +0.045 [+0.020, +0.069] |

*Supplementary Table R3.9c.* The seed-matched three-reference comparison behind "35 of 36" (R3.1): chr1 SNV precision, recall and F1 after population filtering of the in-sample models on the GRCh38-only graph and both HPRC graphs, three seeds each, with base and mapping quality standardised by the training statistics (the setting with GRCh38-only seeds; COLO829T Illumina on v1.1 at its validation-selected epoch); last column: seeds on which the HPRC model has the higher F1.

| tumor | platform | reference | precision: seeds 0 / 1 / 2 (mean ± SD) | recall | F1 | HPRC > GRCh38-only, same seed |
|---|---|---|---|---|---|---|
| HG008T | Illumina | GRCh38-only | 0.334 / 0.394 / 0.388 (0.372 ± 0.033) | 0.786 / 0.783 / 0.789 (0.786 ± 0.003) | 0.469 / 0.524 / 0.520 (0.504 ± 0.030) | - |
| HG008T | Illumina | HPRC v1.1 (d9) | 0.452 / 0.458 / 0.460 (0.457 ± 0.005) | 0.831 / 0.819 / 0.801 (0.817 ± 0.015) | 0.585 / 0.587 / 0.585 (0.586 ± 0.001) | 3 / 3 |
| HG008T | Illumina | HPRC v2.1 (d46) | 0.429 / 0.457 / 0.432 (0.439 ± 0.016) | 0.852 / 0.828 / 0.815 (0.832 ± 0.019) | 0.570 / 0.589 / 0.565 (0.575 ± 0.013) | 3 / 3 |
| HG008T | PacBio | GRCh38-only | 0.372 / 0.472 / 0.463 (0.436 ± 0.055) | 0.826 / 0.792 / 0.805 (0.808 ± 0.017) | 0.513 / 0.592 / 0.588 (0.564 ± 0.044) | - |
| HG008T | PacBio | HPRC v1.1 (d9) | 0.479 / 0.443 / 0.464 (0.462 ± 0.018) | 0.829 / 0.808 / 0.831 (0.823 ± 0.013) | 0.608 / 0.572 / 0.595 (0.592 ± 0.018) | 2 / 3 |
| HG008T | PacBio | HPRC v2.1 (d46) | 0.473 / 0.477 / 0.487 (0.479 ± 0.007) | 0.839 / 0.822 / 0.821 (0.827 ± 0.010) | 0.605 / 0.604 / 0.611 (0.607 ± 0.004) | 3 / 3 |
| HG008T | ONT | GRCh38-only | 0.394 / 0.368 / 0.394 (0.385 ± 0.015) | 0.842 / 0.826 / 0.780 (0.816 ± 0.032) | 0.537 / 0.509 / 0.524 (0.523 ± 0.014) | - |
| HG008T | ONT | HPRC v1.1 (d9) | 0.473 / 0.451 / 0.447 (0.457 ± 0.014) | 0.821 / 0.816 / 0.822 (0.820 ± 0.003) | 0.600 / 0.581 / 0.579 (0.587 ± 0.012) | 3 / 3 |
| HG008T | ONT | HPRC v2.1 (d46) | 0.432 / 0.467 / 0.474 (0.457 ± 0.022) | 0.859 / 0.824 / 0.809 (0.831 ± 0.026) | 0.575 / 0.596 / 0.598 (0.589 ± 0.013) | 3 / 3 |
| COLO829T | Illumina | GRCh38-only | 0.391 / 0.289 / 0.499 (0.393 ± 0.105) | 0.583 / 0.583 / 0.564 (0.576 ± 0.011) | 0.468 / 0.386 / 0.530 (0.461 ± 0.072) | - |
| COLO829T | Illumina | HPRC v1.1 (d9) | 0.484 / 0.562 / 0.475 (0.507 ± 0.048) | 0.698 / 0.705 / 0.702 (0.702 ± 0.004) | 0.572 / 0.626 / 0.566 (0.588 ± 0.033) | 3 / 3 |
| COLO829T | Illumina | HPRC v2.1 (d46) | 0.867 / 0.704 / 0.820 (0.797 ± 0.084) | 0.715 / 0.705 / 0.717 (0.712 ± 0.007) | 0.783 / 0.704 / 0.765 (0.751 ± 0.042) | 3 / 3 |
| COLO829T | fiberseq | GRCh38-only | 0.745 / 0.841 / 0.868 (0.818 ± 0.065) | 0.625 / 0.600 / 0.576 (0.600 ± 0.024) | 0.680 / 0.700 / 0.692 (0.691 ± 0.010) | - |
| COLO829T | fiberseq | HPRC v1.1 (d9) | 0.785 / 0.876 / 0.845 (0.835 ± 0.046) | 0.723 / 0.707 / 0.709 (0.713 ± 0.009) | 0.752 / 0.782 / 0.771 (0.768 ± 0.015) | 3 / 3 |
| COLO829T | fiberseq | HPRC v2.1 (d46) | 0.866 / 0.846 / 0.849 (0.854 ± 0.011) | 0.707 / 0.713 / 0.737 (0.719 ± 0.016) | 0.779 / 0.774 / 0.789 (0.780 ± 0.008) | 3 / 3 |
| COLO829T | ONT | GRCh38-only | 0.340 / 0.353 / 0.475 (0.390 ± 0.075) | 0.543 / 0.540 / 0.569 (0.551 ± 0.016) | 0.418 / 0.427 / 0.518 (0.454 ± 0.055) | - |
| COLO829T | ONT | HPRC v1.1 (d9) | 0.793 / 0.816 / 0.833 (0.814 ± 0.020) | 0.701 / 0.693 / 0.673 (0.689 ± 0.014) | 0.744 / 0.749 / 0.745 (0.746 ± 0.003) | 3 / 3 |
| COLO829T | ONT | HPRC v2.1 (d46) | 0.788 / 0.734 / 0.744 (0.755 ± 0.029) | 0.692 / 0.703 / 0.718 (0.704 ± 0.013) | 0.737 / 0.718 / 0.730 (0.729 ± 0.009) | 3 / 3 |

Data and code: `supplementary/R3.9/seed_table.tsv`, `make_seed_table.py` (from the RIS summary of 0a07f91 / b662f00), `three_ref_table.tsv`, `make_three_ref_table.py`, `paired_table.tsv`, `make_paired_table.py` (from `experiments/E16_bootstrap/paired.tsv`).

### R3.10 Low-VAF and INDEL claims

**Comment（要点）。**

- HapMap 的低 VAF 结果要加限定：那是 germline 细胞系，不是 somatic lineage。
- SNV 和 INDEL 分开报告。Pansoma 的优势在 SNV；INDEL 准确度一般，在部分 ONT 设置下还低于 ClairS-TO。整体结论要收窄。
- 用真正的低 VAF somatic benchmark，或者对已验证 somatic 变异做 in-silico dilution（COLO829T 或 HG008T），来补充 HapMap。
- INDEL 按 insertion / deletion、长度、序列 context（homopolymer、tandem repeat）和平台分层，好把低 VAF 生物学、平台错误模式和 INDEL 类别区分开。

**Response.**

We agree with all four points and have revised the manuscript accordingly.

1. We moved the HapMap mixture from the main benchmark to the Supplement. Its positives are germline alleles of other individuals, so it is reported only as detection by allele fraction in that mixture, not as somatic accuracy (R2.1, R2.4).
2. The revised manuscript reports SNVs and INDELs separately throughout and limits any claimed advantage to SNVs. For GRCh38-positioned calls in the benchmark regions after the common population filter, Pansoma's SNV F1 exceeds that of ClairS-TO and DeepSomatic-TO in all six tested settings; for INDELs the results are mixed and we claim no advantage: with the final INDEL models (trained on six cancer lines with chr1 held out, three seeds) at the operating threshold, Pansoma's INDEL F1 (0.09–0.25) is higher than both comparators in three of the six settings (HG008T PacBio and ONT, COLO829T Fiber-seq) and lower in the other three (DeepSomatic-TO on Illumina for both tumors, ClairS-TO on COLO829T ONT). The Abstract and Discussion state this.
3. For somatic variants at reduced purity, we added in-silico dilutions of COLO829T with its matched normal at 75, 50 and 25% purity on three platforms, with recall reported by truth VAF (R2.4).
4. A new Supplementary Table stratifies Pansoma's INDEL performance by insertion or deletion, length, repeat context and platform for both tumors. INDEL calling is limited mainly by population filtering of truth INDELs, errors concentrated in long homopolymers (most false calls on the HG008T long-read platforms) and the choice of threshold; the dilutions contain too few true INDELs to assess the effect of allele fraction (R1.6).

### R3.m1 Parameter discrepancies between Methods and Supplement

**Comment（要点）。**

- 主 Methods（AF > 10 %、2 条 supporting reads、排除 MAPQ ≤ 5、80 epochs）和 Supplementary Methods（AF ≥ 0.08、≥ 3 个 ALT observation、mean ALT BQ ≥ 10、MAPQ ≥ 10、INDEL ≤ 50 bp、80–120 epochs）不一致。
- 上游（candidate）和下游（VCF、过滤、评估）步骤写得比网络部分简略得多。

（审稿人原文不放进 repo：GitHub 仓库是 public 的；原文在 Pansoma_Resubmit/Reviewer_Letter_Response.docx。）

**Response.**

We agree that the inconsistent parameter descriptions hindered reproducibility. We have reconciled the main and Supplementary Methods against the configuration files used for the revised analyses; both now give the values below, which supersede those of the original submission, and describe the upstream and downstream steps in the same detail as the network. The table summarizes the standard settings; Supplementary Table [X] records the settings and exceptions of each reported model (for example input normalization, training data and the purity dilutions).

| Step | Standard setting |
|---|---|
| Node selection | Alignments with MAPQ > 5; a node is selected when more than 5% of them carry a mismatch or indel on it |
| Candidates | Reads with MAPQ > 10; ≥ 3 ALT reads, each with base quality ≥ 10 at the allele (per read, not averaged over reads); allele fraction ≥ 0.08 (SNVs and INDELs); INDELs of 1–50 bp; at most 800 reads per node |
| Tensors | 8 channels × up to 200 reads × 101 columns along each read's graph walk (R1.1) |
| Labels and split | Somatic, germline or other (R1.8); chr2–22 for training, with 5% of 20,000-node blocks held out for validation; chr1 of every evaluated tumor held out for testing (R1.7) |
| Training | ConvNeXt: SNV 50.4 M parameters, 20 training epochs; INDEL 8.7 M, 12 training epochs. AdamW, learning rate 10⁻⁴, weight decay 0.05, one warm-up epoch, then cosine decay; batch 1,024 per GPU; cross-entropy with inverse-frequency class weights to the power 0.75 (SNV) or 0.5 (INDEL) |
| Checkpoint | The epoch with the best validation somatic F1 |
| Calling | PASS at the recall-0.9 validation threshold; INDELs also at the validation best-F1 threshold (R1.6); the optional off-reference rescue for SNVs is off by default; SNV results are reported with and without it (R1.3, R3.4) |
| Population filtering (R2.6) | Exact allele match. SNVs: gnomAD or CoLoRSdb AF ≥ 10⁻⁴, non-somatic dbSNP 138 records, GATK 1000 Genomes panel of normals. INDELs: gnomAD or CoLoRSdb only, AF ≥ 0.01 in homopolymers ≥ 8 bp, otherwise ≥ 10⁻⁴ |
| Evaluation | rtg vcfeval on chr1 within each benchmark's regions, after GRCh38 projection; graph-coordinate calls reported separately (R3.4) |

### R3.m2 Comparator setup

**Comment（要点）。**

- 像描述 Pansoma 一样详细地记录比较器：软件版本、aligner 和版本、预训练 checkpoint、确切的 tumor-only 配置、PASS / filter 标准、参考资源、population / PoN 资源、命令行。
- Pansoma 的主要优势是 precision，而 precision 对过滤配置很敏感。

（审稿人原文不放进 repo：GitHub 仓库是 public 的；原文在 Pansoma_Resubmit/Reviewer_Letter_Response.docx。）

**Response.**

We thank the reviewer and agree. The revised Methods describe both comparators in the same detail as Pansoma, and a new Supplementary Note gives the exact command line of every run.

| Platform | Input GRCh38 alignment (HG008T / COLO829T / HCC1395) | DeepSomatic-TO model | ClairS-TO model |
|---|---|---|---|
| Illumina | BWA-MEM 0.7.17 / Sentieon BWA-MEM 202308.01 / BWA-MEM 0.7.17 | `WGS_TUMOR_ONLY` | `ilmn_ssrs` (COLO829T: `ilmn`, the synthetic-data `ss` model) |
| PacBio HiFi / Fiber-seq | pbmm2 1.10.0 / 1.13.0 / 1.10.0 | `PACBIO_TUMOR_ONLY` | `hifi_revio_ssrs` |
| ONT | minimap2 2.26 | `ONT_TUMOR_ONLY` | `ont_r10_dorado_sup_5khz_ssrs` |

- **Software and inputs.** DeepSomatic 1.9.0 and ClairS-TO 0.4.2 (HG008T, HCC1395) or 0.4.1 (COLO829T) were run from the developers' container images on the providers' GRCh38 alignments (GIAB, SMaHT, CASTLE) of the reads Pansoma used; we aligned the SEQC2 Illumina reads ourselves. The reference was GRCh38 (the no-alt analysis set for HCC1395), and calls were evaluated on chr1. ClairS-TO has no dorado SUP v5 model, so HG008T ONT used its SUP v4.2 model. The HG008T ONT comparison was rerun with DeepSomatic's `ONT_TUMOR_ONLY` model, replacing the original run that used `PACBIO_TUMOR_ONLY`.
- **Population resources and PASS.** DeepSomatic-TO receives population allele frequencies as an input channel (DeepVariant 1000 Genomes calls for Illumina, CoLoRSdb for long reads; AF > 0.05), and we enabled its default panel-of-normals filter (`--use_default_pon_filtering=true`). ClairS-TO used its default tagging (gnomAD and CoLoRSdb AF ≥ 0.001, dbSNP 138 non-somatic sites, the GATK 1000 Genomes panel of normals) and its default quality filters. All other options were left at their defaults. For both tools, PASS means FILTER = PASS.
- **Scoring.** All call sets are scored with Pansoma's protocol (same truth, benchmark regions and rtg vcfeval options). We report native PASS calls and PASS calls after applying the same additional population filter to all callers; precision–recall curves are generated from all records after this common filtering step, ranked by QUAL (R1.6).
- **Training data.** Pansoma's within-tumor models use chr2–22 of the evaluated tumor for training. These results are distinguished from the cross-tumor evaluations in R1.7, where the evaluated tumor is excluded from training. On HCC1395, DeepSomatic's tumor-only models, like Pansoma's multi-tumor model, were trained on this line with chr1 held out; ClairS-TO's `ssrs` models exclude HCC1395.

### R3.m3 Composition of the pangenome reference

**Comment（要点）。**

- 写明 HPRC 的 release 和版本、包含的 haplotypes、allele-frequency 过滤标准、structural variants 的处理、图的构建方法。
- 说明 benchmark 个体或其近亲是否进了评估用的图。

（审稿人原文不放进 repo：GitHub 仓库是 public 的；原文在 Pansoma_Resubmit/Reviewer_Letter_Response.docx。）

**Response.**

We thank the reviewer. The original Methods named the reference only as an "allele-frequency-filtered HPRC pangenome graph"; the revised Methods specify each graph, and Supplementary Table [X] lists the individuals and their ancestry.

- **HPRC v1.1 graph (d9 filter):** HPRC release 1, built with Minigraph-Cactus on GRCh38 (Liao et al., 2023); 90 haplotypes: two from each of 44 individuals, plus GRCh38 and CHM13.
- **HPRC v2.1 graph (d46 filter):** HPRC release 2, built the same way (Lucas et al., 2026); 464 haplotypes: two from each of 231 individuals, plus GRCh38 and CHM13.
- **Frequency filter:** both graphs are used as the HPRC distributes them. The filter removes nodes and edges supported by fewer than 9 of the 90 (d9) or 46 of the 464 (d46) haplotypes, about 10%; the GRCh38 and CHM13 reference paths are never filtered. Allele-carrier counts in this letter (R1.m1) are over the 463 haplotypes other than GRCh38; a non-reference allele kept by the d46 filter is carried by at least 46 of them.
- **GRCh38-only graph** (R1.2, R3.1): built with vg autoindex from the GRCh38 no-alt analysis set.

**Structural variants** pass the same filter, and those kept are graph paths that reads can follow. Pansoma calls only SNVs and INDELs of 1–50 bp; longer events are not called.

**Benchmark individuals.** None of the six cancer lines (HG008, COLO829, HCC1395, HCC1954, H1437, H2009) appears in either graph's sample list. We additionally examined rare-allele sharing for the HG008, COLO829 and HCC1395 normals: of their germline SNVs carried by only one or two haplotypes of a release, the HPRC individual sharing the most carries 8.2% (v1.1) and 2.9% (v2.1). These measurements are supplementary evidence, not a formal exclusion of relatedness; the other three lines were checked by sample list only. The HapMap mixture donors, in contrast, are HPRC individuals; that experiment uses a graph without them (R1.5).

"AF-HPRC" and "AF-Graph" (original Fig. 2e and Fig. 5) are now "HPRC v1.1 graph (d9 filter)" throughout.

### R3.m4 Terminology and language

**Comment（要点）。**

- 术语不一致，有若干语法问题。
- "panel of normals" 对 gnomAD、dbSNP、1000G、CoLoRSdb 来说范围太宽。
- 指出的几个具体句子要改。
- NPU 的全称要统一。
- Pansoma 和 Pansoma-TO 混用；"Supplement Figure" 和 "Supplementary Figure" 混用。
- 最后做一次语言编辑。

（审稿人原文不放进 repo：GitHub 仓库是 public 的；原文在 Pansoma_Resubmit/Reviewer_Letter_Response.docx。）

**Response.**

We thank the reviewer and have standardized terminology throughout the manuscript and Supplementary Information. The four background resources are now described collectively as "population filtering resources," with "panel of normals" reserved for the GATK 1000 Genomes resource and for the comparators' own panel-of-normals filters. Each resource and its filtering rule are specified in Methods.

The revised pipeline constructs tensors directly from sorted, indexed GAM alignments and no longer produces NPU files. We have therefore removed the NPU terminology and the obsolete format description (original Supplementary Table 1).

We consistently use "Pansoma," "Supplementary Fig." and "Supplementary Table," and have standardized reference names, platform labels and candidate terminology. The grammatical errors identified by the reviewer have been corrected, including "Pansoma converts," "technical details" and "read pileup." The manuscript and Supplementary Information have also been edited for language and consistency.
