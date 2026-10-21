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

We agree with the first point. On the second, there was no leakage: the submitted HapMap results were produced on a graph built without any of the six donors, with the same Minigraph-Cactus construction and d9 filter but from 40 HPRC year-1 individuals plus GRCh38 and CHM13 (82 haplotypes; HG00438, HG02257, HG02486, HG02622, HG002, HG005 and NA19240 excluded). The Methods now state this. The HAPMAP6 positives are germline alleles of the five minor donors (VAF 0.25–16.5%), not acquired mutations. Four of the donors are in the HPRC v1.1 graph and all six in the HPRC v2.1 graph (d46 filter), as our companion study noted, which is why the donor-free graph was used. Even with the donors in a graph, the overlap could not inflate the results: apart from the GRCh38 and CHM13 reference paths, our graphs keep only nodes on about 10% or more of their haplotypes, so a donor's private alleles are filtered out, and a donor allele the graph carries is absorbed and removes a positive. The donors' rare alleles are also found in CoLoRSdb, one of our population filtering resources, far more often than those of our cancer benchmarks' normals (chr20 SNVs: 14.5–25.9% against 4.0–7.9%). The benchmark therefore asks the caller to report alleles that its population prior, used by the trained model (R2.1), is designed to suppress.

We removed the HapMap mixture from the main benchmark and withdrew the HapMap-trained low-VAF model; no claim of the revised manuscript rests on HapMap results. Low-VAF performance is now measured on in-silico dilutions of COLO829T with COLO829BL on three platforms (nominal purity 75, 50 and 25%; measured 0.66 / 0.40 / 0.20 for Illumina and Fiber-seq, 0.70 / 0.44 / 0.22 for ONT). One model per platform and variant type is trained on chr2–22 of all four purities and tested on chr1 of each; the tumor itself is not held out.

- Fiber-seq SNV F1 after population filtering, at one recall-0.9 validation threshold pooled over all purities, is 0.779 / 0.687 / 0.665 / 0.720 at 100 / 75 / 50 / 25% purity (95% block-bootstrap interval at 50%: 0.648–0.682; 547–705 graph-coordinate calls per purity are not scored; R3.4). The share of truth SNVs that any model can reach falls from 0.821 to 0.711; on chr1, 16.1% of truth SNVs at full purity and 26.9% at 25% purity lie below the 0.08 candidate allele-fraction floor.
- The dip at 75 and 50% comes from the tumor's loss of heterozygosity: where the tumor lost the ALT copy of a germline heterozygous SNV, the admixed normal reads put that allele at (1 − purity)/2, inside the somatic allele-fraction range of those mixtures. The false positives at these two purities have a germline-like spectrum (cosine to the germline spectrum 0.788 and 0.768, against 0.536 at full purity).
- In the Illumina and ONT mixtures, tumor- and normal-derived reads differ in base quality, a cue that a real tumor lacks, so those scores are likely optimistic. For Illumina INDELs, removing the cue (per-read base-quality centring) lowers the dilutions' test-tuned best F1 (upper bound) to 0.133–0.146, below the 100% sample (0.182), as expected if the cue had inflated them. For Illumina SNVs, the dilutions' PASS F1 falls from 0.804–0.840 to 0.678–0.709 (precision from 0.97 to 0.68–0.74) [待补: ONT after removing this cue].

[待补: dilutions of a second tumor not used in training (HCC1395)]

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
   for HG008T SNVs (three platforms each; the three-seed mean is higher on every platform), and for INDELs from the
   13–19% (HG008T) and 4–22% (COLO829T) of the original submission to 24–26% and 23–26%.
2. *Recall at equal precision.* Under the identical protocol (R3.m2), at each comparator's post-filter precision,
   Pansoma's SNV recall on HG008T is 0.856–0.898 (v2.1, three-seed means), against 0.806–0.889 for DeepSomatic-TO and
   ClairS-TO. It is lower only than DeepSomatic-TO's on Illumina (HG008T 0.871 against 0.889; COLO829T 0.781 against
   0.835).
3. *Where recall is lost.* A per-allele funnel places each missed truth variant at its first loss. For HG008T SNVs the
   largest loss is the threshold; on COLO829T, 15–16% of truth SNVs are subclonal mutations below the 0.08 candidate
   allele-fraction floor. Most truth INDELs are population alleles: under the same population-filtering policy, the
   maximum attainable recall is 0.342 (HG008T) and 0.382 (COLO829T), and after that filtering the comparators reach
   0.015–0.183.

**Supplementary (R1.6).**

*Recall before and after the revision* (chr1, after population filtering; original submission: the chr1 results of
the submitted models on the HPRC v1.1 graph (d9 filter); revision: PASS at the recall-0.9 validation threshold, ranges
over three platforms and, for SNVs, three training seeds):

| | HG008T SNV | COLO829T SNV | HG008T INDEL | COLO829T INDEL |
|---|---|---|---|---|
| original submission (v1.1) | 0.759–0.828 | 0.575–0.596 | 0.13–0.19 | 0.04–0.22 |
| revision, HPRC v1.1 | 0.803–0.831 | 0.675–0.723 | 0.228–0.256 | 0.244–0.275 |
| revision, HPRC v2.1 (d46 filter) | 0.814–0.859 | 0.700–0.734 | 0.237–0.256 | 0.229–0.260 |

*Where recall is lost* (chr1, benchmark regions, HPRC v2.1 graph (d46 filter), models trained on chr2–22 of the same
tumor; ranges over three platforms and, for SNVs, three training seeds; INDEL candidates include partial allele
matches):

| Recall | SNV HG008T | SNV COLO829T | INDEL HG008T | INDEL COLO829T |
|---|---|---|---|---|
| truth allele reaches a candidate | 0.948–0.954 | 0.812–0.817 | 0.611–0.757 | 0.565–0.695 |
| PASS | 0.841–0.900 | 0.716–0.751 | 0.375–0.431 | 0.382–0.443 |
| PASS, after population filtering | 0.814–0.859 | 0.700–0.734 | 0.237–0.256 | 0.229–0.260 |
| any threshold, after population filtering | 0.889–0.900 | 0.788–0.795 | 0.237–0.257 | 0.244–0.267 |
| maximum attainable under the same population-filtering policy (whole chr1 truth) | 0.913 | 0.966 | 0.342 | 0.382 |

- *Operating point.* PASS is the recall-0.9 validation threshold, set on held-out validation blocks before population
  filtering. Records down to p = 0.01 are kept as LowQual and population matches are tagged, not deleted, so users can
  move along the precision–recall curve.
- *Funnel.* On all three references (GRCh38-only, HPRC v1.1, HPRC v2.1 with the d46 filter), HG008T SNVs are called at
  81.8%, 82.0% and 83.2% (HPRC graphs: mean over three seeds), the threshold being the largest loss on every reference
  and seed; candidate-stage losses are mostly graph alleles and repeat absorption that filtering would tag anyway
  (R1.4). On COLO829T, 17–18% of truth SNVs get no candidate, mostly subclonal mutations below the 0.08 allele-fraction
  floor (15.1–16.1% of chr1 truth SNVs), whereas recall at allele fraction 0.05–0.2 is 0.850–0.966 once a candidate
  exists (R2.3, R2.4).
- *INDELs.* 418 of the 635 HG008T chr1 truth INDELs (81 of 131 on COLO829T) are exact population alleles at the
  filtering frequencies, which sets the maximum attainable recall above; our post-filter recall at the recall-0.9
  validation threshold is close to the curve's maximum. A GRCh38-only graph gives more HG008T truth INDELs a candidate
  (0.75–0.87 against 0.61–0.76 on v2.1), but these are almost all population alleles, so 24–27% of truth INDELs are
  called on every reference and both tumors. At the validation best-F1 threshold, long-read INDEL F1 rises (HG008T
  PacBio, v1.1: 0.092 to 0.227) while recall falls to 0.073–0.158 on HG008T (v1.1) and 0.015–0.061 on COLO829T (v2.1);
  we report both operating points. On HCC1395 (SEQC2), only 7.9% of truth INDELs are population-tagged and 0.76–0.85
  reach a candidate (v2.1; HG008T: 0.43–0.50). The Abstract and Discussion state the INDEL limitation (R3.10).
- *Comparators* (R3.m2). At their own PASS, DeepSomatic-TO and ClairS-TO reach HG008T SNV recall 0.816–0.911 and
  0.845–0.859 (precision 0.167–0.359), and after the same population filtering 0.806–0.889 and 0.831–0.845, against
  Pansoma's 0.814–0.859 (v2.1); on COLO829T they reach 0.555–0.835 and Pansoma 0.700–0.734, DeepSomatic-TO on Illumina
  being the highest (R2.3). For INDELs their post-filter recall is 0.019–0.141 on HG008T and 0.015–0.183 on COLO829T.
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
- *Augmentation.* Joint COLO829T models add in-silico purity dilutions (75, 50, 25%; R1.5). [待补: ONT after removing a
  read-source base-quality cue (R1.5; done for Illumina)]

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

**Response.**

We added HCC1395 (SEQC2 truth v1.2.1, 41,072 somatic alleles) and three CASTLE lines, HCC1954, H1437 and H2009 (25,014, 84,744 and 161,091 somatic alleles on chr1–22), each sequenced with Illumina, PacBio HiFi and ONT and processed with the identical pipeline on the HPRC v2.1 graph (d46 filter).

*HCC1395 chr1 with models that never saw HCC1395* (preliminary: one training run, earlier quality normalization; 3,440 somatic SNVs in the benchmark region; F1 after population filtering at the recall-0.9 validation threshold of the training tumor; R1.7): models trained on HG008T reach 0.790 (Illumina), 0.698 (PacBio) and 0.689 (ONT); models trained on COLO829T (Fiber-seq for PacBio) reach 0.269, 0.741 and 0.621. Their test-tuned best F1 (upper bound) is 0.810, 0.819 and 0.691, so most of their shortfall comes from the transferred threshold, and part from a UV preference learned on COLO829T.

*CASTLE lines, chr1* (preliminary; same models and protocol): models trained on HG008T reach 0.787–0.937 on H1437 and H2009, and models trained on COLO829T 0.454–0.809, again losing mostly at the transferred threshold (test-tuned best F1 0.828–0.950). On HCC1954 every model has low precision (0.22–0.37): 93–99% of these false positives are germline alleles of its matched normal, while its truth set shows no germline content (2.4% common population alleles, no loss of heterozygosity in the normal), so this is germline leakage by our models, not a truth artefact. The H1437 truth, in contrast, includes about 30% germline alleles that its B-lymphoblastoid normal lost through loss of heterozygosity (chr9q, chr1:9–12 Mb), which explains the line's candidate coverage of about 70%; on chr1, our test chromosome, they are negligible (R2.2). A multi-tumor model trained on HCC1395 (without chr1 and chr20) and the three CASTLE lines reaches SNV F1 0.845 / 0.820 / 0.658 on HCC1395 chr1 (Illumina / PacBio / ONT; 0.850 / 0.825 / 0.730 with HG008T and COLO829T also in training), above every single-tumor model of another line, and generalises to HG008T, which it never saw (0.605 / 0.614 / 0.587; R1.7). Under the same protocol on HCC1395 chr1, DeepSomatic-TO, whose models were trained on HCC1395 with chr1 held out, reaches 0.848 / 0.839 / 0.760, and ClairS-TO, which never saw HCC1395, 0.832 / 0.814 / 0.760: Pansoma is level with DeepSomatic-TO on Illumina and 0.014 / 0.030 below it on PacBio / ONT, and above ClairS-TO except on ONT. On HG008T, which none of the three was trained on, the multi-tumor model is above both (DeepSomatic-TO 0.457–0.557, ClairS-TO 0.552–0.591). For INDELs DeepSomatic-TO leads on HCC1395 (best F1 0.519–0.672 against 0.255–0.662).

How we use these lines follows from how the comparators were trained.

- DeepSomatic's multi-cancer models were trained on HCC1395, HCC1937, HCC1954, H1437 and H2009 with chr1 held out, our test chromosome too, and the CASTLE truth sets of HCC1954, H1437 and H2009 were derived from DeepSomatic's calls (biased toward it, as its authors note). ClairS-TO's real-tumor model was fine-tuned on HCC1937, HCC1954, H1437 and H2009 (R3.m2). We therefore use these three lines for training diversity and held-out-line tests of Pansoma only, and do not compare the tools on them.
- HCC1395 has a truth set independent of both tools (SEQC2). ClairS-TO excluded it from training; DeepSomatic held out its chr1 only, whereas our models above never saw HCC1395. HCC1395 chr1 is therefore our main external comparison, together with HG008T, on which neither tool was trained.
- The SEQC2 and CASTLE truth sets likely omit many INDELs in long homopolymers and STRs (15.5% of SEQC2 INDELs lie in homopolymers of ≥ 7 bp, against 74.9% for HG008T), so unmatched INDEL candidates in such repeats are not used as negatives on these lines (R2.2) and their INDEL results are reported separately.

We have not added HCC1937 and Hs578T, the other two lines of the DeepSomatic study: both are breast cancer lines, a type already represented by HCC1395 and HCC1954, and each further line adds three platforms of data to every training run.

### R1.m2 Tensor windows truncated at node boundaries

**Comment（要点）。**

- 窗口限制在当前节点内、在节点边界截断，相邻节点和跨路径的证据丢掉了。
- 短节点、单碱基节点很难处理；很多 alternative 节点比窗口还短。
- 可信的 graph-native 方法需要跨相邻节点的 path-aware tensor，并显式建模 edges、bubbles、alternative paths 和 competing path support。

**Response.**

The reviewer describes the first version correctly; the revised tensors are not truncated at node boundaries.

- Each tensor (8 channels × 200 rows × 101 columns) is centred on the candidate site, and each row follows its read's alignment walk across node boundaries, on reference and alternative nodes alike. A short or single-base node thus keeps the context of its neighbours.
- The graph base (channel 5) and the number of haplotype paths through the node (channel 6, relative to the graph's haplotype count) are filled per column from the node the read traverses there.
- Rows are grouped by allele (A1, A2, …, REF, OTHER) and clustered within a group by their events, which include every node transition, so reads taking the same alternative path sit together.

The classifier uses this path information.

- Alone, the site column's path count separates somatic from germline chr1 SNV test tensors with an AUC of 0.831–0.865 on the HPRC v1.1 graph (d9 filter) and 0.885–0.915 on the HPRC v2.1 graph (d46 filter), and the number of distinct codes within ±10 columns, a proxy for node changes, with 0.803–0.842 (all 12 graph-based test sets; 0.500 on the GRCh38-only graph).
- We edited channel 6 at inference, without retraining, in the same 12 SNV models (HG008T and COLO829T, three platforms, two graphs). With the plane set to its training mean, somatic average precision falls from 0.199–0.352 to 0.044–0.066 (HG008T) and from 0.639–0.829 to 0.329–0.520 (COLO829T). Setting every node to "all haplotypes" has a similar effect and raises the germline tensors called somatic 3.4- to 7.6-fold. On the GRCh38-only graph, where the plane is constant, nothing changes.
- After population filtering (a per-tensor estimate), these two changes to the plane alter F1 by only −0.066 to +0.007, except HG008T ONT on v1.1 with the training mean (0.596 to 0.145); replacing the plane with another tensor's lowers F1 to 0.212–0.551. Whether the plane improves the final calls is the retraining ablation of R1.1.
- The plane also suppresses somatic SNVs off GRCh38: of the assembly-derived HG008T somatic SNVs whose allele lies on an off-reference node, the models call 0 of 112 (PacBio) and 1 of 89 (ONT) on v1.1, against 61 and 49 with every node set to all haplotypes (v2.1, three platforms: 1–4 of 67–135, against 49–97), although 57–79% per dataset have read support as somatic (R1.1). An optional off-reference rescue recovers 22–66 of the read-supported ones per dataset at tensor level (R3.m1); we report every SNV result with and without the rescue (R3.4).

Supplementary Fig. R1.m2 shows one such tensor: at a somatic SNV, its rows cross several nodes, and the ALT reads follow a node carried by 10–12 of the 90 haplotypes about 15 bp upstream, while the REF reads follow one carried by 89.

We agree that the tensors do not model bubbles explicitly: competing path support is represented through the reads' walks and per-node haplotype counts, not as a bubble-level object, and we state this as a limitation (R1.1).

**Supplementary (R1.m2).**

![A chr1 SNV tensor whose rows follow their reads across node boundaries](supplementary/R1.m2/fig_R1.m2_tensor_across_nodes.png)

*R1.m2,* the eight channels of one chr1 SNV tensor (HG008T PacBio HiFi, HPRC v1.1 graph with the d9 filter; somatic truth
SNV chr1:242,680,268 G>A, allele fraction 0.50; 49 ALT, 47 REF and 2 other reads). Rows are reads, grouped by the
allele they carry at the site (A1 = ALT, REF, OTHER); columns are the 101 positions centred on the site (column 50).
Each row follows its read's alignment walk, so the window spans several graph nodes. Panel 6 gives the graph base of the
node the read traverses at each column, panel 7 the number of haplotype paths through that node (coded relative to the
graph's 90 haplotypes: 99 = 89 haplotypes, 10 = 10). In columns 30–37 the ALT reads traverse a node carried by 10–12
haplotypes and the REF reads one carried by 89, so the tensor records that the somatic allele lies on the haplotype with
the less common upstream allele. Rendered with `scripts/visualize_tensor.py` (`supplementary/R1.m2/make_figure.sh`).

## Reviewer 2

### R2.1 Why include the HapMap benchmark

**Comment（要点）。**

- HapMap 的 positives 是 germline variant，它们的 donor 在 HPRC 图里，所以理论上 pangenome 比对会把它们去掉。
- 请讨论为什么纳入这个 benchmark，以及它在模拟 somatic mutation 上的局限。

**Response.**

We included the SMaHT HAPMAP6 mixture because it offers a large truth set at very low allele fractions (0.25–16.5%). We agree that it has three limits for a pangenome-based tumor-only caller:

1. Its positives are germline alleles of the five minor donors that the HG005 background lacks. They are not acquired mutations, so they carry no somatic mutational process, clonal structure or tumor copy-number context.
2. Four of the donors are in the HPRC v1.1 graph (d9 filter), and all six are in the HPRC v2.1 graph (d46 filter). Apart from their GRCh38 and CHM13 reference paths, both graphs keep only nodes carried by about 10% or more of their haplotypes, so they absorb the positives that are common among those haplotypes; a donor's private alleles are filtered out even when the donor is in the graph. The trained classifier relies on this prior: on our tumor benchmarks, when the haplotype-count channel is set at inference so that every node appears to carry all haplotypes, 3.4–7.6 times as many germline tensors pass the calling threshold (12 HPRC-graph SNV models, chr1), although population filtering removes most of them again (R1.1).
3. The donors are public reference samples. Their rare alleles are found in CoLoRSdb, one of our population filtering resources, far more often than those of our cancer benchmarks' normals (chr20 SNVs: 14.5–25.9% against 4.0–7.9%; R2.6), so filtering also works against these positives.

The mixture therefore measures sensitivity to low-fraction germline alleles. It does not measure tumor-only somatic calling.

We removed it from the main benchmark, and no claim of the revised manuscript rests on HapMap results. In its place we use tumor/normal dilutions that contain true somatic mutations: COLO829T at 75, 50 and 25% purity on three platforms, with one model per platform trained on chr2–22 of all purities and tested per purity on chr1; the tumor itself is not held out (R1.5, R2.4). [待补: Illumina SNVs and ONT after removing a read-source base-quality cue (R1.5; done for Illumina INDELs); dilutions of a second tumor not used in training (HCC1395)]

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
serve as cross-sample tests (R1.7).

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

*Recall under a common protocol.* We re-scored the original DeepSomatic-TO and ClairS-TO call sets with our protocol
(chr1, benchmark regions, the same population filter; R3.m2). On HG008T, SNV recall is now comparable: 0.814–0.859 for
Pansoma on the HPRC v2.1 graph (d46 filter; three training seeds; v1.1: R1.6), against 0.806–0.889 for DeepSomatic-TO
and 0.831–0.845 for ClairS-TO, with Pansoma's models trained on chr2–22 of the same tumor. On COLO829T, Pansoma reaches
0.700–0.734 and the comparators 0.555–0.835, and a substantial deficit remains against DeepSomatic-TO on Illumina
(0.705–0.717 against 0.835). This deficit lies both among truth SNVs with an allele fraction below 0.08 (30% recovered
against 6–7%; Pansoma's floor applies to the allele fraction in the reads, which exceeds 0.08 at a few of these sites)
and above it (94% against 83–85%), so the candidate floor does not explain it fully; we address low-allele-fraction
detection in R2.4. INDEL recall remains low in absolute terms, although under the same protocol it is higher than the
comparators' on both tumors (0.229–0.260 against 0.015–0.183).

*Where recall is lost.* A per-allele funnel assigns each missed truth allele to the first stage at which it is lost
(R1.6). On HG008T the largest loss on every reference is the classification threshold (8.9–10.9% of truth SNVs, against
2.1–5.2% without a candidate and 3.2–5.2% first lost at filtering). The population filter tags 61 of the 702 truth SNVs
(8.7%), which any caller using the same filter loses; the tagged alleles not counted at the filtering stage are already
lost before it. On COLO829T, candidate generation is the main limitation: 17–18% of truth SNVs get no candidate
(GRCh38-only graph and HPRC v2.1 graph), mainly because 15.1–16.1% have a truth allele fraction below the 0.08 candidate
floor.

*The two categories.* We regret that the original labels were not explained. "Variant node not ranked #1" and "first
aligned node mismatch" described implementation-specific failure modes of the original pipeline's node-level candidate
selection, which the revised pipeline replaces. We have therefore replaced them in Supplementary Fig. 2 with one
read-level class for every truth allele without a candidate, on six datasets: no read with MAPQ > 10; fewer than three
ALT reads or an allele fraction below the threshold; the ALT is an existing graph allele; repeat absorption or
respelling; an indel longer than 50 bp. Losses after candidate generation, at the threshold and at filtering, are
reported separately (R1.6). For example, on the HPRC v1.1 graph, 2,739 (92.0%) of the 2,976 HG008T PacBio truth INDELs
without a candidate (chr1–22) are fully or partly encoded in the graph. [待补: the same classes on the HPRC v2.1 tensor
sets (E1.7, running)]

*The normal at missed sites on graph nodes.* We examined the matched normal directly, using its haplotype-resolved
assemblies (HG008-N v6.2; COLO829BL). [待补: E20, COLO829BL ONT reads at the absorbed sites] On the HPRC v2.1 graph (d46
filter), 1,607 HG008T and 401 COLO829T truth alleles get no candidate on at least one platform although a graph path
spells their ALT and the GRCh38-only graph gives them a candidate. Only 110 (6.8%) and 27 (6.7%) are carried by the
matched normal; 1,481 (92.2%) and 371 (92.5%) are absent from both normal haplotypes, and their ALT paths run through
nodes walked by other HPRC haplotypes (median 138 and 110). At the 1,351 HG008T INDEL loci absent from the normal, the
HG008-T tumor assembly carries a sequence found in neither normal haplotype at 1,347 and exactly the truth allele at
1,153, and the tumor's Illumina reads show the expected allele at 97.4% of them (95.7% at random non-absorbed truth
INDELs). The HPRC v1.1 graph of the original submission gives the same picture (Supplementary block). At these HG008T
loci the absorbed alleles are therefore supported somatic events that escape candidate generation because the graph
already represents them, and the manuscript now describes graph absorption as a cause of missed somatic alleles rather
than as evidence of benchmark error. Most of these alleles also match population resources (the population filter tags
97.0% of the HG008T and 98.5% of the COLO829T ones) and would be removed by population filtering on any reference (R1.4;
in the chr1 benchmark regions, the GRCh38-only graph calls at most two of the 39–93 HG008T truth INDELs per platform
that v2.1 absorbs), but this overlap does not establish that they are germline in the patient. COLO829T has no tumor
assembly; its matched normal carries the exact ALT at 26 (7.5%) of the 347 absorbed INDELs and an ALT of the same length
with other bases at 22 more, so up to 48 (13.8%) may be germline-like, and we state this uncertainty explicitly.

**Supplementary (R2.3).**

- *The HPRC v1.1 graph (original submission).* By the read-level definition used originally (tumor reads align perfectly
  to an existing graph path), 2,478 HG008T and 497 COLO829T truth alleles are absorbed; 110 (4.4%) and 34 (6.8%) are
  carried by the matched normal and 2,085 (84.1%) and 430 (86.5%) are absent from both normal haplotypes. At the 1,935
  HG008T INDEL loci absent from the normal, the tumor assembly carries a sequence in neither normal haplotype at 1,931
  and exactly the truth allele at 1,652, and the tumor's Illumina reads show the expected allele at 95.4%. With the v2.1
  definition above, v1.1 absorbs 2,002 HG008T and 446 COLO829T truth alleles, of which 104 (5.2%) and 26 (5.8%) are
  carried by the normal; the GRCh38-only graph calls at most one of the 38–114 HG008T truth INDELs per platform that
  v1.1 absorbs (chr1, benchmark regions).
- *Definitions.* Absorbed on v2.1: no candidate on at least one platform, a tensor on the GRCh38-only graph, and an HPRC
  v2.1 path that spells the ALT exactly. Carried by the normal: a normal haplotype has the ALT at the locus.
  Germline-like (COLO829T): no normal haplotype carries the ALT, but one has an allele of exactly the ALT length with
  other bases.

### R2.4 Low-VAF detection and the low-VAF model

**Comment（要点）。**

- 10 % 的 read-support 阈值会限制低 purity、高异质性真实样本里的低 VAF 检测：真正的低 VAF SNV 应该怎么测？
- 专门的 low-VAF model 描述不清：它和主模型有什么不同？
- 它能不能找回低 VAF 变异？主模型在多低的 VAF 开始变差？要求做 VAF 递减的 sensitivity 分析。
- 两个模型各自在什么情况下用？

**Response.**

*Thresholds.* A candidate now needs an allele fraction of at least 0.08 and at least three ALT-supporting reads (base quality ≥ 10), for SNVs and INDELs alike. Below this floor a variant never reaches the classifier; on chr1 of undiluted COLO829T this concerns 15.1–16.1% of the truth SNVs.

*Sensitivity by VAF.* COLO829T, chr1, benchmark regions (for COLO829T the whole genome), HPRC v2.1 graph (d46 filter), three platforms, PASS calls after population filtering. A truth SNV's VAF is its allele fraction in the tumor reads, or the SMaHT truth VAF when it has no candidate:

| truth VAF | truth SNVs | recall | recall of those with a candidate |
|---|---|---|---|
| < 0.05 | 284–301 | 0 | – |
| 0.05–0.10 | 139–176 | 0.230–0.348 | 0.853–0.941 |
| 0.10–0.20 | 149–177 | 0.672–0.897 | 0.850–0.966 |
| 0.20–0.30 | 59–88 | 0.763–0.853 | 0.804–0.889 |
| 0.30–0.50 | 266–391 | 0.808–0.857 | 0.840–0.861 |
| ≥ 0.50 | 1,356–1,507 | 0.836–0.857 | 0.841–0.867 |

Recall falls below 50% in the 0.05–0.10 bin, where only 34–68 truth SNVs have a candidate; once a candidate exists, recall is 0.804–0.966 in every bin. HG008T has too few low-VAF truth SNVs on chr1 for this analysis.

*Decreasing purity.* For each platform, one model was trained on chr2–22 of COLO829T and its in-silico dilutions with COLO829BL (nominal purity 75, 50, 25%) and tested on chr1 of each purity at one recall-0.9 validation threshold shared by all purities; training and test share the tumor (one seed). On Fiber-seq, after population filtering and in the benchmark regions (graph-coordinate calls not scored; R3.4), SNV recall is 0.698 / 0.639 / 0.621 / 0.663 from 100% to 25% purity, at precision 0.881 / 0.741 / 0.715 / 0.786, while the share of truth SNVs whose expected VAF is below the 0.08 floor rises from 16.1% to 26.9%. At 25% purity, recall is close to the highest the curve reaches after filtering (0.663 against 0.688), so candidate generation is the main limit; by VAF it is 0.899 at 0.10–0.20 (829 of 922 truth SNVs) and 0.420 at 0.05–0.10. At 75% and 50% the classifier loses more (curve maxima 0.767 / 0.740) and precision falls: 86.1% and 86.5% of the false calls carry a germline allele of the normal (68.4% undiluted), as expected when germline SNVs in the tumor's regions of loss of heterozygosity enter the somatic allele-fraction range (R1.5). [待补: Illumina SNVs and ONT after removing a read-source base-quality cue that makes their current scores optimistic (R1.5; done for Illumina INDELs); dilutions of a second tumor not used in training (HCC1395)]

*The low-VAF model.* The original low-VAF model was the same kind of network trained on the HapMap mixture, with candidates built at an allele fraction of 0.015 instead of 0.10. Its positives were germline alleles, so it is withdrawn with that benchmark (R1.5).

*Usage.* We recommend the default model for expected somatic VAFs of about 0.10 and above, where call-level recall in the 0.10–0.20 bin is 0.672–0.897 in undiluted COLO829T and 0.899 at 25% purity. No truth SNV below VAF 0.05 was called; such variants need a lower candidate floor or a dedicated ultra-low-VAF method.

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

**Comment（要点）。**（原审稿文字不在 repo 里；下面按详版的逐条 Response 反推。）

1. Figure 2d、2f 画的是 truth-set variants 还是 Pansoma 的输出，要说清楚。
2. 图的引用顺序（2e 在 2f 之前引用）。
3. Supplementary Fig. 2 的 COLO829T panel 没有标注。
4. Supplementary Fig. 1 的点标签压在点上。
5. 第 13 页的 "truth overlapping" 没有定义。
6. Figure 4 里把 truth-set 的 spectrum 和各 caller 放在一起（同 R3.5）。
7. Figure 4 的 panel 引用顺序和正文不一致。
8. 定义 candidate 和 candidate site。

**Response.**

We thank the reviewer and made the following changes.

1. Figures 2d and 2f show truth-set variants, not Pansoma output. The legends now say so, and both panels moved to the
   Supplement.
2. Figures are now cited in order (2e before 2f).
3. The COLO829T panel of Supplementary Fig. 2 is labelled.
4. Point labels in Supplementary Fig. 1 were moved off the points.
5. "Truth overlapping" was replaced by defined categories, for example a call within 10 bp of a truth allele but with a
   different allele.
6. The SBS96 spectrum of each truth set is now shown beside each caller's in Figure 4 (R3.5). On chr1, in the
   benchmark regions, Pansoma's calls after population filtering match their truth spectrum (cosine 0.996–0.997 on
   COLO829T and 0.940–0.958 on HG008T, HPRC graphs). [待补: comparator panels]
7. Figure 4 panels are cited in the order of the text (see item 2).
8. A *candidate* is one allele at a graph position: node, start on the node's forward strand, REF, ALT and type (SNV,
   insertion or deletion). A *candidate site* is all candidates that pass the filters at one node position and variant
   class (SNV or INDEL); they are encoded together in one tensor, each read labelled with the allele it carries (or as
   REF or OTHER).

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

We agree. We no longer present graph-coordinate calls as a strength, and every run also reports a complete output counting them as false positives (R3.4). We validated the graph-coordinate SNV PASS calls on chr1 of nine final runs: HG008T on the HPRC v1.1 graph (d9 filter) and the HPRC v2.1 graph (d46 filter), three platforms each, and COLO829T on v2.1 (115–4,569 calls per run; for COLO829T ONT a random 3,000) [待补: COLO829T on v1.1].

1. *Origin.* Most are produced by an off-reference rescue step; without it, 6–162 per run remain (R3.4). The classifier treats the low haplotype counts of nodes off GRCh38 as evidence against a somatic call (R1.1). For 57–79% per dataset of the assembly-derived HG008T off-reference somatic SNVs (all autosomes), the consensus of their ALT reads is in the HG008-T tumor assembly and in neither HG008-N haplotype. Of these read-supported SNVs the models call 0–4 per dataset without the rescue and 22–66 with it (tensor level), and 63–87% of the rescue's calls on these sites are read-supported.
2. *Placement.* Via the GRCh38 nodes flanking their bubble, 2.2–19.1% per run become a GRCh38 allele, adding seven true positives after population filtering across the twelve runs (0–2 per run; four distinct truth SNVs). Of 89–1,718 calls per run placed only to a GRCh38 interval, one matches a missed truth variant (R3.4).
3. *Assemblies.* We aligned each call's ALT context to the HG008-T tumor assembly and both HG008-N haplotypes (COLO829T: the two COLO829BL haplotypes only, so somatic status cannot be confirmed). In controls (150 per group and run), 145–149 placed true positives are tumor-only on HG008T (COLO829T: 146–150 absent from the normal), and 135–146 (HG008T) and 92–134 (COLO829T) placed false positives after filtering are in the normal. In exclusive classes, the graph-coordinate calls are population alleles (ALT context spelled by GRCh38 or an HPRC haplotype or walk) in 49.8–91.4%, germline in 4.2–45.4%, and tumor-only in 0.2–0.9% on HG008T (COLO829T: 2.2–5.8% absent from the normal). By the assemblies alone, 1.1–5.2% per HG008T run are tumor-only, most of them also population alleles.
4. *Graph-native population filter.* It tags 0–1 of 150 placed true positives per run, but also 74.2–88.3% of the assembly-derived off-reference somatic SNVs (all chromosomes, six HG008T sets), and 20 of the 26 PASS calls of such SNVs on chr1 (six runs). It describes calls and cannot certify them.
5. *Cross-platform.* On v2.1, 81 calls (same node, offset and allele) are PASS on all three HG008T platforms, 68–69 of them population alleles.

So on HG008T at most 1.1–5.2% of graph-coordinate calls per run have assembly support as somatic; the rest are mostly population or germline alleles that a GRCh38-based filter cannot reach, and counting them as false positives (R3.4) is close to their true status. Because the rescue also recovers read-supported off-reference somatic SNVs (point 1), we report every SNV result with and without the rescue; that training with assembly-derived off-reference labels would address the cause is stated as a limitation.

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

   | | with the off-reference rescue | without it |
   |---|---|---|
   | graph-coordinate share of PASS calls | 8.0–68.0% (highest for ONT) | 0.5–10.1% |
   | complete-output F1 (GRCh38-only graph: 0.197–0.680) | 0.342–0.736 | 0.524–0.780 |
   | graph-coordinate calls with a CHM13 allele | 36.7–65.9% per run | 76.6% (663 calls) |

   With the rescue, the HPRC graphs fall below the GRCh38-only graph in four of the twelve comparisons, all long-read
   (HG008T ONT on both graphs, COLO829T ONT on v2.1, COLO829T Fiber-seq on v1.1). Without it, all twelve are above it
   except HG008T ONT on v2.1 inside the benchmark regions (0.533 against 0.537).
2. *Why the rescue exists.* Our truth sets are on GRCh38, so the labelled training examples lie almost entirely on nodes
   of the GRCh38 path, which nearly every haplotype traverses (93.6–99.7% of the somatic examples; all germline and
   non-variant ones). Through the path-count channel, the classifier therefore scores variants on nodes carried by few
   haplotypes low (R1.1). The optional rescue re-scores tensors on nodes off GRCh38 with the path count of a node that
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
   somatic calls: with or without one, about 95% of the SNV calls are HPRC or germline alleles, and 1–2% are tumor-only
   or absent from the normal in the assemblies (R3.3). A Supplementary Table lists each call's node, class, GRCh38
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

*Low VAF.* We agree: the HapMap mixture, whose positives are germline alleles of other individuals, is no longer a main benchmark (R1.5). It is replaced by in-silico dilutions of COLO829T, whose somatic variants are validated, at 75, 50 and 25% purity on three platforms (R2.4). On Fiber-seq, one model gives SNV F1 0.665–0.779 across purities and SNV recall 0.779–0.899 at truth allele fractions of 0.1–0.2, 0.212–0.420 at 0.05–0.1 and zero below 0.05, where no candidate is built. The Illumina and ONT dilutions are reported only after removing a base-quality cue that separates their tumor and normal reads (R1.5) [待补].

*SNVs and INDELs are reported separately* throughout, and any claimed advantage is limited to SNVs where the re-scored comparators support it. Below: chr1, HPRC v1.1 graph (d9 filter) unless stated, benchmark regions, after population filtering, same-tumor models (graph-coordinate calls: R3.4).

- SNV F1 at the recall-0.9 validation threshold is 0.585 / 0.621 / 0.596 on HG008T (Illumina / PacBio / ONT) and 0.778 / 0.781 / 0.706 on COLO829T (Illumina / Fiber-seq / ONT).
- INDEL F1 is far lower: 0.251 / 0.092 / 0.015 on HG008T at that threshold and 0.270 / 0.262 / 0.144 as test-tuned best F1 (upper bound); on COLO829T 0.104 / 0.024 / 0.003 and 0.161 / 0.142 / 0.095. At the validation best-F1 threshold it is 0.218 / 0.227 / 0.117 on HG008T; on COLO829T that threshold passes only 14–32 calls (F1 0.041–0.080).
- ClairS-TO and DeepSomatic-TO are re-scored under the same protocol (R3.m2) [待补: numbers]; we state every setting where ClairS-TO is better, including ONT INDELs.

*INDEL stratification.* A Supplementary Table stratifies INDEL performance by insertion or deletion, length (1, 2–4, 5–10, > 10 bp), repeat context and platform, for every Pansoma model and both comparators [待补: comparators]:

- *Classes.* On HG008T chr1, 74.1% of truth INDELs lie in homopolymers of ≥ 7 bp and 11.1% in STRs. Population filtering removes 418 of the 635, capping recall at 0.342 for any caller using it (COLO829T: 0.382; R1.6); of the truth INDELs it keeps, 64.0–72.2% are called on every reference. On the HPRC graphs recall is lower for insertions than deletions (0.20–0.23 against 0.26–0.30), lowest for 2–4-bp changes (0.05–0.08) and near zero in STRs (0.00–0.04), where filtering removes most truth (57 of 65 on HG008T).
- *Platform error modes* (HG008T). On ONT, 15,797 of 17,023 false INDEL calls lie in homopolymers of ≥ 7 bp, where only 92 of 15,889 calls are true; at the validation best-F1 threshold, precision there rises to 0.25 (ONT) and 0.40 (PacBio). On Illumina, homopolymer precision (0.26) is not far below that outside repeats (0.35).
- *Allele fraction.* In the Fiber-seq dilutions, test-tuned best INDEL F1 is 0.049–0.074 at every purity (7–17 true calls), too few to resolve a trend.

INDEL accuracy is thus limited by overlap with population alleles, long-read homopolymer errors and, on long reads, the threshold, rather than by allele fraction.

### R3.m1 Parameter discrepancies between Methods and Supplement

**Comment（要点）。**

- 主 Methods（AF > 10 %、2 条 supporting reads、排除 MAPQ ≤ 5、80 epochs）和 Supplementary Methods（AF ≥ 0.08、≥ 3 个 ALT
  observation、mean ALT BQ ≥ 10、MAPQ ≥ 10、INDEL ≤ 50 bp、80–120 epochs）不一致。
- 上游（candidate）和下游（VCF、过滤、评估）步骤写得比网络部分简略得多。

**Response.**

We agree. One parameter table (Supplementary Table [待补]) now serves the main and Supplementary Methods, and every value is recorded in each run's manifest and checkpoint; re-running prediction and evaluation of 30 in-sample checkpoints on a second cluster reproduces the reported chr1 F1 within 0.002 (R3.9). Where the texts disagreed, the values below were used for every reported result; the others (including 80–120 epochs) are superseded. "MAPQ ≤ 5" and "MAPQ ≥ 10" were two different steps, now both stated:

- *Alignment:* vg Giraffe 1.66 for the GRCh38-only graph and, on HG008T, the HPRC v1.1 graph (d9 filter); on COLO829T v1.1, 1.64–1.66 (Fiber-seq, ONT) and ≤ 1.62 (Illumina), dated from the alignment files; 1.77 for the HPRC v2.1 graph (d46 filter).
- *Node selection:* alignments with MAPQ > 5; a node is a target when more than 5% of the reads aligned to it carry a mismatch, insertion or deletion on it.
- *Candidates:* reads with MAPQ ≤ 10 are dropped; at least 3 ALT reads, each ALT observation with base quality ≥ 10 (per observation, not a mean); allele fraction ≥ 0.08 for SNVs and INDELs; INDELs up to 50 bp; at most 800 reads per node (deterministic subsampling).
- *Tensors:* 8 channels × 200 rows × 101 columns (R1.m2).
- *Training:* chromosomes 2–22, 5% of their ~1-Mb blocks held out for validation; chromosome 1 is the test chromosome. SNV: 50.4 M-parameter ConvNeXt, 20 epochs, class weights pow 0.75. INDEL: 8.7 M parameters, 12 epochs, square-root class weights. Both: learning rate 1e-4, batch 1,024 per GPU. The table lists the GPU count of every model and the deviations of eight HPRC v1.1 models, which were also retrained with the matched recipes (R3.1); SNV models of the reference comparison were trained with three seeds (R3.9). The checkpoint is the epoch with the best thresholded validation somatic F1.
- *Normalisation:* every reported model z-scores base and mapping quality and the haplotype-count channel with its training statistics. Mapping quality comes from our own alignment and has the same scale across datasets; base quality depends on instrument and basecaller, so the code can instead standardise it per dataset, over all chromosomes of that dataset. [待补: whether the final models use this]
- *Calls:* PASS at the recall-0.9 validation threshold; for INDELs we also report the validation best-F1 threshold (R3.10) [待补: which of the two is the INDEL default]. An optional off-reference rescue calls SNVs on off-reference nodes by argmax after their path counts in the 21 central columns are raised to the all-haplotype code; we report every SNV result with and without the rescue (R3.4).
- *Population filtering resources (R2.6):* exact-allele matches. SNVs: gnomAD or CoLoRSdb at AF ≥ 1e-4, dbSNP records not flagged somatic, the 1000 Genomes panel. INDELs: gnomAD or CoLoRSdb at AF ≥ 0.01 in homopolymers of ≥ 8 bp, else ≥ 1e-4.
- *Evaluation:* rtg vcfeval against each truth set in its benchmark regions, and on the whole chromosome. Graph-coordinate calls go to a separate VCF that vcfeval does not score; a complete-output metric counts them as false positives (R3.4).

The Methods now describe these upstream and downstream steps in the same detail as the network.

### R3.m2 Comparator setup

**Comment（要点）。**

- 像描述 Pansoma 一样详细地记录比较器：软件版本、aligner 和版本、预训练 checkpoint、确切的 tumor-only 配置、PASS / filter 标准、参考资源、population / PoN 资源、命令行。
- Pansoma 的主要优势是 precision，而 precision 对过滤配置很敏感。

**Response.**

We agree. A new Supplementary Table documents both comparators at the level of detail given for Pansoma: software version, model checkpoint, tumor-only configuration, input BAM with its aligner and version, reference, population and panel-of-normals resources, PASS criteria and exact command lines.

*Software and models.* DeepSomatic 1.9.0 in tumor-only mode (`--model_type` WGS_TUMOR_ONLY, PACBIO_TUMOR_ONLY or ONT_TUMOR_ONLY; the HG008T ONT run, first made with the PacBio model, was repeated with the ONT model) and ClairS-TO 0.4.2 (HG008T) and 0.4.1 (COLO829T) with its platform models (Illumina `ilmn_ssrs` on HG008T and `ilmn` on COLO829T, PacBio `hifi_revio_ssrs`, ONT `ont_r10_dorado_sup_5khz_ssrs`). ClairS-TO has no model for the dorado SUP v5 basecalls of the HG008T ONT data, so its SUP 5 kHz model was used. All runs used the GRCh38 alignments released with the data (GIAB for HG008T, SMaHT for COLO829T), from the same reads as Pansoma's graph alignments, on chr1.

*Population resources and PASS criteria.* DeepSomatic tumor-only uses a population allele-frequency input channel (DeepVariant 1000 Genomes calls for Illumina, CoLoRSdb for long reads, AF > 0.05); its panel-of-normals filter applies only with `--use_default_pon_filtering=true`, which all our runs set. ClairS-TO tags non-somatic variants with four resources by default: gnomAD (AF ≥ 0.001), dbSNP 138 non-somatic sites, the GATK 1000 Genomes panel of normals and CoLoRSdb (AF ≥ 0.001), matching the allele for gnomAD and dbSNP and only the position for the other two. Our population filtering resources are the same four source files, with exact-allele matching throughout (R2.6). The two rules differ mainly for INDELs: on the whole truth sets (no benchmark regions), ClairS-TO's default rule tags 87.1% of the HG008T somatic INDEL alleles, 74.0% on COLO829T and 13.8% on HCC1395, against 62.7%, 53.0% and 7.9% for ours; for SNVs they agree within about one percentage point on six truth sets.

*Scoring.* The existing ClairS-TO and DeepSomatic call sets are re-scored with exactly Pansoma's protocol: the same truth, benchmark regions, chromosome and evaluation (rtg vcfeval 3.12.1, `--squash-ploidy --sample ALT,ALT`), in three settings: (1) the tool's own PASS; (2) its PASS after our population filtering; (3) all its records after our population filtering, ranked by its QUAL, as a precision–recall curve. Because precision depends on filtering, all three are reported; R1.6 and R2.3 summarise them.

*Training data.* Neither tool was trained on HG008; COLO829 is a ClairS-TO benchmark sample, not a training sample, and is not named in the DeepSomatic paper. DeepSomatic holds out chr1 in every model; the tumor-only models of its paper were trained on HCC1395, and from release 1.9.0 its WGS_TUMOR_ONLY model is trained on all of the paper's datasets. ClairS-TO's README states that all models exclude chr20, whereas its paper states chr1; its `ssrs` models were fine-tuned on HCC1937, HCC1954, H1437 and H2009. Pansoma's models in the main comparison are trained on chr2–22 of the same tumor; this asymmetry favours Pansoma and is stated beside every table, together with the cross-tumor results of R1.7.

### R3.m3 Composition of the pangenome reference

**Comment（要点）。**

- Document the HPRC release and version, the included haplotypes, the allele-frequency filtering, the treatment of structural variants and the graph construction.
- State whether benchmark individuals or their close relatives contributed to the graph.

**Response.**

The Methods now specify each reference.

- *HPRC v1.1 graph (d9 filter; release 1).* The Minigraph-Cactus graph with GRCh38 as the reference, in the frequency-filtered form distributed by the HPRC. Nodes and edges carried by fewer than 9 of the 90 haplotypes (about 10%) are removed; the two reference paths, GRCh38 and CHM13, are exempt.
  - The 90 haplotypes are 44 individuals × 2, plus GRCh38 and CHM13. The HPRC held out HG002, HG005 and NA19240.
  - The 44 individuals are African 23, American 16, East Asian 4, South Asian 1 and European 0.
- *HPRC v2.1 graph (d46 filter; release 2).* The same construction, filtered at 46 of 464 haplotypes.
  - The 464 haplotypes are 231 individuals × 2, plus GRCh38 and CHM13. The HPRC held out HG00272.
  - The individuals are African 70, East Asian 50, American 44, South Asian 36, European 29, and two without a population label (HG002, HG005).
  - We used the distributed graph and rebuilt its distance and minimizer indexes with vg 1.77.
- *GRCh38-only graph.* Built by vg autoindex (vg 1.66) from the GRCh38 no-alt analysis set that all truth sets use. Every node carries the same 17 paths, so its haplotype-count channel is constant.

*Structural variants.* The graphs keep the structural variants that pass the frequency filter, and reads align along them. Pansoma's candidates, however, are limited to indels of 1–50 bp; longer events are logged and not called.

*Benchmark individuals.* None of the cancer cell lines used here, nor their matched normals, is in either graph: HG008, COLO829, HCC1395, HCC1954, H1437 and H2009 are absent from both sample lists. These lines are anonymised, so we also tested for close relatives genetically. Of the germline SNVs of the HG008, COLO829 and HCC1395 normals that are rare in each release (on one or two of its non-reference haplotypes), no individual carries more than 8.2% in v1.1 or 2.9% in v2.1, whereas a first-degree relative would carry about half. In addition, because both graphs keep only non-reference nodes on at least 9 of 90 or 46 of 464 haplotypes, alleles private to one individual or a small family cannot enter either graph.

The HapMap-mixture donors are in the graphs: four of the six in v1.1 and all six in v2.1. R1.5 explains how we handle this.

The HG008 donor is of European ancestry, a group absent from v1.1 and present in v2.1. After population filtering, the results on the two graphs are close (R1.2).

### R3.m4 Terminology and language

**Comment（要点）。**

- 术语不一致，有若干语法问题。
- "panel of normals" 对 gnomAD、dbSNP、1000G、CoLoRSdb 来说范围太宽。
- 指出的几个具体句子要改。
- NPU 要给出全称。
- Pansoma 和 Pansoma-TO 混用；"Supplement Figure" 和 "Supplementary Figure" 混用。
- 最后做一次语言编辑。

**Response.**

We thank the reviewer and standardised the terminology throughout.

- The method is called Pansoma everywhere (not Pansoma-TO), and figures are cited as "Supplementary Figure".
- NPU is expanded at first use in the main text and in the Supplementary Methods ([待补: full name]).
- "Panel of normals" is no longer used as the collective name of gnomAD, dbSNP, CoLoRSdb and the GATK 1000 Genomes
  panel (only the last is a panel of normals in the usual sense). We now call the four "population filtering resources"
  (abbreviated PoN in figures and tables) and the step "population filtering"; each resource and its exact matching
  rule are listed in Methods (R2.6, R3.m1).
- Terms that recur across the paper (candidate, candidate site, off-reference node, graph-coordinate call, graph-absorbed
  allele, and the two operating thresholds, "the recall-0.9 validation threshold" and "the validation best-F1
  threshold") are defined once in Methods and used consistently (R2.m, item 8).
- The quoted sentences were corrected, for example "One key challenge in calling variants on graph alignments is that
  they differ fundamentally …", "Pansoma converts the variants …", "with the technical details described in Methods",
  and "read pileup" in Figure 5.
- The manuscript was edited for grammar and language by [待补: who did the language edit].
