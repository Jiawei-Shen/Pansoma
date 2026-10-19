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

We agree that "fully exploit the richness of pangenome data" overstated the original design and that node identifiers alone do not make a method pangenome-native; the phrase is removed and Pansoma is now described as pangenome-aware. The revision uses the graph more fully in two ways and tests its effect directly.

1. *Path-following tensors.* In the first version a tensor was confined to the node carrying the candidate, so short nodes truncated the read context (R1.m2). Each tensor row now follows its read's walk through the graph: the 101-column window around a candidate extends across every node the read traverses, on reference and alternative nodes alike, and each column records the base of the node the read passes there (channel 5). The context is therefore no longer bounded by node length: on chr1, the window of 57–76% of somatic SNV candidates and 92–97% of somatic INDEL candidates reaches beyond the candidate's own node on the HPRC graphs.

2. *A graph-native channel, with ablation tests.* Channel 6 gives, for every column, the number of haplotype paths through the node the read traverses, relative to the graph's haplotype count, which only a pangenome provides. Alone, the site's path count separates somatic from germline chr1 SNV tensors with an AUC of 0.831–0.865 on the HPRC v1.1 graph (d9 filter) and 0.885–0.915 on the HPRC v2.1 graph (d46 filter), against 0.500 on a GRCh38-only graph. Ablating channel 6 at inference (training-mean or all-haplotype value) in the 12 HPRC SNV models lowers somatic average precision from 0.199–0.352 to 0.044–0.087 (HG008T) and from 0.639–0.829 to 0.329–0.520 (COLO829T); setting it to the all-haplotype count lets 3.4–7.6 times as many germline tensors pass the threshold. Retraining without channel 6 on both HPRC graphs confirms this: chr1 somatic average precision falls from 0.199–0.353 to 0.065–0.089 on HG008T (COLO829T, v2.1: from 0.755–0.829 to 0.484–0.592), and SNV F1 after population filtering falls in 11 of 12 models (median 0.03, at most 0.31); on v2.1, where the models with channel 6 have three training seeds, all six models without it fall below that seed range by test-tuned best F1. Population filtering removes most of these extra germline calls, so channel 6 and the filter partly overlap. Channel 6 also suppresses somatic SNVs on off-reference nodes; we report every SNV result with and without an off-reference rescue (R3.3, R3.4).

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

*v1.1 versus v2.1.* The larger v2.1 graph produces a further 6–11% fewer germline candidates and keeps candidates for as many somatic alleles. Its gain is mainly before population filtering (SNV F1 +0.079 to +0.222, through precision); afterwards the two graphs are equal within the training-seed range on five datasets (three-seed means differ by −0.018 to +0.015), and indel F1 differs by at most 0.006. On COLO829T Illumina, v1.1 is 0.163 lower, through false calls in two segments of chr1 on every seed.

*Mapping ambiguity.* Graph alignment with Giraffe is as accurate as linear mapping and reduces reference bias (Sirén et al., 2021). On the same 10%-filtered HPRC v1.1 graph, Giraffe-based calling made 34% fewer small-variant errors than GRCh38, with the largest gains in difficult regions (challenging medically relevant genes, regions where GRCh38 is erroneous, L1HS repeats), and mapping to GRCh38 left 6.6–8.5 times more reads unmapped in VNTRs (Liao et al., 2023).

*Absorbed somatic alleles.* The concern holds for indels. Relative to GRCh38-only, v2.1 leaves 24–42% of HG008T and COLO829T truth indels without a candidate; at 85–96% of these loci the alternative reads follow graph nodes, entirely (43–59%) or with a residual edit (30–52%), mostly in repeats. Population filtering removes 76–86% of these indels anyway, so the net loss is 3.5–7.0% of truth indels. Population filtering itself removes 53–63% of truth indels, and comparable panels are the standard germline filter of tumor-only callers (Mutect2, ClairS-TO, DeepSomatic), so the graph adds little to a loss that every tumor-only caller shares (R1.4).

### R1.3 Difficult regions and ONT ambiguity

**Comment（要点）。**

- 在重复、paralog、结构复杂和低 mappability 区域，pangenome 比对仍可能有歧义，ONT 这类高噪声 reads 尤其如此；haplotype 越多，竞争比对位置可能越多。
- 汇总统计（更多 perfect reads、更少 mismatch / indel / soft clip）不能说明改善发生在 pangenome 应该起作用的地方。
- ONT 的 graph-coordinate predictions 更多：是分辨率提高，还是额外的歧义？
- 需要按区域的 benchmark 和错误分析。

**Response.**

We agree that aggregate statistics do not show where alignment improves, and added region- and call-level analyses.

1. *Alignment by region* (same reads, six datasets, SMaHT region classes). In easy regions the graphs give cleaner alignments at the same locus: on the HPRC v1.1 graph (d9 filter), HG008T PacBio mismatched bases per kb fall from 1.02 to 0.16, and 99.93–99.99% of long reads keep their placement. In extreme regions they add ambiguity: on v1.1, 6.26–7.41% of long reads lose MAPQ ≥ 30 and 3.05–4.93% gain it; on the HPRC v2.1 graph (d46 filter), aligned with a newer aligner release, 13.97–17.30% lose it. For short reads, lower MAPQ is not more misplacement: of simulated HG008-N Illumina reads of known origin, v1.1 and v2.1 misplace 0.98% and 0.93% of placed reads against 1.24% on a GRCh38-only graph (R3.6). Simulated ONT reads differ: v1.1 misplaces slightly fewer (0.61% against 0.77%), but v2.1 about twice as many (1.48%), almost all in extreme regions, in line with the more frequent confident disagreements of real ONT reads on v2.1 (HG008T: 1,953 against 1,019 reads on v1.1).

2. *Calls by region* (chr1 SNVs, PASS after population filtering, benchmark regions). On HG008T the arms are close in easy regions (F1 0.627–0.642 on GRCh38-only, 0.636–0.659 on the HPRC graphs). In extreme regions the GRCh38-only graph makes 351–2,210 PASS calls per dataset at precision 0.056–0.253, against 102–281 at 0.161–0.764 on the HPRC graphs. Placed ONT SNVs lose no more recall in extreme regions (HG008T, v2.1: 0.451 against 0.341–0.374), but on COLO829T their precision there is lower (0.452 against 0.693–0.764). [待补: per-stratum intervals]

3. *ONT graph-coordinate calls.* ONT makes the most graph-coordinate calls (chr1: HG008T 1,409 / 1,261 on v1.1 / v2.1, against 116–754 for Illumina and 211–244 for PacBio; COLO829T 1,675 / 4,566). Most are produced by the off-reference rescue; without it, 6–162 per run remain (R3.4). In the nine graph runs available (R3.3):
   - 2.8–19.1% can be written as a GRCh38 variant; after filtering this adds six true positives across the nine runs (0–2 per run; three distinct truth SNVs; R3.4).
   - In exclusive priority order, 49.8–91.4% match an HPRC allele and 4.2–45.4% are patient germline by assembly; 1–4 calls per HG008T run are tumor-only (4–24 by the assemblies alone). The HPRC-allele test is not specific: it tags 72 of 97 assembly-derived off-reference somatic SNVs (HG008T ONT, v2.1; R3.3).

   With the rescue, HG008T ONT complete-output F1 (no benchmark regions) is 0.342 / 0.354, below GRCh38-only (0.365); without it, 0.579 / 0.524, and all nine HPRC runs re-called locally are above GRCh38-only; inside the benchmark regions HG008T ONT on v2.1 stays 0.004 below (0.533 against 0.537). We report every SNV result with and without the rescue. The extra ONT calls reflect added ambiguity, mostly population and germline alleles, not resolution.

### R1.4 Candidate generation removes somatic variants already encoded in the graph

**Comment（要点）。**

- Candidates come only from mismatches / indels against the traversed node. A somatic allele that is already an alternative path, with tumor reads aligning perfectly to it, gives no candidate, tensor or call.
- Called a "fatal methodological issue".
- Supplementary Fig. 3 wrongly describes these false negatives as GRCh38 reference-bias artefacts.
- An allele on an alternative path is not necessarily germline: the same allele can be germline in one person and somatic in another.

**Response.**

The reviewer is correct: an allele on a graph path is not necessarily germline, and reads that align to it without edit give no candidate. We withdrew the reference-bias reading of Supplementary Fig. 3, state the mechanism as a limitation, and quantified the loss.

*Read level (HPRC v1.1 graph, d9 filter; patient assemblies; R3.7).* ALT reads of 2,478 HG008T truth alleles (201 SNVs, 2,277 INDELs; 26.8% of truth INDELs) and 497 COLO829T alleles (67 SNVs, 430 INDELs; 22.5%) align without edit to an existing path. Only 110 (4.4%) and 34 (6.8%) are the patient's germline allele; 2,085 (84.1%) and 430 (86.5%) are carried only by other HPRC individuals.

*Population filtering.* Tumor-only pipelines remove such alleles at filtering; ours (R2.6) tags 98.5–99.0% of these SNVs and 95.0–97.7% of these INDELs, against 3.6–10.2% of all truth SNVs and 53.0–62.7% of all truth INDELs. Absorbed and untagged remain 2 SNVs and 113 INDELs for HG008T (1.2% of truth INDELs) and 1 SNV and 10 INDELs for COLO829T (0.5%); on the HPRC v2.1 graph (d46 filter), 0.2–0.7% of truth INDELs per platform.

*Against a GRCh38-only graph* (same reads and code; chr1–22; tensor level). Net of gains, the HPRC graphs give a tensor to 0.2–0.6% fewer truth SNVs and 1.7–7.7% fewer truth INDELs that filtering would keep. Absorption is a small part: of the untagged INDELs lost on v1.1 (HG008T Illumina / PacBio / ONT: 461 / 550 / 276), 35 / 30 / 16 are absorbed by a strict spelling test; most others get no candidate without an exactly matching graph path or fall below the candidate allele-fraction cut. In the scored region (chr1, benchmark BED), 0–2 absorbed, untagged INDELs are lost per dataset and graph, and the GRCh38-only arm calls most of the untagged INDELs that the graphs lose (HG008T: 16–26 of 16–29; COLO829T: 4–12 of 10–13). Its final INDEL recall is similar: in a per-allele recall funnel, about two thirds of truth INDELs are population alleles, most of which it calls and filtering removes, and 24–27% of truth INDELs are called on every reference, slightly more on GRCh38-only (R1.6). But it also makes the patient's germline INDELs candidates and ranks somatic ones worse: its chr1 INDEL F1 after filtering is lower on all six datasets at the validation best-F1 threshold and by test-tuned best F1 (upper bound), though higher at the recall-0.9 validation threshold on the four long-read datasets, where every arm scores 0.133 or less (R3.1).

Linear tumor-only callers drop such alleles too: an earlier DeepSomatic tumor-only run labels 72 of 131 COLO829T Fiber-seq chr1 truth INDELs germline. [待补: ClairS-TO and DeepSomatic at the absorbed loci]

The classifier applies the population prior a second time through the haplotype-count channel: of the assembly-derived HG008T somatic SNVs on off-reference nodes, 57–79% per dataset have read support as somatic, yet the models call only 0–4 of these (R1.1, R3.3). The Discussion states both as limitations of population-based tumor-only calling.

### R1.5 HapMap mixture as a somatic benchmark, and data leakage

**Comment（要点）。**

- HapMap 的 "somatic-like" positives 是少数个体（minor donors）的 germline allele；Pansoma 的前提是图会压掉 population germline allele，两者冲突，这个 benchmark 可能系统性地对 pangenome caller 不利，也不是获得性低 VAF 突变的合适检验。
- 数据泄漏：HPRC v1.1 含 6 个 donor 中的 4 个（HG00438、HG02257、HG02486、HG02622）；companion paper（Fu, Xin, Miao et al., bioRxiv 2026.01.04.697580）提过这个问题，但没有排除这些 donor。

**Response.**

We agree with both points. The HAPMAP6 positives are germline alleles of the five minor donors (VAF 0.25–16.5%), not acquired mutations. Four of the six donors (HG00438, HG02257, HG02486, HG02622) are in the HPRC v1.1 graph (d9 filter), as our companion study noted, and all six, including HG002 and the 83.5% background HG005, are in the HPRC v2.1 graph (d46 filter). Apart from the GRCh38 and CHM13 reference paths, our graphs keep only nodes on about 10% or more of their haplotypes, so they absorb the positives that are common among those haplotypes whether or not the donors are in the graph; a donor's private alleles are filtered out either way. The donors' rare alleles are also found in CoLoRSdb, one of our population filtering resources, far more often than those of our cancer benchmarks' normals (chr20 SNVs: 14.5–25.9% against 4.0–7.9%). The benchmark therefore asks the caller to report alleles that its population prior, used by the trained model (R2.1), is designed to suppress.

We removed the HapMap mixture from the main benchmark and withdrew the HapMap-trained low-VAF model; no claim of the revised manuscript rests on HapMap results. Low-VAF performance is now measured on in-silico dilutions of COLO829T with COLO829BL on three platforms (nominal purity 75, 50 and 25%; measured 0.66 / 0.40 / 0.20 for Illumina and Fiber-seq, 0.70 / 0.44 / 0.22 for ONT). One model per platform and variant type is trained on chr2–22 of all four purities and tested on chr1 of each; the tumor itself is not held out.

- Fiber-seq SNV F1 after population filtering, at one recall-0.9 validation threshold pooled over all purities, is 0.779 / 0.687 / 0.665 / 0.720 at 100 / 75 / 50 / 25% purity (95% block-bootstrap interval at 50%: 0.648–0.682; 547–705 graph-coordinate calls per purity are not scored; R3.4). The share of truth SNVs that any model can reach falls from 0.821 to 0.711; on chr1, 16.1% of truth SNVs at full purity and 26.9% at 25% purity lie below the 0.08 candidate allele-fraction floor.
- The dip at 75 and 50% comes from the tumor's loss of heterozygosity: where the tumor lost the ALT copy of a germline heterozygous SNV, the admixed normal reads put that allele at (1 − purity)/2, inside the somatic allele-fraction range of those mixtures. The false positives at these two purities have a germline-like spectrum (cosine to the germline spectrum 0.788 and 0.768, against 0.536 at full purity).
- In the Illumina and ONT mixtures, tumor- and normal-derived reads differ in base quality, a cue that a real tumor lacks, so those scores are likely optimistic. For Illumina INDELs, removing the cue (per-read base-quality centring) lowers the dilutions' test-tuned best F1 (upper bound) to 0.133–0.146, below the 100% sample (0.182), as expected if the cue had inflated them [待补: Illumina SNVs and ONT after removing this cue].

[待补: dilutions of a second tumor not used in training (HCC1395)]

### R1.6 Low sensitivity

**Comment（要点）。**

- F1 的提升来自 precision；Supplementary Fig. 1 和 Tables 2–3 里 COLO829T、HG008T 的 SNV recall 低于 ClairS-TO 和 DeepSomatic-TO。
- INDEL recall 只有约 4–22 %（COLO829T）和 13–19 %（HG008T）。
- tumor-only 没有 matched normal 或第二阶段来补回漏掉的 variant，所以 sensitivity 很重要，F1 不能替代它。

**Response.**

We agree. Sensitivity is essential in tumor-only calling, and in the original submission Pansoma's SNV recall was below that of ClairS-TO and DeepSomatic-TO on both tumors, and its INDEL recall was low. The revision changes the operating point and reports recall separately from F1.

*Operating point.* PASS is the recall-0.9 validation threshold, set on held-out validation blocks before population filtering. Records down to p = 0.01 are kept as LowQual and population matches are tagged, not deleted, so users can move along the precision–recall curve, now shown for every caller.

*Where recall is lost* (chr1, benchmark regions, HPRC v1.1 graph with the d9 filter, models trained on chr2–22 of the same tumor; ranges over three platforms; HG008T INDEL candidates include partial allele matches; Supplementary Table [待补]):

| Recall | SNV HG008T | SNV COLO829T | INDEL HG008T | INDEL COLO829T |
|---|---|---|---|---|
| truth allele reaches a candidate | 0.944–0.951 | 0.812–0.816 | 0.539–0.702 | 0.397–0.458 |
| PASS | 0.855–0.862 | 0.705–0.736 | 0.359–0.431 | 0.321–0.420 |
| PASS, after population filtering | 0.826–0.831 | 0.686–0.717 | 0.228–0.256 | 0.237–0.267 |
| any threshold, after population filtering | 0.888–0.897 | 0.790–0.795 | 0.238–0.257 | 0.237–0.275 |
| bound for any caller using the same filter (whole chr1 truth) | 0.913 | 0.966 | 0.342 | 0.382 |

A per-allele funnel on all three references, GRCh38-only, v1.1 and the HPRC v2.1 graph (d46 filter), places each missed truth variant at its first loss. HG008T SNVs are called at 81.8%, 82.8% and 84.9%, the threshold being the largest loss on every reference; candidate-stage losses are mostly graph alleles and repeat absorption that filtering would tag anyway (R1.4). On COLO829T, 17–18% of truth SNVs get no candidate, mostly subclonal mutations below the 0.08 allele-fraction floor (15.1–16.1% of chr1 truth SNVs), whereas recall at allele fraction 0.05–0.2 is 0.850–0.966 once a candidate exists (R2.3, R2.4).

*INDELs.* The Abstract and Discussion now state this limitation (R3.10). 418 of the 635 HG008T chr1 truth INDELs (81 of 131 on COLO829T) are exact population alleles at the filtering frequencies, so no caller filtered this way exceeds a recall of 0.342 (0.382), and our post-filter recall at the recall-0.9 validation threshold is close to the curve's maximum. The graph also absorbs repeat-length INDELs (R1.4): a GRCh38-only graph gives more HG008T truth INDELs a candidate (0.75–0.87 against 0.54–0.70), but these are almost all population alleles, so 24–27% of truth INDELs are called on every reference and both tumors (funnel). At the validation best-F1 threshold, long-read INDEL F1 rises (HG008T PacBio, v1.1: 0.092 to 0.227), but recall falls to 0.073–0.158 on HG008T (v1.1) and 0.015–0.061 on COLO829T (v2.1), so we report both operating points. Both limits depend on the truth set: on HCC1395 (SEQC2), only 7.9% of truth INDELs are population-tagged and 0.76–0.85 reach a candidate (v2.1; HG008T: 0.43–0.50).

*Other tumors* (preliminary, HPRC v2.1 graph; R1.7). Without retraining, HG008T SNV models reach a post-filter recall of 0.640 / 0.701 / 0.589 on COLO829T (Illumina / Fiber-seq / ONT), against 0.715 / 0.707 / 0.692 for COLO829T's own models, and 0.773–0.930 on the lung lines H1437 and H2009; COLO829T models lose recall through threshold transfer (Illumina: 0.157 on HCC1395 and 0.296–0.305 on H1437 and H2009, at precision 0.908–0.973). A model retrained on four other lines (R1.7) reaches post-filter recall 0.666 / 0.703 on HG008T (Illumina / PacBio) and 0.527 / 0.447 on COLO829T, neither seen in training, at higher precision than the single-tumor models on HG008T [待补: ONT].

*Comparators*, re-scored under the identical protocol (R3.m2). At their own PASS, DeepSomatic tumor-only and ClairS-TO reach HG008T SNV recall 0.816–0.911 and 0.845–0.859 (precision 0.167–0.359), and after the same population filtering 0.806–0.889 and 0.831–0.845, against Pansoma's 0.826–0.831. At equal precision (Pansoma's curve at each comparator's post-filter precision) Pansoma's recall is 0.858–0.891 on HG008T, against the comparators' 0.806–0.889, and it is lower only than DeepSomatic on Illumina (HG008T 0.871 against 0.889; COLO829T, v2.1: 0.791 against 0.835). For indels the comparators' post-filter recall is 0.019–0.141 on HG008T and 0.015–0.183 on COLO829T (R2.3).

### R1.7 Training scale and tumor diversity

**Comment（要点）。**

- 约 200 M 参数的 ConvNeXt 主要用 COLO829（40–50 k truth variants）训练，数据量够不够？
- 报告确切参数量，以及每个平台、每种 variant type 的正 / 负训练样本数，并说明数据量撑得住模型容量。
- COLO829 是单一的 UV 主导的 melanoma，模型可能学到 tumor-specific 特征；DeepSomatic 用了多个细胞系并做了大量 augmentation。
- 把训练扩展到不同 mutational process 的肿瘤，并直接证明模型不是 tumor-specific。

**Response.**

We agree that one melanoma line was too narrow a basis, and we report the requested numbers.

*Model size.* The revised SNV network (ConvNeXt) has 50.4 M parameters, the INDEL network 8.7 M (Supplementary Methods); the original 199.6 M configuration overfitted (HG008T PacBio validation loss 0.51 to 2.16).

*Training examples.* Supplementary Table [待补: number] lists the examples per class for every model. On the pangenome graphs, SNV models see 7,098–7,274 (HG008T) or 33,472–34,110 (COLO829T) somatic, 256,583–381,159 germline and 67,735–1,620,829 non-variant examples; INDEL models see 3,121–3,761 (HG008T) and 971–1,141 (COLO829T) somatic examples. Every positive is an exact allele match, and each SNV positive is a different truth allele.

*Capacity.* With about 7,000 parameters per HG008T somatic SNV example, the network can fit its training set (final training loss 0.002–0.139), so generalisation rests on early stopping on held-out node blocks (best epochs 4–12 of 20). On HG008T PacBio (HPRC v2.1), F1 after population filtering rises from 0.574 with 10% of the training blocks to 0.597 with half and 0.607 with all (two replicates per fraction), and 8.7 M and 199.6 M networks reach 0.616 and 0.603, against 0.604–0.611 for 50.4 M over three seeds: accuracy gains little from more data and nothing from capacity.

*Tumor diversity.* COLO829T positives are 71% C>T, 65% at dipyrimidines (UV); HG008T positives are 41% C>T, 19% C>A and 8.4% C>G at TCW. We added two breast (HCC1395, HCC1954) and two lung adenocarcinoma lines (H1437, H2009), three platforms each, processed identically (R1.m1). They add APOBEC-like (HCC1395: 24–25% C>G) and tobacco-like processes (H1437, H2009: 42% and 38% C>A).

*Tumor specificity, a first direct test (preliminary).* We applied each HPRC v2.1 SNV model to chr1 of tumors it was not trained on, at the recall-0.9 validation threshold of its own training tumor (F1 after population filtering, benchmark BED, one training run; Supplementary Table [待补]):

- Models trained on HG008T (pancreatic; about 7,200 somatic examples) reach 0.708, 0.754 and 0.630 on COLO829T (Illumina, PacBio HiFi / Fiber-seq, ONT), 0.024–0.106 below COLO829T's own models (paired block-bootstrap intervals exclude zero), 0.689–0.790 on HCC1395 and 0.787–0.937 on the lung lines H1437 and H2009. Their true positives keep each target's spectrum (R3.5).
- Models trained on COLO829T reach 0.359, 0.557 and 0.487 on HG008T (its own models: 0.570, 0.605, 0.575). They lose partly through the transferred threshold (HCC1395 Illumina: recall 0.157 at precision 0.934; test-tuned best F1 0.810) and partly, it appears, through a UV preference: on HG008T they recall 0.493–0.841 of UV-context C>T but 0.016–0.129 of T>C (HG008T's own models: 0.905–0.925 and 0.613–0.645).
- On HCC1954 every model has low precision (0.22–0.37): 93–99% of the false positives are germline alleles of its matched normal, while its truth set shows no germline content, so the models let germline variants of this line through (R1.m1).

These models standardise base and mapping qualities by the training tumor's statistics. [待补: the same tests with base quality standardised per dataset (R3.m1)]

*Multi-tumor models* (preliminary: one run each, HPRC v2.1, Illumina and PacBio; Supplementary Table [待补]). Trained on HCC1395 (without chr1 and chr20) and the three CASTLE lines, a model reaches SNV F1 0.837 / 0.804 on HCC1395 chr1 (Illumina / PacBio), against at most 0.790 / 0.741 for a single-tumor model of another line, and 0.608 / 0.615 on the unseen HG008T (its own model: 0.570 / 0.603; higher precision, lower recall). On the unseen, UV-driven COLO829T it reaches 0.643 / 0.594 through lower recall; adding HG008T and COLO829T (without chr1) to training restores 0.742 / 0.774 (its own model: 0.783 / 0.779). [待补: ONT after a normalisation fix at prediction; leave-one-line-out of the CASTLE lines]

*Augmentation.* Joint COLO829T models add in-silico purity dilutions (75, 50, 25%; R1.5). [待补: Illumina SNVs and ONT after removing a read-source base-quality cue (R1.5; done for Illumina INDELs)]

### R1.8 Construction of the training set

**Comment（要点）。**

- 正样本：GRCh38 上的 truth 投影到 graph 坐标并按 exact allele 匹配；truth allele 本身就是 graph path、投影有歧义、graph 和 linear 表示不同（尤其 INDEL）时，正样本怎么定？
- graph-only 位置通过相邻的 1-bp SNV node 间接得到 GRCh38 坐标。
- 被吸收的位点永远成不了正样本（Supp. Fig. 3）。
- 要明确写出标签和投影过程，包括 graph 一侧是否做 normalisation，并说明正样本集合没有被系统性扭曲。
- 负样本：什么叫 "determinable"？germline 是否显式包含？artefact 怎么定义？VAF 和 purity 怎么取样？

**Response.**

We agree that the description was insufficient. Labelling is now described in full (Methods, Supplementary Note; each label manifest records the rule file's checksum).

*Positives and projection.* Each truth ALT is split, trimmed and, in repeats, expanded to every equivalent placement (up to 5,000). Each placement is converted along the graph's GRCh38 path into a node key, so no separate normalisation is needed on either side. A tensor is somatic only when its representative allele equals a truth key. Where graph and linear representations differ, matches with allele or read-haplotype overlap above 0.45 are recorded as partial and excluded from training, so every training positive is an exact match (in an earlier model version this raised HG008T Illumina's test-tuned best F1 (upper bound) from 0.623 to 0.665).

*Graph-only positions* are no longer placed through neighbouring nodes; without a truth match, branch-node tensors are excluded from training. The model thus never sees a positive off GRCh38 and scores a low haplotype path count as not somatic: of the assembly-derived HG008T somatic SNVs on off-reference nodes of the HPRC v1.1 graph (d9 filter), 0 of 112 (PacBio) and 1 of 89 (ONT) pass the threshold (HPRC v2.1 graph: 1–4 of 67–135), although 57–79% per dataset have read support as somatic. An optional off-reference rescue re-scores off-reference SNV tensors with the full haplotype count at the site and recovers 22–66 read-supported ones per dataset; it produces most graph-coordinate SNV PASS calls, and we report every SNV result with and without the rescue (R3.3, R3.4).

*Is the positive set distorted?* Alleles that are already graph paths produce no candidate, but on v1.1 95.0–99.0% of these absorbed truth alleles are population-tagged anyway; the graph-specific loss is 1.2% of HG008T and 0.5% of COLO829T truth INDELs (R1.4). In the training region, 94.8–95.9% (HG008T) and 89.5–90.5% (COLO829T) of truth SNVs stand behind a training example on either HPRC graph, against 98.5–98.7% and 90.1–91.4% on the GRCh38-only graph; the COLO829T shortfall is mostly subclonal SNVs below the candidate floor (15.1–16.1% of its chr1 truth SNVs lie below the 0.08 floor; R2.4). INDELs are affected more (42.6–66.7%, against 57.4–96.7% on the GRCh38-only graph). Supplementary Table [待补] gives every truth allele's status.

*Negatives.* Germline is an explicit class from the matched normal's truth (HG008-N and COLO829BL assembly calls; DeepVariant calls with QUAL ≥ 10 for the SEQC2 and CASTLE normals): 256,583–381,159 germline SNV examples per pangenome model. Non-variant examples (artefacts) are tensors on a GRCh38 node inside the confident region without a matching truth allele, including errors next to true variants. "Determinable" means exactly this: inside the confident region, at a unique GRCh38 position, and either an exact truth allele or none; all other tensors are ignored.

*VAF and purity* are not resampled: each line keeps its natural distribution (median somatic SNV allele fraction 0.53–0.54 for HG008T, 0.46–0.51 for COLO829T), and joint COLO829T models add purity dilutions (R1.5).

### R1.m1 Additional benchmarks

**Comment（要点）。**

- 在 SEQC2 的 HCC1395 tumor / normal 细胞系上做 benchmark。
- 在 DeepSomatic 研究里另外五个细胞系（CASTLE）上做 benchmark。

**Response.**

We added HCC1395 (SEQC2 truth v1.2.1, 41,072 somatic alleles) and three CASTLE lines, HCC1954, H1437 and H2009 (25,014, 84,744 and 161,091 somatic alleles on chr1–22), each sequenced with Illumina, PacBio HiFi and ONT and processed with the identical pipeline on the HPRC v2.1 graph (d46 filter).

*HCC1395 chr1 with models that never saw HCC1395* (preliminary: one training run, earlier quality normalization; 3,440 somatic SNVs in the benchmark region; F1 after population filtering at the recall-0.9 validation threshold of the training tumor; R1.7): models trained on HG008T reach 0.790 (Illumina), 0.698 (PacBio) and 0.689 (ONT); models trained on COLO829T (Fiber-seq for PacBio) reach 0.269, 0.741 and 0.621. Their test-tuned best F1 (upper bound) is 0.810, 0.819 and 0.691, so most of their shortfall comes from the transferred threshold, and part from a UV preference learned on COLO829T.

*CASTLE lines, chr1* (preliminary; same models and protocol): models trained on HG008T reach 0.787–0.937 on H1437 and H2009, and models trained on COLO829T 0.454–0.809, again losing mostly at the transferred threshold (test-tuned best F1 0.828–0.950). On HCC1954 every model has low precision (0.22–0.37): 93–99% of these false positives are germline alleles of its matched normal, while its truth set shows no germline content (2.4% common population alleles, no loss of heterozygosity in the normal), so this is germline leakage by our models, not a truth artefact. The H1437 truth, in contrast, includes about 30% germline alleles that its B-lymphoblastoid normal lost through loss of heterozygosity (chr9q, chr1:9–12 Mb), which explains the line's candidate coverage of about 70%; on chr1, our test chromosome, they are negligible (R2.2). A multi-tumor model trained on HCC1395 (without chr1 and chr20) and the three CASTLE lines reaches SNV F1 0.837 / 0.804 on HCC1395 chr1 (Illumina / PacBio), above every single-tumor model of another line, and generalises to HG008T, which it never saw (0.608 / 0.615; R1.7). [待补: ONT after a normalisation fix at prediction; leave-one-line-out of the CASTLE lines; ClairS-TO and DeepSomatic on these lines under the same protocol]

How we use these lines follows from how the comparators were trained.

- DeepSomatic's multi-cancer models were trained on HCC1395, HCC1937, HCC1954, H1437 and H2009 with chr1 held out, our test chromosome too, and the CASTLE truth sets of HCC1954, H1437 and H2009 were derived from DeepSomatic's calls (biased toward it, as its authors note). ClairS-TO's real-tumor model was fine-tuned on HCC1937, HCC1954, H1437 and H2009 (R3.m2). We therefore use these three lines mainly for training diversity and held-out-line tests of Pansoma, and report comparator results on them with this caveat.
- HCC1395 has a truth set independent of both tools (SEQC2). ClairS-TO excluded it from training; DeepSomatic held out its chr1 only, whereas our models above never saw HCC1395. HCC1395 chr1 is therefore our main external comparison, together with HG008T, on which neither tool was trained.
- The SEQC2 and CASTLE truth sets likely omit many INDELs in long homopolymers and STRs (15.5% of SEQC2 INDELs lie in homopolymers of ≥ 7 bp, against 74.9% for HG008T), so unmatched INDEL candidates in such repeats are not used as negatives on these lines (R2.2) and their INDEL results are reported separately.

We have not added HCC1937 and Hs578T, the other two lines of the DeepSomatic study. [待补: the reason]

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
- After population filtering (a per-tensor estimate), these two edits change F1 by only −0.066 to +0.007, except HG008T ONT on v1.1 with the training mean (0.596 to 0.145); replacing the plane with another tensor's lowers F1 to 0.212–0.551. Whether the plane improves the final calls is the retraining ablation of R1.1.
- The plane also suppresses somatic SNVs off GRCh38: of the assembly-derived HG008T somatic SNVs whose allele lies on an off-reference node, the models call 0 of 112 (PacBio) and 1 of 89 (ONT) on v1.1, against 61 and 49 with every node set to all haplotypes (v2.1, three platforms: 1–4 of 67–135, against 49–97), although 57–79% per dataset have read support as somatic (R1.1). An optional off-reference rescue recovers 22–66 of the read-supported ones per dataset at tensor level (R3.m1); we report every SNV result with and without the rescue (R3.4).

A schematic of a tensor crossing node boundaries is added as Supplementary Fig. [待补].

We agree that the tensors do not model bubbles explicitly: competing path support is represented through the reads' walks and per-node haplotype counts, not as a bubble-level object, and we state this as a limitation (R1.1).

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

*Recall difference to DeepSomatic.* We withdrew the truth-error explanation: on HG008T the tumor assembly shows that most graph-absorbed truth alleles are genuine somatic events (below; R1.4, R3.7). We re-scored the original DeepSomatic tumor-only and ClairS-TO call sets with our protocol (chr1, benchmark regions, the same population filter; R3.m2). On HG008T the SNV recall gap is then small: 0.806–0.889 (DeepSomatic) and 0.831–0.845 (ClairS-TO), against 0.826–0.831 for Pansoma on the HPRC v1.1 graph (d9 filter; 0.839–0.859 on the HPRC v2.1 graph, d46 filter), with models trained on chr2–22 of the same tumor. On COLO829T the comparators reach 0.555–0.835 and Pansoma 0.686–0.717. Pansoma's F1 is the highest on five of the six datasets (on all six with v2.1). The largest remaining gap, DeepSomatic on COLO829T Illumina (0.835 against 0.711–0.715), lies both below an allele fraction of 0.08, where Pansoma builds no candidate (30% of these truth SNVs recovered against 6%), and above it (94% against 84%). For indels Pansoma's recall is the highest on both tumors (0.228–0.260 against 0.015–0.183). A per-allele funnel (R1.6) shows where Pansoma loses the rest. On HG008T the largest loss on every reference is the threshold (6.8–10.9% of truth SNVs, against 2.1–5.2% without a candidate and 3.0–5.2% removed by filtering); 61 of the 702 truth SNVs are population alleles that any caller using the same filter loses. On COLO829T (GRCh38-only graph and HPRC v2.1 graph, d46 filter), 17–18% of truth SNVs get no candidate, mainly because 15.1–16.1% have a tumor VAF below the 0.08 candidate floor (R2.4).

*The two categories.* "Variant node not ranked #1" and "first aligned node mismatch" belonged to the earlier pipeline, which we replaced. Every truth allele without a candidate now gets one read-level class, on six datasets (later losses: R1.6): no read with MAPQ > 10; fewer than three ALT reads or an allele fraction below the threshold; the ALT is an existing graph allele; repeat absorption or respelling; an indel longer than 50 bp. For example, of 2,976 HG008T PacBio truth INDELs without a candidate (chr1–22), 2,739 (92.0%) are fully or partly encoded in the graph.

*The normal at missed sites on graph nodes.* Instead of realigning normal reads, we used the haplotype-resolved matched-normal assemblies (HG008-N v6.2; COLO829BL), which show the normal's sequence directly. Of the truth alleles whose tumor reads align without edit to an existing v1.1 graph path (HG008T 2,478; COLO829T 497), only 110 (4.4%) and 34 (6.8%) are carried by the matched normal; 2,085 (84.1%) and 430 (86.5%) are absent from both normal haplotypes and lie on nodes walked by other HPRC individuals. At the 1,935 HG008T INDEL loci absent from the normal, the HG008-T tumor assembly carries a sequence found in neither normal haplotype at 1,931 and exactly the truth allele at 1,652, and the tumor's Illumina reads show the expected allele at 96.2% of them (95.7% at random non-absorbed truth INDELs). On HG008T these are real somatic events hidden by the graph, and the manuscript now says so; most also coincide with population alleles that the filter removes on any reference (R1.4): in the chr1 benchmark regions, the GRCh38-only arm calls at most one of the 38–114 HG008T truth INDELs that v1.1 absorbs. COLO829T has no tumor assembly, and up to 15.8% of its absorbed INDELs may be germline-like.

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

We compared three references built from the same reads, tensor code and labels: a GRCh38-only graph and the HPRC v1.1
and v2.1 graphs (R1.2, R3.1). The precision gain comes mainly from fewer germline candidates, and it is largest in
difficult regions and where germline variants cluster.

*Fewer germline candidates.* At candidate generation the HPRC v1.1 and v2.1 graphs produce 84–90% and 85–91% fewer
germline candidates than the GRCh38-only graph, because reads carrying a germline allele already in the graph follow its
path without a mismatch. Before population filtering, the HG008T models call 117,205–144,251 germline tensors somatic on
GRCh38-only, against 2,528–5,147 on the HPRC graphs. After filtering, in the benchmark regions, the HPRC graphs lead by
+0.028 to +0.292 SNV F1 (means over three training seeds; single-run paired block-bootstrap intervals exclude zero in all twelve comparisons; R3.9).

*Difficult regions.* We stratified chr1 by the SMaHT easy, difficult and extreme regions. They are defined on GRCh38
from the 1000 Genomes accessibility mask and a pangenome mask, not from any sample, so they apply to HG008T as well as
COLO829T (R1.3). The gain grows with difficulty: in extreme regions the HPRC graphs raise SNV precision from 0.056–0.253
to 0.161–0.764 (F1 +0.099 to +0.531 in all nine graph–dataset comparisons; COLO829T Illumina: 39 instead of 2,087 false
calls); in difficult regions F1 rises in eight of nine (−0.005 to +0.259); in easy regions HG008T changes by −0.003 to
+0.020.

*Clustered germline variation.* We grouped calls by the number of the normal's germline variants within 50 bp (0, 1,
2+). For HG008T SNVs the F1 gain is largest with two or more nearby (0.107–0.204, against 0.022–0.087 with none), through
precision; recall stays 0.256–0.308 on every reference. The reads show why. At truth SNVs with two or more germline
variants nearby, reads differ from the reference at 1.8–3.0 positions each on GRCh38-only, against 0.4–1.8 on the HPRC
graphs, and the share of reads with at least two differences falls from 53–74% to 5–30%; with none nearby the
references give the same reads, and mapping quality is the same throughout (mean 58–60). The GRCh38-only false calls at
such sites (125–236 per HG008T platform) are 82–97% germline alleles of the normal or calls within 10 bp of one; 86–94%
of them are not called on the HPRC graphs, and 41–87% are not even candidates there. COLO829T agrees (94–97% of 23–49
false calls removed).

*Limits.* For indels the graphs lose where germline variants cluster: only 5.9–21.6% of the truth indels in that
stratum get an exact candidate, against 62.7–94.1% on GRCh38-only (R1.4). The false calls left after filtering are
mostly the patient's rare germline alleles, which neither the graph nor the filter contains (R2.6), and v2.1 adds no
precision over v1.1 after filtering (R1.2).

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

We agree that lines 292–294 overstated the effect, and have rewritten them. Supplementary Fig. 4 is replaced by counts
and proportions for the same reads on four references (the GRCh38 BAM, a GRCh38-only graph, and the HPRC v1.1 and v2.1
graphs), for both tumors and all three platforms (Supplementary Table and Figure below). In summary:

- *Illumina:* full-length edit-free alignments of the same reads rise by about 10 percentage points on the HPRC graphs
  (HG008T 69% to 79%, COLO829T 71% to 81–82%). The GRCh38 BAM and the GRCh38-only graph give the same rate (HG008T
  68.5% and 68.8%, COLO829T 70.6% and 70.7%), so the gain comes from the graph's content, not from the aligner.
- *Long reads* are rarely edit-free end to end. On HPRC v1.1 their mismatched bases per kb fall (HiFi 1.57 to 0.63, ONT
  5.22 to 4.26), while indel rates, soft clipping and unmapped reads change little. HPRC v2.1 aligns nearly every HiFi
  read, but places the additional ones in difficult regions at lower mapping quality.
- *More edit-free alignments are not more correct alignments.* On reads simulated from the HG008-N assembly, the share
  placed correctly is the same on every reference (Illumina 98.6–98.7%), while confident misplacements fall for short
  reads and HiFi; ONT reads are the exception on v2.1 (R3.6).

*Published benchmarks agree.* On the 10%-filtered HPRC release 1 graph (the filter of our v1.1 graph), 78.1–78.9% of
Illumina reads aligned perfectly, against 68.7% with BWA-MEM on GRCh38, at a slight cost in mapping confidence (Hickey et
al., 2024), and Giraffe-based calling made 34% fewer small-variant errors than GRCh38-based workflows (Liao et al.,
2023). Giraffe is as accurate as linear mappers on simulated short reads (Sirén et al., 2021) and, on the HPRC release
2 graph with the same d46 filter as our v2.1 graph, mapped simulated HiFi reads with the highest true-positive rate and
the fewest confident errors, with similar or slightly better variant calls than linear mappers (Chang et al., 2025,
preprint).

The revised text says that the graph mainly removes edits at alleles that the sample shares with HPRC haplotypes.

**Supplementary (R2.7).**

![Illumina, the same reads on four references](supplementary/R2.7/fig_R2.7a.png)

*a,* Illumina reads in 48 GRCh38 windows (24 random, 24 at segmental duplications): the share aligned end to end with no
mismatch, indel or soft clip.

![Long reads, mismatched bases and indel events per kb](supplementary/R2.7/fig_R2.7b.png)

*b,* Long reads (a 10% genome-wide sample of read names, the same reads on every graph): mismatched bases and indel
events per aligned kb. The higher HiFi and Fiber-seq rates on v2.1 come from the reads it newly aligns (unmapped 2.6–3.4%
on v1.1, 0.01% on v2.1), mostly in difficult regions.

![Simulated reads, correct placements and confident misplacements](supplementary/R2.7/fig_R2.7c.png)

*c,* Reads simulated from the HG008-N assembly (Illumina 5,000,000; HiFi-like and ONT-like 40,000 each): the share placed
within 50 bp of the true origin, and the misplaced reads with MAPQ ≥ 30 per 1,000 reads.

Supplementary Table (same samples as *a* and *b*; MAPQ as reported by each aligner, so not comparable between the BAM and
the graphs):

| dataset | reference | unmapped % | MAPQ 60 % | soft-clipped reads % | edit-free reads % | mismatched bases / kb | indel events / kb |
|---|---|---:|---:|---:|---:|---:|---:|
| HG008T Illumina | GRCh38 BAM (BWA-MEM) | 0.00 | 88.2 | 2.8 | 68.50 | 4.73 | 0.21 |
| HG008T Illumina | GRCh38-only graph | 0.00 | 86.9 | 2.2 | 68.85 | 4.63 | 0.21 |
| HG008T Illumina | HPRC v1.1 graph (d9) | 0.00 | 86.2 | 1.9 | 79.15 | 3.69 | 0.05 |
| HG008T Illumina | HPRC v2.1 graph (d46) | 0.00 | 86.3 | 1.9 | 79.24 | 3.67 | 0.05 |
| COLO829T Illumina | GRCh38 BAM (BWA-MEM) | 0.00 | 86.5 | 3.3 | 70.57 | 3.93 | 0.20 |
| COLO829T Illumina | GRCh38-only graph | 0.00 | 85.6 | 2.6 | 70.72 | 3.95 | 0.21 |
| COLO829T Illumina | HPRC v1.1 graph (d9) | 0.00 | 85.1 | 2.6 | 81.25 | 2.97 | 0.04 |
| COLO829T Illumina | HPRC v2.1 graph (d46) | 0.00 | 85.5 | 2.4 | 81.61 | 2.94 | 0.04 |
| HG008T PacBio HiFi | GRCh38-only graph | 2.79 | 95.9 | 2.6 | 0.15 | 1.57 | 1.58 |
| HG008T PacBio HiFi | HPRC v1.1 graph (d9) | 2.64 | 95.9 | 2.0 | 1.39 | 0.63 | 1.39 |
| HG008T PacBio HiFi | HPRC v2.1 graph (d46) | 0.01 | 83.3 | 3.8 | 1.42 | 1.72 | 1.56 |
| COLO829T Fiber-seq | GRCh38-only graph | 3.37 | 93.1 | 3.3 | 0.18 | 1.80 | 1.82 |
| COLO829T Fiber-seq | HPRC v1.1 graph (d9) | 3.23 | 93.3 | 2.6 | 1.50 | 0.84 | 1.62 |
| COLO829T Fiber-seq | HPRC v2.1 graph (d46) | 0.01 | 79.3 | 5.0 | 1.54 | 2.19 | 1.87 |
| HG008T ONT | GRCh38-only graph | 7.86 | 95.2 | 30.4 | 1.32 | 5.22 | 4.27 |
| HG008T ONT | HPRC v1.1 graph (d9) | 7.75 | 95.0 | 29.9 | 2.54 | 4.26 | 4.12 |
| HG008T ONT | HPRC v2.1 graph (d46) | 8.01 | 94.2 | 29.6 | 2.50 | 4.35 | 4.33 |
| COLO829T ONT | GRCh38-only graph | 13.53 | 93.9 | 99.0 | 0.00 | 9.16 | 8.59 |
| COLO829T ONT | HPRC v1.1 graph (d9) | 13.23 | 93.5 | 99.0 | 0.01 | 8.36 | 8.55 |
| COLO829T ONT | HPRC v2.1 graph (d46) | 13.61 | 89.5 | 99.0 | 0.01 | 8.93 | 9.12 |

Data: `supplementary/R2.7/alignment_stats.tsv`, `supplementary/R2.7/simulated_placement.tsv` (made by
`supplementary/R2.7/make_figures.py` from `experiments/E11a_alignment_stats/results/*.json` and
`experiments/E11b_simulated_alignment/results/summary.tsv`).

*Methods.*
1. *Alignments.* Each dataset's tumor reads were aligned with vg Giraffe to the GRCh38-only graph and to the HPRC v1.1
   (d9) and v2.1 (d46) graphs (versions as in R3.m1); the GRCh38 BAMs are the data providers' BWA-MEM 0.7.17 alignments
   (HG008T: GIAB; COLO829T: SMaHT, Sentieon 202308.01 with duplicate marking, indel realignment and base recalibration).
2. *Same reads.* Reads are selected by a hash of the read name, so every reference sees the same reads. Illumina: 10% of
   the reads in the 48 windows (± 1 kb) that any reference places there; long reads: 10% of all read names, whole files
   scanned.
3. *Measures.* Edit-free: mapped with no mismatch, indel or soft clip. Mismatched bases and indel events per kb of
   aligned bases, from each alignment's edits; unmapped, soft-clipped and MAPQ shares per read.
4. *Simulation.* Reads drawn from one-to-one segments of the HG008-N assembly on chr1 and chr20, with known GRCh38
   origin (Illumina 2 × 150 bp, 0.2% substitutions; HiFi-like ~15 kb, 0.1% errors; ONT-like with the pbsim3 R10.4
   model), aligned to each graph and projected to GRCh38. Correct: the right chromosome and a read end within 50 bp of
   the truth.

*References.*
- Sirén J. et al. Pangenomics enables genotyping of known structural variants in 5202 diverse genomes. *Science* 374,
  abg8871 (2021).
- Liao W.-W. et al. A draft human pangenome reference. *Nature* 617, 312–324 (2023).
- Hickey G. et al. Pangenome graph construction from genome alignments with Minigraph-Cactus. *Nat. Biotechnol.* 42,
  663–673 (2024).
- Chang X. et al. Rapid, accurate long- and short-read mapping to large pangenome graphs with vg Giraffe. Research
  Square preprint, doi:10.21203/rs.3.rs-7747403/v1 (2025; not peer reviewed).

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

We ran this ablation on HG008T and COLO829T, three platforms each. The same reads were aligned with vg Giraffe to a GRCh38-only graph built from the same FASTA, the HPRC v1.1 graph (90 haplotypes; d9 filter) and the HPRC v2.1 graph (464 haplotypes; d46 filter); everything downstream is identical (candidates, tensors, labels, network, training on chr2–22, threshold, population filtering), and models are tested on chr1 of the same tumor. SNV results are three-seed means after population filtering, in the benchmark regions.

1. *Both HPRC graphs outperform GRCh38-only.* SNV F1 is 0.586–0.768 (v1.1) and 0.575–0.780 (v2.1), against 0.454–0.691, with higher precision and recall on all six datasets (F1 +0.028 to +0.292; all paired block-bootstrap intervals exclude zero, R3.9). INDEL F1 at the validation best-F1 threshold is higher in all twelve comparisons; at the recall-0.9 threshold, GRCh38-only is higher on the four long-read datasets, where all arms are poor (R3.10).

2. *v2.1 against v1.1: a gain before population filtering, little after it.* Before filtering, v2.1 is higher on all six datasets (SNV F1 +0.079 to +0.222), mostly through precision (+0.089 to +0.280; recall within 0.016). After filtering, which removes most of the same germline variants, the two graphs are equal within the seed range on five datasets (−0.018 to +0.015); on COLO829T Illumina, v1.1 stays 0.163 lower through false calls in two segments of chr1 (R1.2).

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
2. *Placement.* Via the GRCh38 nodes flanking their bubble, 2.2–19.1% per run become a GRCh38 allele, adding six true positives after population filtering across the twelve runs (0–2 per run; three distinct truth SNVs). Of 89–1,718 calls per run placed only to a GRCh38 interval, one matches a missed truth variant (R3.4).
3. *Assemblies.* We aligned each call's ALT context to the HG008-T tumor assembly and both HG008-N haplotypes (COLO829T: the two COLO829BL haplotypes only, so somatic status cannot be confirmed). In controls (150 per group and run), 145–149 placed true positives are tumor-only on HG008T (COLO829T: 148–150 absent from the normal), and 135–146 (HG008T) and 110–131 (COLO829T) placed false positives after filtering are in the normal. In exclusive classes, the graph-coordinate calls are population alleles (ALT context spelled by GRCh38 or an HPRC haplotype or walk) in 49.8–91.4%, germline in 4.2–45.4%, and tumor-only in 0.2–0.9% on HG008T (COLO829T: 2.4–5.8% absent from the normal). By the assemblies alone, 1.1–5.2% per HG008T run are tumor-only, most of them also population alleles.
4. *Graph-native population filter.* It tags 0–1 of 150 placed true positives per run, but also 74.2–88.3% of the assembly-derived off-reference somatic SNVs (all chromosomes, six HG008T sets), and 20 of the 26 PASS calls of such SNVs on chr1 (six runs). It describes calls and cannot certify them.
5. *Cross-platform.* On v2.1, 81 calls (same node, offset and allele) are PASS on all three HG008T platforms, 68–69 of them population alleles.

So on HG008T at most 1.1–5.2% of graph-coordinate calls per run have assembly support as somatic; the rest are mostly population or germline alleles that a GRCh38-based filter cannot reach, and counting them as false positives (R3.4) is close to their true status. Because the rescue also recovers read-supported off-reference somatic SNVs (point 1), we report every SNV result with and without the rescue; that training with assembly-derived off-reference labels would address the cause is stated as a limitation.

### R3.4 Metrics exclude graph-coordinate calls; "about 2%"; CHM13

**Comment（要点）。**

- TP / FP / FN are computed only for calls that project to GRCh38, so precision, recall and F1 describe a subset, not the complete output.
- The graph-coordinate fraction depends on the platform and is large: above 10 %, and sometimes above 20 % for ONT (Supplementary Table 6). This contradicts the Discussion's "about 2%".
- Are graph-coordinate calls truly graph-dependent, or just not projectable to GRCh38? What fraction can be represented on T2T-CHM13, and what remains graph-specific?

**Response.**

We agree, and the revised evaluation covers the complete output.

1. *"About 2%" was wrong and has been removed.* On chr1, after population filtering and without benchmark regions, graph-coordinate calls are 8.1–52.4% of the PASS SNV output on the HPRC v1.1 graph (d9 filter) and 11.4–68.0% on the HPRC v2.1 graph (d46 filter), highest for ONT throughout. Every benchmark table now reports this fraction, replacing Supplementary Table 6.
2. *Complete-output metrics.* Every run reports F1 for the GRCh38-evaluable calls and for the complete output, counting every graph-coordinate PASS call as a false positive (a worst case). HG008T SNV, chr1, after population filtering, no benchmark regions, HPRC v1.1 / v2.1 / GRCh38-only graph:

   | | Illumina | PacBio | ONT |
   |---|---|---|---|
   | evaluable calls | 0.577 / 0.561 / 0.197 | 0.614 / 0.600 / 0.359 | 0.586 / 0.563 / 0.365 |
   | complete output | 0.546 / 0.414 / 0.197 | 0.553 / 0.533 / 0.359 | 0.342 / 0.354 / 0.365 |

   Under this worst case the HPRC graphs fall below the GRCh38-only graph in four of the twelve SNV comparisons, all long-read: HG008T ONT on both graphs, COLO829T ONT on v2.1 (0.369 against 0.418) and COLO829T Fiber-seq on v1.1 (0.572 against 0.680).
3. *Graph-specific or merely unprojected?* In the twelve graph runs, the GRCh38 nodes flanking each call's bubble give a GRCh38 allele for 2.2–19.1% of calls. After population filtering this adds six true positives across the twelve runs (0–2 per run; three distinct truth SNVs), and of the 89–1,718 calls per run placed only to a GRCh38 interval, one matches a missed truth variant. The rest have no GRCh38 allele (inserted sequence, the GRCh38 base on a branch node, no consistent anchors); most are population or germline alleles, and on HG008T only 1.1–5.2% per run have tumor-only support in the assemblies (R3.3). Projection moves complete-output F1 by at most 0.009. A Supplementary Table lists each call's node, offset, class, GRCh38 interval and assembly class.
4. *Origin and trade-off.* Most graph-coordinate SNV PASS calls are produced by an off-reference rescue step that raises the haplotype count on nodes off GRCh38 (R3.3); without it, 6–162 per run remain. The rescue recovers off-reference somatic SNVs that rtg cannot credit: at the assembly-derived HG008T sites, it calls 22–66 read-supported somatic SNVs per dataset (0–4 without it; tensor level), and 63–87% of its calls there are read-supported (R3.3). On chr1, however, that is 2–8 such SNVs per HG008T run against 210–1,410 graph-coordinate PASS calls, and the rescue causes the long-read losses above. Without it, the placed calls are unchanged and all twelve HPRC runs are above the GRCh38-only graph without benchmark regions (complete-output F1 0.524–0.780 against 0.197–0.680); inside the benchmark regions, HG008T ONT on v2.1 stays 0.004 below (0.533 against 0.537). We report every SNV result with and without the rescue.
5. *CHM13.* We did not compute a CHM13 projection. CHM13 is an unfiltered reference path in both graphs, but our truth sets are on GRCh38, so CHM13 coordinates would not change any score; point 3 separates the calls that GRCh38 can represent from those it cannot.

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

1. *What the gain is.* For the same HG008T Illumina reads (48 GRCh38 windows; R2.7), 79.36% align without any edit on the HPRC v1.1 graph (d9 filter) against 69.12% on the GRCh38-only graph, and 99.74% of the reads that become edit-free keep their GRCh38 position. The graph removes edits at alleles that the sample shares with HPRC haplotypes; this alone is not evidence of more correct placement.

2. *MAPQ, same reads.* Over whole files, the share of mapped reads with MAPQ 60 on the GRCh38-only, v1.1 and v2.1 (d46 filter) graphs is 95.88%, 95.88% and 83.32% for HG008T PacBio and 92.38%, 92.63% and 76.87% for COLO829T Fiber-seq. The v2.1 drop is smaller on ONT and absent in Illumina windows. v2.1 was aligned with a newer Giraffe version, so graph and aligner effects are not separated. By region, both HPRC graphs add some ambiguity in SMaHT extreme regions and segmental duplications (R1.3).

3. *Correctness with known origin.* We simulated reads from the HG008-N v6.2 diploid assembly (in neither HPRC graph) on chr1 and chr20, where it aligns one-to-one to GRCh38. For 5 million Illumina reads, v1.1 and v2.1 give more edit-free alignments than the GRCh38-only graph (72.94% and 73.06% against 63.69%) and fewer misplaced reads (0.983% and 0.927% of placed reads against 1.235%). "Perfect but wrong" alignments do not increase (33,204 and 32,025 against 33,847, nearly all in SMaHT extreme regions), and the error rate among placements with MAPQ ≥ 30 falls from 0.093% to 0.019% and 0.013%; MAPQ is over-confident on all three. The gain concentrates in reads spanning a structural variant (misplaced 51.0% and 46.0% against 71.7%), although fewer of them are placed correctly (21.2% and 21.7% against 25.1%); far from any HG008-N variant the HPRC graphs misplace slightly more (0.807% against 0.697%). HiFi-like reads go the same way (misplaced 0.328% and 0.383% against 0.495%). Noisy ONT-like reads do not: v1.1 misplaces slightly fewer (0.614% against 0.773%), but v2.1 about twice as many (1.483%; at MAPQ ≥ 30, 0.764% against 0.433%), almost all in SMaHT extreme regions, as the confident disagreements of real ONT reads on v2.1 suggested (R1.3); its newer Giraffe release may contribute. The hardest regions are under-represented, and the v1.1 arm used Giraffe 1.65.

4. *Phased haplotypes on real reads.* We compare the graph placement of real HG008-T reads in difficult regions with their placement on the phased HG008-N assembly [待补].

5. *Completeness versus multi-haplotype modelling.* In the reference ablation (R3.1), the HPRC graphs raise chr1 SNV F1 after population filtering by +0.028 to +0.292 on all six datasets (means over three training seeds). Both HPRC graphs contain CHM13 as an unfiltered reference path, so this gain, like the simulation above, combines reference completeness and haplotype diversity. We did not add a T2T-CHM13 arm, because our truth sets are on GRCh38, and state this as a limitation.

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
  without edits to an existing graph path. Only 110 (4.4%) are carried by the matched normal, and in 109 of them the
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

We keep the scope to tumor-only somatic calling in cancer samples. The title, abstract, Results and Discussion now state
this consistently, and broader claims about "somatic variation" in general were removed.

The method is not built for mosaicism in normal tissue, and the revised Methods say why:

- A candidate needs at least three ALT reads and an allele fraction of at least 0.08 for SNVs and INDELs, so variants at
  the much lower allele fractions typical of mosaicism in normal tissue are not evaluated.
- The training positives are tumor variants, mostly clonal: the somatic SNV training tensors of the single-sample
  models have a median allele fraction of 0.53–0.54 (HG008T) and 0.46–0.51 (COLO829T), and only 0.4–1.2% (HG008T) and
  6.5–7.3% (COLO829T) lie below 0.20.
- The classifier is trained to separate a tumor's somatic variants from the germline variants and sequencing errors of
  the same sample.

Validating on an independently characterised mosaicism dataset such as the SMaHT tissues would need a low-allele-fraction
candidate configuration, new training labels and new models. We describe this in the Discussion as future work rather
than claim it.

Within the cancer scope we extend the validation beyond HG008T and COLO829T to further cancer cell lines (R1.7,
R1.m1) [待补: the lines with final results], and the HapMap mixture, whose positives are germline alleles of other individuals, is no longer used as evidence of
somatic performance (R1.5).

### R3.9 Model stability and statistical uncertainty

**Comment（要点）。**

- 结果都是点估计；用多 seed 重训、基因组区间的 block bootstrap 置信区间或染色体间方差评估稳定性。
- 特别是：多 seed 重训，报告 held-out 染色体和独立 HG008T benchmark 上 precision / recall / F1 的分布。
- 用染色体或 block 级重抽样给 Pansoma 与比较器之间差值的置信区间。
- 报告对固定 0.5 阈值的敏感性。

**Response.**

We agree that single point estimates are not enough, and we now report uncertainty for every headline number.

1. *Repeated training.* Training is deterministic for a given seed and dataset (three models retrained with the same seed reproduced their chr1 calls exactly), and inference is stable across hardware: re-running 30 in-sample checkpoints on another cluster reproduced their chr1 F1 after population filtering within 0.002. The variance we report therefore comes from the seed, which also redraws the validation blocks. Every in-sample SNV model of the three-reference comparison (six datasets, three references) was retrained with two further seeds. Across the three seeds of a model, chr1 somatic average precision varies by at most 0.046 on the HPRC graphs, but F1 after filtering at the PASS threshold by 0.003–0.080 (median 0.022; up to 0.144 on GRCh38-only), because the recall-0.9 validation threshold moves with the redrawn validation blocks; the test-tuned best F1 varies by at most 0.065. Differences below about 0.05 between single runs are therefore within training variance. We report the mean and range over the three seeds; both HPRC graphs stay above GRCh38-only in 35 of 36 seed-matched comparisons (R3.1) [待补: further seeds for INDEL models and held-out tumors]. HG008T held out from training replaces the original COLO829T-to-HG008T benchmark; a first single-seed test on the HPRC v2.1 graph (d46 filter), with models that still standardised base and mapping qualities by the training tumor's statistics, loses 0.046–0.211 F1 against HG008T's own models (R1.7).

2. *Block bootstrap.* 95% intervals come from 2,000 resamples of the 1-Mb blocks of chr1 over rtg's true-positive, false-positive and false-negative records; a paired difference uses the same resampled blocks for both members. The test sets are small (697 somatic SNVs on HG008T chr1, 131 somatic INDELs on COLO829T chr1): SNV F1 intervals are 0.033–0.074 wide on the HPRC graphs and up to 0.341 on the GRCh38-only graph, whose false positives cluster in a few blocks; 5-Mb blocks widen them by a median factor of 1.16. Paired intervals are much narrower. Every HPRC − GRCh38-only SNV difference excludes zero (twelve pairs; R3.1), and so does every cross-tumor SNV loss. The two HPRC graphs differ by only 0.015–0.021 on HG008T, with intervals that just exclude zero, yet their three-seed means there differ by −0.011 to +0.015: the intervals cover the sampling of chr1, not training variance. Graph-coordinate calls lie outside rtg and these intervals (R3.4). Intervals for the differences to each comparator follow once they are re-scored [待补].

3. *Operating point.* The fixed 0.5 threshold is no longer used. PASS is the recall-0.9 validation threshold (90% somatic recall on the validation blocks), never tuned on the test chromosome. For SNVs it lands within 0.002–0.043 F1 of the test-tuned best F1 (upper bound) in all twelve in-sample HPRC models; across tumors it transfers less well (the COLO829T Illumina model reaches recall 0.157 on HCC1395 at its own threshold, while its curve peaks at F1 0.810; R1.7). For INDELs on long reads it falls far short of the test-tuned best F1; the validation best-F1 threshold instead raises long-read HG008T INDEL F1 by 0.102–0.139 (paired intervals exclude zero) but lowers it on Illumina and leaves only 0–41 calls on COLO829T (R3.10). We therefore show full curves for every model, report test-tuned best F1 only as a labelled upper bound, and give metrics at validation recall 0.95 and at argmax [待补].

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
- *Node selection:* alignments with MAPQ > 5; a node is a target when more than 5% of its alignments carry an edit.
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
