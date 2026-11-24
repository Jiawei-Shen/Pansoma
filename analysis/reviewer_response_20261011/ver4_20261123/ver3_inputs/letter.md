**<span class="underline">General response</span>**

We are extremely grateful for the reviews’ thoughtful comments and suggestions. This motivates us to improve Pansoma for better performance and somatic variant calling. Through additional work and a more thorough presentation of the Pansoma software, we have addressed all the issues raised. Below, we have provided our point-by-point responses.

**<span class="underline">Point-by-point responses</span>**

**Reviewer \#1 (Remarks to the Author):**  
  
This manuscript presents Pansoma, a machine-learning framework for tumor-only somatic variant calling directly from pangenome graph alignments. Pansoma aligns sequencing reads to a human pangenome graph, summarizes the resulting alignments into node-level pileups, and uses image-like tensors derived from these pileups for somatic variant prediction. The idea of incorporating population-scale pangenome information into tumor-only somatic variant calling is potentially interesting because pangenome references may help reduce reference bias and provide additional information for distinguishing germline from somatic variation, with the pangenome effectively serving as a population-level “generic normal” prior.  
  
However, in its current form, I do not think the method demonstrates what it claims. I have serious and fundamental concerns about the current study design and how the pangenome graph is utilized. Pansoma claims to be a pangenome-native somatic variant caller but does not exploit graph-specific information beyond alignment stage. Graph topology, edges, bubbles, cross-node connectivity, and broader haplotype/path information, is all discarded during training and variant inference. As a result, the downstream model resembles a conventional pileup-based deep-learning architecture, such as DeepVariant or DeepSomatic (Pansoma is not very much different from what DeepVariant and DeepSomatic have done), rather than a method specifically designed to model the pangenome itself. The benchmark results alone do not demonstrate that Pansoma meaningfully exploits graph-native information during variant inference. The manuscript therefore overstates its main methodological novelty.  
  
Another critical concern is that Pansoma does not convincingly demonstrate the specific advantages of using a pangenome graph. A major motivation for pangenome references is that population haplotypes and alternative paths can improve read alignment in repetitive, low-mappability, and structurally complex regions where a linear reference is limited. However, the manuscript provides no region-stratified analysis showing that Pansoma actually improves somatic variant calling in these regions, nor does it provide an ablation that separates the contribution of pangenome-based alignment from that of the downstream node-pileup model.  
  
Also, several aspects of the methodological design and evaluation, including the questionable choice of the HapMap mixture as a somatic benchmark and the relatively low recall. More seriously, the candidate-generation strategy may completely miss true somatic variants already represented on alternative graph paths, preventing them from ever reaching model inference. This exposes a fundamental conceptual problem in the current design. All these require substantial methodological reconsideration and raise concerns about whether the central claims of the manuscript are sufficiently supported.

**Response:** We thank the reviewer for the detailed and constructive comments. In response, we have substantially revised Pansoma's tensor representation to incorporate pangenome-native information and enable tensors to span adjacent graph nodes, allowing the model to capture graph topology and haplotype context. We have also expanded our evaluation to assess Pansoma's performance across different genomic regions and to better understand the sources of its improvements.

Specifically, we evaluated the contributions of the pangenome beyond graph alignment alone (R1.1–R1.3), quantified the somatic alleles absorbed by alternative graph paths and compared this effect with standard population filtering (R1.4), and moved the HapMap mixture experiments from the main benchmarking results to the Supplementary Information (R1.5). We now report sensitivity separately from F1-score (R1.6), reduced model complexity and expanded the training dataset with four additional tumor cell lines, including evaluations on tumors excluded from training (R1.7), and provided a comprehensive description of the training dataset and labeling procedure (R1.8).

Detailed point-by-point responses are provided below. 

Major concerns:

1\. The central claim that Pansoma is “pangenome-native” is not convincingly supported by its method design. The authors state that Pansoma can “fully exploit the richness of pangenome data” (line 89), yet in practice the pangenome is used primarily during vg Giraffe alignment, while the downstream representation remains essentially a local pileup. The released implementation does not explicitly model graph edges, bubbles/superbubbles, cross-node context, or long-range haplotype structure in either the input representation or its neural network. As a result, the method largely amounts to pangenome alignment followed by a standard pileup-based neural network, similar to existing DL variant callers such as DeepVariant, DeepSomatic, and ClairS. Replacing linear genomic coordinates with graph node IDs alone is insufficient to establish a method as pangenome-native. Moreover, haplotype- and path-level information, which represents one of the key advantages of a pangenome graph, is largely discarded. The authors should therefore provide substantial additional evidence or re-design their architecture to demonstrate that Pansoma derives measurable benefit from pangenome-specific information beyond the alignment stage itself.

**Response:** We agree that "fully exploit the richness of pangenome data" overstated the original design and that node identifiers alone do not make a method pangenome-native. We respectfully disagree, however, that Pansoma therefore resembles DeepSomatic or ClairS. Those callers take reads aligned to GRCh38 and build pileups in linear coordinates. They cannot produce Pansoma's tensors, even for a candidate confined to a single node, due to the great difference between graph and linear read alignments.

To make Pansoma pangenome-native, we revised the tensor format so that it encodes graph topology and haplotype context. Each row now follows its read's walk through the graph across node boundaries, and a dedicated channel records, at every position, how many haplotype paths pass through the node the read traverses.

1.  **Path-following context.** We revised our pipeline from single-node pileups to cross-node pileups. Each tensor row follows a read across nodes, preserving the flanking sequence within its 101-column window (Response Figure 1.1).

> We evaluated this change by masking this context at inference, which reduced SNV recall in all twelve comparisons (Response Figure 1.2).

2.  **Haplotype-count channel.** We added a new reference channel to encode the fraction of graph haplotypes that traverse each node on a read’s path (Response Figure 1.1). To evaluate its contribution, we retrained the SNV models with and without this channel using the same dataset and pipeline. Removing channel 6 reduced SNV F1 before population filtering in all twelve comparisons and after filtering in eleven of twelve; the only exception was HG008T ONT on v2.1. (Response Figure 1.3)

Response Figure 1.1

![Image preview](./media/image1.png)

**Response Figure 1.1.** The eight channels represent (1) read bases, (2) graph-reference bases corresponding to each read alignment, (3) the number of distinct GBWT haplotype paths supporting each traversed node, (4) base qualities, (5) the site allele supported by each read, (6) mapping qualities, (7) alignment operations and (8) read strand relative to the forward orientation of the candidate node. Black vertical lines delimit the candidate interval. **a.** An INDEL tensor whose rows follow their reads across node boundaries. An 8-bp deletion (TITATATA\>T) at chr1:23,017,177 (GRCh38) is shown using a candidate-centered tensor constructed from read alignments to the HPRC v2.1 graph. Reads are grouped by their allele support (A1, REF or OTHER), with 16 of 48 reads supporting the deletion (allele fraction = 0.333). **b.** A somatic G\>A substitution at chr1:242,529,453 (GRCh38) is shown using a candidate-centered tensor constructed from read alignments to the HPRC v1.1 graph.

Notably, channel 2 encodes the haplotype bases across the reads aligned reference base pairs, and channel 3 encodes node structure through path counts. Because every position within a node is traversed by the same haplotype paths, the path count is constant within a node, and a change in path count marks a boundary between nodes. Also, a bubble structure is visible in channel 3 in panel b.

Response Figure 1.2

![](./media/image2.png)

**Response Figure 1.2** **Evaluation of** **path-following context by** **masking graph context outside the candidate node at inference.** This evaluation is done by single tumor model respectfully (COLO829T and HG008T model) on HPRC v2.1 graph chromosome 1 SNVs in benchmark regions after population filtering using four resource VCFs. (a) Precision, recall and F1-scores. Blue circles represent the full input tensor, and orange diamonds represent tensors with all eight channels zeroed outside each read's traversal of the candidate node. (b) Differences in precision (P), recall (R) and F1-score between the two conditions, calculated as full tensor minus candidate-node-only input. Positive values indicate improved performance with the full tensor. The main improvement is on the recall.

Response Figure 1.3

![](./media/image3.png)

Response Figure 1.3 Effect of the haplotype-count **channel**. The evaluation was performed separately using single tumor models trained on COLO829T and HG008T, with chromosome 1 SNVs restricted to benchmark regions. For each tumor dataset, we trained two models differing only in the inclusion of the haplotype-count channel with identical configurations. (a) Precision, recall and F1-score after population filtering using four resource VCFs, together with F1-score before population filtering. Blue circles indicate models trained with the haplotype-count channel, and orange diamonds indicate models trained without it. Paired evaluations used the same Pansoma pipeline. (b) Differences in each metric shown in a, calculated as the value obtained with the haplotype-count channel minus that obtained without it. 

1.  
2.  
3.  
2\. Follow the first comment, a fundamental issue is that the authors do not provide a systematic analysis of how the composition and size of the pangenome reference affect somatic variant calling. It remains unclear how the number and diversity of reference haplotypes, haplotype walks, and alternative paths influence node-pileup construction and germline identification. This is crucial in understanding pangenome, because the graph itself determines which alleles are already represented as population variation and therefore which variants may be suppressed before reaching the classifier. This becomes important as newer references such as HPRC2 incorporate many more haplotypes. While a larger graph may improve population representation, it may also increase mapping ambiguity or cause true somatic mutations that overlap population alleles to be absorbed into existing graph paths and missed.

**Response:** We thank the reviewer for these comments. We acknowledge the reviewer's concern that increasingly complex pangenome graphs may absorb more somatic variants into alternative paths, potentially increasing false negatives. However, our quantitative analysis suggests that this effect is limited and does not represent an inherent limitation of the approach. Specially,

1.  **Net benefit in the three-reference ablation.** We have directly evaluated reference composition through a three-reference analyses of germline suppression and Pansoma performance. Even HPRC v2.1 contains more than five times as many haplotypes as HPRC v1.1, because of that nodes supported by fewer than 10% of haplotype paths are excluded by allele-frequency filtering, the node number increasement is limited.

The HPRC graphs reduced matched-normal germline candidates by 84–91% while largely retaining somatic SNV candidates (Response Figure 1.4). And we can see the performance of Pansoma on the three graphs; the benefit is especially evident before population filtering (Response Tables 1.1 and 1.2).

2.  **Somatic variants** **absorption.** We agree that pangenome graphs can absorb some somatic variants, but do not consider this a fundamental or critical issue in pangenome-based somatic calling. In the tumor-only setting studied, most of the existing methods adopt population filtering, and the existing population resource VCF files also overlaps the somatic tumors (Response Table 1.3). This represents a sensitivity–precision trade-off that must be assessed against the reduction in germline false positives and the resulting calling performance.

3.  **Mapping ambiguity**. Pansoma focuses on tumor-only somatic variant calling from graph alignments generated by upstream mapper vg Giraffe. Read placement and mapping ambiguity are primarily addressed by the upstream mapper, whose accuracy has been benchmarked and whose capabilities continue to develop (\[Sirén et al., 2021\] \[Chang et al., 2025, preprint\]). Furthermore, as shown in points 1 and 2, the number of somatic variants absorbed by alternative graph paths did not increase substantially, despite the increase in the number of haplotypes from 90 to 464.

Response Figure 1.4

![](./media/image4.png)

**Response Figure 1.4** We evaluated matched-normal germline alleles (5,061,892 for HG008T, 5,062,184 for COLO829T and 4,828,024 for HCC1395T) and somatic truth SNVs (8,690, 42,035 and 39,447, respectively) across chromosomes 1–22. A truth allele was considered retained if it was represented by a candidate with a matching representative allele. All truth records were included, irrespective of whether they fell within benchmark regions. The same sequencing reads, candidate-generation criteria and tensor construction code were used across all three graph references.

Response Table 1.1

| **chr1 SNVs, no population filter (P / R / F1)** | **Pansoma GRCh38-only** | **Pansoma HPRC v1.1** | **Pansoma HPRC v2.1**  | **ClairS-TO**         | **DeepSomatic-TO**    |
| ------------------------------------------------ | ----------------------- | --------------------- | ---------------------- | --------------------- | --------------------- |
| HG008T Illumina                                  | 0.005 / 0.838 / 0.010   | 0.106 / 0.848 / 0.189 | 0.194 / 0.864 / 0.316  | 0.003 / 0.902 / 0.006 | 0.029 / 0.935 / 0.056 |
| HG008T PacBio HiFi                               | 0.005 / 0.862 / 0.010   | 0.098 / 0.856 / 0.175 | 0.207 / 0.852 / 0.332  | 0.003 / 0.894 / 0.006 | 0.059 / 0.829 / 0.109 |
| HG008T ONT                                       | 0.004 / 0.867 / 0.009   | 0.094 / 0.855 / 0.170 | 0.190 / 0.878 / 0.312  | 0.003 / 0.884 / 0.006 | 0.047 / 0.910 / 0.090 |
| COLO829T Illumina                                | 0.044 / 0.593 / 0.081   | 0.297 / 0.722 / 0.421 | 0.580 / 0.728 / 0.644  | 0.010 / 0.796 / 0.020 | 0.097 / 0.858 / 0.173 |
| COLO829T Fiber-seq                               | 0.052 / 0.616 / 0.096   | 0.435 / 0.733 / 0.545 | 0.613 / 0.734 / 0.668  | 0.011 / 0.771 / 0.021 | 0.255 / 0.569 / 0.352 |
| COLO829T ONT                                     | 0.038 / 0.565 / 0.071   | 0.416 / 0.708 / 0.524 | 0.532 / 0.722 / 0.612  | 0.010 / 0.709 / 0.020 | 0.152 / 0.714 / 0.251 |
| HCC1395T Illumina                                | –                       | –                     | 0.706 / 0.830 / 0.763ᵃ | 0.014 / 0.884 / 0.028 | 0.139 / 0.944 / 0.242 |
| HCC1395T PacBio HiFi                             | –                       | –                     | 0.639 / 0.817 / 0.717ᵃ | 0.014 / 0.845 / 0.028 | 0.257 / 0.944 / 0.404 |
| HCC1395T ONT                                     | –                       | –                     | 0.560 / 0.704 / 0.624ᵃ | 0.013 / 0.752 / 0.025 | 0.217 / 0.781 / 0.340 |

**Response Table 1.1.** SNVs without any population filter on the three references, with the comparators (chr1, benchmark regions, PASS calls; precision / recall / F1). Pansoma before our population filter (it has none of its own): HG008T and COLO829T, the single-tumor models of each reference (trained on that tumor's chr2–22), mean of three seeds (GRCh38-only: training-set statistics; HPRC v1.1 and v2.1: base quality standardised per dataset). ᵃ HCC1395T has tensors only on HPRC v2.1 and no single-tumor model: the six-tumor model (R1.7; HCC1395T chr1 not in training, one run). ClairS-TO rerun without its panel-of-normals step; DeepSomatic-TO with its PON filter lifted (it still takes population allele frequencies as input).

Response Table 1.2

| **chr1 INDELs, no population filter (P / R / F1)** | **Pansoma GRCh38-only** | **Pansoma HPRC v1.1** | **Pansoma HPRC v2.1**  | **ClairS-TO**         | **DeepSomatic-TO**    |
| -------------------------------------------------- | ----------------------- | --------------------- | ---------------------- | --------------------- | --------------------- |
| HG008T Illumina                                    | 0.262 / 0.245 / 0.242   | 0.083 / 0.297 / 0.129 | 0.092 / 0.333 / 0.144  | 0.004 / 0.101 / 0.007 | 0.017 / 0.206 / 0.032 |
| HG008T PacBio HiFi                                 | 0.213 / 0.258 / 0.213   | 0.350 / 0.229 / 0.274 | 0.395 / 0.215 / 0.278  | 0.004 / 0.092 / 0.007 | 0.080 / 0.037 / 0.051 |
| HG008T ONT                                         | 0.019 / 0.361 / 0.035   | 0.025 / 0.233 / 0.044 | 0.011 / 0.300 / 0.020  | 0.004 / 0.087 / 0.007 | 0.043 / 0.019 / 0.026 |
| COLO829T Illumina                                  | 0.175 / 0.209 / 0.182   | 0.039 / 0.293 / 0.066 | 0.023 / 0.354 / 0.043  | 0.002 / 0.191 / 0.003 | 0.010 / 0.244 / 0.019 |
| COLO829T Fiber-seq                                 | 0.158 / 0.209 / 0.180   | 0.064 / 0.295 / 0.093 | 0.074 / 0.237 / 0.111  | 0.002 / 0.160 / 0.004 | 0.124 / 0.092 / 0.105 |
| COLO829T ONT                                       | 0.005 / 0.366 / 0.010   | 0.004 / 0.257 / 0.008 | 0.002 / 0.308 / 0.004  | 0.002 / 0.137 / 0.004 | 0.016 / 0.015 / 0.015 |
| HCC1395T Illumina                                  | –                       | –                     | 0.417 / 0.744 / 0.498ᵃ | 0.006 / 0.699 / 0.012 | 0.022 / 0.835 / 0.042 |
| HCC1395T PacBio HiFi                               | –                       | –                     | 0.515 / 0.699 / 0.589ᵃ | 0.005 / 0.526 / 0.011 | 0.271 / 0.647 / 0.382 |
| HCC1395T ONT                                       | –                       | –                     | 0.165 / 0.642 / 0.260ᵃ | 0.006 / 0.564 / 0.012 | 0.249 / 0.481 / 0.328 |

**Response Table 1.2.** INDELs without any population filter on the three references, with the comparators, as in Response Table 1.1. Pansoma: HG008T and COLO829T, single-tumor models with the final INDEL recipe on each reference, mean of three seeds; these models take gnomAD and CoLoRSdb allele frequencies as input, so, unlike the SNV models, their values without the filter are not free of population information. ᵃ HCC1395T has tensors only on HPRC v2.1: the final six-tumor INDEL models (HCC1395T chr1 and chr20 not in training), three seeds.

Response Table 1.3 Somatic candidate loss on the HPRC graph

<table>
<thead>
<tr class="header">
<th><strong>% of all truth alleles (ranges over three platforms)</strong></th>
<th><strong>HG008T<br />
INDELs</strong></th>
<th><strong>COLO829T<br />
INDELs</strong></th>
<th><strong>HCC1395T<br />
INDELs</strong></th>
<th><strong>HG008T<br />
SNVs</strong></th>
<th><strong>COLO829T<br />
SNVs</strong></th>
<th><strong>HCC1395T<br />
SNVs</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>truth alleles, chr1–22</td>
<td>8,496</td>
<td>1,912</td>
<td>1,625</td>
<td>8,690</td>
<td>42,035</td>
<td>39,447</td>
</tr>
<tr class="even">
<td>somatic truth alleles tagged by population filtering</td>
<td>62.7%</td>
<td>52.6%</td>
<td>7.9%</td>
<td>10.5%</td>
<td>3.4%</td>
<td>2.8%</td>
</tr>
<tr class="odd">
<td>candidate on the GRCh38-only graph but not on HPRC v1.1</td>
<td>31.1–44.8%</td>
<td>25.8–35.1%</td>
<td>4.1–5.5%</td>
<td>3.4–3.9%</td>
<td>0.5–0.9%</td>
<td>0.6–0.8%</td>
</tr>
<tr class="even">
<td><strong>graph-specific loss left after filtering, HPRC v1.1</strong></td>
<td>3.8–7.2%</td>
<td>5.0–6.8%</td>
<td>1.3–1.9%</td>
<td>0.3–0.4%</td>
<td>0.2–0.6%</td>
<td>0.1–0.2%</td>
</tr>
<tr class="odd">
<td>candidate on the GRCh38-only graph but not on HPRC v2.1</td>
<td>29.4–42.1%</td>
<td>23.9–32.0%</td>
<td>4.0–4.9%</td>
<td>3.0–3.6%</td>
<td>0.5–0.9%</td>
<td>0.8–0.9%</td>
</tr>
<tr class="even">
<td><strong>graph-specific loss left after filtering, HPRC v2.1</strong></td>
<td>3.5–7.0%</td>
<td>5.4–6.9%</td>
<td>1.1–1.4%</td>
<td>0.2–0.4%</td>
<td>0.2–0.7%</td>
<td>0.1–0.3%</td>
</tr>
</tbody>
</table>

**Response Table 1.3.** Somatic SNVs and INDELs lost during graph-based candidate generation and their overlap with population filtering. Population filtering follows our resource VCFs: four for SNVs; gnomAD/CoLoRSdb only for INDELs: SNVs use gnomAD, dbSNP 138 (non-somatic records), the GATK 1000 Genomes PoN and CoLoRSdb; INDELs use gnomAD and CoLoRSdb only. A filter tag indicates that a somatic truth allele meets the filtering criteria, not that it is germline in the patient. Candidate presence requires an exact match to any candidate allele, including a non-representative allele. Percentages use the complete chr1–22 truth counts. Ranges span three platforms.

1.  
2.  
3.  
4.  
3\. The manuscript fails to demonstrate one of the primary practical motivations for using a pangenome reference. Pangenome alignment can remain highly ambiguous in repetitive, paralogous, structurally complex, and low-mappability regions, especially for noisy ONT reads. This concern is also related to graph composition, since adding more haplotypes and alternative paths may introduce additional competing mappings in difficult regions. The authors report that graph alignment increases perfectly aligned reads and reduces mismatches, indels, and soft clipping relative to GRCh38. However, these aggregate statistics do not demonstrate that the improvement occurs in the regions where pangenome references are expected to provide the greatest advantage. 

This is important for ONT, for which the authors observe more graph-coordinate predictions and attribute this partly to higher sequencing error and more complex alignment profiles. Without region-specific benchmarking and error analysis, it remains unclear whether these graph-based signals represent improved resolution or simply additional alignment ambiguity. The current results therefore do not establish that Pansoma meaningfully addresses the challenging regions that motivate the use of a pangenome in the first place.

**Response:** We thank the reviewer for this comment and agree that the original manuscript lacked region-specific benchmarking. We have therefore added region-specific analyses of somatic variant calling performance, using the same sequencing reads and analysis pipeline across the GRCh38-only and HPRC graphs. These analyses cover easy, difficult and extreme genomic regions defined by the SMaHT consortium and provided as BED files (Response Table 1.4 SMaHT), as well as the GIAB v3.5 genomic stratifications. We have separately addressed the concern regarding alignment ambiguity in point 3 of our response to R1.2. The results are presented below:

**Calling in difficult regions.** Before the population filtering, SNV F1 was higher on the HPRC graphs in all twelve dataset–graph comparisons in the GIAB union of difficult regions. Gains extended to low-mappability, segmental-duplication and repeat regions, with no significantly lower SNV F1 in any difficult stratum. In SMaHT extreme regions, HPRC v2.1 reduced false-positive calls and improved F1 across all datasets. INDELs showed a different pattern: in \>=7bp homopolymers regions, Illumina gained, whereas long reads lost precision (Response Table 1.5).

**ONT graph-coordinate calls.** We admit that the ONT produced the most graph-coordinate PASS calls cannot be verified as the true tumor somatic variants. And we have narrow text in the manuscript for the calls on the off reference.

Response Table 1.4 SMaHT

| **Category** | **Definition**                                                                                                   | **Interpretation**                                                      |
| ------------ | ---------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------- |
| Easy         | Genomic regions outside the GIAB v3.3 difficult-region union                                                     | Regions expected to be relatively straightforward for variant detection |
| Difficult    | Subsets of GIAB-defined difficult regions, divided into four non-overlapping tiers according to Umap mappability | Challenging regions that retain some degree of short-read mappability   |
| Extreme      | Regions outside both Easy and Difficult                                                                          | The most challenging regions, characterized by poor mappability         |

Response Table 1.4 SMaHT. SMaHT consortium provides explicit, reproducible definitions of Easy, Difficult, and Extreme genomic regions.

Response Table 1.5 Regional variant calling before population filtering

<table>
<thead>
<tr class="header">
<th><strong>chr1 dataset</strong></th>
<th><strong>HPRC graph</strong></th>
<th><strong>SMaHT extreme</strong><br />
<strong>PASS calls (precision)</strong></th>
<th><strong>SMaHT extreme</strong><br />
<strong>ΔF1</strong></th>
<th><strong>GIAB all difficult</strong><br />
<strong>ΔF1</strong></th>
<th><strong>GIAB low map or segdup</strong><br />
<strong>ΔF1</strong></th>
<th><strong>GIAB tandem</strong><br />
<strong>repeats and</strong><br />
<strong>homopolymers</strong><br />
<strong>ΔF1</strong></th>
<th><strong>GIAB outside difficult</strong><br />
<strong>ΔF1</strong></th>
<th><strong>INDEL F1</strong><br />
<strong>homopolymers</strong><br />
<strong>≥ 7 bp</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>HG008T<br />
Illumina</td>
<td>v1.1</td>
<td>15,955 (0.003)<br />
→ 601 (0.045)</td>
<td>+0.072<br />
[0.046, 0.099]</td>
<td>+0.147<br />
[0.120, 0.177]</td>
<td>+0.143<br />
[0.097, 0.193]</td>
<td>+0.103<br />
[0.073, 0.137]</td>
<td>+0.188<br />
[0.167, 0.210]</td>
<td>0.027<br />
→ 0.121</td>
</tr>
<tr class="even">
<td>HG008T<br />
Illumina</td>
<td>v2.1</td>
<td>15,955 (0.003)<br />
→ 599 (0.065)</td>
<td>+0.107<br />
[0.073, 0.144]</td>
<td>+0.207<br />
[0.172, 0.247]</td>
<td>+0.214<br />
[0.156, 0.279]</td>
<td>+0.133<br />
[0.094, 0.178]</td>
<td>+0.315<br />
[0.286, 0.344]</td>
<td>0.027<br />
→ 0.138</td>
</tr>
<tr class="odd">
<td>HG008T<br />
PacBio HiFi</td>
<td>v1.1</td>
<td>17,989 (0.003)<br />
→ 489 (0.070)</td>
<td>+0.111<br />
[0.077, 0.143]</td>
<td>+0.167<br />
[0.137, 0.198]</td>
<td>+0.144<br />
[0.100, 0.187]</td>
<td>+0.142<br />
[0.099, 0.186]</td>
<td>+0.186<br />
[0.165, 0.208]</td>
<td>0.027<br />
→ 0.046</td>
</tr>
<tr class="even">
<td>HG008T<br />
PacBio HiFi</td>
<td>v2.1</td>
<td>17,989 (0.003)<br />
→ 330 (0.112)</td>
<td>+0.169<br />
[0.121, 0.218]</td>
<td>+0.263<br />
[0.221, 0.308]</td>
<td>+0.198<br />
[0.139, 0.261]</td>
<td>+0.221<br />
[0.166, 0.279]</td>
<td>+0.313<br />
[0.285, 0.341]</td>
<td>0.027<br />
→ 0.056</td>
</tr>
<tr class="odd">
<td>HG008T<br />
ONT</td>
<td>v1.1</td>
<td>19,676 (0.003)<br />
→ 571 (0.058)</td>
<td>+0.094<br />
[0.064, 0.126]</td>
<td>+0.153<br />
[0.126, 0.184]</td>
<td>+0.128<br />
[0.088, 0.170]</td>
<td>+0.130<br />
[0.093, 0.171]</td>
<td>+0.179<br />
[0.159, 0.201]</td>
<td>0.010<br />
→ 0.010</td>
</tr>
<tr class="even">
<td>HG008T<br />
ONT</td>
<td>v2.1</td>
<td>19,676 (0.003)<br />
→ 669 (0.070)</td>
<td>+0.118<br />
[0.083, 0.159]</td>
<td>+0.195<br />
[0.159, 0.236]</td>
<td>+0.160<br />
[0.109, 0.227]</td>
<td>+0.148<br />
[0.107, 0.192]</td>
<td>+0.313<br />
[0.284, 0.343]</td>
<td>0.010<br />
→ 0.009</td>
</tr>
<tr class="odd">
<td>COLO829T<br />
Illumina</td>
<td>v1.1</td>
<td>8,143 (0.016)<br />
→ 3,102 (0.045)</td>
<td>+0.053<br />
[0.023, 0.165]</td>
<td>+0.101<br />
[0.052, 0.243]</td>
<td>+0.044<br />
[0.019, 0.158]</td>
<td>+0.042<br />
[0.019, 0.147]</td>
<td>+0.535<br />
[0.510, 0.557]</td>
<td>0.007<br />
→ 0.032</td>
</tr>
<tr class="even">
<td>COLO829T<br />
Illumina</td>
<td>v2.1</td>
<td>8,143 (0.016)<br />
→ 296 (0.449)</td>
<td>+0.471<br />
[0.399, 0.541]</td>
<td>+0.508<br />
[0.460, 0.554]</td>
<td>+0.463<br />
[0.381, 0.546]</td>
<td>+0.473<br />
[0.393, 0.542]</td>
<td>+0.629<br />
[0.613, 0.644]</td>
<td>0.007<br />
→ 0.041</td>
</tr>
<tr class="odd">
<td>COLO829T<br />
Fiber-seq</td>
<td>v1.1</td>
<td>5,684 (0.023)<br />
→ 844 (0.172)</td>
<td>+0.224<br />
[0.123, 0.369]</td>
<td>+0.286<br />
[0.193, 0.390]</td>
<td>+0.196<br />
[0.097, 0.362]</td>
<td>+0.323<br />
[0.268, 0.377]</td>
<td>+0.501<br />
[0.472, 0.526]</td>
<td>0.005<br />
→ 0.010</td>
</tr>
<tr class="even">
<td>COLO829T<br />
Fiber-seq</td>
<td>v2.1</td>
<td>5,684 (0.023)<br />
→ 321 (0.455)</td>
<td>+0.481<br />
[0.393, 0.570]</td>
<td>+0.503<br />
[0.447, 0.559]</td>
<td>+0.445<br />
[0.343, 0.561]</td>
<td>+0.415<br />
[0.323, 0.510]</td>
<td>+0.631<br />
[0.615, 0.646]</td>
<td>0.005<br />
→ 0.008</td>
</tr>
<tr class="odd">
<td>COLO829T<br />
ONT</td>
<td>v1.1</td>
<td>10,110 (0.012)<br />
→ 815 (0.171)</td>
<td>+0.241<br />
[0.168, 0.338]</td>
<td>+0.293<br />
[0.221, 0.375]</td>
<td>+0.217<br />
[0.134, 0.342]</td>
<td>+0.276<br />
[0.208, 0.349]</td>
<td>+0.478<br />
[0.454, 0.501]</td>
<td>0.002<br />
→ 0.002</td>
</tr>
<tr class="even">
<td>COLO829T<br />
ONT</td>
<td>v2.1</td>
<td>10,110 (0.012)<br />
→ 492 (0.274)</td>
<td>+0.348<br />
[0.254, 0.482]</td>
<td>+0.407<br />
[0.321, 0.500]</td>
<td>+0.303<br />
[0.207, 0.458]</td>
<td>+0.452<br />
[0.391, 0.511]</td>
<td>+0.595<br />
[0.580, 0.610]</td>
<td>0.002<br />
→ 0.001</td>
</tr>
</tbody>
</table>

**Response Table 1.5.** SNV calling by region on the HPRC graphs against the GRCh38-only graph (chr1, benchmark regions, before population filtering (no filtered by the four resource VCFs)). SMaHT extreme regions: PASS calls and precision on each graph and the paired ΔF1 (HPRC − GRCh38-only) with its 95% interval (2,000 resamples of 1-Mb blocks). GIAB v3.5 stratifications: all difficult regions (46.8 Mb of chr1), low mappability or segmental duplications, tandem repeats and homopolymers, and the complement of all difficult regions, with the same paired intervals. Last column: INDEL F1 in homopolymers of at least 7 bp, from the earlier INDEL models of the three-reference comparison. HP denotes homopolymers. In the ΔF1 columns, the second line gives the 95% CI. Arrows show GRCh38-only → HPRC.

4.Candidate generation in Pansoma can delete true somatic variants that are encoded in the graph. Pansoma defines candidates only through mismatches and indels against the traversed node. If a true somatic allele is already encoded in an alternative path and tumor reads align perfectly to it, Pansoma treats it as the local graph reference and eliminates it before inference: no candidate, tensor, or call is produced. This is a fatal methodological issue. Supplementary Figure 3 exposes this failure, yet the authors incorrectly describe these false negatives as GRCh38 reference-bias artifacts. The presence of an allele on an alternative path in the pangenome does not mean that it is a germline variant. The same allele can be germline in one individual and somatic in another.

**Response:** We agree that this mechanism exists: an allele on a graph path is not necessarily germline in the patient, and tumor reads that align to it without a mismatch produce no candidate. We have withdrawn the reference-bias interpretation of Supplementary Fig. 3 and matched-normal reads and assemblies show that most of these alleles are genuine somatic events. However, We do not consider this a fatal methodological issue. As discussed in R1.2, graph-mediated suppression can be viewed as analogous to population filtering applied during alignment. And the quantify analysis shows that Pangenome still contribute positively to the somatic variant calling. In our benchmarks, both DeepSomatic-TO and ClairS-TO used population-based filtering, and these resources also contain alleles present in the somatic ground-truth sets. To assess how many somatic variants are lost specifically through graph absorption, we quantified the graph-specific loss across the complete truth sets and its overlap with population filtering (Response Table 1.3).

**Tumor-only population filtering context.** As discussed in R1.2, population filtering involves a trade-off: it is widely used to remove false positive when a matched normal is unavailable but can also exclude true somatic variants. Graph-mediated suppression uses overlapping population information at an earlier stage. We therefore quantified graph-specific absorption of somatic SNVs and INDELs and its overlap with population filtering (Response Table 1.3).

**Somatic loss and germline suppression***.* Graph-based candidate generation does suppress true somatic SNVs and INDELs. However, this loss substantially overlaps population-based filtering, the population resource VCF files are also used by DeepSomatic-TO and ClairS-TO. After applying the same population filter across references, the additional graph-specific loss is 0.2–0.7% of truth SNVs and 3.5–7.0% of truth INDELs (Response Table 1.3). At the same candidate-generation stage, the HPRC graphs v2.1 reduced matched-normal germline candidates by 84–91% relative to GRCh38-only (Response Figure 1.4).

5\. The HapMap low-VAF experiment raises a fundamental concern about whether the proposed pangenome-based mechanism is being evaluated with an appropriate somatic benchmark. The manuscript explicitly states that the “somatic-like” variants in the HapMap mixture are actually germline variants from the minor HapMap individuals, which are used to low VAFs by mixing multiple cell lines. However, a central rationale of Pansoma is that germline alleles already represented in the pangenome can align cleanly to existing graph paths and therefore be suppressed as somatic candidates. The authors themselves report that a large fraction of the germline variants in the benchmark normal samples are already present in the HPRC graph, and they argue that this representation helps suppress germline‑like signals. These two assumptions are potentially contradictory: the positive variants in the HapMap benchmark are themselves population germline variants, whereas true acquired somatic mutations would generally not be expected to follow the same population representation. This raises a major concern about the biological validity of the training and evaluation strategy using pangenome graph using synthetic HapMap results. As a result, the HapMap mixture may systematically disadvantage a pangenome-based caller and does not provide a biologically appropriate setting for evaluating its ability to detect true acquired somatic mutations at low VAF.  
  
The HapMap experiment can also suffer from data leakage. HPRC v1.1 contains four of the six HapMap donors used in the mixture HG00438, HG02257, HG02486, and HG02622. The authors themselves also raise this point in a companion paper (Fu, Xin, Miao et al., "Leveraging Human Pangenome for Improved Somatic Variant Detection," bioRxiv 2026.01.04.697580). Yet these donors are not excluded when benchmarking Pansoma.

**Response:** We thank the reviewer for raising these concerns. As a synthetic mixture of germline variation from multiple individuals, HapMap is less suitable as a primary benchmark for Pansoma. We have therefore moved it to the Supplement as an ultra-deep dataset for low-allele-fraction detection (about 500×; truth VAF 1–16.5%). Somatic calling at reduced purity is assessed using COLO829T tumor/normal dilutions. The original manuscript may not have clearly described our design, but the donor-related leakage suggested by the reviewer did not occur.

**Data leakage.** We acknowledge that the original manuscript did not explain this design clearly enough; however, our experiment did not involve the donor-related data leakage raised by the reviewer. The HapMap experiment used a separately constructed graph rather than either standard HPRC graph, built from 40 HPRC version 1 individuals plus GRCh38 and CHM13 (82 haplotypes), using the same Minigraph-Cactus construction and d9 filter. We excluded all six mixture individuals: HG005, the background individual, and the five minor donors HG002, HG00438, HG02257, HG02486 and HG02622 (Fu, Xin, Miao et al., 2026, Cell Genomics). Our low-VAF model used this graph for both training (chr2–22) and evaluation (held-out chr1 of the mixture); neither stage used a graph containing the mixture individuals. The Methods section now clarify this design.

Positive alleles can still be graphing paths if another individual carries them (44.8% of chr1 SNV and 20.8% of INDEL positives). Such reads can align without a mismatch, rarely becoming candidates (1.3% of these SNVs reach a tensor). A shared evaluation mask excludes graph-carried alleles, and their flanking reference spans for both Pansoma and Mutect2 (Response Table 1.6), so excluded alleles do not contribute directly to precision or recall. Mutect2 has higher SNV recall, whereas Pansoma has higher SNV precision and INDEL recall.

Response Table 1.6

|                                                   | SNVs            | INDELs         |
| ------------------------------------------------- | --------------- | -------------- |
| truth alleles                                     | 392,567         | 124,499        |
| carried by the graph (left out of the evaluation) | 175,773 (44.8%) | 25,938 (20.8%) |
| scored                                            | 215,094         | 92,079         |

**Response Table 1.6.** HapMap mixture without the graph-carried truth alleles (chr1, benchmark regions, truth VAF 1–16.5%).

6\. Although Pansoma achieves higher F1 scores in several benchmark settings, these gains appear to be driven primarily by improved precision rather than improved sensitivity. Pansoma does not demonstrate a convincing sensitivity advantage for tumor-only somatic variant calling, which limits its practical value for tumor-only downstream interpretation. Supplementary Figure 1 and Tables 2–3 show that Pansoma generally has lower SNV recall than ClairS-TO and DeepSomatic-TO across COLO829T and HG008T. The situation is more concerning for INDELs, where Pansoma recall is only \~4–22% in COLO829T and \~13–19% in HG008T. In my understanding, sensitivity is of great importance because there is no matched normal sample or second calling stage that can readily rescue somatic variants missed during candidate generation or classification. Missing such a large fraction of true variants is a serious limitation for tumor-only analysis, where false negatives may be permanently excluded from downstream interpretation in many discovery and clinical-prioritization settings. A higher F1 score alone does not compensate for this loss of sensitivity in many practical tumor-only applications.

**Response:** We agree that sensitivity is essential in tumor-only calling and that F1 alone cannot fully capture it. We revised the overall pipeline and increased the loss weight assigned to true somatic examples during training, improving SNV recall. INDEL sensitivity and sensitivity below VAF 0.08 remain the main limitations of Pansoma, but it is competitive compared to the compared tools.

We explored recall analysis (Response Table 1.7, Response Figure 1.6), with full precision–recall curves for every caller (Response Figure 1.5). For each caller, we compared them with the truth set inside the benchmark regions using RTG vcfeval and traced the precision–recall curve by sweeping a threshold on the QUAL value of the caller's VCF ​respectively. The circle on each curve marks that caller's PASS call set after the same filter. The precision–recall curves show that Pansoma does not achieve a high F1 score by trading precision for recall. 

We further examined the truth variants that receive no candidate (Response Figure 1.6). The extent to which somatic variants are absorbed by the graph varies considerably across datasets. For INDELs, we acknowledge that the graph absorbs a substantial fraction of candidates in COLO829T and HG008T (Response Figure 1.6b). However, as displayed in the figures and tables below, this does not critically limit performance, and Pansoma shows no significant loss of recall (Response Figure 1.5). Additionally, the majority of the absorbed INDEL are homopolymers and tandem repeats (Response Table 1.7).

Response Figure 1.5

![](./media/image5.png)

**Response Figure 1.5. Precision–recall curves on chr1 for Pansoma, DeepSomatic-TO and ClairS-TO.** (a) SNVs and (b) INDELs in HCC1395T, COLO829T and HG008T (rows) on three platforms (columns), in the benchmark regions after the same population filter for every caller. Curves rank calls by each caller’s score (Pansoma: somatic probability; DeepSomatic-TO and ClairS-TO: QUAL, all records except RefCall and NoCall). Dots mark each caller’s PASS calls after the filter. Pansoma uses the six-tumor HPRC v2.1 models; chr1 of these three tumors was not used in training. Recall axes in b are cut to the range of the curves.

Response Figure 1.6

![](./media/image6.png)

Response Figure 1.6 | Why truth variants receive no candidate. Truth variants on chr1–22 without a candidate on the HPRC v2.1 graph (all truth alleles, inside and outside the benchmark regions), classified from the reads at each site. One pie per dataset; n, number of truth variants without a candidate; slices of 8% or more are labelled. Colors are defined separately for a and b. a, SNVs. S3: reads follow the graph branch carrying the ALT, leaving no edit. S4, repeats: (a) the SNV is written as an edit on another branch node; (b) the ALT haplotype is a graph path; (c) reads show no clear ALT. b, INDELs. I1: ALT reads follow a graph branch without an edit. I2: ALT reads follow a similar graph allele, retaining a residual edit with a different spelling.

Response Table 1.7

| **Truth INDELs (chr1–22)**                     | **HG008T**        | **COLO829T**    | **HCC1395T**   |
| ---------------------------------------------- | ----------------- | --------------- | -------------- |
| Absorbed by the graph (union across platforms) | 2,277             | 430             | 52             |
| Homopolymer ≥ 7 bp                             | 1,853 (81.4%)     | 241 (56.0%)     | 30 (57.7%)     |
| Short tandem repeat (unit 2–6 bp, ≥ 3 copies)  | 292 (12.8%)       | 83 (19.3%)      | 5 (9.6%)       |
| **Homopolymer + short tandem repeat, total**   | **2,145 (94.2%)** | **324 (75.3%)** | **35 (67.3%)** |
| Other repeats                                  | 6 (0.3%)          | 20 (4.7%)       | 1 (1.9%)       |
| Substitution within repetitive sequence        | 104 (4.6%)        | 16 (3.7%)       | 2 (3.8%)       |
| Others                                         | 21 (0.9%)         | 61 (14.2%)      | 12 (23.1%)     |
| Non-repetitive sequence                        | 1 (0.0%)          | 9 (2.1%)        | 2 (3.8%)       |

**Response Table 1.7.** Sequence context of somatic truth INDELs absorbed by the HPRC v2.1 graph (chr1–22). An INDEL is absorbed when, on at least one platform, its ALT reads align to a graph path without any edit, leaving no candidate. First row: absorbed INDELs and their percentage of all truth INDELs; other rows: percentages of the absorbed INDELs. Classes describe the somatic change relative to the matched-normal haplotype. Other repeats: shorter homopolymers (4–6 bp), longer-unit tandem repeats, or INDELs that do not change a repeat by whole units. Substitution within repetitive sequence: a base substitution on a germline length difference. Others: allele already on the normal haplotype, complex change, or unresolved normal haplotype.

7\. The training data may be insufficient in scale for a model of this size, and insufficient in tumor diversity for a general-purpose tumor-only somatic variant caller. Pansoma uses a ConvNeXt backbone with around 200M parameters, but Pansoma seems to rely primarily on COLO829 for model training, which provides only on 40k-50k truth somatic variants. It is therefore unclear whether this is sufficient to train such a large neural network. The authors should report the exact parameter count and the number of positive and negative training examples per platform and variant type and provide evidence that the data volume supports the chosen capacity. Moreover, COLO829 is also a single melanoma line dominated by UV signatures, so the model may be learning tumor-specific rather than general somatic features. By comparison, DeepSomatic was trained across multiple cancer cell lines with extensive augmentation to build a broader training distribution. Authors should consider extending training to more types of tumors with a distinct mutational process and demonstrate directly that the current model is not tumor specific.

**Response:** We thank the reviewer for this comment. We reduced model size, expanded training to additional tumor types and tested the models on tumors excluded from training. Response Table 1.8 reports training counts.

**Model size and training data. **We understand the reviewer’s concern that a large model may memorize the training data and overfit. The checkpoints we use do not show this: each is selected on validation data held out from training. We tested the smaller model on the toy sets and reduced the SNV and INDEL networks from 199.6 M to 50.4 M and 8.7 M parameters, respectively, because the smaller networks achieve similar F1 in our tests. As the reviewer requested, Response Table 1.8 reports training counts by model and class. More generally, unlike traditional statistical models, modern deep neural networks are often highly over-parameterized and can still generalize well even when the number of parameters far exceeds the number of training samples.

**Tumor diversity,** **multi-tumor training and evaluation.** We added HCC1395 and HCC1954 (breast), and H1437 and H2009 (lung adenocarcinoma), each on three platforms, broadening the mutational spectra beyond COLO829T's UV signature. We trained four-tumor models using HCC1395 and HCC1954 (breast) and H1437 and H2009 (lung adenocarcinoma), and six-tumor models additionally including chr2–22 of COLO829T and HG008T. We evaluated these models across tumor datasets. Four-tumor models achieve HG008T SNV F1 comparable to its own models and above both comparators (DeepSomatic-TO and ClairS-TO) but remain below both on COLO829T Fiber-seq and ONT. Six-tumor models achieve higher SNV F1 than both comparators on HG008T and COLO829T; these are within-tumor tests with chr1 held out.

Response Table 1.8.

<table>
<thead>
<tr class="header">
<th><strong>Model</strong><br />
<strong>(graph)</strong></th>
<th><strong>Somatic</strong><br />
<strong>(train)</strong></th>
<th><strong>Somatic</strong><br />
<strong>(validation)</strong></th>
<th><strong>Germline</strong><br />
<strong>(train)</strong></th>
<th><strong>Other</strong><br />
<strong>(train)</strong></th>
<th><strong>Other used</strong><br />
<strong>per epoch</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Single-tumor SNV<br />
(HPRC v1.1 and v2.1)</td>
<td>7,098–35,196</td>
<td>350–2,127</td>
<td>256,583–381,159</td>
<td>67,735–2,078,603</td>
<td>67,735–1,620,829</td>
</tr>
<tr class="even">
<td>Single-tumor INDEL<br />
(HPRC v1.1 and v2.1)</td>
<td>971–3,761</td>
<td>52–194</td>
<td>55,260–220,610</td>
<td>20,551–3,289,717</td>
<td>20,551–3,289,717</td>
</tr>
<tr class="odd">
<td>COLO829T four-purity SNV<br />
(HPRC v1.1)</td>
<td>124,070–125,916</td>
<td>7,614–7,713</td>
<td>1,356,321–1,665,260</td>
<td>421,331–3,047,060</td>
<td>421,331–3,047,060</td>
</tr>
<tr class="even">
<td>Four-tumor SNV<br />
(HPRC v2.1)</td>
<td>237,577–245,732</td>
<td>11,592–12,005</td>
<td>1,389,181–1,590,228</td>
<td>314,119–21,105,224</td>
<td>219,951–7,759,510</td>
</tr>
<tr class="odd">
<td>Six-tumor SNV<br />
(HPRC v2.1)</td>
<td>278,615–288,202</td>
<td>13,664–14,144</td>
<td>1,939,846–2,248,674</td>
<td>589,993–24,803,310</td>
<td>430,892–9,682,334</td>
</tr>
<tr class="even">
<td>Four-tumor INDEL<br />
(HPRC v2.1)</td>
<td>9,717–9,930</td>
<td>431–462</td>
<td>324,181–827,146</td>
<td>408,782–10,465,888</td>
<td>408,782–10,465,888</td>
</tr>
<tr class="odd">
<td>Six-tumor INDEL<br />
(HPRC v2.1)</td>
<td>14,132–14,854</td>
<td>620–682</td>
<td>440,263–1,215,616</td>
<td>467,482–15,606,962</td>
<td>467,482–15,606,962</td>
</tr>
</tbody>
</table>

**Response Table 1.8.** Training examples per model (summary; the full table lists every model by class, tumor, platform and variant type). Four-tumor models: HCC1395, HCC1954, H1437 and H2009; six-tumor models: the same four lines, HG008T and COLO829T. SNV models: one run (seed 0) each; four-tumor and six-tumor INDEL models (the final INDEL models): three seeds each. COLO829T four-purity sets count each dataset separately, so a truth allele can appear up to four times (125,916 Illumina SNV positives are 33,797 distinct alleles). Validation “other” also contains off-reference tensors without a truth match, which training does not use. In the current SNV pipeline (the Illumina and PacBio HiFi / Fiber-seq models of COLO829T and the four-tumor and six-tumor models), candidates of the four added lines and of COLO829T are built with an allele-frequency floor of 0.05 (HG008T and all ONT sets: 0.08), and “other” examples with allele frequency below 0.1 are sampled at 25% in each epoch (last column); all other models use every example in every epoch. The four-tumor and six-tumor SNV models use all autosomes of HCC1954, H1437 and H2009 and chr2–19, 21, 22 of HCC1395; chr1 of HG008T, COLO829T and HCC1395 is never in training. The four-tumor and six-tumor INDEL models hold out chr1 of every tumor (and chr20 of HCC1395).

8\. The construction of the training set is not also described in enough detail to be evaluated, and the choices that were made are also not justified. For positive examples, truth variants defined on GRCh38 are projected onto graph coordinates and labeled by exact allele match. However, the same biological variant can be represented differently across graph paths, especially for indels, and graph-only positions are assigned GRCh38 coordinates indirectly through adjacent 1-bp SNV nodes. It is therefore unclear how a positive label is assigned when the truth allele is already represented as a graph path, when a projection is ambiguous, or when the graph and linear representations of the same variant differ. The authors' own false-negative analysis (Supplementary Fig. 3) shows that some truth variants have their apparent allele fraction reduced to nearly zero once reads align cleanly to an alternative path, so these sites never reach candidate generation and cannot become positives at all. The authors should describe the labeling and projection procedure explicitly, including whether normalization was applied on the graph side, and justify why this procedure does not systematically distort the positive set.  
  
For negative examples, the only definition given is that determinable non-matches were labeled negative, with no explanation of what makes a candidate determinable. This matters more for a tumor-only caller than for a paired one, because the model must learn to reject germline alleles and alignment artifacts without a matched normal, and it can only do so if these classes are actually represented among the negatives. The authors should explain how the negative set was assembled, including whether germline variants were explicitly included, how artifacts were defined, and how VAF and tumor-purity conditions were sampled.

**Response:** We thank the reviewer and agree that the original description of the training set was insufficient. The revised Methods and a Supplementary Note describe it in full. We made each choice explicit and conservative and anchored the labels to the GRCh38 path, on which the truth sets are defined. In the absence of an established convention for labelling pangenome graph alignments, some choices differ from those used with linear references. These differences reflect graph-based representation and do not, by themselves, indicate labelling errors.

1.  **Positions***.* Our truth VCFs use GRCh38 coordinates. Each HPRC graph also retains GRCh38 as a reference path, allowing truth loci to be converted into graph node-and-offset coordinates where the path placement is unique. These coordinates locate the truth somatic variants on the graph; ambiguous placements are not forced into exact labels.

2.  **Labelling***.* Candidates are discovered independently from graph-aligned read pileups and compared with the projected truth alleles (R1.minor2). The representative allele, supported by the most ALT reads, determines each site's label: an exact somatic truth match receives label 1. A truth match confined to another allele does not automatically make the site positive. If no read pileup supports an eligible candidate, the truth variant is missed. INDELs are left-aligned and compared at repeat-equivalent truth placements. In repeats, part of an INDEL may match ALT nodes, leaving a partial edit. The original labelling rule assigns label 1 to partial matches with more than 45% allele or read-haplotype overlap; these partial annotations are retained, but the current training recipe excludes them. The current pipeline also assigns germline label 2 using matched-normal truth VCFs and the same projection procedure. Unmatched candidates outside the confident region (the intersection of somatic and germline benchmark regions), and off-reference candidates without a truth match, are excluded with label −1; remaining eligible unmatched candidates receive label 0 (other). Unlabeled off-reference candidates count as validation negatives, a model-selection assumption rather than proof that they are non-somatic. They remain eligible for inference, and output can include calls in graph coordinates.

3.  **VAF and purity. **Candidate VAF is estimated from the supporting reads in the graph-aligned pileup (R1.m2). Training retains each tumor's observed VAF distribution; joint COLO829T models additionally use in-silico tumor/normal dilutions at 75%, 50% and 25% purity.

The situations described in the comment involve graph absorption or differences in variant representation:

  - **Different representations across paths.** This occurs particularly for INDELs. In tandem repeats, homopolymers and other short repeats, reads may partly match existing graph sequence, leaving a residual INDEL edit with a different spelling (Response Figure 1.6). Related SNV respelling is less common (Response Figure 1.6). Pansoma records partial matches, but a residual INDEL alone is not credited as the complete GRCh38 truth allele during evaluation.

  - **Graph-only positions.** Borrowing coordinates from neighboring 1-bp nodes was an original projection convention that has been removed. The current pipeline follows each read across the corresponding graph-reference bases and left-aligns candidate INDELs in node-forward coordinates. In most cases, Giraffe places a read base on a single-base ALT node only when the base matches that node exactly. A truth SNV carried by such a node can therefore be absorbed rather than emitted as a mismatch; these SNV cases are uncommon (Response Figure 1.6). Path changes can also split an INDEL into multiple residual edits, but such cases are extremely rare. We record them as partial INDEL matches and exclude them from current training, equivalent to label −1.

  - **Truth alleles absorbed by graph paths.** Reads matching an ALT path without an edit yield no candidate (Response Figure 1.6). This mechanism is discussed in the previous response.

These cases therefore fall into two categories: absorption of an allele already represented by the graph, and a graph-relative edit that differs from the complete linear-reference allele. The latter representation is inherited from upstream alignments to ALT paths and need not indicate a mapping error; partial somatic matches are nevertheless excluded from current training. Pansoma supports analysis directly in graph coordinates, providing a foundation for future analysis and benchmarking with graph references. These graph-relative candidate definitions, together with the fixed-alignment and fixed-candidate ablation in R1.1, contradict the reviewer’s statement that “Pansoma is not very much different from what DeepVariant and DeepSomatic have done”. Similarity in the neural-network architecture does not establish equivalence in the variant-calling method.

Minor concerns:  
1\. Beyond COLO829 and HG008, the authors should benchmark Pansoma on additional well-established datasets such as the SEQC2 HCC1395 cell lines, which provides a high-confidence somatic variant truth set, as well as the five additional cell lines evaluated in the DeepSomatic study.

We added four cancer cell lines: HCC1395 with the SEQC2 v1.2.1 truth set, and HCC1954, H1437 and H2009 from CASTLE, covering breast cancer and lung adenocarcinoma. All four have Illumina, PacBio HiFi and ONT data processed on HPRC v2.1. In the revised manuscript, we benchmarked both the four-tumor model (trained on HCC1395 and three CASTLE tumor datasets) and the six-tumor model (trained on the same four datasets plus HG008T and COLO829T).

2\. Pansoma’s tensor windows are restricted to the single current node and truncated at node boundaries. Adjacent-node and path-spanning evidence has already been discarded. This contradicts the authors’ claim that Pansoma can “fully exploit the richness of pangenome data.” This design also makes it difficult to deal with short nodes or even single-base node. In fact, a large fraction of alternative nodes has the length shorter than the length of Pansoma’s tensor window. A credible graph-native method requires path-aware tensors spanning neighboring nodes rather than arbitrarily terminating the context at node boundaries and explicitly model graph edges, bubbles, alternative-path relationships, cross-node context, or competing path support.

As noted in R1.1, we agree that “fully exploit the richness of pangenome data” overstated the original design. The revised pipeline removes the node-boundary truncation described in the original Methods. At each column, the tensor records the traversed node’s graph base and haplotype-path count (panels 2 and 3 of Response Figure 1.1). Reads are grouped by candidate allele and ordered by alignment events, including node transitions, to organize evidence from different graph walks. Two examples show how the graph-reference and path-count channels (panels 2 and 3 of Response Figure 1.1) carry graph structure into the tensor. In panel a, a somatic 8-bp deletion spans six graph nodes of 1–2 bp, carried by 19 to 78 of the 464 haplotypes; the path count changes from node to node across the deletion. In panel b, 21 bases from a somatic SNV, all 38 ALT-supporting reads traverse a bubble branch carried by 18 of 90 haplotypes, while all 40 REF-supporting reads traverse the GRCh38 branch carried by 72. The original representation would have truncated this context. The ablations of both features are addressed in R1.1.

Reviewer \#1 (Remarks on code availability):  
  
README is pretty comprehensive.  
  
  
Reviewer \#2 (Remarks to the Author):  
  
Shen et al. developed Pansoma, a novel pangenome-aware, machine learning-based algorithm for detecting somatic single nucleotide variants (SNVs) and short insertions/deletions (indels) from tumor-only sequencing data. As tumor sequencing becomes increasingly common in clinical practice, matched normal samples are not always available, making it challenging to distinguish somatic mutations from germline variants. Pangenome and other graph-based reference genomes improve the representation of human germline variation, facilitating more accurate identification of true somatic mutations. Pansoma is the first somatic variant caller designed to operate directly on variation graph-aligned sequencing data. The method identifies candidate mismatch sites from graph-aligned reads, stores them in a novel node pileup (NPU) format, generates feature tensors containing sequence and alignment information, and uses a deep convolutional neural network to classify candidate variants as somatic SNVs or indels.  
  
The authors evaluated Pansoma using three benchmark datasets: COLO829T (a melanoma cell line), GIEB HG008T (a pancreatic tumor model), and the SMaHT Consortium HapMap Mixture benchmark, which was designed to assess performance on low variant allele frequency (VAF) variants. Across these datasets, Pansoma outperformed two commonly used tumor-only somatic variant callers, DeepSomatic-TO and ClairS-TO, primarily through improved precision, although this was accompanied by reduced recall. Additionally, a low-VAF version of Pansoma outperformed Mutect2 on the HapMap Mixture benchmark, highlighting its potential for improved detection of low-frequency somatic variants in tumor-only sequencing data.  
  
As the field continues towards adoption of graph-based references, the field needs graph-aware tools for somatic mutation detection. Pansoma is the first to help fill this gap. The manuscript is well written and clear. The manuscript would be strengthened by more discussion on some benchmarking choices as well as an expanded investigation on the reasons for the increased precision and decreased recall of Pansoma. More specific comments included below:   
  
Major:   
1\. The HapMap benchmark is a mixture of multiple (germline) samples at various frequencies, allowing for the detection of low-VAF events (as low as 0.5%). However, these are all germline variants, which are in theory removed by the pangenome alignments given that all of the HapMap components are included in the HPRC reference genome. The manuscript would benefit from more discussion as to why the HapMap benchmark is being included in the study as well as the limitations to this benchmarking approach for modelling somatic mutations.

**Response:** We thank the reviewer for this helpful comment. We apologize for not explaining the design and scope of the HapMap analysis clearly in the original manuscript. Although the standard HPRC graphs contain four (HPRC v1.1) or all six (HPRC v2.1) mixture donors, our analysis used a separate graph constructed without any of them. Alleles shared with other graph individuals can still be represented as paths, accounting for 44.8% of chr1 SNV and 20.8% of INDEL positives. Reads carrying these alleles can align without an edit, so they rarely become candidates. The same evaluation mask removes graph-carried alleles for both Pansoma and Mutect2. The retained benchmark focuses on alleles absent from the graph, which, like somatic mutations absent from the reference, must be detected as mismatches or INDELs against the graph.

We included this deeply sequenced mixture (about 500×) to examine low-allele-fraction detection. Although the mixture contains lower-VAF events, we restricted evaluation to truth VAF 1–16.5% because Pansoma's performance below 1% is limited; we now state this limitation explicitly. We also recognize that the positives are germline alleles of the minor donors and that mixing cell lines does not reproduce tumor mutational processes, clonal evolution or copy-number change. Filtering against our four resource VCFs would tag 95.7% of the chr1 SNV positives, so we do not apply this filter. The analysis therefore measures low-fraction allele detection rather than the separation of somatic from germline variants. We have moved the results to the Supplement, where we will provide additional technical details on the mixture, graph construction and evaluation. Somatic performance at lower purity is assessed using COLO829T tumor/normal dilutions.

Response Table 2.1

|                                                   | SNVs            | INDELs         |
| ------------------------------------------------- | --------------- | -------------- |
| truth alleles                                     | 392,567         | 124,499        |
| carried by the graph (left out of the evaluation) | 175,773 (44.8%) | 25,938 (20.8%) |
| scored                                            | 215,094         | 92,079         |

**Response Table 2.1.** HapMap mixture without the graph-carried truth alleles (chr1, benchmark regions, truth VAF 1–16.5%).

2\. While it is appreciated that the benchmarks used as quite standard, the manuscript would be strengthened by a brief discussion on how the ground truth datasets are generated for each benchmark, particularly, to strengthen the authors claim that the ground truth may contain germline variants and sequencing artifacts.

**Response:** We thank the reviewer for this helpful suggestion. The revised Methods explicitly describe the construction of the benchmark truth sets (Response Table 2.2) and their use in labelling graph-derived tensors. We have withdrawn the statement that the truth sets may contain germline variants and sequencing artifacts, because graph alignment alone cannot establish such errors (discussed in R2.3). Instead, we have added a brief discussion of truth alleles shared with matched-normal assemblies and how their interpretation depends on haplotype phase.

**Allele overlap with the matched normal.** We applied the same strict definition to all PASS truth alleles on chr1–22 in both tumors. Using both normal haplotypes reconstructed from assembly-based dipcall calls, we required one to equal GRCh38 plus the truth allele across the variant site and any tandem repeat containing it. This identified 144 HG008T alleles (8 SNVs, 136 INDELs; 0.8% of its truth alleles) and 45 COLO829T INDELs (2.4% of its truth INDELs; no SNVs). All but one were heterozygous in the normal.

At these heterozygous sites, one normal haplotype already carries the allele while the other does not. A change on the latter can therefore remain a somatic event under a haplotype-based truth definition, even though the resulting allele is already present on the other haplotype. Without phase information, the same allele can appear germline-like. We therefore describe this as an ambiguity of allele interpretation, rather than evidence of a truth-set error.

Response Table 2.2

| **Tumor**                   | **Somatic truth**                                                        | **How it was made and validated**                                                                              | **Regions** | **Reference**                                                               |
| --------------------------- | ------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------- | ----------- | --------------------------------------------------------------------------- |
| HG008T (pancreatic)         | GIAB HG008-T draft v0.2; 9,609 SNVs, 9,343 INDELs                        | comparison of the tumor (HG008-T) and normal (HG008-N v6.2) assemblies                                         | 2.81 Gb     | [Wagner et al., 2026, preprint](https://doi.org/10.64898/2026.05.01.722316) |
| COLO829T (melanoma)         | SMaHT COLO829 reference set (Ha et al., 2026); 44,005 SNVs, 2,059 INDELs | calls of four Illumina somatic callers, kept when PacBio reads show the ALT in the tumor and not in the normal | genome-wide | [Ha et al., 2026](https://doi.org/10.1016/j.xgen.2026.101390)               |
| HCC1395 (breast)            | SEQC2 v1.2.1, high and medium confidence; 41,072                         | 378 Illumina call sets (replicates × aligners × callers), checked with other technologies                      | 2.49 Gb     | [Fang et al., 2021](https://doi.org/10.1038/s41587-021-00993-6)             |
| HCC1954, H1437, H2009       | CASTLE; 26,578, 87,410, 168,397                                          | DeepSomatic calls on tumor–normal data of several technologies                                                 | 2.51 Gb     | [Park et al., 2026](https://doi.org/10.1038/s41587-025-02839-x)             |
| HapMap mixture (Supplement) | Kong et al. v1.6                                                         | minor donors' alleles absent from HG005, from donor assemblies                                                 | 2.59 Gb     | [Kong et al., 2026](https://doi.org/10.1016/j.xgen.2026.101358)             |

**Response Table 2.2.** Somatic benchmark truth sets, their construction, validation and references. 

3\. The authors show that Pansoma out-performs both ClairS-TO and DeepSomatic-TO by increasing precision which leads to improved F1-score overall. However, Pansoma has worse recall. The decrease in recall is \~20% in COLO829T and \~10% in HG008T from DeepSomatic. The authors suggest that some of this reduced recall is due to errors in the ground truth datasets. However, the explanation of the reduced recall could be expanded.  
\- The authors investigate the reasons behind false negatives for both HG008T and COLO829T in Supplementary Figure 2. They indicate some variants can be explained by limited allele frequency support in either the linear or the graph alignments, however, a large proportion of false negatives were labeled as “variant node not ranked \#1” and “first aligned node mismatch”. These labels need more explanation.  
\- Similarly, the authors propose that many false negatives overlap existing nodes in the pangenome reference suggesting they may correspond to germline alleles rather than somatic SNVs. To further show this, the authors could assess the alignment of the normal samples at the same sites and confirm that reads are mapping to the alternative paths preferentially. 

**Response:** We thank the reviewer for this important point. We agree that the recall deficit needed a clearer explanation and that graph alignment alone cannot establish errors in the truth sets. We withdrew that interpretation, re-scored all callers under one protocol and replaced the original false-negative categories with an allele-level analysis. We also use the matched normal resources to explore the calls off the reference path.

**Recall under a common protocol**. We revised the overall Pansoma pipeline and increased the loss weight assigned to true somatic examples during training, improving SNV recall. Although Pansoma's SNV recall is not consistently the highest across datasets, it is now broadly comparable to that of the comparators. We explored recall analysis (Response Figure 2.1, Response Figure 2.2), with full precision–recall curves for every caller (Response Figure 2.1). For each caller, we compared them with the truth set inside the benchmark regions using RTG vcfeval and traced the precision–recall curve by sweeping a threshold on the QUAL value of the caller's VCF ​respectively. The circle on each curve marks that caller's PASS call set after the same filter. The precision–recall curves show that Pansoma does not achieve a high F1 score by trading precision for recall. 

**The lost candidates.** We withdrew that interpretation and further examined the truth variants that receive no candidate (Response Figure 2.2). The extent to which somatic variants are absorbed by the graph varies considerably across datasets. For INDELs, we acknowledge that the graph absorbs a substantial fraction of candidates in COLO829T and HG008T (Response Figure 2.2b). However, as displayed in the below figures and tables, this does not critically limit performance, and Pansoma shows no significant loss of recall (Response Figure 2.1). Additionally, the majority of the absorbed INDEL are homopolymers and tandem repeats (Response Table 2.3).

**The normal at missed sites on graph nodes.** Additionally, we carefully examined truth alleles affected by complete or partial graph absorption. The results show that most of these false negatives are genuine somatic variants absorbed by the graph, which can be one reason for the low recall. We also compared this loss with existing population variant resources (four resource VCFs applied by Pansoma, DeepSomatic and ClairS) and found that the graph-specific loss is small: most of the loss is shared between the pangenome graph and the population filtering resources (Response Table 2.4). This explains why Pansoma still has a leading F1 score despite its comparatively low recall.

For the matched-normal analysis, we focused on the strict subset whose ALT is spelled exactly by an HPRC v2.1 path and which lose their candidate on at least one platform despite having one on the GRCh38-only graph (1,607 HG008T and 401 COLO829T alleles, chr1–22).

1.  **Normal reads, as the reviewer suggested.** We realigned COLO829BL normal ONT reads using the tumor's mapping settings. Only 1.1% of REF-like reads followed an ALT-like path, compared with 95% of ALT-like reads; REF-like normal reads were therefore rarely redirected onto ALT-like paths. At absorbed SNV loci, 2.2% of normal reads carried the ALT, compared with 45% of tumor reads. Read-level classification at repeat INDELs remained inconclusive.

2.  **Phased normal assemblies (HG008N, COLO829BL).** In this subset, the matched normal carries 110 HG008T alleles (6.8%) and 27 COLO829T alleles (6.7%); about 92% are absent from both normal haplotypes. At 1,153 of the 1,351 HG008T INDEL loci absent from the normal (85.3%), the HG008-T tumor assembly carries exactly the benchmark allele. Of the 401 COLO829T alleles, 347 are INDELs; 26 are carried by the normal, alongside the one normal-carried SNV. At 22 further INDELs, a normal haplotype carries an allele of the same length but different bases. At these 22 sites, the tumor hifiasm assembly carries the truth allele at 7 and only normal-haplotype sequences at 11; 4 remain unresolved. Counting the 11 normal-only cases and up to 4 unresolved cases as germline-like gives 37–41 of 347 INDELs (10.7–11.8%). This is an assembly-based classification, not a confirmed germline fraction: failure to recover a somatic allele in the assembly does not establish its absence from the tumor. The assembly uses the analyzed HiFi and ONT reads and therefore provides consistency evidence rather than independent validation

Response Figure 2.1

![](./media/image5.png)

**Response Figure 2.1. Precision–recall curves on chr1 for Pansoma, DeepSomatic-TO and ClairS-TO.** (a) SNVs and (b) INDELs in HCC1395T, COLO829T and HG008T (rows) on three platforms (columns), in the benchmark regions after the same population filter for every caller. Curves rank calls by each caller’s score (Pansoma: somatic probability; DeepSomatic-TO and ClairS-TO: QUAL, all records except RefCall and NoCall). Dots mark each caller’s PASS calls after the filter. Pansoma uses the six-tumor HPRC v2.1 models; chr1 of these three tumors was not used in training. Recall axes in b are cut to the range of the curves.

Response Figure 2.2

![](./media/image6.png)

Response Figure 2.2 | Why truth variants receive no candidate. Truth variants on chr1–22 without a candidate on the HPRC v2.1 graph (all truth alleles, inside and outside the benchmark regions), classified from the reads at each site. One pie per dataset; n, number of truth variants without a candidate; slices of 8% or more are labelled. Colors are defined separately for a and b. a, SNVs. S3: reads follow the graph branch carrying the ALT, leaving no edit. S4, repeats: (a) the SNV is written as an edit on another branch node; (b) the ALT haplotype is a graph path; (c) reads show no clear ALT. b, INDELs. I1: ALT reads follow a graph branch without an edit. I2: ALT reads follow a similar graph allele, retaining a residual edit with a different spelling.

Response Table 2.3

| **Truth INDELs (chr1–22)**                     | **HG008T**        | **COLO829T**    | **HCC1395T**   |
| ---------------------------------------------- | ----------------- | --------------- | -------------- |
| Absorbed by the graph (union across platforms) | 2,277             | 430             | 52             |
| Homopolymer ≥ 7 bp                             | 1,853 (81.4%)     | 241 (56.0%)     | 30 (57.7%)     |
| Short tandem repeat (unit 2–6 bp, ≥ 3 copies)  | 292 (12.8%)       | 83 (19.3%)      | 5 (9.6%)       |
| **Homopolymer + short tandem repeat, total**   | **2,145 (94.2%)** | **324 (75.3%)** | **35 (67.3%)** |
| Other repeats                                  | 6 (0.3%)          | 20 (4.7%)       | 1 (1.9%)       |
| Substitution within repetitive sequence        | 104 (4.6%)        | 16 (3.7%)       | 2 (3.8%)       |
| Others                                         | 21 (0.9%)         | 61 (14.2%)      | 12 (23.1%)     |
| Non-repetitive sequence                        | 1 (0.0%)          | 9 (2.1%)        | 2 (3.8%)       |

**Response Table 2.3.** Sequence context of somatic truth INDELs absorbed by the HPRC v2.1 graph (chr1–22). An INDEL is absorbed when, on at least one platform, its ALT reads align to a graph path without any edit, leaving no candidate. First row: absorbed INDELs and their percentage of all truth INDELs; other rows: percentages of the absorbed INDELs. Classes describe the somatic change relative to the matched-normal haplotype. Other repeats: shorter homopolymers (4–6 bp), longer-unit tandem repeats, or INDELs that do not change a repeat by whole units. Substitution within repetitive sequence: a base substitution on a germline length difference. Others: allele already on the normal haplotype, complex change, or unresolved normal haplotype.

Response Table 2.4 Somatic candidate loss on the HPRC graph

<table>
<thead>
<tr class="header">
<th><strong>% of all truth alleles (ranges over three platforms)</strong></th>
<th><strong>HG008T<br />
INDELs</strong></th>
<th><strong>COLO829T<br />
INDELs</strong></th>
<th><strong>HCC1395T<br />
INDELs</strong></th>
<th><strong>HG008T<br />
SNVs</strong></th>
<th><strong>COLO829T<br />
SNVs</strong></th>
<th><strong>HCC1395T<br />
SNVs</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>truth alleles, chr1–22</td>
<td>8,496</td>
<td>1,912</td>
<td>1,625</td>
<td>8,690</td>
<td>42,035</td>
<td>39,447</td>
</tr>
<tr class="even">
<td>somatic truth alleles tagged by population filtering</td>
<td>62.7%</td>
<td>52.6%</td>
<td>7.9%</td>
<td>10.5%</td>
<td>3.4%</td>
<td>2.8%</td>
</tr>
<tr class="odd">
<td>candidate on the GRCh38-only graph but not on HPRC v1.1</td>
<td>31.1–44.8%</td>
<td>25.8–35.1%</td>
<td>4.1–5.5%</td>
<td>3.4–3.9%</td>
<td>0.5–0.9%</td>
<td>0.6–0.8%</td>
</tr>
<tr class="even">
<td><strong>graph-specific loss left after filtering, HPRC v1.1</strong></td>
<td>3.8–7.2%</td>
<td>5.0–6.8%</td>
<td>1.3–1.9%</td>
<td>0.3–0.4%</td>
<td>0.2–0.6%</td>
<td>0.1–0.2%</td>
</tr>
<tr class="odd">
<td>candidate on the GRCh38-only graph but not on HPRC v2.1</td>
<td>29.4–42.1%</td>
<td>23.9–32.0%</td>
<td>4.0–4.9%</td>
<td>3.0–3.6%</td>
<td>0.5–0.9%</td>
<td>0.8–0.9%</td>
</tr>
<tr class="even">
<td><strong>graph-specific loss left after filtering, HPRC v2.1</strong></td>
<td>3.5–7.0%</td>
<td>5.4–6.9%</td>
<td>1.1–1.4%</td>
<td>0.2–0.4%</td>
<td>0.2–0.7%</td>
<td>0.1–0.3%</td>
</tr>
</tbody>
</table>

**Response Table 2.4.** Somatic SNVs and INDELs lost during graph-based candidate generation and their overlap with population filtering. Population filtering follows our resource VCFs: four for SNVs; gnomAD/CoLoRSdb only for INDELs: SNVs use gnomAD, dbSNP 138 (non-somatic records), the GATK 1000 Genomes PoN and CoLoRSdb; INDELs use gnomAD and CoLoRSdb only. A filter tag indicates that a somatic truth allele meets the filtering criteria, not that it is germline in the patient. Candidate presence requires an exact match to any candidate allele, including a non-representative allele. Percentages use the complete chr1–22 truth counts. Ranges span three platforms.

4\. In a similar vein, the authors suggest the majority of false-negatives are a result of intentionally excluding low-VAF variants (as shown in Supplementary Figure 2). Thresholding candidates at 10% read support may limit the detection of true somatic SNV with low VAF as many real-world somatic-mutation analyses will deal with low-VAF given low tumor purity / high heterogeneity.   
\- How do the authors propose to handle the detection of true low VAF somatic SNVs?   
\- The authors also developed a specialized low VAF model, however, there is limited description of how this model differs from their main model.   
\- Does the specialized low VAF model help recover low VAF variants? Is there a VAF threshold that the main model stops performing well on? It would be helpful if the authors conducted sensitivity analyses on decreasing VAF to demonstrate how sensitive their model is to VAF.  
\- The manuscript would also benefit from an expanded discussion on how they propose to use each model and when.

**Response:** We agree that low-VAF detection is a limitation of Pansoma, and we have narrowed the scope to focus on tumors rather than somatic mosaicism. The revised Pansoma model lowered the VAF to 0.05 as same as the compared tools (DeepSomatic-TO and ClairS-TO), and we also present the Pansoma evaluation in different VAF. We answer the reviewer’s questions in turn.

1.  **How are true low-VAF SNVs handled?** The revised manuscript states that Pansoma does not currently support somatic SNVs below an allele fraction of about 0.05: a candidate needs at least three ALT reads and an observed allele fraction of at least 0.05 (lowered from the submitted 0.10).

2.  **How does the low-VAF model differ?** The low-VAF HapMap model used tensor of up to 500 reads (we deeper the tensor height) with a candidate floor of 0.015 allele frequency on a graph excluding all mixture donor. Its positives are the minor donors' germline alleles, not somatic mutations.

3.  **Does it recover low-VAF somatic variants, and where does the main model fail?** The main model fails below the 0.05 floor (point 1). The low-VAF model was an exploratory test of detection in ultra-deep data (about 500×). Since the HapMap is a germline variant simulation dataset, it does not have the somatic related features, its performance on the mixture thus reflects low-fraction detection, not somatic–germline separation, so we report the HapMap analysis in the Supplement as exploratory rather than in the main benchmark.

4.  **Which model and when?** In the revised manuscript, the release models are the six-tumor models for both SNVs and INDELs. The HapMap model is an experimental option for cell-line mixtures only. As the reviewer suggested, we have expanded the discussion to describe when each model should be used in the revision.

5\. The authors investigate reasons for the increase in false negatives in Pansoma; however, it would be equally interesting to better understand the decrease in false positives. What is leading to the increased precision of Pansoma? Can the majority be explained by reduced germline contamination or are their other factors?  
\- A primary motivator for graph-aware somatic mutation detection is distinguishing germline variants from somatic mutations. However, beyond more accurate modeling of germline variation at the site itself, does graph-aware somatic mutation detection improve somatic mutation calling if there are multiple germline variants in the region? 

**Response:** We thank the reviewer for this suggestion. Our three-reference comparison supports three contributors to higher SNV precision: reduced germline contamination, fewer false calls in difficult regions, and improved calling near clustered germline variants. Their relative importance varies by dataset, and these categories overlap.

1.  **Reduced germline contamination***.* With the same reads and candidate rules, the HPRC graphs reduce matched-normal germline candidates by 84–91% relative to GRCh38-only while largely retaining somatic SNV candidates (Response Figure 2.3). This reduces germline contamination before population filtering.

2.  **Fewer false calls in difficult regions.** On COLO829T, the GRCh38-only graph's excess false calls are mainly from SMaHT extreme regions, where HPRC v2.1 reduces false positives and improves precision across all three platforms (Response Table 2.5).

3.  **The false positives of the GRCh38-only Pansoma model represented in the graph.** In the revision, we trained Pansoma using the same pipeline on the GRCh38-only graph (which contains only the GRCh38 path). We then compared the false positives of the GRCh38-only Pansoma model with the matched-normal germline variants and the population variants in the HPRC v2.1 graph. The results show that over 70% of the SNV false positives of the GRCh38-only Pansoma model are germline alleles represented in the HPRC v2.1 graph (Response Table 2.6).

4.  **Benefit near clustered germline variants.** For each candidate site, we counted other germline variants identified in the matched normal within 50 bp on either side, excluding any germline variant overlapping the candidate itself. We then grouped sites into those with no nearby germline variants, one, or at least two. On HG008T, the HPRC graphs give the largest SNV F1 gains where at least two germline variants are nearby (+0.107 to +0.204), mainly through higher precision (Response Table 2.7). The benefit therefore extends beyond germline variation at the candidate site itself. For INDELs, exact-allele candidate availability is lower on the graphs where germline variants cluster, and the results are more variable (Response Table 2.8).

Response Figure 2.3

![](./media/image4.png)

**Response Figure 2.3** We evaluated matched-normal germline alleles (5,061,892 for HG008T, 5,062,184 for COLO829T and 4,828,024 for HCC1395T) and somatic truth SNVs (8,690, 42,035 and 39,447, respectively) across chromosomes 1–22. A truth allele was considered retained if it was represented by a candidate with a matching representative allele. All truth records were included, irrespective of whether they fell within benchmark regions. The same sequencing reads, candidate-generation criteria and tensor construction code were used across all three graph references.

Response Table 2.5 Regional variant calling before population filtering

<table>
<thead>
<tr class="header">
<th><strong>chr1 dataset</strong></th>
<th><strong>HPRC graph</strong></th>
<th><strong>SMaHT extreme</strong><br />
<strong>PASS calls (precision)</strong></th>
<th><strong>SMaHT extreme</strong><br />
<strong>ΔF1</strong></th>
<th><strong>GIAB all difficult</strong><br />
<strong>ΔF1</strong></th>
<th><strong>GIAB low map or segdup</strong><br />
<strong>ΔF1</strong></th>
<th><strong>GIAB tandem</strong><br />
<strong>repeats and</strong><br />
<strong>homopolymers</strong><br />
<strong>ΔF1</strong></th>
<th><strong>GIAB outside difficult</strong><br />
<strong>ΔF1</strong></th>
<th><strong>INDEL F1</strong><br />
<strong>homopolymers</strong><br />
<strong>≥ 7 bp</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>HG008T<br />
Illumina</td>
<td>v1.1</td>
<td>15,955 (0.003)<br />
→ 601 (0.045)</td>
<td>+0.072<br />
[0.046, 0.099]</td>
<td>+0.147<br />
[0.120, 0.177]</td>
<td>+0.143<br />
[0.097, 0.193]</td>
<td>+0.103<br />
[0.073, 0.137]</td>
<td>+0.188<br />
[0.167, 0.210]</td>
<td>0.027<br />
→ 0.121</td>
</tr>
<tr class="even">
<td>HG008T<br />
Illumina</td>
<td>v2.1</td>
<td>15,955 (0.003)<br />
→ 599 (0.065)</td>
<td>+0.107<br />
[0.073, 0.144]</td>
<td>+0.207<br />
[0.172, 0.247]</td>
<td>+0.214<br />
[0.156, 0.279]</td>
<td>+0.133<br />
[0.094, 0.178]</td>
<td>+0.315<br />
[0.286, 0.344]</td>
<td>0.027<br />
→ 0.138</td>
</tr>
<tr class="odd">
<td>HG008T<br />
PacBio HiFi</td>
<td>v1.1</td>
<td>17,989 (0.003)<br />
→ 489 (0.070)</td>
<td>+0.111<br />
[0.077, 0.143]</td>
<td>+0.167<br />
[0.137, 0.198]</td>
<td>+0.144<br />
[0.100, 0.187]</td>
<td>+0.142<br />
[0.099, 0.186]</td>
<td>+0.186<br />
[0.165, 0.208]</td>
<td>0.027<br />
→ 0.046</td>
</tr>
<tr class="even">
<td>HG008T<br />
PacBio HiFi</td>
<td>v2.1</td>
<td>17,989 (0.003)<br />
→ 330 (0.112)</td>
<td>+0.169<br />
[0.121, 0.218]</td>
<td>+0.263<br />
[0.221, 0.308]</td>
<td>+0.198<br />
[0.139, 0.261]</td>
<td>+0.221<br />
[0.166, 0.279]</td>
<td>+0.313<br />
[0.285, 0.341]</td>
<td>0.027<br />
→ 0.056</td>
</tr>
<tr class="odd">
<td>HG008T<br />
ONT</td>
<td>v1.1</td>
<td>19,676 (0.003)<br />
→ 571 (0.058)</td>
<td>+0.094<br />
[0.064, 0.126]</td>
<td>+0.153<br />
[0.126, 0.184]</td>
<td>+0.128<br />
[0.088, 0.170]</td>
<td>+0.130<br />
[0.093, 0.171]</td>
<td>+0.179<br />
[0.159, 0.201]</td>
<td>0.010<br />
→ 0.010</td>
</tr>
<tr class="even">
<td>HG008T<br />
ONT</td>
<td>v2.1</td>
<td>19,676 (0.003)<br />
→ 669 (0.070)</td>
<td>+0.118<br />
[0.083, 0.159]</td>
<td>+0.195<br />
[0.159, 0.236]</td>
<td>+0.160<br />
[0.109, 0.227]</td>
<td>+0.148<br />
[0.107, 0.192]</td>
<td>+0.313<br />
[0.284, 0.343]</td>
<td>0.010<br />
→ 0.009</td>
</tr>
<tr class="odd">
<td>COLO829T<br />
Illumina</td>
<td>v1.1</td>
<td>8,143 (0.016)<br />
→ 3,102 (0.045)</td>
<td>+0.053<br />
[0.023, 0.165]</td>
<td>+0.101<br />
[0.052, 0.243]</td>
<td>+0.044<br />
[0.019, 0.158]</td>
<td>+0.042<br />
[0.019, 0.147]</td>
<td>+0.535<br />
[0.510, 0.557]</td>
<td>0.007<br />
→ 0.032</td>
</tr>
<tr class="even">
<td>COLO829T<br />
Illumina</td>
<td>v2.1</td>
<td>8,143 (0.016)<br />
→ 296 (0.449)</td>
<td>+0.471<br />
[0.399, 0.541]</td>
<td>+0.508<br />
[0.460, 0.554]</td>
<td>+0.463<br />
[0.381, 0.546]</td>
<td>+0.473<br />
[0.393, 0.542]</td>
<td>+0.629<br />
[0.613, 0.644]</td>
<td>0.007<br />
→ 0.041</td>
</tr>
<tr class="odd">
<td>COLO829T<br />
Fiber-seq</td>
<td>v1.1</td>
<td>5,684 (0.023)<br />
→ 844 (0.172)</td>
<td>+0.224<br />
[0.123, 0.369]</td>
<td>+0.286<br />
[0.193, 0.390]</td>
<td>+0.196<br />
[0.097, 0.362]</td>
<td>+0.323<br />
[0.268, 0.377]</td>
<td>+0.501<br />
[0.472, 0.526]</td>
<td>0.005<br />
→ 0.010</td>
</tr>
<tr class="even">
<td>COLO829T<br />
Fiber-seq</td>
<td>v2.1</td>
<td>5,684 (0.023)<br />
→ 321 (0.455)</td>
<td>+0.481<br />
[0.393, 0.570]</td>
<td>+0.503<br />
[0.447, 0.559]</td>
<td>+0.445<br />
[0.343, 0.561]</td>
<td>+0.415<br />
[0.323, 0.510]</td>
<td>+0.631<br />
[0.615, 0.646]</td>
<td>0.005<br />
→ 0.008</td>
</tr>
<tr class="odd">
<td>COLO829T<br />
ONT</td>
<td>v1.1</td>
<td>10,110 (0.012)<br />
→ 815 (0.171)</td>
<td>+0.241<br />
[0.168, 0.338]</td>
<td>+0.293<br />
[0.221, 0.375]</td>
<td>+0.217<br />
[0.134, 0.342]</td>
<td>+0.276<br />
[0.208, 0.349]</td>
<td>+0.478<br />
[0.454, 0.501]</td>
<td>0.002<br />
→ 0.002</td>
</tr>
<tr class="even">
<td>COLO829T<br />
ONT</td>
<td>v2.1</td>
<td>10,110 (0.012)<br />
→ 492 (0.274)</td>
<td>+0.348<br />
[0.254, 0.482]</td>
<td>+0.407<br />
[0.321, 0.500]</td>
<td>+0.303<br />
[0.207, 0.458]</td>
<td>+0.452<br />
[0.391, 0.511]</td>
<td>+0.595<br />
[0.580, 0.610]</td>
<td>0.002<br />
→ 0.001</td>
</tr>
</tbody>
</table>

**Response Table 2.5.** SNV calling by region on the HPRC graphs against the GRCh38-only graph (chr1, benchmark regions, before population filtering (no filtered by the four resource VCFs)). SMaHT extreme regions: PASS calls and precision on each graph and the paired ΔF1 (HPRC − GRCh38-only) with its 95% interval (2,000 resamples of 1-Mb blocks). GIAB v3.5 stratifications: all difficult regions (46.8 Mb of chr1), low mappability or segmental duplications, tandem repeats and homopolymers, and the complement of all difficult regions, with the same paired intervals. Last column: INDEL F1 in homopolymers of at least 7 bp, from the earlier INDEL models of the three-reference comparison. HP denotes homopolymers. In the ΔF1 columns, the second line gives the 95% CI. Arrows show GRCh38-only → HPRC.

Response Table 2.6

<table>
<thead>
<tr class="header">
<th><strong>Dataset</strong></th>
<th><strong>GRCh38-graph</strong><br />
<strong>FP before population filtering</strong></th>
<th><strong>Matched-normal germline</strong><br />
<strong>in HPRC v2.1 graph</strong><br />
<strong>(share; candidates retained on HPRC v2.1 model)</strong></th>
<th><strong>Matched-normal germline</strong><br />
<strong>not in HPRC v2.1 graph</strong><br />
<strong>(share; candidates retained on HPRC v2.1 model)</strong></th>
<th><strong>Non-germline FP</strong><br />
<strong>(% of GRCh38-</strong><br />
<strong>graph FP)</strong></th>
<th><strong>HPRC v2.1</strong><br />
<strong>FP before population filtering</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>HG008T Illumina</td>
<td>121,139</td>
<td>84%; 0.3%</td>
<td>13%; 85%</td>
<td>3%</td>
<td>2,870</td>
</tr>
<tr class="even">
<td>HG008T PacBio HiFi</td>
<td>130,720</td>
<td>85%; 0.4%</td>
<td>14%; 83%</td>
<td>1%</td>
<td>2,590</td>
</tr>
<tr class="odd">
<td>HG008T ONT</td>
<td>147,531</td>
<td>86%; 0.5%</td>
<td>13%; 85%</td>
<td>2%</td>
<td>3,038</td>
</tr>
<tr class="even">
<td>COLO829T Illumina</td>
<td>36,134</td>
<td>78%; 0.2%</td>
<td>9%; 87%</td>
<td>13%</td>
<td>986</td>
</tr>
<tr class="odd">
<td>COLO829T Fiber-seq</td>
<td>37,622</td>
<td>84%; 0.6%</td>
<td>10%; 88%</td>
<td>6%</td>
<td>969</td>
</tr>
<tr class="even">
<td>COLO829T ONT</td>
<td>37,199</td>
<td>71%; 0.6%</td>
<td>10%; 87%</td>
<td>18%</td>
<td>1,354</td>
</tr>
</tbody>
</table>

**Response Table 2.6.** SNV false positives before population filtering on GRCh38-graph and HPRC v2.1 (chr1, benchmark regions). The two germline columns group matched-normal germline alleles by whether they are represented in HPRC v2.1. Each cell gives the share of GRCh38-graph FP (%) and the percentage of those same alleles retained as candidates (forming tensors) on HPRC v2.1. Graph representation is inferred from carriage in at least 46 of 463 HPRC haplotypes. Non-germline gives the share of GRCh38-graph FP without a normal germline match; the first and final FP columns give the totals for each reference. Graph-represented germline alleles account for 71–86% of GRCh38-graph FP, but only 0.2–0.6% retain candidates, compared with 83–88% of graph-absent germline alleles.

Response Table 2.7

<table>
<thead>
<tr class="header">
<th><strong>Platform</strong></th>
<th><strong>Reference</strong></th>
<th><strong>False</strong><br />
<strong>positives</strong><br />
<strong>(0 / 1 / ≥2)</strong></th>
<th><strong>Precision</strong><br />
<strong>(0 / 1 / ≥2)</strong></th>
<th><strong>Recall</strong><br />
<strong>(0 / 1 / ≥2)</strong></th>
<th><strong>F1</strong><br />
<strong>(0 / 1 / ≥2)</strong></th>
<th><strong>ΔF1, HPRC −</strong><br />
<strong>GRCh38-only</strong><br />
<strong>(0 / 1 / ≥2)</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Illumina</td>
<td>GRCh38-only</td>
<td>693 / 160 / 236</td>
<td>0.403 / 0.307 / 0.041</td>
<td>0.832 / 0.732 / 0.256</td>
<td>0.543 / 0.433 / 0.070</td>
<td>–</td>
</tr>
<tr class="even">
<td></td>
<td>HPRC v1.1</td>
<td>515 / 114 / 73</td>
<td>0.490 / 0.390 / 0.141</td>
<td>0.881 / 0.753 / 0.308</td>
<td>0.629 / 0.514 / 0.194</td>
<td>+0.087 / +0.081 / +0.123</td>
</tr>
<tr class="odd">
<td></td>
<td>HPRC v2.1</td>
<td>590 / 136 / 66</td>
<td>0.460 / 0.367 / 0.154</td>
<td>0.897 / 0.814 / 0.308</td>
<td>0.608 / 0.506 / 0.205</td>
<td>+0.066 / +0.073 / +0.135</td>
</tr>
<tr class="even">
<td>PacBio HiFi</td>
<td>GRCh38-only</td>
<td>645 / 144 / 180</td>
<td>0.433 / 0.336 / 0.058</td>
<td>0.877 / 0.753 / 0.282</td>
<td>0.580 / 0.465 / 0.096</td>
<td>–</td>
</tr>
<tr class="odd">
<td></td>
<td>HPRC v1.1</td>
<td>488 / 77 / 18</td>
<td>0.503 / 0.490 / 0.357</td>
<td>0.879 / 0.763 / 0.256</td>
<td>0.639 / 0.597 / 0.299</td>
<td>+0.060 / +0.132 / +0.203</td>
</tr>
<tr class="even">
<td></td>
<td>HPRC v2.1</td>
<td>525 / 96 / 29</td>
<td>0.485 / 0.439 / 0.293</td>
<td>0.882 / 0.773 / 0.308</td>
<td>0.626 / 0.560 / 0.300</td>
<td>+0.047 / +0.095 / +0.204</td>
</tr>
<tr class="odd">
<td>ONT</td>
<td>GRCh38-only</td>
<td>650 / 129 / 125</td>
<td>0.435 / 0.371 / 0.081</td>
<td>0.891 / 0.784 / 0.282</td>
<td>0.584 / 0.503 / 0.126</td>
<td>–</td>
</tr>
<tr class="even">
<td></td>
<td>HPRC v1.1</td>
<td>540 / 97 / 23</td>
<td>0.478 / 0.419 / 0.343</td>
<td>0.881 / 0.722 / 0.308</td>
<td>0.619 / 0.530 / 0.324</td>
<td>+0.035 / +0.027 / +0.199</td>
</tr>
<tr class="odd">
<td></td>
<td>HPRC v2.1</td>
<td>610 / 126 / 52</td>
<td>0.455 / 0.379 / 0.188</td>
<td>0.909 / 0.794 / 0.308</td>
<td>0.607 / 0.513 / 0.233</td>
<td>+0.022 / +0.010 / +0.107</td>
</tr>
</tbody>
</table>

**Response Table 2.7.** SNV calling on HG008T by the number of nearby matched-normal germline variants (chr1, benchmark regions, PASS calls after population filtering). For each truth SNV and each call, we counted the germline variants of the matched normal (HG008-N, Dipcall) within 50 bp on either side, excluding any germline variant overlapping the site itself. Each cell gives the values for sites with 0 / 1 / ≥2 nearby germline variants. Recall is computed by the group of the truth SNV (561 / 97 / 39 truth SNVs) and precision by the group of the call. ΔF1 is the HPRC graph minus the GRCh38-only graph on the same platform, calculated before rounding. Single-tumor models trained on chr2–22 of HG008T on each reference, one training seed per model; PASS at the validation recall-0.9 threshold.

Response Table 2.8

<table>
<thead>
<tr class="header">
<th><strong>Platform</strong></th>
<th><strong>Reference</strong></th>
<th><strong>False</strong><br />
<strong>positives</strong><br />
<strong>(0 / 1 / ≥2)</strong></th>
<th><strong>Precision</strong><br />
<strong>(0 / 1 / ≥2)</strong></th>
<th><strong>Recall</strong><br />
<strong>(0 / 1 / ≥2)</strong></th>
<th><strong>F1</strong><br />
<strong>(0 / 1 / ≥2)</strong></th>
<th><strong>ΔF1, HPRC −</strong><br />
<strong>GRCh38-only</strong><br />
<strong>(0 / 1 / ≥2)</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Illumina</td>
<td>GRCh38-only</td>
<td>139 / 34 / 17</td>
<td>0.423 / 0.105 / 0.105</td>
<td>0.226 / 0.048 / 0.039</td>
<td>0.294 / 0.066 / 0.057</td>
<td>–</td>
</tr>
<tr class="even">
<td></td>
<td>HPRC v1.1</td>
<td>85 / 13 / 7</td>
<td>0.448 / 0.235 / 0.222</td>
<td>0.153 / 0.048 / 0.039</td>
<td>0.228 / 0.079 / 0.067</td>
<td>−0.067 / +0.014 / +0.010</td>
</tr>
<tr class="odd">
<td></td>
<td>HPRC v2.1</td>
<td>68 / 26 / 14</td>
<td>0.456 / 0.212 / 0.067</td>
<td>0.126 / 0.083 / 0.020</td>
<td>0.198 / 0.120 / 0.030</td>
<td>−0.097 / +0.054 / −0.027</td>
</tr>
<tr class="even">
<td>PacBio HiFi</td>
<td>GRCh38-only</td>
<td>325 / 137 / 79</td>
<td>0.246 / 0.062 / 0.060</td>
<td>0.235 / 0.107 / 0.098</td>
<td>0.240 / 0.078 / 0.074</td>
<td>–</td>
</tr>
<tr class="odd">
<td></td>
<td>HPRC v1.1</td>
<td>107 / 23 / 11</td>
<td>0.473 / 0.258 / 0.083</td>
<td>0.212 / 0.095 / 0.020</td>
<td>0.293 / 0.139 / 0.032</td>
<td>+0.053 / +0.061 / −0.042</td>
</tr>
<tr class="even">
<td></td>
<td>HPRC v2.1</td>
<td>90 / 12 / 5</td>
<td>0.486 / 0.400 / 0.167</td>
<td>0.188 / 0.095 / 0.020</td>
<td>0.271 / 0.154 / 0.035</td>
<td>+0.031 / +0.076 / −0.039</td>
</tr>
<tr class="odd">
<td>ONT</td>
<td>GRCh38-only</td>
<td>197 / 28 / 9</td>
<td>0.140 / 0.034 / 0.100</td>
<td>0.071 / 0.012 / 0.020</td>
<td>0.094 / 0.018 / 0.033</td>
<td>–</td>
</tr>
<tr class="even">
<td></td>
<td>HPRC v1.1</td>
<td>36 / 6 / 1</td>
<td>0.526 / 0.250 / 0.500</td>
<td>0.088 / 0.024 / 0.020</td>
<td>0.152 / 0.043 / 0.038</td>
<td>+0.058 / +0.026 / +0.005</td>
</tr>
<tr class="odd">
<td></td>
<td>HPRC v2.1</td>
<td>34 / 7 / 2</td>
<td>0.564 / 0.222 / 0.333</td>
<td>0.097 / 0.024 / 0.020</td>
<td>0.166 / 0.043 / 0.037</td>
<td>+0.072 / +0.025 / +0.004</td>
</tr>
</tbody>
</table>

**Response Table 2.8.** INDEL calling on HG008T by the number of nearby matched-normal germline variants (chr1, benchmark regions). For each truth INDEL and each call, we counted the germline variants of the matched normal (HG008-N, Dipcall) within 50 bp on either side, excluding any germline variant overlapping the site itself. Each cell gives the values for sites with 0 / 1 / ≥2 nearby germline variants. Recall is computed by the group of the truth INDEL (452 / 84 / 51 truth INDELs) and precision by the group of the call. ΔF1 is the HPRC graph minus the GRCh38-only graph on the same platform, calculated before rounding.

6\. The authors call somatic variants using tumour-only models, supplemented with filtering (using 4 PoNs) to remove known germline variants and sequencing artefacts. This filtering approach is standard to many tumor-only variant-calling approaches. However, the authors state numerous times that alignment to the pangenome (HPRC) should improve read mapping (reduce sequence artefacts) and reduce germline-variant contamination.   
\- Some of the PoNs used were developed using common benchmarking datasets and include variants from the tested samples \[the HapMap benchmark in particular\]. This may have inflated the reduced germline-variant contamination observed here as compared to real-world samples  
\- Fig2e shows that \>20% of all germline variants were not captured by the AF-HPRC graph. Were these variants successfully removed from the somatic variant calls using the PoNs?  
\- The authors only evaluate variant calling without the PoN for COLO829. Could the authors also include the performance of the Pansoma without the PoN for HG008T as well?

**Response:** We thank the reviewer for these questions.

1.  **Overlap between the filter resources and the test samples. **We use the four resources (gnomAD v2 genomes \[Karczewski et al., 2020\], dbSNP 138 \[Sherry et al., 2001\], the GATK 1000 Genomes panel of normals \[GATK resource bundle\] and CoLoRSdb v1.1.0 \[Lake and the CoLoRS Consortium, 2024\]) used by both ClairS-TO and DeepSomatic-TO. None of the six benchmark lines is listed among the samples of the 1000 Genomes panel or CoLoRSdb; gnomAD and dbSNP publish no sample lists. A resource built from an individual would contain that individual's private alleles, yet gnomAD holds only 0.6–2.5% and CoLoRSdb 4.0–7.9% of the benchmark normals' rare chr20 germline SNVs, the lowest of the 12 individuals tested. We find no evidence of overlapping. For HapMap mixtures, no population filtering was applied in the HapMap analysis. We apologize that this was not clearly stated in the original manuscript.

2.  **Germline variants absent from the graph.** About 13% of the matched normal's germline SNVs and 22–25% of its germline INDELs in the benchmark regions are absent from the HPRC graphs; they become candidates and are labelled germline in training. On HG008T chr1, the model rejects 79–88% of these SNV candidates and the filter removes 74–89% of those called somatic. Only 2.4–3.2% of the original germline SNV candidates remain after classification and filtering. These residual calls, representing rare germline alleles absent from both the graph and the resources, account for 83–98% of false calls after filtering.

3.  **Performance without population filtering (Response Tables 2.9 and 2.10).** On HG008T, single-tumor SNV F1 is 0.31–0.33 on HPRC v2.1 and about 0.01 on the GRCh38-only graph at similar recall; on COLO829T it rises from 0.07–0.10 to 0.61–0.67. Without their post-calling population filters, ClairS-TO reaches 0.006–0.028 and DeepSomatic-TO, which retains population allele frequencies as input, 0.056–0.404. Because Pansoma’s SNV calls are compared across references before post-calling population filtering, the observed HPRC benefit cannot be attributed to that filter or to overlap between its resources and the test samples; it is consistent with the preferential reduction of germline candidates on the graphs (Response Figure 2.3). The release six-tumor model reaches SNV F1 0.44–0.45 before filtering on HG008T, with higher precision but lower recall than the single-tumor models on the same HPRC v2.1 reference. INDELs show no consistent difference between references, and because their models take population allele frequencies as input, their unfiltered results retain population information.

Response Table 2.9

| **chr1 SNVs, no population filter (P / R / F1)** | **Pansoma GRCh38-only** | **Pansoma HPRC v1.1** | **Pansoma HPRC v2.1**  | **ClairS-TO**         | **DeepSomatic-TO**    |
| ------------------------------------------------ | ----------------------- | --------------------- | ---------------------- | --------------------- | --------------------- |
| HG008T Illumina                                  | 0.005 / 0.838 / 0.010   | 0.106 / 0.848 / 0.189 | 0.194 / 0.864 / 0.316  | 0.003 / 0.902 / 0.006 | 0.029 / 0.935 / 0.056 |
| HG008T PacBio HiFi                               | 0.005 / 0.862 / 0.010   | 0.098 / 0.856 / 0.175 | 0.207 / 0.852 / 0.332  | 0.003 / 0.894 / 0.006 | 0.059 / 0.829 / 0.109 |
| HG008T ONT                                       | 0.004 / 0.867 / 0.009   | 0.094 / 0.855 / 0.170 | 0.190 / 0.878 / 0.312  | 0.003 / 0.884 / 0.006 | 0.047 / 0.910 / 0.090 |
| COLO829T Illumina                                | 0.044 / 0.593 / 0.081   | 0.297 / 0.722 / 0.421 | 0.580 / 0.728 / 0.644  | 0.010 / 0.796 / 0.020 | 0.097 / 0.858 / 0.173 |
| COLO829T Fiber-seq                               | 0.052 / 0.616 / 0.096   | 0.435 / 0.733 / 0.545 | 0.613 / 0.734 / 0.668  | 0.011 / 0.771 / 0.021 | 0.255 / 0.569 / 0.352 |
| COLO829T ONT                                     | 0.038 / 0.565 / 0.071   | 0.416 / 0.708 / 0.524 | 0.532 / 0.722 / 0.612  | 0.010 / 0.709 / 0.020 | 0.152 / 0.714 / 0.251 |
| HCC1395T Illumina                                | –                       | –                     | 0.706 / 0.830 / 0.763ᵃ | 0.014 / 0.884 / 0.028 | 0.139 / 0.944 / 0.242 |
| HCC1395T PacBio HiFi                             | –                       | –                     | 0.639 / 0.817 / 0.717ᵃ | 0.014 / 0.845 / 0.028 | 0.257 / 0.944 / 0.404 |
| HCC1395T ONT                                     | –                       | –                     | 0.560 / 0.704 / 0.624ᵃ | 0.013 / 0.752 / 0.025 | 0.217 / 0.781 / 0.340 |

**Response Table 2.9.** SNVs without any population filter on the three references, with the comparators (chr1, benchmark regions, PASS calls; precision / recall / F1). Pansoma before our population filter (it has none of its own): HG008T and COLO829T, the single-tumor models of each reference (trained on that tumor's chr2–22), mean of three seeds (GRCh38-only: training-set statistics; HPRC v1.1 and v2.1: base quality standardised per dataset). ᵃ HCC1395T has tensors only on HPRC v2.1 and no single-tumor model: the six-tumor model (HCC1395T chr1 not in training, one run). ClairS-TO rerun without its panel-of-normals step; DeepSomatic-TO with its PON filter lifted (it still takes population allele frequencies as input).

Response Table 2.10

| **chr1 INDELs, no population filter (P / R / F1)** | **Pansoma GRCh38-only** | **Pansoma HPRC v1.1** | **Pansoma HPRC v2.1**  | **ClairS-TO**         | **DeepSomatic-TO**    |
| -------------------------------------------------- | ----------------------- | --------------------- | ---------------------- | --------------------- | --------------------- |
| HG008T Illumina                                    | 0.262 / 0.245 / 0.242   | 0.083 / 0.297 / 0.129 | 0.092 / 0.333 / 0.144  | 0.004 / 0.101 / 0.007 | 0.017 / 0.206 / 0.032 |
| HG008T PacBio HiFi                                 | 0.213 / 0.258 / 0.213   | 0.350 / 0.229 / 0.274 | 0.395 / 0.215 / 0.278  | 0.004 / 0.092 / 0.007 | 0.080 / 0.037 / 0.051 |
| HG008T ONT                                         | 0.019 / 0.361 / 0.035   | 0.025 / 0.233 / 0.044 | 0.011 / 0.300 / 0.020  | 0.004 / 0.087 / 0.007 | 0.043 / 0.019 / 0.026 |
| COLO829T Illumina                                  | 0.175 / 0.209 / 0.182   | 0.039 / 0.293 / 0.066 | 0.023 / 0.354 / 0.043  | 0.002 / 0.191 / 0.003 | 0.010 / 0.244 / 0.019 |
| COLO829T Fiber-seq                                 | 0.158 / 0.209 / 0.180   | 0.064 / 0.295 / 0.093 | 0.074 / 0.237 / 0.111  | 0.002 / 0.160 / 0.004 | 0.124 / 0.092 / 0.105 |
| COLO829T ONT                                       | 0.005 / 0.366 / 0.010   | 0.004 / 0.257 / 0.008 | 0.002 / 0.308 / 0.004  | 0.002 / 0.137 / 0.004 | 0.016 / 0.015 / 0.015 |
| HCC1395T Illumina                                  | –                       | –                     | 0.417 / 0.744 / 0.498ᵃ | 0.006 / 0.699 / 0.012 | 0.022 / 0.835 / 0.042 |
| HCC1395T PacBio HiFi                               | –                       | –                     | 0.515 / 0.699 / 0.589ᵃ | 0.005 / 0.526 / 0.011 | 0.271 / 0.647 / 0.382 |
| HCC1395T ONT                                       | –                       | –                     | 0.165 / 0.642 / 0.260ᵃ | 0.006 / 0.564 / 0.012 | 0.249 / 0.481 / 0.328 |

**Response Table 2.10.** INDELs without any population filter on the three references, with the comparators, as in Response Table 2.9. Pansoma: HG008T and COLO829T, single-tumor models with the final INDEL recipe on each reference, mean of three seeds; these models take gnomAD and CoLoRSdb allele frequencies as input, so, unlike the SNV models, their values without the filter are not free of population information. ᵃ HCC1395T has tensors only on HPRC v2.1: the final six-tumor INDEL models (HCC1395T chr1 and chr20 not in training), three seeds.

7\. Lines 292-294 state that graph-based alignment increased perfect aligned reads and reduce soft-clipping. Based on the referenced figure (Supplementary Figure 4), this appears overstated. Please provide actual read counts or show proportion of reads to improve clarity.

We agree that lines 292–294 overstated the alignment improvements. The original text read: “Across both datasets, graph alignment increased perfect alignments by up to 11% in WGS data and by more than four-fold in long-read datasets, while substantially reducing soft-clipped bases, substitutions and indel events (Supplementary Fig. 4).”

**Response:** We agree that reporting relative changes without absolute numbers was unclear. Supplementary Fig. 4 is replaced by absolute read counts and proportions for both tumors on the three platforms and four references (linear GRCh38, a GRCh38-only graph, HPRC v1.1 and v2.1; Response Figure 2.4 and Response Table 2.11).

  - **Illumina.** Perfectly aligned reads rise by about 10 percentage points on the HPRC graphs (HG008T: 1.87 of 2.74 billion reads, 68.3%, on linear GRCh38 against 2.12 of 2.70 billion, 78.5–78.7%; COLO829T: 69.9% against 79.8–80.5%). The GRCh38-only graph gives nearly the linear rates (68.7%, 70.0%), so the gain reflects the additional allelic sequence of the HPRC graphs.

  - **Long reads.** Perfectly aligned reads remain rare (PacBio 0.15–0.18% against 1.4–1.5%; ONT at most 2.8%), so the "four-fold" increase was on a very small base. Reads with a substitution fall (PacBio 95–96% to 66–75%; ONT 93–94% to 87–90%), but nearly every long read still carries an indel.

  - **Soft clipping** falls mainly between the linear and GRCh38-only alignments, that is, with the aligner rather than with the HPRC sequence; on long reads v2.1 is not consistently lower than v1.1.

The proportion of perfectly aligned reads can serve as one indicator of alignment quality, reflecting agreement with the traversed reference sequence. However, perfect alignment alone does not establish correct read placement or improved variant calling. We therefore complemented these statistics with analyses of read placement and downstream calling performance using simulated reads (Response Table 2.12).

Response Figure 2.4 Reads on chr1–22 on four references

![Supplementary Figure R2.7: Reads on chr1–22 on four references](./media/image7.png)

**Response Figure 2.4.** Reads aligned to chr1–22 for each dataset and reference: the share perfectly aligned (no substitution, indel or soft clipping) and the shares with at least one substitution, indel or soft-clipped end. Rows: Illumina; PacBio (HG008T HiFi, COLO829T Fiber-seq); ONT. Note the different y-axis ranges. Proportions are computed among the reads placed on chr1–22 on each reference. These totals differ by up to 7% between references, because each reference leaves a different number of reads unmapped, on other contigs or only on non-reference nodes; all totals are listed in Response Table 2.11. The three graphs are aligned with vg Giraffe; linear GRCh38 uses the data providers' alignments.

Response Table 2.11

| **dataset**               | **reference**   | **reads on chr1–22** | **perfectly aligned**  | **with substitution** | **with indel**      | **with soft clipping** |
| ------------------------- | --------------- | -------------------- | ---------------------- | --------------------- | ------------------- | ---------------------- |
| HG008T Illumina           | GRCh38 (linear) | 2,739,599,264        | 1,869,903,341 (68.25%) | 810,918,802 (29.60%)  | 78,905,626 (2.88%)  | 96,334,574 (3.52%)     |
| HG008T Illumina           | GRCh38 graph    | 2,710,092,615        | 1,861,828,942 (68.70%) | 800,943,101 (29.55%)  | 79,136,287 (2.92%)  | 74,218,393 (2.74%)     |
| HG008T Illumina           | HPRC v1.1 graph | 2,704,033,209        | 2,121,622,695 (78.46%) | 560,922,772 (20.74%)  | 21,774,798 (0.81%)  | 66,387,175 (2.46%)     |
| HG008T Illumina           | HPRC v2.1 graph | 2,695,745,607        | 2,121,713,347 (78.71%) | 552,985,555 (20.51%)  | 19,851,148 (0.74%)  | 63,829,495 (2.37%)     |
| COLO829T Illumina         | GRCh38 (linear) | 2,296,226,777        | 1,604,645,813 (69.88%) | 627,386,769 (27.32%)  | 72,684,106 (3.17%)  | 95,094,133 (4.14%)     |
| COLO829T Illumina         | GRCh38 graph    | 2,283,443,818        | 1,599,146,760 (70.03%) | 629,366,832 (27.56%)  | 71,551,146 (3.13%)  | 78,056,743 (3.42%)     |
| COLO829T Illumina         | HPRC v1.1 graph | 2,285,859,727        | 1,823,482,724 (79.77%) | 423,169,962 (18.51%)  | 23,805,599 (1.04%)  | 78,789,217 (3.45%)     |
| COLO829T Illumina         | HPRC v2.1 graph | 2,265,562,422        | 1,822,701,086 (80.45%) | 410,679,337 (18.13%)  | 20,250,428 (0.89%)  | 66,561,782 (2.94%)     |
| HG008T PacBio HiFi        | GRCh38 (linear) | 12,931,641           | 18,818 (0.15%)         | 12,264,874 (94.84%)   | 12,862,468 (99.47%) | 961,155 (7.43%)        |
| HG008T PacBio HiFi        | GRCh38 graph    | 12,675,374           | 18,819 (0.15%)         | 12,020,862 (94.84%)   | 12,606,854 (99.46%) | 314,848 (2.48%)        |
| HG008T PacBio HiFi        | HPRC v1.1 graph | 12,689,479           | 179,824 (1.42%)        | 8,558,774 (67.45%)    | 12,377,219 (97.54%) | 228,857 (1.80%)        |
| HG008T PacBio HiFi        | HPRC v2.1 graph | 12,977,109           | 184,297 (1.42%)        | 8,605,339 (66.31%)    | 12,666,519 (97.61%) | 422,918 (3.26%)        |
| COLO829T PacBio Fiber-seq | GRCh38 (linear) | 25,308,248           | 44,846 (0.18%)         | 24,325,270 (96.12%)   | 25,138,324 (99.33%) | 2,529,764 (10.00%)     |
| COLO829T PacBio Fiber-seq | GRCh38 graph    | 25,026,903           | 44,873 (0.18%)         | 24,011,605 (95.94%)   | 24,858,107 (99.33%) | 936,884 (3.74%)        |
| COLO829T PacBio Fiber-seq | HPRC v1.1 graph | 25,070,732           | 386,300 (1.54%)        | 18,732,761 (74.72%)   | 24,353,478 (97.14%) | 717,920 (2.86%)        |
| COLO829T PacBio Fiber-seq | HPRC v2.1 graph | 26,066,512           | 397,244 (1.52%)        | 19,343,766 (74.21%)   | 25,350,225 (97.25%) | 1,538,959 (5.90%)      |
| HG008T ONT                | GRCh38 (linear) | 3,208,600            | 43,815 (1.37%)         | 3,009,593 (93.80%)    | 3,084,419 (96.13%)  | 1,215,595 (37.89%)     |
| HG008T ONT                | GRCh38 graph    | 3,100,754            | 43,748 (1.41%)         | 2,869,956 (92.56%)    | 2,950,998 (95.17%)  | 930,380 (30.00%)       |
| HG008T ONT                | HPRC v1.1 graph | 3,099,600            | 85,847 (2.77%)         | 2,691,454 (86.83%)    | 2,929,825 (94.52%)  | 913,057 (29.46%)       |
| HG008T ONT                | HPRC v2.1 graph | 3,087,083            | 84,623 (2.74%)         | 2,704,935 (87.62%)    | 2,944,975 (95.40%)  | 899,235 (29.13%)       |
| COLO829T ONT              | GRCh38 (linear) | 12,636,129           | 723 (0.01%)            | 11,794,165 (93.34%)   | 11,928,948 (94.40%) | 12,590,223 (99.64%)    |
| COLO829T ONT              | GRCh38 graph    | 11,839,534           | 778 (0.01%)            | 10,956,152 (92.54%)   | 11,047,335 (93.31%) | 11,723,760 (99.02%)    |
| COLO829T ONT              | HPRC v1.1 graph | 11,877,368           | 1,121 (0.01%)          | 10,662,112 (89.77%)   | 11,030,919 (92.87%) | 11,760,246 (99.01%)    |
| COLO829T ONT              | HPRC v2.1 graph | 11,789,372           | 1,098 (0.01%)          | 10,645,364 (90.30%)   | 11,021,767 (93.49%) | 11,674,637 (99.03%)    |

**Response Table 2.11.** Counts and proportions of reads aligned to chr1–22 on four references.

*Methods.*

1\. *Alignments.* Each dataset's tumor reads were aligned with vg Giraffe to the GRCh38-only graph and to the HPRC v1.1 (d9) and v2.1 (d46) graphs. Linear GRCh38 is the data providers' alignment: for HG008T (GIAB) BWA-MEM 0.7.17, pbmm2 1.10.0 and minimap2 2.26; for COLO829T (SMaHT) BWA-MEM 0.7.17 in Sentieon 202308.01 (with duplicate marking, indel realignment and base recalibration), pbmm2 1.13.0 and minimap2 2.26.

2\. *Reads counted.* Every read of the dataset, each mate of a pair separately, whose alignment is on chr1–22: on linear GRCh38 the primary alignment (secondary and supplementary alignments are not counted); on a graph, the GRCh38 chromosome holding most of the read's reference-path bases. Unmapped reads, reads on other contigs and reads aligned only to non-reference nodes are not counted.

3\. *Measures.* Perfectly aligned: aligned end to end with no substitution, indel or soft clipping. Substitution: at least one mismatched base (linear: NM minus inserted and deleted bases; graph: mismatches against the sequence of the read's graph path, so a match to a non-reference allele is not a substitution). Indel: at least one inserted or deleted base. Soft clipping: read bases at either end not aligned to the reference (linear: CIGAR S in the primary alignment; hard clipping occurs only in supplementary alignments, which are not counted; graph: read bases at an alignment end not aligned to the path).

4\. On linear GRCh38 a long read split by the aligner keeps its remainder as a supplementary alignment, so its primary alignment is soft-clipped; Giraffe writes one alignment per read.

Response Table 2.12

|                                           | **GRCh38 (linear ref)** | **GRCh38-graph** | **CHM13-graph** | **HPRC v1.1** | **HPRC v2.1** |
| ----------------------------------------- | ----------------------- | ---------------- | --------------- | ------------- | ------------- |
| **HG008-N**                               |                         |                  |                 |               |               |
| Illumina, perfectly aligned               | 63.7%                   | 63.7%            | 64.8%           | 73.0%         | 73.1%         |
| Illumina, misplaced (% of placed)         | 1.23%                   | 1.19%            | 1.12%           | 0.95%         | 0.90%         |
| Illumina, misplaced at MAPQ ≥ 30 (reads)  | 5,555                   | 4,240            | 4,574           | 823           | 567           |
| HiFi-like, misplaced                      | 0.51%                   | 0.47%            | 0.41%           | 0.31%         | 0.36%         |
| HiFi-like, misplaced at MAPQ ≥ 30 (reads) | 57                      | 11               | 13              | 3             | 1             |
| ONT-like, misplaced                       | 0.52%                   | 0.75%            | 0.68%           | 0.60%         | 1.45%         |
| ONT-like, misplaced at MAPQ ≥ 30 (reads)  | 45                      | 164              | 146             | 127           | 294           |
| **COLO829BL**                             |                         |                  |                 |               |               |
| Illumina, perfectly aligned               | 63.8%                   | 63.8%            | 64.8%           | 73.0%         | 73.1%         |
| Illumina, misplaced (% of placed)         | 1.13%                   | 1.09%            | 1.06%           | 0.88%         | 0.85%         |
| Illumina, misplaced at MAPQ ≥ 30 (reads)  | 5,491                   | 4,386            | 5,625           | 1,368         | 1,216         |
| HiFi-like, misplaced                      | 0.40%                   | 0.38%            | 0.42%           | 0.28%         | 0.34%         |
| HiFi-like, misplaced at MAPQ ≥ 30 (reads) | 31                      | 11               | 31              | 4             | 4             |
| ONT-like, misplaced                       | 0.42%                   | 0.64%            | 0.66%           | 0.56%         | 1.50%         |
| ONT-like, misplaced at MAPQ ≥ 30 (reads)  | 36                      | 153              | 145             | 129           | 315           |

**Response Table 2.12.** Placement accuracy of reads simulated from chr1 and chr20 of the phased diploid assemblies HG008N v6.2 and COLO829BL (assembled with Verkko v2.1). Both samples used identical simulation settings, generating 2.5 million Illumina read pairs (2 × 150 bp), 40,000 HiFi-like reads and 40,000 ONT-like reads per sample. Long-read template lengths were normally distributed (mean 15 kb; SD 3 kb for HiFi-like and 4 kb for ONT-like), truncated to 5–30 kb and 3–30 kb, respectively, before sequencing errors were simulated. Templates were sampled from one-to-one assembly-to-GRCh38 intervals; evaluation used reads with unambiguous truth positions on both GRCh38 and CHM13. Linear GRCh38 uses BWA-MEM for Illumina and minimap2 for long reads. GRCh38-graph and CHM13-graph were built from the corresponding FASTA sequences using vg autoindex; all graph arms use vg Giraffe. Misplacement means an incorrect chromosome or both alignment endpoints more than 50 bp from their corresponding true endpoints. Misplacement percentages use placed reads as the denominator; MAPQ ≥ 30 rows give read counts. In both samples, HPRC improves Illumina and HiFi-like placement relative to all three single-genome controls.

Minor:   
1\. Are the variants shown in Figure 2 (panels d and f) generated from Pansoma or from previous linear-reference analyses? If these are based off the ground-truth variant calls, it is not clear they are added new information and could be moved to supplementary.

Neither panel uses Pansoma output; both summarize the published benchmark truth sets. Figure 2d shows the allele-fraction distributions of truth SNVs and INDELs from COLO829T, HG008T and the HapMap mixture, while Figure 2f shows the SBS96 spectra and COSMIC signature decomposition of COLO829T and HG008T truth SNVs. As suggested, both panels have been moved to Supplementary, whose legend identifies the truth sets and their sources.

2\. Figure 2f is referred to before Figure 2e in the main text 

After panels d and f were moved to the Supplementary Information, the remaining Figure 2 panels were relabeled in the order of their first citation. The text and legend have been updated accordingly.

3\. In Supplementary Figure 2, the graph corresponding to COLO829T is not clearly labeled

In the revised Supplementary Fig. 2, which classifies truth alleles lost during candidate generation, each chart is explicitly labeled with its tumor, sequencing platform and variant type.

4\. The purity and VAF labels in Supplementary Figure 1 are often overlapping the points which make it hard to read

In the revised Supplementary Fig. 1, purity and allele-fraction information appears in the legend rather than beside individual points, eliminating overlapping annotations.

5\. What do the authors mean by “truth overlapping” on page 13

“Truth-overlapping calls” (page 13, line 349) meant PASS SNV calls matching a benchmark truth allele, rather than calls merely sharing its genomic position. The “shared core” refers to the subset of these true positives detected by all three callers. The revised text uses “true positives” and describes the shared core separately.

6\. In Figure 4a, the authors show the cosine similarity between the signature decomposition of each tool with the signature decomposition of the benchmark truth set. It would be helpful to also include the signature composition of the truth set in the middle bar plots as well. 

We agree that the truth should be displayed directly alongside the callers. The revised Figure 4a replaces the signature-decomposition bars with SBS96 spectra and also include the signature composition of the truth set. For each tumor and sequencing platform, it shows the chr1 somatic truth spectrum beside those of Pansoma, DeepSomatic-TO and ClairS-TO, together with the matched-normal germline spectrum. Each caller’s cosine similarity to the truth is calculated from these SBS96 profiles. This presentation allows direct comparison across all 96 substitution contexts and shows which contexts contribute to agreement or disagreement.

7\. The order of the figures in Figure 4 does not match how they are referenced in the main text which at times makes it confusing to read

The revised Figure 4 panels are labeled in the order of their first citation, with corresponding updates to the text and legend.

8\. In Supplementary Figure 5, what is the difference between "candidate sites" and "candidates"?

In the original Supplementary Fig. 5, a *candidate site* was a position containing at least one passing variant allele, while a *candidate* was an individual passing allele. Thus, a position containing an SNV and a deletion counted as one site and two candidates. The revised figure follows the revised pipeline’s grouping: each candidate is an allele defined by its graph node, position, REF and ALT; each candidate site groups passing alleles at the same node position within one variant class, SNV or INDEL, and produces one tensor. Consequently, SNV and INDEL candidates at the same position are counted as separate candidate sites. The Methods and legend now explicitly define these counting units.

Reviewer \#2 (Remarks on code availability):  
  
The model has a series of dependencies that have prevented full testing, however, we are continuing to troubleshoot to get the model running on our platform.  
  
  
Reviewer \#3 (Remarks to the Author):  
  
This manuscript presents Pansoma, a machine learning-based framework for tumor-only somatic variant calling that operates directly on pangenome graph alignments. The authors introduce a node-centric pileup representation, construct graph-derived tensors for machine-learning classification, and evaluate the framework using COLO829T, HG008T and SMaHT HapMap mixture datasets across short- and long-read sequencing platforms.  
Addressing the primary limitation of tumor-unmatched somatic variant detection, specifically, the discrimination of true somatic events from co-occurring germline alleles and alignment noise, this investigation evaluates the utility of a pangenome graph reference architecture. Incorporating population genetics directly into the alignment workflow rather than relying solely on downstream filtering rules provides a principled computational strategy. The authors substantiate this approach with concrete performance metrics: a sizable subset of benchmark germline variation is directly represented in the AF-HPRC graph, graph-guided alignment minimizes read-mapping errors and suppresses false-positive candidate calls, and the Pansoma workflow consistently outperforms existing methods in SNV calling precision and overall F1 performance across multiple validation benchmarks.  
Overall, I find Pansoma technically interesting and potentially important. To my knowledge, while graph-based alignment methods are well established, an end-to-end graph-native framework specifically developed for somatic variant calling represents a meaningful methodological advance.  
My main concern is not whether Pansoma performs well; the benchmarks indicate that it does. It is whether the manuscript demonstrates that this advantage is specifically attributable to the pangenome-native representation, the study's central conceptual advance. For this reason, I recommend major revision.  
  
Major comments  
1\. The specific contribution of the pangenome warrants more direct isolation  
The study articulates a logical mechanistic framework wherein the pangenome structure accommodates baseline germline variants via alternative haplotypes. This graph-guided alignment minimizes non-biological discrepancies and mapping errors, thereby suppressing upstream candidate noise and translating directly into improved variant precision. Several components of this model are well supported by the current data-graph-based alignment increases the proportion of perfectly aligned reads and reduces substitutions, indels, and soft-clipping, while graph-informed candidate generation substantially lowers the number of putative variants.  
Nevertheless, the end-to-end benchmarking does not cleanly disentangle the contribution of the pangenome per se. Pansoma diverges from ClairS-TO and DeepSomatic-TO not only in its reference representation but also in alignment strategy, candidate generation, tensorization, model architecture, and post-processing. As the Methods note, the comparator callers employ their standard GRCh38 workflows, whereas Pansoma relies on vg graph alignment coupled to a bespoke graph-native pipeline. To definitively validate the manuscript's primary claim, a rigorous reference-level ablation study is strongly recommended. Specifically, the authors should compare the Pansoma pipeline's performance using a GRCh38 linear reference (or single-path graph equivalent) directly against the full HPRC pangenome graph, ensuring that candidate extraction, tensor generation, deep learning architecture, training setup, and filtering thresholds remain identical across both arms. Such an analysis is critical to elucidate graph-specific benefits from algorithmic optimizations, thereby demonstrating how much of the performance advance is directly driven by the pangenome architecture.

**Response:** We thank the reviewer and agree that comparisons with other callers alone cannot isolate the contribution of the pangenome. We therefore performed the three-reference comparison, comparing a GRCh38-only graph with the allele-filtered HPRC v1.1 (90 haplotypes) and HPRC v2.1 (464 haplotypes) graphs. For HG008T and COLO829T, each across three sequencing platforms, the same reads were aligned with vg Giraffe to each reference. We used the exact same Pansoma pipeline across references, except that the haplotype-count channel was scaled to each graph's haplotype number and was constant on the GRCh38-only graph. Models were trained separately for each reference with three seeds, using chr2–22 of the same tumor for training and holding out chr1 for testing. To further assess the impact of the pangenome graph, we examined performance before applying the post-calling population filter, with SNV and INDEL results reported in Response Tables 3.3 and 3.4, respectively. SNV results after population filtering are reported in Response Table 3.1, with paired confidence intervals in Response Table 3.2. After the same population filter, both HPRC graphs outperform GRCh38-only on all six datasets, increasing mean SNV F1 by 0.028–0.292, with higher precision and recall (Response Table 3.1).

1.  **Improvement before population filtering***.* Both HPRC graphs also outperform GRCh38-only before filtering, with HPRC v2.1 giving the highest SNV F1 on all six datasets. HPRC v2.1 raises F1 from approximately 0.01 to 0.31–0.33 on HG008T and from 0.07–0.10 to 0.61–0.67 on COLO829T. The advantage therefore cannot be attributed to downstream population filtering (Response Table 3.3).

2.  **Germline candidate reduction.** Relative to GRCh38-only, germline candidates decrease by 84–90% on HPRC v1.1 and 85–91% on HPRC v2.1, while most somatic SNV candidates are retained (Response Figure 3.1). This reduction occurs before model classification and population filtering, providing a mechanistic link between graph representation and improved SNV precision. The advantage also persists after removing the haplotype-count channel.

Response Figure 3.1

![](./media/image4.png)

**Response Figure 3.1** We evaluated matched-normal germline alleles (5,061,892 for HG008T, 5,062,184 for COLO829T and 4,828,024 for HCC1395T) and somatic truth SNVs (8,690, 42,035 and 39,447, respectively) across chromosomes 1–22. A truth allele was considered retained if it was represented by a candidate with a matching representative allele. All truth records were included, irrespective of whether they fell within benchmark regions. The same sequencing reads, candidate-generation criteria and tensor construction code were used across all three graph references.

Response Table 3.1

| **tumor** | **platform** | **reference**                               | **precision**                         | **recall**                            | **F1**                                | **HPRC \> GRCh38-only, same seed** |
| --------- | ------------ | ------------------------------------------- | ------------------------------------- | ------------------------------------- | ------------------------------------- | ---------------------------------- |
| HG008T    | Illumina     | GRCh38-only                                 | 0.334 / 0.394 / 0.388 (0.372 ± 0.033) | 0.786 / 0.783 / 0.789 (0.786 ± 0.003) | 0.469 / 0.524 / 0.520 (0.504 ± 0.030) | \-                                 |
| HG008T    | Illumina     | HPRC v1.1 (d9)                              | 0.452 / 0.458 / 0.460 (0.457 ± 0.005) | 0.831 / 0.819 / 0.801 (0.817 ± 0.015) | 0.585 / 0.587 / 0.585 (0.586 ± 0.001) | 3 / 3                              |
| HG008T    | Illumina     | HPRC v2.1 (d46)                             | 0.429 / 0.457 / 0.432 (0.439 ± 0.016) | 0.852 / 0.828 / 0.815 (0.832 ± 0.019) | 0.570 / 0.589 / 0.565 (0.575 ± 0.013) | 3 / 3                              |
| HG008T    | PacBio HiFi  | GRCh38-only                                 | 0.372 / 0.472 / 0.463 (0.436 ± 0.055) | 0.826 / 0.792 / 0.805 (0.808 ± 0.017) | 0.513 / 0.592 / 0.588 (0.564 ± 0.044) | \-                                 |
| HG008T    | PacBio HiFi  | HPRC v1.1 (d9)                              | 0.479 / 0.443 / 0.464 (0.462 ± 0.018) | 0.829 / 0.808 / 0.831 (0.823 ± 0.013) | 0.608 / 0.572 / 0.595 (0.592 ± 0.018) | 2 / 3                              |
| HG008T    | PacBio HiFi  | HPRC v2.1 (d46)                             | 0.473 / 0.477 / 0.487 (0.479 ± 0.007) | 0.839 / 0.822 / 0.821 (0.827 ± 0.010) | 0.605 / 0.604 / 0.611 (0.607 ± 0.004) | 3 / 3                              |
| HG008T    | ONT          | GRCh38-only                                 | 0.394 / 0.368 / 0.394 (0.385 ± 0.015) | 0.842 / 0.826 / 0.780 (0.816 ± 0.032) | 0.537 / 0.509 / 0.524 (0.523 ± 0.014) | \-                                 |
| HG008T    | ONT          | HPRC v1.1 (d9)                              | 0.473 / 0.451 / 0.447 (0.457 ± 0.014) | 0.821 / 0.816 / 0.822 (0.820 ± 0.003) | 0.600 / 0.581 / 0.579 (0.587 ± 0.012) | 3 / 3                              |
| HG008T    | ONT          | HPRC v2.1 (d46)                             | 0.432 / 0.467 / 0.474 (0.457 ± 0.022) | 0.859 / 0.824 / 0.809 (0.831 ± 0.026) | 0.575 / 0.596 / 0.598 (0.589 ± 0.013) | 3 / 3                              |
| COLO829T  | Illumina     | GRCh38-only                                 | 0.391 / 0.289 / 0.499 (0.393 ± 0.105) | 0.583 / 0.583 / 0.564 (0.576 ± 0.011) | 0.468 / 0.386 / 0.530 (0.461 ± 0.072) | \-                                 |
| COLO829T  | Illumina     | HPRC v1.1 (d9), validation-selected epoch 8 | 0.484 / 0.562 / 0.475 (0.507 ± 0.048) | 0.698 / 0.705 / 0.702 (0.702 ± 0.004) | 0.572 / 0.626 / 0.566 (0.588 ± 0.033) | 3 / 3                              |
| COLO829T  | Illumina     | HPRC v2.1 (d46)                             | 0.867 / 0.704 / 0.820 (0.797 ± 0.084) | 0.715 / 0.705 / 0.717 (0.712 ± 0.007) | 0.783 / 0.704 / 0.765 (0.751 ± 0.042) | 3 / 3                              |
| COLO829T  | Fiber-seq    | GRCh38-only                                 | 0.745 / 0.841 / 0.868 (0.818 ± 0.065) | 0.625 / 0.600 / 0.576 (0.600 ± 0.024) | 0.680 / 0.700 / 0.692 (0.691 ± 0.010) | \-                                 |
| COLO829T  | Fiber-seq    | HPRC v1.1 (d9)                              | 0.785 / 0.876 / 0.845 (0.835 ± 0.046) | 0.723 / 0.707 / 0.709 (0.713 ± 0.009) | 0.752 / 0.782 / 0.771 (0.768 ± 0.015) | 3 / 3                              |
| COLO829T  | Fiber-seq    | HPRC v2.1 (d46)                             | 0.866 / 0.846 / 0.849 (0.854 ± 0.011) | 0.707 / 0.713 / 0.737 (0.719 ± 0.016) | 0.779 / 0.774 / 0.789 (0.780 ± 0.008) | 3 / 3                              |
| COLO829T  | ONT          | GRCh38-only                                 | 0.340 / 0.353 / 0.475 (0.390 ± 0.075) | 0.543 / 0.540 / 0.569 (0.551 ± 0.016) | 0.418 / 0.427 / 0.518 (0.454 ± 0.055) | \-                                 |
| COLO829T  | ONT          | HPRC v1.1 (d9)                              | 0.793 / 0.816 / 0.833 (0.814 ± 0.020) | 0.701 / 0.693 / 0.673 (0.689 ± 0.014) | 0.744 / 0.749 / 0.745 (0.746 ± 0.003) | 3 / 3                              |
| COLO829T  | ONT          | HPRC v2.1 (d46)                             | 0.788 / 0.734 / 0.744 (0.755 ± 0.029) | 0.692 / 0.703 / 0.718 (0.704 ± 0.013) | 0.737 / 0.718 / 0.730 (0.729 ± 0.009) | 3 / 3                              |

**Response Table 3.1.** Per-seed chr1 SNV results of the three-reference comparison.

Precision, recall and F1 in benchmark regions after population filtering (our four resource VCFs), at validation-selected checkpoints and PASS thresholds (validation recall 0.9 before filtering). Each metric gives seeds 0 / 1 / 2 followed by mean ± sample standard deviation (SD; n − 1 denominator). The final column counts seeds in which the HPRC model has higher F1 than the GRCh38-only model of the same seed.

*Models and training.* Models are those of the three-reference comparison in R3.1, with chr2–22 of the same tumor used for training, chr1 held out, and base and mapping qualities standardised using the training statistics. The HPRC v2.1 models used for comparisons with other callers (R3.5) standardise base quality per dataset; their three-seed mean F1 differs from the corresponding v2.1 models in this table by at most 0.008.

*Exploratory, not used in any comparison:* COLO829T Illumina on HPRC v1.1 at the earlier epochs 6 / 5 / 6, examined after the chr1 results were seen (base quality standardised per dataset): precision 0.832 / 0.705 / 0.782, recall 0.723 / 0.705 / 0.698, F1 0.773 / 0.705 / 0.738 (0.739 ± 0.034). From epoch 7 (seed 1: epoch 6) these models call unlabelled candidates in chr1:121–126 and 143–144 Mb (1,075–1,639 false calls at the validation-selected epoch 8; 15–502 at epochs 6 / 5 / 6) without any change in validation F1.

*Reproducibility.* Three SNV models retrained with the same seed reproduced their chr1 calls exactly. Training is not bit-reproducible in general: of nine final INDEL models retrained on the same GPU type, three reproduced the original exactly and six differed at the level of floating-point rounding; at the same epoch, chr1 F1 differed by up to 0.041. Re-running prediction and evaluation of 30 checkpoints on a second cluster reproduced their chr1 F1 after population filtering within 0.002.

Response Table 3.2

| **tumor** | **platform** | **training set** | **precision**                         | **recall**                            | **F1**                                |
| --------- | ------------ | ---------------- | ------------------------------------- | ------------------------------------- | ------------------------------------- |
| HG008T    | Illumina     | 6 tumors         | 0.376 / 0.379 / 0.460 (0.405 ± 0.048) | 0.111 / 0.138 / 0.148 (0.132 ± 0.019) | 0.171 / 0.202 / 0.224 (0.199 ± 0.027) |
| HG008T    | Illumina     | 4 tumors         | 0.484 / 0.552 / 0.534 (0.524 ± 0.035) | 0.080 / 0.063 / 0.066 (0.070 ± 0.009) | 0.137 / 0.113 / 0.118 (0.123 ± 0.013) |
| HG008T    | PacBio HiFi  | 6 tumors         | 0.468 / 0.448 / 0.510 (0.475 ± 0.031) | 0.160 / 0.170 / 0.177 (0.169 ± 0.009) | 0.239 / 0.247 / 0.263 (0.249 ± 0.012) |
| HG008T    | PacBio HiFi  | 4 tumors         | 0.431 / 0.485 / 0.520 (0.479 ± 0.045) | 0.080 / 0.082 / 0.090 (0.084 ± 0.005) | 0.135 / 0.140 / 0.154 (0.143 ± 0.010) |
| HG008T    | ONT          | 6 tumors         | 0.192 / 0.212 / 0.153 (0.185 ± 0.030) | 0.150 / 0.153 / 0.162 (0.155 ± 0.006) | 0.168 / 0.178 / 0.157 (0.168 ± 0.010) |
| HG008T    | ONT          | 4 tumors         | 0.419 / 0.439 / 0.416 (0.424 ± 0.012) | 0.084 / 0.085 / 0.089 (0.086 ± 0.003) | 0.139 / 0.143 / 0.146 (0.143 ± 0.003) |
| COLO829T  | Illumina     | 6 tumors         | 0.262 / 0.176 / 0.354 (0.264 ± 0.089) | 0.130 / 0.145 / 0.130 (0.135 ± 0.009) | 0.173 / 0.159 / 0.190 (0.174 ± 0.015) |
| COLO829T  | Illumina     | 4 tumors         | 0.405 / 0.464 / 0.448 (0.439 ± 0.030) | 0.115 / 0.099 / 0.099 (0.104 ± 0.009) | 0.179 / 0.164 / 0.163 (0.168 ± 0.009) |
| COLO829T  | Fiber-seq    | 6 tumors         | 0.350 / 0.327 / 0.275 (0.317 ± 0.039) | 0.107 / 0.137 / 0.107 (0.117 ± 0.018) | 0.164 / 0.194 / 0.154 (0.170 ± 0.021) |
| COLO829T  | Fiber-seq    | 4 tumors         | 0.231 / 0.154 / 0.177 (0.187 ± 0.039) | 0.115 / 0.092 / 0.107 (0.104 ± 0.012) | 0.153 / 0.115 / 0.133 (0.134 ± 0.019) |
| COLO829T  | ONT          | 6 tumors         | 0.082 / 0.060 / 0.048 (0.064 ± 0.017) | 0.160 / 0.160 / 0.160 (0.160 ± 0.000) | 0.109 / 0.088 / 0.074 (0.090 ± 0.017) |
| COLO829T  | ONT          | 4 tumors         | 0.047 / 0.036 / 0.066 (0.050 ± 0.015) | 0.130 / 0.107 / 0.122 (0.120 ± 0.012) | 0.069 / 0.054 / 0.085 (0.070 ± 0.016) |

**Response Table 3.2.** Per-seed chr1 INDEL results of the final models.

Precision, recall and F1 in benchmark regions after gnomAD/CoLoRSdb filtering, at validation-selected checkpoints and PASS thresholds (validation recall 0.8 after filtering). Each metric gives seeds 0 / 1 / 2 followed by mean ± sample standard deviation (SD; n − 1 denominator).

*Models and training.* All models use HPRC v2.1. Six-tumor training includes HCC1395T, HCC1954T, H1437T, H2009T, HG008T and COLO829T; four-tumor training uses the first four and excludes HG008T and COLO829T. Chr1 is held out for every training tumor, with chr20 also held out for HCC1395T. Population filtering uses AF ≥ 0.01 in homopolymers ≥ 8 bp and AF ≥ 1e-4 elsewhere.

Response Table 3.3

| **chr1 SNVs, no population filter (P / R / F1)** | **Pansoma GRCh38-only** | **Pansoma HPRC v1.1** | **Pansoma HPRC v2.1**  | **ClairS-TO**         | **DeepSomatic-TO**    |
| ------------------------------------------------ | ----------------------- | --------------------- | ---------------------- | --------------------- | --------------------- |
| HG008T Illumina                                  | 0.005 / 0.838 / 0.010   | 0.106 / 0.848 / 0.189 | 0.194 / 0.864 / 0.316  | 0.003 / 0.902 / 0.006 | 0.029 / 0.935 / 0.056 |
| HG008T PacBio HiFi                               | 0.005 / 0.862 / 0.010   | 0.098 / 0.856 / 0.175 | 0.207 / 0.852 / 0.332  | 0.003 / 0.894 / 0.006 | 0.059 / 0.829 / 0.109 |
| HG008T ONT                                       | 0.004 / 0.867 / 0.009   | 0.094 / 0.855 / 0.170 | 0.190 / 0.878 / 0.312  | 0.003 / 0.884 / 0.006 | 0.047 / 0.910 / 0.090 |
| COLO829T Illumina                                | 0.044 / 0.593 / 0.081   | 0.297 / 0.722 / 0.421 | 0.580 / 0.728 / 0.644  | 0.010 / 0.796 / 0.020 | 0.097 / 0.858 / 0.173 |
| COLO829T Fiber-seq                               | 0.052 / 0.616 / 0.096   | 0.435 / 0.733 / 0.545 | 0.613 / 0.734 / 0.668  | 0.011 / 0.771 / 0.021 | 0.255 / 0.569 / 0.352 |
| COLO829T ONT                                     | 0.038 / 0.565 / 0.071   | 0.416 / 0.708 / 0.524 | 0.532 / 0.722 / 0.612  | 0.010 / 0.709 / 0.020 | 0.152 / 0.714 / 0.251 |
| HCC1395T Illumina                                | –                       | –                     | 0.706 / 0.830 / 0.763ᵃ | 0.014 / 0.884 / 0.028 | 0.139 / 0.944 / 0.242 |
| HCC1395T PacBio HiFi                             | –                       | –                     | 0.639 / 0.817 / 0.717ᵃ | 0.014 / 0.845 / 0.028 | 0.257 / 0.944 / 0.404 |
| HCC1395T ONT                                     | –                       | –                     | 0.560 / 0.704 / 0.624ᵃ | 0.013 / 0.752 / 0.025 | 0.217 / 0.781 / 0.340 |

**Response Table 3.3.** SNVs without any population filter on the three references, with the comparators (chr1, benchmark regions, PASS calls; precision / recall / F1). Pansoma before our population filter (it has none of its own): HG008T and COLO829T, the single-tumor models of each reference (trained on that tumor's chr2–22), mean of three seeds (GRCh38-only: training-set statistics, R3.1; HPRC v1.1 and v2.1: base quality standardised per dataset). ᵃ HCC1395T has tensors only on HPRC v2.1 and no single-tumor model: the six-tumor model (HCC1395T chr1 not in training, one run). ClairS-TO rerun without its panel-of-normals step; DeepSomatic-TO with its PON filter lifted (it still takes population allele frequencies as input). After filtering: Response Table 3.1 (single-tumor models).

Response Table 3.4

| **chr1 INDELs, no population filter (P / R / F1)** | **Pansoma GRCh38-only** | **Pansoma HPRC v1.1** | **Pansoma HPRC v2.1**  | **ClairS-TO**         | **DeepSomatic-TO**    |
| -------------------------------------------------- | ----------------------- | --------------------- | ---------------------- | --------------------- | --------------------- |
| HG008T Illumina                                    | 0.262 / 0.245 / 0.242   | 0.083 / 0.297 / 0.129 | 0.092 / 0.333 / 0.144  | 0.004 / 0.101 / 0.007 | 0.017 / 0.206 / 0.032 |
| HG008T PacBio HiFi                                 | 0.213 / 0.258 / 0.213   | 0.350 / 0.229 / 0.274 | 0.395 / 0.215 / 0.278  | 0.004 / 0.092 / 0.007 | 0.080 / 0.037 / 0.051 |
| HG008T ONT                                         | 0.019 / 0.361 / 0.035   | 0.025 / 0.233 / 0.044 | 0.011 / 0.300 / 0.020  | 0.004 / 0.087 / 0.007 | 0.043 / 0.019 / 0.026 |
| COLO829T Illumina                                  | 0.175 / 0.209 / 0.182   | 0.039 / 0.293 / 0.066 | 0.023 / 0.354 / 0.043  | 0.002 / 0.191 / 0.003 | 0.010 / 0.244 / 0.019 |
| COLO829T Fiber-seq                                 | 0.158 / 0.209 / 0.180   | 0.064 / 0.295 / 0.093 | 0.074 / 0.237 / 0.111  | 0.002 / 0.160 / 0.004 | 0.124 / 0.092 / 0.105 |
| COLO829T ONT                                       | 0.005 / 0.366 / 0.010   | 0.004 / 0.257 / 0.008 | 0.002 / 0.308 / 0.004  | 0.002 / 0.137 / 0.004 | 0.016 / 0.015 / 0.015 |
| HCC1395T Illumina                                  | –                       | –                     | 0.417 / 0.744 / 0.498ᵃ | 0.006 / 0.699 / 0.012 | 0.022 / 0.835 / 0.042 |
| HCC1395T PacBio HiFi                               | –                       | –                     | 0.515 / 0.699 / 0.589ᵃ | 0.005 / 0.526 / 0.011 | 0.271 / 0.647 / 0.382 |
| HCC1395T ONT                                       | –                       | –                     | 0.165 / 0.642 / 0.260ᵃ | 0.006 / 0.564 / 0.012 | 0.249 / 0.481 / 0.328 |

**Response Table 3.4.** INDELs without any population filter on the three references, with the comparators, as in Response Table 3.3. Pansoma: HG008T and COLO829T, single-tumor models with the final INDEL recipe on each reference, mean of three seeds; these models take gnomAD and CoLoRSdb allele frequencies as input, so, unlike the SNV models, their values without the filter are not free of population information. ᵃ HCC1395T has tensors only on HPRC v2.1: the final six-tumor INDEL models (HCC1395T chr1 and chr20 not in training), three seeds. After filtering: R3.1 (point 3, single-tumor models).

2\. Candidate reduction should be evaluated together with retention of bona fide somatic variants  
The reduction in candidate burden is one of the strongest mechanistic results. However, fewer candidates are beneficial only if germline variants and technical artifacts are preferentially removed while true somatic variants are retained.  
The authors currently emphasize overall candidate reduction, but the composition of removed candidates is not sufficiently resolved.  
I suggest comparing linear-reference and graph-based candidate generation for at least:  
1\) benchmark somatic truth variants,  
2\) matched-normal germline variants,  
3\) background/non-truth candidates.  
To mechanistic link graph representation with downstream accuracy, the authors should consider stratifying false-positive reductions based on whether the underlying germline allele is encoded within the HPRC pangenome graph. Such an analysis would directly connect the observation that over 79% of normal-sample SNVs and 73% of INDELs overlap with the AF-HPRC graph to the observed gains in precision, thereby demonstrating that graph-level allele representation is the primary driver of germline variant suppression.

**Response:** We thank the reviewer for these insights. Using the same reads and candidate-generation code, we compared somatic truth alleles, matched-normal germline alleles and background candidates between the GRCh38-graph and HPRC v2.1 graphs (Response Table 3.5).

**Candidate composition.** On HPRC v2.1, approximately 89–96% of somatic truth SNVs remained candidates, while the proportion of germline SNVs that became candidates fell from 82–88% on GRCh38-graph to 8–9%. Non-truth SNV candidates also decreased by 5–59% (Response Table 3.5). We additionally observed a substantial reduction in somatic INDEL candidates in both HG008T and COLO829T. To investigate this loss, we examined benchmark somatic truth variants absorbed by the graph and found that most were located in tandem repeats, homopolymers, and other short repetitive regions (Response Table 3.7). We further compared these absorbed somatic variants against existing population resources, including our four resource VCFs, and found substantial overlap (Response Table 3.6).

**Graph representation and precision.** To explore the relationship between allele, overlap and precision improvement, we examined Pansoma calls before population filtering on the GRCh38-graph. Using HPRC haplotype support as a proxy for graph representation, graph-represented germline alleles accounted for 71–86% of SNV false positives. With the same reads aligned to HPRC v2.1, only 0.2–0.6% of these alleles retained tensors, whereas 83–88% of graph-absent germline alleles still formed tensors (Response Table 3.8). This selective reduction before population filtering supports a direct link between graph representation and improved SNV precision. For graph-represented germline INDELs, tensor retention was 0.1–0.3% on Illumina and 15–22% on long reads (Response Table 3.9).

Response Table 3.5

| **class (HG008T; COLO829T)**                            | **HG008T GRCh38-graph** | **HG008T HPRC v2.1** | **COLO829T GRCh38-graph** | **COLO829T HPRC v2.1** |
| ------------------------------------------------------- | ----------------------- | -------------------- | ------------------------- | ---------------------- |
| somatic SNVs (8,543; 42,035)                            | 98.5 / 98.4 / 98.7      | 95.8 / 95.2 / 95.6   | 90.9 / 90.3 / 89.6        | 90.1 / 89.8 / 89.2     |
| somatic INDELs (7,933; 1,912)                           | 80.1 / 84.3 / 58.0      | 45.9 / 50.6 / 43.7   | 91.4 / 96.3 / 75.4        | 60.4 / 65.9 / 58.4     |
| germline SNVs (3,737,022; 3,729,123)                    | 82.4 / 82.8 / 83.4      | 7.5 / 7.6 / 8.0      | 86.5 / 87.2 / 88.3        | 8.0 / 8.4 / 8.9        |
| germline INDELs (964,490; 977,293)                      | 70.6 / 75.7 / 64.4      | 5.4 / 14.9 / 14.7    | 72.6 / 78.7 / 65.2        | 5.7 / 15.8 / 17.3      |
| background SNV candidates (label 0), change in number   | –                       | −12 / −59 / −39      | –                         | −28 / −35 / −5         |
| background INDEL candidates (label 0), change in number | –                       | −84 / −1 / +1        | –                         | −81 / +4 / +2          |

**Response Table 3.5.** Candidate retention (%) and background candidate count changes (%) on GRCh38-graph and HPRC v2.1 within the benchmark regions on chr1–22. Platform order: Illumina / PacBio HiFi / ONT (HG008T); Illumina / Fiber-seq / ONT (COLO829T). Parentheses give truth counts (HG008T; COLO829T). Background candidates match neither a somatic nor a germline truth allele (label 0).

Response Table 3.6

| **Truth INDELs (chr1–22)**                     | **HG008T**        | **COLO829T**    | **HCC1395T**   |
| ---------------------------------------------- | ----------------- | --------------- | -------------- |
| Absorbed by the graph (union across platforms) | 2,277             | 430             | 52             |
| Homopolymer ≥ 7 bp                             | 1,853 (81.4%)     | 241 (56.0%)     | 30 (57.7%)     |
| Short tandem repeat (unit 2–6 bp, ≥ 3 copies)  | 292 (12.8%)       | 83 (19.3%)      | 5 (9.6%)       |
| **Homopolymer + short tandem repeat, total**   | **2,145 (94.2%)** | **324 (75.3%)** | **35 (67.3%)** |
| Other repeats                                  | 6 (0.3%)          | 20 (4.7%)       | 1 (1.9%)       |
| Substitution within repetitive sequence        | 104 (4.6%)        | 16 (3.7%)       | 2 (3.8%)       |
| Others                                         | 21 (0.9%)         | 61 (14.2%)      | 12 (23.1%)     |
| Non-repetitive sequence                        | 1 (0.0%)          | 9 (2.1%)        | 2 (3.8%)       |

**Response Table 3.6.** Sequence context of somatic truth INDELs absorbed by the HPRC v2.1 graph (chr1–22). An INDEL is absorbed when, on at least one platform, its ALT reads align to a graph path without any edit, leaving no candidate. First row: absorbed INDELs and their percentage of all truth INDELs; other rows: percentages of the absorbed INDELs. Classes describe the somatic change relative to the matched-normal haplotype. Other repeats: shorter homopolymers (4–6 bp), longer-unit tandem repeats, or INDELs that do not change a repeat by whole units. Substitution within repetitive sequence: a base substitution on a germline length difference. Others: allele already on the normal haplotype, complex change, or unresolved normal haplotype.

Response Table 3.7 Somatic candidate loss on the HPRC graph

<table>
<thead>
<tr class="header">
<th><strong>% of all truth alleles (ranges over three platforms)</strong></th>
<th><strong>HG008T<br />
INDELs</strong></th>
<th><strong>COLO829T<br />
INDELs</strong></th>
<th><strong>HCC1395T<br />
INDELs</strong></th>
<th><strong>HG008T<br />
SNVs</strong></th>
<th><strong>COLO829T<br />
SNVs</strong></th>
<th><strong>HCC1395T<br />
SNVs</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>truth alleles, chr1–22</td>
<td>8,496</td>
<td>1,912</td>
<td>1,625</td>
<td>8,690</td>
<td>42,035</td>
<td>39,447</td>
</tr>
<tr class="even">
<td>somatic truth alleles tagged by population filtering</td>
<td>62.7%</td>
<td>52.6%</td>
<td>7.9%</td>
<td>10.5%</td>
<td>3.4%</td>
<td>2.8%</td>
</tr>
<tr class="odd">
<td>candidate on the GRCh38-only graph but not on HPRC v1.1</td>
<td>31.1–44.8%</td>
<td>25.8–35.1%</td>
<td>4.1–5.5%</td>
<td>3.4–3.9%</td>
<td>0.5–0.9%</td>
<td>0.6–0.8%</td>
</tr>
<tr class="even">
<td><strong>graph-specific loss left after filtering, HPRC v1.1</strong></td>
<td>3.8–7.2%</td>
<td>5.0–6.8%</td>
<td>1.3–1.9%</td>
<td>0.3–0.4%</td>
<td>0.2–0.6%</td>
<td>0.1–0.2%</td>
</tr>
<tr class="odd">
<td>candidate on the GRCh38-only graph but not on HPRC v2.1</td>
<td>29.4–42.1%</td>
<td>23.9–32.0%</td>
<td>4.0–4.9%</td>
<td>3.0–3.6%</td>
<td>0.5–0.9%</td>
<td>0.8–0.9%</td>
</tr>
<tr class="even">
<td><strong>graph-specific loss left after filtering, HPRC v2.1</strong></td>
<td>3.5–7.0%</td>
<td>5.4–6.9%</td>
<td>1.1–1.4%</td>
<td>0.2–0.4%</td>
<td>0.2–0.7%</td>
<td>0.1–0.3%</td>
</tr>
</tbody>
</table>

**Response Table 3.7.** Somatic SNVs and INDELs lost during graph-based candidate generation and their overlap with population filtering. Population filtering follows our resource VCFs: four for SNVs; gnomAD/CoLoRSdb only for INDELs: SNVs use gnomAD, dbSNP 138 (non-somatic records), the GATK 1000 Genomes PoN and CoLoRSdb; INDELs use gnomAD and CoLoRSdb only. A filter tag indicates that a somatic truth allele meets the filtering criteria, not that it is germline in the patient. Candidate presence requires an exact match to any candidate allele, including a non-representative allele. Percentages use the complete chr1–22 truth counts. Ranges span three platforms.

Response Table 3.8

| **dataset**        | **GRCh38-graph FP before population filtering** | **HPRC v2.1 Graph ∩ matched normal germline (share; candidates retained in V2.1)** | **non HPRC v2.1 Graph but matched normal germline (share; candidates retained in V2.1)** | **non-germline FP (% of GRCh38-graph FP)** | **total HPRC v2.1 FP before population filtering** |
| ------------------ | ----------------------------------------------- | ---------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------- | ------------------------------------------ | -------------------------------------------------- |
| HG008T Illumina    | 121,139                                         | 84%; 0.3%                                                                          | 13%; 85%                                                                                 | 3%                                         | 2,870                                              |
| HG008T PacBio HiFi | 130,720                                         | 85%; 0.4%                                                                          | 14%; 83%                                                                                 | 1%                                         | 2,590                                              |
| HG008T ONT         | 147,531                                         | 86%; 0.5%                                                                          | 13%; 85%                                                                                 | 2%                                         | 3,038                                              |
| COLO829T Illumina  | 36,134                                          | 78%; 0.2%                                                                          | 9%; 87%                                                                                  | 13%                                        | 986                                                |
| COLO829T Fiber-seq | 37,622                                          | 84%; 0.6%                                                                          | 10%; 88%                                                                                 | 6%                                         | 969                                                |
| COLO829T ONT       | 37,199                                          | 71%; 0.6%                                                                          | 10%; 87%                                                                                 | 18%                                        | 1,354                                              |

**Response Table 3.8.** SNV false positives before population filtering on GRCh38-graph and HPRC v2.1 (chr1, benchmark regions). The two germline columns group matched-normal germline alleles by whether they are represented in HPRC v2.1. Each cell gives the share of GRCh38-graph FP (%) and the percentage of those same alleles retained as candidates (forming tensors) on HPRC v2.1. Graph representation is inferred from carriage in at least 46 of 463 HPRC haplotypes. Non-germline gives the share of GRCh38-graph FP without a normal germline match; the first and final FP columns give the totals for each reference. Graph-represented germline alleles account for 71–86% of GRCh38-graph FP, but only 0.2–0.6% retain candidates, compared with 83–88% of graph-absent germline alleles.

Response Table 3.9

| **dataset**        | **GRCh38-graph FP before population filtering** | **HPRC v2.1 Graph ∩ matched normal germline (share; candidates retained in V2.1)** | **non HPRC v2.1 Graph but matched normal germline (share; candidates retained in V2.1)** | **non-germline FP (% of GRCh38-graph FP)** | **total HPRC v2.1 FP before population filtering** |
| ------------------ | ----------------------------------------------- | ---------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------- | ------------------------------------------ | -------------------------------------------------- |
| HG008T Illumina    | 39,930                                          | 79%; 0.3%                                                                          | 16%; 44%                                                                                 | 5%                                         | 2,867                                              |
| HG008T PacBio HiFi | 36,354                                          | 71%; 20%                                                                           | 18%; 37%                                                                                 | 11%                                        | 10,035                                             |
| HG008T ONT         | 65,734                                          | 46%; 15%                                                                           | 11%; 39%                                                                                 | 43%                                        | 41,068                                             |
| COLO829T Illumina  | 22,260                                          | 70%; 0.1%                                                                          | 14%; 36%                                                                                 | 16%                                        | 2,576                                              |
| COLO829T Fiber-seq | 37,523                                          | 71%; 18%                                                                           | 16%; 35%                                                                                 | 13%                                        | 11,714                                             |
| COLO829T ONT       | 70,514                                          | 35%; 22%                                                                           | 8%; 42%                                                                                  | 57%                                        | 57,713                                             |

**Response Table 3.9.** INDEL false positives before population filtering on GRCh38-graph and HPRC v2.1 (chr1, benchmark regions). In each germline column, the two values give the share of GRCh38-graph FP (%) and the percentage of those same alleles retained as candidates on HPRC v2.1. Germline and graph-representation groups follow Response Table 3.8. Non-germline gives the share without a normal germline match, and the first and final FP columns give the totals for each reference. Graph-represented germline INDEL candidate retention is 0.1–0.3% on Illumina and 15–22% on long reads; non-germline calls account for 43–57% of GRCh38-graph ONT FP.

3\. Variants residing in graph coordinates require more rigorous validation  
Among Pansoma's most consequential capabilities is its capacity to retain variants that cannot be unambiguously encoded on GRCh38. The authors report that a subset of final predictions exists exclusively in graph coordinates and offer mutational-signature analyses together with representative read-level examples as evidence of biological plausibility. This is a compelling and potentially important result; however, the supporting evidence remains largely indirect.  
The authors are encouraged to perform orthogonal validation on a representative sample of graph-exclusive variants by integrating complementary data sources, including matched-normal alignment evidence, phased haplotype assemblies, or cross-platform validation datasets. Demonstrating the authenticity of these unique calls through independent confirmation would significantly strengthen the manuscript's thesis, proving that pangenome reference architectures reliably extend variant detection into genomic regions inaccessible to conventional linear alignment workflows.

Response. We thank the reviewer for suggesting these additional validation analyses. We examined off-reference SNV PASS calls from HG008T, COLO829T and HCC1395T across three sequencing platforms using matched-normal assemblies, tumor assemblies and matched-normal reads (HPRC v1.1 and v2.1; HCC1395T: v2.1 only; Response Table 3.10).

1.  Phased haplotype assemblies. Each call's ALT context was aligned to both matched-normal haplotypes and a tumor assembly.

2.  Matched-normal reads. We used HG008-N PacBio HiFi, COLO829BL ONT and HCC1395-BL PacBio HiFi reads, with tumor reads of the same platforms.

In the default output, 34 calls from the HPRC v2.1 model in chr1 were supported by tumor-only assemblies, consistent with potential somatic origin. However, the majority of off-reference calls lacked sufficient supporting evidence. (Response Table 3.10).

We hypothesize that this limited off-reference performance partly reflects the training distribution: all somatic training positives lie on the GRCh38 path. The revised Pansoma pipeline also incorporates graph-derived features, including the haplotype-count channel, which may reinforce a learned preference for reference-path predictions and limit off-reference prediction. We have therefore narrowed our claims about off-reference detection and treat these predictions as exploratory candidates requiring further validation.

Response Table 3.10

<table>
<thead>
<tr class="header">
<th><strong>Dataset</strong><br />
<strong>(chr1 SNV)</strong></th>
<th><strong>Graph-</strong><br />
<strong>coordinate</strong><br />
<strong>PASS calls</strong></th>
<th><strong>In a normal</strong><br />
<strong>assembly</strong></th>
<th><strong>Tumor-only</strong></th>
<th><strong>Neither</strong></th>
<th><strong>ALT in ≥2</strong><br />
<strong>normal</strong><br />
<strong>reads</strong></th>
<th><strong>No normal</strong><br />
<strong>ALT read,</strong><br />
<strong>ALT in</strong><br />
<strong>tumor reads</strong></th>
<th><strong>Tumor-only</strong><br />
<strong>and no</strong><br />
<strong>normal</strong><br />
<strong>ALT read</strong></th>
<th><strong>Also PASS</strong><br />
<strong>on both</strong><br />
<strong>other</strong><br />
<strong>platforms</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>HG008T Illumina</td>
<td>29</td>
<td>25 (86.2%)</td>
<td>3 (10.3%)</td>
<td>1 (3.4%)</td>
<td>27 (93.1%)</td>
<td>2 (6.9%)</td>
<td>2</td>
<td>3</td>
</tr>
<tr class="even">
<td>HG008T PacBio HiFi</td>
<td>6</td>
<td>5 (83.3%)</td>
<td>1 (16.7%)</td>
<td>0 (0.0%)</td>
<td>5 (83.3%)</td>
<td>1 (16.7%)</td>
<td>1</td>
<td>3</td>
</tr>
<tr class="odd">
<td>HG008T ONT</td>
<td>162</td>
<td>158 (97.5%)</td>
<td>3 (1.9%)</td>
<td>1 (0.6%)</td>
<td>66 (40.7%)</td>
<td>7 (4.3%)</td>
<td>2</td>
<td>3</td>
</tr>
<tr class="even">
<td>COLO829T Illumina</td>
<td>18</td>
<td>14 (77.8%)</td>
<td>4 (22.2%)</td>
<td>0 (0.0%)</td>
<td>10 (55.6%)</td>
<td>3 (16.7%)</td>
<td>3</td>
<td>1</td>
</tr>
<tr class="odd">
<td>COLO829T Fiber-seq</td>
<td>78</td>
<td>72 (92.3%)</td>
<td>4 (5.1%)</td>
<td>2 (2.6%)</td>
<td>59 (75.6%)</td>
<td>2 (2.6%)</td>
<td>1</td>
<td>1</td>
</tr>
<tr class="even">
<td>COLO829T ONT</td>
<td>98</td>
<td>90 (91.8%)</td>
<td>4 (4.1%)</td>
<td>4 (4.1%)</td>
<td>46 (46.9%)</td>
<td>6 (6.1%)</td>
<td>3</td>
<td>1</td>
</tr>
<tr class="odd">
<td>HCC1395T Illumina</td>
<td>13</td>
<td>9 (69.2%)</td>
<td>3 (23.1%)</td>
<td>1 (7.7%)</td>
<td>7 (53.8%)</td>
<td>2 (15.4%)</td>
<td>2</td>
<td>4</td>
</tr>
<tr class="even">
<td>HCC1395T PacBio HiFi</td>
<td>56</td>
<td>50 (89.3%)</td>
<td>6 (10.7%)</td>
<td>0 (0.0%)</td>
<td>43 (76.8%)</td>
<td>5 (8.9%)</td>
<td>5</td>
<td>4</td>
</tr>
<tr class="odd">
<td>HCC1395T ONT</td>
<td>81</td>
<td>73 (90.1%)</td>
<td>6 (7.4%)</td>
<td>2 (2.5%)</td>
<td>40 (49.4%)</td>
<td>6 (7.4%)</td>
<td>5</td>
<td>4</td>
</tr>
</tbody>
</table>

**Response Table 3.10.** Validation of graph-coordinate SNV PASS calls on HPRC v2.1 (chr1, default output without off-reference rescue). Every call was evaluated. Assembly classes are mutually exclusive: in a normal assembly (ALT context matches a matched-normal haplotype), tumor-only (matches the tumor assembly only) and neither (matches neither assembly or cannot be aligned uniquely). “ALT in ≥ 2 normal reads” requires at least two normal ALT reads; “no normal ALT read, ALT in tumor reads” requires no normal ALT read, at least five normal REF reads and at least two tumor ALT reads. The next column requires both tumor-only assembly support and this read criterion. The final column counts calls with the same node, offset and allele also PASS on both other platforms of the same tumor and graph. Assemblies: HG008-N v6.2 and HG008-T v3.2; the COLO829BL donor-specific assembly and a hifiasm assembly of COLO829T HiFi/ONT reads; hifiasm assemblies of HCC1395-BL and HCC1395 HiFi reads and HCC1395BL v1.0. COLO829T and HCC1395T tumor assemblies reuse the analysed long reads and provide consistency rather than independent validation. Normal/tumor read platforms are PacBio HiFi for HG008 and HCC1395, and ONT for COLO829. HCC1395T uses six-tumor models. Of 541 calls, 34 (6.3%) have tumor-only assembly support.

1.  
2.  
3.  
4\. The reported benchmark metrics do not evaluate the complete Pansoma output  
Conventional TP/FP/FN evaluation is performed only on Pansoma calls that can be projected back to GRCh38. Calls that remain exclusively in graph coordinates are analyzed separately. Thus, the reported precision, recall and F1 values represent the GRCh38-projectable subset of Pansoma calls rather than the complete output of the method.  
This distinction is important because the graph-coordinate fraction is platform dependent and can be significant, particularly for ONT data. Supplementary Table 6 reports graph-coordinate fractions exceeding 10% and, in some settings, 20%. In addition, the statement in the Discussion that only “about 2%” of variants remain on alternative graph paths appears inconsistent with the higher fractions reported in Supplementary Table 6 and should be re-evaluated.  
An additional consideration is whether calls classified as graph-coordinate are truly dependent on graph representation or are simply not projectable to GRCh38. The authors note that the projection framework can be extended to CHM13. It would therefore be informative to determine what fraction of GRCh38-unprojectable calls can be represented on the T2T-CHM13 reference, and what fraction remains genuinely graph-specific even when a more complete linear reference is considered. Such an analysis would help distinguish limitations of GRCh38 from advantages that are intrinsic to the pangenome graph representation.

Response. We thank the reviewer for these insights and suggestions. We projected GRCh38-unprojectable SNV PASS calls to the embedded T2T-CHM13 path. In the default output, 396 of 541 v2.1 calls (73.2%) have a CHM13 allele (Response Table 3.11). Pansoma can also export these calls in CHM13 coordinates. We hypothesize that this high CHM13 representation partly reflects the training bias discussed in R3.3: all somatic training positives lie on the GRCh38 path, while the haplotype-count channel records how many haplotype paths traverse each node. This combination may favor positive predictions on nodes traversed by more paths.

We also evaluated off-reference calls using matched-normal assemblies and reads; the results are presented in R3.3 (Response Table 3.10). The inaccurate “about 2%” statement has been removed: the graph-coordinate fraction in the default output varies from 0.5–10.1% for HG008T/COLO829T and 0.4–2.2% for HCC1395T.

However, we respectfully differ on combining all Pansoma calls into a single TP/FP/FN evaluation under the available benchmarks. The benchmark datasets define truth alleles and confident evaluation regions in GRCh38 coordinates. Calls without an unambiguous GRCh38 projection cannot be assessed directly against these benchmarks; insufficient truth evidence prevents reliable assignment of TP or FP status. These datasets also lack a somatic truth set defined on the corresponding HPRC graph release, so we cannot directly quantify precision and recall across the complete graph-coordinate output. Additionally, the mapping between matched-normal assemblies and the HPRC graph contains ambiguities and gaps, preventing us from reliably assigning FN status and establishing a complete graph-wide somatic benchmark.

We therefore report precision, recall and F1 for the GRCh38-evaluable calls within benchmark regions and assess off-reference calls separately. These metrics describe the evaluable subset and are not presented as the accuracy of Pansoma's complete output.

Response Table 3.11

| **Tumor**                | **Platform** | **Calls**       | **CHM13 allele** |
| ------------------------ | ------------ | --------------- | ---------------- |
| HG008T                   | Illumina     | 29              | 17 (58.6%)       |
| HG008T                   | PacBio HiFi  | 6               | 6 (100%)         |
| HG008T                   | ONT          | 162             | 147 (90.7%)      |
| COLO829T                 | Illumina     | 18              | 13 (72.2%)       |
| COLO829T                 | Fiber-seq    | 78              | 66 (84.6%)       |
| COLO829T                 | ONT          | 98              | 31 (31.6%)       |
| HCC1395T                 | Illumina     | 13              | 11 (84.6%)       |
| HCC1395T                 | PacBio HiFi  | 56              | 51 (91.1%)       |
| HCC1395T                 | ONT          | 81              | 54 (66.7%)       |
| **Total**                | **541**      | **396 (73.2%)** |                  |
| Total excluding HCC1395T | 391          | 280 (71.6%)     |                  |

**Response Table 3.11.** CHM13 representation of graph-coordinate SNV PASS calls on HPRC v2.1 (chr1, default output; nine runs). HCC1395T uses six-tumor models. Overall, 396 of 541 calls (73.2%) have a CHM13 allele, including 116 of 150 HCC1395T calls (77.3%). Excluding HCC1395T, the total is 280 of 391 calls (71.6%). A CHM13 coordinate establishes representability, not somatic status.

5\. The Figure 4 mutational-signature analysis should incorporate the truth set directly  
In Figure 4a, the overlap analysis includes the benchmark truth set, yet the mutational-signature stacked bars presented immediately below encompass only Pansoma, ClairS-TO, and DeepSomatic-TO, with the ground-truth signature entering only indirectly as the reference for cosine-similarity calculation. I recommend displaying the spectra side by side, ground truth alongside Pansoma, DeepSomatic-TO, and ClairS-TO, so that the specific signatures underlying inter-caller agreement or disagreement can be readily identified. This is particularly relevant for HG008T, whose mutational profile is heterogeneous and for which the reported cosine similarities are comparatively modest.   
To strengthen the claim that graph-coordinate variants reflect true biological events through mutational-signature concordance, the authors should provide the underlying 96-channel single base substitution (SBS96) profiles in the Supplementary Information.

**Response:** We thank the reviewer for these helpful suggestions. We have revised Figure 4a to show the benchmark truth spectrum alongside those of Pansoma, DeepSomatic-TO and ClairS-TO, with the matched-normal germline spectrum included for reference (Response Figure 3.2). We have also added the SBS96 profiles of graph-coordinate SNV calls and their summary statistics to the Supplementary Information (Response Figure 3.3).

Response Figure 3.2

![](./media/image8.png)

**Response Figure 3.2.** Revised Figure 4a showing SBS96 spectra (percentage of SNVs in each channel) for the chr1 somatic truth set, matched-normal germline SNVs in the same benchmark regions (HG008-N or COLO829BL, Dipcall), and PASS SNV calls after population filtering from Pansoma (HPRC v2.1, d46), DeepSomatic-TO (v1.9.0), and ClairS-TO (v0.4.1/v0.4.2), shown separately for each sequencing platform.

Response Figure 3.3

![](./media/image9.png)

**Response Figure 3.3.** SBS96 spectra of Pansoma graph-coordinate SNV PASS calls on chr1, shown separately for each tumor and sequencing platform. Because these calls have no GRCh38 coordinate, the standard GRCh38-coordinate population filter and benchmark-region restriction do not apply. Trinucleotide context was derived directly from the graph path using the reference context surrounding each call; 7 of 6,025 calls without a usable context were excluded. For COLO829T ONT, a random subset of 3,000 of 4,569 calls was plotted.

6\. More perfect alignment does not necessarily imply more correct alignment  
The manuscript interprets increased perfect alignment and reduced mismatches/soft clipping as evidence of improved alignment to the pangenome. This interpretation is plausible.  
However, because the pangenome contains more alternative sequence paths than GRCh38, an increase in perfect alignments does not by itself prove that reads are mapped more correctly.  
To isolate the true mechanism behind the observed performance improvements, a control comparison evaluating alignment quality across GRCh38, T2T-CHM13, and the full HPRC graph is strongly recommended. Such an analysis is crucial to delineate improvements stemming simply from reference completeness (i.e., linear T2T assembly) versus those uniquely conferred by graph-based multi-haplotype modeling.  
If feasible, the authors could provide additional evidence using:  
MAPQ distributions,  
phased haplotype concordance,  
assembly-based truth,  
repetitive or segmentally duplicated regions,  
other difficult genomic regions.  
This would help distinguish “more perfectly aligned” from “more accurately aligned.”

**Response:** We agree that perfect alignment alone does not establish correct read placement. We revised the claims and report absolute alignment counts and proportions in Response Figure 3.4 and Response Table 3.12. We address the suggested evidence directly:

1.  **MAPQ distributions. **We examined MAPQ distributions and observed misplacement for simulated HG008N reads (Response Table 3.14). Response Table 3.13 reports high-confidence misplacements (MAPQ ≥ 30) for both samples and all five references.

2.  **Phased haplotype concordance. **We used phased diploid assemblies to establish simulated-read origins but did not perform a separate read-level haplotype-concordance analysis due to a lack of relevant resources.

3.  **Assembly-based truth.** Reads were simulated from HG008N v6.2 and COLO829BL on chr1 and chr20 and aligned to linear GRCh38, GRCh38-graph, CHM13-graph and both HPRC graphs. Response Table 3.13 compares reads with one-to-one truth positions on both GRCh38 and CHM13.

4.  **Repetitive and segmentally duplicated regions.** We stratified variant calling performance by low-mappability/segmental-duplication regions, tandem repeats and homopolymers (Response Table 3.15).

5.  **Other difficult regions**. SMaHT difficult/extreme-region analyses cover calling performance, real-read placement agreement and MAPQ changes (Response Table 3.15). The excess simulated ONT-like misplacement on HPRC v2.1 is concentrated in extreme regions.

Together, these experiments support improved read placement with vg Giraffe on HPRC graphs for Illumina and HiFi-like reads within the evaluable regions. We also observed increased misplacement of ONT-like reads on HPRC v2.1. The HPRC2 study reports greater tandem-repeat allele diversity, while the Giraffe study highlights how long-read length and sequencing errors increase the graph-alignment search space \[[HPRC2, 2026](https://doi.org/10.64898/2026.07.21.739710); [Chang et al., 2025](https://doi.org/10.1101/2025.09.29.678807)\]. We hypothesize that the expanded repeat-allele diversity in v2.1 increases alignment ambiguity for noisy reads and contributes to the excess misplacements concentrated in SMaHT extreme regions.

Response Figure 3.4 Reads on chr1–22 on four references

![Supplementary Figure R2.7: Reads on chr1–22 on four references](./media/image7.png)

**Response Figure 3.4.** Reads aligned to chr1–22 for each dataset and reference: the share perfectly aligned (no substitution, indel or soft clipping) and the shares with at least one substitution, indel or soft-clipped end. Rows: Illumina; PacBio (HG008T HiFi, COLO829T Fiber-seq); ONT. Note the different y-axis ranges. Proportions are computed among the reads placed on chr1–22 on each reference. These totals differ by up to 7% between references, because each reference leaves a different number of reads unmapped, on other contigs or only on non-reference nodes; all totals are listed in Response Table 3.12. The three graphs are aligned with vg Giraffe; linear GRCh38 uses the data providers' alignments.

Response Table 3.12

| **dataset**               | **reference**   | **reads on chr1–22** | **perfectly aligned**  | **with substitution** | **with indel**      | **with soft clipping** |
| ------------------------- | --------------- | -------------------- | ---------------------- | --------------------- | ------------------- | ---------------------- |
| HG008T Illumina           | GRCh38 (linear) | 2,739,599,264        | 1,869,903,341 (68.25%) | 810,918,802 (29.60%)  | 78,905,626 (2.88%)  | 96,334,574 (3.52%)     |
| HG008T Illumina           | GRCh38 graph    | 2,710,092,615        | 1,861,828,942 (68.70%) | 800,943,101 (29.55%)  | 79,136,287 (2.92%)  | 74,218,393 (2.74%)     |
| HG008T Illumina           | HPRC v1.1 graph | 2,704,033,209        | 2,121,622,695 (78.46%) | 560,922,772 (20.74%)  | 21,774,798 (0.81%)  | 66,387,175 (2.46%)     |
| HG008T Illumina           | HPRC v2.1 graph | 2,695,745,607        | 2,121,713,347 (78.71%) | 552,985,555 (20.51%)  | 19,851,148 (0.74%)  | 63,829,495 (2.37%)     |
| COLO829T Illumina         | GRCh38 (linear) | 2,296,226,777        | 1,604,645,813 (69.88%) | 627,386,769 (27.32%)  | 72,684,106 (3.17%)  | 95,094,133 (4.14%)     |
| COLO829T Illumina         | GRCh38 graph    | 2,283,443,818        | 1,599,146,760 (70.03%) | 629,366,832 (27.56%)  | 71,551,146 (3.13%)  | 78,056,743 (3.42%)     |
| COLO829T Illumina         | HPRC v1.1 graph | 2,285,859,727        | 1,823,482,724 (79.77%) | 423,169,962 (18.51%)  | 23,805,599 (1.04%)  | 78,789,217 (3.45%)     |
| COLO829T Illumina         | HPRC v2.1 graph | 2,265,562,422        | 1,822,701,086 (80.45%) | 410,679,337 (18.13%)  | 20,250,428 (0.89%)  | 66,561,782 (2.94%)     |
| HG008T PacBio HiFi        | GRCh38 (linear) | 12,931,641           | 18,818 (0.15%)         | 12,264,874 (94.84%)   | 12,862,468 (99.47%) | 961,155 (7.43%)        |
| HG008T PacBio HiFi        | GRCh38 graph    | 12,675,374           | 18,819 (0.15%)         | 12,020,862 (94.84%)   | 12,606,854 (99.46%) | 314,848 (2.48%)        |
| HG008T PacBio HiFi        | HPRC v1.1 graph | 12,689,479           | 179,824 (1.42%)        | 8,558,774 (67.45%)    | 12,377,219 (97.54%) | 228,857 (1.80%)        |
| HG008T PacBio HiFi        | HPRC v2.1 graph | 12,977,109           | 184,297 (1.42%)        | 8,605,339 (66.31%)    | 12,666,519 (97.61%) | 422,918 (3.26%)        |
| COLO829T PacBio Fiber-seq | GRCh38 (linear) | 25,308,248           | 44,846 (0.18%)         | 24,325,270 (96.12%)   | 25,138,324 (99.33%) | 2,529,764 (10.00%)     |
| COLO829T PacBio Fiber-seq | GRCh38 graph    | 25,026,903           | 44,873 (0.18%)         | 24,011,605 (95.94%)   | 24,858,107 (99.33%) | 936,884 (3.74%)        |
| COLO829T PacBio Fiber-seq | HPRC v1.1 graph | 25,070,732           | 386,300 (1.54%)        | 18,732,761 (74.72%)   | 24,353,478 (97.14%) | 717,920 (2.86%)        |
| COLO829T PacBio Fiber-seq | HPRC v2.1 graph | 26,066,512           | 397,244 (1.52%)        | 19,343,766 (74.21%)   | 25,350,225 (97.25%) | 1,538,959 (5.90%)      |
| HG008T ONT                | GRCh38 (linear) | 3,208,600            | 43,815 (1.37%)         | 3,009,593 (93.80%)    | 3,084,419 (96.13%)  | 1,215,595 (37.89%)     |
| HG008T ONT                | GRCh38 graph    | 3,100,754            | 43,748 (1.41%)         | 2,869,956 (92.56%)    | 2,950,998 (95.17%)  | 930,380 (30.00%)       |
| HG008T ONT                | HPRC v1.1 graph | 3,099,600            | 85,847 (2.77%)         | 2,691,454 (86.83%)    | 2,929,825 (94.52%)  | 913,057 (29.46%)       |
| HG008T ONT                | HPRC v2.1 graph | 3,087,083            | 84,623 (2.74%)         | 2,704,935 (87.62%)    | 2,944,975 (95.40%)  | 899,235 (29.13%)       |
| COLO829T ONT              | GRCh38 (linear) | 12,636,129           | 723 (0.01%)            | 11,794,165 (93.34%)   | 11,928,948 (94.40%) | 12,590,223 (99.64%)    |
| COLO829T ONT              | GRCh38 graph    | 11,839,534           | 778 (0.01%)            | 10,956,152 (92.54%)   | 11,047,335 (93.31%) | 11,723,760 (99.02%)    |
| COLO829T ONT              | HPRC v1.1 graph | 11,877,368           | 1,121 (0.01%)          | 10,662,112 (89.77%)   | 11,030,919 (92.87%) | 11,760,246 (99.01%)    |
| COLO829T ONT              | HPRC v2.1 graph | 11,789,372           | 1,098 (0.01%)          | 10,645,364 (90.30%)   | 11,021,767 (93.49%) | 11,674,637 (99.03%)    |

**Response Table 3.12.** Counts and proportions of reads aligned to chr1–22 on four references.

*Methods.*

1\. *Alignments.* Each dataset's tumor reads were aligned with vg Giraffe to the GRCh38-only graph and to the HPRC v1.1 (d9) and v2.1 (d46) graphs (versions as in R3.m1). Linear GRCh38 is the data providers' alignment: for HG008T (GIAB) BWA-MEM 0.7.17, pbmm2 1.10.0 and minimap2 2.26; for COLO829T (SMaHT) BWA-MEM 0.7.17 in Sentieon 202308.01 (with duplicate marking, indel realignment and base recalibration), pbmm2 1.13.0 and minimap2 2.26.

2\. *Reads counted.* Every read of the dataset, each mate of a pair separately, whose alignment is on chr1–22: on linear GRCh38 the primary alignment (secondary and supplementary alignments are not counted); on a graph, the GRCh38 chromosome holding most of the read's reference-path bases. Unmapped reads, reads on other contigs and reads aligned only to non-reference nodes are not counted.

3\. *Measures.* Perfectly aligned: aligned end to end with no substitution, indel or soft clipping. Substitution: at least one mismatched base (linear: NM minus inserted and deleted bases; graph: mismatches against the sequence of the read's graph path, so a match to a non-reference allele is not a substitution). Indel: at least one inserted or deleted base. Soft clipping: read bases at either end not aligned to the reference (linear: CIGAR S in the primary alignment; hard clipping occurs only in supplementary alignments, which are not counted; graph: read bases at an alignment end not aligned to the path).

4\. On linear GRCh38 a long read split by the aligner keeps its remainder as a supplementary alignment, so its primary alignment is soft-clipped; Giraffe writes one alignment per read.

Response Table 3.13

|                                           | **GRCh38 (linear ref)** | **GRCh38-graph** | **CHM13-graph** | **HPRC v1.1** | **HPRC v2.1** |
| ----------------------------------------- | ----------------------- | ---------------- | --------------- | ------------- | ------------- |
| **HG008-N**                               |                         |                  |                 |               |               |
| Illumina, perfectly aligned               | 63.7%                   | 63.7%            | 64.8%           | 73.0%         | 73.1%         |
| Illumina, misplaced (% of placed)         | 1.23%                   | 1.19%            | 1.12%           | 0.95%         | 0.90%         |
| Illumina, misplaced at MAPQ ≥ 30 (reads)  | 5,555                   | 4,240            | 4,574           | 823           | 567           |
| HiFi-like, misplaced                      | 0.51%                   | 0.47%            | 0.41%           | 0.31%         | 0.36%         |
| HiFi-like, misplaced at MAPQ ≥ 30 (reads) | 57                      | 11               | 13              | 3             | 1             |
| ONT-like, misplaced                       | 0.52%                   | 0.75%            | 0.68%           | 0.60%         | 1.45%         |
| ONT-like, misplaced at MAPQ ≥ 30 (reads)  | 45                      | 164              | 146             | 127           | 294           |
| **COLO829BL**                             |                         |                  |                 |               |               |
| Illumina, perfectly aligned               | 63.8%                   | 63.8%            | 64.8%           | 73.0%         | 73.1%         |
| Illumina, misplaced (% of placed)         | 1.13%                   | 1.09%            | 1.06%           | 0.88%         | 0.85%         |
| Illumina, misplaced at MAPQ ≥ 30 (reads)  | 5,491                   | 4,386            | 5,625           | 1,368         | 1,216         |
| HiFi-like, misplaced                      | 0.40%                   | 0.38%            | 0.42%           | 0.28%         | 0.34%         |
| HiFi-like, misplaced at MAPQ ≥ 30 (reads) | 31                      | 11               | 31              | 4             | 4             |
| ONT-like, misplaced                       | 0.42%                   | 0.64%            | 0.66%           | 0.56%         | 1.50%         |
| ONT-like, misplaced at MAPQ ≥ 30 (reads)  | 36                      | 153              | 145             | 129           | 315           |

**Response Table 3.13.** Placement accuracy of reads simulated from chr1 and chr20 of the phased diploid assemblies HG008N v6.2 and COLO829BL (assembled with Verkko v2.1). Both samples used identical simulation settings, generating 2.5 million Illumina read pairs (2 × 150 bp), 40,000 HiFi-like reads and 40,000 ONT-like reads per sample. Long-read template lengths were normally distributed (mean 15 kb; SD 3 kb for HiFi-like and 4 kb for ONT-like), truncated to 5–30 kb and 3–30 kb, respectively, before sequencing errors were simulated. Templates were sampled from one-to-one assembly-to-GRCh38 intervals; evaluation used reads with unambiguous truth positions on both GRCh38 and CHM13. Linear GRCh38 uses BWA-MEM for Illumina and minimap2 for long reads. GRCh38-graph and CHM13-graph were built from the corresponding FASTA sequences using vg autoindex; all graph arms use vg Giraffe. Misplacement means an incorrect chromosome or both alignment endpoints more than 50 bp from their corresponding true endpoints. Misplacement percentages use placed reads as the denominator; MAPQ ≥ 30 rows give read counts. In both samples, HPRC improves Illumina and HiFi-like placement relative to all three single-genome controls.

Response Table 3.14

| **read type** | **MAPQ** | **GRCh38-graph**    | **CHM13 graph**     | **HPRC v1.1**        | **HPRC v2.1**        |
| ------------- | -------- | ------------------- | ------------------- | -------------------- | -------------------- |
| Illumina      | 0        | 41,978 (70.3%)      | 41,091 (69.4%)      | 45,494 (61.2%)       | 44,253 (59.7%)       |
| Illumina      | 1–9      | 43,093 (48.8%)      | 37,091 (47%)        | 49,532 (37.1%)       | 46,648 (38.5%)       |
| Illumina      | 10–19    | 28,286 (18.9%)      | 27,747 (15.5%)      | 27,169 (5.56%)       | 27,301 (3.33%)       |
| Illumina      | 20–29    | 10,165 (12.9%)      | 9,830 (9.65%)       | 10,156 (3.57%)       | 10,220 (2.29%)       |
| Illumina      | 30–39    | 17,540 (5.87%)      | 17,028 (5.16%)      | 16,416 (1.55%)       | 16,238 (0.918%)      |
| Illumina      | 40–49    | 16,423 (4.83%)      | 15,929 (3.93%)      | 16,206 (1.37%)       | 16,175 (0.89%)       |
| Illumina      | 50–59    | 7,390 (4.14%)       | 7,498 (3.84%)       | 7,109 (1.18%)        | 7,312 (0.602%)       |
| Illumina      | 60       | 4,831,721 (0.0498%) | 4,829,923 (0.0576%) | 4,808,982 (0.00715%) | 4,810,960 (0.00611%) |
| HiFi-like     | 0        | 3 (100%)            | 11 (90.9%)          | 7 (42.9%)            | 0 (–)                |
| HiFi-like     | 1–9      | 612 (27.3%)         | 604 (21.5%)         | 707 (15.6%)          | 775 (18.8%)          |
| HiFi-like     | 10–19    | 35 (25.7%)          | 36 (19.4%)          | 46 (21.7%)           | 172 (2.33%)          |
| HiFi-like     | 20–29    | 50 (12%)            | 35 (11.4%)          | 45 (4.44%)           | 243 (0%)             |
| HiFi-like     | 30–39    | 145 (3.45%)         | 133 (3.01%)         | 125 (1.6%)           | 719 (0%)             |
| HiFi-like     | 40–49    | 474 (1.48%)         | 450 (0.667%)        | 409 (0.244%)         | 1,971 (0.0507%)      |
| HiFi-like     | 50–59    | 1,515 (0%)          | 1,479 (0%)          | 1,387 (0.0721%)      | 4,766 (0.042%)       |
| HiFi-like     | 60       | 37,146 (0.00269%)   | 36,992 (0.0162%)    | 37,226 (0.00537%)    | 31,303 (0%)          |
| ONT-like      | 0        | 1 (100%)            | 4 (100%)            | 9 (77.8%)            | 0 (–)                |
| ONT-like      | 1–9      | 259 (43.2%)         | 239 (41.8%)         | 336 (28.6%)          | 329 (52.6%)          |
| ONT-like      | 10–19    | 110 (15.5%)         | 115 (12.2%)         | 86 (8.14%)           | 133 (47.4%)          |
| ONT-like      | 20–29    | 127 (5.51%)         | 137 (3.65%)         | 173 (1.73%)          | 128 (42.2%)          |
| ONT-like      | 30–39    | 112 (1.79%)         | 118 (1.69%)         | 90 (1.11%)           | 88 (52.3%)           |
| ONT-like      | 40–49    | 64 (6.25%)          | 58 (5.17%)          | 53 (3.77%)           | 88 (45.5%)           |
| ONT-like      | 50–59    | 46 (13%)            | 38 (10.5%)          | 43 (6.98%)           | 67 (47.8%)           |
| ONT-like      | 60       | 38,999 (0.405%)     | 38,770 (0.353%)     | 38,930 (0.321%)      | 38,883 (0.465%)      |

**Response Table 3.14.** MAPQ distributions and observed misplacement for simulated reads from HG008N v6.2 on chr1 and chr20. Each cell gives the number of placed reads in the MAPQ bin, followed by the percentage misplaced within that bin. A dash indicates an empty bin. This analysis covers the four vg Giraffe graph arms; it does not include the linear GRCh38 mapper arm or COLO829BL. Each arm uses its reference-evaluable placed reads, whereas Response Table 3.13 uses the subset with one-to-one truth on both GRCh38 and CHM13. Misplacement means placement on the wrong chromosome, or both alignment endpoints lying more than 50 bp from their corresponding true endpoints. High MAPQ does not eliminate errors; Illumina misplacement is lower on the HPRC graphs across the reported MAPQ bins.

Response Table 3.15 Regional variant calling before population filtering

<table>
<thead>
<tr class="header">
<th><strong>chr1 dataset</strong></th>
<th><strong>HPRC graph</strong></th>
<th><strong>SMaHT extreme</strong><br />
<strong>PASS calls (precision)</strong></th>
<th><strong>SMaHT extreme</strong><br />
<strong>ΔF1</strong></th>
<th><strong>GIAB all difficult</strong><br />
<strong>ΔF1</strong></th>
<th><strong>GIAB low map or segdup</strong><br />
<strong>ΔF1</strong></th>
<th><strong>GIAB tandem</strong><br />
<strong>repeats and</strong><br />
<strong>homopolymers</strong><br />
<strong>ΔF1</strong></th>
<th><strong>GIAB outside difficult</strong><br />
<strong>ΔF1</strong></th>
<th><strong>INDEL F1</strong><br />
<strong>homopolymers</strong><br />
<strong>≥ 7 bp</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>HG008T<br />
Illumina</td>
<td>v1.1</td>
<td>15,955 (0.003)<br />
→ 601 (0.045)</td>
<td>+0.072<br />
[0.046, 0.099]</td>
<td>+0.147<br />
[0.120, 0.177]</td>
<td>+0.143<br />
[0.097, 0.193]</td>
<td>+0.103<br />
[0.073, 0.137]</td>
<td>+0.188<br />
[0.167, 0.210]</td>
<td>0.027<br />
→ 0.121</td>
</tr>
<tr class="even">
<td>HG008T<br />
Illumina</td>
<td>v2.1</td>
<td>15,955 (0.003)<br />
→ 599 (0.065)</td>
<td>+0.107<br />
[0.073, 0.144]</td>
<td>+0.207<br />
[0.172, 0.247]</td>
<td>+0.214<br />
[0.156, 0.279]</td>
<td>+0.133<br />
[0.094, 0.178]</td>
<td>+0.315<br />
[0.286, 0.344]</td>
<td>0.027<br />
→ 0.138</td>
</tr>
<tr class="odd">
<td>HG008T<br />
PacBio HiFi</td>
<td>v1.1</td>
<td>17,989 (0.003)<br />
→ 489 (0.070)</td>
<td>+0.111<br />
[0.077, 0.143]</td>
<td>+0.167<br />
[0.137, 0.198]</td>
<td>+0.144<br />
[0.100, 0.187]</td>
<td>+0.142<br />
[0.099, 0.186]</td>
<td>+0.186<br />
[0.165, 0.208]</td>
<td>0.027<br />
→ 0.046</td>
</tr>
<tr class="even">
<td>HG008T<br />
PacBio HiFi</td>
<td>v2.1</td>
<td>17,989 (0.003)<br />
→ 330 (0.112)</td>
<td>+0.169<br />
[0.121, 0.218]</td>
<td>+0.263<br />
[0.221, 0.308]</td>
<td>+0.198<br />
[0.139, 0.261]</td>
<td>+0.221<br />
[0.166, 0.279]</td>
<td>+0.313<br />
[0.285, 0.341]</td>
<td>0.027<br />
→ 0.056</td>
</tr>
<tr class="odd">
<td>HG008T<br />
ONT</td>
<td>v1.1</td>
<td>19,676 (0.003)<br />
→ 571 (0.058)</td>
<td>+0.094<br />
[0.064, 0.126]</td>
<td>+0.153<br />
[0.126, 0.184]</td>
<td>+0.128<br />
[0.088, 0.170]</td>
<td>+0.130<br />
[0.093, 0.171]</td>
<td>+0.179<br />
[0.159, 0.201]</td>
<td>0.010<br />
→ 0.010</td>
</tr>
<tr class="even">
<td>HG008T<br />
ONT</td>
<td>v2.1</td>
<td>19,676 (0.003)<br />
→ 669 (0.070)</td>
<td>+0.118<br />
[0.083, 0.159]</td>
<td>+0.195<br />
[0.159, 0.236]</td>
<td>+0.160<br />
[0.109, 0.227]</td>
<td>+0.148<br />
[0.107, 0.192]</td>
<td>+0.313<br />
[0.284, 0.343]</td>
<td>0.010<br />
→ 0.009</td>
</tr>
<tr class="odd">
<td>COLO829T<br />
Illumina</td>
<td>v1.1</td>
<td>8,143 (0.016)<br />
→ 3,102 (0.045)</td>
<td>+0.053<br />
[0.023, 0.165]</td>
<td>+0.101<br />
[0.052, 0.243]</td>
<td>+0.044<br />
[0.019, 0.158]</td>
<td>+0.042<br />
[0.019, 0.147]</td>
<td>+0.535<br />
[0.510, 0.557]</td>
<td>0.007<br />
→ 0.032</td>
</tr>
<tr class="even">
<td>COLO829T<br />
Illumina</td>
<td>v2.1</td>
<td>8,143 (0.016)<br />
→ 296 (0.449)</td>
<td>+0.471<br />
[0.399, 0.541]</td>
<td>+0.508<br />
[0.460, 0.554]</td>
<td>+0.463<br />
[0.381, 0.546]</td>
<td>+0.473<br />
[0.393, 0.542]</td>
<td>+0.629<br />
[0.613, 0.644]</td>
<td>0.007<br />
→ 0.041</td>
</tr>
<tr class="odd">
<td>COLO829T<br />
Fiber-seq</td>
<td>v1.1</td>
<td>5,684 (0.023)<br />
→ 844 (0.172)</td>
<td>+0.224<br />
[0.123, 0.369]</td>
<td>+0.286<br />
[0.193, 0.390]</td>
<td>+0.196<br />
[0.097, 0.362]</td>
<td>+0.323<br />
[0.268, 0.377]</td>
<td>+0.501<br />
[0.472, 0.526]</td>
<td>0.005<br />
→ 0.010</td>
</tr>
<tr class="even">
<td>COLO829T<br />
Fiber-seq</td>
<td>v2.1</td>
<td>5,684 (0.023)<br />
→ 321 (0.455)</td>
<td>+0.481<br />
[0.393, 0.570]</td>
<td>+0.503<br />
[0.447, 0.559]</td>
<td>+0.445<br />
[0.343, 0.561]</td>
<td>+0.415<br />
[0.323, 0.510]</td>
<td>+0.631<br />
[0.615, 0.646]</td>
<td>0.005<br />
→ 0.008</td>
</tr>
<tr class="odd">
<td>COLO829T<br />
ONT</td>
<td>v1.1</td>
<td>10,110 (0.012)<br />
→ 815 (0.171)</td>
<td>+0.241<br />
[0.168, 0.338]</td>
<td>+0.293<br />
[0.221, 0.375]</td>
<td>+0.217<br />
[0.134, 0.342]</td>
<td>+0.276<br />
[0.208, 0.349]</td>
<td>+0.478<br />
[0.454, 0.501]</td>
<td>0.002<br />
→ 0.002</td>
</tr>
<tr class="even">
<td>COLO829T<br />
ONT</td>
<td>v2.1</td>
<td>10,110 (0.012)<br />
→ 492 (0.274)</td>
<td>+0.348<br />
[0.254, 0.482]</td>
<td>+0.407<br />
[0.321, 0.500]</td>
<td>+0.303<br />
[0.207, 0.458]</td>
<td>+0.452<br />
[0.391, 0.511]</td>
<td>+0.595<br />
[0.580, 0.610]</td>
<td>0.002<br />
→ 0.001</td>
</tr>
</tbody>
</table>

**Response Table 3.15.** SNV calling by region on the HPRC graphs against the GRCh38-only graph (chr1, benchmark regions, before population filtering (no filtered by the four resource VCFs)). SMaHT extreme regions: PASS calls and precision on each graph and the paired ΔF1 (HPRC − GRCh38-only) with its 95% interval (2,000 resamples of 1-Mb blocks). GIAB v3.5 stratifications: all difficult regions (46.8 Mb of chr1), low mappability or segmental duplications, tandem repeats and homopolymers, and the complement of all difficult regions, with the same paired intervals. Last column: INDEL F1 in homopolymers of at least 7 bp, from the earlier INDEL models of the three-reference comparison (R3.1). HP denotes homopolymers. In the ΔF1 columns, the second line gives the 95% CI. Arrows show GRCh38-only → HPRC.

7\. The conclusions drawn from specific benchmark false negatives should be stated with greater caution, as the current interpretation may be overly definitive.  
Supplementary Figure 3 characterizes several benchmark false-negative variants and concludes that they reflect GRCh38 reference-bias artifacts rather than genuinely missed variants, on the grounds that the corresponding reads align perfectly to alternative graph paths. While this interpretation is plausible, graph-alignment behavior alone may be insufficient to overturn established benchmark truth annotations. Reclassifying these events requires orthogonal validation, such as matched-normal alignment evidence or haplotype-resolved assembly data. In the absence of such supporting evidence, the authors are advised to adopt more cautious phrasing and qualify their claims accordingly.

**Response:** We agree and have withdrawn the interpretation that these false negatives are reference-bias artefacts or benchmark errors. The original Supplementary Fig. 3 used HPRC v1.1; reassessment of the alleles defined by its original read-level criterion classified most alleles as absent from the matched normal. Here, we extend this mechanistic analysis to a strictly defined set of HPRC v2.1 alleles. We tested 1,607 HG008T (140 SNVs and 1,467 INDELs) and 401 COLO829T (54 SNVs and 347 INDELs) absorbed truth alleles using phased matched-normal assemblies (HG008N; COLO829BL), tumor assemblies, matched-normal reads and HPRC haplotype paths (Response Table 3.16). These alleles are spelled exactly by an HPRC v2.1 path but lose their candidate on at least one platform despite having one on the GRCh38-only graph. Only 6.8% and 6.7%, respectively, are carried by the matched normal, while 92.2% and 92.5% are absent from both normal haplotypes. Of the 1,351 HG008T and 319 COLO829T INDEL loci absent from the normal, the tumor assemblies carry exactly the benchmark allele at 1,153 (85.3%) and 226 (70.8%), respectively; the COLO829T assembly uses the same tumor long reads analyzed here and provides consistency rather than independent validation.

These results are consistent with somatic variants lost through graph absorption and do not justify reclassifying these alleles as truth-set errors. Supplementary Fig. 3 is now titled "Benchmark somatic truth alleles missed through graph absorption". We therefore retain these alleles in the benchmark truth set and count them as false negatives wherever Pansoma misses them within the evaluation regions.

Response Table 3.16

| **Evidence**                                                                                                       | **HG008T**                                    | **COLO829T**                                                 |
| ------------------------------------------------------------------------------------------------------------------ | --------------------------------------------- | ------------------------------------------------------------ |
| Absorbed truth alleles (SNVs / INDELs)                                                                             | 1,607 (140 / 1,467)                           | 401 (54 / 347)                                               |
| ALT carried by a matched-normal haplotype                                                                          | 120 (7.5%)                                    | 35 (8.7%)                                                    |
| ALT absent from both normal haplotypes                                                                             | 1,461 (90.9%)                                 | 363 (90.5%)                                                  |
| Normal status ambiguous                                                                                            | 26 (1.6%)                                     | 3 (0.7%)                                                     |
| Normal-absent INDELs: tumor assembly carries exactly the truth ALT                                                 | 1,148 / 1,346 (85.3%)                         | 226 / 311 (72.7%)ᵃ                                           |
| Normal-absent INDELs: tumor reads show the expected allele, % of evaluable loci (random non-absorbed truth INDELs) | Illumina 97.4 (95.7); PacBio HiFi 97.5 (96.4) | Illumina 97.6 (98.4); Fiber-seq 97.3 (97.6); ONT 79.8 (82.1) |
| Normal reads: REF-like reads on the ALT pathᵇ                                                                      | 427 / 58,736 (0.7%)                           | 367 / 33,418 (1.1%)                                          |
| Normal reads: ALT-like reads on the ALT pathᵇ                                                                      | 42,753 / 44,207 (96.7%)                       | 20,288 / 21,285 (95.3%)                                      |
| Normal-absent SNVs: reads carrying the ALT, normal / tumorᵇ                                                        | 24.7% / 59.1%ᶜ                                | 2.0% / 45.5%                                                 |

**Response Table 3.16.** Assembly and read evidence for truth alleles absorbed by HPRC v2.1 paths (chr1–22). An allele is absorbed when an HPRC v2.1 path spells its ALT exactly and it has a candidate on the GRCh38-only graph but none on HPRC v2.1 on at least one platform. Normal-status percentages use all absorbed alleles; tumor-assembly and tumor-read rows use the normal-absent INDELs, with random non-absorbed truth INDELs as controls in parentheses. Normal-read rows give pooled read counts. Details are given in Supplementary Note R2.3.

ᵃ Built from the analysed HiFi and ONT tumor reads, so it shows consistency rather than independent support. ᵇ PacBio HiFi for HG008-N and ONT for COLO829BL, realigned with Giraffe; tumor reads were processed identically. ᶜ The HG008T truth set was called against the HG008-N assembly and projected to GRCh38, so in repeats a truth SNV can have an ALT base that the normal already carries; at 300 random non-absorbed truth SNVs, 0.2% of normal reads carry the ALT. The COLO829T truth set was called on GRCh38 and is not affected.

8\. The validation is predominantly tumor-centric; applicability to somatic mosaicism in normal tissue remains unresolved  
The manuscript introduces Pansoma broadly as a somatic variant caller, yet its biological validation rests on tumor-derived samples and a synthetic low-VAF admixture. COLO829T is a melanoma-derived benchmark and HG008T a pancreatic cancer sample, while the low-VAF HapMap benchmark comprises germline variants from minor cell-line constituents employed as controlled proxies for somatic-like events. No authentic normal-tissue somatic-mosaicism dataset such as SMaHT tissues appears to be evaluated. This distinction is material, as somatic variants arising in normal tissues occupy a substantially different biological context and frequently present at low variant allele fractions.  
If Pansoma is intended specifically for tumor-only somatic variant calling, the manuscript should articulate this scope consistently and refrain from broader claims. If, conversely, the authors intend general applicability to somatic variation beyond malignancy, validation on an independently characterized normal-tissue somatic-mosaicism dataset would considerably strengthen the work. 

**Response:** We agree and have narrowed the scope to tumor-only somatic variant calling in cancer samples; the title, Abstract, Results and Discussion were revised and claims about normal-tissue mosaicism removed. Considering the related concerns raised by Reviewers 1 and 2, we agree that the HapMap mixture, whose positives are minor donors' germline alleles, is less suited to Pansoma's tumor-only somatic-calling scope. We have therefore moved this analysis to the Supplement as an exploratory low-VAF experiment.

We expanded the training set by adding HCC1395T and HCC1954T (breast cancer), and H1437T and H2009T (lung adenocarcinoma). Our cancer benchmarks cover HG008T (pancreatic cancer), COLO829T (melanoma) and HCC1395T (breast cancer). Pansoma has not been trained or validated for somatic mosaicism in normal tissue. Its main models are trained on tumor somatic variants and require at least three ALT reads and an allele fraction of at least 0.08, so the much lower allele fractions typical of mosaicism lie outside their scope. The exploratory low-VAF model (candidate floor 0.015) is reported only as allele detection in the HapMap mixture. Lower tumor fractions are tested only within the cancer scope, on in-silico purity dilutions of COLO829T, and are not extrapolated to normal tissue. The Discussion names validation on characterized normal-tissue data as future work.

9\. Assessment of Model Stability and Statistical Uncertainty  
The selected performance metrics are appropriate; however, the presentation relies almost exclusively on point estimates. Given the inherent susceptibility of machine learning models to variation during selection and training, a formal assessment of model stability is necessary. The authors are encouraged to evaluate uncertainty using approaches such as multi-seed retrainings, block-bootstrap confidence intervals across genomic intervals, or chromosome-wide variance estimates to ensure the robustness of the reported performance metrics.  
A particularly informative analysis would be to retrain Pansoma across multiple random seeds and report the distribution of precision, recall and F1 on both the held-out chromosome and the independent HG008T benchmark. In addition, chromosome-level or genomic-block re-sampling could be used to provide confidence intervals for the performance difference between Pansoma and the comparator callers, rather than reporting only point estimates. Sensitivity of the results to the fixed 0.5 classification threshold would also help establish that the reported advantage is not dependent on a single operating point.

**Response:** We agree that single point estimates are insufficient. The revision assesses variation across training seeds, performance across genomic regions and sensitivity to the classification threshold.

1.  *Training seeds.* SNV reference models and final INDEL models were trained with three seeds and tested on held-out chr1 of HG008T and COLO829T (Response Tables 3.1 and 3.2). The standard deviation of F1 across seeds was 0.001–0.042 (median 0.012) for HPRC SNV models and 0.003–0.027 (median 0.014) for final INDEL models. Across the three seeds, SNV and INDEL performance was broadly stable, with no large shifts observed (Response Tables 3.1 and 3.2).

2.  *Performance across genomic regions.* To explore the performance difference between Pansoma and the comparator callers, we evaluated Pansoma across different genomic regions, including low-mappability and segmentally duplicated regions, tandem repeats, homopolymers and SMaHT extreme regions, to assess its performance in specific genomic contexts (Response Table 3.15).

3.  *Tumors excluded from training.* To assess performance independently of the training tumor, we evaluated Pansoma models across tumor datasets, including tumors excluded from training.

4.  *Operating point.* We added germline variants as a third label in Pansoma's training sets and replaced the fixed 0.5 threshold with validation-selected PASS thresholds. We also generated precision–recall curves by varying the classification threshold to assess its effect on performance. \[Figure\] 

10\. Need for Qualified Interpretation of Low-VAF and Insertion/Deletion Benchmarks  
The utility of the HapMap benchmark for low-VAF assessment is acknowledged; nevertheless, because these synthetic constructs rely on mixed germline lines rather than genuine somatic lineages, the generalizability of the low-VAF results should be stated with appropriate caution. Additionally, the manuscript should clearly distinguish SNV performance from INDEL accuracy. Pansoma's competitive edge is substantially more pronounced in SNV calling, whereas INDEL accuracy remains modest and is surpassed by ClairS-TO in certain ONT settings. Overall performance claims should be refined to accurately reflect these boundaries.  
To better define these limits, the authors could complement the HapMap mixture with either a bona fide low-VAF somatic benchmark or an in silico dilution series based on validated somatic variants from COLO829T or HG008T. In addition, INDEL performance should be stratified by insertion/deletion class, event length, sequence context such as homopolymers or tandem repeats, and sequencing platform. These analyses would help determine whether the observed limitations reflect low-VAF biology, platform-specific error modes, or particular classes of indels.

Response. We agree with all four points and have revised the manuscript accordingly.

1.  The HapMap mixture moved from the main benchmark to the Supplement and is reported only as detection by allele fraction in that mixture, not as somatic accuracy.

2.  SNVs and INDELs are reported separately throughout. The revised pipeline substantially improves INDEL performance in several settings as shown in the manuscript. The six-tumor models also outperform both comparator callers in three of the six HG008T and COLO829T tumor–platform settings. Nevertheless, performance remains platform-dependent, and we do not claim a consistent INDEL advantage over other callers. The revision reflects these improvements and limitations.

3.  For somatic variants at reduced purity, we added in-silico dilutions of COLO829T with its matched normal (75, 50 and 25%) on three platforms, with recall by expected VAF.

4.  A new Supplementary Table stratifies INDEL calls by insertion or deletion, length, repeat context and platform for both tumors, with truth, TP, FP and FN per stratum and both comparators. INDEL calling is limited mainly by three factors: a. population filtering (gnomAD and CoLoRSdb VCFs only) of truth INDELs (384 of 587 HG008T chr1 truth INDELs in the benchmark regions are population alleles, 294 of them in homopolymers of at least 7 bp), b. errors in long homopolymers (44–70% of HG008T false calls), and c. the threshold. The dilutions contain too few true INDELs (131 per purity) to assess allele-fraction effects.

Minor and technical comments  
1\) Additional methodological details and explicit reconciliation are required  
The methodological description appears disproportionately focused on the neural-network architecture, whereas several upstream and downstream steps that directly influence calling performance are comparatively condensed.  
More importantly, several key parameters differ between the main Methods and Supplementary Methods.  
The main Methods describe:  
non-reference AF \>10%,  
2 supporting reads,  
exclusion of reads with MAPQ ≤5.  
The Supplementary Methods instead describe:  
AF ≥0.08,  
≥3 alternative observations,  
mean alternative BQ ≥10,  
MAPQ ≥10,  
INDEL length ≤50 bp.  
Similarly, the main Methods describe 80 training epochs, whereas the Supplementary Methods state 80–120 epochs.  
Resolving these discrepancies as well as more detailed description are essential, as they directly impact the study’s reproducibility and the overall interpretation of model performance.

We agree that the inconsistent parameter descriptions hindered reproducibility. The two sets of values described different analyses: the main Methods listed the thresholds of the original candidate-count comparison between GRCh38 and the graph (original Supplementary Fig. 5), whereas the Supplementary Methods listed those of tensor generation. We have reconciled the main and Supplementary Methods against the configuration files used for the revised analyses and expanded the descriptions of candidate generation, filtering and evaluation. Both now give the values in Response Table 3.17, which supersede those of the original submission.

Response Table 3.17

| **Step**             | **Standard setting**                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| -------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Node selection       | Alignments with MAPQ \>5; a node is selected when more than 5% contain a mismatch or INDEL on that node.                                                                                                                                                                                                                                                                                                                                             |
| Candidates           | Reads with MAPQ \>10; at least three ALT-supporting observations with base quality ≥10 per observation; allele fraction ≥0.05 for SNVs (≥0.08 for HG008T and all ONT datasets) and ≥0.08 for INDELs; INDEL length 1–50 bp; at most 800 reads per node.                                                                                                                                                                                               |
| Tensors              | Eight channels × up to 200 reads × 101 columns following each read’s graph walk. SNV models also take population site features of the full graph as the input.                                                                                                                                                                                                                                                                                       |
| Labels and splits    | Somatic, germline and other classes. The six-tumor models train on chr2–22 of HG008T and COLO829T, chr2–19 and chr21–22 of HCC1395, and all autosomes of the three CASTLE training lines. Five percent of 20,000-node blocks are reserved for validation.                                                                                                                                                                                            |
| Training             | ConvNeXt: 50.4 M parameters and 20 training epochs for SNVs; 8.7 M parameters and 12 epochs for INDELs. AdamW, learning rate 10⁻⁴, weight decay 0.05, one warm-up epoch followed by cosine decay, and batch size 1,024 per GPU. Cross-entropy class weights are proportional to inverse class frequency raised to 0.75 for SNVs or 0.5 for INDELs. For SNV models, “other” examples with allele fraction below 0.1 are sampled at 25% in each epoch. |
| Checkpoint           | The standard rule selects the epoch with the highest validation somatic F1.                                                                                                                                                                                                                                                                                                                                                                          |
| Population filtering | Exact allele matching. SNVs: gnomAD or CoLoRSdb AF ≥10⁻⁴, non-somatic dbSNP 138 records, and the GATK 1000 Genomes panel of normals. INDELs: gnomAD or CoLoRSdb AF ≥0.01 in homopolymers ≥8 bp, otherwise ≥10⁻⁴.                                                                                                                                                                                                                                     |
| Evaluation           | GRCh38-projectable calls evaluated with rtg vcfeval on chr1 within benchmark regions; graph-coordinate calls assessed separately.                                                                                                                                                                                                                                                                                                                    |

**Response Table 3.17.** Standard settings of the revised Pansoma pipeline.

2\) The setup and parameters for comparator workflows should be documented in greater detail to ensure fair benchmarking.  
Pansoma is described in considerable technical detail, including the vg giraffe version, but the comparator workflows are largely described as their “standard GRCh38-based alignment and calling workflows.”  
For reproducibility and fair benchmarking, the authors should provide:  
software versions,  
aligners and versions,  
pretrained checkpoints/models,  
exact tumor-only configurations,  
PASS/filtering criteria,  
reference resources,  
population/PoN resources,  
command-line parameters.  
This is particularly important because the main advantage reported for Pansoma is increased precision, which can be sensitive to filtering configuration.

The revised Methods document comparator versions, input alignments, reference resources, pretrained models and filtering settings. Response Table 3.18 lists the alignment software and the pretrained model for each dataset, and the paragraphs that follow give the software versions, the options we set, the population resources and the scoring.

Response Table 3.18

| **Platform**            | **Alignment software: HG008T / COLO829T / HCC1395**          | **DeepSomatic-TO model** | **ClairS-TO model**                                           |
| ----------------------- | ------------------------------------------------------------ | ------------------------ | ------------------------------------------------------------- |
| Illumina                | BWA-MEM 0.7.17 / Sentieon BWA-MEM 202308.01 / BWA-MEM 0.7.17 | WGS\_TUMOR\_ONLY         | ilmn\_ssrs; COLO829T used ilmn, corresponding to the SS model |
| PacBio HiFi / Fiber-seq | pbmm2 1.10.0 / 1.13.0 / 1.10.0                               | PACBIO\_TUMOR\_ONLY      | hifi\_revio\_ssrs                                             |
| ONT                     | minimap2 2.26                                                | ONT\_TUMOR\_ONLY         | ont\_r10\_dorado\_sup\_5khz\_ssrs                             |

**Response Table 3.18.** Alignment software and pretrained comparator models for each dataset.

**Software and inputs.** We used DeepSomatic 1.9.0 and ClairS-TO 0.4.2 (0.4.1 for COLO829T), run from the developers’ container images on the providers’ GRCh38 alignments of the reads used by Pansoma; we aligned the SEQC2 Illumina reads ourselves. The calling reference was GRCh38 (the no-alt analysis set for HCC1395), and calls were evaluated on chr1. HG008T ONT reads were basecalled with dorado SUP v5.0; ClairS-TO provides no SUP v5 model, so its SUP v4.2 model was used. For HG008T ONT, DeepSomatic-TO was rerun with ONT\_TUMOR\_ONLY, replacing the original submission’s run with PACBIO\_TUMOR\_ONLY.

**Population resources and scoring.** DeepSomatic-TO used its population-frequency input resources, and we enabled its panel-of-normals filtering (--use\_default\_pon\_filtering=true). ClairS-TO used its default quality filters and population tagging: gnomAD and CoLoRSdb AF ≥0.001, non-somatic dbSNP 138 records, and the GATK 1000 Genomes panel of normals. All other options were left at their defaults. PASS denotes FILTER=PASS.

All callers were evaluated against the same truth sets and benchmark regions with the same rtg vcfeval options, and we report PASS calls after the same additional population filter. Comparator precision–recall curves rank by QUAL all records that remain after this filter, including non-PASS variants and excluding reference/no-call records.

**Training data.** Pansoma’s within-tumor HG008T and COLO829T models use chr2–22 of the evaluated tumor; cross-tumor evaluations exclude the evaluated tumor from training. On HCC1395, DeepSomatic-TO and Pansoma’s multi-tumor models include the line in training with chr1 held out; ClairS-TO’s SSRS models exclude HCC1395 entirely.

3\) The composition of the pangenome reference should be specified in greater detail  
The manuscript refers throughout to an allele-frequency-filtered HPRC graph, and Figure 5 specifies the HPRC v1.1 AF-Graph. Because the pangenome constitutes the central methodological element, the Methods should document the HPRC release and version, the haplotypes included, the allele-frequency filtering criteria, the treatment of structural variants, and the graph-construction procedure. It should further be stated whether any benchmark individuals or haplotypes closely related to them, contributed to the graph used for evaluation. Because the manuscript notes that diagnostic efficacy hinges on graph composition and haplotype representation, disclosing this information is essential to assessing the external validity and translational utility of the proposed pipeline.

The revised Methods specify each graph’s release, construction procedure, haplotype composition and filtering rules.

  - **HPRC v1.1, d9:** the release-1 Minigraph-Cactus graph based on GRCh38 (Liao et al., 2023), containing 88 sample haplotypes from 44 individuals plus GRCh38 and CHM13, giving 90 in total.

  - **HPRC v2.1, d46:** the corresponding release-2 graph (Lucas et al., 2026), containing 462 sample haplotypes from 231 individuals plus GRCh38 and CHM13, giving 464 in total.

  - **Frequency filtering:** we used the graphs as distributed by HPRC, with support thresholds of nine and 46 haplotypes, respectively. The GRCh38 and CHM13 reference paths are retained irrespective of these thresholds. The allele-carrier analysis excludes GRCh38 and therefore uses 463 haplotypes for v2.1.

  - **GRCh38-only control:** constructed with vg autoindex from the GRCh38 no-alt analysis set (R3.1).

Structural-variant sequences are subject to the same graph-filtering procedure and remain available for alignment when retained. Pansoma itself calls only SNVs and INDELs of 1–50 bp.

None of the six cancer lines appears in either graph’s sample list. We additionally examined rare-allele sharing for the HG008, COLO829 and HCC1395 normals. Among their germline SNVs carried by one or two graph haplotypes, the maximum fraction shared with an HPRC individual was 8.2% for v1.1 and 2.9% for v2.1. These measurements provide supplementary evidence but do not formally exclude relatedness. The remaining three lines were checked by sample list only. The HapMap experiment used a separate graph excluding all six mixture donors.

Graph names have been standardized throughout the manuscript to identify the release and filtering threshold explicitly.

4\) Terminology and Language Consistency  
Throughout the manuscript and Supplementary Information, terminology is not always used consistently, and several minor grammatical and stylistic issues are present. The authors are encouraged to perform a careful language and consistency check, particularly for method names, abbreviations, figure/supplementary nomenclature, and technical terminology.  
The term “panel of normals” is used broadly for resources that include gnomAD, dbSNP, a 1000 Genomes PoN and CoLoRSdb. These are not all conventional PoNs; “population/background filtering resources” may be more precise.  
“One key challenge in performing variant calling on graph alignments is it differs fundamentally…” is grammatically strange and incorrect.   
“Pansoma convert the variants…”   
“with the technique details described in Methods”   
“read pile up” in Figure 5  
The expansion and terminology for NPU should be standardized between the main text and Supplementary Methods.  
Pansoma versus Pansoma-TO, and “Supplement Figure” versus “Supplementary Figure”, should be standardized throughout.  
The manuscript would benefit from a careful final language edit

We have standardized terminology throughout the manuscript and Supplementary Information. The four background resources are now described collectively as “population filtering resources.” “Panel of normals” is used for the GATK 1000 Genomes resource and, where appropriate, the comparators’ own panel-of-normals filters. Each resource and its filtering rule are specified in Methods.

The revised pipeline constructs tensors directly from sorted, indexed GAM alignments and no longer produces NPU files. We have therefore removed the NPU terminology and obsolete format description from the original Supplementary Table 1.

We consistently use “Pansoma,” “Supplementary Fig.” and “Supplementary Table,” and have standardized reference names, platform labels and candidate terminology. The grammatical errors identified by the reviewer have been corrected, including “Pansoma converts,” “technical details” and “read pileup.” The sentence concerning differences between graph and linear-reference alignments has also been rewritten for clarity. The manuscript and Supplementary Information have been edited throughout for language and consistency.

Overall assessment  
In summary, I consider Pansoma a novel and potentially valuable framework for graph-native somatic variant calling. The improvement in tumor-only SNV precision is encouraging, and several analyses support the proposed benefit of pangenome-aware read representation.  
The primary issue is whether the manuscript currently demonstrates sufficiently and directly that this advantage is specifically attributable to the pangenome itself, rather than to the broader Pansoma architecture and workflow.  
A direct reference-level ablation, a clearer accounting of which candidates are removed by the graph, stronger validation of graph-coordinate calls, and reconciliation of the methodological inconsistencies would substantially strengthen the study.  
In addition, the current validation is primarily cancer-centered. Either the scope should remain tumor-only, or broader applicability to somatic mosaicism should be supported by additional data.  
I therefore recommend major revision.
