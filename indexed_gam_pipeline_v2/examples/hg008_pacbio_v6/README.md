# HG008 PacBio examples in tensor format v6

The v6 counterparts of the 100 v5 examples in `../hg008_pacbio_v5/`: each of the 98
production batches holding a v5 example was rebuilt with the v6 code (`--debug-rows`,
Slurm 363665; every output passed `validate_examples.py`, 196/196) and the site now holding
the example was rendered as `<v5 name>_v6.png`, so each image can be compared with its v5
original side by side. Regenerate with `make_examples.py <rebuilt batches> <supplement outputs...>`.

v6 changes visible here: indels left-normalized (both strands put a repeat indel at the same
place), one tensor per site with every passing allele, channel 2 (panel 3) spells the allele
each row carries (A1, A2, ... or REF; blank = OTHER), rows in blocks A1.., REF, OTHER (labels
and separators on panels 1 and 3) ordered by similarity inside each block, and path counts
exact up to 100 (panel 7).

`index.tsv`: the match (`same position`; `left-normalized` = the indel moved left on its node;
`left-normalized to another node` = onto a neighbouring node, `(supplement run)` when that node
is not a target and the site was built by the supplement run), the v6 site with all its alleles
and AFs, site coverage and block counts, the strand of the A1 rows (`a1_one_strand` = all of
>= 3 A1 rows on one strand) and the v5 values. The rebuild: 98 batches (Slurm 363680, main run
byte-identical to the previous v6 build plus `displaced_nodes.tsv`), then a supplement run over
the 4,098 non-target nodes that normalized indels landed on (576 INDEL sites; 44/44 outputs
pass `validate_examples.py`); no site appears twice (7,589 INDEL site ids).

98 of the 100 examples have a v6 site. The two without (`082`, `086`, both 1-bp insertions in
T homopolymers chopped into short nodes) lost their site for real: in `082` reads on different
graph branches normalize the same insertion to two places (6 and 3 reads; neither passes the
filters), in `086` the insertion's AF on its normalized site is 6/110 = 0.055, below the INDEL
threshold 0.08 (v5: 10/110). Over the 98 batches 111 of 6,028 v5 INDEL sites (1.8 %) have no
equivalent allele in v6, mostly at AF < 0.2 (see the v2 README, "supplement run").
