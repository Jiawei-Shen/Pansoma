# HG008 PacBio examples in tensor format v6

The v6 counterparts of the 100 v5 examples in `../hg008_pacbio_v5/`: each of the 98
production batches holding a v5 example was rebuilt with the v6 code (`--debug-rows`,
Slurm 363665; every output passed `validate_examples.py`, 196/196) and the site now holding
the example was rendered as `<v5 name>_v6.png`, so each image can be compared with its v5
original side by side. Regenerate with `make_examples.py <rebuilt batches>`.

v6 changes visible here: indels left-normalized (both strands put a repeat indel at the same
place), one tensor per site with every passing allele, channel 2 (panel 3) spells the allele
each row carries (A1, A2, ... or REF; blank = OTHER), rows in blocks A1.., REF, OTHER (labels
and separators on panels 1 and 3) ordered by similarity inside each block, and path counts
exact up to 100 (panel 7).

`index.tsv`: the match (`same position`, `left-normalized` = the indel moved left, `nearest
site`), the v6 site with all its alleles and AFs, site coverage and block counts, the strand
of the A1 rows (`a1_one_strand` = all of >= 3 A1 rows on one strand) and the v5 values.

**Known issue — 23 of the 100 examples have no v6 image (`match = none`), all indels.** Tensors
are only built for discovery's target nodes, which were chosen from vg's *raw* indel
placement. Left-normalization moves some indels onto a neighbouring node that is not a target
(e.g. `034`: DEL T at node 41415778 moves across the single-base T nodes 41415777/41415776,
neither a target), and those candidates are dropped; conversely indels that vg placed on
non-target nodes can move into a target node and appear as new sites. The target-node
selection has to be made consistent with normalization before regenerating the data.
