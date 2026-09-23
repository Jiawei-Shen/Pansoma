# 100 HG008 PacBio v5 production tensors

Random sample (seed 20260923) of the tensors written by job 363651
(`indexed-gam-candidate-v5`, site units, 800-record cap, SNV AF ≥ 0.06, INDEL AF ≥ 0.08):
34 SNP (`000`–`033`), 33 DEL (`034`–`066`), 33 INS (`067`–`099`). Regenerate with
`make_examples.py` (matplotlib interpreter, see the v2 README, "Rendering PNGs").

`index.tsv` columns: candidate identity, task, coverage / ALT / REF / other counts, AF,
selected rows, and the strand of the rows whose pixels show the candidate ALT
(`alt_pixel_rows`, `alt_pixel_rows_reverse`, `alt_reverse_fraction`; SNP = X op with the ALT
base and BQ ≥ 10 at the anchor, INS = I ops spelling the ALT over the candidate columns,
DEL = D ops over them), plus the reverse fraction of all selected rows.

Known issue visible here: in 14 of the 33 INS examples every ALT row comes from one strand
(0 of 34 SNP, 0 of 33 DEL). vg places an indel at one end of a repeat in *read* orientation,
so in forward node coordinates forward and reverse reads put the same repeat indel at
opposite ends (sometimes on different nodes of a repeat chopped into short nodes) and the
pipeline, which has no indel normalization, keeps them as two strand-specific candidates.
The old `.dat/.idx` pipeline has no normalization either.
