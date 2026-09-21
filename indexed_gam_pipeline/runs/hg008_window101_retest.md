# HG008 anchor-window regression and width-101 test

User-approved representation: seven channels, 200 rows, **101 columns**. Anchor
column is zero-based 50. Each row uses its original GAM edits and oriented graph
path. Keep 50 alignment columns before the anchor, the anchor, and 50 following
columns. Long I/D context is cropped at the window edge. Missing coverage is
all-zero; a deletion is an explicit D with reference base, MAPQ and GBWT count.

Rows are ranked by substitution + insertion + deletion bases in the displayed
window, stably grouped by visible oriented node path in first-occurrence order,
then uniformly sampled in order to at most 200. Coverage/support/AF are computed
from all eligible records. Revisions are saved separately from the v4 channel
schema: `anchor-centered-columns-v1` and `window-edit-bp-group-uniform-v1`.

## Verification

- 40 regression tests passed, including native GBWT integration.
- Synthetic cases cover 100-bp DEL display (the production candidate filter
  remains <=50 bp), 150-bp insertion, background insertion inside a deletion,
  reverse orientation, partial coverage, 50/anchor/50 coordinates, and 401-to-200
  stable grouped sampling across all seven channels.
- Slurm **362251**, source **d4a9c91**, 2 CPUs, 64 GiB requested RAM, 4-hour limit.
- Outputs and per-step timings:
  `/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/run/window101_retest_20260921T172347Z`

`retest_windows.py` runs the original batch 111 against its original source
snapshot to capture the exact failure candidate and insertion lengths. An
expected exit status 1 with a matching diagnostic is required. It then builds
and audits all passing candidates from this batch using the corrected window,
requires the original failure candidate to appear, and proceeds to 100 examples
from the first 1,024 discovery nodes. The NPY shard setting remains 2,048, so
these bounded tests produce partial final shards. Every source-audited test
output is also checked structurally against `(7,200,101)` and `int32`.

The old graph index/discovery and failed production output are preserved. This
is a bounded regression and HG008 test, not a claim of full production completion.
See the run's `status.json`, `*.resources.txt`, `*_audit.json`, and tensor
`batch_timing.ndjson` for actual outcomes and timings.

## Earlier row-order-only test

Job 362249 (100-column legacy window) generated 95 tensors from its 640-node
subset and passed both source and structural audits. Its controller reported
failure because it expected exactly 100, but the subset had only 95 passing
candidates. Build: 336.78 s; audit: 44.46 s. The width-101 test uses 1,024 early
nodes to provide enough candidates. Earlier outputs are not mixed with the new
window representation.

## Reproduced failure and corrected batch

Original candidate: `22739:0:DEL:AACCCCGGGAACCGCCT>` (17-bp deletion).
Read `m84039_240114_012401_s1/77992551/ccs` contributes 91 insertion bases
inside that reference interval, so the former shared block required 108 columns.
Original failure reproduced in 105.65 s. The corrected batch generated **32**
tensors in **44.05 s** and passed its source audit in **30.86 s**, including the
failure candidate. Its coverage/ALT/REF/other counts remain **51/20/28/3**, exactly
matching the pre-error original computation. The corrected batch audit checked
1,103 selected rows and 97,554 reference-base columns.

The 100-example stage also passed: build **149.48 s**, source audit **99.94 s**;
82 DEL, 14 INS, 4 SNP. Output shard shape `(100,7,200,101)`, int32. It used
66 cached GAM group hits and 66 misses, decoding 66,000 indexed records. Peak
build RSS reported by GNU time: 6,506,416 KiB. These measurements include the
bounded run's startup/I/O conditions and do not establish whole-genome throughput.

`continue_hg008.py` resumes full tensor building only after validating that the
regression passed, core source files match the tested snapshot, the original
input fingerprints are unchanged, and the full graph/discovery completed. It
uses 512-node batches with an 8 GiB GAM cache as in the verified 100-example test,
2,048-example shards, width 101, no max-tensors/read/node limit. It writes to
`.building_tensors_window101` and validates before publishing into the original
output directory. Previous failures and tests are retained separately. Automatic
restart from a partial tensor build is still not implemented.

## Full continuation submission

Full build resumed as Slurm **362252** on `tsingtao`, source commit **91ece2c**,
4 CPUs, 128 GiB requested RAM, 14-day limit. At startup vg version verification
completed and full tensor building entered `running`. This is submission/startup
confirmation, not a completion claim.

Current run directory:
`/scratch/jshen/data/HG008_GIAB/pansoma_v2_tensors/Liss_lab_PacBio_Revio_20240125/run/full_window101_20260921T173318Z`

`run/active_run.json` points to this job and its current `status.json`; the old
`run/status.json` remains the historical failed job record. New tensors are
built in `.building_tensors_window101/` and are only published to the requested
output directory after full structural validation. The full shard shape is
`(2048,7,200,101)`, with a smaller final shard permitted. Full discovery is reused;
all 16,560,350 selected nodes are passed to the builder with no max-tensors limit.
