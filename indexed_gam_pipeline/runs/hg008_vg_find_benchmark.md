# HG008 vg find versus Python GAM reader — 2026-09-21

Job 362263 compared `/scratch/jshen/bin/vg_v1.77.0 find` against the current
`IndexedGam` on the same first 1,024 discovery nodes, in two 512-node batches.
Execution order: vg twice (one query per batch), Python twice with one persistent
8 GiB group cache, then vg twice again. One allocation, two CPUs, 24 GiB memory;
OMP_NUM_THREADS=1 and OPENBLAS_NUM_THREADS=1. OS/storage cache was uncontrolled.
Production code was not changed and full job 362252 was not resumed.

## Method

VG queries use `find -l INPUT.gam -o MIN:MAX`, writing a temporary GAM. Python
then decodes that result and filters membership against the exact target node
set, because the ID interval can include additional nodes. No MAPQ filter is
applied before completeness comparison. Counts after MAPQ >10 are also recorded.

Both readers' records are serialized deterministically and checked using:

- SHA-256 multiset with occurrence counts (preserves duplicate records).
- Length-delimited ordered SHA-256 stream (checks record order).

VG totals include process startup, querying, compressed GAM output, output
parsing, exact-node filtering and the comparison hashes. Python totals include
fetching, equivalent membership checks and hashing. Python reader construction
is outside its per-batch timing; VG startup is inside. Graph queries, edit
expansion and tensor building are not part of this benchmark. No genome-wide
completeness or speedup claim is made from these two batches.

## Results

| Reader | Batch 1 | Batch 2 | Total |
|---|---:|---:|---:|
| vg before | 5.574 s | 2.859 s | 8.433 s |
| Python, persistent cache | 98.054 s | 0.397 s | 98.451 s |
| vg after | 8.361 s | 3.046 s | 11.407 s |

VG command-only times were 5.098/2.433 seconds before and 7.883/2.410 seconds
after. Peak VG subprocess RSS was about 64–66 MiB; this excludes the Python
controller/output parser and is not an end-to-end memory peak. Python's cache
accounted for 5,204,194,223 bytes; accounted cache size is not process RSS.

Every query returned the same records, in the same order and with the same
duplicate multiplicities:

| Batch | Node interval | Exact matched records | MAPQ >10 |
|---|---|---:|---:|
| 1 | 128–2430 | 612 | 612 |
| 2 | 2431–4960 | 549 | 190 |

For these particular ranges no extra records remained to be removed by exact
node filtering; this must not be assumed for other batches. Repeated VG output
files differ at the compressed-byte level, but parsed deterministic records are
identical.

## Interpretation

VG was **8.6–11.7x faster for the combined two-batch query** in this test.
Its advantage is in avoiding the current reader's costly initialization;
Python's already-cached second batch is faster than invoking VG again.
The Python reader indexed 66,000 records across 66 GAM groups before selecting
612 records in batch 1, and reused all 66 groups in batch 2. This benchmark does
not instrument VG's internal group/record visitation, so the relative effects
of its query pruning and native implementation are not established here.

These results justify prioritizing a VG-backed fetch prototype over implementing
a complex multiprocessing replacement for the Python reader immediately.
Before production replacement, test more consecutive and separated/high-depth
regions, retain exact-node filtering, verify record order/multiplicity and final
tensor/metadata equality, and compare sustained query time with Python's cache.
No production reader replacement or full-run restart was performed.

## Records

- [Commands, timing, counts and complete comparison hashes](hg008_vg_find_benchmark/status.json)
- [Benchmark script](hg008_vg_find_benchmark/compare.py)
- VG GNU time resource records are in the same directory.
- Temporary GAM outputs: `tmp/vg_find_benchmark/` (not committed).
