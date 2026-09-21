# HG008 GAM reader diagnosis — 2026-09-21

Production reader and tensor code are unchanged. The full tensor job 362252
remains stopped. These are read-only Slurm diagnostics, not tensor builds.

## Application-cache initialization (job 362261)

The reader queried the first 512 discovery nodes, then the next 512 using the
same reader and its existing 8 GiB GAM-group cache. An in-memory copy of
`_indexed_group` added phase timers; the report records the source SHA-256.
This instrumentation adds small unquantified overhead. OS/storage cache state
was not controlled.

| First-batch activity | Wall seconds |
|---|---:|
| BGZF group reading, decompression and message framing | 68.279 |
| Protobuf alignment decoding | 8.242 |
| Extract node memberships and build postings | 20.192 |
| Convert posting lists to tuples | 0.609 |
| Calculate cache size | 3.965 |
| Whole fetch plus diagnostic consumer | 101.867 |

There were 47 merged ranges, 66 groups, 66,000 source records and 612 selected
records. Reading the source records generated 3,829,019,863 raw message bytes;
the node index contained 36,573,703 postings and 521,717 node keys across groups.
The cache accounted for 5,204,194,223 bytes (a conservative estimate, not RSS).

Second batch: 66 group-cache hits, 549 selected records (190 pass MAPQ >10),
0.183 seconds including 0.100 seconds of diagnostic serialization/hashing.
No cold-cache work repeated. Ordered record hashes are saved in the JSON.

The builder's `gam_fetch_seconds` also includes caller-side MAPQ/chromosome
filtering and context-node collection inside its fetch loop; this reader-only
diagnostic does not perform that entire caller loop.

## BGZF-only follow-up (job 362262)

Same ranges, without parsing alignments or constructing postings:

- BGZF grouping: 49.433 seconds wall, 23.020 seconds process CPU.
- 66 groups / 66,000 records / 3,829,019,863 raw message bytes.
- Subsequent raw compressed-span read: 934,883,767 bytes in 3.517 seconds wall,
  0.581 seconds CPU. This read follows the BGZF pass and can benefit from cache.
  Compressed span endpoints use virtual-offset compressed block starts; this
  probe does not include the trailing partial BGZF block of each range.

These are different passes, so do not subtract raw-read time from BGZF time to
claim an exact decompression cost. Process CPU establishes material execution
cost; wall minus CPU includes waiting and descheduling, not only disk latency.
The two jobs may run on different cluster nodes with different storage/cache
conditions. No measured threaded-decompression speedup is claimed.

## Proposed optimization order (not implemented)

1. Evaluate two-thread native BGZF decompression with bounded input/output
   buffers. Preserve GAI virtual offsets, GAM group boundaries, original record
   order and duplicate records. Threaded decompression must use an interface
   that actually supports BGZF virtual seeking; the current reader does not
   pass a thread setting. A whole-file ordinary gzip replacement is unsuitable
   for indexed random retrieval. Measure CPU/wall and memory before adopting.
2. Optimize node-postings construction, ideally keeping the record/node scan in
   native code and using compact record-index storage. Current Python traversal
   costs about 20 seconds during this initialization. Preserve repeated-visit
   deduplication per record, distinct duplicate records, and long-read coverage.
   Simply changing tuple to array does not establish a speedup; benchmark it.
3. Replace per-posting `sys.getsizeof(int)` calls with conservative arithmetic
   based on record-index bounds, keeping container/raw-byte accounting and the
   existing byte cap. This targets ~4 seconds here, so it is a secondary gain.
4. Keep the current group cache and contiguous processing. Consider bounded
   read-ahead only after measuring the new reader. Do not first copy the full
   176 GB GAM or increase cache limits without evidence of benefit.

This does not change the separately agreed 1,000-node candidate FIFO limit.
Do not change GAI query pruning as part of the first optimization: reducing
fetched ranges requires independent completeness tests, especially for long
reads spanning many nodes. Future acceptance should compare returned records
in exact order/multiplicity, then tensors/metadata, and separate startup from
steady-state timings over more than two batches.

## Reproducibility

- [Phase measurements](hg008_gam_reader_profile/profile.json)
- [BGZF measurements](hg008_gam_reader_profile/bgzf.json)
- [Phase diagnostic source](hg008_gam_reader_profile/profile_reader.py)
- [BGZF diagnostic source](hg008_gam_reader_profile/profile_bgzf.py)
