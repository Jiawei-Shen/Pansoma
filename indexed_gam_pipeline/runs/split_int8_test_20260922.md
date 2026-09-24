# Separate SNV and INDEL tensor generation

User-confirmed AF thresholds: SNV >=0.06, INDEL >=0.08. ALT support >=3.
512 node tasks, 32 concurrent builders, 8192 MiB GAM cache per builder.
Each task runs one shared SNV/INDEL builder, with outputs under SNV/task_NNNN
and INDEL/task_NNNN. Each batch fetches/decodes once for both types. The builder
exits before both validations and before the next task starts. Graph index and discovery are reused.
Full decoding remains selected, matching the previous production job.
Storage: int8, missing quality -1, channel 7 min(raw_count // 4,127).

Validation: 86 pipeline tests run, 84 passed, 2 native-tool tests skipped.
New synthetic end-to-end test checks SNP 3/50 passes at AF .06, DEL 3/50
fails at AF .08, DEL 4/50 passes, separated summaries/shards, and int8 count
encoding (331 ->82). Existing tests cover dynamic refill, cancellation,
partition coverage, storage boundaries and mixed-format merge rejection.

Shared-pass regression: 250 synthetic records fetch once/decode 250 times,
versus two independent type builds fetching twice/decoding 500 times. Separated
NPY files and variant summaries match those two independent builds byte for
byte, including 1-example shards. Cases cover DEL and INS threshold boundaries,
empty outputs, a combined max-tensors limit, and two candidate workers.
Production job 362291 retains its frozen two-pass code; this update does not
modify or restart that job.
