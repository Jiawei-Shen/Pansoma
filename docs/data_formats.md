# Data Formats

## Chromosome Node Filters

Text files with one integer node ID per non-empty line. Multiple files can be passed to `--chr_nodes`; the tensor generator unions them.

## Latest `.idx` Format

The current tensor generator expects a 4-byte little-endian block count followed by 30-byte entries:

```text
<I Q I I H I I>
node_id, offset, block_size, n_records, flags, max_read_len, max_cigar_len
```

## Candidate Node JSON

Supported structures:

```json
[
  {"node_id": "123", "sequence": "ACGT"}
]
```

or:

```json
{
  "nodes": [
    {"node_id": "123", "sequence": "ACGT"}
  ]
}
```

