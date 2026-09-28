"""Small helpers shared across the package: atomic JSON, file stamps, format names, node lists, batching."""
import hashlib
import json
from pathlib import Path

# Files made before this package name the same formats with a "-v1" suffix (the HPRC graph index, the
# reference-path directory, merged tensor sets); readers accept both, writers use the plain names.
EARLIER_FORMAT_NAMES = {"gbz-all-nodes-v1": "gbz-all-nodes", "distinct-gbwt-paths-v1": "distinct-gbwt-paths",
                        "gfa-reference-path-v1": "gfa-reference-path", "chromosome-shards-v1": "chromosome-shards"}


def format_name(value):
    """A format name as this package writes it (see EARLIER_FORMAT_NAMES); a non-string value is returned as is,
    so the reader refuses it like any other unknown name."""
    return EARLIER_FORMAT_NAMES.get(value, value) if isinstance(value, str) else value


def write_json(path, value):
    """Write atomically so a reader never observes a half-written file."""
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def read_json(path):
    return json.loads(Path(path).read_text())


def stamp(path):
    """Cheap identity of a large input: resolved path, size and mtime."""
    path = Path(path).resolve()
    s = path.stat()
    return dict(path=str(path), size=s.st_size, mtime_ns=s.st_mtime_ns)


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def new_output(path):
    """Create an output directory; refuse to reuse one that already has content."""
    path = Path(path)
    if path.exists() and any(path.iterdir()):
        raise ValueError(f"Output directory must be empty: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_nodes(path):
    """Sorted, unique, positive node IDs from a one-per-line text file."""
    result = sorted({int(line) for line in Path(path).read_text().split()})
    if not result or result[0] <= 0:
        raise ValueError("Node list must contain positive node IDs")
    return result


def batches(nodes, size, max_span):
    """Contiguous batches bounded by node count and by node-ID span."""
    batch = []
    for nid in sorted(nodes):
        if batch and (len(batch) >= size or nid - batch[0] > max_span):
            yield batch
            batch = []
        batch.append(nid)
    if batch:
        yield batch

