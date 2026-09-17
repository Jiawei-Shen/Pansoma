"""GAM mappings to unpadded segments, then node-forward pileup records."""
from collections import defaultdict


def raw_segments(alignment, wanted, min_mapq=10):
    # Preserve the original NPU writer's >10 gate plus tensor's >= threshold.
    if alignment.mapping_quality <= 10 or alignment.mapping_quality < min_mapq:
        return
    read_offset = 0
    for mapping in alignment.path.mapping:
        length = sum(e.to_length for e in mapping.edit)
        start = read_offset
        read_offset += length  # Includes mappings to non-target nodes.
        nid = mapping.position.node_id
        if nid not in wanted:
            continue
        ops = []
        for edit in mapping.edit:
            f, t = edit.from_length, edit.to_length
            if f == t:
                ops.append((f, "X" if edit.sequence else "M"))
            elif f and not t:
                ops.append((f, "D"))
            elif t and not f:
                ops.append((t, "I"))
            else:
                # The old writer silently omits these edits. Do not make a
                # malformed CIGAR or silently change their meaning here.
                raise ValueError(f"Unsupported complex edit {f}->{t} at node {nid}")
        seq = alignment.sequence[start:read_offset].upper()
        qual = list(alignment.quality[start:read_offset])
        if len(seq) != length or len(qual) != length:
            raise ValueError(f"Sequence/quality length mismatch at node {nid}")
        if not ops:
            continue
        yield nid, dict(offset_on_node=int(mapping.position.offset),
                        read_sequence=seq, processed_quality_values=qual,
                        cigar_ops=ops,
                        original_cigar_str="".join(f"{n}{op}" for n, op in ops),
                        strand="-" if mapping.position.is_reverse else "+",
                        mapping_quality=int(alignment.mapping_quality))


def orient(segment, node_length):
    result = dict(segment)
    span = sum(n for n, op in segment["cigar_ops"] if op in "MX=D")
    offset = segment["offset_on_node"]
    if offset < 0 or offset + span > node_length:
        raise ValueError("Mapping exceeds node sequence; check graph/GAM compatibility")
    if segment["strand"] == "-":
        result["read_sequence"] = segment["read_sequence"].translate(
            str.maketrans("ACGTN", "TGCAN"))[::-1]
        result["processed_quality_values"] = segment["processed_quality_values"][::-1]
        result["cigar_ops"] = segment["cigar_ops"][::-1]
        result["offset_on_node"] = node_length - span - offset
    return result


def collect(alignments, nodes, min_mapq=10):
    result = defaultdict(list)
    for alignment in alignments:
        for nid, segment in raw_segments(alignment, nodes, min_mapq):
            result[nid].append(segment)
    return result
