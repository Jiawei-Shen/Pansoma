"""(HG008-T ONT-UL, v3 tensors) SNV misses in detail: are the reads a simple base change at this position, or a different local haplotype?

For every SNV no_candidate locus (and SNV controls), over MAPQ>10 reads anchored on GRCh38 20 bp
outside the site (nearest matched GRCh38 base at <= pos-20 and >= pos+20, up to 2 kb away):
  same_length   the read's local sequence has the REF length (no indel between the anchors)
  base          for same-length reads, the read base at the SNV position (REF / ALT / other)
  via           the node the read uses at that position: GRCh38 node or not
Output: snv_detail.jsonl (one record per locus).
"""
import json, os, sys
from bisect import bisect_left, bisect_right
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reads as R  # noqa: E402

HERE = Path(__file__).resolve().parent
FLANK, FAR = 20, 2000


def analyse(chunk):
    path = R.ReferencePath(R.REFPATH)
    fasta = pysam.FastaFile(R.FASTA)
    reader = R.IndexedGam(R.GAM, cache_bytes=512 << 20)
    visits, chrom_of, start0, rev_of, lengths = (np.asarray(path.visits), np.asarray(path.chrom), np.asarray(path.start0),
                                                 np.asarray(path.reverse), np.asarray(path.lengths))
    out, genome_chrom = [], None
    with R.GraphIndex(R.GRAPH) as graph:
        for m in chunk:
            if m["chrom"] != genome_chrom:
                genome_chrom, genome = m["chrom"], fasta.fetch(m["chrom"]).upper()
                contig = path.contig_index[m["chrom"]]
                part = path.contig_slice(m["chrom"])
                wn, ws = np.asarray(path.path_nodes[part]), np.asarray(path.path_starts[part])
            p = m["pos0"]
            left, right = p - FLANK, p + FLANK
            a = max(0, int(np.searchsorted(ws, left - 300, side="right")) - 1)
            b = int(np.searchsorted(ws, right + 300, side="right"))
            nodes = {int(n) for n in wn[a:b]}
            info = dict(truth_id=m["truth_id"], status=m["status"], chrom=m["chrom"], vcf_pos=int(m["vcf_pos"]),
                        ref=m["ref"], alt=m["alt"], records=0, low_mapq=0, spanning=0, same_length=0,
                        indel_length_hist=Counter(), base=Counter(), via=Counter(), alt_via=Counter())
            records = list(reader.fetch(nodes, {}))
            info["records"] = len(records)
            good = [x for x in records if x.mapping_quality > R.MIN_MAPQ]
            info["low_mapq"] = len(records) - len(good)
            subs = [s for s in (R.trim(x, nodes) for x in good) if s]
            seqs = {n: r["sequence"] for n, r in graph.get_nodes({mp.position.node_id for s in subs for mp in s.path.mapping}).items()} if subs else {}
            for s in subs:
                read, _ = R.dec.decode_alignment(s, seqs, 50, None)
                cols = read.columns
                lin = []
                for c in cols:
                    on = visits[c.node] == 1 and chrom_of[c.node] == contig and not c.boundary
                    lin.append(int(start0[c.node] + (lengths[c.node] - 1 - c.pos if rev_of[c.node] else c.pos)) if on else None)
                matched = {x: i for i, x in enumerate(lin) if x is not None and cols[i].op == "M" and cols[i].read == cols[i].ref}
                keys = sorted(matched)
                k, j = bisect_right(keys, left) - 1, bisect_left(keys, right)
                if k < 0 or j >= len(keys) or left - keys[k] > FAR or keys[j] - right > FAR:
                    continue
                pl, pr = keys[k], keys[j]
                il, ir = matched[pl], matched[pr]
                flip = cols[il].reverse != bool(rev_of[cols[il].node])
                i0, i1 = sorted((il, ir))
                between = cols[i0 + 1:i1]
                seq = "".join(c.read for c in between if c.read != "-")
                if flip:
                    seq = R.rc(seq)
                info["spanning"] += 1
                diff = len(seq) - (pr - pl - 1)
                info["indel_length_hist"][diff if abs(diff) <= 5 else (">5" if diff > 0 else "<-5")] += 1
                # the node the read uses at the SNV position: a GRCh38 column with linear position p, or not
                at = [i for i in range(i0 + 1, i1) if lin[i] == p]
                via = "grch38_node" if at else "other_node_or_skip"
                info["via"][via] += 1
                if diff != 0:
                    continue
                info["same_length"] += 1
                base = seq[p - pl - 1]
                kind = "REF" if base == m["ref"] else "ALT" if base == m["alt"] else "other"
                info["base"][kind] += 1
                if kind == "ALT":
                    info["alt_via"][via] += 1
            for key in ("indel_length_hist", "base", "via", "alt_via"):
                info[key] = {str(k): v for k, v in info[key].items()}
            out.append(info)
    return out


def main():
    rows = [r for r in R.loci() if r["kind"] == "SNP"]
    rows.sort(key=lambda m: (m["chrom"], m["pos0"]))
    n = int(os.environ.get("SLURM_CPUS_PER_TASK", 8))
    chunks = [rows[i:i + 5] for i in range(0, len(rows), 5)]
    with ProcessPoolExecutor(n) as pool, open(HERE / "snv_detail.jsonl", "w") as f:
        for part in pool.map(analyse, chunks):
            for r in part:
                f.write(json.dumps(r) + "\n")
            f.flush()


if __name__ == "__main__":
    main()
