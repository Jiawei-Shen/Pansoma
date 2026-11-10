"""R1.m2 INDEL examples: chr1 somatic truth INDEL tensors whose window spans several graph nodes.

Candidates: chr1 INDEL tensors of the twelve HPRC sets (HG008T / COLO829T x Illumina, PacBio HiFi or Fiber-seq, ONT x
HPRC v1.1 d9 / v2.1 d46) labelled 1 with reason representative_allele_is_somatic_truth (exact match to a truth INDEL),
truth FILTER PASS and inside the benchmark BED. Per tensor (read-only):
- node_len, start, end: the candidate's node and site offsets (chr1_variant_summary.ndjson, graph index sqlite);
- allele_nodes: nodes the allele itself occupies (1 + len(summary path));
- ref_nodes: GRCh38-path nodes overlapping the 101-column window [pos0 - 50, pos0 + span + 50) (grch38_path arrays);
- branch_cols: window columns where some read row sits on a node carried by fewer than all haplotypes (channel 6 code
  1..99, padding 0 excluded);
- linked_cols: of those, columns where the ALT (A1) and REF rows sit on different nodes in the same way the SNV example
  does: >= 80% of ALT rows share one channel-6 code, >= 80% of REF rows share another, outside the site columns.
Score = linked columns first, then allele_nodes, ref_nodes, branch_cols; ALT and REF rows >= 5 each.
usage: python select_examples.py OUT.tsv"""
import json
import sqlite3
import sys
from pathlib import Path

import numpy as np

D = Path("/scratch/jshen/data/pansoma_v2_tensors")
SETS = ["HG008T_PacBio", "HG008T_PacBio_d46", "HG008T_ONT", "HG008T_ONT_d46", "HG008T_Illumina", "HG008T_Illumina_d46",
        "COLO829T_fiberseq", "COLO829T_fiberseq_d46", "COLO829T_ONT", "COLO829T_ONT_d46", "COLO829T_Illumina",
        "COLO829T_Illumina_d46"]
GRCH38 = {90: Path("/scratch/jshen/data/pansoma_v2_tensors/graph_index/hprc-v1.1-mc-grch38.d9.grch38_path"),
          464: Path("/scratch/jshen/data/AF-Filtered_VG_Indexes_HPRC_V2/graph_index/hprc-v2.1-mc-grch38.d46.grch38_path")}


def chr1_path(h):
    p = GRCH38[h]
    c = next(c for c in json.load(open(p / "meta.json"))["contigs"] if c["name"] == "chr1")
    o, n = c["offset"], c["nodes"]
    return (np.load(p / "path_starts.npy", mmap_mode="r")[o:o + n],
            np.load(p / "path_nodes.npy", mmap_mode="r")[o:o + n])


def dominant(codes):
    if len(codes) == 0:
        return None, 0.0
    v, c = np.unique(codes, return_counts=True)
    i = int(np.argmax(c))
    return int(v[i]), c[i] / len(codes)


def main():
    out = open(sys.argv[1], "w")
    cols = ["set", "candidate_id", "shard_file", "index", "truth", "event", "len", "af", "alt_rows", "ref_rows",
            "node_len", "start", "end", "allele_nodes", "ref_nodes", "branch_cols", "linked_cols", "linked_detail"]
    print("\t".join(cols), file=out)
    for s in SETS:
        kd = D / s / "tensors" / "INDEL"
        man = json.load(open(kd / "manifest.json"))
        h = man["parameters"]["haplotypes"]
        starts, nodes = chr1_path(h)
        con = sqlite3.connect(f"file:{man['graph_index']['path']}?mode=ro", uri=True)
        lab = {}
        with open(kd / "chr1_labels.ndjson") as f:
            for line in f:
                if "representative_allele_is_somatic_truth" not in line:
                    continue
                r = json.loads(line)
                t = [x for x in r["somatic"] if x.get("representative")]
                if r["label"] == 1 and t and "PASS" in t[0]["filters"] and t[0]["in_bed"]:
                    lab[r["candidate_id"]] = (r, t[0])
        summ = {}
        with open(kd / "chr1_variant_summary.ndjson") as f:
            for line in f:
                cid = line[18:line.index('"', 18)]
                if cid in lab:
                    summ[cid] = json.loads(line)
        shards = {}
        n = 0
        for cid, (r, t) in lab.items():
            sm = summ[cid]
            sf = r["shard_file"]
            if sf not in shards:
                shards[sf] = np.load(kd / sf, mmap_mode="r")
            x = np.asarray(shards[sf][r["index_within_shard"]])
            g = {b["allele"]: (b["start_row"], b["end_row"]) for b in sm["row_groups"]}
            if "A1" not in g or "REF" not in g:
                continue
            a0, a1 = g["A1"]; r0, r1 = g["REF"]
            if a1 - a0 < 5 or r1 - r0 < 5:
                continue
            L = con.execute("select length(seq) from nodes where node_id=?", (sm["node_id"],)).fetchone()[0]
            ch6 = x[6].astype(int)
            site = set(range(sm["candidate_columns"][0], sm["candidate_columns"][1] + 1))
            branch = [c for c in range(ch6.shape[1]) if ((ch6[:, c] > 0) & (ch6[:, c] < 100)).any()]
            linked = []
            for c in branch:
                if c in site:
                    continue
                va, fa = dominant(ch6[a0:a1, c][ch6[a0:a1, c] > 0])
                vr, fr = dominant(ch6[r0:r1, c][ch6[r0:r1, c] > 0])
                if va is not None and vr is not None and va != vr and fa >= 0.8 and fr >= 0.8:
                    linked.append(f"{c}:{va}/{vr}")
            gr = r["grch38"]
            span = max(len(gr["ref"]), 1)
            lo, hi = gr["pos0"] - 50, gr["pos0"] + span + 50
            i0 = max(int(np.searchsorted(starts, lo, side="right")) - 1, 0)
            i1 = int(np.searchsorted(starts, hi, side="left"))
            row = [s, cid, sf, r["index_within_shard"],
                   f"{t['chrom']}:{t['vcf_pos']}:{t['vcf_ref']}>{t['vcf_alt']}", sm["event_type"], sm["event_length"],
                   f"{sm['af']:.3f}", a1 - a0, r1 - r0, L, sm["start"], sm["end"], 1 + len(sm["path"]), i1 - i0,
                   len(branch), len(linked), ",".join(linked)]
            print("\t".join(map(str, row)), file=out)
            n += 1
        out.flush()
        print(f"{s}: {len(lab)} exact PASS in-BED truth INDEL tensors, {n} with >= 5 ALT and REF rows",
              file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
