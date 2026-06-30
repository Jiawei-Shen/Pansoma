#!/usr/bin/env python3

import csv
import json
import subprocess
import sys
from pathlib import Path


XG = "/scratch/jshen/data/HG008_GIAB/indexes_HG008N_curatedv6_250714_polished6.2_rebuild_hprc_d9_vg1660/hprc-v1.1-mc-grch38.d9.plus_HG008N_curatedv6_250714_polished6_2.xg"
INPUT_TSV = "/scratch/jshen/data/HG008_GIAB/changed_repr_records/changed_representation_variants.tsv"
OUTPUT_TSV = "/scratch/jshen/data/HG008_GIAB/changed_repr_records/changed_representation_variants.with_nodes.tsv"


def run_cmd(cmd):
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed:\n{' '.join(cmd)}\n\nSTDERR:\n{result.stderr}"
        )
    return result.stdout


def extract_node_ids_from_vg_find(xg_path, path_name, start_pos, end_pos):
    """
    Query vg find on a path interval and return sorted unique node IDs.

    start_pos/end_pos are passed directly into vg find as:
        path:start-end
    """
    region = f"{path_name}:{start_pos}-{end_pos}"

    # vg find -> protobuf graph; vg view -j -> JSON
    p1 = subprocess.Popen(
        ["vg", "find", "-x", xg_path, "-p", region],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=False,
    )
    p2 = subprocess.Popen(
        ["vg", "view", "-j", "-"],
        stdin=p1.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    p1.stdout.close()

    out2, err2 = p2.communicate()
    _, err1 = p1.communicate()

    if p1.returncode != 0:
        raise RuntimeError(f"vg find failed for region {region}\nSTDERR:\n{err1.decode() if isinstance(err1, bytes) else err1}")
    if p2.returncode != 0:
        raise RuntimeError(f"vg view failed for region {region}\nSTDERR:\n{err2}")

    if not out2.strip():
        return []

    graph_json = json.loads(out2)
    nodes = graph_json.get("node", [])
    node_ids = sorted({int(n["id"]) for n in nodes if "id" in n})
    return node_ids


def choose_query_interval(pos, ref, alt, grch38_type):
    """
    Decide which GRCh38 interval to query.

    Rules:
    - SNV/MNV: span the REF allele
    - DEL: span the deleted reference bases
    - INS: use the anchor base at POS
    - fallback: span REF length if possible, else single base
    """
    pos = int(pos)
    ref_len = len(ref) if ref else 1

    t = (grch38_type or "").upper()

    if t in {"SNV", "SNP"}:
        return pos, pos
    elif t == "MNV":
        return pos, pos + ref_len - 1
    elif t == "DEL":
        return pos, pos + ref_len - 1
    elif t == "INS":
        return pos, pos
    else:
        return pos, pos + ref_len - 1


def main():
    input_path = Path(INPUT_TSV)
    output_path = Path(OUTPUT_TSV)

    if not input_path.exists():
        print(f"ERROR: input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)
    if not Path(XG).exists():
        print(f"ERROR: xg file not found: {XG}", file=sys.stderr)
        sys.exit(1)

    with input_path.open() as fin, output_path.open("w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        fieldnames = reader.fieldnames + [
            "PATH_NAME",
            "QUERY_START",
            "QUERY_END",
            "NODE_COUNT",
            "NODE_IDS",
        ]
        writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for i, row in enumerate(reader, 1):
            chrom = row["#CHROM"]
            pos = int(row["POS"])
            ref = row["REF"]
            alt = row["ALT"]
            grch38_type = row.get("GRCh38_TYPE", "")

            path_name = f"GRCh38#0#{chrom}"
            qstart, qend = choose_query_interval(pos, ref, alt, grch38_type)

            try:
                node_ids = extract_node_ids_from_vg_find(
                    XG, path_name, qstart, qend
                )
            except Exception as e:
                print(f"[WARN] line {i}: {chrom}:{pos} failed: {e}", file=sys.stderr)
                node_ids = []

            out_row = dict(row)
            out_row["PATH_NAME"] = path_name
            out_row["QUERY_START"] = qstart
            out_row["QUERY_END"] = qend
            out_row["NODE_COUNT"] = len(node_ids)
            out_row["NODE_IDS"] = ",".join(map(str, node_ids))
            writer.writerow(out_row)

            if i % 100 == 0:
                print(f"[INFO] processed {i} variants", file=sys.stderr)

    print(f"[DONE] wrote: {output_path}")


if __name__ == "__main__":
    main()