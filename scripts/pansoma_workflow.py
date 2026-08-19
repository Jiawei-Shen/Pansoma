#!/usr/bin/env python3
"""Container-oriented end-to-end Pansoma training and inference workflows."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Iterator, Optional, Sequence, TextIO

import numpy as np
import pysam


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ML_ROOT = PROJECT_ROOT / "machine_learning" / "pansoma_net"
PYTHON = sys.executable


def log(message: str) -> None:
    print(f"[pansoma] {message}", flush=True)


def run(
    command: Sequence[str],
    *,
    cwd: Path = PROJECT_ROOT,
    env: Optional[dict] = None,
    stdout_path: Optional[Path] = None,
) -> None:
    rendered = shlex.join(str(part) for part in command)
    log(f"run: {rendered}")
    if stdout_path is None:
        subprocess.run(
            [str(part) for part in command],
            cwd=str(cwd),
            env=env,
            check=True,
        )
        return

    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = stdout_path.with_suffix(stdout_path.suffix + ".tmp")
    temporary.unlink(missing_ok=True)
    try:
        with temporary.open("wb") as output:
            subprocess.run(
                [str(part) for part in command],
                cwd=str(cwd),
                env=env,
                stdout=output,
                check=True,
            )
        temporary.replace(stdout_path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def require_file(path_text: str, label: str) -> Path:
    path = Path(path_text).expanduser().resolve()
    if not path.is_file():
        raise SystemExit(f"ERROR: {label} not found: {path}")
    return path


def require_directory(path_text: str, label: str) -> Path:
    path = Path(path_text).expanduser().resolve()
    if not path.is_dir():
        raise SystemExit(f"ERROR: {label} not found: {path}")
    return path


def is_complete(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def file_signature(path_text: Optional[str]) -> Optional[dict[str, object]]:
    if not path_text:
        return None
    path = Path(path_text).expanduser().resolve()
    signature: dict[str, object] = {"path": str(path)}
    if path.is_file():
        stat = path.stat()
        signature.update({"size": stat.st_size, "mtime_ns": stat.st_mtime_ns})
    else:
        signature["missing"] = True
    return signature


def inference_configuration(args: argparse.Namespace, checkpoint: Path) -> dict[str, object]:
    """Return the result-affecting inputs used to validate inference caches."""
    return {
        "schema_version": 1,
        "sample": args.sample,
        "checkpoint": file_signature(str(checkpoint)),
        "fastq1": file_signature(args.fastq1),
        "fastq2": file_signature(args.fastq2),
        "gbz": file_signature(args.gbz),
        "min_index": file_signature(args.min_index),
        "dist_index": file_signature(args.dist_index),
        "zipcode_index": file_signature(args.zipcode_index),
        "gfa": file_signature(args.gfa),
        "chromosomes": list(args.chromosomes),
        "mapper_preset": args.mapper_preset,
        "variant_type": args.variant_type,
        "min_af": args.min_af,
        "min_variants": args.min_variants,
        "min_allele_bq": args.min_allele_bq,
        "min_mapq": args.min_mapq,
        "max_indel_len": args.max_indel_len,
        "shard_size": args.shard_size,
        "batch_size": args.batch_size,
        "device": args.device,
        "min_true_prob": args.min_true_prob,
        "min_true_prob_no_anchor": args.min_true_prob_no_anchor,
    }


def indexed_vcf_matches_manifest(
    vcf: Path,
    manifest: Path,
    expected: dict[str, object],
) -> bool:
    index_ready = is_complete(Path(str(vcf) + ".tbi")) or is_complete(
        Path(str(vcf) + ".csi")
    )
    if not is_complete(vcf) or not index_ready or not is_complete(manifest):
        return False
    try:
        with manifest.open("r", encoding="utf-8") as handle:
            observed = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return False
    return observed == expected


def clear_indexed_vcf(vcf: Path, manifest: Path) -> None:
    for path in (vcf, Path(str(vcf) + ".tbi"), Path(str(vcf) + ".csi"), manifest):
        path.unlink(missing_ok=True)


def run_if_missing(
    name: str,
    command: Sequence[str],
    outputs: Iterable[Path],
    *,
    cwd: Path = PROJECT_ROOT,
    env: Optional[dict] = None,
    stdout_path: Optional[Path] = None,
) -> None:
    expected = list(outputs)
    if expected and all(is_complete(path) for path in expected):
        log(f"reuse {name}: {', '.join(str(path) for path in expected)}")
        return
    run(command, cwd=cwd, env=env, stdout_path=stdout_path)
    missing = [str(path) for path in expected if not is_complete(path)]
    if missing:
        raise SystemExit(f"ERROR: {name} finished without expected output(s): {', '.join(missing)}")


def atomic_json_dump(data: object, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, ensure_ascii=False)
    temporary.replace(output)


def open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def graph_fingerprint(paths: Sequence[Path]) -> str:
    digest = hashlib.sha256()
    for path in paths:
        stat = path.stat()
        digest.update(str(path).encode("utf-8"))
        digest.update(str(stat.st_size).encode("ascii"))
        digest.update(str(stat.st_mtime_ns).encode("ascii"))
    return digest.hexdigest()[:16]


def stage_truth_vcf(vcf: Path, output_dir: Path) -> Path:
    has_index = Path(str(vcf) + ".tbi").is_file() or Path(str(vcf) + ".csi").is_file()
    if vcf.name.endswith(".vcf.gz") and has_index:
        return vcf
    if not (vcf.name.endswith(".vcf") or vcf.name.endswith(".vcf.gz")):
        raise SystemExit("ERROR: --truth-vcf must end with .vcf or .vcf.gz")

    output_dir.mkdir(parents=True, exist_ok=True)
    truth_id = graph_fingerprint([vcf])
    staged = output_dir / f"truth.{truth_id}.vcf.gz"
    index = Path(str(staged) + ".tbi")
    if is_complete(staged) and is_complete(index):
        log(f"reuse indexed truth VCF: {staged}")
        return staged

    plain = output_dir / f"truth.{truth_id}.vcf"
    log(f"stage and index truth VCF: {vcf}")
    if vcf.name.endswith(".gz"):
        with gzip.open(vcf, "rb") as source, plain.open("wb") as target:
            shutil.copyfileobj(source, target)
    else:
        shutil.copyfile(vcf, plain)

    try:
        pysam.tabix_compress(str(plain), str(staged), force=True)
        pysam.tabix_index(str(staged), preset="vcf", force=True)
    finally:
        plain.unlink(missing_ok=True)

    if not is_complete(staged) or not is_complete(index):
        raise SystemExit(f"ERROR: failed to create indexed truth VCF: {staged}")
    return staged


def validate_truth_contigs(vcf: Path, chromosomes: Sequence[str]) -> None:
    with pysam.VariantFile(str(vcf)) as handle:
        available = set(handle.header.contigs)
    missing = set(chromosomes).difference(available)
    if missing:
        raise SystemExit(
            f"ERROR: truth VCF is missing requested contigs: {sorted(missing)}. "
            "VCF and graph chromosome names must match."
        )


def labels_match_truth(tensor_dir: Path, truth_id: str) -> bool:
    manifest = tensor_dir / "labels.manifest.json"
    labels = list(tensor_dir.glob("*_labels.npy"))
    if not is_complete(manifest) or not labels:
        return False
    try:
        with manifest.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return False
    return data.get("truth_fingerprint") == truth_id


def clear_stale_labels(tensor_dir: Path) -> None:
    for pattern in ("*_labels.npy", "variant_summary_classified.ndjson", "labels.manifest.json"):
        for path in tensor_dir.glob(pattern):
            path.unlink()


def chromosome_filters(base: Path, chromosome: str) -> tuple[Path, Path]:
    chrom_dir = base / chromosome
    component = chrom_dir / f"{chromosome}.component.nodes.raw.txt"
    grch38 = chrom_dir / f"{chromosome}.GRCh38_path.nodes.raw.txt"
    for path in (component, grch38):
        if not is_complete(path):
            raise SystemExit(f"ERROR: chromosome node filter not found or empty: {path}")
    return component, grch38


def common_paths(args: argparse.Namespace) -> dict[str, Path]:
    work = Path(args.work_dir).expanduser().resolve()
    paths = {
        "work": work,
        "align": work / "alignment",
        "store": work / "alignment_store",
        "nodes": work / "nodes",
        "truth": work / "truth",
        "tensors": work / "tensors",
        "models": work / "models",
        "results": work / "results",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


LONG_READ_PRESETS = {"hifi", "r10"}


def resolve_mapper_preset(args: argparse.Namespace) -> str:
    preset = args.mapper_preset
    if preset in LONG_READ_PRESETS and args.fastq2:
        raise SystemExit(
            f"ERROR: vg giraffe preset {preset!r} accepts one FASTQ; omit --fastq2"
        )
    return preset


WALK_TOKEN = re.compile(r"([<>])([^<>]+)")


def parse_w_walk(walk: str) -> Iterator[tuple[str, str]]:
    for orientation, node_id in WALK_TOKEN.findall(walk):
        yield node_id, "+" if orientation == ">" else "-"


def parse_p_walk(walk: str) -> Iterator[tuple[str, str]]:
    for token in walk.split(","):
        if token and token[-1] in "+-":
            yield token[:-1], token[-1]


def segment_length(fields: list[str]) -> int:
    sequence = fields[2] if len(fields) > 2 else "*"
    if sequence != "*":
        return len(sequence)
    for field in fields[3:]:
        if field.startswith("LN:i:"):
            return int(field[5:])
    raise ValueError(f"GFA segment {fields[1]} has no sequence or LN:i length")


def build_coordinate_tables(gfa: Path, output_dir: Path, chromosomes: Sequence[str]) -> None:
    expected = [output_dir / f"{chrom}.GRCh38.nodes.tsv" for chrom in chromosomes]
    if all(is_complete(path) for path in expected):
        log(f"reuse GRCh38 coordinate tables: {output_dir}")
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    targets = set(chromosomes)
    walk_paths = {chrom: output_dir / f".{chrom}.walk.tmp" for chrom in chromosomes}
    starts: dict[str, int] = {}
    found_w: set[str] = set()
    found_p: set[str] = set()
    lengths: dict[str, int] = {}

    for path in walk_paths.values():
        path.unlink(missing_ok=True)

    log(f"scan GFA for segment lengths and GRCh38 paths: {gfa}")
    try:
        with open_text(gfa) as handle:
            for line in handle:
                if line.startswith("S\t"):
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) >= 3:
                        lengths[fields[1]] = segment_length(fields)
                    continue

                if line.startswith("W\t"):
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 7 or fields[1] != "GRCh38" or fields[2] != "0":
                        continue
                    chromosome = fields[3]
                    if chromosome not in targets or chromosome in found_w:
                        continue
                    starts[chromosome] = int(fields[4]) + 1
                    with walk_paths[chromosome].open("w", encoding="utf-8") as walk_handle:
                        for node_id, strand in parse_w_walk(fields[6]):
                            walk_handle.write(f"{node_id}\t{strand}\n")
                    found_w.add(chromosome)
                    continue

                if line.startswith("P\t"):
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 3:
                        continue
                    for chromosome in targets.difference(found_w).difference(found_p):
                        if fields[1] != f"GRCh38#0#{chromosome}":
                            continue
                        starts[chromosome] = 1
                        with walk_paths[chromosome].open("w", encoding="utf-8") as walk_handle:
                            for node_id, strand in parse_p_walk(fields[2]):
                                walk_handle.write(f"{node_id}\t{strand}\n")
                        found_p.add(chromosome)
                        break

        missing_paths = targets.difference(found_w.union(found_p))
        if missing_paths:
            raise SystemExit(f"ERROR: GRCh38 paths not found in GFA: {sorted(missing_paths)}")

        for chromosome, output in zip(chromosomes, expected):
            temporary = output.with_suffix(output.suffix + ".tmp")
            current_position = starts[chromosome]
            count = 0
            with walk_paths[chromosome].open("r", encoding="utf-8") as walk_handle, temporary.open(
                "w", encoding="utf-8"
            ) as target:
                target.write("node_id\tgrch38_position_start\tstrand_in_path\tlength\tchrom\n")
                for line in walk_handle:
                    node_id, strand = line.rstrip("\n").split("\t")
                    length = lengths.get(node_id)
                    if length is None:
                        raise SystemExit(
                            f"ERROR: GRCh38 path node {node_id} has no segment length in GFA"
                        )
                    target.write(
                        f"{node_id}\t{current_position}\t{strand}\t{length}\t{chromosome}\n"
                    )
                    current_position += length
                    count += 1
            temporary.replace(output)
            log(f"GRCh38 coordinate table {chromosome}: {count:,} nodes")
    finally:
        for path in walk_paths.values():
            path.unlink(missing_ok=True)


def build_graph_resources(
    args: argparse.Namespace,
    paths: dict[str, Path],
    gbz: Path,
    gfa: Path,
) -> tuple[Path, Path, str]:
    fingerprint = graph_fingerprint([gbz, gfa])
    cache_root = (
        Path(args.resource_cache).expanduser().resolve()
        if args.resource_cache
        else paths["work"] / "graph_resources"
    )
    graph_root = cache_root / fingerprint
    filter_base = graph_root / "chromosome_filters"
    coordinate_dir = graph_root / "coordinates"

    filters_ready = all(
        is_complete(filter_base / chrom / f"{chrom}.component.nodes.raw.txt")
        and is_complete(filter_base / chrom / f"{chrom}.GRCh38_path.nodes.raw.txt")
        for chrom in args.chromosomes
    )
    if filters_ready:
        log(f"reuse chromosome node filters: {filter_base}")
    else:
        env = os.environ.copy()
        env.update(
            {
                "GBZ": str(gbz),
                "GFA": str(gfa),
                "OUTDIR": str(filter_base),
                "CHROMOSOMES": " ".join(args.chromosomes),
            }
        )
        run(["bash", "scripts/build_chr_node_filters.sh"], env=env)

    for chromosome in args.chromosomes:
        chromosome_filters(filter_base, chromosome)
    build_coordinate_tables(gfa, coordinate_dir, args.chromosomes)

    manifest = graph_root / "manifest.json"
    atomic_json_dump(
        {
            "fingerprint": fingerprint,
            "gbz": str(gbz),
            "gfa": str(gfa),
            "chromosomes": list(args.chromosomes),
        },
        manifest,
    )
    return filter_base, coordinate_dir, fingerprint


def enrich_candidate_node_map(
    raw_json: Path,
    output_json: Path,
    coordinate_dir: Path,
    chromosomes: Sequence[str],
) -> None:
    if is_complete(output_json):
        log(f"reuse coordinate-aware candidate node map: {output_json}")
        return

    with raw_json.open("r", encoding="utf-8") as handle:
        records = json.load(handle)
    if not isinstance(records, list):
        raise SystemExit(f"ERROR: candidate node JSON must contain a list: {raw_json}")

    records_by_id = {
        str(record["node_id"]): record
        for record in records
        if isinstance(record, dict) and record.get("node_id") is not None
    }
    anchored = 0
    ambiguous: set[str] = set()
    for chromosome in chromosomes:
        table = coordinate_dir / f"{chromosome}.GRCh38.nodes.tsv"
        with table.open("r", encoding="utf-8") as handle:
            next(handle, None)
            for line in handle:
                node_id, position, strand, length, chrom = line.rstrip("\n").split("\t")
                record = records_by_id.get(node_id)
                if record is None:
                    continue
                if "grch38_position_start" in record:
                    ambiguous.add(node_id)
                    continue
                record.update(
                    {
                        "grch38_position_start": int(position),
                        "strand_in_path": strand,
                        "length": int(length),
                        "chrom": chrom,
                    }
                )
                anchored += 1

    for node_id in ambiguous:
        record = records_by_id[node_id]
        for key in ("grch38_position_start", "strand_in_path", "chrom"):
            record.pop(key, None)
    anchored -= len(ambiguous)
    atomic_json_dump(records, output_json)
    log(
        f"candidate node map: total={len(records):,}, GRCh38 anchors={anchored:,}, "
        f"ambiguous={len(ambiguous):,}"
    )


def prepare_sample(args: argparse.Namespace) -> tuple[dict[str, Path], Path, dict[str, Path]]:
    fastq1 = require_file(args.fastq1, "FASTQ 1")
    fastq2 = require_file(args.fastq2, "FASTQ 2") if args.fastq2 else None
    mapper_preset = resolve_mapper_preset(args)
    gbz = require_file(args.gbz, "GBZ graph")
    min_index = require_file(args.min_index, "minimizer index")
    dist_index = require_file(args.dist_index, "distance index")
    zipcode_index = (
        require_file(args.zipcode_index, "zipcode index") if args.zipcode_index else None
    )
    if mapper_preset in LONG_READ_PRESETS and zipcode_index is None:
        raise SystemExit(
            "ERROR: PacBio HiFi and ONT R10 require --zipcode-index together "
            "with a long-read minimizer index"
        )
    gfa = require_file(args.gfa, "GFA graph")
    paths = common_paths(args)
    filter_base, coordinate_dir, graph_id = build_graph_resources(args, paths, gbz, gfa)
    gam = paths["align"] / f"{args.sample}.gam"
    stats = paths["store"] / f"{args.sample}.unperfect_nodes.pkl"
    store_prefix = paths["store"] / args.sample
    dat = store_prefix.with_suffix(".dat")
    idx = store_prefix.with_suffix(".idx")
    raw_candidate_nodes = paths["nodes"] / f"{args.sample}.{graph_id}.candidate_nodes.raw.json"
    candidate_nodes = paths["nodes"] / f"{args.sample}.{graph_id}.candidate_nodes.json"

    align_command = [
        "vg",
        "giraffe",
        "-Z",
        str(gbz),
        "-m",
        str(min_index),
        "-d",
        str(dist_index),
        "-f",
        str(fastq1),
        "-b",
        mapper_preset,
        "-t",
        str(args.threads),
        "-p",
        "-o",
        "gam",
    ]
    if fastq2:
        align_command.extend(["-f", str(fastq2)])
    if zipcode_index:
        align_command.extend(["-z", str(zipcode_index)])

    run_if_missing(
        "FASTQ alignment",
        align_command,
        [gam],
        stdout_path=gam,
    )
    run_if_missing(
        "imperfect-node scan",
        [
            PYTHON,
            "-u",
            "scripts/find_unperfect_nodes.py",
            str(gam),
            "--output",
            str(stats),
            "--output_format",
            "pickle",
            "--milestone",
            str(args.scan_milestone),
            "--threads",
            str(args.threads),
        ],
        [stats],
    )
    run_if_missing(
        "alignment store",
        [
            PYTHON,
            "-u",
            "scripts/build_dat_idx.py",
            str(gam),
            str(stats),
            str(store_prefix),
            "--milestone",
            str(args.store_milestone),
            "--threads",
            str(args.threads),
        ],
        [dat, idx],
    )
    run_if_missing(
        "raw candidate node map",
        [
            PYTHON,
            "-u",
            "scripts/build_node_json.py",
            "--gfa",
            str(gfa),
            "--idx",
            str(idx),
            "--out",
            str(raw_candidate_nodes),
        ],
        [raw_candidate_nodes],
    )
    enrich_candidate_node_map(
        raw_candidate_nodes,
        candidate_nodes,
        coordinate_dir,
        args.chromosomes,
    )

    tensor_dirs: dict[str, Path] = {}
    for chromosome in args.chromosomes:
        component, grch38 = chromosome_filters(filter_base, chromosome)
        tensor_dir = paths["tensors"] / chromosome
        tensor_dir.mkdir(parents=True, exist_ok=True)
        summary = tensor_dir / "variant_summary.ndjson"
        data_files = sorted(tensor_dir.glob("*_data.npy"))
        if is_complete(summary) and data_files:
            log(f"reuse tensors for {chromosome}: {tensor_dir}")
        else:
            run(
                [
                    PYTHON,
                    "-u",
                    "scripts/generate_testing_tensors.py",
                    str(dat),
                    str(idx),
                    str(tensor_dir),
                    str(candidate_nodes),
                    "--chr_nodes",
                    str(component),
                    str(grch38),
                    "--num_workers",
                    str(args.tensor_workers),
                    "--variant_type",
                    args.variant_type,
                    "--min_af",
                    str(args.min_af),
                    "--min_variants",
                    str(args.min_variants),
                    "--min_allele_bq",
                    str(args.min_allele_bq),
                    "--min_mapq",
                    str(args.min_mapq),
                    "--max_indel_len",
                    str(args.max_indel_len),
                    "--shard_size",
                    str(args.shard_size),
                ]
            )
            if not is_complete(summary) or not list(tensor_dir.glob("*_data.npy")):
                raise SystemExit(f"ERROR: tensor generation produced no usable shards for {chromosome}")
        tensor_dirs[chromosome] = tensor_dir

    return paths, candidate_nodes, tensor_dirs


def link_or_filter_shard(
    data_path: Path,
    labels_path: Path,
    out_data: Path,
    out_labels: Path,
) -> int:
    out_data.parent.mkdir(parents=True, exist_ok=True)
    if is_complete(out_data) and is_complete(out_labels):
        labels = np.load(out_labels, mmap_mode="r")
        return int(labels.shape[0])

    data = np.load(data_path, mmap_mode="r")
    labels = np.load(labels_path, mmap_mode="r")
    if data.shape[0] != labels.shape[0]:
        raise SystemExit(
            f"ERROR: shard/label length mismatch: {data_path}={data.shape[0]}, "
            f"{labels_path}={labels.shape[0]}"
        )
    keep = np.isin(labels, (0, 1))
    kept = int(keep.sum())
    if kept == 0:
        log(f"skip shard with no known labels: {data_path}")
        return 0

    if kept == int(labels.shape[0]):
        for source, target in ((data_path, out_data), (labels_path, out_labels)):
            if target.is_symlink() and not target.exists():
                target.unlink()
            if not target.exists():
                target.symlink_to(source)
    else:
        indices = np.flatnonzero(keep)
        filtered_data = np.lib.format.open_memmap(
            out_data,
            mode="w+",
            dtype=data.dtype,
            shape=(kept, *data.shape[1:]),
        )
        filtered_labels = np.lib.format.open_memmap(
            out_labels,
            mode="w+",
            dtype=np.int8,
            shape=(kept,),
        )
        chunk_size = 256
        for start in range(0, kept, chunk_size):
            end = min(start + chunk_size, kept)
            source_indices = indices[start:end]
            filtered_data[start:end] = data[source_indices]
            filtered_labels[start:end] = labels[source_indices]
        filtered_data.flush()
        filtered_labels.flush()
        del filtered_data, filtered_labels
    return kept


def assemble_training_dataset(
    tensor_dirs: dict[str, Path],
    dataset_root: Path,
    validation_chromosomes: set[str],
) -> None:
    total = {"train": 0, "val": 0}
    for chromosome, tensor_dir in tensor_dirs.items():
        split = "val" if chromosome in validation_chromosomes else "train"
        split_dir = dataset_root / split
        for data_path in sorted(tensor_dir.glob("*_data.npy")):
            labels_path = Path(str(data_path).replace("_data.npy", "_labels.npy"))
            if not is_complete(labels_path):
                raise SystemExit(f"ERROR: missing labels for shard: {data_path}")
            renamed_data = split_dir / f"{chromosome}_{data_path.name}"
            renamed_labels = Path(str(renamed_data).replace("_data.npy", "_labels.npy"))
            if is_complete(renamed_data) and is_complete(renamed_labels):
                total[split] += int(np.load(renamed_labels, mmap_mode="r").shape[0])
                continue

            kept = link_or_filter_shard(
                data_path,
                labels_path,
                renamed_data,
                renamed_labels,
            )
            total[split] += kept

    if total["train"] == 0 or total["val"] == 0:
        raise SystemExit(
            "ERROR: training requires non-empty train and validation splits. "
            "Provide at least one validation chromosome and one training chromosome."
        )
    log(f"training dataset ready: train={total['train']:,}, val={total['val']:,}")


def train(args: argparse.Namespace) -> None:
    truth_vcf_source = require_file(args.truth_vcf, "truth VCF")
    truth_id = graph_fingerprint([truth_vcf_source])
    paths = common_paths(args)
    truth_vcf = stage_truth_vcf(truth_vcf_source, paths["truth"])
    validate_truth_contigs(truth_vcf, args.chromosomes)
    if args.epochs < 5:
        raise SystemExit("ERROR: --epochs must be at least 5 because checkpoints start at epoch 5")

    paths, candidate_nodes, tensor_dirs = prepare_sample(args)
    validation = set(args.val_chromosomes)
    unknown = validation.difference(args.chromosomes)
    if unknown:
        raise SystemExit(f"ERROR: validation chromosomes are not in --chromosomes: {sorted(unknown)}")

    for chromosome, tensor_dir in tensor_dirs.items():
        if labels_match_truth(tensor_dir, truth_id):
            log(f"reuse tensor labels for {chromosome}: {tensor_dir}")
        else:
            clear_stale_labels(tensor_dir)
            run(
                [
                    PYTHON,
                    "-u",
                    "scripts/label_tensors.py",
                    str(tensor_dir / "variant_summary.ndjson"),
                    str(candidate_nodes),
                    str(truth_vcf),
                    "--chr",
                    chromosome,
                    "--data-dir",
                    str(tensor_dir),
                ]
            )
            atomic_json_dump(
                {
                    "truth_fingerprint": truth_id,
                    "truth_vcf": str(truth_vcf_source),
                    "chromosome": chromosome,
                },
                tensor_dir / "labels.manifest.json",
            )
        if not list(tensor_dir.glob("*_labels.npy")):
            raise SystemExit(f"ERROR: labeling produced no label shards for {chromosome}")

    dataset_root = paths["work"] / "training_dataset" / truth_id
    assemble_training_dataset(tensor_dirs, dataset_root, validation)

    model_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else paths["models"]
    model_dir.mkdir(parents=True, exist_ok=True)
    stable_model = model_dir / "best_model.pth"
    if is_complete(stable_model):
        log(f"reuse trained model: {stable_model}")
        return

    trainer = "scripts/train_5channels_npy_pansoma.py"
    train_args = [
        trainer,
        "--data_paths",
        str(dataset_root),
        "--output_path",
        str(model_dir),
        "--epochs",
        str(args.epochs),
        "--lr",
        str(args.learning_rate),
        "--batch_size",
        str(args.batch_size),
        "--num_workers",
        str(args.model_workers),
        "--loss_type",
        args.loss_type,
        "--pos_weight",
        str(args.pos_weight),
        "--gamma",
        str(args.gamma),
        "--training_data_ratio",
        str(args.training_data_ratio),
    ]
    if args.nproc_per_node > 1:
        run(
            [
                "torchrun",
                "--standalone",
                f"--nproc_per_node={args.nproc_per_node}",
                *train_args,
                "--ddp",
            ],
            cwd=ML_ROOT,
        )
    else:
        run([PYTHON, *train_args], cwd=ML_ROOT)

    checkpoints = sorted(model_dir.glob("model_e*_f1_*.pth"), key=lambda path: path.stat().st_mtime)
    if not checkpoints:
        checkpoints = sorted(model_dir.glob("*.pth"), key=lambda path: path.stat().st_mtime)
    if not checkpoints:
        raise SystemExit(f"ERROR: training completed without a checkpoint in {model_dir}")
    if stable_model.exists() or stable_model.is_symlink():
        stable_model.unlink()
    stable_model.symlink_to(checkpoints[-1].name)
    log(f"trained model: {stable_model} -> {checkpoints[-1].name}")


def infer(args: argparse.Namespace) -> None:
    checkpoint = require_file(args.checkpoint, "model checkpoint")
    inference_config = inference_configuration(args, checkpoint)
    output_vcf = (
        Path(args.output_vcf).expanduser().resolve()
        if args.output_vcf
        else Path(args.work_dir).expanduser().resolve() / "results" / f"{args.sample}.vcf.gz"
    )
    if not str(output_vcf).endswith(".vcf.gz"):
        raise SystemExit("ERROR: --output-vcf must end with .vcf.gz")
    output_manifest = Path(str(output_vcf) + ".manifest.json")
    final_expected = {"scope": "merged", "inference": inference_config}
    if indexed_vcf_matches_manifest(output_vcf, output_manifest, final_expected):
        log(f"reuse inference result: {output_vcf}")
        return

    paths, candidate_nodes, tensor_dirs = prepare_sample(args)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    per_chromosome: list[Path] = []
    for chromosome, tensor_dir in tensor_dirs.items():
        out_prefix = paths["results"] / f"{args.sample}.{chromosome}"
        chrom_vcf = Path(str(out_prefix) + ".linear.vcf.gz")
        chrom_manifest = Path(str(chrom_vcf) + ".manifest.json")
        chrom_expected = {
            "scope": "chromosome",
            "chromosome": chromosome,
            "inference": inference_config,
        }
        if not indexed_vcf_matches_manifest(chrom_vcf, chrom_manifest, chrom_expected):
            clear_indexed_vcf(chrom_vcf, chrom_manifest)
            command = [
                PYTHON,
                "scripts/test_5channels_npy_pansoma.py",
                "--input_dir",
                str(tensor_dir),
                "--ckpt",
                str(checkpoint),
                "--out_prefix",
                str(out_prefix),
                "--input_mode",
                "shard",
                "--num_workers",
                str(args.model_workers),
                "--map_json",
                str(candidate_nodes),
                "--variant_summary",
                str(tensor_dir / "variant_summary.ndjson"),
                "--batch_size",
                str(args.batch_size),
                "--device",
                args.device,
                "--emit",
                "linear",
                "--normalize",
                "--min_true_prob",
                str(args.min_true_prob),
                "--min_true_prob_no_anchor",
                str(args.min_true_prob_no_anchor),
            ]
            run(command, cwd=ML_ROOT)
            if not is_complete(chrom_vcf) or not is_complete(Path(str(chrom_vcf) + ".tbi")):
                raise SystemExit(
                    f"ERROR: inference produced no indexed VCF for {chromosome}"
                )
            atomic_json_dump(chrom_expected, chrom_manifest)
        per_chromosome.append(chrom_vcf)

    clear_indexed_vcf(output_vcf, output_manifest)
    merge_linear_vcfs(per_chromosome, output_vcf)
    if not is_complete(output_vcf) or not is_complete(Path(str(output_vcf) + ".tbi")):
        raise SystemExit(f"ERROR: inference merge produced no indexed VCF: {output_vcf}")
    atomic_json_dump(final_expected, output_manifest)
    log(f"inference VCF: {output_vcf}")


def build_node_filters(args: argparse.Namespace) -> None:
    gbz = require_file(args.gbz, "GBZ graph")
    gfa = require_file(args.gfa, "GFA graph")
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(
        {
            "GBZ": str(gbz),
            "GFA": str(gfa),
            "OUTDIR": str(output_dir),
            "CHROMOSOMES": " ".join(args.chromosomes),
        }
    )
    run(["bash", "scripts/build_chr_node_filters.sh"], env=env)
    for chromosome in args.chromosomes:
        chromosome_filters(output_dir, chromosome)
    log(f"chromosome node filters: {output_dir}")


def merge_linear_vcfs(input_vcfs: Sequence[Path], output_vcf: Path) -> None:
    metadata: list[str] = []
    metadata_seen: set[str] = set()
    contigs: list[str] = []
    contig_seen: set[str] = set()
    column_header = None

    for path in input_vcfs:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("##contig="):
                    if line not in contig_seen:
                        contig_seen.add(line)
                        contigs.append(line)
                elif line.startswith("##"):
                    if line not in metadata_seen:
                        metadata_seen.add(line)
                        metadata.append(line)
                elif line.startswith("#CHROM") and column_header is None:
                    column_header = line
                elif not line.startswith("#"):
                    break

    if not metadata or column_header is None:
        raise SystemExit("ERROR: per-chromosome VCF headers are incomplete")

    plain_vcf = Path(str(output_vcf)[:-3])
    with plain_vcf.open("w", encoding="utf-8") as output:
        output.writelines(metadata)
        output.writelines(contigs)
        output.write(column_header)
        for path in input_vcfs:
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                for line in handle:
                    if not line.startswith("#"):
                        output.write(line)

    produced = Path(
        pysam.tabix_index(
            str(plain_vcf),
            preset="vcf",
            force=True,
            keep_original=False,
        )
    )
    if produced.resolve() != output_vcf.resolve():
        produced.replace(output_vcf)
        produced_index = Path(str(produced) + ".tbi")
        if produced_index.exists():
            produced_index.replace(Path(str(output_vcf) + ".tbi"))


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--sample", required=True, help="sample identifier used in output names")
    parser.add_argument("--fastq1", required=True, help="read 1 FASTQ/FASTQ.GZ")
    parser.add_argument("--fastq2", help="read 2 FASTQ/FASTQ.GZ for paired-end data")
    parser.add_argument("--gbz", required=True, help="vg GBZ graph")
    parser.add_argument("--min-index", required=True, help="vg minimizer index (.min)")
    parser.add_argument("--dist-index", required=True, help="vg distance index (.dist)")
    parser.add_argument(
        "--zipcode-index",
        help="vg zipcode index (.zipcodes); required for PacBio HiFi and ONT R10",
    )
    parser.add_argument("--gfa", required=True, help="matching graph GFA with node sequences")
    parser.add_argument(
        "--resource-cache",
        help="optional reusable graph-resource cache; defaults to WORK_DIR/graph_resources",
    )
    parser.add_argument("--work-dir", required=True, help="persistent directory for all intermediates")
    parser.add_argument("--chromosomes", nargs="+", default=[f"chr{i}" for i in range(1, 23)])
    parser.add_argument(
        "--mapper-preset",
        choices=["default", "chaining-sr", "fast", "srold", "hifi", "r10"],
        default="default",
        help="vg giraffe -b preset (default: default)",
    )
    parser.add_argument("--variant-type", choices=["snp", "indel", "all"], default="snp")
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--tensor-workers", type=int, default=8)
    parser.add_argument("--scan-milestone", type=int, default=10_000_000)
    parser.add_argument("--store-milestone", type=int, default=1_000_000)
    parser.add_argument("--min-af", type=float, default=0.08)
    parser.add_argument("--min-variants", type=int, default=3)
    parser.add_argument("--min-allele-bq", type=float, default=10.0)
    parser.add_argument("--min-mapq", type=int, default=10)
    parser.add_argument("--max-indel-len", type=int, default=50)
    parser.add_argument("--shard-size", type=int, default=32768)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pansoma",
        description="Run end-to-end Pansoma workflows from FASTQ files.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    doctor_parser = subparsers.add_parser("doctor", help="validate the container runtime")
    doctor_parser.set_defaults(
        handler=lambda _args: run(
            [PYTHON, "scripts/check_environment.py", "--strict-external", "--ml"]
        )
    )

    filters_parser = subparsers.add_parser(
        "build-node-filters",
        help="build reusable chromosome component/GRCh38 node filters",
    )
    filters_parser.add_argument("--gbz", required=True)
    filters_parser.add_argument("--gfa", required=True)
    filters_parser.add_argument("--output-dir", required=True)
    filters_parser.add_argument("--chromosomes", nargs="+", default=[f"chr{i}" for i in range(1, 23)])
    filters_parser.set_defaults(handler=build_node_filters)

    train_parser = subparsers.add_parser("train", help="FASTQ + truth VCF -> trained PansomaNet model")
    add_common_arguments(train_parser)
    train_parser.add_argument(
        "--truth-vcf",
        required=True,
        help="truth variants in .vcf or .vcf.gz format; compression/indexing is automatic",
    )
    train_parser.add_argument("--val-chromosomes", nargs="+", default=["chr1"])
    train_parser.add_argument("--output-dir", help="model output directory; defaults to WORK_DIR/models")
    train_parser.add_argument("--epochs", type=int, default=70)
    train_parser.add_argument("--learning-rate", type=float, default=1e-4)
    train_parser.add_argument("--batch-size", type=int, default=32)
    train_parser.add_argument("--model-workers", type=int, default=4)
    train_parser.add_argument("--nproc-per-node", type=int, default=1, help="number of GPUs for DDP training")
    train_parser.add_argument(
        "--loss-type",
        choices=["focal", "weighted_ce"],
        default="weighted_ce",
    )
    train_parser.add_argument("--pos-weight", type=float, default=88.0)
    train_parser.add_argument("--gamma", type=float, default=2.0)
    train_parser.add_argument("--training-data-ratio", type=float, default=1.0)
    train_parser.set_defaults(handler=train)

    infer_parser = subparsers.add_parser("infer", help="FASTQ + checkpoint -> indexed linear VCF")
    add_common_arguments(infer_parser)
    infer_parser.add_argument("--checkpoint", required=True)
    infer_parser.add_argument("--output-vcf", help="final .vcf.gz path; defaults to WORK_DIR/results/SAMPLE.vcf.gz")
    infer_parser.add_argument("--batch-size", type=int, default=64)
    infer_parser.add_argument("--model-workers", type=int, default=4)
    infer_parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    infer_parser.add_argument("--min-true-prob", type=float, default=0.5)
    infer_parser.add_argument("--min-true-prob-no-anchor", type=float, default=0.2)
    infer_parser.set_defaults(handler=infer)

    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        args.handler(args)
    except subprocess.CalledProcessError as exc:
        raise SystemExit(f"ERROR: command failed with exit code {exc.returncode}") from exc
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
