#!/usr/bin/env python3
"""Tag or remove Pansoma VCF records found in panels of normals."""

from __future__ import annotations

import argparse
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable

import pysam


DEFAULT_MATCHING = ("allele", "allele", "position", "position")
PON_NAMES = ("PoN1_gnomAD", "PoN2_dbSNP", "PoN3_1000G", "PoN4_CoLoRSdb")


def chromosome_aliases(chromosome: str) -> tuple[str, ...]:
    """Return common aliases while preserving the requested name first."""
    aliases = [chromosome]
    if chromosome.startswith("chr"):
        aliases.append(chromosome[3:])
    else:
        aliases.append(f"chr{chromosome}")

    if chromosome in {"M", "MT", "chrM", "chrMT"}:
        aliases.extend(("M", "MT", "chrM", "chrMT"))
    return tuple(dict.fromkeys(aliases))


def resolve_contig(vcf: pysam.VariantFile, chromosome: str) -> str | None:
    contigs = vcf.header.contigs
    return next((alias for alias in chromosome_aliases(chromosome) if alias in contigs), None)


def record_matches(
    pon: pysam.VariantFile,
    contig: str,
    position: int,
    ref: str,
    alts: Iterable[str],
    matching: str,
) -> bool:
    """Query one indexed PoN for an exact position or exact allele match."""
    input_alts = {alt.upper() for alt in alts if alt}
    try:
        records = pon.fetch(contig, max(0, position - 1), position)
    except (OSError, ValueError) as exc:
        raise RuntimeError(
            f"Could not query {pon.filename!r}. Confirm that the VCF is BGZF "
            "compressed and has a .tbi or .csi index."
        ) from exc

    for candidate in records:
        if candidate.pos != position:
            continue
        if matching == "position":
            return True
        candidate_alts = {alt.upper() for alt in (candidate.alts or ()) if alt}
        if (candidate.ref or "").upper() == ref.upper() and input_alts & candidate_alts:
            return True
    return False


def add_output_header(header: pysam.VariantHeader) -> None:
    if "PanelOfNormals" not in header.filters:
        header.filters.add(
            "PanelOfNormals",
            None,
            None,
            "Variant matched at least one configured panel of normals",
        )
    if "PANSOMA_PON" not in header.info:
        header.info.add(
            "PANSOMA_PON",
            ".",
            "String",
            "Panel(s) of normals matched by Pansoma",
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Tag Pansoma VCF records found in four indexed PoNs. By default all "
            "records are retained; use --drop-matched to remove matched records."
        )
    )
    parser.add_argument("input_vcf", help="Input .vcf.gz from Pansoma inference")
    parser.add_argument("output_vcf", help="Output BGZF-compressed .vcf.gz")
    parser.add_argument(
        "--pon",
        nargs=4,
        required=True,
        metavar=("PON1", "PON2", "PON3", "PON4"),
        help="Indexed gnomAD, dbSNP, 1000G, and CoLoRSdb PoN VCFs, in that order",
    )
    parser.add_argument(
        "--drop-matched",
        action="store_true",
        help="Remove PoN-matched records instead of retaining them with a FILTER tag",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_path = Path(args.input_vcf).expanduser().resolve()
    output_path = Path(args.output_vcf).expanduser().resolve()
    pon_paths = [Path(path).expanduser().resolve() for path in args.pon]

    if input_path == output_path:
        sys.exit("ERROR: input_vcf and output_vcf must be different paths")
    for path in (input_path, *pon_paths):
        if not path.is_file():
            sys.exit(f"ERROR: VCF not found: {path}")
    if not str(output_path).endswith(".vcf.gz"):
        sys.exit("ERROR: output_vcf must end with .vcf.gz")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    input_vcf = pysam.VariantFile(str(input_path))
    pons = [pysam.VariantFile(str(path)) for path in pon_paths]
    header = input_vcf.header.copy()
    add_output_header(header)

    total = 0
    matched = 0
    written = 0
    matched_by_pon: Counter[str] = Counter()
    contig_cache: dict[tuple[int, str], str | None] = {}
    missing_contigs: set[tuple[str, str]] = set()

    try:
        with pysam.VariantFile(str(output_path), "wz", header=header) as output_vcf:
            for record in input_vcf:
                total += 1
                matches: list[str] = []

                for index, (pon, matching) in enumerate(zip(pons, DEFAULT_MATCHING)):
                    cache_key = (index, record.contig)
                    if cache_key not in contig_cache:
                        contig_cache[cache_key] = resolve_contig(pon, record.contig)
                    pon_contig = contig_cache[cache_key]
                    if pon_contig is None:
                        missing_contigs.add((PON_NAMES[index], record.contig))
                        continue

                    if record_matches(
                        pon,
                        pon_contig,
                        record.pos,
                        record.ref or "",
                        record.alts or (),
                        matching,
                    ):
                        matches.append(PON_NAMES[index])
                        matched_by_pon[PON_NAMES[index]] += 1

                if matches:
                    matched += 1
                    if args.drop_matched:
                        continue
                    record.translate(header)
                    if set(record.filter.keys()) == {"PASS"}:
                        record.filter.clear()
                    record.filter.add("PanelOfNormals")
                    record.info["PANSOMA_PON"] = tuple(matches)
                else:
                    record.translate(header)

                output_vcf.write(record)
                written += 1
    finally:
        input_vcf.close()
        for pon in pons:
            pon.close()

    try:
        pysam.tabix_index(str(output_path), preset="vcf", force=True)
    except (OSError, RuntimeError) as exc:
        sys.exit(f"ERROR: wrote {output_path}, but Tabix indexing failed: {exc}")

    action = "dropped" if args.drop_matched else "tagged"
    print(f"Input records: {total:,}")
    print(f"PoN-matched records {action}: {matched:,}")
    print(f"Output records: {written:,}")
    for name in PON_NAMES:
        print(f"  {name}: {matched_by_pon[name]:,}")
    for name, chromosome in sorted(missing_contigs):
        print(
            f"WARNING: {name} has no contig matching {chromosome!r}; "
            "records on this contig could not be checked",
            file=sys.stderr,
        )
    print(f"Output VCF: {output_path}")
    print(f"Tabix index: {output_path}.tbi")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
