#!/usr/bin/env python3
"""Tag or remove Pansoma VCF records found in panels of normals.

A record matches a PoN when a PoN record has its position, REF and one of its ALTs, under the rule of
the record's kind:
- SNV (every ALT of REF's length): gnomAD and CoLoRSdb with that ALT's AF >= 0.0001; dbSNP when the
  record is not somatic (SAO != 2); 1000G any such record.
- INDEL (an ALT of another length): only gnomAD and CoLoRSdb, with that ALT's AF >= 0.01. dbSNP and
  1000G are not used: most somatic INDELs are homopolymer / STR length changes that exist as population
  alleles at AF 0.001-0.05, where the individual's germline INDELs sit mostly at AF >= 0.05
  (analysis/tensor_recall_20260930/pon/indel_rules)."""

from __future__ import annotations

import argparse
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable

import pysam


PON_NAMES = ("PoN1_gnomAD", "PoN2_dbSNP", "PoN3_1000G", "PoN4_CoLoRSdb")
SNV_MIN_AF = 0.0001
INDEL_MIN_AF = 0.01
PON_RULES = {  # per kind, per PoN in PON_NAMES order; None = the PoN is not used for the kind
    "SNV": (dict(min_af=SNV_MIN_AF), dict(non_somatic=True), dict(), dict(min_af=SNV_MIN_AF)),
    "INDEL": (dict(min_af=INDEL_MIN_AF), None, None, dict(min_af=INDEL_MIN_AF)),
}


def record_kind(ref: str, alts: Iterable[str]) -> str:
    """INDEL when any (non-symbolic) ALT differs from REF in length, else SNV."""
    for alt in alts:
        if not alt or alt.startswith("<") or alt in ("*", "."):
            continue
        if len(alt) != len(ref):
            return "INDEL"
    return "SNV"


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
    min_af: float | None = None,
    non_somatic: bool = False,
) -> bool:
    """Query one indexed PoN for an exact allele match (with the ALT's AF >= min_af; not SAO=2 if non_somatic)."""
    input_alts = {alt.upper() for alt in alts if alt}
    try:
        records = pon.fetch(contig, max(0, position - 1), position)
    except (OSError, ValueError) as exc:
        raise RuntimeError(
            f"Could not query {pon.filename!r}. Confirm that the VCF is BGZF "
            "compressed and has a .tbi or .csi index."
        ) from exc

    for candidate in records:
        if candidate.pos != position or (candidate.ref or "").upper() != ref.upper():
            continue
        if non_somatic and candidate.info.get("SAO") == 2:
            continue
        candidate_alts = [(alt or "").upper() for alt in (candidate.alts or ())]
        hits = [k for k, alt in enumerate(candidate_alts) if alt in input_alts]
        if not hits:
            continue
        if min_af is None:
            return True
        af = candidate.info.get("AF")
        af = af if isinstance(af, tuple) else (af,) * len(candidate_alts)
        # VCF Float is 32-bit: an AF written as 0.01 reads back as 0.0099999998, so compare with a tolerance
        if any(af[k] is not None and af[k] >= min_af * (1 - 1e-6) for k in hits):
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
            "Tag Pansoma VCF records found in four indexed PoNs (SNV: gnomAD / CoLoRSdb AF >= "
            f"{SNV_MIN_AF}, dbSNP non-somatic, 1000G; INDEL: gnomAD / CoLoRSdb AF >= {INDEL_MIN_AF} only). "
            "By default all records are retained; use --drop-matched to remove matched records."
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
    by_kind: Counter[str] = Counter()
    matched_by_kind: Counter[str] = Counter()
    matched_by_pon: Counter[str] = Counter()
    contig_cache: dict[tuple[int, str], str | None] = {}
    missing_contigs: set[tuple[str, str]] = set()

    try:
        with pysam.VariantFile(str(output_path), "wz", header=header) as output_vcf:
            for record in input_vcf:
                total += 1
                kind = record_kind(record.ref or "", record.alts or ())
                by_kind[kind] += 1
                matches: list[str] = []

                for index, (pon, rule) in enumerate(zip(pons, PON_RULES[kind])):
                    if rule is None:
                        continue
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
                        **rule,
                    ):
                        matches.append(PON_NAMES[index])
                        matched_by_pon[PON_NAMES[index]] += 1

                if matches:
                    matched += 1
                    matched_by_kind[kind] += 1
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
    print(f"Input records: {total:,} (SNV {by_kind['SNV']:,}, INDEL {by_kind['INDEL']:,})")
    print(f"PoN-matched records {action}: {matched:,} (SNV {matched_by_kind['SNV']:,}, INDEL {matched_by_kind['INDEL']:,})")
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
