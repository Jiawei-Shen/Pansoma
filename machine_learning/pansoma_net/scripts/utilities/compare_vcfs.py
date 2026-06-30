#!/usr/bin/env python3
"""
compare_vcfs.py
---------------
Compare two VCF-like files and report how many variant records occur in both.

Features:
- Uses pysam.VariantFile for standard VCF/BCF when possible.
- Falls back to a tolerant text parser for "sites" VCFs that only have 5 columns
  (#CHROM POS ID REF ALT) or have whitespace issues.
- Supports allele-exact (CHROM, POS, REF, ALT) or position-only (CHROM, POS) overlap.
- Expands multi-allelic records (ALT with commas) into separate records.
- Optionally writes the intersection keys to a TSV.

Notes:
- No reference normalization/left-trimming is performed.
- Chromosomes are canonicalized by stripping a leading 'chr' (case-insensitive).
"""

import argparse
import sys
import re
from typing import Iterable, Iterator, Set, Tuple, Optional

try:
    import pysam  # type: ignore
except Exception:
    pysam = None  # We'll fall back to text parsing if not available

KeyAllele = Tuple[str, int, str, str]
KeyPos = Tuple[str, int]


def canon_chrom(chrom: str) -> str:
    """Canonicalize chromosome by stripping leading 'chr' (common cross-dataset discrepancy)."""
    if chrom.lower().startswith("chr"):
        return chrom[3:]
    return chrom


def iter_keys_pysam(path: str, pos_only: bool) -> Iterator[Tuple]:
    """
    Yield keys from a VCF/BCF using pysam.
    Keys are either (chrom, pos) or (chrom, pos, ref, alt), with chrom canonicalized and pos 1-based.
    """
    if pysam is None:
        raise RuntimeError("pysam is not available")

    try:
        vf = pysam.VariantFile(path)
    except Exception as e:
        raise RuntimeError(f"pysam failed to open '{path}': {e}")

    for rec in vf.fetch():  # streaming; index not strictly required for full scan
        chrom = canon_chrom(str(rec.chrom))
        pos = int(rec.pos)  # already 1-based
        ref = str(rec.ref) if rec.ref is not None else ""
        alts = rec.alts if rec.alts is not None else ()
        if pos_only:
            yield (chrom, pos)
        else:
            if not alts:
                # No ALT; skip
                continue
            for alt in alts:
                if alt is None or alt == ".":
                    continue
                yield (chrom, pos, ref, str(alt))


def looks_like_sites_header(header_line: str) -> bool:
    """
    Determine if a header line corresponds to a 5-column sites header:
    '#CHROM POS ID REF ALT'
    """
    # Split by any whitespace to handle tabs/spaces confusion
    cols = re.split(r"\s+", header_line.strip())
    return (
        len(cols) >= 5
        and cols[0].lstrip("#").upper() == "CHROM"
        and cols[1].upper() == "POS"
        and cols[2].upper() == "ID"
        and cols[3].upper() == "REF"
        and cols[4].upper() == "ALT"
    )


def iter_keys_text(path: str, pos_only: bool) -> Iterator[Tuple]:
    """
    Tolerant text-based iterator:
    - Accepts gz or plain text (pysam handles gz transparently; here we just open plain).
    - Skips meta lines starting with '##'.
    - Detects header and allows both 5-col sites header and 8+ col standard header.
    - Splits on arbitrary whitespace (tabs or spaces).
    - Expands multi-allelic ALTs.
    """
    # We won't use gzip here to keep deps minimal; pysam handles gz best.
    # But the user files are typically .vcf.gz; if opening fails, suggest using pysam or zcat.
    opener = open
    if path.endswith(".gz"):
        # Try via pysam's TabixFile to stream lines robustly without parsing VCF fields
        if pysam is None:
            raise RuntimeError(
                f"Cannot read compressed file '{path}' without pysam. "
                f"Install pysam or provide an uncompressed file."
            )
        try:
            tbx = pysam.TabixFile(path)
        except Exception as e:
            raise RuntimeError(f"Failed to tabix-open '{path}': {e}")

        header_processed = False
        sites_mode = False
        for line in tbx.header:  # header lines only
            line = line.rstrip("\n")
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header_processed = True
                if looks_like_sites_header(line):
                    sites_mode = True
                else:
                    # Standard VCF header assumed (8+ columns)
                    sites_mode = False
        if not header_processed:
            # No header lines exposed by Tabix (some site files might not be tabix-indexed properly)
            # We'll still try to parse records below.
            pass

        # iterate records
        for rec in tbx.fetch():
            if not rec or rec.startswith("##"):
                continue
            if rec.startswith("#CHROM"):
                # already handled
                continue
            cols = re.split(r"\s+", rec.strip())
            if len(cols) < 5:
                # malformed; skip
                continue
            chrom = canon_chrom(cols[0])
            try:
                pos = int(cols[1])
            except Exception:
                continue
            if pos_only:
                yield (chrom, pos)
            else:
                if sites_mode or len(cols) == 5:
                    # sites-like: CHROM POS ID REF ALT
                    ref = cols[3]
                    for alt in cols[4].split(","):
                        alt = alt.strip()
                        if alt and alt != ".":
                            yield (chrom, pos, ref, alt)
                else:
                    # standard-ish: CHROM POS ID REF ALT QUAL FILTER INFO ...
                    if len(cols) < 8:
                        # not enough columns; treat as sites-like fallback
                        ref = cols[3]
                        for alt in cols[4].split(","):
                            alt = alt.strip()
                            if alt and alt != ".":
                                yield (chrom, pos, ref, alt)
                    else:
                        ref = cols[3]
                        for alt in cols[4].split(","):
                            alt = alt.strip()
                            if alt and alt != ".":
                                yield (chrom, pos, ref, alt)
        return  # end gz path

    # Plain text path
    header_processed = False
    sites_mode = False
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header_processed = True
                sites_mode = looks_like_sites_header(line)
                break
        # Now parse rest of file (or entire file if header missing)
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            if raw.startswith("##"):
                continue
            if raw.startswith("#CHROM"):
                # skip header row when encountered
                continue
            rec = raw.strip()
            if not rec:
                continue
            cols = re.split(r"\s+", rec)
            if len(cols) < 5:
                continue
            chrom = canon_chrom(cols[0])
            try:
                pos = int(cols[1])
            except Exception:
                continue
            if pos_only:
                yield (chrom, pos)
            else:
                if sites_mode or len(cols) == 5:
                    ref = cols[3]
                    for alt in cols[4].split(","):
                        alt = alt.strip()
                        if alt and alt != ".":
                            yield (chrom, pos, ref, alt)
                else:
                    # standard-ish (8+ cols) or partial; be tolerant
                    ref = cols[3]
                    for alt in cols[4].split(","):
                        alt = alt.strip()
                        if alt and alt != ".":
                            yield (chrom, pos, ref, alt)


def collect_keys(path: str, pos_only: bool) -> Set[Tuple]:
    """
    Try pysam first; if it fails due to header/format, fall back to tolerant text mode.
    """
    # First attempt: pysam VariantFile
    if pysam is not None:
        try:
            return set(iter_keys_pysam(path, pos_only))
        except Exception as e:
            # Fall through to text mode
            sys.stderr.write(f"[info] Falling back to tolerant text parser for '{path}': {e}\n")
    else:
        sys.stderr.write("[info] pysam not available; using tolerant text parser.\n")

    # Fallback: tolerant text mode
    return set(iter_keys_text(path, pos_only))


def main():
    ap = argparse.ArgumentParser(description="Compare two VCF-like files and count overlaps.")
    ap.add_argument("vcf1", help="First VCF/BCF/VCF.gz (standard or sites-like)")
    ap.add_argument("vcf2", help="Second VCF/BCF/VCF.gz (standard or sites-like)")
    ap.add_argument("--pos-only", action="store_true",
                    help="Compare using position only (CHROM, POS) instead of allele-exact")
    ap.add_argument("-o", "--output-tsv", default=None,
                    help="Optional path to write intersect keys as TSV")
    args = ap.parse_args()

    pos_only = args.pos_only

    s1 = collect_keys(args.vcf1, pos_only)
    s2 = collect_keys(args.vcf2, pos_only)

    inter = s1 & s2
    only1 = s1 - s2
    only2 = s2 - s1

    print("=== Comparison summary ===")
    print(f"File 1: {args.vcf1}")
    print(f"File 2: {args.vcf2}")
    print(f"Mode:   {'POSITION-ONLY (CHROM,POS)' if pos_only else 'ALLELE-EXACT (CHROM,POS,REF,ALT)'}")
    print()
    print(f"File1 variants: {len(s1)}")
    print(f"File2 variants: {len(s2)}")
    print(f"Intersection:   {len(inter)}")
    print(f"Only in File1:  {len(only1)}")
    print(f"Only in File2:  {len(only2)}")

    if args.output_tsv:
        with open(args.output_tsv, "w", encoding="utf-8") as out:
            if pos_only:
                out.write("CHROM\tPOS\n")
                for chrom, pos in sorted(inter, key=lambda x: (x[0], x[1])):
                    out.write(f"{chrom}\t{pos}\n")
            else:
                out.write("CHROM\tPOS\tREF\tALT\n")
                for chrom, pos, ref, alt in sorted(inter, key=lambda x: (x[0], x[1], x[2], x[3])):
                    out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\n")
        print(f"\nWrote intersection to: {args.output_tsv}")


if __name__ == "__main__":
    main()
