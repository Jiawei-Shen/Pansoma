#!/usr/bin/env python3
import argparse
import sys
from typing import List, Optional
import math

import matplotlib.pyplot as plt
import pysam


def is_symbolic(allele: Optional[str]) -> bool:
    if allele is None:
        return True
    return allele.startswith("<") and allele.endswith(">")

def classify_variant(record) -> str:
    """Return 'snv' if all non-symbolic ALT alleles are single-base changes,
    'indel' if any length differs, otherwise 'other' (e.g., symbolic)."""
    ref = record.ref
    alts = record.alts or ()
    if any(is_symbolic(a) for a in alts):
        return "other"
    if all(len(ref) == 1 and len(a) == 1 for a in alts):
        return "snv"
    if any(len(a) != len(ref) for a in alts):
        return "indel"
    return "other"

def extract_vafs_from_sample(record, sample_name, per_alt_mode="max"):
    """
    Try FORMAT/AF first; if absent, fall back to FORMAT/AD.
    per_alt_mode: 'max' | 'first' | 'all'
    Returns a list of VAFs for this record based on the mode.
    """
    vafs: List[float] = []

    if sample_name not in record.samples:
        return vafs

    sm = record.samples[sample_name]

    # 1) FORMAT/AF
    af = sm.get("AF")
    if af is not None:
        # af can be a float or list/tuple
        vals = af if isinstance(af, (list, tuple)) else [af]
        vals = [float(x) for x in vals if x is not None and not math.isnan(float(x))]
        vafs = vals

    # 2) FORMAT/AD (ref,alt1,alt2,...)
    elif sm.get("AD") is not None:
        ad = sm.get("AD")
        # pysam returns a tuple of ints for AD (ref, alt1, alt2, ...)
        if ad and len(ad) >= 2:
            ref_dp = ad[0] if ad[0] is not None else 0
            for i in range(1, len(ad)):
                alt_dp = ad[i] if ad[i] is not None else 0
                dp = ref_dp + alt_dp
                if dp > 0:
                    vafs.append(alt_dp / dp)

    # Collapse per_alt_mode
    if not vafs:
        return []

    if per_alt_mode == "first":
        return [vafs[0]]
    elif per_alt_mode == "max":
        return [max(vafs)]
    elif per_alt_mode == "all":
        return vafs
    else:
        return [max(vafs)]

def extract_vafs_site_level(record, per_alt_mode="max"):
    """
    Use INFO/AF if present (site-level). Return per_alt_mode-applied list.
    """
    v = record.info.get("AF")
    if v is None:
        return []
    vals = v if isinstance(v, (list, tuple)) else [v]
    vals = [float(x) for x in vals if x is not None and not math.isnan(float(x))]
    if not vals:
        return []
    if per_alt_mode == "first":
        return [vals[0]]
    elif per_alt_mode == "max":
        return [max(vals)]
    elif per_alt_mode == "all":
        return vals
    else:
        return [max(vals)]

def passes_min_dp(record, sample_name: Optional[str], min_dp: int) -> bool:
    if min_dp <= 0:
        return True
    # Prefer per-sample DP if sample provided; else site-level DP if present
    if sample_name and sample_name in record.samples:
        sm = record.samples[sample_name]
        dp = sm.get("DP")
        if dp is not None:
            try:
                return int(dp) >= min_dp
            except Exception:
                pass
    # Site-level DP as fallback
    site_dp = record.info.get("DP")
    if site_dp is not None:
        try:
            return int(site_dp) >= min_dp
        except Exception:
            pass
    # If no DP fields, don't filter it out based on DP
    return True

def main():
    p = argparse.ArgumentParser(
        description="Plot VAF distribution from a VCF (supports INFO/AF, FORMAT/AF, or FORMAT/AD)."
    )
    p.add_argument("vcf", help="Input VCF/BCF (.vcf, .vcf.gz, or .bcf)")
    p.add_argument("-s", "--sample", help="Sample name (for per-sample AF/AD). If omitted, tries INFO/AF or first sample.")
    p.add_argument("-t", "--type", choices=["snv", "indel", "both"], default="both",
                   help="Variant type to include (default: both)")
    p.add_argument("--pass-only", action="store_true", help="Keep only PASS variants")
    p.add_argument("--min-dp", type=int, default=0, help="Minimum DP (sample DP preferred, else site DP). Default: 0")
    p.add_argument("--per-alt", choices=["max", "first", "all"], default="max",
                   help="How to aggregate multi-allelic VAFs (default: max)")
    p.add_argument("--bins", type=int, default=60, help="Histogram bins (default: 60)")
    p.add_argument("-o", "--out", default="vaf_hist.png", help="Output PNG file (default: vaf_hist.png)")
    p.add_argument("--title", default=None, help="Custom plot title")
    p.add_argument("--xmax", type=float, default=1.0, help="Max X (default: 1.0)")
    p.add_argument("--xmin", type=float, default=0.0, help="Min X (default: 0.0)")
    args = p.parse_args()

    # Open VCF/BCF
    try:
        vcf = pysam.VariantFile(args.vcf)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to open {args.vcf}: {e}\n")
        sys.exit(1)

    # Resolve sample
    chosen_sample = args.sample
    if chosen_sample is None and len(vcf.header.samples) > 0:
        # If user didn't specify, prefer per-sample fields by taking the first sample
        chosen_sample = vcf.header.samples[0]

    vafs: List[float] = []
    total = 0
    kept = 0

    for rec in vcf.fetch():
        total += 1

        if args.pass_only and rec.filter.keys() not in (set(), {"PASS"}):
            # pysam represents filters as a dict-like; PASS often shows as empty set() or {"PASS"}
            continue

        vtype = classify_variant(rec)
        if args.type == "snv" and vtype != "snv":
            continue
        if args.type == "indel" and vtype != "indel":
            continue
        if vtype == "other":  # skip symbolic or unusual types
            continue

        if not passes_min_dp(rec, chosen_sample, args.min_dp):
            continue

        # Prefer per-sample (if we have a sample) over site-level
        vals = []
        if chosen_sample and chosen_sample in rec.samples:
            vals = extract_vafs_from_sample(rec, chosen_sample, per_alt_mode=args.per_alt)

        if not vals:
            vals = extract_vafs_site_level(rec, per_alt_mode=args.per_alt)

        # Keep only valid [0,1] VAFs
        vals = [x for x in vals if x is not None and 0.0 <= x <= 1.0 and not math.isnan(x)]
        if not vals:
            continue

        kept += 1
        vafs.extend(vals)

    # Summary
    sys.stderr.write(f"[INFO] Records scanned: {total}\n")
    sys.stderr.write(f"[INFO] Records kept after filters: {kept}\n")
    sys.stderr.write(f"[INFO] VAF values collected: {len(vafs)}\n")

    if len(vafs) == 0:
        sys.stderr.write("[WARN] No VAF values found after filtering. Nothing to plot.\n")
        sys.exit(2)

    # Plot
    plt.figure()
    plt.hist(vafs, bins=args.bins, range=(args.xmin, args.xmax))
    ttl = args.title
    if not ttl:
        base = f"VAF Distribution ({args.type.upper()})"
        src = "FORMAT/AF|AD or INFO/AF"
        if args.pass_only:
            base += " | PASS"
        if args.min_dp > 0:
            base += f" | DP≥{args.min_dp}"
        ttl = f"{base} | {src}"
        if chosen_sample:
            ttl += f" | sample={chosen_sample}"
    plt.title(ttl)
    plt.xlabel("VAF")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(args.out, dpi=200)
    sys.stderr.write(f"[INFO] Saved plot to {args.out}\n")


if __name__ == "__main__":
    main()
