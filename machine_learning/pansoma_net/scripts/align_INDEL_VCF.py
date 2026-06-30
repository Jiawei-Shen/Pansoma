#!/usr/bin/env python3
"""
align_INDEL_VCF.py

Strictly anchor Pansoma-style INDEL VCF against a linear FASTA.

Updated:
  - Sorts all output records by contig order in VCF header, then by position,
    BEFORE writing final VCF and creating tabix index.
  - This prevents index creation failure caused by reordered deletion records
    (e.g. POS -> POS-1).

Your requested behavior (exact):
  1) TYPE=I (insertion):
       - ALT must start with REF
       - If ALT is missing the anchor prefix: ALT = REF + ALT
       - If ALT == REF (lost inserted bases): ALT = REF + REF  (e.g., A->AA)
       - POS does NOT change
       - IMPORTANT: Every TYPE=I is considered "fixed" (processed), even if already OK.

  2) TYPE=D (deletion):
       - Add exactly ONE left-flanking base to BOTH REF and ALT
       - NEW_POS = POS - 1
       - NEW_REF = left_base + deleted_seq
       - NEW_ALT = left_base
       - deleted_seq is fetched from FASTA at POS with length = len(original REF)
       - IMPORTANT: Every TYPE=D is considered "fixed" (processed).

Other:
  - Do NOT merge anything. One input record => one output record (for I/D).
  - Preserve INFO; fix single-valued INFO/AD (e.g. AD=59) by moving it to INFO/ADALT.
  - Adds INFO/FIXED + INFO/FIXNOTE for ALL TYPE=I and TYPE=D records (even if unchanged).
  - Writes bgzip VCF + tabix index unless --no-index.

Example:
  python scripts/align_INDEL_VCF.py \
    --vcf  /path/in.vcf.gz \
    --fasta /path/ref.fa \
    --out  /path/out.realigned.vcf.gz
"""

import argparse
import sys
from typing import Dict, List, Any

import pysam


def build_contig_mapper(vcf_contigs: List[str], fasta_contigs: List[str]) -> Dict[str, str]:
    """Map VCF contig names to FASTA contigs using exact/chr conversions."""
    fasta_set = set(fasta_contigs)
    m: Dict[str, str] = {}
    for c in vcf_contigs:
        if c in fasta_set:
            m[c] = c
        elif c.startswith("chr") and c[3:] in fasta_set:
            m[c] = c[3:]
        elif ("chr" + c) in fasta_set:
            m[c] = "chr" + c
        else:
            m[c] = c
    return m


def fetch_ref(fa: pysam.FastaFile, contig: str, pos_1based: int, length: int) -> str:
    """Fetch reference substring from FASTA using 1-based position."""
    if pos_1based < 1:
        raise ValueError(f"Invalid position {pos_1based} (must be >=1)")
    start0 = pos_1based - 1
    return fa.fetch(contig, start0, start0 + length).upper()


def make_header(in_header: pysam.VariantHeader) -> pysam.VariantHeader:
    out = in_header.copy()

    if "FIXED" not in out.info:
        out.add_line('##INFO=<ID=FIXED,Number=0,Type=Flag,Description="Record processed/fixed against FASTA by align_INDEL_VCF.py">')
    if "FIXNOTE" not in out.info:
        out.add_line('##INFO=<ID=FIXNOTE,Number=1,Type=String,Description="Fix/process note">')

    if "ADALT" not in out.info:
        out.add_line('##INFO=<ID=ADALT,Number=1,Type=Integer,Description="Alt allele depth (from original INFO/AD when AD was single-valued)">')

    return out


def _as_list_of_ints(x: Any) -> List[int]:
    if x is None:
        return []
    if isinstance(x, (tuple, list)):
        out = []
        for v in x:
            try:
                out.append(int(v))
            except Exception:
                pass
        return out
    if isinstance(x, (int, float)):
        return [int(x)]
    try:
        return [int(str(x))]
    except Exception:
        return []


def sanitize_info(rec: pysam.VariantRecord) -> Dict[str, Any]:
    """Copy INFO safely; move single-valued INFO/AD to INFO/ADALT."""
    info: Dict[str, Any] = {}
    for k in rec.info.keys():
        try:
            info[k] = rec.info[k]
        except Exception:
            continue

    if "AD" in info:
        ad_list = _as_list_of_ints(info["AD"])
        if len(ad_list) == 1:
            info["ADALT"] = int(ad_list[0])
            info.pop("AD", None)

    return info


def serialize_record_fields(rec: pysam.VariantRecord) -> Dict[str, Any]:
    """Copy a record into plain Python fields so it can be sorted and re-written later."""
    return {
        "contig": rec.contig,
        "start": rec.start,
        "stop": rec.stop,
        "id": rec.id,
        "qual": rec.qual,
        "alleles": tuple(rec.alleles),
        "filter": list(rec.filter.keys()),
        "info": dict(rec.info),
        "samples": {
            s: {k: rec.samples[s][k] for k in rec.samples[s].keys()}
            for s in rec.samples
        },
    }


def main():
    ap = argparse.ArgumentParser(description="Anchor Pansoma INDEL VCF against FASTA (strict rules).")
    ap.add_argument("--vcf", required=True, help="Input VCF(.gz)")
    ap.add_argument("--fasta", required=True, help="Reference FASTA (must have .fai)")
    ap.add_argument("--out", required=True, help="Output VCF.GZ")
    ap.add_argument("--log", default=None, help="TSV log (default: <out>.fix_log.tsv)")
    ap.add_argument("--no-index", action="store_true", help="Do not create .tbi index")
    ap.add_argument("--force-ref-from-fasta", action="store_true",
                    help="For TYPE=I, overwrite REF with FASTA at POS (len=original REF).")

    args = ap.parse_args()
    log_path = args.log if args.log else (args.out + ".fix_log.tsv")

    invcf = pysam.VariantFile(args.vcf)
    fa = pysam.FastaFile(args.fasta)

    out_header = make_header(invcf.header)
    contig_map = build_contig_mapper(list(invcf.header.contigs), list(fa.references))

    contig_order = {name: i for i, name in enumerate(out_header.contigs)}
    default_contig_rank = len(contig_order)

    out_records: List[Dict[str, Any]] = []

    fixed = 0
    kept = 0
    dropped = 0

    with open(log_path, "w", encoding="utf-8") as logf:
        logf.write("#chrom\told_pos\told_ref\told_alt\tnew_pos\tnew_ref\tnew_alt\taction\tnote\n")

        for rec in invcf.fetch():
            chrom = rec.contig
            fasta_chrom = contig_map.get(chrom, chrom)
            if fasta_chrom not in fa.references:
                raise SystemExit(f"ERROR: contig '{chrom}' not in FASTA (tried '{fasta_chrom}')")

            if rec.alts is None or len(rec.alts) != 1:
                dropped += 1
                logf.write(f"{chrom}\t{rec.pos}\t{rec.ref}\t{rec.alts}\t.\t.\t.\tDROP\tmissing_or_multiallelic\n")
                continue

            alt_in = rec.alts[0]
            ref_in = rec.ref
            vtype = str(rec.info.get("TYPE", "")).upper()
            info = sanitize_info(rec)

            # ---------------- TYPE=I: insertion ----------------
            if vtype == "I":
                ref = ref_in.upper()
                alt = alt_in.upper()

                if args.force_ref_from_fasta:
                    try:
                        ref = fetch_ref(fa, fasta_chrom, rec.pos, len(ref))
                    except Exception as e:
                        logf.write(f"{chrom}\t{rec.pos}\t{ref_in}\t{alt_in}\t.\t.\t.\tWARN\tforce_ref_failed:{e}\n")

                note = "ins_ok"
                if not alt.startswith(ref):
                    alt = ref + alt
                    note = "ins_prefix_ref_to_alt"
                else:
                    if len(alt) == len(ref):
                        alt = ref + alt
                        note = "ins_alt_eq_ref_make_refref"

                new = {
                    "contig": rec.contig,
                    "start": rec.start,
                    "stop": rec.start + len(ref),
                    "id": rec.id,
                    "qual": rec.qual,
                    "alleles": (ref, alt),
                    "filter": list(rec.filter.keys()),
                    "info": dict(info),
                    "samples": {
                        s: {k: rec.samples[s][k] for k in rec.samples[s].keys()}
                        for s in rec.samples
                    },
                }
                new["info"]["FIXED"] = True
                new["info"]["FIXNOTE"] = note

                out_records.append(new)

                fixed += 1
                kept += 1
                logf.write(f"{chrom}\t{rec.pos}\t{ref_in}\t{alt_in}\t{rec.pos}\t{ref}\t{alt}\tFIX_INS\t{note}\n")
                continue

            # ---------------- TYPE=D: deletion ----------------
            if vtype == "D":
                pos = rec.pos
                if pos <= 1:
                    dropped += 1
                    logf.write(f"{chrom}\t{pos}\t{ref_in}\t{alt_in}\t.\t.\t.\tDROP\tDEL_no_left_anchor\n")
                    continue

                try:
                    left_base = fetch_ref(fa, fasta_chrom, pos - 1, 1)
                    del_len = len(ref_in)
                    deleted_seq = fetch_ref(fa, fasta_chrom, pos, del_len)
                except Exception as e:
                    dropped += 1
                    logf.write(f"{chrom}\t{pos}\t{ref_in}\t{alt_in}\t.\t.\t.\tDROP\tDEL_fasta_fetch_failed:{e}\n")
                    continue

                new_pos = pos - 1
                new_ref = (left_base + deleted_seq).upper()
                new_alt = left_base.upper()

                start0 = new_pos - 1
                stop0 = start0 + len(new_ref)

                new = {
                    "contig": rec.contig,
                    "start": start0,
                    "stop": stop0,
                    "id": rec.id,
                    "qual": rec.qual,
                    "alleles": (new_ref, new_alt),
                    "filter": list(rec.filter.keys()),
                    "info": dict(info),
                    "samples": {
                        s: {k: rec.samples[s][k] for k in rec.samples[s].keys()}
                        for s in rec.samples
                    },
                }
                new["info"]["FIXED"] = True
                new["info"]["FIXNOTE"] = "del_add_left_base"

                out_records.append(new)

                fixed += 1
                kept += 1
                logf.write(f"{chrom}\t{pos}\t{ref_in}\t{alt_in}\t{new_pos}\t{new_ref}\t{new_alt}\tFIX_DEL\tdel_add_left_base\n")
                continue

            # ---------------- Other types: pass through ----------------
            new = {
                "contig": rec.contig,
                "start": rec.start,
                "stop": rec.stop,
                "id": rec.id,
                "qual": rec.qual,
                "alleles": tuple(rec.alleles),
                "filter": list(rec.filter.keys()),
                "info": dict(info),
                "samples": {
                    s: {k: rec.samples[s][k] for k in rec.samples[s].keys()}
                    for s in rec.samples
                },
            }

            out_records.append(new)
            kept += 1
            logf.write(f"{chrom}\t{rec.pos}\t{ref_in}\t{alt_in}\t{rec.pos}\t{rec.ref}\t{alt_in}\tKEEP\tunhandled_type:{vtype}\n")

    invcf.close()
    fa.close()

    # -------- sort before write --------
    out_records.sort(
        key=lambda r: (
            contig_order.get(r["contig"], default_contig_rank),
            r["start"],
            r["stop"],
            r["alleles"][0],
            r["alleles"][1] if len(r["alleles"]) > 1 else ""
        )
    )

    outvcf = pysam.VariantFile(args.out, "wz", header=out_header)

    for r in out_records:
        new_rec = outvcf.new_record(
            contig=r["contig"],
            start=r["start"],
            stop=r["stop"],
            id=r["id"],
            qual=r["qual"],
            alleles=r["alleles"],
            filter=r["filter"],
            info=r["info"],
        )
        for s, sample_data in r["samples"].items():
            new_rec.samples[s].update(sample_data)
        outvcf.write(new_rec)

    outvcf.close()

    if not args.no_index:
        pysam.tabix_index(args.out, preset="vcf", force=True)

    print(f"[done] out={args.out}", file=sys.stderr)
    print(f"[done] log={log_path}", file=sys.stderr)
    print(f"[stats] kept={kept} fixed={fixed} dropped={dropped}", file=sys.stderr)
    if not args.no_index:
        print(f"[done] index={args.out}.tbi", file=sys.stderr)


if __name__ == "__main__":
    main()