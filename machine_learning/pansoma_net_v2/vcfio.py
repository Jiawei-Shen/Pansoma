"""Plain-text VCF helpers for graph_vcf, linear_vcf and vcfeval: a record is a list of its tab-separated columns,
INFO a dict in column order, and every output is bgzipped and tabix-indexed with pysam."""
import gzip
import os
import sys
from pathlib import Path

import pysam

REPO = Path(__file__).resolve().parents[2]


def import_pipeline():
    """Make indexed_gam_pipeline_v4 importable (the CLIs run from machine_learning/); its ReferencePath and
    GraphIndex are the ones that built and merged the tensors."""
    if str(REPO) not in sys.path:
        sys.path.append(str(REPO))


def parse_info(text):
    if text == ".":
        return {}
    info = {}
    for item in text.split(";"):
        key, _, value = item.partition("=")
        info[key] = value if _ else True
    return info


def format_info(info):
    items = [key if value is True else f"{key}={value}" for key, value in info.items()
             if value is not None and value is not False]
    return ";".join(items) or "."


def read_vcf(path):
    """(header lines without the #CHROM line, the #CHROM line, iterator over records as column lists)."""
    stream = gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)
    header = []
    for line in stream:
        if line.startswith("#CHROM"):
            columns = line.rstrip("\n")
            break
        header.append(line.rstrip("\n"))
    else:
        stream.close()
        raise ValueError(f"{path}: no #CHROM line")

    def records():
        with stream:
            for line in stream:
                if line.strip():
                    yield line.rstrip("\n").split("\t")
    return header, columns, records()


def meta(header, key):
    """Value of a ##<key>=value header line, or None."""
    prefix = f"##{key}="
    return next((line[len(prefix):] for line in header if line.startswith(prefix)), None)


def write_vcf(path, header, columns, records):
    """Write records (column lists, already sorted) to path (.vcf.gz), bgzipped, with a .tbi."""
    path = str(path)
    if not path.endswith(".vcf.gz"):
        raise ValueError(f"{path}: the output must end in .vcf.gz")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    plain = path[:-3]
    with open(plain, "w") as f:
        for line in header:
            f.write(line + "\n")
        f.write(columns + "\n")
        for record in records:
            f.write("\t".join(record) + "\n")
    for old in (path, path + ".tbi"):
        if os.path.exists(old):
            os.remove(old)
    pysam.tabix_index(plain, preset="vcf", force=True)  # compresses plain to path and removes plain
    return path


def stats_path(vcf):
    """<name>.stats.json next to <name>.vcf.gz."""
    vcf = str(vcf)
    return Path((vcf[:-len(".vcf.gz")] if vcf.endswith(".vcf.gz") else vcf) + ".stats.json")
