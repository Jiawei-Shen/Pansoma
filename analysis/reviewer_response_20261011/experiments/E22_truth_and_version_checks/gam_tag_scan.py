#!/usr/bin/env python3
"""E22 check 2: read-only scan of the first groups of a GAM (BGZF-compressed vg type-tagged protobuf stream).

Prints the tag of every non-"GAM" group (giraffe >= some version writes one PARAMS_JSON group with its parameters
early in the stream) and its payload, and stops after --groups groups. Nothing is written; the file is opened
read-only.
usage: gam_tag_scan.py GAM [--groups N]
"""
import argparse
import gzip


def varint(fh):
    shift = result = 0
    while True:
        b = fh.read(1)
        if not b:
            raise EOFError
        x = b[0]
        result |= (x & 0x7F) << shift
        shift += 7
        if x < 0x80:
            return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gam")
    ap.add_argument("--groups", type=int, default=200)
    a = ap.parse_args()
    seen, first_other = {}, None
    reads = 0
    with open(a.gam, "rb") as raw, gzip.GzipFile(fileobj=raw) as fh:
        g = 0
        try:
            while g < a.groups:
                cnt = varint(fh)
                tag = fh.read(varint(fh)).decode("ascii", "replace")
                payload = None
                for k in range(cnt - 1):
                    n = varint(fh)
                    if tag != "GAM" and k == 0:
                        payload = fh.read(n)
                    else:
                        fh.seek(n, 1)
                seen[tag] = seen.get(tag, 0) + 1
                if tag == "GAM":
                    reads += cnt - 1
                else:
                    print(f"group {g} (after {reads} alignments): tag {tag!r}, {cnt - 1} message(s)")
                    print("  payload:", payload.decode("utf-8", "replace") if payload else None)
                    if first_other is None:
                        first_other = g
                g += 1
        except EOFError:
            print("end of file")
    print(f"{a.gam}: scanned {g} groups, {reads} alignments; tags {seen}")


if __name__ == "__main__":
    main()
