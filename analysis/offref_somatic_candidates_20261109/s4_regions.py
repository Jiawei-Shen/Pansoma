"""s4d: merged fetch regions of the loci that one label has not counted yet (e.g. loci added to loci.tsv later), for a
remote fetch of the new regions only (s4_fetch_hg008n.py LABEL REGIONS_BED TAG).

usage: s4_regions.py TUMOR LABEL OUT_BED [SELECT]      (SELECT all / cand / germ, as in s4_count)
"""
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA  # noqa: E402


def main():
    tumor, label, out = sys.argv[1], sys.argv[2], sys.argv[3]
    select = sys.argv[4] if len(sys.argv) > 4 else 'all'
    d = f'{DATA}/s4/{tumor}'
    done = set()
    for p in glob.glob(f'{d}/counts/{label}.*.tsv'):
        with open(p) as f:
            next(f)
            done.update(line.split('\t', 1)[0] for line in f)
    iv = []
    for r in csv.DictReader(open(f'{d}/loci.tsv'), delimiter='\t'):
        if r['locus_id'] in done or (select != 'all' and r.get(select) != '1'):
            continue
        iv.append((r['chrom'], int(r['start']), int(r['end'])))
    iv.sort()
    n = 0
    with open(out + '.tmp', 'w') as w:
        cur = None
        for c, s, e in iv:
            if cur and cur[0] == c and s <= cur[2] + 100:
                cur[2] = max(cur[2], e)
            else:
                if cur:
                    w.write(f'{cur[0]}\t{cur[1]}\t{cur[2]}\n')
                    n += 1
                cur = [c, s, e]
        if cur:
            w.write(f'{cur[0]}\t{cur[1]}\t{cur[2]}\n')
            n += 1
    os.replace(out + '.tmp', out)
    print(tumor, label, select, 'loci not counted', len(iv), 'merged regions', n, '->', out)


if __name__ == '__main__':
    main()
