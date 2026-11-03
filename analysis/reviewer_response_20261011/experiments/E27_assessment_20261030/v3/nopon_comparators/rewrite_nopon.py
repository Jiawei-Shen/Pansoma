"""Turn a comparator VCF into its 'no population filtering' version: records whose only reason for failing is the tool's
own panel-of-normals / population-database filter become PASS; everything else is unchanged.
  ClairS-TO:   FILTER == NonSomatic, at least one PoN_* INFO flag and no Verdict_Germline flag -> PASS
               (NonSomatic from the verdict module, or combined with LowQual etc., stays filtered)
  DeepSomatic: FILTER == PON -> PASS (GERMLINE is the model's own class and stays)
usage: rewrite_nopon.py CALLER IN.vcf.gz OUT.vcf (plain text; bgzip afterwards)"""
import collections, gzip, sys
caller, src, dst = sys.argv[1:4]
n = collections.Counter()
with gzip.open(src, 'rt') as i, open(dst, 'w') as o:
    for line in i:
        if line[0] == '#':
            o.write(line); continue
        f = line.rstrip('\n').split('\t')
        flt, info = f[6], f[7].split(';')
        if caller == 'clairs_to' and flt == 'NonSomatic':
            pon = any(x.startswith('PoN_') for x in info); verdict = 'Verdict_Germline' in info
            if pon and not verdict:
                f[6] = 'PASS'; n['NonSomatic(PoN only)->PASS'] += 1
            else:
                n['NonSomatic kept (verdict)' if verdict else 'NonSomatic kept (no PoN flag)'] += 1
        elif caller == 'deepsomatic' and flt == 'PON':
            f[6] = 'PASS'; n['PON->PASS'] += 1
        n['records'] += 1
        o.write('\t'.join(f) + '\n')
print(src, dict(n), flush=True)
