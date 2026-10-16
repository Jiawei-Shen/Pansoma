"""E7 step 0 (v2.1 only): are the genotypes of the decomposed (wave) VCF complete, so that s1 can skip the v1.1 repair?

For every top-level (LV=0) record of the raw deconstruct VCF in REGION with REF <= 100 kb (the bubbles vcfbub -l 0
-a 100000 keeps whole) and every haplotype with a called allele:
  expected non-ref = the raw allele is not 0;
  observed non-ref = the haplotype carries a non-0 allele in at least one wave record that comes from this raw record
                     (wave ID = raw ID, or raw ID + '_<n>' after vcfwave; IDs of records merged by
                     `bcftools norm -m +any` are joined with ';').
A wave VCF with incomplete genotypes (the v1.1 bug, see s1) shows 'missing' = expected but not observed. Run it on v1.1
too (E7_GRAPH=v1.1) as the positive control.
Usage: E7_GRAPH=v2.1 s0_check_wave_gt.py RAW_VCF REGION
"""
import collections, json, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import BCFTOOLS, GRAPH, HPRC_VCF

raw_vcf, region = sys.argv[1], sys.argv[2]


def query(vcf, fmt):
    p = subprocess.Popen([BCFTOOLS, 'query', '-r', region, '-f', fmt, vcf], stdout=subprocess.PIPE, text=True,
                         bufsize=1 << 20)
    for line in p.stdout:
        yield line.rstrip('\n').split('\t')
    assert p.wait() == 0, (vcf, region)


def haps(gts):
    """[(haplotype index, allele string)] for the GT strings (column 0 = CHM13, haploid)."""
    out = []
    for i, g in enumerate(gts):
        a = g.replace('/', '|').split('|')
        bits = [0] if i == 0 else [2 * i - 1, 2 * i]
        out += list(zip(bits, a))
    return out


# wave: raw ID -> set of haplotypes non-ref in any wave record of that bubble; obs1 uses only the records with a single
# ID (a merged record cannot tell which of its ALT alleles came from which bubble, so it is used only for 'missing')
obs, obs1 = collections.defaultdict(set), collections.defaultdict(set)
nwave = nmerged = 0
for f in query(HPRC_VCF, '%ID[\t%GT]\n'):
    nwave += 1
    nonref = {h for h, a in haps(f[1:]) if a not in ('.', '0')}
    ids = f[0].split(';')
    nmerged += len(ids) > 1
    for x in ids:
        base, _, suf = x.rpartition('_')
        rid = base if (base and suf.isdigit()) else x
        obs[rid] |= nonref
        if len(ids) == 1:
            obs1[rid] |= nonref
examples = []

c = collections.Counter()
for f in query(raw_vcf, '%ID\t%REF\t%INFO/LV[\t%GT]\n'):
    rid, ref, lv, gts = f[0], f[1], f[2], f[3:]
    if lv != '0' or len(ref) > 100000:
        continue
    c['raw_records'] += 1
    o = obs.get(rid)
    if o is None:
        c['raw_records_without_wave'] += 1
        continue
    for h, a in haps(gts):
        if a == '.':
            continue
        e = a != '0'
        seen = h in o
        c['hap_calls'] += 1
        c['expected_nonref'] += e
        c['missing (expected, not in wave)'] += e and not seen
        c['spurious (raw ref, non-ref in a single-ID wave record)'] += (not e) and h in obs1.get(rid, ())
        if e and not seen and len(examples) < 5:
            examples.append((rid, len(ref), h, a))
c['wave_records'] = nwave
c['wave_records_merged'] = nmerged
res = dict(graph=GRAPH, region=region, **c)
res['missing_rate'] = round(c['missing (expected, not in wave)'] / max(1, c['expected_nonref']), 6)
res['missing_examples'] = examples
print(json.dumps(res))
