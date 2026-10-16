"""Reference point for the graph PoN's false-tag rate: does the repository's GRCh38 PoN (scripts/filter_panel_of_normals.py,
unmodified, run as a CLI) tag the GRCh38 truth records of the same gold events?

Every gold event (assembly_offref_20260929/offref_somatic.tsv) lists its HG008-T truth records (grch38_records,
chrom:pos:REF>ALT, VCF coordinates). They are written to a small VCF, tagged with the four PoNs, and an event counts
as tagged when any of its records is tagged.
Output: DATA/gold_grch38_pon/{gold_truth.vcf.gz, gold_truth.pon.vcf.gz, events.tsv}
"""
import csv
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e9c_lib as L  # noqa: E402

PON_DIR = '/scratch/jshen/data/Pansoma/panel_of_normal_VCFs'
PONS = [f'{PON_DIR}/{n}' for n in ('af-only-gnomad.hg38.vcf.gz', 'Homo_sapiens_assembly38.dbsnp138.vcf.gz',
                                   '1000g_pon.hg38.vcf.gz', 'CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz')]
PY = '/wanglab/jshen/anaconda3/bin/python'


def records(field):
    """'chrom:pos:REF>ALT,...' -> [(chrom, pos, ref, alt)]; a piece without ':' is a further ALT of the previous record."""
    out = []
    for x in field.split(','):
        if ':' in x:
            chrom, pos, allele = x.split(':')
            ref, alt = allele.split('>')
            out.append((chrom, int(pos), ref, alt))
        elif out:
            c, p, r, a = out[-1]
            out.append((c, p, r, x))
    return out


def rid(t):
    return f'{t[0]}:{t[1]}:{t[2]}>{t[3]}'


def main():
    out = f'{L.DATA}/gold_grch38_pon'
    os.makedirs(out, exist_ok=True)
    gold = list(csv.DictReader(open(L.GOLD), delimiter='\t'))
    recs = set()
    for r in gold:
        recs.update(records(r['grch38_records']))
    order = {f'chr{i}': i for i in range(1, 23)}
    order.update(chrX=23, chrY=24)
    import pysam
    fa = pysam.FastaFile(L.FASTA)
    with open(f'{out}/gold_truth.vcf', 'w') as w:
        w.write('##fileformat=VCFv4.2\n')
        for c in sorted(order, key=order.get):
            w.write(f'##contig=<ID={c},length={fa.get_reference_length(c)}>\n')
        w.write(f'##reference=file://{L.FASTA}\n')
        w.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for chrom, pos, ref, alt in sorted(recs, key=lambda t: (order.get(t[0], 99), t[1], t[2], t[3])):
            w.write(f'{chrom}\t{pos}\t{chrom}:{pos}:{ref}>{alt}\t{ref}\t{alt}\t.\tPASS\t.\n')
    subprocess.run(f'module load samtools/1.18 2>/dev/null; bgzip -f {out}/gold_truth.vcf && tabix -f -p vcf {out}/gold_truth.vcf.gz',
                   shell=True, check=True, executable='/bin/bash')
    subprocess.run([PY, f'{L.REPO}/scripts/filter_panel_of_normals.py', f'{out}/gold_truth.vcf.gz', f'{out}/gold_truth.pon.vcf.gz',
                    '--pon', *PONS, '--fasta', L.FASTA], check=True)
    tagged = {}
    for f in L.vcf_records(f'{out}/gold_truth.pon.vcf.gz'):
        tagged[f[2]] = 'PanelOfNormals' in f[6]
    with open(f'{out}/events.tsv', 'w') as w:
        w.write('event\ttensor\ttensor_exact\tchrom\tgrch38_records\tgrch38_pon_tagged\n')
        for r in gold:
            t = any(tagged.get(rid(x), False) for x in records(r['grch38_records']))
            w.write(f"{r['event']}\t{r['tensor']}\t{r['tensor_exact']}\t{r['contig'].split('_')[0]}\t{r['grch38_records']}\t{int(t)}\n")
    n = len(gold)
    ex = [r for r in gold if r['tensor_exact'] == 'True']
    hit = lambda r: any(tagged.get(rid(x), False) for x in records(r['grch38_records']))  # noqa: E731
    print(f'gold events {n}: GRCh38 PoN tags {sum(map(hit, gold))}; exact-key events {len(ex)}: tagged {sum(map(hit, ex))}')


if __name__ == '__main__':
    main()
