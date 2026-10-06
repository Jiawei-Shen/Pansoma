"""Step 5b: population metadata of the HPRC v1.1 samples (the 44 sample columns of the d9 deconstruct VCF + CHM13 + GRCh38).

Sources (downloaded once into $D/metadata_sources/, then read from there):
  HPRC  human-pangenomics/HPP_Year1_Data_Freeze_v1.0 sample_metadata/hprc_year1_sample_metadata.txt
        (Sex, Subpopulation = 1000G population code, Superpopulation, Notes)       -> primary source
  IGSR  internationalgenome.org api/beta/sample/<sample> (one JSON per sample: sex, population code / name /
        superpopulation)                                                            -> cross-check
  1000G phase3 20131219.populations.tsv                                            -> population names (description)
  Coriell catalogue page (Product=DNA)                                              -> cross-check of samples IGSR lacks
Rules: population_code / superpopulation from HPRC, checked against IGSR (agree column); population_name from the 1000G
table. NA21309 has no code / superpopulation in HPRC (Notes 'HAPMAP - MAASAI IN KINYAWA, KENYA') and is in neither IGSR
nor the Coriell page: code MKK (the HapMap 3 code of that population) and superpopulation AFR are written with a note
(the population lives in Kenya; not an explicit field of any source). Nothing else is filled in by hand.
The HG008 donor (GIAB) ancestry is recorded from the GIAB HG008 paper (McDaniel et al. 2025, Sci Data 12:1195, doi 10.1038/s41597-025-05438-2, PMC12267650).
Outputs: $D/hprc_sample_metadata.tsv, $D/hg008_donor_ancestry.tsv; prints the superpopulation composition.
"""
import collections, csv, html, json, os, re, urllib.request
import pysam
D = '/scratch/jshen/data/pansoma_net_v2_runs/graph_absorbed_somatic_20261001'
S = f'{D}/metadata_sources'
U_HPRC = 'https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0/main/sample_metadata/hprc_year1_sample_metadata.txt'
U_POP = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv'
U_IGSR = 'https://www.internationalgenome.org/api/beta/sample/{}'
U_COR = 'https://www.coriell.org/0/Sections/Search/Sample_Detail.aspx?Ref={}&Product=DNA'
U_HG008 = 'https://pmc.ncbi.nlm.nih.gov/articles/PMC12267650/'
SUPER = {'AFR': 'African', 'AMR': 'Admixed American', 'EAS': 'East Asian', 'SAS': 'South Asian', 'EUR': 'European'}
os.makedirs(S, exist_ok=True)


def get(url, name):
    p = f'{S}/{name}'
    if not os.path.exists(p):
        try:
            urllib.request.urlretrieve(url, p)
        except urllib.error.HTTPError as e:           # IGSR answers 404 for a sample it does not have
            open(p, 'w').write('{}' if e.code == 404 and name.endswith('.json') else '')
    return open(p, errors='ignore').read()


samples = list(pysam.VariantFile('/scratch/jshen/data/AF-Filtered_VG_Indexes/hprc-v1.1-mc-grch38.d9.vcf.gz').header.samples)
hprc = {r['Sample']: r for r in csv.DictReader(get(U_HPRC, 'hprc_year1_sample_metadata.txt').splitlines(), delimiter='\t')}
popname = {r['Population Code']: r['Population Description'] for r in
           csv.DictReader(get(U_POP, '1000G_20131219.populations.tsv').splitlines(), delimiter='\t') if r['Population Code']}
rows = [dict(sample='GRCh38', population_code='', population_name='reference', superpopulation='reference', sex='',
             source_url='', igsr_check='', note='graph reference path; allele 0 of every VCF record, no GT column')]
for s in samples:
    if s == 'CHM13':
        rows.append(dict(sample=s, population_code='', population_name='reference (T2T CHM13, hydatidiform mole)',
                         superpopulation='reference', sex='', source_url='', igsr_check='',
                         note='haploid (one haplotype, GT has one allele)'))
        continue
    h = hprc.get(s)
    if h is None:
        rows.append(dict(sample=s, note='NOT FOUND in HPRC metadata')); continue
    j = json.loads(get(U_IGSR.format(s), f'igsr_{s}.json')); j = j.get('_source', j)
    ip = j.get('populations') or []
    code, sup, note = h['Subpopulation'], h['Superpopulation'], h['Notes']
    if ip:
        ok = ip[0]['code'] == code and ip[0].get('superpopulationCode') == sup and (j.get('sex') or '?')[0] == h['Sex'][0].upper()
        check = f"IGSR {ip[0]['code']}/{ip[0].get('superpopulationCode')}/{j.get('sex')}: {'agree' if ok else 'DISAGREE'}"
        url = f'{U_HPRC} ; {U_IGSR.format(s)}'
    else:
        t = get(U_COR.format(s), f'coriell_{s}.html')
        t = re.sub(r'\s+', ' ', html.unescape(re.sub(r'<[^>]+>', ' ', re.sub(r'<(script|style).*?</\1>', '', t, flags=re.S))))
        m = re.search(r'Description: (.*?) Affected:.*?Sex: (\w+)', t)
        check = f'not in IGSR; Coriell: {m.group(1)} / {m.group(2)}' if m else 'not in IGSR; not on the Coriell page'
        url = f'{U_HPRC} ; {U_COR.format(s)}' if m else U_HPRC
    if not code and 'MAASAI IN KINYAWA' in note.upper():
        code, sup = 'MKK', 'AFR'; popname['MKK'] = 'Maasai in Kinyawa, Kenya (HapMap 3)'
        note += '; code MKK = HapMap 3 Maasai in Kinyawa, Kenya; superpopulation AFR from the population location (not a source field)'
    rows.append(dict(sample=s, population_code=code, population_name=popname.get(code, ''),
                     superpopulation=sup, superpopulation_name=SUPER.get(sup, ''), sex=h['Sex'], cohort=h['Cohort'],
                     family_id=h['FamilyID'], source_url=url, igsr_check=check, note=note))
cols = ['sample', 'population_code', 'population_name', 'superpopulation', 'superpopulation_name', 'sex', 'cohort', 'family_id',
        'source_url', 'igsr_check', 'note']
with open(f'{D}/hprc_sample_metadata.tsv', 'w') as w:
    w.write('\t'.join(cols) + '\n')
    for r in rows:
        w.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
with open(f'{D}/hg008_donor_ancestry.tsv', 'w') as w:
    w.write('sample\treported_ancestry\tsex\tsource\tsource_url\n')
    w.write('HG008\tEuropean genetic ancestry (paper: "derived from a female of European genetic ancestry")\tfemale, 61 years, XX\t'
            'McDaniel et al. 2025 Sci Data 12:1195 (doi 10.1038/s41597-025-05438-2) "Development and extensive sequencing of a broadly-consented Genome in a Bottle '
            f'matched tumor-normal pair" (GIAB HG008, PDAC)\t{U_HG008}\n')
hs = [r for r in rows if r['superpopulation'] != 'reference']
print(len(hs), 'samples', len(rows), 'rows;', sum('DISAGREE' in r.get('igsr_check', '') for r in hs), 'disagree with IGSR;',
      sum('not in IGSR' in r.get('igsr_check', '') for r in hs), 'not in IGSR')
print('superpopulation', dict(collections.Counter(r['superpopulation'] for r in hs)))
print('population', dict(collections.Counter(r['population_code'] for r in hs)))
