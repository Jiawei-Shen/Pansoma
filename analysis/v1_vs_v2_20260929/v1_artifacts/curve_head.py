import json, gzip, bisect, collections, sys
import numpy as np, pysam
V2='/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_BCM_Illumina-WGS_20240313/v3_tensors'
TRUTH='/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_tumorvariants.vcf.gz'
SBED='/scratch/jshen/data/HG008_GIAB/draft_v02_benchmark/HG008-T_somatic_smvar_benchmark_v0.2_all.bed'
GBED='/scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.bed'
GVCF='/scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.vcf.gz'
VAL='/scratch/jshen/data/Pansoma/HG008_GIAB/AF_HPRC/5ch_training_data_SNV/val/SNV_chr1_5chan_tensor_dataset/variant_summary_classified.ndjson'
TEST='/scratch/jshen/data/Pansoma/HG008_GIAB/AF_HPRC/5ch_testing_data_SNV/5ch_testing_data_SNV_chr1/variant_summary.ndjson'
R2='/scratch/jshen/Pansoma_testing_results_V2/HG008_GIAB/AF_HPRC'
VCFS={
 'e053_linear_May14': R2+'/pansoma_HG008T_WGS_ALL_chr/pansoma_HG008T_WGS_chr1/pansoma-to_SNV_pansoma_HG008T_WGS_chr1.linear.vcf.gz',
 'e053_linear_PoN_Apr9': R2+'/pansoma_HG008T_WGS_ALL_chr/pansoma_HG008T_WGS_chr1/pansoma-to_SNV_pansoma_HG008T_WGS_chr1.linear.filtered_PoN.vcf.gz',
 'e053_linear_Mar20': R2+'/pansoma_HG008T_WGS_Chr1/pansoma-to_SNV_pansoma_HG008T_WGS_Chr1.linear.vcf.gz',
 'e053_linear_PoN_Mar19': R2+'/pansoma_HG008T_WGS_Chr1/pansoma-to_SNV_pansoma_HG008T_WGS_Chr1.linear.filtered_PoN.vcf.gz',
 'run1_linear_Mar9': '/scratch/jshen/tmp/HG008T_WGS_AF-HPRC_pansoma.test.linear.vcf.gz',
 'run1_linear_PoN_Mar9': '/scratch/jshen/tmp/HG008T_WGS_AF-HPRC_pansoma.test.linear.filtered_PoN.vcf.gz',
 'HG008Nadded_linear_PoN': '/scratch/jshen/Pansoma_testing_results_V2/HG008_GIAB/AF_HPRC_HG008N_added/pansoma_HG008T_WGS_ALL_chr/pansoma_HG008T_WGS_chr1/pansoma-to_SNV_pansoma_HG008T_WGS_chr1.linear.filtered_PoN.vcf.gz',
}
def bed(path, chrom='chr1'):
    iv=[]
    op=gzip.open if path.endswith('.gz') else open
    for l in op(path,'rt'):
        if l.startswith(('#','track','browser')): continue
        f=l.split('\t')
        if f[0]==chrom: iv.append((int(f[1]),int(f[2])))
    iv.sort(); m=[]
    for s,e in iv:
        if m and s<=m[-1][1]: m[-1]=(m[-1][0],max(m[-1][1],e))
        else: m.append((s,e))
    return m
def intersect(a,b):
    i=j=0; out=[]
    while i<len(a) and j<len(b):
        s=max(a[i][0],b[j][0]); e=min(a[i][1],b[j][1])
        if s<e: out.append((s,e))
        if a[i][1]<b[j][1]: i+=1
        else: j+=1
    return out
conf=intersect(bed(SBED),bed(GBED))
starts=[s for s,e in conf]
def inbed(pos1):
    p=pos1-1; k=bisect.bisect_right(starts,p)-1
    return k>=0 and conf[k][0]<=p<conf[k][1]
print('confident bp chr1', sum(e-s for s,e in conf))
# 697 truths
rows=[l.rstrip('\n').split('\t') for l in open(V2+'/somatic.recall.tsv')]
hdr=rows[0]; rows=[dict(zip(hdr,r)) for r in rows[1:]]
t697=[r for r in rows if r['chrom']=='chr1' and r['kind']=='SNP' and r['passed']=='True' and r['in_bed']=='True']
print('n697', len(t697), collections.Counter(r['status'] for r in t697))
tv=pysam.VariantFile(TRUTH)
truth={}  # (pos,alt)->truth_id
allsom_pos=set()
for r in tv.fetch('chr1'):
    allsom_pos.add(r.pos)
pos697={int(r['vcf_pos']):r for r in t697}
for r in tv.fetch('chr1'):
    if r.pos in pos697 and len(r.ref)==1 and all(len(a)==1 for a in r.alts):
        for a in r.alts: truth[(r.pos,r.ref.upper(),a.upper())]=pos697[r.pos]['truth_id']
print('truth keys', len(truth), 'distinct ids', len(set(truth.values())), 'inbed check', sum(inbed(p) for p in pos697), 'not in my BED:', [p for p in pos697 if not inbed(p)])
v2status={r['truth_id']:r['status'] for r in t697}
# germline chr1 SNP alleles
gv=pysam.VariantFile(GVCF); germ=set()
for r in gv:
    if r.chrom!='chr1': continue
    if len(r.ref)==1:
        for a in r.alts or []:
            if len(a)==1: germ.add((r.pos,r.ref.upper(),a.upper()))
print('germline chr1 SNP alleles', len(germ))

# VCF evaluation
def evalvcf(path):
    calls={}
    for r in pysam.VariantFile(path).fetch('chr1') if path.endswith('.gz') else []:
        for a in r.alts:
            k=(r.pos,r.ref.upper(),a.upper())
            p=float(r.info.get('PROB',1.0))
            calls[k]=max(calls.get(k,0),p)
    items=[]
    for k,p in calls.items():
        tid=truth.get(k)
        if tid is not None: items.append((p,tid,'TP'))
        elif inbed(k[0]): items.append((p,k,'FP_germ' if k in germ else ('FP_somother' if k[0] in allsom_pos else 'FP')))
    return calls, items
