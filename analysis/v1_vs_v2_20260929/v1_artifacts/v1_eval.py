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
# v1 val candidates
cand=collections.Counter(); cand_ids=set(); valkeys=set(); lab=collections.Counter()
val_inbed=collections.Counter()
for l in open(VAL):
    r=json.loads(l)
    valkeys.add((r['node_id'],r['variant_key']))
    gp=r['genomic_position']; k=(gp,r['v_ref'].upper(),r['v_alt'].upper())
    ib=inbed(gp); tid=truth.get(k)
    if tid is not None: cand_ids.add(tid)
    g = k in germ
    val_inbed[(r['classification'], 'inBED' if ib else 'outBED', 'truth697' if tid is not None else ('germline' if g else 'other'))]+=1
print('v1 val tensors by (label, bed, match):'); [print('  ',k,v) for k,v in sorted(val_inbed.items())]
print('v1 ceiling (697 with a v1 val/ref-path candidate):', len(cand_ids), len(cand_ids)/697)
print('  v2 status of the v1-found:', collections.Counter(v2status[t] for t in cand_ids))
print('  v2 status of the v1-missed:', collections.Counter(v2status[t] for t in set(v2status)-cand_ids))
print('  v1-missed but v2 has tensor:', sorted(int([r for r in t697 if r['truth_id']==t][0]['vcf_pos']) for t in set(v2status)-cand_ids if v2status[t]=='tensor_representative')[:40])
print('  v2-missed but v1 has cand:', sorted(int([r for r in t697 if r['truth_id']==t][0]['vcf_pos']) for t in cand_ids if v2status[t]!='tensor_representative'))
# test set vs val
testkeys=set(); ntest=0
for l in open(TEST):
    r=json.loads(l); ntest+=1; testkeys.add((r['node_id'],r['variant_key']))
print('test n', ntest, 'distinct', len(testkeys), 'val n', len(valkeys), 'val-test', len(valkeys-testkeys), 'test-val', len(testkeys-valkeys))
json.dump(sorted(cand_ids, key=int), open('/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/v1_artifacts/v1_candidate_truth_ids.json','w'))
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
res={}
for name,path in VCFS.items():
    try:
        calls,items=evalvcf(path)
    except Exception as e:
        print(name,'ERR',e); continue
    items.sort(key=lambda x:-x[0])
    tps=set(); fp=0; fpg=0; curve=[]
    for p,x,t in items:
        if t=='TP': tps.add(x)
        else:
            fp+=1
            if t=='FP_germ': fpg+=1
        curve.append((p,len(tps),fp))
    ntp=len(tps); P=ntp/max(1,ntp+fp); R=ntp/697
    cnt=collections.Counter(t for _,_,t in items)
    print(f'{name}: calls={len(calls)} evaluated_inBED_or_TP={len(items)} TP={ntp} FP_inBED={fp} (of which HG008N-germline allele {fpg}, somatic-other-kind pos {cnt["FP_somother"]}) P={P:.4f} R={R:.4f} F1={2*P*R/max(1e-9,P+R):.4f} minPROB={min(p for p,_,_ in items):.3f}')
    res[name]=curve
# curve points for e053 May14
def points(curve):
    # thresholds unique
    out={}
    best=(0,None)
    pts=[]
    last={}
    for p,tp,fp in curve: last[p]=(tp,fp)
    ths=sorted(last, reverse=True)
    for th in ths:
        tp,fp=last[th]; P=tp/(tp+fp); R=tp/697; F=2*P*R/(P+R) if P+R>0 else 0
        pts.append((th,tp,fp,P,R,F))
        if F>best[0]: best=(F,(th,tp,fp,P,R))
    out['maxF1']=best
    for pr in [0.05,0.10,0.15,0.20]:
        c=[x for x in pts if x[3]>=pr]; out[f'recall@P>={pr}']=max(c,key=lambda x:x[4])[4] if c else None
    for rc in [0.5,0.7,0.8,0.9]:
        c=[x for x in pts if x[4]>=rc]; out[f'P@R>={rc}']=max(c,key=lambda x:x[3])[3] if c else None
    out['at_min_threshold']=pts[-1]
    return out, pts
for name in ['e053_linear_May14','e053_linear_PoN_Apr9','run1_linear_Mar9','run1_linear_PoN_Mar9','HG008Nadded_linear_PoN']:
    if name in res:
        o,pts=points(res[name]); print(name, json.dumps(o, default=str))
        if name=='e053_linear_May14':
            for th in [0.5,0.6,0.7,0.8,0.9,0.95,0.99]:
                c=[x for x in pts if x[0]>=th]
                if c: x=c[-1]; print(f'   thr>={th}: TP={x[1]} FP={x[2]} P={x[3]:.4f} R={x[4]:.4f} F1={x[5]:.4f}')
