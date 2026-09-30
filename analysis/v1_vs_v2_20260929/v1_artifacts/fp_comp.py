import json, gzip, collections, pysam, numpy as np
exec(open('curve_head.py').read())
from pon_check import pon_hit
P='/scratch/jshen/data/pansoma_net_v2_runs/HG008_Illumina_SNV_base/test_chr1/Liss_lab_BCM_Illumina-WGS_20240313.v3_tensors.SNV.predictions.ndjson.gz'
L='/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_BCM_Illumina-WGS_20240313/v3_tensors/SNV/chr1_labels.ndjson'
VS='/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_BCM_Illumina-WGS_20240313/v3_tensors/SNV/chr1_variant_summary.ndjson'
af={}
for l in open(VS):
    r=json.loads(l); af[r['candidate_id']]=(r.get('af'), r.get('coverage'))
def binaf(a):
    for lo,hi in [(0,0.1),(0.1,0.2),(0.2,0.35),(0.35,0.65),(0.65,0.9),(0.9,1.01)]:
        if lo<=a<hi: return f'{lo}-{hi}'
v2=[]
for lp, ll in zip(gzip.open(P,'rt'), open(L)):
    p=json.loads(lp)
    if not p['in_test'] or p['p_somatic']<0.1388: continue
    l=json.loads(ll); g=l['grch38']
    if not g: continue
    k=(g['pos0']+1,g['ref'].upper(),g['alt'].upper())
    if k in truth or not inbed(k[0]): continue
    v2.append((k, af[p['candidate_id']][0], af[p['candidate_id']][1]))
v1=[]
for r in pysam.VariantFile(VCFS['e053_linear_May14']).fetch('chr1'):
    k=(r.pos,r.ref.upper(),r.alts[0].upper())
    if k in truth or not inbed(k[0]): continue
    v1.append((k, float(r.info['AF']), int(r.info['DP'])))
v2k={x[0] for x in v2}; v1k={x[0] for x in v1}
print('pre-PoN FP alleles: v1', len(v1k), 'v2', len(v2k), 'shared', len(v1k&v2k))
for name, lst in [('v1',v1),('v2',v2)]:
    c=collections.Counter(); a=collections.Counter()
    for k,f,d in lst:
        pn=pon_hit(k[0],k[1],k[2],'strict'); gm=k in germ
        c[('PoN' if pn else 'noPoN', 'germ' if gm else 'nongerm')]+=1
        if not pn: a[('germ' if gm else 'nongerm', binaf(f))]+=1
    print(name, dict(c))
    print('   post-PoN FP by (germline?, AF bin):', sorted(a.items()))
# attribution: are v2's post-PoN non-germline FPs v1 candidates? what did v1 do with them?
v1cand={}
for l in open('/scratch/jshen/data/Pansoma/HG008_GIAB/AF_HPRC/5ch_training_data_SNV/val/SNV_chr1_5chan_tensor_dataset/variant_summary_classified.ndjson'):
    r=json.loads(l); v1cand[(r['genomic_position'],r['v_ref'].upper(),r['v_alt'].upper())]=(r['alt_allele_frequency'], r['coverage_at_locus'], r['alt_allele_count'])
v1call={}
for r in pysam.VariantFile(VCFS['e053_linear_May14']).fetch('chr1'):
    v1call[(r.pos,r.ref.upper(),r.alts[0].upper())]=float(r.info['PROB'])
grp=collections.Counter(); afpair=[]
for k,f,d in v2:
    if pon_hit(k[0],k[1],k[2],'strict') or k in germ: continue
    if k in v1cand:
        grp['v1_candidate_called' if k in v1call else 'v1_candidate_not_called(<0.5)']+=1
        afpair.append((f, v1cand[k][0]))
    else: grp['not_a_v1_candidate']+=1
print('v2 post-PoN non-germline FPs vs v1:', dict(grp))
if afpair:
    a=np.array(afpair); print('  v2 af median', np.median(a[:,0]), 'v1 af median (same sites)', np.median(a[:,1]))
# v2 TP recall by v1-candidate AF for truths: same 608?
