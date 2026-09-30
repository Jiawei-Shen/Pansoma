import json, gzip, collections
exec(open('curve_head.py').read())
from pon_check import pon_hit
P='/scratch/jshen/data/pansoma_net_v2_runs/HG008_Illumina_SNV_base/test_chr1/Liss_lab_BCM_Illumina-WGS_20240313.v3_tensors.SNV.predictions.ndjson.gz'
L='/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_BCM_Illumina-WGS_20240313/v3_tensors/SNV/chr1_labels.ndjson'
calls={}; partial=set(); nin=0; nin_mine=0; mism=0
for lp, ll in zip(gzip.open(P,'rt'), open(L)):
    p=json.loads(lp); l=json.loads(ll)
    if p['candidate_id']!=l['candidate_id']: mism+=1; continue
    g=l.get('grch38')
    if not g: continue
    k=(g['pos0']+1, g['ref'].upper(), g['alt'].upper())
    if p['in_test']: nin+=1
    if not (k in truth or inbed(k[0])): continue
    nin_mine+=1
    s=p['p_somatic']
    if s>calls.get(k,-1): calls[k]=s
    if p['truth_ids'] and k not in truth: partial.add(k)
print('mismatch', mism, 'v2 in_test tensors', nin, 'tensors in my BED or TP', nin_mine, 'distinct alleles', len(calls), 'partial-indel alleles', len(partial))
items=sorted(calls.items(), key=lambda x:-x[1])
TOP=30000
pon={}
for k,s in items[:TOP]:
    pon[k]=pon_hit(k[0],k[1],k[2],'strict')
def curve(use_pon, excl_partial=False):
    tps=set(); fp=0; pts=[]
    for k,s in items[:TOP]:
        if use_pon and pon[k]: continue
        if k in truth: tps.add(truth[k])
        elif excl_partial and k in partial: pass
        else: fp+=1
        pts.append((s,len(tps),fp))
    # collapse to thresholds
    last={}
    for s,tp,f in pts: last[s]=(tp,f)
    out=[(s,)+last[s] for s in sorted(last, reverse=True)]
    return out
def summarize(name, pts):
    best=max(pts, key=lambda x: (2*x[1]/(x[1]+x[2])*x[1]/697/((x[1]/(x[1]+x[2]))+x[1]/697)) if x[1] else 0)
    def PR(x): return x[1]/(x[1]+x[2]), x[1]/697
    b=PR(best); f1=2*b[0]*b[1]/(b[0]+b[1])
    s=f'{name}: maxF1={f1:.3f} (P={b[0]:.3f} R={b[1]:.3f} thr={best[0]:.4f} TP={best[1]} FP={best[2]})'
    for pr in [0.05,0.10,0.15,0.20,0.30,0.50]:
        c=[PR(x)[1] for x in pts if PR(x)[0]>=pr]; s+=f' | R@P>={pr}: {max(c):.3f}' if c else f' | R@P>={pr}: NA'
    for rc in [0.5,0.7,0.8,0.87,0.9]:
        c=[PR(x)[0] for x in pts if PR(x)[1]>=rc]; s+=f' | P@R>={rc}: {max(c):.3f}' if c else f' | P@R>={rc}: NA'
    s+=f' | last: thr={pts[-1][0]:.4f} TP={pts[-1][1]} FP={pts[-1][2]}'
    print(s)
for nm,up,ex in [('v2_base_noPoN',False,False),('v2_base_noPoN_exclPartial',False,True),('v2_base_4PoN',True,False),('v2_base_4PoN_exclPartial',True,True)]:
    summarize(nm, curve(up,ex))
# PoN effect on TPs/FPs at the threshold giving v2 recall ~0.80/0.87 pre-PoN
pre=curve(False)
for target in [0.805,0.872]:
    x=[p for p in pre if p[1]/697>=target][0]; thr=x[0]
    sub=[(k,s) for k,s in items[:TOP] if s>=thr]
    tp=len({truth[k] for k,s in sub if k in truth}); fp=sum(1 for k,s in sub if k not in truth)
    tp2=len({truth[k] for k,s in sub if k in truth and not pon[k]}); fp2=sum(1 for k,s in sub if k not in truth and not pon[k])
    germ_fp=sum(1 for k,s in sub if k not in truth and k in germ)
    print(f'v2 at thr {thr:.4f}: pre-PoN TP={tp} FP={fp} (HG008N germline alleles {germ_fp}) P={tp/(tp+fp):.3f} R={tp/697:.3f} | post-PoN TP={tp2} FP={fp2} P={tp2/(tp2+fp2):.3f} R={tp2/697:.3f}')
