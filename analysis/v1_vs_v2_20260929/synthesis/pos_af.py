"""AF of positive (somatic) training labels: v2 chr2-22 (exact vs residual partial) and v1 AF_HPRC chr1 val true. Read-only."""
import json, pickle, collections, time, numpy as np
T0=time.time()
V2='/scratch/jshen/data/pansoma_v2_tensors/HG008T_Illumina/tensors/SNV/'
pos={}
for c in [f'chr{i}' for i in range(2,23)]:
    for l in open(V2+f'{c}_labels.ndjson'):
        if '"label": 1,' not in l: continue
        r=json.loads(l)
        if r['label']==1: pos[r['candidate_id'],c]=r['reason']
print('v2 chr2-22 label-1 tensors',len(pos),collections.Counter(pos.values()).most_common(6),'t=%.0f'%(time.time()-T0))
af=collections.defaultdict(list)
for c in [f'chr{i}' for i in range(2,23)]:
    for l in open(V2+f'{c}_variant_summary.ndjson'):
        cid=l[18:l.index("\"",18)]
        k=(cid,c)
        if k in pos:
            r=json.loads(l); af['partial' if 'partial' in pos[k] else 'exact'].append(r['af'])
for k,v in af.items():
    v=np.array(v); print('v2',k,len(v),'median %.3f'%np.median(v),'AF<0.2: %d (%.1f%%)'%((v<0.2).sum(),100*(v<0.2).mean()))
allv=np.concatenate([np.array(v) for v in af.values()]); lo=(allv<0.2).sum()
print('v2 all positives',len(allv),'AF<0.2',lo, 'of which partial', (np.array(af['partial'])<0.2).sum())
d=pickle.load(open('../head_to_head/v1_val_recs.pkl','rb'))
t=np.array([r[7] for r in d if r[2]==1]); print('v1 chr1 val true',len(t),'median %.3f'%np.median(t),'AF<0.2 %d (%.1f%%)'%((t<0.2).sum(),100*(t<0.2).mean()))
print('t=%.0f'%(time.time()-T0))
