import json, gzip, collections
exec(open('curve_head.py').read())
P='/scratch/jshen/data/pansoma_net_v2_runs/HG008_Illumina_SNV_base/test_chr1/Liss_lab_BCM_Illumina-WGS_20240313.v3_tensors.SNV.predictions.ndjson.gz'
L='/scratch/jshen/data/pansoma_v2_tensors/Liss_lab_BCM_Illumina-WGS_20240313/v3_tensors/SNV/chr1_labels.ndjson'
c=collections.Counter()
for lp, ll in zip(gzip.open(P,'rt'), open(L)):
    p=json.loads(lp); l=json.loads(ll); g=l.get('grch38')
    ib = bool(g) and inbed(g['pos0']+1)
    c[(p['in_test'], 'myBED' if ib else ('offref' if not g else 'outBED'), p['reason'])]+=1
for k,v in sorted(c.items(), key=lambda x:(-x[1])): print(v,k)
