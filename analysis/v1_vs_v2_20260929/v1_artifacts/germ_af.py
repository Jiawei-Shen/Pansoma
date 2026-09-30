import json, pysam, numpy as np, collections
GVCF='/scratch/jshen/data/HG008_GIAB/dipcall_HG008N_GRCh38/HG008N_GRCh38_dipcall.dip.vcf.gz'
germ={}
for r in pysam.VariantFile(GVCF):
    if r.chrom!='chr1':
        if germ: break
        continue
    if len(r.ref)==1:
        gt=[s['GT'] for s in r.samples.values()][0]
        for i,a in enumerate(r.alts or [], start=1):
            if len(a)==1: germ[(r.pos,r.ref.upper(),a.upper())]='hom' if gt.count(i)==2 else 'het'
P='/scratch/jshen/data/Pansoma/HG008_GIAB/'
bins=[0,0.1,0.2,0.35,0.65,0.9,1.01]
for d in ['AF_HPRC','AF_HPRC_HG008N_added']:
    af=collections.defaultdict(list)
    for l in open(P+d+'/5ch_training_data_SNV/val/SNV_chr1_5chan_tensor_dataset/variant_summary_classified.ndjson'):
        r=json.loads(l); k=(r['genomic_position'],r['v_ref'].upper(),r['v_alt'].upper())
        cls='true' if r['classification']=='true' else ('germ_'+germ[k] if k in germ else 'false_other')
        af[cls].append(r['alt_allele_frequency'])
    print(d)
    for cls,v in sorted(af.items()):
        v=np.array(v); h=np.histogram(v,bins=bins)[0]
        print(f'  {cls:12s} n={len(v):7d} median={np.median(v):.3f} hist[0-.1,.1-.2,.2-.35,.35-.65,.65-.9,.9-1]={list(h)}')
