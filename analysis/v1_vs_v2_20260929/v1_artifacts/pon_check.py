import pysam, sys
PD='/scratch/jshen/data/Panels_of_Normal/'
PONS=[pysam.VariantFile(PD+f) for f in ['1000g_pon.hg38.vcf.gz','af-only-gnomad.hg38.vcf.gz','CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.vcf.gz','Homo_sapiens_assembly38.dbsnp138.vcf.gz']]
def pon_hit(pos, ref, alt, mode='any'):
    for v in PONS:
        for r in v.fetch('chr1', pos-1, pos):
            if r.pos!=pos: continue
            if mode=='any':
                if r.ref.upper()==ref and alt in [a.upper() for a in (r.alts or [])]: return True
            else:
                if r.ref.upper()==ref and [a.upper() for a in (r.alts or [])]==[alt]: return True
    return False
if __name__=='__main__':
    f='/scratch/jshen/Pansoma_testing_results_V2/HG008_GIAB/AF_HPRC/pansoma_HG008T_WGS_ALL_chr/pansoma_HG008T_WGS_chr1/pansoma-to_SNV_pansoma_HG008T_WGS_chr1.linear.vcf.gz'
    keep={'any':0,'strict':0}; n=0
    for r in pysam.VariantFile(f).fetch('chr1'):
        n+=1
        for m in keep:
            if not pon_hit(r.pos, r.ref.upper(), r.alts[0].upper(), m): keep[m]+=1
    print('n',n,'kept',keep, '(bcftools result: 1114)')
