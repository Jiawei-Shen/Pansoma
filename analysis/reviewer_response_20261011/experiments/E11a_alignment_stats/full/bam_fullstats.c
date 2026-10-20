/* E11a full-data scan of one contig of an indexed GRCh38 BAM: per-read alignment counts, the counterpart of
gam_fullstats for the linear alignments.

usage: bam_fullstats BAM CONTIG
Records: primary only (flags 0x100 secondary and 0x800 supplementary skipped and counted); one record = one read or
one mate; duplicates (0x400) are kept, as in the GAMs, and counted.
Per mapped primary record (flag 0x4 records placed at their mate are counted as unmapped):
  substitution  NM - inserted - deleted bases > 0 (NM counts mismatches, inserted and deleted bases)
  indel      a CIGAR I or D operation
  soft clip  a CIGAR S or H operation
  perfect    mapped, no substitution, no indel, no soft clip
Output (stdout, TSV): "CONTIG records perfect sub indel clip read_bases aligned_bases sub_bases ins_bases del_bases
indel_events clip_bases", then "#other unmapped_placed secondary supplementary duplicates no_nm".
*/
#include <stdio.h>
#include <stdlib.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/sam.h>

extern int hfile_set_blksize(hFILE *fp, size_t bufsiz); /* exported by htslib, declared only inside hfile.h's hwrite */

int main(int argc, char **argv) {
    if (argc != 3) { fprintf(stderr, "usage: bam_fullstats BAM CONTIG\n"); return 2; }
    samFile *in = sam_open(argv[1], "r");
    if (!in) { fprintf(stderr, "bam_fullstats: open %s\n", argv[1]); return 2; }
    if (in->format.compression == bgzf) hfile_set_blksize(in->fp.bgzf->fp, 16 << 20); /* 16 MB reads: small reads stall on a busy BeeGFS */
    hts_set_threads(in, 2);
    sam_hdr_t *h = sam_hdr_read(in);
    hts_idx_t *idx = sam_index_load(in, argv[1]);
    if (!h || !idx) { fprintf(stderr, "bam_fullstats: header or index\n"); return 2; }
    hts_itr_t *it = sam_itr_querys(idx, h, argv[2]);
    if (!it) { fprintf(stderr, "bam_fullstats: region %s\n", argv[2]); return 2; }
    int tid = it->tid;
    bam1_t *b = bam_init1();
    long long n[12] = {0}, unm = 0, sec = 0, sup = 0, dup = 0, nonm = 0;
    int r;
    while ((r = sam_itr_next(in, it, b)) >= 0) {
        uint16_t fl = b->core.flag;
        if (b->core.tid != tid) continue;
        if (fl & BAM_FSECONDARY) { sec++; continue; }
        if (fl & BAM_FSUPPLEMENTARY) { sup++; continue; }
        if (fl & BAM_FDUP) dup++;
        if (fl & BAM_FUNMAP) { unm++; continue; }
        uint32_t *cig = bam_get_cigar(b);
        long long ins = 0, del = 0, clip = 0, aln = 0, ev = 0, xb = 0;
        int has_eqx = 0;
        for (uint32_t i = 0; i < b->core.n_cigar; i++) {
            int op = bam_cigar_op(cig[i]);
            long long len = bam_cigar_oplen(cig[i]);
            switch (op) {
            case BAM_CMATCH: aln += len; break;
            case BAM_CEQUAL: aln += len; has_eqx = 1; break;
            case BAM_CDIFF: aln += len; xb += len; has_eqx = 1; break;
            case BAM_CINS: ins += len; ev++; break;
            case BAM_CDEL: del += len; ev++; break;
            case BAM_CSOFT_CLIP: case BAM_CHARD_CLIP: clip += len; break;
            default: break;
            }
        }
        long long sub;
        uint8_t *nm = bam_aux_get(b, "NM");
        if (nm) {
            sub = bam_aux2i(nm) - ins - del;
            if (sub < 0) sub = 0;
        } else {
            nonm++;
            sub = has_eqx ? xb : 0;
        }
        n[0]++;
        if (!sub && !ins && !del && !clip) n[1]++;
        if (sub) n[2]++;
        if (ins || del) n[3]++;
        if (clip) n[4]++;
        n[5] += b->core.l_qseq; n[6] += aln; n[7] += sub; n[8] += ins; n[9] += del; n[10] += ev; n[11] += clip;
    }
    if (r < -1) { fprintf(stderr, "bam_fullstats: read error %d in %s\n", r, argv[2]); return 2; }
    printf("%s", argv[2]);
    for (int k = 0; k < 12; k++) printf("\t%lld", n[k]);
    printf("\n#other\t%lld\t%lld\t%lld\t%lld\t%lld\n", unm, sec, sup, dup, nonm);
    bam_destroy1(b); hts_itr_destroy(it); hts_idx_destroy(idx); sam_hdr_destroy(h); sam_close(in);
    return 0;
}
