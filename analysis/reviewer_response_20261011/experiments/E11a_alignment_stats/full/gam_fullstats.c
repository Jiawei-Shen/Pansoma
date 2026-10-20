/* E11a full-data scan of one segment of a sorted GAM: per-read alignment counts by GRCh38 placement contig.

usage: gam_fullstats GAM CHROM_NPY NCONTIG START END
  GAM        sorted GAM (BGZF; groups of [count, "GAM" tag, count - 1 Alignment messages])
  CHROM_NPY  the pipeline's grch38_path/chrom.npy (int16 per node ID; -1 = node off the GRCh38 path)
  START END  virtual offsets of the segment (group starts; END normalized, -1 = to EOF)
Per Alignment (one record = one read or one mate), the same rules as scan_gam.read_stats:
  mapped     the path has >= 1 mapping
  substitution  an edit with from_length == to_length > 0 and a sequence, or a complex edit (from != to, both > 0)
  indel      an edit with from_length == 0 < to_length or to_length == 0 < from_length, or a complex edit; an
             insertion that is the first edit of the first mapping or the last edit of the last mapping is a soft clip
  soft clip  as above (not an indel)
  perfect    mapped, no substitution, no indel, no soft clip
  placement  the contig holding most anchored bases (sum of from_length + 1 over the mappings on GRCh38-path nodes);
             ties -> lower contig index; -1 = mapped but no GRCh38-path node
Output (stdout, TSV): one line per placement with records: contig index (-2 = unmapped, -1 = unplaced), records,
perfect, with substitution, with indel, with soft clip, read bases, aligned bases (from == to), substituted bases,
inserted bases, deleted bases, indel events, soft-clipped bases; then a line "#total records groups secondary".
*/
#include <fcntl.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>

extern int hfile_set_blksize(hFILE *fp, size_t bufsiz); /* exported by htslib, declared only inside hfile.h's hwrite */

#define NF 12
static int16_t *CH;
static int64_t NNODE;
static int NC;
static int64_t (*CNT)[NF]; /* [NC + 2][NF]: index 0 = unmapped, 1 = unplaced, 2 + c = contig c */

static void die(const char *m) { fprintf(stderr, "gam_fullstats: %s\n", m); exit(2); }

static int bvarint(BGZF *fp, uint64_t *v) {
    uint64_t x = 0;
    for (int s = 0; s < 64; s += 7) {
        int c = bgzf_getc(fp);
        if (c < 0) return c == -1 && s == 0 ? 0 : -1;
        x |= (uint64_t)(c & 0x7f) << s;
        if (!(c & 0x80)) { *v = x; return 1; }
    }
    return -1;
}

static inline uint64_t mvarint(const uint8_t **p, const uint8_t *e) {
    uint64_t x = 0;
    for (int s = 0; s < 64 && *p < e; s += 7) {
        uint8_t c = *(*p)++;
        x |= (uint64_t)(c & 0x7f) << s;
        if (!(c & 0x80)) return x;
    }
    die("bad varint");
    return 0;
}

static inline void skip(const uint8_t **p, const uint8_t *e, int wt) {
    if (wt == 0) mvarint(p, e);
    else if (wt == 1) *p += 8;
    else if (wt == 5) *p += 4;
    else if (wt == 2) { uint64_t n = mvarint(p, e); *p += n; }
    else die("bad wire type");
    if (*p > e) die("field past message end");
}

typedef struct { uint32_t f, t; uint8_t seq; } Edit;
static Edit *ED;
static size_t EDCAP;
static int64_t *TALLY;
static int *TOUCH;

static void edit_push(size_t k, Edit e) {
    if (k >= EDCAP) { EDCAP = EDCAP ? EDCAP * 2 : 4096; ED = realloc(ED, EDCAP * sizeof(Edit)); if (!ED) die("oom"); }
    ED[k] = e;
}

static void alignment(const uint8_t *p, const uint8_t *e, int64_t *secondary) {
    const uint8_t *path = NULL, *pathe = NULL;
    uint64_t read_len = 0;
    while (p < e) {
        uint64_t key = mvarint(&p, e);
        int f = key >> 3, wt = key & 7;
        if (wt == 2 && (f == 1 || f == 2)) {
            uint64_t n = mvarint(&p, e);
            if (f == 1) read_len = n;
            else { path = p; pathe = p + n; }
            p += n;
        } else if (wt == 0 && f == 15) {
            if (mvarint(&p, e)) (*secondary)++;
        } else skip(&p, e, wt);
    }
    /* path: mappings (field 2) -> position (1: node_id 1), edits (2: from 1, to 2, sequence 3) */
    size_t ne = 0, nmap = 0, first_map_edits = 0, last_map_edits = 0;
    int ntouch = 0;
    while (path && path < pathe) {
        uint64_t key = mvarint(&path, pathe);
        int f = key >> 3, wt = key & 7;
        if (!(f == 2 && wt == 2)) { skip(&path, pathe, wt); continue; }
        uint64_t n = mvarint(&path, pathe);
        const uint8_t *m = path, *me = path + n;
        path = me;
        int64_t node = 0;
        uint64_t fsum = 0;
        size_t ne0 = ne;
        while (m < me) {
            uint64_t k2 = mvarint(&m, me);
            int f2 = k2 >> 3, w2 = k2 & 7;
            if (f2 == 1 && w2 == 2) {
                uint64_t n2 = mvarint(&m, me);
                const uint8_t *q = m, *qe = m + n2;
                m = qe;
                while (q < qe) {
                    uint64_t k3 = mvarint(&q, qe);
                    if ((k3 >> 3) == 1 && (k3 & 7) == 0) node = (int64_t)mvarint(&q, qe);
                    else skip(&q, qe, k3 & 7);
                }
            } else if (f2 == 2 && w2 == 2) {
                uint64_t n2 = mvarint(&m, me);
                const uint8_t *q = m, *qe = m + n2;
                m = qe;
                Edit ed = {0, 0, 0};
                while (q < qe) {
                    uint64_t k3 = mvarint(&q, qe);
                    int f3 = k3 >> 3, w3 = k3 & 7;
                    if (f3 == 1 && w3 == 0) ed.f = mvarint(&q, qe);
                    else if (f3 == 2 && w3 == 0) ed.t = mvarint(&q, qe);
                    else if (f3 == 3 && w3 == 2) { uint64_t l = mvarint(&q, qe); ed.seq = l > 0; q += l; }
                    else skip(&q, qe, w3);
                }
                edit_push(ne++, ed);
                fsum += ed.f;
            } else skip(&m, me, w2);
        }
        if (nmap == 0) first_map_edits = ne - ne0;
        last_map_edits = ne - ne0;
        nmap++;
        int c = (node >= 0 && node < NNODE) ? CH[node] : -1;
        if (c >= 0 && c < NC) {
            if (!TALLY[c]) TOUCH[ntouch++] = c;
            TALLY[c] += (int64_t)fsum + 1;
        }
    }
    int64_t *row;
    if (!nmap) { row = CNT[0]; row[0]++; row[6] += read_len; return; }
    int best = -1;
    int64_t bv = -1;
    for (int i = 0; i < ntouch; i++) {
        int c = TOUCH[i];
        if (TALLY[c] > bv || (TALLY[c] == bv && c < best)) { bv = TALLY[c]; best = c; }
        TALLY[c] = 0;
    }
    row = CNT[best < 0 ? 1 : best + 2];
    uint64_t match = 0, sub = 0, ins = 0, del = 0, ev = 0, clip = 0;
    int prev = 0;
    size_t sc_l = 0, sc_r = 0;
    if (ne && first_map_edits && ED[0].f == 0 && ED[0].t > 0) sc_l = 1;
    if (ne && last_map_edits && ED[ne - 1].f == 0 && ED[ne - 1].t > 0 && !(nmap == 1 && ne == 1)) sc_r = 1;
    for (size_t i = 0; i < ne; i++) {
        Edit d = ED[i];
        if ((i == 0 && sc_l) || (i == ne - 1 && sc_r)) { clip += d.t; prev = 5; continue; }
        if (d.f == d.t) {
            if (d.seq) { sub += d.f; prev = 1; } else { match += d.f; prev = 0; }
        } else if (d.t == 0) {
            del += d.f; if (prev != 3) ev++; prev = 3;
        } else if (d.f == 0) {
            ins += d.t; if (prev != 2) ev++; prev = 2;
        } else {
            uint32_t lo = d.f < d.t ? d.f : d.t;
            sub += lo;
            if (d.f > d.t) del += d.f - d.t; else ins += d.t - d.f;
            ev++; prev = 4;
        }
    }
    row[0]++;
    if (!sub && !ins && !del && !clip) row[1]++;
    if (sub) row[2]++;
    if (ins || del) row[3]++;
    if (clip) row[4]++;
    row[5] += read_len; row[6] += match + sub; row[7] += sub; row[8] += ins; row[9] += del; row[10] += ev; row[11] += clip;
}

int main(int argc, char **argv) {
    if (argc != 6) die("usage: gam_fullstats GAM CHROM_NPY NCONTIG START END");
    NC = atoi(argv[3]);
    int64_t start = atoll(argv[4]), end = atoll(argv[5]);
    int fd = open(argv[2], O_RDONLY);
    if (fd < 0) die("open chrom.npy");
    struct stat st;
    fstat(fd, &st);
    uint8_t *npy = mmap(NULL, st.st_size, PROT_READ, MAP_SHARED, fd, 0);
    if (npy == MAP_FAILED || memcmp(npy, "\x93NUMPY", 6)) die("not a .npy");
    size_t hl = npy[6] == 1 ? (size_t)(npy[8] | npy[9] << 8) + 10 : (size_t)(npy[8] | npy[9] << 8 | npy[10] << 16 | (size_t)npy[11] << 24) + 12;
    if (!strstr((char *)npy + 10, "'<i2'")) die("chrom.npy is not int16");
    CH = (int16_t *)(npy + hl);
    NNODE = (st.st_size - hl) / 2;
    CNT = calloc(NC + 2, sizeof *CNT);
    TALLY = calloc(NC, sizeof *TALLY);
    TOUCH = calloc(NC, sizeof *TOUCH);
    BGZF *fp = bgzf_open(argv[1], "r");
    if (!fp) die("open GAM");
    hfile_set_blksize(fp->fp, 16 << 20); /* 16 MB reads: the default 32 kB reads are latency-bound on BeeGFS */
    if (bgzf_seek(fp, start, SEEK_SET) < 0) die("seek");
    uint8_t *buf = NULL;
    size_t cap = 0;
    int64_t records = 0, groups = 0, secondary = 0;
    for (;;) {
        int64_t here = bgzf_tell(fp);
        if (end >= 0 && here == end) break;
        if (end >= 0 && here > end) die("walked past the segment end");
        uint64_t count, len;
        int r = bvarint(fp, &count);
        if (r == 0) break;
        if (r < 0) die("truncated group count");
        groups++;
        for (uint64_t i = 0; i < count; i++) {
            if (bvarint(fp, &len) != 1) die("truncated message length");
            if (len > cap) { cap = len * 2; buf = realloc(buf, cap); if (!buf) die("oom"); }
            if (bgzf_read(fp, buf, len) != (ssize_t)len) die("truncated message");
            if (i == 0) {
                if (len != 3 || memcmp(buf, "GAM", 3)) die("group tag is not GAM");
                continue;
            }
            records++;
            alignment(buf, buf + len, &secondary);
        }
    }
    if (end < 0 && bgzf_check_EOF(fp) != 1) fprintf(stderr, "gam_fullstats: warning: no BGZF EOF marker\n");
    bgzf_close(fp);
    for (int i = 0; i < NC + 2; i++) {
        if (!CNT[i][0]) continue;
        printf("%d", i - 2);
        for (int k = 0; k < NF; k++) printf("\t%lld", (long long)CNT[i][k]);
        printf("\n");
    }
    printf("#total\t%lld\t%lld\t%lld\n", (long long)records, (long long)groups, (long long)secondary);
    return 0;
}
