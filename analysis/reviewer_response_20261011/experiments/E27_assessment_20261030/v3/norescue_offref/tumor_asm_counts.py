import csv, glob, json, os, sys
sys.path.insert(0, '/scratch/jshen/Github/Pansoma/tmp/assessment_20261030/e9c_norescue/code')
import collect as C, e9c_lib as L
asm = {s: C.load_asm(s) for s in C.CATS}
print('run | unplaced PASS | ALT in tumor asm | in tumor asm, not in normal asm | exclusive tumor-only (not population)')
for p in sorted(glob.glob(f'{L.DATA}/runs/*/calls.tsv')):
    d = os.path.dirname(p); rj = json.load(open(f'{d}/run.json')); s = rj['sample']
    hits, mapped = asm[s]
    n = in_t = t_only = excl = 0
    for r in csv.DictReader(open(p), delimiter='\t'):
        if r['group'] != 'unplaced_pass': continue
        n += 1
        if r.get('tensor') != 'True': continue
        alt = 'alt_read' if len(r.get('alt_read', '')) >= 40 else 'alt_graph'
        if not r.get(alt) or not r.get(f'{alt}_site'): continue
        qa = C.qid(r[alt], r[f'{alt}_site'])
        ha = hits.get(qa, {})
        if ha.get('tumor', 0) > 0: in_t += 1
        c, _ = C.asm_class(r, s, asm[s])
        if c == 'tumor_only':
            t_only += 1
            if not C.pop_tag(r): excl += 1
    print(os.path.basename(d).replace('rcpu__', '').replace('__norescue__vcf_chr1', ''), n, in_t, t_only, excl, sep=' | ')
