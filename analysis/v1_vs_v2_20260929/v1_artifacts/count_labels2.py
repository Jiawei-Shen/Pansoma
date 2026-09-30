import numpy as np, glob, os, collections, sys, re
def summ(files):
    c=collections.Counter(); n=0; dt=set(); shp=set()
    for f in files:
        a=np.load(f, mmap_mode='r'); u,k=np.unique(np.asarray(a),return_counts=True); dt.add(str(a.dtype))
        for x,y in zip(u,k): c[int(x)]+=int(y)
    return c, dt
for base in sys.argv[1:]:
    tot=collections.Counter()
    percrom=collections.defaultdict(collections.Counter)
    for f in sorted(glob.glob(base+'/**/*_labels.npy', recursive=True)):
        m=re.search(r'(chr\w+?)_\d+_labels',os.path.basename(f)) or re.search(r'(chr\w+)',f)
        a=np.load(f, mmap_mode='r'); u,k=np.unique(np.asarray(a),return_counts=True)
        for x,y in zip(u,k): percrom[m.group(1)][int(x)]+=int(y); tot[int(x)]+=int(y)
    print('BASE',base, 'dtype', a.dtype, 'nfiles', len(glob.glob(base+'/**/*_labels.npy', recursive=True)))
    for ch in sorted(percrom, key=lambda s:int(s[3:]) if s[3:].isdigit() else 99): print('  ',ch, dict(percrom[ch]))
    print('  TOTAL', dict(tot), sum(tot.values()))
