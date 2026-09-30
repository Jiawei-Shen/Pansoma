import numpy as np, glob, os, collections, sys
for base in sys.argv[1:]:
    tot = collections.Counter()
    for d in sorted(glob.glob(base+'/*')):
        c = collections.Counter()
        for f in sorted(glob.glob(d+'/*_labels.npy')):
            a = np.load(f, mmap_mode='r')
            u, n = np.unique(np.asarray(a), return_counts=True)
            for x, y in zip(u, n): c[int(x)] += int(y)
        nd = sum(np.load(f, mmap_mode='r').shape[0] for f in glob.glob(d+'/*_data.npy'))
        print(os.path.basename(base), os.path.basename(d), dict(c), 'data_n', nd, 'nshards', len(glob.glob(d+'/*_labels.npy')))
        tot.update(c)
    print('TOTAL', base, dict(tot))
