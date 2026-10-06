"""Read-level estimate: what fraction of GAM alignments (reads) touch a branch node, and how many of those touch the
branch nodes without any edit? (s1 counts node visits exactly; this samples whole reads.)

Sample: 400 random group starts of the sorted GAM's GAI (seed 20261002), from each the next `PER_START` records; records
decoded with the pipeline's vg_pb2 (indexed_gam_pipeline_v4). Reads with MAPQ > 5 (discovery's rule) are counted.
Per read: visits (mappings with a node), branch visits (node not on the GRCh38 path, classify.py ref() rule), and edits
per visit as vg wrote them (an edit with from_length != to_length or a sequence = mismatch / indel / soft clip; no
left-normalization, unlike discovery). Classes:
  grch38_only          - every visit on GRCh38 nodes
  branch_all_perfect   - >= 1 branch visit, none of the branch visits has an edit
  branch_some_edited   - >= 1 branch visit with an edit
Also: reads with no edit at all (whole alignment perfect). Output <set>.read_sample.json in this folder.
Run: python s2_read_sample.py SET [PER_START]   (default PER_START 250 short reads / 25 long reads)
"""
import json, random, sys
from pathlib import Path
import numpy as np
import pysam
sys.path.insert(0, '/scratch/jshen/Github/Pansoma')
from indexed_gam_pipeline_v4.gam_reader import IndexedGam, group, decode   # noqa: E402
T = Path('/scratch/jshen/data/pansoma_v2_tensors')
GP = T / 'graph_index/hprc-v1.1-mc-grch38.d9.grch38_path'
OUT = Path(__file__).resolve().parent
s = sys.argv[1]
gam = json.load(open(T / s / 'discovery/discovery_report.json'))['gam']
per = int(sys.argv[2]) if len(sys.argv) > 2 else (250 if 'Illumina' in s else 25)
chrom, visits = np.load(GP / 'chrom.npy'), np.load(GP / 'visits.npy')
REF = (chrom >= 0) & (visits == 1)
reader = IndexedGam(gam, gam + '.gai', cache_bytes=1)
starts = sorted({st for _, _, runs in reader.bins for st, _ in runs})
pick = sorted(random.Random(20261002).sample(starts, min(400, len(starts))))
C = dict(records=0, mapq_gt5=0, grch38_only=0, branch_all_perfect=0, branch_some_edited=0, whole_read_perfect=0,
         visits=0, branch_visits=0, branch_visits_perfect=0, grch38_visits=0, grch38_visits_perfect=0)
with pysam.BGZFile(gam, 'rb') as stream:
    for st in pick:
        stream.seek(st)
        got = 0
        while got < per:
            msgs = group(stream)
            if msgs is None:
                break
            for raw in msgs:
                a = decode(raw)
                C['records'] += 1; got += 1
                if a.mapping_quality <= 5:
                    continue
                C['mapq_gt5'] += 1
                nb = nbe = 0; any_edit = False
                for m in a.path.mapping:
                    n = m.position.node_id
                    if not n:
                        continue
                    ed = any(e.from_length != e.to_length or e.sequence for e in m.edit)
                    any_edit |= ed
                    C['visits'] += 1
                    if REF[n]:
                        C['grch38_visits'] += 1; C['grch38_visits_perfect'] += not ed
                    else:
                        nb += 1; nbe += ed
                        C['branch_visits'] += 1; C['branch_visits_perfect'] += not ed
                C['whole_read_perfect'] += not any_edit
                C['grch38_only' if nb == 0 else 'branch_all_perfect' if nbe == 0 else 'branch_some_edited'] += 1
                if got >= per:
                    break
m = C['mapq_gt5']
C.update(set=s, gam=gam, starts_sampled=len(pick), per_start=per,
         frac_reads_touching_branch=(C['branch_all_perfect'] + C['branch_some_edited']) / m,
         frac_of_branch_reads_all_perfect=C['branch_all_perfect'] / max(1, C['branch_all_perfect'] + C['branch_some_edited']),
         frac_branch_visits_perfect=C['branch_visits_perfect'] / max(1, C['branch_visits']),
         frac_grch38_visits_perfect=C['grch38_visits_perfect'] / max(1, C['grch38_visits']),
         frac_whole_read_perfect=C['whole_read_perfect'] / m)
json.dump(C, open(OUT / f'{s}.read_sample.json', 'w'), indent=1)
print(json.dumps(C, indent=1))
