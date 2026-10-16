"""E11a cross-check: whole-GAM `vg stats -a` reports that already exist next to some of the GAMs used here (written in
2026-03..05 from the unsorted GAMs that were later sorted, same records). Writes results/vg_stats_crosscheck.md.
vg stats definitions: aligned = has a path; perfect = every edit a match (no soft clip); gapless = no indel
(soft clips allowed); softclip events = reads with a soft clip at an end, counted per end."""
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
FILES = [
    ('HG008T_PacBio', 'd9', '/scratch/jshen/data/HG008_GIAB/raw_sequencing_data/Liss_lab_PacBio_Revio_20240125/HG008_PacBio_gam_stats.txt'),
    ('HG008T_ONT', 'd9', '/scratch/jshen/data/HG008_GIAB/raw_sequencing_data/Liss_lab_Northeastern-ONT-UL-20241216/HG008_ONT_gam_stats.txt'),
    ('HG008T_Illumina', 'd9', '/scratch/jshen/data/HG008_GIAB/vg_giraffe_GAM/Liss_lab_BCM_Illumina-WGS_20240313/HG008-T_p23_BCM_IlluminaWGS_20240313.vg.stats'),
    ('HG008T_Illumina', 'linear', '/scratch/jshen/data/HG008_GIAB/raw_sequencing_data/Liss_lab_BCM_Illumina-WGS_20240313/vg_linear_hg38/HG008T_BCM_IlluminaWGS_hg38_linear_vg.gam.stats.txt'),
    ('COLO829T_fiberseq', 'linear', '/scratch/jshen/data/COLO829T/COLO829T_fiberseq/vg_linear_hg38/SMHTCOLO829T_fiberseq_hg38_linear_vg.gam.stats.txt'),
    ('COLO829T_ONT', 'd9', '/scratch/jshen/data/COLO829T/COLO829T_ONT/COLO829T_ONT_merged.vg_stats.txt'),
    ('COLO829T_ONT', 'linear', '/scratch/jshen/data/COLO829T/COLO829T_ONT/vg_linear_hg38/COLO829T_ONT_hg38_linear_vg.gam.stats.txt'),
]


def parse(path):
    t = Path(path).read_text()
    g = lambda k: int(re.search(rf'^{k}: (\d+)', t, re.M).group(1))
    sc = re.search(r'^Softclips: (\d+) bp in (\d+) read events', t, re.M)
    mq = re.search(r'^Mapping quality: mean ([\d.]+)', t, re.M)
    return dict(total=g('Total alignments'), aligned=g('Total aligned'), perfect=g('Total perfect'),
                gapless=g('Total gapless \\(softclips allowed\\)'), softclip_events=int(sc.group(2)) if sc else None,
                mapq_mean=float(mq.group(1)) if mq else None)


L = ['# Existing whole-GAM `vg stats -a` reports (cross-check of the sampled numbers)', '',
     '| set | arm | records | unmapped % | perfect % of all | gapless % of all | soft-clip events per 100 mapped reads | mean MAPQ | file |',
     '|---|---|---:|---:|---:|---:|---:|---:|---|']
for s, a, p in FILES:
    d = parse(p)
    L.append(f"| {s} | {a} | {d['total']:,} | {100 * (d['total'] - d['aligned']) / d['total']:.2f} | "
             f"{100 * d['perfect'] / d['total']:.2f} | {100 * d['gapless'] / d['total']:.2f} | "
             f"{100 * d['softclip_events'] / d['aligned']:.2f} | {d['mapq_mean']:.2f} | `{p}` |")
(HERE / 'results' / 'vg_stats_crosscheck.md').write_text('\n'.join(L) + '\n')
print('\n'.join(L))
