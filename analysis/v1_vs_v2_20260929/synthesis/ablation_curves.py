import contextlib, io, sys, pickle
sys.path.insert(0, "/scratch/jshen/data/pansoma_net_v2_runs/v1_vs_v2_20260929/synthesis")
with contextlib.redirect_stdout(io.StringIO()):
    import my_check as M
for r in ("HG008_Illumina_SNV_keephard_w100_nopartial",):
    M.V2R[r] = M.read_v2(r)
PONF = {k: any(v) for k, v in pickle.load(open(M.CACHE, "rb")).items()}
for r in ("HG008_Illumina_SNV_base", "HG008_Illumina_SNV_keephard_w100", "HG008_Illumina_SNV_allnon_w100", "HG008_Illumina_SNV_keephard_w100_nopartial"):
    s = M.V2R[r][0]
    for tag, c in (("noPoN", M.curve(s)), ("PoN", M.curve(s, pon=PONF))):
        m = M.summ(c); print(f"{r.replace('HG008_Illumina_SNV_',''):<22} {tag:<5}", {k: m[k] for k in ("R@P0.10", "R@P0.15", "R@P0.20", "P@R0.5", "P@R0.7", "P@R0.8", "P@R0.9", "maxF1")})
for tag, c in (("noPoN", M.curve(M.V1["v1_e053"])), ("ownPoN", M.curve(M.V1["v1_e053_ownPoN"]))):
    m = M.summ(c); print(f"{'v1 e053':<22} {tag:<5}", {k: m[k] for k in ("R@P0.10", "R@P0.15", "R@P0.20", "P@R0.5", "P@R0.7", "P@R0.8", "maxF1")})
