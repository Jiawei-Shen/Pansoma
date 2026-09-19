# Source this file from Bash before running Pansoma commands on this machine:
#   source scripts/use_vg.sh
# Select the user-specified vg instead of an older executable from Conda/PATH.
export PANSOMA_VG=/scratch/jshen/bin/vg_v1.77.0
if [[ ! -x "$PANSOMA_VG" ]]; then
    echo "Missing executable: $PANSOMA_VG" >&2
    return 1 2>/dev/null || exit 1
fi
export PATH="/scratch/jshen/bin:$PATH"
hash -r
# vg_v1.77.0 currently links to /scratch/jshen/bin/vg. Verify bare vg resolves
# to the requested executable, since existing scripts invoke it by that name.
if [[ ! "$(command -v vg)" -ef "$PANSOMA_VG" ]]; then
    echo "vg does not resolve to $PANSOMA_VG; check shell aliases/functions" >&2
    return 1 2>/dev/null || exit 1
fi
