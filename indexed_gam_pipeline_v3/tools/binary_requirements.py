"""objdump report of a binary: the newest GLIBC / GLIBCXX / CXXABI symbol versions it needs and the
CPU extensions its code uses (AVX / AVX-512 registers, BMI instructions).

    python -m indexed_gam_pipeline_v3.tools.binary_requirements indexed_gam_pipeline_v3/_fastdecode*.so [BINARY ...]

Prints {binary: report}. Used by tools.graph_index_build compile for the graph-index builder and on
demand for the native decoder module (`native compile` no longer reports it).
"""
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess


def binary_requirements(path):
    """Newest GLIBC / GLIBCXX / CXXABI symbol versions a binary needs and the CPU extensions
    its code uses (AVX/AVX-512 registers, BMI instructions), from objdump; {} without objdump."""
    objdump = shutil.which("objdump")
    if not objdump:
        return {}
    symbols = subprocess.run([objdump, "-T", str(path)], capture_output=True, text=True).stdout
    versions = {}
    for prefix in ("GLIBC", "GLIBCXX", "CXXABI"):
        found = set(re.findall(rf"\b{prefix}_([0-9][0-9.]*[0-9])", symbols))
        if found:
            versions[prefix] = max(found, key=lambda v: tuple(int(x) for x in v.split(".")))
    extensions = dict(avx_ymm=0, avx512_zmm=0, bmi=0)
    bmi = {"pdep", "pext", "bzhi", "shlx", "shrx", "sarx", "rorx", "mulx", "andn", "blsi", "blsr", "blsmsk"}
    with subprocess.Popen([objdump, "-d", "--no-show-raw-insn", str(path)], stdout=subprocess.PIPE, text=True,
                          stderr=subprocess.DEVNULL) as dump:
        for line in dump.stdout:
            if "%ymm" in line:
                extensions["avx_ymm"] += 1
            if "%zmm" in line:
                extensions["avx512_zmm"] += 1
            fields = line.rsplit("\t", 1)[-1].split()
            if fields and fields[0] in bmi:
                extensions["bmi"] += 1
    return dict(symbol_versions=versions, cpu_extensions=extensions)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("binaries", nargs="+", help="shared objects or executables")
    args = parser.parse_args(argv)
    missing = [b for b in args.binaries if not Path(b).is_file()]
    if missing:
        parser.error("not a file: " + ", ".join(missing))
    if not shutil.which("objdump"):
        parser.error("objdump not found on PATH")
    print(json.dumps({b: binary_requirements(b) for b in args.binaries}, indent=2))


if __name__ == "__main__":
    main()
