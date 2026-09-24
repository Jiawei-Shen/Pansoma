"""Optional native (C++) decoder: candidates.decode_alignment, about 30x faster.

candidates.decode_alignment stays the reference and the fallback. The compiled module
`_fastdecode` (built from fastdecode.cpp) is used only when
  * it imports (it was built for this Python and platform),
  * it was compiled from the fastdecode.cpp next to it (the source SHA-256 is compiled in), and
  * it reproduces the Python decoder exactly on a fixed set of synthetic records (checked once
    per process);
otherwise the builder decodes in Python and its manifest says why (`decoder`). A record the
native decoder raises on is decoded again in Python (any error is then the reference's), so a
task completes whenever it would complete in pure Python. The output is identical either way
(tests/test_native_decoder.py; README "native decoder").

Decoded columns stay in the C++ struct array (ColumnArray, ~24 bytes per column instead of a
~120-byte Column object); candidates.Column objects are made only where downstream code reads
them (the visit windows of candidate nodes), in cached 128-column blocks.

Build once per checkout and Python version, before `orchestrate prepare` (which freezes the
package, compiled module included):

    python -m indexed_gam_pipeline_v2.native compile [--cxx g++] [--no-portable]
    python -m indexed_gam_pipeline_v2.native check

Portability: standard C++17 plus the pybind11 headers, no -march / fast-math flags. On Linux,
libstdc++ and libgcc are linked statically when the toolchain has their static archives
(otherwise dynamically), so the module needs little beyond the C library; `compile` and `check`
report the glibc / libstdc++ symbol versions and CPU extensions a binary requires
(binary_requirements, also used for the graph index builder).
"""
import argparse
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys
import sysconfig
import tempfile

from indexed_gam_pipeline_v2 import vg_pb2
from indexed_gam_pipeline_v2.candidates import Candidate, Column, Observation, Read, Visit, decode_alignment, rc
from indexed_gam_pipeline_v2.common import sha256_file

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "fastdecode.cpp"
MODULE = "_fastdecode"
BUILD_INFO = HERE / f"{MODULE}.build.json"
CHOICES = ("auto", "native", "python")
BLOCK = 128
SELF_TEST_CASES = 300


# --- decoded reads ------------------------------------------------------------------

class ColumnArray:
    """Read.columns backed by the C++ struct array; candidates.Column objects on access.

    Supports len(), indexing, slicing, iteration and == like the list it replaces. Column
    objects are built per 128-column block and cached, so overlapping visit windows share them.
    """
    __slots__ = ("array", "n", "blocks")

    def __init__(self, array):
        self.array, self.n, self.blocks = array, len(array), {}

    def __len__(self):
        return self.n

    def _block(self, b):
        cols = self.blocks.get(b)
        if cols is None:
            a = self.array[b * BLOCK:(b + 1) * BLOCK]
            cols = self.blocks[b] = [
                Column(r, f, q, o, node, pos, rev, visit, bound) for r, f, q, o, node, pos, rev, visit, bound in zip(
                    a["read"].tobytes().decode("latin-1"), a["ref"].tobytes().decode("latin-1"), a["quality"].tolist(),
                    a["op"].tobytes().decode("latin-1"), a["node"].tolist(), a["pos"].tolist(),
                    a["reverse"].astype(bool).tolist(), a["visit"].tolist(), a["boundary"].astype(bool).tolist())]
        return cols

    def __getitem__(self, key):
        if isinstance(key, slice):
            start, stop, step = key.indices(self.n)
            if step != 1:
                return [self[i] for i in range(start, stop, step)]
            out = []
            for b in range(start // BLOCK, (stop - 1) // BLOCK + 1) if start < stop else ():
                block = self._block(b)
                out.extend(block[max(start - b * BLOCK, 0):min(stop - b * BLOCK, len(block))])
            return out
        i = key + self.n if key < 0 else key
        if not 0 <= i < self.n:
            raise IndexError("column index out of range")
        return self._block(i // BLOCK)[i % BLOCK]

    def __iter__(self):
        for b in range((self.n + BLOCK - 1) // BLOCK):
            yield from self._block(b)

    def __eq__(self, other):
        return list(self) == list(other)

    def __repr__(self):
        return f"ColumnArray({self.n} columns)"


def _unsupported_event(item):
    """The reference's unsupported-event dict (same keys in the same order) from the C++ tuple."""
    if item[0] == 0:
        _, node, start, ref, alt, kind, mapping, edit = item
        return dict(Candidate(node, start, ref, alt, kind).metadata(), mapping_index=mapping, edit_index=edit,
                    reason="indel_exceeds_limit")
    if item[0] == 1:
        _, node, start, ref, alt, _, mapping, edit = item
        return dict(node_id=node, start=start, ref=ref, alt=alt, mapping_index=mapping, edit_index=edit,
                    reason="complex_replacement_not_supported")
    _, node, start, kind, length, mapping, mappings = item
    return dict(node_id=node, start=start, event_type=kind, event_length=length, mapping_index=mapping,
                mappings=mappings, reason="indel_exceeds_limit_across_mappings")


def decode_native(module, alignment, sequences, max_indel=50, target_nodes=None, left_align=True):
    """decode_alignment's (Read, unsupported) computed by the native module."""
    serialized = alignment.SerializeToString(deterministic=True)  # the digest's bytes; C++ parses the same
    columns, visits, observations, moves, unsupported = module.decode(serialized, sequences, target_nodes,
                                                                      max_indel, left_align)
    visit_objects = [Visit(*v) for v in zip(visits["node"].tolist(), visits["start"].tolist(), visits["end"].tolist(),
                                            visits["reverse"].astype(bool).tolist(), visits["index"].tolist(),
                                            visits["first"].tolist(), visits["last"].tolist(),
                                            visits["node_length"].tolist())]
    obs = [Observation(Candidate(node, start, ref, alt, kind, path), visit, quality)
           for kind, node, start, ref, alt, path, visit, quality in observations]
    read = Read(alignment.name, hashlib.sha256(serialized).hexdigest(), alignment.mapping_quality,
                ColumnArray(columns), visit_objects, obs, moves)
    return read, [_unsupported_event(u) for u in unsupported]


class NativeDecoder:
    """decode_alignment through the native module; a record it raises on is decoded in Python."""

    def __init__(self, module, python=decode_alignment):
        self.module, self.python, self.fallbacks = module, python, 0

    def __call__(self, alignment, sequences, max_indel=50, target_nodes=None, left_align=True):
        try:
            return decode_native(self.module, alignment, sequences, max_indel, target_nodes, left_align)
        except Exception:  # noqa: BLE001 -- the Python decoder decides (and raises its own error, if any)
            pass
        result = self.python(alignment, sequences, max_indel, target_nodes=target_nodes, left_align=left_align)
        self.fallbacks += 1
        return result


# --- equivalence self-test ------------------------------------------------------------

def snapshot(result):
    """Everything decode_alignment returns, with value types, as comparable tuples."""
    read, unsupported = result
    return (tuple(tuple((type(getattr(c, f)).__name__, getattr(c, f)) for f in Column.__slots__) for c in read.columns),
            tuple((v.node, v.start, v.end, v.reverse, v.index, v.first, v.last, v.node_length) for v in read.visits),
            tuple((o.candidate, o.visit, o.quality, type(o.quality).__name__) for o in read.observations),
            tuple(read.moves), tuple(tuple(d.items()) for d in unsupported), (read.name, read.digest, read.mapq))


def _outcome(function, case):
    alignment, sequences, max_indel, targets, left = case
    try:
        return "ok", snapshot(function(alignment, sequences, max_indel, target_nodes=targets, left_align=left))
    except Exception as exc:  # noqa: BLE001 -- the exception type must match
        return "error", type(exc).__name__


def _random_edits(rng, length, offset):
    edits, cursor = [], offset
    for _ in range(rng.randint(1, 6)):
        room, kind = length - cursor, rng.choice("MMXXIDC")
        if kind == "I":
            inserted = "".join(rng.choice("ACGTN") for _ in range(rng.randint(1, 4)))
            edits.append((0, len(inserted), inserted))
        elif room == 0:
            continue
        elif kind == "M":
            f = rng.randint(1, room)
            edits.append((f, f, ""))
        elif kind == "X":
            f = rng.randint(1, min(room, 3))
            edits.append((f, f, "".join(rng.choice("ACGTN") for _ in range(f))))
        elif kind == "D":
            edits.append((rng.randint(1, min(room, 4)), 0, ""))
        else:
            f = rng.randint(1, min(room, 3))
            t = rng.choice([x for x in range(1, 5) if x != f])
            edits.append((f, t, "".join(rng.choice("ACGT") for _ in range(t))))
        cursor += edits[-1][0]
    return edits or [(0, 1, rng.choice("ACGT"))]


def _scattered(rng):
    """Mappings at random nodes and offsets, both orientations, every edit kind, N bases."""
    sequences = {n: "".join(rng.choice("ACGT" + ("N" if rng.random() < 0.2 else "")) for _ in range(rng.randint(1, 9)))
                 for n in range(1, rng.randint(2, 5))}
    specs = []
    for _ in range(rng.randint(1, 5)):
        node = rng.choice(sorted(sequences))
        offset = rng.randint(0, len(sequences[node]))
        specs.append((node, offset, rng.random() < 0.5, _random_edits(rng, len(sequences[node]), offset)))
    return sequences, specs


def _chain(rng):
    """A read along a chain of short low-complexity nodes, forward or reverse, with repeat indels
    (left-normalization across nodes, merging, the max_indel limit)."""
    motif = rng.choice(["A", "AT", "CAG", "AAC", "ACGT", "GC"])
    text = []
    while len(text) < rng.randint(8, 60):
        text.append(rng.choice("ACGTN") if rng.random() < 0.06 else motif[len(text) % len(motif)])
    text, sequences, i = "".join(text), {}, 0
    while i < len(text):
        size = rng.choice([1, 1, 2, 3, 4, 6])
        sequences[len(sequences) + 1], i = text[i:i + size], i + size
    reverse = rng.random() < 0.5
    order = sorted(sequences)[::-1] if reverse else sorted(sequences)
    first = rng.randint(0, len(order) - 1)
    last = rng.randint(first, len(order) - 1)
    specs = []
    for position, node in enumerate(order[first:last + 1]):
        ref = rc(sequences[node]) if reverse else sequences[node]
        offset = rng.randint(0, len(ref) - 1) if position == 0 and rng.random() < 0.3 else 0
        end = len(ref) if node != order[last] or rng.random() < 0.7 else rng.randint(offset, len(ref))
        edits, cursor = [], offset
        while cursor < end or (not edits and rng.random() < 0.5):
            room = end - cursor
            kind = rng.choice("MMMMXIIDDC") if room else "I"
            if kind == "I":
                base = ref[cursor - 1] if cursor else ref[0]
                length = rng.choice([1, 1, 2, 3, 4])
                edits.append((0, length, "".join(base if rng.random() < 0.8 else rng.choice("ACGTN") for _ in range(length))))
                if not room:
                    break
                continue
            if kind == "M":
                f = rng.randint(1, room)
                edits.append((f, f, ""))
            elif kind == "X":
                f = 1
                edits.append((1, 1, rng.choice([b for b in "ACGTN" if b != ref[cursor]])))
            elif kind == "D":
                f = rng.randint(1, min(room, 4))
                edits.append((f, 0, ""))
            else:
                f = rng.randint(1, min(room, 2))
                t = rng.choice([x for x in (1, 2, 3) if x != f])
                edits.append((f, t, "".join(rng.choice("ACGT") for _ in range(t))))
            cursor += f
        specs.append((node, offset, reverse, edits or [(0, 1, rng.choice("ACGT"))]))
    return sequences, specs


def synthetic_case(rng, index):
    """One random decode case: (alignment, sequences, max_indel, target_nodes, left_align)."""
    sequences, specs = (_chain if index % 2 else _scattered)(rng)
    a = vg_pb2.Alignment(name=f"synthetic{index}", mapping_quality=rng.choice([0, 20, 60, 255]))
    for node, offset, reverse, edits in specs:
        m = a.path.mapping.add()
        m.position.node_id, m.position.offset, m.position.is_reverse = node, offset, reverse
        ref = rc(sequences[node]) if reverse else sequences[node]
        cursor = offset
        for f, t, seq in edits:
            m.edit.add(from_length=f, to_length=t, sequence=seq)
            a.sequence += seq if seq else ref[cursor:cursor + f] if t else ""
            cursor += f
    choice = rng.random()
    a.quality = (b"" if choice < 0.15 else bytes(rng.randint(0, 255) for _ in a.sequence) if choice < 0.6
                 else bytes([30] * len(a.sequence)))
    if index % 97 == 0 and a.sequence:  # a malformed record: quality length != sequence length
        a.quality = a.quality + b"\x1e" if a.quality else b"\x1e\x1e"
    targets = rng.choice([None, set(sequences), set(rng.sample(sorted(sequences), max(1, len(sequences) // 2)))])
    return a, sequences, rng.choice([2, 3, 50]), targets, rng.random() < 0.9


def self_test(module, cases=SELF_TEST_CASES, seed=20260924):
    """None if `module` reproduces decode_alignment on `cases` synthetic records, else the first mismatch."""
    rng = random.Random(seed)
    for index in range(cases):
        case = synthetic_case(rng, index)
        expected = _outcome(decode_alignment, case)
        actual = _outcome(lambda *a, **kw: decode_native(module, *a, **kw), case)
        if expected != actual:
            return f"native decoder differs from the Python decoder on synthetic record {index}"
    return None


# --- loading and selection --------------------------------------------------------------

_LOADED = None


def load():
    """(module or None, info) -- imported, hash-checked and self-tested once per process."""
    global _LOADED
    if _LOADED is None:
        _LOADED = _load()
    return _LOADED


def _load():
    command = f"python -m {__package__}.native compile"
    try:
        module = importlib.import_module(f"{__package__}.{MODULE}")
    except ImportError as exc:
        return None, dict(reason=f"not built or not loadable here ({exc}); build it with `{command}`")
    expected = sha256_file(SOURCE) if SOURCE.exists() else None
    if getattr(module, "source_sha256", None) != expected:
        return None, dict(reason=f"stale build (compiled from another fastdecode.cpp); rebuild with `{command}`")
    try:
        problem = self_test(module)
    except Exception as exc:  # noqa: BLE001 -- any failure means: do not use it
        problem = f"self-test failed: {type(exc).__name__}: {exc}"
    if problem:
        return None, dict(reason=problem)
    info = dict(module=Path(module.__file__).name, source_sha256=expected, self_test_records=SELF_TEST_CASES)
    if BUILD_INFO.exists():
        build = json.loads(BUILD_INFO.read_text())
        info["build"] = {k: build.get(k) for k in ("compiler_version", "static_runtime", "python", "platform")}
    return module, info


def select_decoder(choice="auto", python=decode_alignment):
    """(decode function, manifest info) for --decoder auto|native|python.

    auto: native when available, else `python`; the PANSOMA_DECODER environment variable
    (native|python) overrides auto only. native: the native decoder or a ValueError.
    """
    requested = choice or "auto"
    if requested not in CHOICES:
        raise ValueError(f"--decoder must be one of {', '.join(CHOICES)}")
    effective = requested
    if requested == "auto" and os.environ.get("PANSOMA_DECODER") in ("native", "python"):
        effective = os.environ["PANSOMA_DECODER"]
    info = dict(requested=requested, **({"environment": f"PANSOMA_DECODER={effective}"} if effective != requested else {}))
    if effective == "python":
        return python, dict(info, used="python")
    module, details = load()
    if module is None:
        if effective == "native":
            raise ValueError(f"native decoder unavailable: {details['reason']}")
        return python, dict(info, used="python", **details)
    return NativeDecoder(module, python), dict(info, used="native", **details)


# --- building -------------------------------------------------------------------------

def _version(cxx):
    try:
        return subprocess.run([cxx, "--version"], capture_output=True, text=True).stdout.splitlines()[0]
    except (OSError, IndexError):
        return None


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


SELF_TEST_SNIPPET = """
import importlib.util, sys
from pathlib import Path
from indexed_gam_pipeline_v2 import native
path = next(Path(sys.argv[1]).glob("_fastdecode*"))
spec = importlib.util.spec_from_file_location("_fastdecode", path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
problem = native.self_test(module, cases=int(sys.argv[2]))
sys.exit(problem)
"""


def compile_decoder(cxx=None, portable=True, self_test_cases=5000):
    """Compile fastdecode.cpp into this package as _fastdecode<EXT_SUFFIX>, after a self-test.

    The module is compiled into a temporary directory, checked against the Python decoder on
    `self_test_cases` synthetic records in a fresh interpreter, and only then moved into place
    (with _fastdecode.build.json: compiler, flags, requirements). portable: link libstdc++ and
    libgcc statically when possible. Returns (path, build info).
    """
    import pybind11
    cxx = cxx or os.environ.get("CXX") or "g++"
    digest = sha256_file(SOURCE)
    output = HERE / (MODULE + (sysconfig.get_config_var("EXT_SUFFIX") or ".so"))
    flags = ["-O3", "-std=c++17", "-fPIC", "-shared", "-fvisibility=hidden", "-DNDEBUG",
             f'-DFASTDECODE_SOURCE_SHA256="{digest}"']
    includes = [f"-I{pybind11.get_include()}", f"-I{sysconfig.get_paths()['include']}"]
    if sys.platform == "darwin":
        flags += ["-undefined", "dynamic_lookup"]
    attempts = ([["-static-libstdc++", "-static-libgcc"]] if portable and sys.platform.startswith("linux") else []) + [[]]
    with tempfile.TemporaryDirectory(prefix=".fastdecode-build-", dir=HERE) as temp:
        target = Path(temp) / output.name
        errors = []
        for extra in attempts:
            result = subprocess.run([cxx, *flags, *includes, str(SOURCE), *extra, "-o", str(target)],
                                    capture_output=True, text=True)
            if result.returncode == 0:
                break
            errors.append(result.stderr[-3000:])
        else:
            raise RuntimeError("compiling the native decoder failed:\n" + "\n".join(errors))
        check = subprocess.run([sys.executable, "-c", SELF_TEST_SNIPPET, temp, str(self_test_cases)],
                               capture_output=True, text=True, cwd=HERE.parent)
        if check.returncode:
            raise RuntimeError(f"the compiled native decoder failed its self-test: {check.stdout}{check.stderr}")
        os.replace(target, output)
    info = dict(source_sha256=digest, compiler=cxx, compiler_version=_version(cxx), flags=flags + extra,
                static_runtime=bool(extra), python=sys.version.split()[0], platform=sysconfig.get_platform(),
                pybind11=pybind11.__version__, self_test_records=self_test_cases,
                requirements=binary_requirements(output))
    BUILD_INFO.write_text(json.dumps(info, indent=2) + "\n")
    return output, info


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    c = commands.add_parser("compile", help="compile fastdecode.cpp into this package (after a self-test)")
    c.add_argument("--cxx", help="C++17 compiler (default: $CXX, else g++)")
    c.add_argument("--no-portable", action="store_true", help="link libstdc++/libgcc dynamically")
    commands.add_parser("check", help="will the builder use the native decoder? (and why not)")
    args = parser.parse_args(argv)
    if args.command == "compile":
        output, info = compile_decoder(args.cxx, portable=not args.no_portable)
        print(json.dumps(dict(output=str(output), **info), indent=2))
        return
    module, info = load()
    report = dict(native_decoder="available" if module else "unavailable (the builder decodes in Python)", **info)
    if module is not None:
        report["requirements"] = binary_requirements(module.__file__)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
