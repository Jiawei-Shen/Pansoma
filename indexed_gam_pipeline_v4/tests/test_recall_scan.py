"""tensor_postprocessing.recall_scan: the parallel scan of the filtered streams gives truth_labels.scan_filtered's dict
(items and order) whatever the byte ranges, read blocks and workers, and recall_scan.label_run (truth_labels.label_run
verbatim, the scan swapped) writes truth_labels.label_run's labels, manifests, truth tables and recall reports."""
import contextlib
import hashlib
import inspect
import io
import json
import os
from pathlib import Path
import random
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np
import pysam

from . import golden
from .. import orchestrate
from ..tensor_postprocessing import recall_scan, truth_labels
from ..tensor_postprocessing.reference_path import ReferencePath

# SHA-256 of inspect.getsource(truth_labels.label_run). recall_scan.label_run is a verbatim copy: when label_run
# changes there, copy it again and pin the new hash here.
LABEL_RUN_SHA256 = "24e053b628be97419493331e60af6314b66b6eb69cc8dc177300c9fb2627a818"


def random_stream(rng, keys, others, final_newline):
    """Filtered-candidate ndjson in the shapes the scan meets, and some it never does: ids in and out of the truth
    keys, repeated with other reasons, reasons missing or empty, candidate_id not the first key or twice (the first
    counts) or empty (the next one counts), lines without an id, empty lines, CRLF, non-ASCII ids."""
    lines = []
    for _ in range(rng.randint(0, 60)):
        cid = rng.choice(keys + others)
        record = {"candidate_id": cid, "node_id": rng.randint(1, 10 ** 9), "ref": "", "alt": "ACGT"[rng.randint(0, 3)]}
        shape = rng.random()
        if shape < 0.1:
            record = {"site_id": "s", "candidate_id": cid}
        elif shape < 0.15:
            record = {"site_id": "s"}
        elif shape < 0.2:
            record = {"candidate_id": cid, "alleles": [{"candidate_id": rng.choice(keys)}]}
        elif shape < 0.25:
            record = {"candidate_id": "", "alleles": [{"candidate_id": cid}]}
        reasons = rng.choice([None, [], ["min_af"], ["min_variants", "min_af"], ["max_indel_len"], ["a b", "c"]])
        if reasons is not None and "candidate_id" in record:
            record["reasons"] = reasons
        line = json.dumps(record, ensure_ascii=rng.random() < 0.5)
        if rng.random() < 0.05:
            line = ""
        lines.append(line + ("\r" if rng.random() < 0.05 else ""))
    text = "\n".join(lines) + ("\n" if final_newline and lines else "")
    return text.encode()


def id_pools(rng):
    keys = [f"{rng.randint(1, 10 ** 6)}:{rng.randint(0, 99)}:SNP:A>C" for _ in range(12)] + ["5:1:INS:>ÉA", "7:0:DEL:A>"]
    others = [f"{rng.randint(1, 10 ** 6)}:0:INS:>GG" for _ in range(12)] + ["8:2:SNP:C>ü"]
    return keys, others


class LeakyKeys(set):
    """A key set whose iteration (the hash table) also yields ids that `in` rejects: hash hits that are not truth
    keys, as a 64-bit collision would give."""

    def __init__(self, keys, extra):
        super().__init__(keys)
        self.extra = list(extra)

    def __iter__(self):
        yield from set.__iter__(self)
        yield from self.extra

    def __len__(self):
        return set.__len__(self) + len(self.extra)


class ScanTest(unittest.TestCase):
    def test_ranges_start_at_line_starts_and_cover_the_file(self):
        rng = random.Random(1)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "filtered_candidates.ndjson"
            for trial in range(20):
                keys, others = id_pools(rng)
                data = random_stream(rng, keys, others, final_newline=trial % 2 == 0)
                path.write_bytes(data)
                for parts in range(1, len(data) + 3, max(1, len(data) // 97)):
                    cuts = recall_scan.byte_ranges(path, parts)
                    self.assertLessEqual(len(cuts), parts)
                    self.assertEqual([a for a, _ in cuts[1:]], [b for _, b in cuts[:-1]])
                    self.assertEqual((cuts[0][0], cuts[-1][1]) if cuts else (0, 0), (0, len(data)))
                    self.assertTrue(all(data[a - 1:a] == b"\n" for a, _ in cuts[1:]))

    def test_every_range_and_block_size_gives_the_lines_of_one_pass(self):
        """Ranges of 1-64 bytes, reads of 1, 2, 5, 64 bytes and BLOCK: the hits of the ranges in order are those of
        one range read at once (a hash table holding the keys and one other id)."""
        rng = random.Random(2)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "filtered_candidates.ndjson"
            for trial in range(6):
                keys, others = id_pools(rng)
                data = random_stream(rng, keys, others, final_newline=trial % 2 == 1)
                path.write_bytes(data)
                table = np.unique(np.array([recall_scan.key_hash(k.encode()) for k in keys + others[:1]], dtype=np.uint64))
                with mock.patch.object(recall_scan, "_TABLE", table):
                    whole = recall_scan.scan_range(path, 0, len(data), block=len(data) + 1)
                    for chunk in range(1, 65):
                        cuts = recall_scan.byte_ranges(path, -(-len(data) // chunk))
                        for block in (1, 2, 5, 64, recall_scan.BLOCK):
                            hits = [h for a, b in cuts for h in recall_scan.scan_range(path, a, b, block)]
                            self.assertEqual(hits, whole, (trial, chunk, block))
                self.assertTrue(data.count(b"\n") < 5 or whole)

    def test_parallel_scan_equals_the_serial_scan(self):
        """1, 3 and 8 workers over ranges of 3 and 7 bytes and the default ranges: the same dict, items and order
        (the last line of an id wins, the first gives its place); hash hits that are not keys are dropped."""
        rng = random.Random(3)
        with tempfile.TemporaryDirectory() as tmp:
            for trial in range(3):
                keys, others = id_pools(rng)
                data = random_stream(rng, keys, others, final_newline=trial % 2 == 0)
                (Path(tmp) / "filtered_candidates.ndjson").write_bytes(data)
                expected = list(truth_labels.scan_filtered(tmp, set(keys)).items())
                for workers, chunk in ((1, 3), (3, 7), (8, None)):
                    parts = -(-len(data) // chunk) if chunk else None
                    with self.subTest(trial=trial, workers=workers, chunk=chunk):
                        found = recall_scan.scan_filtered(tmp, set(keys), workers=workers, parts=parts, block=5)
                        self.assertEqual(list(found.items()), expected)
                leaky = recall_scan.scan_filtered(tmp, LeakyKeys(keys, others), workers=3)
                self.assertEqual(list(leaky.items()), expected)
            self.assertEqual(recall_scan.scan_filtered(Path(tmp) / "missing", set(keys)), {})
            (Path(tmp) / "filtered_candidates.ndjson").write_bytes(b"")
            self.assertEqual(recall_scan.scan_filtered(tmp, set(keys)), {})

    def test_invalid_utf8_ids_raise_as_in_the_serial_scan(self):
        with tempfile.TemporaryDirectory() as tmp:
            line = json.dumps({"candidate_id": "1:0:SNP:A>C", "reasons": ["min_af"]}).encode() + b"\n"
            (Path(tmp) / "filtered_candidates.ndjson").write_bytes(line * 50 + b'{"candidate_id": "9:\xff:SNP"}\n' + line)
            with self.assertRaises(UnicodeDecodeError) as serial:
                truth_labels.scan_filtered(tmp, {"1:0:SNP:A>C"})
            with self.assertRaises(UnicodeDecodeError) as parallel:
                recall_scan.scan_filtered(tmp, {"1:0:SNP:A>C"}, workers=3, parts=40)
            self.assertEqual(str(parallel.exception), str(serial.exception))


class LabelRunTest(unittest.TestCase):
    def test_label_run_is_truth_labels_label_run_with_the_scan_swapped(self):
        source = inspect.getsource(truth_labels.label_run)
        self.assertEqual(hashlib.sha256(source.encode()).hexdigest(), LABEL_RUN_SHA256,
                         "truth_labels.label_run changed: copy it into recall_scan.label_run, then pin its new SHA-256")
        self.assertEqual(inspect.getsource(recall_scan.label_run), source)

        def names(code):
            found = set(code.co_names)
            for constant in code.co_consts:
                if inspect.iscode(constant):
                    found |= names(constant)
            return found

        used = names(truth_labels.label_run.__code__) & set(vars(truth_labels))
        self.assertIn("scan_filtered", used)
        for name in sorted(used - {"scan_filtered"}):
            self.assertIs(vars(recall_scan).get(name), vars(truth_labels)[name], name)
        self.assertIs(recall_scan.label_run.__globals__["scan_filtered"], recall_scan.scan_filtered)

    def test_labels_and_recall_equal_truth_labels_on_the_o1_world(self):
        """O1's world (golden.mini_world), finalized by orchestrate (recall_scan.label_run), then labelled again by
        truth_labels.label_run and recall_scan.label_run on two copies, with truth alleles that have no tensor and
        whose keys the filtered streams name (repeated, other reasons, no reasons, an id in both kinds, no final
        newline): every file the same, labels.manifest.json apart from `created`, the same return value."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs, run = root / "inputs", root / "run"
            inputs.mkdir()
            argv = golden.mini_world(inputs, run)
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                orchestrate.main(["prepare", *argv])
            with mock.patch.dict(os.environ, dict(SLURM_CPUS_PER_TASK="2")), contextlib.redirect_stdout(io.StringIO()):
                orchestrate.run(run)
            config = json.loads((run / "config.json").read_text())
            post, labels, kinds = config["postprocess"], config["postprocess"]["labels"], list(config["variant_outputs"])
            self.assertTrue(json.loads((run / "status.json").read_text())["labeled"])
            fasta = pysam.FastaFile(labels["reference_fasta"])
            c1, c2 = fasta.fetch("chr1"), fasta.fetch("chr2")
            fasta.close()

            def alt(base):
                return {"A": "C", "C": "G", "G": "T", "T": "A"}[base]

            # no read carries a variant at chr1 offsets 6 and 31 or chr2 offset 2 (golden.mini_world)
            somatic = root / "somatic.vcf"
            somatic.write_text(Path(labels["somatic_vcf"]).read_text() + f"chr1\t7\t.\t{c1[6]}\t{alt(c1[6])}\t.\tPASS\t.\tGT\t0|1\n"
                               f"chr1\t31\t.\t{c1[30:32]}\t{c1[30]}\t.\tPASS\t.\tGT\t0|1\n")
            germline = root / "germline.vcf"
            germline.write_text(Path(labels["germline_vcf"]).read_text() + f"chr2\t3\t.\t{c2[2]}\t{alt(c2[2])}\t.\tPASS\t.\tGT\t0|1\n")
            locator = truth_labels.Locator(ReferencePath(post["reference_path"]))
            (snv, deletion), (germline_snv,) = ([a["keys"] for a in truth_labels.TruthSet(name, vcf, bed, labels["reference_fasta"],
                                                                                         locator).alleles[-n:]]
                                                for name, vcf, bed, n in (("somatic", somatic, labels["somatic_bed"], 2),
                                                                          ("germline", germline, labels["germline_bed"], 1)))
            self.assertTrue(snv and deletion and germline_snv)
            streams = {"SNV": [(snv[0], ["min_af"]), ("1:4:SNP:A>T", ["min_af"]), (germline_snv[0], None),
                               (snv[0], ["min_variants", "min_af"])],
                       "INDEL": [(deletion[0], []), (deletion[-1], ["max_indel_len"]), (snv[0], ["early_af"]),
                                 (germline_snv[0], ["min_variants"])]}
            copies = {name: root / name for name in ("serial", "parallel")}
            for copy in copies.values():
                shutil.copytree(config["tensors"], copy)
                for kind, items in streams.items():
                    lines = [json.dumps(dict(candidate_id=cid, node_id=1, **({} if r is None else dict(reasons=r))))
                             for cid, r in items]
                    with (copy / kind / "filtered_candidates.ndjson").open("a") as out:
                        out.write("\n".join(lines) + ("\n" if kind == "SNV" else ""))
            arguments = (post["reference_path"], labels["reference_fasta"], somatic, labels["somatic_bed"], germline,
                         labels["germline_bed"], root / "truth")
            serial = truth_labels.label_run(copies["serial"], kinds, *arguments)
            tables = {p.name: p.read_bytes() for p in (root / "truth").iterdir()}
            parallel = recall_scan.label_run(copies["parallel"], kinds, *arguments)
            self.assertEqual(parallel, serial)
            self.assertEqual({p.name: p.read_bytes() for p in (root / "truth").iterdir()}, tables)
            # the added alleles (the last of each table) are filtered with the reasons of their ids' last lines
            status = {name: [line.split("\t", 6)[-1] for line in (copies["serial"] / f"{name}.recall.tsv").read_text()
                             .splitlines()] for name in ("somatic", "germline")}
            self.assertEqual(status["somatic"][-2:], ["filtered\tearly_af", "filtered\tmax_indel_len"])
            self.assertEqual(status["germline"][-1:], ["filtered\tmin_variants"])
            files = {name: sorted(p.relative_to(copy) for p in copy.rglob("*") if p.is_file()) for name, copy in copies.items()}
            self.assertEqual(files["parallel"], files["serial"])
            for relative in files["serial"]:
                a, b = (copies[name] / relative for name in ("serial", "parallel"))
                if relative.name == "labels.manifest.json":
                    self.assertEqual(dict(json.loads(b.read_text()), created=None), dict(json.loads(a.read_text()), created=None))
                    self.assertEqual(list(json.loads(b.read_text())), list(json.loads(a.read_text())))
                else:
                    self.assertEqual(b.read_bytes(), a.read_bytes(), str(relative))


if __name__ == "__main__":
    unittest.main()
