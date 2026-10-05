"""gam_reader.BgzfReader against pysam.BGZFile (htslib): seek/tell/read on BGZF files with tiny, empty and
multi-group blocks; read_group against group() on random, truncated and malformed group streams; position();
corrupt blocks raise OSError like pysam."""
import io
from pathlib import Path
import random
import struct
import tempfile
import unittest
import zlib

import pysam

from .fixtures import encode_varint
from ..gam_reader import BgzfReader, group

EOF_BLOCK = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def bgzf_block(data, level=6):
    """One BGZF block holding `data` (at most 0xff00 bytes; b"" gives the EOF marker block)."""
    packer = zlib.compressobj(level, zlib.DEFLATED, -15)
    deflated = packer.compress(data) + packer.flush()
    size = 18 + len(deflated) + 8
    return (b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", size - 1)
            + deflated + struct.pack("<II", zlib.crc32(data), len(data)))


def bgzf_bytes(chunks, eof=True):
    """A BGZF file of one block per chunk (b"" = an empty block), with the EOF marker block unless not `eof`;
    returns (bytes, [(block address, data)])."""
    out, blocks = bytearray(), []
    for chunk in chunks:
        blocks.append((len(out), chunk))
        out += bgzf_block(chunk)
    if eof:
        blocks.append((len(out), b""))
        out += EOF_BLOCK
    return bytes(out), blocks


def random_cuts(rng, payload, boundaries=(), empty=0.1, largest=0xff00):
    """`payload` cut into block chunks of random sizes (some cuts on `boundaries`), with random empty blocks."""
    cuts, at = [], 0
    marks = sorted(set(boundaries))
    while at < len(payload):
        if rng.random() < empty:
            cuts.append(b"")
        kind = rng.random()
        if kind < 0.3:
            size = rng.randint(1, 12)
        elif kind < 0.5 and marks:
            later = [m for m in marks if m > at]
            size = (rng.choice(later[:4]) - at) if later else rng.randint(1, 4000)
        elif kind < 0.9:
            size = rng.randint(1, 4000)
        else:
            size = rng.randint(1, largest)
        cuts.append(payload[at:at + min(size, largest)])
        at += min(size, largest)
    if rng.random() < empty:
        cuts.append(b"")
    return cuts


def group_bytes(messages, tag=b"GAM"):
    return encode_varint(len(messages) + 1) + encode_varint(len(tag)) + tag + b"".join(
        encode_varint(len(m)) + m for m in messages)


def random_groups(rng, count):
    """A group stream: GAM groups of 0-6 messages (some empty), foreign-tagged and count-0 groups; and the
    payload offsets where groups and messages start."""
    payload, marks = bytearray(), []
    for _ in range(count):
        marks.append(len(payload))
        kind = rng.random()
        if kind < 0.05:
            payload += encode_varint(0)
            continue
        tag = b"GAM" if kind < 0.85 else rng.choice([b"PARAMS_JSON", b"X", b""])
        messages = [bytes(rng.getrandbits(8) for _ in range(rng.choice([0, 1, 5, 130, rng.randint(0, 3000)])))
                    for _ in range(rng.randint(0, 6))]
        payload += group_bytes(messages, tag)
    return bytes(payload), marks


def outcomes(read, stream, limit=10000):
    """[group values or ('error', type, message)] read until EOF or the first error."""
    found = []
    for _ in range(limit):
        try:
            value = read(stream)
        except (ValueError, OSError, MemoryError) as error:
            found.append(("error", type(error).__name__, str(error)))
            return found
        found.append(value)
        if value is None:
            return found
    return found


class BgzfReaderTest(unittest.TestCase):
    def write(self, directory, name, data):
        path = Path(directory) / name
        path.write_bytes(data)
        return path

    def test_seek_tell_read_match_pysam(self):
        rng = random.Random(11)
        with tempfile.TemporaryDirectory() as directory:
            for case in range(60):
                payload = bytes(rng.getrandbits(8) for _ in range(rng.randint(0, 30000)))
                data, blocks = bgzf_bytes(random_cuts(rng, payload, empty=0.2), eof=rng.random() < 0.9)
                path = self.write(directory, f"r{case}.gz", data)
                offsets = [(address << 16) | at for address, chunk in blocks for at in range(0, len(chunk) + 1)
                           if at in (0, len(chunk)) or rng.random() < 0.05] + [len(data) << 16]
                with pysam.BGZFile(str(path), "rb") as reference, BgzfReader(path) as reader:
                    for step in range(80):
                        if rng.random() < 0.3:
                            offset = rng.choice(offsets)
                            reference.seek(offset)
                            reader.seek(offset)
                            if rng.random() < 0.5:
                                reader.limit(rng.choice([None, offset, offset + (rng.randint(0, 5) << 16)]))
                        size = rng.choice([0, 1, 2, 7, 100, 5000, 70000])
                        expected = reference.read(size)
                        self.assertEqual(reader.read(size), expected, (case, step))
                        self.assertEqual(reader.tell(), reference.tell(), (case, step))
                # Offsets past the end of their block fail in both; after an empty block the offset applies to the
                # next data block (htslib keeps a seek's offset), past EOF it fails.
                odd = [(address << 16) | (len(chunk) + rng.randint(1, 3)) for address, chunk in blocks]
                for offset in rng.sample(odd, min(4, len(odd))):
                    found = []
                    for stream in (pysam.BGZFile(str(path), "rb"), BgzfReader(path)):
                        stream.seek(offset)
                        try:
                            result = (stream.read(3), stream.tell())
                        except OSError as error:
                            result = str(error)
                        try:
                            stream.close()
                            found.append((result, None))
                        except OSError as error:
                            found.append((result, str(error)))
                    self.assertEqual(found[1], found[0], (case, offset))

    def test_position_names_where_a_read_continues(self):
        rng = random.Random(12)
        with tempfile.TemporaryDirectory() as directory:
            for case in range(40):
                payload = bytes(rng.getrandbits(8) for _ in range(rng.randint(1, 20000)))
                data, blocks = bgzf_bytes(random_cuts(rng, payload, empty=0.3))
                path = self.write(directory, f"p{case}.gz", data)
                flat, before = {}, 0
                for address, chunk in blocks:
                    flat[address] = before
                    before += len(chunk)
                with pysam.BGZFile(str(path), "rb") as reference, BgzfReader(path) as reader:
                    for address, chunk in blocks:
                        for at in {0, len(chunk), rng.randint(0, len(chunk))}:
                            offset = (address << 16) | at
                            reference.seek(offset)
                            rest = reference.read(1 << 20)
                            found = reader.position(offset)
                            if found is None:
                                self.assertEqual(rest, b"")
                            else:
                                self.assertEqual(payload[flat[found[0]] + found[1]:], rest)
                                self.assertLess(found[1], len(dict(blocks)[found[0]]))
                            reader.seek(offset)
                            self.assertEqual(reader.position(), found)  # unloaded: computed from tell()
                            if found is not None:
                                reader.read(1)
                                reader.seek(offset)
                                self.assertEqual(reader.read(len(rest)), rest)
                                self.assertEqual(reader.position(), None)

    def test_read_group_matches_group_on_random_and_truncated_streams(self):
        rng = random.Random(13)
        with tempfile.TemporaryDirectory() as directory:
            for case in range(150):
                payload, marks = random_groups(rng, rng.randint(0, 25))
                kind = rng.random()
                if kind < 0.3 and payload:
                    payload = payload[:rng.randint(0, len(payload) - 1)]  # truncated anywhere
                elif kind < 0.4:
                    payload += rng.choice([b"\x80", b"\xff" * 9 + b"\x02", b"\xff" * 10 + b"\x01", b"\x81\x80"])
                elif kind < 0.5:  # a length beyond EOF; one no allocation can hold (pysam: MemoryError)
                    payload += encode_varint(3) + encode_varint(3) + b"GAM" + encode_varint(rng.choice([999, 1 << 62]))
                data, _ = bgzf_bytes(random_cuts(rng, payload, boundaries=marks, empty=0.15))
                path = self.write(directory, f"g{case}.gam", data)
                with pysam.BGZFile(str(path), "rb") as reference:
                    expected = outcomes(group, reference)
                if not payload.endswith(encode_varint(1 << 62)):
                    self.assertEqual(outcomes(group, io.BytesIO(payload)), expected, case)
                with BgzfReader(path) as reader:
                    self.assertEqual(outcomes(BgzfReader.read_group, reader), expected, case)
                with pysam.BGZFile(str(path), "rb") as reference, BgzfReader(path) as reader:
                    while True:  # the same values and tell() after every group as group() over pysam
                        try:
                            value = group(reference)
                        except (ValueError, OSError, MemoryError):
                            break
                        positions = []
                        self.assertEqual(reader.read_group(positions), value, case)
                        self.assertEqual(reader.tell(), reference.tell(), case)
                        self.assertEqual(len(positions), len(value or []))
                        if value is None:
                            break
                        resume = reference.tell()
                        for (start, end), message in zip(positions, value):
                            reference.seek(start)
                            self.assertEqual(reference.read(len(encode_varint(len(message))) + len(message)),
                                             encode_varint(len(message)) + message)
                            self.assertEqual(reference.tell(), end)
                        reference.seek(resume)

    def test_corrupt_blocks_raise_like_pysam(self):
        """Each corruption htslib refuses raises OSError('Error reading from BGZFile'), and leaving the `with`
        block OSError('Error closing BGZFile object'), as with pysam; a wrong ISIZE (htslib does not check it)
        reads on, except in a strict reader."""
        payload, marks = random_groups(random.Random(14), 30)
        chunks = random_cuts(random.Random(15), payload, boundaries=marks, empty=0.0)
        data, blocks = bgzf_bytes(chunks)
        middle = blocks[len(blocks) // 2][0]
        size = struct.unpack_from("<H", data, middle + 16)[0] + 1

        def flip(at):
            return data[:at] + bytes([data[at] ^ 1]) + data[at + 1:]

        corrupt = {
            "crc": flip(middle + size - 8), "deflate": data[:middle + 18] + b"\xff\xff\xff" + data[middle + 21:],
            "magic": flip(middle), "size": data[:middle + 16] + b"\x05\x00" + data[middle + 18:],
            "truncated block": data[:middle + size - 3], "truncated header": data[:middle + 5],
            "over 64 KiB": data[:middle] + bgzf_block(b"x" * 70000) + data[middle:],
        }

        def attempt(make, action):
            """(the action's result or error, close()'s error)."""
            stream = make()
            try:
                result = action(stream)
            except (ValueError, OSError) as error:
                result = ("error", type(error).__name__, str(error))
            try:
                stream.close()
                return result, None
            except OSError as error:
                return result, str(error)

        failed = ("error", "OSError", "Error reading from BGZFile")
        with tempfile.TemporaryDirectory() as directory:
            for name, broken in corrupt.items():
                path = self.write(directory, name.replace(" ", "_") + ".gam", broken)
                for action in (lambda s: s.read(len(payload) + 1), lambda s: outcomes(group, s)):
                    expected = attempt(lambda: pysam.BGZFile(str(path), "rb"), action)
                    self.assertEqual(expected[1], "Error closing BGZFile object", name)
                    self.assertTrue(expected[0] == failed or expected[0][-1] == failed, name)
                    self.assertEqual(attempt(lambda: BgzfReader(path), action), expected, name)
                self.assertEqual(attempt(lambda: BgzfReader(path), lambda s: outcomes(BgzfReader.read_group, s)),
                                 attempt(lambda: pysam.BGZFile(str(path), "rb"), lambda s: outcomes(group, s)), name)
            path = self.write(directory, "isize.gam", flip(middle + size - 4))
            with pysam.BGZFile(str(path), "rb") as reference, BgzfReader(path) as reader:
                self.assertEqual(reader.read(len(payload) + 1), reference.read(len(payload) + 1))
            self.assertEqual(attempt(lambda: BgzfReader(path, strict=True), lambda s: s.read(len(payload) + 1)),
                             (failed, "Error closing BGZFile object"))


if __name__ == "__main__":
    unittest.main()
