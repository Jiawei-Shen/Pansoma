#!/usr/bin/env python3
import argparse, gzip, sqlite3, pickle, time, gc
from collections import defaultdict
import vg_pb2

def read_varint(stream):
    v, shift = 0, 0
    while True:
        b = stream.read(1)
        if not b:
            raise EOFError
        v |= (b[0] & 0x7F) << shift
        if not (b[0] & 0x80):
            break
        shift += 7
    return v

def file_is_gzip(path):
    with open(path,"rb") as f:
        return f.read(2)==b"\x1f\x8b"

def gam_record_iter(path, tag="GAM"):
    openf = gzip.open if file_is_gzip(path) else open
    with openf(path,"rb") as f:
        while True:
            try:
                n = read_varint(f)
                if   n==0: continue
                tlen = read_varint(f)
                t    = f.read(tlen).decode()
                if t!=tag:
                    for _ in range(n-1):
                        skip = read_varint(f); f.seek(skip,1)
                    continue
                for _ in range(n-1):
                    sz = read_varint(f)
                    yield f.read(sz)
            except EOFError:
                break

def process_alignment(msg, wanted, chrom):
    a = vg_pb2.Alignment(); a.ParseFromString(msg)
    if chrom and not any(p.name==chrom for p in a.refpos):
        return []
    seq, qual, mq = a.sequence, a.quality, a.mapping_quality
    off=0; rows=[]
    for m in a.path.mapping:
        nid = m.position.node_id
        if nid not in wanted:
            for e in m.edit:
                off += max(e.from_length, len(e.sequence))
            continue
        node_off = m.position.offset
        strand   = "-" if m.position.is_reverse else "+"
        parts, qb = [], []
        for e in m.edit:
            if e.from_length:
                frag = seq[off:off+e.from_length]
                parts.append(frag.lower() if e.sequence else frag)
                qb.extend(qual[off:off+e.from_length])
                off += e.from_length
            elif e.sequence:
                L = len(e.sequence)
                parts.append(seq[off:off+L].lower())
                qb.extend(qual[off:off+L])
                off += L
        rows.append((nid, node_off, strand, "".join(parts), bytes(qb).hex(), mq))
    return rows

def flush(conn, buf):
    if not buf: return
    conn.execute("BEGIN")
    conn.executemany(
        "INSERT INTO segments VALUES (?,?,?,?,?,?)", buf
    )
    conn.execute("COMMIT")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("gam"); p.add_argument("stats")
    p.add_argument("sqlite"); p.add_argument("--milestone",type=int,default=1_000_000)
    p.add_argument("--chr",default="")
    args = p.parse_args()

    # load and filter nodes
    print(f"\nload and filter nodes")
    with open(args.stats,"rb") as f: stats=pickle.load(f)
    wanted = {int(n) for n,s in stats.items() if s["not_perfect"]>1 and s["not_perfect"]/(s["perfect"]+s["not_perfect"])>0.10}
    tot, unp = len(stats), sum(1 for s in stats.values() if s["not_perfect"]>0)
    print(f"\nTotal nodes: {tot}")
    print(f"Nodes with ≥1 un‑perfect: {unp} ({unp/tot*100:.2f}%)")
    print(f"Nodes passing filter: {len(wanted)} ({len(wanted)/tot*100:.2f}%)\n")
    del stats; gc.collect()

    # init SQLite
    conn = sqlite3.connect(args.sqlite)
    # 1) speed up journaling
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=OFF;")
    conn.execute("PRAGMA cache_size = -1000000;")
    conn.execute("PRAGMA temp_store=MEMORY;")
    # 2) wrap all inserts in one big transaction
    conn.execute("BEGIN")
    conn.execute("DROP TABLE IF EXISTS segments")
    conn.execute("""CREATE TABLE segments (
        node_id INTEGER,
        offset  INTEGER,
        strand  TEXT,
        seq     TEXT,
        bq      TEXT,
        rq      INTEGER
    )""")
    # conn.execute("CREATE INDEX idx_node ON segments(node_id)")

    buf = []
    reads = 0
    next_m = args.milestone
    t0 = time.perf_counter()


    for raw in gam_record_iter(args.gam):
        rows = process_alignment(raw, wanted, args.chr)
        if rows:
            buf.extend(rows)
        reads += 1

        if reads >= next_m:
            # flush(conn, buf)
            if not buf: break
            # conn.execute("BEGIN")
            conn.executemany(
                "INSERT INTO segments VALUES (?,?,?,?,?,?)", buf
            )
            # conn.execute("COMMIT")
            buf.clear()
            dt = time.perf_counter() - t0
            print(f"{reads} reads | {dt:.1f}s")
            next_m += args.milestone

    flush(conn, buf)

    conn.execute("COMMIT")
    conn.execute("CREATE INDEX idx_node ON segments(node_id)")
    conn.close()
    print(f"{reads} reads | {dt:.1f}s")
    print(f"\n✅ Done. Total reads: {reads}")

if __name__=="__main__":
    main()
