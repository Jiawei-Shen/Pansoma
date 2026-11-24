#!/usr/bin/env python3
"""Sentence and caption clause of the letter for Response Figure 1.2, from results.json (collect.py).
Writes letter_facts.json only when all nine comparisons are finished. usage: python3 letter_facts.py"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
J = json.loads((HERE / "results.json").read_text())
rows = J["rows"]
done = [r for r in rows if r.get("own_node") and r["own_node"].get("TP") is not None]
assert len(rows) == 9, len(rows)
if len(done) < 9:
    raise SystemExit(f"{len(done)} of 9 comparisons finished; letter_facts.json not written")


def d(r, k):
    return r["full"][f"{k}_exact"] - r["own_node"][f"{k}_exact"]


dR, dP, dF = [d(r, "R") for r in rows], [d(r, "P") for r in rows], [d(r, "F1") for r in rows]
nR, nP, nF = sum(x > 0 for x in dR), sum(x > 0 for x in dP), sum(x > 0 for x in dF)
assert nR == 9, nR                    # the sentence below says "all nine"
rng = lambda v: f"{min(v):.3f}–{max(v):.3f}"
f_clause = (f"and F1 in all nine (by {rng(dF)})" if nF == 9 else
            f"and F1 in {nF} of nine")
sentence = ("We evaluated this change by masking this context at inference in the six-tumor models, which reduced SNV "
            f"recall in all nine comparisons (by {rng(dR)}) {f_clause} (Response Figure 1.2).")
WORD = ["none", "one", "two", "three", "four", "five", "six", "seven", "eight", "all nine"]
part = lambda n: WORD[n] if n == 9 else f"{WORD[n]} of nine"
caption_tail = (f"Recall is higher with the full tensor in all nine comparisons, precision in {part(nP)} and F1 in "
                f"{part(nF)}. Six-tumor models, seed 0; PASS at the threshold that recalls 95% of the validation somatic variants after population filtering.")
out = dict(sentence=sentence, caption_tail=caption_tail, n_recall_lower=nR, n_precision_lower=nP, n_f1_lower=nF,
           delta_recall=[round(x, 4) for x in dR], delta_precision=[round(x, 4) for x in dP], delta_f1=[round(x, 4) for x in dF],
           datasets=[f"{r['tumor']} {r['platform']}" for r in rows])
(HERE / "letter_facts.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
print(json.dumps(out, indent=1, ensure_ascii=False))
