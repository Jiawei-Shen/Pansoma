"""Build, from overview_<DATE>.md and points/<R-ID>.md in reviewer order:
  response_point_by_point_<DATE>.md  overview + every point file in full (状态 / Comment / Response / 证据 / 还缺什么)
  response_letter_<DATE>.md          the letter draft: per point the heading, the reviewer's asks and the English Response
usage: python assemble_point_by_point.py [DATE]   (default 20261015). Standard library only; rerun after editing points."""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ORDER = {
    "Reviewer 1": ["R1.1", "R1.2", "R1.3", "R1.4", "R1.5", "R1.6", "R1.7", "R1.8", "R1.m1", "R1.m2"],
    "Reviewer 2": ["R2.1", "R2.2", "R2.3", "R2.4", "R2.5", "R2.6", "R2.7", "R2.m"],
    "Reviewer 3": ["R3.1", "R3.2", "R3.3", "R3.4", "R3.5", "R3.6", "R3.7", "R3.8", "R3.9", "R3.10",
                   "R3.m1", "R3.m2", "R3.m3", "R3.m4"],
}


def section(text, start, stops):
    """The text from the line starting with `start` up to the next line starting with any of `stops`."""
    i = text.find(start)
    if i < 0:
        return ""
    j = min([k for k in (text.find(s, i + len(start)) for s in stops) if k >= 0] or [len(text)])
    return text[i:j].strip()


def main():
    date = sys.argv[1] if len(sys.argv) > 1 else "20261015"
    full = [(HERE / f"overview_{date}.md").read_text().rstrip() + "\n"]
    letter = [f"# Pansoma: response to the reviewers (draft {date})\n\nThe English responses of "
              f"`response_point_by_point_{date}.md`, without the internal evidence. `[待补: …]` marks numbers still "
              "pending.\n"]
    words = []
    for n, (reviewer, points) in enumerate(ORDER.items(), 1):
        full.append(f"\n---\n\n## 6.{n} {reviewer}\n")
        letter.append(f"\n## {reviewer}\n")
        for p in points:
            text = (HERE / "points" / f"{p}.md").read_text().strip()
            full.append("\n" + text + "\n")
            head = text.splitlines()[0]
            comment = section(text, "**Comment", ["**Response", "**证据"])
            resp = section(text, "**Response.**", ["**证据", "**还缺什么"])
            body = resp.replace("**Response.**", "**Response.**", 1)
            letter.append(f"\n{head}\n\n{comment}\n\n{body}\n")
            prose = re.sub(r"^\|.*$", "", resp, flags=re.M)  # word count without tables
            words.append((p, len(prose.replace("**Response.**", "").split())))
    (HERE / f"response_point_by_point_{date}.md").write_text("".join(full))
    (HERE / f"response_letter_{date}.md").write_text("".join(letter))
    print(f"{sum(len(v) for v in ORDER.values())} points -> response_point_by_point_{date}.md, response_letter_{date}.md")
    print("Response words (tables excluded):", ", ".join(f"{p} {w}" for p, w in words))


if __name__ == "__main__":
    main()
