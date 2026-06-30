#!/usr/bin/env python3
"""
Convert a Python pickle dict to JSON.

Example:
    python3 pkl_to_json.py stats.pkl stats.json --pretty --sort-keys
"""

import argparse
import json
import pickle
import sys
from typing import Any


def load_pickle(path: str) -> Any:
    with open(path, "rb") as f:
        return pickle.load(f)


def save_json(obj: Any, path: str, indent: int | None, sort_keys: bool, ensure_ascii: bool) -> None:
    # Your original note: keys are already strings; keep as-is.
    with open(path, "w", encoding="utf-8") as out:
        json.dump(obj, out, indent=indent, sort_keys=sort_keys, ensure_ascii=ensure_ascii)


def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Convert a pickle file to JSON. Assumes the pickle contains JSON-serializable data."
    )
    p.add_argument("input", help="Input .pkl file")
    p.add_argument("output", help="Output .json file")
    p.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print JSON (indent=2).",
    )
    p.add_argument(
        "--sort-keys",
        action="store_true",
        help="Sort JSON object keys.",
    )
    p.add_argument(
        "--ensure-ascii",
        action="store_true",
        help="Escape non-ASCII characters (default: False).",
    )
    return p.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    try:
        obj = load_pickle(args.input)
        indent = 2 if args.pretty else None
        save_json(obj, args.output, indent=indent, sort_keys=args.sort_keys, ensure_ascii=args.ensure_ascii)
    except FileNotFoundError as e:
        print(f"File not found: {e}", file=sys.stderr)
        return 1
    except (pickle.UnpicklingError, json.JSONDecodeError) as e:
        print(f"Serialization error: {e}", file=sys.stderr)
        return 1
    except TypeError as e:
        # If something inside the pickle isn't JSON-serializable, tell the user which option might help.
        print(
            f"TypeError while converting to JSON: {e}\n"
            f"Tip: ensure the pickle contains only JSON-serializable types (dict/list/str/int/float/bool/None).",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))