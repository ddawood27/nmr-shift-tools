#!/usr/bin/env python3
"""
ss_snippet_fullclass.py

Classifies full AA-snippet secondary structure.

For each exact occurrence of a sequence snippet (2–10 aa):
- If ALL residues share the same SS → classify as that SS (H, E, .)
- Otherwise → 'mixed ss'

Defaults are hard-coded for merged_shifts_dssp.csv:
  entry_id, CHAINID, label_seq_id, RESNAME, secondary_structure
"""

from __future__ import annotations

import argparse
from collections import Counter
import sys
import pandas as pd


# ---- FIXED COLUMN NAMES (do not change unless file format changes) ----
COL_ID = "entry_id"
COL_CHAIN = "CHAINID"
COL_RESNUM = "label_seq_id"
COL_AA = "RESNAME"
COL_SS = "secondary_structure"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="merged_shifts_dssp.csv")
    p.add_argument("--snippet", required=True, help="Sequence snippet (2–10 aa)")
    p.add_argument(
        "--include-unknown",
        action="store_true",
        help="Include occurrences with missing SS as 'unknown ss'",
    )
    return p.parse_args()


def validate_snippet(snippet: str) -> str:
    s = snippet.strip().upper()
    if not (2 <= len(s) <= 15):
        raise ValueError("Snippet length must be between 2 and 10 aa.")
    allowed = set("ACDEFGHIKLMNPQRSTVWYBXZ")
    bad = [c for c in s if c not in allowed]
    if bad:
        raise ValueError(f"Invalid amino-acid codes: {bad}")
    return s


def find_all_occurrences(seq: str, sub: str) -> list[int]:
    starts = []
    i = 0
    while True:
        j = seq.find(sub, i)
        if j == -1:
            break
        starts.append(j)
        i = j + 1  # allow overlaps
    return starts


def normalize_ss(raw: str) -> str | None:
    """
    Normalize SS to single-character.
    Coil is '.', as in DSSP-derived data.
    """
    if raw is None:
        return None
    s = str(raw)
    if s.lower() == "nan":
        return None
    c = s[0]
    return c


def main() -> int:
    args = parse_args()
    try:
        snippet = validate_snippet(args.snippet)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2

    try:
        df = pd.read_csv(args.input)
    except Exception as e:
        print(f"Error reading CSV: {e}", file=sys.stderr)
        return 2

    required = [COL_ID, COL_CHAIN, COL_RESNUM, COL_AA, COL_SS]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"Error: missing required columns: {missing}", file=sys.stderr)
        return 2

    df = df.copy()
    df[COL_RESNUM] = pd.to_numeric(df[COL_RESNUM], errors="coerce")
    df = df.dropna(subset=[COL_RESNUM, COL_AA])

    df[COL_AA] = df[COL_AA].astype(str).str.upper().str.strip()

    class_counts = Counter()
    total_occurrences = 0

    for (pid, chain), g in df.groupby([COL_ID, COL_CHAIN], sort=False):
        g = g.sort_values(COL_RESNUM, kind="mergesort")

        seq = "".join(g[COL_AA].tolist())
        ss_raw = g[COL_SS].tolist()

        for start in find_all_occurrences(seq, snippet):
            seg_ss_raw = ss_raw[start : start + len(snippet)]
            seg_ss = []

            missing = False
            for raw in seg_ss_raw:
                c = normalize_ss(raw)
                if c is None:
                    missing = True
                    break
                seg_ss.append(c)

            if missing:
                if args.include_unknown:
                    class_counts["unknown ss"] += 1
                    total_occurrences += 1
                continue

            if all(x == seg_ss[0] for x in seg_ss):
                label = seg_ss[0]
            else:
                label = "mixed ss"

            class_counts[label] += 1
            total_occurrences += 1

    if total_occurrences == 0:
        print(f"No occurrences found for snippet '{snippet}'.")
        return 0

    print(f"Snippet: {snippet}")
    print(f"Occurrences counted: {total_occurrences}")

    for label, count in class_counts.most_common():
        pct = 100.0 * count / total_occurrences
        print(f"{pct:5.1f}% {label}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

