#!/usr/bin/env python3
"""
snippet_shift_avg.py

Find exact occurrences of an AA snippet within each chain, then for a specified
position within that snippet, extract chemical shifts and compute mean/SD.

Supports:
  - ShiftX2-style .cs files (multi-header sections; whitespace-delimited)
  - CSV "wide" tables (one row per residue, atom columns present)

Behavior:
  - Snippet matching is always done within each chain separately.
  - At the target position within each matched snippet, we collect shifts for
    ALL available atom columns and report mean/SD per atom.

Notes for .cs:
  - .cs files often have multiple sections. This parser merges them into one
    row per (chain, resnum) by taking the first non-missing value per atom.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

AA_ALLOWED = set("ACDEFGHIKLMNPQRSTVWYBXZ")

# Atoms we will report (if present)
ATOM_ORDER = ["CA", "CB", "CO", "N", "H", "HA", "HB"]

# Column names used internally
COL_ID = "entry_id"
COL_CHAIN = "CHAINID"
COL_RESNUM = "label_seq_id"
COL_AA = "RESNAME"

# Some .cs/exports use alternate labels; normalize to our canonical set
ATOM_ALIAS = {
    "C": "CO",
    "C'": "CO",
    "CPRIME": "CO",
    "HN": "H",
    "H1": "H",
    "HA1": "HA",
    "HA2": "HA",
    "HA3": "HA",
    "HB1": "HB",
    "HB2": "HB",
    "HB3": "HB",
}


def normalize_atom(a: str) -> str:
    x = str(a).strip().upper().replace(" ", "")
    return ATOM_ALIAS.get(x, x)


def validate_snippet(snippet: str) -> str:
    s = snippet.strip().upper()
    if not (2 <= len(s) <= 80):
        raise ValueError("Snippet length must be between 2 and 80 aa.")
    bad = [c for c in s if c not in AA_ALLOWED]
    if bad:
        raise ValueError(f"Invalid amino-acid codes in snippet: {bad}")
    return s


def find_all_occurrences(seq: str, sub: str) -> List[int]:
    """Return start indices (0-based), allowing overlaps."""
    starts: List[int] = []
    i = 0
    while True:
        j = seq.find(sub, i)
        if j == -1:
            break
        starts.append(j)
        i = j + 1
    return starts


def _to_float(x: str) -> Optional[float]:
    x = str(x).strip()
    if x == "" or x == "****":
        return None
    try:
        return float(x)
    except Exception:
        return None


def mean_sd(values: List[float]) -> Tuple[Optional[float], Optional[float], int]:
    vals = []
    for v in values:
        if v is None:
            continue
        try:
            fv = float(v)
        except Exception:
            continue
        if math.isnan(fv):
            continue
        vals.append(fv)

    n = len(vals)
    if n == 0:
        return None, None, 0
    m = sum(vals) / n
    if n == 1:
        return m, 0.0, 1
    var = sum((x - m) ** 2 for x in vals) / (n - 1)
    return m, math.sqrt(var), n


def parse_cs_file(path: str) -> pd.DataFrame:
    """
    Parse ShiftX2-style .cs file into a DataFrame with columns:
      entry_id, CHAINID, label_seq_id, RESNAME, and any atom columns present.

    Handles multiple header blocks by letting duplicates through; we collapse later.
    """
    rows: List[Dict[str, Any]] = []
    current_header: Optional[List[str]] = None

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()

            # header line in .cs contains Num and RES (case-sensitive in your files)
            if "Num" in parts and "RES" in parts:
                current_header = parts
                continue

            if current_header is None:
                continue

            if len(parts) < 3:
                continue

            chain = parts[0].strip()
            try:
                resnum = int(parts[1])
            except Exception:
                continue
            aa = parts[2].strip().upper()
            if aa not in AA_ALLOWED:
                continue

            rec: Dict[str, Any] = {
                COL_ID: os.path.basename(path),
                COL_CHAIN: chain,
                COL_RESNUM: resnum,
                COL_AA: aa,
            }

            # Map header columns to values; normalize atom names
            # Note: in .cs, columns are aligned with parts indices
            for idx, col in enumerate(current_header):
                atom = normalize_atom(col)
                if atom in ATOM_ORDER and idx < len(parts):
                    rec[atom] = _to_float(parts[idx])

            rows.append(rec)

    if not rows:
        return pd.DataFrame(columns=[COL_ID, COL_CHAIN, COL_RESNUM, COL_AA] + ATOM_ORDER)

    return pd.DataFrame(rows)


def collapse_to_one_row_per_residue(df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse duplicates into one row per (entry_id, chain, resnum).
    For each atom column, take the first non-missing value across duplicate rows.
    """
    df = df.copy()
    df[COL_RESNUM] = pd.to_numeric(df[COL_RESNUM], errors="coerce")
    df = df.dropna(subset=[COL_RESNUM, COL_AA])
    df[COL_RESNUM] = df[COL_RESNUM].astype(int)
    df[COL_AA] = df[COL_AA].astype(str).str.upper().str.strip()

    # Ensure atom columns exist (important for CSV inputs)
    for atom in ATOM_ORDER:
        if atom not in df.columns:
            df[atom] = pd.NA

    def first_nonnull(series: pd.Series):
        for v in series:
            if pd.notna(v):
                return v
        return pd.NA

    agg: Dict[str, Any] = {COL_AA: "first"}
    for atom in ATOM_ORDER:
        agg[atom] = first_nonnull

    out = (
        df.sort_values([COL_ID, COL_CHAIN, COL_RESNUM], kind="mergesort")
          .groupby([COL_ID, COL_CHAIN, COL_RESNUM], as_index=False)
          .agg(agg)
    )
    return out


def load_input(path: str) -> pd.DataFrame:
    if path.lower().endswith(".cs"):
        df = parse_cs_file(path)
        return collapse_to_one_row_per_residue(df)

    # CSV assumed: must contain at least chain/resnum/aa columns; atom cols optional
    df = pd.read_csv(path)
    # Require these core columns (user can rename in their CSV to match these if needed)
    missing = [c for c in [COL_ID, COL_CHAIN, COL_RESNUM, COL_AA] if c not in df.columns]
    if missing:
        raise ValueError(
            f"CSV missing required columns {missing}. "
            f"Expected at least: {COL_ID}, {COL_CHAIN}, {COL_RESNUM}, {COL_AA}."
        )
    return collapse_to_one_row_per_residue(df)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="Input file: .cs or CSV")
    p.add_argument("--snippet", required=True, help="Sequence snippet to search (e.g., AAAAAGA)")
    p.add_argument(
        "--position",
        type=int,
        required=True,
        help="1-based position within the snippet to sample (1..len(snippet))",
    )
    p.add_argument(
        "--multiple-chain",
        action="store_true",
        help="Analyze all chains. If not set, only the first chain encountered is used.",
    )
    p.add_argument(
        "--debug",
        action="store_true",
        help="Print occurrence counts per chain to stderr.",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()

    try:
        snippet = validate_snippet(args.snippet)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2

    if args.position < 1 or args.position > len(snippet):
        print(f"Error: --position must be between 1 and {len(snippet)} for snippet '{snippet}'.", file=sys.stderr)
        return 2

    try:
        df = load_input(args.input)
    except Exception as e:
        print(f"Error reading input: {e}", file=sys.stderr)
        return 2

    if df.empty:
        print("No residue data found in input.", file=sys.stderr)
        return 2

    # If not multiple-chain, restrict to first chain (stable)
    if not args.multiple_chain:
        first_chain = str(df[COL_CHAIN].iloc[0])
        df = df[df[COL_CHAIN].astype(str) == first_chain].copy()

    # Collect values per atom for the target position in each snippet match
    values_by_atom: Dict[str, List[float]] = {a: [] for a in ATOM_ORDER}
    total_matches = 0  # matches of the sequence snippet (regardless of missing shifts)

    for (pid, chain), g in df.groupby([COL_ID, COL_CHAIN], sort=False):
        g = g.sort_values(COL_RESNUM, kind="mergesort")
        seq = "".join(g[COL_AA].tolist())

        starts = find_all_occurrences(seq, snippet)
        if args.debug:
            print(f"[debug] {pid} chain {chain}: occurrences={len(starts)}", file=sys.stderr)

        if not starts:
            continue

        # we will index into g rows by sequence index
        # convert atom columns to numeric lists aligned with sequence order
        atom_arrays = {}
        for atom in ATOM_ORDER:
            atom_arrays[atom] = pd.to_numeric(g[atom], errors="coerce").tolist()

        for start in starts:
            total_matches += 1
            idx0 = start + (args.position - 1)
            if idx0 < 0 or idx0 >= len(seq):
                continue
            for atom in ATOM_ORDER:
                v = atom_arrays[atom][idx0]
                # keep only usable numeric values
                try:
                    fv = float(v)
                except Exception:
                    continue
                if math.isnan(fv):
                    continue
                values_by_atom[atom].append(fv)

    if total_matches == 0:
        print(f"No occurrences found for snippet '{snippet}'.")
        return 0

    print(f"Snippet: {snippet}")
    print(f"Position within snippet: {args.position} (1-based)")
    print(f"Sequence matches found: {total_matches}")
    print("\nAtom\tN\tMean\tSD")

    any_reported = False
    for atom in ATOM_ORDER:
        m, sd, n = mean_sd(values_by_atom[atom])
        if n == 0:
            continue
        any_reported = True
        print(f"{atom}\t{n}\t{m:.4f}\t{sd:.4f}")

    if not any_reported:
        print("(No usable chemical shift values found at that position for any atom columns.)")
        print("Tip: try a different position, or confirm the .cs section contains those atom shifts.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
