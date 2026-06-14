#!/usr/bin/env python3
"""
Compute normalized position of an amino acid within contiguous secondary-structure segments.

Reads a DSSP summary CSV and prints statistics only (no files written).

Normalized position definition:
  - 0.0 = N-terminal edge of SS segment
  - 1.0 = C-terminal edge of SS segment
"""

import argparse
from pathlib import Path
import sys

import pandas as pd


def build_segments(df, chain_col, resid_col, ss_col, ss_target):
    df = df.copy()
    df = df.sort_values([chain_col, resid_col]).reset_index(drop=True)

    is_target = df[ss_col] == ss_target
    df["_is_target"] = is_target

    prev_chain = df[chain_col].shift()
    prev_is = df["_is_target"].shift()
    prev_res = df[resid_col].shift()

    is_consecutive = df[resid_col] == (prev_res + 1)

    segment_break = (
        (df[chain_col] != prev_chain) |
        (df["_is_target"] != prev_is) |
        (df["_is_target"] & (~is_consecutive))
    )

    df["_seg_id"] = segment_break.cumsum()
    df["segment_id"] = df["_seg_id"].where(df["_is_target"], pd.NA)

    df["segment_len"] = (
        df.groupby("segment_id")[resid_col]
        .transform("size")
        .where(df["_is_target"], pd.NA)
    )

    seg_min = (
        df.groupby("segment_id")[resid_col]
        .transform("min")
    )
    df["offset"] = (df[resid_col] - seg_min).where(df["_is_target"], pd.NA)

    def norm_pos(row):
        if pd.isna(row["segment_len"]):
            return pd.NA
        L = int(row["segment_len"])
        if L <= 1:
            return 0.0
        return float(row["offset"]) / float(L - 1)

    df["norm_pos"] = df.apply(norm_pos, axis=1)

    return df.drop(columns=["_is_target", "_seg_id"])


def main():
    ap = argparse.ArgumentParser(
        description="Compute normalized location of an amino acid within a secondary structure (no output files)."
    )
    ap.add_argument("--dssp", required=True, help="DSSP summary CSV")
    ap.add_argument("--aa", required=True, help="Amino acid (one-letter code, e.g. P)")
    ap.add_argument("--ss", required=True, help="Secondary structure code (e.g. E, H)")
    ap.add_argument("--chain-col", default="label_asym_id")
    ap.add_argument("--resid-col", default="label_seq_id")
    ap.add_argument("--aa-col", default="label_comp_id")
    ap.add_argument("--ss-col", default="secondary_structure")
    ap.add_argument("--preview", type=int, default=0,
                    help="Print first N matching rows for inspection (default: 0)")

    args = ap.parse_args()

    path = Path(args.dssp)
    if not path.is_file():
        print(f"Error: file not found: {path}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(path)

    for c in (args.chain_col, args.resid_col, args.aa_col, args.ss_col):
        if c not in df.columns:
            print(f"Error: missing column '{c}'", file=sys.stderr)
            sys.exit(1)

    df[args.chain_col] = df[args.chain_col].fillna("").astype(str)
    df[args.aa_col] = df[args.aa_col].fillna("").astype(str).str.upper().str.strip()
    df[args.ss_col] = df[args.ss_col].fillna("").astype(str).str.strip()
    df[args.resid_col] = pd.to_numeric(df[args.resid_col], errors="coerce")
    df = df.dropna(subset=[args.resid_col])
    df[args.resid_col] = df[args.resid_col].astype(int)

    aa = args.aa.upper()
    ss = args.ss

    df2 = build_segments(df, args.chain_col, args.resid_col, args.ss_col, ss)

    hits = df2[
        (df2[args.aa_col] == aa) &
        (df2[args.ss_col] == ss)
    ]

    if hits.empty:
        print(f"No residues found for AA='{aa}' in SS='{ss}'.")
        return

    norm = pd.to_numeric(hits["norm_pos"], errors="coerce").dropna()
    seglen = pd.to_numeric(hits["segment_len"], errors="coerce").dropna()

    print(f"\nAA = '{aa}' in SS = '{ss}'")
    print(f"Total matches: {len(hits)}")

    print("\nNormalized position (0 = N-edge, 1 = C-edge):")
    print(f"  mean   : {norm.mean():.3f}")
    print(f"  median : {norm.median():.3f}")
    print(f"  std    : {norm.std():.3f}")
    print(f"  min    : {norm.min():.3f}")
    print(f"  max    : {norm.max():.3f}")

    print("\nSegment length:")
    print(f"  mean   : {seglen.mean():.2f}")
    print(f"  median : {seglen.median():.2f}")
    print(f"  min    : {seglen.min():.0f}")
    print(f"  max    : {seglen.max():.0f}")

    if args.preview > 0:
        print(f"\nPreview (first {args.preview} hits):")
        cols = [
            args.chain_col,
            args.resid_col,
            args.aa_col,
            "segment_len",
            "offset",
            "norm_pos",
        ]
        print(hits[cols].head(args.preview).to_string(index=False))


if __name__ == "__main__":
    main()

