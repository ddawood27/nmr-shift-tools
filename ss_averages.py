#!/usr/bin/env python3
"""
ss_averages.py

Purpose:
    Average a numeric column (e.g. CA, CB, CO, N, H, HA, phi, psi, accessibility)
    for a specified amino acid, split by secondary-structure state.

    For a given amino acid and secondary-structure symbol, this script computes:
      1) mean over all residues of that amino acid
      2) mean over residues of that amino acid IN the specified SS
      3) mean over residues of that amino acid NOT in the specified SS

    Designed to work directly with merged_shifts_dssp.csv produced by cif_to_csv.py,
    but will accept any CSV that has:
      - an amino-acid column
      - a secondary-structure column
      - a numeric value column to average
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Average a numeric value by amino acid and secondary structure."
    )

    p.add_argument(
        "csv_file",
        help="Input CSV file (e.g. merged_shifts_dssp.csv)",
    )

    p.add_argument(
        "--target-aa",
        required=True,
        help=(
            "Amino acid to analyze (e.g. ALA, GLY, A, G). "
            "Will be compared case-insensitively to the amino-acid column."
        ),
    )

    p.add_argument(
        "--target-ss",
        required=True,
        help=(
            "Secondary-structure state to filter on (e.g. H, E, T, .). "
            "Will be compared case-insensitively to the secondary-structure column."
        ),
    )

    p.add_argument(
        "--value-col",
        required=True,
        help=(
            "Name of the numeric column to average "
            "(e.g. CA, CB, CO, N, H, HA, phi, psi, accessibility)."
        ),
    )

    p.add_argument(
        "--aa-col",
        default="RESNAME",
        help=(
            "Name of the amino-acid column (default: RESNAME). "
            "In merged_shifts_dssp.csv this is one-letter codes."
        ),
    )

    p.add_argument(
        "--ss-col",
        default="secondary_structure",
        help="Name of the secondary-structure column (default: secondary_structure).",
    )

    p.add_argument(
        "--out-csv",
        help="Optional: write the summary (3 rows) to this CSV file.",
    )

    p.add_argument(
        "--debug",
        action="store_true",
        help="Print extra info about filters and counts.",
    )

    return p.parse_args()


def main():
    args = parse_args()

    csv_path = Path(args.csv_file)
    if not csv_path.is_file():
        print(f"Error: input file not found: {csv_path}", file=sys.stderr)
        sys.exit(1)

    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"Error reading CSV: {e}", file=sys.stderr)
        sys.exit(1)

    # Check columns
    for col in [args.aa_col, args.ss_col, args.value_col]:
        if col not in df.columns:
            print(f"Error: column '{col}' not found in CSV.", file=sys.stderr)
            print(f"Available columns: {list(df.columns)}", file=sys.stderr)
            sys.exit(1)

    aa_col = args.aa_col
    ss_col = args.ss_col
    val_col = args.value_col

    # Normalize case for comparisons
    target_aa = str(args.target_aa).upper()
    target_ss = str(args.target_ss).upper()

    aa_series = df[aa_col].astype(str).str.upper()
    ss_series = df[ss_col].astype(str).str.upper()

    # Coerce numeric column
    df[val_col] = pd.to_numeric(df[val_col], errors="coerce")

    # Masks
    mask_aa = aa_series == target_aa
    mask_ss = ss_series == target_ss

    df_aa = df[mask_aa].copy()
    df_aa_ss = df_aa[mask_ss[mask_aa]].copy()
    df_aa_not_ss = df_aa[~mask_ss[mask_aa]].copy()

    # Drop NaN in the value column
    vals_all = df_aa[val_col].dropna()
    vals_ss = df_aa_ss[val_col].dropna()
    vals_not_ss = df_aa_not_ss[val_col].dropna()

    if args.debug:
        print(f"Total rows in file: {len(df)}")
        print(f"Rows with target AA ({target_aa}): {len(df_aa)}")
        print(f"Rows with target AA in SS '{target_ss}': {len(df_aa_ss)}")
        print(f"Rows with target AA NOT in SS '{target_ss}': {len(df_aa_not_ss)}")

    def safe_mean(x):
        return float(x.mean()) if len(x) > 0 else float("nan")

    def safe_std(x):
        return float(x.std(ddof=1)) if len(x) > 1 else float("nan")

    summary_rows = [
        {
            "subset": "all_target_aa",
            "target_aa": target_aa,
            "target_ss": target_ss,
            "value_col": val_col,
            "n": int(len(vals_all)),
            "mean": safe_mean(vals_all),
            "std": safe_std(vals_all),
        },
        {
            "subset": "in_target_ss",
            "target_aa": target_aa,
            "target_ss": target_ss,
            "value_col": val_col,
            "n": int(len(vals_ss)),
            "mean": safe_mean(vals_ss),
            "std": safe_std(vals_ss),
        },
        {
            "subset": "not_target_ss",
            "target_aa": target_aa,
            "target_ss": target_ss,
            "value_col": val_col,
            "n": int(len(vals_not_ss)),
            "mean": safe_mean(vals_not_ss),
            "std": safe_std(vals_not_ss),
        },
    ]

    summary_df = pd.DataFrame(summary_rows)

    # Pretty print
    print()
    print(f"Summary for AA={target_aa}, SS='{target_ss}', value={val_col}")
    print(summary_df.to_string(index=False))
    print()

    # Optional CSV
    if args.out_csv:
        out_path = Path(args.out_csv)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(out_path, index=False)
        print(f"Wrote summary to {out_path}")


if __name__ == "__main__":
    main()

