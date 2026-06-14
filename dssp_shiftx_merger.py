#!/usr/bin/env python3
"""
dssp_shiftx_merger.py

Purpose:
    Merge DSSP output (.dssp) with SHIFTX2 backbone-prediction files.

    Input:
      1) DSSP file from xssp (text .dssp, DSSP v3.x fixed-width format)
      2) SHIFTX2 backbone output file:
           - "BACKBONE ATOMS" block format, OR
           - whitespace table format starting with "NUM RES ...", OR
           - a real CSV/TSV

    Output:
      - A DSSP CSV: <dssp_stem>_dssp_summary.csv
      - A merged SHIFTX2 + DSSP CSV: merged_shifts_dssp.csv (or custom name)

    Merge key:
      - Residue index (DSSP col 1–5) <-> SHIFTX Num/NUM/RESID
      - Chain ID (DSSP CHAIN) <-> SHIFTX CHAINID (added if missing)

DSSP columns emitted:
  - entry_id
  - label_asym_id
  - label_seq_id
  - auth_seq_id
  - label_comp_id
  - secondary_structure
  - accessibility
  - phi
  - psi
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


# ---------- DSSP (.dssp) parsing ----------

def parse_dssp_to_df(dssp_path: Path) -> pd.DataFrame:
    """
    Parse a DSSP v3.x .dssp file into a DataFrame.

    Fixed-column layout (standard DSSP v3.x / mkdssp):
      cols  1- 5 : residue index (sequential)                -> label_seq_id
      cols  6-10 : residue number (author numbering)         -> auth_seq_id
      col     12 : chain identifier                          -> label_asym_id
      col     14 : amino-acid code (one-letter)              -> label_comp_id
      col     17 : secondary-structure code (H,E,G,I,B,P,T,S or space/.)
      cols  35-38: ACC (accessible surface area)
      cols ~105-110: PHI (backbone dihedral, degrees)
      cols ~111-116: PSI (backbone dihedral, degrees)

    NOTE: DSSP is fixed-width. We slice PHI/PSI using the same offsets that match
    the mkdssp/xssp output seen in your example file.
    """
    text = dssp_path.read_text()
    header_marker = "#  RESIDUE AA STRUCTURE"
    idx = text.find(header_marker)
    if idx == -1:
        raise RuntimeError(f"DSSP header '{header_marker}' not found in {dssp_path}")

    lines = text[idx:].splitlines()
    data = []
    entry_id = dssp_path.stem

    for line in lines[1:]:
        if not line.strip() or line.startswith("#"):
            continue
        if len(line) < 60:
            continue

        # Parse indices
        try:
            seq_id = int(line[0:5])        # sequential residue index
            auth_seq_id = int(line[5:10])  # residue number (often same)
        except ValueError:
            continue

        chain = line[11].strip() or ""
        aa = line[13].strip() or "X"

        ss_char = line[16].strip()
        secondary_structure = ss_char if ss_char else "."

        try:
            acc = int(line[34:38])
        except ValueError:
            acc = None

        # PHI/PSI: fixed-width slices tuned to mkdssp/xssp format
        # In the example file:
        #   phi is around line[104:110], psi around line[110:116]
        phi = None
        psi = None
        if len(line) >= 116:
            try:
                phi = float(line[104:110].strip())
            except ValueError:
                phi = None
            try:
                psi = float(line[110:116].strip())
            except ValueError:
                psi = None

        data.append(
            (entry_id, chain, seq_id, auth_seq_id, aa, secondary_structure, acc, phi, psi)
        )

    if not data:
        raise RuntimeError(f"No residue records parsed from {dssp_path}")

    return pd.DataFrame(
        data,
        columns=[
            "entry_id",
            "label_asym_id",
            "label_seq_id",
            "auth_seq_id",
            "label_comp_id",
            "secondary_structure",
            "accessibility",
            "phi",
            "psi",
        ],
    )


# ---------- SHIFTX2 parsing ----------

def load_shiftx(path: Path, assume_chain: str = "A", debug: bool = False) -> pd.DataFrame:
    """
    Load a SHIFTX/SHiFTX2 output file.

    Supports:
      1) BACKBONE ATOMS block format:
           BACKBONE ATOMS
           Num  RES   CA   CB   CO   N   H   HA
           ...

      2) SHIFTX2 whitespace table format (looks like CSV but isn't):
           NUM RES   HA     H       N        CA      CB       C
           --- --- ------ ------ -------- ------- ------- --------
           1   M   ...

         This function parses ONLY the first such table block.

      3) Generic delimited CSV/TSV:
         Uses pandas.read_csv(sep=None, engine="python") to infer delimiter.

    Output columns are standardized where possible:
      - RESID (Int64), RESNAME (string), CHAINID (string)
      - plus any shift columns present (CA/CB/CO/N/H/HA etc.)
    """
    lines = path.read_text().splitlines()
    if not lines:
        raise RuntimeError(f"SHIFTX file {path} is empty.")

    # --- Case 1: SHIFTX BACKBONE ATOMS block ---
    if lines[0].strip() == "BACKBONE ATOMS":
        if debug:
            print("[DEBUG] Detected SHIFTX BACKBONE ATOMS format.")

        header_line_idx = 1
        header = lines[header_line_idx].split()

        rows = []
        for ln in lines[header_line_idx + 1:]:
            s = ln.strip()
            if not s:
                break
            parts = s.split()
            if len(parts) != len(header):
                break
            rows.append(parts)

        if not rows:
            raise RuntimeError("BACKBONE ATOMS format detected but no backbone rows parsed.")

        df = pd.DataFrame(rows, columns=header)

        # Convert numeric columns; keep residue name as string
        for col in header:
            if col not in ("RES",):
                df[col] = pd.to_numeric(df[col], errors="coerce")

        df = df.rename(columns={"Num": "RESID", "RES": "RESNAME"})
        if "RESID" not in df.columns or "RESNAME" not in df.columns:
            raise RuntimeError("BACKBONE ATOMS parsing failed: missing Num/RES columns.")

        df["RESID"] = pd.to_numeric(df["RESID"], errors="coerce").astype("Int64")
        df["CHAINID"] = assume_chain

        if debug:
            print(f"[DEBUG] Final SHIFTX columns: {list(df.columns)}")

        return df

    # --- Case 2: “fake csv” whitespace table (NUM/RES or NUM/AA header + dashed separator) ---
    header0 = lines[0].strip()
    header0_upper = header0.upper()
    toks0 = header0_upper.split()
    # SHIFTX2 "chemical shift" tables sometimes label residue as RES, sometimes as AA
    has_num_res_header = ("NUM" in toks0) and (("RES" in toks0) or ("AA" in toks0))
    has_dash_separator = len(lines) > 1 and lines[1].strip().startswith("---")

    if has_num_res_header and has_dash_separator:
        if debug:
            print("[DEBUG] Detected SHIFTX2 whitespace-table format (NUM/RES header + dashed separator).")

        header = lines[0].split()
        rows = []

        # Parse only the first block; later sections may change columns
        for ln in lines[2:]:
            s = ln.strip()
            if not s:
                break

            # If a new section header starts, stop
            # (Conservative: if first char is alpha and line isn't a data row)
            if s and s[0].isalpha():
                break

            parts = s.split()
            if len(parts) != len(header):
                break
            rows.append(parts)

        if not rows:
            raise RuntimeError("Whitespace-table detected but no data rows were parsed.")

        df = pd.DataFrame(rows, columns=header)

        rename_map = {
            "NUM": "RESID",
            "Num": "RESID",
            "RES": "RESNAME",
            "Res": "RESNAME",
            "AA": "RESNAME",  # SHIFTX2 table format often uses AA (1-letter)
            "C": "CO",  # carbonyl often labeled "C" in this format
        }
        df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

        if "RESID" not in df.columns or "RESNAME" not in df.columns:
            raise RuntimeError("Whitespace-table parsing failed: missing NUM/RES (RESID/RESNAME).")

        # Convert numeric columns
        for col in df.columns:
            if col not in ("RESNAME",):
                df[col] = pd.to_numeric(df[col], errors="coerce")

        df["RESID"] = pd.to_numeric(df["RESID"], errors="coerce").astype("Int64")
        df["CHAINID"] = assume_chain

        if debug:
            print(f"[DEBUG] Final SHIFTX columns: {list(df.columns)}")

        return df

    # --- Case 3: Generic delimited CSV/TSV ---
    if debug:
        print("[DEBUG] Using pandas.read_csv (auto-sep) for SHIFTX file (generic delimited).")

    df = pd.read_csv(path, sep=None, engine="python")
    df.columns = [str(c).strip() for c in df.columns]

    # Standardize common column names if present
    if "NUM" in df.columns and "RESID" not in df.columns:
        df = df.rename(columns={"NUM": "RESID"})
    if "Num" in df.columns and "RESID" not in df.columns:
        df = df.rename(columns={"Num": "RESID"})
    if "RES" in df.columns and "RESNAME" not in df.columns:
        df = df.rename(columns={"RES": "RESNAME"})
    if "Res" in df.columns and "RESNAME" not in df.columns:
        df = df.rename(columns={"Res": "RESNAME"})
    if "C" in df.columns and "CO" not in df.columns:
        df = df.rename(columns={"C": "CO"})

    if "CHAINID" not in df.columns:
        df["CHAINID"] = assume_chain

    if "RESID" in df.columns:
        df["RESID"] = pd.to_numeric(df["RESID"], errors="coerce").astype("Int64")

    return df


# ---------- CLI + merge logic ----------

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Merge DSSP (.dssp) secondary-structure data with SHIFTX backbone "
            "chemical-shift predictions."
        )
    )

    p.add_argument("dssp_file", help="Input DSSP file (e.g. model.dssp)")
    p.add_argument("shiftx_file", help="SHIFTX2 output file (BACKBONE ATOMS / whitespace table / CSV)")

    p.add_argument(
        "--dssp-out",
        help="Optional DSSP-only CSV output file. Default: <dssp_stem>_dssp_summary.csv",
    )
    p.add_argument(
        "--merged-out",
        help="Optional merged SHIFTX+DSSP CSV output. Default: merged_shifts_dssp.csv",
    )

    p.add_argument(
        "--dssp-chain-col",
        default="label_asym_id",
        help="DSSP chain column name (default: label_asym_id).",
    )
    p.add_argument(
        "--dssp-res-col",
        default="label_seq_id",
        help="DSSP residue index column name (default: label_seq_id).",
    )
    p.add_argument(
        "--shift-chain-col",
        default="CHAINID",
        help="SHIFTX chain column name (default: CHAINID).",
    )
    p.add_argument(
        "--shift-res-col",
        default="RESID",
        help="SHIFTX residue index column name (default: RESID).",
    )

    p.add_argument(
        "--assume-chain",
        default="A",
        help="If SHIFTX output lacks a chain column, use this chain ID (default: A).",
    )

    p.add_argument(
        "--join-how",
        choices=["left", "inner"],
        default="left",
        help="Merge type: left = keep all SHIFTX rows (default), inner = overlaps only.",
    )
    p.add_argument(
        "--debug",
        action="store_true",
        help="Print extra information about parsing and merging.",
    )

    return p.parse_args()


def main():
    args = parse_args()

    dssp_path = Path(args.dssp_file)
    shift_path = Path(args.shiftx_file)

    if not dssp_path.is_file():
        print(f"Error: DSSP file not found: {dssp_path}", file=sys.stderr)
        sys.exit(1)
    if not shift_path.is_file():
        print(f"Error: SHIFTX file not found: {shift_path}", file=sys.stderr)
        sys.exit(1)

    dssp_out = Path(args.dssp_out) if args.dssp_out else dssp_path.with_name(
        dssp_path.stem + "_dssp_summary.csv"
    )
    merged_out = Path(args.merged_out) if args.merged_out else Path("merged_shifts_dssp.csv")

    # 1) Parse DSSP
    print(f"Parsing DSSP from {dssp_path} ...")
    try:
        dssp_df = parse_dssp_to_df(dssp_path)
    except RuntimeError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    dssp_out.parent.mkdir(parents=True, exist_ok=True)
    dssp_df.to_csv(dssp_out, index=False)
    print(f"Wrote DSSP CSV to {dssp_out} ({len(dssp_df)} rows)")

    # 2) Load SHIFTX
    print(f"Reading SHIFTX file from {shift_path} ...")
    try:
        shift_df = load_shiftx(shift_path, assume_chain=args.assume_chain, debug=args.debug)
    except Exception as e:
        print(f"Error reading SHIFTX file: {e}", file=sys.stderr)
        sys.exit(1)

    # Check required columns exist
    for col in (args.dssp_chain_col, args.dssp_res_col):
        if col not in dssp_df.columns:
            print(
                f"Error: DSSP column '{col}' not found. Available: {list(dssp_df.columns)}",
                file=sys.stderr,
            )
            sys.exit(1)

    for col in (args.shift_chain_col, args.shift_res_col):
        if col not in shift_df.columns:
            print(
                f"Error: SHIFTX column '{col}' not found. Available: {list(shift_df.columns)}",
                file=sys.stderr,
            )
            sys.exit(1)

    # Coerce residue indices to numeric
    dssp_df[args.dssp_res_col] = pd.to_numeric(dssp_df[args.dssp_res_col], errors="coerce")
    shift_df[args.shift_res_col] = pd.to_numeric(shift_df[args.shift_res_col], errors="coerce")

    if args.debug:
        print("\n[DEBUG] Column mapping:")
        print(f"  DSSP chain:   {args.dssp_chain_col}")
        print(f"  DSSP resid:   {args.dssp_res_col}")
        print(f"  SHIFTX chain: {args.shift_chain_col}")
        print(f"  SHIFTX resid: {args.shift_res_col}")
        print(f"[DEBUG] DSSP rows:   {len(dssp_df)}")
        print(f"[DEBUG] SHIFTX rows: {len(shift_df)}\n")
        print(f"[DEBUG] SHIFTX columns: {list(shift_df.columns)}\n")

    # 3) Merge
    merged_df = pd.merge(
        shift_df,
        dssp_df,
        left_on=[args.shift_chain_col, args.shift_res_col],
        right_on=[args.dssp_chain_col, args.dssp_res_col],
        how=args.join_how,
        suffixes=("_shiftx", "_dssp"),
    )

    merged_out.parent.mkdir(parents=True, exist_ok=True)
    merged_df.to_csv(merged_out, index=False)
    print(
        f"Wrote merged SHIFTX+DSSP CSV to {merged_out} "
        f"({len(merged_df)} rows, join='{args.join_how}')"
    )


if __name__ == "__main__":
    main()

