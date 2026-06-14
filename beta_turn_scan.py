#!/usr/bin/env python3
"""beta_turn_scan.py

Scan for beta-turn types using PHI/PSI in 4-residue windows (i..i+3).

Input formats supported:
  1) CSI-like whitespace table:
     RESNUM AA3 SST SSS RCI PHI PSI

  2) Merged CSV (your current format):
     RESID, RESNAME, ..., CHAINID, secondary_structure, accessibility, phi, psi

Turn type is determined by PHI/PSI of residues i+1 and i+2, matched to
canonical targets within an angle tolerance (default ±30°).

Canonical targets (deg):
  Type I:   i+1 (-60, -30),  i+2 (-90,   0)
  Type II:  i+1 (-60, 120),  i+2 ( 80,   0)
  Type I':  i+1 ( 60,  30),  i+2 ( 90,   0)
  Type II': i+1 ( 60,-120),  i+2 (-80,   0)

Notes:
  * DSSP may output 360.0 for undefined dihedrals (often at termini). This
    script treats 360.0 as missing (NA).

Usage:
  python beta_turn_scan.py merged_shifts_dssp.csv
  python beta_turn_scan.py merged_shifts_dssp.csv --tol 30 --out beta_turns.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple


@dataclass
class Residue:
    resnum: int
    resname: str
    chain: str
    phi: Optional[float]
    psi: Optional[float]
    dssp: Optional[str] = None
    accessibility: Optional[float] = None


PATTERNS: Dict[str, Tuple[float, float, float, float]] = {
    "Type I":   (-60.0,  -30.0,  -90.0,   0.0),
    "Type II":  (-60.0,  120.0,   80.0,   0.0),
    "Type I'":  ( 60.0,   30.0,   90.0,   0.0),
    "Type II'": ( 60.0, -120.0,  -80.0,   0.0),
}


def parse_float(x: object) -> Optional[float]:
    if x is None:
        return None
    s = str(x).strip()
    if not s or s.upper() == "NA":
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    # DSSP "missing" convention for dihedrals is often 360.0
    if v == 360.0:
        return None
    return v


def norm_angle(a: float) -> float:
    return (a + 180.0) % 360.0 - 180.0


def ang_close(obs: float, ref: float, tol: float) -> bool:
    o = norm_angle(obs)
    r = norm_angle(ref)
    diff = abs(o - r)
    diff = min(diff, 360.0 - diff)
    return diff <= tol


def wrapped_diff(obs: float, ref: float) -> float:
    o = norm_angle(obs)
    r = norm_angle(ref)
    d = abs(o - r)
    return min(d, 360.0 - d)


def classify_turn(phi1: float, psi1: float, phi2: float, psi2: float, tol: float) -> str:
    matches: List[Tuple[str, float]] = []
    for name, (r_phi1, r_psi1, r_phi2, r_psi2) in PATTERNS.items():
        ok = (
            ang_close(phi1, r_phi1, tol)
            and ang_close(psi1, r_psi1, tol)
            and ang_close(phi2, r_phi2, tol)
            and ang_close(psi2, r_psi2, tol)
        )
        if ok:
            score = max(
                wrapped_diff(phi1, r_phi1),
                wrapped_diff(psi1, r_psi1),
                wrapped_diff(phi2, r_phi2),
                wrapped_diff(psi2, r_psi2),
            )
            matches.append((name, score))

    if not matches:
        return "no beta turn"

    matches.sort(key=lambda x: x[1])
    return matches[0][0]


def detect_format(path: str) -> str:
    if path.lower().endswith(".csv"):
        return "csv"
    return "csi"


def load_csi_table(lines: Iterable[str]) -> List[Residue]:
    # RESNUM AA3 SST SSS RCI PHI PSI
    by_key: Dict[Tuple[str, int], Residue] = {}
    for ln in lines:
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split()
        if len(parts) < 7:
            continue
        try:
            resnum = int(parts[0])
        except ValueError:
            continue
        resname = parts[1]
        phi = parse_float(parts[5])
        psi = parse_float(parts[6])
        chain = "?"
        by_key[(chain, resnum)] = Residue(resnum=resnum, resname=resname, chain=chain, phi=phi, psi=psi)

    residues = list(by_key.values())
    residues.sort(key=lambda r: (r.chain, r.resnum))
    return residues


def _find_col(fieldnames: List[str], candidates: List[str]) -> Optional[str]:
    lower = {f.lower(): f for f in fieldnames}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    return None


def load_merged_csv(path: str) -> List[Residue]:
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError("CSV has no header row.")

        fn = reader.fieldnames
        col_resid = _find_col(fn, ["RESID", "resid", "residue", "resnum"])
        col_resname = _find_col(fn, ["RESNAME", "resname", "aa", "amino_acid"])
        col_chain = _find_col(fn, ["CHAINID", "chainid", "chain", "label_asym_id"])
        # Optional dihedrals. Some merged files may not include these.
        col_phi = _find_col(fn, ["PHI", "phi"])
        col_psi = _find_col(fn, ["PSI", "psi"])
        col_ss = _find_col(fn, ["secondary_structure", "dssp", "ss"])
        col_acc = _find_col(fn, ["accessibility", "acc", "asa"])

        if not col_resid or not col_resname:
            raise ValueError(f"CSV missing required columns. Found: {fn}")

        if not col_chain:
            col_chain = None

        # Stash presence for downstream messaging
        load_merged_csv.has_phi_psi = bool(col_phi and col_psi)  # type: ignore[attr-defined]

        # de-duplicate by (chain, resid): keep last occurrence
        by_key: Dict[Tuple[str, int], Residue] = {}

        for row in reader:
            try:
                resnum = int(row[col_resid])
            except Exception:
                continue

            resname = str(row[col_resname]).strip()
            chain = str(row[col_chain]).strip() if col_chain and row.get(col_chain) else "?"

            phi = parse_float(row[col_phi]) if col_phi and (col_phi in row) else None
            psi = parse_float(row[col_psi]) if col_psi and (col_psi in row) else None
            dssp = str(row[col_ss]).strip() if col_ss and row.get(col_ss) is not None else None
            acc = parse_float(row[col_acc]) if col_acc else None

            by_key[(chain, resnum)] = Residue(
                resnum=resnum,
                resname=resname,
                chain=chain,
                phi=phi,
                psi=psi,
                dssp=dssp,
                accessibility=acc,
            )

        residues = list(by_key.values())
        residues.sort(key=lambda r: (r.chain, r.resnum))
        return residues


def format_window(w: List[Residue]) -> str:
    return "-".join(f"{r.resname}{r.resnum}" for r in w)


def scan_chain(chain_res: List[Residue], tol: float) -> List[Dict[str, object]]:
    out: List[Dict[str, object]] = []
    if len(chain_res) < 4:
        return out

    for k in range(len(chain_res) - 3):
        w = chain_res[k:k+4]
        r_i1 = w[1]
        r_i2 = w[2]

        if r_i1.phi is None or r_i1.psi is None or r_i2.phi is None or r_i2.psi is None:
            call = "no beta turn"
        else:
            call = classify_turn(r_i1.phi, r_i1.psi, r_i2.phi, r_i2.psi, tol)

        out.append({
            "chain": w[0].chain,
            "window": format_window(w),
            "i": w[0].resnum,
            "i+1": w[1].resnum,
            "i+2": w[2].resnum,
            "i+3": w[3].resnum,
            "phi(i+1)": r_i1.phi,
            "psi(i+1)": r_i1.psi,
            "phi(i+2)": r_i2.phi,
            "psi(i+2)": r_i2.psi,
            "ss(i)": w[0].dssp,
            "ss(i+1)": w[1].dssp,
            "ss(i+2)": w[2].dssp,
            "ss(i+3)": w[3].dssp,
            "call": call,
        })

    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="Merged CSV (recommended) or CSI-like table")
    ap.add_argument("--tol", type=float, default=30.0, help="Angle tolerance in degrees (default 30)")
    ap.add_argument("--out", default=None, help="Optional output CSV path")
    args = ap.parse_args()

    fmt = detect_format(args.input)

    if args.input == "-":
        residues = load_csi_table(sys.stdin.read().splitlines())
    elif fmt == "csv":
        residues = load_merged_csv(args.input)
        has_phi_psi = getattr(load_merged_csv, "has_phi_psi", False)
    else:
        with open(args.input, "r", encoding="utf-8", errors="replace") as f:
            residues = load_csi_table(f.readlines())

    # group by chain
    chains: Dict[str, List[Residue]] = {}
    for r in residues:
        chains.setdefault(r.chain, []).append(r)

    # output rows
    all_rows: List[Dict[str, object]] = []

    print(f"# Beta-turn scan (tolerance = ±{args.tol:.0f}°)")
    print("# Turn type determined by PHI/PSI of residues i+1 and i+2 in each 4-residue window")
    if fmt == "csv" and not has_phi_psi:
        print("# NOTE: This CSV does not include phi/psi columns, so windows will be labeled 'no beta turn'.")
        print("#       Re-run your merger with phi/psi enabled (DSSP dihedrals) to classify turn types.")
    print()

    print(
        f"{'chain':<5}  {'i..i+3 window':<36}  {'i':>4}  {'i+1':>6}  {'i+2':>6}  {'i+3':>6}  "
        f"{'phi(i+1)':>9}  {'psi(i+1)':>9}  {'phi(i+2)':>9}  {'psi(i+2)':>9}  "
        f"{'ss(i)':>5}  {'ss(i+1)':>7}  {'ss(i+2)':>7}  {'ss(i+3)':>7}  call"
    )
    print("-" * 155)

    def fmtf(x: Optional[float]) -> str:
        return "NA" if x is None else f"{x:.1f}"

    for chain_id, chain_res in sorted(chains.items(), key=lambda x: x[0]):
        chain_res.sort(key=lambda r: r.resnum)
        rows = scan_chain(chain_res, args.tol)
        all_rows.extend(rows)

        for row in rows:
            print(
                f"{row['chain']:<5}  {row['window']:<36}  {row['i']:>4}  {row['i+1']:>6}  {row['i+2']:>6}  {row['i+3']:>6}  "
                f"{fmtf(row['phi(i+1)']):>9}  {fmtf(row['psi(i+1)']):>9}  {fmtf(row['phi(i+2)']):>9}  {fmtf(row['psi(i+2)']):>9}  "
                f"{(row['ss(i)'] or 'NA'):>5}  {(row['ss(i+1)'] or 'NA'):>7}  {(row['ss(i+2)'] or 'NA'):>7}  {(row['ss(i+3)'] or 'NA'):>7}  "
                f"{row['call']}"
            )

    if args.out:
        fieldnames = [
            "chain", "window", "i", "i+1", "i+2", "i+3",
            "phi(i+1)", "psi(i+1)", "phi(i+2)", "psi(i+2)",
            "ss(i)", "ss(i+1)", "ss(i+2)", "ss(i+3)",
            "call",
        ]
        with open(args.out, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for r in all_rows:
                w.writerow(r)
        print(f"\n# Wrote: {args.out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
