#!/usr/bin/env python3
"""
Analyze CSI output to summarize 4-residue turn motifs (e.g., T1, T2).

Given a CSI output file with columns like:
#RES SEQ SST SSS RCI PHI PSI
153 GLU C C 0.02 NA NA
...

This script finds 4-residue occurrences where the SSS column equals the
requested turn type (e.g., T1) for 4 consecutive residues, extracts the
4-amino-acid sequence, counts occurrences, and reports percentages.

Default behavior treats each turn as a *non-overlapping* block of 4 residues
within a contiguous run of the requested turn type.
Use --overlap to count overlapping 4-residue windows within a run.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
}

def parse_csi(path: str):
    """Yield tuples (res_id:int, aa1:str, sss:str) in file order."""
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                res_id = int(parts[0])
            except ValueError:
                continue
            aa3 = parts[1].upper()
            sss = parts[3]  # turn label column in CSI output
            aa1 = AA3_TO_1.get(aa3, "X")
            yield (res_id, aa1, sss)

def find_turn_motifs(records, turn_type: str, overlap: bool = False):
    """Return list of 4AA motif strings for the given turn_type."""
    motifs: list[str] = []
    recs = list(records)
    i = 0
    n = len(recs)

    while i < n:
        res_id, aa1, sss = recs[i]
        if sss != turn_type:
            i += 1
            continue

        # Start of a contiguous run (must be consecutive residue numbers)
        run_res = [recs[i]]
        j = i + 1
        while j < n:
            prev_res_id = run_res[-1][0]
            res_id_j, aa1_j, sss_j = recs[j]
            if sss_j == turn_type and res_id_j == prev_res_id + 1:
                run_res.append(recs[j])
                j += 1
            else:
                break

        run_len = len(run_res)
        if run_len >= 4:
            if overlap:
                for k in range(0, run_len - 3):
                    motifs.append("".join(a for _, a, _ in run_res[k:k+4]))
            else:
                k = 0
                while k + 4 <= run_len:
                    motifs.append("".join(a for _, a, _ in run_res[k:k+4]))
                    k += 4

        i = j  # continue after the run

    return motifs

def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(
        description="Summarize 4-residue amino-acid motifs for a given CSI turn type (e.g., T1, T2).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("csi_file", help="Path to CSI output file (.out)")
    ap.add_argument("turn_type", help="Turn type to analyze (e.g., T1, T2)")
    ap.add_argument(
        "--overlap",
        action="store_true",
        help="Count overlapping 4-residue windows within a contiguous turn run (sliding window).",
    )
    ap.add_argument(
        "--min-count",
        type=int,
        default=1,
        help="Only print motifs with at least this many occurrences.",
    )
    args = ap.parse_args(argv)

    turn_type = args.turn_type.strip()
    motifs = find_turn_motifs(parse_csi(args.csi_file), turn_type=turn_type, overlap=args.overlap)

    if not motifs:
        print(f"No 4-residue motifs found for turn type '{turn_type}'.")
        return 0

    counts = Counter(motifs)
    total = sum(counts.values())
    mode = "overlapping" if args.overlap else "non-overlapping"

    print("CSI turn motif summary")
    print(f"File       : {args.csi_file}")
    print(f"Turn type  : {turn_type}")
    print(f"Counting   : {mode} 4-residue blocks")
    print(f"Total turns: {total}\n")

    print(f"{'4AA motif':<10} {'Count':>7} {'Percent':>9}")
    print("-" * 28)
    for motif, c in counts.most_common():
        if c < args.min_count:
            continue
        pct = (c / total) * 100.0
        print(f"{motif:<10} {c:>7d} {pct:>8.2f}%")

    return 0

if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
