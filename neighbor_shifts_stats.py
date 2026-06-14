#!/usr/bin/env python3

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path
from statistics import mean, stdev
from typing import Dict, Tuple, List, Optional


AA3_TO_AA1 = {
    "ALA": "A","ARG": "R","ASN": "N","ASP": "D","CYS": "C",
    "GLN": "Q","GLU": "E","GLY": "G","HIS": "H","ILE": "I",
    "LEU": "L","LYS": "K","MET": "M","PHE": "F","PRO": "P",
    "SER": "S","THR": "T","TRP": "W","TYR": "Y","VAL": "V",
}
AA1_SET = set(AA3_TO_AA1.values())
ATOM_RE = re.compile(r"^[A-Z][A-Z0-9]{0,2}$")


def normalize_res(name: str) -> str:
    n = name.strip().upper()
    if len(n) == 1:
        if n not in AA1_SET:
            raise ValueError(f"Unknown residue: {name}")
        return n
    if len(n) == 3:
        return AA3_TO_AA1.get(n, None)
    raise ValueError(f"Invalid residue format: {name}")


def read_merged_csv(path: str) -> Tuple[Dict[int,str], Dict[int,Dict[str,float]]]:
    seq_map = {}
    shift_map = defaultdict(dict)

    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        headers_lower = {h.lower(): h for h in reader.fieldnames}

        idx_col = headers_lower.get("resid") or headers_lower.get("res_id")
        aa_col  = headers_lower.get("resname") or headers_lower.get("res_name")

        if not idx_col or not aa_col:
            raise ValueError("Could not detect RESID/RESNAME columns.")

        atom_cols = [h for h in reader.fieldnames if ATOM_RE.fullmatch(h.upper())]

        for row in reader:
            try:
                rid = int(row[idx_col])
            except:
                continue

            seq_map[rid] = normalize_res(row[aa_col])

            for col in atom_cols:
                val = row[col].strip()
                if not val:
                    continue
                try:
                    shift_map[rid][col.upper()] = float(val)
                except:
                    continue

    return seq_map, shift_map


def read_shiftx_table(path: str) -> Tuple[Dict[int,str], Dict[int,Dict[str,float]]]:
    seq_map = {}
    shift_map = defaultdict(dict)

    with open(path) as f:
        lines = [l for l in f if l.strip() and not l.startswith("#")]

    header = None
    for i, line in enumerate(lines[:50]):
        fields = re.split(r"\s+", line.strip())
        if "NUM" in fields and "RES" in fields:
            header = fields
            start = i + 1
            break

    if header is None:
        raise ValueError("Could not detect SHIFTX header.")

    rid_col = header.index("NUM")
    aa_col  = header.index("RES")

    atom_cols = [(i,h) for i,h in enumerate(header) if ATOM_RE.fullmatch(h) and h not in {"NUM","RES"}]

    for line in lines[start:]:
        fields = re.split(r"\s+", line.strip())
        try:
            rid = int(fields[rid_col])
        except:
            continue

        seq_map[rid] = normalize_res(fields[aa_col])

        for idx, atom in atom_cols:
            if idx >= len(fields):
                continue
            try:
                shift_map[rid][atom] = float(fields[idx])
            except:
                continue

    return seq_map, shift_map


def find_matches(seq_map, target, partner, gap, order):
    matches = []
    ids = sorted(seq_map.keys())
    id_set = set(ids)

    for rid in ids:
        rid2 = rid + gap + 1
        if rid2 not in id_set:
            continue

        a = seq_map[rid]
        b = seq_map[rid2]

        if order == "forward" and a == target and b == partner:
            matches.append(rid)

        elif order == "reverse" and a == partner and b == target:
            matches.append(rid2)

        elif order == "either":
            if a == target and b == partner:
                matches.append(rid)
            elif a == partner and b == target:
                matches.append(rid2)

    return matches


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input_file")
    ap.add_argument("-t","--target", required=True)
    ap.add_argument("-p","--partner", required=True)
    ap.add_argument("-g","--gap", required=True, type=int)
    ap.add_argument("-o","--order", required=True, choices=["forward","reverse","either"])
    args = ap.parse_args()

    target = normalize_res(args.target)
    partner = normalize_res(args.partner)

    if args.input_file.endswith(".csv"):
        seq_map, shift_map = read_merged_csv(args.input_file)
    else:
        seq_map, shift_map = read_shiftx_table(args.input_file)

    matches = find_matches(seq_map, target, partner, args.gap, args.order)

    if not matches:
        print("No motif matches found.")
        return

    values = defaultdict(list)

    for rid in matches:
        for atom, shift in shift_map[rid].items():
            values[atom].append(shift)

    print(f"\nSummary for motif {target}-{args.gap}-{partner} ({args.order})")
    print(f"Matched residues: {len(matches)}\n")

    for atom in sorted(values.keys()):
        vals = values[atom]
        if len(vals) == 1:
            print(f"{atom}: n=1 value={vals[0]:.3f}")
        else:
            print(f"{atom}: n={len(vals)} mean={mean(vals):.3f} stdev={stdev(vals):.3f}")


if __name__ == "__main__":
    main()
