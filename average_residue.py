#!/usr/bin/env python3
"""average_residue_modified.py

Compute per-atom mean and SD for a chosen residue type from:
  1) ShiftX2-style WIDE tables with repeated header blocks (original behavior)
  2) A normal single-header CSV like merged_shifts_dssp.csv (e.g., RESID/RESNAME/CHAINID)

Upgrades:
- Accepts merged-column aliases: RESID/RESNAME/CHAINID (and related mmCIF-ish fields)
- Supports both single-header and multi-header formats automatically
- Optional --multiple-chain to keep chain-specific residue pools separate
"""

import argparse
import csv
import glob
import math
import re
from collections import OrderedDict

AA20 = set("ACDEFGHIKLMNPQRSTVWY")
NUM_RE = re.compile(r"^-?\d+$")

ATOM_ORDER = ["CA", "CB", "CO", "N", "H", "HA", "HB"]
CANONICAL = set(ATOM_ORDER)

HEADER_ALIAS = {
    "CA": "CA",
    "CB": "CB",
    "CO": "CO",
    "C": "CO",
    "C'": "CO",
    "CPRIME": "CO",
    "N": "N",
    "H": "H",
    "HN": "H",
    "H1": "H",
    "HA": "HA",
    "HA1": "HA",
    "HA2": "HA",
    "HA3": "HA",
    "HB": "HB",
    "HB1": "HB",
    "HB2": "HB",
    "HB3": "HB",
}

NUM_COL_ALIASES = [
    "NUM", "RESID", "RESNUM", "RESSEQ", "SEQID",
    "LABEL_SEQ_ID", "AUTH_SEQ_ID",
]
RES_COL_ALIASES = [
    "RES", "RESNAME", "AA", "RESIDUE",
    "LABEL_COMP_ID", "AUTH_COMP_ID",
]
CHAIN_COL_ALIASES = [
    "CHAIN", "CHAINID", "CHAIN_ID",
    "LABEL_ASYM_ID", "AUTH_ASYM_ID",
]


def try_float(x):
    try:
        if x is None:
            return None
        s = str(x).strip()
        if s == "" or s.upper() in {"NA", "NAN", "NULL", "NONE"}:
            return None
        return float(s)
    except Exception:
        return None


def normalize_colname(s: str) -> str:
    s = str(s).strip().lstrip("\ufeff")
    s = s.replace("’", "'").replace("`", "'").replace("′", "'").replace('"', "")
    return s.upper().replace(" ", "").replace("-", "").replace("_", "")


def detect_delim(sample_line: str):
    if sample_line.count(",") >= max(sample_line.count("\t"), sample_line.count(";"), 1):
        return ","
    if sample_line.count("\t") > 1:
        return "\t"
    if sample_line.count(";") > 1:
        return ";"
    return None


def _pick_first_present(norm_headers, aliases):
    for a in aliases:
        na = normalize_colname(a)
        if na in norm_headers:
            return norm_headers.index(na)
    return None


def parse_single_header_table(path, multiple_chain=False):
    """Parse a standard delimited table with ONE header row (covers merged_shifts_dssp.csv)."""
    with open(path, "r", encoding="utf-8", errors="ignore", newline="") as f:
        first = f.readline()
        if not first:
            return [], []
        delim = detect_delim(first) or ","
        f.seek(0)

        reader = csv.reader(f, delimiter=delim)
        header = next(reader, None)
        if not header:
            return [], []

        norm_header = [normalize_colname(h) for h in header]

        idx_num = _pick_first_present(norm_header, NUM_COL_ALIASES)
        idx_res = _pick_first_present(norm_header, RES_COL_ALIASES)
        idx_chain = _pick_first_present(norm_header, CHAIN_COL_ALIASES) if multiple_chain else None

        idx_to_atom = {}
        atoms_seen = set()
        for j, name in enumerate(header):
            canon = HEADER_ALIAS.get(normalize_colname(name))
            if canon in CANONICAL:
                idx_to_atom[j] = canon

        if idx_num is None or idx_res is None:
            # Not a compatible single-header file
            return None, None

        merged = OrderedDict()

        for row in reader:
            if not row or len(row) <= max(idx_num, idx_res):
                continue

            num_raw = row[idx_num]
            if not NUM_RE.match(str(num_raw).strip()):
                continue
            resnum = int(str(num_raw).strip())

            aa = str(row[idx_res]).strip().upper()
            if aa not in AA20:
                continue

            if multiple_chain and idx_chain is not None and idx_chain < len(row):
                chain = str(row[idx_chain]).strip() or "?"
                key = (chain, resnum, aa)
                rec = merged.setdefault(key, {"Chain": chain, "Num": resnum, "RES": aa})
            else:
                key = (resnum, aa)
                rec = merged.setdefault(key, {"Num": resnum, "RES": aa})

            for j, canon in idx_to_atom.items():
                if j < len(row):
                    v = try_float(row[j])
                    if v is not None:
                        rec[canon] = v
                        atoms_seen.add(canon)

        rows = list(merged.values())
        atoms_present = [a for a in ATOM_ORDER if a in atoms_seen]
        return rows, atoms_present


def tokenize(line, delim):
    line = line.lstrip("\ufeff")
    if delim is None:
        return line.split()
    return next(csv.reader([line], delimiter=delim))


def has_cols(tokens, *cols):
    n = [normalize_colname(t) for t in tokens]
    return all(normalize_colname(c) in n for c in cols)


def index_of(tokens, col):
    return [normalize_colname(t) for t in tokens].index(normalize_colname(col))


def parse_multi_header_file(path, multiple_chain=False):
    """Original behavior: scan for ShiftX2-style repeated header blocks."""
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        raw = [ln.rstrip("\n") for ln in f]

    merged = OrderedDict()
    atoms_seen = set()
    i = 0
    n = len(raw)

    while i < n:
        line = raw[i]
        header, delim = None, None

        parts_ws = line.split()
        if has_cols(parts_ws, "NUM") and (has_cols(parts_ws, "RES") or has_cols(parts_ws, "RESNAME")):
            header, delim = parts_ws, None
        else:
            d = detect_delim(line)
            if d:
                parts = tokenize(line, d)
                if has_cols(parts, "NUM") and (has_cols(parts, "RES") or has_cols(parts, "RESNAME")):
                    header, delim = parts, d

        if header is None:
            i += 1
            continue

        i += 1

        idx_num = index_of(header, "NUM")
        idx_res = index_of(header, "RES") if has_cols(header, "RES") else index_of(header, "RESNAME")

        idx_chain = None
        if multiple_chain and idx_num > 0:
            idx_chain = idx_num - 1

        idx_to_atom = {}
        for j, name in enumerate(header):
            if j <= idx_res:
                continue
            canon = HEADER_ALIAS.get(normalize_colname(name))
            if canon in CANONICAL:
                idx_to_atom[j] = canon

        while i < n:
            line = raw[i]
            tokens = tokenize(line, delim)

            if tokens and (has_cols(tokens, "NUM") and (has_cols(tokens, "RES") or has_cols(tokens, "RESNAME"))):
                break

            if not tokens or len(tokens) <= max(idx_num, idx_res):
                i += 1
                continue

            if not NUM_RE.match(tokens[idx_num]):
                i += 1
                continue

            resnum = int(tokens[idx_num])
            aa = str(tokens[idx_res]).strip().upper()
            if aa not in AA20:
                i += 1
                continue

            if multiple_chain and idx_chain is not None and idx_chain < len(tokens):
                chain = str(tokens[idx_chain]).strip() or "?"
                key = (chain, resnum, aa)
                rec = merged.setdefault(key, {"Chain": chain, "Num": resnum, "RES": aa})
            else:
                key = (resnum, aa)
                rec = merged.setdefault(key, {"Num": resnum, "RES": aa})

            for j, canon in idx_to_atom.items():
                if j < len(tokens):
                    v = try_float(tokens[j])
                    if v is not None:
                        atoms_seen.add(canon)
                        rec[canon] = v

            i += 1

    rows = list(merged.values())
    atoms_present = [a for a in ATOM_ORDER if a in atoms_seen]
    return rows, atoms_present


def clean_vals(values):
    return [v for v in values if v is not None]


def stats(values):
    vals = clean_vals(values)
    n = len(vals)
    if n == 0:
        return 0, None, None
    m = sum(vals) / n
    if n == 1:
        return 1, m, 0.0
    var = sum((x - m) ** 2 for x in vals) / (n - 1)
    return n, m, math.sqrt(var)


def compute_residue_stats(rows, atoms, residue, chain=None):
    out = {}
    for a in atoms:
        if chain is None:
            series = [r.get(a) for r in rows if r.get("RES") == residue]
        else:
            series = [r.get(a) for r in rows if r.get("RES") == residue and r.get("Chain") == chain]
        n, m, sd = stats(series)
        out[a] = {"n": n, "mean": m, "sd": sd}
    return out


def parse_any(path, multiple_chain=False):
    rows, atoms = parse_single_header_table(path, multiple_chain=multiple_chain)
    if rows is not None:
        return rows, atoms
    return parse_multi_header_file(path, multiple_chain=multiple_chain)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--residue", "-r", required=True, help="1-letter residue code, e.g. M")
    ap.add_argument(
        "--multiple-chain",
        action="store_true",
        help="Treat each chain independently instead of merging chains.",
    )
    args = ap.parse_args()

    res = args.residue.strip().upper()
    if res not in AA20:
        raise SystemExit(f"Residue must be a 1-letter AA code (got {res!r}).")

    filelist = []
    for pat in args.files:
        filelist.extend(glob.glob(pat))
    if not filelist:
        raise SystemExit("No input files matched.")

    for fname in filelist:
        print(f"\n===== Processing: {fname} =====")
        rows, atoms = parse_any(fname, multiple_chain=args.multiple_chain)

        if not rows:
            print("No rows parsed (file format/headers not recognized or empty).")
            continue
        if not atoms:
            print("Parsed residues, but found no recognizable atom shift columns.")
            continue

        if args.multiple_chain:
            chains = sorted({r.get("Chain", "?") for r in rows if r.get("Chain") is not None})
            for ch in chains:
                n_res = sum(1 for r in rows if r.get("RES") == res and r.get("Chain") == ch)
                print(f"\nChain {ch}: (Number of {res} residues averaged: {n_res})")
                stats_dict = compute_residue_stats(rows, atoms, res, chain=ch)
                print(f"Stats for residue {res} (chain {ch})")
                for atom in atoms:
                    n = stats_dict[atom]["n"]
                    m = stats_dict[atom]["mean"]
                    sd = stats_dict[atom]["sd"]
                    if m is None:
                        print(f"  {atom}: mean=NA sd=NA (n=0)")
                    else:
                        print(f"  {atom}: mean={m:.3f} sd={sd:.3f} (n={n})")
        else:
            n_residues = sum(1 for r in rows if r.get("RES") == res)
            print(f"(Number of {res} residues averaged: {n_residues})")
            stats_dict = compute_residue_stats(rows, atoms, res)
            print(f"\nStats for residue {res}")
            for atom in atoms:
                n = stats_dict[atom]["n"]
                m = stats_dict[atom]["mean"]
                sd = stats_dict[atom]["sd"]
                if m is None:
                    print(f"  {atom}: mean=NA sd=NA (n=0)")
                else:
                    print(f"  {atom}: mean={m:.3f} sd={sd:.3f} (n={n})")


if __name__ == "__main__":
    main()
