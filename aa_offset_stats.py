#!/usr/bin/env python3
"""
aa_offset_stats_single.py

Compute the distribution of amino acids found at a signed offset n
away from a target amino acid in a sequence (from tables containing NUM and RES).

- ONE input file only (by design)
- Prints full list to terminal
- No CSV outputs
- Supports repeated header blocks
- Supports CSV / TSV / semicolon / whitespace delim
- Optional: also compute the opposite offset (-n) via --both-directions
"""

import argparse
import csv
import re
from collections import Counter, OrderedDict

AA20 = set("ACDEFGHIKLMNPQRSTVWY")
NUM_RE = re.compile(r"^-?\d+$")


def normalize_colname(s: str) -> str:
    s = str(s).strip().lstrip("\ufeff")
    return s.upper().replace(" ", "").replace("-", "").replace("_", "")


def detect_delim(line: str):
    # Prefer commas if clearly present; else tabs; else semicolons; else whitespace
    if line.count(",") >= max(line.count("\t"), line.count(";"), 1):
        return ","
    if line.count("\t") > 1:
        return "\t"
    if line.count(";") > 1:
        return ";"
    return None  # whitespace


def tokenize(line: str, delim):
    line = line.lstrip("\ufeff")
    if delim is None:
        return line.split()
    return next(csv.reader([line], delimiter=delim))


def has_cols(tokens, *cols):
    n = [normalize_colname(t) for t in tokens]
    return all(c in n for c in cols)


def index_of(tokens, col):
    return [normalize_colname(t) for t in tokens].index(col)


def parse_sequence_multi_header(path: str):
    """
    Parse file potentially containing repeated header blocks.
    Returns an ordered list of (Num, RES) in the order encountered,
    de-duplicated by residue number (keeping first seen).
    """
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        raw = [ln.rstrip("\n") for ln in f]

    seq_by_num = OrderedDict()

    i, n = 0, len(raw)
    while i < n:
        line = raw[i]
        header, delim = None, None

        # detect header with NUM & RES
        parts_ws = line.split()
        if has_cols(parts_ws, "NUM", "RES"):
            header, delim = parts_ws, None
        else:
            d = detect_delim(line)
            if d:
                parts = tokenize(line, d)
                if has_cols(parts, "NUM", "RES"):
                    header, delim = parts, d

        if header is None:
            i += 1
            continue

        # header found
        i += 1
        try:
            idx_num = index_of(header, "NUM")
            idx_res = index_of(header, "RES")
        except ValueError:
            continue

        # read rows until next header
        while i < n:
            line = raw[i]
            tokens = tokenize(line, delim)

            # break if next header encountered
            if tokens and has_cols(tokens, "NUM", "RES") and tokens is not header:
                break

            if not tokens:
                i += 1
                continue

            if len(tokens) <= max(idx_num, idx_res):
                i += 1
                continue

            if not NUM_RE.match(tokens[idx_num]):
                i += 1
                continue

            try:
                resnum = int(tokens[idx_num])
            except Exception:
                i += 1
                continue

            aa = str(tokens[idx_res]).strip().upper()
            if aa not in AA20:
                i += 1
                continue

            # keep first occurrence per residue number
            if resnum not in seq_by_num:
                seq_by_num[resnum] = aa

            i += 1

    return [(num, aa) for num, aa in seq_by_num.items()]


def compute_offset_distribution(seq, target_aa: str, offset: int, require_consecutive: bool):
    """
    seq: list of (Num, RES) ordered by first-seen Num
    offset: signed int, e.g. +2 or -1
    require_consecutive:
      - True: only count if Num(target)+offset exists exactly
      - False: use parsed order adjacency (i+offset)
    """
    counts = Counter()
    total = 0

    if not seq:
        return counts, total

    if require_consecutive:
        num_to_idx = {num: i for i, (num, _) in enumerate(seq)}
        for num, aa in seq:
            if aa != target_aa:
                continue
            partner_num = num + offset
            j = num_to_idx.get(partner_num)
            if j is None:
                continue
            partner_aa = seq[j][1]
            counts[partner_aa] += 1
            total += 1
    else:
        for i, (_, aa) in enumerate(seq):
            if aa != target_aa:
                continue
            j = i + offset
            if j < 0 or j >= len(seq):
                continue
            partner_aa = seq[j][1]
            counts[partner_aa] += 1
            total += 1

    return counts, total


def print_distribution(counts: Counter, total: int, target: str, offset: int):
    direction = f"{offset:+d}"
    print(f"\n== Distribution at offset {direction} from target '{target}' ==")
    print(f"Total matched target sites with a valid partner at offset {direction}: {total}")
    if total == 0:
        print("(No matches under current settings.)")
        return

    items = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    print("\nAA\tCount\tPercent")
    for aa, c in items:
        pct = (c / total) * 100.0
        print(f"{aa}\t{c}\t{pct:6.2f}%")


def main():
    ap = argparse.ArgumentParser(
        description="AA distribution n residues away from a target AA (ONE file only)."
    )
    ap.add_argument("file", help="Single input file (CSV/TSV/whitespace) containing NUM and RES columns.")
    ap.add_argument("--aa", required=True, help="Target amino acid (one-letter), e.g. P")
    ap.add_argument("--n", type=int, required=True,
                    help="Signed offset (e.g. 1 means i+1 toward C-terminus, -1 means i-1 toward N-terminus).")
    ap.add_argument("--both-directions", action="store_true",
                    help="Also compute the opposite offset (-n) in addition to n.")
    ap.add_argument("--no-index-check", action="store_true",
                    help="Do NOT require residue numbers to be consecutive; use parsed row order instead.")
    args = ap.parse_args()

    target = args.aa.strip().upper()
    if target not in AA20:
        raise SystemExit("ERROR: --aa must be one of ACDEFGHIKLMNPQRSTVWY")

    seq = parse_sequence_multi_header(args.file)
    if not seq:
        raise SystemExit("ERROR: No sequence rows detected (need a header containing NUM and RES).")

    require_consecutive = not args.no_index_check

    offsets = [args.n]
    if args.both_directions and args.n != 0:
        offsets.append(-args.n)

    print(f"File: {args.file}")
    print(f"Target AA: {target}")
    print(f"Offsets: {', '.join([f'{o:+d}' for o in offsets])}")
    print(f"Mode: {'NUM-based (require consecutive residue numbers)' if require_consecutive else 'Row-order adjacency (no index check)'}")

    for off in offsets:
        counts, total = compute_offset_distribution(seq, target, off, require_consecutive=require_consecutive)
        print_distribution(counts, total, target, off)


if __name__ == "__main__":
    main()

