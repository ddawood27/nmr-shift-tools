#!/usr/bin/env python3
import re
import sys
from pathlib import Path

OUT_HEADER = "#NUM AA HA CA CB CO N HN"

SIDECHAIN_MARKERS = ("SIDECHAIN PROTON", "SIDECHAIN CARBON")
BACKBONE_MARKER = "BACKBONE ATOMS"

def replace_stars_with_zeros(s: str) -> str:
    # Preserve width: "****" -> "0000", "***" -> "000", etc.
    return re.sub(r"\*+", lambda m: "0" * len(m.group(0)), s)

def strip_backbone_and_sidechains(lines):
    # Remove BACKBONE ATOMS line(s)
    cleaned = [ln for ln in lines if ln.strip().upper() != BACKBONE_MARKER]
    # Cut off at first SIDECHAIN marker
    for i, ln in enumerate(cleaned):
        up = ln.strip().upper()
        if any(up.startswith(m) for m in SIDECHAIN_MARKERS):
            return cleaned[:i]
    return cleaned

def find_backbone_header(lines):
    """
    Finds the header line for the backbone table.
    Accepts:
      - 'Num RES CA CB CO N H HA'
      - 'C Num RES CA CB CO N H HA'  (leading C is chain column label)
    """
    for i, ln in enumerate(lines):
        up = ln.strip().upper()
        if not up:
            continue

        # Must contain these backbone columns
        if ("NUM" in up and "RES" in up and "CA" in up and "CB" in up and "CO" in up and "N" in up and "H" in up and "HA" in up):
            # avoid matching sidechain headers etc by requiring Num + RES
            if re.search(r"\bNUM\b", up) and re.search(r"\bRES\b", up):
                return i
    return None

def is_chain(tok: str) -> bool:
    return re.fullmatch(r"[A-Za-z0-9]", tok) is not None

def is_int(tok: str) -> bool:
    return re.fullmatch(r"\d+", tok) is not None

def format_aligned(rows):
    # rows are [NUM, AA, HA, CA, CB, CO, N, HN]
    headers = ["#NUM", "AA", "HA", "CA", "CB", "CO", "N", "HN"]
    widths = [len(h) for h in headers]
    for r in rows:
        for j, v in enumerate(r):
            widths[j] = max(widths[j], len(v))

    out = [OUT_HEADER]
    for r in rows:
        out.append(" ".join(r[j].rjust(widths[j]) for j in range(8)))
    return "\n".join(out) + "\n"

def process_file(in_path: Path, out_path: Path):
    with open(in_path, "r", newline="") as f:
        lines = f.readlines()

    # Replace stars everywhere first
    lines = [replace_stars_with_zeros(ln) for ln in lines]

    # Remove backbone marker + cut off sidechain sections
    lines = strip_backbone_and_sidechains(lines)

    hdr_idx = find_backbone_header(lines)
    if hdr_idx is None:
        raise ValueError("Could not find backbone header line (Num/RES/CA/CB/CO/N/H/HA).")

    # Parse data lines after header:
    # Each row is either:
    #   A  1  A  CA  CB  CO  N  H  HA   (9 tokens, chain included)
    # or:
    #   1  A  CA  CB  CO  N  H  HA      (8 tokens, no chain)
    rows = []
    for ln in lines[hdr_idx + 1:]:
        if not ln.strip():
            continue

        # stop if we hit any other section header
        up = ln.strip().upper()
        if up.startswith("SIDECHAIN PROTON") or up.startswith("SIDECHAIN CARBON"):
            break

        toks = ln.split()
        if not toks:
            continue

        # with chain
        if len(toks) >= 9 and is_chain(toks[0]) and is_int(toks[1]):
            # drop chain token
            toks = toks[1:]

        # now expect: Num RES CA CB CO N H HA  (8 tokens)
        if len(toks) < 8 or not is_int(toks[0]):
            # ignore non-data lines
            continue

        toks = toks[:8]
        num, aa, ca, cb, co, n, hn, ha = toks

        # output needed: NUM AA HA CA CB CO N HN
        rows.append([num, aa, ha, ca, cb, co, n, hn])

    if not rows:
        raise ValueError("Found backbone header but parsed 0 data rows.")

    out_text = format_aligned(rows)

    with open(out_path, "w", newline="\n") as f:
        f.write(out_text)

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {Path(sys.argv[0]).name} <input.cs> <output.shifty>", file=sys.stderr)
        sys.exit(2)

    process_file(Path(sys.argv[1]), Path(sys.argv[2]))

if __name__ == "__main__":
    main()
