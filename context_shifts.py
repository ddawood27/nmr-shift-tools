#!/usr/bin/env python3
"""
context_shifts.py — averages for poly-A runs and their flanking alanines.

UPDATED:
- Added --multiple-chain flag.
  * Default (off): old behavior, merges chains (Num, RES).
  * On: treat each chain independently and run the run/flank logic per chain,
        then pool results across chains.
"""

import argparse, glob, math, re, csv
from collections import defaultdict, OrderedDict

# --- Basic constants & helpers ---

AA20 = set("ACDEFGHIKLMNPQRSTVWY")
NUM_RE = re.compile(r"^-?\d+$")
ATOM_ORDER = ["CA", "CB", "CO", "N", "H", "HA", "HB"]
CANONICAL = set(ATOM_ORDER)

# Map various atom labels to canonical names
HEADER_ALIAS = {
    "CA": "CA", "CB": "CB", "CO": "CO",
    "C": "CO", "C'": "CO", "CPRIME": "CO",
    "H": "H", "HN": "H",
    "HA": "HA", "HA1": "HA", "HA2": "HA", "HA3": "HA",
    "HB": "HB", "HB1": "HB", "HB2": "HB", "HB3": "HB",
    "N": "N",
}

def try_float(x):
    try:
        return float(x)
    except Exception:
        return None

def normalize_colname(s: str) -> str:
    s = s.strip().replace("’","'").replace("`","'").replace("′","'").replace('"',"")
    return s.upper().replace(" ", "").replace("-", "").replace("_", "")

def norm_list(tokens):
    return [normalize_colname(str(t)) for t in tokens]

def has_cols(tokens, *cols):
    n = norm_list(tokens)
    return all(c in n for c in cols)

def index_of(tokens, col):
    return norm_list(tokens).index(col)

# --- Delimiter & tokenization helpers ---

def detect_delim(line):
    """
    Guess a delimiter: prefer comma, then tab, then semicolon; else whitespace.
    """
    if line.count(",") >= max(line.count("\t"), line.count(";"), 1):
        return ","
    if line.count("\t") > 1:
        return "\t"
    if line.count(";") > 1:
        return ";"
    return None  # fall back to .split()

def tokenize(line, delim):
    if delim is None:
        return line.split()
    return next(csv.reader([line], delimiter=delim))

# --- Header mapping for WIDE tables ---

def build_mapping_from_header(header_tokens, idx_res, exact_columns=None):
    """
    For WIDE tables: map column indices to canonical atom names.
    If exact_columns is given, only those names are allowed and required.
    """
    idx_to_atom, atom_sources = {}, defaultdict(list)
    if exact_columns:
        wanted = [normalize_colname(c) for c in exact_columns]
        for j, name in enumerate(header_tokens):
            if j <= idx_res:
                continue
            nm = normalize_colname(name)
            if nm in wanted:
                if nm not in CANONICAL:
                    raise SystemExit(f"--exact-columns includes unknown atom '{nm}'.")
                idx_to_atom[j] = nm
                atom_sources[nm].append(nm)
        missing = [c for c in wanted if c not in atom_sources]
        if missing:
            raise SystemExit(f"Missing required columns in this block: {missing}")
    else:
        for j, name in enumerate(header_tokens):
            if j <= idx_res:
                continue
            canon = HEADER_ALIAS.get(normalize_colname(name))
            if canon in CANONICAL:
                idx_to_atom[j] = canon
                atom_sources[canon].append(name)
    return idx_to_atom, atom_sources

def detect_long_cols(header_tokens):
    """
    Detect LONG-table columns: ATOMNAME/ATOM and SHIFT.
    Returns (idx_atomname, idx_shift) or (None, None) if not found.
    """
    norm = norm_list(header_tokens)
    idx_atom = None
    for cand in ("ATOMNAME", "ATOM"):
        if cand in norm:
            idx_atom = norm.index(cand)
            break
    idx_shift = norm.index("SHIFT") if "SHIFT" in norm else None
    return idx_atom, idx_shift

# --- Parser that supports multi-headers, WIDE/LONG, mixed delimiters ---

def parse_multi_header_file(path, debug_columns=False, exact_columns=None, multiple_chain=False):
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        raw = [ln.rstrip("\n") for ln in f]

    rows_by_key = OrderedDict()
    atoms_seen = set()
    i = 0
    block_id = 0
    n = len(raw)

    while i < n:
        line = raw[i]

        # Look for a header line containing NUM and RES (case-insensitive)
        header, delim = None, None
        parts_ws = line.split()
        if has_cols(parts_ws, "NUM", "RES"):
            header, delim = parts_ws, None
        else:
            d = detect_delim(line)
            if d is not None:
                parts = tokenize(line, d)
                if has_cols(parts, "NUM", "RES"):
                    header, delim = parts, d

        if header is None:
            i += 1
            continue

        block_id += 1
        i += 1
        try:
            idx_num = index_of(header, "NUM")
            idx_res = index_of(header, "RES")
        except ValueError:
            continue

        # If multiple_chain, assume "chain" is the column immediately before NUM
        idx_chain = (idx_num - 1) if (multiple_chain and idx_num > 0) else None

        # Try WIDE
        idx_to_atom, atom_sources = build_mapping_from_header(header, idx_res, exact_columns)
        wide_mode = len(idx_to_atom) > 0

        # If not WIDE, try LONG
        idx_atom_long, idx_shift_long = (None, None)
        long_mode = False
        if not wide_mode and not exact_columns:
            idx_atom_long, idx_shift_long = detect_long_cols(header)
            long_mode = (idx_atom_long is not None and idx_shift_long is not None)

        if not wide_mode and not long_mode:
            if debug_columns:
                print(f"[Block {block_id}] No usable atom columns after RES and not a LONG table.")
            # Skip until the next possible header
            while i < n:
                maybe = raw[i]
                d2 = detect_delim(maybe)
                toks = tokenize(maybe, d2) if d2 else maybe.split()
                if has_cols(toks, "NUM", "RES"):
                    break
                i += 1
            continue

        if debug_columns:
            if wide_mode:
                print(f"== Block {block_id}: WIDE mapping (header -> canonical) ==")
                for j in sorted(idx_to_atom):
                    print(f"  {header[j]}  ->  {idx_to_atom[j]}")
                print()
            else:
                print(f"== Block {block_id}: LONG mode ATOM={header[idx_atom_long]} SHIFT={header[idx_shift_long]} ==\n")

        # Parse this block's rows
        while i < n:
            line = raw[i]
            tokens = tokenize(line, delim)
            # New header?
            if tokens and has_cols(tokens, "NUM", "RES") and tokens is not header:
                break
            if not tokens:
                i += 1
                continue

            try:
                if len(tokens) <= max(idx_num, idx_res) or not NUM_RE.match(tokens[idx_num]):
                    i += 1
                    continue
            except Exception:
                i += 1
                continue

            try:
                resnum = int(tokens[idx_num])
            except Exception:
                i += 1
                continue

            aa = str(tokens[idx_res]).upper()
            if aa not in AA20:
                i += 1
                continue

            if multiple_chain and idx_chain is not None and idx_chain < len(tokens):
                chain = str(tokens[idx_chain]).strip()
                key = (chain, resnum, aa)
                if key not in rows_by_key:
                    rows_by_key[key] = {"Chain": chain, "Num": resnum, "RES": aa}
            else:
                key = (resnum, aa)
                if key not in rows_by_key:
                    rows_by_key[key] = {"Num": resnum, "RES": aa}

            rec = rows_by_key[key]

            if wide_mode:
                per_atom = defaultdict(list)
                for j, canon in idx_to_atom.items():
                    if j < len(tokens):
                        v = try_float(tokens[j])
                        if v is not None:
                            per_atom[canon].append(v)
                if exact_columns:
                    for canon, vals in per_atom.items():
                        if len(vals) > 1:
                            raise SystemExit(
                                f"[{path}] Multiple columns found for atom '{canon}' in exact mode; "
                                f"columns = {atom_sources[canon]}"
                            )
                for canon, vals in per_atom.items():
                    atoms_seen.add(canon)
                    val = sum(vals)/len(vals) if vals else None
                    if rec.get(canon) is not None:
                        rec[canon] = (rec[canon] + val) / 2.0 if val is not None else rec[canon]
                    elif val is not None:
                        rec[canon] = val
            else:
                # LONG mode: pivot ATOMNAME/SHIFT rows into per-residue wide fields
                if idx_atom_long >= len(tokens) or idx_shift_long >= len(tokens):
                    i += 1
                    continue
                atom_name = normalize_colname(tokens[idx_atom_long])
                canon = HEADER_ALIAS.get(atom_name)
                if canon in CANONICAL:
                    v = try_float(tokens[idx_shift_long])
                    if v is not None:
                        atoms_seen.add(canon)
                        if rec.get(canon) is None:
                            rec[canon] = v
                        else:
                            rec[canon] = (rec[canon] + v) / 2.0

            i += 1

    rows = list(rows_by_key.values())
    atoms_present = [a for a in ATOM_ORDER if a in atoms_seen]

    # If multiple_chain, keep deterministic order within each chain by Num
    if multiple_chain:
        rows.sort(key=lambda r: (r.get("Chain",""), r["Num"]))

    return rows, atoms_present

# --- Run finding & flanking logic ---

def find_runs(residues, run_aa, min_len):
    """
    Find contiguous runs of run_aa with length >= min_len.
    Returns a list of (start_index, end_index) in 0-based sequence index.
    """
    runs = []
    i = 0
    n = len(residues)
    while i < n:
        if residues[i] != run_aa:
            i += 1
            continue
        j = i
        while j < n and residues[j] == run_aa:
            j += 1
        if j - i >= min_len:
            runs.append((i, j - 1))
        i = j
    return runs

def flanks_any_side_by_distance(residues, runs, flank_aa, flank_dist):
    """
    For each run, compute its *qualifying* flanking positions at distance flank_dist:
      left_index  = s - flank_dist
      right_index = e + flank_dist

    A run is kept if EITHER side exists AND is flank_aa:
      residues[left_index] == flank_aa  (in range)  OR
      residues[right_index] == flank_aa (in range)

    Returns:
      kept_runs               : list of (s, e)
      flank_positions_per_run : list of lists of qualifying flank indices per run
    """
    if flank_dist <= 0:
        return runs, [[] for _ in runs]

    n = len(residues)
    kept_runs = []
    flank_positions_per_run = []

    for s, e in runs:
        left_index = s - flank_dist
        right_index = e + flank_dist
        this_flanks = []

        if 0 <= left_index < n and residues[left_index] == flank_aa:
            this_flanks.append(left_index)
        if 0 <= right_index < n and residues[right_index] == flank_aa:
            this_flanks.append(right_index)

        if this_flanks:
            kept_runs.append((s, e))
            flank_positions_per_run.append(this_flanks)

    return kept_runs, flank_positions_per_run

# --- Stats & formatting ---

def average(vals):
    clean = [
        x for x in vals
        if x is not None and not (isinstance(x, float) and math.isnan(x))
    ]
    return (sum(clean) / len(clean)) if clean else None

def fmt_pool_table(title, data, counts):
    """
    data: dict pool_name -> {atom: avg}
    counts: dict pool_name -> count
    """
    lines = [title]
    in_run = counts.get("IN_RUN", 0)
    flank = counts.get("FLANK", 0)
    runflank = counts.get("RUN_FLANK", 0)
    lines.append(
        f"(counts: IN_RUN={in_run}, FLANK={flank}, RUN_FLANK={runflank})"
    )

    for pool in ["IN_RUN", "FLANK", "RUN_FLANK"]:
        if pool not in data:
            continue
        lines.append(f"\n[{pool}]")
        for atom in sorted(data[pool].keys()):
            v = data[pool][atom]
            s = f"{v:.3f}" if isinstance(v, (float, int)) and v is not None else "NA"
            lines.append(f"  {atom:>3}: {s}")
    return "\n".join(lines)

# --- Main CLI ---

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Averages for target residues inside poly-runs, their flanking positions, "
            "and the combined pool.\n"
            "Runs: poly-<run_aa> stretches of length >= run_length.\n\n"
            "Flanking behavior:\n"
            "  * flank-dist == 1 : no flanking filter; use all runs; FLANK empty, RUN_FLANK = IN_RUN.\n"
            "  * flank-dist >= 2 : run is kept if <flank_aa> is present at distance flank_dist on\n"
            "                      at least one side (left or right).\n"
            "Pools:\n"
            "  IN_RUN    : target_aa residues inside those runs\n"
            "  FLANK     : target_aa residues at qualifying flank positions\n"
            "  RUN_FLANK : union of IN_RUN and FLANK\n"
        )
    )
    ap.add_argument("files", nargs="+", help="Shift tables (txt/csv; wildcards allowed).")
    ap.add_argument("--target-aa", default="A", type=lambda s: s.strip().upper(),
                    help="One-letter code of amino acid to average (default A).")
    ap.add_argument("--run-aa", default=None, type=lambda s: s.strip().upper(),
                    help="One-letter code defining the poly runs (default: same as --target-aa).")
    ap.add_argument("--run-length", type=int, default=3,
                    help="Minimum length to call a run (default 3).")
    ap.add_argument("--flank-aa", default=None, type=lambda s: s.strip().upper(),
                    help="Amino acid used for flanking filter when flank-dist >= 2 (default: run-aa).")
    ap.add_argument("--flank-dist", type=int, default=1,
                    help="If 1: no flanking filter (all runs kept). "
                         "If >=2: distance from run boundaries to look for flank-aa.")
    ap.add_argument("--exact-columns", nargs="+", metavar="COL",
                    help="Use ONLY these exact header names per block (e.g., CA CB CO N H HA). (WIDE only)")
    ap.add_argument("--debug-columns", action="store_true",
                    help="Print header→atom mapping or LONG column detection per block.")
    ap.add_argument("--show-values", nargs="*", metavar="ATOM",
                    help="Print raw values for listed atoms in IN_RUN, FLANK and RUN_FLANK pools. "
                         "If provided without atoms, defaults to CA CB CO.")
    ap.add_argument("--multiple-chain", action="store_true",
                    help="Treat each chain independently (run/flank logic per chain), then pool across chains.")
    args = ap.parse_args()

    # Defaults for AA choices
    if args.run_aa is None:
        args.run_aa = args.target_aa
    if args.flank_aa is None:
        args.flank_aa = args.run_aa

    # Basic sanity checks
    for label, aa in [("target-aa", args.target_aa),
                      ("run-aa", args.run_aa),
                      ("flank-aa", args.flank_aa)]:
        if aa not in AA20:
            raise SystemExit(f"Error: --{label} must be a standard one-letter amino acid.")

    if args.flank_dist < 1:
        raise SystemExit("Error: --flank-dist must be >= 1.")

    # Expand wildcards
    filelist = []
    for f in args.files:
        filelist.extend(glob.glob(f))
    if not filelist:
        raise SystemExit("No matching files found.")

    for fname in filelist:
        print(f"\n===== Processing: {fname} =====")
        rows, atoms_present = parse_multi_header_file(
            fname,
            debug_columns=args.debug_columns,
            exact_columns=args.exact_columns,
            multiple_chain=args.multiple_chain,
        )
        if not rows:
            print("No residue rows parsed. Is the file empty or headers missing?")
            continue

        atoms_to_use = [a for a in ATOM_ORDER if a in atoms_present]

        # --- NEW: if multiple-chain, split rows by chain and process each chain independently ---
        if args.multiple_chain:
            chains = OrderedDict()
            for r in rows:
                ch = r.get("Chain", "")
                chains.setdefault(ch, []).append(r)

            # global pools across all chains
            in_run_pool_all = []
            flank_pool_all = []

            total_runs_found = 0
            total_runs_kept = 0

            for ch, chain_rows in chains.items():
                # ensure sorted by Num
                chain_rows = sorted(chain_rows, key=lambda r: r["Num"])
                residues = [r["RES"] for r in chain_rows]

                all_runs = find_runs(residues, args.run_aa, args.run_length)
                total_runs_found += len(all_runs)
                if not all_runs:
                    continue

                if args.flank_dist == 1:
                    kept_runs = all_runs
                    flank_positions_per_run = [[] for _ in all_runs]
                else:
                    kept_runs, flank_positions_per_run = flanks_any_side_by_distance(
                        residues, all_runs, args.flank_aa, args.flank_dist
                    )
                    if not kept_runs:
                        continue

                total_runs_kept += len(kept_runs)

                # pools for this chain
                in_run_pool = []
                flank_pool = []
                used_flank_indices = set()  # per chain

                for (s, e), flank_indices in zip(kept_runs, flank_positions_per_run):
                    for idx in range(s, e + 1):
                        if chain_rows[idx]["RES"] == args.target_aa:
                            in_run_pool.append(chain_rows[idx])

                    for idx in flank_indices:
                        if idx not in used_flank_indices and chain_rows[idx]["RES"] == args.target_aa:
                            flank_pool.append(chain_rows[idx])
                            used_flank_indices.add(idx)

                in_run_pool_all.extend(in_run_pool)
                flank_pool_all.extend(flank_pool)

            run_flank_pool_all = in_run_pool_all + flank_pool_all

            pools = {
                "IN_RUN": in_run_pool_all,
                "FLANK": [] if args.flank_dist == 1 else flank_pool_all,
                "RUN_FLANK": run_flank_pool_all if args.flank_dist != 1 else in_run_pool_all,
            }
            counts = {name: len(lst) for name, lst in pools.items()}

            avgs = {}
            for name, lst in pools.items():
                avgs[name] = {atom: average([r.get(atom) for r in lst]) for atom in atoms_to_use}

            if args.flank_dist == 1:
                desc = f"poly-{args.run_aa} runs (len≥{args.run_length}), no flanking filter (per chain)"
            else:
                desc = (f"poly-{args.run_aa} runs (len≥{args.run_length}) with {args.flank_aa} at distance "
                        f"{args.flank_dist} on at least one side (per chain)")

            title = f"Averages for {args.target_aa}: {desc}"
            print(fmt_pool_table(title, avgs, counts))
            print(f"\n(chain mode) runs found total={total_runs_found}, kept total={total_runs_kept}")

            # Optional raw values
            if args.show_values is not None:
                if len(args.show_values) == 0:
                    which_atoms = ["CA", "CB", "CO"]
                else:
                    which_atoms = [a.upper() for a in args.show_values]
                which_atoms = [a for a in which_atoms if a in atoms_to_use]

                for pool_name, pool in pools.items():
                    print(f"\n== Raw {pool_name} values ==")
                    for atom in which_atoms:
                        print(f"-- {atom} --")
                        print("Chain\tNum\tValue")
                        for r in pool:
                            v = r.get(atom)
                            s = "NA" if v is None else f"{v:.3f}"
                            print(f"{r.get('Chain','')}\t{r['Num']}\t{s}")
                print()

            continue  # done with this file

        # --- Original single-sequence behavior (unchanged) ---
        residues = [r["RES"] for r in rows]

        # 1) Find poly-runs of run-aa
        all_runs = find_runs(residues, args.run_aa, args.run_length)
        if not all_runs:
            print("No poly-runs found with the given run-aa and run-length.")
            continue

        # 2) Determine which runs to keep and their flanking positions
        if args.flank_dist == 1:
            kept_runs = all_runs
            flank_positions_per_run = [[] for _ in all_runs]
        else:
            kept_runs, flank_positions_per_run = flanks_any_side_by_distance(
                residues, all_runs, args.flank_aa, args.flank_dist
            )
            if not kept_runs:
                print("No runs satisfied the flanking criteria (any-side) at this distance.")
                continue

        # 3) Build pools: IN_RUN, FLANK, RUN_FLANK, restricted to target-aa
        in_run_pool = []
        flank_pool = []
        used_flank_indices = set()

        for (s, e), flank_indices in zip(kept_runs, flank_positions_per_run):
            for idx in range(s, e + 1):
                if rows[idx]["RES"] == args.target_aa:
                    in_run_pool.append(rows[idx])

            for idx in flank_indices:
                if idx not in used_flank_indices and rows[idx]["RES"] == args.target_aa:
                    flank_pool.append(rows[idx])
                    used_flank_indices.add(idx)

        run_flank_pool = in_run_pool + flank_pool

        pools = {
            "IN_RUN": in_run_pool,
            "FLANK": flank_pool,
            "RUN_FLANK": run_flank_pool,
        }
        counts = {name: len(lst) for name, lst in pools.items()}

        avgs = {}
        for name, lst in pools.items():
            avgs[name] = {
                atom: average([r.get(atom) for r in lst])
                for atom in atoms_to_use
            }

        # 4) Print summary table
        if args.flank_dist == 1:
            desc = f"poly-{args.run_aa} runs (len≥{args.run_length}), no flanking filter"
        else:
            desc = (
                f"poly-{args.run_aa} runs (len≥{args.run_length}) with {args.flank_aa} at distance "
                f"{args.flank_dist} on at least one side"
            )
        title = f"Averages for {args.target_aa}: {desc}"
        print(fmt_pool_table(title, avgs, counts))

        # 5) Optionally show raw values
        if args.show_values is not None:
            if len(args.show_values) == 0:
                which_atoms = ["CA", "CB", "CO"]
            else:
                which_atoms = [a.upper() for a in args.show_values]
            which_atoms = [a for a in which_atoms if a in atoms_to_use]

            for pool_name, pool in pools.items():
                print(f"\n== Raw {pool_name} values ==")
                for atom in which_atoms:
                    print(f"-- {atom} --")
                    print("Num\tValue")
                    for r in pool:
                        v = r.get(atom)
                        s = "NA" if v is None else f"{v:.3f}"
                        print(f"{r['Num']}\t{s}")
            print()

if __name__ == "__main__":
    main()
