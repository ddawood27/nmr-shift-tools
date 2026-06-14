# nmr-shift-tools
A lightweight Python toolkit for merging DSSP secondary structure assignments with SHIFTX2 backbone chemical shift predictions, then analyzing the result by residue type, sequence motif, secondary structure state, and position. Ten command-line scripts, minimal dependencies.



## Overview

These scripts address a common workflow in protein NMR analysis: you have backbone chemical shift predictions (from SHIFTX2) and secondary structure assignments (from DSSP), and you want to ask questions like:

- What is the average CA shift for all leucines in beta-sheets?
- How often does proline appear two residues before glycine?
- What secondary structure does this sequence motif adopt?
- Where in a helix does proline tend to sit?

The scripts are lightweight, dependency-minimal, and designed to be run directly from the terminal.

---

## Recommended Workflow

```
DSSP file (.dssp)  +  SHIFTX2 output (.cs / .csv)
              │
              ▼
    dssp_shiftx_merger.py
              │
              ▼
    merged_shifts_dssp.csv   ◄──── used by most other scripts
```

Most scripts accept `merged_shifts_dssp.csv` as their primary input. A few (like `average_residue.py` and `neighbor_shifts_stats.py`) also accept raw SHIFTX2 tables directly.

---

## Scripts

### `dssp_shiftx_merger.py`
Merges a DSSP `.dssp` file with SHIFTX2 chemical shift predictions into a single CSV. Produces two outputs: a DSSP-only summary and a full merged table. Supports SHIFTX2 BACKBONE ATOMS format, whitespace-table format, and generic CSV/TSV.

```bash
python3 dssp_shiftx_merger.py model.dssp shiftx2_output.csv
```

Key flags: `--assume-chain A`, `--join-how {left,inner}`, `--merged-out FILE`, `--debug`

---

### `average_residue.py`
Computes per-atom mean and standard deviation of chemical shifts for a specified amino acid across one or more input files.

```bash
python3 average_residue.py merged_shifts_dssp.csv -r K --show-values CA CB CO N H HA
```

Key flags: `-r / --residue`, `--multiple-chain`, `--show-values [ATOMS]`

---

### `ss_averages.py`
Averages a numeric column (chemical shift, phi, psi, accessibility, etc.) for a given amino acid, split by secondary structure state. Reports three groups: all residues of that amino acid, those in the target SS, and those not in it.

```bash
python3 ss_averages.py merged_shifts_dssp.csv --target-aa Q --target-ss H --value-col CA
```

Key flags: `--target-aa`, `--target-ss`, `--value-col`, `--out-csv`, `--debug`

---

### `ss_location.py`
Computes the normalized position of an amino acid within contiguous secondary-structure segments (0.0 = N-terminal edge, 1.0 = C-terminal edge). Useful for asking, for example, whether proline preferentially appears at the ends of helices.

```bash
python3 ss_location.py --dssp dssp_summary.csv --aa P --ss H --preview 10
```

Key flags: `--aa`, `--ss`, `--preview N`

---

### `ss_snippet_fullclass.py`
Classifies the secondary structure of every occurrence of a short amino acid sequence snippet (2–15 residues). Each occurrence is labeled as a single SS type (H, E, .) if all residues agree, or "mixed ss" if they don't.

```bash
python3 ss_snippet_fullclass.py --input merged_shifts_dssp.csv --snippet GPG --include-unknown
```

Key flags: `--snippet`, `--include-unknown`

---

### `snippet_shift_avg.py`
Finds all occurrences of a sequence snippet and computes mean/SD chemical shifts at a specified position within it. Works with both SHIFTX2 `.cs` files (multi-header) and merged CSVs.

```bash
python3 snippet_shift_avg.py --input merged_shifts_dssp.csv --snippet AAGLY --position 3
```

Key flags: `--snippet`, `--position` (1-based), `--multiple-chain`, `--debug`

---

### `neighbor_shifts_stats.py`
Computes mean chemical shifts for a target amino acid when it appears in a specific sequence motif relative to a partner residue. Supports forward, reverse, and either-direction matching.

```bash
python3 neighbor_shifts_stats.py merged_shifts_dssp.csv --target G --partner G --gap 1 --order forward
```

Key flags: `-t / --target`, `-p / --partner`, `-g / --gap`, `-o / --order {forward,reverse,either}`

---

### `aa_offset_stats.py`
Computes the amino acid frequency distribution at a signed sequence offset (±n) from a target residue. Useful for characterizing what tends to appear before or after a given amino acid.

```bash
python3 aa_offset_stats.py filename.csv --aa P --n 2
python3 aa_offset_stats.py filename.csv --aa P --n 2 --both-directions
```

Key flags: `--aa`, `--n`, `--both-directions`, `--no-index-check`

---

### `beta_turn_scan.py`
Scans for beta-turn types (I, II, I', II') using phi/psi angles in 4-residue windows. Classifies each window by matching against canonical dihedral targets within a configurable angle tolerance.

```bash
python3 beta_turn_scan.py merged_shifts_dssp.csv --tol 30 --out beta_turns.csv
```

Key flags: `--tol` (angle tolerance in degrees, default 30), `--out`

---

### `context_shifts.py`
Averages shifts for a target residue, bucketed by its distance to poly-runs of a specified amino acid. Useful for quantifying context effects (e.g., how do alanine shifts differ inside vs. near vs. far from poly-Ala stretches?).

```bash
python3 context_shifts.py filename.csv --target-aa A --run-aa A --run-length 5 --flank-dist 2
```

Key flags: `--target-aa`, `--run-aa`, `--run-length`, `--flank-aa`, `--flank-dist`, `--multiple-chain`, `--show-values`

---

## Input File Formats

Most scripts accept the merged CSV produced by `dssp_shiftx_merger.py`. The expected columns are:

| Column | Description |
|---|---|
| `entry_id` | Protein/structure identifier |
| `CHAINID` | Chain identifier |
| `label_seq_id` / `RESID` | Residue sequence number |
| `RESNAME` | One-letter amino acid code |
| `secondary_structure` | DSSP secondary structure code (H, E, T, S, ., etc.) |
| `accessibility` | Solvent-accessible surface area |
| `phi`, `psi` | Backbone dihedral angles (degrees) |
| `CA`, `CB`, `CO`, `N`, `H`, `HA` | Backbone chemical shifts (ppm) |

Scripts that accept raw SHIFTX2 files (`.cs`, whitespace tables) handle repeated header blocks automatically.

---

## Dependencies

- Python 3.8+
- `pandas` (required by `dssp_shiftx_merger.py`, `ss_averages.py`, `ss_location.py`, `snippet_shift_avg.py`, `ss_snippet_fullclass.py`, `beta_turn_scan.py`)
- Standard library only for `aa_offset_stats.py`, `average_residue.py`, `neighbor_shifts_stats.py`, `context_shifts.py`

Install pandas if needed:

```bash
pip install pandas
```

---

## Secondary Structure Codes (DSSP)

| Code | Meaning |
|---|---|
| H | Alpha helix |
| E | Beta strand |
| T | Hydrogen-bonded turn |
| S | Bend |
| G | 3/10 helix |
| I | Pi helix |
| B | Isolated beta bridge |
| . | Coil / unassigned |
