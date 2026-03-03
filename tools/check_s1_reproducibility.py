# -*- coding: utf-8 -*-
"""
check_s1_reproducibility.py
---------------------------
Utility to compare two Supplementary_Data_S1_master.xlsx files (or two outputs/tables directories)
and report whether they are effectively identical.

Usage examples:
    python tools/check_s1_reproducibility.py --old ./outputs/tables --new /path/to/other/outputs/tables
    python tools/check_s1_reproducibility.py --old old.xlsx --new new.xlsx

Arguments:
    --old   Path to outputs/tables directory OR a .xlsx file (Supplementary_Data_S1_master.xlsx)
    --new   Path to outputs/tables directory OR a .xlsx file (Supplementary_Data_S1_master.xlsx)

Exit codes:
    0 = match (within tolerance)
    1 = mismatch
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_SHEETS = [
    "phenotypes_master",
    "pareto_front",
    "entropy_data",
    "control_shift",
    "diff_correlation",
    "diff_corr_pvalues",
    "convergence_stats",
]


def resolve_s1(p: Path) -> Path:
    if p.is_dir():
        cand = p / "Supplementary_Data_S1_master.xlsx"
        if cand.exists():
            return cand
        # sometimes user points to outputs instead of outputs/tables
        cand2 = p / "tables" / "Supplementary_Data_S1_master.xlsx"
        if cand2.exists():
            return cand2
        raise FileNotFoundError(f"Could not find Supplementary_Data_S1_master.xlsx under: {p}")
    if p.is_file():
        return p
    raise FileNotFoundError(str(p))


def compare_frames(a: pd.DataFrame, b: pd.DataFrame, atol: float, rtol: float) -> tuple[bool, str]:
    if a.shape != b.shape:
        return False, f"Shape differs: {a.shape} vs {b.shape}"

    if list(a.columns) != list(b.columns):
        return False, "Column names/order differs."

    # Try numeric comparison where possible; otherwise fallback to string comparison.
    diffs = []
    for col in a.columns:
        sa, sb = a[col], b[col]
        # Normalize NaNs
        if pd.api.types.is_numeric_dtype(sa) and pd.api.types.is_numeric_dtype(sb):
            va = pd.to_numeric(sa, errors="coerce").to_numpy()
            vb = pd.to_numeric(sb, errors="coerce").to_numpy()
            ok = np.allclose(va, vb, atol=atol, rtol=rtol, equal_nan=True)
            if not ok:
                # show worst offender
                delta = np.nanmax(np.abs(va - vb))
                diffs.append(f"{col}: max|Δ|={delta}")
        else:
            ta = sa.astype(str).fillna("NaN").to_numpy()
            tb = sb.astype(str).fillna("NaN").to_numpy()
            ok = np.array_equal(ta, tb)
            if not ok:
                diffs.append(f"{col}: non-numeric mismatch")
    if diffs:
        return False, " ; ".join(diffs[:8]) + (" ..." if len(diffs) > 8 else "")
    return True, "OK"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True, help="Old baseline: xlsx file or outputs/tables directory")
    ap.add_argument("--new", required=True, help="New rerun: xlsx file or outputs/tables directory")
    ap.add_argument("--sheets", nargs="*", default=DEFAULT_SHEETS, help="Sheets to compare (default: key sheets)")
    ap.add_argument("--atol", type=float, default=0.0, help="Absolute tolerance for numeric cells")
    ap.add_argument("--rtol", type=float, default=0.0, help="Relative tolerance for numeric cells")
    args = ap.parse_args()

    old_path = resolve_s1(Path(args.old))
    new_path = resolve_s1(Path(args.new))

    print(f"[INFO] old: {old_path}")
    print(f"[INFO] new: {new_path}")

    ok_all = True
    for sheet in args.sheets:
        try:
            a = pd.read_excel(old_path, sheet_name=sheet)
            b = pd.read_excel(new_path, sheet_name=sheet)
        except ValueError as e:
            print(f"[WARN] sheet '{sheet}' not found in one of the files: {e}")
            ok_all = False
            continue

        ok, msg = compare_frames(a, b, atol=args.atol, rtol=args.rtol)
        status = "PASS" if ok else "FAIL"
        print(f"[{status}] {sheet}: {msg}")
        ok_all = ok_all and ok

    if ok_all:
        print("✅ Supplementary Data S1 appears consistent.")
        return 0
    print("❌ Mismatch detected. If mismatch is expected, adjust sheets or tolerance.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
