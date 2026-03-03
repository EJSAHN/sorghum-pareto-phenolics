# -*- coding: utf-8 -*-
"""
run_reproduce_and_robustness.py
-------------------------------
One-command runner for (reproduce + robustness):

  (1) Re-run the main pipeline (sorghum_analysis_pipeline.py)
  (2) Optionally compare the new Supplementary_Data_S1_master.xlsx to a baseline
  (3) Run robustness/sensitivity add-on analyses (tools/sorghum_robustness_suite.py)

Designed workflow:
  baseline folder (data + previous outputs) -> rerun pipeline -> verify S1 -> robustness suite

Usage examples:
    python tools/run_reproduce_and_robustness.py --repo . --data /path/to/tables --baseline /path/to/baseline/outputs/tables
    python tools/run_reproduce_and_robustness.py --repo .   (if inputs already exist in repo root or repo/data)

Arguments:
    --repo       Repo root containing sorghum_analysis_pipeline.py
    --data       Folder containing Table S2.xlsx ... Table S6.xlsx (optional)
    --baseline   Baseline outputs/tables folder OR baseline Supplementary_Data_S1_master.xlsx (optional)
    --skip-pipeline  Skip re-running the pipeline and only run the robustness suite (optional)
    --topk       k for Top-k yield comparisons in the robustness suite (optional; default 8)

Notes:
- This script does NOT modify your manuscript. It only runs analysis and creates tables.
- It copies Table S2..S6 into repo/data ONLY if they are missing in repo root or repo/data.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


INPUT_FILES = ["Table S2.xlsx", "Table S3.xlsx", "Table S4.xlsx", "Table S5.xlsx", "Table S6.xlsx"]


def find_inputs(repo: Path) -> dict[str, Path]:
    """Return mapping filename -> path if found in repo root or repo/data."""
    found = {}
    for fn in INPUT_FILES:
        p1 = repo / fn
        p2 = repo / "data" / fn
        if p1.exists():
            found[fn] = p1
        elif p2.exists():
            found[fn] = p2
    return found


def ensure_inputs(repo: Path, data_root: Path | None) -> None:
    """Ensure Table S2..S6 exist in repo (root or repo/data). Copy from data_root if missing."""
    repo = repo.resolve()
    existing = find_inputs(repo)
    missing = [fn for fn in INPUT_FILES if fn not in existing]

    if not missing:
        print("[INFO] Inputs already present in repo (root or data/).")
        return

    if data_root is None:
        raise FileNotFoundError(
            "Missing input files in repo and no --data folder provided. Missing: " + ", ".join(missing)
        )

    data_root = data_root.resolve()
    dest_dir = repo / "data"
    dest_dir.mkdir(parents=True, exist_ok=True)

    for fn in missing:
        src = data_root / fn
        if not src.exists():
            raise FileNotFoundError(f"Could not find {fn} under data folder: {data_root}")
        dst = dest_dir / fn
        print(f"[INFO] Copying {src} -> {dst}")
        shutil.copy2(src, dst)


def run_py(script: Path, cwd: Path) -> None:
    """Run a python script with the current interpreter."""
    cmd = [sys.executable, str(script)]
    print("[RUN]", " ".join(cmd))
    subprocess.run(cmd, cwd=str(cwd), check=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True, help="Path to GitHub repo containing sorghum_analysis_pipeline.py")
    ap.add_argument("--data", default=None, help="Folder containing Table S2..S6 (e.g., Sorghum_Mutant_Geometry2)")
    ap.add_argument("--baseline", default=None, help="Baseline outputs/tables dir OR baseline S1 xlsx file for comparison")
    ap.add_argument("--skip-pipeline", action="store_true", help="Skip re-running pipeline and only run add-on analyses")
    ap.add_argument("--topk", type=int, default=8, help="k for Top-k yield comparisons in add-on analyses")
    args = ap.parse_args()

    repo = Path(args.repo).resolve()
    if not repo.exists():
        raise FileNotFoundError(str(repo))

    data_root = Path(args.data).resolve() if args.data else None
    pipeline = repo / "sorghum_analysis_pipeline.py"
    if not pipeline.exists():
        raise FileNotFoundError(f"Could not find pipeline script: {pipeline}")

    # 1) Inputs
    ensure_inputs(repo, data_root)

    # 2) Run pipeline
    if not args.skip_pipeline:
        run_py(pipeline, cwd=repo)
    else:
        print("[INFO] Skipping pipeline run (per --skip-pipeline).")

    # 3) Optional: compare S1 to baseline
    if args.baseline:
        checker = Path(__file__).parent / "check_s1_reproducibility.py"
        if not checker.exists():
            raise FileNotFoundError(f"Missing helper: {checker}")
        cmd = [sys.executable, str(checker), "--old", args.baseline, "--new", str(repo / "outputs" / "tables")]
        print("[RUN]", " ".join(cmd))
        subprocess.run(cmd, check=False)  # do not raise; we still want add-on tables even if mismatch

    # 4) Add-on reviewer analyses
    addon = Path(__file__).parent / "sorghum_robustness_suite.py"
    if not addon.exists():
        raise FileNotFoundError(f"Missing helper: {addon}")
    cmd = [sys.executable, str(addon), "--root", str(repo), "--topk", str(args.topk)]
    print("[RUN]", " ".join(cmd))
    subprocess.run(cmd, check=True)

    print("✅ Done. Check outputs under:")
    print(f"   {repo / 'outputs' / 'tables'}")
    print(f"   {repo / 'outputs' / 'tables' / 'reviewer_response'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
