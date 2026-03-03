# -*- coding: utf-8 -*-
"""
sorghum_robustness_suite.py
---------------------------
Robustness and sensitivity analyses for the sorghum Pareto–phenolics pipeline.

What it does (from Supplementary_Data_S1_master.xlsx):
  1) Overall + within-stratum (types) correlations between entropy and yield metrics
  2) Sensitivity of Pareto set to choice of SNP proxy (total SNP, heterozy count, heterozy ratio)
  3) Jackknife stability of Pareto membership and knee selection (default proxy: total SNP)
  4) Top-k vs Pareto comparison table (quantifies overlap)
  5) Objective-switch demo table (top-5 under different objectives)

Outputs:
  ./outputs/tables/reviewer_response/
    correlations_entropy_targets.csv
    correlations_by_type.csv (if types present)
    pareto_sensitivity_sets.csv
    pareto_sensitivity_summary.csv
    jackknife_stability.csv
    topk_vs_pareto.csv
    objective_switch_top5.csv

Usage examples:
    python tools/sorghum_robustness_suite.py --root .
    python tools/sorghum_robustness_suite.py --s1 /path/to/Supplementary_Data_S1_master.xlsx

Arguments:
    --root   Project root containing outputs/tables/Supplementary_Data_S1_master.xlsx
    --s1     Direct path to Supplementary_Data_S1_master.xlsx
    --topk   k for Top-k yield comparison (optional; default 8)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PHENO_COLS = [
    "luteolinidin_diglucoside",
    "luteolin_glucoside",
    "apigeninidin_glucoside",
    "luteolinidin",
    "apigeninidin",
    "5_o_me_luteolinidin",
]


def resolve_s1(root: Path | None, s1: Path | None) -> Path:
    if s1 is not None:
        s1p = Path(s1)
        if s1p.is_file():
            return s1p
        raise FileNotFoundError(str(s1p))

    if root is None:
        raise ValueError("Provide either --root or --s1")

    rootp = Path(root)
    cand = rootp / "outputs" / "tables" / "Supplementary_Data_S1_master.xlsx"
    if cand.exists():
        return cand
    cand2 = rootp / "outputs" / "Supplementary_Data_S1_master.xlsx"
    if cand2.exists():
        return cand2
    raise FileNotFoundError(f"Could not find Supplementary_Data_S1_master.xlsx under {rootp}")


def pareto_mask_two(df: pd.DataFrame, max_col: str, min_col: str) -> np.ndarray:
    """Efficient frontier for (maximize max_col) and (minimize min_col)."""
    y = pd.to_numeric(df[max_col], errors="coerce").to_numpy(float)
    b = pd.to_numeric(df[min_col], errors="coerce").to_numpy(float)
    V = np.vstack([y, -b]).T
    n = V.shape[0]
    eff = np.ones(n, dtype=bool)
    for i in range(n):
        if not eff[i]:
            continue
        # dominated if some point is >= in both and > in at least one
        dominated = np.any(np.all(V >= V[i], axis=1) & np.any(V > V[i], axis=1))
        if dominated:
            eff[i] = False
    return eff


def knee_point_index(y: np.ndarray, b: np.ndarray, mask: np.ndarray) -> int | None:
    """Knee point among Pareto points using distance to line between extremes (in normalized space)."""
    idx = np.where(mask)[0]
    if len(idx) == 0:
        return None

    Y = np.vstack([y[idx], -b[idx]]).T
    Yn = (Y - np.nanmin(Y, axis=0)) / (np.nanmax(Y, axis=0) - np.nanmin(Y, axis=0) + 1e-12)

    A = Yn[np.nanargmin(Yn[:, 0])]
    B = Yn[np.nanargmax(Yn[:, 0])]
    BA = B - A
    denom = np.linalg.norm(BA) + 1e-12

    dists = []
    for p in Yn:
        cross_mag = abs(BA[0] * (p[1] - A[1]) - BA[1] * (p[0] - A[0]))
        dists.append(cross_mag / denom)

    knee_local = int(np.nanargmax(dists))
    return int(idx[knee_local])


def corr_row(x: pd.Series, y: pd.Series, method: str) -> dict:
    x = pd.to_numeric(x, errors="coerce")
    y = pd.to_numeric(y, errors="coerce")
    ok = x.notna() & y.notna()
    if ok.sum() < 3:
        return {"n": int(ok.sum()), "r": np.nan, "p": np.nan}

    if method == "pearson":
        r, p = stats.pearsonr(x[ok], y[ok])
    elif method == "spearman":
        r, p = stats.spearmanr(x[ok], y[ok])
    else:
        raise ValueError(method)

    return {"n": int(ok.sum()), "r": float(r), "p": float(p)}


def drop_duplicate_columns(df_in: pd.DataFrame) -> pd.DataFrame:
    """Remove duplicate column names (keeps first occurrence)."""
    return df_in.loc[:, ~df_in.columns.duplicated()].copy()


def safe_concat(dfs: list[pd.DataFrame]) -> pd.DataFrame:
    """
    Concatenate with strong protections:
      - reset index
      - drop duplicate column names within each df
      - align on the union of columns (missing filled with NaN)
    """
    cleaned = []
    all_cols: list[str] = []
    seen = set()

    # first pass: clean each and collect union columns (preserve order)
    for d in dfs:
        d2 = drop_duplicate_columns(d.reset_index(drop=True))
        cleaned.append(d2)
        for c in d2.columns:
            if c not in seen:
                all_cols.append(c)
                seen.add(c)

    # second pass: reindex to union columns
    aligned = []
    for d in cleaned:
        aligned.append(d.reindex(columns=all_cols))
    return pd.concat(aligned, ignore_index=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=None, help="Project root containing outputs/tables/Supplementary_Data_S1_master.xlsx")
    ap.add_argument("--s1", default=None, help="Direct path to Supplementary_Data_S1_master.xlsx")
    ap.add_argument("--topk", type=int, default=8, help="k for Top-k yield comparison (default 8)")
    args = ap.parse_args()

    s1_path = resolve_s1(Path(args.root) if args.root else None, Path(args.s1) if args.s1 else None)

    df = pd.read_excel(s1_path, sheet_name="phenotypes_master")
    df = drop_duplicate_columns(df)

    # Derived totals
    df["aglycone_total"] = df[["luteolinidin", "apigeninidin", "5_o_me_luteolinidin"]].astype(float).sum(axis=1)
    df["glycoside_total"] = df[["luteolinidin_diglucoside", "luteolin_glucoside", "apigeninidin_glucoside"]].astype(float).sum(axis=1)
    df["phenolic_total_6"] = df[PHENO_COLS].astype(float).sum(axis=1)

    out_root = s1_path.parent / "reviewer_response"
    out_root.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------
    # 1) Correlations
    # -------------------------------------------------
    targets = ["aglycone_total", "glycoside_total", "phenolic_total_6"]
    rows = []
    for t in targets:
        for method in ["pearson", "spearman"]:
            rr = corr_row(df["entropy_phenolics"], df[t], method)
            rows.append({"target": t, "method": method, **rr})

    pd.DataFrame(rows).sort_values(["target", "method"]).to_csv(
        out_root / "correlations_entropy_targets.csv", index=False
    )

    # within-type correlations (structure-aware)
    type_rows = []
    if "types" in df.columns:
        for tp, sub in df.groupby("types"):
            for t in targets:
                rr = corr_row(sub["entropy_phenolics"], sub[t], "pearson")
                type_rows.append({"types": tp, "target": t, **rr})

    corr_by_type = pd.DataFrame(type_rows)
    if not corr_by_type.empty:
        corr_by_type.to_csv(out_root / "correlations_by_type.csv", index=False)

    # -------------------------------------------------
    # 2) Pareto sensitivity to burden proxy
    # -------------------------------------------------
    proxies = [
        ("no_of_total_snp", "total_snp"),
        ("no_of_heterozygous_snp", "heterozy_snp"),
        ("heterozygosity_ratio", "heterozy_ratio"),
    ]

    set_rows = []
    summary_rows = []
    base_set = None

    for col, label in proxies:
        if col not in df.columns:
            continue
        mask = pareto_mask_two(df, "aglycone_total", col)
        lines = df.loc[mask, "lines"].astype(str).tolist()
        for ln in lines:
            set_rows.append({"proxy": label, "lines": ln})

        s = set(lines)
        if base_set is None:
            base_set = s

        inter = len(s & base_set)
        union = len(s | base_set)
        summary_rows.append(
            {
                "proxy": label,
                "n_pareto": len(s),
                "intersection_with_first": inter,
                "jaccard_with_first": inter / union if union else np.nan,
            }
        )

    pd.DataFrame(set_rows).to_csv(out_root / "pareto_sensitivity_sets.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(out_root / "pareto_sensitivity_summary.csv", index=False)

    # -------------------------------------------------
    # 3) Jackknife stability (default proxy: total SNP)
    # -------------------------------------------------
    if "no_of_total_snp" in df.columns:
        y = pd.to_numeric(df["aglycone_total"], errors="coerce").to_numpy(float)
        b = pd.to_numeric(df["no_of_total_snp"], errors="coerce").to_numpy(float)
        n = len(df)

        pareto_freq = np.zeros(n, dtype=int)
        knee_freq = np.zeros(n, dtype=int)

        for leave_out in range(n):
            idx = np.array([i for i in range(n) if i != leave_out], dtype=int)
            y2, b2 = y[idx], b[idx]
            mask = pareto_mask_two(pd.DataFrame({"y": y2, "b": b2}), "y", "b")
            pareto_idx = idx[mask]
            pareto_freq[pareto_idx] += 1
            knee = knee_point_index(y2, b2, mask)
            if knee is not None:
                knee_freq[idx[knee]] += 1

        jack = pd.DataFrame(
            {
                "lines": df["lines"].astype(str),
                "pareto_jackknife_freq": pareto_freq / n,
                "knee_jackknife_freq": knee_freq / n,
                "aglycone_total": y,
                "no_of_total_snp": b,
                "types": df["types"] if "types" in df.columns else None,
            }
        ).sort_values(["pareto_jackknife_freq", "knee_jackknife_freq"], ascending=False)

        jack.to_csv(out_root / "jackknife_stability.csv", index=False)

    # -------------------------------------------------
    # 4) Top-k vs Pareto
    # -------------------------------------------------
    k = int(args.topk)
    topk = df.sort_values("aglycone_total", ascending=False).head(k).copy()

    if "no_of_total_snp" in df.columns:
        pareto_default = pareto_mask_two(df, "aglycone_total", "no_of_total_snp")
        pareto_df = df.loc[pareto_default].copy()
        topk["in_pareto_default"] = topk["lines"].astype(str).isin(pareto_df["lines"].astype(str))

        keep = ["lines", "aglycone_total", "entropy_phenolics", "no_of_total_snp", "heterozygosity_ratio", "in_pareto_default"]
        if "types" in topk.columns:
            keep.insert(1, "types")
        topk_out = drop_duplicate_columns(topk[keep])
    else:
        keep = ["lines", "aglycone_total", "entropy_phenolics"]
        if "types" in topk.columns:
            keep.insert(1, "types")
        topk_out = drop_duplicate_columns(topk[keep])

    topk_out.to_csv(out_root / "topk_vs_pareto.csv", index=False)

    # -------------------------------------------------
    # 5) Objective-switch demo (addresses "depends on target trait")
    # -------------------------------------------------
    demos = []

    def add_rank(df_in: pd.DataFrame, metric: str, k2: int, label: str) -> pd.DataFrame:
        df_in = drop_duplicate_columns(df_in)
        sub = df_in.sort_values(metric, ascending=False).head(k2).copy()
        sub["objective"] = label
        sub["rank"] = np.arange(1, len(sub) + 1)

        keep_cols = ["objective", "rank", "lines"]
        if "types" in sub.columns:
            keep_cols.append("types")
        keep_cols += [metric, "entropy_phenolics"]
        for extra in ["no_of_total_snp", "heterozygosity_ratio"]:
            if extra in sub.columns:
                keep_cols.append(extra)

        # ensure all keep cols exist
        for c in keep_cols:
            if c not in sub.columns:
                sub[c] = np.nan

        out = sub[keep_cols].copy()
        out = drop_duplicate_columns(out)
        return out

    k_demo = 5
    demos.append(add_rank(df, "aglycone_total", k_demo, "maximize_aglycone_total"))
    demos.append(add_rank(df, "entropy_phenolics", k_demo, "maximize_entropy"))
    demos.append(add_rank(df, "phenolic_total_6", k_demo, "maximize_total_phenolics"))

    demo_df = safe_concat(demos)
    demo_df.to_csv(out_root / "objective_switch_top5.csv", index=False)

    print(f"✅ Wrote reviewer-response tables to: {out_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())