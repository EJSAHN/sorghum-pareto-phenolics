# -*- coding: utf-8 -*-
"""
Sorghum Mutant Analysis Pipeline (GitHub Release)
=================================================

This script rebuilds ALL results from raw inputs and writes outputs to ./outputs/.
It is designed for reviewers and public users.

INPUTS
------
Place the following Excel files either:
  (A) in the repo root (same folder as this script), OR
  (B) in a subfolder named ./data/

Required:
  - Table S2.xlsx  (Phenotypes; includes the 6 phenolics; row order defines sample_id 1..96)
  - Table S3.xlsx  (SNP burden; includes sample + no_of_total_snp + no_of_heterozygous_snp)
  - Table S4.xlsx  (Genotypes; includes id, pos, ref, Allele, and sample columns 1..96;
                    sometimes stored as "23_1", "23_2" etc — this script canonicalizes them)
  - Table S5.xlsx  (GWAS hits; preferred: chromosome + position; fallback: parse from snp)
  - Table S6.xlsx  (GWAS hits; same)

OUTPUTS
-------
Written to:
  ./outputs/tables/
    phenotypes_master.csv
    pareto_front.csv
    differential_correlation_matrix.csv
    diff_corr_pvalues.csv
    metabolic_control_shift.csv
    hit_genotypes_wide.csv
    hit_snp_altfreq_delta.csv
    key_snps_for_convergence.csv
    genetic_convergence.csv
    convergence_stats.csv
    entropy_vs_yield_data.csv
    run_summary.json
    Supplementary_Data_S1_master.xlsx

Reproducibility Notes
---------------------
- Delta AF uses: alt allele = allele != ref (robust)
- Convergence PCA encodes alt allele as the SECOND allele in "Allele" (split('/')[1])
  to match the original study’s published behavior (≈ 1.07x convergence ratio).
- Random seed fixed at 42.

Dependencies
------------
pip install pandas numpy scipy scikit-learn networkx openpyxl
Optional (UMAP):
pip install umap-learn
"""

from __future__ import annotations

import json
import re
import shutil
import time
from pathlib import Path
from typing import Iterable, List, Optional, Tuple, Dict

import numpy as np
import pandas as pd
import networkx as nx
from scipy.spatial.distance import pdist
from scipy.stats import norm
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler


# =====================================================
# CONFIG (GitHub)
# =====================================================
REPO_ROOT = Path(__file__).resolve().parent
DATA_DIR = (REPO_ROOT / "data") if (REPO_ROOT / "data").exists() else REPO_ROOT
OUT_DIR = REPO_ROOT / "outputs"

# Backup existing outputs before running (default OFF for GitHub cleanliness)
CLEAN_ROOM_BACKUP = False

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

USE_UMAP_IF_AVAILABLE = True

PHENO_COLS = [
    "luteolinidin_diglucoside",
    "luteolin_glucoside",
    "apigeninidin_glucoside",
    "luteolinidin",
    "apigeninidin",
    "5_o_me_luteolinidin",
]

PARETO_OBJECTIVES: List[Tuple[str, str]] = [
    ("aglycone_total", "max"),
    ("no_of_total_snp", "min"),
]

DELTA_THRESHOLD = 0.25
CONVERGENCE_ALT_INDEX = 1  # allele[1] (SECOND allele) for PCA encoding (published behavior)

IUPAC_BASES = {
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
}


# =====================================================
# Utilities
# =====================================================
def now_stamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def backup_outputs(outputs_dir: Path) -> None:
    if outputs_dir.exists():
        dst = outputs_dir.parent / f"{outputs_dir.name}_backup_{now_stamp()}"
        shutil.move(str(outputs_dir), str(dst))
        print(f"ℹ️  Backed up existing outputs -> {dst.name}")

def make_snake(x: object) -> str:
    s = str(x).strip().replace("\n", " ").replace("≤", "<=").replace("≥", ">=")
    s = re.sub(r"\(.*?\)", "", s)
    s = re.sub(r"[^0-9a-zA-Z]+", "_", s).strip("_").lower()
    s = re.sub(r"_+", "_", s)
    return s

def sanitize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    cols = [make_snake(c) for c in df.columns]
    seen: Dict[str, int] = {}
    out = []
    for c in cols:
        if c in seen:
            seen[c] += 1
            out.append(f"{c}_{seen[c]}")
        else:
            seen[c] = 0
            out.append(c)
    df.columns = out
    return df

def read_excel_header_by_tokens(path: Path, tokens_any: List[str], scan_rows: int = 80) -> pd.DataFrame:
    """
    Finds a header row by scanning the first scan_rows rows for presence of ALL tokens in tokens_any
    (after snake_case conversion). If that fails, tries a weaker rule: any token match.
    """
    raw = pd.read_excel(path, header=None, engine="openpyxl")
    want = [make_snake(t) for t in tokens_any]

    header_idx = None
    for i in range(min(scan_rows, len(raw))):
        row = [make_snake(v) for v in raw.iloc[i].tolist()]
        if all(t in row for t in want):
            header_idx = i
            break

    if header_idx is None:
        for i in range(min(scan_rows, len(raw))):
            row = [make_snake(v) for v in raw.iloc[i].tolist()]
            if any(t in row for t in want):
                header_idx = i
                break

    if header_idx is None:
        raise ValueError(f"Could not detect header row in {path.name} using tokens={tokens_any}")

    header = raw.iloc[header_idx].tolist()
    df = raw.iloc[header_idx + 1 :].copy()
    df.columns = [make_snake(h) for h in header]
    df = df.dropna(how="all").reset_index(drop=True)
    return df

def shannon_entropy(values: Iterable[float]) -> float:
    arr = np.array(list(values), dtype=float)
    total = np.nansum(arr)
    if not np.isfinite(total) or total <= 0:
        return np.nan
    p = arr / total
    p = p[(p > 0) & np.isfinite(p)]
    return float(-(p * np.log(p)).sum())

def compute_manifold_coords(X: np.ndarray) -> pd.DataFrame:
    Xs = StandardScaler().fit_transform(X)
    pca = PCA(n_components=5, random_state=RANDOM_SEED)
    Z = pca.fit_transform(Xs)
    out = pd.DataFrame({f"pca{i+1}": Z[:, i] for i in range(5)})

    if USE_UMAP_IF_AVAILABLE:
        try:
            import umap  # type: ignore
            reducer = umap.UMAP(n_components=2, n_neighbors=12, min_dist=0.15, random_state=RANDOM_SEED)
            Zu = reducer.fit_transform(Xs)
            out["manifold1"], out["manifold2"] = Zu[:, 0], Zu[:, 1]
            out.attrs["method"] = "UMAP"
            return out
        except Exception:
            pass

    out["manifold1"], out["manifold2"] = out["pca1"], out["pca2"]
    out.attrs["method"] = "PCA"
    return out

def pareto_mask(df: pd.DataFrame, objectives: List[Tuple[str, str]]) -> np.ndarray:
    V = []
    for col, sense in objectives:
        v = df[col].to_numpy(dtype=float)
        if sense.lower() == "min":
            v = -v
        V.append(v)
    V = np.vstack(V).T

    n = V.shape[0]
    eff = np.ones(n, dtype=bool)
    for i in range(n):
        if not eff[i]:
            continue
        dominated = np.any(np.all(V >= V[i], axis=1) & np.any(V > V[i], axis=1))
        if dominated:
            eff[i] = False
    return eff

def pareto_knee_point(df: pd.DataFrame, objectives: List[Tuple[str, str]], mask: np.ndarray) -> pd.Series:
    if len(objectives) != 2:
        return pd.Series(dtype=float)
    sub = df.loc[mask].copy()
    if sub.empty:
        return pd.Series(dtype=float)

    cols = [o[0] for o in objectives]
    Y = sub[cols].to_numpy(float)
    for j, (_, sense) in enumerate(objectives):
        if sense.lower() == "min":
            Y[:, j] = -Y[:, j]

    Yn = (Y - np.nanmin(Y, axis=0)) / (np.nanmax(Y, axis=0) - np.nanmin(Y, axis=0) + 1e-12)
    A = Yn[np.nanargmin(Yn[:, 0])]
    B = Yn[np.nanargmax(Yn[:, 0])]
    BA = B - A
    denom = np.linalg.norm(BA) + 1e-12

    dists = []
    for p in Yn:
        cross_mag = abs(BA[0] * (p[1] - A[1]) - BA[1] * (p[0] - A[0]))
        dists.append(cross_mag / denom)

    return sub.iloc[int(np.nanargmax(dists))]

def safe_mean_pdist(X: np.ndarray) -> float:
    if X.shape[0] < 2:
        return float("nan")
    return float(np.mean(pdist(X)))


# =====================================================
# GWAS hit pairs (chr,pos)
# =====================================================
def parse_snp_fallback(s: object) -> Tuple[Optional[int], Optional[int]]:
    s = str(s)
    m = re.search(r"(\d+)[^\d]+(\d+)", s)
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2))

def build_hit_pairs(s5: pd.DataFrame, s6: pd.DataFrame) -> pd.DataFrame:
    hits = pd.concat([s5, s6], ignore_index=True).copy()

    if "chromosome" in hits.columns and "position" in hits.columns:
        hp = hits[["chromosome", "position"]].copy()
        hp["chr"] = pd.to_numeric(hp["chromosome"], errors="coerce")
        hp["pos"] = pd.to_numeric(hp["position"], errors="coerce")
        hp = hp[["chr", "pos"]].dropna()
        hp["chr"] = hp["chr"].astype(int)
        hp["pos"] = hp["pos"].astype(int)
        return hp.drop_duplicates().reset_index(drop=True)

    if "snp" not in hits.columns:
        raise KeyError("S5/S6 must contain (chromosome, position) OR a snp column.")
    parsed = hits["snp"].apply(parse_snp_fallback)
    hp = pd.DataFrame(parsed.tolist(), columns=["chr", "pos"]).dropna()
    hp["chr"] = hp["chr"].astype(int)
    hp["pos"] = hp["pos"].astype(int)
    return hp.drop_duplicates().reset_index(drop=True)


# =====================================================
# Table S4 reader (includes the missing-column fix)
# =====================================================
def normalize_s4_header_token(h: object) -> str:
    try:
        if isinstance(h, (int, np.integer)):
            return str(int(h))
        if isinstance(h, (float, np.floating)) and np.isfinite(h) and float(h).is_integer():
            return str(int(h))
    except Exception:
        pass
    s = str(h).strip()
    if re.fullmatch(r"\d+\.0", s):
        return s.split(".")[0]
    return s

def read_table_s4(path: Path) -> pd.DataFrame:
    print("   -> Reading Table S4 (this can take a few minutes)...")
    raw = pd.read_excel(path, header=None, engine="openpyxl")

    header_idx = None
    for i in range(min(120, len(raw))):
        row = [str(x).strip().lower() for x in raw.iloc[i].tolist()]
        if ("id" in row) and ("pos" in row) and ("ref" in row):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError("Could not find S4 header row (need id/pos/ref).")

    header = [normalize_s4_header_token(x) for x in raw.iloc[header_idx].tolist()]
    df = raw.iloc[header_idx + 1 :].copy()
    df.columns = [make_snake(h) for h in header]
    df = df.dropna(how="all").reset_index(drop=True)

    for req in ["id", "pos", "ref"]:
        if req not in df.columns:
            raise ValueError(f"Table S4 missing required column: {req}")

    if "allele" not in df.columns:
        for cand in ["alleles", "alt", "alts"]:
            if cand in df.columns:
                df = df.rename(columns={cand: "allele"})
                break
    if "allele" not in df.columns:
        raise ValueError("Table S4 missing required column: allele (Allele)")

    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df = df.dropna(subset=["pos"]).copy()
    df["pos"] = df["pos"].astype(int)

    def id_to_num(x: object) -> Optional[int]:
        m = re.search(r"(\d+)", str(x))
        return int(m.group(1)) if m else None

    df["_id_num"] = df["id"].apply(id_to_num)
    base = pd.to_numeric(df["_id_num"], errors="coerce").min()
    if not np.isfinite(base):
        raise ValueError("Could not derive numeric base from S4 id column.")
    df["chr"] = (df["_id_num"] - base + 1).astype(int)
    df = df.rename(columns={"id": "contig"}).drop(columns=["_id_num"])

    df.columns = [normalize_s4_header_token(c) for c in df.columns]

    # canonicalize "23_1" -> "23"
    rename = {}
    for c in df.columns:
        s = str(c).strip()

        m = re.fullmatch(r"(\d+)_\d+$", s)
        if m:
            base_num = m.group(1)
            if base_num not in df.columns:
                rename[c] = base_num

        m = re.fullmatch(r"(\d+)_\d+_depth$", s)
        if m:
            base_num = m.group(1)
            target = f"{base_num}_depth"
            if target not in df.columns:
                rename[c] = target

    if rename:
        df = df.rename(columns=rename)

    need = [str(i) for i in range(1, 97)]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise ValueError(f"S4 is missing genotype columns: {missing}")

    return df


# =====================================================
# Delta AF and Convergence encoding
# =====================================================
def pick_alt_not_ref(ref: str, allele_field: str) -> Optional[str]:
    ref_u = str(ref).strip().upper()
    alleles = [a.strip().upper() for a in str(allele_field).split("/")]
    alleles = [a for a in alleles if a not in ("", "NA", "NAN", "NONE")]
    if not alleles:
        return None
    for a in alleles:
        if a != ref_u:
            return a
    return alleles[0]

def pick_allele_by_index(allele_field: str, index: int) -> Optional[str]:
    alleles = [a.strip().upper() for a in str(allele_field).split("/")]
    alleles = [a for a in alleles if a not in ("", "NA", "NAN", "NONE")]
    if not alleles:
        return None
    if len(alleles) > index:
        return alleles[index]
    return alleles[0]

def allele_freq_group(row: pd.Series, group_cols: List[str], alt: str, ref: str) -> float:
    ref_u = str(ref).strip().upper()
    alt_u = str(alt).strip().upper()

    alt_count = 0.0
    denom = 0.0

    for c in group_cols:
        call = str(row[c]).strip().upper()
        if call in ("", ".", "N", "NA", "NAN", "NONE", "0"):
            continue

        if call in ("A", "C", "G", "T"):
            if call == ref_u:
                denom += 2
            elif call == alt_u:
                alt_count += 2
                denom += 2
        elif call in IUPAC_BASES:
            bases = IUPAC_BASES[call]
            if alt_u in bases:
                alt_count += 1
            denom += 2

    return alt_count / denom if denom > 0 else np.nan


# =====================================================
# MAIN
# =====================================================
def main() -> None:
    print("Starting Sorghum Mutant Analysis Pipeline (GitHub)...")
    print(f"Data directory: {DATA_DIR.resolve()}")
    print("Dependencies: pandas numpy scipy scikit-learn networkx openpyxl (optional: umap-learn)")

    required = ["Table S2.xlsx", "Table S3.xlsx", "Table S4.xlsx", "Table S5.xlsx", "Table S6.xlsx"]
    for fn in required:
        if not (DATA_DIR / fn).exists():
            raise FileNotFoundError(
                f"Missing required input: {DATA_DIR / fn}\n"
                "Place Table S2–S6 in the repo root OR in ./data/"
            )

    if CLEAN_ROOM_BACKUP:
        backup_outputs(OUT_DIR)

    ensure_dir(OUT_DIR)
    tab_dir = OUT_DIR / "tables"
    ensure_dir(tab_dir)

    # 1) Load tables
    print("--- 1. Loading phenotype / burden / GWAS tables ---")
    s2 = sanitize_columns(read_excel_header_by_tokens(DATA_DIR / "Table S2.xlsx", ["lines"]))
    s3 = sanitize_columns(read_excel_header_by_tokens(DATA_DIR / "Table S3.xlsx", ["sample"]))
    s5 = sanitize_columns(read_excel_header_by_tokens(DATA_DIR / "Table S5.xlsx", ["snp"]))
    s6 = sanitize_columns(read_excel_header_by_tokens(DATA_DIR / "Table S6.xlsx", ["snp"]))

    for c in PHENO_COLS:
        if c not in s2.columns:
            raise KeyError(f"Table S2 missing phenolic column: {c}")

    s2 = s2.copy()
    s2["sample_id"] = np.arange(1, len(s2) + 1)

    if "lines" not in s2.columns or "sample" not in s3.columns:
        raise KeyError("Missing required identifier columns: 'lines' in S2 and 'sample' in S3.")

    master = s2.merge(s3, left_on="lines", right_on="sample", how="left")

    if "no_of_total_snp" not in master.columns or "no_of_heterozygous_snp" not in master.columns:
        raise KeyError("Table S3 must contain no_of_total_snp and no_of_heterozygous_snp.")

    master["aglycone_total"] = master[["luteolinidin", "apigeninidin", "5_o_me_luteolinidin"]].sum(axis=1)
    master["phenolic_total_6"] = master[PHENO_COLS].sum(axis=1)
    master["entropy_phenolics"] = master[PHENO_COLS].apply(lambda r: shannon_entropy(r.values), axis=1)
    master["heterozygosity_ratio"] = master["no_of_heterozygous_snp"] / master["no_of_total_snp"].replace(0, np.nan)

    # 2) Manifold + Pareto
    print("--- 2. Manifold & Pareto Analysis ---")
    coords = compute_manifold_coords(master[PHENO_COLS].to_numpy(float))
    master = pd.concat([master, coords], axis=1)

    master["pareto"] = pareto_mask(master, PARETO_OBJECTIVES)
    knee = pareto_knee_point(master, PARETO_OBJECTIVES, master["pareto"].to_numpy(dtype=bool))
    master["is_knee"] = False
    knee_name = None
    if not knee.empty:
        master.loc[master.index == knee.name, "is_knee"] = True
        knee_name = str(master.loc[master.index == knee.name, "lines"].iloc[0])

    print(f"Pareto Lines: {int(master['pareto'].sum())}, Knee: {knee_name}")

    # 3) Network + Fisher Z
    print("--- 3. Differential network analysis (Fisher's Z-test) ---")
    pareto_grp = master[master["pareto"] == True]
    non_grp = master[master["pareto"] == False]
    corr_p = pareto_grp[PHENO_COLS].astype(float).corr()
    corr_np = non_grp[PHENO_COLS].astype(float).corr()
    diff_corr = corr_p - corr_np

    n1, n2 = len(pareto_grp), len(non_grp)
    r1 = corr_p.clip(-0.999, 0.999)
    r2 = corr_np.clip(-0.999, 0.999)
    z1 = 0.5 * np.log((1 + r1) / (1 - r1))
    z2 = 0.5 * np.log((1 + r2) / (1 - r2))
    se = np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    p_values = 2 * (1 - norm.cdf(np.abs((z1 - z2) / se)))
    p_val_df = pd.DataFrame(p_values, index=corr_p.index, columns=corr_p.columns)

    # 4) Control shift
    print("--- 4. Metabolic control shift (eigenvector centrality) ---")
    G_p = nx.from_pandas_adjacency(corr_p.abs().fillna(0))
    G_n = nx.from_pandas_adjacency(corr_np.abs().fillna(0))
    cent_p = nx.eigenvector_centrality(G_p, weight="weight", max_iter=1000)
    cent_n = nx.eigenvector_centrality(G_n, weight="weight", max_iter=1000)
    control_shift = pd.DataFrame({
        "Compound": list(cent_p.keys()),
        "Centrality_Elite": list(cent_p.values()),
        "Centrality_NonElite": [cent_n[k] for k in cent_p.keys()],
    })
    control_shift["Shift"] = control_shift["Centrality_Elite"] - control_shift["Centrality_NonElite"]
    control_shift.sort_values("Shift", ascending=False, inplace=True)

    # 5) Genetic convergence
    print("--- 5. Genetic convergence (clean-room rebuild) ---")
    hit_pairs = build_hit_pairs(s5, s6).drop_duplicates().reset_index(drop=True)

    s4 = read_table_s4(DATA_DIR / "Table S4.xlsx")
    s4_hit = s4.merge(hit_pairs, on=["chr", "pos"], how="inner").copy()

    print(f"   -> Found {len(hit_pairs)} unique hit SNPs, {len(s4_hit)} present in Table S4.")
    if s4_hit.empty:
        raise RuntimeError("No hit SNPs matched in S4.")

    sample_cols = [str(i) for i in range(1, 97)]
    elite_ids = set(master.loc[master["pareto"] == True, "sample_id"].astype(int).tolist())
    elite_cols = [str(i) for i in range(1, 97) if i in elite_ids]
    non_cols = [str(i) for i in range(1, 97) if i not in elite_ids]

    hit_genowide = s4_hit[["chr", "contig", "pos", "ref", "allele"] + sample_cols].copy()
    hit_genowide.to_csv(tab_dir / "hit_genotypes_wide.csv", index=False)

    # Delta AF (alt != ref)
    rows = []
    for _, row in hit_genowide.iterrows():
        ref = str(row["ref"]).strip().upper()
        alt = pick_alt_not_ref(ref, str(row["allele"]))
        if alt is None:
            continue
        f_elite = allele_freq_group(row, elite_cols, alt=alt, ref=ref)
        f_non = allele_freq_group(row, non_cols, alt=alt, ref=ref)
        rows.append((int(row["chr"]), int(row["pos"]), f_elite, f_non, f_elite - f_non))

    diffs = pd.DataFrame(rows, columns=["chr", "pos", "alt_freq_elite", "alt_freq_nonelite", "delta"])
    diffs.to_csv(tab_dir / "hit_snp_altfreq_delta.csv", index=False)

    top_snps = diffs.loc[diffs["delta"].abs() > DELTA_THRESHOLD, ["chr", "pos"]].drop_duplicates()
    top_snps.to_csv(tab_dir / "key_snps_for_convergence.csv", index=False)
    print(f"   -> Selected {len(top_snps)} Key SNPs (|ΔAF| > {DELTA_THRESHOLD}).")

    target = pd.merge(hit_genowide, top_snps, on=["chr", "pos"], how="inner")
    if target.empty:
        raise RuntimeError("No SNPs available for convergence PCA after filtering.")

    X_list = []
    for _, row in target.iterrows():
        ref = str(row["ref"]).strip().upper()
        alt = pick_allele_by_index(str(row["allele"]), CONVERGENCE_ALT_INDEX) or ref

        vals = []
        for c in sample_cols:
            call = str(row[c]).strip().upper()
            if call in ("", ".", "N", "NA", "NAN", "NONE", "0"):
                v = np.nan
            elif call == ref:
                v = 0
            elif call == alt:
                v = 2
            elif call in IUPAC_BASES:
                v = 1
            else:
                v = np.nan
            vals.append(v)
        X_list.append(vals)

    X = np.array(X_list, dtype=float).T
    X_imp = SimpleImputer(strategy="most_frequent").fit_transform(X)
    coords_g = PCA(n_components=2, random_state=RANDOM_SEED).fit_transform(X_imp)

    elite_mask = master["pareto"].to_numpy(dtype=bool)
    d_elite = safe_mean_pdist(coords_g[elite_mask])
    d_non = safe_mean_pdist(coords_g[~elite_mask])
    ratio = d_non / d_elite if np.isfinite(d_elite) and d_elite > 0 else np.nan
    print(f"   -> Convergence Ratio: {ratio:.2f}x")

    convergence_df = pd.DataFrame(coords_g, columns=["Genetic_PC1", "Genetic_PC2"])
    convergence_df["sample_id"] = np.arange(1, 97)
    convergence_df["lines"] = master["lines"].astype(str).tolist()
    convergence_df["pareto"] = elite_mask
    convergence_df.to_csv(tab_dir / "genetic_convergence.csv", index=False)

    convergence_stats = pd.DataFrame([{
        "dist_elite": d_elite,
        "dist_non": d_non,
        "convergence_ratio_non_over_elite": ratio,
        "n_hit_pairs": int(len(hit_pairs)),
        "n_matched_s4": int(len(hit_genowide)),
        "n_key_snps": int(len(top_snps)),
        "delta_threshold": float(DELTA_THRESHOLD),
        "convergence_alt_index": int(CONVERGENCE_ALT_INDEX),
        "random_seed": int(RANDOM_SEED),
        "delta_rule": "alt != ref",
        "pca_alt_rule": "allele[1]",
    }])
    convergence_stats.to_csv(tab_dir / "convergence_stats.csv", index=False)

    # 6) Save outputs
    print("--- 6. Saving final outputs ---")
    master.to_csv(tab_dir / "phenotypes_master.csv", index=False)
    master.loc[master["pareto"]].to_csv(tab_dir / "pareto_front.csv", index=False)
    diff_corr.to_csv(tab_dir / "differential_correlation_matrix.csv")
    p_val_df.to_csv(tab_dir / "diff_corr_pvalues.csv")
    control_shift.to_csv(tab_dir / "metabolic_control_shift.csv", index=False)

    entropy_data = master[["sample_id", "lines", "entropy_phenolics", "aglycone_total", "pareto"]].copy()
    entropy_data.to_csv(tab_dir / "entropy_vs_yield_data.csv", index=False)

    out_xlsx = tab_dir / "Supplementary_Data_S1_master.xlsx"
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        master.to_excel(writer, sheet_name="phenotypes_master", index=False)
        master.loc[master["pareto"]].to_excel(writer, sheet_name="pareto_front", index=False)

        diff_corr.to_excel(writer, sheet_name="diff_correlation")
        p_val_df.to_excel(writer, sheet_name="diff_corr_pvalues")
        control_shift.to_excel(writer, sheet_name="control_shift", index=False)
        entropy_data.to_excel(writer, sheet_name="entropy_data", index=False)

        hit_pairs.sort_values(["chr", "pos"]).to_excel(writer, sheet_name="hit_snp_pairs", index=False)
        diffs.to_excel(writer, sheet_name="hit_snp_altfreq_delta", index=False)
        top_snps.sort_values(["chr", "pos"]).to_excel(writer, sheet_name="key_snps_for_convergence", index=False)

        convergence_df.to_excel(writer, sheet_name="genetic_convergence", index=False)
        convergence_stats.to_excel(writer, sheet_name="convergence_stats", index=False)

        master[["sample_id", "lines"]].drop_duplicates().to_excel(writer, sheet_name="sampleid_to_line", index=False)

    summary = {
        "repo_root": str(REPO_ROOT),
        "data_dir": str(DATA_DIR),
        "outputs_dir": str(OUT_DIR),
        "pareto_count": int(master["pareto"].sum()),
        "knee_line": knee_name,
        "hit_pairs": int(len(hit_pairs)),
        "matched_s4_hits": int(len(hit_genowide)),
        "key_snps": int(len(top_snps)),
        "convergence_ratio": float(ratio) if np.isfinite(ratio) else None,
        "delta_rule": "alt != ref",
        "pca_alt_rule": "allele[1]",
        "delta_threshold": float(DELTA_THRESHOLD),
        "convergence_alt_index": int(CONVERGENCE_ALT_INDEX),
        "random_seed": int(RANDOM_SEED),
    }
    (tab_dir / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"✅ DONE. Outputs written to: {OUT_DIR.resolve()}")


if __name__ == "__main__":
    main()
