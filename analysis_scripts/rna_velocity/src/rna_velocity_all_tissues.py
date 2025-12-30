#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RNA velocity example-input pipeline (publication-friendly, minimal, no Hebrew text).

This file is a cleaned re-organization of the original analysis script:
- keeps the same core ideas: neighbors in PCA space, KNN pooling, gamma estimation,
  velocity computation, projection to PCA, and ΔPC magnitude/angle outputs.
- adds a small, explicit spatial proximity layer (boundary distance threshold in microns).

Expected inputs (CSV):
  1) pc_origin_counts_0_*.csv
     Index: cell_id
     Required columns: cell_type, tissue (or group), PC1, PC2, PC3, x_um, y_um, radius_um
  2) spliced_final.csv
     Index: cell_id
     Columns: a dummy first column + gene columns (kept for compatibility)
  3) unspliced_final.csv
     Same format as spliced_final.csv
  4) neighbors.csv (optional; will be computed if missing)
     Each column is a cell_id; each column contains an ordered list of cell_ids by increasing PCA distance.

Outputs (under artifacts/):
  - gammas_table.csv
  - per_cell_delta_pc.csv
  - per_cell_mag_angle.csv
  - per_cell_min_boundary_um.csv
  - per_pair_mwu_magnitude_by_proximality.csv

Notes:
  - Boundary distance (microns): d_boundary = d_centroid - (r_cell + r_neighbor)
  - Proximal threshold: 1.0 um (default), can be changed via CLI.
"""

from __future__ import annotations

import os
import argparse
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy import stats



from scipy.stats import norm

def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction (q-values)."""
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    p_ok = p[ok]
    order = np.argsort(p_ok)
    ranked = p_ok[order]
    q_ok = ranked * (len(ranked) / (np.arange(1, len(ranked) + 1)))
    q_ok = np.minimum.accumulate(q_ok[::-1])[::-1]
    q_ok = np.clip(q_ok, 0.0, 1.0)
    q[ok] = q_ok[np.argsort(order)]
    return q

def circular_mean(angle_rad: np.ndarray) -> float:
    """Circular mean of angles in radians."""
    a = np.asarray(angle_rad, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return float("nan")
    return float(np.arctan2(np.mean(np.sin(a)), np.mean(np.cos(a))))

def circular_diff(a: float, b: float) -> float:
    """Smallest signed difference a-b on [-pi, pi]."""
    if not (np.isfinite(a) and np.isfinite(b)):
        return float("nan")
    d = a - b
    return float((d + np.pi) % (2 * np.pi) - np.pi)

# -----------------------------
# Utilities
# -----------------------------

def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path

def safe_div(a: np.ndarray, b: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    return a / np.maximum(b, eps)

def neighbors_from_pca(pc: pd.DataFrame) -> pd.DataFrame:
    names = pc.index.tolist()
    X = pc[["PC1", "PC2", "PC3"]].to_numpy(float)
    out = {}
    for i, name in enumerate(names):
        d = np.linalg.norm(X - X[i], axis=1)
        order = np.argsort(d)
        out[name] = [names[j] for j in order]
    return pd.DataFrame(out)

def project_to_fixed_pca(S0: pd.DataFrame, n_components: int = 3) -> Tuple[StandardScaler, PCA, pd.DataFrame]:
    """
    Fit scaler+PCA on S0 and return (scaler, pca, pc0_df).
    S0 should be cells x genes (float), already normalized in the way you want.
    """
    scaler = StandardScaler(with_mean=True, with_std=True).fit(S0)
    X0 = scaler.transform(S0)
    pca = PCA(n_components=n_components).fit(X0)
    pc0 = pca.transform(X0)
    pc0_df = pd.DataFrame(pc0, index=S0.index, columns=[f"PC{i+1}" for i in range(n_components)])
    return scaler, pca, pc0_df

def map_velocity_to_pca(S0: pd.DataFrame, V: pd.DataFrame, scaler: StandardScaler, pca: PCA, eps: float = 1e-2) -> pd.DataFrame:
    """
    Finite difference mapping:
      V_pc ≈ (PCA(scaler(S0 + eps*V)) - PCA(scaler(S0))) / eps
    """
    X0 = scaler.transform(S0)
    PC0 = pca.transform(X0)

    X1 = scaler.transform(S0 + eps * V)
    PC1 = pca.transform(X1)

    Vpc = (PC1 - PC0) / eps
    return pd.DataFrame(Vpc, index=S0.index, columns=[f"vPC{i+1}" for i in range(Vpc.shape[1])])


# -----------------------------
# KNN pooling + gamma estimation
# -----------------------------

@dataclass
class VelocityConfig:
    k_neighbors: int = 40
    corner_fraction: float = 0.075  # PERCENTAGE in the original script
    l_filter: float = 0.038         # expression threshold for filtering (max s_knn and max u_knn)
    proximal_threshold_um: float = 1.0
    t_future: float = 3.0

def read_counts(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    if df.shape[1] < 2:
        raise ValueError(f"Counts file must have a dummy first column + gene columns: {path}")
    # keep columns[1:] for compatibility with the original script
    return df.iloc[:, 1:].copy()

def compute_knn_pooled(S: pd.DataFrame, U: pd.DataFrame, neighbors: pd.DataFrame, k: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Pool S and U over each cell's K nearest neighbors (by PCA order in neighbors.csv).
    Pooling follows the original idea: sum neighbor counts and then normalize by total counts across the pool.
    """
    cell_ids = S.index.tolist()
    genes = S.columns.tolist()

    S_pool = np.zeros((len(cell_ids), len(genes)), dtype=float)
    U_pool = np.zeros_like(S_pool)

    # Precompute arrays for speed
    S_arr = S.to_numpy(float)
    U_arr = U.to_numpy(float)
    idx_map = {cid: i for i, cid in enumerate(cell_ids)}

    for i, cid in enumerate(cell_ids):
        nn = neighbors[cid].tolist()[:k]
        nn_idx = [idx_map[n] for n in nn if n in idx_map]
        if len(nn_idx) == 0:
            nn_idx = [i]

        S_sum = S_arr[nn_idx].sum(axis=0)
        U_sum = U_arr[nn_idx].sum(axis=0)

        S_pool[i] = safe_div(S_sum, S_sum.sum())
        U_pool[i] = safe_div(U_sum, U_sum.sum())

    return (
        pd.DataFrame(S_pool, index=cell_ids, columns=genes),
        pd.DataFrame(U_pool, index=cell_ids, columns=genes),
    )

def estimate_gammas(S_knn: pd.DataFrame, U_knn: pd.DataFrame, cfg: VelocityConfig) -> Tuple[pd.Series, pd.Series]:
    """
    Estimate gamma per gene using linear regression without intercept on points near corners:
      - top-right corner (max s_knn, max u_knn)
      - bottom-left corner (0,0)
    A gene is flagged as filtered if max(s_knn) < l_filter AND max(u_knn) < l_filter.
    """
    genes = S_knn.columns
    gammas = {}
    filtered = {}

    n = S_knn.shape[0]
    m = max(int(n * cfg.corner_fraction), 1)

    for g in genes:
        x = S_knn[g].to_numpy(float)
        y = U_knn[g].to_numpy(float)

        filtered[g] = (np.max(x) < cfg.l_filter) and (np.max(y) < cfg.l_filter)

        # Select points near the corners
        max_x, max_y = float(np.max(x)), float(np.max(y))
        # distances to top-right and bottom-left
        d_tr = np.sqrt((x - max_x) ** 2 + (y - max_y) ** 2)
        d_bl = np.sqrt((x - 0.0) ** 2 + (y - 0.0) ** 2)

        tr_idx = np.argsort(d_tr)[:m]
        bl_idx = np.argsort(d_bl)[:m]
        sel = np.unique(np.concatenate([tr_idx, bl_idx]))

        X = x[sel].reshape(-1, 1)
        Y = y[sel]

        reg = LinearRegression(fit_intercept=False)
        reg.fit(X, Y)
        gammas[g] = float(reg.coef_[0])

    return pd.Series(gammas), pd.Series(filtered)

def compute_velocity(S: pd.DataFrame, U: pd.DataFrame, gammas: pd.Series, filtered: pd.Series) -> pd.DataFrame:
    """
    v = u_norm - gamma * s_norm (gene-wise).
    Here we compute using per-cell normalized S and U (sum to 1 per cell).
    Filtered genes are zeroed.
    """
    # Per-cell normalization (same spirit as the original Cell.normalization)
    S_norm = S.div(S.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    U_norm = U.div(U.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)

    V = U_norm - S_norm.mul(gammas, axis=1)
    V.loc[:, filtered[filtered].index] = 0.0
    return V

def compute_future_state(S: pd.DataFrame, V: pd.DataFrame, t: float) -> pd.DataFrame:
    """
    s(t) = s(0) + v * t, clipped at 0.
    """
    S_norm = S.div(S.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    St = (S_norm + V * t).clip(lower=0.0)
    return St


# -----------------------------
# Spatial proximity: boundary distance
# -----------------------------

def min_boundary_distance_by_neighbor_type(meta: pd.DataFrame, threshold_um: float) -> pd.DataFrame:
    """
    For each cell, compute min boundary distance to each neighbor cell type.
    d_boundary = d_centroid - (r_cell + r_neighbor)
    """
    required = {"x_um", "y_um", "radius_um", "cell_type"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(f"Meta is missing required columns for proximity: {sorted(missing)}")

    xy = meta[["x_um", "y_um"]].to_numpy(float)
    r = meta["radius_um"].to_numpy(float)
    ctype = meta["cell_type"].astype(str).to_numpy()

    unique_types = sorted(pd.unique(ctype))
    rows = []
    for i, cell_id in enumerate(meta.index):
        row = {"cell_id": cell_id}
        for nt in unique_types:
            mask = (ctype == nt)
            mask[i] = False
            if not np.any(mask):
                row[f"min_boundary_um__to_{nt}"] = np.nan
                row[f"is_proximal__to_{nt}__thr{threshold_um:g}um"] = False
                continue

            d_centroid = np.sqrt(((xy[mask] - xy[i]) ** 2).sum(axis=1))
            d_boundary = d_centroid - (r[i] + r[mask])
            dmin = float(np.min(d_boundary))
            row[f"min_boundary_um__to_{nt}"] = dmin
            row[f"is_proximal__to_{nt}__thr{threshold_um:g}um"] = (dmin <= threshold_um)
        rows.append(row)

    return pd.DataFrame(rows).set_index("cell_id")


# -----------------------------
# Pairwise tests: magnitude by proximality
# -----------------------------

def mwu_by_proximality(meta: pd.DataFrame, per_cell: pd.DataFrame, min_boundary: pd.DataFrame, cfg: VelocityConfig) -> pd.DataFrame:
    """
    For each ordered pair (primary cell type X, neighbor type Y, X!=Y):
      - define per-cell proximal/distal groups for X based on min boundary distance to Y (<= threshold)
      - within each stratum, compare magnitude between groups (5X vs WT) using MWU

    Returns a long table.
    """
    if "group" not in meta.columns:
        raise ValueError("Meta must include a 'group' column with values like '5X' and 'WT'.")

    merged = meta.join(per_cell, how="inner").join(min_boundary, how="inner")
    if merged.empty:
        raise ValueError("No overlapping cells between meta, per_cell, and min_boundary.")

    cts = sorted(pd.unique(merged["cell_type"].astype(str)))
    rows = []

    for X in cts:
        dfX = merged[merged["cell_type"].astype(str) == X].copy()
        if dfX.empty:
            continue

        for Y in cts:
            if Y == X:
                continue

            prox_col = f"is_proximal__to_{Y}__thr{cfg.proximal_threshold_um:g}um"
            if prox_col not in dfX.columns:
                continue

            for prox_label, df_stratum in [("proximal", dfX[dfX[prox_col]]),
                                           ("distal",   dfX[~dfX[prox_col]])]:
                # Compare 5X vs WT magnitude within stratum
                a = df_stratum[df_stratum["group"] == "5X"]["mag"].to_numpy(float)
                b = df_stratum[df_stratum["group"] == "WT"]["mag"].to_numpy(float)

                if len(a) < 5 or len(b) < 5:
                    continue

                u_stat, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                A12 = u_stat / (len(a) * len(b))
                cliffs = 2 * A12 - 1

                rows.append({
                    "primary_type": X,
                    "neighbor_type": Y,
                    "stratum": prox_label,
                    "n_5X": int(len(a)),
                    "n_WT": int(len(b)),
                    "mwu_p": float(p),
                    "cliffs_delta": float(cliffs),
                    "threshold_um": float(cfg.proximal_threshold_um),
                })

    return pd.DataFrame(rows)


def permutation_test_phase(meta: pd.DataFrame, per_cell: pd.DataFrame, min_boundary: pd.DataFrame, cfg: VelocityConfig, n_perm: int = 2000, seed: int = 7) -> pd.DataFrame:
    """Distance-stratified phase analysis using a label permutation test."""
    if "group" not in meta.columns:
        raise ValueError("Meta must include a 'group' column with values like '5X' and 'WT'.")

    rng = np.random.default_rng(seed)
    merged = meta.join(per_cell, how="inner").join(min_boundary, how="inner")
    cts = sorted(pd.unique(merged["cell_type"].astype(str)))
    rows: list[dict] = []

    for X in cts:
        dfX = merged[merged["cell_type"].astype(str) == X].copy()
        if dfX.empty:
            continue

        for Y in cts:
            if Y == X:
                continue

            prox_col = f"is_proximal__to_{Y}__thr{cfg.proximal_threshold_um:g}um"
            if prox_col not in dfX.columns:
                continue

            for stratum, dfS in [("proximal", dfX[dfX[prox_col]]), ("distal", dfX[~dfX[prox_col]])]:
                ang = dfS["angle_rad"].to_numpy(float)
                grp = dfS["group"].astype(str).to_numpy()
                ok = np.isfinite(ang) & np.isin(grp, ["5X", "WT"])
                ang = ang[ok]
                grp = grp[ok]

                a = ang[grp == "5X"]
                b = ang[grp == "WT"]
                if a.size < 5 or b.size < 5:
                    continue

                diff_obs = circular_diff(circular_mean(a), circular_mean(b))

                diffs = np.zeros(n_perm, dtype=float)
                for i in range(n_perm):
                    perm = grp.copy()
                    rng.shuffle(perm)
                    diffs[i] = circular_diff(circular_mean(ang[perm == "5X"]), circular_mean(ang[perm == "WT"]))

                p = (np.sum(np.abs(diffs) >= np.abs(diff_obs)) + 1.0) / (n_perm + 1.0)

                rows.append({
                    "primary_type": X,
                    "neighbor_type": Y,
                    "stratum": stratum,
                    "n_5X": int((grp == "5X").sum()),
                    "n_WT": int((grp == "WT").sum()),
                    "delta_mean_phase_rad": float(diff_obs),
                    "p_perm": float(p),
                    "q_perm": float("nan"),
                    "n_perm": int(n_perm),
                    "threshold_um": float(cfg.proximal_threshold_um),
                })

    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_perm"] = bh_fdr(out["p_perm"].to_numpy(float))
    return out


def distance_magnitude_correlation(meta: pd.DataFrame, per_cell: pd.DataFrame, min_boundary: pd.DataFrame) -> pd.DataFrame:
    """Continuous distance vs magnitude correlation (WT and 5X) with Fisher r-to-z comparison."""
    merged = meta.join(per_cell, how="inner").join(min_boundary, how="inner")
    cts = sorted(pd.unique(merged["cell_type"].astype(str)))
    rows: list[dict] = []

    for X in cts:
        dfX = merged[merged["cell_type"].astype(str) == X].copy()
        if dfX.empty:
            continue

        for Y in cts:
            if Y == X:
                continue

            dist_col = f"min_boundary_um__to_{Y}"
            if dist_col not in dfX.columns:
                continue

            df_wt = dfX[dfX["group"] == "WT"]
            df_5x = dfX[dfX["group"] == "5X"]

            d_wt = df_wt[dist_col].to_numpy(float)
            m_wt = df_wt["mag"].to_numpy(float)
            ok_wt = np.isfinite(d_wt) & np.isfinite(m_wt)
            d_wt, m_wt = d_wt[ok_wt], m_wt[ok_wt]

            d_5x = df_5x[dist_col].to_numpy(float)
            m_5x = df_5x["mag"].to_numpy(float)
            ok_5x = np.isfinite(d_5x) & np.isfinite(m_5x)
            d_5x, m_5x = d_5x[ok_5x], m_5x[ok_5x]

            if d_wt.size < 10 or d_5x.size < 10:
                continue

            r_wt, _ = stats.pearsonr(d_wt, m_wt)
            r_5x, _ = stats.pearsonr(d_5x, m_5x)

            z_wt = np.arctanh(np.clip(r_wt, -0.999999, 0.999999))
            z_5x = np.arctanh(np.clip(r_5x, -0.999999, 0.999999))
            se = np.sqrt(1.0 / (d_wt.size - 3) + 1.0 / (d_5x.size - 3))
            z = (z_5x - z_wt) / se
            p = 2.0 * (1.0 - norm.cdf(np.abs(z)))

            rows.append({
                "primary_type": X,
                "target_type": Y,
                "n_WT": int(d_wt.size),
                "n_5X": int(d_5x.size),
                "r_WT": float(r_wt),
                "r_5X": float(r_5x),
                "z_fisher": float(z),
                "p_fisher": float(p),
                "q_fisher": float("nan"),
            })

    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_fisher"] = bh_fdr(out["p_fisher"].to_numpy(float))
    return out


def gene_velocity_distance_correlation(meta: pd.DataFrame, V: pd.DataFrame, min_boundary: pd.DataFrame) -> pd.DataFrame:
    """Gene-level correlation between velocity and continuous distance; BH within each (source,target) pair."""
    merged = meta.join(min_boundary, how="inner")
    cts = sorted(pd.unique(merged["cell_type"].astype(str)))
    genes = V.columns.tolist()
    rows: list[dict] = []

    for X in cts:
        idx_X = merged.index[merged["cell_type"].astype(str) == X].intersection(V.index)
        if idx_X.empty:
            continue

        for Y in cts:
            if Y == X:
                continue

            dist_col = f"min_boundary_um__to_{Y}"
            if dist_col not in merged.columns:
                continue

            d = merged.loc[idx_X, dist_col].to_numpy(float)
            ok_d = np.isfinite(d)
            if ok_d.sum() < 10:
                continue

            tmp = []
            pvals = []
            for g in genes:
                v = V.loc[idx_X, g].to_numpy(float)
                ok = ok_d & np.isfinite(v)
                if ok.sum() < 10:
                    continue
                r, p = stats.pearsonr(d[ok], v[ok])
                tmp.append((g, float(r), float(p), int(ok.sum())))
                pvals.append(float(p))

            if not tmp:
                continue

            qvals = bh_fdr(np.array(pvals, dtype=float))
            for i, (g, r, p, n) in enumerate(tmp):
                rows.append({
                    "source_type": X,
                    "target_type": Y,
                    "gene": g,
                    "n_cells": n,
                    "r": r,
                    "p": p,
                    "q": float(qvals[i]),
                })

    return pd.DataFrame(rows)


def mwu_by_proximality_by_region(meta: pd.DataFrame, per_cell: pd.DataFrame, min_boundary: pd.DataFrame, cfg: VelocityConfig) -> pd.DataFrame:
    """
    Region-stratified version of mwu_by_proximality.
    Repeats the distance-stratified magnitude comparison independently within each region.
    """
    if "region" not in meta.columns:
        raise ValueError("Meta must include a 'region' column for region-stratified analysis.")
    out_rows = []
    for region, metaR in meta.groupby(meta["region"].astype(str)):
        cellsR = metaR.index
        perR = per_cell.loc[per_cell.index.intersection(cellsR)]
        minR = min_boundary.loc[min_boundary.index.intersection(cellsR)]
        if perR.empty or minR.empty:
            continue
        tmp = mwu_by_proximality(metaR, perR, minR, cfg)
        if tmp.empty:
            continue
        tmp.insert(0, "region", region)
        if "q_mwu" not in tmp.columns and "mwu_p" in tmp.columns:
            tmp["q_mwu"] = bh_fdr(tmp["mwu_p"].to_numpy(float))
        out_rows.append(tmp)
    if not out_rows:
        return pd.DataFrame()
    return pd.concat(out_rows, ignore_index=True)


def permutation_test_phase_by_region(meta: pd.DataFrame, per_cell: pd.DataFrame, min_boundary: pd.DataFrame, cfg: VelocityConfig, n_perm: int = 2000, seed: int = 7) -> pd.DataFrame:
    """
    Region-stratified version of permutation_test_phase.
    Repeats the distance-stratified phase comparison independently within each region.
    """
    if "region" not in meta.columns:
        raise ValueError("Meta must include a 'region' column for region-stratified analysis.")
    out_rows = []
    for region, metaR in meta.groupby(meta["region"].astype(str)):
        cellsR = metaR.index
        perR = per_cell.loc[per_cell.index.intersection(cellsR)]
        minR = min_boundary.loc[min_boundary.index.intersection(cellsR)]
        if perR.empty or minR.empty:
            continue
        tmp = permutation_test_phase(metaR, perR, minR, cfg, n_perm=n_perm, seed=seed)
        if tmp.empty:
            continue
        tmp.insert(0, "region", region)
        out_rows.append(tmp)
    if not out_rows:
        return pd.DataFrame()
    return pd.concat(out_rows, ignore_index=True)


def gene_velocity_distance_correlation_by_region(meta: pd.DataFrame, V: pd.DataFrame, min_boundary: pd.DataFrame) -> pd.DataFrame:
    """
    Region-stratified version of gene_velocity_distance_correlation.
    Repeats gene-level velocity–distance correlations independently within each region.
    """
    if "region" not in meta.columns:
        raise ValueError("Meta must include a 'region' column for region-stratified analysis.")
    out_rows = []
    for region, metaR in meta.groupby(meta["region"].astype(str)):
        cellsR = metaR.index
        VR = V.loc[V.index.intersection(cellsR)]
        minR = min_boundary.loc[min_boundary.index.intersection(cellsR)]
        if VR.empty or minR.empty:
            continue
        tmp = gene_velocity_distance_correlation(metaR, VR, minR)
        if tmp.empty:
            continue
        tmp.insert(0, "region", region)
        out_rows.append(tmp)
    if not out_rows:
        return pd.DataFrame()
    return pd.concat(out_rows, ignore_index=True)


# -----------------------------
# Main run
# -----------------------------

def run(meta_csv: str, spliced_csv: str, unspliced_csv: str, neighbors_csv: str | None, out_dir: str, cfg: VelocityConfig) -> None:
    out_dir = ensure_dir(out_dir)

    meta = pd.read_csv(meta_csv, index_col=0)
    S = read_counts(spliced_csv)
    U = read_counts(unspliced_csv)

    # Align all indices
    shared = meta.index.intersection(S.index).intersection(U.index)
    meta = meta.loc[shared].copy()
    S = S.loc[shared].copy()
    U = U.loc[shared].copy()

    # Neighbors
    if neighbors_csv and os.path.exists(neighbors_csv):
        neighbors = pd.read_csv(neighbors_csv)
    else:
        neighbors = neighbors_from_pca(meta[["PC1", "PC2", "PC3"]])
        neighbors_path = os.path.join(out_dir, "neighbors.csv")
        neighbors.to_csv(neighbors_path, index=False)
        neighbors_csv = neighbors_path

    # KNN pooling for gamma
    S_knn, U_knn = compute_knn_pooled(S, U, neighbors, k=cfg.k_neighbors)

    gammas, filtered = estimate_gammas(S_knn, U_knn, cfg)
    gammas.name = "gamma"
    filtered.name = "filtered"

    gammas_table = pd.concat([gammas, filtered], axis=1).reset_index().rename(columns={"index": "gene"})
    gammas_table.to_csv(os.path.join(out_dir, "gammas_table.csv"), index=False)

    # Velocity and future state
    V = compute_velocity(S, U, gammas, filtered)
    St = compute_future_state(S, V, t=cfg.t_future)

    # Fixed PCA basis (fit on S0 normalized, then project St and compute ΔPC)
    S0_norm = S.div(S.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    scaler0, pca0, pc0 = project_to_fixed_pca(S0_norm, n_components=3)
    pc_t = pd.DataFrame(pca0.transform(scaler0.transform(St)), index=St.index, columns=["PC1", "PC2", "PC3"])

    delta = pc_t[["PC1", "PC2"]] - pc0[["PC1", "PC2"]]
    per_cell_delta = delta.rename(columns={"PC1": "delta_pc1", "PC2": "delta_pc2"})
    per_cell_delta.to_csv(os.path.join(out_dir, "per_cell_delta_pc.csv"))

    mag = np.hypot(per_cell_delta["delta_pc1"], per_cell_delta["delta_pc2"])
    ang = np.arctan2(per_cell_delta["delta_pc2"], per_cell_delta["delta_pc1"])
    per_cell_mag = pd.DataFrame({"mag": mag, "angle_rad": ang}, index=per_cell_delta.index)
    per_cell_mag.to_csv(os.path.join(out_dir, "per_cell_mag_angle.csv"))

    # Proximity by boundary distance
    min_boundary = min_boundary_distance_by_neighbor_type(meta, threshold_um=cfg.proximal_threshold_um)
    min_boundary.to_csv(os.path.join(out_dir, "per_cell_min_boundary_um.csv"))

    # Pairwise magnitude MWU tests by proximality (+ BH-FDR)
    per_cell = per_cell_delta.join(per_cell_mag)
    tests_mag = mwu_by_proximality(meta, per_cell, min_boundary, cfg)
    if not tests_mag.empty:
        tests_mag["q_mwu"] = bh_fdr(tests_mag["mwu_p"].to_numpy(float))
    tests_mag.to_csv(os.path.join(out_dir, "per_pair_mwu_magnitude_by_proximality.csv"), index=False)

    # Pairwise phase permutation tests by proximality (+ BH-FDR)
    tests_phase = permutation_test_phase(meta, per_cell, min_boundary, cfg, n_perm=2000, seed=7)
    tests_phase.to_csv(os.path.join(out_dir, "per_pair_phase_permutation_by_proximality.csv"), index=False)

    # Continuous distance vs magnitude correlations (WT vs 5X via Fisher r-to-z, + BH-FDR)
    tests_corr = distance_magnitude_correlation(meta, per_cell, min_boundary)
    tests_corr.to_csv(os.path.join(out_dir, "per_pair_distance_magnitude_correlation.csv"), index=False)

    # Gene-level velocity vs distance correlations (+ BH-FDR within each source-target pair)
    tests_gene = gene_velocity_distance_correlation(meta, V, min_boundary)
    tests_gene.to_csv(os.path.join(out_dir, "gene_velocity_distance_correlation.csv"), index=False)

    # Region-stratified repetition of Sections 6–7 and 9 (requires meta["region"])
    if "region" in meta.columns:
        reg_mag = mwu_by_proximality_by_region(meta, per_cell, min_boundary, cfg)
        reg_mag.to_csv(os.path.join(out_dir, "per_region_pair_mwu_magnitude_by_proximality.csv"), index=False)

        reg_phase = permutation_test_phase_by_region(meta, per_cell, min_boundary, cfg, n_perm=2000, seed=7)
        reg_phase.to_csv(os.path.join(out_dir, "per_region_pair_phase_permutation_by_proximality.csv"), index=False)

        reg_gene = gene_velocity_distance_correlation_by_region(meta, V, min_boundary)
        reg_gene.to_csv(os.path.join(out_dir, "per_region_gene_velocity_distance_correlation.csv"), index=False)

# Minimal run summary
    print("[OK] Finished RNA velocity example-input run")
    print(f"     cells={len(shared)}, genes={S.shape[1]}, K={cfg.k_neighbors}, t={cfg.t_future}")
    print(f"     RESULTS -> {os.path.abspath(out_dir)}")


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--meta_csv", required=True, help="pc_origin_counts_0_*.csv (index=cell_id)")
    p.add_argument("--spliced_csv", required=True, help="spliced_final.csv")
    p.add_argument("--unspliced_csv", required=True, help="unspliced_final.csv")
    p.add_argument("--neighbors_csv", default=None, help="neighbors.csv (optional)")
    p.add_argument("--out_dir", default="RESULTS", help="Output directory")
    p.add_argument("--k_neighbors", type=int, default=40)
    p.add_argument("--corner_fraction", type=float, default=0.075)
    p.add_argument("--l_filter", type=float, default=0.038)
    p.add_argument("--proximal_threshold_um", type=float, default=1.0)
    p.add_argument("--t_future", type=float, default=3.0)
    return p

def main() -> None:
    args = build_argparser().parse_args()
    cfg = VelocityConfig(
        k_neighbors=args.k_neighbors,
        corner_fraction=args.corner_fraction,
        l_filter=args.l_filter,
        proximal_threshold_um=args.proximal_threshold_um,
        t_future=args.t_future,
    )
    run(
        meta_csv=args.meta_csv,
        spliced_csv=args.spliced_csv,
        unspliced_csv=args.unspliced_csv,
        neighbors_csv=args.neighbors_csv,
        out_dir=args.out_dir,
        cfg=cfg,
    )

if __name__ == "__main__":
    main()
