# -*- coding: utf-8 -*-
# =========================================================================
# Script: S9_PySuStaln_GMV.py
# -------------------------------------------------------------------------
# Description:
#   Pipeline for SuStaIn (Subtype and Stage Inference) analysis on z-score biomarkers.
#   Process flow:
#   1. Input Excel with z-score biomarkers.
#   2. Sign flip (auto_mean or from_csv).
#   3. Feature aggregation (LH/RH totals for SalVentAttn, DorsAttn, SomMot, Vis).
#   4. SuStaIn (ZscoreSustain) with fixed z-step staging (e.g., 0.5..3.0).
#   5. Fit K=1..N_S_max_search and save per-K results.
#   6. (Optional) Cross-validation CVIC model selection -> subtype & stage assignment.
#
#   Windows-safe (spawn). No control_col.
#
# Inputs:
#   - data/Patient_Combined_Table2_GM.xlsx (or configured input)
#   - [Optional] sign_flip_flags.csv (if SIGN_FLIP_MODE="from_csv")
#
# Outputs:
#   - outputs/S9_PySuStaln/ (Results including assignments, models, and CV stats)
#
# Author: Qirui Zhang, Farber Institute for Neuroscience, Thomas Jefferson University
# Date: 12/24/2025
# =========================================================================

import os
import json
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
import multiprocessing as mp

# ---- pySuStaIn ----
from pySuStaIn.ZscoreSustain import ZscoreSustain

# =========================
# Part 1: Configuration & Paths
# =========================

# Detect Project Root
# Assumes script is in src/pipeline/ or similar depth.
# If "data" folder exists in current working directory, use CWD.
if os.path.exists("data"):
    project_root = os.getcwd()
else:
    # Fallback: climb up 3 levels from this script (src/pipeline/script.py -> ProjectRoot)
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

data_dir = os.path.join(project_root, "data")
result_dir = os.path.join(project_root, "outputs", "PySuStaln","GMV")
os.makedirs(result_dir, exist_ok=True)

# Input Settings
xlsx_filename = "Network_Wscore.xlsx"  # Updated input file
xlsx_path     = os.path.join(data_dir, xlsx_filename)
# If using absolute path for now (comment out above line if needed):
# xlsx_path = r"E:\ALL_Epi\NCT_analysis\Patient_Combined_Table2_GM.xlsx"

sheet_name = 0
id_col     = "SubID"

# Sign flip Mode: "auto_mean" or "from_csv"
SIGN_FLIP_MODE      = "auto_mean"
sign_flip_csv_path  = os.path.join(data_dir, "sign_flip_flags.csv")

# Feature Aggregation Settings
# Combine Left/Right hemisphere networks: SalVentAttn=(A+B), DorsAttn=(A+B), SomMot=(A+B), Vis=(Cent+Peri)
DROP_SOURCE_PARTS = True  # True: Replace A/B with aggregated columns; False: Keep both

COMBO_MAP = {
    # Left Hemisphere
    "GMV_LH_Striatum":    ["GMV_L_Caudate", "GMV_L_Putamen"],
    # Right Hemisphere
    "GMV_RH_Striatum":    ["GMV_R_Caudate", "GMV_R_Putamen"]
}

# Candidate Biomarker Columns (Intersection with Excel columns will be used)
candidate_cols = [
    "GMV_LH_TempPar","GMV_LH_DefaultC","GMV_LH_DefaultB","GMV_LH_DefaultA",
    "GMV_LH_Cont",
    "GMV_LH_Limbic",
    "GMV_LH_SalVentAttn",
    "GMV_LH_DorsAttn",
    "GMV_LH_SomMot",
    "GMV_LH_Vis",
    "GMV_RH_TempPar","GMV_RH_DefaultC","GMV_RH_DefaultB","GMV_RH_DefaultA",
    "GMV_RH_Cont",
    "GMV_RH_Limbic",
    "GMV_RH_SalVentAttn",
    "GMV_RH_DorsAttn",
    "GMV_RH_SomMot",
    "GMV_RH_Vis",
    "GMV_L_Caudate","GMV_L_Hippocampus",
    "GMV_L_Pallidum","GMV_L_Putamen","GMV_L_Thalamus",
    "GMV_R_Caudate","GMV_R_Hippocampus",
    "GMV_R_Pallidum","GMV_R_Putamen","GMV_R_Thalamus",
]

# SuStaIn & CV Parameters
N_S_max_search     = 4
n_startpoints      = 25
n_mcmc_iterations  = 10000              # Recommended >= 100000 for final run
random_state       = 42
dataset_name       = "dataset"

USE_PARALLEL_STARTPOINTS = True
DO_CVIC_SELECTION        = True
n_folds                  = 5

# Z-score Thresholding Parameters
Z_MIN = 1.0
Z_MAX = 3.0
STEP  = 1.0


# =========================
# Part 2: Utility Functions
# =========================

def make_thresholds_05_to_p95(series, step=0.5):
    """
    Generate threshold sequence [0.5, 1.0, ..., <=P95] for a biomarker.
    (Currently unused in main logic, retained for reference)
    """
    arr = pd.to_numeric(series, errors="coerce").to_numpy()
    arr = arr[~np.isnan(arr)]
    if arr.size == 0:
        return np.array([])
    
    max95 = float(max(0.0, np.percentile(arr, 95)))
    if max95 < step:
        return np.array([])
    return np.arange(step, max95 + 1e-9, step)


def build_sustain(X, Z_vals, Z_max, biomarker_labels, N_S_max, output_folder, dataset_label):
    """
    Instantiate ZscoreSustain object.
    """
    return ZscoreSustain(
        X, Z_vals, Z_max, biomarker_labels,
        n_startpoints, N_S_max, n_mcmc_iterations,
        output_folder, dataset_label,
        USE_PARALLEL_STARTPOINTS,
        random_state
    )


def fit_and_save_for_K(K, data_mat, z_vals, Z_max, biomarker_labels, df, id_col_name, base_outdir, dataset_lbl):
    """
    Fit SuStaIn model for a specific K (number of subtypes).
    Saves: assignment, probabilities, samples (if available), and model summary.
    """
    k_dir = os.path.join(base_outdir, f"K{K:02d}")
    os.makedirs(k_dir, exist_ok=True)

    sustain = build_sustain(
        data_mat, z_vals, Z_max, biomarker_labels, K, k_dir, f"{dataset_lbl}_K{K}"
    )

    # Run algorithm (Handle different return signatures in pySuStaIn versions)
    samples_sequence = samples_f = None
    prob_subtype_stage = None
    
    try:
        # Full return signature
        (samples_sequence, samples_f, ml_subtype, prob_ml_subtype,
         ml_stage, prob_ml_stage, prob_subtype_stage) = sustain.run_sustain_algorithm()
    except Exception:
        # Fallback 1
        sustain.run_sustain_algorithm()
        try:
            (ml_subtype, ml_stage, prob_subtype_stage) = sustain.subtype_and_stage_individuals(data_mat)
            prob_ml_subtype = None; prob_ml_stage = None
        except Exception:
            # Fallback 2 (newData)
            (ml_subtype, prob_ml_subtype, ml_stage, prob_ml_stage,
             prob_subtype, prob_stage, prob_subtype_stage) = sustain.subtype_and_stage_individuals_newData(
                    data_mat, None, None, 1000
                )

    # Save Assignments (Subtype/Stage)
    assign = pd.DataFrame({
        "subject_index": np.arange(len(df)),
        "ml_subtype": np.asarray(ml_subtype).flatten(),
        "ml_stage":   np.asarray(ml_stage).flatten()
    })
    if id_col_name and (id_col_name in df.columns):
        assign[id_col_name] = df[id_col_name].values
    assign.to_csv(os.path.join(k_dir, f"assignment_K{K}.csv"), index=False)

    # Save Probabilities (Subject x Subtype x Stage)
    if prob_subtype_stage is not None:
        prob_rows = []
        P = np.asarray(prob_subtype_stage)
        for i in range(P.shape[0]):
            for s in range(P.shape[1]):
                for st in range(P.shape[2]):
                    prob_rows.append({
                        "subject_index": i,
                        "subtype": s + 1,
                        "stage": st + 1,
                        "prob": float(P[i, s, st])
                    })
        pd.DataFrame(prob_rows).to_csv(
            os.path.join(k_dir, f"probabilities_K{K}.csv"), index=False
        )

    # Save MCMC Samples
    try:
        if samples_sequence is not None or samples_f is not None:
            np.savez_compressed(
                os.path.join(k_dir, f"sustain_samples_K{K}.npz"),
                samples_sequence=samples_sequence,
                samples_f=samples_f
            )
    except Exception:
        pass

    # Model Summary JSON
    summary = {
        "K": int(K),
        "n_subjects": int(data_mat.shape[0]),
        "n_biomarkers": int(data_mat.shape[1]),
        "n_startpoints": int(n_startpoints),
        "n_mcmc_iterations": int(n_mcmc_iterations),
        "use_parallel_startpoints": bool(USE_PARALLEL_STARTPOINTS),
        "random_state": int(random_state),
        "dataset_name": dataset_lbl,
        "outdir": k_dir
    }
    with open(os.path.join(k_dir, f"model_summary_K{K}.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    print(f"[SuStaIn] Loop K={K} completed. Results in: {k_dir}")


def _load_sign_flags_from_csv(csv_path):
    """
    Load sign flip flags from CSV. Looks for 'sign' or 'sign_after_flip' column.
    """
    sf = pd.read_csv(csv_path)
    cols = [c.lower() for c in sf.columns]
    colmap = dict(zip(cols, sf.columns))
    
    if "biomarker" not in cols:
        raise ValueError("sign_flip_flags.csv must contain 'biomarker' column.")
    
    if "sign" in cols:
        s_col = colmap["sign"]
    elif "sign_after_flip" in cols:
        s_col = colmap["sign_after_flip"]
    else:
        raise ValueError("sign_flip_flags.csv must contain 'sign' or 'sign_after_flip' column.")
    
    sub = sf[[colmap["biomarker"], s_col]].copy()
    sub.columns = ["biomarker", "sign"]
    sub["sign"] = sub["sign"].astype(float).apply(lambda v: 1.0 if v >= 0 else -1.0)
    return dict(zip(sub["biomarker"], sub["sign"]))


# =========================
# Part 3: Main Execution
# =========================

def main():
    print("=========================================")
    print("      S9 PySuStaln Correspondence        ")
    print("=========================================")
    print(f"Project Root   : {project_root}")
    print(f"Data File      : {xlsx_path}")
    print(f"Output Dir     : {result_dir}")
    print("-----------------------------------------")

    if not os.path.exists(xlsx_path):
        raise FileNotFoundError(f"Data file not found: {xlsx_path}")

    # ------------------------------------
    # 3.1 Load Excel & Feature Aggregation
    # ------------------------------------
    print("Loading data...")
    df = pd.read_excel(xlsx_path, sheet_name=sheet_name)
    
    # --- Subject Filtering ---
    # Exclude subjects with T1notUse = 1
    if "T1notUse" in df.columns:
        n_before = len(df)
        df = df[df["T1notUse"] != 1].copy()
        n_after = len(df)
        if n_before != n_after:
            print(f"[Filter] Excluded {n_before - n_after} subjects with T1notUse=1. Remaining: {n_after}")
    else:
        print("[Filter] Warning: Column 'T1notUse' not found. Skipping exclusion.")
    
    # Exclude subjects with missing IQ (no cognitive data)
    if "IQ" in df.columns:
        n_before = len(df)
        df = df[df["IQ"].notna()].copy()
        n_after = len(df)
        if n_before != n_after:
            print(f"[Filter] Excluded {n_before - n_after} subjects with missing IQ. Remaining: {n_after}")
    else:
        print("[Filter] Warning: Column 'IQ' not found. Skipping exclusion.")

    created_combo_cols = []
    for new_col, parts in COMBO_MAP.items():
        exist = [p for p in parts if p in df.columns]
        # Only aggregate if all parts exist (robustness)
        if len(exist) == len(parts):
            df[new_col] = df[exist].mean(axis=1, skipna=True)
            created_combo_cols.append(new_col)
        else:
            pass # Skip if parts missing

    if len(created_combo_cols) > 0:
        pd.Series(created_combo_cols, name="combo_created").to_csv(
            os.path.join(result_dir, "combo_created.csv"), index=False
        )

    # Filter candidate columns
    if DROP_SOURCE_PARTS:
        # Remove source components (parts) from candidates
        parts_all = set(sum(COMBO_MAP.values(), []))
        candidate_cols_final = [c for c in candidate_cols if c not in parts_all]
        # Add created aggregates
        for newc in created_combo_cols:
            if newc not in candidate_cols_final:
                candidate_cols_final.append(newc)
    else:
        # Keep source components + add aggregates
        candidate_cols_final = list(candidate_cols)
        for newc in created_combo_cols:
            if newc not in candidate_cols_final:
                candidate_cols_final.append(newc)

    # Keep only columns present in DF
    df_cols = [c for c in candidate_cols_final if c in df.columns]
    if not df_cols:
        raise ValueError("No candidate biomarker columns found in Excel. Check column names.")

    pd.Series(df_cols, name="available_biomarkers").to_csv(
        os.path.join(result_dir, "available_biomarkers.csv"), index=False
    )

    # ------------------------------------
    # 3.2 Sign Flip (No t-test)
    # ------------------------------------
    X = df[df_cols].copy()
    flip_flags = {}

    if SIGN_FLIP_MODE.lower() == "from_csv":
        if not os.path.isfile(sign_flip_csv_path):
             # Try checking if it's just 'sign_flip_flags.csv' in result_dir
             pass
        if not os.path.isfile(sign_flip_csv_path):
            raise FileNotFoundError(f"SIGN_FLIP_MODE='from_csv' but file not found: {sign_flip_csv_path}")
        
        csv_flags = _load_sign_flags_from_csv(sign_flip_csv_path)
        for col in df_cols:
            s = csv_flags.get(col, 1.0)
            if s < 0:
                X[col] = -1.0 * X[col]
                flip_flags[col] = -1
            else:
                flip_flags[col] = 1
                
        # Save used flags
        pd.DataFrame({
            "biomarker": df_cols,
            "sign_after_flip": [flip_flags[c] for c in df_cols]
        }).to_csv(os.path.join(result_dir, "sign_flip_flags_used.csv"), index=False)

    elif SIGN_FLIP_MODE.lower() == "auto_mean":
        for col in df_cols:
            mu = pd.to_numeric(df[col], errors="coerce").mean()
            if pd.notna(mu) and (mu < 0):
                X[col] = -1.0 * X[col]
                flip_flags[col] = -1
            else:
                flip_flags[col] = 1
        
        pd.DataFrame({
            "biomarker": df_cols,
            "sign_after_flip": [flip_flags[c] for c in df_cols]
        }).to_csv(os.path.join(result_dir, "sign_flip_flags_auto_mean.csv"), index=False)
    else:
        raise ValueError("SIGN_FLIP_MODE must be 'auto_mean' or 'from_csv'.")

    # ------------------------------------
    # 3.3 Adaptive Z-thresholding
    # ------------------------------------
    # Generate thresholds from Z_MIN to Z_MAX (or P95) with STEP
    if STEP <= 0: raise ValueError("STEP must be > 0.")
    if Z_MAX < Z_MIN: raise ValueError(f"Z_MAX({Z_MAX}) must be >= Z_MIN({Z_MIN}).")

    kept_cols = []
    zlists    = []   # Threshold sequence for each biomarker
    Z_max_vals = []  # Max threshold for each biomarker

    for col in df_cols:
        arr = pd.to_numeric(X[col], errors="coerce").to_numpy()
        arr = arr[~np.isnan(arr)]
        if arr.size == 0: continue
        
        p95 = float(np.percentile(arr, 95))
        p95 = max(0.0, p95)
        col_max = min(Z_MAX, p95)

        if col_max < Z_MIN:
            continue

        n_steps = int(np.floor((col_max - Z_MIN) / STEP)) + 1
        thr_i = Z_MIN + STEP * np.arange(n_steps, dtype=float)
        thr_i = np.round(thr_i, 10)

        if thr_i.size >= 1:
            kept_cols.append(col)
            zlists.append(thr_i)
            Z_max_vals.append(col_max)

    if len(kept_cols) == 0:
        raise RuntimeError("No biomarkers generated >=1 threshold. Check data scale or Z_MIN/STEP.")

    # Pad threshold sequences to matrix form (ragged array -> matrix with zeros)
    max_len = max(len(v) for v in zlists)
    z_vals  = np.zeros((len(zlists), max_len), dtype=float)
    for i, v in enumerate(zlists):
        z_vals[i, :len(v)] = v

    data_mat         = X[kept_cols].to_numpy(float)
    biomarker_labels = kept_cols
    Z_max_arr        = np.asarray(Z_max_vals, dtype=float)

    # Save Intermediate Data
    np.savetxt(os.path.join(result_dir, "z_vals.csv"), z_vals, delimiter=",")
    
    pd.DataFrame({"biomarker": biomarker_labels, "Z_max": Z_max_arr}).to_csv(
        os.path.join(result_dir, "z_max_per_biomarker.csv"), index=False
    )

    long_rows = []
    for name, thr_i in zip(biomarker_labels, zlists):
        for j, t in enumerate(thr_i, start=1):
            long_rows.append({"biomarker": name, "stage_idx": j, "threshold": float(t)})
    pd.DataFrame(long_rows).to_csv(
        os.path.join(result_dir, "z_thresholds_per_biomarker_long.csv"), index=False
    )

    pd.DataFrame(data_mat, columns=biomarker_labels).to_csv(
        os.path.join(result_dir, "data_matrix_after_flip.csv"), index=False
    )
    pd.Series(biomarker_labels, name="kept_for_thresholding").to_csv(
        os.path.join(result_dir, "kept_biomarkers.csv"), index=False
    )

    # ------------------------------------
    # 3.4 SuStaIn Modeling (Loop K)
    # ------------------------------------
    print(f"[SuStaIn] Starting model fitting for K=1..{N_S_max_search}")
    for K in range(1, N_S_max_search + 1):
        fit_and_save_for_K(
            K, data_mat, z_vals, Z_max_arr, biomarker_labels, df, id_col, result_dir, dataset_name
        )
    print(f"[SuStaIn] Fitting completed.")

    # ------------------------------------
    # 3.5 CV Model Selection (Optional)
    # ------------------------------------
    if DO_CVIC_SELECTION:
        print("[SuStaIn] Starting CVIC Cross-Validation...")
        kf = KFold(n_splits=n_folds, shuffle=True, random_state=random_state)

        test_idxs = []
        for _, te in kf.split(data_mat):
            te = np.asarray(te, dtype=int).ravel()
            test_idxs.append(te)
        
        sustain_allK = build_sustain(
            data_mat, z_vals, Z_max_arr, biomarker_labels, N_S_max_search, result_dir, dataset_name
        )
        
        # Run CVIC
        CVIC, loglike_matrix = sustain_allK.cross_validate_sustain_model(test_idxs)

        # Normalize logL matrix shape
        loglike_matrix = np.asarray(loglike_matrix)
        # Expected: (n_folds, N_S_max_search). If transposed, fix it.
        if loglike_matrix.shape == (N_S_max_search, len(test_idxs)):
            loglike_matrix = loglike_matrix.T

        if loglike_matrix.ndim != 2 or loglike_matrix.shape[1] != N_S_max_search:
            raise RuntimeError(f"loglike_matrix shape mismatch: {loglike_matrix.shape}")

        # Export Wide Format
        wide_cols = [f"K{k}" for k in range(1, N_S_max_search + 1)]
        df_logL_wide = pd.DataFrame(loglike_matrix, columns=wide_cols)
        df_logL_wide.insert(0, "fold_idx", np.arange(1, loglike_matrix.shape[0] + 1))
        df_logL_wide.to_csv(os.path.join(result_dir, "cv_logL_wide.csv"), index=False)

        # Export Long Format
        long_rows = []
        n_folds_eff = loglike_matrix.shape[0]
        for fi in range(n_folds_eff):
            for k in range(1, N_S_max_search + 1):
                long_rows.append({
                    "fold_idx": fi + 1,
                    "K": k,
                    "logL": float(loglike_matrix[fi, k - 1])
                })
        pd.DataFrame(long_rows).to_csv(os.path.join(result_dir, "cv_logL_long.csv"), index=False)

        # Summary Stats
        logL_mean = np.mean(loglike_matrix, axis=0)
        logL_std  = np.std(loglike_matrix, axis=0, ddof=1)
        logL_sem  = logL_std / np.sqrt(loglike_matrix.shape[0])
        CVIC      = np.asarray(CVIC).flatten()

        cv_df = pd.DataFrame({
            "K": np.arange(1, N_S_max_search + 1, dtype=int),
            "CVIC": CVIC,
            "logL_mean": logL_mean,
            "logL_std": logL_std,
            "logL_sem": logL_sem
        })
        cv_df.to_csv(os.path.join(result_dir, "cv_summary_tutorial_style.csv"), index=False)

        # Select Best K
        best_K = int(cv_df.loc[cv_df["CVIC"].idxmin(), "K"])
        print(f"[SuStaIn] Best K selected by CVIC: {best_K}")
    else:
        print("[SuStaIn] CVIC Selection skipped.")


# =========================
# Entry Point (Windows Safe)
# =========================

if __name__ == "__main__":
    mp.freeze_support()
    
    # Environment configs for numerical stability
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
    
    np.random.seed(random_state)
    
    main()
