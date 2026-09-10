# =========================================================================
# Script: S8_CCA.R
# -------------------------------------------------------------------------
# Description:
#   Sparse CCA (PMA) Pipeline for Network-Clinical Correspondence.
#
#   1. Loads prepared network/clinical data (from S7).
#   2. Performs Sparse CCA (PMA package) between Imaging (X) and Clinical (Y) sets.
#   3. Supports:
#      - Penalty parameter selection (CCA.permute).
#      - Significance testing (MultiCCA.permute for 1st component).
#      - Variance explained & Redundancy Index calculation.
#      - Visualization of correlations and loadings.
#      - Parallel computation (future.apply).
#
# Inputs:
#   - data/Network_Wscore.xlsx (Created by S7)
#
# Outputs:
#   - outputs/CCA/ (CSV summaries, Loadings, Plots, Matrices)
#
# Author: Qirui Zhang,  Farber Institute for Neuroscience, Department of Neurology, Thomas Jefferson University
# Email: qirui.zhang@jefferson.edu;fmrizhangqr@126.com
# Date: 12/24/2025
# =========================================================================

# =========================
# Part 1: Environment & Libraries
# =========================
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(PMA)
  library(future)
  library(future.apply)
})

# =========================
# Part 2: Configuration & Paths
# =========================
# Determine project root (dynamic)
# Checks if 'data' folder exists in current WD; otherwise assumes we are in src/pipeline
if (dir.exists("data")) {
  project_root <- getwd()
} else {
  # Fallback: assuming script run from src/pipeline or similar depth
  # Adjust as needed for specific environment
  project_root <- normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
}

data_dir    <- file.path(project_root, "data")
out_dir     <- file.path(project_root, "outputs", "CCA")
xlsx_file   <- file.path(data_dir, "Network_Wscore.xlsx")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# User Parameters
sheet                 <- 1
id_col                <- "SubID"

K_components_global   <- 4      # Max canonical components
nperms_select_global  <- 200    # Permutations for penalty selection (CCA.permute)
min_var_global        <- 1e-3   # Near-zero variance threshold

# Permutation Test Settings (MultiCCA.permute)
do_permutation        <- TRUE
nperm_test_global     <- 1000   # For p-value calculation 500-2000
seed_select_global    <- 123
seed_perm_global      <- 456

# Parallel execution
workers_perm <- max(1, parallel::detectCores(logical = TRUE) - 1)
# plan(multisession, workers = workers_perm) # Uncomment if using future_lapply inside functions

# Penalty Grid (Proportion 0~1) - Must be paired for PMA::CCA.permute
penaltyxs_ratio <- seq(0.25, 0.95, length.out = 9)
penaltyzs_ratio <- penaltyxs_ratio

# =========================
# Part 3: Data Definitions & Loading
# =========================
# We use try-catch for file reading to provide clearer error if path is wrong due to root detection
if (!file.exists(xlsx_file)) {
  stop(sprintf("Input file not found: %s\nPlease check working directory or project root logic.", xlsx_file))
}

message(sprintf("Loading data from: %s", xlsx_file))
df <- read_xlsx(xlsx_file, sheet = sheet) |> as.data.frame()

if (!id_col %in% names(df)) stop(sprintf("ID column '%s' not found in input data.", id_col))

# Filter out subjects with missing Cognition
if ("IQ" %in% names(df)) {
  n_total <- nrow(df)
  df <- df %>% filter(!is.na(IQ))
  n_rem <- nrow(df)
  if (n_total - n_rem > 0) {
    message(sprintf("Filtered out %d subjects with missing IQ. Remaining: %d", n_total - n_rem, n_rem))
  }
} else {
  warning("Column 'IQ' not found. Skipping IQ filtering.")
}

# Define Column Sets
# ------------------
# X Sets (Imaging)
x1_cols <- c("Normativity", "Non normativity", names(df)[str_starts(names(df), "Correspondence_")])  # Network Correspondence
x2_cols <- c(names(df)[str_starts(names(df), "Hubness_")])                                         # Hubness

# Y Sets (Clinical) - Reordered/Grouped
y1_cols  <- c("IQ","TMTA","TMTB","CVLT TL","CVLT LDFR","BNT","Letter Fluency",
              "Semantic fluency","Vocabulary","Similarities","Matrix Reasoning",
              "Digit Span","Block Design","Coding","WCST PR","WCST CC",
              "Logical Memory 1","Logical Memory 2","ROCF Copy","Pegboard R","Pegboard L") # Cognition

y2_cols  <- c("EpilepsyType_TLE_EXE","FTBTC_Yes_No","Lateralization_Left","Lateralization_Right","Lateralization_Unclear",
              "Pathology_UHS","Pathology_BHS","Pathology_Lesion","Pathology_Normal", "ASM_count", "AgeOnset") # Epilepsy Clinicals

y3_cols  <- c("IQ","Vocabulary","Similarities")                                        # General Intelligence
y4_cols  <- c("BNT","Letter Fluency","Semantic fluency")                               # Language
y5_cols  <- c("CVLT TL","CVLT LDFR","Logical Memory 1","Logical Memory 2")             # Memory
y6_cols  <- c("Matrix Reasoning","Block Design","WCST PR","WCST CC","ROCF Copy","Digit Span") # Executive
y7_cols  <- c("TMTA","TMTB","Coding","Pegboard R","Pegboard L")                        # Motor Skills
y8_cols  <- c("EpilepsyType_TLE_EXE","FTBTC_Yes_No")                                   # Epilepsy Type
y9_cols  <- c("Lateralization_Left","Lateralization_Right","Lateralization_Unclear")   # Lateralization
y10_cols <- c("Pathology_UHS","Pathology_BHS","Pathology_Lesion","Pathology_Normal")   # Pathology
y11_cols <- "ASM_count"
y12_cols <- c("SeizureDuration")                                          # Seizure Burden AgeOnset

keep_unique <- function(v) v[!duplicated(v)]

x_image_cols    <- keep_unique(c(x1_cols, x2_cols))
y_clinical_cols <- keep_unique(c(y1_cols,y2_cols,y3_cols,y4_cols,y5_cols,y6_cols,y7_cols,y8_cols,y9_cols,y10_cols,y11_cols,y12_cols))

# Validate columns
stopifnot(all(x_image_cols %in% names(df)), all(y_clinical_cols %in% names(df)))

# Named Lists for Iteration
X_sets <- list(
  X1_Network_Topology            = x_image_cols,
  X2_Network_Correspondence = x1_cols,
  X3_Hubness         = x2_cols
)

Y_sets <- list(
  Y1_Clinical          = y_clinical_cols,
  Y2_Cognition         = y1_cols,
  Y3_EpilepsyClinicals = y2_cols,
  Y4_GeneralIntelligence= y3_cols,
  Y5_Language          = y4_cols,
  Y6_Memory            = y5_cols,
  Y7_Executive         = y6_cols,
  Y8_MotorSkills       = y7_cols,
  Y9_EpilepsyType      = y8_cols,
  Y10_Lateralization   = y9_cols,
  Y11_Pathology        = y10_cols,
  Y12_ASM              = y11_cols,
  Y13_SeizureBurden    = y12_cols
)

# =========================
# Part 4: Utility Functions
# =========================
to_numeric_matrix <- function(d, cols) {
  m <- d[, cols, drop = FALSE]
  m[] <- lapply(m, function(x){
    if (is.factor(x)) as.numeric(as.character(x)) else as.numeric(x)
  })
  as.matrix(m)
}

scale_cols <- function(M) {
  s <- apply(M, 2, sd, na.rm = TRUE)
  m <- apply(M, 2, mean, na.rm = TRUE)
  sweep(sweep(M, 2, m, "-"), 2, ifelse(s == 0 | is.na(s), 1, s), "/")
}

# ---------- Added: Comp1 Observed Correlation + Permutation Z/p for Fixed Penalty ----------
perm_z_p_for_penalty <- function(X, Y, penx_ratio, penz_ratio,
                                 nperms = nperms_select_global,
                                 seed = seed_select_global) {
  # X, Y already standardized
  set.seed(seed)
  # Comp1 fit (fixed penalty ratio; PMA::CCA requires [0,1] for "standard")
  fit1 <- PMA::CCA(x = X, z = Y, typex = "standard", typez = "standard",
                   K = 1, standardize = FALSE,
                   penaltyx = penx_ratio, penaltyz = penz_ratio)
  U <- as.vector(X %*% fit1$u[,1])
  V <- as.vector(Y %*% fit1$v[,1])
  r_obs <- suppressWarnings(cor(U, V))
  
  # Permutation (shuffle Y rows)
  r_perm <- numeric(nperms)
  for (b in seq_len(nperms)) {
    idx <- sample.int(nrow(Y))
    fit_b <- PMA::CCA(x = X, z = Y[idx, , drop = FALSE],
                      typex = "standard", typez = "standard",
                      K = 1, standardize = FALSE,
                      penaltyx = penx_ratio, penaltyz = penz_ratio)
    U_b <- as.vector(X %*% fit_b$u[,1])
    V_b <- as.vector(Y[idx, , drop = FALSE] %*% fit_b$v[,1])
    r_perm[b] <- suppressWarnings(cor(U_b, V_b))
  }
  mu  <- mean(r_perm, na.rm = TRUE)
  sdv <- stats::sd(r_perm, na.rm = TRUE)
  z   <- if (isTRUE(is.finite(sdv)) && sdv > 0) (r_obs - mu)/sdv else NA_real_
  p   <- mean(r_perm >= r_obs, na.rm = TRUE)
  
  tibble(
    penaltyx_ratio = penx_ratio,
    penaltyz_ratio = penz_ratio,
    r_obs          = r_obs,
    r_perm_mean    = mu,
    r_perm_sd      = sdv,
    z_value        = z,
    p_value        = p,
    nperm          = nperms
  )
}

# ---------- True Data Fit + Penalty Selection (CCA.permute; Output ratio in (0,1]) ----------
fit_true_and_penalties <- function(
    X, Y,
    k_components   = K_components_global,
    nperms_select  = nperms_select_global,
    min_var        = min_var_global,
    seed_select    = seed_select_global,
    penxs          = penaltyxs_ratio,        # Added: Pass in explicitly
    penzs          = penaltyzs_ratio         # Added
) {
  # Near-zero variance filter
  keepX <- which(apply(X, 2, var, na.rm = TRUE) > min_var)
  keepY <- which(apply(Y, 2, var, na.rm = TRUE) > min_var)
  X <- as.matrix(X[, keepX, drop = FALSE])
  Y <- as.matrix(Y[, keepY, drop = FALSE])
  if (ncol(X) == 0 || ncol(Y) == 0) stop("After variance filter, X or Y has zero columns.")
  
  # Complete cases
  cc <- complete.cases(cbind(X, Y))
  X <- X[cc, , drop = FALSE]
  Y <- Y[cc, , drop = FALSE]
  if (nrow(X) < 5) stop("Too few complete cases after filtering.")
  
  # Standardize
  X <- scale_cols(X); Y <- scale_cols(Y)
  
  # Component count limit
  K_eff <- max(1, min(k_components, nrow(X)-1, ncol(X), ncol(Y)))
  
  # Select penalty parameters (Use CCA.permute for two blocks; use paired penalty ratio vectors)
  set.seed(seed_select)
  per <- PMA::CCA.permute(
    x = X, z = Y, typex = "standard", typez = "standard",
    nperms = nperms_select, trace = FALSE,
    penaltyxs = penxs,    # Use explicit grid
    penaltyzs = penzs     # Use explicit grid (paired with xs)
  )
  
  # True data fit (Sparsity using best ratio from CCA.permute)
  fit <- PMA::CCA(
    x = X, z = Y, typex = "standard", typez = "standard",
    K = K_eff, standardize = FALSE,          # Already standardized
    penaltyx = per$bestpenaltyx,
    penaltyz = per$bestpenaltyz,
    v = per$v.init
  )
  
  # Observed correlation (each component)
  obs_cor <- sapply(seq_len(fit$K), function(i) {
    ux <- as.vector(X %*% fit$u[, i])
    vz <- as.vector(Y %*% fit$v[, i])
    suppressWarnings(cor(ux, vz))
  })
  
  # Added: Penalty-sweep Z/p for given grid (Comp1 only)
  sweep_tbl <- purrr::map2_dfr(
    penxs, penzs,
    ~perm_z_p_for_penalty(X, Y, .x, .y,
                          nperms = nperms_select,  # Use selection permutation count; increase for more stability if needed
                          seed   = seed_select)
  )
  
  list(
    fit       = fit,
    X         = X,
    Y         = Y,
    obs_cor   = obs_cor,
    penalties = c(penaltyx = per$bestpenaltyx, penaltyz = per$bestpenaltyz), # Ratio
    keepX_idx = keepX,
    keepY_idx = keepY,
    cc_idx    = which(cc),    # Row indices of complete cases in original matrix
    sweep_tbl = sweep_tbl     # Return penalty sweep results
  )
}

# ---------- MultiCCA.permute: Return p & z for the first component ----------
# IMPORTANT: MultiCCA.permute requires L1 upper bound for each block (range [1, sqrt(p_k)]),
# while CCA.permute returns a ratio in (0,1]. Thus, scaling and truncation are needed.
multiCCA_pz_for_first_component <- function(
    X, Y, penalties_pair, nperm = nperm_test_global, seed = seed_perm_global
) {
  # At least two columns required for identifiability
  if (ncol(X) < 2 || ncol(Y) < 2) {
    return(list(p_first = NA_real_, z_first = NA_real_))
  }
  
  # CCA.permute ratios -> MultiCCA.permute L1 upper bounds
  px_ccap <- as.numeric(penalties_pair["penaltyx"]) # (0,1]
  pz_ccap <- as.numeric(penalties_pair["penaltyz"]) # (0,1]
  Px <- px_ccap * sqrt(ncol(X))
  Pz <- pz_ccap * sqrt(ncol(Y))
  Px <- min(max(Px, 1), sqrt(ncol(X)))  # Truncate to [1, sqrt(p)]
  Pz <- min(max(Pz, 1), sqrt(ncol(Y)))
  
  pen_mat <- matrix(c(Px, Pz), nrow = 2, ncol = 1)
  
  set.seed(seed)
  per_multi <- PMA::MultiCCA.permute(
    xlist = list(X, Y),
    penalties = pen_mat,
    type = "standard",
    nperms = nperm,
    niter = 3,
    trace = FALSE,
    standardize = TRUE
  )
  
  # Compatibility extraction (mostly scalar; if vector, take the 1st component)
  p_first <- tryCatch(as.numeric(per_multi$pvals[1]),  error = function(e) NA_real_)
  z_first <- tryCatch(as.numeric(per_multi$zstat[1]),  error = function(e) NA_real_)
  list(p_first = p_first, z_first = z_first)
}

# =========================
# Part 5: Main Execution
# =========================
# ---------- Main Loop: Iterate over all (Xk, Ym) ----------
summary_rows <- list()
all_loadings <- list()
all_u1v1     <- list()

for (x_name in names(X_sets)) {
  Xmat_raw <- to_numeric_matrix(df, X_sets[[x_name]])
  
  for (y_name in names(Y_sets)) {
    Ymat_raw <- to_numeric_matrix(df, Y_sets[[y_name]])
    
    message(sprintf("Running TRUE sCCA: %s vs %s ...", x_name, y_name))
    tr <- tryCatch(
      fit_true_and_penalties(Xmat_raw, Ymat_raw),
      error = function(e) e
    )
    if (inherits(tr, "error")) {
      warning(sprintf("sCCA failed for %s vs %s: %s", x_name, y_name, tr$message))
      next
    }
    
    fit <- tr$fit; X <- tr$X; Y <- tr$Y
    obs_cor   <- tr$obs_cor
    penalties <- tr$penalties
    K_eff <- fit$K
    
    # Export penalty sweep Z/p table (for comp1)
    sweep_file <- file.path(out_dir, paste0("penalty_sweep_comp1_", x_name, "_vs_", y_name, ".csv"))
    write.csv(tr$sweep_tbl, file = sweep_file, row.names = FALSE)
    
    # ===== MultiCCA.permute Permutation p & z (Comp1 only) =====
    if (do_permutation) {
      message("  Running MultiCCA.permute for p & z (comp1 only) ...")
      pz <- tryCatch(
        multiCCA_pz_for_first_component(X, Y, penalties,
                                        nperm = nperm_test_global,
                                        seed  = seed_perm_global),
        error = function(e) list(p_first = NA_real_, z_first = NA_real_)
      )
      p_vec <- c(pz$p_first, rep(NA_real_, max(0, K_eff - 1)))
      z_vec <- c(pz$z_first, rep(NA_real_, max(0, K_eff - 1)))
    } else {
      p_vec <- rep(NA_real_, K_eff)
      z_vec <- rep(NA_real_, K_eff)
    }
    
    # ===== Variance Explained Calculation =====
    total_var_X <- sum(apply(X, 2, var))
    total_var_Y <- sum(apply(Y, 2, var))
    var_explained_X <- numeric(K_eff)
    var_explained_Y <- numeric(K_eff)
    VE_X <- numeric(K_eff)
    VE_Y <- numeric(K_eff)
    Redund_X_from_Y <- numeric(K_eff)
    Redund_Y_from_X <- numeric(K_eff)
    
    for (ci in seq_len(K_eff)) {
      U <- as.vector(X %*% fit$u[, ci])
      V <- as.vector(Y %*% fit$v[, ci])
      
      var_explained_X[ci] <- stats::var(U) / total_var_X
      var_explained_Y[ci] <- stats::var(V) / total_var_Y
      
      corr_X_U <- suppressWarnings(stats::cor(X, U)) # pX × 1
      corr_Y_V <- suppressWarnings(stats::cor(Y, V))
      corr_X_V <- suppressWarnings(stats::cor(X, V))
      corr_Y_U <- suppressWarnings(stats::cor(Y, U))
      
      VE_X[ci]             <- mean(corr_X_U^2, na.rm = TRUE)
      VE_Y[ci]             <- mean(corr_Y_V^2, na.rm = TRUE)
      Redund_X_from_Y[ci]  <- mean(corr_X_V^2, na.rm = TRUE)
      Redund_Y_from_X[ci]  <- mean(corr_Y_U^2, na.rm = TRUE)
    }
    
    # ===== Save Non-zero Loadings =====
    u_list <- lapply(seq_len(fit$K), function(i) {
      nz <- which(fit$u[, i] != 0)
      tibble(feature = colnames(X)[nz], weight = fit$u[nz, i], comp = i, side = "X")
    })
    v_list <- lapply(seq_len(fit$K), function(i) {
      nz <- which(fit$v[, i] != 0)
      tibble(feature = colnames(Y)[nz], weight = fit$v[nz, i], comp = i, side = "Y")
    })
    ld <- bind_rows(bind_rows(u_list), bind_rows(v_list)) %>%
      mutate(X_set = x_name, Y_set = y_name)
    all_loadings[[paste(x_name, y_name, sep = "_vs_")]] <- ld
    
    # ===== Summary Table (Added Variance Ratio & Redundancy Index) =====
    sm <- tibble(
      X_set = x_name, Y_set = y_name,
      component = seq_along(obs_cor),
      canonical_correlation = obs_cor,
      p_multicca = p_vec,              # Valid for comp1 only
      z_multicca = z_vec,
      var_explained_X = var_explained_X,
      var_explained_Y = var_explained_Y,
      VE_X = VE_X,
      VE_Y = VE_Y,
      Redund_X_from_Y = Redund_X_from_Y,
      Redund_Y_from_X = Redund_Y_from_X,
      n_used = nrow(X),
      penaltyx_ratio = as.numeric(penalties["penaltyx"]),
      penaltyz_ratio = as.numeric(penalties["penaltyz"])
    )
    summary_rows[[length(summary_rows) + 1]] <- sm
    
    # ===== Export CSV (Loadings / Summary) =====
    write.csv(ld, file = file.path(out_dir, paste0("loadings_", x_name, "_vs_", y_name, ".csv")), row.names = FALSE)
    write.csv(sm, file = file.path(out_dir, paste0("summary_", x_name, "_vs_", y_name, ".csv")), row.names = FALSE)
    
    # ===== Added: Export U1 / V1 for 1st Component (Scores per Subject) =====
    if (fit$K >= 1) {
      U1 <- as.vector(tr$X %*% fit$u[, 1])
      V1 <- as.vector(tr$Y %*% fit$v[, 1])
      u1v1_tbl <- tibble(
        SubID    = df[[id_col]][tr$cc_idx],  # Complete cases only
        X_set    = x_name, Y_set = y_name,
        component= 1L,
        U1 = U1, V1 = V1
      )
      write.csv(u1v1_tbl, file = file.path(out_dir, paste0("U1V1_comp1_", x_name, "_vs_", y_name, ".csv")), row.names = FALSE)
      all_u1v1[[paste(x_name, y_name, sep = "_vs_")]] <- u1v1_tbl
    }
    
    # ===== Visualization: Correlation + p-value (MultiCCA p only) =====
    sm_plot <- sm %>%
      mutate(label_str = ifelse(is.na(p_multicca), "MultiCCA p=NA",
                                sprintf("MultiCCA p=%.3f", p_multicca)))
    ymax <- max(0.05 + suppressWarnings(max(sm_plot$canonical_correlation, na.rm = TRUE)), 0.5)
    if (!is.finite(ymax)) ymax <- 0.5
    
    p_plot <- ggplot(sm_plot, aes(x = factor(component), y = canonical_correlation)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = label_str), vjust = -0.2, size = 3, na.rm = TRUE) +
      ylim(0, ymax) +
      labs(title = paste0("Canonical correlations + MultiCCA p (", x_name, " vs ", y_name, ")"),
           x = "Component", y = "Correlation") +
      theme_minimal(base_size = 12)
    ggsave(filename = file.path(out_dir, paste0("corr_p_", x_name, "_vs_", y_name, ".png")),
           plot = p_plot, width = 6, height = 4, dpi = 150)
    
    # ===== Visualization: Non-zero Loadings per Component (Top 30 features max) =====
    if (nrow(ld) > 0) {
      for (ci in seq_len(K_eff)) {
        ld_ci <- ld %>%
          dplyr::filter(comp == ci) %>% droplevels() %>%
          arrange(desc(abs(weight))) %>%
          slice_head(n = min(30, nrow(.)))
        if (nrow(ld_ci) == 0) next
        p_ci <- ggplot(ld_ci, aes(x = reorder(feature, abs(weight)), y = weight, fill = side)) +
          geom_col(show.legend = FALSE) +
          coord_flip() +
          facet_wrap(~side, scales = "free_y") +
          scale_fill_manual(values = c("X" = "#FF6666", "Y" = "#6699FF")) +
          labs(title = paste0("Top ", nrow(ld_ci), " Non-zero loadings (Comp ", ci, "): ",
                              x_name, " vs ", y_name),
               x = "Feature", y = "Weight") +
          theme_minimal(base_size = 15)
        ggsave(filename = file.path(out_dir, paste0("loadings_comp", ci, "_", x_name, "_vs_", y_name, ".png")),
               plot = p_ci, width = 8, height = 6, dpi = 150)
      }
    }
    
  } # end for y_name
} # end for x_name

# ---------- Added: Supplementary Calculation X*X and Y*Y (Upper triangle only to avoid duplication; same output style as X*Y) ----------
# X*X
x_names <- names(X_sets)
for (i in seq_along(x_names)) {
  for (j in i:length(x_names)) {
    xA <- x_names[i]; xB <- x_names[j]
    message(sprintf("Running TRUE sCCA (X*X): %s vs %s ...", xA, xB))
    Xmat_raw <- to_numeric_matrix(df, X_sets[[xA]])
    Ymat_raw <- to_numeric_matrix(df, X_sets[[xB]])
    
    tr <- tryCatch(fit_true_and_penalties(Xmat_raw, Ymat_raw), error = function(e) e)
    if (inherits(tr, "error")) {
      warning(sprintf("sCCA failed for %s vs %s: %s", xA, xB, tr$message))
      next
    }
    
    fit <- tr$fit; X <- tr$X; Y <- tr$Y
    obs_cor   <- tr$obs_cor
    penalties <- tr$penalties
    K_eff <- fit$K
    
    sweep_file <- file.path(out_dir, paste0("penalty_sweep_comp1_", xA, "_vs_", xB, ".csv"))
    write.csv(tr$sweep_tbl, file = sweep_file, row.names = FALSE)
    
    if (do_permutation) {
      pz <- tryCatch(
        multiCCA_pz_for_first_component(X, Y, penalties,
                                        nperm = nperm_test_global,
                                        seed  = seed_perm_global),
        error = function(e) list(p_first = NA_real_, z_first = NA_real_)
      )
      p_vec <- c(pz$p_first, rep(NA_real_, max(0, K_eff - 1)))
      z_vec <- c(pz$z_first, rep(NA_real_, max(0, K_eff - 1)))
    } else {
      p_vec <- rep(NA_real_, K_eff); z_vec <- rep(NA_real_, K_eff)
    }
    
    total_var_X <- sum(apply(X, 2, var))
    total_var_Y <- sum(apply(Y, 2, var))
    var_explained_X <- numeric(K_eff)
    var_explained_Y <- numeric(K_eff)
    VE_X <- numeric(K_eff); VE_Y <- numeric(K_eff)
    Redund_X_from_Y <- numeric(K_eff); Redund_Y_from_X <- numeric(K_eff)
    
    for (ci in seq_len(K_eff)) {
      U <- as.vector(X %*% fit$u[, ci]); V <- as.vector(Y %*% fit$v[, ci])
      var_explained_X[ci] <- stats::var(U) / total_var_X
      var_explained_Y[ci] <- stats::var(V) / total_var_Y
      corr_X_U <- suppressWarnings(stats::cor(X, U))
      corr_Y_V <- suppressWarnings(stats::cor(Y, V))
      corr_X_V <- suppressWarnings(stats::cor(X, V))
      corr_Y_U <- suppressWarnings(stats::cor(Y, U))
      VE_X[ci]            <- mean(corr_X_U^2, na.rm = TRUE)
      VE_Y[ci]            <- mean(corr_Y_V^2, na.rm = TRUE)
      Redund_X_from_Y[ci] <- mean(corr_X_V^2, na.rm = TRUE)
      Redund_Y_from_X[ci] <- mean(corr_Y_U^2, na.rm = TRUE)
    }
    
    u_list <- lapply(seq_len(fit$K), function(k) {
      nz <- which(fit$u[, k] != 0)
      tibble(feature = colnames(X)[nz], weight = fit$u[nz, k], comp = k, side = "X")
    })
    v_list <- lapply(seq_len(fit$K), function(k) {
      nz <- which(fit$v[, k] != 0)
      tibble(feature = colnames(Y)[nz], weight = fit$v[nz, k], comp = k, side = "Y")
    })
    ld <- bind_rows(bind_rows(u_list), bind_rows(v_list)) %>%
      mutate(X_set = xA, Y_set = xB)
    all_loadings[[paste(xA, xB, sep = "_vs_")]] <- ld
    
    sm <- tibble(
      X_set = xA, Y_set = xB,
      component = seq_along(obs_cor),
      canonical_correlation = obs_cor,
      p_multicca = p_vec, z_multicca = z_vec,
      var_explained_X = var_explained_X, var_explained_Y = var_explained_Y,
      VE_X = VE_X, VE_Y = VE_Y,
      Redund_X_from_Y = Redund_X_from_Y, Redund_Y_from_X = Redund_Y_from_X,
      n_used = nrow(X),
      penaltyx_ratio = as.numeric(penalties["penaltyx"]),
      penaltyz_ratio = as.numeric(penalties["penaltyz"])
    )
    summary_rows[[length(summary_rows) + 1]] <- sm
    
    write.csv(ld, file = file.path(out_dir, paste0("loadings_", xA, "_vs_", xB, ".csv")), row.names = FALSE)
    write.csv(sm, file = file.path(out_dir, paste0("summary_", xA, "_vs_", xB, ".csv")),  row.names = FALSE)
    
    if (fit$K >= 1) {
      U1 <- as.vector(tr$X %*% fit$u[, 1])
      V1 <- as.vector(tr$Y %*% fit$v[, 1])
      u1v1_tbl <- tibble(
        SubID = df[[id_col]][tr$cc_idx],
        X_set = xA, Y_set = xB, component = 1L,
        U1 = U1, V1 = V1
      )
      write.csv(u1v1_tbl, file = file.path(out_dir, paste0("U1V1_comp1_", xA, "_vs_", xB, ".csv")), row.names = FALSE)
      all_u1v1[[paste(xA, xB, sep = "_vs_")]] <- u1v1_tbl
    }
    
    sm_plot <- sm %>%
      mutate(label_str = ifelse(is.na(p_multicca), "MultiCCA p=NA",
                                sprintf("MultiCCA p=%.3f", p_multicca)))
    ymax <- max(0.05 + suppressWarnings(max(sm_plot$canonical_correlation, na.rm = TRUE)), 0.5)
    if (!is.finite(ymax)) ymax <- 0.5
    p_plot <- ggplot(sm_plot, aes(x = factor(component), y = canonical_correlation)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = label_str), vjust = -0.2, size = 3, na.rm = TRUE) +
      ylim(0, ymax) +
      labs(title = paste0("Canonical correlations + MultiCCA p (", xA, " vs ", xB, ")"),
           x = "Component", y = "Correlation") +
      theme_minimal(base_size = 12)
    ggsave(filename = file.path(out_dir, paste0("corr_p_", xA, "_vs_", xB, ".png")),
           plot = p_plot, width = 6, height = 4, dpi = 150)
  }
}

# Y*Y
y_names <- names(Y_sets)
for (i in seq_along(y_names)) {
  for (j in i:length(y_names)) {
    yA <- y_names[i]; yB <- y_names[j]
    message(sprintf("Running TRUE sCCA (Y*Y): %s vs %s ...", yA, yB))
    Xmat_raw <- to_numeric_matrix(df, Y_sets[[yA]])  # Left Block
    Ymat_raw <- to_numeric_matrix(df, Y_sets[[yB]])  # Right Block
    
    tr <- tryCatch(fit_true_and_penalties(Xmat_raw, Ymat_raw), error = function(e) e)
    if (inherits(tr, "error")) {
      warning(sprintf("sCCA failed for %s vs %s: %s", yA, yB, tr$message))
      next
    }
    
    fit <- tr$fit; X <- tr$X; Y <- tr$Y
    obs_cor   <- tr$obs_cor
    penalties <- tr$penalties
    K_eff <- fit$K
    
    sweep_file <- file.path(out_dir, paste0("penalty_sweep_comp1_", yA, "_vs_", yB, ".csv"))
    write.csv(tr$sweep_tbl, file = sweep_file, row.names = FALSE)
    
    if (do_permutation) {
      pz <- tryCatch(
        multiCCA_pz_for_first_component(X, Y, penalties,
                                        nperm = nperm_test_global,
                                        seed  = seed_perm_global),
        error = function(e) list(p_first = NA_real_, z_first = NA_real_)
      )
      p_vec <- c(pz$p_first, rep(NA_real_, max(0, K_eff - 1)))
      z_vec <- c(pz$z_first, rep(NA_real_, max(0, K_eff - 1)))
    } else {
      p_vec <- rep(NA_real_, K_eff); z_vec <- rep(NA_real_, K_eff)
    }
    
    total_var_X <- sum(apply(X, 2, var))
    total_var_Y <- sum(apply(Y, 2, var))
    var_explained_X <- numeric(K_eff)
    var_explained_Y <- numeric(K_eff)
    VE_X <- numeric(K_eff); VE_Y <- numeric(K_eff)
    Redund_X_from_Y <- numeric(K_eff); Redund_Y_from_X <- numeric(K_eff)
    
    for (ci in seq_len(K_eff)) {
      U <- as.vector(X %*% fit$u[, ci]); V <- as.vector(Y %*% fit$v[, ci])
      var_explained_X[ci] <- stats::var(U) / total_var_X
      var_explained_Y[ci] <- stats::var(V) / total_var_Y
      corr_X_U <- suppressWarnings(stats::cor(X, U))
      corr_Y_V <- suppressWarnings(stats::cor(Y, V))
      corr_X_V <- suppressWarnings(stats::cor(X, V))
      corr_Y_U <- suppressWarnings(stats::cor(Y, U))
      VE_X[ci]            <- mean(corr_X_U^2, na.rm = TRUE)
      VE_Y[ci]            <- mean(corr_Y_V^2, na.rm = TRUE)
      Redund_X_from_Y[ci] <- mean(corr_X_V^2, na.rm = TRUE)
      Redund_Y_from_X[ci] <- mean(corr_Y_U^2, na.rm = TRUE)
    }
    
    u_list <- lapply(seq_len(fit$K), function(k) {
      nz <- which(fit$u[, k] != 0)
      tibble(feature = colnames(X)[nz], weight = fit$u[nz, k], comp = k, side = "X")
    })
    v_list <- lapply(seq_len(fit$K), function(k) {
      nz <- which(fit$v[, k] != 0)
      tibble(feature = colnames(Y)[nz], weight = fit$v[nz, k], comp = k, side = "Y")
    })
    ld <- bind_rows(bind_rows(u_list), bind_rows(v_list)) %>%
      mutate(X_set = yA, Y_set = yB)
    all_loadings[[paste(yA, yB, sep = "_vs_")]] <- ld
    
    sm <- tibble(
      X_set = yA, Y_set = yB,
      component = seq_along(obs_cor),
      canonical_correlation = obs_cor,
      p_multicca = p_vec, z_multicca = z_vec,
      var_explained_X = var_explained_X, var_explained_Y = var_explained_Y,
      VE_X = VE_X, VE_Y = VE_Y,
      Redund_X_from_Y = Redund_X_from_Y, Redund_Y_from_X = Redund_Y_from_X,
      n_used = nrow(X),
      penaltyx_ratio = as.numeric(penalties["penaltyx"]),
      penaltyz_ratio = as.numeric(penalties["penaltyz"])
    )
    summary_rows[[length(summary_rows) + 1]] <- sm
    
    write.csv(ld, file = file.path(out_dir, paste0("loadings_", yA, "_vs_", yB, ".csv")), row.names = FALSE)
    write.csv(sm, file = file.path(out_dir, paste0("summary_", yA, "_vs_", yB, ".csv")),  row.names = FALSE)
    
    if (fit$K >= 1) {
      U1 <- as.vector(tr$X %*% fit$u[, 1])
      V1 <- as.vector(tr$Y %*% fit$v[, 1])
      u1v1_tbl <- tibble(
        SubID = df[[id_col]][tr$cc_idx],
        X_set = yA, Y_set = yB, component = 1L,
        U1 = U1, V1 = V1
      )
      write.csv(u1v1_tbl, file = file.path(out_dir, paste0("U1V1_comp1_", yA, "_vs_", yB, ".csv")), row.names = FALSE)
      all_u1v1[[paste(yA, yB, sep = "_vs_")]] <- u1v1_tbl
    }
    
    sm_plot <- sm %>%
      mutate(label_str = ifelse(is.na(p_multicca), "MultiCCA p=NA",
                                sprintf("MultiCCA p=%.3f", p_multicca)))
    ymax <- max(0.05 + suppressWarnings(max(sm_plot$canonical_correlation, na.rm = TRUE)), 0.5)
    if (!is.finite(ymax)) ymax <- 0.5
    p_plot <- ggplot(sm_plot, aes(x = factor(component), y = canonical_correlation)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = label_str), vjust = -0.2, size = 3, na.rm = TRUE) +
      ylim(0, ymax) +
      labs(title = paste0("Canonical correlations + MultiCCA p (", yA, " vs ", yB, ")"),
           x = "Component", y = "Correlation") +
      theme_minimal(base_size = 12)
    ggsave(filename = file.path(out_dir, paste0("corr_p_", yA, "_vs_", yB, ".png")),
           plot = p_plot, width = 6, height = 4, dpi = 150)
  }
}

# ===== Summary Export (Keep original export unchanged) =====
summary_tbl <- bind_rows(summary_rows)
write.csv(summary_tbl,   file = file.path(out_dir, "ALL_pairs_summary.csv"), row.names = FALSE)

all_loadings_tbl <- bind_rows(all_loadings)
write.csv(all_loadings_tbl, file = file.path(out_dir, "ALL_pairs_loadings.csv"), row.names = FALSE)

# Added: Summary of all X*Y U1/V1 (Keep original logic)
if (length(all_u1v1) > 0) {
  U1V1_all_tbl <- dplyr::bind_rows(all_u1v1)
  write.csv(U1V1_all_tbl, file = file.path(out_dir, "ALL_pairs_U1V1_comp1.csv"), row.names = FALSE)
}
message("Done. Results saved to: ", normalizePath(out_dir))

# ===== Compile R-value and P-value Matrix for 1st Component (Keep original X*Y version) =====
X_names <- names(X_sets)
Y_names <- names(Y_sets)

# Initialize empty matrix
R_mat <- matrix(NA_real_, nrow = length(X_names), ncol = length(Y_names),
                dimnames = list(X_names, Y_names))
P_mat <- matrix(NA_real_, nrow = length(X_names), ncol = length(Y_names),
                dimnames = list(X_names, Y_names))

# Extract data from summary_rows to fill (X*Y only)
for (sm in summary_rows) {
  if (!"component" %in% names(sm)) next
  sm1 <- sm %>% filter(component == 1)
  if (nrow(sm1) == 0) next
  xn <- unique(sm1$X_set)
  yn <- unique(sm1$Y_set)
  if (xn %in% X_names && yn %in% Y_names) {
    R_mat[xn, yn] <- sm1$canonical_correlation
    P_mat[xn, yn] <- sm1$p_multicca
  }
}

# Export original matrix (Keep)
R_df <- as.data.frame(R_mat) %>% tibble::rownames_to_column("X_set")
P_df <- as.data.frame(P_mat) %>% tibble::rownames_to_column("X_set")
write.csv(R_df, file = file.path(out_dir, "Matrix_firstComp_R.csv"), row.names = FALSE)
write.csv(P_df, file = file.path(out_dir, "Matrix_firstComp_p.csv"), row.names = FALSE)

# ---------- Added: Export (Xs+Ys)*(Xs+Ys) Square Matrix ----------
Combined_names <- c(names(X_sets), names(Y_sets))
R_full <- matrix(NA_real_, nrow = length(Combined_names), ncol = length(Combined_names),
                 dimnames = list(Combined_names, Combined_names))
P_full <- matrix(NA_real_, nrow = length(Combined_names), ncol = length(Combined_names),
                 dimnames = list(Combined_names, Combined_names))

for (sm in summary_rows) {
  if (!"component" %in% names(sm)) next
  sm1 <- sm %>% filter(component == 1)
  if (nrow(sm1) == 0) next
  xn <- unique(sm1$X_set)
  yn <- unique(sm1$Y_set)
  if (xn %in% Combined_names && yn %in% Combined_names) {
    r <- sm1$canonical_correlation[1]
    p <- sm1$p_multicca[1]
    R_full[xn, yn] <- r
    P_full[xn, yn] <- p
    # Symmetric fill (if full square matrix is needed)
    if (is.na(R_full[yn, xn])) R_full[yn, xn] <- r
    if (is.na(P_full[yn, xn])) P_full[yn, xn] <- p
  }
}
R_full_df <- as.data.frame(R_full) %>% tibble::rownames_to_column("Set")
P_full_df <- as.data.frame(P_full) %>% tibble::rownames_to_column("Set")
write.csv(R_full_df, file = file.path(out_dir, "Matrix_firstComp_R_full.csv"), row.names = FALSE)
write.csv(P_full_df, file = file.path(out_dir, "Matrix_firstComp_p_full.csv"), row.names = FALSE)

# =========================================================================
# Part 6: Bootstrap Loading Stability Analysis (X1 vs Y1 only)
# =========================================================================
message("\n==========================================================")
message("Running sCCA Bootstrap Loading Stability Analysis (X1 vs Y1)...")
message("==========================================================\n")

# Get optimal penalty parameters from the original fit of X1 vs Y1
x_name <- "X1_Network_Topology"
y_name <- "Y1_Clinical"

Xmat_raw <- to_numeric_matrix(df, X_sets[[x_name]])
Ymat_raw <- to_numeric_matrix(df, Y_sets[[y_name]])

# Preprocess identically to fit_true_and_penalties
keepX <- which(apply(Xmat_raw, 2, var, na.rm = TRUE) > min_var_global)
keepY <- which(apply(Ymat_raw, 2, var, na.rm = TRUE) > min_var_global)
X_clean <- as.matrix(Xmat_raw[, keepX, drop = FALSE])
Y_clean <- as.matrix(Ymat_raw[, keepY, drop = FALSE])

cc <- complete.cases(cbind(X_clean, Y_clean))
X_cc <- X_clean[cc, , drop = FALSE]
Y_cc <- Y_clean[cc, , drop = FALSE]

X_scaled <- scale_cols(X_cc)
Y_scaled <- scale_cols(Y_cc)

# Run CCA.permute on the true data to get the optimal penalties
set.seed(seed_select_global)
per_true <- PMA::CCA.permute(
  x = X_scaled, z = Y_scaled, typex = "standard", typez = "standard",
  nperms = nperms_select_global, trace = FALSE,
  penaltyxs = penaltyxs_ratio,
  penaltyzs = penaltyzs_ratio
)

opt_penaltyx <- per_true$bestpenaltyx
opt_penaltyz <- per_true$bestpenaltyz
message(sprintf("Optimal penalties selected: penaltyx = %.3f, penaltyz = %.3f", opt_penaltyx, opt_penaltyz))

# Fit observed (non-bootstrap) sCCA on the true data for all components
K_eff <- max(1, min(K_components_global, nrow(X_scaled)-1, ncol(X_scaled), ncol(Y_scaled)))
fit_obs <- PMA::CCA(
  x = X_scaled, z = Y_scaled, typex = "standard", typez = "standard",
  K = K_eff, standardize = FALSE,
  penaltyx = opt_penaltyx, penaltyz = opt_penaltyz,
  v = per_true$v.init
)

# Extract observed loadings
u_obs <- fit_obs$u # ncol(X_scaled) x K_eff
v_obs <- fit_obs$v # ncol(Y_scaled) x K_eff

# Bootstrap settings
n_boot <- 1000
n_subj <- nrow(X_scaled)

# We will store the loadings for each bootstrap run across all components
boot_u_array <- array(NA_real_, dim = c(n_boot, ncol(X_scaled), K_eff))
boot_v_array <- array(NA_real_, dim = c(n_boot, ncol(Y_scaled), K_eff))

# Match complete components jointly across X and Y, retaining raw estimates.
all_permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(x, function(i) cbind(i, all_permutations(x[x != i]))))
}
component_orders <- all_permutations(seq_len(K_eff))
reference <- rbind(u_obs, v_obs)
reference <- sweep(reference, 2, sqrt(colSums(reference^2)), "/")
boot_u_raw <- boot_u_array
boot_v_raw <- boot_v_array
boot_indices <- matrix(NA_integer_, n_boot, n_subj)
boot_match <- boot_sign <- boot_similarity <- matrix(NA_real_, n_boot, K_eff)
boot_error <- rep("", n_boot)

set.seed(789) # seed for bootstrap resampling

for (b in seq_len(n_boot)) {
  if (b %% 100 == 0) {
    message(sprintf("  Bootstrap iteration %d / %d ...", b, n_boot))
  }
  
  boot_idx <- sample(n_subj, replace = TRUE)
  boot_indices[b, ] <- boot_idx
  X_b <- X_scaled[boot_idx, , drop = FALSE]
  Y_b <- Y_scaled[boot_idx, , drop = FALSE]
  
  # Re-scale within the bootstrap sample
  X_b_scaled <- scale_cols(X_b)
  Y_b_scaled <- scale_cols(Y_b)
  
  fit_b <- tryCatch({
    PMA::CCA(
      x = X_b_scaled, z = Y_b_scaled, typex = "standard", typez = "standard",
      K = K_eff, standardize = FALSE,
      penaltyx = opt_penaltyx, penaltyz = opt_penaltyz, trace = FALSE
    )
  }, error = function(e) { boot_error[b] <<- conditionMessage(e); NULL })

  if (!is.null(fit_b)) {
    if (fit_b$K != K_eff || any(!is.finite(c(fit_b$u, fit_b$v)))) {
      boot_error[b] <- "Incomplete or nonfinite component weights"
      next
    }
    boot_u_raw[b, , ] <- fit_b$u
    boot_v_raw[b, , ] <- fit_b$v
    candidate <- rbind(fit_b$u, fit_b$v)
    similarity <- crossprod(reference, sweep(candidate, 2, sqrt(colSums(candidate^2)), "/"))
    scores <- apply(component_orders, 1, function(p) sum(abs(similarity[cbind(seq_len(K_eff), p)])))
    matched <- component_orders[which.max(scores), ]
    direction <- ifelse(similarity[cbind(seq_len(K_eff), matched)] < 0, -1, 1)
    boot_match[b, ] <- matched
    boot_sign[b, ] <- direction
    boot_similarity[b, ] <- abs(similarity[cbind(seq_len(K_eff), matched)])
    boot_u_array[b, , ] <- sweep(fit_b$u[, matched, drop = FALSE], 2, direction, "*")
    boot_v_array[b, , ] <- sweep(fit_b$v[, matched, drop = FALSE], 2, direction, "*")
  }
}

saveRDS(list(u = boot_u_array, v = boot_v_array, raw_u = boot_u_raw, raw_v = boot_v_raw,
             observed_u = u_obs, observed_v = v_obs, features_X = colnames(X_scaled),
             features_Y = colnames(Y_scaled), participant_ids = df[[id_col]][cc],
             indices = boot_indices, matched_component = boot_match, sign = boot_sign,
             similarity = boot_similarity, error = boot_error, n = n_subj, B = n_boot,
             seed = 789L, penalties = c(opt_penaltyx, opt_penaltyz), session = sessionInfo()),
        file.path(out_dir, "bootstrap_weights_X1_vs_Y1.rds"))
# Combine statistics for all K_eff components
boot_stats_list <- list()

for (k in seq_len(K_eff)) {
  u_obs_k <- u_obs[, k]
  v_obs_k <- v_obs[, k]
  
  boot_u_k <- boot_u_array[, , k]
  boot_v_k <- boot_v_array[, , k]
  
  abs_boot_u_k <- abs(boot_u_k)
  abs_boot_v_k <- abs(boot_v_k)
  
  # Calculate statistics for X features (Imaging)
  stats_u <- tibble(
    Feature = colnames(X_scaled),
    Observed_Loading = u_obs_k,
    Selection_Frequency = colMeans(boot_u_k != 0, na.rm = TRUE),
    Mean_Abs_Loading = colMeans(abs_boot_u_k, na.rm = TRUE),
    Bootstrap_SD = apply(abs_boot_u_k, 2, sd, na.rm = TRUE),
    CI_Lower_Abs = apply(abs_boot_u_k, 2, quantile, probs = 0.025, na.rm = TRUE),
    CI_Upper_Abs = apply(abs_boot_u_k, 2, quantile, probs = 0.975, na.rm = TRUE)
  ) %>%
    mutate(
      Signed_Mean_Loading = sign(Observed_Loading) * Mean_Abs_Loading,
      CI_Lower_Signed = pmin(sign(Observed_Loading) * CI_Lower_Abs, sign(Observed_Loading) * CI_Upper_Abs),
      CI_Upper_Signed = pmax(sign(Observed_Loading) * CI_Lower_Abs, sign(Observed_Loading) * CI_Upper_Abs),
      Side = "X",
      Component = k
    )
  
  # Calculate statistics for Y features (Clinical)
  stats_v <- tibble(
    Feature = colnames(Y_scaled),
    Observed_Loading = v_obs_k,
    Selection_Frequency = colMeans(boot_v_k != 0, na.rm = TRUE),
    Mean_Abs_Loading = colMeans(abs_boot_v_k, na.rm = TRUE),
    Bootstrap_SD = apply(abs_boot_v_k, 2, sd, na.rm = TRUE),
    CI_Lower_Abs = apply(abs_boot_v_k, 2, quantile, probs = 0.025, na.rm = TRUE),
    CI_Upper_Abs = apply(abs_boot_v_k, 2, quantile, probs = 0.975, na.rm = TRUE)
  ) %>%
    mutate(
      Signed_Mean_Loading = sign(Observed_Loading) * Mean_Abs_Loading,
      CI_Lower_Signed = pmin(sign(Observed_Loading) * CI_Lower_Abs, sign(Observed_Loading) * CI_Upper_Abs),
      CI_Upper_Signed = pmax(sign(Observed_Loading) * CI_Lower_Abs, sign(Observed_Loading) * CI_Upper_Abs),
      Side = "Y",
      Component = k
    )
  
  boot_stats_k <- bind_rows(stats_u, stats_v)
  boot_stats_list[[k]] <- boot_stats_k
  
  # ===== Plot stability for component k =====
  # Show only top 30 features (by absolute weight) selected in the observed true sCCA fit
  plot_data <- boot_stats_k %>% 
    filter(Observed_Loading != 0) %>%
    arrange(desc(abs(Observed_Loading))) %>%
    slice_head(n = 30) %>%
    arrange(Side, desc(abs(Observed_Loading)))
  
  if (nrow(plot_data) > 0) {
    # Whiskers end at the most extreme weights within 1.5 IQR of the box.
    bootstrap_long <- bind_rows(
      as_tibble(boot_u_k, .name_repair = ~colnames(X_scaled)) %>% mutate(Side = "X"),
      as_tibble(boot_v_k, .name_repair = ~colnames(Y_scaled)) %>% mutate(Side = "Y")
    ) %>% group_by(Side) %>% mutate(Resample = row_number()) %>% ungroup() %>%
      pivot_longer(-c(Side, Resample), names_to = "Feature", values_to = "Weight") %>%
      filter(is.finite(Weight)) %>%
      inner_join(plot_data %>% select(Feature, Side, Observed_Loading), by = c("Feature", "Side"))
    write.csv(bootstrap_long, file.path(out_dir, sprintf("bootstrap_weights_comp%d_X1_vs_Y1.csv", k)), row.names = FALSE)
    box_data <- bootstrap_long %>% group_by(Feature, Side, Observed_Loading) %>%
      summarise(Q1 = quantile(Weight, 0.25), Median = median(Weight),
                Q3 = quantile(Weight, 0.75),
                Lower_whisker = min(Weight[Weight >= Q1 - 1.5 * (Q3 - Q1)]),
                Upper_whisker = max(Weight[Weight <= Q3 + 1.5 * (Q3 - Q1)]),
                B_valid = n(), .groups = "drop")
    p_k <- ggplot(box_data, aes(x = reorder(Feature, abs(Observed_Loading)), fill = Side)) +
      geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
      geom_boxplot(aes(ymin = Lower_whisker, lower = Q1, middle = Median, upper = Q3, ymax = Upper_whisker),
                   stat = "identity", width = 0.6, staplewidth = 0.4, linewidth = 0.35, show.legend = FALSE) +
      geom_point(aes(y = Observed_Loading), shape = 18, color = "black", size = 2) +
      coord_flip() + facet_wrap(~Side, scales = "free_y") +
      scale_fill_manual(values = c("X" = "#FF6666", "Y" = "#6699FF"), guide = "none") +
      labs(x = NULL, y = "Canonical weight") +
      theme_minimal(base_size = 12) + theme(legend.position = "none", axis.text.y = element_text(size = 10, face = "bold"))

    plot_file <- file.path(out_dir, sprintf("bootstrap_loadings_stability_comp%d_X1_vs_Y1.png", k))
    ggsave(filename = plot_file, plot = p_k, width = 10, height = 7.5, dpi = 300)
    message(sprintf("Saved bootstrap plot for Component %d to: %s", k, plot_file))
  } else {
    message(sprintf("No features were selected (Observed_Loading != 0) for Component %d. Skipping plot.", k))
  }
}

boot_stats_all <- bind_rows(boot_stats_list)
stats_file <- file.path(out_dir, "bootstrap_loading_stability_X1_vs_Y1.csv")
write.csv(boot_stats_all, file = stats_file, row.names = FALSE)
message(sprintf("Saved bootstrap statistics CSV to: %s", stats_file))

