# Epilepsy Network Topology Analysis

This repository contains the complete analysis pipeline for the study on functional-structural network correspondence and hubness comparisons in epilepsy. The project integrates multimodal MRI data (fMRI, sMRI) to explore network topology alterations in Focal Epilepsy (FE) and Temporal Lobe Epilepsy (TLE), using Normative Modeling and Subtype and Stage Inference (SuStaIn).

## Background

The primary objective of this project is to investigate the decoupling between functional and structural connectivity in epilepsy. We employ:

- **Network Correspondence**: Quantifying the overlap between individual functional networks and normative structural atlases.
- **Hubness Mapping**: Identifying critical network hubs and their disruption.
- **Normative Modeling**: Using W-scores to map individual deviations from a healthy control reference model, accounting for age and sex.
- **Multivariate Association**: Sparse Canonical Correlation Analysis (sCCA) to link network metrics with clinical variables.
- **Disease Progression Modeling**: SuStaIn to identify distinct neurodegenerative subtypes and stages.

## Pipeline Overview

The analysis is organized into sequential steps (`S1` to `S12`), categorized by their analytical focus.

### 1. Data Organization & Preprocessing

* **S1_Org_outputs.m**:  
    Organizes raw derivatives (CSV, NIfTI) into a structured project directory. It aggregates subject metadata and computes initial geometric overlaps (Dice coefficients) for subcortical structures.

### 2. Network Correspondence Analysis

* **S2_correspondence_analysis.m**:  
    Calculates the correspondence (overlap) between subject-specific functional networks and standard functional atlases (e.g., Yeo 7/17 networks). Defines metrics like "Normativity" and "Maximum Match".
- **S3_correspondence_normative_modelling.m**:  
    Applies normative modeling (using `PCNtoolkit`) to the correspondence metrics. Computes W-scores (Z-scores adjusted for covariates) to quantify patient-specific deviations.
- **S4_correspondence_statistic.m**:  
    Performs group-level statistical comparisons of W-scores (Patients vs. Controls). Includes False Discovery Rate (FDR) correction and generates visualization plots (Boxplots, Radar plots, Brain surfaces).

### 3. Hubness & Gray Matter Analysis

* **S5_hubmess_pipeline.m**:  
    Analyzes Functional Hubness maps. Similar to S3, it extracts ROI-based hubness metrics, performs normative modeling, and maps statistical deviations on the brain surface.
- **S6_TJU_GM_pipeline.m**:  
    External validation processing for the TJU cohort, aiming to replicate findings using Gray Matter Volume (GMV) or other structural metrics.

### 4. Multivariate & Progression Modeling

* **S7_prepare_data_for_CCA_SusStain.m**:  
    Aggregates all computed features (Correspondence W-scores, Hubness W-scores, GMV, Clinical Demographics) into a single dataset for advanced modeling.
- **S8_CCA.R**:  
    Runs Sparse Canonical Correlation Analysis (sCCA) to identify latent modes of association between brain network deviations and clinical phenotypes (e.g., duration of epilepsy, cognitive scores).
- **S9_PySuStaln_Correspondence.py**:  
    Executes the SuStaIn algorithm (Subtype and Stage Inference). It identifies distinct temporal progression patterns (subtypes) of network correspondence loss across the patient population.
- **S10_plot_pySuStaIn.m** & **S11_pySuStaIn_statisitic.m**:  
    Visualize SuStaIn outputs:
  - Subtype probability maps.
  - Staging progression diagrams.
  - Group-wise statistics of assigned subtypes and stages.

### 5. Validation

* **S12_Validation.m**:  
    Performs cross-cohort validation (e.g., comparing JLH and TJU sites). Includes:
  - Consistency checks of W-scores across sites.
  - Spatial correlation analysis (Spin tests) to verify topological similarity of findings.

## Dependencies

### MATLAB

- **Statistics and Machine Learning Toolbox**
- **Bioinformatics Toolbox**
- **External Toolboxes** (included in `src/utils` or required externally):
  - `PCNtoolkit` (MATLAB wrapper)
  - `ENIGMA Toolbox` (for Spin tests/Surface plotting)
  - `BrainNet Viewer` or `SurfStat` (for visualization)

### Python

- `numpy`, `pandas`, `scipy`
- `pySuStaIn` (for S9)
- `pcn_wscore` (custom wrapper for normative modeling)

### R

- `PMA` (Penalized Multivariate Analysis for CCA)
- `ggplot2` (for visualization)

## Usage

1. **Setup**: Ensure all data is placed in the root directory as specified in `S1`.
2. **Run Sequentially**: Execute scripts `S1` through `S12` in order.
    - MATLAB scripts should be run from the `src/pipeline` directory.
    - Python/R scripts are called for specific modeling steps (S8, S9).
3. **Configuration**: Check the "header" section of each script to adjust paths and parameters (e.g., `Project_Dir`, `Atlas_Name`).

## Citation

If you use this code or pipeline, please cite:
> **[Placeholder: Zhang et al., "Epilepsy Network Topology...", 2026]**
