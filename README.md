# Epilepsy Network Topology Analysis

This repository contains the complete analysis pipeline for the study: **Mapping Individualized Dual-Axis Network Topology in Focal Epilepsy: Divergent Alterations in System Integrity, Integration, and Clinical Correlates**. The project leverages individualized resting-state functional MRI (rs-fMRI) networks to explore complementary, dual-axis network topology alterations across the spectrum of focal epilepsy and other common epilepsy syndromes, using Normative Modeling and Subtype and Stage Inference (SuStaIn).

## Background

The primary objective of this project is to characterize patient-specific alterations in network organization. We decompose system-level functional topology into two complementary axes:

- **System Integrity (Network Correspondence):** Quantifying the alignment (normativity) and idiosyncratic deviations (non-normativity) between individualized data-driven functional networks and canonical intrinsic connectivity systems across multiple parcellation schemes.
- **System Integration (k-hubness):** Quantifying cross-system communication and multi-functionality by measuring the multi-network participation of brain regions across overlapping functional systems.
- **Normative Modeling:** Using W-scores to map individual deviations from a healthy control reference model, adjusting for covariates like age and sex.
- **Multivariate Association:** Sparse Canonical Correlation Analysis (sCCA) to link dual-axis network metrics with clinical and cognitive phenotypes.
- **Disease Progression Modeling:** SuStaIn to identify distinct trajectories and stages of network correspondence disruption.

## Pipeline Overview

The analysis is organized into sequential steps (`S1` to `S12`), categorized by their analytical focus.

### 1. Data Organization & Preprocessing

- **S1_Org_outputs.m**:  
    Organizes raw derivatives (CSV, NIfTI) into a structured project directory. Aggregates subject metadata (including clinical and cognitive variables) and computes initial geometric overlaps (Dice coefficients) for subcortical structures.

### 2. Network Correspondence Analysis

- **S2_correspondence_analysis.m**:  
    Calculates the correspondence (overlap) between subject-specific, overlapping functional networks (derived via SPARK) and standard functional atlases (e.g., Yeo 17, Gordon, HCP-ICA). Defines metrics for Normativity and Non-normativity.
- **S3_correspondence_normative_modelling.m**:  
    Applies normative modeling (using `PCNtoolkit`) to the correspondence metrics. Computes W-scores (Z-scores adjusted for covariates) to quantify patient-specific expected deviations relative to healthy participants.
- **S4_correspondence_statistic.m**:  
    Performs group-level statistical comparisons of correspondence W-scores (Patients vs. Controls). Utilizes permutation-based max-T correction and generates consensus visualization plots across multiple atlases.

### 3. Hubness & Gray Matter Analysis

- **S5_hubmess_pipeline.m**:  
    Analyzes Functional k-hubness maps. It extracts ROI-based hubness metrics, performs normative modeling to generate W-scores, and maps statistical deviations (multi-network integration changes) on the brain surface.
- **S6_TJU_GM_pipeline.m**:  
    Processes structural MRI data for the TJU cohort via CAT12, extracting Gray Matter Volume (GMV) deviations using normative modeling to provide anatomical features for downstream, parallel SuStaIn modeling.

### 4. Multivariate & Progression Modeling

- **S7_prepare_data_for_CCA_SusStain.m**:  
    Aggregates computed features (Correspondence W-scores, Hubness W-scores, GMV W-scores, Clinical/Cognitive Demographics) into unified datasets for advanced multivariate modeling and progression inference.
- **S8_CCA.R**:  
    Runs Sparse Canonical Correlation Analysis (sCCA) to identify dissociable latent modes of association linking dual-axis brain network topology (correspondence and k-hubness) with clinical phenotypes and neurocognitive scores.
- **S9_PySuStaln_Correspondence.py** & **S9_PySuStaln_GMV.py**:  
    Executes the SuStaIn algorithm (Subtype and Stage Inference). Identifies distinct temporal progression trajectories of network correspondence loss and gray matter atrophy across the patient population.
- **S10_plot_pySuStaIn.m** & **S11_pySuStaIn_statisitic.m**:  
    Visualize and statistically evaluate SuStaIn outputs:
  - Subtype probability maps and staging progression diagrams.
  - Group-wise statistics and behavioral comparisons (e.g., PCA of cognition) of assigned subtypes and stages.

### 5. Validation

- **S12_Validation.m**:  
    Performs independent cross-cohort validation (comparing the TJU discovery dataset to the JLH validation dataset). Evaluates spatial reproducibility of correspondence and k-hubness alterations via Spin tests, and expands analysis to characterize distinct topology across other epilepsy syndromes (e.g., GGE, SeLECTS, Absence Epilepsy).

## Dependencies

### MATLAB

- **Statistics and Machine Learning Toolbox**
- **Bioinformatics Toolbox**
- **External Toolboxes** (included in `src/utils` or required externally):
  - `PCNtoolkit` (MATLAB wrapper)
  - `ENIGMA Toolbox` (for Spin tests/Surface plotting)
  - `PERMUTOOLS` (for max-T permutation inference)

### Python
>
> **Note:** When running [PCNtoolkit](https://github.com/amarquand/PCNtoolkit) and [pySuStaIn](https://github.com/ucl-pond/pySuStaIn), please configure the Python environment properly and add the corresponding conda Python interpreter to the MATLAB environment.

### R

- `PMA` (Penalized Multivariate Analysis for sCCA)

## Data Availability

The computation data and pre-calculated results can be downloaded from:

> **Zhang, Q., Arielle, D., Sharifzadeh Javidi, S., Ankeeta, A., Sperling, M., Zhang, Z., & Tracy, J. (2026). Mapping Individualized Dual-Axis Network Topology in Focal Epilepsy: Divergent Alterations in System Integrity, Integration, and Clinical Correlates [Data set]. Zenodo. <https://doi.org/10.5281/zenodo.18808544>**

Please download the dataset and extract it into the root directory of this repository. Note that the Zenodo repository also contains all pre-computed results.

## Usage

1. **Setup**: Ensure all data is organized in the root directory and initial outputs are structured via `S1`.
2. **Run Sequentially**: Execute scripts `S1` through `S12` in order.
    - MATLAB scripts should be run from the `src/pipeline` directory.
    - Python/R scripts are executed for specific multivariate and progression modeling (S8, S9).
3. **Configuration**: Check the parameter section of each script to adjust paths or model settings.

## Citation

If you use this code or pipeline, please cite:
> **Zhang Q., Dascal A., et al., "Mapping Individualized Dual-Axis Network Topology in Focal Epilepsy: Divergent Alterations in System Integrity, Integration, and Clinical Correlates", 2026**
