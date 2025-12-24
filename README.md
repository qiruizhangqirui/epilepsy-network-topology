# Epilepsy Network Correspondence & Normative Modeling Pipeline

This repository contains a MATLAB/Python pipeline for analyzing functional network correspondence and computing W-scores using normative modeling (PCNtoolkit) in epilepsy patients.

## Project Structure

The project is organized as follows:

```text
/
├── data/                 # Input data (Subjects.xlsx, atlases, raw derivatives)
├── outputs/              # Generated results and figures
├── src/
│   ├── pipeline/         # Main analysis scripts (S1 - S5)
│   └── utils/            # Helper functions and Python wrappers
└── README.md             # This file
```

## Setup & Requirements

### Prerequisites
1.  **MATLAB** (with Statistics and Machine Learning Toolbox, Bioinformatics Toolbox).
2.  **Python 3.x** (with `pcntoolkit`, `numpy`, `scipy`, `pandas`).
3.  **PCNtoolkit**: Ensure the python environment where PCNtoolkit is installed is active or accessible by MATLAB.

### Configuration
The pipeline uses **relative paths** automatically. No manual `config.m` is required.
- Ensure all raw data (`Subjects.xlsx`, `NCT_Derivatives`, `SPARK_Derivatives`, `NCT_atlases`) are placed inside the `data/` folder.

## Usage Pipeline

Run the scripts in `src/pipeline/` in the following order:

### 1. Organize Inputs
**Script:** `src/pipeline/S1_Org_outputs.m`
- **Purpose**: Scans `NCT_Derivatives` (correspondence CSVs) and `SPARK_Derivatives` (Hubness NIfTIs).
- **Output**: Aggregates data into `outputs/All_derivatives_struct.mat`.

### 2. Correspondence Analysis
**Script:** `src/pipeline/S2_correspondence_analysis.m`
- **Purpose**: Computes network correspondence measures (Dice coefficients) between subject-specific networks and standard atlases (e.g., Yeo 17, Glasser, HCP ICA).
- **Output**: `outputs/Correspondence_measures.mat`.

### 3. Normative Modeling (Correspondence)
**Script:** `src/pipeline/S3_correspondence_normative_modelling.m`
- **Purpose**: Trains normative models (using Healthy Controls as baseline) and computes **W-scores** for correspondence measures.
- **Method**: Uses `PCNtoolkit` (Bayesian Linear Regression) via Python integration.
- **Output**: `outputs/Correspondence_Wscore.mat`.

### 4. Statistical Analysis
**Script:** `src/pipeline/S4_correspondence_statistic.m`
- **Purpose**: Performs group statistical analysis on W-scores.
- **Features**:
    - Boxplots per group/atlas.
    - Matrix plots of Normativity/Non-normativity.
    - Radar plots for 17-network profiles.
    - Cortical surface visualizations (using Schaefer 200 parcellation).
- **Output**: Figures and stats in `outputs/correspondence_statistic/`.

### 5. Hubness Analysis
**Script:** `src/pipeline/S5_hubmess_pipeline.m`
- **Purpose**: Separate pipeline for "Hubness" map analysis.
- **Steps**:
    1.  Extracts ROI means from Hubness maps.
    2.  Computes W-scores using normative modeling.
    3.  Generates statistical maps (Surface/Volume).
- **Output**: `outputs/hubness_statistic/`.

## Key Utility Functions
Located in `src/utils/`:
- `pcn_wscore.py`: Python wrapper for PCNtoolkit.
- `compute_wscore_pcn.m`: MATLAB interface to the Python wrapper.
- `calc_consensus_normativity.m`: Helper for consensus metric calculation.
- `plot_*.m`: Various plotting utilities for boxplots, surfaces, and matrices.
