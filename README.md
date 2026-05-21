# Testing MONET method:
# Simulation-Based Method Comparison for Detecting Selection in Quantitative Traits

This repository contains simulation scripts, data preparation and analysis workflows, and results from comparing different methods of detecting selection on quantitative traits in structured populations.

## Project Overview

This study evaluates three different statistical methods for detecting selection on quantitative traits using simulated data. We simulated populations with varying structures and selection regimes using Quantinemo, then compared the performance of three methods (including a variation of MONET which considers habitat information):

1. **Qst-Fst** approach
2. **Driftsel** (From Ovaskainen and collaborators 2011 & 2013)
3. **MONET** (Log-Ratio of Ancestral Variances - from do O and collaborators 2025)
4. **MONET w/ habitat**

## Workflow

### 1. Simulations (Quantinemo)
We ran individual-based simulations using Quantinemo with:
- **Population structures**: 
  - Stepping Stone (SS) - 18 populations
  - Island Model (IM) - 18 and 9 populations
  - Hierarchical - 18 populations
- **Selection regimes**: Varying selection coefficients (wdiff), environmental variance (wvar), and correlation between selection and relatedness (correlation).

### 2. Data Preparation
All simulated data were standardized for analysis across methods:
- Generated F1 individuals through controlled breeding designs
- Calculated phenotypes for F1 generation
- Estimated coancestry matrices
- Ensured identical datasets for all three methods

### 3. Method Application
Applied each of the three methods to the exact same prepared datasets to obtain:
- P-values (Qst-Fst and Driftsel)
- S-values (MONET)

### 4. Performance Evaluation
Built ROC curves and calculated AUC (Area Under the Curve) to compare method performance across different scenarios.

## Repository Structure

```
Testing_MONET/
├── simulation_scripts/          # Quantinemo configuration files (.ini) for running simulations
│                                 # Includes all population structures (SS, IM, Hierarchical) with
│                                 # varying selection coefficients and environmental variance
│
├── method_testing_scripts/      # Scripts for data preparation and method application
│   ├── breeding_tests/          # R scripts for creating F1 generations under different
│   │                            # breeding designs (varying pop numbers, F1 numbers, sire/dam ratios)
│   └── testing_scripts/         # R scripts implementing each method (Qst-Fst, MONET, Driftsel)
│                                # and combining results across simulations
│
├── tools/                       # Utility R functions used across all analyses
│                                # (coancestry calculations, F1 creation, method implementations)
│
├── raw_data/                    # Combined results (p-values/s-values) from all method applications
│                                # TSV files containing test statistics for each locus and scenario
│
├── analysis_scripts/            # R scripts for generating ROC curves and calculating performance
│                                # metrics (AUC, TPR at specific FPR thresholds)
│
├── results/                     # Final results: AUC summaries (CSV), ROC 
│                                # and performance tables (CSV/TEX)
│
└── case_study/                  # Empirical re-analysis of Ramirez-Valiente et al. (2023)
    ├── envrionmental_data.r     # Download WorldClim/SoilGrids, compute climate PCs
    ├── comparing_myPCA_theirPCA.r  # Validate our PCs against published Figure 3
    ├── comparing_structures.r   # Compare oak Θ_P to simulated metapopulation structures
    ├── logAV_only.r             # MONET base model (no covariates) for each trait
    ├── logAV_local_environment_analysis.r  # MONET + PC_climate1/2 as fixed effects
    ├── results_logAV_only/      # Summary CSVs and .rds fits from base model
    └── results_logAV_env/       # Summary CSVs and .rds fits from environment model
```

### Case Study: Mediterranean Oaks (`case_study/`)

This directory contains an empirical re-analysis of the Ramirez-Valiente et al. (2023) dataset on two Mediterranean oak species (*Quercus faginea* and *Quercus lusitanica*) grown in common garden experiments.

**Scripts and what they do:**

- `standardized_pedigree_ramirez-valiente2022.R` / `standardized_phenotypes_ramirez-valiente2022.R` — prepare and standardize the raw pedigree and phenotype data from the original dataset into the formats expected by MONET.

- `ThetaP_TheM_ramirez-valiente2022.R` — computes the wild-population coancestry matrix (Θ_P) and the offspring kinship matrix (The.M) from genotype data. Writes `ramirez-valiente2022_{species}_{garden}_ThetaP.rds` and `ramirez-valiente2022_{species}_{garden}_TheM.rds`.

- `envrionmental_data.r` — downloads WorldClim 2.1 bioclimatic variables, SoilGrids 2.0 soil pH, and computes annual and summer moisture indices via the Hargreaves-Samani PET formula. Runs one PCA per species on all 22 environmental variables and writes the first two climate PCs (PC_climate1, PC_climate2) back into `pop_id_withhabitat.csv`. Also writes `pops_with_climate_soilpH_Im.csv` and `pca_climate_per_species.rds`.

- `comparing_myPCA_theirPCA.r` — validates that our recomputed climate PCs reproduce the population ordering reported in Ramirez-Valiente et al. (2023) Figure 3 (read off manually into `rv_coords.csv`). Computes per-axis Spearman rank correlations and a Procrustes correlation on the two-dimensional PC plane. Uses `pops_with_climate_soilpH_Im.csv` and `rv_coords.csv`.

- `comparing_structures.r` — visually compares the Θ_P structure of the oak metapopulations against simulated stepping-stone / island / hierarchical replicates. Produces `ThetaP_comparison_structures.pdf`.

- `logAV_only.r` — fits the base MONET model (no environmental covariates) to each trait, testing only for local vs. global adaptation via the log ratio of ancestral variances. Configured via environment variables (`MONET_SPECIES`, `MONET_GARDEN`, etc.). Results are written to `results_logAV_only/summary_logAV_only_{species}_{garden}.csv` and one `.rds` per trait.

- `logAV_local_environment_analysis.r` — fits the MONET model with PC_climate1 and PC_climate2 as fixed effects, jointly testing for adaptation signal and climate association. Results are written to `results_logAV_env/summary_logAV_env_{species}_{garden}.csv` and one `.rds` per trait.

- `job_logAV_alone_analysis.sh` / `job_logAV_envrionment_analysis.sh` — SLURM array job scripts that launch `logAV_only.r` and `logAV_local_environment_analysis.r` across all (species, garden) combinations on a cluster.

**Key data files:**

| File | Description |
|------|-------------|
| `pop_id_withhabitat.csv` | One row per individual; includes population coordinates and, after running `envrionmental_data.r`, the climate PCs. |
| `ramirez-valiente2022_{species}_{garden}_phenotypes_standardized.csv` | Standardized phenotype data (long format). |
| `ramirez-valiente2022_{species}_{garden}_pedigree_standardized.csv` | Pedigree linking offspring to wild-population parents. |
| `ramirez-valiente2022_{species}_{garden}_ThetaP.rds` | Wild-population coancestry matrix. |
| `ramirez-valiente2022_{species}_{garden}_TheM.rds` | Offspring kinship matrix. |
| `rv_coords.csv` | Per-population PC scores read from Ramirez-Valiente et al. (2023) Figure 3, used for validation only. |

---

### Directory Descriptions

**simulation_scripts/**  
Contains all Quantinemo `.ini` configuration files used to generate simulated datasets. These define population parameters (size, structure, migration rates), selection coefficients (wdiff), environmental variance (wvar), and genetic architecture for quantitative traits.

**method_testing_scripts/**  
Two subdirectories:
- `breeding_tests/`: Scripts that create F1 individuals from simulated parental populations using various breeding designs - `testing_data_preparation.r` is the script we used on the basic breeding desing, which we at times refer to as "full breeding desing".
- `testing_scripts/`: Scripts that apply Qst-Fst, MONET, MONET w/ habitat, and Driftsel methods to the prepared data, ensuring all methods analyze identical datasets

**tools/**  
Reusable R functions including:
- Coancestry matrix estimation from genotype data
- F1 generation from parental populations
- Driftsel implementation
- MONET implementation (with and without habitat information)
- Kinship calculations from pedigree data (unused in final versions)

**raw_data/**  
Tab-separated files containing the raw test statistics (p-values for Qst-Fst and MONET, s-values for Driftsel) for each trait across all simulated scenarios and replicates.

**analysis_scripts/**  
R scripts that generate ROC curves by comparing TPR vs FPR. Calculates Area Under the Curve (AUC) and True Positive Rates at specified False Positive Rate thresholds.

**results/**  
Final outputs including:
- AUC summary tables comparing method performance across scenarios
- ROC curve plots
- Performance metrics tables in CSV and LaTeX formats
