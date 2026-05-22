# Testing MONET method:

# Simulation-Based Method Comparison for Detecting Selection in Quantitative Traits

This repository contains simulation scripts, data preparation and analysis workflows, and results from comparing different methods of detecting selection on quantitative traits in structured populations.

  - [Project Overview](#project-overview)
  - [Workflow](#workflow)
    - [1. Simulations (Quantinemo)](#1-simulations-quantinemo)
    - [2. Data Preparation](#2-data-preparation)
    - [3. Method Application](#3-method-application)
    - [4. Performance Evaluation](#4-performance-evaluation)
  - [Repository Structure](#repository-structure)
    - [Directory Descriptions](#directory-descriptions)
    - [Case Study: Mediterranean Oaks (`case_study/`)](#case-study-mediterranean-oaks-case_study)

## Project Overview

This study evaluates three different statistical methods for detecting selection on quantitative traits using simulated data. We simulated populations with varying structures and selection regimes using Quantinemo ([DOI: 10.1093/bioinformatics/bty737](https://academic.oup.com/bioinformatics/article/35/5/886/5078477?login=true)), then compared the performance of three methods (including a variation of MONET which considers habitat information):

1. **Qst-Fst** approach ([DOI: 10.1534/genetics.108.099812](https://academic.oup.com/genetics/article/183/3/1055/6063141?login=true))
2. **Driftsel** from Ovaskainen and collaborators 2011 ([DOI: 10.1534/genetics.111.129387](https://academic.oup.com/genetics/article/189/2/621/6063860)) & 2013 ([DOI: 10.1111/1755-0998.12111](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12111))
3. **MONET**, Log-Ratio of Ancestral Variances from do O and collaborators 2025 ([
DOI: 10.1371/journal.pgen.1011871](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1011871))
4. **MONET w/ habitat**

## Workflow

### 1. Simulations (Quantinemo)

We ran individual-based simulations using Quantinemo with:

- **Population structures**: 
  - Stepping Stone (SS) - 18 populations
  - Island Model (IM) - 18 and 9 populations
  - Hierarchical - 18 populations
- **Selection regimes**: Varying selection coefficients (wdiff), environmental variance (wvar), and correlation between selection and relatedness (correlation).

The [experiment_config.yml](method_testing_scripts/testing_scripts/experiment_config.yml) file
lists all the parameter combinations defining the simulations. Quantinemo configuration files (`.ini`) are in the [simulation_scripts](simulation_scripts) directory.

### 2. Data Preparation

All simulated data were standardized for analysis across methods:

- Generated F1 individuals through controlled breeding designs ([create_F1.r](tools/create_F1.r))
- Calculated phenotypes for F1 generation ([testing_data_preparation*.r](method_testing_scripts/breeding_tests/))
- Estimated coancestry matrices ([calculate_coancestries.r](tools/calculate_coancestries.r))
- Ensured identical datasets for all three methods

### 3. Method Application

Applied each of the three methods (under [tools](tools)) to the exact same prepared datasets to obtain:

- P-values (Qst-Fst and Driftsel)
- S-values (MONET)

### 4. Performance Evaluation

Built ROC curves and calculated Area Under the Curve (AUC) to compare method performance across different scenarios (see [analysis_scripts](analysis_scripts)).

## Repository Structure

```txt
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

### Directory Descriptions

**simulation_scripts/**  
Contains all Quantinemo `.ini` configuration files used to generate simulated datasets. These define population parameters (size, structure, migration rates), selection coefficients (wdiff), environmental variance (wvar), and genetic architecture for quantitative traits.

**method_testing_scripts/**  
Two subdirectories:

- `breeding_tests/`: Scripts that create F1 individuals from simulated parental populations using various breeding designs
  - [`testing_data_preparation.r`](method_testing_scripts/breeding_tests/testing_data_preparation.r) is the script we used on the basic breeding desing, which we at times refer to as "full breeding desing".
  - `testing_data_preparation_*number.r` represent the scripts creating F1 generations under different breeding designs (varying pop numbers, F1 numbers, sire/dam ratios)

- `testing_scripts/`: Scripts that apply the different methods to the prepared data, ensuring all methods analyze identical datasets
  - Qst-Fst: [`testing_qstfst_analysis.r`](method_testing_scripts/testing_scripts/testing_qstfst_analysis.r)
  - MONET: [`testing_monet_analysis.r`](method_testing_scripts/testing_scripts/testing_monet_analysis.r)
  - MONET w/ habitat: [`testing_monet_whabitat_analysis.r`](method_testing_scripts/testing_scripts/testing_monet_whabitat_analysis.r)
  - Driftsel: [`testing_driftsel_analysis.r`](method_testing_scripts/testing_scripts/testing_driftsel_analysis.r)

**tools/**  
Reusable R functions including:

- Coancestry matrix estimation from genotype data ([`calculate_coancestries.r`](tools/calculate_coancestries.r))
- F1 generation from parental populations ([`create_F1.r`](tools/create_F1.r))
- Driftsel implementation ([`driftsel.r`](tools/driftsel.r))
- MONET implementation, with ([`monet_whabitat.r`](tools/monet_whabitat.r)) and without habitat information ([`monet.r`](tools/monet.r))
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

### Case Study: Mediterranean Oaks (`case_study/`)

This directory contains an empirical re-analysis of the Ramirez-Valiente et al. (2023; [DOI: 10.1111/mec.16816](https://pubmed.ncbi.nlm.nih.gov/36479963/)) dataset on two Mediterranean oak species (*Quercus faginea* and *Quercus lusitanica*) grown in common garden experiments.

**Scripts and what they do:**

- `standardized_pedigree_ramirez-valiente2022.R` / `standardized_phenotypes_ramirez-valiente2022.R` — prepare and standardize the raw pedigree and phenotype data from the original dataset into the formats expected by MONET.

- [`ThetaP_TheM_ramirez-valiente2022.R`](case_study/ThetaP_TheM_ramirez-valiente2022.R) — computes the wild-population coancestry matrix (Θ_P) and the offspring kinship matrix (The.M) from genotype data. Writes `ramirez-valiente2022_{species}_{garden}_ThetaP.rds` and `ramirez-valiente2022_{species}_{garden}_TheM.rds`.

- [`envrionmental_data.r`](case_study/envrionmental_data.r) — downloads WorldClim 2.1 bioclimatic variables, SoilGrids 2.0 soil pH, and computes annual and summer moisture indices via the Hargreaves-Samani PET formula. Runs one PCA per species on all 22 environmental variables and writes the first two climate PCs (PC_climate1, PC_climate2) back into `pop_id_withhabitat.csv`. Also writes `pops_with_climate_soilpH_Im.csv` and `pca_climate_per_species.rds`.

- [`comparing_myPCA_theirPCA.r`](case_study/comparing_myPCA_theirPCA.r) — validates that our recomputed climate PCs reproduce the population ordering reported in Ramirez-Valiente et al. (2023) Figure 3 (read off manually into `rv_coords.csv`). Computes per-axis Spearman rank correlations and a Procrustes correlation on the two-dimensional PC plane. Uses `pops_with_climate_soilpH_Im.csv` and `rv_coords.csv`.

- [`comparing_structures.r`](case_study/comparing_structures.r) — visually compares the Θ_P structure of the oak metapopulations against simulated stepping-stone / island / hierarchical replicates. Produces `ThetaP_comparison_structures.pdf`.

- [`logAV_only.r`](case_study/logAV_only.r) — fits the base MONET model (no environmental covariates) to each trait, testing only for local vs. global adaptation via the log ratio of ancestral variances. Configured via environment variables (`MONET_SPECIES`, `MONET_GARDEN`, etc.). Results are written to `results_logAV_only/summary_logAV_only_{species}_{garden}.csv` and one `.rds` per trait.

- [`logAV_local_environment_analysis.r`](case_study/logAV_local_environment_analysis.r) — fits the MONET model with PC_climate1 and PC_climate2 as fixed effects, jointly testing for adaptation signal and climate association. Results are written to `results_logAV_env/summary_logAV_env_{species}_{garden}.csv` and one `.rds` per trait.

- [`job_logAV_alone_analysis.sh`](case_study/job_logAV_alone_analysis.sh) / [`job_logAV_envrionment_analysis.sh`](case_study/job_logAV_envrionment_analysis.sh) — SLURM array job scripts that launch `logAV_only.r` and `logAV_local_environment_analysis.r` across all (species, garden) combinations on a cluster.

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
