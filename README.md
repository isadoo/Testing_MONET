# Testing MONET:

This repository is related to a paper under review.
In this work we present MONET - MOdel of Neutral Evolution of Traits. 
Check out the model's repository here: [MONET](https://github.com/isadoo/MONET)
This particular repository is for the part of the testing and applying MONET for testing different hypothesis. 

We tested two methods of MONET (logAV and habitat as a covariate) and other comparable methods previously published using simulated data
AND we analyzed two case studies (Meditaranian Oak & *Geranium carolinianum*)

  - [Project Overview](#project-overview)
  - [Repository Structure](#repository-structure)
    - [Directory Descriptions](#directory-descriptions)
    - [Breeding and sampling tests](#breeding-and-sampling-tests)
    - [Case Study: Mediterranean Oaks (`case_study_Oak/`)](#case-study-mediterranean-oaks-case_study_oak)
    - [Case Study: *Geranium carolinianum* (`case_study_Geranium/`)](#case-study-geranium-carolinianum-case_study_geranium)

## Project Overview

This study presents MONET and compares it to two different statistical models for detecting selection on quantitative traits using simulated data. One of the statistical models is based on the classic Qst-Fst, and the other is Driftsel. From Driftsel two test statistics can emerge: S-test and H-test. 

Similarly, from MONET one can get several test statistics for different hypotheses related to trait evolution. In our simulation tests we show two, the logAV (previously presented), and the habitat correlation test. 

1. **Qst-Fst** approach ([DOI: 10.1534/genetics.108.099812](https://academic.oup.com/genetics/article/183/3/1055/6063141?login=true))
2. **Driftsel** from Ovaskainen and collaborators 2011 ([DOI: 10.1534/genetics.111.129387](https://academic.oup.com/genetics/article/189/2/621/6063860)) & 2013 ([DOI: 10.1111/1755-0998.12111](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12111))
3. **Driftsel H-test** from Karhunen and collaborators 2014 ([DOI:10.1111/evo.12268](doi.org/10.1111/evo.12268))
4. **MONET: logAV**, Log-Ratio of Ancestral Variances from do O and collaborators 2025 ([
DOI: 10.1371/journal.pgen.1011871](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1011871))
5. **MONET: habitat** --> presented in this work!

We simulated populations with varying structures and selection regimes using Quantinemo ([DOI: 10.1093/bioinformatics/bty737](https://academic.oup.com/bioinformatics/article/35/5/886/5078477?login=true)).

## Repository Structure

```txt
Testing_MONET/
├── simulation_scripts/          # Quantinemo configuration files (.ini) for running simulations
│                                 # Includes all population structures (SS, IM, Hierarchical) with
│                                 # varying wdiff (Δθ) and wvar (ω_o)
│
├── method_testing_scripts/      # Scripts for data preparation and method application
│   ├── breeding_tests/          # R scripts for creating F1 generations under different
│   │                            # breeding designs (varying pop numbers, F1 numbers, sire/dam ratios)
│   └── testing_scripts/         # R scripts implementing each method (Qst-Fst, MONET, MONET w/ habitat, Driftsel S and H-test)
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
├── results/                     # Tables for the paper and SI (CSV/TEX)
│
├── case_study_Oak/              # Re-analysis of Ramirez-Valiente et al. (2023)
└── case_study_Geranium/         # Re-analysis of Shirk & Hamrick (2014)
```

### Directory Descriptions

**simulation_scripts/**  
Quantinemo `.ini` configuration files used to generate the simulated datasets: stepping stones (18 populations), island model (18 and 9 populations) and hierarchical (18 populations). They define the population parameters (size, structure, migration rates), the genetic architecture of the traits and the selection parameters: wdiff (Δθ, the distance between the most divergent optima: 0, 4.6, 10), wvar (ω_o, the width of the fitness function: 10, 22, 50) and the correlation pattern between optima and population structure (cline or parabola, called `sine` in the files, for stepping stones; grouped or swapped for hierarchical). [`experiment_config.yml`](method_testing_scripts/testing_scripts/experiment_config.yml) lists all the parameter combinations.

**method_testing_scripts/**  
Two subdirectories:

- `breeding_tests/`: Scripts that create F1 individuals from simulated parental populations using various breeding designs
  - [`testing_data_preparation.r`](method_testing_scripts/breeding_tests/testing_data_preparation.r) is the script we used on the basic breeding design, which we at times refer to as "full breeding design".
  - `testing_data_preparation_popnumber.r`, `_F1number.r`, `_familynumber_3by3.r` and `_familynumber_1by9.r` create the F1 generations of the sampling and breeding tests (see [Breeding and sampling tests](#breeding-and-sampling-tests)).

- `testing_scripts/`: Scripts that apply the different methods to the prepared data, ensuring all methods analyze identical datasets
  - Qst-Fst: [`testing_qstfst_analysis.r`](method_testing_scripts/testing_scripts/testing_qstfst_analysis.r)
  - MONET: [`testing_monet_analysis.r`](method_testing_scripts/testing_scripts/testing_monet_analysis.r)
  - MONET w/ habitat: [`testing_monet_whabitat_analysis.r`](method_testing_scripts/testing_scripts/testing_monet_whabitat_analysis.r)
  - Driftsel: [`testing_driftsel_analysis.r`](method_testing_scripts/testing_scripts/testing_driftsel_analysis.r)
  - Driftsel H-test: [`testing_driftsel_H_analysis.r`](method_testing_scripts/testing_scripts/testing_driftsel_H_analysis.r). Runs the Driftsel MCMC and computes both the H-test (`H.test`, with the simulated optima as habitat) and the S-test on the same chain (`S_value_Driftsel_Hrun`). `n_nearPD_H` counts the posterior samples where the drift covariance `2 G ⊗ θ` was not positive definite and was replaced by the nearest positive definite matrix.
  - MONET w/ habitat and the H-test take their habitat values from [`tools/habitat_optima.r`](tools/habitat_optima.r).

**tools/**  
Reusable R functions including:

- Coancestry matrix estimation from genotype data ([`calculate_coancestries.r`](tools/calculate_coancestries.r))
- F1 generation from parental populations ([`create_F1.r`](tools/create_F1.r))
- Driftsel implementation, S-test (`neut.test`) and H-test (`H.test`) ([`driftsel.r`](tools/driftsel.r)). `H.test` uses `Matrix::nearPD` when a posterior sample of `2 G ⊗ θ` is not positive definite and returns how often this happened (attribute `n.nearPD`).
- Habitat values (the simulated selection optima, per population structure, correlation and wdiff) used by MONET w/ habitat and the Driftsel H-test ([`habitat_optima.r`](tools/habitat_optima.r)).
- MONET implementation, with ([`monet_whabitat.r`](tools/monet_whabitat.r)) and without habitat information ([`monet.r`](tools/monet.r))

**raw_data/**  
`full_breeding_combined_review2026.tsv` and `all_methods_full_breeding_review2026.tsv` contain the full-breeding results with all five methods (MONET, MONET w/ habitat, Qst-Fst, Driftsel S and H-test) run on the same, seeded F1 data for SS / IM / hierarchical (18 populations), wdiff 4.6 and 10, wvar 10 / 22 / 50. The combined file adds the earlier results for the scenarios where MONET w/ habitat and the H-test do not apply (wdiff 0, no selection) and for the 9-population island models (column `source`). `3methods_full_breeding_combined_november.tsv` is the earlier November 2025 run (Qst-Fst, MONET and Driftsel S), used for Table 2 and the full-design rows of Table 3.
`breeding_run*.tsv`: results of the sampling and breeding tests (see [Breeding and sampling tests](#breeding-and-sampling-tests)).

**analysis_scripts/**  
`Results_ROC_Figures_monetwhabitat.r` makes the full-breeding ROC figures and AUC table used in the paper (Qst-Fst, MONET, MONET w/ habitat, Driftsel S and H-test), from `raw_data/full_breeding_combined_review2026.tsv`.

How each statistic is turned into a detection:

- Qst-Fst, MONET and MONET w/ habitat give p-values: selected if p < threshold.
- Driftsel S (average over MCMC samples of `pchisq` of the Mahalanobis distance, `mean(neut.test(...))`) and Driftsel H (`H.test`: proportion of MCMC samples where the observed Mantel statistic between population effects and habitat is larger than one drawn from the drift model `2 G ⊗ θ`) are treated the same way: folded, `min(x, 1 - x)`, and the p thresholds are mapped to folded thresholds with a normal of sd 0.10 (`S_fold_from_p`). Under drift both are centred on 0.5; H is wider than S in SS and hierarchical (sd 0.27 / 0.23 vs 0.13 / 0.15) and similar in IM (~0.2).
- The p grid goes down to 1e-7 (log-spaced below 1e-4) so that the mapped folded thresholds reach 0.
- Karhunen et al. (2014) cut-offs S or H < 0.05 / > 0.95 correspond to p < 0.10 two-sided. The calibration check at the end of `Results_ROC_Figures_monetwhabitat.r` also prints the neutral FPR with these cut-offs.

`AUC_total_and_partial.r` computes the partial AUC over FPR 0–0.1 and the total AUC for every scenario (IM 18 pops, SS, hierarchical) and method, from `raw_data/full_breeding_combined_review2026.tsv`. The ROC curves are the exact empirical ones from the raw statistics (`pROC`, selected vs neutral trait; p-values, folded S and folded H, smaller = more extreme), with no threshold grid. The partial AUC is divided by 0.1, so it is the mean power over FPR 0–0.1 (1 = perfect, 0.05 = random); the total AUC is over FPR 0–1 (1 = perfect, 0.5 = random). Outputs (tables for the SI, best method per row in bold): `results/partialAUC_FPR01_review2026.tex`/`.csv` and `results/totalAUC_FPR_review2026.tex`/`.csv` (the total AUC agrees with `AUC_summary_table_whabitat.csv` to within 0.02).

`Table2_TPR_at_FPR005.r` makes Table 2 of the paper (TPR at an FPR of 0.05 for every scenario with 18 subpopulations) from `raw_data/3methods_full_breeding_combined_november.tsv`: ROC curves on a grid of p-value thresholds, with the TPR interpolated at FPR = 0.05. Outputs: `results/TPR_at_FPR_005_all_scenarios_november2025.csv` and `results/TPR_at_FPR_005_table_november2025.tex`. `Table3_sampling_breeding.r` makes Table 3 (see [Breeding and sampling tests](#breeding-and-sampling-tests)). Both are run from `analysis_scripts/`.

**results/**  
The tables written by the scripts in `analysis_scripts/` (CSV and LaTeX).

### Breeding and sampling tests

These tests (Table 3 of the paper) measure how the sampling and common garden design affect power. Each test changes one aspect of the full design (18 subpopulations, 5 × 5 crossing, 50 F1 per subpopulation), for wdiff 4.6 and wvar 10.

| Test | Script (`method_testing_scripts/breeding_tests/`) | Change from the full design | Raw data (`raw_data/`) |
|------|------|------|------|
| Number of subpopulations | `testing_data_preparation_popnumber.r` | 9 consecutive subpopulations: 4 to 12 (stepping stones), 6 to 14 (hierarchical). Island model: the separate 9-subpopulation simulation | `breeding_run1_popnumber.tsv` |
| Number of F1 | `testing_data_preparation_F1number.r` | 25 of the 50 F1 per subpopulation, drawn at random | `breeding_run2_F1number.tsv` |
| 3 × 3 crossing | `testing_data_preparation_familynumber_3by3.r` | 3 dams × 3 sires (6 of the 10 parents), 18 F1 per subpopulation | `breeding_run3_3x3.tsv` |
| 1 × 9 crossing | `testing_data_preparation_familynumber_1by9.r` | 1 dam × 9 sires, 18 F1 per subpopulation; no Qst-Fst (design not supported) | `breeding_run4_1x9.tsv` |

[`Table3_sampling_breeding.r`](analysis_scripts/Table3_sampling_breeding.r) computes the TPR at an FPR of 0.05 as for Table 2, with the full-design rows from `3methods_full_breeding_combined_november.tsv`, and writes `results/Table3_sampling_breeding_TPR_FPR005.csv` and `.tex`.

### Case Study: Mediterranean Oaks (`case_study_Oak/`)

This directory contains a re-analysis of the Ramirez-Valiente et al. (2023; [DOI: 10.1111/mec.16816](https://pubmed.ncbi.nlm.nih.gov/36479963/)) dataset on two Mediterranean oak species (*Quercus faginea* and *Quercus lusitanica*) grown in common garden experiments.

**Scripts and what they do:**

- `standardized_pedigree_ramirez-valiente2022.R` / `standardized_phenotypes_ramirez-valiente2022.R` — prepare and standardize the raw pedigree and phenotype data from the original dataset into the formats expected by MONET.

- [`ThetaP_TheM_ramirez-valiente2022.R`](case_study_Oak/ThetaP_TheM_ramirez-valiente2022.R) — computes the wild-population coancestry matrix (Θ_P) and the offspring kinship matrix (The.M) from genotype data. Writes `ramirez-valiente2022_{species}_{garden}_ThetaP.rds` and `ramirez-valiente2022_{species}_{garden}_TheM.rds`.

- [`envrionmental_data.r`](case_study_Oak/envrionmental_data.r) — downloads WorldClim 2.1 bioclimatic variables, SoilGrids 2.0 soil pH, and computes annual and summer moisture indices via the Hargreaves-Samani PET formula. Runs one PCA per species on all 22 environmental variables and writes the first two climate PCs (PC_climate1, PC_climate2) back into `pop_id_withhabitat.csv`. Also writes `pops_with_climate_soilpH_Im.csv` and `pca_climate_per_species.rds`.

- [`comparing_myPCA_theirPCA.r`](case_study_Oak/comparing_myPCA_theirPCA.r) — validates that our recomputed climate PCs reproduce the population ordering reported in Ramirez-Valiente et al. (2023) Figure 3 (read off manually into `rv_coords.csv`). Computes per-axis Spearman rank correlations and a Procrustes correlation on the two-dimensional PC plane. Uses `pops_with_climate_soilpH_Im.csv` and `rv_coords.csv`.

- [`comparing_structures.r`](case_study_Oak/comparing_structures.r) — visually compares the Θ_P structure of the oak metapopulations against simulated stepping-stone / island / hierarchical replicates. Produces `ThetaP_comparison_structures.pdf`.

- [`logAV_only.r`](case_study_Oak/logAV_only.r) — fits the base MONET model (no environmental covariates) to each trait, testing only for local vs. global adaptation via the log ratio of ancestral variances. Configured via environment variables (`MONET_SPECIES`, `MONET_GARDEN`, etc.). Results are written to `results_logAV_only/summary_logAV_only_{species}_{garden}.csv` and one `.rds` per trait.

- [`logAV_local_environment_analysis.r`](case_study_Oak/logAV_local_environment_analysis.r) — fits the MONET model with PC_climate1 and PC_climate2 as fixed effects, jointly testing for adaptation signal and climate association. Results are written to `results_logAV_env/summary_logAV_env_{species}_{garden}.csv` and one `.rds` per trait.

- [`job_logAV_alone_analysis.sh`](case_study_Oak/job_logAV_alone_analysis.sh) / [`job_logAV_envrionment_analysis.sh`](case_study_Oak/job_logAV_envrionment_analysis.sh) — SLURM array job scripts that launch `logAV_only.r` and `logAV_local_environment_analysis.r` across all (species, garden) combinations on a cluster.

**Key data files:**

| File | Description |
|------|-------------|
| `pop_id_withhabitat.csv` | One row per individual; includes population coordinates and, after running `envrionmental_data.r`, the climate PCs. |
| `ramirez-valiente2022_{species}_{garden}_phenotypes_standardized.csv` | Standardized phenotype data (long format). |
| `ramirez-valiente2022_{species}_{garden}_pedigree_standardized.csv` | Pedigree linking offspring to wild-population parents. |
| `ramirez-valiente2022_{species}_{garden}_ThetaP.rds` | Wild-population coancestry matrix. |
| `ramirez-valiente2022_{species}_{garden}_TheM.rds` | Offspring kinship matrix. |
| `rv_coords.csv` | Per-population PC scores read from Ramirez-Valiente et al. (2023) Figure 3, used for validation only. |

### Case Study: *Geranium carolinianum* (`case_study_Geranium/`)

This directory contains an empirical re-analysis of the Shirk and Hamrick (2014; [DOI: 10.1111/evo.12474](https://doi.org/10.1111/evo.12474)) common-garden dataset on *Geranium carolinianum*, an annual plant native to the United States and introduced to China. The original data are from Dryad ([DOI: 10.5061/dryad.nm41q](https://doi.org/10.5061/dryad.nm41q)).

**Design of the dataset:**

- 20 source populations: 10 from the native range (USA) and 10 from the introduced range (China; population codes starting with `CN-`).
- Seeds were collected per maternal plant in the wild, so the offspring are coded as **maternal half-sib families** (200 families, 1400 plants in experiment 1).
- Wild individuals (782) were genotyped at 28 allozyme loci, used to estimate the population coancestry matrix (Θ_P).
- The data come in two experiments (`exp1`, `exp2`) because the plants measured for morphological / life-history traits and those measured for stable isotopes are named differently or are different individuals:
  - `exp1` (`Gercar qst data.csv`): morphological and life-history traits.
  - `exp2` (`Gercar qst SI data.txt`): stable isotopes and leaf C/N (595 plants).

**Question.** Instead of a local environmental variable, the habitat covariate here is the **country** of origin (China vs USA, USA = reference). Two models are fitted, both with country as a fixed effect.

The **country-only model**, `Y ~ country + (1 | gr(pop, cov = 2Θ_P)) + (1 | gr(ind, cov = M))`, tests jointly:

1. local vs global adaptation: is the log ratio of ancestral variances significantly positive (local) or negative (global)?
2. country association: is the country coefficient (β) significantly different from zero, beyond the neutral expectation given by Θ_P (Q_CT − F_CT analogue)?

The **by-just-pop model**, `Y ~ country + (1 | gr(pop, cov = 2Θ_P, by = country)) + (1 | gr(ind, cov = M))`, estimates the among-population variance (V_AB) separately for each country, with a single within-population variance (V_AW) for both. Population effects keep their coancestry between countries (cov(u_p, u_q) = sd_c(p) · sd_c(q) · 2Θ_P[p, q]). It gives one log ratio per country (`logAV_USA` = log(V_AB_USA / V_AW), `logAV_China` = log(V_AB_China / V_AW)), their difference (`logAV_diff` = log(V_AB_China / V_AB_USA)), and the country β as above.

Each model has a full version (`*_analysis.r`, individual effects kept) and a **marginal version** (`*_marginal.r`). Because M is block-diagonal by population, each population block can be rotated to the eigenbasis of M, where observations are independent with variance sd_A² · λ + σ². This integrates the individual effects out of the sampler. It is exact only for gaussian traits, so the Poisson and Student-t traits are always fitted with the full version.

**Scripts and what they do** (in `Shirk_Geranium_2014_EFLP_pedigree/`):

- [`standardized_pedigree_shirk2014.R`](case_study_Geranium/Shirk_Geranium_2014_EFLP_pedigree/standardized_pedigree_shirk2014.R) — builds the half-sib-dam pedigree (`id`, `sire`, `dam`, `sire_pop`, `dam_pop`; each offspring gets a unique unknown sire) for each experiment. Writes `shirk2014_{exp}_pedigree_standardized.csv`.
- [`standardized_phenotypes_shirk2014.R`](case_study_Geranium/Shirk_Geranium_2014_EFLP_pedigree/standardized_phenotypes_shirk2014.R) — reshapes the phenotypes to long format (`individual`, `block`, `trait`, `trait_id`). Writes `shirk2014_{exp}_phenotypes_standardized.csv`.
- [`dosage_genotypes_shirk2014.R`](case_study_Geranium/Shirk_Geranium_2014_EFLP_pedigree/dosage_genotypes_shirk2014.R) — converts the allozyme genotypes into allele dosages (one column per allele). Writes `shirk2014_genotypes_dosage.csv`.
- [`ThetaP_TheM_shirk2014.R`](case_study_Geranium/Shirk_Geranium_2014_EFLP_pedigree/ThetaP_TheM_shirk2014.R) — computes Θ_P from the wild genotypes and the offspring kinship matrix (The.M) from the pedigree, with `calculate_coancestries()` (see [`tools/`](tools)). Writes `shirk2014_{exp}_ThetaP.rds`, `shirk2014_{exp}_TheM.rds` and the `ThetaP_*.png` / `TheM_*.png` heat maps for visual inspection.
- [`logAV_country_analysis.r`](case_study_Geranium/Shirk_Geranium_2014_EFLP_pedigree/logAV_country_analysis.r) — fits MONET ([`tools/monet.r`](tools/monet.r)) with `country` as a fixed effect, one trait at a time. Configured through environment variables: `MONET_EXPERIMENT` (`exp1` or `exp2`), `MONET_TRAIT` (optional, run a single trait), `MONET_BASEDIR`, `MONET_TOOLSDIR`, `MONET_OUTDIR`. MCMC: 4 chains, 12000 iterations (4000 warm-up), thin 2, `adapt_delta = 0.99`, `max_treedepth = 15`. Θ_P is restricted to the populations observed for each trait.
- [`job_logAV_country_analysis.sh`](case_study_Geranium/Shirk_Geranium_2014_EFLP_pedigree/job_logAV_country_analysis.sh) — SLURM array job, one task per (experiment, trait) pair (17 traits in `exp1`, 5 in `exp2`). Logs go to `country_logs/`.
- [`logAV_country_marginal.r`](case_study_Geranium/Shirk_Geranium_2014_EFLP_pedigree/logAV_country_marginal.r) / [`job_logAV_country_marginal.sh`](case_study_Geranium/Shirk_Geranium_2014_EFLP_pedigree/job_logAV_country_marginal.sh) — the country-only model in the marginal form, with the same priors, MCMC settings and standardization. Run for the traits that did not converge (`dff_days`, `juvlf_len`, `seedwt`) or did not finish in time (`harvest_ht`, `rose_dia`) with `logAV_country_analysis.r`. Logs go to `country_marginal_logs/`.
- [`logAV_country_by_just_pop_marginal.r`](case_study_Geranium/Shirk_Geranium_2014_EFLP_pedigree/logAV_country_by_just_pop_marginal.r) / [`job_logAV_country_by_just_pop_marginal.sh`](case_study_Geranium/Shirk_Geranium_2014_EFLP_pedigree/job_logAV_country_by_just_pop_marginal.sh) — the by-just-pop model in the marginal form, for the 17 gaussian traits. Logs go to `country_by_just_pop_marginal_logs/`.
- [`logAV_country_by_just_pop_analysis.r`](case_study_Geranium/Shirk_Geranium_2014_EFLP_pedigree/logAV_country_by_just_pop_analysis.r) / [`job_logAV_country_by_just_pop_analysis.sh`](case_study_Geranium/Shirk_Geranium_2014_EFLP_pedigree/job_logAV_country_by_just_pop_analysis.sh) — the by-just-pop model in full, for the 5 non-gaussian traits (`boltA`, `boltB`, `boltC`, `totC`, `totN`). For the Poisson traits, V_AW is on the latent (log) scale. Logs go to `country_by_just_pop_logs/`.

All these scripts are configured through the same environment variables and use the same MCMC settings as `logAV_country_analysis.r`.

**Trait preparation** (done inside each analysis script):

- `germinated` and `planted` are protocol dates, not phenotypes (planting is always 2 days after scarification), and are dropped.
- `dff` and `juvlf_date` are stored as calendar-date serials in the Dryad file, not as durations. They are converted to durations following the original paper: `dff_days = dff − germinated` (days to first flower, corrected for germination time) and `juvlf_days = juvlf_date − planted` (juvenile establishment time).
- `repro_ratio = repwt / totalAG` (reproductive ratio) is added because it is the allocation trait tested in the original paper; `repwt` is kept as an additional trait.
- `boltA`, `boltB`, `boltC` (numbers of bolting stems) are counts: fitted with a Poisson family and not standardized.
- `totC` and `totN` (`exp2`) are fitted with a Student-t family: a few plants have %C / %N values far outside the plausible range while their isotope values look normal, which suggests transcription errors. Since this is published data, these plants are kept and the heavy-tailed likelihood down-weights them instead of removing them.
- All other traits use a Gaussian family on standardized values, with exponential(1) priors on `sigma` and `sd`. Replicate measurements of the same plant are averaged.

| Experiment | Traits |
|------------|--------|
| `exp1` | `boltA`, `boltB`, `boltC`, `dff_days`, `harvest_ht`, `juvlf_days`, `juvlf_len`, `lfarea`, `lfmass`, `repro_ratio`, `repwt`, `rootwt`, `rose_dia`, `seedwt`, `SLA`, `totalAG`, `vegwt` |
| `exp2` | `delC13`, `delN15`, `totC`, `totN`, `weight` |

**Outputs, country-only model** (`results_logAV_country/` and `results_logAV_country_marginal/`):

- `monet_{exp}_{trait}.rds` — full MONET result object for each trait (in `results_logAV_country_marginal/`: a data frame of posterior draws with the same columns as `monet()$sampling`, plus `sigma`).
- `summary_logAV_country_{exp}_{trait}.csv` — one row per trait: `n_individuals`, `n_populations`, log ratio (mean, median, 95% CI, p-value), country β (mean, median, 95% CI, p-value), and convergence diagnostics (`n_divergent`, `max_rhat`). Same columns in both directories.
- `summary_logAV_country_{exp}.csv` — combined summary, written only when all traits of an experiment are run in one job (without `MONET_TRAIT`).

**Outputs, by-just-pop model** (`results_logAV_country_by_just_pop/` for the non-gaussian traits and `results_logAV_country_by_just_pop_marginal/` for the gaussian traits, same file names and columns):

- `monet_country_by_just_pop_{exp}_{trait}.rds` — posterior draws: fixed effects, SDs, `V_AB_USA`, `V_AB_China`, `V_AW`, `log_ratio_USA`, `log_ratio_China` and `log_ratio_diff`.
- `summary_logAV_country_by_just_pop_{exp}_{trait}.csv` — one row per trait: `n_individuals`, `n_pop_USA`, `n_pop_China`, mean, median, 95% CI and p-value of `logAV_USA`, `logAV_China`, `logAV_diff` and `country_beta`, and convergence diagnostics (`n_divergent`, `max_rhat`).

**Key data files:**

| File | Description |
|------|-------------|
| `doi_10_5061_dryad_nm41q__v20140605/` | Original Dryad files, with their own READMEs describing each column. |
| `shirk2014_{exp}_pedigree_standardized.csv` | Half-sib-dam pedigree linking offspring to their source population. |
| `shirk2014_{exp}_phenotypes_standardized.csv` | Phenotypes in long format. |
| `shirk2014_genotypes_dosage.csv` | Allele dosages of the genotyped wild individuals. |
| `shirk2014_{exp}_ThetaP.rds` | Wild-population coancestry matrix (20 x 20). |
| `shirk2014_{exp}_TheM.rds` | Offspring kinship matrix. |
