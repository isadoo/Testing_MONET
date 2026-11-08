# Rerunning Traits with Poor MCMC Convergence

## Problem

Some traits from the initial monet analysis showed poor MCMC convergence diagnostics:
- **Rhat > 1.05**: Chains haven't converged (should be ≤ 1.01)
- **ESS < 300**: Insufficient effective sample size (want ≥ 400, ideally ≥ 1000)

These indicate the posterior samples are unreliable and should not be used for inference.

## Solution Strategy

### Improved MCMC Settings

The rerun uses **STRICT** settings that address common hierarchical model convergence issues:

| Parameter | Original | Strict | Reason |
|-----------|----------|--------|--------|
| `warmup` | 2000 | 4000 | More adaptation time for sampler to find optimal step size |
| `iter` | 5000 | 12000 | More total iterations |
| `thin` | 2 | 4 | Keep thinning to manage memory with longer chains |
| `adapt_delta` | 0.8 (default) | 0.99 | Smaller step sizes → fewer divergences (critical for hierarchical models) |
| `max_treedepth` | 10 (default) | 15 | Allow longer trajectories for complex posterior geometry |
| Effective draws | ~1500 | ~8000 | More usable posterior samples |

**Key insight**: The default `adapt_delta=0.8` is often too aggressive for hierarchical animal models with weak likelihoods and funnel-shaped posteriors. Increasing to 0.99 dramatically reduces divergent transitions.

### Why These Settings Work

1. **More warmup (4000)**: The sampler needs time to:
   - Estimate the mass matrix (covariance structure)
   - Find optimal step size
   - Escape problematic regions of parameter space

2. **Higher adapt_delta (0.99)**: 
   - Forces smaller NUTS step sizes
   - Reduces divergent transitions
   - Essential for hierarchical models with weak data

3. **Deeper trees (15)**:
   - Handles complex posterior geometry
   - Important when random effects have low variance
   - Prevents "maximum treedepth" warnings

4. **More iterations (12000)**:
   - Compensates for thinning
   - Provides ~8000 effective draws
   - Improves precision of credible intervals

## Workflow

### Quick Start

```bash
cd /work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/case_study

# Run complete workflow: identify failures → submit rerun jobs
./workflow_rerun_failures.sh
```

### Step-by-Step

#### 1. Identify Failed Traits

```bash
module load r-light
Rscript identify_convergence_failures.r
```

This scans all `.rds` files in `results_logAV_env/` and creates:
- `results_logAV_env/traits_needing_rerun.csv` — list of problematic traits

Criteria:
- max_rhat > 1.05
- min_ess < 300 (checks both Bulk_ESS and Tail_ESS)

#### 2. Submit Rerun Job

```bash
sbatch job_rerun_failed_traits.sh
```

Job specs:
- Array job (tasks 0-2 for faginea_GH, faginea_OC, lusitanica_GH)
- 72 hour time limit (strict MCMC takes ~3-5× longer)
- 32 GB memory (double the original)
- 4 CPUs per task

#### 3. Monitor Progress

```bash
# Check queue
squeue -u $USER

# Watch output (adjust task ID as needed)
tail -f environment_logs/rerun_60985088_0.stdout

# Check for completion
ls -lh results_logAV_env/summary_logAV_env_RERUN_*.csv
```

## Output Files

### Created Files

1. **`traits_needing_rerun.csv`**
   - List of problematic traits with diagnostics
   - Columns: species, garden, trait, max_rhat, min_ess, reason

2. **`monet_<species>_<garden>_<trait>.rds`** (overwritten)
   - Full monet results with strict settings
   - Includes `full_model` (brms object) for detailed diagnostics

3. **`summary_logAV_env_RERUN_<species>_<garden>.csv`**
   - Summary table with convergence diagnostics
   - New columns: `min_bulk_ess`, `min_tail_ess`
   - Compare to original `summary_logAV_env_<species>_<garden>.csv`

### Validation

After the job completes, verify convergence:

```r
library(brms)

# Load a rerun result
fit <- readRDS("results_logAV_env/monet_faginea_GH_AGR.TreatWW.rds")

# Check diagnostics
fit$convergence$max_rhat        # Should be < 1.01
summary(fit$full_model)         # Check all Rhat and ESS values

# Look for warnings
plot(fit$full_model)            # Trace plots should look like hairy caterpillars
pairs(fit$full_model)           # Check for funnel geometry
```

### Expected Results

**Before (original run)**:
- AGR.TreatWW: Rhat = 1.35, ESS as low as 9
- Gsmass.TreatD (lusitanica): Rhat = 1.32
- Thickness.TreatD (lusitanica): Rhat = 1.61

**After (strict settings)**:
- All Rhat ≤ 1.01
- All ESS ≥ 400 (most > 1000)
- No (or very few) divergent transitions
- Reliable credible intervals

## Notes on Alternative Strategies

### What About Non-Centered Parameterization?

If traits **still** fail with strict settings, the model may have funnel geometry (common when group-level variance is small). The brms default uses **centered** parameterization:

```
β_pop ~ Normal(μ, σ_pop)  # centered
```

For funnel problems, **non-centered** works better:
```
z_pop ~ Normal(0, 1)
β_pop = μ + σ_pop * z_pop  # non-centered
```

In brms, force this with:
```r
bf(trait ~ ..., center = FALSE)
```

However, this requires modifying `monet.r` and is only needed if strict MCMC settings don't resolve the issue.

### Why Not Just Drop Thinning?

Thinning (`thin=4`) throws away 75% of samples, which seems wasteful. However:
- **Memory**: Full chains take ~4× more RAM
- **Autocorrelation**: Doesn't help convergence (adaptation phase matters more)
- **ESS**: Effective sample size adjusts for autocorrelation

If memory isn't a constraint, you could try `thin=1` with `iter=5000, warmup=3000` for faster runs with similar effective sample size.

## Files in This Directory

| File | Purpose |
|------|---------|
| `identify_convergence_failures.r` | Scan .rds files, identify bad traits |
| `logAV_rerun_failed_traits.r` | Rerun with strict MCMC settings |
| `job_rerun_failed_traits.sh` | SLURM array job script |
| `workflow_rerun_failures.sh` | Complete workflow (identify → submit) |
| `RERUN_CONVERGENCE.md` | This file |

## References

- Stan Development Team (2024). "Runtime warnings and convergence problems"
  https://mc-stan.org/misc/warnings.html

- Betancourt (2017). "Diagnosing Biased Inference with Divergences"
  https://mc-stan.org/users/documentation/case-studies/divergences_and_bias.html

- McElreath (2020). *Statistical Rethinking*, Chapter 9 (MCMC diagnostics)
