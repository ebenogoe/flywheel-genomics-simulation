# Flywheel Genomics — Breeding Scheme Simulation

Forward-time breeding simulation pipeline comparing Genetic Male Sterility (GMS),
Chemical Sterilization (CS), and Conventional Inbred Development in a recurrent
selection programme. Built with [AlphaSimR](https://cran.r-project.org/package=AlphaSimR).

## Requirements

R packages: `AlphaSimR`, `ggplot2`, `statgenGWAS`, `future`, `future.apply`,
`progressr`, `patchwork`, `dplyr`

## File overview

| File | Purpose |
|------|---------|
| `new_main_run_all_scenarios_unique.R` | Main orchestration script (30 replicates, all schemes) |
| `functions_utils.R` | Shared helpers: GWAS, selection, diversity metrics, crossing |
| `male_sterility_scheme_gwas.R` | Genetic Male Sterility scheme |
| `chem_sterilization_scheme_gwas.R` | Chemical Sterilization scheme |
| `conventional_inbred_scheme_gwas.R` | Conventional Inbred Development scheme |
| `plot_individual_plots_true_mean.R` | Plotting function (called by both scripts above) |
| `quick_test.R` | Smoke test — 2 reps, reduced parameters |

## Usage

**Quick smoke test** (2 reps, ~5–10 min):

```r
source("quick_test.R")
```

**Full run** (30 reps, all founder sizes — hours):

```r
source("new_main_run_all_scenarios_unique.R")
```

Scheme toggles in `params`:

```r
run_ms3  = TRUE   # Genetic Male Sterility
run_chem = TRUE   # Chemical Sterilization
run_conv = TRUE   # Conventional Inbred Development
```

## Key parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `Ne` | 25 | Ancestral effective population size (coalescent) |
| `n.chr` | 10 | Chromosomes |
| `n.sites` | 300 | Segregating sites per chromosome |
| `base_pop_size` | 300 | Founder pool (Ne × 12) |
| `n.reps` | 30 | Independent replicates |
| `n.cycles_ms` | 50 | Recurrent cycles (GMS/CS) |
| `n.cycles_exe_rm` | 10 | Recurrent cycles (Conventional) |

Founder population generated with MaCS: θ = 2.5 × 10⁻⁶/bp, ρ = 1 × 10⁻⁶/bp,
chromosome length 10⁸ bp, `inbred = FALSE` (diverse diploid founders).
