# Flywheel Genomics Simulation

Simulation scripts for two related studies. Set your working directory to the
relevant subfolder before running any scripts.

---

## breeding_simulation

Forward-time breeding simulation comparing three recurrent selection schemes
(Genetic Male Sterility, Chemical Sterilization, and Conventional Inbred Development)
using [AlphaSimR](https://cran.r-project.org/package=AlphaSimR).

### Requirements

```r
install.packages(c("AlphaSimR", "ggplot2", "statgenGWAS",
                   "future", "future.apply", "progressr",
                   "patchwork", "dplyr"))
```

### Files

| File | Purpose |
|------|---------|
| `new_main_run_all_scenarios_unique.R` | Main run script (30 replicates, all schemes) |
| `functions_utils.R` | Shared helper functions |
| `male_sterility_scheme_gwas.R` | Genetic Male Sterility scheme |
| `chem_sterilization_scheme_gwas.R` | Chemical Sterilization scheme |
| `conventional_inbred_scheme_gwas.R` | Conventional Inbred Development scheme |
| `plot_individual_plots_true_mean.R` | Plotting (called by main and quick test scripts) |
| `quick_test.R` | Smoke test (2 replicates, reduced parameters) |

### Usage

**Quick smoke test (2 replicates, around 5 to 10 minutes):**

```r
source("quick_test.R")
```

**Full run (30 replicates, all founder sizes, several hours):**

```r
source("new_main_run_all_scenarios_unique.R")
```

Scheme toggles in `params`:

```r
run_ms3  = TRUE   # Genetic Male Sterility
run_chem = TRUE   # Chemical Sterilization
run_conv = TRUE   # Conventional Inbred Development
```

### Key parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `Ne` | 25 | Ancestral effective population size (coalescent) |
| `n.chr` | 10 | Chromosomes |
| `n.sites` | 300 | Segregating sites per chromosome |
| `base_pop_size` | 300 | Founder pool (Ne x 12) |
| `n.reps` | 30 | Independent replicates |
| `n.cycles_ms` | 50 | Recurrent cycles (Genetic Male Sterility and Chemical Sterilization) |
| `n.cycles_exe_rm` | 10 | Recurrent cycles (Conventional Inbred) |
| `n.gens` | 5 | Generations per cycle (Conventional Inbred only) |

Founder population generated with MaCS: theta = 2.5e-6/bp, rho = 1e-6/bp,
chromosome length 1e8 bp, `inbred = FALSE` (diverse diploid founders).

---

## gwas_simulation

Power and false positive rate simulation comparing three GWAS models across two
plant populations using [simplePHENOTYPES](https://github.com/samuelbfernandes/simplePHENOTYPES)
and [GAPIT](https://github.com/jiabowang/GAPIT).

### Populations

- **SAP**: Sorghum Association Panel
- **HBP**: Haitian Breeding Population

Genotype data for both populations are not included in this repository. The pipeline
expects preprocessed genotype objects and PCs/kinship matrices as inputs (see Usage below).

### Models compared

| Label | Description |
|-------|-------------|
| GLM | Naive model, no population structure or kinship correction |
| MLM + PC | Mixed model with principal components as fixed covariates |
| MLM + PC + K | Full mixed model with principal components and kinship matrix |

### Requirements

```r
# GAPIT (from GitHub)
devtools::install_github("jiabowang/GAPIT", force = TRUE)

# Remaining packages (from CRAN)
install.packages(c("simplePHENOTYPES", "dplyr", "readr",
                   "stringr", "future", "furrr", "progressr"))
```

### Files

| File | Purpose |
|------|---------|
| `pipeline_parallel.R` | Power simulation (100 replicates, parallelized via future/furrr) |
| `pipeline_fpr.R` | False positive rate simulation (100 null replicates) |
| `plot_summary_barplot.R` | Summary plots for power results |
| `plot_summary_barplot_fpr.R` | Summary plots for FPR results |

### Usage

The pipeline expects the following preprocessed inputs in the working directory:

- `sap_geno_df_MAF0.05.RData` and `hbp_geno_df_chibas_MAF0.05.RData`: filtered genotype matrices
- `SAP/SAP_PCs.rds` and `SAP/SAP_kinship.rds`: SAP principal components and kinship matrix
- `HBP/HBP_PCs.rds` and `HBP/HBP_kinship.rds`: HBP principal components and kinship matrix

**Power simulation:**

```r
source("pipeline_parallel.R")
```

Results are written to `gwas_power_results/`. Key outputs:
- `detection_results_all_reps.csv`: per-replicate detection results
- `power_summary.csv`: power by population, model, and window size

**FPR simulation:**

```r
source("pipeline_fpr.R")
```

Results are written to `gwas_fpr_results/`.

**Generate summary figures:**

```r
source("plot_summary_barplot.R")
source("plot_summary_barplot_fpr.R")
```

### Key parameters (pipeline_parallel.R)

| Parameter | Value | Description |
|-----------|-------|-------------|
| `N_REPS` | 100 | Independent replicates |
| `H2` | 0.5 | Simulated heritability |
| `ADD_EFFECT` | 1 | Additive QTN effect size |
| `WINDOW_SIZES` | 10, 25, 100 kb | Detection window sizes around true QTN |
| `N_PC` | 3 | Principal components used as covariates |
| `N_WORKERS` | 8 | Parallel workers (adjust to available cores) |
