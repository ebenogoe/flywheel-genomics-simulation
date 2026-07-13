# Flywheel Genomics Simulation

Simulation scripts for two related studies. Set your working directory to the
relevant subfolder before running any scripts.

## breeding_simulation

Forward-time breeding simulation comparing three recurrent selection schemes
(Genetic Male Sterility, Chemical Sterilization, and Conventional Inbred Development)
using [AlphaSimR](https://cran.r-project.org/package=AlphaSimR).

| File | Purpose |
|------|---------|
| `new_main_run_all_scenarios_unique.R` | Main run script (30 replicates) |
| `functions_utils.R` | Shared helper functions |
| `male_sterility_scheme_gwas.R` | Genetic Male Sterility scheme |
| `chem_sterilization_scheme_gwas.R` | Chemical Sterilization scheme |
| `conventional_inbred_scheme_gwas.R` | Conventional Inbred Development scheme |
| `plot_individual_plots_true_mean.R` | Plotting |
| `quick_test.R` | Smoke test (2 replicates, reduced parameters) |

R packages required: `AlphaSimR`, `ggplot2`, `statgenGWAS`, `future`, `future.apply`, `progressr`, `patchwork`, `dplyr`

## gwas_simulation

Power and false positive rate simulation comparing three GWAS models across two
plant populations (SAP: Sorghum Association Panel, HBP: Haitian Breeding Population)
using [simplePHENOTYPES](https://github.com/samuelbfernandes/simplePHENOTYPES) and
[GAPIT](https://github.com/jiabowang/GAPIT).

| File | Purpose |
|------|---------|
| `pipeline_parallel.R` | Power simulation (100 replicates, parallelized) |
| `pipeline_fpr.R` | False positive rate simulation (100 null replicates) |
| `plot_summary_barplot.R` | Summary plots for power results |
| `plot_summary_barplot_fpr.R` | Summary plots for FPR results |

R packages required: `GAPIT`, `simplePHENOTYPES`, `dplyr`, `readr`, `stringr`, `future`, `furrr`, `progressr`
