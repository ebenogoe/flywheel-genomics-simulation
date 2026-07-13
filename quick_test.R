library(AlphaSimR)
library(ggplot2)
library(statgenGWAS)
library(future)
library(future.apply)
library(progressr)
library(patchwork)
library(dplyr)

source("functions_utils.R")

# ── Reduced parameters for smoke test ────────────────────────────────────────
params <- list(
  Ne = 25, n.chr = 10, n.sites = 300,
  founder_sizes          = c(10),
  n.reps                 = 2,
  target.pop.size        = 2500,
  n.crosses              = 625,
  n.cycles_ms            = 3,
  n_ms_select            = 125,
  n_fertile_select       = 250,
  males_per_female       = 5,
  tracked_trait_id       = 2,
  max_attempts           = 10,
  n.cycles_exe_rm        = 2,
  n.gens                 = 5,
  select_top_fams_conv   = 25,
  n_select_vector  = c(NA, NA, 2, 125, 50, 25),
  n.progeny_vector = c(1, 20, 10, 20, 50, NA),
  use_selection_index      = TRUE,
  selection_index_weights  = c(0.8, 0, 0.2),
  run_gwas          = TRUE,
  gwas_cycle        = 2,
  gwas_threshold    = 0.05,
  gwas_random_size  = 500,
  gwas_distance     = 10000,
  n.QTL.trait3      = 1,
  use_qtl_filtering_t3       = TRUE,
  use_qtl_filtering_founders = FALSE,
  n.top.qtl                  = 1,
  estimate_ne        = TRUE,
  prop.markers       = 0.3,
  n.marker.sample    = 5,
  n.replicates.bgld  = 5,
  maf_threshold_value = 0.2,

  # ── Scheme toggle ─────────────────────────────────────────────────────────
  run_ms3  = TRUE,    # Genetic Male Sterility
  run_chem = TRUE,   # Chemical Sterilization
  run_conv = TRUE     # Conventional Inbred Development
)

initialize_globals <- function() {
  list2env(params, envir = .GlobalEnv)
  base_pop_size <<- Ne * 12
}
initialize_globals()

log_timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
log_filename  <- paste0("quick_test_", log_timestamp, ".txt")
sink(log_filename, split = TRUE, append = TRUE, type = "output")
cat("Quick test started:", format(Sys.time()), "\n\n")
start.time <- Sys.time()

# ── Inline create_base_population & sample_founders (same as main script) ────
create_base_population <- function(n.chr, n.sites, base_pop_size) {
  # Blueprint: 300 inbred lines, 10 chr, 1×10^8 bp/chr (manualGenLen=1 Morgan),
  #            Ne=25, mu=2.5×10^-8/bp, r=1×10^-8/bp (1 cM/Mb).
  # theta = 4·Ne·mu = 4·25·2.5e-8 = 2.5e-6 per bp
  # rho   = 4·Ne·r  = 4·25·1e-8   = 1e-6   per bp
  # -eN: Ne was 4× larger (=100) ten generations ago
  base_pop <- runMacs(
    nInd          = base_pop_size,
    nChr          = n.chr,
    segSites      = n.sites,
    inbred        = FALSE,
    manualCommand = paste(
      "100000000",                    # 1×10^8 bp per chromosome
      "-t", 4 * Ne * 2.5e-8,          # θ = 4·Ne·mu per bp; mu = 2.5×10^-8
      "-r", 4 * Ne * 1e-8,             # ρ = 4·Ne·r  per bp; r  = 1×10^-8
      "-eN", 10 / (4 * Ne), 100 / Ne  # Ne×4 ten generations ago
    ),
    manualGenLen  = 1,
    nThreads      = NULL
  )
  SP_local <- SimParam$new(base_pop)
  SP_local$addTraitADEG(nQtlPerChr = 100, mean = 0, var = 5, meanDD = 0.3, relAA = 0.3, varGxE = 2)
  SP_local$addTraitA(nQtlPerChr = c(0,0,0,0,0,0,0,0,0,1), mean = 0, var = 1)
  SP_local$addTraitAD(nQtlPerChr = 1, mean = 0, var = 1, meanDD = 0.5)
  SP_local$setVarE(h2 = c(0.5, 0.5, 0.5))
  SP_local$addSnpChip(nSnpPerChr = 0.1 * n.sites)
  list(base_pop = base_pop, SP = SP_local)
}

sample_founders <- function(base_pop_data, founder_size, use_qtl_filtering = FALSE, n.top.qtl = 1) {
  base_pop <- base_pop_data$base_pop
  SP       <- base_pop_data$SP
  founders <- base_pop[sample(seq_len(base_pop@nInd), founder_size)]
  initial_parents <- newPop(founders, simParam = SP)
  rm(founders); gc()
  list(population = initial_parents, SP = SP)
}

# ── Run ───────────────────────────────────────────────────────────────────────
plan(multisession, workers = 2)

cat("=== Creating base populations ===\n")
base_populations <- future_lapply(seq_len(params$n.reps), function(rep) {
  create_base_population(params$n.chr, params$n.sites, base_pop_size)
}, future.seed = TRUE)

all_results      <- list()
all_gwas_results <- list()

assemble <- function(results_list, label) {
  breeding <- Filter(Negate(is.null), lapply(results_list, `[[`, "breeding_results"))
  df <- do.call(rbind, breeding); df$scheme <- label
  gwas_list <- lapply(seq_along(results_list), function(i) {
    g <- results_list[[i]]$gwas_results
    if (is.null(g)) return(NULL)
    g$scheme <- label; g$founder_size <- size; g
  })
  list(breeding = df, gwas = Filter(Negate(is.null), gwas_list))
}

for (size in params$founder_sizes) {
  scheme_breeding <- list()

  if (run_ms3) {
    cat("\n=== SCHEME 1a: Genetic Male Sterility ===\n")
    source("male_sterility_scheme_gwas.R")
    ms_results <- future_lapply(seq_len(params$n.reps), function(rep) {
      fd <- sample_founders(base_populations[[rep]], size)
      assign("SP", fd$SP, envir = .GlobalEnv)
      r <- ms3_single_rep(rep, selected_founders = fd$population)
      rm("SP", envir = .GlobalEnv); r
    }, future.seed = TRUE)
    ms_a <- assemble(ms_results, "Genetic Male Sterility")
    all_gwas_results <- c(all_gwas_results, ms_a$gwas)
    scheme_breeding[["ms3"]] <- ms_a$breeding
  }

  if (run_chem) {
    cat("\n=== SCHEME 1b: Chemical Sterilization ===\n")
    source("chem_sterilization_scheme_gwas.R")
    chem_results <- future_lapply(seq_len(params$n.reps), function(rep) {
      fd <- sample_founders(base_populations[[rep]], size)
      assign("SP", fd$SP, envir = .GlobalEnv)
      r <- chem_single_rep(rep, selected_founders = fd$population)
      rm("SP", envir = .GlobalEnv); r
    }, future.seed = TRUE)
    chem_a <- assemble(chem_results, "Chemical Sterilization")
    all_gwas_results <- c(all_gwas_results, chem_a$gwas)
    scheme_breeding[["chem"]] <- chem_a$breeding
  }

  if (run_conv) {
    cat("\n=== SCHEME 2: Conventional Inbred ===\n")
    source("conventional_inbred_scheme_gwas.R")
    conv_results <- future_lapply(seq_len(params$n.reps), function(rep) {
      fd <- sample_founders(base_populations[[rep]], size)
      assign("SP", fd$SP, envir = .GlobalEnv)
      r <- conv_single_rep(rep, selected_founders = fd$population)
      rm("SP", envir = .GlobalEnv); r
    }, future.seed = TRUE)
    conv_a <- assemble(conv_results, "Conventional Inbred")
    all_gwas_results <- c(all_gwas_results, conv_a$gwas)
    scheme_breeding[["conv"]] <- conv_a$breeding
  }

  common_cols <- Reduce(intersect, lapply(scheme_breeding, colnames))
  combined    <- do.call(rbind, lapply(scheme_breeding, function(df) df[, common_cols]))
  combined$founder_size <- size
  all_results[[as.character(size)]] <- combined
}

final_results      <- do.call(rbind, all_results)
final_results_file <<- paste0("quick_test_results_", log_timestamp, ".csv")
write.csv(final_results, final_results_file, row.names = FALSE)
cat("\nBreeding results written to:", final_results_file, "\n")

if (length(all_gwas_results) > 0) {
  gwas_results_file <<- paste0("quick_test_gwas_", log_timestamp, ".csv")
  write.csv(do.call(rbind, all_gwas_results), gwas_results_file, row.names = FALSE)
  cat("GWAS results written to:", gwas_results_file, "\n")
}

cat("\n=== Generating plots ===\n")
source("plot_individual_plots_true_mean.R")
plot_all_schemes_fair_with_lm(final_results_file)
cat("Plots saved.\n")

end.time <- Sys.time()
cat("\nStarted:", format(start.time, "%X"),
    "\nEnded:  ", format(end.time,   "%X"),
    "\nDuration:", format(end.time - start.time), "\n\n")
cat("=== QUICK TEST COMPLETE ===\n")
sink()
cat("Log saved to:", log_filename, "\n")
