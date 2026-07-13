library(AlphaSimR)
library(ggplot2)
library(statgenGWAS)
library(future)
library(future.apply)
library(progressr)
library(patchwork)
library(dplyr)

source("functions_utils.R")


# ============================================================
# PARAMETERS
# ============================================================

params <- list(

  # ── Step 1: Founder pool (RunMacs) ──────────────────────────────────────────
  Ne      = 25,    # Effective population size
  n.chr   = 10,    # Chromosomes
  n.sites = 300,   # Segregating sites per chromosome
  # Traits:
  #   Trait 1 — Yield:           polygenic, 100 QTL/chr, h2 = 0.5
  #   Trait 2 — ms3 locus:       A model, 1 QTL on chr 10 (male sterility gene)
  #   Trait 3 — Stress response: oligogenic, 1 QTL/chr, h2 = 0.5

  # ── Shared ──────────────────────────────────────────────────────────────────
  founder_sizes   = c(10),   # Founder pool size(s) to test
  n.reps          = 30,      # Independent replications
  target.pop.size = 2500,    # Working population size per cycle
  n.crosses       = 625,     # Crosses per cycle (125F × 5M)

  # ── Scheme toggle ────────────────────────────────────────────────────────────
  run_ms3  = TRUE,   # Genetic Male Sterility
  run_chem = TRUE,   # Chemical Sterilization
  run_conv = TRUE,   # Conventional Inbred Development

  # ── Step 3: Genetic Male Sterility (Scheme 1a) ──────────────────────────────
  n.cycles_ms      = 50,    # Recurrent cycles
  n_ms_select      = 125,   # Females (ms/ms) selected per cycle
  n_fertile_select = 250,   # Males (Ms/*) selected per cycle
  males_per_female = 5,     # Each female crosses with 5 males (with replacement)
  tracked_trait_id = 2,     # Trait index for ms3 locus
  max_attempts     = 10,    # Max resampling attempts to find an ms/ms founder

  # ── Step 3: Conventional Inbred Development (Scheme 2) ──────────────────────
  n.cycles_exe_rm      = 10,   # Recurrent cycles
  n.gens               = 5,    # Generations per cycle (F1 → F5)
  select_top_fams_conv = 25,   # Families selected in Step 2 selfing phase

  # Individuals / families kept at each generation within a cycle
  n_select_vector = c(
    NA,   # F0 — not used
    NA,   # F1 — not used
    2,    # F2 → F3: 2 best plants per family  (250 × 10 = 2 500)
    125,  # F3 → F4: best 125 families          (125 × 20 = 2 500)
    50,   # F4 → F5: top 50 families             (50 × 50 = 2 500)
    25    # F5:       1 line from 25 best families
  ),

  # Offspring produced per family at each generation
  n.progeny_vector = c(
    1,    # F0 → F1:  1 progeny/cross  → 125 F1 families
    20,   # F1 → F2:  20 plants/family → 125 × 20 = 2 500
    10,   # F2 → F3:  10 plants/family → 250 × 10 = 2 500
    20,   # F3 → F4:  20 plants/family → 125 × 20 = 2 500
    50,   # F4 → F5:  50 plants/family →  50 × 50 = 2 500
    NA    # F5:       selection only
  ),

  # ── Selection ────────────────────────────────────────────────────────────────
  use_selection_index     = TRUE,
  selection_index_weights = c(0.8, 0, 0.2),   # 80 % Yield, 20 % Stress response

  # ── GWAS ─────────────────────────────────────────────────────────────────────
  run_gwas         = TRUE,
  gwas_cycle       = 10,
  gwas_threshold   = 0.05,
  gwas_random_size = 1000,
  gwas_distance    = 10000,
  n.QTL.trait3     = 1,

  # ── QTL-assisted selection (cycles 11+) ──────────────────────────────────────
  use_qtl_filtering_t3       = TRUE,
  use_qtl_filtering_founders = FALSE,
  n.top.qtl                  = 1,

  # ── Diversity / Ne estimation ─────────────────────────────────────────────────
  estimate_ne         = TRUE,
  prop.markers        = 0.3,
  n.marker.sample     = 5,
  n.replicates.bgld   = 30,
  maf_threshold_value = 0.2
)


initialize_globals <- function() {
  list2env(params, envir = .GlobalEnv)
  base_pop_size <<- Ne * 12   # 300 inbred lines for the founder pool
}
initialize_globals()


# ── Logging ───────────────────────────────────────────────────────────────────
log_timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
log_filename  <- paste0("overall_scenarios_", log_timestamp, ".txt")
sink(log_filename, split = TRUE, append = TRUE, type = "output")
start.time <- Sys.time()


# ============================================================
# STEP 1: CREATE FOUNDER POOL
# RunMacs: 300 inbred lines, 10 chr, 300 seg. sites/chr, Ne = 25
# ============================================================

create_base_population <- function(n.chr, n.sites, base_pop_size) {
  # Manuscript/blueprint: 300 diverse individuals, 10 chr, 1×10^8 bp/chr,
  #   Ne = 25, mu = 2.5×10^-8/bp, r = 1×10^-8/bp (1 cM/Mb → 1 Morgan/chr).
  # theta = 4·Ne·mu = 4·25·2.5e-8 = 2.5e-6 per bp
  # rho   = 4·Ne·r  = 4·25·1e-8   = 1e-6   per bp
  base_pop <- runMacs(
    nInd          = base_pop_size,
    nChr          = n.chr,
    segSites      = n.sites,
    inbred        = FALSE,
    manualCommand = paste(
      "100000000",                    # 1×10^8 bp per chromosome
      "-t", 4 * Ne * 2.5e-8,          # θ per bp
      "-r", 4 * Ne * 1e-8,             # ρ per bp
      "-eN", 10 / (4 * Ne), 100 / Ne  # Ne×4 ten generations ago
    ),
    manualGenLen  = 1,
    nThreads      = NULL
  )

  SP_local <- SimParam$new(base_pop)
  SP_local$addTraitADEG(nQtlPerChr = 100,
                        mean = 0, var = 5, meanDD = 0.3, relAA = 0.3, varGxE = 2)
  SP_local$addTraitA(nQtlPerChr = c(0,0,0,0,0,0,0,0,0,1), mean = 0, var = 1)
  SP_local$addTraitAD(nQtlPerChr = 1, mean = 0, var = 1, meanDD = 0.5)
  SP_local$setVarE(h2 = c(0.5, 0.5, 0.5))
  SP_local$addSnpChip(nSnpPerChr = 0.1 * n.sites)

  list(base_pop = base_pop, SP = SP_local)
}

sample_founders <- function(base_pop_data, founder_size,
                            use_qtl_filtering = FALSE, n.top.qtl = 1) {
  base_pop <- base_pop_data$base_pop
  SP       <- base_pop_data$SP

  if (use_qtl_filtering_founders) {
    qtl_effects   <- SP$traits[[2]]@addEff
    top_idx       <- order(abs(qtl_effects), decreasing = TRUE)[seq_len(n.top.qtl)]
    fav_allele    <- ifelse(qtl_effects[top_idx] > 0, 2L, 0L)
    qtl_geno_base <- pullQtlGeno(base_pop, trait = 2)
    fav_count     <- rowSums(sweep(qtl_geno_base[, top_idx, drop = FALSE], 2, fav_allele, `==`))
    founders      <- base_pop[order(fav_count, decreasing = TRUE)[seq_len(founder_size)]]
  } else {
    founders <- base_pop[sample(seq_len(base_pop@nInd), founder_size)]
  }

  initial_parents <- newPop(founders, simParam = SP)
  rm(founders); gc()
  list(population = initial_parents, SP = SP)
}


# ============================================================
# STEP 3: RUN ALL SCENARIOS
# ============================================================

run_all_scenarios <- function() {
  plan(multisession, workers = 8)
  handlers(global = TRUE)
  handlers("rstudio")

  all_results      <- list()
  all_gwas_results <- list()

  cat("\n=== Creating base populations ===\n")
  base_populations <- future_lapply(seq_len(params$n.reps), function(rep) {
    create_base_population(params$n.chr, params$n.sites, base_pop_size)
  }, future.seed = TRUE)


  assemble <- function(results_list, scheme_label, size) {
    breeding  <- Filter(Negate(is.null), lapply(results_list, `[[`, "breeding_results"))
    df        <- do.call(rbind, breeding)
    df$scheme <- scheme_label
    gwas_list <- lapply(seq_along(results_list), function(i) {
      g <- results_list[[i]]$gwas_results
      if (is.null(g)) return(NULL)
      g$scheme <- scheme_label; g$founder_size <- size; g
    })
    list(breeding = df, gwas = Filter(Negate(is.null), gwas_list))
  }

  for (size in params$founder_sizes) {
    cat("\n\n=== Founder size:", size, "===\n")
    scheme_breeding <- list()

    ## ── Scheme 1a: Genetic Male Sterility ───────────────────────────────────
    if (run_ms3) {
      cat("\n=== SCHEME 1a: Genetic Male Sterility ===\n")
      source("male_sterility_scheme_gwas.R")
      ms_results <- future_lapply(seq_len(params$n.reps), function(rep) {
        fd <- sample_founders(base_populations[[rep]], size,
                              use_qtl_filtering = params$use_qtl_filtering_founders,
                              n.top.qtl         = params$n.top.qtl)
        assign("SP", fd$SP, envir = .GlobalEnv)
        result <- ms3_single_rep(rep, selected_founders = fd$population)
        rm("SP", envir = .GlobalEnv); result
      }, future.seed = TRUE)
      ms_a <- assemble(ms_results, "Genetic Male Sterility", size)
      all_gwas_results <- c(all_gwas_results, ms_a$gwas)
      scheme_breeding[["ms3"]] <- ms_a$breeding
    }

    ## ── Scheme 1b: Chemical Sterilization ───────────────────────────────────
    if (run_chem) {
      cat("\n=== SCHEME 1b: Chemical Sterilization ===\n")
      source("chem_sterilization_scheme_gwas.R")
      chem_results <- future_lapply(seq_len(params$n.reps), function(rep) {
        fd <- sample_founders(base_populations[[rep]], size,
                              use_qtl_filtering = params$use_qtl_filtering_founders,
                              n.top.qtl         = params$n.top.qtl)
        assign("SP", fd$SP, envir = .GlobalEnv)
        result <- chem_single_rep(rep, selected_founders = fd$population)
        rm("SP", envir = .GlobalEnv); result
      }, future.seed = TRUE)
      chem_a <- assemble(chem_results, "Chemical Sterilization", size)
      all_gwas_results <- c(all_gwas_results, chem_a$gwas)
      scheme_breeding[["chem"]] <- chem_a$breeding
    }

    ## ── Scheme 2: Conventional Inbred Development ───────────────────────────
    if (run_conv) {
      cat("\n=== SCHEME 2: Conventional Inbred Development ===\n")
      source("conventional_inbred_scheme_gwas.R")
      conv_results <- future_lapply(seq_len(params$n.reps), function(rep) {
        fd <- sample_founders(base_populations[[rep]], size,
                              use_qtl_filtering = params$use_qtl_filtering_founders,
                              n.top.qtl         = params$n.top.qtl)
        assign("SP", fd$SP, envir = .GlobalEnv)
        result <- conv_single_rep(rep, selected_founders = fd$population)
        rm("SP", envir = .GlobalEnv); result
      }, future.seed = TRUE)
      conv_a <- assemble(conv_results, "Conventional Inbred", size)
      all_gwas_results <- c(all_gwas_results, conv_a$gwas)
      scheme_breeding[["conv"]] <- conv_a$breeding
    }

    ## ── Combine across active schemes ────────────────────────────────────────
    common_cols <- Reduce(intersect, lapply(scheme_breeding, colnames))
    combined    <- do.call(rbind, lapply(scheme_breeding, function(df) df[, common_cols]))
    combined$founder_size <- size
    all_results[[as.character(size)]] <- combined
  }


  ## ── Save outputs ──────────────────────────────────────────────────────────
  final_results      <- do.call(rbind, all_results)
  final_results_file <<- paste0("all_scenarios_results_", log_timestamp, ".csv")
  write.csv(final_results, final_results_file, row.names = FALSE)

  if (length(all_gwas_results) > 0) {
    final_gwas        <- do.call(rbind, all_gwas_results)
    gwas_results_file <<- paste0("all_gwas_results_", log_timestamp, ".csv")
    write.csv(final_gwas, gwas_results_file, row.names = FALSE)
    cat("\n=== GWAS results saved to:", gwas_results_file, "===\n")
  }

  final_results
}


# ── Execute ───────────────────────────────────────────────────────────────────
final_data <- run_all_scenarios()
source("plot_individual_plots_true_mean.R")
plot_all_schemes_fair_with_lm(final_results_file)

end.time <- Sys.time()
cat("\nStarted at", format(start.time, "%X"), "\n")
cat("Ended at",   format(end.time,   "%X"), "\n")
cat("Duration:",  format(end.time - start.time), "\n\n")
sink()
cat("Log saved to:", log_filename, "\n")
