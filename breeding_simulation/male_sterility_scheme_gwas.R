# ============================================================
# SCHEME 1a: RAPID CYCLING — GENETIC MALE STERILITY
# ms/ms  → male sterile (females)
# Ms/ms or Ms/Ms → fertile (males)
# ============================================================

identify_male_sterile <- function(pop) {
  pullQtlGeno(pop, trait = tracked_trait_id)[, 1] == 0
}

identify_fertiles <- function(pop) {
  pullQtlGeno(pop, trait = tracked_trait_id)[, 1] >= 1
}

message("Loading Genetic Male Sterility scheme functions...")


# ============================================================
ms3_single_rep <- function(rep_index, selected_founders) {

  all_results      <- list()
  all_gwas_results <- list()
  founder_size     <- selected_founders@nInd


  # ── Step 2: Create founding S0 ──────────────────────────────────────────────
  # Blueprint: "All founders crossed to male sterile donor"
  # Goal: produce F1s carrying the ms allele, then randCross → 2 500 founding S0.
  #
  # Path A (inbred=FALSE): natural Ms/ms het (dosage=1) exists → use as donor.
  # Path B (inbred=TRUE):  no natural hets → cross ms/ms (dosage=0) × non-ms
  #                        to produce all-Ms/ms F1 donors first.
  #
  # Resample if NEITHER path is viable (no hets AND no ms/ms, or no non-ms).
  ms_dosage   <- pullQtlGeno(selected_founders, trait = 2)[, 1]
  het_indices <- which(ms_dosage == 1)   # Ms/ms heterozygotes
  msms_idx    <- which(ms_dosage == 0)   # ms/ms homozygous
  nonms_idx   <- which(ms_dosage > 0)    # Ms/ms or Ms/Ms

  can_proceed <- length(het_indices) > 0 ||
                 (length(msms_idx) > 0 && length(nonms_idx) > 0)

  attempt <- 1
  while (!can_proceed && attempt <= max_attempts) {
    cat(sprintf("Attempt %d: cannot form Ms/ms donors — resampling...\n", attempt))
    new_base          <- create_base_population(params$n.chr, params$n.sites, base_pop_size)
    fd                <- sample_founders(new_base, founder_size,
                                         use_qtl_filtering = params$use_qtl_filtering_founders,
                                         n.top.qtl         = params$n.top.qtl)
    selected_founders <- fd$population
    SP                <<- fd$SP
    ms_dosage         <- pullQtlGeno(selected_founders, trait = 2, simParam = SP)[, 1]
    het_indices       <- which(ms_dosage == 1)
    msms_idx          <- which(ms_dosage == 0)
    nonms_idx         <- which(ms_dosage > 0)
    can_proceed       <- length(het_indices) > 0 ||
                         (length(msms_idx) > 0 && length(nonms_idx) > 0)
    attempt           <- attempt + 1
  }
  if (!can_proceed) stop("Cannot create Ms/ms donors after max_attempts resamples.")

  # ── Cycle 0: parent baseline (collected after donor validation) ──────────────
  parent_results <- add_ms3_cols(
    collect_pop_metrics(selected_founders, recurrent_cycle = 0, generation = 0, SP),
    selected_founders
  )

  # Build crossing plan to create Ms/ms F1 donors
  if (length(het_indices) > 0) {
    # Path A: use existing Ms/ms het as donor
    cat(sprintf("Path A (het donor): Ms/ms founder at index %d.\n", het_indices[1]))
    other_idx     <- setdiff(seq_len(selected_founders@nInd), het_indices[1])
    cross_plan_s2 <- cbind(rep(het_indices[1], length(other_idx)), other_idx)
  } else {
    # Path B: ms/ms × each non-ms → all F1 are Ms/ms
    cat(sprintf("Path B (inbred donors): ms/ms(%d) × %d non-ms founders.\n",
                msms_idx[1], length(nonms_idx)))
    cross_plan_s2 <- cbind(rep(msms_idx[1], length(nonms_idx)), nonms_idx)
  }

  f1_from_donor <- makeCross(selected_founders, crossPlan = cross_plan_s2,
                              nProgeny = 50, simParam = SP)

  # Randomly intermate F1 offspring: 625 crosses × 4 progeny = 2 500 founding S0
  current_S0 <- randCross(f1_from_donor, nCrosses = n.crosses, nProgeny = 4, simParam = SP)
  cat(sprintf("Step 2 complete: founding S0 = %d individuals\n", current_S0@nInd))
  rm(f1_from_donor); gc()


  # ── Step 3: Recurrent selection cycles ──────────────────────────────────────
  for (recurrent_cycle in seq_len(n.cycles_ms)) {
    cat(sprintf(">>> Rep %d  Cycle %d\n", rep_index, recurrent_cycle))

    current_selection_trait <- if (recurrent_cycle <= 10) 1L else 3L

    # GWAS for Trait 3 on full S0 at cycle gwas_cycle (before any selection filter)
    if (run_gwas && recurrent_cycle == gwas_cycle) {
      all_gwas_results[[1]] <- run_scheme_gwas(current_S0, rep_index,
                                                recurrent_cycle, generation = 1L, SP)
    }

    # Collect metrics on the full S0 (before QTL filter, so they reflect true pop state)
    cycle_results <- add_ms3_cols(
      collect_pop_metrics(current_S0, recurrent_cycle, generation = 1L, SP),
      current_S0
    )
    all_results[[recurrent_cycle]] <- list(parent_results, cycle_results)

    # QTL-assisted pre-selection: retain only homozygous favorable individuals (cycles 11+)
    if (recurrent_cycle > 10 && use_qtl_filtering_t3)
      current_S0 <- apply_qtl_filter(current_S0, SP)

    # Parental selection: 125 ms/ms females + 250 Ms/* males
    is_female <- identify_male_sterile(current_S0)
    is_male   <- identify_fertiles(current_S0)
    females_pool <- current_S0[is_female]
    males_pool   <- current_S0[is_male]

    if (females_pool@nInd == 0 || males_pool@nInd == 0) {
      cat("    Insufficient ms/ms or fertile individuals — stopping.\n"); stop()
    }

    # Cap at pool size so selectInd is never asked for more than available
    n_select_f <- min(n_ms_select,      females_pool@nInd)
    n_select_m <- min(n_fertile_select, males_pool@nInd)

    use_index_this_cycle <- use_selection_index && recurrent_cycle > 10
    selected_females <- select_with_method(females_pool, n_select_f,
                                           current_selection_trait,
                                           use_index     = use_index_this_cycle,
                                           index_weights = selection_index_weights, simParam = SP)
    selected_males   <- select_with_method(males_pool, n_select_m,
                                           current_selection_trait,
                                           use_index     = use_index_this_cycle,
                                           index_weights = selection_index_weights, simParam = SP)
    rm(females_pool, males_pool); gc()

    cat(sprintf("    Selected: %d females (ms/ms) + %d males (Ms/*)\n",
                selected_females@nInd, selected_males@nInd))

    # Crossing: 125F × 5M × 4 offspring = 2 500 S0
    current_S0 <- generate_structured_crosses(females          = selected_females,
                                               males            = selected_males,
                                               simParam         = SP,
                                               targetSize       = target.pop.size,
                                               malesPerFemale   = males_per_female,
                                               nProgenyPerCross = 4)
    rm(selected_females, selected_males); gc()
    cat(sprintf("    New S0: %d individuals\n", current_S0@nInd))
  }

  rm(current_S0); gc()

  # ── Combine and return ───────────────────────────────────────────────────────
  parent_df  <- all_results[[1]][[1]]
  cycle_df   <- do.call(rbind, lapply(all_results, `[[`, 2))
  final_df   <- cbind(rep = rep_index, rbind(parent_df, cycle_df))

  gwas_df <- if (run_gwas && length(all_gwas_results) > 0)
    do.call(rbind, all_gwas_results) else NULL

  list(breeding_results = final_df, gwas_results = gwas_df)
}
