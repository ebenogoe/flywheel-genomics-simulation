# ============================================================
# SCHEME 1b: RAPID CYCLING — CHEMICAL STERILIZATION
# Females are randomly assigned from the phenotypically best
# 375 individuals (top n_ms_select + n_fertile_select).
# Chemical sterilisation is applied randomly within the
# selected pool — not preferentially to the top-ranked plants.
# Same crossing geometry as Scheme 1a: 125F × 5M × 4 offspring = 2 500 S0.
# ============================================================

message("Loading Chemical Sterilization scheme functions...")


# ============================================================
chem_single_rep <- function(rep_index, selected_founders) {

  all_results      <- list()
  all_gwas_results <- list()

  # ── Cycle 0: parent baseline ────────────────────────────────────────────────
  parent_results <- add_ms3_cols(
    collect_pop_metrics(selected_founders, recurrent_cycle = 0, generation = 0, SP),
    selected_founders
  )


  # ── Step 2: Create founding S0 ──────────────────────────────────────────────
  # CS founding: all founders randomly intermated (no ms introgression needed).
  # 625 random crosses × 4 offspring = 2 500 founding S0.
  current_S0 <- randCross(selected_founders, nCrosses = n.crosses, nProgeny = 4, simParam = SP)
  cat(sprintf("Step 2 complete: founding S0 = %d individuals\n", current_S0@nInd))
  rm(selected_founders); gc()


  # ── Step 3: Recurrent selection cycles ──────────────────────────────────────
  for (recurrent_cycle in seq_len(n.cycles_ms)) {
    cat(sprintf(">>> Rep %d  Cycle %d\n", rep_index, recurrent_cycle))

    current_selection_trait <- if (recurrent_cycle <= 10) 1L else 3L
    use_index_this_cycle    <- use_selection_index && recurrent_cycle > 10

    # GWAS for Trait 3 on full S0 at cycle gwas_cycle (before any selection filter)
    if (run_gwas && recurrent_cycle == gwas_cycle) {
      all_gwas_results[[1]] <- run_scheme_gwas(current_S0, rep_index,
                                                recurrent_cycle, generation = 1L, SP)
    }

    # Collect metrics on the full S0 (before QTL filter, reflects true pop state)
    cycle_results <- add_ms3_cols(
      collect_pop_metrics(current_S0, recurrent_cycle, generation = 1L, SP),
      current_S0
    )
    all_results[[recurrent_cycle]] <- list(parent_results, cycle_results)

    # QTL-assisted pre-selection: retain only homozygous favorable individuals (cycles 11+)
    if (recurrent_cycle > 10 && use_qtl_filtering_t3)
      current_S0 <- apply_qtl_filter(current_S0, SP)

    # Parental selection: select the best n_ms_select + n_fertile_select individuals,
    # then randomly assign n_ms_select as females (chemically sterilised) and the
    # remainder as males. Manuscript: "125 females and 250 males were assigned randomly
    # from a pool of the best 375 individuals."
    n_total   <- min(n_ms_select + n_fertile_select, current_S0@nInd)
    top_pool  <- select_with_method(current_S0, n_total,
                                    current_selection_trait,
                                    use_index     = use_index_this_cycle,
                                    index_weights = selection_index_weights,
                                    simParam      = SP)

    female_idx       <- sample(seq_len(top_pool@nInd), min(n_ms_select, top_pool@nInd))
    male_idx         <- setdiff(seq_len(top_pool@nInd), female_idx)
    selected_females <- top_pool[female_idx]
    selected_males   <- top_pool[male_idx]
    rm(top_pool); gc()

    cat(sprintf("    Selected: %d females (chem.) + %d males\n",
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
  parent_df <- all_results[[1]][[1]]
  cycle_df  <- do.call(rbind, lapply(all_results, `[[`, 2))
  final_df  <- cbind(rep = rep_index, rbind(parent_df, cycle_df))

  gwas_df <- if (run_gwas && length(all_gwas_results) > 0)
    do.call(rbind, all_gwas_results) else NULL

  list(breeding_results = final_df, gwas_results = gwas_df)
}
