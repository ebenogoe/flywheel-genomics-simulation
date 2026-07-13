# ============================================================
# SCHEME 2: CONVENTIONAL PURE LINE DEVELOPMENT
# One recurrent cycle = 5 years (F1 → F5)
# ============================================================

message("Loading Conventional Inbred scheme functions...")

conv_single_rep <- function(rep_index, selected_founders) {

  all_results      <- list()
  all_gwas_results <- list()

  # ── Cycle 0: parent baseline ────────────────────────────────────────────────
  parent_results <- collect_pop_metrics(selected_founders, recurrent_cycle = 0,
                                        generation = 0, SP)


  # ── Step 2: Create F5 founder lines ─────────────────────────────────────────
  # 625 random crosses × 4 offspring = 2 500 individuals
  # Phenotypically select 25 best families, 1 line/family → 25 F1 lines
  # Self F1 → F5 (4 rounds of selfing) → 25 inbred founder lines
  initial_cross <- randCross(selected_founders, nCrosses = n.crosses, nProgeny = 4, simParam = SP)
  top_families  <- selectFam(initial_cross, nFam = select_top_fams_conv,
                              use = "pheno", selectTop = TRUE, simParam = SP)
  founder_lines <- selectWithinFam(top_families, nInd = 1, use = "pheno",
                                   selectTop = TRUE, simParam = SP)
  rm(initial_cross, top_families); gc()

  f5_lines <- founder_lines
  for (gen in seq_len(4)) {
    self_plan <- cbind(seq_len(f5_lines@nInd), seq_len(f5_lines@nInd))
    f5_lines  <- makeCross(f5_lines, crossPlan = self_plan, nProgeny = 1, simParam = SP)
  }
  rm(founder_lines, selected_founders); gc()
  cat(sprintf("Step 2 complete: %d F5 founder lines created\n", f5_lines@nInd))

  # f5_lines is the starting parent pool for recurrent cycle 1
  current_parents <- f5_lines
  rm(f5_lines); gc()


  # ── Step 3: Recurrent selection cycles ──────────────────────────────────────
  for (recurrent_cycle in seq_len(n.cycles_exe_rm)) {
    cat(sprintf(">>> Rep %d  Cycle %d\n", rep_index, recurrent_cycle))

    current_selection_trait <- if (recurrent_cycle <= 10) 1L else 3L
    use_index_this_cycle    <- use_selection_index && recurrent_cycle > 10

    # ── Year 1: 125 biparental crosses ────────────────────────────────────────
    # Each parent crosses with 5 random partners (directed, no deduplication)
    # → 25 parents × 5 mates = 125 crosses × 1 progeny = 125 F1 families
    n_parents    <- current_parents@nInd
    n_partners   <- 5L
    cross_plan   <- do.call(rbind, lapply(seq_len(n_parents), function(parent) {
      partners <- sample(setdiff(seq_len(n_parents), parent), n_partners)
      cbind(rep(parent, n_partners), partners)
    }))
    F1 <- makeCross(current_parents, crossPlan = cross_plan,
                    nProgeny = n.progeny_vector[1], simParam = SP)
    rm(cross_plan); gc()
    cat(sprintf("  Year 1: %d F1 families created\n", F1@nInd))

    # Collect metrics on full F1 (before QTL filter, so they reflect true pop state)
    cycle_results <- collect_pop_metrics(F1, recurrent_cycle, generation = 1L, SP)
    all_results[[recurrent_cycle]] <- list(parent_results, cycle_results)

    # QTL-assisted pre-selection: retain only homozygous favorable individuals (cycles 11+)
    if (recurrent_cycle > 10 && use_qtl_filtering_t3)
      F1 <- apply_qtl_filter(F1, SP)


    # ── Year 2: F1 → F2 ───────────────────────────────────────────────────────
    # Bulk-self: 20 plants per family → 125 × 20 = 2 500
    # Within-family selection: keep 2 best plants per family → 250 selected
    self_plan <- cbind(seq_len(F1@nInd), seq_len(F1@nInd))
    F2 <- makeCross(F1, crossPlan = self_plan,
                    nProgeny = n.progeny_vector[2], simParam = SP)
    cat(sprintf("  Year 2: F2 = %d individuals (%d families × %d plants)\n",
                F2@nInd, F1@nInd, n.progeny_vector[2]))

    # GWAS for Trait 3 on F2 at cycle gwas_cycle
    if (run_gwas && recurrent_cycle == gwas_cycle) {
      all_gwas_results[[1]] <- run_scheme_gwas(F2, rep_index,
                                                recurrent_cycle, generation = 2L, SP)
    }

    F2_selected <- selectWithinFam(F2, nInd = n_select_vector[3],
                                   trait = current_selection_trait,
                                   use = "pheno", selectTop = TRUE, simParam = SP)
    rm(F1, F2); gc()
    cat(sprintf("  Year 2: %d F2 plants selected (%d/family)\n",
                F2_selected@nInd, n_select_vector[3]))


    # ── Year 3: F2 → F3 headrows ──────────────────────────────────────────────
    # Self: 10 plants per row → 250 × 10 = 2 500
    # Family selection: keep best 125 families, 1 representative each
    self_plan <- cbind(seq_len(F2_selected@nInd), seq_len(F2_selected@nInd))
    F3 <- makeCross(F2_selected, crossPlan = self_plan,
                    nProgeny = n.progeny_vector[3], simParam = SP)
    rm(F2_selected); gc()
    cat(sprintf("  Year 3: F3 = %d individuals\n", F3@nInd))

    F3_fam_sel <- selectFam(F3, nFam = n_select_vector[4],
                             trait = current_selection_trait,
                             use = "pheno", selectTop = TRUE, simParam = SP)
    F3_reps    <- selectWithinFam(F3_fam_sel, nInd = 1,
                                  trait = current_selection_trait,
                                  use = "pheno", selectTop = TRUE, simParam = SP)
    rm(F3, F3_fam_sel); gc()
    cat(sprintf("  Year 3: %d families selected → 1 rep each\n", n_select_vector[4]))


    # ── Year 4: F3 → F4 ───────────────────────────────────────────────────────
    # Self: 20 plants per family → 125 × 20 = 2 500
    # Family selection: keep top 50 families, 1 representative each
    self_plan <- cbind(seq_len(F3_reps@nInd), seq_len(F3_reps@nInd))
    F4 <- makeCross(F3_reps, crossPlan = self_plan,
                    nProgeny = n.progeny_vector[4], simParam = SP)
    rm(F3_reps); gc()
    cat(sprintf("  Year 4: F4 = %d individuals\n", F4@nInd))

    F4_fam_sel <- selectFam(F4, nFam = n_select_vector[5],
                             trait = current_selection_trait,
                             use = "pheno", selectTop = TRUE, simParam = SP)
    F4_reps    <- selectWithinFam(F4_fam_sel, nInd = 1,
                                  trait = current_selection_trait,
                                  use = "pheno", selectTop = TRUE, simParam = SP)
    rm(F4, F4_fam_sel); gc()
    cat(sprintf("  Year 4: %d families selected → 1 rep each\n", n_select_vector[5]))


    # ── Year 5: F4 → F5 ───────────────────────────────────────────────────────
    # Self: 50 plants per family → 50 × 50 = 2 500
    # Select 1 line from each of the 25 best families → parents for next cycle
    self_plan <- cbind(seq_len(F4_reps@nInd), seq_len(F4_reps@nInd))
    F5 <- makeCross(F4_reps, crossPlan = self_plan,
                    nProgeny = n.progeny_vector[5], simParam = SP)
    rm(F4_reps); gc()
    cat(sprintf("  Year 5: F5 = %d individuals\n", F5@nInd))

    F5_fam_sel      <- selectFam(F5, nFam = n_select_vector[6],
                                  trait = current_selection_trait,
                                  use = "pheno", selectTop = TRUE, simParam = SP)
    current_parents <- selectWithinFam(F5_fam_sel, nInd = 1,
                                       trait = current_selection_trait,
                                       use = "pheno", selectTop = TRUE, simParam = SP)
    rm(F5, F5_fam_sel); gc()
    cat(sprintf("  Year 5: %d F5 lines selected as parents for next cycle\n",
                current_parents@nInd))
  }

  rm(current_parents); gc()

  # ── Combine and return ───────────────────────────────────────────────────────
  parent_df <- all_results[[1]][[1]]
  cycle_df  <- do.call(rbind, lapply(all_results, `[[`, 2))
  final_df  <- cbind(rep = rep_index, rbind(parent_df, cycle_df))

  gwas_df <- if (run_gwas && length(all_gwas_results) > 0)
    do.call(rbind, all_gwas_results) else NULL

  list(breeding_results = final_df, gwas_results = gwas_df)
}
