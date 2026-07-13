# ============================================================
# GWAS
# ============================================================

sim.gwas <- function(SNP, n.ind, random.size, p, distance, pheno, SNP.Map, QTLmap) {

  if (random.size < n.ind) {
    random.indices <- sample(1:n.ind, random.size, replace = FALSE)
    SNP   <- SNP[random.indices, ]
    pheno <- pheno[random.indices, ]
  }

  geno             <- SNP
  colnames(geno)   <- paste("X", colnames(geno), sep = "")
  rownames(SNP.Map) <- SNP.Map$SNP.names
  SNP.Map          <- SNP.Map[, c("chr", "pos")]

  common_markers <- intersect(colnames(geno), rownames(SNP.Map))
  geno    <- geno[, common_markers]
  SNP.Map <- SNP.Map[common_markers, ]

  gData      <- createGData(geno = geno, map = SNP.Map, pheno = pheno)
  gwas_result <- runSingleTraitGwas(gData      = gData,
                                    traits     = "Trait3",
                                    GLSMethod  = "single",
                                    thrType    = "bonf",
                                    alpha      = 0.05)

  sig_snps  <- gwas_result$GWAResult[[1]]
  hit_count <- 0

  if (is.data.frame(sig_snps) && nrow(sig_snps) > 0) {
    sig_snps_filtered <- sig_snps[sig_snps$pValue < p & !is.na(sig_snps$pValue), ]
    cat(" >> Significant SNPs found:", nrow(sig_snps_filtered), "\n")

    if (nrow(sig_snps_filtered) > 0) {
      for (i in seq_len(nrow(sig_snps_filtered))) {
        qtl_same_chr <- QTLmap[QTLmap$chr == sig_snps_filtered$chr[i], ]
        if (nrow(qtl_same_chr) > 0 &&
            any(abs(qtl_same_chr$pos - sig_snps_filtered$pos[i]) <= distance))
          hit_count <- hit_count + 1
      }
    }
  }

  cat(" >> Hit count:", hit_count, "\n")
  return(hit_count)
}


# ============================================================
# SELECTION
# ============================================================

select_with_method <- function(population, n_to_select, selection_trait,
                               use_index = FALSE, index_weights = NULL,
                               family_selection = FALSE, n_families = NULL,
                               simParam = SP) {

  if (population@nInd == 0) return(population)

  if (family_selection) {
    if (use_index)
      return(selectFam(population, nFam = n_families, trait = selIndex,
                       selectTop = TRUE, simParam = simParam, b = index_weights))
    return(selectFam(population, nFam = n_families, trait = selection_trait,
                     selectTop = TRUE, simParam = simParam))
  }

  if (use_index)
    return(selectInd(population, nInd = n_to_select, trait = selIndex,
                     selectTop = TRUE, simParam = simParam, b = index_weights))
  return(selectInd(population, trait = selection_trait, nInd = n_to_select,
                   selectTop = TRUE, simParam = simParam))
}


# ============================================================
# ms3 / CHEMICAL STERILIZATION SHARED HELPERS
# ============================================================

# 125 females × 5 males (with replacement) × 4 offspring = 2 500 S0
generate_structured_crosses <- function(females, males, simParam,
                                        targetSize       = 2500,
                                        malesPerFemale   = 5,
                                        nProgenyPerCross = 4) {
  cross_plan <- do.call(rbind, lapply(seq_len(females@nInd), function(i) {
    cbind(rep(i, malesPerFemale),
          sample(seq_len(males@nInd), malesPerFemale, replace = TRUE) + females@nInd)
  }))
  cross_pop <- mergePops(list(females, males))
  offspring <- makeCross(cross_pop, crossPlan = cross_plan,
                         nProgeny = nProgenyPerCross, simParam = simParam)
  if (offspring@nInd > targetSize)
    offspring <- offspring[sample(seq_len(offspring@nInd), targetSize)]
  offspring
}

# ms3 locus allele frequency summary
get_ms3_allele_freq <- function(pop) {
  dosage             <- pullQtlGeno(pop, trait = tracked_trait_id)[, 1]
  normal_allele_freq <- mean(dosage) / 2
  ms3_allele_freq    <- 1 - normal_allele_freq
  list(ms3_allele_freq    = ms3_allele_freq,
       normal_allele_freq = normal_allele_freq,
       maf  = min(normal_allele_freq, ms3_allele_freq),
       n_00 = sum(dosage == 0),
       n_01 = sum(dosage == 1),
       n_11 = sum(dosage == 2))
}

# Append ms3 genotype counts to a metrics data.frame
# Used by both GMS and CS schemes (neutral marker tracking in CS)
add_ms3_cols <- function(df, pop) {
  s <- get_ms3_allele_freq(pop)
  cbind(df,
        ms3_allele_freq    = round(s$ms3_allele_freq,    4),
        normal_allele_freq = round(s$normal_allele_freq, 4),
        ms3_maf            = round(s$maf,                4),
        n_male_sterile     = s$n_00,
        n_segregating      = s$n_01,
        n_fertile_hom      = s$n_11)
}


# ============================================================
# DIVERSITY METRICS
# ============================================================

calculate_msv <- function(my_pop, selection_trait) {
  f1_pedigree <- getPed(my_pop)
  f1_bvs      <- bv(my_pop)
  family_ids  <- unique(f1_pedigree[, c("mother", "father")])
  family_msvs <- numeric(nrow(family_ids))

  for (i in seq_len(nrow(family_ids))) {
    prog_idx <- which(f1_pedigree[, "mother"] == family_ids[i, "mother"] &
                      f1_pedigree[, "father"] == family_ids[i, "father"])
    if (length(prog_idx) > 1)
      family_msvs[i] <- var(f1_bvs[prog_idx, selection_trait])
  }

  average_msv <- mean(family_msvs, na.rm = TRUE)
  cat(sprintf("    Average Mendelian Sampling Variance: %f\n", average_msv))
  round(average_msv, 2)
}


calculate_pi <- function(mypop_geno, mypop) {
  p                   <- colSums(mypop_geno) / (2 * mypop@nInd)
  nucleotide_diversity <- mean(2 * p * (1 - p), na.rm = TRUE)
  cat(sprintf("    Nucleotide Diversity (Pi): %f\n", nucleotide_diversity))
  nucleotide_diversity
}


calcBackgroundLD <- function(pop, simParam,
                              n_chr            = 10,
                              n_markers_sample = 5,
                              n_replicates     = 100,
                              maf_threshold    = maf_threshold_value,
                              verbose          = TRUE) {

  if (verbose) {
    cat(">>> Starting background LD calculation...\n")
    if (!is.null(maf_threshold))
      cat(sprintf("    MAF Threshold: >= %.2f\n", maf_threshold))
    cat(sprintf("    Chromosomes: %d | Markers per chr: %d | Replicates: %d\n",
                n_chr, n_markers_sample, n_replicates))
  }

  geno_matrix <- pullSegSiteGeno(pop, simParam = simParam)
  marker_map  <- getGenMap(object = simParam)
  sampling_indices_list <- split(seq_len(nrow(marker_map)), f = marker_map$chr)

  if (!is.null(maf_threshold)) {
    p_allele <- colSums(geno_matrix) / (2 * pop@nInd)
    maf      <- pmin(p_allele, 1 - p_allele)
    passing  <- which(maf >= maf_threshold)
    sampling_indices_list <- lapply(sampling_indices_list,
                                    function(x) intersect(x, passing))
    if (verbose)
      cat(sprintf("    Applied MAF filter: %d of %d markers remain.\n",
                  length(unlist(sampling_indices_list)), nrow(marker_map)))
  }

  replicate_avg_r2 <- numeric(n_replicates)

  for (rep in seq_len(n_replicates)) {
    all_pairwise_r2 <- c()

    for (chr_i in seq_len(n_chr)) {
      src_idx <- sampling_indices_list[[chr_i]]
      if (length(src_idx) < n_markers_sample) next
      sampled_src <- sample(src_idx, n_markers_sample)

      for (chr_j in seq_len(n_chr)[-chr_i]) {
        cmp_idx <- sampling_indices_list[[chr_j]]
        if (length(cmp_idx) < n_markers_sample) next
        sampled_cmp <- sample(cmp_idx, n_markers_sample)

        r2_matrix       <- suppressWarnings(
          cor(geno_matrix[, sampled_src], geno_matrix[, sampled_cmp])^2)
        all_pairwise_r2 <- c(all_pairwise_r2, as.vector(r2_matrix))
      }
    }

    replicate_avg_r2[rep] <- if (length(all_pairwise_r2) > 0)
      mean(all_pairwise_r2, na.rm = TRUE) else NA
  }

  background_ld <- mean(replicate_avg_r2, na.rm = TRUE)

  if (verbose) {
    cat(">>> Background LD calculation complete!\n")
    cat(sprintf("    Final average LD (r^2): %.4f\n", background_ld))
  }

  list(mean_bg_ld      = round(background_ld, 4),
       replicate_bg_ld = replicate_avg_r2)
}


# LD-based Ne estimator (mirrors gl.LDNe: singleton removal, inter-chromosomal pairs)
estimate_Ne_like_glLDNe <- function(geno_matrix) {
  n_ind  <- nrow(geno_matrix)
  snps   <- ncol(geno_matrix)

  chrom_vector  <- factor(rep(paste0("Chr", seq_len(n.chr)), length.out = snps))
  allele_counts <- colSums(geno_matrix, na.rm = TRUE)
  non_singletons <- which(allele_counts > 1 & allele_counts < (2 * n_ind - 1))

  geno_matrix  <- geno_matrix[, non_singletons, drop = FALSE]
  chrom_vector <- chrom_vector[non_singletons]

  pair_idx    <- combn(seq_len(ncol(geno_matrix)), 2)
  valid_pairs <- which(chrom_vector[pair_idx[1, ]] != chrom_vector[pair_idx[2, ]])
  pair_idx    <- pair_idx[, valid_pairs]

  r2_vals <- numeric(ncol(pair_idx))
  for (i in seq_len(ncol(pair_idx))) {
    snp1 <- geno_matrix[, pair_idx[1, i]]
    snp2 <- geno_matrix[, pair_idx[2, i]]
    if (var(snp1, na.rm = TRUE) > 0 && var(snp2, na.rm = TRUE) > 0)
      r2_vals[i] <- cor(snp1, snp2, use = "pairwise.complete.obs")^2
  }

  r2_vals  <- r2_vals[!is.na(r2_vals)]
  mean_r2  <- mean(r2_vals)
  Ne_val   <- 1 / (3 * (mean_r2 - (1 / (2 * n_ind))))

  list(Ne = Ne_val, r2 = mean_r2)
}


# ============================================================
# SHARED HELPERS — used by all scheme simulation files
# ============================================================

# Collect standard per-population metrics into a one-row data.frame.
collect_pop_metrics <- function(pop, recurrent_cycle, generation, SP) {
  sel_trait <- if (recurrent_cycle <= 10) 1L else 3L

  qtl_geno     <- pullQtlGeno(pop, trait = 3)
  allele_freqs <- colMeans(qtl_geno, na.rm = TRUE) / 2
  maf_t3       <- pmin(allele_freqs, 1 - allele_freqs)

  geno <- pullSegSiteGeno(pop)
  if (estimate_ne && nrow(geno) >= 2) {
    sampled_cols <- sample(seq_len(ncol(geno)), floor(ncol(geno) * prop.markers))
    est    <- estimate_Ne_like_glLDNe(geno[, sampled_cols])
    ne_val <- round(est$Ne, 4)
    r2_val <- round(est$r2, 4)
  } else {
    ne_val <- r2_val <- NA
  }

  data.frame(
    recurrent_cycle        = recurrent_cycle,
    generation             = generation,
    n_ind                  = pop@nInd,
    pheno_mean             = round(meanP(pop)[[1]], 4),
    pheno_var              = round(varP(pop)[1],    4),
    geno_mean              = round(meanG(pop)[[1]], 4),
    geno_var               = round(varG(pop)[1],    4),
    Ne_like_glLDNe         = ne_val,
    r2_like_glLDNe         = r2_val,
    avg_maf                = mean(maf_t3, na.rm = TRUE),
    n_fixed_qtl            = sum(maf_t3 == 0),
    nucleotide_diversity   = calculate_pi(geno, pop),
    mendelian_sampling_var = calculate_msv(pop, sel_trait),
    background_ld          = calcBackgroundLD(pop, SP,
                               n_chr            = n.chr,
                               n_markers_sample = n.marker.sample,
                               n_replicates     = n.replicates.bgld)$mean_bg_ld
  )
}


# Run GWAS for Trait 3 on the supplied population.
run_scheme_gwas <- function(pop, rep_index, recurrent_cycle, generation, SP) {
  QTLmap  <- getQtlMap(trait = 3, simParam = SP)
  QTLgeno <- pullQtlGeno(pop, trait = 3, simParam = SP)
  ns      <- nrow(QTLgeno)
  maf_qtl <- apply(rbind(0.5 * colSums(QTLgeno) / ns,
                         1 - 0.5 * colSums(QTLgeno) / ns), 2, min)
  qtl_loss_rate <- sum(maf_qtl == 0) / (n.QTL.trait3 * n.chr)

  SNP     <- pullSnpGeno(pop, snpChip = 1, simParam = SP)
  pheno_d <- data.frame(genotype = rownames(SNP), pheno(pop))

  SNP_Map           <- getSnpMap(snpChip = 1, simParam = SP)[, -4]
  colnames(SNP_Map) <- c("SNP.names", "chr", "pos")
  SNP_Map$allele1   <- "A"
  SNP_Map$allele2   <- "T"
  SNP_Map$chr       <- as.integer(SNP_Map$chr)
  SNP_Map$SNP.names <- paste0("X", SNP_Map$SNP.names)

  n_hits <- sim.gwas(SNP = SNP, n.ind = pop@nInd, random.size = gwas_random_size,
                     p = gwas_threshold, distance = gwas_distance,
                     pheno = pheno_d, SNP.Map = SNP_Map, QTLmap = QTLmap)

  rm(QTLgeno, SNP, pheno_d, SNP_Map); gc()

  data.frame(
    rep             = rep_index,
    recurrent_cycle = recurrent_cycle,
    generation      = generation,
    qtl_loss_rate   = qtl_loss_rate,
    hit_rate        = n_hits / (n.QTL.trait3 * n.chr),
    n_qtl_detected  = n_hits,
    n_ind           = pop@nInd
  )
}


# Retain only individuals homozygous for the favorable allele at the top n.top.qtl
# QTLs for Trait 3. Falls back to top carriers if no fully-homozygous individuals exist.
apply_qtl_filter <- function(pop, SP) {
  qtl_effects <- SP$traits[[3]]@addEff
  top_idx     <- order(abs(qtl_effects), decreasing = TRUE)[seq_len(n.top.qtl)]
  fav_allele  <- ifelse(qtl_effects[top_idx] > 0, 1L, 0L)

  qtl_geno  <- pullQtlGeno(pop, trait = 3)[, top_idx, drop = FALSE]
  fav_count <- matrix(0L, nrow(qtl_geno), ncol(qtl_geno))
  for (i in seq_len(ncol(qtl_geno)))
    fav_count[, i] <- if (fav_allele[i] == 1L) qtl_geno[, i] else 2L - qtl_geno[, i]

  is_homozygous_fav <- rowSums(fav_count == 2L) == n.top.qtl
  if (any(is_homozygous_fav)) return(pop[is_homozygous_fav])

  warning("No individuals homozygous for favorable allele at all QTL; keeping top carriers.")
  top_carriers <- order(rowSums(fav_count >= 1L), decreasing = TRUE)[seq_len(min(100L, pop@nInd))]
  pop[top_carriers]
}
