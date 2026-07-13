# GWAS False Positive Rate (FPR) Simulation Pipeline
# Null scenario: pure noise phenotypes (no QTN, no genetic signal)
# Any significant SNP detected = false positive
# 3 models × 2 populations × 100 reps

library(GAPIT)
library(dplyr)
library(readr)

# ======================
# CONFIGURATION PARAMS
# ======================
N_REPS       <- 100
MAF_THRESHOLD <- 0
N_PC         <- 3

# Output directory for results
RESULTS_DIR <- "gwas_fpr_results"
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ======================
# HELPER FUNCTIONS
# ======================
# Align kinship matrix to sample IDs
align_kinship <- function(kinship_mat, sample_ids) {
  if (is.null(rownames(kinship_mat)) || all(rownames(kinship_mat) %in% as.character(1:ncol(kinship_mat)))) {
    if (length(sample_ids) != nrow(kinship_mat)) {
      stop("Sample count mismatch: kinship (", nrow(kinship_mat),
           ") vs samples (", length(sample_ids), ")")
    }
    rownames(kinship_mat) <- sample_ids
    colnames(kinship_mat) <- sample_ids
  } else {
    common_ids <- intersect(rownames(kinship_mat), sample_ids)
    if (length(common_ids) < length(sample_ids)) {
      warning("Kinship matrix missing samples: ",
              setdiff(sample_ids, rownames(kinship_mat)))
    }
    kinship_mat <- kinship_mat[common_ids, common_ids, drop = FALSE]
  }
  
  kinship_df <- data.frame(
    ID = rownames(kinship_mat),
    as.data.frame(kinship_mat),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  return(kinship_df)
}

# Format genotype data for GAPIT (handles both populations)
format_geno_for_gapit <- function(geno_df, pop_type = c("SAP", "HBP")) {
  pop_type  <- match.arg(pop_type)
  snp_names <- rownames(geno_df)
  
  if (pop_type == "HBP") {
    snp_info <- strsplit(snp_names, "_")
    chrom    <- sapply(snp_info, `[`, 1)
    pos      <- sapply(snp_info, `[`, 2)
  } else {  # SAP: "S1_12345678"
    snp_info <- strsplit(gsub("^S", "", snp_names), "_")
    chrom    <- sapply(snp_info, `[`, 1)
    pos      <- sapply(snp_info, `[`, 2)
  }
  
  GM <- data.frame(
    SNP        = snp_names,
    Chromosome = chrom,
    Position   = as.numeric(pos),
    stringsAsFactors = FALSE
  )
  
  GD <- t(geno_df)
  GD <- data.frame(
    Taxa = rownames(GD),
    as.data.frame(GD),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(GD) <- NULL
  
  return(list(GM = GM, GD = GD))
}

# Run all 3 GWAS models for a population
run_gwas_models <- function(Y, GD, GM, pcs_df, kinship_mat, pop_name, rep_id) {
  results <- list()
  
  # GLM
  cat("\n[Rep", rep_id, "] Running GLM for", pop_name, "...\n")
  results$glm <- GAPIT(
    Y                 = Y,
    GD                = GD,
    GM                = GM,
    PCA.total         = 0,
    kinship.algorithm = NULL,
    SNP.MAF           = MAF_THRESHOLD,
    model             = "GLM",
    file.output       = FALSE,
    Geno.View.output  = FALSE,
    PCA.View.output   = FALSE
  )
  
  # GLM + PCs
  cat("[Rep", rep_id, "] Running GLM+PC for", pop_name, "...\n")
  pc_match <- pcs_df[match(Y$Taxa, pcs_df$sample), c("sample", paste0("PC", 1:N_PC))]
  colnames(pc_match)[1] <- "Taxa"
  
  results$glm_pc <- GAPIT(
    Y                 = Y,
    GD                = GD,
    GM                = GM,
    CV                = pc_match,
    PCA.total         = 0,
    kinship.algorithm = NULL,
    SNP.MAF           = MAF_THRESHOLD,
    model             = "GLM",
    file.output       = FALSE,
    Geno.View.output  = FALSE,
    PCA.View.output   = FALSE
  )
  
  
  # MLM + PCs + Kinship
  cat("[Rep", rep_id, "] Running MLM PC+K for", pop_name, "...\n")
  kinship_aligned <- align_kinship(kinship_mat, Y$Taxa)
  
  results$mlm_pk <- GAPIT(
    Y                = Y,
    GD               = GD,
    GM               = GM,
    CV               = pc_match,
    KI               = kinship_aligned,
    SNP.MAF          = MAF_THRESHOLD,
    model            = "MLM",
    file.output      = FALSE,
    Geno.View.output = FALSE,
    PCA.View.output  = FALSE
  )
  
  # MLM + Kinship only
  cat("[Rep", rep_id, "] Running MLM+K for", pop_name, "...\n")
  
  results$mlm_k <- GAPIT(
    Y = Y,
    GD = GD,
    GM = GM,
    KI = kinship_aligned,
    PCA.total = 0,
    kinship.algorithm = NULL,
    SNP.MAF = MAF_THRESHOLD,
    model = "MLM",
    file.output = FALSE,
    Geno.View.output = FALSE,
    PCA.View.output = FALSE
  )
  
  return(results)
}

# Check for false positives (Bonferroni threshold)
# In this null scenario, any significant SNP = false positive
check_false_positives <- function(gwas_result) {
  gwas_df     <- gwas_result$GWAS
  bonf_thresh <- 0.05 / nrow(gwas_df)
  n_sig       <- sum(gwas_df$P.value < bonf_thresh, na.rm = TRUE)
  
  return(data.frame(
    n_sig_snps       = n_sig,
    any_false_positive = n_sig > 0,
    bonf_thresh      = bonf_thresh,
    n_snps_tested    = nrow(gwas_df)
  ))
}

# ======================
# MAIN PIPELINE
# ======================

# Load base genotype data ONCE
cat("Loading base genotype data...\n")
load("sap_geno_df_MAF0.05.RData")
load("hbp_geno_df_chibas_MAF0.05.RData")

# Load population-specific resources
cat("Loading PCs and kinship...\n")
sap_pcs     <- readRDS("SAP/SAP_PCs.rds")
hbp_pcs     <- readRDS("HBP/HBP_PCs.rds")
sap_kinship <- readRDS("SAP/SAP_kinship.rds")
hbp_kinship <- readRDS("HBP/HBP_kinship.rds")

# Pre-format genotype data once (no QTN removal needed in null scenario)
cat("Formatting genotype data for GAPIT...\n")
sap_gapit <- format_geno_for_gapit(sap_geno_df_filtered, pop_type = "SAP")
hbp_gapit <- format_geno_for_gapit(hbp_geno_df_filtered, pop_type = "HBP")

# Get sample IDs from GD (Taxa column)
sap_sample_ids <- sap_gapit$GD$Taxa
hbp_sample_ids <- hbp_gapit$GD$Taxa

# Initialize results storage
all_results <- list()

# ======================
# NULL SIMULATION LOOP
# ======================
for (rep_id in 1:N_REPS) {
  cat("\n========================================\n")
  cat("NULL REPlicate", rep_id, "of", N_REPS, "\n")
  cat("========================================\n")
  
  # --- SAP null phenotype ---
  # Pure noise ie no QTN, no heritability. Any hit is a false positive
  cat("\n--- SAP Population (null) ---\n")
  sap_null_pheno <- data.frame(
    Taxa  = sap_sample_ids,
    Pheno = rnorm(length(sap_sample_ids)), 
    stringsAsFactors = FALSE
  )
  
  sap_gwas_null <- run_gwas_models(
    Y           = sap_null_pheno,
    GD          = sap_gapit$GD,
    GM          = sap_gapit$GM,
    pcs_df      = sap_pcs,
    kinship_mat = sap_kinship,
    pop_name    = "SAP_null",
    rep_id      = rep_id
  )
  
  sap_fp <- lapply(c("glm", "glm_pc", "mlm_k", "mlm_pk"), function(model) {
    check_false_positives(sap_gwas_null[[model]]) %>%
      mutate(model = model, population = "SAP")
  })
  
  # --- HBP null phenotype ---
  cat("\n--- HBP Population (null) ---\n")
  hbp_null_pheno <- data.frame(
    Taxa  = hbp_sample_ids,
    Pheno = rnorm(length(hbp_sample_ids)),
    stringsAsFactors = FALSE
  )
  
  hbp_gwas_null <- run_gwas_models(
    Y           = hbp_null_pheno,
    GD          = hbp_gapit$GD,
    GM          = hbp_gapit$GM,
    pcs_df      = hbp_pcs,
    kinship_mat = hbp_kinship,
    pop_name    = "HBP_null",
    rep_id      = rep_id
  )
  
  hbp_fp <- lapply(c("glm", "glm_pc", "mlm_k", "mlm_pk"), function(model) {
    check_false_positives(hbp_gwas_null[[model]]) %>%
      mutate(model = model, population = "HBP")
  })
  
  # Combine this replicate's results
  rep_df <- bind_rows(c(sap_fp, hbp_fp)) %>%
    mutate(replicate = rep_id)
  
  all_results[[rep_id]] <- rep_df
  
  # Incremental checkpoint save
  if (rep_id %% 10 == 0) {
    saveRDS(bind_rows(all_results), file.path(RESULTS_DIR, "fpr_results_incremental.RDS"))
    cat("\n[Checkpoint] Saved results up to replicate", rep_id, "\n")
  }
}

# ======================
# AGGREGATE & SUMMARIZE
# ======================
cat("\n\n=== FINAL FPR RESULTS ===\n")
final_df <- bind_rows(all_results)

# FPR = proportion of null replicates
fpr_summary <- final_df %>%
  group_by(population, model) %>%
  summarise(
    n_reps = n(),
    n_reps_with_fp = sum(any_false_positive, na.rm = TRUE),
    fpr = n_reps_with_fp / n_reps,
    avg_n_fp_snps = mean(n_sig_snps, na.rm = TRUE),
    avg_n_snps_tested = mean(n_snps_tested, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(population, model)

print(fpr_summary)

# Save outputs
write.csv(final_df, file.path(RESULTS_DIR, "fpr_results_all_reps.csv"),  row.names = FALSE)
write.csv(fpr_summary, file.path(RESULTS_DIR, "fpr_summary.csv"), row.names = FALSE)

cat("\nFPR pipeline complete! Results saved to:", RESULTS_DIR, "\n")