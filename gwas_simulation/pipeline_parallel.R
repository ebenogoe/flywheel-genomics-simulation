# GWAS Simulation Pipeline — Parallelized
# Strategy: future + furrr with per-rep isolated directories
# Each worker gets its own subdirectory to avoid file collisions in simplePHENOTYPES

library(GAPIT)
library(simplePHENOTYPES)
library(dplyr)
library(readr)
library(stringr)
library(future)
library(furrr)       # install.packages("furrr") if needed
library(progressr)   # install.packages("progressr") if needed

# ======================
# CONFIGURATION PARAMS
# ======================
N_REPS        <- 100
WINDOW_SIZES  <- c(10, 25, 100)  # kb
H2            <- 0.5
ADD_EFFECT    <- 1
MAF_THRESHOLD <- 0
N_PC          <- 3

# How many parallel workers to use.
# 6 workers is a safe default for M2 Pro (6 perf cores); drop to 4 if you see memory pressure.
N_WORKERS <- 8

RESULTS_DIR <- "gwas_power_results"
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ======================
# HELPER FUNCTIONS
# ======================

# Create rep-level directory only (let simplePHENOTYPES create population subdirs)
make_rep_dir <- function(base_dir, rep_id) {
  path <- file.path(base_dir, paste0("rep_", rep_id))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path)
}

# Remove rep scratch dir when done (keeps disk clean)
clean_rep_dir <- function(base_dir, rep_id) {
  path <- file.path(base_dir, paste0("rep_", rep_id))
  if (dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
}

# Read QTNs from simplePHENOTYPES output in a per-rep directory
read_qtns <- function(pheno_dir) {
  qtn_file <- file.path(pheno_dir, "Additive_Selected_QTNs.txt")
  if (!file.exists(qtn_file)) stop("QTN file not found: ", qtn_file)
  
  qtns <- read_delim(qtn_file, delim = "\t", col_types = cols())
  
  if ("snp" %in% colnames(qtns)) {
    qtns$snp_id <- qtns$snp
  } else if ("SNP" %in% colnames(qtns)) {
    qtns$snp_id <- qtns$SNP
  }
  
  qtns$chr <- as.character(str_extract(qtns$snp_id, "(?<=S?)[0-9]+(?=_)"))
  qtns$pos <- as.numeric(str_extract(qtns$snp_id, "(?<=_)[0-9]+$"))
  
  qtns %>% select(snp_id, chr, pos)
}

# Remove QTN SNP from genotype matrix
remove_qtn_from_geno <- function(geno_df, qtn_snp, pop_type = c("SAP", "HBP")) {
  pop_type  <- match.arg(pop_type)
  qtn_clean <- gsub("^S", "", qtn_snp)
  
  qtn_to_remove <- if (pop_type == "SAP") paste0("S", qtn_clean) else qtn_clean
  
  if (!(qtn_to_remove %in% rownames(geno_df))) {
    warning("QTN SNP not found in genotype data: ", qtn_to_remove)
    return(geno_df)
  }
  
  geno_filtered <- geno_df[rownames(geno_df) != qtn_to_remove, , drop = FALSE]
  message("  Removed QTN: ", qtn_to_remove, " | Remaining SNPs: ", nrow(geno_filtered))
  geno_filtered
}

# Align kinship matrix to sample IDs
align_kinship <- function(kinship_mat, sample_ids) {
  if (is.null(rownames(kinship_mat)) ||
      all(rownames(kinship_mat) %in% as.character(seq_len(ncol(kinship_mat))))) {
    if (length(sample_ids) != nrow(kinship_mat))
      stop("Sample count mismatch: kinship (", nrow(kinship_mat),
           ") vs samples (", length(sample_ids), ")")
    rownames(kinship_mat) <- sample_ids
    colnames(kinship_mat) <- sample_ids
  } else {
    common_ids   <- intersect(rownames(kinship_mat), sample_ids)
    kinship_mat  <- kinship_mat[common_ids, common_ids, drop = FALSE]
  }
  
  data.frame(
    ID = rownames(kinship_mat),
    as.data.frame(kinship_mat),
    check.names    = FALSE,
    stringsAsFactors = FALSE
  )
}

# Format genotype data for GAPIT
format_geno_for_gapit <- function(geno_df, pop_type = c("SAP", "HBP")) {
  pop_type  <- match.arg(pop_type)
  snp_names <- rownames(geno_df)
  
  if (pop_type == "HBP") {
    snp_info <- strsplit(snp_names, "_")
  } else {
    snp_info <- strsplit(gsub("^S", "", snp_names), "_")
  }
  
  chrom <- sapply(snp_info, `[`, 1)
  pos   <- sapply(snp_info, `[`, 2)
  
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
  
  list(GM = GM, GD = GD)
}

# Run all 3 GWAS models
run_gwas_models <- function(Y, GD, GM, pcs_df, kinship_mat, pop_name, rep_id, n_pc) {
  results <- list()
  
  message("  [Rep ", rep_id, "] GLM — ", pop_name)
  results$glm <- GAPIT(
    Y = Y, GD = GD, GM = GM,
    PCA.total          = 0,
    kinship.algorithm  = NULL,
    SNP.MAF            = MAF_THRESHOLD,
    model              = "GLM",
    file.output        = FALSE,
    Geno.View.output   = FALSE,
    PCA.View.output    = FALSE
  )
  
  message("  [Rep ", rep_id, "] MLM+PC — ", pop_name)
  pc_match <- pcs_df[match(Y$Taxa, pcs_df$sample),
                     c("sample", paste0("PC", seq_len(n_pc)))]
  colnames(pc_match)[1] <- "Taxa"
  
  results$mlm_pc <- GAPIT(
    Y = Y, GD = GD, GM = GM, CV = pc_match,
    PCA.total          = 0,
    kinship.algorithm  = NULL,
    SNP.MAF            = MAF_THRESHOLD,
    model              = "MLM",
    file.output        = FALSE,
    Geno.View.output   = FALSE,
    PCA.View.output    = FALSE
  )
  
  message("  [Rep ", rep_id, "] MLM+PC+K — ", pop_name)
  kinship_aligned <- align_kinship(kinship_mat, Y$Taxa)
  
  results$mlm_pk <- GAPIT(
    Y = Y, GD = GD, GM = GM, CV = pc_match, KI = kinship_aligned,
    SNP.MAF            = MAF_THRESHOLD,
    model              = "MLM",
    file.output        = FALSE,
    Geno.View.output   = FALSE,
    PCA.View.output    = FALSE
  )
  
  results
}

# Check QTN detection across window sizes
check_detection_multiwindow <- function(gwas_result, true_qtn, window_kb_vec) {
  gwas_df          <- gwas_result$GWAS
  gwas_df$Position <- as.numeric(gsub(",", "", as.character(gwas_df$Position)))
  bonf_thresh      <- 0.05 / nrow(gwas_df)
  sig_snps         <- gwas_df[gwas_df$P.value < bonf_thresh, ]
  
  results <- lapply(window_kb_vec, function(window_kb) {
    window_bp <- window_kb * 1000
    
    if (nrow(sig_snps) > 0) {
      chr_sig <- sig_snps[sig_snps$Chromosome == true_qtn$chr, ]
      if (nrow(chr_sig) > 0) {
        distances   <- abs(chr_sig$Position - true_qtn$pos)
        min_dist    <- min(distances)
        detected    <- min_dist <= window_bp
        closest_snp <- chr_sig$SNP[which.min(distances)]
      } else {
        min_dist <- NA; detected <- FALSE; closest_snp <- NA
      }
    } else {
      min_dist <- NA; detected <- FALSE; closest_snp <- NA
    }
    
    data.frame(
      window_kb   = window_kb,
      detected    = detected,
      distance_bp = min_dist,
      closest_snp = closest_snp,
      n_sig_snps  = nrow(sig_snps),
      bonf_thresh = bonf_thresh
    )
  })
  
  bind_rows(results)
}

# ======================
# PER-REP WORKER FUNCTION
# All data objects passed explicitly — workers are isolated processes
# ======================
run_one_rep <- function(
    rep_id,
    scratch_dir,
    # SAP
    sap_geno_df, sap_numeric_format_filtered, sap_pcs, sap_kinship,
    # HBP
    hbp_geno_df, hbp_numeric_format_filtered, hbp_pcs, hbp_kinship,
    # Config
    h2, add_effect, n_pc, window_sizes, maf_threshold
) {
  message("\n=== Rep ", rep_id, " started (PID ", Sys.getpid(), ") ===")
  
  # ---- SAP ----
  message("[Rep ", rep_id, "] SAP: simulating phenotype...")
  rep_dir <- make_rep_dir(scratch_dir, rep_id)
  
  # Pre-create SAP parent directory (simplePHENOTYPES will create SAP_phenotypes)
  dir.create(file.path(rep_dir, "SAP"), showWarnings = FALSE, recursive = TRUE)
  
  sap_pheno_obj <- create_phenotypes(
    geno_obj    = sap_numeric_format_filtered,
    add_QTN_num = 1,
    h2          = h2,
    model       = "A",
    add_effect  = add_effect,
    rep         = 1,
    home_dir    = rep_dir,
    constraints = list(maf_above = 0.05, maf_below = 0.06),
    output_dir  = file.path("SAP", "SAP_phenotypes")
  )
  
  sap_pheno_dir <- file.path(rep_dir, "SAP", "SAP_phenotypes")
  sap_qtns <- read_qtns(sap_pheno_dir)
  message("[Rep ", rep_id, "] SAP QTN: ", sap_qtns$snp_id,
          " chr", sap_qtns$chr, " pos", sap_qtns$pos)
  
  sap_geno_filtered <- remove_qtn_from_geno(sap_geno_df, sap_qtns$snp_id, "SAP")
  sap_gapit         <- format_geno_for_gapit(sap_geno_filtered, "SAP")
  
  sap_pheno_file <- file.path(sap_pheno_dir,
                              "Simulated_Data_1_Reps_Herit_0.5.txt")
  sap_pheno <- read.table(sap_pheno_file, header = TRUE)
  sap_pheno_clean <- data.frame(
    Taxa  = sap_pheno$X.Trait.,
    Pheno = sap_pheno$Pheno,
    stringsAsFactors = FALSE
  )
  
  sap_gwas <- run_gwas_models(
    Y = sap_pheno_clean, GD = sap_gapit$GD, GM = sap_gapit$GM,
    pcs_df = sap_pcs, kinship_mat = sap_kinship,
    pop_name = "SAP", rep_id = rep_id, n_pc = n_pc
  )
  
  sap_detection_df <- bind_rows(lapply(c("glm", "mlm_pc", "mlm_pk"), function(model) {
    check_detection_multiwindow(sap_gwas[[model]], sap_qtns[1, ], window_sizes) %>%
      mutate(model = model)
  })) %>% mutate(population = "SAP")
  
  # ---- HBP ----
  message("[Rep ", rep_id, "] HBP: simulating phenotype...")
  
  # Pre-create HBP parent directory (simplePHENOTYPES will create HBP_phenotypes)
  dir.create(file.path(rep_dir, "HBP"), showWarnings = FALSE, recursive = TRUE)
  
  hbp_pheno_obj <- create_phenotypes(
    geno_obj    = hbp_numeric_format_filtered,
    add_QTN_num = 1,
    h2          = h2,
    model       = "A",
    add_effect  = add_effect,
    rep         = 1,
    home_dir    = rep_dir,
    constraints = list(maf_above = 0.05, maf_below = 0.06),
    output_dir  = file.path("HBP", "HBP_phenotypes")
  )
  
  hbp_pheno_dir <- file.path(rep_dir, "HBP", "HBP_phenotypes")
  hbp_qtns <- read_qtns(hbp_pheno_dir)
  message("[Rep ", rep_id, "] HBP QTN: ", hbp_qtns$snp_id,
          " chr", hbp_qtns$chr, " pos", hbp_qtns$pos)
  
  hbp_geno_filtered <- remove_qtn_from_geno(hbp_geno_df, hbp_qtns$snp_id, "HBP")
  hbp_gapit         <- format_geno_for_gapit(hbp_geno_filtered, "HBP")
  
  hbp_pheno_file <- file.path(hbp_pheno_dir,
                              "Simulated_Data_1_Reps_Herit_0.5.txt")
  hbp_pheno <- read.table(hbp_pheno_file, header = TRUE)
  hbp_pheno_clean <- data.frame(
    Taxa  = hbp_pheno$X.Trait.,
    Pheno = hbp_pheno$Pheno,
    stringsAsFactors = FALSE
  )
  
  hbp_gwas <- run_gwas_models(
    Y = hbp_pheno_clean, GD = hbp_gapit$GD, GM = hbp_gapit$GM,
    pcs_df = hbp_pcs, kinship_mat = hbp_kinship,
    pop_name = "HBP", rep_id = rep_id, n_pc = n_pc
  )
  
  hbp_detection_df <- bind_rows(lapply(c("glm", "mlm_pc", "mlm_pk"), function(model) {
    check_detection_multiwindow(hbp_gwas[[model]], hbp_qtns[1, ], window_sizes) %>%
      mutate(model = model)
  })) %>% mutate(population = "HBP")
  
  # ---- Combine & clean up ----
  clean_rep_dir(scratch_dir, rep_id)
  
  bind_rows(sap_detection_df, hbp_detection_df) %>%
    mutate(replicate = rep_id)
}


# ======================
# MAIN PIPELINE
# ======================

cat("Loading base genotype data...\n")
load("sap_geno_df_MAF0.05.RData")
load("hbp_geno_df_chibas_MAF0.05.RData")

cat("Loading PCs and kinship...\n")
sap_pcs     <- readRDS("SAP/SAP_PCs.rds")
hbp_pcs     <- readRDS("HBP/HBP_PCs.rds")
sap_kinship <- readRDS("SAP/SAP_kinship.rds")
hbp_kinship <- readRDS("HBP/HBP_kinship.rds")

# Scratch dir for per-rep isolated phenotype files (use absolute path)
SCRATCH_DIR <- normalizePath(file.path(RESULTS_DIR, "scratch"), mustWork = FALSE)
dir.create(SCRATCH_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Set up parallel backend ----
# multisession = separate R processes (fork-safe, works on macOS/M2)
# Each worker gets a full copy of the exported data via future's globals mechanism
plan(multisession, workers = N_WORKERS)

# ---- Progress reporting ----
# progressr sends progress updates from workers back to the main process
handlers(global = TRUE)
handlers("cli")   # Uses cli package for a clean progress bar; falls back to "txtprogressbar"

cat(sprintf(
  "\nStarting %d replicates across %d workers...\n\n",
  N_REPS, N_WORKERS
))

# ---- Run in parallel with progress ----
all_results <- with_progress({
  p <- progressor(steps = N_REPS)
  
  future_map(
    seq_len(N_REPS),
    function(rep_id) {
      # Signal progress BEFORE the heavy work so the bar updates when a worker picks up a job
      p(message = sprintf("Rep %d/%d", rep_id, N_REPS))
      
      run_one_rep(
        rep_id = rep_id,
        scratch_dir = SCRATCH_DIR,
        # SAP
        sap_geno_df                  = sap_geno_df,
        sap_numeric_format_filtered  = sap_numeric_format_filtered,
        sap_pcs                      = sap_pcs,
        sap_kinship                  = sap_kinship,
        # HBP
        hbp_geno_df                  = hbp_geno_df,
        hbp_numeric_format_filtered  = hbp_numeric_format_filtered,
        hbp_pcs                      = hbp_pcs,
        hbp_kinship                  = hbp_kinship,
        # Config
        h2            = H2,
        add_effect    = ADD_EFFECT,
        n_pc          = N_PC,
        window_sizes  = WINDOW_SIZES,
        maf_threshold = MAF_THRESHOLD
      )
    },
    .options = furrr_options(
      seed     = TRUE,           # Reproducible RNG across workers
      chunk_size = 1,            # Send reps one at a time for even load balancing
      packages = c(              # Explicitly declare packages each worker needs
        "GAPIT", "simplePHENOTYPES",
        "dplyr", "readr", "stringr"
      )
    )
  )
})

# Reset to sequential
plan(sequential)

# ======================
# CHECKPOINTING
# Incremental saves happen after every 10 completed reps (regardless of order)
# ======================
final_df <- bind_rows(all_results)

# Save checkpoint every 10 reps based on completed replicate numbers
completed_reps <- sort(unique(final_df$replicate))
for (i in seq_along(completed_reps)) {
  if (i %% 10 == 0) {
    saveRDS(
      final_df %>% filter(replicate %in% completed_reps[seq_len(i)]),
      file.path(RESULTS_DIR, sprintf("power_results_checkpoint_%03d.RDS", i))
    )
    cat("[Checkpoint] Saved results for", i, "completed replicates\n")
  }
}

# ======================
# AGGREGATE & SUMMARIZE
# ======================
cat("\n\n=== FINAL RESULTS ===\n")

power_summary <- final_df %>%
  group_by(population, model, window_kb) %>%
  summarise(
    n_detected   = sum(detected, na.rm = TRUE),
    n_total      = n(),
    power        = n_detected / n_total,
    avg_distance = mean(distance_bp, na.rm = TRUE),
    .groups      = "drop"
  ) %>%
  arrange(population, model, window_kb)

print(power_summary)

write.csv(final_df,      file.path(RESULTS_DIR, "detection_results_all_reps.csv"), row.names = FALSE)
write.csv(power_summary, file.path(RESULTS_DIR, "power_summary.csv"),              row.names = FALSE)

cat("\nPipeline complete! Results saved to:", RESULTS_DIR, "\n")