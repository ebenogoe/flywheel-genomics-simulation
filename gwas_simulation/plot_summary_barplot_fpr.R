# ============================================================================
# FPR Plot: Compare models side-by-side, colored by population
# ============================================================================
library(ggplot2)
library(dplyr)

# Load summary results
fpr_summary <- read.csv("gwas_fpr_results/fpr_summary.csv")

# Clean model labels for display
fpr_summary <- fpr_summary %>%
  mutate(
    model_f = factor(model,
                     levels = c("glm", "glm_pc", "mlm_k", "mlm_pk"),
                     labels = c("GLM", "GLM + PC", "MLM + K", "MLM + PC + K")),
    population = factor(population, levels = c("SAP", "HBP"))
  )

colors <- c("SAP" = "#2166AC", "HBP" = "#D6604D")

p_fpr <- ggplot(fpr_summary, aes(x = model_f, y = fpr, fill = population)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),
           color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", linewidth = 0.5) +
  annotate("text", x = 0.5, y = 0.055, label = "Expected FPR (0.05)",
           hjust = 0, vjust = 0, size = 3.2, color = "black") +
  scale_fill_manual(values = colors) +
  scale_y_continuous(limits = c(0, .5), breaks = seq(0, 0.5, 0.1)) +
  labs(
    title = "Genome-Wide False Positive Rate by Model",
    subtitle = "Bonferroni threshold (α = 0.05); dashed line = nominal FPR",
    x     = "Model",
    y     = "False Positive Rate",
    fill  = "Population"
  ) +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "bottom",
    strip.background  = element_rect(fill = "gray90"),
    text              = element_text(family = "Helvetica")
  )

print(p_fpr)

# Save as PDF
cairo_pdf("gwas_fpr_results/fpr_plot.pdf", width = 7, height = 6, family = "Helvetica")
print(p_fpr)
dev.off()

cat("\nFPR plot saved as: gwas_fpr_results/fpr_plot.pdf\n")