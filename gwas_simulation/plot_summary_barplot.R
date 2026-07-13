# GWAS Power Analysis Visualization
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)  # for combining plots

# Read the data
data <- read.csv("gwas_power_results/power_summary.csv")

# Rename HBP to Chibas
data$population <- gsub("HBP", "Chibas", data$population)

# Format model names properly
data$model <- case_when(
  data$model == "glm" ~ "GLM",
  data$model == "glm_pc" ~ "GLM PC",
  data$model == "mlm_k" ~ "MLM K",
  data$model == "mlm_pk" ~ "MLM PC+K",
  TRUE ~ data$model
)

# Create a factor for window_kb for better plotting
data$window_kb_f <- factor(data$window_kb, levels = c(10, 25, 100))

# Create a factor for model with proper order
data$model <- factor(data$model, levels = c("GLM", "GLM PC", "MLM K", "MLM PC+K"))

# Define color palette
colors <- c("Chibas" = "#66C2A5", "SAP" = "#FC8D62")

# Set base theme with Helvetica font and no gridlines
theme_custom <- function() {
  theme_bw(base_size = 12, base_family = "Helvetica") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      strip.background = element_rect(fill = "gray90"),
      text = element_text(family = "Helvetica")
    )
}


# ============================================================================
# Updated Plot: Compare populations side-by-side, faceted by model
# ============================================================================
p_final <- ggplot(data, aes(x = window_kb_f, y = power, fill = population)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", linewidth = 0.3) +
  facet_wrap(~model, ncol = 2) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(title = "Detection Power by Model and Window Size",
       x = "Window (kb)",
       y = "Power",
       fill = "Population") +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90"),
    text = element_text(family = "Helvetica")
  )

# Save as PDF
cairo_pdf("gwas_power_results/power_plot_WITHOUT_maf_constraints_h05.pdf", width = 10, height = 6, family = "Helvetica")
print(p_final)
dev.off()

cat("\nFinal plot saved!\n")



sample_gene <- "DX9009|TCON100011"
strsplit(sample_gene, "\\|")[[1]][1]


# ============================================================================
# Plot 1: Power by Window Size (faceted by model)
# ============================================================================
p1 <- ggplot(data, aes(x = window_kb_f, y = power, color = population, group = population)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  facet_wrap(~model, ncol = 3) +
  scale_color_manual(values = colors) +
  labs(title = "Detection Power by Window Size and Model",
       x = "Window Size (kb)",
       y = "Power",
       color = "Population") +
  theme_custom()

# ============================================================================
# Plot 2: Power comparison across models (faceted by population)
# ============================================================================
p2 <- ggplot(data, aes(x = model, y = power, fill = window_kb_f)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  facet_wrap(~population) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Detection Power by Model and Window Size",
       x = "Model",
       y = "Power",
       fill = "Window (kb)") +
  theme_custom() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================================================
# Plot 3: Average Distance comparison
# ============================================================================
p3 <- ggplot(data, aes(x = window_kb_f, y = avg_distance, color = population, group = population)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  facet_wrap(~model, ncol = 3, scales = "free_y") +
  scale_color_manual(values = colors) +
  labs(title = "Average Distance by Window Size and Model",
       x = "Window Size (kb)",
       y = "Average Distance",
       color = "Population") +
  theme_custom()

# ============================================================================
# Plot 4: Power vs Average Distance scatter plot
# ============================================================================
p4 <- ggplot(data, aes(x = avg_distance, y = power, color = population, shape = model)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(title = "Power vs Average Distance",
       x = "Average Distance",
       y = "Power",
       color = "Population",
       shape = "Model") +
  theme_custom()

# ============================================================================
# Plot 5: Heatmap of power values
# ============================================================================
p5 <- ggplot(data, aes(x = window_kb_f, y = model, fill = power)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", power)), color = "black", size = 4, family = "Helvetica") +
  facet_wrap(~population) +
  scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850", 
                       midpoint = 0.85, limits = c(0.6, 1)) +
  labs(title = "Power Heatmap",
       x = "Window Size (kb)",
       y = "Model",
       fill = "Power") +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    text = element_text(family = "Helvetica")
  )

# ============================================================================
# Save individual plots as PDF using Cairo
# ============================================================================
cairo_pdf("gwas_power_results/plot1_power_by_window.pdf", width = 10, height = 6, family = "Helvetica")
print(p1)
dev.off()

cairo_pdf("gwas_power_results/plot2_power_by_model.pdf", width = 10, height = 6, family = "Helvetica")
print(p2)
dev.off()

cairo_pdf("gwas_power_results/plot3_avg_distance.pdf", width = 10, height = 6, family = "Helvetica")
print(p3)
dev.off()

cairo_pdf("gwas_power_results/plot4_power_vs_distance.pdf", width = 8, height = 6, family = "Helvetica")
print(p4)
dev.off()

cairo_pdf("gwas_power_results/plot5_heatmap.pdf", width = 10, height = 5, family = "Helvetica")
print(p5)
dev.off()

# ============================================================================
# Create combined summary plot
# ============================================================================
combined_plot <- (p1 / p2) | p5
cairo_pdf("gwas_power_results/combined_summary.pdf", width = 16, height = 10, family = "Helvetica")
print(combined_plot)
dev.off()

# ============================================================================
# Print summary statistics
# ============================================================================
cat("\n=== Summary Statistics ===\n")
summary_stats <- data %>%
  group_by(population, model) %>%
  summarise(
    mean_power = mean(power),
    max_power = max(power),
    min_power = min(power),
    mean_distance = mean(avg_distance),
    .groups = "drop"
  )
print(summary_stats)

cat("\n=== Best performing combinations (highest power) ===\n")
best <- data %>%
  arrange(desc(power)) %>%
  head(5)
print(best)

cat("\nPlots saved as PDF files to gwas_power_results/ directory\n")


# ============================================================================
# Updated Plot: Power by Window Size, grouped by Model, faceted by Population
# ============================================================================
p_final <- ggplot(data, aes(x = window_kb_f, y = power, fill = model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", linewidth = 0.3) +
  facet_wrap(~population, ncol = 2) +
  scale_fill_manual(values = c("GLM" = "#66C2A5", "MLM PC" = "#FC8D62", "MLM PC+K" = "#8DA0CB")) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(title = "Detection Power by Model and Window Size",
       x = "Window (kb)",
       y = "Power",
       fill = "Model") +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90"),
    text = element_text(family = "Helvetica")
  )

# Save as PDF
cairo_pdf("gwas_power_results/final_power_plot.pdf", width = 10, height = 6, family = "Helvetica")
print(p_final)
dev.off()

cat("\nFinal plot saved as: gwas_power_results/final_power_plot.pdf\n")
