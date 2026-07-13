plot_all_schemes_fair_with_lm <- function(csv_file) {
  data <- read.csv(csv_file)
  
  # Convert founder_size to factor for proper legend
  data$founder_size <- factor(data$founder_size)

  # Scheme aesthetics — only apply levels/colors present in the data
  scheme_levels  <- c("Genetic Male Sterility", "Chemical Sterilization", "Conventional Inbred")
  scheme_colors  <- c("Genetic Male Sterility" = "#E41A1C",
                      "Chemical Sterilization"  = "#377EB8",
                      "Conventional Inbred"     = "#984EA3")
  scheme_alphas  <- c("Genetic Male Sterility" = 1.0,
                      "Chemical Sterilization"  = 1.0,
                      "Conventional Inbred"     = 0.5)
  active_schemes <- intersect(scheme_levels, unique(data$scheme))
  data$scheme    <- factor(data$scheme, levels = active_schemes)

  # Rapid cycling = 1 generation/cycle; conventional = n.gens (5) generations/cycle
  data <- data %>%
    mutate(
      total_generations = ifelse(scheme %in% c("Genetic Male Sterility", "Chemical Sterilization"),
                                 recurrent_cycle,
                                 recurrent_cycle * n.gens)
    )
  
  # Modified function to create individual plots for each metric
  make_individual_plot <- function(df, yvar, ylab, title, filename_suffix) {
    # Summarize across replicates so we get the true mean at each generation
    summary_data <- df %>%
      dplyr::group_by(scheme, total_generations) %>%
      dplyr::summarise(
        mean_y = mean(.data[[yvar]], na.rm = TRUE),
        se = sd(.data[[yvar]], na.rm = TRUE) / sqrt(dplyr::n()),
        .groups = "drop"
      )
    
    p <- ggplot(summary_data, aes(x = total_generations, y = mean_y, color = scheme, group = scheme)) +
      # true means with lines + points
       geom_line(aes(alpha = scheme), linewidth = 1.2) +
       geom_point(aes(alpha = scheme), size = 3) +
      #geom_line(linewidth = 1.2) +
      #geom_point(aes(shape = scheme), size = 3, fill = "white") +
      #scale_shape_manual(values = c(16, 17, 15, 18)) +  # circle, triangle, square, diamond
      # shaded SE ribbon
      geom_ribbon(aes(ymin = mean_y - se, ymax = mean_y + se, fill = scheme), alpha = 0.15, color = NA) +
      # optional smooth overlay for trend
      #geom_smooth(se = FALSE, method = "loess", formula = y ~ x, span = 0.75, linewidth = 1, linetype = "dashed") +
      
      theme_minimal() +
      labs(x = "Generation", y = ylab, title = title) +
      scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +
      #scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
      theme(
        text = element_text(family = "Helvetica"),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(15, 70, 15, 30) # Add margins: top, right, bottom, left
      ) +
      scale_color_manual(values = scheme_colors[active_schemes]) +
      scale_fill_manual(values  = scheme_colors[active_schemes]) +
      scale_alpha_manual(values = scheme_alphas[active_schemes]) +
      
      guides(alpha = "none")  # hides alpha from legend
    
    # Save individual plot
    # plot_timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    # ggsave(paste0(filename_suffix, "_", plot_timestamp, ".jpg"), 
    #        p, width = 10, height = 8, dpi = 600)
    
    plot_timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    ggsave(
      paste0(filename_suffix, "_", plot_timestamp, ".pdf"),
      width = 10, height = 8, dpi = 600,
      device = cairo_pdf,
      family = "Arial"
    )
    
    return(p)
  }
  
  # Create and save individual plots
  gain_plot <- make_individual_plot(data, "geno_mean", "Phenotypic Mean", "", "phenotypic_mean")
  var_plot <- make_individual_plot(data, "geno_var", "Genetic Variance", "", "genetic_variance")
  maf_plot <- make_individual_plot(data, "avg_maf", "Average MAF", "", "average_maf")
  fixed_qtl_plot <- make_individual_plot(data, "n_fixed_qtl", "Fixed QTL", "", "fixed_qtl")
  ne_plot <- make_individual_plot(data, "Ne_like_glLDNe", "Effective Population Size (Ne)", "", "effective_pop_size")
  pi_plot <- make_individual_plot(data, "nucleotide_diversity", "Nucleotide Diversity", "", "nucleo_diversity")
  bgld_plot <- make_individual_plot(data, "background_ld", "Background LD", "", "background_ld")
  
  
  # Create combined plot for comparison (optional)
  combined_plot <- (gain_plot + var_plot) / (maf_plot + fixed_qtl_plot) / (ne_plot + pi_plot) / (bgld_plot + plot_spacer()) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  # Add overall title to combined plot
  combined_plot <- combined_plot +
    plot_annotation(title = "Comparison of Breeding Schemes (Adjusted for Generations)",
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 26)))
  
  # Save combined plot
  # plot_timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  # ggsave(paste0("fair_comparison_generations_", plot_timestamp, ".jpg"), 
  #        combined_plot, width = 18, height = 18, dpi = 600)
  
  plot_timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  ggsave(
    paste0("fair_comparison_generations_", plot_timestamp, ".pdf"),
    width = 18, height = 18, dpi = 600,
    device = cairo_pdf,
    family = "Arial"
  )
  
  # Return the combined plot for display
  return(combined_plot)
}