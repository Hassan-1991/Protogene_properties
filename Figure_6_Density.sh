library(ggplot2)
library(dplyr)
library(gridExtra)

###DENSITY PLOT###

# Step 1: Load data
data <- read.csv("/stor/work/Ochman/hassan/protogene_extension/expression_location_properties/sequence_properties/sequenceproperties.marked.conservation.final.plotting.csv")

# Step 2: Define feature-column mapping with appropriate ylim values
features <- list(
  list(name = "GC", ylim = c(0, 0.15)),
  list(name = "polarAA", ylim = c(0, 0.125)),
  list(name = "cai", ylim = c(0, 10))
)

# Step 3: Define the annotation pairs (these define the rows now)
pairs <- list(
  c("annotated_conserved", "control_control"),
  c("annotated_ORFan", "control_control"),
  c("novel_conserved", "control_control"),
  c("novel_ORFan", "control_control"),
  c("novel_total", "control_control")
)

# Step 4: Function to generate plots for a given feature
generate_plots <- function(feature_name, ylim_vals) {
  # Filter data for species "Ecoli" and the given feature
  ecoli_data <- data %>%
    filter(species == "Mycobacterium", feature == feature_name)
  
  # Ensure control_control is always present
  control_data <- ecoli_data %>%
    filter(annot_status_conservation == "control_control")
  
  # Generate density plots
  plots <- lapply(pairs, function(groups) {
    # Filter data for the current pair of groups
    filtered_data <- ecoli_data %>%
      filter(annot_status_conservation %in% groups)
    
    # Make sure control group is included if missing
    if (!any(filtered_data$annot_status_conservation == "control_control")) {
      filtered_data <- bind_rows(filtered_data, control_data)
    }
    
    # Explicitly set the color mapping to avoid ggplot's automatic reassignment
    filtered_data$annot_status_conservation <- factor(
      filtered_data$annot_status_conservation, 
      levels = c("control_control", "annotated_conserved", "annotated_ORFan", "novel_conserved", "novel_ORFan", "novel_total")
    )
    
    # Create the density plot
    ggplot(filtered_data, aes(x = value, color = annot_status_conservation, fill = annot_status_conservation)) +
      geom_density(alpha = 0.4, adjust = 1) +
      labs(x = feature_name, y = "Density") +
      scale_color_manual(values = c(
        "control_control" = "green", 
        "annotated_conserved" = "darkmagenta", 
        "annotated_ORFan" = "pink", 
        "novel_conserved" = "blue", 
        "novel_ORFan" = "skyblue", 
        "novel_total" = "blueviolet"
      )) +
      scale_fill_manual(values = c(
        "control_control" = "green", 
        "annotated_conserved" = "darkmagenta", 
        "annotated_ORFan" = "pink", 
        "novel_conserved" = "blue", 
        "novel_ORFan" = "skyblue", 
        "novel_total" = "blueviolet"
      )) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.y = element_line(size = 0.25),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1)
      ) +
      ylim(ylim_vals)  # Set y-axis limit
  })
  
  return(plots)
}

# Step 5: Generate plots for each feature
all_plots <- lapply(features, function(f) generate_plots(f$name, f$ylim))

# Step 6: Transpose the layout - extract each row from different features
transposed_plots <- list()
for (i in 1:length(pairs)) {
  transposed_plots <- append(transposed_plots, list(all_plots[[1]][[i]], all_plots[[2]][[i]], all_plots[[3]][[i]]))
}

# Step 7: Print the final 5x3 collage (5 rows, 3 columns)
p1 <- grid.arrange(grobs = transposed_plots, ncol = 3)

#ggsave("/stor/scratch/Ochman/hassan/112724_protogene_extension/consolidating_datasets/Ecoli_density.pdf", p1, width = 24, height = 30, units = "in")
ggsave("/stor/scratch/Ochman/hassan/112724_protogene_extension/consolidating_datasets/Mycobacterium_density.pdf", p1, width = 24, height = 30, units = "in")

###VIOLIN PLOT###
# Function to generate plot for a dataset
generate_plot <- function(data, variable) {
  data %>%
    filter(feature == variable) %>%
    ggplot(aes(x = annot_status_conservation, y = value)) +
    geom_violin(aes(fill = annot_status_conservation), color = "black", alpha = 1, width = 0.9) +
    stat_summary(fun.data = mean_cl_boot, geom = "point", shape = 18, size = 3, color = "yellow") +
    scale_fill_manual(values = c("maroon", "pink", "mediumblue", "lightskyblue2", "blueviolet","grey")) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),         # Remove X-axis labels
      axis.text.x = element_blank(),          # Remove X-axis text
      axis.text.y = element_text(size = 20),  # Increase Y-axis text size
      axis.ticks.x = element_blank(),         # Remove X-axis ticks
      axis.ticks.y = element_line(size = 1.5),# Make Y-axis ticks longer
      axis.ticks.length = unit(0.3, "cm"),    # Increase tick length
      panel.border = element_rect(color = "black", fill = NA, size = 1.5), # Box around the plot
      panel.grid.major.x = element_blank(),   # Remove vertical gridlines
      panel.grid.minor.x = element_blank(),
      legend.position = "none"
    )
}

# Read the datasets
df_ecoli <- read.csv("/stor/work/Ochman/hassan/protogene_extension/expression_location_properties/sequence_properties/sequenceproperties.marked.conservation.final.plotting.csv") %>% filter(species=="Ecoli") %>% filter(annot_status_conservation!="novel_others") %>% filter(annot_status_conservation!="annotated_others")
df_mt <- read.csv("/stor/work/Ochman/hassan/protogene_extension/expression_location_properties/sequence_properties/sequenceproperties.marked.conservation.final.plotting.csv") %>% filter(species=="Mycobacterium") %>% filter(annot_status_conservation!="novel_others") %>% filter(annot_status_conservation!="annotated_others")

# Define levels for factor variable
df_mt$annot_status_conservation <- factor(df_mt$annot_status_conservation, levels = c("annotated_conserved", "annotated_ORFan", "novel_conserved", "novel_ORFan", "novel_total", "control_control"))
df_ecoli$annot_status_conservation <- factor(df_ecoli$annot_status_conservation, levels = c("annotated_conserved", "annotated_ORFan", "novel_conserved", "novel_ORFan", "novel_total", "control_control"))

# Generate plots for MTb dataset
plot_mt_gc3rd <- generate_plot(df_mt, "GC")
plot_mt_compbias <- generate_plot(df_mt, "polarAA")
plot_mt_cai <- generate_plot(df_mt, "cai")

# Generate plots for Ecoli dataset
plot_ecoli_gc3rd <- generate_plot(df_ecoli, "GC")
plot_ecoli_compbias <- generate_plot(df_ecoli, "polarAA")
plot_ecoli_cai <- generate_plot(df_ecoli, "cai")

# Arrange plots into a 2x2 grid
combined_plots <- grid.arrange(
  plot_ecoli_gc3rd, plot_ecoli_compbias, plot_ecoli_cai,
  plot_mt_gc3rd, plot_mt_compbias, plot_mt_cai,
  ncol = 3, top = "Comparison of GC_3rd and CAI across E. coli and MTb"
)

# Print the collage
print(combined_plots)

ggsave("/stor/scratch/Ochman/hassan/112724_protogene_extension/consolidating_datasets/combined_violin_plots.pdf", plot = combined_plots, width = 18, height = 10, units = "in", device = "pdf")

###SIGNIFICANCE VALUES###
compute_wilcoxon <- function(data, feature) {
  filtered_data <- data %>% filter(feature == !!feature)
  
  # Perform pairwise Wilcoxon test
  test_results <- pairwise.wilcox.test(filtered_data$value, filtered_data$annot_status_conservation, 
                                       p.adjust.method = "bonferroni")  # Adjust for multiple comparisons
  
  # Convert to a tidy dataframe
  result_df <- as.data.frame(as.table(test_results$p.value)) %>%
    rename(Category1 = Var1, Category2 = Var2, p_value = Freq) %>%
    mutate(Feature = feature)
  
  return(result_df)
}

# Compute pairwise Wilcoxon test results for each feature in E. coli
wilcoxon_gc <- compute_wilcoxon(df_ecoli, "GC")
wilcoxon_polarAA <- compute_wilcoxon(df_ecoli, "polarAA")
wilcoxon_cai <- compute_wilcoxon(df_ecoli, "cai")

# Combine results into one dataframe
wilcoxon_results <- bind_rows(wilcoxon_gc, wilcoxon_polarAA, wilcoxon_cai)

# Print or save the results
print(wilcoxon_results)
write.csv(wilcoxon_results, "/stor/scratch/Ochman/hassan/112724_protogene_extension/consolidating_datasets/wilcoxon_p_values.csv", row.names = FALSE)
