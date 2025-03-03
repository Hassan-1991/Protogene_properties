#Dumping all code here, will re-arrange/format later:

#Expression, cumdist

p1 <- read.csv("Ecoli.REL606.1tpm.conditions.conservation.tsv",sep='\t')
p2 <- read.csv("Ecoli.REL606.0.3tpm.conditions.conservation.tsv",sep='\t')

generate_plot_data <- function(df) {
  cutoffs <- seq(-0.001, 41, by = 1)
  plot_data_list <- list()
  
  for (line_value in unique(df$category)) {
    line_data <- df[df$category == line_value, ]
    fraction_passed <- numeric(length(cutoffs))
    
    for (i in seq_along(cutoffs)) {
      fraction_passed[i] <- sum(line_data$conditions > cutoffs[i]) / nrow(line_data)
    }
    
    plot_data <- data.frame(Cutoff = cutoffs, Fraction_Passed = fraction_passed, Line = line_value)
    plot_data_list[[line_value]] <- plot_data
  }
  
  combined_plot_data <- do.call(rbind, plot_data_list)
  
  return(combined_plot_data)
}

# Generate plot data for p1 and p2
plot_data_p1 <- generate_plot_data(p1)
plot_data_p2 <- generate_plot_data(p2)

# Define line colors with new labels
#line_colors <- c(
#  "EC_Ndah" = "",
#  "MS" = "magenta",
#  "Nakahigashi" = "#4f83cc",
#  "Stringer" = "#8da0cb",
#  "VanOrsdel" = "green",
#  "Weaver" = "#08306b",
#  "annotated" = "darkmagenta",
#  "control" = "grey"
#)

line_colors <- c(
  "annotated_conserved" = "darkmagenta",
  "annotated_ORFan" = "magenta",
  "novel_conserved" = "#08306b",
  "novel_ORFan" = "#8da0cb",
  "control" = "grey"
)

# Create ggplot cumulative distribution plot for p1
plot_p1 <- ggplot(plot_data_p1 %>% filter(Line!="MS") %>% filter(Line!="EC_Ndah") %>% filter(Line!="VanOrsdel"), aes(x = Cutoff, y = Fraction_Passed, color = Line)) +
  geom_line(size = 1) +
  labs(
    x = "Number of unique conditions in which genes are expressed",
    y = "Fraction of genes in each category"
  ) +
  scale_x_continuous(breaks = seq(0, 42, by = 3)) +
  scale_color_manual(values = line_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20),  # Adjust the size as needed
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24), # Adjust the size as needed
    axis.title.y = element_text(size = 24), # Adjust the size as needed
    legend.text = element_text(size = 12),  # Adjust the size as needed
    legend.title = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks = element_line(color = "black", size = 0.7)
  )

# Create ggplot cumulative distribution plot for p2
plot_p2 <- ggplot(plot_data_p2 %>% filter(Line!="MS") %>% filter(Line!="EC_Ndah") %>% filter(Line!="VanOrsdel"), aes(x = Cutoff, y = Fraction_Passed, color = Line)) +
  geom_line(size = 1) +
  labs(
    x = "Number of unique conditions in which genes are expressed",
    y = "Fraction of genes in each category"
  ) +
  scale_x_continuous(breaks = seq(0, 42, by = 3)) +
  scale_color_manual(values = line_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20),  # Adjust the size as needed
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24), # Adjust the size as needed
    axis.title.y = element_text(size = 24), # Adjust the size as needed
    legend.position = "none",  # Remove the legend from the second plot
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks = element_line(color = "black", size = 0.7)
  )

# Combine the two plots side-by-side with a shared legend
cumdist_combined_plot <- plot_p1 + plot_p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Print the combined plot
print(cumdist_combined_plot)

ggsave("/stor/scratch/Ochman/hassan/112724_protogene_extension/consolidating_datasets/Expression_cumdist.pdf", cumdist_combined_plot, width = 18, height = 10, units = "in")

###TPM count###

library(ggplot2)
library(dplyr)

setwd("/stor/work/Ochman/hassan/protogene_extension/expression_location_properties")

p1 <- read.csv("Ecoli.REL606.tpm.numbers.conservation",sep='\t')

dataset_order <- c("annotated_conserved", "annotated_ORFan", "novel_conserved", "novel_ORFan")

# Normalize counts relative to the 'total' values
p1_normalized <- p1 %>%
  group_by(dataset) %>%
  mutate(fraction = ifelse(category != "total", count / count[category == "total"], NA)) %>%
  filter(category %in% c("1tpm", "5tpm"))  # Exclude the "total" rows from plotting

p1_normalized$dataset <- factor(p1_normalized$dataset, levels = dataset_order)

# Create grouped bar plot
ggplot(p1_normalized, aes(x = dataset, y = fraction, fill = category)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  labs(x = "Dataset", y = "Fraction", fill = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability


ggplot(p1_normalized, aes(x = dataset, y = fraction, fill = category)) +
  geom_col(data = p1_normalized %>% filter(category == "1tpm"), position = "identity", fill = "steelblue") +  # Full 1tpm bars
  geom_col(data = p1_normalized %>% filter(category == "5tpm"), position = "identity", fill = "orange") +  # Smaller 5tpm bars
  theme_minimal() +
  labs(x = "Dataset", y = "Fraction", fill = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

ggplot(p1_normalized, aes(x = dataset, y = fraction, fill = category)) +
  geom_col(data = p1_normalized %>% filter(category == "1tpm"), position = "identity", 
           fill = "lightblue", color = "black", width = 0.7) +  # Light blue with black border
  geom_col(data = p1_normalized %>% filter(category == "5tpm"), position = "identity", 
           fill = "dodgerblue4", color = "black", width = 0.7) +  # Darker blue with black border
  theme_minimal(base_size = 14) +
  labs(x = "Dataset", y = "Fraction", fill = "Category") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Box around the plot
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.ticks = element_line(color = "black"),  # Add tick marks
    axis.line.x = element_line(color = "black", size = 0.25),  # Line along the base of bars
    panel.grid.major.y = element_line(color = "gray80", linetype = "solid"),  # Horizontal gridlines
    panel.grid.minor = element_blank()  # Keep minor gridlines off
  )

ggsave("/stor/scratch/Ochman/hassan/112724_protogene_extension/consolidating_datasets/Expression_TPMcount.pdf", width = 12, height = 8, units = "in")

#MeanTPM violin plot

library(ggplot2)
library(readr)
library(ggsignif)  # For significance bars

# Load the data
df <- read_tsv("Ecoli.REL606.tpm.conservation.meanvalues.tsv", col_names = TRUE)

# Convert category to a factor for proper ordering
df$category <- factor(df$category, levels = c("annotated_conserved", "annotated_ORFan", "novel_conserved","novel_ORFan"))

# Define custom colors
custom_colors <- c("annotated_conserved" = "darkmagenta",
                   "annotated_ORFan" = "salmon",
                   "novel_ORFan" = "lightblue",
                   "novel_conserved" = "blue")

# Define pairwise comparisons
comparisons <- list(
  c("annotated_ORFan", "annotated_conserved"),
  c("annotated_ORFan", "novel_ORFan"),
  c("annotated_ORFan", "novel_conserved"),
  c("annotated_conserved", "novel_ORFan"),
  c("annotated_conserved", "novel_conserved"),
  c("novel_ORFan", "novel_conserved")
)

# Create the violin plot
ggplot(df, aes(x = category, y = log(meantpm + 1), fill = category)) +
  geom_violin(trim = FALSE, alpha = 0.8) +  # Violin plot without trimming
  scale_fill_manual(values = custom_colors) +  # Set custom colors
  geom_signif(comparisons = comparisons, step_increase = 0.1, tip_length = 0.02) +  # Add significance bars
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),  # Larger Y-axis label
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1),  # Larger X-axis text
        axis.text.y = element_text(size = 16),  # Larger Y-axis text
        axis.ticks.y = element_line(),  # Add Y-axis ticks
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),  # Box around plot
        legend.position = "none") +
  ylab("Log(Mean TPM + 1)")


# Save the plot (optional)
ggsave("violin_plot.png", width = 6, height = 4, dpi = 300)
