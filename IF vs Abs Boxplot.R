# Load required libraries
library(ggplot2)
library(ggdist)
library(dplyr)
library(glue)
library(ggtext)

# Define the color palette and fonts
okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                       "#0072B2", "#D55E00", "#CC79A7", "#999999")
bg_color <- "grey97"
font_family <- "Fira Sans"

# Example Absorbance and Impact Factor Data
absorbance_data <- list(
  Nature = list(A260_A280 = c(2.00, 1.95, 2.10, 1.99, 1.97), A260_A230 = c(2.20, 2.15, 2.18, 2.21, 2.19)),
  Science = list(A260_A280 = c(2.00, 1.95, 2.10, 1.99, 1.97), A260_A230 = c(2.19, 2.14, 2.22, 2.20, 2.18)),
  "Nature Biotechnology" = list(A260_A280 = c(2.02, 1.93, 1.97, 2.00, 1.98), A260_A230 = c(2.12, 2.10, 2.18, 2.15, 2.17)),
  "Genome Research" = list(A260_A280 = c(1.94, 1.92, 2.01, 1.96, 1.98), A260_A230 = c(2.15, 2.12, 2.16, 2.14, 2.19)),
  Lancet = list(A260_A280 = c(2.10, 2.05, 2.15, 2.11, 2.08), A260_A230 = c(2.25, 2.20, 2.23, 2.26, 2.24)),
  Cell = list(A260_A280 = c(2.05, 2.01, 2.12, 2.10, 2.06), A260_A230 = c(2.22, 2.19, 2.21, 2.23, 2.20)),
  "New England Journal of Medicine" = list(A260_A280 = c(2.11, 2.08, 2.14, 2.12, 2.10), A260_A230 = c(2.30, 2.25, 2.28, 2.31, 2.29)),
  "Nature Reviews Molecular Cell Biology" = list(A260_A280 = c(2.12, 2.09, 2.18, 2.15, 2.13), A260_A230 = c(2.32, 2.28, 2.30, 2.34, 2.31)),
  "Nature Reviews Drug Discovery" = list(A260_A280 = c(2.13, 2.10, 2.16, 2.14, 2.12), A260_A230 = c(2.33, 2.29, 2.31, 2.35, 2.32)),
  "Nature Reviews Cancer" = list(A260_A280 = c(2.14, 2.11, 2.17, 2.15, 2.13), A260_A230 = c(2.34, 2.30, 2.32, 2.36, 2.33)),
  "Nature Energy" = list(A260_A280 = c(2.09, 2.06, 2.14, 2.12, 2.10), A260_A230 = c(2.29, 2.26, 2.28, 2.32, 2.30)),
  "Nature Materials" = list(A260_A280 = c(2.07, 2.03, 2.11, 2.09, 2.05), A260_A230 = c(2.27, 2.24, 2.26, 2.30, 2.28)),
  "Nature Medicine" = list(A260_A280 = c(2.08, 2.04, 2.13, 2.11, 2.07), A260_A230 = c(2.28, 2.25, 2.27, 2.31, 2.29)),
  "Nature Communications" = list(A260_A280 = c(2.00, 1.98, 2.05, 2.03, 2.01), A260_A230 = c(2.20, 2.18, 2.22, 2.24, 2.21)),
  "Journal of the American Medical Association" = list(A260_A280 = c(2.06, 2.02, 2.10, 2.08, 2.04), A260_A230 = c(2.26, 2.22, 2.25, 2.29, 2.27)),
  "Chemical Reviews" = list(A260_A280 = c(2.04, 2.00, 2.08, 2.06, 2.02), A260_A230 = c(2.24, 2.20, 2.23, 2.27, 2.25)),
  "Nature Methods" = list(A260_A280 = c(2.03, 2.00, 2.07, 2.05, 2.01), A260_A230 = c(2.23, 2.19, 2.22, 2.26, 2.24)),
  "Nature Genetics" = list(A260_A280 = c(2.02, 1.99, 2.06, 2.04, 2.00), A260_A230 = c(2.22, 2.18, 2.21, 2.25, 2.23)),
  "Nature Physics" = list(A260_A280 = c(2.01, 1.97, 2.05, 2.03, 1.99), A260_A230 = c(2.21, 2.17, 2.20, 2.24, 2.22)),
  "Nature Chemistry" = list(A260_A280 = c(2.00, 1.96, 2.04, 2.02, 1.98), A260_A230 = c(2.20, 2.16, 2.19, 2.23, 2.21)),
  "Nature Immunology" = list(A260_A280 = c(1.99, 1.95, 2.03, 2.01, 1.97), A260_A230 = c(2.19, 2.15, 2.18, 2.22, 2.20))
)


impact_factors <- data.frame(
  Journal = c("Nature", "Science", "Nature Biotechnology", "Lancet", "Cell",
              "New England Journal of Medicine", "Nature Reviews Molecular Cell Biology",
              "Nature Reviews Drug Discovery", "Nature Reviews Cancer", "Nature Energy",
              "Nature Materials", "Nature Medicine", "Nature Communications",
              "Journal of the American Medical Association", "Chemical Reviews",
              "Nature Methods", "Nature Genetics", "Nature Physics", "Nature Chemistry",
              "Nature Immunology", "Genome Research"),
  '2019' = c(42.78, 41.84, 36.56, 60.39, 38.64, 74.70, 94.44, 84.70, 60.72, 60.86,
           43.84, 36.13, 12.12, 45.54, 52.76, 47.99, 38.33, 19.68, 30.36, 25.60, 10.78),
  '2020' = c(49.96, 47.73, 31.86, 79.32, 41.58, 91.25, 94.44, 84.70, 60.72, 60.86,
           43.84, 36.13, 14.92, 56.27, 54.30, 47.99, 38.33, 19.68, 30.36, 25.60, 8.27),
  '2021' = c(49.96, 47.73, 54.90, 79.32, 41.58, 91.25, 94.44, 84.70, 60.72, 60.86,
           43.84, 36.13, 14.92, 56.27, 54.30, 47.99, 38.33, 19.68, 30.36, 25.60, 8.96),
  '2022' = c(69.50, 63.70, 54.90, 202.73, 66.85, 176.08, 94.44, 84.70, 60.72, 60.86,
           43.84, 87.42, 17.69, 157.34, 72.09, 47.99, 38.33, 19.68, 30.36, 25.60, 6.70),
  '2023' = c(69.50, 63.70, 54.90, 202.73, 66.85, 176.08, 94.44, 84.70, 60.72, 60.86,
           43.84, 87.42, 17.69, 157.34, 72.09, 47.99, 38.33, 19.68, 30.36, 25.60, 6.13)
)

# Prepare Absorbance Data
long_data <- do.call(rbind, lapply(names(absorbance_data), function(journal) {
  data.frame(
    Journal = journal,
    Absorbance_Type = rep(c("A260/A280", "A260/A230"), each = 500),
    Absorbance = c(
      absorbance_data[[journal]]$A260_A280,
      absorbance_data[[journal]]$A260_A230
    )
  )
}))

colnames(impact_factors) <- sub("^X", "", colnames(impact_factors))

library(tidyr)

# Reshape Impact Factor Data
impact_factors_long <- impact_factors %>%
  pivot_longer(cols = 2019:2023, names_to = "Year", values_to = "Impact_Factor") %>%
  mutate(Year = as.numeric(Year))

library(tidyr)
library(glue)

# Normalize Impact Factor for Barplots
impact_factors$Impact_Factor_Normalized <- impact_factors$`2023` / max(impact_factors$`2023`) * 2.5

# Subtitle for the Plot
plot_subtitle <- glue(
  "Combined visualization of absorbance ratios (A260/A280, A260/A230) and 2023 Impact Factors (IF) for journals."
)

# Reorder Journals by 2023 Impact Factor (Highest to Lowest)
long_data$Journal <- factor(
  long_data$Journal,
  levels = impact_factors$Journal[order(-impact_factors$`2023`)]
)


# Reorder impact_factors data for barplots
impact_factors <- impact_factors %>%
  arrange(desc('2023'))

library(tidyr)

impact_factors_long <- impact_factors %>%
  pivot_longer(
    cols = `2019`:`2023`, # Use backticks for numeric column names
    names_to = "Year",
    values_to = "Impact_Factor"
  ) %>%
  mutate(Year = as.numeric(Year))




# Log-transform 2023 Impact Factors for normalization
impact_factors$Log_Impact_Factor <- log10(impact_factors$'2023')
impact_factors$Impact_Factor_Normalized <- impact_factors$Log_Impact_Factor / max(impact_factors$Log_Impact_Factor) * 2.5
########
library(ggplot2)
library(patchwork)
library(dplyr)

# Reorder Journals by 2023 Impact Factor
impact_factors <- impact_factors %>%
  mutate(Journal = factor(Journal, levels = Journal[order(-`2023`)]))

long_data <- long_data %>%
  mutate(Journal = factor(Journal, levels = impact_factors$Journal))

# Half-Eye Plot for Absorbance Ratios
half_eye_plot <- ggplot(long_data, aes(x = Absorbance, y = Journal, fill = Absorbance_Type)) +
  stat_halfeye(aes(color = Absorbance_Type), width = 0.6, point_size = 1.5, alpha = 0.4) +
  scale_fill_manual(values = okabe_ito_palette[1:2], name = "Absorbance Type") +
  scale_color_manual(values = c("A260/A280" = "blue", "A260/A230" = "green"), name = "Absorbance Type") +
  labs(
    title = "Absorbance Ratios (A260/A280, A260/A230)",
    x = "Absorbance Ratio",
    y = "Journal"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(hjust = 0),
    plot.title = element_text(size = 14)
  )

# Horizontal Bar Plot for Impact Factors
bar_plot <- ggplot(impact_factors, aes(x = Log_Impact_Factor, y = Journal)) +
  geom_bar(stat = "identity", fill = "grey60", color = "black", alpha = 0.8, width = 0.6) +
  scale_x_continuous(name = "Impact Factor (Log-Scaled)", expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "2023 Impact Factors",
    x = "Log Impact Factor",
    y = NULL  # Keep y-axis consistent with the other plot
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),  # Hide y-axis labels to avoid duplication
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 14)
  )

# Combine Plots Using Patchwork
combined_plot <- half_eye_plot | bar_plot +
  plot_annotation(
    title = "Combined Visualization of Absorbance Ratios and Impact Factors",
    subtitle = "Left: Absorbance Ratios (A260/A280, A260/A230). Right: Log-Scaled Impact Factors for 2023.",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, margin = margin(b = 10))
    )
  )

# Print the Combined Plot
print(combined_plot)

