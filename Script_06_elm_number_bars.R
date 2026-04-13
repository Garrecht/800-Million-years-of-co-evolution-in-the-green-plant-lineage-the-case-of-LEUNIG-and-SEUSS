# Load libraries
library(dplyr)
library(tidyr)
library(readr)

# Read the two TSV files
df_counts <- read_tsv("Elm_counts.txt")   # gene_name, Elm_name, count
df_elm    <- read_tsv("Elm_descr.txt")   # Elm_name, interaction_domain, functional_site_class

# Inspect structure (optional)
# head(df_counts)
# head(df_elm)

# Convert counts table to wide format (one column per gene)
df_wide <- df_counts %>%
  pivot_wider(
    names_from = gene_name,
    values_from = count,
    values_fill = list(count = 0)
  )

# Join with Elm annotation table
df_final <- df_elm %>%
  left_join(df_wide, by = "Elm_name")

# Write output
write_tsv(df_final, "expanded_output.tsv")

# Print result
print(df_final)





library(tidyverse)

# Read data
df <- read_tsv("Elm_matrix.txt")

# Reshape (wide → long)
df_long <- df %>%
  pivot_longer(
    cols = -Elm_name,
    names_to = "gene_name",
    values_to = "count"
  ) %>%
  filter(count > 0)   # optional: remove zeros for cleaner plots

# Create faceted bar plots (one per gene)
p <- ggplot(df_long, aes(x = Elm_name, y = count, fill = Elm_name)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ gene_name, ncol = 1, scales = "free_y") +
  labs(
    title = "Elm Motif Counts per Gene",
    x = "Elm motif",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

# Show plot
print(p)