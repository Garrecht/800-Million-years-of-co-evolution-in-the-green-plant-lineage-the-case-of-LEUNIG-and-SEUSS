# ---- Load libraries ----
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)

# ---- Input ----
input_file <- "elm_data2.txt"

# ---- Read data ----
df <- read_tsv(input_file, col_types = cols(
  gene_name = col_character(),
  Elm_name  = col_character(),
  position  = col_character(),
  sequence  = col_character()
))

# ---- Clean data ----
df <- df %>%
  mutate(
    gene_name = trimws(gene_name),
    Elm_name  = trimws(Elm_name),
    position  = trimws(position)
  )

# =========================================================
# 1. Extract all unique Elm names
# =========================================================

all_elms <- df %>%
  distinct(Elm_name)

write_tsv(all_elms, "all_elm_names.tsv")

# =========================================================
# 2. Find Elm names shared by ALL genes
# =========================================================

n_genes <- df %>%
  distinct(gene_name) %>%
  nrow()

shared_elms <- df %>%
  distinct(gene_name, Elm_name) %>%
  group_by(Elm_name) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n == n_genes) %>%
  select(Elm_name)

write_tsv(shared_elms, "shared_elm_names.tsv")

print(shared_elms)

# =========================================================
# 3. Count occurrences per gene
# =========================================================

counts <- df %>%
  group_by(gene_name, Elm_name) %>%
  summarise(count = n(), .groups = "drop")

write_tsv(counts, "elm_counts_per_gene.tsv")

# =========================================================
# 4. Extract start and end positions
# =========================================================

df_pos <- df %>%
  mutate(
    start = as.integer(str_extract(position, "^[0-9]+")),
    end   = as.integer(str_extract(position, "(?<=-)[0-9]+"))
  ) %>%
  filter(!is.na(start) & !is.na(end))

write_tsv(
  df_pos %>% select(gene_name, Elm_name, start, end),
  "elm_positions.tsv"
)
# =========================================================
# 5. Plot motif occupancy (dot plot)
# =========================================================

# Expand positions
df_expanded <- df_pos %>%
  rowwise() %>%
  mutate(pos = list(seq(start, end))) %>%
  unnest(pos)

# Order Elm names by frequency (optional)
elm_order <- df %>%
  count(Elm_name, sort = TRUE) %>%
  pull(Elm_name)

df_expanded$Elm_name <- factor(df_expanded$Elm_name, levels = elm_order)

n_genes <- df %>% distinct(gene_name) %>% nrow()

shared_elms <- df %>%
  distinct(gene_name, Elm_name) %>%
  group_by(Elm_name) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n == n_genes) %>%
  pull(Elm_name)

# ---- Mark shared motifs in the expanded data ----
df_expanded <- df_pos %>%
  rowwise() %>%
  mutate(pos = list(seq(start, end))) %>%
  unnest(pos) %>%
  mutate(shared = if_else(Elm_name %in% shared_elms, "Shared", "Other"))

# ---- Order Elm names (optional) ----
elm_order <- df %>%
  count(Elm_name, sort = TRUE) %>%
  pull(Elm_name)

df_expanded$Elm_name <- factor(df_expanded$Elm_name, levels = elm_order)

# ---- Plot ----
pos_plot <- ggplot(df_expanded, aes(x = Elm_name, y = pos, color = shared)) +
  geom_point(alpha = 0.7, size = 1.8) +
  facet_wrap(~ gene_name, ncol = 1, scales = "free_y", strip.position = "right") +
  theme_minimal() +
  labs(
    title = "ELM motif occupancy across protein sequences",
    x = "ELM motif",
    y = "Position (amino acid)",
    color = "Motif type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold"),
    legend.position = "none" 
  ) +
  scale_color_manual(values = c("Shared" = "steelblue", "Other" = "darkgrey"))

print(pos_plot)






# ---- Identify motifs shared by all genes ----
n_genes <- df %>% distinct(gene_name) %>% nrow()

shared_all <- df %>%
  distinct(gene_name, Elm_name) %>%
  group_by(Elm_name) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n == n_genes) %>%
  pull(Elm_name)

# ---- Filter expanded data to only shared motifs ----
df_shared <- df_pos %>%
  rowwise() %>%
  mutate(pos = list(seq(start, end))) %>%
  unnest(pos) %>%
  filter(Elm_name %in% shared_all)

# ---- Order Elm_names (optional) ----
df_shared$Elm_name <- factor(df_shared$Elm_name, levels = sort(unique(df_shared$Elm_name)))

my_palette <- c(
  "MOD_PIKK_1" = "#648FFF",
  "MOD_PKA_2" = "#785EF0",
  "MOD_Plk_1" = "#DC267F",
  "MOD_ProDKin_1" = "#FE6100",
  "TRG_ER_diArg_1" = "#FFB000"
  # add more Elm_names as needed
)

gene_order <- c("AtSEU", "Zm00001d049125", "Amtr071_SEU", "CYCAS_010520",
                "AtSLK2","Zm00001d035990", "AmTr069_SLK", "CYCAS_008694",
                "Ceric414", "Ceric673", "Ceric23s", "MaPoSEU1", "Knit",
                "CYCAS_024791", "MaPoSEU2" )  # replace with your gene names in order

# ---- Convert gene_name to factor with this order ----
df_shared$gene_name <- factor(df_shared$gene_name, levels = gene_order)

# ---- Plot with custom colors ----
pos_plot <- ggplot(df_shared, aes(y = Elm_name, x = pos, color = Elm_name)) +
  geom_point(alpha = 0.8, size = 2) +
  facet_wrap(~ gene_name, ncol = 1, scales = "free_y", strip.position = "right") +
  theme_minimal() +
  labs(
    title = "Motifs shared by all genes",
    y = "ELM motif",
    x = "Position (amino acid)",
    color = "ELM motif"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold"),
    legend.position = "none"
  ) +
  scale_color_manual(values = my_palette)

print(pos_plot)
