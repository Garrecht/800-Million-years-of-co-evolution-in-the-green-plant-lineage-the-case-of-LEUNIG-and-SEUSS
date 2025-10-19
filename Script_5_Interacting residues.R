# Load necessary libraries
library(ggplot2)
library(dplyr)
library(patchwork)

##### IMPORTANT ####
# The sequences in "Sequences_R.fa" MUST be sorted like the alphabetically sorted .cif files. For each .cif file, two sequences are present, the one mentioned first as the first one, then the second one (see file names)

##### Intersecting close AS and Protein Sequences ####

# Import and read files
fasta_file <- "Sequences_R.fa" # File containing SORTED fasta sequences
interaction_file_50 <- "Interacting residues_short_50.txt" # File containing interaction residues in short format (pIDDT50)
interaction_file_70 <- "Interacting residues_short_70.txt" # File containing pIDDT >70 residues

fasta_lines <- readLines(fasta_file)
interaction_lines_50 <- readLines(interaction_file_50)
interaction_lines_70 <- readLines(interaction_file_70)

# Filter sequences and names from fasta file and removes ">"
sequence_names <- fasta_lines[seq(1, length(fasta_lines), by = 2)]
sequences <- fasta_lines[seq(2, length(fasta_lines), by = 2)]
sequence_names <- gsub("^>", "", sequence_names)

# Filter interaction lines to exclude lines with names (i.e.: Residues in...)
interaction_lines_50 <- interaction_lines_50[!grepl("^Residues", interaction_lines_50)]
interaction_lines_70 <- interaction_lines_70[!grepl("^Residues", interaction_lines_70)]


all_results <- data.frame()

# Process each sequence and its interaction data
for (i in seq_along(sequences)) {
   sequence <- sequences[i]
   sequence_name <- sequence_names[i]

   interaction_positions_50 <- as.numeric(unlist(strsplit(interaction_lines_50[i], ",")))
   interaction_positions_70 <- as.numeric(unlist(strsplit(interaction_lines_70[i], ",")))

   residues <- strsplit(sequence, "")[[1]]
   residue_positions <- 1:length(residues)
   sequence_df <- data.frame(
      Sequence_Name = sequence_name,
      Position = residue_positions,
      Residue = residues
   )

   closeAS_50 <- data.frame(Position_50 = interaction_positions_50)
   closeAS_70 <- data.frame(Position_70 = interaction_positions_70)

   sequence_df <- sequence_df %>%
      mutate(Interacting = if_else(Position %in% closeAS_70$Position, "70",
         if_else(Position %in% closeAS_50$Position, "50", "0")
      ))

   all_results <- bind_rows(all_results, sequence_df)
}

#### Plotting of results: ####
max_length <- max(all_results$Position)
plot_list <- list()

# Loop through each unique sequence name
for (sequence_name in unique(all_results$Sequence_Name)) {
   sequence_data <- all_results %>%
      filter(Sequence_Name == sequence_name)

   p <- ggplot(sequence_data, aes(x = Position, y = 1, fill = Interacting)) +
      geom_tile(height = 1) +
      scale_fill_manual(
         values = c("70" = "#038f0a", "50" = "#19059c", "0" = "#d9d5d4"),
         name = "Residues",
         labels = c("70" = "≤ 4.5 Å & pIDDT > 70", "50" = "≤ 4.5 Å & 50 < pIDDT ≤ 70", "0" = "> 4.5 Å | pIDDT ≤ 50")
      ) +
      scale_x_continuous(limits = c(0, max_length)) +
      labs(
         x = "Residue Position",
         y = sequence_name
      ) +
      theme_minimal() +
      theme(
         aspect.ratio = 0.02,
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.title.y = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.margin = margin(0, 0, 0, 0)
      )

   plot_list[[sequence_name]] <- p
}

# Remove X-axis labels for all but the last plot (to keep as legend)
plot_list <- lapply(seq_along(plot_list), function(i) {
   if (i < length(plot_list)) {
      plot_list[[i]] + theme(
         axis.title.x = element_blank(),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank()
      )
   } else {
      plot_list[[i]]
   }
})

# Remove legends for plots without any residues lower than pIDDT 70 or 50. (For combining legend via patchwork)
remove_short_legends <- function(plot, legend_column) {
   # Extract the data used in the legend
   data <- ggplot_build(plot)$data[[1]]

   # Count the number of unique entries in the specified column
   num_unique <- length(unique(data[[legend_column]]))

   if (num_unique < 3) {
      plot <- plot + guides(color = "none", fill = "none")
   }

   return(plot)
}
plot_list <- lapply(plot_list, function(plot) {
   remove_short_legends(plot, "group")
})

# Combine plots in pairs to signify interactions, stack them vertically, and add a combined title
combined_plots <- lapply(seq(1, length(plot_list), by = 2), function(i) {
   plot_pair <- plot_list[c(i, i + 1)]

   plot_pair <- lapply(plot_pair, function(p) {
      p + theme(axis.title.y = element_blank())
   })

   y_titles <- sapply(plot_pair, function(p) {
      title <- ggplot_build(p)$plot$labels$y
      # Remove unwanted suffixes
      gsub("_p[0-9]+", "", title)
   })
   combined_title <- paste(y_titles[1], "+", y_titles[2])

   # Add title to the first plot of the plot pair
   plot_pair <- lapply(seq_along(plot_pair), function(i) {
      if (i < length(plot_pair)) {
         plot_pair[[i]] + labs(title = combined_title)
      } else {
         plot_pair[[i]]
      }
   })

   combined_plot <- wrap_plots(plot_pair, ncol = 1)
   return(combined_plot)
})

# Combine all paired plots together in a single layout (all stacked vertically)
final_combined_plot <- wrap_plots(combined_plots, ncol = 1) +
   plot_layout(guides = "collect")

ggsave("Interacting_residues.pdf", final_combined_plot, width = 10, height = 15)
