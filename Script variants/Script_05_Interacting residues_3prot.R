# Load necessary libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)
library(Biostrings)
library(writexl)

##### IMPORTANT ####
# The sequences in "sequences_R_3prot.fa" MUST be sorted like the alphabetically sorted .cif files. For each .cif file, three sequences are present, the one mentioned first as the first one, then the second one and the final one(see file names)
# The sequences need to be named with the sequence name, the name of the compared sequence, and 
# a number referencing the corresponding partner protein sequence!

# I.e
# Fold_AtLUG_AtHCA19_AtSEU
# -->
# >AtLUG_AtHDA19_p1
# [SEQUENCE]
# >AtLUG_AtSEU_p2
# [SEQUENCE]
# >AtHDA19_AtLUG_p1
# [SEQUENCE]
# AtHDA19_AtSEU_p3
# [SEQUENCE]
# >AtSEU_AtLUG_p2
# [SEQUENCE]
# >AtSEU_AtHDA19_p3
# [SEQUENCE]
# ...

# Even if sequences are present multiple times!
# Otherwise, the output of the script WILL be wrong!

##### Intersecting close AS and Protein sequences ####

# Import and read files
fasta_file_3p <- "sequences_R_3prot.fa" # File containing SORTED fasta sequences
interaction_file_3p_50 <- "../AlphaFold models/Interacting residues_short_3p_50.txt" # File containing interaction residues in short format (pIDDT50)
interaction_file_3p_70 <- "../AlphaFold models/Interacting residues_short_3p_70.txt" # File containing pIDDT >70 residues

interaction_lines_3p_50 <- readLines(interaction_file_3p_50)
interaction_lines_3p_70 <- readLines(interaction_file_3p_70)


# Use Biostirngs to read fasta file
fasta_lines_3p <- readAAStringSet(fasta_file_3p)
sequence_names_3p <- names(fasta_lines_3p)
sequences_3p <- as.character(fasta_lines_3p)


# Filter interaction lines to exclude lines with names (i.e.: Residues in...)
interaction_lines_3p_50 <- interaction_lines_3p_50[!grepl("^Residues", interaction_lines_3p_50)]
interaction_lines_3p_70 <- interaction_lines_3p_70[!grepl("^Residues", interaction_lines_3p_70)]

# Initialize an empty data frame to store all results
all_results_3p <- data.frame()

# Process each sequence and its interaction data
for (i in seq_along(sequences_3p)) {
   # Extract the sequence and its name
   sequence <- sequences_3p[i]
   sequence_name <- sequence_names_3p[i]

   # Extract interaction positions for this sequence
   interaction_positions_3p_50 <- as.numeric(unlist(strsplit(interaction_lines_3p_50[i], ",")))
   interaction_positions_3p_70 <- as.numeric(unlist(strsplit(interaction_lines_3p_70[i], ",")))

   # Split sequence into residues and into the protein number (third element of name)
   residues <- strsplit(sequence, "")[[1]]
   residue_positions <- 1:length(residues)
   ProtNr <- sapply(strsplit(sequence_name,"_"), `[`, 3)
   
   # Create a dataframe
   sequence_df <- data.frame(
      Sequence_Name = sequence_name,
      ProtNr = ProtNr,
      Position = residue_positions,
      Residue = residues
   )

   # Create interaction data frame
   closeAS_50 <- data.frame(Position_50 = interaction_positions_3p_50)
   closeAS_70 <- data.frame(Position_70 = interaction_positions_3p_70)

   # Mark interacting residues
   sequence_df <- sequence_df %>%
      mutate(Interacting = if_else(Position %in% closeAS_70$Position, "70",
         if_else(Position %in% closeAS_50$Position, "50", "0")
      ))

   # Append results to the combined data frame
   all_results_3p <- bind_rows(all_results_3p, sequence_df)
}


#### OPTIONAL: Create a shorter, sorted result list ####
# Filter desired interactions
#all_results_3p <- all_results_3p %>%
#   filter(grepl("_p1$|_p10$|_p14$|_p15$|_p18$", Sequence_Name))

# Get list of all filtered names and bring them in requested order
#levels(factor(all_results_3p$Sequence_Name))
#order <- c("AtLUG_p1", "AtSEU_p1", "CrLUG_p10", "CrSEU1_p10", "PpLUG1_p18", "PpSEU2_p18", "MpLUG_p15", "MpSEU1_p15", "KnLUG_p14", "KnSEU_p14")

#all_results_3p <- all_results_3p %>%
#   mutate(Sequence_Name = factor(Sequence_Name, levels = order)) %>% 
#   arrange(Sequence_Name)


#### Plotting of results ####

# Find maximum protein length
max_length <- max(all_results_3p$Position)

# Order results according to protein complex number, seuqnece name and residue position
all_results_3p <- all_results_3p %>% 
   arrange(ProtNr, Sequence_Name, Position)

# Initialize a list to store plots
plot_list_3p <- list()

# Loop through each unique sequence name
for (sequence_name in unique(all_results_3p$Sequence_Name)) {
   # Filter data for the current sequence
   sequence_data <- all_results_3p %>%
      filter(Sequence_Name == sequence_name)

   # Create the plot
   p <- ggplot(sequence_data, aes(x = Position, y = 1, fill = Interacting)) +
      geom_tile(height = 1) + # Create squares for each residue
      scale_fill_manual(
         values = c("70" = "#038f0a", "50" = "#19059c", "0" = "#d9d5d4"),
         name = "Residues",
         labels = c("70" = "≤ 4.5 Å & pIDDT > 70", "50" = "≤ 4.5 Å & 50 < pIDDT ≤ 70", "0" = "> 4.5 Å | pIDDT ≤ 50")
      ) +
      scale_x_continuous(limits = c(0, max_length)) + # Apply common x-axis scale
      labs(
         x = "Residue Position",
         y = sequence_name
      ) +
      theme_minimal() +
      theme(
         aspect.ratio = 0.02, # Compresses Y-axis
         axis.text.y = element_blank(), # Remove Y-axis text
         axis.ticks.y = element_blank(), # Remove Y-axis ticks
         axis.title.y = element_blank(), # Removes Y-axis title
         panel.grid.major = element_blank(), # Remove grid lines
         panel.grid.minor = element_blank(),
         plot.margin = margin(0, 0, 0, 0) # Remove margins around the plot
      )

   # Add the plot to the list
   plot_list_3p[[sequence_name]] <- p
}


# Remove X-axis labels for all but the last plot (to keep as legend)
plot_list_3p <- lapply(seq_along(plot_list_3p), function(i) {
   if (i < length(plot_list_3p)) {
      plot_list_3p[[i]] + theme(
         axis.title.x = element_blank(),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank()
      )
   } else {
      plot_list_3p[[i]] # Keep the X-axis for the last plot
   }
})


# Remove legends for plots without any residues lower than pIDDT 70 or 50. (For combining legend via patchwork necessary)
remove_short_legends <- function(plot, legend_column) {
   # Extract the data used in the legend
   data <- ggplot_build(plot)$data[[1]]

   # Count the number of unique entries in the specified column
   num_unique <- length(unique(data[[legend_column]]))

   # Remove the legend if fewer than 3 unique entries
   if (num_unique < 3) {
      plot <- plot + guides(color = "none", fill = "none")
   }

   return(plot)
}

# Apply the function to each plot in the list
plot_list_3p <- lapply(plot_list_3p, function(plot) {
   remove_short_legends(plot, "group")
})


# Combine plots in pairs  to signify complexes, stack them vertically, and add a combined title
combined_plots_3p <- lapply(seq(1, length(plot_list_3p), by = 2), function(i) {
   # Select two plots of a protein complex
   plot_pair <- plot_list_3p[c(i, i + 1)]

   # Remove the Y-axis titles for both plots in the pair
   plot_pair <- lapply(plot_pair, function(p) {
      p + theme(axis.title.y = element_blank())
   })

   # Get the y-axis titles from both plots and use as titles
   y_titles <- sapply(plot_pair, function(p) {
      title <- ggplot_build(p)$plot$labels$y
      # Remove unwanted suffixes
      gsub("_p[0-9]+", "", title)
   })

   # Combine the y-axis titles for the pair (e.g., "Title 1 + Title 2")
   y_title_1 <- sapply(strsplit(y_titles[1],"_"), `[`, 1)
   y_title_2 <- sapply(strsplit(y_titles[2],"_"), `[`, 1)
   
   combined_title <- paste(y_title_1, "+", y_title_2)
   # Add title to the first plot of the plot pair
   plot_pair <- lapply(seq_along(plot_pair), function(i) {
      if (i < length(plot_pair)) {
         plot_pair[[i]] + labs(title = combined_title)
      } else {
         plot_pair[[i]]
      }
   })


   # Combine the two plots vertically
   combined_plot <- wrap_plots(plot_pair, ncol = 1) +
      plot_layout(heights = c(0.5, 0.5))
   return(combined_plot)
})


# Combine all paired plots together in a single layout (all stacked vertically)
final_combined_plot <- wrap_plots(combined_plots_3p, ncol = 1) +
   plot_layout(guides = "collect")
# Save results
ggsave("Interacting_residues_full_3p.pdf", final_combined_plot, width = 10, height = 9, device = cairo_pdf)
ggsave("Interacting_residues_full_3p.svg", final_combined_plot, width = 10, height = 9)
#ggsave("Interacting_residues_short.pdf", final_combined_plot, width = 10, height = 5, device = cairo_pdf)
