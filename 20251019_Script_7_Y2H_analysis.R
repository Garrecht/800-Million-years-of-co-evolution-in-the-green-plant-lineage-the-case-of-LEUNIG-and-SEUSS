## Information & Settings ####
# Goal: Automated Processing (Gitter, rename and analyse) pictures from Y2H-screens

# Use the gitter.batch() function, which allows you to supply a reference screen that is analyzed first and used to create a grid that is overlayed onto the second plate (the plate you want to analyse)
# For this, the images need to be as similar as possible regarding pixel size, position of the plate ect.
# Crop the images to the INSIDE of the plate using an ImageJ script to keep everything consistent
# Check the gitted images visually to confirm that the grids were placed correctly! IMPORTANT!

# Also, after the gitter analysis, check the output images for unwanted artifacts or noise and make sure that all colonies were captured by the analysis
# THIS IS THE ONLY WAY TO ENSURE THAT THE DATA WAS COLLECTED PROPERLY!

# If colonies are not captured correctly, try the following:
# Re-gitter the image. Sometimes this alone gets rid of the noise
# Re-crop the images to make sure your positioning is ok
# Add a small bit of gausian blur (start with 1, slowly increase) in case there is  high noise in the gitter_picture plate
# If there are non-captured spots (black) within colonies, they are caused by the general color of the colonies/color differences within the colonies which is wrongfully assumed to be different by the binary transformation
# Try reducing the background in ImageJ using rolling ball substraction


# Gittered reference image that contains all 96 colonies in a good quality
file_ref <- "./Gitter/Plate_reference.jpg"

# Define save location of files
dat_save <- "./Gitter/Output_file/"
grid_save <- "./Gitter/Output_image/"

# Define your source directory
dir_source <- "./Gitter/Source/"

# Gitter reference file
dat_ref <- gitter(file_ref,
   dat.save = dat_save,
   grid.save = grid_save,
   plate.format = 96,
   verbose = "p"
)


# IMPROTANT:
# The script assumes that your files ar named after the following scheme:
# Plate_[ID]v[REPLICATE]_[Treatment]
# [ID] and [REPLICATE] should be numeric values
# [Treatment] should be one of the following
# LW_Transfer [This is the reference plate]
# LWH
# LWH_3-AT


# ANALYSIS:
# Idea: Calculate sizes of LWH and 3-AT colony growth in relation to the LW-transfer colonies
# --> The smaller the difference, the more likely that there is a genuine protein interaction, as this allows for growth of colonies on the selective media


# Definition of selective thresholds
colony_size_threshold <- 900 # Based on a first look, growing spots seem to have a colony size of >1000

# DEFINITION OF WHEN WE CONSIDER COLONY GRWOTH TO BE A GENUINE PROTEIN INTERACTION
# Currently: If assay colonies have grown to a size of 50% or more in relation to control colony
interaction_true <- 0.5


# LOOKUP-TABLE
# To automatically map the row/col coordinates of the plates to the proteins, you must create a lookup table in excel:
# One sheet consisting of seven columns with the following data:
# [Plate], [row], [row_tag], [row_protein], [column], [col_tag], [col_protein]
# [Plate] is plate number WITHOUT replicates (i.e 1, 2, 3 ect.)
# [row/col] is position of row/col
# [tag] is AD/BD
# [protein] is protein name

# Load lookup-table
library(readxl)
lookup_grid <- read_excel("./Y2H_Lookup_table.xlsx", sheet = "grid")

# Merge the tag and protein columns to one combined column and remove the old ones
lookup_grid <- lookup_grid %>%
   mutate(tag_prot_row = paste0(tag1, "_", protein1)) %>%
   select(!c(tag1, protein1))

lookup_grid <- lookup_grid %>%
   mutate(tag_prot_col = paste0(tag2, "_", protein2)) %>%
   select(!c(tag2, protein2))

# Change plate number from "Double" to "Character" and add "Plate_" for later merge with results
lookup_grid$plate <- as.character(lookup_grid$plate)
lookup_grid$plate <- paste0("Plate_", lookup_grid$plate)


## Setup ####

# Install necessary packages
install.packages("gplots")
install.packages("LSD")
install.packages("RColorBrewer")
install.packages("plyr")
install.packages("stringr")
install.packages("logging")

install.packages("devtools")
library(devtools)
devtools::install_github("omarwagih/gitter")

install.packages("BiocManager")
library(BiocManager)
BiocManager::install("EBImage")
BiocManager::install("ComplexHeatmap")

# Load necessary packages

library(gitter)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(writexl)

## Batch Gitter ####
gitter.batch(dir_source,
   ref.image.file = file_ref,
   dat.save = dat_save,
   grid.save = grid_save,
   plate.format = 96,
   verbose = "p"
)


# This code loads multiple tab-separated gitter data files, appends them row-wise, and adds columns indicating the source file and treatment for each row:
combined_gitter_data <- function(directory) {
   gitter_data_list <- list.files(dat_save, pattern = "\\.dat$", full.names = TRUE)
   gitter_data_inter <- list()

   # The col names are unfortunatly marked as comments, so we need to supply them manually
   colnames_gitter <- c("row", "col", "size", "circularity", "flags")

   for (file in gitter_data_list) {
      data_f <- read.table(file,
         comment.char = "#",
         sep = "\t"
      )
      colnames(data_f) <- colnames_gitter

      # Splits the file name into separet stings at each "_" in the name to precisely name the plates and treatments
      split_name <- strsplit(basename(file), "_")[[1]]
      data_f$plate <- paste(split_name[1], split_name[2], sep = "_")

      # An if loop that creates a treatment column and checks if the fourth string starts with "cropped"
      # Will be the case for LWH-plates, but not LW_Transfer or LWH_3AT
      if (length(split_name) >= 4) {
         if (startsWith(split_name[4], "cropped")) {
            data_f$treatment <- paste(split_name[3])
         } else {
            data_f$treatment <- paste(split_name[3], split_name[4], sep = "_")
         }
      } else {
         data_f$treatment <- paste("Unknown")
      }
      gitter_data_inter[[file]] <- data_f
   }
   gitter_data_inter <- do.call(rbind, gitter_data_inter)
   rownames(gitter_data_inter) <- 1:nrow(gitter_data_inter)
   return(gitter_data_inter)
}
gitter_data <- combined_gitter_data()


## Analysis of results ####

# Sets the names of the treatments to numerical values, Easier for referencing and working with them
gitter_data <- gitter_data %>% mutate(
   treatment =
      recode(treatment,
         "LW_Transfer" = 1,
         "LWH" = 2,
         "LWH_3AT" = 3,
         "Unknown" = 4
      )
)
gitter_data <- filter(gitter_data, treatment != 4)


#### Function to obtain and process data from singular plates ####
process_plate <- function(gitter_data) {
   #### #### Extracting data from singular plates and creating a three-dimensional matrix of position and treatment ####

   unique_rows <- sort(unique(gitter_data$row))
   unique_cols <- sort(unique(gitter_data$col))
   unique_treatments <- sort(unique(gitter_data$treatment))
   plate_matrix <- array(NA,
      dim = c(length(unique_rows), length(unique_cols), length(unique_treatments)),
      dimnames = list(unique_rows, unique_cols, unique_treatments)
   )

   for (treat in unique_treatments) {
      treatment_data <- gitter_data %>% filter(treatment == !!treat)
      for (i in 1:nrow(treatment_data)) {
         row_index <- which(unique_rows == treatment_data$row[i])
         col_index <- which(unique_cols == treatment_data$col[i])
         plate_matrix[row_index, col_index, treat] <- treatment_data$size[i]
      }
   }

   #### #### The analysis that is being done with the data ####

   # Sets the LW-Media as the base for the next analysis steps
   base_treatment <- unique_treatments[1]

   # Create a mask for none-growing colonies on the LW-transfer treatment (1)
   threshold_mask <- abs(plate_matrix[, , base_treatment]) < colony_size_threshold

   treatment_matrix <- array(NA, dim = dim(plate_matrix), dimnames = dimnames(plate_matrix))
   for (treat in unique_treatments) {
      if (treat != base_treatment) {
         treatment_matrix[, , treat] <- plate_matrix[, , treat] / plate_matrix[, , base_treatment]
         treatment_matrix[treatment_matrix > 1] <- 1
      }
   }

   for (treat in unique_treatments) {
      if (treat != base_treatment) {
         treatment_matrix[, , treat][threshold_mask] <- NA
      }
   }
   treatment_df <- as.data.frame(as.table(treatment_matrix))
   colnames(treatment_df) <- c("row", "col", "treatment", "rel.size")

   interaction_treat <- treatment_df %>%
      filter(rel.size >= interaction_true)

   list(
      "Gitter analysis" = plate_matrix,
      "Assay effects" = treatment_matrix,
      "No colony growth" = threshold_mask,
      "True interactions" = interaction_treat
   )
}


#### Function to loop the analysis through all plates ####
analyze_all_plates <- function(gitter_data) {
   plate_list <- split(gitter_data, gitter_data$plate)
   plate_list <- plate_list[sapply(plate_list, function(gitter_data) {
      any(gitter_data$treatment == "1")
   })]
   gitter_results <- list()

   # Process each plate by looping through each plate via their ID
   for (plate_id in names(plate_list)) {
      cat("Processing plate", plate_id, "\n")
      plate_analysis <- plate_list[[plate_id]]
      gitter_results[[plate_id]] <- process_plate(plate_analysis)
   }
   gitter_results
}
results <- analyze_all_plates(gitter_data)


#### Function to obtain a dataframe of all true interactions ####
get_protein_interaction <- function(results) {
   interaction_results <- data.frame()
   for (plate_id in names(results)) {
      plate <- results[[plate_id]]
      if ("True interactions" %in% names(plate)) {
         interaction_df <- plate$`True interactions`
      }
      if (nrow(interaction_df) == 0) {
         next
      }

      interaction_df$plate <- plate_id
      interaction_results <- bind_rows(interaction_results, interaction_df)
   }

   # Modify the interaction_results table to store the replicate number in a different column to allow for merging of the results with the lookup table names

   interaction_results <- interaction_results %>%
      mutate(replicate = gsub("Plate_[0-9](.*)", "\\1", plate)) %>%
      mutate(plate = gsub("v[0-9]+$", "\\1", plate))

   interaction_results$row <- as.numeric(interaction_results$row)
   interaction_results$col <- as.numeric(interaction_results$col)

   merged_results <- interaction_results %>%
      inner_join(final_lookup, by = c("plate", "row", "col")) %>%
      mutate(treatment = case_when(
         treatment == "2" ~ "LWH",
         treatment == "3" ~ "LWH_3-AT",
         TRUE ~ treatment
      )) %>%
      separate(tag_prot_row, into = c("tag_row", "prot_row"), sep = "_", extra = "merge") %>%
      separate(tag_prot_col, into = c("tag_col", "prot_col"), sep = "_", extra = "merge")
   return(merged_results)
}
Assay_results <- get_protein_interaction(results)


# Shorten results: How many times (i.e replicates and both directions of interactions) is a interaction present in the Assay_results

short_and_sort <- function(Assay_results) {
   Assay_results_empty <- Assay_results %>%
      filter(prot_row == "empty" | prot_col == "empty") %>%
      group_by(tag_row, prot_row, tag_col, prot_col) %>%
      summarize(count = n(), .groups = "drop")

   Assay_results_M <- Assay_results %>%
      mutate(tag_prot_row = paste0(tag_row, "_", prot_row)) %>%
      select(!c(tag_row, prot_row))

   Assay_results_M <- Assay_results_M %>%
      mutate(tag_prot_col = paste0(tag_col, "_", prot_col)) %>%
      select(!c(tag_col, prot_col))


   Assay_results_short <- Assay_results_M %>%
      rowwise() %>%
      mutate(
         Sorted_Names = list(c(treatment, sort(c(tag_prot_row, tag_prot_col))))
      ) %>%
      ungroup() %>%
      mutate(
         treatment = sapply(Sorted_Names, `[`, 1),
         protein_1 = sapply(Sorted_Names, `[`, 2),
         protein_2 = sapply(Sorted_Names, `[`, 3)
      ) %>%
      group_by(treatment, protein_1, protein_2) %>%
      summarize(
         count = n(), .groups = "drop"
      ) %>%
      separate(protein_1, into = c("tag_1", "prot_1"), sep = "_", extra = "merge") %>%
      separate(protein_2, into = c("tag_2", "prot_2"), sep = "_", extra = "merge")

   write_xlsx(
      list(
         Assay_results = Assay_results,
         Empty_interactions = Assay_results_empty,
         Assay_results_short = Assay_results_short
      ),
      "./Assay_results.xlsx"
   )
}

short_and_sort(Assay_results)


#### Visualizing Analysis Output ####

generate_bubble <- function(results) {
   # Recode the treatment names for better overview (for use in plot only!)
   treatment_labels <- c(
      "1" = "LW",
      "2" = "LWH",
      "3" = "LWH + 3-AT"
   )
   pdf(
      file = "./Y2H_results_plots.pdf",
      width = 8, height = 6
   )

   for (plate_id in names(results)) {
      plate_results <- results[[plate_id]]
      treatment_results <- plate_results$`Assay effects`
      treatment_name <- dimnames(treatment_results)[[3]]
      for (treatment in treatment_name) {
         # Remove analysis of LW-plate against LW-plate from dataframe (treatment = 1), as the output is just NAs
         bubble_matrix <- treatment_results[, , treatment]
         if (all(is.na(bubble_matrix))) {
            next
         }

         bubble_df <- as.data.frame(as.table(bubble_matrix))
         colnames(bubble_df) <- c("row", "col", "rel.size")


         # Replace NA values for missing colonies with -1, a value that is otherwise impossible
         bubble_df <- bubble_df %>%
            mutate(rel.size = if_else(is.na(rel.size), -1, rel.size))

         bubble_df <- bubble_df %>%
            mutate(threshold = if_else(rel.size < interaction_true, "FALSE", "TRUE"))

         print(
            plot <- ggplot(
               bubble_df,
               aes(x = col, y = row)
            ) +
               geom_point(shape = 21, fill = "darkgrey") +
               geom_point(
                  shape = 21,
                  aes(size = rel.size, fill = threshold)
               ) +
               scale_radius(
                  range = c(3, 15),
                  limits = c(0, 1),
                  breaks = c(0, 0.25, 0.5, 0.75, 1),
                  labels = c("NA", "0.25", "0.5", "0.75", "1")
               ) +
               scale_fill_manual(
                  values = c(
                     "FALSE" = "#F61E1E",
                     "TRUE" = "#4ADB2C"
                  ),
                  name = "Protein interaction"
               ) +
               labs(
                  title = paste("Bubble Plot for Treatment:", plate_id, ",", treatment_labels[treatment]),
                  x = "Column",
                  y = "Row"
               ) +
               guides(
                  fill = "none",
                  size = guide_legend(
                     title = "Rel. size",
                     override.aes = list(
                        size = c(1.5, 6, 9, 12, 15),
                        fill = c("darkgrey", "#F61E1E", "#4ADB2C", "#4ADB2C", "#4ADB2C")
                     )
                  )
               ) +
               scale_y_discrete(limits = rev) +
               theme_minimal() +
               theme(
                  legend.position = "right",
                  panel.grid.major = element_blank(),
                  panel.border = element_rect(fill = NA),
                  axis.ticks = element_line(color = "black")
               ) +
               coord_fixed()
         )

         cat("Processing plate", plate_id, ",", treatment_labels[treatment], "\n")
      }
   }
   dev.off()
}

generate_bubble(results)
