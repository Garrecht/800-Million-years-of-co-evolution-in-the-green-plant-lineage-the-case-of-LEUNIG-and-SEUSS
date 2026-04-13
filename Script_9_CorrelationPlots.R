# simple correlation matrices
library(ComplexHeatmap)
library(grid)

tpm_file   <- read.csv("Marchantia_Expression.csv", header = T, sep = ";") 
top_correl <- "high_correlation_pairs_marchantia.tsv"

#transform data 
rownames(tpm_file) <- tpm_file$Name
tpm_file$Name <- NULL
head(tpm_file)
head(t(tpm_file))

#correlate
gene_cor <- cor(t(tpm_file), use = "pairwise.complete.obs", method = "pearson")


#make list of top correlated genes, and write table of pairs correlating above threshold
threshold <- 0.9

upper_idx <- upper.tri(gene_cor, diag = FALSE)

high_cor_pairs <- data.frame(
  Gene1 = rownames(gene_cor)[row(gene_cor)[upper_idx]],
  Gene2 = colnames(gene_cor)[col(gene_cor)[upper_idx]],
  Correlation = gene_cor[upper_idx]
)

high_cor_pairs <- high_cor_pairs[high_cor_pairs$Correlation > threshold, ]

head(high_cor_pairs)

write.table(
  high_cor_pairs,
  file = top_correl,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# Create an empty character matrix of same dimensions, for white middle line
annot_mat <- matrix("", nrow = nrow(gene_cor), ncol = ncol(gene_cor))
rownames(annot_mat) <- rownames(gene_cor)
colnames(annot_mat) <- colnames(gene_cor)

#make heatmap
Heatmap(gene_cor,
        name = "Pearson correlation",
        
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        
        row_dend_width = unit(1.5, "cm"),
        row_names_side = "left",
        
        row_names_max_width = unit(12.5, "cm"),
        
        column_names_side = "top",
        column_dend_side = "top",
        
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 7),
        
        column_dend_height = unit(1.5, "cm"),
        
        col = colorRampPalette(c("darkorange", "#F0F0F0", "#1BAEA7"))(100),
        
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 6), 
          labels_gp = gpar(fontsize = 6),  
          legend_height = unit(2.5, "cm"),   
          grid_width = unit(1.5, "mm")      
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(i == j) {
            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = "white", col = "grey"))
          }
        }
)
