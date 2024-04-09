# Follow me on computingbiology.blog
#
# Companion code for 

##### LIBRARIES #####
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(biomaRt)
library(oranges)



##### GENE SET OVERLAPS #####
glyc <- grep("glycolysis", oranges::cpdb_data$pathway_info$pathway_name)

O <- matrix(0, nrow = length(glyc), ncol = length(glyc))
glyco_comparison <- vector(mode = "list", length = length(glyc))

for (i in seq(1, length(glyc))) {
  
  gcount <- length(names(which(oranges::cpdb_data$pathway_matrix[glyc[i], ] != 0)))
  psource <- oranges::cpdb_data$pathway_info$pathway_source[glyc[i]]
  
  for (j in seq(1, length(glyc))) {
    
    aub <- length(union(names(which(oranges::cpdb_data$pathway_matrix[glyc[i], ] != 0)), names(which(oranges::cpdb_data$pathway_matrix[glyc[j], ] != 0))))
    aintb <- length(intersect(names(which(oranges::cpdb_data$pathway_matrix[glyc[i], ] != 0)), names(which(oranges::cpdb_data$pathway_matrix[glyc[j], ] != 0))))
    O[i, j] <- aintb / aub
    
  }
  
  glyco_comparison[[i]] <- tibble::tibble(pathway_name = oranges::cpdb_data$pathway_info$pathway_name[glyc[i]],
                                          pathway_source = psource,
                                          number_genes = gcount)
  
}
glyco_comparison <- do.call(what = rbind, glyco_comparison)

# plot differences
glyco_comparison %>% 
  ggplot(aes(x = pathway_source, y = number_genes)) + 
  geom_bar(stat = "identity") + 
  labs(x = NULL, y = "Number of genes in gene set") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 14),
        panel.background = element_blank())

# plot matrix of similarities
rownames(O) <- glyco_comparison$pathway_source
colnames(O) <- glyco_comparison$pathway_source


corrplot::corrplot(O, 
                   method = "circle", 
                   type = "upper", 
                   tl.col = "black", 
                   is.corr = FALSE, 
                   tl.cex = 1)

