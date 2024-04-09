# Hairballs and biological organization
# load("~/work/data/biogrid_human.RData")
library(oranges)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(biomaRt)


# keep only protein coding genes
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'www.ensembl.org')
genes <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name","transcript_biotype"), filters = c("transcript_biotype","chromosome_name"),values = list("protein_coding",c(1:22)), mart = mart)

z <- AnnotationDbi::select(x = org.Hs.eg.db, 
                           keys = unique(c(A$entrez_a, A$entrez_b)), 
                           keytype = "ENTREZID", 
                           columns = c("SYMBOL", "GENENAME"))

biogrid_pcg <- intersect(which(A$gene_a %in% genes$external_gene_name), which(A$gene_b %in% genes$external_gene_name))

B <- A[biogrid_pcg, ]

reactome_glycolysis <- which(cpdb_data$pathway_info$pathway_source == "reactome" & cpdb_data$pathway_info$pathway_name == "glycolysis")
genes_glycolysis <- names(which(cpdb_data$pathway_matrix[reactome_glycolysis, ] == 1))


glycolytic_network <- B %>% 
  dplyr::filter(entrez_a %in% genes_glycolysis | entrez_b %in% genes_glycolysis) %>%
  dplyr::select(biogrid_id, entrez_a, entrez_b, gene_a, gene_b)

gene_annotations <- tibble::tibble(entrez_id = c(glycolytic_network$entrez_a, glycolytic_network$entrez_b),
                                   gene_name = c(glycolytic_network$gene_a, glycolytic_network$gene_b)) %>% 
  dplyr::distinct() %>%
  tibble::add_column(in_glycolysis = 0) %>%
  dplyr::mutate(in_glycolysis = replace(entrez_id %in% genes_glycolysis, 1, 0))
