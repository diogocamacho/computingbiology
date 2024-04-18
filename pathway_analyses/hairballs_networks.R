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


# write out the network
glycolytic_network <- B %>% 
  dplyr::filter(entrez_a %in% genes_glycolysis | entrez_b %in% genes_glycolysis) %>%
  dplyr::select(biogrid_id, entrez_a, entrez_b, gene_a, gene_b)  %>%
  tibble::add_column(in_glycolysis_a = 0, in_glycolysis_b = 0)
glycolytic_network$in_glycolysis_a[glycolytic_network$entrez_a %in% genes_glycolysis] <- 1
glycolytic_network$in_glycolysis_b[glycolytic_network$entrez_b %in% genes_glycolysis] <- 1

interaction_pair <- vector(length = nrow(glycolytic_network))
for(i in seq(1, nrow(glycolytic_network))) {
  
  interaction_pair[i] <- paste(sort(c(glycolytic_network$entrez_a[i], glycolytic_network$entrez_b[i])), collapse = " :: ")
  
}

glycolytic_network$interaction_pair <- interaction_pair

glycolytic_network <- glycolytic_network %>% 
  dplyr::distinct(interaction_pair, .keep_all = TRUE) %>%
  dplyr::filter(gene_a != gene_b)


gene_annotations <- tibble::tibble(entrez_id = c(glycolytic_network$entrez_a, glycolytic_network$entrez_b),
                                   gene_name = c(glycolytic_network$gene_a, glycolytic_network$gene_b),
                                   in_glycolysis = c(glycolytic_network$in_glycolysis_a, glycolytic_network$in_glycolysis_b)) %>% 
  dplyr::distinct()

write.csv(x = glycolytic_network, file = "pathway_analyses/glycolysis_network.csv", row.names = FALSE, quote = FALSE)
write.csv(x = gene_annotations, file = "pathway_analyses/gene_annotations.csv", row.names = FALSE, quote = FALSE)


##### network analytics #####
g <- glycolytic_network %>%
  dplyr::select(gene_a, gene_b) %>%
  igraph::graph_from_data_frame(directed = FALSE)

# 
gene_annotations %>% 
  dplyr::mutate(gene_centrality = degree(g)) %>%
  dplyr::mutate(gene_betweeness = betweenness(g)) %>%
  dplyr::arrange(desc(gene_betweeness))


# a plot of the data frame above
gene_annotations %>% 
  dplyr::mutate(gene_centrality = degree(g)) %>%
  dplyr::mutate(gene_betweeness = betweenness(g)) %>%
  dplyr::mutate(normalized_betweeness = gene_betweeness / max(gene_betweeness)) %>%
  dplyr::arrange(desc(gene_betweeness)) %>% 
  dplyr::slice(1:50) %>%
  ggplot(aes(x = forcats::fct_reorder(gene_name, normalized_betweeness, .desc = TRUE), y = normalized_betweeness)) + 
  geom_point(shape = 21) + 
  labs(x = "Gene name", y = "Betweeness") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


##### community detection #####
leiden_communities <- cluster_leiden(graph = g, 
                                     objective_function = "CPM", 
                                     resolution_parameter = 0.01, 
                                     n_iterations = 1000)


