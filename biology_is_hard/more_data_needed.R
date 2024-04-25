# Show me how biology is hard
# load some libraries first
library(cauldron)
library(oranges)
library(tidyverse)
library(ggplot2)

# get valproic acid data
valproic50 <- which(cauldron::cmap_instances$name == "valproic acid" & 
      cauldron::cmap_instances$cell_line == "MCF7" & 
      cauldron::cmap_instances$concentration == 50 & 
      cauldron::cmap_instances$vehicle == "DMSO")


A <- cauldron::cmap_dataset$M[, valproic50]

# calculate number perturbation impact on genes and samples
number_pgenes <- apply(A, 2, function(x) length(which(x != 0)))
genes_pert <- apply(A, 1, function(x) length(which(x != 0)))

# get set of genes perturbed by valproic acid
va_gene_set <- vector(mode = "character")
for (i in seq(1, ncol(A))) {
  id1 <- which(A[, i] != 0)
  gsym <- cauldron::cmap_dataset$G[id1, 1]
  va_gene_set <- union(va_gene_set, gsym)
}

# do the same but for all data
all_perts <- apply(cauldron::cmap_dataset$M, 2, function(x) length(which(x != 0)))
genes_all <- apply(cauldron::cmap_dataset$M, 1, function(x) length(which(x != 0)))
hist(all_perts)

# plot some things
hist(genes_pert[genes_pert != 0])

### hypergeometric tests
va_enrichment_all <- vector(mode = "list", length = ncol(A))
for (i in seq(1, ncol(A))) {
  id1 <- which(A[, i] != 0)
  gsym <- cauldron::cmap_dataset$G[id1, 2]
  enr <- oranges::oranges(query_entrez = gsym, universe_entrez = cauldron::cmap_dataset$G[, 2])
  enr <- enr %>% dplyr::filter(p_val < 0.05) %>% tibble::add_column(sample_id = i)
  va_enrichment_all[[i]] <- enr
}
va_enrichment_all <- do.call(rbind, va_enrichment_all)

va_enrichment_all %>% 
  dplyr::filter(p_val < 0.001) %>%
  ggplot(aes(x = -log10(p_val), y = name, fill = factor(sample_id))) + 
  geom_point(shape = 21) + 
  labs(x = "-log10(enrichment p-value)", y = "Pathway") + 
  theme_bw() + 
  theme(axis.text.y = element_blank(),
        panel.grid = element_blank())

# are there commonalities?
va_enrichment_all %>% 
  dplyr::count(name) %>% 
  dplyr::filter(n > 2) %>%
  ggplot(aes(x = name, y = n)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5))

# what is covered?
va_enrichment_all %>% 
  dplyr::filter(name == "notch_signaling_pathway") %>%
  ggplot(aes(x = data_source, y = pathway_proportion)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5))