#' ---
#' title: Network modules and drug discovery
#' author: Diogo M. Camacho
#' date: 4/25/2024
#' ---

# load libraries
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(biomaRt)
library(igraph)
library(patchwork)

# load biogrid data
load(file = "biogrid_human.RData")

##### write out network #####
g <- A %>%
  dplyr::select(entrez_a, entrez_b) %>%
  igraph::graph_from_data_frame(directed = FALSE)

##### identify network communities #####
leiden_communities <- igraph::cluster_leiden(graph = g, 
                                     objective_function = "CPM", 
                                     resolution_parameter = .1, 
                                     n_iterations = 1000)

# size of communities
cluster_size <- integer(length = leiden_communities$nb_clusters)
for(i in seq(1, length(cluster_size))) {
  cluster_size[i] <- length(V(g)[which(leiden_communities$membership == i)])
}

leiden_clusters <- tibble::tibble(cluster_name = paste("cluster_", seq(1, leiden_communities$nb_clusters), sep = ""),
                               cluster_size = cluster_size)

leiden_clusters %>% 
  ggplot(aes(x = cluster_name, y = cluster_size)) + 
  geom_point(shape = 21, color = "black", fill = "red", size = 2, alpha = 0.5) +
  labs(x = "Cluster", y = "Cluster size") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# plot(cluster_size)
# leiden_communities$nb_clusters
# leiden_communities$nb_clusters - length(which(cluster_size < 5))


##### summarizing perturbations on communities #####
# before ignoring clusters, let's put them in the context of the drug data
# cauldron::cmap_dataset$G$var2

drug_clusters <- vector(mode = "list", length = leiden_communities$nb_clusters)
for (i in seq(1, length(drug_clusters))) {
  
  gene_members <- intersect(cauldron::cmap_dataset$G$var2, names(V(g)[which(leiden_communities$membership == i)]))
  cluster_size <- length(gene_members)

  drug_clusters[[i]] <- tibble::tibble(cluster_id = i,
                                       cluster_size = cluster_size,
                                       cluster_name = paste("cluster_", i, sep = ""),
                                       gene_members = gene_members)
  
}
drug_clusters <- do.call(rbind, drug_clusters)


# compute a cluster "activity" for each drug
cn <- unique(drug_clusters$cluster_id)
dc_data <- matrix(0, nrow = leiden_communities$nb_clusters, ncol = ncol(cauldron::cmap_dataset$M))

for(i in seq(1, length(cn))) {
  
  tmp1 <- drug_clusters$gene_members[drug_clusters$cluster_id == cn[i]]
  tmp2 <- which(cauldron::cmap_dataset$G$var2 %in% tmp1)
  tmp3 <- cauldron::cmap_dataset$M[tmp2, ]
  if(length(tmp2) > 1) {
    dc_data[i, ] <- Matrix::colSums(tmp3)  
  } else {
    dc_data[i, ] <- tmp3
  }
  
  
}

# example plots
# cid <- sample(x = nrow(dc_data), 1)
cid <- 321
did <- sample(x = ncol(dc_data), 1)

drugvis_example <- tibble::tibble(cluster_effect = dc_data[, did],
                                  cluster_name = paste("cluster_", seq(1, leiden_communities$nb_clusters), sep = ""))

clustervis_example <- tibble::tibble(drug_name = cauldron::cmap_instances$name,
                                     drug_effect = dc_data[cid, ])

p1 <- drugvis_example %>% 
  # dplyr::count(cluster_effect) %>% 
  # ggplot(aes(x = cluster_effect, y = n)) + 
  ggplot(aes(x = cluster_name, y = cluster_effect)) +
  geom_point(shape = 21, color = "black", fill = "red", size = 2, alpha = 0.5) +
  # geom_bar(stat = "identity", color = "black", fill = "red", alpha = 0.5) + 
  labs(x = "Cluster", y = "Drug effect on cluster") + 
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank())

p2 <- clustervis_example %>% 
  # dplyr::count(cluster_effect) %>% 
  # ggplot(aes(x = cluster_effect, y = n)) + 
  ggplot(aes(x = drug_name, y = drug_effect)) +
  geom_point(shape = 21, color = "black", fill = "red", size = 2, alpha = 0.5) +
  # geom_bar(stat = "identity", color = "black", fill = "red", alpha = 0.5) + 
  labs(x = "Drug", y = "Drug effect on cluster") + 
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank())


p1 + p2

# how many drugs does a given cluster respond to? 
num_effects <- tibble::tibble(cluster_name = paste("cluster_", seq(1, leiden_communities$nb_clusters), sep = ""),
                              perturbing_drugs = apply(dc_data, 1, function(y) length(which(y != 0))))

num_effects %>% 
  # ggplot(aes(x = forcats::fct_reorder(cluster_name, perturbing_drugs, .desc = TRUE), y = perturbing_drugs)) + 
  ggplot(aes(x = cluster_name, y = perturbing_drugs)) + 
  geom_point(shape = 21, fill = "red", color = "black", size = 2, alpha = 0.5) +
  labs(x = "Cluster", y = "Number of drug that modulate cluster") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

##### Clean up the community data set #####
nodrugeffect <- apply(dc_data, 1, function(y) length(which(y == 0)))
nix_clusters <- which(nodrugeffect == ncol(dc_data))
dc_data <- dc_data[-nix_clusters, ]

##### DRUID style TF-IDF #####
# split data into up and down effects for clusters
# the TFIDF matrix can only be computed on positive integers, so that's why we do this trick
cnames_up <- paste(unique(drug_clusters$cluster_name), "up")
cnames_down <- paste(unique(drug_clusters$cluster_name), "down")

a1 <- matrix(data = 0, nrow = ncol(dc_data), ncol = length(cnames_up))
a2 <- matrix(data = 0, nrow = ncol(dc_data), ncol = length(cnames_down))
for (i in seq(1, nrow(a1))) {
  x1 <- dc_data[, i]
  upcl <- which(x1 > 0)
  dncl <- which(x1 < 0)
  a1[i, upcl] <- x1[upcl]
  a2[i, dncl] <- abs(x1[dncl])
}
a3 <- cbind(a1, a2)
rownames(a3) <- cauldron::cmap_instances$name
colnames(a3) <- c(cnames_up, cnames_down)

# now calculate the TFIDF with this matrix
# tidy data
x1 <-  tm::as.DocumentTermMatrix(x = a3, weighting = 1) # <-- drugs as documents, clusters as terms
x2 <-  tm::as.DocumentTermMatrix(x = t(a3), weighting = 1) # <-- clusters as documents, drugs as terms
dtm1 <- tidytext::tidy(x1)
dtm2 <- tidytext::tidy(x2)

# count words in both dtms
words1 <- dtm1 %>% 
  dplyr::count(term, document) %>% 
  dplyr::ungroup()

words2 <- dtm2 %>% 
  dplyr::count(term, document) %>% 
  dplyr::ungroup()

# count all terms
terms1 <- words1 %>% 
  dplyr::group_by(term) %>% 
  dplyr::summarize(total = sum(n))

terms2 <- words2 %>% 
  dplyr::group_by(term) %>% 
  dplyr::summarize(total = sum(n))

# calculate tf-idf
tfidf1 <- words1 %>% 
  tidytext::bind_tf_idf(document, term, n)

tfidf2 <- words2 %>% 
  tidytext::bind_tf_idf(document, term, n)

# combine into 1 tfidf matrix
# output matrix is tfidf1 * t(tfidf2)
y1 <- tfidf1 %>% dplyr::arrange(., term)
y2 <- tfidf2 %>% dplyr::arrange(., document)
y3 <- y1$tf_idf * y2$tf_idf

M <- matrix(0, nrow = nrow(a3), ncol = ncol(a3))
for (i in seq(1, nrow(a3))) {
  b1 <- which(y1$document == rownames(a3)[i]) # <-- row ids for drugs
  b2 <- match(y1$term[b1], colnames(a3)) # <-- col ids for clusters
  M[i, b2] <- y3[b1]
}
rownames(M) <- rownames(a3)
colnames(M) <- colnames(a3)

# make sparse matrix
M <- Matrix::Matrix(M, sparse = TRUE)

# compute cross product for TFIDF
cpm <- apply(M, 1, crossprod)



##### Disease data #####
# the disease data will come from an analyses of the differentially expressed data
# because our collapsing of the perturbation data is based on a sum of perturbation across 
# the different genes in a cluster, we will do something similar after we compute
# the differential expression

# load data
load(file = "differential_networks_disease/lung_cancer_gse19804.RData")

# source functions
source("./limma_dge.R")

# simplest case is comparing (pooled) normal samples to (pooled) disease samples
tumor_lung <- grep("tumor",samples_lung$clinical_diagnosis)
tumor_lung <- which(colnames(expression_lung) %in% samples_lung$cel_file[tumor_lung])
healthy_lung <- grep("normal",samples_lung$clinical_diagnosis)
healthy_lung <- which(colnames(expression_lung) %in% samples_lung$cel_file[healthy_lung])

# limma for differential expression
lung_res <- limma_dge(expression_data = expression_lung,caseIds = tumor_lung,ctrIds = healthy_lung)

# binarize differential expression as +1 for up and -1 for down
lung_dge <- lung_res %>% 
  as_tibble %>% 
  tibble::add_column(gene_symbol = genes_lung$SYMBOL, 
                     entrez = genes_lung$ENTREZID, 
                     effect = 0) 
lung_dge$effect[lung_res$logFC > 1 & lung_res$adj.P.Val < 0.05] <- 1
lung_dge$effect[lung_res$logFC < -1 & lung_res$adj.P.Val < 0.05] <- -1

# compute clusters for lung differential expression
cn <- unique(drug_clusters$cluster_id)
lung_clusters <- integer(length = leiden_communities$nb_clusters)

for(i in seq(1, length(cn))) {
  
  tmp1 <- drug_clusters$gene_members[drug_clusters$cluster_id == cn[i]]
  tmp2 <- which(lung_dge$entrez %in% tmp1)
  lung_clusters[i] <- sum(lung_dge$effect[tmp2])

}

lung_clusters <- cbind(lung_clusters, 0.01)

lung_cnames <- paste("cluster_", cn, sep = "")

# now compute the cosine similarity between our disease data and the drug TFIDF
# we will use DRUID's framework
qv <- druid_geneset(dge_matrix = lung_clusters, desired_effect = "neg", entrez = lung_cnames, gene_space = colnames(M))
qv_cs <- cosine_similarity(query_vector = qv, tfidf_matrix = M, tfidf_crossprod_mat = cpm)
prandom <- random_probability(similarity_results = qv_cs, 
                              gs_size = length(which(lung_clusters[, 1] != 0)), 
                              num_sets = 1000, 
                              target_tfidf = M, 
                              tfidf_crossprod_mat = cpm)
prandom[which(t2 < 3)] <- 1

dscore <- DRUID::druid_score(similarity_results = qv_cs, 
                      random_probabilities = prandom, 
                      num_random = 1000)

res <- tibble(cosine_similarity = as.vector(qv_cs),
              probability_random = prandom,
              druid_score = as.vector(dscore))

res <- res %>% 
  tibble::add_column(., query_size = sum(qv), .before = 1) %>%
  tibble::add_column(., number_matches = t2, .before = 2) # %>%
  # tibble::add_column(., percent_matched = t2 / sum(query_vector), .before = 3) %>%
  # tibble::add_column(., matched_genes = b1, .before = 4)

res <- res %>% 
  tibble::add_column(., drug_name = rownames(M), .before = 1) %>%
  # tibble::add_column(., concentration = cauldron::cmap_instances$concentration, .before = 2) %>%
  # tibble::add_column(., cell_line = as.character(cauldron::cmap_instances$cell_line), .before = 3) %>%
  # tibble::add_column(., data_source = names(cauldron::druid_potion)[selection], .before = 1) %>%
  dplyr::arrange(., desc(druid_score)) %>%
  dplyr::filter(., number_matches >= 5) %>% 
  dplyr::distinct()

