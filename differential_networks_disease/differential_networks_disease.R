#' ---
#' title: Differential networks for disease characterization
#' author: Diogo M. Camacho
#' date: 4/25/2024
#' ---

# load data
load(file = "differential_networks_disease/lung_cancer_gse19804.RData")
load(file = "biogrid_human.RData")

# source functions
source("./limma_dge.R")

##### biogrid network #####
g <- A %>%
  dplyr::select(entrez_a, entrez_b) %>%
  igraph::graph_from_data_frame(directed = FALSE)

M <- igraph::as_adjacency_matrix(g)
M[M > 1] <- 1
diag(M) <- 0
GRN <- igraph::graph_from_adjacency_matrix(M, mode = "undirected") %>% 
  igraph::as_data_frame() %>% 
  tibble::as_tibble()

##### subset genes in data that are in network #####
GRN_lung <- GRN %>% dplyr::filter(from %in% genes_lung$ENTREZID, to %in% genes_lung$ENTREZID)

##### calculate differential network: correlation based #####
set.seed(08030727)
exids <- sample(x = union(GRN_lung$from, GRN_lung$to), size = 2500, replace = FALSE)

example_data <- expression_lung[which(genes_lung$ENTREZID %in% exids), ]
example_genes <- genes_lung[which(genes_lung$ENTREZID %in% exids), ]

D <- diffnet::diff_top(data = example_data, 
                       group1 = samples_lung$clinical_diagnosis == "normal", 
                       group2 = samples_lung$clinical_diagnosis == "primary tumor")

D <- D %>% 
  tibble::add_column(entrez_a = example_genes$ENTREZID[D$x],
                     entrez_b = example_genes$ENTREZID[D$y],
                     symbol_a = example_genes$SYMBOL[D$x],
                     symbol_b = example_genes$SYMBOL[D$y], .before = 1) %>%
  dplyr::select(-x, -y)

# add interaction column for comparison
D <- D %>% tibble::add_column(., interaction = paste(D$entrez_a, D$entrez_b, sep = " -- "))
G <- GRN_lung %>% tibble::add_column(., interaction = paste(GRN_lung$from, GRN_lung$to, sep = " -- "))


##### igraph object #####
g2 <- D %>%
  dplyr::select(entrez_a, entrez_b) %>%
  igraph::graph_from_data_frame(directed = FALSE)

graph_cent <- igraph::centr_degree(g2)


##### simple node metric #####
ag <- union(D$symbol_a, D$symbol_b)
int_count <- vector(length = length(ag))
for(i in seq(1, length(ag))) {
  
  int_count[i] <- length(which(D$symbol_a == ag[i])) + length(which(D$symbol_b == ag[i]))
  
}

# make a better data frame
topdiff_ratio <- vector(length = length(ag), mode = "list")
for(i in seq(1, length(ag))) {
  
  x <- D %>% dplyr::filter(symbol_a == ag[i] | symbol_b == ag[i]) %>% dplyr::count(change_type)
  
  y1 <- x$n[x$change_type == "gain of edge"]
  y2 <- x$n[x$change_type == "loss of edge"]
  y3 <- x$n[x$change_type == "edge present, same sign"]
  
  if(length(y1) != 0 & length(y2) != 0) {
    topdiff_ratio[[i]] <- tibble::tibble(gene = ag[i],
                                         change_type = c("gain", "loss"),
                                         count = c(y1, y2))
    } else {
      if(length(y1) == 0 & length(y2) != 0) {
        topdiff_ratio[[i]] <- tibble::tibble(gene = ag[i],
                                             change_type = c("gain", "loss"),
                                             count = c(0, y2))
      } else {
        if(length(y1) != 0 & length(y2) == 0) {
          topdiff_ratio[[i]] <- tibble::tibble(gene = ag[i],
                                               change_type = c("gain", "loss"),
                                               count = c(y1, 0))
        }
      }
    }
  
}
topdiff_ratio <- do.call(rbind, topdiff_ratio)

# plot 
topdiff_ratio %>% 
  dplyr::filter(count > 40) %>% 
  ggplot(aes(y = forcats::fct_reorder(gene, count, .desc = TRUE), x = count, fill = change_type)) + 
  geom_point(shape = 21) + 
  labs(x = "number of edges", y = NULL) +
  theme_bw()


# plot some things
D %>% 
  dplyr::count(symbol_a, change_type) %>% 
  dplyr::filter(change_type == "gain of edge" | change_type == "loss of edge") %>% 
  dplyr::filter(n > 25) %>% 
  ggplot(aes(y = forcats::fct_reorder(symbol_a, n, .desc = TRUE), x = n, fill = change_type)) + 
  geom_point(shape = 21) + 
  labs(x = "number of edges", y = NULL) +
  theme_bw()

D %>% 
  dplyr::count(symbol_b, change_type) %>% 
  dplyr::filter(change_type == "gain of edge" | change_type == "loss of edge") %>% 
  dplyr::filter(n > 25) %>% 
  ggplot(aes(y = forcats::fct_reorder(symbol_b, n, .desc = TRUE), x = n, fill = change_type)) + 
  geom_point(shape = 21) + 
  theme_bw()

D %>% 
  dplyr::count(symbol_a, change_type) %>% 
  dplyr::filter(change_type == "loss of edge") %>% 
  dplyr::filter(n > 30) %>% 
  ggplot(aes(y = forcats::fct_reorder(symbol_b, n, .desc = TRUE), x = n, fill = change_type)) + 
  geom_point(shape = 21) + 
  labs(x = "number of edges", y = NULL) +
  theme_bw()

##### DIFFERENCES TO DGE #####
# simplest case is comparing (pooled) normal samples to (pooled) disease samples
tumor_lung <- grep("tumor",samples_lung$clinical_diagnosis)
tumor_lung <- which(colnames(expression_lung) %in% samples_lung$cel_file[tumor_lung])
healthy_lung <- grep("normal",samples_lung$clinical_diagnosis)
healthy_lung <- which(colnames(expression_lung) %in% samples_lung$cel_file[healthy_lung])

# limma for differential expression
lung_res <- limma_dge(expression_data = expression_lung,caseIds = tumor_lung,ctrIds = healthy_lung)

lung_res <- as_tibble(lung_res) %>%
  tibble::add_column(gene_symbol = genes_lung$SYMBOL, .before = 1)

top_diff_genes <- lung_res %>% dplyr::filter(abs(logFC) > 1, adj.P.Val < 0.05)


##### Complementing network analyses and DGE #####
genes_example <- tibble(ag, int_count) %>% dplyr::arrange(desc(int_count)) %>% dplyr::filter(., int_count >= 100)

diff_genes_net <- vector(mode = "list", length = dim(genes_example)[1])
for(i in seq(1, length(diff_genes_net))) {
  x1 <- D %>% dplyr::filter(symbol_a == genes_example$ag[i] | symbol_b == genes_example$ag[i]) %>% dplyr::select(symbol_a, symbol_b)
  x2 <- setdiff(union(x1$symbol_a, x1$symbol_b), genes_example$ag[i])
  x3 <- intersect(x2, top_diff_genes$gene_symbol)
  
  diff_genes_net[[i]] <- tibble::tibble(gene_interest = genes_example$ag[i],
                                        number_interactions = length(x2),
                                        number_differential_interactions = length(x3))
  
}
diff_genes_net <- do.call(rbind, diff_genes_net)


# ##### subsetting to experimental interactions #####
# D_known <- D[which(D$interaction %in% G$interaction), ]
# 
# ag_known <- union(D_known$symbol_a, D_known$symbol_b)
# int_count_known <- vector(length = length(ag_known))
# for(i in seq(1, length(ag_known))) {
#   
#   int_count_known[i] <- length(which(D_known$symbol_a == ag_known[i])) + length(which(D_known$symbol_b == ag_known[i]))
#   
# }
# 
# 
# topdiff_known <- vector(length = length(ag_known), mode = "list")
# for(i in seq(1, length(ag_known))) {
#   
#   x <- D_known %>% dplyr::filter(symbol_a == ag[i] | symbol_b == ag[i]) %>% dplyr::count(change_type)
#   
#   y1 <- x$n[x$change_type == "gain of edge"]
#   y2 <- x$n[x$change_type == "loss of edge"]
#   y3 <- x$n[x$change_type == "edge present, same sign"]
#   
#   if(length(y1) != 0 & length(y2) != 0) {
#     topdiff_known[[i]] <- tibble::tibble(gene = ag_known[i],
#                                          change_type = c("gain", "loss"),
#                                          count = c(y1, y2))
#   } else {
#     if(length(y1) == 0 & length(y2) != 0) {
#       topdiff_known[[i]] <- tibble::tibble(gene = ag_known[i],
#                                            change_type = c("gain", "loss"),
#                                            count = c(0, y2))
#     } else {
#       if(length(y1) != 0 & length(y2) == 0) {
#         topdiff_known[[i]] <- tibble::tibble(gene = ag_known[i],
#                                              change_type = c("gain", "loss"),
#                                              count = c(y1, 0))
#       }
#     }
#   }
#   
# }
# topdiff_known <- do.call(rbind, topdiff_known)
# 
# # plot 
# topdiff_known %>% 
#   ggplot(aes(y = forcats::fct_reorder(gene, count, .desc = TRUE), x = count, fill = change_type)) + 
#   geom_point(shape = 21) + 
#   labs(x = "number of edges", y = NULL) +
#   theme_bw()
# 
# 
# 
# 
# D_known %>% 
#   dplyr::count(symbol_a, change_type) %>% 
#   dplyr::filter(change_type == "gain of edge" | change_type == "loss of edge") %>% 
#   # dplyr::filter(n > 25) %>% 
#   ggplot(aes(y = forcats::fct_reorder(symbol_a, n, .desc = TRUE), x = n, fill = change_type)) + 
#   geom_point(shape = 21) + 
#   theme_bw()
# 
# D_known %>% 
#   dplyr::count(symbol_b, change_type) %>% 
#   dplyr::filter(change_type == "gain of edge" | change_type == "loss of edge") %>% 
#   # dplyr::filter(n > 25) %>% 
#   ggplot(aes(y = forcats::fct_reorder(symbol_b, n, .desc = TRUE), x = n, fill = change_type)) + 
#   geom_point(shape = 21) + 
#   theme_bw()