# load("~/work/data/lung_cancer_gse19804.RData")

# source("/Volumes/HOME/scripts/r/expression_analysis/limma_dge.R")

# simplest case is comparing (pooled) normal samples to (pooled) disease samples
tumor_lung <- grep("tumor", samples_lung$clinical_diagnosis)
tumor_lung <- which(colnames(expression_lung) %in% samples_lung$cel_file[tumor_lung])
healthy_lung <- grep("normal",samples_lung$clinical_diagnosis)
healthy_lung <- which(colnames(expression_lung) %in% samples_lung$cel_file[healthy_lung])

# limma for differential expression
lung_res <- limma_dge(expression_data = expression_lung,caseIds = tumor_lung,ctrIds = healthy_lung)

##### RUN DRUID ####
nsclc_druid <- DRUID::concoct(dge_matrix = cbind(lung_res$logFC,lung_res$adj.P.Val), 
                       num_random = 10000, 
                       druid_direction = "neg", 
                       fold_thr = 1, 
                       pvalue_thr = 0.01,
                       entrez = genes_lung$ENTREZID)

#### DO THE SAME ANALYSES FOR DISEASES
a1 <- DRUID::druid_geneset(dge_matrix = cbind(lung_res$logFC,lung_res$adj.P.Val), 
                           desired_effect = "pos", 
                           fold_thr = 1, 
                           pvalue_thr = 0.01, 
                           entrez = geo_genes$`GeneID/GSE`, 
                           gene_space = colnames(geo_tfidf))

a2 <- DRUID::cosine_similarity(query_vector = a1, 
                               tfidf_matrix = geo_tfidf, 
                               tfidf_crossprod_mat = geo_cpm)

a3 <- DRUID::random_probability(similarity_results = a2, 
                                gs_size = sum(a1), 
                                num_sets = 10000, 
                                target_tfidf = geo_tfidf, 
                                tfidf_crossprod_mat = geo_cpm)

a4 <- DRUID::druid_score(similarity_results = a2, random_probabilities = a3, num_random = 10000)

# results ----
# message("Building results data frame...")
disease_druid <- tibble::tibble(cosine_similarity = as.vector(a2),
              probability_random = a3,
              druid_score = as.vector(a4))

disease_druid <- disease_druid %>% 
  tibble::add_column(., query_size = sum(a1), .before = 1)

disease_druid <- disease_druid %>% 
  tibble::add_column(., disease = Y$Disease, .before = 1) %>%
  tibble::add_column(., sample_source = Y$sample_source, .before = 2) %>%
  dplyr::arrange(., desc(druid_score))


#### DO THE SAME ANALYSES FOR CCLE
a1 <- DRUID::druid_geneset(dge_matrix = cbind(lung_res$logFC,lung_res$adj.P.Val), 
                           desired_effect = "pos", 
                           fold_thr = 1, 
                           pvalue_thr = 0.01, 
                           entrez = ccle_genes$`GeneID/NA`, 
                           gene_space = colnames(ccle_tfidf))

a2 <- DRUID::cosine_similarity(query_vector = a1, 
                               tfidf_matrix = ccle_tfidf, 
                               tfidf_crossprod_mat = ccle_cpm)

a3 <- DRUID::random_probability(similarity_results = a2, 
                                gs_size = sum(a1), 
                                num_sets = 10000, 
                                target_tfidf = ccle_tfidf, 
                                tfidf_crossprod_mat = ccle_cpm)

a4 <- DRUID::druid_score(similarity_results = a2, 
                         random_probabilities = a3, 
                         num_random = 10000)

# results ----
# message("Building results data frame...")
ccle_druid <- tibble::tibble(cosine_similarity = as.vector(a2),
                                probability_random = a3,
                                druid_score = as.vector(a4))

ccle_druid <- ccle_druid %>% 
  tibble::add_column(., query_size = sum(a1), .before = 1)

ccle_druid <- ccle_druid %>% 
  tibble::add_column(., cell_line = Y$CellLine, .before = 1) %>%
  tibble::add_column(., cell_line_tissue = Y$Tissue, .before = 2) %>%
  dplyr::arrange(., desc(druid_score))



##### DRUG TO DISEASE ####
drug_id <- which(cauldron::cmap_instances$name == "monobenzone")[1]
test_profile <- cauldron::cmap_dataset$M[, drug_id]

a1 <- DRUID::druid_geneset(dge_matrix = cbind(test_profile, 0), 
                           desired_effect = "neg", 
                           fold_thr = 0, 
                           pvalue_thr = 1, 
                           entrez = cauldron::cmap_dataset$G$var2, 
                           gene_space = colnames(tcga_tfidf))

a2 <- DRUID::cosine_similarity(query_vector = a1, 
                               tfidf_matrix = tcga_tfidf, 
                               tfidf_crossprod_mat = tcga_cpm)

a3 <- DRUID::random_probability(similarity_results = a2, 
                                gs_size = sum(a1), 
                                num_sets = 10000, 
                                target_tfidf = tcga_tfidf, 
                                tfidf_crossprod_mat = tcga_cpm)

a4 <- DRUID::druid_score(similarity_results = a2, 
                         random_probabilities = a3, 
                         num_random = 10000)

tcga_drug <- tibble::tibble(cosine_similarity = as.vector(a2),
                            probability_random = a3,
                            druid_score = as.vector(a4))

tcga_drug <- tcga_drug %>% 
  tibble::add_column(., query_size = sum(a1), .before = 1)

tcga_drug <- tcga_drug %>% 
  tibble::add_column(., cancer_type = tcga_metadata$`Cancer Name_Cancer Acronym`, .before = 1) %>%
  # tibble::add_column(., sample_tissue = Y$sample_source, .before = 2) %>%
  dplyr::arrange(., desc(druid_score))

# compare to random
sample_tcga <- sample(x = 1:nrow(tcga_metadata), size = 10000, replace = TRUE)
tcga_metadata[sample_tcga, ] %>% 
  dplyr::count(`Cancer Name_Cancer Acronym`) %>% 
  dplyr::mutate(., prop_random = n / 10000) %>%
  dplyr::left_join(., y = tcga_metadata %>% 
                     dplyr::count(`Cancer Name_Cancer Acronym`) %>% 
                     dplyr::mutate(., prop_orig = n / nrow(tcga_metadata)), by = "Cancer Name_Cancer Acronym") %>%
  dplyr::left_join(., y = tcga_drug %>% 
                     dplyr::slice(1:100) %>% 
                     dplyr::count(., cancer_type) %>% 
                     dplyr::mutate(., prop_druid = n / 100), by = c("Cancer Name_Cancer Acronym" = "cancer_type")) %>%  
  dplyr::select(`Cancer Name_Cancer Acronym`, prop_orig, prop_random, prop_druid) %>%
  dplyr::arrange(., desc(prop_druid))

dplyr::select(`Cancer Name_Cancer Acronym`, prop_orig, prop_random) %>%
  dplyr::arrange(., desc(prop_random))



# Towards an MoA ----------------------------------------------------------
# using the GTEx data to compute tissue specificity
utissues <- unique(sample_mapping$SMTSD)

sel_tissues <- vector(length = 0L)
number_samples <- vector(length = length(utissues))
for(i in seq(1,length(utissues))) {
  
  usamps <- which(sample_mapping$SMTSD == utissues[i])
  nsamps <- length(usamps)
  number_samples[i] <- nsamps
  if(nsamps > 10) {
    rsamps <- sample(usamps, 10, replace = FALSE)
  } else {
    rsamps <- usamps
  }
  
  sel_tissues <- c(sel_tissues, rsamps)
  
}

cids <- which(colnames(M) %in% sample_mapping$SAMPID[sel_tissues])
sampled_expression <- M[, cids]

# keep only protein coding genes
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'www.ensembl.org')
genes <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name","transcript_biotype"), filters = c("transcript_biotype","chromosome_name"),values = list("protein_coding",c(1:22)), mart = mart)
gids <- which(G %in% genes$external_gene_name)
sampled_expression <- sampled_expression[gids, ]
G <- G[gids]

# clean up some stuff
sids <- sample_mapping$SAMPID
sname <- character(length = length(sids))
for (i in seq(1, length(sids))) {
  
  sname[i] <- paste(strsplit(sids[i], split = "-")[[1]][c(1,2)], collapse = "-")
  
}


samples_df <- tibble::tibble(sample_id = colnames(sampled_expression),
                             tissue_source = sample_mapping$SMTS[which(sample_mapping$SAMPID %in% colnames(sampled_expression))],
                             tissue_detail = sample_mapping$SMTSD[which(sample_mapping$SAMPID %in% colnames(sampled_expression))]) #,


### calculate sample specificity
# scale the data for any sample
column_scale <- scale(x = sampled_expression)

# scale the data for any gene
row_scale <- t(scale(t(x = sampled_expression)))


specificity_z <- (row_scale + column_scale) / sqrt(2)

# now look at our genes in the context of the lung tissue
gin <- strsplit(nsclc_druid$matched_genes[1], " \\| ")[[1]]
gin_up <- sapply(gin[grep("up", gin)], function(y) strsplit(y, " up")[[1]][1])
gin_down <- sapply(gin[grep("down", gin)], function(y) strsplit(y, " down")[[1]][1])
mono_genes <- c(gin_down, gin_up)

which(G %in% mono_genes)
