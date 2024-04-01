# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Data plots for selected GEO samples

library(DESeq2)
library(GEOquery)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(biomaRt)


##### DATA PROCESSING #####
# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE179347", "file=GSE179347_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

# get metadata annotations for samples
X <- GEOquery::getGEO(GEO = "GSE179347")

simple_meta <- data.frame(accession = pData(X[[1]])[1:9, 2],
                          group = pData(X[[1]])[1:9, 48])
rownames(simple_meta) <- rownames(pData(X[[1]]))[1:9]
simple_meta$group_name <- c(rep("control", 3),
                            rep("no_glucose", 3),
                            rep("glucose", 3))


# pre-filter low count genes
# keep genes with at least 10 counts in at least 2 samples
keep <- rowSums( tbl >= 10 ) >= 2
tbl <- tbl[keep, ]

# keep only protein coding genes
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'www.ensembl.org')
genes <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name","transcript_biotype"), filters = c("transcript_biotype","chromosome_name"),values = list("protein_coding",c(1:22)), mart = mart)

z <- AnnotationDbi::select(x = org.Hs.eg.db, 
                           keys = as.character(rownames(tbl)), 
                           keytype = "ENTREZID", 
                           columns = c("SYMBOL", "GENENAME"))

gids <- which(z$SYMBOL %in% genes$external_gene_name)
E <- tbl[gids, ] # <-- final expression table
G <- z[gids, ] # <-- gene annotation data frame


##### Differential expression #####
dds <- DESeqDataSetFromMatrix(countData = E[,1:9],
                              colData = simple_meta,
                              design = ~ group_name)

dds <- DESeq(object = dds)

glc_noglc <- results(dds, 
                     contrast = c("group_name", "glucose", "no_glucose"),
                     pAdjustMethod = "fdr",
                     cooksCutoff = FALSE)

res <- tibble::tibble(entrez_id = rownames(glc_noglc),
                        logFC = as.vector(as.numeric(glc_noglc$log2FoldChange)),
                        p_val = as.vector(as.numeric(glc_noglc$pvalue)),
                        q_val = as.vector(as.numeric(glc_noglc$padj)))

diff_genes <- res %>% 
  tibble::add_column(gene_symbol = G$SYMBOL) %>% 
  tibble::add_column(gene_name = G$GENENAME) %>% 
  dplyr::filter(., q_val < 0.05, abs(logFC) > 1)


# plot differential expression results
res %>% 
  tibble::add_column(., color = "gray") %>%
  dplyr::mutate(., color = replace(color, list = which(logFC > 1 & q_val < 0.05), values = "red")) %>%
  dplyr::mutate(., color = replace(color, list = which(logFC < -1 & q_val < 0.05), values = "red")) %>%
  ggplot() + 
  geom_point(aes(x = logFC, y = -log10(q_val), color = color), alpha = 0.5, size = 3) +
  scale_color_manual(values = c("gray", "red")) +
  geom_hline(yintercept = -log10(0.05), lty = 2) +
  geom_vline(xintercept = -1, lty = 2) +
  geom_vline(xintercept = 1, lty = 2) +
  labs(x = "log2(fold change)", y = "-log10(p-value)") +
  theme_bw() +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 24, color = "black"),
        legend.position = "none",
        panel.grid = element_blank())


# plot differentially expressed genes only
diff_genes %>%
  # dplyr::arrange(desc(logFC)) %>%
  # tibble::add_column(rank_order = seq(1, nrow(diff_genes))) %>%
  ggplot(aes(x = logFC, y = fct_reorder(gene_symbol, logFC))) +
  geom_point(alpha = 0.5, size = 3) + 
  labs(x = "log2(fold change)", y = "Gene symbol") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 24, color = "black"),
        legend.position = "none",
        panel.grid = element_blank())





