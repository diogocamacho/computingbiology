##### LIBRARIES #####
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(biomaRt)
library(oranges)



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

##### LOAD GENE SETS #####
kegg_pathways <- oranges::cpdb_data$pathway_info$pathway_name[oranges::cpdb_data$pathway_info$pathway_source == "kegg"]
