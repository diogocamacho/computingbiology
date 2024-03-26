#' Tissue Specificity
#' 
#' Diogo M. Camacho
#' March 2024

# libraries
library(readxl)
library(readr)
library(biomaRt)  
library(ggplot2)

# set seed for reproducibility
set.seed(20240324)

# Data
# We will be using the TPM matrix from GTex V8 which you can download here:
#
# https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression

# metadata
sample_mapping <- readr::read_delim("~/work/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
subject_mapping <- readr::read_delim("~/work/data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

# TPM data
gtex_data <- readr::read_delim("~/work/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", delim = "\t", skip = 2)

# make expression data frame a matrix
G <- gtex_data$Description
M <- as.matrix(gtex_data[, 3:ncol(gtex_data)])

### memory saver
# selecting only a few samples, to make data matrix smaller to be more memory efficient for this example
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


# now clean matrix
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

### example on expression profile
# get 5 random genes
random_genes <- sample(x = 1:nrow(sampled_expression), 5)

# plot the expression of one of these
plot(sampled_expression[random_genes[3], ])
G[random_genes[3]]


### calculate sample specificity
# scale the data for any sample
column_scale <- scale(x = sampled_expression)

# scale the data for any gene
row_scale <- t(scale(t(x = sampled_expression)))


specificity_z <- (row_scale + column_scale) / sqrt(2)

### examples
# gapdh
gapdh_id <- which(G == "GAPDH")
# make a data frame to plot TPM vs z-score data

gapdh_df <- tibble::tibble(sample_name = colnames(sampled_expression),
                           tpm_data = sampled_expression[gapdh_id, ],
                           zscore_data = specificity_z[gapdh_id, ],
                           tissue_source = samples_df$tissue_source)

ggplot(gapdh_df, aes(x = tpm_data, y = zscore_data, color = tissue_source)) +
  geom_point() + 
  geom_hline(yintercept = 3, lty = 2, color = "darkblue") +
  theme_bw()

length(which(specificity_z[gapdh_id, ] > 3)) / ncol(specificity_z)

# GPL1R
glp_id <- which(G == "GLP1R")
# make a data frame to plot TPM vs z-score data

glp_df <- tibble::tibble(sample_name = colnames(sampled_expression),
                           tpm_data = sampled_expression[glp_id, ],
                           zscore_data = specificity_z[glp_id, ],
                           tissue_source = samples_df$tissue_source)

ggplot(glp_df, aes(x = tpm_data, y = zscore_data, color = tissue_source)) +
  geom_point() + 
  geom_hline(yintercept = 3, lty = 2, color = "darkblue") +
  theme_bw()

length(which(specificity_z[glp_id, ] > 3)) / ncol(specificity_z)


# ALK
alk_id <- which(G == "ALK")
# make a data frame to plot TPM vs z-score data

alk_df <- tibble::tibble(sample_name = colnames(sampled_expression),
                         tpm_data = sampled_expression[alk_id, ],
                         zscore_data = specificity_z[alk_id, ],
                         tissue_source = samples_df$tissue_source)

ggplot(alk_df, aes(x = tpm_data, y = zscore_data, color = tissue_source)) +
  geom_point() + 
  geom_hline(yintercept = 3, lty = 2, color = "darkblue") +
  theme_bw()

length(which(specificity_z[alk_id, ] > 3)) / ncol(specificity_z)





