library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(biomaRt)


A <- readr::read_tsv("~/Downloads/BIOGRID-ALL-4.4.232.tab3.txt") 

A <- A %>% 
  dplyr::filter(`Organism ID Interactor A` == 9606, `Organism ID Interactor B` == 9606) %>%
  dplyr::select(biogrid_id = `#BioGRID Interaction ID`,
                entrez_a = `Entrez Gene Interactor A`, 
                entrez_b = `Entrez Gene Interactor B`,
                gene_a = `Official Symbol Interactor A`,
                gene_b = `Official Symbol Interactor B`,
                score = Score, `Experimental System`, `Experimental System Type`, `Publication Source`, Throughput, Score) 

save(file = "~/work/data/biogrid_human.RData", A)

