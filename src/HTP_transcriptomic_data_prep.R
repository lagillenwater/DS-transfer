#### This program preps the HTP data for downstream analysis
## Pipeline:
##    - load count data
##    - log transform
##    - apply variance filter
##    - pivot to wider

library(tidyverse)

## load expression data
expression <- read.delim("../data/HTP_WholeBlood_RNAseq_Counts_Synapse.txt")

## logtranformation
expression$logValue <- log(expression$Value+1)

## apply variance filter to protein coding genes
## for now I'm selecting the top 1000 genes for testing purposes
top_1000 <- expression %>%
    filter(Gene_type == "protein_coding") %>%
    group_by(EnsemblID) %>%
    summarise(variance = var(logValue)) %>%
    arrange(desc(variance)) %>%
    top_n(1000) 


variance_filtered <- expression %>%
    filter(EnsemblID %in% top_1000$EnsemblID) %>%
    as_tibble() 
    



## pivot wider. Have to include the unique identifier for row numbers to avoid any errors. 
htp_expr <- variance_filtered %>%
    select(LabID, Gene_name,logValue) %>%
    pivot_wider(names_from = Gene_name, id_cols = LabID, values_from = logValue)




### get the mappings of gene_name, ensembl, and chr from expression


gene_map <- expression %>%
    select(Gene_name, EnsemblID,Chr) %>%
    filter(EnsemblID %in% top_1000$EnsemblID) %>%
    as_tibble()

expression_list<- list(expression = htp_expr,gene_map = gene_map)
save(expression_list, file = "../data/HTP_transcription_counts_wide_1000.Rdata")


