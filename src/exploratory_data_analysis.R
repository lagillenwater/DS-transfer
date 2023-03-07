## Sample data processing

library(tidyverse)


## Create smaller sample datasets for testing

## Prostate 
meta <- read.csv("./data/tcga/PRAD_metadata.csv")
expr <- read.csv("./data/tcga/PRAD_expression.csv")

## Identify  ID's with primary tumor and normal tissue samples
tcga_ids <- meta %>% select(tcga.cgc_case_id,tcga.cgc_sample_sample_type)  %>% filter(tcga.cgc_sample_sample_type == "Solid Tissue Normal") %>% .$tcga.cgc_case_id

## sample ids are different from tcga ids
sample_ids <- meta %>% filter(tcga.cgc_case_id %in% tcga_ids) %>% .$X

## filter expression  data and metadata by those ids
sample_expression <- expr %>% filter(X %in% sample_ids)
sample_meta <- meta %>% filter(tcga.cgc_case_id %in% tcga_ids)

## Apply variance filter to genes
gene_vars <- sample_expression %>%
    summarise(across(where(is.numeric), var))

high_variance_genes <- gene_vars %>%
    pivot_longer(cols = everything(), names_to = "gene",values_to = "variance") %>%
    arrange(desc(variance)) %>%
    top_n(5000)

