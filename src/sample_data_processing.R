## Sample data processing

library(tidyverse)

## Create smaller sample datasets for testing

## Prostate Cancer Data and gtex data
tcga_meta <- read.csv("./data/tcga/PRAD_metadata.csv")
tcga_expr <- read.csv("./data/tcga/PRAD_expression.csv")
gtex_meta <- read.csv("./data/gtex/PROSTATE_metadata.csv")
gtex_expr <- read.csv("./data/gtex/PROSTATE_expression.csv")


## filter tcga data to tumor data only
tumor_expr <- tcga_expr %>%
    filter(X %in% ( tcga_meta %>% filter(tcga.cgc_sample_sample_type == "Primary Tumor") %>% filter(grepl("01A", tcga.cgc_sample_id)) %>% .$X)) 
## Also filtered based on the sample_id containing "01A" to remove duplicates

any(duplicated(tumor_expr$X)) ## Checking to ensure that there are not duplicated 

## combine the tcga and gtex expression data

identical(names(tumor_expr), names(gtex_expr)) ## TRUE
expr <- rbind(tumor_expr,gtex_expr)

## log transform
expr[,2:ncol(expr)] <- log2(expr[,2:ncol(expr)]+1)


## Apply variance filter to genes
gene_vars <- expr %>%
    summarise(across(where(is.numeric), var))


high_variance_genes <- gene_vars %>%
    pivot_longer(cols = everything(), names_to = "gene",values_to = "variance") %>%
    arrange(desc(variance)) %>%
    filter(variance > .1)

expr_filtered <- expr %>%
    dplyr::select(c('X', high_variance_genes$gene))


## Convert Ensembl IDs to Gene Names
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- names(expr_filtered)[2:ncol(expr_filtered)]
alternate <- gsub(".*[.]","", genes)
genes <- gsub("[.].*","", genes)
G_list <- getBM(filters= "ensembl_gene_id",attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
G_list <- G_list %>% filter(hgnc_symbol != "")

names(expr_filtered) <- gsub("[.].*","", names(expr_filtered))
expr_filtered <- expr_filtered %>%
    dplyr::select(c("X", G_list$ensembl_gene_id)) %>%
    subset(.,select = which(!duplicated(names(.))))
names(expr_filtered)[2:ncol(expr_filtered)] <- G_list$hgnc_symbol


## Export the expression data based on applying a variance filter
save(expr_filtered, file = "./data/processed/prostate_gene_expression.RData")

## Combine and filter the metadata
meta <- gtex_meta %>%
    full_join(tcga_meta, by = "X", .keep_all = T) %>%
    filter(X %in% expr_filtered$X)

## Save the metadata
save(meta, file = "./data/processed/prostate_metadata.RData")
