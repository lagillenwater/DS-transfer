## GSEA
library(fgsea)
library(tidyverse)
set.seed(42)

myGO <- gmtPathways("./data/Gene Pathways/c2.all.v2023.1.Hs.symbols.gmt")
deg_res <- read.csv("./results/tcga/deg/prostate_deg.csv")

gene_list <- deg_res$logFC
names(gene_list) <- deg_res$ID
gene_list = gene_list[!duplicated(names(gene_list))]

gsea_res <- fgsea(pathways = myGO,
                  stats = gene_list,
                  minSize = 15,
                  maxSize = 10000)

gsea_res <- gsea_res %>%
    arrange(pval)
