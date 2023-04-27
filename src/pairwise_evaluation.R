#### This program is for performing pairwise gene evaluation
library(tidyverse)

load("../data/HTP_transcription_counts_wide_1000.Rdata")

htp_expr <- expression_list$expression

## column medians
gene_medians <- htp_expr %>%
    summarise_if(is.numeric, mean) %>%
    t() %>%
    as.numeric()

## find gene ratios

gene_ratios <- outer(gene_medians ,gene_medians, "/")
colnames(gene_ratios) <- names(htp_expr)[-1]
rownames(gene_ratios) <- names(htp_expr)[-1]

## pivot to longer
test <- gene_ratios %>%
    as.tibble() %>%
    mutate(Gene1 = rownames(gene_ratios))%>%
    pivot_longer(cols = colnames(gene_ratios),  names_to = "Gene2", values_to = "gene_ratio") %>%
    filter(Gene1 != Gene2) %>%
    arrange(desc(gene_ratio))
