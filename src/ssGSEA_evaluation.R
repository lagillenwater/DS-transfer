## This is a program for performing ssGSEA evaluation of the a gene expression dataset

## Install  the ssGSEA package
if (!require("devtools", quietly = TRUE)){
  install.packages("devtools")
}

if (!require("ssGSEA2", quietly = TRUE)){
devtools::install_github("nicolerg/ssGSEA2")
}

## Load  libraries
library(ssGSEA2)
library(tidyverse)
library(cmapR)
library(ComplexHeatmap)

## Load the data
load("../data/HTP_transcription_counts_wide_protein_coding_variance_filtered.Rdata")

htp_expr <- expression_list$expression


### Testing testing with data

## creating GCT object
htp_gct_mat <- htp_expr %>%
    column_to_rownames("LabID") %>%
    t()

htp_gct <- new("GCT", htp_gct_mat)

## writing gct object to disk
write_gct(htp_gct, "../data/HTP_transcription_counts_protein_coding_variance_filtered.gct", appenddim = F)

## Hallmark Gene Sets

## testing

res = run_ssGSEA2("../data/HTP_transcription_counts_protein_coding_variance_filtered.gct",
                  output.prefix = "test_hallmarks",
                  gene.set.databases = "../data/Gene Pathways/h.all.v2023.1.Hs.symbols.gmt",
                  output.directory = "../results/tests",
                  sample.norm.type = "none", 
                  weight = 0.75, 
                  correl.type = "rank", 
                  statistic = "area.under.RES",
                  output.score.type = "NES", 
                  nperm = 10, 
                  min.overlap = 5, 
                  extended.output = TRUE, 
                  global.fdr = FALSE,
                  log.file = "../results/tests/run.log")

scores <- parse_gctx("../results/tests/test_hallmarks-scores.gct")
scores_mat <- mat(scores)

### mapping the data

### loading metadata 
meta <- read.delim("../data/P4C_metadata_021921_Costello.txt")

### processing metadata
meta_heatmap <- meta %>%
    select(LabID, Karyotype) %>%
    filter(LabID %in% ids(scores, dimension = "column")) %>%
    arrange(match(LabID, ids(scores, dimension = "column"))) %>%
    column_to_rownames('LabID')

## fitlering scores matrix
scores_mat <- scores_mat[, colnames(scores_mat) %in% rownames(meta_heatmap)]

identical(rownames(meta_heatmap), colnames(scores_mat))

ha <- HeatmapAnnotation(df = meta_heatmap)

Heatmap(scores_mat,
        show_column_names = F,
        top_annotation = ha,
        row_names_gp = grid::gpar(fontsize = 7))
