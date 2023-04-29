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

## Load the data
load("../data/HTP_transcription_counts_wide_1000.Rdata")

htp_expr <- expression_list$expression

# Download example input
download.file(url = "https://raw.githubusercontent.com/nicolerg/ssGSEA2/master/example/PI3K_pert_logP_n2x23936.gct",              destfile = "/tmp/PI3K_pert_logP_n2x23936.gct")

# Download gene set database 
download.file(url = "https://raw.githubusercontent.com/nicolerg/ssGSEA2/master/example/ptm.sig.db.all.flanking.human.v1.8.1.gmt",destfile  = "/tmp/ptm.sig.db.all.flanking.human.v1.8.1.gmt")

res = run_ssGSEA2("/tmp/PI3K_pert_logP_n2x23936.gct",
                  output.prefix = "example",
                  gene.set.databases = "/tmp/ptm.sig.db.all.flanking.human.v1.8.1.gmt",
                  output.directory = "/tmp",
                  sample.norm.type = "none", 
                  weight = 0.75, 
                  correl.type = "rank", 
                  statistic = "area.under.RES",
                  output.score.type = "NES", 
                  nperm = 1000, 
                  min.overlap = 5, 
                  extended.output = TRUE, 
                  global.fdr = FALSE,
                  log.file = "/tmp/run.log")



## cmap exploration


my_ds <- parse_gctx("/tmp/PI3K_pert_logP_n2x23936.gct")




### Testing testing with data

## creating GCT object
htp_gct_mat <- htp_expr %>%
    column_to_rownames("LabID") %>%
    t()

htp_gct <- new("GCT", htp_gct_mat)

## writing gct object to disk
write_gct(htp_gct, "../data/HTP_transcription_counts.gct", appenddim = F)

## Hallmark Gene Sets

## testing

res = run_ssGSEA2("../data/HTP_transcription_counts.gct",
                  output.prefix = "test_hallmarks",
                  gene.set.databases = "../data/Gene Pathways/h.all.v2023.1.Hs.symbols.gmt",
                  output.directory = "../results/tests",
                  sample.norm.type = "none", 
                  weight = 0.75, 
                  correl.type = "rank", 
                  statistic = "area.under.RES",
                  output.score.type = "NES", 
                  nperm = 1000, 
                  min.overlap = 5, 
                  extended.output = TRUE, 
                  global.fdr = FALSE,
                  log.file = "../results/tests/run.log")
